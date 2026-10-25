// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "HydroMechanicsProcess.h"

#include <cassert>
#include <set>
#include "NumLib/DOF/GlobalMatrixProviders.h"
#include "NumLib/Exceptions.h"
#include "MathLib/LinAlg/FinalizeVectorAssembly.h"

#include "HydroMechanicsFEM.h"
#include "MaterialLib/SolidModels/MFront/MFront.h"
#include "HydroMechanicsProcessData.h"
#include "MeshLib/Elements/Utils.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/Deformation/SolidMaterialInternalToSecondaryVariables.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/Utils/CreateLocalAssemblersTaylorHood.h"
#include "ProcessLib/Utils/SetIPDataInitialConditions.h"

namespace ProcessLib
{
namespace HydroMechanics
{
template <int DisplacementDim>
HydroMechanicsProcess<DisplacementDim>::HydroMechanicsProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    HydroMechanicsProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const use_monolithic_scheme,
    bool const is_linear)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      AssemblyMixin<HydroMechanicsProcess<DisplacementDim>>{
          *_jacobian_assembler, is_linear, use_monolithic_scheme},
      process_data_(std::move(process_data))
{
    updateSafetyMetadata(0.);
    // For numerical Jacobian
    if (this->_jacobian_assembler->needsPicardAssembly() &&
        _use_monolithic_scheme)
    {
        OGS_FATAL(
            "Numerical Jacobian is not supported for the "
            "HydroMechanicsProcess using the monolithic scheme yet.");
    }

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "sigma_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, local_assemblers_, &LocalAssemblerIF::getSigma));

    _integration_point_writer.emplace_back(
        std::make_unique<MeshLib::IntegrationPointWriter>(
            "epsilon_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            integration_order, local_assemblers_,
            &LocalAssemblerIF::getEpsilon));

    if (process_data_.mechanical_drainage)
    {
        for (auto const& [name, getter] : std::array{
            std::pair{"mechanical_pressure_ip", &LocalAssemblerIF::getMechanicalPressure},
            std::pair{"mechanical_excess_pressure_ip", &LocalAssemblerIF::getMechanicalExcess},
            std::pair{"mechanical_fluid_content_change_ip", &LocalAssemblerIF::getMechanicalFluidContent},
            std::pair{"mechanical_reference_strain_ip", &LocalAssemblerIF::getMechanicalReferenceStrain},
            std::pair{"mechanical_reference_pressure_ip", &LocalAssemblerIF::getMechanicalReferencePressure}})
        {
            _integration_point_writer.emplace_back(
                std::make_unique<MeshLib::IntegrationPointWriter>(name, 1,
                    integration_order, local_assemblers_, getter));
        }
    }
    if (process_data_.use_positive_pore_pressure_in_mechanics && !process_data_.mechanical_drainage)
    {
        _integration_point_writer.emplace_back(
            std::make_unique<MeshLib::IntegrationPointWriter>(
                "mechanical_pressure_ip", 1, integration_order,
                local_assemblers_, &LocalAssemblerIF::getMechanicalPressure));
    }
    if (process_data_.use_positive_pore_pressure_in_mechanics)
    {
        _integration_point_writer.emplace_back(
            std::make_unique<MeshLib::IntegrationPointWriter>(
                "applied_biot_pressure_ip", 1, integration_order,
                local_assemblers_, &LocalAssemblerIF::getAppliedBiotPressure));
    }
    if (process_data_.safety && process_data_.safety->usesActualPressure())
    {
        _integration_point_writer.emplace_back(std::make_unique<MeshLib::IntegrationPointWriter>(
            "mechanical_actual_reference_pressure_ip",1,integration_order,
            local_assemblers_,&LocalAssemblerIF::getMechanicalActualReferencePressure));
    }
    if (process_data_.mechanical_drainage &&
        std::any_of(process_data_.mechanical_drainage->materials.begin(),
                    process_data_.mechanical_drainage->materials.end(),
                    [](auto const& entry) { return entry.second.usesTotalStress(); }))
    {
        for (auto const& [name, getter] : std::array{
            std::pair{"mechanical_stress_basis_ip", &LocalAssemblerIF::getMechanicalStressBasis},
            std::pair{"mechanical_pressure_valid_ip", &LocalAssemblerIF::getMechanicalPressureValid}})
        {
            _integration_point_writer.emplace_back(
                std::make_unique<MeshLib::IntegrationPointWriter>(name, 1,
                    integration_order, local_assemblers_, getter));
        }
    }
    if (!_use_monolithic_scheme)
    {
        _integration_point_writer.emplace_back(
            std::make_unique<MeshLib::IntegrationPointWriter>(
                "strain_rate_variable_ip", 1, integration_order,
                local_assemblers_, &LocalAssemblerIF::getStrainRateVariable));
    }
}

template <int DisplacementDim>
bool HydroMechanicsProcess<DisplacementDim>::isLinear() const
{
    return AssemblyMixin<HydroMechanicsProcess<DisplacementDim>>::isLinear();
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
HydroMechanicsProcess<DisplacementDim>::getMatrixSpecifications(
    const int process_id) const
{
    // For the monolithic scheme or the M process (deformation) in the staggered
    // scheme.
    if (process_id == process_data_.mechanics_related_process_id)
    {
        auto const& l = *_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }

    // For staggered scheme and H process (pressure).
    auto const& l = *local_to_global_index_map_with_base_nodes_;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &sparsity_pattern_with_linear_element_};
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes = std::make_unique<MeshLib::MeshSubset>(
        _mesh, _mesh.getNodes(), process_data_.use_taylor_hood_elements);

    // Create single component dof in the mesh's base nodes.
    base_nodes_ = MeshLib::getBaseNodes(_mesh.getElements());
    mesh_subset_base_nodes_ = std::make_unique<MeshLib::MeshSubset>(
        _mesh, base_nodes_, process_data_.use_taylor_hood_elements);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};
    local_to_global_index_map_single_component_ =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    if (process_data_.isMonolithicSchemeUsed())
    {
        // For pressure, which is the first
        std::vector<MeshLib::MeshSubset> all_mesh_subsets{
            *mesh_subset_base_nodes_};

        // For displacement.
        const int monolithic_process_id = 0;
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        getProcessVariables(monolithic_process_id)[1]
                            .get()
                            .getNumberOfGlobalComponents(),
                        [&]() { return *_mesh_subset_all_nodes; });

        std::vector<int> const vec_n_components{1, DisplacementDim};
        _local_to_global_index_map =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);
        assert(_local_to_global_index_map);
    }
    else
    {
        // For displacement equation.
        const int process_id = 1;
        std::vector<MeshLib::MeshSubset> all_mesh_subsets;
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        getProcessVariables(process_id)[0]
                            .get()
                            .getNumberOfGlobalComponents(),
                        [&]() { return *_mesh_subset_all_nodes; });

        std::vector<int> const vec_n_components{DisplacementDim};
        _local_to_global_index_map =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);

        // For pressure equation.
        // Collect the mesh subsets with base nodes in a vector.
        std::vector<MeshLib::MeshSubset> all_mesh_subsets_base_nodes{
            *mesh_subset_base_nodes_};
        local_to_global_index_map_with_base_nodes_ =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets_base_nodes),
                // by location order is needed for output
                NumLib::ComponentOrder::BY_LOCATION);

        sparsity_pattern_with_linear_element_ = NumLib::computeSparsityPattern(
            *local_to_global_index_map_with_base_nodes_, _mesh);

        assert(_local_to_global_index_map);
        assert(local_to_global_index_map_with_base_nodes_);
    }
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::createLocalAssemblersHM<DisplacementDim,
                                        HydroMechanicsLocalAssembler>(
        mesh.getElements(), dof_table, local_assemblers_,
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        process_data_);

    auto add_secondary_variable = [&](std::string const& name,
                                      int const num_components,
                                      auto get_ip_values_function)
    {
        _secondary_variables.addSecondaryVariable(
            name,
            makeExtrapolator(num_components, getExtrapolator(),
                             local_assemblers_,
                             std::move(get_ip_values_function)));
    };

    add_secondary_variable("sigma",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerIF::getIntPtSigma);

    add_secondary_variable("epsilon",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerIF::getIntPtEpsilon);

    add_secondary_variable("velocity", DisplacementDim,
                           &LocalAssemblerIF::getIntPtDarcyVelocity);

    //
    // enable output of internal variables defined by material models
    //
    ProcessLib::Deformation::solidMaterialInternalToSecondaryVariables<
        LocalAssemblerIF>(process_data_.solid_materials,
                          add_secondary_variable);

    if(!process_data_.excavation_release.empty())
    {
        process_data_.excavation_release_force=MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh),"excavation_release_force",MeshLib::MeshItemType::Node,DisplacementDim);
        std::fill(process_data_.excavation_release_force->begin(),process_data_.excavation_release_force->end(),0.);
    }

    if (!process_data_.displacement_reference_times.empty())
    {
        process_data_.displacement_relative = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "displacement_relative",
            MeshLib::MeshItemType::Node, DisplacementDim);
        process_data_.displacement_reference_values.assign(mesh.getNumberOfNodes() * DisplacementDim, 0.);
    }

    process_data_.pressure_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "pressure_interpolated",
            MeshLib::MeshItemType::Node, 1);

    process_data_.principal_stress_vector[0] =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_vector_1",
            MeshLib::MeshItemType::Cell, 3);

    process_data_.principal_stress_vector[1] =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_vector_2",
            MeshLib::MeshItemType::Cell, 3);

    process_data_.principal_stress_vector[2] =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_vector_3",
            MeshLib::MeshItemType::Cell, 3);

    process_data_.principal_stress_values =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "principal_stress_values",
            MeshLib::MeshItemType::Cell, 3);

    process_data_.permeability = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "permeability",
        MeshLib::MeshItemType::Cell,
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim));

    if (process_data_.output_material_state)
    {
        ProcessLib::Deformation::solidMaterialInternalVariablesToIntegrationPointWriter(
            process_data_.solid_materials, local_assemblers_,
            _integration_point_writer, integration_order);
    }

    setIPDataInitialConditions(_integration_point_writer, mesh.getProperties(),
                               local_assemblers_);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(&LocalAssemblerIF::initialize,
                                                local_assemblers_,
                                                *_local_to_global_index_map);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::initializeBoundaryConditions(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    if (process_data_.isMonolithicSchemeUsed())
    {
        const int process_id_of_hydromechanics = 0;
        initializeProcessBoundaryConditionsAndSourceTerms(
            *_local_to_global_index_map, process_id_of_hydromechanics, media);
        return;
    }

    // Staggered scheme:
    // for the equations of pressure
    const int hydraulic_process_id = 0;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *local_to_global_index_map_with_base_nodes_, hydraulic_process_id,
        media);

    // for the equations of deformation.
    const int mechanical_process_id = 1;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map, mechanical_process_id, media);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for HydroMechanics");

    // Note: This assembly function is for the Picard nonlinear solver. Since
    // only the Newton-Raphson method is employed to simulate coupled HM
    // processes in this class, this function is actually not used so far.

    if (this->_jacobian_assembler->needsPicardAssembly() &&
        (!_use_monolithic_scheme))
    {
        if (process_id == process_data_.hydraulic_process_id)
        {
            this->_jacobian_assembler->setNonDeformationComponentIDs(
                {process_id} /* pressure variable id */);
        }
        else
        {
            OGS_FATAL(
                "Numerical Jacobian is only supported for the "
                "liquid fluid process in the staggered HydroMechanicsProcess.");
        }
    }

    AssemblyMixin<HydroMechanicsProcess<DisplacementDim>>::assemble(
        t, dt, x, x_prev, process_id, M, K, b);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::assembleResidual(
    double const t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalVector& b, GlobalMatrix& scratch)
{
    if (!supportsResidualScreening())
    { OGS_FATAL("Residual screening requires ordinary monolithic 2D HM without strict drainage or Safety."); }
    struct ResetFlag
    {
        bool& flag;
        bool previous;
        ~ResetFlag() { flag = previous; }
    } reset{process_data_.residual_screening_assembly,
            process_data_.residual_screening_assembly};
    process_data_.residual_screening_assembly = true;
    assembleWithJacobian(t, dt, x, x_prev, process_id, b, scratch);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalVector& b, GlobalMatrix& Jac)
{
    if (_boundary_conditions[process_id].hasActiveSet())
    { seepage_time_=t;seepage_dt_=dt;seepage_previous_=x_prev; }
    // For the monolithic scheme
    bool const use_monolithic_scheme = process_data_.isMonolithicSchemeUsed();
    if (use_monolithic_scheme)
    {
        DBUG(
            "Assemble the Jacobian of HydroMechanics for the monolithic "
            "scheme.");
    }
    else
    {
        // For the staggered scheme
        if (process_id == process_data_.hydraulic_process_id)
        {
            DBUG(
                "Assemble the Jacobian equations of liquid fluid process in "
                "HydroMechanics for the staggered scheme.");
        }
        else
        {
            DBUG(
                "Assemble the Jacobian equations of mechanical process in "
                "HydroMechanics for the staggered scheme.");
        }
    }

    AssemblyMixin<HydroMechanicsProcess<DisplacementDim>>::assembleWithJacobian(
        t, dt, x, x_prev, process_id, b, Jac);
    applyExcavationRelease(t,b);
}

template <int DisplacementDim>
NumLib::IterationResult HydroMechanicsProcess<DisplacementDim>::postIterationConcreteProcess(GlobalVector const& x)
{
    if (!_boundary_conditions[0].hasEnabledActiveSet(seepage_time_))
    {
        return NumLib::IterationResult::SUCCESS;
    }
    auto const specification=getMatrixSpecifications(0);std::size_t matrix_id=0;
    auto& jac=NumLib::GlobalMatrixProvider::provider.getMatrix(specification,matrix_id);
    auto& rhs=NumLib::GlobalVectorProvider::provider.getVector(specification);
    jac.setZero();rhs.setZero();std::vector<GlobalVector*> current{const_cast<GlobalVector*>(&x)};
    try
    { this->assembleWithJacobian(seepage_time_,seepage_dt_,current,seepage_previous_,0,rhs,jac); }
    catch (NumLib::AssemblyException const& error)
    {
        WARN("Seepage trial reassembly failed; reject global increment: {}",error.what());
        NumLib::GlobalMatrixProvider::provider.releaseMatrix(jac);
        NumLib::GlobalVectorProvider::provider.releaseVector(rhs);
        return NumLib::IterationResult::FAILURE;
    }
    bool const changed=_boundary_conditions[0].updateActiveSet(seepage_time_,x);
    NumLib::GlobalMatrixProvider::provider.releaseMatrix(jac);NumLib::GlobalVectorProvider::provider.releaseVector(rhs);
    if (changed) { INFO("HM seepage active set changed; repeat Newton before acceptance."); }
    return changed?NumLib::IterationResult::REPEAT_ITERATION:NumLib::IterationResult::SUCCESS;
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::updateSafetyMetadata(double const coordinate)
{
    if (!process_data_.steady_state_seepage.intervals.empty())
    {
        auto* field = MeshLib::getOrCreateMeshProperty<double>(
            _mesh, "steady_state_seepage", MeshLib::MeshItemType::Node, 1);
        std::fill(field->begin(), field->end(),
                  process_data_.steady_state_seepage.active(coordinate) ? 1. : 0.);
    }
    if (!process_data_.transient_seepage.intervals.empty())
    {
        auto* field = MeshLib::getOrCreateMeshProperty<double>(
            _mesh, "transient_seepage", MeshLib::MeshItemType::Node, 1);
        std::fill(field->begin(), field->end(),
                  process_data_.transient_seepage.active(coordinate) ? 1. : 0.);
    }
    if (!process_data_.safety) { return; }
    auto const& settings=*process_data_.safety;
    double const rate_valid=!settings.working_hydraulic_rates || settings.pressure_policy=="NotApplicable" ||
        NumLib::Time{coordinate}<=NumLib::Time{settings.initialisation_end} || settings.active(coordinate) || process_data_.steady_state_seepage.active(coordinate) ? 0.:1.;
    for (auto const& [name,value]:std::array{
        std::pair{"physical_time",std::min(coordinate,settings.physical_time)},
        std::pair{"strength_reduction_factor",settings.factor(coordinate)},
        std::pair{"hydraulic_rate_valid",rate_valid}})
    {
        auto* field=MeshLib::getOrCreateMeshProperty<double>(_mesh,name,MeshLib::MeshItemType::Node,1);
        std::fill(field->begin(),field->end(),value);
    }
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep HydroMechanicsProcess.");
    if (process_data_.safety && !sourceTermsSupportSafety(process_id))
    { OGS_FATAL("Structural source requires a qualified explicit Safety strength/resistance policy."); }
    if(process_data_.native_checkpoint && !sourceTermsSupportCheckpoint(process_id))
    { OGS_FATAL("Structural source lacks qualified checkpoint state."); }
    updateSafetyMetadata(t);
    if (process_data_.safety)
    {
        auto const& settings=*process_data_.safety;
        setSourceSafetyState(process_id,settings.active(t),settings.factor(t),settings.physical_time);
        for (auto const& [id,material]:process_data_.solid_materials)
        {
            auto* mfront=dynamic_cast<MaterialLib::Solids::MFront::MFront<DisplacementDim>*>(material.get());
            if (!mfront) { OGS_FATAL("Safety requires MFront MC relations."); }
            mfront->setStrengthReductionFactors(settings.factor(t-dt),settings.factor(t));
        }
        updateSafetyMetadata(t);
    }


    if (hasMechanicalProcess(process_id))
    {
        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &LocalAssemblerIF::preTimestep, local_assemblers_,
            getActiveElementIDs(), *_local_to_global_index_map, *x[process_id],
            t, dt);

        AssemblyMixin<
            HydroMechanicsProcess<DisplacementDim>>::updateActiveElements();
    if (process_data_.adaptive_rollback)
    {
        adaptive_ip_snapshot_.clear();
        for (auto const& assembler:local_assemblers_) { adaptive_ip_snapshot_.push_back(assembler->captureIntegrationPointCheckpoint()); }
        adaptive_bc_snapshot_=_boundary_conditions[process_id].checkpointState();adaptive_snapshot_time_=t-dt;
    }
    prepareExcavationRelease(x,t,dt,process_id);
    }
}

template <int DisplacementDim>
std::vector<std::vector<std::string>>
HydroMechanicsProcess<DisplacementDim>::initializeAssemblyOnSubmeshes(
    std::vector<std::reference_wrapper<MeshLib::Mesh>> const& meshes)
{
    INFO("HydroMechanicsProcess initializeSubmeshOutput().");
    std::vector<std::vector<std::string>> per_process_residuum_names;
    if (_process_variables.size() == 1)  // monolithic
    {
        per_process_residuum_names = {{"MassFlowRate", "NodalForces"}};
    }
    else  // staggered
    {
        per_process_residuum_names = {{"MassFlowRate"}, {"NodalForces"}};
    }

    AssemblyMixin<HydroMechanicsProcess<DisplacementDim>>::
        initializeAssemblyOnSubmeshes(meshes, per_process_residuum_names);

    return per_process_residuum_names;
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::postAcceptedTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, int const process_id)
{
    if(!process_data_.excavation_release.empty())
    {
        excavation_accepted_time_=t;
        excavation_accepted_solution_.resize(x[process_id]->size());
        for(GlobalIndexType i=0;i<x[process_id]->size();++i)
        { excavation_accepted_solution_[i]=x[process_id]->get(i); }
        for(auto const& interval:process_data_.excavation_release)
        {
            if(std::abs(t-interval.start)<1e-10)
            {
                excavation_accepted_residual_=excavationReferenceResidual(*x[process_id],t,1.,process_id);
                break;
            }
        }
    }
    if (process_data_.displacement_relative &&
        process_data_.displacement_reference_next < process_data_.displacement_reference_times.size())
    {
        double const reference = process_data_.displacement_reference_times[process_data_.displacement_reference_next];
        double const tolerance = 1e-12 * std::max(1., reference);
        if (t > reference + tolerance)
        { OGS_FATAL("Accepted time grid skipped displacement reference t={}", reference); }
        if (std::abs(t - reference) <= tolerance)
        {
            for (auto const* node : _mesh.getNodes())
            {
                for (int c = 0; c < DisplacementDim; ++c)
                {
                    auto const index = getDOFTable(process_id).getGlobalIndex(
                        {_mesh.getID(), MeshLib::MeshItemType::Node, node->getID()}, 1, c);
                    if (index < 0 || index == NumLib::MeshComponentMap::nop)
                    { OGS_FATAL("Displacement references require complete serial displacement DOFs."); }
                    process_data_.displacement_reference_values[node->getID() * DisplacementDim + c] = x[process_id]->get(index);
                    // Secondary variables were computed before global acceptance;
                    // update the boundary-time output after the accepted capture.
                    (*process_data_.displacement_relative)[node->getID() * DisplacementDim + c] = 0.;
                }
            }
            ++process_data_.displacement_reference_next;
            INFO("HM captured accepted displacement reference at t={}", t);
        }
    }
    if (process_data_.mechanical_drainage)
    {
        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &LocalAssemblerIF::acceptMechanicalDrainage, local_assemblers_,
            getActiveElementIDs(), t, *x[process_id], getDOFTable(process_id));
    }
    if (process_data_.native_checkpoint)
    {
        auto& settings=*process_data_.native_checkpoint;
        auto values=[&](GlobalVector const& v)
        {
            if (v.getRangeBegin()!=0 || v.getRangeEnd()!=v.size())
            { OGS_FATAL("Native HM checkpoint requires serial complete vectors."); }
            std::vector<double> result(v.size());
            for (GlobalIndexType i=0;i<v.size();++i) { result[i]=v.get(i); }
            return result;
        };
        if (std::abs(t-settings.initialisation_time)<=1e-12*std::max(1.,t))
        { settings.initialisation_solution=values(*x[process_id]); }
        for (auto const& point:settings.points)
        {
            if (std::abs(t-point.time)>1e-12*std::max(1.,point.time)) { continue; }
            NativeCheckpointState state;
            state.time=t;state.initialisation_time=settings.initialisation_time;
            state.id=point.id;state.model_identity=point.model_identity;
            state.runtime_identity=NativeCheckpointDetail::runtimeIdentity();
            state.solution=values(*x[process_id]);state.initialisation_solution=settings.initialisation_solution;
            if (state.initialisation_solution.size()!=state.solution.size())
            { OGS_FATAL("Native HM checkpoint has no accepted initialisation reference."); }
            state.displacement_reference=process_data_.displacement_reference_values;
            state.displacement_reference_next=process_data_.displacement_reference_next;
            auto property=[&](char const* name)
            {
                if (!_mesh.getProperties().hasPropertyVector(name)) { return std::vector<double>{}; }
                auto const* p=_mesh.getProperties().template getPropertyVector<double>(name);
                return std::vector<double>(p->begin(),p->end());
            };
            state.nodal_forces=property("NodalForces");state.mass_flow_rate=property("MassFlowRate");
            for (auto const& assembler:local_assemblers_)
            { state.elements.push_back(assembler->captureIntegrationPointCheckpoint()); }
            state.boundary_states=_boundary_conditions[process_id].checkpointState();
            state.source_states=sourceCheckpointState(process_id);
            if (!process_data_.excavation_release.empty())
            {
                double const cursor = excavation_release_captured_ == std::numeric_limits<std::size_t>::max()
                    ? -1. : static_cast<double>(excavation_release_captured_);
                state.excavation_states = {{excavation_accepted_time_, cursor},
                    excavation_accepted_solution_, excavation_accepted_residual_, excavation_release_values_};
            }
            writeNativeCheckpoint(point.file,state);
            INFO("HM wrote accepted checkpoint {} at t={} to {}",point.id,point.time,point.file);
        }
    }
}

template <int DisplacementDim>
bool HydroMechanicsProcess<DisplacementDim>::restoreInitialConditionsConcreteProcess(
    std::vector<GlobalVector*> const& x, std::vector<GlobalVector*> const& x_prev,
    double const t, int const process_id)
{
    if (!process_data_.native_checkpoint || process_data_.native_checkpoint->restart_file.empty()) { return false; }
#ifdef USE_PETSC
    OGS_FATAL("Native HM checkpoint first scope requires serial Eigen vectors.");
#else
    auto& settings=*process_data_.native_checkpoint;
    auto const state=decodeNativeCheckpoint(NativeCheckpointDetail::fileBytes(settings.restart_file));
    if (settings.restored || state.time!=t || state.initialisation_time!=settings.initialisation_time ||
        state.model_identity!=settings.restart_identity ||
        state.runtime_identity!=NativeCheckpointDetail::runtimeIdentity() ||
        state.solution.size()!=static_cast<std::size_t>(x[process_id]->size()) ||
        state.initialisation_solution.size()!=state.solution.size() ||
        state.elements.size()!=local_assemblers_.size() ||
        state.displacement_reference.size()!=process_data_.displacement_reference_values.size() ||
        state.displacement_reference_next>process_data_.displacement_reference_times.size())
    { OGS_FATAL("Native HM checkpoint time/model/runtime/vector/reference layout mismatch."); }
    std::size_t expected_references=0;
    for (double const reference:process_data_.displacement_reference_times)
    { if (reference<=t+1e-12*std::max(1.,t)) { ++expected_references; } }
    if (state.displacement_reference_next!=expected_references)
    { OGS_FATAL("Native HM checkpoint displacement-reference cursor mismatch."); }
    auto check_property=[&](std::vector<double> const& values, std::size_t const components)
    { if (!values.empty() && values.size()!=_mesh.getNumberOfNodes()*components)
      { OGS_FATAL("Native HM checkpoint nodal diagnostic layout mismatch."); } };
    check_property(state.nodal_forces,DisplacementDim);check_property(state.mass_flow_rate,1);
    if (!process_data_.excavation_release.empty())
    {
        auto const& fields = state.excavation_states;
        if (fields.size()!=4 || fields[0].size()!=2 || fields[0][0]!=t ||
            fields[1].size()!=state.solution.size() ||
            (!fields[2].empty() && fields[2].size()!=state.solution.size()) ||
            (!fields[3].empty() && fields[3].size()!=state.solution.size()) ||
            fields[0][1]<-1. || std::floor(fields[0][1])!=fields[0][1] ||
            fields[0][1]>=static_cast<double>(process_data_.excavation_release.size()))
        { OGS_FATAL("Native HM checkpoint excavation state/layout mismatch."); }
    }
    else if (!state.excavation_states.empty())
    { OGS_FATAL("Excavation checkpoint requires the original release schedule."); }
    _boundary_conditions[process_id].validateCheckpointState(state.boundary_states,t);
    // All element payloads must validate before any global/IP state changes.
    for (std::size_t i=0;i<local_assemblers_.size();++i)
    { local_assemblers_[i]->validateIntegrationPointCheckpoint(state.elements[i]); }
    for (std::size_t i=0;i<local_assemblers_.size();++i)
    { local_assemblers_[i]->restoreIntegrationPointCheckpoint(state.elements[i]); }
    for (GlobalIndexType i=0;i<x[process_id]->size();++i)
    { x[process_id]->set(i,state.solution[i]);x_prev[process_id]->set(i,state.solution[i]); }
    MathLib::LinAlg::finalizeAssembly(*x[process_id]);MathLib::LinAlg::finalizeAssembly(*x_prev[process_id]);
    settings.initialisation_solution=state.initialisation_solution;
    process_data_.displacement_reference_values=state.displacement_reference;
    process_data_.displacement_reference_next=state.displacement_reference_next;
    GlobalVector initialisation_reference(x[process_id]->size());
    for (GlobalIndexType i=0;i<initialisation_reference.size();++i)
    { initialisation_reference.set(i,state.initialisation_solution[i]); }
    MathLib::LinAlg::finalizeAssembly(initialisation_reference);
    MathLib::LinAlg::setLocalAccessibleVector(initialisation_reference);
    std::vector<GlobalVector*> references{&initialisation_reference};
    _boundary_conditions[process_id].postAcceptedTimestep(state.initialisation_time,references,process_id);
    _boundary_conditions[process_id].restoreCheckpointState(state.boundary_states,t);
    restoreSourceCheckpointState(process_id,state.source_states);
    if(process_data_.safety)
    {
        auto const& safety=*process_data_.safety;
        setSourceSafetyState(process_id,safety.active(t),safety.factor(t),safety.physical_time);
    }
    auto restore_property=[&](char const* name,std::vector<double> const& values,int const components)
    {
        if (values.empty()) { return; }
        auto* property=MeshLib::getOrCreateMeshProperty<double>(_mesh,name,MeshLib::MeshItemType::Node,components);
        std::copy(values.begin(),values.end(),property->begin());
    };
    restore_property("NodalForces",state.nodal_forces,DisplacementDim);
    restore_property("MassFlowRate",state.mass_flow_rate,1);
    updateSafetyMetadata(t);
    if (!process_data_.excavation_release.empty())
    {
        auto const& fields = state.excavation_states;
        excavation_accepted_time_=fields[0][0];
        excavation_release_captured_=fields[0][1]<0. ? std::numeric_limits<std::size_t>::max()
            : static_cast<std::size_t>(fields[0][1]);
        excavation_accepted_solution_=fields[1];
        excavation_accepted_residual_=fields[2];
        excavation_release_values_=fields[3];
    }
    settings.restored=true;
    INFO("HM restored complete accepted checkpoint {} at physical t={}",state.id,t);
    return true;
#endif
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, double const t, double const dt,
    const int process_id)
{
    if (process_id != process_data_.hydraulic_process_id)
    {
        return;
    }

    DBUG("PostTimestep HydroMechanicsProcess.");

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::postTimestep, local_assemblers_,
        getActiveElementIDs(), getDOFTables(x.size()), x, x_prev, t, dt,
        process_id);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::postNonLinearSolverConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, const double t, double const dt,
    const int process_id)
{
    DBUG("PostNonLinearSolver HydroMechanicsProcess.");

    // Calculate strain, stress or other internal variables of mechanics.
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::postNonLinearSolver, local_assemblers_,
        getActiveElementIDs(), getDOFTables(x.size()), x, x_prev, t, dt,
        process_id);
    if (process_data_.initialisation_gate &&
        (NumLib::Time{t} == NumLib::Time{process_data_.initialisation_gate->reference_time} ||
         (process_data_.safety && process_data_.safety->verify_equilibrium && process_data_.safety->active(t))))
    {
        auto const& gate = *process_data_.initialisation_gate;
        auto const specification = getMatrixSpecifications(process_id);
        std::size_t matrix_id = 0;
        auto& jac = NumLib::GlobalMatrixProvider::provider.getMatrix(specification, matrix_id);
        auto& residual = NumLib::GlobalVectorProvider::provider.getVector(specification);
        jac.setZero();
        residual.setZero();
        // Reassemble at the final nonlinear solution, before history commit.
        // This includes natural tractions and sources, without suppressing supports.
        this->assembleWithJacobian(t, dt, x, x_prev, process_id, residual, jac);
        MathLib::finalizeVectorAssembly(residual);
        std::set<GlobalIndexType> essential;
        for (auto const& bc : *getKnownSolutions(t, *x[process_id], process_id))
        { essential.insert(bc.ids.begin(), bc.ids.end()); }
        double norm2 = 0.;
        for (auto const* node : _mesh.getNodes())
        {
            for (int c = 0; c < DisplacementDim; ++c)
            {
                auto const index = getDOFTable(process_id).getGlobalIndex(
                    {_mesh.getID(), MeshLib::MeshItemType::Node, node->getID()}, 1, c);
                if (index < 0 || index == NumLib::MeshComponentMap::nop)
                { OGS_FATAL("Initialisation gate currently requires complete serial mechanical DOFs."); }
                if (!essential.contains(index)) { norm2 += residual.get(index)*residual.get(index); }
            }
        }
        double worst_yield = -std::numeric_limits<double>::infinity();
        for (auto const& assembler : local_assemblers_)
        { worst_yield = std::max(worst_yield, assembler->getInitialisationYieldResidual(process_data_.safety?process_data_.safety->factor(t):1.)); }
        NumLib::GlobalMatrixProvider::provider.releaseMatrix(jac);
        NumLib::GlobalVectorProvider::provider.releaseVector(residual);
        double const norm = std::sqrt(norm2);
        INFO("HM initialisation gate at t={}: free_force={} N/m, MC_yield={} Pa.", t, norm, worst_yield);
        if (!std::isfinite(norm) || !std::isfinite(worst_yield) ||
            norm > gate.force_tolerance || worst_yield > gate.yield_tolerance)
        {
            if (process_data_.safety && process_data_.safety->active(t)) { OGS_FATAL("HM safety gate rejected equilibrium/admissibility."); }
            OGS_FATAL("HM initialisation gate rejected equilibrium/admissibility before working phases.");
        }
        INFO("HM initialisation gate passed before working phases.");
    }
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::
    setInitialConditionsConcreteProcess(std::vector<GlobalVector*> const& x,
                                        double const t,
                                        int const process_id)
{
    // So far, this function only sets the initial stress using the input data.
    if (process_id != process_data_.mechanics_related_process_id)
    {
        return;
    }

    DBUG("Set initial conditions of HydroMechanicsProcess.");
    postAcceptedTimestepConcreteProcess(x, t, process_id);

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::setInitialConditions, local_assemblers_,
        getActiveElementIDs(), getDOFTables(x.size()), x, t, process_id);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::computeSecondaryVariableConcrete(
    double const t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& x_prev, const int process_id)
{
    if (process_id != process_data_.hydraulic_process_id)
    {
        return;
    }

    DBUG("Compute the secondary variables for HydroMechanicsProcess.");
    if (process_data_.displacement_relative)
    {
        for (auto const* node : _mesh.getNodes())
        {
            for (int c = 0; c < DisplacementDim; ++c)
            {
                auto const index = getDOFTable(process_id).getGlobalIndex(
                    {_mesh.getID(), MeshLib::MeshItemType::Node, node->getID()}, 1, c);
                if (index < 0 || index == NumLib::MeshComponentMap::nop)
                { OGS_FATAL("Displacement references require complete serial displacement DOFs."); }
                auto const slot = node->getID() * DisplacementDim + c;
                (*process_data_.displacement_relative)[slot] = x[process_id]->get(index) -
                    process_data_.displacement_reference_values[slot];
            }
        }
    }


    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerIF::computeSecondaryVariable, local_assemblers_,
        getActiveElementIDs(), getDOFTables(x.size()), t, dt, x, x_prev,
        process_id);
}

template <int DisplacementDim>
std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
HydroMechanicsProcess<DisplacementDim>::getDOFTableForExtrapolatorData() const
{
    const bool manage_storage = false;
    return std::make_tuple(local_to_global_index_map_single_component_.get(),
                           manage_storage);
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
HydroMechanicsProcess<DisplacementDim>::getDOFTable(const int process_id) const
{
    if (hasMechanicalProcess(process_id))
    {
        return *_local_to_global_index_map;
    }

    // For the equation of pressure
    return *local_to_global_index_map_with_base_nodes_;
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::rollbackRejectedTimestep(int const process_id)
{
    if (!process_data_.adaptive_rollback) { return; }
    if (adaptive_ip_snapshot_.size()!=local_assemblers_.size()) { OGS_FATAL("Missing accepted adaptive snapshot."); }
    for (std::size_t i=0;i<local_assemblers_.size();++i) { local_assemblers_[i]->validateIntegrationPointCheckpoint(adaptive_ip_snapshot_[i]); }
    _boundary_conditions[process_id].validateCheckpointState(adaptive_bc_snapshot_,adaptive_snapshot_time_);
    for (std::size_t i=0;i<local_assemblers_.size();++i) { local_assemblers_[i]->restoreIntegrationPointCheckpoint(adaptive_ip_snapshot_[i]); }
    _boundary_conditions[process_id].restoreCheckpointState(adaptive_bc_snapshot_,adaptive_snapshot_time_);
    for (std::size_t i=0;i<local_assemblers_.size();++i)
    { if (local_assemblers_[i]->captureIntegrationPointCheckpoint()!=adaptive_ip_snapshot_[i]) { OGS_FATAL("Adaptive integration-point rollback is not byte-identical."); } }
    if (_boundary_conditions[process_id].checkpointState()!=adaptive_bc_snapshot_) { OGS_FATAL("Adaptive boundary rollback is not exact."); }
    INFO("HM restored accepted integration-point/MFront and boundary state after rejected global increment (byte-identical).");
}

template <int DisplacementDim>
std::vector<double> HydroMechanicsProcess<DisplacementDim>::excavationReferenceResidual(
    GlobalVector const& reference,double const t,double const dt,int const process_id)
{
    if(!sourceTermsSupportCheckpoint(process_id))
    { OGS_FATAL("Ramped excavation requires transactional structural source histories."); }
    GlobalVector copy{reference};std::vector<GlobalVector*> current{&copy};
    auto const source_snapshot=sourceCheckpointState(process_id);
    auto const boundary_snapshot=_boundary_conditions[process_id].checkpointState();
    std::vector<std::vector<std::string>> snapshot;
    for(auto const& assembler:local_assemblers_)
    { snapshot.push_back(assembler->captureIntegrationPointCheckpoint()); }
    auto const specification=getMatrixSpecifications(process_id);std::size_t matrix_id=0;
    auto& jac=NumLib::GlobalMatrixProvider::provider.getMatrix(specification,matrix_id);
    auto& rhs=NumLib::GlobalVectorProvider::provider.getVector(specification);
    jac.setZero();rhs.setZero();std::vector<double> result(reference.size());
    auto restore=[&]()
    {
        process_data_.excavation_reference_assembly=false;
        for(std::size_t i=0;i<snapshot.size();++i)
        { local_assemblers_[i]->restoreIntegrationPointCheckpoint(snapshot[i]); }
        _boundary_conditions[process_id].restoreCheckpointState(boundary_snapshot,excavation_accepted_time_);
        restoreSourceCheckpointState(process_id,source_snapshot);
        for(std::size_t i=0;i<snapshot.size();++i)
        { if(local_assemblers_[i]->captureIntegrationPointCheckpoint()!=snapshot[i])
          { OGS_FATAL("Excavation reference assembly changed accepted integration-point history."); } }
        if(sourceCheckpointState(process_id)!=source_snapshot)
        { OGS_FATAL("Excavation reference assembly changed accepted structural history."); }
        NumLib::GlobalMatrixProvider::provider.releaseMatrix(jac);
        NumLib::GlobalVectorProvider::provider.releaseVector(rhs);
    };
    process_data_.excavation_reference_assembly=true;
    try
    {
        this->assembleWithJacobian(t,dt,current,current,process_id,rhs,jac);
        MathLib::finalizeVectorAssembly(rhs);
        for(GlobalIndexType i=0;i<reference.size();++i) { result[i]=-rhs.get(i); }
    }
    catch(...) { restore();throw; }
    restore();return result;
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::prepareExcavationRelease(
    std::vector<GlobalVector*> const& x,double const t,double const dt,int const process_id)
{
    if(process_data_.excavation_release.empty()) { return; }
    std::size_t phase=0;
    for(auto const& interval:process_data_.excavation_release)
    {
        double const event=interval.start+1e-8*std::max(1.,std::abs(interval.start)/2.);
        if(t>=event && t<=interval.end+1e-10)
        {
            if(excavation_release_captured_==phase) { return; }
            if(std::abs(excavation_accepted_time_-interval.start)>1e-9 ||
               excavation_accepted_solution_.size()!=static_cast<std::size_t>(x[process_id]->size()) ||
               excavation_accepted_residual_.size()!=excavation_accepted_solution_.size())
            { OGS_FATAL("Excavation release requires the accepted phase-start solution and residual."); }
            GlobalVector reference{*x[process_id]};
            for(GlobalIndexType i=0;i<reference.size();++i)
            { reference.set(i,excavation_accepted_solution_[i]); }
            MathLib::LinAlg::finalizeAssembly(reference);
            auto const current=excavationReferenceResidual(reference,t,dt,process_id);
            std::set<GlobalIndexType> essential;
            for(auto const& bc:*getKnownSolutions(t,reference,process_id))
            { essential.insert(bc.ids.begin(),bc.ids.end()); }
            excavation_release_values_.assign(reference.size(),0.);
            double norm2=0.,old_norm2=0.;
            for(auto const* node:_mesh.getNodes())
            {
                for(int c=0;c<DisplacementDim;++c)
                {
                    auto const index=getDOFTable(process_id).getGlobalIndex(
                        {_mesh.getID(),MeshLib::MeshItemType::Node,node->getID()},1,c);
                    if(index<0 || index==NumLib::MeshComponentMap::nop)
                    { OGS_FATAL("Excavation release requires complete serial mechanical DOFs."); }
                    if(essential.contains(index)) { continue; }
                    double const force=current[index]-excavation_accepted_residual_[index];
                    excavation_release_values_[index]=force;norm2+=force*force;
                    old_norm2+=excavation_accepted_residual_[index]*excavation_accepted_residual_[index];
                    if(std::abs(force)>1e-10)
                    { INFO("HM excavation release node {} component {}: start={}, end={}, force={}",node->getID(),c,interval.start,interval.end,force); }
                }
            }
            excavation_release_captured_=phase;
            INFO("HM excavation release captured accepted t={}: change_norm={} N/m, prior_free_residual={} N/m; hydraulic residual is not compensated.",interval.start,std::sqrt(norm2),std::sqrt(old_norm2));
            return;
        }
        ++phase;
    }
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::applyExcavationRelease(double const t,GlobalVector& rhs)
{
    if(process_data_.excavation_reference_assembly || process_data_.excavation_release.empty()) { return; }
    auto* const output=process_data_.excavation_release_force;
    std::fill(output->begin(),output->end(),0.);
    if(excavation_release_captured_>=process_data_.excavation_release.size()) { return; }
    auto const& interval=process_data_.excavation_release[excavation_release_captured_];
    double const fraction=t>interval.start?interval.remaining(t):0.;
    if(fraction==0.) { return; }
    for(auto const* node:_mesh.getNodes())
    {
        for(int c=0;c<DisplacementDim;++c)
        {
            auto const index=getDOFTable(0).getGlobalIndex(
                {_mesh.getID(),MeshLib::MeshItemType::Node,node->getID()},1,c);
            double const force=fraction*excavation_release_values_.at(index);
            rhs.add(index,force);(*output)[node->getID()*DisplacementDim+c]=force;
        }
    }
}

template class HydroMechanicsProcess<2>;
template class HydroMechanicsProcess<3>;

}  // namespace HydroMechanics
}  // namespace ProcessLib
