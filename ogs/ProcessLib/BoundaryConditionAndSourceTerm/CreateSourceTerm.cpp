// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateSourceTerm.h"

#include <cassert>
#include <numeric>

#include "CreateAnchorTerm.h"
#include "Plate2D.h"
#include "Interface2D.h"
#include "EmbeddedPile2D.h"
#include "ParameterLib/Utils.h"
#include "CreateEmbeddedAnchor.h"
#include "CreateNodalSourceTerm.h"
#include "CreateVolumetricSourceTerm.h"
#include "Python/CreatePythonSourceTerm.h"
#include "SourceTerm.h"
#include "SourceTermConfig.h"
#include "ProcessLib/ProcessVariable.h"

namespace ProcessLib
{
std::unique_ptr<SourceTermBase> createSourceTerm(
    const SourceTermConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table_bulk,
    const MeshLib::Mesh& source_term_mesh, const int variable_id,
    const unsigned integration_order, const unsigned shapefunction_order,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    [[maybe_unused]] std::vector<std::reference_wrapper<ProcessVariable>> const&
        all_process_variables_for_this_process,
    const MeshLib::Mesh& bulk_mesh)
{
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
    auto const type = config.config.peekConfigParameter<std::string>("type");

    // check basic data consistency
    if (variable_id >=
            static_cast<int>(dof_table_bulk.getNumberOfVariables()) ||
        config.component_id >=
            dof_table_bulk.getNumberOfVariableComponents(variable_id))
    {
        OGS_FATAL(
            "Variable id or component id too high. Actual values: ({:d}, "
            "{:d}), maximum values: ({:d}, {:d}).",
            variable_id, config.component_id,
            dof_table_bulk.getNumberOfVariables(),
            dof_table_bulk.getNumberOfVariableComponents(variable_id));
    }

    if (!source_term_mesh.getProperties()
             .template existsPropertyVector<std::size_t>(
                 MeshLib::getBulkIDString(MeshLib::MeshItemType::Node)) &&
        (type != "EmbeddedAnchor" && type != "EmbeddedPile2D"))
    {
        OGS_FATAL(
            "The required bulk node ids map does not exist in the source term "
            "mesh '{:s}'.",
            source_term_mesh.getName());
    }
    std::vector<MeshLib::Node*> const& source_term_nodes =
        source_term_mesh.getNodes();
    DBUG(
        "Found {:d} nodes for source term at mesh '{:s}' for the variable {:d} "
        "and component {:d}",
        source_term_nodes.size(), source_term_mesh.getName(), variable_id,
        config.component_id);

    MeshLib::MeshSubset source_term_mesh_subset(source_term_mesh,
                                                source_term_nodes);

    if (type == "Nodal")
    {
        auto dof_table_source_term =
            dof_table_bulk.deriveBoundaryConstrainedMap(
                variable_id, {config.component_id},
                std::move(source_term_mesh_subset));
        return ProcessLib::createNodalSourceTerm(
            config.config, config.mesh, std::move(dof_table_source_term),
            source_term_mesh.getID(), variable_id, config.component_id,
            parameters);
    }
    if(type=="EmbeddedPile2D")
    {
        config.config.checkConfigParameter("type","EmbeddedPile2D");
        if(config.component_id!=0 || dof_table_bulk.getNumberOfVariableComponents(variable_id)!=2 || bulk_mesh.getDimension()!=2)
        { OGS_FATAL("EmbeddedPile2D requires component0 on planar displacement[2]."); }
        auto const reduction=config.config.getConfigParameter<std::string>("resistance_reduction_policy","Unspecified");
        if(reduction!="Unspecified" && reduction!="Fixed" && reduction!="ReduceWithSoil")
        { OGS_FATAL("Embedded resistance_reduction_policy must be Fixed or ReduceWithSoil."); }
        auto parameter=[&](char const* key)->ParameterLib::Parameter<double> const*
        {
            auto const name=config.config.getConfigParameterOptional<std::string>(key);
            return name?&ParameterLib::findParameter<double>(*name,parameters,1,&source_term_mesh):nullptr;
        };
        auto const connection=config.config.getConfigParameter<std::string>("head_connection","Independent");
        if(connection!="Independent" && connection!="BondedToSoil" && connection!="BondedToBulkNode")
        { OGS_FATAL("Embedded head_connection must be Independent, BondedToSoil or BondedToBulkNode."); }
        auto const head_bulk=config.config.getConfigParameterOptional<std::size_t>("head_bulk_node_id");
        if((connection=="BondedToBulkNode")!=static_cast<bool>(head_bulk))
        { OGS_FATAL("head_bulk_node_id is required only for BondedToBulkNode."); }
        bool const interaction=config.config.getConfigParameter<bool>("use_interaction_factors",false);
        bool const stiffness=config.config.getConfigParameter<bool>("use_member_stiffness",false);
        StructuralElements::PileProperties props{
            config.config.getConfigParameter<double>("axial_stiffness"),config.config.getConfigParameter<double>("bending_stiffness"),
            config.config.getConfigParameter<double>("shaft_stiffness"),config.config.getConfigParameter<double>("shaft_resistance"),
            config.config.getConfigParameter<double>("lateral_stiffness"),config.config.getConfigParameter<double>("lateral_resistance",std::numeric_limits<double>::infinity()),
            config.config.getConfigParameter<double>("base_stiffness",0.),config.config.getConfigParameter<double>("base_resistance",0.)};
        return std::make_unique<EmbeddedPile2D>(bulk_mesh,source_term_mesh,dof_table_bulk,variable_id,props,
            config.config.getConfigParameter<std::size_t>("head_node",0),config.config.getConfigParameter<std::size_t>("tip_node",source_term_mesh.getNumberOfNodes()-1),
            config.config.getConfigParameter<double>("spacing",1.),config.config.getConfigParameter<double>("activation_time",0.),reduction,
            parameter("head_axial_force_parameter"),parameter("head_lateral_force_parameter"),
            parameter("head_axial_displacement_parameter"),parameter("head_lateral_displacement_parameter"),connection!="Independent",interaction,stiffness,
            config.config.getConfigParameter<double>("prestress_force",0.),
            config.config.getConfigParameter<double>("prestress_end_time",-1.),
            config.config.getConfigParameter<bool>("prestress_ramp",false),head_bulk.value_or(std::numeric_limits<std::size_t>::max()));
    }
    if(type=="Interface2D")
    {
        config.config.checkConfigParameter("type","Interface2D");
        if(config.component_id!=0 || dof_table_bulk.getNumberOfVariableComponents(variable_id)!=2 || bulk_mesh.getDimension()!=2)
        { OGS_FATAL("Interface2D requires component0 on planar displacement[2]."); }
        auto map=dof_table_bulk.deriveBoundaryConstrainedMap(variable_id,{0,1},std::move(source_term_mesh_subset));
        auto const normal=config.config.getConfigParameter<std::vector<double>>("normal");
        if(normal.size()!=2) { OGS_FATAL("Interface2D normal requires two values."); }
        auto const reduction=config.config.getConfigParameter<std::string>("strength_reduction_policy","Unspecified");
        if(reduction!="Unspecified" && reduction!="Fixed" && reduction!="ReduceWithSoil")
        { OGS_FATAL("Interface strength_reduction_policy must be Fixed or ReduceWithSoil."); }
        auto const opening=config.config.getConfigParameter<std::string>("opening_shear_policy","ReleaseAllShear");
        if (opening!="ReleaseAllShear" && opening!="RetainCohesion")
        { OGS_FATAL("Interface opening_shear_policy must be ReleaseAllShear or RetainCohesion."); }
        int pressure_variable=-1;
        if(config.config.getConfigParameter<bool>("use_hydraulic_pressure",false))
        {
            if(all_process_variables_for_this_process.size()!=2 || variable_id!=1 ||
               dof_table_bulk.getNumberOfVariableComponents(0)!=1)
            { OGS_FATAL("Wet Interface2D requires monolithic HM with pressure first and displacement second."); }
            pressure_variable=0;
        }
        auto const preload_name=config.config.getConfigParameterOptional<std::string>("initial_normal_pressure_parameter");
        ParameterLib::Parameter<double> const* preload=preload_name
            ? &ParameterLib::findParameter<double>(*preload_name,parameters,1) : nullptr;
        return std::make_unique<Interface2D>(std::move(map),source_term_mesh,Eigen::Vector2d(normal[0],normal[1]),
            config.config.getConfigParameter<double>("normal_stiffness"),config.config.getConfigParameter<double>("shear_stiffness"),
            config.config.getConfigParameter<double>("cohesion"),config.config.getConfigParameter<double>("friction_coefficient"),reduction,
            bulk_mesh,dof_table_bulk,pressure_variable,config.config.getConfigParameter<unsigned>("soil_side",0),preload,opening=="RetainCohesion");
    }
    if (type == "Plate2D")
    {
        config.config.checkConfigParameter("type", "Plate2D");
        if (config.component_id!=0 ||
            dof_table_bulk.getNumberOfVariableComponents(variable_id)!=2 ||
            bulk_mesh.getDimension()!=2)
        { OGS_FATAL("Plate2D requires planar displacement with two components."); }
        auto map=dof_table_bulk.deriveBoundaryConstrainedMap(
            variable_id,{0,1},std::move(source_term_mesh_subset));
        double const EA=config.config.getConfigParameter<double>("axial_stiffness");
        double const EI=config.config.getConfigParameter<double>("bending_stiffness");
        double const GA=config.config.getConfigParameter<double>("shear_stiffness",std::numeric_limits<double>::infinity());
        double Np=0.,Mp=0.;
        auto const law=config.config.getConfigParameterOptional<std::string>("section_law");
        if(law)
        {
            if(*law!="DiamondPerfectPlastic") { OGS_FATAL("Unknown plate section_law {}",*law); }
            Np=config.config.getConfigParameter<double>("axial_capacity");
            Mp=config.config.getConfigParameter<double>("moment_capacity");
            if(!std::isfinite(Np) || Np<=0. || !std::isfinite(Mp) || Mp<=0. || std::isfinite(GA))
            { OGS_FATAL("Plastic plates need positive capacities and omitted shear_stiffness."); }
        }
        auto const distributed=config.config.getConfigParameter<std::vector<double>>("distributed_load",std::vector<double>{0.,0.});
        if(distributed.size()!=2) { OGS_FATAL("Plate2D distributed_load requires two global components."); }
        std::optional<StructuralElements::DistributedLoadRamp> load_ramp;
        if(auto ramp=config.config.getConfigSubtreeOptional("distributed_load_ramp"))
        {
            load_ramp.emplace(ramp->getConfigParameter<double>("start_time"),
                              ramp->getConfigParameter<double>("end_time"));
        }
        return std::make_unique<Plate2D>(std::move(map),source_term_mesh,EA,EI,GA,Np,Mp,Eigen::Vector2d(distributed[0],distributed[1]),load_ramp);
    }
    if (type == "Anchor")
    {
        const int number_of_components =
            dof_table_bulk.getNumberOfVariableComponents(variable_id);
        std::vector<int> component_ids(number_of_components);
        std::iota(std::begin(component_ids), std::end(component_ids), 0);
        auto dof_table_source_term =
            dof_table_bulk.deriveBoundaryConstrainedMap(
                variable_id, component_ids, std::move(source_term_mesh_subset));
        const int bulk_mesh_dimension =
            dof_table_bulk.getMeshSubset(variable_id, 0)
                .getMesh()
                .getDimension();
        if (bulk_mesh_dimension != number_of_components)
        {
            OGS_FATAL(
                "For the Anchor source term type,"
                "the bulk mesh dimension needs to be the same "
                "as the number of process variable components.");
        }
        switch (bulk_mesh_dimension)
        {
            case 2:
                return ProcessLib::createAnchorTerm<2>(
                    config.config, config.mesh,
                    std::move(dof_table_source_term), source_term_mesh.getID(),
                    variable_id, parameters);
            case 3:
                return ProcessLib::createAnchorTerm<3>(
                    config.config, config.mesh,
                    std::move(dof_table_source_term), source_term_mesh.getID(),
                    variable_id, parameters);
            default:
                OGS_FATAL(
                    "Anchor can not be instantiated "
                    "for mesh dimensions other than two or three. "
                    "{}-dimensional mesh was given.",
                    bulk_mesh_dimension);
        }
    }
    if (type == "EmbeddedAnchor")
    {
        const int number_of_components =
            dof_table_bulk.getNumberOfVariableComponents(variable_id);

        const int bulk_mesh_dimension = bulk_mesh.getDimension();
        if (bulk_mesh_dimension != number_of_components)
        {
            OGS_FATAL(
                "For the EmbeddedAnchor source term type,"
                "the bulk mesh dimension needs to be the same "
                "as the number of process variable components.");
        }
        switch (bulk_mesh_dimension)
        {
            case 2:
                return ProcessLib::createEmbeddedAnchor<2>(
                    config.config, config.mesh, bulk_mesh, dof_table_bulk,
                    source_term_mesh.getID(), variable_id);
            case 3:
                return ProcessLib::createEmbeddedAnchor<3>(
                    config.config, config.mesh, bulk_mesh, dof_table_bulk,
                    source_term_mesh.getID(), variable_id);
            default:
                OGS_FATAL(
                    "Anchor can not be instantiated "
                    "for mesh dimensions other than two or three. "
                    "{}-dimensional mesh was given.",
                    bulk_mesh_dimension);
        }
    }
    if (type == "Line" || type == "Volumetric")
    {
        auto dof_table_source_term =
            dof_table_bulk.deriveBoundaryConstrainedMap(
                variable_id, {config.component_id},
                std::move(source_term_mesh_subset));
        auto const& bulk_mesh_dimension =
            dof_table_bulk.getMeshSubset(variable_id, config.component_id)
                .getMesh()
                .getDimension();
        return ProcessLib::createVolumetricSourceTerm(
            config.config, bulk_mesh_dimension, config.mesh,
            std::move(dof_table_source_term), parameters, integration_order,
            shapefunction_order);
    }
    if (type == "Python")
    {
        auto dof_table_source_term =
            dof_table_bulk.deriveBoundaryConstrainedMap(
                std::move(source_term_mesh_subset));

        return ProcessLib::createPythonSourceTerm(
            config.config, config.mesh, std::move(dof_table_source_term),
            variable_id, config.component_id, integration_order,
            shapefunction_order, source_term_mesh.getDimension(),
            all_process_variables_for_this_process);
    }

    OGS_FATAL("Unknown source term type: `{:s}'.", type);
}
}  // namespace ProcessLib
