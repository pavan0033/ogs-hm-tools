// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateHydroMechanicsProcess.h"

#include <cassert>
#include <range/v3/range.hpp>
#include <string>

#include "HydroMechanicsProcess.h"
#include "HydroMechanicsProcessData.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/Constant.h"
#include "MaterialLib/MPL/Properties/Linear.h"
#include "MaterialLib/MPL/Utils/CheckMPLPhasesForSinglePhaseFlow.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MeshLib/Utils/Is2DMeshOnRotatedVerticalPlane.h"
#include "ParameterLib/Utils.h"
#include "ParameterLib/FunctionParameter.h"
#include "ProcessLib/Common/HydroMechanics/CreateInitialStress.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
namespace HydroMechanics
{

CouplingScheme parseCouplingScheme(
    std::optional<BaseLib::ConfigTree> const& config)
{
    if (!config)
    {
        return Monolithic{};
    }

    auto const coupling_scheme_type =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__coupling_scheme__type}
        config->getConfigParameter<std::string>("type");

    if (coupling_scheme_type == "monolithic")
    {
        return Monolithic{};
    }

    // Default value is 0.5, which is recommended by [Mikelic & Wheeler].
    double const fixed_stress_stabilization_parameter =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__coupling_scheme__fixed_stress_stabilization_parameter}
        config->getConfigParameter<double>(
            "fixed_stress_stabilization_parameter", 0.5);

    DBUG("Using value {:g} for coupling parameter of staggered scheme.",
         fixed_stress_stabilization_parameter);

    {  // Check parameter value. Optimum is not a-priori known, but within
       // certain interval [Storvik & Nordbotten]
        double const csp_min = 1.0 / 6.0;
        double const csp_max = 1.0;
        if (fixed_stress_stabilization_parameter < csp_min ||
            fixed_stress_stabilization_parameter > csp_max)
        {
            WARN(
                "Value of coupling scheme parameter = {:g} is out of "
                "reasonable range ({:g}, {:g}).",
                fixed_stress_stabilization_parameter, csp_min, csp_max);
        }
    }

    bool const fixed_stress_over_time_step =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__coupling_scheme__fixed_stress_over_time_step}
        config->getConfigParameter<std::string>("fixed_stress_over_time_step",
                                                "false") == "true";

    return Staggered{fixed_stress_stabilization_parameter,
                     fixed_stress_over_time_step};
}

template <int DisplacementDim>
std::unique_ptr<Process> createHydroMechanicsProcess(
    std::string const& name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order, BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "HYDRO_MECHANICS");
    DBUG("Create HydroMechanicsProcess.");

    if (DisplacementDim == 2)
    {
        if (mesh.isAxiallySymmetric() &&
            MeshLib::is2DMeshOnRotatedVerticalPlane(mesh))
        {
            OGS_FATAL(
                "Mesh {:s} is on a plane rotated around the vertical axis. The "
                "axisymmetric problem can not use such mesh.",
                mesh.getName());
        }
    }

    auto const coupling_scheme = parseCouplingScheme(
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__coupling_scheme}
        config.getConfigSubtreeOptional("coupling_scheme"));

    /// \section processvariableshm Process Variables

    //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_p = nullptr;
    ProcessVariable* variable_u = nullptr;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;

    int const hydraulic_process_id = 0;
    int mechanics_related_process_id = 0;

    if (std::holds_alternative<Monolithic>(coupling_scheme))
    {
        /// Primary process variables as they appear in the global component
        /// vector:
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__process_variables__pressure}
             "pressure",
             //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__process_variables__displacement}
             "displacement"});
        variable_p = &per_process_variables[0].get();
        variable_u = &per_process_variables[1].get();
        process_variables.push_back(std::move(per_process_variables));
    }

    if (std::holds_alternative<Staggered>(coupling_scheme))
    {
        using namespace std::string_literals;
        for (auto const& variable_name : {"pressure"s, "displacement"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        mechanics_related_process_id = 1;
        variable_p = &process_variables[hydraulic_process_id][0].get();
        variable_u = &process_variables[mechanics_related_process_id][0].get();
    }

    DBUG("Associate displacement with process variable '{:s}'.",
         variable_u->getName());

    if (variable_u->getNumberOfGlobalComponents() != DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '{:s}' is different "
            "from the displacement dimension: got {:d}, expected {:d}",
            variable_u->getName(),
            variable_u->getNumberOfGlobalComponents(),
            DisplacementDim);
    }

    DBUG("Associate pressure with process variable '{:s}'.",
         variable_p->getName());
    if (variable_p->getNumberOfGlobalComponents() != 1)
    {
        OGS_FATAL(
            "Pressure process variable '{:s}' is not a scalar variable but has "
            "{:d} components.",
            variable_p->getName(),
            variable_p->getNumberOfGlobalComponents());
    }

    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
            parameters, local_coordinate_system, materialIDs(mesh), config);

    /// \section parametershm Process Parameters
    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (b.size() != DisplacementDim)
        {
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is {:d}, displacement "
                "dimension is {:d}",
                b.size(), DisplacementDim);
        }

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__mass_lumping}
    auto mass_lumping = config.getConfigParameter<bool>("mass_lumping", false);

    auto const is_linear =
        //! \ogs_file_param{prj__processes__process__linear}
        config.getConfigParameter("linear", false);

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    std::array const requiredMediumProperties = {
        MaterialPropertyLib::reference_temperature,
        MaterialPropertyLib::permeability, MaterialPropertyLib::porosity,
        MaterialPropertyLib::biot_coefficient};
    std::array const requiredFluidProperties = {MaterialPropertyLib::viscosity,
                                                MaterialPropertyLib::density};
    std::array const requiredSolidProperties = {MaterialPropertyLib::density};

    MaterialPropertyLib::checkMPLPhasesForSinglePhaseFlow(mesh, media_map);

    // The uniqueness of phase has already been checked in
    // `checkMPLPhasesForSinglePhaseFlow`.
    MaterialPropertyLib::Variable const phase_variable =
        (*ranges::begin(media_map.media()))
                ->hasPhase(MaterialPropertyLib::PhaseName::Gas)
            ? MaterialPropertyLib::Variable::gas_phase_pressure
            : MaterialPropertyLib::Variable::liquid_phase_pressure;
    for (auto const& medium : media_map.media())
    {
        checkRequiredProperties(*medium, requiredMediumProperties);
        checkRequiredProperties(fluidPhase(*medium), requiredFluidProperties);
        checkRequiredProperties(
            medium->phase(MaterialPropertyLib::PhaseName::Solid),
            requiredSolidProperties);
    }
    DBUG("Media properties verified.");

    // Initial stress conditions
    auto initial_stress = ProcessLib::createInitialStress<DisplacementDim>(
        config, parameters, mesh);

    const bool use_taylor_hood_elements = variable_p->getShapeFunctionOrder() !=
                                          variable_u->getShapeFunctionOrder();

    HydroMechanicsProcessData<DisplacementDim> process_data{
        materialIDs(mesh),
        std::move(media_map),
        std::move(solid_constitutive_relations),
        initial_stress,
        specific_body_force,
        mass_lumping,
        coupling_scheme,
        hydraulic_process_id,
        mechanics_related_process_id,
        use_taylor_hood_elements,
        phase_variable};

    if (auto const safety=config.getConfigSubtreeOptional("native_safety"))
    {
        if (DisplacementDim!=2 || mesh.isAxiallySymmetric() ||
            !process_data.isMonolithicSchemeUsed() || is_linear)
        { OGS_FATAL("Safety requires nonlinear monolithic plane-strain HM."); }
        SafetyContinuation settings;
        settings.physical_time=safety->getConfigParameter<double>("physical_time");
        settings.initialisation_end=safety->getConfigParameter<double>("initialisation_end");
        settings.pressure_policy=safety->getConfigParameter<std::string>("pressure_policy");
        settings.working_hydraulic_rates=safety->getConfigParameter<bool>("working_hydraulic_rates",true);
        settings.verify_equilibrium=safety->getConfigParameter<bool>("verify_equilibrium",false);
        for (auto const& point:safety->getConfigSubtreeList("point"))
        {
            settings.coordinates.push_back(point.getConfigParameter<double>("coordinate"));
            settings.factors.push_back(point.getConfigParameter<double>("factor"));
        }
        if (settings.coordinates.size()<2 || settings.coordinates.front()!=settings.physical_time)
        { OGS_FATAL("Invalid native safety continuation coordinates."); }
        process_data.safety=std::move(settings);
    }

    if (auto const checkpoint = config.getConfigSubtreeOptional("native_checkpoint"))
    {
        if constexpr (DisplacementDim != 2)
        { OGS_FATAL("Native HM checkpoint first scope requires plane strain."); }
        if (!process_data.isMonolithicSchemeUsed())
        { OGS_FATAL("Native HM checkpoint requires monolithic HM."); }
        NativeCheckpointSettings settings;
        settings.initialisation_time=checkpoint->getConfigParameter<double>("initialisation_time");
        settings.restart_file=checkpoint->getConfigParameter<std::string>("restart_file", "");
        settings.restart_identity=checkpoint->getConfigParameter<std::string>("restart_identity", "");
        for (auto const& entry:checkpoint->getConfigSubtreeList("point"))
        {
            settings.points.push_back({entry.getConfigParameter<std::string>("id"),
                entry.getConfigParameter<std::string>("file"),
                entry.getConfigParameter<std::string>("model_identity"),
                entry.getConfigParameter<double>("time")});
        }
        process_data.native_checkpoint=std::move(settings);
    }

    if (auto births=config.getConfigSubtreeOptional("construction_births"))
    {
        for (auto const& birth:births->getConfigSubtreeList("birth"))
        {
            process_data.construction_births.push_back({birth.getConfigParameter<int>("material_id"),
                birth.getConfigParameter<double>("time"),birth.getConfigParameter<bool>("reset")});
        }
    }
    if (auto const seepage = config.getConfigSubtreeOptional("steady_state_seepage"))
    {
        if (DisplacementDim != 2 || mesh.isAxiallySymmetric() ||
            !process_data.isMonolithicSchemeUsed() || is_linear ||
            process_data.mechanical_drainage)
        { OGS_FATAL("SteadyStateSeepage requires nonlinear monolithic plane-strain ordinary HM."); }
        double previous_end = -1.;
        for (auto const& interval : seepage->getConfigSubtreeList("interval"))
        {
            double const start = interval.getConfigParameter<double>("start");
            double const end = interval.getConfigParameter<double>("end");
            if (!std::isfinite(start) || !std::isfinite(end) ||
                start < 0. || end <= start || start < previous_end)
            { OGS_FATAL("Invalid SteadyStateSeepage interval."); }
            process_data.steady_state_seepage.intervals.emplace_back(start, end);
            previous_end = end;
        }
    }
    if (auto const seepage = config.getConfigSubtreeOptional("transient_seepage"))
    {
        if (DisplacementDim != 2 || mesh.isAxiallySymmetric() ||
            !process_data.isMonolithicSchemeUsed() || is_linear || process_data.mechanical_drainage)
        { OGS_FATAL("TransientSeepage requires nonlinear monolithic plane-strain ordinary HM."); }
        double previous_end = -1.;
        for (auto const& interval : seepage->getConfigSubtreeList("interval"))
        {
            double const start = interval.getConfigParameter<double>("start");
            double const end = interval.getConfigParameter<double>("end");
            if (!std::isfinite(start) || !std::isfinite(end) || start < 0. ||
                end <= start || start < previous_end)
            { OGS_FATAL("Invalid TransientSeepage interval."); }
            process_data.transient_seepage.intervals.emplace_back(start, end);
            previous_end = end;
        }
    }
    if(auto const release=config.getConfigSubtreeOptional("excavation_release"))
    {
#ifdef USE_PETSC
        OGS_FATAL("Ramped excavation first scope requires serial Eigen.");
#endif
        if(DisplacementDim!=2 || mesh.isAxiallySymmetric() ||
           !process_data.isMonolithicSchemeUsed() || is_linear ||
           process_data.mechanical_drainage || process_data.safety ||
           !process_data.construction_births.empty())
        { OGS_FATAL("Ramped excavation requires ordinary nonlinear monolithic plane-strain MC without Safety/construction."); }
        double previous=-1.;
        for(auto const& interval:release->getConfigSubtreeList("interval"))
        {
            double const start=interval.getConfigParameter<double>("start");
            double const end=interval.getConfigParameter<double>("end");
            if(!std::isfinite(start) || !std::isfinite(end) || start<=0. || end<=start || start<previous)
            { OGS_FATAL("Invalid excavation release interval."); }
            process_data.excavation_release.push_back({start,end});previous=end;
        }
    }
    process_data.adaptive_rollback=config.getConfigParameter<bool>("adaptive_rollback",false);
    process_data.initialisation_load_ramp_end =
        config.getConfigParameter<double>("initialisation_load_ramp_end", 0.);
    if (!std::isfinite(process_data.initialisation_load_ramp_end) ||
        process_data.initialisation_load_ramp_end < 0. ||
        (process_data.initialisation_load_ramp_end > 0. &&
         !process_data.isMonolithicSchemeUsed()))
    {
        OGS_FATAL("Synchronous initialisation currently requires monolithic HM and finite nonnegative ramp end.");
    }

    if (auto const references = config.getConfigSubtreeOptional("displacement_reference_times"))
    {
        if (!process_data.isMonolithicSchemeUsed())
        { OGS_FATAL("Displacement references currently require monolithic HM."); }
        for (auto const& entry : references->getConfigSubtreeList("time"))
        {
            double const time = entry.getValue<double>();
            if (!std::isfinite(time) || time <= 0. ||
                (!process_data.displacement_reference_times.empty() &&
                 time <= process_data.displacement_reference_times.back()))
            { OGS_FATAL("Displacement reference times must be positive, finite and increasing."); }
            process_data.displacement_reference_times.push_back(time);
        }
    }

    if (auto const gate = config.getConfigSubtreeOptional("initialisation_gate"))
    {
        if (DisplacementDim != 2 || mesh.isAxiallySymmetric() ||
            !process_data.isMonolithicSchemeUsed() || is_linear)
        { OGS_FATAL("Initialisation gate currently requires nonlinear monolithic plane-strain 2D HM."); }
        InitialisationGate settings{gate->getConfigParameter<double>("reference_time"),
            gate->getConfigParameter<double>("force_tolerance"),
            gate->getConfigParameter<double>("yield_tolerance"), {}};
        for (auto const material : gate->getConfigSubtreeList("material"))
        {
            int const id = material.getConfigAttribute<int>("id");
            if (auto const name = material.getConfigParameterOptional<std::string>("cohesion_parameter"))
            {
                auto const& parameter = ParameterLib::findParameter<double>(*name, parameters, 1);
                auto const* function = dynamic_cast<ParameterLib::FunctionParameter<double> const*>(&parameter);
                if (parameter.isTimeDependent() && (!function || function->usesTimeSymbol()))
                { OGS_FATAL("Safety cohesion must be time independent."); }
                settings.cohesion_parameters.emplace(id, &parameter);
            }
            settings.materials.emplace(id, MCAdmissibilityParameters{
                material.getConfigParameter<double>("cohesion"),
                material.getConfigParameter<double>("friction_angle"),
                material.getConfigParameter<double>("transition_angle"),
                material.getConfigParameter<double>("cutoff")});
        }
        for (auto const& [id, medium] : media)
        {
            if (!settings.materials.contains(id)) { OGS_FATAL("Initialisation gate requires exact constitutive/medium ID coverage."); }
        }
        if (settings.materials.size() != media.size()) { OGS_FATAL("Initialisation gate has uncovered constitutive IDs."); }
        process_data.initialisation_gate = std::move(settings);
    }

    auto const material_state_output =
        config.getConfigParameter<std::string>("material_state_output", "None");
    if (material_state_output != "None" &&
        material_state_output != "IntegrationPoints")
    {
        OGS_FATAL("material_state_output must be None or IntegrationPoints");
    }
    process_data.output_material_state = material_state_output == "IntegrationPoints";
    // A single source of fluid compressibility: the density EOS derivative.
    // Reject obsolete inputs rather than silently dropping their storage.
    for (auto const& [id, medium] : media)
    {
        if (fluidPhase(*medium).hasProperty(MaterialPropertyLib::PropertyType::bulk_modulus))
        { OGS_FATAL("HM fluid bulk_modulus input has been removed. Define pressure-dependent fluid density instead; beta_f=(1/rho)*d(rho)/d(p). Constant density means zero fluid compressibility."); }
    }
    auto densityDerivedModulus = [&](MaterialPropertyLib::Medium const& medium) -> std::optional<double>
    {
        namespace MPL = MaterialPropertyLib;
        auto const& density=fluidPhase(medium).property(MPL::PropertyType::density);
        if (!dynamic_cast<MPL::Constant const*>(&density) && !dynamic_cast<MPL::Linear const*>(&density))
        { OGS_FATAL("Strict mechanical storage currently supports Constant or pressure-Linear density; general EOS remains available to ordinary HM."); }
        MPL::VariableArray vars;
        vars.gas_phase_pressure=0.;
        vars.liquid_phase_pressure=0.;
        vars.temperature=medium.property(MPL::PropertyType::reference_temperature)
            .template value<double>(MPL::EmptyVariableArray,ParameterLib::SpatialPosition{},0.,0.);
        double const rho=density.template value<double>(vars,ParameterLib::SpatialPosition{},0.,0.);
        double const derivative=density.template dValue<double>(vars,phase_variable,ParameterLib::SpatialPosition{},0.,0.);
        if (!std::isfinite(rho) || rho<=0. || !std::isfinite(derivative) || derivative<0.)
        { OGS_FATAL("Density-derived mechanical storage needs positive density and finite nonnegative pressure derivative at reference pressure zero."); }
        return derivative>0.?std::optional<double>{rho/derivative}:std::nullopt;
    };
    if (auto const drainage = config.getConfigSubtreeOptional("mechanical_drainage"))
    {
        if (DisplacementDim != 2 || mesh.isAxiallySymmetric() ||
            !std::holds_alternative<Monolithic>(coupling_scheme) || is_linear ||
            variable_p->getShapeFunctionOrder() != 1 ||
            variable_u->getShapeFunctionOrder() != 2 || !process_data.output_material_state)
        {
            OGS_FATAL("Mechanical drainage currently requires nonlinear monolithic plane-strain P1/P2 HM with IntegrationPoints material-state output.");
        }
        for (auto const check : drainage->getConfigSubtreeList("static_parameter"))
        {
            auto const name = check.getConfigParameter<std::string>("parameter");
            auto const& parameter = ParameterLib::findParameter<double>(name, parameters, 1);
            auto const* function = dynamic_cast<ParameterLib::FunctionParameter<double> const*>(&parameter);
            if (!function || function->usesTimeSymbol())
            { OGS_FATAL("Mechanical drainage spatial stiffness/strength must be a static scalar Function."); }
        }
        MechanicalDrainageSettings settings{
            drainage->getConfigParameter<double>("reference_time"), {}};
        if (!std::isfinite(settings.reference_time) || settings.reference_time <= 0)
        {
            OGS_FATAL("Mechanical drainage requires positive finite reference time.");
        }
        for (auto const& entry : drainage->getConfigSubtreeList("material"))
        {
            int const id = entry.getConfigAttribute<int>("id");
            auto const type = entry.getConfigParameter<std::string>("type");
            auto const explicit_modulus = entry.getConfigParameterOptional<double>("biot_modulus");
            std::optional<double> fluid_modulus;
            if (id >= 0 && media.contains(id) && type != "UndrainedTotalStress")
            { fluid_modulus=densityDerivedModulus(*media.at(id)); }
            auto const pressure_coupling = entry.getConfigParameter<std::string>("pressure_coupling", "ClipActual");
            if (pressure_coupling != "ClipActual" && pressure_coupling != "ReferenceAndExcess")
            { OGS_FATAL("Unknown mechanical drainage pressure_coupling."); }
            if (type == "UndrainedTotalStress" && pressure_coupling != "ClipActual")
            { OGS_FATAL("Undrained C has no Biot pressure coupling."); }
            double const modulus = explicit_modulus.value_or(0.);
            bool const total_stress = type == "UndrainedTotalStress";
            if (total_stress && (explicit_modulus || fluid_modulus))
            {
                OGS_FATAL("UndrainedTotalStress uses calibrated total-stress stiffness, not Biot or fluid bulk storage.");
            }
            if ((explicit_modulus && fluid_modulus) ||
                (explicit_modulus && (!std::isfinite(modulus) || modulus <= 0.)) ||
                (fluid_modulus && (!std::isfinite(*fluid_modulus) || *fluid_modulus <= 0.)) ||
                (type == "UndrainedEffectiveStress" && !explicit_modulus && !fluid_modulus))
            {
                OGS_FATAL("Undrained mechanical drainage requires positive finite biot_modulus or compressible density EOS, never both. Drained material may omit storage.");
            }
            if (id < 0 || !media.contains(id) ||
                (type != "DrainedEffectiveStress" && type != "UndrainedEffectiveStress" && !total_stress))
            {
                OGS_FATAL("Mechanical drainage requires a known MaterialID and a supported effective-stress or UndrainedTotalStress mode.");
            }
            if (!settings.materials.emplace(id, MechanicalDrainageMaterial{
                total_stress ? MechanicalDrainageType::UndrainedTotalStress :
                (type == "DrainedEffectiveStress" ? MechanicalDrainageType::DrainedEffectiveStress : MechanicalDrainageType::UndrainedEffectiveStress),
                modulus, fluid_modulus.value_or(0.), pressure_coupling == "ReferenceAndExcess"}).second)
            {
                OGS_FATAL("Duplicate mechanical drainage MaterialID.");
            }
        }
        for (auto const& [id, medium] : media)
        {
            if (!settings.materials.contains(id))
            {
                OGS_FATAL("Mechanical drainage must explicitly cover every medium MaterialID.");
            }
            for (auto const property : {MaterialPropertyLib::PropertyType::biot_coefficient, MaterialPropertyLib::PropertyType::porosity,
                MaterialPropertyLib::PropertyType::reference_temperature, MaterialPropertyLib::PropertyType::permeability})
            {
                if (!dynamic_cast<MaterialPropertyLib::Constant const*>(&medium->property(property)))
                {
                    OGS_FATAL("Mechanical drainage requires Constant medium properties in this first scope.");
                }
            }
            double const alpha = medium->property(MaterialPropertyLib::PropertyType::biot_coefficient)
                .template value<double>(MaterialPropertyLib::EmptyVariableArray, ParameterLib::SpatialPosition{}, 0., 0.);
            if (!std::isfinite(alpha) || alpha <= 0 || alpha > 1)
            {
                OGS_FATAL("Mechanical drainage requires 0 < alpha <= 1.");
            }
            auto const& storage = settings.materials.at(id);
            if (storage.fluid_bulk_modulus != 0.)
            {
                double const porosity = medium->property(MaterialPropertyLib::PropertyType::porosity)
                    .template value<double>(MaterialPropertyLib::EmptyVariableArray, ParameterLib::SpatialPosition{}, 0., 0.);
                if (!std::isfinite(porosity) || porosity <= 0. || porosity >= 1. || alpha < porosity)
                {
                    OGS_FATAL("Derived mechanical storage requires 0 < porosity <= alpha <= 1 and porosity < 1.");
                }
                INFO("Mechanical drainage MaterialID {} derives linearised Biot storage from density-EOS reference bulk stiffness={} Pa, porosity={}, alpha={} and the drained elastic bulk modulus.", id, storage.fluid_bulk_modulus, porosity, alpha);
            }
            auto const& fluid = fluidPhase(*medium);
            auto const& solid = medium->phase(MaterialPropertyLib::PhaseName::Solid);
            if ((!dynamic_cast<MaterialPropertyLib::Constant const*>(&fluid.property(MaterialPropertyLib::PropertyType::density)) &&
                 !dynamic_cast<MaterialPropertyLib::Linear const*>(&fluid.property(MaterialPropertyLib::PropertyType::density))) ||
                !dynamic_cast<MaterialPropertyLib::Constant const*>(&fluid.property(MaterialPropertyLib::PropertyType::viscosity)) ||
                !dynamic_cast<MaterialPropertyLib::Constant const*>(&solid.property(MaterialPropertyLib::PropertyType::density)))
            {
                OGS_FATAL("Mechanical drainage uses Constant/Linear liquid density and Constant solid density/viscosity; local storage is explicit or linearised from the density EOS at reference pressure zero.");
            }
        }
        process_data.mechanical_drainage = std::move(settings);
    }


    if (process_data.safety && process_data.safety->usesActualPressure())
    {
        auto const& safety=*process_data.safety;
        if (!process_data.mechanical_drainage)
        {
            process_data.mechanical_drainage=MechanicalDrainageSettings{safety.physical_time,{}};
            for (auto const& [id,medium]:media)
            {
                double const kf=densityDerivedModulus(*medium).value_or(0.);
                process_data.mechanical_drainage->materials.emplace(id,MechanicalDrainageMaterial{
                    MechanicalDrainageType::DrainedEffectiveStress,0.,kf});
            }
        }
        auto& settings=*process_data.mechanical_drainage;
        settings.safety_reference_time=safety.physical_time;
        for (auto const& [id,medium]:media)
        {
            auto material=settings.materials.at(id);
            if (material.usesTotalStress())
            {
                if (safety.pressure_policy!="MaterialDrainage")
                { OGS_FATAL("Actual-pressure safety cannot convert TSA stress basis."); }
                settings.safety_materials.emplace(id,material);
                continue;
            }
            if (safety.pressure_policy!="MaterialDrainage")
            {
                material.type=safety.pressure_policy=="AdditionalUndrained" ?
                    MechanicalDrainageType::UndrainedEffectiveStress : MechanicalDrainageType::DrainedEffectiveStress;
            }
            if (material.type==MechanicalDrainageType::UndrainedEffectiveStress &&
                material.biot_modulus==0. && material.fluid_bulk_modulus==0.)
            {
                double const alpha=medium->property(MaterialPropertyLib::PropertyType::biot_coefficient)
                    .template value<double>(MaterialPropertyLib::EmptyVariableArray,ParameterLib::SpatialPosition{},0.,0.);
                double const n=medium->property(MaterialPropertyLib::PropertyType::porosity)
                    .template value<double>(MaterialPropertyLib::EmptyVariableArray,ParameterLib::SpatialPosition{},0.,0.);
                if (!(alpha>n && alpha<1.))
                { OGS_FATAL("AdditionalUndrained needs finite positive flow-compatible storage: supply pressure-dependent fluid density when grain compliance is zero."); }
            }
            settings.safety_materials.emplace(id,material);
        }
    }

    process_data.use_positive_pore_pressure_in_mechanics =
        config.getConfigParameter<bool>("use_positive_pore_pressure_in_mechanics", false);
    if (process_data.use_positive_pore_pressure_in_mechanics &&
        (DisplacementDim != 2 || mesh.isAxiallySymmetric() ||
         !process_data.isMonolithicSchemeUsed() || is_linear ||
         (process_data.safety && process_data.safety->pressure_policy != "FrozenAcceptedField" &&
             !process_data.safety->usesActualPressure())))
    {
        OGS_FATAL("Positive mechanical pore pressure currently requires nonlinear monolithic plane-strain HM; ordinary Safety/checkpoints use FrozenAcceptedField; actual-pressure Safety captures the accepted signed fluid state and requires effective-stress materials.");
    }

    if (process_data.use_positive_pore_pressure_in_mechanics && process_data.mechanical_drainage)
    {
        for (auto const& [id, material] : process_data.mechanical_drainage->materials)
        {
            if (material.usesTotalStress())
            { OGS_FATAL("Positive Biot pressure requires effective-stress materials; TSA has no Biot pressure."); }
        }
    }

    // These engineering profiles are qualified only from a virgin state.
    // A result VTU is not a complete constitutive/host/BC restart contract.
    bool requires_virgin_history = process_data.initialisation_gate.has_value() || process_data.native_checkpoint.has_value() || process_data.safety.has_value();
    if (process_data.mechanical_drainage)
    {
        for (auto const& [id, material] : process_data.mechanical_drainage->materials)
        {
            requires_virgin_history |= material.usesTotalStress();
        }
    }
    if (requires_virgin_history)
    {
        for (auto const& name : mesh.getProperties().getPropertyVectorNames())
        {
            if (name == "sigma_ip" || name == "epsilon_ip" ||
                name == "strain_rate_variable_ip" ||
                name == "mechanical_pressure_ip" ||
                name.starts_with("mechanical_reference_") ||
                name.starts_with("material_state_variable_"))
            {
                OGS_FATAL("Virgin-gravity TSA/equilibrium_check cannot import mechanical history '{}' from the bulk mesh. Complete accepted-state restart is not implemented; use an unseeded mesh.", name);
            }
        }
    }

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<HydroMechanicsProcess<DisplacementDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        std::holds_alternative<Monolithic>(coupling_scheme), is_linear);
}

template std::unique_ptr<Process> createHydroMechanicsProcess<2>(
    std::string const& name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

template std::unique_ptr<Process> createHydroMechanicsProcess<3>(
    std::string const& name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

}  // namespace HydroMechanics
}  // namespace ProcessLib
