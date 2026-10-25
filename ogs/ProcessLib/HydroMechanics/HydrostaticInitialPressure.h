// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <boost/property_tree/xml_parser.hpp>
#include <algorithm>
#include <initializer_list>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <set>
#include <regex>
#include <vector>

#include "MeshLib/Elements/Utils.h"
#include "MeshLib/Mesh.h"

namespace ProcessLib::HydroMechanics
{
namespace InitialGroundwaterDetail
{
using boost::property_tree::ptree;

[[noreturn]] inline void fail(std::string const& message)
{
    throw std::runtime_error("HM_INITIAL_GROUNDWATER: " + message);
}

inline void keys(ptree const& tree,
                 std::initializer_list<std::string> allowed,
                 std::initializer_list<std::string> repeated = {})
{
    for (auto const& [key, unused] : tree)
    {
        if (std::find(allowed.begin(), allowed.end(), key) == allowed.end())
        {
            fail("unsupported entry " + key);
        }
        if (tree.count(key) != 1 &&
            std::find(repeated.begin(), repeated.end(), key) == repeated.end())
        {
            fail("duplicate entry " + key);
        }
    }
}

inline double finiteNumber(ptree const& tree, std::string const& key)
{
    auto const text = tree.get<std::string>(key);
    std::size_t end = 0;
    double value;
    try
    {
        value = std::stod(text, &end);
    }
    catch (std::exception const&)
    {
        fail("invalid number for " + key);
    }
    if (end != text.size() || !std::isfinite(value))
    {
        fail("finite number required for " + key);
    }
    return value;
}

// Boundary water levels use the EOS reference density (p_ref=0 gauge).
// This is a boundary-pressure convention, not a constant-density flow model.
inline double referenceLiquidDensity(ptree const& property)
{
    auto const type=property.get<std::string>("type");
    if (type=="Constant") { return finiteNumber(property,"value"); }
    if (type!="Linear" || property.count("independent_variable")!=1)
    { fail("native water levels require Constant or pressure-Linear liquid density"); }
    keys(property,{"name","type","reference_value","independent_variable"});
    auto const& variable=property.get_child("independent_variable");
    keys(variable,{"variable_name","reference_condition","slope"});
    if (variable.get<std::string>("variable_name")!="liquid_phase_pressure" ||
        finiteNumber(variable,"reference_condition")!=0. || finiteNumber(variable,"slope")<0.)
    { fail("native water levels require Linear liquid_phase_pressure density referenced at zero gauge pressure with nonnegative slope"); }
    return finiteNumber(property,"reference_value");
}

inline std::string number(double const value)
{
    std::ostringstream text;
    text << std::setprecision(std::numeric_limits<double>::max_digits10)
         << value;
    return text.str();
}
}  // namespace InitialGroundwaterDetail

}  // namespace ProcessLib::HydroMechanics
#include "EngineeringPhases.h"
#include "EngineeringHydraulics.h"
#include "AutomaticReservoirSurface.h"
#include "ReservoirCatalogue.h"
#include "EngineeringSeepage.h"
#include "EngineeringSafety.h"
#include "EngineeringCheckpoints.h"
#include "EngineeringExcavation.h"
#include "EngineeringConstruction.h"
namespace ProcessLib::HydroMechanics
{

struct GravityInitialisation
{
    std::string bulk_mesh;
    std::string pressure_mesh;
};

/// Check actual loaded mesh coverage; a boundary subset cannot implement a
/// fixed initial field throughout the domain.
inline void validateGravityPressureMesh(GravityInitialisation const& gravity,
                                       MeshLib::Mesh const& bulk,
                                       MeshLib::Mesh const& pressure)
{
    using namespace InitialGroundwaterDetail;
    auto const* ids = pressure.getProperties().getPropertyVector<std::size_t>(
        "bulk_node_ids", MeshLib::MeshItemType::Node, 1);
    if (!ids || ids->size() != pressure.getNumberOfNodes())
    {
        fail("gravity pressure mesh requires UInt64 bulk_node_ids for every node");
    }
    std::set<std::size_t> covered;
    for (std::size_t i = 0; i < ids->size(); ++i)
    {
        auto const id = (*ids)[i];
        if (id >= bulk.getNumberOfNodes() || !covered.insert(id).second)
        {
            fail("gravity pressure mesh has duplicate or invalid bulk_node_ids");
        }
        if ((pressure.getNode(i)->asEigenVector3d() -
             bulk.getNode(id)->asEigenVector3d()).norm() > 1e-12)
        {
            fail("gravity pressure mesh node coordinates disagree with bulk map");
        }
    }
    for (auto const* node : MeshLib::getBaseNodes(bulk.getElements()))
    {
        if (!covered.contains(node->getID()))
        {
            fail("FixedInitialField pressure mesh must cover every primary pressure node");
        }
    }
    if (bulk.getProperties().hasPropertyVector("sigma_ip") ||
        bulk.getProperties().hasPropertyVector("epsilon_ip"))
    {
        fail("Fixed-field initialisation initialisation starts from virgin material, not an IP restart");
    }
}

/// Lower named hydrostatic initial and boundary fields to stock parameters.
/// Initial and optional boundary bindings are explicit; constitutive state is unchanged.
inline bool prepareHydrostaticInitialPressure(std::stringstream& input,
                                             std::string const& output_directory,
                                             std::optional<GravityInitialisation>& gravity,
                                             std::filesystem::path const& project_directory)
{
    using namespace InitialGroundwaterDetail;
    if (input.str().find("<native_checkpoint") == std::string::npos &&
        input.str().find("<native_safety") == std::string::npos &&
        input.str().find("<checkpoints") == std::string::npos &&
        input.str().find("<restart") == std::string::npos &&
        input.str().find("<groundwater_states") == std::string::npos &&
        input.str().find("<initial_pressure_guess") == std::string::npos &&
        input.str().find("<hydraulic_boundaries") == std::string::npos &&
        input.str().find("<surface_water_") == std::string::npos &&
        input.str().find("UndrainedTotalStress") == std::string::npos &&
        input.str().find("<mechanical_drainage") == std::string::npos &&
        input.str().find("<displacement_reference_times") == std::string::npos &&
        input.str().find("<adaptive_rollback") == std::string::npos &&
        input.str().find("<steady_state_seepage") == std::string::npos &&
        input.str().find("<transient_seepage") == std::string::npos &&
        input.str().find("<initialisation_gate") == std::string::npos &&
        input.str().find("<initialisation_load_ramp_end") == std::string::npos)
    {
        return false;
    }
    std::stringstream parsed(input.str());
    ptree document;
    boost::property_tree::read_xml(
        parsed, document,
        boost::property_tree::xml_parser::trim_whitespace |
            boost::property_tree::xml_parser::no_comments);
    auto& root = document.get_child("OpenGeoSysProject");
    auto& processes = root.get_child("processes");
    bool has_catalogue = false;
    for (auto const& [tag, candidate] : processes)
    {
        if (tag == "process" && candidate.get<std::string>("type", "") == "HYDRO_MECHANICS" &&
            (candidate.count("steady_state_seepage") || candidate.count("transient_seepage") || candidate.count("adaptive_rollback") || candidate.count("construction_births") || candidate.count("native_safety") || candidate.count("native_checkpoint") || candidate.count("displacement_reference_times") || candidate.count("initialisation_gate") || candidate.count("initialisation_load_ramp_end")))
        {
            fail("Reserved HM lowering keys must not be supplied directly; use initialisation/equilibrium_check and load_ramp.");
        }
        if (tag == "process" && candidate.get_child_optional("initialisation.initial_pressure_guess"))
        { has_catalogue = true; }
        if (tag == "process")
        {
            if (auto phases = candidate.get_child_optional("phases"))
            {
                for (auto const& [phase_tag, phase] : *phases)
                { has_catalogue |= phase.count("initial_pressure_guess") != 0; }
            }
        }
        if (tag == "process" && (candidate.count("groundwater_states") ||
                                candidate.count("hydraulic_boundaries") ||
                                candidate.count("surface_water_states") ||
                                candidate.count("surface_water_loads") ||
                                candidate.count("mechanical_drainage")))
        {
            has_catalogue = true;
        }
    }
    if (!has_catalogue)
    {
        return false;
    }
    if (processes.count("process") != 1)
    {
        fail("one HM process required for this feature");
    }
    auto& process = processes.get_child("process");
    // Phase-owned initialisation is canonicalised before checkpoint identity and
    // all existing physical lowering. The legacy process-level form remains valid.
    if (auto phases = process.get_child_optional("phases"))
    {
        bool first = true;
        for (auto& [tag, phase] : *phases)
        {
            bool owns_initialisation = false;
            for (auto const key : {"method", "groundwater", "initial_pressure_guess", "load_ramp",
                                   "pressure_mesh", "equilibrium_check",
                                   "stress_basis", "k0",
                                   "vertical_effective_stress_parameter",
                                   "stress_parameter"})
            {
                owns_initialisation |= phase.count(key) != 0;
            }
            if (owns_initialisation)
            {
                if (tag != "phase" || !first ||
                    phase.get<std::string>("type", "") != "Initialisation")
                {
                    fail("initialisation settings belong only to the first Initialisation phase");
                }
                if (process.count("initialisation"))
                {
                    fail("choose phase-owned or process-level initialisation, not both");
                }
                if (phase.count("method") != 1)
                {
                    fail("phase-owned initialisation requires exactly one method");
                }
                ptree initial;
                for (auto const key : {"method", "groundwater", "hydraulic_policy",
                                       "initial_pressure_guess", "load_ramp", "pressure_mesh",
                                       "equilibrium_check", "stress_basis", "k0",
                                       "vertical_effective_stress_parameter",
                                       "stress_parameter"})
                {
                    if (phase.count(key) > 1)
                    {
                        fail("duplicate phase-owned initialisation setting " + std::string(key));
                    }
                    if (auto setting = phase.get_child_optional(key))
                    {
                        initial.put_child(key, *setting);
                        phase.erase(key);
                    }
                }
                process.put_child("initialisation", initial);
            }
            first = false;
        }
    }
    bool const zero_initial_guess = static_cast<bool>(
        process.get_child_optional("initialisation.initial_pressure_guess"));
    if (zero_initial_guess)
    {
        auto& initial = process.get_child("initialisation");
        keys(initial.get_child("initial_pressure_guess"), {});
        if (initial.count("initial_pressure_guess") != 1 ||
            initial.get<std::string>("initial_pressure_guess") != "Zero" ||
            initial.get<std::string>("method", "") != "Gravity" ||
            (initial.get<std::string>("hydraulic_policy", "") != "SteadyStateSeepage" &&
             initial.get<std::string>("hydraulic_policy", "") != "FixedInitialField"))
        { fail("initial_pressure_guess Zero requires Gravity with SteadyStateSeepage or FixedInitialField"); }
        if (initial.count("groundwater") || process.count("groundwater_states"))
        { fail("choose Zero pressure guess or a groundwater-level guess, not both"); }
        // Reuse existing catalogue lowering privately; the generated initial
        // parameter below is literal zero, never rho*g*(0-y).
        ptree state, states;
        state.put("<xmlattr>.id", "InternalZeroPressureGuess");
        state.put("type", "Hydrostatic");
        state.put("water_level", "0");
        state.put("pressure_policy", "Signed");
        states.push_back({"state", state});
        process.put_child("groundwater_states", states);
        initial.erase("initial_pressure_guess");
        initial.put("groundwater.<xmlattr>.ref", "InternalZeroPressureGuess");
    }
    // User-facing drainage names are lowered before any stress-basis routing.
    // A/B share the local water-content law, but B has calibrated su strength.
    if (auto materials = process.get_child_optional("mechanical_drainage"))
    {
        for (auto& [tag, material] : *materials)
        {
            if (tag != "material") { continue; }
            auto const type = material.get<std::string>("type");
            if (type == "UndrainedB")
            {
                int const id = material.get<int>("<xmlattr>.id");
                bool found_relation = false;
                for (auto const& [relation_tag, relation] : process)
                {
                    if (relation_tag != "constitutive_relation") { continue; }
                    auto const relation_id = relation.get_optional<int>("<xmlattr>.id");
                    if (relation_id && *relation_id != id) { continue; }
                    found_relation = true;
                    for (auto const& [property_tag, property] : relation.get_child("material_properties"))
                    {
                        auto const name = property.get<std::string>("<xmlattr>.name");
                        if (name != "FrictionAngle" && name != "DilatancyAngle") { continue; }
                        auto const parameter_name = property.get<std::string>("<xmlattr>.parameter");
                        bool zero = false;
                        for (auto const& [parameter_tag, parameter] : root.get_child("parameters"))
                        {
                            if (parameter.get<std::string>("name") != parameter_name) { continue; }
                            if (parameter.get<std::string>("type") == "Constant")
                            {
                                zero = finiteNumber(parameter, parameter.count("values") ? "values" : "value") == 0.;
                            }
                        }
                        if (!zero)
                        { fail("UndrainedB requires Constant zero FrictionAngle and DilatancyAngle; Cohesion is calibrated su and stiffness is effective E/nu."); }
                    }
                }
                if (!found_relation) { fail("UndrainedB has no matching constitutive relation"); }
            }
            if (type == "UndrainedA" || type == "UndrainedB")
            {
                if (material.get<std::string>("pressure_coupling", "ReferenceAndExcess") != "ReferenceAndExcess")
                { fail("Undrained A/B require ReferenceAndExcess coupling to retain water stiffness above the phreatic surface"); }
                material.put("pressure_coupling", "ReferenceAndExcess");
                material.put("type", "UndrainedEffectiveStress");
            }
            else if (type == "UndrainedC")
            { material.put("type", "UndrainedTotalStress"); }
            else if (type == "Drained")
            { material.put("type", "DrainedEffectiveStress"); }
        }
    }
    bool has_total_stress = false;
    bool has_other_basis = false;
    if (auto const materials = process.get_child_optional("mechanical_drainage"))
    {
        for (auto const& [tag, material] : *materials)
        {
            if (tag != "material") { continue; }
            bool const tsa = material.get<std::string>("type", "") == "UndrainedTotalStress";
            has_total_stress |= tsa;
            has_other_basis |= !tsa;
        }
    }
    bool const homogeneous_tsa = has_total_stress && !has_other_basis;
    bool const tsa_placeholder = homogeneous_tsa && !process.count("groundwater_states");
    if (homogeneous_tsa && process.count("initialisation") &&
        !process.get_child_optional("initialisation.hydraulic_policy"))
    {
        process.put("initialisation.hydraulic_policy", "FixedInitialField");
    }
    if (tsa_placeholder)
    {
        if (process.get_child_optional("initialisation.groundwater"))
        {
            fail("TSA groundwater reference has no matching catalogue; omit both for the automatic inactive zero-pressure field.");
        }
        ptree state;
        state.put("<xmlattr>.id", "InternalTsaPlaceholder");
        state.put("type", "Hydrostatic");
        state.put("water_level", "0");
        state.put("pressure_policy", "Signed");
        ptree states;
        states.push_back({"state", state});
        process.put_child("groundwater_states", states);
        process.put("initialisation.groundwater.<xmlattr>.ref", "InternalTsaPlaceholder");
    }
    if (process.get<std::string>("type") != "HYDRO_MECHANICS" ||
        process.count("groundwater_states") != 1 ||
        process.count("initialisation") != 1)
    {
        fail("one HM groundwater catalogue and explicit initialisation required");
    }
    prepareReservoirCatalogue(process);
    prepareEngineeringConstruction(root,process);
    auto const checkpoint_schedule = prepareEngineeringCheckpoints(root,process,
        project_directory,output_directory.empty()?std::filesystem::path("."):std::filesystem::path(output_directory));
    auto const excavation_water_tractions=prepareEngineeringExcavation(root,process);
    prepareAutomaticReservoirSurfaces(root,process,project_directory,
        output_directory.empty()?std::filesystem::path("."):std::filesystem::path(output_directory));
    std::vector<std::string> drainage_phase_ids;
    if (auto phases = process.get_child_optional("phases"))
    {
        for (auto const& [tag, phase] : *phases)
        { if (tag == "phase") { drainage_phase_ids.push_back(phase.get<std::string>("<xmlattr>.id")); } }
    }
    prepareEngineeringSafety(root,process);
    auto const phase_schedule = prepareEngineeringPhases(root, process, true);
    prepareEngineeringSeepage(root, process, phase_schedule);
    if (auto drainage = process.get_child_optional("mechanical_drainage"))
    {
        for(auto const& [vtag,variable]:root.get_child("process_variables"))
        {
            if(auto sources=variable.get_child_optional("source_terms"))
            {
                for(auto const& [stag,source]:*sources)
                {
                    if(source.get<std::string>("type","")=="Interface2D" &&
                       source.get<bool>("use_hydraulic_pressure",false))
                    { fail("Wet Interface2D currently traces ordinary HM primary pressure; strict material actual-IP pressure is unsupported."); }
                }
            }
        }
        for (auto const& [tag, unused] : *drainage)
        {
            if (tag != "reference" && tag != "material")
            {
                fail("mechanical_drainage requires reference and material entries");
            }
        }
        if (drainage->count("reference") != 1 || phase_schedule.times.size() < 2 ||
            !process.get_child_optional("initialisation.method"))
        { fail("mechanical_drainage requires an accepted phase reference and Gravity initialisation"); }
        auto const reference = drainage->get<std::string>("reference");
        std::size_t reference_index = 1;
        if (reference != "InitialisationEnd")
        {
            auto const found = std::find(drainage_phase_ids.begin(), drainage_phase_ids.end(), reference);
            if (found == drainage_phase_ids.end()) { fail("Unknown mechanical_drainage reference phase " + reference); }
            reference_index = std::distance(drainage_phase_ids.begin(), found) + 1;
            if (reference_index >= phase_schedule.times.size()) { fail("Drainage reference must precede Safety continuation"); }
        }
        if (phase_schedule.hydraulic_phases)
        {
            for (std::size_t i = reference_index; i < phase_schedule.hydraulic_policies.size(); ++i)
            {
                if (phase_schedule.hydraulic_policies[i] != "FrozenAcceptedField")
                { fail("After a material drainage reference, hydraulic phases must use FrozenAcceptedField; Coupled/Steady flow is a different formulation"); }
            }
        }
        if (auto safety = process.get_optional<double>("native_safety.physical_time"))
        {
            if (phase_schedule.times.at(reference_index) > *safety)
            { fail("Mechanical drainage reference must be an accepted working phase before Safety"); }
        }
        std::set<std::string> static_parameters;
        bool total_stress = false;
        bool effective_stress = false;
        for (auto const& [tag, material] : *drainage)
        {
            if (tag != "material") { continue; }
            bool const tsa = material.get<std::string>("type") == "UndrainedTotalStress";
            total_stress |= tsa;
            effective_stress |= !tsa;
        }
        if (total_stress && (reference != "InitialisationEnd" || phase_schedule.hydraulic_phases))
        { fail("UndrainedC/TSA requires virgin InitialisationEnd reference without preceding seepage; accepted stress-basis conversion is unsupported"); }
        if (total_stress && (process.get<std::string>("initialisation.method") != "Gravity" ||
            process.count("initial_stress")))
        {
            fail("TSA requires virgin Gravity initialisation and no assigned initial stress; accepted stress-basis conversion is not supported.");
        }
        bool const mixed_basis = total_stress && effective_stress;
        std::map<int, bool> basis_by_material;
        for (auto const& [tag, material] : *drainage)
        {
            if (tag != "material") { continue; }
            int const id = material.get<int>("<xmlattr>.id");
            if (id < 0 || !basis_by_material.emplace(id,
                material.get<std::string>("type") == "UndrainedTotalStress").second)
            {
                fail("mechanical_drainage requires distinct nonnegative material IDs");
            }
        }
        if (mixed_basis || (total_stress && process.count("constitutive_relation") > 1))
        {
            std::set<int> constitutive_ids;
            for (auto const& [tag, relation] : process)
            {
                if (tag != "constitutive_relation") { continue; }
                auto const id = relation.get_optional<int>("<xmlattr>.id");
                if (!id || !basis_by_material.contains(*id) || !constitutive_ids.insert(*id).second)
                {
                    fail("Mixed stress bases require one explicit, distinct constitutive_relation id matching every mechanical material ID; no shared-law fallback.");
                }
            }
            if (constitutive_ids.size() != basis_by_material.size())
            {
                fail("Mixed stress bases require exact constitutive/material ID coverage.");
            }
        }
        keys(drainage->get_child("reference"), {});
        for (auto const& [tag, variable] : root.get_child("process_variables"))
        {
            if (variable.count("deactivated_subdomains"))
            {
                fail("mechanical_drainage does not yet support element activation/deactivation");
            }
        }
        for (auto const& [tag, relation] : process)
        {
            if (tag == "constitutive_relation" &&
                (relation.get<std::string>("type") != "MFront" ||
                 relation.get<std::string>("behaviour") != "MohrCoulombAbboSloan"))
            {
                fail("mechanical_drainage first scope supports MFront MohrCoulombAbboSloan only");
            }
        }
        for (auto const& [tag, relation] : process)
        {
            if (tag != "constitutive_relation") { continue; }
            bool const relation_is_tsa = total_stress &&
                (!mixed_basis || basis_by_material.at(relation.get<int>("<xmlattr>.id")));
            for (auto const& [property_tag, property] : relation.get_child("material_properties"))
            {
                auto const parameter_name = property.get<std::string>("<xmlattr>.parameter");
                bool constant = false;
                for (auto const& [parameter_tag, parameter] : root.get_child("parameters"))
                {
                    if (parameter.get<std::string>("name") == parameter_name)
                    {
                        constant = parameter.get<std::string>("type") == "Constant";
                        auto const property_name = property.get<std::string>("<xmlattr>.name");
                        if (!relation_is_tsa && parameter.get<std::string>("type") == "Function" &&
                            (property_name == "YoungModulus" || property_name == "Cohesion") &&
                            parameter.count("expression") == 1 &&
                            !std::regex_search(parameter.get<std::string>("expression"), std::regex("(^|[^a-zA-Z0-9_])t([^a-zA-Z0-9_]|$)")))
                        {
                            constant = true; // A runtime time-symbol check is also mandatory.
                            static_parameters.insert(parameter_name);
                        }
                        if (relation_is_tsa && constant)
                        {
                            std::string const key = parameter.count("values") ? "values" : "value";
                            double const value = finiteNumber(parameter, key);
                            if ((property_name == "YoungModulus" && value <= 0.) ||
                                (property_name == "PoissonRatio" && (value <= -1. || value >= .5)) ||
                                (property_name == "Cohesion" && value <= 0.))
                            {
                                fail("TSA requires finite Eu > 0, -1 < nuu < 0.5 and calibrated su > 0.");
                            }
                        }
                        if (relation_is_tsa && constant && (property_name == "FrictionAngle" || property_name == "DilatancyAngle"))
                        {
                            std::string const key = parameter.count("values") ? "values" : "value";
                            if (finiteNumber(parameter, key) != 0.)
                            {
                                fail("UndrainedTotalStress requires zero FrictionAngle and DilatancyAngle; Cohesion represents calibrated su and elastic parameters are total-stress Eu/nuu.");
                            }
                        }
                    }
                }
                if (!constant)
                {
                    fail("mechanical_drainage first scope requires Constant MC material parameters");
                }
            }
        }
        drainage->erase("reference");
        drainage->put("reference_time", number(phase_schedule.times.at(reference_index)));
        for (auto const& name : static_parameters)
        { ptree parameter; parameter.put("parameter", name); drainage->push_back({"static_parameter", parameter}); }
    }

    if (auto const check = process.get_child_optional("initialisation.equilibrium_check"))
    {
        keys(*check, {"force_tolerance", "yield_tolerance"});
        if (process.get<std::string>("initialisation.method", "") != "Gravity" ||
            phase_schedule.times.size() < 2 || process.count("initialisation_gate"))
        {
            fail("equilibrium_check requires a native virgin Gravity phase; initialisation_gate is reserved");
        }
        double const force = finiteNumber(*check, "force_tolerance");
        double const yield = finiteNumber(*check, "yield_tolerance");
        if (force <= 0. || yield <= 0.) { fail("equilibrium_check tolerances must be positive and finite"); }
        ptree gate;
        gate.put("reference_time", number(phase_schedule.times.at(1)));
        gate.put("force_tolerance", number(force));
        gate.put("yield_tolerance", number(yield));
        std::set<int> ids;
        for (auto const& [tag, relation] : process)
        {
            if (tag != "constitutive_relation") { continue; }
            if (relation.get<std::string>("type") != "MFront" ||
                relation.get<std::string>("behaviour") != "MohrCoulombAbboSloan")
            { fail("equilibrium_check requires Constant-property MohrCoulombAbboSloan"); }
            int const id = relation.get<int>("<xmlattr>.id", 0);
            if (id < 0 || !ids.insert(id).second) { fail("equilibrium_check requires distinct nonnegative constitutive IDs"); }
            std::map<std::string, double> values;
            std::optional<std::string> spatial_cohesion;
            for (auto const& [unused, property] : relation.get_child("material_properties"))
            {
                auto const name = property.get<std::string>("<xmlattr>.parameter");
                bool found = false;
                for (auto const& [unused_parameter, parameter] : root.get_child("parameters"))
                {
                    if (parameter.get<std::string>("name") != name) { continue; }
                    auto const property_name=property.get<std::string>("<xmlattr>.name");
                    bool const spatial=parameter.get<std::string>("type")=="Function" &&
                        (property_name=="Cohesion" || property_name=="YoungModulus") &&
                        parameter.count("expression")==1 &&
                        !std::regex_search(parameter.get<std::string>("expression"),std::regex("(^|[^a-zA-Z0-9_])t([^a-zA-Z0-9_]|$)"));
                    if (found || (parameter.get<std::string>("type") != "Constant" && !spatial))
                    { fail("equilibrium_check requires unique Constant MC parameters or static scalar Function E/cohesion"); }
                    found = true;
                    if (spatial)
                    {
                        if (property_name=="Cohesion") { spatial_cohesion=name; }
                        values.emplace(property_name,1.); // Cohesion is validated at each IP in the runtime gate.
                    }
                    else { values.emplace(property_name,finiteNumber(parameter, parameter.count("values") ? "values" : "value")); }
                }
                if (!found) { fail("equilibrium_check cannot resolve MC parameter " + name); }
            }
            ptree material;
            material.put("<xmlattr>.id", id);
            if (spatial_cohesion) { material.put("cohesion_parameter",*spatial_cohesion); }
            for (auto const& [property, key] : std::map<std::string, std::string>{
                {"Cohesion", "cohesion"}, {"FrictionAngle", "friction_angle"},
                {"TransitionAngle", "transition_angle"}, {"TensionCutOffParameter", "cutoff"}})
            {
                if (!values.contains(property)) { fail("equilibrium_check missing MC property " + property); }
                material.put(key, number(values.at(property)));
            }
            double const c = values.at("Cohesion");
            double const phi = values.at("FrictionAngle");
            double const trans = values.at("TransitionAngle");
            double const cutoff = values.at("TensionCutOffParameter");
            double const cap = phi > 0. ? std::min(100., .1*c/std::tan(phi*std::numbers::pi/180.)) : 100.;
            if (c <= 0. || phi < 0. || phi >= 90. || trans <= 0. || trans >= 30. || cutoff < 0. || (!spatial_cohesion && cutoff > cap))
            { fail("equilibrium_check requires admissible MC parameters and the declared cutoff cap"); }
            gate.push_back({"material", material});
        }
        for (auto const& [tag, variable] : root.get_child("process_variables"))
        {
            if (variable.count("deactivated_subdomains")) { fail("equilibrium_check does not yet support activation/deactivation"); }
        }
        process.put_child("initialisation_gate", gate);
    }
    double ramp_end = 0.;
    if (auto const ramp = process.get_child_optional("initialisation.load_ramp"))
    {
        keys(*ramp, {});
        if (ramp->data() != "Synchronous" ||
            process.get<std::string>("initialisation.method", "") != "Gravity" ||
            phase_schedule.times.size() < 2)
        {
            fail("load_ramp Synchronous requires native virgin Gravity initialisation");
        }
        if (process.count("initialisation_load_ramp_end"))
        {
            fail("reserved initialisation_load_ramp_end collision");
        }
        ramp_end = phase_schedule.times.at(1);
        process.put("initialisation_load_ramp_end", number(ramp_end));
    }
    bool const steady_initial = process.get<std::string>("initialisation.hydraulic_policy", "") == "SteadyStateSeepage";
    if (steady_initial && process.get_child("initialisation").count("pressure_mesh"))
    { fail("SteadyStateSeepage uses boundary heads, not an all-domain pressure_mesh."); }
    if (steady_initial && (ramp_end <= 0. || !phase_schedule.hydraulic_phases ||
        phase_schedule.hydraulic_policies.front() != "SteadyStateSeepage"))
    { fail("SteadyStateSeepage initialisation requires load_ramp Synchronous and an initial SteadyStateSeepage hydraulic phase."); }
    if (phase_schedule.hydraulic_phases &&
        phase_schedule.hydraulic_policies.front() != process.get<std::string>("initialisation.hydraulic_policy"))
    { fail("Initialisation and first phase hydraulic_policy must match."); }
    std::string const ramp_factor = ramp_end > 0.
        ? "min(max(t/" + number(ramp_end) + ",0),1)*(" : "";
    std::string const ramp_close = ramp_end > 0. ? ")" : "";
    auto const& initial = process.get_child("initialisation");
    bool const fixed_field_mode = initial.count("method") != 0;
    std::string method;
    if (fixed_field_mode)
    {
        method = initial.get<std::string>("method");
        if (method == "Gravity")
        {
            keys(initial, {"method", "groundwater", "hydraulic_policy", "pressure_mesh", "load_ramp", "equilibrium_check"});
        }
        else if (method == "K0")
        {
            keys(initial, {"method", "groundwater", "hydraulic_policy", "pressure_mesh",
                           "stress_basis", "k0", "vertical_effective_stress_parameter"});
            keys(initial.get_child("k0"), {});
            keys(initial.get_child("vertical_effective_stress_parameter"), {});
            if (initial.get<std::string>("stress_basis") != "Effective" ||
                finiteNumber(initial, "k0") < 0)
            {
                fail("K0 requires explicit Effective stress basis and finite nonnegative k0");
            }
        }
        else if (method == "AssignedStress")
        {
            keys(initial, {"method", "groundwater", "hydraulic_policy", "pressure_mesh",
                           "stress_basis", "stress_parameter"});
            keys(initial.get_child("stress_parameter"), {});
            auto const basis = initial.get<std::string>("stress_basis");
            if (basis != "Effective" && basis != "Total")
            {
                fail("AssignedStress requires explicit Effective or Total stress_basis");
            }
        }
        else
        {
            fail("unsupported initialisation method " + method);
        }
        if (method != "Gravity")
        {
            keys(initial.get_child("stress_basis"), {});
            // Every region must support a transactional, admissible stress seed.
            // Coverage and property validation remain owned by process creation.
            for (auto const& [tag, law] : process)
            {
                if (tag != "constitutive_relation") { continue; }
                auto const behaviour = law.get<std::string>("behaviour", "");
                bool const elastic = behaviour == "StandardElasticityBrick";
                bool const mc = (behaviour == "MohrCoulombAbboSloan" ||
                                 behaviour == "MohrCoulombAbboSloanTension") &&
                                law.get_child_optional("substepping");
                bool const hs = behaviour == "PublishedHS" &&
                                law.get_child_optional("substepping");
                if (law.get<std::string>("type", "") != "MFront" ||
                    (!elastic && !mc && !hs))
                {
                    fail("stress seed requires MFront elasticity, transactional MC or PublishedHS in every material region");
                }
            }
        }
        keys(initial.get_child("method"), {});
        keys(initial.get_child("hydraulic_policy"), {});
        if (initial.count("pressure_mesh"))
        {
            keys(initial.get_child("pressure_mesh"), {});
        }
        if (initial.get<std::string>("hydraulic_policy") != "FixedInitialField" && !steady_initial)
        {
            fail("initialisation requires FixedInitialField or SteadyStateSeepage");
        }
        if (process.count("initial_stress") || process.count("hydraulic_boundaries"))
        {
            fail("initialisation owns pressure and stress seeds; remove competing initial_stress or hydraulic_boundaries");
        }
        auto const pressure_mesh = initial.get<std::string>("pressure_mesh", "");
        bool exists = false;
        std::string bulk_mesh;
        for (auto const& [tag, entry] : root.get_child("meshes"))
        {
            if (tag != "mesh")
            {
                continue;
            }
            auto const name = std::filesystem::path(entry.data()).stem().string();
            if (bulk_mesh.empty())
            {
                bulk_mesh = name;
            }
            exists = exists || name == pressure_mesh;
        }
        if (!pressure_mesh.empty() && (!exists || pressure_mesh == bulk_mesh))
        {
            fail("Fixed-field initialisation pressure_mesh must name a separate declared mapped mesh");
        }
        gravity = GravityInitialisation{bulk_mesh, pressure_mesh};
        if (root.get<std::string>("time_loop.processes.process.time_stepping.type")=="PhaseAdaptive")
        {
            if (phase_schedule.times.size()<3 && !steady_initial) { fail("Adaptive initialisation needs a following working/hold phase unless it is standalone steady seepage."); }
        }
        else
        {
        auto const& time_process = root.get_child("time_loop.processes.process");
        if (root.get_child("time_loop.processes").count("process") != 1 ||
            time_process.get<std::string>("time_stepping.type") != "FixedTimeStepping")
        {
            fail("Fixed-field initialisation verification requires one FixedTimeStepping process");
        }
        auto const start = finiteNumber(time_process, "time_stepping.t_initial");
        auto const finish = finiteNumber(time_process, "time_stepping.t_end");
        double elapsed = 0;
        std::size_t steps = 0;
        for (auto const& [tag, pair] : time_process.get_child("time_stepping.timesteps"))
        {
            if (tag != "pair")
            {
                fail("Fixed-field initialisation requires explicit timestep pairs");
            }
            auto const repeat = pair.get<int>("repeat");
            double const dt = finiteNumber(pair, "delta_t");
            if (repeat <= 0 || dt <= 0 || start + elapsed + dt >= finish)
            {
                if (steps == 0)
                {
                    fail("Fixed-field initialisation needs an equilibrium step followed by an unchanged hold");
                }
            }
            if (repeat <= 0 || dt <= 0)
            {
                fail("positive gravity timestep count and size required");
            }
            steps += repeat;
            elapsed += repeat * dt;
        }
        if (steps < 2 || finish <= start || elapsed < finish - start)
        {
            fail("Fixed-field initialisation requires a complete schedule with equilibrium and hold steps");
        }
        }
    }
    else
    {
        keys(initial, {"groundwater"});
    }
    auto const& binding = initial.get_child("groundwater");
    keys(binding, {"<xmlattr>"});
    keys(binding.get_child("<xmlattr>"), {"ref"});
    auto const reference = binding.get<std::string>("<xmlattr>.ref");
    std::map<std::string, double> levels;
    std::map<std::string, std::string> elevations;
    std::map<std::string, bool> clipped;
    for (auto const& [tag, state] : process.get_child("groundwater_states"))
    {
        if (tag != "state") { fail("catalogue requires state entries"); }
        keys(state, {"<xmlattr>", "type", "water_level", "pressure_policy", "points", "outside_extent"});
        keys(state.get_child("<xmlattr>"), {"id"});
        keys(state.get_child("type"), {});
        keys(state.get_child("pressure_policy"), {});
        auto const id = state.get<std::string>("<xmlattr>.id");
        if (id.empty() || id.find_first_not_of(
            "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_") != std::string::npos ||
            elevations.contains(id))
        { fail("groundwater state requires a unique nonempty identifier"); }
        auto const type = state.get<std::string>("type");
        auto const policy = state.get<std::string>("pressure_policy");
        if (policy != "Signed" && policy != "CompressionOnly")
        { fail("groundwater pressure_policy must be Signed or CompressionOnly"); }
        clipped[id] = policy == "CompressionOnly";
        if (type == "Hydrostatic")
        {
            if (state.count("points") || state.count("outside_extent"))
            { fail("Hydrostatic uses water_level, not phreatic points"); }
            keys(state.get_child("water_level"), {});
            double const level = finiteNumber(state, "water_level");
            levels.emplace(id, level);
            elevations.emplace(id, number(level));
        }
        else if (type == "PhreaticLine")
        {
            if (state.count("water_level") || state.get<std::string>("outside_extent", "") != "Constant")
            { fail("PhreaticLine requires points and explicit outside_extent Constant, not water_level"); }
            keys(state.get_child("outside_extent"), {});
            std::vector<std::pair<double, double>> points;
            for (auto const& [point_tag, point] : state.get_child("points"))
            {
                if (point_tag != "point") { fail("PhreaticLine points requires point entries"); }
                keys(point, {"x", "y"});
                keys(point.get_child("x"), {});keys(point.get_child("y"), {});
                double const x = finiteNumber(point, "x"), y = finiteNumber(point, "y");
                if (!points.empty() && (x <= points.back().first ||
                    !std::isfinite(x - points.back().first) || !std::isfinite(y - points.back().second)))
                { fail("PhreaticLine x coordinates must be strictly increasing with finite differences"); }
                points.emplace_back(x, y);
            }
            if (points.size() < 2) { fail("PhreaticLine needs at least two coordinate points"); }
            std::string expression = number(points.front().second);
            for (std::size_t i = 1; i < points.size(); ++i)
            {
                auto const [x0, y0] = points[i-1];auto const [x1, y1] = points[i];
                expression += "+(" + number(y1-y0) + ")*min(max((x-(" + number(x0) + "))/(" + number(x1-x0) + "),0),1)";
            }
            elevations.emplace(id, std::move(expression));
        }
        else { fail("groundwater type must be Hydrostatic or PhreaticLine"); }
    }
    if (!elevations.contains(reference)) { fail("unknown groundwater reference " + reference); }
    if (steady_initial && (clipped.at(reference) || !levels.contains(reference)))
    { fail("SteadyStateSeepage uses a Hydrostatic/Signed initial guess; computed boundary-driven pressure replaces that guess"); }
    std::istringstream gravity_text(
        process.get<std::string>("specific_body_force"));
    double gx, gy;
    std::string extra;
    if (!(gravity_text >> gx >> gy) || (gravity_text >> extra) ||
        !std::isfinite(gx) || !std::isfinite(gy) || gx != 0 || gy > 0 || (gy == 0 && !tsa_placeholder))
    {
        fail("this feature requires 2D vertical gravity: 0 negative_g");
    }
    std::optional<double> density;
    for (auto const& [tag, medium] : root.get_child("media"))
    {
        if (tag != "medium")
        {
            fail("media requires medium entries");
        }
        std::size_t found = 0;
        for (auto const& [phase_tag, phase] : medium.get_child("phases"))
        {
            if (phase_tag != "phase" ||
                phase.get<std::string>("type") != "AqueousLiquid")
            {
                continue;
            }
            for (auto const& [property_tag, property] :
                 phase.get_child("properties"))
            {
                if (property_tag != "property" ||
                    property.get<std::string>("name") != "density")
                {
                    continue;
                }
                ++found;
                double const rho = referenceLiquidDensity(property);
                if (rho <= 0 || (density && *density != rho))
                {
                    fail("positive common liquid density required across media");
                }
                density = rho;
            }
        }
        if (found != 1)
        {
            fail("each medium needs exactly one AqueousLiquid density");
        }
    }
    if (!density || !std::isfinite(*density * -gy) ||
        !std::isfinite(*density * -gy))
    {
        fail("finite liquid unit weight required");
    }
    auto const pressure_variable =
        process.get<std::string>("process_variables.pressure");
    ptree* variable = nullptr;
    for (auto& [tag, candidate] : root.get_child("process_variables"))
    {
        if (tag == "process_variable" &&
            candidate.get<std::string>("name") == pressure_variable)
        {
            if (variable)
            {
                fail("duplicate pressure variable");
            }
            variable = &candidate;
        }
    }
    if (!variable || variable->count("initial_condition"))
    {
        fail("pressure must have exactly one initial owner: remove competing initial_condition");
    }
    // Define each catalogue entry once, then bind initial and boundary uses
    // independently. All references share the same density/gravity conversion.
    auto parameterName = [](std::string const& id)
    { return "__hm_initial_groundwater_" + id; };
    auto& parameters = root.get_child("parameters");
    for (auto const& [id, elevation] : elevations)
    {
        if (levels.contains(id) && !std::isfinite(*density * -gy * levels.at(id)))
        {
            fail("finite hydrostatic pressure offset required for every state");
        }
        auto const name = parameterName(id);
        for (auto const& [tag, parameter] : parameters)
        {
            if (tag == "parameter" &&
                parameter.get<std::string>("name") == name)
            {
                fail("reserved generated parameter collision");
            }
        }
        ptree parameter;
        parameter.put("name", name);
        parameter.put("type", "Function");
        std::string const pressure_expression = number(*density * -gy) + "*((" + elevation + ")-y)";
        parameter.put("expression", (tsa_placeholder || (zero_initial_guess && id == reference)) ? std::string{"0"} :
            ramp_factor + (clipped.at(id) ? "max(0," + pressure_expression + ")" : pressure_expression) + ramp_close);
        parameters.push_back({"parameter", parameter});
    }
    variable->put("initial_condition", parameterName(reference));
    if (fixed_field_mode)
    {
        ptree boundaries;
        if (variable->get<int>("order") != 1)
        { fail("Fixed-field initialisation requires P1 pressure."); }
        if (auto existing=variable->get_child_optional("boundary_conditions"))
        {
            for (auto const& [tag, boundary] : *existing)
            {
                auto const type=boundary.get<std::string>("type", "");
                if (type!="PumpingWell" && type!="SeepageFace")
                { fail("Fixed-field initialisation FixedInitialField owns all pressure constraints; remove stock pressure BCs"); }
                double const start=finiteNumber(boundary,"start_time");
                double const end=finiteNumber(boundary,"end_time");
                if (!phase_schedule.hydraulic_phases ||
                    start<phase_schedule.times.at(1) || end<=start ||
                    end>phase_schedule.times.back())
                { fail("PumpingWell/SeepageFace requires an interval in later hydraulic phases."); }
                for (std::size_t i=0;i<phase_schedule.hydraulic_policies.size();++i)
                {
                    auto const& policy=phase_schedule.hydraulic_policies[i];
                    bool const overlap=start<phase_schedule.times[i+1] && end>phase_schedule.times[i];
                    if (overlap && (policy=="FixedInitialField" || policy=="FrozenAcceptedField"))
                    { fail("PumpingWell/SeepageFace cannot overlap a prescribed interior pressure field."); }
                }
                boundaries.push_back({tag,boundary});
            }
        }
        ptree bc;
        bool const automatic_pressure_nodes = gravity->pressure_mesh.empty();
        bc.put("mesh", automatic_pressure_nodes ? gravity->bulk_mesh : gravity->pressure_mesh);
        bc.put("type", automatic_pressure_nodes ? "InitialFieldDirichlet" : "Dirichlet");
        if (automatic_pressure_nodes && ramp_end > 0.)
        {
            bc.put("evaluate_at_current_time", true);
        }
        bc.put("component", 0);
        bc.put("parameter", parameterName(reference));
        if (phase_schedule.hydraulic_phases)
        {
            std::size_t prefix = 0;
            while (prefix < phase_schedule.hydraulic_policies.size() &&
                   phase_schedule.hydraulic_policies[prefix] == "FixedInitialField")
            {
                ++prefix;
            }
            if (automatic_pressure_nodes)
            {
                bc.put("end_time", number(phase_schedule.times.at(prefix)));
            }
            else
            {
                bc.put("type", "DirichletWithinTimeInterval");
                bc.put("time_interval.start", 0);
                bc.put("time_interval.end", number(phase_schedule.times.at(prefix)));
            }
        }
        if (!steady_initial) { boundaries.push_back({"boundary_condition", bc}); }
        for (std::size_t i=0;i<phase_schedule.hydraulic_policies.size();++i)
        {
            if (phase_schedule.hydraulic_policies[i]!="FrozenAcceptedField") { continue; }
            ptree frozen;
            frozen.put("mesh",gravity->bulk_mesh);frozen.put("type","AcceptedFieldDirichlet");
            frozen.put("component",0);frozen.put("parameter",parameterName(reference));
            frozen.put("capture_time",number(phase_schedule.times.at(i)));
            frozen.put("hold_until",number(phase_schedule.times.at(i+1)));
            boundaries.push_back({"boundary_condition",frozen});
        }
        variable->put_child("boundary_conditions", boundaries);
        auto const displacement_name = process.get<std::string>("process_variables.displacement");
        auto zeroConstant = [&](std::string const& name)
        {
            for (auto const& [tag, parameter] : parameters)
            {
                if (tag != "parameter" || parameter.get<std::string>("name") != name)
                {
                    continue;
                }
                if (parameter.get<std::string>("type") != "Constant")
                {
                    return false;
                }
                std::istringstream values(parameter.get<std::string>("values"));
                double value;
                bool found = false;
                while (values >> value)
                {
                    found = true;
                    if (value != 0 || !std::isfinite(value))
                    {
                        return false;
                    }
                }
                return found && values.eof();
            }
            return false;
        };
        for (auto const& [tag, candidate] : root.get_child("process_variables"))
        {
            if (tag != "process_variable" || candidate.get<std::string>("name") != displacement_name)
            {
                continue;
            }
            if (candidate.get<int>("order") != 2 ||
                !zeroConstant(candidate.get<std::string>("initial_condition")))
            {
                fail("Fixed-field initialisation requires quadratic displacement with zero initial displacement");
            }
            for (auto const& [bc_tag, boundary] : candidate.get_child("boundary_conditions"))
            {
                if (bc_tag != "boundary_condition" ||
                    boundary.get<std::string>("type") != "Dirichlet" ||
                    !zeroConstant(boundary.get<std::string>("parameter")))
                {
                    fail("Fixed-field initialisation allows only fixed zero displacement supports and explicit constant surface water");
                }
            }
        }
    }

    if (fixed_field_mode && method != "Gravity")
    {
        auto const source_name = initial.get<std::string>(method == "K0" ?
            "vertical_effective_stress_parameter" : "stress_parameter");
        ptree const* source = nullptr;
        for (auto const& [tag, parameter] : parameters)
        {
            if (tag == "parameter" && parameter.get<std::string>("name") == source_name)
            {
                if (source)
                {
                    fail("duplicate stress-seed parameter");
                }
                source = &parameter;
            }
        }
        if (!source)
        {
            fail("unknown stress-seed parameter " + source_name);
        }
        std::vector<std::string> expressions;
        auto const type = source->get<std::string>("type");
        if (type == "Constant")
        {
            std::istringstream values(source->get<std::string>("values"));
            std::string token;
            while (values >> token)
            {
                ptree number_tree;
                number_tree.put("value", token);
                expressions.push_back(number(finiteNumber(number_tree, "value")));
            }
        }
        else if (type == "Function")
        {
            for (auto const& [tag, entry] : *source)
            {
                if (tag == "expression")
                {
                    if (std::regex_search(entry.data(),
                        std::regex("(^|[^a-zA-Z0-9_])t([^a-zA-Z0-9_]|$)")))
                    {
                        fail("initial stress expressions must be time independent");
                    }
                    expressions.push_back(entry.data());
                }
            }
        }
        else
        {
            fail("stress seed supports Constant or time-independent Function parameters only");
        }
        if (expressions.size() != (method == "K0" ? 1u : 4u))
        {
            fail("K0 vertical stress needs one component; AssignedStress needs xx yy zz xy");
        }
        auto stress_name = source_name;
        if (method == "K0")
        {
            stress_name = "__hm_k0_initial_stress";
            for (auto const& [tag, parameter] : parameters)
            {
                if (tag == "parameter" && parameter.get<std::string>("name") == stress_name)
                {
                    fail("reserved K0 parameter collision");
                }
            }
            auto const vertical = expressions.front();
            auto const horizontal = number(finiteNumber(initial, "k0")) + "*(" + vertical + ")";
            ptree parameter;
            parameter.put("name", stress_name);
            parameter.put("type", "Function");
            for (auto const& expression : {horizontal, vertical, horizontal, std::string("0")})
            {
                parameter.add("expression", expression);
            }
            parameters.push_back({"parameter", parameter});
        }
        ptree stress;
        stress.put("<xmlattr>.type", initial.get<std::string>("stress_basis") == "Effective" ? "effective" : "total");
        stress.data() = stress_name;
        process.put_child("initial_stress", stress);
    }

    if (process.count("hydraulic_boundaries") > 1)
    {
        fail("duplicate hydraulic_boundaries catalogue");
    }
    if (auto const bindings = process.get_child_optional("hydraulic_boundaries"))
    {
        if (bindings->empty())
        {
            fail("hydraulic_boundaries requires at least one boundary");
        }
        if (!variable->get_child_optional("boundary_conditions"))
        {
            variable->add_child("boundary_conditions", ptree{});
        }
        auto& boundaries = variable->get_child("boundary_conditions");
        for (auto const& [tag, boundary] : *bindings)
        {
            if (tag != "boundary")
            {
                fail("hydraulic_boundaries requires boundary entries");
            }
            keys(boundary, {"mesh", "groundwater"});
            keys(boundary.get_child("mesh"), {});
            auto const mesh = boundary.get<std::string>("mesh");
            auto const& groundwater = boundary.get_child("groundwater");
            keys(groundwater, {"<xmlattr>"});
            if (!groundwater.data().empty())
            {
                fail("groundwater binding accepts a ref attribute only");
            }
            keys(groundwater.get_child("<xmlattr>"), {"ref"});
            auto const state = groundwater.get<std::string>("<xmlattr>.ref");
            if (!elevations.contains(state))
            {
                fail("unknown boundary groundwater reference " + state);
            }
            bool mesh_exists = false;
            for (auto const& [mesh_tag, entry] : root.get_child("meshes"))
            {
                if (mesh_tag == "mesh" &&
                    std::filesystem::path(entry.data()).stem() == mesh)
                {
                    mesh_exists = true;
                }
            }
            if (mesh.empty() || !mesh_exists)
            {
                fail("boundary mesh must name a declared mesh: " + mesh);
            }
            for (auto const& [bc_tag, bc] : boundaries)
            {
                if (bc_tag == "boundary_condition" &&
                    bc.get<std::string>("mesh", "") == mesh)
                {
                    fail("competing pressure boundary on mesh " + mesh);
                }
            }
            ptree bc;
            bc.put("mesh", mesh);
            bc.put("type", "Dirichlet");
            bc.put("component", 0);
            bc.put("parameter", parameterName(state));
            boundaries.push_back({"boundary_condition", bc});
        }
        process.erase("hydraulic_boundaries");
    }
    // Exterior water is a mechanical normal pressure, independent of the
    // groundwater initial field and hydraulic Dirichlet conditions.
    bool const has_surface_states = process.count("surface_water_states") != 0;
    bool const has_surface_loads = process.count("surface_water_loads") != 0;
    if (has_surface_states || has_surface_loads)
    {
        if (process.count("surface_water_states") != 1 ||
            process.count("surface_water_loads") != 1)
        {
            fail("one surface_water_states catalogue and surface_water_loads binding list required");
        }
        std::map<std::string, std::string> surface_parameters;
        for (auto const& [tag, state] : process.get_child("surface_water_states"))
        {
            if (tag != "state")
            {
                fail("surface water catalogue requires state entries");
            }
            keys(state, {"<xmlattr>", "type", "water_level", "pressure_policy"});
            keys(state.get_child("<xmlattr>"), {"id"});
            keys(state.get_child("type"), {});
            keys(state.get_child("water_level"), {});
            keys(state.get_child("pressure_policy"), {});
            auto const id = state.get<std::string>("<xmlattr>.id");
            if (id.empty() || id.find_first_not_of(
                    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_") !=
                    std::string::npos)
            {
                fail("surface water state id must be a nonempty identifier");
            }
            if (state.get<std::string>("type") != "Hydrostatic" ||
                state.get<std::string>("pressure_policy") != "CompressionOnly")
            {
                fail("surface water requires explicit Hydrostatic/CompressionOnly");
            }
            double const level = finiteNumber(state, "water_level");
            if (!std::isfinite(*density * -gy * level))
            {
                fail("finite surface water pressure offset required");
            }
            auto const name = "__hm_surface_water_" + id;
            if (!surface_parameters.emplace(id, name).second)
            {
                fail("duplicate surface water state " + id);
            }
            for (auto const& [parameter_tag, parameter] : parameters)
            {
                if (parameter_tag == "parameter" &&
                    parameter.get<std::string>("name") == name)
                {
                    fail("reserved surface water parameter collision");
                }
            }
            ptree parameter;
            parameter.put("name", name);
            parameter.put("type", "Function");
            parameter.put("expression", ramp_factor + "max(0," + number(*density * -gy) +
                "*((" + number(level) + ")-y))" + ramp_close);
            parameters.push_back({"parameter", parameter});
        }
        auto const displacement_name =
            process.get<std::string>("process_variables.displacement");
        ptree* displacement = nullptr;
        for (auto& [tag, candidate] : root.get_child("process_variables"))
        {
            if (tag == "process_variable" &&
                candidate.get<std::string>("name") == displacement_name)
            {
                if (displacement)
                {
                    fail("duplicate displacement variable");
                }
                displacement = &candidate;
            }
        }
        if (!displacement || displacement->get<int>("components") != 2)
        {
            fail("surface water requires the 2D displacement variable");
        }
        if (!displacement->get_child_optional("boundary_conditions"))
        {
            displacement->add_child("boundary_conditions", ptree{});
        }
        auto& mechanical_boundaries = displacement->get_child("boundary_conditions");
        auto const& loads = process.get_child("surface_water_loads");
        if (loads.empty())
        {
            fail("surface_water_loads requires at least one load");
        }
        for (auto const& [tag, load] : loads)
        {
            if (tag != "load")
            {
                fail("surface_water_loads requires load entries");
            }
            keys(load, {"mesh", "surface_water"});
            keys(load.get_child("mesh"), {});
            auto const mesh = load.get<std::string>("mesh");
            bool exists = false;
            for (auto const& [mesh_tag, entry] : root.get_child("meshes"))
            {
                if (mesh_tag == "mesh" &&
                    std::filesystem::path(entry.data()).stem() == mesh)
                {
                    exists = true;
                }
            }
            if (mesh.empty() || !exists)
            {
                fail("surface water mesh must name a declared mesh: " + mesh);
            }
            auto const& binding = load.get_child("surface_water");
            keys(binding, {"<xmlattr>"});
            keys(binding.get_child("<xmlattr>"), {"ref"});
            if (!binding.data().empty())
            {
                fail("surface_water binding accepts a ref attribute only");
            }
            auto const state = binding.get<std::string>("<xmlattr>.ref");
            if (!surface_parameters.contains(state))
            {
                fail("unknown surface water reference " + state);
            }
            for (auto const& [bc_tag, bc] : mechanical_boundaries)
            {
                if (bc_tag == "boundary_condition" &&
                    bc.get<std::string>("mesh", "") == mesh)
                {
                    fail("competing mechanical boundary on surface water mesh " + mesh);
                }
            }
            ptree bc;
            bc.put("mesh", mesh);
            bc.put("type", "NormalTraction");
            bc.put("parameter", surface_parameters.at(state));
            mechanical_boundaries.push_back({"boundary_condition", bc});
        }
        process.erase("surface_water_states");
        process.erase("surface_water_loads");
    }
    if (phase_schedule.hydraulic_phases && !fixed_field_mode)
    {
        fail("hydraulic phases require an explicit initialisation method");
    }
    for (auto& [tag,variable]:root.get_child("process_variables"))
    {
        if (variable.get<std::string>("name","")==process.get<std::string>("process_variables.displacement"))
        {
            for (auto const& boundary:excavation_water_tractions)
            { variable.get_child("boundary_conditions").push_back({"boundary_condition",boundary}); }
        }
    }
    bindEngineeringLoads(root, process, phase_schedule);
    bindEngineeringHydraulics(root, process, phase_schedule,
                             levels.contains(reference) && !clipped.at(reference) ? levels.at(reference) : std::numeric_limits<double>::quiet_NaN(), *density * -gy);
    process.erase("groundwater_states");
    process.erase("initialisation");
    cropCheckpointTimeStepping(root,checkpoint_schedule);
    input.str("");
    input.clear();
    boost::property_tree::write_xml(input, document);
    std::filesystem::path const directory =
        output_directory.empty() ? "." : output_directory;
    std::filesystem::create_directories(directory);
    boost::property_tree::write_xml(
        (directory / "hm_initial_groundwater_expanded.prj").string(),
        document);
    return true;
}
}  // namespace ProcessLib::HydroMechanics
