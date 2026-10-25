// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <boost/property_tree/ptree.hpp>
#include <cmath>
#include <iomanip>
#include <map>
#include <set>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ProcessLib::HydroMechanics
{
struct EngineeringLoad
{
    std::string mesh;
    std::vector<double> values{0.};
    std::string type = "NormalPressure";
    int component = -1;
    std::string spatial_scale_parameter;
};

struct EngineeringHead
{
    std::string mesh;
    bool surface_water;
    std::vector<double> values;
    std::vector<bool> active;
    bool seepage=false;
    std::string pressure_mesh;
    bool water_at_ips=false;
};

struct EngineeringPhaseSchedule
{
    bool hydraulic_phases = false;
    std::vector<std::string> hydraulic_policies;
    std::map<std::string, EngineeringHead> heads;
    std::vector<double> times{0.};
    std::map<std::string, EngineeringLoad> loads;
};

inline EngineeringPhaseSchedule prepareEngineeringPhases(
    boost::property_tree::ptree& root,
    boost::property_tree::ptree& process,
    bool const allow_empty_mechanical_loads = false)
{
    using namespace InitialGroundwaterDetail;
    EngineeringPhaseSchedule schedule;
    if (!process.count("phases") && !process.count("mechanical_loads") && !process.count("hydraulic_controls"))
    {
        return schedule;
    }
    schedule.hydraulic_phases = process.count("hydraulic_controls") != 0;
    if (process.count("phases") != 1 || (!schedule.hydraulic_phases && !allow_empty_mechanical_loads && process.count("mechanical_loads") != 1) || process.count("mechanical_loads") > 1 ||
        !process.get_child_optional("initialisation.method"))
    {
        fail("phases require one mechanical_loads catalogue and explicit initialisation method");
    }
    std::set<std::string> meshes;
    for (auto const& [tag, entry] : root.get_child("meshes"))
    {
        if (tag == "mesh")
        {
            meshes.insert(std::filesystem::path(entry.data()).stem().string());
        }
    }
    if (schedule.hydraulic_phases)
    {
        std::set<std::string> controlled_meshes;
        for (auto const& [tag, control] : process.get_child("hydraulic_controls"))
        {
            if (tag != "boundary")
            {
                fail("hydraulic_controls requires boundary entries");
            }
            keys(control, {"<xmlattr>", "type", "mesh", "initial_head", "pressure_policy", "apply_surface_water", "automatic_pressure_mesh", "water_at_integration_points"});
            keys(control.get_child("<xmlattr>"), {"id"});
            for (auto const key : {"type", "mesh", "initial_head", "pressure_policy", "apply_surface_water"})
            {
                keys(control.get_child(key), {});
            }
            auto const id = control.get<std::string>("<xmlattr>.id");
            auto const mesh = control.get<std::string>("mesh");
            auto const water = control.get<std::string>("apply_surface_water");
            if (id.empty() || id.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_") != std::string::npos ||
                !meshes.contains(mesh) || !controlled_meshes.insert(mesh).second ||
                (control.get<std::string>("type") != "HydrostaticHead" && control.get<std::string>("type") != "ReservoirSeepage") ||
                control.get<std::string>("pressure_policy") != "Signed" ||
                (water != "true" && water != "false"))
            {
                fail("invalid hydraulic control ID/type/mesh/Signed policy/surface-water flag");
            }
            double const level = finiteNumber(control, "initial_head");
            if (!schedule.heads.emplace(id, EngineeringHead{mesh, water == "true", {level}, {}}).second)
            {
                fail("duplicate hydraulic control ID");
            }
            schedule.heads.at(id).pressure_mesh=control.get<std::string>("automatic_pressure_mesh",mesh);
            if (!meshes.contains(schedule.heads.at(id).pressure_mesh)) { fail("unknown automatic pressure mesh"); }
            schedule.heads.at(id).water_at_ips=control.get<bool>("water_at_integration_points",false);
            schedule.heads.at(id).seepage=control.get<std::string>("type")=="ReservoirSeepage";
        }
    }
    auto const mechanical = process.get_child_optional("mechanical_loads");
    boost::property_tree::ptree const empty;
    for (auto const& [tag, load] : (mechanical ? *mechanical : empty))
    {
        if (tag != "load")
        {
            fail("mechanical_loads requires load entries");
        }
        keys(load, {"<xmlattr>", "type", "mesh", "component", "reference", "spatial_scale_parameter"});
        keys(load.get_child("<xmlattr>"), {"id"});
        keys(load.get_child("type"), {});
        keys(load.get_child("mesh"), {});
        auto const id = load.get<std::string>("<xmlattr>.id");
        auto const mesh = load.get<std::string>("mesh");
        auto const type = load.get<std::string>("type");
        bool const displacement = type == "PrescribedDisplacement";
        int component = -1;
        if (displacement)
        {
            keys(load.get_child("component"), {});
            keys(load.get_child("reference"), {});
            component = load.get<int>("component");
            if ((component != 0 && component != 1) ||
                load.get<std::string>("reference") != "InitialisationEnd")
            {
                fail("PrescribedDisplacement requires component 0/1 and reference InitialisationEnd");
            }
        }
        else if (load.count("component") || load.count("reference"))
        {
            fail("NormalPressure must not specify displacement component/reference");
        }
        if (id.empty() || id.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_") != std::string::npos ||
            (type != "NormalPressure" && !displacement) || !meshes.contains(mesh))
        {
            fail("mechanical load requires unique identifier, declared mesh and NormalPressure/PrescribedDisplacement type");
        }
        std::string const spatial=load.get<std::string>("spatial_scale_parameter","");
        if(!spatial.empty())
        {
            if(displacement) { fail("spatial_scale_parameter currently applies only to NormalPressure"); }
            bool found=false;
            for(auto const& [tag,p]:root.get_child("parameters"))
            {
                if(tag!="parameter" || p.get<std::string>("name")!=spatial) { continue; }
                auto const kind=p.get<std::string>("type");
                if(kind=="Function")
                {
                    if(p.count("expression")!=1 || std::regex_search(p.get<std::string>("expression"),std::regex(R"(\bt\b)")))
                    { fail("spatial pressure shape must be a scalar time-independent Function"); }
                }
                else if(kind=="Constant")
                {
                    std::istringstream values(p.get<std::string>("values"));double v;std::string extra;
                    if(!(values>>v) || !std::isfinite(v) || v<0. || (values>>extra))
                    { fail("spatial pressure shape must be scalar finite and nonnegative"); }
                }
                else { fail("spatial pressure shape supports only Constant or static Function"); }
                found=true;
            }
            if(!found) { fail("unknown spatial_scale_parameter"); }
        }
        if (!schedule.loads.emplace(id, EngineeringLoad{mesh, {0.}, type, component, spatial}).second)
        {
            fail("duplicate mechanical load " + id);
        }
    }
    if (schedule.loads.empty() && !schedule.hydraulic_phases && !allow_empty_mechanical_loads)
    {
        fail("mechanical_loads must not be empty");
    }
    auto& time_process = root.get_child("time_loop.processes.process");
    if (time_process.count("time_stepping"))
    {
        fail("phases own time_stepping; remove the competing stock schedule");
    }
    boost::property_tree::ptree stepping;
    stepping.put("type", "FixedTimeStepping");
    stepping.put("t_initial", 0);
    bool adaptive_schedule=false; boost::property_tree::ptree intervals;
    std::set<std::string> phase_ids;
    std::size_t index = 0;
    for (auto const& [tag, phase] : process.get_child("phases"))
    {
        if (tag != "phase")
        {
            fail("phases requires ordered phase entries");
        }
        keys(phase, {"<xmlattr>", "type", "duration", "steps", "stepping", "load", "hydraulic_policy", "heads", "active_heads", "reset_displacement_reference"}, {"load"});
        keys(phase.get_child("<xmlattr>"), {"id"});
        keys(phase.get_child("type"), {});
        keys(phase.get_child("duration"), {});
        if (phase.count("steps")) { keys(phase.get_child("steps"), {}); }
        auto const id = phase.get<std::string>("<xmlattr>.id");
        auto const type = phase.get<std::string>("type");
        if (id.empty() || !phase_ids.insert(id).second ||
            (index == 0 ? type != "Initialisation" : (type != "MechanicalLoading" && !(schedule.hydraulic_phases && type == "Loading"))))
        {
            fail("first phase must be Initialisation; subsequent phases MechanicalLoading; IDs unique");
        }
        if (schedule.hydraulic_phases)
        {
            if (phase.count("hydraulic_policy"))
            {
                keys(phase.get_child("hydraulic_policy"), {});
            }
            // Only the initial phase inherits its declared initialisation policy.
            // Working phases must still state their hydraulic calculation type.
            auto const policy = index == 0
                ? phase.get<std::string>("hydraulic_policy", process.get<std::string>("initialisation.hydraulic_policy"))
                : phase.get<std::string>("hydraulic_policy");
            if ((policy != "FixedInitialField" && policy != "Coupled" && policy != "TransientSeepage" && policy != "Undrained" && policy != "FrozenAcceptedField" && policy != "SteadyStateSeepage") ||
                (index == 0 && policy != "FixedInitialField" && policy != "SteadyStateSeepage"))
            {
                fail("initial hydraulic phase must be FixedInitialField or SteadyStateSeepage; later policies also include Coupled/Undrained/FrozenAcceptedField");
            }
            if (policy == "FixedInitialField")
            {
                for (auto const& previous : schedule.hydraulic_policies)
                {
                    if (previous != "FixedInitialField")
                    {
                        fail("cannot re-enter FixedInitialField after pressure release");
                    }
                }
            }
            if (policy == "FrozenAcceptedField" && process.count("mechanical_drainage"))
            {
                // Freeze the primary reference only. Effective-stress material
                // laws still calculate their actual IP pressure independently.
                for (auto const& [tag, material] : process.get_child("mechanical_drainage"))
                { if (tag == "material" && material.get<std::string>("type") == "UndrainedTotalStress")
                  { fail("FrozenAcceptedField cannot prescribe an actual fluid state for TSA materials."); } }
            }
            schedule.hydraulic_policies.push_back(policy);
            for (auto& [id, head] : schedule.heads)
            {
                head.values.push_back(head.values.back());
                head.active.push_back(policy == "Coupled" || policy == "TransientSeepage" || policy == "SteadyStateSeepage");
            }
            if (auto const selected = phase.get_child_optional("active_heads"))
            {
                if (!selected->data().empty())
                {
                    fail("active_heads requires head entries, not text");
                }
                if (policy != "Coupled" && policy != "TransientSeepage" && policy != "SteadyStateSeepage")
                {
                    fail("active_heads requires Coupled or SteadyStateSeepage");
                }
                std::set<std::string> activated;
                for (auto& [name, head] : schedule.heads)
                {
                    head.active.back() = false;
                }
                for (auto const& [tag, selection] : *selected)
                {
                    if (tag != "head")
                    {
                        fail("active_heads requires head entries");
                    }
                    keys(selection, {"<xmlattr>"});
                    keys(selection.get_child("<xmlattr>"), {"ref"});
                    auto const ref = selection.get<std::string>("<xmlattr>.ref");
                    if (!selection.data().empty() || !schedule.heads.contains(ref) ||
                        !activated.insert(ref).second)
                    {
                        fail("active_heads requires known unique head references");
                    }
                    schedule.heads.at(ref).active.back() = true;
                }
            }
            if (auto const actions = phase.get_child_optional("heads"))
            {
                if (index == 0 || (policy != "Coupled" && policy != "TransientSeepage" && policy != "SteadyStateSeepage"))
                {
                    fail("head changes require a working Coupled or SteadyStateSeepage phase");
                }
                std::set<std::string> changed;
                for (auto const& [tag, action] : *actions)
                {
                    if (tag != "head")
                    {
                        fail("heads requires head entries");
                    }
                    keys(action, {"<xmlattr>", "target_head"});
                    keys(action.get_child("<xmlattr>"), {"ref"});
                    keys(action.get_child("target_head"), {});
                    auto const ref = action.get<std::string>("<xmlattr>.ref");
                    if (!schedule.heads.contains(ref) || !changed.insert(ref).second)
                    {
                        fail("phase head requires known unique reference");
                    }
                    if (!schedule.heads.at(ref).active.back())
                    {
                        fail("cannot change target_head for an inactive hydraulic boundary");
                    }
                    schedule.heads.at(ref).values.back() = finiteNumber(action, "target_head");
                }
            }
        }
        else if (phase.count("hydraulic_policy") || phase.count("heads") || phase.count("active_heads"))
        {
            fail("phase hydraulic actions require explicit hydraulic_controls catalogue");
        }
        double const duration = finiteNumber(phase, "duration");
        if ((phase.count("steps")+phase.count("stepping"))!=1)
        { fail("Choose exactly one of fixed steps or adaptive stepping per phase."); }
        bool const adaptive=phase.count("stepping")!=0;
        int const steps = phase.get<int>("steps",40);
        if (duration<=0. || steps<(index==0?2:1)) { fail("Positive duration and valid fixed step count required."); }
        ptree interval;interval.put("start",number(schedule.times.back()));interval.put("end",number(schedule.times.back()+duration));interval.put("adaptive",adaptive);
        double initial=duration/steps,minimum=initial,maximum=initial;
        if (adaptive)
        {
            auto const& a=phase.get_child("stepping");keys(a,{"type","initial_increment","minimum_increment","maximum_increment","grow_below_iterations","reduce_above_iterations"});
            if (a.get<std::string>("type")!="Adaptive") { fail("stepping requires Adaptive."); }
            double const lo=finiteNumber(a,"minimum_increment"),hi=finiteNumber(a,"maximum_increment"),first=finiteNumber(a,"initial_increment");
            if (!(0.<lo && lo<=first && first<=hi && hi<=1.)) { fail("Adaptive increment fractions require 0<minimum<=initial<=maximum<=1."); }
            int const grow=a.get<int>("grow_below_iterations",4),reduce=a.get<int>("reduce_above_iterations",8);
            if (grow<0 || reduce<=grow) { fail("Adaptive effort range requires 0<=grow_below_iterations<reduce_above_iterations"); }
            interval.put("grow_below_iterations",grow);interval.put("reduce_above_iterations",reduce);
            initial=duration*first;minimum=duration*lo;maximum=duration*hi;adaptive_schedule=true;
        }
        interval.put("initial",number(initial));interval.put("minimum",number(minimum));interval.put("maximum",number(maximum));intervals.push_back({"interval",interval});
        if (duration <= 0 || steps < (index == 0 ? 2 : 1) ||
            !std::isfinite(schedule.times.back() + duration))
        {
            fail("phase duration/steps must be positive; initialisation needs equilibrium plus hold");
        }
        for (auto& [name, load] : schedule.loads)
        {
            load.values.push_back(load.values.back());
        }
        std::set<std::string> changed_loads;
        for (auto const& [tag, action] : phase)
        {
            if (tag != "load") { continue; }
            if (index == 0)
            {
                fail("initialisation phase must not apply a working load");
            }
            keys(action, {"<xmlattr>", "target_pressure", "target_displacement"});
            keys(action.get_child("<xmlattr>"), {"ref"});
            auto const reference = action.get<std::string>("<xmlattr>.ref");
            if (!schedule.loads.contains(reference))
            {
                fail("phase load requires a known reference");
            }
            if (!changed_loads.insert(reference).second)
            {
                fail("duplicate mechanical load action for " + reference);
            }
            auto& load = schedule.loads.at(reference);
            bool const displacement = load.type == "PrescribedDisplacement";
            auto const target = displacement ? "target_displacement" : "target_pressure";
            auto const other = displacement ? "target_pressure" : "target_displacement";
            if (action.count(other))
            {
                fail("phase target does not match mechanical load type");
            }
            keys(action.get_child(target), {});
            double const value = finiteNumber(action, target);
            if (!displacement && value < 0)
            {
                fail("NormalPressure requires nonnegative target_pressure");
            }
            load.values.back() = value;
        }
        if (phase.count("reset_displacement_reference"))
        {
            keys(phase.get_child("reset_displacement_reference"), {});
            auto const reset = phase.get<std::string>("reset_displacement_reference");
            if (index == 0 || (reset != "true" && reset != "false"))
            {
                fail("reset_displacement_reference requires true/false on a working phase");
            }
            if (reset == "true")
            {
                process.add("displacement_reference_times.time", number(schedule.times.back()));
            }
        }
        schedule.times.push_back(schedule.times.back() + duration);
        boost::property_tree::ptree pair;
        pair.put("repeat", steps);
        pair.put("delta_t", number(duration / steps));
        stepping.add_child("timesteps.pair", pair);
        ++index;
    }
    bool const standalone_steady = index == 1 && schedule.hydraulic_phases &&
        schedule.hydraulic_policies.front() == "SteadyStateSeepage";
    if (index < 2 && !standalone_steady)
    {
        fail("phases require initialisation followed by a working phase, or one standalone SteadyStateSeepage initialisation");
    }
    if (adaptive_schedule)
    {
        if (process.get<std::string>("coupling_scheme","monolithic")!="monolithic" || process.get<bool>("linear",false))
        { fail("Adaptive native phases require nonlinear monolithic HM."); }
        for (auto const& [tag, law] : process)
        {
            if (tag != "constitutive_relation") { continue; }
            auto const behaviour = law.get<std::string>("behaviour", "");
            if (law.get<std::string>("type", "") != "MFront" ||
                (behaviour != "MohrCoulombAbboSloan" &&
                 behaviour != "MohrCoulombAbboSloanTension" &&
                 behaviour != "PublishedHS"))
            { fail("Adaptive rollback currently requires MFront MC or PublishedHS."); }
            if (behaviour == "PublishedHS" && !law.count("substepping"))
            { fail("PublishedHS adaptive phases require explicit transactional substepping."); }
        }
        stepping.erase("timesteps");stepping.put("type","PhaseAdaptive");
        for (auto const& entry:intervals) { stepping.push_back(entry); }
        process.put("adaptive_rollback",true);
        auto& output=root.get_child("time_loop.output");auto suffix=output.get<std::string>("suffix","_ts_{:timestep}_t_{:time}");
        if (suffix.find(":converged")==std::string::npos) { suffix+="{:converged}"; }output.put("suffix",suffix);
    }
    stepping.put("t_end", number(schedule.times.back()));
    time_process.add_child("time_stepping", stepping);
    return schedule;
}

inline void bindEngineeringLoads(boost::property_tree::ptree& root,
                                 boost::property_tree::ptree& process,
                                 EngineeringPhaseSchedule const& schedule)
{
    using namespace InitialGroundwaterDetail;
    if (schedule.loads.empty())
    {
        if (!schedule.hydraulic_phases && schedule.times.size() > 1)
        {
            process.erase("mechanical_loads");
            process.erase("phases");
        }
        return;
    }
    auto& parameters = root.get_child("parameters");
    if (!root.get_child_optional("curves"))
    {
        root.add_child("curves", boost::property_tree::ptree{});
    }
    auto& curves = root.get_child("curves");
    auto const displacement_name = process.get<std::string>("process_variables.displacement");
    for (auto& [tag, variable] : root.get_child("process_variables"))
    {
        if (tag != "process_variable" || variable.get<std::string>("name") != displacement_name)
        {
            continue;
        }
        auto& boundaries = variable.get_child("boundary_conditions");
        for (auto const& [id, load] : schedule.loads)
        {
            for (auto const& [bc_tag, bc] : boundaries)
            {
                if (bc_tag == "boundary_condition" && bc.get<std::string>("mesh") == load.mesh)
                {
                    fail("phase load has competing mechanical condition on " + load.mesh);
                }
            }
            auto const name = "__hm_phase_" + id;
            for (auto const& [tag, parameter] : parameters)
            {
                if (parameter.get<std::string>("name") == name || parameter.get<std::string>("name") == name + "_unit")
                {
                    fail("reserved phase parameter collision");
                }
            }
            for (auto const& [tag, curve] : curves)
            {
                if (curve.get<std::string>("name") == name)
                {
                    fail("reserved phase curve collision");
                }
            }
            std::string coordinates, values;
            for (std::size_t i = 0; i < schedule.times.size(); ++i)
            {
                coordinates += number(schedule.times[i]) + " ";
                values += number(load.values[i]) + " ";
            }
            boost::property_tree::ptree curve;
            curve.put("name", name);
            curve.put("coords", coordinates);
            curve.put("values", values);
            curves.push_back({"curve", curve});
            boost::property_tree::ptree unit;
            unit.put("name", name + "_unit");
            unit.put("type", "Constant");
            unit.put("values", 1);
            parameters.push_back({"parameter", unit});
            boost::property_tree::ptree parameter;
            parameter.put("name", name);
            parameter.put("type", "CurveScaled");
            parameter.put("curve", name);
            parameter.put("parameter", load.spatial_scale_parameter.empty()?name + "_unit":load.spatial_scale_parameter);
            parameters.push_back({"parameter", parameter});
            boost::property_tree::ptree bc;
            bc.put("mesh", load.mesh);
            bc.put("type", load.type == "PrescribedDisplacement" ?
                               "RelativeDirichlet" : "NormalTraction");
            if (load.type == "PrescribedDisplacement")
            {
                bc.put("component", load.component);
                bc.put("reference_time", number(schedule.times.at(1)));
            }
            bc.put("parameter", name);
            boundaries.push_back({"boundary_condition", bc});
        }
    }
    process.erase("mechanical_loads");
    process.erase("phases");
}
}  // namespace ProcessLib::HydroMechanics
