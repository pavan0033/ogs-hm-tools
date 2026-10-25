// SPDX-License-Identifier: BSD-3-Clause
#pragma once

namespace ProcessLib::HydroMechanics
{
inline void bindEngineeringHydraulics(
    boost::property_tree::ptree& root,
    boost::property_tree::ptree& process,
    EngineeringPhaseSchedule const& schedule,
    double const initial_level, double const water_weight)
{
    using namespace InitialGroundwaterDetail;
    if (!schedule.hydraulic_phases)
    {
        return;
    }
    if (process.count("hydraulic_boundaries") || process.count("surface_water_states") ||
        process.count("surface_water_loads"))
    {
        fail("hydraulic_controls owns phase heads and companion surface water; remove competing catalogues");
    }
    auto& parameters = root.get_child("parameters");
    auto const pressure = process.get<std::string>("process_variables.pressure");
    auto const displacement = process.get<std::string>("process_variables.displacement");
    for (auto const& [id, head] : schedule.heads)
    {
        if (!head.seepage && !process.get_child_optional("steady_state_seepage") && std::isfinite(initial_level) && std::abs(head.values.front() - initial_level) >
            1e-12 * std::max(1., std::abs(initial_level)))
        {
            fail("hydraulic control initial_head must match the initial hydrostatic groundwater state");
        }
        auto const name = "__hm_hyd_" + id;
        for (auto const& [tag, parameter] : parameters)
        {
            if (parameter.get<std::string>("name") == name + "_pressure" ||
                parameter.get<std::string>("name") == name + "_water")
            {
                fail("reserved hydraulic parameter collision");
            }
        }
        std::string level_expression = number(head.values.front());
        for (std::size_t i = 0; i + 1 < schedule.times.size(); ++i)
        {
            if (!std::isfinite(water_weight * head.values.at(i + 1)))
            {
                fail("hydraulic head yields nonfinite pressure");
            }
            double const change = head.values[i + 1] - head.values[i];
            if (change != 0.)
            {
                level_expression += "+(" + number(change) + ")*min(max((t-(" +
                    number(schedule.times[i]) + "))/(" +
                    number(schedule.times[i + 1] - schedule.times[i]) + "),0),1)";
            }
        }
        auto const expression = number(water_weight) + "*((" + level_expression + ")-y)";
        ptree p;
        p.put("name", name + "_pressure");
        p.put("type", "Function");
        bool const steady_initial = process.get<std::string>("initialisation.hydraulic_policy", "") == "SteadyStateSeepage";
        double const ramp_end = process.get<double>("initialisation_load_ramp_end", 0.);
        p.put("expression", steady_initial ?
            "min(max(t/(" + number(ramp_end) + "),0),1)*(" + expression + ")" : expression);
        parameters.push_back({"parameter", p});
        if (head.surface_water)
        {
            p.put("name", name + "_water");
            std::string water_expression="max(0," + expression + ")";
            double const ramp_end=process.get<double>("initialisation_load_ramp_end",0.);
            if (ramp_end>0.)
            {
                water_expression="min(max(t/(" + number(ramp_end) + "),0),1)*(" + water_expression + ")";
            }
            p.put("expression",water_expression);
            parameters.push_back({"parameter", p});
        }
        for (auto& [tag, variable] : root.get_child("process_variables"))
        {
            if (tag != "process_variable")
            {
                continue;
            }
            auto const variable_name = variable.get<std::string>("name");
            if (variable_name == pressure)
            {
                for (std::size_t i = 0; i < schedule.hydraulic_policies.size(); ++i)
                {
                    if (!head.active.at(i))
                    {
                        continue;
                    }
                    ptree bc;
                    bc.put("mesh", head.pressure_mesh.empty()?head.mesh:head.pressure_mesh);
                    bc.put("type", head.seepage ? "ReservoirSeepage" : "PhaseDirichlet");
                    bc.put("component", 0);
                    bc.put("parameter", name + "_pressure");
                    bc.put("start_time", number(schedule.times[i]));
                    bc.put("end_time", number(schedule.times[i + 1]));
                    variable.get_child("boundary_conditions").push_back({"boundary_condition", bc});
                }
            }
            if (variable_name == displacement && head.surface_water)
            {
                auto& boundaries = variable.get_child("boundary_conditions");
                for (auto const& [bc_tag, bc] : boundaries)
                {
                    if (bc_tag == "boundary_condition" && bc.get<std::string>("mesh") == head.mesh)
                    {
                        // A declared surcharge and exterior water pressure are
                        // additive natural loads, not competing constraints.
                        bool additive=false;
                        for (auto const& [load_id,load] : schedule.loads)
                        {
                            additive |= load.type=="NormalPressure" && load.mesh==head.mesh &&
                                bc.get<std::string>("type")=="NormalTraction" &&
                                bc.get<std::string>("parameter","")=="__hm_phase_"+load_id;
                        }
                        if (!additive)
                        { fail("companion water load has competing mechanical condition on " + head.mesh); }
                    }
                }
                ptree bc;
                bc.put("mesh", head.mesh);
                bc.put("type", "NormalTraction");
                bc.put("parameter", name + "_water");
                if (head.water_at_ips) { bc.put("parameter_evaluation","IntegrationPoints"); }
                boundaries.push_back({"boundary_condition", bc});
            }
        }
    }
    process.erase("hydraulic_controls");
    process.erase("phases");
}
}  // namespace ProcessLib::HydroMechanics
