// SPDX-License-Identifier: BSD-3-Clause
#pragma once

namespace ProcessLib::HydroMechanics
{
/// Lower only explicit native policies; the internal interval key is reserved.
inline void prepareEngineeringSeepage(
    boost::property_tree::ptree const& root,
    boost::property_tree::ptree& process,
    EngineeringPhaseSchedule const& schedule)
{
    using namespace InitialGroundwaterDetail;
    ptree transient_intervals;
    for (std::size_t i = 0; i < schedule.hydraulic_policies.size(); ++i)
    {
        if (schedule.hydraulic_policies[i] != "TransientSeepage") { continue; }
        if (process.get<std::string>("coupling_scheme", "monolithic") != "monolithic" ||
            process.get<bool>("linear", false) || process.count("construction") ||
            process.count("construction_births") || process.count("deactivated_subdomains"))
        { fail("TransientSeepage requires nonlinear monolithic fixed-mesh HM"); }
        ptree interval;
        interval.put("start", number(schedule.times.at(i)));
        interval.put("end", number(schedule.times.at(i + 1)));
        transient_intervals.push_back({"interval", interval});
    }
    if (!transient_intervals.empty())
    { process.put_child("transient_seepage", transient_intervals); }
    ptree intervals;
    for (std::size_t i = 0; i < schedule.hydraulic_policies.size(); ++i)
    {
        if (schedule.hydraulic_policies[i] != "SteadyStateSeepage") { continue; }
        if (!std::any_of(schedule.heads.begin(), schedule.heads.end(),
                         [i](auto const& item) { return item.second.active.at(i); }))
        {
            fail("SteadyStateSeepage needs an active prescribed-head boundary; an all-Neumann pressure problem has no unique datum.");
        }
        ptree interval;
        interval.put("start", number(schedule.times.at(i)));
        interval.put("end", number(schedule.times.at(i + 1)));
        intervals.push_back({"interval", interval});
    }
    if (intervals.empty() && transient_intervals.empty()) { return; }
    // Effective-stress drainage may start at a later accepted phase. Its
    // lowering validates that every steady interval ends by that reference.
    if (auto materials = process.get_child_optional("mechanical_drainage"))
    {
        for (auto const& [tag, material] : *materials)
        { if (tag == "material" && material.get<std::string>("type") == "UndrainedTotalStress")
          { fail("Steady seepage cannot initialise a total-stress TSA basis"); } }
    }
    auto const method=process.get<std::string>("initialisation.method", "");
    bool const initial_steady=!schedule.hydraulic_policies.empty() &&
        schedule.hydraulic_policies.front()=="SteadyStateSeepage";
    if (process.get<std::string>("coupling_scheme", "monolithic") != "monolithic" ||
        process.get<bool>("linear", false) ||
        process.count("construction") || process.count("construction_births") ||
        process.count("deactivated_subdomains") ||
        (method!="Gravity" && method!="K0" && method!="AssignedStress") ||
        (initial_steady && method!="Gravity"))
    {
        fail("SteadyStateSeepage requires nonlinear monolithic HM; initial seepage requires Gravity, later seepage permits Gravity/K0/AssignedStress. Construction birth is unsupported.");
    }
    for (auto const& [tag, law] : process)
    {
        if (tag == "constitutive_relation" &&
            (law.get<std::string>("type", "") != "MFront" ||
             (law.get<std::string>("behaviour", "") != "MohrCoulombAbboSloan" &&
              law.get<std::string>("behaviour", "") != "MohrCoulombAbboSloanTension" &&
              law.get<std::string>("behaviour", "") != "PublishedHS")))
        { fail("SteadyStateSeepage currently requires MFront MC or PublishedHS."); }
        if (tag == "constitutive_relation" &&
            law.get<std::string>("behaviour", "") == "PublishedHS" &&
            !law.count("substepping"))
        { fail("PublishedHS seepage requires explicit transactional substepping."); }
    }
    // Steady seepage uses the liquid reference density at zero gauge pressure,
    // including when a pressure-Linear EOS is configured for later Coupled
    // phases. Constant transport restores exact p(lambda)=lambda*p(full).
    auto constant = [&](auto const& self, ptree const& node, bool liquid) -> void
    {
        for (auto const& [tag, child] : node)
        {
            if (tag == "property" && child.get<std::string>("type", "") != "Constant")
            {
                if (!(liquid && child.get<std::string>("name", "")=="density" && child.get<std::string>("type", "")=="Linear"))
                { fail("SteadyStateSeepage requires Constant transport properties except pressure-Linear liquid density."); }
                if (!(referenceLiquidDensity(child)>0.)) { fail("positive liquid reference density required"); }
            }
            self(self, child, tag=="phase"?child.template get<std::string>("type", "")=="AqueousLiquid":liquid);
        }
    };
    constant(constant, root.get_child("media"), false);
    if (!intervals.empty()) { process.put_child("steady_state_seepage", intervals); }
}
}  // namespace ProcessLib::HydroMechanics
