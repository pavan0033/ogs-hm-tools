// SPDX-License-Identifier: BSD-3-Clause
#pragma once
namespace ProcessLib::HydroMechanics
{
inline bool prepareEngineeringConstruction(boost::property_tree::ptree& root,
                                           boost::property_tree::ptree& process)
{
    using namespace InitialGroundwaterDetail;
    if (!process.count("phases")) { return false; }
    bool requested=false;
    for (auto const& [tag,p]:process.get_child("phases")) { requested |= p.count("activate")!=0; }
    if (!requested) { return false; }
    if (process.count("checkpoints") || process.count("restart") || process.count("mechanical_drainage") ||
        process.get<std::string>("coupling_scheme","monolithic")!="monolithic" ||
        process.get_child_optional("initialisation.equilibrium_check"))
    { fail("Construction first scope requires ordinary monolithic native fixed-step HM without checkpoint, strict drainage, initialisation gate."); }
    double end=0.;for (auto const& [tag,p]:process.get_child("phases")) { end+=finiteNumber(p,"duration"); }
    struct Region { double inactive_start=0.; bool active=false; bool born=false; };
    std::map<int, Region> regions;
    double start=0.;ptree masks,births;
    for (auto& [tag,p]:process.get_child("phases"))
    {
        for (auto const action_name:{"activate","excavate"})
        {
            auto action=p.get_child_optional(action_name);if (!action) { continue; }
            if (start<=0. || p.get<std::string>("type","")!="Loading" || p.get<std::string>("hydraulic_policy","")!="Coupled")
            { fail("Construction actions require later coupled Loading phases."); }
            if (std::string(action_name)=="activate")
            {
                keys(*action,{"material_id","state_policy"});int const id=action->get<int>("material_id");
                if (id<0 || regions[id].active) { fail("Construction activation needs a nonnegative inactive material ID."); }
                auto& region=regions[id];auto const policy=action->get<std::string>("state_policy");
                if ((policy!="StressFreeBirth" && policy!="PreserveHistory") || (!region.born && policy!="StressFreeBirth"))
                { fail("First activation needs StressFreeBirth; later activation explicitly chooses StressFreeBirth or PreserveHistory."); }
                ptree mask;mask.put("time_interval.start",number(region.inactive_start));mask.put("time_interval.end",number(start));mask.put("material_ids",id);masks.push_back({"deactivated_subdomain",mask});
                ptree birth;birth.put("material_id",id);birth.put("time",number(start));birth.put("reset",policy=="StressFreeBirth");births.push_back({"birth",birth});region.born=true;region.active=true;
            }
            else
            {
                keys(*action,{"material_ids"});
                int const id=action->get<int>("material_ids");
                auto found=regions.find(id);
                if (found==regions.end() || !found->second.active) { fail("Construction removal must match an active material ID."); }
                found->second.active=false;found->second.inactive_start=start+1e-8*std::max(1.,std::abs(start)/2.);
            }
            p.erase(action_name);
        }
        start+=finiteNumber(p,"duration");
    }
    for (auto const& [material,region]:regions)
    {
        if (region.active) { continue; }
        ptree mask;mask.put("time_interval.start",number(region.inactive_start));mask.put("time_interval.end",number(end));mask.put("material_ids",material);masks.push_back({"deactivated_subdomain",mask});
    }
    for (auto& [tag,v]:root.get_child("process_variables"))
    {
        if (v.get<std::string>("name","")!=process.get<std::string>("process_variables.pressure") && v.get<std::string>("name","")!=process.get<std::string>("process_variables.displacement")) { continue; }
        if (v.count("deactivated_subdomains")) { fail("Construction owns both active-domain masks."); }
        v.put_child("deactivated_subdomains",masks);
    }
    process.put_child("construction_births",births);return true;
}
}
