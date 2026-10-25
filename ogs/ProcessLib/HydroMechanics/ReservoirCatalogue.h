// SPDX-License-Identifier: BSD-3-Clause
#pragma once

namespace ProcessLib::HydroMechanics
{
/// User-facing reservoir level, shared by exposed ground and far-field head.
inline void prepareReservoirCatalogue(boost::property_tree::ptree& process)
{
    using namespace InitialGroundwaterDetail;
    auto reservoirs=process.get_child_optional("reservoirs");
    if (!reservoirs) { return; }
    if (reservoirs->empty() || process.count("hydraulic_controls"))
    { fail("reservoirs requires entries and replaces hydraulic_controls"); }
    ptree controls;
    std::map<std::string,std::string> companions;
    std::set<std::string> ids;
    for (auto const& [tag,reservoir]:*reservoirs)
    {
        if (tag!="reservoir") { fail("reservoirs requires reservoir entries"); }
        keys(reservoir,{"<xmlattr>","water_level","surface","far_field_mesh","apply_water_traction"});
        keys(reservoir.get_child("<xmlattr>"),{"id"});
        keys(reservoir.get_child("water_level"),{});
        auto const id=reservoir.get<std::string>("<xmlattr>.id");
        if (id.empty() || id.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_")!=std::string::npos || !ids.insert(id).second)
        { fail("reservoir requires unique simple ID"); }
        double const level=finiteNumber(reservoir,"water_level");
        auto const traction = reservoir.get<std::string>("apply_water_traction", "true");
        if (traction != "true" && traction != "false")
        { fail("reservoir apply_water_traction must be true or false"); }
        if (reservoir.count("apply_water_traction"))
        { keys(reservoir.get_child("apply_water_traction"), {}); }
        ptree boundary;
        boundary.put("<xmlattr>.id",id);
        boundary.put("type","ReservoirSeepage");
        boundary.put("initial_head",number(level));
        boundary.put("pressure_policy","Signed");
        boundary.put("apply_surface_water",traction);
        boundary.put("water_at_integration_points",true);
        auto surface=reservoir.get_child("surface");
        if (surface.count("mesh"))
        {
            keys(surface,{"mesh"});keys(surface.get_child("mesh"),{});
            surface.put("type","Mesh");
        }
        boundary.put_child("exterior_surface",surface);
        controls.push_back({"boundary",boundary});
        if (auto far=reservoir.get_optional<std::string>("far_field_mesh"))
        {
            keys(reservoir.get_child("far_field_mesh"),{});
            if (far->empty()) { fail("far_field_mesh must name a boundary mesh"); }
            auto const companion=id+"FarField";
            if (!ids.insert(companion).second) { fail("reservoir far-field ID collision"); }
            companions.emplace(id,companion);
            boundary.erase("exterior_surface");
            boundary.put("<xmlattr>.id",companion);
            boundary.put("type","HydrostaticHead");
            boundary.put("mesh",*far);
            boundary.put("apply_surface_water","false");
            controls.push_back({"boundary",boundary});
        }
    }
    for (auto& [tag,phase]:process.get_child("phases"))
    {
        auto heads=phase.get_child_optional("heads");
        if (!heads) { continue; }
        ptree expanded;
        for (auto const& [head_tag,head]:*heads)
        {
            auto const id=head.get<std::string>("<xmlattr>.ref");
            // Companion actions are internal; one reservoir action changes both.
            if (!ids.contains(id) || std::any_of(companions.begin(),companions.end(),[&](auto const& item){ return item.second==id; }))
            { fail("reservoir phase head must reference a declared reservoir ID"); }
            expanded.push_back({head_tag,head});
            if (companions.contains(id))
            { auto copy=head;copy.put("<xmlattr>.ref",companions.at(id));expanded.push_back({head_tag,copy}); }
        }
        phase.put_child("heads",expanded);
    }
    process.put_child("hydraulic_controls",controls);
    process.erase("reservoirs");
}
} // namespace ProcessLib::HydroMechanics
