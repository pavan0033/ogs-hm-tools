// SPDX-License-Identifier: BSD-3-Clause
#pragma once
namespace ProcessLib::HydroMechanics
{
/// A Safety interval is dimensionless numerical continuation, never one second.
inline double phaseAnalysisDuration(boost::property_tree::ptree const& phase)
{
    return phase.get<std::string>("type")=="Safety" ? 1. :
        InitialGroundwaterDetail::finiteNumber(phase,"duration");
}
inline void prepareEngineeringSafety(boost::property_tree::ptree& root,
                                     boost::property_tree::ptree& process)
{
    using namespace InitialGroundwaterDetail;
    if (!process.count("phases")) { return; }
    bool requested=false;
    for (auto const& [tag,phase]:process.get_child("phases"))
    { requested |= tag=="phase" && phase.get<std::string>("type","")=="Safety"; }
    if (!requested) { return; }
    if (process.get<std::string>("coupling_scheme","monolithic")!="monolithic" ||
        process.get<std::string>("initialisation.method","")!="Gravity" ||
        process.get<bool>("linear",false) || root.get_child("processes").size()!=1 ||
        root.get<std::string>("time_loop.processes.process.time_discretization.type","")!="BackwardEuler")
    { fail("Safety requires one nonlinear monolithic Gravity/BackwardEuler native HM process."); }
    for (auto const& [tag,parameter]:root.get_child("parameters"))
    {
        if (tag!="parameter" || parameter.get<std::string>("type","")=="Constant") { continue; }
        auto const name=parameter.get<std::string>("name");
        bool used=false;
        for (auto const& [relation_tag,relation]:process)
        {
            if (relation_tag!="constitutive_relation") { continue; }
            for (auto const& [property_tag,property]:relation.get_child("material_properties"))
            {
                if (property.get<std::string>("<xmlattr>.parameter")!=name) { continue; }
                auto const property_name=property.get<std::string>("<xmlattr>.name");
                if (property_name!="YoungModulus" && property_name!="Cohesion")
                { fail("Safety spatial Function scope permits YoungModulus and Cohesion only."); }
                used=true;
            }
        }
        // A static spatial load shape is independent of the Safety coordinate.
        // Its phase amplitude remains frozen by the engineering load schedule.
        if (auto loads=process.get_child_optional("mechanical_loads"))
        {
            for (auto const& [load_tag,load]:*loads)
            {
                if (load_tag=="load" &&
                    load.get<std::string>("type","")=="NormalPressure" &&
                    load.get<std::string>("spatial_scale_parameter","")==name)
                { used=true; }
            }
        }
        if (!used || parameter.get<std::string>("type","")!="Function" || parameter.count("expression")!=1 ||
            std::regex_search(parameter.get<std::string>("expression",""),std::regex("(^|[^a-zA-Z0-9_])t([^a-zA-Z0-9_]|$)")))
        { fail("Safety requires Constant parameters or time-independent scalar Functions for YoungModulus/Cohesion or declared spatial load shapes; RF is not physical time."); }
    }
    auto constant_properties=[&](auto const& self,ptree const& node, bool liquid)->void
    {
        liquid = liquid || node.get<std::string>("type", "") == "AqueousLiquid";
        for (auto const& [tag,child]:node)
        {
            bool const reference_eos = liquid && child.get<std::string>("name", "") == "density" &&
                child.get<std::string>("type", "") == "Linear" &&
                child.get<std::string>("independent_variable.variable_name", "") == "liquid_phase_pressure" &&
                child.get<double>("independent_variable.reference_condition", 1.) == 0. &&
                child.get<double>("independent_variable.slope", -1.) >= 0. &&
                child.get<double>("reference_value", 0.) > 0. &&
                child.count("independent_variable") == 1;
            if (tag=="property" && child.get<std::string>("type","")!="Constant" && !reference_eos)
            { fail("Safety requires Constant MPL properties or a pressure-Linear liquid density EOS referenced at zero gauge pressure."); }
            self(self,child,liquid);
        }
    };
    constant_properties(constant_properties,root.get_child("media"),false);
    for (auto const& [tag,relation]:process)
    { if (tag=="constitutive_relation" && (relation.get<std::string>("type","")!="MFront" ||
        relation.get<std::string>("behaviour","")!="MohrCoulombAbboSloan"))
      { fail("Safety currently supports MohrCoulombAbboSloan only."); } }
    bool tsa=false,strict_esa=false;
    if (auto materials=process.get_child_optional("mechanical_drainage"))
    {
        for (auto const& [tag,material]:*materials)
        {
            if (tag!="material") { continue; }
            bool const total=material.get<std::string>("type","")=="UndrainedTotalStress";
            tsa|=total;strict_esa|=!total;
        }
    }
    bool const flow_history = process.count("hydraulic_controls") != 0;
    if ((tsa && flow_history) || (!tsa && !strict_esa && process.count("hydraulic_controls") != 1))
    { fail("Safety needs hydraulic_controls for ordinary HM; TSA cannot convert a preceding flow stress basis."); }
    ptree settings;double coordinate=0.,previous=1.;bool safety_started=false;std::size_t index=0;
    for (auto& [tag,phase]:process.get_child("phases"))
    {
        if (tag!="phase") { fail("Safety requires native phase entries."); }
        auto const type=phase.get<std::string>("type");
        if (type!="Safety")
        {
            if (safety_started) { fail("Working phases cannot follow Safety continuation."); }
            coordinate+=finiteNumber(phase,"duration");
            if (index==0) { settings.put("initialisation_end",number(coordinate)); }
        }
        else
        {
            keys(phase,{"<xmlattr>","type","steps","stepping","target_factor","pore_pressure_policy","reset_displacement_reference","verify_equilibrium"});
            if (index==0) { fail("Safety requires an accepted initialisation/working state."); }
            auto const policy=phase.get<std::string>("pore_pressure_policy");
            bool const actual=policy=="FrozenAcceptedActualPressure" || policy=="AdditionalUndrained";
            bool const material_policy=policy=="MaterialDrainage";
            if ((tsa && strict_esa) ? !material_policy :
                (tsa ? policy!="NotApplicable" :
                !(actual || (strict_esa && material_policy) || (!strict_esa && policy=="FrozenAcceptedField"))))
            { fail("Safety requires NotApplicable for TSA, an explicit actual-pressure policy for strict ESA, or FrozenAcceptedField for ordinary HM."); }
            if (safety_started && policy!=settings.get<std::string>("pressure_policy"))
            { fail("All Safety phases must retain the same pressure policy in this scope."); }
            double const target=finiteNumber(phase,"target_factor");
            if (phase.get<int>("steps",1)<1 || target<previous || target<1.)
            { fail("Safety target_factor must be finite, nondecreasing and >=1; steps positive."); }
            if (!safety_started)
            {
                settings.put("physical_time",number(coordinate));settings.put("pressure_policy",policy);
                settings.put("working_hydraulic_rates",!tsa && (!strict_esa || flow_history));
                settings.put("verify_equilibrium",phase.get<bool>("verify_equilibrium",false));
                ptree point;point.put("coordinate",number(coordinate));point.put("factor",1.);settings.push_back({"point",point});
            }
            safety_started=true;coordinate+=1.;
            ptree point;point.put("coordinate",number(coordinate));point.put("factor",number(target));settings.push_back({"point",point});
            previous=target;phase.erase("target_factor");phase.erase("pore_pressure_policy");phase.erase("verify_equilibrium");
            phase.put("type",(tsa || (strict_esa && !flow_history))?"MechanicalLoading":"Loading");phase.put("duration",1.);
            if (!tsa && (!strict_esa || flow_history)) { phase.put("hydraulic_policy","FrozenAcceptedField"); }
        }
        ++index;
    }
    auto& output=root.get_child("time_loop.output");
    auto suffix=output.get<std::string>("suffix","_ts_{:timestep}_t_{:time}");
    if (suffix.find(":converged")==std::string::npos) { suffix+="{:converged}"; }
    output.put("suffix",suffix);
    process.put_child("native_safety",settings);
}
} // namespace ProcessLib::HydroMechanics
