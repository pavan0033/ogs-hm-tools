// SPDX-License-Identifier: BSD-3-Clause
#pragma once
namespace ProcessLib::HydroMechanics
{
/// Compile a native removal action into identical pressure/displacement masks.
/// Exposed-face pressure and water-traction parameters remain explicit inputs.
inline std::vector<boost::property_tree::ptree> prepareEngineeringExcavation(boost::property_tree::ptree& root,
                                         boost::property_tree::ptree& process)
{
    using namespace InitialGroundwaterDetail;
    if (!process.count("phases")) { return {}; }
    bool requested=false;
    for (auto const& [tag,phase]:process.get_child("phases"))
    { requested |= phase.count("excavate")!=0; }
    if (!requested) { return {}; }
    if (process.count("mechanical_drainage") ||
        process.get<std::string>("coupling_scheme","monolithic")!="monolithic")
    { fail("Native excavation currently requires ordinary monolithic coupled HM without strict mechanical drainage."); }
    for (auto const& [tag,variable]:root.get_child("process_variables"))
    {
        if(variable.count("deactivated_subdomains"))
        { fail("Native excavation owns the deactivated_subdomains masks; remove competing stock masks."); }
    }
    double end=0.;for (auto const& [tag,phase]:process.get_child("phases")) { end+=finiteNumber(phase,"duration"); }
    double start=0.;std::set<int> removed;
    for (auto& [tag,phase]:process.get_child("phases"))
    {
        if (auto action=phase.get_child_optional("excavate"))
        {
            if (start<=0. || phase.get<std::string>("type","")!="Loading" ||
                (phase.get<std::string>("hydraulic_policy","")!="Coupled" &&
                 phase.get<std::string>("hydraulic_policy","")!="SteadyStateSeepage"))
            { fail("excavate requires a later Coupled or SteadyStateSeepage Loading phase."); }
            keys(*action,{"material_ids","pressure_parameter","stress_release","hydraulic_boundary"});
            auto const hydraulic=action->get<std::string>("hydraulic_boundary","PrescribedPressure");
            if(hydraulic!="PrescribedPressure" && hydraulic!="Natural")
            { fail("excavation hydraulic_boundary must be PrescribedPressure or Natural."); }
            auto const release=action->get<std::string>("stress_release","Instantaneous");
            if(release!="Instantaneous" && release!="Ramp")
            { fail("excavation stress_release must be Instantaneous or Ramp."); }
            if(release=="Ramp")
            {
                if(phase.count("load") || phase.count("heads"))
                { fail("Ramped excavation first scope separates mechanical/head changes into other phases."); }
                if(root.get_child("time_loop.processes").size()!=1)
                { fail("Ramped excavation first scope requires a single HM equation system."); }
                double const duration=finiteNumber(phase,"duration");
                double const minimum=phase.count("stepping") ?
                    duration*finiteNumber(phase.get_child("stepping"),"minimum_increment") :
                    duration/finiteNumber(phase,"steps");
                if(minimum<=5e-8*std::max(1.,std::abs(start)/2.))
                { fail("Ramped excavation minimum time increment must resolve the phase-start removal event."); }
                for(auto const& [tag,relation]:process)
                {
                    if(tag!="constitutive_relation") { continue; }
                    auto const behaviour=relation.get<std::string>("behaviour","");
                    if(relation.get<std::string>("type")!="MFront" ||
                       (behaviour!="MohrCoulombAbboSloan" && behaviour!="MohrCoulombAbboSloanTension" && behaviour!="PublishedHS"))
                    { fail("Ramped excavation requires transactional MFront MC soil or PublishedHS."); }
                    if(behaviour=="PublishedHS" && !relation.get_child_optional("substepping"))
                    { fail("PublishedHS excavation requires explicit transactional substepping."); }
                    for(auto const& [property_tag,property]:relation.get_child("material_properties"))
                    {
                        auto const name=property.get<std::string>("<xmlattr>.parameter");
                        for(auto const& [parameter_tag,parameter]:root.get_child("parameters"))
                        {
                            if(parameter.get<std::string>("name","")!=name) { continue; }
                            auto const type=parameter.get<std::string>("type");
                            if(type!="Constant" && type!="Function")
                            { fail("Ramped excavation material parameters must be Constant or static Function."); }
                            for(auto const& [key,value]:parameter)
                            { if(key=="expression" && std::regex_search(value.data(),std::regex(R"(\bt\b)")))
                              { fail("Ramped excavation material parameters must be time independent."); } }
                        }
                    }
                }
                for(auto const& [variable_tag,variable]:root.get_child("process_variables"))
                {
                    if(!variable.count("source_terms")) { continue; }
                    for(auto const& [source_tag,source]:variable.get_child("source_terms"))
                    {
                        if(auto ramp=source.get_child_optional("distributed_load_ramp"))
                        {
                            double const begin=finiteNumber(*ramp,"start_time");
                            double const finish=finiteNumber(*ramp,"end_time");
                            if(begin<start+duration-1e-10 && finish>start+1e-10)
                            { fail("Ramped excavation separates plate distributed-load ramps into other phases."); }
                        }
                        double const activation=source.get<double>("activation_time",0.);
                        double const lockoff=source.get<double>("prestress_end_time",-1.);
                        if((activation>=start-1e-10 && activation<start+duration-1e-10) ||
                           (lockoff>start+1e-10 && activation<start+duration-1e-10))
                        { fail("Ramped excavation separates structural activation/jacking into other phases."); }
                    }
                }
                auto const solver=root.get_child("time_loop.processes").front().second.get<std::string>("nonlinear_solver");
                bool newton=false;
                for(auto const& [tag,setting]:root.get_child("nonlinear_solvers"))
                { newton |= setting.get<std::string>("name")==solver && setting.get<std::string>("type")=="Newton"; }
                if(!newton) { fail("Ramped excavation requires Newton assembly."); }
                ptree interval;interval.put("start",number(start));interval.put("end",number(start+finiteNumber(phase,"duration")));
                if(!process.count("excavation_release")) { process.put_child("excavation_release",ptree{}); }
                process.get_child("excavation_release").push_back({"interval",interval});
            }
            std::istringstream ids(action->get<std::string>("material_ids"));int id;bool any=false;
            while (ids>>id) { if (id<0 || !removed.insert(id).second) { fail("excavation material IDs must be unique nonnegative removals."); } any=true; }
            if (!any || !ids.eof()) { fail("excavation material_ids must be nonempty integers."); }
            auto const parameter=action->get<std::string>("pressure_parameter");bool found=false;
            for (auto const& [tag,p]:root.get_child("parameters")) { found |= p.get<std::string>("name","")==parameter; }
            if (!found) { fail("excavation pressure_parameter must name an existing exposed-face pressure parameter."); }
            if(release=="Ramp")
            {
                for(auto const& [tag,p]:root.get_child("parameters"))
                {
                    if(p.get<std::string>("name","")!=parameter) { continue; }
                    auto const type=p.get<std::string>("type");
                    if(type!="Constant" && type!="Function")
                    { fail("Ramped excavation pressure_parameter must be Constant or Function."); }
                }
            }
            // An accepted phase-end belongs to the old active domain. The
            // removal is first assembled in the following step, as in the
            // independently verified stock excavation control.
            double const event=start+1e-8*std::max(1.,std::abs(start)/2.);
            for (auto& [vtag,variable]:root.get_child("process_variables"))
            {
                bool const pressure=variable.get<std::string>("name","")==process.get<std::string>("process_variables.pressure");
                bool const displacement=variable.get<std::string>("name","")==process.get<std::string>("process_variables.displacement");
                if (!pressure && !displacement) { continue; }
                ptree sub;sub.put("time_interval.start",number(event));sub.put("time_interval.end",number(end));sub.put("material_ids",action->get<std::string>("material_ids"));
                if (pressure && hydraulic=="PrescribedPressure") { sub.put("boundary_parameter",parameter); }
                if (!variable.count("deactivated_subdomains"))
                { variable.put_child("deactivated_subdomains",ptree{}); }
                variable.get_child("deactivated_subdomains").push_back({"deactivated_subdomain",sub});
            }
            phase.erase("excavate");
        }
        start+=finiteNumber(phase,"duration");
    }
    std::vector<ptree> water_tractions;
    for (auto& [tag,variable]:root.get_child("process_variables"))
    {
        if (variable.get<std::string>("name","")!=process.get<std::string>("process_variables.displacement")) { continue; }
        auto& boundaries=variable.get_child("boundary_conditions");
        for (auto it=boundaries.begin();it!=boundaries.end();)
        {
            if (it->second.get<std::string>("type","")=="NormalTraction")
            {
                keys(it->second,{"mesh","type","parameter"});
                water_tractions.push_back(it->second);it=boundaries.erase(it);
            }
            else { ++it; }
        }
    }
    return water_tractions;
}
}
