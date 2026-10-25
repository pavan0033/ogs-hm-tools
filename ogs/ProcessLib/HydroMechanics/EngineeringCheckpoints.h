// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include "NativeCheckpoint.h"

namespace ProcessLib::HydroMechanics
{
struct EngineeringCheckpointSchedule
{
    double restart_time=0.;
    std::optional<double> stop_time;
};
inline std::string checkpointModelIdentity(
    boost::property_tree::ptree root, std::size_t const phase_count,
    std::filesystem::path const& project_directory)
{
    using namespace NativeCheckpointDetail;
    auto& process=root.get_child("processes.process");
    process.erase("checkpoints");process.erase("restart");process.erase("native_checkpoint");
    auto& phases=process.get_child("phases");
    std::size_t count=0;
    for (auto i=phases.begin();i!=phases.end();)
    { if (++count>phase_count) { i=phases.erase(i); } else { ++i; } }
    root.get_child("time_loop").erase("output");root.erase("test_definition");
    // Solver algorithms do not define the accepted physical state. Retain all
    // material, mesh, loading, boundary and accepted-phase-prefix checks while
    // permitting a controlled solver comparison from one identical capsule.
    root.erase("linear_solvers");root.erase("nonlinear_solvers");
    std::ostringstream canonical;boost::property_tree::write_xml(canonical,root);
    std::string identity=canonical.str();
    for (auto const& [tag,entry]:root.get_child("meshes"))
    { if (tag=="mesh") { identity+=entry.data()+":"+hex(hash(fileBytes(project_directory/entry.data())))+"\n"; } }
    return hex(hash(identity));
}
inline EngineeringCheckpointSchedule prepareEngineeringCheckpoints(
    boost::property_tree::ptree const& root,
    boost::property_tree::ptree& process,
    std::filesystem::path const& project_directory,
    std::filesystem::path const& output_directory)
{
    using namespace InitialGroundwaterDetail;
    EngineeringCheckpointSchedule result;
    if (!process.count("checkpoints") && !process.count("restart")) { return result; }
    if (process.get<std::string>("coupling_scheme","monolithic")!="monolithic" ||
        (process.get<std::string>("initialisation.method","")!="Gravity" &&
         process.get<std::string>("initialisation.method","")!="AssignedStress") ||
        root.get_child("processes").size()!=1 ||
        root.get<std::string>("time_loop.processes.process.time_discretization.type","")!="BackwardEuler")
    { fail("checkpoints require one monolithic Gravity or AssignedStress HM process, native phases and BackwardEuler"); }
    std::ostringstream raw;boost::property_tree::write_xml(raw,root);
    if (raw.str().find("deactivated_subdomains")!=std::string::npos ||
        raw.str().find("fixed_output_times")!=std::string::npos ||
        raw.str().find("RelativeDirichlet")!=std::string::npos ||
        raw.str().find("compensate_non_equilibrium_initial_residuum")!=std::string::npos ||
        root.get_child_optional("time_loop.output.meshes") ||
        process.get<bool>("linear",false) || process.count("initial_stress"))
    { fail("checkpoint first scope excludes activation, assigned stress, manual relative BCs, residual compensation, linear caching and submesh output"); }
    std::map<std::string,std::pair<std::size_t,double>> ends;
    double time=0.;std::size_t count=0;double initialisation_time=0.;
    for (auto const& [tag,phase]:process.get_child("phases"))
    {
        if (tag!="phase") { fail("checkpoint requires native phase entries"); }
        time+=phaseAnalysisDuration(phase);
        if (++count==1) { initialisation_time=time; }
        if (!ends.emplace(phase.get<std::string>("<xmlattr>.id"),std::make_pair(count,time)).second)
        { fail("checkpoint phase IDs must be unique"); }
    }
    ptree settings;settings.put("initialisation_time",number(initialisation_time));
    if (auto entries=process.get_child_optional("checkpoints"))
    {
        std::set<std::string> identifiers;std::set<double> times;
        if (entries->empty()) { fail("checkpoints must contain checkpoint entries"); }
        for (auto const& [tag,entry]:*entries)
        {
            if (tag!="checkpoint") { fail("checkpoints requires checkpoint entries"); }
            keys(entry,{"<xmlattr>","after_phase","stop_after"});keys(entry.get_child("<xmlattr>"),{"id"});
            keys(entry.get_child("after_phase"),{});
            auto const id=entry.get<std::string>("<xmlattr>.id");auto const phase=entry.get<std::string>("after_phase");
            if (id.empty() || id.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_")!=std::string::npos ||
                !identifiers.insert(id).second || !ends.contains(phase) || !times.insert(ends.at(phase).second).second)
            { fail("checkpoint requires unique simple id and unique known after_phase"); }
            bool stop=false;
            if (entry.count("stop_after"))
            { keys(entry.get_child("stop_after"),{});auto const value=entry.get<std::string>("stop_after");
              if (value!="true" && value!="false") { fail("checkpoint stop_after requires true/false"); }
              stop=value=="true"; }
            if (stop)
            { if (result.stop_time) { fail("only one stop_after checkpoint is supported"); }
              result.stop_time=ends.at(phase).second; }
            ptree point;point.put("id",id);point.put("time",number(ends.at(phase).second));
            point.put("file",std::filesystem::absolute(output_directory/("checkpoint_"+id+".hmcp")).string());
            point.put("model_identity",checkpointModelIdentity(root,ends.at(phase).first,project_directory));
            settings.push_back({"point",point});
        }
    }
    if (auto restart=process.get_child_optional("restart"))
    {
        keys(*restart,{"file"});keys(restart->get_child("file"),{});
        auto const path=std::filesystem::absolute(project_directory/restart->get<std::string>("file"));
        auto const state=decodeNativeCheckpoint(NativeCheckpointDetail::fileBytes(path));
        bool found=false;
        for (auto const& [id,end]:ends)
        {
            if (std::abs(state.time-end.second)<=1e-12*std::max(1.,std::abs(end.second)))
            {
                auto const identity=checkpointModelIdentity(root,end.first,project_directory);
                if (state.initialisation_time!=initialisation_time || state.model_identity!=identity)
                { fail("checkpoint model/mesh/accepted-phase-prefix identity mismatch"); }
                settings.put("restart_identity",identity);found=true;break;
            }
        }
        if (!found || state.time>=time-1e-12*std::max(1.,std::abs(time)) || (result.stop_time && *result.stop_time<=state.time+1e-12*std::max(1.,std::abs(state.time))))
        { fail("restart must match a native phase end with remaining working steps"); }
        result.restart_time=state.time;settings.put("restart_file",path.string());
    }
    process.erase("checkpoints");process.erase("restart");process.put_child("native_checkpoint",settings);
    return result;
}
inline void cropCheckpointTimeStepping(boost::property_tree::ptree& root,
                                      EngineeringCheckpointSchedule const& schedule)
{
    using namespace InitialGroundwaterDetail;
    if (schedule.restart_time==0. && !schedule.stop_time) { return; }
    auto& stepper=root.get_child("time_loop.processes.process.time_stepping");
    double const end=schedule.stop_time.value_or(stepper.get<double>("t_end"));
    if (stepper.get<std::string>("type")=="PhaseAdaptive")
    { stepper.put("t_initial",number(schedule.restart_time));stepper.put("t_end",number(end));return; }
    ptree pairs;double time=0.;std::string previous_increments;
    for (auto const& [tag,pair]:stepper.get_child("timesteps"))
    {
        double const next=time+pair.get<int>("repeat")*pair.get<double>("delta_t");
        if (time>=schedule.restart_time-1e-10 && next<=end+1e-10) { pairs.push_back({tag,pair}); }
        if (next<=schedule.restart_time+1e-10)
        {
            for (int i=0;i<pair.get<int>("repeat");++i)
            { previous_increments+=number(pair.get<double>("delta_t"))+" "; }
        }
        time=next;
    }
    if (pairs.empty()) { fail("checkpoint crop produced no working steps"); }
    if (schedule.restart_time>0.)
    { stepper.put("initial_time_origin",0.);stepper.put("initial_time_increments",previous_increments); }
    stepper.put("t_initial",number(schedule.restart_time));stepper.put("t_end",number(end));stepper.put_child("timesteps",pairs);
}
} // namespace ProcessLib::HydroMechanics
