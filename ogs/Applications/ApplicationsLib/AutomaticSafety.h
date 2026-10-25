// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <boost/property_tree/xml_parser.hpp>
#include <filesystem>
#include <fstream>
#include <functional>
#include <cmath>
#include <optional>
#include "ProcessLib/HydroMechanics/NativeCheckpoint.h"
namespace ApplicationsLib
{
struct AutomaticSafetyPlan
{
    boost::property_tree::ptree project;
    std::filesystem::path directory,checkpoint;
    double tolerance,maximum;
    std::string phase_id;
};
inline void writeAutomaticProject(std::filesystem::path const& file,boost::property_tree::ptree const& project)
{
    std::ofstream out(file);
    if (!out) { throw std::runtime_error("Cannot create automatic safety project."); }
    boost::property_tree::write_xml(out,project);
}
inline std::unique_ptr<AutomaticSafetyPlan> prepareAutomaticSafety(std::stringstream& input,
    std::filesystem::path const& project_file,std::string& output)
{
    if (input.str().find("<search")==std::string::npos) { return {}; }
    using namespace ProcessLib::HydroMechanics::InitialGroundwaterDetail;
    ptree document;std::stringstream parsed(input.str());boost::property_tree::read_xml(parsed,document,boost::property_tree::xml_parser::trim_whitespace|boost::property_tree::xml_parser::no_comments);
    auto& root=document.get_child("OpenGeoSysProject");
    if (root.get_child("processes").size()!=1 || root.get<std::string>("processes.process.type","")!="HYDRO_MECHANICS")
    { fail("Automatic safety requires one native HM process."); }
    auto& p=root.get_child("processes.process");
    if (!p.count("phases") || p.count("checkpoints") || p.count("restart"))
    { fail("Automatic safety owns its working checkpoint; remove manual checkpoint/restart controls."); }
    auto& phases=p.get_child("phases");
    if (phases.size()<2) { fail("Automatic safety requires prior working/initialisation phases."); }
    auto& last=phases.back().second;
    if (last.get<std::string>("type","")!="Safety" || last.get<std::string>("search","")!="Automatic")
    { fail("Automatic search must be the final Safety phase."); }
    for (auto i=phases.begin();std::next(i)!=phases.end();++i)
    { if (i->second.get<std::string>("type","")=="Safety") { fail("Automatic search currently starts after working phases, not after another Safety phase."); } }
    keys(last,{"<xmlattr>","type","search","factor_tolerance","maximum_factor","pore_pressure_policy","reset_displacement_reference"});
    double const tolerance=finiteNumber(last,"factor_tolerance"),maximum=finiteNumber(last,"maximum_factor");
    if (!(tolerance>0. && maximum>1. && tolerance<maximum-1. && maximum<=100.))
    { fail("Automatic safety requires 0<tolerance<maximum_factor-1 and 1<maximum_factor<=100."); }
#if defined(USE_PETSC) || defined(USE_LIS)
    fail("Automatic safety currently requires the serial Eigen build.");
#endif
    for (auto& [tag,mesh]:root.get_child("meshes"))
    { if (tag=="mesh") { mesh.data()=std::filesystem::weakly_canonical(project_file.parent_path()/mesh.data()).string(); } }
    last.erase("search");last.erase("factor_tolerance");last.erase("maximum_factor");last.put("steps",1);last.put("target_factor",1.);last.put("verify_equilibrium",true);
    // A default force/yield gate is explicit in every retained generated input.
    if (phases.front().second.count("method"))
    {
        auto& initial=phases.front().second;
        if (!initial.count("equilibrium_check"))
        { initial.put("equilibrium_check.force_tolerance",.001);initial.put("equilibrium_check.yield_tolerance",.1); }
    }
    else if (!p.get_child_optional("initialisation.equilibrium_check"))
    { p.put("initialisation.equilibrium_check.force_tolerance",.001);p.put("initialisation.equilibrium_check.yield_tolerance",.1); }
    auto previous=phases.end();--previous;--previous;
    ptree point;point.put("<xmlattr>.id","AutomaticWorking");point.put("after_phase",previous->second.get<std::string>("<xmlattr>.id"));
    ptree checkpoints;checkpoints.push_back({"checkpoint",point});p.put_child("checkpoints",checkpoints);
    auto plan=std::make_unique<AutomaticSafetyPlan>();plan->project=document;plan->tolerance=tolerance;plan->maximum=maximum;plan->phase_id=last.get<std::string>("<xmlattr>.id");
    plan->directory=std::filesystem::absolute(output.empty()?".":output)/"automatic_safety";
    if (!std::filesystem::create_directory(plan->directory)) { fail("automatic_safety output already exists; use a fresh output directory."); }
    std::filesystem::create_directory(plan->directory/"working");
    output=(plan->directory/"working").string();plan->checkpoint=plan->directory/"working/checkpoint_AutomaticWorking.hmcp";
    document.put("OpenGeoSysProject.processes.process.checkpoints.checkpoint.stop_after",true);
    writeAutomaticProject(plan->directory/"working/input.prj",document);
    input.str("");input.clear();boost::property_tree::write_xml(input,document);input.seekg(0);
    return plan;
}
struct AutomaticSafetyTrial { double target;int steps;bool accepted;std::string directory; };
inline bool executeAutomaticSafety(AutomaticSafetyPlan const& plan,
    std::function<bool(std::string const&,std::string const&)> const& run)
{
    using namespace ProcessLib::HydroMechanics::NativeCheckpointDetail;
    auto const original=hash(fileBytes(plan.checkpoint));std::vector<AutomaticSafetyTrial> trials;
    double lower=1.,upper=plan.maximum;std::string status="searching";
    auto report=[&]()
    {
        std::ofstream out(plan.directory/"result.json");out<<std::setprecision(17);
        out<<"{\n\"status\":\""<<status<<"\",\n\"accepted_lower_factor\":"<<lower<<",\n\"candidate_upper_factor\":";
        if (status=="ceiling_reached") { out<<"null"; } else { out<<upper; }
        out<<",\n\"factor_tolerance\":"<<plan.tolerance<<",\n\"interval_width\":"<<(upper-lower)<<",\n\"interpretation\":\"Numerical accepted-equilibrium/failure interval; failure alone is not a mathematical collapse proof\",\n\"trials\":[";
        bool first=true;for (auto const& t:trials)
        { if (!first) { out<<","; }first=false;out<<"{\"target\":"<<t.target<<",\"steps\":"<<t.steps<<",\"accepted\":"<<(t.accepted?"true":"false")<<",\"directory\":\""<<t.directory<<"\"}"; }
        out<<"]\n}\n";
    };
    auto trial=[&](double target,double step)
    {
        if (trials.size()>=80) { status="trial_limit_unresolved";report();throw std::runtime_error("Automatic safety trial limit reached; unresolved result retained."); }
        auto document=plan.project;auto& p=document.get_child("OpenGeoSysProject.processes.process");
        p.put("restart.file",plan.checkpoint.string());auto& phase=p.get_child("phases").back().second;
        int const steps=std::max(2,static_cast<int>(std::ceil((target-1.)/step)));
        phase.put("target_factor",target);phase.put("steps",steps);
        auto const folder=plan.directory/("trial_"+std::to_string(trials.size()));std::filesystem::create_directory(folder);
        writeAutomaticProject(folder/"input.prj",document);
        INFO("Automatic safety trial {}: RF={}, continuation steps={}.",trials.size(),target,steps);
        bool accepted=run((folder/"input.prj").string(),folder.string());
        if (hash(fileBytes(plan.checkpoint))!=original) { throw std::runtime_error("Automatic safety input checkpoint was modified."); }
        trials.push_back({target,steps,accepted,folder.filename().string()});report();return accepted;
    };
    if (!trial(1.,.01)) { status="reference_state_unresolved";report();return false; }
    // A second, finer continuation pass independently checks RF-increment
    // sensitivity. It reuses the same original accepted working archive.
    for (int pass=0;pass<2;++pass)
    {
        double const step=pass==0?.01:.005;
        double start=pass==0?1.:std::max(1.,lower-.05);
        double high=pass==0?std::min(plan.maximum,1.25):std::min(plan.maximum,upper+.05);
        lower=1.;upper=plan.maximum;
        if (start>1. && trial(start,step)) { lower=start; }
        while (trial(high,step))
        {
            lower=high;
            if (high==plan.maximum) { status="ceiling_reached";report();INFO("Automatic safety ceiling reached: RF >= {} (no limiting FoS found).",lower);return true; }
            high=std::min(plan.maximum,high*1.25);
        }
        upper=high;
        while (upper-lower>plan.tolerance)
        { double const mid=lower+(upper-lower)/2.;if (trial(mid,step)) { lower=mid; } else { upper=mid; } }
    }
    status="numerical_interval_resolved";report();
    INFO("Automatic safety result: accepted RF={}, failure candidate RF={}, width={} (requested {}).",lower,upper,upper-lower,plan.tolerance);
    return true;
}
}
