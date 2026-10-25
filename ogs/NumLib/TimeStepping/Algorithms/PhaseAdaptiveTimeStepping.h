// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include "TimeStepAlgorithm.h"
#include "BaseLib/ConfigTree.h"
#include <algorithm>
#include <limits>
namespace NumLib
{
struct AdaptivePhaseInterval
{
    double start,end,initial,minimum,maximum;
    bool adaptive;
    int grow_below_iterations=4;
    int reduce_above_iterations=8;
};
/// Nonlinear-iteration adaptation; endpoints are mandatory and never crossed.
/// This controls solver robustness, not an estimate of transient truncation error.
class PhaseAdaptiveTimeStepping final : public TimeStepAlgorithm
{
public:
    PhaseAdaptiveTimeStepping(double start,double end,std::vector<AdaptivePhaseInterval> phases)
        :TimeStepAlgorithm(start,end),phases_(std::move(phases))
    {
        if (phases_.empty() || !std::isfinite(start) || !std::isfinite(end) || end<=start) { OGS_FATAL("Invalid adaptive phase interval bounds."); }
        double previous=phases_.front().start;
        for (auto const& p:phases_)
        {
            if (!std::isfinite(p.start) || !std::isfinite(p.end) || !(p.end>p.start) || Time{previous}!=Time{p.start} ||
                !std::isfinite(p.initial) || !std::isfinite(p.minimum) || !std::isfinite(p.maximum) || !(0<p.minimum && p.minimum<=p.initial && p.initial<=p.maximum) || p.grow_below_iterations<0 || p.reduce_above_iterations<=p.grow_below_iterations)
            { OGS_FATAL("Invalid adaptive phase increments or noncontiguous schedule."); }
            previous=p.end;
        }
        if (Time{start}<Time{phases_.front().start} || Time{end}>Time{phases_.back().end}) { OGS_FATAL("Adaptive interval does not cover the requested run."); }
    }
    AdaptivePhaseInterval const& phase(double t) const
    {
        for (auto const& p:phases_) { if (Time{t}<Time{p.end}) { return p; } }
        return phases_.back();
    }
    double next(double,int iterations,TimeStep& previous,TimeStep& current) override
    {
        bool const failed=!current.isAccepted();
        double const t=failed?current.previous()():current.current()();
        auto const& p=phase(t);
        double dt=p.initial;
        if (failed)
        {
            if (!p.adaptive || current.dt()<=p.minimum*(1+1e-12))
            { OGS_FATAL("Adaptive phase reached its minimum increment without convergence."); }
            dt=std::max(p.minimum,current.dt()*.5);
        }
        else if (current.timeStepNumber()!=0 && Time{t}>Time{p.start})
        {
            dt=p.adaptive?std::clamp(current.dt()*(iterations<=p.grow_below_iterations?1.5:iterations>p.reduce_above_iterations?.5:1.),p.minimum,p.maximum):p.initial;
        }
        double const remaining=std::max(0.,p.end-t);
        // A rounded nominal increment must land on the mandatory endpoint,
        // rather than leave a second, sub-minimum interval of a few ulps.
        // This only adjusts the proposed increment, never an accepted state.
        double const roundoff=64*std::numeric_limits<double>::epsilon()*
            std::max({1.,std::abs(t),std::abs(p.end),std::abs(dt)});
        dt=std::abs(dt-remaining)<=roundoff?remaining:std::min(dt,remaining);
        if (failed)
        {
            auto const n=current.timeStepNumber();
            previous=TimeStep{Time{t},Time{t}+dt,n};
            current=previous;current.setAccepted(false);
        }
        return dt;
    }
    bool canReduceTimestepSize(TimeStep const& a,TimeStep const& b) const override
    {
        auto const& attempt=a.current()>b.current()?a:b;
        auto const& p=phase(attempt.previous()());
        return p.adaptive && attempt.dt()>p.minimum*(1+1e-12);
    }
private:
    std::vector<AdaptivePhaseInterval> phases_;
};
inline std::unique_ptr<TimeStepAlgorithm> createPhaseAdaptiveTimeStepper(BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type","PhaseAdaptive");
    double const start=config.getConfigParameter<double>("t_initial");
    double const end=config.getConfigParameter<double>("t_end");
    std::vector<AdaptivePhaseInterval> phases;
    for (auto const& p:config.getConfigSubtreeList("interval"))
    { phases.push_back({p.getConfigParameter<double>("start"),p.getConfigParameter<double>("end"),
        p.getConfigParameter<double>("initial"),p.getConfigParameter<double>("minimum"),p.getConfigParameter<double>("maximum"),p.getConfigParameter<bool>("adaptive"),p.getConfigParameter<int>("grow_below_iterations",4),p.getConfigParameter<int>("reduce_above_iterations",8)}); }
    if (phases.empty() || !std::isfinite(start) || !std::isfinite(end) || end<=start)
    { OGS_FATAL("Invalid adaptive phase schedule."); }
    return std::make_unique<PhaseAdaptiveTimeStepping>(start,end,std::move(phases));
}
}
