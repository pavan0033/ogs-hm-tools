// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <algorithm>
#include <cmath>
#include <optional>
#include <utility>
#include "NewtonStepStrategy.h"

namespace NumLib
{
/// Evaluate every proposed step; never accept an unchecked minimum step.
template <typename Evaluate>
std::optional<std::pair<double,double>> backtrackArmijo(
    double const initial, double const contraction, double const minimum,
    double const armijo, unsigned const max_trials, double const base_merit,
    double const merit_floor, Evaluate&& evaluate)
{
    double step=initial;
    for (unsigned trial=0;trial<max_trials;++trial)
    {
        auto const merit=evaluate(step);
        if (merit && std::isfinite(*merit) &&
            (*merit<=(1.-armijo*step)*base_merit || *merit<=merit_floor))
        { return std::pair{step,*merit}; }
        if (step<=minimum) { break; }
        step=std::max(minimum,step*contraction);
    }
    return std::nullopt;
}

/// Serial monolithic backtracking with freshly assembled, constrained residuals.
class BacktrackingNewtonStrategy final : public NewtonStepStrategy
{
public:
    BacktrackingNewtonStrategy(double initial,double contraction,double minimum,
        double armijo,unsigned max_trials,double merit_tolerance, bool residual_screening=false);
    StepResult applyStep(GlobalVector const& x,GlobalVector const& direction,
        GlobalVector const& residual,GlobalMatrix const& jacobian,
        GlobalVector& next,NewtonStepContext& context,int iteration) override;
private:
    double initial_,contraction_,minimum_,armijo_,merit_tolerance_;
    unsigned max_trials_;
    bool residual_screening_;
};
} // namespace NumLib
