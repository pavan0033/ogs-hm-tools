// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace ProcessLib::StructuralElements
{
/// Analysis-time amplitude of an elastic plate's complete consistent load.
struct DistributedLoadRamp
{
    double start;
    double end;
    DistributedLoadRamp(double begin, double finish) : start(begin), end(finish)
    {
        if (!std::isfinite(start) || !std::isfinite(end) || end<=start)
        { throw std::invalid_argument("Plate distributed load ramp needs finite increasing start/end times"); }
    }
    double factor(double time) const
    { return std::clamp((time-start)/(end-start),0.,1.); }
};
}
