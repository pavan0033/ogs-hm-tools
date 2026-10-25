// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <algorithm>

namespace ProcessLib::HydroMechanics
{
/// Temporary mechanical equilibrium compensation after soil removal.
/// Hydraulic equations are never compensated. At phase end the load is zero.
struct ExcavationRelease
{
    double start;
    double end;

    double remaining(double const t) const
    {
        return std::clamp((end-t)/(end-start),0.,1.);
    }
};
}
