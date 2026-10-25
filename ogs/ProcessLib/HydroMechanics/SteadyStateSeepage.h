// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <utility>
#include <vector>
#include "NumLib/TimeStepping/Time.h"

namespace ProcessLib::HydroMechanics
{
/// The incoming step owns its end point; the next policy starts after it.
struct SteadyStateSeepage
{
    std::vector<std::pair<double, double>> intervals;

    bool active(double const time) const
    {
        for (auto const& [start, end] : intervals)
        {
            if (NumLib::Time{time} > NumLib::Time{start} &&
                NumLib::Time{time} <= NumLib::Time{end})
            {
                return true;
            }
        }
        return false;
    }
};
}  // namespace ProcessLib::HydroMechanics
