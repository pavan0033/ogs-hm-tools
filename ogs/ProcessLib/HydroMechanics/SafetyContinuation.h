// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include "BaseLib/Error.h"
#include "NumLib/TimeStepping/Time.h"
namespace ProcessLib::HydroMechanics
{
struct SafetyContinuation
{
    double physical_time=0.,initialisation_end=0.;
    std::string pressure_policy;
    std::vector<double> coordinates,factors;
    bool working_hydraulic_rates=true;
    bool verify_equilibrium=false;
    bool usesActualPressure() const
    { return pressure_policy=="FrozenAcceptedActualPressure" || pressure_policy=="AdditionalUndrained" || pressure_policy=="MaterialDrainage"; }
    bool active(double const coordinate) const
    { return NumLib::Time{coordinate}>NumLib::Time{physical_time}; }
    double factor(double const coordinate) const
    {
        if (coordinate<=coordinates.front()) { return 1.; }
        for (std::size_t i=1;i<coordinates.size();++i)
        {
            if (coordinate<=coordinates[i])
            {
                double const fraction=(coordinate-coordinates[i-1])/(coordinates[i]-coordinates[i-1]);
                return factors[i-1]+fraction*(factors[i]-factors[i-1]);
            }
        }
        return factors.back();
    }
};
} // namespace ProcessLib::HydroMechanics
