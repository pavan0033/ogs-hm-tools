// SPDX-License-Identifier: BSD-3-Clause
#pragma once

namespace ProcessLib::HydroMechanics
{
/// Local screened-node complementarity. Rate is a mass extraction magnitude;
/// reaction is positive INTO the domain, as in HM MassFlowRate.
enum class PumpingState { Off = 0, Rate = 1, MinimumHead = 2 };

inline PumpingState pumpingState(PumpingState const old,
                                 double const pressure,
                                 double const minimum_pressure,
                                 double const reaction,
                                 double const maximum_mass_rate,
                                 double const pressure_tolerance,
                                 double const mass_tolerance)
{
    if (old == PumpingState::Rate)
    {
        return pressure < minimum_pressure - pressure_tolerance
                   ? PumpingState::MinimumHead
                   : PumpingState::Rate;
    }
    if (old == PumpingState::MinimumHead)
    {
        if (reaction > mass_tolerance) { return PumpingState::Off; }
        if (reaction < -maximum_mass_rate - mass_tolerance)
        { return PumpingState::Rate; }
        return PumpingState::MinimumHead;
    }
    return pressure > minimum_pressure + pressure_tolerance
               ? PumpingState::Rate
               : PumpingState::Off;
}
}  // namespace ProcessLib::HydroMechanics
