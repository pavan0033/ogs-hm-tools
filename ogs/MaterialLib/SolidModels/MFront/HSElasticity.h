// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace MaterialLib::Solids::MFront
{
/// Read-only unloading modulus of the imported PublishedHS model, in Pa.
/// The 1 kPa numerator floor is inherited from its Matriz6 routine; it is
/// stiffness regularisation, not added cohesion or tensile strength.
inline double hsUnloadingYoungModulus(double const reference_modulus,
                                     double const exponent,
                                     double const cohesion,
                                     double const friction_degrees,
                                     double const reference_pressure,
                                     double const least_compression)
{
    if (!std::isfinite(reference_modulus) || reference_modulus <= 0. ||
        !std::isfinite(exponent) || exponent < 0. ||
        !std::isfinite(cohesion) || cohesion < 0. ||
        !std::isfinite(friction_degrees) || friction_degrees <= 0. ||
        friction_degrees >= 90. || !std::isfinite(reference_pressure) ||
        reference_pressure <= 0. || !std::isfinite(least_compression))
    {
        throw std::invalid_argument("Invalid PublishedHS unloading parameters");
    }
    double const phi = friction_degrees * std::numbers::pi / 180.;
    double const shift = cohesion * std::cos(phi);
    double const numerator =
        std::max({shift + least_compression * std::sin(phi), 1000., shift});
    return reference_modulus *
           std::pow(numerator /
                        (shift + reference_pressure * std::sin(phi)),
                    exponent);
}
}  // namespace MaterialLib::Solids::MFront
