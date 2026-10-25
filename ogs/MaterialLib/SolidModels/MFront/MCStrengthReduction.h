// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/BehaviourData.hxx>
#include <MGIS/Behaviour/Variable.hxx>
#include <algorithm>
#include <cmath>
#include <numbers>
#include "BaseLib/Error.h"

namespace MaterialLib::Solids::MFront
{
/// Scale a detached trial's freshly evaluated base MC properties. RF is a
/// dimensionless continuation coordinate, independent of physical time.
/// Caller supplies the previous/current accepted/trial factors explicitly.
inline void applyMCStrengthReduction(mgis::behaviour::BehaviourData& data,
                                     double const previous_factor,
                                     double const current_factor)
{
    if (data.s0.b.behaviour!="MohrCoulombAbboSloan" ||
        !std::isfinite(previous_factor) || !std::isfinite(current_factor) ||
        previous_factor<1. || current_factor<previous_factor)
    { OGS_FATAL("MC strength reduction requires MC and finite nondecreasing factors >=1."); }
    auto reduce=[](auto const& state, double const factor)
    {
        auto properties=state.material_properties;
        auto offset=[&](char const* name) {
            return mgis::behaviour::getVariableOffset(state.b.mps,name,state.b.hypothesis);
        };
        auto& cohesion=properties.at(offset("Cohesion"));
        auto& phi=properties.at(offset("FrictionAngle"));
        auto& psi=properties.at(offset("DilatancyAngle"));
        double const cutoff=properties.at(offset("TensionCutOffParameter"));
        if (!std::isfinite(cohesion) || cohesion<=0. || !std::isfinite(phi) ||
            phi<0. || phi>=90. || !std::isfinite(psi) || psi<0. || psi>phi ||
            !std::isfinite(cutoff) || cutoff<0.)
        { OGS_FATAL("MC strength reduction needs positive cohesion, 0<=psi<=phi<90 degrees and finite nonnegative cutoff."); }
        constexpr double radians=std::numbers::pi/180.;
        double const cap=phi==0.?100.:std::min(100.,.1*cohesion/std::tan(phi*radians));
        if (cutoff>cap*(1.+1e-12)) { OGS_FATAL("MC strength reduction cutoff exceeds user cap in Pa."); }
        if (factor==1.) { return properties; }
        cohesion/=factor;
        phi=std::atan(std::tan(phi*radians)/factor)/radians;
        psi=std::atan(std::tan(psi*radians)/factor)/radians;
        return properties;
    };
    // Validate and compute both states before modifying either trial buffer.
    auto previous=reduce(data.s0,previous_factor);
    auto current=reduce(data.s1,current_factor);
    data.s0.material_properties=std::move(previous);
    data.s1.material_properties=std::move(current);
}
} // namespace MaterialLib::Solids::MFront
