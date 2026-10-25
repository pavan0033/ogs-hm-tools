// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace ProcessLib::StructuralElements
{
struct ContactState { double plastic_slip=0.; };
struct ContactResponse
{
    Eigen::Vector2d traction;
    Eigen::Matrix2d tangent;
    ContactState state;
};
/// Closure is positive in compression; slip is signed. Open contact transmits
/// neither normal nor shear traction by default. RetainCohesion instead keeps
/// the cohesive shear spring/cap at zero normal pressure while open. It has no
/// damage or fracture-energy law. The Coulomb tangent is nonsymmetric.
inline ContactResponse contactTrial(Eigen::Vector2d const& deformation,
                                    ContactState const& accepted,
                                    double kn,double kt,double cohesion,double mu,
                                    bool retain_cohesion_when_open=false)
{
    if(!deformation.allFinite() || !std::isfinite(accepted.plastic_slip) ||
       !std::isfinite(kn) || kn<=0. || !std::isfinite(kt) || kt<=0. ||
       !std::isfinite(cohesion) || cohesion<0. || !std::isfinite(mu) || mu<0.)
    { throw std::invalid_argument("Invalid contact properties/deformation"); }
    if(deformation[0]<0. && (!retain_cohesion_when_open || cohesion==0.))
    { return {Eigen::Vector2d::Zero(),Eigen::Matrix2d::Zero(),{deformation[1]}}; }
    double const normal=kn*std::max(deformation[0],0.);
    double const normal_tangent=deformation[0]>=0.?kn:0.;
    double const trial=kt*(deformation[1]-accepted.plastic_slip);
    double const limit=cohesion+mu*normal;
    Eigen::Matrix2d tangent=Eigen::Matrix2d::Zero();tangent(0,0)=normal_tangent;
    // A frictionless, cohesionless interface has identically zero shear
    // traction, including at zero slip. Its derivative must also be zero.
    if(limit==0.)
    { return {{normal,0.},tangent,{deformation[1]}}; }
    if(std::abs(trial)<=limit)
    { tangent(1,1)=kt;return {{normal,trial},tangent,accepted}; }
    double const sign=trial<0.?-1.:1.;
    tangent(1,0)=sign*mu*normal_tangent;
    return {{normal,sign*limit},tangent,{deformation[1]-sign*limit/kt}};
}
}
