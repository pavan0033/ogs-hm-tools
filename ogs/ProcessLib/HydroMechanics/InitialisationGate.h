// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <Eigen/Core>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include "BaseLib/Error.h"
#include "ParameterLib/Parameter.h"
namespace ProcessLib::HydroMechanics
{
struct MCAdmissibilityParameters
{
    double cohesion, friction_angle, transition_angle, cutoff;
};
struct InitialisationGate
{
    double reference_time, force_tolerance, yield_tolerance;
    std::map<int, MCAdmissibilityParameters> materials;
    std::map<int, ParameterLib::Parameter<double> const*> cohesion_parameters;
};
/// Abbo-Sloan invariants check; input shear uses Kelvin convention.
inline double mcAdmissibilityResidual(Eigen::Vector4d const& sigma,
                                    MCAdmissibilityParameters const& p)
{
    if (!sigma.allFinite()) { OGS_FATAL("Initialisation gate: nonfinite MC stress."); }
    double const phi = p.friction_angle * std::numbers::pi / 180.;
    double const trans = p.transition_angle * std::numbers::pi / 180.;
    double const mean = sigma.head<3>().sum() / 3.;
    Eigen::Vector4d dev = sigma;
    dev.head<3>().array() -= mean;
    double const j2 = std::max(.5 * dev.squaredNorm(), 1e-14);
    Eigen::Matrix3d tensor = Eigen::Matrix3d::Zero();
    tensor.diagonal() = dev.head<3>();
    tensor(0,1) = tensor(1,0) = dev[3] / std::sqrt(2.);
    double const arg = std::clamp(-3.*std::sqrt(3.)*tensor.determinant() /
        (2.*std::pow(j2, 1.5)), -1.+1e-14, 1.-1e-14);
    double const theta = std::asin(arg) / 3.;
    double k;
    if (std::abs(theta) < trans)
    {
        k = std::cos(theta) - std::sin(phi)*std::sin(theta)/std::sqrt(3.);
    }
    else
    {
        double const sign = theta >= 0. ? 1. : -1.;
        double const term1 = std::cos(trans)-std::sin(phi)*std::sin(trans)/std::sqrt(3.);
        double const term2 = sign*std::sin(trans)+std::sin(phi)*std::cos(trans)/std::sqrt(3.);
        double const denominator = 18.*std::pow(std::cos(3.*trans),3);
        double const b = (sign*std::sin(6.*trans)*term1-6.*std::cos(6.*trans)*term2)/denominator;
        double const c = (-std::cos(3.*trans)*term1-3.*sign*std::sin(3.*trans)*term2)/denominator;
        double const a = -std::sin(phi)*sign*std::sin(trans)/std::sqrt(3.) -
            b*sign*std::sin(3.*trans)-c*std::pow(std::sin(3.*trans),2)+std::cos(trans);
        k = a+b*arg+c*arg*arg;
    }
    double const value = mean*std::sin(phi)+std::sqrt(j2*k*k+
        p.cutoff*p.cutoff*std::pow(std::sin(phi),2))-p.cohesion*std::cos(phi);
    if (!std::isfinite(value)) { OGS_FATAL("Initialisation gate: nonfinite yield residual."); }
    return value;
}
} // namespace ProcessLib::HydroMechanics
