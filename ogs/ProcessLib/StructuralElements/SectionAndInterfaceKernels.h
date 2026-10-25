// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <Eigen/Core>
#include <cmath>
#include <stdexcept>

#include "StructuralKernels.h"

namespace ProcessLib::StructuralElements
{
struct SectionState
{
    Eigen::Vector2d plastic_deformation = Eigen::Vector2d::Zero();
};
struct SectionResponse
{
    Eigen::Vector2d forces;
    Eigen::Matrix2d tangent;
    SectionState trial_state;
};
/// Associated perfect-plastic section: |N|/Np + |M|/Mp <= 1.
/// Deformation = [axial strain, curvature], stiffness = [EA, EI].
inline SectionResponse sectionTrial(Eigen::Vector2d const& deformation,
                                    SectionState const& accepted,
                                    Eigen::Vector2d const& stiffness,
                                    Eigen::Vector2d const& capacities)
{
    if (!deformation.allFinite() || !accepted.plastic_deformation.allFinite() ||
        !stiffness.allFinite() || !capacities.allFinite() ||
        (stiffness.array() <= 0.).any() || (capacities.array() <= 0.).any())
    {
        throw std::invalid_argument("Invalid section material input");
    }
    Eigen::Vector2d const trial =
        stiffness.cwiseProduct(deformation - accepted.plastic_deformation);
    if ((trial.cwiseAbs().cwiseQuotient(capacities)).sum() <= 1.)
    {
        return {trial, stiffness.asDiagonal(), accepted};
    }
    Eigen::Vector2d const magnitude = trial.cwiseAbs();
    Eigen::Vector2d const weight = stiffness.cwiseQuotient(capacities);
    double lower = 0.;
    double upper = magnitude.cwiseQuotient(weight).maxCoeff();
    for (int i = 0; i < 100; ++i)
    {
        double const lambda = (lower + upper) / 2.;
        double const yield = (magnitude - lambda * weight)
                                 .cwiseMax(0.).cwiseQuotient(capacities).sum();
        if (yield > 1.) { lower = lambda; } else { upper = lambda; }
    }
    Eigen::Vector2d const returned = (magnitude - upper * weight).cwiseMax(0.);
    Eigen::Vector2d normal = Eigen::Vector2d::Zero();
    Eigen::Matrix2d tangent = Eigen::Matrix2d::Zero();
    Eigen::Vector2d force;
    for (int i = 0; i < 2; ++i)
    {
        double const sign = trial[i] < 0. ? -1. : 1.;
        force[i] = sign * returned[i];
        if (returned[i] > 0.)
        {
            normal[i] = sign / capacities[i];
            tangent(i, i) = stiffness[i];
        }
    }
    Eigen::Vector2d const direction = tangent * normal;
    tangent -= direction * direction.transpose() / normal.dot(direction);
    return {force, tangent,
            SectionState{deformation - force.cwiseQuotient(stiffness)}};
}

struct InterfaceState { double plastic_slip = 0.; };
struct InterfaceResponse
{
    double resistance;
    double tangent;
    InterfaceState trial_state;
};
/// Relative displacement law; stiffness/strength units must match force or
/// distributed resistance. Positive slip is compressive penetration for a base.
inline InterfaceResponse slipTrial(double slip, InterfaceState const& accepted,
                                   double stiffness, double resistance_limit,
                                   bool compression_only)
{
    auto const response = axialTrial(slip, AxialState{accepted.plastic_slip},
                                     stiffness, resistance_limit,
                                     resistance_limit, compression_only);
    return {response.force, response.tangent,
            InterfaceState{response.trial_state.plastic_strain}};
}
}  // namespace ProcessLib::StructuralElements
