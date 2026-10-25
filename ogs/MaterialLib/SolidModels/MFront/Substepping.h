// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <MGIS/Behaviour/BehaviourData.hxx>
#include <algorithm>
#include <exception>
#include <functional>
#include <optional>
#include <utility>
#include <vector>

#include "BaseLib/ConfigTree-fwd.h"

namespace MaterialLib::Solids::MFront
{
struct SubsteppingSettings
{
    unsigned initial_substeps = 1;
    unsigned max_substeps = 128;
    unsigned growth_factor = 2;
    // Absolute perturbation of a Kelvin strain component (dimensionless).
    double tangent_perturbation = 1e-7;
    // Opt-in directional derivative if exactly one whole-path probe fails.
    bool allow_one_sided_tangent = false;
};

std::optional<SubsteppingSettings> createSubsteppingSettings(
    BaseLib::ConfigTree const& config);
void validateSubsteppingSettings(SubsteppingSettings const& settings);
void checkSubsteppingBehaviour(mgis::behaviour::Behaviour const& behaviour);

/// Transactional integration: on failure, data is unchanged. On success s0
/// remains the GLOBAL starting state. Each retry and tangent probe starts from
/// that state. Numerical tangents differentiate the complete, fixed-N path.
/// The callback is injectable for testing; it returns the MGIS status code.
std::optional<unsigned> integrateWithSubstepping(
    mgis::behaviour::BehaviourData& data,
    SubsteppingSettings const& settings,
    std::function<int(mgis::behaviour::BehaviourData&)> const& integrate,
    bool compute_tangent = true);
}  // namespace MaterialLib::Solids::MFront
