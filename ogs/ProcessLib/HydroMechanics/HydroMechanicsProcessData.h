// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <utility>
#include <variant>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/VariableType.h"
#include "ParameterLib/Parameter.h"
#include "MechanicalDrainage.h"
#include "NativeCheckpoint.h"
#include "SafetyContinuation.h"
#include "InitialisationGate.h"
#include "SteadyStateSeepage.h"
#include "ExcavationRelease.h"
#include "ProcessLib/Common/HydroMechanics/InitialStress.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;
}
}  // namespace MaterialLib
namespace ProcessLib
{
namespace HydroMechanics
{

struct Monolithic
{
};

struct Staggered
{
    /// An optional input to set an algorithmic parameter of the staggered
    /// scheme. In the HydroMechanics process the fixed-stress split has been
    /// implemented as staggered scheme with a stabilization parameter to be
    /// set. For more details see [user guide -
    /// conventions](https://www.opengeosys.org/docs/userguide/basics/conventions/#fixed-stress-split-for-hydro-mechanical-processes).
    double const fixed_stress_stabilization_parameter;

    /// If set, the volume stress rate is fixed over time step, e.g.
    /// \f[\dot{\sigma}_v^{n} = \dot{\sigma}_v^{n-1}, \f]
    /// where \f$n\f$ is the time step index.
    /// Otherwise, the volume stress rate is fixed over coupling iteration, e.g.
    /// \f[\dot{\sigma}_v^{n, k} = \dot{\sigma}_v^{n, k-1}, \f]
    /// where \f$k\f$ is the coupling iteration index, and
    /// \f[ \dot{()}^{n, k} = \left(()^{n, k} - ()^{n-1}\right)/dt. \f]
    bool const fixed_stress_over_time_step;
};

using CouplingScheme = std::variant<Monolithic, Staggered>;

template <int DisplacementDim>
struct HydroMechanicsProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;

    /// The constitutive relation for the mechanical part.
    std::map<int, std::shared_ptr<
                      MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;

    InitialStress const initial_stress;

    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;

    /// If set mass lumping will be applied to the pressure equation.
    bool const mass_lumping = false;

    CouplingScheme coupling_scheme;

    /// ID of hydraulic process.
    int const hydraulic_process_id = 0;

    /// ID of the processes that contains mechanical process.
    int const mechanics_related_process_id = 0;

    const bool use_taylor_hood_elements = false;

    MaterialPropertyLib::Variable const phase_variable{};

    MeshLib::PropertyVector<double>* pressure_interpolated = nullptr;
    std::array<MeshLib::PropertyVector<double>*, 3> principal_stress_vector = {
        nullptr, nullptr, nullptr};
    MeshLib::PropertyVector<double>* principal_stress_values = nullptr;

    /// Total permeability as a symmetric tensor of length 4 or 6
    /// with elements in the order k_xx, k_yy, k_zz, k_xy, k_yz, k_xz
    MeshLib::PropertyVector<double>* permeability = nullptr;

    /// Optional exact material-state diagnostics without nodal extrapolation.
    bool output_material_state = false;

    struct ConstructionBirth { int material_id; double time; bool reset; };
    std::vector<ConstructionBirth> construction_births;
    bool adaptive_rollback=false;
    std::vector<ExcavationRelease> excavation_release;
    bool excavation_reference_assembly=false;
    bool residual_screening_assembly=false;
    MeshLib::PropertyVector<double>* excavation_release_force=nullptr;
    bool use_positive_pore_pressure_in_mechanics = false;
    SteadyStateSeepage steady_state_seepage;
    SteadyStateSeepage transient_seepage;
    std::optional<NativeCheckpointSettings> native_checkpoint;
    std::optional<SafetyContinuation> safety;
    std::optional<MechanicalDrainageSettings> mechanical_drainage;
    std::optional<InitialisationGate> initialisation_gate;

    std::vector<double> displacement_reference_times;
    std::size_t displacement_reference_next = 0;
    std::vector<double> displacement_reference_values;
    MeshLib::PropertyVector<double>* displacement_relative = nullptr;

    /// Numerical virgin-gravity ramp, shared with hydrostatic pressure/tractions.
    double initialisation_load_ramp_end = 0.;
    double initialisationLoadFactor(double const t) const
    {
        if (initialisation_load_ramp_end <= 0.) { return 1.; }
        return std::clamp(t / initialisation_load_ramp_end, 0., 1.);
    }

    bool isMonolithicSchemeUsed() const
    {
        return std::holds_alternative<Monolithic>(coupling_scheme);
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace HydroMechanics
}  // namespace ProcessLib
