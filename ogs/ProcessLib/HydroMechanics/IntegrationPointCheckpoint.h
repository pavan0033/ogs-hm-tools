// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "MaterialLib/SolidModels/MFront/CheckpointData.h"
#include "MaterialLib/SolidModels/MFront/MFrontGeneric.h"
#include "MathLib/KelvinVector.h"
#include "MechanicalDrainage.h"

namespace ProcessLib::HydroMechanics
{
namespace IPCheckpointDetail
{
using namespace MaterialLib::Solids::MFront::CheckpointDetail;
inline std::string readBytes(std::istream& in)
{
    auto const size = read<std::uint64_t>(in);
    if (size > 16 * 1024 * 1024)
    { OGS_FATAL("HM checkpoint nested capsule exceeds its size limit."); }
    std::string bytes(size, '\0');
    if (!in.read(bytes.data(), size))
    { OGS_FATAL("Truncated HM checkpoint nested capsule."); }
    return bytes;
}
template <typename Vector> void writeKelvin(std::ostream& out, Vector const& vector)
{
    writeVector(out, std::vector<double>(vector.data(), vector.data() + vector.size()));
}
template <typename Vector> void readKelvin(std::istream& in, Vector& vector)
{
    std::vector<double> values(vector.size());
    readVector(in, values);
    vector = Eigen::Map<Vector const>(values.data());
    if (!vector.allFinite()) { OGS_FATAL("Nonfinite HM checkpoint stress/strain."); }
}
} // namespace IPCheckpointDetail

template <int Dim> struct DecodedHMIntegrationPoint
{
    using KV = MathLib::KelvinVector::KelvinVectorType<Dim>;
    KV stress, previous_stress, strain, previous_strain;
    double strain_rate, pressure, excess, fluid_content;
    MechanicalDrainageState drainage;
    std::unique_ptr<MaterialLib::Solids::MFront::MaterialStateVariablesMFront<Dim>> material;
};

/// Host IP capsule only. The process owner must supply accepted-step semantics,
/// model/runtime identity, nodal state, clock and boundary-reference persistence.
template <int Dim, typename IP>
std::string captureHMIntegrationPoint(IP const& ip, bool const total_stress)
{
    using namespace IPCheckpointDetail;
    using MSV = MaterialLib::Solids::MFront::MaterialStateVariablesMFront<Dim>;
    auto const* state = dynamic_cast<MSV const*>(ip.material_state_variables.get());
    if (!state || (state->_behaviour_data.s0.b.behaviour != "MohrCoulombAbboSloan" && state->_behaviour_data.s0.b.behaviour != "MohrCoulombAbboSloanTension" && state->_behaviour_data.s0.b.behaviour != "PublishedHS"))
    { OGS_FATAL("HM checkpoint currently requires MohrCoulombAbboSloan or PublishedHS MFront history."); }
    if (total_stress && state->_behaviour_data.s0.b.behaviour == "PublishedHS")
    { OGS_FATAL("PublishedHS checkpoint requires an effective-stress basis."); }
    std::ostringstream out(std::ios::out | std::ios::binary);
    write(out, static_cast<std::uint32_t>(Dim));
    write(out, static_cast<std::uint32_t>(total_stress));
    write(out, ip.integration_weight);
    writeKelvin(out, ip.sigma_eff); writeKelvin(out, ip.sigma_eff_prev);
    writeKelvin(out, ip.eps); writeKelvin(out, ip.eps_prev);
    write(out, ip.strain_rate_variable);
    write(out, ip.mechanical_pressure); write(out, ip.mechanical_excess);
    write(out, ip.mechanical_fluid_content);
    write(out, static_cast<std::uint32_t>(ip.mechanical_drainage.hasReference()));
    write(out, ip.mechanical_drainage.referenceStrain());
    write(out, ip.mechanical_drainage.referencePressure());
    write(out, ip.mechanical_drainage.referenceActualPressure());
    write(out, static_cast<std::int32_t>(state->equivalent_plastic_strain_offset_));
    write(out, state->last_local_substeps);
    writeString(out, MaterialLib::Solids::MFront::captureMFrontCheckpoint(state->_behaviour_data));
    std::string const bytes = out.str();
    std::ostringstream envelope(std::ios::out | std::ios::binary);
    envelope.write("OGSHMIP2", 8);
    write(envelope, static_cast<std::uint64_t>(bytes.size()));
    write(envelope, checksum(bytes));
    envelope.write(bytes.data(), bytes.size());
    return envelope.str();
}

/// Validate and decode into detached storage. No receiver mutation on rejection.
template <int Dim, typename IP>
DecodedHMIntegrationPoint<Dim> decodeHMIntegrationPoint(
    IP const& target, std::string_view const capsule, bool const total_stress)
{
    using namespace IPCheckpointDetail;
    using MSV = MaterialLib::Solids::MFront::MaterialStateVariablesMFront<Dim>;
    if (capsule.size() < 24 || capsule.size() > 32 * 1024 * 1024 ||
        capsule.substr(0, 8) != "OGSHMIP2")
    { OGS_FATAL("Invalid HM integration-point checkpoint header/size."); }
    std::istringstream envelope(std::string(capsule.substr(8)), std::ios::in | std::ios::binary);
    auto const size = read<std::uint64_t>(envelope);
    auto const sum = read<std::uint64_t>(envelope);
    auto const bytes = capsule.substr(24);
    if (size != bytes.size() || sum != checksum(bytes))
    { OGS_FATAL("HM integration-point checkpoint length/checksum mismatch."); }
    std::istringstream in(std::string(bytes), std::ios::in | std::ios::binary);
    if (read<std::uint32_t>(in) != Dim ||
        read<std::uint32_t>(in) != static_cast<std::uint32_t>(total_stress) ||
        read<double>(in) != target.integration_weight)
    { OGS_FATAL("HM integration-point checkpoint dimension/basis/quadrature mismatch."); }
    DecodedHMIntegrationPoint<Dim> result;
    readKelvin(in, result.stress); readKelvin(in, result.previous_stress);
    readKelvin(in, result.strain); readKelvin(in, result.previous_strain);
    result.strain_rate = read<double>(in);
    result.pressure = read<double>(in); result.excess = read<double>(in);
    result.fluid_content = read<double>(in);
    auto const present = read<std::uint32_t>(in);
    double const strain = read<double>(in);
    double const pressure = read<double>(in);
    double const actual_pressure = read<double>(in);
    if (present > 1)
    { OGS_FATAL("Invalid HM checkpoint drainage-reference flag."); }
    result.drainage = MechanicalDrainageState::fromCheckpoint(present, strain, pressure, actual_pressure);
    auto const* current = dynamic_cast<MSV const*>(target.material_state_variables.get());
    if (total_stress && current && current->_behaviour_data.s0.b.behaviour == "PublishedHS")
    { OGS_FATAL("PublishedHS checkpoint requires an effective-stress basis."); }
    auto const offset = read<std::int32_t>(in);
    double const substeps = read<double>(in);
    if (!current || (current->_behaviour_data.s0.b.behaviour != "MohrCoulombAbboSloan" && current->_behaviour_data.s0.b.behaviour != "MohrCoulombAbboSloanTension" && current->_behaviour_data.s0.b.behaviour != "PublishedHS") ||
        offset != current->equivalent_plastic_strain_offset_ ||
        !std::isfinite(substeps) || substeps < 0. || std::floor(substeps) != substeps ||
        !std::isfinite(result.strain_rate))
    { OGS_FATAL("Invalid HM checkpoint material-wrapper identity/diagnostics."); }
    for (double const value : {result.pressure, result.excess, result.fluid_content})
    {
        if (!std::isfinite(value) && !(total_stress && std::isnan(value)))
        { OGS_FATAL("Invalid HM checkpoint pressure/content for declared stress basis."); }
    }
    result.material = std::make_unique<MSV>(*current);
    MaterialLib::Solids::MFront::restoreMFrontCheckpoint(
        result.material->_behaviour_data, readBytes(in));
    result.material->last_local_substeps = substeps;
    if (in.peek() != std::char_traits<char>::eof())
    { OGS_FATAL("Trailing HM integration-point checkpoint data."); }
    return result;
}

template <int Dim, typename IP>
void commitHMIntegrationPoint(IP& target, DecodedHMIntegrationPoint<Dim>&& state)
{
    target.sigma_eff = state.stress; target.sigma_eff_prev = state.previous_stress;
    target.eps = state.strain; target.eps_prev = state.previous_strain;
    target.strain_rate_variable = state.strain_rate;
    target.mechanical_pressure = state.pressure; target.mechanical_excess = state.excess;
    target.mechanical_fluid_content = state.fluid_content;
    target.mechanical_drainage = state.drainage;
    target.material_state_variables = std::move(state.material);
}
} // namespace ProcessLib::HydroMechanics
