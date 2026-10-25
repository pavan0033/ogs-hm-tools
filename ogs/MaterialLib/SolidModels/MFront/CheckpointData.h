// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/BehaviourData.hxx>
#include <bit>
#include <cstdint>
#include <cstring>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include "BaseLib/Error.h"

namespace MaterialLib::Solids::MFront
{
namespace CheckpointDetail
{
inline std::uint64_t checksum(std::string_view const bytes)
{
    std::uint64_t result = 14695981039346656037ULL;
    for (unsigned char const byte : bytes) { result = (result ^ byte) * 1099511628211ULL; }
    return result;
}
template <typename T> void write(std::ostream& out, T const value)
{
    static_assert(std::is_trivially_copyable_v<T>);
    out.write(reinterpret_cast<char const*>(&value), sizeof(T));
}
template <typename T> T read(std::istream& in)
{
    static_assert(std::is_trivially_copyable_v<T>);
    T value;
    if (!in.read(reinterpret_cast<char*>(&value), sizeof(T)))
    { OGS_FATAL("Truncated MFront checkpoint data."); }
    return value;
}
inline void writeString(std::ostream& out, std::string const& value)
{
    write(out, static_cast<std::uint64_t>(value.size()));
    out.write(value.data(), value.size());
}
inline void matchString(std::istream& in, std::string const& expected)
{
    auto const size = read<std::uint64_t>(in);
    if (size != expected.size()) { OGS_FATAL("MFront checkpoint behaviour/build identity mismatch."); }
    std::string value(size, '\0');
    if (!in.read(value.data(), size) || value != expected)
    { OGS_FATAL("MFront checkpoint behaviour/build identity mismatch."); }
}
inline void writeVector(std::ostream& out, std::vector<double> const& values)
{
    write(out, static_cast<std::uint64_t>(values.size()));
    if (!values.empty()) { out.write(reinterpret_cast<char const*>(values.data()), values.size()*sizeof(double)); }
}
inline void readVector(std::istream& in, std::vector<double>& values)
{
    if (read<std::uint64_t>(in) != values.size())
    { OGS_FATAL("MFront checkpoint vector layout mismatch."); }
    if (!values.empty() && !in.read(reinterpret_cast<char*>(values.data()), values.size()*sizeof(double)))
    { OGS_FATAL("Truncated MFront checkpoint vector."); }
}
inline void writeState(std::ostream& out, mgis::behaviour::State const& state)
{
    write(out, state.stored_energy); write(out, state.dissipated_energy); write(out, state.mass_density);
    writeVector(out, state.gradients); writeVector(out, state.thermodynamic_forces);
    writeVector(out, state.material_properties); writeVector(out, state.internal_state_variables);
    writeVector(out, state.external_state_variables);
}
inline void readState(std::istream& in, mgis::behaviour::State& state)
{
    state.stored_energy=read<double>(in); state.dissipated_energy=read<double>(in); state.mass_density=read<double>(in);
    readVector(in, state.gradients); readVector(in, state.thermodynamic_forces);
    readVector(in, state.material_properties); readVector(in, state.internal_state_variables);
    readVector(in, state.external_state_variables);
}
} // namespace CheckpointDetail

/// Exact MGIS data capsule for the same binary/schema/platform only. The host
/// must additionally fingerprint its loaded libraries and capture FE/BC history.
/// The transient error-buffer pointer is deliberately owned by the receiver.
inline std::string captureMFrontCheckpoint(mgis::behaviour::BehaviourData const& data)
{
    using namespace CheckpointDetail;
    static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559);
    std::ostringstream payload(std::ios::out|std::ios::binary);
    auto const& b=data.s0.b;
    write(payload, static_cast<std::uint32_t>(std::endian::native == std::endian::little ? 1 : 2));
    writeString(payload,b.library);writeString(payload,b.behaviour);writeString(payload,b.build_id);
    writeString(payload,b.tfel_version);writeString(payload,b.function);
    write(payload,static_cast<std::int32_t>(b.hypothesis));
    write(payload,data.dt);write(payload,data.rdt);write(payload,data.speed_of_sound);writeVector(payload,data.K);
    writeState(payload,data.s0);writeState(payload,data.s1);
    std::string const bytes=payload.str();
    std::ostringstream capsule(std::ios::out|std::ios::binary);
    capsule.write("OGSMGIS1",8);write(capsule,static_cast<std::uint64_t>(bytes.size()));write(capsule,checksum(bytes));
    capsule.write(bytes.data(),bytes.size());return capsule.str();
}

/// Parse into a detached trial; rejection never mutates the receiving state.
inline void restoreMFrontCheckpoint(mgis::behaviour::BehaviourData& target,
                                   std::string_view const capsule)
{
    using namespace CheckpointDetail;
    if (capsule.size()<24 || capsule.substr(0,8)!="OGSMGIS1")
    { OGS_FATAL("Invalid MFront checkpoint header."); }
    std::istringstream envelope(std::string(capsule.substr(8)),std::ios::in|std::ios::binary);
    auto const size=read<std::uint64_t>(envelope);auto const sum=read<std::uint64_t>(envelope);
    auto const bytes=capsule.substr(24);
    if (size!=bytes.size() || size>16*1024*1024 || checksum(bytes)!=sum)
    { OGS_FATAL("MFront checkpoint length/checksum mismatch."); }
    std::istringstream in(std::string(bytes),std::ios::in|std::ios::binary);
    if (read<std::uint32_t>(in)!=(std::endian::native==std::endian::little ? 1u : 2u))
    { OGS_FATAL("MFront checkpoint platform mismatch."); }
    auto const& b=target.s0.b;
    matchString(in,b.library);matchString(in,b.behaviour);matchString(in,b.build_id);
    matchString(in,b.tfel_version);matchString(in,b.function);
    if (read<std::int32_t>(in)!=static_cast<std::int32_t>(b.hypothesis))
    { OGS_FATAL("MFront checkpoint hypothesis mismatch."); }
    auto trial=target;
    trial.dt=read<double>(in);trial.rdt=read<double>(in);trial.speed_of_sound=read<double>(in);readVector(in,trial.K);
    readState(in,trial.s0);readState(in,trial.s1);
    if (in.peek()!=std::char_traits<char>::eof()) { OGS_FATAL("MFront checkpoint has trailing data."); }
    char* const error_buffer=target.error_message;
    target=trial;
    // Supported MGIS assignment omits State::mass_density.
    target.s0.mass_density=trial.s0.mass_density;target.s1.mass_density=trial.s1.mass_density;
    target.error_message=error_buffer;
}
} // namespace MaterialLib::Solids::MFront
