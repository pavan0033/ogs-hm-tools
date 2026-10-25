// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "BaseLib/Error.h"

namespace ProcessLib::HydroMechanics
{
namespace NativeCheckpointDetail
{
constexpr std::uint64_t limit = 256 * 1024 * 1024;
inline std::uint64_t hash(std::string_view const bytes)
{
    std::uint64_t value = 14695981039346656037ULL;
    for (unsigned char const byte : bytes) { value = (value ^ byte) * 1099511628211ULL; }
    return value;
}
inline std::string hex(std::uint64_t const value)
{ std::ostringstream out; out << std::hex << value; return out.str(); }
inline std::string fileBytes(std::filesystem::path const& path)
{
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in || in.tellg() < 0 || static_cast<std::uint64_t>(in.tellg()) > limit)
    { OGS_FATAL("Cannot read bounded HM checkpoint/model file {}", path.string()); }
    std::string value(static_cast<std::size_t>(in.tellg()), '\0');
    in.seekg(0);
    if (!in.read(value.data(), value.size()))
    { OGS_FATAL("Cannot read HM checkpoint/model file {}", path.string()); }
    return value;
}
template <typename T> void write(std::ostream& out, T const value)
{ out.write(reinterpret_cast<char const*>(&value), sizeof(value)); }
template <typename T> T read(std::istream& in)
{
    T value;
    if (!in.read(reinterpret_cast<char*>(&value), sizeof(value)))
    { OGS_FATAL("Truncated native HM checkpoint."); }
    return value;
}
inline void bytes(std::ostream& out, std::string const& value)
{ write(out, static_cast<std::uint64_t>(value.size()));out.write(value.data(), value.size()); }
inline std::string bytes(std::istream& in)
{
    auto const size=read<std::uint64_t>(in);
    if (size>limit) { OGS_FATAL("Native HM checkpoint field too large."); }
    std::string value(size,'\0');
    if (!in.read(value.data(),size)) { OGS_FATAL("Truncated native HM checkpoint field."); }
    return value;
}
inline void vector(std::ostream& out, std::vector<double> const& v)
{
    write(out, static_cast<std::uint64_t>(v.size()));
    if (!v.empty()) { out.write(reinterpret_cast<char const*>(v.data()),v.size()*sizeof(double)); }
}
inline std::vector<double> vector(std::istream& in)
{
    auto const size=read<std::uint64_t>(in);
    if (size>limit/sizeof(double)) { OGS_FATAL("Native HM checkpoint vector too large."); }
    std::vector<double> v(size);
    if (!v.empty() && !in.read(reinterpret_cast<char*>(v.data()),v.size()*sizeof(double)))
    { OGS_FATAL("Truncated native HM checkpoint vector."); }
    for (double const value:v) { if (!std::isfinite(value)) { OGS_FATAL("Nonfinite native HM checkpoint nodal/reference field."); } }
    return v;
}
inline std::string runtimeIdentity()
{
    // Linux/serial same-build first scope. Paths and file contents, never addresses.
    std::ifstream maps("/proc/self/maps");
    if (!maps) { OGS_FATAL("Native HM checkpoint requires Linux runtime identity support."); }
    std::set<std::string> paths;
    paths.insert(std::filesystem::canonical("/proc/self/exe").string());
    for (std::string line;std::getline(maps,line);)
    {
        auto const start=line.find('/');
        if (start==std::string::npos) { continue; }
        auto const path=line.substr(start);
        if (path.find(".so")!=std::string::npos && std::filesystem::is_regular_file(path))
        { paths.insert(std::filesystem::canonical(path).string()); }
    }
    std::string identity;
    for (auto const* key : {"OMP_NUM_THREADS","OGS_ASM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"})
    { auto const* value=std::getenv(key);identity+=std::string(key)+"="+(value?value:"")+"\n"; }
    for (auto const& path:paths) { identity+=path+":"+hex(hash(fileBytes(path)))+"\n"; }
    return hex(hash(identity));
}
} // namespace NativeCheckpointDetail

struct NativeCheckpointState
{
    double time=0., initialisation_time=0.;
    std::string id, model_identity, runtime_identity;
    std::vector<double> solution, initialisation_solution;
    std::vector<double> displacement_reference;
    std::uint64_t displacement_reference_next=0;
    std::vector<double> nodal_forces, mass_flow_rate;
    std::vector<std::vector<std::string>> elements;
    std::vector<std::vector<double>> boundary_states;
    std::vector<std::vector<double>> source_states;
    // Accepted ramped-excavation state; absent in legacy capsules.
    std::vector<std::vector<double>> excavation_states;
};
struct NativeCheckpointPoint
{ std::string id, file, model_identity; double time; };
struct NativeCheckpointSettings
{
    double initialisation_time;
    std::vector<NativeCheckpointPoint> points;
    std::string restart_file, restart_identity;
    std::vector<double> initialisation_solution;
    bool restored=false;
};
inline std::string encodeNativeCheckpoint(NativeCheckpointState const& state)
{
    using namespace NativeCheckpointDetail;
    std::ostringstream out(std::ios::out|std::ios::binary);
    write(out,state.time);write(out,state.initialisation_time);
    bytes(out,state.id);bytes(out,state.model_identity);bytes(out,state.runtime_identity);
    vector(out,state.solution);vector(out,state.initialisation_solution);
    vector(out,state.displacement_reference);write(out,state.displacement_reference_next);
    vector(out,state.nodal_forces);vector(out,state.mass_flow_rate);
    write(out,static_cast<std::uint64_t>(state.elements.size()));
    for (auto const& element:state.elements)
    {
        write(out,static_cast<std::uint64_t>(element.size()));
        for (auto const& ip:element) { bytes(out,ip); }
    }
    write(out,static_cast<std::uint64_t>(state.boundary_states.size()));
    for (auto const& boundary:state.boundary_states) { vector(out,boundary); }
    if(!state.source_states.empty() || !state.excavation_states.empty())
    {
        write(out,static_cast<std::uint64_t>(state.source_states.size()));
        for(auto const& source:state.source_states) { vector(out,source); }
    }
    if (!state.excavation_states.empty())
    {
        write(out, static_cast<std::uint64_t>(state.excavation_states.size()));
        for (auto const& field : state.excavation_states) { vector(out, field); }
    }
    std::string const body=out.str();
    std::ostringstream envelope(std::ios::out|std::ios::binary);
    envelope.write(!state.excavation_states.empty()?"OGSHMCP4":(state.source_states.empty()?"OGSHMCP2":"OGSHMCP3"),8);write(envelope,static_cast<std::uint64_t>(body.size()));write(envelope,hash(body));
    envelope.write(body.data(),body.size());return envelope.str();
}
inline NativeCheckpointState decodeNativeCheckpoint(std::string_view const capsule)
{
    using namespace NativeCheckpointDetail;
    if (capsule.size()<24 || capsule.size()>limit || (capsule.substr(0,8)!="OGSHMCP2" && capsule.substr(0,8)!="OGSHMCP3" && capsule.substr(0,8)!="OGSHMCP4"))
    { OGS_FATAL("Invalid native HM checkpoint header/size."); }
    std::istringstream envelope(std::string(capsule.substr(8)),std::ios::in|std::ios::binary);
    auto const size=read<std::uint64_t>(envelope);auto const checksum=read<std::uint64_t>(envelope);
    auto const body=capsule.substr(24);
    if (size!=body.size() || checksum!=hash(body)) { OGS_FATAL("Native HM checkpoint length/checksum mismatch."); }
    std::istringstream in(std::string(body),std::ios::in|std::ios::binary);
    NativeCheckpointState state;
    state.time=read<double>(in);state.initialisation_time=read<double>(in);
    if (!std::isfinite(state.time) || !std::isfinite(state.initialisation_time) ||
        state.initialisation_time<=0. ||
        state.time < state.initialisation_time - 1e-12*std::max(1.,state.initialisation_time))
    { OGS_FATAL("Invalid native HM checkpoint physical time."); }
    state.id=bytes(in);state.model_identity=bytes(in);state.runtime_identity=bytes(in);
    state.solution=vector(in);state.initialisation_solution=vector(in);
    state.displacement_reference=vector(in);state.displacement_reference_next=read<std::uint64_t>(in);
    state.nodal_forces=vector(in);state.mass_flow_rate=vector(in);
    auto const elements=read<std::uint64_t>(in);
    if (elements>body.size()/8) { OGS_FATAL("Invalid native HM checkpoint element count."); }
    state.elements.resize(elements);
    for (auto& element:state.elements)
    {
        auto const ips=read<std::uint64_t>(in);
        if (ips>body.size()/8) { OGS_FATAL("Invalid native HM checkpoint IP count."); }
        element.reserve(ips);
        for (std::uint64_t i=0;i<ips;++i) { element.push_back(bytes(in)); }
    }
    auto const boundaries=read<std::uint64_t>(in);
    if (boundaries>body.size()/8) { OGS_FATAL("Invalid boundary checkpoint count."); }
    for (std::uint64_t i=0;i<boundaries;++i) { state.boundary_states.push_back(vector(in)); }
    if(capsule.substr(0,8)=="OGSHMCP3" || capsule.substr(0,8)=="OGSHMCP4")
    {
        auto const sources=read<std::uint64_t>(in);
        if(sources>body.size()/8) { OGS_FATAL("Invalid structural checkpoint source count."); }
        for(std::uint64_t i=0;i<sources;++i) { state.source_states.push_back(vector(in)); }
    }
    if (capsule.substr(0,8)=="OGSHMCP4")
    {
        auto const fields=read<std::uint64_t>(in);
        if (fields!=4) { OGS_FATAL("Invalid excavation checkpoint field count."); }
        for (std::uint64_t i=0;i<fields;++i) { state.excavation_states.push_back(vector(in)); }
    }
    if (in.peek()!=std::char_traits<char>::eof()) { OGS_FATAL("Trailing native HM checkpoint data."); }
    return state;
}
inline void writeNativeCheckpoint(std::filesystem::path const& file, NativeCheckpointState const& state)
{
    if (std::filesystem::exists(file)) { OGS_FATAL("Refusing to overwrite HM checkpoint {}",file.string()); }
    auto const bytes=encodeNativeCheckpoint(state);
    auto const temporary=file.string()+".partial";
    if (std::filesystem::exists(temporary)) { OGS_FATAL("Refusing to overwrite partial HM checkpoint."); }
    { std::ofstream out(temporary,std::ios::binary);out.write(bytes.data(),bytes.size());out.close();
      if (!out) { OGS_FATAL("Failed to write HM checkpoint {}",file.string()); } }
    std::filesystem::rename(temporary,file);
}
} // namespace ProcessLib::HydroMechanics
