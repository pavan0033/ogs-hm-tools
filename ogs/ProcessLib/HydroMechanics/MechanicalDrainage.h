// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <algorithm>
#include <cmath>
#include <map>
#include <limits>
#include <optional>
#include "BaseLib/Error.h"
#include "NumLib/TimeStepping/Time.h"

namespace ProcessLib::HydroMechanics
{
enum class MechanicalDrainageType { DrainedEffectiveStress, UndrainedEffectiveStress, UndrainedTotalStress };
struct MechanicalDrainageMaterial
{
    MechanicalDrainageType type;
    double biot_modulus;
    double fluid_bulk_modulus = 0.;
    bool reference_and_excess = false;

    double appliedPressure(double actual, double excess, double reference,
                           bool positive_reference) const
    {
        if (!positive_reference) { return actual; }
        return reference_and_excess ? std::max(reference, 0.) + excess
                                    : std::max(actual, 0.);
    }
    double pressureDerivative(double actual, bool positive_reference) const
    {
        return !positive_reference || reference_and_excess || actual > 0. ? 1. : 0.;
    }


    bool usesTotalStress() const { return type == MechanicalDrainageType::UndrainedTotalStress; }

    MechanicalDrainageMaterial resolved(double const alpha, double const porosity,
                                       double const drained_bulk_modulus) const
    {
        if (usesTotalStress()) { OGS_FATAL("TSA has no Biot storage response."); }
        if (biot_modulus != 0.) { return *this; }
        auto result = *this;
        if (type == MechanicalDrainageType::DrainedEffectiveStress && fluid_bulk_modulus == 0.)
        {
            // Pressure is fixed: compressibility does not enter the response.
            result.biot_modulus = std::numeric_limits<double>::infinity();
            return result;
        }
        if (!std::isfinite(fluid_bulk_modulus) || fluid_bulk_modulus < 0. ||
            !std::isfinite(porosity) || porosity <= 0. || porosity >= 1. ||
            !std::isfinite(alpha) || alpha < porosity || alpha > 1. ||
            !std::isfinite(drained_bulk_modulus) || drained_bulk_modulus <= 0.)
        {
            OGS_FATAL("Derived mechanical storage needs positive finite fluid_bulk_modulus and drained elastic bulk modulus, and 0 < porosity <= alpha <= 1 with porosity < 1.");
        }
        // Ks = Kd/(1-alpha), avoiding infinity when alpha=1.
        double const compliance = (fluid_bulk_modulus == 0. ? 0. : porosity / fluid_bulk_modulus) +
            (alpha - porosity) * (1. - alpha) / drained_bulk_modulus;
        result.biot_modulus = 1. / compliance;
        if (!std::isfinite(result.biot_modulus) || result.biot_modulus <= 0.)
        {
            OGS_FATAL("Derived Biot modulus must be finite and positive.");
        }
        return result;
    }
};
struct MechanicalDrainageSettings
{
    double reference_time;
    std::map<int, MechanicalDrainageMaterial> materials;
    double safety_reference_time=0.;
    std::map<int,MechanicalDrainageMaterial> safety_materials;
    MechanicalDrainageMaterial const& materialAt(int const id,double const t) const
    {
        return safety_reference_time>0. && NumLib::Time{t}>NumLib::Time{safety_reference_time}
            ? safety_materials.at(id) : materials.at(id);
    }
    bool isSafetyReference(double const t) const
    { return safety_reference_time>0. && NumLib::Time{t}==NumLib::Time{safety_reference_time}; }
    bool active(double const t) const
    {
        return NumLib::Time{t} > NumLib::Time{reference_time};
    }
};
struct MechanicalDrainageResponse
{
    double pressure;
    double excess;
    double fluid_content;
    double volumetric_stiffness;
};
/// Immutable accepted reference; evaluating/rejecting trials never recaptures it.
class MechanicalDrainageState
{
public:
    void accept(double const t, double const reference_time,
                double const volumetric_strain, double const pressure)
    {
        if (!std::isfinite(reference_time) || reference_time <= 0 ||
            !std::isfinite(t) || !std::isfinite(volumetric_strain) || !std::isfinite(pressure))
        {
            OGS_FATAL("Mechanical drainage requires finite state and positive reference time.");
        }
        if (_reference || NumLib::Time{t} < NumLib::Time{reference_time})
        {
            return;
        }
        if (NumLib::Time{t} != NumLib::Time{reference_time})
        {
            OGS_FATAL("Mechanical drainage skipped the accepted reference time.");
        }
        _reference = Reference{volumetric_strain, pressure, pressure};
    }
    /// Caller is the process accepted-step hook, never a trial evaluation.
    void rebaseAccepted(double const volumetric_strain,double const primary_pressure,
                         double const actual_pressure)
    {
        if (!std::isfinite(volumetric_strain) || !std::isfinite(primary_pressure) || !std::isfinite(actual_pressure))
        { OGS_FATAL("Nonfinite accepted mechanical pressure reference."); }
        _reference=Reference{volumetric_strain,primary_pressure,actual_pressure};
    }
    MechanicalDrainageResponse evaluate(double const reference_pressure,
                                        double const volumetric_strain,
                                        double const alpha,
                                        MechanicalDrainageMaterial const& material) const
    {
        if (material.usesTotalStress()) { OGS_FATAL("TSA does not predict fluid pressure or fluid content."); }
        double const M = material.biot_modulus;
        if (!_reference || !std::isfinite(reference_pressure) ||
            !std::isfinite(volumetric_strain) || !std::isfinite(alpha) ||
            alpha <= 0 || alpha > 1 || M <= 0 ||
            (!std::isfinite(M) && !(material.type == MechanicalDrainageType::DrainedEffectiveStress && M == std::numeric_limits<double>::infinity())))
        {
            OGS_FATAL("Mechanical drainage needs an accepted reference, 0 < alpha <= 1 and finite positive Biot modulus.");
        }
        if (std::abs(reference_pressure - _reference->pressure) > 1e-7)
        {
            OGS_FATAL("Mechanical drainage primary pressure must remain at its accepted reference field.");
        }
        double const change = volumetric_strain - _reference->strain;
        bool const undrained = material.type == MechanicalDrainageType::UndrainedEffectiveStress;
        double const excess = undrained ? -alpha * M * change : 0.;
        // Preserve legacy arithmetic exactly when the two references coincide.
        bool const legacy=_reference->actual_pressure==_reference->pressure;
        double const pressure = legacy ? reference_pressure + excess : _reference->actual_pressure + excess;
        double const correction = legacy ? excess : (_reference->actual_pressure-reference_pressure)+excess;
        double const content = legacy ? (reference_pressure - _reference->pressure + excess) / M + alpha * change
                                      : excess / M + alpha * change;
        double const stiffness = undrained ? alpha * alpha * M : 0.;
        if (!std::isfinite(pressure) || !std::isfinite(content) || !std::isfinite(stiffness))
        {
            OGS_FATAL("Mechanical drainage response overflow.");
        }
        return {pressure, correction, content, stiffness};
    }
    bool hasReference() const { return _reference.has_value(); }
    static MechanicalDrainageState fromCheckpoint(
        bool const present, double const strain, double const pressure)
    { return fromCheckpoint(present,strain,pressure,pressure); }
    static MechanicalDrainageState fromCheckpoint(
        bool const present,double const strain,double const pressure,double const actual_pressure)
    {
        MechanicalDrainageState state;
        if (!std::isfinite(strain) || !std::isfinite(pressure) || !std::isfinite(actual_pressure) ||
            (!present && (strain != 0. || pressure != 0. || actual_pressure != 0.)))
        { OGS_FATAL("Invalid checkpoint mechanical drainage reference."); }
        if (present) { state._reference = Reference{strain, pressure, actual_pressure}; }
        return state;
    }
    double referenceStrain() const { return _reference ? _reference->strain : 0.; }
    double referencePressure() const { return _reference ? _reference->pressure : 0.; }
    double referenceActualPressure() const { return _reference ? _reference->actual_pressure : 0.; }
private:
    struct Reference { double strain; double pressure; double actual_pressure; };
    std::optional<Reference> _reference;
};
} // namespace ProcessLib::HydroMechanics
