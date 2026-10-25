// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Substepping.h"

#include <spdlog/fmt/ranges.h>

#include <Eigen/Core>
#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/Integrate.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"

namespace MaterialLib::Solids::MFront
{
void validateSubsteppingSettings(SubsteppingSettings const& s)
{
    if (s.initial_substeps == 0 || s.max_substeps < s.initial_substeps ||
        s.growth_factor < 2 || !std::isfinite(s.tangent_perturbation) ||
        s.tangent_perturbation <= 0)
    {
        OGS_FATAL(
            "Invalid MFront substepping settings: require 1 <= "
            "initial_substeps <= max_substeps, growth_factor >= 2, "
            "and a finite positive tangent_perturbation.");
    }
}

std::optional<SubsteppingSettings> createSubsteppingSettings(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__solid__constitutive_relation__MFront__substepping}
    auto const tree = config.getConfigSubtreeOptional("substepping");
    if (!tree)
    {
        return std::nullopt;
    }
    SubsteppingSettings s;
    //! \ogs_file_param{material__solid__constitutive_relation__MFront__substepping__initial_substeps}
    auto const initial = tree->getConfigParameter<int>("initial_substeps", 1);
    //! \ogs_file_param{material__solid__constitutive_relation__MFront__substepping__max_substeps}
    auto const maximum = tree->getConfigParameter<int>("max_substeps", 128);
    //! \ogs_file_param{material__solid__constitutive_relation__MFront__substepping__growth_factor}
    auto const growth = tree->getConfigParameter<int>("growth_factor", 2);
    if (initial < 1 || maximum < initial || growth < 2)
    {
        OGS_FATAL("Invalid MFront substepping counts.");
    }
    s.initial_substeps = initial;
    s.max_substeps = maximum;
    s.growth_factor = growth;
    //! \ogs_file_param{material__solid__constitutive_relation__MFront__substepping__tangent_perturbation}
    s.tangent_perturbation = tree->getConfigParameter<double>(
        "tangent_perturbation", s.tangent_perturbation);
    s.allow_one_sided_tangent =
        tree->getConfigParameter<bool>("allow_one_sided_tangent", false);
    validateSubsteppingSettings(s);
    return s;
}

void checkSubsteppingBehaviour(mgis::behaviour::Behaviour const& b)
{
    using B = mgis::behaviour::Behaviour;
    using V = mgis::behaviour::Variable;
    if (b.btype != B::STANDARDSTRAINBASEDBEHAVIOUR ||
        b.kinematic != B::SMALLSTRAINKINEMATIC ||
        (b.hypothesis != mgis::behaviour::Hypothesis::PLANESTRAIN &&
         b.hypothesis != mgis::behaviour::Hypothesis::TRIDIMENSIONAL) ||
        b.gradients.size() != 1 || b.thermodynamic_forces.size() != 1 ||
        b.to_blocks.size() != 1 || b.gradients[0].name != "Strain" ||
        b.gradients[0].type != V::STENSOR ||
        b.thermodynamic_forces[0].name != "Stress" ||
        b.thermodynamic_forces[0].type != V::STENSOR ||
        b.to_blocks[0].first.name != "Stress" ||
        b.to_blocks[0].second.name != "Strain" ||
        b.to_blocks[0].first.type != V::STENSOR ||
        b.to_blocks[0].second.type != V::STENSOR || b.esvs.size() > 1 ||
        (!b.esvs.empty() &&
         (b.esvs[0].name != "Temperature" || b.esvs[0].type != V::SCALAR)))
    {
        OGS_FATAL(
            "MFront substepping supports only small-strain Strain/Stress "
            "behaviours in plane strain or 3D, with a single mechanical "
            "tangent "
            "block and at most a scalar Temperature external variable. "
            "Finite-strain and coupled tangent blocks are unsupported.");
    }
}

namespace
{
using Data = mgis::behaviour::BehaviourData;
using Vector = Eigen::VectorXd;
using Matrix =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

auto view(std::vector<double> const& values)
{
    return Eigen::Map<Vector const>(values.data(), values.size());
}
auto view(std::vector<double>& values)
{
    return Eigen::Map<Vector>(values.data(), values.size());
}

void interpolate(std::vector<double> const& a,
                 std::vector<double> const& b,
                 double const fraction,
                 std::vector<double>& out)
{
    if (fraction == 1.)
    {
        out = b;
        return;
    }
    view(out) = (1. - fraction) * view(a) + fraction * view(b);
}

bool admissibleMCStart(Data const& data)
{
    auto const& b = data.s0.b;
    if (b.behaviour != "MohrCoulombAbboSloan" && b.behaviour != "MohrCoulombAbboSloanTension")
    {
        return true;
    }
    auto property = [&](std::vector<double> const& values, char const* name)
    {
        auto const offset = mgis::behaviour::getVariableOffset(b.mps, name, b.hypothesis);
        return values.at(offset);
    };
    for (auto const* state : {&data.s0, &data.s1})
    {
        auto const& values = state->material_properties;
        double const c = property(values, "Cohesion");
        double const phi = property(values, "FrictionAngle");
        double const a = property(values, "TensionCutOffParameter");
        double const transition = property(values, "TransitionAngle");
        if (!std::isfinite(c) || c < 0 || !std::isfinite(phi) || phi < 0 || phi >= 90 ||
            !std::isfinite(transition) || transition <= 0 || transition >= 30 ||
            !std::isfinite(a) || a < 0)
        {
            OGS_FATAL("Invalid Abbo-Sloan material properties (SI stress units, angle degrees required)");
        }
        if (b.behaviour == "MohrCoulombAbboSloanTension")
        {
            double const ft=property(values,"TensileStrength");
            if (!std::isfinite(ft) || ft<0.) { OGS_FATAL("TensileStrength must be finite and nonnegative (Pa)."); }
        }
        double const cap = phi == 0 ? 100. : std::max(100., .1*c/std::tan(phi*std::numbers::pi/180.));
        if (a > cap*(1+1e-12))
        {
            OGS_FATAL("Abbo-Sloan TensionCutOffParameter a={} Pa exceeds max(100 Pa, 0.1*c*cot(phi))={} Pa", a, cap);
        }
    }
    // The initial-admissibility guard is qualified for the plane-strain model.
    // The constitutive integrator remains responsible for all other hypotheses.
    if (data.s0.thermodynamic_forces.size() != 4)
    {
        return true;
    }
    auto const& f = data.s0.thermodynamic_forces;
    Eigen::Matrix3d tensor = Eigen::Matrix3d::Zero();
    tensor.diagonal() = Eigen::Vector3d(f[0], f[1], f[2]);
    tensor(0,1) = tensor(1,0) = f[3]/std::sqrt(2.);
    if (b.behaviour == "MohrCoulombAbboSloanTension")
    {
        double const ft=property(data.s0.material_properties,"TensileStrength");
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(tensor);
        double const tolerance=1e-8*std::max({1.,ft,tensor.norm(),property(data.s0.material_properties,"Cohesion")});
        if (eigen.info()!=Eigen::Success || eigen.eigenvalues().maxCoeff()>ft+tolerance) { return false; }
    }
    double const mean = tensor.trace()/3.;
    Eigen::Matrix3d const dev = tensor - mean*Eigen::Matrix3d::Identity();
    double const physical_J2 = .5*dev.squaredNorm();
    // Regularise the Lode-angle division only. A numerical J2 floor must not
    // manufacture shear strength demand at the zero-stress, zero-cohesion apex.
    double const J2 = std::max(physical_J2, 1e-14);
    double const J3 = dev.determinant();
    double const arg = std::clamp(-3*std::sqrt(3.)*J3/(2*J2*std::sqrt(J2)), -1.+1e-14, 1.-1e-14);
    double const theta = std::asin(arg)/3.;
    auto const& values = data.s0.material_properties;
    double const phi = property(values,"FrictionAngle")*std::numbers::pi/180.;
    double const tr = property(values,"TransitionAngle")*std::numbers::pi/180.;
    double const c = property(values,"Cohesion");
    double const a = property(values,"TensionCutOffParameter");
    double K;
    if (std::abs(theta) < tr)
    {
        K = std::cos(theta)-std::sin(phi)*std::sin(theta)/std::sqrt(3.);
    }
    else
    {
        double const sign = theta < 0 ? -1. : 1.;
        double const t1 = std::cos(tr)-std::sin(phi)*std::sin(tr)/std::sqrt(3.);
        double const t2 = sign*std::sin(tr)+std::sin(phi)*std::cos(tr)/std::sqrt(3.);
        double const den = 18*std::pow(std::cos(3*tr),3);
        double const B = (sign*std::sin(6*tr)*t1-6*std::cos(6*tr)*t2)/den;
        double const C = (-std::cos(3*tr)*t1-3*sign*std::sin(3*tr)*t2)/den;
        double const A = -std::sin(phi)*sign*std::sin(tr)/std::sqrt(3.)-B*sign*std::sin(3*tr)-C*std::pow(std::sin(3*tr),2)+std::cos(tr);
        K = A+B*arg+C*arg*arg;
    }
    double const yield = mean*std::sin(phi)+std::sqrt(physical_J2*K*K+a*a*std::pow(std::sin(phi),2))-c*std::cos(phi);
    double const tolerance = 1e-8*std::max({1., c, tensor.norm()});
    if (!std::isfinite(yield) || yield > tolerance)
    {
        WARN("MC_INITIAL_STATE_INADMISSIBLE: yield={} Pa, tolerance={} Pa. Substepping cannot repair an inadmissible global starting stress.",yield,tolerance);
        return false;
    }
    return true;
}

// A successful return code does not enforce plastic complementarity. Reject
// decreasing accumulated plastic flow and a purported elastic step that changes
// plastic strain. Retrying starts from the unchanged accepted state.
bool admissibleMCIncrement(Data const& data)
{
    auto const& b = data.s0.b;
    if (b.behaviour == "MohrCoulombAbboSloanTension")
    {
        for (auto const* name : {"EquivalentPlasticStrain","TensileMultiplier0","TensileMultiplier1","TensileMultiplier2"})
        {
            auto const offset=mgis::behaviour::getVariableOffset(b.isvs,name,b.hypothesis);
            if (data.s1.internal_state_variables.at(offset)-data.s0.internal_state_variables.at(offset)<-1e-12) { return false; }
        }
        return true;
    }
    if (b.behaviour != "MohrCoulombAbboSloan")
    {
        return true;
    }
    auto const lambda = mgis::behaviour::getVariableOffset(
        b.isvs, "EquivalentPlasticStrain", b.hypothesis);
    double const increment = data.s1.internal_state_variables.at(lambda) -
                             data.s0.internal_state_variables.at(lambda);
    if (increment < -1e-12)
    {
        return false;
    }
    // Thermal/property increments may contain nonmechanical strain. This
    // kinematic check is restricted to unchanged properties and temperature.
    if (std::abs(increment) <= 1e-14 &&
        data.s0.material_properties == data.s1.material_properties &&
        data.s0.external_state_variables == data.s1.external_state_variables)
    {
        auto const young_offset = mgis::behaviour::getVariableOffset(
            b.mps, "YoungModulus", b.hypothesis);
        auto const poisson_offset = mgis::behaviour::getVariableOffset(
            b.mps, "PoissonRatio", b.hypothesis);
        double const young = data.s0.material_properties.at(young_offset);
        double const poisson = data.s0.material_properties.at(poisson_offset);
        double const shear = young / (2. * (1. + poisson));
        double const lame = young * poisson / ((1. + poisson) * (1. - 2. * poisson));
        Vector const increment_strain = view(data.s1.gradients) - view(data.s0.gradients);
        Vector expected = 2. * shear * increment_strain;
        expected.head<3>().array() += lame * increment_strain.head<3>().sum();
        Vector const actual = view(data.s1.thermodynamic_forces) -
                              view(data.s0.thermodynamic_forces);
        double const tolerance = 1e-5 + 1e-8 * std::max({
            expected.norm(), view(data.s0.thermodynamic_forces).norm(),
            view(data.s1.thermodynamic_forces).norm()});
        if ((actual - expected).lpNorm<Eigen::Infinity>() > tolerance)
        {
            return false;
        }
    }
    return true;
}

std::optional<Data> integratePath(Data const& start, unsigned const n,
                                  std::function<int(Data&)> const& integrate,
                                  bool const compute_tangent = true)
{
    Data trial = start;
    // MGIS's default error buffer is shared. Never use it for local retries.
    std::array<char, 512> error{};
    trial.error_message = error.data();
    trial.dt = start.dt / n;
    // Convert each recommendation (relative to dt/n) back to the global dt.
    // OGS currently does not consume rdt, but do not erase material advice.
    double recommended_factor = std::numeric_limits<double>::infinity();
    for (unsigned i = 0; i < n; ++i)
    {
        if (i != 0)
        {
            trial.s0 = trial.s1;
            // MGIS versions supported by OGS omit density in State assignment.
            trial.s0.mass_density = trial.s1.mass_density;
        }
        // Do not use a previous global Newton iteration's trial history.
        trial.s1 = trial.s0;
        double const fraction = (static_cast<double>(i) + 1.) / n;
        interpolate(start.s0.gradients, start.s1.gradients, fraction,
                    trial.s1.gradients);
        interpolate(start.s0.material_properties, start.s1.material_properties,
                    fraction, trial.s1.material_properties);
        interpolate(start.s0.external_state_variables,
                    start.s1.external_state_variables, fraction,
                    trial.s1.external_state_variables);
        trial.s1.mass_density = (1. - fraction) * start.s0.mass_density +
                                fraction * start.s1.mass_density;
        view(trial.K).setZero();
        // For a subdivided path the caller differentiates the complete stress
        // update. Intermediate local tangents are unused and may not exist at
        // a surface corner, even when the physical path and its probes succeed.
        trial.K[0] = static_cast<double>(
            (compute_tangent && n == 1 && start.s0.b.behaviour != "MohrCoulombAbboSloanTension" &&
             start.s0.b.behaviour != "PublishedHS")
                   ? mgis::behaviour::IntegrationType::
                         INTEGRATION_CONSISTENT_TANGENT_OPERATOR
                   : mgis::behaviour::IntegrationType::
                         INTEGRATION_NO_TANGENT_OPERATOR);
        // MGIS treats rdt as input/output; initialise each independent call.
        trial.rdt = 1.;
        if (integrate(trial) != 1 || !std::isfinite(trial.rdt) ||
            trial.rdt <= 0. ||
            !view(trial.s1.thermodynamic_forces).allFinite() ||
            !view(trial.s1.internal_state_variables).allFinite() ||
            !view(trial.K).allFinite())
        {
            return std::nullopt;
        }
        if (!admissibleMCIncrement(trial))
        {
            return std::nullopt;
        }
        recommended_factor = std::min(recommended_factor, trial.rdt / n);
    }
    // The host owns the global state commit. Never return the last substep's
    // s0.
    trial.s0 = start.s0;
    trial.s0.mass_density = start.s0.mass_density;
    trial.dt = start.dt;
    trial.rdt = recommended_factor;
    trial.error_message = start.error_message;
    return trial;
}
}  // namespace

std::optional<unsigned> integrateWithSubstepping(
    Data& data, SubsteppingSettings const& settings,
    std::function<int(Data&)> const& integrate,
    bool const compute_tangent)
{
    validateSubsteppingSettings(settings);
    auto const size = data.s1.gradients.size();
    if (size == 0 || data.s0.gradients.size() != size ||
        data.s0.thermodynamic_forces.size() != size ||
        data.s1.thermodynamic_forces.size() != size ||
        data.s0.material_properties.size() !=
            data.s1.material_properties.size() ||
        data.s0.external_state_variables.size() !=
            data.s1.external_state_variables.size() ||
        data.s0.internal_state_variables.size() !=
            data.s1.internal_state_variables.size() ||
        data.K.size() != size * size)
    {
        OGS_FATAL("Unexpected MFront mechanical tangent dimensions.");
    }
    if (!std::isfinite(data.dt) || data.dt < 0.)
    {
        OGS_FATAL(
            "MFront substepping requires a finite nonnegative time increment.");
    }
    if (!admissibleMCStart(data))
    {
        return std::nullopt;
    }
    for (unsigned n = settings.initial_substeps;;)
    {
        auto result = integratePath(data, n, integrate, compute_tangent);
        if (result && compute_tangent && (n > 1 ||
                       data.s0.b.behaviour == "MohrCoulombAbboSloanTension" ||
                       data.s0.b.behaviour == "PublishedHS"))
        {
            Matrix tangent(size, size);
            bool tangent_ok = true;
            for (std::size_t j = 0; j < size; ++j)
            {
                Data plus = data;
                Data minus = data;
                double const h = settings.tangent_perturbation;
                plus.s1.gradients[j] += h;
                minus.s1.gradients[j] -= h;
                if (plus.s1.gradients[j] == data.s1.gradients[j] ||
                    minus.s1.gradients[j] == data.s1.gradients[j])
                {
                    OGS_FATAL(
                        "MFront tangent_perturbation is too small to change "
                        "strain component {}.",
                        j);
                }
                auto const upper = integratePath(plus, n, integrate);
                auto const lower = integratePath(minus, n, integrate);
                if (upper && lower)
                {
                    tangent.col(j) = (view(upper->s1.thermodynamic_forces) -
                                      view(lower->s1.thermodynamic_forces)) /
                                     (2. * h);
                }
                else if (settings.allow_one_sided_tangent && upper)
                {
                    tangent.col(j) = (view(upper->s1.thermodynamic_forces) -
                                      view(result->s1.thermodynamic_forces)) /
                                     h;
                }
                else if (settings.allow_one_sided_tangent && lower)
                {
                    tangent.col(j) = (view(result->s1.thermodynamic_forces) -
                                      view(lower->s1.thermodynamic_forces)) /
                                     h;
                }
                else
                {
                    tangent_ok = false;
                    break;
                }
            }
            if (tangent_ok && tangent.allFinite())
            {
                Eigen::Map<Matrix>(result->K.data(), size, size) = tangent;
            }
            else
            {
                result.reset();
            }
        }
        if (result)
        {
            double const density0 = result->s0.mass_density;
            double const density1 = result->s1.mass_density;
            data = std::move(*result);
            data.s0.mass_density = density0;
            data.s1.mass_density = density1;
            if (n>1) DBUG("MFront selected local substeps={} dt={}",n,data.dt);
            return n;
        }
        if (n == settings.max_substeps)
        {
            WARN(
                "MFront substepping exhausted: dt={} strain0=[{}] strain1=[{}] "
                "stress0=[{}] internal0=[{}] properties0=[{}]",
                data.dt, fmt::join(data.s0.gradients, ","),
                fmt::join(data.s1.gradients, ","),
                fmt::join(data.s0.thermodynamic_forces, ","),
                fmt::join(data.s0.internal_state_variables, ","),
                fmt::join(data.s0.material_properties, ","));
            return std::nullopt;
        }
        // Include a non-power-of-growth maximum and avoid integer overflow.
        n = n > settings.max_substeps / settings.growth_factor
                ? settings.max_substeps
                : n * settings.growth_factor;
    }
}
}  // namespace MaterialLib::Solids::MFront
