// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MaterialLib/SolidModels/MFront/Substepping.h"

#include <gtest/gtest.h>

#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/Integrate.hxx>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>

#include "BaseLib/ConfigTree.h"
#include "ProcessLib/HydroMechanics/ExcavationRelease.h"
#include "ParameterLib/ConstantParameter.h"
#include "MaterialLib/SolidModels/MFront/MFrontGeneric.h"
#include "Tests/TestTools.h"
#include "BaseLib/Logging.h"
#include "MaterialLib/SolidModels/MFront/HSElasticity.h"
#include "ProcessLib/StructuralElements/Contact2D.h"
#include "ProcessLib/StructuralElements/DistributedLoadRamp.h"
#include "MathLib/LinAlg/Eigen/EigenLinearSolver.h"
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"

TEST(HardeningSoilMC, UnloadingModulusAtReferenceConfinement)
{
    using MaterialLib::Solids::MFront::hsUnloadingYoungModulus;
    EXPECT_DOUBLE_EQ(hsUnloadingYoungModulus(9e7, .5, 0., 34., 1e5, 1e5), 9e7);
    EXPECT_NEAR(hsUnloadingYoungModulus(9e7, .5, 0., 34., 1e5, 4e5), 1.8e8, 1e-7);
}

TEST(HardeningSoilMC, StiffnessFloorDoesNotRequireCohesion)
{
    using MaterialLib::Solids::MFront::hsUnloadingYoungModulus;
    EXPECT_GT(hsUnloadingYoungModulus(9e7, .5, 0., 34., 1e5, 0.), 0.);
    EXPECT_DOUBLE_EQ(hsUnloadingYoungModulus(9e7, .5, 0., 34., 1e5, 0.),
                     hsUnloadingYoungModulus(9e7, .5, 0., 34., 1e5, -100.));
    EXPECT_THROW(hsUnloadingYoungModulus(9e7, .5, -1., 34., 1e5, 0.), std::invalid_argument);
}

TEST(HMLinearSolver, SparseLUAcceptsConstrainedZeroCouplingPattern)
{
    MathLib::EigenMatrix matrix(3);
    for (int i = 0; i < 3; ++i) { matrix.setValue(i, i, 1.); }
    matrix.setValue(0, 2, 0.);
    matrix.setValue(2, 1, 0.);
    MathLib::EigenVector rhs(3), solution(3);
    for (int i = 0; i < 3; ++i) { rhs[i] = i + 1.; }
    solution.setZero();
    MathLib::EigenOption option;
    option.solver_type = MathLib::EigenOption::SolverType::SparseLU;
    MathLib::EigenLinearSolver solver("zero_coupling", option);
    ASSERT_TRUE(solver.solve(matrix, rhs, solution));
    for (int i = 0; i < 3; ++i) { EXPECT_DOUBLE_EQ(solution[i], rhs[i]); }
}

TEST(StructuralContact, FrictionlessZeroSlipHasZeroShearTangent)
{
    using ProcessLib::StructuralElements::contactTrial;
    Eigen::Vector2d const u(1e-4, 0.);
    auto const response = contactTrial(u, {}, 2e6, 1e6, 0., 0.);
    EXPECT_DOUBLE_EQ(response.traction[0], 200.);
    EXPECT_DOUBLE_EQ(response.traction[1], 0.);
    EXPECT_DOUBLE_EQ(response.tangent(1, 1), 0.);
    double const h = 1e-7;
    auto plus = u;
    auto minus = u;
    plus[1] += h;
    minus[1] -= h;
    Eigen::Vector2d const derivative =
        (contactTrial(plus, {}, 2e6, 1e6, 0., 0.).traction -
         contactTrial(minus, {}, 2e6, 1e6, 0., 0.).traction) /
        (2 * h);
    EXPECT_DOUBLE_EQ(derivative[1], 0.);
}

namespace MFrontSubsteppingTests
{
using namespace MaterialLib::Solids::MFront;
using Data = mgis::behaviour::BehaviourData;

struct MFrontSubstepping : testing::Test
{
    mgis::behaviour::Behaviour behaviour = mgis::behaviour::load(
        "libOgsMFrontBehaviour", "StandardElasticityBrick",
        mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data data{behaviour};
    SubsteppingSettings settings;

    MFrontSubstepping()
    {
        data.dt = 1.;
        data.s0.gradients.assign(4, 0.);
        data.s1.gradients = {.8, .2, 0., 0.};
        data.s0.thermodynamic_forces.assign(4, 0.);
        data.s1.thermodynamic_forces.assign(4, -999.);
        data.s0.internal_state_variables = {0.};
        data.s1.internal_state_variables = {-999.};
        settings.max_substeps = 8;
    }

    // History-dependent nonlinear update, with a nonsymmetric tangent.
    // Failure deliberately poisons trial data to test transactional rollback.
    static int update(Data& d, double limit = .25)
    {
        double const de = d.s1.gradients[0] - d.s0.gradients[0];
        EXPECT_TRUE(d.K[0] == 0. || d.K[0] == 4.);
        EXPECT_DOUBLE_EQ(d.s1.internal_state_variables[0],
                         d.s0.internal_state_variables[0]);
        if (std::abs(de) > limit)
        {
            d.s0.internal_state_variables[0] = -100.;
            d.s1.internal_state_variables[0] = -100.;
            d.s1.thermodynamic_forces[0] = -100.;
            return -1;
        }
        double const q = d.s0.internal_state_variables[0];
        double const next = q + de * (1. + q) * (1. + q);
        d.s1.internal_state_variables[0] = next;
        d.s1.thermodynamic_forces = {next, 2. * next + 3. * d.s1.gradients[1],
                                     d.s1.gradients[2], d.s1.gradients[3]};
        std::fill(d.K.begin(), d.K.end(), 0.);
        d.K[0] = (1. + q) * (1. + q);
        d.K[4] = 2. * d.K[0];
        d.K[5] = 3.;
        d.K[10] = d.K[15] = 1.;
        return 1;
    }
};

TEST_F(MFrontSubstepping, RollbackAndWholeIncrementTangent)
{
    auto const count = integrateWithSubstepping(
        data, settings, [](Data& d) { return update(d); });
    ASSERT_EQ(count, 4u);
    double q = 0.;
    double tangent = 0.;
    double last = 0.;
    for (unsigned i = 0; i < *count; ++i)
    {
        last = (1. + q) * (1. + q);
        tangent = (1. + 2. * .2 * (1. + q)) * tangent + last / *count;
        q += .2 * last;
    }
    EXPECT_NEAR(data.s1.internal_state_variables[0], q, 1e-13);
    EXPECT_NEAR(data.K[0], tangent, 1e-7);
    EXPECT_NEAR(data.K[4], 2. * tangent, 1e-7);
    EXPECT_NEAR(data.K[5], 3., 1e-7);
    EXPECT_DOUBLE_EQ(data.K[1], 0.);
    EXPECT_GT(std::abs(data.K[0] - last), .1);
    EXPECT_DOUBLE_EQ(data.s0.internal_state_variables[0], 0.);
    EXPECT_DOUBLE_EQ(data.s0.gradients[0], 0.);
    EXPECT_DOUBLE_EQ(data.dt, 1.);
}

TEST_F(MFrontSubstepping, WholePathDoesNotRequireIntermediateTangents)
{
    auto stress_only = [](Data& d)
    {
        if (d.K[0] == 4.)
        {
            d.s0.internal_state_variables[0] = -100.;
            return -1;
        }
        return update(d, 1.);
    };
    ASSERT_EQ(integrateWithSubstepping(data, settings, stress_only), 2u);
    EXPECT_NEAR(data.s1.internal_state_variables[0], 1.184, 1e-13);
    EXPECT_NEAR(data.K[0], 2.04, 1e-7);
    EXPECT_NEAR(data.K[4], 4.08, 1e-7);
    EXPECT_NEAR(data.K[5], 3., 1e-7);
    EXPECT_DOUBLE_EQ(data.s0.internal_state_variables[0], 0.);
}

TEST_F(MFrontSubstepping, SingleStepStillRequestsLocalTangent)
{
    data.s1.gradients[0] = .1;
    auto local = [](Data& d)
    {
        EXPECT_DOUBLE_EQ(d.K[0], 4.);
        return update(d);
    };
    ASSERT_EQ(integrateWithSubstepping(data, settings, local), 1u);
    EXPECT_DOUBLE_EQ(data.K[0], 1.);
}

TEST_F(MFrontSubstepping, FailureLeavesCallerUnchanged)
{
    Data const saved = data;
    settings.max_substeps = 2;
    EXPECT_FALSE(integrateWithSubstepping(data, settings,
                                          [](Data& d) { return update(d); }));
    EXPECT_EQ(data.s0.internal_state_variables,
              saved.s0.internal_state_variables);
    EXPECT_EQ(data.s1.internal_state_variables,
              saved.s1.internal_state_variables);
    EXPECT_EQ(data.s1.thermodynamic_forces, saved.s1.thermodynamic_forces);
    EXPECT_EQ(data.K, saved.K);
    EXPECT_EQ(data.dt, saved.dt);
}

TEST_F(MFrontSubstepping, FailedPerturbationRetriesEntirePath)
{
    settings.initial_substeps = 4;
    auto const count = integrateWithSubstepping(
        data, settings, [](Data& d) { return update(d, .20000000001); });
    EXPECT_EQ(count, 8u);
    EXPECT_DOUBLE_EQ(data.s0.internal_state_variables[0], 0.);
}

TEST_F(MFrontSubstepping, NonPowerMaximumIsTried)
{
    settings.max_substeps = 3;
    EXPECT_EQ(integrateWithSubstepping(data, settings,
                                       [](Data& d) { return update(d, .3); }),
              3u);
}

TEST_F(MFrontSubstepping, OneStepKeepsNativeTangent)
{
    data.s1.gradients[0] = .1;
    EXPECT_EQ(integrateWithSubstepping(data, settings,
                                       [](Data& d) { return update(d); }),
              1u);
    EXPECT_DOUBLE_EQ(data.K[0], 1.);
    EXPECT_DOUBLE_EQ(data.K[4], 2.);
}

TEST_F(MFrontSubstepping, RealElasticBehaviourAndDriverInterpolation)
{
    checkSubsteppingBehaviour(behaviour);
    data.s0.internal_state_variables =
        mgis::behaviour::State(behaviour).internal_state_variables;
    data.s1.internal_state_variables = data.s0.internal_state_variables;
    data.s0.material_properties = {1e7, .3};
    data.s1.material_properties = {1e7, .3};
    data.s0.external_state_variables = {293.15};
    data.s1.external_state_variables = {303.15};
    Data direct = data;
    direct.s1 = direct.s0;
    direct.s1.gradients = data.s1.gradients;
    direct.s1.external_state_variables = data.s1.external_state_variables;
    direct.K[0] = 4.;
    auto integrate = [this](Data& d)
    {
        auto v = mgis::behaviour::make_view(d);
        return mgis::behaviour::integrate(v, behaviour);
    };
    ASSERT_EQ(integrate(direct), 1);
    settings.initial_substeps = 4;
    ASSERT_EQ(integrateWithSubstepping(data, settings, integrate), 4u);
    for (std::size_t i = 0; i < data.K.size(); ++i)
    {
        EXPECT_NEAR(data.K[i], direct.K[i], .1);
    }
    EXPECT_EQ(data.s0.external_state_variables[0], 293.15);
    EXPECT_EQ(data.s1.external_state_variables[0], 303.15);
    EXPECT_NEAR(data.s1.thermodynamic_forces[0],
                direct.s1.thermodynamic_forces[0], 1e-7);
}

TEST_F(MFrontSubstepping, RejectCoupledAndFiniteStrainBehaviours)
{
    auto coupled = behaviour;
    coupled.to_blocks.push_back(coupled.to_blocks[0]);
    EXPECT_ANY_THROW(checkSubsteppingBehaviour(coupled));
    auto finite = behaviour;
    finite.btype = mgis::behaviour::Behaviour::STANDARDFINITESTRAINBEHAVIOUR;
    EXPECT_ANY_THROW(checkSubsteppingBehaviour(finite));
}

TEST(MFrontSubsteppingConfig, ParseAndValidate)
{
    auto tree = Tests::readXml(R"(<substepping>
      <initial_substeps>2</initial_substeps><max_substeps>16</max_substeps>
      <growth_factor>3</growth_factor><tangent_perturbation>2e-7</tangent_perturbation>
    </substepping>)");
    BaseLib::ConfigTree config(std::move(tree), "",
                               BaseLib::ConfigTree::onerror,
                               BaseLib::ConfigTree::onwarning);
    auto const s = createSubsteppingSettings(config);
    ASSERT_TRUE(s);
    EXPECT_EQ(s->initial_substeps, 2u);
    EXPECT_EQ(s->max_substeps, 16u);
    EXPECT_EQ(s->growth_factor, 3u);
    EXPECT_EQ(s->tangent_perturbation, 2e-7);
    auto invalid = *s;
    invalid.growth_factor = 1;
    EXPECT_ANY_THROW(validateSubsteppingSettings(invalid));
    invalid = *s;
    invalid.tangent_perturbation = 0;
    EXPECT_ANY_THROW(validateSubsteppingSettings(invalid));
}


TEST_F(MFrontSubstepping, InterpolationHistoryDensityAndTimeRecommendation)
{
    settings.initial_substeps = settings.max_substeps = 4;
    data.s0.material_properties = {10., 20.};
    data.s1.material_properties = {30., 60.};
    data.s0.external_state_variables = {290.};
    data.s1.external_state_variables = {310.};
    data.s0.mass_density = 100.;
    data.s1.mass_density = 200.;
    data.s0.stored_energy = 3.;
    data.s0.dissipated_energy = 4.;
    unsigned calls = 0;
    auto const count = integrateWithSubstepping(
        data, settings,
        [&](Data& d)
        {
            unsigned const i = calls++ % 4;
            EXPECT_DOUBLE_EQ(d.dt, .25);
            EXPECT_DOUBLE_EQ(d.rdt, 1.);
            EXPECT_DOUBLE_EQ(d.s0.mass_density, 100. + 25. * i);
            EXPECT_DOUBLE_EQ(d.s1.mass_density, 100. + 25. * (i + 1));
            EXPECT_DOUBLE_EQ(d.s0.material_properties[0], 10. + 5. * i);
            EXPECT_DOUBLE_EQ(d.s1.material_properties[0], 15. + 5. * i);
            EXPECT_DOUBLE_EQ(d.s1.material_properties[1], 30. + 10. * i);
            EXPECT_DOUBLE_EQ(d.s0.external_state_variables[0], 290. + 5. * i);
            EXPECT_DOUBLE_EQ(d.s1.external_state_variables[0], 295. + 5. * i);
            EXPECT_DOUBLE_EQ(d.s0.internal_state_variables[0], i);
            EXPECT_DOUBLE_EQ(d.s1.internal_state_variables[0], i);
            EXPECT_DOUBLE_EQ(d.s0.stored_energy, 3. + i);
            EXPECT_DOUBLE_EQ(d.s0.dissipated_energy, 4. + 2. * i);
            d.s1.internal_state_variables[0] = i + 1;
            d.s1.stored_energy = 4. + i;
            d.s1.dissipated_energy = 6. + 2. * i;
            d.s1.thermodynamic_forces = d.s1.gradients;
            std::fill(d.K.begin(), d.K.end(), 0.);
            d.rdt = 2. + i;
            return 1;
        });
    ASSERT_EQ(count, 4u);
    EXPECT_EQ(calls, 36u);
    EXPECT_DOUBLE_EQ(data.s0.mass_density, 100.);
    EXPECT_DOUBLE_EQ(data.s1.mass_density, 200.);
    EXPECT_DOUBLE_EQ(data.s0.internal_state_variables[0], 0.);
    EXPECT_DOUBLE_EQ(data.s1.internal_state_variables[0], 4.);
    EXPECT_DOUBLE_EQ(data.rdt, .5);
    EXPECT_DOUBLE_EQ(data.dt, 1.);
}

TEST_F(MFrontSubstepping, SingleStepRecommendationAndUnreliableResult)
{
    data.s1.gradients[0] = .1;
    ASSERT_EQ(integrateWithSubstepping(data, settings,
                                       [](Data& d)
                                       {
                                           int const status = update(d);
                                           d.rdt = 1.5;
                                           return status;
                                       }),
              1u);
    EXPECT_DOUBLE_EQ(data.rdt, 1.5);
    Data const original = data;
    unsigned calls = 0;
    EXPECT_FALSE(integrateWithSubstepping(
        data, settings,
        [&](Data& d)
        {
            ++calls;
            d.rdt = .5;
            d.s1.internal_state_variables[0] = -42.;
            return 0;  // Successful but unreliable: must not be accepted.
        }));
    EXPECT_EQ(calls, 4u);  // N=1,2,4,8: first substep fails each time.
    EXPECT_EQ(data.s1.internal_state_variables,
              original.s1.internal_state_variables);
    EXPECT_EQ(data.rdt, original.rdt);
}

TEST_F(MFrontSubstepping, NonFiniteStressAndRecommendationFailTransactionally)
{
    for (bool const bad_stress : {false, true})
    {
        Data const original = data;
        EXPECT_FALSE(integrateWithSubstepping(
            data, settings,
            [&](Data& d)
            {
                if (bad_stress)
                    d.s1.thermodynamic_forces[0] =
                        std::numeric_limits<double>::quiet_NaN();
                else
                    d.rdt = std::numeric_limits<double>::quiet_NaN();
                return 1;
            }));
        EXPECT_EQ(data.s1.thermodynamic_forces,
                  original.s1.thermodynamic_forces);
        EXPECT_EQ(data.K, original.K);
    }
}

TEST(MFrontSubsteppingElastic, NativeEquivalenceAndKelvinShearIn2DAnd3D)
{
    for (auto const hypothesis : {mgis::behaviour::Hypothesis::PLANESTRAIN,
                                  mgis::behaviour::Hypothesis::TRIDIMENSIONAL})
    {
        auto b = mgis::behaviour::load("libOgsMFrontBehaviour",
                                       "StandardElasticityBrick", hypothesis);
        checkSubsteppingBehaviour(b);
        Data original{b};
        original.dt = .25;
        original.s0.material_properties =
            original.s1.material_properties = {1e7, .3};
        original.s0.external_state_variables =
            original.s1.external_state_variables = {293.15};
        for (std::size_t j = 0; j < original.s1.gradients.size(); ++j)
            original.s1.gradients[j] = (j + 1.) * 1e-4;
        original.K[0] =
            static_cast<double>(mgis::behaviour::IntegrationType::
                                    INTEGRATION_CONSISTENT_TANGENT_OPERATOR);
        auto integrate = [&b](Data& d)
        {
            auto v = mgis::behaviour::make_view(d);
            return mgis::behaviour::integrate(v, b);
        };
        Data direct = original;
        ASSERT_EQ(integrate(direct), 1);
        for (unsigned const n : {1u, 4u})
        {
            Data wrapped = original;
            SubsteppingSettings settings;
            settings.initial_substeps = settings.max_substeps = n;
            ASSERT_EQ(integrateWithSubstepping(wrapped, settings, integrate),
                      n);
            for (std::size_t i = 0; i < direct.K.size(); ++i)
                EXPECT_NEAR(wrapped.K[i], direct.K[i], 1e-3);
            for (std::size_t i = 0; i < direct.s1.gradients.size(); ++i)
            {
                EXPECT_NEAR(wrapped.s1.thermodynamic_forces[i],
                            direct.s1.thermodynamic_forces[i], 1e-8);
                if (i >= 3)
                {
                    // Kelvin strain=sqrt(2)*epsilon_ij; Kelvin
                    // stress=sqrt(2)*sigma_ij.
                    EXPECT_NEAR(wrapped.K[i * direct.s1.gradients.size() + i],
                                1e7 / 1.3, 1e-3);
                }
            }
            if (n == 1)
            {
                EXPECT_EQ(wrapped.K, direct.K);
                EXPECT_EQ(wrapped.s1.internal_state_variables,
                          direct.s1.internal_state_variables);
                EXPECT_EQ(wrapped.rdt, direct.rdt);
            }
        }
    }
}

TEST_F(MFrontSubstepping, RejectUnsupportedHypothesisAndExternalVariables)
{
    auto b = behaviour;
    b.hypothesis = mgis::behaviour::Hypothesis::PLANESTRESS;
    EXPECT_ANY_THROW(checkSubsteppingBehaviour(b));
    b = behaviour;
    b.esvs.push_back(b.esvs.front());
    EXPECT_ANY_THROW(checkSubsteppingBehaviour(b));
    b = behaviour;
    b.esvs.front().name = "UnknownDriver";
    EXPECT_ANY_THROW(checkSubsteppingBehaviour(b));
    b = behaviour;
    b.to_blocks[0].first.type = mgis::behaviour::Variable::SCALAR;
    EXPECT_ANY_THROW(checkSubsteppingBehaviour(b));
}

TEST(MFrontSubsteppingConfig, InvalidLimitsAndPerturbations)
{
    for (double const h : {0., -1e-7, std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::quiet_NaN()})
    {
        SubsteppingSettings s;
        s.tangent_perturbation = h;
        EXPECT_ANY_THROW(validateSubsteppingSettings(s));
    }
    SubsteppingSettings s;
    s.initial_substeps = 0;
    EXPECT_ANY_THROW(validateSubsteppingSettings(s));
    s.initial_substeps = 4;
    s.max_substeps = 2;
    EXPECT_ANY_THROW(validateSubsteppingSettings(s));
    for (auto const xml :
         {"<substepping><initial_substeps>-1</initial_substeps></substepping>",
          "<substepping><max_substeps>0</max_substeps></substepping>",
          "<substepping><growth_factor>1</growth_factor></substepping>"})
    {
        auto tree = Tests::readXml(xml);
        BaseLib::ConfigTree config(std::move(tree), "",
                                   BaseLib::ConfigTree::onerror,
                                   BaseLib::ConfigTree::onwarning);
        EXPECT_ANY_THROW(createSubsteppingSettings(config));
    }
    auto tree = Tests::readXml("<unrelated>1</unrelated>");
    BaseLib::ConfigTree config(std::move(tree), "",
                               BaseLib::ConfigTree::onerror,
                               BaseLib::ConfigTree::onwarning);
    config.ignoreConfigParameter("unrelated");
    EXPECT_FALSE(createSubsteppingSettings(config));
}

TEST_F(MFrontSubstepping, RejectMalformedDataAndUnrepresentablePerturbation)
{
    auto bad = data;
    bad.s0.gradients.pop_back();
    EXPECT_ANY_THROW(integrateWithSubstepping(
        bad, settings, [](Data& d) { return update(d); }));
    bad = data;
    bad.dt = -1.;
    EXPECT_ANY_THROW(integrateWithSubstepping(
        bad, settings, [](Data& d) { return update(d); }));
    settings.initial_substeps = 4;
    settings.tangent_perturbation = 1e-30;
    EXPECT_ANY_THROW(integrateWithSubstepping(
        data, settings, [](Data& d) { return update(d); }));
}
TEST_F(MFrontSubstepping, FailureAfterSuccessfulSubstepRestartsOriginalHistory)
{
    settings.initial_substeps = 2;
    auto const count = integrateWithSubstepping(
        data, settings,
        [](Data& d)
        {
            if (d.dt == .5 && d.s0.gradients[0] > 0.)
            {
                d.s0.internal_state_variables[0] = -100.;
                d.s1.internal_state_variables[0] = -100.;
                return -1;
            }
            return update(d, 1.);
        });
    ASSERT_EQ(count, 4u);
    double q = 0.;
    for (unsigned i = 0; i < 4; ++i)
        q += .2 * (1. + q) * (1. + q);
    EXPECT_NEAR(data.s1.internal_state_variables[0], q, 1e-13);
    EXPECT_DOUBLE_EQ(data.s0.internal_state_variables[0], 0.);
}

TEST_F(MFrontSubstepping, OneSidedWholePathDerivativeIsOptIn)
{
    auto corner = [](Data& d)
    {
        if (d.s1.gradients[0] > .79 && d.s1.gradients[1] < .2)
            return -1;
        return update(d);
    };
    Data const original = data;
    EXPECT_FALSE(integrateWithSubstepping(data, settings, corner));
    EXPECT_EQ(data.s1.thermodynamic_forces, original.s1.thermodynamic_forces);
    settings.allow_one_sided_tangent = true;
    ASSERT_EQ(integrateWithSubstepping(data, settings, corner), 4u);
    EXPECT_NEAR(data.K[5], 3., 1e-7);
    EXPECT_DOUBLE_EQ(data.s0.internal_state_variables[0], 0.);
}

TEST_F(MFrontSubstepping, LowerWholePathDerivativeHasCorrectSign)
{
    settings.allow_one_sided_tangent = true;
    auto corner = [](Data& d)
    {
        if (d.s1.gradients[0] > .79 && d.s1.gradients[1] > .2)
            return -1;
        return update(d);
    };
    ASSERT_EQ(integrateWithSubstepping(data, settings, corner), 4u);
    EXPECT_NEAR(data.K[5], 3., 1e-7);
}

TEST_F(MFrontSubstepping, TwoFailedWholePathProbesStillReject)
{
    settings.allow_one_sided_tangent = true;
    Data const original = data;
    auto corner = [](Data& d)
    {
        if (d.s1.gradients[0] > .79 && d.s1.gradients[1] != .2)
        {
            d.s0.internal_state_variables[0] = -100.;
            d.s1.thermodynamic_forces[0] = -100.;
            return -1;
        }
        return update(d);
    };
    EXPECT_FALSE(integrateWithSubstepping(data, settings, corner));
    EXPECT_EQ(data.s0.internal_state_variables,
              original.s0.internal_state_variables);
    EXPECT_EQ(data.s1.thermodynamic_forces, original.s1.thermodynamic_forces);
}

TEST(MFrontMCGuard, CapAndInadmissibleStartAreRejectedBeforeTrials)
{
    auto b = mgis::behaviour::load("libOgsMFrontBehaviour", "MohrCoulombAbboSloan",
                                 mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data data{b};
    data.dt = 1.;
    data.s0.material_properties = data.s1.material_properties = {1e7,.3,100,30,0,29,101};
    data.s0.external_state_variables = data.s1.external_state_variables = {293.15};
    SubsteppingSettings settings;
    int calls = 0;
    auto callback = [&](Data&) { ++calls; return -1; };
    EXPECT_ANY_THROW(integrateWithSubstepping(data, settings, callback));
    EXPECT_EQ(calls, 0);
    data.s0.material_properties.back() = data.s1.material_properties.back() = 10;
    data.s0.thermodynamic_forces = {-1050,-10500,-1050,0};
    data.s1.thermodynamic_forces = data.s0.thermodynamic_forces;
    auto const saved = data;
    EXPECT_FALSE(integrateWithSubstepping(data,settings,callback));
    EXPECT_EQ(calls,0);
    EXPECT_EQ(data.s0.thermodynamic_forces,saved.s0.thermodynamic_forces);
    EXPECT_EQ(data.s1.internal_state_variables,saved.s1.internal_state_variables);
}

TEST(MFrontMCGuard, AdmissiblePrestressConditionsTheRealMaterial)
{
    auto b = mgis::behaviour::load("libOgsMFrontBehaviour", "MohrCoulombAbboSloan",
                                 mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data data{b}; data.dt=1.;
    data.s0.material_properties = data.s1.material_properties = {1e7,.3,100,30,0,29,10};
    data.s0.external_state_variables = data.s1.external_state_variables = {293.15};
    data.s0.thermodynamic_forces = data.s1.thermodynamic_forces = {-6300,-10500,-6300,0};
    auto const original=data;
    auto integrate=[&](Data& d){auto view=mgis::behaviour::make_view(d);return mgis::behaviour::integrate(view,b);};
    SubsteppingSettings settings;
    ASSERT_EQ(integrateWithSubstepping(data,settings,integrate),1u);
    EXPECT_EQ(data.s0.internal_state_variables,original.s0.internal_state_variables);
    EXPECT_NEAR(data.s1.internal_state_variables.back(),0.,1e-12);
    for(int i=0;i<3;++i)
    {
        double const expected=(1.3*original.s0.thermodynamic_forces[i]-.3*(-23100))/1e7;
        EXPECT_NEAR(data.s1.internal_state_variables[i],expected,1e-12);
        EXPECT_NEAR(data.s1.thermodynamic_forces[i],original.s0.thermodynamic_forces[i],1e-6);
    }
}

TEST(MFrontMCRecovery, RealPlasticPointAndWholePathTangentSensitivity)
{
    auto b = mgis::behaviour::load("libOgsMFrontBehaviour", "MohrCoulombAbboSloan",
                                 mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data start{b};start.dt=.1;
    start.s0.material_properties=start.s1.material_properties={1e7,.3,100,20,0,29,10};
    start.s0.external_state_variables=start.s1.external_state_variables={293.15};
    // Captured accepted HM element 0/IP 0 target that required eight local substeps.
    start.s1.gradients={-1.1113657405904537e-08,-0.0007528362520822507,0.,-5.3244349381400745e-09};
    auto integrate=[&](Data& d){auto view=mgis::behaviour::make_view(d);return mgis::behaviour::integrate(view,b);};
    Data one=start;one.K[0]=4.;
    EXPECT_NE(integrate(one),1);
    Data result=start;
    SubsteppingSettings settings;
    auto selected=integrateWithSubstepping(result,settings,integrate);
    ASSERT_TRUE(selected);EXPECT_GT(*selected,1u);
    EXPECT_GT(result.s1.internal_state_variables.back(),1e-7);
    EXPECT_EQ(result.s0.internal_state_variables,start.s0.internal_state_variables);
    double const trace_elastic=result.s1.internal_state_variables[0]+result.s1.internal_state_variables[1]+result.s1.internal_state_variables[2];
    EXPECT_NEAR(trace_elastic,start.s1.gradients[0]+start.s1.gradients[1]+start.s1.gradients[2],1e-10);
    settings.initial_substeps=settings.max_substeps=*selected;
    settings.tangent_perturbation*=.5;
    Data probe=start;ASSERT_TRUE(integrateWithSubstepping(probe,settings,integrate));
    double scale=1.;
    for(auto const value : result.K) { scale=std::max(scale,std::abs(value)); }
    double maximum_difference=0.;
    for(std::size_t i=0;i<result.K.size();++i)
    {
        maximum_difference=std::max(maximum_difference,std::abs(result.K[i]-probe.K[i]));
        // FD cancellation at nominally zero cross derivatives needs an absolute
        // floor relative to the stiffness scale, not a 1 Pa/strain component.
        EXPECT_NEAR(result.K[i],probe.K[i],1e-8*scale+1e-3*std::abs(result.K[i]));
    }
    EXPECT_LT(maximum_difference/scale,1e-6);
    RecordProperty("whole_path_tangent_relative_difference",maximum_difference/scale);
}

TEST_F(MFrontSubstepping, HostCommitPreservesDensity)
{
    MaterialStateVariablesMFront<2> state{0, behaviour};
    state._behaviour_data.s0.mass_density = 100.;
    state._behaviour_data.s1.mass_density = 200.;
    state._behaviour_data.s1.gradients[0] = .01;
    state.pushBackState();
    EXPECT_DOUBLE_EQ(state._behaviour_data.s0.mass_density, 200.);
    EXPECT_DOUBLE_EQ(state._behaviour_data.s0.gradients[0], .01);
}

} // namespace MFrontSubsteppingTests

int main(int argc, char** argv)
{
 BaseLib::initOGSLogger("info");
 ::testing::InitGoogleTest(&argc, argv);
 return RUN_ALL_TESTS();
}

// Reference capture is tested independently of mesh assembly using a prescribed
// increment boundary. The production factory supplies a stock Dirichlet BC.
#include "ProcessLib/BoundaryConditionAndSourceTerm/RelativeDirichletBoundaryCondition.h"

namespace
{
class TestIncrementBC final : public ProcessLib::BoundaryCondition
{
public:
    double at_reference = 0.;
    bool reverse_ids = false;
    void getEssentialBCValues(
        double const t, GlobalVector const&,
        NumLib::IndexValueVector<GlobalIndexType>& values) const override
    {
        values.ids = reverse_ids ? std::vector<GlobalIndexType>{1, 0}
                                 : std::vector<GlobalIndexType>{0, 1};
        values.values = {at_reference - .01 * std::max(0., t - 1.),
                         at_reference - .01 * std::max(0., t - 1.)};
    }
};

TEST(HMRelativeDisplacement, InactiveUntilAcceptedReference)
{
    auto increment = std::make_unique<TestIncrementBC>();
    ProcessLib::RelativeDirichletBoundaryCondition bc(std::move(increment), 1.);
    GlobalVector x(2);
    x[0] = -.002;
    x[1] = -.003;
    NumLib::IndexValueVector<GlobalIndexType> values;
    bc.getEssentialBCValues(.5, x, values);
    EXPECT_TRUE(values.ids.empty());
    bc.getEssentialBCValues(1., x, values);
    EXPECT_TRUE(values.ids.empty());
    // A converged trial is not global acceptance.
    bc.postTimestep(1., {&x}, 0);
    EXPECT_THROW(bc.getEssentialBCValues(1.1, x, values), std::runtime_error);
    bc.postAcceptedTimestep(1., {&x}, 0);
    bc.getEssentialBCValues(1.1, x, values);
    EXPECT_DOUBLE_EQ(values.values[0], -.003000000000000001);
    EXPECT_NEAR(values.values[1], -.004, 1e-15);
}

TEST(HMRelativeDisplacement, RejectedTrialAndLaterAcceptanceCannotRecapture)
{
    ProcessLib::RelativeDirichletBoundaryCondition bc(
        std::make_unique<TestIncrementBC>(), 1.);
    GlobalVector x(2);
    x[0] = -.002;
    x[1] = -.003;
    bc.postAcceptedTimestep(1., {&x}, 0);
    x[0] = 999.;
    x[1] = 888.;
    bc.postTimestep(1., {&x}, 0);
    bc.postAcceptedTimestep(2., {&x}, 0);
    NumLib::IndexValueVector<GlobalIndexType> values;
    bc.getEssentialBCValues(1.2, x, values);
    EXPECT_NEAR(values.values[0], -.004, 1e-15);
    EXPECT_NEAR(values.values[1], -.005, 1e-15);
    EXPECT_DOUBLE_EQ(x[0], 999.);
}

TEST(HMRelativeDisplacement, SkippedReferenceIsRejected)
{
    ProcessLib::RelativeDirichletBoundaryCondition bc(
        std::make_unique<TestIncrementBC>(), 1.);
    GlobalVector x(2);
    EXPECT_THROW(bc.postAcceptedTimestep(1.1, {&x}, 0), std::runtime_error);
}

TEST(HMRelativeDisplacement, NonzeroInitialIncrementIsRejected)
{
    auto increment = std::make_unique<TestIncrementBC>();
    increment->at_reference = -.1;
    ProcessLib::RelativeDirichletBoundaryCondition bc(std::move(increment), 1.);
    GlobalVector x(2);
    EXPECT_THROW(bc.postAcceptedTimestep(1., {&x}, 0), std::runtime_error);
}

TEST(HMRelativeDisplacement, ChangedDofMapIsRejected)
{
    auto increment = std::make_unique<TestIncrementBC>();
    auto* editable = increment.get();
    ProcessLib::RelativeDirichletBoundaryCondition bc(std::move(increment), 1.);
    GlobalVector x(2);
    x.setZero();
    bc.postAcceptedTimestep(1., {&x}, 0);
    editable->reverse_ids = true;
    NumLib::IndexValueVector<GlobalIndexType> values;
    EXPECT_THROW(bc.getEssentialBCValues(1.1, x, values), std::runtime_error);
}

TEST(HMRelativeDisplacement, InvalidReferenceTimeIsRejected)
{
    EXPECT_THROW(ProcessLib::RelativeDirichletBoundaryCondition(
                     std::make_unique<TestIncrementBC>(), 0.),
                 std::runtime_error);
}
}  // namespace

#include "ProcessLib/BoundaryConditionAndSourceTerm/PhaseDirichletBoundaryCondition.h"
namespace
{
TEST(HMPhaseHydraulics, OpenLeftClosedRightOwnership)
{
    ProcessLib::PhaseDirichletBoundaryCondition bc(std::make_unique<TestIncrementBC>(), 1., 2.);
    GlobalVector x(2);
    x.setZero();
    NumLib::IndexValueVector<GlobalIndexType> values;
    bc.getEssentialBCValues(1., x, values);
    EXPECT_TRUE(values.ids.empty());
    bc.getEssentialBCValues(1.1, x, values);
    EXPECT_EQ(values.ids.size(), 2);
    bc.getEssentialBCValues(2., x, values);
    EXPECT_EQ(values.ids.size(), 2);
    bc.getEssentialBCValues(2.1, x, values);
    EXPECT_TRUE(values.ids.empty());
}
TEST(HMPhaseHydraulics, ConsecutivePhasesCannotBothOwnSharedEndpoint)
{
    ProcessLib::PhaseDirichletBoundaryCondition first(std::make_unique<TestIncrementBC>(), 1., 2.);
    ProcessLib::PhaseDirichletBoundaryCondition second(std::make_unique<TestIncrementBC>(), 2., 3.);
    GlobalVector x(2);
    x.setZero();
    NumLib::IndexValueVector<GlobalIndexType> a, b;
    first.getEssentialBCValues(2., x, a);
    second.getEssentialBCValues(2., x, b);
    EXPECT_FALSE(a.ids.empty());
    EXPECT_TRUE(b.ids.empty());
}
TEST(HMPhaseHydraulics, InvalidIntervalsRejected)
{
    EXPECT_THROW(ProcessLib::PhaseDirichletBoundaryCondition(std::make_unique<TestIncrementBC>(), 2., 1.), std::runtime_error);
    EXPECT_THROW(ProcessLib::PhaseDirichletBoundaryCondition(std::make_unique<TestIncrementBC>(), 1., 1.), std::runtime_error);
}
}  // namespace

#include "ProcessLib/HydroMechanics/MechanicalDrainage.h"
namespace
{
using namespace ProcessLib::HydroMechanics;
TEST(HMMechanicalDrainage, UndrainedConservesLocalContentAndCompressionGeneratesPressure)
{
    MechanicalDrainageState state;
    state.accept(1., 1., -.001, 10000.);
    MechanicalDrainageMaterial material{MechanicalDrainageType::UndrainedEffectiveStress, 8e7};
    auto const r = state.evaluate(10000., -.0013, .8, material);
    EXPECT_NEAR(r.excess, 19200., 1e-10);
    EXPECT_NEAR(r.pressure, 29200., 1e-10);
    EXPECT_NEAR(r.fluid_content, 0., 1e-18);
    EXPECT_DOUBLE_EQ(r.volumetric_stiffness, .8*.8*8e7);
}
TEST(HMMechanicalDrainage, DrainedRetainsNonzeroReferenceAndAccountsForWater)
{
    MechanicalDrainageState state;
    state.accept(1., 1., -.001, 12345.);
    auto const r = state.evaluate(12345., -.0013, .8,
        {MechanicalDrainageType::DrainedEffectiveStress, 8e7});
    EXPECT_DOUBLE_EQ(r.pressure, 12345.);
    EXPECT_DOUBLE_EQ(r.excess, 0.);
    EXPECT_NEAR(r.fluid_content, -.00024, 1e-17);
    EXPECT_DOUBLE_EQ(r.volumetric_stiffness, 0.);
}
TEST(HMMechanicalDrainage, ReferenceOnlyAtAcceptedEndpoint)
{
    MechanicalDrainageState state;
    MechanicalDrainageMaterial material{MechanicalDrainageType::UndrainedEffectiveStress, 8e7};
    state.accept(.9, 1., -.001, 10000.);
    EXPECT_THROW(state.evaluate(10000., -.001, .8, material), std::runtime_error);
    EXPECT_THROW(state.accept(1.1, 1., -.001, 10000.), std::runtime_error);
    state.accept(1., 1., -.001, 10000.);
    state.accept(2., 1., 999., 888.);
    EXPECT_DOUBLE_EQ(state.referenceStrain(), -.001);
    EXPECT_DOUBLE_EQ(state.referencePressure(), 10000.);
}
TEST(HMMechanicalDrainage, RejectedAndRepeatedTrialsDoNotChangeReference)
{
    MechanicalDrainageState state;
    state.accept(1., 1., -.001, 10000.);
    MechanicalDrainageMaterial material{MechanicalDrainageType::UndrainedEffectiveStress, 8e7};
    auto const initial = state.evaluate(10000., -.001, .8, material);
    state.evaluate(10000., -.01, .8, material); // discarded trial
    EXPECT_THROW(state.evaluate(10001., -.01, .8, material), std::runtime_error);
    auto const after = state.evaluate(10000., -.001, .8, material);
    EXPECT_DOUBLE_EQ(after.pressure, initial.pressure);
    EXPECT_DOUBLE_EQ(after.fluid_content, initial.fluid_content);
}
TEST(HMMechanicalDrainage, EffectiveToTotalStressVolumetricTangentMatchesFiniteDifference)
{
    MechanicalDrainageState state;
    state.accept(1., 1., -.001, 10000.);
    MechanicalDrainageMaterial material{MechanicalDrainageType::UndrainedEffectiveStress, 8e7};
    double const alpha = .8, e = -.0013, h = 1e-9;
    auto const upper = state.evaluate(10000., e+h, alpha, material);
    auto const lower = state.evaluate(10000., e-h, alpha, material);
    auto const r = state.evaluate(10000., e, alpha, material);
    double const derivative = -alpha*(upper.pressure-lower.pressure)/(2*h);
    EXPECT_NEAR(derivative, r.volumetric_stiffness, 1e-2);
}
TEST(HMMechanicalDrainage, UnloadingAndHoldUseOriginalReference)
{
    MechanicalDrainageState state;
    state.accept(1., 1., -.001, 10000.);
    MechanicalDrainageMaterial material{MechanicalDrainageType::UndrainedEffectiveStress, 8e7};
    auto const first = state.evaluate(10000., -.0013, .8, material);
    auto const partial = state.evaluate(10000., -.0012, .8, material);
    EXPECT_LT(partial.pressure, first.pressure);
    EXPECT_NEAR(partial.excess, 12800., 1e-10);
    auto const hold = state.evaluate(10000., -.0012, .8, material);
    EXPECT_DOUBLE_EQ(hold.pressure, partial.pressure);
}
TEST(HMMechanicalDrainage, InvalidStorageCouplingAndTrialStatesRejected)
{
    MechanicalDrainageState state;
    state.accept(1., 1., -.001, 10000.);
    for (double const m : {0., -1., std::numeric_limits<double>::infinity()})
    {
        EXPECT_THROW(state.evaluate(10000., -.001, .8,
            {MechanicalDrainageType::UndrainedEffectiveStress, m}), std::runtime_error);
    }
    EXPECT_THROW(state.evaluate(10000., -.001, 0.,
        {MechanicalDrainageType::UndrainedEffectiveStress, 8e7}), std::runtime_error);
    EXPECT_THROW(state.evaluate(10000., std::numeric_limits<double>::quiet_NaN(), .8,
        {MechanicalDrainageType::UndrainedEffectiveStress, 8e7}), std::runtime_error);
}
TEST(HMMechanicalDrainage, StrictStageStartsAfterReferenceNotAtIt)
{
    MechanicalDrainageSettings settings{1., {}};
    EXPECT_FALSE(settings.active(1.));
    EXPECT_TRUE(settings.active(1.1));
}
} // namespace

namespace
{
TEST(HMMechanicalDrainage, RealPlasticMCCombinedTotalStressTangent)
{
    using namespace MFrontSubsteppingTests;
    auto behaviour = mgis::behaviour::load("libOgsMFrontBehaviour", "MohrCoulombAbboSloan",
        mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data initial{behaviour};
    initial.dt = .1;
    initial.s0.material_properties = initial.s1.material_properties = {1e7,.3,100,20,0,29,10};
    initial.s0.external_state_variables = initial.s1.external_state_variables = {293.15};
    initial.s1.gradients = {-1.1113657405904537e-08,-0.0007528362520822507,0.,-5.3244349381400745e-09};
    auto integrate = [&](Data& data) { auto view=mgis::behaviour::make_view(data); return mgis::behaviour::integrate(view,behaviour); };
    SubsteppingSettings settings;
    Data result=initial;
    ASSERT_TRUE(integrateWithSubstepping(result,settings,integrate));
    ASSERT_GT(result.s1.internal_state_variables.back(),1e-7);
    MechanicalDrainageState pressure_state;
    pressure_state.accept(1.,1.,0.,10000.);
    MechanicalDrainageMaterial material{MechanicalDrainageType::UndrainedEffectiveStress,8e7};
    double const alpha=.8,h=1e-9;
    auto total_stress=[&](Data data)
    {
        if(!integrateWithSubstepping(data,settings,integrate)) { throw std::runtime_error("MC probe failed"); }
        double const volume=data.s1.gradients[0]+data.s1.gradients[1]+data.s1.gradients[2];
        auto const p=pressure_state.evaluate(10000.,volume,alpha,material).pressure;
        for(int i=0;i<3;++i) { data.s1.thermodynamic_forces[i]-=alpha*p; }
        return data.s1.thermodynamic_forces;
    };
    for(int column=0;column<4;++column)
    {
        Data up=initial,down=initial;
        up.s1.gradients[column]+=h;down.s1.gradients[column]-=h;
        auto const upper=total_stress(up),lower=total_stress(down);
        for(int row=0;row<4;++row)
        {
            double const expected=result.K[row*4+column]+((row<3 && column<3)?alpha*alpha*material.biot_modulus:0.);
            double const measured=(upper[row]-lower[row])/(2*h);
            EXPECT_NEAR(measured,expected,2e-4*std::abs(expected)+10.);
        }
    }
}
} // namespace


namespace
{
// Constitutive prerequisite only: this does not enable a TSA HM process mode.
TEST(MFrontTotalStress, PureShearStrengthHistoryAndConfinementIndependence)
{
    using namespace MFrontSubsteppingTests;
    auto behaviour = mgis::behaviour::load(
        "libOgsMFrontBehaviour", "MohrCoulombAbboSloan",
        mgis::behaviour::Hypothesis::PLANESTRAIN);
    constexpr double young = 1e7, poisson = .49, strength = 1000.;
    constexpr double shear = young / (2. * (1. + poisson));
    auto integrate = [&](Data& d)
    {
        auto view = mgis::behaviour::make_view(d);
        return mgis::behaviour::integrate(view, behaviour);
    };
    std::ofstream output;
    if (auto const path = std::getenv("OGS_TSA_MATERIAL_CSV"))
    {
        output.open(path);
        ASSERT_TRUE(output.good());
        output << std::setprecision(17)
               << "confinement,step,gamma,tau,mean_stress,plastic_multiplier,plastic_volume,substeps\n";
    }
    std::vector<double> reference;
    double max_error = 0., max_volume = 0.;
    for (double const confinement : {0., 10000., 100000.})
    {
        Data data{behaviour};
        data.dt = .1;
        data.s0.material_properties = data.s1.material_properties =
            {young, poisson, strength, 0., 0., 29., 10.};
        data.s0.external_state_variables = data.s1.external_state_variables = {293.15};
        data.s0.thermodynamic_forces = data.s1.thermodynamic_forces =
            {-confinement, -confinement, -confinement, 0.};
        SubsteppingSettings settings;
        double previous_gamma = 0., expected_tau = 0., previous_lambda = 0.;
        // Load beyond yield, unload through zero, reverse to negative yield,
        // then hold. Engineering shear gamma = sqrt(2)*Kelvin strain.
        for (int step = 0; step <= 120; ++step)
        {
            double const gamma = step <= 40 ? step * .00005
                : (step <= 100 ? .002 - (step - 40) * .00005 : -.001);
            data.s1 = data.s0;
            data.s1.gradients[3] = gamma / std::sqrt(2.);
            auto const count = integrateWithSubstepping(data, settings, integrate);
            ASSERT_TRUE(count) << "step=" << step << " confinement=" << confinement;
            expected_tau = std::clamp(expected_tau + shear * (gamma - previous_gamma),
                                      -strength, strength);
            double const tau = data.s1.thermodynamic_forces[3] / std::sqrt(2.);
            double const mean = (data.s1.thermodynamic_forces[0] +
                data.s1.thermodynamic_forces[1] + data.s1.thermodynamic_forces[2]) / 3.;
            auto const& elastic = data.s1.internal_state_variables;
            double const elastic_volume = elastic[0] + elastic[1] + elastic[2];
            double const plastic_volume = -elastic_volume -
                3. * confinement * (1. - 2. * poisson) / young;
            max_error = std::max(max_error, std::abs(tau - expected_tau));
            max_volume = std::max(max_volume, std::abs(plastic_volume));
            EXPECT_NEAR(tau, expected_tau, 1e-5);
            EXPECT_NEAR(mean, -confinement, 1e-5);
            EXPECT_NEAR(plastic_volume, 0., 1e-12);
            EXPECT_GE(elastic.back() + 1e-12, previous_lambda);
            if (confinement == 0.)
                reference.push_back(tau);
            else
                EXPECT_NEAR(tau, reference[step], 1e-5);
            if (output)
                output << confinement << ',' << step << ',' << gamma << ',' << tau
                       << ',' << mean << ',' << elastic.back() << ','
                       << plastic_volume << ',' << *count << '\n';
            previous_lambda = elastic.back();
            previous_gamma = gamma;
            data.s0 = data.s1;
        }
        EXPECT_GT(previous_lambda, 1e-4);
    }
    RecordProperty("maximum_shear_error_Pa", max_error);
    RecordProperty("maximum_plastic_volume", max_volume);
}
}  // namespace

namespace
{
TEST(MFrontTotalStress, CutoffIsInactiveAtZeroFrictionAndEffectiveStressControlDiffers)
{
    using namespace MFrontSubsteppingTests;
    auto behaviour = mgis::behaviour::load(
        "libOgsMFrontBehaviour", "MohrCoulombAbboSloan",
        mgis::behaviour::Hypothesis::PLANESTRAIN);
    auto endpoint = [&](double phi, double cutoff, double confinement)
    {
        Data data{behaviour};
        data.dt = .1;
        data.s0.material_properties = data.s1.material_properties =
            {1e7, .49, 1000., phi, 0., 29., cutoff};
        data.s0.external_state_variables = data.s1.external_state_variables = {293.15};
        data.s0.thermodynamic_forces = data.s1.thermodynamic_forces =
            {-confinement, -confinement, -confinement, 0.};
        auto integrate = [&](Data& d)
        {
            auto view = mgis::behaviour::make_view(d);
            return mgis::behaviour::integrate(view, behaviour);
        };
        SubsteppingSettings settings;
        for (int i = 1; i <= 40; ++i)
        {
            data.s1 = data.s0;
            data.s1.gradients[3] = i * .00005 / std::sqrt(2.);
            if (!integrateWithSubstepping(data, settings, integrate))
                throw std::runtime_error("TSA material control failed");
            data.s0 = data.s1;
        }
        return data.s1.thermodynamic_forces[3] / std::sqrt(2.);
    };
    for (double const a : {0., 10., 100.})
        EXPECT_NEAR(endpoint(0., a, 10000.), 1000., 1e-5);
    double const unconfined = endpoint(20., 10., 0.);
    double const confined = endpoint(20., 10., 10000.);
    EXPECT_GT(confined - unconfined, 3000.);
    RecordProperty("effective_stress_control_confinement_gain_Pa", confined-unconfined);
}
}  // namespace

namespace
{
TEST(HMMechanicalStorage, DerivedWaterStorageAndCompressibleGrains)
{
    MechanicalDrainageMaterial const water{MechanicalDrainageType::UndrainedEffectiveStress, 0., 2.2e9};
    EXPECT_DOUBLE_EQ(water.resolved(1., .3, 1e7).biot_modulus, 2.2e9/.3);
    double const kd = 1e7/(3.*(1.-2.*.3));
    double const expected = 1./(.3/2.2e9+(.8-.3)*(1.-.8)/kd);
    EXPECT_DOUBLE_EQ(water.resolved(.8,.3,kd).biot_modulus, expected);
    MechanicalDrainageMaterial const explicit_value{MechanicalDrainageType::UndrainedEffectiveStress, 8e7};
    EXPECT_DOUBLE_EQ(explicit_value.resolved(.8,.3,kd).biot_modulus,8e7);
}
TEST(HMMechanicalStorage, DrainedResponseNeedsNoStorage)
{
    MechanicalDrainageMaterial const material{MechanicalDrainageType::DrainedEffectiveStress, 0.};
    auto const resolved = material.resolved(1.,.3,1e7);
    MechanicalDrainageState state;
    state.accept(1.,1.,-.001,10000.);
    auto const result = state.evaluate(10000.,-.0013,1.,resolved);
    EXPECT_DOUBLE_EQ(result.pressure,10000.);
    EXPECT_NEAR(result.fluid_content,-.0003,1e-17);
    EXPECT_DOUBLE_EQ(result.volumetric_stiffness,0.);
}
TEST(HMMechanicalStorage, InvalidDerivedInputsAreRejected)
{
    MechanicalDrainageMaterial const water{MechanicalDrainageType::UndrainedEffectiveStress,0.,2.2e9};
    for(double const n : {0.,-1.,1.,std::numeric_limits<double>::quiet_NaN()})
        EXPECT_THROW(water.resolved(1.,n,1e7),std::runtime_error);
    EXPECT_THROW(water.resolved(.2,.3,1e7),std::runtime_error);
    EXPECT_THROW(water.resolved(1.,.3,0.),std::runtime_error);
    MechanicalDrainageMaterial const missing{MechanicalDrainageType::UndrainedEffectiveStress,0.};
    EXPECT_THROW(missing.resolved(1.,.3,1e7),std::runtime_error);
}
}  // namespace


namespace
{
TEST(HMTotalStress, CannotAccidentallyUseEffectiveStressStorage)
{
    MechanicalDrainageMaterial const material{MechanicalDrainageType::UndrainedTotalStress,0.};
    EXPECT_TRUE(material.usesTotalStress());
    EXPECT_THROW(material.resolved(1.,.3,1e7),std::runtime_error);
    MechanicalDrainageState state;
    state.accept(1.,1.,-.001,10000.);
    EXPECT_THROW(state.evaluate(10000.,-.0011,1.,material),std::runtime_error);
}
}  // namespace

#include "ProcessLib/HydroMechanics/InitialisationGate.h"
TEST(HMInitialisationGate, IndependentAbboSloanInvariantReferences)
{
    using namespace ProcessLib::HydroMechanics;
    EXPECT_NEAR(mcAdmissibilityResidual(Eigen::Vector4d{0,0,0,0}, MCAdmissibilityParameters{100,20,29,10}), -90.549060645334151, 1e-7);
    EXPECT_NEAR(mcAdmissibilityResidual(Eigen::Vector4d{-10000,-20000,-14000,565.68542494923804}, MCAdmissibilityParameters{100,20,29,10}), -208.29573719834477, 1e-7);
    EXPECT_NEAR(mcAdmissibilityResidual(Eigen::Vector4d{-12000,-4000,-4000,0}, MCAdmissibilityParameters{2000,30,29,10}), -1709.3776597370768, 1e-7);
    EXPECT_NEAR(mcAdmissibilityResidual(Eigen::Vector4d{-30000,-28000,-28000,989.94949366116657}, MCAdmissibilityParameters{1000,0,29,10}), 220.65556157336755, 1e-7);
    EXPECT_NEAR(mcAdmissibilityResidual(Eigen::Vector4d{100,200,0,424.26406871192853}, MCAdmissibilityParameters{100,20,29,10}), 261.49009269320766, 1e-7);
    EXPECT_NEAR(mcAdmissibilityResidual(Eigen::Vector4d{-100,100,0,0}, MCAdmissibilityParameters{100,20,20,10}), 6.0892097158756684, 1e-7);
    EXPECT_NEAR(mcAdmissibilityResidual(Eigen::Vector4d{-100,-100,200,0}, MCAdmissibilityParameters{100,20,29,10}), 73.540829285660124, 1e-7);
    EXPECT_THROW(mcAdmissibilityResidual(Eigen::Vector4d::Constant(std::numeric_limits<double>::quiet_NaN()), MCAdmissibilityParameters{100,20,29,10}), std::runtime_error);
}

#include "MaterialLib/SolidModels/MFront/CheckpointData.h"
#include <array>
#include <iterator>
namespace MFrontSubsteppingTests
{
TEST_F(MFrontSubstepping, CheckpointRoundTripPreservesAllMGISDataAndReceiverBuffer)
{
    // The substepping fixture replaces the real elastic ISV layout with a
    // synthetic history variable. Use the actual behaviour schema for restart.
    Data actual{behaviour};data=actual;
    data.s0.thermodynamic_forces={-1.,-2.,-3.,.4};
    data.s1.thermodynamic_forces={-4.,-5.,-6.,.7};
    data.s0.material_properties={1e7,.3};data.s1.material_properties={2e7,.25};
    data.s0.external_state_variables={293.15};data.s1.external_state_variables={294.15};
    data.dt=.125;data.rdt=.75;data.speed_of_sound=1234.;
    data.s0.mass_density=987.;data.s1.mass_density=1234.;
    data.s0.stored_energy=1.;data.s1.stored_energy=2.;
    data.s0.dissipated_energy=3.;data.s1.dissipated_energy=4.;
    data.K.assign(data.K.size(), 9876.);
    data.s0.gradients[0]=.0001;data.s1.gradients[0]=.002;
    auto const saved=captureMFrontCheckpoint(data);
    Data receiver{behaviour};std::array<char,512> error{};receiver.error_message=error.data();
    restoreMFrontCheckpoint(receiver,saved);
    EXPECT_EQ(captureMFrontCheckpoint(receiver),saved);
    EXPECT_EQ(receiver.error_message,error.data());
    EXPECT_DOUBLE_EQ(receiver.s0.mass_density,987.);EXPECT_DOUBLE_EQ(receiver.s1.mass_density,1234.);
    receiver.s1.gradients[0]=999.;EXPECT_EQ(captureMFrontCheckpoint(data),saved);
}
TEST_F(MFrontSubstepping, CheckpointCorruptionAndWrongBehaviourLeaveReceiverUnchanged)
{
    auto const saved=captureMFrontCheckpoint(data);Data receiver{behaviour};
    auto const original=captureMFrontCheckpoint(receiver);
    EXPECT_THROW(restoreMFrontCheckpoint(receiver,saved),std::runtime_error);
    EXPECT_EQ(captureMFrontCheckpoint(receiver),original);
    for (auto corrupt : {saved.substr(0,10),saved.substr(0,saved.size()-1),std::string("broken")})
    { EXPECT_THROW(restoreMFrontCheckpoint(receiver,corrupt),std::runtime_error);EXPECT_EQ(captureMFrontCheckpoint(receiver),original); }
    auto corrupt=saved;corrupt.back()^=1;
    EXPECT_THROW(restoreMFrontCheckpoint(receiver,corrupt),std::runtime_error);
    EXPECT_EQ(captureMFrontCheckpoint(receiver),original);
    auto const mc=mgis::behaviour::load("libOgsMFrontBehaviour","MohrCoulombAbboSloan",mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data incompatible{mc};auto const untouched=captureMFrontCheckpoint(incompatible);
    EXPECT_THROW(restoreMFrontCheckpoint(incompatible,saved),std::runtime_error);
    EXPECT_EQ(captureMFrontCheckpoint(incompatible),untouched);
}
TEST(HMMFrontCheckpoint, PlasticMCBranchesMatchUninterruptedMaterialHistory)
{
    auto const b=mgis::behaviour::load("libOgsMFrontBehaviour","MohrCoulombAbboSloan",mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data state{b};state.dt=1.;
    state.s0.material_properties=state.s1.material_properties={1e7,.49,1000,0,0,29,10};
    state.s0.external_state_variables=state.s1.external_state_variables={293.15};
    state.s0.gradients.assign(4,0.);state.s1.gradients=state.s0.gradients;state.s1.gradients[3]=.002/std::sqrt(2.);
    auto integrate=[&](Data& d){auto view=mgis::behaviour::make_view(d);return mgis::behaviour::integrate(view,b);};
    SubsteppingSettings settings;ASSERT_TRUE(integrateWithSubstepping(state,settings,integrate));
    ASSERT_GT(state.s1.internal_state_variables.back(),1e-5);
    mgis::behaviour::update(state);state.s0.mass_density=state.s1.mass_density;
    auto const checkpoint=captureMFrontCheckpoint(state);
    RecordProperty("checkpoint_bytes",static_cast<int>(checkpoint.size()));
    RecordProperty("checkpoint_plastic_multiplier",std::to_string(state.s1.internal_state_variables.back()));
    RecordProperty("independent_branches",3);
    if (auto const path=std::getenv("HM_CHECKPOINT_TEST_OUTPUT"))
    {
        { std::ofstream out(path,std::ios::binary);ASSERT_TRUE(out.good());out.write(checkpoint.data(),checkpoint.size());ASSERT_TRUE(out.good()); }
        std::ifstream in(path,std::ios::binary);
        std::string const persisted((std::istreambuf_iterator<char>(in)),std::istreambuf_iterator<char>());
        ASSERT_EQ(persisted,checkpoint);
        Data loaded{b};restoreMFrontCheckpoint(loaded,persisted);
        EXPECT_EQ(captureMFrontCheckpoint(loaded),checkpoint);
    }
    int checkpoint_branch_index=0;
    for (double const gamma : {.0015,-.001,.003})
    {
        Data reference=state;Data branch{b};restoreMFrontCheckpoint(branch,checkpoint);
        reference.s1.gradients[3]=branch.s1.gradients[3]=gamma/std::sqrt(2.);
        if (gamma < 0.)
        {
            Data rejected = branch;
            auto const unchanged = captureMFrontCheckpoint(rejected);
            SubsteppingSettings single_step;
            single_step.max_substeps = 1;
            EXPECT_FALSE(integrateWithSubstepping(rejected, single_step, integrate));
            EXPECT_EQ(captureMFrontCheckpoint(rejected), unchanged);
        }
        auto const reference_count = integrateWithSubstepping(reference,settings,integrate);
        auto const branch_count = integrateWithSubstepping(branch,settings,integrate);
        ASSERT_TRUE(reference_count);
        ASSERT_TRUE(branch_count);
        EXPECT_EQ(reference_count, branch_count);
        if (gamma < 0.) { EXPECT_GT(*branch_count, 1u); }
        RecordProperty("branch_substeps_" + std::to_string(checkpoint_branch_index),
                       static_cast<int>(*branch_count));
        EXPECT_EQ(captureMFrontCheckpoint(reference),captureMFrontCheckpoint(branch));
        // Replay equivalence alone cannot establish physical admissibility.
        double const shear_modulus = 1e7 / (2. * (1. + .49));
        double const initial_tau = state.s0.thermodynamic_forces[3] / std::sqrt(2.);
        double const initial_gamma = state.s0.gradients[3] * std::sqrt(2.);
        double const expected_tau = std::clamp(
            initial_tau + shear_modulus * (gamma - initial_gamma), -1000., 1000.);
        double const expected_lambda = state.s0.internal_state_variables.back() +
            std::abs(gamma - initial_gamma - (expected_tau - initial_tau) / shear_modulus);
        EXPECT_NEAR(branch.s1.thermodynamic_forces[3] / std::sqrt(2.), expected_tau, 1e-6);
        EXPECT_NEAR(branch.s1.internal_state_variables.back(), expected_lambda, 1e-10);
        if (auto const path=std::getenv("HM_CHECKPOINT_TEST_OUTPUT"))
        {
            auto const bytes=captureMFrontCheckpoint(branch);
            std::ofstream out(std::string(path)+".branch_"+std::to_string(checkpoint_branch_index),std::ios::binary);
            ASSERT_TRUE(out.good());out.write(bytes.data(),bytes.size());ASSERT_TRUE(out.good());
        }
        ++checkpoint_branch_index;
        EXPECT_EQ(captureMFrontCheckpoint(state),checkpoint);
    }
}
} // namespace MFrontSubsteppingTests

#include "ProcessLib/HydroMechanics/IntegrationPointCheckpoint.h"

namespace MFrontSubsteppingTests
{
namespace HMIP = ProcessLib::HydroMechanics;
struct CheckpointHostIP
{
    using KV = MathLib::KelvinVector::KelvinVectorType<2>;
    KV sigma_eff = KV::Zero(), sigma_eff_prev = KV::Zero();
    KV eps = KV::Zero(), eps_prev = KV::Zero();
    double integration_weight = .125;
    double strain_rate_variable = .375;
    HMIP::MechanicalDrainageState mechanical_drainage;
    double mechanical_pressure = 0., mechanical_excess = 0., mechanical_fluid_content = 0.;
    std::unique_ptr<MaterialLib::Solids::MechanicsBase<2>::MaterialStateVariables> material_state_variables;
    explicit CheckpointHostIP(mgis::behaviour::Behaviour const& b, int offset = 4)
        : material_state_variables(std::make_unique<MaterialStateVariablesMFront<2>>(offset, b)) {}
    auto& material()
    { return static_cast<MaterialStateVariablesMFront<2>&>(*material_state_variables); }
};

TEST(HMIntegrationPointCheckpoint, CompleteHostStatePlasticBranchesAndAtomicRejection)
{
    auto const b = mgis::behaviour::load("libOgsMFrontBehaviour", "MohrCoulombAbboSloan", mgis::behaviour::Hypothesis::PLANESTRAIN);
    CheckpointHostIP source{b};
    auto& data = source.material()._behaviour_data;
    data.dt = 1.;
    data.s0.material_properties = data.s1.material_properties = {1e7,.49,1000,0,0,29,10};
    data.s0.external_state_variables = data.s1.external_state_variables = {293.15};
    data.s1.gradients = {-.0001,-.0001,-.0001,.002/std::sqrt(2.)};
    auto integrate = [&](Data& d) { auto view = mgis::behaviour::make_view(d); return mgis::behaviour::integrate(view,b); };
    SubsteppingSettings settings;
    auto const count = integrateWithSubstepping(data,settings,integrate);
    ASSERT_TRUE(count);
    source.material().last_local_substeps = *count;
    mgis::behaviour::update(data);
    data.s0.mass_density = data.s1.mass_density = 2050.;
    source.sigma_eff = Eigen::Map<CheckpointHostIP::KV const>(data.s1.thermodynamic_forces.data());
    source.eps = Eigen::Map<CheckpointHostIP::KV const>(data.s1.gradients.data());
    // Distinct previous values explicitly test that no field is reconstructed
    // or accidentally replaced by the current state during roundtrip.
    source.sigma_eff_prev = source.sigma_eff * .75;
    source.eps_prev = source.eps * .5;
    source.mechanical_drainage.accept(1.,1.,0.,10000.);
    auto response = source.mechanical_drainage.evaluate(10000.,source.eps.head<3>().sum(),1.,
        {HMIP::MechanicalDrainageType::UndrainedEffectiveStress,1e7});
    source.mechanical_pressure=response.pressure;source.mechanical_excess=response.excess;source.mechanical_fluid_content=response.fluid_content;
    auto const capsule = HMIP::captureHMIntegrationPoint<2>(source,false);
    CheckpointHostIP restored{b};
    auto parsed = HMIP::decodeHMIntegrationPoint<2>(restored,capsule,false);
    HMIP::commitHMIntegrationPoint<2>(restored,std::move(parsed));
    EXPECT_EQ(HMIP::captureHMIntegrationPoint<2>(restored,false),capsule);
    EXPECT_TRUE((restored.sigma_eff.array()==source.sigma_eff.array()).all());
    EXPECT_TRUE((restored.sigma_eff_prev.array()==source.sigma_eff_prev.array()).all());
    EXPECT_TRUE((restored.eps.array()==source.eps.array()).all());
    EXPECT_TRUE((restored.eps_prev.array()==source.eps_prev.array()).all());
    EXPECT_DOUBLE_EQ(restored.strain_rate_variable,source.strain_rate_variable);
    EXPECT_DOUBLE_EQ(restored.mechanical_pressure,source.mechanical_pressure);
    EXPECT_DOUBLE_EQ(restored.mechanical_excess,source.mechanical_excess);
    EXPECT_DOUBLE_EQ(restored.mechanical_fluid_content,source.mechanical_fluid_content);
    EXPECT_DOUBLE_EQ(restored.mechanical_drainage.referenceStrain(),source.mechanical_drainage.referenceStrain());
    EXPECT_DOUBLE_EQ(restored.mechanical_drainage.referencePressure(),source.mechanical_drainage.referencePressure());
    EXPECT_DOUBLE_EQ(restored.material().last_local_substeps,source.material().last_local_substeps);
    EXPECT_EQ(captureMFrontCheckpoint(restored.material()._behaviour_data),captureMFrontCheckpoint(source.material()._behaviour_data));
    EXPECT_NE(restored.mechanical_pressure,10000.);
    EXPECT_DOUBLE_EQ(restored.material()._behaviour_data.s0.mass_density,2050.);
    EXPECT_TRUE(restored.mechanical_drainage.hasReference());
    if (auto const path=std::getenv("HM_IP_CHECKPOINT_TEST_OUTPUT"))
    { std::ofstream out(path,std::ios::binary);out.write(capsule.data(),capsule.size());ASSERT_TRUE(out.good()); }
    int branch_index=0;
    for (double const gamma : {.0015,-.001,.003})
    {
        CheckpointHostIP direct{b}, branch{b};
        // Uninterrupted control is a native copy of the source, not a decode.
        direct.sigma_eff=source.sigma_eff;direct.sigma_eff_prev=source.sigma_eff_prev;
        direct.eps=source.eps;direct.eps_prev=source.eps_prev;
        direct.strain_rate_variable=source.strain_rate_variable;
        direct.mechanical_drainage=source.mechanical_drainage;
        direct.mechanical_pressure=source.mechanical_pressure;direct.mechanical_excess=source.mechanical_excess;
        direct.mechanical_fluid_content=source.mechanical_fluid_content;
        direct.material_state_variables=std::make_unique<MaterialStateVariablesMFront<2>>(source.material());
        HMIP::commitHMIntegrationPoint<2>(branch,HMIP::decodeHMIntegrationPoint<2>(restored,capsule,false));
        for (auto* ip : {&direct,&branch})
        {
            auto& d=ip->material()._behaviour_data;
            d.s1.gradients[3]=gamma/std::sqrt(2.);
            auto const n=integrateWithSubstepping(d,settings,integrate);ASSERT_TRUE(n);
            ip->material().last_local_substeps=*n;
            ip->sigma_eff=Eigen::Map<CheckpointHostIP::KV const>(d.s1.thermodynamic_forces.data());
            ip->eps=Eigen::Map<CheckpointHostIP::KV const>(d.s1.gradients.data());
        }
        double const G=1e7/(2*(1+.49));
        double const tau0=source.sigma_eff[3]/std::sqrt(2.);
        double const tau=std::clamp(tau0+G*(gamma-.002),-1000.,1000.);
        double const lambda=source.material().getEquivalentPlasticStrain()+std::abs(gamma-.002-(tau-tau0)/G);
        EXPECT_NEAR(branch.sigma_eff[3]/std::sqrt(2.),tau,1e-6);
        EXPECT_NEAR(branch.material().getEquivalentPlasticStrain(),lambda,1e-10);
        EXPECT_EQ(HMIP::captureHMIntegrationPoint<2>(branch,false),HMIP::captureHMIntegrationPoint<2>(direct,false));
        EXPECT_EQ(HMIP::captureHMIntegrationPoint<2>(restored,false),capsule);
        if (auto const path=std::getenv("HM_IP_CHECKPOINT_TEST_OUTPUT"))
        {
            auto const bytes=HMIP::captureHMIntegrationPoint<2>(branch,false);
            std::ofstream out(std::string(path)+".branch_"+std::to_string(branch_index),std::ios::binary);
            out.write(bytes.data(),bytes.size());ASSERT_TRUE(out.good());
        }
        ++branch_index;
    }
    std::vector<std::string> broken{capsule.substr(0,12),capsule+"x",capsule};
    broken.back()[30]^=1;
    for (auto const& bytes:broken)
    {
        EXPECT_ANY_THROW(HMIP::decodeHMIntegrationPoint<2>(restored,bytes,false));
        EXPECT_EQ(HMIP::captureHMIntegrationPoint<2>(restored,false),capsule);
    }
    EXPECT_ANY_THROW(HMIP::decodeHMIntegrationPoint<2>(restored,capsule,true));
    restored.integration_weight=.25;
    EXPECT_ANY_THROW(HMIP::decodeHMIntegrationPoint<2>(restored,capsule,false));
    restored.integration_weight=.125;
    EXPECT_EQ(HMIP::captureHMIntegrationPoint<2>(restored,false),capsule);
    CheckpointHostIP wrong_offset{b,3};
    EXPECT_ANY_THROW(HMIP::decodeHMIntegrationPoint<2>(wrong_offset,capsule,false));
}

TEST(HMIntegrationPointCheckpoint, TotalStressNaNAndAbsentReferenceStayDistinct)
{
    auto const b=mgis::behaviour::load("libOgsMFrontBehaviour","MohrCoulombAbboSloan",mgis::behaviour::Hypothesis::PLANESTRAIN);
    CheckpointHostIP source{b}, restored{b};
    double const nan=std::numeric_limits<double>::quiet_NaN();
    source.mechanical_pressure=source.mechanical_excess=source.mechanical_fluid_content=nan;
    auto const capsule=HMIP::captureHMIntegrationPoint<2>(source,true);
    HMIP::commitHMIntegrationPoint<2>(restored,HMIP::decodeHMIntegrationPoint<2>(restored,capsule,true));
    EXPECT_EQ(HMIP::captureHMIntegrationPoint<2>(restored,true),capsule);
    EXPECT_TRUE(std::isnan(restored.mechanical_pressure));
    EXPECT_FALSE(restored.mechanical_drainage.hasReference());
    EXPECT_ANY_THROW(HMIP::MechanicalDrainageState::fromCheckpoint(false,1.,0.));
    EXPECT_ANY_THROW(HMIP::MechanicalDrainageState::fromCheckpoint(true,nan,0.));
    auto const zero_reference=HMIP::MechanicalDrainageState::fromCheckpoint(true,0.,0.);
    EXPECT_TRUE(zero_reference.hasReference());
    EXPECT_EQ(zero_reference.referencePressure(),0.);
    EXPECT_ANY_THROW(HMIP::decodeHMIntegrationPoint<2>(restored,capsule,false));
    // Production TSA retains accepted mechanical reference strain, although
    // its pressure diagnostics remain invalid. Preserve this legitimate state.
    source.mechanical_drainage=HMIP::MechanicalDrainageState::fromCheckpoint(true,-.00006,0.);
    auto const mechanical_reference_capsule=HMIP::captureHMIntegrationPoint<2>(source,true);
    HMIP::commitHMIntegrationPoint<2>(restored,HMIP::decodeHMIntegrationPoint<2>(restored,mechanical_reference_capsule,true));
    EXPECT_TRUE(restored.mechanical_drainage.hasReference());
    EXPECT_EQ(restored.mechanical_drainage.referenceStrain(),-.00006);
    EXPECT_TRUE(std::isnan(restored.mechanical_pressure));
    EXPECT_EQ(HMIP::captureHMIntegrationPoint<2>(restored,true),mechanical_reference_capsule);
}
} // namespace MFrontSubsteppingTests

#include "ProcessLib/HydroMechanics/NativeCheckpoint.h"
namespace MFrontSubsteppingTests
{
TEST(HMNativeCheckpoint, EnvelopePreservesEveryProcessField)
{
    HMIP::NativeCheckpointState source;
    source.time=2.;source.initialisation_time=1.;source.id="Loaded";
    source.model_identity="model";source.runtime_identity="runtime";
    source.solution={1.,2.};source.initialisation_solution={3.,4.};
    source.displacement_reference={5.,6.};source.displacement_reference_next=2;
    source.nodal_forces={7.,8.};source.mass_flow_rate={9.};
    source.elements={{"first","second"},{"third"}};
    auto const bytes=HMIP::encodeNativeCheckpoint(source);
    auto const copy=HMIP::decodeNativeCheckpoint(bytes);
    EXPECT_EQ(copy.time,2.);EXPECT_EQ(copy.initialisation_time,1.);
    EXPECT_EQ(copy.id,"Loaded");EXPECT_EQ(copy.model_identity,"model");EXPECT_EQ(copy.runtime_identity,"runtime");
    EXPECT_EQ(copy.solution,source.solution);EXPECT_EQ(copy.initialisation_solution,source.initialisation_solution);
    EXPECT_EQ(copy.displacement_reference,source.displacement_reference);EXPECT_EQ(copy.displacement_reference_next,2u);
    EXPECT_EQ(copy.nodal_forces,source.nodal_forces);EXPECT_EQ(copy.mass_flow_rate,source.mass_flow_rate);
    EXPECT_EQ(copy.elements,source.elements);EXPECT_EQ(HMIP::encodeNativeCheckpoint(copy),bytes);
}
TEST(HMNativeCheckpoint, StructuralPayloadVersionAndLegacyCompatibility)
{
    HMIP::NativeCheckpointState source;
    source.time=2.;source.initialisation_time=1.;source.solution={1.};
    auto const legacy=HMIP::encodeNativeCheckpoint(source);
    EXPECT_EQ(legacy.substr(0,8),"OGSHMCP2");
    EXPECT_TRUE(HMIP::decodeNativeCheckpoint(legacy).source_states.empty());
    source.source_states={{1.,.25,.125},{},{1.,2.,3.}};
    auto const modern=HMIP::encodeNativeCheckpoint(source);
    EXPECT_EQ(modern.substr(0,8),"OGSHMCP3");
    auto const copy=HMIP::decodeNativeCheckpoint(modern);
    EXPECT_EQ(copy.source_states,source.source_states);
    EXPECT_EQ(HMIP::encodeNativeCheckpoint(copy),modern);
    source.source_states[0][1]=std::numeric_limits<double>::quiet_NaN();
    EXPECT_ANY_THROW(HMIP::decodeNativeCheckpoint(HMIP::encodeNativeCheckpoint(source)));
}
TEST(HMNativeCheckpoint, ExcavationPayloadPreservesReleaseAndLegacySources)
{
    HMIP::NativeCheckpointState source;
    source.time=7.;source.initialisation_time=1.;source.solution={1.,2.};
    source.excavation_states={{7.,1.},{1.,2.},{3.,4.},{5.,6.}};
    auto const bytes=HMIP::encodeNativeCheckpoint(source);
    EXPECT_EQ(bytes.substr(0,8),"OGSHMCP4");
    auto const copy=HMIP::decodeNativeCheckpoint(bytes);
    EXPECT_EQ(copy.excavation_states,source.excavation_states);
    EXPECT_TRUE(copy.source_states.empty());
    EXPECT_EQ(HMIP::encodeNativeCheckpoint(copy),bytes);
    source.source_states={{.125,42.}};
    EXPECT_EQ(HMIP::decodeNativeCheckpoint(HMIP::encodeNativeCheckpoint(source)).source_states,source.source_states);
    source.excavation_states.pop_back();
    EXPECT_ANY_THROW(HMIP::decodeNativeCheckpoint(HMIP::encodeNativeCheckpoint(source)));
}
TEST(HMNativeCheckpoint, EnvelopeRejectsCorruptionAndNonfiniteProcessFields)
{
    HMIP::NativeCheckpointState source;source.time=2.;source.initialisation_time=1.;source.solution={1.};
    auto const bytes=HMIP::encodeNativeCheckpoint(source);
    EXPECT_ANY_THROW(HMIP::decodeNativeCheckpoint(bytes.substr(0,10)));
    EXPECT_ANY_THROW(HMIP::decodeNativeCheckpoint(bytes+"trailing"));
    auto corrupt=bytes;corrupt.back()^=1;
    EXPECT_ANY_THROW(HMIP::decodeNativeCheckpoint(corrupt));
    source.time=std::numeric_limits<double>::quiet_NaN();
    EXPECT_ANY_THROW(HMIP::decodeNativeCheckpoint(HMIP::encodeNativeCheckpoint(source)));
    source.time=2.;source.solution[0]=std::numeric_limits<double>::infinity();
    EXPECT_ANY_THROW(HMIP::decodeNativeCheckpoint(HMIP::encodeNativeCheckpoint(source)));
}
}

#include "ProcessLib/BoundaryConditionAndSourceTerm/AcceptedFieldDirichletBoundaryCondition.h"
namespace MFrontSubsteppingTests
{
TEST(HMAcceptedPressure, AcceptedCaptureImmutableAcrossQueriesAndRestart)
{
    auto make=[] {
        auto bc=std::make_unique<ProcessLib::AcceptedFieldDirichletBoundaryCondition>(
            std::make_unique<TestIncrementBC>(),2.,4.);
        bc->setIndices({0,1});return bc;
    };
    auto bc=make();GlobalVector x(2);x.set(0,1700.);x.set(1,9500.);
    NumLib::IndexValueVector<GlobalIndexType> values;
    bc->getEssentialBCValues(2.,x,values);EXPECT_TRUE(values.ids.empty());
    EXPECT_ANY_THROW(bc->getEssentialBCValues(2.1,x,values));
    bc->postAcceptedTimestep(1.9,{&x},0);
    EXPECT_ANY_THROW(bc->getEssentialBCValues(2.1,x,values));
    bc->postAcceptedTimestep(2.,{&x},0);
    auto state=bc->checkpointState();
    x.set(0,-999.);bc->getEssentialBCValues(3.,x,values);
    EXPECT_EQ(values.values,(std::vector<double>{1700.,9500.}));
    bc->postAcceptedTimestep(3.,{&x},0);EXPECT_EQ(bc->checkpointState(),state);
    auto restored=make();restored->restoreCheckpointState(state,3.);
    restored->getEssentialBCValues(4.,x,values);EXPECT_EQ(values.values,(std::vector<double>{1700.,9500.}));
    restored->getEssentialBCValues(4.1,x,values);EXPECT_TRUE(values.ids.empty());
    auto corrupted=state;corrupted[3]=7.;
    EXPECT_ANY_THROW(restored->restoreCheckpointState(corrupted,3.));
    EXPECT_EQ(restored->checkpointState(),state);
    EXPECT_ANY_THROW(restored->validateCheckpointState(state,1.));
    EXPECT_ANY_THROW(make()->postAcceptedTimestep(2.1,{&x},0));
}
}

#include "MaterialLib/SolidModels/MFront/MCStrengthReduction.h"
namespace MFrontSubsteppingTests
{
TEST(MFrontStrengthReduction, BothEndpointPropertiesAndAtomicInvalidInput)
{
    auto b=mgis::behaviour::load("libOgsMFrontBehaviour","MohrCoulombAbboSloan",mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data data{b};data.s0.material_properties=data.s1.material_properties={1e7,.3,1000.,30.,5.,29.,10.};
    auto const original=data;
    applyMCStrengthReduction(data,1.5,3.);
    EXPECT_DOUBLE_EQ(data.s0.material_properties[2],1000./1.5);
    EXPECT_DOUBLE_EQ(data.s1.material_properties[2],1000./3.);
    constexpr double rad=std::numbers::pi/180.;
    EXPECT_NEAR(std::tan(data.s0.material_properties[3]*rad),std::tan(30.*rad)/1.5,1e-15);
    EXPECT_NEAR(std::tan(data.s1.material_properties[3]*rad),std::tan(30.*rad)/3.,1e-15);
    EXPECT_NEAR(std::tan(data.s1.material_properties[4]*rad),std::tan(5.*rad)/3.,1e-15);
    for (auto i:{0u,1u,5u,6u}) { EXPECT_EQ(data.s1.material_properties[i],original.s1.material_properties[i]); }
    auto const before=data;
    EXPECT_ANY_THROW(applyMCStrengthReduction(data,3.,2.));
    EXPECT_EQ(data.s0.material_properties,before.s0.material_properties);
    EXPECT_EQ(data.s1.material_properties,before.s1.material_properties);
    data=original;data.s1.material_properties[6]=101.;auto const bad=data;
    EXPECT_ANY_THROW(applyMCStrengthReduction(data,1.,2.));
    EXPECT_EQ(data.s0.material_properties,bad.s0.material_properties);
    EXPECT_EQ(data.s1.material_properties,bad.s1.material_properties);
}
TEST(MFrontStrengthReduction, FixedShearStrainReturnsToReducedStrengthWithoutReset)
{
    auto b=mgis::behaviour::load("libOgsMFrontBehaviour","MohrCoulombAbboSloan",mgis::behaviour::Hypothesis::PLANESTRAIN);
    auto integrate=[&](Data& d) { auto view=mgis::behaviour::make_view(d);return mgis::behaviour::integrate(view,b); };
    std::ofstream out;
    if (auto const path=std::getenv("HM_STRENGTH_REDUCTION_CSV"))
    { out.open(path);ASSERT_TRUE(out.good());out<<std::setprecision(17)<<"substeps,RF,gamma,tau,expected_tau,plastic_multiplier\n"; }
    double const young=1e7,nu=.49,gamma=.0002,shear=young/(2*(1+nu));
    for (unsigned const n:{1u,8u})
    {
        Data data{b};data.dt=1.;
        data.s0.material_properties=data.s1.material_properties={young,nu,1000.,0.,0.,29.,10.};
        data.s0.external_state_variables=data.s1.external_state_variables={293.15};
        data.s1.gradients[3]=gamma/std::sqrt(2.);
        SubsteppingSettings settings;settings.initial_substeps=settings.max_substeps=n;
        ASSERT_TRUE(integrateWithSubstepping(data,settings,integrate));mgis::behaviour::update(data);
        double previous=1.,last_lambda=0.;
        for (double const rf:{1.,1.5,2.,3.,4.})
        {
            data.s1=data.s0;
            data.s0.material_properties=data.s1.material_properties={young,nu,1000.,0.,0.,29.,10.};
            applyMCStrengthReduction(data,previous,rf);
            ASSERT_TRUE(integrateWithSubstepping(data,settings,integrate))<<n<<" "<<rf;
            double const tau=data.s1.thermodynamic_forces[3]/std::sqrt(2.);
            double const expected=std::min(shear*gamma,1000./rf);
            EXPECT_NEAR(tau,expected,1e-5);
            EXPECT_NEAR(data.s1.gradients[3],gamma/std::sqrt(2.),1e-15);
            EXPECT_GE(data.s1.internal_state_variables.back()+1e-12,last_lambda);
            EXPECT_NEAR(data.s1.internal_state_variables.back(),gamma-expected/shear,1e-9);
            if (out) { out<<n<<','<<rf<<','<<gamma<<','<<tau<<','<<expected<<','<<data.s1.internal_state_variables.back()<<'\n'; }
            last_lambda=data.s1.internal_state_variables.back();previous=rf;mgis::behaviour::update(data);
        }
    }
}
}

#include "ProcessLib/HydroMechanics/SafetyContinuation.h"
namespace MFrontSubsteppingTests
{
TEST(MFrontStrengthReduction, UnitFactorIsExactlyTheBasePropertyState)
{
    auto b=mgis::behaviour::load("libOgsMFrontBehaviour","MohrCoulombAbboSloan",mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data data{b};data.s0.material_properties=data.s1.material_properties={1e7,.3,1000.,30.,5.,29.,10.};
    auto const before=data;applyMCStrengthReduction(data,1.,1.);
    EXPECT_EQ(data.s0.material_properties,before.s0.material_properties);
    EXPECT_EQ(data.s1.material_properties,before.s1.material_properties);
}
TEST(HMSafetyContinuation, DimensionlessFactorIsSeparateFromFrozenPhysicalTime)
{
    HMIP::SafetyContinuation settings{2.,1.,"NotApplicable",{2.,3.,4.},{1.,1.5,2.}};
    EXPECT_FALSE(settings.active(2.));EXPECT_TRUE(settings.active(2.1));
    EXPECT_EQ(settings.factor(0.),1.);EXPECT_EQ(settings.factor(2.),1.);
    EXPECT_EQ(settings.factor(2.5),1.25);EXPECT_EQ(settings.factor(3.5),1.75);
    EXPECT_EQ(settings.factor(4.),2.);EXPECT_EQ(settings.physical_time,2.);
}
}

namespace MFrontSubsteppingTests
{
TEST(HMActualPressureSafety, UnequalPrimaryActualReferenceAndLocalConservation)
{
    HMIP::MechanicalDrainageState state;
    state.rebaseAccepted(-.001,10000.,13000.);
    auto const frozen=state.evaluate(10000.,-.00125,.8,
        {HMIP::MechanicalDrainageType::DrainedEffectiveStress,8e7});
    EXPECT_EQ(frozen.pressure,13000.);EXPECT_EQ(frozen.excess,3000.);
    EXPECT_NEAR(frozen.fluid_content,-.0002,1e-15);EXPECT_EQ(frozen.volumetric_stiffness,0.);
    auto const undrained=state.evaluate(10000.,-.00125,.8,
        {HMIP::MechanicalDrainageType::UndrainedEffectiveStress,8e7});
    EXPECT_NEAR(undrained.pressure,29000.,1e-10);EXPECT_NEAR(undrained.excess,19000.,1e-10);
    EXPECT_NEAR(undrained.fluid_content,0.,1e-15);EXPECT_NEAR(undrained.volumetric_stiffness,5.12e7,1e-7);
    EXPECT_EQ(state.referencePressure(),10000.);EXPECT_EQ(state.referenceActualPressure(),13000.);
    EXPECT_ANY_THROW(state.evaluate(13000.,-.00125,.8,
        {HMIP::MechanicalDrainageType::UndrainedEffectiveStress,8e7}));
    EXPECT_ANY_THROW(state.rebaseAccepted(-.001,10000.,std::numeric_limits<double>::quiet_NaN()));
    EXPECT_EQ(state.referenceActualPressure(),13000.);
}
TEST(HMActualPressureSafety, StorageMatchesConstantDensityFlowCompliance)
{
    using Type=HMIP::MechanicalDrainageType;
    double const bulk=1e7/(3.*(1.-2.*.3));
    auto const grain=HMIP::MechanicalDrainageMaterial{Type::UndrainedEffectiveStress,0.,0.}.resolved(.8,.3,bulk);
    EXPECT_NEAR(1./grain.biot_modulus,(.8-.3)*(1.-.8)/bulk,1e-22);
    auto const water=HMIP::MechanicalDrainageMaterial{Type::UndrainedEffectiveStress,0.,2.2e9}.resolved(.8,.3,bulk);
    EXPECT_NEAR(1./water.biot_modulus,.3/2.2e9+(.8-.3)*(1.-.8)/bulk,1e-22);
    EXPECT_ANY_THROW((HMIP::MechanicalDrainageMaterial{Type::UndrainedEffectiveStress,0.,0.}.resolved(1.,.3,bulk)));
    auto const alpha1=HMIP::MechanicalDrainageMaterial{Type::UndrainedEffectiveStress,0.,2.2e9}.resolved(1.,.3,bulk);
    EXPECT_NEAR(alpha1.biot_modulus,2.2e9/.3,1e-5);
}
TEST(HMActualPressureSafety, WorkingModeRemainsUntilAcceptedSafetyBoundary)
{
    HMIP::MechanicalDrainageSettings settings{1.,{{0,{HMIP::MechanicalDrainageType::UndrainedEffectiveStress,8e7}}}};
    settings.safety_reference_time=2.;settings.safety_materials={{0,{HMIP::MechanicalDrainageType::DrainedEffectiveStress,8e7}}};
    EXPECT_EQ(settings.materialAt(0,2.).type,HMIP::MechanicalDrainageType::UndrainedEffectiveStress);
    EXPECT_EQ(settings.materialAt(0,2.1).type,HMIP::MechanicalDrainageType::DrainedEffectiveStress);
    EXPECT_TRUE(settings.isSafetyReference(2.));EXPECT_FALSE(settings.isSafetyReference(2.1));
}
TEST(HMActualPressureSafety, UnequalReferenceCheckpointRoundtripAndValidation)
{
    auto const b=mgis::behaviour::load("libOgsMFrontBehaviour","MohrCoulombAbboSloan",mgis::behaviour::Hypothesis::PLANESTRAIN);
    CheckpointHostIP source{b},target{b};source.mechanical_pressure=13000.;source.mechanical_excess=3000.;
    source.mechanical_drainage.rebaseAccepted(-.001,10000.,13000.);
    auto const bytes=HMIP::captureHMIntegrationPoint<2>(source,false);
    auto decoded=HMIP::decodeHMIntegrationPoint<2>(target,bytes,false);
    EXPECT_EQ(decoded.drainage.referenceActualPressure(),13000.);
    HMIP::commitHMIntegrationPoint<2>(target,std::move(decoded));
    EXPECT_EQ(HMIP::captureHMIntegrationPoint<2>(target,false),bytes);
    EXPECT_EQ(target.mechanical_drainage.referencePressure(),10000.);
    EXPECT_ANY_THROW(HMIP::MechanicalDrainageState::fromCheckpoint(false,0.,0.,13000.));
    EXPECT_ANY_THROW(HMIP::MechanicalDrainageState::fromCheckpoint(true,0.,0.,std::numeric_limits<double>::quiet_NaN()));
    EXPECT_EQ(HMIP::captureHMIntegrationPoint<2>(target,false),bytes);
}
}

namespace MFrontSubsteppingTests
{
TEST(HMActualPressureSafety, MaterialModesRetainDistinctSafetyResponses)
{
    using Type=HMIP::MechanicalDrainageType;
    HMIP::MechanicalDrainageSettings settings{1.,{
        {0,{Type::DrainedEffectiveStress,8e7}},
        {1,{Type::UndrainedEffectiveStress,8e7}},
        {2,{Type::UndrainedTotalStress,0.}}}};
    settings.safety_reference_time=4.;settings.safety_materials=settings.materials;
    HMIP::MechanicalDrainageState state;state.rebaseAccepted(-.001,10000.,13000.);
    auto const drained=state.evaluate(10000.,-.00125,.8,settings.materialAt(0,4.1));
    auto const undrained=state.evaluate(10000.,-.00125,.8,settings.materialAt(1,4.1));
    EXPECT_DOUBLE_EQ(drained.pressure,13000.);
    EXPECT_NEAR(undrained.pressure,29000.,1e-10);
    EXPECT_NEAR(undrained.fluid_content,0.,1e-15);
    EXPECT_TRUE(settings.materialAt(2,4.1).usesTotalStress());
}
}

#include "NumLib/TimeStepping/Algorithms/PhaseAdaptiveTimeStepping.h"
TEST(HMAdaptivePhases, GrowthBoundaryAndNewPhaseInitialIncrement)
{
    NumLib::PhaseAdaptiveTimeStepping stepper(0.,3.,{{0.,1.,.025,.0001,.1,true},{1.,3.,.2,.02,.4,true}});
    NumLib::TimeStep previous{NumLib::Time{0.}},current{NumLib::Time{0.}};
    EXPECT_DOUBLE_EQ(stepper.next(0.,0,previous,current),.025);
    NumLib::updateTimeSteps(.025,previous,current);
    EXPECT_DOUBLE_EQ(stepper.next(0.,2,previous,current),.0375);
    current=NumLib::TimeStep{NumLib::Time{.95},NumLib::Time{1.},10};
    EXPECT_DOUBLE_EQ(stepper.next(0.,2,previous,current),.2);
    current=NumLib::TimeStep{NumLib::Time{.9},NumLib::Time{.975},9};
    EXPECT_NEAR(stepper.next(0.,2,previous,current),.025,1e-15);
}
TEST(HMAdaptivePhases, RejectionRestoresClockAndHonoursMinimum)
{
    NumLib::PhaseAdaptiveTimeStepping stepper(0.,1.,{{0.,1.,.8,.1,.8,true}});
    NumLib::TimeStep previous{NumLib::Time{0.}},current{NumLib::Time{0.},NumLib::Time{.8},1};current.setAccepted(false);
    EXPECT_TRUE(stepper.canReduceTimestepSize(current,previous));
    EXPECT_DOUBLE_EQ(stepper.next(0.,10,previous,current),.4);
    EXPECT_DOUBLE_EQ(current.previous()(),0.);EXPECT_DOUBLE_EQ(current.current()(),.4);
    previous=NumLib::TimeStep{NumLib::Time{0.}};
    current=NumLib::TimeStep{NumLib::Time{0.},NumLib::Time{.1},1};current.setAccepted(false);
    EXPECT_FALSE(stepper.canReduceTimestepSize(current,previous));
    EXPECT_ANY_THROW(stepper.next(0.,10,previous,current));
}

#include "ProcessLib/HydroMechanics/SteadyStateSeepage.h"
TEST(HMSteadySeepage, IncomingIntervalOwnsBoundaryAndTransientRelease)
{
    ProcessLib::HydroMechanics::SteadyStateSeepage policy{{{0.,1.},{3.,4.}}};
    EXPECT_FALSE(policy.active(0.));
    EXPECT_TRUE(policy.active(.01));
    EXPECT_TRUE(policy.active(1.));
    EXPECT_FALSE(policy.active(1.01));
    EXPECT_FALSE(policy.active(3.));
    EXPECT_TRUE(policy.active(3.01));
    EXPECT_TRUE(policy.active(4.));
    EXPECT_FALSE(policy.active(4.01));
}

TEST(HMAdaptivePhases, CustomEffortRangeRetainsConvergedPlasticSteps)
{
    NumLib::PhaseAdaptiveTimeStepping stepper(0.,1.,{{0.,1.,.1,.001,.2,true,6,20}});
    NumLib::TimeStep previous{NumLib::Time{0.}},current{NumLib::Time{0.},NumLib::Time{.1},1};
    EXPECT_DOUBLE_EQ(stepper.next(0.,10,previous,current),.1);
    EXPECT_DOUBLE_EQ(stepper.next(0.,21,previous,current),.05);
    EXPECT_NEAR(stepper.next(0.,6,previous,current),.15,1e-15);
    current.setAccepted(false);
    EXPECT_DOUBLE_EQ(stepper.next(0.,10,previous,current),.05);
}

#include "NumLib/ODESolver/BacktrackingNewtonStrategy.h"
TEST(HMBacktracking, CubicOvershootFindsDecreasingStep)
{
    double const x=.1;
    double const f=x*x*x-1.;
    double const direction=f/(3*x*x);
    auto const step=NumLib::backtrackArmijo(1.,.5,1e-6,1e-4,21,.5*f*f,1e-24,
        [&](double a)->std::optional<double>{double const y=x-a*direction;double const r=y*y*y-1.;return .5*r*r;});
    ASSERT_TRUE(step);EXPECT_LT(step->first,1.);EXPECT_LT(step->second,.5*f*f);
}
TEST(HMBacktracking, InvalidTrialsNeverBecomeAcceptedMinimumStep)
{
    auto result=NumLib::backtrackArmijo(1.,.5,.01,1e-4,21,1.,1e-24,
        [](double)->std::optional<double>{return std::nullopt;});
    EXPECT_FALSE(result);
    result=NumLib::backtrackArmijo(1.,.5,.01,1e-4,21,1.,1e-24,
        [](double)->std::optional<double>{return 2.;});
    EXPECT_FALSE(result);
}

#include "ParameterLib/FunctionParameter.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

TEST(HMSafetyStaticField, TimeSymbolMetadataMatchesActualFunctionResponse)
{
    ParameterLib::SpatialPosition pos;
    pos.setCoordinates(MathLib::Point3d{{0., .25, 0.}});
    ParameterLib::FunctionParameter<double> spatial(
        "spatial_cohesion", {"100+50*(1-y)"}, {});
    ParameterLib::FunctionParameter<double> temporal(
        "temporal_cohesion", {"100+50*(1-y)+t"}, {});
    EXPECT_FALSE(spatial.usesTimeSymbol());
    EXPECT_TRUE(temporal.usesTimeSymbol());
    EXPECT_DOUBLE_EQ(spatial(0., pos)[0], 137.5);
    EXPECT_DOUBLE_EQ(spatial(123., pos)[0], 137.5);
    EXPECT_DOUBLE_EQ(temporal(123., pos)[0] - temporal(0., pos)[0], 123.);
    // Preserve the old conservative API; only the explicit Safety check uses
    // expression-symbol metadata.
    EXPECT_TRUE(spatial.isTimeDependent());
    EXPECT_TRUE(temporal.isTimeDependent());
}


TEST(HMMechanicalDrainage, ReferenceAndExcessRetainsUndrainedStiffnessAcrossPhreaticSurface)
{
    using namespace ProcessLib::HydroMechanics;
    constexpr double M = 7.333333333333333e9;
    MechanicalDrainageMaterial material{MechanicalDrainageType::UndrainedEffectiveStress, M, 0., true};
    for (double const initial : {-1.e6, 1.e4})
    {
        MechanicalDrainageState state;
        state.accept(1., 1., 0., initial);
        for (double const strain : {-1.e-7, 1.e-7})
        {
            auto const response = state.evaluate(initial, strain, 1., material);
            double const applied = material.appliedPressure(response.pressure, response.excess, initial, true);
            EXPECT_NEAR(applied, std::max(initial, 0.) - M * strain, 1.e-7);
            EXPECT_NEAR(response.fluid_content, 0., 1.e-20);
            EXPECT_EQ(material.pressureDerivative(response.pressure, true), 1.);
            constexpr double h = 1.e-10;
            auto const plus = state.evaluate(initial, strain+h, 1., material);
            auto const minus = state.evaluate(initial, strain-h, 1., material);
            double const difference = (material.appliedPressure(plus.pressure, plus.excess, initial, true) -
                material.appliedPressure(minus.pressure, minus.excess, initial, true))/(2*h);
            EXPECT_NEAR(difference, -M, M*1.e-7);
        }
    }
    material.reference_and_excess = false;
    EXPECT_EQ(material.appliedPressure(-1.e6, 1000., -1001000., true), 0.);
    EXPECT_EQ(material.pressureDerivative(-1.e6, true), 0.);
}

namespace MFrontTensionTests
{
using namespace MaterialLib::Solids::MFront;
using Data = mgis::behaviour::BehaviourData;
TEST(MFrontTension, HydrostaticTensionAndElasticUnload)
{
    auto b = mgis::behaviour::load("libOgsMFrontBehaviour", "MohrCoulombAbboSloanTension",
                                 mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data data{b}; data.dt=1.;
    data.s0.material_properties=data.s1.material_properties={1.e7,.3,250.,32.,0.,29.,100.,250.};
    data.s0.external_state_variables=data.s1.external_state_variables={293.15};
    data.s1.gradients={.001,.001,.001,0.};
    auto integrate=[&](Data& d){auto view=mgis::behaviour::make_view(d);return mgis::behaviour::integrate(view,b);};
    SubsteppingSettings settings;settings.initial_substeps=2;settings.max_substeps=128;settings.tangent_perturbation=1.e-8;
    ASSERT_TRUE(integrateWithSubstepping(data,settings,integrate));
    for(int i=0;i<3;++i) EXPECT_NEAR(data.s1.thermodynamic_forces[i],250.,1.e-3);
    EXPECT_NEAR(data.s1.thermodynamic_forces[3],0.,1.e-6);
    auto const saved=data.s1.internal_state_variables;
    data.s0=data.s1;data.s1=data.s0;
    for(int i=0;i<3;++i)data.s1.gradients[i]-=5.e-6;
    ASSERT_TRUE(integrateWithSubstepping(data,settings,integrate));
    for(int i=0;i<3;++i)EXPECT_NEAR(data.s1.thermodynamic_forces[i],125.,1.e-3);
    for(std::size_t i=4;i<saved.size();++i)EXPECT_NEAR(data.s1.internal_state_variables[i],saved[i],1.e-9);
}
TEST(MFrontTension, ZeroStrengthHydrostaticOpening)
{
    auto b=mgis::behaviour::load("libOgsMFrontBehaviour","MohrCoulombAbboSloanTension",mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data data{b};data.dt=1.;data.s0.material_properties=data.s1.material_properties={1.e7,.3,250.,32.,0.,29.,100.,0.};
    data.s0.external_state_variables=data.s1.external_state_variables={293.15};data.s1.gradients={1.e-4,1.e-4,1.e-4,0.};
    auto integrate=[&](Data& d){auto view=mgis::behaviour::make_view(d);return mgis::behaviour::integrate(view,b);};
    SubsteppingSettings settings;settings.initial_substeps=2;settings.max_substeps=128;settings.tangent_perturbation=1.e-8;
    ASSERT_TRUE(integrateWithSubstepping(data,settings,integrate));
    for(auto stress:data.s1.thermodynamic_forces)EXPECT_NEAR(stress,0.,1.e-3);
}
TEST(MFrontTension, ZeroCohesionSandCompressionAndOpening)
{
    auto const behaviour=mgis::behaviour::load(
        "libOgsMFrontBehaviour","MohrCoulombAbboSloanTension",
        mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data data{behaviour};data.dt=1.;
    data.s0.material_properties=data.s1.material_properties=
        {3.e7,.2,0.,34.,4.,29.,0.,0.};
    data.s0.external_state_variables=data.s1.external_state_variables={293.15};
    data.s1.gradients={-1.e-4,-1.e-4,0.,0.};
    auto integrate=[&](Data& d)
    { auto view=mgis::behaviour::make_view(d);return mgis::behaviour::integrate(view,behaviour); };
    SubsteppingSettings settings;settings.initial_substeps=2;
    settings.max_substeps=128;settings.tangent_perturbation=1.e-8;
    ASSERT_TRUE(integrateWithSubstepping(data,settings,integrate));
    EXPECT_NEAR(data.s1.thermodynamic_forces[0],-4166.666666666667,1.e-4);
    EXPECT_NEAR(data.s1.thermodynamic_forces[1],-4166.666666666667,1.e-4);
    EXPECT_NEAR(data.s1.thermodynamic_forces[2],-1666.666666666667,1.e-4);
    data.s0=data.s1;data.s1=data.s0;
    data.s1.gradients={1.e-4,1.e-4,0.,0.};
    ASSERT_TRUE(integrateWithSubstepping(data,settings,integrate));
    for(auto stress:data.s1.thermodynamic_forces)
    { EXPECT_NEAR(stress,0.,1.e-3); }
    for(auto value:data.s1.internal_state_variables)
    { EXPECT_TRUE(std::isfinite(value)); }
}
}
namespace MFrontTensionTests
{
TEST(MFrontTension, RotatedRankineReturnAndWholePathTangent)
{
    auto b=mgis::behaviour::load("libOgsMFrontBehaviour","MohrCoulombAbboSloanTension",mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data start{b};start.dt=1.;start.s0.material_properties=start.s1.material_properties={1.e7,.3,10000.,30.,0.,29.,100.,0.};
    start.s0.external_state_variables=start.s1.external_state_variables={293.15};
    start.s1.gradients={.000065,.000065,-.00013,.000195*std::sqrt(2.)};
    auto integrate=[&](Data& d){auto v=mgis::behaviour::make_view(d);return mgis::behaviour::integrate(v,b);};
    SubsteppingSettings settings;settings.initial_substeps=2;settings.max_substeps=128;settings.tangent_perturbation=1.e-8;
    Data data=start;ASSERT_TRUE(integrateWithSubstepping(data,settings,integrate));
    EXPECT_NEAR(data.s1.thermodynamic_forces[0],-928.571428571,1.e-3);
    EXPECT_NEAR(data.s1.thermodynamic_forces[1],-928.571428571,1.e-3);
    EXPECT_NEAR(data.s1.thermodynamic_forces[2],-1857.142857142,1.e-3);
    EXPECT_NEAR(data.s1.thermodynamic_forces[3],928.571428571*std::sqrt(2.),1.e-3);
    settings.tangent_perturbation=5.e-9;Data probe=start;ASSERT_TRUE(integrateWithSubstepping(probe,settings,integrate));
    for(std::size_t i=0;i<data.K.size();++i)EXPECT_NEAR(data.K[i],probe.K[i],5.);
}
TEST(MFrontTension, InactiveCapMatchesLegacyPlasticShear)
{
    SubsteppingSettings settings;settings.initial_substeps=8;settings.max_substeps=128;settings.tangent_perturbation=1.e-8;
    std::vector<double> reference_stress;double reference_lambda=0.;
    for(bool tension:{false,true})
    {
        auto b=mgis::behaviour::load("libOgsMFrontBehaviour",tension?"MohrCoulombAbboSloanTension":"MohrCoulombAbboSloan",mgis::behaviour::Hypothesis::PLANESTRAIN);
        Data data{b};data.dt=1.;std::vector<double> properties={1.e7,.3,10000.,0.,0.,29.,100.};if(tension)properties.push_back(0.);
        data.s0.material_properties=data.s1.material_properties=properties;
        data.s0.external_state_variables=data.s1.external_state_variables={293.15};data.s0.thermodynamic_forces=data.s1.thermodynamic_forces={-20000.,-20000.,-20000.,0.};data.s1.gradients={0.,0.,0.,.01};
        auto integrate=[&](Data& d){auto v=mgis::behaviour::make_view(d);return mgis::behaviour::integrate(v,b);};
        ASSERT_TRUE(integrateWithSubstepping(data,settings,integrate));
        EXPECT_NEAR(data.s1.thermodynamic_forces[3],10000.*std::sqrt(2.),1.e-3);
        if(!tension){reference_stress=data.s1.thermodynamic_forces;reference_lambda=data.s1.internal_state_variables[4];}
        else{for(int i=0;i<4;++i)EXPECT_NEAR(data.s1.thermodynamic_forces[i],reference_stress[i],1.e-3);EXPECT_NEAR(data.s1.internal_state_variables[4],reference_lambda,1.e-9);}
    }
}
TEST(MFrontTension, SimultaneousShearAndTensileReturn)
{
    auto b=mgis::behaviour::load("libOgsMFrontBehaviour","MohrCoulombAbboSloanTension",mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data data{b};data.dt=1.;data.s0.material_properties=data.s1.material_properties={1.e7,.3,10000.,0.,0.,29.,100.,0.};data.s0.external_state_variables=data.s1.external_state_variables={293.15};
    data.s1.gradients={.0026,-.00195,-.00065,0.};
    auto integrate=[&](Data& d){auto v=mgis::behaviour::make_view(d);return mgis::behaviour::integrate(v,b);};
    SubsteppingSettings settings;settings.initial_substeps=8;settings.max_substeps=128;settings.tangent_perturbation=1.e-8;
    ASSERT_TRUE(integrateWithSubstepping(data,settings,integrate));
    auto const& s=data.s1.thermodynamic_forces;
    EXPECT_LE(std::max({s[0],s[1],s[2]}),1.e-3);
    EXPECT_LE(std::max({s[0],s[1],s[2]})-std::min({s[0],s[1],s[2]}),20000.+1.e-2);
    EXPECT_GT(data.s1.internal_state_variables[4],0.);
    EXPECT_GT(data.s1.internal_state_variables[5]+data.s1.internal_state_variables[6]+data.s1.internal_state_variables[7],0.);
}
}
namespace MFrontTensionTests
{
TEST(MFrontTension, CompleteHostCapsulePreservesTensileHistoryAndRollback)
{
    using namespace MFrontSubsteppingTests;
    auto const b=mgis::behaviour::load("libOgsMFrontBehaviour","MohrCoulombAbboSloanTension",mgis::behaviour::Hypothesis::PLANESTRAIN);
    CheckpointHostIP source{b};auto& data=source.material()._behaviour_data;
    data.dt=1.;data.s0.material_properties=data.s1.material_properties={1.e7,.3,250.,32.,0.,29.,100.,250.};data.s0.external_state_variables=data.s1.external_state_variables={293.15};data.s1.gradients={.001,.001,.001,0.};
    auto integrate=[&](Data& d){auto v=mgis::behaviour::make_view(d);return mgis::behaviour::integrate(v,b);};
    SubsteppingSettings settings;settings.initial_substeps=2;settings.tangent_perturbation=1.e-8;
    ASSERT_TRUE(integrateWithSubstepping(data,settings,integrate));source.material().pushBackState();
    source.sigma_eff=Eigen::Map<CheckpointHostIP::KV const>(data.s1.thermodynamic_forces.data());source.eps=Eigen::Map<CheckpointHostIP::KV const>(data.s1.gradients.data());
    auto const capsule=HMIP::captureHMIntegrationPoint<2>(source,false);
    auto const saved=data.s1.internal_state_variables;
    for(int i=0;i<3;++i)data.s1.gradients[i]+=.0001;
    ASSERT_TRUE(integrateWithSubstepping(data,settings,integrate));source.material().pushBackState();
    auto parsed=HMIP::decodeHMIntegrationPoint<2>(source,capsule,false);HMIP::commitHMIntegrationPoint<2>(source,std::move(parsed));
    EXPECT_EQ(HMIP::captureHMIntegrationPoint<2>(source,false),capsule);
    EXPECT_EQ(source.material()._behaviour_data.s1.internal_state_variables,saved);
}
}

namespace MCElasticStiffnessTests
{
TEST(MFrontElasticStiffness, CohesionlessMCUsesDeclaredElasticModuli)
{
    using namespace MaterialLib::Solids::MFront;
    ParameterLib::ConstantParameter<double> young("E", 9.e7);
    ParameterLib::ConstantParameter<double> poisson("nu", .2);
    ParameterLib::ConstantParameter<double> zero("zero", 0.);
    ParameterLib::ConstantParameter<double> phi("phi", 34.);
    ParameterLib::ConstantParameter<double> psi("psi", 4.);
    ParameterLib::ConstantParameter<double> transition("transition", 29.);
    auto behaviour = mgis::behaviour::load(
        "libOgsMFrontBehaviour", "MohrCoulombAbboSloanTension",
        mgis::behaviour::Hypothesis::PLANESTRAIN);
    MFrontGeneric<2, boost::mp11::mp_list<Strain>,
        boost::mp11::mp_list<Stress>, boost::mp11::mp_list<Temperature>>
        material(std::move(behaviour),
            {&young, &poisson, &zero, &phi, &psi, &transition, &zero, &zero},
            {}, std::nullopt);
    ParameterLib::SpatialPosition position;
    auto const elastic = material.getMCElasticStiffness(1., position);
    ASSERT_TRUE(elastic);
    EXPECT_NEAR((*elastic)(0,0), 1.e8, 1.e-6);
    EXPECT_NEAR((*elastic)(0,1), 2.5e7, 1.e-6);
    EXPECT_NEAR((*elastic)(3,3), 7.5e7, 1.e-6);
    EXPECT_NEAR(material.getBulkModulus(1., position, &*elastic), 5.e7, 1.e-6);
}
}

namespace MFrontTensionTests
{
TEST(MFrontTension, CohesionlessRankineApexIsObjective)
{
    auto const behaviour = mgis::behaviour::load(
        "libOgsMFrontBehaviour", "MohrCoulombAbboSloanTension",
        mgis::behaviour::Hypothesis::PLANESTRAIN);
    
    double const fraction = 1.;
    for (double const angle : {0., .37})
    {
        Data data{behaviour}; data.dt = 1.;
        data.s0.material_properties = data.s1.material_properties =
            {3.e7, .2, 0., 34., 4., 29., 0., 0.};
        data.s0.external_state_variables = data.s1.external_state_variables = {293.15};
        double const c = std::cos(angle), s = std::sin(angle);
        data.s1.gradients = {3.e-4*c*c+2.e-4*s*s,
            3.e-4*s*s+2.e-4*c*c, 1.e-4, std::sqrt(2.)*1.e-4*s*c};
        auto integrate = [&](Data& d)
        { auto view=mgis::behaviour::make_view(d);return mgis::behaviour::integrate(view,behaviour); };
        SubsteppingSettings settings;settings.initial_substeps=2;
        settings.max_substeps=128;settings.tangent_perturbation=1.e-8;
        ASSERT_TRUE(integrateWithSubstepping(data,settings,integrate));
        for (auto const stress : data.s1.thermodynamic_forces)
        { EXPECT_NEAR(stress, 0., 1.e-3); }
        auto const lambda = mgis::behaviour::getVariableOffset(
            behaviour.isvs, "EquivalentPlasticStrain", behaviour.hypothesis);
        auto const tensile = mgis::behaviour::getVariableOffset(
            behaviour.isvs, "TensilePlasticStrain", behaviour.hypothesis);
        EXPECT_NEAR(data.s1.internal_state_variables[lambda], 0., 1.e-9);
        for (std::size_t i=0;i<4;++i)
        {
            double const expected = data.s1.gradients[i]*fraction;
            EXPECT_NEAR(data.s1.internal_state_variables[tensile+i],expected,1.e-9);
        }
    }
}
}

namespace MFrontTensionTests
{
TEST(MFrontTension, SingleIncrementRetainsRotatedWholePathTangent)
{
    auto const b=mgis::behaviour::load("libOgsMFrontBehaviour",
        "MohrCoulombAbboSloanTension",mgis::behaviour::Hypothesis::PLANESTRAIN);
    std::vector<double> reference;
    for(unsigned const n:{1u,2u})
    {
        Data data{b};data.dt=1.;
        data.s0.material_properties=data.s1.material_properties={1.e7,.3,10000.,30.,0.,29.,100.,0.};
        data.s0.external_state_variables=data.s1.external_state_variables={293.15};
        data.s1.gradients={.000065,.000065,-.00013,.000195*std::sqrt(2.)};
        unsigned calls=0;
        auto integrate=[&](Data& d){++calls;auto v=mgis::behaviour::make_view(d);return mgis::behaviour::integrate(v,b);};
        SubsteppingSettings settings;settings.initial_substeps=n;
        settings.max_substeps=n;settings.tangent_perturbation=1.e-8;
        ASSERT_TRUE(integrateWithSubstepping(data,settings,integrate));
        EXPECT_GE(calls,9u*n);
        EXPECT_NEAR(data.s1.thermodynamic_forces[0],-928.571428571,1.e-3);
        EXPECT_NEAR(data.s1.thermodynamic_forces[2],-1857.142857142,1.e-3);
        EXPECT_NEAR(data.s1.thermodynamic_forces[3],928.571428571*std::sqrt(2.),1.e-3);
        if(reference.empty()){reference=data.K;}
        else{for(std::size_t i=0;i<reference.size();++i)EXPECT_NEAR(data.K[i],reference[i],5.);}
    }
}
TEST(MFrontTension, SiltIncrementRetriesIfOneStepTangentProbeFails)
{
    auto const b=mgis::behaviour::load("libOgsMFrontBehaviour",
        "MohrCoulombAbboSloanTension",mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data start{b};start.dt=1.;
    start.s0.external_state_variables=start.s1.external_state_variables={293.15};
    start.s0.material_properties={60000000,0.20000000000000001,1000,30,0,29,10,0};
    start.s0.internal_state_variables={-0.00020445281926400512,-0.00012430535852950701,-5.4925756581060382e-05,0.00016617325205118866,0,0,0,0,0,0,0,0};
    start.s0.thermodynamic_forces={-16617.373202776467,-12610.000166051559,-9141.0200686292283,8308.6626025594323};
    start.s1.material_properties={60000000,0.20000000000000001,1000,30,0,29,10,0};
    start.s1.internal_state_variables={-0.00020445281926400512,-0.00012430535852950701,-5.4925756581060382e-05,0.00016617325205118866,0,0,0,0,0,0,0,0};
    start.s1.thermodynamic_forces={-16617.373202776467,-12610.000166051559,-9141.0200686292283,8308.6626025594323};
    start.s0.gradients={-0.0001495270626829452,9.539766779473411e-05,-0,0.00016617325205118855};
    start.s1.gradients={-0.00043716586244874403,0.00069255102825309368,0,0.00029367869835328277};
    unsigned calls=0;
    auto integrate=[&](Data& d){++calls;auto v=mgis::behaviour::make_view(d);return mgis::behaviour::integrate(v,b);};
    SubsteppingSettings settings;settings.initial_substeps=1;settings.max_substeps=1;
    settings.tangent_perturbation=1.e-8;Data result=start;
    // The primal one-step return succeeds, but one whole-path tangent probe
    // fails. Never publish that incomplete derivative or any trial history.
    EXPECT_FALSE(integrateWithSubstepping(result,settings,integrate));
    EXPECT_EQ(result.s1.gradients,start.s1.gradients);
    EXPECT_EQ(result.s1.thermodynamic_forces,start.s1.thermodynamic_forces);
    EXPECT_EQ(result.s1.internal_state_variables,start.s1.internal_state_variables);
    EXPECT_EQ(result.K,start.K);
    settings.max_substeps=512;
    auto const selected=integrateWithSubstepping(result,settings,integrate);
    ASSERT_TRUE(selected);EXPECT_GT(*selected,1u);
    for(auto const value:result.K)EXPECT_TRUE(std::isfinite(value));
    for(auto const value:result.s1.thermodynamic_forces)EXPECT_TRUE(std::isfinite(value));

}
}

TEST(HMExcavationRelease, TemporaryForceVanishesAtPhaseEnd)
{
    ProcessLib::HydroMechanics::ExcavationRelease const release{2.,4.};
    EXPECT_DOUBLE_EQ(release.remaining(2.),1.);
    EXPECT_DOUBLE_EQ(release.remaining(2.5),.75);
    EXPECT_DOUBLE_EQ(release.remaining(3.),.5);
    EXPECT_DOUBLE_EQ(release.remaining(4.),0.);
    EXPECT_DOUBLE_EQ(release.remaining(5.),0.);
}


TEST(HMInterface, CohesiveOpeningContinuityTangentAndHistory)
{
    using ProcessLib::StructuralElements::contactTrial;
    using ProcessLib::StructuralElements::ContactState;
    double constexpr kn=2e6,kt=1e6,c=50.,mu=.3;
    ContactState const virgin{};
    for (double const slip : {1e-5,1e-4,-1e-4,2e-4,-2e-4})
    {
        auto const open=contactTrial({-1e-12,slip},virgin,kn,kt,c,mu,true);
        auto const closed=contactTrial({1e-12,slip},virgin,kn,kt,c,mu,true);
        EXPECT_NEAR(open.traction[1],closed.traction[1],1e-6);
        EXPECT_DOUBLE_EQ(open.traction[0],0.);
        for(double const gap : {-1e-4,1e-4})
        {
            Eigen::Vector2d const u(gap,slip);
            auto const r=contactTrial(u,virgin,kn,kt,c,mu,true);
            for(int j=0;j<2;++j)
            {
                auto plus=u.eval(),minus=u.eval();plus[j]+=1e-9;minus[j]-=1e-9;
                Eigen::Vector2d const fd=(contactTrial(plus,virgin,kn,kt,c,mu,true).traction-
                    contactTrial(minus,virgin,kn,kt,c,mu,true).traction)/(2e-9);
                EXPECT_LT((fd-r.tangent.col(j)).norm(),1e-3);
            }
        }
    }
    auto const loaded=contactTrial({-1e-4,1e-4},virgin,kn,kt,c,mu,true);
    EXPECT_DOUBLE_EQ(loaded.traction[1],50.);
    EXPECT_DOUBLE_EQ(virgin.plastic_slip,0.);
    auto const unload=contactTrial({-1e-4,8e-5},loaded.state,kn,kt,c,mu,true);
    EXPECT_NEAR(unload.traction[1],30.,1e-12);
    auto const reclose=contactTrial({1e-4,8e-5},loaded.state,kn,kt,c,mu,true);
    EXPECT_NEAR(reclose.traction[1],30.,1e-12);
    for(double const gap : {-1e-4,0.,1e-4})
    {
        auto const legacy=contactTrial({gap,1e-4},loaded.state,kn,kt,0.,mu);
        auto const selected=contactTrial({gap,1e-4},loaded.state,kn,kt,0.,mu,true);
        EXPECT_EQ(legacy.traction,selected.traction);
        EXPECT_EQ(legacy.tangent,selected.tangent);
        EXPECT_DOUBLE_EQ(legacy.state.plastic_slip,selected.state.plastic_slip);
    }
}

namespace MFrontTensionTests
{
TEST(MFrontTension, CohesionlessSpectralReturnMatchesIndependentLodeRoot)
{
 auto const b=mgis::behaviour::load("libOgsMFrontBehaviour","MohrCoulombAbboSloanTension",mgis::behaviour::Hypothesis::PLANESTRAIN);
 Data data{b};data.dt=0.024218749999999734;
 data.s0.external_state_variables=data.s1.external_state_variables={293.15};
 data.s0.material_properties=data.s1.material_properties={90000000,0.20000000000000001,0,34,4,29,0,0};
 data.s0.thermodynamic_forces=data.s1.thermodynamic_forces={-5351.4646999258166,-1513.6555449126358,-5313.5019163882171,-22.185936274226044};
 data.s0.internal_state_variables=data.s1.internal_state_variables={-4.4289257862951622e-05,6.8815308705574566e-06,-4.3783087415783634e-05,-2.9581248365634726e-07,0.0026511710868914106,0,0,0,0,0,0,0};
 data.s0.gradients={-0.0011501372394887212,0.0018775939165045718,0,3.6575267530388991e-05};data.s1.gradients={-0.0013407785855562094,0.0021882580145790487,0,4.5893694801558128e-05};
 auto integrate=[&](Data& d){auto v=mgis::behaviour::make_view(d);return mgis::behaviour::integrate(v,b);};
 SubsteppingSettings settings;settings.initial_substeps=1;settings.max_substeps=1;settings.tangent_perturbation=1.e-8;
 ASSERT_TRUE(integrateWithSubstepping(data,settings,integrate));
 auto const& s=data.s1.thermodynamic_forces;
 double const delta=std::hypot(s[0]-s[1],std::sqrt(2.)*s[3]);
 std::array<double,3> principal{(s[0]+s[1]-delta)/2.,(s[0]+s[1]+delta)/2.,s[2]};
 std::sort(principal.begin(),principal.end());
 EXPECT_NEAR(principal[0],-163.99493169564562,1.e-5);
 EXPECT_NEAR(principal[1],-162.80172984646472,1.e-5);
 EXPECT_NEAR(principal[2],-46.38199658136328,1.e-5);
 EXPECT_NEAR(data.s1.internal_state_variables[4]-data.s0.internal_state_variables[4],.0005923435217095208,1.e-10);
 for(auto const v:data.K) EXPECT_TRUE(std::isfinite(v));
 // Elastic reversal must leave accumulated shear history unchanged.
 Data reverse=data;mgis::behaviour::update(reverse);
 auto const strain=reverse.s0.gradients;for(unsigned i=0;i<4;++i)reverse.s1.gradients[i]=strain[i]-1.e-10*(data.s1.gradients[i]-data.s0.gradients[i]);
 ASSERT_TRUE(integrateWithSubstepping(reverse,settings,integrate));
 EXPECT_NEAR(reverse.s1.internal_state_variables[4],reverse.s0.internal_state_variables[4],1.e-10);
}
}

namespace MFrontTensionTests
{
TEST(MFrontTension, MixedCohesionlessApexHasAdmissibleHistoryAndCompleteTangent)
{
 auto const b=mgis::behaviour::load("libOgsMFrontBehaviour","MohrCoulombAbboSloanTension",mgis::behaviour::Hypothesis::PLANESTRAIN);
 Data data{b};data.dt=7.3242187500044409e-05;
 data.s0.external_state_variables=data.s1.external_state_variables={293.15};
 data.s0.material_properties=data.s1.material_properties={90000000,0.20000000000000001,0,34,4,29,0,0};
 data.s0.thermodynamic_forces=data.s1.thermodynamic_forces={-5.3746035994832315,-1.5250739441613683,-5.3403288305866186,-0.17248773356598193};
 data.s0.internal_state_variables=data.s1.internal_state_variables={-4.4461367161484818e-08,6.8656949094733526e-09,-4.4004370242863319e-08,-2.2998364475464259e-09,0.0032017801536933792,0,0,0,0,0,0,0};
 data.s0.gradients={-0.001317632736060676,0.0021646071743574818,0,4.1980650038778774e-05};data.s1.gradients={-0.0013178277875678058,0.0021649389453224019,0,4.1976587248397526e-05};
 unsigned calls=0;auto integrate=[&](Data& d){++calls;auto v=mgis::behaviour::make_view(d);return mgis::behaviour::integrate(v,b);};
 SubsteppingSettings settings;settings.initial_substeps=1;settings.max_substeps=1;settings.tangent_perturbation=1.e-11;
 ASSERT_TRUE(integrateWithSubstepping(data,settings,integrate));EXPECT_EQ(calls,9u);
 for(auto const stress:data.s1.thermodynamic_forces)EXPECT_NEAR(stress,0.,1.e-7);
 for(auto const tangent:data.K)EXPECT_NEAR(tangent,0.,1.e-3);
 auto const& a=data.s0.internal_state_variables;auto const& z=data.s1.internal_state_variables;
 double const lambda=z[4]-a[4];ASSERT_GT(lambda,0.);
 double trace=0.;for(unsigned i=0;i<3;++i)trace+=data.s1.gradients[i]-data.s0.gradients[i]+a[i]-z[i]-(z[8+i]-a[8+i]);
 EXPECT_NEAR(trace/lambda,std::sin(4.*std::acos(-1.)/180.),1.e-9);
 for(unsigned i=5;i<8;++i)EXPECT_GE(z[i]-a[i],-1.e-18);
 EXPECT_GE(z[8]-a[8],-1.e-18);EXPECT_GE(z[9]-a[9],-1.e-18);EXPECT_GE(z[10]-a[10],-1.e-18);
 EXPECT_GE((z[8]-a[8])*(z[9]-a[9])-.5*(z[11]-a[11])*(z[11]-a[11]),-1.e-30);
 Data reload=data;mgis::behaviour::update(reload);reload.s1.gradients[0]-=1.e-7;reload.s1.gradients[1]-=1.e-7;
 ASSERT_TRUE(integrateWithSubstepping(reload,settings,integrate));
 EXPECT_NEAR(reload.s1.internal_state_variables[4],reload.s0.internal_state_variables[4],1.e-12);
 for(unsigned i=8;i<12;++i)EXPECT_NEAR(reload.s1.internal_state_variables[i],reload.s0.internal_state_variables[i],1.e-12);
 EXPECT_NEAR(reload.s1.thermodynamic_forces[0],-12.5,1.e-5);
 EXPECT_NEAR(reload.s1.thermodynamic_forces[2],-5.,1.e-5);
}
}

TEST(StructuralPlate, DistributedLoadRampHasBoundedAmplitude)
{
    using ProcessLib::StructuralElements::DistributedLoadRamp;
    DistributedLoadRamp const ramp(1.,2.);
    EXPECT_DOUBLE_EQ(ramp.factor(0.),0.);EXPECT_DOUBLE_EQ(ramp.factor(1.),0.);
    EXPECT_DOUBLE_EQ(ramp.factor(1.25),.25);EXPECT_DOUBLE_EQ(ramp.factor(1.5),.5);
    EXPECT_DOUBLE_EQ(ramp.factor(2.),1.);EXPECT_DOUBLE_EQ(ramp.factor(3.),1.);
    EXPECT_THROW(DistributedLoadRamp(1.,1.),std::invalid_argument);
    EXPECT_THROW(DistributedLoadRamp(2.,1.),std::invalid_argument);
    EXPECT_THROW(DistributedLoadRamp(0.,std::numeric_limits<double>::infinity()),std::invalid_argument);
}

TEST(HMAdaptivePhases, CapturedRoundoffRemainderLandsOnEndpoint)
{
    NumLib::PhaseAdaptiveTimeStepping stepper(5.,7.,{
        {5.,6.,.025,.000001,.1,true,6,20},
        {6.,7.,.025,.000001,.1,true,6,20}});
    double const t=5.993749999999981;
    NumLib::TimeStep previous{NumLib::Time{5.987499999999982}};
    NumLib::TimeStep current{NumLib::Time{5.987499999999982},NumLib::Time{t},131};
    double const dt=stepper.next(0.,9,previous,current);
    EXPECT_DOUBLE_EQ(dt,6.-t);
    NumLib::updateTimeSteps(dt,previous,current);
    EXPECT_EQ(current.current(),NumLib::Time{6.});
    EXPECT_DOUBLE_EQ(stepper.next(0.,9,previous,current),.025);
}
TEST(HMAdaptivePhases, EndpointSnapDoesNotSwallowResolvedConstructionOffset)
{
    NumLib::PhaseAdaptiveTimeStepping stepper(5.,7.,{
        {5.,6.,.025,.000001,.1,true,6,20},
        {6.,7.,.025,.000001,.1,true,6,20}});
    double const t=5.99;
    NumLib::TimeStep previous{NumLib::Time{5.98}};
    NumLib::TimeStep current{NumLib::Time{t-(.01-1e-8)},NumLib::Time{t},2};
    double const nominal=current.dt();
    double const dt=stepper.next(0.,9,previous,current);
    EXPECT_DOUBLE_EQ(dt,nominal);
    EXPECT_GT(6.-(t+dt),.9e-8);
}

namespace MFrontSubsteppingTests
{
TEST(HardeningSoilMC, ActualPlasticHistoryCapsuleAndRejectedTrialAreAtomic)
{
    auto const* path=std::getenv("HS_GENERIC_LIBRARY");
    if (!path) { GTEST_SKIP()<<"Set HS_GENERIC_LIBRARY and HS_KERNEL_LIBRARY for imported HS qualification."; }
    auto const b=mgis::behaviour::load(path,"PublishedHS",mgis::behaviour::Hypothesis::PLANESTRAIN);
    CheckpointHostIP source{b,-1};
    auto& d=source.material()._behaviour_data;
    d.dt=1.;
    d.s0.material_properties=d.s1.material_properties={30e6,30e6,90e6,.5,0.,34.,4.,.2,1e5,.4408070965292531,.9,0.,60e6/1.1,1.3,1.85,1.};
    d.s0.external_state_variables=d.s1.external_state_variables={293.15};
    d.s0.thermodynamic_forces=d.s1.thermodynamic_forces={-1e5,-1e5,-1e5,0.};
    d.s0.gradients=d.s1.gradients={0.,0.,0.,0.};
    d.s1.gradients={.00005,-.0002,.00005,0.};
    auto integrate=[&](Data& data){auto view=mgis::behaviour::make_view(data);return mgis::behaviour::integrate(view,b);};
    SubsteppingSettings settings; settings.max_substeps=128;
    auto const count=integrateWithSubstepping(d,settings,integrate);
    ASSERT_TRUE(count); source.material().last_local_substeps=*count;
    ASSERT_EQ(d.s1.internal_state_variables.size(),3u);
    ASSERT_EQ(d.s1.internal_state_variables[2],1.);
    ASSERT_LT(d.s1.internal_state_variables[1],0.);
    mgis::behaviour::update(d);
    source.sigma_eff=Eigen::Map<CheckpointHostIP::KV const>(d.s1.thermodynamic_forces.data());
    source.sigma_eff_prev=source.sigma_eff*.8;
    source.eps=Eigen::Map<CheckpointHostIP::KV const>(d.s1.gradients.data());
    source.eps_prev=source.eps*.7;
    auto const saved=HMIP::captureHMIntegrationPoint<2>(source,false);
    CheckpointHostIP restored{b,-1};
    HMIP::commitHMIntegrationPoint<2>(restored,HMIP::decodeHMIntegrationPoint<2>(restored,saved,false));
    EXPECT_EQ(HMIP::captureHMIntegrationPoint<2>(restored,false),saved);
    for (double strain : {-.0001,-.0004})
    {
        Data direct=d, resumed=restored.material()._behaviour_data;
        direct.s1.gradients[1]=resumed.s1.gradients[1]=strain;
        ASSERT_TRUE(integrateWithSubstepping(direct,settings,integrate));
        ASSERT_TRUE(integrateWithSubstepping(resumed,settings,integrate));
        EXPECT_EQ(captureMFrontCheckpoint(direct),captureMFrontCheckpoint(resumed));
    }
    Data rejected=d; auto const before=captureMFrontCheckpoint(rejected);
    auto poison=[](Data& bad){bad.s0.internal_state_variables[0]=999.;bad.s1.internal_state_variables[1]=999.;return -1;};
    EXPECT_FALSE(integrateWithSubstepping(rejected,settings,poison));
    EXPECT_EQ(captureMFrontCheckpoint(rejected),before);
    auto broken=saved; broken.back()^=1;
    EXPECT_ANY_THROW(HMIP::decodeHMIntegrationPoint<2>(restored,broken,false));
    EXPECT_EQ(HMIP::captureHMIntegrationPoint<2>(restored,false),saved);
    EXPECT_ANY_THROW(HMIP::captureHMIntegrationPoint<2>(source,true));
}
} // namespace MFrontSubsteppingTests

namespace MFrontSubsteppingTests
{
TEST(HardeningSoilMC, ConfiguredSingleStepTangentHasIndependentElasticReference)
{
    auto const* path=std::getenv("HS_GENERIC_LIBRARY");
    if (!path) { GTEST_SKIP()<<"Imported HS library required."; }
    auto const b=mgis::behaviour::load(path,"PublishedHS",mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data d{b};d.dt=1.;
    d.s0.material_properties=d.s1.material_properties={30e6,30e6,90e6,.5,0.,34.,4.,.2,1e5,.4408070965292531,.9,0.,60e6/1.1,1.3,1.85,1.};
    d.s0.external_state_variables=d.s1.external_state_variables={293.15};
    d.s0.thermodynamic_forces=d.s1.thermodynamic_forces={-1e5,-1e5,-1e5,0.};
    d.s0.internal_state_variables=d.s1.internal_state_variables={.01,-1e9,1.};
    d.s0.gradients=d.s1.gradients={0.,0.,0.,0.};
    SubsteppingSettings settings;settings.initial_substeps=1;settings.max_substeps=1;
    settings.tangent_perturbation=1e-9;
    bool requested_intermediate_tangent=false;
    auto integrate=[&](Data& data)
    {
        requested_intermediate_tangent|=data.K[0]!=0.;
        auto view=mgis::behaviour::make_view(data);
        return mgis::behaviour::integrate(view,b);
    };
    ASSERT_TRUE(integrateWithSubstepping(d,settings,integrate));
    EXPECT_FALSE(requested_intermediate_tangent);
    double const E=90e6,nu=.2,G=E/(2*(1+nu)),lambda=E*nu/((1+nu)*(1-2*nu));
    for(int i=0;i<4;++i)
    {
        for(int j=0;j<4;++j)
        {
            double const expected=i<3&&j<3 ? lambda+(i==j?2*G:0.) : (i==j?2*G:0.);
            // At isotropic stress, the least-principal-stress stiffness is
            // directionally nonsmooth. Finite-h curvature is O(E^2*h/pref),
            // rather than roundoff alone. Require one part per million.
            EXPECT_NEAR(d.K[i*4+j],expected,1e-6*std::max(E,std::abs(expected)));
        }
    }
}
} // namespace MFrontSubsteppingTests

#include "NumLib/ODESolver/NonlinearSystem.h"
#include "NumLib/Exceptions.h"
namespace
{
class BacktrackingAssemblyCounter final
    : public NumLib::NonlinearSystem<NumLib::NonlinearSolverTag::Newton>
{
public:
    explicit BacktrackingAssemblyCounter(bool const change_boundary,
                                         bool const screening = false,
                                         bool const reject_validation = false)
        : change_boundary_(change_boundary), screening_(screening),
          reject_validation_(reject_validation) {}
    unsigned assemblies = 0;
    unsigned residual_assemblies = 0;
    double assembled_free = 0.;
    bool isLinear() const override { return false; }
    bool requiresNormalization() const override { return false; }
    MathLib::MatrixSpecifications getMatrixSpecifications(int) const override
    { return {2, 2, nullptr, nullptr}; }
    void assemble(std::vector<GlobalVector*> const& x,
                  std::vector<GlobalVector*> const&, int) override
    {
        ++assemblies;
        assembled_free = (*x[0])[0];
        if (reject_validation_ && residual_assemblies && assembled_free > 1.2)
        { throw NumLib::AssemblyException("Injected full tangent failure"); }
    }
    bool supportsResidualScreening() const override { return screening_; }
    void assembleResidual(std::vector<GlobalVector*> const& x,
                          std::vector<GlobalVector*> const&, int) override
    { ++residual_assemblies; assembled_free = (*x[0])[0]; }
    std::vector<GlobalIndexType> getIndicesOfResiduumWithoutInitialCompensation() const override
    { return {}; }
    void setReleaseNodalForces(GlobalVector const*, int) override {}
    void getResidual(GlobalVector const&, GlobalVector const&,
                     GlobalVector& residual) const override
    { residual.set(0, assembled_free * assembled_free - 1.); residual.set(1, 0.); }
    void getJacobian(GlobalMatrix& matrix) const override
    { matrix.setZero(); matrix.getRawMatrix().coeffRef(0, 0) = 2. * assembled_free; matrix.getRawMatrix().coeffRef(1, 1) = 1.; }
    void computeKnownSolutions(GlobalVector const&, int) override {}
    void applyKnownSolutions(GlobalVector& x) const override
    { if (change_boundary_) { x.set(1, 5.); } }
    void applyKnownSolutionsNewton(GlobalMatrix&, GlobalVector& residual,
                                  GlobalVector const&, GlobalVector& delta) const override
    { residual.set(1, 0.); delta.set(1, 0.); }
    void applyKnownSolutionsPETScSNES(GlobalMatrix&, GlobalVector&,
                                     GlobalVector&) const override {}
    void updateConstraints(GlobalVector&, GlobalVector&, int) override {}
private:
    bool const change_boundary_;
    bool const screening_;
    bool const reject_validation_;
};

void checkBacktrackingAssemblyReuse(bool const changed_boundary)
{
    BacktrackingAssemblyCounter system(changed_boundary);
    GlobalVector x(2), previous(2), direction(2), residual(2), next(2);
    x.set(0, .5); x.set(1, 0.); previous.setZero(); direction.setZero();
    direction.set(0, -.75);
    std::vector<GlobalVector*> current{&x}, old{&previous};
    system.assemble(current, old, 0);
    GlobalMatrix equilibrated(2);
    equilibrated.getRawMatrix().coeffRef(0, 0) = 1234.; equilibrated.getRawMatrix().coeffRef(1, 1) = 1.;
    NumLib::NewtonStepContext context{system, old, 0, &current, nullptr};
    NumLib::BacktrackingNewtonStrategy strategy(1., .5, 1e-6, 1e-4, 21, 1e-9);
    auto const result = strategy.applyStep(x, direction, residual, equilibrated,
                                          next, context, 1);
    ASSERT_TRUE(result.success);
    EXPECT_DOUBLE_EQ(result.step_length, 1.);
    EXPECT_DOUBLE_EQ(next[0], 1.25);
    EXPECT_DOUBLE_EQ(next[1], changed_boundary ? 5. : 0.);
    // One original Newton assembly, one genuinely new trial; a changed
    // essential base additionally needs its own assembly.
    EXPECT_EQ(system.assemblies, changed_boundary ? 3u : 2u);
}
}  // namespace
TEST(HMBacktracking, UnchangedBaseReusesPhysicalAssembly)
{ checkBacktrackingAssemblyReuse(false); }
TEST(HMBacktracking, ChangedEssentialBaseIsReassembled)
{ checkBacktrackingAssemblyReuse(true); }

TEST(HMBacktracking, ResidualScreenRejectsOvershootAndFullyValidatesCandidate)
{
    BacktrackingAssemblyCounter system(false, true);
    GlobalVector x(2), previous(2), direction(2), residual(2), next(2);
    x.set(0, .1); x.set(1, 0.); previous.setZero(); direction.setZero();
    direction.set(0, -4.95);
    std::vector<GlobalVector*> current{&x}, old{&previous};
    system.assemble(current, old, 0);
    GlobalMatrix matrix(2);
    NumLib::NewtonStepContext context{system, old, 0, &current, nullptr};
    NumLib::BacktrackingNewtonStrategy strategy(1., .5, 1e-6, 1e-4, 21, 1e-9, true);
    auto const result = strategy.applyStep(x, direction, residual, matrix, next, context, 1);
    ASSERT_TRUE(result.success);
    EXPECT_DOUBLE_EQ(result.step_length, .25);
    EXPECT_EQ(system.residual_assemblies, 3u);
    EXPECT_EQ(system.assemblies, 2u);
    EXPECT_DOUBLE_EQ(system.assembled_free, next[0]);
}
TEST(HMBacktracking, ResidualScreenCannotAcceptFailedFullTangentValidation)
{
    BacktrackingAssemblyCounter system(false, true, true);
    GlobalVector x(2), previous(2), direction(2), residual(2), next(2);
    x.set(0, .1); x.set(1, 0.); previous.setZero(); direction.setZero();
    direction.set(0, -4.95);
    std::vector<GlobalVector*> current{&x}, old{&previous};
    system.assemble(current, old, 0);
    GlobalMatrix matrix(2);
    NumLib::NewtonStepContext context{system, old, 0, &current, nullptr};
    NumLib::BacktrackingNewtonStrategy strategy(1., .5, 1e-6, 1e-4, 21, 1e-9, true);
    auto const result = strategy.applyStep(x, direction, residual, matrix, next, context, 1);
    ASSERT_TRUE(result.success);
    EXPECT_DOUBLE_EQ(result.step_length, .125);
    EXPECT_EQ(system.residual_assemblies, 4u);
    EXPECT_EQ(system.assemblies, 3u);
    EXPECT_DOUBLE_EQ(system.assembled_free, next[0]);
}
namespace MFrontSubsteppingTests
{
TEST_F(MFrontSubstepping, ResidualPathRetainsStressHistoryAndTransactionalRetries)
{
    Data full = data, screened = data;
    unsigned full_calls = 0, residual_calls = 0;
    auto const a = integrateWithSubstepping(full, settings,
        [&](Data& d) { ++full_calls; return update(d); });
    auto const b = integrateWithSubstepping(screened, settings,
        [&](Data& d) { ++residual_calls; return update(d); }, false);
    ASSERT_TRUE(a); ASSERT_EQ(a, b);
    EXPECT_EQ(full.s1.thermodynamic_forces, screened.s1.thermodynamic_forces);
    EXPECT_EQ(full.s1.internal_state_variables, screened.s1.internal_state_variables);
    EXPECT_EQ(full.s0.gradients, screened.s0.gradients);
    EXPECT_LT(residual_calls, full_calls);
    Data rejected = data; auto const original = captureMFrontCheckpoint(rejected);
    EXPECT_FALSE(integrateWithSubstepping(rejected, settings,
        [](Data& d) { d.s0.gradients[0] = 999.; return -1; }, false));
    EXPECT_EQ(captureMFrontCheckpoint(rejected), original);
}
} // namespace MFrontSubsteppingTests

namespace MFrontSubsteppingTests
{
TEST(HardeningSoilMC, ResidualScreenKeepsActualPlasticStateWithoutTangentProbes)
{
    auto const* path = std::getenv("HS_GENERIC_LIBRARY");
    if (!path) { GTEST_SKIP() << "Imported HS library required."; }
    auto const b = mgis::behaviour::load(path, "PublishedHS",
                                       mgis::behaviour::Hypothesis::PLANESTRAIN);
    Data full{b}; full.dt = 1.;
    full.s0.material_properties = full.s1.material_properties =
        {30e6,30e6,90e6,.5,0.,34.,4.,.2,1e5,.4408070965292531,.9,0.,60e6/1.1,1.3,1.85,1.};
    full.s0.external_state_variables = full.s1.external_state_variables = {293.15};
    full.s0.thermodynamic_forces = full.s1.thermodynamic_forces = {-1e5,-1e5,-1e5,0.};
    full.s0.gradients = {0.,0.,0.,0.};
    full.s1.gradients = {.00005,-.0002,.00005,0.};
    Data screened = full;
    SubsteppingSettings settings;
    unsigned full_calls = 0, residual_calls = 0;
    auto call = [&](Data& d) { auto view = mgis::behaviour::make_view(d); return mgis::behaviour::integrate(view,b); };
    auto const a = integrateWithSubstepping(full, settings,
        [&](Data& d) { ++full_calls; return call(d); });
    auto const c = integrateWithSubstepping(screened, settings,
        [&](Data& d) { ++residual_calls; return call(d); }, false);
    ASSERT_TRUE(a); ASSERT_EQ(a,c);
    EXPECT_EQ(full.s1.thermodynamic_forces,screened.s1.thermodynamic_forces);
    EXPECT_EQ(full.s1.internal_state_variables,screened.s1.internal_state_variables);
    EXPECT_EQ(full.s0.gradients,screened.s0.gradients);
    EXPECT_LT(residual_calls,full_calls);
    EXPECT_EQ(screened.s1.internal_state_variables[2],1.);
}
} // namespace MFrontSubsteppingTests
TEST(HMBacktracking, ResidualScreeningRejectsUnsupportedEquationSystem)
{
    BacktrackingAssemblyCounter system(false);
    GlobalVector x(2), previous(2), direction(2), residual(2), next(2);
    x.setZero(); previous.setZero(); direction.setZero(); residual.setZero();
    std::vector<GlobalVector*> current{&x}, old{&previous};
    GlobalMatrix matrix(2);
    NumLib::NewtonStepContext context{system, old, 0, &current, nullptr};
    NumLib::BacktrackingNewtonStrategy strategy(1., .5, 1e-6, 1e-4, 21, 1e-9, true);
    EXPECT_ANY_THROW(strategy.applyStep(x, direction, residual, matrix, next, context, 1));
    EXPECT_EQ(system.assemblies, 0u);
    EXPECT_EQ(system.residual_assemblies, 0u);
}

#include "ProcessLib/StructuralElements/EmbeddedPileNetwork.h"

TEST(EmbeddedPileAcceptedState, RoundoffPerturbationsPreserveLocalEquilibrium)
{
using namespace ProcessLib::StructuralElements;
EmbeddedPileNetwork net({0.0,10.816653826391969,11.267347735824966,11.718041645257966,12.168735554690965,12.619429464123963,13.070123373556962,13.520817282989961,13.971511192422959,14.422205101855958},{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {8, 9}},0,9,{500000000.0,2800000.0,100000000.0,400000.0,1000000000.0,INFINITY,0.,0.},{0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0},{{500000000.0,0.0},{499748851.36979634,2811087.288955104},{499748851.36979634,2811087.288955104},{499748851.36979634,2811087.288955104},{499748851.36979634,2811087.288955104},{499748851.36979634,2811087.288955104},{499748851.36979634,2811087.288955104},{499748851.36979634,2811087.288955104},{499748851.36979634,2811087.288955104}});
auto state=net.initialState();
std::vector<double> data={2.0,10.0,8.000000000000002,-0.011449271336601384,0.0166693140486074,0.0,-0.017751419572049593,0.011236782635320082,-0.0007475258083181732,-0.01735844457612094,0.010882483520765313,-0.000764914420765777,-0.017034179837204275,0.010528339687313733,-0.0007667529213703612,-0.016759355858568,0.010213519758540296,-0.0005880527375004215,-0.01652654654188498,0.01001907357269146,-0.0003323109589328738,-0.016331293543998038,0.00989987593520605,-0.00028550149500631097,-0.01617562803254951,0.009759798059836835,-0.00038742976965998424,-0.016068601015520295,0.009574000226665111,-0.00048523589703931803,-0.016026305891263637,0.009355889085102733,-0.000526895992379028,-0.0029674750912184224,-0.0021475666305890936,0.0,0.00032144664226228695,-4.01095117036151e-05,0.0,0.0003390142060737781,-5.4890983397529896e-05,0.0,0.00035269184678970213,-6.140403192296983e-05,0.0,0.00036288292940848114,-5.9835564055955744e-05,0.0,0.00036986351031599455,-5.218063127930053e-05,0.0,0.0003745874199642929,-4.255016303495303e-05,0.0,0.00037739999030492996,-3.176350997355308e-05,0.0,0.00037830122133790574,-1.9820672095100673e-05,0.0,0.00037729111306322017,-6.721649399595817e-06,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.001817351115156172};
std::vector<double> sv={-0.011449271336601384,0.0166693140486074,0.0,-0.015574321847119597,0.011234646866788213,0.0,-0.015763668580052535,0.010885574265999362,0.0,-0.015903567664567048,0.010535381360713048,0.0,-0.015804656313494758,0.010200122033784044,0.0,-0.015658068594300703,0.01001490882435384,0.0,-0.015404346528483666,0.009907737404076356,0.0,-0.015025998883265688,0.009761324152365618,0.0,-0.014523025658646774,0.009575669069221636,0.0,-0.01389542685462692,0.009350772154644398,0.0};
for(int i=0;i<30;++i){state.q[i]=data[3+i];state.reference[i]=data[3+30+i];}
for(int i=0;i<10;++i){state.shaft[i].plastic_slip=data[3+60+i];state.lateral[i].plastic_slip=data[3+70+i];}
state.base.plastic_slip=data[83];state.prestress_strain=data[84];
Eigen::VectorXd soil=Eigen::Map<Eigen::VectorXd>(sv.data(),30);
PileControls c;c.external_force=Eigen::VectorXd::Zero(soil.size());c.bond_head_to_soil=true;

for(double delta : {0.,1e-15,-1e-15,1e-13,-1e-13,1e-11,-1e-11,1e-9,-1e-9})
{
auto trial=soil;trial[0]+=delta;
EXPECT_NO_THROW((void)net.trial(trial,state,c));
}
}
