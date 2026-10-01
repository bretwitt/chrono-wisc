// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2026 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: bgwitt
// =============================================================================
//
// Tests of SCMSindyResidual: the polynomial library, sparse recovery of a known
// model, saving and loading, the training-range clamps, the wheel-ground state
// on a slope, and SCM applying a registered force correction.
//
// =============================================================================

#include <cmath>
#include <cstdio>
#include <random>
#include <set>
#include <string>

#include "gtest/gtest.h"

#include "chrono/core/ChDataPath.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemSMC.h"

#include "chrono_vehicle/terrain/SCMSindyResidual.h"
#include "chrono_vehicle/terrain/SCMSindyStiffness.h"
#include "chrono_vehicle/terrain/SCMTerrain.h"

using namespace chrono;
using namespace chrono::vehicle;

namespace {

// Samples of a known sparse model, with a little noise
void AddSyntheticSamples(SCMSindyResidual& model, int n, double noise) {
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> slip(-0.2, 0.9), lat(-0.3, 0.3), slope(-0.1, 0.35), sink(0.0, 0.3);
    std::normal_distribution<double> eps(0, noise);
    for (int i = 0; i < n; ++i) {
        SCMSindyResidual::Sample s;
        s.x = {slip(rng), lat(rng), slope(rng), sink(rng)};
        const double S = s.x[0], B = s.x[1], T = s.x[2], Z = s.x[3];
        s.y[0] = -0.4 * S * T - 0.2 * T + eps(rng);  // lost traction on slopes, worse with slip
        s.y[1] = -0.3 * B + eps(rng);                // lateral resistance
        s.y[2] = 0.5 * Z * Z + eps(rng);             // stiffer soil at depth
        s.y[3] = 0.1 + eps(rng);                     // constant extra rolling resistance
        s.load = 100;
        model.AddSample(s);
    }
}

class SlopeFunctor : public SCMTerrain::HeightFunctor {
  public:
    explicit SlopeFunctor(double grade) : m_grade(grade) {}
    virtual double GetInitHeight(const ChVector2i& loc, double delta) override { return loc.x() * delta * m_grade; }

  private:
    double m_grade;
};

}  // namespace

TEST(SCMSindyResidual, Library) {
    SCMSindyResidual model(3);
    ASSERT_EQ(model.GetNumTerms(), 84);  // C(6 + 3, 3)
    EXPECT_EQ(model.GetTermName(0), "1");
    std::set<std::string> names;
    for (int k = 0; k < model.GetNumTerms(); ++k)
        names.insert(model.GetTermName(k));
    EXPECT_EQ(names.size(), 84u);
    EXPECT_TRUE(names.count("s*theta^2"));

    std::vector<double> theta;
    model.EvaluateLibrary({2, 3, 5, 7, 11, 13}, theta);
    for (int k = 0; k < model.GetNumTerms(); ++k)
        if (model.GetTermName(k) == "s*theta^2")
            EXPECT_DOUBLE_EQ(theta[k], 50.0);
}

TEST(SCMSindyResidual, SparseRecovery) {
    SCMSindyResidual model(3);
    AddSyntheticSamples(model, 4000, 0.002);
    const auto report = model.Fit();

    const int expected_terms[] = {2, 1, 1, 1};
    for (int i = 0; i < SCMSindyResidual::kNumOutputs; ++i) {
        EXPECT_EQ(report[i].num_terms, expected_terms[i]) << SCMSindyResidual::GetOutputName(i) << "\n" << model.GetEquations();
        EXPECT_LT(report[i].rmse, 0.003);
    }
    auto coef = [&](int i, const std::string& term) {
        for (int k = 0; k < model.GetNumTerms(); ++k)
            if (model.GetTermName(k) == term)
                return model.GetCoefficient(i, k);
        return std::nan("");
    };
    EXPECT_NEAR(coef(0, "s*theta"), -0.4, 0.01);
    EXPECT_NEAR(coef(0, "theta"), -0.2, 0.01);
    EXPECT_NEAR(coef(1, "b"), -0.3, 0.01);
    EXPECT_NEAR(coef(2, "z^2"), 0.5, 0.01);
    EXPECT_NEAR(coef(3, "1"), 0.1, 0.01);
}

TEST(SCMSindyResidual, SaveLoadAndClamp) {
    SCMSindyResidual model(2);
    AddSyntheticSamples(model, 2000, 0.0);
    model.Fit();
    const std::string file = GetChronoOutputPath() + "utest_scm_sindy_residual.json";
    model.Save(file);

    SCMSindyResidual loaded(2);
    loaded.Load(file);
    ASSERT_TRUE(loaded.IsTrained());
    const SCMSindyResidual::Features x = {0.5, 0.1, 0.2, 0.1};
    const auto y0 = model.Evaluate(x);
    const auto y1 = loaded.Evaluate(x);
    for (int i = 0; i < SCMSindyResidual::kNumOutputs; ++i)
        EXPECT_NEAR(y0[i], y1[i], 1e-12);

    // A model of another library degree is refused
    SCMSindyResidual other(3);
    EXPECT_THROW(other.Load(file), std::runtime_error);
    std::remove(file.c_str());

    // Far outside the training range, features are clamped: no blow-up of the polynomial
    const auto far = model.Evaluate({50, -50, 5, 10});
    const auto edge = model.Evaluate({0.9, -0.3, 0.35, 0.3});
    for (int i = 0; i < SCMSindyResidual::kNumOutputs; ++i)
        EXPECT_NEAR(far[i], edge[i], 0.02);
}

TEST(SCMSindyResidual, WheelState) {
    // Ground rising 10 degrees along +x, a wheel on it rolling uphill with 30% slip and 4 cm sinkage
    const double alpha = 10 * CH_DEG_TO_RAD, r = 0.25, sink = 0.04, v = 0.1, slip = 0.3;
    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -1.62));
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    auto scm = chrono_types::make_shared<SCMTerrain>(&sys, false);
    scm->Initialize(chrono_types::make_shared<SlopeFunctor>(std::tan(alpha)), 0.02);

    auto wheel = chrono_types::make_shared<ChBody>();
    wheel->SetPos(ChVector3d(0, 0, (r - sink) / std::cos(alpha)));
    const ChVector3d up_slope(std::cos(alpha), 0, std::sin(alpha));
    wheel->SetPosDt(v * up_slope);
    wheel->SetAngVelParent(ChVector3d(0, v / (1 - slip) / r, 0));
    sys.AddBody(wheel);

    SCMSindyResidual model(2);
    model.SetTerrain(scm.get());
    model.AddWheel(wheel, r);

    const ChVector3d n(-std::sin(alpha), 0, std::cos(alpha));
    SCMSindyResidual::State s;
    ASSERT_TRUE(model.GetState(*wheel, 200.0 * n, s));
    EXPECT_NEAR(s.x[0], slip, 1e-9);
    EXPECT_NEAR(s.x[1], 0, 1e-9);
    EXPECT_NEAR(s.x[2], alpha, 1e-3);
    EXPECT_NEAR(s.x[3], sink / r, 1e-3);
    EXPECT_NEAR(s.load, 200, 0.5);
    EXPECT_NEAR(s.l.Dot(up_slope), 1, 1e-6);

    // A trained model's correction is along the wheel-ground frame, in units of the force scale
    AddSyntheticSamples(model, 1000, 0.0);
    model.Fit();
    model.SetForceScale(150);
    ChVector3d df, dt;
    model.Correct(*wheel, 200.0 * n, VNULL, df, dt);
    const auto y = model.Evaluate({slip, 0, alpha, sink / r});
    EXPECT_NEAR(df.Dot(up_slope), 150 * y[0], 1e-6);
    EXPECT_NEAR(df.Dot(n), 150 * y[2], 1e-6);
    EXPECT_NEAR(dt.y(), 150 * r * y[3], 1e-6);
    EXPECT_LT(df.Dot(up_slope), 0);  // traction lost on the slope

    // Above the ground, the correction ramps out over the contact gap (2 cm), whether SCM touches the wheel or not
    wheel->SetPos(ChVector3d(0, 0, (r + 0.01) / std::cos(alpha)));
    ChVector3d df_half, dt_half;
    model.Correct(*wheel, VNULL, VNULL, df_half, dt_half);
    const auto y_half = model.Evaluate({slip, 0, alpha, -0.01 / r, 0});
    EXPECT_NEAR(df_half.Dot(up_slope), 0.5 * 150 * y_half[0], 1e-6);
    wheel->SetPos(ChVector3d(0, 0, (r + 0.03) / std::cos(alpha)));
    EXPECT_FALSE(model.GetState(*wheel, VNULL, s));
    model.Correct(*wheel, VNULL, VNULL, df_half, dt_half);
    EXPECT_TRUE(df_half.IsNull());
}

TEST(SCMSindyResidual, CorrectionOffContact) {
    // A body SCM does not touch, named by the correction: SCM still applies the correction to it
    class Lift : public SCMTerrain::ContactForceCorrection {
      public:
        explicit Lift(ChBody* body) : m_body(body) {}
        virtual void Correct(ChBody& body, const ChVector3d& force, const ChVector3d&, ChVector3d& dforce, ChVector3d& dtorque) override {
            dforce = ChVector3d(10, 0, 0);
            dtorque = VNULL;
            touched = touched || !force.IsNull();
        }
        virtual std::vector<ChBody*> GetBodies() const override { return {m_body}; }
        bool touched = false;

      private:
        ChBody* m_body;
    };

    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(VNULL);
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    auto block = chrono_types::make_shared<ChBodyEasyBox>(0.2, 0.2, 0.2, 1000, true, true, chrono_types::make_shared<ChContactMaterialSMC>());
    block->SetPos(ChVector3d(0, 0, 0.5));  // well above the ground
    sys.AddBody(block);
    SCMTerrain scm(&sys, false);
    scm.Initialize(4, 4, 0.02);
    auto lift = chrono_types::make_shared<Lift>(block.get());
    scm.RegisterContactForceCorrection(lift);
    for (int i = 0; i < 100; ++i)
        sys.DoStepDynamics(1e-3);
    EXPECT_FALSE(lift->touched);
    EXPECT_NEAR(block->GetPosDt().x(), 10.0 / block->GetMass() * 0.1, 1e-3);  // a = F / m over 0.1 s
}

TEST(SCMSindyResidual, CorrectionApplied) {
    // A block resting on SCM soil, with a correction pushing it along +x: it slides, while GetContactForceBody still
    // reports SCM's own force
    class Push : public SCMTerrain::ContactForceCorrection {
      public:
        virtual void Correct(ChBody& body, const ChVector3d& force, const ChVector3d&, ChVector3d& dforce, ChVector3d& dtorque) override {
            dforce = ChVector3d(2 * force.z(), 0, 0);
            dtorque = VNULL;
            ++calls;
        }
        int calls = 0;
    };

    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -9.81));
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    auto block = chrono_types::make_shared<ChBodyEasyBox>(0.4, 0.4, 0.2, 1000, true, true, chrono_types::make_shared<ChContactMaterialSMC>());
    block->SetPos(ChVector3d(0, 0, 0.1));
    sys.AddBody(block);

    SCMTerrain scm(&sys, false);
    scm.SetSoilParameters(2e6, 0, 1.1, 0, 30, 0.01, 4e7, 3e4);
    scm.Initialize(4, 4, 0.02);
    auto push = chrono_types::make_shared<Push>();
    scm.RegisterContactForceCorrection(push);

    for (int i = 0; i < 300; ++i)
        sys.DoStepDynamics(1e-3);

    EXPECT_GT(push->calls, 0);
    EXPECT_GT(block->GetPos().x(), 0.05);
    ChVector3d f, t;
    ASSERT_TRUE(scm.GetContactForceBody(block, f, t));
    EXPECT_LT(std::abs(f.x()), f.z());  // the reported force is SCM's, without the push
}

TEST(SCMSindyResidual, GravitySlopeOnFlatGround) {
    // Flat ground under gravity tilted back by 10 degrees, as CRM slope tests do: a wheel rolling along +x climbs
    const double alpha = 10 * CH_DEG_TO_RAD, r = 0.25;
    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(1.62 * ChVector3d(-std::sin(alpha), 0, -std::cos(alpha)));
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    auto scm = chrono_types::make_shared<SCMTerrain>(&sys, false);
    scm->Initialize(4, 4, 0.02);

    auto wheel = chrono_types::make_shared<ChBody>();
    wheel->SetPos(ChVector3d(0, 0, r - 0.02));
    sys.AddBody(wheel);
    SCMSindyResidual model;
    model.SetTerrain(scm.get());
    model.AddWheel(wheel, r);
    SCMSindyResidual::State s;
    ASSERT_TRUE(model.GetState(*wheel, ChVector3d(0, 0, 100.0), s));
    EXPECT_NEAR(s.x[2], alpha, 1e-9);
    EXPECT_NEAR(s.x[3], 0.02 / r, 1e-9);
}

TEST(SCMSindyResidual, FeatureMask) {
    // Without lateral slip (bit 1) or rut depth (bit 4): monomials of degree <= 2 in 3 features
    const unsigned int mask = 0xF & ~2u;
    SCMSindyResidual model(2, mask);
    ASSERT_EQ(model.GetNumTerms(), 10);
    for (int k = 0; k < model.GetNumTerms(); ++k)
        EXPECT_EQ(model.GetTermName(k).find('b'), std::string::npos) << model.GetTermName(k);

    AddSyntheticSamples(model, 500, 0.0);
    model.Fit();
    const std::string file = GetChronoOutputPath() + "utest_scm_sindy_mask.json";
    model.Save(file);
    SCMSindyResidual all(2);
    EXPECT_THROW(all.Load(file), std::runtime_error);
    auto loaded = SCMSindyResidual::FromFile(file);
    EXPECT_EQ(loaded->GetFeatureMask(), mask);
    const auto y0 = model.Evaluate({0.5, 0, 0.2, 0.1});
    const auto y1 = loaded->Evaluate({0.5, 0, 0.2, 0.1});
    for (int i = 0; i < SCMSindyResidual::kNumOutputs; ++i)
        EXPECT_NEAR(y0[i], y1[i], 1e-12);
    std::remove(file.c_str());
}

TEST(SCMSindyResidual, MirroredWheel) {
    // Two wheels in the same state, one mounted with its axle reversed (as Viper's left wheels): same features, and
    // the same correction in the absolute frame
    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -1.62));
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    auto scm = chrono_types::make_shared<SCMTerrain>(&sys, false);
    scm->Initialize(4, 4, 0.02);

    const double r = 0.25, v = 0.1, slip = 0.4;
    auto make_wheel = [&](const ChQuaternion<>& rot) {
        auto w = chrono_types::make_shared<ChBody>();
        w->SetPos(ChVector3d(0, 0, r - 0.03));
        w->SetRot(rot);
        w->SetPosDt(ChVector3d(v, 0, 0));
        w->SetAngVelParent(ChVector3d(0, v / (1 - slip) / r, 0));
        sys.AddBody(w);
        return w;
    };
    auto right = make_wheel(QUNIT);
    auto left = make_wheel(QuatFromAngleZ(CH_PI));

    SCMSindyResidual model(2);
    model.SetTerrain(scm.get());
    model.AddWheel(right, r);
    model.AddWheel(left, r);
    SCMSindyResidual::State sr, sl;
    ASSERT_TRUE(model.GetState(*right, ChVector3d(0, 0, 100.0), sr));
    ASSERT_TRUE(model.GetState(*left, ChVector3d(0, 0, 100.0), sl));
    for (int j = 0; j < SCMSindyResidual::kNumFeatures; ++j)
        EXPECT_NEAR(sr.x[j], sl.x[j], 1e-9) << SCMSindyResidual::GetFeatureName(j);
    EXPECT_NEAR(sr.x[0], slip, 1e-9);
    EXPECT_NEAR(sl.l.x(), 1, 1e-9);

    AddSyntheticSamples(model, 500, 0.0);
    model.Fit();
    ChVector3d fr, tr, fl, tl;
    model.Correct(*right, ChVector3d(0, 0, 100.0), VNULL, fr, tr);
    model.Correct(*left, ChVector3d(0, 0, 100.0), VNULL, fl, tl);
    EXPECT_NEAR((fr - fl).Length(), 0, 1e-9);
    EXPECT_NEAR((tr - tl).Length(), 0, 1e-9);
}

TEST(SCMSindyResidual, SlidingBackWhileSpinning) {
    // A wheel spinning forward while sliding back faster than its rim moves: slip past 1, not flipped to a skid
    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -1.62));
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    auto scm = chrono_types::make_shared<SCMTerrain>(&sys, false);
    scm->Initialize(4, 4, 0.02);
    const double r = 0.25;
    auto wheel = chrono_types::make_shared<ChBody>();
    wheel->SetPos(ChVector3d(0, 0, r - 0.03));
    wheel->SetPosDt(ChVector3d(-0.15, 0, 0));
    wheel->SetAngVelParent(ChVector3d(0, 0.1 / r, 0));
    sys.AddBody(wheel);
    SCMSindyResidual model(2);
    model.SetTerrain(scm.get());
    model.AddWheel(wheel, r);
    SCMSindyResidual::State s;
    ASSERT_TRUE(model.GetState(*wheel, ChVector3d(0, 0, 100.0), s));
    EXPECT_NEAR(s.x[0], (0.1 + 0.15) / 0.15, 1e-9);
    EXPECT_NEAR(s.l.x(), 1, 1e-9);
}

TEST(SCMSindyResidual, RutAhead) {
    // A 5 cm rut a wheel radius ahead of the wheel, as a rear wheel finds behind a front one
    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -1.62));
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    SCMTerrain scm(&sys, false);
    scm.Initialize(4, 4, 0.02);
    std::vector<SCMTerrain::NodeLevel> rut;
    for (int i = 5; i <= 30; ++i)
        for (int j = -8; j <= 8; ++j)
            rut.push_back({ChVector2i(i, j), -0.05});
    scm.SetModifiedNodes(rut);

    const double r = 0.25;
    auto wheel = chrono_types::make_shared<ChBody>();
    wheel->SetPos(ChVector3d(0, 0, r - 0.03));
    wheel->SetPosDt(ChVector3d(0.1, 0, 0));
    wheel->SetAngVelParent(ChVector3d(0, 0.1 / r, 0));
    sys.AddBody(wheel);
    SCMSindyResidual model(2);
    model.SetTerrain(&scm);
    model.AddWheel(wheel, r);
    SCMSindyResidual::State s;
    ASSERT_TRUE(model.GetState(*wheel, ChVector3d(0, 0, 100.0), s));
    EXPECT_NEAR(s.x[4], 0.05 / r, 1e-6);
    EXPECT_NEAR(s.x[3], 0.03 / r, 1e-6);

    // Travelling the other way, the ground ahead is fresh
    wheel->SetPosDt(ChVector3d(-0.1, 0, 0));
    wheel->SetAngVelParent(ChVector3d(0, -0.1 / r, 0));
    ASSERT_TRUE(model.GetState(*wheel, ChVector3d(0, 0, 100.0), s));
    EXPECT_NEAR(s.x[4], 0, 1e-6);
}

TEST(SCMSindyResidual, PhysicsTerms) {
    SCMSindyResidual model(1, 0x1F & ~2u);  // s, theta, z, d
    const int base = model.GetNumTerms();
    model.AddLibraryTerm("sat(s/0.1)");
    model.AddLibraryTerm("z^0.5");
    model.AddLibraryTerm("relu(s-1)");
    model.AddLibraryTerm("sat(s/0.3)*z");
    ASSERT_EQ(model.GetNumTerms(), base + 4);
    EXPECT_EQ(model.GetTermName(base), "sat(s/0.1)");
    EXPECT_EQ(model.GetTermName(base + 1), "z^0.5");
    EXPECT_EQ(model.GetTermName(base + 2), "relu(s-1)");
    EXPECT_EQ(model.GetTermName(base + 3), "sat(s/0.3)*z");

    std::vector<double> th;
    model.EvaluateLibrary({0.2, 0, 0, 0.09, 0}, th);
    EXPECT_NEAR(th[base], 1 - std::exp(-2.0), 1e-12);
    EXPECT_NEAR(th[base + 1], 0.3, 1e-12);
    EXPECT_NEAR(th[base + 2], 0, 1e-12);
    EXPECT_NEAR(th[base + 3], (1 - std::exp(-0.2 / 0.3)) * 0.09, 1e-12);
    model.EvaluateLibrary({-0.2, 0, 0, -0.01, 0}, th);  // odd saturation; fractional power of the positive part
    EXPECT_NEAR(th[base], -(1 - std::exp(-2.0)), 1e-12);
    EXPECT_NEAR(th[base + 1], 0, 1e-12);
    model.EvaluateLibrary({1.5, 0, 0, 0, 0}, th);
    EXPECT_NEAR(th[base + 2], 0.5, 1e-12);

    EXPECT_THROW(model.AddLibraryTerm("sat(s/0.1)"), std::invalid_argument);  // already there
    EXPECT_THROW(model.AddLibraryTerm("b"), std::invalid_argument);           // not among the features
    EXPECT_THROW(model.AddLibraryTerm("sat(s)"), std::invalid_argument);      // malformed
    EXPECT_THROW(model.AddLibraryTerm("q^2"), std::invalid_argument);         // unknown feature
}

TEST(SCMSindyResidual, PhysicsTermsParsimony) {
    // A saturating traction residual: a physics term captures it in one term where monomials need many
    auto make = [](bool physics) {
        auto m = std::make_shared<SCMSindyResidual>(3, 0x1F & ~2u);
        if (physics)
            m->AddLibraryTerm("sat(s/0.1)");
        std::mt19937 rng(11);
        std::uniform_real_distribution<double> s(-0.2, 1.0), z(0, 0.3);
        for (int i = 0; i < 2000; ++i) {
            SCMSindyResidual::Sample smp{};
            smp.x = {s(rng), 0, 0, z(rng), 0};
            smp.y = {0.4 * std::copysign(1 - std::exp(-std::abs(smp.x[0]) / 0.1), smp.x[0]) - 0.2 * smp.x[3], 0, 0, 0};
            smp.load = 100;
            m->AddSample(smp);
        }
        return m;
    };
    SCMSindyResidual::FitSettings fs;
    fs.threshold = 0.005;
    auto poly = make(false);
    auto phys = make(true);
    const auto rp = poly->Fit(fs);
    const auto rq = phys->Fit(fs);
    EXPECT_EQ(rq[0].num_terms, 2) << phys->GetEquations();
    EXPECT_LT(rq[0].rmse, 1e-6);
    EXPECT_GT(rp[0].num_terms, rq[0].num_terms);
    EXPECT_GT(rp[0].rmse, 10 * rq[0].rmse);

    // Saved and loaded with its library
    const std::string file = GetChronoOutputPath() + "utest_scm_sindy_physics.json";
    phys->Save(file);
    auto loaded = SCMSindyResidual::FromFile(file);
    ASSERT_EQ(loaded->GetNumTerms(), phys->GetNumTerms());
    const auto y0 = phys->Evaluate({0.3, 0, 0, 0.1, 0});
    const auto y1 = loaded->Evaluate({0.3, 0, 0, 0.1, 0});
    EXPECT_NEAR(y0[0], y1[0], 1e-12);
    std::remove(file.c_str());
}

TEST(SCMSindyResidual, SlipHistory) {
    // h follows slip through a first-order lag in the system's time
    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -1.62));
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    auto scm = chrono_types::make_shared<SCMTerrain>(&sys, false);
    scm->Initialize(4, 4, 0.02);
    const double r = 0.25;
    auto wheel = chrono_types::make_shared<ChBody>();
    wheel->SetPos(ChVector3d(0, 0, r - 0.02));
    wheel->SetPosDt(ChVector3d(0.1, 0, 0));
    wheel->SetAngVelParent(ChVector3d(0, 0.1 / r, 0));  // rolling: slip 0
    wheel->SetFixed(true);
    sys.AddBody(wheel);
    SCMSindyResidual model(2);
    model.SetTerrain(scm.get());
    model.AddWheel(wheel, r);
    model.SetHistoryTime(0.5);
    SCMSindyResidual::State s;
    ASSERT_TRUE(model.GetState(*wheel, VNULL, s));
    EXPECT_NEAR(s.x[5], 0, 1e-12);
    wheel->SetAngVelParent(ChVector3d(0, 0.2 / r, 0));  // slip 0.5 from now on
    sys.SetChTime(0.5);
    ASSERT_TRUE(model.GetState(*wheel, VNULL, s));
    EXPECT_NEAR(s.x[0], 0.5, 1e-9);
    EXPECT_NEAR(s.x[5], 0.5 * (1 - std::exp(-1.0)), 1e-9);
    ASSERT_TRUE(model.GetState(*wheel, VNULL, s));  // same time: no further update
    EXPECT_NEAR(s.x[5], 0.5 * (1 - std::exp(-1.0)), 1e-9);
}

TEST(SCMSindyResidual, OutputTermsWeightsEnsemble) {
    // Output 2 restricted to terms with d; the others free
    auto make = [](bool restrict_normal) {
        auto m = std::make_shared<SCMSindyResidual>(1, 0x3F & ~2u);
        if (restrict_normal)
            m->SetOutputTerms(2, [](const std::string& t) { return t.find('d') != std::string::npos; });
        std::mt19937 rng(5);
        std::uniform_real_distribution<double> u(0, 0.5);
        for (int i = 0; i < 400; ++i) {
            SCMSindyResidual::Sample smp{};
            smp.x = {u(rng), 0, u(rng), u(rng), u(rng), u(rng)};
            smp.y = {0.3 * smp.x[0], 0, 0.5 * smp.x[3] + 0.2 * smp.x[4], 0};
            smp.load = 100;
            m->AddSample(smp);
        }
        return m;
    };
    auto m = make(true);
    m->Fit();
    for (int k = 0; k < m->GetNumTerms(); ++k) {
        if (m->GetTermName(k).find('d') == std::string::npos)
            EXPECT_EQ(m->GetCoefficient(2, k), 0.0) << m->GetTermName(k);
    }
    EXPECT_NE(m->GetCoefficient(0, 1), 0.0);  // "s" free on output 0

    // A weight of 2 is the same as the sample twice
    SCMSindyResidual a(1, 0x3F & ~2u), b(1, 0x3F & ~2u);
    std::mt19937 rng(9);
    std::uniform_real_distribution<double> u(0, 1);
    for (int i = 0; i < 50; ++i) {
        SCMSindyResidual::Sample smp{};
        smp.x = {u(rng), 0, u(rng), u(rng), u(rng), u(rng)};
        smp.y = {smp.x[0] + 0.1 * u(rng), 0, 0, 0};
        smp.load = 100;
        const int copies = i % 3 == 0 ? 2 : 1;
        smp.weight = copies;
        a.AddSample(smp);
        smp.weight = 1;
        for (int c = 0; c < copies; ++c)
            b.AddSample(smp);
    }
    SCMSindyResidual::FitSettings fs;
    fs.threshold = 1e-6;
    a.Fit(fs);
    b.Fit(fs);
    for (int k = 0; k < a.GetNumTerms(); ++k)
        EXPECT_NEAR(a.GetCoefficient(0, k), b.GetCoefficient(0, k), 1e-9);

    // Ensemble: spread from members, saved and loaded with the allowed terms
    auto e = make(true);
    e->Fit();
    auto c0 = e->GetCoefficients(), c1 = c0;
    c1[0][1] += 0.2;  // member 1 differs in "s" on output 0
    e->SetEnsemble({c0, c1});
    const SCMSindyResidual::Features x = {0.25, 0, 0.25, 0.25, 0.25, 0.25};
    EXPECT_NEAR(e->EvaluateSpread(x)[0], 0.1 * 0.25, 1e-9);
    const std::string file = GetChronoOutputPath() + "utest_scm_sindy_ensemble.json";
    e->Save(file);
    auto l = SCMSindyResidual::FromFile(file);
    EXPECT_NEAR(l->EvaluateSpread(x)[0], 0.1 * 0.25, 1e-9);
    for (int k = 0; k < l->GetNumTerms(); ++k)
        EXPECT_EQ(l->IsTermAllowed(2, k), e->IsTermAllowed(2, k));
    std::remove(file.c_str());
}

TEST(SCMSindyStiffness, MultiplierUnderWheel) {
    // A map of constant ln multiplier ln 3: soil under the wheel is 3x as stiff, soil away from it keeps the base
    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -1.62));
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    SCMTerrain scm(&sys, false);
    scm.Initialize(4, 4, 0.02);
    const double r = 0.25;
    auto wheel = chrono_types::make_shared<ChBody>();
    wheel->SetPos(ChVector3d(0.5, 0, r - 0.02));
    sys.AddBody(wheel);

    auto map = std::make_shared<SCMSindyResidual>(1, 0x3F & ~2u);
    for (int i = 0; i < 40; ++i) {
        SCMSindyResidual::Sample smp{};
        smp.x = {0.01 * i, 0, 0, 0.005 * i, 0, 0.01 * i};
        smp.y = {0, 0, std::log(3.0), 0};
        smp.load = 100;
        map->AddSample(smp);
    }
    map->Fit();
    map->SetTerrain(&scm);
    map->AddWheel(wheel, r);
    SCMSindyStiffness stiffness(map, &scm, {1e5, 0, 1.1, 0, 30, 0.01, 2e7, 3e4});
    EXPECT_NEAR(stiffness.GetMultiplier(*wheel), 3.0, 1e-4);  // to the fit's ridge regularization

    double kphi, kc, n, c, phi, j, k, d;
    stiffness.Set(ChVector3d(0.5, 0.05, 0), kphi, kc, n, c, phi, j, k, d);  // under the wheel
    EXPECT_NEAR(kphi, 3e5, 30.0);
    EXPECT_NEAR(phi, 30, 1e-12);
    stiffness.Set(ChVector3d(-1.0, 0, 0), kphi, kc, n, c, phi, j, k, d);  // away from it
    EXPECT_NEAR(kphi, 1e5, 1e-6);
}
