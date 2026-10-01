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
// Sparse (SINDy) regression model of the soil forces on a wheel that SCM misses,
// trained on CRM, added to SCM's forces as a correction.
//
// =============================================================================

#ifndef SCM_SINDY_RESIDUAL_H
#define SCM_SINDY_RESIDUAL_H

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "chrono/physics/ChBody.h"

#include "chrono_vehicle/ChApiVehicle.h"
#include "chrono_vehicle/terrain/SCMTerrain.h"

namespace chrono {
namespace vehicle {

/// @addtogroup vehicle_terrain
/// @{

/// Residual model of the soil forces on a wheel: what a reference soil model (such as CRM) gives beyond what SCM
/// gives, as a sparse combination of polynomial terms in the wheel's state, found by sparse identification (SINDy:
/// sequentially thresholded least squares). SCM has no soil flow, so it misses how a slipping wheel digs in on a
/// slope and how soil pushed aside resists it; the residual carries those effects back into SCM.
///
/// The model works in a wheel-ground frame: n the ground's normal under the wheel (SCM's undeformed ground), a the
/// wheel's axle, l = a x n the direction the wheel rolls, t = n x l. The axle is taken the way the wheel spins (or,
/// barely spinning, the way it travels), so mirrored wheels and wheels driving in reverse see the same features. Its
/// features, all dimensionless:
/// - s: longitudinal slip, (w r - v_l) / max(|w r|, |v_l|, v_eps): from -1 (locked, skidding) through 0 (rolling) and
///   1 (spinning in place) to 2 (spinning forward while sliding back as fast)
/// - b: lateral slip, v_t / max(|w r|, |v_l|, v_eps), clamped to [-1, 1]
/// - theta: uphill slope of the ground along l, against gravity (rad): slopes made by tilting gravity over flat ground
///   count as slopes
/// - z: sinkage below the undeformed ground over the wheel radius
/// - d: depth of the rut in SCM's ground a wheel radius ahead of the wheel, over the wheel radius: a wheel running in
///   another's rut (a rover's rear wheels) sits deep below the undeformed ground but on soil already compacted
/// - h: slip low-pass filtered over a history time (SetHistoryTime, default 0.5 s): shear builds in the soil with the
///   displacement a slipping wheel accumulates (Janosi-Hanamoto), so traction lags slip; h - s carries that history
/// Its outputs, per unit of a constant force scale F (SetForceScale; the nominal wheel load, say):
/// - dFl/F, dFt/F, dFn/F: residual force along l, t and n
/// - dTa/(F r): residual torque about the axle
/// The residual is a function of the wheel's state alone. It is not scaled by SCM's load: in training, SCM sees the
/// reference model's sinkage, where its load can be far from the wheel's, and in closed loop the corrected normal force
/// settles the wheel at the reference model's sinkage. The correction applies while the wheel is at or below the
/// undeformed ground, and ramps out as it rises a contact gap (SetContactGap) above it; SCM need not touch the wheel. Where
/// SCM's ground is lower than the reference model's (SCM compacts a rut the reference model's soil dilates into), the
/// correction is the reference model's whole force.
///
/// Each output can be given its own candidate terms (SetOutputTerms), samples can be weighted (Sample::weight), and a
/// bootstrap ensemble of fits can be stored with the model (SetEnsemble): its spread, and how far a wheel's state lies
/// outside the training range, tell how far to trust the correction (GetWheelTrust).
///
/// Train with AddSample (PlanetHybridTerrain::SetResidualRecorder collects them while CRM carries the wheels) and Fit,
/// save with Save; apply by registering the model with the SCM terrain (SCMTerrain::RegisterContactForceCorrection)
/// and adding the wheels. Outside the range it was trained over the model is not evaluated: features are clamped to
/// the training range and each output to the largest magnitude seen in training.
class CH_VEHICLE_API SCMSindyResidual : public SCMTerrain::ContactForceCorrection {
  public:
    static constexpr int kNumFeatures = 6;
    static constexpr int kNumOutputs = 4;
    using Features = std::array<double, kNumFeatures>;
    using Outputs = std::array<double, kNumOutputs>;

    /// State of a wheel in its wheel-ground frame.
    struct State {
        Features x;        ///< features: s, b, theta, z, d, h
        double load;       ///< SCM normal load (N)
        double radius;     ///< wheel radius (m)
        ChVector3d l, t, n, a;  ///< wheel-ground frame (absolute)
    };

    /// A training sample: features, residual outputs (see class description), and SCM's normal load (N).
    struct Sample {
        Features x;
        Outputs y;
        double load;
        double weight = 1;  ///< weight in the fit (say, one over the number of samples from the same test)
    };

    /// Settings of the sparse regression.
    struct FitSettings {
        double threshold = 0.01;  ///< terms contributing less than this (in output units, RMS over the samples) are dropped
        double ridge = 1e-6;      ///< Tikhonov regularization of the normalized least-squares problem
        int max_iterations = 20;  ///< thresholding iterations
    };

    /// Quality of the fit of one output, over the training samples.
    struct FitReport {
        int num_terms;  ///< terms kept
        double r2;      ///< coefficient of determination
        double rmse;    ///< root mean square error (output units)
    };

    /// Model with a library of all monomials in the features up to total degree `degree` (1 to 5), in the features
    /// flagged in `features` (bit j for feature j; default: all). Leave out features the training data does not vary,
    /// such as lateral slip on a single-wheel rig: their terms would only fit noise.
    explicit SCMSindyResidual(int degree = 3, unsigned int features = 0x3F);

    /// Add a term to the library (before Fit): a product of factors separated by '*', each one of
    /// - `x` or `x^p`: a power of feature x (a non-integer p applies to max(x, 0), as a Bekker sinkage exponent),
    /// - `sat(x/c)`: sign(x) (1 - exp(-|x| / c)), a response saturating over a scale c (Janosi-Hanamoto shear),
    /// - `relu(x-a)`: max(x - a, 0), a hinge where a regime starts (a wheel spinning past slip 1),
    /// with x a feature name (s, b, theta, z, d, h) among the model's features, such as "sat(s/0.1)*z^0.5". Terms from
    /// physics can describe the residual with far fewer terms than monomials. Throws on a malformed term.
    void AddLibraryTerm(const std::string& term);

    /// Restrict output `i` to the library terms `allow` accepts, by name (default: all terms; saved with the model). An
    /// output's physics can call for its own terms: the normal force in a rut, say, from rut terms only.
    void SetOutputTerms(int i, const std::function<bool(const std::string& term)>& allow);
    bool IsTermAllowed(int i, int k) const { return m_allowed[i][k]; }

    /// Time constant of the slip history feature h (s, default: 0.5).
    void SetHistoryTime(double tau) { m_history_time = tau; }

    /// A bootstrap ensemble of fits (coefficients per output, per member), stored with the model; Evaluate uses the model's
    /// own coefficients, EvaluateSpread the members' standard deviation.
    using Coefficients = std::array<std::vector<double>, kNumOutputs>;
    void SetEnsemble(const std::vector<Coefficients>& members) { m_ensemble = members; }
    const std::vector<Coefficients>& GetEnsemble() const { return m_ensemble; }
    Coefficients GetCoefficients() const { return m_coefs; }

    /// Standard deviation of the ensemble members' outputs at `x` (zero without an ensemble).
    Outputs EvaluateSpread(const Features& x) const;

    /// How far to trust the correction on a wheel, as of its last correction: the largest ensemble spread over the
    /// applied outputs (output units), and how far its state lies outside the training range (the largest excess over a
    /// feature's range, in units of that range). Returns false if the wheel was not corrected.
    bool GetWheelTrust(const ChBody& wheel, double& spread, double& outside) const;

    /// Evaluate corrections on the wheels of `terrain`, which must outlive the model's use (it is not owned: the
    /// terrain holds the model once registered).
    void SetTerrain(const SCMTerrain* terrain) { m_terrain = terrain; }

    /// A wheel to correct and sample: its radius (m) and axle, in its reference frame.
    void AddWheel(std::shared_ptr<ChBody> wheel, double radius, const ChVector3d& axis = ChVector3d(0, 1, 0));

    /// Forget the wheels added (to reuse the model in another simulation).
    void ClearWheels() {
        m_wheels.clear();
        m_wheel_bodies.clear();
    }

    /// Force the outputs are scaled by (N, default: 100); saved with the model. Set before training.
    void SetForceScale(double force) { m_force_scale = force; }
    double GetForceScale() const { return m_force_scale; }

    /// Height of the wheel's lowest point above the undeformed ground over which the correction ramps out, and above
    /// which the wheel is neither corrected nor sampled (m, default: 0.02).
    void SetContactGap(double gap) { m_contact_gap = gap; }
    double GetContactGap() const { return m_contact_gap; }

    /// Speed below which slips are not resolved (m/s, default: 0.02).
    void SetSpeedEps(double eps) { m_speed_eps = eps; }

    /// Scale of the correction applied (default: 1; 0 turns it off).
    void SetGain(double gain) { m_gain = gain; }
    double GetGain() const { return m_gain; }

    /// Whether output `i` is applied (default: all; saved with the model). An output can be fitted but left out, such as
    /// the normal force where SCM's own is calibrated to the reference model and the residual is not a function of the
    /// wheel's state (at one sinkage, a stiff reference soil carries whatever load the wheel has).
    void SetOutputApplied(int i, bool applied) { m_applied[i] = applied; }
    bool IsOutputApplied(int i) const { return m_applied[i]; }

    /// State of a wheel (added with AddWheel) given SCM's force on it (absolute; zero if SCM does not touch it). Returns
    /// false if the wheel is not known, SCM's terrain is not set, or the wheel is more than the contact gap above the
    /// undeformed ground.
    bool GetState(const ChBody& wheel, const ChVector3d& scm_force, State& state) const;

    /// Add a sample from SCM's and the reference model's wrenches on a wheel (absolute frame, force at the center of
    /// mass, torque about it). Returns false if no sample was taken (see GetState).
    bool AddSample(const ChBody& wheel,
                   const ChVector3d& scm_force,
                   const ChVector3d& scm_torque,
                   const ChVector3d& ref_force,
                   const ChVector3d& ref_torque);

    /// The sample AddSample would add, without adding it (to average samples over time before adding them).
    bool MakeSample(const ChBody& wheel,
                    const ChVector3d& scm_force,
                    const ChVector3d& scm_torque,
                    const ChVector3d& ref_force,
                    const ChVector3d& ref_torque,
                    Sample& sample) const;

    /// Add a sample directly.
    void AddSample(const Sample& sample) { m_samples.push_back(sample); }

    const std::vector<Sample>& GetSamples() const { return m_samples; }
    void ClearSamples() { m_samples.clear(); }

    /// Write the samples as CSV (s, b, theta, z, d, load, dFl/F, dFt/F, dFn/F, dTa/(F r)), for offline study.
    void WriteSamples(const std::string& filename) const;

    /// Fit the model to the samples. Also sets the training range the model is clamped to.
    std::array<FitReport, kNumOutputs> Fit(const FitSettings& settings);
    std::array<FitReport, kNumOutputs> Fit() { return Fit(FitSettings()); }

    /// Whether the model has coefficients (from Fit or Load).
    bool IsTrained() const { return m_trained; }

    /// Library: number of terms, and the name of a term (such as "s*theta^2").
    int GetNumTerms() const { return static_cast<int>(m_terms.size()); }
    std::string GetTermName(int k) const;
    int GetDegree() const { return m_degree; }
    unsigned int GetFeatureMask() const { return m_features; }

    /// Evaluate the library at `x` (GetNumTerms values).
    void EvaluateLibrary(const Features& x, std::vector<double>& theta) const;

    /// Model outputs at `x`, with the training-range clamps.
    Outputs Evaluate(const Features& x) const;

    /// Coefficient of term `k` for output `i`.
    double GetCoefficient(int i, int k) const { return m_coefs[i][k]; }

    /// Human-readable equations of the fitted model.
    std::string GetEquations() const;

    /// Save the fitted model to a JSON file, and load one (which must have the same library degree and features).
    void Save(const std::string& filename) const;
    void Load(const std::string& filename);

    /// A model of the library degree and features saved in a JSON file, loaded from it.
    static std::shared_ptr<SCMSindyResidual> FromFile(const std::string& filename);

    static const char* GetFeatureName(int j);
    static const char* GetOutputName(int i);

    /// SCMTerrain::ContactForceCorrection
    virtual void Correct(ChBody& body,
                         const ChVector3d& force,
                         const ChVector3d& torque,
                         ChVector3d& dforce,
                         ChVector3d& dtorque) override;

    /// SCMTerrain::ContactForceCorrection: the wheels, corrected whether SCM touches them or not.
    virtual std::vector<ChBody*> GetBodies() const override;

  private:
    struct Wheel {
        double radius;
        ChVector3d axis;
        // Slip history (feature h), and trust as of the last correction; updated from const queries
        mutable double history = 0;
        mutable double history_time = -1;
        mutable double spread = 0;
        mutable double outside = 0;
        mutable bool corrected = false;
    };

    // A library term: a product of factors
    struct Factor {
        enum Kind { POW, SAT, RELU } kind;
        int feature;
        double param;  // exponent, saturation scale, or hinge offset
    };
    using Term = std::vector<Factor>;

    void BuildLibrary();
    Term ParseTerm(const std::string& text) const;
    static std::string TermName(const Term& term);

    int m_degree;
    unsigned int m_features;
    std::vector<Term> m_terms;  // the library
    std::array<std::vector<bool>, kNumOutputs> m_allowed;  // terms each output may use
    std::vector<Coefficients> m_ensemble;
    double m_history_time = 0.5;
    std::array<std::vector<double>, kNumOutputs> m_coefs;
    Features m_x_min, m_x_max;  // training range of the features
    Outputs m_y_max;            // largest output magnitude in training
    bool m_trained = false;

    const SCMTerrain* m_terrain = nullptr;
    std::unordered_map<const ChBody*, Wheel> m_wheels;
    std::vector<std::shared_ptr<ChBody>> m_wheel_bodies;  // keep the keys alive
    double m_contact_gap = 0.02;
    double m_force_scale = 100;
    double m_speed_eps = 0.02;
    double m_gain = 1;
    std::array<bool, kNumOutputs> m_applied = {true, true, true, true};
    std::vector<Sample> m_samples;
};

/// @} vehicle_terrain

}  // namespace vehicle
}  // namespace chrono

#endif
