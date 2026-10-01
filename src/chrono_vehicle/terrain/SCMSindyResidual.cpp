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

#include "chrono_vehicle/terrain/SCMSindyResidual.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include <Eigen/Dense>

#include "chrono/input_output/ChUtilsJSON.h"
#include "chrono/physics/ChSystem.h"
#include "chrono_vehicle/ChWorldFrame.h"

#include "chrono_thirdparty/rapidjson/document.h"

namespace chrono {
namespace vehicle {

namespace {
const char* kFeatureNames[SCMSindyResidual::kNumFeatures] = {"s", "b", "theta", "z", "d", "h"};
const char* kOutputNames[SCMSindyResidual::kNumOutputs] = {"dFl/F", "dFt/F", "dFn/F", "dTa/(F r)"};
}  // namespace

SCMSindyResidual::SCMSindyResidual(int degree, unsigned int features) : m_degree(degree), m_features(features & 0x3F) {
    if (degree < 1 || degree > 5)
        throw std::invalid_argument("SCMSindyResidual: library degree must be 1 to 5");
    if (m_features == 0)
        throw std::invalid_argument("SCMSindyResidual: no features");
    BuildLibrary();
    for (auto& c : m_coefs)
        c.assign(m_terms.size(), 0.0);
    for (auto& a : m_allowed)
        a.assign(m_terms.size(), true);
    m_x_min.fill(0);
    m_x_max.fill(0);
    m_y_max.fill(0);
}

const char* SCMSindyResidual::GetFeatureName(int j) {
    return kFeatureNames[j];
}

const char* SCMSindyResidual::GetOutputName(int i) {
    return kOutputNames[i];
}

void SCMSindyResidual::BuildLibrary() {
    // All monomials up to total degree m_degree, by degree, then lexicographically
    m_terms.clear();
    for (int d = 0; d <= m_degree; ++d) {
        std::array<int, kNumFeatures> e{};
        // Enumerate exponent vectors summing to d
        std::vector<std::array<int, kNumFeatures>> terms;
        auto recurse = [&](auto&& self, int j, int left) -> void {
            if (j == kNumFeatures - 1) {
                e[j] = left;
                terms.push_back(e);
                return;
            }
            for (int p = left; p >= 0; --p) {
                e[j] = p;
                self(self, j + 1, left - p);
            }
        };
        recurse(recurse, 0, d);
        for (const auto& term : terms) {
            bool used = true;
            for (int j = 0; j < kNumFeatures; ++j)
                used = used && (term[j] == 0 || (m_features >> j & 1));
            if (!used)
                continue;
            Term t;
            for (int j = 0; j < kNumFeatures; ++j)
                if (term[j] > 0)
                    t.push_back({Factor::POW, j, double(term[j])});
            m_terms.push_back(t);
        }
    }
}

std::string SCMSindyResidual::TermName(const Term& term) {
    std::ostringstream os;
    for (size_t f = 0; f < term.size(); ++f) {
        const Factor& fa = term[f];
        if (f)
            os << "*";
        const char* x = kFeatureNames[fa.feature];
        switch (fa.kind) {
            case Factor::POW:
                os << x;
                if (fa.param != 1)
                    os << "^" << fa.param;
                break;
            case Factor::SAT:
                os << "sat(" << x << "/" << fa.param << ")";
                break;
            case Factor::RELU:
                os << "relu(" << x << "-" << fa.param << ")";
                break;
        }
    }
    const std::string name = os.str();
    return name.empty() ? "1" : name;
}

std::string SCMSindyResidual::GetTermName(int k) const {
    return TermName(m_terms[k]);
}

SCMSindyResidual::Term SCMSindyResidual::ParseTerm(const std::string& text) const {
    auto fail = [&](const std::string& why) { throw std::invalid_argument("SCMSindyResidual: term \"" + text + "\": " + why); };
    auto feature = [&](const std::string& name) {
        for (int j = 0; j < kNumFeatures; ++j)
            if (name == kFeatureNames[j]) {
                if (!(m_features >> j & 1))
                    fail("feature " + name + " is not among the model's");
                return j;
            }
        fail("unknown feature " + name);
        return -1;
    };
    auto number = [&](const std::string& v) {
        size_t used = 0;
        double x = 0;
        try {
            x = std::stod(v, &used);
        } catch (const std::exception&) {
            used = 0;
        }
        if (used != v.size() || v.empty())
            fail("bad number " + v);
        return x;
    };
    Term term;
    if (text == "1")
        return term;
    std::istringstream is(text);
    std::string f;
    while (std::getline(is, f, '*')) {
        if (f.rfind("sat(", 0) == 0 && f.back() == ')') {
            const std::string in = f.substr(4, f.size() - 5);
            const size_t slash = in.find('/');
            if (slash == std::string::npos)
                fail("sat needs x/c");
            const double c = number(in.substr(slash + 1));
            if (!(c > 0))
                fail("sat scale must be positive");
            term.push_back({Factor::SAT, feature(in.substr(0, slash)), c});
        } else if (f.rfind("relu(", 0) == 0 && f.back() == ')') {
            const std::string in = f.substr(5, f.size() - 6);
            const size_t minus = in.find('-');
            if (minus == std::string::npos)
                fail("relu needs x-a");
            term.push_back({Factor::RELU, feature(in.substr(0, minus)), number(in.substr(minus + 1))});
        } else {
            const size_t caret = f.find('^');
            const double p = caret == std::string::npos ? 1.0 : number(f.substr(caret + 1));
            if (!(p > 0))
                fail("exponent must be positive");
            term.push_back({Factor::POW, feature(f.substr(0, caret)), p});
        }
    }
    if (term.empty())
        fail("empty");
    return term;
}

void SCMSindyResidual::AddLibraryTerm(const std::string& text) {
    Term term = ParseTerm(text);
    const std::string name = TermName(term);
    for (const auto& t : m_terms)
        if (TermName(t) == name)
            throw std::invalid_argument("SCMSindyResidual: term " + name + " is already in the library");
    m_terms.push_back(term);
    for (auto& c : m_coefs)
        c.push_back(0.0);
    for (auto& a : m_allowed)
        a.push_back(true);
    m_trained = false;
}

void SCMSindyResidual::SetOutputTerms(int i, const std::function<bool(const std::string&)>& allow) {
    for (size_t k = 0; k < m_terms.size(); ++k)
        m_allowed[i][k] = allow(TermName(m_terms[k]));
    m_trained = false;
}

SCMSindyResidual::Outputs SCMSindyResidual::EvaluateSpread(const Features& x) const {
    Outputs spread{};
    if (m_ensemble.size() < 2)
        return spread;
    Features xc;
    for (int j = 0; j < kNumFeatures; ++j)
        xc[j] = std::clamp(x[j], m_x_min[j], m_x_max[j]);
    std::vector<double> theta;
    EvaluateLibrary(xc, theta);
    for (int i = 0; i < kNumOutputs; ++i) {
        double sum = 0, sum2 = 0;
        for (const auto& m : m_ensemble) {
            double v = 0;
            for (size_t k = 0; k < theta.size() && k < m[i].size(); ++k)
                v += m[i][k] * theta[k];
            v = std::clamp(v, -m_y_max[i], m_y_max[i]);
            sum += v;
            sum2 += v * v;
        }
        const double n = double(m_ensemble.size());
        spread[i] = std::sqrt(std::max(0.0, sum2 / n - (sum / n) * (sum / n)));
    }
    return spread;
}

bool SCMSindyResidual::GetWheelTrust(const ChBody& wheel, double& spread, double& outside) const {
    auto itr = m_wheels.find(&wheel);
    if (itr == m_wheels.end() || !itr->second.corrected)
        return false;
    spread = itr->second.spread;
    outside = itr->second.outside;
    return true;
}

void SCMSindyResidual::EvaluateLibrary(const Features& x, std::vector<double>& theta) const {
    theta.resize(m_terms.size());
    for (size_t k = 0; k < m_terms.size(); ++k) {
        double v = 1;
        for (const Factor& f : m_terms[k]) {
            const double u = x[f.feature];
            switch (f.kind) {
                case Factor::POW: {
                    const int n = static_cast<int>(f.param);
                    if (n == f.param) {
                        double p = 1;
                        for (int i = 0; i < n; ++i)
                            p *= u;
                        v *= p;
                    } else {
                        v *= std::pow(std::max(u, 0.0), f.param);
                    }
                    break;
                }
                case Factor::SAT:
                    v *= std::copysign(1 - std::exp(-std::abs(u) / f.param), u);
                    break;
                case Factor::RELU:
                    v *= std::max(u - f.param, 0.0);
                    break;
            }
        }
        theta[k] = v;
    }
}

SCMSindyResidual::Outputs SCMSindyResidual::Evaluate(const Features& x) const {
    Features xc;
    for (int j = 0; j < kNumFeatures; ++j)
        xc[j] = std::clamp(x[j], m_x_min[j], m_x_max[j]);
    std::vector<double> theta;
    EvaluateLibrary(xc, theta);
    Outputs y;
    for (int i = 0; i < kNumOutputs; ++i) {
        double v = 0;
        for (size_t k = 0; k < theta.size(); ++k)
            v += m_coefs[i][k] * theta[k];
        y[i] = std::clamp(v, -m_y_max[i], m_y_max[i]);
    }
    return y;
}

void SCMSindyResidual::AddWheel(std::shared_ptr<ChBody> wheel, double radius, const ChVector3d& axis) {
    if (!wheel || !(radius > 0) || axis.Length() < 1e-9)
        throw std::invalid_argument("SCMSindyResidual::AddWheel: null wheel, bad radius or axis");
    m_wheels[wheel.get()] = {radius, axis.GetNormalized()};
    m_wheel_bodies.push_back(std::move(wheel));
}

bool SCMSindyResidual::GetState(const ChBody& wheel, const ChVector3d& scm_force, State& state) const {
    auto itr = m_wheels.find(&wheel);
    if (itr == m_wheels.end() || !m_terrain)
        return false;
    const double r = itr->second.radius;
    const ChFrameMoving<>& frame = wheel.GetFrameRefToAbs();
    const ChVector3d& center = frame.GetPos();
    const ChVector3d& vertical = ChWorldFrame::Vertical();
    // Up is against gravity, if there is any
    ChVector3d up = vertical;
    if (const ChSystem* sys = wheel.GetSystem()) {
        const ChVector3d& g = sys->GetGravitationalAcceleration();
        if (g.Length() > 1e-9)
            up = -g.GetNormalized();
    }

    // Wheel-ground frame
    const ChVector3d n = m_terrain->GetInitNormal(center);
    ChVector3d a = frame.GetRotMat() * itr->second.axis;
    ChVector3d l = Vcross(a, n);
    if (l.Length() < 1e-6)
        return false;
    l.Normalize();

    const double load = scm_force.Dot(n);

    // Slips: rim speed against the center's speed along and across the rolling direction. The axle is taken the way the
    // wheel spins (a wheel mounted mirrored, or driving in reverse, is the same wheel), or, barely spinning, the way it
    // travels. Not by the sum of the two: a wheel spinning forward while sliding back faster would flip to a skid.
    const ChVector3d v = frame.GetPosDt();
    double rim = wheel.GetAngVelParent().Dot(a) * r;
    double vl = v.Dot(l);
    if (std::abs(rim) > m_speed_eps ? rim < 0 : vl < 0) {
        a = -a;
        l = -l;
        rim = -rim;
        vl = -vl;
    }
    const ChVector3d t = Vcross(n, l);
    const double den = std::max({std::abs(rim), std::abs(vl), m_speed_eps});

    state.x[0] = (rim - vl) / den;
    state.x[1] = std::clamp(v.Dot(t) / den, -1.0, 1.0);
    // Slip history: slip low-pass filtered in the system's time, updated once per time
    const Wheel& w = itr->second;
    const double now = wheel.GetSystem() ? wheel.GetSystem()->GetChTime() : 0.0;
    if (w.history_time < 0) {
        w.history = state.x[0];
        w.history_time = now;
    } else if (now > w.history_time) {
        w.history += (1 - std::exp(-(now - w.history_time) / m_history_time)) * (state.x[0] - w.history);
        w.history_time = now;
    }
    state.x[5] = w.history;
    state.x[2] = std::asin(std::clamp(l.Dot(up), -1.0, 1.0));
    // Sinkage: how far the wheel's lowest point is below the undeformed ground, along its normal
    const double clearance = (ChWorldFrame::Height(center) - m_terrain->GetInitHeight(center)) * n.Dot(vertical);
    if (clearance - r > m_contact_gap)
        return false;
    state.x[3] = (r - clearance) / r;
    // Rut ahead: how far SCM's ground a radius ahead lies below the undeformed ground
    const ChVector3d ahead = center + r * l;
    state.x[4] = (m_terrain->GetInitHeight(ahead) - m_terrain->GetHeight(ahead)) * n.Dot(vertical) / r;
    state.load = load;
    state.radius = r;
    state.l = l;
    state.t = t;
    state.n = n;
    state.a = a;
    return true;
}

bool SCMSindyResidual::MakeSample(const ChBody& wheel,
                                  const ChVector3d& scm_force,
                                  const ChVector3d& scm_torque,
                                  const ChVector3d& ref_force,
                                  const ChVector3d& ref_torque,
                                  Sample& sample) const {
    State s;
    if (!GetState(wheel, scm_force, s))
        return false;
    const ChVector3d df = ref_force - scm_force;
    const ChVector3d dt = ref_torque - scm_torque;
    sample.x = s.x;
    const double F = m_force_scale;
    sample.y = {df.Dot(s.l) / F, df.Dot(s.t) / F, df.Dot(s.n) / F, dt.Dot(s.a) / (F * s.radius)};
    sample.load = s.load;
    return true;
}

bool SCMSindyResidual::AddSample(const ChBody& wheel,
                                 const ChVector3d& scm_force,
                                 const ChVector3d& scm_torque,
                                 const ChVector3d& ref_force,
                                 const ChVector3d& ref_torque) {
    Sample sample;
    if (!MakeSample(wheel, scm_force, scm_torque, ref_force, ref_torque, sample))
        return false;
    m_samples.push_back(sample);
    return true;
}

void SCMSindyResidual::WriteSamples(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error("SCMSindyResidual: cannot write " + filename);
    out << "s,b,theta,z,d,h,load,dFl_F,dFt_F,dFn_F,dTa_Fr\n";
    out << std::setprecision(9);
    for (const auto& s : m_samples) {
        for (double x : s.x)
            out << x << ",";
        out << s.load;
        for (double y : s.y)
            out << "," << y;
        out << "\n";
    }
}

std::array<SCMSindyResidual::FitReport, SCMSindyResidual::kNumOutputs> SCMSindyResidual::Fit(const FitSettings& settings) {
    const int N = static_cast<int>(m_samples.size());
    const int K = GetNumTerms();
    if (N < 2 * K)
        throw std::runtime_error("SCMSindyResidual::Fit: " + std::to_string(N) + " samples for " + std::to_string(K) +
                                 " library terms; need at least twice as many");

    // Training range
    for (int j = 0; j < kNumFeatures; ++j) {
        m_x_min[j] = m_x_max[j] = m_samples[0].x[j];
        for (const auto& s : m_samples) {
            m_x_min[j] = std::min(m_x_min[j], s.x[j]);
            m_x_max[j] = std::max(m_x_max[j], s.x[j]);
        }
    }
    for (int i = 0; i < kNumOutputs; ++i) {
        m_y_max[i] = 0;
        for (const auto& s : m_samples)
            m_y_max[i] = std::max(m_y_max[i], std::abs(s.y[i]));
    }

    // Weights, normalized to a mean of one; rows of the least-squares problem scaled by their square roots
    Eigen::VectorXd sw(N);
    double wsum = 0;
    for (int r = 0; r < N; ++r)
        wsum += m_samples[r].weight;
    for (int r = 0; r < N; ++r)
        sw[r] = std::sqrt(m_samples[r].weight * N / wsum);

    // Library matrix, columns normalized to unit (weighted) RMS so the threshold is in output units for every term
    Eigen::MatrixXd A(N, K);
    std::vector<double> theta;
    for (int r = 0; r < N; ++r) {
        EvaluateLibrary(m_samples[r].x, theta);
        for (int k = 0; k < K; ++k)
            A(r, k) = theta[k] * sw[r];
    }
    Eigen::VectorXd scale(K);
    for (int k = 0; k < K; ++k) {
        scale[k] = std::sqrt(A.col(k).squaredNorm() / N);
        if (scale[k] > 1e-12)
            A.col(k) /= scale[k];
    }

    std::array<FitReport, kNumOutputs> report;
    for (int i = 0; i < kNumOutputs; ++i) {
        Eigen::VectorXd y(N);
        for (int r = 0; r < N; ++r)
            y[r] = m_samples[r].y[i] * sw[r];

        // Sequentially thresholded least squares, over the terms this output may use
        std::vector<bool> active(K);
        for (int k = 0; k < K; ++k)
            active[k] = scale[k] > 1e-12 && m_allowed[i][k];
        Eigen::VectorXd xi = Eigen::VectorXd::Zero(K);
        for (int it = 0; it < settings.max_iterations; ++it) {
            std::vector<int> cols;
            for (int k = 0; k < K; ++k)
                if (active[k])
                    cols.push_back(k);
            xi.setZero();
            if (cols.empty())
                break;
            Eigen::MatrixXd As(N, cols.size());
            for (size_t c = 0; c < cols.size(); ++c)
                As.col(c) = A.col(cols[c]);
            Eigen::MatrixXd G = As.transpose() * As / N;
            G.diagonal().array() += settings.ridge;
            const Eigen::VectorXd sol = G.ldlt().solve(As.transpose() * y / N);
            bool changed = false;
            for (size_t c = 0; c < cols.size(); ++c) {
                if (std::abs(sol[c]) < settings.threshold) {
                    active[cols[c]] = false;
                    changed = true;
                } else {
                    xi[cols[c]] = sol[c];
                }
            }
            if (!changed)
                break;
        }

        // Unscaled coefficients and fit quality
        const Eigen::VectorXd pred = A * xi;
        const double mean = y.dot(sw) / sw.squaredNorm();
        const double ss_res = (y - pred).squaredNorm();
        const double ss_tot = (y - mean * sw).squaredNorm();
        report[i].num_terms = 0;
        for (int k = 0; k < K; ++k) {
            m_coefs[i][k] = xi[k] != 0 ? xi[k] / scale[k] : 0;
            report[i].num_terms += xi[k] != 0;
        }
        report[i].r2 = ss_tot > 0 ? 1 - ss_res / ss_tot : 1;
        report[i].rmse = std::sqrt(ss_res / N);
    }
    m_trained = true;
    return report;
}

std::string SCMSindyResidual::GetEquations() const {
    std::ostringstream os;
    os << std::showpos << std::setprecision(4);
    for (int i = 0; i < kNumOutputs; ++i) {
        os << kOutputNames[i] << " =";
        bool any = false;
        for (int k = 0; k < GetNumTerms(); ++k) {
            if (m_coefs[i][k] == 0)
                continue;
            os << " " << m_coefs[i][k] << " " << GetTermName(k);
            any = true;
        }
        if (!any)
            os << " 0";
        os << "\n";
    }
    return os.str();
}

void SCMSindyResidual::Save(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error("SCMSindyResidual: cannot write " + filename);
    auto array = [&](const auto& v) {
        out << "[";
        for (size_t k = 0; k < v.size(); ++k)
            out << (k ? ", " : "") << v[k];
        out << "]";
    };
    out << std::setprecision(17);
    out << "{\n  \"Type\": \"SCMSindyResidual\",\n  \"Degree\": " << m_degree << ",\n";
    out << "  \"Features\": [";
    bool first = true;
    for (int j = 0; j < kNumFeatures; ++j) {
        if (m_features >> j & 1) {
            out << (first ? "" : ", ") << "\"" << kFeatureNames[j] << "\"";
            first = false;
        }
    }
    out << "],\n  \"Outputs\": [";
    for (int i = 0; i < kNumOutputs; ++i)
        out << (i ? ", " : "") << "\"" << kOutputNames[i] << "\"";
    out << "],\n  \"Terms\": [";
    for (int k = 0; k < GetNumTerms(); ++k)
        out << (k ? ", " : "") << "\"" << GetTermName(k) << "\"";
    out << "],\n  \"Feature Min\": ";
    array(m_x_min);
    out << ",\n  \"Feature Max\": ";
    array(m_x_max);
    out << ",\n  \"Output Max\": ";
    array(m_y_max);
    out << ",\n  \"Applied Outputs\": [";
    for (int i = 0; i < kNumOutputs; ++i)
        out << (i ? ", " : "") << (m_applied[i] ? "true" : "false");
    out << "]";
    bool restricted = false;
    for (const auto& a : m_allowed)
        for (bool b : a)
            restricted = restricted || !b;
    if (restricted) {
        out << ",\n  \"Allowed Terms\": [\n";
        for (int i = 0; i < kNumOutputs; ++i) {
            out << "    [";
            for (size_t k = 0; k < m_allowed[i].size(); ++k)
                out << (k ? ", " : "") << (m_allowed[i][k] ? 1 : 0);
            out << (i + 1 < kNumOutputs ? "],\n" : "]\n");
        }
        out << "  ]";
    }
    if (!m_ensemble.empty()) {
        out << ",\n  \"Ensemble\": [\n";
        for (size_t e = 0; e < m_ensemble.size(); ++e) {
            out << "    [";
            for (int i = 0; i < kNumOutputs; ++i) {
                out << (i ? ", " : "");
                array(m_ensemble[e][i]);
            }
            out << (e + 1 < m_ensemble.size() ? "],\n" : "]\n");
        }
        out << "  ]";
    }
    out << ",\n  \"History Time\": " << m_history_time;
    out << ",\n  \"Force Scale\": " << m_force_scale << ",\n  \"Contact Gap\": " << m_contact_gap
        << ",\n  \"Speed Eps\": " << m_speed_eps << ",\n  \"Coefficients\": [\n";
    for (int i = 0; i < kNumOutputs; ++i) {
        out << "    ";
        array(m_coefs[i]);
        out << (i + 1 < kNumOutputs ? ",\n" : "\n");
    }
    out << "  ]\n}\n";
}

std::shared_ptr<SCMSindyResidual> SCMSindyResidual::FromFile(const std::string& filename) {
    rapidjson::Document d;
    ReadFileJSON(filename, d);
    if (!d.IsObject() || !d.HasMember("Degree") || !d.HasMember("Features") || !d["Features"].IsArray())
        throw std::runtime_error("SCMSindyResidual: " + filename + ": not an SCMSindyResidual model");
    unsigned int features = 0;
    for (const auto& f : d["Features"].GetArray())
        for (int j = 0; j < kNumFeatures; ++j)
            if (std::string(f.GetString()) == kFeatureNames[j])
                features |= 1u << j;
    auto model = std::make_shared<SCMSindyResidual>(d["Degree"].GetInt(), features);
    model->Load(filename);
    return model;
}

void SCMSindyResidual::Load(const std::string& filename) {
    rapidjson::Document d;
    ReadFileJSON(filename, d);
    auto fail = [&](const std::string& what) { throw std::runtime_error("SCMSindyResidual: " + filename + ": " + what); };
    if (!d.IsObject() || !d.HasMember("Type") || std::string(d["Type"].GetString()) != "SCMSindyResidual")
        fail("not an SCMSindyResidual model");
    if (!d.HasMember("Degree") || d["Degree"].GetInt() != m_degree)
        fail("library degree differs from the model's (" + std::to_string(m_degree) + ")");
    unsigned int features = 0;
    if (d.HasMember("Features") && d["Features"].IsArray()) {
        for (const auto& f : d["Features"].GetArray())
            for (int j = 0; j < kNumFeatures; ++j)
                if (std::string(f.GetString()) == kFeatureNames[j])
                    features |= 1u << j;
    }
    if (features != m_features)
        fail("library features differ from the model's");
    auto read = [&](const char* name, auto& v, size_t n, size_t n_min) {
        if (!d.HasMember(name) || !d[name].IsArray() || d[name].Size() > n || d[name].Size() < n_min)
            fail(std::string("bad or missing \"") + name + "\"");
        std::fill(v.begin(), v.end(), 0.0);
        for (rapidjson::SizeType k = 0; k < d[name].Size(); ++k)
            v[k] = d[name][k].GetDouble();
    };
    // Models saved before the history feature have one feature fewer (and do not use it)
    read("Feature Min", m_x_min, kNumFeatures, kNumFeatures - 1);
    read("Feature Max", m_x_max, kNumFeatures, kNumFeatures - 1);
    read("Output Max", m_y_max, kNumOutputs, kNumOutputs);
    // The library, as saved: monomials and any terms added to them
    if (!d.HasMember("Terms") || !d["Terms"].IsArray())
        fail("bad or missing \"Terms\"");
    std::vector<Term> terms;
    try {
        for (const auto& t : d["Terms"].GetArray())
            terms.push_back(ParseTerm(t.GetString()));
    } catch (const std::invalid_argument& e) {
        fail(e.what());
    }
    m_terms = terms;
    if (!d.HasMember("Coefficients") || !d["Coefficients"].IsArray() || d["Coefficients"].Size() != kNumOutputs)
        fail("bad or missing \"Coefficients\"");
    for (int i = 0; i < kNumOutputs; ++i) {
        const auto& row = d["Coefficients"][i];
        m_coefs[i].assign(m_terms.size(), 0.0);
        if (!row.IsArray() || row.Size() != m_terms.size())
            fail("coefficient row " + std::to_string(i) + " does not match the library");
        for (rapidjson::SizeType k = 0; k < row.Size(); ++k)
            m_coefs[i][k] = row[k].GetDouble();
    }
    for (auto& a : m_allowed)
        a.assign(m_terms.size(), true);
    if (d.HasMember("Allowed Terms") && d["Allowed Terms"].IsArray() && d["Allowed Terms"].Size() == kNumOutputs) {
        for (int i = 0; i < kNumOutputs; ++i) {
            const auto& row = d["Allowed Terms"][i];
            if (!row.IsArray() || row.Size() != m_terms.size())
                fail("\"Allowed Terms\" row " + std::to_string(i) + " does not match the library");
            for (rapidjson::SizeType k = 0; k < row.Size(); ++k)
                m_allowed[i][k] = row[k].GetInt() != 0;
        }
    }
    m_ensemble.clear();
    if (d.HasMember("Ensemble") && d["Ensemble"].IsArray()) {
        for (const auto& mem : d["Ensemble"].GetArray()) {
            if (!mem.IsArray() || mem.Size() != kNumOutputs)
                fail("bad \"Ensemble\" member");
            Coefficients c;
            for (int i = 0; i < kNumOutputs; ++i) {
                if (!mem[i].IsArray() || mem[i].Size() != m_terms.size())
                    fail("\"Ensemble\" member does not match the library");
                for (const auto& v : mem[i].GetArray())
                    c[i].push_back(v.GetDouble());
            }
            m_ensemble.push_back(c);
        }
    }
    if (d.HasMember("History Time"))
        m_history_time = d["History Time"].GetDouble();
    if (!d.HasMember("Force Scale"))
        fail("missing \"Force Scale\" (a model saved before outputs were scaled by a constant force)");
    m_force_scale = d["Force Scale"].GetDouble();
    if (d.HasMember("Applied Outputs") && d["Applied Outputs"].IsArray() && d["Applied Outputs"].Size() == kNumOutputs) {
        for (rapidjson::SizeType i = 0; i < kNumOutputs; ++i)
            m_applied[i] = d["Applied Outputs"][i].GetBool();
    }
    if (d.HasMember("Contact Gap"))
        m_contact_gap = d["Contact Gap"].GetDouble();
    if (d.HasMember("Speed Eps"))
        m_speed_eps = d["Speed Eps"].GetDouble();
    m_trained = true;
}

void SCMSindyResidual::Correct(ChBody& body, const ChVector3d& force, const ChVector3d& torque, ChVector3d& dforce, ChVector3d& dtorque) {
    dforce = VNULL;
    dtorque = VNULL;
    State s;
    if (!m_trained || m_gain == 0 || !GetState(body, force, s))
        return;
    Outputs y = Evaluate(s.x);
    for (int i = 0; i < kNumOutputs; ++i)
        if (!m_applied[i])
            y[i] = 0;

    // Trust: ensemble spread over the applied outputs, and how far outside the training range the state lies
    const Wheel& w_rec = m_wheels.find(&body)->second;
    const Outputs spread = EvaluateSpread(s.x);
    w_rec.spread = 0;
    for (int i = 0; i < kNumOutputs; ++i)
        if (m_applied[i])
            w_rec.spread = std::max(w_rec.spread, spread[i]);
    w_rec.outside = 0;
    for (int j = 0; j < kNumFeatures; ++j) {
        if (!(m_features >> j & 1))
            continue;
        const double range = std::max(m_x_max[j] - m_x_min[j], 1e-9);
        w_rec.outside = std::max(w_rec.outside, std::max(s.x[j] - m_x_max[j], m_x_min[j] - s.x[j]) / range);
    }
    w_rec.corrected = true;

    // Full at the undeformed ground and below, none a contact gap above it
    const double ramp = m_contact_gap > 0 ? std::clamp(1 + s.x[3] * s.radius / m_contact_gap, 0.0, 1.0) : 1.0;
    const double w = m_gain * ramp * m_force_scale;
    dforce = w * (y[0] * s.l + y[1] * s.t + y[2] * s.n);
    dtorque = (w * s.radius * y[3]) * s.a;
}

std::vector<ChBody*> SCMSindyResidual::GetBodies() const {
    std::vector<ChBody*> bodies;
    for (const auto& b : m_wheel_bodies)
        bodies.push_back(b.get());
    return bodies;
}

}  // namespace vehicle
}  // namespace chrono
