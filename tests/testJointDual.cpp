// Standalone checks of the joint dual objectives (H, H_alpha, F): gradient
// and Hessian vs. finite differences, and the Hessian null-spaces of
// latex/contributions.tex (eq:joint-nullspace).

#include "core/domain.h"
#include "core/jointDual.h"
#include "core/shiftedDiagram.h"
#include "testUtils.h"

#include <Eigen/Dense>

#include <cmath>
#include <string>
#include <vector>

namespace
{

using PP::Domain;
using PP::JointDual;
using PP::ShiftedDiagram;
using PPTest::Check;
using Vector = JointDual::Vector;

std::string
_GetDualName(const JointDual& dual)
{
    if (not dual.HasShifts()) {
        return "F";
    }
    if (not dual.HasWeights()) {
        return "G";
    }
    return std::isfinite(dual.GetAlpha()) ? "H_alpha" : "H";
}

double
_Evaluate(const JointDual& dual, const Vector& z)
{
    ShiftedDiagram d;
    std::string err;
    dual.Build(z, d, err);
    return dual.GetValue(z, d);
}

Vector
_EvaluateGradient(const JointDual& dual, const Vector& z)
{
    ShiftedDiagram d;
    std::string err;
    dual.Build(z, d, err);
    Vector g;
    dual.GetGradient(z, d, g);
    return g;
}

// -Hessian at z, dense, without regularization.
Eigen::MatrixXd
_GetNegatedHessian(const JointDual& dual, const Vector& z)
{
    ShiftedDiagram d;
    std::string err;
    dual.Build(z, d, err);
    JointDual::SparseMatrix N;
    dual.GetNegatedHessian(d, 0.0, N);
    return Eigen::MatrixXd(N);
}

// Small random multipliers (shifts and weights).
Vector
_MakeRandomMultipliers(const JointDual& dual, double scale, unsigned seed)
{
    const std::vector<double> r = PPTest::RandomValues(dual.GetNumVariables(), -scale, scale, seed);
    return Eigen::Map<const Vector>(r.data(), Eigen::Index(r.size()));
}

void
_TestDerivatives(const JointDual& dual, const std::string& name)
{
    const std::string what = name + " " + _GetDualName(dual);
    const Vector z = _MakeRandomMultipliers(dual, 0.002, 7);

    ShiftedDiagram d;
    std::string err;
    const bool ok = dual.Build(z, d, err);
    Check(ok and d.GetNumEmptyCells() == 0, what + ": power diagram without empty cells");
    if (not ok) {
        return;
    }
    Vector g;
    dual.GetGradient(z, d, g);
    const Eigen::MatrixXd N = _GetNegatedHessian(dual, z);

    const double h = 1e-6;
    double gradErr = 0.0, hessErr = 0.0;
    for (int k = 0; k < dual.GetNumVariables(); ++k) {
        Vector zp = z, zm = z;
        zp[k] += h;
        zm[k] -= h;
        gradErr = std::max(gradErr, std::fabs((_Evaluate(dual, zp) - _Evaluate(dual, zm)) / (2 * h) - g[k]));
        const Vector col = (_EvaluateGradient(dual, zp) - _EvaluateGradient(dual, zm)) / (2 * h);
        hessErr = std::max(hessErr, (col + N.col(k)).lpNorm<Eigen::Infinity>());
    }
    const double gradScale = std::max(1e-3, g.lpNorm<Eigen::Infinity>());
    const double hessScale = N.cwiseAbs().maxCoeff();
    Check(gradErr < 1e-5 * gradScale,
          what + ": gradient matches finite differences of the value (rel err " +
              PPTest::FormatValue(gradErr / gradScale) + ")");
    Check(hessErr < 1e-4 * hessScale,
          what + ": Hessian matches finite differences of the gradient (rel err " +
              PPTest::FormatValue(hessErr / hessScale) + ")");
    Check((N - N.transpose()).cwiseAbs().maxCoeff() < 1e-12 * hessScale, what + ": Hessian is symmetric");
}

// Null-space modes of the Hessian of H at z (eq:joint-nullspace):
//   du_i = mu p_i + t,   dw_i = mu gamma_i + 2 t^T cbar_i + c,
// gamma_i = |cbar_i|^2 - (|cbar_i - q_i|^2 - w_i) (eq:pivot), with the
// generators p_i and the original weights w_i of z.
std::vector<Vector>
_MakeNullModes(const JointDual& dual, const Vector& z, const ShiftedDiagram& d)
{
    const int dim = dual.GetDimension();
    const int s = dual.GetStride();
    const int n = dual.GetNumSites();
    const std::vector<double>& q = dual.GetSites();
    const std::vector<double>& cbar = dual.GetTargets();
    std::vector<Vector> modes(dim + 2, Vector::Zero(dual.GetNumVariables()));
    for (int i = 0; i < n; ++i) {
        double gamma = z[s * i + dim];
        for (int c = 0; c < dim; ++c) {
            const double cb = cbar[dim * i + c];
            gamma += cb * cb - (cb - q[dim * i + c]) * (cb - q[dim * i + c]);
        }
        for (int c = 0; c < dim; ++c) {
            modes[0][s * i + c] = d.GetGenerator(i)[c];       // dilation mu
            modes[1 + c][s * i + c] = 1.0;                    // translation t_c
            modes[1 + c][s * i + dim] = 2.0 * cbar[dim * i + c];
        }
        modes[0][s * i + dim] = gamma;
        modes[dim + 1][s * i + dim] = 1.0;                    // constant c
    }
    return modes;
}

void
_TestNullSpace(const Domain& D, const std::string& name, int n)
{
    const int dim = D.GetDimension();
    const std::vector<double> sites = PPTest::RandomSites(n, dim, 0.05, 0.95, 3);
    const std::vector<double> targets = PPTest::RandomSites(n, dim, 0.05, 0.95, 4);
    const std::vector<double> volumes(n, D.GetVolume() / n);

    const JointDual H(D, sites, targets, volumes, HUGE_VAL, JointDual::Unknowns::ShiftsAndWeights);
    const Vector z = _MakeRandomMultipliers(H, 0.01, 9);
    ShiftedDiagram d;
    std::string err;
    H.Build(z, d, err);
    const Eigen::MatrixXd N = _GetNegatedHessian(H, z);
    const double scale = N.cwiseAbs().maxCoeff();

    double modeErr = 0.0;
    for (const Vector& m : _MakeNullModes(H, z, d)) {
        modeErr = std::max(modeErr, (N * m).lpNorm<Eigen::Infinity>() / m.lpNorm<Eigen::Infinity>());
    }
    Check(modeErr < 1e-10 * scale, name + " H: the d + 2 homothety modes are in the Hessian null-space (rel " +
                                       PPTest::FormatValue(modeErr / scale) + ")");

    // N is positive semi-definite with exactly d + 2 zero eigenvalues.
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(N);
    const Eigen::VectorXd ev = eig.eigenvalues();
    int zeros = 0;
    for (int k = 0; k < ev.size(); ++k) {
        zeros += std::fabs(ev[k]) < 1e-9 * scale;
    }
    Check(ev[0] > -1e-9 * scale and zeros == dim + 2,
          name + " H: -Hessian is PSD with d + 2 = " + std::to_string(dim + 2) + " zero eigenvalues (found " +
              std::to_string(zeros) + ")");

    // H_alpha keeps only the constant on the weights.
    const JointDual Ha(D, sites, targets, volumes, 0.5, JointDual::Unknowns::ShiftsAndWeights);
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigA(_GetNegatedHessian(Ha, z));
    const Eigen::VectorXd evA = eigA.eigenvalues();
    int zerosA = 0;
    for (int k = 0; k < evA.size(); ++k) {
        zerosA += std::fabs(evA[k]) < 1e-9 * scale;
    }
    Check(evA[0] > -1e-9 * scale and zerosA == 1, name + " H_alpha: -Hessian is PSD with a single zero eigenvalue");
}

// The shifted sites reproduce the power distance identity eq:complete-square:
// |x - q_i|^2 - w_i - 2 u_i^T (x - cbar_i) = |x - p_i|^2 - lambda_i.
void
_TestShiftedSites(const Domain& D, const std::string& name)
{
    const int dim = D.GetDimension();
    const int n = 5;
    const std::vector<double> q = PPTest::RandomSites(n, dim, 0.1, 0.9, 1);
    const std::vector<double> cbar = PPTest::RandomSites(n, dim, 0.1, 0.9, 2);
    const std::vector<double> u = PPTest::RandomSites(n, dim, -0.05, 0.05, 3);
    const std::vector<double> w = PPTest::RandomValues(n, -0.01, 0.01, 4);
    ShiftedDiagram d;
    std::string err;
    d.Build(D, q, cbar, u, w, err);
    const std::vector<double> x = PPTest::RandomSites(1, dim, 0.0, 1.0, 5);
    double maxErr = 0.0;
    for (int i = 0; i < n; ++i) {
        double lhs = -w[i], rhs = -d.GetShiftedWeight(i);
        for (int c = 0; c < dim; ++c) {
            const double k = dim * i + c;
            lhs += (x[c] - q[k]) * (x[c] - q[k]) - 2.0 * u[k] * (x[c] - cbar[k]);
            rhs += (x[c] - d.GetGenerator(i)[c]) * (x[c] - d.GetGenerator(i)[c]);
        }
        maxErr = std::max(maxErr, std::fabs(lhs - rhs));
    }
    Check(maxErr < 1e-14, name + ": shifted sites satisfy |x-q|^2 - w - 2u.(x-cbar) = |x-p|^2 - lambda");
}

void
_TestDomain(const Domain& D, const std::string& name, int n)
{
    const int dim = D.GetDimension();
    const std::vector<double> sites = PPTest::RandomSites(n, dim, 0.05, 0.95, 11);
    const std::vector<double> targets = PPTest::RandomSites(n, dim, 0.05, 0.95, 12);
    const std::vector<double> volumes = PPTest::ScaleFractions(PPTest::RandomFractions(n, 13), D.GetVolume());

    _TestShiftedSites(D, name);
    _TestDerivatives(JointDual(D, sites, targets, volumes, HUGE_VAL, JointDual::Unknowns::ShiftsAndWeights), name);
    _TestDerivatives(JointDual(D, sites, targets, volumes, 0.3, JointDual::Unknowns::ShiftsAndWeights), name);
    _TestDerivatives(JointDual(D, sites, targets, volumes, HUGE_VAL, JointDual::Unknowns::Weights), name);
    JointDual G(D, sites, targets, volumes, HUGE_VAL, JointDual::Unknowns::Shifts);
    G.SetFixedWeights(PPTest::RandomValues(n, -0.002, 0.002, 14));
    _TestDerivatives(G, name);
    _TestNullSpace(D, name, n);
}

} // namespace

int
main()
{
    Domain square, L, cube;
    PPTest::MakeUnitSquare(square);
    PPTest::MakeLShape(L);
    PPTest::MakeUnitCube(cube);

    _TestDomain(square, "2D square", 12);
    _TestDomain(L, "2D L-shape", 12);
    _TestDomain(cube, "3D cube", 10);
    return PPTest::Report();
}
