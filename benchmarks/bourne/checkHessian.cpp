// Compares the assembled Hessian of H / H_alpha with central finite
// differences of the gradient at given multipliers, to test whether the
// Hessian is still accurate near degenerate facets (benchmarks/bourne).
//
//   checkHessian <case.txt> <z.txt> [h ...]
//   checkHessian <case.txt> <z.txt> eig [k]
//
// The second form prints the k smallest eigenvalues of the reduced Hessian
// of H in the shifts (the weights solved out: the Schur complement
// S = N_uu - N_uw N_ww^-1 N_wu, with eta = 0), after removing the gauge
// modes (translations and the dilation u_i = p_i - c), normalized by the
// mass scale 2 Vbar_i, with a smoothness ratio per eigenvector: the mean
// |du_i - du_j|^2 over facet neighbors divided by the mean |du_i|^2 (small:
// a smooth displacement field).
//
// case.txt is a benchBourne case file (dim, box, n, solver, alpha, volumes,
// targets, sites); z.txt holds the multipliers z (stride d + 1 per site:
// the shift u, then w). For each step h, prints the relative error per block
// (shift-shift, shift-weight, weight-weight), the worst column, and the
// smallest generator distance and facet length among its neighbors.

#include "core/domain.h"
#include "core/jointDual.h"
#include "core/shiftedDiagram.h"

#include "../common/boxDomain.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace
{

using PP::JointDual;
using Vector = JointDual::Vector;

bool
_ReadValues(std::istream& in, size_t count, std::vector<double>& out)
{
    out.resize(count);
    for (double& v : out) {
        if (not (in >> v)) {
            return false;
        }
    }
    return true;
}

struct Case
{
    std::vector<double> _box, _volumes, _targets, _sites;
    int                 _dim = 2;
    int                 _n = 0;
    std::string         _solver = "H";
    double              _alpha = HUGE_VAL;
};

bool
_ReadCase(const std::string& path, Case& c)
{
    std::ifstream in(path);
    std::string key;
    while (in >> key) {
        if (key == "dim") {
            in >> c._dim;
        } else if (key == "box") {
            _ReadValues(in, size_t(2 * c._dim), c._box);
        } else if (key == "n") {
            in >> c._n;
        } else if (key == "solver") {
            in >> c._solver;
        } else if (key == "alpha") {
            in >> c._alpha;
        } else if (key == "volumes") {
            _ReadValues(in, size_t(c._n), c._volumes);
        } else if (key == "targets") {
            _ReadValues(in, size_t(c._dim * c._n), c._targets);
        } else if (key == "sites") {
            _ReadValues(in, size_t(c._dim * c._n), c._sites);
        } else {
            std::string rest;
            std::getline(in, rest);  // other benchBourne keys
        }
    }
    return (c._dim == 2 or c._dim == 3) and c._box.size() == size_t(2 * c._dim) and not c._volumes.empty();
}

Vector
_Gradient(const JointDual& dual, const Vector& z, bool& ok)
{
    PP::ShiftedDiagram d;
    std::string err;
    ok = dual.Build(z, d, err) and d.GetNumEmptyCells() == 0;
    Vector g;
    dual.GetGradient(z, d, g);
    return g;
}

} // namespace

int
main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "usage: checkHessian <case.txt> <z.txt> [h ...]" << std::endl;
        return 2;
    }
    Case c;
    if (not _ReadCase(argv[1], c)) {
        std::cerr << "bad case file" << std::endl;
        return 2;
    }
    std::vector<double> zv;
    {
        std::ifstream in(argv[2]);
        double x;
        while (in >> x) {
            zv.push_back(x);
        }
    }
    const bool eig = argc >= 4 and std::string(argv[3]) == "eig";
    const int numEig = (eig and argc >= 5) ? std::atoi(argv[4]) : 8;
    std::vector<double> steps;
    for (int k = eig ? argc : 3; k < argc; ++k) {
        steps.push_back(std::atof(argv[k]));
    }
    if (steps.empty()) {
        steps = {1e-5, 1e-7, 1e-9};
    }

    std::string err;
    PP::Domain domain;
    if (not PPBench::MakeBoxDomain(c._dim, c._box, domain, err)) {
        std::cerr << err << std::endl;
        return 2;
    }
    const double alpha = c._solver == "Halpha" ? c._alpha : HUGE_VAL;
    const JointDual dual(domain, c._sites, c._targets, c._volumes, alpha, JointDual::Unknowns::ShiftsAndWeights);
    if (zv.size() != size_t(dual.GetNumVariables())) {
        std::cerr << "z has " << zv.size() << " values, expected " << dual.GetNumVariables() << std::endl;
        return 2;
    }
    const Vector z = Eigen::Map<const Vector>(zv.data(), Eigen::Index(zv.size()));

    PP::ShiftedDiagram d;
    if (not dual.Build(z, d, err)) {
        std::cerr << err << std::endl;
        return 2;
    }
    JointDual::SparseMatrix Ns;
    dual.GetNegatedHessian(d, 0.0, Ns);
    const Eigen::MatrixXd N(Ns);
    if (eig) {
        const int n = dual.GetNumSites(), s = dual.GetStride(), dim = dual.GetDimension();
        std::vector<int> iu, iw;
        for (int i = 0; i < n; ++i) {
            for (int c = 0; c < dim; ++c) {
                iu.push_back(s * i + c);
            }
            iw.push_back(s * i + dim);
        }
        auto sub = [&](const std::vector<int>& r, const std::vector<int>& c) {
            Eigen::MatrixXd M(r.size(), c.size());
            for (size_t a = 0; a < r.size(); ++a) {
                for (size_t b = 0; b < c.size(); ++b) {
                    M(a, b) = N(r[a], c[b]);
                }
            }
            return M;
        };
        const Eigen::MatrixXd Nuu = sub(iu, iu), Nuw = sub(iu, iw), Nww = sub(iw, iw);
        // N_ww is singular along the constant weight: use the pseudo-inverse
        // (its range is orthogonal to the constant, as is N_wu's).
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ew(Nww);
        Eigen::VectorXd inv = ew.eigenvalues();
        const double tol = 1e-12 * inv.cwiseAbs().maxCoeff();
        for (int k = 0; k < inv.size(); ++k) {
            inv[k] = std::fabs(inv[k]) > tol ? 1.0 / inv[k] : 0.0;
        }
        const Eigen::MatrixXd NwwPinv = ew.eigenvectors() * inv.asDiagonal() * ew.eigenvectors().transpose();
        Eigen::MatrixXd S = Nuu - Nuw * NwwPinv * Nuw.transpose();
        // Mass normalization: S' = D^-1/2 S D^-1/2 with D = 2 Vbar_i per
        // shift (the scale of the centroid term), so eigenvalues are
        // dimensionless.
        Eigen::VectorXd dinv(dim * n);
        for (int i = 0; i < n; ++i) {
            dinv.segment(dim * i, dim).setConstant(1.0 / std::sqrt(2.0 * c._volumes[i]));
        }
        S = dinv.asDiagonal() * S * dinv.asDiagonal();
        // Gauge modes (translations, dilation about the origin of the
        // generators), orthonormalized in the scaled metric, projected out.
        Eigen::MatrixXd G(dim * n, dim + 1);
        G.setZero();
        for (int i = 0; i < n; ++i) {
            for (int c2 = 0; c2 < dim; ++c2) {
                G(dim * i + c2, c2) = 1.0 / dinv[dim * i];
                G(dim * i + c2, dim) = d.GetGenerator(i)[c2] / dinv[dim * i];
            }
        }
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(G);
        const Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(dim * n, dim + 1);
        const Eigen::MatrixXd P = Eigen::MatrixXd::Identity(dim * n, dim * n) - Q * Q.transpose();
        const Eigen::MatrixXd Sp = P * S * P;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Sp);
        const Eigen::VectorXd lam = es.eigenvalues();
        // Skip the dim + 1 projected (zero) modes: take the eigenvalues above
        // the gauge cutoff, sorted ascending.
        const double big = lam.cwiseAbs().maxCoeff();
        std::printf("n %d  facets %zu  reduced Hessian (scaled): max eigenvalue %.3e\n", n, d.GetNumFacets(), big);
        int printed = 0;
        for (int k = 0; k < lam.size() and printed < numEig; ++k) {
            if (std::fabs(lam[k]) < 1e-10 * big) {
                continue;  // gauge
            }
            const Eigen::VectorXd v = dinv.asDiagonal() * es.eigenvectors().col(k);
            double num = 0.0, den = 0.0;
            for (size_t f = 0; f < d.GetNumFacets(); ++f) {
                const auto& F = d.GetFacet(f);
                num += (v.segment(dim * F._i, dim) - v.segment(dim * F._j, dim)).squaredNorm();
            }
            num /= std::max<size_t>(1, d.GetNumFacets());
            den = v.squaredNorm() / n;
            std::printf("  eig %2d: %.3e  (ratio to max %.1e)  smoothness %.3f\n", printed, lam[k], lam[k] / big, num / den);
            ++printed;
        }
        return 0;
    }
    const int s = dual.GetStride(), n = dual.GetNumSites();

    // Smallest generator distance and facet length per site.
    std::vector<double> minDist(size_t(n), HUGE_VAL), minArea(size_t(n), HUGE_VAL);
    double globalDist = HUGE_VAL, globalArea = HUGE_VAL;
    for (size_t f = 0; f < d.GetNumFacets(); ++f) {
        const auto& F = d.GetFacet(f);
        const double dist = d.GetGeneratorDistance(f);
        for (int i : {F._i, F._j}) {
            minDist[i] = std::min(minDist[i], dist);
            minArea[i] = std::min(minArea[i], F._area);
        }
        globalDist = std::min(globalDist, dist);
        globalArea = std::min(globalArea, F._area);
    }
    const double cell = std::sqrt(domain.GetVolume() / n);
    std::printf("n %d  facets %zu  min generator distance %.2e h  min facet length %.2e h\n", n, d.GetNumFacets(),
                globalDist / cell, globalArea / cell);

    for (double h : steps) {
        // Block maxima of |N_fd - N| and of |N|.
        double err[3] = {0, 0, 0}, scale[3] = {0, 0, 0};
        double worst = 0.0;
        int worstCol = -1;
        bool allOk = true;
        for (int k = 0; k < dual.GetNumVariables(); ++k) {
            Vector zp = z, zm = z;
            zp[k] += h;
            zm[k] -= h;
            bool okp, okm;
            const Vector col = (_Gradient(dual, zm, okm) - _Gradient(dual, zp, okp)) / (2 * h);
            allOk = allOk and okp and okm;
            for (int r = 0; r < dual.GetNumVariables(); ++r) {
                const int b = (r % s == s - 1 ? 1 : 0) + (k % s == s - 1 ? 1 : 0);  // 0 uu, 1 uw, 2 ww
                err[b] = std::max(err[b], std::fabs(col[r] - N(r, k)));
                scale[b] = std::max(scale[b], std::fabs(N(r, k)));
            }
            const double e = (col - N.col(k)).lpNorm<Eigen::Infinity>() / std::max(1e-300, N.col(k).lpNorm<Eigen::Infinity>());
            if (e > worst) {
                worst = e;
                worstCol = k;
            }
        }
        const int site = worstCol / s;
        std::printf("h %.0e%s: rel err uu %.2e  uw %.2e  ww %.2e | worst column %d (site %d, %s): %.2e, "
                    "site min gen dist %.2e h, min facet %.2e h\n",
                    h, allOk ? "" : " (empty cell at a probe)", err[0] / scale[0], err[1] / scale[1],
                    err[2] / scale[2], worstCol, site, worstCol % s == s - 1 ? "w" : "u", worst,
                    minDist[site] / cell, minArea[site] / cell);
    }
    return 0;
}
