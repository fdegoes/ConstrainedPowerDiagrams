#include "core/jointDual.h"
#include "core/domain.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace PP
{

JointDual::JointDual(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& targets,
    const std::vector<double>& volumes,
    double alpha,
    Unknowns unknowns)
    : _domain(domain),
      _sites(sites),
      _targets(unknowns == Unknowns::Weights ? sites : targets),
      _volumes(volumes),
      _alpha(alpha),
      _unknowns(unknowns),
      _dim(domain.GetDimension())
{
}

void
JointDual::Unpack(const Vector& z, std::vector<double>& shifts, std::vector<double>& weights) const
{
    const int n = GetNumSites();
    shifts.assign(_dim * n, 0.0);
    weights.resize(n);
    for (int i = 0; i < n; ++i) {
        if (HasShifts()) {
            for (int c = 0; c < _dim; ++c) {
                shifts[_dim * i + c] = z[GetStride() * i + c];
            }
        }
        weights[i] = HasWeights() ? z[_GetWeightIndex(i)] : _fixedWeights[i];
    }
}

JointDual::Vector
JointDual::Pack(const std::vector<double>& shifts, const std::vector<double>& weights) const
{
    Vector z = Vector::Zero(GetNumVariables());
    for (int i = 0; i < GetNumSites(); ++i) {
        if (HasShifts()) {
            for (int c = 0; c < _dim; ++c) {
                z[GetStride() * i + c] = shifts[_dim * i + c];
            }
        }
        if (HasWeights()) {
            z[_GetWeightIndex(i)] = weights[i];
        }
    }
    return z;
}

bool
JointDual::Build(const Vector& z, ShiftedDiagram& diagram, std::string& err) const
{
    if (z.size() != GetNumVariables()) {
        err = "Wrong number of unknowns.";
        return false;
    }
    if (not HasWeights() and _fixedWeights.size() != size_t(GetNumSites())) {
        err = "The fixed weights of G are not set.";
        return false;
    }
    std::vector<double> shifts, weights;
    Unpack(z, shifts, weights);
    return diagram.Build(_domain, _sites, _targets, shifts, weights, err);
}

double
JointDual::GetValue(const Vector& z, const ShiftedDiagram& diagram) const
{
    double h = diagram.GetEnergy();
    if (HasWeights()) {
        for (int i = 0; i < GetNumSites(); ++i) {
            h += z[_GetWeightIndex(i)] * _volumes[i];
        }
    }
    if (HasShifts() and std::isfinite(_alpha)) {
        for (int i = 0; i < GetNumSites(); ++i) {
            const double u2 = z.segment(GetStride() * i, _dim).squaredNorm();
            h -= _volumes[i] * u2 / _alpha;
        }
    }
    return h;
}

void
JointDual::GetGradient(const Vector& z, const ShiftedDiagram& diagram, Vector& g) const
{
    g.resize(GetNumVariables());
    for (int i = 0; i < GetNumSites(); ++i) {
        const double V = diagram.GetVolume(i);
        if (HasWeights()) {
            g[_GetWeightIndex(i)] = _volumes[i] - V;
        }
        if (not HasShifts()) {
            continue;
        }
        // 2 V_i (cbar_i - c_i) - 2 / alpha Vbar_i u_i.
        const double* c = diagram.GetCentroid(i);
        const bool soft = std::isfinite(_alpha);
        for (int k = 0; k < _dim; ++k) {
            double gk = 2.0 * V * (_targets[_dim * i + k] - c[k]);
            if (soft) {
                gk -= 2.0 / _alpha * _volumes[i] * z[GetStride() * i + k];
            }
            g[GetStride() * i + k] = gk;
        }
    }
}

JointDual::Vector
JointDual::GetMassDiagonal() const
{
    Vector m(GetNumVariables());
    for (int i = 0; i < GetNumSites(); ++i) {
        const double r = std::pow(_volumes[i], 1.0 / _dim);
        if (HasShifts()) {
            m.segment(GetStride() * i, _dim).setConstant(_volumes[i] * r * r);
        }
        if (HasWeights()) {
            m[_GetWeightIndex(i)] = _volumes[i];
        }
    }
    return m;
}

void
JointDual::GetResiduals(const Vector& g, double& volume, double& shift) const
{
    volume = 0.0;
    shift = 0.0;
    const double p = (_dim + 1.0) / _dim;
    for (int i = 0; i < GetNumSites(); ++i) {
        if (HasWeights()) {
            volume = std::max(volume, std::fabs(g[_GetWeightIndex(i)]) / _volumes[i]);
        }
        if (HasShifts()) {
            const double gu = g.segment(GetStride() * i, _dim).norm();
            shift = std::max(shift, 0.5 * gu / std::pow(_volumes[i], p));
        }
    }
}

double
JointDual::GetResidual(const Vector& g) const
{
    double volume, shift;
    GetResiduals(g, volume, shift);
    return std::max(volume / _volumeTolerance, shift / _shiftTolerance);
}

double
JointDual::_GetCellSize() const
{
    return std::pow(_domain.GetVolume() / double(GetNumSites()), 1.0 / _dim);
}

// Block-sparse assembly of N, written in place into its compressed column
// storage. Block column a holds the blocks of a and of its neighbors (the
// sites sharing a facet of positive measure) in increasing order; column k
// of block column a starts at s^2 off_a + k s deg_a, and block row b at
// list position p holds rows s b .. s b + s - 1 from p s on.
struct JointDual::_Assembly
{
    _Assembly(const ShiftedDiagram& diagram, int n, int s, SparseMatrix& N) : _s(s), _N(N)
    {
        // Neighbor lists, with the site itself.
        std::vector<int> deg(size_t(n), 1);
        for (size_t f = 0; f < diagram.GetNumFacets(); ++f) {
            if (diagram.GetGeneratorDistance(f) > 0.0) {
                ++deg[size_t(diagram.GetFacet(f)._i)];
                ++deg[size_t(diagram.GetFacet(f)._j)];
            }
        }
        _off.assign(size_t(n) + 1, 0);
        for (int a = 0; a < n; ++a) {
            _off[size_t(a) + 1] = _off[size_t(a)] + deg[size_t(a)];
        }
        std::vector<int> lists(static_cast<size_t>(_off[size_t(n)]));
        std::vector<int> fill(_off.begin(), _off.end() - 1);
        for (int a = 0; a < n; ++a) {
            lists[size_t(fill[size_t(a)]++)] = a;
        }
        for (size_t f = 0; f < diagram.GetNumFacets(); ++f) {
            if (diagram.GetGeneratorDistance(f) > 0.0) {
                const int i = diagram.GetFacet(f)._i, j = diagram.GetFacet(f)._j;
                lists[size_t(fill[size_t(i)]++)] = j;
                lists[size_t(fill[size_t(j)]++)] = i;
            }
        }
        _self.resize(size_t(n));
        for (int a = 0; a < n; ++a) {
            std::sort(lists.begin() + _off[size_t(a)], lists.begin() + _off[size_t(a) + 1]);
            _self[size_t(a)] = _Find(lists, a, a);
        }
        // Positions of the off-diagonal blocks of each facet.
        _facetPos.assign(diagram.GetNumFacets(), {0, 0});
        for (size_t f = 0; f < diagram.GetNumFacets(); ++f) {
            if (diagram.GetGeneratorDistance(f) > 0.0) {
                const int i = diagram.GetFacet(f)._i, j = diagram.GetFacet(f)._j;
                _facetPos[f] = {_Find(lists, j, i), _Find(lists, i, j)};
            }
        }

        const Eigen::Index size = Eigen::Index(s) * n;
        const Eigen::Index nnz = Eigen::Index(s) * s * _off[size_t(n)];
        N.resize(size, size);
        N.resizeNonZeros(nnz);
        int* outer = N.outerIndexPtr();
        int* inner = N.innerIndexPtr();
        for (int a = 0; a < n; ++a) {
            const int d = deg[size_t(a)];
            for (int k = 0; k < s; ++k) {
                const int start = _GetColumnStart(a, k);
                outer[s * a + k] = start;
                for (int p = 0; p < d; ++p) {
                    for (int r = 0; r < s; ++r) {
                        inner[start + p * s + r] = s * lists[size_t(_off[size_t(a)] + p)] + r;
                    }
                }
            }
        }
        outer[size] = int(nnz);
        std::fill(N.valuePtr(), N.valuePtr() + nnz, 0.0);
    }

    // Entry (r, c) of the block of rows of site b in block column a, b at
    // position p of the list of a.
    double& At(int a, int p, int r, int c)
    {
        return _N.valuePtr()[_GetColumnStart(a, c) + p * _s + r];
    }
    double& Diagonal(int a, int r, int c) { return At(a, _self[size_t(a)], r, c); }

    int                              _s;
    SparseMatrix&                    _N;
    std::vector<int>                 _off;
    std::vector<int>                 _self;      // position of a in its list
    std::vector<std::array<int, 2>>  _facetPos;  // (i in list of j, j in list of i)

private:
    int _GetColumnStart(int a, int k) const
    {
        const int d = _off[size_t(a) + 1] - _off[size_t(a)];
        return _s * _s * _off[size_t(a)] + k * _s * d;
    }
    // Position of b in the list of a.
    int _Find(const std::vector<int>& lists, int a, int b) const
    {
        const auto first = lists.begin() + _off[size_t(a)];
        const auto last = lists.begin() + _off[size_t(a) + 1];
        return int(std::lower_bound(first, last, b) - first);
    }
};

void
JointDual::_AddPairRegularization(const ShiftedDiagram& diagram, _Assembly& N) const
{
    // Shift-shift diagonal blocks of the Hessian assembled so far (the facet
    // blocks), copied before any pair term is added.
    const int n = GetNumSites();
    const int dd = _dim * _dim;
    std::vector<double> blocks(size_t(n) * size_t(dd));
    for (int a = 0; a < n; ++a) {
        for (int r = 0; r < _dim; ++r) {
            for (int c = 0; c < _dim; ++c) {
                blocks[size_t(a) * dd + size_t(_dim * r + c)] = N.Diagonal(a, r, c);
            }
        }
    }
    // Curvature e^T B_a e of site a along the unit direction e.
    auto curvature = [&](int a, const double* e) {
        const double* B = &blocks[size_t(a) * dd];
        double k = 0.0;
        for (int r = 0; r < _dim; ++r) {
            for (int c = 0; c < _dim; ++c) {
                k += e[r] * B[_dim * r + c] * e[c];
            }
        }
        return k;
    };
    const double cell = _GetCellSize();
    for (size_t f = 0; f < diagram.GetNumFacets(); ++f) {
        const double d = diagram.GetGeneratorDistance(f);
        if (not (d > 0.0)) {
            continue;  // coincident (no direction)
        }
        const ShiftedDiagram::Facet& F = diagram.GetFacet(f);
        const double radius = _pairActivity * _pairRadiusRatio * cell;
        if (d >= radius) {
            continue;  // far enough apart
        }
        const double* pi = diagram.GetGenerator(F._i);
        const double* pj = diagram.GetGenerator(F._j);
        // Unit direction of the pair.
        const double inv = 1.0 / d;
        double e[3];
        for (int c = 0; c < _dim; ++c) {
            e[c] = (pj[c] - pi[c]) * inv;
        }
        // The stiffness of the pair, from the curvature of its own sites
        // along its direction.
        const double rhoMax = _pairScale * std::max(curvature(F._i, e), curvature(F._j, e));
        if (not (rhoMax > 0.0)) {
            continue;
        }
        // rhoMax at d = 0, falling to 0 with zero slope at the radius.
        const double x = 1.0 - d / radius;
        const double rho = rhoMax * (x * x) * (x * x);
        for (int a : {F._i, F._j}) {
            for (int r = 0; r < _dim; ++r) {
                for (int c = 0; c < _dim; ++c) {
                    N.Diagonal(a, r, c) += rho * e[r] * e[c];
                }
            }
        }
    }
}

void
JointDual::_AddFacetBlocks(
    const ShiftedDiagram& diagram,
    size_t f,
    double minDistance,
    _Assembly& N) const
{
    // The Hessian blocks of the facet are 1/(2 d_ij) ∫ phi_a phi_b^T with
    // phi_a = [2 (x - cbar_a), 1] (eq:appendix-block): +-M_ab for the
    // pair (a, b), off-diagonal with sign -1 in N and diagonal with sign +1.
    const ShiftedDiagram::Facet& F = diagram.GetFacet(f);
    const int ends[2] = {F._i, F._j};
    const int s = GetStride();
    const std::array<int, 2>& pos = N._facetPos[f];
    const double ww = diagram.GetWeightWeightCoef(f, minDistance);
    double B[9], ba[3], bb[3];
    for (int ka = 0; ka < 2; ++ka) {
        for (int kb = 0; kb < 2; ++kb) {
            const int a = ends[ka];
            const int b = ends[kb];
            const double sign = (a == b) ? 1.0 : -1.0;
            // The diagonal block of cell a gathers the (a, a) terms; the
            // block (a, b) is in block column b, at the position of a.
            const int column = a == b ? a : b;
            const int p = a == b ? N._self[size_t(a)] : pos[ka];
            if (HasWeights()) {
                N.At(column, p, s - 1, s - 1) += sign * ww;
            }
            if (not HasShifts()) {
                continue;
            }
            diagram.GetShiftShiftBlock(f, a, b, B, minDistance);
            for (int r = 0; r < _dim; ++r) {
                for (int c = 0; c < _dim; ++c) {
                    N.At(column, p, r, c) += sign * B[_dim * r + c];
                }
            }
            if (not HasWeights()) {
                continue;
            }
            diagram.GetShiftWeightBlock(f, a, ba, minDistance);
            diagram.GetShiftWeightBlock(f, b, bb, minDistance);
            for (int r = 0; r < _dim; ++r) {
                N.At(column, p, r, _dim) += sign * ba[r];
                N.At(column, p, _dim, r) += sign * bb[r];
            }
        }
    }
}

void
JointDual::GetNegatedHessian(const ShiftedDiagram& diagram, double eta, SparseMatrix& N) const
{
    // Hessian distance clamp, relative to the cell size.
    const double minDistance = _hessianMinDistanceRatio > 0.0 ? _hessianMinDistanceRatio * _GetCellSize() : 0.0;
    const int n = GetNumSites();
    const int s = GetStride();
    _Assembly A(diagram, n, s, N);
    for (size_t f = 0; f < diagram.GetNumFacets(); ++f) {
        if (not (diagram.GetGeneratorDistance(f) > 0.0)) {
            continue;
        }  // coincident generators share no facet of positive measure
        _AddFacetBlocks(diagram, f, minDistance, A);
    }
    if (HasShifts() and _pairRadiusRatio > 0.0 and _pairActivity > 0.0) {
        _AddPairRegularization(diagram, A);
    }
    const bool soft = HasShifts() and std::isfinite(_alpha);
    const Vector M = GetMassDiagonal();
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < s; ++k) {
            double d = eta * M[s * i + k];
            if (soft and k < _dim) {
                d += 2.0 / _alpha * _volumes[i];
            }
            A.Diagonal(i, k, k) += d;
        }
    }
}

} // namespace PP
