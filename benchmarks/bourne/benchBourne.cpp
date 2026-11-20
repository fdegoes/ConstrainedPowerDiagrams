// Runs one dual solve of the Bourne et al. benchmarks (benchmarks/bourne)
// and writes the result and the trajectory of iterates as JSON.
//
//   benchBourne <case.txt> <out.json>
//
// The case file is whitespace-separated "key value..." records:
//   dim      2 or 3 (default 2; must precede box, targets and sites)
//   box      the lower then the upper corner (2 dim values), e.g.
//            xmin ymin xmax ymax in 2D
//   n        number of cells
//   solver   H | Halpha | G
//   alpha    stiffness of Halpha
//   start    s of MakeJointStart (0: from the sites; not used by G, which
//            starts from zero shifts)
//   volumetol, shifttol, maxiter
//   eta             DualSolver::Parameters::_regularization (default if absent)
//   noclamp         DualSolver::Parameters::_hessianClamp = false
//   pairreg R s     DualSolver::Parameters::_pairRegularization
//                   Radius, Scale
//   notrajectory    do not record the iterates (the result only)
//   lighttrajectory record the iterates without z, the shifted weights and
//                   the close pairs (for the stopping rules of other methods)
//   restarts k g    DualSolver::Parameters::_pairRestarts, _pairRestartGrowth
//   stretchlimit t  DualSolver::Parameters::_stretchLimit
//   volumes  n values
//   targets  dim n values (cbar)
//   sites    dim n values (q)
//   weights  n values, the fixed weights w of G (required by solver G)

#include "core/domain.h"
#include "core/dualSolver.h"
#include "core/jointDual.h"
#include "core/powerCells.h"

#include "../common/boxDomain.h"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{

using PP::DualSolver;
using Clock = std::chrono::steady_clock;

struct Case
{
    int                 _dim = 2;
    std::vector<double> _box;
    int                 _n = 0;
    std::string         _solver = "H";
    double              _alpha = HUGE_VAL;
    double              _start = 0.0;
    double              _volumeTol = 1e-8;
    double              _shiftTol = 1e-8;
    int                 _maxIterations = 100;
    double              _eta = -1.0;                 // < 0: the solver default
    bool                _clamp = true;
    double              _pairRadius = -1.0;  // < 0: the solver defaults (pairreg)
    double              _pairScale = 2.0;
    int                 _restarts = -1;  // < 0: the solver default
    double              _restartGrowth = 2.0;
    double              _stretchLimit = -1.0;  // < 0: the solver default
    bool                _trajectory = true;
    bool                _lightTrajectory = false;
    std::vector<double> _volumes, _targets, _sites;
    std::vector<double> _weights;  // fixed weights of G
};

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

bool
_ReadCase(const std::string& path, Case& c, std::string& err)
{
    std::ifstream in(path);
    if (not in) {
        err = "cannot open " + path;
        return false;
    }
    std::string key;
    while (in >> key) {
        bool ok = true;
        if (key == "dim") {
            ok = bool(in >> c._dim) and (c._dim == 2 or c._dim == 3);
        } else if (key == "box") {
            ok = _ReadValues(in, size_t(2 * c._dim), c._box);
        } else if (key == "n") {
            ok = bool(in >> c._n);
        } else if (key == "solver") {
            ok = bool(in >> c._solver);
        } else if (key == "alpha") {
            ok = bool(in >> c._alpha);
        } else if (key == "start") {
            ok = bool(in >> c._start);
        } else if (key == "volumetol") {
            ok = bool(in >> c._volumeTol);
        } else if (key == "shifttol") {
            ok = bool(in >> c._shiftTol);
        } else if (key == "maxiter") {
            ok = bool(in >> c._maxIterations);
        } else if (key == "eta") {
            ok = bool(in >> c._eta);
        } else if (key == "noclamp") {
            c._clamp = false;
        } else if (key == "notrajectory") {
            c._trajectory = false;
        } else if (key == "lighttrajectory") {
            c._lightTrajectory = true;
        } else if (key == "restarts") {
            ok = bool(in >> c._restarts >> c._restartGrowth);
        } else if (key == "pairreg") {
            ok = bool(in >> c._pairRadius >> c._pairScale);
        } else if (key == "stretchlimit") {
            ok = bool(in >> c._stretchLimit);
        } else if (key == "volumes") {
            ok = _ReadValues(in, size_t(c._n), c._volumes);
        } else if (key == "targets") {
            ok = _ReadValues(in, size_t(c._dim * c._n), c._targets);
        } else if (key == "sites") {
            ok = _ReadValues(in, size_t(c._dim * c._n), c._sites);
        } else if (key == "weights") {
            ok = _ReadValues(in, size_t(c._n), c._weights);
        } else {
            err = "unknown key " + key;
            return false;
        }
        if (not ok) {
            err = "bad value for " + key;
            return false;
        }
    }
    if (c._box.size() != size_t(2 * c._dim) or c._volumes.empty() or c._targets.empty() or c._sites.empty()) {
        err = "missing box, volumes, targets or sites";
        return false;
    }
    if (c._solver == "G" and c._weights.empty()) {
        err = "solver G needs the fixed weights";
        return false;
    }
    return true;
}

// A JSON number: null for nan or inf, which JSON cannot represent.
struct _Num
{
    double _v;
};

static std::ostream&
operator<<(std::ostream& out, _Num n)
{
    if (std::isfinite(n._v)) {
        out << n._v;
    } else {
        out << "null";
    }
    return out;
}

// A JSON string, with its quotes and the escapes of quotes, backslashes and
// control characters.
struct _Str
{
    const std::string& _s;
};

static std::ostream&
operator<<(std::ostream& out, _Str s)
{
    out << '"';
    for (const char ch : s._s) {
        if (ch == '"' or ch == '\\') {
            out << '\\' << ch;
        } else if (static_cast<unsigned char>(ch) < 0x20) {
            out << ' ';
        } else {
            out << ch;
        }
    }
    return out << '"';
}

void
_WriteArray(std::ostream& out, const std::vector<double>& v)
{
    out << "[";
    for (size_t k = 0; k < v.size(); ++k) {
        out << (k ? "," : "") << _Num{v[k]};
    }
    out << "]";
}

struct Iterate
{
    double              _seconds;
    int                 _numBuilds;
    double              _value;
    double              _volumeResidual;
    double              _shiftResidual;
    std::vector<double> _generators;
    std::vector<double> _shiftedWeights;
    std::vector<double> _volumes;
    std::vector<double> _centroids;
    std::vector<double> _z;
    std::vector<int>    _closePairs;  // i, j of the facets with d_ij < R (pair regularization)
};

void
_GetCells(const PP::ShiftedDiagram& d, std::vector<double>& volumes, std::vector<double>& centroids)
{
    const int n = d.GetNumSites();
    const int dim = d.GetDimension();
    volumes.resize(size_t(n));
    centroids.resize(size_t(dim * n));
    for (int i = 0; i < n; ++i) {
        volumes[i] = d.GetVolume(i);
        for (int k = 0; k < dim; ++k) {
            centroids[size_t(dim * i + k)] = d.GetCentroid(i)[k];
        }
    }
}

} // namespace

int
main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "usage: benchBourne <case.txt> <out.json>" << std::endl;
        return 2;
    }
    Case c;
    std::string err;
    if (not _ReadCase(argv[1], c, err)) {
        std::cerr << err << std::endl;
        return 2;
    }

    PP::Domain domain;
    if (not PPBench::MakeBoxDomain(c._dim, c._box, domain, err)) {
        std::cerr << err << std::endl;
        return 2;
    }

    DualSolver::Parameters params;
    params._volumeTolerance = c._volumeTol;
    params._shiftTolerance = c._shiftTol;
    params._maxIterations = c._maxIterations;
    if (c._eta >= 0.0) {
        params._regularization = c._eta;
    }
    params._hessianClamp = c._clamp;
    if (c._restarts >= 0) {
        params._pairRestarts = c._restarts;
        params._pairRestartGrowth = c._restartGrowth;
    }
    if (c._pairRadius >= 0.0) {
        params._pairRegularizationRadius = c._pairRadius;
        params._pairRegularizationScale = c._pairScale;
    }
    if (c._stretchLimit >= 0.0) {
        params._stretchLimit = c._stretchLimit;
    }

    // Radius of the pair regularization, to record the regularized pairs.
    // The cell size (|domain| / n)^(1/d), as in JointDual::_GetCellSize().
    const double cellSize = std::pow(domain.GetVolume() / double(c._n), 1.0 / double(c._dim));
    const double pairRadius = (c._pairRadius >= 0.0 ? c._pairRadius : params._pairRegularizationRadius) * cellSize;

    // The trajectory is recorded outside the timed solve: the observer
    // only copies the iterate, and its cost is subtracted.
    std::vector<Iterate> trajectory;
    double observerSeconds = 0.0;
    const Clock::time_point t0 = Clock::now();
    int startBuilds = 0;
    if (c._trajectory) params._onIterate = [&](const PP::JointDual&, const DualSolver::Vector& z, const DualSolver::Result& r) {
        const Clock::time_point s = Clock::now();
        Iterate it;
        it._seconds = std::chrono::duration<double>(s - t0).count() - observerSeconds;
        it._numBuilds = startBuilds + r._previousBuilds + r._numBuilds;
        it._value = r._value;
        it._volumeResidual = r._volumeResidual;
        it._shiftResidual = r._shiftResidual;
        it._generators = r._diagram.GetGenerators();
        _GetCells(r._diagram, it._volumes, it._centroids);
        if (c._lightTrajectory) {
            trajectory.push_back(std::move(it));
            observerSeconds += std::chrono::duration<double>(Clock::now() - s).count();
            return;
        }
        it._shiftedWeights = r._diagram.GetShiftedWeights();
        it._z.assign(z.data(), z.data() + z.size());
        for (size_t f = 0; f < r._diagram.GetNumFacets(); ++f) {
            const double d = r._diagram.GetGeneratorDistance(f);
            if (pairRadius > 0.0 and d > 0.0 and d < pairRadius) {
                it._closePairs.push_back(r._diagram.GetFacet(f)._i);
                it._closePairs.push_back(r._diagram.GetFacet(f)._j);
            }
        }
        trajectory.push_back(std::move(it));
        observerSeconds += std::chrono::duration<double>(Clock::now() - s).count();
    };

    // Start: the volume solve F at the generators (1 - s) q + s cbar. G
    // has fixed weights and no start solve: it starts from zero shifts.
    const bool fixedWeights = c._solver == "G";
    std::vector<double> shifts, weights;
    DualSolver::Result start;
    // Empty cells of the start diagram: those of (q, w) for G (zero
    // shifts), those of the F solve otherwise (-1: no start diagram).
    int startEmptyCells = -1;
    if (fixedWeights) {
        start._status = DualSolver::Status::Converged;
        shifts.assign(c._sites.size(), 0.0);
        weights = c._weights;
        PP::PowerCells cells;
        if (not PP::ComputePowerCells(domain, c._sites, c._weights, cells, err)) {
            std::cerr << err << std::endl;
            return 2;
        }
        startEmptyCells = cells.GetNumEmptyCells();
    } else {
        start = PP::MakeJointStart(domain, c._sites, c._targets, c._volumes, c._start, shifts, weights);
        if (start._diagram.GetNumSites() == c._n) {
            startEmptyCells = start._diagram.GetNumEmptyCells();
        }
    }
    startBuilds = start._numBuilds;
    const double startSeconds = std::chrono::duration<double>(Clock::now() - t0).count();

    // A start stopped by round-off (e.g. cocircular sites) close to the
    // prescribed volumes is still used; startConverged records it.
    // MakeJointStart returns before writing the shifts when the F solve does
    // not converge; with s = 0 they are zero and the weights are lambda_F.
    const bool usable = start.HasConverged() or
                        (start._status == DualSolver::Status::LineSearchFailed and start._volumeResidual <= 1e-6);
    if (usable and shifts.size() != c._sites.size()) {
        if (c._start != 0.0) {
            std::cerr << "unconverged start with s != 0" << std::endl;
            return 2;
        }
        shifts.assign(c._sites.size(), 0.0);
    }

    DualSolver::Result res = start;
    if (usable) {
        if (c._solver == "Halpha") {
            res = PP::SoftenedSolver(params).Solve(domain, c._sites, c._targets, c._volumes, c._alpha, shifts, weights);
        } else if (fixedWeights) {
            res = PP::CentroidSolver(params).Solve(domain, c._sites, c._targets, weights, c._volumes, shifts);
        } else {
            res = PP::JointSolver(params).Solve(domain, c._sites, c._targets, c._volumes, shifts, weights);
        }
    }
    const double seconds = std::chrono::duration<double>(Clock::now() - t0).count() - observerSeconds;

    std::vector<double> volumes, centroids;
    _GetCells(res._diagram, volumes, centroids);

    std::ofstream out(argv[2]);
    out.precision(17);
    out << "{\"dim\":" << c._dim << ",\"status\":\"" << DualSolver::GetStatusString(res._status) << "\""
        << ",\"converged\":" << (res.HasConverged() ? "true" : "false")
        << ",\"startConverged\":" << (start.HasConverged() ? "true" : "false")
        << ",\"startVolumeResidual\":" << _Num{start._volumeResidual} << ",\"startEmptyCells\":" << startEmptyCells
        << ",\"startIterations\":" << start._iterations
        << ",\"startBuilds\":" << start._numBuilds << ",\"startSeconds\":" << startSeconds
        << ",\"iterations\":" << res._iterations << ",\"attempts\":" << res._attempts
        << ",\"attemptIterations\":" << res._attemptIterations
        << ",\"numBuilds\":" << startBuilds + (usable ? res._numBuilds : 0) << ",\"seconds\":" << seconds
        << ",\"buildSeconds\":" << res._buildSeconds << ",\"hessianSeconds\":" << res._hessianSeconds
        << ",\"factorSeconds\":" << res._factorSeconds << ",\"volumeResidual\":" << _Num{res._volumeResidual}
        << ",\"shiftResidual\":" << _Num{res._shiftResidual} << ",\"value\":" << _Num{res._value} << ",\"error\":"
        << _Str{res._error};
    out << ",\"shifts\":";
    _WriteArray(out, shifts);
    out << ",\"weights\":";
    _WriteArray(out, weights);
    out << ",\"generators\":";
    _WriteArray(out, res._diagram.GetGenerators());
    out << ",\"shiftedWeights\":";
    _WriteArray(out, res._diagram.GetShiftedWeights());
    out << ",\"volumes\":";
    _WriteArray(out, volumes);
    out << ",\"centroids\":";
    _WriteArray(out, centroids);
    out << ",\"trajectory\":[";
    for (size_t k = 0; k < trajectory.size(); ++k) {
        const Iterate& it = trajectory[k];
        out << (k ? "," : "") << "{\"seconds\":" << it._seconds << ",\"numBuilds\":" << it._numBuilds
            << ",\"value\":" << _Num{it._value} << ",\"volumeResidual\":" << _Num{it._volumeResidual}
            << ",\"shiftResidual\":" << _Num{it._shiftResidual} << ",\"generators\":";
        _WriteArray(out, it._generators);
        out << ",\"shiftedWeights\":";
        _WriteArray(out, it._shiftedWeights);
        out << ",\"volumes\":";
        _WriteArray(out, it._volumes);
        out << ",\"centroids\":";
        _WriteArray(out, it._centroids);
        out << ",\"z\":";
        _WriteArray(out, it._z);
        out << ",\"closePairs\":[";
        for (size_t k = 0; k < it._closePairs.size(); ++k) {
            out << (k ? "," : "") << it._closePairs[k];
        }
        out << "]}";
    }
    out << "]}" << std::endl;
    return 0;
}
