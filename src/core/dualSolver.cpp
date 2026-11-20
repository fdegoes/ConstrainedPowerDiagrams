#include "core/dualSolver.h"
#include "core/domain.h"
#include "core/geometry.h"

#include <Eigen/Dense>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <ostream>

namespace PP
{

namespace
{

// Appends a sentence to a (possibly empty) message.
void
_AppendMessage(std::string& message, const std::string& sentence)
{
    if (not message.empty()) {
        message += " ";
    }
    message += sentence;
}

// True if two generators of the diagram are equal: their cells cannot both
// be non-empty (e.g. duplicate sites with equal shifts and weights).
bool
_HasCoincidentGenerators(const ShiftedDiagram& diagram)
{
    const int dim = diagram.GetDimension();
    std::vector<int> order(diagram.GetNumSites());
    std::iota(order.begin(), order.end(), 0);
    const auto less = [&](int a, int b) {
        const double* pa = diagram.GetGenerator(a);
        const double* pb = diagram.GetGenerator(b);
        return std::lexicographical_compare(pa, pa + dim, pb, pb + dim);
    };
    std::sort(order.begin(), order.end(), less);
    for (size_t k = 1; k < order.size(); ++k) {
        if (not less(order[k - 1], order[k])) {
            return true;
        }
    }
    return false;
}

bool
_CheckDualInputs(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& targets,
    const std::vector<double>& volumes,
    std::string& err)
{
    if (not domain.IsValid()) {
        err = "Invalid clipping domain.";
        return false;
    }
    const size_t n = volumes.size();
    const size_t dim = size_t(domain.GetDimension());
    if (n == 0 or sites.size() != dim * n or targets.size() != dim * n) {
        err = "Number of sites, targets and target volumes do not match.";
        return false;
    }
    if (not std::all_of(volumes.begin(), volumes.end(), [](double v) { return v > 0.0 and std::isfinite(v); })) {
        err = "Target volumes must be positive.";
        return false;
    }
    if (not AreAllFinite(sites) or not AreAllFinite(targets)) {
        err = "Sites or targets are not finite.";
        return false;
    }
    return true;
}

// Multipliers z from shifts and weights, zeros where the sizes do not
// match the dual.
JointDual::Vector
_PackInitial(const JointDual& dual, std::vector<double>& shifts, std::vector<double>& weights)
{
    const size_t n = size_t(dual.GetNumSites());
    if (shifts.size() != size_t(dual.GetDimension()) * n) {
        shifts.assign(dual.GetDimension() * n, 0.0);
    }
    if (weights.size() != n) {
        weights.assign(n, 0.0);
    }
    return dual.Pack(shifts, weights);
}

double
_Seconds(std::chrono::steady_clock::time_point start)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

// Objective changes below this are indistinguishable from round-off.
double
_ComputeRoundoff(double value)
{
    return 16.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::fabs(value));
}

DualSolver::Result
_MakeInputError(const std::string& err)
{
    DualSolver::Result res;
    res._error = err;
    return res;
}

} // namespace

const char*
DualSolver::GetStatusString(Status s)
{
    switch (s) {
    case Status::Converged:           return "converged";
    case Status::MaxIterations:       return "max iterations";
    case Status::LineSearchFailed:    return "line search failed";
    case Status::DampingFailed:       return "damping failed";
    case Status::FactorizationFailed: return "factorization failed";
    case Status::InvalidInput:        return "invalid input";
    }
    return "unknown";
}

bool
DualSolver::_Start(const JointDual& dual, const Vector& z, Result& res)
{
    const auto t0 = std::chrono::steady_clock::now();
    const bool built = dual.Build(z, res._diagram, res._error);
    res._buildSeconds += _Seconds(t0);
    if (not built) {
        res._status = Status::InvalidInput;
        return false;
    }
    if (_HasCoincidentGenerators(res._diagram)) {
        res._status = Status::InvalidInput;
        res._error = "Some generators coincide at the start (e.g. duplicate sites), so their cells cannot both be non-empty.";
        return false;
    }
    const int empty = res._diagram.GetNumEmptyCells();
    if (empty > 0) {
        _AppendMessage(res._warning, std::to_string(empty) + " cell(s) are empty at the start.");
    }
    res._value = dual.GetValue(z, res._diagram);
    if (not std::isfinite(res._value)) {
        res._status = Status::InvalidInput;
        res._error = "The objective is not finite at the start.";
        return false;
    }
    return true;
}

bool
DualSolver::_Factorize(const JointDual& dual, const ShiftedDiagram& diagram, Result& res)
{
    const auto t0 = std::chrono::steady_clock::now();
    dual.GetNegatedHessian(diagram, _params._regularization, _N);
    res._hessianSeconds += _Seconds(t0);
    _borderRows.clear();
    const auto t1 = std::chrono::steady_clock::now();
    // The damped matrices N + mu D^2 share the pattern of N (its diagonal
    // is stored): one symbolic analysis per iteration.
    // N is made of (d + 1) x (d + 1) blocks, one per pair of neighbor sites.
    _llt.SetBlockSize(dual.GetStride());
    _llt.analyzePattern(_N);
    _damped = _N;
    _dampedDiagonal.resize(size_t(_N.rows()));
    for (Eigen::Index k = 0; k < _N.outerSize(); ++k) {
        SparseMatrix::InnerIterator it(_damped, k);
        while (it and it.row() < k) {
            ++it;
        }
        if (not it or it.row() != k) {
            return false;  // not reached: GetNegatedHessian() stores the diagonal
        }
        _dampedDiagonal[size_t(k)] = &it.valueRef();
    }
    bool ok = true;
    if (_UsesUndampedStep()) {
        _llt.factorize(_N);
        ok = _llt.info() == Eigen::Success;
    }
    res._factorSeconds += _Seconds(t1);
    return ok;
}

bool
DualSolver::_ComputeStep(
    const JointDual& dual,
    const Vector&,
    const ShiftedDiagram&,
    const Vector& g,
    Vector& dz,
    Vector& stationarity)
{
    dz = _SolveFactorized(g);
    stationarity = g;
    return _llt.info() == Eigen::Success;
}

void
DualSolver::_Log(const Result& res, double stepLength) const
{
    if (_params._log == nullptr) {
        return;
    }
    *_params._log << "  iter " << res._iterations << "  value " << res._value << "  volume residual "
                  << res._volumeResidual << "  shift residual " << res._shiftResidual << "  step " << stepLength
                  << std::endl;
}

bool
DualSolver::_TakeLineSearchStep(const JointDual& dual, const Vector& g, const Vector& dz, Vector& z, Result& res)
{
    const auto t0 = std::chrono::steady_clock::now();
    const LineSearch::Result ls =
        LineSearch(_params._lineSearch).Search(dual, res._value, g, res._residual, dz, z, res._diagram);
    res._buildSeconds += _Seconds(t0);
    res._numBuilds += ls._halvings + 1;
    if (not ls._accepted) {
        // ls._error: why the last trial diagram could not be built, if so.
        res._status = Status::LineSearchFailed;
        res._error = ls._error;
        return false;
    }
    res._value = ls._value;
    ++res._iterations;
    _Log(res, ls._stepLength);
    return true;
}

double
DualSolver::_ComputeTrialResidual(const JointDual& dual, const Vector& g) const
{
    if (_borderRows.empty()) {
        return dual.GetResidual(g);
    }
    // Residual of g - 2 C^T nu with nu minimizing it in the least-squares
    // sense, as at the current iterate (_ComputeStationarity()).
    Eigen::VectorXd nu;
    Vector stationarity;
    if (not _ProjectOutRows(_borderRows, _borderMetric, g, nu, stationarity)) {
        return HUGE_VAL;
    }
    return dual.GetResidual(stationarity);
}

bool
DualSolver::_ProjectOutRows(
    const std::vector<Vector>& rows,
    const Vector& metric,
    const Vector& g,
    Eigen::VectorXd& nu,
    Vector& stationarity) const
{
    const int m = int(rows.size());
    Eigen::MatrixXd G(m, m);
    Eigen::VectorXd r(m);
    for (int k = 0; k < m; ++k) {
        const Vector wk = metric.cwiseProduct(rows[k]);
        r[k] = wk.dot(g);
        for (int l = 0; l < m; ++l) {
            G(k, l) = 2.0 * wk.dot(rows[l]);
        }
    }
    const Eigen::LDLT<Eigen::MatrixXd> ldlt(G);
    if (ldlt.info() != Eigen::Success) {
        return false;
    }
    nu = ldlt.solve(r);
    stationarity = g;
    for (int k = 0; k < m; ++k) {
        stationarity -= 2.0 * nu[k] * rows[k];
    }
    return true;
}

bool
DualSolver::_ComputeUndampedMultipliers(const std::vector<Vector>& rows, const Vector& g, Eigen::VectorXd& nu)
{
    _llt.factorize(_N);
    if (_llt.info() != Eigen::Success) {
        return false;
    }
    Vector dz, stationarity;
    return _ComputeBorderedStep(rows, g, dz, nu, stationarity);
}

bool
DualSolver::_ComputeStationarity(const JointDual&, const Vector& g, Vector& stationarity)
{
    stationarity = g;
    return true;
}

double
DualSolver::_ComputeModelIncrease(const Vector& g, const Vector& p) const
{
    return g.dot(p) - 0.5 * p.dot(_N * p);
}

bool
DualSolver::_ComputeDampedStep(const JointDual& dual, const Vector& g, double mu, Vector& p)
{
    std::copy(_N.valuePtr(), _N.valuePtr() + _N.nonZeros(), _damped.valuePtr());
    for (size_t k = 0; k < _dampedDiagonal.size(); ++k) {
        *_dampedDiagonal[k] += mu * _scale2[Eigen::Index(k)];
    }
    _llt.factorize(_damped);
    if (_llt.info() != Eigen::Success) {
        return false;
    }
    if (_borderRows.empty()) {
        p = _SolveFactorized(g);
        return true;
    }
    const std::vector<Vector> rows = _borderRows;
    Eigen::VectorXd nu;
    Vector stationarity;
    return _ComputeBorderedStep(rows, g, p, nu, stationarity);
}

bool
DualSolver::_TakeDampedStep(const JointDual& dual, const Vector& g, Vector& z, Result& res)
{
    const Damping& lm = _params._damping;
    _scale2 = _N.diagonal().cwiseMax(1e-300);
    if (_resetDamping) {
        _mu = lm._initialDamping;
        _muGrowth = 2.0;
        _resetDamping = false;
    }
    ShiftedDiagram trial;
    Vector p;
    std::string trialError;
    for (int k = 0; k <= lm._maxRejections; ++k) {
        if (_mu > lm._maxDamping) {
            break;
        }
        const auto t0 = std::chrono::steady_clock::now();
        const bool stepped = _ComputeDampedStep(dual, g, _mu, p);
        res._factorSeconds += _Seconds(t0);
        if (not stepped) {
            res._status = Status::FactorizationFailed;
            return false;
        }
        if (_params._stretchLimit > 0.0 and dual.HasFrozenEdges()) {
            // Fraction to the boundary: shorten the step so that no frozen
            // edge shrinks below the limit, rather than rejecting it (which
            // would also raise the damping of the next steps).
            p *= dual.GetMaxStep(z, p, _params._stretchLimit);
        }
        const double predicted = _ComputeModelIncrease(g, p);
        const Vector zTrial = z + p;
        ++res._numBuilds;
        const auto t1 = std::chrono::steady_clock::now();
        // A trial whose diagram cannot be built (e.g. an undefined power
        // vertex of nearly collinear generators) is rejected like any other.
        const bool built = dual.Build(zTrial, trial, trialError);
        res._buildSeconds += _Seconds(t1);
        double value = res._value;
        double rho = -1.0;
        if (built and _AcceptTrial(dual, zTrial, trial, predicted, lm._acceptRatio, res, value, rho)) {
            // Nielsen: mu *= max(1/3, 1 - (2 rho - 1)^3).
            const double t = 2.0 * rho - 1.0;
            _mu = std::max(lm._minDamping, _mu * std::max(1.0 / 3.0, 1.0 - t * t * t));
            _muGrowth = 2.0;
            z = zTrial;
            res._diagram = std::move(trial);
            res._value = value;
            ++res._iterations;
            _Log(res, _mu);
            return true;
        }
        _mu = std::max(_mu, 1e-300) * _muGrowth;
        _muGrowth *= 2.0;
    }
    res._status = Status::DampingFailed;
    if (not trialError.empty()) {
        res._error = trialError;
    }
    return false;
}

bool
DualSolver::_AcceptTrial(
    const JointDual& dual,
    const Vector& zTrial,
    const ShiftedDiagram& trial,
    double predicted,
    double acceptRatio,
    const Result& res,
    double& value,
    double& rho) const
{
    rho = -1.0;
    if (not (trial.GetMinVolume() > _params._lineSearch._minVolume)) {
        return false;
    }
    if (_params._stretchLimit > 0.0 and dual.HasFrozenEdges() and
        not (dual.GetMinStretch(zTrial) >= _params._stretchLimit * (1.0 - 1e-12))) {
        return false;
    }
    value = dual.GetValue(zTrial, trial);
    const double roundoff = _ComputeRoundoff(res._value);
    if (predicted <= roundoff) {
        Vector gTrial;
        dual.GetGradient(zTrial, trial, gTrial);
        // The increase is not measurable: accept if the residual decreases
        // (the same residual as res._residual: of the stationarity condition
        // with the border rows, when there are any) and the objective does
        // not decrease beyond round-off.
        rho = (value >= res._value - roundoff and _ComputeTrialResidual(dual, gTrial) < res._residual) ? 1.0 : -1.0;
        return rho > 0.0;
    }
    rho = (value - res._value) / predicted;
    return rho > acceptRatio;
}

DualSolver::Result
DualSolver::_Solve(JointDual& dual, Vector& z, bool useLineSearch)
{
    _useLineSearch = useLineSearch;
    Result res;
    res._previousBuilds = _buildOffset;
    if (not _Start(dual, z, res)) {
        return res;
    }
    res._numBuilds = 1;

    dual.SetTolerances(_params._volumeTolerance, _params._shiftTolerance);
    dual.SetHessianMinDistanceRatio(_params._hessianClamp ? _params._regularization : 0.0);
    dual.SetPairRegularization(_params._pairRegularizationRadius, _params._pairRegularizationScale);
    dual.ClearFrozenEdges();
    const bool freeze = dual.HasShifts() and not _useLineSearch and _params._stretchLimit > 0.0;
    double startShiftResidual = -1.0;
    _resetDamping = true;
    Vector g, dz, stationarity;
    for (;;) {
        dual.GetGradient(z, res._diagram, g);
        if (_params._pairRegularizationRadius > 0.0 and dual.HasShifts()) {
            // Activity of the pair regularization from the shift residual of
            // the gradient: 1 at the start, 0 at convergence.
            double volume, shift;
            dual.GetResiduals(g, volume, shift);
            if (startShiftResidual < 0.0) {
                startShiftResidual = shift;
            }
            dual.SetPairActivity(startShiftResidual > 0.0 ? std::min(1.0, shift / startShiftResidual) : 0.0);
        }
        if (freeze) {
            // Edges of this diagram, whose shrinking the steps limit.
            dual.FreezeEdges(res._diagram);
        }
        if (not _Factorize(dual, res._diagram, res)) {
            res._status = Status::FactorizationFailed;
            return res;
        }
        const auto t0 = std::chrono::steady_clock::now();
        bool stepped = true;
        if (_UsesUndampedStep()) {
            stepped = _ComputeStep(dual, z, res._diagram, g, dz, stationarity);
        } else {
            stepped = _ComputeStationarity(dual, g, stationarity);  // the damped steps alone use N
        }
        res._factorSeconds += _Seconds(t0);
        if (not stepped) {
            res._status = Status::FactorizationFailed;
            return res;
        }
        dual.GetResiduals(stationarity, res._volumeResidual, res._shiftResidual);
        res._residual = dual.GetResidual(stationarity);
        if (_params._onIterate) {
            _params._onIterate(dual, z, res);
        }
        if (res._residual <= 1.0) {
            res._status = Status::Converged;
            return res;
        }
        if (res._iterations >= _params._maxIterations) {
            res._status = Status::MaxIterations;
            return res;
        }
        const bool ok = useLineSearch ? _TakeLineSearchStep(dual, g, dz, z, res)
                                      : _TakeDampedStep(dual, g, z, res);
        if (not ok) {
            return res;
        }
    }
}

DualSolver::Result
SoftenedSolver::Solve(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& targets,
    const std::vector<double>& volumes,
    double alpha,
    std::vector<double>& shifts,
    std::vector<double>& weights)
{
    std::string err;
    if (not _CheckDualInputs(domain, sites, targets, volumes, err)) {
        return _MakeInputError(err);
    }
    if (not (alpha > 0.0) or not std::isfinite(alpha)) {
        return _MakeInputError("The stiffness alpha must be positive and finite.");
    }
    JointDual dual(domain, sites, targets, volumes, alpha, JointDual::Unknowns::ShiftsAndWeights);
    Vector z = _PackInitial(dual, shifts, weights);
    Result res = _Solve(dual, z);
    dual.Unpack(z, shifts, weights);
    return res;
}

DualSolver::Result
CentroidSolver::Solve(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& targets,
    const std::vector<double>& weights,
    const std::vector<double>& volumes,
    std::vector<double>& shifts)
{
    std::string err;
    if (not _CheckDualInputs(domain, sites, targets, volumes, err)) {
        return _MakeInputError(err);
    }
    if (weights.size() != volumes.size() or not AreAllFinite(weights)) {
        return _MakeInputError("The fixed weights must be one finite value per site.");
    }
    JointDual dual(domain, sites, targets, volumes, HUGE_VAL, JointDual::Unknowns::Shifts);
    dual.SetFixedWeights(weights);
    std::vector<double> w = weights;
    Vector z = _PackInitial(dual, shifts, w);
    Result res = _Solve(dual, z);
    dual.Unpack(z, shifts, w);
    return res;
}

DualSolver::Result
JointSolver::Solve(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& targets,
    const std::vector<double>& volumes,
    std::vector<double>& shifts,
    std::vector<double>& weights)
{
    std::string err;
    if (not _CheckDualInputs(domain, sites, targets, volumes, err)) {
        return _MakeInputError(err);
    }
    JointDual dual(domain, sites, targets, volumes, HUGE_VAL, JointDual::Unknowns::ShiftsAndWeights);
    const Vector z0 = _PackInitial(dual, shifts, weights);
    Vector z = z0;
    _buildOffset = 0;
    Result res = _Solve(dual, z);
    // Restarts from the start z0 with a wider pair regularization (see Parameters::_pairRestarts): it changes the early
    // iterations, which decide whether close generator pairs collapse. The
    // result is the converged attempt, or the one with the smallest
    // residual; counts and timers add up over the attempts.
    const Parameters saved = GetParameters();
    Result best = res;
    Vector zBest = z;
    int iterations = res._iterations, builds = res._numBuilds;
    double buildSeconds = res._buildSeconds, hessianSeconds = res._hessianSeconds;
    double factorSeconds = res._factorSeconds;
    std::string warning = res._warning;
    int attempts = 1;
    best._attemptIterations = res._iterations;
    for (int k = 0; k < saved._pairRestarts and saved._pairRegularizationRadius > 0.0; ++k) {
        if (res._status != Status::MaxIterations and res._status != Status::DampingFailed) {
            break;
        }
        GetParameters()._pairRegularizationRadius *= saved._pairRestartGrowth;
        _buildOffset = builds;
        z = z0;
        res = _Solve(dual, z);
        ++attempts;
        iterations += res._iterations;
        builds += res._numBuilds;
        buildSeconds += res._buildSeconds;
        hessianSeconds += res._hessianSeconds;
        factorSeconds += res._factorSeconds;
        if (not res._warning.empty()) {
            warning += (warning.empty() ? "" : " ") + res._warning;
        }
        // Only attempts that ran their iterations have a residual (a failed
        // factorization or start keeps the default 0).
        const bool ran = res._status == Status::MaxIterations or res._status == Status::DampingFailed;
        if (res.HasConverged() or (ran and res._residual < best._residual)) {
            best = res;
            best._attemptIterations = res._iterations;
            zBest = z;
        }
    }
    GetParameters() = saved;
    _buildOffset = 0;
    best._iterations = iterations;
    best._numBuilds = builds;
    best._buildSeconds = buildSeconds;
    best._hessianSeconds = hessianSeconds;
    best._factorSeconds = factorSeconds;
    best._warning = warning;
    best._attempts = attempts;
    dual.Unpack(zBest, shifts, weights);
    return best;
}

DualSolver::Result
MakeJointStart(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& targets,
    const std::vector<double>& volumes,
    double s,
    std::vector<double>& shifts,
    std::vector<double>& weights)
{
    if (sites.size() != targets.size()) {
        return _MakeInputError("Number of sites and targets do not match.");
    }
    const int dim = domain.GetDimension();
    std::vector<double> generators(sites.size());
    for (size_t k = 0; k < sites.size(); ++k) {
        generators[k] = (1.0 - s) * sites[k] + s * targets[k];
    }
    DualSolver::Result res = VolumeSolver().Solve(domain, generators, volumes, weights);
    if (not res.HasConverged()) {
        return res;
    }
    // weights holds lambda_F; convert to the weights of H.
    shifts.resize(sites.size());
    for (size_t i = 0; i < weights.size(); ++i) {
        double w = weights[i];
        for (int c = 0; c < dim; ++c) {
            const size_t k = dim * i + c;
            shifts[k] = generators[k] - sites[k];
            w -= shifts[k] * shifts[k] + 2.0 * shifts[k] * (sites[k] - targets[k]);
        }
        weights[i] = w;
    }
    return res;
}

void
ConstrainedJointSolver::_BuildConstraints(
    const JointDual& dual,
    const std::vector<double>& direction,
    const std::vector<CutCellSide>& sides)
{
    // a = ({Vbar_i p_i}, 0), b_j = ({Vbar_i e_j}, 0) and
    // c_B = ({Vbar_i (p_i - p_{c_B})}, 0), zero on the weights.
    const int dim = dual.GetDimension();
    const int s = dual.GetStride();
    const std::vector<double>& volumes = dual.GetVolumes();
    _rows.assign(dim + 1 + sides.size(), Vector::Zero(dual.GetNumVariables()));
    for (int i = 0; i < dual.GetNumSites(); ++i) {
        for (int c = 0; c < dim; ++c) {
            _rows[0][s * i + c] = volumes[i] * direction[dim * i + c];
            _rows[1 + c][s * i + c] = volumes[i];
        }
    }
    for (size_t l = 0; l < sides.size(); ++l) {
        Vector& row = _rows[dim + 1 + l];
        const int cut = sides[l]._cut;
        for (int i : sides[l]._side) {
            for (int c = 0; c < dim; ++c) {
                row[s * i + c] = volumes[i] * (direction[dim * i + c] - direction[dim * cut + c]);
            }
        }
    }
}

bool
DualSolver::_ComputeBorderedStep(
    const std::vector<Vector>& rows,
    const Vector& g,
    Vector& dz,
    Eigen::VectorXd& nu,
    Vector& stationarity)
{
    _borderRows = rows;
    // One solve per row plus one for g, sharing the factorization of N
    // (eq:flow-constrained-elimination-solve).
    const int m = int(rows.size());
    const Vector z0 = _SolveFactorized(g);
    std::vector<Vector> zc(m);
    for (int k = 0; k < m; ++k) {
        zc[k] = _SolveFactorized(2.0 * rows[k]);
    }

    // Multipliers from C dz = 0 (eq:flow-constrained-elimination-gauge).
    Eigen::MatrixXd G(m, m);
    Eigen::VectorXd r(m);
    for (int k = 0; k < m; ++k) {
        r[k] = rows[k].dot(z0);
        for (int l = 0; l < m; ++l) {
            G(k, l) = rows[k].dot(zc[l]);
        }
    }
    const Eigen::LDLT<Eigen::MatrixXd> ldlt(G);
    if (ldlt.info() != Eigen::Success) {
        return false;
    }
    nu = ldlt.solve(r);

    // dz = z0 - sum_k nu_k zc_k (eq:flow-constrained-elimination-update)
    // and the stationarity residual g - 2 C^T nu.
    dz = z0;
    stationarity = g;
    for (int k = 0; k < m; ++k) {
        dz -= nu[k] * zc[k];
        stationarity -= 2.0 * nu[k] * rows[k];
    }
    return true;
}

bool
ConstrainedJointSolver::_ComputeStationarity(const JointDual& dual, const Vector& g, Vector& stationarity)
{
    // The damped steps that follow are bordered by C; the multipliers of the
    // stationarity condition, by least squares, give the gauge, which is the
    // one of the solution at convergence (g = 2 C^T nu there).
    // Metric 1 / Vbar_i on the shifts: the residual g_u_i / (2 Vbar_i) is the
    // gap between the centroid and the particle of the correction
    // (eq:flow-shake-step-corr), so that sum_i Vbar_i times this gap, i.e.
    // the momentum change, vanishes with the translation rows. The weights
    // have no rows.
    Vector metric = Vector::Ones(dual.GetNumVariables());
    const int s = dual.GetStride(), dim = dual.GetDimension();
    for (int i = 0; i < dual.GetNumSites(); ++i) {
        for (int c = 0; c < dim; ++c) {
            metric[s * i + c] = 1.0 / dual.GetVolumes()[size_t(i)];
        }
    }
    _SetBorderRows(_rows, metric);
    Eigen::VectorXd nu;
    if (not _ProjectOutRows(_rows, metric, g, nu, stationarity)) {
        return false;
    }
    _gauge._dilation = nu[0];
    _gauge._translation.assign(nu.data() + 1, nu.data() + 1 + dim);
    _gauge._arms.assign(nu.data() + 1 + dim, nu.data() + nu.size());
    return true;
}

DualSolver::Result
ConstrainedJointSolver::Solve(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& targets,
    const std::vector<double>& volumes,
    const std::vector<double>& direction,
    const std::vector<CutCellSide>& sides,
    std::vector<double>& shifts,
    std::vector<double>& weights,
    Gauge& gauge)
{
    std::string err;
    if (not _CheckDualInputs(domain, sites, targets, volumes, err)) {
        return _MakeInputError(err);
    }
    if (direction.size() != sites.size() or not AreAllFinite(direction)) {
        return _MakeInputError("The direction must have one finite point per site.");
    }
    const int numSites = int(volumes.size());
    for (const CutCellSide& side : sides) {
        if (side._cut < 0 or side._cut >= numSites or side._side.empty()) {
            return _MakeInputError("Every cut cell must be a site and have a non-empty side.");
        }
        for (int i : side._side) {
            if (i < 0 or i >= numSites or i == side._cut) {
                return _MakeInputError("The side of a cut cell must hold other sites.");
            }
        }
    }
    JointDual dual(domain, sites, targets, volumes, HUGE_VAL, JointDual::Unknowns::ShiftsAndWeights);
    Vector z = _PackInitial(dual, shifts, weights);
    _BuildConstraints(dual, direction, sides);

    // The Newton steps keep C z constant, so the start must satisfy C z = 0.
    double violation = 0.0;
    for (const Vector& row : _rows) {
        violation = std::max(violation, std::fabs(row.dot(z)));
    }
    if (violation > 1e-12 * domain.GetVolume() * (1.0 + z.lpNorm<Eigen::Infinity>())) {
        return _MakeInputError("The initial shifts do not satisfy the linear constraints.");
    }

    _gauge = Gauge();
    Result res = _Solve(dual, z);
    if (res.HasConverged()) {
        // The gauge of the correction from the undamped bordered Newton step
        // at the solution (eq:flow-constrained-elimination-gauge), once per
        // solve: unlike the least-squares multipliers of the iterations, it
        // carries the second-order term of the step, so that the correction
        // keeps the momentum to the accuracy of the step rather than to the
        // tolerance of the solve.
        Vector g;
        dual.GetGradient(z, res._diagram, g);
        Eigen::VectorXd nu;
        if (_ComputeUndampedMultipliers(_rows, g, nu)) {
            const int dim = dual.GetDimension();
            _gauge._dilation = nu[0];
            _gauge._translation.assign(nu.data() + 1, nu.data() + 1 + dim);
            _gauge._arms.assign(nu.data() + 1 + dim, nu.data() + nu.size());
        } else {
            res._status = Status::FactorizationFailed;
        }
    }
    dual.Unpack(z, shifts, weights);
    gauge = _gauge;
    return res;
}

DualSolver::Result
VolumeSolver::Solve(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& volumes,
    std::vector<double>& weights)
{
    std::string err;
    if (not _CheckDualInputs(domain, sites, sites, volumes, err)) {
        return _MakeInputError(err);
    }
    JointDual dual(domain, sites, sites, volumes, HUGE_VAL, JointDual::Unknowns::Weights);
    std::vector<double> shifts;
    Vector z = _PackInitial(dual, shifts, weights);
    Result res = _Solve(dual, z, true);
    dual.Unpack(z, shifts, weights);

    const double D = domain.GetVolume();
    const double sum = std::accumulate(volumes.begin(), volumes.end(), 0.0);
    if (std::fabs(sum - D) > 1e-6 * D) {
        _AppendMessage(res._warning, "Target volumes do not sum to the domain volume; no exact solution exists.");
    }
    return res;
}

} // namespace PP
