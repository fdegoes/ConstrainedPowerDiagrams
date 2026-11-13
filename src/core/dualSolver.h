#ifndef PP_CORE_DUALSOLVER_H
#define PP_CORE_DUALSOLVER_H

#include "core/jointDual.h"
#include "core/lineSearch.h"
#include "core/shiftedDiagram.h"
#include "core/sparseCholesky.h"

#include <Eigen/Core>

#include <cmath>
#include <functional>
#include <iosfwd>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace PP
{

class Domain;

// Newton's method on a concave dual JointDual (latex/appendix-algorithms.tex,
// Algorithm "newton"): from the multipliers z, repeat
//   g = gradient at z
//   N = -Hessian at z + eta M (JointDual::GetNegatedHessian(), with the
//       distance clamp and the pair regularization of the parameters)
//   step from N and g, rebuilding the power diagram at each trial point
// until the gradient residuals are below the tolerances.
//
// The step is Levenberg-Marquardt: (N + mu D^2) p = g, with Marquardt
// scaling D^2 = diag(N), bordered by linear rows C (C p = 0) for
// ConstrainedJointSolver. A trial is accepted when no cell is empty and the
// actual increase of the objective is at least acceptRatio times the
// increase g^T p - p^T N p / 2 of the model (when that increase is below
// the round-off of the objective, when the residual decreases instead). mu
// decreases after accepted steps and grows after rejected ones (Nielsen's
// update), so a rejection only refactorizes, without rebuilding the
// Hessian; mu persists across iterations. VolumeSolver (F) instead takes the
// Newton step with the backtracking line search of Algorithm "line-search".
//
// Termination only depends on the gradient residuals: the solve converges
// when both the volume and the shift residuals are at most their
// tolerances, and fails otherwise (no accepted step, max iterations,
// factorization).
class DualSolver
{
public:
    using Vector = JointDual::Vector;
    using SparseMatrix = JointDual::SparseMatrix;

    struct Result;

    enum class Status
    {
        Converged,            // residuals <= tolerances (the only success)
        MaxIterations,        // maxIterations Newton steps taken
        LineSearchFailed,     // no increase after maxHalvings (F)
        DampingFailed,        // no accepted step after maxRejections, or mu
                              // above maxDamping
        FactorizationFailed,  // Cholesky of N failed
        InvalidInput,         // bad sizes or values, or geogram failure
    };

    struct Damping
    {
        double _initialDamping = 1e-3;  // mu at the first iteration
        double _minDamping = 0.0;       // mu never drops below this
        double _maxDamping = 1e12;      // fail when mu exceeds this
        double _acceptRatio = 1e-4;     // accept the step if rho > acceptRatio
        int    _maxRejections = 60;     // per Newton iteration
    };

    struct Parameters
    {
        int    _maxIterations = 100;
        // Converged when both residuals (JointDual::GetResiduals() of the
        // gradient, or of the constrained stationarity condition) are at
        // most these: every |V_i - Vbar_i| <= volumeTolerance Vbar_i, and
        // every |g_u_i| <= 2 shiftTolerance Vbar_i^((d+1)/d), i.e. a centroid
        // error of about shiftTolerance times the cell size.
        double _volumeTolerance = 1e-8;
        double _shiftTolerance = 1e-8;
        // eta, the single constant of the regularizations of N, which change
        // N but not the gradient, hence not the solution:
        //  - N = -Hessian + eta M, with the per-block mass M of
        //    JointDual::GetMassDiagonal() (units L^-2);
        //  - JointDual::SetHessianMinDistanceRatio(): generator distances in
        //    the Hessian facet coefficients (1 / d_ij) are bounded below by
        //    eta times the cell size (|domain| / n)^(1/d), so that close or merging
        //    generators cannot make N singular;
        // Each part is disabled by setting its own parameter below to 0
        // (the clamp: minDistance; the pair term: radius).
        double _regularization = 1e-3;
        bool   _hessianClamp = true;
        // JointDual::SetPairRegularization(): radius R of the pair
        // regularization as a ratio of the cell size (|domain| / n)^(1/d)
        // (R = 0 disables it), and its maximum rhoMax_ij as a multiple of the
        // larger curvature of the two sites along the pair. R is scaled at every
        // iteration by min(1, r / r0), with r the shift residual and r0 its
        // value at the start, so that the regularization fades out as the
        // solve converges.
        double _pairRegularizationRadius = 0.3;
        double _pairRegularizationScale = 2.0;
        // JointSolver restarts: when a solve stops without converging (max
        // iterations or no accepted step), solve again from the start, up to
        // pairRestarts times (0: no restart), with the radius of the pair
        // regularization multiplied by pairRestartGrowth each time.
        int    _pairRestarts = 3;
        double _pairRestartGrowth = 2.0;
        LineSearch::Parameters _lineSearch;  // VolumeSolver
        Damping                _damping;     // the other solvers
        // Optional progress log (one line per iteration).
        std::ostream* _log = nullptr;
        // Optional observer of every iterate (the start and each accepted
        // step), called with the multipliers z and the result so far (its
        // diagram and residuals are those of z).
        std::function<void(const JointDual&, const Vector&, const Result&)> _onIterate;
    };

    struct Result
    {
        Status         _status = Status::InvalidInput;
        int            _iterations = 0;
        int            _numBuilds = 0;  // power diagrams built (trial points)
        // Power diagrams built by earlier attempts of the same solve
        // (JointSolver restarts), for observers that count builds.
        int            _previousBuilds = 0;
        // Solves run (1 + the JointSolver restarts) and the iterations of
        // the returned one.
        int            _attempts = 1;
        int            _attemptIterations = 0;
        double         _residual = 0.0;        // JointDual::GetResidual(), <= 1 when converged
        double         _volumeResidual = 0.0;  // max_i |g_w_i| / Vbar_i
        double         _shiftResidual = 0.0;   // max_i |g_u_i| / (2 Vbar_i^((d+1)/d))
        // Wall time (seconds) of the power diagram builds (including the
        // line search), the Hessian assembly, and the factorizations and
        // solves.
        double         _buildSeconds = 0.0;
        double         _hessianSeconds = 0.0;
        double         _factorSeconds = 0.0;
        double         _value = 0.0;   // objective at the returned multipliers
        ShiftedDiagram _diagram;       // power diagram at the returned multipliers
        std::string    _warning;
        std::string    _error;

        bool HasConverged() const { return _status == Status::Converged; }
    };

    DualSolver() = default;
    explicit DualSolver(const Parameters& params) : _params(params) { }
    virtual ~DualSolver() = default;

    Parameters& GetParameters() { return _params; }
    const Parameters& GetParameters() const { return _params; }

    static const char* GetStatusString(Status s);

protected:
    int          _buildOffset = 0;           // Result::_previousBuilds of the next solve
    using Cholesky = SparseCholesky;

    // Maximizes dual from z (the solution, or the last accepted iterate, on
    // return). useLineSearch selects the line search instead of
    // Levenberg-Marquardt (VolumeSolver).
    // The tolerances of dual are set from the parameters.
    Result _Solve(JointDual& dual, Vector& z, bool useLineSearch = false);

    // Step dz and the stationarity vector whose residuals are tested (g
    // itself for the plain step) from N (factorized by the caller with
    // _Factorize()) and the gradient g at z. Returns false if a solve
    // failed.
    virtual bool _ComputeStep(
        const JointDual& dual,
        const Vector& z,
        const ShiftedDiagram& diagram,
        const Vector& g,
        Vector& dz,
        Vector& stationarity);

    // Solves N x = b with the factorization of the current iteration.
    Vector _SolveFactorized(const Vector& b) const { return _llt.solve(b); }
    // Whether the iterations need the undamped step of N (_ComputeStep()):
    // for the line search, and for the stationarity residual of the
    // bordered step. The plain Levenberg-Marquardt step only factorizes the
    // damped matrices.
    virtual bool _UsesUndampedStep() const { return _useLineSearch; }

    // Bordered Newton step (latex/appendix-bordered.tex) with the
    // constraint rows C: solves [N C^T; C 0] [dz; 2 nu] = [g; 0] by block
    // elimination, with one solve with N per row plus one for g. On return
    // nu holds one multiplier per row and stationarity is
    // g - 2 C^T nu.
    // The rows are kept for the damped steps of the iteration.
    bool _ComputeBorderedStep(
        const JointDual& dual,
        const std::vector<Vector>& rows,
        const Vector& g,
        Vector& dz,
        Eigen::VectorXd& nu,
        Vector& stationarity);

private:
    bool _Start(const JointDual& dual, const Vector& z, Result& res);
    bool _Factorize(const JointDual& dual, const ShiftedDiagram& diagram, Result& res);
    void _Log(const Result& res, double stepLength) const;

    // Steps of _Solve() (line search for F, damped otherwise); return false to stop with
    // res._status set.
    bool _TakeLineSearchStep(const JointDual& dual, const Vector& g, const Vector& dz, Vector& z, Result& res);
    // Ratio test of the damped step for a trial z + p with predicted
    // increase predicted: accepts when no cell is empty (volume at most
    // _lineSearch._minVolume) and the actual
    // increase is at least acceptRatio times the predicted one. When the
    // predicted increase is within round-off of the objective, the trial is
    // accepted instead if it decreases the gradient residual; rho is then 1
    // or -1. Fills value with the objective at the trial.
    bool _AcceptTrial(
        const JointDual& dual,
        const Vector& zTrial,
        const ShiftedDiagram& trial,
        double predicted,
        double acceptRatio,
        const Result& res,
        double& value,
        double& rho) const;
    bool _TakeDampedStep(const JointDual& dual, const Vector& g, const Vector& dz, Vector& z, Result& res);
    // Step of (N + mu D^2) p = g with the border rows of the current
    // iteration; refactorizes. Returns false if the factorization fails.
    bool _ComputeDampedStep(const JointDual& dual, const Vector& g, double mu, Vector& p);

    double _ComputeModelIncrease(const Vector& g, const Vector& p) const;
    // Residual of a trial gradient comparable with Result::_residual: of g
    // itself, or of the stationarity condition g - 2 C^T nu with the border
    // rows of the current iteration (nu by least squares).
    double _ComputeTrialResidual(const JointDual& dual, const Vector& g) const;

    Parameters   _params;
    Cholesky     _llt;
    SparseMatrix _N;
    SparseMatrix _damped;                  // N + mu D^2, the pattern of N
    std::vector<double*> _dampedDiagonal;  // its diagonal entries
    bool         _useLineSearch = false;
    Vector       _scale2;           // D^2 = diag(N)
    double       _mu = 0.0;         // damping (Levenberg-Marquardt)
    double       _muGrowth = 2.0;   // Nielsen's nu
    bool         _resetDamping = true;  // set mu to its initial value at the next step
    // Border rows of the current step (empty for the plain step).
    std::vector<Vector> _borderRows;
};

// Maximizes the softened joint dual H_alpha (latex/contributions.tex,
// eq:regularized-dual) over the shifts and weights.
class SoftenedSolver : public DualSolver
{
public:
    using DualSolver::DualSolver;

    // sites, targets  dim doubles per site, q and cbar
    // volumes         target volumes Vbar
    // alpha           stiffness (> 0)
    // shifts, weights initial multipliers (zeros if the sizes do not match);
    //                 the solution on return
    Result Solve(
        const Domain& domain,
        const std::vector<double>& sites,
        const std::vector<double>& targets,
        const std::vector<double>& volumes,
        double alpha,
        std::vector<double>& shifts,
        std::vector<double>& weights);
};

// Maximizes the centroid dual G (latex/contributions.tex, eq:centroid-dual)
// over the shifts, with fixed weights: cells with the prescribed centroids,
// when they exist. Its Hessian is negative definite unless the pivots of
// eq:pivot are coplanar, so the plain step is used.
class CentroidSolver : public DualSolver
{
public:
    using DualSolver::DualSolver;

    // sites, targets  dim doubles per site, q and cbar
    // weights         the fixed weights w
    // volumes         positive scales of the cells (regularization and
    //                 residual only, see JointDual)
    // shifts          initial shifts (zeros if the size does not match);
    //                 the solution on return
    Result Solve(
        const Domain& domain,
        const std::vector<double>& sites,
        const std::vector<double>& targets,
        const std::vector<double>& weights,
        const std::vector<double>& volumes,
        std::vector<double>& shifts);
};

// Maximizes the joint dual H (latex/contributions.tex, eq:joint-dual): cells
// with the prescribed volumes and centroids, when they exist.
//
// The Hessian of H vanishes along the d + 2 homotheties of eq:joint-nullspace
// (and along one more mode per extra connected component of the interior
// facets, e.g. a cell with a single interior facet; notes/
// joint_nullspace_proof.md). No border is used: the Levenberg-Marquardt
// damping mu diag(N) gives the step curvature along these modes.
//
// When a solve stops on MaxIterations or DampingFailed, it restarts from the
// start up to pairRestarts times with a stronger pair regularization (see
// Parameters::_pairRestarts), and returns the converged or the best attempt;
// the counts and timers of the Result add up all attempts.
class JointSolver : public DualSolver
{
public:
    using DualSolver::DualSolver;

    // sites, targets  dim doubles per site, q and cbar
    // volumes         target volumes Vbar
    // shifts, weights initial multipliers (zeros if the sizes do not match);
    //                 the solution on return.
    //
    // Expected start: MakeJointStart() (an F solve at generators between the
    // sites and the targets), so that the start already has the prescribed
    // volumes and H only corrects the centroids.
    Result Solve(
        const Domain& domain,
        const std::vector<double>& sites,
        const std::vector<double>& targets,
        const std::vector<double>& volumes,
        std::vector<double>& shifts,
        std::vector<double>& weights);
};

// Start of the H solves: the power diagram of the generators
//   p = (1 - s) q + s cbar
// with the prescribed volumes (VolumeSolver, F, at p), written as multipliers
// of H for the sites q and targets cbar:
//   u = p - q,   w = lambda_F - |u|^2 - 2 u^T (q - cbar)
// (eq:shifted-site, so that the shifted weights are lambda_F).
//
// s = 0 starts from the sites (u = 0, w = F at q): best when the sites are
// already close to the solution, and required when the shifts must satisfy
// linear constraints at the start (ConstrainedJointSolver). s = 1 starts
// from the targets, and the H solve then no longer depends on q. In the
// static benchmarks (sites = reference generators + noise), s = 0.5 was the
// most robust: 36-40/40 solves in 2D for displacements up to the cell size,
// and 40/40 in 3D.
//
// sites, targets  dim doubles per site, q and cbar
// volumes         target volumes Vbar
// weights         initial weights of the F solve (zeros if the size does not
//                 match)
// shifts, weights on return, the start of JointSolver::Solve()
// Returns the result of the F solve.
DualSolver::Result MakeJointStart(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& targets,
    const std::vector<double>& volumes,
    double s,
    std::vector<double>& shifts,
    std::vector<double>& weights);

// Maximizes H subject to the dilation and translation constraints of the
// constrained integrator (latex/flow.tex, eq:flow-line-constraint):
//   sum_i Vbar_i u_i^T p_i = 0,   sum_i Vbar_i u_i = 0,
// for a direction p (the generators of the previous step). Each Newton step
// is the bordered solve of latex/appendix-bordered.tex: d + 2 solves with N
// and a dense (d + 1) x (d + 1) system for the gauge (s, t), the
// multipliers of the constraints (scaled by 1/2).
class ConstrainedJointSolver : public DualSolver
{
public:
    using DualSolver::DualSolver;

    struct Gauge
    {
        double              _dilation = 0.0;  // s
        std::vector<double> _translation;     // t, dim values
    };

    // direction  dim doubles per site, p
    // shifts     must satisfy the constraints (e.g. zeros, the default
    //            when the size does not match)
    // gauge      (s, t) of the last Newton step on return
    Result Solve(
        const Domain& domain,
        const std::vector<double>& sites,
        const std::vector<double>& targets,
        const std::vector<double>& volumes,
        const std::vector<double>& direction,
        std::vector<double>& shifts,
        std::vector<double>& weights,
        Gauge& gauge);

protected:
    bool _ComputeStep(
        const JointDual& dual,
        const Vector& z,
        const ShiftedDiagram& diagram,
        const Vector& g,
        Vector& dz,
        Vector& stationarity) override;
    // The bordered step sets the border rows and the gauge.
    bool _UsesUndampedStep() const override { return true; }

private:
    // Rows a, b_1, ..., b_d of the constraint matrix C
    // (eq:appendix-constraint-matrix).
    void _BuildConstraints(const JointDual& dual, const std::vector<double>& direction);

    std::vector<Vector> _rows;
    Gauge               _gauge;
};

// Maximizes the volume dual F (latex/background.tex, eq:volume-dual) over the
// weights: cells with the prescribed volumes.
class VolumeSolver : public DualSolver
{
public:
    using DualSolver::DualSolver;

    // sites    dim doubles per site
    // volumes  one positive target volume per site, summing to
    //          domain.GetVolume()
    // weights  initial weights (zeros if the size does not match); the
    //          solution on return
    Result Solve(
        const Domain& domain,
        const std::vector<double>& sites,
        const std::vector<double>& volumes,
        std::vector<double>& weights);
};

} // namespace PP

#endif // PP_CORE_DUALSOLVER_H
