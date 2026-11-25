#ifndef PP_CORE_INTEGRATOR_H
#define PP_CORE_INTEGRATOR_H

#include "core/dualSolver.h"

#include <memory>
#include <string>
#include <vector>

namespace PP
{

class Domain;

// State of a particle fluid: positions q and velocities v of the particles,
// and the power diagram of the step as shifts u and weights w, i.e. the
// generators p = q + u and power weights lambda = w + |u|^2 (the convention
// of the SOPs, which read P + u and w + |u|^2). dim doubles per particle for
// the vectors, one per particle for the scalars. Empty velocities, shifts or
// weights are taken as zeros.
struct ParticleState
{
    std::vector<double> _positions;   // q
    std::vector<double> _velocities;  // v
    std::vector<double> _shifts;      // u
    std::vector<double> _weights;     // w
    // Power Particles only: the pressures of the previous step (unused by
    // the other integrators).
    std::vector<double> _pressures;
};

// Energy and momentum of a particle fluid, to monitor a simulation over
// time. Energies follow the Lagrangian of latex/flow.tex
// (eq:flow-lagrangian), L = sum_i Vbar_i |v_i|^2 - U(q), i.e. masses 2 Vbar_i
// (halve them for the usual 1/2 m |v|^2 with m = Vbar):
//   kinetic           K = sum_i Vbar_i |v_i|^2
//   quantization      Q = sum_i ∫_{V_i} |x - c_i|^2 dx    (cells about their centroids)
//   centroid residual R = sum_i Vbar_i |c_i - q_i|^2
//   gravity           G = -2 sum_i Vbar_i g . q_i        (so that q'' = g)
//   potential         U, the integrator's potential (Integrator::GetPotential)
//   total             E = K + U + G, conserved by the variational schemes up to
//                     their time discretization
// of the power diagram of the state, and the momenta (masses Vbar_i):
//   momentum          P = sum_i Vbar_i v_i                 (dim values)
//   angular momentum  sum_i Vbar_i (q_i - o) x v_i about the domain centroid
//                     o (1 value in 2D, 3 in 3D)
struct FlowMeasures
{
    double              _kinetic = 0.0;
    double              _quantization = 0.0;
    double              _centroidResidual = 0.0;
    double              _gravity = 0.0;
    double              _potential = 0.0;
    double              _total = 0.0;
    std::vector<double> _momentum;
    std::vector<double> _angularMomentum;

    // One CSV line per time step, for plotting: the header, then the values
    // at time t.
    static std::string GetCSVHeader(int dim);
    std::string GetCSVRow(double t) const;
};

// Options of the specific integrators (MakeIntegrator); each uses its own.
struct IntegratorOptions
{
    // Softened: stiffness alpha > 0 of the centroid penalty.
    double _alpha = 1.0;
    // GM: scale kappa of the restoring force. 1 is the limit alpha -> 0 of
    // the softened integrator, 1 + alpha the strength of its centroid
    // residual (latex/flow.tex, Link to Gallouet & Merigot), 1 / eps^2 the
    // stiffness of Gallouet & Merigot (2018).
    // Constrained: scale kappa of the potential (the Lagrangian sum V |dq/dt|^2 - kappa Phi), whose
    // force kappa (q - p) pulls the particles towards their generators; 1 is latex/flow.tex.
    double _stiffness = 1.0;
    // Power Particles: fraction beta of the move to the centroid in the
    // advection (0.5 in the 2015 code), Newton steps of the weights
    // (cap_max_iters) and their tolerance (cap_tolerance), and the diagonal
    // regularization eps of the Laplacian and Hessian (< 0: the 2015 values,
    // 1e-8 in 2D and 1e-6 in 3D).
    double _centroidBlend = 0.5;
    int    _weightIterations = 2;
    double _weightTolerance = 0.0;
    double _regularization = -1.0;
    // Constrained: projection of the velocities onto the tangent space of
    // the holonomic constraint (latex/flow.tex, the RATTLE paragraph of
    // sec:flow-discussion), i.e. Vbar-orthogonal to the generators, the
    // translations and the localized homotheties of the diagram: never
    // (SHAKE), on the initial velocities only, or after every step (RATTLE).
    enum class VelocityProjection
    {
        None,
        Initial,
        EveryStep
    };
    VelocityProjection _velocityProjection = VelocityProjection::None;
};

// One time step of a particle fluid whose particles carry the prescribed
// volumes Vbar (their masses, at unit density). The integrators work in 2D
// and 3D (the dimension of the domain). The fluid fills the whole domain:
// the volumes sum to the domain volume and every cell is fluid (free
// surfaces, e.g. with air particles, are left for later). Step() advances the state by dt and
// leaves it unchanged if the step fails; at the end of a step the state
// holds the power diagram of the step re-expressed at the new particles
// (same p and lambda, u = p - q, w = lambda - |u|^2), and the result its
// cells.
class Integrator
{
public:
    struct Parameters
    {
        // The integrators solve to a looser tolerance than the static
        // solves (DualSolver::Parameters, 1e-8), the same for the volumes
        // and the shifts (their updates use the shifts directly).
        Parameters() { _solver._volumeTolerance = _solver._shiftTolerance = 1e-4; }

        // Thresholds, iterations and regularization of the dual solves.
        DualSolver::Parameters _solver;
        // Gravity: the acceleration of the external force (dim values, e.g.
        // {0, -9.8} in 2D; empty: none), applied as a velocity kick
        // v += dt g at the start of each step, before the solve.
        std::vector<double>    _gravity;
        // Initialize(): convergence of the Lloyd iterations of the
        // centroidal diagrams (CentroidalDiagram::Parameters).
        double _lloydTolerance = 1e-9;
        int    _lloydMaxIterations = 1000;
    };

    struct Result
    {
        bool               _converged = false;
        std::string        _status;  // of the last solve
        std::string        _error;
        DualSolver::Result _solve;  // the last dual solve (empty for Power Particles)
        int                _iterations = 0;  // Newton iterations of the step
        // Of the power diagram of the state on return: centroids (dim per
        // particle) and volumes of the cells.
        std::vector<double> _centroids;
        std::vector<double> _cellVolumes;
        double _maxVolumeError = 0.0;  // max_i |V_i - Vbar_i| / Vbar_i
        double _maxCentroidGap = 0.0;  // max_i |q_i - c_i|
        int    _numProjected = 0;      // Power Particles: particles pushed back into the domain
        // Energy and momentum of the state on return. For the softened and GM
        // integrators, its potential is measured on the diagram solved at the
        // start of the step, i.e. at the previous positions: use
        // _startMeasures, or Measure(), to pair the kinetic and potential
        // energies of one time level.
        FlowMeasures _measures;
        // Energy and momentum of the state at the start of the step (before
        // the gravity kick), on the diagram the step solves for it, so that
        // the kinetic and potential energies are those of the same time
        // level, at no extra solve. Empty for Power Particles.
        FlowMeasures _startMeasures;

        bool HasConverged() const { return _converged; }
    };

    // domain must outlive the integrator; volumes (one per particle, summing
    // to the domain volume) are copied.
    Integrator(
        const Domain& domain,
        std::vector<double> volumes,
        const Parameters& params,
        const IntegratorOptions& options);
    virtual ~Integrator() = default;

    Integrator(const Integrator&) = delete;
    Integrator& operator=(const Integrator&) = delete;

    const std::vector<double>& GetVolumes() const { return _volumes; }
    Parameters& GetParameters() { return _params; }
    const Parameters& GetParameters() const { return _params; }

    virtual const char* GetName() const = 0;

    // Makes the power diagram of state consistent with this integrator at
    // the start of a simulation (the positions may move, the velocities are
    // kept); state is unchanged if it fails.
    virtual Result Initialize(ParticleState& state) const = 0;

    // Advances state by dt; state is unchanged if the step fails.
    virtual Result Step(ParticleState& state, double dt) const = 0;

    // The potential U of this integrator from the quantization error Q and
    // the centroid residual R of a diagram (FlowMeasures):
    //   Constrained:     Q + R            (R = 0 at the solution)
    //   Softened:        Q + (1 + alpha) R   (eq:flow-softened-potential)
    //   GM:              kappa (Q + R)
    //   Power Particles: 0 (the pressure is not a potential force)
    virtual double GetPotential(double quantization, double centroidResidual) const = 0;

    // Energy and momentum of state, on its power diagram: built from the
    // state for the constrained integrator and Power Particles, and solved
    // at its positions for the softened and GM integrators (one dual solve),
    // as at the start of their steps.
    virtual bool Measure(const ParticleState& state, FlowMeasures& measures, std::string& err) const;

protected:
    int _GetDimension() const;

    // Copies in to out with the empty fields filled with zeros; false (and
    // res._error) if the sizes do not match the volumes.
    bool _Prepare(const ParticleState& in, ParticleState& out, Result& res) const;

    // v += dt g.
    void _ApplyGravity(std::vector<double>& velocities, double dt) const;

    // Generators p = q + u and power weights lambda = w + |u|^2 of state.
    void _GetDiagram(const ParticleState& state, std::vector<double>& p, std::vector<double>& lambda) const;
    // Sets the shifts and weights of state so that its diagram is (p, lambda)
    // at its positions.
    void _SetDiagram(ParticleState& state, const std::vector<double>& p, const std::vector<double>& lambda) const;

    // Power diagram of state (targets at the positions, so that its power
    // weights are w + |u|^2).
    bool _BuildDiagram(const ParticleState& state, ShiftedDiagram& diagram, std::string& err) const;

    // Fills the cells and diagnostics of res from the diagram of state.
    void _Finish(const ShiftedDiagram& diagram, const ParticleState& state, Result& res) const;

    // Energy and momentum of state, whose power diagram is diagram.
    FlowMeasures _Measure(const ShiftedDiagram& diagram, const ParticleState& state) const;

    // Initial state from the centroidal diagram of the positions: q = p =
    // the centroidal sites, u = 0, w its weights.
    bool _MakeCentroidal(ParticleState& state, Result& res) const;

    const Domain&       _domain;
    std::vector<double> _volumes;
    Parameters          _params;
    IntegratorOptions   _options;
};

enum class IntegratorKind
{
    Constrained,     // ConstrainedIntegrator
    Softened,        // SoftenedIntegrator
    GM,              // GMIntegrator
    PowerParticles,  // PowerParticlesIntegrator
};

// The integrator of the given kind (nullptr for an unknown kind).
std::unique_ptr<Integrator> MakeIntegrator(
    IntegratorKind kind,
    const Domain& domain,
    std::vector<double> volumes,
    const Integrator::Parameters& params = Integrator::Parameters(),
    const IntegratorOptions& options = IntegratorOptions());

// Time step of the CFL condition of the 2015 Power Particles code:
// cfl * r_min / max_i max_k |v_ik|, with r_min the radius of the ball of the
// smallest volume; 0 if all velocities are zero.
double ComputeCFLStep(const std::vector<double>& velocities, const std::vector<double>& volumes, int dim, double cfl);

} // namespace PP

#endif // PP_CORE_INTEGRATOR_H
