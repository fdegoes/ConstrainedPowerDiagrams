#ifndef PP_CORE_CONSTRAINEDINTEGRATOR_H
#define PP_CORE_CONSTRAINEDINTEGRATOR_H

#include "core/integrator.h"

namespace PP
{

// Time step of the constrained integrator (latex/flow.tex,
// sec:flow-constrained; latex/appendix-algorithms.tex, Algorithm
// "flow-constrained"): the particles are the centroids of the power
// diagram of the generators, with prescribed volumes.
//   qtilde  = (1 + dt^2) q^k + dt v^k
//   (c_B, B_B) = the cut cells of the dual graph of p^k and their sides
//   (u, w, s, t, s_B) = ConstrainedJointSolver with sites q^k, targets
//                  qtilde, direction p^k and the sides (c_B, B_B), from u = 0
//                  and the weights of F at q^k (MakeJointStart with s = 0:
//                  the shifts must satisfy the constraints at the start)
//   q^{k+1}_i = qtilde_i - s p^k_i - t - sum_{B_B containing cell i} s_B (p^k_i - p^k_{c_B})
//   v^{k+1} = (q^{k+1} - q^k) / dt
//   p^{k+1} = q^k + u
//   lambda^{k+1} = w + |u|^2 + 2 u^T (q^k - qtilde)
// Initialize() moves the particles to a centroidal diagram (p = q).
class ConstrainedIntegrator : public Integrator
{
public:
    using Integrator::Integrator;

    const char* GetName() const override { return "Constrained"; }
    Result Initialize(ParticleState& state) const override;
    Result Step(ParticleState& state, double dt) const override;
    double GetPotential(double quantization, double centroidResidual) const override
    {
        return quantization + centroidResidual;
    }

    // The gauge (dilation s, translation t, branch dilations s_B) of the last
    // solve of Step().
    struct StepResult : Result
    {
        ConstrainedJointSolver::Gauge _gauge;
    };
    StepResult StepWithGauge(ParticleState& state, double dt) const;

private:
    // The cut cells of the dual graph of a diagram and their sides, from its
    // facets of positive measure between distinct generators.
    static std::vector<CutCellSide> _FindCutCellSides(const ShiftedDiagram& diagram);
    // Projects the velocities of state onto the tangent space of the
    // constraint at the diagram of its generators: removes their
    // Vbar-weighted components along the generators p_i, the d translations
    // and the localized homotheties p_i - p_c on each side (latex/flow.tex,
    // the RATTLE paragraph). False if the projection system is singular.
    bool _ProjectVelocities(const ShiftedDiagram& diagram, ParticleState& state) const;
};

} // namespace PP

#endif // PP_CORE_CONSTRAINEDINTEGRATOR_H
