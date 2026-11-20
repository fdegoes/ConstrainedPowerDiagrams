#ifndef PP_CORE_GMINTEGRATOR_H
#define PP_CORE_GMINTEGRATOR_H

#include "core/integrator.h"

namespace PP
{

// Time step of the scheme of Gallouet & Merigot (2018), the limit alpha -> 0
// of the softened integrator (latex/flow.tex, sec:flow-discussion): the
// particles are the sites of a power diagram with the prescribed volumes and
// no shifts, pulled towards the centroids of their cells.
//   w       = argmax F at q^k (VolumeSolver), warm-started from lambda^k
//   v^{k+1} = v^k + kappa dt (c - q^k)
//   q^{k+1} = q^k + dt v^{k+1}
//   p^{k+1} = q^k,  lambda^{k+1} = w
// with kappa = IntegratorOptions::_stiffness: 1 is the limit alpha -> 0 of the
// softened integrator, 1 + alpha the strength of its centroid residual
// (latex/flow.tex, Link to Gallouet & Merigot), 1 / eps^2 the stiffness of
// Gallouet & Merigot.
// Initialize() solves the weights for the volumes at the positions.
class GMIntegrator : public Integrator
{
public:
    using Integrator::Integrator;

    const char* GetName() const override { return "GM"; }
    Result Initialize(ParticleState& state) const override;
    Result Step(ParticleState& state, double dt) const override;
    bool Measure(const ParticleState& state, FlowMeasures& measures, std::string& err) const override;
    double GetPotential(double quantization, double centroidResidual) const override
    {
        return _options._stiffness * (quantization + centroidResidual);
    }

private:
    // The volume solve at the positions of cur, warm-started from its
    // diagram: weights on return.
    bool _Solve(const ParticleState& cur, std::vector<double>& weights, Result& res) const;
};

} // namespace PP

#endif // PP_CORE_GMINTEGRATOR_H
