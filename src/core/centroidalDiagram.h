#ifndef PP_CORE_CENTROIDALDIAGRAM_H
#define PP_CORE_CENTROIDALDIAGRAM_H

#include "core/dualSolver.h"
#include "core/integrator.h"

#include <string>
#include <vector>

namespace PP
{

class Domain;

// Centroidal power diagram with prescribed volumes, by Lloyd iterations:
// solve the weights for the volumes (VolumeSolver, warm-started), move the
// sites to the centroids of their cells, and repeat until the sites are
// the centroids. Its sites are both the particles and the generators of
// the initial state of the integrators (latex/flow.tex).
class CentroidalDiagram
{
public:
    struct Parameters
    {
        DualSolver::Parameters _volumeSolver;
        int    _maxIterations = 1000;
        // Converged when max_i |c_i - q_i| <= tolerance * diam(domain),
        // with diam the bounding-box diagonal.
        double _tolerance = 1e-8;
    };

    struct Result
    {
        bool        _converged = false;
        int         _iterations = 0;
        double      _residual = 0.0;  // max_i |c_i - q_i| / diam(domain)
        std::string _error;
    };

    CentroidalDiagram() = default;
    explicit CentroidalDiagram(const Parameters& params) : _params(params) { }

    Parameters& GetParameters() { return _params; }
    const Parameters& GetParameters() const { return _params; }

    // sites    dim doubles per site, the initial guess; the centroidal
    //          sites on return
    // volumes  target volumes, summing to domain.GetVolume()
    // weights  initial weights (zeros if the size does not match); the
    //          weights of the centroidal diagram on return
    Result Solve(
        const Domain& domain,
        const std::vector<double>& volumes,
        std::vector<double>& sites,
        std::vector<double>& weights) const;

    // Initial state of the integrators from a centroidal diagram:
    // particles at the sites, zero velocities and shifts, and the weights.
    static ParticleState MakeInitialState(const std::vector<double>& sites, const std::vector<double>& weights);

private:
    Parameters _params;
};

} // namespace PP

#endif // PP_CORE_CENTROIDALDIAGRAM_H
