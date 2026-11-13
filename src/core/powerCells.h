#ifndef PP_CORE_POWERCELLS_H
#define PP_CORE_POWERCELLS_H

#include <string>
#include <vector>

namespace PP
{

class Domain;

// Integrated quantities of the power cells V_i = cell_i ∩ domain.
struct PowerCells
{
    // Facet A_ij shared by two neighbor cells (a segment in 2D, a polygon
    // in 3D).
    struct Facet
    {
        int                 _i = 0;     // first site, i < j
        int                 _j = 0;     // second site
        double              _area = 0.0;  // measure A_ij (length in 2D)
        std::vector<double> _centroid;    // c_ij = 1/A_ij ∫_{A_ij} x dA, dim values
        // Central second moment ∫_{A_ij} (x - c_ij)(x - c_ij)^T dA,
        // dim x dim row-major. With the centroid it gives every facet
        // integral of the centroid-constraint Hessians, e.g.
        //   ∫ (x - a)(x - b)^T = secondMoment + A (c - a)(c - b)^T.
        std::vector<double> _secondMoment;
    };

    int                 _dim = 0;
    std::vector<double> _volumes;        // V_i (area in 2D)
    std::vector<double> _centroids;      // c_i, dim per site (site position if V_i = 0)
    std::vector<double> _secondMoments;  // M_i = ∫_{V_i} |x - q_i|^2 dx
    std::vector<Facet>  _facets;         // pairs of neighbor cells
    // False when the cells could not be clipped (an empty cell would seed
    // geogram's clipping walk): all cells are then reported empty.
    bool                _complete = true;

    int GetNumEmptyCells() const;
};

// Computes the power cells of weighted sites, clipped to domain, and
// their volumes, centroids, second moments and shared facets. The power
// distance of x to site i is |x - sites_i|^2 - weights_i.
//
// sites   dim doubles per site (dim = domain.GetDimension())
// weights 1 double per site
bool ComputePowerCells(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& weights,
    PowerCells& out,
    std::string& err);

} // namespace PP

#endif // PP_CORE_POWERCELLS_H
