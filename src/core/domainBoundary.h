#ifndef PP_CORE_DOMAINBOUNDARY_H
#define PP_CORE_DOMAINBOUNDARY_H

#include <cstddef>
#include <vector>

namespace PP
{

class Domain;

// The boundary of a domain as a level set, for the advection of Power
// Particles (the LevelSet of the 2015 code): the signed distance to the
// boundary (negative inside) and the closest boundary point. The boundary is
// made of the element faces with no neighbor: segments in 2D, triangles in
// 3D. Queries are brute force over the boundary faces.
class DomainBoundary
{
public:
    DomainBoundary() = default;
    explicit DomainBoundary(const Domain& domain);

    int GetDimension() const { return _dim; }
    size_t GetNumFaces() const { return _faces.size() / size_t(_dim); }

    // Signed distance of x to the boundary: < 0 inside the domain, > 0
    // outside.
    double Evaluate(const double* x) const;

    // Closest point of the boundary to x.
    void Project(const double* x, double* closest) const;

private:
    // Inside the domain: x lies in some element.
    bool _IsInside(const double* x) const;
    // Squared distance from x to boundary face f, and its closest point.
    double _FaceDistance2(size_t f, const double* x, double* closest) const;

    int                 _dim = 0;
    std::vector<double> _points;    // dim per vertex
    std::vector<int>    _elements;  // dim + 1 per element
    std::vector<int>    _faces;     // dim vertices per boundary face
};

} // namespace PP

#endif // PP_CORE_DOMAINBOUNDARY_H
