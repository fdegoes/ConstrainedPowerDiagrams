#ifndef PP_CORE_CELLMOMENTS_H
#define PP_CORE_CELLMOMENTS_H

// Moments of power cells and facets, accumulated over the fragments of the
// cells clipped to the domain elements (internal to the core library).
// Shared by ComputePowerCells(), powerCells2D.h and powerCells3D.h.

#include "core/powerCells.h"

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

namespace PP
{

// Largest dimension of the fixed-size moment arrays (cellBackend.h checks it).
constexpr int kMaxDimension = 3;

// Moments of a measure (a cell, or a facet) in dim dimensions, about a
// reference point: mass ∫ dx, first moment ∫ r dx and second moment
// ∫ r r^T dx with r = x - reference.
struct Moments
{
    explicit Moments(int dim = 0) : _dim(dim) { }

    // Adds a simplex with k + 1 vertices (relative to the reference) and
    // measure measure (signed for cells, positive for facets):
    //   ∫ r      = measure / (k + 1) * sum_v v
    //   ∫ r r^T  = measure / ((k + 1)(k + 2)) * (sum_v v v^T + S S^T),
    // with S = sum_v v.
    void AddSimplex(const double* const* v, int nbVertices, double measure)
    {
        double S[kMaxDimension] = {0.0, 0.0, 0.0};
        for (int a = 0; a < nbVertices; ++a) {
            for (int c = 0; c < _dim; ++c) {
                S[c] += v[a][c];
            }
        }
        _mass += measure;
        for (int c = 0; c < _dim; ++c) {
            _first[c] += measure / nbVertices * S[c];
        }
        const double k = measure / (nbVertices * (nbVertices + 1));
        for (int r = 0; r < _dim; ++r) {
            for (int c = 0; c < _dim; ++c) {
                double vv = S[r] * S[c];
                for (int a = 0; a < nbVertices; ++a) {
                    vv += v[a][r] * v[a][c];
                }
                _second[_dim * r + c] += k * vv;
            }
        }
    }

    void Add(const Moments& m, double sign)
    {
        _mass += sign * m._mass;
        for (int c = 0; c < _dim; ++c) {
            _first[c] += sign * m._first[c];
        }
        for (int k = 0; k < _dim * _dim; ++k) {
            _second[k] += sign * m._second[k];
        }
    }

    // Trace of the second moment: ∫ |r|^2.
    double GetTraceSecond() const
    {
        double t = 0.0;
        for (int c = 0; c < _dim; ++c) {
            t += _second[_dim * c + c];
        }
        return t;
    }

    // Centroid, given the reference point.
    std::vector<double> GetCentroid(const double* reference) const
    {
        std::vector<double> c(reference, reference + _dim);
        if (_mass > 0.0) {
            for (int k = 0; k < _dim; ++k) {
                c[k] += _first[k] / _mass;
            }
        }
        return c;
    }

    // Second moment about the centroid: ∫ r r^T - (∫ r)(∫ r)^T / mass.
    std::vector<double> GetCentralSecond() const
    {
        std::vector<double> m(_second, _second + _dim * _dim);
        if (_mass > 0.0) {
            for (int r = 0; r < _dim; ++r) {
                for (int c = 0; c < _dim; ++c) {
                    m[_dim * r + c] -= _first[r] * _first[c] / _mass;
                }
            }
        }
        return m;
    }

    int    _dim;
    double _mass = 0.0;
    double _first[kMaxDimension] = {0.0, 0.0, 0.0};
    double _second[kMaxDimension * kMaxDimension] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
};

// Point p relative to reference, first dim coordinates.
inline void
MakeRelative(const double* p, const double* reference, int dim, double* r)
{
    for (int c = 0; c < dim; ++c) {
        r[c] = p[c] - reference[c];
    }
}

// Accumulates cell and facet moments over the cell fragments.
class Accumulator
{
public:
    Accumulator(int dim, const std::vector<double>& sites)
        : _dim(dim), _n(int(sites.size() / dim)), _sites(sites), _cells(_n, Moments(dim))
    { }

    int GetNumSites() const { return _n; }
    const double* GetSite(int i) const { return &_sites[_dim * i]; }

    // Adds a fragment of cell i (moments about site i). The fragment
    // orientation is arbitrary: its signed mass fixes the sign.
    void AddCellFragment(int i, const Moments& m)
    {
        _cells[i].Add(m, m._mass < 0.0 ? -1.0 : 1.0);
    }

    // Adds a piece of the facet between sites i and j (moments about
    // site i). Each pair is recorded once, from the smaller index; the
    // pieces of a pair are summed by Fill().
    void AddFacetPiece(int i, int j, const Moments& m)
    {
        if (j <= i or j >= _n or not (m._mass > 0.0)) {
            return;
        }
        _facets.emplace_back(_MakeKey(i, j), m);
    }

    // Adds the moments of other (the fragments of another part of the
    // domain).
    void Merge(const Accumulator& other)
    {
        for (int i = 0; i < _n; ++i) {
            _cells[i].Add(other._cells[i], 1.0);
        }
        _facets.insert(_facets.end(), other._facets.begin(), other._facets.end());
    }

    void Fill(PowerCells& out) const
    {
        out._dim = _dim;
        out._volumes.resize(_n);
        out._secondMoments.resize(_n);
        out._centroids.clear();
        for (int i = 0; i < _n; ++i) {
            out._volumes[i] = _cells[i]._mass;
            out._secondMoments[i] = _cells[i].GetTraceSecond();
            const std::vector<double> c = _cells[i].GetCentroid(GetSite(i));
            out._centroids.insert(out._centroids.end(), c.begin(), c.end());
        }
        // The pieces of each pair, in (i, j) order (a stable sort, so that
        // the sums do not depend on the sorting).
        std::vector<std::pair<uint64_t, Moments>> pieces = _facets;
        std::stable_sort(pieces.begin(), pieces.end(),
                         [](const auto& a, const auto& b) { return a.first < b.first; });
        out._facets.clear();
        for (size_t k = 0; k < pieces.size();) {
            Moments m = pieces[k].second;
            size_t l = k + 1;
            for (; l < pieces.size() and pieces[l].first == pieces[k].first; ++l) {
                m.Add(pieces[l].second, 1.0);
            }
            PowerCells::Facet F;
            F._i = int(pieces[k].first >> 32);
            F._j = int(pieces[k].first & 0xffffffffu);
            F._area = m._mass;
            F._centroid = m.GetCentroid(GetSite(F._i));
            F._secondMoment = m.GetCentralSecond();
            out._facets.push_back(std::move(F));
            k = l;
        }
    }

private:
    static uint64_t _MakeKey(int i, int j) { return (uint64_t(i) << 32) | uint64_t(j); }

    int                                   _dim;
    int                                   _n;
    const std::vector<double>&            _sites;
    std::vector<Moments>                  _cells;
    std::vector<std::pair<uint64_t, Moments>> _facets;  // pieces
};

} // namespace PP

#endif // PP_CORE_CELLMOMENTS_H
