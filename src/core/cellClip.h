#ifndef PP_CORE_CELLCLIP_H
#define PP_CORE_CELLCLIP_H

// Shared helpers of the 2D and 3D power cells (powerCells2D.h,
// powerCells3D.h), internal to the core library: a uniform grid of the
// bounding boxes of the domain elements, and a loop over the sites in
// chunks, one per thread.

#include "core/cellMoments.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <thread>
#include <vector>

namespace PP
{

// Uniform grid over the bounding boxes of items in dim dimensions (2 or 3),
// with about one item per grid cell.
template <int dim>
class BoxGrid
{
public:
    // lo, hi: dim values per item.
    BoxGrid(const std::vector<double>& lo, const std::vector<double>& hi) : _itemLo(lo), _itemHi(hi)
    {
        const size_t n = lo.size() / dim;
        for (int c = 0; c < dim; ++c) {
            _lo[c] = HUGE_VAL;
            _hi[c] = -HUGE_VAL;
        }
        for (size_t k = 0; k < n; ++k) {
            for (int c = 0; c < dim; ++c) {
                _lo[c] = std::min(_lo[c], lo[dim * k + c]);
                _hi[c] = std::max(_hi[c], hi[dim * k + c]);
            }
        }
        _res = std::max(1, int(std::pow(double(n), 1.0 / dim)));
        size_t cells = 1;
        for (int c = 0; c < dim; ++c) {
            cells *= size_t(_res);
        }
        _grid.assign(cells, {});
        for (size_t k = 0; k < n; ++k) {
            int r0[dim], r1[dim];
            _GetRange(&lo[dim * k], &hi[dim * k], r0, r1);
            _ForEachCell(r0, r1, [&](size_t g) { _grid[g].push_back(int(k)); });
        }
    }

    // The items whose bounding box overlaps [lo, hi], each once: stamp holds
    // one entry per item (initially not tag), marked with tag.
    void Find(const double* lo, const double* hi, int tag, std::vector<int>& stamp, std::vector<int>& out) const
    {
        out.clear();
        int r0[dim], r1[dim];
        if (not _GetRange(lo, hi, r0, r1)) {
            return;
        }
        _ForEachCell(r0, r1, [&](size_t g) {
            for (int k : _grid[g]) {
                if (stamp[size_t(k)] == tag) {
                    continue;
                }
                bool overlaps = true;
                for (int c = 0; c < dim; ++c) {
                    overlaps = overlaps and _itemHi[dim * size_t(k) + c] >= lo[c] and _itemLo[dim * size_t(k) + c] <= hi[c];
                }
                if (overlaps) {
                    stamp[size_t(k)] = tag;
                    out.push_back(k);
                }
            }
        });
    }

private:
    // Grid cells overlapping [lo, hi]; false if none.
    bool _GetRange(const double* lo, const double* hi, int* r0, int* r1) const
    {
        for (int c = 0; c < dim; ++c) {
            if (hi[c] < _lo[c] or lo[c] > _hi[c]) {
                return false;
            }
            const double w = _hi[c] - _lo[c];
            auto cell = [&](double x) {
                const int k = w > 0.0 ? int((x - _lo[c]) / w * _res) : 0;
                return std::max(0, std::min(_res - 1, k));
            };
            r0[c] = cell(lo[c]);
            r1[c] = cell(hi[c]);
        }
        return true;
    }

    template <class F>
    void _ForEachCell(const int* r0, const int* r1, const F& f) const
    {
        int k[dim];
        std::copy(r0, r0 + dim, k);
        for (;;) {
            size_t g = 0;
            for (int c = dim - 1; c >= 0; --c) {
                g = g * size_t(_res) + size_t(k[c]);
            }
            f(g);
            int c = 0;
            for (; c < dim and k[c] == r1[c]; ++c) {
                k[c] = r0[c];
            }
            if (c == dim) {
                return;
            }
            ++k[c];
        }
    }

    std::vector<double>           _itemLo, _itemHi;
    double                        _lo[dim], _hi[dim];
    int                           _res = 1;
    std::vector<std::vector<int>> _grid;
};

// Runs work(first, last, acc) over the sites [0, n) of acc, in one chunk per
// thread with parallel (each with its own accumulator, merged in chunk
// order, so that the result does not depend on the threads), else in one
// chunk. False if a chunk returns false.
template <class Work>
bool
ForEachSiteChunk(int n, bool parallel, const std::vector<double>& sites, int dim, Accumulator& acc, const Work& work)
{
    if (not parallel) {
        return work(0, n, acc);
    }
    const int chunks = int(std::max(1u, std::min(std::thread::hardware_concurrency(), 64u)));
    std::vector<Accumulator> parts(size_t(chunks) - 1, Accumulator(dim, sites));
    std::vector<char> ok(size_t(chunks), 1);
    std::vector<std::thread> threads;
    auto bound = [&](int c) { return int(int64_t(n) * c / chunks); };
    for (int c = 1; c < chunks; ++c) {
        threads.emplace_back([&, c]() { ok[size_t(c)] = work(bound(c), bound(c + 1), parts[size_t(c) - 1]); });
    }
    ok[0] = work(0, bound(1), acc);
    for (std::thread& t : threads) {
        t.join();
    }
    for (const Accumulator& part : parts) {
        acc.Merge(part);
    }
    return std::all_of(ok.begin(), ok.end(), [](char x) { return x != 0; });
}

} // namespace PP

#endif // PP_CORE_CELLCLIP_H
