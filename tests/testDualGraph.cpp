// Cut cells of the dual graph and their sides (core/dualGraph.h).

#include "core/dualGraph.h"
#include "testUtils.h"

#include <set>
#include <string>
#include <vector>

using namespace PP;
using Edges = std::vector<std::pair<int, int>>;

namespace
{

// Each side must be a connected component of the graph without its cut cell,
// i.e. connected, and linked to the other cells only through the cut cell.
bool
_IsComponentWithoutCut(int n, const Edges& edges, const CutCellSide& side)
{
    const std::set<int> S(side._side.begin(), side._side.end());
    if (S.empty() or S.count(side._cut)) {
        return false;
    }
    std::vector<std::vector<int>> adj(n);
    for (const auto& [i, j] : edges) {
        if ((S.count(i) and not S.count(j) and j != side._cut) or (S.count(j) and not S.count(i) and i != side._cut)) {
            return false;
        }
        if (S.count(i) and S.count(j)) {
            adj[i].push_back(j);
            adj[j].push_back(i);
        }
    }
    std::set<int> seen{*S.begin()};
    std::vector<int> stack{*S.begin()};
    while (not stack.empty()) {
        const int u = stack.back();
        stack.pop_back();
        for (int v : adj[u]) {
            if (seen.insert(v).second) {
                stack.push_back(v);
            }
        }
    }
    return seen.size() == S.size();
}

void
_TestGraph(const std::string& name, int n, const Edges& edges, size_t expected)
{
    const std::vector<CutCellSide> sides = FindCutCellSides(n, edges);
    PPTest::Check(sides.size() == expected,
                  name + ": " + std::to_string(sides.size()) + " sides (b - 1 = " + std::to_string(expected) + ")");
    bool ok = true;
    for (const CutCellSide& side : sides) {
        ok = ok and _IsComponentWithoutCut(n, edges, side);
    }
    PPTest::Check(ok, name + ": every side is a component without its cut cell");
}

}  // namespace

int
main()
{
    // A path of 6 nodes: 5 blocks (bridges), 4 cut cells.
    _TestGraph("path", 6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}}, 4);
    // Two triangles sharing node 2: 2 blocks.
    _TestGraph("bowtie", 5, {{0, 1}, {1, 2}, {2, 0}, {2, 3}, {3, 4}, {4, 2}}, 1);
    // A cycle and a triangulated 3 x 3 grid: a single block.
    _TestGraph("cycle", 6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}}, 0);
    Edges grid;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            const int k = 3 * r + c;
            if (c < 2) grid.emplace_back(k, k + 1);
            if (r < 2) grid.emplace_back(k, k + 3);
            if (c < 2 and r < 2) grid.emplace_back(k, k + 4);
        }
    }
    _TestGraph("grid", 9, grid, 0);
    // A star with 4 leaves: 4 blocks at the center.
    _TestGraph("star", 5, {{0, 1}, {0, 2}, {0, 3}, {0, 4}}, 3);
    // A triangle with a dangling edge at each corner, listed with duplicates
    // and a self-loop: 4 blocks.
    _TestGraph("triangle with tails", 6, {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 4}, {2, 5}, {1, 0}, {3, 3}}, 3);
    // The search rooted inside a block or at a leaf gives the same count.
    _TestGraph("path rooted at its middle", 5, {{2, 1}, {2, 3}, {1, 0}, {3, 4}}, 3);
    return PPTest::Report();
}
