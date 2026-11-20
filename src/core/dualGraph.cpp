#include "core/dualGraph.h"

#include <algorithm>

namespace PP
{

std::vector<CutCellSide>
FindCutCellSides(int numCells, const std::vector<std::pair<int, int>>& edges)
{
    // Adjacency lists without duplicates or self-loops.
    std::vector<std::vector<int>> adj(numCells);
    for (const auto& [i, j] : edges) {
        if (i != j and 0 <= i and i < numCells and 0 <= j and j < numCells) {
            adj[i].push_back(j);
            adj[j].push_back(i);
        }
    }
    for (std::vector<int>& a : adj) {
        std::sort(a.begin(), a.end());
        a.erase(std::unique(a.begin(), a.end()), a.end());
    }

    // Iterative depth-first search: discovery time, low point, and the last
    // discovery time in the subtree, so that the subtree of v is the range
    // order[tin[v] .. tout[v]] of the nodes in discovery order.
    std::vector<int> tin(numCells, -1), low(numCells), tout(numCells), parent(numCells, -1), next(numCells, 0);
    std::vector<int> order;
    order.reserve(numCells);
    std::vector<std::pair<int, int>> closures;  // (u, v): a block closed at u by its child v
    for (int root = 0; root < numCells; ++root) {
        if (tin[root] >= 0) {
            continue;
        }
        const size_t firstClosure = closures.size();
        std::vector<int> stack{root};
        tin[root] = low[root] = int(order.size());
        order.push_back(root);
        while (not stack.empty()) {
            const int u = stack.back();
            if (next[u] < int(adj[u].size())) {
                const int v = adj[u][next[u]++];
                if (tin[v] < 0) {
                    parent[v] = u;
                    tin[v] = low[v] = int(order.size());
                    order.push_back(v);
                    stack.push_back(v);
                } else if (v != parent[u]) {
                    low[u] = std::min(low[u], tin[v]);
                }
                continue;
            }
            stack.pop_back();
            tout[u] = int(order.size()) - 1;
            const int p = parent[u];
            if (p >= 0) {
                low[p] = std::min(low[p], low[u]);
                if (low[u] >= tin[p]) {
                    closures.emplace_back(p, u);
                }
            }
        }
        // Drop one block at the root: the dilations of all its sides about
        // the root sum to a homothety.
        for (size_t k = firstClosure; k < closures.size(); ++k) {
            if (closures[k].first == root) {
                closures.erase(closures.begin() + long(k));
                break;
            }
        }
    }

    std::vector<CutCellSide> sides;
    sides.reserve(closures.size());
    for (const auto& [u, v] : closures) {
        CutCellSide side;
        side._cut = u;
        side._side.assign(order.begin() + tin[v], order.begin() + tout[v] + 1);
        std::sort(side._side.begin(), side._side.end());
        sides.push_back(std::move(side));
    }
    return sides;
}

}  // namespace PP
