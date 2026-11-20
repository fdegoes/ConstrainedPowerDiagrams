#ifndef PP_CORE_DUALGRAPH_H
#define PP_CORE_DUALGRAPH_H

#include <utility>
#include <vector>

namespace PP
{

// A side of a cut cell of the dual graph (latex/appendix-nullspace.tex,
// "Cut cells & branches"; the branches B_B of latex/flow.tex): the cells S that removing the cut cell c
// separates from the rest, excluding c. Each side carries the localized
// homothety that dilates the generators of S about the generator of c.
struct CutCellSide
{
    int              _cut = -1;  // c
    std::vector<int> _side;      // S, sorted
};

// The b - 1 sides of the localized homotheties of a dual graph with b blocks
// (latex/flow.tex, eq:flow-line-constraint-arm): one side per block, found
// by the depth-first search of Hopcroft and Tarjan, except one block at the
// root of each connected component, since the dilations of all the sides at
// the root sum to a homothety. Empty when the graph has no cut cell.
//
// numCells  number of nodes (cells)
// edges     dual edges (i, j); duplicates and self-loops are ignored
std::vector<CutCellSide> FindCutCellSides(int numCells, const std::vector<std::pair<int, int>>& edges);

}  // namespace PP

#endif
