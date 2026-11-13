# Null space of the Hessian of the joint dual H (draft)

Working note, not part of the paper. Notation as in `latex/alias.tex` and
`latex/appendix-joint.tex`. Generators are assumed to be in **general
position**: no d + 2 lifted generators lie on a common non-vertical
hyperplane, so every interior vertex of the power diagram is shared by exactly
d + 1 cells, and every interior ridge (codimension-2 face) by exactly 3 cells.

## Statement

Let V_1, ..., V_n be the cells of the power diagram of the shifted sites
(p_i, lambda_i), clipped to the domain Omega, all non-empty. Call a facet
A_ij = V_i ∩ V_j (of positive (d-1)-measure) *interior*, and an interior ridge
a (d-2)-face of the diagram that lies in the interior of Omega. Define the
**facet graph** G: one node per interior facet, and an edge between two facets
that share an interior ridge. Let C be its number of connected components.

**Proposition.** Under general position,

    ker ∇²H  ≅  { piecewise-affine psi on the diagram, continuous across every interior facet },

and

    dim ker ∇²H = (d + 1) + C.

For a connected facet graph (C = 1) this is the (d + 2)-dimensional family of
homotheties of eq:joint-nullspace. Every additional component adds exactly
one mode, a local "hinge" or "fold" of that component.

## Step 1: the kernel is the space of continuous piecewise-affine fields

From appendix-joint, with delta psi_i(x) = 2 delta u_i^T (x - cbar_i) + delta w_i,

    delta z^T ∇²H delta z = - sum_{interior A_ij} 1 / (2 |p_i - p_j|) ∫_{A_ij} (delta psi_i - delta psi_j)^2 dA.

Boundary facets do not appear. Every term is non-positive, so delta z is in the
kernel iff delta psi_i = delta psi_j on every interior facet A_ij (both are
affine and the facet has positive measure, so equality a.e. on A_ij is
equality on its affine hull).

The map delta z = (delta u_i, delta w_i)_i -> (delta psi_i)_i is a linear
bijection between R^{(d+1) n} and the space of collections of affine functions
(one per cell): delta u_i is half the gradient of delta psi_i, and
delta w_i = delta psi_i(cbar_i). Hence

    ker ∇²H  ≅  L := { (psi_i)_i affine, psi_i = psi_j on A_ij for every interior facet }.

## Step 2: splitting off the affine part

Global affine functions psi_i = a^T x + b for all i are in L. This gives d + 1
dimensions; they are the translations t and the constant c of
eq:joint-nullspace. Write L = Aff ⊕ L_0, with L_0 a complement (e.g. fields
that vanish identically on one chosen cell).

## Step 3: gradient jumps (Maxwell-Cremona)

For psi in L and an interior facet A_ij with unit normal n_ij (from V_i to
V_j), psi_j - psi_i is affine and vanishes on the affine hull of A_ij, hence

    ∇psi_j - ∇psi_i = kappa_ij n_ij      for a scalar kappa_ij,

and psi_j - psi_i = kappa_ij n_ij^T (x - x_ij) with x_ij on A_ij. Conversely,
on a connected set of cells, the gradients and the field are determined from
one cell by integrating the jumps across facets.

**Closure at a ridge.** Around an interior ridge R, shared (general position)
by exactly three cells V_a, V_b, V_c and three facets A_ab, A_bc, A_ca, going
around R returns to the starting cell, so

    kappa_ab n_ab + kappa_bc n_bc + kappa_ca n_ca = 0.            (*)

The three normals lie in the 2-plane orthogonal to R and are pairwise
non-parallel (the three facets are distinct and meet transversally), so they
span that plane and satisfy exactly one linear relation. Therefore (*) fixes
the three kappa's up to one common scalar: **the kappa's of the facets around
a ridge are proportional with fixed ratios.** (The affine parts are consistent
automatically: all three differences vanish on R.)

This is the Maxwell-Cremona correspondence: the kappa's are a self-stress of
the reciprocal framework (the edges p_i p_j of the regular triangulation),
since n_ij is parallel to p_j - p_i.

## Step 4: one scalar per component of the facet graph

Propagating (*) along the edges of G, all kappa's on one connected component
K of G are fixed up to a single scalar s_K. Conversely, every choice of the
scalars (s_K)_K is realized by an element of L_0:

* **Existence of a nonzero solution on each component.** The power diagram
  itself provides one: its lifting psi_i^0(x) = 2 p_i^T x - |p_i|^2 + lambda_i
  is continuous across all facets (it is the max of the lifted affine
  functions, which agree on the facets), with jumps kappa^0_ij = 2 |p_j - p_i|
  != 0 on every facet. So the ratios imposed by (*) are those of kappa^0, and
  s_K kappa^0 restricted to K is a valid jump pattern for every s_K.
  (On the whole diagram, psi^0 is the dilation mode mu.)

* **Independence across components, and integrability.** Different
  components share no ridge, so their constraints are independent. To build
  psi from the jumps s_K kappa^0, integrate from a reference cell along paths
  of cells crossing facets. The result is path independent iff it closes
  around every interior ridge (the only places where cycles of cells close up
  in the interior); this is exactly (*). Cycles around boundary pieces impose
  nothing, since the boundary carries no facet constraint. (In a non-simply
  connected domain, cycles of cells around a hole also close trivially here,
  because the jumps along such a cycle sum to the jumps around the ridges it
  encloses... TODO: check; all numerical tests are simply connected or the
  L-shape.)

Therefore L_0 ≅ R^C (one scalar per component), and

    dim ker ∇²H = dim Aff + dim L_0 = (d + 1) + C.        ∎ (modulo TODO)

## Consequences

* **C = 1 (connected facet graph).** The kernel is exactly the d + 2
  homotheties: d + 1 affine modes plus the global scale of the lifting
  (the dilation mu). This is the generic case and the statement of
  eq:joint-nullspace.

* **Single-facet cell.** A cell V_i with exactly one interior facet A_ij:
  the facet has no interior ridge (a third cell at a ridge of A_ij would be a
  second interior neighbour of V_i), so it is an isolated node of G and adds
  one mode. Explicitly, delta psi_i = n_ij^T (x - x_ij) on V_i and 0 elsewhere:
  the lifted plane of V_i rotates about the facet. This is *sufficient but not
  necessary* for extra modes.

* **Cells with several facets can also create modes.** Any component of
  interior facets that touches the rest of the skeleton only through the
  domain boundary adds a mode, e.g. a chain of facets cutting across the
  domain from boundary to boundary. In the 2D tests, 8 of 59 diagrams with
  extra modes had no single-facet cell (components of 2-5 facets over 3-4
  cells). So **no bound on the number of facets of a cell characterizes the
  extra modes**; the condition is the connectivity of the facet graph.

* **The diagram is still unique.** All kernel directions leave every cell
  unchanged (they are liftings of the same partition), so the non-uniqueness
  is only in the multipliers, as for the homotheties.

## Numerical check

`dim ker = d + 1 + C` matched the numerical nullity of the assembled Hessian
in 480/480 random 2D diagrams (unit square and L-shape, n = 4 ... 80; C from
1 to 4). Not yet checked in 3D (needs ridge adjacency from the cell
polyhedra); in all 3D benchmarks so far the nullity was d + 2.

## Open points

1. The non-simply-connected case (TODO in Step 4).
2. 3D numerical check.
3. Citation: Maxwell-Cremona in the plane (Hopcroft-Kahn 1992; Whiteley
   1982); higher dimensions (Rybnikov; Aurenhammer, "A criterion for the
   affine equivalence of cell complexes in R^d and convex polyhedra in
   R^{d+1}", 1987, which also covers power diagrams as projections of
   polyhedra).
