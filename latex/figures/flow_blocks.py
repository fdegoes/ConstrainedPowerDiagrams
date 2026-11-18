#!/usr/bin/env python3
"""Draws the blocks figure of the constrained dynamics (flow.tex, Sec. 5.1): a power diagram in
general position on an N-shaped domain, whose two bars are joined by a diagonal. The cells at the
two inner corners of the N are cut cells, which split the dual graph into a chain of three blocks:
the lower part of the left bar (root), the diagonal with the upper part of the left bar and the
lower part of the right bar (B_1), and the upper part of the right bar (B_2).
  figures/flow-blocks.tex   the cells tinted by block (root in gray), the cut cells in orange
  figures/flow-dfs.tex      the depth-first search of Algorithm 3 from the bottom of the left bar:
                            every cell points to its parent along a dual edge, coloured by block,
                            the search enters every block through the red cells, and the facets
                            that cut off each branch are drawn thick
Needs numpy (and matplotlib, for the geometry helpers of nullspace_reducible.py):
  external/bourne/.venv/bin/python latex/figures/flow_blocks.py [--seed 1]
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nullspace_reducible as nr  # noqa: E402

import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
S, H = 0.36, 0.5                                      # width of the bars, height of the diagonal
BOX = (0.0, 0.0, 1.0, 1.0)
PIECES = [np.array(P, float) for P in (                # convex pieces of the N, counter-clockwise
    [(0, 0), (S, 0), (S, 1), (0, 1)],
    [(S, 1 - H), (1 - S, 0), (1 - S, H), (S, 1)],
    [(1 - S, 0), (1, 0), (1, 1), (1 - S, 1)])]
OUTLINE = [(0, 0), (S, 0), (S, 1 - H), (1 - S, 0), (1, 0), (1, 1), (1 - S, 1), (1 - S, H), (S, 1), (0, 1)]
CORNERS = [np.array((S, 1 - H)), np.array((1 - S, H))]   # the inner corners of the N
FILL = ["nsGray!22", "nsBlue!28", "nsGreen!32"]            # root block, B_1, B_2
EDGE = ["black!55", "nsBlue", "nsGreen!75!black"]
LTOL = 1e-7


def poly_clip(P, Q):
    """The part of the convex polygon P inside the convex counter-clockwise polygon Q."""
    for A, B in zip(Q, np.roll(Q, -1, axis=0)):
        if len(P) < 3:
            return []
        n = np.array([B[1] - A[1], A[0] - B[0]])         # outward normal of the edge AB
        P = nr.clip(P, n, n @ A)
    return P if len(P) >= 3 else []


def seg_clip(A, B, Q):
    """The part of the segment AB inside the convex counter-clockwise polygon Q (Cyrus-Beck)."""
    d = B - A
    t0, t1 = 0.0, 1.0
    for E, F in zip(Q, np.roll(Q, -1, axis=0)):
        n = np.array([F[1] - E[1], E[0] - F[0]])
        num, den = n @ (E - A), n @ d
        if abs(den) < 1e-15:
            if num < 0:
                return None
            continue
        t = num / den
        if den > 0:
            t1 = min(t1, t)
        else:
            t0 = max(t0, t)
    return (A + t0 * d, A + t1 * d) if t1 - t0 > LTOL else None


def diagram(Y, w):
    """Cells, as pieces clipped to the convex pieces of the N, and facets, as clipped segments."""
    conv = nr.voronoi(Y, w, BOX)
    pieces = [[np.array(q) for q in (poly_clip([tuple(p) for p in P], Q) for Q in PIECES) if q] for P in conv]
    facets = {}
    for i, P in enumerate(conv):
        for j in range(i + 1, len(Y)):
            a = 2.0 * (Y[j] - Y[i])
            b = Y[j] @ Y[j] - Y[i] @ Y[i] - w[j] + w[i]
            on = [k for k in range(len(P)) if abs(a @ P[k] - b) < 1e-9 * (1 + np.linalg.norm(a))]
            if len(on) < 2:
                continue
            d = np.array([-a[1], a[0]])                # extreme points along the line, since the
            s = [d @ P[k] for k in on]                  # clipping may repeat a vertex
            A, B = P[on[int(np.argmin(s))]], P[on[int(np.argmax(s))]]
            segs = [g for g in (seg_clip(A, B, Q) for Q in PIECES) if g is not None]
            if sum(np.linalg.norm(g[1] - g[0]) for g in segs) > 1e-6:
                facets[(i, j)] = segs
    return pieces, facets


def centroid(ps):
    a = np.array([nr.area(P) for P in ps])
    return sum(ai * nr.centroid(P) for ai, P in zip(a, ps)) / a.sum()


def touching(P, Q):
    """Whether two convex polygons share a boundary segment of positive length."""
    for A, B in zip(P, np.roll(P, -1, axis=0)):
        d = (B - A) / np.linalg.norm(B - A)
        n = np.array([-d[1], d[0]])
        for E, F in zip(Q, np.roll(Q, -1, axis=0)):
            if abs(n @ (E - A)) < LTOL and abs(n @ (F - A)) < LTOL:
                lo, hi = sorted((d @ (E - A), d @ (F - A)))
                if min(hi, d @ (B - A)) - max(lo, 0.0) > 1e-6:
                    return True
    return False


def connected(ps):
    seen, todo = {0}, [0]
    while todo:
        u = todo.pop()
        for v in range(len(ps)):
            if v not in seen and touching(ps[u], ps[v]):
                seen.add(v)
                todo.append(v)
    return len(seen) == len(ps)


def on_lines(p, polys):
    for Q in polys:
        Q = np.vstack([Q, Q[:1]])
        for A, B in zip(Q[:-1], Q[1:]):
            d = B - A
            t = np.clip((p - A) @ d / (d @ d), 0, 1)
            if np.linalg.norm(A + t * d - p) < 1e-7:
                return True
    return False


def blocks_dfs(n, edges, r):
    """Algorithm 3: the DFS from r, with the parents and the blocks, each closed at a cell u."""
    nbr = [[] for _ in range(n)]
    for i, j in edges:
        nbr[i].append(j)
        nbr[j].append(i)
    t, low, parent = [-1] * n, [0] * n, [-1] * n
    stack, blocks, clock = [], [], [0]

    def visit(u):
        t[u] = low[u] = clock[0]
        clock[0] += 1
        for v in nbr[u]:
            if t[v] < 0:
                parent[v] = u
                stack.append((u, v))
                visit(v)
                low[u] = min(low[u], low[v])
                if low[v] >= t[u]:               # u cuts off the subtree of v
                    es = set()
                    while True:
                        e = stack.pop()
                        es.add(frozenset(e))
                        if e == (u, v):
                            break
                    blocks.append((u, v, es))
            elif v != parent[u] and t[v] < t[u]:
                stack.append((u, v))
                low[u] = min(low[u], t[v])

    sys.setrecursionlimit(10000)
    visit(r)
    return parent, blocks


def subtree(parent, v):
    kids = {}
    for c, p in enumerate(parent):
        kids.setdefault(p, []).append(c)
    out, todo = set(), [v]
    while todo:
        u = todo.pop()
        out.add(u)
        todo += kids.get(u, [])
    return out


def in_domain(p):
    return any(nr.inside(Q, p) for Q in PIECES)


def example(seed, n, cut_offset, w_cut, gap, iters, jitter):
    rng = np.random.default_rng(seed)
    # the cut sites sit inside the bars, level with the inner corners
    C = np.array([CORNERS[0] - (cut_offset, 0.0), CORNERS[1] + (cut_offset, 0.0)])
    Y = []
    while len(Y) < n:
        p = rng.random(2)
        if in_domain(p) and np.min(np.linalg.norm(C - p, axis=1)) > gap:
            Y.append(p)
    Y = np.vstack([np.array(Y), C])
    m = len(Y)
    cuts = (m - 2, m - 1)
    w = np.zeros(m)
    w[list(cuts)] = w_cut
    for _ in range(iters):                            # Lloyd, with the cut sites fixed
        pieces, _ = diagram(Y, w)
        Y[:n] = [centroid(ps) for ps in pieces[:n]]
    Y[:n] += jitter * rng.standard_normal((n, 2))
    pieces, facets = diagram(Y, w)

    assert all(pieces), "empty cell"
    assert all(connected(ps) for ps in pieces), "a cell reaches across a notch"
    assert all(any(nr.inside(P, y) for P in ps) for ps, y in zip(pieces, Y)), "a generator lies outside its cell"
    ends = {}
    for segs in facets.values():
        for A, B in segs:
            for p in (A, B):
                if not on_lines(p, PIECES):
                    key = tuple(np.round(p, 6))
                    ends[key] = ends.get(key, 0) + 1
    assert all(c == 3 for c in ends.values()), "not in general position: %s" % sorted(set(ends.values()))
    r = min(range(n), key=lambda k: Y[k][1] + Y[k][0])  # bottom of the left bar
    parent, blocks = blocks_dfs(m, list(facets), r)
    assert all(p >= 0 for k, p in enumerate(parent) if k != r), "the dual graph is not connected"
    found = {u for u, _, _ in blocks if u != r}
    assert found == set(cuts), "cut cells %s, expected %s" % (sorted(found), list(cuts))
    assert len(blocks) == 3, "%d blocks" % len(blocks)
    blocks = [next(b for b in blocks if b[0] == u) for u in (r,) + cuts]
    branch = [None] + [subtree(parent, v) - {u} for u, v, _ in blocks[1:]]
    R = [[b for b in (1, 2) if k in branch[b]] for k in range(m)]
    print("blocks: n = %d, %d dual edges, |R_i| = %s" % (m, len(facets), np.bincount([len(x) for x in R])))
    return Y, pieces, facets, r, parent, blocks, cuts, branch, R


def fmt(p):
    return "(%.4f,%.4f)" % (p[0], p[1])


def path(P):
    return " -- ".join(fmt(p) for p in P) + " -- cycle"


def picture(body):
    head = ("%% generated by figures/flow_blocks.py; do not edit by hand\n\\begin{tikzpicture}[x=4cm,y=4cm]\n"
            "  \\useasboundingbox (0,0) rectangle (1,1);\n")
    return head + body + "\\end{tikzpicture}\n"


def frame():
    return "  \\draw[line width=0.8pt, line join=miter] %s;\n" % path(OUTLINE)


def facet_lines(facets, style, skip=()):
    return "".join("  \\draw[%s] %s -- %s;\n" % (style, fmt(A), fmt(B))
                   for e, segs in facets.items() if e not in skip for A, B in segs)


def seg_dist(p, A, B):
    d = B - A
    t = np.clip((p - A) @ d / max(d @ d, 1e-15), 0, 1)
    return np.linalg.norm(A + t * d - p)


def free_spot(candidates, dots, segs):
    """The candidate farthest from the dots and the segments."""
    def clearance(p):                                 # labels are about twice as wide as tall
        return min([np.linalg.norm((p - q) * (0.5, 1.0)) for q in dots] + [seg_dist(p, A, B) for A, B in segs])
    return max(candidates, key=clearance)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n", type=int, default=54, help="sites besides the two cut sites")
    ap.add_argument("--cut-offset", type=float, default=0.18, help="distance of the cut sites from the corners")
    ap.add_argument("--w-cut", type=float, default=0.04, help="weight of the two cut sites")
    ap.add_argument("--gap", type=float, default=0.2, help="initial clearance around the cut sites")
    ap.add_argument("--iters", type=int, default=100)
    ap.add_argument("--jitter", type=float, default=0.004, help="perturbation of the sites after Lloyd")
    args = ap.parse_args()
    Y, pieces, facets, r, parent, blocks, cuts, branch, R = example(
        args.seed, args.n, args.cut_offset, args.w_cut, args.gap, args.iters, args.jitter)
    m = len(Y)
    block_of_edge = {e: b for b, (_, _, es) in enumerate(blocks) for e in es}
    home = [next(b for b, (_, _, es) in enumerate(blocks) if any(k in e for e in es)) for k in range(m)]
    i = max((k for k in range(m) if len(R[k]) == 2), key=lambda k: Y[k][1] + Y[k][0])      # top right
    j = min((k for k in range(m) if R[k] == [1] and k not in cuts),
            key=lambda k: np.linalg.norm(Y[k] - (0.5, 0.5)))                               # centre
    root_cells = [k for k in range(m) if not R[k] and k not in cuts and k != r]
    k_ = min(root_cells, key=lambda k: np.linalg.norm(Y[k] - centroid([P for q in root_cells for P in pieces[q]])))
    dot = "  \\fill[black] %s circle (0.012);\n"
    big = "  \\fill[nsRed] %s circle (0.022);\n"
    border = np.array(OUTLINE + OUTLINE[:1], float)
    arrows = [(Y[k], Y[p]) for k, p in enumerate(parent) if p >= 0]
    grid = [np.array([x, y]) for x in np.arange(0.04, 0.97, 0.01) for y in np.arange(0.04, 0.97, 0.01)]

    def cells_of(ks):
        return [g for g in grid if any(nr.inside(P, g) for k in ks for P in pieces[k])]

    # left: the cells tinted by block, the cut cells in orange
    body = "".join("  \\fill[%s] %s;\n" % ("nsCut" if k in cuts else FILL[home[k]], path(P))
                   for k, ps in enumerate(pieces) for P in ps)
    body += facet_lines(facets, "black!80, line width=0.35pt")
    body += frame()
    facet_segs = [g for segs in facets.values() for g in segs] + list(zip(border[:-1], border[1:]))
    for b, name in ((0, "$B_0$"), (1, "$B_1$"), (2, "$B_2$")):
        ks = [k for k in range(m) if home[k] == b and k not in cuts]
        c = centroid([P for k in ks for P in pieces[k]])  # the centroid of the block, moved
        spots = cells_of(ks)                               # into it if it falls outside
        body += "  \\node[nsLabel, fill=white, rounded corners=1pt] at %s {%s};\n" % (
            fmt(min(spots, key=lambda g: np.linalg.norm(g - c))), name)
    for c, name in zip(cuts, ("$c_{B_1}$", "$c_{B_2}$")):
        body += "  \\node[nsLabel] at %s {%s};\n" % (fmt(centroid(pieces[c])), name)
    with open(os.path.join(HERE, "flow-blocks.tex"), "w") as f:
        f.write(picture(body))

    # right: the DFS tree, every cell pointing to its parent, coloured by block, the facets that
    # cut off each branch drawn thick in the colour of its block
    cutoff = {b: [e for e in facets if blocks[b][0] in e and (set(e) - {blocks[b][0]}) <= branch[b]] for b in (1, 2)}
    body = "".join("  \\fill[yellow!45] %s;\n" % path(P) for k in (i, j, k_) for P in pieces[k])
    body += facet_lines(facets, "nsGray!70, line width=0.35pt", skip=cutoff[1] + cutoff[2])
    for b in (1, 2):
        body += facet_lines({e: facets[e] for e in cutoff[b]}, "%s, line width=1.4pt, line cap=round" % EDGE[b])
    body += frame()
    for k, p in enumerate(parent):
        if p >= 0:
            body += "  \\draw[->, >=stealth, %s, line width=0.6pt, shorten <=1.3pt, shorten >=2.2pt] %s -- %s;\n" % (
                EDGE[block_of_edge[frozenset((k, p))]], fmt(Y[k]), fmt(Y[p]))
    body += "".join((big if k == r or k in cuts else dot) % fmt(y) for k, y in enumerate(Y))
    thick = [g for b in (1, 2) for e in cutoff[b] for g in facets[e]]
    keep = arrows + thick + list(zip(border[:-1], border[1:]))
    halo = "nsLabel, fill=white, inner sep=0.5pt, rounded corners=1pt"
    placed = []
    for k, name in ((i, "$i$"), (j, "$j$"), (k_, "$k$")):
        ring = [Y[k] + 0.04 * np.array([np.cos(t), np.sin(t)]) for t in np.linspace(0, 2 * np.pi, 16, endpoint=False)]
        placed.append(free_spot(ring, [q for q in Y if q is not Y[k]], keep))
        body += "  \\node[%s] at %s {%s};\n" % (halo, fmt(placed[-1]), name)
    # the branch labels go in the notches of the N, next to the edge of their branch: the
    # diagonal for B_1 (bottom notch) and the upper right bar for B_2 (top notch)
    def diagonal_x(y, upper):
        y0 = 1.0 if upper else 1 - H
        return S + (y0 - y) / (1 - H) * (1 - 2 * S)
    spots = {1: (diagonal_x(0.1, False) - 0.065, 0.1), 2: (1 - S - 0.065, 0.86)}
    for b, name in ((1, "$\\branch_{B_1}$"), (2, "$\\branch_{B_2}$")):
        body += "  \\node[nsLabel, text=%s] at %s {%s};\n" % (EDGE[b], fmt(spots[b]), name)
    with open(os.path.join(HERE, "flow-dfs.tex"), "w") as f:
        f.write(picture(body))
    print("i = %d with R_i = %s, j = %d with R_j = %s, k = %d with R_k = %s" % (i, R[i], j, R[j], k_, R[k_]))


if __name__ == "__main__":
    main()
