#!/usr/bin/env python3
"""Tables of the H and G runs on the synthetic benchmark (the report in
benchmarks/synthetic/REPORT.md), from <results>/ours.csv (run_synthetic.py)
and <results>/oursG.csv (run_synthetic_g.py).

Usage:
  external/bourne/.venv/bin/python benchmarks/synthetic/report.py [--results benchmarks/synthetic/results]
"""

import argparse
import collections
import csv
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NOISE_RHO = {0: 0.8, 1: 0.5, 2: 0.25}


def load(path):
    return {r["case"]: r for r in csv.DictReader(open(path))}


def is_true(s):
    return s == "True"


def fmt_range(xs, f="%.1e"):
    xs = [x for x in xs if np.isfinite(x)]
    if not xs:
        return "-"
    return (f + " – " + f) % (min(xs), max(xs))


def table_groups(H, G):
    print("| noise, volumes | H converged | G converged | G invalid start | G failed | G converged / valid |")
    print("|---|---|---|---|---|---|")
    groups = collections.defaultdict(list)
    for k, r in H.items():
        groups[(int(r["noise"]), r["volumes"])].append(k)
    for (noise, vol), ks in sorted(groups.items()):
        h = sum(is_true(H[k]["converged"]) for k in ks)
        o = collections.Counter(G[k]["outcome"] for k in ks if k in G)
        valid = len(ks) - o["invalid start"]
        print("| %d (ρ %.2f), %s | %d/%d | %d/%d | %d | %d | %d/%d |" % (
            noise, NOISE_RHO[noise], vol, h, len(ks), o["converged"], len(ks), o["invalid start"], o["failed"],
            o["converged"], valid))
    hs = sum(is_true(r["converged"]) for r in H.values())
    o = collections.Counter(r["outcome"] for r in G.values())
    print("| **all** | %d/%d | %d/%d | %d | %d | %d/%d |" % (
        hs, len(H), o["converged"], len(G), o["invalid start"], o["failed"], o["converged"],
        len(G) - o["invalid start"]))


def table_invalid_by_n(G):
    print("| n | " + " | ".join("noise %d %s" % (z, v) for z in (0, 1, 2) for v in ("uniform", "voronoi")) + " |")
    print("|---|" + "---|" * 6)
    byn = collections.defaultdict(collections.Counter)
    for r in G.values():
        if r["outcome"] == "invalid start":
            byn[int(r["n"])][(int(r["noise"]), r["volumes"])] += 1
    ns = sorted({int(r["n"]) for r in G.values()})
    for n in ns:
        if byn[n]:
            print("| %d | " % n + " | ".join(str(byn[n][(z, v)]) for z in (0, 1, 2) for v in ("uniform", "voronoi"))
                  + " |")


def table_accuracy(H, G):
    print("| solver | cases | iterations | max volume residual | max shift residual | generator error / h | diagram mismatch |")
    print("|---|---|---|---|---|---|---|")
    hc = [r for r in H.values() if is_true(r["converged"])]
    gc = [r for r in G.values() if r["outcome"] == "converged"]
    for name, rows in (("H", hc), ("G", gc)):
        f = lambda key: [float(r[key]) for r in rows]
        print("| %s | %d | %s | %s | %s | %s | %s |" % (
            name, len(rows), fmt_range(f("iterations"), "%d"), fmt_range(f("volumeResidual")),
            fmt_range(f("shiftResidual")), fmt_range(f("generatorError")), fmt_range(f("diagramMismatch"))))


def list_failures(H, G):
    print("H not converged:")
    for k, r in sorted(H.items(), key=lambda kv: (kv[1]["volumes"], int(kv[1]["n"]))):
        if not is_true(r["converged"]):
            print("  %-26s %-15s %4s it %s att  res %.1e / %.1e   G: %s" % (
                k, r["status"], r["iterations"], r["attempts"], float(r["volumeResidual"]),
                float(r["shiftResidual"]), G[k]["outcome"] if k in G else "-"))
    print("G failed with a valid start:")
    for k, r in sorted(G.items(), key=lambda kv: (kv[1]["volumes"], int(kv[1]["n"]))):
        if r["outcome"] == "failed":
            print("  %-26s %-15s %4s it  generator error %.2f h   H: %s" % (
                k, r["status"], r["iterations"], float(r["generatorError"]),
                "converged" if is_true(H[k]["converged"]) else H[k]["status"]))
    e = sorted(int(r["startEmptyCells"]) for r in G.values() if r["outcome"] == "invalid start")
    if e:
        print("empty cells of the invalid starts: min %d, median %d, max %d" % (e[0], e[len(e) // 2], e[-1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    args = ap.parse_args()
    H = load(os.path.join(args.results, "ours.csv"))
    G = load(os.path.join(args.results, "oursG.csv"))
    print("## Convergence per noise level and volume distribution\n")
    table_groups(H, G)
    print("\n## Invalid starts of G per n\n")
    table_invalid_by_n(G)
    print("\n## Accuracy of the converged solves\n")
    table_accuracy(H, G)
    print("\n## Failures\n")
    list_failures(H, G)


if __name__ == "__main__":
    main()
