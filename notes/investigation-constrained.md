# Investigation: making the constrained integrator work

Opened 2026-11-24, after the relaxed integrator turned out to equal GM at equal
force and a PP-like non-symplectic step was ruled out (the method must be
symplectic). The paper (Sec. 5, flow results) is paused meanwhile.

## The missing element: the strength of the potential

The Lagrangian of latex/flow.tex is sum V |dq/dt|^2 - Phi(q), with Phi the
quantization error of the cells about their centroids, the particles. Its force
is -u = q - p (eq:flow-potential-gradient): each particle is pulled towards its
generator. This is the dual of GM:
- GM: generators = particles, centroids off by the gap c - q, force kappa (c - q);
- constrained: centroids = particles, generators off by u = p - q, force -u,
  with a unit strength fixed by the choice of the Lagrangian.
GM at kappa ~ 1 lets the particles cluster (gap ~ 1.4 h at T = 1); the
constrained integrator at unit strength lets the generators drift (u growing
roughly as t^2) until the diagram degenerates at t ~ 0.02 (n = 5000). The flow
shears the particle arrangement on a time h / v ~ 0.14, the unit potential
restores it on a time ~ 1.

Scaling the potential, L = sum V |dq/dt|^2 - kappa Phi (IntegratorOptions::
_stiffness for the constrained integrator; the prediction becomes
q + dt v + kappa dt^2 (q - p^k), the same discrete equation with the known part
of the force taken at the current generators):

| kappa | dt 0.001, to t = 0.1 |
|---|---|
| 0 | fails at step 19 |
| 1 | fails at step 21 |
| 100 | runs |
| 10^4 | fails at step 1 |

kappa = 100, dt = 0.001, uniform four-vortex case, n = 5000: up to t = 0.23 the
generators stay within 0.11 h of the particles (rms 0.02-0.06 h), the kinetic
energy is kept to 99.6 % (0.00250 -> 0.00249), the closest pair is 0.83 h (PP
at T = 1: 0.70 h), and volumes and centroids are exact (1e-8). At t = 0.236 all
generators jump by about h within a few steps (rms |u| 0.06 h -> 1.1 h) and
the step fails.

kappa is again a stiffness. Unlike GM's 1/eps^2 it is not tied to dt: its role
is that the regularization acts faster than the flow shears the particles,
1/sqrt(kappa) < h / v, which suggests kappa ~ (v / h)^2 (here (0.1 / 0.014)^2
~ 50; kappa = 100 works, 1 does not). To be checked.

## The second missing element: a stiff soft mode

Both the late failure and the step-size limit are one instability. Before the
failure at t ~ 0.236 (kappa = 100, dt = 0.001), the generator offsets grow as a
global anisotropic strain of the generators, u ~ S q (S11 - S22 dominant; the
dilation, rotation and translation stay negligible), which alternates in sign
every step and grows by a factor 1.09 (step 227) to 1.74 (step 236): the
signature of Stormer-Verlet beyond its stability bound.

This strain is a soft mode of the centroid map: at step 226 the generator
strain swings by 3.8e-3 per step while the particle strain changes by about
dt x 9e-5 = 9e-8, an amplification A ~ 4e4. Through the force kappa (q - p) of
the potential, the particles feel a stiffness omega^2 ~ kappa A, and Verlet is
stable only for dt omega < 2:

| kappa | dt | dt sqrt(kappa A) | observed |
|---|---|---|---|
| 100 | 0.001 | ~ 2 | stable until A grows, then the alternating mode |
| 100 | 0.002 | ~ 4 | fails at step 4 |
| 1 | 0.01 | ~ 2 | fails at steps 2-3 |
| 1 | 0.001 | ~ 0.2 | stable, but u drifts (weak potential) until t ~ 0.02 |

So the constrained system is stiff: a weak potential lets the generators
drift, a strong one makes the soft generator modes stiff for the particles.
Candidate fixes, both symplectic:
1. treat the potential implicitly (implicit midpoint, or a symplectic
   splitting with the stiff part solved exactly), which removes the dt bound
   for linear stiff modes; each step is already a nonlinear solve;
2. remove the soft modes from the force, as the homotheties are removed by
   the gauge, if they can be characterized (the anisotropic strain here).

## Open problems

1. The scaling of A with n, the flow and the domain, and whether the soft
   modes are only the affine strains.
2. The required kappa (regularization faster than the shear time h / v).
3. Which fix keeps the scheme variational and cheap.

Correction: an earlier facet script dropped the facets dual to convex-hull
edges of the generators (their cells are unbounded before clipping), and
reported corner cells with a single interior facet and cut cells next to them.
They do not exist: the integrator's dual graph has no cut cells in these runs
(sides 0 at every step), and the corner cells have two facets of about h.

## Ruled out earlier (with kappa = 1)

Start and step tolerances, a normal velocity component (SHAKE and RATTLE fail
alike), facet flips alone, a dilation gauge drift (the generators of the
integrator are those of a static H solve on the same particles, up to a
homothety), and the cut-cell detection (Hopcroft-Tarjan, b - 1 sides).
Stale velocities in the first BNOT inputs made the breakdown 8 times earlier
(fixed: bnot_centroidal.py --velocity 4vortices).

## Overnight report (2026-11-25)

See experiments/pp2d/report/REPORT.md: the explicit scheme with the softened potential does not work; Newmark (theta = 1/4) runs at dt = 0.01 until t ~ 1 (four-vortex, corner flips) or 1.5-3 (one-vortex, Newton stalls).
