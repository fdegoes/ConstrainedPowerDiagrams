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

## Open problems

1. The jump at t ~ 0.236. Not the depth-first root of FindCutCellSides (cell 0
   is in the main block, every side is a single corner cell). A cut cell at the
   top-left corner disappears around step 229; the jump comes 6 steps later.
   Next: log the solver and the gauge (s, t, s_B) over steps 225-237.
2. The step size. With kappa = 100, dt = 0.002 fails at step 4 and 0.005 at
   step 1; with kappa = 1, dt = 0.01 fails at step 3 at n = 5000. The SHAKE
   correction moves along the old generators p^k, which turn fast (the
   generators move ~10-200 times faster than the particles), so the line misses
   the constraint set. The failures depend on kappa dt^2 and on dt.
3. The scaling of the required kappa with n and v.

## Ruled out earlier (with kappa = 1)

Start and step tolerances, a normal velocity component (SHAKE and RATTLE fail
alike), facet flips alone, a dilation gauge drift (the generators of the
integrator are those of a static H solve on the same particles, up to a
homothety), and the cut-cell detection (Hopcroft-Tarjan, b - 1 sides).
Stale velocities in the first BNOT inputs made the breakdown 8 times earlier
(fixed: bnot_centroidal.py --velocity 4vortices).
