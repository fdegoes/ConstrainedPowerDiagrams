# Investigation C: a variational Power Particles

Opened 2026-11-24. The paper (Sec. 5 and the flow results) is paused until
this and `investigation-relaxed.md` give a definite answer.

## Idea

GM (Gallouet & Merigot 2018, eq. 5.4) with the stiffness kappa = 1/eps^2 is
  v^{k+1} = v^k + kappa dt (c^k - q^k),   q^{k+1} = q^k + dt v^{k+1},
with c^k the centroids of the cells of q^k with the prescribed volumes. With
eps = dt this gives
  q^{k+1} = c^k + dt v^k,
which is the site update of Power Particles (a Lloyd step onto the centroids,
then advection; latex/flow.tex, Link to deGoes2015PowerParticles), except that
PP advects with its projected (discretely divergence-free) velocities and stores
them as the new velocities, while GM's velocities absorb the Lloyd correction
(c - q) / dt.

So Power Particles is close to GM at eps = dt plus a velocity projection. If
the projection can be written as a variational step (a constraint on the
velocities, as RATTLE projects onto the tangent space of a holonomic
constraint, or a discrete Lagrange-d'Alembert step for the nonholonomic
condition that the velocities preserve the volumes to first order), the
result is a variational scheme whose stiffness is tied to dt by construction
rather than tuned.

## Caveats

- eps = dt is outside GM's convergence regime (their Theorem 1.4 needs
  tau = o(eps^2); their experiments use tau / eps between 0.2 and 0.5). The
  scheme is stable (dt < 2 eps), but its energy behaviour is untested.
- The Lloyd correction (c - q)/dt in GM's velocity is a large, non-smooth
  term; PP discards it from the stored velocity. Whether a variational
  projection keeps PP's good behaviour is the question.

## Plan

1. Check numerically that GM with kappa = 1/dt^2 and PP have the same site
   update on the four-vortex inputs (one step, same state), and compare their
   runs to T = 1: energy, centroid gap, volume error, closest pair.
2. Write PP's velocity projection as the minimizer of a discrete action or as
   a projection onto the tangent space of the volume constraint, and derive
   the step that replaces it (notes, then a prototype integrator).
3. Prototype the integrator in `src/core` next to `gmIntegrator.cpp`, with a
   test that it reduces to PP's site update.

## Outcome

To be written.
