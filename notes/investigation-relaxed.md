# Investigation A: the relaxed integrator as the variational method

Opened 2026-11-24. The paper (Sec. 5 and the flow results) is paused until
this and `investigation-variational-pp.md` give a definite answer.

## Setting

Relaxed (softened) integrator, `src/core/softenedIntegrator.cpp`: the particle
q is the site of the joint dual, the targets are q, and the power cell of
particle i has the generator p_i = q_i + u_i. At the maximizer of H_alpha,
u_i = alpha (q_i - c_i), so the centroid, the particle and the generator are
collinear, with q between them and |p - q| = alpha |q - c|:
- alpha -> 0: p = q, GM with kappa = 1;
- alpha -> inf: c = q, the constrained integrator.
The force is (1 + alpha)(c - q), with c the centroid of the shifted diagram.

## What is known (four-vortex inputs, n = 5000, corrected velocities)

- The stability bound dt < 2/sqrt(1 + alpha) is not the binding one.
- The binding limit is the distinctness of the generators q + u. At dt = 0.01,
  T = 1, uniform case: alpha = 0.1, 0.5, 1 run; 2 fails at step 97, 3 at 87,
  5 at 76, 10 at 58, 100 at 20. With the stale velocities of the first BNOT
  inputs the limit was about 30 times lower (0.1 failed).
- The shift grows as alpha times the residual; its mechanism is open. The
  estimate alpha ~ (h / 2 v_n)^2, based on a first-order drift off the
  realizable set, was withdrawn: the velocity is tangent (cos(v, q) = 1e-3).
- GM at eps = 3 dt (kappa 1111 and 4444) runs at both dt. PP, GM and the relaxed
  runs at alpha 0.028 and 0.003 run to T = 1 on all three cases at dt 0.01 and
  0.005, with step tolerance 1e-8.

## Questions

1. Is the relaxed scheme better than GM at equal force, 1 + alpha = kappa?
   Metrics per run: total energy drift, centroid gap max |c - q|, closest
   particle pair (clustering), volume error, Newton iterations and time.
2. How does the largest usable alpha scale with n (1000, 5000, 20000), the
   flow speed (0.05, 0.1, 0.2) and dt (0.01, 0.005, 0.0025)? Is it a property
   of the resolution (like GM's eps), and can it be predicted?
3. What makes the shifts grow: the residual |q - c| along the run, or the
   shift for a fixed residual? Track both.

## Plan

1. A runner for the GM-versus-relaxed comparison at equal force, on the
   uniform case first: alpha in {0.1, 0.5, 1}, kappa = 1 + alpha, both dt.
2. A summary script over `energy.csv` and the frame CSVs for the metrics of
   question 1, plus the closest pair per frame.
3. The scaling study of question 2, with BNOT inputs at n = 1000 and 20000
   (`experiments/pp2d/bnot_centroidal.py`, field `--velocity 4vortices`).

## Outcome

To be written: whether the relaxed scheme is a better variational method than
GM, and what fixes its alpha.
