# Viscosity: to revisit

None of our integrators (constrained, softened, GM, Power Particles) has a
viscosity step yet (benchmarks/pp2d skips the `viscosity:` parameter).

The 2015 Power Particles code (Before2015/pp, pp2/SimDEC.cpp `diffuse`,
pp2/SimDyn.cpp `apply_viscosisty`) applies it after the pressure step as an
implicit diffusion per velocity component,

    (diag(V) + nu dt L_star) v* = diag(V) v,

with L_star the one-ring star Laplacian of the cells (star_ij = A_ij / l_ij,
domain-boundary and air facets skipped, no density scaling). It is used by
`rt` (0.002) and the free-surface examples (circle, layers, splash, wave).

Open question: how to add viscosity to our variational formulation
(latex/flow.tex), e.g. as a Rayleigh dissipation term in the discrete
Lagrangian, a split implicit step as in Power Particles, or through the
facet Laplacian of H. To decide before running rt or the free-surface cases.
