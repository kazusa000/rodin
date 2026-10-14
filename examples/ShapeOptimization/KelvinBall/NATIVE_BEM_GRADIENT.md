# Native Rodin BEM gradient: first verification stage

On `codex/p03-bem-gradient`, numerical states, derivatives and verification updates
are C++/Rodin. The first stage verifies six states, a fixed-topology coordinate
gradient and a real update. Stability sweeps and cost comparisons are deferred.

## Build and run

Use the existing Rodin CMake/dependency configuration, with examples enabled.
Build outside the source checkout and execute in a unique results directory:

```sh
cmake --build "$BUILD" --target KelvinBallBEMValidation KelvinBallBEMMMGValidation KelvinBallQuotientProjectionTest KelvinBall
cd "$RESULTS"
OMP_NUM_THREADS=1 "$BUILD/examples/ShapeOptimization/KelvinBall/KelvinBallBEMValidation"
```

`KelvinBallBEMValidation` checks the kernel formula, a closed 512-panel unit sphere,
a fixed O24 chiral surface, centered differences with three decreasing steps,
and the dilation derivative. It then backtracks an actual coordinate update,
rescales to the original enclosed volume, and solves all six states again.
Exit 0 requires both finite-difference agreement and an improved objective.
It writes `bem-validation.json`, `bem-initial.mesh` and `bem-updated.mesh`.
This is a fixed-topology surface update; it is distinct from a level-set/MMG update.

`KelvinBallQuotientProjectionTest` checks the coordinate pullback against the
existing averaged vector transport, projection idempotence, and rotated scalar
and vector traces on a canonical two-tetrahedron chamber. Averaging coincident representatives
does not repair a nonmatching surface triangulation.

`KelvinBallBEMMMGValidation` is the full-volume verification route. It uses a
13-point grid on the cube `[-2,2]^3`, a fixed `0.02 H3 H6` radial seed, the same
BEM class, Rodin H1/volume identification, Eikonal redistancing, RK4 nodal
transport and MMG level-set reconstruction. There are no artificial quotient
cuts in this fixture. It saves the two actual volume meshes and
`bem-mmg-validation.json`; both designs must pass closed-surface state checks.
It checks a positive first-order direction and reports the fresh objective
change after remeshing, including whether it improved. It does not enforce O24
symmetry during meshing/update and is not the production quotient optimizer.

The existing chamber optimizer has an opt-in native entry:

```sh
KelvinBall --bem-gradient --bem-order=3 --bem-epsilon-factor=0.6 \
  --n=5 --hmin-factor=0.8 --hmax-factor=1.2 --seed-amplitude=0.02 \
  --iterations=2 --step=0.01 --reconstruction=mmg --save-mesh
```

Two evaluated designs request one transport/reconstruction update. Every design
must have a closed sewn body before it is solved. Native mode uses the BEM load
for the Rodin H1 metric and volume projection, averages coincident quotient
representatives, and protects cut vertices/triangles in the MMG route.
It saves `bem-state-N.json`, `bem-chamber-N.mesh` and `candidate.mesh`.
It does not emit placeholder FEM velocity or pressure fields. An actual new
design's objective must be read from its fresh BEM solve, independently of the
reported first-order prediction. Remeshing changes topology and epsilon.

The successful/failed chamber cases, exact source SHA and commands belong in the
workspace's dated run record. Closed-surface verification alone does not certify
the MMG/SWIFT chamber routes or a long optimization loop. The default FEM path
and the frozen October 11 baseline retain their existing behavior.

## Discrete model

`BEM.cpp` assembles `A = Integral(Potential(kernel, trial), test)` in a vector P1
space on an embedded triangular surface. `DenseProblem` and Rodin's
`PartialPivLU` solve `A F = B`; six rigid-velocity right-hand sides share the
factorization. `R = B.transpose() * F` gives the resistance matrix.

The kernel is the three-dimensional regularized Stokeslet, including
`1/(8*pi*mu)`. Epsilon is `epsilonFactor * sqrt(totalArea/panels)`.
`RegularizedStokeslet` opts the P1 global integrator into tensor surface
quadrature using Rodin's permutation-symmetric positive Witherden--Vincent rule.
The marker leaves the previous path for other kernels unchanged. This is a
regularized Galerkin model; singular-kernel convergence is not claimed.

`k` and `q` are the mean traces of the translation and rotation diagonal blocks;
`c` averages the traces of the two coupling blocks. The objective is
`rho = abs(c)/sqrt(k*q)`. With the symmetric matrix `W = d(rho)/dR`,

```text
d(rho) = 2 trace(W F^T dB) - trace(W F^T dA F).
```

The derivative includes mapped quadrature positions, panel measures, rotational
loads, and the global area dependence of epsilon. Reference P1 coordinates and
topology are fixed during this differentiation. The implementation uses space
DOF maps, not an assumed global component numbering. It is an exact derivative
of this discrete model, not a proof of a continuous Hadamard density.
At unresolved zero coupling, `abs(c)` has no unique derivative and the gradient
request fails; sphere state verification therefore requests no derivative.

State checks reject nonfinite solutions, relative residual above `1e-9`,
operator asymmetry above `1e-12`, resistance reciprocity error above `1e-9`,
nonpositive `k/q`, open edges and nonmanifold edges. Sphere verification reports
the actual errors against `6*pi` and `8*pi` (initial coarse accuracy threshold
20%); it does not establish mesh independence.

The current support claim is local CPU, FP64, affine triangular surfaces and
dense storage. In native chamber mode, stage 3 includes state and derivative
evaluation. Its timer is therefore not directly comparable to a FEM state-only
timer. Controlled stability and timing studies require their own protocol.
