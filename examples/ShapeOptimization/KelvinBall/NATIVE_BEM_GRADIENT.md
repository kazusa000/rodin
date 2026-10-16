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

The full-volume fixture also accepts `--seed-amplitude=A` (default `0.02`) and
`--step-factors=0,0.01,0.03,0.1` (default `0.01`). Actual transport time is
the factor times the grid spacing. Multiple factors share one initial shape,
state, H1 direction and signed-distance field; candidates are independent,
not successive optimization iterations. A zero factor copies the signed-distance
field without transport and isolates redistancing/reconstruction changes.
Each candidate has a separate `step-N/` mesh and JSON record. Exit zero means
the requested reconstruction/state checks ran successfully; `improved` records
the objective comparison and does not certify volume acceptance or optimization.
The seed coefficient multiplies `H3 H6`; its maximum radial deformation is
approximately `0.2963 A`, not `A`.

`scripts/bem-gradient-seed-step-sweep.sh EXECUTABLE NEW_RESULTS_DIRECTORY`
runs amplitudes `0.02`, `0.2`, `0.6`, each with the four factors above, and
records real exit codes. Existing results directories are refused. The BEM
quadrature, epsilon policy, grid, smoothing and reconstruction settings remain
fixed. All numerical calculations still use the native C++ executable.

## Consecutive accepted full-volume updates

```sh
KelvinBallBEMMMGValidation --seed-amplitude=0.6 --step-factors=0.03 \
  --updates=3 --volume-tolerance=0.005 --max-backtracks=6
```

With `--updates>0`, the fixture recomputes six states and the discrete derivative
on each accepted shape, rebuilds its H1 system, and starts each update at the
specified factor. Rejected trials halve the factor, with at most six halvings
by default. Candidates are accepted only after closed-surface state checks,
strict rho improvement and weak Armijo gain
`rho_new >= rho_old + 1e-4 * dt * slope`. Total volume drift is bounded relative
to the original obstacle volume, not just the preceding iteration. The default
0.5% volume tolerance is a configurable coarse-mesh experimental policy, not
exact volume conservation or an accuracy estimate. No post-hoc rescaling or
geometric snapping is used to obtain acceptance.

Trial fields borrow an independent working copy; rejection preserves the accepted
mesh. Every attempt retains its candidate mesh and JSON, including volume-only
rejections where the expensive BEM state is deliberately skipped. Checkpoints
`accepted-N.mesh` and `bem-mmg-history.json` record accepted progress. Exit zero
requires all requested updates; exit two means bounded backtracking exhausted
without accepting a candidate. Other failures exit one. Existing accepted histories
are refused. `--updates=0` preserves the independent seed/step diagnostic mode.

`scripts/bem-gradient-short-run.sh EXECUTABLE NEW_RESULTS_DIRECTORY` runs this
three-update configuration. The route still does not enforce O24 symmetry;
the production quotient seam issue and stability/cost studies remain separate.

Long background campaigns may request up to 1000 updates. To continue an accepted
checkpoint without refitting it, supply `--initial-mesh=ACCEPTED.mesh`,
`--target-volume=ORIGINAL_VOLUME` and `--expected-rho=ACCEPTED_RHO`. The native
MMG loader restores mesh metadata; the first fresh BEM solve must reproduce
the accepted rho to relative `1e-10`. Volume is checked against the supplied
original target, so restarting does not reset the allowed cumulative drift.
All other numerical settings and acceptance rules are unchanged. Checkpoints
in a resumed campaign use local numbering starting at zero; record the prior
step offset in the campaign record. The update limit is an operational bound,
not a claim of convergence. Line-search exhaustion still terminates the run.

`scripts/bem-gradient-long-run.sh EXECUTABLE NEW_RESULTS_DIRECTORY ACCEPTED_MESH
ORIGINAL_VOLUME ACCEPTED_RHO` requests 100 further updates with seed provenance
`0.6`, factor `0.03`, volume tolerance `0.005`, and six backtracks. The shell
only launches the native executable and records its real exit status.

To retain the H1 volume-nullspace projection while disabling volume-based
rejection, append `no-volume-rejection` to that script invocation. This selects
`--no-volume-rejection` and factor `0.10` (physical transport time `0.10*h`,
versus the previous `0.03*h`). The native flag disables only the starting-mesh
and candidate cumulative volume tolerance guards. The projected H1 direction
and its two solves are unchanged; no volume penalty, correction or rescaling
is added. The original target remains a reporting reference, with
`volume_projection_enabled=true`, `volume_rejection_enabled=false` and
`volume_tolerance=null` in the history. Positive finite obstacle volume and
closed-surface/state validity remain necessary to evaluate the objective.
Strict rho increase, weak Armijo and six halvings remain active. The previous
volume-rejection mode remains available for replaying earlier experiments.

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
