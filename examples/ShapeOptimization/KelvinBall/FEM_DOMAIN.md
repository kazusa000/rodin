# L15 graded FEM campaign

The source branch is `codex/p03-fem-domain`, based on the compiler-only Carlos
reference `3ae2c041`. The selected initial sphere has radius one and the full
outer cube has edge 30 (`outer-radius=15`). The geometry and FEM optimization
method are inherited; no BEM, new thickness correction or seed selection is used.

The reference spacing remains exactly `h0=1/6`, matching n13/L2. The same
minimum MMG size `h0/10`, far size `10h0`, Hausdorff tolerance `0.1h0^2`,
regularization `4h0` and step `0.1h0` remain fixed. A new opt-in coordinate
grid is uniform through coordinate radius 2 and then grows by 1.3, capped
at the far size. The identical monotone map on all axes preserves chamber
cut planes. It avoids allocating an n91 uniform cube before adaptation.
This is a seed grading policy, not a guarantee that every adapted edge has
exactly the requested length. Existing distance-weighted MMG adaptation
targets the same near size around the evolving interface and coarse far field.

Use the preset `reproduction/fem-domain-L15.json`. The native diagnostic checks
the original near coordinates, unchanged MMG size targets and monotone finite
graded intervals. The CPU-only runner verifies the unchanged n13/L2 one-update
FEM values against the frozen baseline, and validates L15 geometry before
starting the optimization. All stages execute on 28h in separate output folders.

The actual optimizer requests 21 evaluated states: initial state plus 20 real
level-set/MMG updates. Completed-state CSV, transport increments, reconstruction
timings, commands, source/submodule/binary hashes and actual exit codes are saved.
MMG retries are explicitly disabled to avoid silently halving the mesh scale.
An exception stops the campaign and retains inputs; no tolerance is relaxed.
The existing failed-input exports remain enabled. Full fields stay on 28h.

The launcher `scripts/run-fem-domain.py` requires an exact clean pushed commit,
new build/output directories and the original 28h baseline build and campaign
records. It inherits only the verified dependency paths/configuration, not old
binaries. CPU execution uses four OpenMP threads and one BLAS thread by default.

The historical regression uses two OpenMP threads to match its frozen campaign;
production uses four. On 28h the frozen and new binaries agree exactly within
each thread setting, but the first reconstructed shape differs across settings.
This is retained as an execution-setting dependence, not a numerical-source change.
The first validation-only failure is retained. An explicit `--reuse-native-campaign`
may reuse its successful native build only after clean source, unchanged numerical
files, dependency IDs, build path and executable hash are checked. The controller
and actually compiled numerical commits are recorded separately.

## Fixed chiral initialization

The opt-in `--seed-amplitude=0.02` uses the same fixed reflection-odd O24
radial seed as the earlier coupled research branch. For unit direction
`(x,y,z)`, `H3=3 sqrt(3) xyz`,
`H6=6 sqrt(3) (x²-y²)(y²-z²)(z²-x²)`, and the initial level set is
`|X| - (1 + amplitude H3 H6)`. Zero remains the original sphere. This is
a prescribed perturbation, not an initialization selected by measured scores.
MMG and SWIFT initialization use the same level-set function; initial geometry
errors are measured against the chosen seed, rather than a unit sphere.

Use `--preset reproduction/fem-domain-L15-chiral.json` with the runner for
the chiral repeat. Only initialization differs from the original L15 preset;
box size, near-field sizes, gradation, step, regularization and 20 real updates
stay fixed. The native seed diagnostic checks all 24 rotations, reflection
oddness, origin finiteness and the unchanged zero-amplitude sphere. Original
one-update regression is still performed with amplitude zero and two threads.
The existing fixed-step optimizer does not reject objective decreases; the
completion certificate counts transported/reconstructed updates, not accepted
Armijo steps or scientific convergence. Periodic trace and mesh warnings remain
part of the numerical assessment, and are not fixed by adding a seed.

## Historical parameter control

`reproduction/fem-domain-L15-chiral-success-params.json` restores the successful
historical `regularization=100` and `step=0.5`, while retaining L15, fixed
near-field spacing, the 0.02 O24 seed and the current MMG/transport workflow.
The physical regularization length is `100/6` and time step `0.5/6`.
This isolates the two optimizer parameters; it does not reproduce the old
L2 domain, n20 grid, thickness constraint, periodic distance projection or
L2 transport. Original one-update regression still uses baseline defaults.
Numerical sources are unchanged, so the runner can reuse the verified chiral
native build after checking source/dependency identity and its binary hash.
The allowed reuse differences are restricted to the runner, this document and
the specifically named parameter-control preset. The compiled and controller
commits are both recorded; a previous user stop does not invalidate its passed
native build and preflight checks.
