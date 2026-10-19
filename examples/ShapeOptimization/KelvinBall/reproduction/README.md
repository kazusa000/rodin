# FEM old: Cubic using the archived Terra numerical pipeline

Base: `e004696706c169a775f039dc0782b07205c370da`, the Rodin source of the
2026-10-01 Terra result preserved in `out/project-03/success/terra-rho-unknown`.
The user explicitly selected a fixed amplitude-0.02 Cubic chiral seed:
`r = 1 + 0.02*H3*H6`. There is no seed search. The existing Terra seed is unchanged.

The old solver, level-set transport, MMG reconstruction and shape-update code
remain unchanged. The only change inside the optimization implementation is the
printed seed description. `SymmetryTest.cpp` checks every proper rotation and
also checks that the Cubic seed breaks reflection symmetry and uses the specified
amplitude. Both Cubic and Terra geometry tests must pass before optimization.

Run only on 28h. `run_fem_old.py` uses the archived parameters:
`--n=20 --outer-radius=2 --regularization=100 --step=0.5
--advection-quadrature=3 --reconstruction=mmg --mmg-adapt --mmg-snap=0.005
--thickness-min=0 --iterations=20 --save-mesh`, OMP6/BLAS1.
Twenty solved states are nineteen real updates. There is no BEM boundary feedback,
new transport, volume rescaling or thickness constraint.

The standalone CMake wrapper changes build integration only and records the
actual configuration, compiler exits, binary SHA, native exit and complete CSV
rows. It uses the existing 28h dependency prefix. Outputs/builds must be fresh and
outside the source tree. The runner refuses non-28h hosts and dirty/unpinned source.

After optimization, export the original saved HDF5 doubles with the already
validated frozen exporter `2b3466878cd8c063a5c09ef5c5cd17a4bcf1ba75`; retain the
last saved state's relation to the last solved state. Wait for this task's Terra
BEM process to exit before starting one independent Cubic BEM measurement. The
unchanged verifier is `12f38a9a33dbbb0a1b7c0d0affcfff91620b927a`: CPU quotient-CG,
all original panels, epsilon factor .6, subdivision0, mu1, center0, residual gate
1e-6 for the original system and both RHS. No smoothing, welding or panel reduction.
FEM rho is an internal diagnostic; neither a normal exit nor a converged BEM solve
proves geometric closure or continuum accuracy. No automatic numerical restarts.

`evaluate_saved_step.py` evaluates an explicitly selected, completely solved
saved state without any optimizer entry point. It verifies the original FEM
commit and HDF5 SHA, uses the same frozen double-precision exporter and BEM
acceptance gate, and waits for the prior Terra BEM's actual exit. The `--step`
selection and `--input-sha256` are mandatory; `--prepare-only` exports and checks
the geometry without launching a solver. Outputs must be fresh. This is the
replacement for the cancelled automatic final-state evaluation when the user
requests an earlier saved state.
