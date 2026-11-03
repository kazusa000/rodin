# Carlos L10 reproduction, 3 November 2026

The experiment uses upstream `d948427a50c00264c81ca40b2f616220581c7a26`,
including the bounded thickness correction published on 3 November. Carlos's
email specifies `--outer-radius=10 --n=60`. The remaining numerical options
follow the working SWIFT command in the upstream README: 200 evaluated states,
regularization factor 20, step factor 0.2, absolute minimum thickness 0.05,
MMG adaptation, size factors 0.1/10, MUMPS with four threads, and SWIFT budgets
30 outer / 15 inner iterations. The unmodified initial condition is a unit sphere.

The initial background is the upstream uniform grid, with spacing `10/59`.
The physical step is `2/59`, and the regularization length is `200/59`.
There is no graded-domain extension. All numerical sources are upstream;
the only C++ edit supplies a dependent-type `typename` required by GCC.
The standalone build also links the packaged OpenMPI C++ compatibility library.

`run_carlos_l10.py` executes only on 28h. It requires fresh external build and
output directories, clean exact source commits, and initialized pinned MMG and
termcolor submodules. Configure and compile each have a four-hour budget;
optimization and the final independent BEM each have a 72-hour budget.
No small-box FEM control, restart, or automatic parameter change is performed.

The final shape is the last complete, contiguous FEM CSV state with its
matching original HDF5 geometry. A later saved but unsolved shape is excluded.
If the native run stops early, its actual exit is retained and the final complete
state is still assessed; completion of the requested budget remains false.
200 evaluated states correspond to 199 real updates.

The fixed utilities commit is `b8f3c8c5929c21aa9bd72236c0c1e843e0a72a64`,
and the BEM verifier is `12f38a9a33dbbb0a1b7c0d0affcfff91620b927a`.
The original full HDF5 interface is exported with double precision and unchanged
panels. The protocol uses cubic O24 symmetry, epsilon factor 0.6, subdivision 0,
viscosity 1, center 0, CPU quotient CG, relative tolerance `1e-6`, and at most
6000 iterations. The actual solver exit, convergence, total residual, both RHS
residuals, source, mesh hash, and full-panel count are checked before acceptance.
BEM results are retained separately from the internal FEM score.

Run the controller selection tests on 28h with:

```sh
cd /path/to/source/examples/ShapeOptimization/KelvinBall/reproduction
python -B -m unittest -v test_run_carlos_l10.py
```

The launch record supplies `--source`, `--commit`, `--build`, `--output`,
`--prefix`, `--utilities`, and `--verifier` as absolute paths. Each stage records
its command, actual process exit, and log. Original geometry and open symmetry
seams are preserved. The experiment measures a fixed discretization and BEM
regularization; the iteration budget is an operational limit.
