# L10 recovery after the 2026-11-06 host reboot

The parent run preserves CSV states 0..23 and complete chamber, fluid, distance,
velocity, pressure and direction files for state 24. It died during reconstruction
towards 25. Its independent pre-fit SWIFT background was never written.

`run_l10_continue.py` recovers the original state-24 tetrahedron order, materials
and double coordinates, using the matching last native fluid mesh for boundary
labels. The solid cuts use the original vertex references, with original label
histograms required. Recovery does not move, refit, smooth, weld or simplify the
saved geometry. MEDIT supplies topology; coordinates are overwritten directly
from the original binary HDF5 dataset before evaluating anything.

The missing independent SWIFT background is explicitly reinitialized as a copy
of the saved fitted chamber. This is a new continuation trajectory from the saved
shape. It is not a replay of the missing pre-fit background. Parameters, original
target volume, framework and numerical update loop remain unchanged.

Before updates, a single-state reevaluation must reproduce exact chamber/fluid
coordinates, ordered tetrahedra and materials, lower-dimensional label
histograms, every saved velocity/pressure field and the complete Theta and
Distance fields. Field tolerances are 1e-7/1e-9 (rtol/atol), with 1e-12/1e-12 for
Distance. Missing fields, changed labels, changed coordinates or failed comparisons
stop the controller. The recovered state's new full-precision k/c/q/rho are saved;
state 24 had no complete parent CSV row.

The continuation evaluates local 0..76, corresponding to global 24..100, with 76
new updates after the recovered state. Local 0 repeats the accepted recovery state.
There is no BEM, wall-clock cutoff, recurring monitor or automatic retry.

Every iteration starts by saving full native chamber/background topology,
double coordinates, attributes, transformations and six MMG marker sets.
`checkpoint-NNNNNN/READY` marks completed persistence, not completed evaluation.
The binary `/State` vector is: format version (1), global index, original target
volume, background h, previous rho, predicted rho change, previous volume,
previous thickness penalty, predicted thickness change, predicted volume change,
initial h0, regularization factor, step factor, minimum thickness, hmin, hmax.
Absent prior scalars use NaN. Checkpoint I/O must pass a synthetic native dataset
roundtrip and reject out-of-range MMG marker indices before recovery evaluation.
