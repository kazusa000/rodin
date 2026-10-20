# October 11 baseline with a graded initial exterior

The source baseline is `3ae2c041`. Only Configuration and the background coordinates in Sphere are changed; the optimization, NodalTransport, Metrics, thickness, trace assembly and framework remain identical to that baseline.

`--graded-chamber --h=0.16666666666666666 --outer-radius=15 --near-field-radius=2 --chamber-gradation=1.3` fixes the algorithmic reference spacing while expanding only the initial exterior grid. The core retains the original coordinate multiplication. Exterior intervals grow by 1.3 up to the configured hmax (10h by default).

MMG retains the original edge-length-driven optimization and 0.1h/10h bounds; this run does not enable `--mmg-adapt` or its different Welsch size map. Grading is prescribed on the initial background; later meshes inherit their edge sizes through the baseline MMG path. It is not an enforced analytic size field at every iteration, and an unchanged h0 does not promise identical surface panels in L2 and L15.

The controller first reproduces both original L2 evaluated states at step 0.1, including a uniform/graded limiting mesh and numerical comparison. Only after these controls pass does it start L15 from the unit sphere at the user-requested doubled step 0.2, regularization 4, thickness 0, MMG, snap 0 and 20 evaluated states (19 real updates). It performs no BEM evaluation or feedback and never restarts a failed stage.

All compilation and computation run on 28h. Configure/build logs, native exit codes, source and binary hashes, control results, complete CSV state counts and the scientific limits are recorded in the external campaign directory. FEM rho is an internal diagnostic.

Use the historical baseline's OpenMP2/BLAS1 for both controls and L15. On 28h the retained original baseline binary at OpenMP6 reproduces the changed second state of the first preflight, while OpenMP2 exactly restores the historical two-state result. The uniform and limiting graded paths agree at either thread count; this is sensitivity to parallel roundoff and subsequent remeshing, not evidence of a changed derivative. Keep this failed preflight as evidence.

`--reuse-build-record` accepts a previously successful compilation only after verifying the numerical source commit, original CMake source path and binary SHA256. Controller-only changes can then reuse that binary without touching its source or build. The controller and numerical source commits are recorded separately, and the full controls are rerun before optimization.


## Saved-state continuation

`--restart-from <completed-L15-campaign> --updates 20` resumes global19 and requests 20 new real updates, through global39. Local0 is a repeated evaluation of global19 and does not count as an update. L15, h0, MMG, step0.2, regularization4 and OpenMP2/BLAS1 remain fixed; no BEM is launched.

The saved XDMF mesh has double coordinates and tetrahedra, but does not persist lower-dimensional connectivity. `restart_snapshot.py` retains every vertex and original tetrahedron/material in order, recovers the existing solid/fluid interface and preserves the original matching native fluid boundary labels; missing solid-cut labels use original vertex references or a whole-triangle halfspace with the fixed-plane roundoff tolerance. It requires all label counts and ridge counts to match the snapshot, and serializes coordinates with 17 significant digits and an exact float64 round-trip check. It never moves, welds, remeshes or simplifies a point. Missing or ambiguous boundaries are errors. The native MMG entry bypasses sphere initialization and reads this volume mesh; required MMG features are recovered by the unchanged normal reconstruction path.

Before continuing, the controller repeats the small-box regression, then fully evaluates the imported global19 state and requires k/c/q/volume to match the prior final row. The original Theta direction and Distance fields must also match their saved counterparts. Its exported vertices, tetrahedron ordering and materials must be bit-identical to the old HDF5 snapshot; vertex/edge/face label histograms must match. The original initial target volume remains the CSV reporting reference. The direction still enforces zero first-order volume change; it does not correct accumulated finite-step volume drift. The numerical update loop from the prior baseline is checked byte-for-byte unchanged.
