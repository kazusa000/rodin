# October 11 baseline with a graded initial exterior

The source baseline is `3ae2c041`. Only Configuration and the background coordinates in Sphere are changed; the optimization, NodalTransport, Metrics, thickness, trace assembly and framework remain identical to that baseline.

`--graded-chamber --h=0.16666666666666666 --outer-radius=15 --near-field-radius=2 --chamber-gradation=1.3` fixes the algorithmic reference spacing while expanding only the initial exterior grid. The core retains the original coordinate multiplication. Exterior intervals grow by 1.3 up to the configured hmax (10h by default).

MMG retains the original edge-length-driven optimization and 0.1h/10h bounds; this run does not enable `--mmg-adapt` or its different Welsch size map. Grading is prescribed on the initial background; later meshes inherit their edge sizes through the baseline MMG path. It is not an enforced analytic size field at every iteration, and an unchanged h0 does not promise identical surface panels in L2 and L15.

The controller first reproduces both original L2 evaluated states at step 0.1, including a uniform/graded limiting mesh and numerical comparison. Only after these controls pass does it start L15 from the unit sphere at the user-requested doubled step 0.2, regularization 4, thickness 0, MMG, snap 0 and 20 evaluated states (19 real updates). It performs no BEM evaluation or feedback and never restarts a failed stage.

All compilation and computation run on 28h. Configure/build logs, native exit codes, source and binary hashes, control results, complete CSV state counts and the scientific limits are recorded in the external campaign directory. FEM rho is an internal diagnostic.

Use the historical baseline's OpenMP2/BLAS1 for both controls and L15. On 28h the retained original baseline binary at OpenMP6 reproduces the changed second state of the first preflight, while OpenMP2 exactly restores the historical two-state result. The uniform and limiting graded paths agree at either thread count; this is sensitivity to parallel roundoff and subsequent remeshing, not evidence of a changed derivative. Keep this failed preflight as evidence.

`--reuse-build-record` accepts a previously successful compilation only after verifying the numerical source commit, original CMake source path and binary SHA256. Controller-only changes can then reuse that binary without touching its source or build. The controller and numerical source commits are recorded separately, and the full controls are rerun before optimization.
