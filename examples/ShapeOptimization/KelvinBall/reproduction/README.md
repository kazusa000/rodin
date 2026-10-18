# Original-core L15 control

This branch starts from the complete Rodin `9b3dfd7e` used in the successful
small-box experiment. `src/`, the Stokes and gradient assembly, distance
projection, transport, thickness load, and reconstruction algorithm are
unchanged. Only `Configuration` and the sphere chamber/size-map construction
gain an opt-in graded exterior; the old uniform path remains the default.

Use `--graded-chamber --h=0.10526315789473684 --outer-radius=15`.
The old core `[0,2]` has the same 20 axis coordinates, computed as `i*h`.
Only the exterior intervals grow by 1.3, capped at `10h`. The MMG target is
`h` throughout the old core, then increases linearly to `10h` outside.
The reference `h=2/19` is independent of the outer extent, including MMG retries.

`run_original_l15.py` runs only on 28h. It builds from this committed source,
runs the configuration test, and checks two L2 initial states (the default
uniform path and the graded path with no exterior). Both must reproduce the
historical k/c/q/volume and identical HDF5 mesh datasets. It then runs the
L15 case from the original analytic unit sphere for 20 evaluated states,
with reg100h, dt0.5h, thickness2h, quadrature3, snap0.005, and the original MMG
defaults. Outputs and actual exit codes are external to the source tree.
There is no BEM evaluation or boundary feedback, restart mesh, added seed,
surface smoothing, new locator, or coverage projection in this entry.

The enlarged-domain initial MMG mesh is allowed to differ from L2: the sphere
and near-field target resolution are fixed, but global remeshing is not a
guarantee of identical initial interface triangles. The historic locator's
relative acceptance rule is unchanged. FEM rho remains an internal diagnostic,
and a completed finite iteration budget is not a convergence or independent
BEM certificate. All former results and source variants remain preserved.
