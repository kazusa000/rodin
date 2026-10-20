# Restarted Cubic with baseline thickness

`run_thickness_resume.py` is the 28h-only, frozen-commit entry point for the
2026-10-19 restart of the original Cubic global10. It does not generate a new
seed. It exports the original double-precision volume mesh, checks exact mesh
coordinates/tetrahedra/vertex and cell references after loading, and requires the unmodified
state solve to recover global10 before adding thickness. The observed reindexing
discrepancy is recorded explicitly (rho about 0.094%); the guard for old FEM
k/c/q/rho is 0.2%, while geometric volume/reference remain within 1e-14. This is
not a bitwise FEM reproduction. Enabling thickness must then leave this
reconstructed state unchanged to relative 1e-8. Coordinates/tetrahedra still
require exact equality. A previous native state-only exit0 may be reused only
with identical input/geometry checks and binary SHA. XDMF snapshots omit
lower-dimensional connectivity: boundary labels are recovered uniquely from
tetrahedral material adjacency and the fixed chamber planes, with every face
and edge label count checked against the original HDF5. These entities are
reindexed; no coordinates or tetrahedra are changed. Feature edges/corners
follow the existing Sphere boundary-feature rules.

`--thickness-method=baseline --thickness-min=2` ports Carlos baseline
`3ae2c0414ffebf11f5814a46510ddf4d8375ba0a`: first outward crossing along inward
geometric normals on all 24 rotated chamber surfaces, nominal distance 2h,
guard distance 2h+2dt=3h, and its bounded paired-point repair. The ray-based load
is a local surrogate, not an exact derivative of the measured penalty. The rho
and thickness directions are separately regularized and projected onto the
volume nullspace. The adaptive multiplier seeks the baseline aggregate repair
rate; its repair-only fallback remains volume-projected. This is a soft
aggregate constraint, not a guarantee of pointwise minimum thickness.

The original h=2/19, dt=.5h, regularization=100h, MMG and transport stay fixed.
The original chamber target volume .17391249230727751 is retained. Six solved
states local0..5 are global10..15, hence five real updates. FEM rho remains an
internal diagnostic. The final complete state receives the existing FP64 CPU
BEM evaluation (full original panels, epsilon factor .6, rtol 1e-6), including
final residual verification with the original NumPy operator. No BEM feedback
enters the updates. Original results and success packages are immutable.

The geometry test covers known unit-cube first-exit distances, nominal/guard
consistency, paired load, inactive thickness, and invalid nominal distance.
The controller also validates separate volume projections and thickness
surrogate descent on the actual restarted global10 before any real update.

## Saved-field MMG diagnosis

`KelvinBall --replay-mmg=<old-optimization-directory> --replay-mmg-count=3`
with the original restart mesh/snapshot and original MMG parameters replays
only the saved `KelvinBall.Chamber.Advected.000000..2.h5` fields. It returns
before the optimization loop: no Stokes solve, gradient, transport update, or
BEM evaluation is started. Before each field is used, coordinates and ordered
tetrahedra must exactly equal its original chamber snapshot; a mismatch stops
the replay instead of transferring the field to another mesh.

For diagnosis only, snapshots and boundary CSVs are written at input,
before-cut, after-cut, before-adapt and after-adapt. The final fixed-plane
check retains its original 1e-6 threshold. Normal invocations have no stage
recording, parameter changes, or relaxed checks. Use a fresh scratch/output
directory and preserve the original run and binary.

`--replay-mmg-lock-cuts` is an isolated causal control: the existing
`Sphere::protectFixedGeometry(..., true)` also requires the cut-plane triangles
during adaptation. It is rejected outside replay. It is not a production fix
or an optimization setting: freezing all cut triangles can retain thin cells,
and does not demonstrate that a long optimization is stable.
