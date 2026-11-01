# Carlos upstream SWIFT, direct L15 experiment

Base: `a9f1115ae71525978f086917a73958f95d77d1bf` (`example/KelvinBall`, 2026-10-31).
This entry runs only on the authorized host `28h`. No small-box FEM controls are run.
The synthetic checkpoint fixture is an I/O and coordinate-generation test, not an optimization experiment.

The numerical changes from upstream are confined to the optional graded initial
background in Configuration/Sphere. The optimizer additionally supports opt-in
native checkpoint I/O. One missing `typename` before `decltype(classifier)::Parameters`
is a GCC compatibility correction, not a change to the SWIFT algorithm. The standalone
build explicitly links the packaged OpenMPI C++ compatibility library required by
its MUMPS headers; the numerical executable still runs serial MPI / OpenMP4. Framework
`src/`, gradients, Stokes solves, rotated coupling, transport, and fitting remain upstream.

## Fixed experiment

- Unit sphere; box halfwidth 15, full side 30.
- Reference h0=1/6; inner halfwidth 2 retains coordinates i*h0.
- Initial outer intervals grow by 1.3, capped at 10h0. Subsequent SWIFT/MMG
  adaptation uses the upstream algorithm, not a newly imposed analytic size field.
- SWIFT with `--mmg-adapt`, hmin=.1h0, hmax=10h0; regularization=20h0,
  step=.2h0, absolute minimum thickness=.05.
- MUMPS, SWIFT outer/inner limits 30/15, OpenMP4, BLAS1.
- First six evaluated states 0..5 are five real updates. A separate one-state
  readback of state5 adds no update. Only after checkpoint, full-field consistency,
  independent BEM acceptance, and positive net FEM diagnostic change does the
  controller resume for 15 further updates (local0=global5; local15=global20).
  Net FEM change is only a conservative workflow gate, not physical validation.
- Independent BEM on saved steps5/20: full original interface, unchanged double
  coordinates, cubic O24, epsilon_factor=.6, subdivision0, mu1, center0,
  CPU quotient-CG, rtol1e-6 (total and both RHS), maxiter6000, block256, row_chunk32.
  The immutable exporter/validator and verifier are pinned by the controller.
  BEM does not feed back into FEM. Fixed regularized BEM is not continuum truth.

## Persistence and failure behavior

Each `checkpoint-NNNNNN` contains native full topology, coordinates, materials,
all-dimensional labels, and explicit MMG corner/ridge/required entity tags for
**both** the fitted chamber and the evolving SWIFT background. Binary scalar state
retains the original volume reference. `READY` means persistence completed; only a
full CSV row proves that a state was solved. XDMF alone is not used to restart.
The readback gate compares every checkpoint dataset exactly, k/c/q/rho/volume and
original target volume at rtol1e-8/atol1e-10, full Theta at 1e-7/1e-9 and full Distance
at 1e-12/1e-12. A mismatch stops the queue without relaxing tolerances.

Every stage has its real child exit, command, timestamps and log. Compilation
budget is 4h, pilot/continuation each 8h, restart validation 2h, each BEM72h.
No automatic retry, parameter change, geometry repair or small-box fallback occurs.
Old results remain untouched. The controller refuses existing build/output paths.

Run `run_l15_swift.py --help` with the scientific venv on 28h; supply clean pinned
source/commit, fresh external build/output paths, the recorded dependency prefix,
immutable utility checkout b8f3c8c and verifier12f38a9. Launch through a unique tmux
wrapper recording its own exit independently of the campaign exit.
