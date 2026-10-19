# Native infinite-exterior BEM on the original O24 chamber

The design mesh is the Carlos chamber `0 <= |z| <= y <= x <= L`. The frozen
baseline is `3ae2c0414ffebf11f5814a46510ddf4d8375ba0a`. Only the state model
changes from finite-box Stokes FEM to the free-space regularized single layer.
The outer box is a fixed boundary of the H1/level-set extension, not a no-slip
fluid wall in this BEM model. No Python state solver or FEM proxy gradient is used.

## State space and native assembly

`BEM::evaluateChamber` extracts only Gamma within the chamber. It does not sew
volume copies and does not call the full-surface solver. `BEMQuotientKernel`
acts on 3x3 responses:

```
K_Q[F](x,y) = sum_{g in O24} S_epsilon(x-g*y) g F(y) g^T.
```

The load index rotates as well as the physical-vector index. Translation and
rotation are two matrix-valued families (nine components each), assembled with
Rodin `P1(surface,3,3)`, `Integral(Potential(kernel,trial),test)` and `DenseProblem`.
`ChamberSymmetry::tensorBasis` constructs E from vertex orbits and stabilizers,
including symmetry axes, and the native solver uses `E^T A E` and `E^T b`.
Residuals are measured in the admissible test space. No image unknowns are stored.

If B and X contain these two families, C=B^T X. The three isotropic resistance
coefficients are `k=8 C00`, `q=8 C11`, `c=4(C01+C10)`, accounting for 24 patches
and three load columns. `rho=abs(c)/sqrt(k*q)`. Zero coupling has no absolute-value
gradient and is rejected when a derivative is requested. The epsilon model is
`0.6*sqrt(chamber_area/chamber_panels)`, equal to the full-image area/panel ratio.

## Coordinate derivative and chamber update

The discrete differential is `2 <dB,XW> - <XW,dA X>`. It includes both integration
measures, both image coordinates (`dr=dx-g*dy`), the rotational load, and the
area-dependent epsilon. It is deposited directly into chamber Gamma vertices
and then into the original chamber shape DOFs. No full-body coordinate gradient
or finite-difference direction is used by the optimizer.

The Carlos H1 metric, rotated Nitsche coupling, volume-nullspace projection,
Eikonal and rotated RK4 characteristic continuation remain in use. Scalar/vector
representative projections now use chamber orbits, without constructing a sewn
volume. Volume is reported and projected to first order; there is no volume-drift
acceptance cap, correction or post-hoc rescaling.

For the native MMG path, transport and cutting use a fixed copy of the original
uniform chamber background. Both paired traces therefore share the same scalar
nodal data and initial trace triangulation. Distance and direction are evaluated
spatially from the fitted chamber; nodal copying between the two meshes is not
used. Each MMG trial cuts a fresh background copy. The original MMG cutter,
material checks and fixed-plane projection are retained. This is an explicit
change to the reconstruction input, introduced to avoid independent remeshed
cut traces. It does not freeze the physical Gamma/cut intersection across updates.

Virtual image edge incidence must be two everywhere before BEM assembly.
The audit stores coordinate keys and edges, not image meshes. Open edges are
errors; tolerance is not enlarged, faces are not filled, and the check is never
bypassed. Exported sewn fields are visualization only.

## Verification and bounded use

`KelvinBallBEMChamberValidation --reference` checks states, an independent
full-surface-only reference on the same rotated patch, and centered finite
differences on rotation-compatible chamber coordinate directions. The reference
is a test oracle, not the production state path and never a full-volume optimizer.
The derivative test uses three halved steps with a final relative-error target
of 0.2%; the reference target is 1e-8, constrained state residual 1e-9, and zero
unmatched image edges. These are discrete consistency gates, not mesh convergence.

```
KelvinBall --bem-gradient --n=13 --seed-amplitude=0.1 --iterations=2 \
  --step=0.10 --reconstruction=mmg --mmg-snap=0.001 --mmg-retries=0 --save-mesh
```

Two evaluations request one actual update. `dt=0.10*L/(n-1)`, not the full-box
fixture spacing. The bounded validation must also inspect fresh BEM rho after
MMG, compare it with the Armijo prediction, and verify the chamber checkpoints.
No long campaign is launched by these checks. SWIFT, stability/cost studies and
long optimization remain separate validation tasks. Historical full-volume
fixtures retain their old purpose and must not be used as chamber acceptance.

The bounded update check enables the existing Carlos near-vertex snapping
parameter at 0.001 of a crossed edge, to avoid tiny tetrahedra when the zero
contour nearly hits a background node. This is a recorded reconstruction
parameter, not a volume constraint or seam-welding tolerance. A reconstruction
failure ends that check immediately (`mmg-retries=0`); it is never accepted by
relaxing the fixed-plane or positive-Jacobian checks.

During the native MMG level-set cut only, `setPreserveVolumeBoundary` enables
MMG3D nosurf/nomove. Required input faces intersected by the level set may be
replaced; without these temporary policies MMG can subsequently round the
chamber planes. Level-set intersections are still generated afresh each update.
After cutting, the existing optimizer can move interior vertices while required
cut faces and vertices are protected. Fixed-plane/Jacobian and image-edge audits
remain mandatory. The default FEM cutter and 2D/surface kernels retain their
existing policies. Exact cut inputs and pre-projection candidates are retained.

## Verified checkpoint (2026-10-18)

On 28h, numerical commit `379bd613f6b97d65aca4d81b0deaa4288908245c`
completed `scripts/bem-gradient-chamber-check.sh` with build, state/FD, native
update and combined acceptance exits all zero. The coarse reference/FD check
has 36 chamber panels, 243 matrix-state DOFs and 164 independent DOFs. Its
maximum resistance/rho relative difference from the full-surface reference is
`5.618068023815625e-13`; the final centered-FD relative error is
`3.191422722907536e-7`, decreasing by approximately four on each halved step.

The n=13 update uses the command above. Rho changes from
`6.0418514161355301e-10` to `3.4841093740257715e-8`; both states have zero
unmatched image edges and residuals below `3e-16`. The projected volume slope
is `-8.6736173798840355e-18`. The chamber volume changes by `-0.8463987521%`;
it is recorded, without a drift rejection cap. The mesh passes the existing
Jacobian/fixed-plane checks; final minimum tetrahedron quality is only
`0.004457135467049566`, so this is not evidence of robust long optimization.

The first transported scalar has rotated jump `1.11e-15` and zero fallback
vertices. After reconstruction, Eikonal traces on the fitted, nonmatching cut
triangulations have jump `0.04942269014`; the next transport would be evaluated
and projected onto the common chamber background. No second update was run.
The first-order rho prediction is `4.177957982374644e-7`, whereas the actual
increase is `3.4236908598644163e-8`. Weak Armijo passes, but remeshing changes
both panels and area-based epsilon. These tests establish discrete consistency
and one actual update, not a mesh-independent physical improvement, convergence,
second-step seam accuracy, numerical stability or computational cost.

Remote evidence (retained in place):
`/home/profil/p90318/work/main/m2_internship/out/project-03/bem-gradient/native-chamber-20261018-379bd613`.
Earlier failed candidates and diagnostic snapshots are preserved separately.

## Bounded continuation

`scripts/bem-gradient-chamber-run.py --executable ABSOLUTE_KELVINBALL
--results NEW_ABSOLUTE_RESULTS --checkpoint ACCEPTED_CHAMBER_MESH
--expected-rho ACCEPTED_RHO --updates 20` continues from an accepted chamber
checkpoint. Each child process performs exactly one native update with the
validated n=13 and MMG/snap parameters. Step and H1 length factors are
explicit `--step` and `--regularization` options (defaults 0.10 and 4).
Both must be finite and positive, are recorded in history and passed unchanged
to the native child. The step=1.0, H1=4 continuation starts from the last
accepted chamber mesh, not the separately diagnosed larger-step candidate. The Python driver only launches
native processes, checks their JSON/CSV outputs and copies accepted meshes.
It neither solves states nor constructs gradients, transport or geometry.

The same output checker is used by the one-step check. Every restart must
reproduce the preceding accepted rho. No candidate is copied to accepted-N
until all checks pass; rejected candidates remain in their trial directory.
There is no volume-drift rejection and no failed-campaign automatic restart.
The driver records a bounded continuation result, not a convergence claim.

Optional `--max-backtracks 2` requests the initial step and at most two halvings:
for `--step 1.0`, the factors are exactly 1.0, 0.5 and 0.25. Each candidate
starts from the same accepted checkpoint and recomputes the native direction;
a rejected mesh never becomes the next trial input. Objective rejection or a
nonzero native candidate result permits the next smaller step. Native return
codes remain in each trial; if every new candidate has a native failure, the
last actual native failure is propagated. Broken restart, projection or output
contracts stop immediately. Only a successful candidate advances the accepted
count. Three rejected candidates exhaust that iteration and stop the campaign;
this does not establish convergence. The next accepted iteration starts again
at the requested initial step. The default of zero halvings preserves the
previous fixed-step behavior and output paths.

`--first-backtrack 1` may skip an already documented initial-step rejection
on the first iteration only. The launch record must identify the prior failed
trial and verify the checkpoint, parameters and numerical binary hashes.
It records no new measurement at the skipped step. Later iterations restart
at trial zero. In line-search mode, files are under
`iteration-N/trial-M/update`, with validation alongside each update directory.
