# Tetrahedral extension of the September 22 implementation

`KelvinBallTetra` reuses the same six-state, full-rho derivative, null-space
volume projection, scalar trace identification, characteristic transport and
MMG reconstruction as `KelvinBall`. No assembly or derivative term is changed.
The original cubic executable remains the default, including its sphere seed.

The new chamber is `x >= |z|, y >= |z|`, truncated by the same cube of half-width
L. Twelve determinant-one signed permutation matrices with positive sign
product reconstruct the full domain. The cut labels are specific to this mode:

| Label | Plane | Pair |
| --- | --- | --- |
| 31 | y=z | 31 -> 32 via (z,x,y) |
| 32 | x=z | inverse of the preceding map |
| 33 | x=-z | 33 -> 35 via (y,-z,-x) |
| 35 | y=-z | inverse of the preceding map |

There is no x=y boundary or self-paired half-plane. The two outer cube faces
share the no-slip label 7, while their intersection is marked as a ridge by
distinct local feature identifiers. Those identifiers do not alter mesh labels.
The fixed-boundary tolerance remains 1e-6.

The initial radial function on unit directions is
`r=1+0.02*(H3+H6)/sqrt(2)`, with `H3=3sqrt(3)xyz` and
`H6=6sqrt(3)(x²-y²)(y²-z²)(z²-x²)`. This is one fixed sphere perturbation,
not a score-selected initial design. Its group invariance, breaking of an
additional cubic quarter-turn and reflection, and cut maps are checked by
`KelvinBallTetraGeometryTest`; `KelvinBallCubicGeometryTest` checks the default.
Only the MMG route is supported for the tetrahedral extension.

For the requested experiment use `--thickness-min=0`: no extra thickness
constraint. FEM values remain bounded-domain diagnostics; free-space BEM and
seam diagnostics must be recorded separately before physical conclusions.
