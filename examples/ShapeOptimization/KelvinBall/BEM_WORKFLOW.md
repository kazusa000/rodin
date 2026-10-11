# Single-checkout KelvinBall + BEM workflow

All source belongs to the Rodin fork's `codex/clement-swift-bem` branch.
The native KelvinBall solver, BEM/CUDA operator, six outer-wall velocity
families, acceptance controller, diagnostics and viewers are distributed in
this example. No Project03 checkout, nested BEM repository or binary-only
dependency is required.

## Source and unchanged method

The native solver descends from Carlos's `example/KelvinBall@79718dc`;
the existing integrated native source is unchanged by this migration.
Python workflows come from Project03 `067f8e9`, while the complete BEM and
CUDA operator come from the actually used verifier `12f38a9`, not the
different historical scorer at Project03 HEAD. Source/function checksums are
in [migration-provenance.json](reproduction/migration-provenance.json).
109 original numerical definitions retain identical ASTs after removal of
import declarations; the CUDA arithmetic body is unchanged. Only package
imports, source/output paths, deployment and source-version bookkeeping change.

The method remains n13 cubic O24, a fixed `0.02 H3 H6` perturbed sphere,
six BEM wall traces, FEM internal states and strain direction, volume-nullspace
projection, rotated characteristic transport and MMG/SWIFT reconstruction.
Soft thickness is `0.06 R0`, weight 1, with the existing exact discrete
derivative. Candidates require the previous positive-area full-solid
connectivity gate and actual BEM `rho-P` Armijo acceptance. Carlos's new
thickness correction is not enabled simultaneously. No additional shape
constraint, seed selection, seam repair, panel reduction or MoveMesh is added.

The retained protocol is all original O24 panels, epsilon factor 0.6, FP64,
block 256, 2000 iterations and one hour per stage. BEM uses `1e-6`, with at
most one `1e-9` refinement when the weak-coupling threshold `1e-4` is met.
Both CPU and GPU true residuals must meet the effective recorded tolerance.
Numerical parameters are in [clement-swift-bem.json](reproduction/clement-swift-bem.json).
Node configuration cannot override them.

The FEM proxy direction is not an exact BEM derivative: the previous centered
BEM FD check had the same sign but approximately 69 times smaller magnitude.
This limitation remains unchanged. Open interface exports are original-panel
BEM diagnostics, not closed-body physical or printing acceptance.

## Layout

```text
KelvinBall/
  *.cpp, *.h                  Existing native Carlos solver and BEM interfaces
  python/kelvin_bem/          BEM, IO, geometry, acceptance and orchestration
  cuda/bem_quotient_cuda.cu   Optional FP64 GPU operator source
  scripts/kelvin-bem.py       One entry point, usable from any working directory
  scripts/templates/         Lossless full-resolution rotating/slider viewer
  tests/                     Migrated tests and single-checkout regression checks
  input/                     Small versioned sphere test mesh
  reproduction/              Numerical preset and migration provenance
```

Builds, fields, logs, shared libraries and generated HTML belong outside the
source checkout. The existing workspace can keep using `tmp/project-03/`
for builds and `out/project-03/clement-swift/` for results. Old Project03
repositories, worktrees and experiments remain retained and unmodified.

## Installation

Clone the existing fork, select the branch and initialize the required native
submodules. Native build dependencies follow the Rodin installation guide:
C++20 compiler, CMake, Eigen, Boost, HDF5, MUMPS and its support libraries.
The validated deployment uses MUMPS 5.6.2; migration does not upgrade it.
Python 3.11 or newer is required (`hashlib.file_digest` is used).
Create an isolated environment and install this example's requirements.
CUDA is optional for CPU tests and explicitly required for GPU measurement.

```sh
git clone --branch codex/clement-swift-bem https://github.com/kazusa000/rodin.git
cd rodin
git submodule update --init third-party/mmg third-party/termcolor
python3 -m venv /path/outside/rodin/venv
/path/outside/rodin/venv/bin/python -m pip install -r examples/ShapeOptimization/KelvinBall/requirements.txt
python3 examples/ShapeOptimization/KelvinBall/scripts/kelvin-bem.py --help
```

The controller and CPU/GPU checkouts must use the same clean commit.
Native build records include the source commit, submodule versions, commands,
binary hashes and runtime library paths. BEM records include that same source
commit, input and boundary-table hashes, CUDA source/library hashes and both
true residual certificates. The CUDA library is built from the bundled source
and verifies its embedded source hash; the legacy binary hash is attribution,
not an installation requirement.

## Native and CUDA builds

Commands below use generic paths; configure dependencies for the local machine
through `--cmake-arg` and its existing isolated environment. Outputs must be
new directories, never an old result/build tree. Native source defaults to
this checkout; another native repository is rejected.

```sh
python examples/ShapeOptimization/KelvinBall/scripts/kelvin-bem.py native build \
  --build-dir /path/outside/rodin/native-build --jobs 2 \
  --cmake-arg=-DCMAKE_PREFIX_PATH=/path/to/dependency-prefix \
  --cmake-arg=-DMUMPS_ROOT=/path/to/dependency-prefix
python examples/ShapeOptimization/KelvinBall/scripts/kelvin-bem.py native checks \
  --build-dir /path/outside/rodin/native-build --output /path/outside/rodin/native-checks
python examples/ShapeOptimization/KelvinBall/scripts/kelvin-bem.py build-cuda \
  --nvcc /path/to/cuda/bin/nvcc --output /path/outside/rodin/libbem_quotient_cuda.so
```

CUDA compilation retains compute capability 7.5, FP64, `--fmad=false` and
the original compiler flags. GPU selection never silently falls back to CPU.

## Measurement and paired acceptance

```sh
KELVIN_SYMMETRY=cubic python examples/ShapeOptimization/KelvinBall/scripts/kelvin-bem.py bem \
  --mesh /path/to/accepted.mesh --output /path/outside/rodin/new-measurement \
  --cuda-library /path/outside/rodin/libbem_quotient_cuda.so
```

For the distributed workflow, copy [nodes.example.json](reproduction/nodes.example.json)
to an external deployment configuration and set hosts, interpreter/checkouts,
build and unique output paths. Explicit CLI deployment flags take precedence.
The fixed preset continues to govern the numerical method.

```sh
KELVIN_SYMMETRY=cubic python examples/ShapeOptimization/KelvinBall/scripts/kelvin-bem.py paired-check \
  --nodes-config /path/outside/rodin/nodes.json --output /path/outside/rodin/controller-results
```

`paired-check` retains the validated bounded campaign: one BEM-accepted update
for each backend from a shared input, centered FD diagnostics and same-target
reconstruction. It is not a long optimization or a T12 migration. `--mesh`
can supply an existing accepted checkpoint; omitting it generates the fixed
perturbed sphere. The `native step` and `bem` stages are the same building
blocks for subsequent updates. No historical 03 optimization driver is
transplanted as a second native solver.

The controller checks its own and both workers' commit equality, keeps actual
failure/exit evidence and backtracks only the previously identified
candidate-local fixed-boundary or BEM iteration-limit failures. Unknown
source/import/transport failures abort; no tolerance relaxation or infinite
restart is introduced. Unique outputs retain rejected candidates as evidence.

## Diagnostics and viewer

`export-chamber` restores exact saved chamber coordinates, `export-target`
preserves the original saved nodal target, and `audit` compares reconstruction
against that target. `viewer` retains the rotating, iteration-slider HTML
format and lossless full-panel data. Accepted history chains and old result
metadata remain readable without importing code from their original 03 tree.
For remote viewer export, `--reader-source` points to this Rodin checkout,
not to an old Project03 reader. No HDF5 is downloaded by that viewer workflow.

Each command provides `--help`. Tests run from the single checkout:

```sh
KELVIN_SYMMETRY=cubic python -m pytest -q examples/ShapeOptimization/KelvinBall/tests
KELVIN_SYMMETRY=cubic KELVIN_RUN_CUDA=1 BEM_CUDA_LIBRARY=/path/to/libbem_quotient_cuda.so \
  python -m pytest -q examples/ShapeOptimization/KelvinBall/tests/test_bem_cuda.py
```

GPU tests are explicitly skipped without hardware opt-in, not marked as passed.
Migration acceptance also requires a same-input native/BEM campaign and old
saved stopping-case replay, in addition to these tests. Validation results
and the exact executing commit are recorded separately from subsequent
documentation-only commits.

The completed migration checks and retained-data comparison are documented in
[MIGRATION_VALIDATION.md](MIGRATION_VALIDATION.md), with compact public evidence
in [migration-validation.json](reproduction/migration-validation.json).
