#!/usr/bin/env bash
# Build and run the tetrahedral KelvinBall from a frozen checkout on 28h.
set -eu
root=${1:?workspace root required}
src=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
sha=$(git -C "$src" rev-parse HEAD)
test -z "$(git -C "$src" status --porcelain)"
tmp=$root/tmp/project-03
name=L10-n50-tetra-20261122-${sha:0:8}
resume=()
if [[ $# -eq 2 || ( $# -eq 3 && $3 == amd ) ]]; then
  parent=$(realpath "$2")
  test -r "$parent/kelvin-ball.csv"
  name=L10-n50-tetra-continue-${sha:0:8}
  resume=("--continue-from=$parent")
  if [[ $# -eq 3 ]]; then
    resume+=("--swift-linear-mumps-amd=1")
  fi
elif [[ $# -ne 1 ]]; then
  echo "Usage: $0 workspace-root [parent-optimization-directory [amd]]" >&2
  exit 2
fi
out=$root/out/project-03/fem-domain/$name
build=$tmp/$name-build
mkdir "$out"
trap 'echo "$?" > "$out/exit-code.txt"' EXIT
printf '%s\n' "$sha" > "$out/source-commit.txt"
dep=$tmp/carlos-bem-parallel-20261004-deps/prefix/usr
lib=$dep/lib/x86_64-linux-gnu
blas=$tmp/l10-openblas-lp64-20261105/prefix/usr/lib/x86_64-linux-gnu/openblas-pthread
export LD_LIBRARY_PATH="$blas:$lib:$lib/openmpi/lib:$lib/hdf5/serial${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=7
stage() (
  name=$1
  shift
  mkdir "$out/$name"
  cd "$out/$name"
  trap 'echo "$?" > exit-code.txt' EXIT
  printf '%q ' "$@" > command.txt
  printf '\n' >> command.txt
  "$@" > run.log 2>&1 &
  pid=$!
  echo "$pid" > pid.txt
  wait "$pid"
)
stage configure cmake -S "$src/examples/ShapeOptimization/KelvinBall/reproduction" -B "$build" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$dep" -DMUMPS_ROOT="$dep" \
  "-DMPI_C_COMPILER_INCLUDE_DIRS=$lib/openmpi/include;$lib/openmpi/include/openmpi" \
  -DMPI_C_HEADER_DIR="$lib/openmpi/include" -DMPI_C_LIB_NAMES=mpi \
  -DMPI_mpi_LIBRARY="$lib/openmpi/lib/libmpi.so" -DKELVIN_TETRAHEDRAL=ON
stage compile cmake --build "$build" --target KelvinBall KelvinBallTetraTest -j4
stage symmetry "$build/kelvin/KelvinBallTetraTest"
sha256sum "$build/bin/KelvinBall" > "$out/binary.sha256"
ldd "$build/bin/KelvinBall" > "$out/libraries.txt"
stage optimization "$build/bin/KelvinBall" \
  --n=50 --iterations=200 --outer-radius=10 \
  --hmin-factor=0.05 --hmax-factor=50 \
  --regularization=20 --step=0.2 \
  --reconstruction=swift --mmg-adapt --thickness-min=0.05 \
  --swift-linear-solver=mumps --swift-linear-threads=7 \
  --swift-model-fit=10 \
  --swift-convergence-iterations-outer=30 --swift-convergence-iterations-inner=15 "${resume[@]}"
