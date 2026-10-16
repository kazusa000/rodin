#!/usr/bin/env bash
# Resume native BEM from an accepted mesh; preserve original volume and real exit.
set -u
if [[ $# != 5 || "$1" != /* || ! -x "$1" || "$2" != /* || -e "$2" || "$3" != /* || ! -f "$3" ]]; then
  printf 'Usage: %s ABSOLUTE_EXECUTABLE NEW_ABSOLUTE_RESULTS ACCEPTED_MESH ORIGINAL_VOLUME ACCEPTED_RHO\n' "$0" >&2
  exit 2
fi
executable=$1
results=$2
checkpoint=$3
target_volume=$4
expected_rho=$5
mkdir -p "$results" || exit 2
arguments=(--seed-amplitude=0.6 --step-factors=0.03 --updates=100 --volume-tolerance=0.005 --max-backtracks=6
  "--initial-mesh=$checkpoint" "--target-volume=$target_volume" "--expected-rho=$expected_rho")
printf '%q ' "$executable" "${arguments[@]}" > "$results/command.txt"
printf '\n' >> "$results/command.txt"
(
  cd "$results" || exit 2
  OMP_NUM_THREADS=1 "$executable" "${arguments[@]}" > run.log 2>&1
)
result=$?
printf '%s\n' "$result" > "$results/run.exit"
exit "$result"
