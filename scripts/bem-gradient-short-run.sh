#!/usr/bin/env bash
# Native C++ short trajectory; reject existing output and retain actual exit status.
set -u
if [[ $# != 2 || "$1" != /* || ! -x "$1" || "$2" != /* || -e "$2" ]]; then
  printf 'Usage: %s ABSOLUTE_EXECUTABLE NEW_ABSOLUTE_RESULTS_DIRECTORY\n' "$0" >&2
  exit 2
fi
executable=$1
results=$2
mkdir -p "$results" || exit 2
printf '%s\n' "$executable --seed-amplitude=0.6 --step-factors=0.03 --updates=3 --volume-tolerance=0.005 --max-backtracks=6" > "$results/command.txt"
(
  cd "$results" || exit 2
  OMP_NUM_THREADS=1 "$executable" --seed-amplitude=0.6 --step-factors=0.03 --updates=3 --volume-tolerance=0.005 --max-backtracks=6 > run.log 2>&1
)
result=$?
printf '%s\n' "$result" > "$results/run.exit"
exit "$result"
