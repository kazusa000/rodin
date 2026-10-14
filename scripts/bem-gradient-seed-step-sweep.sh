#!/usr/bin/env bash
# Controlled native BEM seed/step experiment; each trial starts at the same seed.
set -u
if [[ $# != 2 || "$1" != /* || ! -x "$1" || "$2" != /* || -e "$2" ]]; then
  printf 'Usage: %s ABSOLUTE_EXECUTABLE NEW_ABSOLUTE_RESULTS_DIRECTORY\n' "$0" >&2
  exit 2
fi
executable=$1
results=$2
mkdir -p "$results" || exit 2
failed=0
for amplitude in 0.02 0.2 0.6; do
  directory="$results/seed-$amplitude"
  mkdir "$directory" || exit 2
  printf '%s\n' "$executable --seed-amplitude=$amplitude --step-factors=0,0.01,0.03,0.1" > "$directory/command.txt"
  (
    cd "$directory" || exit 2
    OMP_NUM_THREADS=1 "$executable" "--seed-amplitude=$amplitude" --step-factors=0,0.01,0.03,0.1 > run.log 2>&1
  )
  result=$?
  printf '%s\n' "$result" > "$directory/run.exit"
  if [[ $result != 0 ]]; then failed=1; fi
done
exit "$failed"
