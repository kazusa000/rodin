#!/usr/bin/env bash
# One chamber state/gradient validation and one chamber MMG update; never a long run.
set -u
if [[ $# != 2 || "$1" != /* || ! -x "$1" || "$2" != /* || -e "$2" ]]; then
  printf 'Usage: %s ABSOLUTE_KELVINBALL NEW_ABSOLUTE_RESULTS\n' "$0" >&2
  exit 2
fi
executable=$1
results=$2
validator="$(dirname "$executable")/KelvinBallBEMChamberValidation"
[[ -x "$validator" ]] || exit 2
mkdir -p "$results/validation" "$results/update" || exit 2
export OMP_NUM_THREADS=1
(cd "$results/validation" && "$validator" --reference > run.log 2>&1)
result=$?
printf '%s\n' "$result" > "$results/validation/run.exit"
if [[ $result != 0 ]]; then exit "$result"; fi
arguments=(--bem-gradient --n=13 --seed-amplitude=0.1 --iterations=2
  --step=0.10 --reconstruction=mmg --mmg-snap=0.001 --mmg-retries=0 --save-mesh)
printf '%q ' "$executable" "${arguments[@]}" > "$results/update/command.txt"
printf '\n' >> "$results/update/command.txt"
(cd "$results/update" && "$executable" "${arguments[@]}" > run.log 2>&1)
result=$?
printf '%s\n' "$result" > "$results/update/run.exit"
if [[ $result != 0 ]]; then exit "$result"; fi
# Read native C++ outputs; Python performs no state solve or direction construction.
python3 "$(dirname "$0")/bem-gradient-chamber-run.py" --check-update "$results"
result=$?
printf '%s\n' "$result" > "$results/validation.exit"
exit "$result"
