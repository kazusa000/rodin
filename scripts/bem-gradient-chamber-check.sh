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
python3 - "$results" <<'CHECK'
import csv,json,math,pathlib,sys
root=pathlib.Path(sys.argv[1]);u=root/'update'
states=[json.loads((u/f'bem-state-{i}.json').read_text()) for i in [0,1]]
rows=list(csv.DictReader((u/'kelvin-ball.csv').open()))
assert len(rows)==2,'Expected exactly two chamber evaluations'
old,new=[s['rho'] for s in states]
predicted=float(rows[1]['predicted_delta_rho'])
checks={'strict_rho_increase':new>old,'weak_armijo':predicted>0 and new>=old+1e-4*predicted,
 'closed_rotated_interface':all(s['open_edges']==0 for s in states),
 'native_state_residual':all(math.isfinite(s['residual']) and s['residual']<=1e-9 for s in states),
 'chamber_state_space':all(s['design_domain']=='O24-chamber' and s['independent_dofs']<=s['state_dofs'] for s in states),
 'volume_projection':abs(float(rows[0]['d_volume_theta']))<=1e-10,
 'real_transport':float(rows[0]['advection_increment'])>1e-12,
 'chamber_checkpoints':all((u/f'bem-chamber-{i}.mesh').exists() for i in [0,1])}
record={'passed':all(checks.values()),'design_domain':'O24-chamber','gradient_domain':'O24-chamber',
 'reconstruction_domain':'O24-chamber','bem_domain':'infinite-exterior','checks':checks,
 'initial_rho':old,'updated_rho':new,'predicted_delta_rho':predicted,'state_evaluations':2,'real_updates':1,
 'volume_projection_enabled':True,'volume_rejection_enabled':False,'volume_tolerance':None,
 'relative_chamber_volume_change':float(rows[1]['volume'])/float(rows[0]['volume'])-1}
(root/'chamber-update-validation.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
raise SystemExit(0 if record['passed'] else 2)
CHECK
result=$?
printf '%s\n' "$result" > "$results/validation.exit"
exit "$result"
