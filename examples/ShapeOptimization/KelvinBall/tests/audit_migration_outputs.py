#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Read-only migration audit; write only its requested external report."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

def read(path):
    return json.loads(Path(path).read_text())

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def check_bem(value, mesh, pin):
    assert value["source"]["commit"] == pin and value["source"]["status"] == ""
    assert sha(mesh) == value["mesh_sha256"]
    assert value["group_order"] == 24 and not value["subsampled"]
    assert value["full_panels"] == 24 * value["chamber_panels"]
    assert value["epsilon_factor"] == 0.6 and value["solver_maxiter"] == 2000
    rtol = value["solver_rtol"]
    assert rtol in (1e-6, 1e-9)
    d = value["diagnostics"]
    assert d["solver_converged"] and d["cpu_residual_verified"]
    assert d["solver_rtol"] == rtol
    assert all(0 <= x <= rtol for x in d["rhs_relative_residuals"] + d["cpu_rhs_relative_residuals"])
    assert d["cuda"]["dtype"] == "float64" and not d["cuda"]["fast_math"] and not d["cuda"]["fmad"]
    assert math.isclose(value["rho"], abs(value["c"]) / math.sqrt(value["k"]*value["q"]), rel_tol=1e-13)
    return dict(rtol=rtol, panels=value["full_panels"], rho=value["rho"],
                cpu_residuals=d["cpu_rhs_relative_residuals"],
                gpu_residuals=d["rhs_relative_residuals"])

def main():
    p=argparse.ArgumentParser()
    for key in ("old", "new", "replay", "old-replay", "checks", "pin", "output"):
        p.add_argument("--"+key, required=True)
    a=p.parse_args()
    old,new,replay,old_replay=map(Path,(a.old,a.new,a.replay,a.old_replay))
    campaign=read(new/"campaign.json")
    assert campaign["success"] and campaign["status"] == "both_backends_accepted"
    assert campaign["source"]["commit"] == a.pin and campaign["worker_commit"] == a.pin
    checks=read(Path(a.checks)/"run-record.json")
    assert checks["success"] and len(checks["cases"]) == 4
    assert checks["project"]["commit"] == a.pin
    assert all(c["success"] and c["exit_code"] == 0 and not c["timed_out"] for c in checks["cases"])
    measurements=[]
    for stage in ("initial", "mmg-trial-0", "swift-trial-0", "fd-0.01-1", "fd-0.01--1", "fd-0.005-1", "fd-0.005--1"):
        directory=stage+"-bem"
        prev,value=read(old/directory/"bem.json"),read(new/directory/"bem.json")
        audit=check_bem(value,new/directory/"input.mesh",a.pin)
        assert value["mesh_sha256"] == prev["mesh_sha256"], (stage,"candidate hash changed")
        assert value["table_sha256"] == prev["table_sha256"], (stage,"outer trace changed")
        assert all(value[k] == prev[k] for k in ("k","c","q","rho")), (stage,"coefficient changed")
        worker=read(new/directory/"worker.json")
        assert worker["success"] and worker["source"]["commit"] == a.pin
        measurements.append(dict(stage=stage,mesh_sha256=value["mesh_sha256"],coefficients_bitwise_equal=True,outer_trace_bitwise_equal=True,**audit))
    steps=[]
    old_campaign=read(old/"campaign.json")
    for trial in campaign["trials"]:
        assert trial["trial"] == 0 and trial["accepted"]
        label=trial["method"]+"-trial-0"
        previous=next(x for x in old_campaign["trials"] if x["method"] == trial["method"] and x["trial"] == 0)
        for key in ("direction","thickness","geometry","accepted","acceptance","predicted_delta","actual_delta"):
            assert trial[key] == previous[key], (label,key)
        assert trial["geometry"]["solid_connectivity"]["full_components"] == 1
        log=(new/label/"run.log").read_text()
        residuals=[float(x) for x in re.findall(r"Chamber Stokes loadOffset=\d+ residual=([^ ]+)",log)]
        assert len(residuals) == 2 and max(residuals) <= 1e-8
        steps.append(dict(method=trial["method"],rho=trial["bem"]["rho"],penalty=trial["thickness"]["penalty"],
                          full_panels=trial["bem"]["full_panels"],original_stokes_residuals=residuals,
                          acceptance=trial["acceptance"],same_as_before=True))
    assert campaign["finite_differences"] == old_campaign["finite_differences"]
    assert all(x["success"] for x in campaign["same_target"])
    for x in campaign["same_target"]:
        assert x["run"]["project"]["commit"] == a.pin
    replay_checks=[]
    for method in ("mmg","swift"):
        current=read(replay/method/"audit.json")
        previous=read(old_replay/method/"audit.json")
        assert current == previous,(method,"replay audit changed")
        run=read(replay/method/"run-record.json")
        assert run["success"] and run["project"]["commit"] == a.pin
        replay_checks.append(dict(method=method,candidate_sha256=current["candidate_sha256"],
                                  quality_min=current["quality_min"],same_as_before=True))
    report=dict(passed=True,executing_commit=a.pin,native_checks=4,measurements=measurements,
                accepted_steps=steps,finite_differences_unchanged=True,
                native_proxy_is_not_exact_bem_gradient=True,saved_stopping_replay=replay_checks)
    output=Path(a.output)
    assert not output.exists()
    output.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))

if __name__ == "__main__":
    main()
