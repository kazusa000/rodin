from pathlib import Path
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
import pytest
from kelvin_bem import compare_backends as pair


def data():
    b0, b1 = dict(rho=.001), dict(rho=.0011)
    t = dict(group_order=24, passed=True, minimum=.06, penalty=1e-5, deepest=.01)
    g = dict(solid_connectivity=dict(group_order=24, full_components=1))
    d = dict(rho=.001, penalty=1e-5, dt=.01, merit_slope=.002)
    return b0, b1, t, t, g, g, d


def test_actual_bem_merit_accepts_not_native_prediction():
    values = data()
    assert pair.decide(*values)[0]
    values[1]["rho"] = .0009
    assert not pair.decide(*values)[0]


def test_fragments_rejected_even_with_large_rho():
    values = list(data())
    values[5] = dict(solid_connectivity=dict(group_order=24, full_components=2))
    assert not pair.decide(*values)[0]


def test_wrong_native_coefficients_fail():
    values = data()
    values[-1]["rho"] *= 2
    with pytest.raises(ValueError, match="consume"):
        pair.decide(*values)


def test_safe_shell_quoting():
    assert pair.shell(["/a b", "x;y"]) == "'/a b' 'x;y'"


def test_only_identified_candidate_failures_are_backtracked():
    assert pair.known_rejection("native",-6,False,
        "Mesh reconstruction changed the fixed boundary geometry.")
    assert pair.known_rejection("bem",1,False,
        "RuntimeError: quotient-CG did not converge in 2000 iterations;")
    assert pair.known_rejection("bem",1,False,"wrong BEM verifier module") is None
    assert pair.known_rejection("native",1,False,"source must be clean and pinned") is None
    assert pair.known_rejection("bem",1,True,"quotient-CG did not converge in 2000 iterations;") is None


def test_native_failure_keeps_actual_record_and_log(tmp_path):
    import argparse, json
    c = object.__new__(pair.Controller)
    c.args = argparse.Namespace(cpu_host="cpu",cpu_python="/python",
        cpu_source="/source",build="/build")
    c.cpu, c.output = Path("/remote"), tmp_path
    def ssh(*args):
        raise pair.CommandFailure(1,False,tmp_path/"launch.log")
    c.ssh = ssh
    def get(host,remote,local):
        local.parent.mkdir(parents=True,exist_ok=True)
        if str(remote).endswith("run-record.json"):
            local.write_text(json.dumps(dict(success=False,cases=[
                dict(exit_code=-6,timed_out=False)])))
        else:
            local.write_text("Mesh reconstruction changed the fixed boundary geometry.")
        return local
    c.get = get
    with pytest.raises(pair.CandidateFailure) as caught:
        c.native("candidate","step")
    assert caught.value.code == -6
    assert caught.value.rejection() == "invalid fixed-boundary geometry"


@pytest.mark.parametrize("online", [True, False])
def test_node_readiness_uses_configured_deployment_not_personal_ips(monkeypatch, online):
    import argparse, json
    controller = object.__new__(pair.Controller)
    controller.args = argparse.Namespace(
        tailscale="tailscale", cpu_tailscale_ip="100.1.2.3", gpu_tailscale_ip="100.4.5.6",
        cpu_host="other-cpu", gpu_host="other-gpu", cpu_source="/cpu/rodin",
        gpu_source="/gpu/rodin", worker_commit="same-sha")
    calls = []
    def output(command, **kwargs):
        calls.append(command)
        if command[0] == "tailscale":
            return json.dumps(dict(BackendState="Running", Self=dict(Online=True), Peer={
                "cpu": dict(TailscaleIPs=["100.1.2.3"], Online=True),
                "gpu": dict(TailscaleIPs=["100.4.5.6"], Online=online)}))
        assert command[0] == "ssh" and command[3] in ("other-cpu", "other-gpu")
        return "same-sha" if "rev-parse" in command[-1] else ""
    monkeypatch.setattr(pair.subprocess, "check_output", output)
    if online:
        controller.prepare()
        assert len(calls) == 5
    else:
        with pytest.raises(RuntimeError, match="Tailscale-online"):
            controller.prepare()
        assert len(calls) == 1
