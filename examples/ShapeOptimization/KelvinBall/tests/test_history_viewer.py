import base64
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
import gzip
import json
from pathlib import Path
import re
import struct

import pytest

ROOT = Path(__file__).resolve().parents[1]
from kelvin_bem import history_viewer as viewer


def score():
    return dict(mesh_sha256="original", group_order=12, subsampled=False,
                full_panels=12, chamber_panels=1, k=1., c=-.002, q=4., rho=.001,
                diagnostics=dict(solver_converged=True, cpu_residual_verified=True,
                                 rhs_relative_residuals=[9e-11, 8e-11],
                                 cpu_rhs_relative_residuals=[9e-11, 8e-11]))


def test_rejects_failed_measurement_and_wrong_input():
    with pytest.raises(ValueError, match="input hash"):
        viewer.validate_score(score(), "different")
    s = score()
    s["diagnostics"]["cpu_rhs_relative_residuals"][1] = 5e-9
    with pytest.raises(ValueError, match="residual acceptance"):
        viewer.validate_score(s, "original")


def test_rejects_cubic_subsampling_and_wrong_coefficients():
    for key, value in [("group_order", 24), ("subsampled", True), ("full_panels", 6)]:
        s = score()
        s[key] = value
        with pytest.raises(ValueError, match="all original T12"):
            viewer.validate_score(s, "original")
    s = score()
    s["rho"] = .02
    with pytest.raises(ValueError, match="k/c/q"):
        viewer.validate_score(s, "original")


def test_lossless_payload_retains_mesh_and_measured_values(tmp_path):
    f = dict(step=13, points=[.12345678912345678, .25, -.5], triangles=[0, 0, 0],
             rho=.001168267727287026)
    data = dict(frames=[f], rotations=[[[1, 0, 0], [0, 1, 0], [0, 0, 1]]])
    output = tmp_path / "viewer.html"
    record = viewer.render(data, output)
    encoded = re.search(r'application/octet-stream">([^<]+)</script>', output.read_text())[1]
    raw = gzip.decompress(base64.b64decode(encoded))
    assert raw[:4] == b"MSH3"
    length = struct.unpack_from("<I", raw, 4)[0]
    decoded = json.loads(raw[8:8+length])
    frame = decoded["frames"][0]
    offset = 8 + length
    count = frame.pop("pointValues")
    packed = bytes(raw[offset+byte*count+i] for i in range(count) for byte in range(8))
    frame["points"] = list(struct.unpack(f"<{count}d", packed))
    offset += count * 8
    count = frame.pop("triangleValues")
    width = frame.pop("indexWidth")
    packed = bytes(raw[offset+byte*count+i] for i in range(count) for byte in range(width))
    frame["triangles"] = list(struct.unpack(f"<{count}{'H' if width == 2 else 'I'}", packed))
    previous = 0
    for i in range(0, count, 3):
        first = previous + frame["triangles"][i]
        frame["triangles"][i] = first
        frame["triangles"][i+1] += first
        frame["triangles"][i+2] += first
        previous = first
    assert offset + count * width == len(raw)
    assert decoded == data
    assert record["geometry_modified"] is False
    assert record["bytes"] < 1_000_000
    with pytest.raises(FileExistsError):
        viewer.render(data, output)


@pytest.mark.parametrize("rtol", [1e-9, 1e-6])
def test_effective_tolerance_certificate(rtol):
    s = score()
    s["solver_rtol"] = s["diagnostics"]["solver_rtol"] = rtol
    s["diagnostics"]["rhs_relative_residuals"] = [rtol * .9] * 2
    s["diagnostics"]["cpu_rhs_relative_residuals"] = [rtol * .91] * 2
    viewer.validate_score(s, "original")
    s["diagnostics"]["cpu_rhs_relative_residuals"][0] = rtol * 1.01
    with pytest.raises(ValueError, match="residual acceptance"):
        viewer.validate_score(s, "original")


def test_rejects_inconsistent_or_unknown_tolerance():
    for rtol, diagnostic in [(1e-5, 1e-5), (1e-6, 1e-9)]:
        s = score()
        s["solver_rtol"] = rtol
        s["diagnostics"]["solver_rtol"] = diagnostic
        with pytest.raises(ValueError, match="tolerance certificate"):
            viewer.validate_score(s, "original")


def campaign(start, updates):
    history = [dict(update=i, mesh=f"mesh-{start+i}",
                    bem=dict(mesh_sha256=f"hash-{start+i}"), geometry={})
               for i in range(updates + 1)]
    return dict(symmetry="terra", history=history, accepted_updates=updates,
                final=history[-1]["bem"], final_mesh=history[-1]["mesh"],
                trials=[dict(update=updates+1, accepted=False)])


def test_current_chain_contains_exactly_accepted_steps_0_to_22():
    campaigns = [campaign(0, 9), campaign(9, 11), campaign(20, 2)]
    entries = viewer.chain_entries(campaigns)
    assert [e[1] for e in entries] == list(range(23))
    assert len({e[3]["mesh_sha256"] for e in entries}) == 23
    assert entries[0][0] == "初始扰动球"
    campaigns[1]["history"][0]["bem"]["mesh_sha256"] = "wrong"
    with pytest.raises(ValueError, match="seed differs"):
        viewer.chain_entries(campaigns)


def test_rejects_history_with_missing_checkpoint():
    c = campaign(0, 9)
    del c["history"][4]
    with pytest.raises(ValueError, match="checkpoint is inconsistent"):
        viewer.chain_entries([c])


def test_cubic_uses_full24_and_never_relables_a_t12_measurement():
    s=score()
    with pytest.raises(ValueError,match='all original O24'):
        viewer.validate_score(s,'original','cubic')
    s['group_order']=s['full_panels']=24
    viewer.validate_score(s,'original','cubic')
    with pytest.raises(ValueError,match='all original T12'):
        viewer.validate_score(s,'original')
    campaigns=[campaign(0,1),campaign(1,16)]
    for c in campaigns: c['symmetry']='cubic'
    entries=viewer.chain_entries(campaigns,'cubic')
    assert [e[1] for e in entries]==list(range(18))
    assert len({e[3]['mesh_sha256'] for e in entries})==18
    with pytest.raises(ValueError,match='selected symmetry'):
        viewer.chain_entries(campaigns)


@pytest.mark.parametrize('symmetry,group',[('terra',12),('cubic',24)])
def test_frames_use_selected_rotation_count_and_volume(tmp_path,symmetry,group):
    from types import SimpleNamespace
    import hashlib
    import numpy as np
    mesh=tmp_path/'accepted.mesh';mesh.write_text('raw original interface')
    digest=hashlib.sha256(mesh.read_bytes()).hexdigest()
    s=dict(score(),mesh_sha256=digest,group_order=group,full_panels=group)
    points=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
    faces=np.array([[0,1,2]])
    boundary=SimpleNamespace(read_mesh=lambda path:(points,faces,np.array([13]),[],[]))
    rotations12=['T']*12;rotations24=['O']*24
    exporter=SimpleNamespace(material_interface=lambda *args:(faces,1.),geometry=SimpleNamespace(
        proper_tetrahedral_rotations=lambda:rotations12,proper_octahedral_rotations=lambda:rotations24))
    geometry=dict(mesh_sha256=digest,volume_full=group,**{f'full{group}':dict(boundary_edges=9)})
    data=viewer.collect_frames([('shape',0,mesh,s,geometry)],[dict(status='done',source=dict(commit='frozen'))],
        tmp_path,boundary,exporter,symmetry)
    assert data['rotations']==(rotations12 if group==12 else rotations24)
    assert data['frames'][0]['volume']==group and data['frames'][0]['fullTriangles']==group


def test_compression_reindexing_retains_exact_oriented_geometry():
    points = [1.234567891234567, 0., 0., 0., 1., 0., 0., 0., 1., -1., 0., 0.]
    faces = [0, 2, 3, 0, 1, 2, 3, 2, 1]
    reordered_points, reordered_faces = viewer.locality_order(points, faces)
    def geometry(p, triangles):
        oriented = []
        for i in range(0, len(triangles), 3):
            face = tuple(tuple(p[3*v:3*v+3]) for v in triangles[i:i+3])
            oriented.append(min(face, face[1:]+face[:1], face[2:]+face[:2]))
        return sorted(oriented)
    assert geometry(points, faces) == geometry(reordered_points, reordered_faces)
