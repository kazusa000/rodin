"""Recover original tetrahedra and chamber face labels from saved XDMF doubles."""
from pathlib import Path
import hashlib
import json

import h5py
import numpy as np


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def export_restart(snapshot, output, outer_radius):
    """Write a lossless MEDIT volume checkpoint, rejecting ambiguous boundaries."""
    snapshot, output = Path(snapshot), Path(output)
    if output.exists():
        raise ValueError('Require a fresh restart export')
    source_hash = digest(snapshot)
    with h5py.File(snapshot) as data:
        vertices = data['Mesh/Geometry/Vertices'][:]
        stream = data['Mesh/XDMF/Topology'][:]
        attributes = [data[f'Mesh/Attributes/{dim}'][:] for dim in range(4)]
    tets = stream.reshape(-1, 5)
    if not np.all(tets[:, 0] == 6) or not np.isfinite(vertices).all():
        raise ValueError('Require finite original linear tetrahedra')
    tets = tets[:, 1:].astype(np.int64)
    labels = attributes[3]
    if len(tets) != len(labels) or set(np.unique(labels)) != {2, 3}:
        raise ValueError('Require original solid/fluid materials')
    # The paired x=y labels can move away from z=0 during upstream MMG.
    # Recover their original mapping from the last preserved native fluid
    # coverage mesh, never replace it by a new sign-based partition.
    native = snapshot.parent / "kelvin-cut-coverage-failure.mesh"
    fluid_snapshot = snapshot.with_name(snapshot.name.replace("Chamber", "Fluid"))
    native_hash = digest(native)
    lines = native.read_text().splitlines()
    def section(name):
        index = next(i for i, line in enumerate(lines) if line.strip() == name) + 1
        while not lines[index].strip():
            index += 1
        count = int(lines[index])
        return lines[index+1:index+1+count]
    ascii_vertices = np.array([[float(x) for x in line.split()[:3]] for line in section("Vertices")])
    with h5py.File(fluid_snapshot) as data:
        fluid_vertices = data['Mesh/Geometry/Vertices'][:]
        fluid_stream = data['Mesh/XDMF/Topology'][:].reshape(-1, 5)
        fluid_labels = data['Mesh/Attributes/3'][:]
    if ascii_vertices.shape != fluid_vertices.shape or not np.allclose(ascii_vertices, fluid_vertices, rtol=6e-6, atol=1e-9):
        raise ValueError('The native coverage mesh is not the matching saved fluid state')
    ascii_tets = np.array([[int(x) for x in line.split()] for line in section("Tetrahedra")])
    if not np.array_equal(ascii_tets[:, :4] - 1, fluid_stream[:, 1:]) or not np.array_equal(ascii_tets[:, 4], fluid_labels):
        raise ValueError('The native coverage mesh has different fluid tetrahedra or materials')
    lookup = {tuple(point): i for i, point in enumerate(vertices)}
    if len(lookup) != len(vertices):
        raise ValueError('Coincident original vertices require an explicit parent map')
    try:
        fluid_to_full = np.array([lookup[tuple(point)] for point in fluid_vertices])
    except KeyError as error:
        raise ValueError('A preserved fluid coordinate is absent from the complete volume') from error
    native_faces = {}
    for line in section("Triangles"):
        a, b, c, label = map(int, line.split())
        if label in (7, 13, 31, 32, 33, 35):
            native_faces[tuple(sorted(fluid_to_full[np.array([a, b, c])-1]))] = label
    templates = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]])
    faces = np.sort(tets[:, templates].reshape(-1, 3), axis=1)
    unique, inverse, counts = np.unique(faces, axis=0, return_inverse=True, return_counts=True)
    if counts.max() > 2:
        raise ValueError('Nonmanifold tetrahedral faces')
    solid_owners = np.bincount(inverse, weights=np.repeat(labels == 2, 4))
    gamma = (counts == 2) & (solid_owners == 1)
    selected = (counts == 1) | gamma
    triangles = unique[selected]
    triangle_labels = np.zeros(len(triangles), dtype=np.int64)
    triangle_labels[gamma[selected]] = 13
    # The upstream projection fixes these planes. This tolerance classifies
    # coordinates only; it never moves or merges any point.
    plane_tolerance = 1e-10
    for i, face in enumerate(triangles):
        if triangle_labels[i]:
            continue
        points = vertices[face]
        preserved = native_faces.get(tuple(face))
        if preserved is not None:
            triangle_labels[i] = preserved
            continue
        if np.max(np.abs(points[:, 0] - outer_radius)) <= plane_tolerance:
            triangle_labels[i] = 7
        elif np.max(np.abs(points[:, 0] - points[:, 1])) <= plane_tolerance:
            references = set(int(v) for v in attributes[0][face]) & {33, 35}
            if len(references) == 1:
                triangle_labels[i] = next(iter(references))
            elif len(references) == 0 and np.all(points[:, 2] >= 0):
                triangle_labels[i] = 33
            elif len(references) == 0 and np.all(points[:, 2] <= 0):
                triangle_labels[i] = 35
            else:
                raise ValueError(f'A solid-cut triangle has ambiguous original references: {face}, {references}')
        elif np.max(np.abs(points[:, 1] - points[:, 2])) <= plane_tolerance:
            triangle_labels[i] = 31
        elif np.max(np.abs(points[:, 1] + points[:, 2])) <= plane_tolerance:
            triangle_labels[i] = 32
        else:
            raise ValueError(f'Unclassified original boundary triangle {face}')
    null = np.iinfo(np.uint64).max
    expected_values, expected_counts = np.unique(attributes[2], return_counts=True)
    expected = {int(v): int(n) for v, n in zip(expected_values, expected_counts) if v != null}
    values, nums = np.unique(triangle_labels, return_counts=True)
    actual = {int(v): int(n) for v, n in zip(values, nums)}
    if actual != expected or len(unique) != len(attributes[2]):
        raise ValueError(f'Recovered face labels differ from original: {actual} != {expected}')
    edge_labels = {}
    for face, label in zip(triangles, triangle_labels):
        for edge in ((face[0], face[1]), (face[0], face[2]), (face[1], face[2])):
            edge_labels.setdefault(edge, set()).add(int(label))
    fixed = {7, 31, 32, 33, 35}
    ridges = [edge for edge, refs in edge_labels.items() if len(refs) >= 2 and refs & fixed]
    expected_edges = {int(v): int(n) for v, n in zip(*np.unique(attributes[1], return_counts=True)) if v != null}
    if expected_edges != {40: len(ridges)}:
        raise ValueError(f'Ridge labels differ from original: {expected_edges} vs {len(ridges)}')
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('w') as file:
        file.write('MeshVersionFormatted 2\nDimension 3\nVertices\n'+str(len(vertices))+'\n')
        for point, label in zip(vertices, attributes[0]):
            file.write(' '.join(format(float(x), '.17g') for x in point)+f' {int(label)}\n')
        file.write('Edges\n'+str(len(ridges))+'\n')
        for a, b in ridges:
            file.write(f'{a+1} {b+1} 40\n')
        file.write('Triangles\n'+str(len(triangles))+'\n')
        for face, label in zip(triangles, triangle_labels):
            file.write(' '.join(str(int(v)+1) for v in face)+f' {label}\n')
        file.write('Tetrahedra\n'+str(len(tets))+'\n')
        for tet, label in zip(tets, labels):
            file.write(' '.join(str(int(v)+1) for v in tet)+f' {label}\n')
        file.write('End\n')
    # Verify decimal serialization restores the exact float64 values.
    with output.open() as file:
        for _ in range(4):
            next(file)
        restored = np.array([[float(x) for x in next(file).split()[:3]] for _ in vertices])
    if not np.array_equal(restored, vertices) or digest(snapshot) != source_hash:
        raise ValueError('Original double coordinates or source changed')
    record = dict(source=str(snapshot), source_sha256=source_hash, output=str(output),
                  output_sha256=digest(output), vertices=len(vertices), tetrahedra=len(tets),
                  gamma_panels=int(gamma.sum()), face_label_counts=actual, ridge_edges=len(ridges),
                  original_tet_order_and_materials=True, exact_double_roundtrip=True,
                  geometry_modified=False, boundary_label_method='matching original native fluid faces plus original solid-cut vertex references; original label counts required',
                  native_fluid=str(native), native_fluid_sha256=native_hash, solid_cut_fallback='whole-triangle hemisphere only when no original vertex cut reference exists; full direction regression required', fluid_snapshot=str(fluid_snapshot),
                  fluid_snapshot_sha256=digest(fluid_snapshot),
                  plane_classification_tolerance=plane_tolerance)
    output.with_suffix('.json').write_text(json.dumps(record, indent=2)+'\n')
    return record
