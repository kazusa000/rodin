#!/usr/bin/env python3
"""Resume the original Cubic global10 with baseline thickness, then independent BEM."""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time

from run_fem_old import EXAMPLE, EXPORTER, VERIFIER, complete_rows, digest, git, pinned, save

BASELINE = '3ae2c0414ffebf11f5814a46510ddf4d8375ba0a'
OLD = 'fed5689233fc035e99f9feb1fe4bad8d9e681eef'
CPU = 'fa595d6f9f7f002b23127d202a9608f4b009d09c'
LIBRARY = 'cb4260eff5d79bbc276c623d95f0119d755042de8b0e79e691959cfcdb527a66'
INPUT = '5f2e8ab6bb31180c7447a55625d2c9abee55f3d96c12045739301bb9a0ff30d2'
PARAMETERS = ['--n=20', '--outer-radius=2', '--regularization=100', '--step=0.5',
    '--advection-quadrature=3', '--reconstruction=mmg', '--mmg-adapt',
    '--mmg-snap=0.005', '--thickness-weight=1', '--save-mesh']


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def restart_mesh(snapshot, target):
    """Recover labelled cut/outer faces from the unchanged tetrahedral complex.

    XDMF snapshots omit lower-dimensional connectivity, so their entity numbers
    cannot be restored. Face geometry/material adjacency uniquely recover the
    original boundary conditions; histograms must match every stored label.
    """
    import h5py
    import numpy as np
    from collections import Counter, defaultdict
    with h5py.File(snapshot) as f:
        vertices = f['Mesh/Geometry/Vertices'][:]
        tets = f['Mesh/XDMF/Topology'][:].reshape(-1, 5)[:, 1:]
        attrs = [f['Mesh/Attributes/'+str(i)][:] for i in range(4)]
    null = np.iinfo(np.uint64).max
    owners = defaultdict(list)
    for i, tet in enumerate(tets):
        for local in [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]:
            owners[tuple(sorted(int(tet[k]) for k in local))].append(i)
    faces = []
    histogram = Counter()
    edge_labels = defaultdict(set)
    vertex_labels = defaultdict(set)
    for face, cells in owners.items():
        if len(cells) == 2:
            label = 13 if attrs[3][cells[0]] != attrs[3][cells[1]] else int(null)
        elif len(cells) == 1:
            x = vertices[list(face)]
            if np.max(np.abs(x[:, 0]-2)) <= 1e-8:
                label = 7
            elif np.max(np.abs(x[:, 1]-x[:, 2])) <= 1e-8:
                label = 31
            elif np.max(np.abs(x[:, 1]+x[:, 2])) <= 1e-8:
                label = 32
            elif np.max(np.abs(x[:, 0]-x[:, 1])) <= 1e-8:
                label = 33 if x[:, 2].mean() >= 0 else 35
            else:
                raise ValueError('Unclassified original chamber boundary')
        else:
            raise ValueError('Nonmanifold volume face')
        histogram[label] += 1
        if label != null:
            faces.append((*face, label))
            for v in face:
                vertex_labels[v].add(label)
            for pair in [(0, 1), (1, 2), (0, 2)]:
                edge_labels[tuple(sorted(face[k] for k in pair))].add(label)
    expected = Counter(map(int, attrs[2]))
    if histogram != expected:
        raise ValueError('Recovered face labels differ from original HDF5 counts')
    fixed = {7, 31, 32, 33, 35}
    features = [(edge, labels) for edge, labels in edge_labels.items()
                if len(labels) >= 2 and labels & fixed]
    if sum(attrs[1] == 40) != len(features) or set(map(int, attrs[1])) != {40, int(null)}:
        raise ValueError('Recovered feature edges differ from original counts')
    ridges = [i+1 for i, (_, labels) in enumerate(features)
              if len({33 if label == 35 else label for label in labels}) >= 2]
    corners = [v+1 for v, labels in vertex_labels.items() if len(labels) >= 3 and labels & fixed]
    with target.open('w') as stream:
        stream.write('MeshVersionFormatted 2\nDimension 3\nVertices\n'+str(len(vertices))+'\n')
        for x, ref in zip(vertices, attrs[0]):
            if ref == null:
                raise ValueError('Unexpected unset source vertex reference')
            stream.write(' '.join(format(float(c), '.17g') for c in x)+' '+str(ref)+'\n')
        stream.write('Edges\n'+str(len(features))+'\n')
        for (a, b), _ in features:
            stream.write(f'{a+1} {b+1} 40\n')
        stream.write('Triangles\n'+str(len(faces))+'\n')
        for a, b, c, label in faces:
            stream.write(f'{a+1} {b+1} {c+1} {label}\n')
        stream.write('Tetrahedra\n'+str(len(tets))+'\n')
        for tet, label in zip(tets, attrs[3]):
            stream.write(' '.join(str(int(v)+1) for v in tet)+' '+str(label)+'\n')
        for name, entities in [('Ridges', ridges), ('Corners', corners)]:
            stream.write(name+'\n'+str(len(entities))+'\n')
            for entity in entities:
                stream.write(str(entity)+'\n')
        stream.write('End\n')
    return dict(geometry_modified=False, tetrahedra_unchanged=True,
                face_label_counts=dict(histogram), feature_edges=len(features),
                ridges=len(ridges), corners=len(corners), lower_dimensional_entities_reindexed=True)


def same_geometry(original, loaded):
    import h5py
    import numpy as np
    from collections import Counter
    with h5py.File(original) as a, h5py.File(loaded) as b:
        names = ['Geometry/Vertices', 'XDMF/Topology', 'Attributes/0', 'Attributes/3']
        checks = {n: np.array_equal(a['Mesh/'+n][:], b['Mesh/'+n][:]) for n in names}
        for d in [1, 2]:
            checks['Attributes/'+str(d)+'_counts'] = (Counter(map(int, a['Mesh/Attributes/'+str(d)][:]))
                == Counter(map(int, b['Mesh/Attributes/'+str(d)][:])))
    if not all(checks.values()):
        raise ValueError('Restart altered coordinates/cells/labels: '+str(checks))
    return checks


def check_direction(row):
    for key in ['thickness_penalty', 'thickness_nominal_penalty', 'thickness_minimum_exit',
                'thickness_maximum_deficit', 'thickness_repair_d_volume', 'd_thickness_theta',
                'd_volume_theta', 'thickness_guard_distance']:
        if not math.isfinite(float(row[key])):
            raise ValueError('Missing finite thickness diagnostic: '+key)
    if abs(float(row['d_volume_theta'])) > 1e-9 or abs(float(row['thickness_repair_d_volume'])) > 1e-9:
        raise ValueError('Independent volume-nullspace projection failed')
    if float(row['d_thickness_theta']) > 1e-10:
        raise ValueError('Selected direction worsens thickness surrogate')
    if not math.isclose(float(row['thickness_guard_distance']), 3*2/19, rel_tol=1e-12):
        raise ValueError('Wrong guard distance')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ['source', 'prefix', 'build', 'output', 'input', 'exporter', 'verifier', 'cpu_source', 'library']:
        parser.add_argument('--'+key.replace('_', '-'), type=Path, required=True)
    parser.add_argument('--commit', required=True)
    parser.add_argument('--reuse-build', action='store_true',
                        help='Reuse an inactive build configured for this exact source path')
    args = parser.parse_args()
    if platform.node() != '28h':
        raise ValueError('All computation is authorized only on 28h')
    for source, commit in [(args.source, args.commit), (args.exporter, EXPORTER),
                           (args.verifier, VERIFIER), (args.cpu_source, CPU)]:
        pinned(source, commit)
    if digest(args.library) != LIBRARY:
        raise ValueError('Unverified native BEM library')
    unchanged = ['src', EXAMPLE+'Sphere.cpp', EXAMPLE+'Configuration.cpp', EXAMPLE+'Common.cpp',
                 EXAMPLE+'Metrics.cpp', EXAMPLE+'RotatedNitscheIntegrator.cpp', EXAMPLE+'SewedOutput.cpp']
    if git(args.source, 'diff', OLD, args.commit, '--', *unchanged):
        raise ValueError('Old numerical core changed outside restart/thickness')
    manifest = json.loads((args.input/'campaign.json').read_text())
    if manifest['commit'] != OLD or manifest['symmetry'] != 'cubic':
        raise ValueError('Wrong original experiment')
    old_rows = complete_rows(args.input/'optimization/kelvin-ball.csv')
    original = old_rows[10]
    target = float(old_rows[0]['target_volume'])
    snapshot = args.input/'optimization/KelvinBall.Chamber.mesh.000010.h5'
    if digest(snapshot) != INPUT or not math.isclose(target, .17391249230727751, rel_tol=1e-14):
        raise ValueError('Original global10 or target volume changed')
    if args.output.exists() or args.output.resolve().is_relative_to(args.source.resolve()):
        raise ValueError('Require fresh output outside source')
    if args.build.resolve().is_relative_to(args.source.resolve()):
        raise ValueError('Build must be outside source')
    if args.build.exists():
        expected = 'CMAKE_HOME_DIRECTORY:INTERNAL='+str(args.source/EXAMPLE/'reproduction')
        if not args.reuse_build or expected not in (args.build/'CMakeCache.txt').read_text().splitlines():
            raise ValueError('Existing build is not the explicitly requested source-matched cache')
    args.output.mkdir(parents=True)
    record = dict(state='preparing', commit=args.commit, thickness_baseline_commit=BASELINE,
        original_commit=OLD, input_snapshot=str(snapshot), input_sha256=INPUT,
        original_target_volume=target, local0_global_step=10, requested_real_updates=5,
        requested_solved_states=6, final_global_step=15, bem_feedback=False,
        nominal_thickness=2*2/19, guard_distance=3*2/19, started_at=time.time(), stages=[])
    environment = dict(os.environ, KELVIN_SYMMETRY='cubic', OMP_NUM_THREADS='6',
        OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    lib = args.prefix/'lib/x86_64-linux-gnu'
    environment['LD_LIBRARY_PATH'] = ':'.join(map(str, [lib, lib/'openmpi/lib', lib/'hdf5/serial']))
    record['reused_build'] = args.reuse_build
    record['environment'] = {k: environment[k] for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'LD_LIBRARY_PATH']}
    proc = None
    checkpoint = lambda: save(args.output/'campaign.json', record)

    def stop(signum, _frame):
        if proc is not None and proc.poll() is None:
            os.killpg(proc.pid, signal.SIGTERM)
        raise InterruptedError('Stopped by signal '+str(signum))
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)

    def run(name, command, hours):
        nonlocal proc
        directory = args.output/name
        directory.mkdir(parents=True)
        item = dict(name=name, command=list(map(str, command)), started_at=time.time())
        record['stages'].append(item)
        with (directory/'run.log').open('w') as log:
            proc = subprocess.Popen(command, cwd=directory, env=environment,
                stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            item['pid'] = proc.pid
            checkpoint()
            try:
                actual = proc.wait(timeout=hours*3600)
            except (subprocess.TimeoutExpired, InterruptedError) as error:
                item['interrupted'] = repr(error)
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    actual = proc.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    actual = proc.wait()
        item.update(exit_code=actual, finished_at=time.time())
        (directory/'exit-code.txt').write_text(str(actual)+'\n')
        checkpoint()
        if item.get('interrupted'):
            raise RuntimeError('Stage interrupted: '+name)
        return actual

    code = 1
    try:
        sys.path.insert(0, str(args.exporter/'scripts'))
        os.environ['KELVIN_SYMMETRY'] = 'cubic'
        export = module('restart_export', args.exporter/'scripts/reproduce_smallbox_original.py')
        directory = args.output/'restart'
        directory.mkdir()
        native = directory/'native-export.mesh'
        mesh = directory/'original-double.mesh'
        export.export_snapshot(snapshot, native)
        record['restart_certificate'] = export.original_double_mesh(snapshot, native, mesh)
        if record['restart_certificate']['gamma_panels'] != int(original['interface_triangles']):
            raise ValueError('Restart lost original panels')
        mesh = directory/'restart-labelled.mesh'
        record['restart_label_recovery'] = restart_mesh(snapshot, mesh)
        record['restart_mesh_sha256'] = digest(mesh)
        configure = ['cmake', '-S', str(args.source/EXAMPLE/'reproduction'), '-B', str(args.build),
            '-DCMAKE_BUILD_TYPE=Release', f'-DCMAKE_PREFIX_PATH={args.prefix}', f'-DMUMPS_ROOT={args.prefix}',
            f'-DMPI_C_HEADER_DIR={lib/"openmpi/include"}',
            f'-DMPI_C_COMPILER_INCLUDE_DIRS={lib/"openmpi/include"};{lib/"openmpi/include/openmpi"}',
            '-DMPI_C_LIB_NAMES=mpi', f'-DMPI_mpi_LIBRARY={lib/"openmpi/lib/libmpi.so"}',
            f'-DCARLOS_MUMPS_SUPPORT_LIBRARIES={lib/"openmpi/lib/libmpi_cxx.so"}']
        if run('configure', configure, 1) or run('compile', ['cmake', '--build', str(args.build),
                '--target', 'KelvinBall', 'KelvinBallCubicGeometryTest', 'KelvinBallTetraGeometryTest',
                'KelvinBallThicknessTest', '-j', '2'], 4):
            raise RuntimeError('Build failed')
        binary = args.build/'bin/KelvinBall'
        record['binary_sha256'] = digest(binary)
        for name in ['CubicGeometry', 'TetraGeometry', 'Thickness']:
            if run('validation/'+name, [str(args.build/f'bin/KelvinBall{name}Test')], .25):
                raise RuntimeError(name+' control failed')
        restart = [f'--restart-mesh={mesh}', f'--restart-target-volume={target:.17g}', '--initial-global-step=10']
        common = [str(binary), *PARAMETERS, *restart]
        if run('validation/geometry', [*common, '--thickness-min=0', '--geometry-only'], .25):
            raise RuntimeError('Restart mesh load failed')
        loaded = args.output/'validation/geometry/RestartGeometry.Chamber.mesh.000000.h5'
        record['restart_geometry_exact'] = same_geometry(snapshot, loaded)
        if run('validation/state', [*common, '--thickness-min=0', '--state-only', '--iterations=1'], 1):
            raise RuntimeError('Restart state solve failed')
        row = complete_rows(args.output/'validation/state/kelvin-ball.csv')[0]
        record['initial_state_regression'] = {key: dict(original=float(original[key]), restart=float(row[key]))
            for key in ['k', 'c', 'q', 'rho', 'volume', 'target_volume']}
        if not all(math.isclose(float(row[key]), float(original[key]), rel_tol=1e-8, abs_tol=1e-11)
                   for key in record['initial_state_regression']):
            raise ValueError('Restart did not recover the original global10 FEM state')
        thick = ['--thickness-min=2', '--thickness-method=baseline']
        if run('validation/direction', [*common, *thick, '--iterations=1'], 1):
            raise RuntimeError('Thickness identification failed')
        direction = complete_rows(args.output/'validation/direction/kelvin-ball.csv')[0]
        check_direction(direction)
        record['direction_validation'] = {key: direction[key] for key in direction
            if key.startswith('thickness_') or key in ['d_volume_theta', 'd_thickness_theta']}
        record['state'] = 'optimizing'
        native_exit = run('optimization', [*common, *thick, '--iterations=6'], 12)
        record['native_exit_code'] = native_exit
        rows = complete_rows(args.output/'optimization/kelvin-ball.csv')
        for i, row in enumerate(rows):
            check_direction(row)
            if int(row['global_iteration']) != i+10 or float(row['target_volume']) != target:
                raise ValueError('Restart global numbering or original volume reference changed')
        record.update(fully_evaluated_states=len(rows), completed_real_updates=max(0, len(rows)-1),
                      budget_completed=native_exit == 0 and len(rows) == 6)
        if not rows:
            raise ValueError('No complete state available for final independent BEM')
        local = len(rows)-1
        final = args.output/'optimization'/f'KelvinBall.Chamber.mesh.{local:06d}.h5'
        record.update(selected_bem_local_step=local, selected_bem_global_step=10+local,
            selected_bem_geometry_fem_evaluated=True, state='exporting_final')
        directory = args.output/'final-geometry'
        directory.mkdir()
        native = directory/'native-export.mesh'
        mesh = directory/'original-double.mesh'
        export.export_snapshot(final, native)
        record['final_certificate'] = export.original_double_mesh(final, native, mesh)
        panels = int(rows[-1]['interface_triangles'])
        if record['final_certificate']['gamma_panels'] != panels:
            raise ValueError('Final export lost original panels')
        jobs = [dict(name='cubic-final', symmetry='cubic', input=str(mesh), sha256=digest(mesh),
                     panels=panels, full_panels=24*panels, saved_step=10+local)]
        jobs_file = args.output/'bem-jobs.json'
        save(jobs_file, jobs)
        record['state'] = 'independent_bem'
        command = [sys.executable, '-B', str(args.cpu_source/'scripts/run_saved_cpu_queue.py'),
            '--jobs', str(jobs_file), '--output', str(args.output/'bem'), '--verifier', str(args.verifier),
            '--library', str(args.library), '--commit', CPU, '--threads', '12']
        if run('bem-wrapper', command, 73):
            raise RuntimeError('Independent final BEM failed')
        accepted = module('cpu_acceptance', args.cpu_source/'scripts/run_saved_cpu_queue.py')
        record['independent_bem'] = accepted.accepted(args.output/'bem/cubic-final', jobs[0], args.library, CPU)
        if digest(snapshot) != INPUT or digest(args.output/'restart/restart-labelled.mesh') != record['restart_mesh_sha256']:
            raise ValueError('Original restart input changed during run')
        pinned(args.source, args.commit)
        record['state'] = 'completed' if record['budget_completed'] else 'native_failed_bem_completed'
        code = 0 if record['budget_completed'] else 1
    except Exception as error:
        record.update(state='failed', error=repr(error))
    finally:
        record.update(exit_code=code, finished_at=time.time())
        checkpoint()
        (args.output/'exit-code.txt').write_text(str(code)+'\n')
    return code


if __name__ == '__main__':
    raise SystemExit(main())
