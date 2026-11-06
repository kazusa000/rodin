"""Compare recovered fields at original vertex indices before allowing new updates."""
from pathlib import Path
import h5py
import numpy as np


def compare_fields(original, recovered, step=24):
    original, recovered = Path(original), Path(recovered)
    evidence = {}
    for grid in ('Chamber', 'Fluid'):
        old = original / f'KelvinBall.{grid}.mesh.{step:06d}.h5'
        new = recovered / f'KelvinBall.{grid}.mesh.000000.h5'
        with h5py.File(old) as a, h5py.File(new) as b:
            for name in ('Mesh/Geometry/Vertices', 'Mesh/XDMF/Topology', 'Mesh/Attributes/3'):
                if not np.array_equal(a[name][:], b[name][:]):
                    raise ValueError(f'Original {grid} coordinates, tetrahedron order or materials changed: {name}')
            for dim in range(3):
                if not np.array_equal(np.sort(a[f'Mesh/Attributes/{dim}'][:]),
                                      np.sort(b[f'Mesh/Attributes/{dim}'][:])):
                    raise ValueError(f'Original {grid} labels changed in dimension {dim}')
        evidence[grid + '_geometry'] = 'exact coordinates, ordered tetrahedra, materials; label histograms match'
    names = ['Chamber.Theta', 'Chamber.Distance']
    names += [f'Fluid.{pressure}{kind}_{i}' for pressure in ('', 'Pressure_')
              for kind in ('Translation', 'Rotation') for i in range(3)]
    for name in names:
        # Upstream field names are collected by the controller from the original XDMF files.
        old = original / f'KelvinBall.{name}.{step:06d}.h5'
        new = recovered / f'KelvinBall.{name}.000000.h5'
        if not old.exists() or not new.exists():
            raise ValueError(f'Missing complete saved field: {name}')
        with h5py.File(old) as a, h5py.File(new) as b:
            x, y = a['GridFunction/Values/Data'][:], b['GridFunction/Values/Data'][:]
        rtol, atol = ((1e-12, 1e-12) if name.endswith('Distance') else (1e-7, 1e-9))
        if x.shape != y.shape or not np.isfinite(y).all() or not np.allclose(x, y, rtol=rtol, atol=atol):
            maximum = float(np.max(np.abs(x-y))) if x.shape == y.shape else None
            raise ValueError(f'Restored {name} differs; maximum absolute difference={maximum}')
        evidence[name] = dict(max_abs=float(np.max(np.abs(x-y))), rtol=rtol, atol=atol)
    return evidence


def compare_native(a, b):
    with h5py.File(a) as x, h5py.File(b) as y:
        def datasets(file):
            names = []
            file.visititems(lambda name, item: names.append(name) if isinstance(item, h5py.Dataset) else None)
            return names
        names = datasets(x)
        if names != datasets(y) or any(not np.array_equal(x[n][()], y[n][()]) for n in names):
            raise ValueError('Full native checkpoint roundtrip changed a dataset')
    return len(names)
