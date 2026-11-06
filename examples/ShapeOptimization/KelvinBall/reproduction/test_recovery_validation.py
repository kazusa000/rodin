"""Regression guards for refusing changed or incomplete restart evidence."""
import tempfile
import unittest
from pathlib import Path
import h5py
import numpy as np
from recovery_validation import compare_fields


class RecoveryChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.old, self.new = self.root/'old', self.root/'new'
        for root, step in ((self.old, 24), (self.new, 0)):
            root.mkdir()
            for grid in ('Chamber', 'Fluid'):
                with h5py.File(root/f'KelvinBall.{grid}.mesh.{step:06d}.h5', 'w') as f:
                    f['Mesh/Geometry/Vertices'] = np.array([[1., 0., 0.], [1., 1., 0.], [1., 0., 1.], [0., 0., 0.]])
                    f['Mesh/XDMF/Topology'] = np.array([6, 0, 1, 2, 3], dtype='u8')
                    for dim in range(4):
                        f[f'Mesh/Attributes/{dim}'] = np.array([2], dtype='u8')
            names = ['Chamber.Theta', 'Chamber.Distance']
            names += [f'Fluid.{p}{k}_{i}' for p in ('', 'Pressure_') for k in ('Translation', 'Rotation') for i in range(3)]
            for name in names:
                with h5py.File(root/f'KelvinBall.{name}.{step:06d}.h5', 'w') as f:
                    f['GridFunction/Values/Data'] = np.ones((4, 3))

    def change(self, file, dataset, index, value):
        with h5py.File(self.new/file, 'r+') as f:
            f[dataset][index] = value
        with self.assertRaises(ValueError):
            compare_fields(self.old, self.new)

    def test_identical_saved_fields_pass(self):
        self.assertIn('Chamber.Theta', compare_fields(self.old, self.new))

    def test_one_ulp_coordinate_change_rejected(self):
        self.change('KelvinBall.Chamber.mesh.000000.h5', 'Mesh/Geometry/Vertices', (0,0), np.nextafter(1., 2.))

    def test_reordered_tetrahedron_rejected(self):
        self.change('KelvinBall.Chamber.mesh.000000.h5', 'Mesh/XDMF/Topology', 1, 1)

    def test_material_change_rejected(self):
        self.change('KelvinBall.Fluid.mesh.000000.h5', 'Mesh/Attributes/3', 0, 3)

    def test_face_label_change_rejected(self):
        self.change('KelvinBall.Chamber.mesh.000000.h5', 'Mesh/Attributes/2', 0, 35)

    def test_changed_theta_rejected(self):
        self.change('KelvinBall.Chamber.Theta.000000.h5', 'GridFunction/Values/Data', (0,0), 1.01)

    def test_changed_distance_rejected(self):
        self.change('KelvinBall.Chamber.Distance.000000.h5', 'GridFunction/Values/Data', (0,0), 1.00000001)

    def test_missing_flow_field_rejected(self):
        (self.new/'KelvinBall.Fluid.Pressure_Rotation_2.000000.h5').unlink()
        with self.assertRaises(ValueError):
            compare_fields(self.old, self.new)


if __name__ == '__main__':
    unittest.main()
