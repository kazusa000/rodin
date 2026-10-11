from pathlib import Path
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
import sys
import numpy as np
from kelvin_bem.audit_clement_reconstruction import quality, evaluate


def test_regular_tetrahedron_quality_is_one():
    x = np.array([[0,0,0],[1,0,0],[.5,np.sqrt(3)/2,0],[.5,np.sqrt(3)/6,np.sqrt(2/3)]])
    np.testing.assert_allclose(quality(x,np.array([[0,1,2,3]])),1,rtol=1e-14)


def test_original_p1_interpolation_and_no_nearest_extrapolation():
    x = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float)
    values = 2*x[:,0]-3*x[:,1]+4*x[:,2]+.1
    result = evaluate(x,np.array([[0,1,2,3]]),values,np.array([[.2,.1,.3],[2,2,2]]))
    np.testing.assert_allclose(result[0],1.4,rtol=1e-14)
    assert np.isnan(result[1])
