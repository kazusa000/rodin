import sys
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
from pathlib import Path
import numpy as np
import pytest

from kelvin_bem.tetra_solid_connectivity import solid_connectivity, triangle_overlap_area
from kelvin_bem import export_carlos_interface as exporter


def chamber():
    x=np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,1],[1,1,-1]],float)
    t=np.array([[0,1,2,3],[0,1,2,4]])
    return x,t


def test_overlap_handles_nonmatching_faces_but_not_point_or_edge_contacts():
    a=np.array([[0,0],[1,0],[0,1]],float)
    assert triangle_overlap_area(a,a[::-1])==pytest.approx(.5)
    assert triangle_overlap_area(a,np.array([[0,0],[.5,0],[0,.5]]))==pytest.approx(.125)
    assert triangle_overlap_area(a,np.array([[1,0],[2,0],[1,1]]))==0
    assert triangle_overlap_area(a,np.array([[0,0],[1,0],[0,-1]]))==0


def test_rotated_cuts_connect_twelve_chambers_without_welding():
    x,t=chamber(); original=x.copy()
    d=solid_connectivity(x,t,np.full(len(t),2),exporter.geometry.proper_tetrahedral_rotations())
    assert d['chamber_components']==1 and d['full_components']==1
    assert d['positive_area_cut_joins']==11
    assert sum(d['component_volumes_full'])==pytest.approx(4)
    assert np.array_equal(x,original)


@pytest.mark.parametrize('touch',[False,True])
def test_internal_detached_or_vertex_touching_tetra_is_a_fragment(touch):
    x,t=chamber()
    if touch:
        x=np.vstack((x,[[1.2,1.2,1.1],[1.3,1.2,1.1],[1.2,1.3,1.1]]))
        t=np.vstack((t,[3,5,6,7]))
    else:
        x=np.vstack((x,[[1.5,1.5,0],[1.6,1.5,0],[1.5,1.6,0],[1.5,1.5,.1]]))
        t=np.vstack((t,[5,6,7,8]))
    d=solid_connectivity(x,t,np.full(len(t),2),exporter.geometry.proper_tetrahedral_rotations())
    assert d['chamber_components']==2 and d['full_components']==13


def test_nonmanifold_solid_and_empty_solid_are_invalid():
    x,t=chamber(); rotations=exporter.geometry.proper_tetrahedral_rotations()
    with pytest.raises(ValueError,match='empty'):
        solid_connectivity(x,t,np.full(len(t),3),rotations)
    with pytest.raises(ValueError,match='nonmanifold'):
        solid_connectivity(x,np.vstack((t,t[0])),np.full(3,2),rotations)


def test_cubic_self_paired_cut_connects_all_24_and_detects_island():
    x=np.array([[0,0,0],[1,0,0],[1,1,1],[1,1,-1]],float)
    t=np.array([[0,1,2,3]])
    rotations=exporter.geometry.proper_octahedral_rotations()
    d=solid_connectivity(x,t,np.full(1,2),rotations,symmetry='cubic')
    assert d['full_components']==1 and d['positive_area_cut_joins']==23
    assert sum(d['component_volumes_full'])==pytest.approx(8)
    x=np.vstack((x,[[.7,.3,0],[.8,.3,0],[.7,.4,0],[.7,.3,.1]]))
    t=np.vstack((t,[4,5,6,7]))
    broken=solid_connectivity(x,t,np.full(2,2),rotations,symmetry='cubic')
    assert broken['chamber_components']==2 and broken['full_components']==25
    with pytest.raises(ValueError,match='symmetry/group mismatch'):
        solid_connectivity(x,t,np.full(2,2),rotations,symmetry='terra')


def test_fixed_cubic_seed_is_proper_rotation_invariant_not_mirror_invariant():
    def radius(x):
        a,b,c=np.asarray(x)/np.linalg.norm(x)
        h3=3*np.sqrt(3)*a*b*c
        h6=6*np.sqrt(3)*(a*a-b*b)*(b*b-c*c)*(c*c-a*a)
        return 1+.02*h3*h6
    x=np.array([.21,.53,.82])
    for r in exporter.geometry.proper_octahedral_rotations():
        assert radius(np.asarray(r)@x)==pytest.approx(radius(x),abs=1e-15)
    assert abs(radius(x)-radius(x*np.array([-1,1,1])))>1e-4
