import pytest
import numpy as np

from models import Satellite

def test_zero_angle_position():
    """
    Zero: at angle = 0 the satellite should lie on the +x axis at (r, 0, 0).
    """
    r = 1e7
    sat = Satellite(name="Zero", distance_from_sun=r)
    pos = sat.position(0.0)
    assert np.allclose(pos, [r, 0.0, 0.0], atol=1e-8)


def test_one_quarter_turn_position():
    """
    One: at angle = \pi/2 the satellite should lie on the +y axis at (0, r, 0).
    """
    r = 2e7
    sat = Satellite(name="OneQ", distance_from_sun=r)
    pos = sat.position(np.pi/2)
    assert np.allclose(pos, [0.0, r, 0.0], atol=1e-8)


def test_lots_of_angles_position():
    """
    Lots: sampling many angles should return an (N × 3) array of positions,
    each exactly at distance r from the origin.
    """
    r = 5.5e6
    sat = Satellite(name="Lots", distance_from_sun=r)
    angles = np.linspace(0, 2*np.pi, 100)
    positions = np.array([sat.position(a) for a in angles])
    assert positions.shape == (100, 3)
    # Check every point lies on the circle of radius r
    dists = np.linalg.norm(positions, axis=1)
    assert np.allclose(dists, r, rtol=1e-12)


def test_many_full_rotations_position():
    """
    Many: after any integer number of full 2*\pi rotations, position should
    return to its starting point (within floating-point tolerance).
    """
    r = 3e7
    sat = Satellite(name="Many", distance_from_sun=r)
    start = sat.position(0.37)
    for k in range(1, 6):
        pos_k = sat.position(0.37 + 2*k*np.pi)
        assert np.allclose(pos_k, start, rtol=1e-12)


def test_opps_invalid_angle_type():
    """
    Opps: non-numeric angles should raise a TypeError.
    """
    sat = Satellite(name="Oops", distance_from_sun=1e7)
    with pytest.raises(TypeError):
        sat.position("not a number")
