import pytest
import numpy as np

from models import Satellite
from constants import G, solar_mass

def test_zero_angle_velocity():
    """
    Zero: at angle = 0, vx = 0 and vy = r * omega, where
    omega = sqrt(G*M / r^3). Velocity should point purely along +y.
    """
    r = 1e7
    sat = Satellite(name="ZeroVel", distance_from_sun=r)
    omega = np.sqrt(G.value * solar_mass.value / r**3)
    v = sat.velocity(0.0)
    expected = np.array([0.0, r * omega, 0.0])
    assert np.allclose(v, expected, rtol=1e-12)


def test_one_quarter_turn_velocity():
    """
    One: at angle = pi/2, vx = -r * omega and vy = 0. Velocity should point purely
    in the -x direction.
    """
    r = 2e7
    sat = Satellite(name="OneQVel", distance_from_sun=r)
    omega = np.sqrt(G.value * solar_mass.value / r**3)
    v = sat.velocity(np.pi/2)
    expected = np.array([-r * omega, 0.0, 0.0])
    assert np.allclose(v, expected, rtol=1e-12)


def test_lots_of_angles_velocity():
    """
    Lots: sampling many angles should return a 3‐vector each time and the speed
    magnitude abs(v) should remain constant = r * omega.
    """
    r = 5e6
    sat = Satellite(name="LotsVel", distance_from_sun=r)
    omega = np.sqrt(G.value * solar_mass.value / r**3)
    angles = np.linspace(0, 2*np.pi, 50)
    vels = np.array([sat.velocity(a) for a in angles])
    assert vels.shape == (50, 3)
    speeds = np.linalg.norm(vels, axis=1)
    assert np.allclose(speeds, r * omega, rtol=1e-12)


def test_many_full_rotations_velocity():
    """
    Many: after any integer number of full 2*\pi rotations, velocity should
    return to its value at the base angle.
    """
    r = 3e7
    sat = Satellite(name="ManyVel", distance_from_sun=r)
    base_angle = 1.23
    v0 = sat.velocity(base_angle)
    for k in range(1, 6):
        v_k = sat.velocity(base_angle + 2*k*np.pi)
        assert np.allclose(v_k, v0, rtol=1e-12)


def test_opps_invalid_angle_type():
    """
    Opps: non-numeric angles should raise a TypeError.
    """
    sat = Satellite(name="OopsVel", distance_from_sun=1e7)
    with pytest.raises(TypeError):
        sat.velocity("not a number")
