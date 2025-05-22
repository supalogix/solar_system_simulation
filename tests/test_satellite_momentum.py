import pytest
import numpy as np

from models import Satellite
from constants import G, solar_mass

def test_zero_angle_momentum():
    """
    Zero: at angle = 0, velocity = [0, r*omega, 0], so momentum = m * [0, r*omega, 0].
    """
    r = 1e7
    m = 5.0e3
    sat = Satellite(name="ZeroMom", distance_from_sun=r, mass=m)
    omega = np.sqrt(G.value * solar_mass.value / r**3)
    mom = sat.momentum(0.0)
    expected = np.array([0.0, m * r * omega, 0.0])
    assert np.allclose(mom, expected, rtol=1e-12)


def test_one_quarter_turn_momentum():
    """
    One: at angle = \pi/2, momentum should point almost entirely in −x.
    We check the x‐component directly, and assert the y‐component is negligible.
    """
    r = 2e7
    m = 4.2e3
    sat = Satellite(name="OneMom", distance_from_sun=r, mass=m)
    omega = np.sqrt(G.value * solar_mass.value / r**3)

    mom = sat.momentum(np.pi/2)
    px_expected = -m * r * omega

    # x should match to high relative precision
    assert np.isclose(mom[0], px_expected, rtol=1e-12)

    # y‐component should be many orders smaller than |px|
    assert abs(mom[1]) < abs(px_expected) * 1e-12

    # z remains exactly zero
    assert mom[2] == 0.0


def test_lots_of_angles_momentum():
    """
    Lots: sampling many angles yields a (N×3) array whose magnitudes all equal m*r*omega.
    """
    r = 3.3e7
    m = 7.7e3
    sat = Satellite(name="LotsMom", distance_from_sun=r, mass=m)
    omega = np.sqrt(G.value * solar_mass.value / r**3)
    angles = np.linspace(0, 2*np.pi, 50)
    moms = np.array([sat.momentum(a) for a in angles])
    assert moms.shape == (50, 3)
    speeds = np.linalg.norm(moms, axis=1)
    assert np.allclose(speeds, m * r * omega, rtol=1e-12)


def test_many_full_rotations_momentum_repeats():
    """
    Many: after any integer number of 2*\pi rotations, momentum returns to initial value.
    """
    r = 4e7
    m = 1.1e4
    sat = Satellite(name="ManyMom", distance_from_sun=r, mass=m)
    base_angle = 0.73
    mom0 = sat.momentum(base_angle)
    for k in range(1, 6):
        mom_k = sat.momentum(base_angle + 2*k*np.pi)
        assert np.allclose(mom_k, mom0, rtol=1e-12)


def test_opps_invalid_angle_type():
    """
    Opps: non-numeric angle should raise a TypeError downstream in sin/cos.
    """
    sat = Satellite(name="OopsMom", distance_from_sun=1e7, mass=1.0)
    with pytest.raises(TypeError):
        sat.momentum("not an angle")
