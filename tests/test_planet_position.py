import pytest
import numpy as np

from models import Planet
from constants import G, solar_mass

def test_zero_time_position():
    """
    Zero: at t=0 the eccentric anomaly E=0, so the planet sits at
    periapsis on the x-axis (distance = a*(1−e)).
    """
    p = Planet(name="TestZero", a=1e7, e=0.3, i=0.0, Omega=0.0, omega=0.0, nu=0.0, mass=1e23)
    pos = p.position(0.0)
    expected_x = p.a * (1 - p.e)
    assert np.allclose(pos, [expected_x, 0.0, 0.0], atol=1e-8)


def test_one_quarter_orbit():
    """
    One: for a circular orbit (e=0), after one quarter of the period T/4,
    the planet should be at y = +a, x = 0.
    """
    p = Planet(name="TestOne", a=1e7, e=0.0, i=0.0, Omega=0.0, omega=0.0, nu=0.0, mass=1e23)
    # period T = 2 * \pi * sqrt(a^3/( G* M))
    T = 2 * np.pi * np.sqrt(p.a**3 / (G.value * solar_mass.value))
    pos = p.position(T / 4)
    assert np.allclose(pos, [0.0, p.a, 0.0], atol=1e-6)


def test_lots_of_time_samples():
    """
    Lots: computing position at many time samples should return an array of
    the correct shape (N×3).
    """
    p = Planet(name="TestLots", a=1e7, e=0.2, i=10.0, Omega=20.0, omega=30.0, nu=40.0, mass=1e23)
    times = np.linspace(0, 1e5, 50)
    positions = np.array([p.position(t) for t in times])
    assert positions.shape == (50, 3)


def test_many_orbits_returns_initial():
    """
    Many: after an integer number of full periods, the planet should return
    (approximately) to its starting position.
    """
    p = Planet(name="TestMany", a=1e7, e=0.0, i=0.0, Omega=0.0, omega=0.0, nu=0.0, mass=1e23)
    T = 2 * np.pi * np.sqrt(p.a**3 / (G.value * solar_mass.value))
    start = p.position(0.0)
    end   = p.position(10 * T)
    assert np.allclose(end, start, atol=1e-6)


def test_opps_invalid_eccentricity():
    """
    Opps: eccentricity >= 1 is not a valid ellipse.
    We expect the x and y components to be NaN, and z to remain 0.
    """
    p = Planet(name="TestOops",
               a=1e7, e=1.5,
               i=0.0, Omega=0.0,
               omega=0.0, nu=0.0,
               mass=1e23)

    pos = p.position(0.0)
    # x & y should be NaN
    assert np.isnan(pos[0]) and np.isnan(pos[1])
    # z should remain exactly zero (we never alter the z-component)
    assert pos[2] == 0.0

