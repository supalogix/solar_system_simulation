import pytest
import numpy as np

from models import Planet
from constants import G, solar_mass

def test_zero_velocity():
    """
    Zero: For a circular orbit (e=0) at t=0, there is no radial component and
    the tangential speed is v = a * n = a * sqrt(G M / a^3) = sqrt(G M / a).
    At nu=0 this should point entirely along +y.
    """
    a = 1e7
    p = Planet(name="Zero", a=a, e=0.0, i=0.0, Omega=0.0, omega=0.0, nu=0.0, mass=1e23)
    v = p.velocity(0.0)
    expected_speed = np.sqrt(G.value * solar_mass.value / a)
    assert np.allclose(v, [0.0, expected_speed, 0.0], atol=1e-8)


def test_one_quarter_period_velocity():
    """
    One: For a circular orbit after T/4, the planet should be at $\nu=\pi/2$,
    so velocity vector = [-v_tan, 0, 0].
    """
    a = 1e7
    p = Planet(name="One", a=a, e=0.0, i=0.0, Omega=0.0, omega=0.0, nu=0.0, mass=1e23)
    T = 2 * np.pi * np.sqrt(a**3 / (G.value * solar_mass.value))
    v = p.velocity(T/4)
    expected_speed = np.sqrt(G.value * solar_mass.value / a)
    assert np.allclose(v, [-expected_speed, 0.0, 0.0], atol=1e-6)


def test_lots_of_velocity_samples():
    """
    Lots: Sampling velocity at many times should return a 3‐vector each time
    and remain finite for a valid elliptical orbit (e<1).
    """
    p = Planet(name="Lots", a=1e7, e=0.3, i=10.0, Omega=20.0, omega=30.0, nu=40.0, mass=1e23)
    times = np.linspace(0, 1e5, 50)
    vels = np.array([p.velocity(t) for t in times])
    # shape must be (N,3)
    assert vels.shape == (50, 3)
    # all entries finite
    assert np.isfinite(vels).all()


def test_many_orbits_velocity_repeats():
    """
    Many: After an integer number of full periods, the velocity should return
    (approximately) to its initial value.
    """
    a = 1e7
    p = Planet(name="Many", a=a, e=0.0, i=0.0, Omega=0.0, omega=0.0, nu=0.0, mass=1e23)
    T = 2 * np.pi * np.sqrt(a**3 / (G.value * solar_mass.value))
    v0 = p.velocity(0.0)
    vN = p.velocity(5 * T)
    assert np.allclose(vN, v0, atol=1e-6)


def test_opps_invalid_eccentricity_velocity():
    """
    Opps: For e >= 1 the orbit is invalid. We expect x and y to become NaN
    (due to sqrt of negative), and z to remain zero.
    """
    p = Planet(name="Oops", a=1e7, e=1.5, i=0.0, Omega=0.0, omega=0.0, nu=0.0, mass=1e23)
    v = p.velocity(0.0)
    assert np.isnan(v[0]) and np.isnan(v[1])
    assert v[2] == 0.0
