import pytest
import numpy as np

from models import Planet
from constants import G, solar_mass

def test_zero_momentum():
    """
    Zero: For a circular orbit (e=0) at t=0, the planet has zero radial speed
    and tangential speed v = sqrt(G * M / a), so momentum = m * v along +y.
    """
    a = 1e7
    m = 2.0e23
    p = Planet("ZeroMom", a=a, e=0.0, i=0, Omega=0, omega=0, nu=0, mass=m)
    mom = p.momentum(0.0)
    expected_speed = np.sqrt(G.value * solar_mass.value / a)
    expected = np.array([0.0, m * expected_speed, 0.0])

    # allow a small relative tolerance
    assert np.allclose(mom, expected, rtol=1e-12)


def test_one_quarter_period_momentum():
    """
    One: For a circular orbit after T/4, momentum should lie almost entirely
    along −x. We check the x‐component directly, and only assert that the tiny
    y‐leakage is negligible compared to abs(px).
    """
    a = 1e7
    m = 3.5e23
    p = Planet("OneMom", a=a, e=0.0, i=0, Omega=0, omega=0, nu=0, mass=m)
    T = 2 * np.pi * np.sqrt(a**3 / (G.value * solar_mass.value))

    mom = p.momentum(T/4)
    # expected main component
    expected_speed = np.sqrt(G.value * solar_mass.value / a)
    px_expected = -m * expected_speed

    # x should match within relative tolerance
    assert np.isclose(mom[0], px_expected, rtol=1e-12)

    # y should be many orders smaller than abs(px)
    assert abs(mom[1]) < abs(px_expected) * 1e-12

    # z is exactly zero
    assert mom[2] == 0.0


def test_lots_of_momentum_samples():
    """
    Lots: Sampling momentum at many times should return a 3‐vector each time
    and remain finite when e<1.
    """
    p = Planet("LotsMom", a=1e7, e=0.4, i=5, Omega=15, omega=25, nu=35, mass=4.2e23)
    times = np.linspace(0, 2e5, 60)
    moms = np.array([p.momentum(t) for t in times])
    assert moms.shape == (60, 3)
    assert np.isfinite(moms).all()


def test_many_orbits_momentum_repeats():
    """
    Many: After an integer number of full periods, momentum returns to its
    initial direction and magnitude. We compare the dominant y‐component,
    and ensure any x‐leakage remains negligible.
    """
    a = 1e7
    m = 5.1e23
    p = Planet("ManyMom", a=a, e=0.0, i=0, Omega=0, omega=0, nu=0, mass=m)
    T = 2 * np.pi * np.sqrt(a**3 / (G.value * solar_mass.value))

    mom0 = p.momentum(0.0)
    momN = p.momentum(7 * T)

    # y‐component should match within relative tolerance
    assert np.isclose(momN[1], mom0[1], rtol=1e-12)

    # any x‐component must be negligible compared to abs(py)
    assert abs(momN[0]) < abs(mom0[1]) * 1e-12

    # z remains zero
    assert momN[2] == 0.0


def test_opps_invalid_eccentricity_momentum():
    """
    Opps: For e >= 1 the orbit is invalid. x and y momentum should be NaN,
    z remains zero.
    """
    p = Planet("OopsMom", a=1e7, e=1.2, i=0, Omega=0, omega=0, nu=0, mass=6.3e23)
    mom = p.momentum(0.0)
    assert np.isnan(mom[0]) and np.isnan(mom[1])
    assert mom[2] == 0.0
