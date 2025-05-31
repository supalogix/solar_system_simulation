r"""
This module provides support functions to calculate relativistic time dilation
effects for a given set of massive bodies and observer positions. It exports
utilities that, given the mass, position, and momentum of each planet, compute
how much slower a clock ticks compared to one far from any gravitational source,
and returns the result in attoseconds (10^{-18} s).

Time dilation arises from two distinct influences:

1. Gravitational time dilation  
   A mass “bends” spacetime, causing clocks closer to it to tick more slowly.
   In the weak-field approximation, the fractional slowdown is
   $\frac{G\,m}{r\,c^2}$, where
   $m$ is the body’s mass, $r$ is the distance to the clock, $G$ is
   Newton’s constant, and $c$ is the speed of light.

2. Kinematic time dilation  
   A moving clock also ticks more slowly than one at rest.  For speeds
   $v \ll c$, the tiny fractional difference per second is
   $\Delta t - \tau \approx \frac{1}{2}\,\frac{v^2}{c^2}$,
   where $v$ is the object’s speed.  Note that even if a body is very far
   away (so its gravitational effect is negligible), its high speed will
   still produce a kinematic slowdown detectable here.

By summing these two contributions for all planets and then converting the
combined fractional rate into attoseconds, the functions in this file deliver
a convenient, high-precision measure of time dilation experienced by an
observer in any specified location.
"""

import mpmath as mp
import numpy as np
from constants import G_mp, c_mp, afs_conversion

# Use high precision so tiny relativistic effects are captured
mp.dps = 50

def calculate_time_dilation_per_planet(planets, satellites, t=0):
    """
    Return gravitational and kinematic time dilation in attoseconds for each planet
    at each satellite observation point.

    Args:
        planets (list of Planet): Planet objects providing mass, position(t), and momentum(t).
        satellites (list of Planet): Planet objects providing mass, position(t), and momentum(t).
        t (float, optional): Time in seconds at which to sample each planet's state.
            Defaults to 0.

    Returns:
        dict of str -> dict of str -> float:
            Nested dictionary where each satellite key maps to a dict of planet names
            and their corresponding time dilation (in attoseconds).
    """

    sat_keys = [s.name for s in satellites]
    names    = [p.name for p in planets]

    # Prepare a container to accumulate each satellite’s time-dilation contributions
    # for every planet, initializing all values to zero before we start summing.
    result = {
        key: {name: 0.0 for name in names}
        for key in sat_keys
    }


    for planet in planets:
        # High-precision mass m
        m = mp.mpf(planet.mass)

        # Planet position for gravitational term
        planet_pos = planet.position(t)

        # We need the planet’s squared speed to know how its motion slows its clock.
        # Since we already have its momentum (mass * speed), we derive v^2 from that.
        # Computing it once here lets us reuse the value for each satellite point,
        # saving work and keeping the code efficient.
        p_vec = planet.momentum(t)
        p_sq  = mp.mpf(np.dot(p_vec, p_vec))
        v_sq  = p_sq / (m**2)


        for satellite in satellites:
            # Distance for gravitational dilation
            d = np.linalg.norm(satellite.position(t) - planet_pos)

            if d == 0:
                continue

            # Gravitational slowdown: mass warps spacetime and causes time to flow more slowly
            # closer to the planet. Clocks near the planet tick slightly slower than those
            # farther away. In the weak gravitational field of a planet, we can approximate
            # this tiny effect with a simple formula.
            td_grav = (G_mp * m) / (d * c_mp**2)

            # Kinematic slowdown: according to relativity, a moving clock ticks more slowly
            # compared to one at rest. Even at everyday speeds, this tiny delay grows with
            # the square of the speed, so we include a small correction for the planet’s motion.
            td_kin = v_sq / (2 * c_mp**2)

            # Convert to attoseconds and store
            result[satellite.name][planet.name] = float((td_grav + td_kin) * afs_conversion)

    return result

def generate_simulation(planets, satellites, t=0):
    """
    Generate a snapshot log of all planets and satellites at time t, including positions,
    velocities, momenta, and time‐dilation info for each satellite.

    Parameters
    ----------
    planets : list of Planet
        A list of Planet objects representing the primary bodies (Sun, Mercury, Venus, etc.).
    satellites : list of Planet
        A list of Planet objects representing test‐particles or Trojan/Lagrange‐point satellites.
    t : float, optional
        Simulation time in seconds at which to evaluate each body's state. Default is 0.

    Returns
    -------
    dict
        A nested dictionary with the following structure:

        {
            "time": float,  # equals t
            "planets": {
                "<planet_name>": {
                    "position": { "px": float, "py": float, "pz": float },
                    "velocity": { "vx": float, "vy": float, "vz": float },
                    "momentum": { "mx": float, "my": float, "mz": float }
                },
                ...
            },
            "satellites": {
                "<satellite_name>": {
                    "position": { "px": float, "py": float, "pz": float },
                    "velocity": { "vx": float, "vy": float, "vz": float },
                    "momentum": { "mx": float, "my": float, "mz": float },
                    "time_dilation_per_planet": {
                        "<planet_name>": float,  # time dilation contribution in seconds
                        ...
                    },
                    "time_dilation_total": float  # sum of all per‐planet time dilations
                },
                ...
            }
        }
    """

    log = {}
    log['time'] = t
    log['planets'] = {}
    log['satellites'] = {}

    for p in planets:
        px, py, pz = p.position(t)
        vx, vy, vz = p.velocity(t)
        mx, my, mz = p.momentum(t)
        log['planets'][p.name] = {}
        log['planets'][p.name]['position'] = {'px': px, 'py': py, 'pz': pz}
        log['planets'][p.name]['velocity'] = {'vx': vx, 'vy': vy, 'vz': vz}
        log['planets'][p.name]['momentum'] = {'mx': mx, 'my': my, 'mz': mz}

    for s in satellites:
        px, py, pz = s.position(t)
        vx, vy, vz = s.velocity(t)
        mx, my, mz = s.momentum(t)
        log['satellites'][s.name] = {}
        log['satellites'][s.name]['position'] = {'px': px, 'py': py, 'pz': pz}
        log['satellites'][s.name]['velocity'] = {'vx': vx, 'vy': vy, 'vz': vz}
        log['satellites'][s.name]['momentum'] = {'mx': mx, 'my': my, 'mz': mz}

    time_dilations_per_planet = calculate_time_dilation_per_planet(
        planets=planets, 
        satellites=satellites, 
        t=t
    )

    sums = { L: sum(vals.values()) for L, vals in time_dilations_per_planet.items() }

    for name, val in time_dilations_per_planet.items():
        log['satellites'][name]['time_dilation_per_planet'] = val
        log['satellites'][name]['time_dilation_total'] = sums[name]

    return log
