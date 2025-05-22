r"""
This module provides support functions to calculate relativistic time dilation
effects for a given set of massive bodies and observer positions. It exports
utilities that, given the mass, position, and momentum of each planet, compute
how much slower a clock ticks compared to one far from any gravitational source,
returning the result in attoseconds (10^{-18} s).

Time dilation arises from two distinct influences:

1. Gravitational time dilation  
   A mass “bends” spacetime, causing clocks closer to it to tick more slowly.
   In the weak-field approximation, the fractional slowdown is
   $\displaystyle \frac{G\,m}{r\,c^2}$, where
   $m$ is the body’s mass, $r$ is the distance to the clock, $G$ is
   Newton’s constant, and $c$ is the speed of light.

2. Kinematic time dilation  
   A moving clock also ticks more slowly than one at rest.  For speeds
   $v \ll c$, the tiny fractional difference per second is
   $\displaystyle \Delta t - \tau \approx \frac{1}{2}\,\frac{v^2}{c^2}$,
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
        satellites (dict of str -> array-like): Mapping of satellite point names (e.g., "L1")
            to their 3D observer positions.
        t (float, optional): Time in seconds at which to sample each planet's state.
            Defaults to 0.

    Returns:
        dict of str -> dict of str -> float:
            Nested dictionary where each satellite key maps to a dict of planet names
            and their corresponding time dilation (in attoseconds).
    """
    sat_keys = list(satellites.keys())
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


        for key, obs_pos in satellites.items():
            # Distance for gravitational dilation
            d = np.linalg.norm(obs_pos - planet_pos)

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
            result[key][planet.name] = float((td_grav + td_kin) * afs_conversion)

    return result

def calculate_time_dilation(planets, position, t=0):
    """
    Return the combined gravitational and kinematic time dilation at the given
    observer position and time, in attoseconds, by summing each planet’s
    contribution computed by calculate_time_dilation_per_planet.

    Args:
        planets (list of Planet): Planet objects with mass, position(), and momentum().
        position (array-like): 3D coordinates of the observer.
        t (float, optional): Time in seconds to sample each planet’s state. Defaults to 0.

    Returns:
        float: Net time dilation offset per second, in attoseconds.
    """
    # Compute per‐planet dilation at a single “dummy” satellite key
    per_planet = calculate_time_dilation_per_planet(
        planets,
        satellites={ "_observer": position },
        t=t
    )
    # Extract the dict for our dummy key and sum all planet contributions
    return sum(per_planet["_observer"].values())