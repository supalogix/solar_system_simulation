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


# Gravitational slowdown: mass warps spacetime and causes time to flow more slowly
# closer to the planet. Clocks near the planet tick slightly slower than those
# farther away. In the weak gravitational field of a planet, we can approximate
# this tiny effect with a simple formula.
def compute_gravitation_slowdown(G_mp, m, d, c_mp):
    td_grav = (G_mp * m) / (d * c_mp**2)

    return td_grav

# Kinematic slowdown: according to relativity, a moving clock ticks more slowly
# compared to one at rest. Even at everyday speeds, this tiny delay grows with
# the square of the speed, so we include a small correction for the planet’s motion.
def compute_kinematic_slowdown(v_sq, c_mp):
    td_kin = v_sq / (2 * c_mp**2)

    return td_kin

# ---- Newton–Raphson root finder ----
def newton_raphson(F, x0, dx=1e-6, max_iter=20, tol=1e-12):
    """
    Generic 1D Newton–Raphson solver.
    F: function of a single variable.
    x0: initial guess.
    dx: finite‐difference step for derivative.
    """
    x = x0
    for _ in range(max_iter):
        f = F(x)
        # approximate derivative
        df = (F(x + dx) - f) / dx
        if abs(df) < 1e-16:
            break
        dx_x = f / df
        x -= dx_x
        if abs(dx_x) < tol:
            break
    return x

# ---- Lagrange‐point finder ----
def find_collinear_lagrange(planet_name, planets, t=0):
    """
    Compute instantaneous L1, L2, L3 for the given planet.

    Args:
        planet_name (str): Name of the secondary body (e.g. 'Jupiter').
        planets (list of Planet): Must include 'Sun' and the named planet.
        t (float): Epoch time for sampling positions.

    Returns:
        dict: {'L1': np.array([x,y,z]), 'L2': ..., 'L3': ...}
    """
    # Extract Sun and target planet
    sun = next(p for p in planets if p.name == 'Sun')
    sec = next(p for p in planets if p.name == planet_name)

    # Masses
    M_S = mp.mpf(sun.mass)
    M_P = mp.mpf(sec.mass)

    # Positions
    r_S = np.array(sun.position(t))
    r_P = np.array(sec.position(t))

    # Relative vector and distance
    R_vec = r_P - r_S
    R = mp.sqrt((R_vec**2).sum())

    # Instantaneous mean motion
    G = mp.mpf(G_mp)
    mu_T = G * (M_S + M_P)
    Omega = mp.sqrt(mu_T / R**3)

    # Define the force‐balance functions
    def F12(r):
        return G*M_S/r**2 - G*M_P/(R - r)**2 - (Omega**2)*r

    def F3(r):
        return G*M_S/r**2 + G*M_P/(R + r)**2 - (Omega**2)*r

    # Initial guesses (collinear approx)
    mu = M_P/(M_S + M_P)
    r1_guess = R * (1 - (mu/3)**(1/3))
    r2_guess = R * (1 + (mu/3)**(1/3))
    r3_guess = R * (1 + 7*mu/12)

    # Solve for the three roots
    r1 = newton_raphson(F12, r1_guess)
    r2 = newton_raphson(F12, r2_guess)
    r3 = newton_raphson(F3,  r3_guess)

    # Convert to Cartesian positions
    unit = np.array(R_vec / float(R))
    L1 = r_S + float(r1) * unit
    L2 = r_S + float(r2) * unit
    L3 = r_S - float(r3) * unit

    return {'L1': L1, 'L2': L2, 'L3': L3}

def find_equilateral_lagrange(planet_name, planets, t=0):
    """
    Compute instantaneous L4 and L5 for the given planet by 
    phase‐shifting its Sun‐orbit by ±60° in the instantaneous orbital plane.

    Args:
        planet_name (str): Name of the secondary body (e.g. 'Jupiter').
        planets (list of Planet): Must include 'Sun' and the named planet.
        t (float): Epoch time for sampling positions and velocities.

    Returns:
        dict: {'L4': np.array([x,y,z]), 'L5': np.array([x,y,z])}
    """
    # 1) find Sun and secondary
    sun = next(p for p in planets if p.name == 'Sun')
    sec = next(p for p in planets if p.name == planet_name)

    # 2) get their state vectors at t
    r_S = np.array(sun.position(t))
    r_P = np.array(sec.position(t))
    v_P = np.array(sec.velocity(t))

    # 3) define the instantaneous orbital plane axes
    R_vec = r_P - r_S
    R = np.linalg.norm(R_vec)
    x_hat = R_vec / R

    # angular‐momentum vector gives plane normal
    L_vec = np.cross(R_vec, v_P)
    norm_L = np.linalg.norm(L_vec)
    if norm_L < 1e-12:
        raise ValueError("Degenerate orbit: can't define plane normal")
    z_hat = L_vec / norm_L

    # in‐plane orthonormal
    y_hat = np.cross(z_hat, x_hat)

    # 4) rotate by ±60° in that plane
    cos60, sin60 = 0.5, np.sqrt(3)/2

    L4_dir =  x_hat * cos60 + y_hat * sin60
    L5_dir =  x_hat * cos60 - y_hat * sin60

    # 5) position = Sun + R * direction
    L4 = r_S + R * L4_dir
    L5 = r_S + R * L5_dir

    return {'L4': L4, 'L5': L5}

def find_all_lagrange_points(planet_name, planets, t=0):
    """
    Compute all five Sun–planet Lagrange points (L1–L5) at time t.

    Args:
        planet_name (str): Name of the secondary body, e.g. 'Jupiter'.
        planets (list of Planet): Must include the Sun and the named planet.
        t (float): Epoch time for sampling positions/velocities.

    Returns:
        dict: {
            'L1': np.array([x,y,z]),
            'L2': np.array([x,y,z]),
            'L3': np.array([x,y,z]),
            'L4': np.array([x,y,z]),
            'L5': np.array([x,y,z]),
        }
    """
    # collinear points
    collinear = find_collinear_lagrange(planet_name, planets, t)
    # equilateral points
    equilateral = find_equilateral_lagrange(planet_name, planets, t)

    # merge and return
    return {
        'L1': collinear['L1'],
        'L2': collinear['L2'],
        'L3': collinear['L3'],
        'L4': equilateral['L4'],
        'L5': equilateral['L5'],
    }



def calculate_time_dilation_per_planet(planets, planet_name, t=0):
    """
    Return gravitational and kinematic time dilation in attoseconds
    for each Sun–planet Lagrange point of the given planet.

    Args:
        planets (list of Planet): Planet objects providing mass, position(t), velocity(t), and momentum(t).
        planet_name (str): Name of the primary planet (e.g. 'Jupiter').
        t (float, optional): Epoch time in seconds at which to sample each planet's state. Defaults to 0.

    Returns:
        dict of str -> dict of str -> float:
            Outer keys are 'L1'–'L5'; inner dict maps each planet name
            to its time dilation (in attoseconds) at that Lagrange point.
    """
    # find the five Lagrange‐point positions at time t
    lag_points = find_all_lagrange_points(planet_name, planets, t)
    point_keys = list(lag_points.keys())           # ['L1','L2','L3','L4','L5']
    planet_names = [p.name for p in planets]

    # initialize result structure
    result = {
        pt: {name: 0.0 for name in planet_names}
        for pt in point_keys
    }

    for p in planets:
        # high‐precision mass
        m = mp.mpf(p.mass)

        # planet’s instantaneous position and momentum
        pos_p = np.array(p.position(t))
        p_vec = np.array(p.momentum(t), dtype=float)
        p_sq  = mp.mpf(np.dot(p_vec, p_vec))
        v_sq  = p_sq / (m**2)

        for pt in point_keys:
            # distance from planet to this Lagrange point
            lp_pos = np.array(lag_points[pt], dtype=float)
            d = np.linalg.norm(lp_pos - pos_p)
            if d == 0:
                continue

            # gravitational slowdown
            td_grav = compute_gravitation_slowdown(
                G_mp=G_mp,
                m=m,
                d=d,
                c_mp=c_mp
            )

            # kinematic slowdown
            td_kin = compute_kinematic_slowdown(
                v_sq=v_sq,
                c_mp=c_mp
            )

            # store (convert to attoseconds)
            result[pt][p.name] = float((td_grav + td_kin) * afs_conversion)

    return result

def generate_simulation(planets, planet_name, t=0):
    """
    Generate a snapshot log of all planets and the five Lagrange‐point satellites
    for the given planet at time t, including positions and time‐dilation info.

    Parameters
    ----------
    planets : list of Planet
        A list of Planet objects representing the bodies (Sun, Mercury, Venus, etc.).
    planet_name : str
        Name of the planet whose Lagrange points we want (e.g. 'Jupiter').
    t : float, optional
        Simulation time in seconds at which to evaluate each body's state. Default is 0.

    Returns
    -------
    dict
        A nested dictionary with the following structure:

        {
            "time": float,  # equals t
            "planets": { ... },
            "satellites": {
                "L1": {
                    "position": { "px": float, "py": float, "pz": float },
                    "time_dilation_per_planet": { "<planet_name>": float, ... },
                    "time_dilation_total": float
                },
                ...
                "L5": { ... }
            }
        }
    """
    log = {
        'time': t,
        'planets': {},
        'satellites': {}
    }

    # --- record all planet states ---
    for p in planets:
        px, py, pz = p.position(t)
        vx, vy, vz = p.velocity(t)
        mx, my, mz = p.momentum(t)
        log['planets'][p.name] = {
            'position': {'px': px, 'py': py, 'pz': pz},
            'velocity': {'vx': vx, 'vy': vy, 'vz': vz},
            'momentum': {'mx': mx, 'my': my, 'mz': mz}
        }

    # --- compute L1–L5 positions ---
    lag_points = find_all_lagrange_points(planet_name, planets, t)
    # --- compute time dilations at those points ---
    td_per_planet = calculate_time_dilation_per_planet(planets, planet_name, t)

    # --- assemble satellites entries ---
    for Lname, pos in lag_points.items():
        px, py, pz = pos
        td_dict = td_per_planet[Lname]
        total_td = sum(td_dict.values())

        log['satellites'][Lname] = {
            'position': {'px': px, 'py': py, 'pz': pz},
            'time_dilation_per_planet': td_dict,
            'time_dilation_total': total_td
        }

    return log