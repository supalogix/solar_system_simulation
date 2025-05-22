#!/usr/bin/env python3
"""
This script drives our end-to-end time dilation study for solar system Lagrange satellites.

We define each planet’s orbit around the Sun (actually computed from the solar system’s
barycenter) and then place “satellites” at the five Jupiter Lagrange points.  For each day
of a year, we compute how much proper time accumulates at each point—and then compare
everything to L4.  Choosing L4 as our zero point is purely arbitrary; we could just as
easily have used L1, L2, L3, or L5.  In principle, the most natural reference would be the
system barycenter, since that’s the frame in which we calculate the orbits.  But in reality
no probe ever sits at the barycenter: our only real measurements come from points we can
actually place satellites.  By anchoring to a true Lagrange-point location, we ground
our offsets in data that could be gathered by an actual mission.
"""

import numpy as np

from models import Planet, Satellite
from time_dilation import calculate_time_dilation
from plotting import plot_all_offsets, plot_individual, plot_satellite_contributions_over_orbit

def main():

    # Rather than using abstract units, we anchor each planet’s orbital elements and mass to Earth’s 
    # baseline so it’s easy to see how each world differs:
    #  - semi-major axis (a) is given in meters but described as a multiple of 1 AU (Earth’s orbit)
    #  - eccentricity (e) measures how stretched the orbit is (Earth’s e≈0.0167 is nearly circular)
    #  - inclination (i) is the tilt of the orbital plane relative to Earth’s (i=0 by definition)
    #  - mass is quoted in kilograms but noted as a multiple of Earth’s mass
    # Choosing Earth as the reference makes these comparisons intuitive: we know what 1 AU,
    # e≈0.017, i=0° and 1 Earth-mass feel like, so we can immediately grasp how Mercury’s
    # 0.39 AU or Jupiter’s 318× Earth-mass stand out.
    planets = [
        Planet(
            name="Mercury",
            a=5.791e10,    # semi-major axis: 0.39 AU (~0.39 x Earth's 1 AU)
            e=0.2056,      # eccentricity: 0.2056 (~12.3 x Earth's e=0.0167)
            i=7.0,         # inclination: 7.0 deg vs Earth's 0 deg
            Omega=48.3,    # longitude of ascending node: 48.3 deg vs 0 deg
            omega=29.1,    # argument of periapsis: 29.1 deg vs Earth's 102.9 deg
            nu=0,          # true anomaly at epoch: 0 deg
            mass=3.285e23  # mass: 0.055 x Earth mass
        ),

        Planet(
            name="Venus",
            a=1.082e11,    # semi-major axis: 0.72 AU (~0.72 x Earth's 1 AU)
            e=0.0067,      # eccentricity: 0.0067 (~0.4 x Earth's e)
            i=3.4,         # inclination: 3.4 deg vs Earth's 0 deg
            Omega=76.7,    # longitude of ascending node: 76.7 deg vs 0 deg
            omega=54.9,    # argument of periapsis: 54.9 deg vs Earth's 102.9 deg
            nu=0,          # true anomaly at epoch: 0 deg
            mass=4.867e24  # mass: 0.82 x Earth mass
        ),

        Planet(
            name="Earth",
            a=1.496e11,    # semi-major axis: 1.00 AU (baseline)
            e=0.0167,      # eccentricity: 0.0167 (baseline)
            i=0.0,         # inclination: 0 deg (baseline)
            Omega=0.0,     # longitude of ascending node: 0 deg (baseline)
            omega=102.9,   # argument of periapsis: 102.9 deg
            nu=0,          # true anomaly at epoch: 0 deg
            mass=5.972e24  # mass: 1.00 x Earth mass (baseline)
        ),

        Planet(
            name="Mars",
            a=2.279e11,    # semi-major axis: 1.52 AU (~1.52 x Earth's 1 AU)
            e=0.0935,      # eccentricity: 0.0935 (~5.6 x Earth's e)
            i=1.9,         # inclination: 1.9 deg vs 0 deg
            Omega=49.6,    # longitude of ascending node: 49.6 deg vs 0 deg
            omega=286.5,   # argument of periapsis: 286.5 deg
            nu=0,          # true anomaly at epoch: 0 deg
            mass=6.417e23  # mass: 0.11 x Earth mass
        ),

        Planet(
            name="Jupiter",
            a=7.785e11,    # semi-major axis: 5.20 AU (~5.20 x Earth's 1 AU)
            e=0.0484,      # eccentricity: 0.0484 (~2.9 x Earth's e)
            i=1.3,         # inclination: 1.3 deg vs 0 deg
            Omega=100.5,   # longitude of ascending node: 100.5 deg vs 0 deg
            omega=273.6,   # argument of periapsis: 273.6 deg
            nu=0,          # true anomaly at epoch: 0 deg
            mass=1.898e27  # mass: 318 x Earth mass
        ),

        Planet(
            name="Saturn",
            a=1.429e12,    # semi-major axis: 9.58 AU (~9.58 x Earth's 1 AU)
            e=0.0565,      # eccentricity: 0.0565 (~3.4 x Earth's e)
            i=2.5,         # inclination: 2.5 deg vs 0 deg
            Omega=113.6,   # longitude of ascending node: 113.6 deg vs 0 deg
            omega=339.6,   # argument of periapsis: 339.6 deg
            nu=0,          # true anomaly at epoch: 0 deg
            mass=5.683e26  # mass: 95 x Earth mass
        ),

        Planet(
            name="Uranus",
            a=2.871e12,    # semi-major axis: 19.19 AU (~19.19 x Earth's 1 AU)
            e=0.0463,      # eccentricity: 0.0463 (~2.8 x Earth's e)
            i=0.8,         # inclination: 0.8 deg vs 0 deg
            Omega=74.0,    # longitude of ascending node: 74.0 deg vs 0 deg
            omega=96.6,    # argument of periapsis: 96.6 deg
            nu=0,          # true anomaly at epoch: 0 deg
            mass=8.681e25  # mass: 14 x Earth mass
        ),

        Planet(
            name="Neptune",
            a=4.495e12,    # semi-major axis: 30.07 AU (~30.07 x Earth's 1 AU)
            e=0.0086,      # eccentricity: 0.0086 (~0.5 x Earth's e)
            i=1.8,         # inclination: 1.8 deg vs 0 deg
            Omega=131.8,   # longitude of ascending node: 131.8 deg vs 0 deg
            omega=276.5,   # argument of periapsis: 276.5 deg
            nu=0,          # true anomaly at epoch: 0 deg
            mass=1.024e26  # mass: 17 x Earth mass
        ),
    ]

    # Find Jupiter for reference by name
    jupiter = next(p for p in planets if p.name == "Jupiter")

    # Initialize satellites
    satellites = {
        "L1": Satellite(name="L1 Satellite", distance_from_sun=jupiter.a * 0.99),
        "L2": Satellite(name="L2 Satellite", distance_from_sun=jupiter.a * 1.01),
        "L3": Satellite(name="L3 Satellite", distance_from_sun=jupiter.a),
        "L4": Satellite(name="L4 Satellite", distance_from_sun=jupiter.a),
        "L5": Satellite(name="L5 Satellite", distance_from_sun=jupiter.a),
    }

    # Generate an hourly time grid over one full Earth year (365 days × 24 hours),
    # giving us 8,760 evenly spaced time steps so we can track how the proper‐time
    # offsets evolve throughout the orbit.
    days_in_year = 365 * 24
    times = np.arange(days_in_year)

    # Prepare storage for proper-time offsets
    colors = ["red", "green", "purple", "orange", "blue"]
    proper_time_offset_data = {key: [] for key in satellites}

    # We step through each simulated day to see how the Lagrange points move over time:
    #   - Converting days to seconds keeps our inputs consistent with the physics math.
    #   - Fetching Jupiter’s position defines the axis along which its Lagrange points lie.
    #   - Extracting the orbital angle pinpoints exactly where each Lagrange point sits for
    #     the proper-time calculation.
    for day in times:

        # Convert our day-based timestep into seconds, since all orbital and
        # time-dilation calculations expect time in SI units (seconds).
        t = day * 86400

        # Determine Jupiter’s location in space at time t so we can
        # anchor the Lagrange points along the Sun–Jupiter line for this moment.
        jpos = jupiter.position(t)

        # Compute Jupiter’s true bearing around the Sun (0 to 2*\pi), taking into account
        # the signs of both x and y. This ensures we place each Lagrange point at the
        # correct angular offset relative to Jupiter.
        angle = np.arctan2(jpos[1], jpos[0])

        # Compute each Lagrange point’s position by offsetting Jupiter’s orbital angle:
        # - L1 and L2 lie directly along the Sun–Jupiter line (zero offset)
        # - L3 sits opposite Jupiter (offset \pi)
        # - L4 and L5 form equilateral triangles with Jupiter (offset by 120 degrees)
        positions = {
            key: sat.position(angle + offset)
            for (key, sat), offset in zip(
                satellites.items(), [0, 0, np.pi, 2*np.pi/3, -2*np.pi/3]
            )
        }

        # Compute the proper‐time slowdown at each Lagrange point by feeding its
        # current position into our time‐dilation function, so we can compare
        # how much each clock runs slow relative to our chosen reference.
        tds = {
            key: calculate_time_dilation(planets=planets, position=pos)
            for key, pos in positions.items()
        }

        ref = tds["L4"]

        # Record each point’s time‐dilation difference from our reference (L4)
        # so we track how much faster or slower each clock runs compared to L4 over time.
        for key, td in tds.items():
            proper_time_offset_data[key].append(td - ref)


    # Plot combined and individual proper-time offsets
    plot_all_offsets(times=times, data=proper_time_offset_data, colors=colors)
    plot_individual(times=times, data=proper_time_offset_data, colors=colors)

    # Plot satellite contributions over a full orbit for each satellite
    for sat in satellites.values():
        plot_satellite_contributions_over_orbit(planets=planets, satellite=sat, steps=360)

if __name__ == "__main__":
    main()
