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
from config import planets_data, colors
from plotting import plot_all_offsets, plot_individual, plot_satellite_contributions_over_orbit

def generate_positions(satellites, angle):
    # Compute each Lagrange point’s position by offsetting the main planet’s orbital angle:
    # - L1 and L2 lie directly along the Sun–Planet line (zero offset)
    # - L3 sits opposite the Planet (offset \pi)
    # - L4 and L5 form equilateral triangles with the Planet (offset by 120 degrees)
    positions = {
        key: sat.position(angle + offset)
        for (key, sat), offset in zip(
            satellites.items(), [0, 0, np.pi, 2*np.pi/3, -2*np.pi/3]
        )
    }

    return positions

def generate_time_dilations(planets, positions):
    # Compute the proper‐time slowdown at each Lagrange point by feeding its
    # current position into our time‐dilation function, so we can compare
    # how much each clock runs slow relative to our chosen reference.
    tds = {
        key: calculate_time_dilation(planets, position=pos)
        for key, pos in positions.items()
    }

    return tds

def generate_proper_time_offset(planets, satellites, main_planet, reference_satellite):
    # Generate an hourly time grid over one full Earth year (365 days × 24 hours),
    # giving us 8,760 evenly spaced time steps so we can track how the proper‐time
    # offsets evolve throughout the orbit.
    days_in_year = 365 * 24
    times = np.arange(days_in_year)

    # Prepare storage for proper-time offsets
    proper_time_offset_data = {key: [] for key in satellites}

    # We step through each simulated day to see how the Lagrange points move over time:
    #   - Converting days to seconds keeps our inputs consistent with the physics math.
    #   - Fetching the planets’s position defines the axis along which its Lagrange points lie.
    #   - Extracting the orbital angle pinpoints exactly where each Lagrange point sits for
    #     the proper-time calculation.
    for day in times:

        # Convert our day-based timestep into seconds, since all orbital and
        # time-dilation calculations expect time in SI units (seconds).
        t = day * 86400

        # Determine the main planet’s location in space at time t so we can
        # anchor the Lagrange points along the Sun–Planet line for this moment.
        main_planet_position = main_planet.position(t)

        # Compute Jupiter’s true bearing around the Sun (0 to 2*\pi), taking into account
        # the signs of both x and y. This ensures we place each Lagrange point at the
        # correct angular offset relative to Jupiter.
        angle = np.arctan2(main_planet_position[1], main_planet_position[0])

        # Compute each Lagrange point’s position by offsetting Jupiter’s orbital angle:
        # - L1 and L2 lie directly along the Sun–Jupiter line (zero offset)
        # - L3 sits opposite Jupiter (offset \pi)
        # - L4 and L5 form equilateral triangles with Jupiter (offset by 120 degrees)
        positions = generate_positions(satellites=satellites, angle=angle)

        # Compute the proper‐time slowdown at each Lagrange point by feeding its
        # current position into our time‐dilation function, so we can compare
        # how much each clock runs slow relative to our chosen reference.
        tds = generate_time_dilations(planets=planets, positions=positions)

        ref = tds[reference_satellite]

        # Record each point’s time‐dilation difference from our reference (L4)
        # so we track how much faster or slower each clock runs compared to L4 over time.
        for key, td in tds.items():
            proper_time_offset_data[key].append(td - ref)
        
    return times, proper_time_offset_data

def generate_planets(planets_data):
    planets = [Planet(**pdata) for pdata in planets_data]

    return planets

def generate_satellites(planets, main_planet_name):
    # Find Jupiter for reference by name
    main_planet = next(p for p in planets if p.name == main_planet_name)

    # Initialize satellites
    satellites = {
        "L1": Satellite(name="L1", distance_from_sun=main_planet.a * 0.99),
        "L2": Satellite(name="L2", distance_from_sun=main_planet.a * 1.01),
        "L3": Satellite(name="L3", distance_from_sun=main_planet.a),
        "L4": Satellite(name="L4", distance_from_sun=main_planet.a),
        "L5": Satellite(name="L5", distance_from_sun=main_planet.a),
    }

    return main_planet, satellites

def generate_plots(times, proper_time_offset_data, colors):
    # Plot combined and individual proper-time offsets
    plot_all_offsets(times=times, data=proper_time_offset_data, colors=colors)
    plot_individual(times=times, data=proper_time_offset_data, colors=colors)

def main():
    planets = generate_planets(
        planets_data=planets_data
    )

    main_planet, satellites = generate_satellites(
        planets=planets,
        main_planet_name="Jupiter"
    )

    times, proper_time_offset_data = generate_proper_time_offset(
        planets=planets, 
        satellites=satellites, 
        main_planet=main_planet,
        reference_satellite="L1"
    )

    generate_plots(
        times=times,
        proper_time_offset_data=proper_time_offset_data,
        colors=colors
    )

    # Plot satellite contributions over a full orbit for each satellite
    #for sat in satellites.values():
    #    plot_satellite_contributions_over_orbit(planets=planets, satellite=sat, steps=360)

if __name__ == "__main__":
    main()
