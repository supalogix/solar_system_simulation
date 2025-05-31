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
import json

from models import Planet
from time_dilation import generate_simulation
from config import planets_data, planets_data_2, colors
from constants import solar_mass
from helpers import generate_lagrange_orbital_params

def generate_planets(planets_data):
    """
    Generate a list of Planet objects (including the Sun) from raw orbital data.

    This function takes a list of dictionaries describing each planet’s orbital elements
    and mass, computes the Sun’s barycentric semi-major axis based on those masses and
    semi-major axes, inserts a Sun entry at the beginning of the list, and then
    constructs a Planet instance for each entry.

    Parameters
    ----------
    planets_data : list of dict
        A list where each dict corresponds to a planet and must contain the keys:
            - "name"   (str): Identifier for the body (e.g., "Mercury", "Earth").
            - "type"   (str): Body type (e.g., "PLANET").
            - "a"      (float): Semi-major axis of the orbit in meters.
            - "e"      (float): Orbital eccentricity (0 <= e < 1).
            - "i"      (float): Inclination in degrees.
            - "Omega"  (float): Longitude of ascending node in degrees.
            - "omega"  (float): Argument of periapsis in degrees.
            - "nu"     (float): True anomaly at epoch in degrees.
            - "mass"   (float): Mass of the body in kilograms.

    Returns
    -------
    list of Planet
        A list of Planet instances, in which the first element is the Sun (with
        its semi-major axis set to the barycentric value computed from the input
        planets), followed by the planets in the same order as `planets_data`.
    """

    # Compute the Sun’s reflex semi‐major axis about the system barycenter
    # based on: Murray, C. D., & Dermott, S. F. (2000). Solar System Dynamics. 
    M_sun = 1.9885e30  # kg
    sun_a = sum(p["mass"] * p["a"] for p in planets_data) / M_sun

    # Insert the Sun at the start of the list:
    planets_data.insert(0, {
        "name": "Sun",
        "type": "PLANET",
        "a": sun_a,     
        "e": 0.0,       
        "i": 0.0,       # define the reference plane
        "Omega": 0.0,
        "omega": 0.0,
        "nu": 0.0,
        "mass": M_sun
    })

    planets = [Planet(**pdata) for pdata in planets_data]

    return planets

def to_serializable(val):
    # NumPy scalar?
    if isinstance(val, np.generic):
        return val.item()

    raise TypeError(f"Type {type(val)} not serializable")


def main():
    config = next(p for p in planets_data_2 if p['name'] == "Jupiter")

    planets = generate_planets(planets_data=planets_data_2)

    satellites_data = generate_lagrange_orbital_params(config, solar_mass.value)

    satellites = [Planet(**pdata) for pdata in satellites_data]

    days_in_year = 365 * 24
    times = np.arange(days_in_year)

    for day in times:
        t = day * 86400

        log = generate_simulation(
            planets=planets, 
            satellites=satellites,
            t=t
        )

        L1 = log['satellites']['L1']['time_dilation_total']
        L2 = log['satellites']['L2']['time_dilation_total']
        L3 = log['satellites']['L3']['time_dilation_total']
        L4 = log['satellites']['L4']['time_dilation_total']
        L5 = log['satellites']['L5']['time_dilation_total']

        log['day'] = day
        log['time_dilations'] = {}
        log['time_dilations']['L4'] = {}
        log['time_dilations']['L4']['L1'] = L1 - L4
        log['time_dilations']['L4']['L2'] = L2 - L4
        log['time_dilations']['L4']['L3'] = L3 - L4
        log['time_dilations']['L4']['L4'] = L4 - L4
        log['time_dilations']['L4']['L5'] = L5 - L4
        log['time_dilations']['L5'] = {}
        log['time_dilations']['L5']['L1'] = L1 - L5
        log['time_dilations']['L5']['L2'] = L2 - L5
        log['time_dilations']['L5']['L3'] = L3 - L5
        log['time_dilations']['L5']['L4'] = L4 - L5
        log['time_dilations']['L5']['L5'] = L5 - L5
        json_pretty = json.dumps(log, indent=2, default=to_serializable)
        print(json_pretty)

if __name__ == "__main__":
    main()
