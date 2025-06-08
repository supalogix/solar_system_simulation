#!/usr/bin/env python3
"""
This script drives our end-to-end time dilation study for solar system Lagrange satellites,
and dumps both the configuration and the streaming simulation data to a timestamped YAML
file, with full support for numpy types via custom PyYAML representers.
"""

import os
import asyncio
import datetime
import random
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import yaml
import aiofiles

from models import Planet
from time_dilation import generate_simulation
from config import planets_data_2
from constants import solar_mass
from helpers import generate_lagrange_orbital_params

# ---- Register NumPy representers with PyYAML ----
# PyYAML cannot serialize NumPy scalars and arrays
# Teach SafeDumper to convert numpy types into Python floats/ints and lists
# Avoids RepresenterError during safe_dump calls
# Ensures valid, human-readable YAML output when streaming large numpy-based simulation data

def _np_scalar_representer(dumper, data):
    """
    Represent any numpy scalar (e.g. np.float64, np.int32) as a native YAML float/int.
    """

    # Use the float tag; PyYAML will coerce to int if no decimal point
    return dumper.represent_scalar(
        u'tag:yaml.org,2002:float',
        str(data.item())
    )

def _np_array_representer(dumper, data):
    """
    Represent numpy.ndarray as a YAML sequence (list).
    """

    return dumper.represent_sequence(
        u'tag:yaml.org,2002:seq',
        data.tolist()
    )

# Attach to SafeDumper so that yaml.safe_dump will use these
yaml.SafeDumper.add_multi_representer(np.generic, _np_scalar_representer)
yaml.SafeDumper.add_multi_representer(np.ndarray, _np_array_representer)


# ---- YAML INITIALIZER ----
# Prepare a fresh YAML file for streaming 
# Creates/truncates simulation_dump/file 
# Writes top-level 'configuration:' block 
# Starts empty 'simulation_data:' sequence
# Allows appending each day's log under it 
# Maintains clear, human-readable structure

async def init_yaml_file(path: str,
                         config_planets: list[dict],
                         config_sats:   list[dict]):
    """
    Create/truncate the YAML file at `path`, write:
      configuration:
        planets:   <config_planets>
        satellites:<config_sats>
      simulation_data:
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    async with aiofiles.open(path, 'w') as f:
        # Top-level "configuration:" key
        await f.write("configuration:\n")
        # Dump the planets + satellites under that, indented 2 spaces
        block = yaml.safe_dump(
            {'planets': config_planets, 'satellites': config_sats},
            sort_keys=False,
            default_flow_style=False
        )
        for line in block.splitlines():
            await f.write(f"  {line}\n")
        # Start the empty simulation_data list
        await f.write("simulation_data:\n")


# ---- YAML APPENDER ----
# Append logs under top-level 'simulation_data:'
# Prefix each entry with '- ' and indent nested
# Ensures valid YAML sequence formatting
# Streams entries without holding all in memory
# Maintains clear, human-readable output


async def append_log(path: str, log: dict):
    """
    Append one `log` dict under simulation_data as a new sequence item:
      simulation_data:
        - day: 0
          ...
        - day: 1
          ...
    """
    yaml_str = yaml.safe_dump(log, sort_keys=False, default_flow_style=False)
    lines    = yaml_str.splitlines()
    # Prefix first line with "  - ", then indent the rest 4 spaces
    chunk    = ["  - " + lines[0]] + [f"    {ln}" for ln in lines[1:]]
    text     = "\n".join(chunk) + "\n"

    async with aiofiles.open(path, 'a') as f:
        await f.write(text)


# ---- PLANET / SATELLITE UTILS ----
# Generate Planet instances from raw orbital data
# Compute Sun’s barycentric semi-major axis
# Insert Sun entry and preserve planet order
# Encapsulate setup for planets and satellites
# Keeps simulation initialization concise


def generate_planets(planets_data):
    """
    Insert the Sun (with barycentric semi-major axis) into planets_data and
    return a list of Planet instances — but first randomize each planet’s
    starting true anomaly so they all begin at a random point in their orbits.
    """

    # 1) Randomize starting angle for each planet in degrees
    for p in planets_data:
        p['nu'] = random.uniform(0.0, 360.0)

    # 2) Compute the Sun’s reflex semi-major axis about the system barycenter
    M_sun = solar_mass.value
    sun_a = sum(p["mass"] * p["a"] for p in planets_data) / M_sun

    # 3) Insert the Sun at the start (we keep its nu = 0)
    planets_data.insert(0, {
        "name":  "Sun",
        "type":  "PLANET",
        "a":      sun_a,
        "e":      0.0,
        "i":      0.0,
        "Omega":  0.0,
        "omega":  0.0,
        "nu":     0.0,
        "mass":   M_sun
    })

    # 4) Build Planet objects
    return [Planet(**pdata) for pdata in planets_data]



def simulate_day(day: int,
                 planets_raw: list[dict],
                 planet_name: str) -> dict:
    """
    Worker function that lives in a subprocess.
    Reconstructs Planet objects, runs generate_simulation for the five
    Lagrange points of `planet_name` at time t, then computes and returns the log.
    """
    # Rebuild fresh Planet instances in this process
    planets = generate_planets([dict(p) for p in planets_raw])

    t = day * 86400
    # New signature: just planets + planet_name + time
    log = generate_simulation(planets, planet_name, t)

    # Compute pairwise time‐dilation offsets relative to L4 and L5
    totals = {pt: log['satellites'][pt]['time_dilation_total']
              for pt in ("L1","L2","L3","L4","L5")}
    log['day'] = day
    log['time_dilations'] = {
        'L4': {pt: totals[pt] - totals['L4'] for pt in totals},
        'L5': {pt: totals[pt] - totals['L5'] for pt in totals},
    }

    return log


async def run_simulation():
    now      = datetime.datetime.now()
    stamp    = now.strftime("%Y%m%d_%H%M%S")
    out_path = os.path.join("simulation_dump",
                            f"simulation_{stamp}.yaml")

    # Prepare raw planet configs
    planets_raw = [dict(p) for p in planets_data_2]
    # Identify the planet whose Lagrange points we want
    jup_cfg     = next(p for p in planets_data_2 if p['name']=="Jupiter")
    planet_name = jup_cfg['name']

    # Initialize YAML with configuration
    await init_yaml_file(out_path, planets_raw, None)

    #days = 365 * 24
    days = 1
    loop = asyncio.get_running_loop()

    with ProcessPoolExecutor() as pool:
        tasks = [
            loop.run_in_executor(
                pool,
                simulate_day,
                day,
                planets_raw,
                planet_name
            )
            for day in range(days)
        ]

        for completed in asyncio.as_completed(tasks):
            log = await completed
            await append_log(out_path, log)

    print(f"Wrote simulation to '{out_path}'")


if __name__ == "__main__":
    asyncio.run(run_simulation())