#!/usr/bin/env python3

import yaml
import numpy as np
import matplotlib.pyplot as plt
import argparse

# --------------------------
# Load YAML data
# --------------------------
def load_simulation_data(yaml_path):
    with open(yaml_path, 'r') as f:
        data = yaml.safe_load(f)
    return data['simulation_data']

# --------------------------
# Extract total time dilations for each satellite
# --------------------------
def extract_total_dilations(sim_data):
    satellites = ['L0', 'L1', 'L2', 'L3', 'L4', 'L5']
    logs = {sat: [] for sat in satellites}
    days = []

    for entry in sim_data:
        day = entry['day']
        days.append(day)
        for sat in satellites:
            td_total = entry['satellites'][sat]['time_dilation_total']
            logs[sat].append(td_total)

    return np.array(days), logs

# --------------------------
# Plot total time dilations
# --------------------------
def plot_total_dilations(days, logs):
    plt.figure(figsize=(12, 8))
    for sat, values in logs.items():
        plt.plot(days, np.array(values)/1e9, label=sat)  # convert to seconds x10⁹

    plt.xlabel("Day")
    plt.ylabel("Total Time Dilation (seconds × 10⁹)")
    plt.title("Total Time Dilation per Satellite")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

# --------------------------
# Main entrypoint
# --------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("yaml_file", help="Path to YAML simulation file")
    args = parser.parse_args()

    sim_data = load_simulation_data(args.yaml_file)
    days, logs = extract_total_dilations(sim_data)
    plot_total_dilations(days, logs)
