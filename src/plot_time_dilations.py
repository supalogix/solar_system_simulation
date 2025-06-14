#!/usr/bin/env python3

import yaml
import numpy as np
import matplotlib.pyplot as plt
import itertools
import argparse

# --------------------------
# Helper: Load YAML
# --------------------------
def load_simulation_data(yaml_path):
    with open(yaml_path, 'r') as f:
        data = yaml.safe_load(f)
    return data['simulation_data']

# --------------------------
# Helper: Extract Time Dilations
# --------------------------
def extract_pairwise_dilations(sim_data):
    # Generate a dictionary where each key is (anchor, target)
    pairwise_logs = {}

    for entry in sim_data:
        day = entry['day']
        td = entry['time_dilations']
        for anchor in td:
            for target in td[anchor]:
                pair = (anchor, target)
                if pair not in pairwise_logs:
                    pairwise_logs[pair] = []
                pairwise_logs[pair].append( (day, td[anchor][target]) )

    return pairwise_logs

# --------------------------
# Plotting function
# --------------------------
def plot_pairwise_dilations(pairwise_logs):
    # Sort pairs to have nice grouping
    pairs = sorted(pairwise_logs.keys())

    for anchor in ['L0','L1','L2','L3','L4','L5']:
        for target in ['L0','L1','L2','L3','L4','L5']:
            if anchor == target:
                continue
            pair = (anchor, target)
            if pair not in pairwise_logs:
                continue

            data = pairwise_logs[pair]
            days, dilations = zip(*data)

            plt.figure(figsize=(10,6))
            plt.plot(days, np.array(dilations)/1e9, label=f'{anchor} - {target}')
            plt.xlabel("Day")
            plt.ylabel("Time Dilation Difference (seconds × 10⁹)")
            plt.title(f"Time Dilation: {anchor} relative to {target}")
            plt.legend()
            plt.grid()
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
    pairwise_logs = extract_pairwise_dilations(sim_data)
    plot_pairwise_dilations(pairwise_logs)
