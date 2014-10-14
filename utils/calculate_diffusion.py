#!/bin/sh env python

import matplotlib.pyplot as plt
import numpy as np

# Get energy grid
with open('energy.out', 'r') as fh:
    lines = fh.read().splitlines()
energy = np.array([float(x) for x in lines])

# Get cumulative migration area rate
with open('migration_rate.out', 'r') as fh:
    lines = fh.read().splitlines()
mig_rate = np.array([float(x) for x in lines])

# Get cumulate kill rate
with open('kill_rate.out', 'r') as fh:
    lines = fh.read().splitlines()
kill_rate = np.array([float(x) for x in lines])

# Calculate cumulate migration area
mig_area = 1.0/6.0*mig_rate/kill_rate

# Make plots
mig_rate_plt = plt.figure()
mig_rate_ax = plt.gca()
mig_rate_ax.set_yscale('log')
mig_rate_ax.set_xscale('log')
mig_rate_ax.set_ylabel('Cumulative migration rate')
mig_rate_ax.set_xlabel('Energy [MeV]')
mig_rate_ax.plot(energy, mig_rate)

kill_rate_plt = plt.figure()
kill_rate_ax = plt.gca()
kill_rate_ax.set_yscale('log')
kill_rate_ax.set_xscale('log')
kill_rate_ax.set_ylabel('Cumulative kill rate')
kill_rate_ax.set_xlabel('Energy [MeV]')
kill_rate_ax.plot(energy, kill_rate)

mig_rate_plt = plt.figure()
mig_rate_ax = plt.gca()
mig_rate_ax.set_yscale('log')
mig_rate_ax.set_xscale('log')
mig_rate_ax.set_ylabel('Cumulative migration area')
mig_rate_ax.set_xlabel('Energy [MeV]')
mig_rate_ax.plot(energy, mig_area)
plt.show()
