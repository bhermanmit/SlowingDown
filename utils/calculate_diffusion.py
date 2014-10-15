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

# Get flux
with open('flux.out', 'r') as fh:
    lines = fh.read().splitlines()
flux = np.array([float(x) for x in lines])

# Get removal rate
with open('removal_rate.out', 'r') as fh:
    lines = fh.read().splitlines()
removal_rate = np.array([float(x) for x in lines])

# Calculate cumulate migration area
mig_area = 1.0/6.0*mig_rate/kill_rate

# Calculate removal xs
removal_xs = removal_rate / flux

# Calculate cumulative diffusion coefficient
difftog = mig_area*removal_xs

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

flux_plt = plt.figure()
flux_ax = plt.gca()
flux_ax.set_yscale('log')
flux_ax.set_xscale('log')
flux_ax.set_ylabel('Flux')
flux_ax.set_xlabel('Energy [MeV]')
flux_ax.plot(energy, flux)

removal_xs_plt = plt.figure()
removal_xs_ax = plt.gca()
removal_xs_ax.set_yscale('log')
removal_xs_ax.set_xscale('log')
removal_xs_ax.set_ylabel('Removal XS')
removal_xs_ax.set_xlabel('Energy [MeV]')
removal_xs_ax.plot(energy, removal_xs)

difftog_plt = plt.figure()
difftog_ax = plt.gca()
difftog_ax.set_yscale('log')
difftog_ax.set_xscale('log')
difftog_ax.set_ylabel('Cumlative D to group')
difftog_ax.set_xlabel('Energy [MeV]')
difftog_ax.plot(energy, difftog)

plt.show()
