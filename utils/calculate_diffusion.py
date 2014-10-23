#!/bin/sh env python

import matplotlib.pyplot as plt
import numpy as np

def plot(energy, yaxis, ylabel):
    total_plt = plt.figure()
    total_ax = plt.gca()
    total_ax.set_yscale('log')
    total_ax.set_xscale('log')
    total_ax.set_ylabel(ylabel)
    total_ax.set_xlabel('Energy [MeV]')
    total_ax.plot(energy, yaxis)

def readin(filename):
    with open(filename, 'r') as fh:
        lines = fh.read().splitlines()
    data = np.array([float(x) for x in lines])
    data = data[1:] # remove most thermal group
    data = data[::-1] # reverse area for easy computation later
    return data

# Get energy grid
energy = readin("energy.out")
ng = energy.shape[0]

# Get cumulative r^2 rate
r2c_rate = readin("r2c.out")

# Get cumulate weight for r^2 rate  
wc_rate = readin("wc.out") 

# Get flux
flux = readin("flux.out") 

# Get scattering rate
scatt_rate = readin("scattering.out")

# Get absortpion rate
abs_rate = readin("absorption.out")

# Get cumulative outscatter rate
outscatt_rate = readin("outscatterc.out")

# Get within group scattering rate
winscatt_rate = readin("withinscatterc.out")

# Calculate xs
xs_a = abs_rate / flux
xs_s = scatt_rate / flux
flux = flux/flux.sum()

# Calculate cumulative flux, scattering and absoption
cflux = np.cumsum(flux)
cabs = np.cumsum(xs_a*flux) 
cscatt = np.cumsum(xs_s*flux)

# Calculate cumulative xs
xsc_a = cabs / cflux
xsc_s = cscatt / cflux
probc = outscatt_rate / (outscatt_rate + winscatt_rate)
xsc_r = xsc_a + xsc_s*probc

# Calculate cumulate migration area
mig_area = 1.0/6.0*r2c_rate/wc_rate

# Calculate diffusion coefficient
diffc = mig_area*xsc_r

# Calculate group diffusion coefficient
diff = np.zeros(diffc.shape)
diff[0]  = diffc[0]
for i in range(diff.shape[0]-1):
    ii = i + 1
    diffrate = 0.0
    fluxsum = flux[ii]
    for j in range(ii):
        diffrate += diff[j]*flux[j]
        fluxsum += flux[j]
        diff[ii] = (diffc[ii]*fluxsum - diffrate)/flux[ii]

# Plots
plot(energy, flux, "Flux")
plot(energy, cflux, "Cumulative Flux")
plot(energy, xsc_a, "Cumuative Absorption XS")
plot(energy, xsc_s, "Cumulative Scattering XS")
plot(energy, probc, "Outscattering probability")
plot(energy, xsc_r, "Cumulative Removal XS")
plot(energy, mig_area, "Cumulative Migration Area")
plot(energy, diffc, "Cumulative Diff Coef")
plt.show()
