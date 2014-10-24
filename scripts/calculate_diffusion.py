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

def plot_fast(energy, yaxis, ylabel):
    total_plt = plt.figure()
    total_ax = plt.gca()
    total_ax.set_yscale('log')
    total_ax.set_xscale('log')
    total_ax.set_ylabel(ylabel)
    total_ax.set_xlabel('Energy [MeV]')
    total_ax.plot(energy, yaxis)
    plt.xlim([0.625e-6, 20.0])

def readin(filename):
    with open(filename, 'r') as fh:
        lines = fh.read().splitlines()
    data = np.array([float(x) for x in lines])
    data = data[1:] # remove most thermal group
    data = data[::-1] # reverse area for easy computation later
    return data

# True answer
energy_ratio = readin("energy_ratio.dat")
ratio = readin("ratio.dat")

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

# Get p1 scattering rate
p1_scatt_rate = readin("p1_scattering.out")

# Get cumulative outscatter rate
outscatt_rate = readin("outscatterc.out")

# Get within group scattering rate
winscatt_rate = readin("withinscatterc.out")

# Calculate xs
xs_a = abs_rate / flux
xs_s = scatt_rate / flux
p1_xs_s = p1_scatt_rate / flux
flux = flux/flux.sum()
xs_t = xs_a + xs_s

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

# Calculate transport xs based on diffusion coeff
xs_tr = 1.0/(3.0*diff) 

# Calculate transport xs based on outscatter approx
xs_tr_out = xs_t - p1_xs_s

# Plots
plot(energy, flux, "Flux")
plot(energy, xs_t, "Total XS")
plot(energy, cflux, "Cumulative Flux")
plot(energy, xsc_a, "Cumuative Absorption XS")
plot(energy, xsc_s, "Cumulative Scattering XS")
plot(energy, probc, "Outscattering probability")
plot(energy, xsc_r, "Cumulative Removal XS")
plot(energy, mig_area, "Cumulative Migration Area")
plot(energy, diffc, "Cumulative Diff Coef")
plot(energy, p1_xs_s / xs_s, "Mu_bar")

# Custom plot for diffusion coefficient
plot_fast(energy, diff, "Diff Coef")
plot_fast(energy, xs_tr, "Transport XS")
ratio_plt = plt.figure()
ratio_ax = plt.gca()
ratio_ax.set_xscale('log')
ratio_ax.set_ylabel("Transport-to-Total XS")
ratio_ax.set_xlabel('Energy [MeV]')
ratio_ax.plot(energy_ratio, ratio, 'r-', label="Exact from P1 Theory", linewidth=3.0)
ratio_ax.plot(energy, xs_tr/xs_t, 'b--', label="New tally", linewidth=3.0) 
ratio_ax.plot(energy, xs_tr_out/xs_t, 'g.', label="Outscatter Approx.", linewidth=3.0)
legend = ratio_ax.legend(loc='upper center', shadow=True)
frame = legend.get_frame()
frame.set_facecolor('0.90')
plt.xlim([0.625e-6, 20.0])

plt.show()
