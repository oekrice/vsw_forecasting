"""
    VIZ/DIAGNOSTICS.PY
    Script to plot global diagnostics for one or more dumfric-fcast runs.

    Use this to monitor progress or compare runs.
    
    A Yeates 2025-Apr
"""
import numpy as np
from scipy.io import netcdf_file
import os
import sys
if (len(sys.argv) > 1):    # to avoid needing xwindows (allows ssh remote use)
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt
import datetime

# SELECT RUNS AND COLOURS FOR GRAPHS:
path = '/Volumes/bmjg46/prs/dumfric/'
# runs = ['fcast-repl2_s20250401.12', 'fcast-over2_s20250401.12', 'fcast-repl2_u20250409.12', 'fcast-over2_u20250409.12' ]  
# runs = ['fcast-repl2', 'fcast-over2']
# runs = ['fcast-repl4', 'fcast-repl6']
runs = ['fcast-repl7', 'fcast-repl8']
colors = ['tab:blue', 'tab:red']

# Read in diagnostic files:
t, dtime, pflux, oflux, energy, current, poynting, runtype, labels = [], [], [], [], [], [], [], [], []


pflux_max, oflux_max, energy_max, current_max, poynting_min = 0, 0, 0, 0, 0
for i in range(len(runs)):

    # List updates of this run:
    for file in os.listdir(path):
        if file.startswith(runs[i]+"_s"):
            srun = file
            break
    runs1 = []
    for file in os.listdir(path):
        if file.startswith(runs[i]+"_u"):
            runs1.append(file)
    runs1.sort()
    runs1 = [srun] + runs1[:-1] # leave off last run

    labels1 = [runs[i].replace('_', '\_')]

    t1, dtime1, pflux1, oflux1, energy1, current1, poynting1 = [], [], [], [], [], [], []

    for run in runs1:
        fh = netcdf_file(path+run+'/diagnostics.nc', 'r', mmap=False)
        t1.append(fh.variables['t'][:]/86400.0)
        # Convert time to datetime object:
        fid = open(path+run+'/prepared.nml')
        lines = fid.readlines()
        fid.close()
        starttime = datetime.datetime.strptime(lines[5].split(' ')[-1].split(',')[0], '%Y%m%d.%H')
        dtime1.append( [starttime + datetime.timedelta(days=t00) for t00 in t1[-1]] )
        pflux1.append(fh.variables['phot_flux'][:])
        oflux1.append(fh.variables['open_flux'][:])
        energy1.append(fh.variables['energy'][:])
        current1.append(fh.variables['mean_current'][:])
        poynting1.append(fh.variables['poynting_outer'][:])
        fh.close()

        pflux_max = max(pflux_max, max(pflux1[-1]))
        oflux_max = max(oflux_max, max(oflux1[-1]))
        energy_max = max(energy_max, max(energy1[-1]))
        current_max = max(current_max, max(current1[-1]))
        poynting_min = min(poynting_min, min(poynting1[-1]))

        labels1.append('')

    t.append(t1)
    dtime.append(dtime1)
    pflux.append(pflux1)
    oflux.append(oflux1)
    energy.append(energy1)
    current.append(current1)
    poynting.append(poynting1)
    labels.append(labels1)


plt.figure(figsize=(6,10))
plt.rc("text", usetex=True)
plt.rc("font", family="serif")

plt.subplot(511)
for k, run in enumerate(runs):
    for j in range(len(t[k])):
        plt.plot(dtime[k][j], pflux[k][j]/1e24, label=labels[k][j], color=colors[k])
        plt.plot(dtime[k][j][-1], pflux[k][j][-1]/1e24, '.', color=colors[k])
        plt.ylabel(r'$\Phi_{\rm phot}$ [$10^{24}$ Mx]',)
        plt.xlim(dtime[k][0][0], dtime[k][-1][-1] + datetime.timedelta(days=2))
        plt.ylim(0, 1.1*pflux_max/1e24)
        plt.legend()
        plt.gcf().autofmt_xdate()

plt.subplot(512)
for k, run in enumerate(runs):
    for j in range(len(t[k])):
        plt.plot(dtime[k][j], oflux[k][j]/1e23, color=colors[k])
        plt.plot(dtime[k][j][-1], oflux[k][j][-1]/1e23, '.', color=colors[k])
        plt.ylabel(r'$\Phi_{\rm open}$ [$10^{23}$ Mx]')
        plt.xlim(dtime[k][0][0], dtime[k][-1][-1] + datetime.timedelta(days=2))
        plt.ylim(0, 1.1*oflux_max/1e23)
        plt.gcf().autofmt_xdate()

plt.subplot(513)
for k, run in enumerate(runs):
    for j in range(len(t[k])):
        plt.plot(dtime[k][j], energy[k][j]/1e33, color=colors[k])
        plt.plot(dtime[k][j][-1], energy[k][j][-1]/1e33, '.', color=colors[k])
        plt.ylabel(r'$E_{\rm mag}$ [$10^{33}$ erg]')
        plt.xlim(dtime[k][0][0], dtime[k][-1][-1] + datetime.timedelta(days=2))
        plt.ylim(0, 1.1*energy_max/1e33)
        plt.gcf().autofmt_xdate()

plt.subplot(514)
for k, run in enumerate(runs):
    for j in range(len(t[k])):
        plt.plot(dtime[k][j], current[k][j]/1e-10, color=colors[k])
        plt.plot(dtime[k][j][-1], current[k][j][-1]/1e-10, '.', color=colors[k])
        plt.ylabel(r'$\langle|\nabla\times\textbf{B}|\rangle$ [$10^{-10}$ Gcm$^{-1}$]')
        plt.xlim(dtime[k][0][0], dtime[k][-1][-1] + datetime.timedelta(days=2))
        plt.ylim(0, 1.1*current_max/1e-10)
        plt.gcf().autofmt_xdate()

plt.subplot(515)
for k, run in enumerate(runs):
    for j in range(len(t[k])):
        plt.plot(dtime[k][j], poynting[k][j]/1e26, color=colors[k])
        plt.plot(dtime[k][j][-1], poynting[k][j][-1]/1e26, '.', color=colors[k])
        plt.ylabel(r'$S_{\rm outer}$ [$10^{26}$ erg$\,$s$^{-1}$]')
        plt.xlim(dtime[k][0][0], dtime[k][-1][-1] + datetime.timedelta(days=2))
        plt.ylim(1.1*poynting_min/1e26, -0.1*poynting_min/1e26)
        plt.gcf().autofmt_xdate()

plt.tight_layout()

# Output:
if (len(sys.argv) > 1):
    plt.savefig('diagnostics.png', bbox_inches='tight')
else:
    plt.show()
