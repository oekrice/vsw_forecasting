"""
    VIZ/MAPS.PY
    Script to compare predicted lat-long maps of magnetic field from different dumfric-fcast runs.
    
    A Yeates 2025-Apr
"""
import os
import numpy as np
from scipy.io import netcdf_file, FortranFile
import sys
if (len(sys.argv) > 1):    # to avoid needing xwindows (allows ssh remote use)
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt
import datetime
sys.path.append("./viz/tools")
import utils
import wind

path = '/Volumes/bmjg46/prs/dumfric/'
# runs = ['fcast-repl2_s20250401.12', 'fcast-over2_s20250401.12']  
# runs = ['fcast-repl2_u20250411.12', 'fcast-over2_u20250411.12']  
# runs = ['fcast-repl3_s20250401.12', 'fcast-repl3_u20250409.12', 'fcast-repl3_u20250411.12']  
# runs = ['fcast-repl5_s20250428.12']
# runs=['fcast-repl4_s20250301.12']
runs=['fcast-repl6_s20250301.12']
colors = ['tab:blue', 'tab:brown', 'tab:orange']

# Location of heliospheric boundary [outer boundary of Schatten extension, in R_sun]:
r_hb = 21.5

# Whether to re-trace field line maps/coronal holes:
retrace = True

# For purposes of labelling:
runs1 = [run.replace('_', '\_') for run in runs]
nruns = len(runs)

# GET FINAL SNAPSHOT FOR EACH RUN
# --------------------------------
snaps = []
for i in range(nruns):
    snaps1 = utils.list_snaps(path+runs[i]+"/", prefix="b")
    snaps.append(snaps1[-1])
    print(path+runs[i], snaps[i])

# POSTPROCESSING TO COMPUTE SCHATTEN EXTENSION AND FIELD LINE MAPPINGS
# --------------------------------------------------------------------
if retrace:
    for i in range(nruns):
        wind.windbnd(snaps[i], r_hb, path=path+runs[i]+"/")


# PLOTS
# -----
plt.figure(figsize=(20,nruns*3))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

for i in range(nruns):

    # PHOTOSPHERIC FIELD
    # ------------------
    fh = netcdf_file(path+runs[i]+"/"+snaps[i], "r", mmap=False)
    r = fh.variables["r"][:]
    th = fh.variables["th"][:]
    ph = fh.variables["ph"][:]
    br0 = fh.variables["br"][1:-1,1:-1,0]
    fh.close()
    plt.subplot(nruns, 4, 4*i+1)
    plt.pcolormesh(np.rad2deg(ph), np.cos(th), br0.T, vmin=-50, vmax=50, cmap='bwr')
    plt.title(runs[i])
    plt.xlabel('Carrington longitude')
    plt.ylabel('Sine latitude')

    # CORONAL HOLE MAP
    # ----------------
    fid = FortranFile(path+runs[i]+"/chmap_"+snaps[i]+'.unf', 'r')
    s_chmap = fid.read_reals(dtype=np.float64)
    ph_chmap = fid.read_reals(dtype=np.float64)
    chmap = fid.read_ints(dtype=np.int32).reshape((np.size(ph_chmap,0)-1, np.size(s_chmap,0)-1))
    chmap = np.swapaxes(chmap, 0, 1)
    fid.close()
    plt.subplot(nruns, 4, 4*i+2)
    plt.pcolormesh(np.rad2deg(ph_chmap), s_chmap, chmap, vmin=-1, vmax=1, cmap='bwr')
    plt.pcolormesh(np.rad2deg(ph), np.cos(th), br0.T, vmin=-50, vmax=50, cmap='bwr', alpha=0.2)
    plt.xlabel('Carrington longitude')
    # plt.ylabel('Sine latitude')
    plt.title('Coronal holes (phot bdry)')

    # HELIOSPHERIC BOUNDARY
    # ---------------------
    fh = netcdf_file(path+runs[i]+"/windbound_"+snaps[i], "r", mmap=False)
    fs = fh.variables["expansionFactor"][:]
    chd = fh.variables["CHBDistance"][:]
    vr = fh.variables["vr"][:]
    br_bnd = fh.variables["br"][:]
    date = fh.date.decode("utf-8")
    fh.close()
    plt.subplot(nruns, 4, 4*i+3)
    plt.pcolormesh(np.rad2deg(ph), np.cos(th), br_bnd, vmin=-1e-3, vmax=1e-3, cmap='bwr')
    plt.xlabel('Carrington longitude')
    # plt.ylabel('Sine latitude')
    plt.title(r'$B_r$ (heliospheric bdry)')

    plt.subplot(nruns, 4, 4*i+4)
    plt.pcolormesh(np.rad2deg(ph), np.cos(th), vr, vmin=0, vmax=1000, cmap='cubehelix')
    plt.colorbar()
    plt.xlabel('Carrington longitude')
    # plt.ylabel('Sine latitude')
    plt.title(r'$v_r$ [km$\,$s$^{-1}$] (heliospheric bdry)')

plt.tight_layout()

# Output:
if (len(sys.argv) > 1):
    plt.savefig('maps.png', bbox_inches='tight')
else:
    plt.show()
