"""
    VIZ/HUXT-FCAST.PY
    Script to run HUXt for one or more dumfric-fcast models, and plot predicted vsw vs time.
    
    [need to run maps.py first to generate wind boundary file]

    A Yeates 2025-Apr
"""
import os
import numpy as np
from scipy.io import netcdf_file, FortranFile
import astropy.units as u
from sunpy.coordinates.sun import B0
import sys
if (len(sys.argv) > 1):    # to avoid needing xwindows (allows ssh remote use)
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt
import datetime
sys.path.append("./viz/tools")
import utils
import wind
sys.path.append('./viz/HUXt-master/code')
import huxt as H
import huxt_inputs as Hin
import huxt_analysis as HA

path = '/Volumes/bmjg46/prs/dumfric/'
# runs = ['fcast-repl2_s20250401.12', 'fcast-over2_s20250401.12']  
runs = ['fcast-repl2_s20250401.12', 'fcast-over2_u20250409.12', 'fcast-over2_u20250411.12']  
# colors = ['tab:blue', 'tab:brown']
colors = ['tab:blue', 'tab:brown', 'tab:red']

labels = [run.replace('_', '\_') for run in runs]
nruns = len(runs)

# GET FINAL SNAPSHOT FOR EACH RUN
# --------------------------------
dtime_huxt = []
vsw_huxt = []
bpol_huxt = []
for i in range(nruns):
    wfiles1 = utils.list_snaps(path+runs[i]+"/", prefix="windbound_b")
    wfile = wfiles1[-1]

    # Get radius of heliospheric boundary:
    fh = netcdf_file(path+runs[i]+'/schat_'+wfile[10:], "r", mmap=False)
    r = fh.variables["r"][:]
    fh.close()
    r_hb = r[-1]
    
    # Read model output on heliospheric boundary:
    fh = netcdf_file(path+runs[i]+"/"+wfile, "r", mmap=False)
    fs = fh.variables["expansionFactor"][:]
    chd = fh.variables["CHBDistance"][:]
    vr_bnd = fh.variables["vr"][:]
    br_bnd = fh.variables["br"][:]
    date = fh.date.decode("utf-8")
    fh.close()
    
    dtime_snap = datetime.datetime.strptime(date, "%Y-%m-%d %H:%M")
    simtime = 27.27*u.day

    # Latitude of Earth at simulation time:
    _, vr_longs, vr_lats, br_map, br_longs, br_lats = Hin.get_PFSS_maps(path+runs[i]+"/"+wfile)
    E_lat = np.deg2rad(B0(dtime_snap))
    iE_lat = np.argmin(abs(vr_lats - E_lat))
    
    # Run HUXt:
    fcast_length = 7  # in days
    t_init = dtime_snap
    cr, cr_lon_init = Hin.datetime2huxtinputs(t_init)
    vr_long = vr_bnd[iE_lat,:]*u.km/u.s
    br_long = br_bnd[iE_lat,:]
    model = H.HUXt(v_boundary=vr_long, b_boundary=br_long, cr_num=cr, cr_lon_init=cr_lon_init, latitude=E_lat, lon_out=0.0*u.deg, simtime=fcast_length*u.day, dt_scale=4, frame = 'synodic', r_min = r_hb*u.solRad)
    model.solve([]) 

    # Extract time series at Earth:
    earth_series = HA.get_observer_timeseries(model, observer = 'Earth')
    dtime_huxt.append( earth_series['time'] )
    # mjd_huxt = earth_series['mjd']
    vsw_huxt.append( earth_series['vsw'] )
    bpol_huxt.append( earth_series['bpol'] )

    plt.plot(dtime_huxt[-1], vsw_huxt[-1], color=colors[i], label=labels[i])

plt.legend()
plt.show()