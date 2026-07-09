import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from datetime import datetime
import csv
import wind_forecast as fcast  #This should now contain everything we need...

#%%
#==================================================
# reading in PFSS+SCS model boundary data and reconstructed observations at 21.5 solar radii
#==================================================

# read in hovmuller plot of pfss output along sub earth path across time
rhb=21.5 # outer boundary height [rS]
rss = 2.5 # source surface [rS]
hov_file = f'./data/shared_data/boundary_distribution.nc'
input_file = os.path.abspath(hov_file)
loaded_ds = xr.open_dataset(input_file)
loaded_ds = loaded_ds.sortby('mjd')

# Assign fields to arrays and coordinates
raw_surface_velocities = loaded_ds['v_21p5'].to_numpy() # target velocities <==== THESE ARE THE BACKMAPPED VELOCITY BOUNDARIES
observation_times = loaded_ds['mjd']

loaded_ds.close() # close file to keep file explorer happy

#Convert to datetime objects to make comparisons easy

raw_boundary_times = Time(observation_times, format='mjd').to_datetime()

#Compare the two datasets, checking for overlaps in log times
sources = ["p2g"]#, "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]
target_nlon = np.shape(raw_surface_velocities)[1]
source_coordinates = np.linspace(0.0,2*np.pi,361)
target_coordinates = np.linspace(0.0,2*np.pi,target_nlon)

for si, source in enumerate(sources):
    batch_id = si
    if (batch_id//2)%2 == 0:
        is_pfss = True
        model = "pfss"
    else:
        is_pfss = False
        model = "outflow"
    if (batch_id%2) == 0:
        rss = 2.5
    else:
        rss = 5.0
    if (batch_id//4) == 0:
        data_source = "gong"
    else:
        data_source = "hmi"

    nicetitle = f"{model}, rss = {rss}, source = {data_source}"

    #Do some comparisons here.
    #Establish the permissible ranges, and remove cme times or thereabouts (using established code)
    log_fname = f'./data/{source}/directory.csv'
    source_times = []; source_ids = []
    with open(log_fname, "r", encoding="utf-8") as f:
        data = csv.reader(f)
        for row in data:
            if len(row) > 3:
                source_times.append(datetime.strptime(row[2], "%Y-%m-%d %H:%M:%S"))
                source_ids.append(int(row[1]))

    source_times = np.array(source_times); source_ids = np.array(source_ids)
    #Filter to make the two datasets line up. Hopefully everything is a daily cadence.
    min_time = max(sorted(raw_boundary_times)[0], sorted(source_times)[0])
    max_time = min(sorted(raw_boundary_times)[-1], sorted(source_times)[-1])

    surface_velocities = raw_surface_velocities[(raw_boundary_times >= min_time) & (raw_boundary_times <= max_time),:]
    boundary_times = raw_boundary_times[(raw_boundary_times >= min_time) & (raw_boundary_times <= max_time)]
    source_times = source_times[(source_times >= min_time) & (source_times <= max_time)]

    print(np.shape(surface_velocities))


    #Run through source times and check they exist in boundary times. If not, don't include them
    source_ids = []
    scount = 0
    for bi in range(len(boundary_times)):
        found = False
        while not found:
            if source_times[scount] == boundary_times[bi]:
                source_ids.append(scount)
                scount += 1
                found = True
            else:
                scount += 1

    source_ids = np.array(source_ids)   #These are the snaps from the source file which need to be compared against. Good good.
    chd_data = np.zeros((len(source_ids),target_nlon))
    fs_data = np.zeros((len(source_ids),target_nlon))

    for i, source_id in enumerate(source_ids[:]):
        print(i/len(source_ids))
        s0, ph0, br0, fs, chd = fcast.data_functions.load_chb_distances(source ,source_id)
        fs = np.mean(fs, axis=0)
        fs = np.interp(target_coordinates,source_coordinates,fs)

        chd = np.mean(chd, axis=0)
        chd = np.interp(target_coordinates,source_coordinates,chd)

        chd_data[i] = chd
        fs_data[i] = fs


        #cmap = fcast.model_functions.compute_vr(source_id, source, method="wsa", params = None, doplot=False, iteration=0, huxt_name="test1", output_cmaps=False)
        #Sum this in the relevant direction and interpolate onto the correct grid size
        # cmap = np.mean(cmap, axis=0)
        # cmap = np.interp(target_coordinates,source_coordinates,cmap)
        # hovmoller_data[i] = cmap
        #

    print(np.shape(surface_velocities))
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12,6), sharey=True)#, dpi=300)

    im0 = axes[0].imshow(chd_data, extent=[0.0, 2*np.pi, np.min(source_times), np.max(source_times)], aspect='auto', cmap='inferno', origin='lower', vmin=0, vmax=12)
    fig.colorbar(im0, ax = axes[0], label='DCHB [deg]', aspect=20, orientation='horizontal', location='top')

    im1 = axes[1].imshow(fs_data, extent=[0.0, 2*np.pi, np.min(source_times), np.max(source_times)], aspect='auto', cmap='inferno', vmax=1000, origin='lower')
    fig.colorbar(im1, ax = axes[1], label='fxp', aspect=20, orientation='horizontal', location='top')

    im1 = axes[2].imshow(surface_velocities, extent=[0.0, 2*np.pi, np.min(source_times), np.max(source_times)], aspect='auto', cmap='inferno', vmax=600, origin='lower')
    fig.colorbar(im1, ax = axes[2], label='Vsw @ 21.5rS [km/s]', aspect=20, orientation='horizontal', location='top')

    axes[0].set_ylabel('MJD [Days]')
    axes[0].set_xlabel('Carrington Lon. [Rad.]')
    axes[1].set_xlabel('Carrington Lon. [Rad.]')
    axes[2].set_xlabel('Carrington Lon. [Rad.]')

    plt.savefig('./plots/boundary_comparison/calculated_reference.png')
    plt.show()

# ## Plot arrays to see that things are behaving
# fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12,6), sharey=True)#, dpi=300)
#
# im0 = axes[0].imshow(dchb_array, extent=[sub_earth_longs[0], sub_earth_longs[-1], mjd_times[0], mjd_times[-1]], aspect='auto', cmap='inferno', origin='lower')
# fig.colorbar(im0, ax = axes[0], label='DCHB [deg]', aspect=20, orientation='horizontal', location='top')
#
# im1 = axes[1].imshow(fxp_array, extent=[sub_earth_longs[0], sub_earth_longs[-1], mjd_times[0], mjd_times[-1]], aspect='auto', cmap='inferno', vmax=1000, origin='lower')
# fig.colorbar(im1, ax = axes[1], label='fxp', aspect=20, orientation='horizontal', location='top')
#
# im1 = axes[2].imshow(v_21p5, extent=[sub_earth_longs[0], sub_earth_longs[-1], mjd_times[0], mjd_times[-1]], aspect='auto', cmap='inferno', vmax=600, origin='lower')
# fig.colorbar(im1, ax = axes[2], label='Vsw @ 21.5rS [km/s]', aspect=20, orientation='horizontal', location='top')
#
# axes[0].set_ylabel('MJD [Days]')
# axes[0].set_xlabel('Carrington Lon. [Rad.]')
# axes[1].set_xlabel('Carrington Lon. [Rad.]')
# axes[2].set_xlabel('Carrington Lon. [Rad.]')
#
# plt.show()
