import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

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
dchb_array = loaded_ds['dchb_array'].to_numpy() # input param 1
fxp_array = loaded_ds['fxp_array'].to_numpy() # input param 2
v_21p5 = loaded_ds['v_21p5'].to_numpy() # target velocities <==== THESE ARE THE BACKMAPPED VELOCITY BOUNDARIES
sub_earth_longs = np.radians(loaded_ds['long'])
mjd_times = loaded_ds['mjd']

loaded_ds.close() # close file to keep file explorer happy

print(np.min(dchb_array), np.max(dchb_array))
print(np.min(fxp_array), np.max(fxp_array))
## Plot arrays to see that things are behaving
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12,6), sharey=True)#, dpi=300)

im0 = axes[0].imshow(dchb_array, extent=[sub_earth_longs[0], sub_earth_longs[-1], mjd_times[0], mjd_times[-1]], aspect='auto', cmap='inferno', origin='lower')
fig.colorbar(im0, ax = axes[0], label='DCHB [deg]', aspect=20, orientation='horizontal', location='top')

im1 = axes[1].imshow(fxp_array, extent=[sub_earth_longs[0], sub_earth_longs[-1], mjd_times[0], mjd_times[-1]], aspect='auto', cmap='inferno', vmax=1000, origin='lower')
fig.colorbar(im1, ax = axes[1], label='fxp', aspect=20, orientation='horizontal', location='top')

im1 = axes[2].imshow(v_21p5, extent=[sub_earth_longs[0], sub_earth_longs[-1], mjd_times[0], mjd_times[-1]], aspect='auto', cmap='inferno', vmax=600, origin='lower')
fig.colorbar(im1, ax = axes[2], label='Vsw @ 21.5rS [km/s]', aspect=20, orientation='horizontal', location='top')

axes[0].set_ylabel('MJD [Days]')
axes[0].set_xlabel('Carrington Lon. [Rad.]')
axes[1].set_xlabel('Carrington Lon. [Rad.]')
axes[2].set_xlabel('Carrington Lon. [Rad.]')

plt.savefig('./plots/boundary_comparison/reference.png')
