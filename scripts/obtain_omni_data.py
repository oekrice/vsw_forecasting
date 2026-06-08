#This is a testbed for making the forecast scripts actually nice, and doing it all properly and things. HA, that went well!

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import cma
import csv
from datetime import datetime

this_directory = os.getcwd() + "/"
sys.path.append(this_directory+'/viz/HUXt-master/code')

import huxt_inputs as Hin
import huxt as H
import huxt_analysis as HA

import wind_forecast as fcast  #This should now contain everything we need...

#Just a script to download OMNI measurements and save them in a reasonable format.
#Will eventually look into removing CMEs and other snazzy things.
#But this will require some nuance.

print('Downloading OMNI data')
path = os.getcwd()

dtime_min = datetime(2005, 1, 1)
dtime_max = datetime(2026, 1, 1)

data_omni = Hin.get_omni(dtime_min, dtime_max)
all_dtime_omni = data_omni['datetime']
all_vsw_omni = data_omni['V'].values

omni_data = []

for i in range(len(all_dtime_omni)):
    omni_data.append([all_dtime_omni[i], all_vsw_omni[i]])

omni_fname = './data/shared_data/omni.csv'
with open(omni_fname, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(omni_data)

# dtime_omni, vsw_omni = fcast.get_average_speeds(all_dtime_omni, all_vsw_omni, spinup_time = 0.0, cadence=24, target_times=target_times, verbose=True)
#
# np.savetxt(('./data/' + run + '/times_omni.txt'), dtime_omni.astype("datetime64[s]").astype(str), fmt="%s")
# np.savetxt(('./data/' + run + '/vs_omni.txt'), vsw_omni)


