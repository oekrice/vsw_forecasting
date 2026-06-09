#This script is just for creating metadata for base data which were created before I had done the metadata stuff properly.
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import cma
import csv
from datetime import datetime, timedelta

import wind_forecast as fcast  #This should now contain everything we need...
from dtaidistance import dtw

start = datetime(2010, 1, 1) #This CAN'T change for a given run name. BE CAREFUL

obs_times = [start + timedelta(days=i) for i in range(5478)]

run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

for batch_id in range(7,8):

    #Get the model setup depending on the batch numbers
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
        source = "gong"
    else:
        source = "hmi"

    run_name = run_names[batch_id]


    run_parameters = {"observation_time": obs_times,
                    "base_name": run_names[batch_id],
                    "run_name": "test_run",
                    "model_type": model,
                    "calculate_base_model": True,
                    "overwrite_base_model": False,
                    "calculate_huxt": False,
                    "r_ss": rss,
                    "WSA_type": "standard",
                    "WSA_parameters": None,
                    "verbose": True,
                    "data_source": source,
                    "resolutions": [120,180,360],
                    "r_hb": 21.5,
                    "match_flag": False,
                    "velocity_type": "wsa",
                    "spinup_time": 5,
                    "forecast_length": 5,
                    "verbose": True,
                    "optimisation_type": "wasserstein",
                    "do_plots": True}

    snap_subset = np.arange(len(obs_times))
    for snap_id in snap_subset:
        print(batch_id, snap_id/len(snap_subset))
        obs_time = run_parameters["observation_time"][snap_id]
        resolutions =  [120,180,360]

        #Check these files exist, then update the directory
        target_fname = f'./data/{run_names[batch_id]}/chb_{snap_id:09d}.nc'

        if os.path.exists(target_fname):
            fcast.update_directory("base", run_names[batch_id], snap_id, [obs_time, rss, model, source, resolutions])
            fcast.update_directory("chbetc", run_names[batch_id], snap_id, [21.5])  #This should update the log of what has been calculated already. Will be tricksy, I think, to make sure the data stays uncorrupted.
        else:
            print('File not found')

