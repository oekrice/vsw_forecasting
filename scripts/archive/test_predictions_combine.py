#This script is for taking an optimised distribution (or otherwise) and outputting (saving) a timeseries of the predicited velocities.

#Initially this script is going to test the three types of raw WSA parameters suggested in the 2010 paper. And then perhaps my own versions based on the relationships accounting for changes in CHB and expansion factors.

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import cma
import csv
import multiprocessing as mp
from datetime import datetime, timedelta

import wind_forecast as fcast  #This should now contain everything we need...
from dtaidistance import dtw
import random
import time

import matplotlib
#matplotlib.use('Agg')
#This script should just run the base model and HuxT, at a low resolution.
#Will automatically create a run ID with parameters encoded into the outputs, one hopes.

if len(sys.argv) > 1:
    batch_select = int(sys.argv[1])
else:
    raise Exception('Specify batch number.')

if len(sys.argv) > 2:
    base_selection = int(sys.argv[2])
else:
    raise Exception('Specify base selection.')

start = datetime(2010, 1, 1) #This CAN'T change for a given run name. BE CAREFUL
#obs_times = [start + timedelta(days=i) for i in range(5478)]
obs_times = [start + timedelta(days=i) for i in range(0, 5478)]

n_cores = 8

plot_specific = -1   #Just evalulate a specific snap. Set to -1 for the latest one
plot_continuous = True   #Will wait for outputs and keep up (if possible)
find_min_sigma = True
use_neural_net = True

#Specify input parameters as a dictionary, which can be embiggened or ensmallened as necessary.
#Will check against whether sufficient data exists which matches what has been asked for, and will recalculate if necessary.
#Let's specify literally everything here, all the parameters which can happen.
#Will need a lookup table or equivalent to find data which matches things as they should.
#Can specify file name to look up WSA parameters? Yeah, probably.
run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

batch_bases = ["wsa", "wsa_nocmes", "rms_nocmes", "corr_nocmes"]
titles1 = ["Default WSA", "Velocities Optimised for Distributions", "Velocities Optimised for RMS", "Velocities Optimised for Correlation"]

batch_base = batch_bases[base_selection]
title1 = titles1[base_selection]

snap_subset = np.arange(len(obs_times)) #This should just start from the start now

print('Running model with base velocities', title1, run_names[batch_select])
for plot_num, batch_id in enumerate(np.arange(batch_select,batch_select+1)):

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

    nicetitle = f"{model}, rss = {rss}, source = {source}"

    if batch_base == "wsa":
        print('Using default wsa parameters as a reference for later')
        batch_name = f"combine_raw _{batch_id}"
        velocity_type = "wsa"

    nicetitle = f"{model}, rss = {rss}, source = {source}"
    test_parameters = {"observation_time": obs_times,
                    "base_name": run_names[batch_id],
                    "run_name": batch_name,
                    "model_type": model,
                    "calculate_base_model": False,
                    "overwrite_base_model": False,
                    "calculate_huxt": True,
                    "r_ss": rss,
                    "WSA_type": "standard",
                    "WSA_parameters": None,
                    "data_source": source,
                    "resolutions": [120,180,360],
                    "r_hb": 21.5,
                    "match_flag": False,
                    "velocity_type": velocity_type,
                    "spinup_time": 5,
                    "forecast_length": 5,
                    "verbose": True,
                    "optimisation_type": "least_squares",
                    "do_plots": False}

    if not os.path.exists(f"./data/{test_parameters['run_name']}"):
        os.mkdir(f"./data/{test_parameters['run_name']}")

    if not os.path.exists(f"./plots/{test_parameters['run_name']}"):
        os.mkdir(f"./plots/{test_parameters['run_name']}")

    if find_min_sigma:   #Use the parameters at the point at which the solution appears to have converged the best
        sigma_index = fcast.stats_functions.find_ideal_sigma_index(sigmas)
        theta = thetas[-1]
    else:
        theta = thetas[-1]
    nthetas = np.shape(thetas)[0]

    iteration = 0
    print('Current theta', theta)

    fcast.model_functions.run_model(test_parameters, theta=theta, snap_subset=snap_subset, iteration=iteration, save_speeds=True)

