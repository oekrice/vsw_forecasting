#This is a testbed for making the forecast scripts actually nice, and doing it all properly and things. HA, that went well!

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
#This script should just run the base model and HuxT, at a low resolution.
#Will automatically create a run ID with parameters encoded into the outputs, one hopes.

start = datetime(2010, 1, 1) #This CAN'T change for a given run name. BE CAREFUL
obs_times = [start + timedelta(days=i) for i in range(5478)]
n_cores = 8

plot_specific = -1   #Just evalulate a specific snap. Set to -1 for the latest one
find_min_sigma = True

use_neural_net = False
#Specify input parameters as a dictionary, which can be embiggened or ensmallened as necessary.
#Will check against whether sufficient data exists which matches what has been asked for, and will recalculate if necessary.
#Let's specify literally everything here, all the parameters which can happen.
#Will need a lookup table or equivalent to find data which matches things as they should.
#Can specify file name to look up WSA parameters? Yeah, probably.
run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

base_selection = 3
batch_bases = ["wsa", "wsa_nocmes", "rms_nocmes", "corr_nocmes"]
titles1 = ["Default WSA", "Velocities Optimised for Distributions", "Velocities Optimised for RMS", "Velocities Optimised for Correlation"]

batch_base = batch_bases[base_selection]
title1 = titles1[base_selection]

fig1, axs1 = plt.subplots(2,4, figsize=(12,6))

for plot_num, batch_id in enumerate(np.arange(0,8)):


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

    batch_base = "combined_dynamic"

    #Have a flag to plot the default parameters if just 'wsa' is selected

    if batch_base == "wsa":
        print('Using default wsa parameters as a reference for later')
        batch_name = f"wsa_nocmes_{batch_id}"
        velocity_type = "wsa"

    else:
        batch_name = f"{batch_base}_{batch_id}"
        velocity_type = "wsa_scaled"

    # batch_name = f"net_test_{batch_id}"
    # velocity_type = "neural_net"
    # batch_base = "net_test"

    nicetitle = f"{model}, rss = {rss}, source = {source}"
    test_parameters = {"observation_time": obs_times,
                    "base_name": run_names[batch_id],
                    "run_name": batch_name,
                    "model_type": model,
                    "calculate_base_model": False,
                    "overwrite_base_model": False,
                    "calculate_huxt": False,
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
                    "optimisation_type": "distribution",
                    "do_plots": False}

    if not os.path.exists(f"./data/{test_parameters['run_name']}"):
        os.mkdir(f"./data/{test_parameters['run_name']}")

    if not os.path.exists(f"./plots/{test_parameters['run_name']}"):
        os.mkdir(f"./plots/{test_parameters['run_name']}")

    if len(sys.argv) > 2:
        if sys.argv[2] == 'ham8':
            print('Copying log file from Hamilton')
            os.system(f"scp -r vgjn10@hamilton8.dur.ac.uk:/nobackup/vgjn10/projects/vsw_forecasting/data/{test_parameters["run_name"]}/log.csv ./data/{test_parameters["run_name"]}/")
        else:
            raise Exception("Second argument not recognised. Use 'ham8' to copy Hamilton data")

    def load_directory():
        directory_fname = f'./data/{test_parameters["run_name"]}/log.csv'
        if os.path.exists(directory_fname):
            #This directory already exists. Hopefully with proper header information etc
            directory_data = []
            with open(directory_fname, "r", encoding="utf-8") as f:
                data = csv.reader(f)
                for row in data:
                    directory_data.append(row)
        else:
            print(f'Directory data not found with fname {directory_fname}')

        #Run through loaded directory and deduce information
        scores = []
        sigmas = []
        thetas = []
        for row in directory_data[1:]:
            scores.append(float(row[1]))
            sigmas.append(float(row[2]))
            thetas.append(row[3:])

        scores = np.array(scores)
        sigmas = np.array(sigmas)
        thetas = np.array(thetas, dtype='float')

        return scores, sigmas, thetas

    nsamples = 1
    valid_snaps = np.arange(5478)
    random.shuffle(valid_snaps)
    snap_subset = np.array([0] + list(valid_snaps[:nsamples-1]))

    try:
        scores, sigmas, thetas = load_directory()
    except:
        print('Directory not found...')
        continue


    if find_min_sigma:   #Use the parameters at the point at which the solution appears to have converged the best
        sigma_index = fcast.stats_functions.find_ideal_sigma_index(sigmas, overall_min=True, rawscores = scores)
        theta = thetas[-1]
    else:
        theta = thetas[-1]


    nthetas = np.shape(thetas)[0]

    if plot_specific < 0:
        i = len(thetas) - 1
    else:
        i = plot_specific

    iteration = 0
    print('Current theta', theta)
    if batch_base == "wsa":
        _, cmaps = fcast.model_functions.compute_vr(0, run_name, method=test_parameters["velocity_type"], params = None, doplot=test_parameters["do_plots"], iteration=iteration, huxt_name=test_parameters["run_name"], output_cmaps=True)
    else:
        _, cmaps = fcast.model_functions.compute_vr(0, run_name, method=test_parameters["velocity_type"], params = theta, doplot=test_parameters["do_plots"], iteration=iteration, huxt_name=test_parameters["run_name"], output_cmaps=True)
    vmesh, vr = cmaps


    ax = axs1[plot_num//4, plot_num%4]
    ax.pcolormesh(vr[30:-30,:], vmax=800, vmin=200)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(nicetitle)

plt.suptitle(title1)
plt.tight_layout()
plt.savefig(f'./plots/velocities_{batch_base}.png')
plt.close()




