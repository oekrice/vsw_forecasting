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
matplotlib.use('Agg')
#This script should just run the base model and HuxT, at a low resolution.
#Will automatically create a run ID with parameters encoded into the outputs, one hopes.

start = datetime(2010, 1, 1) #This CAN'T change for a given run name. BE CAREFUL
obs_times = [start + timedelta(days=i) for i in range(5478)]
n_cores = 8

plot_specific = -1   #Just evalulate a specific snap. Set to -1 for the latest one
plot_continuous = True   #Will wait for outputs and keep up (if possible)

#Specify input parameters as a dictionary, which can be embiggened or ensmallened as necessary.
#Will check against whether sufficient data exists which matches what has been asked for, and will recalculate if necessary.
#Let's specify literally everything here, all the parameters which can happen.
#Will need a lookup table or equivalent to find data which matches things as they should.
#Can specify file name to look up WSA parameters? Yeah, probably.
run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

nsamples = 250
valid_snaps = np.arange(5478)
random.shuffle(valid_snaps)
snap_subset = np.array([0] + list(valid_snaps[:nsamples-1]))

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

    nicetitle = f"{model}, rss = {rss}, source = {source}"

    test_parameters = {"observation_time": obs_times,
                    "base_name": run_names[batch_id],
                    "run_name": f"combined_dynamic_{batch_id}",
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
                    "velocity_type": "wsa",
                    "spinup_time": 5,
                    "forecast_length": 5,
                    "verbose": True,
                    "optimisation_type": "distribution",
                    "do_plots": True}

    if not os.path.exists(f"./data/{test_parameters['run_name']}"):
        os.mkdir(f"./data/{test_parameters['run_name']}")

    if not os.path.exists(f"./plots/{test_parameters['run_name']}"):
        os.mkdir(f"./plots/{test_parameters['run_name']}")

    print('Copying log file from Hamilton')
    os.system(f"scp -r vgjn10@hamilton8.dur.ac.uk:/nobackup/vgjn10/projects/vsw_forecasting/data/{test_parameters["run_name"]}/log.csv ./data/{test_parameters["run_name"]}/")

    def evaluate_theta(theta, snap_subset, iteration):
        """
        Will carry on even if there are errors.
        """

        skillscores = fcast.run_model(test_parameters, theta=None, snap_subset=snap_subset, iteration=iteration)

        minimiser = np.mean(skillscores)
        return minimiser

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

    scores, sigmas, thetas = load_directory()

    sigma_index = fcast.stats_functions.find_ideal_sigma_index(sigmas)
    fig, axs = plt.subplots(2, figsize = (10,7))
    axs[0].plot(scores)
    axs[0].scatter([sigma_index],[scores[sigma_index]], c = 'red')
    axs[1].plot(sigmas)
    axs[1].set_ylim(0.0,0.5)
    axs[1].scatter([sigma_index],[sigmas[sigma_index]], c = 'red')
    plt.savefig(f'./plots/converges/converge_{test_parameters["run_name"]}.png' % (batch_id))
    plt.close()

