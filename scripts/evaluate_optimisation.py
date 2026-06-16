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

run_name = "optimise_run_2"
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

if len(sys.argv) > 1:
    batch_id = int(sys.argv[1])
else:
    raise Exception('Specify batch number.')

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

test_parameters = {"observation_time": obs_times,
                "base_name": run_names[batch_id],
                "run_name": f"optimise_run_{batch_id}",
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

def evaluate_theta(theta, snap_subset, iteration):
    """
    Will carry on even if there are errors.
    """

    skillscores = fcast.run_model(test_parameters, theta=theta, snap_subset=snap_subset, iteration=iteration)

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

nsamples = 25
valid_snaps = np.arange(5478)
random.shuffle(valid_snaps)
snap_subset = np.array([0] + list(valid_snaps[:nsamples-1]))

if not plot_continuous:

    scores, sigmas, thetas = load_directory()

    fig, axs = plt.subplots(2, figsize = (10,7))
    axs[0].plot(scores)
    axs[1].plot(sigmas)
    plt.savefig('./plots/%s/converge.png' % test_parameters["run_name"])
    plt.close()

    nthetas = np.shape(thetas)[0]

    if plot_specific < 0:
        i = len(thetas) - 1
    else:
        i = plot_specific

    print('Current theta', thetas[i])
    skillscore = evaluate_theta(thetas[i], snap_subset=snap_subset, iteration=i)
    print('Current skillscore', skillscore)

else:
    print('Waiting for continuous updates of performance')
    counter = 0
    while True:
        scores, sigmas, thetas = load_directory()
        nthetas = np.shape(thetas)[0]
        if nthetas > counter:
            i = nthetas - 1
            fig, axs = plt.subplots(2, figsize = (10,7))
            axs[0].plot(scores)
            axs[1].plot(sigmas)
            plt.savefig('./plots/%s/converge.png' % test_parameters["run_name"])
            plt.close()

            print('Current theta', thetas[i])
            skillscore = evaluate_theta(thetas[i], snap_subset=snap_subset, iteration=i)
            print('Current skillscore', skillscore)
            counter = nthetas
        else:
            time.sleep(5.0)




