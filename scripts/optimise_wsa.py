#This is a testbed for making the forecast scripts actually nice, and doing it all properly and things. HA, that went well!

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import cma
import csv
import multiprocessing as mp
from datetime import datetime, timedelta

import wind_forecast as wf  #This should now contain everything we need...
from dtaidistance import dtw
import random

import matplotlib
matplotlib.use('Agg')
#This script should just run the base model and HuxT, at a low resolution.
#Will automatically create a run ID with parameters encoded into the outputs, one hopes.

start = datetime(2010, 1, 1) #This CAN'T change for a given run name. BE CAREFUL
obs_times = [start + timedelta(days=i) for i in range(5478)]
test_single =  False

if "SLURM_JOB_ID" in os.environ:
    n_cores = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
    print('Number of slurm-allocated cores:', n_cores)
else:
    print('Running locally (not on slurm)')
    n_cores = 8

nsamples = 50
extend_current_run = True
use_neural_net = True

if not use_neural_net:
    theta_size = 8
else:
    theta_size = 41

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

if use_neural_net:
    batch_name = f"net_test_{batch_id}"
    velocity_type = "neural_net"
else:
    batch_name = f"corr_nocmes_{batch_id}"
    velocity_type = "wsa_scaled"

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
                "verbose": False,
                "optimisation_type": "correlation",
                "do_plots": False,
                "filter_cmes": True}

if not os.path.exists(f"./data/{test_parameters['run_name']}"):
    os.mkdir(f"./data/{test_parameters['run_name']}")

#Save a log to let the thing know it's started, for logging purposes
np.savetxt(f"./data/{test_parameters['run_name']}/start.dat", [n_cores])
print('Running job with name', test_parameters['run_name'], 'using data', test_parameters['base_name'])

def evaluate_with_timeout(pool, theta, timeout=60):
    result = pool.apply_async(evaluate_theta_safe, (theta,))
    try:
        return result.get(timeout=timeout)
    except Exception:
        return 1e12

def evaluate_theta(theta, snap_subset):
    """
    Will carry on even if there are errors.
    """

    skillscores = wf.model_functions.run_model(test_parameters, theta=theta, snap_subset=snap_subset)

    minimiser = np.mean(skillscores)

    np.savetxt(f"./data/{test_parameters['run_name']}/start_theta.dat", [theta])

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

def run_cma_mp(n_cores=None):

    if n_cores is None:
        n_cores = int(os.environ.get("SLURM_CPUS_PER_TASK", mp.cpu_count()))

    pool = mp.Pool(processes=n_cores)
    best_loss = float("inf")
    best_theta = None

    popsize = n_cores
    while popsize < 16:
         popsize += n_cores

    print('Ncores:', n_cores, 'Population size', popsize)


    if not extend_current_run:
        if os.path.exists(f'data/{test_parameters["run_name"]}/log.csv'):
            os.remove(f'data/{test_parameters["run_name"]}/log.csv')

        es = cma.CMAEvolutionStrategy(np.zeros(theta_size), 0.1, {'verb_disp': 1, 'popsize': popsize})
        best_losses = []
        sigmas = []

    else:
        if os.path.exists(f'data/{test_parameters["run_name"]}/log.csv'):
            scores, sigmas, thetas = load_directory()
            sigmas = list(sigmas)
            best_losses = list(scores)
            es = cma.CMAEvolutionStrategy(thetas[-1], sigmas[-1], {'verb_disp': 1, 'popsize': popsize})
            print('Using existing run, initial conditions', thetas[-1], sigmas[-1])
        else:
            best_losses = []
            sigmas = []
            es = cma.CMAEvolutionStrategy(np.zeros(theta_size), 0.1, {'verb_disp': 1, 'popsize': popsize})

    valid_snaps = np.arange(len(obs_times))#[1-cme_mask]

    if test_parameters["filter_cmes"]:
        cme_mask = wf.data_functions.get_cme_times(obs_times)
        valid_times = np.where(cme_mask == 0)[0]
        valid_snaps = valid_snaps[valid_times]
    random.shuffle(valid_snaps)

    snap_subset = valid_snaps[:nsamples].copy()

    with mp.Pool(processes=n_cores) as pool:
        while not es.stop():

            if (len(sigmas)%25) == 0:  #I've not really tested whether this makes any meaningful difference...
                random.shuffle(valid_snaps)
                snap_subset = valid_snaps[:nsamples].copy()

            solutions = es.ask()

            results = [
                pool.apply_async(evaluate_theta, (theta,snap_subset))
                for theta in solutions
            ]

            losses = []
            for r in results:
                try:
                    losses.append(r.get(timeout=60.0))
                except Exception:
                    losses.append(1e12)

            es.tell(solutions, losses)
            best_loss = 1e6
            best_theta = None
            for theta, loss in zip(solutions, losses):
                print("Current score", loss)

                if loss < best_loss:
                    best_loss = loss
                    best_theta = theta.copy()

            #Run the model with the best ones to do a plot of progress?
            #fcast.run_model(test_parameters, best_theta, iteration = len(best_losses), snap_subset=[0])

            best_losses.append(best_loss)
            sigmas.append(es.sigma)

            wf.data_functions.update_theta_record(test_parameters, best_loss, es.sigma, best_theta)  #This keeps a record of which thetas are good, and the sigma at that time.

            if not os.path.exists('plots'):
                os.mkdir('plots')
            if not os.path.exists(f'plots/{test_parameters["run_name"]}'):
                os.mkdir(f'plots/{test_parameters["run_name"]}')

            if False:
                fig, axs = plt.subplots(2, figsize = (10,7))
                axs[0].plot(best_losses)
                axs[1].plot(sigmas)
                plt.savefig('./plots/%s/converge.png' % test_parameters["run_name"])
                plt.close()

    return

if not test_single:
    if __name__ == "__main__":
        run_cma_mp(n_cores=n_cores)

else:
    valid_snaps = np.arange(len(obs_times))#[1-cme_mask]
    filter_for_cmes = True
    if filter_for_cmes:
        cme_mask = wf.data_functions.get_cme_times(obs_times)
        valid_times = np.where(cme_mask == 0)[0]
        valid_snaps = valid_snaps[valid_times]

    random.shuffle(valid_snaps)
    snap_subset = valid_snaps[:nsamples]

    skillscore = evaluate_theta(np.zeros(theta_size), snap_subset=snap_subset)
    print('Current skillscore', skillscore)

