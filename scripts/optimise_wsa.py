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

import matplotlib
matplotlib.use('Agg')
#This script should just run the base model and HuxT, at a low resolution.
#Will automatically create a run ID with parameters encoded into the outputs, one hopes.

start = datetime(2010, 1, 1) #This CAN'T change for a given run name. BE CAREFUL
obs_times = [start + timedelta(days=i) for i in range(5478)]
test_single =  False
n_cores = 8
nsamples = 50

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
                "verbose": False,
                "optimisation_type": "distribution",
                "do_plots": False}

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

    skillscores = fcast.run_model(test_parameters, theta=theta, snap_subset=snap_subset)

    minimiser = np.mean(skillscores)
    return minimiser

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

    es = cma.CMAEvolutionStrategy(np.zeros(8), 0.1, {'verb_disp': 1, 'popsize': popsize})
    if os.path.exists(f'data/{test_parameters["run_name"]}/log.csv'):
        os.remove(f'data/{test_parameters["run_name"]}/log.csv')

    best_losses = []
    sigmas = []
    with mp.Pool(processes=n_cores) as pool:
        while not es.stop():


            nsamples = 25
            valid_snaps = np.arange(5478)
            random.shuffle(valid_snaps)
            snap_subset = valid_snaps[:nsamples]

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

            fcast.update_theta_record(test_parameters, best_loss, es.sigma, best_theta)  #This keeps a record of which thetas are good, and the sigma at that time.

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
    valid_snaps = np.arange(5478)
    random.shuffle(valid_snaps)
    snap_subset = valid_snaps[:nsamples]

    skillscore = evaluate_theta(np.zeros(8), snap_subset=snap_subset)
    print('Current skillscore', skillscore)

