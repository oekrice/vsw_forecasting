#Script for plotting colourmaps of the chb, expansion factors and velocities for a given snapshot (will use halfway through the 'good' rotation?)_plot

#This is a testbed for making the forecast scripts actually nice, and doing it all properly and things. HA, that went well!

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import cma
import csv
import time
import multiprocessing as mp
from datetime import datetime, timedelta

import wind_forecast as wf  #This should now contain everything we need...
from dtaidistance import dtw
import random
from scipy.optimize import minimize
from scipy.stats import pearsonr
import cmocean

from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

import matplotlib as mpl

#Matplotlib preamble (should use on all plots really for consistency)
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 12,        # Default font size
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
})
import matplotlib
#matplotlib.use('Agg')
#This script should just run the base model and HuxT, at a low resolution.
#Will automatically create a run ID with parameters encoded into the outputs, one hopes.

start = datetime(2010, 1, 1) #This CAN'T change for a given run name. BE CAREFUL
#obs_times = [start + timedelta(days=i) for i in range(5478)]
obs_times = [start + timedelta(days=i) for i in range(5478)]

fig_width = 443.57848/72

if "SLURM_JOB_ID" in os.environ:
    n_cores = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
    print('Number of slurm-allocated cores:', n_cores)
else:
    print('Running locally (not on slurm)')
    n_cores = 8


set_type = 'good'

def run_model_combined(run_parameters, theta=None, snap_subset=None, iteration=0, output_distributions=False, save_speeds=False, use_old_chb_formula=False):
    global global_vr
    global global_r
    global global_rms
    global ax_row
    #Alas I think this needs to be here now, as it's getting far too complicated.
            #Cycle through the reuqested observation times. Can be just one, or several
    if "observation_time" in run_parameters:
        if not isinstance(run_parameters["observation_time"], (list, np.ndarray)):
            run_parameters["observation_time"] = [run_parameters["observation_time"]]
    else:
        raise Exception('Observation time not provided.')

    if run_parameters["base_name"] is not None:
        if "match_flag" in run_parameters:
            data_lookup = wf.data_functions.check_existing_data(run_parameters, match_flag = run_parameters["match_flag"]) #Will eventually output flags as to whether these data have been already succesfully calculated with the given inputs.
        else:
            data_lookup = wf.data_functions.check_existing_data(run_parameters, match_flag = True)
    else:
        data_lookup = [0] * len(run_parameters["observation_time"])

    if run_parameters["model_type"] == "outflow":
        is_pfss = False
    elif run_parameters["model_type"] == "pfss":
        is_pfss = True
    else:
        raise Exception("Model type not recognised. Need 'outflow' or 'pfss'")

    if run_parameters["base_name"] is not None:
        run_name = run_parameters["base_name"]
    else:
        run_name = "tmp"

    if snap_subset is None:
        snap_subset = np.arange(len(run_parameters["observation_time"]))

    for snap_id in snap_subset:
        obs_time = run_parameters["observation_time"][snap_id]
        #This will do the base calculations (the ones which take some time). I think it makes sense to put these in this separate loop, at least for now.
        if run_parameters["calculate_base_model"]:  #Do the PFSS/Outflow calculation
            if run_parameters["verbose"]:
                print(f'Running base model at time {obs_time}')
            #Check if base model already exists for these data, and the metadata all match (should put this check in an extra function. Each 'data' file should have a lookup table for it, I think)
            if data_lookup[snap_id] < 1 or run_parameters["overwrite_base_model"]: #The base model data doesn't exist, so recalculate
                calculate_outflow(snap_id,  obs_time, rss = run_parameters["r_ss"],
                                        overwrite=True, is_pfss=is_pfss, output_directory=run_name, save_snap=True,
                                        source= run_parameters["data_source"],resolutions=run_parameters["resolutions"])
            if data_lookup[snap_id] < 2 or run_parameters["overwrite_base_model"]: #The CHB and expansion factor data doesn't exist, so do that.
                calculate_chb_exp(snap_id,  run_name, r_hb = run_parameters["r_hb"], purge_data=True, use_old_formula=use_old_chb_formula)

    if run_parameters["verbose"]:
        print('Coronal hole distances and expansion factors computed for all requested snaps')

    min_omni_time = min(run_parameters["observation_time"]) - timedelta(days=30)
    max_omni_time = max(run_parameters["observation_time"]) + timedelta(days=30)
    omni_fname = './data/shared_data/omni.csv'
    omni_data = []

    with open(omni_fname, "r", encoding="utf-8") as f:
        data = csv.reader(f)
        for row in data:
            if datetime.strptime(row[0], "%Y-%m-%d %H:%M:%S") >= min_omni_time and datetime.strptime(row[0], "%Y-%m-%d %H:%M:%S") <= max_omni_time:
                omni_data.append([datetime.strptime(row[0], "%Y-%m-%d %H:%M:%S"), float(row[1])])

    omni_data = np.array(omni_data)

    #OK, do the data processing here for shifting and things
    omni_dates = omni_data[:,0]
    omni_vs = omni_data[:,1]


    if run_parameters["verbose"]:
        print('OMNI Data Obtained')

    skillscores = []

    alltimes = []
    allspeeds = []
    allspeeds_ref = []


    if run_parameters["calculate_huxt"]:
        for si, snap_id in enumerate(snap_subset):
            obs_time = run_parameters["observation_time"][snap_id]
            if run_parameters["verbose"]:
                print(f'Running HuXT forecast model at time {obs_time}')

            if run_parameters["velocity_type"] == "wsa" or run_parameters["velocity_type"] == "wsa_scaled" or run_parameters["velocity_type"] == "wsa_combined":
                #This is the polynomial expression
                if theta is not None:
                    if si == 0:
                        vr, chb, fs = wf.field_calculations.compute_vr(snap_id, run_name, run_parameters["velocity_type"], params = theta, doplot=run_parameters["do_plots"], iteration=iteration, huxt_name=run_parameters["run_name"], output_cmaps=True)
                    else:
                        vr, chb, fs =  wf.field_calculations.compute_vr(snap_id, run_name, run_parameters["velocity_type"], params = theta, doplot=False, iteration=iteration, huxt_name=run_parameters["run_name"], output_cmaps=True)
                else:
                    if run_parameters["verbose"]:
                        print('Using default WSA parameters')
                    if si == 0:
                        vr, chb, fs =  wf.field_calculations.compute_vr(snap_id, run_name, method="wsa", params = None, doplot=run_parameters["do_plots"], iteration=iteration, huxt_name=run_parameters["run_name"], output_cmaps=True)
                    else:
                        vr, chb, fs =  wf.field_calculations.compute_vr(snap_id, run_name, method="wsa", params = None, doplot=False, iteration=iteration, huxt_name=run_parameters["run_name"], output_cmaps=True)

            elif run_parameters["velocity_type"] == "neural_net":
                if run_parameters["verbose"]:
                    print('Using neural net parameters', theta)
                if si == 0:
                    vr, chb, fs =  wf.field_calculations.compute_vr(snap_id, run_name, method="neural_net", params = theta, doplot=run_parameters["do_plots"], iteration=iteration, huxt_name=run_parameters["run_name"], output_cmaps=True)
                else:
                    vr, chb, fs =  wf.field_calculations.compute_vr(snap_id, run_name, method="neural_net", params = theta, doplot=False, iteration=iteration, huxt_name=run_parameters["run_name"], output_cmaps=True)

            else:
                raise Exception("Velocity calculation type not recognised...")

            if run_parameters["verbose"]:
                print('Calculating VSW with HuxT...')

            #Using this velocity map (and nothing else?), run HuxT
            #ax.plot(times, model_speeds, linewidth=0.2, c = 'green')

            #Use the map from the middle of the rotatio
            global_vr = vr
            #ax_row[0].pcolormesh(vr)

        return vr, chb, fs

def evaluate_theta(theta, snap_subset):
    """
    Will carry on even if there are errors.
    """


    return None

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

nsamples = len(obs_times)
extend_current_run = False
use_neural_net = True

if not use_neural_net:
    theta_size = 9
else:
    theta_size = 41

#Specify input parameters as a dictionary, which can be embiggened or ensmallened as necessary.
#Will check against whether sufficient data exists which matches what has been asked for, and will recalculate if necessary.
#Let's specify literally everything here, all the parameters which can happen.
#Will need a lookup table or equivalent to find data which matches things as they should.
#Can specify file name to look up WSA parameters? Yeah, probably.

iteration_number = 0
mappables = [None, None, None]

fig, axs = plt.subplots(2,3, figsize=(fig_width,fig_width*0.5), constrained_layout=True)  #This one is now for colourmaps

#Run through vanilla WSA, Otpimised WSA and Neural Net. To make the point.

for bi, batch_id in enumerate([0,7]):
    run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

    ax_row = axs[bi, :]

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

    colours = ['red', 'green', 'blue']
    titles =  ['Default WSA Parameters', 'Optimised WSA Parameters', 'Neural Net']
    global_vr = np.zeros((180,360))  #Put this here and it might work out in the end?
    global_r = 0.
    global_rms = 0.

    model_type = 0

    plot_colour = colours[model_type]
    batch_name = f"overfit_{set_type}_{batch_id}"
    velocity_type = "wsa"   #To be used for the physics-informed parameter-changing. Just a select few of them. See if there are any mad patterns.

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
                    "optimisation_type": "correlation",
                    "do_plots": False,
                    "filter_cmes": True}

    if not os.path.exists(f"./data/{test_parameters['run_name']}"):
        os.mkdir(f"./data/{test_parameters['run_name']}")

    #Save a log to let the thing know it's started, for logging purposes
    np.savetxt(f"./data/{test_parameters['run_name']}/start.dat", [n_cores])
    print('Running job with name', test_parameters['run_name'], 'using data', test_parameters['base_name'])

    selected_snap = 123

    valid_snaps = [selected_snap]
    snap_subset = valid_snaps

    vr, chbs, fs = run_model_combined(test_parameters, snap_subset=snap_subset)

    print(np.shape(vr), np.shape(chbs), np.shape(fs))
    xs = np.linspace(-180,180,361)
    ys = np.linspace(-90,90,181)
    toplots = [chbs, fs, vr]
    vmins = [0,0,250]
    vmaxs = [10,300,850]
    for i in range(3):
        ax_row[i].pcolormesh(xs, ys, toplots[i], vmin=vmins[i], vmax=vmaxs[i], rasterized=True, cmap=cmocean.cm.solar)

    for i in range(1,3):
        ax_row[i].set_yticks([])

    if bi == 0:
        for i in range(3):
            ax_row[i].set_xticks([])

    if bi == 0:
        ax_row[0].set_title('CH Boundary Distance')
        ax_row[1].set_title('Expansion Factor')
        ax_row[2].set_title('WSA Velocity')



for i in range(3):
    norm = mpl.colors.Normalize(vmin=vmins[i], vmax=vmaxs[i])

    sm = mpl.cm.ScalarMappable(
        norm=norm,
        cmap=cmocean.cm.solar
    )
    sm.set_array([])

    cbar = fig.colorbar(
        sm,
        ax=axs[:, i],
        orientation='horizontal',
        location='bottom',
        aspect=20,
        pad = 0.01
    )

# ax.set_xlabel('Time')
# ax.set_ylabel('Wind Speed')
# plt.setp(ax.get_xticklabels(), rotation=30, ha = "right")

# ax.set_xticks([])
# ax.set_yticks([])

#plt.legend(fontsize=12)
plt.savefig(f'./paper/plots/1_plot_cmaps.pdf', dpi=300, bbox_inches="tight")
plt.show()
plt.close()
