#This script is for comparing the overall distributions of coronal hole distances and expansion factors

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
from scipy.stats import pearsonr

import matplotlib

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 12,        # Default font size
    "axes.labelsize": 12,
    "axes.titlesize": 8,
    "xtick.labelsize": 12,
    "ytick.labelsize": 8,
})

fig_width = 443.57848/72

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

def make_nicetitle(id):
    if (id//2)%2 == 0:
        is_pfss = True
        model = "PFSS"
    else:
        is_pfss = False
        model = "Outflow"
    if (id%2) == 0:
        rss = 2.5
    else:
        rss = 5.0
    if (id//4) == 0:
        source = "GONG"
    else:
        source = "HMI"

    rss_string = "r_{ss}"
    nicetitle = f"{model}, ${rss_string} = {rss}$, {source}"
    return nicetitle

if False:   #Compile and save the data

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

        print(run_name)
        allchds = []
        allfs = []
        nsnaps = 100 #len(obs_times)

        for snap_id in np.linspace(0,len(obs_times)-1, nsnaps):
            snap_id = int(snap_id)
            #print(snap_id/len(obs_times))
            s0, ph0, br0, fs, chd = fcast.data_functions.load_chb_distances(run_name, snap_id)
            allchds = allchds + list(chd[60:120,:].flatten()[::100])
            allfs = allfs + list(fs[60:120,:].flatten()[::100])
            del s0, ph0, br0, fs, chd

        print('Data length', len(allchds))
        np.save(f'./data/factor_comparison/chds_{plot_num}.npy', allchds)
        np.save(f'./data/factor_comparison/fs_{plot_num}.npy', allfs)

else:   #Analyse the data
    chb_ref = np.load(f'./data/factor_comparison/chds_0.npy')
    fs_ref = np.load(f'./data/factor_comparison/fs_0.npy')

    chb_biases = []
    print(len(chb_ref))
    fig1, axs1 = plt.subplots(2,4, figsize=(fig_width, fig_width/2))
    for plot_num, batch_id in enumerate(np.arange(0,8)):


        chb = np.load(f'./paper/data/factor_comparison/chds_{plot_num}.npy')

        ax = axs1[plot_num//4, plot_num%4]

        ax.set_xticks([])
        ax.set_yticks([])
        r, _ = pearsonr(chb, chb_ref)
        bias = np.mean(chb)/np.mean(chb_ref)
        ax.set_title(f"{make_nicetitle(plot_num)}, \n r = {r:.3f}, bias = {bias:.3f}")
        ax.scatter(chb_ref, chb, alpha = 0.1, c = 'blue', s = 0.1, rasterized=True)
        chb_biases.append(bias)

    plt.tight_layout()
    plt.savefig('./paper/plots/4_compare_chbs.pdf', dpi=600)
    plt.show()
    plt.close()

    # np.savetxt('./data/shared_data/chb_biases.txt', chb_biases, delimiter = ',')
    #
    # fs_biases = []
    #
    # fig1, axs1 = plt.subplots(2,4, figsize=(12,6))
    # for plot_num, batch_id in enumerate(np.arange(0,8)):
    #
    #
    #     fs = np.load(f'./paper/data/factor_comparison/fs_{plot_num}.npy')
    #
    #     ax = axs1[plot_num//4, plot_num%4]
    #     ax.set_xticks([])
    #     ax.set_yticks([])
    #     r, _ = pearsonr(fs, fs_ref)
    #     bias = np.mean(fs)/np.mean(fs_ref)
    #     ax.set_title(f"{make_nicetitle(plot_num)}, \n r = {r:.3f}, bias = {bias:.3f}")
    #     ax.scatter(fs_ref, fs, alpha = 0.1, c = 'blue', s = 0.1)
    #     fs_biases.append(bias)
    #
    # plt.suptitle('Expansion Factors')
    # plt.tight_layout()
    # plt.savefig('./plots/fs_corr.png')
    # plt.show()
    # plt.close()
    #
    # np.savetxt('./data/shared_data/fs_biases.txt', fs_biases, delimiter = ',')
