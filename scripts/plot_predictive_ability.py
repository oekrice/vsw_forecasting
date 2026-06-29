#This script requires test_predictions to have been run beforehand. Just loads in text files of speeds, which are hopefully aligned, and calculates/plots various things
#Seems to be best to do things separately like this to allow for a lot of flexibility.
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import wind_forecast as fcast
from scipy.ndimage import gaussian_filter1d
from datetime import datetime, timedelta

filter_for_cmes = True

cmap = plt.get_cmap("tab10")

if len(sys.argv) > 1:
    plot_type = int(sys.argv[1])
else:
    plot_type = -1

def scale_function(m, c, series):
    """
    Just does mx+c on the timeseries. Can optimise the skill scores based on that, hopefully
    """
    return series*m + c*100

#Load in the scaling (generated using scale_for_persistence.py)

unoptimised_scales = []
optimised_scales = []

do_ss_scaling = True
if not do_ss_scaling:  #Do linear scaling based on rms, not skillscore
    for i in range(8):
        opt_title = 'raw'
        batch_name = f'wsa_nocmes_{i}'
        if os.path.exists(f'./data/scaling_data/{batch_name}_{opt_title}_rms.txt'):
            scale_data = np.loadtxt(f'./data/scaling_data/{batch_name}_{opt_title}_rms.txt')
            unoptimised_scales.append(scale_data[1:])
        else:
            unoptimised_scales.append(None)

        opt_title = 'dist'

        if os.path.exists(f'./data/scaling_data/{batch_name}_{opt_title}_rms.txt'):
            scale_data = np.loadtxt(f'./data/scaling_data/{batch_name}_{opt_title}_rms.txt')
            optimised_scales.append(scale_data[1:])
        else:
            optimised_scales.append(None)

else:
    for i in range(8):
        opt_title = 'raw'
        batch_name = f'wsa_nocmes_{i}'
        if os.path.exists(f'./data/scaling_data/{batch_name}_{opt_title}.txt'):
            scale_data = np.loadtxt(f'./data/scaling_data/{batch_name}_{opt_title}.txt')
            unoptimised_scales.append(scale_data[1:])
        else:
            unoptimised_scales.append(None)

        opt_title = 'dist'

        if os.path.exists(f'./data/scaling_data/{batch_name}_{opt_title}.txt'):
            scale_data = np.loadtxt(f'./data/scaling_data/{batch_name}_{opt_title}.txt')
            optimised_scales.append(scale_data[1:])
        else:
            optimised_scales.append(None)

# unoptimised_scales = [None, None, None, None, None, None, None, [1.08, -0.156]]
# optimised_scales = [None, None, None, None, None, None, None, [ 0.8411262,  -0.05513745]]


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

    nicetitle = f"{model}, r_ss = {rss}, {source}"
    return nicetitle


if plot_type == -1 or plot_type == 0: #Do timeseries and print out RMS values. Alas these appear to be consistently worse once optimised. Bugger. Yes.
    batch_names = []

    do_scaled = True
    do_optimised = True

    for i in range(8):
        batch_names.append(f'wsa_nocmes_{i}')

    fig = plt.figure(figsize=(12,6))
    for i, batch_name in enumerate(batch_names):

        #Hopefully all things should be arranged nicely time-wise, but do need to check as much
        wsa_fname = f'./data/raw_speeds/{batch_name}_wsa_speeds.txt'
        optimised_fname = f'./data/raw_speeds/{batch_name}_wsa_scaled_speeds.txt'
        ref_fname_0 = f'./data/raw_speeds/{batch_name}_wsa_speeds_ref.txt'
        ref_fname_1 = f'./data/raw_speeds/{batch_name}_wsa_scaled_speeds_ref.txt'


        if not(os.path.exists(wsa_fname) and os.path.exists(optimised_fname) and os.path.exists(ref_fname_0) and os.path.exists(ref_fname_1)):
            print('Files not found...', wsa_fname, optimised_fname, ref_fname_0, ref_fname_1)
            continue

        wsa = np.loadtxt(wsa_fname, delimiter = ',')
        optimised_wsa = np.loadtxt(optimised_fname, delimiter = ',')
        omni_ref_0 = np.loadtxt(ref_fname_0, delimiter = ',')
        omni_ref_1 = np.loadtxt(ref_fname_1, delimiter = ',')

        timeseries = np.loadtxt(f'./data/raw_speeds/{batch_name}_wsa_scaled_times.txt', dtype='datetime64[s]', delimiter = ',')

        filter_for_cmes = False
        if filter_for_cmes:
            cme_mask = fcast.data_functions.get_cme_times(timeseries)
            invalid_times = np.where(cme_mask == 1)[0]
            wsa[invalid_times] = np.nan
            optimised_wsa[invalid_times] = np.nan
            omni_ref_0[invalid_times] = np.nan
            omni_ref_1[invalid_times] = np.nan

        if do_scaled:  #Apply scaling as determined by scale_for_persistence.py
            if unoptimised_scales[i] is not None and optimised_scales[i] is not None:
                wsa = scale_function(unoptimised_scales[i][0], unoptimised_scales[i][1], wsa)
                optimised_wsa = scale_function(optimised_scales[i][0], optimised_scales[i][1], optimised_wsa)

        data_lengths = [len(wsa), len(optimised_wsa), len(omni_ref_0), len(omni_ref_1)]

        if not min(data_lengths) == max(data_lengths):
            raise Exception("Timeseries data lengths don't match. Not sure what to do...")

        if do_optimised:
            abs_difference = np.abs(optimised_wsa-omni_ref_0)
        else:
            abs_difference = np.abs(wsa-omni_ref_0)
        #Split this up into some time chunks and average. Otherwise is quite a mess. Current timeseries is every day.
        #Need to do this and THEN get the RMS values once CMEs are taken into account
        difference_cadence = 30*24 #Number of hours to split

        time_slices = []
        diff_slices = []
        i_min = 0; i_max = difference_cadence
        while i_max < len(abs_difference):
            time_slices.append(timeseries[i_min] + 0.5*(timeseries[i_max]-timeseries[i_min]))
            diff_slices.append(np.nanmean(abs_difference[i_min:i_max]))
            i_min = i_min + difference_cadence
            i_max = i_max + difference_cadence

        #Plot the differences NOT removing the CMEs

        # Then get the actual RMS figures using the CME filter. This is all awfully complicated.
        filter_for_cmes = True
        if filter_for_cmes:
            cme_mask = fcast.data_functions.get_cme_times(timeseries)
            invalid_times = np.where(cme_mask == 1)[0]
            wsa[invalid_times] = np.nan
            optimised_wsa[invalid_times] = np.nan
            omni_ref_0[invalid_times] = np.nan
            omni_ref_1[invalid_times] = np.nan

        if do_optimised:
            rms = np.sqrt(np.nanmean((optimised_wsa - omni_ref_0)**2))
        else:
            rms = np.sqrt(np.nanmean((wsa - omni_ref_0)**2))

        # print('Standard STD for ', batch_name, np.sqrt(np.nanmean((wsa-omni_ref_0)**2)))
        # print('Optimised STD for', batch_name, np.sqrt(np.nanmean((optimised_wsa-omni_ref_1)**2)))
        plt.plot(time_slices, diff_slices, label = f'{make_nicetitle(i)}, rms = {rms:.0f}km/s')

    plt.legend(fontsize=10)
    plt.title('Mean absolute wind speed error, optimised for distributions and scaled for Skill Score')
    plt.ylim(0,500)
    plt.tight_layout()
    plt.savefig('./plots/errors_dist_ss.png')
    plt.show()

if plot_type == -1 or plot_type == 1: #Do histogram comparison
    batch_names = []
    for i in range(8):
        batch_names.append(f'wsa_nocmes_{i}')

    fig = plt.figure(figsize=(12,6))
    for i, batch_name in enumerate(batch_names):
        wsa_fname = f'./data/raw_speeds/{batch_name}_wsa_speeds.txt'
        optimised_fname = f'./data/raw_speeds/{batch_name}_wsa_scaled_speeds.txt'
        ref_fname_0 = f'./data/raw_speeds/{batch_name}_wsa_speeds_ref.txt'
        ref_fname_1 = f'./data/raw_speeds/{batch_name}_wsa_scaled_speeds_ref.txt'

        if not(os.path.exists(wsa_fname) and os.path.exists(optimised_fname) and os.path.exists(ref_fname_0) and os.path.exists(ref_fname_1)):
            continue

        wsa = np.loadtxt(wsa_fname, delimiter = ',')
        optimised_wsa = np.loadtxt(optimised_fname, delimiter = ',')
        omni_ref_0 = np.loadtxt(ref_fname_0, delimiter = ',')
        omni_ref_1 = np.loadtxt(ref_fname_1, delimiter = ',')

        timeseries = np.loadtxt(f'./data/raw_speeds/{batch_name}_wsa_scaled_times.txt', dtype='datetime64[s]', delimiter = ',')

        if filter_for_cmes:
            cme_mask = fcast.data_functions.get_cme_times(timeseries)
            invalid_times = np.where(cme_mask == 1)[0]
            wsa[invalid_times] = np.nan
            optimised_wsa[invalid_times] = np.nan
            omni_ref_0[invalid_times] = np.nan
            omni_ref_1[invalid_times] = np.nan



        nbins = 101
        hist_wsa, _ = np.histogram(wsa, bins=nbins, range=(0.0,1000.0))
        hist_optimised_wsa, _ = np.histogram(optimised_wsa, bins=nbins, range=(0.0,1000.0))
        hist_ref_0, _ = np.histogram(omni_ref_0, bins=nbins, range=(0.0,1000.0))
        hist_ref_1, _ = np.histogram(omni_ref_1, bins=nbins, range=(0.0,1000.0))

        hist_wsa = hist_wsa/np.sum(hist_wsa)
        hist_optimised_wsa = hist_optimised_wsa/np.sum(hist_optimised_wsa)
        hist_ref = hist_ref_0/np.sum(hist_ref_0)
        distance0 = np.sum((hist_wsa - hist_ref)**2)
        distance1 = np.sum((hist_optimised_wsa - hist_ref)**2)

        print('Standard distance:', distance0)
        print('Optimised distance:', distance1)
        plt.plot(hist_wsa)
        plt.plot(hist_optimised_wsa)
        plt.plot(hist_ref, c = 'black')

        plt.show()

if plot_type == -1 or plot_type == 2: #Do 'persistence metric' or equivalent, for a variety of thresholds.

    run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

    doruns = np.arange(8)
    batch_names = []
    for i in range(8):
        batch_names.append(f'wsa_nocmes_{i}')

    for i in doruns:
        fig = plt.figure(figsize=(12,6))

        batch_name = batch_names[i]
        wsa_fname = f'./data/raw_speeds/{batch_name}_wsa_speeds.txt'
        optimised_fname = f'./data/raw_speeds/{batch_name}_wsa_scaled_speeds.txt'
        ref_fname_0 = f'./data/raw_speeds/{batch_name}_wsa_speeds_ref.txt'
        ref_fname_1 = f'./data/raw_speeds/{batch_name}_wsa_scaled_speeds_ref.txt'

        timeseries = np.loadtxt(f'./data/raw_speeds/{batch_name}_wsa_scaled_times.txt', dtype='datetime64[s]', delimiter = ',')

        if not(os.path.exists(wsa_fname) and os.path.exists(optimised_fname) and os.path.exists(ref_fname_0) and os.path.exists(ref_fname_1)):
            continue

        wsa = np.loadtxt(wsa_fname, delimiter = ',')
        optimised_wsa = np.loadtxt(optimised_fname, delimiter = ',')
        omni_ref_0 = np.loadtxt(ref_fname_0, delimiter = ',')
        omni_ref_1 = np.loadtxt(ref_fname_1, delimiter = ',')

        if filter_for_cmes:
            cme_mask = fcast.data_functions.get_cme_times(timeseries, also_filter_persistence=False)
            invalid_times = np.where(cme_mask == 1)[0]
            wsa[invalid_times] = np.nan
            optimised_wsa[invalid_times] = np.nan
            omni_ref_0[invalid_times] = np.nan
            omni_ref_1[invalid_times] = np.nan

        if True:  #Apply scaling as determined by scale_for_persistence.py
            if unoptimised_scales[i] is not None and optimised_scales[i] is not None:
                wsa = scale_function(unoptimised_scales[i][0], unoptimised_scales[i][1], wsa)
                optimised_wsa = scale_function(optimised_scales[i][0], optimised_scales[i][1], optimised_wsa)

        thresholds = np.arange(450,600,5)
        scores = fcast.stats_functions.do_met_stats(None, wsa, omni_ref_0)

        plt.plot(thresholds, scores, label = f'{make_nicetitle(i)}, default parameters')

        scores = fcast.stats_functions.do_met_stats(None, optimised_wsa, omni_ref_0)
        plt.plot(thresholds, scores, label = f'{make_nicetitle(i)}, optimised parameters')

        plt.legend()
        plt.title('Raw persistence skill scores')

        plt.tight_layout()
        plt.show()

if plot_type == -1 or plot_type == 3: #Do 'persistence metric' or equivalent, for a variety of thresholds.

    run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

    doruns = np.arange(8)
    batch_names = []
    for i in range(8):
        batch_names.append(f'wsa_nocmes_{i}')

    fig = plt.figure(figsize=(12,6))

    for i in doruns:

        batch_name = batch_names[i]
        wsa_fname = f'./data/raw_speeds/{batch_name}_wsa_speeds.txt'
        optimised_fname = f'./data/raw_speeds/{batch_name}_wsa_scaled_speeds.txt'
        ref_fname_0 = f'./data/raw_speeds/{batch_name}_wsa_speeds_ref.txt'
        ref_fname_1 = f'./data/raw_speeds/{batch_name}_wsa_scaled_speeds_ref.txt'

        timeseries = np.loadtxt(f'./data/raw_speeds/{batch_name}_wsa_scaled_times.txt', dtype='datetime64[s]', delimiter = ',')

        if not(os.path.exists(wsa_fname) and os.path.exists(optimised_fname) and os.path.exists(ref_fname_0) and os.path.exists(ref_fname_1)):
            continue

        wsa = np.loadtxt(wsa_fname, delimiter = ',')
        optimised_wsa = np.loadtxt(optimised_fname, delimiter = ',')
        omni_ref_0 = np.loadtxt(ref_fname_0, delimiter = ',')
        omni_ref_1 = np.loadtxt(ref_fname_1, delimiter = ',')

        if filter_for_cmes:
            cme_mask = fcast.data_functions.get_cme_times(timeseries, also_filter_persistence=False)
            invalid_times = np.where(cme_mask == 1)[0]
            wsa[invalid_times] = np.nan
            optimised_wsa[invalid_times] = np.nan
            omni_ref_0[invalid_times] = np.nan
            omni_ref_1[invalid_times] = np.nan

        if True:  #Apply scaling as determined by "scale_for_persistence.py"
            if unoptimised_scales[i] is not None and optimised_scales[i] is not None:
                wsa = scale_function(unoptimised_scales[i][0], unoptimised_scales[i][1], wsa)
                optimised_wsa = scale_function(optimised_scales[i][0], optimised_scales[i][1], optimised_wsa)

        thresholds = np.arange(450,600,5)
        scores = fcast.stats_functions.do_met_stats(None, wsa, omni_ref_0)

        plt.plot(thresholds, scores, label = f'{make_nicetitle(i)}, default parameters')

        scores = fcast.stats_functions.do_met_stats(None, optimised_wsa, omni_ref_0)
        #plt.plot(thresholds, scores, label = f'{make_nicetitle(i)}')#, optimised parameters')

    plt.ylim(-0.1,0.4)
    plt.axhline(c='black', linestyle='dashed')
    plt.legend(fontsize=10)
    #plt.title('Skill scores by threshold, optimised for distributions')
    plt.title('Skill scores by threshold, default scaled WSA')

    plt.tight_layout()
    plt.savefig('./plots/skillscores_wsa_scaled.png')
    plt.show()
