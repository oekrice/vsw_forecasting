#This script requires test_predictions to have been run beforehand. Just loads in text files of speeds, which are hopefully aligned, and calculates/plots various things
#Seems to be best to do things separately like this to allow for a lot of flexibility.
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import wind_forecast as fcast

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

unoptimised_scales = [None, None, None, None, None, None, None, [1.08, -0.156]]
optimised_scales = [None, None, None, None, None, None, None, [ 0.8411262,  -0.05513745]]

def make_nicetitle(id):
    if (id//2)%2 == 0:
        is_pfss = True
        model = "pfss"
    else:
        is_pfss = False
        model = "outflow"
    if (id%2) == 0:
        rss = 2.5
    else:
        rss = 5.0
    if (id//4) == 0:
        source = "gong"
    else:
        source = "hmi"

    nicetitle = f"{model}, rss = {rss}, source = {source}"
    return nicetitle

if plot_type == -1 or plot_type == 0: #Do timeseries and print out RMS values. Alas these appear to be consistently worse once optimised. Bugger.
    batch_names = []
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


        if filter_for_cmes:
            cme_mask = fcast.data_functions.get_cme_times(timeseries)
            invalid_times = np.where(cme_mask == 1)[0]
            wsa[invalid_times] = np.nan
            optimised_wsa[invalid_times] = np.nan
            omni_ref_0[invalid_times] = np.nan
            omni_ref_1[invalid_times] = np.nan

        data_lengths = [len(wsa), len(optimised_wsa), len(omni_ref_0), len(omni_ref_1)]

        if not min(data_lengths) == max(data_lengths):
            raise Exception("Timeseries data lengths don't match. Not sure what to do...")

        plt.plot(timeseries, wsa, c = cmap(i), linestyle = 'dashed')
        plt.plot(timeseries, optimised_wsa, c = cmap(i), linestyle='solid')
        plt.plot(timeseries, omni_ref_1, c = 'black')

        print('Standard STD for ', batch_name, np.sqrt(np.nanmean((wsa-omni_ref_0)**2)))
        print('Optimised STD for', batch_name, np.sqrt(np.nanmean((optimised_wsa-omni_ref_1)**2)))
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

        if False:  #Apply scaling as determined by scale_for_persistence.py
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
