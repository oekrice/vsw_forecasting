#This script requires test_predictions to have been run beforehand. Just loads in text files of speeds, which are hopefully aligned, and calculates/plots various things
#Seems to be best to do things separately like this to allow for a lot of flexibility.
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import wind_forecast as fcast

cmap = plt.get_cmap("tab10")

if len(sys.argv) > 1:
    plot_type = int(sys.argv[1])
else:
    plot_type = -1

if plot_type == -1 or plot_type == 0: #Do timeseries and print out RMS values. Alas these appear to be consistently worse once optimised. Bugger.
    batch_names = []
    for i in range(7):
        batch_names.append(f'optimise_run_{i}')

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

        plt.plot(wsa, c = cmap(i), linestyle = 'dashed')
        plt.plot(optimised_wsa, c = cmap(i), linestyle='solid')
        plt.plot(omni_ref_1, c = 'black')

        print('Standard STD for ', batch_name, np.sqrt(np.nanmean((wsa-omni_ref_0)**2)))
        print('Optimised STD for', batch_name, np.sqrt(np.nanmean((optimised_wsa-omni_ref_1)**2)))
    plt.show()

if plot_type == -1 or plot_type == 1: #Do histogram comparison
    batch_names = []
    for i in range(8):
        batch_names.append(f'optimise_run_{i}')

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

    batch_names = []
    for i in range(7):
        batch_names.append(f'optimise_run_{i}')
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
        thresholds = np.arange(450,600,5)
        scores = fcast.do_met_stats(None, wsa, omni_ref_0)

        plt.plot(thresholds, scores, label = run_names[i])

    plt.legend()
    plt.title('Raw persistence skill scores')

    plt.tight_layout()
    plt.show()
