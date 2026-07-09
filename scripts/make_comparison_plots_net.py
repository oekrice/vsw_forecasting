#This has now become (as predicted) an unweildy mess.
#I'm going to attempt to write what all the options and processes are. Bear with.

"""
There are two stages for processing the data.

The first, which we shall call 'optimisation', is when the parameters in the WSA model (or potentially a neural net) are modified to fit either a distribution, rms, or something else entirely.

The second, which we shall call 'scaling', is how the output velocities can be scaled to minimise something else. This (so far) can be rms, skill score, or potentially distributions.

Today, I'd like to generalise this so it works nicely. Alas the 'unscaled' raw wsa speeds don't follow the same pattern as everything else, for reasons. Perhaps we should change that.

Either way, the data can all be read in AT THE START, and cehcekd for consistency etc., before various plots can be made resulting from them.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import wind_forecast as fcast
from scipy.ndimage import gaussian_filter1d
from datetime import datetime, timedelta
from scipy.stats import pearsonr

parameter_sources = ["net_test"]
scale_sources = ["raw", "ss", "rms", "dist", "ss_raw"]

for a in range(1):
    #bs = [0,1,2,3,4]
    bs = [0]
    for b in bs:

        #Pick the desired combination here. The titles above should be kept consistent, but can obviously be added to if desired.
        #This code should be copied into 'scale_velocites' for consistency and neatness.
        parameter_select = a; scale_select = b

        parameter_source = parameter_sources[parameter_select]
        scale_source = scale_sources[scale_select]

        parameter_shortnames = ["net"]
        parameter_shortname = parameter_shortnames[parameter_select]


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

        scales = []

        for i in range(8):
            if scale_source != "raw":
                if os.path.exists(f'./data/scaling_data/{parameter_shortname}_{scale_source}_{i}.txt'):
                    scale_data = np.loadtxt(f'./data/scaling_data/{parameter_shortname}_{scale_source}_{i}.txt')
                    scales.append(scale_data[1:])
                else:
                    raise Exception("Scaling source not found. Try running the scaling script for this combination of sources/parameters.")
            else:
                scales.append(None)

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

        def make_suptitle(parameter_source, scale_source):
            root = "net_test"

            if scale_source == "raw":
                root += ""
            elif parameter_source != "raw":
                root += " and "
            else:
                root += ", "

            if scale_source == "raw":
                end = ""
            elif scale_source == "ss":
                end = "Scaled for Persistence Skill Score"
            elif scale_source == "rms":
                end = "Scaled for RMS"
            elif scale_source == "dist":
                end = "Scaled for Speed Distributions"
            elif scale_source == "ss_raw":
                end = "Scaled for Raw Skill Score"
            else:
                raise Exception("Scale source not regonised")

            return root+end

        suptitle = make_suptitle(parameter_source, scale_source)

        print('Doing plots for:', suptitle)

        if plot_type == -1 or plot_type == 0: #Do timeseries and print out RMS values. Alas these appear to be consistently worse once optimised. Bugger. Yes.
            batch_names = []

            for i in range(8):
                batch_names.append(f'{parameter_source}_{i}')

            omni_fname = f'./data/raw_speeds/net_test_7_neural_net_speeds_ref.txt'
            fig = plt.figure(figsize=(12,6))
            for i, batch_name in enumerate(batch_names):

                #Hopefully all things should be arranged nicely time-wise, but do need to check as much
                if parameter_source == "raw":
                    data_fname = f'./data/raw_speeds/wsa_nocmes_{i}_neural_net_speeds.txt'

                else:
                    data_fname = f'./data/raw_speeds/{batch_name}_neural_net_speeds.txt'

                if not(os.path.exists(omni_fname) and os.path.exists(data_fname)):
                    print('Files not found...', omni_fname, data_fname)
                    continue

                wsa = np.loadtxt(data_fname, delimiter = ',')
                omni = np.loadtxt(omni_fname, delimiter = ',')

                if parameter_source == "raw":
                    timeseries = np.loadtxt(f'./data/raw_speeds/wsa_nocmes_{i}_neural_net_times.txt', dtype='datetime64[s]', delimiter = ',')
                else:
                    timeseries = np.loadtxt(f'./data/raw_speeds/{batch_name}_neural_net_times.txt', dtype='datetime64[s]', delimiter = ',')

                "Always filter for CMEs, but save this data separately"
                cme_mask = fcast.data_functions.get_cme_times(timeseries)
                invalid_times = np.where(cme_mask == 1)[0]

                wsa_filtered = wsa.copy()
                omni_filtered = omni.copy()
                wsa_filtered[invalid_times] = np.nan
                omni_filtered[invalid_times] = np.nan

                if scales[i] is not None:
                    wsa = scale_function(scales[i][0], scales[i][1], wsa)
                    wsa_filtered = scale_function(scales[i][0], scales[i][1], wsa_filtered)

                data_lengths = [len(wsa), len(omni)]

                if not min(data_lengths) == max(data_lengths):
                    raise Exception("Timeseries data lengths don't match. Not sure what to do...")

                abs_difference = np.abs(wsa-omni)
                #Split this up into some time chunks and average. Otherwise is quite a mess. Current timeseries is every day.
                #Need to do this and THEN get the RMS values once CMEs are taken into account
                difference_cadence = 30*24 #Number of hours to split
                #Plot the differences NOT removing the CMEs

                time_slices = []
                diff_slices = []
                i_min = 0; i_max = difference_cadence
                while i_max < len(abs_difference):
                    time_slices.append(timeseries[i_min] + 0.5*(timeseries[i_max]-timeseries[i_min]))
                    diff_slices.append(np.nanmean(abs_difference[i_min:i_max]))
                    i_min = i_min + difference_cadence
                    i_max = i_max + difference_cadence

                rms = np.sqrt(np.nanmean((wsa_filtered - omni)**2))

                # print('Standard STD for ', batch_name, np.sqrt(np.nanmean((wsa-omni_ref_0)**2)))
                # print('Optimised STD for', batch_name, np.sqrt(np.nanmean((optimised_neural_net-omni_ref_1)**2)))
                plt.plot(time_slices, diff_slices, label = f'{make_nicetitle(i)}, rms = {rms:.0f}km/s')

            plt.legend(fontsize=10)
            plt.title(suptitle)
            plt.ylim(0,300)
            plt.tight_layout()
            plt.savefig(f'./plots/data_comparison/errors_{parameter_shortname}_{scale_source}.png')
            plt.close()

        #Now the histograms of distributions.

        if plot_type == -1 or plot_type == 1: #Do timeseries and print out RMS values. Alas these appear to be consistently worse once optimised. Bugger. Yes.
            batch_names = []

            for i in range(8):
                batch_names.append(f'{parameter_source}_{i}')

            omni_fname = f'./data/raw_speeds/net_test_7_neural_net_speeds_ref.txt'
            fig1, axs1 = plt.subplots(2,4, figsize=(12,6))

            for i, batch_name in enumerate(batch_names):

                #Hopefully all things should be arranged nicely time-wise, but do need to check as much
                if parameter_source == "raw":
                    data_fname = f'./data/raw_speeds/wsa_nocmes_{i}_neural_net_speeds.txt'

                else:
                    data_fname = f'./data/raw_speeds/{batch_name}_neural_net_speeds.txt'

                if not(os.path.exists(omni_fname) and os.path.exists(data_fname)):
                    print('Files not found...', omni_fname, data_fname)
                    continue

                wsa = np.loadtxt(data_fname, delimiter = ',')
                omni = np.loadtxt(omni_fname, delimiter = ',')

                if parameter_source == "raw":
                    timeseries = np.loadtxt(f'./data/raw_speeds/wsa_nocmes_{i}_neural_net_times.txt', dtype='datetime64[s]', delimiter = ',')
                else:
                    timeseries = np.loadtxt(f'./data/raw_speeds/{batch_name}_neural_net_times.txt', dtype='datetime64[s]', delimiter = ',')

                "Always filter for CMEs, but save this data separately"
                cme_mask = fcast.data_functions.get_cme_times(timeseries)
                invalid_times = np.where(cme_mask == 1)[0]

                wsa_filtered = wsa.copy()
                omni_filtered = omni.copy()
                wsa_filtered[invalid_times] = np.nan
                omni_filtered[invalid_times] = np.nan

                if scales[i] is not None:
                    wsa = scale_function(scales[i][0], scales[i][1], wsa)
                    wsa_filtered = scale_function(scales[i][0], scales[i][1], wsa_filtered)

                data_lengths = [len(wsa), len(omni)]

                if not min(data_lengths) == max(data_lengths):
                    raise Exception("Timeseries data lengths don't match. Not sure what to do...")

                abs_difference = np.abs(wsa-omni)
                #Split this up into some time chunks and average. Otherwise is quite a mess. Current timeseries is every day.
                #Need to do this and THEN get the RMS values once CMEs are taken into account
                difference_cadence = 30*24 #Number of hours to split
                #Plot the differences NOT removing the CMEs

                time_slices = []
                diff_slices = []
                i_min = 0; i_max = difference_cadence
                while i_max < len(abs_difference):
                    time_slices.append(timeseries[i_min] + 0.5*(timeseries[i_max]-timeseries[i_min]))
                    diff_slices.append(np.nanmean(abs_difference[i_min:i_max]))
                    i_min = i_min + difference_cadence
                    i_max = i_max + difference_cadence

                rms = np.sqrt(np.nanmean((wsa_filtered - omni)**2))

                nbins = 101
                hist_neural_net, _ = np.histogram(wsa_filtered, bins=nbins, range=(0.0,1000.0))
                hist_ref, _ = np.histogram(omni_filtered, bins=nbins, range=(0.0,1000.0))

                hist_neural_net = hist_neural_net/np.sum(hist_neural_net)
                hist_ref = hist_ref/np.sum(hist_ref)
                distance = np.sum((hist_neural_net - hist_ref)**2)

                print('Histogram Distance:', distance)

                ax = axs1[i//4, i%4]
                ax.plot(hist_neural_net)
                ax.plot(hist_ref, c = 'black', linestyle = 'dashed')
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_title(make_nicetitle(i), fontsize = 10)
                ax.set_ylim(-0.005,0.15)

            plt.suptitle(suptitle)
            plt.tight_layout()
            plt.savefig(f'./plots/data_comparison/distributions_{parameter_shortname}_{scale_source}.png')

            plt.close()

        #Finally do the skill scores along the lines of above. One plot is probably sufficient. Unless we do it by year? Could be complicated though.
        if plot_type == -1 or plot_type == 2: #Do 'persistence metric' or equivalent, for a variety of thresholds.

            run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

            batch_names = []

            for i in range(8):
                batch_names.append(f'{parameter_source}_{i}')

            omni_fname = f'./data/raw_speeds/net_test_7_neural_net_speeds_ref.txt'
            fig = plt.figure(figsize=(12,6))
            for i, batch_name in enumerate(batch_names):

                #Hopefully all things should be arranged nicely time-wise, but do need to check as much
                if parameter_source == "raw":
                    data_fname = f'./data/raw_speeds/wsa_nocmes_{i}_neural_net_speeds.txt'

                else:
                    data_fname = f'./data/raw_speeds/{batch_name}_neural_net_speeds.txt'

                if not(os.path.exists(omni_fname) and os.path.exists(data_fname)):
                    print('Files not found...', omni_fname, data_fname)
                    continue

                wsa = np.loadtxt(data_fname, delimiter = ',')
                omni = np.loadtxt(omni_fname, delimiter = ',')

                if parameter_source == "raw":
                    timeseries = np.loadtxt(f'./data/raw_speeds/wsa_nocmes_{i}_neural_net_times.txt', dtype='datetime64[s]', delimiter = ',')
                else:
                    timeseries = np.loadtxt(f'./data/raw_speeds/{batch_name}_neural_net_times.txt', dtype='datetime64[s]', delimiter = ',')

                "Always filter for CMEs, but save this data separately"
                cme_mask = fcast.data_functions.get_cme_times(timeseries)
                invalid_times = np.where(cme_mask == 1)[0]

                wsa_filtered = wsa.copy()
                omni_filtered = omni.copy()
                wsa_filtered[invalid_times] = np.nan
                omni_filtered[invalid_times] = np.nan

                if scales[i] is not None:
                    wsa = scale_function(scales[i][0], scales[i][1], wsa)
                    wsa_filtered = scale_function(scales[i][0], scales[i][1], wsa_filtered)

                data_lengths = [len(wsa), len(omni)]

                if not min(data_lengths) == max(data_lengths):
                    raise Exception("Timeseries data lengths don't match. Not sure what to do...")


                thresholds = np.arange(450,600,5)
                scores = fcast.stats_functions.do_met_stats(None, wsa_filtered, omni_filtered)

                plt.plot(thresholds, scores, label = f'{make_nicetitle(i)}')

            plt.axhline(0.0, c = 'black', linestyle='dashed')
            plt.ylim(-0.35,0.35)
            plt.legend(fontsize = 10)
            plt.suptitle(f"Skill Scores, {suptitle}")

            plt.tight_layout()
            plt.savefig(f'./plots/data_comparison/skillscores_{parameter_shortname}_{scale_source}.png')

            plt.close()


        #Actually finally, plot the predicted speeds against the actual ones, and do some correlations (maybe)
        if plot_type == -1 or plot_type == 3:

            run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

            batch_names = []

            for i in range(8):
                batch_names.append(f'{parameter_source}_{i}')

            omni_fname = f'./data/raw_speeds/net_test_7_neural_net_speeds_ref.txt'
            fig1, axs1 = plt.subplots(2,4, figsize=(12,6))
            for i, batch_name in enumerate(batch_names):

                #Hopefully all things should be arranged nicely time-wise, but do need to check as much
                if parameter_source == "raw":
                    data_fname = f'./data/raw_speeds/wsa_nocmes_{i}_neural_net_speeds.txt'

                else:
                    data_fname = f'./data/raw_speeds/{batch_name}_neural_net_speeds.txt'

                if not(os.path.exists(omni_fname) and os.path.exists(data_fname)):
                    print('Files not found...', omni_fname, data_fname)
                    continue

                wsa = np.loadtxt(data_fname, delimiter = ',')
                omni = np.loadtxt(omni_fname, delimiter = ',')

                if parameter_source == "raw":
                    timeseries = np.loadtxt(f'./data/raw_speeds/wsa_nocmes_{i}_neural_net_times.txt', dtype='datetime64[s]', delimiter = ',')
                else:
                    timeseries = np.loadtxt(f'./data/raw_speeds/{batch_name}_neural_net_times.txt', dtype='datetime64[s]', delimiter = ',')

                "Always filter for CMEs, but save this data separately"
                cme_mask = fcast.data_functions.get_cme_times(timeseries)
                invalid_times = np.where(cme_mask == 1)[0]

                wsa_filtered = wsa.copy()
                omni_filtered = omni.copy()
                wsa_filtered[invalid_times] = np.nan
                omni_filtered[invalid_times] = np.nan

                if scales[i] is not None:
                    wsa = scale_function(scales[i][0], scales[i][1], wsa)
                    wsa_filtered = scale_function(scales[i][0], scales[i][1], wsa_filtered)

                data_lengths = [len(wsa), len(omni)]

                if not min(data_lengths) == max(data_lengths):
                    raise Exception("Timeseries data lengths don't match. Not sure what to do...")

                nas = np.logical_or(np.isnan(wsa_filtered), np.isnan(omni_filtered))
                r, _ = pearsonr(wsa_filtered[~nas], omni_filtered[~nas])
                ax = axs1[i//4, i%4]
                ax.scatter(wsa_filtered, omni_filtered, c = 'black', s = 0.1)
                #ax.plot(hist_ref, c = 'black', linestyle = 'dashed')
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xlim(200,800)
                ax.set_ylim(200,800)
                ax.set_title(f"{make_nicetitle(i)}, r = {r:.3f}", fontsize = 10)
                #ax.set_ylim(-0.005,0.15)

            plt.suptitle(f"{suptitle}")

            plt.tight_layout()
            plt.savefig(f'./plots/data_comparison/correlation_{parameter_shortname}_{scale_source}.png')

            plt.close()
