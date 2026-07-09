#This has now become (as predicted) an unweildy mess.
#I'm going to attempt to write what all the options and processes are. Bear with.

"""
There are two stages for processing the data.

The first, which we shall call 'optimisation', is when the parameters in the WSA model (or potentially a neural net) are modified to fit either a distribution, rms, or something else entirely.

The second, which we shall call 'scaling', is how the output velocities can be scaled to minimise something else. This (so far) can be rms, skill score, or potentially distributions.

Today, I'd like to generalise this so it works nicely. Alas the 'unscaled' raw wsa speeds don't follow the same pattern as everything else, for reasons. Perhaps we should change that.

Either way, the data can all be read in AT THE START, and cehcekd for consistency etc., before various plots can be made resulting from them.

I'm adding to this file to look into comparisons with the persistence model. If it is well-correlated, then it might be really good to take a combination of persistence PLUS any information given by the VSW predictions (difference in one month to the next, perhaps, and scaled?)
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import wind_forecast as fcast
from scipy.ndimage import gaussian_filter1d
from datetime import datetime, timedelta
from scipy.stats import pearsonr
from scipy.optimize import minimize

parameter_sources = ["raw", "wsa_nocmes", "rms_nocmes", "corr_nocmes"]
scale_sources = ["raw", "ss", "rms", "dist", "ss_raw"]

for a in [0]:
    bs = [0]
    #bs = [0,2,3]
    for b in bs:

        #Pick the desired combination here. The titles above should be kept consistent, but can obviously be added to if desired.
        #This code should be copied into 'scale_velocites' for consistency and neatness.
        parameter_select = a; scale_select = b

        parameter_source = parameter_sources[parameter_select]
        scale_source = scale_sources[scale_select]

        parameter_shortnames = ["raw", "dist", "rms", "corr"]
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
            if parameter_source == "raw":
                root = "Default WSA Parameters"
            elif parameter_source == "wsa_nocmes":
                root = "WSA Optimised for Speed Distributions"
            elif parameter_source == "rms_nocmes":
                root = "WSA Optimised for RMS"
            elif parameter_source == "corr_nocmes":
                root = "WSA Optimised for Correlation"
            else:
                raise Exception("Parameter source not regonised")

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

        omni_fname = f'./data/raw_speeds/wsa_nocmes_0_wsa_speeds_ref.txt'

        persistence_time = 27.27 #Time in days for persistence model checks
        persistence_time_int = int(persistence_time*24)  #The cut in hours
        omni = np.loadtxt(omni_fname, delimiter = ',')
        persist_data = np.nan*omni
        persist_data[persistence_time_int:] = omni[:-persistence_time_int]

        timeseries = np.loadtxt(f'./data/raw_speeds/wsa_nocmes_0_wsa_scaled_times.txt', dtype='datetime64[s]', delimiter = ',')  #Just pick any of the timeseries -- they all should match

        # time_windows = np.linspace(0,10,25)  #This is for testing the different CME filtering timescales

        # xs = []; ys = []; y2s = []
        # for window in time_windows:
        #     xs.append(window); ys.append(r)
        #     print(window, r)
        #
        # plt.plot(xs, ys/np.max(ys))
        # plt.plot(xs, y2s/np.max(y2s))
        # plt.show()

        cme_mask = fcast.data_functions.get_cme_times(timeseries)

        invalid_times = np.where(cme_mask == 1)[0]

        persist_data_filtered = persist_data.copy()
        omni_filtered = omni.copy()
        persist_data_filtered[invalid_times] = np.nan
        omni_filtered[invalid_times] = np.nan

        data_lengths = [len(persist_data_filtered), len(omni_filtered)]

        if not min(data_lengths) == max(data_lengths):
            raise Exception("Timeseries data lengths don't match. Not sure what to do...")

        nas = np.logical_or(np.isnan(omni_filtered), np.isnan(persist_data_filtered))

        r, _ = pearsonr(omni_filtered[~nas], persist_data_filtered[~nas])

        cmap, xs, ys = fcast.stats_functions.find_data_colourmap(omni_filtered[~nas], persist_data_filtered[~nas], 300)

        # plt.pcolormesh(xs, ys, cmap.T, vmax=np.percentile(cmap,99.5))
        # plt.show()

        print('Persistence correlation:', r)
        #Actually finally, plot the predicted speeds against the actual ones, and do some correlations (maybe)
        if plot_type == -1 or plot_type == 3:

            run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

            batch_names = []

            for i in range(8):
                batch_names.append(f'{parameter_source}_{i}')

            fig1, axs1 = plt.subplots(2,4, figsize=(12,6))
            for i, batch_name in enumerate(batch_names[:]):

                #Hopefully all things should be arranged nicely time-wise, but do need to check as much
                if parameter_source == "raw":
                    data_fname = f'./data/raw_speeds/wsa_nocmes_{i}_wsa_speeds.txt'

                else:
                    data_fname = f'./data/raw_speeds/{batch_name}_wsa_scaled_speeds.txt'

                if not(os.path.exists(omni_fname) and os.path.exists(data_fname)):
                    print('Files not found...', omni_fname, data_fname)
                    continue

                wsa = np.loadtxt(data_fname, delimiter = ',')


                if parameter_source == "raw":
                    timeseries = np.loadtxt(f'./data/raw_speeds/wsa_nocmes_{i}_wsa_scaled_times.txt', dtype='datetime64[s]', delimiter = ',')
                else:
                    timeseries = np.loadtxt(f'./data/raw_speeds/{batch_name}_wsa_scaled_times.txt', dtype='datetime64[s]', delimiter = ',')

                omni_shift = np.nan*omni
                omni_shift[persistence_time_int:] = omni[:-persistence_time_int]

                wsa_shift = np.nan*omni
                wsa_shift[persistence_time_int:] = wsa[:-persistence_time_int]

                cme_mask = fcast.data_functions.get_cme_times(timeseries)
                invalid_times = np.where(cme_mask == 1)[0]

                wsa_filtered = wsa.copy()
                omni_filtered = omni.copy()
                wsa_filtered[invalid_times] = np.nan
                omni_filtered[invalid_times] = np.nan

                wsa_shift_filtered = wsa_shift.copy()
                omni_shift_filtered = omni_shift.copy()
                wsa_shift_filtered[invalid_times] = np.nan
                omni_shift_filtered[invalid_times] = np.nan

                if scales[i] is not None:
                    wsa = scale_function(scales[i][0], scales[i][1], wsa)
                    wsa_filtered = scale_function(scales[i][0], scales[i][1], wsa_filtered)

                data_lengths = [len(wsa), len(omni)]

                if not min(data_lengths) == max(data_lengths):
                    raise Exception("Timeseries data lengths don't match. Not sure what to do...")

                wsa_diff = wsa_filtered - wsa_shift_filtered   #Difference in WSA prediction since the last month

                xdata = omni_shift_filtered
                ydata = omni_filtered

                nas = np.logical_or(np.isnan(xdata), np.isnan(ydata))
                r, _ = pearsonr(xdata[~nas], ydata[~nas])
                rms = np.sqrt(np.nanmean((xdata - ydata)**2))

                #print('Persistence correlation and RMS:',r, rms)

                xdata = wsa_filtered
                ydata = omni_filtered

                nas = np.logical_or(np.isnan(xdata), np.isnan(ydata))
                r, _ = pearsonr(xdata[~nas], ydata[~nas])
                rms = np.sqrt(np.nanmean((xdata - ydata)**2))
                #print('Raw HuxT correlation and RMS:',r, rms)

                def check_combination_factor(factor):
                    #Gives 1-the correation, as we want something that can be minimised
                    combined_metric = omni_shift_filtered + factor*wsa_diff
                    nas = np.logical_or(np.isnan(combined_metric), np.isnan(omni_filtered))
                    r, _ = pearsonr(combined_metric[~nas], omni_filtered[~nas])
                    return 1.0-r

                optimum_factor = minimize(check_combination_factor, x0 = 1.0)
                #print('Optimum factor', optimum_factor.x)

                best_metric = omni_shift_filtered + optimum_factor.x*wsa_diff
                xdata = best_metric
                ydata = omni_filtered

                nas = np.logical_or(np.isnan(xdata), np.isnan(ydata))
                r, _ = pearsonr(xdata[~nas], ydata[~nas])
                rms = np.sqrt(np.nanmean((xdata - ydata)**2))
                print('Combined correlation and RMS:',r, rms)

                ax = axs1[i//4, i%4]

                cmap, xs, ys = fcast.stats_functions.find_data_colourmap(best_metric[~nas], omni_filtered[~nas], 300, xmin=200, xmax=800, ymin=200, ymax=800)

                ax.pcolormesh(xs, ys, cmap.T, vmax=np.percentile(cmap,99.5))

                #ax.scatter(wsa_filtered, omni_filtered, c = 'black', s = 0.1)
                #ax.plot(hist_ref, c = 'black', linestyle = 'dashed')
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xlim(200,800)
                ax.set_ylim(200,800)
                ax.set_title(f"{make_nicetitle(i)}, r = {r:.3f}", fontsize = 10)
                #ax.set_ylim(-0.005,0.15)

            plt.suptitle(f"{suptitle}")

            plt.tight_layout()
            plt.savefig(f'./plots/combined_test.png')

            plt.show()
