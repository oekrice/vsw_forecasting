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


if len(sys.argv) > 1:
    parameter_set = int(sys.argv[1])
else:
    raise Exception('Specify parameter set.')

parameter_sources = ["raw", "wsa_nocmes", "rms_nocmes", "corr_nocmes"]
scale_sources = ["raw", "ss", "rms", "dist", "ss_raw"]

comparison_type = 0   #Types are 0, 1,2,3
# 0 = Raw WSA
# 1 = Direct combination
# 2 = Optimised combination
# 3 = OMNI Reference

type_titles = ["Raw WSA", "Direct Persistence Combination", "Optimised Persistence Combination", "OMNI Reference"]

for a in [0]:
    bs = [0]
    #bs = [0,2,3]
    for b in bs:

        #Pick the desired combination here. The titles above should be kept consistent, but can obviously be added to if desired.
        #This code should be copied into 'scale_velocites' for consistency and neatness.
        parameter_select = a; scale_select = b

        parameter_source = parameter_sources[parameter_select]
        scale_source = scale_sources[scale_select]

        parameter_shortnames = ["Set 1", "Set 2", "Set 3"]
        parameter_shortname = parameter_shortnames[parameter_select]

        filter_for_cmes = True

        cmap = plt.get_cmap("tab10")

        #Load in the scaling (generated using scale_for_persistence.py)

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


        print('Doing plots for:', parameter_shortnames[parameter_select])

        omni_fname = f'./data/raw_speeds/wsa_nocmes_0_wsa_speeds_ref.txt'

        persistence_time = 27.27 #Time in days for persistence model checks
        persistence_time_int = int(persistence_time*24)  #The cut in hours
        omni = np.loadtxt(omni_fname, delimiter = ',')
        persist_data = np.nan*omni
        persist_data[persistence_time_int:] = omni[:-persistence_time_int]

        timeseries = np.loadtxt(f'./data/raw_speeds/wsa_nocmes_0_wsa_scaled_times.txt', dtype='datetime64[s]', delimiter = ',')  #Just pick any of the timeseries -- they all should match

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
        rms = np.sqrt(np.nanmean((omni_filtered[~nas] - persist_data_filtered[~nas])**2))
        mae = np.nanmean(np.abs(omni_filtered[~nas]-persist_data_filtered[~nas]))
        diffsim, _ = fcast.stats_functions.get_distribution_similarity(omni_filtered[~nas], persist_data_filtered[~nas], iteration=0, doplots=False)

        #Actually finally, plot the predicted speeds against the actual ones, and do some correlations (maybe)

        run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

        batch_names = []

        for i in range(8):
            batch_names.append(f'{parameter_source}_{i}')

        fig1, axs1 = plt.subplots(2,4, figsize=(12,6))
        for i, batch_name in enumerate(batch_names[:]):

            data_fname = f'./data/raw_speeds/combine_raw_{parameter_set}_{i}_wsa_speeds.txt'

            if not(os.path.exists(omni_fname) and os.path.exists(data_fname)):
                print('Files not found...', omni_fname, data_fname)
                continue

            wsa = np.loadtxt(data_fname, delimiter = ',')

            timeseries = np.loadtxt(f'./data/raw_speeds/wsa_nocmes_0_wsa_times_ref.txt', dtype='datetime64[s]', delimiter = ',')

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

            def check_combination_factor_simple(factor):
                #Gives 1-the correation, as we want something that can be minimised
                combined_metric = omni_shift_filtered + factor*wsa_diff
                nas = np.logical_or(np.isnan(combined_metric), np.isnan(omni_filtered))
                r, _ = pearsonr(combined_metric[~nas], omni_filtered[~nas])
                return 1.0-r

            def check_combination_factor(factor):
                #Gives the geometic mean of the two things that can be minimised.
                combined_metric = omni_shift_filtered + factor*wsa_diff
                nas = np.logical_or(np.isnan(combined_metric), np.isnan(omni_filtered))
                r, _ = pearsonr(combined_metric[~nas], omni_filtered[~nas])
                rms = np.sqrt(np.nanmean((combined_metric[~nas] - omni_filtered[~nas])**2))
                ref_r, _ = pearsonr(omni_shift_filtered[~nas], omni_filtered[~nas])
                ref_rms = np.sqrt(np.nanmean((omni_shift_filtered[~nas] - omni_filtered[~nas])**2))

                correlation_improvement = (1.0-r)/(1.0-ref_r)
                rms_improvement = rms/ref_rms

                if correlation_improvement >= 1.0:
                    #Introduce a harsh penalty for getting worse. Needs to be continuous though.
                    correlation_improvement = correlation_improvement + 10*(correlation_improvement - 1.0)
                if rms_improvement >= 1.0:
                    rms_improvement = rms_improvement + 10*(rms_improvement - 1.0)


                res = np.sqrt(correlation_improvement*rms_improvement)

                #print('r/reference, rms/reference. result', r, ref_r, rms, ref_rms, res)

                return res

            def check_raw_factors(factor):
                #Scales the wsa result to minimise the RMS/Correlation combination (will be shifted such that the mean is correct?)

                combined_metric = factor*wsa_filtered

                combined_metric += np.nanmean(omni_filtered - combined_metric)

                nas = np.logical_or(np.logical_or(np.isnan(combined_metric), np.isnan(omni_filtered)), np.isnan(omni_shift_filtered))

                r, _ = pearsonr(combined_metric[~nas], omni_filtered[~nas])
                rms = np.sqrt(np.nanmean((combined_metric[~nas] - omni_filtered[~nas])**2))

                ref_r, _ = pearsonr(omni_shift_filtered[~nas], omni_filtered[~nas])
                ref_rms = np.sqrt(np.nanmean((omni_shift_filtered[~nas] - omni_filtered[~nas])**2))

                correlation_improvement = (1.0-r)/(1.0-ref_r)
                rms_improvement = rms/ref_rms

                if correlation_improvement >= 1.0:
                    #Introduce a harsh penalty for getting worse. Needs to be continuous though.
                    correlation_improvement = correlation_improvement + 10*(correlation_improvement - 1.0)
                if rms_improvement >= 1.0:
                    rms_improvement = rms_improvement + 10*(rms_improvement - 1.0)

                #res = np.sqrt(correlation_improvement*rms_improvement)
                res = rms_improvement
                #print('r/reference, rms/reference. result', r, ref_r, rms, ref_rms, res)

                return res

            #optimum_factor = minimize(check_combination_factor, x0 = 1.0)
            #best_metric = omni_shift_filtered + optimum_factor.x*wsa_diff

            optimum_factor = minimize(check_raw_factors, x0 = 1.0)
            best_metric = optimum_factor.x*wsa_filtered
            best_metric += np.nanmean(omni_filtered - best_metric)

            #print('Optimum factor', optimum_factor.x)

            #best_metric =  wsa_filtered#omni_shift_filtered
            xdata = best_metric
            ydata = omni_filtered

            diffsim, _ = fcast.stats_functions.get_distribution_similarity(xdata, ydata, iteration=0, doplots=False)

            nas = np.logical_or(np.isnan(xdata), np.isnan(ydata))
            r, _ = pearsonr(xdata[~nas], ydata[~nas])
            rms = np.sqrt(np.nanmean((xdata - ydata)**2))
            mae = np.nanmean(np.abs(xdata-ydata))

            thresholds = np.arange(450,600,5)
            ss_persist = fcast.stats_functions.do_met_stats(timeseries, model_speeds=xdata, reference_speeds=ydata, compare_to_persist=True, persistence_reference = omni_shift_filtered, persistence_cadence=int(24*27.7), thresholds=thresholds)
            ss_persist = np.mean(ss_persist)
            ss_raw = fcast.stats_functions.do_met_stats(timeseries, model_speeds=xdata, reference_speeds=ydata, compare_to_persist=False, thresholds=thresholds)
            ss_raw = np.mean(ss_raw)

            print('Combined correlation, RMS, MAE, diff, and skillscores:',r, rms, mae, diffsim, ss_persist, ss_raw)

            #print('Combined correlation, RMS, MAE and distribution similarity:',r, rms, mae, diffsim)

            ax = axs1[i//4, i%4]

            cmap, xs, ys = fcast.stats_functions.find_data_colourmap(best_metric[~nas], omni_filtered[~nas], 300, xmin=200, xmax=800, ymin=200, ymax=800)

            ax.pcolormesh(xs, ys, cmap.T, vmax=np.percentile(cmap,99.5))

            #ax.scatter(wsa_filtered, omni_filtered, c = 'black', s = 0.1)
            #ax.plot(hist_ref, c = 'black', linestyle = 'dashed')
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlim(200,800)
            ax.set_ylim(200,800)
            ax.set_title(f"{make_nicetitle(i)} \n r = {r:.3f}, rms = {rms:.0f}, mae = {mae:.0f}, dist = {diffsim:.3f} \n ss_persist = {ss_persist:.2f}, ss_raw = {ss_raw:.2f}", fontsize = 8)
            #ax.set_ylim(-0.005,0.15)

        plt.suptitle(f"Scaled WSA")

        plt.tight_layout()
        plt.savefig(f'./plots/combine_plots/{parameter_set}_{comparison_type}.png')

        plt.close()
