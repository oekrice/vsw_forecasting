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
import drms
import wind_forecast as fcast
from scipy.ndimage import gaussian_filter1d
from datetime import datetime, timedelta
from scipy.stats import pearsonr
from scipy.optimize import minimize

parameter_sources = ["combine", "combine"]
scale_sources = ["WSA", "combine", "combine_optimised", "reference"]

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

def get_crot_limits():
    #Obtains the times for the Carrington rotations. Bastardised code from elsewhere, but it should do.

    try:
        c = drms.Client()
        #Find the correct Carrington Rotation for this date.
        crot_times_mdi = c.query(('mdi.synoptic_mr_polfil_96m'), key = ["T_START","T_STOP","CAR_ROT"])
        crot_times_hmi = c.query(('hmi.synoptic_mr_polfil_720s'), key = ["T_START","T_STOP","CAR_ROT"])
    except:
        raise Exception("Failed to find the Carrington Rotation database")

    start_times = []
    end_times = []
    centre_times = []
    crot_numbers = []
    for source, crot_times in enumerate([crot_times_mdi, crot_times_hmi]):
        start_times_raw = list(crot_times.pop("T_START"))
        for i in range(len(start_times_raw)):
            if start_times_raw[i][-6:-4] == "60":
                start_times_raw[i] = start_times_raw[i][:-6] + "00" + start_times_raw[i][-4:]
        end_times_raw = list(crot_times.pop("T_STOP"))
        for i in range(len(end_times_raw)):
            if end_times_raw[i][-6:-4] == "60":
                end_times_raw[i] = end_times_raw[i][:-6] + "00" + end_times_raw[i][-4:]

        for si, crot_number in enumerate(crot_times.pop("CAR_ROT")):
            if (crot_number < 2098 and source == 0) or (crot_number >= 2098 and source == 1):
                start_times.append(datetime.strptime(start_times_raw[si].split('_TAI')[0], "%Y.%m.%d_%H:%M:%S"))
                end_times.append(datetime.strptime(end_times_raw[si].split('_TAI')[0], "%Y.%m.%d_%H:%M:%S"))
                centre_times.append(0.5*(start_times[-1] - end_times[-1]) + start_times[-1])
                crot_numbers.append(crot_number)

    return np.array(crot_numbers), np.array(start_times)

#To copy raw speeds from Hamilton

# scp -r vgjn10@hamilton8.dur.ac.uk:/nobackup/vgjn10/projects/vsw_forecasting/paper/data/raw_speeds/combine_* ./paper/data/raw_speeds
crot_numbers, crot_starts = get_crot_limits()  #Just do this once, whatever.

for a in [1]:
    #bs = [0,1,2]
    bs = [0]
    #bs = [0,2,3]
    if a == 0:
        source_title = "Default Parameters"
    elif a == 1:
        source_title = "Optimised Parameters"
    else:
        raise Exception('Source not recognised')
    for b in bs:

        #Pick the desired combination here. The titles above should be kept consistent, but can obviously be added to if desired.
        #This code should be copied into 'scale_velocites' for consistency and neatness.
        parameter_source = parameter_sources[a]
        scale_source = scale_sources[b]

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

        ref_batch_name = f'{parameter_source}_{0}_{0}'
        omni_fname = f'./paper/data/raw_speeds/{ref_batch_name}_wsa_combined_speeds_ref.txt'

        persistence_time = 27.27 #Time in days for persistence model checks
        persistence_time_int = int(persistence_time*24)  #The cut in hours
        omni = np.loadtxt(omni_fname, delimiter = ',')
        persist_data = np.nan*omni
        persist_data[persistence_time_int:] = omni[:-persistence_time_int]

        timeseries = np.loadtxt(f'./paper/data/raw_speeds/{ref_batch_name}_wsa_combined_times.txt', dtype='datetime64[s]', delimiter = ',')  #Just pick any of the timeseries -- they all should match

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


        run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

        batch_names = []

        for i in range(0,8):
            batch_names.append(f'{parameter_source}_{a}_{i}')

        plotted_omni = False
        fig, axs = plt.subplots(2, figsize = (fig_width,0.6*fig_width))
        for i in [4,5,6,7]:

            batch_name = batch_names[i]
            #Hopefully all things should be arranged nicely time-wise, but do need to check as much
            if a == 1:
                print("Using 'scaled' parameters")
                batch_name = f"cma_{i}"
                #Hopefully all things should be arranged nicely time-wise, but do need to check as much
                data_fname = f'./paper/data/raw_speeds/{batch_name}_wsa_combined_speeds.txt'

            else:
                #Hopefully all things should be arranged nicely time-wise, but do need to check as much
                data_fname = f'./paper/data/raw_speeds/{batch_name}_wsa_combined_speeds.txt'

            if not(os.path.exists(omni_fname) and os.path.exists(data_fname)):
                print('Files not found...', omni_fname, data_fname)
                continue

            wsa = np.loadtxt(data_fname, delimiter = ',')


            timeseries = np.loadtxt(f'./paper/data/raw_speeds/{batch_name}_wsa_combined_times.txt', dtype='datetime64[s]', delimiter = ',')

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

            def check_combination_factor(factor, filter_firsthalf=False):
                #Gives the geometic mean of the two things that can be minimised.


                combined_metric = omni_shift_filtered + factor*wsa_diff
                nas = np.logical_or(np.isnan(combined_metric), np.isnan(omni_filtered))

                if filter_firsthalf: #Just use the first half of each year to do the optimisation (roughly, I'm not caring about leap years or anything here...)
                    #print(len(combined_metric), np.sum(nas))
                    #Disallow the second half of each year, so it's not over-optimised
                    cadence = int(24*365.25)
                    n_min = 0; n_max = cadence//2
                    while n_max < len(combined_metric):
                        nas[n_min:n_max] = True
                        n_min += cadence; n_max += cadence
                    #print(len(combined_metric), np.sum(nas))

                r, _ = pearsonr(combined_metric[~nas], omni_filtered[~nas])
                rms = np.sqrt(np.nanmean((combined_metric[~nas] - omni_filtered[~nas])**2))
                ref_r, _ = pearsonr(omni_shift_filtered[~nas], omni_filtered[~nas])
                ref_rms = np.sqrt(np.nanmean((omni_shift_filtered[~nas] - omni_filtered[~nas])**2))

                # correlation_improvement = (1.0-r)/(1.0-ref_r)
                # rms_improvement = rms/ref_rms

                correlation_improvement = (1.0-r)
                rms_improvement = rms/100.0

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

            def do_crot_stats(counter, metric, omni, timeseries, crot_number, crot_times):
                #Will filter the metric further and do rms/correlation on them. Might need to be smarter with the filtering but we'll see
                #Find min and max indices of the times
                start_time = crot_times[counter]; end_time = crot_times[counter + 1]
                min_index = np.searchsorted(timeseries, start_time)
                max_index = np.searchsorted(timeseries, end_time)
                xdata = metric[min_index:max_index]
                time_slice = timeseries[min_index:max_index]
                ydata = omni[min_index:max_index]

                nas = np.logical_or(np.isnan(xdata), np.isnan(ydata))
                if len(xdata) > 24*25:  #Make sure at least half the rotation is logged
                    if np.sum(nas)/len(xdata) < 0.25:  #Too many nans, don't bother with this one (hopefully won't happen that often...)
                        r, _ = pearsonr(xdata[~nas], ydata[~nas])
                        rms = np.sqrt(np.nanmean((xdata - ydata)**2))
                        #print(crot_number, r, rms, np.sum(nas)/len(xdata))
                    else:
                        r = np.nan
                        rms = np.nan
                else:
                    r = np.nan
                    rms = np.nan

                return r, rms

            def do_crot_stats_average(counter, metric, omni, timeseries, crot_number, crot_times):
                #Will filter the metric further and do rms/correlation on them. Might need to be smarter with the filtering but we'll see
                #Find min and max indices of the times
                start_time = crot_times[counter-5]; end_time = crot_times[counter + 6]
                min_index = np.searchsorted(timeseries, start_time)
                max_index = np.searchsorted(timeseries, end_time)
                xdata = metric[min_index:max_index]
                time_slice = timeseries[min_index:max_index]
                ydata = omni[min_index:max_index]

                nas = np.logical_or(np.isnan(xdata), np.isnan(ydata))
                if len(xdata) > 10*27:  #Make sure at least half the rotation is logged? Eh, this'll do.
                    if np.sum(nas)/len(xdata) < 0.75:  #Too many nans, don't bother with this one (hopefully won't happen that often...)
                        r, _ = pearsonr(xdata[~nas], ydata[~nas])
                        rms = np.sqrt(np.nanmean((xdata - ydata)**2))
                        print(crot_number, r, rms, np.sum(nas)/len(xdata))
                    else:
                        r = np.nan
                        rms = np.nan
                else:
                    r = np.nan
                    rms = np.nan

                return r, rms

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

            if scale_source == "WSA":
                best_metric =  wsa_filtered
            elif scale_source == "combine":
                best_metric = omni_shift_filtered + 1.0*wsa_diff
            elif scale_source == "combine_optimised":
                optimum_factor = minimize(check_combination_factor, x0 = 1.0, args=(True))
                best_metric = omni_shift_filtered + optimum_factor.x*wsa_diff
            elif scale_source == "reference":
                best_metric = omni_shift_filtered
            else:
                raise Exception('Scaling source not found')

            start_plot_cut = 188
            end_plot_cut = 30
            if not plotted_omni: #The OMNI Reference
                all_rmss = np.nan*crot_numbers
                all_rs = np.nan*crot_numbers
                for ci, crot in enumerate(crot_numbers[:-6]):
                    r, rms = do_crot_stats_average(ci, omni_shift_filtered, omni_filtered, timeseries, crot, crot_starts)
                    all_rmss[ci] = rms
                    all_rs[ci] = r
                axs[0].plot(crot_numbers[start_plot_cut:-end_plot_cut], all_rmss[start_plot_cut:-end_plot_cut], c = 'black', linestyle='dashed',linewidth=1.0)
                axs[1].plot(crot_numbers[start_plot_cut:-end_plot_cut], all_rs[start_plot_cut:-end_plot_cut], label = 'Persistence', linestyle='dashed', c = 'black', linewidth=1.0)
                plotted_omni=True
            #Best_metric is what we're comparing against, no matter what. Need to now do stats on each Carrington rotation.
            #I'll try to do this without being clever, but if it's too slow might have to be clever. Let's see.
            all_rmss = np.nan*crot_numbers
            all_rs = np.nan*crot_numbers
            for ci, crot in enumerate(crot_numbers[:-6]):
                r, rms = do_crot_stats_average(ci, best_metric, omni_filtered, timeseries, crot, crot_starts)
                all_rmss[ci] = rms
                all_rs[ci] = r

            axs[0].plot(crot_numbers[start_plot_cut:-end_plot_cut], all_rmss[start_plot_cut:-end_plot_cut], linewidth=1.0)
            axs[0].set_ylabel('RMS')
            axs[0].set_xticks([])
            axs[1].plot(crot_numbers[start_plot_cut:-end_plot_cut], all_rs[start_plot_cut:-end_plot_cut], label = make_nicetitle(i), linewidth=1.0)
            axs[1].set_ylabel('Correlation $r$')
            axs[1].set_xlabel('Carrington Rotation')

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


        handles, labels = axs[1].get_legend_handles_labels()
        fig.legend(handles, labels,
           loc="lower center",
           ncol=3,                  # adjust as needed
           bbox_to_anchor=(0.5, -0.0), fontsize=8)

        plt.tight_layout(rect=[0, 0.08, 1, 1])
        plt.savefig(f'./paper/plots/7_optimised_time_plot.pdf')
        plt.show()
        plt.close()
