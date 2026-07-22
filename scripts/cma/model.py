#This is a testbed for making the forecast scripts actually nice, and doing it all properly and things. HA, that went well!

import os
import sys
import numpy as np
import csv
from datetime import datetime, timedelta
#
import wind_forecast as wf  #This should now contain everything we need...
from scipy.stats import pearsonr

def run_model(run_parameters, theta=None, snap_subset=None, iteration=0, output_distributions=False, save_speeds=False, use_old_chb_formula=False):
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

    run_name = run_parameters["base_name"]

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
                if si == 0:
                    vr = wf.field_calculations.compute_vr(snap_id, run_name, run_parameters["velocity_type"], params = theta, doplot=run_parameters["do_plots"], iteration=iteration, huxt_name=run_parameters["run_name"])
                else:
                    vr = wf.field_calculations.compute_vr(snap_id, run_name, run_parameters["velocity_type"], params = theta, doplot=False, iteration=iteration, huxt_name=run_parameters["run_name"])

            else:
                raise Exception("Velocity calculation type not recognised...")

            if run_parameters["verbose"]:
                print('Calculating VSW with HuxT...')

            #Using this velocity map (and nothing else?), run HuxT
            times, model_speeds, _ = wf.field_calculations.get_vsw(snap_id, run_parameters["base_name"], obs_time, vr, spinup_cme_days=run_parameters["spinup_time"], fcast_length=run_parameters["forecast_length"])
            omni_times, omni_speeds = wf.stats_functions.get_average_speeds([omni_dates], [omni_vs], spinup_time = run_parameters["spinup_time"], cadence = 24, target_times = times, verbose=run_parameters["verbose"])

            alltimes.append(times)
            allspeeds.append(model_speeds)
            allspeeds_ref.append(omni_speeds)

        persistence_time = 27.27 #Time in days for persistence model checks
        persistence_time_int = int(persistence_time*24)  #The cut in hours

        #Combine into one set of values, mainly consisting of NaNs (which is what we want really)
        times_avg, wsa = wf.stats_functions.get_average_speeds(alltimes, allspeeds, spinup_time=run_parameters["spinup_time"], cadence = 24, verbose=False)
        _, omni = wf.stats_functions.get_average_speeds(alltimes, allspeeds_ref, spinup_time=run_parameters["spinup_time"], cadence = 24, verbose=False)

        omni_shift = np.nan*omni
        omni_shift[persistence_time_int:] = omni[:-persistence_time_int]

        wsa_shift = np.nan*omni
        wsa_shift[persistence_time_int:] = wsa[:-persistence_time_int]

        cme_mask = wf.data_functions.get_cme_times(times_avg)
        invalid_times = np.where(cme_mask == 1)[0]

        wsa_filtered = wsa.copy()
        omni_filtered = omni.copy()
        wsa_filtered[invalid_times] = np.nan
        omni_filtered[invalid_times] = np.nan

        wsa_shift_filtered = wsa_shift.copy()
        omni_shift_filtered = omni_shift.copy()
        wsa_shift_filtered[invalid_times] = np.nan
        omni_shift_filtered[invalid_times] = np.nan

        #Can change this is necessary. But should be fine.
        best_metric = wsa_filtered

        nas = np.logical_or(np.isnan(wsa_filtered), np.isnan(omni_filtered))
        bestr, _ = pearsonr(best_metric[~nas], omni_filtered[~nas])
        best_rms = np.sqrt(np.nanmean((best_metric[~nas] - omni_filtered[~nas])**2))


        correlation_improvement = (1.0-bestr)
        rms_improvement = best_rms/100.0

        if correlation_improvement >= 1.0:
            #Introduce a harsh penalty for getting worse. Needs to be continuous though.
            correlation_improvement = correlation_improvement + 10*(correlation_improvement - 1.0)
        if rms_improvement >= 1.0:
            rms_improvement = rms_improvement + 10*(rms_improvement - 1.0)

        res = np.sqrt(correlation_improvement*rms_improvement)

        print('r, rms, res',  bestr, best_rms, res)
        result = res #CHANGE THIS TO RES IN DUE COURSE

        if np.isnan(result):
            result = 1e6
        return result


