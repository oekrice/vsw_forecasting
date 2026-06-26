'#These scripts are for the higher-level stuff, like calculating the entire HuxT prediction'

import os
import sys
from pathlib import Path

import numpy as np
import datetime
from scipy.io import netcdf_file

from datetime import datetime, timedelta
import csv

from .data_functions import check_existing_data
from .field_calculations import compute_vr, get_vsw
from .stats_functions import get_average_speeds, get_distribution_similarity, get_wasserstein_distance

def run_model(run_parameters, theta=np.zeros(8), snap_subset=None, iteration=0, output_distributions=False, save_speeds=False):
    """
    Using the parameter disctionary, will run the base model AND HuxT. If a run_name is provided, will save out data as it goes.
    So many variations need to be tested here, but I think I can do it...
    This will output a timeseries with speeds. That is all.
    """

    #Cycle through the reuqested observation times. Can be just one, or several
    if "observation_time" in run_parameters:
        if not isinstance(run_parameters["observation_time"], (list, np.ndarray)):
            run_parameters["observation_time"] = [run_parameters["observation_time"]]
    else:
        raise Exception('Observation time not provided.')

    if run_parameters["base_name"] is not None:
        if "match_flag" in run_parameters:
            data_lookup = check_existing_data(run_parameters, match_flag = run_parameters["match_flag"]) #Will eventually output flags as to whether these data have been already succesfully calculated with the given inputs.
        else:
            data_lookup = check_existing_data(run_parameters, match_flag = True)
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
            if data_lookup[snap_id] < 1: #The base model data doesn't exist, so recalculate
                fcast.calculate_outflow(snap_id,  obs_time, rss = run_parameters["r_ss"],
                                        overwrite=True, is_pfss=is_pfss, output_directory=run_name, save_snap=True,
                                        source= run_parameters["data_source"],resolutions=run_parameters["resolutions"])
            if data_lookup[snap_id] < 2: #The CHB and expansion factor data doesn't exist, so do that.
                fcast.calculate_chb_exp(snap_id,  run_name, r_hb = run_parameters["r_hb"], purge_data=True)

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
    omni_dates = omni_data[:,0]
    omni_vs = omni_data[:,1]

    if run_parameters["verbose"]:
        print('OMNI Data Obtained')

    skillscores = []

    alltimes = []
    allspeeds = []
    allspeeds_ref = []
    for si, snap_id in enumerate(snap_subset):
        obs_time = run_parameters["observation_time"][snap_id]
        if run_parameters["verbose"]:
            print(f'Running HuXT forecast model at time {obs_time}')

        if run_parameters["velocity_type"] == "wsa" or run_parameters["velocity_type"] == "wsa_scaled" :
            #This is the polynomial expression
            if theta is not None:
                if si == 0:
                    vr = compute_vr(snap_id, run_name, method="wsa_scaled", params = theta, doplot=run_parameters["do_plots"], iteration=iteration, huxt_name=run_parameters["run_name"])
                else:
                    vr = compute_vr(snap_id, run_name, method="wsa_scaled", params = theta, doplot=False, iteration=iteration, huxt_name=run_parameters["run_name"])
            else:
                if run_parameters["verbose"]:
                    print('Using default WSA parameters')
                if si == 0:
                    vr = compute_vr(snap_id, run_name, method="wsa", params = None, doplot=run_parameters["do_plots"], iteration=iteration, huxt_name=run_parameters["run_name"])
                else:
                    vr = compute_vr(snap_id, run_name, method="wsa", params = None, doplot=False, iteration=iteration, huxt_name=run_parameters["run_name"])

        elif run_parameters["velocity_type"] == "neural_net":
            if run_parameters["verbose"]:
                print('Using neural net parameters', theta)
            if si == 0:
                vr = compute_vr(snap_id, run_name, method="neural_net", params = theta, doplot=run_parameters["do_plots"], iteration=iteration, huxt_name=run_parameters["run_name"])
            else:
                vr = compute_vr(snap_id, run_name, method="neural_net", params = theta, doplot=False, iteration=iteration, huxt_name=run_parameters["run_name"])

        else:
            raise Exception("Velocity calculation type not recognised...")

        if run_parameters["verbose"]:
            print('Calculating VSW with HuxT...')

        #Using this velocity map (and nothing else?), run HuxT
        times, model_speeds, _ = get_vsw(snap_id, run_parameters["base_name"], obs_time, vr, spinup_cme_days=run_parameters["spinup_time"], fcast_length=run_parameters["forecast_length"])
        omni_times, omni_speeds = get_average_speeds([omni_dates], [omni_vs], spinup_time = 0, cadence = 24, target_times = times, verbose=run_parameters["verbose"])

        alltimes.append(times)
        allspeeds.append(model_speeds)
        allspeeds_ref.append(omni_speeds)

        if run_parameters["do_plots"] and si == 0:
            import matplotlib.pyplot as plt

            if not os.path.exists('plots'):
                os.mkdir('plots')
            if not os.path.exists(f'plots/{run_parameters["run_name"]}'):
                os.mkdir(f'plots/{run_parameters["run_name"]}')

            fig = plt.figure(figsize=(10,7))
            plt.plot(omni_times, model_speeds)
            plt.plot(omni_times, omni_speeds)
            #plt.title(f'Distance metric: {skillscores[-1]}')
            plt.savefig('./plots/%s/timeseries_%05d.png' % (run_parameters["run_name"], iteration))
            print(f'Plot saved to {'./plots/%s/timeseries_%05d.png' % (run_parameters["run_name"], iteration)}')
            #plt.show()
            plt.close()

    if run_parameters["optimisation_type"] == "dtw":
        from dtaidistance import dtw

        #For DTW, need to compare each set individually, which I concede is a bit of a pain.
        for si in range(len(allspeeds)):
            omni_speeds = allspeeds_ref[si]
            model_speeds = allspeeds[si]
            dtw_distance = dtw.distance(omni_speeds, model_speeds)
            skillscores.append(dtw_distance)
    elif run_parameters["optimisation_type"] == "least_squares":

        times_avg, speeds_avg = get_average_speeds(alltimes, allspeeds, spinup_time = 0, cadence = 24, verbose=True)
        times_ref_avg, speeds_ref_avg = get_average_speeds(alltimes, allspeeds_ref, spinup_time = 0, cadence = 24, verbose=True)

        if save_speeds:
            if run_parameters["verbose"] == True:
                print('Saving out raw speed data...')

            times_avg = np.array(times_avg, dtype='datetime64[s]')
            times_ref_avg = np.array(times_ref_avg, dtype='datetime64[s]')
            #Save out the speeds to a normal txt file, so analysis on them is easy. Do need all the information though.
            if not os.path.exists('./data/raw_speeds/'):
                os.mkdir('./data/raw_speeds/')
            np.savetxt(f'./data/raw_speeds/{run_parameters["run_name"]}_{run_parameters["velocity_type"]}_speeds.txt', speeds_avg, delimiter = ',')
            np.savetxt(f'./data/raw_speeds/{run_parameters["run_name"]}_{run_parameters["velocity_type"]}_speeds_ref.txt', speeds_ref_avg, delimiter = ',')
            np.savetxt(f'./data/raw_speeds/{run_parameters["run_name"]}_{run_parameters["velocity_type"]}_times.txt', times_avg, fmt='%s', delimiter = ',')
            np.savetxt(f'./data/raw_speeds/{run_parameters["run_name"]}_{run_parameters["velocity_type"]}_times_ref.txt', times_ref_avg, fmt='%s', delimiter = ',')

            #print('Overall maximum speed (and ref):', np.max(allspeeds), np.max(allspeeds_ref))

            if run_parameters["verbose"] == True:
                print(f'Raw speed data saved with root {run_parameters["run_name"]}_{run_parameters["velocity_type"]}')
        speeds = np.concatenate(allspeeds)
        speeds_ref = np.concatenate(allspeeds_ref)

        leastsquares_distance = np.sqrt(np.mean((speeds - speeds_ref)**2))
        skillscores.append(leastsquares_distance)
    elif run_parameters["optimisation_type"] == "wasserstein":
        speeds = np.concatenate(allspeeds)
        speeds_ref = np.concatenate(allspeeds_ref)
        wasserstein_distance = get_wasserstein_distance(speeds, speeds_ref, huxt_name=run_parameters["run_name"], iteration=iteration, doplots=run_parameters["do_plots"])
        skillscores.append(wasserstein_distance)
    elif run_parameters["optimisation_type"] == "distribution":
        speeds = np.concatenate(allspeeds)
        speeds_ref = np.concatenate(allspeeds_ref)
        distribution_similarity, dists = get_distribution_similarity(speeds, speeds_ref, huxt_name=run_parameters["run_name"], iteration=iteration, doplots=run_parameters["do_plots"])
        skillscores.append(distribution_similarity)
        print('Overall maximum speed (and ref):', np.max(speeds), np.max(speeds_ref))
    else:
        raise Exception('Optimisation type not recognised')

    if run_parameters["verbose"]:
        print('Skillscore', skillscores[0])

    if not output_distributions:
        return np.mean(skillscores)
    else:
        return np.mean(skillscores), dists

def do_model_statistics(run_name, compare_to_persist=True):
    import matplotlib.pyplot as plt


    mean_times = np.loadtxt('./data/' + run_name + '/times.txt', dtype='datetime64[s]', delimiter = ',')
    mean_speeds = np.loadtxt('./data/' + run_name + '/vs.txt', delimiter = ',')

    times_omni, speeds_omni = fcast.obtain_omni_data(run_name, mean_times, overwrite=False)

    fig = plt.figure(figsize = (10,7))

    plt.plot(mean_times, mean_speeds, label = run_name)
    plt.plot(times_omni, speeds_omni, label = 'OMNI Data', c = 'black')

    plt.legend()
    plt.savefig(f'./plots/speeds_{run_name}.png')
    plt.close()

    thresholds = np.arange(450,600,5)
    skillscores = fcast.do_met_stats(times_omni, model_speeds=mean_speeds, reference_speeds=speeds_omni, compare_to_persist=compare_to_persist, persistence_cadence=int(24*27.7), thresholds=thresholds)

    fig = plt.figure(figsize = (10,5))

    plt.plot(thresholds, skillscores, label = run_name)
    plt.title(f'Mean skillscore = {np.mean(skillscores)}')
    plt.legend()
    if compare_to_persist:
        plt.savefig(f'./plots/skillscores_relative_{run_name}.png')
    else:
        plt.savefig(f'./plots/skillscores_raw_{run_name}.png')

    plt.close()

    return skillscores


def calculate_base_model(run_name, obs_times, overwrite=False, is_pfss=False, rss=2.5, source="hmi", snap_subset=None):

    if snap_subset is None:
        snap_subset = np.arange(len(obs_times))

    alltimes = []; allspeeds = []
    for snap_id in snap_subset:
        obs_time = obs_times[snap_id]
        #This is a pretty nice way of doing it I think
        print(f'Running snap {snap_id} from batch "{run_name}""')
        print(f'Is PFSS: {is_pfss}, rss = {rss}, Data source:', source)
        data_fname = f'./data/{run_name}/chb_{snap_id:09d}.nc'
        if not os.path.exists(data_fname) or overwrite:
            #Calculate the outflow fields
            fcast.calculate_outflow(snap_id,  obs_time, rss = rss, overwrite=False, is_pfss=is_pfss, output_directory=run_name, save_snap=True, source=source)

            #Calculate the schatten extension, expansion factors and coronal hole distances. I think putting this as a separate function is just fine. Doesn't require any information beyond the file names.
            fcast.calculate_chb_exp(snap_id,  run_name, purge_data=True)

    return None
