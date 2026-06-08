'#These scripts are for the higher-level stuff, like calculating the entire HuxT prediction'

import os
import sys
from pathlib import Path

import numpy as np
import datetime
from scipy.io import netcdf_file

this_directory = os.getcwd() + "/"
sys.path.append(this_directory+"viz/tools")
import wind
import wind_forecast as fcast

import matplotlib.pyplot as plt

#from field_calculations import calculate_outflow, compute_coronal_hole_distances, compute_vr, get_vsw, get_vsw_dynamic, get_average_speeds

def do_model_run(run_name, obs_times, velocity_parameters=[285, 625+285, 0.22222, 1, 0.8, 2, 2, 3], overwrite=False, is_pfss=False, rss=2.5, source="hmi", snap_subset=None, Net=None, doplot=False, verbose=False):

    if snap_subset is None:
        snap_subset = np.arange(len(obs_times))

    alltimes = []; allspeeds = []
    for snap_id in snap_subset:
        if verbose:
            print(f'Doing snap {snap_id} of {len(snap_subset)}')
        obs_time = obs_times[snap_id]
        #This is a pretty nice way of doing it I think
        # print(f'Running snap {snap_id} from batch "{run_name}"')
        # print(f'Is PFSS: {is_pfss}, rss = {rss}, Parameters: {velocity_parameters}', 'Data source:', source)
        data_fname = f'./data/{run_name}/chb_{snap_id:09d}.nc'
        if not os.path.exists(data_fname) or overwrite:
            #Calculate the outflow fields
            fcast.calculate_outflow(snap_id,  obs_time, rss = rss, overwrite=False, is_pfss=is_pfss, output_directory=run_name, save_snap=True, source=source)

            #Calculate the schatten extension, expansion factors and coronal hole distances. I think putting this as a separate function is just fine. Doesn't require any information beyond the file names.
            fcast.calculate_chb_exp(snap_id,  run_name, purge_data=True)

        #Compute velocity map (this requires use of variable parameters)
        if Net is None:
            vr = fcast.compute_vr(snap_id,  run_name, params = velocity_parameters)
        else:
            vr = fcast.compute_vr_net(snap_id,  run_name, Net)

        #Using this velocity map (and nothing else?), run HuxT
        times, speeds, _ = fcast.get_vsw(snap_id, run_name, obs_time, vr, spinup_cme_days=5, fcast_length=5)

        alltimes.append(times)
        allspeeds.append(speeds)

        if snap_id%100 == 0 and snap_id != 0:  #Do helpful things as it goes along as Hamilton doesn't tell any progress otherwise

            mean_times, mean_speeds = fcast.get_average_speeds(alltimes, allspeeds, spinup_time = 5, cadence = 24)

            np.savetxt('./data/' + run_name + '/times.txt', mean_times.astype("datetime64[s]").astype(str), fmt="%s")
            np.savetxt('./data/' + run_name + '/vs.txt', mean_speeds, delimiter = ',')


    mean_times, mean_speeds = fcast.get_average_speeds(alltimes, allspeeds, spinup_time = 5, cadence = 24)

    np.savetxt('./data/' + run_name + '/times.txt', mean_times.astype("datetime64[s]").astype(str), fmt="%s")
    np.savetxt('./data/' + run_name + '/vs.txt', mean_speeds, delimiter = ',')

    times_omni, speeds_omni = fcast.obtain_omni_data(run_name, mean_times, overwrite=False)

    thresholds = np.arange(400,600,5)
    skillscores = fcast.do_met_stats(times_omni, model_speeds=mean_speeds, reference_speeds=speeds_omni, compare_to_persist=compare_to_persist, persistence_cadence=int(24*27.7), thresholds=thresholds)

    if doplot:
        fig = plt.figure(figsize = (10,7))

        plt.plot(mean_times, mean_speeds, label = run_name)
        plt.plot(times_omni, speeds_omni, label = 'OMNI Data', c = 'black')

        plt.title(f'Avg. skillscore: {np.mean(skillscores):.2f}')
        plt.legend()
        plt.savefig(f'./plots/speeds_{run_name}.png')
        plt.close()


    if Net is not None:
        distance_metric = np.nanmean((mean_speeds - speeds_omni)**2)
        return skillscores, distance_metric
    else:
        return skillscores

def do_model_statistics(run_name, compare_to_persist=True):


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
