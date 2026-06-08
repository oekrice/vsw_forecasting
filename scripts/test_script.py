#This is a testbed for making the forecast scripts actually nice, and doing it all properly and things. HA, that went well!

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import cma
import csv
from datetime import datetime, timedelta

import wind_forecast as fcast  #This should now contain everything we need...

#This script should just run the base model and HuxT, at a low resolution.
#Will automatically create a run ID with parameters encoded into the outputs, one hopes.

obs_time = [datetime(2010, 1, 1), datetime(2010, 1, 20)]#, datetime(2010, 1, 2), datetime(2010, 1, 3)]

#Specify input parameters as a dictionary, which can be embiggened or ensmallened as necessary.
#Will check against whether sufficient data exists which matches what has been asked for, and will recalculate if necessary.
#Let's specify literally everything here, all the parameters which can happen.
#Will need a lookup table or equivalent to find data which matches things as they should.
#Can specify file name to look up WSA parameters? Yeah, probably.

test_parameters = {"observation_time": obs_time,
                  "base_name": "test1",
                  "model_type": "pfss",
                  "calculate_base_model": True,
                  "overwrite_base_model": False,
                  "calculate_huxt": True,
                  "r_ss": 2.5,
                  "WSA_type": "standard",
                  "WSA_parameters": None,
                  "verbose": True,
                  "data_source": "hmi",
                  "resolutions": [60,90,180],
                  "r_hb": 21.5,
                  "match_flag": False,
                  "velocity_type": "wsa",
                  "spinup_time": 5,
                  "forecast_length": 5}


def run_model(run_parameters, snap_subset=None):
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
            data_lookup = fcast.check_existing_data(run_parameters, match_flag = run_parameters["match_flag"]) #Will eventually output flags as to whether these data have been already succesfully calculated with the given inputs.
        else:
            data_lookup = fcast.check_existing_data(run_parameters, match_flag = True)
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

            #Check if base model already exists for these data, and the metadata all match (should put this check in an extra function. Each 'data' file should have a lookup table for it, I think)
            if data_lookup[snap_id] < 1: #The base model data doesn't exist, so recalculate
                fcast.calculate_outflow(snap_id,  obs_time, rss = run_parameters["r_ss"],
                                        overwrite=True, is_pfss=is_pfss, output_directory=run_name, save_snap=True,
                                        source= run_parameters["data_source"],resolutions=run_parameters["resolutions"])
            if data_lookup[snap_id] < 2: #The CHB and expansion factor data doesn't exist, so do that.
                fcast.calculate_chb_exp(snap_id,  run_name, r_hb = run_parameters["r_hb"], purge_data=True)

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
    print('OMNI Data Obtained')

    if False:
        alltimes = []
        allspeeds = []
        for snap_id in snap_subset:
            obs_time = run_parameters["observation_time"][snap_id]

            if run_parameters["velocity_type"] == "wsa":
                #This is the polynomial expression
                vr = fcast.compute_vr(snap_id, run_name, method="wsa_scaled", params = None, doplot=False)
            else:
                raise Exception("Velocity calculation type not recognised...")

            #Using this velocity map (and nothing else?), run HuxT
            times, speeds, _ = fcast.get_vsw(snap_id, run_name, obs_time, vr, spinup_cme_days=run_parameters["spinup_time"], fcast_length=run_parameters["forecast_length"])
            alltimes.append(times)
            allspeeds.append(speeds)


        mean_times, mean_speeds = fcast.get_average_speeds(alltimes, allspeeds, spinup_time = 5, cadence = 24)

        np.savetxt('./data/' + run_name + '/times.txt', mean_times.astype("datetime64[s]").astype(str), fmt="%s")
        np.savetxt('./data/' + run_name + '/vs.txt', mean_speeds, delimiter = ',')

    mean_times = np.loadtxt('./data/' + run_name + '/times.txt', dtype='datetime64[s]', delimiter = ',')
    mean_speeds = np.loadtxt('./data/' + run_name + '/vs.txt', delimiter = ',')

    mean_omni_times, mean_omni_speeds = fcast.get_average_speeds([omni_dates], [omni_vs], spinup_time = 0, cadence = 24, target_times = mean_times)
    plt.plot(mean_omni_times, mean_omni_speeds)
    plt.plot(mean_times, mean_speeds)
    plt.show()

    sys.exit()

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

run_model(test_parameters)
