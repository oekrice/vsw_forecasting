#This is a testbed for making the forecast scripts actually nice, and doing it all properly and things. HA, that went well!

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import cma
import csv
from datetime import datetime

import wind_forecast as fcast  #This should now contain everything we need...

#This script should just run the base model and HuxT, at a low resolution.
#Will automatically create a run ID with parameters encoded into the outputs, one hopes.

obs_time = [datetime(2010, 1, 1)]#, datetime(2010, 1, 2), datetime(2010, 1, 3)]

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
                  "r_ss": 5.0,
                  "WSA_type": "standard",
                  "WSA_parameters": None,
                  "verbose": True,
                  "data_source": "hmi",
                  "resolutions": [60,90,180],
                  "r_hb": 21.5,
                  "match_flag": True}

def check_existing_data(run_parameters, match_flag=True):
    """
    If a run name is given, check against the lookup table as to whether the lower boundary data exists for these parameters.
    Will raise appropriate errors, if appropriate.
    If not, will provide a database of what has been calculated for this 'name'

    Code:
    -1 for file exists but metadata is different. Flag if necessary.
    0 for file not existing
    1 for base model done
    2 for chbmap done
    """
    fname = run_parameters["base_name"]
    directory_fname = f'./data/{fname}/directory.csv'
    if os.path.exists(directory_fname):
        directory_data = []
        with open(directory_fname, "r", encoding="utf-8") as f:
            data = csv.reader(f)
            for row in data:
                directory_data.append(row)
    else:
        directory_data = []

    #Run through each of the desired snaps and see what exists. Need integer codes for this really. They are now defined above
    check_codes = [0]*len(run_parameters["observation_time"])
    for row in directory_data:
        if row[0] == "base":
            snap_id = int(row[1])
            expected_row = ["base", str(snap_id), str(run_parameters["observation_time"][snap_id]), str(run_parameters["r_ss"]),
                            run_parameters["model_type"], run_parameters["data_source"],
                            str(run_parameters["resolutions"][0]), str(run_parameters["resolutions"][1]), str(run_parameters["resolutions"][2])]
            if row == expected_row and check_codes[snap_id] < 1:
                check_codes[snap_id] = 1
            else:
                check_codes[snap_id] = -1
                if match_flag:
                    print("Existing parameter set:", row)
                    print("Specified parameter set", expected_row)
                    raise Exception("Existing data doesn't match these parameters in this run. Aborting... To overwrite with these new parameters put 'match_flag=False''")
                else:
                    print("New parameters do not match existing ones, but proceeding anyway. Set 'match_flag=True' to stop this")
        elif row[0] == "chbetc":
            snap_id = int(row[1])
            expected_row = ["chbetc", str(snap_id), str(run_parameters["r_hb"])]
            if row == expected_row:
                check_codes[snap_id] = 2
            else:
                check_codes[snap_id] = -1
                if match_flag:
                    print("Existing parameter set:", row)
                    print("Specified parameter set", expected_row)
                    raise Exception("Existing data doesn't match these parameters in this run. Aborting... To overwrite with these new parameters put 'match_flag=False''")
                else:
                    print("New parameters do not match existing ones, but proceeding anyway. Set 'match_flag=True' to stop this")
        else:
            raise Exception('Not done this bit yet')

    return check_codes

def run_model(run_parameters):
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

    for snap_id, obs_time in enumerate(run_parameters["observation_time"]):
        if run_parameters["calculate_base_model"]:  #Do the PFSS/Outflow calculation

            #Check if base model already exists for these data, and the metadata all match (should put this check in an extra function. Each 'data' file should have a lookup table for it, I think)
            if data_lookup[snap_id] < 1: #The base model data doesn't exist, so recalculate
                fcast.calculate_outflow(snap_id,  obs_time, rss = run_parameters["r_ss"],
                                        overwrite=True, is_pfss=is_pfss, output_directory=run_name, save_snap=True,
                                        source= run_parameters["data_source"],resolutions=run_parameters["resolutions"])
            if data_lookup[snap_id] < 2: #The CHB and expansion factor data doesn't exist, so do that.
                fcast.calculate_chb_exp(snap_id,  run_name, r_hb = run_parameters["r_hb"], purge_data=True)

    sys.exit()

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

run_model(test_parameters)
