#This is a testbed for making the forecast scripts actually nice, and doing it all properly and things. HA, that went well!

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import cma
from datetime import datetime

import wind_forecast as fcast  #This should now contain everything we need...

#This script should just run the base model and HuxT, at a low resolution.
#Will automatically create a run ID with parameters encoded into the outputs, one hopes.
print('running')

obs_time = datetime(2010, 1, 1)

#Specify input parameters as a dictionary, which can be embiggened or ensmallened as necessary.
#Will check against whether sufficient data exists which matches what has been asked for, and will recalculate if necessary.
#Let's specify literally everything here, all the parameters which can happen.
#Will need a lookup table or equivalent to find data which matches things as they should.
#Can specify file name to look up WSA parameters? Yeah, probably.

test_parameters = {"observation_time": obs_time,
                  "run_name": None,
                  "model_type": "outflow",
                  "calculate_base_model": True,
                  "calculate_huxt": True,
                  "r_ss": 2.5,
                  "WSA_type": "standard",
                  "WSA_parameters": None,
                  "verbose": True}


def run_model(run_parameters):
    """
    Using the parameter disctionary, will run the base model AND HuxT. If a run_name is provided, will save out data as it goes.
    So many variations need to be tested here, but I think I can do it...
    This will output a timeseries with speeds. That is all.
    """

    #Cycle through the reuqested observation times. Can be just one, or several


    sys.exit()
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

#run_model(test_parameters)
