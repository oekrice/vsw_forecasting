#This script requires test_predictions to have been run beforehand. Just loads in text files of speeds, which are hopefully aligned, and calculates/plots various things
#Seems to be best to do things separately like this to allow for a lot of flexibility.
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import wind_forecast as fcast
import cma

#Hamilton copy command: (to get the speeds and times)
# scp -r vgjn10@hamilton8.dur.ac.uk:/nobackup/vgjn10/projects/vsw_forecasting/data/raw_speeds ./data

for a in range(4):
    bs = [1,2,3,4]
    for b in bs:
        parameter_sources = ["raw", "wsa_nocmes", "rms_nocmes", "corr_nocmes"]
        scale_sources = ["raw", "ss", "rms", "dist", "ss_raw"]

        #Pick the desired combination here. The titles above should be kept consistent, but can obviously be added to if desired.
        #This code should be copied into 'scale_velocites' for consistency and neatness.
        parameter_select = a; scale_select = b

        parameter_source = parameter_sources[parameter_select]
        scale_source = scale_sources[scale_select]

        parameter_shortnames = ["raw", "dist", "rms", "corr"]
        parameter_shortname = parameter_shortnames[parameter_select]


        def scale_function(m, c, series):
            """
            Just does mx+c on the timeseries. Can optimise the skill scores based on that, hopefully
            """
            return series*m + c*100

        #Need to run through and do optimised scales and unoptimised scales. Automation is somewhat key as I can't keep track of 16 things.
        run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

        batch_names = []

        for i in range(8):
            batch_names.append(f'{parameter_source}_{i}')
        #fig = plt.figure(figsize=(12,6))

        if not os.path.exists('./data/scaling_data/'):
            os.mkdir('./data/scaling_data/')


        for i, batch_name in enumerate(batch_names):
            omni_fname = f'./data/raw_speeds/wsa_nocmes_0_wsa_speeds_ref.txt'

            if parameter_source == "raw":
                data_fname = f'./data/raw_speeds/wsa_nocmes_{i}_wsa_speeds.txt'
            else:
                data_fname = f'./data/raw_speeds/{batch_name}_wsa_scaled_speeds.txt'

            if not(os.path.exists(omni_fname) and os.path.exists(data_fname)):
                print('Files not found...', omni_fname, data_fname)
                continue

            wsa = np.loadtxt(data_fname, delimiter = ',')
            omni = np.loadtxt(omni_fname, delimiter = ',')

            if parameter_source == "raw":
                timeseries = np.loadtxt(f'./data/raw_speeds/wsa_nocmes_{i}_wsa_scaled_times.txt', dtype='datetime64[s]', delimiter = ',')
            else:
                timeseries = np.loadtxt(f'./data/raw_speeds/{batch_name}_wsa_scaled_times.txt', dtype='datetime64[s]', delimiter = ',')

            "Always filter for CMEs, but save this data separately"
            cme_mask = fcast.data_functions.get_cme_times(timeseries)
            invalid_times = np.where(cme_mask == 1)[0]

            wsa_filtered = wsa.copy()
            omni_filtered = omni.copy()
            wsa_filtered[invalid_times] = np.nan
            omni_filtered[invalid_times] = np.nan

            def test_scale_paras(scalers, series):
                series_scale = scale_function(scalers[0], scalers[1], series)
                if scale_source == "ss":
                    scores = fcast.stats_functions.do_met_stats(None, series_scale , omni_filtered, compare_to_persist=True)
                    score = 1.0 - np.mean(scores)
                elif scale_source == "ss_raw":
                    scores = fcast.stats_functions.do_met_stats(None, series_scale , omni_filtered, compare_to_persist=False)
                    score = 1.0 - np.mean(scores)
                elif scale_source == "raw":
                    raise Exception("Trying to optimise nothing. Don't do that")
                elif scale_source == "rms":
                    score = np.sqrt(np.nanmean((series_scale - omni_filtered)**2))
                elif scale_source == "dist":
                    score, _ = fcast.stats_functions.get_distribution_similarity(series_scale, omni_filtered, doplots=False)
                else:
                    raise Exception("Not done this bit yet")

                print(scalers, score)
                np.savetxt(f'./data/scaling_data/{parameter_shortname}_{scale_source}_{i}.txt', [score, scalers[0], scalers[1]])
                return score

            thresholds = np.arange(450,600,5)
            # scores = fcast.stats_functions.do_met_stats(None, wsa, omni_ref_0)
            score = test_scale_paras([1.0,0.0],wsa_filtered)
            print('Score', score)

            es = cma.CMAEvolutionStrategy([1.0,0.0], 0.10, {'verb_disp': 1, 'popsize': 10, 'maxiter': 50})
            es.optimize(test_scale_paras, args=([wsa]))
            es.result_pretty()

