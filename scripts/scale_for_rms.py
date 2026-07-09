#This script requires test_predictions to have been run beforehand. Just loads in text files of speeds, which are hopefully aligned, and calculates/plots various things
#Seems to be best to do things separately like this to allow for a lot of flexibility.
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import wind_forecast as fcast
import cma

filter_for_cmes = True

cmap = plt.get_cmap("tab10")


def make_nicetitle(id):
    if (id//2)%2 == 0:
        is_pfss = True
        model = "pfss"
    else:
        is_pfss = False
        model = "outflow"
    if (id%2) == 0:
        rss = 2.5
    else:
        rss = 5.0
    if (id//4) == 0:
        source = "gong"
    else:
        source = "hmi"

    nicetitle = f"{model}, rss = {rss}, source = {source}"
    return nicetitle

def scale_function(m, c, series):
    """
    Just does mx+c on the timeseries. Can optimise the skill scores based on that, hopefully
    """
    return series*m + c*100

#Need to run through and do optimised scales and unoptimised scales. Automation is somewhat key as I can't keep track of 16 things.
run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

unoptimised_scales = [None, None, None, None, None, None, [1.08, -0.156]]
optimised_scales = [None, None, None, None, None, None, None]

doruns = np.arange(0,8)
batch_names = []
for i in range(8):
    batch_names.append(f'wsa_nocmes_{i}')
#fig = plt.figure(figsize=(12,6))

if not os.path.exists('./data/scaling_data/'):
    os.mkdir('./data/scaling_data/')

for i in doruns:

    for is_optimised in [False, True]:  #Do the optimised run and the unoptimised run

        if is_optimised:
            opt_title = 'dist'
        else:
            opt_title = 'raw'

        batch_name = batch_names[i]
        wsa_fname = f'./data/raw_speeds/{batch_name}_wsa_speeds.txt'
        optimised_fname = f'./data/raw_speeds/{batch_name}_wsa_scaled_speeds.txt'
        ref_fname_0 = f'./data/raw_speeds/{batch_name}_wsa_speeds_ref.txt'
        ref_fname_1 = f'./data/raw_speeds/{batch_name}_wsa_scaled_speeds_ref.txt'

        timeseries = np.loadtxt(f'./data/raw_speeds/{batch_name}_wsa_scaled_times.txt', dtype='datetime64[s]', delimiter = ',')

        if not(os.path.exists(wsa_fname) and os.path.exists(optimised_fname) and os.path.exists(ref_fname_0) and os.path.exists(ref_fname_1)):
            continue

        wsa = np.loadtxt(wsa_fname, delimiter = ',')
        optimised_wsa = np.loadtxt(optimised_fname, delimiter = ',')
        omni_ref_0 = np.loadtxt(ref_fname_0, delimiter = ',')
        omni_ref_1 = np.loadtxt(ref_fname_1, delimiter = ',')

        if filter_for_cmes:
            cme_mask = fcast.data_functions.get_cme_times(timeseries, also_filter_persistence=False)
            invalid_times = np.where(cme_mask == 1)[0]
            wsa[invalid_times] = np.nan
            optimised_wsa[invalid_times] = np.nan
            omni_ref_0[invalid_times] = np.nan
            omni_ref_1[invalid_times] = np.nan

        def test_scale_paras(scalers, series):

            series_scale = scale_function(scalers[0], scalers[1], series)

            rms = np.sqrt(np.nanmean((series_scale - omni_ref_0)**2))

            np.savetxt(f'./data/scaling_data/{batch_name}_{opt_title}_rms.txt', [rms, scalers[0], scalers[1]])
            return rms

        thresholds = np.arange(450,600,5)
        # scores = fcast.stats_functions.do_met_stats(None, wsa, omni_ref_0)
        score = test_scale_paras([1.0,0.0],wsa)
        print('Score', score)

        es = cma.CMAEvolutionStrategy([1.0,0.0], 0.10, {'verb_disp': 1, 'popsize': 10, 'maxiter': 50})
        if is_optimised:
            es.optimize(test_scale_paras, args=([optimised_wsa]))
        else:
            es.optimize(test_scale_paras, args=([wsa]))
        es.result_pretty()

    # plt.plot(thresholds, scores, label = f'{make_nicetitle(i)}, default parameters')

    #scores = fcast.stats_functions.do_met_stats(None, optimised_wsa, omni_ref_0)
    # plt.plot(thresholds, scores, label = f'{make_nicetitle(i)}, optimised parameters')

# plt.legend()
# plt.title('Raw persistence skill scores')
#
# plt.tight_layout()
# plt.show()
