#This is a script designed to test a single parameter based on the entire test length.
#Each iteration will take a couple of hours, so it's important these are done in paralell.
#Could try testing all 9 parameters over a reasonable range, and iterating based on that? Is a bit bodgy, but I am rather grasping at straws here...

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import cma
import csv
import multiprocessing as mp
from datetime import datetime, timedelta

import wind_forecast as fcast  #This should now contain everything we need...
from dtaidistance import dtw
import random
import time

import matplotlib
#matplotlib.use('Agg')
#This script should just run the base model and HuxT, at a low resolution.
#Will automatically create a run ID with parameters encoded into the outputs, one hopes.

ntests = 10

if len(sys.argv) > 1:
    batch_select = int(sys.argv[1])
else:
    raise Exception('Specify batch number.')


if len(sys.argv) > 2:
    sweep_index = int(sys.argv[2])
else:
    raise Exception('Specify sweep index.')

parameter_to_change = sweep_index//ntests
parameter_pick = sweep_index%ntests

parameter_fname = f'./data/shared_data/sweep_paras_{batch_select}.csv'

if os.path.exists(parameter_fname):
    print('Alread-optimised parameters exist. Using them...')
    with open(parameter_fname, "r", encoding="utf-8") as f:
        data = csv.reader(f)
        for ri, row in enumerate(data):
            parameter_set = np.array(row).astype('float')
        #Just automate it. Otherwise we'll just get stuck...
        scale_limits = np.nan*np.ones((9,2))
        scale_limits[:,0] = parameter_set/1.5
        scale_limits[:,1] = parameter_set*1.5

else:
    parameter_set = [285,910,2/9,1.0,0.8,2  ,2,3,1]

    scale_limits = np.loadtxt('./data/shared_data/wsa_limits.dat', delimiter = ',')

select_scale = np.linspace(scale_limits[parameter_to_change][0], scale_limits[parameter_to_change][1], 10)
parameter_set[parameter_to_change] = select_scale[parameter_pick]
print('Changing parameter:', parameter_to_change)
print('To value:', select_scale[parameter_pick])

start = datetime(2010, 1, 1) #This CAN'T change for a given run name. BE CAREFUL
obs_times = [start + timedelta(days=i) for i in range(5478)]
#obs_times = [start + timedelta(days=i) for i in range(0, 10)]

plot_specific = -1   #Just evalulate a specific snap. Set to -1 for the latest one
plot_continuous = True   #Will wait for outputs and keep up (if possible)
find_min_sigma = True
use_neural_net = True

#Specify input parameters as a dictionary, which can be embiggened or ensmallened as necessary.
#Will check against whether sufficient data exists which matches what has been asked for, and will recalculate if necessary.
#Let's specify literally everything here, all the parameters which can happen.
#Will need a lookup table or equivalent to find data which matches things as they should.
#Can specify file name to look up WSA parameters? Yeah, probably.
run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

batch_id = batch_select

snap_subset = np.arange(len(obs_times)) #This should just start from the start now


#Get the model setup depending on the batch numbers
if (batch_id//2)%2 == 0:
    is_pfss = True
    model = "pfss"
else:
    is_pfss = False
    model = "outflow"
if (batch_id%2) == 0:
    rss = 2.5
else:
    rss = 5.0
if (batch_id//4) == 0:
    source = "gong"
else:
    source = "hmi"
run_name = run_names[batch_id]

nicetitle = f"{model}, rss = {rss}, source = {source}"

print('Using default wsa parameters as a reference for later')
batch_name = f"sweep_{batch_id}_{parameter_to_change}_{parameter_pick}"
velocity_type = "wsa_combined"

nicetitle = f"{model}, rss = {rss}, source = {source}"
test_parameters = {"observation_time": obs_times,
                "base_name": run_names[batch_id],
                "run_name": batch_name,
                "model_type": model,
                "calculate_base_model": False,
                "overwrite_base_model": False,
                "calculate_huxt": True,
                "r_ss": rss,
                "WSA_type": "standard",
                "WSA_parameters": None,
                "data_source": source,
                "resolutions": [120,180,360],
                "r_hb": 21.5,
                "match_flag": False,
                "velocity_type": velocity_type,
                "spinup_time": 5,
                "forecast_length": 5,
                "verbose": True,
                "optimisation_type": "least_squares",
                "do_plots": False}

theta = parameter_set

fcast.model_functions.run_model(test_parameters, theta=theta, snap_subset=snap_subset, save_speeds=True)

#Save out this specific parameter combination
with open(f'./data/raw_speeds/{batch_name}.csv', "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows([parameter_set])

