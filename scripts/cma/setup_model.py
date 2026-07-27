#MPI (ish) CMA test. As the standard sweep approach is terrible, but this might be a decent compromise. Can have quite big parameter sets.

import os, sys
import numpy as np
from datetime import datetime, timedelta

from model import run_model

start = datetime(2010, 1, 1) #This CAN'T change for a given run name. BE CAREFUL
obs_times = [start + timedelta(days=i) for i in range(5478)]
#obs_times = [start + timedelta(days=i) for i in range(60)]

test_single = False

nsamples = len(obs_times)
extend_current_run = False

theta_size = 9

#Specify input parameters as a dictionary, which can be embiggened or ensmallened as necessary.
#Will check against whether sufficient data exists which matches what has been asked for, and will recalculate if necessary.
#Let's specify literally everything here, all the parameters which can happen.
#Will need a lookup table or equivalent to find data which matches things as they should.
#Can specify file name to look up WSA parameters? Yeah, probably.

run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

if len(sys.argv) > 1:
    batch_id = int(sys.argv[1])
else:
    raise Exception('Specify batch number.')

if len(sys.argv) > 2:
    process = int(sys.argv[2])
else:
    raise Exception('Specify process ID.')

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

batch_name = f"cma_{batch_id}"
velocity_type = "wsa_scaled"   #To be used for the physics-informed parameter-changing. Just a select few of them. See if there are any mad patterns.


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
                "verbose": False,
                "optimisation_type": "correlation",
                "do_plots": False,
                "filter_cmes": True}

#Load in theta from the process ID.
input_data = np.loadtxt('./data/cma_data/parameters_%05d_%05d.txt' % (batch_id, process), delimiter = ',')
theta = input_data

print('Theta', theta)
snap_subset = np.arange(len(obs_times))
skillscore = run_model(test_parameters, theta=theta, snap_subset=snap_subset)

#Save out results. This is all. The model should wait for the rest.
np.savetxt('./data/cma_data/results_%05d_%05d.txt' % (batch_id, process), [skillscore])


