#This is a testbed for making the forecast scripts actually nice, and doing it all properly and things. HA, that went well!

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import cma
import csv
from datetime import datetime, timedelta

import wind_forecast as fcast  #This should now contain everything we need...
from dtaidistance import dtw


#This script should just run the base model and HuxT, at a low resolution.
#Will automatically create a run ID with parameters encoded into the outputs, one hopes.

obs_time = [datetime(2010, 1, 1), datetime(2010, 1, 2), datetime(2010, 1, 3)]

#Specify input parameters as a dictionary, which can be embiggened or ensmallened as necessary.
#Will check against whether sufficient data exists which matches what has been asked for, and will recalculate if necessary.
#Let's specify literally everything here, all the parameters which can happen.
#Will need a lookup table or equivalent to find data which matches things as they should.
#Can specify file name to look up WSA parameters? Yeah, probably.

test_parameters = {"observation_time": obs_time,
                  "base_name": "test2",
                  "run_name": "test_run",
                  "model_type": "outflow",
                  "calculate_base_model": True,
                  "overwrite_base_model": False,
                  "calculate_huxt": True,
                  "r_ss": 2.5,
                  "WSA_type": "standard",
                  "WSA_parameters": None,
                  "verbose": True,
                  "data_source": "hmi",
                  "resolutions": [120,180,360],
                  "r_hb": 21.5,
                  "match_flag": False,
                  "velocity_type": "wsa",
                  "spinup_time": 5,
                  "forecast_length": 5,
                  "verbose": True,
                  "optimisation_type": "wasserstein",
                  "do_plots": True}

fcast.run_model(test_parameters)
