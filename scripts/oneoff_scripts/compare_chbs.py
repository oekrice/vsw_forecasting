#One-off script to test the different methods of plotting coronal hole boundaries (and to see if the data actually make any sense...)
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import wind_forecast as fcast

test_parameters = {"observation_time": datetime(2025, 3, 1),
                "base_name": "chb_test_base",
                "run_name": "chb_test",
                "model_type": "pfss",
                "calculate_base_model": True,
                "overwrite_base_model": True,
                "calculate_huxt": False,
                "r_ss": 2.5,
                "WSA_type": "standard",
                "WSA_parameters": None,
                "data_source": "gong",
                "resolutions": [120,180,360],
                #"resolutions": [60,90,180],
                "r_hb": 21.5,
                "match_flag": False,
                "velocity_type": "wsa",
                "spinup_time": 5,
                "forecast_length": 5,
                "verbose": True,
                "optimisation_type": "least_squares",
                "do_plots": False}

fcast.model_functions.run_model(test_parameters, use_old_chb_formula=True)
s0, ph0, br0, fs, chd = fcast.data_functions.load_chb_distances(test_parameters["base_name"] ,0)
old_chd = chd.flatten()

fcast.model_functions.run_model(test_parameters, use_old_chb_formula=False)
s0, ph0, br0, fs, chd = fcast.data_functions.load_chb_distances(test_parameters["base_name"] ,0)
new_chd = chd.flatten()

fig = plt.figure(figsize=(12,6))
plt.scatter(new_chd, old_chd,s=0.1)
plt.savefig('./plots/boundary_detection_comparison.png')
plt.show()
