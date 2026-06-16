#This script is just for copying the directory data over to Hamilton, so it doesn't need to be calculated again.
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import cma
import csv
from datetime import datetime, timedelta

import wind_forecast as fcast  #This should now contain everything we need...
from dtaidistance import dtw

start = datetime(2010, 1, 1) #This CAN'T change for a given run name. BE CAREFUL

obs_times = [start + timedelta(days=i) for i in range(5478)]

run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]

for batch_id in range(0,1):

    source_location = f'./data/{run_names[batch_id]}/directory.csv'
    destination_location = f'/nobackup/vgjn10/projects/vsw_forecasting/data/{run_names[batch_id]}'

    command = 'scp ' + source_location + ' vgjn10@hamilton8.dur.ac.uk:' + destination_location
    print(source_location)
    print(destination_location)
    print(command)

    os.system(command)
