#This script runs the base models (outflow/pfss and schatten) and is designed to be done in multiple threads, as it can be very slow.

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import cma
from datetime import datetime, timedelta

import wind_forecast as fcast  #This should now contain everything we need...

start = datetime(2010, 1, 1) #This CAN'T change for a given run name. BE CAREFUL

obs_times = [start + timedelta(days=i) for i in range(5478)]

net_names = ["pfss_2.5_gong", "pfss_5.0_gong", "outflow_2.5_gong", "outflow_5.0_gong", "pfss_2.5_hmi", "pfss_5.0_hmi", "outflow_2.5_hmi", "outflow_5.0_hmi"]

run_names = ["p2g", "p5g", "o2g", "o5g", "p2h", "p5h", "o2h", "o5h"]


if len(sys.argv) > 1:
    batch_id = int(sys.argv[1])
else:
    raise Exception('Specify batch.')

if len(sys.argv) > 2:
    nprocs = int(sys.argv[2])
    print(f'Running as a batch subset with {nprocs} threads')
else:
    nprocs = 1

if len(sys.argv) > 3:
    thread_number = int(sys.argv[3])
    print(f'This is thread number {thread_number}')
else:
    thread_number = 0

#Get the model setup depending on the batch numbers
if (batch_id//2)%2 == 0:
    is_pfss = True
else:
    is_pfss = False

if (batch_id%2) == 0:
    rss = 2.5
else:
    rss = 5.0

if (batch_id//4) == 0:
    source = "gong"
else:
    source = "hmi"

run_name = run_names[batch_id]

if nprocs > 1:
    #Need to establish the subsets for multiprocessing
    set_starts = np.linspace(0, len(obs_times), nprocs + 1).astype('int')
    snap_subset = np.arange(set_starts[thread_number], set_starts[thread_number+1])
else:
    snap_subset = None

fcast.calculate_base_model(run_name, obs_times, is_pfss=is_pfss, rss=rss, source=source, snap_subset=snap_subset)







