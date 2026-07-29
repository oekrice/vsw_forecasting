#This is a very one-off script to take 'optimum' parameters from a set of Hamilton runs and save them out as raw parameter values.
#Probably lots of scope for this to go wrong but meh.

import numpy as np
import os, sys
import csv

do_combine = True
base_name = 'crot_fit'
run_title = 'cma_data'

if not do_combine:
    destination_fname = './data/shared_data/optimum_parameters.csv'
    log_fname = f"./data/{run_title}/log.csv"
else:
    destination_fname = './data/shared_data/optimum_parameters_combine.csv'
    log_fname = f"./data/{run_title}/log.csv"

if os.path.exists(destination_fname):
    os.remove(destination_fname)

for run_num in range(8):

    if not do_combine:
        #Copy log files from Hamilton
        if not os.path.exists(f"./data/{run_title}/{run_num}_log.csv"):
            copy_command = f'scp -r vgjn10@hamilton8.dur.ac.uk:/nobackup/vgjn10/projects/vsw_forecasting/data/{run_title}/{run_num}_log.csv ./data/{run_title}/'
            if not os.path.exists(f"./data/{run_title}/"):
                os.mkdir(f"./data/{run_title}/")
            os.system(copy_command)

        bestscore = 1e6
        bestrow = None
        #Find base parameters giving the best scores
        with open(f"./data/{run_title}/{run_num}_log.csv", "r", encoding="utf-8") as f:
            log_data = csv.reader(f)
            for row in log_data:
                if not row[0].isnumeric():
                    continue
                score = float(row[1])
                if score < bestscore:
                    bestscore = score
                    bestrow = row
        bestrow = np.array(bestrow, dtype='float')

    else:
        #Copy log files from Hamilton
        if not os.path.exists(f"./data/{run_title}/{run_num}_combine_log.csv"):
            copy_command = f'scp -r vgjn10@hamilton8.dur.ac.uk:/nobackup/vgjn10/projects/vsw_forecasting/data/{run_title}/{run_num}_combine_log.csv ./data/{run_title}/'
            if not os.path.exists(f"./data/{run_title}/"):
                os.mkdir(f"./data/{run_title}/")
            os.system(copy_command)

        bestscore = 1e6
        bestrow = None
        #Find base parameters giving the best scores
        with open(f"./data/{run_title}/{run_num}_combine_log.csv", "r", encoding="utf-8") as f:
            log_data = csv.reader(f)
            for row in log_data:
                if not row[0].isnumeric():
                    continue
                score = float(row[1])
                if score < bestscore:
                    bestscore = score
                    bestrow = row
        bestrow = np.array(bestrow, dtype='float')
    #Convert these parameters into 'raw parameter space'.
    #Requires the limits to be consistet throughout, but this can be stolen from compute_vr, I think.
    parameter_set = np.zeros(9)
    params = bestrow[2:]
    scale_limits= np.loadtxt('./data/shared_data/wsa_limits.dat', delimiter = ',')

    def scale_parameter(i, x):
        return 0.5*(1.0 + np.tanh(x))*(scale_limits[i][1] - scale_limits[i][0]) + scale_limits[i][0]

    parameter_set[0] = scale_parameter(0, params[0])
    parameter_set[1] = scale_parameter(1, params[1])
    parameter_set[2] = scale_parameter(2, params[2])
    parameter_set[3] = scale_parameter(3, params[3])
    parameter_set[4] = scale_parameter(4, params[4])
    parameter_set[5] = scale_parameter(5, params[5])
    parameter_set[6] = scale_parameter(6, params[6])
    parameter_set[7] = scale_parameter(7, params[7])
    parameter_set[8] = scale_parameter(8, params[8])

    #Save these out in a sensible place, to be read in by the 'speeds' script
    print(parameter_set)

    if os.path.exists(destination_fname):
        #This directory already exists. Hopefully with proper header information etc
        directory_data = []
        with open(destination_fname, "r", encoding="utf-8") as f:
            data = csv.reader(f)
            for row in data:
                directory_data.append(row)
    else:
        directory_data = []

    new_row_data = []
    for i in range(len(parameter_set)):
        new_row_data.append(parameter_set[i])
    directory_data.append(new_row_data)

    with open(destination_fname, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(directory_data)



