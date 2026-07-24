#Overall wrapper for the CMA using many processes.
#Can update this to use history is possible/necessary?

#Load in an existing CMA log if it exists?
import os, sys
import cma
import csv
import numpy as np
import subprocess
import time

do_combine = True

if len(sys.argv) > 1:
    batch_id = int(sys.argv[1])
else:
    raise Exception('Specify batch number.')

if not os.path.exists('./data/cma_data'):
    os.mkdir('./data/cma_data')

overall_log = []

if os.path.exists('./data/cma_data/overall_log.csv'):
    with open('./data/cma_data/overall_log.csv', "r", encoding="utf-8") as f:
        data = csv.reader(f)
        for row in data:
            overall_log.append(row)
else:
    overall_log = []

print(overall_log)

def await_solutions(batch_id, popsize):
    #Runs the processes (this can be tweaked for Hamilton?) and await the solutions
    #Wait (for a certain amount of time...) for the optimum parameters to be saved out. Bit wasteful just to have one file, but meh.

    if "SLURM_JOB_ID" in os.environ:# or True:
        #raise Exception('This is running on slurm. Not done that bit yet...)')
        subprocess.Popen(["sbatch",  "scripts/cma/cma_array.sh", str(batch_id)])

    else:
        print('Running locally (not on slurm).')
        for i in range(popsize):
            subprocess.Popen(["python",  "scripts/cma/setup_model.py", str(batch_id), str(i)])

    #Wait for result files to be made.
    results = np.full(popsize, fill_value=-1.0)
    time_start = time.time()
    timeout = 3600*6 # Three hours for timeout, eventually?
    while np.nanmin(results) < 0 and time.time() - time_start < timeout:
        for i in range(popsize):
            if results[i] > 0 or np.isnan(results[i]):
                continue
            if os.path.exists('./data/cma_data/results_%05d_%05d.txt' % (batch_id, i)):
                time.sleep(1.0) # Let it write properly, just in case
                results[i] = np.loadtxt('./data/cma_data/results_%05d_%05d.txt' % (batch_id, i))
                print('Updated!', results)

        print(time.time() - time_start, results, np.nanmin(results) )
        time.sleep(1.0)

    results[results < 0.0] = 1e6
    return results


data_length = 0
if do combine:
    directory_fname = f'./data/cma_data/{batch_id}_log.csv'
else:
    directory_fname = f'./data/cma_data/{batch_id}_combine_log.csv'

if os.path.exists(directory_fname):
    #This directory already exists. Hopefully with proper header information etc
    directory_data = []
    with open(directory_fname, "r", encoding="utf-8") as f:
        data = csv.reader(f)
        for row in data:
            directory_data.append(row)
    directory_data = np.array(directory_data, dtype='float')
    all_parameters = directory_data[:,2:]
    all_solutions = directory_data[:,1]

    best_parameters = all_parameters[np.where(all_solutions == np.nanmin(all_solutions))[0][0]]
else:
    best_parameters = np.zeros(9)

print('Best', best_parameters)
popsize = 16
es = cma.CMAEvolutionStrategy(best_parameters, 0.1, {'verb_disp': 1, 'popsize': popsize})


while not es.stop():

    #If CMA data already exists, can just use that until it's all caught up? No, apparently not... Can use it as a starting point though.

    parameters = es.ask()
    #Save out the parameter sets into temporary files here. Make sure to remove any previous ones as that could cause confusion
    for i in range(len(parameters)):
       if os.path.exists('./data/cma_data/parameters_%05d_%05d.txt' % (batch_id,i)):
            os.remove('./data/cma_data/parameters_%05d_%05d.txt' % (batch_id,i))
       if os.path.exists('./data/cma_data/results_%05d_%05d.txt' % (batch_id,i)):
            os.remove('./data/cma_data/results_%05d_%05d.txt' % (batch_id,i))
       np.savetxt('./data/cma_data/parameters_%05d_%05d.txt' % (batch_id,i), parameters[i], delimiter = ',')

    solutions = await_solutions(batch_id, popsize)

    #Save out to a CMA data file
    if os.path.exists(directory_fname):
        #This directory already exists. Hopefully with proper header information etc
        directory_data = []
        with open(directory_fname, "r", encoding="utf-8") as f:
            data = csv.reader(f)
            for row in data:
                directory_data.append(row)
    else:
        directory_data = []

    for i in range(len(solutions)):
        if np.isnan(solutions[i]):
            solutions[i] = 1e10

        snap_id = len(directory_data)
        new_row_data = [snap_id, solutions[i]] + list(parameters[i])
        directory_data.append(new_row_data)

    with open(directory_fname, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(directory_data)

    es.tell(parameters, solutions)
