from datetime import datetime
import os
import drms
import numpy as np
from scipy.io import netcdf_file
import astropy.units as u
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.stats import wasserstein_distance

class VelocityNet():
    """
    Class for the neural net used for determining velocities from expansion factors and coronal hole distances (and maybe others?).
    If possible, put all the weighting functions and things in here
    """


    def __init__(self, n_nodes=10, n_inputs=2, max_velocity=1000, min_velocity=0):
        #Generate arraysfor the weights and things. Then broadcast to these arrays so the dimensions don't get mixed up...
        self.n_nodes = n_nodes
        self.n_inputs = n_inputs
        self.nparas = n_inputs*n_nodes + n_nodes*2 + 1
        self.max_velocity = max_velocity
        self.min_velocity = min_velocity

        self.weights_in = np.zeros((n_inputs, n_nodes))
        self.biases_in = np.zeros((n_nodes))
        self.weights_out = np.zeros((n_nodes))
        self.biases_out = np.zeros((1))
        self.parameter_set = np.zeros((self.nparas))

    def generate_random_seed(self, sigma=0.1):
        """
        Does a seed with Gaussian weightings of sigma
        """
        self.parameter_set = np.random.normal(scale=sigma, size = self.nparas)
        self.weights_in[:,:] = np.reshape(self.parameter_set[0:self.n_inputs*self.n_nodes], shape = np.shape(self.weights_in))
        self.biases_in[:] = self.parameter_set[self.n_inputs*self.n_nodes:self.n_inputs*self.n_nodes+self.n_nodes]
        self.weights_out[:] = self.parameter_set[self.n_inputs*self.n_nodes+self.n_nodes:self.n_inputs*self.n_nodes+self.n_nodes*2]
        self.biases_out[:] = self.parameter_set[self.n_inputs*self.n_nodes+self.n_nodes*2:self.n_inputs*self.n_nodes+self.n_nodes*2+1]

    def update_network(self, parameter_set):
        """
        For a given parameter set, updates the network biases etc.
        """
        self.parameter_set = parameter_set.copy()
        self.weights_in[:,:] = np.reshape(parameter_set[0:self.n_inputs*self.n_nodes], shape = np.shape(self.weights_in))
        self.biases_in[:] = parameter_set[self.n_inputs*self.n_nodes:self.n_inputs*self.n_nodes+self.n_nodes]
        self.weights_out[:] = parameter_set[self.n_inputs*self.n_nodes+self.n_nodes:self.n_inputs*self.n_nodes+self.n_nodes*2]
        self.biases_out[:] = parameter_set[self.n_inputs*self.n_nodes+self.n_nodes*2:self.n_inputs*self.n_nodes+self.n_nodes*2+1]

    def sigmoid(self, x):
        """
        Returns the sigmoid function of the input x
        """
        return 1.0/(1.0 + np.exp(-x))

    def velocity(self, inputs):
        """
        For given input arrays, runs the neural net to find the expected output (number between 0 and 1)
        """
        inputs = np.array(inputs)
        node_activations = np.zeros((self.n_nodes, np.shape(inputs[0])[0], np.shape(inputs[0])[1]))

        for n in range(self.n_nodes):
            node_activations[n] += np.sum(self.weights_in[:,n, np.newaxis, np.newaxis]*inputs[:], axis=0) + self.biases_in[n]
        node_activations = np.clip(node_activations, a_min = 1e-3, a_max = 1e3)  #Stop over and underflow in the exponentials
        node_activations = self.sigmoid(node_activations)
        output = np.sum(node_activations*self.weights_out[:,np.newaxis,np.newaxis], axis=0) + self.biases_out
        output = self.sigmoid(output)

        return self.min_velocity + output*(self.max_velocity - self.min_velocity)

    def save_current_state(self, minimiser, fname='net_state.txt'):
        """
        Appends to the filename the current ability of the net in question. SO can restart a training run with impunity etc.
        """
        step = 0
        if os.path.exists(fname):
            with open(fname) as f:
                for line in f.readlines():
                    step += 1

        save_line = [step, minimiser] + self.parameter_set.tolist()
        with open(fname, "a") as f:
            f.write(" ".join(f"{x:.6f}" for x in save_line) + "\n")
        return

    def load_best_state(self, fname='net_state.txt', override_nnodes=False):

        #Determine the correct number of parameters for this best state
        if os.path.exists(fname):
            print('Using bespoke best state')
            best_score = 1e6; best_id = 0
            best_parameters = []
            with open(fname, "r") as f:
                data = f.readlines()
                cut = len(data)
                for li, line in enumerate(data[-cut:]):
                    if float(line.split(' ')[1]) < best_score:
                        best_score = float(line.split(' ')[1])
                        best_id = len(data) - cut + li

            for val in data[best_id].split(' ')[2:]:
                best_parameters.append(float(val))
            print('Best score', best_score)
        elif os.path.exists('./nets/default.txt'):
            print('Using default state')
            best_score = 1e6; best_id = 0
            best_parameters = []
            with open('./nets/default.txt', "r") as f:
                data = f.readlines()
                cut = len(data)
                for li, line in enumerate(data[-cut:]):
                    if float(line.split(' ')[1]) < best_score:
                        best_score = float(line.split(' ')[1])
                        best_id = len(data) - cut + li
            for val in data[best_id].split(' ')[2:]:
                best_parameters.append(float(val))
        else:
            raise Exception('Log file not found...')

        target_nparas = len(best_parameters)

        n_best = 0
        for n in range(100):
            if target_nparas == self.n_inputs*n + n*2 + 1:
                n_best = n
                print(f'Best parameters are for {n} nodes')
                break

        if n_best != self.n_nodes:
            print("Best parameters are not for the correct amount of nodes")

        if n_best < self.n_nodes and override_nnodes:
            print('Using fewer nodes and will extend as appropriate in due course')
            self.n_nodes = n_best
            n_nodes = n_best
            n_inputs = self.n_inputs
            self.nparas = n_inputs*n_nodes + n_nodes*2 + 1

            self.weights_in = np.zeros((n_inputs, n_nodes))
            self.biases_in = np.zeros((n_nodes))
            self.weights_out = np.zeros((n_nodes))
            self.biases_out = np.zeros((1))
            self.parameter_set = np.zeros((self.nparas))

        elif n_best == self.n_nodes:
            print('Best scores match. Carry on.')
        else:
            raise Exception("Cannot use best parameters. Sort it out.")

        self.parameter_set[:] = np.array(best_parameters)
        self.weights_in[:,:] = np.reshape(self.parameter_set[0:self.n_inputs*self.n_nodes], shape = np.shape(self.weights_in))
        self.biases_in[:] = self.parameter_set[self.n_inputs*self.n_nodes:self.n_inputs*self.n_nodes+self.n_nodes]
        self.weights_out[:] = self.parameter_set[self.n_inputs*self.n_nodes+self.n_nodes:self.n_inputs*self.n_nodes+self.n_nodes*2]
        self.biases_out[:] = self.parameter_set[self.n_inputs*self.n_nodes+self.n_nodes*2:self.n_inputs*self.n_nodes+self.n_nodes*2+1]
        return

    def load_latest_state(self, fname='net_state.txt'):
        """
        Loads the last state found in the log, not necessarily the best
        """
        if os.path.exists(fname):
            best_score = 1e6; best_id = 0
            best_parameters = []
            with open(fname, "r") as f:
                data = f.readlines()
                cut = len(data)
                line = data[-1]

            for val in data[-1].split(' ')[2:]:
                best_parameters.append(float(val))
            print('Latest score', best_score)
        else:
            raise Exception('Log file not found...')

        self.parameter_set[:] = np.array(best_parameters)
        self.weights_in[:,:] = np.reshape(self.parameter_set[0:self.n_inputs*self.n_nodes], shape = np.shape(self.weights_in))
        self.biases_in[:] = self.parameter_set[self.n_inputs*self.n_nodes:self.n_inputs*self.n_nodes+self.n_nodes]
        self.weights_out[:] = self.parameter_set[self.n_inputs*self.n_nodes+self.n_nodes:self.n_inputs*self.n_nodes+self.n_nodes*2]
        self.biases_out[:] = self.parameter_set[self.n_inputs*self.n_nodes+self.n_nodes*2:self.n_inputs*self.n_nodes+self.n_nodes*2+1]
        return

    def extend_net(self, n_nodes_target):
        """
        Adds more neurons to an existing net, up to n_nodes_target of them.
        """
        nparas_new = self.n_inputs*n_nodes_target + n_nodes_target*2 + 1
        new_parameter_set = np.zeros(nparas_new)

        new_weights = np.zeros((self.n_inputs, n_nodes_target))
        new_weights[:,:self.n_nodes] = self.weights_in

        new_parameter_set[0:self.n_inputs*n_nodes_target] = np.reshape(new_weights, shape = (self.n_inputs*n_nodes_target))
        new_parameter_set[self.n_inputs*n_nodes_target:self.n_inputs*n_nodes_target+self.n_nodes] = self.biases_in[:]
        new_parameter_set[self.n_inputs*n_nodes_target+n_nodes_target:self.n_inputs*n_nodes_target+n_nodes_target+self.n_nodes] = self.weights_out[:]
        new_parameter_set[self.n_inputs*n_nodes_target+n_nodes_target*2:self.n_inputs*n_nodes_target+n_nodes_target*2+1] = self.biases_out[:]

        #Update Net metadata
        self.parameter_set = new_parameter_set
        self.n_nodes = n_nodes_target
        self.nparas = self.n_inputs*self.n_nodes + self.n_nodes*2 + 1

        self.weights_in = np.zeros((self.n_inputs, n_nodes_target))
        self.biases_in = np.zeros((n_nodes_target))
        self.weights_out = np.zeros((n_nodes_target))
        self.biases_out = np.zeros((1))

        self.weights_in[:,:] = np.reshape(self.parameter_set[0:self.n_inputs*self.n_nodes], shape = (self.n_inputs, self.n_nodes))
        self.biases_in[:] = self.parameter_set[self.n_inputs*self.n_nodes:self.n_inputs*self.n_nodes+self.n_nodes]
        self.weights_out[:] = self.parameter_set[self.n_inputs*self.n_nodes+self.n_nodes:self.n_inputs*self.n_nodes+self.n_nodes*2]
        self.biases_out[:] = self.parameter_set[self.n_inputs*self.n_nodes+self.n_nodes*2:self.n_inputs*self.n_nodes+self.n_nodes*2+1]

        return


def get_wasserstein_distance(speeds1, speeds2, nbins=101, doplots=False, huxt_name=None, iteration=0):
    """
    Given two (aligned) distributions of speeds, returns the Wasserstein distance between them.
    Will bin everything between 0 and 1000 km/s. nbins is to be determined empirically?
    """
    hist1, _ = np.histogram(speeds1, bins=nbins, range=(0.0,1000.0))
    hist2, _ = np.histogram(speeds2, bins=nbins, range=(0.0,1000.0))
    hist1 = hist1/np.sum(hist1)
    hist2 = hist2/np.sum(hist2)
    distance = wasserstein_distance(hist1, hist2)*1000
    if doplots:
        fig = plt.figure(figsize=(10,7))
        plt.plot(np.linspace(0,1000,len(hist1)), hist1)
        plt.plot(np.linspace(0,1000,len(hist2)), hist2)
        plt.title(f'Wasserstein Distance: {distance}')
        plt.tight_layout()
        plt.savefig('./plots/%s/hists_%05d.png' % (huxt_name, iteration))
        print(f'Plot saved to {'./plots/%s/hists_%05d.png' % (huxt_name, iteration)}')
        plt.close()

    return distance

def get_average_speeds(alltimes, allspeeds, spinup_time = 0, cadence=10, weighted_average = True, verbose=False, plot_averaging=False, target_times=None):
    """
    Using the dtimes, speeds and parameters for spinup etc., find the average values over the specified forecast interval.
    Should be able to deal with missing data, but hopefully everything should be on the same cadence so that's fine...

    This one runs through in the other dimension, so hopefully should be considerably faster

    Cadence is the number of outputs per day
    """

    #Run through each time which exists, check all of them off and take means from there.
    #Refactor so you get a list for each time? Yes. But deal with it hourly (round down to nearest hour)

    if target_times is None:
        mins = [min(i for i in row) for row in alltimes]
        maxs = [max(i for i in row )for row in alltimes]

        tmin = min(mins)
        tmax = max(maxs)

        #Convert these (and all) the times to seconds and add on the spinup time
        tmin = int(tmin.timestamp())  + 3600*24*spinup_time
        tmax = int(tmax.timestamp())

        #Figure out the (hourly) times at which speeds can be logged
        potential_times = np.arange(tmin, tmax, step = int(3600*24/cadence))

    else:
        potential_times = target_times.astype("datetime64[s]").astype(np.int64)

    allvalues = [[] for _ in range(len(potential_times))]
    allprops = [[] for _ in range(len(potential_times))]

    if plot_averaging:
        fig = plt.figure(figsize=(10,7))

    #If there's only one set, this can be sped up quite considerably. It's only when using huxt inputs that one needs to do it the slow way
    if len(alltimes) == 1:
        if verbose:
            print('Doing direct interpolation (for OMNI etc.)')
        set_num = 0
        times_secs = alltimes[set_num].astype("datetime64[s]").astype(np.int64)
        interp = interp1d(times_secs, allspeeds[set_num], bounds_error = False, fill_value= "extrapolate")
        for ti, check_time in enumerate(potential_times):
            allvalues[ti] = interp(check_time)
        mean_values = allvalues
    else:

        for set_num in range(len(alltimes)):
            if verbose:
                print(f'Averaging set {set_num} of {len(alltimes)}...')
            #Run through each INDIVIDUAL set
            times_secs = alltimes[set_num].astype("datetime64[s]").astype(np.int64)
            interp = interp1d(times_secs, allspeeds[set_num], bounds_error = True)

            if plot_averaging:
                plt.plot(times_secs, allspeeds[set_num], linewidth = 0.5, c = 'black')
            #Then run through the times which may lie in the range

            #I think this can be sped up considerably by using searchsorted. It's far too slow as things stand
            start_check_index = np.searchsorted(potential_times, np.min(times_secs))

            for ti, check_time in enumerate(potential_times[start_check_index:]):
                if verbose:
                    print(f'Progress: {100*ti/len(potential_times[start_check_index:])}%')

                inrange = True
                if check_time < min(times_secs):
                    inrange = False
                    continue
                if check_time > max(times_secs):
                    inrange = False
                    break  #Need to check this is the correct functionality
                if check_time < times_secs[0] + 3600*24*spinup_time:
                    inrange = False
                    continue

                if inrange:
                    #This check time is in the range. Add on the value*prop, and the prop itself
                    frac = (check_time - times_secs[0] - 3600*24*spinup_time)/(max(times_secs) - times_secs[0] - 3600*24*spinup_time)
                    if weighted_average:
                        prop = 1 - (2*frac - 1)**2
                    else:
                        prop = 1

                    allvalues[start_check_index + ti].append(interp(check_time)*prop)
                    allprops[start_check_index + ti].append(prop)

        #Take the mean for each of the above times
        mean_values = np.zeros(len(potential_times))
        for ti, time in enumerate(potential_times):
            if len(allvalues[ti]) > 0 and np.sum(allprops[ti]) > 0:
                mean_values[ti] = np.sum(allvalues[ti])/np.sum(allprops[ti])
            else:
                mean_values[ti] = np.nan

    final_times = np.array(potential_times).astype("datetime64[s]")
    mean_values = np.array(mean_values)

    if plot_averaging:
        plt.plot(potential_times, mean_values, c = 'green')
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.show()
        plt.close(fig)

    del potential_times, interp, times_secs, allspeeds
    return final_times, mean_values

def fast_fractions(vs, velocity_threshold = 500, cadence_days = 5):
    """
    Does the stats on the fraction of 'fast' wind speed times for each of the models.
    Do need to adjust for nans, probably.
    """
    cadence_int = int(cadence_days*24)
    cut_start = 0
    cut_end = cut_start + cadence_int
    props = []
    while cut_start < len(vs):
        if (np.count_nonzero(~np.isnan(vs[cut_start:cut_end]))) > 0:
            prop = np.sum(vs[cut_start:cut_end] > velocity_threshold) / (np.count_nonzero(~np.isnan(vs[cut_start:cut_end])))
        else:
            prop = np.nan
        props.append(prop)
        cut_start = cut_start + cadence_int
        cut_end = cut_start + cadence_int

    return np.array(props)

def do_met_stats(times_omni, model_speeds, reference_speeds, compare_to_persist=True, persistence_cadence=int(24*27.7), thresholds = np.arange(450,600,5), cut_start=None, cut_end=None):
    """
    Do the met office stats -- the proportion of fast solar wind in 5 day chunks and comparison to persistence model (if appropriate)
    Doesn't do any of the slow calculations, but does require the times and vs to be appropriately saved
    """

    if compare_to_persist:
        persistence_int = persistence_cadence
        vs_persist = reference_speeds.copy()
        vs_persist[persistence_int:] = vs_persist[:-persistence_int]
        vs_persist[:persistence_int] = np.nan


    "#The persistence forecast is meaningless for the first month anyway, so just don't calculate these ones"

    threshold_scores = np.zeros(len(thresholds))

    for ti, threshold in enumerate(thresholds):
        velocity_threshold = threshold

        if cut_start is not None and cut_end is not None:
            model_fractions = fast_fractions(model_speeds[cut_start:cut_end], velocity_threshold = velocity_threshold)
            base_fractions = fast_fractions(reference_speeds[cut_start:cut_end], velocity_threshold = velocity_threshold)
            if compare_to_persist:
                persist_fractions = fast_fractions(vs_persist[cut_start:cut_end], velocity_threshold=velocity_threshold)
        else:
            model_fractions = fast_fractions(model_speeds, velocity_threshold = velocity_threshold)
            base_fractions = fast_fractions(reference_speeds, velocity_threshold = velocity_threshold)
            if compare_to_persist:
                persist_fractions = fast_fractions(vs_persist, velocity_threshold=velocity_threshold)

        if compare_to_persist:
            if np.sum(~np.isnan(persist_fractions)) == 0:
                skillscore = np.nan
            else:
                if np.nanmean((persist_fractions - base_fractions)**2) > 0.0:
                    skillscore = 1.0 - np.nanmean((model_fractions - base_fractions)**2)/np.nanmean((persist_fractions - base_fractions)**2)
                else:
                    skillscore = np.nan
        else:
            skillscore = 1.0 - np.nanmean((model_fractions - base_fractions)**2)

        threshold_scores[ti] = skillscore

    return threshold_scores



