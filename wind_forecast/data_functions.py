#This module is to contain all the functions pertaining to obtaining or wrangling the raw data

from datetime import datetime, timedelta
import os, sys
import drms
import numpy as np
from scipy.io import netcdf_file
import astropy.units as u

import csv

def get_PFSS_maps_local(br_map, vr_map, phi, cotheta):
    """
    Transforms saved data into the correct units etc. to be used by HuxT. Avoids any reading/writing.
    """
    phi = phi * u.rad
    theta = (np.pi / 2 - np.arccos(cotheta)) * u.rad
    vr_lats = theta[:, 0]
    br_lats = vr_lats
    vr_longs = phi[0, :]
    br_longs = vr_longs

    return vr_map, vr_longs, vr_lats, br_map, br_longs, br_lats

def get_cme_times(obs_times, cme_fname='./data/shared_data/cme_list.csv', cme_time_window_days=5, also_filter_persistence=False, persistence_cadence=27.7):
    """
    Given a set of valid snap times, determine whether these correspond to CME arrivals (up to the 5 day bodge period).
    This is more complicated than it seems!
    Don't use the precomputed OMNI file as this is quick enough just to do again. Produces a mask for the allowable subset of times to sample for an optimisation.
    """
    print('Filtering for CMEs')
    if os.path.exists(cme_fname):
        cme_data = []
        with open(cme_fname, "r", encoding="utf-8") as f:
            data = csv.reader(f)
            for row in data:
                try:
                    cme_data.append([datetime.strptime(row[1], "%Y/%m/%d %H%M"), datetime.strptime(row[2], "%Y/%m/%d %H%M")])
                except:
                    pass
                #cme_data.append([row[1],row[2]])
    else:
        raise Exception('CME List not found')

    #This is just a list of the start and end times of each CME.
    cme_data = np.array(cme_data)
    min_i = 0
    cme_flags = np.zeros(len(obs_times))

    #Logs the times of the CMES
    for ci in range(len(cme_data)):
        cme_start = cme_data[ci][0] - timedelta(days=cme_time_window_days)
        cme_end = cme_data[ci][1] + timedelta(days=cme_time_window_days)
        iscme = 0
        i = min_i
        go = True
        while go and i < len(obs_times):
            if cme_start < obs_times[i] and cme_end > obs_times[i]:
                cme_flags[i] = 1
                if iscme == 0:
                    min_i = i
                iscme = 1
            if obs_times[i] > cme_end:
                break
            i += 1
    if also_filter_persistence:
        cme_persist_flags = np.zeros(len(obs_times))

        min_i = 0

        #Logs the times of the CMEs, offset by one solar rotation (into the future)
        for ci in range(len(cme_data)):
            cme_start = cme_data[ci][0] - timedelta(days=cme_time_window_days) + timedelta(days=persistence_cadence)
            cme_end = cme_data[ci][1] + timedelta(days=cme_time_window_days) + timedelta(days=persistence_cadence)
            iscme = 0
            i = min_i
            go = True
            while go and i < len(obs_times):
                if cme_start < obs_times[i] and cme_end > obs_times[i]:
                    cme_persist_flags[i] = 1
                    if iscme == 0:
                        min_i = i
                    iscme = 1
                if obs_times[i] > cme_end:
                    break
                i += 1

        cme_flags = np.maximum(cme_flags, cme_persist_flags)

    return np.array(cme_flags).astype('int')


def load_chb_distances(run_name, snap_id):
    """
    Loads in the coronal boundary distances and expansion factors. Just for plotting etc. for now I think.
    I'd like to avoid all readin/writing after this point for speed reasons, but that may be tricky


    Parameters
    ----------
    run_name : string
        Name of the batch, used for the folder to which the upper boundaries are saved
    snap_id: int
        Number of the time snap for this particular calculation

    Returns
    -------
    distances : array
        Distances to coronal hole boundaries
    """
    path = os.getcwd() + '/' + "data" + '/' + run_name
    snap_fname = 'chb_%09d.nc' % snap_id

    chb_fname = path + '/' + snap_fname

    if not os.path.exists(chb_fname):
        raise Exception('Coronal hole boundary file not found')
    fid = netcdf_file(chb_fname, "r")
    s0 = fid.variables["cos(th)"][:].copy()
    ph0 = fid.variables["ph"][:].copy()
    br0 = fid.variables["br"][:].copy()
    fs = fid.variables["expansionFactor"][:].copy()
    chd = fid.variables["CHBDistance"][:].copy()
    fid.close()

    return s0, ph0, br0, fs, chd

def get_source_times(src_folder):
    '''
    #This gets the available GONG input data times
    '''
    times = []
    path = os.getcwd()

    if not os.path.exists(path + "/" + src_folder):
        raise Exception(f'Data folder not found at {path + "/" + src_folder}. Check root directory location?')

    for fname in os.listdir(path + "/" + src_folder):
        if fname.startswith("evo.Earth"):
            times.append(datetime.strptime(fname[-13:-3], "%Y%m%d%H"))

    return sorted(times)

def get_crot_times(data_dir=None):
    """
    Obtain the Carrington rotation times. To be used for modifying the HMI/MDI maps to be in line with the GONG maps
    Uses the code I wrote for outflowpy
    """

    if data_dir is not None:
        if os.path.exists(f'{data_dir}/crot_times.npy'):
            try:
                centre_times = np.load(f'{data_dir}/shared_data/crot_times.npy', allow_pickle=True)
                return centre_times
            except:
                print("Carrington rotation data not found, trying to download it...")

    try:
        c = drms.Client()
        #Find the correct Carrington Rotation for this date.
        crot_times_mdi = c.query(('mdi.synoptic_mr_polfil_96m'), key = ["T_START","T_STOP","CAR_ROT"])
        crot_times_hmi = c.query(('hmi.synoptic_mr_polfil_720s'), key = ["T_START","T_STOP","CAR_ROT"])
    except:
        raise Exception("Failed to find the Carrington Rotation database")

    start_times = []
    end_times = []
    centre_times = []

    for source, crot_times in enumerate([crot_times_mdi, crot_times_hmi]):
        start_times_raw = list(crot_times.pop("T_START"))
        for i in range(len(start_times_raw)):
            if start_times_raw[i][-6:-4] == "60":
                start_times_raw[i] = start_times_raw[i][:-6] + "00" + start_times_raw[i][-4:]
        end_times_raw = list(crot_times.pop("T_STOP"))
        for i in range(len(end_times_raw)):
            if end_times_raw[i][-6:-4] == "60":
                end_times_raw[i] = end_times_raw[i][:-6] + "00" + end_times_raw[i][-4:]

        for si, crot_number in enumerate(crot_times.pop("CAR_ROT")):
            if (crot_number < 2098 and source == 0) or (crot_number >= 2098 and source == 1):
                start_times.append(datetime.strptime(start_times_raw[si].split('_TAI')[0], "%Y.%m.%d_%H:%M:%S"))
                end_times.append(datetime.strptime(end_times_raw[si].split('_TAI')[0], "%Y.%m.%d_%H:%M:%S"))
                centre_times.append(0.5*(start_times[-1] - end_times[-1]) + start_times[-1])

    if data_dir is not None:
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
        if not os.path.exists(f'{data_dir}/shared_data'):
            os.mkdir(f'{data_dir}/shared_data')

        np.save(f'{data_dir}/shared_data/crot_times.npy', centre_times)

    return np.array(centre_times)

def obtain_enlil_data(run, target_times = []):
    """
    Obtains the met office enlil predicions and saves out in the same format as those generated by my PFSS or outflow calculations
    """

    path = os.getcwd()
    wsapath = path + "/GONG/"
    all_dtime_wsa = []
    all_vsw_wsa = []

    obs_times = get_observation_times("GONG")

    for time in obs_times[:]:
        wsafile = time.strftime("evo.Earth.%Y%m%d%H.nc")
        fh = netcdf_file(wsapath+wsafile, "r", mmap=False)
        t_wsa = fh.variables["TIME"][:]
        v1_wsa = fh.variables["V1"][:]
        fh.close()
        dtime_wsa0 = datetime.datetime.strptime(wsafile, "evo.Earth.%Y%m%d%H.nc")
        dtime_wsa1 = np.array([dtime_wsa0 + datetime.timedelta(seconds=int(t)) for t in t_wsa])
        # ind = [(t > dtime_wsa0) for t in dtime_wsa1]
        all_dtime_wsa.append(dtime_wsa1[dtime_wsa1 > dtime_wsa0])
        all_vsw_wsa.append(v1_wsa[dtime_wsa1 > dtime_wsa0]/1e3)

    dtime_wsa, vsw_wsa = get_average_speeds(all_dtime_wsa, all_vsw_wsa, spinup_time = 0.0, target_times = target_times)

    dtime_wsa = np.array(dtime_wsa)
    vsw_wsa = np.array(vsw_wsa)

    np.savetxt((path + '/data/' + run + '/times_enlil.txt'), dtime_wsa.astype("datetime64[s]").astype(str), fmt="%s")
    np.savetxt((path + '/data/' + run + '/vs_enlil.txt'), vsw_wsa)

    return dtime_wsa, vsw_wsa

def obtain_omni_data(run, target_times, overwrite=False):
    """
    Downloads the OMNI data and saves in the standard format
    """
    redo = False
    if not os.path.exists('./data/' + run + '/times_omni.txt') or not os.path.exists('./data/' + run + '/vs_omni.txt') or overwrite:
        redo = True
    else:
        #Check time data matches
        dtime_omni = np.loadtxt('./data/' + run + '/times_omni.txt', dtype='datetime64[s]', delimiter = ',')
        if dtime_omni[0] != target_times[0] or dtime_omni[-1] != target_times[-1]:
            redo = True

    if redo:
        print('Downloading OMNI data')
        path = os.getcwd()

        dtime_min = target_times[0]
        dtime_max = target_times[-1]

        data_omni = Hin.get_omni(dtime_min, dtime_max)
        all_dtime_omni = [data_omni['datetime']]
        all_vsw_omni = [data_omni['V'].values]

        dtime_omni, vsw_omni = fcast.get_average_speeds(all_dtime_omni, all_vsw_omni, spinup_time = 0.0, cadence=24, target_times=target_times, verbose=True)

        np.savetxt(('./data/' + run + '/times_omni.txt'), dtime_omni.astype("datetime64[s]").astype(str), fmt="%s")
        np.savetxt(('./data/' + run + '/vs_omni.txt'), vsw_omni)

    else:
        dtime_omni = np.loadtxt('./data/' + run + '/times_omni.txt', dtype='datetime64[s]', delimiter = ',')
        vsw_omni = np.loadtxt('./data/' + run + '/vs_omni.txt', delimiter = ',')

    return dtime_omni, vsw_omni

def update_theta_record(test_parameters, best_loss, sigma, best_theta):
    """
    This should be a relatively simple one to update a theta record during an optimisation run.
    Will hopefully then allow for fully automatic evalulation of the ability of the converged parameters.
    """

    if not os.path.exists('data'):
        os.mkdir('data')
    if not os.path.exists(f'data/{test_parameters["run_name"]}'):
        os.mkdir(f'data/{test_parameters["run_name"]}')

    directory_fname = f'./data/{test_parameters["run_name"]}/log.csv'
    if os.path.exists(directory_fname):
        #This directory already exists. Hopefully with proper header information etc
        directory_data = []
        with open(directory_fname, "r", encoding="utf-8") as f:
            data = csv.reader(f)
            for row in data:
                directory_data.append(row)
    else:
        header_row = [test_parameters["base_name"], test_parameters["optimisation_type"]]
        directory_data = [header_row]

    snap_id = len(directory_data) - 1
    new_row_data = [snap_id, best_loss, sigma]
    for i in range(len(best_theta)):
        new_row_data.append(best_theta[i])
    directory_data.append(new_row_data)

    with open(directory_fname, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(directory_data)

    print(f'Record updated for iteration {snap_id}')
    return True

def update_directory(update_type, fname, snap_id, args):
    """
    Updates the directory for the base directory, to check whether things need to be redone or not.
    Exact formatting etc. needs to be determined, but a .csv is probably the best way forward?
    """
    directory_fname = f'./data/{fname}/directory.csv'
    if update_type == "base":
        if os.path.exists(directory_fname):
            directory_data = []
            with open(directory_fname, "r", encoding="utf-8") as f:
                data = csv.reader(f)
                for row in data:
                    directory_data.append(row)
        else:
            directory_data = []

        new_row_data = ["base", snap_id, args[0], args[1], args[2], args[3], args[4][0], args[4][1], args[4][2]]
        #Check for an ID in the directory. If it exists, replace it. If not,add it.
        data_added = False
        for ri, row in enumerate(directory_data):
            if snap_id == int(row[1]):
                directory_data[ri] = new_row_data
                data_added = True
                break
        if not data_added:
            directory_data.append(new_row_data)
        with open(directory_fname, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(directory_data)
    elif update_type == "chbetc":
        if os.path.exists(directory_fname):
            directory_data = []
            with open(directory_fname, "r", encoding="utf-8") as f:
                data = csv.reader(f)
                for row in data:
                    directory_data.append(row)
        else:
            directory_data = []

        new_row_data = ["chbetc", snap_id, args[0]]
        #Check for an ID in the directory. If it exists, replace it. If not,add it.
        data_added = False
        for ri, row in enumerate(directory_data):
            if snap_id == row[1]:
                directory_data[row] = new_row_data
                data_added = True
                break
        if not data_added:
            directory_data.append(new_row_data)
        with open(directory_fname, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(directory_data)
    else:
        raise Exception('Update type not recognised.')


def check_existing_data(run_parameters, match_flag=True):
    """
    If a run name is given, check against the lookup table as to whether the lower boundary data exists for these parameters.
    Will raise appropriate errors, if appropriate.
    If not, will provide a database of what has been calculated for this 'name'

    Code:
    -1 for file exists but metadata is different. Flag if necessary.
    0 for file not existing
    1 for base model done
    2 for chbmap done
    """
    fname = run_parameters["base_name"]
    directory_fname = f'./data/{fname}/directory.csv'
    if os.path.exists(directory_fname):
        directory_data = []
        with open(directory_fname, "r", encoding="utf-8") as f:
            data = csv.reader(f)
            for row in data:
                directory_data.append(row)
    else:
        directory_data = []

    #Run through each of the desired snaps and see what exists. Need integer codes for this really. They are now defined above
    check_codes = [0]*len(run_parameters["observation_time"])
    for row in directory_data:
        if row[0] == "base":
            snap_id = int(row[1])
            if snap_id < len(check_codes):
                expected_row = ["base", str(snap_id), str(run_parameters["observation_time"][snap_id]), str(run_parameters["r_ss"]),
                                run_parameters["model_type"], run_parameters["data_source"],
                                str(run_parameters["resolutions"][0]), str(run_parameters["resolutions"][1]), str(run_parameters["resolutions"][2])]
                if row == expected_row and check_codes[snap_id] < 1:
                    check_codes[snap_id] = 1
                else:
                    check_codes[snap_id] = -1
                    if match_flag:
                        print("Existing parameter set:", row)
                        print("Specified parameter set", expected_row)
                        raise Exception("Existing data doesn't match these parameters in this run. Aborting... To overwrite with these new parameters put 'match_flag=False''")
                    else:
                        print("Existing parameter set:", row)
                        print("Specified parameter set", expected_row)
                        print("New parameters do not match existing ones, but proceeding anyway. Set 'match_flag=True' to stop this")
        elif row[0] == "chbetc":
            snap_id = int(row[1])
            if snap_id < len(check_codes):
                expected_row = ["chbetc", str(snap_id), str(run_parameters["r_hb"])]

                if row == expected_row and check_codes[snap_id] > -1:
                    check_codes[snap_id] = 2
                else:
                    check_codes[snap_id] = -1
                    if match_flag:
                        print("Existing parameter set:", row)
                        print("Specified parameter set", expected_row)
                        raise Exception("Existing data doesn't match these parameters in this run. Aborting... To overwrite with these new parameters put 'match_flag=False''")
                    else:
                        print("Existing parameter set:", row)
                        print("Specified parameter set", expected_row)
                        print("New parameters do not match existing ones, but proceeding anyway. Set 'match_flag=True' to stop this")
        else:
            raise Exception('Not done this bit yet')

    return check_codes
