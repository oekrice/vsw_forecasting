

import os
import numpy as np
from scipy.io import netcdf_file, FortranFile
import sys
import matplotlib.pyplot as plt
import datetime
from datetime import timedelta
import outflowpy
import sunpy
from scipy.interpolate import interp1d
import astropy.units as u
from sunpy.coordinates.sun import B0

this_directory = os.getcwd() + "/"
sys.path.append(this_directory +"prepare")

import pfss
import data_gong
import output_netcdf
import copy

sys.path.append(this_directory+"viz/tools")
import utils
import wind
sys.path.append(this_directory+'/viz/HUXt-master/code')
import huxt_inputs as Hin
import huxt as H
import huxt_analysis as HA

from astropy.time import Time

import wind_forecast as fcast
import astropy.units as u
import csv

def get_cme_fname(src_folder, tmatch):
    """
    Gets the file name of any relevant CME files
    """
    time_string = tmatch.strftime("%Y%m%d")
    fnames = []
    path = os.getcwd()
    for fname in os.listdir(path + "/" + src_folder):
        if fname.startswith("cone2bc_" + time_string):
            fnames.append(fname)
    if len(fnames) > 0:
        return path + "/" + src_folder + "/" + fnames[0]
    else:
        return None


def compute_vr(snap_id, run_name, method="wsa", params=[285, 625+285, 0.22222, 1, 0.8, 2, 2, 3], doplot=False):
    """
        Compute map of the solar wind speed v_r given the coronal hole boundary distance (chb, in degrees) and flux tube expansion factor (fs).
        ary -- 2019/09/13
    """

    s0, ph0, br0, fs, chd = fcast.load_chb_distances(run_name, snap_id)
    if params is None and method == "wsa_scaled":
        params = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    elif params is None and method == "wsa":
        params = [285, 625+285, 0.22222, 1, 0.8, 2, 2, 3]

    params = np.abs(params)

    fs = fs.copy()
    chb = chd.copy()


    if method == "marion":
        # e.g. params=[200.0, 700.0, 7.0, 1/2.5]
        vslow = params[0]
        vfast = params[1]
        w = params[2]
        d = params[3]
        vr = vslow + (vfast - vslow) * (np.deg2rad(chb) * w) ** d
        vr[vr < vslow] = vslow
        vr[vr > vfast] = vfast

    if method == "wsa":

        # e.g. params=[285, 625+285, 0.22222, 1, 0.8, 2, 2, 3]
        vslow = params[0]
        vfast = params[1]
        a = params[2]
        b = params[3]
        g = params[4]
        w = params[5]
        d = params[6]
        i = params[7]
        fs[fs < 0] = 1  # numerical error leading to negative fs
        vr = (
            vslow
            + ((vfast - vslow) / (1.0 + fs) ** a) * (b - g * np.exp(-(chb / w) ** d)) ** i
        )

    if method == "wsa_scaled":
        #This uses the limit data file to get each parameter while keeping them reasonable.
        scale_limits= np.loadtxt('./data/shared_data/wsa_limits.dat', delimiter = ',')
        def scale_parameter(i, x):
            return 0.5*(1.0 + np.tanh(x))*(scale_limits[i][1] - scale_limits[i][0]) + scale_limits[i][0]
        vslow = scale_parameter(0, params[0])
        vfast = scale_parameter(1, params[1])
        a = scale_parameter(2, params[2])
        b = scale_parameter(3, params[3])
        g = scale_parameter(4, params[4])
        w = scale_parameter(5, params[5])
        d = scale_parameter(6, params[6])
        i = scale_parameter(7, params[7])
        fs[fs < 0] = 1  # numerical error leading to negative fs
        vr = (
            vslow
            + ((vfast - vslow) / (1.0 + fs) ** a) * (b - g * np.exp(-(chb / w) ** d)) ** i
        )

    if method == "riley":
        # e.g. params=[0.1051, 0.0101, 333.0, 631.0]
        ep = params[0]
        w = params[1]
        vslow = params[2]
        vfast = params[3]
        vr = vslow + 0.5*(vfast - vslow) * (1.0 + np.tanh((np.deg2rad(chb) - ep) / w))

    if doplot:

        #Calculate meshgrid of reasonable speeds
        fss = np.linspace(0,300,200)
        chbs = np.linspace(0,30,200)

        fss, chbs = np.meshgrid(fss, chbs)
        vmesh = (vslow  + ((vfast - vslow) / (1.0 + fss) ** a) * (b - g * np.exp(-(chbs / w) ** d)) ** i)

        #Do a plot of the chbs, expansions, velocity map and resulting pattern
        fig, axs = plt.subplots(2,2, figsize = (10,7))

        im = axs[0,0].pcolormesh(fs, vmin = 0, vmax = np.percentile(fs, 99))
        axs[0,0].set_title('Expansion Factors')
        plt.colorbar(im, ax = axs[0,0])

        im = axs[0,1].pcolormesh(chb, vmin = 0, vmax = np.percentile(chb, 99))
        axs[0,1].set_title('Coronal Hole Boundary Distances')
        plt.colorbar(im, ax = axs[0,1])

        im = axs[1,0].pcolormesh(vmesh, vmin = 0, vmax = 1000)
        axs[1,0].set_title('Velocity function')
        plt.colorbar(im, ax = axs[1,0])

        im = axs[1,1].pcolormesh(vr, vmin = 0, vmax = np.percentile(vr, 99))
        axs[1,1].set_title('Velocity map')
        plt.colorbar(im, ax = axs[1,1])

        for i in range(2):
            for ax in axs[i]:
                ax.set_xticks([])
                ax.set_yticks([])
        plt.tight_layout()
        plt.show()
        #plt.savefig(f'./plots/vr_{run_name}_{snap_id}.png')
        plt.close()

    return vr

def compute_vr_net(snap_id, run_name, Net):
    """
        Compute map of the solar wind speed v_r given the coronal hole boundary distance (chb, in degrees) and flux tube expansion factor (fs).
        ary -- 2019/09/13
    """

    s0, ph0, br0, fs, chd = fcast.load_chb_distances(run_name, snap_id)

    vr = Net.velocity([fs, chd])


    return vr

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
                print('Entry updated')
                directory_data[row] = new_row_data
                data_added = True
                break
        if not data_added:
            print('Entry added')
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


def calculate_outflow(snap_id, obs_time, output_directory=None, overwrite=False, rss=2.5, is_pfss=False, save_snap=True, source=None, resolutions=[120,180,360]):
    """
    For a given observation time, calculates the PFSS field in the required format. Saves to an appropriate netcdf file.
    """
    path = os.getcwd() + '/' + "data"

    if save_snap and not output_directory:
        output_directory = "test"
        print(f"Data directory not provided, so saving outputs in folder 'test'. Set save_snap to False if you don't want to do this")

    if source is None or source not in ["gong", "hmi"]:
        raise Exception('Lower boundary source not recognised. Use wither "gong" or "hmi" (the latter of which uses MDI if necessary)')

    dir_fname = path + '/' + output_directory + '/'
    snap_fname = path + '/' + output_directory + '/' + 'outflow_%09d.nc' % snap_id

    if not os.path.isfile(snap_fname) or overwrite:  #If file doesn't exist, do a thing
        print('Computing Outflow model and wind speed at time', obs_time)

        dtime = obs_time
        nr = resolutions[0]
        ns = resolutions[1]
        nph = resolutions[2]
        rss = rss
        nsmooth = 7

        if source == "gong":
            print(f'Using GONG map at time {obs_time}')
            br0, t_actual, mapedge = data_gong.readhrlymap_gong(dtime, ns, nph, nsmooth=nsmooth)
            header = outflowpy.utils.carr_cea_wcs_header(dtime, br0.T.shape)
            input_map = sunpy.map.Map((br0, header))
            roll_data = input_map.data
        elif source == "hmi":
            print(f'Using HMI/MDI map at time {obs_time}')
            #Need to shift this data so it aligns with the GONG equivalent. Need to use the carrington rotation times for this

            crot_times = fcast.get_crot_times("data")
            crot_ind = np.searchsorted(crot_times, dtime)
            time_difference = np.min(dtime - crot_times)
            frac = (obs_time - crot_times[crot_ind -1])/(crot_times[crot_ind] - crot_times[crot_ind -1])
            input_map = outflowpy.obtain_data.prepare_hmi_mdi_time(dtime.strftime("%Y-%m-%dT%H:00:00"), ns, nph, smooth = 1.0*5e-2/nph, use_cached = True, cache_directory='./_download_cache/')
            input_map = sunpy.map.Map(np.roll(input_map.data, -int(frac*nph), axis=1), input_map.meta)

        # PLOT LOWER BOUNDARY MAP
        # fig = plt.figure()
        # plt.pcolormesh(input_map.data, vmin = -300, vmax = 300, cmap = 'seismic')
        # plt.colorbar()
        # plt.tight_layout()
        # plt.savefig(f'./plots/{source}_{obs_time}.png')
        # plt.close()

        if is_pfss:
            outflow_in = outflowpy.Input(input_map, nr, rss, mf_constant=0.0)
        else:
            outflow_in = outflowpy.Input(input_map, nr, rss)

        outflow_out = outflowpy.outflow(outflow_in)

        br, bs, bp = outflow_out.bcx
        bs = -bs  #This convention appears to differ from pfsspy and outflowpy, but is clearly important
        r = outflow_out.grid.rg; th = outflow_out.grid.sg; ph = outflow_out.grid.pg
        rrc = outflow_out.grid.rc; thc = outflow_out.grid.sc; pc = outflow_out.grid.pc

        #Manually extend the cordinate axes to include ghost points
        #And transform to match the original PFSS script
        rrc = np.zeros(len(outflow_out.grid.rc) + 2)
        thc = np.zeros(len(outflow_out.grid.sc) + 2)
        pc = np.zeros(len(outflow_out.grid.pc) + 2)
        rrc[1:-1] = outflow_out.grid.rc
        thc[1:-1] = outflow_out.grid.sc
        pc[1:-1] = outflow_out.grid.pc
        rrc[0] = 2*rrc[1] - rrc[2]; rrc[-1] = 2*rrc[-2] - rrc[-3]
        thc[0] = 2*thc[1] - thc[2]; thc[-1] = 2*thc[-2] - thc[-3]
        pc[0] = 2*pc[1] - pc[2]; pc[-1] = 2*pc[-2] - pc[-3]
        rrc = np.exp(rrc)
        r = np.exp(r)
        th = np.arccos(th)
        thc[1:-1] = np.arccos(thc[1:-1])
        thc[0] = -1; thc[-1] = -1

        if save_snap:
            os.makedirs(dir_fname, exist_ok=True)
            output_netcdf.bc(snap_fname, r, th, ph, rrc, thc, pc, br, bs, bp)
            if is_pfss:
                model = "pfss"
            else:
                model = "outflow"
            update_directory("base", output_directory, snap_id, [obs_time, rss, model, source, resolutions])  #This should update the log of what has been calculated already. Will be tricksy, I think, to make sure the data stays uncorrupted.

        return r, th, ph, rrc, thc, pc, br, bs, bp

    else:
        print('Outflow field already exists for this time. Use flag "overwrite=True" to recalculate.')
        return

def calculate_chb_exp(snap_id, batch_name, r_hb=21.5, overwrite=False, purge_data=False):
    """
    Calculates an outflow/PFSS field, the Schatten extension and saves the coronal hole distances and expansion factors
    Needs the appropriate base field to have been calculated already.

    Parameters
    ----------
    run_name : string
        Name of the batch, used for the folder to which the upper boundaries are saved
    snap_id: int
        Number of the time snap for this particular calculation
    obs_time : datetime object
        Observation time for the lower boundary
    is_pfss (optional) : bool
        If False, does an optimised outflow field. If True, does PFSS
    rss (optional) : float
        Source surface height
    r_hb (optional) : float
        Schatten extension limit

    Returns
    -------
    distances : array
        Distances to coronal hole boundaries
    """
    path = os.getcwd() + '/' + "data"
    base_fname = 'outflow_%09d.nc' % snap_id

    if not os.path.exists(path+'/'+batch_name+"/"+base_fname):
        raise Exception('Base model not yet calculated. Need to run compute_outflow or equivalent.')

    windmap_fname = path+'/'+batch_name+"/" + 'windmap_' + base_fname + '.unf'

    if not os.path.isfile(windmap_fname) or overwrite:  #If file doesn't exist, do a thing
        # '#Calculate the Schatten Extension/coronal holes etc. Requires the above file to be saved as the Fortran reads it in.'
        wind.windmap(base_fname, r_hb, path=path+'/'+batch_name+"/", codepath=os.getcwd() +'/viz/fortran/')
    else:
        print('Schatten file already calculated, so using that. Set overwrite=True to recalculate.')

    chb_fname = path + '/'+ batch_name+ "/" + "chb_" + ('%09d' % snap_id) + ".nc"

    if not os.path.isfile(chb_fname) or overwrite:  #If file doesn't exist, do a thing
        compute_coronal_hole_distances(snap_id, windmap_fname, batch_name, path=path+'/'+batch_name+"/")
        update_directory("chbetc", batch_name, snap_id, [r_hb])  #This should update the log of what has been calculated already. Will be tricksy, I think, to make sure the data stays uncorrupted.
    else:
        print('CHB/Expansion factor file already calculated, so using that. Set overwrite=True to recalculate.')


    if purge_data:
        #Remove the old files which aren't needed and are large. Stops hamilton getting bunged up.
        try:
            os.remove(f'./data/{batch_name}/{base_fname}')
        except:
            pass
        try:
            os.remove(windmap_fname)
        except:
            pass
        try:
            os.remove(f'./data/{batch_name}/chmap_{base_fname}.unf')
        except:
            pass
        try:
            os.remove(f'./data/{batch_name}/schat_{base_fname}')
        except:
            pass
    return

def get_vsw(snap_id, run_name, obs_time, vr_bnd, r_hb=21.5, fcast_length=5, decelerate=True, cone2bcfile='', plot2d=False, spinup_cme_days=0, include_cmes=False, cme_velocity_shift=1.0):
    """
    For given velocity map etc., calculate wind speed at Earth. Includes computing the vsw map on the surface
    - assumes a static heliospheric field
    - length of forecast is fcast_length days.
    - start inserting CMEs spinup_cme_days days before magnetogram time.

    If a cone2bcfile is provided, CMEs are inserted from it. (Format used for Enlil.)

    For testing, set plot2d=True to run HUXt in 2d and plot result instead of returning.
    """

    print('Calculating VSW with HuxT...')
    s0, ph0, br_bnd, fs, chd = fcast.load_chb_distances(run_name, snap_id)

    # Latitude of Earth at simulation time:
    _, vr_longs, vr_lats, br_map, br_longs, br_lats = fcast.get_PFSS_maps_local(br_bnd, vr_bnd, ph0, s0)
    E_lat = np.deg2rad(B0(obs_time))
    iE_lat = np.argmin(abs(vr_lats - E_lat))

    # Run HUXt:
    t_init = obs_time - datetime.timedelta(days=spinup_cme_days)
    cr, cr_lon_init = Hin.datetime2huxtinputs(t_init)
    vr_long = vr_bnd[iE_lat,:]*u.km/u.s
    br_long = br_bnd[iE_lat,:]

    # [optional] Decelerate to compensate for the fact that WSA is designed to work with ballistic mapping:
    if decelerate:
        vr_long, _ = Hin.map_v_inwards(vr_long, 215*u.solRad, vr_longs, r_hb*u.solRad)

    if plot2d:
        model = H.HUXt(v_boundary=vr_long, b_boundary=br_long, cr_num=cr, cr_lon_init=cr_lon_init, latitude=E_lat, simtime=(fcast_length + spinup_cme_days)*u.day, dt_scale=4, frame = 'synodic', r_min = r_hb*u.solRad)
    else:
        model = H.HUXt(v_boundary=vr_long, b_boundary=br_long, cr_num=cr, cr_lon_init=cr_lon_init, latitude=E_lat, lon_out=0.0*u.deg, simtime=(fcast_length + spinup_cme_days)*u.day, dt_scale=4, frame = 'synodic', r_min = r_hb*u.solRad)

    # Add in CMEs from conefile if required:
    # I'm not sure these were actually being used. Let's attempt...

    #Can obtain CMEs at a later date

    # if include_cmes:
    #     #Determine the CMEs which correspond to this input, if any.
    #     print('Including cone CMEs...')
    #     fname = get_cme_fname("GONG", obs_time)
    #     if fname is not None:
    #         cone2bcfile = fname
    #     else:
    #         cone2bcfile = ''
    #
    # if cone2bcfile != '':
    #
    #     #Get cone CMEs in the correct format, including the spinup offset time
    #     cme_list = Hin.ConeFile_to_ConeCME_list(model, cone2bcfile)
    #     cmes_out = []
    #     for cme in cme_list:
    #         cme_replace = copy.deepcopy(cme)
    #         cme_replace.v = cme.v*np.abs(cme_velocity_shift)
    #         if cme_replace.v > model.v_max:
    #             cme_replace.v = model.v_max - 1*(u.km/u.s)
    #         cme_replace.t_launch = cme.t_launch + spinup_cme_days*u.day
    #         if cme_replace.t_launch > 0.0*u.day:
    #             cmes_out.append(cme_replace)
    #
    #     cme_list = cmes_out
    #
    #     model.solve(cme_list, tag='cone_cme_test')
    #
    #     # # Manual CME for testing:
    #     # cme = H.ConeCME(t_launch=1*u.day, longitude=360*u.deg, latitude = 0*u.deg, width=60*u.deg, v=2500*(u.km/u.s), thickness=5*u.solRad, initial_height=r_hb*u.solRad)
    #     # model.solve([cme], tag='cone_cme_test')
    # else:

    model.solve([])

    if plot2d:
        HA.animate(model, tag='cone_cme_test')
        # t_interest = fcast_length//2*u.day
        # fig, ax = HA.plot(model, t_interest)
        # plt.show()
        sys.exit()

    # Extract time series at Earth:
    earth_series = HA.get_observer_timeseries(model, observer = 'Earth')

    dtime = earth_series['time']
    vsw = earth_series['vsw']
    bpol = earth_series['bpol']

    return dtime, vsw, bpol

def expansionfactor(rm, r0, brm, br0):
    """
        Compute map of flux tube expansion factor given br0 at r=r0 and brm at r=rm (radii) of each field line.
    """

    fs = brm * rm ** 2 / br0 / r0 ** 2

    return fs

def compute_coronal_hole_distances(snap_id, windmap_fname, batch_name, path='./', save_file=True):
    """
    This step should be included in the 'slow runs', as it requires no knowledge of the velocity maps.
    Should output an array of coronal hole distances, which can be saved and read-in quickly by the HuxT solver.

    Parameters
    ----------
    snap_id : int
        Snap number (for fname purposes)
    windmap_fname : string
        File name of the 'windmap'
    save_file (optional) : bool
        if True, saves out as an .nc

    Returns
    -------
    distances : array
        Distances to coronal hole boundaries
    """

    # Compute coronal hole boundary distances and flux tube expansion factors
    # -----------------------------------------------------------------------
    if not os.path.exists(windmap_fname):
        raise Exception("Windmap file doesn't exist at ", windmap_fname)

    fid = FortranFile(windmap_fname, "r")
    ns = fid.read_ints(dtype=np.int32)[0]
    nph = fid.read_ints(dtype=np.int32)[0]
    # Outer boundary of Schatten model:
    s0 = fid.read_reals(dtype=np.float64).reshape((ns, nph))
    ph0 = fid.read_reals(dtype=np.float64).reshape((ns, nph))
    br0 = fid.read_reals(dtype=np.float64).reshape((ns, nph))
    # Outer boundary of dumfric model:
    br1 = fid.read_reals(dtype=np.float64).reshape((ns, nph))
    # Photosphere:
    sm = fid.read_reals(dtype=np.float64).reshape((ns, nph))
    phm = fid.read_reals(dtype=np.float64).reshape((ns, nph))
    brm = fid.read_reals(dtype=np.float64).reshape((ns, nph))
    fid.close()

    # Compute map of coronal-hole boundary distances traced down through combined model:
    # - CHB distances:
    chd = sm * 0
    igood = (
        np.abs(sm) <= 1
    )  # points without bad mapping (U-shaped field lines or HCS)

    chd[igood] = chbmap_45(snap_id, sm[igood], phm[igood], path=path, codepath="./viz/fortran/")

    # Compute map of flux-tube expansion factors traced down through combined model:
    fs = sm * 0
    fs0 = expansionfactor(1, 2.5, brm[igood], br1[igood])
    fs[igood] = fs0

    # Fill in points with bad mapping using nearest values:
    for k in range(ns):
        for j in range(nph):
            if igood[k, j] == 0:
                angs = (s0[k, j] - s0) ** 2 + (ph0[k, j] - ph0) ** 2
                angs[~igood] = 999
                kk = np.unravel_index(np.argmin(angs), (ns, nph))
                chd[k, j] = chd[kk]
                fs[k, j] = fs[kk]

    if False:
        fig, axs = plt.subplots(2, figsize = (10,10))

        im = axs[0].pcolormesh(chd, vmin = 0, vmax = np.percentile(chd, 95))
        axs[0].set_title('Coronal Hole Distances')
        plt.colorbar(im, ax = axs[0])

        im = axs[1].pcolormesh(fs, vmin = 0, vmax = np.percentile(fs, 95))
        axs[1].set_title('Expansion Factors')
        plt.colorbar(im, ax = axs[1])

        plt.savefig(f'./plots/chb_exp_{batch_name}_{snap_id}.png')
        plt.close()


    snap = ('%09d' % snap_id) + ".nc"
    chb_fname = path + "chb_" + snap
    print('Saving hole distances to', chb_fname)

    fid = netcdf_file(chb_fname, "w")
    fid.createDimension("nth", ns)
    fid.createDimension("nph", nph)
    vid = fid.createVariable("cos(th)", "d", ("nth", "nph"))
    vid[:] = s0
    vid = fid.createVariable("ph", "d", ("nth", "nph"))
    vid[:] = ph0
    vid = fid.createVariable("br", "d", ("nth", "nph"))
    vid[:] = br0
    vid = fid.createVariable("expansionFactor", "d", ("nth", "nph"))
    vid[:] = fs
    vid = fid.createVariable("CHBDistance", "d", ("nth", "nph"))
    vid[:] = chd
    fid.nth = ns
    fid.nph = nph
    fid.close()

    print('Saved coronal hole distances and expansion factors to file:', chb_fname)

    return

def chbmap_45(snap, sm, pm, path="./", codepath="./fortran/", make_plot=False):

    #Check coronal hole map exists (should already have been calculated)
    ch_fname = os.path.join(path, "chmap_outflow_" + ('%09d' % snap) + ".nc.unf")
    if not os.path.exists(ch_fname):
        raise Exception('Coronal hole map not found at ', ch_fname)


    fid = FortranFile(ch_fname, "r")
    sc = fid.read_reals(dtype=np.float64)
    pc = fid.read_reals(dtype=np.float64)
    chmap = fid.read_ints(dtype=np.int32).reshape(
        (np.size(pc, 0) - 1, np.size(sc, 0) - 1)
    )
    fid.close()


    #### Reworking WSA 4.5 coronal hole boundary function into python for use with PFSS

    footpoint = np.abs(chmap)
    phiSize, thetaSize = footpoint.shape
    hole_bound_i = []
    hole_bound_j = []

    for step_i in range(phiSize): # loop across longitude
        for step_j in range(thetaSize): # loop across latitude
            if footpoint[step_i, step_j] == 0: # check if coordinate is closed field point

                # Get neighboring indices with wrapping and clamping
                i_l = (step_i - 1) % phiSize
                i_r = (step_i + 1) % phiSize
                j_b = max(step_j - 1, 0)
                j_t = min(step_j + 1, thetaSize - 1)

                boundary = 0

                # Check 3x3 neighborhood (excluding center)
                for step_k in range(j_b, j_t + 1):
                    step_h = i_l
                    h_stop = (i_r + 1) % phiSize

                    while True:
                        if footpoint[step_h, step_k] != 0:
                            boundary += 1
                        if step_h == i_r:
                            break
                        step_h = (step_h + 1) % phiSize

                if boundary > 0:
                    hole_bound_i.append(step_i)
                    hole_bound_j.append(step_j)

    # Include coronal hole pixels at poles
    for step_i in range(phiSize):
        if footpoint[step_i, 0] == 0:
            hole_bound_i.append(step_i)
            hole_bound_j.append(0)
        if footpoint[step_i, thetaSize - 1] == 0:
            hole_bound_i.append(step_i)
            hole_bound_j.append(thetaSize - 1)

    # Create array to store boundary coordinates
    ed = np.zeros(footpoint.shape)
    ed[hole_bound_i, hole_bound_j] = 1
    ed = ed.T

    ### Continue as normal ###
    s, p = np.meshgrid(
        0.5 * (sc[1:] + sc[:-1]), 0.5 * (pc[1:] + pc[:-1]), indexing="ij"
    )
    ped = p[ed == 1]
    sed = s[ed == 1]

    # - compute minimum spherical angle from each footpoint to edge list:
    #   (set to 0 for field lines with no photospheric mapping, e.g. open-open)
    chd = pm * 0
    sthm = np.sqrt(1 - sm ** 2)
    sthed = np.sqrt(1 - sed ** 2)
    for i in range(np.size(chd)):
        angs = np.real(np.arccos(sm[i] * sed + sthm[i] * sthed * np.cos(pm[i] - ped)))
        chd[i] = np.min(angs)

    # For plotting purposes, make an alternative map not of the footpoints

    if make_plot:
        fig, axs = plt.subplots(3,1)
        axs[0].pcolormesh(chmap.T)

        xs = []; ys = []; minangs = []
        for i in range(np.size(chd)):
            angs = np.real(np.arccos(sm[i] * sed + sthm[i] * sthed * np.cos(pm[i] - ped)))
            xs.append(sm[i])
            ys.append(pm[i])
            minangs.append(np.min(angs))

        axs[1].pcolormesh(ed)

        axs[2].scatter(ys, xs, s = 0.1, c = minangs)

        for ax in axs:
            ax.set_xticks([])
            ax.set_yticks([])
        plt.tight_layout()
        #plt.savefig('./plots/coronal_holes.png')
        plt.show()

    file_flag = 1 # flag to change file string to know which .nc outputs are using wsa 4.5 edge detection

    # - output distance in degrees:
    return np.rad2deg(chd)

def get_observation_times(src_folder):
    '''
    #This gets the available GONG input data times
    '''
    times = []
    for fname in os.listdir(os.getcwd() + "/" + src_folder):
         if fname.startswith("evo.Earth"):
            times.append(datetime.datetime.strptime(fname[-13:-3], "%Y%m%d%H"))

    return sorted(times)


















