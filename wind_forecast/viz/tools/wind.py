"""
    Routines for solar wind output from DuMFric.

    ary -- 2025 May
"""
import os
this_directory = os.getcwd() + "/"

import sys
import numpy as np
from scipy.io import netcdf_file, FortranFile
from scipy.ndimage import sobel
import datetime
from astropy.time import Time
import astropy.units as u
from sunpy.coordinates.sun import B0
import matplotlib.pyplot as plt

import copy


def windmap(snap, r_hb, path="./", codepath="./viz/fortran/"):
    """
    Using precompiled fortran, computes the Schatten extension, coronal hole map, expansion factors etc.
    Doesn't compute the wind speed, as that has many quick-to-change parameters

    This saves out to a 'windmap' netcdf file
    """
    #Check file exists
    if not os.path.exists(path+snap):
        raise Exception('Outflow/PFSS output not found at location ', path+snap)
    # Compute coronal hole map
    # ------------------------
    os.system(codepath + "bin/tracer " + path + " " + snap + " chmap")
    print('Coronal hole map calculated.')
    # Compute Schatten extension and combined field line mapping
    # ----------------------------------------------------------
    os.system(codepath + "bin/tracer " + path + " " + snap + " wind %g" % r_hb)
    print('Schatten extension calculated.')

    return


def windbnd(snap, r_hb, path="./", codepath="./viz/fortran/"):
    """
    Using the outputs from the computed Schatten field, determine the distribution of radial velocities at r_hb
    """
    # Compute coronal hole boundary distances and flux tube expansion factors
    # -----------------------------------------------------------------------
    fid = FortranFile(path + "windmap_" + snap + ".unf", "r")
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
    chd[igood] = chbmap_45(snap, sm[igood], phm[igood], path=path, codepath=codepath)

    # Compute map of flux-tube expansion factors traced down through combined model:
    fs = sm * 0
    fs0 = expansionfactor(snap, 1, 2.5, brm[igood], br1[igood], path=path)
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

    # Empirical model for vr on outer boundary:
    # vr = compute_vr(chd, fs, method="wsa", params=[250.0, 625.0, 0.25, 1, 0.8, 1, 4, 1])
    vr = compute_vr(chd, fs, method="wsa", params = [285, 625+285, 0.22222, 1, 0.8, 2, 2, 3])

    # Save boundary data to file:
    viewdate = datetime.datetime.strptime(snap[-14:-3], "%Y%m%d.%H")
    fid = netcdf_file(path + "windbound_" + snap, "w")
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
    vid = fid.createVariable("vr", "d", ("nth", "nph"))
    vid[:] = vr
    fid.date = viewdate.strftime("%Y-%m-%d %H:00")
    fid.nth = ns
    fid.nph = nph
    fid.close()
    

def chbmap_old(snap, sm, pm, path="./", codepath="./fortran/"):
    """
        Compute map of coronal hole boundary distances given 2d arrays of points (sm, pm) on photosphere.
    """

    # Compute coronal hole map (if necessary):
    if not os.path.exists(path + "chmap_" + snap + ".unf"):
        utils.compile_f90(codepath=codepath)
        os.system(codepath + "bin/tracer " + path + " " + snap + " chmap")
    fid = FortranFile(path + "chmap_" + snap + ".unf", "r")
    sc = fid.read_reals(dtype=np.float64)
    pc = fid.read_reals(dtype=np.float64)
    chmap = fid.read_ints(dtype=np.int32).reshape(
        (np.size(pc, 0) - 1, np.size(sc, 0) - 1)
    )
    chmap = np.swapaxes(chmap, 0, 1)
    fid.close()

    # - detect edges in map:
    sx = sobel(chmap, axis=0, mode="constant")
    sy = sobel(chmap, axis=1, mode="constant")
    ed = (np.hypot(sx, sy) > 0.005).astype("int")
    ed[0, :] = 0
    ed[-1, :] = 0
    ed[:, 0] = 0
    ed[:, -1] = 0
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

    # - output distance in degrees:
    return np.rad2deg(chd)

def chbmap_45(snap, sm, pm, path="./", codepath="./fortran/", make_plot=False):

    # Compute coronal hole map (if necessary):
    if not os.path.exists(os.path.join(path, "chmap_" + snap + ".unf")):
        utils.compile_f90(codepath=codepath)
        os.system(codepath + "/bin/tracer " + path + "/ " + snap + " chmap")

    fid = FortranFile(os.path.join(path, "chmap_" + snap + ".unf"), "r")
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
        plt.savefig('./plots/coronal_holes.png')
        plt.show()

    file_flag = 1 # flag to change file string to know which .nc outputs are using wsa 4.5 edge detection

    # - output distance in degrees:
    return np.rad2deg(chd)

def expansionfactor(snap, rm, r0, brm, br0, path="./"):
    """
        Compute map of flux tube expansion factor given br0 at r=r0 and brm at r=rm (radii) of each field line.
    """

    fs = brm * rm ** 2 / br0 / r0 ** 2

    return fs


def compute_vr(chb, fs, method="marion", params=[200.0, 700.0, 7.0, 1/2.5]):
    """
        Compute map of the solar wind speed v_r given the coronal hole boundary distance (chb, in degrees) and flux tube expansion factor (fs).
        ary -- 2019/09/13
    """

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

        # To account for change in resolution:
        chb1 = chb + 2

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
            + ((vfast - vslow) / (1.0 + fs) ** a) * (b - g * np.exp(-(chb1 / w) ** d)) ** i
        )

    if method == "riley":
        # e.g. params=[0.1051, 0.0101, 333.0, 631.0]
        ep = params[0]
        w = params[1]
        vslow = params[2]
        vfast = params[3]
        vr = vslow + 0.5*(vfast - vslow) * (1.0 + np.tanh((np.deg2rad(chb) - ep) / w))

    return vr


def get_huxt_vsw(windbound_file, r_hb, fcast_length=7, decelerate=True, cone2bcfile='', plot2d=False, spinup_cme_days=0):
    """
    For given windbound file, run HUXt model and generate forecast vsw time series at 1 AU.
    - assumes a static heliospheric field
    - length of forecast is fcast_length days.
    - start inserting CMEs spinup_cme_days days before magnetogram time.

    If a cone2bcfile is provided, CMEs are inserted from it. (Format used for Enlil.)
    
    For testing, set plot2d=True to run HUXt in 2d and plot result instead of returning.
    """

    # Read model output on heliospheric boundary:
    fh = netcdf_file(windbound_file, "r", mmap=False)
    fs = fh.variables["expansionFactor"][:]
    chd = fh.variables["CHBDistance"][:]
    vr_bnd = fh.variables["vr"][:]
    br_bnd = fh.variables["br"][:]
    date = fh.date.decode("utf-8")
    fh.close()
    
    dtime_snap = datetime.datetime.strptime(date, "%Y-%m-%d %H:%M")
    simtime = 27.27*u.day

    # Latitude of Earth at simulation time:
    _, vr_longs, vr_lats, br_map, br_longs, br_lats = Hin.get_PFSS_maps(windbound_file)
    E_lat = np.deg2rad(B0(dtime_snap))
    iE_lat = np.argmin(abs(vr_lats - E_lat))
    
    # Run HUXt:
    t_init = dtime_snap - datetime.timedelta(days=spinup_cme_days)
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
    if cone2bcfile != '':
        cme_list = Hin.ConeFile_to_ConeCME_list(model, cone2bcfile)
        cmes_out = []
        for cme in cme_list:
            cme_replace = copy.deepcopy(cme)
            if cme.v > model.v_max:
                print(cme.v)
                cme_replace.v = model.v_max - 1*(u.km/u.s)
            cmes_out.append(cme_replace)
        cme_list = cmes_out
        model.solve(cme_list, tag='cone_cme_test')
    
        # # Manual CME for testing:
        # cme = H.ConeCME(t_launch=-0.2*u.day, longitude=360*u.deg, latitude = 0*u.deg, width=60*u.deg, v=1000*(u.km/u.s), thickness=5*u.solRad, initial_height=r_hb*u.solRad)
        # model.solve([cme], tag='cone_cme_test')
    else:
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

    del earth_series; del fs; del chd; del vr_bnd; del br_bnd

    return dtime, vsw, bpol
