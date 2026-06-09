"""
    PREPARE_UPDATE.PY - prepare update of existing DUMFRIC-FCAST run.

    A Yeates -- Apr 2025
"""
# Other modules:
import os
import sys
import datetime
import ftplib
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter as gauss
from scipy.ndimage import label
from scipy.io import netcdf_file
import astropy.units as u
from sunpy.coordinates.sun import carrington_rotation_number, carrington_rotation_time
# Modules in this directory:
import preptools
import data_gong
import uncurl
import pfss
from Parameters import Parameters
import output_netcdf
sys.path.append('./viz/tools/')
import ayplot

# Path to existing run directory [with trailing slash]:
inputpath = "/Volumes/bmjg46/prs/dumfric/fcast-repl2_u20250409.12/"

# Desired finish time:
t_f = "20250411.12"
dtime_f = datetime.datetime.strptime(t_f, "%Y%m%d.%H")
print("Updated forecast time: t_f = "+t_f)

# Create new output directory:
outputpath = inputpath[:-14]+"_u"+t_f+"/"
print("Create output directory") 
os.system("rm -r " + outputpath)
os.system("mkdir -p " + outputpath)

# Parameters for GONG magnetograms:
gong_smooth = 5e-5
gong_hrly_nsmooth = 7
gong_scaling_factor = 1.2
gong_sig = 3
gong_bpar = 23.5
gong_imbalance = 0.5

# Number of smoothing steps in defining local vector potential for emerging regions:
# (two separate steps)
vecpot_nsmooth = 5
vecpot_nsmooth_twist = 4

# How to choose twist parameter ("uniform" or "alpha_obs"):
twist_type = "uniform"

# Twist coefficient [suggested default 0.1 for "uniform" or 10 for "alpha_obs"]:
tau0 = 0.0

# Determine whether regions are to be replaced or overwritten:
replace = (inputpath[-19:-15] == 'repl')


# GET RESTART SNAPSHOT FOR RUN
# ============================
# - time just before earliest tentative emergence in original run.
rtfiles = []
for file in os.listdir(inputpath):
    if file.startswith('rt'):
        rtfiles.append(file)
rtfiles.sort()
t_restart = rtfiles[0][2:13]


# RECREATE SFT MODEL USED FOR PROCESSING EMERGENCE
# ================================================
# Read previous prepared.nml file
with open(inputpath+"prepared.nml") as oldfile:
    flines = oldfile.readlines()
flines[6] = "   end_time = " + t_f + ",\n"
nr = int(flines[7][8:-2])
ns = int(flines[8][8:-2])
nph = int(flines[9][8:-2])
rss = float(flines[10][9:-2])
cadence_em = int(flines[13][16:-2])
# - initial condition:
fh = netcdf_file(inputpath+"b"+t_restart+".nc", "r", mmap=False)
br0 = fh.variables["br"][:,:,0].T
fh.close()
# - surface flows:
fh = netcdf_file(inputpath+"surface_flows.nc", "r", mmap=False)
eta = fh.variables["etasurf"][:]
vs = fh.variables["vs"][:]
vp = fh.variables["vp"][:]
fh.close()
# - initialise SFT model:
sft = preptools.SFT(ns, nph)
dr, rg = preptools.make_grid_r(nr, rss)
sft.prep_sft(vs, vp, eta)
sft.setbr(br0[1:-1,1:-1])
sft.setdt()
t_sft = datetime.datetime.strptime(t_restart, "%Y%m%d.%H")
# - set parameters:
parameters = Parameters(outputpath, '', '', nr, ns, nph, rss, 0,  0, 0, 0, 0, 0, 0, '', '.TRUE.', '.TRUE.', cadence_em)


# COPY ACROSS PARAMETER FILES THAT ARE UNCHANGED, AND RESTART SNAPSHOT
# ====================================================================
os.system("cp "+inputpath+"friction.nc "+outputpath)
os.system("cp "+inputpath+"outflow.nc "+outputpath)
os.system("cp "+inputpath+"surface_flows.nc "+outputpath)
os.system("cp "+inputpath+"b"+t_restart+".nc "+outputpath)


# GENERATE NEW NAMELIST FILES
# ===========================
# Update end time:
with open(outputpath+"prepared.nml", mode="w") as newfile:
    for fline in flines:
        newfile.write(fline)
# Read previos user.nml file and update directory:
with open("user.nml") as oldfile:
    flines = oldfile.readlines()
print(flines)
flines[1] = 'datadir = "' + outputpath[:-1] + '",\n'
with open("user.nml", mode="w") as newfile:
    for fline in flines:
        newfile.write(fline)    


# FIND NEW EMERGING REGIONS
# =========================
bregs = np.empty(shape=(0, ns, nph))
pregs = []
cr_regs = []
type_regs = []

# [1] Extract definitive regions from hourly maps
# -----------------------------------------------
print("Extracting definitive emerging regions from hourly synoptic maps...")
os.system('mkdir '+outputpath+"im-prepare/")
with open(inputpath+"last_def_map.txt") as fid:
    t_last_def = fid.readline()[:-1]
    cr_first_incomplete = int(fid.readline())
dtime_last_def = datetime.datetime.strptime(t_last_def, "%Y%m%d.%H")
if dtime_last_def < dtime_f:
    # - New definitive time period to be added:
    while (dtime_last_def < dtime_f):
        # Step forward 60 degrees in map:
        dtime_last_def += datetime.timedelta(days=4.545)
        print(dtime_last_def)
        cr1 = int(carrington_rotation_number(dtime_last_def)) 
        bregs1, brmap1, pregs1, lonmin, dtime_actual = data_gong.get_fluxregions_hrly(dtime_last_def, ns, nph, sig=gong_sig, bpar=gong_bpar, unbalance=gong_imbalance, nsmooth=gong_hrly_nsmooth, scaling=gong_scaling_factor, plots=True, outpath=outputpath, lontype='definitive', strip=[120, 180], remove_prv=(cr1==cr_first_incomplete))
        if pregs1 != []:
            bregs = np.concatenate((bregs, bregs1), axis=0)
            pregs += pregs1
            cr_regs += [cr1 for p in pregs1]
            type_regs += ['rdh' for p in pregs1]
    # Save time of last definitive map to file:
    with open(outputpath+"last_def_map.txt", mode="w") as fid:
        fid.write(dtime_actual.strftime("%Y%m%d.%H\n"))
        fid.write("%4.4i" % cr_first_incomplete)
else:
    # - No new definitive time period to be added. Get required details from previous definitive map.
    _, _, mapedge = data_gong.readhrlymap_gong(dtime_last_def, ns, nph, nsmooth=gong_hrly_nsmooth)
    lon2 = ((mapedge + 60) % 360) + 360
    lonmin = lon2 + 120
    cr1 = int(carrington_rotation_number(dtime_last_def))

# [2] Extract tentative regions from most recent hourly map
# ---------------------------------------------------------
print("Extracting tentative emerging regions from hourly synoptic map...")
bregs2, brmap2, pregs2, _, _ = data_gong.get_fluxregions_hrly(dtime_f, ns, nph, sig=gong_sig, bpar=gong_bpar, unbalance=gong_imbalance, nsmooth=gong_hrly_nsmooth, scaling=gong_scaling_factor, plots=True, outpath=outputpath, lontype='tentative', incl_lonmin=lonmin, remove_prv=(cr1==cr_first_incomplete))
if pregs2 != []:
    bregs = np.concatenate((bregs, bregs2), axis=0)
    pregs += pregs2
    cr1 = int(carrington_rotation_number(dtime_f)) 
    cr_regs += [cr1 for p in pregs2]
    type_regs += ['rt' for p in pregs2]

# Get emergence time for each region and sort into time order
# -----------------------------------------------------------
if pregs != []:
    t_ems = []
    for k, preg in enumerate(pregs):
        # - get emergence time from carrington longitude:
        pcen = pregs[k] % (2*np.pi)
        if ((type_regs[k] == "rt") | (type_regs[k] == "rdh")) & (pregs[k] > 4*np.pi):
            t_em = carrington_rotation_time(cr_regs[k]-1, longitude=pcen*u.rad).datetime
        else:
            t_em = carrington_rotation_time(cr_regs[k], longitude=pcen*u.rad).datetime
        # - for tentative regions, ensure they have finished emerging by the end of the simulation:
        if (type_regs[k] == "rt") & (t_em > dtime_f):
            t_em = dtime_f
        t_ems.append(t_em)
    idx = sorted(range(len(t_ems)), key=lambda index: t_ems[index])
    t_ems = [t_ems[i] for i in idx]
    pregs = [pregs[i] for i in idx]
    cr_regs = [cr_regs[i] for i in idx]
    type_regs = [type_regs[i] for i in idx]
    bregs = bregs[idx,:,:]

# Create emerging region files
# ----------------------------
if pregs != []:
    for k, preg in enumerate(pregs):
        pcen = pregs[k] % (2*np.pi)
        # - time when emergence starts [force completion time to be noon then subtract emtime]:
        t_em = t_ems[k]
        t_emstart = t_em.replace(hour=12)
        t_emstart -= datetime.timedelta(seconds=3600*cadence_em)
        # [advance SFT model to time t_emstart, if necessary]
        if (t_emstart > t_sft):
            nsteps = (t_emstart - t_sft).days * sft.ndt
            sft.evolve(nsteps)
            t_sft = t_emstart
        # [advance SFT model to t_em0 without this region, if necessary]
        t_emfin = t_em.replace(hour=12)
        if (t_emfin > t_sft):
            nsteps = (t_emfin - t_sft).days * sft.ndt
            sft.evolve(nsteps)
            t_sft = t_emfin
        # - find region number (in case of multiple regions present at same time):
        ks = 1
        while (os.path.exists(parameters.outputpath+type_regs[k][:2]+t_em.strftime('%Y%m%d.%H')+('_%4.4i.nc' % ks))):
            ks += 1
        br1 = bregs[k,:,:]
        # - get mask of insertion region by smoothing |br| [shift to 180 deg first]:
        msk = np.abs(br1.copy())
        msk /= np.max(msk)
        i = np.argmin(np.abs(sft.pc - pcen)) - parameters.nph//2
        if replace:
            # [subtract SFT-evolved br without this region]
            br1 = br1 - sft.br[1:-1,1:-1]
        msk = np.roll(msk, -i, axis=1)
        br1 = np.roll(br1, -i, axis=1)
        # - set nsmooth (amount to smooth msk) to smallest value that avoids any disconnected
        #   components or holes in the region
        gsmooth = 0
        ncpts = 2
        nholes = 1
        while ((ncpts > 1) | (nholes > 0)):
            msk2 = gauss(msk.copy(), gsmooth)
            _, ncpts = label(msk2 > 1e-3)
            _, nholes = label(msk2 <= 1e-3)
            nholes -= 1
            if ((ncpts > 1) | (nholes > 0)):
                gsmooth += 0.1
        msk = msk2
        # - get region:
        msk1 = np.zeros((parameters.ns+2, parameters.nph+2))
        msk1[1:-1,1:-1] = (msk2 > 1e-3).astype('int')
        br1[msk1[1:-1,1:-1] == 0] = 0
        br1 = preptools.correct_flux_multiplicative(br1)
        # - compute vector potential:
        als0, alp0 = uncurl.e_local(parameters.ns, parameters.nph, -br1, msk1, inductive=True)
        # - smooth divergence with dA/dt = grad(div(A)):
        diva = np.zeros((parameters.ns+1, parameters.nph+1))
        dt_smooth = 0.1
        als, alp = als0.copy(), alp0.copy()
        for m in range(vecpot_nsmooth):
            # divergence at grid points:
            diva[1:-1,1:-1] = als[1:,1:-1] - als[:-1,1:-1] + alp[1:-1,1:] - alp[1:-1,:-1]
            als += dt_smooth*(diva[1:,:] - diva[:-1,:])
            alp += dt_smooth*(diva[:,1:] - diva[:,:-1])
        # - add twist to region without changing br:
        s2, ph2 = np.meshgrid(sft.sc, sft.pc, indexing='ij')
        if twist_type == "uniform":
            scen = np.mean(s2[msk1[1:-1,1:-1] == 1])
            tau = 0
            if scen > 0:
                tau = -tau0
            if scen < 0:
                tau = tau0
        # [first make smoother version of br1 for computing twist function]
        br1sm = br1.copy()
        for m in range(vecpot_nsmooth_twist):
            br1sm[1:-1,1:-1] += dt_smooth*(br1sm[1:-1,2:] + br1sm[1:-1,:-2] + br1sm[:-2,1:-1] + br1sm[:-2,1:-1] - 4*br1sm[1:-1,1:-1])
        als, alp, phi = preptools.addtwist(rg, sft.sc, sft.pc, als, alp, br1sm, msk1[1:-1,1:-1], tau, dsmooth=2, plot=False)
        # - recompute br as check:
        br2 = (als[:,1:] - als[:,:-1])/sft.ds/sft.dph - (alp[1:,:] - alp[:-1,:])/sft.ds/sft.dph
        # - find region number (in case multiple regions present on same day):
        ks = 1
        while (os.path.exists(parameters.outputpath+type_regs[k][:2]+t_emstart.strftime('%Y%m%d.%H')+('_%4.4i.nc' % ks))):
            ks += 1
        # - plot:
        axes = ayplot.make_axes([2,2,2,2,2,2], aspect=[2,2,2], height=5, ncols=3)
        ayplot.mesh(axes[0], np.rad2deg(sft.pg), sft.sg, br1, fmax=50, title = (t_emstart.strftime('%Y%m%d.%H')+('\_%4.4i'%ks)+r', $B_r$'), ylabel='Sine Latitude')
        axes[0].contour(np.rad2deg(sft.pc), sft.sc, msk1[1:-1,1:-1], [0.5], colors='k', linewidths=0.75)
        alsmax = np.max(np.abs(als))
        ayplot.mesh(axes[1], np.rad2deg(sft.pc), sft.sg, als0[:,1:-1], ylabel='Sine Latitude', fmax=alsmax, title=r'$A_s$')
        alpmax = np.max(np.abs(alp))
        ayplot.mesh(axes[4], np.rad2deg(sft.pg), sft.sc, alp0[1:-1,:], ylabel='Sine Latitude', fmax=alpmax, title=r'$A_\phi$')
        ayplot.mesh(axes[2], np.rad2deg(sft.pc), sft.sg, als[:,1:-1], ylabel='Sine Latitude', fmax=alsmax, title=r'$A_s$ smoothed+twisted')
        ayplot.mesh(axes[5], np.rad2deg(sft.pg), sft.sc, alp[1:-1,:], ylabel='Sine Latitude', fmax=alpmax, title=r'$A_\phi$ smoothed+twisted')
        lonmin = np.min(np.rad2deg(ph2[msk1[1:-1,1:-1] == 1])) - 10
        lonmax = np.max(np.rad2deg(ph2[msk1[1:-1,1:-1] == 1])) + 10
        smin = np.min(s2[msk1[1:-1,1:-1] == 1]) - 0.1
        smax = np.max(s2[msk1[1:-1,1:-1] == 1]) + 0.1
        for j in range(6):
            axes[j].set_xlim(lonmin, lonmax)
            axes[j].set_ylim(smin, smax)
        # - rotate back to original longitude:
        # (for als0 and als, have to be careful with boundary values)
        als0[:,:-1] = np.roll(als0[:,:-1], i, axis=1)
        als0[:,-1] = als0[:,0]
        alp0 = np.roll(alp0, i, axis=1)
        als[:,:-1] = np.roll(als[:,:-1], i, axis=1)
        als[:,-1] = als[:,0]
        alp = np.roll(alp, i, axis=1)
        br1 = np.roll(br1, i, axis=1)
        br2 = np.roll(br2, i, axis=1)
        msk1 = np.roll(msk1, i, axis=1)
        # - if all is well output region:
        if (np.sum(np.isnan(als0)) + np.sum(np.isnan(alp0))) > 0:
            continue
        if ( np.max(np.abs(br2-br1)) < 0.1*np.max(np.abs(br1)) ) & ((np.sum(np.isnan(als0)) + np.sum(np.isnan(alp0))) == 0):
            output_netcdf.region(parameters.outputpath+type_regs[k][:2]+t_emstart.strftime('%Y%m%d.%H')+('_%4.4i.nc' % ks), als.T, alp.T)
            # [add this region into SFT model]
            sft.addbr(br2)
        ayplot.mesh(axes[3], np.rad2deg(sft.pg), sft.sg, sft.br[1:-1,1:-1], fmax=50)
        axes[3].contour(np.rad2deg(sft.pc), sft.sc, msk1[1:-1,1:-1], [0.5], colors='k', linewidths=0.75)
        # - save plot:
        plt.tight_layout()
        plt.savefig(parameters.outputpath+'im-prepare/'+type_regs[k][:2]+'_a'+t_emstart.strftime('%Y%m%d.%H')+('_%4.4i.png' % ks))
        plt.close()
        del(axes)
