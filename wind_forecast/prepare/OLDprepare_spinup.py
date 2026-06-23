"""
    PREPARE_SPINUP.PY - prepare new DUMFRIC-FCAST run starting from synoptic map, with given spinup time.

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

# Parent directory for output:
root =  "/Volumes/bmjg46/prs/dumfric/"

# Run identifier:
run_id = "2"

# Whether to insert new regions by overwriting or overlaying:
replace = False
if replace:
    outputpath = root+"fcast-repl"+run_id
else:
    outputpath = root+"fcast-over"+run_id

# Desired finish time:
t_f = "20250401.12"
dtime_f = datetime.datetime.strptime(t_f, "%Y%m%d.%H")
print("Simulation forecast time: t_f = "+t_f)

outputpath += "_s"+t_f+"/"

# Length of spinup in days [needs to be long enough that a full carrington map is available]:
dt_spinup = 60
dtime_0 = dtime_f - datetime.timedelta(days=dt_spinup)
t_0 = dtime_0.strftime("%Y%m%d.%H")
print("Simulation start time:    t_0 = "+t_0)

# Grid sizes in rho, s, phi:
nr = 60
ns = 180
nph = 360 # must be even (for polar boundary condition)

# Outer boundary radius [in Rsun]:
rss = 2.5

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


# CREATE PARAMETER FILES
# ======================
cflfact=0.1 
snapshot_cadence=168 
diagnostic_cadence=10 
nu0=0.36e-5 
ohmic_ratio=0.0
ohmic_enhanced=0.0
hyper_ratio=1.2e-7
snapshot_variables="B"
driver_v='.TRUE.'
driver_em='.TRUE.'
cadence_em=24
starttime = t_0
endtime = t_f

parameters = Parameters(outputpath, starttime, endtime, nr, ns, nph, rss, cflfact,  snapshot_cadence, diagnostic_cadence, nu0, ohmic_ratio, ohmic_enhanced, hyper_ratio, snapshot_variables, driver_v, driver_em, cadence_em)
parameters.write_prepared_file()
parameters.write_user_file('user.nml')

sft = preptools.SFT(parameters.ns, parameters.nph)
dr, rg = preptools.make_grid_r(parameters.nr, parameters.rss)
#-----------------------------------------------------------
# Output surface flows:
RSUN = 6.96e10
SDAY = 86400.0
# (1) Supergranular diffusivity (constant):
eta = 350*1e10/RSUN**2
# (2) Meridional flow velocity (vs on p-ribs):
v0 = 0.015*1e5/RSUN
p0 = 2.33
Du = v0*(1+p0)**(0.5*(p0+1))/p0**(0.5*p0)
vs = Du*sft.sg*(np.sqrt(1 - sft.sg**2))**p0
# (3) Differential rotation velocity (vp on s-ribs):
omA = np.deg2rad(0.18)/SDAY
omB = -np.deg2rad(2.396)/SDAY
omC = -np.deg2rad(1.787)/SDAY
om = omA + omB*sft.sc**2 + omC*sft.sc**4
vp = om*np.sqrt(1-sft.sc**2)
output_netcdf.surfaceFlows(parameters.outputpath+'surface_flows.nc', eta, vs, vp)
#-----------------------------------------------------------
# Output friction coefficient:
# [nu = nu0*r^2*sin^2(th)]:
# - the strength parameter nu0 [in s^-1] is set in user.nml.
rg2, sg2 = np.meshgrid(rg, sft.sg, indexing='ij')
nug = np.exp(2*rg2)*(1 - sg2**2)
nug[0,:] = 0.0  # set to zero at photosphere
output_netcdf.friction(parameters.outputpath+'friction.nc', nug)
#-----------------------------------------------------------
# Output radial outflow speed:
vout = 100*1e5/RSUN
p = 11.5
vr = vout*(np.exp(rg)/parameters.rss)**p
output_netcdf.outflow(parameters.outputpath+'outflow.nc', vr)


# INITIAL CONDITION
# =================
print("Generating initial condition...")
# Use "integrated" (full carrington) synoptic map for the corresponding carrington rotation:
cr_0 = int(carrington_rotation_number(dtime_0))
print('cr_0 = %4.4i' % cr_0)
try:
    br0 = data_gong.readcrmap_gong(cr_0, ns, nph, smooth=gong_smooth, scaling=gong_scaling_factor)[0]    
except:
    print("ERROR: No integral synoptic map available for chosen start date. Extend dt_spinup.")
    sys.exit()
pfss.pfss(br0, parameters.nr, parameters.ns, parameters.nph, parameters.rss, filename=parameters.outputpath+'a_pfss'+parameters.starttime+'.nc')


# EMERGING REGIONS
# ================
cr_f = int(carrington_rotation_number(dtime_f))
# [1] Extract regions from any complete integrated synoptic maps AFTER cr_0
# -------------------------------------------------------------------------
print("Extracting emerging regions from integrated synoptic maps...")
os.system('mkdir '+outputpath+"im-prepare/")
bregs = np.empty(shape=(0, ns, nph))
pregs = []
cr_regs = []
type_regs = []
for cr in range(cr_0+1, cr_f+1):
    print("- CR%4.4i" % cr)
    try:
        br0 = data_gong.readcrmap_gong(cr, parameters.ns, parameters.nph, smooth=gong_smooth, scaling=gong_scaling_factor)[0] 
    except:
        print("[Integrated map not available]")
        break
    bregs1, brmap1, pregs1 = data_gong.get_fluxregions(cr, parameters.ns, parameters.nph, sig=gong_sig, bpar=gong_bpar, unbalance=gong_imbalance, smooth=gong_smooth, scaling=gong_scaling_factor, plots=True, outpath=outputpath)
    if pregs1 != []:
        bregs = np.concatenate((bregs, bregs1), axis=0)
        pregs += pregs1
        cr_regs += [cr for p in pregs1]
        type_regs += ['rd' for p in pregs1]
cr_first_incomplete = cr

# # FUDGE FOR TESTING [create worst case where two carrington maps missing]
# cr -= 1

# [2] Extract definitive regions from hourly maps
# -----------------------------------------------
# i.e. in a strip a suitable buffer zone from most recent update.
# - need to keep track of whether the previous rotation was covered by hourly maps
#  (so fully included already in step [1]) or a finalised carrington map (so not fully included)
print("Extracting definitive emerging regions from hourly synoptic maps...")
dtime = carrington_rotation_time(cr_first_incomplete, 360*u.deg)
while (dtime < dtime_f):
    # Step forward 60 degrees in map:
    dtime += datetime.timedelta(days=4.545)
    cr1 = int(carrington_rotation_number(dtime)) 
    bregs1, brmap1, pregs1, lonmin = data_gong.get_fluxregions_hrly(dtime, parameters.ns, parameters.nph, sig=gong_sig, bpar=gong_bpar, unbalance=gong_imbalance, nsmooth=gong_hrly_nsmooth, scaling=gong_scaling_factor, plots=True, outpath=outputpath, lontype='definitive', strip=[120, 180], remove_prv=(cr1==cr_first_incomplete))
    if pregs1 != []:
        bregs = np.concatenate((bregs, bregs1), axis=0)
        pregs += pregs1
        cr_regs += [cr1 for p in pregs1]
        type_regs += ['rdh' for p in pregs1]
# Save time of last definitive map, and last complete rotation to file:
with open(outputpath+"last_def_map.txt", mode="w") as fid:
    fid.write(dtime.strftime("%Y%m%d.%H\n"))
    fid.write("%4.4i" % cr)

# [3] Extract tentative regions from most recent hourly map
# ---------------------------------------------------------
print("Extracting tentative emerging regions from hourly synoptic map...")
bregs2, brmap2, pregs2, _ = data_gong.get_fluxregions_hrly(dtime_f, parameters.ns, parameters.nph, sig=gong_sig, bpar=gong_bpar, unbalance=gong_imbalance, nsmooth=gong_hrly_nsmooth, scaling=gong_scaling_factor, plots=True, outpath=outputpath, lontype='tentative', incl_lonmin=lonmin, remove_prv=(cr1==cr_first_incomplete))
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
    # [Initialise SFT model]
    sft.prep_sft(vs, vp, eta)
    sft.setbr(br0)
    sft.setdt()
    t_sft = dtime_0
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