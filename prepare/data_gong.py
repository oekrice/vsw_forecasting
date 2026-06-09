"""
    Python tools for reading NSO GONG synoptic maps and mapping to dumfric grid.
    
    Includes routines to read in (i) Carrington rotation maps and (ii) hourly maps.
    
    For details of data: https://gong.nso.edu/data/magmap/
    https://magmap.nso.edu/PRODUCTS
    
    A Yeates - Oct 2022
"""
import numpy as np
import matplotlib.pyplot as plt
import ftplib
from astropy.io import fits
import astropy.units as units
from sunpy.coordinates.sun import carrington_rotation_time, L0
from sunpy.io._fits import get_header
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage import gaussian_filter
from scipy.ndimage.measurements import label
import preptools
import datetime
import os, sys
try:
    sys.path.append('./viz/tools/')
    import ayplot
except:
    pass

#--------------------------------------------------------------------------------
def readhrlymap_gong(t, ns, nph, nsmooth=0):
    """
        Download the most recent hourly-updated GONG synoptic map for time t, map to the DUMFRIC grid, and correct flux balance.
        
        ARGUMENTS:
            t is a datetime object
            ns and nph define the required grid (e.g. 180 and 360)
            nsmooth [optional] controls the number of steps of gaussian smoothing (default 0 is no smoothing)
            
        Note that the spherical harmonic smoothing is not applied here because in the intended use there will be no neighbouring map later in time. Instead, a simple gaussian smoothing is applied.
    """
    
    # (1) READ IN DATA FROM FTP
    # -------------------------

    from astropy.config import set_temp_cache, get_cache_dir
    with set_temp_cache(path=os.getcwd() + "/", delete=False):
        ftp = ftplib.FTP('gong2.nso.edu')
        ftp.login()

        file_t, t_actual = get_gonghrly_filename(t, ftp)

        try:
            brm = (fits.open('ftp://gong2.nso.edu/'+file_t))[0].data
            print('FOUND MOST RECENT GONG HOURLY-UPDATED MAP FOR TIME '+t.strftime('%Y-%m-%d %H:%M')+' AT '+t_actual.strftime('%Y-%m-%d %H:%M'))
            hdr = get_header(fits.open('ftp://gong2.nso.edu/'+file_t))
            mapedge = hdr[0]['MAPEDGE']
        except Exception as error_message:
            print('! FAILED TO LOAD RECENT MAP FOR TIME '+t.strftime('%Y-%m-%d %H:%M'))
            print('Error message:', error_message)
            brm = np.zeros((ns, nph))
            mapedge = 0

    ftp.close()

    nsm = np.size(brm, axis=0)
    npm = np.size(brm, axis=1)
    dsm = 2.0/nsm
    dpm = 2*np.pi/npm
    scm = np.linspace(-1 + 0.5*dsm, 1 - 0.5*dsm, nsm)
    pcm = np.linspace(0.5*dpm, 2*np.pi - 0.5*dpm, npm)
        
    # Remove NaNs:
    brm = np.nan_to_num(brm)

    # (2) APPLY SIMPLE GAUSSIAN SMOOTHING
    # -----------------------------------
    if (nsmooth > 0):
        dt = 0.1
        for k in range(nsmooth):
            brm[1:-1,:] += dt*(np.roll(brm[1:-1,:],-1,axis=1) + np.roll(brm[1:-1,:],1,axis=1) + brm[:-2,:] + brm[2:,:] - 4*brm[1:-1,:])
        
    # (3) INTERPOLATE TO COMPUTATIONAL GRID
    # -------------------------------------
    # Form computational grid arrays:
    ds = 2.0/ns
    dph = 2*np.pi/nph
    sc = np.linspace(-1 + 0.5*ds, 1 - 0.5*ds, ns)
    pc = np.linspace( 0.5*dph, 2*np.pi - 0.5*dph, nph)

    # Note: in the hourly-updated maps, Earth position is always at 60 degrees.
    # Account for this by changing pc.
    # dms = L0(time=t_actual).dms
    # print((dms[0]+dms[1]/60+dms[2]/3600))
    # clon = np.deg2rad((dms[0]+dms[1]/60+dms[2]/3600))
    clon = np.deg2rad(mapedge)
    pc = np.mod(pc - clon, 2.*np.pi)
    
    # Interpolate to the computational grid:
    bri = RectBivariateSpline(pcm, scm, brm.T)
    br = np.zeros((ns, nph))
    for i in range(nph):
        br[:,i] = bri(pc[i], sc).T.flatten()

    del(brm)

    # (4) CORRECT FLUX BALANCE
    # ------------------------
    br = preptools.correct_flux_multiplicative(br)
    
    return br, t_actual, mapedge


#--------------------------------------------------------------------------------
def readcrmap_gong(rot, ns, nph, smooth=0, scaling=1):
    """
        Download the GONG synoptic map for Carrington rotation rot, map to the DUMFRIC grid, and correct flux balance.
        Also read in the neighbouring maps, and put them together for smoothing.
        
        ARGUMENTS:
            rot is the number of the required Carrington rotation (e.g. 2190)
            ns and nph define the required grid (e.g. 180 and 360)
            smooth [optional] controls the strength of smoothing (default 0 is no smoothing)
        
        Sets br=0 if no synoptic map is found (for either the main rotation or the previous/next ones).
    """

    # (1) READ IN DATA FROM FTP AND STITCH TOGETHER 3 ROTATIONS
    # ---------------------------------------------------------
    ftp = ftplib.FTP('gong2.nso.edu')
    ftp.login()
    
    file_rot = get_gongcr_filename(rot, ftp)
    try:
        brm = (fits.open('ftp://gong2.nso.edu/'+file_rot))[0].data
        print('FOUND GONG CR MAP FOR CR%4.4i' % rot)

        file_rot_l = get_gongcr_filename(rot+1, ftp)
        brm_l = (fits.open('ftp://gong2.nso.edu/'+file_rot_l))[0].data
        print('FOUND GONG CR MAP FOR CR%4.4i' % (rot+1))

        file_rot_r = get_gongcr_filename(rot-1, ftp)
        brm_r = (fits.open('ftp://gong2.nso.edu/'+file_rot_r))[0].data
        print('FOUND GONG CR MAP FOR CR%4.4i' % (rot-1))
    except:
        print('! FAILED TO LOAD MAP FOR CR%4.4i AND NEIGHBOURS' % rot)
        brm = np.zeros((ns, nph))

        
    nsm = np.size(brm, axis=0)
    npm = np.size(brm, axis=1)
    dsm = 2.0/nsm
    dpm = 2*np.pi/npm
    scm = np.linspace(-1 + 0.5*dsm, 1 - 0.5*dsm, nsm)
    pcm = np.linspace(0.5*dpm, 2*np.pi - 0.5*dpm, npm)
    
    # Stitch together:
    brm3 = np.concatenate((brm_l, brm, brm_r), axis=1)
    del(brm, brm_l, brm_r)

    # Remove NaNs:
    brm3 = np.nan_to_num(brm3)

    # Coordinates of combined map (pretend it goes only once around Sun in longitude!):
    nsm = np.size(brm3, axis=0)
    npm = np.size(brm3, axis=1)
    dsm = 2.0/nsm
    dpm = 2*np.pi/npm
    scm = np.linspace(-1 + 0.5*dsm, 1 - 0.5*dsm, nsm)
    pcm = np.linspace(0.5*dpm, 2*np.pi - 0.5*dpm, npm)

    # (2) SMOOTH COMBINED MAP WITH SPHERICAL HARMONIC FILTER
    # ------------------------------------------------------
    if (smooth > 0):
        brm3 = preptools.sh_smooth(brm3, smooth)

    # (3) INTERPOLATE CENTRAL MAP TO COMPUTATIONAL GRID
    # -------------------------------------------------
    # Form computational grid arrays:
    ds = 2.0/ns
    dph = 2*np.pi/nph
    sc = np.linspace(-1 + 0.5*ds, 1 - 0.5*ds, ns)
    pc1 = np.linspace( 0.5*dph, 2*np.pi - 0.5*dph, nph)
    pc = pc1/3 + 2*np.pi/3  # coordinate on the stitched grid

    # Interpolate to the computational grid:
    bri = RectBivariateSpline(pcm, scm, brm3.T)
    br = np.zeros((ns, nph))
    for i in range(ns):
        br[i,:] = bri(pc, sc[i]).T.flatten()

    # (4) INTERPOLATE LEFT AND RIGHT MAPS TO COMPUTATIONAL GRID
    # ---------------------------------------------------------
    brl = np.zeros((ns, nph))
    brr = np.zeros((ns, nph))
    for i in range(ns):
        brl[i,:] = bri(pc - 2*np.pi/3, sc[i]).T.flatten()
        brr[i,:] = bri(pc + 2*np.pi/3, sc[i]).T.flatten()

    del(brm3, bri)

    # (5) CORRECT FLUX BALANCE
    # ------------------------
    br = preptools.correct_flux_multiplicative(br)

    # (6) APPLY SCALING
    # -----------------
    br *= scaling
    
    return br, brl, brr

#--------------------------------------------------------------------------------
def get_gonghrly_filename(t, ftp, backdays=28):
    """
    Identify filename for most recent hourly-updated GONG map on FTP server, for given time t. If no map is found in preceding "backdays" days, return ''.
    """
    t0 = t
    # Look back for maps:
    for i in range(backdays):
        # subdir = t0.strftime('QR/bqs/%Y%m/mrbqs%y%m%d/')
        subdir = t0.strftime('QR/zqs/%Y%m/mrzqs%y%m%d/')
        try:
            # List files for that day:
            ftp.cwd(subdir)
            try:
                gongfiles = ftp.nlst()
                ftp.cwd('~/')
                # Look backward through files and find first one that is before given time (if any):
                for f1 in gongfiles[::-1]:
                    t1 = datetime.datetime.strptime(f1[5:16], '%y%m%dt%H%M')
                    if (t1 <= t):
                        return subdir+f1, t1
            except:
                pass
        except:
            pass
        t0 -= datetime.timedelta(days=1)
    
    return '', t

#--------------------------------------------------------------------------------
def get_gongcr_filename(rot, ftp):
    """
    Identify filename for a GONG map on FTP server, otherwise return ''.
    """
    # Identify date corresponding to 180 Carrington longitude (middle of map).
    # [This seems to be how GONG CR maps are labelled.]
    t0 = carrington_rotation_time(rot, longitude=180*units.deg)
    t0.format = 'datetime'
    
    # Get ftp subdirectory [YYYYmm/]:
    # subdir = t0.strftime('QR/mqs/%Y%m/mrmqs%y%m%d/')
    subdir = t0.strftime('QR/nqs/%Y%m/mrnqs%y%m%d/')
    try:
        ftp.cwd(subdir)
    except:
        return ''

    # List files in the directory:
    gongfiles = ftp.nlst()
    ftp.cwd('~/')

    for file in gongfiles:
        if file[17:21] == ('%4.4i' % rot):
            return subdir+file
    return ''

def get_fluxregions(rot, ns, nph, sig=3, bpar=39.8, unbalance=0.5, smooth=0, scaling=1, plots=False, outpath='./'):
    """
        Determine strong flux regions from GONG integrated (full) synoptic map for rotation rot.
        
        Returns arrays on s, phi grid (corrected for flux balance), and longitude centroid.
    """
    
    # Get map and neighbours smoothed and interpolated on dumfric grid.
    # Map itself is corrected for flux balance but neighbours are not.
    br, brl, brr = readcrmap_gong(rot, ns, nph, smooth=smooth, scaling=scaling)
    br3 = np.concatenate((brl, br, brr), axis=1)
    del(brl, brr)
       
    # Coordinates on triple map:
    np3 = np.size(br3, axis=1)
    dp3 = 6*np.pi/np3
    ds = 2./ns
    pc3 = np.linspace(0.5*dp3, 6*np.pi - 0.5*dp3, np3)
    sc3 = np.linspace(-1 + 0.5*ds, 1 - 0.5*ds, ns) 
    
    if (plots):
        axes = ayplot.make_axes([2,2,2,2], ncols=1, aspect=[2,4,1], height=10)
        ayplot.mesh(axes[0], np.rad2deg(pc3), sc3, br3.T, fmax=50)
        axes[0].plot([360, 360], [-1,1], 'k--', linewidth=0.5)
        axes[0].plot([720, 720], [-1,1], 'k--', linewidth=0.5)
        axes[0].text(150, 1.1, '[CR%i]' % (rot+1))
        axes[0].text(520, 1.1, 'CR%i' % rot)
        axes[0].text(850, 1.1, '[CR%i]' % (rot-1))

    # Identify strong flux regions by smoothing absolute flux then contouring.
    m3 = np.abs(br3)
    m3 = gaussian_filter(m3, sig)
    labels, nregs = label(m3 > bpar)
    if (plots):
        ayplot.mesh(axes[1], np.rad2deg(pc3), sc3, m3.T, fmax=bpar)
        axes[1].plot([360, 360], [-1,1], 'k--', linewidth=0.5)
        axes[1].plot([720, 720], [-1,1], 'k--', linewidth=0.5) 
        axes[1].text(20, 0.8, 'SIG = %i,  BPAR = %g' % (sig, bpar))
        ayplot.mesh(axes[2], np.rad2deg(pc3), sc3, labels.T, cmap='nipy_spectral', fmin=0)
        axes[2].plot([360, 360], [-1,1], 'w--', linewidth=0.5)
        axes[2].plot([720, 720], [-1,1], 'w--', linewidth=0.5) 
        
    # Remove regions with too large a flux imbalance, otherwise correct for flux balance:
    ratio = np.zeros(nregs)
    for i in range(1,nregs+1):
        netflux = np.sum(br3[labels==i])
        absflux = np.sum(np.abs(br3[labels==i]))
        if (np.abs(netflux)/absflux > unbalance):
            labels[labels==i] = 0
        else:
            br3[labels==i] = preptools.correct_flux_multiplicative(br3[labels==i])

    # Remove regions with longitude centroid (area not flux) outside original map, 
    # or regions centred too near the pole (due to artefacts - e.g. CR2017 in SOLIS).
    sc, pc = np.meshgrid(sc3, pc3, indexing='ij')
    bregs = []
    pregs = []
    for i in range(1,nregs+1):
        if (np.sum(labels==i) > 0):
            pcen = np.mean(pc[labels==i])
            scen = np.mean(sc[labels==i])
            if ((pcen < 2*np.pi) | (pcen > 4*np.pi) | (abs(scen) > 0.9)):
                labels[labels==i] = 0
            else:
                # Create list of corrected Br arrays for good regions (ensure to include wrap-arounds).
                m1 = br3.copy()
                m1[labels!=i] = 0
                bregs.append(m1[:,:nph] + m1[:, nph:2*nph] + m1[:,2*nph:])
                pregs.append(pcen)
    
    bregs = np.array(bregs)
    if (plots):
        if (np.size(bregs, axis=0) == 0):
            bregs = np.zeros((1, ns, nph))
        ayplot.mesh(axes[3], np.rad2deg(pc3[nph:2*nph]), sc3, np.sum(bregs,axis=0).T, fmax=50, xlim=[0, 1080])
        axes[3].plot([360, 360], [-1,1], 'k--', linewidth=0.5)
        axes[3].plot([720, 720], [-1,1], 'k--', linewidth=0.5) 
        axes[3].text(20, 0.8, 'UNBALANCE = %g' % unbalance)
        plt.savefig(outpath+'im-prepare/regions%4.4i.png' % rot, bbox_inches='tight')
        plt.close()
        del(axes)
        
    return bregs, br, pregs

def get_fluxregions_hrly(t, ns, nph, sig=3, bpar=39.8, unbalance=0.5, nsmooth=0, scaling=1, plots=False, outpath='./', lontype='all', def_strip=[120, 180], remove_prv=True):
    """
        Determine strong flux regions from GONG hourly synoptic map.
        
        Returns arrays on s, phi grid (corrected for flux balance), and longitude centroid.
    """
    
    # Get map and neighbours smoothed and interpolated on dumfric grid.
    br, t_actual, mapedge = readhrlymap_gong(t, ns, nph, nsmooth=nsmooth)
    # - concatenate three copies into triple map (to deal with boundaries):
    br3 = np.concatenate((br, br, br), axis=1)
       
    # Coordinates on triple map:
    np3 = np.size(br3, axis=1)
    dp3 = 6*np.pi/np3
    ds = 2./ns
    pc3 = np.linspace(0.5*dp3, 6*np.pi - 0.5*dp3, np3)
    sc3 = np.linspace(-1 + 0.5*ds, 1 - 0.5*ds, ns) 

    lon1 = ((mapedge) % 360) + 360
    lon2 = ((mapedge + 60) % 360) + 360
    lon3 = ((mapedge + 120) % 360) + 360

    if (plots):
        axes = ayplot.make_axes([2,2,2,2], ncols=1, aspect=[2,4,1], height=10)
        ayplot.mesh(axes[0], np.rad2deg(pc3), sc3, br3.T, fmax=50)
        axes[0].plot([360, 360], [-1,1], 'k--', linewidth=0.5)
        axes[0].plot([720, 720], [-1,1], 'k--', linewidth=0.5)
        axes[0].text(150+360, 1.1, t_actual.strftime("%Y%m%d.%H"))
        axes[0].plot([lon1, lon1], [-1, 1], 'k', linestyle=(0,(5,10)), linewidth=0.75)
        axes[0].plot([lon2, lon2], [-1, 1], 'k--', linewidth=0.75)
        lon3 = ((mapedge + 120) % 360) + 360
        axes[0].plot([lon3, lon3], [-1, 1], 'k', linestyle=(0,(5,10)), linewidth=0.75)

    # Identify strong flux regions by smoothing absolute flux then contouring.
    m3 = np.abs(br3)
    m3 = gaussian_filter(m3, sig)
    labels, nregs = label(m3 > bpar)
    if (plots):
        ayplot.mesh(axes[1], np.rad2deg(pc3), sc3, m3.T, fmax=bpar)
        axes[1].plot([360, 360], [-1,1], 'k--', linewidth=0.5)
        axes[1].plot([720, 720], [-1,1], 'k--', linewidth=0.5)
        axes[1].text(20, 0.8, 'SIG = %i,  BPAR = %g' % (sig, bpar))
        axes[1].plot([lon1, lon1], [-1, 1], 'k', linestyle=(0,(5,10)), linewidth=0.75)
        axes[1].plot([lon2, lon2], [-1, 1], 'k--', linewidth=0.75)
        axes[1].plot([lon3, lon3], [-1, 1],  'k', linestyle=(0,(5,10)),linewidth=0.75)
        ayplot.mesh(axes[2], np.rad2deg(pc3), sc3, labels.T, cmap='nipy_spectral', fmin=0)
        axes[2].plot([360, 360], [-1,1], 'w--', linewidth=0.5)
        axes[2].plot([720, 720], [-1,1], 'w--', linewidth=0.5)
        axes[2].plot([lon1, lon1], [-1, 1],  'w', linestyle=(0,(5,10)), linewidth=0.75)
        axes[2].plot([lon2, lon2], [-1, 1], 'w--', linewidth=0.75)
        axes[2].plot([lon3, lon3], [-1, 1],  'w', linestyle=(0,(5,10)), linewidth=0.75)
                
    # Remove regions with too large a flux imbalance, otherwise correct for flux balance:
    ratio = np.zeros(nregs)
    for i in range(1,nregs+1):
        netflux = np.sum(br3[labels==i])
        absflux = np.sum(np.abs(br3[labels==i]))
        if (np.abs(netflux)/absflux > unbalance):
            labels[labels==i] = 0
        else:
            br3[labels==i] = preptools.correct_flux_multiplicative(br3[labels==i])

    # Remove regions with longitude centroid (area not flux) outside original map, 
    # or regions centred too near the pole (due to artefacts - e.g. CR2017 in SOLIS).
    sc, pc = np.meshgrid(sc3, pc3, indexing='ij')
    bregs = []
    pregs = []
    for i in range(1,nregs+1):
        if (np.sum(labels==i) > 0):
            pcen = np.mean(pc[labels==i])
            scen = np.mean(sc[labels==i])
            # Filter out regions depending on lontype argument:
            if lontype == "all":
                if ((pcen < 2*np.pi) | (pcen > 4*np.pi) | (abs(scen) > 0.9)):
                    labels[labels==i] = 0
                else:
                    # Create list of corrected Br arrays for good regions (ensure to include wrap-arounds).
                    m1 = br3.copy()
                    m1[labels!=i] = 0
                    bregs.append(m1[:,:nph] + m1[:, nph:2*nph] + m1[:,2*nph:])
                    pregs.append(pcen)
            elif lontype == "definitive":
                p0 = np.deg2rad((lon2 + def_strip[0]))
                p1 = np.deg2rad((lon2 + def_strip[1]))
                if remove_prv:
                    # - don't include regions with centres in previous carrington rotation:
                    p1 = min(p1, 4*np.pi)
                if ((pcen < p0) | (pcen > p1) | (abs(scen) > 0.9)):
                    labels[labels==i] = 0
                else:
                    # Create list of corrected Br arrays for good regions (ensure to include wrap-arounds).
                    m1 = br3.copy()
                    m1[labels!=i] = 0
                    bregs.append(m1[:,:nph] + m1[:, nph:2*nph] + m1[:,2*nph:])
                    pregs.append(pcen)
            elif lontype == "tentative":
                # include everything in updated region, up to CM + incl_lonmin:
                p0 = np.deg2rad(lon2 - 60)
                p1 = np.deg2rad(lon2 + 120)
                if remove_prv:
                    # - don't include regions with centres in previous carrington rotation:
                    p1 = min(p1, 4*np.pi)
                if ((pcen < p0) | (pcen > p1) | (abs(scen) > 0.9)):
                    labels[labels==i] = 0
                else:
                    # Create list of corrected Br arrays for good regions (ensure to include wrap-arounds).
                    m1 = br3.copy()
                    m1[labels!=i] = 0
                    bregs.append(m1[:,:nph] + m1[:, nph:2*nph] + m1[:,2*nph:])
                    pregs.append(pcen)
    bregs = np.array(bregs)

    if (plots):
        if (np.size(bregs, axis=0) == 0):
            bregs = np.zeros((1, ns, nph))
        ayplot.mesh(axes[3], np.rad2deg(pc3[nph:2*nph]), sc3, np.sum(bregs,axis=0).T, fmax=50, xlim=[0, 1080])
        axes[3].plot([360, 360], [-1,1], 'k--', linewidth=0.5)
        axes[3].plot([720, 720], [-1,1], 'k--', linewidth=0.5) 
        axes[3].text(20, 0.8, 'UNBALANCE = %g' % unbalance)
        axes[3].plot([lon2, lon2], [-1, 1], 'k--', linewidth=0.75)
        axes[3].plot([lon3, lon3], [-1, 1],  'k', linestyle=(0,(5,10)),linewidth=0.75)
        if lontype == "definitive":
            lonstrip0 = ((lon2 + def_strip[0]) % 360) + 360
            lonstrip1 = ((lon2 + def_strip[1]) % 360) + 360
            axes[3].plot([lonstrip0, lonstrip0], [-1, 1],  color='b', linestyle=(0,(5,10)),linewidth=0.75)
            axes[3].plot([lonstrip1, lonstrip1], [-1, 1],  color='b', linestyle=(0,(5,10)),linewidth=0.75)
            axes[3].plot([lon1, lon1], [-1, 1], 'k', linestyle=(0,(5,10)), linewidth=0.75)
        if lontype == "tentative":
            lonstrip0 = ((lon2 + 120) % 360) + 360
            axes[3].plot([lonstrip0, lonstrip0], [-1, 1],  color='r', linestyle=(0,(5,10)),linewidth=0.75)
            axes[3].plot([lon1, lon1], [-1, 1],  color='r', linestyle=(0,(5,10)),linewidth=0.75)

        if lontype == "definitive":
            plt.savefig(outpath+t_actual.strftime("im-prepare/hrly_def_regions_%Y%m%d.%H.png"), bbox_inches='tight')
        elif lontype == "tentative":
            plt.savefig(outpath+t_actual.strftime("im-prepare/hrly_ten_regions_%Y%m%d.%H.png"), bbox_inches='tight')
        plt.close()
        del(axes)

    return bregs, br, pregs
