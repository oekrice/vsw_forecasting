"""
    Various subroutines used by other prepare scripts.
    
    ARY 19/12/17
"""
import numpy as np
import scipy.linalg as la
from scipy.io import FortranFile


def make_grid_r(nr, rss):
    """
        
    """
    dr = np.log(rss) / nr
    rg = np.linspace(0, np.log(rss), nr + 1)

    return dr, rg


def correct_flux_multiplicative(f):
    """
        Correct the flux balance in the map f (assumes that cells have equal area).
    """

    # Compute positive and negative fluxes:
    ipos = f > 0
    ineg = f < 0
    fluxp = np.abs(np.sum(f[ipos]))
    fluxn = np.abs(np.sum(f[ineg]))

    # Rescale both polarities to mean:
    fluxmn = 0.5 * (fluxn + fluxp)
    f1 = f.copy()
    f1[ineg] *= fluxmn / fluxn
    f1[ipos] *= fluxmn / fluxp

    return f1


def sh_smooth(f, smooth, cutoff=10):
    """
    Smooth array f [on dumfric grid, cell centres, no ghost cells] with spherical harmonic filter exp( -smooth*l*(l+1) ) in spectral space.

    This implementation uses discrete eigenfunctions (in latitude) instead of Plm.

    Parameters:
    smooth -- coefficient of filter exp(-smooth * lam) [set cutoff=0 to include all eigenvalues. This is quicker for small matrices.]
    cutoff -- largest value of smooth*lam to include [so 10 means ignore blm multiplied by exp(-10)]
    """

    nsm = np.size(f, axis=0)
    npm = np.size(f, axis=1)
    dsm = 2.0/nsm
    dpm = 2*np.pi/npm
    scm = np.linspace(-1 + 0.5*dsm, 1 - 0.5*dsm, nsm)
    sgm = np.linspace(-1, 1, nsm+1) 

    # Prepare tridiagonal matrix:
    Fp = sgm * 0  # Lp/Ls on p-ribs
    Fp[1:-1] = np.sqrt(1 - sgm[1:-1] ** 2) / (np.arcsin(scm[1:]) - np.arcsin(scm[:-1])) * dpm
    Vg = Fp / dsm / dpm
    Fs = ((np.arcsin(sgm[1:]) - np.arcsin(sgm[:-1])) / np.sqrt(1 - scm ** 2) / dpm)  # Ls/Lp on s-ribs
    Uc = Fs / dsm / dpm
    # - create off-diagonal part of the matrix:
    off_diag = -Vg[1:nsm]
    # - terms required for m-dependent part of matrix:
    mu = np.fft.fftfreq(npm)
    mu = 4 * np.sin(np.pi * mu) ** 2
    diag1 = Vg[:nsm] + Vg[1:nsm+1]

    # FFT in phi of photospheric distribution at each latitude:
    fhat = np.fft.rfft(f, axis=1)

    # Loop over azimuthal modes (positive m):
    nm = npm//2 + 1
    blm = np.zeros((nsm), dtype="complex")
    fhat1 = np.zeros((nsm, nm), dtype="complex")
    for m in range(nm):
        if (m%100) == 0:
            print(m, nm)
        # - set diagonal terms of matrix:
        diag = diag1 + Uc[:nsm] * mu[m]
        # - compute eigenvectors Q_{lm} and eigenvalues lam_{lm}:
        #   (note that matrix is symmetric tridiag so use special solver)
        if cutoff > 0:
            # - ignore contributions with eigenvalues too large to contribute after smoothing:
            lamax = cutoff/smooth
            lam, Q = la.eigh_tridiagonal(diag, off_diag, select="v", select_range=(0,lamax))
            nsm1 = len(lam) # [full length would be nsm]
        else:
            lam, Q = la.eigh_tridiagonal(diag, off_diag)
            nsm1 = nsm
        # - find coefficients of eigenfunction expansion:
        for l in range(nsm1): 
            blm[l] = np.dot(Q[:,l], fhat[:,m])
            # - apply filter [the eigenvalues should be a numerical approx of lam = l*(l+1)]:
            blm[l] *= np.exp(-smooth*lam[l])
        # - invert the latitudinal transform:
        fhat1[:,m] = np.dot(blm[:nsm1], Q.T)

    # Invert the FFT in longitude:
    f_out = np.real(np.fft.irfft(fhat1, axis=1))
    
    return f_out

# def plgndr(m, x, lmax):
#     """
#         Evaluate associated Legendre polynomials P_lm(x) for given (positive)
#         m, from l=0,lmax, with spherical harmonic normalization included.
#         Only elements l=m:lmax are non-zero.
        
#         Similar to scipy.special.lpmv except that function only works for 
#         small l due to overflow, because it doesn't include the normalizationp.
#     """

#     nx = np.size(x)
#     plm = np.zeros((nx, lmax + 1))
#     pmm = 1
#     if m > 0:
#         somx2 = (1 - x) * (1 + x)
#         fact = 1.0
#         for i in range(1, m + 1):
#             pmm *= somx2 * fact / (fact + 1)
#             fact += 2

#     pmm = np.sqrt((m + 0.5) * pmm)
#     pmm *= (-1) ** m
#     plm[:, m] = pmm
#     if m < lmax:
#         pmmp1 = x * np.sqrt(2 * m + 3) * pmm
#         plm[:, m + 1] = pmmp1
#         if m < lmax - 1:
#             for l in range(m + 2, lmax + 1):
#                 fact1 = np.sqrt(
#                     ((l - 1.0) ** 2 - m ** 2) / (4.0 * (l - 1.0) ** 2 - 1.0)
#                 )
#                 fact = np.sqrt((4.0 * l ** 2 - 1.0) / (l ** 2 - m ** 2))
#                 pll = (x * pmmp1 - pmm * fact1) * fact
#                 pmm = pmmp1
#                 pmmp1 = pll
#                 plm[:, l] = pll
#     return plm


def bmr2d(sc, pc, lat0, lon0, sep0, tilt0, B0, xithresh=9):
    """
        Return magnetic field for 2d tilted magnetic bipole with peak strength B0.
        **Note that all input angles should be in radians.
        
        In an untilted frame where the bipole is at the equator it has the form
            Br = -B0*(lon/sep0)*exp[ 0.5*( 1 - xi ]
            where xi = (lon**2 + 2*lat**2)/sep0**2
        where the strength B1 is chosen so that the maximum of Br is B0.
        
        Note: cutoff based on threshold of xi. 
    """

    ns = np.size(sc)
    nph = np.size(pc)
    ds = sc[1] - sc[0]
    dp = pc[1] - pc[0]
    sc2, pc2 = np.meshgrid(sc, pc, indexing="ij")
    sinc2 = np.sqrt(1 - sc2 ** 2)

    # Cartesian coordinates:
    x = np.cos(pc2) * sinc2
    y = np.sin(pc2) * sinc2
    z = sc2

    # Rotate to frame where BMR is on equator and untilted:
    xb = (
        x * np.cos(lat0) * np.cos(lon0)
        + y * np.cos(lat0) * np.sin(lon0)
        + z * np.sin(lat0)
    )
    yb = (
        x
        * (-np.cos(tilt0) * np.sin(lon0) + np.sin(tilt0) * np.sin(lat0) * np.cos(lon0))
        + y
        * (np.cos(tilt0) * np.cos(lon0) + np.sin(tilt0) * np.sin(lat0) * np.sin(lon0))
        - z * np.sin(tilt0) * np.cos(lat0)
    )
    zb = (
        x
        * (-np.sin(tilt0) * np.sin(lon0) - np.cos(tilt0) * np.sin(lat0) * np.cos(lon0))
        + y
        * (np.sin(tilt0) * np.cos(lon0) - np.cos(tilt0) * np.sin(lat0) * np.sin(lon0))
        + z * np.cos(tilt0) * np.cos(lat0)
    )
    zb[zb > 1] = 1
    zb[zb < -1] = -1

    # Magnetic field of BMR in this frame:
    thb = np.arccos(zb)
    phb = np.arctan2(yb, xb)
    xi = (phb ** 2 + 2 * (0.5 * np.pi - thb) ** 2) / sep0 ** 2
    brb = -B0 * phb / sep0 * np.exp(-xi)

    # Cutoff at threshold:
    brb[xi > xithresh] = 0
    msk = (xi <= xithresh).astype("int")

    # Remove any flux imbalance:
    brb = correct_flux_multiplicative(brb)

    return brb, msk


def bmr2d_a(sg, pg, sc, pc, lat0, lon0, sep0, tilt0, B0, xithresh=9):
    """
        Return analytical (2d) vector potential for tilted magnetic bipole with peak strength B0 (including edge lengths).
        **Note that all input angles should be in radians.
        
        Also returns binary "mask" array for sweeping, with threshold xithresh in terms of "radial" xi coordinate.
        
        See manual for details of the functional form.
    """

    ns = np.size(sc)
    nph = np.size(pc)

    dp = 2 * np.pi / nph

    # - rotation matrix to bmr coordinates (where bmr is on
    #   equator and untilted):
    M = np.zeros((3, 3))
    M[0, :] = np.array(
        [np.cos(lat0) * np.cos(lon0), np.cos(lat0) * np.sin(lon0), np.sin(lat0)]
    )
    M[1, :] = np.array(
        [
            -np.cos(tilt0) * np.sin(lon0) + np.sin(tilt0) * np.sin(lat0) * np.cos(lon0),
            np.cos(tilt0) * np.cos(lon0) + np.sin(tilt0) * np.sin(lat0) * np.sin(lon0),
            -np.sin(tilt0) * np.cos(lat0),
        ]
    )
    M[2, :] = np.array(
        [
            -np.sin(tilt0) * np.sin(lon0) - np.cos(tilt0) * np.sin(lat0) * np.cos(lon0),
            np.sin(tilt0) * np.cos(lon0) - np.cos(tilt0) * np.sin(lat0) * np.sin(lon0),
            np.cos(tilt0) * np.cos(lat0),
        ]
    )
    # - inverse:
    Minv = np.linalg.inv(M)

    # - constant amplitude:
    A0 = 0.5 * sep0 * B0

    # [1] As
    sc2, pg2 = np.meshgrid(sc, pg, indexing="ij")
    lc2 = np.sqrt(1 - sc2 ** 2)
    # - cartesian coordinates:
    x = np.cos(pg2) * lc2
    y = np.sin(pg2) * lc2
    z = sc2
    # - rotate to frame where bmr is on equator and untilted:
    xb = M[0, 0] * x + M[0, 1] * y + M[0, 2] * z
    yb = M[1, 0] * x + M[1, 1] * y + M[1, 2] * z
    zb = M[2, 0] * x + M[2, 1] * y + M[2, 2] * z
    dy = np.arcsin(zb)  # latitude in bmr frame
    dx = np.arctan2(yb, xb)  # longitude in bmr frame
    ay1 = -A0 * np.cos(dy) * np.exp(-(dx ** 2 + 2 * dy ** 2) / sep0 ** 2)
    # - convert to Cartesian components:
    axb = np.cos(dx) * np.sin(dy) * ay1
    ayb = np.sin(dx) * np.sin(dy) * ay1
    azb = -np.cos(dy) * ay1
    # - rotate back to global frame:
    ax = Minv[0, 0] * axb + Minv[0, 1] * ayb + Minv[0, 2] * azb
    ay = Minv[1, 0] * axb + Minv[1, 1] * ayb + Minv[1, 2] * azb
    az = Minv[2, 0] * axb + Minv[2, 1] * ayb + Minv[2, 2] * azb
    # - extract as spherical component:
    als = -np.cos(pg2) * sc2 * ax - np.sin(pg2) * sc2 * ay + lc2 * az
    # - multiply by edge lengths:
    for j in range(ns):
        als[j, :] *= np.arcsin(sg[j + 1]) - np.arcsin(sg[j])

    # [2] Ap
    sg2, pc2 = np.meshgrid(sg, pc, indexing="ij")
    lg2 = np.sqrt(1 - sg2 ** 2)
    # - cartesian coordinates:
    x = np.cos(pc2) * lg2
    y = np.sin(pc2) * lg2
    z = sg2
    # - rotate to frame where bmr is on equator and untilted:
    xb = M[0, 0] * x + M[0, 1] * y + M[0, 2] * z
    yb = M[1, 0] * x + M[1, 1] * y + M[1, 2] * z
    zb = M[2, 0] * x + M[2, 1] * y + M[2, 2] * z
    dy = np.arcsin(zb)  # latitude in bmr frame
    dx = np.arctan2(yb, xb)  # longitude in bmr frame
    ay1 = -A0 * np.cos(dy) * np.exp(-(dx ** 2 + 2 * dy ** 2) / sep0 ** 2)
    # - convert to Cartesian components:
    axb = np.cos(dx) * np.sin(dy) * ay1
    ayb = np.sin(dx) * np.sin(dy) * ay1
    azb = -np.cos(dy) * ay1
    # - rotate back to global frame:
    ax = Minv[0, 0] * axb + Minv[0, 1] * ayb + Minv[0, 2] * azb
    ay = Minv[1, 0] * axb + Minv[1, 1] * ayb + Minv[1, 2] * azb
    az = Minv[2, 0] * axb + Minv[2, 1] * ayb + Minv[2, 2] * azb
    # - extract as spherical component:
    alp = -np.sin(pc2) * ax + np.cos(pc2) * ay
    # - multiply by edge lengths:
    for j in range(ns + 1):
        alp[j, :] *= dp * np.sqrt(1 - sg[j] ** 2)

    # [3] Mask array for sweeping
    sc2, pc2 = np.meshgrid(sc, pc, indexing="ij")
    lc2 = np.sqrt(1 - sc2 ** 2)
    # - cartesian coordinates:
    x = np.cos(pc2) * lc2
    y = np.sin(pc2) * lc2
    z = sc2
    # - rotate to frame where bmr is on equator and untilted:
    xb = M[0, 0] * x + M[0, 1] * y + M[0, 2] * z
    yb = M[1, 0] * x + M[1, 1] * y + M[1, 2] * z
    zb = M[2, 0] * x + M[2, 1] * y + M[2, 2] * z
    dy = np.arcsin(zb)  # latitude in bmr frame
    dx = np.arctan2(yb, xb)  # longitude in bmr frame
    xi = (dx ** 2 + 2 * dy ** 2) / sep0 ** 2

    msk = (xi <= xithresh).astype("int")

    return als, alp, msk


def addtwist(rg, sc, pc, als, alp, brb, msk, tau, dsmooth=2, plot=True):
    """
        Take 2d vector potential, Br, and mask for emerging region, and add horizontal gradient
        term to als and alp, corresponding to twisting with strength tau.
        
        Sign of tau corresponds with sign of helicity.
    """
    ns = np.size(sc)
    nph = np.size(pc)
    ds = sc[1] - sc[0]
    dp = pc[1] - pc[0]

    sc2, pcs = np.meshgrid(sc, pc, indexing="ij")
    Lp = 2 * dp * np.sqrt(1 - sc2 ** 2)  # 2 cells
    Ls = np.arcsin(sc2[2:, 1:-1]) - np.arcsin(sc2[:-2, 1:-1])

    # Gradient of signed Br:
    dbp = np.zeros((ns, nph))
    dbp[1:-1, 1:-1] = 0.5 * (brb[1:-1, 2:] - brb[1:-1, :-2]) / Lp[1:-1, 1:-1]
    dbp[1:-1, 0] = 0.5 * (brb[1:-1, 1] - brb[1:-1, -1]) / Lp[1:-1, 0]
    dbp[1:-1, -1] = 0.5 * (brb[1:-1, 0] - brb[1:-1, -2]) / Lp[1:-1, -1]
    dbs = np.zeros((ns, nph))
    dbs[1:-1, 1:-1] = 0.5 * (brb[2:, 1:-1] + brb[:-2, 1:-1]) / Ls
    g = np.sqrt(dbp ** 2 + dbs ** 2)

    # Gradient of unsigned Br:
    br11 = np.abs(brb)
    dbp = np.zeros((ns, nph))
    dbp[1:-1, 1:-1] = 0.5 * (br11[1:-1, 2:] - br11[1:-1, :-2]) / Lp[1:-1, 1:-1]
    dbp[1:-1, 0] = 0.5 * (br11[1:-1, 1] - br11[1:-1, -1]) / Lp[1:-1, 0]
    dbp[1:-1, -1] = 0.5 * (br11[1:-1, 0] - br11[1:-1, -2]) / Lp[1:-1, -1]
    dbs = np.zeros((ns, nph))
    dbs[1:-1, 1:-1] = 0.5 * (br11[2:, 1:-1] + br11[:-2, 1:-1]) / Ls
    ug = np.sqrt(dbp ** 2 + dbs ** 2)

    # Identify pil:
    pil = (g - ug) ** 2

    # Define potential by smoothing:
    nsmoo = int(np.deg2rad(dsmooth) / 0.2 / dp)
    phic = np.zeros((ns + 2, nph + 2))
    phic[1:-1, 1:-1] = pil
    # - set to zero outside of emergence region [smoothing will then extend it]:
    phic[1:-1, 1:-1] *= msk
    dtdiff = 0.2
    for step in range(nsmoo):
        # - apply periodic boundary conditions in phi:
        phic[:, 0] = phic[:, -2]
        phic[:, -1] = phic[:, 1]
        # - apply diffusion:
        phic[1:-1, 1:-1] += dtdiff * (
            -4 * phic[1:-1, 1:-1]
            + phic[:-2, 1:-1]
            + phic[2:, 1:-1]
            + phic[1:-1, :-2]
            + phic[1:-1, 2:]
        )
    pils = phic[1:-1, 1:-1].copy()
    # - multiply by Br:
    phic[1:-1, 1:-1] *= brb
    # - average to grid points:
    phic[:, 0] = phic[:, -2]
    phic[:, -1] = phic[:, 1]
    phi = 0.25 * (phic[1:, 1:] + phic[:-1, 1:] + phic[1:, :-1] + phic[:-1, :-1])
    # - normalize:
    dr1 = np.exp(rg[1]) - np.exp(rg[0])
    b1max = np.max(
        np.sqrt(
            0.25 * (phi[1:-1, 2:] - phi[1:-1, :-2]) ** 2
            + 0.25 * (phi[2:, 1:-1] - phi[:-2, 1:-1]) ** 2
        )
        / dr1
    )
    phi *= tau * np.max(np.abs(brb)) / b1max

    if plot:
        import matplotlib.pyplot as plt
        import matplotlib as mpl

        mpl.use("Agg")

        s2, ph2 = np.meshgrid(sc, pc, indexing="ij")
        lonmin = np.min(np.rad2deg(ph2[msk == 1])) - 10
        lonmax = np.max(np.rad2deg(ph2[msk == 1])) + 10
        smin = np.min(s2[msk == 1]) - 0.1
        smax = np.max(s2[msk == 1]) + 0.1

        plt.rc("text", usetex=True)
        plt.rc("font", family="serif")
        plt.figure(figsize=(7, 9))
        ax = plt.subplot(411)
        pm = ax.pcolormesh(np.rad2deg(pc), sc, brb, cmap="bwr")
        bmax = np.max(np.abs(brb))
        pm.set_clim(vmin=-bmax, vmax=bmax)
        plt.colorbar(pm)
        ax.set_xlim(0, 360)
        ax.set_ylim(-1, 1)
        ax.set_title(r"$B_r$")
        ax.set_ylabel("s")
        ax.set_xlim(lonmin, lonmax)
        ax.set_ylim(smin, smax)

        ax = plt.subplot(412)
        pm = ax.pcolormesh(np.rad2deg(pc), sc, pil, cmap="bwr")
        bmax = 500  # np.max(np.abs(pil))
        pm.set_clim(vmin=-bmax, vmax=bmax)
        plt.colorbar(pm)
        ax.set_xlim(0, 360)
        ax.set_ylim(-1, 1)
        ax.set_title(r"$f_{\rm pil}$")
        ax.set_ylabel(r"$s$")
        ax.set_xlim(lonmin, lonmax)
        ax.set_ylim(smin, smax)

        ax = plt.subplot(413)
        pm = ax.pcolormesh(np.rad2deg(pc), sc, pils, cmap="bwr")
        bmax = np.max(np.abs(pils))
        pm.set_clim(vmin=-bmax, vmax=bmax)
        plt.colorbar(pm)
        ax.set_xlim(0, 360)
        ax.set_ylim(-1, 1)
        ax.set_title(r"$\langle f_{\rm pil} \rangle$")
        ax.set_ylabel(r"$s$")
        ax.set_xlim(lonmin, lonmax)
        ax.set_ylim(smin, smax)

        ax = plt.subplot(414)
        sg = np.linspace(-1, 1, ns + 1)
        pg = np.linspace(0, 2 * np.pi, nph + 1)
        pm = ax.pcolormesh(np.rad2deg(pg), sg, phi, cmap="bwr")
        bmax = np.max(np.abs(phi))
        pm.set_clim(vmin=-bmax, vmax=bmax)
        plt.colorbar(pm)
        ax.set_xlim(0, 360)
        ax.set_ylim(-1, 1)
        ax.set_title(r"$\Phi$")
        ax.set_ylabel("s")
        ax.set_xlabel(r"$\phi$ [degrees]")
        ax.set_xlim(lonmin, lonmax)
        ax.set_ylim(smin, smax)

        plt.savefig("demo_twist.png", bbox_inches="tight")
        plt.close()

    # Compute vector potential by grad(Phi) - edge lengths are included:
    als += phi[1:, :] - phi[:-1, :]
    alp += phi[:, 1:] - phi[:, :-1]
    return als, alp, phi


def lapeigf(msk, ndiff=100):
    """
    Given a mask function (1 inside domain, 0 outside), compute first Laplace eigenfunction on
    this domainp.
    """

    nx = np.size(msk, axis=0)  # not including ghost layer
    ny = np.size(msk, axis=1)
    imsk = np.nonzero(msk == 1)
    nmsk = np.size(imsk, axis=1)

    # ID of each grid cell (-1 if not in mask):
    idd = np.zeros((nx, ny)) - 1
    idd[imsk] = np.linspace(0, nmsk - 1, nmsk)
    # Add surrounding zero cells:
    idds = np.zeros((nx + 2, ny + 2)) - 1
    idds[1:-1, 1:-1] = idd
    # Periodic boundary in x:
    idds[0, :] = idds[-2, :]
    idds[-1, :] = idds[1, :]

    # Construct matrix A (with step-size h=1):
    # -- number of rows and columns is nmsk
    # -- these are psi in row order starting from bottom left, on interior of mask region only.
    id_dn = idds[1 : nx + 1, 0:ny]
    id_up = idds[1 : nx + 1, 2 : ny + 2]
    id_lft = idds[0:nx, 1 : ny + 1]
    id_rt = idds[2 : nx + 2, 1 : ny + 1]
    A = np.zeros((nmsk, nmsk))
    for i in range(nx):
        for j in range(ny):
            idd1 = int(idd[i, j])
            if idd1 > -1:
                id_up1 = int(id_up[i, j])
                id_dn1 = int(id_dn[i, j])
                id_lft1 = int(id_lft[i, j])
                id_rt1 = int(id_rt[i, j])
                A[idd1, idd1] = 4
                if id_up1 > -1:
                    A[idd1, id_up1] = -1
                else:
                    A[idd1, idd1] += 1
                if id_dn1 > -1:
                    A[idd1, id_dn1] = -1
                else:
                    A[idd1, idd1] += 1
                if id_lft1 > -1:
                    A[idd1, id_lft1] = -1
                else:
                    A[idd1, idd1] += 1
                if id_rt1 > -1:
                    A[idd1, id_rt1] = -1
                else:
                    A[idd1, idd1] += 1

    # Find the first eigenvalue and corresponding eigenfunction:
    # - I think the eigenvalues should be real since the matrix A is symmetric
    lam, w = np.linalg.eig(A)
    i0 = np.argmin(lam)
    psi = np.zeros((nx, ny))
    psi[imsk] = w[:, i0]

    # Normalize to be negative:
    psi = -np.abs(psi)

    # Smooth it slightly by diffusion (periodic boundaries):
    psig = np.zeros((nx + 2, ny + 2))
    psig[1:-1, 1:-1] = psi
    dtdiff = 0.2
    for step in range(ndiff):
        # - apply periodic boundary conditions in x:
        psig[0, :] = psig[-2, :]
        psig[-1, :] = psig[1, :]
        # - apply diffusion:
        psig[1:-1, 1:-1] += dtdiff * (
            -4 * psig[1:-1, 1:-1]
            + psig[:-2, 1:-1]
            + psig[2:, 1:-1]
            + psig[1:-1, :-2]
            + psig[1:-1, 2:]
        )

    del psi

    psig[0, :] = psig[-2, :]
    psig[-1, :] = psig[1, :]
    return psig


class SFT:
    """
        Surface flux transport model.
    """

    def __init__(self, ns, nph):
        self.ns = ns
        self.nph = nph
        self.ds = 2.0 / ns
        self.dph = 2 * np.pi / nph
        self.sg = np.linspace(-1, 1, ns + 1)
        self.sc = np.linspace(-1 + 0.5 * self.ds, 1 - 0.5 * self.ds, ns)
        self.pg = np.linspace(0, 2 * np.pi, nph + 1)
        self.pc = np.linspace(0.5 * self.dph, 2 * np.pi - 0.5 * self.dph, nph)

    def prep_sft(self, vs1, vp1, eta):
        """
            Prepare for surface flux transport.
            Input should be 1D arrays of vs and vp, and a constant eta.
        """
        self.eta = eta

        self.es = np.zeros((self.ns, self.nph + 1))
        self.ep = np.zeros((self.ns + 1, self.nph))
        self.Ns, _ = np.meshgrid(
            np.sqrt(1 - self.sc ** 2) * self.dph, self.pg, indexing="ij"
        )
        self.Np, _ = np.meshgrid(
            np.arcsin(self.sc[1:]) - np.arcsin(self.sc[:-1]), self.pc, indexing="ij"
        )
        self.Ls, _ = np.meshgrid(
            np.arcsin(self.sg[1:]) - np.arcsin(self.sg[:-1]), self.pg, indexing="ij"
        )
        self.Lp, _ = np.meshgrid(
            np.sqrt(1 - self.sg ** 2) * self.dph, self.pc, indexing="ij"
        )
        self.vs, _ = np.meshgrid(vs1, self.pc, indexing="ij")
        self.vp, _ = np.meshgrid(vp1, self.pg, indexing="ij")

    def setdt(self):
        """
            Set the timestep by the CFL conditionp.
        """
        hpmin = np.min(np.abs(self.Ns))
        hsmin = np.min(np.abs(self.Ls))
        dt_eta = min([hpmin ** 2 / self.eta, hsmin ** 2 / self.eta])
        dt_mf = np.min(np.abs(self.Np[:, :] / self.vs[1:-1, :]))
        dt_om = np.min(np.abs(self.Ns[:, :] / self.vp[:, :]))
        dt = 0.2 * min([dt_eta, dt_mf, dt_om])
        # - modify to fit exactly in one day:
        self.ndt = int(86400.0 / dt)  # -number of steps in a day
        self.dt = 86400.0 / self.ndt
        print("Timestep for SFT precomputation is %g secs" % self.dt)

    def setbr(self, br0):
        """
            Set br.
        """
        self.br = np.zeros((self.ns + 2, self.nph + 2))
        self.br[1:-1, 1:-1] = br0.copy()

    def addbr(self, br0):
        """
            Add a given array to br.
        """
        self.br[1:-1, 1:-1] += br0
        self.br[:, 0] = self.br[:, -2]
        self.br[:, -1] = self.br[:, 1]

    def evolve(self, nsteps):
        """
            Evolve br by flux transport evolution for nsteps time steps.
            Input map br0 should have no ghost cells.
        """

        for step in range(nsteps):
            # (i) apply periodic boundary conditions in phi:
            self.br[:, 0] = self.br[:, -2]
            self.br[:, -1] = self.br[:, 1]
            # (ii) average br to ribs:
            brs = 0.5 * (self.br[1:-1, :-1] + self.br[1:-1, 1:])
            brp = 0.5 * (self.br[:-1, 1:-1] + self.br[1:, 1:-1])
            # (iii) emf on s-edges from differential rotation:
            self.es = self.vp * brs
            # (iv) emf on p-edges from meridional flow:
            self.ep = -self.vs * brp
            # (v) emf on ribs from supergranular diffusion:
            self.es -= self.eta * (self.br[1:-1, 1:] - self.br[1:-1, :-1]) / self.Ns
            self.ep[1:-1, :] += (
                self.eta * (self.br[2:-1, 1:-1] - self.br[1:-2, 1:-1]) / self.Np
            )
            # (vi) update br by Faraday+Stokes:
            self.br[1:-1, 1:-1] -= (
                self.dt
                * (
                    self.es[:, 1:] * self.Ls[:, 1:]
                    - self.es[:, :-1] * self.Ls[:, :-1]
                    - self.ep[1:, :] * self.Lp[1:, :]
                    + self.ep[:-1, :] * self.Lp[:-1, :]
                )
                / self.ds
                / self.dph
            )

    def flux(self, br1):
        """
            Compute flux of array br1 on SFT grid.
        """

        return np.sum(np.abs(br1)) * self.ds * self.dph * (6.96e10) ** 2
