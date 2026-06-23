import numpy as np
import scipy.linalg as la
from . import output_netcdf


def pfss(
    br0,
    nr,
    ns,
    nph,
    rss,
    r0=1,
    filename="",
    output="a",
    testQ=False,
    outerBC="radial",
    br1=[],
):
    """
        Compute 3D potential field using eigenfunction method in r,s,p coordinates, on the dumfric grid (equally spaced in
        rho=ln(r/rsun), s=cos(theta0), and p=phi).
       
        The output should have zero current to machine precision,
        when computed with the DuMFriC staggered discretization.
       
        Outer boundary condition is controlled by the flag "outerBC":
        outerBC='radial': Bs = Bp = 0 on r=rss
        outerBC='br': Br = br1 on r=rss
        outerBC='inf': semi-infinite solution
       
        Output depends on the flag 'output':
         output='none': as it says
         output='a': ar*Lr, as*Ls, ap*Lp on cell edges.
         output='bc': br, bs, bp on the centres of the cell faces.
         output='bg': br, bs, bp (weighted) averaged to grid points.
        
        Set testQ=True to compare the discrete eigenfunctions Qj_{lm}  to Plm(cos(th)).
        
    Copyright (C) Anthony R. Yeates, Durham University 29/8/17
    ** modified 1/2/19 to add option of br outer boundary condition.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    """

    # Coordinates:
    ds = 2.0 / ns
    dp = 2 * np.pi / nph
    dr = (np.log(rss) - np.log(r0)) / nr
    rg = np.linspace(np.log(r0), np.log(rss), nr + 1)
    rc = np.linspace(np.log(r0) + 0.5 * dr, np.log(rss) - 0.5 * dr, nr)
    sg = np.linspace(-1, 1, ns + 1)
    sc = np.linspace(-1 + 0.5 * ds, 1 - 0.5 * ds, ns)

    k = np.linspace(0, nr, nr + 1)

    Fp = sg * 0  # Lp/Ls on p-ribs
    Fp[1:-1] = np.sqrt(1 - sg[1:-1] ** 2) / (np.arcsin(sc[1:]) - np.arcsin(sc[:-1])) * dp
    Vg = Fp / ds / dp
    Fs = (
        (np.arcsin(sg[1:]) - np.arcsin(sg[:-1])) / np.sqrt(1 - sc ** 2) / dp
    )  # Ls/Lp on s-ribs
    Uc = Fs / ds / dp

    # FFT in phi of photospheric distribution at each latitude:
    brt = np.fft.rfft(br0, axis=1)

    # FFT in phi of outer distribution at each latitude (if required):
    if outerBC == "br":
        brt1 = np.fft.rfft(br1, axis=1)

    # Prepare tridiagonal matrix:
    # - create off-diagonal part of the matrix:
    A = np.zeros((ns, ns))
    for j in range(ns - 1):
        A[j, j + 1] = -Vg[j + 1]
        A[j + 1, j] = A[j, j + 1]
    # - term required for m-dependent part of matrix:
    mu = np.fft.fftfreq(nph)
    mu = 4 * np.sin(np.pi * mu) ** 2
    # - initialise:
    psir = np.zeros((nr + 1, ns), dtype="complex128")
    psi = np.zeros((nr + 1, ns, nph), dtype="complex")
    e1 = np.exp(dr)
    fact = np.sinh(dr) * (e1 - 1)

    if testQ:
        import scipy.special as sp
        import matplotlib.pyplot as plt

        plt.figure()

    # Loop over azimuthal modes (positive m):
    for m in range(nph // 2 + 1):
        # - set diagonal terms of matrix:
        for j in range(ns):
            A[j, j] = Vg[j] + Vg[j + 1] + Uc[j] * mu[m]
        # - compute eigenvectors Q_{lm} and eigenvalues lam_{lm}:
        #   (note that A is symmetric so use special solver)
        lam, Q = la.eigh(A)
        # - solve quadratic:
        Flm = 0.5 * (1 + e1 + lam * fact)
        ffp = Flm + np.sqrt(Flm ** 2 - e1)
        if outerBC == "inf":
            ffm = Flm - np.sqrt(Flm * Flm - e1)
        else:
            ffm = e1 / ffp
        # - compute radial term for each l (for this m):
        for l in range(ns):
            # - sum c_{lm} + d_{lm} [from BC at photosphere]:
            cdlm = np.dot(Q[:, l], brt[:, m]) / lam[l]
            if outerBC == "radial":
                # - ratio c_{lm}/d_{lm} [numerically safer this way up]
                ratio = (1 - ffm[l]) / (ffp[l] - 1)
                ratio *= (ffm[l] / ffp[l]) ** (nr - 1)
                dlm = cdlm / (1.0 + ratio)
                clm = ratio * dlm
            if outerBC == "br":
                # - sum c_{lm}*ffp**nr + d_{lm}*ffm**nr [from BC at photosphere]:
                cdlm1 = np.dot(Q[:, l], brt1[:, m]) / lam[l] * rss ** 2
                # - solve simultaneously to get clm and dlm:
                clm = (cdlm1 - ffm[l] ** nr * cdlm) / (ffp[l] ** nr - ffm[l] ** nr)
                dlm = (cdlm1 - ffp[l] ** nr * cdlm) / (ffm[l] ** nr - ffp[l] ** nr)
            if outerBC == "inf":
                dlm = cdlm
                clm = cdlm * 0
            # Set psir[:,l] = clm*ffp[l]**k + dlm*ffm[l]**k
            with np.errstate(
                over="raise"
            ):  # catch overflow for high res and split product
                for kdiv in range(1, nr + 1):
                    kd = k // kdiv
                    psir[:, l] = clm
                    try:
                        for k1 in range(kdiv - 1):
                            psir[:, l] *= ffp[l] ** kd
                        psir[:, l] *= ffp[l] ** (k - (kdiv - 1) * kd)
                        break
                    except FloatingPointError:
                        continue
            psir[:, l] += dlm * ffm[l] ** k

        # - compute entry for this m in psit = Sum_l c_{lm}Q_{lm}**j
        psi[:, :, m] = np.dot(psir, Q.T)
        if m > 0:
            psi[:, :, nph - m] = np.conj(psi[:, :, m])

        if testQ & (m == 6):
            isrt = np.argsort(lam, axis=0)  # sort eigenvalues
            lam = lam[isrt]
            istat = np.indices((ns, ns))
            Q = Q[istat[0], isrt]
            plt.clf()
            for l in range(5):
                plm = sp.lpmv(m, m + l, sc)
                Ql = Q[:, l] * Q[1, l] / np.abs(Q[1, l])  # normalise and match sign
                plt.plot(sc, Ql / np.max(np.abs(Ql)), "ko")
                plm = plm * plm[1] / np.abs(plm[1])
                plt.plot(sc, plm / np.max(np.abs(plm)), label="l=%i" % l)
                plt.xlabel(r"$\cos(\theta)$")
            plt.title("m = %i" % m)
            plt.legend()
            plt.savefig("Q.png", bbox_inches="tight")
            plt.show()

    del (psir, mu, A)

    # Compute psi by inverse fft:
    psi = np.real(np.fft.ifft(psi, axis=2))

    # Hence compute vector potential [note index order, for netcdf]:
    alr = np.zeros((nph + 1, ns + 1, nr))
    als = np.zeros((nph + 1, ns, nr + 1))
    alp = np.zeros((nph, ns + 1, nr + 1))

    for j in range(nr + 1):
        for i in range(nph + 1):
            als[i, :, j] = Fs * (psi[j, :, ((i - 1) % nph)] - psi[j, :, ((i) % nph)])
        for i in range(nph):
            alp[i, 1:-1, j] = Fp[1:-1] * (psi[j, 1:, i] - psi[j, :-1, i])

    # Output to netcdf file:
    r = np.exp(rg)
    th = np.arccos(sg)
    ph = np.linspace(0, 2 * np.pi, nph + 1)

    if output == "none":
        return alr, als, alp

    if output == "a":
        output_netcdf.a(filename, r, th, ph, alr, als, alp)

    if (output == "bc") | (output == "bg"):
        rc = np.linspace(-0.5 * dr, np.log(rss) + 0.5 * dr, nr + 2)
        rrc = np.exp(rc)
        thc = np.zeros(ns + 2) - 1
        thc[1:-1] = np.arccos(sc)
        pc = np.linspace(-0.5 * dp, 2 * np.pi + 0.5 * dp, nph + 2)

        # Required face normals:
        dnph = np.zeros((ns + 2, 2))
        dns = np.zeros((ns + 1, 2))
        dnr = np.zeros(ns + 2)
        for k in range(2):
            for j in range(1, ns + 1):
                dnph[j, k] = rrc[k] * np.sqrt(1 - sc[j - 1] ** 2) * dp
            dnph[0, k] = dnph[1, k]
            dnph[-1, k] = dnph[-2, k]
            for j in range(1, ns):
                dns[j, k] = rrc[k] * (np.arcsin(sc[j]) - np.arcsin(sc[j - 1]))
            dns[0, k] = dns[1, k]
            dns[-1, k] = dns[-2, k]
        for j in range(ns + 2):
            dnr[j] = rrc[0] * (np.exp(dr) - 1)
        dnr[0] = -dnr[0]
        dnr[-1] = -dnr[-1]

        # Required area factors:
        Sbr = np.zeros((ns + 2, nr + 1))
        for k in range(nr + 1):
            Sbr[1:-1, k] = np.exp(2 * rg[k]) * ds * dp
            Sbr[0, k] = Sbr[1, k]
            Sbr[-1, k] = Sbr[-2, k]
        Sbs = np.zeros((ns + 1, nr + 2))
        for k in range(nr + 2):
            for j in range(1, ns):
                Sbs[j, k] = (
                    0.5
                    * np.exp(2 * rc[k] - dr)
                    * dp
                    * (np.exp(2 * dr) - 1)
                    * np.sqrt(1 - sg[j] ** 2)
                )
            Sbs[0, k] = Sbs[1, k]
            Sbs[-1, k] = Sbs[-2, k]
        Sbp = np.zeros((ns + 2, nr + 2))
        for k in range(nr + 2):
            for j in range(1, ns + 1):
                Sbp[j, k] = (
                    0.5
                    * np.exp(2 * rc[k] - dr)
                    * (np.exp(2 * dr) - 1)
                    * (np.arcsin(sg[j]) - np.arcsin(sg[j - 1]))
                )
            Sbp[0, k] = Sbp[1, k]
            Sbp[-1, k] = Sbp[-2, k]

        # Compute br*Sbr, bs*Sbs, bp*Sbp at cell centres by Stokes theorem:
        br = np.zeros((nph + 2, ns + 2, nr + 1))
        bs = np.zeros((nph + 2, ns + 1, nr + 2))
        bp = np.zeros((nph + 1, ns + 2, nr + 2))
        br[1:-1, 1:-1, :] = (
            als[1:, :, :] - als[:-1, :, :] + alp[:, :-1, :] - alp[:, 1:, :]
        )
        bs[1:-1, :, 1:-1] = alp[:, :, 1:] - alp[:, :, :-1]
        bp[:, 1:-1, 1:-1] = als[:, :, :-1] - als[:, :, 1:]

        del (alr, als, alp)

        # Fill ghost values with boundary conditions:
        # - constant gradient at outer boundary:
        bs[1:-1, :, -1] = 2 * bs[1:-1, :, -2] - bs[1:-1, :, -3]
        bp[:, 1:-1, -1] = 2 * bp[:, 1:-1, -2] - bp[:, 1:-1, -3]
        # - periodic in phi:
        bs[0, :, :] = bs[-2, :, :]
        bs[-1, :, :] = bs[1, :, :]
        br[0, :, :] = br[-2, :, :]
        br[-1, :, :] = br[1, :, :]
        # js = jp = 0 at photosphere:
        for i in range(nph + 1):
            bp[i, :, 0] = (
                Sbp[:, 0]
                / dnph[:, 0]
                * (
                    bp[i, :, 1] * dnph[:, 1] / Sbp[:, 1]
                    + br[i, :, 0] * dnr[:] / Sbr[:, 0]
                    - br[i + 1, :, 0] * dnr[:] / Sbr[:, 0]
                )
            )
        for i in range(nph + 2):
            bs[i, :, 0] = (
                Sbs[:, 0]
                / dns[:, 0]
                * (
                    bs[i, :, 1] * dns[:, 1] / Sbs[:, 1]
                    + br[i, :-1, 0] * dnr[:-1] / Sbr[:-1, 0]
                    - br[i, 1:, 0] * dnr[1:] / Sbr[1:, 0]
                )
            )
        # - polar boundaries as in dumfric:
        for i in range(nph + 2):
            i1 = (i + nph // 2) % nph
            br[i, -1, :] = br[i1, -2, :]
            br[i, 0, :] = br[i1, 1, :]
            bs[i, -1, :] = 0.5 * (bs[i, -2, :] - bs[i1, -2, :])
            bs[i, 0, :] = 0.5 * (bs[i, 1, :] - bs[i1, 1, :])
        for i in range(nph + 1):
            i1 = (i + nph // 2) % nph
            bp[i, -1, :] = -bp[i1, -2, :]
            bp[i, 0, :] = -bp[i1, 1, :]

        if output == "bc":
            # Remove area factors:
            for i in range(nph + 2):
                br[i, :, :] = br[i, :, :] / Sbr
                bs[i, :, :] = bs[i, :, :] / Sbs
            for i in range(nph + 1):
                bp[i, :, :] = bp[i, :, :] / Sbp

            return r, th, ph, rrc, thc, pc, br, bs, bp

        if output == "bg":
            # Weighted average to grid points:
            brg = br[:-1, :-1, :] + br[1:, :-1, :] + br[1:, 1:, :] + br[:-1, 1:, :]
            bsg = bs[:-1, :, :-1] + bs[1:, :, :-1] + bs[1:, :, 1:] + bs[:-1, :, 1:]
            bpg = bp[:, :-1, :-1] + bp[:, 1:, :-1] + bp[:, 1:, 1:] + bp[:, :-1, 1:]
            for i in range(nph + 1):
                brg[i, :, :] /= 2 * (Sbr[:-1, :] + Sbr[1:, :])
                bsg[i, :, :] /= 2 * (Sbs[:, :-1] + Sbs[:, 1:])
            for i in range(nph + 1):
                bpg[i, :, :] /= (
                    Sbp[:-1, :-1] + Sbp[1:, :-1] + Sbp[1:, 1:] + Sbp[:-1, 1:]
                )

            return r, th, ph, brg, bsg, bpg
