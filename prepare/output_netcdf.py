"""
    Routines for writing netcdf files, as required for DuMFriC.
    
    ARY 24/8/17
"""

from scipy.io import netcdf
import numpy as n


def surfaceFlows(filename, eta, vs, vp):
    """
        Surface flow profiles (eta is a constant, vs on the p ribs
        and vp on the s ribs, both functions of s only).
    """
    ns = n.size(vp)

    fid = netcdf.netcdf_file(filename, "w")
    fid.createDimension("edim", 1)
    eid = fid.createVariable("etasurf", "d", ("edim",))
    eid[:] = eta
    fid.createDimension("sgdim", ns + 1)
    fid.createDimension("scdim", ns)
    vsid = fid.createVariable("vs", "d", ("sgdim",))
    vsid[:] = vs
    vpid = fid.createVariable("vp", "d", ("scdim",))
    vpid[:] = vp
    fid.close()


def surfaceFlows2d(filename, vs, vp):
    """
        Surface flow profiles (vs on the p ribs and vp on the s ribs).
    """
    ns = n.size(vp, axis=0)
    np = n.size(vs, axis=1)

    fid = netcdf.netcdf_file(filename, "w")
    fid.createDimension("scdim", ns)
    fid.createDimension("sgdim", ns + 1)
    fid.createDimension("pcdim", np)
    fid.createDimension("pgdim", np + 1)
    vid = fid.createVariable("vs", "d", ("pcdim", "sgdim"))
    vid[:] = vs.T
    vid = fid.createVariable("vp", "d", ("pgdim", "scdim"))
    vid[:] = vp.T
    fid.close()


def friction(filename, nug):
    """
        Friction coefficient as a function of r and s.
    """

    nr = n.size(nug, axis=0) - 1
    ns = n.size(nug, axis=1) - 1

    fid = netcdf.netcdf_file(filename, "w")
    fid.createDimension("rgdim", nr + 1)
    fid.createDimension("sgdim", ns + 1)
    nuid = fid.createVariable("nug", "d", ("sgdim", "rgdim"))
    nuid[:] = nug.T
    fid.close()


def outflow(filename, vr):
    """
        Radial outflow speed as a function of r [in Rsun/s].
    """
    nr = n.size(vr) - 1
    fid = netcdf.netcdf_file(filename, "w")
    fid.createDimension("rgdim", nr + 1)
    vrid = fid.createVariable("vrg", "d", ("rgdim",))
    vrid[:] = vr
    fid.close()


def a(filename, r, th, ph, apr, aps, app):
    """
        Vector potential * edge lengths on cell edges.
    """

    nr = n.size(r) - 1
    ns = n.size(th) - 1
    np = n.size(ph) - 1

    fid = netcdf.netcdf_file(filename, "w")
    fid.createDimension("rc", nr)
    fid.createDimension("r", nr + 1)
    fid.createDimension("thc", ns)
    fid.createDimension("th", ns + 1)
    fid.createDimension("phc", np)
    fid.createDimension("ph", np + 1)
    vid = fid.createVariable("r", "d", ("r",))
    vid[:] = r
    vid = fid.createVariable("th", "d", ("th",))
    vid[:] = th
    vid = fid.createVariable("ph", "d", ("ph",))
    vid[:] = ph
    vid = fid.createVariable("ar", "d", ("ph", "th", "rc"))
    vid[:] = apr
    vid = fid.createVariable("as", "d", ("ph", "thc", "r"))
    vid[:] = aps
    vid = fid.createVariable("ap", "d", ("phc", "th", "r"))
    vid[:] = app
    fid.close()
    print("Wrote A*L to file " + filename)


def bc(filename, r, th, ph, rc, thc, phc, br, bs, bp):
    """
        Magnetic field components on cell faces, including ghost cells.
    """

    nr = n.size(r) - 1
    ns = n.size(th) - 1
    np = n.size(ph) - 1

    fid = netcdf.netcdf_file(filename, "w")
    fid.createDimension("rc", nr + 2)
    fid.createDimension("r", nr + 1)
    fid.createDimension("thc", ns + 2)
    fid.createDimension("th", ns + 1)
    fid.createDimension("phc", np + 2)
    fid.createDimension("ph", np + 1)
    vid = fid.createVariable("r", "d", ("r",))
    vid[:] = r
    vid = fid.createVariable("th", "d", ("th",))
    vid[:] = th
    vid = fid.createVariable("ph", "d", ("ph",))
    vid[:] = ph
    vid = fid.createVariable("rc", "d", ("rc",))
    vid[:] = rc
    vid = fid.createVariable("thc", "d", ("thc",))
    vid[:] = thc
    vid = fid.createVariable("phc", "d", ("phc",))
    vid[:] = phc
    vid = fid.createVariable("br", "d", ("phc", "thc", "r"))
    vid[:] = br
    vid = fid.createVariable("bth", "d", ("phc", "th", "rc"))
    vid[:] = -bs
    vid = fid.createVariable("bph", "d", ("ph", "thc", "rc"))
    vid[:] = bp
    fid.close()
    print("Wrote B on faces to file " + filename)


def bg(filename, r, th, ph, brg, bsg, bpg):
    """
        Magnetic field components co-located at grid points.
    """

    nr = n.size(r) - 1
    ns = n.size(th) - 1
    np = n.size(ph) - 1

    fid = netcdf.netcdf_file(filename, "w")
    fid.createDimension("r", nr + 1)
    fid.createDimension("th", ns + 1)
    fid.createDimension("ph", np + 1)
    vid = fid.createVariable("r", "d", ("r",))
    vid[:] = r
    vid = fid.createVariable("th", "d", ("th",))
    vid[:] = th
    vid = fid.createVariable("ph", "d", ("ph",))
    vid[:] = ph
    vid = fid.createVariable("br", "d", ("ph", "th", "r"))
    vid[:] = brg
    vid = fid.createVariable("bth", "d", ("ph", "th", "r"))
    vid[:] = -bsg
    vid = fid.createVariable("bph", "d", ("ph", "th", "r"))
    vid[:] = bpg
    fid.close()
    print("Wrote B at grid points to file " + filename)


def esurf(filename, es, ep):
    """
        Electric field * edge lengths on 2D lower boundary.
    """

    ns = n.size(es, axis=0)
    np = n.size(ep, axis=1)

    fid = netcdf.netcdf_file(filename, "w")
    fid.createDimension("scdim", ns)
    fid.createDimension("sgdim", ns + 1)
    fid.createDimension("pcdim", np)
    fid.createDimension("pgdim", np + 1)
    vid = fid.createVariable("es", "d", ("pgdim", "scdim"))
    vid[:] = es.T
    vid = fid.createVariable("ep", "d", ("pcdim", "sgdim"))
    vid[:] = ep.T
    fid.close()
    print("Wrote E*L on lower boundary to file " + filename)


def brsurf(filename, br):
    """
        Br at interior cell centres on 2D lower boundary.
    """

    ns = n.size(br, axis=0)
    np = n.size(br, axis=1)

    fid = netcdf.netcdf_file(filename, "w")
    fid.createDimension("scdim", ns)
    fid.createDimension("pcdim", np)
    vid = fid.createVariable("br", "d", ("pcdim", "scdim"))
    vid[:] = br.T
    fid.close()
    print("Wrote Br on lower boundary to file " + filename)


def region(filename, als, alp):
    """
        Vector potential * edge lengths of emerging 2d region.
    """

    ns = n.size(als, axis=1)
    np = n.size(alp, axis=0)

    fid = netcdf.netcdf_file(filename, "w")
    fid.createDimension("scdim", ns)
    fid.createDimension("sgdim", ns + 1)
    fid.createDimension("pcdim", np)
    fid.createDimension("pgdim", np + 1)
    vid = fid.createVariable("as", "d", ("pgdim", "scdim"))
    vid[:] = als
    vid = fid.createVariable("ap", "d", ("pcdim", "sgdim"))
    vid[:] = alp
    fid.close()
    print("Wrote A*L for emerging region to file " + filename)
