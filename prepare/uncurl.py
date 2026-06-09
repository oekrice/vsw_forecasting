import numpy as np
from scipy import sparse
import scipy.sparse.linalg as la
from scipy.ndimage import binary_dilation

def shrinkage(a, kappa):
    return np.maximum(0, a-kappa) - np.maximum(0, -a-kappa)

def basisPursuit(A, b, rho, alpha, max_iter=1000, abstol=1e-6, reltol=1e-6, lsqtol=1e-4, verbose=True):
    """
    Solves the following problem via ADMM:
        minimize    ||x||_1
        subject to  Ax = b
    The solution is returned in the vector x.
    
    rho is the augmented Lagrangian parameter.
    
    alpha is the over-relaxation parameter (typically 1 to 1.8).
    
    Method described in the paper Boyd et al. (2010).
    
    Ported to Python by ARY 22/8/17 (I use least-squares for the projection so that the method will work for larger systems).
    """

    m = np.size(A,0); n = np.size(A,1)
    
    x = np.zeros(n)
    z = np.zeros(n)
    u = np.zeros(n)
    
    if (verbose):
        print('iter   r norm   eps pri   s norm   eps dual   objective')
   
    abstol = np.sqrt(n)*abstol
   
    for k in range(max_iter):
        # x-update
        x = z - u
        x += la.lsqr(A, b - A.dot(x), atol=lsqtol, btol=lsqtol)[0]

        # z-update with relaxation
        zold = z.copy()
        x_hat = alpha*x + (1 - alpha)*zold
        z = shrinkage(x_hat + u, 1.0/rho)
    
        u += x_hat - z
        
        # termination checks:
        objval = np.linalg.norm(x,1)
        
        r_norm = np.linalg.norm(x-z,2)
        s_norm = np.linalg.norm(-rho*(z - zold),2)
        
        eps_pri = abstol + reltol*max(np.linalg.norm(x,2), np.linalg.norm(-z,2))
        eps_dual = abstol + reltol*np.linalg.norm(rho*u,2)
        
        if (verbose & (k%100 == 0)):
            print(k, r_norm, eps_pri, s_norm, eps_dual, objval)
        
        if ((r_norm < eps_pri) & (s_norm < eps_dual)):
            break
    
    return x

def e_global(ns, nph, dbr, inductive=False, rho=1e5, alpha=1.35):
    """
    Compute the 2d electric field on the dumfric grid given a map of dbr at
    cell centres.
    
    Computes L1-minimum E unless inductive=True, in which case it computes the
    L2-minimum.
    
    Parameters rho and alpha are for the basis pursuit algorithm.
    
    Output is es*Ls and ep*Lp on edges.
    """

    # Set up right-hand side (-dbr*Sr):
    ds = 2.0/ns
    dp = 2*np.pi/nph
    dbr = dbr.T
    b = -dbr.flatten(order='F')*ds*dp

    # Set up matrix:
    # - the matrix A has ns*np rows [grid cells] and 
    #   ns*(np+1) + np*(ns+1) columns [electric field components on edges].
    # - ep is listed first, in column major order starting from bottom left.
    # - the same matrix would be defined by
    #       A = n.zeros((ns*np, ns*(np+1) + np*(ns+1)))
    #       for i in range(ns): # s index
    #           for j in range(np): # p index
    #               A[(i)*np + j, i*(np) + j] = 1
    #               A[(i)*np + j, (i+1)*(np) + j] = -1
    #               A[(i)*np + j, i*(np+1) + j+1 + np*(ns+1)] = 1
    #               A[(i)*np + j, i*(np+1) + j + np*(ns+1)] = -1
    #
    I = sparse.eye(nph+1)
    Jp = (I.A[1:,:] - I.A[:-1,:])
    Ip = sparse.eye(nph)
    I = sparse.eye(ns+1)
    Js = I.A[1:,:] - I.A[:-1,:]
    del(I)
    Is = sparse.eye(ns)
    A = sparse.hstack([-sparse.kron(Js,Ip), sparse.kron(Is,Jp)])
    del(Ip,Is,Jp,Js)
   
    # Solve for E:
    if (inductive):
        x = la.lsqr(A, b)[0]
    else:
        x = basisPursuit(A, b, rho, alpha)

    # Repackage variables: 
    elp = np.reshape(x[:(nph*(ns+1))], (ns+1, nph))
    els = np.reshape(x[(nph*(ns+1)):], (ns, nph+1))

    return els, elp


def e_local(ns, nph, dbr, msk, inductive=False, rho=1e5, alpha=1.35):
    """
    Compute the 2d electric field on the dumfric grid given a map of dbr at
    cell centres and a mask array msk (0 or 1).

    The electric field is computed only on regions where msk=1, with Exn==0 on the
    boundary of each such region.

    Computes L1-minimum E unless inductive=True, in which case it computes the
    L2-minimum.

    Parameters rho and alpha are for the basis pursuit algorithm.

    Output is es*Ls and ep*Lp on edges of full grid.
    """

    # Label each grid cell with index:
    ncl = nph*ns
    icl = np.zeros((ns, nph), dtype='int')
    icl[icl == 0] = np.linspace(0, ncl-1, ncl)

    # Label each s rib and p rib with index:
    nrs0 = (nph+1)*(ns)
    i0rs = np.zeros((ns,nph+1), dtype='int')
    i0rs[i0rs == 0] = np.linspace(0, nrs0-1, nrs0)

    nrp0 = nph*(ns+1)
    i0rp = np.zeros((ns+1,nph), dtype='int')
    i0rp[i0rp == 0] = np.linspace(0, nrp0-1, nrp0)

    # Identify s and p ribs we want to solve for, and assign them "unknown id's":
    rs_slv = (msk[1:-1,:-1]*msk[1:-1,1:]) > 0
    rp_slv = (msk[:-1,1:-1]*msk[1:,1:-1]) > 0
    nrs = np.sum(rs_slv)
    nrp = np.sum(rp_slv)
    rp_uid = np.zeros((ns+1,nph), dtype='int') - 1
    rp_uid[rp_slv] = np.linspace(0, nrp-1, nrp)
    rs_uid = np.zeros((ns,nph+1), dtype='int') - 1
    rs_uid[rs_slv] = np.linspace(nrp, nrp+nrs-1, nrs)

    # Identify grid cells within mask region, and assign them "equation id's":
    eqid = np.zeros((ns, nph), dtype='int') - 1
    imsk = np.nonzero(msk[1:-1,1:-1] == 1)
    nmsk = np.size(imsk, axis=1)
    eqid[imsk] = np.linspace(0, nmsk-1, nmsk)

    # Initialize matrix:
    data = np.zeros((4*nmsk))
    row_ind = data.copy().astype('int')
    col_ind = data.copy().astype('int')

    # Add entries in matrix corresponding to interior cells:
    k = 0
    for i in range(1,ns-1):
        for j in range(1,nph-1):
            if (msk[i+1,j+1] == 1):
                if (msk[i,j+1] == 1):
                    row_ind[k] = eqid[i,j]
                    col_ind[k] = rp_uid[i,j]
                    data[k] = 1
                    k += 1
                if (msk[i+2,j+1] == 1):
                    row_ind[k] = eqid[i,j]
                    col_ind[k] = rp_uid[i+1,j]
                    data[k] = -1
                    k += 1
                if (msk[i+1,j] == 1):
                    row_ind[k] = eqid[i,j]
                    col_ind[k] = rs_uid[i,j]
                    data[k] = -1
                    k += 1
                if (msk[i+1,j+2] == 1):
                    row_ind[k] = eqid[i,j]
                    col_ind[k] = rs_uid[i,j+1]
                    data[k] = 1
                    k += 1

    # Form into sparse matrix:
    A = sparse.csc_matrix((data, (row_ind, col_ind)))

    # Set up right-hand side (-dbr*Sr):
    ds = 2.0/ns
    dp = 2*np.pi/nph
    b = np.zeros((nmsk))
    b[0:nmsk] = -dbr[imsk]*ds*dp

    # Solve for E:
    if (inductive):
        x = la.lsqr(A, b)[0]
    else:
        x = basisPursuit(A, b, rho, alpha)

    # Repackage variables in global arrays:
    elp = np.zeros((ns+1, nph))
    for i in range(ns+1):
        for j in range(nph):
            if (rp_uid[i,j] >= 0):
                elp[i,j] = x[rp_uid[i,j]]
    els = np.zeros((ns, nph+1))
    for i in range(ns):
        for j in range(nph+1):
            if (rs_uid[i,j] >= 0):
                els[i,j] = x[rs_uid[i,j]]

    return els, elp

def noninductive_local(ns, nph, f, msk):
    """
    Given an array f at cell centres and a mask array msk (0 or 1), determine a non-inductive electric field contribution grad(u), by solving the poisson equation
            lap_h(u) = f
    with Dirichlet boundary conditions u=0 on the boundary of msk.
    
    Output is es*Ls and ep*Lp on edges of full grid.
    """
    
    ds = 2.0/ns
    dp = 2*np.pi/nph
    
    # Add ghost points in phi for msk array (don't need them at right-hand end):
    msk1 = np.zeros((ns+2, nph+1))
    msk1[1:-1,1:] = msk
    msk1[1:-1,0] = msk[:,-1]
    
    # Add ghost points to f at cell centres (don't need them at right-hand end):
    f1 = np.zeros((ns+2, nph+1))
    f1[1:-1,1:] = f
    f1[1:-1,0] = f1[1:-1,-1]

    # Average f to grid points (not right-most):
    fg = 0.25*(f1[:-1,:-1] + f1[:-1,1:] + f1[1:,:-1] + f1[1:,1:])
    fg[:,0] = 0.25*(f1[:-1,-1] + f1[:-1,0] + f1[1:,-1] + f1[1:,0])
    
    # Arrays of 1-s^2:
    s = np.linspace(-1, 1, ns+1)
    # - at grid points in s:
    sig = 1 - s**2
    # - at midway points in s:
    sih = 1 - (0.5*(s[:-1] + s[1:]))**2

    # Identify particular grid points we want to solve for, and assign them "unknown id's":
    # (omit right-most column in phi and fill at the end by periodicity)
    g_slv = (msk1[:-1,:-1]*msk1[:-1,1:]*msk1[1:,:-1]*msk1[1:,1:]) > 0
    ngs = np.sum(g_slv)
    g_uid = np.zeros((ns+1,nph), dtype='int') - 1   # at grid points
    g_uid[g_slv] = np.linspace(0, ngs-1, ngs)

    # Initialize matrix:
    data = np.zeros((5*ngs))
    row_ind = data.copy().astype('int')
    col_ind = data.copy().astype('int')
    
    # Add entries in matrix corresponding to unknown grid points:
    k = 0
    for i in range(1,ns):
        for j in range(0,nph):
            if (g_uid[i,j] > -1):
                row_ind[k] = g_uid[i,j]
                col_ind[k] = g_uid[i,j]
                data[k] = -sig[i]*(sih[i-1] + sih[i])/ds**2 - 2/dp**2
                k += 1
                if (g_uid[i+1,j] > -1):
                    row_ind[k] = g_uid[i,j]
                    col_ind[k] = g_uid[i+1,j]
                    data[k] = sig[i]*sih[i]/ds**2
                    k += 1
                if (g_uid[i-1,j] > -1):
                    row_ind[k] = g_uid[i,j]
                    col_ind[k] = g_uid[i-1,j]
                    data[k] = sig[i]*sih[i-1]/ds**2
                    k += 1
                jp = (j + 1) % nph
                if (g_uid[i,jp] > -1):
                    row_ind[k] = g_uid[i,j]
                    col_ind[k] = g_uid[i,jp]
                    data[k] = 1/dp**2
                    k += 1
                jm = (j + nph - 1) % nph
                if (g_uid[i,jm] > -1):
                    row_ind[k] = g_uid[i,j]
                    col_ind[k] = g_uid[i,jm]
                    data[k] = 1/dp**2
                    k += 1

    # Form into sparse matrix:
    A = sparse.csc_matrix((data, (row_ind, col_ind)))

    # Set up right-hand side (-jz):
    b = np.zeros((ngs))
    k = 0
    for i in range(1,ns):
        for j in range(0,nph):
            if (g_uid[i,j] > -1):
                b[k] = sig[i]*fg[i,j]
                k += 1
                
    # Solve for psi at interior grid points:
    x = la.lsqr(A, b)[0]
          
    # Repackage in global array:
    psi = np.zeros((ns+1, nph+1))
    for i in range(ns+1):
        for j in range(nph):
            if (g_uid[i,j] > -1):
                psi[i,j] = x[g_uid[i,j]]

    # Add right-most values by periodicity
    for i in range(ns+1):
        if (g_uid[i,0] >= 0):
            psi[i,-1] = x[g_uid[i,0]]
            
    # Return components of grad(psi) on edges, multiplied by edge lengths:
    dpsis = psi[1:,:] - psi[:-1,:]
    dpsip = psi[:,1:] - psi[:,:-1]
    
    return dpsis, dpsip, psi

