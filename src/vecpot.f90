!*******************************************************************************
MODULE vecpot
!*******************************************************************************
! This module implements the magnetic vector potential.
!*******************************************************************************
    USE ffts
    USE params
    USE mpitools, ONLY: control, mpi_loc, mpi_dims, comm, mpireal, stat, &
        err, outer, surface, north, south, mpi_r_nxt, mpi_r_prv, mpi_p_nxt, &
        mpi_p_prv
    USE grid, ONLY: nr, ns, np, nsg, npg, dr, dp, ds, dns, dnp, Sbr, Sbs, Sbp, &
        dlr, dls, dlp, rc, sc, pc, rg, sg, pg, rss
    USE ncio, ONLY: ncdf_read_3darray

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: alr, als, alp, read_vector_potential, compute_a_from_b
    PUBLIC:: arefs, arefp, zero_aref_outer, compute_aref_outer
    PUBLIC:: set_a_exact

!*******************************************************************************
    REAL(d), ALLOCATABLE:: alr(:,:,:)   ! A_rho*L_rho on edges (1:nr,0:ns,0:np)
    REAL(d), ALLOCATABLE:: als(:,:,:)   ! A_s*L_s on edges (0:nr,1:ns,0:np)
    REAL(d), ALLOCATABLE:: alp(:,:,:)   ! A_phi*L_phi on edges (0:nr,0:ns,1:np)
    REAL(d), ALLOCATABLE:: arefs(:,:)   ! A^ref_s on outer boundary edges
                                        ! (1:ns,0:np)
    REAL(d), ALLOCATABLE:: arefp(:,:)   ! A^ref_phi on outer boundary edges
                                        ! (0:ns,1:np)
    REAL(d), ALLOCATABLE:: u(:,:)
    REAL(d), ALLOCATABLE:: ug(:,:)

    INTEGER, PARAMETER:: tag=1

!*******************************************************************************
CONTAINS
!*******************************************************************************

!===============================================================================
SUBROUTINE read_vector_potential(netcdf_file)
! Read components of A*L on cell edges from file.
    CHARACTER(*):: netcdf_file     ! netcdf file containing alr, als, alp.
    
    ! Allocate arrays:
    IF (.NOT. ALLOCATED(alr)) ALLOCATE(alr(1:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(als)) ALLOCATE(als(0:nr,1:ns,0:np))
    IF (.NOT. ALLOCATED(alp)) ALLOCATE(alp(0:nr,0:ns,1:np))

    ! Read arrays:
    CALL ncdf_read_3darray(netcdf_file, 'ar', alr(1:nr,0:ns,0:np))
    CALL ncdf_read_3darray(netcdf_file, 'as', als(0:nr,1:ns,0:np))
    CALL ncdf_read_3darray(netcdf_file, 'ap', alp(0:nr,0:ns,1:np))

    IF (control) PRINT*,'Read vector potential from file '//netcdf_file

END SUBROUTINE read_vector_potential

!===============================================================================
SUBROUTINE compute_a_from_b(br, bs, bp)
    ! Given existing arrays of br, bs, bp (excluding face areas), compute A*L on
    ! cell edges.
    !
    ! Used when restarting the code from a magnetic snapshot.
    !
    ! Uses the DeVore gauge where ar = 0, with "line integral gauge" on
    ! lower boundary.
    ! - on lower boundary:
    !     as0(s,p) = int_p0^p br
    !     ap0(s,p) = 0
    ! - in volume:
    !     as(r,s,p) = as0(s,p) - int_r0^r bp
    !     ap(r,s,p) = ap0(s,p) + int_r0^r bs
    REAL(d), INTENT(IN):: br(0:nr,0:ns+1,0:np+1)   ! B_rho on faces
    REAL(d), INTENT(IN):: bs(0:nr+1,0:ns,0:np+1)   ! B_s on faces
    REAL(d), INTENT(IN):: bp(0:nr+1,0:ns+1,0:np)   ! B_phi on faces
    
    ! Local variables:
    INTEGER :: i,j,k
    
    ! Allocate arrays:
    IF (.NOT. ALLOCATED(alr)) ALLOCATE(alr(1:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(als)) ALLOCATE(als(0:nr,1:ns,0:np))
    IF (.NOT. ALLOCATED(alp)) ALLOCATE(alp(0:nr,0:ns,1:np))
    
    alr(1:nr,0:ns,0:np) = 0.0_d
    als(0:nr,1:ns,0:np) = 0.0_d
    alp(0:nr,0:ns,1:np) = 0.0_d
    
    ! Compute as on lower boundary edges by integration in p direction:
    DO j = 0, mpi_dims(3) - 2
        IF (surface) THEN
            IF (mpi_loc(3) == j) THEN
                DO k = 1, np
                    als(0,1:ns,k) = als(0,1:ns,k-1) &
                        + br(0,1:ns,k) * Sbr(0,1:ns,k)
                END DO
                CALL MPI_SEND(als(0,1:ns,np), ns, mpireal, mpi_p_nxt, &
                    tag + j, comm, err)
            END IF
            IF (mpi_loc(3) == (j + 1)) THEN
                CALL MPI_RECV(als(0,1:ns,0), ns, mpireal, mpi_p_prv, &
                    tag + j, comm, stat, err)
                DO k = 1, np
                    als(0,1:ns,k) = als(0,1:ns,k-1) &
                        + br(0,1:ns,k) * Sbr(0,1:ns,k)
                END DO
            END IF
        END IF
    END DO

    ! Loop through MPI r-layers and do integration in order:
    DO i = 0, mpi_dims(1)-1
        IF (mpi_loc(1) == i) THEN
            DO j = 1, nr
                als(j,1:ns,0:np) = als(j-1,1:ns,0:np) &
                    - bp(j,1:ns,0:np) * Sbp(j,1:ns,0:np)
                alp(j,0:ns,1:np) = alp(j-1,0:ns,1:np) &
                    + bs(j,0:ns,1:np) * Sbs(j,0:ns,1:np)
                ! Deal with polar boundary:
                IF (south) alp(j,0,1:np) = alp(j-1,0,1:np) &
                    + bs(j,1,1:np) * Sbs(j,1,1:np) + bp(j,1,1:np) * Sbp(j,1,1:np) &
                    - bp(j,1,0:np-1) * Sbp(j,1,0:np-1) + br(j,1,1:np) * Sbr(j,1,1:np) &
                    - br(j-1,1,1:np) * Sbr(j-1,1,1:np)
                IF (north) alp(j,ns,1:np) = alp(j-1,ns,1:np) &
                    + bs(j,ns-1,1:np) * Sbs(j,ns-1,1:np) &
                    - bp(j,ns,1:np) * Sbp(j,ns,1:np) &
                    + bp(j,ns,0:np-1) * Sbp(j,ns,0:np-1) &
                    - br(j,ns,1:np) * Sbr(j,ns,1:np) &
                    + br(j-1,ns,1:np) * Sbr(j-1,ns,1:np)
            END DO
            IF (.NOT. outer) THEN
                CALL MPI_SEND(als(nr,1:ns,0:np), ns * (np + 1), mpireal, &
                    mpi_r_nxt, tag + i, comm, err)
                CALL MPI_SEND(alp(nr,0:ns,1:np), (ns + 1) * np, &
                    mpireal, mpi_r_nxt, tag + i, comm, err)
            END IF
        END IF
        IF (mpi_loc(1) == (i+1)) THEN
            CALL MPI_RECV(als(0,1:ns,0:np), ns * (np + 1), mpireal, &
                mpi_r_prv, tag + i, comm, stat, err)
            CALL MPI_RECV(alp(0,0:ns,1:np),(ns + 1) * np, mpireal, mpi_r_prv, &
                tag + i, comm, stat, err)
        END IF
    END DO
    
END SUBROUTINE compute_a_from_b

!===============================================================================
SUBROUTINE zero_aref_outer()
! Allocate A_s and A_p on outer boundary if necessary, and set to zero.
! For code reasons, need to allocate even for processes not on outer boundary,
! but these are just left at zero.

    IF (.NOT. ALLOCATED(arefs)) ALLOCATE(arefs(1:ns,0:np))
    IF (.NOT. ALLOCATED(arefp)) ALLOCATE(arefp(0:ns,1:np))
    arefs(1:ns,0:np) = 0.0_d
    arefp(0:ns,1:np) = 0.0_d

END SUBROUTINE zero_aref_outer

!===============================================================================
SUBROUTINE compute_aref_outer(br)
! Compute reference A_s and A_p on outer boundary from B_rho. In minimal gauge.
! Uses FFT method [currently single process gathers all arrays and does the
! computation].
    REAL(d), INTENT(IN):: br(0:nr,0:ns+1,0:np+1)    ! B_rho on faces.

    ! Local variables:
    REAL(d):: e1
    INTEGER:: i, j, i0, j0, i1, j1, srank

    
    ! Control gathers -br into global array ug
    i0 = 0
    j0 = 0
    IF (control) THEN
        IF (.NOT. ALLOCATED(ug)) ALLOCATE(ug(0:nsg+1, 0:npg+1))
        ug = 0.0_d
        IF (outer) THEN
            i0 = mpi_loc(2) * ns
            j0 = mpi_loc(3) * np
            ug(i0:i0+ns+1,j0:j0+np+1) = -br(nr,0:ns+1,0:np+1) &
                * Sbr(nr,0:ns+1,0:np+1)
        END IF
    END IF
    DO i = 0, mpi_dims(2) - 1
        DO j = 0, mpi_dims(3) - 1
            IF ((outer) .AND. (mpi_loc(2) == i) .AND. (mpi_loc(3) == j) &
                .AND. (.NOT. control)) &
                CALL MPI_SEND(-br(nr,0:ns+1,0:np+1) * Sbr(nr,0:ns+1,0:np+1), &
                    (ns+2) * (np+2), mpireal, 0, tag + i * j, comm, err)
            IF (control) THEN
                CALL MPI_CART_RANK(comm, (/mpi_dims(1) - 1, i, j/), &
                    srank, err)
                IF (srank /= 0) THEN
                    i1 = i * ns
                    j1 = j * np
                    CALL MPI_RECV(ug(i1:i1+ns+1,j1:j1+np+1), (ns+2) * (np+2), &
                        mpireal, srank, tag + i * j, comm, stat, err)
                END IF
            END IF
        END DO
    END DO

    ! Control solves 2D Poisson equation lap_h(u) = -br at cell centres,
    ! result in ug.
    IF (control) CALL poisson2d(ug)

    ! Control scatters ug to u, including ghost values.
    IF (.NOT. ALLOCATED(u)) ALLOCATE(u(0:ns+1,0:np+1))
    IF (control .AND. outer) u(0:ns+1,0:np+1) = ug(i0:i0+ns+1,j0:j0+np+1)
    DO i = 0, mpi_dims(2) - 1
        DO j = 0, mpi_dims(3) - 1
            IF (control) THEN
                CALL MPI_CART_RANK(comm, (/mpi_dims(1) - 1, i, j/), srank, err)
                IF (srank /= 0) THEN
                    i1 = i * ns
                    j1 = j * np
                    CALL MPI_SEND(ug(i1:i1+ns+1,j1:j1+np+1), (ns+2) * (np+2), &
                        mpireal, srank, tag + i * j, comm, err)
                END IF
            END IF
            IF (outer .AND. (mpi_loc(2) == i) .AND. (mpi_loc(3) == j) .AND. &
                (.NOT. control)) CALL MPI_RECV(u(0:ns+1,0:np+1), &
                    (ns+2) * (np+2), mpireal, 0, tag + i * j, comm, stat, err)
        END DO
    END DO

    ! Locals compute arefs and arefp by A^ref = curl(u*e_r)
    IF (outer) THEN
        e1 = EXP(0.5_d * dr)
        arefs(1:ns,0:np) = -(u(1:ns,1:np+1) - u(1:ns,0:np)) &
            / dnp(nr,1:ns,0:np) / e1
        arefp(0:ns,1:np) = (u(1:ns+1,1:np) - u(0:ns,1:np)) &
            / dns(nr,0:ns,1:np) / e1
    END IF
    
END SUBROUTINE compute_aref_outer

!===============================================================================
SUBROUTINE poisson2d(ug)
! Solve Poisson equation in-place for ug at cell centres on a spherical surface.
! i.e. ug starts out as RHS and ends as solution, both including ghost cells.
!
! Uses FFT method [currently single process gathers all arrays and does
! computation]
    REAL(d), INTENT(INOUT):: ug(0:nsg+1,0:npg+1)
    
    ! Local variables:
    REAL(d):: sg(0:nsg), Vg(0:nsg)
    REAL(d):: sgc(1:nsg), Uc(1:nsg), ev(1:nsg), Q(1:nsg,1:nsg), lam(1:nsg)
    REAL(d):: mu(1:npg)
    COMPLEX(d):: ugz(1:nsg,1:npg), clm(1:nsg)
    INTEGER:: i, m

    ! Global coordinate factors:
    sg(0:nsg) = (/ (-1.0_d + DBLE(i) * ds, i = 0, nsg) /)
    sgc(1:nsg) = (/ (-1.0_d + (DBLE(i) - 0.5_d) * ds, i = 1, nsg) /)
    Vg(0:nsg) = 0.0_d
    Vg(1:nsg-1) = dp * SQRT(1.0_d - sg(1:nsg-1)**2) / (ASIN(sgc(2:nsg)) &
        - ASIN(sgc(1:nsg-1)))
    Uc(1:nsg) = (ASIN(sg(1:nsg)) - ASIN(sg(0:nsg-1))) / SQRT(1.0_d - sgc**2) / dp

    ! FFT of array in p direction:
    ugz(1:nsg,1:npg) = CMPLX(-ug(1:nsg,1:npg), KIND=d)
    CALL fftn(ugz, SHAPE(ugz), DIM=(/2/))
    ugz(1:nsg,1:npg) = ugz/SQRT(DBLE(npg))

    ! Frequencies [+,-] for FFT:
    mu(1:npg) = (/ (DBLE(i), i=0,npg-1) /)
    mu(1:npg) = mu(1:npg) / DBLE(npg)
    mu(npg/2+1:npg) = mu(npg/2+1:npg) - 1.0_d

    ! Loop over azimuthal modes (positive m):
    mu(1:npg) = 4.0_d * SIN(PI * mu(1:npg))**2
    DO m = 1, npg/2 + 1
        ! - prepare tridiagonal matrix (ev = off-diagonal and lam diagonal).
        ! - also prepare identity matrix Q
        ev(1:nsg) = -Vg(0:nsg-1)
        Q(1:nsg,1:nsg) = 0.0_d
        DO i = 1, nsg
            lam(i) = Vg(i-1) + Vg(i) + Uc(i) * mu(m)
            Q(i,i) = 1.0_d
        END DO
        ! - compute eigenvectors Q_{lm} and eigenvalues lam_{lm}
        CALL tri_eig(lam, ev, Q)
        ! - compute array of c_{lm} for each l
        DO i = 1,nsg
            clm(i) = DOT_PRODUCT(Q(1:nsg,i), ugz(1:nsg,m)) / lam(i)
        END DO
        ! - compute entry for this m in psit = Sum_l c_{lm}Q_{lm}^j
        ugz(1:nsg,m) = MATMUL(Q, clm)
        IF (m > 1) ugz(1:nsg,npg+2-m) = CONJG(ugz(1:nsg,m))
    END DO

    ! Inverse FFT:
    CALL fftn(ugz, SHAPE(ugz), DIM=(/2/), INV=.TRUE.)
    ug(1:nsg,1:npg) = DBLE(ugz) * SQRT(DBLE(npg))

    ! Global boundary conditions [0 at poles, periodic in phi]:
    ug(0,0:npg+1) = 0.0_d
    ug(nsg+1,0:npg+1) = 0.0_d
    ug(0:nsg+1,0) = ug(0:nsg+1,npg)
    ug(0:nsg+1,npg+1) = ug(0:nsg+1,1)

END SUBROUTINE poisson2d

!===============================================================================
FUNCTION pythag(a, b)
    REAL(d):: a, b, pythag
    ! ---
    ! (from Numerical Recipes)
    ! Compute (a^2 + b^2)^1/2 without destructive underflow or overflow.
    ! ---
    REAL(d):: absa, absb
    
    absa = ABS(a)
    absb = ABS(b)
    IF (absa > absb) THEN
        pythag = absa * SQRT(1.0_d + (absb/absa)**2)
    ELSE
        IF (absb < EPSILON(1.0_d)) THEN
            pythag = 0.0_d
        ELSE
            pythag = absb * SQRT(1.0_d + (absa/absb)**2)
        END IF
    END IF
    
END FUNCTION pythag

!===============================================================================
SUBROUTINE tri_eig(dv, e, z)
! (tqli from Numerical Recipes - this one based on C version to avoid
!  named DO loops)
! Compute eigenvalues and eigenvectors of a symm. tridiagonal matrix A, using
! the QL algorithm. To find only eigenvalues, omit optional argument z.
! - input: dv(n) - a vector with the diagonal elements of A
!          e(n) - a vector with the subdiagonal of A (ignores e(1))
!          z(n,n) - the identity matrix [optional - if eigenvectors reqd]
! - output: dv(n) - a vector with the eigenvalues
!           z(n,n) - the matrix of eigenvectors - kth column is
!                    normalized eigenvector corresponding to dv(k).
    REAL(d), INTENT(INOUT):: dv(:), e(:)
    REAL(d), OPTIONAL, INTENT(INOUT):: z(:,:)
    
    ! Local variables:
    INTEGER:: i, iter, l, m, n, ndum
    REAL(d):: b, c, dd, f, g, p, r, s
    REAL(d):: ff(SIZE(e))
    LOGICAL:: cflag
    
    n = SIZE(dv)
    IF (PRESENT(z)) ndum = SIZE(z, 1)
    e(:) = EOSHIFT(e(:), 1)  ! convenient to renumber the elements of e
    DO l = 1, n
        iter = 0
        DO
            DO m = l, n - 1  ! Look for a single small subdiagonal element
                ! to split the matrix
                dd = ABS(dv(m)) + ABS(dv(m+1))
                IF (ABS(e(m)) < EPSILON(1.0_d) * dd) EXIT
            END DO
            IF (m == l) EXIT
            IF (iter == 30) THEN
                PRINT*,'Error: too many iterations in tqli'
                STOP
            END IF
            iter = iter+1
            g = (dv(l+1) - dv(l)) / (2.0_d * e(l))  ! form shift
            r = pythag(g, 1.0_d)
            g = dv(m) - dv(l) + e(l) / (g + SIGN(r, g))  ! This is dv_m - k_s
            s = 1.0_d
            c = 1.0_d
            p = 0.0_d
            cflag = .FALSE.
            DO i = m-1, l, -1  ! A plane rotation as in original QL, followed by
                ! Givens rotations to restore tridiagonal form.
                f = s * e(i)
                b = c * e(i)
                r = pythag(f, g)
                e(i+1) = r
                IF (ABS(r) < EPSILON(1.0_d)) THEN  ! recover from underflow
                    dv(i+1) = dv(i+1) - p
                    e(m) = 0.0_d
                    cflag = .TRUE.
                    EXIT
                END IF
                s = f / r
                c = g / r
                g = dv(i+1) - p
                r = (dv(i) - g) * s + 2.0_d * c * b
                p = s * r
                dv(i+1) = g + p
                g = c * r - b
                IF (PRESENT(z)) THEN  ! for eigenvectors
                    ff(1:n) = z(1:n,i+1)
                    z(1:n,i+1) = s * z(1:n,i) + c * ff(1:n)
                    z(1:n,i) = c * z(1:n,i) - s * ff(1:n)
                END IF
            END DO
            IF (cflag) CYCLE
            dv(l) = dv(l) - p
            e(l) = g
            e(m) = 0.0_d
        END DO
    END DO

END SUBROUTINE tri_eig

!===============================================================================
SUBROUTINE a_exact(r, s, p, ar0, as0, ap0, test_case)
! Return exact solution for test field, ar, as, ap components (no weighting).
! Integer test_case specifies which exact solution.
    REAL(d), INTENT(IN):: r, s, p
    REAL(d), INTENT(OUT):: ar0, as0, ap0
    INTEGER, INTENT(IN):: test_case

    ! Local variables:
    REAL(d):: rad, sintheta, cos2theta

    rad = EXP(r)
    
    IF (test_case == 0) THEN    ! axisymmetric PFSS dipole
        ar0 = 0.0_d
        as0 = 0.0_d
        sintheta = SQRT(1.0_d - s**2)
        cos2theta = 2*s**2 - 1
        ap0 = rad * sintheta * (1.0_d + 2.0_d * rss**3 / rad**3) &
            / (1 + 2 * rss**3)
    END IF

END SUBROUTINE a_exact

!===============================================================================
SUBROUTINE set_a_exact(test_case)
! Set alr, als, alp to analytical test field. Including edge lengths.
    INTEGER, INTENT(IN):: test_case

    ! Local variables:
    INTEGER:: i, j, k
    REAL(d):: sintheta, rad, cos2theta, ar0, as0, ap0

    IF (.NOT. ALLOCATED(alr)) ALLOCATE(alr(1:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(als)) ALLOCATE(als(0:nr,1:ns,0:np))
    IF (.NOT. ALLOCATED(alp)) ALLOCATE(alp(0:nr,0:ns,1:np))
    
    alr(1:nr,0:ns,0:np) = 0.0_d
    als(0:nr,1:ns,0:np) = 0.0_d
    alp(0:nr,0:ns,1:np) = 0.0_d

    DO k = 0, np
        DO j = 0, ns
            DO i = 1, nr
                CALL a_exact(rc(i), sg(j), pg(k), ar0, as0, ap0, test_case)
                alr(i,j,k) = ar0 * dlr(i,j,k)
            END DO
        END DO
    END DO

    DO k = 0, np
        DO j = 1, ns
            DO i = 0, nr
                CALL a_exact(rg(i), sc(j), pg(k), ar0, as0, ap0, test_case)
                als(i,j,k) = as0 * dls(i,j,k)
            END DO
        END DO
    END DO

    DO k = 1, np
        DO j = 0, ns
            DO i = 0, nr
                CALL a_exact(rg(i), sg(j), pc(k), ar0, as0, ap0, test_case)
                alp(i,j,k) = ap0 * dlp(i,j,k)
            END DO
        END DO
    END DO
    IF (south) alp(0:nr,0,1:np) = 0.0_d
    IF (north) alp(0:nr,ns,1:np) = 0.0_d

END SUBROUTINE set_a_exact

!*******************************************************************************
END MODULE vecpot
!*******************************************************************************
