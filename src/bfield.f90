!*******************************************************************************
MODULE bfield
!*******************************************************************************
! This module implements the magnetic field.
!*******************************************************************************
    USE params
    USE mpitools, ONLY: control, mpireal, comm, stat, err, mpi_opp, &
        mpi_r_prv, mpi_r_nxt, mpi_s_prv, mpi_s_nxt, mpi_p_prv, mpi_p_nxt, &
        surface, outer, south, north, MPI_MAX
    USE grid, ONLY: nr, ns, np, dns, dnr, dnp, Sbr, Sbs, Sbp, Sbrg, Sbsg, &
        Sbpg, rss, rg, rc, sc, sg, pc, pg
    USE hyper, ONLY: hyper_flag
    USE ncio, ONLY: ncdf_read_3darray, ncdf_init_snapshot_file, &
        ncdf_write_3darray

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: br, bs, bp, bgr, bgs, bgp, bbg, bbgs
    PUBLIC:: read_magnetic_field, compute_b_stokes
    PUBLIC:: test_compute_b_stokes, b_exact, set_b_exact

!*******************************************************************************
    REAL(d), ALLOCATABLE:: br(:,:,:)   ! B_rho on faces (0:nr,0:ns+1,0:np+1)
    REAL(d), ALLOCATABLE:: bs(:,:,:)   ! B_s on faces (0:nr+1,0:ns,0:np+1)
    REAL(d), ALLOCATABLE:: bp(:,:,:)   ! B_phi on faces (0:nr+1,0:ns+1,0:np)
    REAL(d), ALLOCATABLE:: bgr(:,:,:)   ! B_rho at grid pts (0:nr,0:ns,0:np)
    REAL(d), ALLOCATABLE:: bgs(:,:,:)   ! B_s at grid pts (0:nr,0:ns,0:np)
    REAL(d), ALLOCATABLE:: bgp(:,:,:)   ! B_phi at grid pts (0:nr,0:ns,0:np)
    REAL(d), ALLOCATABLE:: bbg(:,:,:)   ! B^2 at grid pts (-1:nr+1,0:ns,0:np)
    REAL(d), ALLOCATABLE:: bbgs(:,:,:)   ! Softened/bounded value of B^2 at
                                         ! grid pts (-1:nr+1,0:ns,0:np)

    INTEGER, PARAMETER:: tag=1

!*******************************************************************************
CONTAINS
!*******************************************************************************

!===============================================================================
SUBROUTINE read_magnetic_field(netcdf_file)
! Read face-centered components B_r, B_theta, B_phi from file and assign to
! B_rho, B_s, B_phi arrays [called br, bs, bp].
    CHARACTER(*):: netcdf_file     ! netcdf file containing br, bth, bph
    
    IF (.NOT. ALLOCATED(br)) ALLOCATE(br(0:nr,0:ns+1,0:np+1))
    IF (.NOT. ALLOCATED(bs)) ALLOCATE(bs(0:nr+1,0:ns,0:np+1))
    IF (.NOT. ALLOCATED(bp)) ALLOCATE(bp(0:nr+1,0:ns+1,0:np))

    ! Read arrays:
    CALL ncdf_read_3darray(netcdf_file, 'br', br(0:nr,0:ns+1,0:np+1))
    CALL ncdf_read_3darray(netcdf_file, 'bth', bs(0:nr+1,0:ns,0:np+1))
    CALL ncdf_read_3darray(netcdf_file, 'bph', bp(0:nr+1,0:ns+1,0:np))

    ! Convert from B_theta to B_s:
    bs(0:nr+1,0:ns,0:np+1) = -bs(0:nr+1,0:ns,0:np+1)

    IF (control) PRINT*,'Read magnetic field from file '//netcdf_file

END SUBROUTINE read_magnetic_field

!===============================================================================
SUBROUTINE compute_b_stokes(alr, als, alp)
! Compute bsr, bss, bsp at cell centres from alr, als, alp using Stokes theorem.
    REAL(d), INTENT(IN):: alr(1:nr,0:ns,0:np)   ! A_rho*L_rho on (interior) edges
    REAL(d), INTENT(IN):: als(0:nr,1:ns,0:np)   ! A_s*L_s on (interior) edges
    REAL(d), INTENT(IN):: alp(0:nr,0:ns,1:np)   ! A_phi*L_phi on (interior) edges

    IF (.NOT. ALLOCATED(br)) ALLOCATE(br(0:nr,0:ns+1,0:np+1))
    IF (.NOT. ALLOCATED(bs)) ALLOCATE(bs(0:nr+1,0:ns,0:np+1))
    IF (.NOT. ALLOCATED(bp)) ALLOCATE(bp(0:nr+1,0:ns+1,0:np))
    IF (.NOT. ALLOCATED(bgr)) ALLOCATE(bgr(0:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(bgs)) ALLOCATE(bgs(0:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(bgp)) ALLOCATE(bgp(0:nr,0:ns,0:np))

    ! Compute interior values (with face area weighting included):
    br(0:nr,1:ns,1:np) = als(0:nr,1:ns,1:np) - als(0:nr,1:ns,0:np-1) &
        + alp(0:nr,0:ns-1,1:np) - alp(0:nr,1:ns,1:np)

    bs(1:nr,0:ns,1:np) = alp(1:nr,0:ns,1:np) - alp(0:nr-1,0:ns,1:np) &
        + alr(1:nr,0:ns,0:np-1) - alr(1:nr,0:ns,1:np)

    bp(1:nr,1:ns,0:np) = alr(1:nr,1:ns,0:np) - alr(1:nr,0:ns-1,0:np) &
        + als(0:nr-1,1:ns,0:np) - als(1:nr,1:ns,0:np)

    ! Compute ghost values:
    CALL bs_ghost_values(br, bs, bp)

    ! Compute B at grid points:
    CALL compute_bg_average(br, bs, bp)

    ! Remove face area weighting:
    br(0:nr,0:ns+1,0:np+1) = br(0:nr,0:ns+1,0:np+1) / Sbr(0:nr,0:ns+1,0:np+1)
    bs(0:nr+1,0:ns,0:np+1) = bs(0:nr+1,0:ns,0:np+1) / Sbs(0:nr+1,0:ns,0:np+1)
    bp(0:nr+1,0:ns+1,0:np) = bp(0:nr+1,0:ns+1,0:np) / Sbp(0:nr+1,0:ns+1,0:np)

END SUBROUTINE compute_b_stokes

!===============================================================================
SUBROUTINE bs_ghost_values(bsr, bss, bsp)
! Compute ghost values of bsr, bss, bsp using combination of MPI transfer and
! global boundary conditions. Order matters.
    REAL(d), INTENT(INOUT):: bsr(0:nr,0:ns+1,0:np+1)   ! B_rho*S_rho on faces
    REAL(d), INTENT(INOUT):: bss(0:nr+1,0:ns,0:np+1)   ! B_s*S_s on faces
    REAL(d), INTENT(INOUT):: bsp(0:nr+1,0:ns+1,0:np)   ! B_phi*S_phi on faces

    ! MPI transfers in p-direction:
    ! - bsr up
    CALL MPI_SENDRECV(bsr(0:nr,0:ns+1,np), (nr+1) * (ns+2), mpireal, &
        mpi_p_nxt, tag, bsr(0:nr,0:ns+1,0), (nr+1) * (ns+2), mpireal, &
        mpi_p_prv, tag, comm, stat, err)
    ! - bss up
    CALL MPI_SENDRECV(bss(0:nr+1,0:ns,np), (nr+2) * (ns+1), mpireal, &
        mpi_p_nxt, tag, bss(0:nr+1,0:ns,0), (nr+2) * (ns+1), mpireal, &
        mpi_p_prv, tag, comm, stat, err)
    ! - bsr down
    CALL MPI_SENDRECV(bsr(0:nr,0:ns+1,1), (nr+1) * (ns+2), mpireal, &
        mpi_p_prv, tag, bsr(0:nr,0:ns+1,np+1), (nr+1) * (ns+2), mpireal, &
        mpi_p_nxt, tag, comm, stat, err)
    ! - bss down
    CALL MPI_SENDRECV(bss(0:nr+1,0:ns,1), (nr+2) * (ns+1), mpireal, &
        mpi_p_prv, tag, bss(0:nr+1,0:ns,np+1), (nr+2) * (ns+1), mpireal, &
        mpi_p_nxt, tag, comm, stat, err)

    ! MPI transfers in s-direction:
    ! - bsr up
    IF (.NOT. north) &
        CALL MPI_SEND(bsr(0:nr,ns,0:np+1), (nr+1) * (np+2), mpireal, &
            mpi_s_nxt, tag, comm, err)
    IF (.NOT. south) &
        CALL MPI_RECV(bsr(0:nr,0,0:np+1), (nr+1) * (np+2), mpireal, &
            mpi_s_prv, tag, comm, stat, err)
    ! - bsp up
    IF (.NOT. north) &
        CALL MPI_SEND(bsp(0:nr+1,ns,0:np), (nr+2) * (np+1), mpireal, &
            mpi_s_nxt, tag, comm, err)
    IF (.NOT. south) &
        CALL MPI_RECV(bsp(0:nr+1,0,0:np), (nr+2) * (np+1), mpireal, &
            mpi_s_prv, tag, comm, stat, err)
    ! - bsr down
    IF (.NOT. south) &
        CALL MPI_SEND(bsr(0:nr,1,0:np+1), (nr+1) * (np+2), mpireal, &
            mpi_s_prv, tag, comm, err)
    IF (.NOT. north) &
        CALL MPI_RECV(bsr(0:nr,ns+1,0:np+1), (nr+1) * (np+2), mpireal, &
            mpi_s_nxt, tag, comm, stat, err)
    ! - bsp down
    IF (.NOT. south) &
        CALL MPI_SEND(bsp(0:nr+1,1,0:np), (nr+2) * (np+1), mpireal, &
            mpi_s_prv, tag, comm, err)
    IF (.NOT. north) &
        CALL MPI_RECV(bsp(0:nr+1,ns+1,0:np), (nr+2) * (np+1), mpireal, &
            mpi_s_nxt, tag, comm, stat, err)

    ! MPI transfers in r-direction:
    ! - bss up
    IF (.NOT. outer) &
        CALL MPI_SEND(bss(nr,0:ns,0:np+1), (ns+1) * (np+2), mpireal, &
            mpi_r_nxt, tag, comm, err)
    IF (.NOT. surface) &
        CALL MPI_RECV(bss(0,0:ns,0:np+1), (ns+1) * (np+2), mpireal, &
            mpi_r_prv, tag, comm, stat, err)
    ! - bsp up
    IF (.NOT. outer) &
        CALL MPI_SEND(bsp(nr,0:ns+1,0:np), (ns+2) * (np+1), mpireal, &
            mpi_r_nxt, tag, comm, err)
    IF (.NOT. surface) &
        CALL MPI_RECV(bsp(0,0:ns+1,0:np), (ns+2) * (np+1), mpireal, &
            mpi_r_prv, tag, comm, stat, err)
    ! - bss down
    IF (.NOT. surface) &
        CALL MPI_SEND(bss(1,0:ns,0:np+1), (ns+1) * (np+2), mpireal, &
            mpi_r_prv, tag, comm, err)
    IF (.NOT. outer) &
        CALL MPI_RECV(bss(nr+1,0:ns,0:np+1), (ns+1) * (np+2), mpireal, &
            mpi_r_nxt, tag, comm, stat, err)
    ! - bsp down
    IF (.NOT. surface) &
        CALL MPI_SEND(bsp(1,0:ns+1,0:np), (ns+2) * (np+1), mpireal, &
                mpi_r_prv, tag, comm, err)
    IF (.NOT. outer) &
        CALL mpi_recv(bsp(nr+1,0:ns+1,0:np), (ns+2) * (np+1), mpireal, &
            mpi_r_nxt, tag, comm, stat, err)
        
    ! BCs at inner boundary (js = jp = 0)
    IF (surface) THEN
        bss(0,0:ns,0:np+1) = Sbs(0,0:ns,0:np+1) / dns(0,0:ns,0:np+1) * ( &
            bss(1,0:ns,0:np+1) * dns(1,0:ns,0:np+1) / Sbs(1,0:ns,0:np+1) &
            + bsr(0,0:ns,0:np+1) * dnr(0,0:ns,0:np+1) / Sbr(0,0:ns,0:np+1) &
            - bsr(0,1:ns+1,0:np+1) * dnr(0,1:ns+1,0:np+1) &
            / Sbr(0,1:ns+1,0:np+1))
        bsp(0,0:ns+1,0:np) = Sbp(0,0:ns+1,0:np) / dnp(0,0:ns+1,0:np) * ( &
            bsp(1,0:ns+1,0:np) * dnp(1,0:ns+1,0:np) / Sbp(1,0:ns+1,0:np) &
            + bsr(0,0:ns+1,0:np) * dnr(0,0:ns+1,0:np) / Sbr(0,0:ns+1,0:np) &
            - bsr(0,0:ns+1,1:np+1) * dnr(0,0:ns+1,1:np+1) &
            / Sbr(0,0:ns+1,1:np+1))
    END IF
    
    ! BCs at outer boundary (js = jp = 0)
    IF (outer) THEN
        bss(nr+1,0:ns,0:np+1) = Sbs(nr+1,0:ns,0:np+1) / &
            dns(nr+1,0:ns,0:np+1) * ( &
            bss(nr,0:ns,0:np+1) * dns(nr,0:ns,0:np+1) / Sbs(nr,0:ns,0:np+1) &
            - bsr(nr,0:ns,0:np+1) * dnr(nr,0:ns,0:np+1) &
            / Sbr(nr,0:ns,0:np+1) &
            + bsr(nr,1:ns+1,0:np+1) * dnr(nr,1:ns+1,0:np+1) &
            / Sbr(nr,1:ns+1,0:np+1))
        bsp(nr+1,0:ns+1,0:np) = Sbp(nr+1,0:ns+1,0:np) &
            / dnp(nr+1,0:ns+1,0:np) * ( &
            bsp(nr,0:ns+1,0:np) * dnp(nr,0:ns+1,0:np) / Sbp(nr,0:ns+1,0:np) &
            - bsr(nr,0:ns+1,0:np) * dnr(nr,0:ns+1,0:np) &
            / Sbr(nr,0:ns+1,0:np) &
            + bsr(nr,0:ns+1,1:np+1) * dnr(nr,0:ns+1,1:np+1) &
            / Sbr(nr,0:ns+1,1:np+1))
    END IF

    ! Polar BCs
    IF (south) THEN
        CALL MPI_SENDRECV(bsr(0:nr,1,0:np+1), (nr+1) * (np+2), mpireal, &
            mpi_opp, tag, bsr(0:nr,0,0:np+1), (nr+1) * (np+2), mpireal, &
            mpi_opp, tag, comm, stat, err)
        CALL MPI_SENDRECV(bss(0:nr+1,1,0:np+1), (nr+2) * (np+2), mpireal, &
            mpi_opp, tag, bss(0:nr+1,0,0:np+1), (nr+2) * (np+2), mpireal, &
            mpi_opp, tag, comm, stat,err)
        bss(0:nr+1,0,0:np+1) = 0.5_d*(bss(0:nr+1,1,0:np+1) &
            - bss(0:nr+1,0,0:np+1))
        CALL MPI_SENDRECV(bsp(0:nr+1,1,0:np), (nr+2) * (np+1), mpireal, &
            mpi_opp, tag, bsp(0:nr+1,0,0:np), (nr+2) * (np+1), mpireal, &
            mpi_opp, tag, comm, stat,err)
        bsp(0:nr+1,0,0:np) = -bsp(0:nr+1,0,0:np)
    END IF

    IF (north) THEN
        CALL MPI_SENDRECV(bsr(0:nr,ns,0:np+1), (nr+1) * (np+2), mpireal, &
            mpi_opp, tag, bsr(0:nr,ns+1,0:np+1), (nr+1) * (np+2), mpireal, &
            mpi_opp, tag, comm, stat,err)
        CALL MPI_SENDRECV(bss(0:nr+1,ns-1,0:np+1), (nr+2) * (np+2), mpireal, &
            mpi_opp, tag, bss(0:nr+1,ns,0:np+1), (nr+2) * (np+2), mpireal, &
            mpi_opp, tag, comm,stat,err)
        bss(0:nr+1,ns,0:np+1) = 0.5_d*(bss(0:nr+1,ns-1,0:np+1) &
            - bss(0:nr+1,ns,0:np+1))
        CALL MPI_SENDRECV(bsp(0:nr+1,ns,0:np), (nr+2) * (np+1), mpireal, &
            mpi_opp, tag, bsp(0:nr+1,ns+1,0:np), (nr+2) * (np+1), mpireal, &
            mpi_opp, tag, comm, stat, err)
        bsp(0:nr+1,ns+1,0:np) = -bsp(0:nr+1,ns+1,0:np)
    END IF

END SUBROUTINE bs_ghost_values

!===============================================================================
SUBROUTINE compute_bg_average(bsr, bss, bsp)
! Compute bgr, bgs, bgp by area-weighted average of bsr, bss, bsp to grid pts.
! Also computes bbg, including ghost values in rho if hyperdiffusion is in use.
    REAL(d), INTENT(IN):: bsr(0:nr,0:ns+1,0:np+1)   ! B_rho*S_rho on faces
    REAL(d), INTENT(IN):: bss(0:nr+1,0:ns,0:np+1)   ! B_s*S_s on faces
    REAL(d), INTENT(IN):: bsp(0:nr+1,0:ns+1,0:np)   ! B_phi*S_phi on faces

    IF (.NOT. ALLOCATED(bgr)) ALLOCATE(bgr(0:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(bgs)) ALLOCATE(bgs(0:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(bgp)) ALLOCATE(bgp(0:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(bbg)) ALLOCATE(bbg(0:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(bbgs)) ALLOCATE(bbgs(-1:nr+1,0:ns,0:np))

    bgr(0:nr,0:ns,0:np) = (bsr(0:nr,0:ns,0:np) + bsr(0:nr,1:ns+1,0:np) &
            + bsr(0:nr,0:ns,1:np+1) + bsr(0:nr,1:ns+1,1:np+1)) &
            / Sbrg(0:nr,0:ns,0:np)
    bgs(0:nr,0:ns,0:np) = (bss(0:nr,0:ns,0:np) + bss(1:nr+1,0:ns,0:np) &
            + bss(0:nr,0:ns,1:np+1) + bss(1:nr+1,0:ns,1:np+1)) &
            / Sbsg(0:nr,0:ns,0:np)
    bgp(0:nr,0:ns,0:np) = (bsp(0:nr,0:ns,0:np) + bsp(0:nr,1:ns+1,0:np) &
            + bsp(1:nr+1,0:ns,0:np) + bsp(1:nr+1,1:ns+1,0:np)) &
            / Sbpg(0:nr,0:ns,0:np)

    bbg(0:nr,0:ns,0:np) = bgr(0:nr,0:ns,0:np) * bgr(0:nr,0:ns,0:np) &
        + bgs(0:nr,0:ns,0:np) * bgs(0:nr,0:ns,0:np) &
        + bgp(0:nr,0:ns,0:np) * bgp(0:nr,0:ns,0:np)
        
    bbgs(0:nr,0:ns,0:np) = soften(bbg(0:nr,0:ns,0:np))

    IF (hyper_flag) THEN
        IF (.NOT. outer) &
            CALL MPI_SEND(bbgs(nr-1,0:ns,0:np), (ns+1) * (np+1), mpireal, &
                mpi_r_nxt, tag, comm, err)
        IF (.NOT.surface) &
            CALL MPI_RECV(bbgs(-1,0:ns,0:np), (ns+1) * (np+1), mpireal, &
                mpi_r_prv, tag, comm, stat, err)
        IF (.NOT.surface) &
            CALL MPI_SEND(bbgs(1,0:ns,0:np), (ns+1) * (np+1), mpireal, &
                mpi_r_prv, tag, comm, err)
        IF (.NOT.outer) &
            CALL MPI_RECV(bbgs(nr+1,0:ns,0:np) , (ns+1) * (np+1), mpireal, &
                mpi_r_nxt, tag, comm, stat, err)
        IF (outer) bbgs(nr+1,0:ns,0:np) = bbgs(nr,0:ns,0:np)
        IF (surface) bbgs(-1,0:ns,0:np) = bbgs(0,0:ns,0:np)
    END IF

END SUBROUTINE compute_bg_average

!===============================================================================
FUNCTION soften(bbg)
! Apply softening function to limit size of bbg.
    REAL(d), INTENT(IN):: bbg(:,:,:)
    REAL(d):: soften(SIZE(bbg, 1),SIZE(bbg, 2),SIZE(bbg, 3))
    REAL(d), PARAMETER:: eps = 1e-4_d
    
    soften = bbg
    WHERE(soften < eps) soften = eps

END FUNCTION soften

!===============================================================================
SUBROUTINE test_compute_b_stokes(alr, als, alp, test_case, output_netcdf)
! Test the computation of B on cell faces against exact B from b_exact routine.
! The vector potential components are supplied as input.
    REAL(d), INTENT(IN):: alr(1:nr,0:ns,0:np)   ! A_rho*L_rho on (interior) edges
    REAL(d), INTENT(IN):: als(0:nr,1:ns,0:np)   ! A_s*L_s on (interior) edges
    REAL(d), INTENT(IN):: alp(0:nr,0:ns,1:np)   ! A_phi*L_phi on (interior) edges
    INTEGER, INTENT(IN):: test_case
    LOGICAL, INTENT(IN):: output_netcdf

    ! Local variables:
    REAL(d), ALLOCATABLE:: error(:,:,:)
    REAL(d):: local_max, error_br, error_bs, error_bp
    REAL(d):: br0, bs0, bp0
    INTEGER:: i, j, k

    CALL compute_b_stokes(alr, als, alp)
    
    IF (output_netcdf) THEN
        ! Output B on cell faces as 3D arrays:
        CALL ncdf_init_snapshot_file('test_b.nc', .TRUE.)
        CALL ncdf_write_3darray('test_b.nc', 'br', 'r', 'thc', 'phc', &
            br(0:nr,0:ns+1,0:np+1))
        CALL ncdf_write_3darray('test_b.nc', 'bth', 'rc', 'th', 'phc', &
            -bs(0:nr+1,0:ns,0:np+1))
        CALL ncdf_write_3darray('test_b.nc', 'bph', 'rc', 'thc', 'ph', &
            bp(0:nr+1,0:ns+1,0:np))
    END IF
    
    ! Maximum error in br:
    ALLOCATE(error(0:nr,1:ns,1:np))
    error(0:nr,1:ns,1:np) = 0.0_d
    DO i = 0, nr
        DO j = 1, ns
            DO k = 1, np
                CALL b_exact(rg(i), sc(j), pc(k), br0, bs0, bp0, test_case)
                error(i,j,k) = br(i,j,k) - br0
            END DO
        END DO
    END DO
    local_max = MAXVAL(ABS(error))
    CALL MPI_ALLREDUCE(local_max, error_br, 1, mpireal, MPI_MAX, &
       comm, err)
       
    ! Maximum error in bs:
    DEALLOCATE(error)
    ALLOCATE(error(1:nr,0:ns,1:np))
    error(1:nr,0:ns,1:np) = 0.0_d
    DO i = 1, nr
        DO j = 0, ns
            DO k = 1, np
                CALL b_exact(rc(i), sg(j), pc(k), br0, bs0, bp0, test_case)
                error(i,j,k) = bs(i,j,k) - bs0
            END DO
        END DO
    END DO
    local_max = MAXVAL(ABS(error))
    CALL MPI_ALLREDUCE(local_max, error_bs, 1, mpireal, MPI_MAX, &
       comm, err)

    ! Maximum error in bp:
    DEALLOCATE(error)
    ALLOCATE(error(1:nr,1:ns,0:np))
    error(1:nr,1:ns,0:np) = 0.0_d
    DO i = 1, nr
        DO j = 1, ns
            DO k = 0, np
                CALL b_exact(rc(i), sc(j), pg(k), br0, bs0, bp0, test_case)
                error(i,j,k) = bp(i,j,k) - bp0
            END DO
        END DO
    END DO
    local_max = MAXVAL(ABS(error))
    CALL MPI_ALLREDUCE(local_max, error_bp, 1, mpireal, MPI_MAX, &
       comm, err)
   
    IF (control) THEN
        PRINT*, 'Max |error| in Br: ', error_br
        PRINT*, 'Max |error| in Bs: ', error_bs
        PRINT*, 'Max |error| in Bp: ', error_bp
    END IF

END SUBROUTINE test_compute_b_stokes

!===============================================================================
SUBROUTINE b_exact(r, s, p, br0, bs0, bp0, test_case)
! Return exact solution for test field, br, bs, bp components.
! Integer test_case specifies which exact solution.
    REAL(d), INTENT(IN):: r, s, p
    REAL(d), INTENT(OUT):: br0, bs0, bp0
    INTEGER, INTENT(IN):: test_case

    ! Local variables:
    REAL(d):: rad

    rad = EXP(r)
    
    IF (test_case == 0) THEN    ! axisymmetric PFSS dipole
        br0 = 2.0_d * (1 + 2 * rss**3/rad**3) / (1 + 2 * rss**3) * s
        bs0 = 2.0_d * (1 - rss**3/rad**3) / (1 + 2 * rss**3) * &
            SQRT(1.0_d - s**2)
        bp0 = 0.0_d
    END IF

END SUBROUTINE b_exact

!===============================================================================
SUBROUTINE set_b_exact(test_case)
! Set br, bs, bp to analytical test field.
    INTEGER, INTENT(IN):: test_case

    ! Local variables:
    INTEGER:: i, j, k
    REAL(d):: sintheta, rad, cos2theta, br0, bs0, bp0

    IF (.NOT. ALLOCATED(br)) ALLOCATE(br(0:nr,0:ns+1,0:np+1))
    IF (.NOT. ALLOCATED(bs)) ALLOCATE(bs(0:nr+1,0:ns,0:np+1))
    IF (.NOT. ALLOCATED(bp)) ALLOCATE(bp(0:nr+1,0:ns+1,0:np))
    
    br(0:nr,0:ns+1,0:np+1) = 0.0_d
    bs(0:nr+1,0:ns,0:np+1) = 0.0_d
    bp(0:nr+1,0:ns+1,0:np) = 0.0_d

    DO k = 1, np
        DO j = 1, ns
            DO i = 0, nr
                CALL b_exact(rg(i), sc(j), pc(k), br0, bs0, bp0, test_case)
                br(i,j,k) = br0
            END DO
        END DO
    END DO

    DO k = 1, np
        DO j = 0, ns
            DO i = 1, nr
                CALL b_exact(rc(i), sg(j), pc(k), br0, bs0, bp0, test_case)
                bs(i,j,k) = bs0
            END DO
        END DO
    END DO

    DO k = 1, np
        DO j = 1, ns
            DO i = 0, nr
                CALL b_exact(rc(i), sc(j), pg(k), br0, bs0, bp0, test_case)
                bp(i,j,k) = bp0
            END DO
        END DO
    END DO

    ! Add face area weightings and compute ghost values:
    br(0:nr,0:ns+1,0:np+1) = br(0:nr,0:ns+1,0:np+1) * Sbr(0:nr,0:ns+1,0:np+1)
    bs(0:nr+1,0:ns,0:np+1) = bs(0:nr+1,0:ns,0:np+1) * Sbs(0:nr+1,0:ns,0:np+1)
    bp(0:nr+1,0:ns+1,0:np) = bp(0:nr+1,0:ns+1,0:np) * Sbp(0:nr+1,0:ns+1,0:np)
    CALL bs_ghost_values(br, bs, bp)
    
    ! Average to grid points:
    CALL compute_bg_average(br, bs, bp)

    ! Remove face area weightings:
    br(0:nr,0:ns+1,0:np+1) = br(0:nr,0:ns+1,0:np+1) / Sbr(0:nr,0:ns+1,0:np+1)
    bs(0:nr+1,0:ns,0:np+1) = bs(0:nr+1,0:ns,0:np+1) / Sbs(0:nr+1,0:ns,0:np+1)
    bp(0:nr+1,0:ns+1,0:np) = bp(0:nr+1,0:ns+1,0:np) / Sbp(0:nr+1,0:ns+1,0:np)

END SUBROUTINE set_b_exact

!*******************************************************************************
END MODULE bfield
!*******************************************************************************
