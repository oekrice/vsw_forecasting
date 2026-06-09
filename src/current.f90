!*******************************************************************************
MODULE current
!*******************************************************************************
! This module implements the electric current density.
!*******************************************************************************
    USE params
    USE mpitools, ONLY: control, north, south, mpireal, mpi_opp, comm, stat, &
        err, MPI_MAX
    USE grid, ONLY: nr, ns, np, Sjr, Sjs, Sjp, dnr, dns, dnp, rc, sc, pc, &
        rg, sg, pg
    USE ncio, ONLY: ncdf_init_snapshot_file, ncdf_write_3darray

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: jr, js, jp, jgr, jgs, jgp, compute_j_stokes, test_compute_j_stokes

!*******************************************************************************
    REAL(d), ALLOCATABLE:: jr(:,:,:)   ! J_rho on staggered faces
                                       ! (0:nr+1,0:ns,0:np)
    REAL(d), ALLOCATABLE:: js(:,:,:)   ! J_s on staggered faces
                                       ! (0:nr,0:ns+1,0:np)
    REAL(d), ALLOCATABLE:: jp(:,:,:)   ! J_phi on staggered faces
                                       ! (0:nr,0:ns,0:np+1)
    REAL(d), ALLOCATABLE:: jgr(:,:,:)   ! J_rho at grid pts (0:nr,0:ns,0:np)
    REAL(d), ALLOCATABLE:: jgs(:,:,:)   ! J_s at grid pts (0:nr,0:ns,0:np)
    REAL(d), ALLOCATABLE:: jgp(:,:,:)   ! J_phi at grid pts (0:nr,0:ns,0:np)
    INTEGER, PARAMETER:: tag=1

!*******************************************************************************
CONTAINS
!*******************************************************************************

!===============================================================================
SUBROUTINE compute_j_stokes(br, bs, bp)
! Compute jr, js, jp on staggered faces from br, bs, bp at faces, using
! Stokes theorem.
    REAL(d), INTENT(IN):: br(0:nr,0:ns+1,0:np+1)   ! B_rho on faces
    REAL(d), INTENT(IN):: bs(0:nr+1,0:ns,0:np+1)   ! B_s on faces
    REAL(d), INTENT(IN):: bp(0:nr+1,0:ns+1,0:np)   ! B_phi on faces
       
    IF (.NOT. ALLOCATED(jr)) ALLOCATE(jr(0:nr+1,0:ns,0:np))
    IF (.NOT. ALLOCATED(js)) ALLOCATE(js(0:nr,0:ns+1,0:np))
    IF (.NOT. ALLOCATED(jp)) ALLOCATE(jp(0:nr,0:ns,0:np+1))

    ! Compute j*Area on staggered faces:
    jr(0:nr+1,0:ns,0:np) = (bs(0:nr+1,0:ns,1:np+1) * dns(0:nr+1,0:ns,1:np+1)) &
        - (bs(0:nr+1,0:ns,0:np) * dns(0:nr+1,0:ns,0:np)) &
        + (bp(0:nr+1,0:ns,0:np) * dnp(0:nr+1,0:ns,0:np)) &
        - (bp(0:nr+1,1:ns+1,0:np) * dnp(0:nr+1,1:ns+1,0:np))

    js(0:nr,0:ns+1,0:np) = (bp(1:nr+1,0:ns+1,0:np) * dnp(1:nr+1,0:ns+1,0:np)) &
        - (bp(0:nr,0:ns+1,0:np) * dnp(0:nr,0:ns+1,0:np)) &
        + (br(0:nr,0:ns+1,0:np) * dnr(0:nr,0:ns+1,0:np)) &
        - (br(0:nr,0:ns+1,1:np+1) * dnr(0:nr,0:ns+1,1:np+1))

    jp(0:nr,0:ns,0:np+1) = (br(0:nr,1:ns+1,0:np+1) * dnr(0:nr,1:ns+1,0:np+1)) &
        - (br(0:nr,0:ns,0:np+1) * dnr(0:nr,0:ns,0:np+1)) &
        + (bs(0:nr,0:ns,0:np+1) * dns(0:nr,0:ns,0:np+1)) &
        - (bs(1:nr+1,0:ns,0:np+1) * dns(1:nr+1,0:ns,0:np+1))

    ! Get values on polar boundaries by boundary condition:
    CALL j_polar_boundary

    ! Remove area weighting:
    jr(0:nr+1,0:ns,0:np) = jr(0:nr+1,0:ns,0:np) / Sjr(0:nr+1,0:ns,0:np)
    js(0:nr,0:ns+1,0:np) = js(0:nr,0:ns+1,0:np) / Sjs(0:nr,0:ns+1,0:np)
    jp(0:nr,0:ns,0:np+1) = jp(0:nr,0:ns,0:np+1) / Sjp(0:nr,0:ns,0:np+1)

    ! Average to get values at grid points:
    CALL compute_jg_average

END SUBROUTINE compute_j_stokes

!===============================================================================
SUBROUTINE j_polar_boundary()
! Set values of jr, jp on polar boundary, using neighbouring values.

    IF (south) THEN
        CALL MPI_SENDRECV(jr(0:nr+1,1,0:np), (nr+2) * (np+1), mpireal, &
            mpi_opp, tag, jr(0:nr+1,0,0:np), (nr+2) * (np+1), mpireal, &
            mpi_opp, tag, comm, stat, err)
        jr(0:nr+1,0,0:np) = 0.5_d * (jr(0:nr+1,1,0:np) + jr(0:nr+1,0,0:np))
        CALL MPI_SENDRECV(jp(0:nr,1,0:np+1), (nr+1) * (np+2), mpireal, &
            mpi_opp, tag, jp(0:nr,0,0:np+1), (nr+1) * (np+2), mpireal, &
            mpi_opp, tag, comm, stat, err)
        jp(0:nr,0,0:np+1) = 0.5_d * (jp(0:nr,1,0:np+1) - jp(0:nr,0,0:np+1))
    END IF
    
    IF (north) THEN
        CALL MPI_SENDRECV(jr(0:nr+1,ns-1,0:np), (nr+2) * (np+1), mpireal, &
            mpi_opp, tag, jr(0:nr+1,ns,0:np), (nr+2) * (np+1), mpireal, &
            mpi_opp, tag, comm, stat, err)
        jr(0:nr+1,ns,0:np) = 0.5_d * (jr(0:nr+1,ns-1,0:np) + jr(0:nr+1,ns,0:np))
        CALL MPI_SENDRECV(jp(0:nr,ns-1,0:np+1), (nr+1) * (np+2), mpireal, &
            mpi_opp, tag, jp(0:nr,ns,0:np+1), (nr+1) * (np+2), mpireal, &
            mpi_opp, tag, comm, stat, err)
        jp(0:nr,ns,0:np+1) = 0.5_d * (jp(0:nr,ns-1,0:np+1) - jp(0:nr,ns,0:np+1))
    END IF

END SUBROUTINE j_polar_boundary

!===============================================================================
SUBROUTINE compute_jg_average()
! Average jr, js, jp from staggered face centres to grid pts of unstaggered
! grid.

    IF (.NOT. ALLOCATED(jgr)) ALLOCATE(jgr(0:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(jgs)) ALLOCATE(jgs(0:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(jgp)) ALLOCATE(jgp(0:nr,0:ns,0:np))
    
    jgr(0:nr,0:ns,0:np) = 0.5_d * (jr(0:nr,0:ns,0:np) + jr(1:nr+1,0:ns,0:np))
    jgs(0:nr,0:ns,0:np) = 0.5_d * (js(0:nr,0:ns,0:np) + js(0:nr,1:ns+1,0:np))
    jgp(0:nr,0:ns,0:np) = 0.5_d * (jp(0:nr,0:ns,0:np) + jp(0:nr,0:ns,1:np+1))
    
END SUBROUTINE compute_jg_average

!===============================================================================
SUBROUTINE test_compute_j_stokes(br, bs, bp, test_case, output_netcdf)
! Test the computation of J on edges against exact J from j_exact routine.
! The magnetic field components are supplied as input.
    REAL(d), INTENT(IN):: br(0:nr,0:ns+1,0:np+1)   ! B_rho on faces
    REAL(d), INTENT(IN):: bs(0:nr+1,0:ns,0:np+1)   ! B_s on faces
    REAL(d), INTENT(IN):: bp(0:nr+1,0:ns+1,0:np)   ! B_phi on faces
    INTEGER, INTENT(IN):: test_case
    LOGICAL, INTENT(IN):: output_netcdf

    ! Local variables:
    REAL(d), ALLOCATABLE:: error(:,:,:)
    REAL(d):: local_max, error_jr, error_js, error_jp
    REAL(d):: jr0, js0, jp0
    INTEGER:: i, j, k

    CALL compute_j_stokes(br, bs, bp)

    IF (output_netcdf) THEN
        ! Output B on cell faces as 3D arrays:
        CALL ncdf_write_3darray('test_b.nc', 'jr', 'r', 'th', 'ph', &
            jgr(0:nr,0:ns,0:np))
        CALL ncdf_write_3darray('test_b.nc', 'jth', 'r', 'th', 'ph', &
            -jgs(0:nr,0:ns,0:np))
        CALL ncdf_write_3darray('test_b.nc', 'jph', 'r', 'th', 'ph', &
            jgp(0:nr,0:ns,0:np))
    END IF

    ! Maximum error in jr:
    ALLOCATE(error(1:nr,0:ns,0:np))
    error(1:nr,0:ns,0:np) = 0.0_d
    DO i = 1, nr
        DO j = 0, ns
            DO k = 0, np
                CALL j_exact(rc(i), sg(j), pg(k), jr0, js0, jp0, test_case)
                error(i,j,k) = jr(i,j,k) - jr0
            END DO
        END DO
    END DO
    local_max = MAXVAL(ABS(error))
    CALL MPI_ALLREDUCE(local_max, error_jr, 1, mpireal, MPI_MAX, &
       comm, err)

    ! Maximum error in js:
    DEALLOCATE(error)
    ALLOCATE(error(0:nr,1:ns,0:np))
    error(0:nr,1:ns,0:np) = 0.0_d
    DO i = 0, nr
        DO j = 1, ns
            DO k = 0, np
                CALL j_exact(rg(i), sc(j), pg(k), jr0, js0, jp0, test_case)
                error(i,j,k) = js(i,j,k) - js0
            END DO
        END DO
    END DO
    local_max = MAXVAL(ABS(error))
    CALL MPI_ALLREDUCE(local_max, error_js, 1, mpireal, MPI_MAX, &
       comm, err)

    ! Maximum error in jp:
    DEALLOCATE(error)
    ALLOCATE(error(0:nr,0:ns,1:np))
    error(0:nr,0:ns,1:np) = 0.0_d
    DO i = 0, nr
        DO j = 0, ns
            DO k = 1, np
                CALL j_exact(rg(i), sg(j), pc(k), jr0, js0, jp0, test_case)
                error(i,j,k) = jp(i,j,k) - jp0
            END DO
        END DO
    END DO
    local_max = MAXVAL(ABS(error))
    CALL MPI_ALLREDUCE(local_max, error_jp, 1, mpireal, MPI_MAX, &
       comm, err)

    IF (control) THEN
        PRINT*, 'Max |error| in Jr: ', error_jr
        PRINT*, 'Max |error| in Js: ', error_js
        PRINT*, 'Max |error| in Jp: ', error_jp
    END IF

END SUBROUTINE test_compute_j_stokes

!===============================================================================
SUBROUTINE j_exact(r, s, p, jr0, js0, jp0, test_case)
! Return exact solution for test field, jr, js, jp components.
! Integer test_case specifies which exact solution.
    REAL(d), INTENT(IN):: r, s, p
    REAL(d), INTENT(OUT):: jr0, js0, jp0
    INTEGER, INTENT(IN):: test_case

    ! Local variables:
    REAL(d):: rad

    rad = EXP(r)
    
    IF (test_case == 0) THEN    ! axisymmetric PFSS dipole
        jr0 = 0.0_d
        js0 = 0.0_d
        jp0 = 0.0_d
    END IF

END SUBROUTINE j_exact

!*******************************************************************************
END MODULE current
!*******************************************************************************
