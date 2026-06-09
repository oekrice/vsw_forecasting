!*******************************************************************************
MODULE outflow
!*******************************************************************************
! This module implements the radial outflow velocity.
! (Note: surface flows are dealt with in the driver.f90 module.)
!*******************************************************************************
    USE params
    USE mpitools, ONLY: control, surface, mpireal, comm, err, MPI_MIN
    USE grid, ONLY: nr, ns, np, dnr
    USE ncio, ONLY: ncdf_read_1darray_r

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: read_outflow, dt_max_outflow, vr

!*******************************************************************************
    REAL(d), ALLOCATABLE:: vr(:)    ! Outflow speed at grid pts (0:nr) [units?]
    REAL(d):: dt_max_outflow=1e10_d   ! Maximum timestep for outflow [s]

!*******************************************************************************
CONTAINS
!*******************************************************************************

!===============================================================================
SUBROUTINE read_outflow(netcdf_file)
! Read 1d profile of outflow velocity in rho and set 1d array.
    CHARACTER(*), INTENT(IN):: netcdf_file     ! netcdf file containing nu.
    
    ! Local variables:
    INTEGER:: i
    REAL(d):: dt0, local_dt, global_dt
    
    ! Allocate array:
    IF (ALLOCATED(vr)) DEALLOCATE(vr)
    ALLOCATE(vr(0:nr))
    vr(0:nr) = 0.0_d
    
    ! Read array:
    CALL ncdf_read_1darray_r(netcdf_file, 'vrg', vr(0:nr))

    ! Set to zero on inner boundary:
    IF (surface) vr(0) = 0.0_d
    
    ! Determine maximum timestep for outflow:
    local_dt = 1.0d10
    DO i = 0, nr
        IF (vr(i) > 1d-10) THEN
            dt0 = MINVAL(dnr(i,1:ns,1:np)) / vr(i)
            IF (dt0 < local_dt) local_dt = dt0
        END IF
    END DO
    CALL MPI_ALLREDUCE(local_dt, global_dt, 1, mpireal, MPI_MIN, comm, err)

    dt_max_outflow = cflfact * global_dt
    
    IF (control) PRINT*,'Read outflow profile from file '//netcdf_file// &
        ' --- dt_max [s] =', dt_max_outflow

END SUBROUTINE read_outflow

!*******************************************************************************
END MODULE outflow
!*******************************************************************************
