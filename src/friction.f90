!*******************************************************************************
MODULE friction
!*******************************************************************************
! This module implements the friction.
!*******************************************************************************
    USE params
    USE mpitools, ONLY: control, surface, mpireal, comm, err, MPI_MIN
    USE grid, ONLY: nr, ns, np, dlp
    USE ncio, ONLY: ncdf_read_2darray_rs

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: nu0, nu, read_friction, dt_max_fric

!*******************************************************************************
    REAL(d):: nu0   ! coefficient of friction [units?]
    REAL(d), ALLOCATABLE:: nu(:,:,:)   ! Friction coeff at grid pts
                                       ! (0:nr,0:ns,0:np)  [units?]
    REAL(d):: dt_max_fric=1e10_d   ! Maximum timestep for friction [s]

!*******************************************************************************
CONTAINS
!*******************************************************************************

!===============================================================================
SUBROUTINE read_friction(netcdf_file)
! Read 2d profile of friction coefficient in (rho, s) and set 3d array.
    CHARACTER(*), INTENT(IN):: netcdf_file     ! netcdf file containing nu.

    ! Local variables:
    INTEGER:: i
    REAL(d):: local_dt, global_dt
        
    ! Allocate array:
    IF (ALLOCATED(nu)) DEALLOCATE(nu)
    ALLOCATE(nu(0:nr,0:ns,0:np))
    nu(0:nr,0:ns,0:np) = 0.0_d

    ! Read array:
    CALL ncdf_read_2darray_rs(netcdf_file, 'nug', nu(0:nr,0:ns,0))

    ! Include strength parameter:
    nu = nu*nu0
    
    ! Copy to other points in p-direction:
    DO i = 1, np
        nu(:,:,i) = nu(:,:,0)
    END DO

    ! Determine maximum timestep for diffusive friction term:
    local_dt = MINVAL(ABS(dlp(1:nr,1:ns-1,0) ** 2 / &
        (1e-10_d + nu(1:nr,1:ns-1,0))))
    CALL MPI_ALLREDUCE(local_dt, global_dt, 1, mpireal, MPI_MIN, comm, err)
    dt_max_fric = cflfact * global_dt
    
    IF (control) PRINT*,'Read friction profile from file '//netcdf_file// &
        ' --- dt_max [s] =', dt_max_fric

END SUBROUTINE read_friction

!*******************************************************************************
END MODULE friction
!*******************************************************************************
