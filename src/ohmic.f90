!*******************************************************************************
MODULE ohmic
!*******************************************************************************
! This module implements ohmic diffusion.
!*******************************************************************************
    USE params
    USE mpitools, ONLY: control, mpireal, MPI_MIN, comm, err
    USE grid, ONLY: nr, ns, np, dlp
    USE friction, ONLY: nu0

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: ohmic_flag, ohmic_ratio, ohmic_eta, ohmic_enhanced, dt_max_ohmic
    PUBLIC:: initialize_ohmic

!*******************************************************************************
    ! Specified in user namelist:
    REAL(d):: ohmic_ratio   ! ratio eta/(R_sun**2 * nu0)  [dimensionless]
    REAL(d):: ohmic_enhanced    ! coeff of enhanced diffusion [dimensionless]

    ! Derived:
    LOGICAL:: ohmic_flag     ! flag for whether ohmic diffusion is turned on
    REAL(d):: ohmic_eta     ! coefficient of ohmic diffusivity [R_sun**2/s]
    REAL(d):: dt_max_ohmic=1e10_d   ! Maximum timestep for ohmic diffusion [s]

!*******************************************************************************
CONTAINS
!*******************************************************************************

!===============================================================================
SUBROUTINE initialize_ohmic()
! Set ohmic diffusion coefficient and calculate maximum timestep.
    
    ! Local variables:
    REAL(d):: local_dt, global_dt

    ! Set ohmic_eta coefficient:
    IF (ohmic_ratio < 1.0d-20) THEN
       ohmic_flag = .FALSE.
       ohmic_eta = 0.0_d
    ELSE
       ohmic_flag = .TRUE.
       ohmic_eta = ohmic_ratio * nu0
    END IF

    ! Calculate maximum timestep:
    IF (ohmic_flag) THEN
        local_dt = MINVAL(ABS(dlp(1:nr,1:ns-1,0) ** 2 / &
            (ohmic_eta)))
        CALL MPI_ALLREDUCE(local_dt, global_dt, 1, mpireal, MPI_MIN, comm, err)
        dt_max_ohmic = cflfact * global_dt
        IF (control) PRINT*,'Ohmic diffusion dt_max [s] =', dt_max_ohmic
    ELSE
        dt_max_ohmic = 1e10_d
    END IF

END SUBROUTINE initialize_ohmic

!*******************************************************************************
END MODULE ohmic
!*******************************************************************************
