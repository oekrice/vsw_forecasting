!*******************************************************************************
MODULE driver
!*******************************************************************************
! This module implements the boundary driving.
!*******************************************************************************
    USE params
    USE mpitools, ONLY: control, surface, comm, err, MPI_MIN, mpireal
    USE grid, ONLY: nr, ns, np, dlp, dns, dnp
    USE ncio, ONLY: ncdf_read_surface_scalar, ncdf_read_surface_1darray_s, &
        ncdf_read_2darray_sp
    USE time, ONLY: time_now, time_to_str

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: cadence_em, driver_em, driver_v, vbs, vbp, eta_surf
    PUBLIC:: dt_max_driver, read_surface_flows, init_emergence, begin_emergence
    PUBLIC:: als_em, alp_em, t_em

!*******************************************************************************
    ! Specified in user namelist:
    INTEGER:: cadence_em   ! time between boundary driver updates [hours]
    LOGICAL:: driver_v  ! turn on boundary driving velocity?
    LOGICAL:: driver_em  ! turn on active region emergence?
    REAL(d):: eta_surf  ! diffusivity on inner boundary [R_sun**2/s]
    REAL(d), ALLOCATABLE:: vbs0(:)  ! steady v_s(s) on inner boundary p-edges
                                    ! [R_sun/s] (0:ns)
    REAL(d), ALLOCATABLE:: vbp0(:)  ! steady v_p(s) on inner boundary s-edges
                                ! [R_sun/s] (1:ns)
    REAL(d), ALLOCATABLE:: vbs(:,:)  ! total v_s on inner boundary p-edges
                                     ! [R_sun/s] (0:ns,1:np)
    REAL(d), ALLOCATABLE:: vbp(:,:)     ! total v_phi on inner boundary s-edges
                                        ! [R_sun/s] (1:ns,0:np)
    REAL(d), ALLOCATABLE:: als_em(:,:)  ! A_s*L_s on inner boundary s-edges
                                        ! (1:ns,0:np)
    REAL(d), ALLOCATABLE:: alp_em(:,:)  ! A_phi*L_phi on inner boundary p-edges
                                        ! (0:ns,1:np)
    REAL(d), ALLOCATABLE:: als0(:,:)  ! A_s*L_s on inner boundary s-edges
                                        ! (1:ns,0:np) [temporary]
    REAL(d), ALLOCATABLE:: alp0(:,:)  ! A_phi*L_phi on inner boundary p-edges
                                        ! (0:ns,1:np) [temporary]s
    REAL(d):: dt_max_driver=1e10_d
    REAL(d):: t_em  ! cadence_em in seconds

!*******************************************************************************
CONTAINS
!*******************************************************************************

!===============================================================================
SUBROUTINE read_surface_flows(netcdf_file)
! Allocate E*L on cell edges if necessary, and set to zero.
    CHARACTER(*), INTENT(IN):: netcdf_file     ! netcdf file containing flows.

    ! Local variables:
    INTEGER:: i
    REAL(d):: local_dt, global_dt

    IF (.NOT. ALLOCATED(vbs0)) ALLOCATE(vbs0(0:ns))
    IF (.NOT. ALLOCATED(vbp0)) ALLOCATE(vbp0(1:ns))
    IF (.NOT. ALLOCATED(vbs)) ALLOCATE(vbs(0:ns,1:np))
    IF (.NOT. ALLOCATED(vbp)) ALLOCATE(vbp(1:ns,0:np))
    vbs0(0:ns) = 0.0_d
    vbp0(1:ns) = 0.0_d
    vbs(0:ns,1:np) = 0.0_d
    vbp(1:ns,0:np) = 0.0_d
    eta_surf = 0.0_d

    CALL ncdf_read_surface_scalar(netcdf_file, 'etasurf', eta_surf)
    CALL ncdf_read_surface_1darray_s(netcdf_file, 'vs', vbs0(0:ns))
    CALL ncdf_read_surface_1darray_s(netcdf_file, 'vp', vbp0(1:ns))

    IF (.NOT. surface) DEALLOCATE(vbs0, vbp0, vbs, vbp)

    IF (surface) THEN
        DO i = 1, np
            vbs(0:ns,i) = vbs0(0:ns)
        END DO
        DO i = 0, np
            vbp(1:ns,i) = vbp0(1:ns)
        END DO
    END IF

    ! Determine maximum timestep for supergranular diffusion:
    IF (surface) THEN
        local_dt = MINVAL(ABS(dlp(0,1:ns-1,0)**2 / (eta_surf + 1e-8_d)))
    ELSE
        local_dt = 1e10_d
    END IF
    CALL MPI_ALLREDUCE(local_dt, global_dt, 1, mpireal, MPI_MIN, comm, err)
    dt_max_driver = global_dt
    
    ! Determine maximum timestep for vbs0:
    IF (surface) THEN
        local_dt = MINVAL(ABS(dns(0,1:ns-1,0) / (vbs0(1:ns-1) + 1e-8_d)))
    ELSE
        local_dt = 1e10_d
    END IF
    CALL MPI_ALLREDUCE(local_dt, global_dt, 1, mpireal, MPI_MIN, comm, err)
    dt_max_driver = MIN(dt_max_driver, global_dt)

    ! Determine maximum timestep for vbp0:
    IF (surface) THEN
        local_dt = MINVAL(ABS(dnp(0,1:ns,0) / (vbp0(1:ns) + 1e-8_d)))
    ELSE
        local_dt = 1e10_d
    END IF
    CALL MPI_ALLREDUCE(local_dt, global_dt, 1, mpireal, MPI_MIN, comm, err)
    dt_max_driver = MIN(dt_max_driver, global_dt)

    dt_max_driver = cflfact * dt_max_driver
    
    IF (control) PRINT*,'Read surface flow profiles from file '//netcdf_file// &
        ' --- dt_max [s] =', dt_max_driver
    
END SUBROUTINE read_surface_flows

!===============================================================================
SUBROUTINE init_emergence()
! Declare arrays needed for region emergence and set emergence time in secs.

    ALLOCATE(als_em(1:ns,0:np))
    ALLOCATE(alp_em(0:ns,1:np))
    ALLOCATE(als0(1:ns,0:np))
    ALLOCATE(alp0(0:ns,1:np))
    
    t_em = DBLE(cadence_em) * 3600.0_d
    
END SUBROUTINE init_emergence

!===============================================================================
SUBROUTINE begin_emergence()
! Look for emerging region files for this time. If any are present, read them in
! and combine their vector potentials in als_em, alp_em.
! Otherwise, set als_em and alp_em to zero.

    ! Local variables:
    INTEGER:: i
    CHARACTER(4):: nreg
    LOGICAL:: exists

    als_em(1:ns,0:np) = 0.0_d
    alp_em(0:ns,1:np) = 0.0_d
    
    ! Read definitive regions:
    i = 1
    DO
        WRITE(nreg, '(I4.4)') i
        INQUIRE(FILE=TRIM(datadir)//'/rd'//time_to_str(time_now)//'_'//nreg// &
            '.nc', EXIST=exists)
        IF (exists) THEN

            ! Read in data to als0 and alp0:
            CALL ncdf_read_2darray_sp('/rd'//time_to_str(time_now)//'_'//nreg// &
                '.nc', 'as', als0(1:ns,0:np))
            CALL ncdf_read_2darray_sp('/rd'//time_to_str(time_now)//'_'//nreg// &
                '.nc', 'ap', alp0(0:ns,1:np))
            
            ! Add to total:
            als_em(1:ns,0:np) = als_em(1:ns,0:np) + als0(1:ns,0:np)
            alp_em(0:ns,1:np) = alp_em(0:ns,1:np) + alp0(0:ns,1:np)

            IF (control) PRINT*,'Emerging rd'//time_to_str(time_now)//'_'// &
                nreg//'.nc'
        ELSE
            EXIT
        END IF
        i = i + 1
    END DO

    ! Read tentative regions [treat the same]:
    i = 1
    DO
        WRITE(nreg, '(I4.4)') i
        INQUIRE(FILE=TRIM(datadir)//'/rt'//time_to_str(time_now)//'_'//nreg// &
            '.nc', EXIST=exists)
        IF (exists) THEN

            ! Read in data to als0 and alp0:
            CALL ncdf_read_2darray_sp('/rt'//time_to_str(time_now)//'_'//nreg// &
                '.nc', 'as', als0(1:ns,0:np))
            CALL ncdf_read_2darray_sp('/rt'//time_to_str(time_now)//'_'//nreg// &
                '.nc', 'ap', alp0(0:ns,1:np))
            
            ! Add to total:
            als_em(1:ns,0:np) = als_em(1:ns,0:np) + als0(1:ns,0:np)
            alp_em(0:ns,1:np) = alp_em(0:ns,1:np) + alp0(0:ns,1:np)

            IF (control) PRINT*,'Emerging rt'//time_to_str(time_now)//'_'// &
                nreg//'.nc'
        ELSE
            EXIT
        END IF
        i = i + 1
    END DO

END SUBROUTINE begin_emergence

!*******************************************************************************
END MODULE driver
!*******************************************************************************
