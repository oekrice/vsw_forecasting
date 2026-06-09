!*******************************************************************************
MODULE init
!*******************************************************************************
! Routines for reading input parameters and initializing the simulation.
!*******************************************************************************
    USE params
    USE mpitools, ONLY: control, comm, err, n_procs, mpireal, MPI_MIN, &
        start_mpi, finish, subdivide_mpi
    USE time, ONLY: start_time, end_time, get_restart_time, initialize_time, &
        restart_time
    USE grid, ONLY: nr, ns, np, nrg, nsg, npg, rss, dlp, &
        create_grid_arrays
    USE bfield, ONLY: read_magnetic_field, br, bs, bp
    USE vecpot, ONLY: read_vector_potential, compute_a_from_b
    USE friction, ONLY: nu0, read_friction
    USE outflow, ONLY: read_outflow
    USE ohmic, ONLY: ohmic_flag, ohmic_ratio, ohmic_enhanced, initialize_ohmic
    USE hyper, ONLY: hyper_flag, hyper_ratio, initialize_hyper
    USE driver, ONLY: cadence_em, driver_em, driver_v, read_surface_flows, &
        init_emergence
    USE diagnostics, ONLY: diagnostic_cadence
    
    IMPLICIT NONE
    PRIVATE

    PUBLIC:: initialize

!*******************************************************************************

!*******************************************************************************
CONTAINS
!*******************************************************************************
    
!===============================================================================
SUBROUTINE initialize()
! Read in the simulation parameters to initialize a simulation.

    CALL start_mpi
    
    CALL read_user_namelist
    CALL read_setup_namelist
    CALL determine_restart
    CALL save_code

    CALL create_grid_arrays

    IF (restart) THEN
        CALL get_restart_time
        CALL read_magnetic_field('b'//restart_time//'.nc')
        CALL compute_a_from_b(br, bs, bp)
    ELSE
        CALL read_vector_potential('a_pfss'//start_time//'.nc')
    END IF

    CALL read_friction('friction.nc')
    CALL read_outflow('outflow.nc')
    CALL initialize_ohmic
    CALL initialize_hyper

    IF (driver_v) CALL read_surface_flows('surface_flows.nc')
    IF (driver_em) CALL init_emergence()

    CALL initialize_time

END SUBROUTINE initialize

!===============================================================================
SUBROUTINE read_user_namelist()
! Read code parameters from user namelist supplied at command line.

    ! Local variables:
    CHARACTER(str_max_len):: user_nml_file
    INTEGER:: tmp_rank, i
    NAMELIST /user/ datadir, cflfact, snapshot_cadence, snapshot_variables, &
        diagnostic_cadence, nu0, ohmic_ratio, ohmic_enhanced, hyper_ratio

    ! Read in name of namelist file from command line:
    CALL GET_COMMAND_ARGUMENT(1, user_nml_file)

    ! Each process reads namelist file in turn, to avoid conflicts:
    CALL MPI_COMM_RANK(comm, tmp_rank, err)
    DO i = 0, n_procs - 1
        IF (tmp_rank == i) THEN

            ! Read in user-defined namelist:
            OPEN(1, FILE=TRIM(user_nml_file), STATUS='old', IOSTAT=err)
            IF (err /= 0) THEN
                PRINT*, 'ERROR: could not open user namelist file '// &
                    TRIM(user_nml_file)
                CALL finish
            ELSE
                READ(1, NML=user)
                CLOSE(1)
            END IF
        END IF
    END DO

    IF (control) PRINT*,'Read user namelist: ', TRIM(user_nml_file)
    IF (control) PRINT*,'Data directory: ', TRIM(datadir)

END SUBROUTINE read_user_namelist

!===============================================================================
SUBROUTINE read_setup_namelist()
! Read code parameters from setup namelist 'prepared.nml' in data directory.
! This must be pre-generated, usually from a python prepare script.

    ! Local variables:
    INTEGER:: tmp_rank, i
    NAMELIST /prepared/ start_time, end_time, cadence_em, nr, ns, np, rss, &
        driver_v, driver_em

    ! Each process reads namelist file in turn, to avoid conflicts:
    CALL MPI_COMM_RANK(comm, tmp_rank, err)
    DO i = 0, n_procs - 1
        IF (tmp_rank == i) THEN

            ! Read in setup namelist:
            OPEN(1, FILE=TRIM(datadir)//'/prepared.nml', STATUS='old', &
                IOSTAT=err)
            IF (err /= 0) THEN
                PRINT*, 'ERROR: problem with data directory.'
                CALL finish
            ELSE
                READ(1, NML=prepared)
                CLOSE(1)
            END IF
        END IF
    END DO
    
    ! Global number of grid cells:
    nrg = nr
    nsg = ns
    npg = np
    
    ! Subdivide grid between MPI processes and set nr, ns, np to local
    ! number of grid cells:
    CALL subdivide_mpi(nrg, nsg, npg, nr, ns, np)

    IF (control) PRINT*,'Global grid resolution [rho, s, phi]:', nrg, nsg, npg
    IF (control) PRINT*,'Outer boundary r_ss [R_sun] =',rss
    IF (control) PRINT*,'Start time: ',start_time,' End time: ',end_time
    IF (control) &
        PRINT*,'==============================================================='

END SUBROUTINE read_setup_namelist

!===============================================================================
SUBROUTINE determine_restart()
! Check for command line argument "-restart".

    ! Local variables:
    CHARACTER(str_max_len):: rstr
    
    CALL GET_COMMAND_ARGUMENT(2, rstr)
    IF (rstr(1:8) == '-restart') restart = .TRUE.

END SUBROUTINE determine_restart

!===============================================================================
SUBROUTINE save_code()
    ! Save copy of source code and user namelist file.
    
    ! Local variables:
    CHARACTER(str_max_len):: user_nml_file

    ! Read in name of namelist file from command line:
    CALL GET_COMMAND_ARGUMENT(1, user_nml_file)

    IF (control) THEN
        CALL SYSTEM('rm -rf '//TRIM(datadir)//'/src')
        CALL SYSTEM('cp -r src '//TRIM(datadir)//'/')
        CALL SYSTEM('cp Makefile '//TRIM(datadir)//'/')
        CALL SYSTEM('cp '//TRIM(user_nml_file)//' '//TRIM(datadir)//'/user.nml')
    END IF
        
END SUBROUTINE save_code

!*******************************************************************************
END MODULE init
!*******************************************************************************
