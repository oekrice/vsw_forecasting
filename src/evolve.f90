!*******************************************************************************
MODULE evolve
!*******************************************************************************
! This module implements the time evolution.
!*******************************************************************************
    USE mpitools, ONLY: control
    USE params
    USE grid, ONLY: nr, ns, np, Sbr, Sbs, Sbp
    USE time
    USE friction, ONLY: dt_max_fric
    USE outflow, ONLY: dt_max_outflow, vr
    USE ohmic, ONLY: dt_max_ohmic, ohmic_flag
    USE hyper, ONLY: dt_max_hyper, hyper_flag, compute_alpha
    USE vecpot, ONLY: alr, als, alp
    USE efield, ONLY: elr, els, elp, egr, egs, egp, add_outflow, add_friction, &
        add_ohmic, add_hyper, add_surface_flows, add_surface_diffusion, &
        add_emergence, zero_efield, mpi_transfer_e, compute_eg_average
    USE bfield, ONLY: br, bs, bp, bgr, bgs, bgp, bbgs, compute_b_stokes
    USE current, ONLY: jr, js, jp, jgr, jgs, jgp, compute_j_stokes
    USE driver, ONLY: driver_em, driver_v, dt_max_driver, cadence_em, &
        begin_emergence
    USE ncio, ONLY: ncdf_init_snapshot_file, ncdf_write_3darray
    USE diagnostics, ONLY: diagnostic_cadence, write_diagnostics

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: run_simulation

!*******************************************************************************
    REAL(d):: dt    ! timestep [s]
    REAL(d), ALLOCATABLE:: k1alr(:,:,:), k1als(:,:,:), k1alp(:,:,:)
    REAL(d), ALLOCATABLE:: k2alr(:,:,:), k2als(:,:,:), k2alp(:,:,:)
    INTEGER, PARAMETER:: TIMESTEPPING_SCHEME = 2

!*******************************************************************************
CONTAINS
!*******************************************************************************

!===============================================================================
SUBROUTINE run_simulation()
! Run the complete simulation. This is the main loop.

    ! Local variables:
    INTEGER:: hours_since_driver
    INTEGER:: hours_since_snap
    INTEGER:: hours_to_evol
    INTEGER:: steps_since_diagnostics
    REAL(d):: seconds_target
    REAL(d):: dt_lim

    hours_since_snap = 0
    hours_since_driver = cadence_em
    steps_since_diagnostics = diagnostic_cadence

    ! Allocate additional arrays depending on timestepping scheme:
    SELECT CASE (TIMESTEPPING_SCHEME)
    CASE (2) ! Explicit Midpoint
        ALLOCATE(k1alr(1:nr,0:ns,0:np))
        ALLOCATE(k1als(0:nr,1:ns,0:np))
        ALLOCATE(k1alp(0:nr,0:ns,1:np))
    CASE (3) ! Ralston
        ALLOCATE(k1alr(1:nr,0:ns,0:np))
        ALLOCATE(k1als(0:nr,1:ns,0:np))
        ALLOCATE(k1alp(0:nr,0:ns,1:np))
        ALLOCATE(k2alr(1:nr,0:ns,0:np))
        ALLOCATE(k2als(0:nr,1:ns,0:np))
        ALLOCATE(k2alp(0:nr,0:ns,1:np))
    END SELECT

    ! Write out initial snapshot:
    IF (.NOT.restart) THEN
        CALL zero_efield
        CALL compute_efield(alr, als, alp)
        CALL write_3dsnapshot('b'//time_to_str(time_now)//'.nc', &
        snapshot_variables)
    END IF
    
    DO
        ! Apply driver (if due):
        IF (driver_em .AND. (hours_since_driver == cadence_em)) THEN
            CALL begin_emergence()
            hours_since_driver = 0
        END IF

        ! Decide how long to evolve for:
        IF (driver_em) THEN
            hours_to_evol = MIN((snapshot_cadence - hours_since_snap), &
                (cadence_em - hours_since_driver))
        ELSE
            hours_to_evol = snapshot_cadence - hours_since_snap
        END IF

        IF ((hours + hours_to_evol) > hours_total) &
            hours_to_evol = hours_total - hours
        
        ! Evolve by repeated timesteps:
        seconds_target = seconds + REAL(hours_to_evol) * 3600.0_d
        DO
            ! Choose timestep length:
            dt_lim = seconds_target - seconds
            IF (ABS(dt_lim) < 1.0d-8) EXIT
            CALL get_dt(dt_lim)

            ! Carry out single timestep:
            CALL single_timestep(dt, TIMESTEPPING_SCHEME)

            ! Write diagnostics:
            IF (steps_since_diagnostics == diagnostic_cadence) THEN
                CALL write_diagnostics('diagnostics.nc', seconds)
                steps_since_diagnostics = 0
            END IF

            ! Advance timer:
            seconds = seconds + dt
            steps_since_diagnostics = steps_since_diagnostics + 1
        END DO

        ! Update the time/date:
        hours_since_snap = hours_since_snap + hours_to_evol
        hours_since_driver = hours_since_driver + hours_to_evol
        hours = hours + hours_to_evol
        CALL increase_time(time_now, hours_to_evol)

        ! Output:
        IF ((hours_since_snap == snapshot_cadence) &
            .OR. (time_now == time_end)) THEN
            CALL write_3dsnapshot('b'//time_to_str(time_now)//'.nc', &
                snapshot_variables)
            hours_since_snap = 0
        END IF
        
        ! End of simulation:
        IF (time_now == time_end) EXIT
    END DO

END SUBROUTINE run_simulation

!===============================================================================
SUBROUTINE get_dt(dt_lim)
! Determine timestep by CFL condition (cap maximum at dt_lim).
    REAL(d), INTENT(IN):: dt_lim

    ! Maximum timestep allowed by friction and outflow:
    dt = MIN(dt_max_fric, dt_max_outflow)

    ! Maximum timestep allowed by surface flows:
    dt = MIN(dt, dt_max_driver)

    ! Maximum timestep allowed by ohmic diffusion:
    dt = MIN(dt, dt_max_ohmic)
    
    ! Maximum timestep allowed by hyperdiffusion:
    dt = MIN(dt, dt_max_hyper)

    dt = MIN(dt, dt_lim)

END SUBROUTINE get_dt

!===============================================================================
SUBROUTINE single_timestep(dt, scheme)
! Update the vector potential by a single timestep size dt, using an explicit
! timestepping method.
    REAL(d), INTENT(IN):: dt    ! Time stepsize [s]
    INTEGER:: scheme     ! Timestepping scheme [1=FwdEuler, 
                         !    2=ExplMidpoint, 3=Ralston]

    SELECT CASE (scheme)
    CASE (1) ! Forward Euler
        CALL zero_efield
        CALL compute_efield(alr, als, alp)
        alr(1:nr,0:ns,0:np) = alr(1:nr,0:ns,0:np) - dt*elr(1:nr,0:ns,0:np)
        als(0:nr,1:ns,0:np) = als(0:nr,1:ns,0:np) - dt*els(0:nr,1:ns,0:np)
        alp(0:nr,0:ns,1:np) = alp(0:nr,0:ns,1:np) - dt*elp(0:nr,0:ns,1:np)
    CASE (2) ! Explicit midpoint method (second-order RK)
        CALL zero_efield
        CALL compute_efield(alr, als, alp)
        k1alr(1:nr,0:ns,0:np) = alr(1:nr,0:ns,0:np) - 0.5_d*dt*elr(1:nr,0:ns,0:np)
        k1als(0:nr,1:ns,0:np) = als(0:nr,1:ns,0:np) - 0.5_d*dt*els(0:nr,1:ns,0:np)
        k1alp(0:nr,0:ns,1:np) = alp(0:nr,0:ns,1:np) - 0.5_d*dt*elp(0:nr,0:ns,1:np)
        CALL zero_efield
        CALL compute_efield(k1alr, k1als, k1alp)
        alr(1:nr,0:ns,0:np) = alr(1:nr,0:ns,0:np) - dt*elr(1:nr,0:ns,0:np)
        als(0:nr,1:ns,0:np) = als(0:nr,1:ns,0:np) - dt*els(0:nr,1:ns,0:np)
        alp(0:nr,0:ns,1:np) = alp(0:nr,0:ns,1:np) - dt*elp(0:nr,0:ns,1:np)
    CASE (3) ! Ralston's method (second-order RK)
        CALL zero_efield
        CALL compute_efield(alr, als, alp)
        k1alr(1:nr,0:ns,0:np) = elr(1:nr,0:ns,0:np)
        k1als(0:nr,1:ns,0:np) = els(0:nr,1:ns,0:np)
        k1alp(0:nr,0:ns,1:np) = elp(0:nr,0:ns,1:np)
        k2alr(1:nr,0:ns,0:np) = alr(1:nr,0:ns,0:np) - 2.0_d/3.0_d*dt*elr(1:nr,0:ns,0:np)
        k2als(0:nr,1:ns,0:np) = als(0:nr,1:ns,0:np) - 2.0_d/3.0_d*dt*els(0:nr,1:ns,0:np)
        k2alp(0:nr,0:ns,1:np) = alp(0:nr,0:ns,1:np) - 2.0_d/3.0_d*dt*elp(0:nr,0:ns,1:np)
        CALL zero_efield
        CALL compute_efield(k2alr, k2als, k2alp)
        alr(1:nr,0:ns,0:np) = alr(1:nr,0:ns,0:np) - 0.25_d*dt*k1alr(1:nr,0:ns,0:np) - 0.75_d*dt*elr(1:nr,0:ns,0:np)
        als(0:nr,1:ns,0:np) = als(0:nr,1:ns,0:np) - 0.25_d*dt*k1als(0:nr,1:ns,0:np) - 0.75_d*dt*els(0:nr,1:ns,0:np)
        alp(0:nr,0:ns,1:np) = alp(0:nr,0:ns,1:np) - 0.25_d*dt*k1alp(0:nr,0:ns,1:np) - 0.75_d*dt*elp(0:nr,0:ns,1:np)
    CASE DEFAULT
        PRINT*, 'Error, unrecognised value scheme=',scheme, &
            ' in evolve.f90/single_timestep'
    END SELECT

END SUBROUTINE single_timestep

!===============================================================================
SUBROUTINE compute_efield(alr, als, alp)
! Compute electric field from given vector potential.
    REAL(d), INTENT(IN):: alr(1:nr,0:ns,0:np)   ! A_rho*L_rho on (interior) edges
    REAL(d), INTENT(IN):: als(0:nr,1:ns,0:np)   ! A_s*L_s on (interior) edges
    REAL(d), INTENT(IN):: alp(0:nr,0:ns,1:np)   ! A_phi*L_phi on (interior) edges

    CALL compute_b_stokes(alr, als, alp)
    CALL compute_j_stokes(br, bs, bp)
    CALL add_outflow(bs, bp, vr)
    CALL add_friction(bgr, bgs, bgp, bbgs, jgr, jgs, jgp)
    IF (hyper_flag) THEN
        CALL compute_alpha(bgr, bgs, bgp, bbgs, jgr, jgs, jgp)
        CALL add_hyper(bgr, bgs, bgp, bbgs)
    END IF
    IF (ohmic_flag) CALL add_ohmic(jr, js, jp, jgr, jgs, jgp, bbgs)
    IF (driver_v) THEN
        CALL add_surface_flows(br)
        CALL add_surface_diffusion(br)
    END IF
    IF (driver_em) CALL add_emergence()
    CALL mpi_transfer_e
    
END SUBROUTINE compute_efield

!===============================================================================
SUBROUTINE write_3dsnapshot(netcdf_file, variables)
! Write snapshot to netcdf file.
    CHARACTER(*), INTENT(IN):: netcdf_file     ! netcdf file to create/overwrite
    CHARACTER(*), INTENT(IN):: variables    ! string specifying variables

    ! Local variables:
    LOGICAL:: b_flag, j_flag, e_flag
    
    b_flag = (INDEX(variables, 'B') > 0) .OR. (INDEX(variables, 'b') > 0)
    j_flag = (INDEX(variables, 'J') > 0) .OR. (INDEX(variables, 'j') > 0)
    e_flag = (INDEX(variables, 'E') > 0) .OR. (INDEX(variables, 'e') > 0)

    ! Create netcdf file:
    CALL ncdf_init_snapshot_file(netcdf_file, b_flag)
    
    IF (b_flag) THEN
        ! Write magnetic field on cell faces:
        CALL ncdf_write_3darray(netcdf_file, 'br', 'r', 'thc', 'phc', &
            br(0:nr,0:ns+1,0:np+1))
        CALL ncdf_write_3darray(netcdf_file, 'bth', 'rc', 'th', 'phc', &
            -bs(0:nr+1,0:ns,0:np+1))
        CALL ncdf_write_3darray(netcdf_file, 'bph', 'rc', 'thc', 'ph', &
            bp(0:nr+1,0:ns+1,0:np))
    END IF
    
    IF (j_flag) THEN
        ! Write current density at grid points:
        CALL ncdf_write_3darray(netcdf_file, 'jr', 'r', 'th', 'ph', &
            jgr(0:nr,0:ns,0:np))
        CALL ncdf_write_3darray(netcdf_file, 'jth', 'r', 'th', 'ph', &
            -jgs(0:nr,0:ns,0:np))
        CALL ncdf_write_3darray(netcdf_file, 'jph', 'r', 'th', 'ph', &
            jgp(0:nr,0:ns,0:np))
    END IF
    
    IF (e_flag) THEN
        ! Compute total electric field at grid points (unweighted):
        CALL compute_eg_average()
        ! Write electric field at grid points:
        CALL ncdf_write_3darray(netcdf_file, 'er', 'r', 'th', 'ph', &
            egr(0:nr,0:ns,0:np))
        CALL ncdf_write_3darray(netcdf_file, 'eth', 'r', 'th', 'ph', &
            -egs(0:nr,0:ns,0:np))
        CALL ncdf_write_3darray(netcdf_file, 'eph', 'r', 'th', 'ph', &
            egp(0:nr,0:ns,0:np))
    END IF
    
    IF (control) PRINT*,'Wrote snapshot to file '//TRIM(netcdf_file)
    
END SUBROUTINE write_3dsnapshot

!*******************************************************************************
END MODULE evolve
!*******************************************************************************
