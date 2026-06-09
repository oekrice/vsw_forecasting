!*******************************************************************************
MODULE diagnostics
!*******************************************************************************
! This module implements the diagnostics file.
!*******************************************************************************
    USE netcdf
    USE params
    USE mpitools, ONLY: control, surface, outer, mpireal, comm, err, MPI_SUM, &
        MPI_INTEGER, MPI_MAX, north, south, rank
    USE grid, ONLY: nr, ns, np, dr, ds, dp, rg, sg, sc, pc, Sbr, Vg, V_total, dls, &
        dlp, rss, Sjr, Sjs, Sjp
    USE ncio, ONLY: ncdf_open_file, ncdf_close_file, &
        ncdf_create_diagnostics_file, ncdf_restart_diagnostics_file, &
        ncdf_write_diagnostic_entry
    USE time, ONLY: start_time, end_time
    USE driver, ONLY: cadence_em
    USE bfield, ONLY: br, bs, bp, bbg, bgr, bgs, bgp, bbgs
    USE current, ONLY: jgr, jgs, jgp
    USE efield, ONLY: els, elp, fgr, fgs, fgp, compute_eg_average, egs, egp
    USE vecpot, ONLY: zero_aref_outer, compute_aref_outer, arefs, arefp, alr, als, alp
    USE hyper, ONLY: hyper_eta4, hr, hs, hp, alg, compute_alpha, compute_hyper_div
    USE friction, ONLY: nu
    USE outflow, ONLY: vr

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: diagnostic_cadence, write_diagnostics, magnetic_energy, &
        unsigned_flux_outer, test_diagnostics

!*******************************************************************************
    INTEGER:: diagnostic_cadence    ! cadence to output diagnostics [time steps]
    REAL(d), PARAMETER:: RSUN = 6.96e10_d     ! radius of Sun [cm]
    REAL(d), PARAMETER:: MU0 = 1.2566370614e-8_d   ! vacuum permeability [H/cm]
    
    LOGICAL, PARAMETER:: unsigned_flux_inner_flag = .TRUE.
    LOGICAL, PARAMETER:: unsigned_flux_outer_flag = .TRUE.
    LOGICAL, PARAMETER:: magnetic_energy_flag = .TRUE.
    LOGICAL, PARAMETER:: mean_current_flag = .TRUE.
    LOGICAL, PARAMETER:: outer_bh_flag = .FALSE.
    LOGICAL, PARAMETER:: helicity_flux_flag = .FALSE.
    LOGICAL, PARAMETER:: energy_terms_flag = .TRUE.
    LOGICAL, PARAMETER:: eruption_diagnostics_flag = .TRUE.
    
    REAL(d), PARAMETER:: MAX_TENSION_RADIUS = 2.0_d
    REAL(d), PARAMETER:: BH_MIN = 0.02_d

    REAL(d), ALLOCATABLE:: darefs(:,:)
    REAL(d), ALLOCATABLE:: darefp(:,:)

    INTEGER, PARAMETER:: n_caps = 18
    REAL(d), PARAMETER:: cap_lats(n_caps)=(/90.0_d, 35.0_d, 35.0_d, 35.0_d, 35.0_d, 0.0_d, 0.0_d,&
                0.0_d, 0.0_d, 0.0_d, 0.0_d, 0.0_d, 0.0_d, -35.0_d,&
                -35.0_d, -35.0_d, -35.0_d, -90.0_d/)
    REAL(d), PARAMETER:: cap_lons(n_caps)=(/0.0_d, 0.0_d, 90.0_d, 180.0_d, 270.0_d, 0.0_d, 45.0_d,&
                90.0_d, 135.0_d, 180.0_d, 225.0_d, 270.0_d, 315.0_d,&
                0.0_d, 90.0_d, 180.0_d, 270.0_d, 0.0_d/)
    REAL(d), PARAMETER:: cap_half_angle = 0.65_d
    REAL(d), ALLOCATABLE:: cap_partial_sum(:), local_helicity(:,:)
    REAL(d), ALLOCATABLE:: cap_mask(:,:,:), cap_helicity_flux(:)

!*******************************************************************************
CONTAINS
!*******************************************************************************

!===============================================================================
SUBROUTINE write_diagnostics(diag_file, t)
! Write diagnostics for current time to file.
    CHARACTER(*), INTENT(IN):: diag_file    ! diagnostics file
    REAL(d), INTENT(IN):: t     ! current time [s]
    
    ! Local variables:
    INTEGER, SAVE:: position = 1   ! index to entry in the file
    REAL(d):: flux_inner, flux_outer, energy, current, bsn, bss, bpn, bps, hflux
    REAL(d):: poynting_inner, poynting_outer, maxten, nptsh
    REAL(d):: energy_flux_hyper_inner, energy_flux_hyper_outer
    REAL(d):: diss_friction, work_outflow, diss_hyper
    INTEGER:: nc_id, i
    LOGICAL:: file_exists
    CHARACTER(40):: cap_varname

    ! Initialize the diagnostics file if this is the first call:
    IF (position == 1) THEN
        IF (control) THEN
            INQUIRE(FILE=TRIM(datadir)//'/'//diag_file, EXIST=file_exists)
            ! If file already exists, then update it, otherwise create a new file:
            IF (restart .AND. file_exists) THEN
                CALL ncdf_restart_diagnostics_file(diag_file, t, position)
            ELSE
                CALL ncdf_create_diagnostics_file(diag_file, start_time, end_time, &
                    cadence_em, unsigned_flux_inner_flag, unsigned_flux_outer_flag, &
                    magnetic_energy_flag, mean_current_flag, outer_bh_flag, &
                    helicity_flux_flag, energy_terms_flag, eruption_diagnostics_flag, &
                    n_caps)
            END IF
        END IF
        ! Initialize reference A (for helicity flux):
        IF (helicity_flux_flag) CALL zero_aref_outer()
   END IF

    ! Compute the diagnostics:
    IF (unsigned_flux_inner_flag) flux_inner = unsigned_flux_inner(br)
    IF (unsigned_flux_outer_flag) flux_outer = unsigned_flux_outer(br)
    IF (magnetic_energy_flag) energy = magnetic_energy(bbg)
    IF (mean_current_flag) current = mean_current(jgr, jgs, jgp)
    IF (outer_bh_flag) THEN
        bsn = outer_bs_north(bs)
        bss = outer_bs_south(bs)
        bpn = outer_bp_north(bp)
        bps = outer_bp_south(bp)
    END IF
    IF (helicity_flux_flag) hflux = helicity_flux(els, elp, arefs, arefp, t)
    IF (energy_terms_flag) THEN
        CALL compute_eg_average()     
        ! Use B from inside domain for averages on outer bds:
        bgp(nr,0:ns,0:np) = 0.5_d * (bp(nr,0:ns,0:np) + bp(nr,1:ns+1,0:np))
        bgs(nr,0:ns,0:np) = 0.5_d * (bs(nr,0:ns,0:np) + bs(nr,0:ns,1:np+1))
        egs(0:nr,0:ns,0:np) = egs(0:nr,0:ns,0:np) * bgp(0:nr,0:ns,0:np)
        egp(0:nr,0:ns,0:np) = egp(0:nr,0:ns,0:np) * bgs(0:nr,0:ns,0:np)
        poynting_inner = poynting_flux(egs, egp, .TRUE.)
        poynting_outer = poynting_flux(egs, egp, .FALSE.)
        CALL compute_alpha(bgr, bgs, bgp, bbgs, jgr, jgs, jgp)
        CALL compute_hyper_div(bbgs)
        hr(0:nr+1,0:ns,0:np) = hr(0:nr+1,0:ns,0:np) / Sjr(0:nr+1,0:ns,0:np)
        hs(0:nr,0:ns+1,0:np) = hs(0:nr,0:ns+1,0:np) / Sjs(0:nr,0:ns+1,0:np)
        hp(0:nr,0:ns,0:np+1) = hp(0:nr,0:ns,0:np+1) / Sjp(0:nr,0:ns,0:np+1)
        energy_flux_hyper_inner = energy_flux_hyper(.TRUE.)
        energy_flux_hyper_outer = energy_flux_hyper(.FALSE.)
        diss_friction = dissipation_friction()
        work_outflow = workdone_outflow()
        diss_hyper = dissipation_hyper()
    END IF
    IF (eruption_diagnostics_flag) THEN
        maxten = max_tension(MAX_TENSION_RADIUS)
        nptsh = npoints_bh(BH_MIN)
    END IF

    ! Write to file:
    IF (control) THEN
        nc_id = ncdf_open_file(diag_file)
        CALL ncdf_write_diagnostic_entry(nc_id, 't', t, position)
        IF (unsigned_flux_inner_flag) &
            CALL ncdf_write_diagnostic_entry(nc_id, 'phot_flux', &
                flux_inner, position)
        IF (unsigned_flux_outer_flag) &
                    CALL ncdf_write_diagnostic_entry(nc_id, 'open_flux', &
                        flux_outer, position)
        IF (magnetic_energy_flag) CALL ncdf_write_diagnostic_entry(nc_id, &
            'energy', energy, position)
        IF (mean_current_flag) CALL ncdf_write_diagnostic_entry(nc_id, &
                'mean_current', current, position)
        IF (outer_bh_flag) THEN
            CALL ncdf_write_diagnostic_entry(nc_id, 'bs_nor', bsn, position)
            CALL ncdf_write_diagnostic_entry(nc_id, 'bs_sou', bss, position)
            CALL ncdf_write_diagnostic_entry(nc_id, 'bp_nor', bpn, position)
            CALL ncdf_write_diagnostic_entry(nc_id, 'bp_sou', bps, position)
        END IF
        IF (helicity_flux_flag) THEN
            CALL ncdf_write_diagnostic_entry(nc_id, 'hr_flux', hflux, position)
            DO i=1,n_caps
                WRITE(cap_varname,"('hr_flux_',I0.3)") i
                CALL ncdf_write_diagnostic_entry(nc_id, TRIM(cap_varname), &
                    cap_helicity_flux(i), position)
            END DO
        END IF
        IF (energy_terms_flag) THEN
            CALL ncdf_write_diagnostic_entry(nc_id, 'poynting_inner', &
                poynting_inner, position)
            CALL ncdf_write_diagnostic_entry(nc_id, 'poynting_outer', &
                poynting_outer, position)
            CALL ncdf_write_diagnostic_entry(nc_id, 'energy_flux_hyper_inner', &
                energy_flux_hyper_inner, position)
            CALL ncdf_write_diagnostic_entry(nc_id, 'energy_flux_hyper_outer', &
                energy_flux_hyper_outer, position)
            CALL ncdf_write_diagnostic_entry(nc_id, 'dissipation_friction', &
                diss_friction, position)      
            CALL ncdf_write_diagnostic_entry(nc_id, 'work_outflow', &
                work_outflow, position)
            CALL ncdf_write_diagnostic_entry(nc_id, 'dissipation_hyper', &
                diss_hyper, position)                         
        END IF
        IF (eruption_diagnostics_flag) THEN
            CALL ncdf_write_diagnostic_entry(nc_id, 'max_tension', maxten, position)
            CALL ncdf_write_diagnostic_entry(nc_id, 'npoints_bh', nptsh, position)
        END IF
        CALL ncdf_close_file(nc_id)
        position = position + 1
    END IF
    
END SUBROUTINE write_diagnostics

!===============================================================================
FUNCTION unsigned_flux_inner(br)
! Compute unsigned flux at inner boundary.
    REAL(d):: unsigned_flux_inner   ! Total unsigned flux [Mx]
    REAL(d), INTENT(IN):: br(0:nr,0:ns+1,0:np+1)    ! B_rho on faces [G]

    ! Local variables:
    REAL(d):: local_tot

    local_tot = 0.0_d
    IF (surface) local_tot = SUM(ABS(br(0,1:ns,1:np))*Sbr(0,1:ns,1:np))
    CALL MPI_ALLREDUCE(local_tot, unsigned_flux_inner, 1, mpireal, MPI_SUM, &
        comm, err)
    
    unsigned_flux_inner = unsigned_flux_inner * RSUN * RSUN

END FUNCTION unsigned_flux_inner

!===============================================================================
FUNCTION unsigned_flux_outer(br)
! Compute unsigned flux at outer boundary.
    REAL(d):: unsigned_flux_outer    ! Total unsigned flux [Mx]
    REAL(d), INTENT(IN):: br(0:nr,0:ns+1,0:np+1)    ! B_rho on faces [G]

    ! Local variables:
    REAL(d):: local_tot

    local_tot = 0.0_d
    IF (outer) local_tot = SUM(ABS(br(nr,1:ns,1:np)*Sbr(nr,1:ns,1:np)))
    CALL MPI_ALLREDUCE(local_tot, unsigned_flux_outer, 1, mpireal, MPI_SUM, &
        comm, err)

    unsigned_flux_outer = unsigned_flux_outer * RSUN * RSUN

END FUNCTION unsigned_flux_outer

!===============================================================================
FUNCTION magnetic_energy(bbg)
! Compute total magnetic energy in simulation volume.
    REAL(d):: magnetic_energy    ! Magnetic energy [ergs]
    REAL(d), INTENT(IN):: bbg(0:nr,0:ns,0:np)    ! B**2 at grid pts [G]

    ! Local variables:
    REAL(d):: local_tot

    local_tot = 0.125_d * SUM( &
        bbg(0:nr-1,0:ns-1,0:np-1) * Vg(0:nr-1,0:ns-1,0:np-1) &
        + bbg(0:nr-1,0:ns-1,1:np) * Vg(0:nr-1,0:ns-1,1:np) &
        + bbg(0:nr-1,1:ns,0:np-1) * Vg(0:nr-1,1:ns,0:np-1)  &
        + bbg(0:nr-1,1:ns,1:np) * Vg(0:nr-1,1:ns,1:np) &
        + bbg(1:nr,0:ns-1,0:np-1) * Vg(1:nr,0:ns-1,0:np-1) &
        + bbg(1:nr,0:ns-1,1:np) * Vg(1:nr,0:ns-1,1:np) &
        + bbg(1:nr,1:ns,0:np-1) * Vg(1:nr,1:ns,0:np-1) &
        + bbg(1:nr,1:ns,1:np) * Vg(1:nr,1:ns,1:np))
    CALL MPI_ALLREDUCE(local_tot, magnetic_energy, 1, mpireal, MPI_SUM, &
        comm, err)

    magnetic_energy = magnetic_energy * RSUN**3 / (8.0_d * PI)
    
END FUNCTION magnetic_energy

!===============================================================================
FUNCTION mean_current(jgr, jgs, jgp)
! Compute total magnetic energy in simulation volume.
    REAL(d):: mean_current    ! mean of |J| [A/m**2]
    REAL(d), INTENT(IN):: jgr(0:nr,0:ns,0:np)    ! J_rho at grid pts [G/R_sun]
    REAL(d), INTENT(IN):: jgs(0:nr,0:ns,0:np)    ! J_s at grid pts [G/R_sun]
    REAL(d), INTENT(IN):: jgp(0:nr,0:ns,0:np)    ! J_phi at grid pts [G/R_sun]

    ! Local variables:
    REAL(d):: local_tot
    REAL(d) :: jg(0:nr,0:ns,0:np)    ! |J| at grid pts [G/R_sun]

    ! Compute |J| at grid points weighted by cell volumes:
    jg(0:nr,0:ns,0:np) = SQRT(jgr(0:nr,0:ns,0:np) * jgr(0:nr,0:ns,0:np) &
        + jgs(0:nr,0:ns,0:np) * jgs(0:nr,0:ns,0:np) &
        + jgp(0:nr,0:ns,0:np) * jgp(0:nr,0:ns,0:np)) * Vg

    ! Volume mean:
    local_tot = 0.125_d * SUM(jg(0:nr-1,0:ns-1,0:np-1) &
        + jg(0:nr-1,0:ns-1,1:np) + jg(0:nr-1,1:ns,0:np-1) &
        + jg(0:nr-1,1:ns,1:np) + jg(1:nr,0:ns-1,0:np-1) &
        + jg(1:nr,0:ns-1,1:np) + jg(1:nr,1:ns,0:np-1) + jg(1:nr,1:ns,1:np))
    CALL MPI_ALLREDUCE(local_tot, mean_current, 1, mpireal, MPI_SUM, &
        comm, err)
    mean_current = mean_current / V_total
    
    ! Convert to A/m**2 [note that G/H = 1e-4 A/m**2]:
!    mean_current = mean_current / MU0/ RSUN * 1e-4_d
    mean_current = mean_current /RSUN ! for comparison with old code

END FUNCTION mean_current

!===============================================================================
FUNCTION outer_bs_north(bs)
! Compute average B_s at half a grid cell below outer boundary, in Northern
! hemisphere.
    REAL(d):: outer_bs_north  ! average B_s [G]
    REAL(d), INTENT(IN):: bs(0:nr+1,0:ns,0:np+1)    ! B_s on faces [G]

    ! Local variables:
    INTEGER:: j
    REAL(d):: local_tot, j0, glob_n

    j0 = 0
    IF (outer) j0 = DBLE(ns*np)
    CALL MPI_ALLREDUCE(j0, glob_n, 1, mpireal, MPI_SUM, comm, err)
    local_tot = 0.0_d
    IF (outer) THEN
        DO j = 0, ns-1
            IF (sg(j) > 0) local_tot = local_tot + SUM(bs(nr,j,1:np))
        END DO
    END IF
    CALL MPI_ALLREDUCE(local_tot, outer_bs_north, 1, mpireal, MPI_SUM, comm, &
        err)
    outer_bs_north = outer_bs_north / glob_n
    
END FUNCTION outer_bs_north

!===============================================================================
FUNCTION outer_bs_south(bs)
! Compute average B_s at half a grid cell below outer boundary, in Southern
! hemisphere.
    REAL(d):: outer_bs_south  ! average B_s [G]
    REAL(d), INTENT(IN):: bs(0:nr+1,0:ns,0:np+1)    ! B_s on faces [G]

    ! Local variables:
    INTEGER:: j
    REAL(d):: local_tot, j0, glob_n

    j0 = 0
    IF (outer) j0 = DBLE(ns*np)
    CALL MPI_ALLREDUCE(j0, glob_n, 1, mpireal, MPI_SUM, comm, err)
    local_tot = 0.0_d
    IF (outer) THEN
        DO j = 0, ns-1
            IF (sg(j) < 0) local_tot = local_tot + SUM(bs(nr,j,1:np))
        END DO
    END IF
    CALL MPI_ALLREDUCE(local_tot, outer_bs_south, 1, mpireal, MPI_SUM, comm, &
        err)
    outer_bs_south = outer_bs_south / glob_n
    
END FUNCTION outer_bs_south

!===============================================================================
FUNCTION outer_bp_north(bp)
! Compute average B_phi at half a grid cell below outer boundary, in Northern
! hemisphere.
    REAL(d):: outer_bp_north  ! average B_phi [G]
    REAL(d), INTENT(IN):: bp(0:nr+1,0:ns+1,0:np)    ! B_phi on faces [G]

    ! Local variables:
    INTEGER:: j
    REAL(d):: local_tot, j0, glob_n

    j0 = 0
    IF (outer) j0 = DBLE(ns*np)
    CALL MPI_ALLREDUCE(j0, glob_n, 1, mpireal, MPI_SUM, comm, err)
    local_tot = 0.0_d
    IF (outer) THEN
        DO j = 0, ns
            IF (sc(j) > 0) local_tot = local_tot + SUM(bp(nr,j,0:np-1))
        END DO
    END IF
    CALL MPI_ALLREDUCE(local_tot, outer_bp_north, 1, mpireal, MPI_SUM, comm, &
        err)
    outer_bp_north = outer_bp_north / glob_n
    
END FUNCTION outer_bp_north

!===============================================================================
FUNCTION outer_bp_south(bp)
! Compute average B_phi at half a grid cell below outer boundary, in Southern
! hemisphere.
    REAL(d):: outer_bp_south ! average B_phi [G]
    REAL(d), INTENT(IN):: bp(0:nr+1,0:ns+1,0:np)    ! B_phi on faces [G]

    ! Local variables:
    INTEGER:: j
    REAL(d):: local_tot, j0, glob_n

    j0 = 0
    IF (outer) j0 = DBLE(ns*np)
    CALL MPI_ALLREDUCE(j0, glob_n, 1, mpireal, MPI_SUM, comm, err)
    local_tot = 0.0_d
    IF (outer) THEN
        DO j = 0, ns
            IF (sc(j) < 0) local_tot = local_tot + SUM(bp(nr,j,0:np-1))
        END DO
    END IF
    CALL MPI_ALLREDUCE(local_tot, outer_bp_south, 1, mpireal, MPI_SUM, comm, &
        err)
    outer_bp_south = outer_bp_south / glob_n
    
END FUNCTION outer_bp_south

!===============================================================================
FUNCTION helicity_flux(els, elp, arefs, arefp, t)
! Compute flux of relative helicity through outer boundary, in P-T gauge.
! -- computes both the total, and separate values in each "cap".
    REAL(d):: helicity_flux     ! units [Mx**2/s]
    REAL(d), INTENT(INOUT):: els(0:nr,1:ns,0:np)   ! E_s*L_s on edges
    REAL(d), INTENT(INOUT):: elp(0:nr,0:ns,1:np)   ! E_phi*L_phi on edges
    REAL(d), INTENT(IN):: arefs(1:ns,0:np)   ! A^ref_s*L_s on outer bndry edges
    REAL(d), INTENT(IN):: arefp(0:ns,1:np)   ! A^ref_phi*L_phi on outer bd edges
    REAL(d), INTENT(IN):: t     ! current time [s]

    ! Local variables:
    REAL(d), SAVE:: tprev = -1.0_d
    REAL(d):: local_tot, angle
    INTEGER:: j, k, i, j0, j1
    CHARACTER(40):: test_filename

    IF (.NOT. ALLOCATED(darefs)) ALLOCATE(darefs(1:ns,1:np))
    IF (.NOT. ALLOCATED(darefp)) ALLOCATE(darefp(1:ns,1:np))
    IF (.NOT. ALLOCATED(cap_partial_sum)) ALLOCATE(cap_partial_sum(1:n_caps))
    IF (.NOT. ALLOCATED(local_helicity)) ALLOCATE(local_helicity(1:ns,1:np))
    IF (.NOT. ALLOCATED(cap_mask)) ALLOCATE(cap_mask(1:n_caps,1:ns,1:np))
    IF (.NOT. ALLOCATED(cap_helicity_flux)) ALLOCATE(cap_helicity_flux(1:n_caps))

    IF (tprev >= 0.0_d) THEN
        local_tot = 0.0_d
        IF (outer) THEN
            ! Previous Ap, averaged to r face centres:
            darefs(1:ns,1:np) = 0.5_d * (arefs(1:ns,0:np-1) + arefs(1:ns,1:np))
            darefp(1:ns,1:np) = 0.5_d * (arefp(0:ns-1,1:np) + arefp(1:ns,1:np))
        END IF
        CALL compute_aref_outer(br)
        IF (outer) THEN
            els(0:nr,1:ns,0:np) = els(0:nr,1:ns,0:np) / dls(0:nr,1:ns,0:np)
            j0 = 0
            j1 = ns
            IF (south) j0 = 1
            IF (north) j1 = ns-1
            elp(0:nr,j0:j1,1:np) = elp(0:nr,j0:j1,1:np)/dlp(0:nr,j0:j1,1:np)
            ! (a) Add (2E + dAp/dt)
            darefs(1:ns,1:np) = (0.5_d * (arefs(1:ns,0:np-1) &
                + arefs(1:ns,1:np)) - darefs(1:ns,1:np)) / (t - tprev)
            darefs(1:ns,1:np) = darefs(1:ns,1:np) + els(nr,1:ns,0:np-1) &
                + els(nr,1:ns,1:np)
            darefp(1:ns,1:np) = (0.5_d * (arefp(0:ns-1,1:np) &
                + arefp(1:ns,1:np)) - darefp(1:ns,1:np)) / (t - tprev)
            darefp(1:ns,1:np) = darefp(1:ns,1:np) + elp(nr,0:ns-1,1:np) &
                + elp(nr,1:ns,1:np)
            ! (b) Compute integral of [Ap x (2E + dAp/dt)].e_r
!            do j=1,ns
!                if (sc(j) > 0) local_tot = local_tot + sum((0.5_d*(arefp(j-1,:) + arefp(j,:))*darefs(j,:) &
!                - 0.5_d*(arefs(j,0:np-1) + arefs(j,1:np))*darefp(j,:)) * Sbr(nr,j,1:np))
!             end do
             
            local_helicity(1:ns,1:np) = (0.5_d * (arefp(0:ns-1,1:np) &
                + arefp(1:ns,1:np)) * darefs(1:ns,1:np) - 0.5_d*(arefs(1:ns,0:np-1) &
                + arefs(1:ns,1:np)) * darefp(1:ns,1:np)) * Sbr(nr,1:ns,1:np)
            local_tot = SUM(local_helicity)
            DO i=1, n_caps
                cap_partial_sum(i) = SUM(cap_mask(i,:,:)*local_helicity)
            END DO
            els(0:nr,1:ns,0:np) = els(0:nr,1:ns,0:np) * dls(0:nr,1:ns,0:np)
            elp(0:nr,j0:j1,1:np) = elp(0:nr,j0:j1,1:np) * dlp(0:nr,j0:j1,1:np)
        END IF
        CALL MPI_ALLREDUCE(local_tot, helicity_flux, 1, mpireal, MPI_SUM, &
            comm, err)
        helicity_flux = helicity_flux * RSUN**4
        DO i = 1,n_caps
            CALL MPI_ALLREDUCE(cap_partial_sum(i), cap_helicity_flux(i), 1, &
                mpireal, MPI_SUM, comm, err)
            cap_helicity_flux(i) = cap_helicity_flux(i) * RSUN**4
        END DO
    ELSE
        CALL compute_aref_outer(br)
        helicity_flux = 0.0_d
        cap_helicity_flux(1:n_caps) = 0.0_d
        ! Precompute masks for each of the caps:
        DO i=1,n_caps
            DO k=1,np
                DO j=1,ns
                    ! Angle on unit sphere between this point and cap centre:
                    angle = ACOS( sc(j)*SIN(cap_lats(i)*PI/180.0_d) &
                         + SQRT(1-sc(j)**2)*COS(cap_lats(i)*PI/180.0_d) &
                         * COS(pc(k) - cap_lons(i)*PI/180.0_d) )
                    IF (angle <= cap_half_angle) THEN
                        cap_mask(i,j,k) = 1.0_d
                    ELSE
                        cap_mask(i,j,k) = 0.0_d
                    END IF
                END DO
            END DO
            WRITE(test_filename,"('test',I0.6,'.',I0.6,'.dat')") rank, i
            OPEN (unit=915+rank,file=test_filename,form='unformatted')
            WRITE(915+rank) sc(:)
            WRITE(915+rank) pc(:)
            WRITE(915+rank) cap_mask(i,:,:)
            CLOSE(915+rank)
        END DO
    END IF
    
    tprev = t

END FUNCTION helicity_flux

!===============================================================================
FUNCTION npoints_bh(threshold)
! On spherical surface half a grid cell below outer boundary, compute number
! of grid points where |B_h| = sqrt(Bs**2 + Bph**2) > threshold.
    REAL(d):: npoints_bh  ! number of grid points (as float for convenience)
    REAL(d), INTENT(IN):: threshold    ! minimum value of |B_h| to count [G]

    ! Local variables:
    REAL(d):: local_tot

    local_tot = 0.0_d
    IF (outer) &
        local_tot = COUNT( bgs(nr,0:ns,0:np-1)**2 + bgp(nr,0:ns,0:np-1)**2 &
            > (threshold**2))

    CALL MPI_ALLREDUCE(local_tot, npoints_bh, 1, mpireal, MPI_SUM, comm, err)      
        
END FUNCTION npoints_bh

!===============================================================================
FUNCTION max_tension(radius)
    ! Compute maximum radial tension force at grid point nearest to
    ! given radius. (Returns zero if the maximum is negative.)
    REAL(d):: max_tension  ! maximum vale of T_r at given radius
    REAL(d), INTENT(IN):: radius    ! radius to compute T_r [Rsun]

    ! Local variables:
    REAL(d):: local_max
    INTEGER:: irmin(1), i

    local_max = 0.0_d
    ! Check whether radius is in this process:
    IF (MINVAL(ABS(rg-LOG(radius))) < dr) THEN
        irmin = MINLOC(ABS(rg-LOG(radius)))
        ! Avoid double counting when nearest point is exactly at a process boundary:
        IF (irmin(1) < nr) THEN
            i = irmin(1)
            ! Compute maximum tension force in this process:
            ! T_r = RSUN * e_r . (B.grad(B)) / B**2
            local_max = MAXVAL( ( fgr(i,:,:) + 0.5_d*EXP(-rg(i))*(bbg(i+1,:,:) &
                - bbg(i,:,:))/dr )/bbg(i,:,:) )
        END IF
    END IF
    
    CALL MPI_ALLREDUCE(local_max, max_tension, 1, mpireal, MPI_MAX, comm, err) 

END FUNCTION max_tension

!===============================================================================
FUNCTION poynting_flux(esbpg, epbsg, inner_flag)
    ! Compute  -1/(4*pi) * int ( E x B ) . er dS on inner boundary
    ! (if inner_flag=.TRUE.) or outer boundary (if inner_flag=.FALSE.).
        REAL(d):: poynting_flux     ! units [erg/s]
        LOGICAL, INTENT(IN):: inner_flag ! .TRUE. for inner ,.FALSE. for outer bdy
        REAL(d), INTENT(INOUT):: esbpg(0:nr,0:ns,0:np)   ! E_s*B_phi at grid pts
        REAL(d), INTENT(INOUT):: epbsg(0:nr,0:ns,0:np)   ! E_phi*B_s at grid pts

        ! Local variables:
        REAL(d):: local_tot
        ! REAL(d):: vr_flux
        REAL(d):: poy(1:ns,1:np)    ! Radial poynting vector on faces
        INTEGER:: i

        poy = 0.0_d
        local_tot = 0.0_d
        IF ((inner_flag).AND.(surface)) THEN
            poy(1:ns,1:np) = esbpg(0,0:ns-1,0:np-1) + esbpg(0,0:ns-1,1:np) &
                + esbpg(0,1:ns,0:np-1) + esbpg(0,1:ns,1:np) &
                - epbsg(0,0:ns-1,0:np-1) - epbsg(0,0:ns-1,1:np) &
                - epbsg(0,1:ns,0:np-1) - epbsg(0,1:ns,1:np)
            local_tot = SUM(poy * Sbr(0,1:ns,1:np)) * 0.25_d/4.0_d/PI
        END IF
        IF ((.NOT.(inner_flag)).AND.(outer)) THEN
            poy(1:ns,1:np) = esbpg(nr,0:ns-1,0:np-1) + esbpg(nr,0:ns-1,1:np) &
                + esbpg(nr,1:ns,0:np-1) + esbpg(nr,1:ns,1:np) &
                - epbsg(nr,0:ns-1,0:np-1) - epbsg(nr,0:ns-1,1:np) &
                - epbsg(nr,1:ns,0:np-1) - epbsg(nr,1:ns,1:np)
            local_tot = SUM(poy * Sbr(nr,1:ns,1:np)) * 0.25_d/4.0_d/PI
        END IF
        
        CALL MPI_ALLREDUCE(local_tot, poynting_flux, 1, mpireal, &
             MPI_SUM, comm, err)
             
        poynting_flux = poynting_flux * RSUN**3

        ! ! Check that we get -1/(4*pi) * int v_out*Bh^2 dS on outer bdy:
        ! IF (.NOT.(inner_flag)) THEN
        !     poy = 0.0_d
        !     local_tot = 0.0_d
        !     IF (outer) THEN
        !         poy(1:ns,1:np) = vr(nr) *  &
        !             (bgs(nr,0:ns-1,0:np-1)**2 + bgp(nr,0:ns-1,0:np-1)**2 &
        !             + bgs(nr,0:ns-1,1:np)**2 + bgp(nr,0:ns-1,1:np)**2 &
        !             + bgs(nr,1:ns,0:np-1)**2 + bgp(nr,1:ns,0:np-1)**2 &
        !             + bgs(nr,1:ns,1:np)**2 + bgp(nr,1:ns,1:np)**2 )
        !         local_tot = SUM(poy * Sbr(nr,1:ns,1:np)) * 0.25_d/4.0_d/PI
        !     ENDIF
        !     CALL MPI_ALLREDUCE(local_tot, vr_flux, 1, mpireal, &
        !     MPI_SUM, comm, err)
        !     vr_flux = -vr_flux * RSUN**3
        !     IF (control) print*, poynting_flux, vr_flux
        !     poynting_flux = vr_flux
        ! END IF


END FUNCTION poynting_flux

!===============================================================================
FUNCTION energy_flux_hyper(inner_flag)
    ! Compute  1/(4*pi) * int alpha * hr dS on inner boundary
    ! (if inner_flag=.TRUE.) or outer boundary (if inner_flag=.FALSE.).
        REAL(d):: energy_flux_hyper     ! units [erg/s]
        LOGICAL, INTENT(IN):: inner_flag ! .TRUE. for inner ,.FALSE. for outer bdy

        ! Local variables:
        REAL(d):: hyp(1:ns,1:np)    ! integrand on horizontal cell faces
        REAL(d):: local_tot

        hyp = 0.0_d
        local_tot = 0.0_d
        IF ((inner_flag).AND.(surface)) THEN
            hyp(1:ns,1:np) = 0.5_d * hyper_eta4 * ( &
                alg(0,0:ns-1,0:np-1) * (hr(0,0:ns-1,0:np-1) + hr(1,0:ns-1,0:np-1)) &
                + alg(0,0:ns-1,1:np) * (hr(0,0:ns-1,1:np) + hr(1,0:ns-1,1:np)) &
                + alg(0,1:ns,0:np-1) * (hr(0,1:ns,0:np-1) + hr(1,1:ns,0:np-1)) &
                + alg(0,1:ns,1:np) * (hr(0,1:ns,1:np) + hr(1,1:ns,1:np)) &
                )
            local_tot = SUM(hyp * Sbr(0,1:ns,1:np)) * 0.125_d/4.0_d/PI
        END IF
        IF (.NOT.(inner_flag).AND.(outer)) THEN
            hyp(1:ns,1:np) = 0.5_d * hyper_eta4 * ( &
                alg(nr,0:ns-1,0:np-1) * (hr(nr,0:ns-1,0:np-1) + hr(nr+1,0:ns-1,0:np-1)) &
                + alg(nr,0:ns-1,1:np) * (hr(nr,0:ns-1,1:np) + hr(nr+1,0:ns-1,1:np)) &
                + alg(nr,1:ns,0:np-1) * (hr(nr,1:ns,0:np-1) + hr(nr+1,1:ns,0:np-1)) &
                + alg(nr,1:ns,1:np) * (hr(nr,1:ns,1:np) + hr(nr+1,1:ns,1:np)) &
                )
            local_tot = SUM(hyp * Sbr(nr,1:ns,1:np)) * 0.125_d/4.0_d/PI
        END IF

        CALL MPI_ALLREDUCE(local_tot, energy_flux_hyper, 1, mpireal, &
        MPI_SUM, comm, err)
        
        energy_flux_hyper = energy_flux_hyper * RSUN**3

END FUNCTION energy_flux_hyper

!===============================================================================
FUNCTION dissipation_friction()
    ! Compute -1/(4*pi) * int |J x B|^2/(nu * B^2) dV.
    ! i.e. volume dissipation of magnetic energy due to friction.
    REAL(d):: dissipation_friction     ! units [erg/s]

    ! Local variables:
    REAL(d):: local_tot
    REAL(d):: fg(0:nr,0:ns,0:np)    ! integrand at grid pts.

    fg(0:nr,0:ns,0:np) = nu(0:nr,0:ns,0:np) * (fgr(0:nr,0:ns,0:np)**2 + &
         fgs(0:nr,0:ns,0:np)**2 + fgp(0:nr,0:ns,0:np)**2) &
         /bbgs(0:nr,0:ns,0:np)

    local_tot = 0.125_d * SUM( &
        fg(0:nr-1,0:ns-1,0:np-1) * Vg(0:nr-1,0:ns-1,0:np-1) &
        + fg(0:nr-1,0:ns-1,1:np) * Vg(0:nr-1,0:ns-1,1:np) &
        + fg(0:nr-1,1:ns,0:np-1) * Vg(0:nr-1,1:ns,0:np-1)  &
        + fg(0:nr-1,1:ns,1:np) * Vg(0:nr-1,1:ns,1:np) &
        + fg(1:nr,0:ns-1,0:np-1) * Vg(1:nr,0:ns-1,0:np-1) &
        + fg(1:nr,0:ns-1,1:np) * Vg(1:nr,0:ns-1,1:np) &
        + fg(1:nr,1:ns,0:np-1) * Vg(1:nr,1:ns,0:np-1) &
        + fg(1:nr,1:ns,1:np) * Vg(1:nr,1:ns,1:np))
    CALL MPI_ALLREDUCE(local_tot, dissipation_friction, 1, mpireal, MPI_SUM, &
        comm, err)

    dissipation_friction = -dissipation_friction * RSUN**3 / (4.0_d * PI)

END FUNCTION dissipation_friction

!===============================================================================
FUNCTION workdone_outflow()
    ! Compute -1/(4*pi) * int v_out e_r.(J x B) dV.
    ! i.e. work done by outflow against Lorentz force in volume.
    REAL(d):: workdone_outflow     ! units [erg/s]

    ! Local variables:
    REAL(d):: local_tot
    REAL(d):: fg(0:nr,0:ns,0:np)    ! integrand at grid pts.
    INTEGER:: i

    DO i=0, nr
        fg(i,0:ns,0:np) = vr(i) * fgr(i,0:ns,0:np)
    END DO

    local_tot = 0.125_d * SUM( &
    fg(0:nr-1,0:ns-1,0:np-1) * Vg(0:nr-1,0:ns-1,0:np-1) &
    + fg(0:nr-1,0:ns-1,1:np) * Vg(0:nr-1,0:ns-1,1:np) &
    + fg(0:nr-1,1:ns,0:np-1) * Vg(0:nr-1,1:ns,0:np-1)  &
    + fg(0:nr-1,1:ns,1:np) * Vg(0:nr-1,1:ns,1:np) &
    + fg(1:nr,0:ns-1,0:np-1) * Vg(1:nr,0:ns-1,0:np-1) &
    + fg(1:nr,0:ns-1,1:np) * Vg(1:nr,0:ns-1,1:np) &
    + fg(1:nr,1:ns,0:np-1) * Vg(1:nr,1:ns,0:np-1) &
    + fg(1:nr,1:ns,1:np) * Vg(1:nr,1:ns,1:np))
    CALL MPI_ALLREDUCE(local_tot, workdone_outflow, 1, mpireal, MPI_SUM, &
    comm, err)

    workdone_outflow = -workdone_outflow * RSUN**3 / (4.0_d * PI)

END FUNCTION workdone_outflow

!===============================================================================
FUNCTION dissipation_hyper()
    ! Compute -1/(4*pi) * int eta4 * B^2 * |grad(alpha)|^2 dV.
    ! i.e. volume dissipation of magnetic energy due to hyperdiffusion.
    REAL(d):: dissipation_hyper     ! units [erg/s]

    ! Local variables:
    REAL(d):: local_tot
    REAL(d):: fg(0:nr,0:ns,0:np)    ! integrand at grid pts.

    fg(0:nr,0:ns,0:np) = 0.0625_d * hyper_eta4 * ( &
        (hr(0:nr,0:ns,0:np) + hr(1:nr+1,0:ns,0:np))**2 &
        + (hs(0:nr,0:ns,0:np) + hs(0:nr,1:ns+1,0:np))**2 &
        + (hp(0:nr,0:ns,0:np) + hp(0:nr,0:ns,1:np+1))**2 ) &
        / bbgs(0:nr,0:ns,0:np)

    local_tot = 0.125_d * SUM( &
        fg(0:nr-1,0:ns-1,0:np-1) * Vg(0:nr-1,0:ns-1,0:np-1) &
        + fg(0:nr-1,0:ns-1,1:np) * Vg(0:nr-1,0:ns-1,1:np) &
        + fg(0:nr-1,1:ns,0:np-1) * Vg(0:nr-1,1:ns,0:np-1)  &
        + fg(0:nr-1,1:ns,1:np) * Vg(0:nr-1,1:ns,1:np) &
        + fg(1:nr,0:ns-1,0:np-1) * Vg(1:nr,0:ns-1,0:np-1) &
        + fg(1:nr,0:ns-1,1:np) * Vg(1:nr,0:ns-1,1:np) &
        + fg(1:nr,1:ns,0:np-1) * Vg(1:nr,1:ns,0:np-1) &
        + fg(1:nr,1:ns,1:np) * Vg(1:nr,1:ns,1:np))

    CALL MPI_ALLREDUCE(local_tot, dissipation_hyper, 1, mpireal, MPI_SUM, &
        comm, err)

    dissipation_hyper = -dissipation_hyper * RSUN**3 / (4.0_d * PI)

END FUNCTION dissipation_hyper

!===============================================================================
SUBROUTINE test_diagnostics(br, bs, bp, bbg, test_case)
! Test (some of) the diagnostics against exact solution.
    REAL(d), INTENT(IN):: br(0:nr,0:ns+1,0:np+1)   ! B_rho on faces
    REAL(d), INTENT(IN):: bs(0:nr+1,0:ns,0:np+1)   ! B_s on faces
    REAL(d), INTENT(IN):: bp(0:nr+1,0:ns+1,0:np)   ! B_phi on faces
    REAL(d), INTENT(IN):: bbg(0:nr,0:ns,0:np)   ! B^2 at grid pts
    INTEGER, INTENT(IN):: test_case
    
    ! Local variables:
    REAL(d):: flux_inner0, flux_outer0, energy0
    REAL(d):: flux_inner, flux_outer, energy

    ! Test computation of reference vector potential:
    CALL zero_aref_outer()
    CALL compute_aref_outer(br)
    IF (control) PRINT*, 'Computed aref on outer boundary.'

    CALL diagnostics_exact(flux_inner0, flux_outer0, energy0, test_case)
    
    energy = magnetic_energy(bbg)
    flux_inner = unsigned_flux_inner(br)
    flux_outer = unsigned_flux_outer(br)

    IF (control) PRINT*, 'energy: ', energy, ' [exact: ', &
        energy0, ']'
    IF (control) PRINT*, 'flux_inner: ', flux_inner, ' [exact: ', &
        flux_inner0, ']'
    IF (control) PRINT*, 'flux_outer: ', flux_outer, ' [exact: ', &
        flux_outer0, ']'

END SUBROUTINE test_diagnostics

!===============================================================================
SUBROUTINE diagnostics_exact(flux_inner0, flux_outer0, energy0, test_case)
! Return exact diagnostics for test cases.
    INTEGER, INTENT(IN):: test_case
    REAL(d), INTENT(INOUT):: flux_inner0, flux_outer0, energy0

    IF (test_case == 0) THEN    ! axisymmetric PFSS dipole
        energy0 = 2.0_d * (rss**3 - 1.0) / 3.0_d / &
            (1 + 2 * rss**3) * 6.96e10_d**3
        flux_inner0 = 4.0_d * PI * 6.96e10**2
        flux_outer0 = 12.0_d * PI * rss**2 / &
            (1.0_d + 2.0_d * rss**3) * 6.96e10**2
    END IF

END SUBROUTINE diagnostics_exact

!*******************************************************************************
END MODULE diagnostics
!*******************************************************************************
