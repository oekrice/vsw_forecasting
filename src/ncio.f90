!*******************************************************************************
MODULE ncio
!*******************************************************************************
! This module implements the magnetic vector potential.
!*******************************************************************************
    USE netcdf
    USE params
    USE mpitools
    USE grid, ONLY: nr, ns, np, nrg, nsg, npg, rg, rc, sg, sc, pg, pc

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: ncdf_open_file, ncdf_close_file
    PUBLIC:: ncdf_read_3darray, ncdf_read_2darray_rs, ncdf_read_2darray_sp
    PUBLIC:: ncdf_read_1darray_r
    PUBLIC:: ncdf_read_surface_1darray_s, ncdf_read_surface_scalar
    PUBLIC:: ncdf_init_snapshot_file, ncdf_write_3darray
    PUBLIC:: ncdf_create_diagnostics_file, ncdf_restart_diagnostics_file
    PUBLIC:: ncdf_write_diagnostic_entry

!*******************************************************************************
    INTEGER, PARAMETER:: nfreal = NF90_DOUBLE  ! type of output variables

!*******************************************************************************
CONTAINS
!*******************************************************************************

!===============================================================================
SUBROUTINE try(status)
! Catch error in reading netcdf fild.
    INTEGER, INTENT(IN):: status
    
    IF (status /= NF90_noerr) THEN
        PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
        CALL finish
    END IF
    
END SUBROUTINE try

!===============================================================================
FUNCTION ncdf_open_file(netcdf_file)
! Open netcdf file and return id.
    INTEGER:: ncdf_open_file
    CHARACTER(*), INTENT(IN):: netcdf_file     ! netcdf file

    CALL try(NF90_OPEN(TRIM(datadir)//'/'//netcdf_file, NF90_WRITE, &
        ncdf_open_file))

END FUNCTION ncdf_open_file

!===============================================================================
SUBROUTINE ncdf_close_file(nc_id)
! Open netcdf file.
    INTEGER, INTENT(IN):: nc_id     ! id of open netcdf file

    CALL try(NF90_CLOSE(nc_id))

END SUBROUTINE ncdf_close_file

!===============================================================================
SUBROUTINE ncdf_read_3darray(netcdf_file, var_name, arr)
! Read local portion of given 3d array from netcdf file.
    CHARACTER(*), INTENT(IN):: netcdf_file     ! netcdf file containing array arr
    CHARACTER(*), INTENT(IN):: var_name
    REAL(d), INTENT(INOUT):: arr(:,:,:)

    ! Local variables:
    INTEGER:: i, nc_id, var_id

    DO i = 0, n_procs - 1
        IF (rank == i) THEN
            CALL try(NF90_OPEN(TRIM(datadir)//'/'//netcdf_file, NF90_NOWRITE, nc_id))
            CALL try(NF90_INQ_VARID(nc_id, var_name, var_id))
            CALL try(NF90_GET_VAR(nc_id, var_id, arr, &
                start=(/mpi_loc(1) * nr + 1, mpi_loc(2) * ns + 1, &
                    mpi_loc(3) * np + 1/), &
                count=SHAPE(arr)))
            CALL try(NF90_CLOSE(nc_id))
        END IF
        CALL MPI_BARRIER(comm, err)
    END DO
    
END SUBROUTINE ncdf_read_3darray

!===============================================================================
SUBROUTINE ncdf_read_2darray_rs(netcdf_file, var_name, arr)
! Read local portion of given 2d array in (rho, s) from netcdf file.
    CHARACTER(*), INTENT(IN):: netcdf_file     ! netcdf file containing array arr
    CHARACTER(*), INTENT(IN):: var_name
    REAL(d), INTENT(INOUT):: arr(:,:)

    ! Local variables:
    INTEGER:: i, nc_id, var_id

    DO i = 0, n_procs - 1
        IF (rank == i) THEN
            CALL try(NF90_OPEN(TRIM(datadir)//'/'//netcdf_file, NF90_NOWRITE, &
                nc_id))
            CALL try(NF90_INQ_VARID(nc_id, var_name, var_id))
            CALL try(NF90_GET_VAR(nc_id, var_id, arr, &
                start=(/mpi_loc(1) * nr + 1, mpi_loc(2) * ns + 1/), &
                count=SHAPE(arr)))
            CALL try(NF90_CLOSE(nc_id))
        END IF
        CALL MPI_BARRIER(comm, err)
    END DO

END SUBROUTINE ncdf_read_2darray_rs

!===============================================================================
SUBROUTINE ncdf_read_2darray_sp(netcdf_file, var_name, arr)
! Read local portion of given 2d array in (s, phi) from netcdf file.
    CHARACTER(*), INTENT(IN):: netcdf_file     ! netcdf file containing array arr
    CHARACTER(*), INTENT(IN):: var_name
    REAL(d), INTENT(INOUT):: arr(:,:)

    ! Local variables:
    INTEGER:: i, nc_id, var_id

    DO i = 0, n_procs - 1
        IF (rank == i) THEN
            CALL try(NF90_OPEN(TRIM(datadir)//'/'//netcdf_file, NF90_NOWRITE, &
                nc_id))
            CALL try(NF90_INQ_VARID(nc_id, var_name, var_id))
            CALL try(NF90_GET_VAR(nc_id, var_id, arr, &
                start=(/mpi_loc(2) * ns + 1, mpi_loc(3) * np + 1/), &
                count=SHAPE(arr)))
            CALL try(NF90_CLOSE(nc_id))
        END IF
        CALL MPI_BARRIER(comm, err)
    END DO

END SUBROUTINE ncdf_read_2darray_sp

!===============================================================================
SUBROUTINE ncdf_read_1darray_r(netcdf_file, var_name, arr)
! Read local portion of given 1d array in rho from netcdf file.
    CHARACTER(*), INTENT(IN):: netcdf_file     ! netcdf file containing array arr
    CHARACTER(*), INTENT(IN):: var_name
    REAL(d), INTENT(INOUT):: arr(:)

    ! Local variables:
    INTEGER:: i, nc_id, var_id

    DO i = 0, n_procs - 1
        IF (rank == i) THEN
            CALL try(NF90_OPEN(TRIM(datadir)//'/'//netcdf_file, NF90_NOWRITE, &
                nc_id))
            CALL try(NF90_INQ_VARID(nc_id, var_name, var_id))
            CALL try(NF90_GET_VAR(nc_id, var_id, arr, &
                start=(/mpi_loc(1) * nr + 1/), count=SHAPE(arr)))
            CALL try(NF90_CLOSE(nc_id))
        END IF
        CALL MPI_BARRIER(comm, err)
    END DO

END SUBROUTINE ncdf_read_1darray_r

!===============================================================================
SUBROUTINE ncdf_read_surface_1darray_s(netcdf_file, var_name, arr)
! Read local portion of given 1d array in s from netcdf file.
    CHARACTER(*), INTENT(IN):: netcdf_file     ! netcdf file containing array arr
    CHARACTER(*), INTENT(IN):: var_name
    REAL(d), INTENT(INOUT):: arr(:)

    ! Local variables:
    INTEGER:: i, nc_id, var_id

    DO i = 0, n_procs - 1
        IF ((rank == i) .AND. surface) THEN
            CALL try(NF90_OPEN(TRIM(datadir)//'/'//netcdf_file, NF90_NOWRITE, &
                nc_id))
            CALL try(NF90_INQ_VARID(nc_id, var_name, var_id))
            CALL try(NF90_GET_VAR(nc_id, var_id, arr, &
                start=(/mpi_loc(2) * ns + 1/), count=SHAPE(arr)))
            CALL try(NF90_CLOSE(nc_id))
        END IF
        CALL MPI_BARRIER(comm, err)
    END DO

END SUBROUTINE ncdf_read_surface_1darray_s

!===============================================================================
SUBROUTINE ncdf_read_surface_scalar(netcdf_file, var_name, val)
! Read scalar from netcdf file.
    CHARACTER(*), INTENT(IN):: netcdf_file     ! netcdf file containing val
    CHARACTER(*), INTENT(IN):: var_name
    REAL(d), INTENT(INOUT):: val

    ! Local variables:
    INTEGER:: i, nc_id, var_id

    DO i = 0, n_procs - 1
        IF ((rank == i) .AND. surface) THEN
            CALL try(NF90_OPEN(TRIM(datadir)//'/'//netcdf_file, NF90_NOWRITE, &
                nc_id))
            CALL try(NF90_INQ_VARID(nc_id, var_name, var_id))
            CALL try(NF90_GET_VAR(nc_id, var_id, val))
            CALL try(NF90_CLOSE(nc_id))
        END IF
        CALL MPI_BARRIER(comm, err)
    END DO

END SUBROUTINE ncdf_read_surface_scalar

!===============================================================================
SUBROUTINE ncdf_init_snapshot_file(netcdf_file, centres_flag)
! Create (or overwrite) netcdf file for output snapshot, create dimensions and
! write grid arrays.
    CHARACTER(*), INTENT(IN):: netcdf_file     ! netcdf file to create
    LOGICAL, INTENT(IN):: centres_flag  ! whether to include cell centre arrays.
    
    ! Local variables:
    INTEGER:: nc_id, r_dim_id, rc_dim_id, s_dim_id, sc_dim_id
    INTEGER:: p_dim_id, pc_dim_id
    
    IF (control) THEN
        CALL try(NF90_CREATE(TRIM(datadir)//'/'//netcdf_file, NF90_CLOBBER, &
            nc_id))
        CALL try(NF90_DEF_DIM(nc_id, 'r', nrg+1, r_dim_id))
        IF (centres_flag) CALL try(NF90_DEF_DIM(nc_id, 'rc', nrg+2, rc_dim_id))
        call try(NF90_DEF_DIM(nc_id, 'th', nsg+1, s_dim_id))
        IF (centres_flag) CALL try(NF90_DEF_DIM(nc_id, 'thc', nsg+2, sc_dim_id))
        call try(NF90_DEF_DIM(nc_id, 'ph', npg+1, p_dim_id))
        IF (centres_flag) CALL try(NF90_DEF_DIM(nc_id, 'phc', npg+2, pc_dim_id))
        call try(NF90_ENDDEF(nc_id))
        call try(NF90_CLOSE(nc_id))
    END IF
    
    ! Write 1d grid arrays:
    CALL ncdf_write_1darray_r(netcdf_file, 'r', 'r', EXP(rg(0:nr)))
    IF (centres_flag) CALL ncdf_write_1darray_r(netcdf_file, 'rc', 'rc', &
        EXP(rc(0:nr+1)))
    CALL ncdf_write_1darray_s(netcdf_file, 'th', 'th', ACOS(sg(0:ns)))
    IF (centres_flag) CALL ncdf_write_1darray_s(netcdf_file, 'thc', 'thc', &
        ACOS(sc(0:ns+1)))
    CALL ncdf_write_1darray_p(netcdf_file, 'ph', 'ph', pg(0:np))
    IF (centres_flag) CALL ncdf_write_1darray_p(netcdf_file, 'phc', 'phc', &
        pc(0:np+1))
    
END SUBROUTINE ncdf_init_snapshot_file

!===============================================================================
SUBROUTINE ncdf_write_3darray(netcdf_file, var_name, r_dim, s_dim, p_dim, arr)
! Write 3d array to netcdf file. File should already exist.
    CHARACTER(*), INTENT(IN):: netcdf_file     ! netcdf file
    CHARACTER(*), INTENT(IN):: var_name    ! variable name (will be created)
    CHARACTER(*), INTENT(IN):: r_dim, s_dim, p_dim  ! names of dimensions
                                                    ! (should already exist)
    REAL(d), INTENT(IN):: arr(:,:,:)

    ! Local variables:
    INTEGER:: i, nc_id, r_dim_id, s_dim_id, p_dim_id, var_id

    ! Create the variable:
    IF (control) THEN
        CALL try(NF90_OPEN(TRIM(datadir)//'/'//netcdf_file, NF90_WRITE, nc_id))
        CALL try(NF90_INQ_DIMID(nc_id, r_dim, r_dim_id))
        CALL try(NF90_INQ_DIMID(nc_id, s_dim, s_dim_id))
        CALL try(NF90_INQ_DIMID(nc_id, p_dim, p_dim_id))
        CALL try(NF90_REDEF(nc_id))
        CALL try(NF90_DEF_VAR(nc_id, var_name, nfreal, (/r_dim_id, s_dim_id, &
            p_dim_id/), var_id))
        CALL try(NF90_ENDDEF(nc_id))
        CALL try(NF90_CLOSE(nc_id))
    END IF

    CALL MPI_BARRIER(comm, err)
    
    ! Each MPI process writes local data to file in turn:
    DO i = 0, n_procs - 1
        IF (rank == i) THEN
            CALL try(NF90_OPEN(TRIM(datadir)//'/'//netcdf_file, NF90_WRITE, &
                nc_id))
            CALL try(NF90_INQ_VARID(nc_id, var_name, var_id))
            CALL try(NF90_PUT_VAR(nc_id, var_id, arr, &
                start=(/mpi_loc(1) * nr + 1, mpi_loc(2) * ns + 1, &
                mpi_loc(3) * np + 1/), count=SHAPE(arr)))
            CALL try(NF90_CLOSE(nc_id))
        END IF
        CALL MPI_BARRIER(comm, err)
    END DO
    
END SUBROUTINE ncdf_write_3darray

!===============================================================================
SUBROUTINE ncdf_write_1darray_r(netcdf_file, var_name, r_dim, arr)
! Write 1d array in rho to netcdf file. File should already exist.
    CHARACTER(*), INTENT(IN):: netcdf_file     ! netcdf file
    CHARACTER(*), INTENT(IN):: var_name    ! variable name (will be created)
    CHARACTER(*), INTENT(IN):: r_dim  ! names of dimension (should already exist)
    REAL(d), INTENT(IN):: arr(:)
    
    ! Local variables:
    INTEGER:: i, nc_id, r_dim_id, var_id

    ! Create the variable:
    IF (control) THEN
        CALL try(NF90_OPEN(TRIM(datadir)//'/'//netcdf_file, NF90_WRITE, nc_id))
        CALL try(NF90_INQ_DIMID(nc_id, r_dim, r_dim_id))
        CALL try(NF90_REDEF(nc_id))
        CALL try(NF90_DEF_VAR(nc_id, var_name, nfreal, (/r_dim_id/), var_id))
        CALL try(NF90_ENDDEF(nc_id))
        CALL try(NF90_CLOSE(nc_id))
    END IF

    CALL MPI_BARRIER(comm, err)
    
    ! Each MPI process writes local data to file in turn:
    DO i = 0, n_procs - 1
        IF (rank == i) THEN
            CALL try(NF90_OPEN(TRIM(datadir)//'/'//netcdf_file, NF90_WRITE, &
                nc_id))
            CALL try(NF90_INQ_VARID(nc_id, var_name, var_id))
            CALL try(NF90_PUT_VAR(nc_id, var_id, arr, &
                start=(/mpi_loc(1) * nr + 1/), count=SHAPE(arr)))
            CALL try(NF90_CLOSE(nc_id))
        END IF
        CALL MPI_BARRIER(comm, err)
    END DO
    
END SUBROUTINE ncdf_write_1darray_r

!===============================================================================
SUBROUTINE ncdf_write_1darray_s(netcdf_file, var_name, s_dim, arr)
! Write 1d array in s to netcdf file. File should already exist.
    CHARACTER(*), INTENT(IN):: netcdf_file     ! netcdf file
    CHARACTER(*), INTENT(IN):: var_name    ! variable name (will be created)
    CHARACTER(*), INTENT(IN):: s_dim  ! names of dimension (should already exist)
    REAL(d), INTENT(IN):: arr(:)
    
    ! Local variables:
    INTEGER:: i, nc_id, s_dim_id, var_id

    ! Create the variable:
    IF (control) THEN
        CALL try(NF90_OPEN(TRIM(datadir)//'/'//netcdf_file, NF90_WRITE, nc_id))
        CALL try(NF90_INQ_DIMID(nc_id, s_dim, s_dim_id))
        CALL try(NF90_REDEF(nc_id))
        CALL try(NF90_DEF_VAR(nc_id, var_name, nfreal, (/s_dim_id/), var_id))
        CALL try(NF90_ENDDEF(nc_id))
        CALL try(NF90_CLOSE(nc_id))
    END IF

    CALL MPI_BARRIER(comm, err)
    
    ! Each MPI process writes local data to file in turn:
    DO i = 0, n_procs - 1
        IF (rank == i) THEN
            CALL try(NF90_OPEN(TRIM(datadir)//'/'//netcdf_file, NF90_WRITE, &
                nc_id))
            CALL try(NF90_INQ_VARID(nc_id, var_name, var_id))
            CALL try(NF90_PUT_VAR(nc_id, var_id, arr, &
                start=(/mpi_loc(2) * ns + 1/), count=SHAPE(arr)))
            CALL try(NF90_CLOSE(nc_id))
        END IF
        CALL MPI_BARRIER(comm, err)
    END DO
    
END SUBROUTINE ncdf_write_1darray_s

!===============================================================================
SUBROUTINE ncdf_write_1darray_p(netcdf_file, var_name, p_dim, arr)
! Write 1d array in phi to netcdf file. File should already exist.
    CHARACTER(*), INTENT(IN):: netcdf_file     ! netcdf file
    CHARACTER(*), INTENT(IN):: var_name    ! variable name (will be created)
    CHARACTER(*), INTENT(IN):: p_dim  ! names of dimension (should already exist)
    REAL(d), INTENT(IN):: arr(:)
    
    ! Local variables:
    INTEGER:: i, nc_id, p_dim_id, var_id

    ! Create the variable:
    IF (control) THEN
        CALL try(NF90_OPEN(TRIM(datadir)//'/'//netcdf_file, NF90_WRITE, nc_id))
        CALL try(NF90_INQ_DIMID(nc_id, p_dim, p_dim_id))
        CALL try(NF90_REDEF(nc_id))
        CALL try(NF90_DEF_VAR(nc_id, var_name, nfreal, (/p_dim_id/), var_id))
        CALL try(NF90_ENDDEF(nc_id))
        CALL try(NF90_CLOSE(nc_id))
    END IF

    CALL MPI_BARRIER(comm, err)
    
    ! Each MPI process writes local data to file in turn:
    DO i = 0, n_procs - 1
        IF (rank == i) THEN
            CALL try(NF90_OPEN(TRIM(datadir)//'/'//netcdf_file, NF90_WRITE, &
                nc_id))
            CALL try(NF90_INQ_VARID(nc_id, var_name, var_id))
            CALL try(NF90_PUT_VAR(nc_id, var_id, arr, &
                start=(/mpi_loc(3) * np + 1/), count=SHAPE(arr)))
            CALL try(NF90_CLOSE(nc_id))
        END IF
        CALL MPI_BARRIER(comm, err)
    END DO
    
END SUBROUTINE ncdf_write_1darray_p

!===============================================================================
SUBROUTINE ncdf_create_diagnostics_file(diag_file, start_time, end_time, &
    cadence_em, unsigned_flux_inner_flag, unsigned_flux_outer_flag, &
    magnetic_energy_flag, average_current_flag, outer_bh_flag, &
    helicity_flux_flag, energy_terms_flag, eruption_diagnostics_flag, n_caps)
! Create diagnostics file, write attributes and initialize required variables.
    CHARACTER(*), INTENT(IN):: diag_file    ! diagnostics file
    CHARACTER(*), INTENT(IN):: start_time, end_time
    INTEGER, INTENT(IN):: cadence_em
    LOGICAL, INTENT(IN):: unsigned_flux_inner_flag, unsigned_flux_outer_flag
    LOGICAL, INTENT(IN):: magnetic_energy_flag, average_current_flag
    LOGICAL, INTENT(IN):: outer_bh_flag, helicity_flux_flag
    LOGICAL, INTENT(IN):: energy_terms_flag, eruption_diagnostics_flag
    INTEGER, INTENT(IN):: n_caps

    ! Local variables:
    INTEGER:: nc_id, t_dim_id, t_id, f_id, o_id, e_id, j_id, bsn_id, bss_id
    INTEGER:: bpn_id, bps_id, hf_id, poyi_id, poyo_id, mt_id, npts_id
    INTEGER:: ehypi_id, ehypo_id, dissf_id, work_id, dissh_id, i
    CHARACTER(40):: cap_varname

    CALL try(NF90_CREATE(TRIM(datadir)//'/'//diag_file, NF90_CLOBBER, nc_id))
    CALL try(NF90_PUT_ATT(nc_id, NF90_GLOBAL, 'start_time', start_time))
    CALL try(NF90_PUT_ATT(nc_id, NF90_GLOBAL, 'end_time', end_time))
    CALL try(NF90_PUT_ATT(nc_id, NF90_GLOBAL, 'cadence_em', cadence_em))
    CALL try(NF90_DEF_DIM(nc_id, 't', NF90_UNLIMITED , t_dim_id))
    CALL try(NF90_DEF_VAR(nc_id, 't', nfreal, t_dim_id, t_id))
    IF (unsigned_flux_inner_flag) CALL try(NF90_DEF_VAR(nc_id, &
        'phot_flux', nfreal, t_dim_id, f_id))
    IF (unsigned_flux_outer_flag) CALL try(NF90_DEF_VAR(nc_id, &
        'open_flux', nfreal, t_dim_id, o_id))
    IF (magnetic_energy_flag) CALL try(NF90_DEF_VAR(nc_id, 'energy', &
        nfreal, t_dim_id, e_id))
    IF (average_current_flag) CALL try(NF90_DEF_VAR(nc_id, 'mean_current', &
        nfreal, t_dim_id, j_id))
    IF (outer_bh_flag) THEN
        CALL try(NF90_DEF_VAR(nc_id, 'bs_nor', nfreal, t_dim_id, bsn_id))
        CALL try(NF90_DEF_VAR(nc_id, 'bs_sou', nfreal, t_dim_id, bss_id))
        CALL try(NF90_DEF_VAR(nc_id, 'bp_nor', nfreal, t_dim_id, bpn_id))
        CALL try(NF90_DEF_VAR(nc_id, 'bp_sou', nfreal, t_dim_id, bps_id))
    END IF
    IF (helicity_flux_flag) THEN
        CALL try(NF90_DEF_VAR(nc_id, 'hr_flux', &
            nfreal, t_dim_id, hf_id))
        DO i=1,n_caps
            write(cap_varname,"('hr_flux_',I0.3)") i
            CALL try(NF90_DEF_VAR(nc_id, TRIM(cap_varname),nfreal, &
                         t_dim_id, hf_id))
        END DO
    END IF
    IF (energy_terms_flag) THEN
        CALL try(NF90_DEF_VAR(nc_id, 'poynting_inner', nfreal, &
            t_dim_id, poyi_id))
        CALL try(NF90_DEF_VAR(nc_id, 'poynting_outer', nfreal, &
            t_dim_id, poyo_id))
        CALL try(NF90_DEF_VAR(nc_id, 'energy_flux_hyper_inner', nfreal, &
            t_dim_id, ehypi_id))
        CALL try(NF90_DEF_VAR(nc_id, 'energy_flux_hyper_outer', nfreal, &
            t_dim_id, ehypo_id))
        CALL try(NF90_DEF_VAR(nc_id, 'dissipation_friction', nfreal, &
        t_dim_id, dissf_id))
        CALL try(NF90_DEF_VAR(nc_id, 'work_outflow', nfreal, &
        t_dim_id, work_id))
        CALL try(NF90_DEF_VAR(nc_id, 'dissipation_hyper', nfreal, &
        t_dim_id, dissh_id))
    END IF
    IF (eruption_diagnostics_flag) THEN
        CALL try(NF90_DEF_VAR(nc_id, 'max_tension', nfreal, t_dim_id, mt_id))
        CALL try(NF90_DEF_VAR(nc_id, 'npoints_bh', nfreal, t_dim_id, npts_id))
    END IF
    CALL try(NF90_ENDDEF(nc_id))
    CALL try(NF90_CLOSE(nc_id))
    
END SUBROUTINE ncdf_create_diagnostics_file

!===============================================================================
SUBROUTINE ncdf_restart_diagnostics_file(diag_file, t, position)
! Identify last position in file before current time, and start overwriting
! after that.
    CHARACTER(*), INTENT(IN):: diag_file    ! diagnostics file
    REAL(d), INTENT(IN):: t     ! time that we are restarting from
    INTEGER, INTENT(INOUT):: position     ! position in file
    
    ! Local variables:
    INTEGER:: nc_id, t_dim_id, t_id, nt
    CHARACTER*(1):: nm
    REAL(d), ALLOCATABLE:: tf(:)
    
    CALL try(NF90_OPEN(TRIM(datadir)//'/'//diag_file, NF90_NOWRITE, nc_id))
    CALL try(NF90_INQ_DIMID(nc_id, 't', t_dim_id))
    CALL try(NF90_INQUIRE_DIMENSION(nc_id, t_dim_id, nm, nt))
    ALLOCATE(tf(nt))
    CALL try(NF90_INQ_VARID(nc_id, 't', t_id))
    CALL try(NF90_GET_VAR(nc_id, t_id, tf))
    CALL try(NF90_CLOSE(nc_id))
    DO WHILE(position < nt)
        IF (tf(position) > t) EXIT
        position = position + 1
    END DO
    DEALLOCATE(tf)
    
END SUBROUTINE ncdf_restart_diagnostics_file

!===============================================================================
SUBROUTINE ncdf_write_diagnostic_entry(nc_id, var_name, val, position)
! Write individual entry to diagnostic file.
    INTEGER, INTENT(IN):: nc_id     ! id of already open netcdf file
    CHARACTER(*), INTENT(IN):: var_name    ! variable name
    REAL(d), INTENT(IN):: val   ! value to write
    INTEGER, INTENT(IN):: position  ! position in file
    
    ! Local variables:
    INTEGER:: var_id

    CALL try(NF90_INQ_VARID(nc_id, var_name, var_id))
    CALL try(NF90_PUT_VAR(nc_id, var_id, val, start=(/position/)))

END SUBROUTINE ncdf_write_diagnostic_entry

!*******************************************************************************
END MODULE ncio
!*******************************************************************************
