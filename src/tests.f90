!*******************************************************************************
MODULE tests
!*******************************************************************************
! This module implements unit tests targeted at individual subroutines.
!
! The routines here are designed to be run in isolation so as to bypass the
! usual code initialisation. This ensures they are isolating the effect of
! particular subroutines.
!*******************************************************************************
    USE params
    USE mpitools, ONLY: start_mpi, subdivide_mpi, finish, control
    USE grid, ONLY: nrg, nsg, npg, rss, nr, ns, np, create_grid_arrays !, &
    USE vecpot, ONLY: alr, als, alp, set_a_exact
    USE bfield, ONLY: br, bs, bp, bbg, test_compute_b_stokes, &
        set_b_exact
    USE diagnostics, ONLY: test_diagnostics
    USE current, ONLY: test_compute_j_stokes

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: unit_tests

!*******************************************************************************
    INTEGER, PARAMETER:: test_case=0       ! 0 for PFSS-dipole
    LOGICAL, PARAMETER:: output_netcdf=.TRUE.

!*******************************************************************************
CONTAINS
!*******************************************************************************

!===============================================================================
SUBROUTINE unit_tests()
! Read 2d profile of friction coefficient in (rho, s) and set 3d array.

    ! Local variables:
    CHARACTER(str_max_len):: command_line

    CALL GET_COMMAND_ARGUMENT(1, command_line)
    IF ((command_line /= '').AND.(command_line /= '-test')) RETURN

    datadir = './'

    CALL start_mpi

    IF (control) PRINT*, '***********************************************'
    IF (control) PRINT*, '********  RUNNING UNIT TESTS ONLY...  *********'
    IF (control) PRINT*, '***********************************************'

    ! Set up grid:
    nrg = 30
    nsg = 90
    npg = 180
    CALL subdivide_mpi(nrg, nsg, npg, nr, ns, np)
    rss = 2.5_d
    CALL create_grid_arrays
    
    ! Test computation of B from exact A:
    CALL set_a_exact(test_case)
    CALL test_compute_b_stokes(alr, als, alp, test_case, output_netcdf)

    ! Test computation of J from exact B:
    CALL set_b_exact(test_case)
    CALL test_compute_j_stokes(br, bs, bp, test_case, output_netcdf)

    ! Test diagnostics:
    CALL set_b_exact(test_case)
    CALL test_diagnostics(br, bs, bp, bbg, test_case)

    CALL finish

END SUBROUTINE unit_tests

!*******************************************************************************
END MODULE tests
!*******************************************************************************
