!*******************************************************************************
MODULE mpitools
!*******************************************************************************
! This module deals with MPI.
!*******************************************************************************
    USE params

    IMPLICIT NONE

    INCLUDE 'mpif.h'
    
    PRIVATE
    
    PUBLIC:: control, surface, outer, south, north
    PUBLIC:: mpireal, err, comm, stat, n_procs, rank, mpi_dims, mpi_loc, mpi_opp
    PUBLIC:: mpi_r_prv, mpi_r_nxt, mpi_s_prv, mpi_s_nxt, mpi_p_prv, mpi_p_nxt
    PUBLIC:: MPI_INTEGER, MPI_CHAR, MPI_MIN, MPI_MAX, MPI_SUM, i1r, i0s, i1s
    PUBLIC:: start_mpi, finish, subdivide_mpi
    
!*******************************************************************************
    LOGICAL:: control
    LOGICAL:: surface
    LOGICAL:: outer
    LOGICAL:: south
    LOGICAL:: north
    INTEGER:: mpireal=MPI_DOUBLE_PRECISION
    INTEGER:: err
    INTEGER:: n_procs
    INTEGER:: comm
    INTEGER:: rank
    INTEGER:: mpi_dims(3)
    INTEGER:: mpi_loc(3)
    INTEGER:: stat(MPI_STATUS_SIZE)
    INTEGER:: mpi_r_prv, mpi_r_nxt
    INTEGER:: mpi_s_prv, mpi_s_nxt
    INTEGER:: mpi_p_prv, mpi_p_nxt
    INTEGER:: mpi_opp
    INTEGER, PARAMETER:: tag=1

    INTEGER:: i1r   ! Radial index of first non-surface layer.
    INTEGER:: i0s, i1s   ! Latitudinal index of first/last interior grid pts.

!*******************************************************************************
CONTAINS
!*******************************************************************************

!===============================================================================
SUBROUTINE start_mpi()
! Initialize MPI.
    
    ! Local variables:
    INTEGER:: tmp_rank

    CALL MPI_INIT(err)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, n_procs, err)
    comm = MPI_COMM_WORLD
    
    ! Determine control process:
    CALL MPI_COMM_RANK(comm, tmp_rank, err)
    control = (tmp_rank == 0)
    
    IF (control) THEN
        PRINT*, &
            '================================================================='
        PRINT*, &
            '==      DUMFRIC-FCAST -- Durham Magneto-frictional code        =='
        PRINT*, &
            '==      in spherical shell using (rho, s, phi) grid.           =='
            PRINT*, &
            '==      [ FORECASTING VERSION ]                                =='            
        PRINT*, &
            '==                                                             =='
        PRINT*, &
            '==       Contact: A.Yeates (anthony.yeates@durham.ac.uk)       =='
        PRINT*, &
            '================================================================='
        PRINT*, &
            'Number of MPI processes:', n_procs
    END IF
    
END SUBROUTINE start_mpi

!===============================================================================
SUBROUTINE finish()
! Terminate MPI and stop the code.

    CALL MPI_FINALIZE(err)
    STOP

END SUBROUTINE finish

!===============================================================================
SUBROUTINE subdivide_mpi(nrg, nsg, npg, nr, ns, np)
! Subdivide the grid between MPI processes.
    INTEGER, INTENT(IN):: nrg   ! global number of grid cells in rho
    INTEGER, INTENT(IN):: nsg   ! global number of grid cells in s
    INTEGER, INTENT(IN):: npg   ! global number of grid cells in phi
    INTEGER, INTENT(OUT):: nr   ! local number of grid cells in rho
    INTEGER, INTENT(OUT):: ns   ! local number of grid cells in s
    INTEGER, INTENT(OUT):: np   ! local number of grid cells in phi

    ! Local variables:
    LOGICAL:: mpi_periodic(3)
    INTEGER:: i
    INTEGER:: mpi_dims_best(3)
    REAL(d):: mean, diff, diff_best
    
    ! Check that nprocs is even
    IF (MOD(n_procs, 2) /= 0) THEN
        IF (control) PRINT*,'ERROR: TOTAL NUMBER OF MPI PROCS MUST BE EVEN'
        CALL finish
    END IF

    ! Choose optimum division of procs that fits grid dimensions:
    diff_best = REAL(n_procs)
    mpi_dims_best = (/-1,-1,-1/)
    DO i = 2, n_procs, 2
        IF (MOD(n_procs, i) == 0) THEN
            mpi_dims = (/0, 0, i/)
            ! Find optimum decomposition with i points in p:
            CALL MPI_DIMS_CREATE(n_procs, 3, mpi_dims, err)
            ! Check whether this is allowed:
            nr = nrg / mpi_dims(1)
            ns = nsg / mpi_dims(2)
            np = npg / mpi_dims(3)
            IF (((nr * mpi_dims(1)) == nrg) &
                .AND. ((ns * mpi_dims(2)) == nsg) &
                .AND. ((np * mpi_dims(3)) == npg)) THEN
                mean = SUM(REAL(mpi_dims))/3.0_d
                diff = MAXVAL(REAL(mpi_dims) - mean)
                IF (diff < diff_best) THEN
                    diff_best = diff
                    mpi_dims_best = mpi_dims
                END IF
            END IF
        END IF
    END DO
    IF (mpi_dims_best(1) * mpi_dims_best(2) * mpi_dims_best(3) == n_procs) THEN
        mpi_dims = mpi_dims_best
        nr = nrg / mpi_dims(1)
        ns = nsg / mpi_dims(2)
        np = npg / mpi_dims(3)
    ELSE
        PRINT*,'ERROR: THIS NUMBER OF MPI PROCS DOES NOT FIT THE GRID'
        CALL finish
    END IF

    mpi_periodic=(/.FALSE.,.FALSE.,.TRUE./)

    CALL MPI_CART_CREATE(MPI_COMM_WORLD, 3, mpi_dims, mpi_periodic, .TRUE., &
        comm, err)

    ! Get my rank:
    CALL MPI_COMM_RANK(comm, rank, err)
    control = (rank == 0)
    
    ! Get coordinates of me and my neighbours:
    CALL MPI_CART_COORDS(comm, rank, 3, mpi_loc, err)
    CALL MPI_CART_SHIFT(comm, 0, 1, mpi_r_prv, mpi_r_nxt, err)
    CALL MPI_CART_SHIFT(comm, 1, 1, mpi_s_prv, mpi_s_nxt, err)
    CALL MPI_CART_SHIFT(comm, 2, 1, mpi_p_prv, mpi_p_nxt, err)

    ! Set surface and outer flags:
    surface = (mpi_loc(1) == 0)
    outer = (mpi_loc(1) == mpi_dims(1) - 1)
    south = (mpi_loc(2) == 0)
    north = (mpi_loc(2) == mpi_dims(2) - 1)

    ! Radial index of first non-surface layer:
    i1r = 0
    IF (surface) i1r = 1

    ! Latitudinal indices of first and last interior grid point:
    i0s = 0
    i1s = ns
    IF (south) i0s = 1
    IF (north) i1s = ns-1

    ! Get rank of opposite process in longitude (for polar BCs):
    CALL MPI_CART_SHIFT(comm, 2, mpi_dims(3) / 2, mpi_opp, i, err)

END SUBROUTINE subdivide_mpi

!*******************************************************************************
END MODULE mpitools
!*******************************************************************************
