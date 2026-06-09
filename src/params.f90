!*******************************************************************************
MODULE params
!*******************************************************************************
! Declare global parameters.
!*******************************************************************************

    IMPLICIT NONE
    
!*******************************************************************************
    INTEGER, PARAMETER:: d = KIND(0.d0) ! precision for floats
    REAL(d), PARAMETER:: PI = 4.0_d * atan(1.0_d)   ! pi

    LOGICAL:: restart   ! flag for whether simulation has been restarted

    INTEGER, PARAMETER:: str_max_len=80  ! max length of strings for i/o
    CHARACTER(str_max_len):: datadir     ! path to data directory
    
    INTEGER:: snapshot_cadence  ! cadence to output snaphots [hrs]
    CHARACTER(5):: snapshot_variables   ! which variables to include in snapshot

    REAL(d):: cflfact   ! CFL parameter to control timestep
    
!*******************************************************************************
END MODULE params
!*******************************************************************************
