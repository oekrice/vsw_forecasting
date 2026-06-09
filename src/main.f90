!*******************************************************************************
PROGRAM main
!*******************************************************************************
! DuMFric main program file. To run, give the user namelist file as a command-line argument. If no command-line argument is given (or "-test") then the code runs the unit tests.
!*******************************************************************************
    USE params
    USE mpitools, ONLY: finish
    USE init, ONLY: initialize
    USE evolve, ONLY: run_simulation
    USE tests, ONLY: unit_tests
    
    IMPLICIT NONE

    CALL unit_tests
    
    CALL initialize

    CALL run_simulation

    CALL finish
    
END PROGRAM main
