!*******************************************************************************
MODULE time
!*******************************************************************************
! This module contains variables and routines for dealing with dates/times
! in the code. These work down to the level of hours.
!*******************************************************************************
    USE params
    USE mpitools, ONLY: control, comm, err, MPI_CHAR

    IMPLICIT NONE
    
    PRIVATE
    
    PUBLIC:: OPERATOR (==)
    PUBLIC:: start_time, end_time, restart_time, time_now, time_end, hours_total
    PUBLIC:: seconds, hours
    PUBLIC:: time_equality, str_to_time, time_to_str, increase_time
    PUBLIC:: elapsed_time, get_restart_time, initialize_time
    
!*******************************************************************************
    CHARACTER(11):: start_time   ! Start time of simulation in "yyyymmdd.hh" fmt.
    CHARACTER(11):: end_time  ! End time of simulation in "yyyymmdd.hh" fmt.
    CHARACTER(11):: restart_time  ! Restart time of simulation in "yyyymmdd.hh" fmt.
    
    REAL(d):: seconds   ! current time counter in seconds
    INTEGER:: hours     ! current time counter in hours
    INTEGER:: hours_total   ! Total length of simulation [hrs]
    
    TYPE :: time_type
        INTEGER:: year
        INTEGER:: month
        INTEGER:: day
        INTEGER:: hour
    END TYPE time_type
    TYPE(time_type):: time_now  ! Current date/time in time_type format.
    TYPE(time_type):: time_end  ! End date/time in time_type format.
    INTERFACE OPERATOR (==)
        MODULE PROCEDURE time_equality
    END INTERFACE

!*******************************************************************************
CONTAINS
!*******************************************************************************

!===============================================================================
FUNCTION time_equality(t1, t2)
! Overload == operator for time_type objects.
    LOGICAL:: time_equality
    TYPE(time_type), INTENT(IN):: t1, t2
    
    time_equality = .FALSE.
    IF ((t1%year == t2%year) .AND. (t1%month == t2%month) &
        .AND. (t1%day == t2%day) .AND. (t1%hour == t2%hour)) &
        time_equality = .TRUE.
    
END FUNCTION time_equality

!===============================================================================
FUNCTION str_to_time(str)
! Create time object from date/time string in the format "yyyymmdd.hh"
    TYPE(time_type):: str_to_time
    CHARACTER(11), INTENT(IN):: str

    READ(str,'(I4.4, I2.2, I2.2, 1X, I2.2)') str_to_time%year, &
        str_to_time%month, str_to_time%day, str_to_time%hour

END FUNCTION str_to_time

!===============================================================================
FUNCTION time_to_str(t1)
    ! Create string in the format "yyyymmdd.hh" from time object.
    CHARACTER(11):: time_to_str
    TYPE(time_type), INTENT(IN):: t1
    
    WRITE(time_to_str,'(I4.4, I2.2, I2.2, A, I2.2)') t1%year, t1%month, &
        t1%day, '.', t1%hour

END FUNCTION time_to_str

!===============================================================================
SUBROUTINE increase_time(t1, dhours)
    ! Advance time object by dhours hours.
    TYPE(time_type), INTENT(INOUT):: t1
    INTEGER, INTENT(IN):: dhours

    ! Local variables:
    INTEGER:: i, dleap

    DO i = 1, dhours
        t1%hour = t1%hour + 1
        IF (t1%hour == 24) THEN
            t1%hour = 0
            t1%day = t1%day + 1
            dleap = 0
            IF (MOD(t1%year, 4) == 0) THEN
                dleap = 1
                IF ((MOD(t1%year, 100) == 0) .AND. (MOD(t1%year, 400) /= 0)) &
                    dleap = 0
            END IF
            IF ((t1%day > 31) .AND. ((t1%month == 1) .OR. (t1%month == 3) &
                .OR. (t1%month == 5) .OR. (t1%month == 7) .OR. (t1%month == 8) &
                .OR. (t1%month == 10) .OR. (t1%month == 12))) THEN
                t1%day = 1
                t1%month = t1%month + 1
                IF (t1%month == 13) THEN
                    t1%month = 1
                    t1%year = t1%year + 1
                END IF
            END IF
            IF ((t1%day > 30) .AND. ((t1%month == 4) .OR. (t1%month == 6) &
                .OR.(t1%month == 9) .OR. (t1%month == 11))) THEN
                t1%day = 1
                t1%month = t1%month + 1
            END IF
            IF ((t1%day > (28 + dleap)) .AND. (t1%month == 2)) THEN
                t1%day = 1
                t1%month = 3
            END IF
        END IF
    END DO

END SUBROUTINE increase_time

!===============================================================================
SUBROUTINE elapsed_time(t1, t2, t_elapsed)
    ! Return elapsed time between two time objects (in hours).
    INTEGER, INTENT(OUT):: t_elapsed
    TYPE(time_type), INTENT(IN):: t1, t2
    
    ! Local variables:
    TYPE(time_type):: t0
                
    t0 = t1
    t_elapsed = 0
    DO
        t_elapsed = t_elapsed + 1
        CALL increase_time(t0, 1)
        IF (t0 == t2) EXIT
    END DO
    
END SUBROUTINE elapsed_time

!===============================================================================
SUBROUTINE get_restart_time()
    ! Find the last magnetic (bc) snapshot in data directory, and set
    ! restart_time to this.
    
    ! Local variables:
    CHARACTER(15):: snap0
    INTEGER:: io
    
    IF (control) THEN
        CALL SYSTEM('ls -d "'//TRIM(datadir)//'"/b????????.??.nc | xargs -n 1 basename > list.tmp')
        OPEN(1, file = 'list.tmp')
        DO
            READ(1,'(15a)', iostat=io) snap0
            IF (io /= 0) EXIT
        END DO
        CLOSE(1)
        CALL SYSTEM('rm list.tmp')
        restart_time(1:11) = snap0(2:12)
        PRINT*,'--- Restarting from ',restart_time,' ---'
    END IF
        
    CALL MPI_BCAST(restart_time, LEN(restart_time), MPI_CHAR, 0, comm, err)
        
END SUBROUTINE get_restart_time

!===============================================================================
SUBROUTINE initialize_time()
! Set current time correctly for start of simulation and get length of
! simulation in hours.
    
    ! Local variables:
    TYPE(time_type):: time_restart  ! Restart time of simulation in time_type &
                                    ! format.
    TYPE(time_type):: time_start  ! Start time of simulation in time_type format.

    ! Convert from strings to time objects:
    IF (restart) time_restart = str_to_time(restart_time)
    time_start = str_to_time(start_time)
    time_end = str_to_time(end_time)
    
    ! Set current time to start or restart time as appropriate:
    IF (restart) THEN
        time_now = time_restart
        IF (start_time == restart_time) THEN
            seconds = 0.0_d
        ELSE
            CALL elapsed_time(time_start, time_restart, hours)
        END IF
        seconds = REAL(hours)*3600.0_d
    ELSE
        time_now = time_start
        hours = 0
        seconds = 0.0_d
    END IF
    
    ! Get total length of simulation in hours:
    CALL elapsed_time(time_start, time_end, hours_total)
    
END SUBROUTINE initialize_time

!*******************************************************************************
END MODULE time
!*******************************************************************************
