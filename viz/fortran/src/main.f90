!****************************************************************
program main
!----------------------------------------------------------------
! Main program: field line tracing and helicity computation
!               for a dumfric snapshot
!----------------------------------------------------------------

use shared
use readnc
use trace
use xtrapol
use outeqm
use fluxes
use omp_lib

implicit none

character(200) :: wkpath, snap, mode, snap_nxt, dhours_str, ir_str, r_hb_str, v1_str
integer :: omp_threads, omp_thread1, nca, dhours, ir
real(d) :: r_hb, v1

!$omp parallel
omp_threads = omp_get_num_threads()
omp_thread1 = omp_get_thread_num()
if (VERBOSE .and. (omp_thread1.eq.0)) print*,'Number of OpenMp threads:',omp_threads
!$omp end parallel

! CHECK COMMAND LINE ARGUMENTS:
nca = command_argument_count()
if (nca.lt.3) then
    print*,'Error! Need at least 2 command line arguments'
end if

! READ COMMAND LINE ARGUMENTS:
call get_command_argument(1, wkpath)    ! working directory
call get_command_argument(2, snap)      ! simulation snapshot
call get_command_argument(3, mode)      ! what to calculate

select case (trim(mode))
    case ('al')     ! compute field-line averaged alpha (FF param) for list of field lines
        call readB(trim(wkpath)//trim(snap))
        call readJ(trim(wkpath)//trim(snap))
        call compute_al(trim(wkpath), 'x0.unf')
    case ('almap')     ! compute map of field-line averaged alpha on lower boundary
        call readB(trim(wkpath)//trim(snap))
        call readJ(trim(wkpath)//trim(snap))   
        call compute_almap(trim(wkpath), trim(snap), ns, np)             
    case ('chmap')  ! "coronal hole" map (open/closed footpoints)
        call readB(trim(wkpath)//trim(snap))
        call compute_chmap(trim(wkpath), trim(snap), ns, np)
    case ('fl')     ! trace list of field lines
        call readB(trim(wkpath)//trim(snap))
        call compute_fl(trim(wkpath), 'x0.unf')
    case ('flh')    ! compute FL helicity for list of field lines
        call readB(trim(wkpath)//trim(snap), .true.)
        call compute_flh(trim(wkpath), 'x0.unf')
    case ('flhmap')    ! compute map of FL helicity on lower boundary
        if (nca.gt.3) then
            call get_command_argument(4, ir_str)
            read(ir_str,*) ir
        else
            ir = 1
        end if
        call readB(trim(wkpath)//trim(snap), .true., .true.)
        call compute_flhmap(trim(wkpath), trim(snap), 2*ns, 2*np, ir)
    case ('hr')    ! compute relative helicity
        call readB(trim(wkpath)//trim(snap), .true., .true.)
        call compute_hr(trim(wkpath), trim(snap))
    case ('pfrf')   ! compute reference potential field
        call readB(trim(wkpath)//trim(snap))
        call compute_pfrf_inner(trim(wkpath), trim(snap))           
    case ('pfss')   ! compute PFSS
        call readB(trim(wkpath)//trim(snap))
        call compute_pfss_inner(trim(wkpath), trim(snap))
    case ('outeqm')   ! compute outflow equilibrium
    ! [requires v1 parameter in km/s as 4th command line argument]
        call get_command_argument(4, v1_str)
        read(v1_str,*) v1
        call readB(trim(wkpath)//trim(snap))
        call compute_outeqm_inner(trim(wkpath), trim(snap), v1)
    case ('wind')   ! add Schatten (PFSS) extension and compute connectivity map
        call get_command_argument(4, r_hb_str)
        read(r_hb_str,*) r_hb
        call readB(trim(wkpath)//trim(snap))
        call compute_schatten(trim(wkpath), trim(snap), r_hb)
        call compute_windmap(trim(wkpath), trim(snap), ns, np)
    case ('hfluxes') ! compute total helicity, dissipation and boundary fluxes
    ! [note: this requires electric field in the snapshot files;
    !        also requires a 4th command-line argument giving next snapshot file
    !        and a 5th argument giving the time difference in hours]
        call get_command_argument(4, snap_nxt)
        call get_command_argument(5, dhours_str)
        read(dhours_str,*) dhours
        call compute_hfluxes(wkpath, snap, snap_nxt, dhours)
    case ('pol') ! compute poloidal potential
        call readB(trim(wkpath)//trim(snap), .false., .false., .true.)
        call writeP(trim(wkpath)//'poloidal.unf')
end select

end program main
