!****************************************************************
module trace
!****************************************************************

    use shared
    use readnc
    
    implicit none

    private
    public :: compute_chmap, compute_fl, compute_windmap, compute_flhmap, compute_flh
    public :: compute_almap, compute_al, compute_hr
    
    contains   
  
    !****************************************************************
    subroutine compute_al(wkpath, infile)
        character*(*), intent(in) :: wkpath, infile
        character*(*), parameter :: ALFILE='al.unf'
        real(d), parameter :: MAXERROR=1e-1_d
        real(d), parameter :: MINB=1.0e-4_d
        integer, parameter :: NMAX=1000        
        ! ---
        ! Compute average alpha = j.B/B^2 along field lines
        ! for given field line startpoints
        ! in infile, and output the result to a new file.
        ! ---
        integer :: nfl, i, j, k
        real(d), dimension(:,:), allocatable :: x0
        real(d), dimension(:,:,:), allocatable :: xl
        integer, dimension(:), allocatable :: flag
        real(d), dimension(:), allocatable :: al, l, vl
        real(d) :: b1(3), j1(3), dd
               
        ! (1) Read starting points from file:
        ! -----------------------------------
        open(1, file=wkpath//infile, form='unformatted')
        read(1) nfl
        allocate(x0(nfl,3))
        read(1) x0
        close(1)
        
        if (VERBOSE) print*,'Computing ', nfl, ' field lines...'
        
        ! (2) Trace fieldlines:
        ! ---------------------
        allocate(xl(nfl,NMAX,3), flag(nfl))
        flag = 0
        call fieldline(x0, MAXERROR, MINB, xl, flag)
        
        ! (3) Compute averaged alpha for each curve:
        ! ------------------------------------------
        allocate(al(nfl), l(NMAX), vl(NMAX))
        al = 0.0_d
        !$omp parallel private(j,k,l,vl,b1,j1)
        !$omp do
        do i=1,nfl
            k = 0
            l = 0.0_d
            vl = 0.0_d
            do j=1,nmax
                if (xl(i,j,3).gt.-50.0_d) then
                    call interpB(xl(i,j,:), b1)
                    call interpJ(xl(i,j,:), j1)
                    vl(j) = sum(b1*j1)/sum(b1*b1)
                    if (k.gt.0) then
                        dd = sum((xl(i,j,:) - xl(i,j-1,:))**2) + 1.0d-12
                        l(j) = l(j-1) + dsqrt(dd)
                        al(i) = al(i) + 0.5_d*(l(j)-l(j-1))*(vl(j) + vl(j-1))
                    end if
                    k = k+1
                end if
            end do
            al(i) = al(i)/maxval(l)
        end do
        !$omp end do
        !$omp end parallel               
        
        ! (4) Output alpha:
        ! -------------------------------
        open(unit=1, file=wkpath//ALFILE, form='unformatted')
        write(1) al/RSUN
        close(1)
            
    end subroutine compute_al
 
    !****************************************************************
    subroutine compute_almap(wkpath, snap,  nsm, npm)
        character*(*), intent(in) :: wkpath, snap
        integer, intent(in) :: nsm, npm
        real(d), parameter :: MAXERROR=1e-1_d
        real(d), parameter :: MINB=1.0e-4_d
        integer, parameter :: NMAX=1000        
        ! ---
        ! Compute map of field-line averaged alpha=j.B/B^2 on lower boundary.
        ! ---
        integer :: i, j, k
        real(d), dimension(:,:), allocatable :: x0
        real(d), dimension(:,:,:), allocatable :: xl
        integer, dimension(:), allocatable :: flag
        real(d), dimension(:), allocatable :: al, l, vl, sg, pg
        real(d) :: b1(3), j1(3), dd, dsm, dpm
               
        ! (1) Generate starting points:
        ! -----------------------------
        if (VERBOSE) print*,'Tracing ', nsm*npm,' field lines upwards [almap]'
        allocate(x0(nsm*npm,3))
        call makeGridSphere(exp(r(1)), nsm, npm, x0)
        
        ! (2) Trace fieldlines:
        ! ---------------------
        allocate(xl(nsm*npm,NMAX,3), flag(nsm*npm))
        flag = 0
        call fieldline(x0, MAXERROR, MINB, xl, flag)
        
        ! (3) Compute field line helicity for each curve:
        ! -----------------------------------------------
        allocate(al(nsm*npm), l(NMAX), vl(NMAX))
        al = 0.0_d
        !$omp parallel private(j,k,l,vl,b1,j1)
        !$omp do
        do i=1,nsm*npm
            k = 0
            l = 0.0_d
            vl = 0.0_d
            do j=1,nmax
                if (xl(i,j,3).gt.-50.0_d) then
                    call interpB(xl(i,j,:), b1)
                    call interpJ(xl(i,j,:), j1)
                    vl(j) = sum(b1*j1)/sum(b1*b1)
                    if (k.gt.0) then
                        dd = sum((xl(i,j,:) - xl(i,j-1,:))**2) + 1.0d-12
                        l(j) = l(j-1) + dsqrt(dd)
                        al(i) = al(i) + 0.5_d*(l(j)-l(j-1))*(vl(j) + vl(j-1))
                    end if
                    k = k+1
                end if
            end do
            al(i) = al(i)/maxval(l)           
        end do
        !$omp end do
        !$omp end parallel               

        ! Make arrays of cell edges for use by meshgrid plotting:
        allocate(sg(nsm+1), pg(npm+1))
        dsm = 2.0_d/float(nsm)
        dpm = TWOPI/float(npm)
        sg = (/( dble(i-1)/dble(nsm-1)*(s(ns)-s(1)-dsm) + s(1), i=1,nsm+1 )/)   
        sg(1) = -1.0d0
        sg(nsm+1) = 1.0d0
        pg = (/( dble(i-1)/dble(npm-1)*(p(np) - p(1)-dpm), i=1,npm+1 )/)

        ! (4) Output field line helicity:
        ! -------------------------------
        open(unit=1, file=wkpath//'almap_'//snap//'.unf', form='unformatted')
        write(1) sg
        write(1) pg
        write(1) al/RSUN
        close(1)
            
    end subroutine compute_almap     
 
    !****************************************************************
    subroutine compute_chmap(wkpath, snap, nsm, npm)
        character*(*), intent(in) :: wkpath, snap
        integer, intent(in) :: nsm, npm
        real(d), parameter :: MAXERROR=1e-1_d
        real(d), parameter :: MINB=1.0e-4_d
        ! ---
        ! Compute "coronal hole" map of open field-line footpoints on lower boundary.
        ! ---
        real(d), dimension(:,:), allocatable :: x0
        real(d), dimension(:,:,:), allocatable :: xl
        integer, dimension(:), allocatable :: flag, cmplt
        integer, dimension(:,:), allocatable :: chmap
        real(d), dimension(:), allocatable :: sm, pm, sg, pg
        integer :: i, ism, ipm
        real(d) :: s0, p0, dsm, dpm
        
        ! (1) Trace field lines upward from uniform grid:
        ! -----------------------------------------------
        if (VERBOSE) print*,'Tracing ', nsm*npm,' field lines upwards [chmap]'
        allocate(x0(nsm*npm,3))
        call makeGridSphere(exp(r(1)), nsm, npm, x0)
        allocate(xl(nsm*npm,2,3), flag(nsm*npm), cmplt(nsm*npm))
        flag = 0
        call fieldline(x0, MAXERROR, MINB, xl, flag, .true.)
        allocate(chmap(nsm,npm))
        cmplt = 0.0d0
        where (flag.eq.6) cmplt = -1.0d0
        where (flag.eq.9) cmplt = 1.0d0
        chmap = reshape(int(cmplt), (/nsm, npm/), order=(/2, 1/))
        
        ! (2) Trace field lines downward from uniform grid on outer boundary:
        ! -------------------------------------------------------------------
        if (VERBOSE) print*,'Tracing ', nsm*npm,' field lines downwards [chmap]'
        call makeGridSphere(exp(r(nr)), nsm, npm, x0)
        call fieldline(x0, MAXERROR, MINB, xl, flag, .true.)
        ! - If end-points are on solar surface, change nearest pixel of chmap:
        allocate(sm(nsm), pm(npm))
        dsm = 2.0_d/float(nsm)
        dpm = TWOPI/float(npm)
        sm = (/( dble(i-1)/dble(nsm-1)*(s(ns)-s(1)-dsm) + s(1) + 0.5_d*dsm, i=1,nsm )/)   
        pm = (/( dble(i-1)/dble(npm-1)*(p(np) - p(1)-dpm) + 0.5*dpm, i=1,npm )/)
        do i=1,nsm*npm
            if (flag(i).eq.6) then
                s0 = xl(i,2,3)/sqrt(xl(i,2,1)**2 + xl(i,2,2)**2 + xl(i,2,3)**2)
                p0 = mod(atan2(xl(i,2,2),xl(i,2,1)) + TWOPI, TWOPI)
                ism = minloc(abs(s0 - sm), 1)
                ipm = minloc(abs(p0 - pm), 1)
                chmap(ism,ipm) = -1
            end if
            if (flag(i).eq.9) then
                s0 = xl(i,1,3)/sqrt(xl(i,1,1)**2 + xl(i,1,2)**2 + xl(i,1,3)**2)
                p0 = mod(atan2(xl(i,1,2),xl(i,1,1)) + TWOPI, TWOPI)
                ism = minloc(abs(s0 - sm), 1)
                ipm = minloc(abs(p0 - pm), 1)
                chmap(ism,ipm) = 1
            end if
        end do

        ! Make arrays of cell edges for use by meshgrid plotting:
        allocate(sg(nsm+1), pg(npm+1))
        sg = (/( dble(i-1)/dble(nsm-1)*(s(ns)-s(1)-dsm) + s(1), i=1,nsm+1 )/)   
        sg(1) = -1.0d0
        sg(nsm+1) = 1.0d0
        pg = (/( dble(i-1)/dble(npm-1)*(p(np) - p(1)-dpm), i=1,npm+1 )/)
        
        ! (3) Output coronal hole map to file:
        ! ------------------------------------
        open(unit=1, file=wkpath//'chmap_'//snap//'.unf', form='unformatted')
        write(1) sg
        write(1) pg
        write(1) chmap
        close(1)
            
    end subroutine compute_chmap
  
    !****************************************************************
    subroutine compute_fl(wkpath, infile)
        character*(*), intent(in) :: wkpath, infile
        real(d), parameter :: MAXERROR=1e-1_d
        real(d), parameter :: MINB=1.0e-4_d
        integer, parameter :: NMAX=5000
        character*(*), parameter :: XLFILE='xl.unf'
        ! ---
        ! Compute field lines from startpoints in infile and output full
        ! field line coordinates to file.
        ! ---
        integer :: nfl
        real(d), dimension(:,:), allocatable :: x0
        real(d), dimension(:,:,:), allocatable :: xl
        integer, dimension(:), allocatable :: flag

        ! (1) Read starting points from file:
        ! -----------------------------------
        open(1, file=wkpath//infile, form='unformatted')
        read(1) nfl
        allocate(x0(nfl,3))
        read(1) x0
        close(1)

        if (VERBOSE) print*,'Computing ', nfl, ' field lines...'

        ! (2) Trace fieldlines:
        ! ---------------------
        allocate(xl(nfl,NMAX,3), flag(nfl))
        flag = 0
        call fieldline(x0, MAXERROR, MINB, xl, flag)

        ! (3) Output fieldlines:
        ! ----------------------
        open(unit=1, file=wkpath//XLFILE, form='unformatted')
        write(1) NMAX
        write(1) xl
        close(1)
            
    end subroutine compute_fl

    !****************************************************************
    subroutine compute_flh(wkpath, infile)
        character*(*), intent(in) :: wkpath, infile
        character*(*), parameter :: FLHFILE='flh.unf'
        real(d), parameter :: MAXERROR=1e-8_d
        real(d), parameter :: MINB=1.0e-4_d
        integer, parameter :: NMAX=4000
        real(d), parameter :: MAXSTEP=0.5_d
        ! ---
        ! Compute field line helicity for given field line startpoints
        ! in infile, and output the result to a new file.
        ! ---
        integer :: nfl, i, j, k, j1
        real(d), dimension(:,:), allocatable :: x0
        real(d), dimension(:,:,:), allocatable :: xl
        integer, dimension(:), allocatable :: flag
        real(d), dimension(:), allocatable :: flh, l, vl
        real(d), dimension(:), allocatable :: rad0, rad1, brm0, brm1, br0
        real(d) :: a1(3), b1(3), dd, sm, pm
               
        ! (1) Read starting points from file:
        ! -----------------------------------
        open(1, file=wkpath//infile, form='unformatted')
        read(1) nfl
        allocate(x0(nfl,3))
        read(1) x0
        close(1)
        
        if (VERBOSE) print*,'Computing ', nfl, ' field lines...'

        ! (2) Trace fieldlines:
        ! ---------------------
        allocate(xl(nfl,NMAX,3), flag(nfl))
        flag = 0

        call fieldline(x0, MAXERROR, MINB, xl, flag, .false., MAXSTEP)

        ! (3) Compute field line helicity for each curve:
        ! -----------------------------------------------
        allocate(flh(nfl), l(NMAX), vl(NMAX))
        allocate(brm0(nfl), brm1(nfl), rad0(nfl), rad1(nfl), br0(nfl))
        flh = 0.0_d
        !$omp parallel private(j,k,l,vl,a1,b1,dd,j1,sm,pm)
        !$omp do
        do i=1,nfl
            k = 0
            l = 0.0_d
            vl = 0.0_d
            do j=1,nmax
                if (xl(i,j,3).gt.-50.0_d) then
                    if (k.eq.0) then
                        ! Save radius and Br of startpoint:
                        rad0(i) = sqrt(xl(i,j,1)**2 + xl(i,j,2)**2 + xl(i,j,3)**2)
                        sm = xl(i,j,3)/rad0(i)
                        pm = mod(atan2(xl(i,j,2),xl(i,j,1)) + TWOPI, TWOPI)
                        call interpB(xl(i,j,:), b1)
                        brm0(i) = sqrt(1.0_d - sm**2)*(cos(pm)*b1(1) + sin(pm)*b1(2)) + sm*b1(3)
                    end if
                    call interpA(xl(i,j,:), a1)
                    call interpB(xl(i,j,:), b1)
                    vl(j) = sum(b1*a1)/dsqrt(sum(b1*b1))
                    if (k.gt.0) then
                        dd = sum((xl(i,j,:) - xl(i,j-1,:))**2) + 1.0d-12
                        l(j) = l(j-1) + dsqrt(dd)
                        flh(i) = flh(i) + 0.5_d*(l(j)-l(j-1))*(vl(j) + vl(j-1))
                    end if
                    j1 = j
                    k = k+1
                end if
            end do
            ! Save radius and Br of endpoint:
            rad1(i) = sqrt(xl(i,j1,1)**2 + xl(i,j1,2)**2 + xl(i,j1,3)**2)
            sm = xl(i,j1,3)/rad1(i)
            pm = mod(atan2(xl(i,j1,2),xl(i,j1,1)) + TWOPI, TWOPI)
            call interpB(xl(i,j1,:), b1)
            brm1(i) = sqrt(1.0_d - sm**2)*(cos(pm)*b1(1) + sin(pm)*b1(2)) + sm*b1(3)
            ! Make array of br at seed points:
            sm = x0(i,3)/sqrt(x0(i,1)**2 + x0(i,2)**2 + x0(i,3)**2)
            pm = mod(atan2(x0(i,2),x0(i,1)) + TWOPI, TWOPI)
            call interpB(x0(i,:), b1)
            br0(i) = sqrt(1.0_d - sm**2)*(cos(pm)*b1(1) + sin(pm)*b1(2)) + sm*b1(3)
        end do
        !$omp end do
        !$omp end parallel

        ! (4) Output field line helicity:
        ! -------------------------------
        open(unit=1, file=wkpath//FLHFILE, form='unformatted')
        write(1) flh*RSUN*RSUN
        write(1) rad0
        write(1) rad1
        write(1) brm0
        write(1) brm1
        write(1) br0
        close(1)
            
    end subroutine compute_flh

    !****************************************************************
    subroutine compute_flhmap(wkpath, snap,  nsm, npm, ir)
        character*(*), intent(in) :: wkpath, snap
        integer, intent(in) :: nsm, npm
        integer, intent(in) :: ir
        real(d), parameter :: MAXERROR=1e-8_d
        real(d), parameter :: MINB=1.0e-4_d
        integer, parameter :: NMAX=4000
        real(d), parameter :: MAXSTEP=0.5_d
        ! ---
        ! Compute map of field line helicity on lower boundary.
        ! ---
        integer :: i, j, k
        real(d), dimension(:,:), allocatable :: x0
        real(d), dimension(:,:,:), allocatable :: xl
        integer, dimension(:), allocatable :: flag
        real(d), dimension(:), allocatable :: flh, l, vl, sg, pg
        real(d) :: a1(3), b1(3), dd, dsm, dpm, hintegral, sm, pm
        real(d), dimension(:), allocatable ::  brm

        ! (1) Generate starting points:
        ! -----------------------------
        if (VERBOSE) print*,'Tracing ', nsm*npm,' field lines [flhmap]'
        allocate(x0(nsm*npm,3))
        call makeGridSphere(exp(r(ir)), nsm, npm, x0)
                
        ! (2) Trace fieldlines:
        ! ---------------------
        allocate(xl(nsm*npm,NMAX,3), flag(nsm*npm))
        flag = 0
        call fieldline(x0, MAXERROR, MINB, xl, flag, .false., MAXSTEP)
        
        ! (3) Compute field line helicity for each curve:
        ! -----------------------------------------------

        allocate(flh(nsm*npm), l(NMAX), vl(NMAX))
        flh = 0.0_d
        !$omp parallel private(j,k,l,vl,a1,b1,dd)
        !$omp do
        do i=1,nsm*npm
            k = 0
            l = 0.0_d
            vl = 0.0_d
            do j=1,nmax
                if (xl(i,j,3).gt.-50.0_d) then
                    call interpA(xl(i,j,:), a1)
                    call interpB(xl(i,j,:), b1)
                    vl(j) = sum(b1*a1)/dsqrt(sum(b1*b1))
                    if (k.gt.0) then
                        dd = sum((xl(i,j,:) - xl(i,j-1,:))**2) + 1.0d-12
                        l(j) = l(j-1) + dsqrt(dd)
                        flh(i) = flh(i) + 0.5_d*(l(j)-l(j-1))*(vl(j) + vl(j-1))
                    end if
                    k = k+1
                end if
            end do
        end do
        !$omp end do
        !$omp end parallel               

        ! (4) Make array of Br at same points as FLH
        ! ----------------------------------------
        allocate(brm(nsm*npm))
        do i=1,nsm*npm
            sm = x0(i,3)/sqrt(x0(i,1)**2 + x0(i,2)**2 + x0(i,3)**2)
            pm = mod(atan2(x0(i,2),x0(i,1)) + TWOPI, TWOPI)
            call interpB(x0(i,:), b1)
            brm(i) = sqrt(1.0_d - sm**2)*(cos(pm)*b1(1) + sin(pm)*b1(2)) + sm*b1(3)
        end do

        ! (5) Make arrays of cell edges for use by meshgrid plotting:
        ! ----------------------------------------
        allocate(sg(nsm+1), pg(npm+1))
        dsm = 2.0_d/float(nsm)
        dpm = TWOPI/float(npm)
        sg = (/( dble(i-1)/dble(nsm-1)*(s(ns)-s(1)-dsm) + s(1), i=1,nsm+1 )/)   
        sg(1) = -1.0d0
        sg(nsm+1) = 1.0d0
        pg = (/( dble(i-1)/dble(npm-1)*(p(np) - p(1)-dpm), i=1,npm+1 )/)

        ! (6) Compute helicity by volume integral:
        ! ----------------------------------------
        hintegral = trap3d(arg*brg + asg*bsg + apg*bpg)

        ! (4) Output field line helicity:
        ! -------------------------------
        open(unit=1, file=wkpath//'flhmap_'//snap//'.unf', form='unformatted')
        write(1) sg
        write(1) pg
        write(1) flh*RSUN*RSUN
        write(1) brm
        write(1) hintegral*RSUN**4
        write(1) exp(r(ir))
        close(1)
            
    end subroutine compute_flhmap


    !****************************************************************
    subroutine compute_hr(wkpath, snap)
        character*(*), intent(in) :: wkpath, snap
        ! ---
        ! Compute volume integrated helicity.
        ! ---
        real(d) :: hintegral

        ! Compute helicity by volume integral:
        ! ----------------------------------------
        hintegral = trap3d(arg*brg + asg*bsg + apg*bpg)

        ! Output helicity (single number) to file:
        ! ----------------------------------------
        open(unit=1, file=wkpath//'hr_'//snap//'.unf', form='unformatted')
        write(1) hintegral*RSUN**4
        close(1)
            
    end subroutine compute_hr

    
    !****************************************************************
    subroutine compute_windmap(wkpath, snap, nsm, npm)
        character*(*), intent(in) :: wkpath, snap
        integer, intent(in) :: nsm, npm
        real(d), parameter :: MAXERROR=1e-1_d
        real(d), parameter :: MINB=1.0e-4_d
        ! ---
        ! Compute solar wind connectivity map down from Schatten outer boundary 
        ! continue through original model to photosphere.
        ! Output the mapping and br, bth, bph values (bth, bph only at outer boundary).
        ! ---
        real(d), dimension(:,:), allocatable :: x0
        real(d), dimension(:,:,:), allocatable :: xl
        integer, dimension(:), allocatable :: flag
        real(d), dimension(:), allocatable :: sm, pm, rm, brm, brm1
        real(d) :: b1(3), bsign
        integer :: i
                          
        ! (1) Define uniform grid on outer (Schatten) boundary and save to file:
        ! ----------------------------------------------------------------------
        call readB(trim(wkpath)//'schat_'//trim(snap))   
        allocate(x0(nsm*npm,3), sm(nsm*npm), pm(nsm*npm), rm(nsm*npm))
        allocate(brm(nsm*npm), brm1(nsm*npm))
        call makeGridSphere(exp(r(nr)), nsm, npm, x0)
        sm = x0(:,3)/sqrt(x0(:,1)**2 + x0(:,2)**2 + x0(:,3)**2)
        pm = mod(atan2(x0(:,2),x0(:,1)) + TWOPI, TWOPI)
        ! - get br on outer boundary (outward):
        do i=1,nsm*npm
            call interpB(x0(i,:), b1)
            brm(i) = sqrt(1.0_d - sm(i)**2)*(cos(pm(i))*b1(1) + sin(pm(i))*b1(2)) + sm(i)*b1(3)
        end do
        
        open(unit=2, file=wkpath//'windmap_'//snap//'.unf', form='unformatted')
        write(2) nsm
        write(2) npm   
        write(2) sm
        write(2) pm
                                 
        ! (2) Trace downward through Schatten from uniform grid on outer boundary:
        ! ------------------------------------------------------------------------
        if (VERBOSE) print*,'Tracing ', nsm*npm,' field lines downwards through Schatten [wind]'
        allocate(xl(nsm*npm,2,3), flag(nsm*npm))
        flag = 0
        call fieldline(x0, MAXERROR, MINB, xl, flag, .true.)
        ! Correct radius to exact boundary:
        rm = sqrt(xl(:,1,1)**2 + xl(:,1,2)**2 + xl(:,1,3)**2)
        x0(:,1) = xl(:,1,1)/rm*exp(r(1))
        x0(:,2) = xl(:,1,2)/rm*exp(r(1))
        x0(:,3) = xl(:,1,3)/rm*exp(r(1))

        ! (3) Read in inner model and use to restore sign of B on outer boundary:
        ! ------------------------------------------------------------------------
        call readB(trim(wkpath)//trim(snap))  
        do i=1,nsm*npm
            call interpB(x0(i,:), b1)
            brm1(i) = sqrt(1.0_d - sm(i)**2)*(cos(pm(i))*b1(1) + sin(pm(i))*b1(2)) + sm(i)*b1(3)
            bsign = brm1(i)/abs(brm1(i))
            brm(i) = brm(i)*bsign
        end do        
        write(2) brm
        write(2) brm1

        ! (4) Continue downward through inner model:
        ! ------------------------------------------
        if (VERBOSE) print*,'Tracing ', nsm*npm,' field lines downwards [wind]'
        flag = 0
        call fieldline(x0, MAXERROR, MINB, xl, flag, .true.) 
        do i=1,nsm*npm
            select case (flag(i))
                case (6) 
                    sm(i) = xl(i,2,3)/sqrt(xl(i,2,1)**2 + xl(i,2,2)**2 + xl(i,2,3)**2)
                    pm(i) = mod(atan2(xl(i,2,2),xl(i,2,1)) + TWOPI, TWOPI)              
                    call interpB(xl(i,2,:), b1)
                    brm(i) = sqrt(1.0_d - sm(i)**2)*(cos(pm(i))*b1(1) + sin(pm(i))*b1(2)) + sm(i)*b1(3)
                case (9) 
                    sm(i) = xl(i,1,3)/sqrt(xl(i,1,1)**2 + xl(i,1,2)**2 + xl(i,1,3)**2)
                    pm(i) = mod(atan2(xl(i,1,2),xl(i,1,1)) + TWOPI, TWOPI) 
                    call interpB(xl(i,1,:), b1)
                    brm(i) = sqrt(1.0_d - sm(i)**2)*(cos(pm(i))*b1(1) + sin(pm(i))*b1(2)) + sm(i)*b1(3)                    
                case default     ! mark field lines that don't reach photosphere
                    sm(i) = -999.0_d
                    pm(i) = -999.0_d
                    brm(i) = -999.0_d
            end select
        end do
        write(2) sm
        write(2) pm
        write(2) brm
        close(2)
             
    end subroutine compute_windmap
    
!-----------------------------------------------------------------------
  
    !****************************************************************
    subroutine fieldline(x0, maxError, minB, xl, endflag, mapOnly, maxStep)
        real(d), intent(in), dimension(:,:) :: x0
        real(d), intent(in) :: maxError
        real(d), intent(in) :: minB
        logical, intent(in), optional :: mapOnly
        real(d), intent(in), optional :: maxStep
        real(d), intent(inout), dimension(:,:,:) :: xl
        integer, intent(inout), dimension(:), optional :: endflag
        ! ---
        ! Trace a bunch of fieldlines in both directions from given
        ! startpoints x0.
        ! Output full field line arrays xl, unless mapOnly=.true. in which
        ! case only endpoints are output.
        !
        ! Codes for endflag is as follows:
        !  0 -- neither end reached boundary
        !  +1 -- start closed
        !  +2 -- end closed
        !  +4 -- start open
        !  +8 -- end open
        !
        ! Optional argument maxStep is maximum stepsize to take, as a fraction of the
        ! radial cell height at the top of the domain. Default is 0.5 but should reduce it
        ! when computing field-line integrals.
        ! ---
        integer :: nfl0, nmax
        logical :: mapOnlyFlag, fullFlag
        real(d) :: rMin, rMax
        real(d) :: ddirn, maxdl, dl, r1, r2, dl_dt, error, k1r, rl
        real(d), dimension(3) :: x1, x2, dx1, dx2, k1, k2
        integer :: nxt, cntr, dirn, i

        if (present(mapOnly)) then
            mapOnlyFlag = mapOnly
        else
            mapOnlyFlag = .false.
        end if
        fullFlag = (.not.mapOnlyFlag)

        nfl0 = size(xl,1)

        if (mapOnlyFlag) then
            nmax = 20003
        else
            nmax = size(xl,2)
        end if

        ! Extents of domain:
        rMin = minval(r)
        rMax = maxval(r)

        ! Convert from rho to r:
        rMin = dexp(rMin)
        rMax = dexp(rMax)

        ! Maximum allowed step-size:
        if (present(maxStep)) then
            maxdl = maxStep*(dexp(r(nr))-dexp(r(nr-1)))
        else
            maxdl = 0.5_d*(dexp(r(nr))-dexp(r(nr-1)))
        end if

        ! Initialize all field lines to "incomplete" status:
        endflag = 0
        
        ! Main loop:
        
        !$omp parallel private(dirn,ddirn,dl,nxt,cntr,k1,k2,dl_dt,x1,x2,r2,dx1,dx2,error)
        !$omp do
        do i=1,nfl0
            ! Initialise variables:
            xl(i,:,:) = 0.0_d
            xl(i,:,3) = -99.0_d

            do dirn=-1,1,2
                ! Reverse direction of field for backward and forward tracing:
                ddirn = dble(dirn)
                dl=maxdl

                if (mapOnlyFlag) then
                    nxt = (dirn+3)/2
                else
                    nxt = nmax/2
                end if

                cntr = 0

                ! Add startpoint to output array:
                xl(i,nxt,:) = x0(i,:)

                ! Interpolate k1:
                call interpB(xl(i,nxt,:), k1)
                dl_dt = dsqrt(sum(k1*k1))
                if (dl_dt < minB) exit   ! stop if null is reached
                dl_dt = dl_dt*ddirn
                k1 = k1/dl_dt
                do
                    ! Compute midpoint:
                    x2 = xl(i,nxt,:) + dl*k1
                    r2 = dsqrt(sum(x2*x2))
                    ! If outside boundary, do Euler step to boundary and stop:
                    if ((r2.lt.rMin).or.(r2.gt.rMax)) then
                        rl = dsqrt(sum(xl(i,nxt,:)*xl(i,nxt,:)))
                        k1r = (xl(i,nxt,1)*k1(1) + xl(i,nxt,2)*k1(2) &
                            + xl(i,nxt,3)*k1(3))/rl
                        if (r2.lt.rMin) then
                            dl = (rMin - rl)/k1r
                        else
                            dl = (rMax - rl)/k1r
                        end if
                        if (fullFlag) nxt = nxt + dirn
                        cntr = cntr + 1
                        if (cntr .ge. nmax/2) then
                            print*,'WARNING: Field line length exceeds NMAX'
                            exit
                        end if
                        if (mapOnlyFlag) then
                            xl(i,nxt,:) = xl(i,nxt,:) + dl*k1
                        else
                            xl(i,nxt,:) = xl(i,nxt-dirn,:) + dl*k1
                        end if
                        if (r2.lt.rMin) then
                            if (dirn.eq.-1) then
                                endflag(i)=endflag(i)+1
                            else
                                endflag(i)=endflag(i)+2
                            end if
                        else
                            if (dirn.eq.-1) then
                                endflag(i)=endflag(i)+4
                            else
                                endflag(i)=endflag(i)+8
                            end if
                        end if
                        exit
                    end if
                    ! Interpolate k2:
                    call interpB(x2, k2)
                    dl_dt = dsqrt(sum(k2*k2))
                    dl_dt = dl_dt*ddirn
                
                    ! Compute first and second-order update:
                    k2 = 0.5_d*(k1 + k2/dl_dt)

                    if (sum(k2*k2) < minB*minB) exit   ! stop if null is reached
                    dx1 = dl*k1
                    dx2 = dl*k2

                    ! Estimate error from difference:
                    error = sum((dx1-dx2)*(dx1-dx2))

                    ! Modify step size depending on error:
                    if (error.lt.maxerror) then
                        dl = maxdl
                    else
                        dl = min(maxdl, 0.85_d*dabs(dl)*(maxerror/error)**0.25_d)
                    end if

                    ! Update if error is small enough:
                    if (error.le.maxerror) then
                        x1 = xl(i,nxt,:) + dx2
                        r1 = dsqrt(sum(x1*x1))
                        ! Return midpoint if full step leaves domain:
                        if ((r1.lt.rMin).or.(r1.gt.rMax)) x1 = x2
                        if (fullFlag) nxt = nxt + dirn
                        cntr = cntr + 1
                        if (cntr .ge. nmax/2) then
                            print*,'WARNING: Field line length exceeds NMAX'
                            exit
                        end if
                        xl(i,nxt,:) = x1
                        ! Interpolate k1 at next point:
                        call interpB(x1, k1)
                        dl_dt = dsqrt(sum(k1**2))
                        if (dl_dt < minB) exit   ! stop if null is reached
                        dl_dt = dl_dt*ddirn
                        k1 = k1/dl_dt
                    end if
            end do
        end do
        end do
        !$omp end do
        !$omp end parallel

    end subroutine fieldline

end module trace
