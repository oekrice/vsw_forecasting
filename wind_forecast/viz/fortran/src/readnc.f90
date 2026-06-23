!****************************************************************
module readnc
!****************************************************************

    use shared
    use netcdf
    use ffts

!----------------------------------------------------------------
! Read netcdf files and set up interpolators.
!----------------------------------------------------------------
    implicit none
    
    private
    public :: readB, readJ, interpA, interpB, interpJ, computeBfromA
    public :: readvec_gpts, writeP
    public :: br0, br1
    public :: arg, asg, apg, brg, bsg, bpg

    logical :: store_glb=.false.
    real(d), dimension(:,:,:,:,:), allocatable:: db, da, dj ! interpolation data
    real(d), dimension(:,:), allocatable:: br0, br1         ! Br on lower boundary
    real(d), dimension(:,:,:), allocatable:: arg, asg, apg, brg, bsg, bpg
    real(d), dimension(:,:,:), allocatable:: pol  ! poloidal potential

    real(d), dimension(64,64) :: m                          ! interpolation coefficients

    contains
    
    !****************************************************************
    subroutine readB(filename, vectorpot, storeflag, polflag)
        character*(*), intent(in) :: filename
        logical, intent(in), optional :: vectorpot, storeflag, polflag
        character*(*), parameter :: GAUGE='minimal'
        ! ---
        ! Read magnetic field from a DuMFric netcdf file, and create interpolator.
        ! ---
        real(d), dimension(:,:,:), allocatable :: br,bs,bp

        if (allocated(db)) then
            deallocate(db)
            deallocate(br0)
            deallocate(br1)
        end if
        if (allocated(da)) deallocate(da)

        if (present(storeflag)) store_glb = storeflag

        if (VERBOSE) print*,'Reading '//filename
        
        ! Read grid (store in shared variables):
        call readDims(filename)
        call readCoords(filename)

        ! Read magnetic field and create interpolator (private to this module):
        allocate(br(nr,ns+1,np+1),bs(nr+1,ns,np+1),bp(nr+1,ns+1,np))
        call readvec_faces(filename,br,bs,bp)

        allocate(db(8,nr,ns,np,3))
        call prepareInterpB(br, bs, bp)   

        ! Compute vector potential:   
        if (present(vectorpot)) then
            if ((vectorpot).and.VERBOSE) &
                print*,'Vector potential in ', GAUGE, ' gauge'
            if (vectorpot) call computeA(br, bs, bp, GAUGE)
        end if

        ! Compute poloidal potential:
        if (present(polflag).and.VERBOSE) print*,'Poloidal potential'
        if (present(polflag)) call computeP(br)

        ! Store Br arrays at cell faces on inner and outer boundaries:
        allocate(br0(ns+1,np+1), br1(ns+1,np+1))
        br0 = br(1,:,:)
        br1 = br(nr,:,:)
        
        deallocate(br, bs, bp)
    
    end subroutine readB
    
    !****************************************************************
    subroutine readJ(filename)
        character*(*), intent(in) :: filename
        ! ---
        ! Read current density at grid points from a DuMFric netcdf file, and create interpolator.
        ! ---
        real(d), dimension(:,:,:), allocatable :: jr, js, jp

        if (allocated(dj)) deallocate(dj)

        if (VERBOSE) print*,'Reading '//filename
        
        ! Read magnetic field and create interpolator (private to this module):
        allocate(jr(nr,ns,np),js(nr,ns,np),jp(nr,ns,np))
        call readvec_gpts(filename, 'j', jr, js, jp)
        allocate(dj(8,nr,ns,np,3))
        call prepareInterpJ(jr, js, jp)   
        
        deallocate(jr, js, jp)
    
    end subroutine readJ    

    !****************************************************************
    subroutine interpA(x1, a1)
        real(d), dimension(3), intent(in) :: x1
        real(d), dimension(3), intent(out) :: a1
        ! ---
        ! C1 tricubic interpolation for vector on regular grid
        ! - using method of Lekien & Marsden 2005.
        ! ---
        real(d) :: r1, s1, p1
        real(d) :: fr, fs, fp
        integer :: ir, is, ip

        real(d), dimension(3) :: al
        real(d) :: fact
        integer :: i,j

        ! Identify required cell.
        ! -----------------------
        r1 = dsqrt(sum(x1*x1))
        s1 = x1(3)/r1
        p1 = mod(atan2(x1(2),x1(1)) + TWOPI, TWOPI)
        r1 = log(r1)

        ir = floor((r1-r(1))/dr) + 1
        is = floor(abs(s1 + 1.0_d)/ds) + 1
        ip = floor(p1/dp) + 1

        fr = (r1 - r(ir))/dr
        fs = (s1 - s(is))/ds
        fp = (p1 - p(ip))/dp

        if (ir.eq.nr) then
            ir = ir-1
            fr = 1.0_d
        end if
        if (is.eq.ns) then
            is = is-1
            fs = 1.0_d
        end if
        if (ir.eq.0) then
            ir = 1
            fr = 0.0_d
        end if
        
        ! Compute interpolation coefficients.
        ! -----------------------------------
        ! Loop through each combination of powers
        a1 = 0.0_d
        do i = 1,64
            al = 0.0_d
            do j = 0,7
                al = al + m(i,8*j+1)*da(j+1,ir,is,ip,:) &
                    + m(i,8*j+2)*da(j+1,ir+1,is,ip,:) &
                    + m(i,8*j+3)*da(j+1,ir,is+1,ip,:) &
                    + m(i,8*j+4)*da(j+1,ir+1,is+1,ip,:) &
                    + m(i,8*j+5)*da(j+1,ir,is,ip+1,:) &
                    + m(i,8*j+6)*da(j+1,ir+1,is,ip+1,:) &
                    + m(i,8*j+7)*da(j+1,ir,is+1,ip+1,:) &
                    + m(i,8*j+8)*da(j+1,ir+1,is+1,ip+1,:)
            end do
            fact = fr**mod(i-1,4)*fs**mod((i-1)/4,4)*fp**((i-1)/16)
            a1 = a1 + al*fact
        end do

    end subroutine interpA   
    
    !****************************************************************
    subroutine interpB(x1, b1)
        real(d), dimension(3), intent(in) :: x1
        real(d), dimension(3), intent(out) :: b1
        ! ---
        ! C1 tricubic interpolation for vector on regular grid
        ! - using method of Lekien & Marsden 2005.
        ! ---
        real(d) :: r1, s1, p1
        real(d) :: fr, fs, fp
        integer :: ir, is, ip

        real(d), dimension(3) :: al
        real(d) :: fact
        integer :: i,j

        ! Identify required cell.
        ! -----------------------
        r1 = dsqrt(sum(x1*x1))
        s1 = x1(3)/r1
        p1 = mod(atan2(x1(2),x1(1)) + TWOPI, TWOPI)
        r1 = log(r1)
        if (r1 < r(1)) r1=r(1)

        ir = floor((r1-r(1))/dr) + 1
        is = floor(abs(s1 + 1.0_d)/ds) + 1
        ip = floor(p1/dp) + 1

        fr = (r1 - r(ir))/dr
        fs = (s1 - s(is))/ds
        fp = (p1 - p(ip))/dp

        if (ir.ge.nr) then
            ir = ir-1
            fr = 1.0_d
        end if
        if (is.eq.ns) then
            is = is-1
            fs = 1.0_d
        end if
        if (ir.eq.0) then
            ir = 1
            fr = 0.0_d
        end if
        
        ! Compute interpolation coefficients.
        ! -----------------------------------
        ! Loop through each combination of powers
        b1 = 0.0_d
        do i = 1,64
            al = 0.0_d
            do j = 0,7
                al = al + m(i,8*j+1)*db(j+1,ir,is,ip,:) &
                    + m(i,8*j+2)*db(j+1,ir+1,is,ip,:) &
                    + m(i,8*j+3)*db(j+1,ir,is+1,ip,:) &
                    + m(i,8*j+4)*db(j+1,ir+1,is+1,ip,:) &
                    + m(i,8*j+5)*db(j+1,ir,is,ip+1,:) &
                    + m(i,8*j+6)*db(j+1,ir+1,is,ip+1,:) &
                    + m(i,8*j+7)*db(j+1,ir,is+1,ip+1,:) &
                    + m(i,8*j+8)*db(j+1,ir+1,is+1,ip+1,:)
            end do
            fact = fr**mod(i-1,4)*fs**mod((i-1)/4,4)*fp**((i-1)/16)
            b1 = b1 + al*fact
        end do

    end subroutine interpB

     !****************************************************************
    subroutine interpJ(x1, j1)
        real(d), dimension(3), intent(in) :: x1
        real(d), dimension(3), intent(out) :: j1
        ! ---
        ! C1 tricubic interpolation for vector on regular grid
        ! - using method of Lekien & Marsden 2005.
        ! ---
        real(d) :: r1, s1, p1
        real(d) :: fr, fs, fp
        integer :: ir, is, ip

        real(d), dimension(3) :: al
        real(d) :: fact
        integer :: i,j

        ! Identify required cell.
        ! -----------------------
        r1 = dsqrt(sum(x1*x1))
        s1 = x1(3)/r1
        p1 = mod(atan2(x1(2),x1(1)) + TWOPI, TWOPI)
        r1 = log(r1)

        ir = floor((r1-r(1))/dr) + 1
        is = floor(abs(s1 + 1.0_d)/ds) + 1
        ip = floor(p1/dp) + 1

        fr = (r1 - r(ir))/dr
        fs = (s1 - s(is))/ds
        fp = (p1 - p(ip))/dp

        if (ir.eq.nr) then
            ir = ir-1
            fr = 1.0_d
        end if
        if (is.eq.ns) then
        is = is-1
            fs = 1.0_d
        end if

        ! Compute interpolation coefficients.
        ! -----------------------------------
        ! Loop through each combination of powers
        j1 = 0.0_d
        do i = 1,64
            al = 0.0_d
            do j = 0,7
                al = al + m(i,8*j+1)*dj(j+1,ir,is,ip,:) &
                    + m(i,8*j+2)*dj(j+1,ir+1,is,ip,:) &
                    + m(i,8*j+3)*dj(j+1,ir,is+1,ip,:) &
                    + m(i,8*j+4)*dj(j+1,ir+1,is+1,ip,:) &
                    + m(i,8*j+5)*dj(j+1,ir,is,ip+1,:) &
                    + m(i,8*j+6)*dj(j+1,ir+1,is,ip+1,:) &
                    + m(i,8*j+7)*dj(j+1,ir,is+1,ip+1,:) &
                    + m(i,8*j+8)*dj(j+1,ir+1,is+1,ip+1,:)
            end do
            fact = fr**mod(i-1,4)*fs**mod((i-1)/4,4)*fp**((i-1)/16)
            j1 = j1 + al*fact
        end do

    end subroutine interpJ
    
    !****************************************************************
    subroutine computeBfromA(ar, as, ap, br, bs, bp)
      real(d), dimension(:,:,:), intent(inout) :: ar, as, ap
      real(d), dimension(:,:,:), intent(inout) :: br, bs, bp
      ! ---
      ! Compute B on cell faces from A on cell edges (mirrors calculation
      ! in DuMFriC code - use for testing computation of A).
      ! ---
      integer :: i, j
      integer, dimension(size(p)+1) :: opp
      real(d) :: e1, e2, fs0, fs1, fr0, ff
      real(d), dimension(size(s)) :: asins, sins
      real(d), dimension(size(r)) :: e2r

      br = 0.0_d
      bs = 0.0_d
      bp = 0.0_d

      !
      ! Prepare
      e1 = exp(dr)
      sins = sqrt(1.0_d - s**2)
      asins = asin(s)
      e2r = exp(2.0_d*r)
      !
      ! Weight ar, as, ap by edge lengths
      if (size(ar).gt.1) then
        do i=1,nr-1
            ar(i,:,:) = ar(i,:,:)*(exp(r(i+1)) - exp(r(i)))
        end do
      end if
      do i=1,nr
         e1 = exp(r(i))
         do j=1,ns-1
            as(i,j,:) = as(i,j,:)*e1*(asins(j+1)-asins(j))
         end do
         do j=1,ns
            ap(i,j,:) = ap(i,j,:)*dp*e1*sins(j)
         end do
      end do
      !
      ! Compute br in interior:
      br(:,2:ns,2:np) = ap(:,1:ns-1,:) - ap(:,2:ns,:) &
           + as(:,:,2:np) - as(:,:,1:np-1)
      br(:,2:ns,2:np) = br(:,2:ns,2:np)/ds/dp
      do i=1,nr
         br(i,2:ns,2:np) = br(i,2:ns,2:np)/e2r(i)
      end do
      !
      ! Compute bs in interior:
      bs(2:nr,:,2:np) = ap(2:nr,:,:) - ap(1:nr-1,:,:)
      if (size(ar).gt.1) bs(2:nr,:,2:np) = bs(2:nr,:,2:np) + ar(:,:,1:np-1) - ar(:,:,2:np) 
      e2 = 0.5_d*dp*(exp(2.0_d*dr) - 1.0_d)
      do i=2,nr
         do j=1,ns
            bs(i,j,:) = bs(i,j,:)/e2/e2r(i-1)/sins(j)
         end do
      end do
      !
      ! Compute bp in interior:
      bp(2:nr,2:ns,:) = as(1:nr-1,:,:) - as(2:nr,:,:)
      if (size(ar).gt.1) bp(2:nr,2:ns,:) =  bp(2:nr,2:ns,:) + ar(:,2:ns,:) - ar(:,1:ns-1,:)      
      e2 = 0.5_d*(exp(2.0_d*dr) - 1.0_d)
      do i=2,nr
         do j=2,ns
            bp(i,j,:) = bp(i,j,:)/e2/e2r(i-1)/(asins(j)-asins(j-1))
         end do
      end do
      !
      ! Outer boundary (open):
      bs(nr+1,:,2:np) = 2.0_d*bs(nr,:,2:np) - bs(nr-1,:,2:np)
      bp(nr+1,2:ns,:) = 2.0_d*bp(nr,2:ns,:) - bp(nr-1,2:ns,:)
      !
      ! Periodic boundaries in p:
      br(:,2:ns,1) = br(:,2:ns,np)
      br(:,2:ns,np+1) = br(:,2:ns,2)
      bs(2:nr+1,:,1) = bs(2:nr+1,:,np)
      bs(2:nr+1,:,np+1) = bs(2:nr+1,:,2)
      !
      ! Boundary conditions at lower boundary (js = jp = 0):
      fs0 = exp(-r(1))*exp(0.5_d*dr)
      fs1 = exp(r(1))*exp(0.5_d*dr)
      fr0 = exp(r(1))*(exp(0.5_d*dr) - exp(-0.5_d*dr))
      do j=1,ns
         ff = asin(s(j)+0.5_d*ds) - asin(s(j)-0.5_d*ds)
         bs(1,j,:) = fs0/ff*(bs(2,j,:)*fs1*ff &
              + (br(1,j,:) - br(1,j+1,:))*fr0)
      end do
      fs0 = fs0/dp
      fs1 = fs1*dp
      fr0 = exp(r(1))*(exp(0.5_d*dr) - exp(-0.5_d*dr))
      do j=2,ns
         ff = sqrt(1.0_d - (s(j)-0.5_d*ds)**2)
         do i=1,np
            bp(1,j,i) = fs0/ff*(bp(2,j,i)*fs1*ff &
                 + (br(1,j,i) - br(1,j,i+1))*fr0)
         end do
      end do
      !
      ! Boundary conditions at poles:
      opp(1:(np+1)/2) = (/ ( (np-1)/2 + i, i=1,(np+1)/2 ) /)
      opp((np+1)/2+1:np+1) = (/ ( i+1, i=1,(np+1)/2 ) /)
      br(:,1,:) = br(:,2,opp)
      br(:,ns+1,:) = br(:,ns,opp)
      bs(:,1,:) = 0.5*(bs(:,2,:) - bs(:,2,opp))
      bs(:,ns,:) = 0.5*(bs(:,ns-1,:) - bs(:,ns-1,opp))
      opp(1:(np+1)/2) = (/ ( (np-1)/2 + i, i=1,(np+1)/2 ) /)
      opp((np+1)/2+1:np) = (/ ( 1 + i, i=1,(np-1)/2 ) /)
      bp(:,1,:) = -bp(:,2,opp(1:np))
      bp(:,ns+1,:) = -bp(:,ns,opp(1:np))
      !
      ! Remove weighting again from ar, as, ap
      if (size(ar).gt.1) then
        do i=1,nr-1
            ar(i,:,:) = ar(i,:,:)/(exp(r(i+1)) - exp(r(i)))
        end do
      end if
      do i=1,nr
         e1 = exp(r(i))
         do j=1,ns-1
            as(i,j,:) = as(i,j,:)/e1/(asins(j+1)-asins(j))
         end do
         do j=2,ns-1
            ap(i,j,:) = ap(i,j,:)/dp/e1/sqrt(1-s(j)**2)
         end do
      end do
      ap(:,1,:) = 0.0_d
      ap(:,ns,:) = 0.0_d

    end subroutine computeBfromA    
    
    !----------------------------------------------------------------------

    !****************************************************************
    subroutine computeA(br, bs, bp, gauge)
        real(d), intent(in), dimension(:,:,:) :: br, bs, bp
        logical, parameter :: CHECK=.true.
        character*(*), intent(in) :: gauge
        ! ---
        ! Compute vector potential.
        ! If gauge='devore' then use DeVore-Coulomb gauge vector potential (as used by Yeates & Hornig, A&A [2016] or Lowder & Yeates, ApJ [2017].
        ! If gauge='minimal' then change to "minimal" gauge.
        ! ---
        real(d) :: ar(nr-1,ns,np), as(nr,ns-1,np), ap(nr,ns,np-1)
        real(d), dimension(:,:,:), allocatable :: br2, bs2, bp2, g
        real(d), dimension(:,:), allocatable :: aps1, Q
        real(d), dimension(:), allocatable :: lam, ev, Lpc, Vc, Lsc, Ls, Ug
        complex(d), dimension(:,:), allocatable :: rhs
        complex(d), dimension(:), allocatable :: clm
        real(d) :: f
        integer :: i, j, m
        real(d) :: Vg(ns), sc(ns-1), Uc(ns-1)
        real(d) :: er(nr), dasins(ns-1), dasinsc(ns-2), sqsc(ns-1), sqs(ns)
        complex(d) :: br00(ns-1,np-1)
        real(d) :: psi(ns-1,np-1), mu(np-1)
    
        ar = 0.0_d
        as = 0.0_d
        ap = 0.0_d
        
        ! (1) Compute Devore-Coulomb (upward) vector potential
        ! ====================================================
        ! (a) Compute contribution from radial integrals
        ! ----------------------------------------------
        ! as = -1/Ls*Sum(Sp*bp)   Ls = e^r(dasin(s+) - dasin(s-))
        !       Sp = 0.5*e^{2rg}*(e^{2dr}-1)(dasin(s+) - dasin(s-))
        ! ap = 1/Lp*Sum(Ss*bs)    Lp = e^r*dp*dsqrt(1-s^2)
        !       Ss = 0.5*e^{2rg}*(e^{2dr}-1)*dp*dsqrt(1-s^2)
        !
        ! so  as = -0.5*(e^{2dr}-1)/e^r*Sum( e^{2rg}*bp )
        ! and ap = 0.5*(e^{2dr}-1)/e^r*Sum( e^{2rg}*bs )
        er = exp(r)       
        f = 0.5_d*(exp(2.0_d*dr)-1.0)
        do i=2,nr
            as(i,:,:) = as(i-1,:,:) - f*er(i-1)**2*bp(i,2:ns,:)
            ap(i,:,:) = ap(i-1,:,:) + f*er(i-1)**2*bs(i,:,2:np)
        end do
        !$omp parallel
        !$omp do
        do i=2,nr
            as(i,:,:) = as(i,:,:)/er(i)
            ap(i,:,:) = ap(i,:,:)/er(i)
        end do
        !$omp end do
        !$omp end parallel
        
        ! (b) Compute contribution from lower boundary
        ! --------------------------------------------
        ! -- inductive A = curl(psi*e_r)
        allocate(Q(ns-1,ns-1), lam(ns-1), ev(ns-1), clm(ns-1))
        ! - arrays of coordinate factors
        sc = (/ (s(1) + (dble(i) + 0.5_d)*ds, i=0,ns-2) /)
        dasins = asin(s(2:ns)) - asin(s(1:ns-1))        
        dasinsc = asin(sc(2:ns-1)) - asin(sc(1:ns-2))
        sqsc = sqrt(1.0_d - sc**2)
        sqs = sqrt(1.0_d - s**2)
        Vg = 0.0_d
        Vg(2:ns-1) = sqs(2:ns-1)/dasinsc(1:ns-2)/ds
        Uc = dasins(1:ns-1)/sqsc/ds/dp**2
        ! - array of br*r**2 on lower boundary
        br00 = cmplx(br(1,2:ns,2:np), kind=d) * er(1)**2
        ! - fft of br00 in p direction
        call fftn(br00, shape(br00), dim=(/2/))
        br00 = br00/sqrt(dble(np-1))
        ! - term required for m-dependent part of matrix
        mu(1:np-1) = (/ (dble(i), i=0,np-2) /)/dble(np-1)
        ! - order frequencies [+,-] for FFT
        mu((np-1)/2+1:np-1) = mu((np-1)/2+1:np-1) - 1.0_d
        mu = 4.0_d*sin(PI*mu)**2
        ! - loop over azimuthal modes (positive m)
        do m=1,(np-1)/2+1
            ! - prepare tridiagonal matrix (ev = off-diagonal and lam diagonal).
            ! - also prepare identity matrix Q
            ev = -Vg(1:ns-1)
            Q = 0.0_d
            do i=1,ns-1
                lam(i) = Vg(i) + Vg(i+1) + Uc(i)*mu(m)
                Q(i,i) = 1.0_d
            end do
            ! - compute eigenvectors Q_{lm} and eigenvalues lam_{lm}
            call trieig(lam, ev, Q)
            ! - compute array of c_{lm} for each l
            do i=1,ns-1
                clm(i) = dot_product(Q(:,i), br00(:,m))/lam(i)
            end do
            ! - compute entry for this m in psit = Sum_l c_{lm}Q_{lm}^j
            br00(:,m) = matmul(Q,clm)
            if (m > 1) br00(:,np+1-m) = conjg(br00(:,m))
        end do
        ! - compute psi by inverse FFT
        call fftn(br00, shape(br00), dim=(/2/), inv=.true.)
        psi = dble(br00)*sqrt(dble(np-1))
        ! - hence compute as and ap on lower boundary (A0 term)
        !$omp parallel
        !$omp do        
        do i=1,np
            as(1,:,i) = (psi(:,mod(i-3+np,np-1)+1) - psi(:,mod(i-2+np,np-1)+1))/sqsc/dp/er(1)
        end do
        !$omp end do  
        !$omp do        
        do i=1,np-1
            ap(1,2:ns-1,i) = (psi(2:ns-1,i) - psi(1:ns-2,i))/dasinsc(1:ns-2)/er(1)
        end do
        !$omp end do
        !$omp end parallel        
        ! - finally update A0 term at other heights
        !   (ratio of Ls0/Ls, Lp0/Lp = e^r0/e^r
        !$omp parallel
        !$omp do
        do i=2,nr
            as(i,:,:) = as(i,:,:) + as(1,:,:)*er(1)/er(i)
            ap(i,:,:) = ap(i,:,:) + ap(1,:,:)*er(1)/er(i)
        end do
        !$omp end do
        !$omp end parallel

        if (gauge.eq.'minimal') then
            print*,'COMPUTING MINIMAL GAUGE'
            ! (2) Compute Aps on outer boundary (needed to fix boundary condition at pole)
            ! ===========================================================================
            ! [note: gauge is already correct at lower boundary from DeVore-Coulomb]
            ! - array of br on outer boundary
            br00 = cmplx(br(nr,2:ns,2:np)*exp(r(nr)), kind=d)
            ! - fft of br00 in p direction
            call fftn(br00, shape(br00), dim=(/2/))
            br00 = br00/sqrt(dble(np-1))
            ! - loop over azimuthal modes (positive m)
            do m=1,(np-1)/2+1
                ! - prepare tridiagonal matrix (ev = off-diagonal and lam diagonal).
                ! - also prepare identity matrix Q
                ev = -Vg(1:ns-1)
                Q = 0.0_d
                do i=1,ns-1
                    lam(i) = Vg(i) + Vg(i+1) + Uc(i)*mu(m)
                    Q(i,i) = 1.0_d
                end do
                ! - compute eigenvectors Q_{lm} and eigenvalues lam_{lm}
                call trieig(lam, ev, Q)
                ! - compute array of c_{lm} for each l
                do i=1,ns-1
                    clm(i) = dot_product(Q(:,i), br00(:,m))/lam(i)
                end do
                ! - compute entry for this m in psit = Sum_l c_{lm}Q_{lm}^j
                br00(:,m) = matmul(Q,clm)
                if (m > 1) br00(:,np+1-m) = conjg(br00(:,m))
            end do
            ! - compute psi by inverse FFT
            call fftn(br00, shape(br00), dim=(/2/), inv=.true.)
            psi = dble(br00)*sqrt(dble(np-1))
            ! - hence compute as on outer boundary
            allocate(aps1(ns-1,np))
            !$omp parallel
            !$omp do              
            do i=1,np
                aps1(:,i) = (psi(:,mod(i-3+np,np-1)+1) - psi(:,mod(i-2+np,np-1)+1))/sqsc/dp
            end do
            !$omp end do
            !$omp end parallel
            deallocate(Q, lam, ev, clm)
            
            ! (3) Compute gauge function g such that lap_h(g) = -div_h(A) on outer boundary
            ! ============================================================================
            ! [with polar boundary condition dg/ds = Aps - As]
            ! Initialise:
            allocate(g(nr,ns,np))   ! gauge function at grid points
            g = 0.0_d
            allocate(Q(ns-2,ns-2), lam(ns-2), ev(ns-2), clm(ns-2), rhs(ns-2,np-1))
            allocate(Lpc(ns-1), Vc(ns-1), Lsc(ns-1), Ls(ns-2), Ug(ns-2))
            Vc = sqsc*dp/dasins(1:ns-1)
            Ug = dasinsc(1:ns-2)/sqs(2:ns-1)/dp
            Ls = er(nr)*dasinsc(1:ns-2)
            Lpc = er(nr)*sqsc*dp
            Lsc = er(nr)*dasins(1:ns-1)
            ! Solve Poisson equation on top boundary:
            ! - rhs [-div(A)]
            !$omp parallel
            !$omp do               
            do i=2,np-1
                rhs(:,i) = Ls*(-ap(nr,2:ns-1,i) + ap(nr,2:ns-1,i-1)) - Lpc(2:ns-1)*as(nr,2:ns-1,i) + Lpc(1:ns-2)*as(nr,1:ns-2,i)
            end do
            !$omp end do
            !$omp end parallel            
            rhs(:,1) = Ls*(-ap(nr,2:ns-1,1) + ap(nr,2:ns-1,np-1)) - Lpc(2:ns-1)*as(nr,2:ns-1,1) + Lpc(1:ns-2)*as(nr,1:ns-2,1)
            ! - incorporate Neumann boundary conditions [dg/ds = Aps - As]:
            rhs(1,:) = rhs(1,:) + Vc(1)*Lsc(1)*(aps1(1,1:np-1) - as(nr,1,1:np-1))
            rhs(ns-2,:) = rhs(ns-2,:) - Vc(ns-1)*Lsc(ns-1)*(aps1(ns-1,1:np-1) - as(nr,ns-1,1:np-1))
            ! - fft of rhs in p direction
            call fftn(rhs, shape(rhs), dim=(/2/))
            rhs = -rhs/sqrt(dble(np-1))
            ! - term required for m-dependent part of matrix
            mu(1:np-1) = (/ (dble(i), i=0,np-2) /)/dble(np-1)
            ! - order frequencies [+,-] for FFT
            mu((np-1)/2+1:np-1) = mu((np-1)/2+1:np-1) - 1.0_d
            mu = 4.0_d*sin(PI*mu)**2
            ! - loop over azimuthal modes (positive m)
            do m=1,(np-1)/2+1
                ! - prepare tridiagonal matrix (ev = off-diagonal and lam diagonal).
                ! - also prepare identity matrix Q
                ev = -Vc(1:ns-2)
                Q = 0.0_d
                do j=2,ns-3
                    lam(j) = Vc(j) + Vc(j+1) + Ug(j)*mu(m)
                    Q(j,j) = 1.0_d
                end do
                ! - modify boundary terms because of Neumann BC above:
                lam(1) = Vc(2) + Ug(1)*mu(m)
                Q(1,1) = 1.0_d
                lam(ns-2) = Vc(ns-2) + Ug(ns-2)*mu(m)
                Q(ns-2,ns-2) = 1.0_d
                ! - compute eigenvectors Q_{lm} and eigenvalues lam_{lm}
                call trieig(lam, ev, Q)
                ! - compute array of c_{lm} for each l
                do j=1,ns-2
                    clm(j) = dot_product(Q(:,j), rhs(:,m))/lam(j)
                end do
                ! - compute entry for this m in psit = Sum_l c_{lm}Q_{lm}^j
                rhs(:,m) = matmul(Q,clm)
                if (m > 1) rhs(:,np+1-m) = conjg(rhs(:,m))
            end do
            ! - compute psi by inverse FFT
            call fftn(rhs, shape(rhs), dim=(/2/), inv=.true.)
            g(nr,2:ns-1,1:np-1) = dble(rhs)*dsqrt(dble(np-1))
            g(nr,2:ns-1,np) = g(nr,2:ns-1,1)
            ! - fill polar values by matching Aps in pole-most grid point:
            g(nr,1,:) = g(nr,2,:) - Lsc(1)*(aps1(1,:) - as(nr,1,:))
            g(nr,ns,:) = g(nr,ns-1,:) + Lsc(ns-1)*(aps1(ns-1,:) - as(nr,ns-1,:))
            ! - add a constant to fix zero mean [only affects ar, not as or ap]:
            g(nr,:,:) = g(nr,:,:) - sum(g(nr,:,:))/dble(ns*np)
            ! Fill in other heights by linear interpolation:
            !$omp parallel
            !$omp do
            do i=2,nr-1
                g(i,:,:) = (r(i) - r(1))/(r(nr) - r(1))*g(nr,:,:)
            end do
            !$omp end do
            !$omp end parallel
            
            ! (4) Change gauge to A + grad(g)
            ! ==============================
            !$omp parallel private(j)   
            !$omp do
            do i=1,nr-1
                ar(i,:,:) =  ar(i,:,:) + (g(i+1,:,:) - g(i,:,:))/(er(i+1) - er(i))
            end do
            !$omp end do
            !$omp do
            do i=1,nr
                do j=1,ns-1
                    as(i,j,:) = as(i,j,:) + (g(i,j+1,:) - g(i,j,:))/er(i)/dasins(j)
                end do
                do j=2,ns-1
                    ap(i,j,:) = ap(i,j,:) + (g(i,j,2:np) - g(i,j,1:np-1))/er(i)/dp/sqs(j)
                end do
            end do    
            !$omp end do
            !$omp end parallel
        end if
                    
        ! (5) Prepare for interpolation
        ! ============================
        allocate(da(8,nr,ns,np,3))
        call prepareInterpA(ar, as, ap) 
        
        ! (6) Check that A curls to give B (for debugging)
        ! ===============================================
        if (CHECK) then
            allocate(br2(nr,ns+1,np+1), bs2(nr+1,ns,np+1), bp2(nr+1,ns+1,np))
            call computeBfromA(ar, as, ap, br2, bs2, bp2)
            print*,'MAX ERROR br =', maxval(abs(br2-br))
            print*,'MAX ERROR bs =', maxval(abs(bs2-bs))
            print*,'MAX ERROR bp =', maxval(abs(bp2-bp))
            open(unit=1, file='test.unf', form='unformatted')
            write(1) bp(nr,:,:)
            write(1) bp2(nr,:,:)
            close(1)           
        end if

    end subroutine computeA

    !****************************************************************
    subroutine computeP(br)
        real(d), intent(in), dimension(:,:,:) :: br
        ! ---
        ! Compute poloidal potential at all radii, by solving Poisson
        ! equation. Located at face centres.
        ! ---
        real(d), dimension(:,:), allocatable :: Q
        real(d), dimension(:), allocatable :: lam, ev
        complex(d), dimension(:), allocatable :: clm
        integer :: i, j, m
        real(d) :: Vg(ns), sc(ns-1), Uc(ns-1)
        real(d) :: dasins(ns-1), dasinsc(ns-2), sqsc(ns-1), sqs(ns)
        complex(d) :: br00(ns-1,np-1)
        real(d) :: mu(np-1)

        allocate(Q(ns-1,ns-1), lam(ns-1), ev(ns-1), clm(ns-1))
        ! - arrays of coordinate factors
        sc = (/ (s(1) + (dble(i) + 0.5_d)*ds, i=0,ns-2) /)
        dasins = asin(s(2:ns)) - asin(s(1:ns-1))
        dasinsc = asin(sc(2:ns-1)) - asin(sc(1:ns-2))
        sqsc = sqrt(1.0_d - sc**2)
        sqs = sqrt(1.0_d - s**2)
        Vg = 0.0_d
        Vg(2:ns-1) = sqs(2:ns-1)/dasinsc(1:ns-2)/ds
        Uc = dasins(1:ns-1)/sqsc/ds/dp**2
        ! - term required for m-dependent part of matrix
        mu(1:np-1) = (/ (dble(i), i=0,np-2) /)/dble(np-1)
        ! - order frequencies [+,-] for FFT
        mu((np-1)/2+1:np-1) = mu((np-1)/2+1:np-1) - 1.0_d
        mu = 4.0_d*sin(PI*mu)**2
        ! - initialise output array
        allocate(pol(nr,ns-1,np-1))
        pol = 0.0_d
        do j=1,nr
            ! - 2d array of br
            br00 = cmplx(br(j,2:ns,2:np), kind=d)
            ! - fft of br00 in p direction
            call fftn(br00, shape(br00), dim=(/2/))
            br00 = br00/sqrt(dble(np-1))
            ! - loop over azimuthal modes (positive m)
            !$omp parallel private(ev, Q, i, lam, clm)
            !$omp do
            do m=1,(np-1)/2+1
                ! - prepare tridiagonal matrix (ev = off-diagonal and lam diagonal).
                ! - also prepare identity matrix Q
                ev = -Vg(1:ns-1)
                Q = 0.0_d
                do i=1,ns-1
                    lam(i) = Vg(i) + Vg(i+1) + Uc(i)*mu(m)
                    Q(i,i) = 1.0_d
                end do
                ! - compute eigenvectors Q_{lm} and eigenvalues lam_{lm}
                call trieig(lam, ev, Q)
                ! - compute array of c_{lm} for each l
                do i=1,ns-1
                    clm(i) = dot_product(Q(:,i), br00(:,m))/lam(i)
                end do
                ! - compute entry for this m in psit = Sum_l c_{lm}Q_{lm}^j
                br00(:,m) = matmul(Q,clm)
                if (m > 1) br00(:,np+1-m) = conjg(br00(:,m))
            end do
            !$omp end do
            !$omp end parallel
            ! - compute P by inverse FFT
            call fftn(br00, shape(br00), dim=(/2/), inv=.true.)
            pol(j,:,:) = dble(br00)*sqrt(dble(np-1))
            pol(j,:,:) = pol(j,:,:)*exp(r(j))**2
        end do

    end subroutine computeP

    !****************************************************************
    subroutine prepareInterpA(ar, as, ap)
        real(d), intent(in), dimension(:,:,:) :: ar, as, ap
        ! ---
        ! Convert input vector field on cell edges to cartesian components at
        ! grid points. Then precompute derivatives and read in matrix for
        ! tricubic interpolation.
        ! - inputs: ar(nr,ns,np), as(nr,ns,np), ap(nr,ns,np) on edges
        ! - store in module-wide: da(8,nr,ns,np,3) -- ax, ay, az and derivatives, at grid pts.
        !                         m(64,64) -- interpolation matrix
        ! ---
        integer :: nr,ns,np,i,j
        real(d):: dr, ds, dp

        nr = size(r)
        ns = size(s)
        np = size(p)

        dr = r(2) - r(1)
        ds = s(2) - s(1)
        dp = p(2) - p(1)

        ! Average to grid pts.
        ! --------------------
        if (.not.allocated(arg)) &
            allocate(arg(nr,ns,np),asg(nr,ns,np),apg(nr,ns,np))
        arg = 0.0_d
        asg = 0.0_d
        apg = 0.0_d
        arg(2:nr-1,:,:) = 0.5_d*(ar(1:nr-2,:,:) + ar(2:nr-1,:,:))
        asg(:,2:ns-1,:) = 0.5_d*(as(:,1:ns-2,:) + as(:,2:ns-1,:))
        apg(:,:,2:np-1) = 0.5_d*(ap(:,:,1:np-2) + ap(:,:,2:np-1))
        apg(:,:,1) = 0.5_d*(ap(:,:,np-1) + ap(:,:,1))
        apg(:,:,np) = apg(:,:,1)
        arg(1,:,:) = 2.0_d*arg(2,:,:) - arg(3,:,:)
        arg(nr,:,:) = 2.0_d*arg(nr-1,:,:) - arg(nr-2,:,:)

        ! Convert to cartesian cmpts, and estimate derivatives at grid points.
        ! --------------------------------------------------------------------
        da = 0.0_d
        ! Function itself:
        do i=1,np
            do j=1,ns
                da(1,:,j,i,1) = dsqrt(1.0_d - s(j)**2)*cos(p(i))*arg(:,j,i) - s(j)*cos(p(i))*asg(:,j,i) - sin(p(i))*apg(:,j,i)
                da(1,:,j,i,2) = dsqrt(1.0_d - s(j)**2)*sin(p(i))*arg(:,j,i) - s(j)*sin(p(i))*asg(:,j,i) + cos(p(i))*apg(:,j,i)
                da(1,:,j,i,3) = s(j)*arg(:,j,i) + dsqrt(1.0_d - s(j)**2)*asg(:,j,i)
            end do
        end do
        if (.not.store_glb) deallocate(arg,asg,apg)
        ! d/dr: (1-sided for boundaries)
        da(2,2:nr-1,:,:,:) = 0.5_d*(da(1,3:nr,:,:,:) - da(1,1:nr-2,:,:,:))
        da(2,1,:,:,:) = 0.5_d*(-3.0_d*da(1,1,:,:,:) + 4.0_d*da(1,2,:,:,:) - da(1,3,:,:,:))
        da(2,nr,:,:,:) = 0.5_d*(3.0_d*da(1,nr,:,:,:) - 4.0_d*da(1,nr-1,:,:,:) + da(1,nr-2,:,:,:))
        ! d/ds: (1-sided for boundaries)
        da(3,:,2:ns-1,:,:) = 0.5_d*(da(1,:,3:ns,:,:) - da(1,:,1:ns-2,:,:))
        da(3,:,1,:,:) = 0.5_d*(-3.0_d*da(1,:,1,:,:) + 4.0_d*da(1,:,2,:,:) - da(1,:,3,:,:))
        da(3,:,ns,:,:) = 0.5_d*(3.0_d*da(1,:,ns,:,:) - 4.0_d*da(1,:,ns-1,:,:) + da(1,:,ns-2,:,:))
        ! d/dp: (periodic for boundaries)
        da(4,:,:,2:np-1,:) = 0.5_d*(da(1,:,:,3:np,:) - da(1,:,:,1:np-2,:))
        da(4,:,:,1,:) = 0.5_d*(da(1,:,:,2,:) - da(1,:,:,np-1,:))
        da(4,:,:,np,:) = 0.5_d*(da(1,:,:,2,:) - da(1,:,:,np-1,:))
        ! d^2/drds: (1-sided on boundaries)
        da(5,2:nr-1,:,:,:) = 0.5_d*(da(3,3:nr,:,:,:) - da(3,1:nr-2,:,:,:))
        da(5,1,:,:,:) = 0.5_d*(-3.0_d*da(3,1,:,:,:) + 4.0_d*da(3,2,:,:,:) - da(3,3,:,:,:))
        da(5,nr,:,:,:) = 0.5_d*(3.0_d*da(3,nr,:,:,:) - 4.0_d*da(3,nr-1,:,:,:) + da(3,nr-2,:,:,:))
        ! d^2/drdp: (1-sided on boundaries)
        da(6,2:nr-1,:,:,:) = 0.5_d*(da(4,3:nr,:,:,:) - da(4,1:nr-2,:,:,:))
        da(6,1,:,:,:) = 0.5_d*(-3.0_d*da(4,1,:,:,:) + 4.0_d*da(4,2,:,:,:) - da(4,3,:,:,:))
        da(6,nr,:,:,:) = 0.5_d*(3.0_d*da(4,nr,:,:,:) - 4.0_d*da(4,nr-1,:,:,:) + da(4,nr-2,:,:,:))
        ! d^2/dsdp: (1-sided on boundaries)
        da(7,:,2:ns-1,:,:) = 0.5_d*(da(4,:,3:ns,:,:) - da(4,:,1:ns-2,:,:))
        da(7,:,1,:,:) = 0.5_d*(-3.0_d*da(4,:,1,:,:) + 4.0_d*da(4,:,2,:,:) - da(4,:,3,:,:))
        da(7,:,ns,:,:) = 0.5_d*(3.0_d*da(4,:,ns,:,:) - 4.0_d*da(4,:,ns-1,:,:) + da(4,:,ns-2,:,:))
        ! d^3/drdsdp: (1-sided on boundaries)
        da(8,2:nr-1,:,:,:) = 0.5_d*(da(7,3:nr,:,:,:) - da(7,1:nr-2,:,:,:))
        da(8,1,:,:,:) = 0.5_d*(-3.0_d*da(7,1,:,:,:) + 4.0_d*da(7,2,:,:,:) - da(7,3,:,:,:))
        da(8,nr,:,:,:) = 0.5_d*(3.0_d*da(7,nr,:,:,:) - 4.0_d*da(7,nr-1,:,:,:) + da(7,nr-2,:,:,:))
        
        call interpMatrix
        
    end subroutine prepareInterpA
    
    !****************************************************************
    subroutine prepareInterpB(br, bs, bp)
        real(d), intent(in), dimension(:,:,:) :: br, bs, bp
        ! ---
        ! Convert input vector field on cell faces to cartesian components at
        ! grid points. Then precompute derivatives and read in matrix for
        ! tricubic interpolation.
        ! - inputs: br(nr,ns+1,np+1), bs(nr+1,ns,np+1), bp(nr+1,ns+1,np)
        !             on cell faces
        ! - store in module-wide: db(8,nr,ns,np,3) -- bx, by, bz and derivatives, at grid pts.
        !                         m(64,64) -- interpolation matrix
        ! ---
        real(d), dimension(:,:), allocatable :: atmp
        real(d), dimension(:), allocatable :: rc, sc
        integer :: i,j

        ! Weighted average to grid points
        ! -------------------------------
        ! Coordinates at face centres:
        allocate(rc(nr+1),sc(ns+1))
        rc = (/ (r(1) + (dble(i) - 0.5_d)*dr, i=0,nr) /)
        sc = (/ (s(1) + (dble(i) - 0.5_d)*ds, i=0,ns) /)
        if (.not.allocated(brg)) &
            allocate(brg(nr,ns,np),bsg(nr,ns,np),bpg(nr,ns,np))
        ! br (no need for weighting since cells have equal area)
        brg = 0.25_d*(br(:,1:ns,1:np) + br(:,1:ns,2:np+1) &
            + br(:,2:ns+1,1:np) + br(:,2:ns+1,2:np+1))
        ! bs
        allocate(atmp(nr+1,ns))
        do i=1,nr+1
        atmp(i,:) = 0.5_d*(dexp(2.0_d*dr) - 1.0_d)*dp*dexp(2.0_d*rc(i) - dr) &
                *dsqrt(1.0_d - s**2)
        end do
        atmp(:,1) = atmp(:,2)
        atmp(:,ns) = atmp(:,ns-1)
        do i=1,np
        bsg(:,:,i) = (bs(1:nr,:,i) + bs(1:nr,:,i+1))*atmp(1:nr,:) &
                + (bs(2:nr+1,:,i) + bs(2:nr+1,:,i+1))*atmp(2:nr+1,:)
        bsg(:,:,i) = 0.5*bsg(:,:,i)/(atmp(1:nr,:) + atmp(2:nr+1,:))
        end do 
        deallocate(atmp)
        ! bp
        allocate(atmp(nr+1,ns+1))
        do i=1,nr+1
        atmp(i,2:ns) = 0.5_d*(dexp(2.0_d*dr) - 1.0_d)*dexp(2.0_d*rc(i) - dr) &
                *(dasin(s(2:ns)) - dasin(s(1:ns-1)))
        end do
        atmp(:,1) = atmp(:,2)
        atmp(:,ns+1) = atmp(:,ns)
        do i=1,np
        bpg(:,:,i) = bp(1:nr,1:ns,i)*atmp(1:nr,1:ns) &
                + bp(1:nr,2:ns+1,i)*atmp(1:nr,2:ns+1) &
                + bp(2:nr+1,1:ns,i)*atmp(2:nr+1,1:ns) &
                + bp(2:nr+1,2:ns+1,i)*atmp(2:nr+1,2:ns+1)
        bpg(:,:,i) = bpg(:,:,i)/(atmp(1:nr,1:ns) + atmp(1:nr,2:ns+1) &
                + atmp(2:nr+1,1:ns) + atmp(2:nr+1,2:ns+1))
        end do
        deallocate(atmp)

        ! Convert to cartesian cmpts at grid points, and estimate derivatives
        ! -------------------------------------------------------------------
        db = 0.0_d
        ! Function itself:
        do i=1,np
            do j=1,ns
                db(1,:,j,i,1) = dsqrt(1.0_d - s(j)**2)*cos(p(i))*brg(:,j,i) - s(j)*cos(p(i))*bsg(:,j,i) - sin(p(i))*bpg(:,j,i)
                db(1,:,j,i,2) = dsqrt(1.0_d - s(j)**2)*sin(p(i))*brg(:,j,i) - s(j)*sin(p(i))*bsg(:,j,i) + cos(p(i))*bpg(:,j,i)
                db(1,:,j,i,3) = s(j)*brg(:,j,i) + dsqrt(1.0_d - s(j)**2)*bsg(:,j,i)
            end do
        end do
        if (.not.store_glb) deallocate(brg,bsg,bpg)
        ! d/dr: (1-sided for boundaries)
        db(2,2:nr-1,:,:,:) = 0.5_d*(db(1,3:nr,:,:,:) - db(1,1:nr-2,:,:,:))
        db(2,1,:,:,:) = 0.5_d*(-3.0_d*db(1,1,:,:,:) + 4.0_d*db(1,2,:,:,:) - db(1,3,:,:,:))
        db(2,nr,:,:,:) = 0.5_d*(3.0_d*db(1,nr,:,:,:) - 4.0_d*db(1,nr-1,:,:,:) + db(1,nr-2,:,:,:))
        ! d/ds: (1-sided for boundaries)
        db(3,:,2:ns-1,:,:) = 0.5_d*(db(1,:,3:ns,:,:) - db(1,:,1:ns-2,:,:))
        db(3,:,1,:,:) = 0.5_d*(-3.0_d*db(1,:,1,:,:) + 4.0_d*db(1,:,2,:,:) - db(1,:,3,:,:))
        db(3,:,ns,:,:) = 0.5_d*(3.0_d*db(1,:,ns,:,:) - 4.0_d*db(1,:,ns-1,:,:) + db(1,:,ns-2,:,:))
        ! d/dp: (periodic for boundaries)
        db(4,:,:,2:np-1,:) = 0.5_d*(db(1,:,:,3:np,:) - db(1,:,:,1:np-2,:))
        db(4,:,:,1,:) = 0.5_d*(db(1,:,:,2,:) - db(1,:,:,np-1,:))
        db(4,:,:,np,:) = 0.5_d*(db(1,:,:,2,:) - db(1,:,:,np-1,:))
        ! d^2/drds: (1-sided on boundaries)
        db(5,2:nr-1,:,:,:) = 0.5_d*(db(3,3:nr,:,:,:) - db(3,1:nr-2,:,:,:))
        db(5,1,:,:,:) = 0.5_d*(-3.0_d*db(3,1,:,:,:) + 4.0_d*db(3,2,:,:,:) - db(3,3,:,:,:))
        db(5,nr,:,:,:) = 0.5_d*(3.0_d*db(3,nr,:,:,:) - 4.0_d*db(3,nr-1,:,:,:) + db(3,nr-2,:,:,:))
        ! d^2/drdp: (1-sided on boundaries)
        db(6,2:nr-1,:,:,:) = 0.5_d*(db(4,3:nr,:,:,:) - db(4,1:nr-2,:,:,:))
        db(6,1,:,:,:) = 0.5_d*(-3.0_d*db(4,1,:,:,:) + 4.0_d*db(4,2,:,:,:) - db(4,3,:,:,:))
        db(6,nr,:,:,:) = 0.5_d*(3.0_d*db(4,nr,:,:,:) - 4.0_d*db(4,nr-1,:,:,:) + db(4,nr-2,:,:,:))
        ! d^2/dsdp: (1-sided on boundaries)
        db(7,:,2:ns-1,:,:) = 0.5_d*(db(4,:,3:ns,:,:) - db(4,:,1:ns-2,:,:))
        db(7,:,1,:,:) = 0.5_d*(-3.0_d*db(4,:,1,:,:) + 4.0_d*db(4,:,2,:,:) - db(4,:,3,:,:))
        db(7,:,ns,:,:) = 0.5_d*(3.0_d*db(4,:,ns,:,:) - 4.0_d*db(4,:,ns-1,:,:) + db(4,:,ns-2,:,:))
        ! d^3/drdsdp: (1-sided on boundaries)
        db(8,2:nr-1,:,:,:) = 0.5_d*(db(7,3:nr,:,:,:) - db(7,1:nr-2,:,:,:))
        db(8,1,:,:,:) = 0.5_d*(-3.0_d*db(7,1,:,:,:) + 4.0_d*db(7,2,:,:,:) - db(7,3,:,:,:))
        db(8,nr,:,:,:) = 0.5_d*(3.0_d*db(7,nr,:,:,:) - 4.0_d*db(7,nr-1,:,:,:) + db(7,nr-2,:,:,:))

        call interpMatrix
        
    end subroutine prepareInterpB   

     !****************************************************************
    subroutine prepareInterpJ(jr, js, jp)
        real(d), intent(in), dimension(:,:,:) :: jr, js, jp
        ! ---
        ! For input vector field at grid points, precompute derivatives
        ! and read in matrix for tricubic interpolation.
        ! - inputs: jr(nr,ns,np), js(nr,ns,np), jp(nr,ns,np)
        !             on cell faces
        ! - store in module-wide: dj(8,nr,ns,np,3) -- jx, jy, jz and derivatives, at grid pts.
        !                         m(64,64) -- interpolation matrix
        ! ---
        integer :: i,j

        ! Convert to cartesian cmpts at grid points, and estimate derivatives
        ! -------------------------------------------------------------------
        dj = 0.0_d
        ! Function itself:
        do i=1,np
            do j=1,ns
                dj(1,:,j,i,1) = dsqrt(1.0_d - s(j)**2)*cos(p(i))*jr(:,j,i) - s(j)*cos(p(i))*js(:,j,i) - sin(p(i))*jp(:,j,i)
                dj(1,:,j,i,2) = dsqrt(1.0_d - s(j)**2)*sin(p(i))*jr(:,j,i) - s(j)*sin(p(i))*js(:,j,i) + cos(p(i))*jp(:,j,i)
                dj(1,:,j,i,3) = s(j)*jr(:,j,i) + dsqrt(1.0_d - s(j)**2)*js(:,j,i)
            end do
        end do
        ! d/dr: (1-sided for boundaries)
        dj(2,2:nr-1,:,:,:) = 0.5_d*(dj(1,3:nr,:,:,:) - dj(1,1:nr-2,:,:,:))
        dj(2,1,:,:,:) = 0.5_d*(-3.0_d*dj(1,1,:,:,:) + 4.0_d*dj(1,2,:,:,:) - dj(1,3,:,:,:))
        dj(2,nr,:,:,:) = 0.5_d*(3.0_d*dj(1,nr,:,:,:) - 4.0_d*dj(1,nr-1,:,:,:) + dj(1,nr-2,:,:,:))
        ! d/ds: (1-sided for boundaries)
        dj(3,:,2:ns-1,:,:) = 0.5_d*(dj(1,:,3:ns,:,:) - dj(1,:,1:ns-2,:,:))
        dj(3,:,1,:,:) = 0.5_d*(-3.0_d*dj(1,:,1,:,:) + 4.0_d*dj(1,:,2,:,:) - dj(1,:,3,:,:))
        dj(3,:,ns,:,:) = 0.5_d*(3.0_d*dj(1,:,ns,:,:) - 4.0_d*dj(1,:,ns-1,:,:) + dj(1,:,ns-2,:,:))
        ! d/dp: (periodic for boundaries)
        dj(4,:,:,2:np-1,:) = 0.5_d*(dj(1,:,:,3:np,:) - dj(1,:,:,1:np-2,:))
        dj(4,:,:,1,:) = 0.5_d*(dj(1,:,:,2,:) - dj(1,:,:,np-1,:))
        dj(4,:,:,np,:) = 0.5_d*(dj(1,:,:,2,:) - dj(1,:,:,np-1,:))
        ! d^2/drds: (1-sided on boundaries)
        dj(5,2:nr-1,:,:,:) = 0.5_d*(dj(3,3:nr,:,:,:) - dj(3,1:nr-2,:,:,:))
        dj(5,1,:,:,:) = 0.5_d*(-3.0_d*dj(3,1,:,:,:) + 4.0_d*dj(3,2,:,:,:) - dj(3,3,:,:,:))
        dj(5,nr,:,:,:) = 0.5_d*(3.0_d*dj(3,nr,:,:,:) - 4.0_d*dj(3,nr-1,:,:,:) + dj(3,nr-2,:,:,:))
        ! d^2/drdp: (1-sided on boundaries)
        dj(6,2:nr-1,:,:,:) = 0.5_d*(dj(4,3:nr,:,:,:) - dj(4,1:nr-2,:,:,:))
        dj(6,1,:,:,:) = 0.5_d*(-3.0_d*dj(4,1,:,:,:) + 4.0_d*dj(4,2,:,:,:) - dj(4,3,:,:,:))
        dj(6,nr,:,:,:) = 0.5_d*(3.0_d*dj(4,nr,:,:,:) - 4.0_d*dj(4,nr-1,:,:,:) + dj(4,nr-2,:,:,:))
        ! d^2/dsdp: (1-sided on boundaries)
        dj(7,:,2:ns-1,:,:) = 0.5_d*(dj(4,:,3:ns,:,:) - dj(4,:,1:ns-2,:,:))
        dj(7,:,1,:,:) = 0.5_d*(-3.0_d*dj(4,:,1,:,:) + 4.0_d*dj(4,:,2,:,:) - dj(4,:,3,:,:))
        dj(7,:,ns,:,:) = 0.5_d*(3.0_d*dj(4,:,ns,:,:) - 4.0_d*dj(4,:,ns-1,:,:) + dj(4,:,ns-2,:,:))
        ! d^3/drdsdp: (1-sided on boundaries)
        dj(8,2:nr-1,:,:,:) = 0.5_d*(dj(7,3:nr,:,:,:) - dj(7,1:nr-2,:,:,:))
        dj(8,1,:,:,:) = 0.5_d*(-3.0_d*dj(7,1,:,:,:) + 4.0_d*dj(7,2,:,:,:) - dj(7,3,:,:,:))
        dj(8,nr,:,:,:) = 0.5_d*(3.0_d*dj(7,nr,:,:,:) - 4.0_d*dj(7,nr-1,:,:,:) + dj(7,nr-2,:,:,:))
        
        call interpMatrix
        
    end subroutine prepareInterpJ
    
    !****************************************************************
    subroutine interpMatrix
        ! ---
        ! Initialize interpolation matrix m for tricubic interpolation.
        ! --

        m(1, :) = (/1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(2, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(3, :) = (/-3, 3, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(4, :) = (/2, -2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(5, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(6, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(7, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, &
        3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, &
        -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(8, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, &
        -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, &
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(9, :) = (/-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, &
        0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(10, :) = (/0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, &
        0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(11, :) = (/9, -9, -9, 9, 0, 0, 0, 0, 6, 3, -6, -3, 0, 0, 0, 0, 6, &
        -6, 3, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, &
        2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(12, :) = (/-6, 6, 6, -6, 0, 0, 0, 0, -3, -3, 3, 3, 0, 0, 0, 0, -4, &
        4, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, &
        -2, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(13, :) = (/2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, &
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(14, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, &
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(15, :) = (/-6, 6, 6, -6, 0, 0, 0, 0, -4, -2, 4, 2, 0, 0, 0, 0, -3, &
        3, -3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, &
        -1, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(16, :) = (/4, -4, -4, 4, 0, 0, 0, 0, 2, 2, -2, -2, 0, 0, 0, 0, 2, &
        -2, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, &
        1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(17, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(18, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(19, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(20, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(21, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(22, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0/)
        m(23, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, &
        3, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0/)
        m(24, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, &
        -2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0/)
        m(25, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, &
        0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(26, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0/)
        m(27, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 9, -9, -9, 9, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 6, 3, -6, -3, 0, 0, 0, 0, 6, &
        -6, 3, -3, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0/)
        m(28, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -6, 6, 6, -6, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -3, -3, 3, 3, 0, 0, 0, 0, -4, &
        4, -2, 2, 0, 0, 0, 0, -2, -2, -1, -1, 0, 0, 0, 0/)
        m(29, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, &
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(30, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0/)
        m(31, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -6, 6, 6, -6, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -4, -2, 4, 2, 0, 0, 0, 0, -3, &
        3, -3, 3, 0, 0, 0, 0, -2, -1, -2, -1, 0, 0, 0, 0/)
        m(32, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 4, -4, -4, 4, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 2, 2, -2, -2, 0, 0, 0, 0, 2, &
        -2, 2, -2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0/)
        m(33, :) = (/-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(34, :) = (/0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0, 0, 3, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(35, :) = (/9, -9, 0, 0, -9, 9, 0, 0, 6, 3, 0, 0, -6, -3, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 3, -3, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(36, :) = (/-6, 6, 0, 0, 6, -6, 0, 0, -3, -3, 0, 0, 3, 3, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -4, 4, 0, 0, -2, 2, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -2, -2, 0, 0, -1, -1, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(37, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, &
        0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, &
        0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(38, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, &
        0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, 0, 0, 0/)
        m(39, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, &
        -9, 0, 0, -9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, &
        3, 0, 0, -6, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, &
        -6, 0, 0, 3, -3, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0/)
        m(40, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, &
        6, 0, 0, 6, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, &
        -3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, &
        4, 0, 0, -2, 2, 0, 0, -2, -2, 0, 0, -1, -1, 0, 0/)
        m(41, :) = (/9, 0, -9, 0, -9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, &
        0, 3, 0, -6, 0, -3, 0, 6, 0, -6, 0, 3, 0, -3, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, &
        0, 2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(42, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 9, 0, -9, 0, -9, 0, 9, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, &
        0, 3, 0, -6, 0, -3, 0, 6, 0, -6, 0, 3, 0, -3, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0/)
        m(43, :) = (/-27, 27, 27, -27, 27, -27, -27, 27, -18, -9, 18, 9, 18, 9, -18, -9, -18, &
        18, -9, 9, 18, -18, 9, -9, -18, 18, 18, -18, -9, 9, 9, -9, -12, &
        -6, -6, -3, 12, 6, 6, 3, -12, -6, 12, 6, -6, -3, 6, 3, -12, &
        12, -6, 6, -6, 6, -3, 3, -8, -4, -4, -2, -4, -2, -2, -1/)
        m(44, :) = (/18, -18, -18, 18, -18, 18, 18, -18, 9, 9, -9, -9, -9, -9, 9, 9, 12, &
        -12, 6, -6, -12, 12, -6, 6, 12, -12, -12, 12, 6, -6, -6, 6, 6, &
        6, 3, 3, -6, -6, -3, -3, 6, 6, -6, -6, 3, 3, -3, -3, 8, &
        -8, 4, -4, 4, -4, 2, -2, 4, 4, 2, 2, 2, 2, 1, 1/)
        m(45, :) = (/-6, 0, 6, 0, 6, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, &
        0, -3, 0, 3, 0, 3, 0, -4, 0, 4, 0, -2, 0, 2, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, &
        0, -2, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(46, :) = (/0, 0, 0, 0, 0, 0, 0, 0, -6, 0, 6, 0, 6, 0, -6, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, &
        0, -3, 0, 3, 0, 3, 0, -4, 0, 4, 0, -2, 0, 2, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -2, 0, -2, 0, -1, 0, -1, 0/)
        m(47, :) = (/18, -18, -18, 18, -18, 18, 18, -18, 12, 6, -12, -6, -12, -6, 12, 6, 9, &
        -9, 9, -9, -9, 9, -9, 9, 12, -12, -12, 12, 6, -6, -6, 6, 6, &
        3, 6, 3, -6, -3, -6, -3, 8, 4, -8, -4, 4, 2, -4, -2, 6, &
        -6, 6, -6, 3, -3, 3, -3, 4, 2, 4, 2, 2, 1, 2, 1/)
        m(48, :) = (/-12, 12, 12, -12, 12, -12, -12, 12, -6, -6, 6, 6, 6, 6, -6, -6, -6, &
        6, -6, 6, 6, -6, 6, -6, -8, 8, 8, -8, -4, 4, 4, -4, -3, &
        -3, -3, -3, 3, 3, 3, 3, -4, -4, 4, 4, -2, -2, 2, 2, -4, &
        4, -4, 4, -2, 2, -2, 2, -2, -2, -2, -2, -1, -1, -1, -1/)
        m(49, :) = (/2, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(50, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, -2, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(51, :) = (/-6, 6, 0, 0, 6, -6, 0, 0, -4, -2, 0, 0, 4, 2, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -3, 3, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, -2, -1, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(52, :) = (/4, -4, 0, 0, -4, 4, 0, 0, 2, 2, 0, 0, -2, -2, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 2, -2, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(53, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, &
        0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, &
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(54, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, &
        0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0/)
        m(55, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, &
        6, 0, 0, 6, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, &
        -2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, &
        3, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0, -2, -1, 0, 0/)
        m(56, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, &
        -4, 0, 0, -4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, &
        2, 0, 0, -2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, &
        -2, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0/)
        m(57, :) = (/-6, 0, 6, 0, 6, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, &
        0, -2, 0, 4, 0, 2, 0, -3, 0, 3, 0, -3, 0, 3, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, &
        0, -1, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(58, :) = (/0, 0, 0, 0, 0, 0, 0, 0, -6, 0, 6, 0, 6, 0, -6, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, &
        0, -2, 0, 4, 0, 2, 0, -3, 0, 3, 0, -3, 0, 3, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, -2, 0, -1, 0/)
        m(59, :) = (/18, -18, -18, 18, -18, 18, 18, -18, 12, 6, -12, -6, -12, -6, 12, 6, 12, &
        -12, 6, -6, -12, 12, -6, 6, 9, -9, -9, 9, 9, -9, -9, 9, 8, &
        4, 4, 2, -8, -4, -4, -2, 6, 3, -6, -3, 6, 3, -6, -3, 6, &
        -6, 3, -3, 6, -6, 3, -3, 4, 2, 2, 1, 4, 2, 2, 1/)
        m(60, :) = (/-12, 12, 12, -12, 12, -12, -12, 12, -6, -6, 6, 6, 6, 6, -6, -6, -8, &
        8, -4, 4, 8, -8, 4, -4, -6, 6, 6, -6, -6, 6, 6, -6, -4, &
        -4, -2, -2, 4, 4, 2, 2, -3, -3, 3, 3, -3, -3, 3, 3, -4, &
        4, -2, 2, -4, 4, -2, 2, -2, -2, -1, -1, -2, -2, -1, -1/)
        m(61, :) = (/4, 0, -4, 0, -4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, &
        0, 2, 0, -2, 0, -2, 0, 2, 0, -2, 0, 2, 0, -2, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, &
        0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        m(62, :) = (/0, 0, 0, 0, 0, 0, 0, 0, 4, 0, -4, 0, -4, 0, 4, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, &
        0, 2, 0, -2, 0, -2, 0, 2, 0, -2, 0, 2, 0, -2, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0/)
        m(63, :) = (/-12, 12, 12, -12, 12, -12, -12, 12, -8, -4, 8, 4, 8, 4, -8, -4, -6, &
        6, -6, 6, 6, -6, 6, -6, -6, 6, 6, -6, -6, 6, 6, -6, -4, &
        -2, -4, -2, 4, 2, 4, 2, -4, -2, 4, 2, -4, -2, 4, 2, -3, &
        3, -3, 3, -3, 3, -3, 3, -2, -1, -2, -1, -2, -1, -2, -1/)
        m(64, :) = (/8, -8, -8, 8, -8, 8, 8, -8, 4, 4, -4, -4, -4, -4, 4, 4, 4, &
        -4, 4, -4, -4, 4, -4, 4, 4, -4, -4, 4, 4, -4, -4, 4, 2, &
        2, 2, 2, -2, -2, -2, -2, 2, 2, -2, -2, 2, 2, -2, -2, 2, &
        -2, 2, -2, 2, -2, 2, -2, 1, 1, 1, 1, 1, 1, 1, 1/)

    end subroutine interpMatrix

    !****************************************************************
    subroutine readDims(filename)
        character*(*), intent(in) :: filename
        ! ---
        ! Get coordinate dimensions from netcdf file.
        ! (these are the number of grid points in each direction)
        ! - input: filename of netcdf file, incl. path.
        ! - store in shared variables: nr, ns, np -> number of grid pts in each direction
        ! ---
        integer :: ncid,dimid
        
        call chk( nf90_open(filename, NF90_NOWRITE, ncid) )
        
        call chk( nf90_inq_dimid(ncid, "r", dimid) )
        call chk( nf90_inquire_dimension(ncid, dimid, len = nr) )
        call chk( nf90_inq_dimid(ncid, "th", dimid) )
        call chk( nf90_inquire_dimension(ncid, dimid, len = ns) )
        call chk( nf90_inq_dimid(ncid, "ph", dimid) )
        call chk( nf90_inquire_dimension(ncid, dimid, len = np) )
        
        call chk( nf90_close(ncid) )
        
    end subroutine readDims
    
    !****************************************************************
    subroutine readCoords(filename)
        character*(*), intent(in) :: filename
        ! ---
        ! Read coordinate arrays (at grid pts) from netcdf file,
        ! and convert to rho, s, phi coordinates.
        ! - input: filename of netcdf file, incl. path.
        ! - store in shared variables: r(nr), s(ns), p(np) -> rho, s=cos(th), phi arrays at grid pts.
        ! ---
        integer :: ncid,varid
        
        if (allocated(r)) deallocate(r)
        if (allocated(s)) deallocate(s)
        if (allocated(p)) deallocate(p)
        allocate(r(nr),s(ns),p(np))
        
        call chk( nf90_open(filename, NF90_NOWRITE, ncid) )
        
        call chk( nf90_inq_varid(ncid, "r", varid) )
        call chk( nf90_get_var(ncid, varid, r) )
        r = log(r)
        
        call chk( nf90_inq_varid(ncid, "th", varid) )
        call chk( nf90_get_var(ncid, varid, s) )
        s = cos(s)
        
        call chk( nf90_inq_varid(ncid, "ph", varid) )
        call chk( nf90_get_var(ncid, varid, p) )
        
        call chk( nf90_close(ncid) )
        
        ! Store grid spacing:
        dr = r(2) - r(1)
        ds = s(2) - s(1)
        dp = p(2) - p(1)
        
    end subroutine readCoords
    
    !****************************************************************
    subroutine readvec_faces(filename,br,bs,bp)
        character*(*), intent(in) :: filename
        real(d), dimension(:,:,:), intent(inout) :: br,bs,bp
        ! ---
        ! Read vector field on cell faces (incl. ghost cells) from netcdf file.
        ! - input: filename of netcdf file (incl. path)
        ! - output: br(nr,ns+1,np+1), bs(nr+1,ns,np+1), bp(nr+1,ns+1,np)
        !           on cell faces
        ! ---      
        integer :: ncid,varid
        
        call chk( nf90_open(filename, NF90_NOWRITE, ncid) )
        
        call chk( nf90_inq_varid(ncid, "br", varid) )
        call chk( nf90_get_var(ncid, varid, br) )
        
        call chk( nf90_inq_varid(ncid, "bth", varid) )
        call chk( nf90_get_var(ncid, varid, bs) )
        bs = -bs
        
        call chk( nf90_inq_varid(ncid, "bph", varid) )
        call chk( nf90_get_var(ncid, varid, bp) )
        
        call chk( nf90_close(ncid) )

    end subroutine readvec_faces
    
    !****************************************************************
    subroutine readvec_gpts(filename,field,vr,vs,vp)
        character*(*), intent(in) :: filename, field
        real(d), dimension(:,:,:), intent(inout) :: vr,vs,vp
        ! ---
        ! Read vector field defined at grid points from netcdf file.
        ! - input: filename of netcdf file (incl. path)
        !          field -- name of variable in file (e.g. 'j' or 'e')
        ! - output: vr(nr,ns,np), vs(nr,ns,np), vp(nr,ns,np) at grid pts
        ! ---     
        integer :: ncid,varid
        
        call chk( nf90_open(filename, NF90_NOWRITE, ncid) )
        
        call chk( nf90_inq_varid(ncid, field//"r", varid) )
        call chk( nf90_get_var(ncid, varid, vr) )
        
        call chk( nf90_inq_varid(ncid, field//"th", varid) )
        call chk( nf90_get_var(ncid, varid, vs) )
        vs = -vs
        
        call chk( nf90_inq_varid(ncid, field//"ph", varid) )
        call chk( nf90_get_var(ncid, varid, vp) )
        
        call chk( nf90_close(ncid) )

    end subroutine readvec_gpts

    !****************************************************************
    subroutine writeP(filename)
        character*(*), intent(in) :: filename
        ! ---
        ! Write poloidal vector potential to unformatted file.
        ! --
        real(d):: sc(ns-1), pc(np-1)
        sc = 0.5*(s(2:ns) + s(1:ns-1))
        pc = 0.5*(p(2:np) + p(1:np-1))

        open(unit=1, file=filename, form='unformatted')
        write(1) r
        write(1) sc
        write(1) pc
        write(1) pol*RSUN*RSUN
        close(1)

    end subroutine writeP

end module readnc
