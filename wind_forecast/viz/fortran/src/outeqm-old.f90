!****************************************************************
module outeqm
!****************************************************************

    use shared
    use readnc
    use ffts
    use netcdf
    
    implicit none

    private
    public :: compute_outeqm_inner
    
    contains
  
    !****************************************************************
    subroutine compute_outeqm_inner(wkpath, snap, v1) 
        character*(*), intent(in) :: wkpath, snap
        real(d), intent(in) :: v1
        ! ---
        ! Compute outflow equilibrium in original domain
        ! (with source surface at top).
        ! Parameter v1 is the speed parameter [km/s]
        ! Save B on cell faces to new netcdf file.
        ! ---       
        real(d) :: br(nr,ns+1,np+1), bs(nr+1,ns,np+1), bp(nr+1,ns+1,np)

        if (VERBOSE) print*,'Computing outflow eqm...'
        call compute_outeqm(br0, v1, br, bs, bp)
        print*, shape(br)
        print*,'oflux=',sum(abs(br(nr,2:ns,2:np)))
        call b_to_netcdf(br, bs, bp, wkpath, 'outeqm_'//snap)
          
    end subroutine compute_outeqm_inner

!-----------------------------------------------------------------------

    !****************************************************************
    subroutine compute_outeqm(brsurf, v1, br, bs, bp) 
        real(d), intent(inout), dimension(:,:,:) :: br, bs, bp
        real(d), intent(in), dimension(:,:) :: brsurf
        real(d), intent(in) :: v1
        ! ---
        ! Compute outflow eqm and return vector potential as, ap on
        ! edges. Component ar is zero in this gauge.
        !   brsurf(ns+1,np+1) is the br array on the lower boundary (including ghost cells)
        !   v1 is the speed parameter (in km/s)
        ! ---
        real(d) :: ms(np-1), trigs(np-1,np-1)
        real(d) :: ls(np-1,ns-1), legs(np-1,ns-1,ns-1)
        ! real(d) :: hc(nr+1), g(nr)
        real(d) :: sig(ns), sc(ns-1), sigc(ns-1)
        real(d) :: cml, fs0, fs1, fr0, ff
        integer :: l, m, i, j, k, opp(np+1)
        real(d) :: rc(nr+1), vc(nr+1), dv(nr+1), radc(nr+1), rrad(nr)
        real(d) :: A, B, C
        real(d) :: hc(nr+1,np-1,ns-1), g(nr,np-1,ns-1)
    
        ! Compute azimuthal eigenvalues m and eigenvectors:
        call findms(ms, trigs)

        ! Find latitudinal eigenvalues l and eigenvectors, for each m:
        do m=1,np-1
            call findls(ms(m), ls(m,:), legs(m,:,:))
        end do

        ! Run through l,m modes and compute expansion coefficients 
        ! + radial modes + magnetic field components:
        br = 0.0_d
        bs = 0.0_d
        bp = 0.0_d

        ! Precompute factors:
        sig = sqrt(1.0_d - s**2)
        sc = (/ (s(1) + (dble(i) - 0.5_d)*ds, i=1,ns-1) /)
        sigc = sqrt(1.0_d - sc**2)
        rc = (/ (r(1) + (dble(i) - 0.5_d)*dr, i=0,nr) /)
        radc = exp(rc)
        rrad = exp(2.0_d*r)
        call vout(rc, v1, vc)
        call voutdiff(rc, v1, dv)
        vc = vc*radc
        dv = dv*radc

        ! Compute c_ml*[radial functions] for each m and l:
        hc = 0.0_d
        hc(nr+1,:,:) = -1.0_d  ! seems to be best numerical option
        g = 0.0_d
        g(1,:,:) = 1.0_d   ! lower boundary condition
        !$omp parallel private(l,k,j,i,cml,A,B,C)
        !$omp do
        do m=1,np-1
            do l=1,ns-1
                ! Compute coefficient cml using orthogonality of eigenvectors:
                cml = 0.0_d
                B = 0.0_d  ! denominator
                do k=2,np
                    do j=2,ns
                        A = trigs(k-1,m)*legs(m,j-1,l)
                        cml= cml + A*brsurf(j,k)
                        B = B + A**2
                    end do
                end do
                if (abs(B) < 1d-10) then
                    cml = 0.0_d
                else
                    cml = cml/B
                end if
        
                if (abs(cml) < 1d-10) then
                    cml = 0.0_d
                else
                    ! Find radial function H_l(rho) at cell centres
                    ! - integrate backwards from the outer boundary:
                    do i=nr-1, 1, -1
                        A = 1.0_d/dr**2
                        B = 0.5_d*(3.0_d - vc(i+1))/dr
                        C = 2.0_d - ls(m,l) - 3.0_d*vc(i+1) - dv(i+1)
                        hc(i,m,l) = hc(i+1,m,l) * (2.0_d*A - C)
                        hc(i,m,l) = hc(i,m,l) + hc(i+2,m,l)*(-A - B)
                        hc(i,m,l) = hc(i,m,l) / (A - B)
                    end do
                    ! - normalize to satisfy lower boundary condition:
                    hc(:,m,l) = hc(:,m,l)*dr/(hc(2,m,l)*radc(2) - hc(1,m,l)*radc(1))

                    ! Find G_l(rho) from H_l(rho) using recurrence.
                    do i=2,nr
                        g(i,m,l) = ( 0.5_d*ls(m,l)*hc(i,m,l)*(rrad(i) - rrad(i-1)) &
                        + g(i-1,m,l)*rrad(i-1) ) / rrad(i)
                    end do

                    hc(:,m,l) = hc(:,m,l)*cml  ! for later convenience
                    g(:,m,l) = g(:,m,l)*cml 
                end if
            end do
        end do
        !$omp end do
        !$omp end parallel

        ! Fill up magnetic field components for each m:
        !$omp parallel private(j,i,m)
        !$omp do
        do k=2,np
            do j=2,ns
                do i=1,nr
                    do m=1,np-1
                        br(i,j,k) = br(i,j,k) + sum(g(i,m,:)*legs(m,j-1,:))*trigs(k-1,m)
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel  
        !$omp parallel private(j,i,m)
        !$omp do
        do k=2,np
            do j=2,ns-1
                do i=2,nr
                    do m=1,np-1
                        bs(i,j,k) = bs(i,j,k) + sig(j)*sum(hc(i,m,:)*(legs(m,j,:) - legs(m,j-1,:)))/ds*trigs(k-1,m)
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel  
        !$omp parallel private(k,i,m)
        !$omp do       
        do j=2,ns
            do k=2,np-1
                do i=2,nr
                    do m=1,np-1
                        bp(i,j,k) = bp(i,j,k) + sum(hc(i,m,:)*legs(m,j-1,:))/sigc(j-1)*(trigs(k,m) - trigs(k-1,m))/dp
                    end do
                end do
            end do
            do i=2,nr
                do m=1,np-1
                    bp(i,j,1) = bp(i,j,1) + sum(hc(i,m,:)*legs(m,j-1,:))/sigc(j-1)*(trigs(1,m) - trigs(np-1,m))/dp
                end do
                bp(i,j,np) = bp(i,j,1)
            end do
        end do
        !$omp end do
        !$omp end parallel         
        
        ! Fill magnetic field boundary values:
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
        !$omp parallel private(ff)
        !$omp do
        do j=1,ns
           ff = asin(s(j)+0.5_d*ds) - asin(s(j)-0.5_d*ds)
           bs(1,j,:) = fs0/ff*(bs(2,j,:)*fs1*ff &
                + (br(1,j,:) - br(1,j+1,:))*fr0)
        end do
        !$omp end do
        !$omp end parallel 
        fs0 = fs0/dp
        fs1 = fs1*dp
        fr0 = exp(r(1))*(exp(0.5_d*dr) - exp(-0.5_d*dr))
        !$omp parallel private(ff, i)
        !$omp do
        do j=2,ns
           ff = sqrt(1.0_d - (s(j)-0.5_d*ds)**2)
           do i=1,np
              bp(1,j,i) = fs0/ff*(bp(2,j,i)*fs1*ff &
                   + (br(1,j,i) - br(1,j,i+1))*fr0)
           end do
        end do
        !$omp end do
        !$omp end parallel 
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

    end subroutine compute_outeqm

    !****************************************************************
    subroutine findms(ms, vs) 
        real(d), intent(inout), dimension(:) :: ms
        real(d), intent(inout), dimension(:,:) :: vs
        ! ---
        ! Calculate azimuthal eigenvalues ms(1:np-1) and 
        ! eigenvectors vs(1:np-1,1:np-1).
        !
        ! (ms would be integers in infinite limit but not necessarily
        ! here.)
        !
        ! We have to combine eigenvectors from both the cosines and 
        ! sines, hence this is more complicated than the latitudinal
        ! direction.
        ! ---
        real(d) :: dvals(np-1), evals(np-1), v1(np-1,np-1)
        integer :: i, j

        vs = 0.0_d
        
        ! (1) Calculate sine eigenvalues:
        ! - diagonal elements
        dvals = 2.0_d
        dvals(1) = 1.0_d
        dvals(np-1) = 1.0_d
        ! - off-diagonal elements (first one is ignored):
        evals = -1.0_d
        ! - initialize v1 to identity matrix:
        v1 = 0.0_d
        do j=1,np-1
            v1(j,j) = 1.0_d
        end do
        call trieig(dvals, evals, v1)  ! dvals will be eigenvalues
        ! Extract only periodic ones (nb. eigenvalues not ordered):
        ms = 0.0_d
        i = 1
        do j=1,np-1
            if (v1(1,j)*v1(np-1,j) > 0) then
                ms(i) = sqrt(abs(dvals(j))/dp**2)
                vs(:,i) = v1(:,j)
                i = i+1
            end if
        end do

        ! (2) Calculate cosine eigenvalues:
        ! - diagonal elements
        dvals = 2.0_d
        dvals(1) = 3.0_d
        dvals(np-1) = 3.0_d
        ! - off-diagonal elements (first one is ignored):
        evals = -1.0_d
        ! - initialize v1 to identity matrix:
        v1 = 0.0_d
        do j=1,np-1
            v1(j,j) = 1.0_d
        end do
        call trieig(dvals, evals, v1)  ! dvals will be eigenvalues
        ! Remove repeated modes (nb. eigenvalues not ordered):
        do j=1,np-1
            if (v1(2,j)*v1(np-2,j) < 0) then
                ms(i) = sqrt(abs(dvals(j))/dp**2)
                vs(:,i) = v1(:,j)
                i = i+1
            end if
        end do

    end subroutine findms

    !****************************************************************
    subroutine findls(m, ls, vs) 
        real(d), intent(in) :: m
        real(d), intent(inout), dimension(:) :: ls
        real(d), intent(inout), dimension(:,:) :: vs
        ! ---
        ! Calculate latitudinal eigenvalues ls(1:ns-1) and 
        ! eigenvectors vs(1:ns-1,1:ns-1) for given m.
        ! ---
        integer :: i, j
        real(d) :: sigc(ns-1), sigs(ns)
        real(d) :: evals(ns-1)

        sigc = (/ (s(1) + (dble(i) + 0.5_d)*ds, i=0,ns-2) /)
        sigc = sqrt(1.0_d - sigc**2)
        sigs = 1.0_d - s**2

        ! Initialize ls to diagonal entries:
        ls = sigs(2:ns) + sigs(1:ns-1) - (ds*m**2/sigc)*(asin(s(1:ns-1)) - asin(s(2:ns)))

        ! Off-diagonal entries:
        evals(2:ns-1) = -sigs(2:ns-1)

        ! Initialize vs to identity matrix:
        vs = 0.0_d
        do j=1,ns-1
            vs(j,j) = 1.0_d
        end do

        ! Compute eigenvalues and eigenvectors:
        call trieig(ls, evals, vs)
        ls = ls/ds**2

    end subroutine findls

    !****************************************************************
    subroutine vout(rs, v1, vr)
        real(d), intent(in), dimension(:) :: rs
        real(d), intent(in) :: v1
        real(d), intent(inout), dimension(:) :: vr
        real(d), parameter :: RCRIT = 10.0_d
        ! ---
        ! Specify v_r function for outflow.
        ! This is (an approximation to) Parker's solar wind solution,
        ! with critical radius RCRIT.
        ! Parameter v1 is the overall magnitude.
        ! (Could be changed to any reasonable equivalent.)
        ! ---

        vr = exp(-2.0_d*RCRIT/exp(rs)) * exp(2.0_d*r(nr))
        vr = vr / exp(-2.0_d*RCRIT/exp(r(nr))) / exp(2.0_d*rs)
        vr = vr * v1

    end subroutine vout

    !****************************************************************
    subroutine voutdiff(rs, v1, dvr)
        real(d), intent(in), dimension(:) :: rs
        real(d), intent(in) :: v1
        real(d), intent(inout), dimension(:) :: dvr
        real(d), parameter :: RCRIT = 10.0_d
        ! ---
        ! Specify v_r'' function for outflow.
        ! This is (an approximation to) Parker's solar wind solution, 
        ! with critical radius RCRIT.
        ! Parameter v1 is the overall magnitude.
        ! (Could be changed to any reasonable equivalent.)
        ! ---

        dvr = v1*exp(2.0_d*r(nr)) * exp(2.0_d*RCRIT/exp(r(nr)))
        dvr = dvr*exp(rs)*2.0_d*exp(-2.0_d*RCRIT/exp(rs))*(RCRIT-exp(rs))*exp(-4.0_d*rs)

    end subroutine voutdiff

    !****************************************************************
    subroutine b_to_netcdf(br, bs, bp, wkpath, snap)
        real(d), intent(inout), dimension(:,:,:) :: br, bs, bp
        character*(*), intent(in) :: wkpath, snap
        integer, parameter:: nfreal=nf90_double
        ! ---
        ! Output br, bs, bp on cell faces to netcdf file.
        ! ---     
        real(d) :: rc(nr+1), sc(ns+1), pc(np+1), thc(ns+1)
        integer :: i
        integer :: ncid, br_vid, bs_vid, bp_vid, r_vid, s_vid, p_vid
        integer :: rc_vid, sc_vid, pc_vid
        integer :: nr_did, nrc_did, ns_did, nsc_did, np_did, npc_did
        
        ! Output to netcdf file:
        if (VERBOSE) print*,'Writing '//wkpath//snap
        call chk(nf90_create(wkpath//snap, nf90_clobber, ncid))
        call chk(nf90_def_dim(ncid,"r",nr,nr_did))
        call chk(nf90_def_dim(ncid,"th",ns,ns_did))
        call chk(nf90_def_dim(ncid,"ph",np,np_did))
        call chk(nf90_def_dim(ncid,"rc",nr+1,nrc_did))
        call chk(nf90_def_dim(ncid,"thc",ns+1,nsc_did))
        call chk(nf90_def_dim(ncid,"phc",np+1,npc_did))

        call chk(nf90_def_var(ncid,"r",nfreal,nr_did,r_vid))
        call chk(nf90_def_var(ncid,"th",nfreal,ns_did,s_vid))
        call chk(nf90_def_var(ncid,"ph",nfreal,np_did,p_vid))
        call chk(nf90_def_var(ncid,"rc",nfreal,nrc_did,rc_vid))
        call chk(nf90_def_var(ncid,"thc",nfreal,nsc_did,sc_vid))
        call chk(nf90_def_var(ncid,"phc",nfreal,npc_did,pc_vid))        
        call chk(nf90_def_var(ncid,"br",nfreal,(/nr_did,nsc_did,npc_did/),br_vid))
        call chk(nf90_def_var(ncid,"bth",nfreal,(/nrc_did,ns_did,npc_did/),bs_vid))
        call chk(nf90_def_var(ncid,"bph",nfreal,(/nrc_did,nsc_did,np_did/),bp_vid))

        call chk(nf90_enddef(ncid))

        call chk(nf90_put_var(ncid,r_vid,exp(r)))
        call chk(nf90_put_var(ncid,s_vid,acos(s)))
        call chk(nf90_put_var(ncid,p_vid,p))
        
        rc = (/ (r(1) + (dble(i) - 0.5_d)*dr, i=0,nr) /)
        sc = (/ (s(1) + (dble(i) - 0.5_d)*ds, i=0,ns) /)
        pc = (/ (p(1) + (dble(i) - 0.5_d)*dp, i=0,np) /)

        thc = acos(sc)
        thc(1) = 2.0_d*thc(2) - thc(3)
        thc(ns+1) = 2.0_d*thc(ns) - thc(ns-1)

        call chk(nf90_put_var(ncid,rc_vid,exp(rc)))
        call chk(nf90_put_var(ncid,sc_vid,thc))
        call chk(nf90_put_var(ncid,pc_vid,pc))        
        
        call chk(nf90_put_var(ncid,br_vid,br))
        call chk(nf90_put_var(ncid,bs_vid,-bs))
        call chk(nf90_put_var(ncid,bp_vid,bp))

        call chk( nf90_close(ncid) )
    
    end subroutine b_to_netcdf
    
end module outeqm
