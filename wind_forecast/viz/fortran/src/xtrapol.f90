!****************************************************************
module xtrapol
!****************************************************************

    use shared
    use readnc
    use ffts
    use netcdf
    
    implicit none

    private
    public :: compute_pfss_inner, compute_pfrf_inner, compute_schatten
    
    contains
  
    !****************************************************************
    subroutine compute_pfss_inner(wkpath, snap) 
        character*(*), intent(in) :: wkpath, snap
        ! ---
        ! Compute PFSS extrapolation in original domain (with source surface at top).
        ! Save B on cell faces to new netcdf file.
        ! ---       
        real(d) :: ar(1,1,1), as(nr,ns-1,np), ap(nr,ns,np-1)
          
        if (VERBOSE) print*,'Computing PFSS...'
        call compute_pf(br0, as, ap)
        ar = 0.0_d
        call a_to_bnetcdf(ar, as, ap, wkpath, 'pfss_'//snap)
          
    end subroutine compute_pfss_inner
    
    !****************************************************************
    subroutine compute_pfrf_inner(wkpath, snap) 
        character*(*), intent(in) :: wkpath, snap
        ! ---
        ! Compute reference potential field in original domain (matching br top and bottom).
        ! Save B on cell faces to new netcdf file.
        ! ---       
        real(d) :: ar(1,1,1), as(nr,ns-1,np), ap(nr,ns,np-1)
          
        if (VERBOSE) print*,'Computing reference potential field...'
        call compute_pf(br0, as, ap, br1)
        ar = 0.0_d
        call a_to_bnetcdf(ar, as, ap, wkpath, 'pfrf_'//snap)
          
    end subroutine compute_pfrf_inner    
          
    !****************************************************************
    subroutine compute_schatten(wkpath, snap, rad2)
        character*(*), intent(in) :: wkpath, snap
        real(d), intent(in) :: rad2
        integer, parameter :: n2=40
        ! ---
        ! Compute (potential) Schatten extension to dumfric datacube.
        ! Input rad2 is the outer source-surface radius (heliospheric boundary).
        ! Parameter:
        !   n2   is the number of grid points in the extension field
        ! Save B on cell faces to new netcdf file [note: output field points
        ! everywhere outward]
        ! ---
        real(d) :: ar(1,1,1), as(n2,ns-1,np), ap(n2,ns,np-1), r1, br00
        integer :: i

        if (VERBOSE) print*,'Computing Schatten extension...'
        
        ! Change r array to correspond to outer region:
        r1 = r(nr)
        deallocate(r)
        allocate(r(n2))
        dr = (log(rad2) - r1)/dble(n2-1)
        r = (/ (r1 + dble(i)*dr, i=0,n2-1) /)
        nr = n2
        
        ! Remove monopole component:
        br00 = sum(abs(br1(2:ns,2:np)))/dble((ns-1)*(np-1))
        
        ! Compute PFSS in outer region:
        call compute_pf(abs(br1)-br00, as, ap)
        ar = 0.0_d 
        
        ! Add back monopole component (do this to Br directly):
        call a_to_bnetcdf(ar, as, ap, wkpath, 'schat_'//snap, monopole=br00)
        
    end subroutine compute_schatten  
    
    
!-----------------------------------------------------------------------

   !****************************************************************
    subroutine compute_pf(brsurf, as, ap, brouter)
        real(d), intent(inout), dimension(:,:,:) :: as, ap
        real(d), intent(in), dimension(:,:) :: brsurf
        real(d), intent(in), dimension(:,:), optional :: brouter
        ! ---
        ! Compute potential field and return vector potential as, ap on
        ! edges. Component ar is zero in this gauge.
        !   brsurf(ns+1,np+1) is the br array on the lower boundary (including ghost cells)
        !
        ! If optional argument brouter(ns+1,np+1) is present than br is set to match this
        ! on the outer boundary, otherwise a PFSS extrapolation is computed.
        ! ---
        real(d) :: e0, e1, ed1, fact, dnr, ratio
        complex(d) :: cdlm0, cdlm1
        integer :: i, j, m
        real(d) :: Vg(ns), sc(ns-1), Uc(ns-1)
        complex(d) :: brt0(ns-1,np-1), psit(nr,ns-1,np-1)
        complex(d), dimension(:,:), allocatable :: brt1
        real(d) :: psi(nr,ns-1,np-1), mu(np-1), Q(ns-1,ns-1)
        real(d) :: lam(ns-1), ev(ns-1), fp(ns-1), fm(ns-1), k(nr)

        ! Precompute:
        e0 = dexp(2.0_d*r(1))
        e1 = dexp(2.0_d*r(nr))
        ed1 = dexp(dr)
        fact = (ed1 - 1.0_d)*(dexp(0.5_d*dr) - dexp(-0.5*dr))
        k = (/ (dble(i), i=0,nr-1) /)
        dnr = dble(nr-1)

        ! Initialise:
        ! - arrays of coordinate factors
        sc = (/ (s(1) + (dble(i) + 0.5_d)*ds, i=0,ns-2) /)
        Vg = 0.0_d
        Vg(2:ns-1) = dsqrt(1.0_d - s(2:ns-1)**2)/ &
            (dasin(sc(2:ns-1)) - dasin(sc(1:ns-2)))/ds
        Uc = (dasin(s(2:ns)) - dasin(s(1:ns-1)))/dsqrt(1.0_d - sc**2)/ds/dp**2

        ! FFT of br0 in p direction
        brt0 = cmplx(brsurf(2:ns,2:np), kind=d)
        call fftn(brt0, shape(brt0), dim=(/2/))
        brt0 = brt0/dsqrt(dble(np-1))

        ! FFT of br1 in p direction [if required]
        if (present(brouter)) then
            allocate(brt1(ns-1,np-1))
            brt1 = cmplx(brouter(2:ns,2:np), kind=d)
            call fftn(brt1, shape(brt1), dim=(/2/))
            brt1 = brt1/dsqrt(dble(np-1))    
        end if
        
        ! Term required for m-dependent part of matrix
        mu(1:np-1) = (/ (dble(i), i=0,np-2) /)/dble(np-1)
        !
        ! - order frequencies [+,-] for FFT
        mu((np-1)/2+1:np-1) = mu((np-1)/2+1:np-1) - 1.0_d
        mu = 4.0_d*dsin(PI*mu)**2

        ! Loop over azimuthal modes (positive m)
        do m=1,(np-1)/2+1
            ! - prepare tridiagonal matrix (ev = off-diagonal and lam diagonal).
            ! - also prepare identity matrix Q
            ev = -Vg(1:ns-1)
            Q = 0.0_d
            do i=1,ns-1
                lam(i) = Vg(i) + Vg(i+1) + Uc(i)*mu(m)
                Q(i,i) = 1.0_d
            end do
            !
            ! - compute eigenvectors Q_{lm} and eigenvalues lam_{lm}
            call trieig(lam, ev, Q)
            !
            ! - solve quadratic
            fp = 0.5_d*(1.0_d + ed1 + lam*fact)
            fp = fp + dsqrt(fp*fp - ed1)
            fm = ed1/fp
            !
            ! - compute radial term for each l (for this m):
            !$omp parallel private(cdlm0, cdlm1, ratio)
            !$omp do
            do i=1,ns-1
                ! - sum c_{lm} + d_{lm} from boundary condition at photosphere
                cdlm0 = dot_product(Q(:,i), brt0(:,m))*e0/lam(i)
                ! - solve simultaneously using outer boundary condition
                if (present(brouter)) then
                    ! - sum c_{lm}*ffp**nr + d_{lm}*ffm**nr
                    cdlm1 = dot_product(Q(:,i), brt1(:,m))*e1/lam(i)
                    ! - add contribution from c_{lm}:
                    psit(:,i,m) = (cdlm1 - fm(i)**dnr*cdlm0)*fp(i)**k/(fp(i)**dnr - fm(i)**dnr)
                    ! - add contribution from d_{lm}:
                    psit(:,i,m) = psit(:,i,m) + &
                            (cdlm1 - fp(i)**dnr*cdlm0)*fm(i)**k/(fm(i)**dnr - fp(i)**dnr)
                else
                    ! - ratio d_{lm}/c_{lm}
                    ratio = (fm(i)**(dnr-1.0_d) - fm(i)**dnr) &
                        /(fp(i)**dnr - fp(i)**(dnr-1.0_d))
                    psit(:,i,m) = cdlm0*(ratio*fp(i)**k + fm(i)**k)/(1.0_d + ratio)
                end if
            end do
            !$omp end do
            !$omp end parallel
            !
            ! - compute entry for this m in DFT of psi:
            psit(:,:,m) = matmul(psit(:,:,m),transpose(Q))
            if (m > 1) psit(:,:,np+1-m) = conjg(psit(:,:,m))
        end do

        ! Compute psi by inverse FFT
        call fftn(psit, shape(psit), dim=(/3/), inv=.true.)
        psi = dble(psit)*dsqrt(dble(np-1))

        ! Hence compute as and ap
        as = 0.0_d
        ap = 0.0_d
        !$omp parallel private(i)
        !$omp do
        do j=1,nr
            do i=1,np
                as(j,:,i) = (psi(j,:,mod(i-2+np-1,np-1)+1) - psi(j,:,mod(i-1+np-1,np-1)+1)) &
                    /dexp(r(j))/dsqrt(1.0_d - sc**2)/dp
            end do
            do i=1,np-1
                ap(j,2:ns-1,i) = (psi(j,2:ns-1,i) - psi(j,1:ns-2,i)) &
                    /dexp(r(j))/(dasin(sc(2:ns-1)) - dasin(sc(1:ns-2)))
            end do
        end do
        !$omp end do
        !$omp end parallel       
            
    end subroutine compute_pf

    !****************************************************************
    subroutine a_to_bnetcdf(ar, as, ap, wkpath, snap, monopole)
        real(d), intent(inout), dimension(:,:,:) :: ar, as, ap
        character*(*), intent(in) :: wkpath, snap
        real(d), intent(in), optional :: monopole
        integer, parameter:: nfreal=nf90_double
        ! ---
        ! Compute br, bs, bp on cell faces from ar, as, ap
        ! and output to netcdf file.
        ! [Note: if ar has only 1 element, then it is assumed zero and ignored.]
        ! Optional monopole component is added to br.
        ! ---     
        real(d) :: br(nr,ns+1,np+1), bs(nr+1,ns,np+1), bp(nr+1,ns+1,np)
        real(d) :: rc(nr+1), sc(ns+1), pc(np+1), thc(ns+1)
        integer :: i
        integer :: ncid, br_vid, bs_vid, bp_vid, r_vid, s_vid, p_vid
        integer :: rc_vid, sc_vid, pc_vid
        integer :: nr_did, nrc_did, ns_did, nsc_did, np_did, npc_did
        
        call computeBfromA(ar, as, ap, br, bs, bp)
         
        ! Add monopole component (used for Schatten solution):
        if (present(monopole)) then
            do i=1,nr
                br(i,:,:) = br(i,:,:) + monopole*(exp(r(1) - r(i)))**2
            end do
        end if
        
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
    
    end subroutine a_to_bnetcdf
    
end module xtrapol
