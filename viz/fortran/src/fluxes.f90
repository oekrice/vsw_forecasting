!****************************************************************
module fluxes
!****************************************************************

    use shared
    use readnc

!----------------------------------------------------------------
! Read netcdf files and set up interpolators.
!----------------------------------------------------------------
    implicit none
    
    private
    public :: compute_hfluxes

    contains

    !****************************************************************
    subroutine compute_hfluxes(wkpath, filename, filename_nxt, dhours)
        !----------------------------------------------------------------
        ! Compute terms in helicity evolution equation, using two snapshots
        ! to estimate time derivatives.
        !----------------------------------------------------------------
        character*(*), intent(in) :: wkpath, filename, filename_nxt
        integer, intent(in) :: dhours
        real(d), dimension(:,:,:), allocatable :: arg_nxt, asg_nxt, apg_nxt
        real(d), dimension(:,:,:), allocatable :: brg_nxt, bsg_nxt, bpg_nxt
        real(d), dimension(:,:,:), allocatable :: erg, esg, epg
        real(d) :: hr, hr_nxt, f_dhdt, f_diss, f_ae0, f_ae1, f_ada0, f_ada1
        real(d) :: dsecs
        integer :: nr,ns,np

        ! Read second snapshot and store A at grid points:
        call readB(trim(wkpath)//trim(filename_nxt), .true., .true.)

        nr = size(r)
        ns = size(s)
        np = size(p)
        allocate(arg_nxt(nr,ns,np),asg_nxt(nr,ns,np),apg_nxt(nr,ns,np))
        arg_nxt = arg
        asg_nxt = asg
        apg_nxt = apg
        brg_nxt = brg
        bsg_nxt = bsg
        bpg_nxt = bpg

        ! Read first snapshot and store both A and B at grid points:
        call readB(trim(wkpath)//trim(filename), .true., .true.)
        
        ! Read electric field from first snapshot:
        allocate(erg(nr,ns,np),esg(nr,ns,np),epg(nr,ns,np))
        call readvec_gpts(trim(wkpath)//trim(filename), 'e', erg, esg, epg)

        ! Compute terms in evolution equation:
        dsecs = dble(dhours)*3600.0_d

        ! (a) d/dt int_V(A.B)dV
        hr = trap3d(arg*brg + asg*bsg + apg*bpg)
        hr_nxt = trap3d(arg_nxt*brg_nxt + asg_nxt*bsg_nxt + apg_nxt*bpg_nxt)
        f_dhdt = (hr_nxt - hr)/dsecs

        ! (b) -2*int_V(E.B)dV
        f_diss = trap3d(erg*brg + esg*bsg + epg*bpg)
        f_diss = -2*f_diss

        ! (c) oint_{S_0} A x [2*E] dS
        f_ae0 = trap2d(apg(1,:,:)*esg(1,:,:) - asg(1,:,:)*epg(1,:,:))
        f_ae0 = 2.0_d*f_ae0*exp(2.0_d*r(1))

        ! (d) oint_{S_1} A x [2*E] dS
        f_ae1 = trap2d(apg(nr,:,:)*esg(nr,:,:) - asg(nr,:,:)*epg(nr,:,:))
        f_ae1 = 2.0_d*f_ae1*exp(2.0_d*r(nr))

        ! (e) oint_{S_0} A x dA/dt dS
        f_ada0 = trap2d(apg(1,:,:)*(asg_nxt(1,:,:) - asg(1,:,:)) &
            - asg(1,:,:)*(apg_nxt(1,:,:) - apg(1,:,:)))
        f_ada0 = f_ada0*exp(2.0_d*r(1))/dsecs

        ! (f) oint_{S_1} A x dA/dt dS
        f_ada1 = trap2d(apg(nr,:,:)*(asg_nxt(nr,:,:) - asg(nr,:,:)) &
            - asg(nr,:,:)*(apg_nxt(nr,:,:) - apg(nr,:,:)))
        f_ada1 = f_ada1*exp(2.0_d*r(nr))/dsecs

        ! Output:
        open(unit=1, file=trim(wkpath)//'hflux_'//trim(filename)//'.unf', form='unformatted')
        write(1) hr*RSUN**4
        write(1) f_dhdt*RSUN**4
        write(1) f_diss*RSUN**4
        write(1) f_ae0*RSUN**4
        write(1) f_ae1*RSUN**4
        write(1) f_ada0*RSUN**4
        write(1) f_ada1*RSUN**4
        close(1)

    end subroutine compute_hfluxes



end module fluxes
