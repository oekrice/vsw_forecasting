!*******************************************************************************
MODULE hyper
!*******************************************************************************
! This module implements 4th-order hyperdiffusion.
!*******************************************************************************
    USE params
    USE mpitools, ONLY: control, north, south, surface, outer, mpireal, &
        MPI_MIN, comm, stat, err, mpi_r_prv, mpi_r_nxt, mpi_s_prv, mpi_s_nxt, &
        mpi_p_prv, mpi_p_nxt, mpi_opp, i0s, i1s
    USE grid, ONLY: nr, ns, np, dlr, dls, dlp, dnr, dns, Sjr, Sjs, Sjp, Vg
    USE friction, ONLY: nu0

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: hyper_ratio, hyper_flag, hyper_eta4, dt_max_hyper, dhg
    PUBLIC:: hr, hs, hp, alg
    PUBLIC:: initialize_hyper, compute_alpha, compute_hyper_div

!*******************************************************************************
    ! Specified in user namelist:
    REAL(d):: hyper_ratio    ! ratio eta4/(R_sun**4 * nu0)  [dimensionless]

    ! Derived:
    LOGICAL:: hyper_flag     ! flag for whether hyperdiffusion is turned on
    REAL(d):: hyper_eta4    ! coefficient of hyperdiffusivity [R_sun**4/s]
    REAL(d), ALLOCATABLE:: alg(:,:,:)   ! alpha = J.B/B**2 at grid pts
                                        ! (0:nr,0:ns,0:np)
    REAL(d), ALLOCATABLE:: hr(:,:,:)    ! eta4*B**2*grad_rho(alpha) on staggered
                                        ! faces (0:nr+1,0:ns,0:np)
    REAL(d), ALLOCATABLE:: hs(:,:,:)    ! eta4*B**2*grad_s(alpha) on staggered
                                        ! faces (0:nr,0:ns+1,0:np)
    REAL(d), ALLOCATABLE:: hp(:,:,:)    ! eta4*B**2*grad_phi(alpha) on staggered
                                        ! faces (0:nr,0:ns,0:np+1)
    REAL(d), ALLOCATABLE:: dhg(:,:,:)   ! div(eta4*B**2*grad(alg))/B**2 at grid
                                        ! pts (0:nr,0:ns,0:np)
    REAL(d):: dt_max_hyper=1e10_d  ! maximum timestep for hyperdiffusion [s]
    INTEGER, PARAMETER:: tag=1

!*******************************************************************************
CONTAINS
!*******************************************************************************

!===============================================================================
SUBROUTINE initialize_hyper()
! Set hyperdiffusion coefficient.

    ! Set hyperdiffusion coefficient:
    IF (hyper_ratio < 1.0d-20) THEN
       hyper_flag = .FALSE.
       hyper_eta4 = 0.0_d
       dt_max_hyper = 1e10_d
    ELSE
       hyper_flag = .TRUE.
       hyper_eta4 = hyper_ratio * nu0
    END IF

END SUBROUTINE initialize_hyper

!===============================================================================
SUBROUTINE compute_alpha(bgr, bgs, bgp, bbgs, jgr, jgs, jgp)
! Compute alg at grid points.
    REAL(d), INTENT(IN):: bgr(0:nr,0:ns,0:np)   ! B_rho at grid pts
    REAL(d), INTENT(IN):: bgs(0:nr,0:ns,0:np)   ! B_s at grid pts
    REAL(d), INTENT(IN):: bgp(0:nr,0:ns,0:np)   ! B_phi at grid pts
    REAL(d), INTENT(IN):: bbgs(-1:nr+1,0:ns,0:np)   ! Softened B**2 at grid pts
    REAL(d), INTENT(IN):: jgr(0:nr,0:ns,0:np)   ! J_rho at grid pts
    REAL(d), INTENT(IN):: jgs(0:nr,0:ns,0:np)   ! J_s at grid pts
    REAL(d), INTENT(IN):: jgp(0:nr,0:ns,0:np)   ! J_phi at grid pts

    IF (.NOT. ALLOCATED(alg)) ALLOCATE(alg(0:nr,0:ns,0:np))

    alg(0:nr,0:ns,0:np) = (bgr(0:nr,0:ns,0:np) * jgr(0:nr,0:ns,0:np) &
        + bgs(0:nr,0:ns,0:np) * jgs(0:nr,0:ns,0:np) &
        + bgp(0:nr,0:ns,0:np) * jgp(0:nr,0:ns,0:np) ) &
        / bbgs(0:nr,0:ns,0:np)

END SUBROUTINE compute_alpha

!===============================================================================
SUBROUTINE compute_hyper_div(bbgs)
! Compute div( eta4 * bbg * grad(alg) )/bbg as required for electric field.
    REAL(d), INTENT(IN):: bbgs(-1:nr+1,0:ns,0:np)   ! Softened B**2 at grid pts

    ! Local variables:
    REAL(d):: local_dt, dt_r, dt_s, dt_p

    IF (.NOT. ALLOCATED(hr)) ALLOCATE(hr(0:nr+1,0:ns,0:np))
    IF (.NOT. ALLOCATED(hs)) ALLOCATE(hs(0:nr,0:ns+1,0:np))
    IF (.NOT. ALLOCATED(hp)) ALLOCATE(hp(0:nr,0:ns,0:np+1))
    IF (.NOT. ALLOCATED(dhg)) ALLOCATE(dhg(0:nr,0:ns,0:np))

    ! Compute hr, hs, hp on interior ribs [eta4*B^2*grad(alpha)]
    ! - note: eta4 and 0.5 are included later
    hr(1:nr,0:ns,0:np) = (bbgs(0:nr-1,0:ns,0:np) + bbgs(1:nr,0:ns,0:np)) &
       * (alg(1:nr,0:ns,0:np) - alg(0:nr-1,0:ns,0:np)) / dlr(1:nr,0:ns,0:np)
    hs(0:nr,1:ns,0:np) = (bbgs(0:nr,0:ns-1,0:np) + bbgs(0:nr,1:ns,0:np)) &
       * (alg(0:nr,1:ns,0:np) - alg(0:nr,0:ns-1,0:np)) / dls(0:nr,1:ns,0:np)
    hp(0:nr,1:ns-1,1:np) = (bbgs(0:nr,1:ns-1,0:np-1) + bbgs(0:nr,1:ns-1,1:np)) &
       * (alg(0:nr,1:ns-1,1:np) - alg(0:nr,1:ns-1,0:np-1)) &
       / dlp(0:nr,1:ns-1,1:np)
    IF (.NOT. south) hp(0:nr,0,1:np) = &
         (bbgs(0:nr,0,0:np-1) + bbgs(0:nr,0,1:np)) &
       * (alg(0:nr,0,1:np) - alg(0:nr,0,0:np-1)) &
       / dlp(0:nr,0,1:np)
    IF (.NOT. north) hp(0:nr,ns,1:np) = &
         (bbgs(0:nr,ns,0:np-1) + bbgs(0:nr,ns,1:np)) &
       * (alg(0:nr,ns,1:np) - alg(0:nr,ns,0:np-1)) &
       / dlp(0:nr,ns,1:np)

    ! MPI transfer hp in p direction
    CALL MPI_SENDRECV(hp(0:nr,0:ns,np), (nr+1) * (ns+1), mpireal, mpi_p_nxt, &
        tag, hp(0:nr,0:ns,0),(nr+1) * (ns+1), mpireal, mpi_p_prv, tag, comm, &
        stat, err)
    CALL MPI_SENDRECV(hp(0:nr,0:ns,1), (nr+1) * (ns+1), mpireal, mpi_p_prv, &
        tag, hp(0:nr,0:ns,np+1), (nr+1) * (ns+1), mpireal, mpi_p_nxt, tag, &
        comm, stat, err)
    
    ! MPI transfer hs at interior boundaries in s direction
    ! - up
    IF (.NOT. north) CALL MPI_SEND(hs(0:nr,ns,0:np), (nr+1) * (np+1), mpireal, &
        mpi_s_nxt, tag, comm, err)
    IF (.NOT. south) CALL MPI_RECV(hs(0:nr,0,0:np), (nr+1) * (np+1), mpireal, &
        mpi_s_prv, tag, comm, stat, err)
    ! - down
    IF (.NOT. south) CALL MPI_SEND(hs(0:nr,1,0:np), (nr+1) * (np+1), mpireal, &
        mpi_s_prv, tag, comm, err)
    IF (.NOT. north) CALL MPI_RECV(hs(0:nr,ns+1,0:np), (nr+1) * (np+1), &
        mpireal, mpi_s_nxt, tag, comm, stat, err)
    
    ! MPI transfer at polar boundaries in s
    IF (south) THEN
        CALL MPI_SENDRECV(hs(0:nr,1,0:np), (nr+1) * (np+1), mpireal, mpi_opp, &
            tag, hs(0:nr,0,0:np), (nr+1) * (np+1), mpireal, mpi_opp, tag, &
            comm, stat, err)
        CALL MPI_SENDRECV(hp(0:nr,1,0:np+1), (nr+1) * (np+2), mpireal, mpi_opp, &
            tag, hp(0:nr,0,0:np+1), (nr+1) * (np+2), mpireal, mpi_opp, tag, &
            comm, stat, err)
        hp(0:nr,0,0:np+1) = 0.5_d * (hp(0:nr,1,0:np+1) - hp(0:nr,0,0:np+1))
    END IF
    IF (north) THEN
        CALL MPI_SENDRECV(hs(0:nr,ns,0:np), (nr+1) * (np+1), mpireal, mpi_opp, &
            tag, hs(0:nr,ns+1,0:np), (nr+1) * (np+1), mpireal, mpi_opp, tag, &
            comm, stat, err)
        CALL MPI_SENDRECV(hp(0:nr,ns-1,0:np+1), (nr+1) * (np+2), mpireal, &
            mpi_opp, tag, hp(0:nr,ns,0:np+1), (nr+1) * (np+2), mpireal, &
            mpi_opp, tag, comm, stat, err)
        hp(0:nr,ns,0:np+1) = 0.5_d * (hp(0:nr,ns-1,0:np+1) - hp(0:nr,ns,0:np+1))
    END IF
    
    ! MPI transfer hr interior boundaries in r direction
    ! - up
    IF (.NOT. outer) CALL MPI_SEND(hr(nr,0:ns,0:np), (ns+1) * (np+1), &
        mpireal, mpi_r_nxt, tag, comm, err)
    IF (.NOT. surface) CALL MPI_RECV(hr(0,0:ns,0:np), (ns+1) * (np+1), &
        mpireal, mpi_r_prv, tag, comm, stat, err)
    ! - down
    IF (.NOT. surface) CALL MPI_SEND(hr(1,0:ns,0:np), (ns+1) * (np+1), &
        mpireal, mpi_r_prv, tag, comm, err)
    IF (.NOT. outer) CALL MPI_RECV(hr(nr+1,0:ns,0:np), (ns+1) * (np+1), &
        mpireal, mpi_r_nxt, tag, comm, stat, err)
    
    ! Apply boundary conditions to hr in r [zero gradient]
    IF (surface) hr(0,0:ns,0:np) = hr(1,0:ns,0:np)
    IF (outer) hr(nr+1,0:ns,0:np) = hr(nr,0:ns,0:np)
    
    ! Compute div(hr,hs,hp) at grid points
    ! [note: Sjs includes opp sign for ghost point hs]
    hr(0:nr+1,0:ns,0:np) = hr(0:nr+1,0:ns,0:np) * Sjr(0:nr+1,0:ns,0:np)
    hs(0:nr,0:ns+1,0:np) = hs(0:nr,0:ns+1,0:np) * Sjs(0:nr,0:ns+1,0:np)
    hp(0:nr,0:ns,0:np+1) = hp(0:nr,0:ns,0:np+1) * Sjp(0:nr,0:ns,0:np+1)
    dhg(0:nr,0:ns,0:np) = (hr(1:nr+1,0:ns,0:np) - hr(0:nr,0:ns,0:np) &
         + hs(0:nr,1:ns+1,0:np) - hs(0:nr,0:ns,0:np) &
         + hp(0:nr,0:ns,1:np+1) - hp(0:nr,0:ns,0:np)) / Vg(0:nr,0:ns,0:np)
    ! - include 0.5 from earlier average
    dhg(0:nr,0:ns,0:np) = 0.5_d * hyper_eta4 * dhg(0:nr,0:ns,0:np) &
        /bbgs(0:nr,0:ns,0:np)
    
    ! Set dhg to zero at pole
    IF (south) dhg(0:nr,0,0:np) = 0.0_d
    IF (north) dhg(0:nr,ns,0:np) = 0.0_d
    
    ! Set maximum timestep:
    local_dt = MINVAL(ABS(dnr(0:nr,0:ns,0:np) / dhg(0:nr,0:ns,0:np)))
    CALL MPI_ALLREDUCE(local_dt, dt_r, 1, mpireal, MPI_MIN, comm, err)
    local_dt = MINVAL(ABS(dns(0:nr,0:ns,0:np) / dhg(0:nr,0:ns,0:np)))
    CALL MPI_ALLREDUCE(local_dt, dt_s, 1, mpireal, MPI_MIN, comm, err)
    local_dt = MINVAL(ABS(dlp(0:nr,i0s:i1s,0:np) / dhg(0:nr,i0s:i1s,0:np)))
    CALL MPI_ALLREDUCE(local_dt, dt_p, 1, mpireal, MPI_MIN, comm, err)
    dt_max_hyper = MINVAL((/dt_r, dt_s, dt_p/)) * cflfact

END SUBROUTINE compute_hyper_div

!*******************************************************************************
END MODULE hyper
!*******************************************************************************
