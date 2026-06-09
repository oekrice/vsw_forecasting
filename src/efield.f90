!*******************************************************************************
MODULE efield
!*******************************************************************************
! This module implements the electric field.
!*******************************************************************************
    USE params
    USE mpitools, ONLY: control, surface, outer, north, south, i1r, mpireal, &
        MPI_MAX, comm, stat, err, mpi_r_prv, mpi_r_nxt, mpi_s_prv, mpi_s_nxt, &
        mpi_p_prv, mpi_p_nxt
    USE grid, ONLY: nr, ns, np, dlr, dls, dlp, dnr, Sbr, Sjs, Sjp, Sbs, Sbp
    USE friction, ONLY: nu
    USE driver, ONLY: vbs, vbp, eta_surf, als_em, alp_em, t_em
    USE ohmic, ONLY: ohmic_enhanced, ohmic_eta
    USE hyper, ONLY: compute_hyper_div, dhg

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: elr, els, elp, egr, egs, egp, fgr, fgs, fgp
    PUBLIC:: zero_efield, add_outflow, add_friction, add_ohmic
    PUBLIC:: add_hyper, add_surface_flows, add_surface_diffusion
    PUBLIC:: add_emergence, mpi_transfer_e, compute_eg_average

!*******************************************************************************
    REAL(d), ALLOCATABLE:: elr(:,:,:)   ! E_rho*L_rho on edges (1:nr,0:ns,0:np)
    REAL(d), ALLOCATABLE:: els(:,:,:)   ! E_s*L_s on edges (0:nr,1:ns,0:np)
    REAL(d), ALLOCATABLE:: elp(:,:,:)   ! E_phi*L_phi on edges (0:nr,0:ns,1:np)
    REAL(d), ALLOCATABLE:: egr(:,:,:)   ! E_rho at grid pts (0:nr,0:ns,0:np)
    REAL(d), ALLOCATABLE:: egs(:,:,:)   ! E_s at grid pts (0:nr,0:ns,0:np)
    REAL(d), ALLOCATABLE:: egp(:,:,:)   ! E_phi at grid pts (0:nr,0:ns,0:np)
    REAL(d), ALLOCATABLE:: fgr(:,:,:)   ! Force F_r at grid pts (0:nr,0:ns,0:np)
    REAL(d), ALLOCATABLE:: fgs(:,:,:)   ! Force F_s at grid pts (0:nr,0:ns,0:np)
    REAL(d), ALLOCATABLE:: fgp(:,:,:) ! Force F_phi at grid pts (0:nr,0:ns,0:np)
    REAL(d), ALLOCATABLE:: jjg(:,:,:) ! |J| at grid pts (0:nr,0:ns,0:np)
    INTEGER, PARAMETER:: tag=1
    
!*******************************************************************************
CONTAINS
!*******************************************************************************

!===============================================================================
SUBROUTINE zero_efield()
! Allocate E*L on cell edges if necessary, and set to zero.
    
    IF (.NOT. ALLOCATED(elr)) ALLOCATE(elr(1:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(els)) ALLOCATE(els(0:nr,1:ns,0:np))
    IF (.NOT. ALLOCATED(elp)) ALLOCATE(elp(0:nr,0:ns,1:np))

    elr(1:nr,0:ns,0:np) = 0.0_d
    els(0:nr,1:ns,0:np) = 0.0_d
    elp(0:nr,0:ns,1:np) = 0.0_d
    
    IF (.NOT. ALLOCATED(egr)) ALLOCATE(egr(0:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(egs)) ALLOCATE(egs(0:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(egp)) ALLOCATE(egp(0:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(fgr)) ALLOCATE(fgr(0:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(fgs)) ALLOCATE(fgs(0:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(fgp)) ALLOCATE(fgp(0:nr,0:ns,0:np))
    IF (.NOT. ALLOCATED(jjg)) ALLOCATE(jjg(0:nr,0:ns,0:np))

END SUBROUTINE zero_efield

!===============================================================================
SUBROUTINE add_outflow(bs, bp, vr)
! Use upwinding to compute -v_r*B_phi on s-edges and v_r*B_s on p-edges.
! Assumes that v_r > 0 (outward flow).
    REAL(d), INTENT(IN):: bs(0:nr+1,0:ns,0:np+1)    ! B_s on faces
    REAL(d), INTENT(IN):: bp(0:nr+1,0:ns+1,0:np)    ! B_phi on faces
    REAL(d), INTENT(IN):: vr(0:nr)  ! outflow velocity profile [R_sun/s]

    ! Local variables:
    INTEGER:: i
    
    DO i = i1r, nr
        els(i,1:ns,0:np) = els(i,1:ns,0:np)  &
            - vr(i)*dls(i,1:ns,0:np) * bp(i,1:ns,0:np)
        elp(i,0:ns,1:np) = elp(i,0:ns,1:np) &
            + vr(i)*dlp(i,0:ns,1:np) * bs(i,0:ns,1:np)
    END DO

END SUBROUTINE add_outflow

!===============================================================================
SUBROUTINE add_friction(bgr, bgs, bgp, bbgs, jgr, jgs, jgp)
! Compute electric field for friction and add to elr, els, elp.
    REAL(d), INTENT(IN):: bgr(0:nr,0:ns,0:np)    ! B_r at grid pts (unweighted)
    REAL(d), INTENT(IN):: bgs(0:nr,0:ns,0:np)    ! B_s at grid pts (unweighted)
    REAL(d), INTENT(IN):: bgp(0:nr,0:ns,0:np)   ! B_phi at grid pts (unweighted)
    REAL(d), INTENT(IN):: bbgs(-1:nr+1,0:ns,0:np)    ! Softened B^2 at grid pts
    REAL(d), INTENT(IN):: jgr(0:nr,0:ns,0:np)    ! J_r at grid pts (unweighted)
    REAL(d), INTENT(IN):: jgs(0:nr,0:ns,0:np)    ! J_s at grid pts (unweighted)
    REAL(d), INTENT(IN):: jgp(0:nr,0:ns,0:np)   ! J_phi at grid pts (unweighted)

    ! Lorentz force at grid pts:
    fgr(0:nr,0:ns,0:np) = jgp(0:nr,0:ns,0:np) * bgs(0:nr,0:ns,0:np) &
        - jgs(0:nr,0:ns,0:np) * bgp(0:nr,0:ns,0:np)
    fgs(0:nr,0:ns,0:np) = jgr(0:nr,0:ns,0:np) * bgp(0:nr,0:ns,0:np) &
            - jgp(0:nr,0:ns,0:np) * bgr(0:nr,0:ns,0:np)
    fgp(0:nr,0:ns,0:np) = jgs(0:nr,0:ns,0:np) * bgr(0:nr,0:ns,0:np) &
        - jgr(0:nr,0:ns,0:np) * bgs(0:nr,0:ns,0:np)
    
    ! Electric field -nu/soften(B**2) * (jxB):
    egr(0:nr,0:ns,0:np) = nu(0:nr,0:ns,0:np) * (fgs(0:nr,0:ns,0:np) * &
        bgp(0:nr,0:ns,0:np) - fgp(0:nr,0:ns,0:np) * bgs(0:nr,0:ns,0:np)) &
        / bbgs(0:nr,0:ns,0:np)
    egs(0:nr,0:ns,0:np) = nu(0:nr,0:ns,0:np) * (fgp(0:nr,0:ns,0:np) * &
        bgr(0:nr,0:ns,0:np) - fgr(0:nr,0:ns,0:np) * bgp(0:nr,0:ns,0:np)) &
        / bbgs(0:nr,0:ns,0:np)
    egp(0:nr,0:ns,0:np) = nu(0:nr,0:ns,0:np) * (fgr(0:nr,0:ns,0:np) * &
        bgs(0:nr,0:ns,0:np) - fgs(0:nr,0:ns,0:np) * bgr(0:nr,0:ns,0:np)) &
        / bbgs(0:nr,0:ns,0:np)

    ! Average to edges and add scale factors:
    elr(1:nr,0:ns,0:np) = elr(1:nr,0:ns,0:np) + 0.5_d * ( &
        egr(0:nr-1,0:ns,0:np) + egr(1:nr,0:ns,0:np)) * dlr(1:nr,0:ns,0:np)
    els(i1r:nr,1:ns,0:np) = els(i1r:nr,1:ns,0:np) + 0.5_d * ( &
        egs(i1r:nr,0:ns-1,0:np) + egs(i1r:nr,1:ns,0:np)) * dls(i1r:nr,1:ns,0:np)
    elp(i1r:nr,0:ns,1:np) = elp(i1r:nr,0:ns,1:np) + 0.5_d * ( &
        egp(i1r:nr,0:ns,0:np-1) + egp(i1r:nr,0:ns,1:np)) * dlp(i1r:nr,0:ns,1:np)
    
END SUBROUTINE add_friction

!===============================================================================
SUBROUTINE add_ohmic(jr, js, jp, jgr, jgs, jgp, bbgs)
! Compute electric field for ohmic diffusion and add to elr, els, elp.
    REAL(d), INTENT(IN):: jr(0:nr+1,0:ns,0:np)  ! J_r on stagd faces (unwtd)
    REAL(d), INTENT(IN):: js(0:nr,0:ns+1,0:np)  ! J_s on stagd faces (unwtd)
    REAL(d), INTENT(IN):: jp(0:nr,0:ns,0:np+1)  ! J_phi on stagd faces (unwtd)
    REAL(d), INTENT(IN):: jgr(0:nr,0:ns,0:np)    ! J_r at grid pts (unweighted)
    REAL(d), INTENT(IN):: jgs(0:nr,0:ns,0:np)    ! J_s at grid pts (unweighted)
    REAL(d), INTENT(IN):: jgp(0:nr,0:ns,0:np)   ! J_phi at grid pts (unweighted)
    REAL(d), INTENT(IN):: bbgs(-1:nr+1,0:ns,0:np)    ! softened B^2 at grid pts

    ! Local variables:
    REAL(d):: local_bmax, glob_bmax
    
    ! Compute maximum |B|:
    local_bmax = MAXVAL(bbgs(0:nr,0:ns,0:np))
    CALL MPI_ALLREDUCE(local_bmax, glob_bmax, 1, mpireal, MPI_MAX, comm, err)
    glob_bmax = 0.5_d * ohmic_enhanced / SQRT(glob_bmax)

    ! Compute |J| at grid points:
    jjg(0:nr,0:ns,0:np) = SQRT(jgr(0:nr,0:ns,0:np) * jgr(0:nr,0:ns,0:np) &
        + jgs(0:nr,0:ns,0:np) * jgs(0:nr,0:ns,0:np) &
        + jgp(0:nr,0:ns,0:np) * jgp(0:nr,0:ns,0:np))

    ! TEMP [COPY GLOVAG]
    WHERE (jjg(0:nr,0:ns,0:np) < 1.0d-3) jjg(0:nr,0:ns,0:np) = 1.0d-3

    ! Add diffusion eta*J on coronal edges:
    elr(1:nr,0:ns,0:np) = elr(1:nr,0:ns,0:np) + ohmic_eta * (1.0_d + &
        glob_bmax * (jjg(0:nr-1,0:ns,0:np) + jjg(1:nr,0:ns,0:np))) &
        * jr(1:nr,0:ns,0:np) * dlr(1:nr,0:ns,0:np)
         
    els(i1r:nr,1:ns,0:np) = els(i1r:nr,1:ns,0:np) + ohmic_eta*(1.0_d + &
        glob_bmax * (jjg(i1r:nr,0:ns-1,0:np) + jjg(i1r:nr,1:ns,0:np))) &
        * js(i1r:nr,1:ns,0:np) * dls(i1r:nr,1:ns,0:np)
         
    elp(i1r:nr,0:ns,1:np) = elp(i1r:nr,0:ns,1:np) + ohmic_eta*(1.0_d + &
        glob_bmax * (jjg(i1r:nr,0:ns,0:np-1) + jjg(i1r:nr,0:ns,1:np))) &
        * jp(i1r:nr,0:ns,1:np) * dlp(i1r:nr,0:ns,1:np)

END SUBROUTINE add_ohmic

!===============================================================================
SUBROUTINE add_hyper(bgr, bgs, bgp, bbgs)
! Compute electric field for hyperdiffusion and add to elr, els, elp.
    REAL(d), INTENT(IN):: bgr(0:nr,0:ns,0:np)  ! B_rho at grid pts (unweighted)
    REAL(d), INTENT(IN):: bgs(0:nr,0:ns,0:np)  ! B_s at grid pts (unwtd)
    REAL(d), INTENT(IN):: bgp(0:nr,0:ns,0:np)  ! B_phi at grid pts (unwtd)
    REAL(d), INTENT(IN):: bbgs(-1:nr+1,0:ns,0:np)  ! Softened B**2 at grid pts

    ! Compute div( eta4*B**2*grad(alg) )/B**2 at grid points:
    CALL compute_hyper_div(bbgs)
    
    ! Compute electric field at grid points:
    egr(0:nr,0:ns,0:np) = -bgr(0:nr,0:ns,0:np)*dhg(0:nr,0:ns,0:np)
    egs(0:nr,0:ns,0:np) = -bgs(0:nr,0:ns,0:np)*dhg(0:nr,0:ns,0:np)
    egp(0:nr,0:ns,0:np) = -bgp(0:nr,0:ns,0:np)*dhg(0:nr,0:ns,0:np)

    ! Average to edges and add scale factors:
    elr(1:nr,0:ns,0:np) = elr(1:nr,0:ns,0:np) + 0.5_d * ( &
        egr(0:nr-1,0:ns,0:np) + egr(1:nr,0:ns,0:np)) * dlr(1:nr,0:ns,0:np)
    els(i1r:nr,1:ns,0:np) = els(i1r:nr,1:ns,0:np) + 0.5_d * ( &
        egs(i1r:nr,0:ns-1,0:np) + egs(i1r:nr,1:ns,0:np)) * dls(i1r:nr,1:ns,0:np)
    elp(i1r:nr,0:ns,1:np) = elp(i1r:nr,0:ns,1:np) + 0.5_d * ( &
        egp(i1r:nr,0:ns,0:np-1) + egp(i1r:nr,0:ns,1:np)) * dlp(i1r:nr,0:ns,1:np)

END SUBROUTINE add_hyper

!===============================================================================
SUBROUTINE add_surface_flows(br)
! Compute electric field for surface flows and add to els, elp.
    REAL(d), INTENT(IN):: br(0:nr,0:ns+1,0:np+1)  ! B_rho on faces

    IF (surface) THEN
        els(0,1:ns,0:np) = els(0,1:ns,0:np) + vbp(1:ns,0:np) * 0.5_d * &
            (br(0,1:ns,0:np) + br(0,1:ns,1:np+1)) * dls(0,1:ns,0:np)
        elp(0,0:ns,1:np) = elp(0,0:ns,1:np) - vbs(0:ns,1:np) * 0.5_d * &
            (br(0,0:ns,1:np) + br(0,1:ns+1,1:np)) * dlp(0,0:ns,1:np)
    END IF

END SUBROUTINE add_surface_flows

!===============================================================================
SUBROUTINE add_surface_diffusion(br)
! Compute electric field for supergranular diffusion and add to els, elp.
    REAL(d), INTENT(IN):: br(0:nr,0:ns+1,0:np+1)  ! B_rho* on faces

    ! Local variables:
    INTEGER:: i, j

    IF (surface) THEN
        els(0,1:ns,0:np) = els(0,1:ns,0:np) + eta_surf * dls(0,1:ns,0:np) &
            * (br(0,1:ns,0:np) * dnr(0,1:ns,0:np) - br(0,1:ns,1:np+1) &
            * dnr(0,1:ns,1:np+1)) / Sjs(0,1:ns,0:np)
        IF (south) THEN
            i = 1
        ELSE
            i = 0
        END IF
        IF (north) THEN
            j = ns-1
        ELSE
            j = ns
        END IF
        elp(0,i:j,1:np) = elp(0,i:j,1:np) + eta_surf * dlp(0,i:j,1:np) &
            * (br(0,i+1:j+1,1:np) * dnr(0,i+1:j+1,1:np) - br(0,i:j,1:np) &
            * dnr(0,i:j,1:np)) / Sjp(0,i:j,1:np)
    END IF
    
END SUBROUTINE add_surface_diffusion

!===============================================================================
SUBROUTINE add_emergence()
! Add surface electric field for new region emergences to els, elp.
! A constant electric field is applied during the emergence.

    IF (surface) THEN
        els(0,1:ns,0:np) = els(0,1:ns,0:np) - als_em(1:ns,0:np)/t_em
        elp(0,0:ns,1:np) = elp(0,0:ns,1:np) - alp_em(0:ns,1:np)/t_em
    END IF

END SUBROUTINE add_emergence

!===============================================================================
SUBROUTINE mpi_transfer_e()
! Transfer values or elr, els, elp on MPI boundaries to the previous MPI process
! and overwrite values that were previously there (to avoid accumulation of
! differences between processes).

    ! Transfer er and es in p direction:
    CALL MPI_SENDRECV(elr(1:nr,0:ns,0), nr * (ns+1), mpireal, mpi_p_prv, tag, &
        elr(1:nr,0:ns,np), nr * (ns+1), mpireal, mpi_p_nxt, tag, comm, stat, &
        err)
    CALL MPI_SENDRECV(els(0:nr,1:ns,0), (nr+1) * ns, mpireal, mpi_p_prv, tag, &
        els(0:nr,1:ns,np), (nr+1) * ns, mpireal, mpi_p_nxt, tag, comm, stat, &
        err)
        
    ! Transfer er and ep in s direction:
    IF (.NOT. south) CALL MPI_SEND(elr(1:nr,0,0:np), nr * (np+1), mpireal, &
        mpi_s_prv, tag, comm, err)
    IF (.NOT. north) CALL MPI_RECV(elr(1:nr,ns,0:np), nr * (np+1), mpireal, &
        mpi_s_nxt, tag, comm, stat, err)
    IF (.NOT. south) CALL MPI_SEND(elp(0:nr,0,1:np), (nr+1) * np, mpireal, &
        mpi_s_prv, tag, comm, err)
    IF (.NOT. north) CALL MPI_RECV(elp(0:nr,ns,1:np), (nr+1) * np, mpireal, &
        mpi_s_nxt, tag, comm, stat, err)
            
    ! Transfer es and ep in r direction:
    IF (.NOT. surface) CALL MPI_SEND(els(0,1:ns,0:np), ns * (np+1), mpireal, &
        mpi_r_prv, tag, comm, err)
    IF (.NOT. outer) CALL MPI_RECV(els(nr,1:ns,0:np), ns * (np+1), mpireal, &
        mpi_r_nxt, tag, comm, stat, err)
    IF (.NOT. surface) CALL MPI_SEND(elp(0,0:ns,1:np), (ns+1) * np, mpireal, &
        mpi_r_prv, tag, comm, err)
    IF (.NOT. outer) CALL MPI_RECV(elp(nr,0:ns,1:np), (ns+1) * np, mpireal, &
        mpi_r_nxt, tag, comm, stat, err)

END SUBROUTINE mpi_transfer_e

!===============================================================================
SUBROUTINE compute_eg_average()
! Average elr, els, elp to grid pts and remove length factors.
! This routine is usually used just for outputting the electric field.

    ! Local variables:
    INTEGER:: i0, i1
    
    ! First and last interior grid points (to avoid dividing by zero dlp):
    IF (south) THEN
        i0 = 1
    ELSE
        i0 = 0
    END IF
    IF (north) THEN
        i1 = ns-1
    ELSE
        i1 = ns
    END IF
    
    ! Average E to interior grid points.
    egr(1:nr-1,0:ns,0:np) = 0.5_d * (elr(1:nr-1,0:ns,0:np) &
        / dlr(1:nr-1,0:ns,0:np) + elr(2:nr,0:ns,0:np) / dlr(2:nr,0:ns,0:np))
    egs(0:nr,1:ns-1,0:np) = 0.5_d * (els(0:nr,1:ns-1,0:np) &
        / dls(0:nr,1:ns-1,0:np) + els(0:nr,2:ns,0:np) / dls(0:nr,2:ns,0:np))
    egp(0:nr,i0:i1,1:np-1) = 0.5_d * (elp(0:nr,i0:i1,1:np-1) &
        / dlp(0:nr,i0:i1,1:np-1) + elp(0:nr,i0:i1,2:np) / dlp(0:nr,i0:i1,2:np))

    ! Fill global boundary values by simple extrapolation:
    IF (surface) egr(0,0:ns,0:np) = 2.0_d * elr(1,0:ns,0:np) &
        / dlr(1,0:ns,0:np) - elr(2,0:ns,0:np) / dlr(2,0:ns,0:np)
    IF (outer) egr(nr,0:ns,0:np) = 2.0_d * elr(nr-1,0:ns,0:np) / &
        dlr(nr-1,0:ns,0:np) - elr(nr-2,0:ns,0:np) / dlr(nr-2,0:ns,0:np)

    ! Use MPI to exchange interior ghost values:
    ! - transfer ep in p direction:
    CALL MPI_SENDRECV(elp(0:nr,i0:i1,1) / dlp(0:nr,i0:i1,1), &
        (nr+1) * (i1-i0+1), mpireal, mpi_p_prv, tag, egp(0:nr,i0:i1,np), &
        (nr+1) * (i1-i0+1), mpireal, mpi_p_nxt, tag, comm, stat, err)
    egp(0:nr,i0:i1,np) = 0.5_d * (egp(0:nr,i0:i1,np) + elp(0:nr,i0:i1,np) &
        / dlp(0:nr,i0:i1,np))
    CALL MPI_SENDRECV(elp(0:nr,i0:i1,np) / dlp(0:nr,i0:i1,np), &
        (nr+1) * (i1-i0+1), mpireal, mpi_p_nxt, tag, egp(0:nr,i0:i1,0), &
        (nr+1) * (i1-i0+1), mpireal, mpi_p_prv, tag, comm, stat, err)
    egp(0:nr,i0:i1,0) = 0.5_d * (egp(0:nr,i0:i1,0) + elp(0:nr,i0:i1,1) &
        / dlp(0:nr,i0:i1,1))
        
    ! Transfer es in s direction:
    IF (.NOT. south) CALL MPI_SEND(els(0:nr,1,0:np) / dls(0:nr,1,0:np), &
        (nr+1) * (np+1), mpireal, mpi_s_prv, tag, comm, err)
    IF (.NOT. north) CALL MPI_RECV(egs(0:nr,ns,0:np), (nr+1) * (np+1), &
        mpireal, mpi_s_nxt, tag, comm, stat, err)
    egs(0:nr,ns,0:np) = 0.5_d * (egs(0:nr,ns,0:np) + els(0:nr,ns,0:np) &
            / dls(0:nr,ns,0:np))
    IF (.NOT. north) CALL MPI_SEND(els(0:nr,ns,0:np) / dls(0:nr,ns,0:np), &
        (nr+1) * (np+1), mpireal, mpi_s_nxt, tag, comm, err)
    IF (.NOT. south) CALL MPI_RECV(egs(0:nr,0,0:np), (nr+1) * (np+1), mpireal, &
            mpi_s_prv, tag, comm, stat, err)
    egs(0:nr,0,0:np) = 0.5_d * (egs(0:nr,0,0:np) + els(0:nr,1,0:np) &
        / dls(0:nr,1,0:np))
            
    ! Transfer er in r direction:
    IF (.NOT. surface) CALL MPI_SEND(elr(1,0:ns,0:np) / dlr(1,0:ns,0:np), &
        (ns+1) * (np+1), mpireal, mpi_r_prv, tag, comm, err)
    IF (.NOT. outer) CALL MPI_RECV(egr(nr,0:ns,0:np), (ns+1) * (np+1), &
        mpireal, mpi_r_nxt, tag, comm, stat, err)
    egr(nr,0:ns,0:np) = 0.5_d * (egr(nr,0:ns,0:np) + elr(nr,0:ns,0:np) &
        / dlr(nr,0:ns,0:np))
    IF (.NOT. outer) CALL MPI_SEND(elr(nr,0:ns,0:np) / dlr(nr,0:ns,0:np), &
        (ns+1) * (np+1), mpireal, mpi_r_nxt, tag, comm, err)
    IF (.NOT. surface) CALL MPI_RECV(egr(0,0:ns,0:np), (ns+1) * (np+1), &
        mpireal, mpi_r_prv, tag, comm, stat, err)
    egr(0,0:ns,0:np) = 0.5_d * (egr(0,0:ns,0:np) + elr(1,0:ns,0:np) &
        / dlr(1,0:ns,0:np))

END SUBROUTINE compute_eg_average

!*******************************************************************************
END MODULE efield
!*******************************************************************************
