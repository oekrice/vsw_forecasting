!*******************************************************************************
MODULE grid
!*******************************************************************************
! This module defines the simulation grid.
!*******************************************************************************
    USE params
    USE mpitools, ONLY: comm, err, mpireal, north, south, mpi_loc, MPI_SUM
    
    IMPLICIT NONE

    PRIVATE
    
    PUBLIC:: nr, ns, np, nrg, nsg, npg, rss, rg, rc, sg, sc, pg, pc, dr, ds, dp
    PUBLIC:: dlr, dls, dlp, dnr, dns, dnp, Sbr, Sbs, Sbp, Sbrg, Sbsg, Sbpg
    PUBLIC:: Sjr, Sjs, Sjp, Vg, V_total
    PUBLIC:: create_grid_arrays
    
!*******************************************************************************
    INTEGER:: nr, ns, np    ! Number of grid cells of this MPI process
                            ! in rho, s, phi directions
    INTEGER:: nrg, nsg, npg  ! Global number of grid cells in r, s, phi dirns
    REAL(d):: rss   ! Outer boundary radius [R_sun]
    REAL(d):: dr, ds, dp    ! Grid spacing in rho, s, phi [dimensionless]
    REAL(d), ALLOCATABLE:: rg(:), sg(:), pg(:)  ! rho, s, phi at grid points
    REAL(d), ALLOCATABLE:: rc(:), sc(:), pc(:)  ! rho, s, phi at cell centres
    REAL(d), ALLOCATABLE:: dlr(:,:,:)   ! Edge lengths in rho
    REAL(d), ALLOCATABLE:: dls(:,:,:)   ! Edge lengths in s
    REAL(d), ALLOCATABLE:: dlp(:,:,:)   ! Edge lengths in phi
    REAL(d), ALLOCATABLE:: dnr(:,:,:)   ! Staggered edge lengths in rho
    REAL(d), ALLOCATABLE:: dns(:,:,:)   ! Staggered edge lengths in s
    REAL(d), ALLOCATABLE:: dnp(:,:,:)   ! Staggered edge lengths in phi
    REAL(d), ALLOCATABLE:: Sbr(:,:,:)   ! Areas of rho-faces.
    REAL(d), ALLOCATABLE:: Sbs(:,:,:)   ! Areas of s-faces.
    REAL(d), ALLOCATABLE:: Sbp(:,:,:)   ! Areas of phi-faces.
    REAL(d), ALLOCATABLE:: Sjr(:,:,:)   ! Areas of staggered rho-faces.
    REAL(d), ALLOCATABLE:: Sjs(:,:,:)   ! Areas of staggered s-faces.
    REAL(d), ALLOCATABLE:: Sjp(:,:,:)   ! Areas of staggered phi-faces.
    REAL(d), ALLOCATABLE:: Sbrg(:,:,:)   ! Combined areas at grid points in rho.
    REAL(d), ALLOCATABLE:: Sbsg(:,:,:)   ! Combined areas at grid points in rho.
    REAL(d), ALLOCATABLE:: Sbpg(:,:,:)   ! Combined areas at grid points in rho.
    REAL(d), ALLOCATABLE:: Vg(:,:,:)   ! Cell volume centered at grid points.
    REAL(d):: V_total    ! Volume of whole domain [R_sun**3]

!*******************************************************************************
CONTAINS
!*******************************************************************************

!===============================================================================
SUBROUTINE create_grid_arrays
! Create local arrays of grid points, lengths, areas and volumes.

    ! Local variables:
    INTEGER:: i, j
    REAL(d):: local_tot

    ! Grid spacings:
    dr = LOG(rss)/REAL(nrg);
    ds = 2.0_d/REAL(nsg);
    dp = 2.0_d*PI/REAL(npg);

    ! Grid point arrays:
    ALLOCATE(rg(0:nr), sg(0:ns), pg(0:np))
    rg(0:nr) = (/ (REAL(mpi_loc(1)*nr + i)*dr, i = 0, nr) /)
    sg(0:ns) = (/ (-1.0_d + REAL(mpi_loc(2)*ns + i)*ds, i = 0, ns) /)
    pg(0:np) = (/ (REAL(mpi_loc(3)*np + i)*dp, i = 0, np) /)

    ! Cell centre arrays: [include ghost points for rho and phi]
    ALLOCATE(rc(0:nr+1), sc(0:ns+1), pc(0:np+1))
    rc(0:nr+1) = (/ ((DBLE(mpi_loc(1)*nr + i) - 0.5_d)*dr, i = 0, nr+1) /)
    sc(0:ns+1) = (/ (-1.0_d + (DBLE(mpi_loc(2)*ns + i)-0.5_d)*ds, i = 0, ns+1) /)
    pc(0:np+1) = (/ ((DBLE(mpi_loc(3)*np + i)-0.5_d)*dp, i = 0, np+1) /)
    IF (south) sc(0) = -1.0_d
    IF (north) sc(ns+1) = 1.0_d

    ! Edge lengths:
    ALLOCATE(dlr(0:nr+1,0:ns,0:np))
    DO i = 0, nr+1
     dlr(i,:,:) = EXP(rc(i)-0.5_d*dr)*(EXP(dr) - 1.0_d)
    END DO
    ALLOCATE(dls(0:nr,0:ns+1,0:np))
    DO i = 0, nr
        DO j = 2, ns-1
            dls(i,j,:) = EXP(rg(i))*(ASIN(sc(j)+0.5_d*ds) &
                - ASIN(sc(j)-0.5_d*ds))
        END DO
        IF (south) THEN
            dls(i,1,:) = EXP(rg(i))*(ASIN(sc(1)+0.5_d*ds) - ASIN(-1.0_d))
            dls(i,0,:) = dls(i,1,:)
        ELSE
            dls(i,1,:) = EXP(rg(i))*(ASIN(sc(1)+0.5_d*ds) - &
                ASIN(sc(1)-0.5_d*ds))
            dls(i,0,:) = EXP(rg(i))*(ASIN(sc(0)+0.5_d*ds) - &
                ASIN(sc(0)-0.5_d*ds))
        END IF
        IF (north) THEN
            dls(i,ns,:) = EXP(rg(i))*(ASIN(1.0_d) - ASIN(sc(ns)-0.5_d*ds))
            dls(i,ns+1,:) = dls(i,ns,:)
        ELSE
            dls(i,ns,:) =  EXP(rg(i))*(ASIN(sc(ns)+0.5_d*ds) - &
                ASIN(sc(ns)-0.5_d*ds))
            dls(i,ns+1,:) = EXP(rg(i))*(ASIN(sc(ns+1)+0.5_d*ds) - &
                ASIN(sc(ns+1)-0.5_d*ds))
        END IF
    END DO

    ALLOCATE(dlp(0:nr,0:ns,0:np+1))
    DO i = 0, nr
        DO j = 1, ns-1
            dlp(i,j,:) = EXP(rg(i))*SQRT(1.0_d - sg(j)**2)*dp
        END DO
        IF (south) THEN
            dlp(i,0,:) = 0.0_d
        ELSE
            dlp(i,0,:) = EXP(rg(i))*SQRT(1.0_d - sg(0)**2)*dp
        END IF
        IF (north) THEN
            dlp(i,ns,:) = 0.0_d
        ELSE
            dlp(i,ns,:) = EXP(rg(i))*SQRT(1.0_d - sg(ns)**2)*dp
        END IF
    END DO

    ! Normal lengths at face centres:
    ALLOCATE(dnr(0:nr,0:ns+1,0:np+1))
    DO i = 0, nr
     dnr(i,:,:) = EXP(rc(i))*(EXP(dr) - 1.0_d)
    END DO
    IF (south) dnr(:,0,:) = -dnr(:,0,:)
    IF (north) dnr(:,ns+1,:) = -dnr(:,ns+1,:)

    ALLOCATE(dns(0:nr+1,0:ns,0:np+1))
    DO i = 0, nr+1
        DO j = 1, ns-1
            dns(i,j,:) = EXP(rc(i))*(ASIN(sc(j+1)) - ASIN(sc(j)))
        END DO
        IF (north) THEN
            dns(i,ns,:) = 1.0_d   ! not used but needed to avoid NaN
        ELSE
            dns(i,ns,:) = exp(rc(i))*(asin(sc(ns+1)) - asin(sc(ns)))
        END IF
        IF (south) THEN
            dns(i,0,:) = 1.0_d
        ELSE
            dns(i,0,:) = EXP(rc(i))*(ASIN(sc(1)) - ASIN(sc(0)))
        END IF
    END DO

    ALLOCATE(dnp(0:nr+1,0:ns+1,0:np))
    DO i = 0, nr+1
        DO j = 1, ns
         dnp(i,j,:) = EXP(rc(i))*SQRT(1.0_d - sc(j)**2)*dp
        END DO
        IF (south) THEN
         dnp(i,0,:) = dnp(i,1,:)
        ELSE
         dnp(i,0,:) = EXP(rc(i))*SQRT(1.0_d - sc(0)**2)*dp
        END IF
        IF (north) THEN
         dnp(i,ns+1,:) = dnp(i,ns,:)
        ELSE
         dnp(i,ns+1,:) = EXP(rc(i))*SQRT(1.0_d - sc(ns+1)**2)*dp
        END IF
    END DO

    ! Face areas:
    ALLOCATE(Sbr(0:nr,0:ns+1,0:np+1))
    DO i = 0, nr
        Sbr(i,:,:) = EXP(2.0_d*rg(i))*ds*dp
    END DO

    ALLOCATE(Sbs(0:nr+1,0:ns,0:np+1))
    DO i = 0, nr+1
        DO j = 1, ns-1
         Sbs(i,j,:) = 0.5_d*EXP(2.0_d*rc(i) - dr)*dp* &
             (EXP(2.0_d*dr) - 1.0_d)*SQRT(1.0_d - sg(j)**2)
        END DO
        IF (south) THEN
         Sbs(i,0,:) = Sbs(i,1,:)
        ELSE
         Sbs(i,0,:) = 0.5_d*EXP(2.0_d*rc(i) - dr)*dp* &
             (EXP(2.0_d*dr) - 1.0_d)*SQRT(1.0_d - sg(0)**2)
        END IF
        IF (north) THEN
         Sbs(i,ns,:) = Sbs(i,ns-1,:)
        ELSE
         Sbs(i,ns,:) = 0.5_d*EXP(2.0_d*rc(i) - dr)*dp* &
             (EXP(2.0_d*dr) - 1.0_d)*SQRT(1.0_d - sg(ns)**2)
        END IF
    END DO

    ALLOCATE(Sbp(0:nr+1,0:ns+1,0:np))
    DO i = 0, nr+1
        DO j= 2, ns-1
            Sbp(i,j,:) = 0.5_d*EXP(2.0_d*rc(i) - dr)*(EXP(2.0_d*dr) - 1.0_d) &
                *(ASIN(sc(j)+0.5_d*ds) - ASIN(sc(j)-0.5*ds))
        END DO
        IF (south) THEN
            Sbp(i,1,:) = 0.5_d*EXP(2.0_d*rc(i) - dr)* &
                (EXP(2.0_d*dr) - 1.0_d)*(ASIN(sc(1)+0.5_d*ds) - ASIN(-1.0_d))
            Sbp(i,0,:) = Sbp(i,1,:)
        ELSE
            Sbp(i,1,:) = 0.5_d*EXP(2.0_d*rc(i) - dr)*(EXP(2.0_d*dr) - 1.0_d) &
                *(ASIN(sc(1)+0.5_d*ds) - ASIN(sc(1)-0.5*ds))
            Sbp(i,0,:) = 0.5_d*EXP(2.0_d*rc(i) - dr)*(EXP(2.0_d*dr) - 1.0_d) &
                *(ASIN(sc(0)+0.5_d*ds) - ASIN(sc(0)-0.5*ds))
        END IF
        IF (north) THEN
            Sbp(i,ns,:) = 0.5_d*EXP(2.0_d*rc(i) - dr)* &
                (EXP(2.0_d*dr) - 1.0_d)*(ASIN(1.0_d) - ASIN(sc(ns)-0.5*ds))
            Sbp(i,ns+1,:) = Sbp(i,ns,:)
        ELSE
            Sbp(i,ns,:) = 0.5_d*EXP(2.0_d*rc(i) - dr)*(EXP(2.0_d*dr) - 1.0_d) &
                *(ASIN(sc(ns)+0.5_d*ds) - ASIN(sc(ns)-0.5*ds))
            Sbp(i,ns+1,:) = 0.5_d*EXP(2.0_d*rc(i) - dr)*(EXP(2.0_d*dr) - 1.0_d) &
                *(ASIN(sc(ns+1)+0.5_d*ds) - ASIN(sc(ns+1)-0.5*ds))
        END IF
    END DO

    ! Areas perp to midpoints of edges:
    ALLOCATE(Sjr(0:nr+1,0:ns,0:np))
    DO i = 0, nr+1
     Sjr(i,:,:) = EXP(2.0_d*rc(i))*ds*dp
    END DO

    ALLOCATE(Sjs(0:nr,0:ns+1,0:np))
    DO i = 0, nr
        DO j = 1, ns
         Sjs(i,j,:) = 0.5_d*EXP(2.0_d*rc(i))*dp* &
             (EXP(2.0_d*dr) - 1.0_d)*SQRT(1.0_d - sc(j)**2)
        END DO
        IF (south) THEN
         Sjs(i,0,:) = -Sjs(i,1,:)
        ELSE
         Sjs(i,0,:) = 0.5_d*EXP(2.0_d*rc(i))*dp*(EXP(2.0_d*dr) - 1.0_d) &
            *SQRT(1.0_d - sc(0)**2)
        END IF
        IF (north) THEN
         Sjs(i,ns+1,:) = -Sjs(i,ns,:)
        ELSE
         Sjs(i,ns+1,:) = 0.5_d*EXP(2.0_d*rc(i))*dp*(EXP(2.0_d*dr) - 1.0_d) &
            *SQRT(1.0_d - sc(ns+1)**2)
        END IF
    END DO

    ALLOCATE(Sjp(0:nr,0:ns,0:np+1))
    DO i = 0, nr
        DO j = 1, ns-1
         Sjp(i,j,:) = 0.5_d*EXP(2.0_d*rc(i))*(EXP(2.0_d*dr) - 1.0_d) &
            *(ASIN(sc(j+1)) - ASIN(sc(j)))
        END DO
        IF (south) THEN
            Sjp(i,0,:) = Sjp(i,1,:)
        ELSE
            Sjp(i,0,:) = 0.5_d*EXP(2.0_d*rc(i))*(exp(2.0_d*dr) - 1.0_d) &
                *(ASIN(sc(1)) - ASIN(sc(0)))
        END IF
        IF (north) THEN
            Sjp(i,ns,:) = Sjp(i,ns-1,:)
        ELSE
            Sjp(i,ns,:) = 0.5_d*EXP(2.0_d*rc(i))*(EXP(2.0_d*dr) - 1.0_d) &
                *(ASIN(sc(ns+1)) - ASIN(sc(ns)))
        END IF
    END DO

    ! Combined areas at grid points (for quick averaging):
    ALLOCATE(Sbrg(0:nr,0:ns,0:np),Sbsg(0:nr,0:ns,0:np),Sbpg(0:nr,0:ns,0:np))
    Sbrg = Sbr(:,0:ns,0:np) + Sbr(:,1:ns+1,0:np) + Sbr(:,0:ns,1:np+1) &
        + Sbr(:,1:ns+1,1:np+1)
    Sbsg = Sbs(0:nr,:,0:np) + Sbs(1:nr+1,:,0:np) + Sbs(0:nr,:,1:np+1) &
        + Sbs(1:nr+1,:,1:np+1)
    Sbpg = Sbp(0:nr,0:ns,:) + Sbp(0:nr,1:ns+1,:) + Sbp(1:nr+1,0:ns,:) &
        + Sbp(1:nr+1,1:ns+1,:)

    ! Cell volume (centered at grid points):
    ALLOCATE(Vg(0:nr,0:ns,0:np))
    DO i = 0, nr
        Vg(i,:,:) = EXP(3.0_d*rc(i))*ds*dp*(EXP(3.0_d*dr) - 1.0_d)/3.0_d
    END DO

    ! Compute volume of whole domain using trapezium rule:
    local_tot = 0.125_d*SUM( &
     Vg(0:nr-1,0:ns-1,0:np-1) + Vg(0:nr-1,0:ns-1,1:np) &
     + Vg(0:nr-1,1:ns,0:np-1) + Vg(0:nr-1,1:ns,1:np) &
     + Vg(1:nr,0:ns-1,0:np-1) + Vg(1:nr,0:ns-1,1:np) &
     + Vg(1:nr,1:ns,0:np-1) + Vg(1:nr,1:ns,1:np))
    CALL MPI_ALLREDUCE(local_tot, V_total, 1, mpireal, MPI_SUM, comm, err)

END SUBROUTINE create_grid_arrays

!*******************************************************************************
END MODULE grid
!*******************************************************************************
