!****************************************************************
module shared

    use netcdf
  
!----------------------------------------------------------------
! Global constants and shared supporting routines.
!----------------------------------------------------------------
    implicit none
  
    logical, parameter :: VERBOSE=.true.

    integer, parameter :: d=kind(0.d0)    ! double precision
    real(d), parameter :: FOURPI = 16.0_d*datan(1.0_d)
    real(d), parameter :: TWOPI=8.0_d*datan(1.0_d)
    real(d), parameter :: PI = 4.0_d*datan(1.0_d)
    real(d), parameter :: RSUN = 6.96d10
    
    integer :: nr, ns, np
    real(d), dimension(:), allocatable :: r, s, p
    real(d) :: dr, ds, dp

    contains
    
    !****************************************************************
    subroutine chk(istatus)
        !----------------------------------------------------------------
        ! Chk (ever so slightly modified from www.unidata.ucar.edu).
        ! For netcdf.
        !----------------------------------------------------------------
        integer, intent (in) :: istatus
        ! ---
        ! Wrapper for netcdf calls (slightly modified from www.unidata.ucar.edu).
        ! ---
        if (istatus /= nf90_noerr) then
        write(*,*) trim(adjustl(nf90_strerror(istatus)))
        end if
        
    end subroutine chk   
    
    !****************************************************************
    subroutine makeGridSphere(radius, nsg, npg, x0)
        integer, intent(in) :: nsg, npg
        real(d), intent(in) :: radius
        real(d), intent(inout), dimension(:,:) :: x0
        ! ---
        ! Generate 2d grid at resolution (nsg, npg) on surface r=radius,
        ! on spherical surface.
        ! Output: gridpoints in Cartesian coordinates.
        ! ---
        real(d) :: dsg, dpg
        real(d), dimension(:), allocatable :: r0, s0, p0
        integer :: i,j

        dsg = 2.0_d/float(nsg)
        dpg = TWOPI/float(npg)
        
        allocate(r0(nsg*npg),s0(nsg*npg),p0(nsg*npg))

        !$omp parallel
        !$omp do
        do i=1,nsg
        do j=1,npg
            r0(j + (i-1)*npg) = radius
            s0(j + (i-1)*npg) = dble(i-1)/dble(nsg-1)*(s(ns)-s(1)-dsg) + s(1) + 0.5_d*dsg
            p0(j + (i-1)*npg) = dble(j-1)/dble(npg-1)*(p(np) - p(1)-dpg) + 0.5*dpg
        end do
        end do
        !$omp end do
        !$omp end parallel

        x0(:,1) = r0*sqrt(1.0_d - s0*s0)*cos(p0)
        x0(:,2) = r0*sqrt(1.0_d - s0*s0)*sin(p0)
        x0(:,3) = r0*s0

    end subroutine makeGridSphere
        
    !****************************************************************
    function pythag(a, b)
        real(d) :: a, b, pythag
        ! ---
        ! (from Numerical Recipes)
        ! Compute (a^2 + b^2)^1/2 without destructive underflow or overflow.
        ! ---
        real(d) :: absa, absb
        
        absa=abs(a)
        absb=abs(b)
        if (absa.gt.absb) then
        pythag=absa*dsqrt(1.0_d+(absb/absa)**2)
        else
        if (absb.eq.0.0_d) then
            pythag=0.0_d
        else
            pythag=absb*dsqrt(1.0_d+(absa/absb)**2)
        end if
        end if
        return
        
    end function pythag
        
    !****************************************************************
    function trap3d(fg)
        real(d), dimension(:,:,:), intent(in) :: fg
        real(d) :: trap3d
        ! ---
        ! Compute volume integral of array given at grid points, using
        ! trapezium rule.
        ! ---
        real(d), dimension(:), allocatable :: fg1
        integer :: i, nr

        nr = size(r)

        ! Do 2d integrals (all cells have equal area):
        allocate(fg1(nr))
        do i=1,nr
            fg1(i) = trap2d(fg(i,:,:))*exp(2.0_d*r(i))
        end do

        ! Do radial integral:
        trap3d = 0.0_d
        do i=1,nr-1
            trap3d = trap3d + 0.5_d*(fg1(i+1) + fg1(i))*(exp(r(i+1)) - exp(r(i)))
        end do

        return
    end function trap3d

    !****************************************************************
    function trap2d(fg)
        real(d), dimension(:,:), intent(in) :: fg
        real(d) :: trap2d
        ! ---
        ! Compute integral of array over constant-r given at grid points, using
        ! trapezium rule.
        ! ---
        integer :: ns, np
        ns = size(s)
        np = size(p)

        trap2d = sum(fg(1:ns-1,1:np-1) + fg(1:ns-1,2:np) &
            + fg(2:ns,1:np-1) + fg(2:ns,2:np))
        trap2d = trap2d*0.25_d*ds*dp

        return
    end function trap2d

    !****************************************************************
    subroutine trieig(dv, e, z)
        real(d), dimension(:), intent(inout) :: dv,e
        real(d), dimension(:,:), optional, intent(inout) :: z
        ! ---
        ! (tqli from Numerical Recipes)
        ! Compute eigenvalues and eigenvectors of a symm. tridiagonal matrix A, using
        ! the QL algorithm. To find only eigenvalues, omit optional argument z.
        ! - input: dv(n) - a vector with the diagonal elements of A
        !          e(n) - a vector with the subdiagonal of A (ignores e(1))
        !          z(n,n) - the identity matrix [optional - if eigenvectors reqd]
        ! - output: dv(n) - a vector with the eigenvalues
        !           z(n,n) - the matrix of eigenvectors - kth column is
        !                    normalized eigenvector corresponding to dv(k).
        ! ---
        integer :: i,iter,l,m,n,ndum
        real(d) :: b,c,dd,f,g,p,r,s
        real(d), dimension(size(e)) :: ff
        
        n = size(dv)
        if (present(z)) ndum = size(z,1)
        e(:) = eoshift(e(:),1)  ! convenient to renumber the elements of e
        do l=1,n
        iter=0
        iterate: do
            do m=l,n-1  ! Look for a single small subdiagonal element to split
                ! the matrix
                dd=abs(dv(m)) + abs(dv(m+1))
                if (abs(e(m))+dd == dd) exit
            end do
            if (m == l) exit iterate
            if (iter == 30) then
                print*,'Error: too many iterations in tqli'
                stop
            end if
            iter=iter+1
            g=(dv(l+1)-dv(l))/(2.0_d*e(l))  ! form shift
            r=pythag(g, 1.0_d)
            g=dv(m)-dv(l)+e(l)/(g+sign(r,g))  ! This is dv_m - k_s
            s=1.0_d
            c=1.0_d
            p=0.0_d
            do i=m-1,l,-1  ! A plane rotation as in original QL, followed by
                ! Givens rotations to restore tridiagonal form.
                f=s*e(i)
                b=c*e(i)
                r=pythag(f,g)
                e(i+1)=r
                if (r == 0.0_d) then  ! recover from underflow
                    dv(i+1) = dv(i+1) - p
                    e(m) = 0.0_d
                    cycle iterate
                end if
                s=f/r
                c=g/r
                g=dv(i+1)-p
                r=(dv(i)-g)*s+2.0_d*c*b
                p=s*r
                dv(i+1)=g+p
                g=c*r-b
                if (present(z)) then  ! for eigenvectors
                    ff(1:n)=z(1:n,i+1)
                    z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
                    z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
                end if
            end do
            dv(l)=dv(l)-p
            e(l)=g
            e(m)=0.0_d
        end do iterate
        end do
        
    end subroutine trieig
    
end module shared
