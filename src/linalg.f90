module linalg
   !! Wrapper module around LAPACK's `dgesv` and lines intersections detector
   use constants, only: pr
   implicit none

   type :: point
      real(pr) :: x
      real(pr) :: y
      integer :: i
      integer :: j
   end type point

   interface intersection
      module procedure :: intersect_two_lines
      module procedure :: intersect_one_line
   end interface

   interface allclose
      module procedure :: allclose_inner
      module procedure :: allclose_outer
   end interface
contains
   function solve_system(a, b) result(x)
      real(pr), intent(in) :: b(:)
      real(pr), intent(in) :: a(size(b), size(b))

      real(pr) :: x(size(b))

      real(8) :: a_lapack(size(b), size(b)), b_lapack(size(b))

      integer :: n, nrhs, lda, ipiv(size(b)), ldb, info

      interface
         subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            integer :: n
            integer :: nrhs
            real(8) :: a(n, n)
            integer :: lda
            integer :: ipiv(n)
            real(8) :: b(n)
            integer :: ldb
            integer :: info
         end subroutine
      end interface

      n = size(a, dim=1)
      nrhs = 1
      lda = n
      ldb = n

      a_lapack = a
      b_lapack = b
      call dgesv(n, nrhs, a_lapack, lda, ipiv, b_lapack, ldb, info)

      x = b_lapack
   end function solve_system

   function intersect_two_lines(l1_x, l1_y, l2_x, l2_y) result(intersections)
      real(pr), intent(in) :: l1_x(:), l1_y(:), l2_x(:), l2_y(:)
      type(point), allocatable :: intersections(:)

      real(pr) :: s, t
      integer :: i, j

      real(pr) :: x, y, xold, yold
      xold = 9999
      yold = 9999

      allocate (intersections(0))

      line1: do i = 2, size(l1_x)
         line2: do j = 2, size(l2_x)
            associate ( &
               x1 => l1_x(i - 1), x2 => l1_x(i), &
               x3 => l2_x(j - 1), x4 => l2_x(j), &
               y1 => l1_y(i - 1), y2 => l1_y(i), &
               y3 => l2_y(j - 1), y4 => l2_y(j))
               call intersects(x1, x2, x3, x4, y1, y2, y3, y4, s, t)

               if (0 <= s .and. s <= 1 .and. 0 <= t .and. t <= 1) then
                  x = s*(x2 - x1) + x1
                  y = s*(y2 - y1) + y1

                  if ( &
                     abs((x - xold)) > 1e-2_pr .and. abs((y - yold)) > 1e-2_pr &
                     ) then
                     xold = x
                     yold = y
                     intersections = [intersections, point(x, y, i, j)]
                     exit line2
                  end if

               end if
            end associate
         end do line2
      end do line1
      if (size(intersections) > 3) then
         deallocate(intersections)
         allocate(intersections(0))
      end if
   end function

   function intersect_one_line(lx, ly) result(intersections)
      real(pr), intent(in) :: lx(:), ly(:)
      type(point), allocatable :: intersections(:)

      real(pr) :: s, t
      integer :: i, j

      real(pr) :: x, y, xold, yold

      xold = 9999
      yold = 9999

      allocate (intersections(0))
      line1: do i = 2, size(lx) - 1
         line2: do j = i + 2, size(lx)
            associate ( &
               x1 => lx(i - 1), x2 => lx(i), &
               x3 => lx(j), x4 => lx(j - 1), &
               y1 => ly(i - 1), y2 => ly(i), &
               y3 => ly(j), y4 => ly(j - 1))

               call intersects(x1, x2, x3, x4, y1, y2, y3, y4, s, t)
               ! print *, "===", i, j,"==="
               ! print *, x1, y1, x2, y2
               ! print *, x3, y3, x4, y4
               ! print *, s, t
               ! print *, "---------------"
               if (0 <= s .and. s <= 1 .and. 0 <= t .and. t <= 1) then
                  x = s*(x2 - x1) + x1
                  y = s*(y2 - y1) + y1
                  if (abs(x - xold) > 1 .or. abs(y - yold) > 1) then
                     xold = x
                     yold = y
                     ! Use earliest point for the "other" line
                     intersections = [intersections, point(x, y, i, j - 1)]
                  end if
               end if
            end associate
         end do line2
      end do line1
      if (size(intersections) > 3) then
         deallocate(intersections)
         allocate(intersections(0))
      end if
   end function

   subroutine intersects(x1, x2, x3, x4, y1, y2, y3, y4, s, t)
      real(pr), intent(in) :: x1, x2, x3, x4, y1, y2, y3, y4
      real(pr), intent(out) :: s, t

      real(pr) :: A(2, 2), b(2), tmp

      A(1, :) = [x2 - x1, x3 - x4]
      A(2, :) = [y2 - y1, y3 - y4]
      b = [x3 - x1, y3 - y1]

      b = solve_system(a, b)
      s = b(1)
      t = b(2)
   end subroutine

   elemental function interpol(x1, x2, y1, y2, x_obj) result(y)
      !! Linear interpolation.
      !!
      !! Calculates the linear interpolation between two points at a desired
      !! x value with the equation:
      !! \[ 
      !!    y = \frac{y_2 - y_1}{x_2 - x_1} \cdot (x_{obj})  - x_1 + y_1
      !! \]
      !! 
      !! Since this function is defined as `elemental` it will also interpolate
      !! a set of vectors.
      !!
      !! Examples of usage:
      !!
      !! ```fortran
      !! x1 = 2
      !! x2 = 5
      !! y1 = 2
      !! y2 = 9
      !! y = interpol(x1, x2, y1, y2, 2.3)
      !! ```
      !!
      !! ```fortran
      !! x1 = 2
      !! x2 = 5
      !! y1 = [2, 6]
      !! y2 = [9, 15]
      !! y = interpol(x1, x2, y1, y2, 2.3)
      !! ```
      real(pr), intent(in) :: x1 !! First point x value
      real(pr), intent(in) :: x2 !! Second point x value
      real(pr), intent(in) :: y1 !! First point y value 
      real(pr), intent(in) :: y2 !! Second point y value
      real(pr), intent(in) :: x_obj !! Desired x value to interpolate
      real(pr) :: y !! y value at `x_obj`
      y = (y2 - y1)/(x2 - x1)*(x_obj - x1) + y1
   end function

   pure logical function allclose_outer(array, array_test, atol)
      real(pr), intent(in) :: array(:)
      real(pr), intent(in) :: array_test(:)
      real(pr), intent(in) :: atol

      allclose_outer = all(abs(array - array_test) < atol)
   end function
   
   pure logical function allclose_inner(array, atol)
      real(pr), intent(in) :: array(:)
      real(pr), intent(in) :: atol

      integer :: i, j
      allclose_inner = .true.

   end function

   subroutine full_newton(fun, iters, X, ns, S, max_iters, F, dF, solvetol)
      !! Subroutine to solve a point in the envelope.
      !!
      !! Procedure that solves a point with the Newton-Raphson method.
      use stdlib_optval, only: optval
      use constants, only: ouput_path
      ! use minpack_module, only: dpmpar, hybrj1
      interface
         subroutine fun(X, ns, S, F, dF)
            !! Function to solve
            import pr
            real(pr), intent(in) :: X(:)
            integer, intent(in) :: ns
            real(pr), intent(in) :: S
            real(pr), intent(out) :: F(size(X))
            real(pr), intent(out) :: dF(size(X), size(X))
         end subroutine
      end interface
      !&<
      integer,  intent(out)    :: iters !! Number of iterations needed
      real(pr), intent(in out) :: X(:)  !! Variables vector
      integer,  intent(in)     :: ns    !! Number of specification
      real(pr), intent(in)     :: S     !! Specification value
      integer, intent(in)      :: max_iters !! Maximum iterations
      real(pr), intent(out)    :: F(size(X)) !! Function values at solved point
      real(pr), intent(out)    :: df(size(X), size(X)) !! Jacobian values
      real(pr), optional, intent(in) :: solvetol
      !&>

      real(pr) :: b(size(X)), A(size(X), size(X))
      real(pr) :: dX(size(X)), tol

      integer :: n, ldfjac, info, lwa, funit
      real(pr) :: fvec(size(x)), fjac(size(x),size(x))
      real(pr), allocatable :: wa(:)

      n = size(X)

      tol = optval(solvetol, 1.e-5_pr)
      ! ldfjac = n
      ! lwa = (n*(n+13))/2
      ! allocate(wa(lwa))
      ! call hybrj1(fcn, n, x, Fvec, Fjac, Ldfjac, Tol, Info, Wa, Lwa)

      ! f = fvec
      ! df = fjac
      ! iters = 3
      ! if (info == 4) iters = max_iters+1

      dX = 20
      b = 500

      open(newunit=funit, file="newton")

      newton: do iters = 1, max_iters
         if (maxval(abs(dx)) < tol .or. maxval(abs(b)) < tol) exit newton
         call fun(X, ns, S, b, a)
         b = -b

         if (any(isnan(b))) dx = dx/maxval(abs(dx))
         dX = solve_system(A, b)

         do while(maxval(abs(dx/x)) > 1)
            dX = dX/2
         end do

         write(funit, *) iters, X, dX

         X = X + dX
      end do newton
      close(funit)

      F = -b
      dF = A

      contains
         subroutine fcn(n, x, fvec, fjac, ldfjac, iflag)
            integer, intent(in) :: N
            real(pr), intent(in) :: x(n)
            real(pr), intent(inout) :: fvec(n)
            real(pr), intent(inout) :: fjac(ldfjac, n)
            integer, intent(in) :: ldfjac
            integer, intent(inout) :: iflag
            call fun(X, ns, S, fvec, fjac)
         end subroutine
   end subroutine
end module linalg


! module polyfit_mod
! 
!   implicit none
!       interface 
!     subroutine sgelss (m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
!       integer :: m
! 		integer :: n
! 		integer :: nrhs
! 		real, dimension( lda, * ) :: a
! 		integer :: lda
! 		real, dimension( ldb, * ) :: b
! 		integer :: ldb
! 		real, dimension( * ) :: s
! 		real :: rcond
! 		integer :: rank
! 		real, dimension( * ) :: work
! 		integer :: lwork
! 		integer :: info 
!     end subroutine
!    end interface
! 
! contains
! 
!   function vander(x,N,increasing) result(v)
!     real, intent(in) :: x(:)
!     integer, intent(in) :: N
!     logical, intent(in), optional :: increasing
!     real :: v(size(x),N)
!     
!     integer :: i
!     logical :: increasing_
! 
!     increasing_ = .false.
!     if (present(increasing)) increasing_ = increasing    
! 
!     if (increasing_) then
!       if (N > 0) v(:,1) = 1
!       if (N > 1) then
!         v(:,2:) = spread(x,dim=2,ncopies=N-1)
!         do i = 2, N-1
!           v(:,i+1) = v(:,i)*v(:,i+1)
!         end do
!       end if
!     else
!       if (N > 0) v(:,N) = 1
!       if (N > 1) then
!         v(:,1:N-1) = spread(x,dim=2,ncopies=N-1)
!         do i = N-1,2,-1
!           v(:,i-1) = v(:,i)*v(:,i-1)
!         end do
!       end if
!     end if
!   end function
! 
!   function polyfit(x,y,deg,rcond,rank,singular_values) result(p)
!     real, intent(in) :: x(:), y(:)
!     integer, intent(in) :: deg
!     real, intent(in), optional :: rcond
!     integer, intent(out), optional :: rank
!     real, intent(out), allocatable, optional :: singular_values(:)
!     real :: p(deg + 1)
! 
!     integer :: order
!     real :: rcond_
!     real, allocatable :: lhs(:,:), rhs(:), scl(:),s(:), work(:)
!     integer :: m, n, rank_, info, lwork
! 
!     order = deg + 1
! 
!     ! check parameters
!     if (deg < 0) error stop "expected deg >= 0"
!     if (size(x) == 0) error stop "expected non-empty vector for x"
!     if (size(x) /= size(y)) error stop "expected x and y to have same length"
! 
!     ! set rcond
!     rcond_ = size(x)*epsilon(x)
!     if (present(rcond)) rcond_ = rcond
! 
!     m = size(x)
!     n = order
!     allocate(lhs(m,n),scl(m))
! 
!     ! set up least squares equation for powers of x  
!     lhs = vander(x,order)
! 
!     call print_mat(lhs)
!     ! scale lhs to improve the condition number
!     scl = sqrt(sum(lhs*lhs,dim=1))
!     print *, "scl = "
!     call print_mat(spread(scl,dim=1,ncopies=m))
!     lhs = lhs/spread(scl,dim=1,ncopies=m)
!     write(*,*)
!     call print_mat(lhs)
! 
!     allocate(rhs,source=y)
! 
!     !
!     ! solve lstsq problem
!     !
!     allocate(s(min(m,n)))
! 
!     ! perform workspace query
!     lwork = -1
!     allocate(work(1))
!     call sgelss(m,n,1,lhs,m,rhs,m,s,rcond_,rank_,work,lwork,info)
!     print *, "info = ", info
! 
!     ! resize workspace
!     lwork = int(work(1))
!     print *, "lwork = ", lwork, work(1)
!     deallocate(work)
!     allocate(work(lwork))
! 
!     print *, rhs
!     ! compute minimum norm solutiong using SVD
!     call sgelss(m,n,1,lhs,m,rhs,m,s,rcond_,rank_,work,lwork,info)
!     if (info /= 0) then
!       if (info < 0) write(*,*) 'sgelss: the ',abs(info),'-th argument had an illegal value'
!       if (info > 0) write(*,*) 'sgelss: the algorithm for computing the SVD failed to converge'
!       error stop
!     end if
!     call print_mat(lhs)
!     print *, rhs
! 
!     if (rank_ /= order) then
!       write(*,*) "polyfit may be poorly conditioned"
!     end if
! 
!     ! broadcast scale coefficients
!     p = rhs/scl
! 
!     if (present(rank)) rank = rank_
!     if (present(singular_values)) then
!       if (allocated(singular_values)) deallocate(singular_values)
!       allocate(singular_values,source=s)
!     end if
!   end function
! 
!   subroutine print_mat(v)
!     real, intent(in) :: v(:,:)
!     integer :: i, j
!     do i = 1, size(v,1)
!       print *, (v(i,j),j=1,size(v,2))
!     end do
!   end subroutine
! 
! end module
! 
! program main
! 
!   use polyfit_mod
!   implicit none
! 
!   real :: x(4), y(4)
!   real :: p(3)
! 
!   x = [1.,2.,3.,4.]
!   y = x**2 + 2*x - 4
! 
!   p = polyfit(x,y,2)
!   print *, p
! 
! end program
! 
