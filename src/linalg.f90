module linalg
   !! Wrapper module around LAPACK's `dgesv`
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
contains
   function solve_system(a, b) result(x)
      real(pr), intent(in out) :: b(:)
      real(pr), intent(in out) :: a(size(b), size(b))

      real(pr) :: x(size(b))

      real(8) :: a_lapack(size(b), size(b)), b_lapack(size(b))

      integer :: n, nrhs, lda, ipiv(size(b)), ldb, info

      interface
         subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            integer :: n
            integer :: nrhs
            real(8) :: a(n,n)
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

      real(pr) :: x, y, xold=9999, yold=9999

      allocate(intersections(0))

      line1: do i=2, size(l1_x)
         line2: do j=2, size(l2_x)
            associate(x1 => l1_x(i-1), x2 => l1_x(i), x3 => l2_x(j-1), x4 => l2_x(j), &
                      y1 => l1_y(i-1), y2 => l1_y(i), y3 => l2_y(j-1), y4 => l2_y(j))
               call intersects(x1, x2, x3, x4, y1, y2, y3, y4, s, t)

               if (0 <= s .and. s <= 1 .and. 0 <= t .and. t <= 1) then
                  x = s * (x2-x1) + x1
                  y = s * (y2-y1) + y1

                  if (abs(x - xold) > 1 .and. abs(y - yold) > 1) then
                     print *, "CROSS", x, y
                     xold = x
                     yold = y
                     intersections = [intersections, point(x, y, i, j)]
                     exit line2
                  end if

               end if
            end associate
         end do line2
      end do line1
   end function
   
   function intersect_one_line(lx, ly) result(intersections)
      real(pr), intent(in) :: lx(:), ly(:)
      type(point), allocatable :: intersections(:)
      character(len=*), parameter :: fmt="(*(G0,:,', '))"

      real(pr) :: s, t
      integer :: i, j

      real(pr) :: x, y, xold=9999, yold=9999

      allocate(intersections(0))
      line1: do i=2, size(lx)-1
         line2: do j=i+15, size(lx)
            associate(x1 => lx(i-1), x2 => lx(i), x3 => lx(j), x4 => lx(j-1), &
                      y1 => ly(i-1), y2 => ly(i), y3 => ly(j), y4 => ly(j-1))
            call intersects(x1, x2, x3, x4, y1, y2, y3, y4, s, t)
            if (0 <= s .and. s <= 1 .and. 0 <= t .and. t <= 1) then
            
               x = s*(x2 - x1) + x1
               y = s*(y2 - y1) + y1

            if (abs(x - xold) > 1 .and. abs(y - yold) > 1) then
               print *, "CROSS"
               print *, i, j, x, y
               write(*, fmt) "x1, y1 = ", x1, y1
               write(*, fmt) "x2, y2 = ", x2, y2
               write(*, fmt) "x3, y3 = ", x3, y3
               write(*, fmt) "x4, y4 = ", x4, y4
               xold = x
               yold = y
               intersections = [intersections, point(x, y, i, j)]
            end if
            end if
            end associate
         end do line2
      end do line1
   end function

   subroutine intersects(x1, x2, x3, x4, y1, y2, y3, y4, s, t)
      real(pr), intent(in) :: x1, x2, x3, x4, y1, y2, y3, y4
      real(pr), intent(out) :: s, t

      real(pr) :: A(2,2), b(2), tmp

      A(1, :) = [x2-x1, x3-x4]
      A(2, :) = [y2-y1, y3-y4]
      b = [x3-x1, y3-y1]

      b = solve_system(a, b)
      s = b(1)
      t = b(2)
   end subroutine

   elemental function interpol(x1, x2, y1, y2, x_obj) result(y)
      real(pr), intent(in) :: x1
      real(pr), intent(in) :: x2
      real(pr), intent(in) :: y1
      real(pr), intent(in) :: y2
      real(pr), intent(in) :: x_obj
      real(pr) :: y
      y = (y2 - y1)/(x2 - x1) * (x_obj - x1) + y1
   end function
end module linalg