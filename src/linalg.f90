module linalg
   !! Wrapper module around LAPACK's `dgesv`
   use constants, only: pr
   implicit none
contains
   function solve_system(a, b) result(x)
      real(pr), intent(in out) :: b(:)
      real(pr), intent(in out) :: a(size(b), size(b))

      real(pr) :: x(size(b))

      real(8) :: a_lapack(size(b), size(b)), b_lapack(size(b))

      integer :: n, nrhs, lda, ipiv(size(b)), ldb, info

      n = size(a, dim=1)
      nrhs = 1
      lda = n
      ldb = n

      a_lapack = a
      b_lapack = b
      call dgesv(n, nrhs, a_lapack, lda, ipiv, b_lapack, ldb, info)

      x = b_lapack
   end function solve_system
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