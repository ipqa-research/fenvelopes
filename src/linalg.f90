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