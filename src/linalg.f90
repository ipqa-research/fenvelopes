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
end module linalg