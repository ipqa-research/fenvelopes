module constants
   use iso_fortran_env, only: real32, real64, real128
   implicit none

   integer, parameter :: pr = real64 !! Precision
   real(pr), parameter :: R = 0.08314472 !! R constant []
   character(len=254) :: ouput_path = "fenvelopes_output/" !! Output path
   character(len=1) :: path_sep = "/" !! Path separator
end module constants
