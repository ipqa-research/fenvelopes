module dtypes
   use constants, only: pr
   use io

   implicit none

   private
   public :: AbsEnvel
   public :: envelope
   public :: critical_point

   type, abstract :: AbsEnvel
   end type
   
   type :: critical_point
      real(pr) :: t
      real(pr) :: p
      real(pr) :: alpha
   end type critical_point

   type :: envelope
      real(pr), allocatable :: vars(:, :)  !! Value of the set of variables at each point
      real(pr), allocatable :: z(:) !! Global composition
      real(pr), allocatable :: t(:) !! Temperature points
      real(pr), allocatable :: p(:) !! Pressure points
      real(pr), allocatable :: logk(:, :) !! ln(K) for each point
      real(pr), allocatable :: logphi(:, :) !! lnphi for each point
      type(critical_point), allocatable :: critical_points(:) !! Critical points
   end type envelope
end module