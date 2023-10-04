module dtypes
   use constants, only: pr
   use io

   implicit none

   private
   public :: AbsEnvel
   public :: envelope
   public :: env3
   public :: print_header
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
   contains
      procedure :: write => write_envel_2
   end type envelope

   type, extends(envelope) :: env3
      real(pr), allocatable :: beta(:) !! Other phase fraction
      real(pr), allocatable :: x(:, :) !! 
      real(pr), allocatable :: y(:, :) !! 
      real(pr), allocatable :: w(:, :) !! 
      real(pr), allocatable :: logks(:, :) !! ln(Ks)
      type(critical_point), allocatable :: ll_critical_points(:)
   contains
      procedure :: write => write_envel_3
   end type env3
contains
   subroutine write_critical_points(self, file_name)
      type(critical_point), intent(in) :: self(:)
      character(len=*), optional, intent(in) :: file_name !! Ouptut file name

      character(len=:), allocatable :: filename
      integer :: file_unit
      integer :: i
      
      if (present(file_name)) then
         filename = file_name
      else 
         filename = "CP"
      end if
      
      open(newunit=file_unit, file=filename)
         write(file_unit, "(A)") "P            T"
         do i = 1, size(self)
            write(file_unit, "(2(E10.5,2x))") self(i)%t,  self(i)%p
         end do
      close(file_unit)
   end subroutine

   subroutine write_envel_2(self, file_name)
      class(envelope), intent(in):: self
      character(len=*), optional, intent(in) :: file_name !! Ouptut file name
      character(len=:), allocatable :: filename
      integer :: i, n, file_unit, n_components

      if (present(file_name)) then
         filename = file_name
      else 
         filename = "envelout-2phase"
      end if

      n = size(self%t)
      n_components = size(self%z)

      associate(t => self%t, p => self%p, &
                logk => self%logk,        &
                z => self%z &
          )
         open(newunit=file_unit, file=filename)
            write(file_unit, *) &
             "P ", "T ", ("K" // str(i), i=1,n_components), &
             ("z" // str(i), i=1,n_components)
            do i=1,n
               write(file_unit, *) p(i), t(i), logk(i, :), z
            end do
         close(file_unit)
      end associate
      
      ! Write Critical Points file
      ! filename = filename // "-CP"

      ! call write_critical_points(self%critical_points, filename)

      deallocate(filename)
   end subroutine write_envel_2

   subroutine write_envel_3(self, file_name)
      class(env3), intent(in):: self
      character(len=*), optional, intent(in) :: file_name
      character(len=:), allocatable :: filename
      integer :: i, n, file_unit, n_components

      if (present(file_name)) then
         filename = file_name
      else 
         filename = "envelout-3phase"
      end if

      n_components = size(self%z)

      associate(t => self%t, p => self%p, &
          logk => self%logk, logks => self%logks, &
          x => self%x, y => self%y, w => self%w, beta => self%beta  &
          )
      n = size(self%t)

      open(newunit=file_unit, file=filename)
      write(file_unit, *) &
          "P ", "T ", "beta ", &
         ("K" // str(i), i=1,n_components), &
         ("KS" // str(i), i=1,n_components), &
         ("x" // str(i), i=1,n_components), &
         ("y" // str(i), i=1,n_components), &
         ("w" // str(i), i=1,n_components)
      do i=1,n
         write(file_unit, *) p(i), t(i), beta(i), logk(i, :), logks(i, :), x(i, :), y(i, :), w(i, :)
      end do
      close(file_unit)
      end associate
      
      ! Write Critical Points file
      filename = filename // "-CP"
      call write_critical_points(self%critical_points, filename)
      deallocate(filename)
   end subroutine write_envel_3

   subroutine print_header(name)
      character(len=250), intent(in) :: name

      print *, "==================================="
      print *, "!", name
      print *, "-----------------------------------"
   end subroutine print_header
end module dtypes