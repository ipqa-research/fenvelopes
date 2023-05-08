module dtypes
   use constants, only: pr
   use io

   implicit none

   private
   public :: envelope
   public :: env3
   public :: point
   public :: kfcross
   public :: print_header
   public :: find_cross
   public :: find_self_cross
   public :: critical_point

   
   type :: critical_point
      real(pr) :: t
      real(pr) :: p
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

   type :: point
      real(pr) :: x
      real(pr) :: y
      integer :: i
      integer :: j
   end type point

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
      filename = filename // "-CP"

      call write_critical_points(self%critical_points, filename)

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

   function kfcross(i, t_values, logk, target_t)
      !! Estimate the Kvalues of an envelope by interpolation from a near point.
      integer, intent(in) :: i
      real(pr), allocatable, intent(in) :: t_values(:) !! Envelope's temperatures
      real(pr), allocatable, intent(in) :: logk(:, :) !! Envelope's kvalues
      real(pr), intent(in) :: target_t !! Target temperature where to interpolate

      real(pr), allocatable :: kfcross(:) !! Kvalues at desired point

      kfcross = (logk(i, :) - logk(i - 1, :)) &
                /&
                (t_values(i) - t_values(i - 1)) &
                * (target_t - t_values(i - 1)) &
                + logk(i - 1, :)
   end function

   subroutine print_header(name)
      character(len=250), intent(in) :: name

      print *, "==================================="
      print *, "!", name
      print *, "-----------------------------------"
   end subroutine print_header
   
   subroutine find_cross(tv1, tv2, pv1, pv2, crossings, crossed)
      !! Find crossings between two given lines
      !!
      !! Returns an array of crossigns, containings the crosses found. Each row
      !! contains the data from each found cross
      !!
      !!  | --------| ------- | ---------------- | ----------------- |
      !!  | x_cross | y_cross | first_line_index | second_line_index |
      !!  | --------| ------- | ---------------- | ----------------- |
      !!

      real(pr), intent(in)  :: tv1(:)  !! First line x values
      real(pr), intent(in)  :: tv2(:)  !! Second line x values
      real(pr), intent(in)  :: pv1(:)  !! First line y values
      real(pr), intent(in)  :: pv2(:)  !! Second line y values
      logical, optional, intent(out) :: crossed

      type(point), allocatable :: crossings(:) !! Array of crossings
      type(point) :: current_cross

      real(pr) :: x11, x12, x21, x22, y11, y12, y21, y22

      real(pr) :: x_cross, y_cross, m1, b1, m2, b2, xlow, xup, ylow, yup
      real(pr), dimension(2) :: xpair_1, xpair_2, ypair_1, ypair_2
      integer :: i, j, n

      if (present(crossed)) then
         crossed = .false.
      end if

      if (allocated(crossings)) then
         deallocate (crossings)
      end if

      allocate (crossings(0))
      n = 0

      do i = 2, size(tv1)
         xpair_1 = tv1(i - 1:i)
         ypair_1 = pv1(i - 1:i)

         x11 = xpair_1(1)
         x12 = xpair_1(2)
         y11 = ypair_1(1)
         y12 = ypair_1(2)

         m1 = (y12 - y11)/(x12 - x11)
         b1 = y11 - m1*x11

         do j = 2, size(tv2)
            xpair_2 = tv2(j - 1:j)
            ypair_2 = pv2(j - 1:j)

            x21 = xpair_2(1)
            x22 = xpair_2(2)
            y21 = ypair_2(1)
            y22 = ypair_2(2)

            m2 = (y22 - y21)/(x22 - x21)
            b2 = y21 - m2*x21

            x_cross = (b1 - b2)/(m2 - m1)
            y_cross = m1*x_cross + b1

            xlow = max(minval(xpair_1), minval(xpair_2))
            xup = min(maxval(xpair_1), maxval(xpair_2))
            ylow = max(minval(ypair_1), minval(ypair_2))
            yup = min(maxval(ypair_1), maxval(ypair_2))

            if ( &
               (xlow <= x_cross) .and. (x_cross <= xup) .and. &
               (ylow <= y_cross) .and. (y_cross <= yup) &
               ) then
               if (present(crossed)) crossed = .true.
               print *, "CROSS:", i, j, x_cross, y_cross

               ! TODO: This should get back, but for some reason now
               ! there is a dimension 0 error that didn't appear before

               ! if ((abs(x_cross - crossings(n)%x) < 0.1) .and. &
               !     (abs(y_cross - crossings(n)%y) < 0.1)) then
               !    print *, "CROSS: Repeated cross, skipping..."
               !    cycle
               ! end if

               current_cross = point(x_cross, y_cross, i, j)
               n = n + 1
               crossings = [crossings, current_cross]

            end if
         end do
      end do
   end subroutine find_cross

   subroutine find_self_cross(array_x, array_y, found_cross, crossed)
      use constants, only: pr
      use array_operations, only: diff, mask

      real(pr), intent(in) :: array_x(:)
      real(pr), intent(in) :: array_y(size(array_x))
      type(point), allocatable, intent(in out) :: found_cross(:)
      logical, optional, intent(out) :: crossed

      logical, allocatable :: filter(:)
      integer, allocatable :: msk(:)
      real(pr) :: min_x, max_x

      integer :: i, idx, idy

      if(present(crossed)) crossed = .false.

      ! All the values with positive delta 
      filter = diff(array_x) > 0

      return

      i = 1
      do while(filter(i))
         ! Find the first ocurrence of a negative delta x
         ! This will give the index of the cricondentherm
         i = i + 1
      end do

      ! if (i < size(array_x)) then
      !    msk = mask(filter(i:)) + i
      !    max_x = maxval(array_x(msk))
      !    min_x = minval(array_x(msk))

      !    ! 
      !    filter = array_x <= max_x - 5 .and. array_x >= min_x - 5 .and. array_y >= 10
      !    msk = mask(filter)

      !    call find_cross(&
      !       array_x(msk), array_x(msk), array_y(msk), array_y(msk), found_cross, crossed &
      !    )

      !    if (size(found_cross) > 1) then
      !       found_cross%i = found_cross%i + msk(1)
      !       found_cross%j = found_cross%j + msk(1)
      !    end if
      ! end if


      ! if (size(found_cross) > 0) then
      !    do i=1,size(found_cross)
      !       ! TODO: This assumes there is only one self-cross, should be better defined
      !       idx = minloc(abs(array_x - found_cross(i)%x), dim=1)
      !       idy = minloc(abs(array_y - found_cross(i)%y), dim=1)

      !       found_cross(i)%i = idx
      !       found_cross(i)%j = idy
      !    end do
      ! end if
   end subroutine find_self_cross
end module dtypes
