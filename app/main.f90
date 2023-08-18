program main
   use envelopes, only: envelope2, max_points, k_wilson_bubble, &
                        max_points, p_wilson, k_wilson
   use dtypes, only: envelope
   use constants, only: pr
   use system, only: z, nc

   implicit none

   real(pr), allocatable :: tv(:) ! Temperatures [K]
   real(pr), allocatable :: pv(:) ! Pressures [bar]
   real(pr), allocatable :: dv(:) ! Densities [mol/L]

   real(pr) :: tcri(4)            ! Critical points temperatures
   real(pr) :: pcri(4)            ! Critical points pressures
   real(pr) :: dcri(4)            ! Critical points densities

   real(pr) :: t, p               ! Temperature and pressure
   real(pr), allocatable :: k(:)  ! K factors

   integer :: n_points, icri(4), ncri, i

   type(envelope) :: bub_env, dew_env

   call setup              !
   call pt_envelopes
   call px_envelopes
contains
   subroutine setup
      use io_nml, only: read_system, write_system
      use system, only: kij
      use inj_envelopes, only: setup_inj => from_nml
      integer :: funit_system
      character(len=254) :: infile
      call get_command_argument(1, value=infile)

      call read_system(trim(infile))
      call setup_inj(trim(infile))

      open (newunit=funit_system, file="systemdata.nml")
      call write_system(funit_system)
      close (funit_system)

      allocate (tv(max_points), pv(max_points), dv(max_points))
      allocate (k(size(z)))
   end subroutine
end program main
