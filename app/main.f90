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

   subroutine pt_envelopes
      !! Calculation of PT envelopes of the main system.
      
      integer :: n

      ! =====================================================================
      !  Bubble envel
      ! ---------------------------------------------------------------------
      call k_wilson_bubble(z, t, p, k)
      call envelope2( &
         1, nc, z, T, P, k, &
         n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, &
         bub_env &
         )
      ! =====================================================================

      ! =====================================================================
      !  Dew/AOP envelopes
      ! ---------------------------------------------------------------------
      t = 315
      p = p_wilson(z, t)
      do while (p > 0.1)
         t = t - 5
         p = p_wilson(z, t)
      end do
      k = 1/k_wilson(t, p)

      call envelope2( &
         2, nc, z, T, P, k, &
         n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, &
         dew_env &
         )

      ! Remove the low pressure parts.
      n = 1
      do i = 2, size(dew_env%t)
         n = n + 1
         if (dew_env%t(i) - dew_env%t(i - 1) < 0) exit
      end do

      if (n /= size(dew_env%t)) then
         dew_env%t = dew_env%t(i:)
         dew_env%p = dew_env%p(i:)
         dew_env%logk = dew_env%logk(i:, :)
      end if
      ! =====================================================================
   end subroutine
end program main
