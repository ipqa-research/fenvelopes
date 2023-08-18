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

   subroutine px_envelopes
      !! Calculation of Px envelopes at selected temperature.
      use inj_envelopes, only: F_injection, full_newton, z_injection, &
                               T_inj => T, injection_envelope, z_0, injection_case, &
                               injelope, funit_output
      use envelopes, only: envelope, k_wilson, p_wilson
      use linalg, only: interpol

      real(pr), allocatable :: X(:), F(:), F2(:), dF(:, :), df_num(:, :)
      real(pr) :: alpha, S
      integer :: ns, i, iters, idx, ti
      integer, allocatable :: i_inj(:)
      real(pr), allocatable :: ts_envel(:), ts(:)
      real(pr) :: t_tol = 2

      type(injelope) :: bub_envels, dew_envels
      allocate (X(nc + 2), F(nc + 2), dF(nc + 2, nc + 2), df_num(nc + 2, nc + 2), F2(nc + 2))

      ! ======================================================================
      !  Setup system
      ! ----------------------------------------------------------------------
      alpha = 0.0
      z_injection = z_injection/sum(z_injection)
      ns = nc + 2
      open (newunit=funit_output, file="px.dat")
      ! ======================================================================

      ! ======================================================================
      !  Bubble envelope
      ! ----------------------------------------------------------------------
      print *, "Running Bubble"
      bubble: block
         real(pr) :: pold
         pold = 0
         ts_envel = pack(bub_env%t, mask=abs(bub_env%t - t_inj) < t_tol)
         do i = 1, size(ts_envel)
            idx = findloc(bub_env%t, value=ts_envel(i), dim=1)
            print *, ts_envel(i)

            p = interpol( &
                bub_env%t(idx), bub_env%t(idx + 1), &
                bub_env%p(idx), bub_env%p(idx + 1), &
                t_inj)

            if (abs(p - pold) < 5) cycle
            pold = p

            k = exp(interpol( &
                    bub_env%t(idx), bub_env%t(idx + 1), &
                    bub_env%logk(idx, :), bub_env%logk(idx + 1, :), &
                    t_inj))

            X(1:nc) = log(K)
            X(nc + 1) = log(P)
            X(nc + 2) = alpha

            call injection_envelope(X, ns, 0.01_pr, bub_envels)
         end do
      end block bubble
      ! ======================================================================

      ! ======================================================================
      !  Dew envelope
      ! ----------------------------------------------------------------------
      print *, "Running Dew"
      dew: block
         real(pr) :: pold
         pold = 0
         ts_envel = pack(dew_env%t, mask=abs(dew_env%t - t_inj) < t_tol)
         do i = 1, size(ts_envel)
            idx = findloc(dew_env%t, value=ts_envel(i), dim=1)

            alpha = 0

            p = interpol( &
                dew_env%t(idx), dew_env%t(idx + 1), &
                dew_env%p(idx), dew_env%p(idx + 1), &
                t_inj)

            if (abs(p - pold) < 5) cycle
            pold = p
            print *, ts_envel(i), p

            k = exp(interpol( &
                    dew_env%t(idx), dew_env%t(idx + 1), &
                    dew_env%logk(idx, :), dew_env%logk(idx + 1, :), &
                    t_inj))

            X(1:nc) = log(K)
            X(nc + 1) = log(P)
            X(nc + 2) = alpha

            call injection_envelope(X, ns, 0.01_pr, dew_envels)
         end do
      end block dew
      ! ======================================================================

      ! ======================================================================
      !  Look for crossings
      ! ----------------------------------------------------------------------
      check_crossings: block
         use linalg, only: point, intersection
         type(point), allocatable :: inter(:)
         inter = intersection( &
                 dew_envels%alpha, dew_envels%p, &
                 bub_envels%alpha, bub_envels%p &
                 )

         do i = 1, size(inter)
            print *, inter(i)
         end do
      end block check_crossings
      ! ======================================================================
      close (funit_output)
      print *, "END"
   end subroutine
end program main
