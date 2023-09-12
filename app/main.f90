program main
   use dtypes, only: envelope
   use constants, only: pr, ouput_path
   use legacy_ar_models, only: nc
   use flap, only: command_line_interface

   implicit none
   real(pr) :: et, st

   type(command_line_interface) :: cli
   integer :: cli_error
   character(len=99) :: cli_string
   
   type(envelope) :: bub_env, dew_env

   call cli%init(progname="envelopes", description="Phase Envelopes")
   call cli%add(&
            switch="--infile", &
            switch_ab="-i", &
            help="Input file", &
            error=cli_error, &
            required=.true.)
   call cli%parse(error=cli_error)

   if (cli_error /= 0) stop

   call system("mkdir -p " // trim(ouput_path))
   call system("rm " // trim(ouput_path) // "*")
   
   call setup ! Setup module variables 

   call cpu_time(st)
   call pt_envelopes ! Calculate PT envelopes at the system's composition
   call cpu_time(et)
   print *, "PT: ", (et-st) * 1000 , "ms"

   call cpu_time(st)
   call px_envelopes ! Calculate Px envelopes
   call cpu_time(et)
   print *, "PX: ", (et-st) * 1000 , "ms"
contains
   subroutine setup
      use io_nml, only: read_system, write_system
      use inj_envelopes, only: setup_inj => from_nml
      integer :: funit_system
      character(len=500) :: infile

      call cli%get(val=infile, switch="--infile", error=cli_error)

      call read_system(trim(infile))
      call setup_inj(trim(infile))

      open (newunit=funit_system, file="systemdata.nml")
      call write_system(funit_system)
      close (funit_system)
   end subroutine

   subroutine pt_envelopes
      use legacy_ar_models, only: z
      use envelopes, only: envelope2, max_points, k_wilson_bubble, &
                        max_points, p_wilson, k_wilson
      !! Calculation of PT envelopes of the main system.
      real(pr), allocatable :: tv(:) ! Temperatures [K]
      real(pr), allocatable :: pv(:) ! Pressures [bar]
      real(pr), allocatable :: dv(:) ! Densities [mol/L]

      real(pr) :: tcri(4)            ! Critical points temperatures
      real(pr) :: pcri(4)            ! Critical points pressures
      real(pr) :: dcri(4)            ! Critical points densities

      real(pr) :: t, p               ! Temperature and pressure
      real(pr), allocatable :: k(:)  ! K factors
      integer :: n_points, icri(4), ncri, i

      integer :: n

      allocate (tv(max_points), pv(max_points), dv(max_points))
      allocate (k(size(z)))

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
      
      ! ======================================================================
      !  Look for crossings
      ! ----------------------------------------------------------------------
      check_crossings: block
         use linalg, only: point, intersection
         type(point), allocatable :: inter(:)
         inter = intersection( &
                 dew_env%t, dew_env%p, &
                 bub_env%t, bub_env%p &
               )
         print *, "Intersections: ", size(inter)
         do i = 1, size(inter)
            print *, inter(i)
         end do
      end block check_crossings
      ! ======================================================================
   end subroutine

   subroutine px_envelopes
      !! Calculation of Px envelopes at selected temperature.
      use inj_envelopes, only: F_injection, full_newton, z_injection, &
                               T_inj => T, injection_envelope, z_0, injection_case, &
                               injelope
      use envelopes, only: envelope, k_wilson, p_wilson
      use linalg, only: interpol

      real(pr), allocatable :: X(:), F(:), F2(:), dF(:, :), df_num(:, :)
      real(pr) :: alpha, S
      integer :: ns, i, iters, idx, ti
      integer, allocatable :: i_inj(:)
      real(pr), allocatable :: ts_envel(:), ts(:)
      real(pr) :: t, p
      real(pr), allocatable :: k(:)
      real(pr) :: t_tol = 2

      type(injelope) :: bub_envels, dew_envels
      allocate (X(nc + 2), F(nc + 2), dF(nc + 2, nc + 2), df_num(nc + 2, nc + 2), F2(nc + 2))

      ! ======================================================================
      !  Setup system
      ! ----------------------------------------------------------------------
      alpha = 0.0
      z_injection = z_injection/sum(z_injection)
      ns = nc + 2
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
            X(nc+1) = log(P)
            X(nc+2) = alpha

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
   end subroutine
end program main
