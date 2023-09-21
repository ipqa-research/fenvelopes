program main
   use dtypes, only: envelope
   use inj_envelopes, only: injelope
   use constants, only: pr, ouput_path
   use legacy_ar_models, only: nc
   use flap, only: command_line_interface

   implicit none
   real(pr) :: et, st

   type(command_line_interface) :: cli
   integer :: cli_error
   character(len=99) :: cli_string

   type(envelope) :: pt_bub, pt_dew
   type(injelope) :: px_bub, px_dew

   call cli%init(progname="envelopes", description="Phase Envelopes")
   call cli%add( &
      switch="--infile", &
      switch_ab="-i", &
      help="Input file", &
      error=cli_error, &
      required=.true.)
   call cli%parse(error=cli_error)

   if (cli_error /= 0) stop

   call system("mkdir -p "//trim(ouput_path))
   call system("rm "//trim(ouput_path)//"*")

   call setup ! Setup module variables

   call cpu_time(st)
   call pt_envelopes ! Calculate PT envelopes at the system's composition
   call cpu_time(et)
   print *, "PT: ", (et - st)*1000, "ms"

   call cpu_time(st)
   call px_envelopes ! Calculate Px envelopes
   call cpu_time(et)
   print *, "PX: ", (et - st)*1000, "ms"
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

      ! ========================================================================
      !  Bubble envel
      ! ------------------------------------------------------------------------
      call k_wilson_bubble(z, t, p, k)
      call envelope2( &
         1, nc, z, T, P, k, &
         n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, &
         pt_bub &
         )
      ! ========================================================================

      ! ========================================================================
      !  Dew/AOP envelopes
      ! ------------------------------------------------------------------------
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
         pt_dew &
         )

      ! Remove the low pressure parts.
      n = 1
      do i = 2, size(pt_dew%t)
         n = n + 1
         if (pt_dew%t(i) - pt_dew%t(i - 1) < 0) exit
      end do

      if (n /= size(pt_dew%t)) then
         pt_dew%t = pt_dew%t(i:)
         pt_dew%p = pt_dew%p(i:)
         pt_dew%logk = pt_dew%logk(i:, :)
      end if
      ! ========================================================================

      ! ========================================================================
      !  Look for crossings
      ! ------------------------------------------------------------------------
      check_crossings: block
         use linalg, only: point, intersection
         type(point), allocatable :: inter(:)
         inter = intersection( &
                 pt_dew%t, pt_dew%p, &
                 pt_bub%t, pt_bub%p &
                 )
         print *, "Intersections: ", size(inter)
         do i = 1, size(inter)
            print *, inter(i)
         end do
      end block check_crossings
      ! ========================================================================
   end subroutine

   subroutine px_envelopes
      !! Calculation of Px envelopes at selected temperature.
      use inj_envelopes, only: full_newton, z_injection, &
                               T_inj => T, injection_envelope, z_0, &
                               injelope, injection_envelope_three_phase
      use envelopes, only: envelope, k_wilson, p_wilson
      use linalg, only: interpol

      real(pr), allocatable :: X(:)
      real(pr) :: alpha, del_S0
      integer :: ns, i, idx
      real(pr), allocatable :: ts_envel(:)
      real(pr) :: p
      real(pr), allocatable :: k(:)
      real(pr) :: t_tol = 2

      ! ========================================================================
      !  Setup system
      ! ------------------------------------------------------------------------
      alpha = 0.0_pr
      z_injection = z_injection/sum(z_injection)
      ns = nc + 2
      del_S0 = 0.1_pr
      ! ========================================================================

      ! ========================================================================
      !  Bubble envelope
      ! ------------------------------------------------------------------------
      print *, "Running Bubble"
      bubble: block
         real(pr) :: pold
         pold = 0
         ts_envel = pack(pt_bub%t, mask=abs(pt_bub%t - t_inj) < t_tol)
         do i = 1, size(ts_envel)
            idx = findloc(pt_bub%t, value=ts_envel(i), dim=1)
            p = interpol( &
                pt_bub%t(idx), pt_bub%t(idx + 1), &
                pt_bub%p(idx), pt_bub%p(idx + 1), &
                t_inj)

            if (abs(p - pold) < 5) cycle
            pold = p

            k = exp(interpol( &
                    pt_bub%t(idx), pt_bub%t(idx + 1), &
                    pt_bub%logk(idx, :), pt_bub%logk(idx + 1, :), &
                    t_inj))

            X = [log(K), log(P), alpha]

            call injection_envelope(X, ns, del_S0, px_bub)
         end do
      end block bubble
      ! ========================================================================

      ! ========================================================================
      !  Dew envelope
      ! ------------------------------------------------------------------------
      print *, "Running Dew"
      dew: block
         real(pr) :: pold
         pold = 0
         ts_envel = pack(pt_dew%t, mask=abs(pt_dew%t - t_inj) < t_tol)
         do i = 1, size(ts_envel)
            idx = findloc(pt_dew%t, value=ts_envel(i), dim=1)
            alpha = 0
            p = interpol( &
                pt_dew%t(idx), pt_dew%t(idx + 1), &
                pt_dew%p(idx), pt_dew%p(idx + 1), &
                t_inj)

            if (abs(p - pold) < 5) cycle
            pold = p
            k = exp(interpol( &
                    pt_dew%t(idx), pt_dew%t(idx + 1), &
                    pt_dew%logk(idx, :), pt_dew%logk(idx + 1, :), &
                    t_inj))

            X = [log(K), log(P), alpha]

            call injection_envelope(X, ns, del_s0, px_dew)
         end do
      end block dew
      ! ========================================================================

      ! ========================================================================
      !  Look for crossings
      ! ------------------------------------------------------------------------
      check_crossings: block
         use linalg, only: point, intersection
         type(point), allocatable :: inter(:), self_inter(:)
         inter = intersection( &
                 px_dew%alpha, px_dew%p, &
                 px_bub%alpha, px_bub%p)
         self_inter = intersection(px_dew%alpha, px_dew%p)

         print *, "Px Intersections:      ", size(inter)
         print *, "Px Self-Intersections: ", size(self_inter)

         do i = 1, size(inter)
            print *, inter(i)
         end do
         ! =====================================================================
         three_phase: block
            integer :: i, j
            ! Variables
            real(pr) :: alpha, beta
            real(pr), allocatable ::  lnKx(:), lnKy(:), X(:)
            
            real(pr) :: phase_x(nc), phase_y(nc), z(nc)
            type(injelope) :: px_bub_3, px_dew_3

            ! =================================================================
            !  Set variables based on intersections
            ! -----------------------------------------------------------------
            if (size(inter) == 0) then
            else
               i = inter(1)%i
               j = inter(1)%j

               alpha = inter(1)%x
               p = inter(1)%y

               lnKx = interpol( &
                      px_dew%alpha(i), px_dew%alpha(i + 1), &
                      px_dew%logk(i, :), px_dew%logk(i + 1, :), &
                      alpha &
                      )

               lnKy = interpol( &
                      px_bub%alpha(j), px_bub%alpha(j + 1), &
                      px_bub%logk(j, :), px_bub%logk(j + 1, :), &
                      alpha &
                      )
            end if
            ! =================================================================

            z = alpha*z_injection + (1 - alpha)*z_0

            ! Bubble line composition
            phase_y = exp(lnKy)*z
            ! Dew line composition
            phase_x = exp(lnKx)*z

            del_S0 = -0.1_pr
            beta = 1
            ns = 2*nc + 3

            ! ==================================================================
            !  Line with incipient phase gas
            ! ------------------------------------------------------------------
            lnKx = log(phase_x/phase_y)
            lnKy = log(z/phase_y)
            X = [lnKx, lnKy, log(p), alpha, beta]
            call injection_envelope_three_phase(X, ns, del_S0, px_bub_3)
            ! ==================================================================

            ! ==================================================================
            !  Line with incipient phase liquid
            ! ------------------------------------------------------------------
            lnKx = log(phase_y/phase_x)
            lnKy = log(z/phase_x)
            X = [lnKx, lnKy, log(p), alpha, beta]
            call injection_envelope_three_phase(X, ns, del_S0, px_dew_3)
            ! ==================================================================
         end block three_phase
      end block check_crossings
   end subroutine
   function px_two_phase(t_inj, pt_env_2, t_tol)
      !! Calculate two phase Px envelopes at a given injection temperature.
      !!
      !! Given an injection temperature `t_inj` and a base PT envelope 
      !! `pt_env_2`, finds all the points on the PT envelope near `t_inj`, based
      !! on an absolute tolerance `t_tol`. These points are used as 
      !! initialization for calculation of Px envelopes.
      
      use linalg, only: interpol
      use inj_envelopes, only: injection_envelope
      
      real(pr), intent(in) :: t_inj !! Injection temperature [K]
      type(envelope), intent(in) :: pt_env_2 !! Base PT envelope
      real(pr), intent(in) :: t_tol !! Absolute temperature tolerance
      type(injelope) :: px_two_phase !! Output Px envelope
      
      real(pr), allocatable :: ts_envel(:) !! Temperatures under tolerance 
      real(pr), allocatable :: k(:) !! K values
      real(pr), allocatable :: X(:) !! Vector of variables
      real(pr) :: alpha !! Amount of injection
      real(pr) :: p !! Pressure of ocurrence
      real(pr) :: pold !! Old pressure, used to assure no repeats

      integer :: i, idx, ns
      real(pr) :: del_S0

      del_S0 = 0.01_pr
      pold = 0

      ts_envel = pack(pt_env_2%t, mask=abs(pt_env_2%t - t_inj) < t_tol)
      do i = 1, size(ts_envel)
         idx = findloc(pt_env_2%t, value=ts_envel(i), dim=1)
         p = interpol( &
               pt_env_2%t(idx), pt_env_2%t(idx + 1), &
               pt_env_2%p(idx), pt_env_2%p(idx + 1), &
               t_inj)

         if (abs(p - pold) < 5) cycle
         pold = p

         k = exp(interpol( &
                  pt_env_2%t(idx), pt_env_2%t(idx + 1), &
                  pt_env_2%logk(idx, :), pt_env_2%logk(idx + 1, :), &
                  t_inj))

         X = [log(K), log(P), alpha]
         ns = size(X)

         call injection_envelope(X, ns, del_S0, px_two_phase)
      end do
   end function
end program main