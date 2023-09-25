program main
   use envelopes, only: PTEnvel3
   use dtypes, only: envelope
   use inj_envelopes, only: injelope
   use constants, only: pr, ouput_path
   use legacy_ar_models, only: nc
   use flap, only: command_line_interface
   use stdlib_ansi, only: blue => fg_color_blue, red => fg_color_red, &
                          operator(//), operator(+), &
                          style_reset, style_blink_fast, style_bold

   implicit none
   real(pr) :: et, st

   type(command_line_interface) :: cli
   integer :: cli_error
   character(len=99) :: cli_string

   type(envelope) :: pt_bub, pt_dew
   type(injelope) :: px_bub, px_dew
   type(PTEnvel3), allocatable :: pt_bub_3(:), pt_dew_3(:)

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
   print *, "PT: ", (et - st)*1000, "cpu ms"

   call cpu_time(st)
   call px_envelopes ! Calculate Px envelopes
   call cpu_time(et)
   print *, "PX: ", (et - st)*1000, "cpu ms"
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
         call exit

         three_phase: block
         use linalg, only: interpol
         use envelopes, only: pt_envelope_three_phase, PTEnvel3
         real(pr), allocatable :: lnKx(:), lnKy(:)
         real(pr), allocatable :: X(:)
         real(pr) :: t, p, beta, del_S0
         real(pr), allocatable :: phase_y(:), phase_x(:)
         integer :: i, j, i_inter=1
         integer :: ns

         allocate(pt_bub_3(size(inter)), pt_dew_3(size(inter)))
         do i_inter=1,size(inter)
            i = inter(i_inter)%i
            j = inter(i_inter)%j

            t = inter(i_inter)%x
            p = inter(i_inter)%y

            lnKx = interpol( pt_dew%t(i), pt_dew%t(i + 1), &
                           pt_dew%logk(i, :), pt_dew%logk(i + 1, :), &
                           t &
                           )

            lnKy = interpol( &
                     pt_bub%t(j), pt_bub%t(j + 1), &
                     pt_bub%logk(j, :), pt_bub%logk(j + 1, :), &
                     t &
                     )

            ! Bubble line composition
            phase_y = exp(lnKy)*z
            ! Dew line composition
            phase_x = exp(lnKx)*z

            del_S0 = -0.01_pr
            beta = 1

            ns = 2*nc + 3

            ! ==================================================================
            !  Line with incipient phase gas
            ! ------------------------------------------------------------------
            print *, "Three Phase: Gas"
            lnKx = log(phase_x/phase_y)
            lnKy = log(z/phase_y)
            X = [lnKx, lnKy, log(p), log(t), beta]
            call pt_envelope_three_phase(X, ns, del_S0, pt_bub_3(i_inter))
            ! ==================================================================
            ! ==================================================================
            !  Line with incipient phase liquid
            ! ------------------------------------------------------------------
            print *, "Three Phase: Liquid"
            lnKx = log(phase_y/phase_x)
            lnKy = log(z/phase_x)
            X = [lnKx, lnKy, log(p), log(t), beta]
            call pt_envelope_three_phase(X, ns, del_S0, pt_dew_3(i_inter))
            ! ==================================================================
         end do
         end block three_phase
      end block check_crossings
      ! ========================================================================
   end subroutine

   subroutine px_envelopes
      !! Calculation of Px envelopes at selected temperature.
      use inj_envelopes, only: full_newton, z_injection, &
                               T_inj => T, injection_envelope, z_0, &
                               injelope, injection_envelope_three_phase, get_z
      use envelopes, only: envelope, k_wilson, p_wilson
      use linalg, only: interpol, point
      
      type(point), allocatable :: inter(:), self_inter(:)

      real(pr) :: del_S0
      integer :: ns, i
      real(pr) :: p
      real(pr) :: t_tol = 2
      real(pr) :: dzda(nc)

      ! ========================================================================
      !  Two phase envelopes
      ! ------------------------------------------------------------------------
      ! print *, red // "Running Bubble" // style_reset
      ! px_bub = px_two_phase(t_inj, pt_bub, t_tol, del_S0=-0.1_pr)
      
      ! print *, blue // "Running Dew" // style_reset
      ! px_dew = px_two_phase(t_inj, pt_dew, t_tol, del_S0=-0.1_pr)

      print *, red // "Running Bubble" // style_reset
      px_bub = px_two_phase(t_inj, pt_bub, t_tol)
      
      print *, blue // "Running Dew" // style_reset
      px_dew = px_two_phase(t_inj, pt_dew, t_tol)
      ! ========================================================================

      ! ========================================================================
      !  Look for crossings
      ! ------------------------------------------------------------------------
      check_crossings: block
         inter = check_intersections(px_dew, px_bub)
         self_inter = check_self_intersections(px_dew)
         print *, style_bold // "Px Intersections:      " // style_reset, size(inter)
         print *, style_bold // "Px Self-Intersections: " // style_reset, size(self_inter)
      end block check_crossings
      ! ========================================================================
      
      ! ========================================================================
      !  Three phase regions
      ! ------------------------------------------------------------------------
      three_phase: block
         use legacy_ar_models, only: TERMO
         integer :: i, j, i_inter
         ! Variables
         real(pr) :: alpha, beta
         real(pr), allocatable ::  lnKx(:), lnKy(:), X(:)
         
         real(pr) :: phase_x(nc), phase_y(nc), z(nc)
         real(pr) :: lnfug_x(nc), dlnphi_dp(nc), dlnphi_dt(nc), dlnphi_dn(nc,nc), v
         real(pr) :: lnfug_y(nc)
         type(injelope) :: px_bub_3, px_dew_3

         ! =====================================================================
         !  Set variables based on intersections
         ! ---------------------------------------------------------------------
         if (size(inter) == 0) then
            px_bub_3 = px_three_phase(t_inj, pt_bub_3, t_tol)
            px_dew_3 = px_three_phase(t_inj, pt_dew_3, t_tol)
         else
            do i_inter=1,size(inter)
               print *, "Intersection: ", inter(i_inter)
            i = inter(i_inter)%i
            j = inter(i_inter)%j

            alpha = inter(i_inter)%x
            p = inter(i_inter)%y

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

            call get_z(alpha, z, dzda)

            ! Bubble line composition
            phase_y = exp(lnKy)*z
            ! Dew line composition
            phase_x = exp(lnKx)*z

            if (i_inter == 1) then
               del_S0 = -0.01_pr
               beta = 1
            else
               del_S0 = 0.01_pr
               beta = 0
            end if

            ns = 2*nc + 3

            ! ==================================================================
            !  Line with incipient phase gas
            ! ------------------------------------------------------------------
            print *, "Three Phase: Gas"
            lnKx = log(phase_x/phase_y)
            lnKy = log(z/phase_y)
            X = [lnKx, lnKy, log(p), alpha, beta]
            call injection_envelope_three_phase(X, ns, del_S0, px_bub_3)
            ! ==================================================================

            ! ==================================================================
            !  Line with incipient phase liquid
            ! ------------------------------------------------------------------
            print *, "Three Phase: Liquid"
            lnKx = log(phase_y/phase_x)
            lnKy = log(z/phase_x)
            X = [lnKx, lnKy, log(p), alpha, beta]
            call injection_envelope_three_phase(X, ns, del_S0, px_dew_3)
            end do
            ! ==================================================================
         end if
      end block three_phase
      ! ========================================================================
   end subroutine

   function check_intersections(this, other) result(intersections)
      !! Find intersections between two injelopes.
      use linalg, only: point, intersection
      type(injelope), intent(in) :: this !! One of the injelopes
      type(injelope), intent(in) :: other !! The other injelope
      type(point), allocatable :: intersections(:) !! Found intersections
      intersections = intersection( &
                 this%alpha, this%p, &
                 other%alpha, other%p  &
      )
   end function

   function check_self_intersections(self) result(self_intersections)
      !! Find self-intersections on a given injelope.
      use linalg, only: point, intersection
      type(injelope), intent(in) :: self !! Envelope to check self-intersections
      type(point), allocatable :: self_intersections(:) !! Found intersections
      self_intersections = intersection(self%alpha, self%p)
   end function

   function px_two_phase(t_inj, pt_env_2, t_tol, del_S0)
      !! Calculate two phase Px envelopes at a given injection temperature.
      !!
      !! Given an injection temperature `t_inj` and a base PT envelope 
      !! `pt_env_2`, finds all the points on the PT envelope near `t_inj`, based
      !! on an absolute tolerance `t_tol`. These points are used as 
      !! initialization for calculation of Px envelopes.
      
      use linalg, only: interpol
      use inj_envelopes, only: injection_envelope
      use stdlib_optval, only: optval
      use envelopes, only: AbsEnvel
      
      real(pr), intent(in) :: t_inj !! Injection temperature [K]
      type(envelope), intent(in) :: pt_env_2 !! Base PT envelope
      real(pr), intent(in) :: t_tol !! Absolute temperature tolerance
      real(pr), optional, intent(in) :: del_S0 !! First point \(\Delta S\)
      type(injelope) :: px_two_phase !! Output Px envelope
      
      real(pr), allocatable :: ts_envel(:) !! Temperatures under tolerance 
      real(pr), allocatable :: k(:) !! K values
      real(pr), allocatable :: X(:) !! Vector of variables
      real(pr) :: alpha !! Amount of injection
      real(pr) :: p !! Pressure of ocurrence
      real(pr) :: pold !! Old pressure, used to assure no repeats

      integer :: i, idx, ns
      real(pr) :: del_S

      del_S = optval(del_S0, 0.1_pr)
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

         alpha = 0
         X = [log(K), log(P), alpha]
         ns = size(X)

         call injection_envelope(X, ns, del_S, px_two_phase)
      end do
   end function

   function px_three_phase(t_inj, pt_env_3, t_tol, del_S0)
      !! Calculate three phase Px envelopes at a given injection temperature.
      !!
      !! Given an injection temperature `t_inj` and a base PT envelope 
      !! `pt_env_2`, finds all the points on the PT envelope near `t_inj`, based
      !! on an absolute tolerance `t_tol`. These points are used as 
      !! initialization for calculation of Px envelopes.
      
      use linalg, only: interpol
      use inj_envelopes, only: injection_envelope_three_phase
      use stdlib_optval, only: optval
      use envelopes, only: AbsEnvel
      
      real(pr), intent(in) :: t_inj !! Injection temperature [K]
      type(PTEnvel3), intent(in) :: pt_env_3(:) !! Base PT envelope
      real(pr), intent(in) :: t_tol !! Absolute temperature tolerance
      real(pr), optional, intent(in) :: del_S0 !! First point \(\Delta S\)
      type(injelope) :: px_three_phase !! Output Px envelope
      
      real(pr), allocatable :: ts_envel(:) !! Temperatures under tolerance 
      real(pr), allocatable :: kx(:), ky(:) !! K values
      real(pr), allocatable :: X(:) !! Vector of variables
      real(pr) :: alpha !! Amount of injection
      real(pr) :: beta  !! Main phases molar fraction
      real(pr) :: p !! Pressure of ocurrence
      real(pr) :: pold !! Old pressure, used to assure no repeats

      integer :: i, idx, ns, i_envel
      real(pr) :: del_S

      del_S = optval(del_S0, 0.05_pr)
      pold = 0

      do i_envel = 1,size(pt_env_3)
      associate(pt => pt_env_3(i_envel))
      ts_envel = pack(pt%t, mask=abs(pt%t - t_inj) < t_tol)
      do i = 1, size(ts_envel)
         idx = findloc(pt%t, value=ts_envel(i), dim=1)
         p = interpol( &
               pt%t(idx), pt%t(idx + 1), &
               pt%p(idx), pt%p(idx + 1), &
               t_inj)

         if (abs(p - pold) < 5) cycle
         pold = p

         kx = exp(interpol( &
                  pt%t(idx), pt%t(idx + 1), &
                  pt%lnkx(idx, :), pt%lnkx(idx + 1, :), &
                  t_inj))
         ky = exp(interpol( &
                  pt%t(idx), pt%t(idx + 1), &
                  pt%lnky(idx, :), pt%lnky(idx + 1, :), &
                  t_inj))
         beta = interpol( &
                  pt%t(idx), pt%t(idx + 1), &
                  pt%beta(idx), pt%beta(idx + 1), &
                  t_inj)
         
         alpha = 0
         X = [log(Kx), log(Ky), log(P), alpha, beta]
         ns = size(X) - 1

         call injection_envelope_three_phase(X, ns, del_S, px_three_phase)
      end do
      end associate
      end do
   end function
end program main