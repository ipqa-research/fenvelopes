program main
   use envelopes, only: PTEnvel3
   use dtypes, only: envelope
   use inj_envelopes, only: injelope
   use constants, only: pr, ouput_path
   use legacy_ar_models, only: nc
   use flap, only: command_line_interface
   use stdlib_ansi, only: blue => fg_color_blue, red => fg_color_red, &
                          operator(//), operator(+), &
                          style_reset, style_blink_fast, style_bold, style_underline

   implicit none
   real(pr) :: et, st

   type(command_line_interface) :: cli
   integer :: cli_error
   character(len=99) :: cli_string

   type(envelope) :: pt_bub, pt_dew, pt_hpl !! Shared 2ph-PT envelopes
   type(PTEnvel3), allocatable :: pt_bub_3(:), pt_dew_3(:) !! Shared 3ph-PT envelopes
   type(injelope) :: px_bub, px_dew, px_hpl !! Shared 2ph-Px envelopes

   ! Setup everything
   call setup

   ! PT Envelopes
   call cpu_time(st)
   call pt_envelopes
   call cpu_time(et)
   print *, "PT: ", (et - st)*1000, "cpu ms"

   ! PX Envelopes
   call cpu_time(st)
   call px_envelopes
   call cpu_time(et)
   print *, "PX: ", (et - st)*1000, "cpu ms"
contains

   subroutine setup_cli
      !! Setup CLI subroutine
      !!
      !! Setup the Command-Line-Interface processor
      call cli%init(progname="envelopes", description="Phase Envelopes")
      call cli%add( &
         switch="--infile", &
         switch_ab="-i", &
         help="Input file", &
         error=cli_error, &
         required=.true.)
      call cli%parse(error=cli_error)

      if (cli_error /= 0) stop
   end subroutine

   subroutine setup
      !! Setup system
      !!
      !! Make output folder (if necessary) and/or clean everyhing in an
      !! existing one. Then read input files to setup needed parameters.
      use io_nml, only: read_system, write_system
      use inj_envelopes, only: setup_inj => from_nml
      integer :: funit_system
      character(len=500) :: infile

      call system("mkdir -p "//trim(ouput_path))
      call system("rm "//trim(ouput_path)//"*")

      call setup_cli
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
                           max_points, p_wilson, k_wilson, find_hpl, get_case
      use linalg, only: point
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

      type(point), allocatable :: intersections(:), self_intersections(:)
      character(len=:), allocatable :: pt_case

      integer :: n

      allocate (tv(max_points), pv(max_points), dv(max_points))
      allocate (k(size(z)))

      print *, style_underline // "PT Regions" // style_reset

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
      !  HPLL Envelope
      ! ------------------------------------------------------------------------
      t = 700.0_pr
      t = pt_bub%t(maxloc(pt_bub%p, dim=1))
      p = maxval(pt_bub%p)*2.5_pr

      call find_hpl(t, p, k)
      k = 1/k
      call envelope2( &
         3, nc, z, T, P, k, &
         n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, &
         pt_hpl &
         )
      ! ========================================================================

      ! ========================================================================
      !  Look for crossings
      ! ------------------------------------------------------------------------
      call get_case(&
         pt_dew, pt_bub, pt_hpl, &
         intersections, self_intersections, pt_case &
      )
      ! print *, style_bold // pt_case // style_reset
      ! ========================================================================

      three_phase: block
         use envelopes, only: pt_three_phase_from_intersection

         allocate(pt_bub_3(size(intersections)), pt_dew_3(size(intersections)))
         select case(pt_case)
         case("2_DEW_BUB")
            call pt_three_phase_from_intersection(&
                  pt_dew, pt_bub, intersections, &
                  pt_bub_3, pt_dew_3 &
            )
         case("1_HPL_DEW")
            call pt_three_phase_from_intersection(&
                  pt_hpl, pt_dew, intersections, &
                  pt_bub_3, pt_dew_3 &
            )
         end select
      end block three_phase
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
      
      print *, style_underline // "----------" // style_reset
      print *, style_underline // "Px Regions" // style_reset
      print *, style_underline // "----------" // style_reset

      ! ========================================================================
      !  Two phase envelopes
      ! ------------------------------------------------------------------------
      print *, red // "Running Bubble" // style_reset
      px_bub = px_two_phase(t_inj, pt_bub, t_tol)
      
      print *, blue // "Running Dew" // style_reset
      px_dew = px_two_phase(t_inj, pt_dew, t_tol)

      ! ========================================================================

      ! ========================================================================
      !  Look for crossings
      ! ------------------------------------------------------------------------
      inter = check_intersections(px_dew, px_bub)
      self_inter = check_self_intersections(px_dew)
      print *, style_bold // "Px Intersections:      " // style_reset, size(inter)
      print *, style_bold // "Px Self-Intersections: " // style_reset, size(self_inter)
      ! ========================================================================
      
      ! ========================================================================
      !  Three phase regions
      ! ------------------------------------------------------------------------
      three_phase: block
         integer :: i
         type(injelope) :: px_bub_3, px_dew_3, px_branch_3(2)

         if (size(inter) == 0) then
            px_bub_3 = px_three_phase(t_inj, pt_bub_3, t_tol)
            px_dew_3 = px_three_phase(t_inj, pt_dew_3, t_tol)
         else
            do i=1,size(inter)
               print *, "Intersection: ", inter(i)
               px_branch_3 = px_three_phase_from_inter(inter(i), px_dew, px_bub)
            end do
         end if

         if (size(self_inter) > 0) then
            do i=1,size(self_inter)
               !TODO: Add a check if one of the previous lines already
               !      in this DSP
               px_branch_3 = px_three_phase_from_inter(&
                  self_inter(i), px_dew, px_dew &
               )
            end do
         end if
      end block three_phase
      ! ========================================================================
   end subroutine

   function px_hpl_line(alpha_0, p)
      ! Find a HPLL PX line at a given pressure, starting from a given alpha
      use legacy_ar_models, only: termo
      use inj_envelopes, only: t
      real(pr), intent(in) :: alpha_0 !! Staring \(\alpha\) to search
      real(pr), intent(in) :: p !! Pressure of HPLL
      type(injelope) :: px_hpl_line !! Resulting HPLL line
      
      real(pr) :: diff
      real(pr) :: lnfug_z(nc), lnfug_y(nc)
      real(pr) :: z(nc), dzda(nc), y(nc), v, k(nc)
      
      real(pr), allocatable :: x(:)
      real(pr) :: del_S0
      integer :: ns, ncomp

      real(pr) :: alpha_in
      
      print *, red // "Running HPLL" // style_reset

      find_hpl: do ncomp = nc, 1, -1
         alpha_in = alpha_0
         y = 0
         y(ncomp) = 1
         diff = -1
         
         do while (diff < 0 .and. alpha_in < 0.9)
            alpha_in = alpha_in + 0.01_pr
            call get_z(alpha_in, z, dzda)
            call termo(nc, 4, 1, t, p, z, v, philog=lnfug_z)
            call termo(nc, 4, 1, t, p, y, v, philog=lnfug_y)
            diff = (log(z(ncomp)) + lnfug_z(ncomp)) - (log(y(ncomp)) + lnfug_y(ncomp))
         end do

         if (alpha_in < 1) exit find_hpl
      end do find_hpl

      k = 1/exp(lnfug_y - lnfug_z)
      X = [log(K), log(P), alpha_in]
      del_S0 = -0.005
      ns = size(X) - 1

      call injection_envelope(X, ns, del_S0, px_hpl_line)
   end function

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

         if (abs(p - pold) < 30) cycle
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
      !! `pt_env_3`, finds all the points on the PT envelope near `t_inj`, based
      !! on an absolute tolerance `t_tol`. These points are used as 
      !! initialization for calculation of Px envelopes.

      use linalg, only: interpol
      use inj_envelopes, only: injection_envelope_three_phase
      use stdlib_optval, only: optval
      use envelopes, only: AbsEnvel

      real(pr), intent(in) :: t_inj !! Injection temperature [K]
      type(PTEnvel3), intent(in) :: pt_env_3(:) !! Base PT envelopes
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

   function px_three_phase_from_inter(&
         inter, px_1, px_2, del_S0, beta0 &
         ) result(px_3)
      use legacy_ar_models, only: nc
      use stdlib_optval, only: optval
      use linalg, only: point, interpol
      use inj_envelopes, only: injelope, injection_envelope_three_phase, get_z
      type(point), intent(in) :: inter
      type(injelope), intent(in) :: px_1, px_2
      type(injelope) :: px_3(2)
      real(pr), optional :: del_S0
      real(pr), optional :: beta0

      integer :: i, j

      real(pr) :: lnKx(nc), lnKy(nc), alpha, p, beta, X(2*nc+3)
      real(pr) :: phase_y(nc), phase_x(nc)
      real(pr) :: del_S
      real(pr) :: z(nc), dzda(nc)
      integer :: ns

      del_S = optval(del_S0, -0.05_pr)
      beta = optval(beta0, 1.0_pr)

      i = inter%i
      j = inter%j

      alpha = inter%x
      p = inter%y

      lnKx = interpol( &
               px_1%alpha(i), px_1%alpha(i + 1), &
               px_1%logk(i, :), px_1%logk(i + 1, :), &
               alpha &
               )

      lnKy = interpol( &
               px_2%alpha(j), px_2%alpha(j + 1), &
               px_2%logk(j, :), px_2%logk(j + 1, :), &
               alpha &
               )

      call get_z(alpha, z, dzda)

      ! Bubble line composition
      phase_y = exp(lnKy)*z
      ! Dew line composition
      phase_x = exp(lnKx)*z

      ns = 2*nc + 3

      ! ==================================================================
      !  Line with incipient phase gas
      ! ------------------------------------------------------------------
      print *, "Three Phase: Gas"
      lnKx = log(phase_x/phase_y)
      lnKy = log(z/phase_y)
      X = [lnKx, lnKy, log(p), alpha, beta]
      call injection_envelope_three_phase(X, ns, del_S, px_3(1))
      ! ==================================================================

      ! ==================================================================
      !  Line with incipient phase liquid
      ! ------------------------------------------------------------------
      print *, "Three Phase: Liquid"
      lnKx = log(phase_y/phase_x)
      lnKy = log(z/phase_x)
      X = [lnKx, lnKy, log(p), alpha, beta]
      call injection_envelope_three_phase(X, ns, del_S, px_3(2))
   end function
end program main