program main
   use dtypes, only: envelope
   use envelopes, only: PTEnvel3
   use inj_envelopes, only: injelope, get_z
   use constants, only: pr, ouput_path
   use legacy_ar_models, only: nc, z
   use stdlib_ansi, only: blue => fg_color_blue, red => fg_color_red, &
      operator(//), operator(+), &
      style_reset, style_blink_fast, style_bold, style_underline
   use flap, only: command_line_interface
   use fenvelopes_cli, only: setup_cli => setup_cli_envelopes
   
   implicit none
   real(pr) :: et, st

   type(command_line_interface) :: cli

   type(envelope) :: pt_bub, pt_dew, pt_hpl !! Shared 2ph-PT envelopes
   type(PTEnvel3), allocatable :: pt_bub_3(:), pt_dew_3(:) !! Shared 3ph-PT envelopes
   type(injelope), allocatable :: px_bub(:), px_dew(:) !! Shared 2ph-Px envelopes
   type(injelope) :: px_hpl(2)

   integer :: cli_error
   logical :: run_px, run_3ph

   real(pr) :: alpha=0.0
   real(pr) :: pt_bub_t0 = 180
   real(pr) :: pt_dew_t0 = 180

   ! Setup everything
   call setup

   call cli%get(alpha, switch="--alpha0", error=cli_error)
   call get_z(alpha, z)
   print "(A)", "=================================================================="
   print "(A)", "Settings"
   print "(A)", "------------------------------------------------------------------"
   print "('𝛼: '(F5.3))", alpha
   if (size(z) > 14) then
      print "(A)", "Z:"
      print "(*(F5.3, 1x)))", z(:14)
      print "(*(F5.3, 1x)))", z(15:)
   else
   print "('Z: '(*(F5.3, 1x)))", z
   end if

   call cli%get(run_3ph, switch="--three_phases", error=cli_error)

   ! PT Envelopes
   call cpu_time(st)
   call pt_envelopes
   call cpu_time(et)
   print *, "PT: ", (et - st)*1000, "cpu ms"

   call cli%get(run_px, switch="--injection", error=cli_error)

   if (run_px) then
      ! PX Envelopes
      call cpu_time(st)
      call px_envelopes
      call cpu_time(et)
      print *, "PX: ", (et - st)*1000, "cpu ms"
   end if
contains
   subroutine setup
      !! Setup system
      !!
      !! Make output folder (if necessary) and/or clean everyhing in an
      !! existing one. Then read input files to setup needed parameters.
      use io_nml, only: read_system, write_system
      use inj_envelopes, only: setup_inj => from_nml
      integer :: funit_system, cli_error
      character(len=500) :: infile

      call system("mkdir -p "//trim(ouput_path))
      call system("rm "//trim(ouput_path)//"*")

      call setup_cli(cli)
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
      use saturation_points, only: bubble_temperature, dew_temperature, hpl_temperature, EquilibriaState
      !! Calculation of PT envelopes of the main system.
      real(pr), allocatable :: tv(:) ! Temperatures [K]
      real(pr), allocatable :: pv(:) ! Pressures [bar]
      real(pr), allocatable :: dv(:) ! Densities [mol/L]

      real(pr) :: tcri(4)            ! Critical points temperatures
      real(pr) :: pcri(4)            ! Critical points pressures
      real(pr) :: dcri(4)            ! Critical points densities

      real(pr) :: t, p               ! Temperature and pressure
      real(pr), allocatable :: k(:)  ! K factors
      integer :: n_points, icri(4), ncri, i, i_max

      type(point), allocatable :: intersections(:), self_intersections(:)
      character(len=:), allocatable :: pt_case

      integer :: n
      type(EquilibriaState) :: bubble, dew

      allocate (tv(max_points), pv(max_points), dv(max_points))
      allocate (k(size(z)))

      print "(A)", "=================================================================="
      print *, style_underline // "PT Regions" // style_reset
      print "(A)", "------------------------------------------------------------------"

      ! ========================================================================
      !  Bubble envel
      ! ------------------------------------------------------------------------
      print *, "Bubble PT"
      bubble = bubble_temperature(z, p=1.0_pr, t0=pt_bub_t0)
      k = bubble%y/z
      t = bubble%t
      p = bubble%p

      call envelope2( &
         1, nc, z, T, P, k, &
         n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, &
         pt_bub &
         )
      print *, "Done"
      ! ========================================================================

      ! ========================================================================
      !  Dew/AOP envelopes
      ! ------------------------------------------------------------------------
      print *, "Dew PT"
      t = pt_dew_t0
      p = p_wilson(z, t)
      do while (p > 0.1)
         t = t - 5
         p = p_wilson(z, t)
      end do

      k = 1/k_wilson(t, p)

      dew = dew_temperature(z, p=0.00001_pr, t0=pt_dew_t0)
      k = dew%y/z
      t = dew%t
      p = dew%p

      call envelope2( &
         2, nc, z, T, P, k, &
         n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, &
         pt_dew &
         )
      print *, "Done"
      ! ========================================================================

      ! ========================================================================
      !  HPLL Envelope
      ! ------------------------------------------------------------------------
      print *, "HPLL PT"
      i_max = maxloc([pt_dew%p], dim=1)

      p = pt_dew%p(i_max) + 50
      t = pt_dew%t(i_max) + 50
      
      dew = hpl_temperature(z, p=p, t0=t, y0=exp(pt_dew%logk(i_max, :)) * z)
      k = dew%y/z
      t = dew%t
      p = dew%p

      call find_hpl(t, p, k)
      call envelope2(&
         3, nc, z, T, P, k, &
         n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, &
         pt_hpl &
         )
      print *, "Done"
      ! ========================================================================

      ! ========================================================================
      !  Look for crossings
      ! ------------------------------------------------------------------------
      call get_case(&
         pt_dew, pt_bub, pt_hpl, &
         intersections, self_intersections, pt_case &
         )
      print *, style_bold // pt_case // style_reset
      ! ========================================================================

      if (run_3ph) then
         three_phase: block
            use linalg, only: allclose
            use envelopes, only: pt_three_phase_from_intersection

            if (size(intersections) >= 1) then
               allocate(pt_bub_3(size(intersections)), pt_dew_3(size(intersections)))
            else
               allocate(pt_bub_3(1), pt_dew_3(1))
            end if

            select case(pt_case)
             case("0")
               isolated: block
                  use legacy_ar_models, only: nc, z
                  use legacy_thermo_properties, only: termo
                  use envelopes, only: pt_envelope_three_phase
                  integer :: ncomp=7
                  real(pr) :: p, v, t
                  real(pr) :: lnphix(nc), lnphiy(nc), lnphiz(nc), &
                     phase_x(nc), phase_y(nc), beta, x(2*nc+3), lnKx(nc), lnKy(nc)

                  type(envelope) :: env2ph

                  !!! OJO: Esta variable se usa como selectora de que linea
                  !!!      utilizar al momento de calcular la linea trifasica,
                  !!!      usualmente usamos pt_bub
                  env2ph = pt_bub

                  beta = z(ncomp)
                  phase_y = 0
                  phase_y(ncomp) = 1
                  p = env2ph%p(1)
                  t = env2ph%t(1)

                  call termo(nc, 1, 4, t, p, phase_y, v, philog=lnphiy)
                  call termo(nc, 1, 4, t, p, z, v, philog=lnphiz)

                  ! Y: Asphaltenes
                  ! X: Vapor
                  ! Z: Main fluid
                  phase_x = exp(env2ph%logk(1, :)) * z
                  call termo(nc, -1, 4, t, p, phase_x, v, philog=lnphix)

                  ! main/vapor
                  lnKx = log(z/phase_x)
                  ! asph/vapor
                  lnKy = lnphix - lnphiy

                  X = [lnKx, lnKy, log(p), log(t), beta]
                  call pt_envelope_three_phase(X, 2*nc+2, 0.1_pr, pt_bub_3(1))
               end block isolated
             case("2_DEW_BUB")
               call pt_three_phase_from_intersection(&
                  pt_dew, pt_bub, intersections, &
                  pt_bub_3, pt_dew_3 &
                  )
             case("2_HPL_BUB_DEW_BUB")
               call pt_three_phase_from_intersection(&
                  pt_hpl, pt_bub, [intersections(1)], &
                  pt_bub_3, pt_dew_3 &
                  )
               call pt_three_phase_from_intersection(&
                  pt_dew, pt_bub, [intersections(2)], &
                  pt_bub_3, pt_dew_3 &
                  )
             case("2_HPL_BUB_HPL_DEW")
               call pt_three_phase_from_intersection(&
                  pt_hpl, pt_dew, [intersections(2)], &
                  pt_bub_3, pt_dew_3 &
                  )
               if (.not. allclose(&
                  [intersections(1)%x, intersections(1)%y], &
                  [intersections(2)%x, intersections(2)%y], &
                  atol=1.0_pr &
               )) call pt_three_phase_from_intersection(&
                  pt_hpl, pt_bub, [intersections(1)], &
                  pt_bub_3, pt_dew_3 &
                  )

             case("2_HPL_BUB")
               call pt_three_phase_from_intersection(&
                  pt_hpl, pt_bub, [intersections(1)], &
                  pt_bub_3, pt_dew_3 &
                  )
               call pt_three_phase_from_intersection(&
                  pt_hpl, pt_bub, [intersections(2)], &
                  pt_bub_3, pt_dew_3 &
                  )
             case("1_HPL_BUB")
               call pt_three_phase_from_intersection(&
                  pt_hpl, pt_bub, intersections, &
                  pt_dew_3, pt_bub_3 &
                  )
             case("1_HPL_DEW")
               call pt_three_phase_from_intersection(&
                  pt_hpl, pt_dew, intersections, &
                  pt_bub_3, pt_dew_3 &
                  )
             case("1_DEW")
               call pt_three_phase_from_intersection(&
                  pt_dew, pt_dew, self_intersections, &
                  pt_dew_3, pt_bub_3 &
                  )
            end select
         end block three_phase
      end if
   end subroutine

   subroutine px_envelopes
      !! Calculation of P\(\alpha\) envelopes at selected temperature.
      use inj_envelopes, only: full_newton, z_injection, &
         T_inj => T, injection_envelope, z_0, &
         injelope, injection_envelope_three_phase, get_z, &
         px_two_phase_from_pt, &
         px_three_phase_from_pt, &
         px_three_phase_from_inter, &
         px_hpl_line
      use envelopes, only: envelope, k_wilson, p_wilson
      use linalg, only: interpol, point, intersection

      real(pr) :: t_tol = 2, p_hpl, a_hpl, z_hpl(nc), y_hpl(nc)
      integer :: max_bub_p

      print *, style_underline // "----------" // style_reset
      print *, style_underline // "Px Regions" // style_reset
      print *, style_underline // "----------" // style_reset

      ! ========================================================================
      !  Two phase envelopes
      ! ------------------------------------------------------------------------
      print *, red // "Running Bubble" // style_reset
      px_bub = px_two_phase_from_pt(t_inj, pt_bub, alpha0=alpha, t_tol=5.0_pr)

      print *, blue // "Running Dew" // style_reset
      px_dew = px_two_phase_from_pt(t_inj, pt_dew, alpha0=alpha, t_tol=5.0_pr)

      print *, blue // "Running HPLL" // style_reset
      manual: block
         real(pr), allocatable :: X(:)
         real(pr) :: k(nc), del_S0, p_0
         integer :: ns
         K = 0.01
         K(nc) = 10
         p_0 = 1100
         X = [log(k), log(p_0), 0.5_pr]
         del_S0 = 0.5_pr
         ns = size(X) - 1
         call injection_envelope(X, ns, del_S0, px_hpl(1))
         
         X(1:nc) = px_hpl(1)%logk(2, :)
         X(nc + 2) = px_hpl(1)%alpha(2)

         call injection_envelope(X, ns, -del_S0, px_hpl(2))
      end block manual
      ! if (allocated(pt_hpl%t)) then
      !     px_hpl = px_two_phase_from_pt(t_inj, pt_hpl, alpha0=alpha, t_tol=5.0_pr)
      ! end if

      ! ========================================================================
      !  Three phase regions
      ! ------------------------------------------------------------------------
      if (run_3ph) then
         three_phase: block
            use dsp_lines, only: injelope, dsp_line_from_dsp_px
            use inj_envelopes, only: PXEnvel3
            type(PXEnvel3) :: px_bub_3, px_dew_3, px_branch_3(2)

            ! =====================================================================
            ! Intersections between lines
            ! ---------------------------------------------------------------------
            print *, "========================================================="
            print *, "DSPs"
            print *, "DEW_BUB"
            call calc_all_dsps(px_dew, px_bub, px_branch_3)
            print *, "BUB_HPL"
            call calc_all_dsps(px_bub, px_hpl, px_branch_3)
            print *, "DEW_HPL"
            call calc_all_dsps(px_dew, px_hpl, px_branch_3)

            print *, "SELF DEW"
            call calc_all_self_dsp(px_dew, px_branch_3)
            print *, "SELF BUB"
            call calc_all_self_dsp(px_bub, px_branch_3)
            print *, "========================================================="
            ! =====================================================================

            ! Isolated lines coming from PT lines
            print *, "Isolated Bub"
            if (allocated(pt_bub_3)) then
               px_bub_3 = px_three_phase_from_pt(t_inj, pt_bub_3, 2*t_tol, alpha0=alpha)
            end if

            if (px_bub_3%beta(size(px_bub_3%beta)) > 1) then
               !! Get back if the 3ph isolated line ended in a DSP
               comeback: block
                  use inj_envelopes, only: get_z
                  integer :: npoints
                  real(pr), allocatable :: xx(:), y(:), w(:)
                  real(pr), allocatable :: X(:), z(:), z_inj(:), lnKx(:), lnKy(:)
                  real(pr) :: p, alpha, beta

                  exit comeback

                  npoints = size(px_bub_3%beta) - 1
                  alpha = px_bub_3%alpha(npoints)
                  beta = 1!px_bub_3%beta(npoints)
                  p = px_bub_3%p(npoints)

                  xx = px_bub_3%x(npoints, :)
                  y = px_bub_3%y(npoints, :)
                  w = px_bub_3%w(npoints, :)

                  z = px_bub_3%z_0
                  z_inj = px_bub_3%z_inj
                  call get_z(alpha, z)

                  lnKx = log(xx/w)
                  lnKy = log(y/w)

                  X = [lnKx, lnKy, log(p), alpha, beta]
                  call injection_envelope_three_phase(X, 2*nc+3, -0.01_pr, px_bub_3)

                  X = [log(xx/z), log(p), alpha]
                  call injection_envelope(X, nc+2, -0.01_pr, px_bub(1))
               end block comeback
            end if

            print *, "Isolated Dew"
            if (allocated(pt_dew_3)) then
               px_dew_3 = px_three_phase_from_pt(t_inj, pt_dew_3, 2*t_tol, alpha0=alpha)
            end if
         end block three_phase
      end if
      ! ========================================================================
   end subroutine

   subroutine calc_all_dsps(px_1, px_2, px_out)
      !! Calculate three-phase lines from DSPs between two lines.
      !!
      !! Find all the intersection points between two phase envelopes and
      !! calculate all the pairs of three phase lines (\(P \alpha\) and DSP
      !! lines)
      use inj_envelopes, only: injelope, px_three_phase_from_inter, PXEnvel3
      use dsp_lines, only: dsp_line_from_dsp_px
      use linalg, only: point, intersection
      type(injelope) :: px_1(:), px_2(:)
      type(PXEnvel3) :: px_out(2)
      type(injelope) :: dsps(2)
      type(point), allocatable :: inter_1_2(:)
      integer :: i, j, k

      do i=1,size(px_1)
         if (.not. allocated(px_1(i)%alpha)) cycle
         do j=1,size(px_2)
            if (.not. allocated(px_2(j)%alpha)) cycle
            print *, i, j
            ! Go through each possible pair of envelopes to find DSPs
            inter_1_2  = intersection(&
               px_1(i)%alpha, px_1(i)%p, &
               px_2(j)%alpha, px_2(j)%p  &
               )
            do k=1,size(inter_1_2)
               ! For each DSP found, calculate the two possible DSP lines
               ! and the two three-phase branchs.
               print *, "Intersection: ", inter_1_2(k)
               ! dsps = dsp_line_from_dsp_px(&
               !    inter_1_2(k), px_1(i), px_2(j) &
               !    )
               px_out = px_three_phase_from_inter(&
                  inter_1_2(k), px_1(i), px_2(j) &
                  )
            end do
         end do
      end do
   end subroutine

   subroutine calc_all_self_dsp(px, px_out)
      !! Calculate three-phase lines from DSPs on a single line.
      !!
      !! From a single envelope find all self-crossing points an calculate
      !! the corresponding three-phase lines for each one of them. Including
      !! both Px as well as DSP lines.
      use inj_envelopes, only: injelope, px_three_phase_from_inter, PXEnvel3
      use dsp_lines, only: dsp_line_from_dsp_px
      use linalg, only: point, intersection
      type(injelope) :: px(:)

      type(PXEnvel3) :: px_out(2)
      type(injelope) :: dsps(2)
      type(point), allocatable :: inter(:)
      integer :: i, j

      do i=1,size(px)
         print *, i
         inter = intersection(px(i)%alpha, px(i)%p)
         if (size(inter) > 0) then
            do j=1,size(inter)
               print *, "Self-intersection: ", inter(j)
               ! dsps = dsp_line_from_dsp_px(inter(j), px(i), px(i))
               px_out = px_three_phase_from_inter(&
                  inter(j), px(i), px(i) &
                  )
            end do
         end if
      end do
   end subroutine
end program main
