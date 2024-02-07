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
   type(injelope) :: px_hpl(1)

   integer :: cli_error
   logical :: run_px, run_3ph

   real(pr) :: alpha=0.0
   real(pr) :: pt_bub_t0 = 180
   real(pr) :: pt_dew_t0 = 180

   ! Setup everything
   call setup
   
   call cli%get(alpha, switch="--alpha0", error=cli_error)
   call get_z(alpha, z)
   print *, "ALPHA: ", alpha
   
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

      print *, style_underline // "PT Regions" // style_reset

      ! ========================================================================
      !  Bubble envel
      ! ------------------------------------------------------------------------
      ! call k_wilson_bubble(z, t_0=pt_bub_t0, p_end=0.5_pr, t=t, p=p, k=k)
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
      ! ========================================================================

      ! ========================================================================
      !  HPLL Envelope
      ! ------------------------------------------------------------------------
      i_max = maxloc([pt_dew%p], dim=1)
      
      p = pt_dew%p(i_max) + 50
      t = pt_dew%t(i_max) + 50

      dew = hpl_temperature(z, p=p, t0=t, y0=exp(pt_dew%logk(i_max, :)) * z)
      k = dew%y/z
      t = dew%t
      p = dew%p

      ! call find_hpl(t, p, k)
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
      print *, style_bold // pt_case // style_reset
      ! ========================================================================

      if (run_3ph) then
      three_phase: block
         use envelopes, only: pt_three_phase_from_intersection

         if (size(intersections) >= 1) then
            allocate(pt_bub_3(size(intersections)), pt_dew_3(size(intersections)))
         else
            allocate(pt_bub_3(1), pt_dew_3(1))
         end if
         
         select case(pt_case)
         case("0")
            isolated: block
               use legacy_ar_models, only: nc, termo, z
               use envelopes, only: pt_envelope_three_phase
               integer :: ncomp=1
               real(pr) :: p, v, t
               real(pr) :: lnphix(nc), lnphiy(nc), &
                           phase_x(nc), phase_y(nc), beta, x(2*nc+3), lnKx(nc), lnKy(nc)

               beta = z(ncomp)
               phase_y = 0
               phase_y(ncomp) = 1
               p = pt_bub%p(1)
               t = pt_bub%t(1)
               
               call termo(nc, 1, 4, t, p, phase_y, v, philog=lnphiy)

               ! Y: Asphaltenes
               ! X: Vapor
               ! Z: Main fluid
               phase_x = exp(pt_bub%logk(1, :)) * z
               call termo(nc, -1, 4, t, p, phase_x, v, philog=lnphix)
               
               ! main/vapor
               lnKx = log(z/phase_x)
               ! asph/vapor
               lnKy = lnphix - lnphiy

               X = [lnKx, lnKy, log(p), log(t), beta]
               call pt_envelope_three_phase(X, 2*nc+2, 0.01_pr, pt_bub_3(1))
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
            call pt_three_phase_from_intersection(&
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
      ! TODO: This is a dirty setup that barely works and should be fixed

      ! if (allocated(px_dew)) then
      ! if (size(px_dew(1)%alpha) > 5) then
      !    max_bub_p = maxloc(px_dew(1)%p, dim=1)
      !    a_hpl = px_dew(1)%alpha(max_bub_p)*0.9
      !    p_hpl = px_dew(1)%p(max_bub_p) * 1.1
      !    call get_z(a_hpl, z_hpl)
      !    y_hpl = exp(px_dew(1)%logk(max_bub_p+2, :)) * z
      !    print *, a_hpl, p_hpl
      !    px_hpl(1) = px_hpl_line(a_hpl, p_hpl, y0=y_hpl)
      ! endif
      ! if (size(px_bub(1)%alpha) > 5) then
      !    max_bub_p = maxloc(px_bub(1)%p, dim=1)
      !    a_hpl = px_bub(1)%alpha(max_bub_p)*1.1
      !    p_hpl = px_bub(1)%p(max_bub_p)
      !    call get_z(a_hpl, z_hpl)
      !    y_hpl = exp(px_bub(1)%logk(max_bub_p+2, :)) * z
      !    px_hpl(1) = px_hpl_line(a_hpl, p_hpl, y0=y_hpl)
      ! end if
      ! end if

      ! ========================================================================

      ! ========================================================================
      !  Three phase regions
      ! ------------------------------------------------------------------------
      if (run_3ph) then
      three_phase: block
         use dsp_lines, only: injelope, dsp_line_from_dsp_px
         type(injelope) :: px_bub_3, px_dew_3, px_branch_3(2)

         ! =====================================================================
         ! Intersections between lines
         ! ---------------------------------------------------------------------
         call calc_all_dsps(px_dew, px_bub, px_branch_3)
         call calc_all_dsps(px_bub, px_hpl, px_branch_3)
         call calc_all_dsps(px_dew, px_hpl, px_branch_3)

         call calc_all_self_dsp(px_dew, px_branch_3)
         call calc_all_self_dsp(px_bub, px_branch_3)
         ! =====================================================================

         ! Isolated lines coming from PT lines
         print *, "Isolated Bub"
         if (allocated(pt_bub_3)) then
            px_bub_3 = px_three_phase_from_pt(t_inj, pt_bub_3, 2*t_tol)
         end if
         
         ! print *, "Isolated Dew"
         if (allocated(pt_dew_3)) then
            px_dew_3 = px_three_phase_from_pt(t_inj, pt_dew_3, t_tol)
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
      use inj_envelopes, only: injelope, px_three_phase_from_inter
      use dsp_lines, only: dsp_line_from_dsp_px
      use linalg, only: point, intersection
      type(injelope) :: px_1(:), px_2(:)
      type(injelope) :: px_out(2)
      type(injelope) :: dsps(2)
      type(point), allocatable :: inter_1_2(:)
      integer :: i, j, k
      
      do i=1,size(px_1)
         do j=1,size(px_2)
            ! Go through each possible pair of envelopes to find DSPs
            inter_1_2  = intersection(&
               px_1(i)%alpha, px_1(i)%p, &
               px_2(j)%alpha, px_2(j)%p  &
            )
            do k=1,size(inter_1_2)
               ! For each DSP found, calculate the two possible DSP lines
               ! and the two three-phase branchs.
               print *, "Intersection: ", inter_1_2(k)
               dsps = dsp_line_from_dsp_px(&
                  inter_1_2(k), px_1(i), px_2(j) &
               )
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
      use inj_envelopes, only: injelope, px_three_phase_from_inter
      use dsp_lines, only: dsp_line_from_dsp_px
      use linalg, only: point, intersection
      type(injelope) :: px(:)
      
      type(injelope) :: px_out(2)
      type(injelope) :: dsps(2)
      type(point), allocatable :: inter(:)
      integer :: i, j

      do i=1,size(px)
         inter = intersection(px(i)%alpha, px(i)%p)
         if (size(inter) > 0) then
            do j=1,size(inter)
               dsps = dsp_line_from_dsp_px(inter(j), px(i), px(i))
               px_out = px_three_phase_from_inter(&
                  inter(j), px(i), px(i) &
               )
            end do
         end if
      end do
   end subroutine
end program main
