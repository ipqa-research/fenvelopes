program main
   use dtypes, only: envelope
   use envelopes, only: PTEnvel3
   use inj_envelopes, only: injelope, get_z
   use constants, only: pr, ouput_path
   use legacy_ar_models, only: nc, z
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
   type(injelope), allocatable :: px_bub(:), px_dew(:), px_hpl !! Shared 2ph-Px envelopes

   ! real(pr) :: alpha=0.95
   real(pr) :: alpha=0.0
   real(pr) :: pt_bub_t0 = 180
   real(pr) :: pt_dew_t0 = 180

   ! Setup everything
   call setup

   call get_z(alpha, z)
   ! PT Envelopes
   call cpu_time(st)
   call pt_envelopes
   call cpu_time(et)
   print *, "PT: ", (et - st)*1000, "cpu ms"
      
   call exit

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
      call cli%add( &
         switch="--injection", &
         switch_ab="-px", &
         help="Trace Px lines", &
         error=cli_error, &
         required=.false.)
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
      call k_wilson_bubble(z, t_0=pt_bub_t0, p_end=0.5_pr, t=t, p=p, k=k)
      call envelope2( &
         1, nc, z, T, P, k, &
         n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, &
         pt_bub &
      )
      ! ========================================================================

      ! ========================================================================
      !  Dew/AOP envelopes
      ! ------------------------------------------------------------------------
      t = pt_dew_t0
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
      ! ========================================================================

      ! ========================================================================
      !  HPLL Envelope
      ! ------------------------------------------------------------------------
      ! t = 700.0_pr
      ! p = maxval([pt_bub%p, pt_dew%p])*1.5_pr

      p = 900.0_pr
      t = pt_dew%t(maxloc(pt_dew%p, dim=1))

      call find_hpl(t, p, k)
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
               use legacy_ar_models, only: nc, z
               use legacy_thermo_properties, only: termo
               use envelopes, only: pt_envelope_three_phase
               real(pr) :: p, v, t
               real(pr) :: lnphix(nc), lnphiy(nc), &
                           phase_x(nc), phase_y(nc), beta, x(2*nc+3), lnKx(nc), lnKy(nc)

               beta = z(nc)
               phase_y = 0
               phase_y(nc) = 1
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
      type(point), allocatable :: inter_dew_bub(:), self_inter_dew(:), self_inter_bub(:)

      real(pr) :: t_tol = 2

      print *, style_underline // "----------" // style_reset
      print *, style_underline // "Px Regions" // style_reset
      print *, style_underline // "----------" // style_reset

      ! ========================================================================
      !  Two phase envelopes
      ! ------------------------------------------------------------------------
      print *, red // "Running Bubble" // style_reset
      px_bub = px_two_phase_from_pt(t_inj, pt_bub, t_tol=5.0_pr)

      print *, blue // "Running Dew" // style_reset
      px_dew = px_two_phase_from_pt(t_inj, pt_dew, t_tol=5.0_pr)

      print *, blue // "Running HPLL" // style_reset
      px_hpl = px_hpl_line(0.99_pr, 1000.0_pr)
      ! ========================================================================

      ! ========================================================================
      !  Three phase regions
      ! ------------------------------------------------------------------------
      three_phase: block
         use dsp_lines, only: injelope, dsp_line_from_dsp_px
         integer :: i, j, k
         type(injelope) :: px_bub_3, px_dew_3, px_branch_3(2)
         type(injelope):: dsps(2)

         ! =====================================================================
         ! Intersections between lines
         ! ---------------------------------------------------------------------
         do i=1,size(px_dew)
            do j=1,size(px_bub)
               ! Go through each possible pair of envelopes to find DSPs
               inter_dew_bub  = intersection(&
                  px_dew(i)%alpha, px_dew(i)%p, &
                  px_bub(j)%alpha, px_bub(j)%p  &
               )
               do k=1,size(inter_dew_bub)
                  ! For each DSP found, calculate the two possible DSP lines
                  ! and the two three-phase branchs.
                  print *, "Intersection: ", inter_dew_bub(k)
                  dsps = dsp_line_from_dsp_px(&
                     inter_dew_bub(k), px_dew(i), px_bub(j) &
                  )
                  px_branch_3 = px_three_phase_from_inter(&
                     inter_dew_bub(k), px_dew(i), px_bub(j) &
                  )
               end do
            end do
         end do
         ! =====================================================================

         ! =====================================================================
         ! Self intersections loops
         ! ---------------------------------------------------------------------
         do i=1,size(px_bub)
            print *, "Checking Bub slef-intersections: ", i
            self_inter_bub = intersection(px_bub(i)%alpha, px_bub(i)%p)
            if (size(self_inter_bub) > 0) then
               do j=1,size(self_inter_bub)
                  dsps = dsp_line_from_dsp_px(self_inter_bub(j), px_bub(i), px_bub(i))
                  px_branch_3 = px_three_phase_from_inter(&
                     self_inter_bub(j), px_bub(i), px_bub(i) &
                  )
               end do
            end if
         end do
         
         do i=1,size(px_dew)
            print *, "Checking Dew self-intersections: ", i
            self_inter_dew = intersection(px_dew(i)%alpha, px_dew(i)%p)
            if (size(self_inter_dew) > 0) then
               do j=1,size(self_inter_dew)
                  dsps = dsp_line_from_dsp_px(self_inter_dew(j), px_dew(i), px_dew(i))
                  px_branch_3 = px_three_phase_from_inter(&
                     self_inter_dew(j), px_dew(i), px_dew(i) &
                  )
               end do
            end if
         end do
         ! =====================================================================
         
         ! Isolated lines coming from PT lines
         print *, "Isolated Bub"
         px_bub_3 = px_three_phase_from_pt(t_inj, pt_bub_3, t_tol)
         print *, "Isolated Dew"
         px_dew_3 = px_three_phase_from_pt(t_inj, pt_dew_3, t_tol)
      end block three_phase
      ! ========================================================================
   end subroutine
end program main
