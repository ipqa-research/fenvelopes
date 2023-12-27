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
   type(injelope) :: px_bub, px_dew, px_hpl !! Shared 2ph-Px envelopes

   real(pr) :: alpha=0.0

   ! Setup everything
   call setup

   call get_z(alpha, z)
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
      call k_wilson_bubble(z, t_0=230.0_pr, p_end=0.5_pr, t=t, p=p, k=k)
      call envelope2( &
         1, nc, z, T, P, k, &
         n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, &
         pt_bub &
      )
      ! ========================================================================

      ! ========================================================================
      !  Dew/AOP envelopes
      ! ------------------------------------------------------------------------
      t = 300
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
      t = 700.0_pr
      t = pt_bub%t(maxloc(pt_bub%p, dim=1))
      p = maxval([pt_bub%p, pt_dew%p])*1.5_pr

      p = 400

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
      print *, style_bold // pt_case // style_reset
      ! ========================================================================

      three_phase: block
         use envelopes, only: pt_three_phase_from_intersection
         allocate(pt_bub_3(size(intersections)), pt_dew_3(size(intersections)))
         select case(pt_case)
         case("2_DEW_BUB")
            dsp_line: block
               use dsp_lines, only: injelope, dsp_line_from_dsp
               type(injelope):: dsps(2)
               integer :: i
               do i=1,size(intersections)
                  dsps = dsp_line_from_dsp(intersections(i), pt_dew, pt_bub)
               end do
            end block dsp_line
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
         case("2_HPL_BUB")
            call pt_three_phase_from_intersection(&
                  pt_hpl, pt_bub, [intersections(1)], &
                  pt_bub_3, pt_dew_3 &
            )
            call pt_three_phase_from_intersection(&
                  pt_hpl, pt_bub, [intersections(2)], &
                  pt_bub_3, pt_dew_3 &
            )
            dsp_line_2hpl_bub: block
               use dsp_lines, only: injelope, dsp_line_from_dsp
               type(injelope):: dsps(2)
               integer :: i
               do i=1,size(intersections)
                  dsps = dsp_line_from_dsp(intersections(i), pt_hpl, pt_bub, alpha0=alpha)
               end do
            end block dsp_line_2hpl_bub
         case("1_HPL_DEW")
            dsp_line_hpl: block
               use dsp_lines, only: injelope, dsp_line_from_dsp
               type(injelope):: dsps(2)
               integer :: i
               do i=1,size(intersections)
                  dsps = dsp_line_from_dsp(intersections(i), pt_hpl, pt_dew, alpha0=alpha)
               end do
            end block dsp_line_hpl
            call pt_three_phase_from_intersection(&
                  pt_hpl, pt_dew, intersections, &
                  pt_bub_3, pt_dew_3 &
            )
         case("1_HPL_BUB")
            dsp_line_hpl_bub: block
               use dsp_lines, only: injelope, dsp_line_from_dsp
               type(injelope):: dsps(2)
               integer :: i
               do i=1,size(intersections)
                  dsps = dsp_line_from_dsp(intersections(i), pt_hpl, pt_bub)
               end do
            end block dsp_line_hpl_bub
            call pt_three_phase_from_intersection(&
                  pt_hpl, pt_bub, intersections, &
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
