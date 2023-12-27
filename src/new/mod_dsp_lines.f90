module dsp_lines
   !! Module to calculate DSP lines
   use constants, only: pr, R
   use dtypes, only: envelope, critical_point
   use inj_envelopes, only: update_spec, injection_case, z_0, z_injection, &
                            get_z, injelope
   use linalg, only: solve_system, interpol, full_newton
   use progress_bar_module, only: progress_bar

   implicit none
   ! ===========================================================================
   !  Parameters
   ! ---------------------------------------------------------------------------
   integer :: env_number = 0 !! Number of calculated envelope
   integer :: max_iters = 10 !! Maximum number of iterations for a newton step
   integer, parameter :: max_points = 1000 !! Maximum number of points for each envelope
   character(len=255) :: FE_LOG

   real(pr) :: del_S_multiplier = 1.0_pr
   real(pr) :: max_da=0.1_pr
   real(pr) :: max_dp=1000.0_pr
   real(pr) :: max_dT=100.0_pr
   ! ===========================================================================
contains
   ! ===========================================================================
   ! Three-phases
   ! ---------------------------------------------------------------------------
   subroutine dsp_line_F(Xvars, ns, S, F, dF)
      !! Function to solve at each point of a DSP line.
      !!
      !! The vector of variables X corresponds to:
      !! \( X = [lnKx_i, lnKy_i lnP, \alpha, T] \)
      !!
      !! While the equations are:
      !!
      !! \( F = [
      !!        lnKx_i - ln \phi_i(z, P, T) + ln \phi_i(x, P, T),
      !!        lnKy_i - ln \phi_i(z, P, T) + ln \phi_i(y, P, T),
      !!        \sum_{i=1}^N (x_i - z_i),
      !!        \sum_{i=1}^N (y_i - z_i),
      !!        X_{ns} - S
      !! ] \)
      use legacy_ar_models, only: TERMO
      use iso_fortran_env, only: error_unit
      real(pr), intent(in)  :: Xvars(:) !! Vector of variables
      integer, intent(in)  :: ns   !! Number of specification
      real(pr), intent(in)  :: S    !! Specification value
      real(pr), intent(out) :: F(size(Xvars)) !! Vector of functions valuated
      real(pr), intent(out) :: df(size(Xvars), size(Xvars)) !! Jacobian matrix

      ! Xvars variables
      real(pr) :: Kx((Size(Xvars)-3)/2)
      real(pr) :: Ky((Size(Xvars)-3)/2)
      real(pr) :: P
      real(pr) :: alpha
      real(pr) :: T

      ! Incipient phase 1 variables
      real(pr) :: Vx
      real(pr), dimension((Size(Xvars)-3)/2) :: x, lnfug_x, dlnphi_dt_x, dlnphi_dp_x
      real(pr), dimension((Size(Xvars)-3)/2, (Size(Xvars)-3)/2) :: dlnphi_dn_x

      ! Incipient phase 2 variables
      real(pr) :: Vy
      real(pr), dimension((Size(Xvars)-3)/2) :: y, lnfug_y, dlnphi_dt_y, dlnphi_dp_y
      real(pr), dimension((Size(Xvars)-3)/2, (Size(Xvars)-3)/2) :: dlnphi_dn_y

      ! Main phase variables
      real(pr) :: Vz
      real(pr), dimension((Size(Xvars)-3)/2) :: z, lnfug_z, dlnphi_dt_z, dlnphi_dp_z
      real(pr), dimension((Size(Xvars)-3)/2, (Size(Xvars)-3)/2) :: dlnphi_dn_z

      ! Derivative of z wrt alpha
      real(pr) :: dzda((Size(Xvars)-3)/2)

      real(pr) :: dzdKx((Size(Xvars)-3)/2), dxdKx((Size(Xvars)-3)/2), dydKx((Size(Xvars)-3)/2)
      real(pr) :: dzdKy((Size(Xvars)-3)/2), dxdKy((Size(Xvars)-3)/2), dydKy((Size(Xvars)-3)/2)

      integer :: i, j, n

      n = (Size(Xvars)-3)/2

      ! Setting variables
      Kx = exp(Xvars(1:n))
      Ky = exp(Xvars(n + 1:2*n))

      P = exp(Xvars(2*n + 1))
      alpha = Xvars(2*n + 2)
      T = exp(Xvars(2*n + 3))

      call get_z(alpha, z, dzda)
      if (any(z < 0)) return

      x = Kx * z
      y = Ky * z

      call TERMO( &
         n, 0, 4, T, P, x, Vx, lnfug_x, dlnphi_dp_x, dlnphi_dt_x, dlnphi_dn_x &
      )
      call TERMO( &
         n, 0, 4, T, P, y, Vy, lnfug_y, dlnphi_dp_y, dlnphi_dt_y, dlnphi_dn_y &
      )
      call TERMO( &
         n, 0, 4, T, P, z, Vz, lnfug_z, dlnphi_dp_z, dlnphi_dt_z, dlnphi_dn_z &
      )

      F(1:n)       = log(Kx) + lnfug_x - lnfug_z
      F(n + 1:2*n) = log(Ky) + lnfug_y - lnfug_z

      F(2*n + 1) = sum(x - z)
      F(2*n + 2) = sum(y - z)
      F(2*n + 3) = Xvars(ns) - S
      df = 0

      do i=1,n
         do j=1,n
            df(i, j)   = x(j) * dlnphi_dn_x(i, j) 
            df(i+n, j+n) = y(j) * dlnphi_dn_y(i, j) 
         end do
         df(i, i) = df(i, i) + 1
         df(i+n,i+n) = df(i+n, i+n) + 1

         ! Derivatives wrt alpha
         df(i,   2*n+2) = sum(dlnphi_dn_x(i, :) * Kx * dzda - dlnphi_dn_z(i, :) * dzda)
         df(i+n, 2*n+2) = sum(dlnphi_dn_y(i, :) * Ky * dzda - dlnphi_dn_z(i, :) * dzda)
      end do

      df(:n, 2*n+1) = P * (dlnphi_dp_x - dlnphi_dp_z)
      df(n+1:2*n, 2*n+1) = P * (dlnphi_dp_y - dlnphi_dp_z)
      
      df(:n, 2*n+3) = T * (dlnphi_dt_x - dlnphi_dt_z)
      df(n+1:2*n, 2*n+3) = T * (dlnphi_dt_y - dlnphi_dt_z)

      df(2*n+1, :n) = x
      df(2*n+2, n+1:2*n) = y

      df(2*n+1, 2*n+2) = sum(Kx * dzda - dzda)
      df(2*n+2, 2*n+2) = sum(Ky * dzda - dzda)

      df(2*n+3, ns) = 1
   end subroutine

   subroutine dsp_line(X0, spec_number, del_S0, envels)
      use constants, only: ouput_path
      use io, only: str
      !! Subroutine to calculate Px phase envelopes via continuation method.
      !! Three phases version.
      real(pr), intent(in) :: X0(:) !! Vector of variables
      integer, intent(in) :: spec_number !! Number of specification
      real(pr), intent(in) :: del_S0 !! \(\Delta S_0\)
      type(injelope), intent(out) :: envels !! Calculated envelopes

      type(critical_point), allocatable :: cps(:)

      real(pr) :: X(size(X0))
      integer :: ns
      real(pr) :: S
      real(pr) :: XS(max_points, size(X0))
      real(pr) :: del_S

      real(pr) :: F(size(X0)), dF(size(X0), size(X0)), dXdS(size(X0))

      integer :: point, iters, n
      integer :: i
      integer :: funit_output
      character(len=254) :: fname_env

      allocate (cps(0))
      X = X0
      n = (size(X0) - 3)/2
      ns = spec_number
      S = X(ns)
      del_S = del_S0

      ! ======================================================================
      !  Output file
      ! ----------------------------------------------------------------------
      env_number = env_number + 1

      write (fname_env, *) env_number
      fname_env = "env-3ph-DSP"//"_"//trim(adjustl(fname_env))
      fname_env = trim(adjustl(ouput_path))//trim(fname_env)//".dat"

      open (newunit=funit_output, file=fname_env)
      write (funit_output, *) "#", env_number
      write (funit_output, *) "STAT", " iters", " ns", " alpha", " P", &
         " T", (" lnKx"//str(i), i=1,n), (" lnKy"//str(i), i=1,n)
      write (funit_output, *) "X0", iters, ns, &
         X(2*n + 2), exp(X(2*n + 1)), exp(X(2*n + 3)), X(:2*n)
      ! ======================================================================

      enveloop: do point = 1, max_points
         call progress_bar(point, max_points, advance=.false.)
         call full_newton(dsp_line_F, iters, X, ns, S, max_iters, F, dF)
         if (iters >= max_iters) then
            call progress_bar(point, max_points, advance=.true.)
            print *, "Breaking: Above max iterations"
            exit enveloop
         end if

         write (funit_output, *) "SOL", iters, ns, &
            X(2*n + 2), exp(X(2*n + 1)), exp(X(2*n + 3)), X(:2*n)
         XS(point, :) = X

         call update_spec(X, ns, del_S, dF, dXdS)
         call fix_step(iters, X, ns, S, del_S, dXdS)

         X = X + dXdS*del_S
         S = X(ns)

         if (any(break_conditions_dsp_line(X, ns, S, del_S)) .and. point > 10) then
            call progress_bar(point, max_points, .true.)
            print *, "Breaking: ", break_conditions_dsp_line(X, ns, S, del_S)
            exit enveloop
         end if

      end do enveloop

      write (funit_output, *) ""
      write (funit_output, *) ""
      write (funit_output, *) "#critical"
      if (size(cps) > 0) then
         do i = 1, size(cps)
            write (funit_output, *) cps(i)%alpha, cps(i)%p
         end do
      else
         write (funit_output, *) "NaN NaN"
      end if

      close (funit_output)
      envels%z = z_0
      envels%z_inj = z_injection
      envels%logk = XS(:point-1, :n)
      envels%alpha = XS(:point-1, n + 2)
      envels%p = exp(XS(:point-1, n + 1))
      envels%critical_points = cps
   end subroutine

   subroutine fix_step(iters, X, ns, S, del_S, dXdS)
         integer, intent(in) :: iters
         real(pr), intent(in) :: X(:)
         integer, intent(in) :: ns
         real(pr), intent(in) :: S
         real(pr), intent(in out) :: del_S
         real(pr), intent(in out) :: dXdS(size(X))
         real(pr) :: Xnew(size(X))
         real(pr) :: da, dP, dT
         integer :: n

         n = (size(X) - 3)/2

         del_S = sign(del_S_multiplier, del_S)*minval([ &
                                            max(abs(sqrt(X(ns))/10), 0.1_pr), &
                                            abs(del_S)*3/iters &
                                            ] &
                                            )

         if (injection_case == "dilution") del_S = 50*del_S

         Xnew = X + dXdS*del_S
         da = (Xnew(2*n + 2)) - (X(2*n + 2))
         dP = exp(Xnew(2*n + 1)) - exp(X(2*n + 1))
         dT = exp(Xnew(2*n + 3)) - exp(X(2*n + 3))

         do while (abs(da) > max_da .or. abs(dP) > max_dP .or. abs(dT) > max_dT)
            dXdS = dXdS/2.0_pr
            Xnew = X + dXdS*del_S
            da = (Xnew(2*n + 2)) - (X(2*n + 2))
            dP = exp(Xnew(2*n + 1)) - exp(X(2*n + 1))
            dT = exp(Xnew(2*n + 3)) - exp(X(2*n + 3))
         end do
   end subroutine
   
   function break_conditions_dsp_line(X, ns, S, del_S)
      !! Set of conditions to break the tracing.
      real(pr) :: X(:) !! Variables vector
      integer :: ns !! Number of specification
      real(pr) :: S !! Value of specification
      real(pr) :: del_S

      integer :: n
      real(pr) :: p, alpha
      logical, allocatable :: break_conditions_dsp_line(:)

      n = (size(X) - 3)/2
      p = exp(X(2*n + 1))
      alpha = X(2*n + 2)

      break_conditions_dsp_line = [ .false. ]
   end function
   ! ===========================================================================
   
   ! ===========================================================================
   ! Initialization procedures
   ! ---------------------------------------------------------------------------
   function dsp_line_from_dsp(&
         dsp, pt_1, pt_2, del_S0, alpha0 &
      ) result(envels)
      !! Calculate a DSP line from an already found double saturation point
      !! (`dsp`) between two PT lines (`pt_1`and `pt_2`).
      !! 
      use legacy_ar_models, only: nc
      use stdlib_optval, only: optval
      use linalg, only: point, interpol
      type(point), intent(in) :: dsp !! DSP point
      type(envelope), intent(in) :: pt_1 !! PT envelope
      type(envelope), intent(in) :: pt_2 !! PT envelope
      type(injelope) :: envels(2) !! Output envelopes
      real(pr), optional :: del_S0 !! Initial step size, defaults to 0.01
      real(pr), optional :: alpha0 !! Initial \(\alpha)\, defaults to 0

      integer :: i, j

      real(pr) :: lnKx(nc), lnKy(nc), t, p, X(2*nc+3), alpha
      real(pr) :: del_S
      real(pr) :: z(nc), dzda(nc)
      integer :: ns

      del_S = optval(del_S0, 0.01_pr)
      alpha = optval(alpha0, 0._pr)

      i = dsp%i
      j = dsp%j

      t = dsp%x
      p = dsp%y

      lnKx = interpol( &
               pt_1%t(i),   pt_1%t(i + 1), &
               pt_1%logk(i, :), pt_1%logk(i + 1, :), &
               t &
            )

      lnKy = interpol( &
               pt_2%t(j), pt_2%t(j + 1), &
               pt_2%logk(j, :), pt_2%logk(j + 1, :), &
               t &
            )

      call get_z(alpha, z, dzda)

      ns = 2*nc + 2

      X = [lnKx, lnKy, log(p), alpha, log(t)]
      call dsp_line(X, ns, del_S, envels(1))
      
      X = [lnKx, lnKy, log(p), alpha, log(t)]
      call dsp_line(X, ns, -del_S, envels(2))
      ! ==================================================================
   end function
   
   function dsp_line_from_dsp_px(&
         dsp, px_1, px_2, del_S0, alpha0 &
      ) result(envels)
      !! Calculate a DSP line from an already found double saturation point
      !! (`dsp`) between two P\(\alpha\) lines (`px_1`and `px_2`).
      !! 
      use legacy_ar_models, only: nc
      use stdlib_optval, only: optval
      use linalg, only: point, interpol
      use inj_envelopes, only: t_inj => T
      type(point), intent(in) :: dsp !! Double Saturation Point
      type(injelope), intent(in) :: px_1 !! P\(\alpha\) line
      type(injelope), intent(in) :: px_2 !! P\(\alpha\) line
      type(injelope) :: envels(2) !! DSP lines
      real(pr), optional :: del_S0 !! Initial step size, defaults to 0.01
      real(pr), optional :: alpha0 !! Initial \(\alpha\), defaults to 0

      integer :: i, j

      real(pr) :: lnKx(nc), lnKy(nc), t, p, X(2*nc+3), alpha
      real(pr) :: del_S
      real(pr) :: z(nc), dzda(nc)
      integer :: ns

      del_S = optval(del_S0, 0.01_pr)
      alpha = optval(alpha0, 0._pr)

      i = dsp%i
      j = dsp%j

      alpha = dsp%x
      p = dsp%y
      t = t_inj

      lnKx = interpol( &
               px_1%alpha(i),   px_1%alpha(i + 1), &
               px_1%logk(i, :), px_1%logk(i + 1, :), &
               alpha &
            )

      lnKy = interpol( &
               px_2%alpha(j), px_2%alpha(j + 1), &
               px_2%logk(j, :), px_2%logk(j + 1, :), &
               alpha &
            )

      call get_z(alpha, z, dzda)

      ns = 2*nc + 3
      print *, "DSP Line", alpha, t, p

      X = [lnKx, lnKy, log(p), alpha, log(t)]
      call dsp_line(X, ns, del_S, envels(1))
      
      X = [lnKx, lnKy, log(p), alpha, log(t)]
      call dsp_line(X, ns, -del_S, envels(2))
      ! ==================================================================
   end function
   ! ===========================================================================
end module