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
   integer :: max_iters = 50 !! Maximum number of iterations for a newton step
   integer, parameter :: max_points = 1000 !! Maximum number of points for each envelope
   character(len=255) :: FE_LOG

   ! Two-phase settings
   real(pr) :: del_S_multiplier = 2.0_pr
   real(pr) :: max_dp=50.0_pr
   real(pr) :: max_dalpha=0.01_pr
   real(pr) :: critical_multiplier = 2.0_pr

   ! Three phase parameters
   real(pr) :: del_S_multiplier_three_phase = 1.7_pr
   real(pr) :: critical_fact = 3.0_pr
   ! ===========================================================================
contains
   subroutine from_nml(filepath)
      use legacy_ar_models, only: nc
      character(len=*), intent(in) :: filepath
      integer :: funit

      namelist /nml_px/ T, z_0, z_injection, injection_case

      allocate (z_0(nc), z_injection(nc))

      open (newunit=funit, file=filepath)
         read (funit, nml=nml_px)
      close (funit)

      z_injection = z_injection/sum(z_injection)
   end subroutine

   ! ===========================================================================
   ! Three-phases
   ! ---------------------------------------------------------------------------
   subroutine F_injection_three_phases(Xvars, ns, S, F, dF)
      !! Function to solve at each point of a three phase envelope.
      !!
      !! The vector of variables X corresponds to:
      !! \( X = [lnKx_i, lnKy_i lnP, \alpha, T] \)
      !!
      !! While the equations are:
      !!
      !! \( F = [
      !!        lnKx_i - ln \phi_i(x, P, T) + ln \phi_i(w, P, T),
      !!        lnKy_i - ln \phi_i(y, P, T) + ln \phi_i(w, P, T),
      !!        \sum_{i=1}^N (w_i) - 1,
      !!        \sum_{i=1}^N (x_i - y_i),
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
      real(pr) :: z((Size(Xvars)-3)/2)
      real(pr) :: Kx((Size(Xvars)-3)/2)
      real(pr) :: Ky((Size(Xvars)-3)/2)
      real(pr) :: P
      real(pr) :: alpha
      real(pr) :: T

      real(pr) :: beta=1

      ! Main phase 1 variables
      real(pr) :: Vx
      real(pr), dimension((Size(Xvars)-3)/2) :: x, lnfug_x, dlnphi_dt_x, dlnphi_dp_x
      real(pr), dimension((Size(Xvars)-3)/2, (Size(Xvars)-3)/2) :: dlnphi_dn_x

      ! Main phase 2 variables
      real(pr) :: Vy
      real(pr), dimension((Size(Xvars)-3)/2) :: y, lnfug_y, dlnphi_dt_y, dlnphi_dp_y
      real(pr), dimension((Size(Xvars)-3)/2, (Size(Xvars)-3)/2) :: dlnphi_dn_y

      ! Incipient phase variables
      real(pr) :: Vw
      real(pr), dimension((Size(Xvars)-3)/2) :: w, lnfug_w, dlnphi_dt_w, dlnphi_dp_w
      real(pr), dimension((Size(Xvars)-3)/2, (Size(Xvars)-3)/2) :: dlnphi_dn_w

      ! Derivative of z wrt alpha
      real(pr) :: dzda((Size(Xvars)-3)/2), dwda((Size(Xvars)-3)/2)

      ! Derivative of w wrt beta
      real(pr) :: dwdb((Size(Xvars)-3)/2)

      real(pr) :: dwdKx((Size(Xvars)-3)/2), dxdKx((Size(Xvars)-3)/2), dydKx((Size(Xvars)-3)/2)
      real(pr) :: dwdKy((Size(Xvars)-3)/2), dxdKy((Size(Xvars)-3)/2), dydKy((Size(Xvars)-3)/2)

      integer :: i, j, n

      n = (Size(Xvars)-3)/2

      ! Setting variables
      Kx = exp(Xvars(1:n))
      Ky = exp(Xvars(n + 1:2*n))
      P = exp(Xvars(2*n + 1))
      alpha = Xvars(2*n + 2)
      T = Xvars(2*n + 3)

      call get_z(alpha, z, dzda)

      w = z/(beta*Ky + (1 - beta)*Kx)
      x = w*Kx
      y = w*Ky

      call TERMO( &
         n, 0, 4, T, P, x, Vx, lnfug_x, dlnphi_dp_x, dlnphi_dt_x, dlnphi_dn_x &
         )
      call TERMO( &
         n, 0, 4, T, P, y, Vy, lnfug_y, dlnphi_dp_y, dlnphi_dt_y, dlnphi_dn_y &
         )
      call TERMO( &
         n, 0, 4, T, P, w, Vw, lnfug_w, dlnphi_dp_w, dlnphi_dt_w, dlnphi_dn_w &
         )

      F(1:n) = Xvars(1:n) + lnfug_x - lnfug_w
      F(n + 1:2*n) = Xvars(n + 1:2*n) + lnfug_y - lnfug_w

      F(2*n + 1) = sum(w) - 1
      F(2*n + 2) = sum(x - y)
      F(2*n + 3) = Xvars(ns) - S

      df = 0
      dwda = 1.0_pr/(beta*Ky + (1 - beta)*Kx)*dzda
      dwdb = z*(Kx - Ky)/((1 - beta)*Kx + beta*Ky)**2

      dwdKx = -z*(1 - beta)/(Ky*beta + (1 - beta)*Kx)**2
      dxdKx = Kx*dwdKx + w
      dydKx = Ky*dwdKx

      dwdKy = -z*(beta)/(Ky*beta + (1 - beta)*Kx)**2
      dxdKy = Kx*dwdKy
      dydKy = Ky*dwdKy + w

      do i = 1, n
         do j = 1, n
            df(i, j) = Kx(j)*(dlnphi_dn_x(i, j)*dxdKx(j) &
                              - dlnphi_dn_w(i, j)*dwdKx(j))
            df(i + n, j) = Kx(j)*(dlnphi_dn_y(i, j)*dydKx(j) &
                                  - dlnphi_dn_w(i, j)*dwdKx(j))

            df(i, j + n) = Ky(j)*(dlnphi_dn_x(i, j)*dxdKy(j) &
                                  - dlnphi_dn_w(i, j)*dwdKy(j))
            df(i + n, j + n) = Ky(j)*(dlnphi_dn_y(i, j)*dydKy(j) &
                                      - dlnphi_dn_w(i, j)*dwdKy(j))
         end do

         df(i, i) = df(i, i) + 1
         df(i + n, i + n) = df(i + n, i + n) + 1
         df(i, 2*n + 2) = sum( &
                          Kx*dlnphi_dn_x(i, :)*dwda - dlnphi_dn_w(i, :)*dwda &
                          )
         df(i + n, 2*n + 2) = sum(Ky*dlnphi_dn_y(i, :)*dwda &
                                  - dlnphi_dn_w(i, :)*dwda)
         df(i, 2*n + 3) = sum(Kx*dlnphi_dn_x(i, :)*dwdb &
                              - dlnphi_dn_w(i, :)*dwdb)
         df(i + n, 2*n + 3) = sum(Ky*dlnphi_dn_y(i, :)*dwdb &
                                  - dlnphi_dn_w(i, :)*dwdb)
         df(2*n + 1, i) = Kx(i)*dwdKx(i)
         df(2*n + 1, i + n) = Ky(i)*dwdKy(i)

         df(2*n + 2, i) = Kx(i)*dxdKx(i) - Kx(i)*dydKx(i)
         df(2*n + 2, i + n) = Ky(i)*dxdKy(i) - Ky(i)*dydKy(i)
      end do

      ! Derivatives wrt P
      df(:n, 2*n + 1) = P*(dlnphi_dp_x - dlnphi_dp_w)
      df(n + 1:2*n, 2*n + 1) = P*(dlnphi_dp_y - dlnphi_dp_w)

      ! Derivatives wrt alpha
      df(2*n + 1, 2*n + 2) = sum(dwda)
      df(2*n + 2, 2*n + 2) = sum(Kx*dwda - Ky*dwda)

      ! Derivatives wrt T
      df(:n, 2*n + 3) = T*(dlnphi_dp_x - dlnphi_dp_w)
      df(n + 1:2*n, 2*n + 1) = T*(dlnphi_dp_y - dlnphi_dp_w)

      ! Derivatives wrt Xs
      df(2*n + 3, :) = 0
      df(2*n + 3, ns) = 1
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
      write (funit_output, *) "#", T
      write (funit_output, *) "STAT", " iters", " ns", " alpha", " P", &
         " T", (" lnKx"//str(i), i=1,n), (" lnKy"//str(i), i=1,n)
      write (funit_output, *) "X0", iters, ns, X(2*n + 2), exp(X(2*n + 1)), &
         X(2*n + 3), X(:2*n)
      ! ======================================================================

      enveloop: do point = 1, max_points
         call progress_bar(point, max_points, advance=.false.)
         call full_newton(F_injection_three_phases, iters, X, ns, S, max_iters, F, dF)
         if (iters >= max_iters) then
            call progress_bar(point, max_points, advance=.true.)
            print *, "Breaking: Above max iterations"
            exit enveloop
         end if

         write (funit_output, *) "SOL", iters, ns, X(2*n + 2), &
            exp(X(2*n + 1)), X(2*n + 3), X(:2*n)
         XS(point, :) = X

         call update_spec(X, ns, del_S, dF, dXdS)

         fix_step: block
         end block fix_step

         detect_critical: block
            real(pr) :: K((size(X0) - 3)/2), Knew((size(X0) - 3)/2), &
                        Xnew(size(X0)), fact
            real(pr) :: pc, alpha_c, dS_c, dXdS_in(size(X0))
            integer :: max_changing, i
            fact = critical_fact

            loop: do i = 0, 1
               Xnew = X + fact*dXdS*del_S

               K = X(i*n + 1:(i + 1)*n)
               Knew = Xnew(i*n + 1:(i + 1)*n)

               max_changing = maxloc(abs(Knew - K), dim=1)

               if (all(K*Knew < 0)) then
                  dS_c = ( &
                         -k(max_changing)*(Xnew(ns) - X(ns)) &
                         /(Knew(max_changing) - K(max_changing)) &
                         )

                  Xnew = X + dXdS*dS_c
                  alpha_c = Xnew(2*n + 2)
                  pc = exp(Xnew(2*n + 1))
                  cps = [cps, critical_point(t, pc, alpha_c)]

                  ! del_S = dS_c + del_S ! * fact
                  ! del_S = del_S * fact
                  del_S = dS_c - sign(1.7_pr, dS_c)*dS_c

                  write (funit_output, *) ""
                  write (funit_output, *) ""
                  exit loop
               end if
            end do loop
         end block detect_critical

         if (x(2*n + 3) > 1 .or. (x(2*n+3) < 0)) then
            call progress_bar(point, max_points, .true.)
            print *, "Breaking: positive 𝜷"
            exit enveloop
         end if

         X = X + dXdS*del_S
         S = X(ns)
         if (any(break_conditions_three_phases(X, ns, S, del_S)) .and. point > 10) then
            call progress_bar(point, max_points, .true.)
            print *, "Breaking: ", break_conditions_three_phases(X, ns, S, del_S)
            exit enveloop
         end if
      end do enveloop

      ! point = point - 1

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
      envels%logk = XS(:point, :n)
      envels%alpha = XS(:point, n + 2)
      envels%p = exp(XS(:point, n + 1))
      envels%critical_points = cps
   end subroutine

   subroutine update_spec(X, ns, S, del_S, dXdS)
         real(pr), intent(in) :: X(:)
         integer, intent(in) :: ns
         real(pr), intent(in) :: S
         real(pr), intent(in) :: del_S
         real(pr), intent(in out) :: dXdS(size(X))
         real(pr) :: Xnew(size(X))
         real(pr) :: dP, dT

         del_S = sign(del_S_multiplier_three_phase, del_S)*minval([ &
                                            max(abs(X(ns)/10), 0.1_pr), &
                                            abs(del_S)*3/iters &
                                            ] &
                                            )

         if (injection_case == "dilution") del_S = 50*del_S

         Xnew = X + dXdS*del_S
         dP = exp(Xnew(2*n + 1)) - exp(X(2*n + 1))
         dalpha = exp(Xnew(2*n + 3)) - exp(X(2*n + 3))

         do while (abs(dP) > max_dp .or. abs(dalpha) > max_dalpha)
            dXdS = dXdS/2.0_pr

            Xnew = X + dXdS*del_S
            dP = exp(Xnew(2*n + 1)) - exp(X(2*n + 1))
            dalpha = Xnew(2*n + 2) - X(2*n + 2)
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
      logical, allocatable :: break_conditions_three_phases(:)

      n = (size(X) - 3)/2
      p = exp(X(2*n + 1))
      alpha = X(2*n + 2)

      break_conditions_three_phases = [ &
                                       ! p < 0.001_pr .or. p > 5000, &
                                       abs(del_S) < 1e-5 &
                                      ]
   end function
   ! ===========================================================================
   
   subroutine get_z(alpha, z, dzda)
      !! Calculate the fluid composition based on an amount of addition
      !! of second fluid.
      !!
      !! The injection can be considered as two kinds of injection:
      !! - Displacement: \( z = \alpha z_i + (1-\alpha) z_0 \)
      !! - Addition:  \( z = \frac{\alpha z_i + (1-\alpha) z_0}{\sum_{i=1}^N \alpha z_i + (1-\alpha) z_0} \)
      real(pr), intent(in)  :: alpha !! Addition percentaje \( \alpha \)
      real(pr), intent(out) :: z(size(z_0)) !! New composition
      real(pr), optional, intent(out) :: dzda(size(z_0)) !! Derivative wrt \(\alpha\)

      select case (injection_case)
      case ("displace")
         z = (z_injection*alpha + (1.0_pr - alpha)*z_0)
         if (present(dzda)) dzda = z_injection - z_0
      case ("dilute")
         z = (z_injection*alpha + z_0)/sum(z_injection*alpha + z_0)
         if (present(dzda)) dzda = -(alpha*z_injection + z_0) &
                *sum(z_injection)/sum(alpha*z_injection + z_0)**2 &
                + z_injection/sum(alpha*z_injection + z_0)
      case default
         z = (z_injection*alpha + (1.0_pr - alpha)*z_0)
         if (present(dzda)) dzda = z_injection - z_0
      end select
   end subroutine

   function get_case(dew_envel, bub_envel) result(n_case)
      type(injelope), intent(in) :: dew_envel
      type(injelope), intent(in) :: bub_envel
      integer :: n_case
   end function

   function remove_duplicates(envels) result(clean_envels)
      !! From a set of envelopes check if they are the same line
      class(injelope) :: envels(:)
      type(injelope), allocatable :: clean_envels(:)

      if (size(envels) /= 1) then
      else
         clean_envels = envels
      end if
   end function

   function same_line(env1, env2)
      !! 
      class(injelope), intent(in) :: env1, env2
      logical :: same_line
   end function

   ! ===========================================================================
   ! Initialization procedures
   ! ---------------------------------------------------------------------------
   function px_three_phase_from_inter(&
         inter, px_1, px_2, del_S0, beta0 &
         ) result(envels)
      use legacy_ar_models, only: nc
      use stdlib_optval, only: optval
      use linalg, only: point, interpol
      type(point), intent(in) :: inter
      type(injelope), intent(in) :: px_1, px_2
      type(injelope) :: envels(2)
      real(pr), optional :: del_S0
      real(pr), optional :: beta0

      integer :: i, j

      real(pr) :: lnKx(nc), lnKy(nc), alpha, p, beta, X(2*nc+3)
      real(pr) :: phase_y(nc), phase_x(nc)
      real(pr) :: del_S
      real(pr) :: z(nc), dzda(nc)
      integer :: ns

      del_S = optval(del_S0, -0.05_pr)
      beta = optval(beta0, 0.99_pr)

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
      call injection_envelope_three_phase(X, ns, del_S, envels(1))
      ! ==================================================================

      ! ==================================================================
      !  Line with incipient phase liquid
      ! ------------------------------------------------------------------
      print *, "Three Phase: Liquid"
      lnKx = log(phase_y/phase_x)
      lnKy = log(z/phase_x)
      X = [lnKx, lnKy, log(p), alpha, beta]
      call injection_envelope_three_phase(X, ns, del_S, envels(2))
   end function
   ! ===========================================================================
end module
