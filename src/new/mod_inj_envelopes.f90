module inj_envelopes
   !! Module to calculate Px phase envelopes 
   use constants, only: pr, R
   use dtypes, only: envelope, critical_point
   use linalg, only: solve_system, interpol

   implicit none

   integer :: env_number = 0

   type, extends(envelope) :: injelope
      real(pr), allocatable :: alpha(:) !! Ammount of injected fluid
      real(pr), allocatable :: z_inj(:) !! Injected fluid composition
      real(pr), allocatable :: z_mix(:, :) !! Composition at each step
   end type

   integer :: max_iters = 500 !! Maximum number of iterations for a newton step
   integer, parameter :: max_points = 800 !! Maximum number of points for each envelope
   real(pr), allocatable :: z_0(:) !! Original fluid composition
   real(pr), allocatable :: z_injection(:) !! Injection fluid composition
   real(pr) :: T !! Temperature of injection
   real(pr) :: del_S = 0.1 !! Specificiation variation
   character(len=10) :: injection_case !! Kind of injection displace|dilute
contains

   subroutine from_nml(filepath)
       ! use system, only: nc
       use legacy_ar_models, only: nc
       character(len=*), intent(in) :: filepath
       integer :: funit

       namelist /nml_px/ T ,z_0, z_injection, injection_case

       allocate(z_0(nc), z_injection(nc))

       open(newunit=funit, file=filepath)
           read(funit, nml=nml_px)
       close(funit)
   end subroutine

   subroutine F_injection(X, ns, S, F, dF)
      !! Function to solve at each point of the phase envelope.
      !!
      !! The vector of variables X corresponds to:
      !! \( X = [lnK_i, lnP, \alpha] \)
      !!
      !! While the equations are:
      !!
      !! \( F = [lnK_i - ln \phi_i(y, P, T) + ln \phi_i(z, P, T), 
      !!         \sum_{i=1}^N, X_{ns} - S] \)
      !!
      !! The injection can be considered as two kinds of injection:
      !! - Displacement: \( z = \alpha z_i + (1-\alpha) z_0 \)
      !! - Addition:  \( z = \frac{\alpha z_i + (1-\alpha) z_0}{\sum_{i=1}^N \alpha z_i + (1-\alpha) z_0} \)
      !!
      use iso_fortran_env, only: error_unit
      use legacy_ar_models, only: TERMO
      real(pr), intent(in)  :: X(:) !! Vector of variables
      integer, intent(in)   :: ns !! Number of specification
      real(pr), intent(in)  :: S !! Specification value
      real(pr), intent(out) :: F(size(X)) !! Vector of functions valuated
      real(pr), intent(out) :: df(size(x), size(X)) !! Jacobian matrix

      ! X variables
      real(pr) :: K(size(X) - 2)
      real(pr) :: alpha
      real(pr) :: P

      ! Main phase variables
      real(pr) :: Vz
      real(pr), dimension(size(X)-2) :: z, lnfug_z, dlnphi_dt_z, dlnphi_dp_z
      real(pr), dimension(size(X)-2, size(X)-2) :: dlnphi_dn_z

      ! Incipient phase variables
      real(pr) :: Vy
      real(pr), dimension(size(X)-2) :: y, lnfug_y, dlnphi_dt_y, dlnphi_dp_y
      real(pr), dimension(size(X)-2, size(X)-2) :: dlnphi_dn_y

      ! Derivative of z wrt alpha
      real(pr) :: dzda(size(X)-2)

      integer :: i, j, n

      n = size(X) - 2
      K = exp(X(1:n))
      P = exp(X(n+1))
      alpha = X(n+2)

      select case(injection_case)
          case("displace")
              z = (z_injection * alpha + (1.0_pr - alpha) * z_0)
              dzda = z_injection - z_0
          case("dilute")
              z = (z_injection * alpha + z_0)/sum(z_injection * alpha + z_0)
              dzda = -(alpha*z_injection + z_0) &
                     * sum(z_injection) / sum(alpha*z_injection + z_0)**2 &
                     + z_injection / sum(alpha*z_injection + z_0)
          case default
              z = (z_injection * alpha + (1.0_pr - alpha) * z_0)
              dzda = z_injection - z_0
      end select
      if (any(z < 0)) z = 0
      
      y = K * z

      call TERMO(n, 0, 4, T, P, y, Vy, lnfug_y, dlnphi_dp_y, dlnphi_dt_y, dlnphi_dn_y)
      call TERMO(n, 0, 4, T, P, z, Vz, lnfug_z, dlnphi_dp_z, dlnphi_dt_z, dlnphi_dn_z)

      F(1:n) = X(:n) + lnfug_y - lnfug_z
      F(n+1) = sum(y - z)
      F(n+2) = X(ns) - S

      df = 0

      do i=1,n
         do j=1,n
            df(i, j) = y(j) * dlnphi_dn_y(i, j)
         end do
         df(i, i) = df(i, i) + 1
         df(i, n+2) = sum(K * dlnphi_dn_y(i, :) * dzda - dlnphi_dn_z(i, :) * dzda)
      end do

      df(:n, n+1) = P * (dlnphi_dp_y - dlnphi_dp_z)
      df(n+1, :n) = y
      df(n+1, n+2) = sum(dzda*(K-1))

      df(n+2, :)  = 0
      df(n+2, ns) = 1
   end subroutine
   
   subroutine F_injection_three_phases(Xvars, ns, S, F, dF)
      use legacy_ar_models, only: TERMO
      !! Function to solve at each point of a three phase envelope.
      !!
      !! The vector of variables X corresponds to:
      !! \( X = [lnKx_i, lnKy_i lnP, \alpha, \beta] \)
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
      use iso_fortran_env, only: error_unit
      real(pr), intent(in)  :: Xvars(:) !! Vector of variables
      integer,  intent(in)  :: ns   !! Number of specification
      real(pr), intent(in)  :: S    !! Specification value
      real(pr), intent(out) :: F(size(Xvars)) !! Vector of functions valuated
      real(pr), intent(out) :: df(size(Xvars), size(Xvars)) !! Jacobian matrix

#define N (size(Xvars) - 3 )/2
      ! Xvars variables
      real(pr) :: z(N)
      real(pr) :: Kx(N)
      real(pr) :: Ky(N)
      real(pr) :: P
      real(pr) :: beta
      real(pr) :: alpha

      ! Main phase 1 variables
      real(pr) :: Vx
      real(pr), dimension(N) :: x, lnfug_x, dlnphi_dt_x, dlnphi_dp_x
      real(pr), dimension(N, N) :: dlnphi_dn_x
 
      ! Main phase 2 variables
      real(pr) :: Vy
      real(pr), dimension(N) :: y, lnfug_y, dlnphi_dt_y, dlnphi_dp_y
      real(pr), dimension(N, N) :: dlnphi_dn_y

      ! Incipient phase variables
      real(pr) :: Vw
      real(pr), dimension(N) :: w, lnfug_w, dlnphi_dt_w, dlnphi_dp_w
      real(pr), dimension(N, N) :: dlnphi_dn_w
 
      ! Derivative of z wrt alpha
      real(pr) :: dzda(N), dwda(N)

      ! Derivative of w wrt beta
      real(pr) :: dwdb(N)

      real(pr) :: dwdKx(N), dxdKx(N), dydKx(N)
      real(pr) :: dwdKy(N), dxdKy(N), dydKy(N)

      integer :: i, j, n
 
      n = N
#undef N

      Kx    = exp(Xvars(1:n))
      Ky    = exp(Xvars(n+1:2*n))
      P     = exp(Xvars(2*n+1))
      alpha = Xvars(2*n+2)
      beta  = Xvars(2*n+3)

      call get_z(alpha, z, dzda)
      if (any(z < 0)) z = 0

      w = z / (beta * Ky + (1-beta) * Kx)
      x = w * Kx
      y = w * Ky

      call TERMO(n, 0, 4, T, P, x, Vx, lnfug_x, dlnphi_dp_x, dlnphi_dt_x, dlnphi_dn_x)
      call TERMO(n, 0, 4, T, P, y, Vy, lnfug_y, dlnphi_dp_y, dlnphi_dt_y, dlnphi_dn_y)
      call TERMO(n, 0, 4, T, P, w, Vw, lnfug_w, dlnphi_dp_w, dlnphi_dt_w, dlnphi_dn_w)

      F(1:n)     = Xvars(1:n)     + lnfug_x - lnfug_w
      F(n+1:2*n) = Xvars(n+1:2*n) + lnfug_y - lnfug_w
 
      F(2*n+1) = sum(w) - 1
      F(2*n+2) = sum(x - y)
      F(2*n+3) = Xvars(ns) - S
 
      df = 0
      dwda = 1.0_pr / (beta * Ky + (1-beta) * Kx) * dzda
      dwdb = z * (Kx - Ky) / ((1 - beta) * Kx + beta * Ky)**2

      dwdKx = -z * (1-beta) / (Ky*beta + (1-beta)*Kx)**2
      dxdKx = Kx * dwdKx + w
      dydKx = Ky * dwdKx

      dwdKy = -z * (beta) / (Ky*beta + (1-beta)*Kx)**2
      dxdKy = Kx * dwdKy
      dydKy = Ky * dwdKy + w

      do i=1,n
         do j=1,n
            df(i, j)   = Kx(j) * (dlnphi_dn_x(i, j)  * dxdKx(j) - dlnphi_dn_w(i,j) * dwdKx(j))
            df(i+n, j) = Kx(j) * (dlnphi_dn_y(i, j)  * dydKx(j) - dlnphi_dn_w(i,j) * dwdKx(j))
            
            df(i, j+n)   = Ky(j) * (dlnphi_dn_x(i, j)  * dxdKy(j) - dlnphi_dn_w(i,j) * dwdKy(j))
            df(i+n, j+n) = Ky(j) * (dlnphi_dn_y(i, j)  * dydKy(j) - dlnphi_dn_w(i,j) * dwdKy(j))
         end do

         df(i, i)     = df(i, i)     + 1
         df(i+n, i+n) = df(i+n, i+n) + 1

         df(i,   2*n+2) = sum(Kx * dlnphi_dn_x(i, :) * dwda - dlnphi_dn_w(i, :) * dwda)
         df(i+n, 2*n+2) = sum(Ky * dlnphi_dn_y(i, :) * dwda - dlnphi_dn_w(i, :) * dwda)
         
         df(i,   2*n+3) = sum(Kx * dlnphi_dn_x(i, :) * dwdb - dlnphi_dn_w(i, :) * dwdb)
         df(i+n, 2*n+3) = sum(Ky * dlnphi_dn_y(i, :) * dwdb - dlnphi_dn_w(i, :) * dwdb)

         df(2*n+1, i)   = Kx(i) * dwdKx(i)
         df(2*n+1, i+n) = Ky(i) * dwdKy(i)

         df(2*n+2, i)   = Kx(i) * dxdKx(i) - Kx(i) * dydKx(i)
         df(2*n+2, i+n) = Ky(i) * dxdKy(i) - Ky(i) * dydKy(i)
      end do

      ! Derivatives wrt P
      df(:n,      2*n+1) = P * (dlnphi_dp_x - dlnphi_dp_w)
      df(n+1:2*n, 2*n+1) = P * (dlnphi_dp_y - dlnphi_dp_w)

      ! Derivatives wrt alpha
      df(2*n+1, 2*n+2) = sum(dwda)
      df(2*n+2, 2*n+2) = sum(Kx * dwda - Ky * dwda)
     
      ! Derivatives wrt beta
      df(2*n+1, 2*n+3) = sum(dwdb)
      df(2*n+2, 2*n+3) = sum(Kx * dwdb - Ky * dwdb)

      ! Derivatives wrt Xs
      df(2*n+3,  :) = 0
      df(2*n+3, ns) = 1
   end subroutine

   subroutine injection_envelope(X0, spec_number, del_S0, envels)
       use constants, only: ouput_path
       !! Subroutine to calculate Px phase envelopes via continuation method
       real(pr), intent(in) :: X0(:) !! Vector of variables
       integer,  intent(in) :: spec_number !! Number of specification
       real(pr), intent(in) :: del_S0 !! \(\Delta S_0\)
       type(injelope), intent(out) :: envels !! Calculated envelopes

       type(critical_point), allocatable :: cps(:)

       real(pr) :: X(size(X0))
       integer :: ns
       real(pr) :: S
       real(pr) :: XS(max_points, size(X0))

       real(pr) :: F(size(X0)), dF(size(X0), size(X0)), dXdS(size(X0))

       integer :: point, iters, n
       integer :: i
       integer :: funit_output
       character(len=254) :: fname_env

       allocate(cps(0))
       X = X0
       n = size(X0) - 2
       ns = spec_number
       S = X(ns)
       del_S = del_S0

       ! ======================================================================
       !  Output file
       ! ----------------------------------------------------------------------
       env_number = env_number + 1
       
       write(fname_env, *) env_number
       fname_env = "env-2ph-PX" // "_" // trim(adjustl(fname_env))
       fname_env = trim(adjustl(ouput_path)) // trim(fname_env) // ".dat"
       
       open(funit_output, file=fname_env)
       write(funit_output, * ) "#", T
       write(funit_output, *) "X0", iters, ns, X(n+2), exp(X(n+1)),  X(:n)
       ! ======================================================================

       enveloop: do point=1, max_points
          call full_newton(f_injection, iters, X, ns, S, F, dF)

          if (iters >= max_iters) then
             exit enveloop
          end if

          write(funit_output, *) "SOL", iters, ns, X(n+2), exp(X(n+1)),  X(:n)
          XS(point, :) = X

          update_spec: block
             real(pr) :: dFdS(size(X0))
             integer  :: ns_new

             dFdS = 0
             dFdS(n+2) = 1

             dXdS = solve_system(dF, dFdS)

             ns_new = maxloc(abs(dXdS), dim=1)

             if (ns_new /= ns) then
                del_S = dXdS(ns_new) * del_S  ! translation of delS to the  new specification variable
                dXdS = dXdS/dXdS(ns_new)
                ns = ns_new
             end if

             del_S = sign(1.0_pr, del_S) * minval( [ &
                max(sqrt(abs(X(ns)))/10, 0.1),       &
                abs(del_S) * 3/iters                 &
               ]                                     &
             )

             if (injection_case == "dilution") del_S = 50*del_S
          end block update_spec

          fix_step: block
             real(pr) :: Xnew(size(X0))
             real(pr) :: dP, dalpha

             Xnew = X + dXdS * del_S
             dP = exp(Xnew(n+1)) - exp(X(n+1))
             dalpha = Xnew(n+2) - X(n+2)

             do while (abs(dP) > 50 .or. abs(dalpha) > 0.03)
                  dXdS = dXdS/10.0_pr

                  Xnew = X + dXdS * del_S
                  dP = exp(Xnew(n+1)) - exp(X(n+1))
                  dalpha = Xnew(n+2) - X(n+2)
               end do
          end block fix_step

          detect_critical: block
             real(pr) :: K(size(X0) - 2), Knew(size(X0) - 2), Xnew(size(X0)), fact
             real(pr) :: pc, alpha_c, dS_c
             integer :: max_changing
             fact = 2.5

             Xnew = X + fact * dXdS * del_S

             K = X(:n)
             Knew = Xnew(:n)

             if (all(K * Knew < 0)) then
                max_changing = maxloc(abs(K - Knew), dim=1)

                dS_c = - k(max_changing) * (Xnew(ns) - X(ns))/(Knew(max_changing) - K(max_changing))
                del_S = dS_c * 1.1

                Xnew = X + dXdS * dS_c
                alpha_c = Xnew(n+2)
                pc = Xnew(n+1)

                cps = [cps, critical_point(t, pc, alpha_c)]
                write(funit_output, *) ""
                write(funit_output, *) ""
             end if
          end block detect_critical

          X = X + dXdS * del_S
          S = X(ns)

          if (any(break_conditions(X, ns, S))) exit enveloop
       end do enveloop

       point = point - 1

       write(funit_output, *) "#critical"
       if (size(cps) > 0 ) then
          do i=1,size(cps)
             write(funit_output, *) cps(i)%t, cps(i)%p
          end do
       else
          write(funit_output, *) "NaN NaN"
       endif
       
       close(funit_output)
       envels%z = z_0
       envels%z_inj = z_injection
       envels%logk = XS(:point, :n)
       envels%alpha = XS(:point, n+2)
       envels%p = exp(XS(:point, n+1))
       envels%critical_points = cps
   end subroutine
   
   subroutine injection_envelope_three_phase(X0, spec_number, del_S0, envels)
       use constants, only: ouput_path
       !! Subroutine to calculate Px phase envelopes via continuation method.
       !! Three phases version.
       real(pr), intent(in) :: X0(:) !! Vector of variables
       integer,  intent(in) :: spec_number !! Number of specification
       real(pr), intent(in) :: del_S0 !! \(\Delta S_0\)
       type(injelope), intent(out) :: envels !! Calculated envelopes

       type(critical_point), allocatable :: cps(:)

       real(pr) :: X(size(X0))
       integer :: ns
       real(pr) :: S
       real(pr) :: XS(max_points, size(X0))

       real(pr) :: F(size(X0)), dF(size(X0), size(X0)), dXdS(size(X0))

       integer :: point, iters, n
       integer :: i
       integer :: funit_output
       character(len=254) :: fname_env

       allocate(cps(0))
       X = X0
       n = (size(X0) - 3)/2
       ns = spec_number
       S = X(ns)
       del_S = del_S0

       ! ======================================================================
       !  Output file
       ! ----------------------------------------------------------------------
       env_number = env_number + 1

       write(fname_env, *) env_number
       fname_env = "env-3ph-PX" // "_" // trim(adjustl(fname_env))
       fname_env = trim(adjustl(ouput_path)) // trim(fname_env) // ".dat"

       open(newunit=funit_output, file=fname_env)
       write(funit_output, * ) "#", T
       write(funit_output, *) "X0", iters, ns, X(2*n+2), exp(X(2*n+1)), X(2*n+3),  X(:2*n)
       ! ======================================================================

       enveloop: do point=1, max_points
          call full_newton(F_injection_three_phases, iters, X, ns, S, F, dF)
          if (iters >= max_iters) then
             exit enveloop
          end if

          write(funit_output, *) "SOL", iters, ns, X(2*n+2), exp(X(2*n+1)), X(2*n+3),  X(:2*n)
          XS(point, :) = X

          update_spec: block
             real(pr) :: dFdS(size(X0))
             integer  :: ns_new

             dFdS = 0
             ! Actually it's -dFdS
             dFdS(2*n+3) = 1

             dXdS = solve_system(dF, dFdS)


             if (maxval(abs(X(:2*n))) < 1) then
                ns_new = maxloc(abs(dXdS(:2*n)), dim=1)  ! T and P not allowed to be chosen close to a critical point
             else
                ns_new = maxloc(abs(dXdS), dim=1)
             end if

             if (ns_new /= ns) then
                del_S = dXdS(ns_new) * del_S  ! translation of delS to the  new specification variable
                dXdS = dXdS/dXdS(ns_new)
                ns = ns_new
             end if

             del_S = sign(1.0_pr, del_S) * minval( [ &
                max(sqrt(abs(X(ns))), 0.1),          &
                abs(del_S) * 3/iters                 &
               ]                                     &
             )

             if (injection_case == "dilution") del_S = 50*del_S
          end block update_spec

          fix_step: block
             real(pr) :: Xnew(size(X0))
             real(pr) :: dP, dalpha

             Xnew = X + dXdS * del_S
             dP = exp(Xnew(2*n+1)) - exp(X(n+1))
             dalpha = Xnew(2*n+2) - X(n+2)

             do while (abs(dP) > 50 .or. abs(dalpha) > 0.03)
               dXdS = dXdS/10.0_pr

               Xnew = X + dXdS * del_S
               dP = exp(Xnew(2*n+1)) - exp(X(2*n+1))
               dalpha = Xnew(2*n+2) - X(2*n+2)
             end do
          end block fix_step

          detect_critical: block
             real(pr) :: K((size(X0) - 3)/2), Knew((size(X0) - 3)/2), Xnew(size(X0)), fact
             real(pr) :: pc, alpha_c, dS_c, dXdS_in(size(X0))
             integer :: max_changing, i
             fact = 15.0_pr

             Xnew = X + fact * dXdS * del_S
             do i=0,1
                K = X(i*n+1:(i+1)*n)
                Knew = Xnew(i*n+1:(i+1)*n)
                max_changing = minloc(abs(Knew - K), dim=1)

                if (all(K * Knew < 0)) then
                   dS_c = - k(max_changing) * (Xnew(ns) - X(ns))/(Knew(max_changing) - K(max_changing))
                   del_S = sign(15.0_pr, dS_c) ! dS_c * 15_pr

                   Xnew = X + dXdS * dS_c
                   alpha_c = Xnew(2*n+2)
                   pc = exp(Xnew(2*n+1))

                   cps = [cps, critical_point(t, pc, alpha_c)]
                   write(funit_output, *) ""
                   write(funit_output, *) ""
                end if
             end do
          end block detect_critical

          if (x(2*n+3) > 1) exit enveloop

          X = X + dXdS * del_S
          S = X(ns)
          if (any(break_conditions_three_phases(X, ns, S))) exit enveloop
       end do enveloop

       point = point - 1

       write(funit_output, *) ""
       write(funit_output, *) ""
       write(funit_output, *) "#critical"
       if (size(cps) > 0) then
          do i=1,size(cps)
             write(funit_output, *) cps(i)%alpha, cps(i)%p
          end do
       else
          write(funit_output, *) "NaN NaN"
       endif

       close(funit_output)
       envels%z = z_0
       envels%z_inj = z_injection
       envels%logk = XS(:point, :n)
       envels%alpha = XS(:point, n+2)
       envels%p = exp(XS(:point, n+1))
       envels%critical_points = cps
   end subroutine

   subroutine full_newton(fun, iters, X, ns, S, F, dF)
      !! Subroutine to solve a point in the envelope.
      interface
         subroutine fun(X, ns, S, F, dF)
            !! Function to solve
            import pr
            real(pr), intent(in) :: X(:)
            integer,  intent(in) :: ns
            real(pr), intent(in) :: S
            real(pr), intent(out) :: F(size(X))
            real(pr), intent(out) :: dF(size(X), size(X))
         end subroutine
      end interface
       integer,  intent(out)    :: iters !! Number of iterations needed
       real(pr), intent(in out) :: X(:)  !! Variables vector
       integer,  intent(in)     :: ns    !! Number of specification
       real(pr), intent(in)     :: S     !! Specification value
       real(pr), intent(out)    :: F(size(X)) !! Function values at solved point
       real(pr), intent(out)    :: df(size(X), size(X)) !! Jacobian values

       real(pr) :: b(size(X)), A(size(X), size(X))

       real(pr) :: dX(size(X)), tol=1e-5

       dX = 20

       newton: do iters=1, max_iters*10
           if (maxval(abs(dx)) < tol) exit newton
           call fun(X, ns, S, b, A)

           b = -b
           dX = solve_system(A, b)

           do while (maxval(abs(dX)) > 1)
              dX = dX/10
           end do

           X = X + dX
       end do newton

       F = b
       dF = A
   end subroutine

   function break_conditions(X, ns, S)
      !! Set of conditions to break the tracing.
      real(pr) :: X(:)
      integer :: ns
      real(pr) :: S

      integer :: n
      real(pr) :: p, alpha
      logical, allocatable :: break_conditions(:)

      n = size(X) - 2
      p = exp(X(n+1))
      alpha = X(n+2)

      break_conditions = [&
          p < 10 .or. p > 1000, &
          abs(del_S) < 1e-8     &
      ]
   end function
end module
