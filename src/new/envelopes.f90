module envelopes
   !! Functions to be used in the different continuation methods to trace
   !! phase envelopes 
   use constants, only: pr
   use linalg, only: solve_system, full_newton
   use dtypes, only: AbsEnvel, envelope, critical_point
   use legacy_ar_models, only: nc
   use legacy_thermo_properties, only: termo
   ! use progress_bar_module, only: progress_bar
   implicit none

   type, extends(AbsEnvel) :: PTEnvel3
      integer :: n
      real(pr), allocatable :: lnKx(:, :)
      real(pr), allocatable :: lnKy(:, :)
      real(pr), allocatable :: T(:)
      real(pr), allocatable :: P(:)
      real(pr), allocatable :: beta(:)
      type(critical_point), allocatable :: critical_points(:)
   end type

   integer, parameter :: max_points = 2000
   integer, parameter :: max_iters = 100
   integer :: env_number = 0

   interface
      function F(X, ns, S)
         import pr
         real(pr), intent(in) :: X(:)
         integer, intent(in) :: ns
         real(pr), intent(in) :: S
         real(pr) :: F
      end function
   end interface
contains
   ! ===========================================================================
   !  Initializators
   ! ---------------------------------------------------------------------------
   subroutine k_wilson_bubble(z, t_0, p_end, t, p, k)
      !! Find the Wilson Kfactors at ~10 bar to initialize a bubble point
      use legacy_ar_models, only: pc, tc, w
      real(pr), intent(in) :: z(:)
      real(pr), intent(in) :: t_0
      real(pr), intent(in) :: p_end
      real(pr), intent(in out) :: p
      real(pr), intent(in out) :: t

      real(pr), intent(out) :: k(size(z))

      P = 100.0
      T = t_0

      do while (P > p_end)
         T = T - 5._pr
         P = 1.0_pr/sum(z * pc*exp(5.373_pr*(1 + w)*(1 - tc/T)))
      end do
      k = k_wilson(t, p)
   end subroutine

   function k_wilson(t, p) result(k)
      ! use system, only: pc, tc, w
      use legacy_ar_models, only: pc, tc, w
      real(pr), intent(in) :: t, p
      real(pr) :: k(size(pc))
      k = pc * exp(5.373_pr * (1.0_pr + w) * (1.0_pr - tc/t))/p
   end function

   function p_wilson(z, t) result(p)
      ! use system, only: pc, tc, w
      use legacy_ar_models, only: pc, tc, w
      real(pr), intent(in) :: t, z(:)
      real(pr) :: p
      P = 1.0_pr/sum(z*pc*exp(5.373_pr*(1 + w)*(1 - tc/T)))
   end function
   ! ===========================================================================
   
   ! ===========================================================================
   ! General routines
   ! ---------------------------------------------------------------------------
   subroutine update_specification(iter, passingcri, X, dF, ns, S, delS, dXdS)
      integer,  intent(in)     :: iter
      logical,  intent(in)     :: passingcri
      real(pr), intent(in)     :: X(nc + 2)
      real(pr), intent(in)     :: dF(nc + 2, nc + 2)
      integer,  intent(in out) :: ns
      real(pr), intent(in out) :: S
      real(pr), intent(in out) :: delS
      real(pr), intent(in out) :: dXdS(nc + 2)

      real(pr) :: dF_dS(nc + 2)
      real(pr) :: bd(nc + 2)
      real(pr) :: AJ(nc + 2, nc + 2)
      real(pr) :: delmax, updel
      integer  :: nsold

      dF_dS = 0
      call dFdS(dF_dS)

      bd = -dF_dS
      AJ = dF
      dXdS = solve_system(AJ, bd)

      ! Selection of (the most changing) variable to be specified for the next point
      nsold = ns

      ns = maxloc(abs(dXdS), dim=1)

      if (maxval(abs(X(:nc))) < 0.2) then
         ns = maxloc(abs(dXdS(:nc)), dim=1)  ! T and P not allowed to be chosen close to a critical point
      end if

      if (ns /= nsold) then
         delS = dXdS(ns)*delS  ! translation of delS to the  new specification variable
         dXdS = dXdS/dXdS(ns)  ! translation of sensitivities
         S = X(ns)             ! update of S
      end if

      ! Setting step in S for the next point to be calculated
      delmax = max(sqrt(abs(X(ns)))/10, 0.1)
      updel = delS*3/iter

      if (passingcri) updel = delS
      if (delS > 0) then
         delS = min(updel, delmax)
      else
         delS = max(updel, -delmax)
      end if
      delS = 3*delS

      S = S + delS
   end subroutine
   ! ===========================================================================

   ! ===========================================================================
   ! Specification function derivatives
   ! ---------------------------------------------------------------------------
   subroutine dFdS(dF_dS)
      ! use system, only: nc
      use legacy_ar_models, only: nc
      real(pr), intent(out) :: dF_dS(nc + 2)

      dF_dS = 0
      dF_dS(nc + 2) = -1
   end subroutine
   ! ---------------------------------------------------------------------------
   
   ! ===========================================================================
   ! Two Phase envelopes
   ! ---------------------------------------------------------------------------
   function X2(kfact, P, T) result(X)
      real(pr), intent(in) :: kfact(nc)
      real(pr), intent(in) :: P
      real(pr), intent(in) :: T

      real(pr) :: X(nc + 2)

      integer :: n

      n = size(kfact)

      X(:n) = log(kfact)
      X(n + 1) = log(T)
      X(n + 2) = log(P)
   end function

   subroutine F2(incipient, z, y, X, S, ns, F, dF)
      character(len=*), intent(in) :: incipient
      real(pr), intent(in) :: z(:)
      real(pr), intent(in) :: X(nc + 2)
      real(pr), intent(in) :: y(nc)
      real(pr), intent(in) :: S
      integer, intent(in) :: ns

      real(pr), intent(out) :: F(nc + 2)
      real(pr), intent(out) :: dF(nc + 2, nc + 2)

      real(pr) :: Vx, Vy, lnfug_x(nc), lnfug_y(nc)
      real(pr) :: dlnphi_dt_x(nc), dlnphi_dt_y(nc)
      real(pr) :: dlnphi_dp_x(nc), dlnphi_dp_y(nc)
      real(pr) :: dlnphi_dn_x(nc, nc), dlnphi_dn_y(nc, nc)

      real(pr) :: T, P

      integer :: ix, iy, n, j

      n = size(z)
      F = 0
      dF = 0

      T = exp(X(n+1))
      P = exp(X(n+2))

      select case(incipient)
      case ("liquid")
         ix = -1
         iy = 1
      case ("vapor")
         ix = 1
         iy = -1
      case ("2ndliquid")
         ix = 1
         iy = 1
      case default
         ix = 0
         iy = 0
      end select

      call TERMO(n, iy, 4, T, P, y, Vy, lnfug_y, dlnphi_dp_y, dlnphi_dt_y, dlnphi_dn_y)
      call TERMO(n, ix, 2, T, P, z, Vx, lnfug_x, dlnphi_dp_x, dlnphi_dt_x, dlnphi_dn_x)

      F(:n) = X(:n) + lnfug_y - lnfug_x  ! X(:n) are LOG_K
      F(n + 1) = sum(y - z)
      F(n + 2) = X(ns) - S

      ! Jacobian Matrix
      do j=1,n
         df(:n, j) = dlnphi_dn_y(:, j) * y(j)
         df(j, j) = dF(j, j) + 1
      end do

      df(:n, n + 1) = T * (dlnphi_dt_y - dlnphi_dt_x)
      df(:n, n + 2) = P * (dlnphi_dp_y - dlnphi_dp_x)

      df(n + 1, :n) = y

      df(n + 2, :) = 0
      df(n + 2, ns) = 1
   end subroutine F2

   subroutine fix_delx(&
         point, iterations, desired_iterations, first_tol, tol, delX &
      )
      integer, intent(in)  :: point
      integer, intent(in)  :: iterations
      integer, intent(in)  :: desired_iterations
      
      real(pr), intent(in) :: first_tol
      real(pr), intent(in) :: tol
      real(pr), intent(in out) :: delX(:)

      if (point == 1) then
         do while (maxval(abs(delX)) > first_tol) 
            ! Too large Newton step --> Reduce it
            delX = delX/2
         end do
      else
         do while (maxval(abs(delX)) > tol)   
            ! Too large Newton step --> Reduce it
            delX = delX/2
         end do
         if (iterations > desired_iterations)  then
            ! too many iterations (sometimes due to oscillatory behavior 
            ! near critical point) --> Reduce it
            delX = delX/2
         endif
      end if
   end subroutine

   subroutine find_hpl(t, p, k)
      !! Find a HPLL initial point at a given pressure
      use legacy_ar_models, only: nc, z
      use legacy_thermo_properties, only: termo
      real(pr), intent(in out) :: t
      real(pr), intent(in) :: p
      real(pr), intent(out) :: k(nc)

      integer :: i, ncomp
      real(pr) :: diff
      real(pr) :: v
      real(pr) :: x(nc), y(nc), lnfug_z(nc), lnfug_y(nc)

      diff = -1

      ncomp = nc

      y = 0
      y(ncomp) = 1

      do while(diff < 0)
         t = t - 1.0_pr
         call termo(nc, 4, 1, t, p, z, v, philog=lnfug_z)
         call termo(nc, 4, 1, t, p, y, v, philog=lnfug_y)
         diff = (log(z(ncomp)) + lnfug_z(ncomp)) - (log(y(ncomp)) + lnfug_y(ncomp))
      end do

      k = 1/exp(lnfug_y - lnfug_z)
      k = 1e-3
      k(ncomp) = 1000
      ! print *, k
      ! print *, 1/k
   end subroutine

   subroutine envelope2(ichoice, n, z, T, P, KFACT, & ! This will probably always exist
                        n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, & ! This shouldnt be here in the future
                        this_envelope) ! This output should encapsulate everything
      use dtypes, only: envelope, critical_point
      use linalg, only: point, solve_system
      use constants, only: ouput_path
      use io, only: str
      implicit none

      ! number of compounds in the system and starting point type
      integer, intent(in) :: n, ichoice

      ! estimated T and P for first point (then used for every point)
      real(pr) :: T, P

      ! Maximun pressure
      real(pr) :: maxP

      ! estimated K factors for first point (then used for every point)
      real(pr), intent(in out) :: KFACT(n)

      ! composition of the system
      real(pr), intent(in) :: z(n)

      ! T, P and Density of the calculated envelope
      real(pr), intent(out) :: Tv(max_points)
      real(pr), intent(out) :: Pv(max_points)
      real(pr), intent(out) :: Dv(max_points)

      ! number of valid elements in Tv, Pv and Dv arrays
      integer, intent(out) :: n_points

      ! positions of the last saturation points before each critical point
      integer, dimension(4), intent(out) :: icri

      ! T, P and Density of critical points
      real(pr), dimension(4), intent(out) :: Tcri(4), Pcri(4), Dcri(4)

      ! number of valid elements in icri, Tcri, Pcri and Dcri arrays
      integer, intent(out) :: ncri

      ! Intermediate variables during calculation process
      real(pr), dimension(n) :: y
      real(pr), dimension(n + 2) :: X, Xold, Xold2, delX, bd, F, dFdS, dXdS
      real(pr), dimension(n + 2, n + 2) :: JAC, AJ
      real(pr) :: Vx
      logical :: run, passingcri, minT, minmaxT

      character(len=:), allocatable :: incipient_phase

      type(envelope), intent(out) :: this_envelope
      real(pr) :: tmp_logk(max_points, n)
      real(pr) :: tmp_logphi(max_points, n)

      ! Extrapolation of variables to detect critical points
      real(pr) :: extra_slope(n + 2)
      real(pr) :: lnK_extrapolated(n)
      real(pr) :: delta_t

      integer :: i
      integer :: iy, ix ! Vapor or liquid selectors

      ! Specification value, delta and index
      real(pr) :: S, delS
      integer :: ns

      real(pr) :: Told2, Told
      real(pr) :: frac

      ! Netwon method
      integer :: iter ! Iteration
      integer :: max_iter

      ! Critical Points
      type(critical_point), allocatable :: critical_points(:)
      integer :: black_i ! Number of steps while trying to escape the CP
      real(pr) :: stepx

      integer :: funit_output
      character(len=254) :: fname_env

      ! ========================================================================
      !  OUTPUT file
      ! ------------------------------------------------------------------------
      env_number = env_number + 1
      write(fname_env, *) env_number
      fname_env = "env-2ph-PT" // "_" // trim(adjustl(fname_env))
      fname_env = trim(adjustl(ouput_path)) // trim(fname_env) // ".dat"

      open(newunit=funit_output, file=fname_env)
      ! ========================================================================

      ! Initialize with zero Tv and Pv
      allocate(this_envelope%vars(max_points-50, n+2))
      Tv = 0
      Pv = 0

      minT = .false.
      minmaxT = .false.
      passingcri = .false.
      Told2 = 0.0
      Told = 10.0
      maxP = 0.d0

      !-----------------------------------------------------------
      ! Continuation method for tracing the envelope starts here
      run = .true.
      i = 0
      ncri = 0
      JAC(n + 1, :) = 0.d0
      X(:n) = log(KFACT)
      X(n + 1) = log(T)
      X(n + 2) = log(P)
      iy = 1
      ix = 1

      select case(ichoice)
      case (1)
         incipient_phase = "vapor"
         iy = -1
      case (2)
         incipient_phase = "liquid"
         ix = -1
      case (3)
         incipient_phase = "2ndliquid"
      end select
      write(funit_output, *) "#", incipient_phase
      write(funit_output, "(*(A,2x))") "STAT", "iter", "ns", "T", "P", &
         ("lnK"//str(i),i=1,n), ("y"//str(i), i=1,n)

      if (ichoice <= 2) then
         ! low T bub (1) or dew (2)
         ! x will be vapor phase during the first part, 
         ! and liquid after a critical point is crossed
         if (ichoice == 1) iy = -1
         if (ichoice == 2) ix = -1
         ns = n + 1
         S = log(T)
         delS = 0.005

         ! Wilson estimate for vapor (or liquid) composition
         y = KFACT*z
      else
         ! (ichoice==3) high P L-L sat
         ! PmaxDewC = maxval(PdewC(1:ilastDewC))
         ns = n + 2
         S = log(P)
         delS = -0.05
         y = kfact * z
         ! y = 0.d0
         ! y(n) = 1.d0
      end if

      Xold = 0.d0
      dFdS = 0.d0
      dFdS(n + 2) = -1.d0

      i=0

      do while (run)
         i = i + 1
         if (i > max_points - 50) then
            exit
         end if
         ! Newton starts here
         delX = 1.0
         iter = 0
         max_iter = 500

         do while (maxval(abs(delX)) > 1.d-9 .and. iter <= max_iter)
            ! Solve point with full Newton method
            call F2(incipient_phase, z, y, X, S, ns, F, JAC)

            iter = iter + 1

            bd = -F
            AJ = JAC
            delX = solve_system(AJ, bd)
            call fix_delX(i, iter, 3, 10.0_pr, 0.08_pr, delX)

            X = X + delX

            if (.not. passingcri .and. i /= 1 &
                .and. iter > 10 &
                .and. maxval(abs(delX)) > 0.001) then 
               ! Too many iterations --> Reduce step to new point

               delS = delS*2.0/4.0
               S = S - delS
               X = Xold + dXdS*delS
            end if

            KFACT = exp(X(:n))
            y = z*KFACT
            T = exp(X(n + 1))
            P = exp(X(n + 2))
         end do

         ! Point converged (unless it jumped out because of high number of iterations)
         if (iter > max_iter) run = .false.
         if (P > maxP) maxP = P

         if (run) write(funit_output, *) "SOL", iter, ns, T, P, X(:n), z*exp(X(:n))

         if (incipient_phase == "liquid" .and. i > 1) then
            ! TODO: If this is the way the low p dew line finishes, 
            ! I think this could be better, like using dPdT
            if (P < Pv(i - 1) .and. P < maxP/5 .and. T > 300) then
               run = .true.  ! to finish envelope going to low T bubble
            end if
         end if

         Tv(i) = T
         Pv(i) = P
         this_envelope%vars(i, :) = X
         tmp_logk(i, :n) = X(:n)

         if (incipient_phase == "2ndliquid" .and. P < 0.1) then
            ! isolated LL line detected. 
            ! Stop and start a new one from low T false bubble point
            run = .false.
         end if

         if (i > max_points - 50) exit

         if (sum(X(:n) * Xold(:n)) < 0) then  ! critical point detected
            ncri = ncri + 1
            icri(ncri) = i - 1
            frac = -Xold(ns)/(X(ns) - Xold(ns))
            Tcri(ncri) = Tv(i - 1) + frac*(T - Tv(i - 1))
            Pcri(ncri) = Pv(i - 1) + frac*(P - Pv(i - 1))
            Dcri(ncri) = Dv(i - 1) + frac*(Dv(i) - Dv(i - 1))

            select case (incipient_phase)
            case("liquid")
               incipient_phase = "vapor"
            case("vapor")
               incipient_phase = "liquid"
            end select

            write(funit_output, *) " "
            write(funit_output, *) " "
            write(funit_output, *) "#", incipient_phase
         end if

         if (run) then
            ! Calculation of sensitivities (dXdS)
            ! dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )

            call update_specification(iter, passingcri, X, JAC, ns, S, delS, dXdS)

            ! Generation of estimates for the next point
            Told2 = Told
            Told = T
            Xold2 = Xold
            Xold = X
            X = Xold + dXdS * delS

            critical_region: block
            black_i = 0

            if (passingcri) passingcri = .false.

            if (i > 10) then
               ! After the 10th step extrapolate the K factors in order
               ! to find a possible critical point (check if all lnK values
               ! change sign).
               extrapolation: block
                  integer :: loc
                  integer :: its
                  real(pr) :: delta
                  real(pr) :: m(size(X))
                  real(pr) :: max_lnK, max_lnK2, delta_lnK
                  real(pr) :: delta_X(size(x))
                  real(pr) :: Xin(size(X))

                  its = 0
                  delta = delS

                  ! Variation of lnK based on deltaS
                  ! m = 37.0_pr * (X - Xold2)/(delta)
                  ! lnK_extrapolated = (delta) * m(:n) + X(:n)
                  lnK_extrapolated = X(:n) + 7 * delS * dXdS(:n)
                  ! print *, X(:n)
                  ! print *, lnK_extrapolated

                  if (all((X(:n) * lnK_extrapolated < 0), dim=1)) then
                     print *, "Extrapol CP"
                     ! All lnK changed sign, so a CP is inminent
                     ! aproach it enough to force the jumping algorithm
                     do while( &
                           maxval(abs(X(:n))) >= 0.02 &
                           .and. all(X(:n)*lnK_extrapolated > 0, dim=1)&
                        )
                        print *, its, "Getting to critical", &
                           exp(X(n+1)), exp(X(n+2)), maxval(abs(X(:n))), &
                           all(X(:n)*lnK_extrapolated > 0, dim=1)
                        its = its + 1
                        delta_x = delta * m
                        X = delta_x + X
                        S = X(ns)
                        if (its > 10) exit
                     end do
                     passingcri = .true.
                  end if
               end block extrapolation
            end if

            do while (maxval(abs(X(:n))) < 0.03)
               ! print *, "Jumping critical"
               ! approaching the black hole... get out of there! (0.03)
               black_i = black_i + 1
               if (black_i > 50) then
                   print *, "Stuck on the black hole"
                   if (black_i > 100) stop
               end if

               stepX = maxval(abs(X(:n) - Xold(:n))) 

               ! the step given by the most changing logK to fall into 
               ! the black hole
               S = S + delS
               X = X + dXdS*delS
               ! passingcri = .true.
               ! if (stepX > 0.07) then
               !    !  half step back
               !    S = S - delS/2
               !    X = X - dXdS*delS/2   
               ! else
               !    ! one more step to jump over the critical point
               !    S = S + delS
               !    X = X + dXdS*delS   
               ! end if
            end do
            end block critical_region

            T = exp(X(n + 1))

            do while (.not. passingcri .and. abs(T - Told) > 15)
               ! Delta T estimations > 7K are not allowed
               delS = delS/2
               S = S - delS
               X = Xold + dXdS*delS
               T = exp(X(n + 1))
            end do

            P = exp(X(n + 2))
            KFACT = exp(X(:n))
            y = z*KFACT

            ! Finish conditions
            if ((dXdS(n + 1)*delS < 0 .and. P < 0.05 .and. T < 50.0) &  ! dew line stops when P<0.1 bar or T<50K
                .or. (P > 1.0 .and. T < 50.0) &   ! bubble line stops when T<50K
                .or. (P > 5000) &
                .or. (abs(dels) < 1.d-10)) then
                run = .false.
            end if
         end if
      end do
      !-----------------------------------------------------------

      n_points = i

      write(funit_output, *) " "
      write(funit_output, *) " "
      write(funit_output, *) "#critical"
      if (ncri == 0) write(funit_output, *) "NaN NaN"
      do i=1, ncri
         write(funit_output, *) Tcri(i), Pcri(i)
      end do

      ! Define envelope values, omit the last point to avoid not really
      ! converged cases
      close(funit_output)
      this_envelope%logk = tmp_logk(:n_points - 1, :)
      this_envelope%logphi = tmp_logphi(:n_points - 1, :)
      this_envelope%t = Tv(:n_points - 1)
      this_envelope%p = Pv(:n_points - 1)
      this_envelope%z = z

      allocate(critical_points(ncri))
      critical_points%t = tcri(:ncri)
      critical_points%p = pcri(:ncri)
      this_envelope%critical_points = critical_points
   end subroutine envelope2
   ! ===========================================================================

   ! ===========================================================================
   ! Three-phase envelopes
   ! ---------------------------------------------------------------------------
   subroutine pt_F_three_phases(Xvars, ns, S, F, dF)
      !! Function to solve at each point of a three phase envelope.
      !!
      !! The vector of variables X corresponds to:
      !! \( X = [lnKx_i, lnKy_i lnP, lnT, \beta] \)
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
      use legacy_ar_models, only: z
      use legacy_thermo_properties, only: termo
      use iso_fortran_env, only: error_unit
      real(pr), intent(in)  :: Xvars(:) !! Vector of variables
      integer, intent(in)  :: ns   !! Number of specification
      real(pr), intent(in)  :: S    !! Specification value
      real(pr), intent(out) :: F(size(Xvars)) !! Vector of functions valuated
      real(pr), intent(out) :: df(size(Xvars), size(Xvars)) !! Jacobian matrix

      ! Xvars variables
      ! real(pr) :: z((Size(Xvars)-3)/2)
      real(pr) :: Kx((Size(Xvars)-3)/2)
      real(pr) :: Ky((Size(Xvars)-3)/2)
      real(pr) :: P
      real(pr) :: T
      real(pr) :: beta

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

      ! Derivative of w wrt beta
      real(pr) :: dwdb((Size(Xvars)-3)/2)

      real(pr) :: dwdKx((Size(Xvars)-3)/2), dxdKx((Size(Xvars)-3)/2), dydKx((Size(Xvars)-3)/2)
      real(pr) :: dwdKy((Size(Xvars)-3)/2), dxdKy((Size(Xvars)-3)/2), dydKy((Size(Xvars)-3)/2)

      integer :: i, j, n

      n = (Size(Xvars)-3)/2

      Kx = exp(Xvars(1:n))
      Ky = exp(Xvars(n + 1:2*n))
      P = exp(Xvars(2*n + 1))
      T = exp(Xvars(2*n + 2))
      beta = Xvars(2*n + 3)

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

      F = 0

      F(1:n) = Xvars(1:n) + lnfug_x - lnfug_w
      F(n + 1:2*n) = Xvars(n + 1:2*n) + lnfug_y - lnfug_w

      F(2*n + 1) = sum(w) - 1
      F(2*n + 2) = sum(x - y)
      F(2*n + 3) = Xvars(ns) - S

      df = 0
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

         ! dlnK_i/dlnK_i
         df(i, i) = df(i, i) + 1
         df(i + n, i + n) = df(i + n, i + n) + 1

         df(i, 2*n + 3) = sum(Kx*dlnphi_dn_x(i, :)*dwdb - dlnphi_dn_w(i, :)*dwdb)
         df(i + n, 2*n + 3) = sum(Ky*dlnphi_dn_y(i, :)*dwdb - dlnphi_dn_w(i, :)*dwdb)

         df(2*n + 1, i) = Kx(i)*dwdKx(i)
         df(2*n + 1, i + n) = Ky(i)*dwdKy(i)

         df(2*n + 2, i) = Kx(i)*dxdKx(i) - Kx(i)*dydKx(i)
         df(2*n + 2, i + n) = Ky(i)*dxdKy(i) - Ky(i)*dydKy(i)
      end do

      ! Derivatives wrt P
      df(:n, 2*n + 1) = P*(dlnphi_dp_x - dlnphi_dp_w)
      df(n + 1:2*n, 2*n + 1) = P*(dlnphi_dp_y - dlnphi_dp_w)

      ! Derivatives wrt T
      df(:n, 2*n + 2) = T*(dlnphi_dt_x - dlnphi_dt_w)
      df(n + 1:2*n, 2*n + 2) = T*(dlnphi_dt_y - dlnphi_dt_w)

      ! Derivatives wrt beta
      df(2*n + 1, 2*n + 3) = sum(dwdb)
      df(2*n + 2, 2*n + 3) = sum(Kx*dwdb - Ky*dwdb)

      ! Derivatives wrt Xs
      df(2*n + 3, :) = 0
      df(2*n + 3, ns) = 1
   end subroutine

   subroutine pt_envelope_three_phase(X0, spec_number, del_S0, envel)
      use constants, only: ouput_path
      use io, only: str
      !! Subroutine to calculate Px phase envelopes via continuation method.
      !! Three phases version.
      real(pr), intent(in) :: X0(:) !! Vector of variables
      integer, intent(in) :: spec_number !! Number of specification
      real(pr), intent(in) :: del_S0 !! \(\Delta S_0\)
      type(PTEnvel3), intent(out) :: envel !! Calculated envelopes

      type(critical_point), allocatable :: cps(:)

      real(pr) :: X(size(X0))
      integer :: ns
      real(pr) :: S
      real(pr) :: del_S
      real(pr) :: XS(max_points, size(X0))

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
      fname_env = "env-3ph-PT"//"_"//trim(adjustl(fname_env))
      fname_env = trim(adjustl(ouput_path))//trim(fname_env)//".dat"

      open (newunit=funit_output, file=fname_env)
      write (funit_output, *) "#"
      write (funit_output, *) "STAT", " iters", " ns", " T", " P", &
         " beta", (" lnKx"//str(i), i=1,n), (" lnKy"//str(i), i=1,n)
      write (funit_output, *) "X0", iters, ns, exp(X(2*n + 2)), exp(X(2*n + 1)), &
         X(2*n + 3), X(:2*n)
      ! ======================================================================

      enveloop: do point = 1, max_points
         ! call progress_bar(point, max_points, advance=.false.)
         call full_newton(pt_F_three_phases, iters, X, ns, S, max_iters, F, dF, solvetol=1e-7_pr)
         if (iters >= max_iters) then
            print *, "Breaking: Above max iterations"
            exit enveloop
         end if

         write (funit_output, *) "SOL", iters, ns, exp(X(2*n + 2)), &
            exp(X(2*n + 1)), X(2*n + 3), X(:2*n)
         XS(point, :) = X

         update_spec: block
            real(pr) :: dFdS(size(X0))
            integer  :: ns_new

            dFdS = 0
            ! Actually it's -dFdS
            dFdS(2*n + 3) = 1

            dXdS = solve_system(dF, dFdS)

            if (maxval(abs(X(:2*n))) < 1) then
               ! T, P and beta not allowed near a CP
               ns_new = maxloc(abs(dXdS(:2*n)), dim=1)
            else
               ns_new = maxloc(abs(dXdS), dim=1)
            end if

            if (ns_new /= ns) then
               ! translation of delS to the  new specification variable
               del_S = dXdS(ns_new)*del_S  
               dXdS = dXdS/dXdS(ns_new)
               ns = ns_new
            end if
         end block update_spec

         fix_step: block
            real(pr) :: Xnew(size(X0))
            real(pr) :: dP, dT

            del_S = sign(1.0_pr, del_S) * minval([ &
                                               max(abs(sqrt(X(ns))/10), 0.1_pr), &
                                               abs(del_S)*3/iters &
                                               ] &
                                               )

            ! Xnew = X + dXdS*del_S
            ! dP = exp(Xnew(2*n + 1)) - exp(X(2*n + 1))
            ! dT = exp(Xnew(2*n + 2)) - exp(X(2*n + 2))

            ! do while (abs(dP) > 50 .or. abs(dT) > 7)
            !    dXdS = dXdS/2.0_pr

            !    Xnew = X + dXdS*del_S
            !    dP = exp(Xnew(2*n + 1)) - exp(X(2*n + 1))
            !    dT = exp(Xnew(2*n + 2)) - exp(X(2*n + 2))
            ! end do
         end block fix_step

         detect_critical: block
            real(pr) :: K((size(X0) - 3)/2), Knew((size(X0) - 3)/2), &
                        Xnew(size(X0)), fact
            real(pr) :: pc, tc, dS_c, dXdS_in(size(X0))
            integer :: max_changing, i
            fact = 4.0_pr

            loop: do i = 0, 1
               Xnew = X + fact*dXdS*del_S

               K = X(i*n + 1:(i + 1)*n)
               Knew = Xnew(i*n + 1:(i + 1)*n)
              
               max_changing = maxloc(abs(Knew - K), dim=1)

               if (all(K*Knew < 0) .or. maxval(abs(K)) < 0.1) then
                  print *, "CP"
                  dS_c = ( &
                         -k(max_changing)*(Xnew(ns) - X(ns)) &
                         /(Knew(max_changing) - K(max_changing)) &
                         )

                  Xnew = X + dXdS*dS_c
                  Tc = exp(Xnew(2*n + 2))
                  pc = exp(Xnew(2*n + 1))
                  cps = [cps, critical_point(tc, pc, 0.0_pr)]

                  del_S = dS_c + 2.5_pr * dS_c

                  write (funit_output, *) ""
                  write (funit_output, *) ""
                  exit loop
               end if
            end do loop
         end block detect_critical

         X = X + dXdS*del_S
         S = X(ns)
         if (any(break_conditions_three_phases(X, ns, S))) then
            print *, "Breaking: ", break_conditions_three_phases(X, ns, S)
            exit enveloop
         end if
      end do enveloop

      write (funit_output, *) ""
      write (funit_output, *) ""
      write (funit_output, *) "#critical"
      if (size(cps) > 0) then
         do i = 1, size(cps)
            write (funit_output, *) cps(i)%t, cps(i)%p
         end do
      else
         write (funit_output, *) "NaN NaN"
      end if

      close (funit_output)
      ! call progress_bar(point, max_points, .true.)
      envel%lnKx = XS(:point, :n)
      envel%lnKy = XS(:point, n+1:2*n)
      envel%P = exp(XS(:point, 2*n + 1))
      envel%T = exp(XS(:point, 2*n + 2))
      envel%beta = XS(:point, 2*n + 3)
      envel%critical_points = cps
   end subroutine

   function break_conditions_three_phases(X, ns, S)
      !! Set of conditions to break the tracing.
      real(pr) :: X(:) !! Variables vector
      integer :: ns !! Number of specification
      real(pr) :: S !! Value of specification

      integer :: n
      real(pr) :: p, t, beta
      logical, allocatable :: break_conditions_three_phases(:)

      n = (size(X) - 3)/2
      p = exp(X(2*n + 1))
      beta = x(2*n+3)
      t = X(2*n + 2)
         
      break_conditions_three_phases = [ &
                                          beta < 0 &
                                          .or. 1 < beta &
                                      ]
   end function
   ! ===========================================================================

   ! ===========================================================================
   !  Intersections and crossings
   ! ---------------------------------------------------------------------------
   subroutine get_case(&
      dew, bub, hpl, intersections, self_intersections, this_case &
   )
      use linalg, only: intersection, point
      type(envelope), intent(in) :: dew
      type(envelope), intent(in) :: bub
      type(envelope), intent(in) :: hpl
      type(point), allocatable, intent(out) :: intersections(:)
      type(point), allocatable, intent(out) :: self_intersections(:)
      character(len=:), allocatable, intent(out) :: this_case

      type(point), allocatable :: inter_dew_bub(:)
      type(point), allocatable :: inter_hpl_bub(:)
      type(point), allocatable :: inter_hpl_dew(:)
      
      type(point), allocatable :: self_inter_dew(:)

      inter_dew_bub = intersection(dew%t, dew%p, bub%t, bub%p)
      inter_hpl_bub = intersection(hpl%t, hpl%p, bub%t, bub%p)
      inter_hpl_dew = intersection(hpl%t, hpl%p, dew%t, dew%p)

      self_inter_dew = intersection(dew%t, dew%p)

      print *, "======INTERSERCTIONS========="
      print *, "DEW_BUB", inter_dew_bub
      print *, "HPL_BUB", inter_hpl_bub
      print *, "HPL_DEW", inter_hpl_dew
      print *, "============================="

      if (size(inter_dew_bub) == 2) then
         this_case = "2_DEW_BUB"
         intersections = inter_dew_bub
      else if (size(inter_hpl_bub) == 1 .and. size(inter_dew_bub) == 1) then
         this_case = "2_HPL_BUB_DEW_BUB"
         intersections = [inter_hpl_bub, inter_dew_bub]
      else if (size(inter_hpl_bub) == 1 .and. size(inter_hpl_dew) == 1) then
         this_case = "2_HPL_BUB_HPL_DEW"
         intersections = [inter_hpl_bub, inter_hpl_dew]
      else if (size(inter_hpl_bub) == 2) then
         this_case = "2_HPL_BUB"
         intersections = [inter_hpl_bub(1), inter_hpl_bub(2)]
      else if (size(inter_hpl_bub) == 1) then
         this_case = "1_HPL_BUB"
         intersections = inter_hpl_bub
      else if (size(inter_hpl_dew) == 1) then
         this_case = "1_HPL_DEW"
         intersections = inter_hpl_dew
      else if (size(self_inter_dew) == 1) then
         this_case = "1_DEW"
         allocate(intersections(0))
         self_intersections = self_inter_dew
      else
         this_case = "0"
         allocate(intersections(0))
      end if
   end subroutine

   subroutine pt_three_phase_from_intersection(&
         pt_x, pt_y, intersections, &
         pt_x_3, pt_y_3 &
      )
      use legacy_ar_models, only: z
      use linalg, only: point, interpol
      type(envelope), intent(in) :: pt_x, pt_y
      type(point), intent(in) :: intersections(:)
      type(PTEnvel3), intent(out) :: pt_x_3(:), pt_y_3(:)
      
      real(pr), allocatable :: lnKx(:), lnKy(:)
      real(pr), allocatable :: X(:)
      real(pr) :: t, p, beta, del_S0
      real(pr), allocatable :: phase_y(:), phase_x(:)
      integer :: i, j, i_inter=1
      integer :: ns

      do i_inter=1,size(intersections)
         i = intersections(i_inter)%i
         j = intersections(i_inter)%j

         t = intersections(i_inter)%x
         p = intersections(i_inter)%y

         lnKx = interpol( pt_x%t(i), pt_x%t(i + 1), &
                          pt_x%logk(i, :), pt_x%logk(i + 1, :), &
                        t &
                        )

         lnKy = interpol( &
                  pt_y%t(j), pt_y%t(j + 1), &
                  pt_y%logk(j, :), pt_y%logk(j + 1, :), &
                  t &
               )

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
         print *, "Three Phase: Gas"
         lnKx = log(phase_x/phase_y)
         lnKy = log(z/phase_y)
         X = [lnKx, lnKy, log(p), log(t), beta]
         call pt_envelope_three_phase(X, ns, del_S0, pt_x_3(i_inter))
         ! ==================================================================

         ! ==================================================================
         !  Line with incipient phase liquid
         ! ------------------------------------------------------------------
         print *, "Three Phase: Liquid"
         lnKx = log(phase_y/phase_x)
         lnKy = log(z/phase_x)
         X = [lnKx, lnKy, log(p), log(t), beta]
         call pt_envelope_three_phase(X, ns, del_S0, pt_y_3(i_inter))
         ! ==================================================================
      end do
      ! ===========================================================================
   end subroutine
end module envelopes
