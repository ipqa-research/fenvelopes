module envelopes
   !! Functions to be used in the different continuation methods to trace
   !! phase envelopes 
   use constants, only: pr
   use linalg, only: solve_system
   ! use system, only: nc
   use dtypes, only: envelope
   use legacy_ar_models, only: nc, termo
   implicit none

   integer, parameter :: max_points = 2000
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
   subroutine k_wilson_bubble(z, t, p, k)
      !! Find the Wilson Kfactors at ~10 bar to initialize a bubble point
      ! use system, only: pc, tc, w
      use legacy_ar_models, only: pc, tc, w
      real(pr), intent(in) :: z(:)
      real(pr), intent(in out) :: p
      real(pr), intent(in out) :: t

      real(pr), intent(out) :: k(size(z))

      P = 100.0
      T = 200.0

      do while (P > 10)
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
      character(len=:), allocatable, intent(in) :: incipient
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

   subroutine envelope2(ichoice, n, z, T, P, KFACT, & ! This will probably always exist
                        n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, & ! This shouldnt be here in the future
                        this_envelope) ! This output should encapsulate everything
      use dtypes, only: envelope, critical_point
      use linalg, only: point, solve_system
      use constants, only: ouput_path
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
         delS = -0.005
         y = 0.d0
         y(n) = 1.d0
      end if

      Xold = 0.d0
      dFdS = 0.d0
      dFdS(n + 2) = -1.d0

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

         if (run) write(funit_output, *) "SOL", iter, ns, T, P, exp(X(:n))

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

                  its = 0
                  delta = delS

                  ! Variation of lnK based on deltaS
                  m = 15.0_pr * (X - Xold2)/(delta)
                  lnK_extrapolated = (delta) * m(:n) + X(:n)

                  if (all((X(:n) * lnK_extrapolated < 0), dim=1)) then
                     ! All lnK changed sign, so a CP is inminent
                     ! aproach it enough to force the jumping algorithm
                     do while( &
                           maxval(abs(X(:n))) > 0.03 &
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
               print *, "Jumping critical"
               ! approaching the black hole... get out of there! (0.03)
               black_i = black_i + 1
               if (black_i > 50) then
                   print *, "Stuck on the black hole"
                   if (black_i > 100) stop
               end if

               stepX = maxval(abs(X(:n) - Xold(:n))) 

               ! the step given by the most changing logK to fall into 
               ! the black hole
               passingcri = .true.
               if (stepX > 0.07) then
                  !  half step back
                  S = S - delS/2
                  X = X - dXdS*delS/2   
               else
                  ! one more step to jump over the critical point
                  S = S + delS
                  X = X + dXdS*delS   
               end if
            end do
            end block critical_region

            T = exp(X(n + 1))

            if (.not. passingcri .and. abs(T - Told) > 7) then 
               ! Delta T estimations > 7K are not allowed
               delS = delS/2
               S = S - delS
               X = Xold + dXdS*delS
               T = exp(X(n + 1))
            end if

            P = exp(X(n + 2))
            KFACT = exp(X(:n))
            y = z*KFACT

            ! Finish conditions
            if ((dXdS(n + 1)*delS < 0 .and. P < 0.1 .or. T < 120.0) &  ! dew line stops when P<0.1 bar or T<150K
                .or. (P > 1.0 .and. T < 150.0) &   ! bubble line stops when T<150K
                .or. (P > 1500) &
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
end module envelopes
