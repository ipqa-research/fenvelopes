module inj_envelopes
   use constants, only: pr, R
   use dtypes, only: envelope
   use linalg, only: solve_system

   implicit none

   integer :: max_iters = 500
   integer, parameter :: max_points = 1000
   real(pr), allocatable :: z_0(:)
   real(pr), allocatable :: z_injection(:)
   real(pr) :: T
   real(pr) :: del_S = 0.1
   character(len=:), allocatable :: injection_case
contains
   ! =============================================================================
   !  Injection envelopes
   ! -----------------------------------------------------------------------------
   subroutine F_injection(X, ns, S, F, dF)
      use iso_fortran_env, only: error_unit
      real(pr), intent(in)  :: X(:)
      integer, intent(in)   :: ns
      real(pr), intent(in)  :: S
      real(pr), intent(out) :: F(size(X))
      real(pr), intent(out) :: df(size(x), size(X))

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
          case("dilution")
              z = (z_injection * alpha + z_0)/sum(z_injection * alpha + z_0)
              dzda = -(alpha*z_injection + z_0) &
                     * sum(z_injection) / sum(alpha*z_injection + z_0)**2 &
                     + z_injection / sum(alpha*z_injection + z_0)
          case default
              z = (z_injection * alpha + (1.0_pr - alpha) * z_0)
              dzda = z_injection - z_0
      end select
      
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

   subroutine injection_envelope(X0, spec_number, envels)
       !! Subroutine to calculate Px phase envelopes
       real(pr), intent(in) :: X0(:) !! Vector of variables
       integer,  intent(in) :: spec_number !! Number of specification
       type(envelope), allocatable, intent(out) :: envels(:) !! Calculated envelopes

       real(pr) :: X(size(X0))
       integer :: ns
       real(pr) :: S
       real(pr) :: XS(max_points, size(X0))

       real(pr) :: F(size(X0)), dF(size(X0), size(X0)), dXdS(size(X0))

       integer :: point, iters, n

       X = X0

       n = size(X0) - 2

       ns = spec_number
       S = X(ns)

       print *, "#", X(n+2)
       print *, "X0", iters, ns, X
       do point=1, max_points
          call full_newton(f_injection, iters, X, ns, S, F, dF)

          if (iters >= max_iters) then
             print *, "Breaking due to not converged point"
             exit
          end if

          print *, "SOL", iters, ns, X

          update_spec: block
             real(pr) :: dFdS(size(X0))
             integer  :: ns_new

             dFdS = dF(n+2, :)
             dXdS = solve_system(dF, dFdS)

             ns_new = maxloc(abs(dXdS), dim=1)
             ns_new = ns

             if (ns_new /= ns) then
                dXdS = dXdS/dXdS(ns_new)
                del_S = dXdS(ns_new) * del_S  ! translation of delS to the  new specification variable
             end if
             ns = ns_new

             del_S = sign(1.0_pr, del_S) * minval( [ &
                max(sqrt(abs(X(ns)))/10, 0.1),       &
                abs(del_S) * 3/iters                 &
               ]                                     &
             )
             ! del_S = del_S*10
          end block update_spec

          fix_step: block
             real(pr) :: Xnew(size(X0))
             real(pr) :: dP, dalpha
             Xnew = X + dXdS * del_S
             dP = exp(Xnew(n+1)) - exp(X(n+1))
             dalpha = exp(Xnew(n+2)) - exp(X(n+2))

             if (&
                     abs(dalpha) > 0.1 &
                     .or. abs(dP) > 50 &
                  ) then

                  Xnew = X + dXdS * del_S
                  dP = exp(Xnew(n+1)) - exp(X(n+1))
                  dalpha = exp(Xnew(n+2)) - exp(X(n+2))

                  dXdS = dXdS/50.0_pr
              end if
          end block fix_step

          detect_critical: block
             real(pr) :: K(size(X0) - 2), Knew(size(X0) - 2), fact
             fact = 50
             K = X(:n)
             Knew = X(:n) + fact * dXdS(:n) * del_S
             ! print *, "EXTRAPOL", ns, point, Knew
             if (all(K * Knew < 0)) then
                dXdS = fact * dXdS
             end if
          end block detect_critical

          X = X + dXdS * del_S
          S = X(ns)

          ! if (any(break_conditions(X, ns, S))) exit
       end do
   end subroutine

   subroutine full_newton(fun, iters, X, ns, S, F, dF)
      interface
         subroutine fun(X, ns, S, F, dF)
            import pr
            real(pr), intent(in) :: X(:)
            integer,  intent(in) :: ns
            real(pr), intent(in) :: S
            real(pr), intent(out) :: F(size(X))
            real(pr), intent(out) :: dF(size(X), size(X))
         end subroutine
      end interface
       integer, intent(out)     :: iters
       real(pr), intent(in out) :: X(:)
       integer, intent(in) :: ns
       real(pr), intent(in) :: S
       real(pr), intent(out)    :: F(size(X))
       real(pr), intent(out)    :: df(size(X), size(X))

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
      real(pr) :: X(:)
      integer :: ns
      real(pr) :: S

      logical :: break_conditions(1)

      break_conditions(1) = (X(size(X)) > 1)
   end function
end module
