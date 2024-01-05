module saturation_points
   use constants, only: pr
   use legacy_ar_models, only: termo, nc
   use envelopes, only: k_wilson
   implicit none

   type :: EquilibriaState
      integer :: iters !! Iterations needed to reach the state
      real(pr), allocatable :: y(:) !! Vapour molar fractions
      real(pr), allocatable :: x(:) !! Liquid molar fractions
      real(pr) :: t !! Temperature [K]
      real(pr) :: p !! Pressure [bar]
   end type

   real(pr) :: tol = 1e-5
   integer :: max_iterations = 100
contains
   type(EquilibriaState) function bubble_pressure(n, t, p0, y0)
      !! Bubble pressure calculation function.
      !!
      !! Calculates the bubble pressure of a multicomponent mixture at a given
      !! composition `n` and temperature `t`
      real(pr), intent(in) :: n(:) !! Composition vector [moles / molar fraction]
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), optional, intent(in) :: p0 !! Initial pressure [bar]
      real(pr), optional, intent(in) :: y0 !! Initial composition

      real(pr) :: p, vy, vz

      real(pr) :: k(size(n)), y(size(n)), z(size(n)), lnk(size(n))
      real(pr) :: lnfug_y(size(n)), dlnphi_dt_y(size(n)), dlnphi_dp_y(size(n))
      real(pr) :: lnfug_z(size(n)), dlnphi_dt_z(size(n)), dlnphi_dp_z(size(n))

      real(pr) :: f, dfdp

      integer :: its

      ! =======================================================================
      ! Handle arguments
      ! -----------------------------------------------------------------------
      z = n/sum(n)
      if (size(n) /= nc) call exit(1)
      if (present (p0)) p = p0
      
      ! Initiliaze with K-wilson factors
      if (present(y0)) then
         y = y0
         k = y/z
      else
         k = k_wilson(t, p)
         y = k * z
      end if
      print *, y
      ! =======================================================================

      ! =======================================================================
      !  Solve point
      ! -----------------------------------------------------------------------
      do its=1,max_iterations
         call termo(nc, -1, 4, t, p, y, vy, philog=lnfug_y, dlphip=dlnphi_dp_y)
         call termo(nc, 1, 4, t, p, z, vz, philog=lnfug_z, dlphip=dlnphi_dp_z)
         do while (any(isnan(lnfug_y)))
         p = p/2.0_pr
            call termo(nc, -1, 4, t, p, y, vy, philog=lnfug_y, dlphip=dlnphi_dp_y)
            call termo(nc, 1, 4, t, p, z, vz, philog=lnfug_z, dlphip=dlnphi_dp_z)
         end do
         lnk = lnfug_z - lnfug_y
         k = exp(lnk)
         f = sum(z*k) - 1
         dfdp = sum(z * k * (dlnphi_dp_z - dlnphi_dp_y))
         p = p - f/dfdp
         y = z*k
         if (abs(f/dfdp) < tol) exit
      end do

      bubble_pressure = EquilibriaState(its, y, z, t, p)
      ! =======================================================================
   end function
   
   type(EquilibriaState) function bubble_temperature(n, p, t0, y0)
      !! Bubble temperature calculation function.
      !!
      !! Calculates the bubble temperature of a multicomponent mixture at a given
      !! composition `n` and temperature `t`
      real(pr), intent(in) :: n(:) !! Composition vector [moles / molar fraction]
      real(pr), intent(in) :: p !! Pressure [bar]
      real(pr), optional, intent(in) :: t0 !! Initial temperature [K]
      real(pr), optional, intent(in) :: y0(:) !! Initial composition

      real(pr) :: t, vy, vz

      real(pr) :: k(size(n)), y(size(n)), z(size(n)), lnk(size(n))
      real(pr) :: lnfug_y(size(n)), dlnphi_dt_y(size(n))
      real(pr) :: lnfug_z(size(n)), dlnphi_dt_z(size(n))

      real(pr) :: f, dfdt

      integer :: its

      ! =======================================================================
      ! Handle arguments
      ! -----------------------------------------------------------------------
      z = n/sum(n)
      if (size(n) /= nc) call exit(1)
      if (present (t0)) t = t0
      
      ! Initiliaze with K-wilson factors
      if (present(y0)) then
         y = y0
         k = y/z
      else
         k = k_wilson(t, p)
         y = k * z
      end if
      ! =======================================================================

      ! =======================================================================
      !  Solve point
      ! -----------------------------------------------------------------------
      do its=1,max_iterations
         call termo(nc, -1, 4, t, p, y, vy, philog=lnfug_y, dlphit=dlnphi_dt_y)
         call termo(nc, 1, 4, t, p, z, vz, philog=lnfug_z, dlphit=dlnphi_dt_z)
         do while(any(isnan(y)) .and. t > 0)
            t = t - 10
            call termo(nc, -1, 4, t, p, y, vy, philog=lnfug_y, dlphit=dlnphi_dt_y)
            call termo(nc, 1, 4, t, p, z, vz, philog=lnfug_z, dlphit=dlnphi_dt_z)
         end do
         lnk = lnfug_z - lnfug_y
         k = exp(lnk)
         f = sum(z*k) - 1
         dfdt = sum(z * k * (dlnphi_dt_z - dlnphi_dt_y))
         t = t - f/dfdt
         y = z*k
         if (abs(f/dfdt) < tol) exit
      end do

      bubble_temperature = EquilibriaState(its, y, z, t, p)
      ! =======================================================================
   end function
end module
