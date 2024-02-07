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

   real(pr) :: tol = 1e-5_pr
   integer :: max_iterations = 100
   real(pr) :: step_tol = 0.1_pr
contains

   type(EquilibriaState) function hpl_pressure(n, t, p0, y0, max_inner_its)
      !! Hpl pressure calculation function.
      !!
      !! Calculates the hpl pressure of a multicomponent mixture.
      use stdlib_optval, only: optval
      real(pr), intent(in) :: n(:) !! Composition vector [moles / molar fraction]
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), optional, intent(in) :: p0 !! Initial pressure [bar]
      real(pr), optional, intent(in) :: y0(:) !! Initial composition
      integer, optional, intent(in) :: max_inner_its(:) !! Inner iterations

      real(pr) :: p, vy, vz

      real(pr) :: k(size(n)), y(size(n)), z(size(n)), lnk(size(n))
      real(pr) :: lnfug_y(size(n)), dlnphi_dp_y(size(n))
      real(pr) :: lnfug_z(size(n)), dlnphi_dp_z(size(n))

      real(pr) :: f, step

      integer :: its, inner_its

      ! =======================================================================
      ! Handle arguments
      ! -----------------------------------------------------------------------
      z = n/sum(n)
      if (size(n) /= nc) call exit(1)
      if (present (p0)) p = p0
      inner_its = optval(inner_its, 50)

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
         call termo(nc, 1, 4, t, p, y, vy, philog=lnfug_y, dlphip=dlnphi_dp_y)
         call termo(nc, 1, 4, t, p, z, vz, philog=lnfug_z, dlphip=dlnphi_dp_z)
         inner_its = 0
         
         do while (any(isnan(lnfug_y)) .and. t > 0)
             inner_its = inner_its + 1
            p = p/2.0_pr
            call termo(nc, 1, 4, t, p, y, vy, philog=lnfug_y, dlphip=dlnphi_dp_y)
            call termo(nc, 1, 4, t, p, z, vz, philog=lnfug_z, dlphip=dlnphi_dp_z)
         end do
         lnk = lnfug_z - lnfug_y
         k = exp(lnk)
         f = sum(z*k) - 1
         step = f/sum(z * k * (dlnphi_dp_z - dlnphi_dp_y))
         do while (abs(step) > abs(step_tol * f))
            step = step/2
         end do
         p = p - step
         y = z*k
         if (abs(step) < tol .and. its > 3) exit
      end do
      hpl_pressure = EquilibriaState(its, y, z, t, p)
      ! =======================================================================
   end function

   type(EquilibriaState) function hpl_temperature(n, p, t0, y0, max_inner_its)
      !! Hpl temperature calculation function.
      !!
      !! Calculates the hpl temperature of a multicomponent mixture.
      use stdlib_optval, only: optval
      real(pr), intent(in) :: n(:) !! Composition vector [moles / molar fraction]
      real(pr), intent(in) :: p !! Temperature [K]
      real(pr), optional, intent(in) :: t0 !! Initial pressure [bar]
      real(pr), optional, intent(in) :: y0(:) !! Initial composition
      integer, optional, intent(in) :: max_inner_its(:) !! Inner iterations

      real(pr) :: t, vy, vz

      real(pr) :: k(size(n)), y(size(n)), z(size(n)), lnk(size(n))
      real(pr) :: lnfug_y(size(n)), dlnphi_dt_y(size(n))
      real(pr) :: lnfug_z(size(n)), dlnphi_dt_z(size(n))

      real(pr) :: f, step

      integer :: its, inner_its

      ! =======================================================================
      ! Handle arguments
      ! -----------------------------------------------------------------------
      z = n/sum(n)
      if (size(n) /= nc) call exit(1)
      if (present (t0)) t = t0
      inner_its = optval(inner_its, 50)

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
         call termo(nc, 1, 4, t, p, y, vy, philog=lnfug_y, dlphit=dlnphi_dt_y)
         call termo(nc, 1, 4, t, p, z, vz, philog=lnfug_z, dlphit=dlnphi_dt_z)
         inner_its = 0
         
         do while (any(isnan(lnfug_y)) .and. t > 0)
             inner_its = inner_its + 1
            t = t - 5.0_pr
            call termo(nc, 1, 4, t, p, y, vy, philog=lnfug_y, dlphit=dlnphi_dt_y)
            call termo(nc, 1, 4, t, p, z, vz, philog=lnfug_z, dlphit=dlnphi_dt_z)
         end do
         lnk = lnfug_z - lnfug_y
         k = exp(lnk)
         f = sum(z*k) - 1
         step = f/sum(z * k * (dlnphi_dt_z - dlnphi_dt_y))
         do while (abs(step) > abs(step_tol * f))
            step = step/2
         end do
         t = t - step
         y = z*k
         if (abs(step) < tol) exit
      end do
      hpl_temperature = EquilibriaState(its, y, z, t, p)
      ! =======================================================================
   end function

   type(EquilibriaState) function bubble_pressure(n, t, p0, y0, max_inner_its)
      !! Bubble pressure calculation function.
      !!
      !! Calculates the bubble pressure of a multicomponent mixture.
      use stdlib_optval, only: optval
      real(pr), intent(in) :: n(:) !! Composition vector [moles / molar fraction]
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), optional, intent(in) :: p0 !! Initial pressure [bar]
      real(pr), optional, intent(in) :: y0(:) !! Initial composition
      integer, optional, intent(in) :: max_inner_its(:) !! Inner iterations

      real(pr) :: p, vy, vz

      real(pr) :: k(size(n)), y(size(n)), z(size(n)), lnk(size(n))
      real(pr) :: lnfug_y(size(n)), dlnphi_dp_y(size(n))
      real(pr) :: lnfug_z(size(n)), dlnphi_dp_z(size(n))

      real(pr) :: f, step

      integer :: its, inner_its

      ! =======================================================================
      ! Handle arguments
      ! -----------------------------------------------------------------------
      z = n/sum(n)
      if (size(n) /= nc) call exit(1)
      if (present (p0)) p = p0
      inner_its = optval(inner_its, 50)

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
         call termo(nc, -1, 4, t, p, y, vy, philog=lnfug_y, dlphip=dlnphi_dp_y)
         call termo(nc, 1, 4, t, p, z, vz, philog=lnfug_z, dlphip=dlnphi_dp_z)
         inner_its = 0
         
         do while (any(isnan(lnfug_y)) .and. t > 0)
             inner_its = inner_its + 1
            p = p/2.0_pr
            call termo(nc, -1, 4, t, p, y, vy, philog=lnfug_y, dlphip=dlnphi_dp_y)
            call termo(nc, 1, 4, t, p, z, vz, philog=lnfug_z, dlphip=dlnphi_dp_z)
         end do
         lnk = lnfug_z - lnfug_y
         k = exp(lnk)
         f = sum(z*k) - 1
         step = f/sum(z * k * (dlnphi_dp_z - dlnphi_dp_y))
         do while (abs(step) > abs(step_tol * f))
            step = step/2
         end do
         p = p - step
         y = z*k
         if (abs(step) < tol) exit
      end do
      bubble_pressure = EquilibriaState(its, y, z, t, p)
      ! =======================================================================
   end function

   type(EquilibriaState) function bubble_temperature(n, p, t0, y0, max_inner_its)
      !! Bubble temperature calculation function.
      !!
      !! Calculates the bubble temperature of a multicomponent mixture.
      use stdlib_optval, only: optval
      real(pr), intent(in) :: n(:) !! Composition vector [moles / molar fraction]
      real(pr), intent(in) :: p !! Temperature [K]
      real(pr), optional, intent(in) :: t0 !! Initial pressure [bar]
      real(pr), optional, intent(in) :: y0(:) !! Initial composition
      integer, optional, intent(in) :: max_inner_its(:) !! Inner iterations

      real(pr) :: t, vy, vz

      real(pr) :: k(size(n)), y(size(n)), z(size(n)), lnk(size(n))
      real(pr) :: lnfug_y(size(n)), dlnphi_dt_y(size(n))
      real(pr) :: lnfug_z(size(n)), dlnphi_dt_z(size(n))

      real(pr) :: f, step

      integer :: its, inner_its

      ! =======================================================================
      ! Handle arguments
      ! -----------------------------------------------------------------------
      z = n/sum(n)
      if (size(n) /= nc) call exit(1)
      if (present (t0)) t = t0
      inner_its = optval(inner_its, 50)

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
         inner_its = 0
         
         do while (any(isnan(lnfug_y)) .and. t > 0)
             inner_its = inner_its + 1
            t = t - 5.0_pr
            call termo(nc, -1, 4, t, p, y, vy, philog=lnfug_y, dlphit=dlnphi_dt_y)
            call termo(nc, 1, 4, t, p, z, vz, philog=lnfug_z, dlphit=dlnphi_dt_z)
         end do
         lnk = lnfug_z - lnfug_y
         k = exp(lnk)
         f = sum(z*k) - 1
         step = f/sum(z * k * (dlnphi_dt_z - dlnphi_dt_y))
         do while (abs(step) > abs(step_tol * f))
            step = step/2
         end do
         t = t - step
         y = z*k
         if (abs(step) < tol) exit
      end do
      bubble_temperature = EquilibriaState(its, y, z, t, p)
      ! =======================================================================
   end function

   type(EquilibriaState) function dew_pressure(n, t, p0, y0, max_inner_its)
      !! Dew pressure calculation function.
      !!
      !! Calculates the dew pressure of a multicomponent mixture.
      use stdlib_optval, only: optval
      real(pr), intent(in) :: n(:) !! Composition vector [moles / molar fraction]
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), optional, intent(in) :: p0 !! Initial pressure [bar]
      real(pr), optional, intent(in) :: y0(:) !! Initial composition
      integer, optional, intent(in) :: max_inner_its(:) !! Inner iterations

      real(pr) :: p, vy, vz

      real(pr) :: k(size(n)), y(size(n)), z(size(n)), lnk(size(n))
      real(pr) :: lnfug_y(size(n)), dlnphi_dp_y(size(n))
      real(pr) :: lnfug_z(size(n)), dlnphi_dp_z(size(n))

      real(pr) :: f, step

      integer :: its, inner_its

      ! =======================================================================
      ! Handle arguments
      ! -----------------------------------------------------------------------
      z = n/sum(n)
      if (size(n) /= nc) call exit(1)
      if (present (p0)) p = p0
      inner_its = optval(inner_its, 50)

      ! Initiliaze with K-wilson factors
      if (present(y0)) then
         y = y0
         k = y/z
      else
         k = k_wilson(t, p)
            y = z / k
      end if
      ! =======================================================================

      ! =======================================================================
      !  Solve point
      ! -----------------------------------------------------------------------
      do its=1,max_iterations
         call termo(nc, 1, 4, t, p, y, vy, philog=lnfug_y, dlphip=dlnphi_dp_y)
         call termo(nc, -1, 4, t, p, z, vz, philog=lnfug_z, dlphip=dlnphi_dp_z)
         inner_its = 0
         
         do while (any(isnan(lnfug_y)) .and. t > 0)
             inner_its = inner_its + 1
            p = p/2.0_pr
            call termo(nc, 1, 4, t, p, y, vy, philog=lnfug_y, dlphip=dlnphi_dp_y)
            call termo(nc, -1, 4, t, p, z, vz, philog=lnfug_z, dlphip=dlnphi_dp_z)
         end do
         lnk = lnfug_z - lnfug_y
         k = exp(lnk)
         f = sum(z*k) - 1
         step = f/sum(z * k * (dlnphi_dp_z - dlnphi_dp_y))
         do while (abs(step) > abs(step_tol * f))
            step = step/2
         end do
         p = p - step
         y = z*k
         if (abs(step) < tol) exit
      end do
      dew_pressure = EquilibriaState(its, y, z, t, p)
      ! =======================================================================
   end function

   type(EquilibriaState) function dew_temperature(n, p, t0, y0, max_inner_its)
      !! Dew temperature calculation function.
      !!
      !! Calculates the dew temperature of a multicomponent mixture.
      use stdlib_optval, only: optval
      real(pr), intent(in) :: n(:) !! Composition vector [moles / molar fraction]
      real(pr), intent(in) :: p !! Temperature [K]
      real(pr), optional, intent(in) :: t0 !! Initial pressure [bar]
      real(pr), optional, intent(in) :: y0(:) !! Initial composition
      integer, optional, intent(in) :: max_inner_its(:) !! Inner iterations

      real(pr) :: t, vy, vz

      real(pr) :: k(size(n)), y(size(n)), z(size(n)), lnk(size(n))
      real(pr) :: lnfug_y(size(n)), dlnphi_dt_y(size(n))
      real(pr) :: lnfug_z(size(n)), dlnphi_dt_z(size(n))

      real(pr) :: f, step

      integer :: its, inner_its

      ! =======================================================================
      ! Handle arguments
      ! -----------------------------------------------------------------------
      z = n/sum(n)
      if (size(n) /= nc) call exit(1)
      if (present (t0)) t = t0
      inner_its = optval(inner_its, 50)

      ! Initiliaze with K-wilson factors
      if (present(y0)) then
         y = y0
         k = y/z
      else
         k = k_wilson(t, p)
            y = z / k
      end if
      ! =======================================================================

      ! =======================================================================
      !  Solve point
      ! -----------------------------------------------------------------------
      do its=1,max_iterations
         call termo(nc, 1, 4, t, p, y, vy, philog=lnfug_y, dlphit=dlnphi_dt_y)
         call termo(nc, -1, 4, t, p, z, vz, philog=lnfug_z, dlphit=dlnphi_dt_z)
         inner_its = 0
         
         do while (any(isnan(lnfug_y)) .and. t > 0)
             inner_its = inner_its + 1
            t = t - 5.0_pr
            call termo(nc, 1, 4, t, p, y, vy, philog=lnfug_y, dlphit=dlnphi_dt_y)
            call termo(nc, -1, 4, t, p, z, vz, philog=lnfug_z, dlphit=dlnphi_dt_z)
         end do
         lnk = lnfug_z - lnfug_y
         k = exp(lnk)
         f = sum(z*k) - 1
         step = f/sum(z * k * (dlnphi_dt_z - dlnphi_dt_y))
         do while (abs(step) > abs(step_tol * f))
            step = step/2
         end do
         t = t - step
         y = z*k
         if (abs(step) < tol) exit
      end do
      dew_temperature = EquilibriaState(its, y, z, t, p)
      ! =======================================================================
   end function
end module
