! module numerical_continuation_mod
! !! Module that implements a generic implementation of the continuation method
! !! to calculate phase envelopes.
!     use yaeos_constants, only: pr
!     implicit none
! 
!     integer, parameter :: max_points=1000
! 
!     type, abstract :: continuation_method
!         real(pr), allocatable :: x(:)
!         integer :: ns
!         real(pr) :: del_S
!     contains
!         procedure, deferred :: f
!         procedure, deferred :: df
!     end type
! 
!     abstract interface
!         subroutine continuation_fun(X, ns, S, F, dF)
!             import pr
!             real(pr), intent(in) :: X
!             integer, intent(in) :: ns
!             real(pr), intent(in) :: S
!             real(pr), intent(out) :: F(size(X))
!             real(pr), intent(out) :: dF(size(X), size(X))
!         end subroutine
!     end interface
! 
! contains
! 
!     subroutine numerical_continuation(&
!         fun, &
!         init_point, init_spec, init_del_spec, &
!         XS &
!     )
!         procedure(continuation_fun) :: fun
!         real(pr), intent(in) :: init_point(:)
!         integer,  intent(in) :: init_spec
!         real(pr), intent(in) :: init_del_spec
!         real(pr), intent(out) :: XS(max_points, size(init_point))
! 
!         ! Inner variables
!         real(pr) :: X(size(init_point))
!         integer :: ns
!         real(pr) :: S
!         real(pr) :: del_S
!         real(pr) :: dXdS(size(init_point))
! 
!         real(pr) :: F(size(init_point)), dF(size(init_point), size(init_point))
! 
!         integer :: point !! Main loop
! 
!         ! Point solving variables
!         integer :: solve_its
! 
!         ! Initialize variables
!         X = init_point
!         ns = init_spec
!         S = X(ns)
!         del_S = init_del_spec
! 
! 
!         main_loop: do point = 1, max_points
!             call solve_point(fun, solve_its, X, ns, S)
!             XS(point, :) = X
!             call update_specification(X, ns, del_S, dF, solve_its, dXdS)
!             call fix_step(X, ns, S, solve_its, del_S, dXdS)
!             call post_process(point, XS, ns, S, del_S, dXdS)
!             if (any(break_conditions(X, ns, S, del_S))) exit main_loop
!             X = X + dXdS * del_S
!             S = X(ns)
!         end do main_loop
!     end subroutine
! end module