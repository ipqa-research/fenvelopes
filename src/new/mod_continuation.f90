!! Module that implements a generic implementation of the continuation method
!! to calculate phase envelopes.
module numerical_continuation_mod
    use yaeos_constants, only: pr
    implicit none

    abstract interface
        subroutine continuation_fun(X, ns, S, F, dF)
            import pr
            real(pr), intent(in) :: X
            integer, intent(in) :: ns
            real(pr), intent(in) :: S
            real(pr), intent(out) :: F(size(X))
            real(pr), intent(out) :: dF(size(X), size(X))
        end subroutine
    end interface

contains

    subroutine numerical_continuation(&
        fun, &
        init_point, init_spec, init_del_spec &
        )
        procedure(continuation_fun) :: fun
        real(pr), intent(in) :: init_point(:)
        integer,  intent(in) :: init_spec
        real(pr), intent(in) :: init_del_spec
    end subroutine
end module