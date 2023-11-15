module pt_2ph
    use constants, only: pr
    real(pr) :: x(16 + 2) = [&
        -1.9505250424097711, &
        -1.1667253258630770, &
        -0.50083650782086886, &
        -1.4403476424022568, &
        -1.0811890481373283, &
        -0.93383702132049806, &
        -0.90796780830843316, &
        -0.78367628417639512, &
        -0.74775960561138821, &
        -0.71024720094329619, &
        -0.69526113124061351, &
        4.0543875795062009E-002, &
        2.7816909475822125, &
        3.9157194615109057, &
        4.5068354677502320, &
        11.360614772439938, &
        6.3025709559471332, & 
        log(232.95614892378970)]
end module

program main
    use legacy_ar_models, only: z
    use pt_2ph, only: x
    use constants, only: pr
    use inj_envelopes, only: from_nml, z_0
    use envelopes, only: F2
    use io_nml, only: read_system
    use fordiff, only: derivative
    implicit none
    integer, parameter :: nvars = 16 + 2
    real(pr) :: y(16)
    character(len=*), parameter :: infile="test/test_f3_pt.nml"
    character(len=:), allocatable :: phase

    real(pr) ::F(nvars),df(nvars, nvars)
    real(pr) :: jac(nvars, nvars), numjac(nvars, nvars), S, diff(nvars, nvars)
    integer :: i, j, ns

    call read_system(infile)
    call from_nml(infile)
    
    ns = nvars
    S = X(ns)

    phase = "liquid"
    call F2(phase, z, y, X, S, ns, F, dF)
    jac = df
    numjac = deriv(x)
    diff = (jac - numjac)/numjac
    
    do ns=1, nvars
        print *, ns
        do i=1, nvars
            print "(I3, x, 3(E15.5, 2x))", i, jac(i,ns), numjac(i,ns), diff(i, ns)
        end do
    end do

    print *, maxloc(abs(diff))


contains
    function fun(x)
        real(pr), intent(in) :: x(:)
        real(pr), allocatable :: fun(:)
        call F2(phase, z, y, X, S, ns, F, dF)
        fun = F
    end function

    function deriv(x)
        real(pr) :: x(:)
        real(pr) :: deriv(size(x), size(x))
        integer :: i
        real(pr) :: dx=1e-2

        real(pr) :: fdx(size(x)), f(size(x))

        do i=1,size(x)
            f = fun(x)
            x(i) = x(i) + dx*abs(x(i))
            fdx = fun(x)
            x(i) = x(i) - dx
            deriv(:, i) = (fdx - f)/dx
        end do
    end function
end program
