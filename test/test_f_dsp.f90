module var
    use constants, only: pr
    real(pr) :: x(16*2+3) = [-1.9505250424097711, &
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
        -0.60278167802332649, &
        -0.33378919660713791, &
        -0.10052077666369962, &
        -0.41963150416019224, &
        -0.22595193414305051, &
        -9.6633793016200953E-002, &
        -9.6893816488812616E-003, &
        3.4503360665734797E-002, &
        0.12298863627470356, &
        0.14608613593991784, &
        0.25088155199109413, &
        0.92977835950126353, &
        2.7781340724066719, &
        3.3719143186445981, &
        3.8606337057567375, &
        5.9948813972520849, &
        log(232.95614892378970), &
        0.0000000000000000, &
        6.3025709559471332 ]


end module

program main
    use var, only: x
    use constants, only: pr
    use inj_envelopes, only: from_nml, z_0
    use dsp_lines, only: dsp_line_F
    use envelopes, only: k_wilson
    use io_nml, only: read_system
    use fordiff, only: derivative
    implicit none
    integer, parameter :: nvars = 16*2 + 3
    character(len=*), parameter :: infile="test/test_f3_pt.nml"

    real(pr) ::F(nvars),df(nvars, nvars)
    real(pr) :: jac(nvars, nvars), numjac(nvars, nvars), S, diff(nvars, nvars)
    integer :: i, j, ns
    
    call read_system(infile)
    call from_nml(infile)
    
    ns = nvars
    S = X(ns)
    call dsp_line_F(X, ns, S, F, dF)
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
    print *, diff(34,16)


contains
    function fun(x)
        real(pr), intent(in) :: x(:)
        real(pr), allocatable :: fun(:)
        call dsp_line_F(X, ns, S, F, dF)
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
            x(i) = x(i) + dx ! *abs(x(i))
            fdx = fun(x)
            x(i) = x(i) - dx
            deriv(:, i) = (fdx - f)/dx
        end do
    end function
end program
