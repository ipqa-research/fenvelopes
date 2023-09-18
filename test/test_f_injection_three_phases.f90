program main
    use constants, only: pr
    use inj_envelopes, only: F_injection_three_phases, from_nml, z_0
    use envelopes, only: k_wilson
    use io_nml, only: read_system
    use fordiff, only: derivative
    integer, parameter :: nvars = 16+3
    character(len=*), parameter :: infile="test/test_f3.nml"

    real(pr) :: X(nvars), F(nvars), Fdx(nvars), df(nvars, nvars), dx
    real(pr) :: jac(nvars, nvars), numjac(nvars, nvars), S
    integer :: i, j, ns
    
    call read_system(infile)
    call from_nml(infile)
    X = [-0.41981837560016311, & 
        -0.79828389545752954, & 
        -0.93607212819652841, & 
        -1.1294431131109890, & 
        -1.2738268650266098, & 
        -0.37650048779143636, & 
        -1.0324858608887026, & 
        4.5238362455206769, & 
        0.29641441051264411, & 
        0.62818018210875071, & 
        7.9138737267353043E-002, & 
        -0.25915674873115258, & 
        -0.47420214799519028, & 
        -0.72823939039186159, & 
        -4.5112757160826309, & 
        -16.235569521702651, & 
        4.8671594176813775, & 
        0.30843430279212669, & 
        1.0000000474974513E-003]

    ns = nvars
    S = X(ns)
    call F_injection_three_phases(X, ns, S, F, dF)
    jac = df

    numjac = derivative(f=fun, x=X, h=1e-15_pr, method='central')
    
    do ns=1, nvars
        print *, ns
        do i=1, nvars
            print "(I3, x, 3(E10.3, 2x))", i, jac(i,ns), numjac(i,ns), (numjac(i, ns) - jac(i, ns)) * 100
        end do
    end do

contains
    function fun(x)
        real(pr), intent(in) :: x(:)
        real(pr), allocatable :: fun(:)
        call F_injection_three_phases(X, ns, S, F, dF)
        fun = F
    end function
end program
