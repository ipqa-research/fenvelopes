program main
    use constants, only: pr
    use envelopes, only: pt_F_three_phases
    use io_nml, only: read_system
    use fordiff, only: derivative
    character(len=*), parameter :: infile="test/test_f3_pt.nml"

    real(pr), allocatable :: X(:), F(:), Fdx(:), df(:, :)
    real(pr), allocatable :: jac(:, :), numjac(:, :)
    real(pr) :: dx, S
    integer :: i, j, ns, len_x
    
    call read_system(infile)
    X=[&
        -1.9513638490072251, &
        -1.1661779375906842, &
        -0.49956967077471320, &
        -1.4404970885168327, &
        -1.0808457388533852, &
        -0.93333583839794887, &
        -0.90754166039592443, &
        -0.78304628593831593, &
        -0.74714660089597174, &
        -0.70952247609720753, &
        -0.69457992914970657, &
        4.3520956423480306E-002, &
        2.7892744922014709, &
        3.9246051762741589, &
        4.5125746361069288, &
        11.378182771852959, &
        -0.60345867918897100, &
        -0.33364362813220155, &
        -9.9834721980533600E-002, &
        -0.41984716093079832, &
        -0.22562098387844687, &
        -9.5953395560971058E-002, &
        -8.8056526984703620E-003, &
        3.5524493620038851E-002, &
        0.12423602823851601, &
        0.14742579483929957, &
        0.25249775329708163, &
        0.93379560927496952, &
        2.7876678324424828, &
        3.3829857970183483, &
        3.8713409641216381, &
        6.0124425598888385, &
        5.4510132751728753, &
        6.3021529774269105, &
        0.6000000000000000 ]

    len_x = size(X)
    allocate(F(len_x), Fdx(len_x), df(len_x, len_x), jac(len_x, len_x), numjac(len_x, len_x) )

    ns = len_x
    S = X(ns)
    call pt_F_three_phases(X, ns, S, F, dF)
    jac = df

    numjac = derivative(f=fun, x=X, h=1e-5_pr, method='central')
    
    do ns=1, len_x
        print *, ns
        do i=1, len_x
            print "(I3, x, 3(E10.3, 2x))", i, jac(i,ns), numjac(i,ns), (numjac(i, ns) - jac(i, ns))
        end do
    end do

    print *, ""
    print *, maxloc(abs(numjac - jac))

    ! if (any(abs(numjac - jac) > 1e-2)) call exit(1)

contains
    function fun(x)
        real(pr), intent(in) :: x(:)
        real(pr), allocatable :: fun(:)
        call pt_F_three_phases(X, ns, S, F, dF)
        fun = F
    end function
end program
