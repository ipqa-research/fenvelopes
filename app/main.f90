program main
    use envelopes, only: envelope2, max_points, k_wilson_bubble
    use dtypes, only: envelope
    use constants, only: pr
    use io_nml, only: read_system, write_system
    use system, only: z, nc
    implicit none

    real(pr), allocatable :: tv(:) ! Temperatures [K]
    real(pr), allocatable :: pv(:) ! Pressures [bar]
    real(pr), allocatable :: dv(:) ! Pressures [bar]
    real(pr) :: tcri(4) ! Critical points temperatures
    real(pr) :: pcri(4) ! Critical points pressures
    real(pr) :: dcri(4) ! Critical points densities
    real(pr) :: t, p ! Temperature and pressure
    real(pr), allocatable :: k(:) ! K factors
    integer :: n_points, icri(4), ncri, funit_system

    type(envelope) :: env

    allocate(tv(max_points), pv(max_points), dv(max_points), k(max_points))
    call read_system("input.nml")

    open(newunit=funit_system, file="systemdata.nml")
    call write_system(funit_system)
    close(funit_system)

    call k_wilson_bubble(z, t, p, k)

    call envelope2(&
       1, nc, z, T, P, k, &
       n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, &
       env &
    )
    call env%write("output.csv")
end program main
