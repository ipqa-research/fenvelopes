program main
   use envelopes, only: envelope2, max_points, k_wilson_bubble, &
                        max_points, p_wilson, k_wilson
   use dtypes, only: envelope
   use constants, only: pr
   use system, only: z, nc

   implicit none

   real(pr), allocatable :: tv(:) ! Temperatures [K]
   real(pr), allocatable :: pv(:) ! Pressures [bar]
   real(pr), allocatable :: dv(:) ! Densities [mol/L]

   real(pr) :: tcri(4)            ! Critical points temperatures
   real(pr) :: pcri(4)            ! Critical points pressures
   real(pr) :: dcri(4)            ! Critical points densities

   real(pr) :: t, p               ! Temperature and pressure
   real(pr), allocatable :: k(:)  ! K factors

   integer :: n_points, icri(4), ncri, i

   type(envelope) :: bub_env, dew_env

   call setup              !
   call pt_envelopes
   call px_envelopes

    call k_wilson_bubble(z, t, p, k)

    call envelope2(&
       1, nc, z, T, P, k, &
       n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, &
       env &
    )
    call env%write("output.csv")
end program main
