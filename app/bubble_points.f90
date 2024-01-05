program main
   !! Program to calculate bubble points.
   use flap, only: command_line_interface
   use fenvelopes_cli, only: setup_cli => setup_cli_bubble

   implicit none

   type(command_line_interface) :: cli

   ! Setup everything
   call setup
   call run
contains
   subroutine setup
      !! Setup system
      !!
      !! Make output folder (if necessary) and/or clean everyhing in an
      !! existing one. Then read input files to setup needed parameters.
      use constants, only: ouput_path
      use io_nml, only: read_system, write_system
      use inj_envelopes, only: setup_inj => from_nml
      integer :: funit_system, cli_error
      character(len=500) :: infile

      call system("mkdir -p "//trim(ouput_path))
      call system("rm "//trim(ouput_path)//"*")

      call setup_cli(cli)
      call cli%get(val=infile, switch="--infile", error=cli_error)
      call read_system(trim(infile))
   end subroutine

   subroutine run
      use constants, only: pr
      use legacy_ar_models, only: z, nc
      use saturation_points, only: bubble_temperature, EquilibriaState
      implicit none
      real(pr) :: t, p, y0(size(z))
      integer :: i

      type(EquilibriaState) :: equi

      do i=1, 100
         p = real(i, pr)
         equi = bubble_temperature(z, p, 300.0_pr)
         print *, equi%t, equi%p, equi%iters
      end do
   end subroutine
end program main
