program main
   use inj_envelopes, only: get_z, setup_inj => from_nml
   use constants, only: pr, ouput_path
   use legacy_ar_models, only: nc, z
   use flap, only: command_line_interface

   implicit none
   real(pr) :: et, st

   type(command_line_interface) :: cli
   integer :: cli_error
   character(len=99) :: cli_string

   real(pr) :: alpha=0.0
   real(pr) :: temperature

   ! Setup everything
   call setup
contains
   subroutine setup_cli
      !! Setup CLI subroutine
      !!
      !! Setup the Command-Line-Interface processor
      call cli%init(progname="isotherm", description="Isotherm calculation")
      call cli%add( &
         switch="--infile", &
         switch_ab="-i",    &
         help="Input file", &
         error=cli_error,   &
         required=.true.    &
      )
      
      call cli%add( &
         switch="--temperature", &
         switch_ab="-t", &
         help="Temperature [K]", &
         error=cli_error, &
         required=.true.)
      call cli%parse(error=cli_error)

      if (cli_error /= 0) stop
   end subroutine

   subroutine setup
      !! Setup system
      !!
      !! Make output folder (if necessary) and/or clean everyhing in an
      !! existing one. Then read input files to setup needed parameters.
      use io_nml, only: read_system, write_system
      use inj_envelopes, only: setup_inj => from_nml
      integer :: funit_system
      character(len=500) :: infile

      call system("mkdir -p "//trim(ouput_path))
      call system("rm "//trim(ouput_path)//"*")

      call setup_cli
      call cli%get(val=infile, switch="--infile", error=cli_error)

      call read_system(trim(infile))
      call setup_inj(trim(infile))
      call get_z(alpha, z)

      open (newunit=funit_system, file="systemdata.nml")
      call write_system(funit_system)
      close (funit_system)
   end subroutine

   subroutine calc_isotherm
      call cli%get(val=temperature, switch="--temperature", error=cli_error)
   end subroutine
end program main
