module fenvelopes_cli
   !! Module for setting up command line interfaces for different programs
   use flap, only: command_line_interface
contains
   subroutine setup_cli_envelopes(cli)
      !! Setup CLI subroutine for phase envelopes calculation
      !!
      !! Setup the Command-Line-Interface processor
      type(command_line_interface), intent(in out) :: cli
      integer :: cli_error

      call cli%init(progname="envelopes", description="Phase Envelopes")
      call cli%add( &
         switch="--infile", &
         switch_ab="-i", &
         help="Input file", &
         error=cli_error, &
         required=.true.)

      call cli%add( &
         switch="--injection", &
         switch_ab="-px", &
         help="Trace Px lines", &
         error=cli_error, &
         required=.false., &
         def=".false.")
      
      call cli%add( &
         switch="--three_phases", &
         switch_ab="-3ph", &
         help="Trace three phase lines", &
         error=cli_error, &
         required=.false., &
         def=".false.")

      call cli%add( &
         switch="--alpha0", &
         switch_ab="-a", &
         help="Initial alpha0, used for PT envelopes", &
         error=cli_error, &
         required=.false., &
         def="0")
      call cli%parse(error=cli_error)
      print *, cli_error

      if (cli_error /= 0) stop
   end subroutine

   subroutine setup_cli_isotherms(cli)
      !! Setup CLI subroutine for isotherms calculations
      !!
      !! Setup the Command-Line-Interface processor
      type(command_line_interface), intent(in out) :: cli
      integer :: cli_error
   end subroutine
   
   subroutine setup_cli_flash(cli)
      !! Setup CLI subroutine for flash calculations
      !!
      !! Setup the Command-Line-Interface processor
      type(command_line_interface), intent(in out) :: cli
      integer :: cli_error
   end subroutine
   
   subroutine setup_cli_bubble(cli)
      !! Setup CLI subroutine for bubble points
      !!
      !! Setup the Command-Line-Interface processor
      type(command_line_interface), intent(in out) :: cli
      integer :: cli_error

      call cli%init(progname="bubble_point", description="Bubble point")
      call cli%add( &
         switch="--infile", &
         switch_ab="-i", &
         help="Input file", &
         error=cli_error, &
         required=.true.)

      if (cli_error /= 0) stop
   end subroutine
end module
