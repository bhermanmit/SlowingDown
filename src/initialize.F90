module initialize

  use input_xml,  only: read_input_xml

  implicit none
  private
  public :: initialize_run

contains

!===============================================================================
! INITIALIZE_RUN
!===============================================================================

  subroutine initialize_run()

    ! Read input
    call read_input_xml()

  end subroutine initialize_run

end module initialize
