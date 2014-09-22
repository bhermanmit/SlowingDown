module initialize

  use input_xml,  only: read_input_xml
  use nuclide_class,  only: n_nuclides
  use particle_class,  only: Particle

  implicit none
  private
  public :: initialize_run

contains

!===============================================================================
! INITIALIZE_RUN
!===============================================================================

  subroutine initialize_run(p)

    type(Particle) :: p

    ! Read input
    call read_input_xml()

    ! Set up arrays in particle
    call p % initialize(n_nuclides)

  end subroutine initialize_run

end module initialize
