module tracking

  use particle_class,  only: Particle

  implicit none
  private
  public :: transport

contains

!===============================================================================
! TRANSPORT
!===============================================================================

  subroutine transport(p)

    type(Particle) :: p

  end subroutine transport

end module tracking 
