module physics

  use particle_class,  only: Particle

  implicit none
  private
  public :: collide

contains

!===============================================================================
! COLLIDE
!===============================================================================

  subroutine collide(p)

    type(Particle) :: p

  end subroutine collide

end module physics
