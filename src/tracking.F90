module tracking

  use particle_class,  only: Particle
  use random,          only: prn

  implicit none
  private
  public :: transport

contains

!===============================================================================
! TRANSPORT samples path length and saves it for later tallying
!===============================================================================

  subroutine transport(p)

    type(Particle) :: p

    real(8) :: dist

    ! Sample path length to collision
    dist = -log(prn())/p % get_macro_total()

    ! Save distance particle traveled
    call p % set_distance(dist)

  end subroutine transport

end module tracking 
