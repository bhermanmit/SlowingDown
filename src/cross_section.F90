module cross_section

  use particle_class,  only: Particle

  implicit none

contains

!===============================================================================
! CALCULATE_XS
!===============================================================================

  subroutine calculate_xs(p)

    type(Particle) :: p

  end subroutine calculate_xs

end module cross_section
