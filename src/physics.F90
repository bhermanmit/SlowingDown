module physics

  use constants,       only: ZERO, REACTION_ABSORBED, REACTION_SCATTERED
  use particle_class,  only: Particle
  use nuclide_class,   only: n_nuclides
  use random,          only: prn

  implicit none
  private
  public :: run_physics

contains

!===============================================================================
! RUN_PHYSICS
!===============================================================================

  subroutine run_physics(p)

    type(Particle), intent(inout) :: p

    ! Determine how far the neutron travels
    call sample_pathlength(p)

    ! Determine what isotope the neutron interacts with
    call sample_isotope(p)

    ! Determine the reaction that occurs
    call sample_reaction(p)

    ! Perform collision physics
    call collision_physics(p)

  end subroutine run_physics

!===============================================================================
! SAMPLE_PATHLENGTH
!===============================================================================

  subroutine sample_pathlength(p)

    type(Particle), intent(inout) :: p

    real(8) :: dist

    ! Sample path length to collision
    dist = -log(prn())/p % get_macro_total()

    ! Save distance particle traveled
    call p % set_distance(dist)

  end subroutine sample_pathlength

!===============================================================================
! SAMPLE_ISOTOPE
!===============================================================================

  subroutine sample_isotope(p)

    type(Particle), intent(inout) :: p

    integer :: nuclide_idx
    real(8) :: cutoff
    real(8) :: prob

    ! Calculate cutoff for random sampling
    cutoff = prn()*p % get_macro_total()

    ! Loop around isotopes
    nuclide_idx = 0
    prob = ZERO
    do while (prob < cutoff)
      nuclide_idx = nuclide_idx + 1 
      prob = prob + p % get_nuclide_macroxs_t(nuclide_idx)
    end do

    ! Save in particle
    call p % set_nuclide_index(nuclide_idx)

  end subroutine sample_isotope

!===============================================================================
! SAMPLE_REACTION
!===============================================================================

  subroutine sample_reaction(p)

    type(Particle), intent(inout) :: p

    integer :: nuclide_idx
    real(8) :: cutoff

    ! Get the nuclide index
    nuclide_idx = p % get_nuclide_index()

    ! Calculate cutoff
    cutoff = prn() * p % get_nuclide_macroxs_t(nuclide_idx)

    ! Check against absorption xs
    if (cutoff < p % get_nuclide_macroxs_a(nuclide_idx)) then
      call p % set_reaction_type(REACTION_ABSORBED)
    else
      call p % set_reaction_type(REACTION_SCATTERED)
    end if

  end subroutine sample_reaction

!===============================================================================
! COLLISION_PHYSICS 
!===============================================================================

  subroutine collision_physics(p)

    type(Particle), intent(inout) :: p

  end subroutine collision_physics

end module physics
