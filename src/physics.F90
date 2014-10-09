module physics

  use constants,       only: ZERO, REACTION_ABSORBED, REACTION_SCATTERED, &
                             ONE, TWO, PI, MASS_NEUTRON
  use particle_class,  only: Particle
  use nuclide_class,   only: n_nuclides, nuclides
  use random,          only: prn
  use tally_class,     only: tal

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

    ! Save tally
    call tal % add_flux_score(p % get_energy(), dist)

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

    integer :: nuclide_idx
    integer :: reaction_type
    real(8) :: A
    real(8) :: E_in
    real(8) :: E_out
    real(8) :: mass_t
    real(8) :: mu_cm
    real(8) :: speed_in
    real(8) :: speed_in_cm
    real(8) :: uvw_out(3)
    real(8) :: uvw_in(3)
    real(8) :: uvw_in_cm(3)
    real(8) :: v_cm(3)
    real(8) :: v_in(3)
    real(8) :: v_in_cm(3)

    ! Get the reaction type
    reaction_type = p % get_reaction_type()

    ! Perform collision physics for each reaction
    select case (reaction_type)

      case (REACTION_SCATTERED)

        ! Get collision nuclide
        nuclide_idx = p % get_nuclide_index()

        ! Get pre-collision direction and energy
        uvw_in = p % get_uvw()
        E_in = p % get_energy()

        ! Get atomic mass 
        mass_t = nuclides(nuclide_idx) % get_mass()

        ! Calculate atomic weight ratio
        A = mass_t/MASS_NEUTRON

        ! Calculate speed
        speed_in = sqrt(TWO*E_in/MASS_NEUTRON)

        ! Calculate velocity
        v_in = speed_in*uvw_in

        ! Calculate center of mass velocity
        v_cm = v_in/(A + ONE)

        ! Calculate neutron velocity in CM
        v_in_cm = v_in - v_cm

        ! Calculate speed of neutron in CM
        speed_in_cm = sqrt(dot_product(v_in_cm, v_in_cm))

        ! Calculate direction of motion in CM
        uvw_in_cm = v_in_cm / speed_in_cm

        ! Sample isotropic mu in CM
        mu_cm = TWO*prn() - ONE 
  
      case (REACTION_ABSORBED)

        call p % set_alive(.false.)

    end select

  end subroutine collision_physics

end module physics
