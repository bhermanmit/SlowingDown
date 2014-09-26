module physics

  use constants,       only: ZERO, REACTION_ABSORBED, REACTION_SCATTERED, &
                             ONE, TWO, PI
  use particle_class,  only: Particle
  use nuclide_class,   only: n_nuclides, nuclides
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

    integer :: nuclide_idx
    integer :: reaction_type
    real(8) :: A
    real(8) :: coszz
    real(8) :: coszt
    real(8) :: Ein
    real(8) :: Eout
    real(8) :: uvw(3)
    real(8) :: uvw0(3)
    real(8) :: phi
    real(8) :: s1
    real(8) :: s2
    real(8) :: term1

    ! Get the reaction type
    reaction_type = p % get_reaction_type()

    ! Perform collision physics for each reaction
    select case (reaction_type)

      case (REACTION_SCATTERED)

        ! Get collision nuclide
        nuclide_idx = p % get_nuclide_index()

        ! Get pre-collision direction
        uvw = p % get_uvw()

        ! Get atomic weight
        A = nuclides(nuclide_idx) % get_A()

        ! Sample isotropic azimuthal angle in COM
        phi = 2.0*PI*prn()

        ! Sample isotropic polar angle in COM
        coszz = ONE - TWO*prn()
        term1 = A**2 + ONE + TWO*A*coszz

        ! Transform scattering polar angle into LAB
        coszt = (ONE + A*coszz)/sqrt(term1)
        s1 = sqrt(ONE - coszt**2)
        s2 = sqrt(ONE - uvw0(1)**2)

        ! Transform cosines relative to incoming flight path
        uvw(1) = coszt*uvw0(1) + (s1*(uvw0(1)*uvw0(3)*cos(phi) - &
             uvw(2)*sin(phi))/s2)
        uvw(2) = coszt*uvw0(2) + (s1*(uvw0(2)*uvw0(3)*cos(phi) + &
             uvw(1)*sin(phi))/s2)
        uvw(3) = coszt*uvw0(3) - s1*s2*cos(phi)

        ! Calculate outgoing energy
        Ein = p % get_energy()
        Eout = Ein*term1/(A + ONE)**2

        ! Set new direction and energy of particle
        call p % set_uvw(uvw)
        call p  % set_energy(Eout)
  
      case (REACTION_ABSORBED)

        call p % set_alive(.false.)

    end select

  end subroutine collision_physics

end module physics
