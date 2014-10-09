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
    real(8) :: mu
    real(8) :: mu_cm
    real(8) :: speed_in
    real(8) :: speed_out
    real(8) :: speed_in_cm
    real(8) :: speed_out_cm
    real(8) :: uvw_in(3)
    real(8) :: uvw_out(3)
    real(8) :: uvw_in_cm(3)
    real(8) :: uvw_out_cm(3)
    real(8) :: v_cm(3)
    real(8) :: v_in(3)
    real(8) :: v_out(3)
    real(8) :: v_in_cm(3)
    real(8) :: v_out_cm(3)

    ! Get the reaction type
    reaction_type = p % get_reaction_type()

    ! Perform collision physics for each reaction
    select case (reaction_type)

      case (REACTION_SCATTERED)

        ! Get collision nuclide
        nuclide_idx = p % get_nuclide_index()

        ! Get pre-collision direction and energy in LAB
        uvw_in = p % get_uvw()
        E_in = p % get_energy()

        ! Get atomic mass 
        mass_t = nuclides(nuclide_idx) % get_mass()

        ! Calculate atomic weight ratio
        A = mass_t/MASS_NEUTRON

        ! Calculate incoming speed in LAB
        speed_in = sqrt(E_in)

        ! Calculate incoming velocity in LAB
        v_in = speed_in*uvw_in

        ! Calculate center of mass velocity
        v_cm = v_in/(A + ONE)

        ! Calculate incoming neutron velocity in CM
        v_in_cm = v_in - v_cm

        ! Calculate incoming speed of neutron in CM
        speed_in_cm = sqrt(dot_product(v_in_cm, v_in_cm))

        ! Calculate incoming direction of motion in CM
        uvw_in_cm = v_in_cm / speed_in_cm

        ! Sample isotropic mu in CM
        mu_cm = TWO*prn() - ONE 

        ! Rotate angle in CM 
        uvw_out_cm = rotate_angle(uvw_in_cm, mu_cm)

        ! In elastic scattering, speed is conserved in CM
        speed_out_cm = speed_in_cm

        ! Calculate outgoing velocity in CM
        v_out_cm = speed_out_cm * uvw_out_cm

        ! Calculate outgoing velocity in LAB
        v_out = v_out_cm + v_cm
 
        ! Calculate outgoing energy in lAB
        E_out = dot_product(v_out, v_out)

        ! Calculate outgoing speed in LAB
        speed_out = sqrt(E_out)

        ! Calculate outgoing direction in LAB
        uvw_out = v_out / speed_out

        ! Calculate mu in LAB
        mu = dot_product(uvw_in, uvw_out)

        ! Back info back in particle
        call p % set_energy(E_out)
        call p % set_uvw(uvw_out)

      case (REACTION_ABSORBED)

        call p % set_alive(.false.)

    end select

  end subroutine collision_physics

!===============================================================================
! ROTATE_ANGLE rotates direction cosines through a polar angle whose cosine is
! mu and through an azimuthal angle sampled uniformly. Note that this is done
! with direct sampling rather than rejection as is done in MCNP and SERPENT.
!===============================================================================

  function rotate_angle(uvw0, mu) result(uvw)

    real(8), intent(in) :: uvw0(3) ! directional cosine
    real(8), intent(in) :: mu      ! cosine of angle in lab or CM
    real(8)             :: uvw(3)  ! rotated directional cosine

    real(8) :: phi    ! azimuthal angle
    real(8) :: sinphi ! sine of azimuthal angle
    real(8) :: cosphi ! cosine of azimuthal angle
    real(8) :: a      ! sqrt(1 - mu^2)
    real(8) :: b      ! sqrt(1 - w^2)
    real(8) :: u0     ! original cosine in x direction
    real(8) :: v0     ! original cosine in y direction
    real(8) :: w0     ! original cosine in z direction

    ! Copy original directional cosines
    u0 = uvw0(1)
    v0 = uvw0(2)
    w0 = uvw0(3)

    ! Sample azimuthal angle in [0,2pi)
    phi = TWO * PI * prn()

    ! Precompute factors to save flops
    sinphi = sin(phi)
    cosphi = cos(phi)
    a = sqrt(max(ZERO, ONE - mu*mu))
    b = sqrt(max(ZERO, ONE - w0*w0))

    ! Need to treat special case where sqrt(1 - w**2) is close to zero by
    ! expanding about the v component rather than the w component
    if (b > 1e-10) then
      uvw(1) = mu*u0 + a*(u0*w0*cosphi - v0*sinphi)/b
      uvw(2) = mu*v0 + a*(v0*w0*cosphi + u0*sinphi)/b
      uvw(3) = mu*w0 - a*b*cosphi
    else
      b = sqrt(ONE - v0*v0)
      uvw(1) = mu*u0 + a*(u0*v0*cosphi + w0*sinphi)/b
      uvw(2) = mu*v0 - a*b*cosphi
      uvw(3) = mu*w0 + a*(v0*w0*cosphi - u0*sinphi)/b
    end if

  end function rotate_angle

end module physics
