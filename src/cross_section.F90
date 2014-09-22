module cross_section

  use constants,       only: ZERO
  use nuclide_class,   only: n_nuclides, nuclides, Nuclide
  use particle_class,  only: Particle

  implicit none

contains

!===============================================================================
! CALCULATE_XS
!===============================================================================

  subroutine calculate_xs(p)

    type(Particle) :: p

    integer :: i
    real(8) :: dens
    real(8) :: E
    real(8) :: macroxs_a
    real(8) :: macroxs_s
    real(8) :: microxs_a
    real(8) :: microxs_s
    type(Nuclide), pointer :: nuc => null()

    ! Get particle energy
    E = p % get_energy()

    ! Reset macros
    macroxs_a = ZERO
    macroxs_s = ZERO

    ! Loop around nuclides
    do i = 1, n_nuclides
      nuc => nuclides(i)

      ! Get nuclide atom density
      dens = nuc % get_density()

      ! Get micro absorption cross section
      microxs_a = nuc % interp_xs_a(E) 

      ! Get micro scattering cross section
      microxs_s = nuc % interp_xs_s(E)

      ! Add micros to particle
      call p % add_micros(i, microxs_a, microxs_s)

      ! Add macros
      macroxs_a = macroxs_a + dens*microxs_a
      macroxs_s = macroxs_s + dens*microxs_s
    end do

    ! Set macros in particle
    call p % add_macros(macroxs_a, macroxs_s)

  end subroutine calculate_xs

end module cross_section
