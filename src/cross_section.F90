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
    real(8) :: macroxs_t
    real(8) :: microxs_a
    real(8) :: microxs_s
    type(Nuclide), pointer :: nuc => null()

    ! Get particle energy
    E = p % get_energy()

    ! Reset macros
    macroxs_t = ZERO

    ! Loop around nuclides
    do i = 1, n_nuclides
      nuc => nuclides(i)

      ! Get nuclide atom density
      dens = nuc % get_density()

      ! Get micro absorption cross section
      microxs_a = nuc % interp_xs_a(E) 

      ! Get micro scattering cross section
      microxs_s = nuc % interp_xs_s(E)

      ! Add macros
      macroxs_a = dens*microxs_a
      macroxs_s = dens*microxs_s

      ! Add macros to particle
      call p % add_macros(i, macroxs_a, macroxs_s)

      ! Accumulate total macroscopic cross section
      macroxs_t = macroxs_t + macroxs_a + macroxs_s
    end do

    ! Set macro total
    call p % set_macro_total(macroxs_t)

  end subroutine calculate_xs

end module cross_section
