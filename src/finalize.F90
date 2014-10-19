module finalize

  use nuclide_class,  only: nuclides, n_nuclides
  use tally_class,    only: flux_tal, abs_tal, scat_tal, r2c_tal, &
                            outscatc_tal, winscatc_tal

  implicit none

contains

!===============================================================================
! FINALIZE_RUN
!===============================================================================

  subroutine finalize_run()

    integer :: i

    ! Write tallies
    call flux_tal % write("flux.out")
    call abs_tal % write("absorption.out")
    call scat_tal % write("scattering.out")
    call r2c_tal % write("r2c.out")
    call outscatc_tal % write("outscatterc.out")
    call winscatc_tal % write("withinscatterc.out")

    ! Free nuclide memory
    do i = 1, n_nuclides
      call nuclides(i) % clear()
    end do

    ! Free tally memory
    call flux_tal % clear() 
    call abs_tal % clear()
    call scat_tal % clear() 
    call r2c_tal % clear()
    call outscatc_tal % clear()
    call winscatc_tal % clear()

  end subroutine finalize_run

end module finalize
