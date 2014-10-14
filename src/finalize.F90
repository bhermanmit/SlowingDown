module finalize

  use nuclide_class,  only: nuclides, n_nuclides
  use tally_class, only: tal

  implicit none

contains

!===============================================================================
! FINALIZE_RUN
!===============================================================================

  subroutine finalize_run()

    integer :: i

    ! Write tallies
    call tal % write()

    ! Free nuclide memory
    do i = 1, n_nuclides
      call nuclides(i) % clear()
    end do

    ! Free tally memory
    call tal % clear()

  end subroutine finalize_run

end module finalize
