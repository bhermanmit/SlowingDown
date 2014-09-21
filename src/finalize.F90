module finalize

  use nuclide_class,  only: nuclides, n_nuclides

  implicit none

contains

!===============================================================================
! FINALIZE_RUN
!===============================================================================

  subroutine finalize_run()

    integer :: i

    ! Free nuclide memory
    do i = 1, n_nuclides
      call nuclides(i) % clear()
    end do

  end subroutine finalize_run

end module finalize
