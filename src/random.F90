module random

  implicit none

contains

!===============================================================================
! PRN
!===============================================================================

  function prn() result(pseudo_rn)

    real(8) :: pseudo_rn

    call random_number(pseudo_rn)

  end function prn

!===============================================================================
! INITIALIZE_PRN
!===============================================================================

  subroutine initialize_prn()

    call random_seed()

  end subroutine initialize_prn

end module random
