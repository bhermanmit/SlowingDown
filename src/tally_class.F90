module tally_class

  implicit none
  private

  type, public :: Tally
    integer :: n_bins
    real(8), allocatable :: bins(:)
    real(8), allocatable :: s1(:)
    real(8), allocatable :: s2(:)
    contains
      procedure, public :: clear => tally_clear
  end type Tally

contains

!===============================================================================
! TALLY_CLEAR
!===============================================================================

  subroutine tally_clear(self)

    class(Tally), intent(inout) :: self

    if (allocated(self % bins)) deallocate(self % bins)
    if (allocated(self % s1)) deallocate(self % s1)
    if (allocated(self % s2)) deallocate(self % s2)

  end subroutine tally_clear

end module tally_class
