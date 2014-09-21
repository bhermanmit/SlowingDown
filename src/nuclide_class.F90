module nuclide_class

  use constants,  only: MAX_WORD_LEN

  implicit none

  type, public :: Nuclide
    private

    ! Nuclide properties
    character(len=MAX_WORD_LEN) :: name ! name of nuclide
    real(8), allocatable :: E(:) ! energy array of nuclide cross section
    real(8), allocatable :: xs_s(:) ! scattering micro xs
    real(8), allocatable :: xs_a(:) ! absorption mircro xs

    ! Nuclide methods
    contains
      procedure, public :: clear => nuclide_clear

  end type Nuclide

contains

!===============================================================================
! NUCLIDE_CLEAR
!===============================================================================

  subroutine nuclide_clear(self)

    class(Nuclide) :: self

    if(allocated(self % E)) deallocate(self % E)
    if(allocated(self % xs_s)) deallocate(self % xs_s)
    if(allocated(self % xs_a)) deallocate(self % xs_a)

  end subroutine nuclide_clear

end module nuclide_class
