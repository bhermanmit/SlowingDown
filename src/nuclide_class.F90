module nuclide_class

  use constants,  only: MAX_WORD_LEN, MAX_LINE_LEN

  implicit none

  type, public :: Nuclide
    private

    ! Nuclide properties
    character(len=MAX_WORD_LEN) :: name ! name of nuclide
    character(len=MAX_LINE_LEN) :: xs_s_path ! path to scattering cross section
    character(len=MAX_LINE_LEN) :: xs_a_path ! path to absorption cross sectoin
    real(8), allocatable :: E(:) ! energy array of nuclide cross section
    real(8), allocatable :: xs_s(:) ! scattering micro xs
    real(8), allocatable :: xs_a(:) ! absorption mircro xs

    ! Nuclide methods
    contains
      procedure, public :: clear => nuclide_clear
      procedure, public :: get_xs_s_path => nuclide_get_xs_s_path
      procedure, public :: set_name => nuclide_set_name
      procedure, public :: set_xs_s_path => nuclide_set_xs_s_path

  end type Nuclide

  ! Saved module variables
  integer, save :: n_nuclides
  type(Nuclide), allocatable, save, target :: nuclides(:)

contains

!===============================================================================
! NUCLIDE_CLEAR
!===============================================================================

  subroutine nuclide_clear(self)

    class(Nuclide), intent(inout) :: self

    if(allocated(self % E)) deallocate(self % E)
    if(allocated(self % xs_s)) deallocate(self % xs_s)
    if(allocated(self % xs_a)) deallocate(self % xs_a)

  end subroutine nuclide_clear

!===============================================================================
! NUCLIDE_GET_XS_S_PATH
!===============================================================================

  function nuclide_get_xs_s_path(self) result(xs_s_path)

    class(Nuclide) :: self
    character(len=MAX_LINE_LEN) :: xs_s_path

    xs_s_path = self % xs_s_path

  end function nuclide_get_xs_s_path

!===============================================================================
! NUCLIDE_SET_NAME
!===============================================================================

  subroutine nuclide_set_name(self, name)

    class(Nuclide), intent(inout) :: self
    character(len=MAX_WORD_LEN), intent(in) :: name

    self % name = name 

  end subroutine nuclide_set_name

!===============================================================================
! NUCLIDE_SET_XS_S_PATH
!===============================================================================

  subroutine nuclide_set_xs_s_path(self, xs_s_path)

    class(Nuclide), intent(inout) :: self
    character(len=MAX_LINE_LEN), intent(in) :: xs_s_path

    self % xs_s_path = xs_s_path

  end subroutine nuclide_set_xs_s_path

end module nuclide_class
