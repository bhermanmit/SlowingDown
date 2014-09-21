module nuclide_class

  use constants,  only: MAX_WORD_LEN, MAX_LINE_LEN

  implicit none
  private

  type, public :: Nuclide
    private

    ! Nuclide properties
    character(len=MAX_WORD_LEN) :: name ! name of nuclide
    character(len=MAX_LINE_LEN) :: xs_file ! cross section file
    real(8), allocatable :: energy(:) ! energy array of nuclide cross section
    real(8), allocatable :: xs_s(:) ! scattering micro xs
    real(8), allocatable :: xs_a(:) ! absorption mircro xs

    ! Nuclide methods
    contains
      procedure, public :: clear => nuclide_clear
      procedure, public :: get_xs_file => nuclide_get_xs_file
      procedure, public :: set_energy => nuclide_set_energy
      procedure, public :: set_name => nuclide_set_name
      procedure, public :: set_xs_s => nuclide_set_xs_s
      procedure, public :: set_xs_file => nuclide_set_xs_file
      procedure, public :: write => nuclide_write

  end type Nuclide

  ! Saved module variables
  integer, public, save :: n_nuclides
  type(Nuclide), public, allocatable, save, target :: nuclides(:)

contains

!===============================================================================
! NUCLIDE_CLEAR
!===============================================================================

  subroutine nuclide_clear(self)

    class(Nuclide), intent(inout) :: self

    if(allocated(self % energy)) deallocate(self % energy)
    if(allocated(self % xs_s)) deallocate(self % xs_s)
    if(allocated(self % xs_a)) deallocate(self % xs_a)

  end subroutine nuclide_clear

!===============================================================================
! NUCLIDE_GET_XS_FILE
!===============================================================================

  function nuclide_get_xs_file(self) result(xs_file)

    class(Nuclide) :: self
    character(len=MAX_LINE_LEN) :: xs_file

    xs_file = self % xs_file

  end function nuclide_get_xs_file

!===============================================================================
! NUCLIDE_SET_ENERGY
!===============================================================================

  subroutine nuclide_set_energy(self, energy)

    class(Nuclide), intent(inout) :: self
    real(8), intent(in) :: energy(:)

    allocate(self % energy(size(energy)))
    self % energy = energy

  end subroutine nuclide_set_energy

!===============================================================================
! NUCLIDE_SET_NAME
!===============================================================================

  subroutine nuclide_set_name(self, name)

    class(Nuclide), intent(inout) :: self
    character(len=MAX_WORD_LEN), intent(in) :: name

    self % name = name 

  end subroutine nuclide_set_name

!===============================================================================
! NUCLIDE_SET_XS_FILE
!===============================================================================

  subroutine nuclide_set_xs_file(self, xs_file)

    class(Nuclide), intent(inout) :: self
    character(len=MAX_LINE_LEN), intent(in) :: xs_file

    self % xs_file = xs_file

  end subroutine nuclide_set_xs_file

!===============================================================================
! NUCLIDE_SET_XS_S
!===============================================================================

  subroutine nuclide_set_xs_s(self, xs_s)

    class(Nuclide), intent(inout) :: self
    real(8), intent(in) :: xs_s(:)

    allocate(self % xs_s(size(xs_s)))
    self % xs_s = xs_s

  end subroutine nuclide_set_xs_s

!===============================================================================
! NUCLIDE_WRITE
!===============================================================================

  subroutine nuclide_write(self)

    class(Nuclide), intent(inout) :: self

    print *, self % name
    print *, self % energy
    print *, self % xs_s

  end subroutine nuclide_write

end module nuclide_class
