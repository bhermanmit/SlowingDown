module tally_class

  use constants,  only: MAX_WORD_LEN, MAX_ENERGY, MIN_ENERGY, EQUAL_LETHARGY, &
                        ZERO
  use math,       only: binary_search

  implicit none
  private

  type, public :: Tally
    character(len=MAX_WORD_LEN) :: score
    integer :: nbins
    integer :: type
    real(8), allocatable :: bins(:)
    real(8), allocatable :: flux_s1(:)
    real(8), allocatable :: flux_s2(:)
    contains
      procedure, public :: add_flux_score => tally_add_flux_score
      procedure, public :: clear => tally_clear
      procedure, public :: get_energy_bin => tally_get_energy_bin
      procedure, public :: get_nbins => tally_get_nbins
      procedure, public :: initialize => tally_initialize
      procedure, public :: set_nbins => tally_set_nbins
      procedure, public :: set_type => tally_set_type
      procedure, public :: write => tally_write
  end type Tally

  type(Tally), public, save :: tal

contains

!===============================================================================
! TALLY_ADD_FLUX_SCORE
!===============================================================================

  subroutine tally_add_flux_score(self, energy, score)

    class(Tally) :: self
    integer :: group
    real(8) :: energy
    real(8) :: score

    ! Get energy group
    group = self % get_energy_bin(energy)

    ! Save tally
    self % flux_s1(group) = self % flux_s1(group) + score
    self % flux_s2(group) = self % flux_s2(group) + score**2

  end subroutine tally_add_flux_score

!===============================================================================
! TALLY_CLEAR
!===============================================================================

  subroutine tally_clear(self)

    class(Tally), intent(inout) :: self

    if (allocated(self % bins)) deallocate(self % bins)
    if (allocated(self % flux_s1)) deallocate(self % flux_s1)
    if (allocated(self % flux_s2)) deallocate(self % flux_s2)

  end subroutine tally_clear

!===============================================================================
! TALLY_GET_ENERGY_BIN
!===============================================================================

  function tally_get_energy_bin(self, energy) result(energy_bin)

    class(Tally) :: self
    real(8) :: energy
    integer :: energy_bin
    
    ! Get energy group
    energy_bin = binary_search(self % bins, self % nbins+1, energy)

  end function tally_get_energy_bin

!===============================================================================
! TALLY_GET_NBINS
!===============================================================================

  function tally_get_nbins(self) result(nbins)

    class(Tally) :: self
    integer :: nbins
    
    nbins = self % nbins

  end function tally_get_nbins

!===============================================================================
! TALLY_INITIALIZE
!===============================================================================

  subroutine tally_initialize(self)

    class(Tally), intent(inout) :: self

    integer :: i
    integer :: nbins
    real(8) :: u_low
    real(8) :: u_hi
    real(8) :: u_width

    ! check for equal lethargy bins
    if (self % type == EQUAL_LETHARGY) then

      ! get number of lethargy bins
      nbins = self % nbins

      ! get lethargy width
      u_low = log(MAX_ENERGY/MAX_ENERGY)
      u_hi = log(MAX_ENERGY/MIN_ENERGY)
      u_width = (u_hi - u_low) / real(nbins, 8)
 
    end if

    ! Adjust bins
    allocate(self % bins(nbins + 1))

    ! Set bins
    self % bins(1) = 0.0
    do i = 1, nbins-1
      self % bins(i+1) = MAX_ENERGY/exp(u_hi - real(i,8)*u_width)
    end do
    self % bins(nbins+1) = MAX_ENERGY

    ! allocate scoring bins
    allocate(self % flux_s1(nbins))
    allocate(self % flux_s2(nbins))
    self % flux_s1 = ZERO
    self % flux_s2 = ZERO

  end subroutine tally_initialize

!===============================================================================
! TALLY_SET_NBINS
!===============================================================================

  subroutine tally_set_nbins(self, nbins)

    class(Tally) :: self
    integer :: nbins

    self % nbins = nbins

  end subroutine tally_set_nbins

!===============================================================================
! TALLY_SET_TYPE
!===============================================================================

  subroutine tally_set_type(self, type)

    class(Tally), intent(inout) :: self
    integer, intent(in) :: type

    self % type = type

  end subroutine tally_set_type

!===============================================================================
! TALLY_WRITE
!===============================================================================

  subroutine tally_write(self)

    class(Tally), intent(inout) :: self

    integer :: i

    ! Write out flux
    open(FILE="flux.out", UNIT=100, ACTION="write")
    do i = 1, self % nbins
      write(100, *) (self % bins(i) + self % bins(i+1))/2.0, self % flux_s1(i)
    end do
    close(100)

  end subroutine tally_write

end module tally_class
