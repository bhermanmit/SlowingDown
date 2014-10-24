module tally_class

  use constants,  only: MAX_WORD_LEN, MAX_ENERGY, MIN_ENERGY, EQUAL_LETHARGY, &
                        ZERO, MAX_FILE_LEN
  use math,       only: binary_search

  implicit none
  private

  type, public :: Tally
    character(len=MAX_WORD_LEN) :: score
    integer :: nbins
    integer :: type
    real(8), allocatable :: bins(:)
    real(8), allocatable :: s1(:)
    real(8), allocatable :: s2(:)
    contains
      procedure, public :: add_score => tally_add_score
      procedure, public :: clear => tally_clear
      procedure, public :: get_energy_bin => tally_get_energy_bin
      procedure, public :: get_nbins => tally_get_nbins
      procedure, public :: initialize => tally_initialize
      procedure, public :: set_nbins => tally_set_nbins
      procedure, public :: set_type => tally_set_type
      procedure, public :: write => tally_write
      procedure, public :: write_energy => tally_write_energy
  end type Tally

  type(Tally), public, save :: flux_tal
  type(Tally), public, save :: abs_tal
  type(Tally), public, save :: scat_tal
  type(Tally), public, save :: r2c_tal
  type(Tally), public, save :: outscatc_tal
  type(Tally), public, save :: winscatc_tal
  type(Tally), public, save :: wc_tal
  type(Tally), public, save :: p1_scat_tal

contains

!===============================================================================
! TALLY_ADD_SCORE
!===============================================================================

  subroutine tally_add_score(self, group, score)

    class(Tally), intent(inout) :: self
    integer, intent(in) :: group
    real(8), intent(in) :: score

    integer :: bin

    ! Change group to bin
    bin = self % get_nbins() - group + 1

    ! Save tally
    self % s1(bin) = self % s1(bin) + score
    self % s2(bin) = self % s2(bin) + score**2

  end subroutine tally_add_score

!===============================================================================
! TALLY_CLEAR
!===============================================================================

  subroutine tally_clear(self)

    class(Tally), intent(inout) :: self

    if (allocated(self % bins)) deallocate(self % bins)
    if (allocated(self % s1)) deallocate(self % s1)
    if (allocated(self % s2)) deallocate(self % s2)

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
    allocate(self % s1(nbins))
    allocate(self % s2(nbins))
    self % s1 = ZERO
    self % s2 = ZERO

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
! TALLY_WRITE_ENERGY
!===============================================================================

  subroutine tally_write_energy(self)

    class(Tally), intent(inout) :: self

    integer :: i

    ! Write out energy pts
    open(FILE="energy.out", UNIT=99, ACTION="write")
    do i = 1, self % nbins
      write(99, *) (self % bins(i) + self % bins(i+1))/2.0
    end do
    close(99)

  end subroutine tally_write_energy

!===============================================================================
! TALLY_WRITE
!===============================================================================

  subroutine tally_write(self, filename)

    class(Tally), intent(inout) :: self
    character(len=*) :: filename

    integer :: i
    integer :: un

    ! Write out energy pts
    open(FILE="energy.out", UNIT=99, ACTION="write")
    do i = 1, self % nbins
      write(99, *) (self % bins(i) + self % bins(i+1))/2.0
    end do
    close(99)

    ! Write out flux
    open(FILE=trim(filename), NEWUNIT=un, ACTION="write")
    do i = 1, self % nbins
      write(un, *) self % s1(i)
    end do
    close(un)

  end subroutine tally_write

end module tally_class
