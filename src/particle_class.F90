module particle_class 

  use constants,  only: ONE, TWO, PI
  use cross_section_header 
  use random,     only: prn

  implicit none

  ! Particle class
  type, public :: Particle
    private

    ! Particle properties
    integer :: n_collisions ! number of collisions
    integer :: nuclide_index ! Index of collision nuclide in macro
    logical :: alive ! is the particle alive?
    real(8) :: dist ! distance particle traveled
    real(8) :: E ! energy
    real(8) :: macro_total
    real(8) :: uvw(3) ! direction of travel
    real(8) :: xyz(3) ! spatial position
    type(MacroXS), allocatable :: macro(:) ! Current macroscopic xs

    ! Particle methods
    contains
      procedure, public :: add_macros => particle_add_macros
      procedure, public :: clear => particle_clear
      procedure, public :: get_alive => particle_get_alive
      procedure, public :: get_distance => particle_get_distance
      procedure, public :: get_energy => particle_get_energy
      procedure, public :: get_macro_total => particle_get_macro_total
      procedure, public :: get_nuclide_macroxs_t => &
                           particle_get_nuclide_macroxs_t
      procedure, public :: get_nuclide_index => particle_get_nuclide_index
      procedure, public :: initialize => particle_initialize
      procedure, public :: set_distance => particle_set_distance
      procedure, public :: set_macro_total => particle_set_macro_total
      procedure, public :: set_nuclide_index => particle_set_nuclide_index
      procedure, public :: set_n_collisions => particle_set_n_collisions
      procedure, public :: start => particle_start

  end type Particle

contains

!===============================================================================
! PARTICLE_ADD_MACROS
!===============================================================================

  subroutine particle_add_macros(self, i, macroxs_a, macroxs_s)

    class(Particle), intent(inout) :: self
    integer, intent(in) :: i
    real(8), intent(in) :: macroxs_a
    real(8), intent(in) :: macroxs_s

    self % macro(i) % xs_a = macroxs_a
    self % macro(i) % xs_s = macroxs_s
    self % macro(i) % xs_t = macroxs_a + macroxs_s

  end subroutine particle_add_macros 

!===============================================================================
! PARTICLE_CLEAR
!===============================================================================

  subroutine particle_clear(self)

    class(Particle), intent(inout) :: self

    if (allocated(self % macro)) deallocate(self % macro)

  end subroutine particle_clear

!===============================================================================
! PARTICLE_GET_ALIVE
!===============================================================================

  function particle_get_alive(self) result(alive)

    class(Particle) :: self
    logical :: alive

    alive = self % alive

  end function particle_get_alive

!===============================================================================
! PARTICLE_GET_DISTANCE
!===============================================================================

  function particle_get_distance(self) result(dist)

    class(Particle) :: self
    real(8) :: dist

    dist = self % dist

  end function particle_get_distance

!===============================================================================
! PARTICLE_GET_ENERGY
!===============================================================================

  function particle_get_energy(self) result(energy)

    class(Particle) :: self
    real(8) :: energy

    energy = self % E

  end function particle_get_energy

!===============================================================================
! PARTICLE_GET_MACRO_TOTAL
!===============================================================================

  function particle_get_macro_total(self) result(macro_total)

    class(Particle) :: self
    real(8) :: macro_total

    macro_total = self % macro_total

  end function particle_get_macro_total

!===============================================================================
! PARTICLE_GET_NUCLIDE_MACROXS_T
!===============================================================================

  function particle_get_nuclide_macroxs_t(self, i) result(macroxs_t)

    class(Particle) :: self
    integer :: i
    real(8) :: macroxs_t

    macroxs_t = self % macro(i) % xs_t

  end function particle_get_nuclide_macroxs_t

!===============================================================================
! PARTICLE_GET_NUCLIDE_INDEX
!===============================================================================

  function particle_get_nuclide_index(self) result(nuclide_index)

    class(Particle) :: self
    integer :: nuclide_index

    nuclide_index = self % nuclide_index

  end function particle_get_nuclide_index

!===============================================================================
! PARTICLE_INITIALIZE
!===============================================================================

  subroutine particle_initialize(self, n)

    class(Particle), intent(inout) :: self
    integer, intent(in) :: n

    allocate(self % macro(n))

  end subroutine particle_initialize

!===============================================================================
! PARTICLE_SET_DISTANCE
!===============================================================================

  subroutine particle_set_distance(self, dist)

    class(Particle), intent(inout) :: self
    real(8), intent(in) :: dist

    self % dist = dist

  end subroutine particle_set_distance

!===============================================================================
! PARTICLE_SET_MACRO_TOTAL
!===============================================================================

  subroutine particle_set_macro_total(self, macro_total)

    class(Particle), intent(inout) :: self
    real(8), intent(in) :: macro_total

    self % macro_total = macro_total

  end subroutine particle_set_macro_total

!===============================================================================
! PARTICLE_SET_NUCLIDE_INDEX
!===============================================================================

  subroutine particle_set_nuclide_index(self, nuclide_index)

    class(Particle), intent(inout) :: self
    integer, intent(in) :: nuclide_index

    self % nuclide_index = nuclide_index

  end subroutine particle_set_nuclide_index

!===============================================================================
! PARTICLE_SET_N_COLLISIONS
!===============================================================================

  subroutine particle_set_n_collisions(self, n_collisions)

    class(Particle), intent(inout) :: self
    integer, intent(in) :: n_collisions

    self % n_collisions = n_collisions

  end subroutine particle_set_n_collisions

!===============================================================================
! PARTICLE_START samples particles' starting properties
!===============================================================================

  subroutine particle_start(self)

    class(Particle) :: self

    real(8) :: cosz ! cosine of z
    real(8) :: phi  ! azimuthal angle
    real(8) :: sinz ! sine of z

    ! Starting location
    self % xyz = (/0.0, 0.0, 0.0/)

    ! Starting direction
    cosz = ONE - TWO*prn() ! z-cosine
    sinz = sqrt(ONE - cosz**2)
    phi = TWO*PI*prn()
    self % uvw(1) = sinz*cos(phi)
    self % uvw(2) = sinz*sin(phi)
    self % uvw(3) = cosz

    ! Starting energy
    self % E = watt_spectrum()

    ! Particle is alive
    self % alive = .true.

  end subroutine particle_start

!===============================================================================
! WATT_SPECTRUM samples the outgoing energy from a Watt energy-dependent fission
! spectrum. Although fitted parameters exist for many nuclides, generally the
! continuous tabular distributions (LAW 4) should be used in lieu of the Watt
! spectrum. This direct sampling scheme is an unpublished scheme based on the
! original Watt spectrum derivation (See F. Brown's MC lectures).
!===============================================================================

  function watt_spectrum() result(E_out)

    real(8)             :: E_out ! energy of emitted neutron

    real(8) :: a ! Watt parameter a
    real(8) :: b ! Watt parameter b
    real(8) :: w ! sampled from Maxwellian

    a     = 0.988
    b     = 2.249
    w     = maxwell_spectrum(a)
    E_out = w + a*a*b/4. + (TWO*prn() - ONE)*sqrt(a*a*b*w)

  end function watt_spectrum

!===============================================================================
! MAXWELL_SPECTRUM samples an energy from the Maxwell fission distribution based
! on a direct sampling scheme. The probability distribution function for a
! Maxwellian is given as p(x) = 2/(T*sqrt(pi))*sqrt(x/T)*exp(-x/T). This PDF can
! be sampled using rule C64 in the Monte Carlo Sampler LA-9721-MS.
!===============================================================================

  function maxwell_spectrum(T) result(E_out)

    real(8), intent(in)  :: T     ! tabulated function of incoming E
    real(8)              :: E_out ! sampled energy

    real(8) :: r1, r2, r3  ! random numbers
    real(8) :: c           ! cosine of pi/2*r3
    
    r1 = prn() 
    r2 = prn()
    r3 = prn()

    ! determine cosine of pi/2*r
    c = cos(PI/TWO*r3)

    ! determine outgoing energy
    E_out = -T*(log(r1) + log(r2)*c*c)

  end function maxwell_spectrum

end module particle_class
