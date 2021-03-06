module particle_class 

  use constants,   only: ONE, TWO, PI, MAX_ENERGY, MIN_ENERGY, &
                         REACTION_SCATTERED
  use cross_section_header 
  use random,      only: prn
  use tally_class,    only: flux_tal, abs_tal, scat_tal, r2c_tal, &
                            outscatc_tal, winscatc_tal, wc_tal, &
                            p1_scat_tal

  implicit none

  ! Particle class
  type, public :: Particle
    private

    ! Particle properties
    integer :: energy_group ! energy group in tally structure
    integer :: energy_group_last ! energy group in tally structure
    integer :: n_collisions ! number of collisions
    integer :: nuclide_index ! Index of collision nuclide in macro
    integer :: reaction_type ! The reaction that occurred
    logical :: alive ! is the particle alive?
    real(8) :: dist ! distance particle traveled
    real(8) :: E ! energy
    real(8) :: E_last ! energy
    real(8) :: macro_total
    real(8) :: mulab ! mu in LAB system
    real(8) :: uvw(3) ! direction of travel
    real(8) :: weight ! statistical weight
    real(8) :: xyz(3) ! spatial position
    real(8) :: xyz_start(3) ! starting spatial position for tallies
    type(MacroXS), allocatable :: macro(:) ! Current macroscopic xs

    ! Particle methods
    contains
      procedure, public :: add_macros => particle_add_macros
      procedure, public :: begin_collision => particle_begin_collision
      procedure, public :: calc_displacement2 => particle_calc_displacement2
      procedure, public :: clear => particle_clear
      procedure, public :: get_alive => particle_get_alive
      procedure, public :: get_distance => particle_get_distance
      procedure, public :: get_energy => particle_get_energy
      procedure, public :: get_energy_last => particle_get_energy_last
      procedure, public :: get_energy_group => particle_get_energy_group
      procedure, public :: get_energy_group_last => &
                           particle_get_energy_group_last
      procedure, public :: get_macro_total => particle_get_macro_total
      procedure, public :: get_mulab => particle_get_mulab
      procedure, public :: get_nuclide_macroxs_a => &
                           particle_get_nuclide_macroxs_a
      procedure, public :: get_nuclide_macroxs_t => &
                           particle_get_nuclide_macroxs_t
      procedure, public :: get_nuclide_index => particle_get_nuclide_index
      procedure, public :: get_reaction_type => particle_get_reaction_type
      procedure, public :: get_uvw => particle_get_uvw
      procedure, public :: get_weight => particle_get_weight
      procedure, public :: initialize => particle_initialize
      procedure, public :: lookup_energy_group => particle_lookup_energy_group
      procedure, public :: mark_position => particle_mark_position
      procedure, public :: move => particle_move
      procedure, public :: set_alive => particle_set_alive
      procedure, public :: set_distance => particle_set_distance
      procedure, public :: set_energy => particle_set_energy
      procedure, public :: set_energy_group => particle_set_energy_group
      procedure, public :: set_macro_total => particle_set_macro_total
      procedure, public :: set_mulab => particle_set_mulab
      procedure, public :: set_nuclide_index => particle_set_nuclide_index
      procedure, public :: set_n_collisions => particle_set_n_collisions
      procedure, public :: set_reaction_type => particle_set_reaction_type
      procedure, public :: set_uvw => particle_set_uvw
      procedure, public :: set_weight => particle_set_weight
      procedure, public :: start => particle_start
      procedure, public :: tally => particle_tally

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
! PARTICLE_BEGIN_COLLISION
!===============================================================================

  subroutine particle_begin_collision(self)

    class(Particle), intent(inout) :: self

    self % E_last = self % E
    self % energy_group_last = self % energy_group

  end subroutine particle_begin_collision

!===============================================================================
! PARTICLE_CALC_DISPLACEMENT
!===============================================================================

  function particle_calc_displacement2(self) result(displacement2)

    class(Particle) :: self
    real(8) :: displacement2

    displacement2 = sum((self % xyz - self % xyz_start)**2)

  end function particle_calc_displacement2

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
! PARTICLE_GET_ENERGY_LAST
!===============================================================================

  function particle_get_energy_last(self) result(energy)

    class(Particle) :: self
    real(8) :: energy

    energy = self % E_last

  end function particle_get_energy_last

!===============================================================================
! PARTICLE_GET_ENERGY_GROUP
!===============================================================================

  function particle_get_energy_group(self) result(energy_group)

    class(Particle) :: self
    integer :: energy_group

    energy_group = self % energy_group

  end function particle_get_energy_group

!===============================================================================
! PARTICLE_GET_ENERGY_GROUP_LAST
!===============================================================================

  function particle_get_energy_group_last(self) result(energy_group)

    class(Particle) :: self
    integer :: energy_group

    energy_group = self % energy_group_last

  end function particle_get_energy_group_last

!===============================================================================
! PARTICLE_GET_MACRO_TOTAL
!===============================================================================

  function particle_get_macro_total(self) result(macro_total)

    class(Particle) :: self
    real(8) :: macro_total

    macro_total = self % macro_total

  end function particle_get_macro_total

!===============================================================================
! PARTICLE_GET_MULAB
!===============================================================================

  function particle_get_mulab(self) result(mulab)

    class(Particle) :: self
    real(8) :: mulab

    mulab = self % mulab

  end function particle_get_mulab

!===============================================================================
! PARTICLE_GET_NUCLIDE_MACROXS_A
!===============================================================================

  function particle_get_nuclide_macroxs_a(self, i) result(macroxs_a)

    class(Particle) :: self
    integer :: i
    real(8) :: macroxs_a

    macroxs_a = self % macro(i) % xs_a

  end function particle_get_nuclide_macroxs_a

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
! PARTICLE_GET_REACTION_TYPE
!===============================================================================

  function particle_get_reaction_type(self) result(reaction_type)

    class(Particle) :: self
    integer :: reaction_type

    reaction_type = self % reaction_type

  end function particle_get_reaction_type

!===============================================================================
! PARTICLE_GET_UVW
!===============================================================================

  function particle_get_uvw(self) result(uvw)

    class(Particle) :: self
    real(8) :: uvw(3)

    uvw = self % uvw

  end function particle_get_uvw

!===============================================================================
! PARTICLE_GET_WEIGHT
!===============================================================================

  function particle_get_weight(self) result(weight)

    class(Particle) :: self
    real(8) :: weight

    weight = self % weight

  end function particle_get_weight

!===============================================================================
! PARTICLE_INITIALIZE
!===============================================================================

  subroutine particle_initialize(self, n)

    class(Particle), intent(inout) :: self
    integer, intent(in) :: n

    allocate(self % macro(n))

  end subroutine particle_initialize

!===============================================================================
! PARTICLE_LOOKUP_ENERGY_GROUP
!===============================================================================

  function particle_lookup_energy_group(self) result(energy_group)

    class(Particle), intent(inout) :: self
    integer :: energy_group

    ! get the energy bin
    energy_group = flux_tal % get_energy_bin(self % E)

    ! convert energy bin to energy group
    energy_group = flux_tal % get_nbins() - energy_group + 1

  end function particle_lookup_energy_group

!===============================================================================
! PARTICLE_MARK_POSITION
!===============================================================================

  subroutine particle_mark_position(self)

    class(Particle), intent(inout) :: self

    self % xyz_start = self % xyz

  end subroutine particle_mark_position

!===============================================================================
! PARTICLE_MOVE
!===============================================================================

  subroutine particle_move(self)

    class(Particle), intent(inout) :: self

    self % xyz = self % xyz + self % dist * self % uvw

  end subroutine particle_move

!===============================================================================
! PARTICLE_SET_ALIVE
!===============================================================================

  subroutine particle_set_alive(self, alive)

    class(Particle), intent(inout) :: self
    logical, intent(in) :: alive

    self % alive = alive

  end subroutine particle_set_alive

!===============================================================================
! PARTICLE_SET_DISTANCE
!===============================================================================

  subroutine particle_set_distance(self, dist)

    class(Particle), intent(inout) :: self
    real(8), intent(in) :: dist

    self % dist = dist

  end subroutine particle_set_distance

!===============================================================================
! PARTICLE_SET_ENERGY
!===============================================================================

  subroutine particle_set_energy(self, energy)

    class(Particle), intent(inout) :: self
    real(8), intent(in) :: energy

    ! set the energy
    self % E = energy 

    ! set the energy group
    call self % set_energy_group(self % lookup_energy_group())

  end subroutine particle_set_energy

!===============================================================================
! PARTICLE_SET_ENERGY_GROUP
!===============================================================================

  subroutine particle_set_energy_group(self, energy_group)

    class(Particle), intent(inout) :: self
    integer, intent(in) :: energy_group

    self % energy_group = energy_group

  end subroutine particle_set_energy_group

!===============================================================================
! PARTICLE_SET_MACRO_TOTAL
!===============================================================================

  subroutine particle_set_macro_total(self, macro_total)

    class(Particle), intent(inout) :: self
    real(8), intent(in) :: macro_total

    self % macro_total = macro_total

  end subroutine particle_set_macro_total

!===============================================================================
! PARTICLE_SET_MULAB
!===============================================================================

  subroutine particle_set_mulab(self, mulab)

    class(Particle), intent(inout) :: self
    real(8), intent(in) :: mulab

    self % mulab = mulab

  end subroutine particle_set_mulab

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
! PARTICLE_SET_REACTION_TYPE
!===============================================================================

  subroutine particle_set_reaction_type(self, reaction_type)

    class(Particle), intent(inout) :: self
    integer, intent(in) :: reaction_type

    self % reaction_type = reaction_type

  end subroutine particle_set_reaction_type

!===============================================================================
! PARTICLE_SET_UVW
!===============================================================================

  subroutine particle_set_uvw(self, uvw)

    class(Particle), intent(inout) :: self
    real(8), intent(in) :: uvw(3)

    self % uvw = uvw

  end subroutine particle_set_uvw

!===============================================================================
! PARTICLE_SET_WEIGHT
!===============================================================================

  subroutine particle_set_weight(self, weight)

    class(Particle), intent(inout) :: self
    real(8), intent(in) :: weight

    self % weight = weight

  end subroutine particle_set_weight

!===============================================================================
! PARTICLE_START samples particles' starting properties
!===============================================================================

  subroutine particle_start(self)

    class(Particle) :: self

    real(8) :: cosz ! cosine of z
    real(8) :: mu   ! polar angle
    real(8) :: phi  ! azimuthal angle
    real(8) :: sinz ! sine of z

    ! Starting weight
    call self % set_weight(ONE)

    ! Starting location
    self % xyz = (/0.0, 0.0, 0.0/)
    call self % mark_position() ! sets this as initial reference position for
                                ! displacement calculation for migration

    ! Starting direction
    phi = TWO*PI*prn()
    mu = TWO*prn() - ONE
    self % uvw(1) = mu
    self % uvw(2) = sqrt(ONE - mu**2) * cos(phi)
    self % uvw(3) = sqrt(ONE - mu**2) * sin(phi)

    ! Starting energy
    self % E = watt_spectrum()

    ! Make sure it is under 20.0 MeV and above 1.e-11 MeV
    self % E = max(self % E, MIN_ENERGY)
    self % E = min(self % E, MAX_ENERGY)

    ! Get neutrons energy group in tally structure
    self % energy_group = self % lookup_energy_group()

    ! Particle is alive
    self % alive = .true.

  end subroutine particle_start

!===============================================================================
! PARTICLE_TALLY
!===============================================================================

  subroutine particle_tally(self)

    class(Particle), intent(inout) :: self

    integer :: energy_group
    integer :: energy_group_last
    integer :: i
    real(8) :: displacement2
    real(8) :: weight

    ! Extract data
    weight = self % get_weight()
    energy_group = self % get_energy_group()
    energy_group_last = self % get_energy_group_last()

    ! Record flux tally at each collision, uses pre-collision energy group
    call flux_tal % add_score(energy_group_last, weight/self % get_macro_total())

    ! Record scattering tallies
    if (self % get_reaction_type() == REACTION_SCATTERED) then
      call scat_tal % add_score(energy_group_last, weight)
      call p1_scat_tal % add_score(energy_group_last, weight*self % get_mulab())
      do i = energy_group, flux_tal % get_nbins()
        call winscatc_tal % add_score(i, weight)
      end do
    else
      call abs_tal % add_score(energy_group_last, weight)
    end if

    ! Record cumulative tallies
    displacement2 = self % calc_displacement2()
    if (self % get_reaction_type() == REACTION_SCATTERED) then

      ! for a scattering out of group h to group g, we
      ! tally the crow flight distance from group g all the
      ! way to the energy group just above. The neutron position
      ! is then marked as the reference position for the next
      ! collision 
      do i = energy_group_last, energy_group - 1
        call r2c_tal % add_score(i, weight*displacement2)
        call wc_tal % add_score(i, weight)
        call outscatc_tal % add_score(i, weight)
      end do
    else

      ! for absorption, we tally in the absorbed groups and all the
      ! groups below in energy
      do i = energy_group, flux_tal % get_nbins()
        call r2c_tal % add_score(i, weight*displacement2)
        call wc_tal % add_score(i, weight)
      end do
    end if

  end subroutine particle_tally

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
