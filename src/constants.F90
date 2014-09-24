module constants

  implicit none

  ! SlowingDown major, minor, and release numbers
  integer, parameter :: VERSION_MAJOR   = 0
  integer, parameter :: VERSION_MINOR   = 0
  integer, parameter :: VERSION_RELEASE = 0

  ! Maximum number of words in a single line, length of line, and length of
  ! single word
  integer, parameter :: MAX_WORDS    = 500
  integer, parameter :: MAX_LINE_LEN = 250
  integer, parameter :: MAX_WORD_LEN = 150
  integer, parameter :: MAX_FILE_LEN = 255

  ! ============================================================================
  ! PHYSICAL CONSTANTS

  ! Values here are from the Committee on Data for Science and Technology
  ! (CODATA) 2010 recommendation (doi:10.1103/RevModPhys.84.1527).

  real(8), parameter ::            &
       PI           = 3.1415926535898_8, & ! pi
       MASS_NEUTRON = 1.008664916,       & ! mass of a neutron in amu
       MASS_PROTON  = 1.007276466812,    & ! mass of a proton in amu
       AMU          = 1.660538921e-27,   & ! 1 amu in kg
       N_AVOGADRO   = 0.602214129,       & ! Avogadro's number in 10^24/mol
       K_BOLTZMANN  = 8.6173324e-11,     & ! Boltzmann constant in MeV/K
       INFINITY     = huge(0.0_8),       & ! positive infinity
       ZERO         = 0.0_8,             &
       ONE          = 1.0_8,             &
       TWO          = 2.0_8

  ! Maximum number of allowed collisions
  integer, parameter :: MAX_COLLISIONS = 1

  ! Reaction that occurred
  integer, parameter :: REACTION_ABSORBED = 1
  integer, parameter :: REACTION_SCATTERED = 2

end module constants
