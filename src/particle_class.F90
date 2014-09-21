module particle_class 

  use constants,  only: ONE, TWO, PI
  use random,     only: prn

  implicit none

  ! Particle class
  type, public :: Particle
    private

    ! Particle properties
    real(8) :: xyz(3) ! spatial position
    real(8) :: uvw(3) ! direction of travel
    real(8) :: E ! energy

    ! Particle methods
    contains
      procedure, public :: start => particle_start

  end type Particle

contains

!======================================================================
! PARTICLE_START samples particles' starting properties
!======================================================================

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