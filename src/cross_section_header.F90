module cross_section_header

  implicit none

  type :: MicroXS
    real(8) :: xs_a
    real(8) :: xs_s
    real(8) :: xs_t
  end type MicroXS

  type :: MacroXS
    real(8) :: xs_a
    real(8) :: xs_s
    real(8) :: xs_t
  end type MacroXS

end module cross_section_header
