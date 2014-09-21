module input_xml

  use output,        only: write_message, fatal_error, MAX_LINE_LEN, message
  use xml_interface

  implicit none

contains

!===============================================================================
! READ_INPUT_XML calls each of the separate subroutines for reading settings,
! geometry, materials, and tallies.
!===============================================================================

  subroutine read_input_xml()

    character(MAX_LINE_LEN) :: filename
    logical :: file_exists
    type(Node), pointer :: doc => null()

    ! Display output message
    message = "Reading input XML file..."
    call write_message()

    ! Check if settings.xml exists
    filename = "input.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      message = "Input file does not exist."
      call fatal_error()
    end if

    ! Parse settings.xml file
    call open_xmldoc(doc, filename)

    ! Close settings XML file
    call close_xmldoc(doc)

  end subroutine read_input_xml

end module input_xml
