module input_xml

  use constants,      only: MAX_WORD_LEN, MAX_LINE_LEN, EQUAL_LETHARGY
  use global
  use nuclide_class,  only: n_nuclides, nuclides, Nuclide
  use output,         only: write_message, fatal_error, message
  use tally_class,    only: flux_tal, abs_tal, scat_tal, r2c_tal, &
                            outscatc_tal, winscatc_tal, wc_tal, &
                            p1_scat_tal
  use xml_interface

  implicit none
  private
  public :: read_input_xml

contains

!===============================================================================
! READ_INPUT_XML calls each of the separate subroutines for reading settings,
! geometry, materials, and tallies.
!===============================================================================

  subroutine read_input_xml()

    character(len=MAX_LINE_LEN) :: filename
    character(len=MAX_WORD_LEN) :: tally_type
    character(len=MAX_LINE_LEN) :: xs_file
    integer :: i
    integer :: tally_nbins
    logical :: file_exists
    real(8) :: density
    type(Node), pointer :: doc => null()
    type(Node), pointer :: node_nuc => null()
    type(Node), pointer :: node_tal => null()
    type(NodeList), pointer :: node_nuc_list => null()
    type(Nuclide), pointer :: nuc => null()

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

    ! Parse input.xml file
    call open_xmldoc(doc, filename)

    ! Number of particles
    if (check_for_node(doc, "n_particles")) then
      call get_node_value(doc, "n_particles", n_particles)
    else
      message = "Number of particles not specified in input."
      call fatal_error()
    end if

    ! Get micro absorption cross section 
    if (check_for_node(doc, "micro_absxs")) then
      call get_node_value(doc, "micro_absxs", micro_absxs)
    else
      message = "Need to specify micro absorption xs in input."
      call fatal_error()
    end if

    ! Get list of nuclides
    call get_node_list(doc, "nuclide", node_nuc_list)

    ! Get the number of nuclides
    n_nuclides = get_list_size(node_nuc_list)
    allocate(nuclides(n_nuclides))

    ! Loop around nuclides and read all input
    do i = 1, n_nuclides

      ! Get pointer to i-th nuclide
      call get_list_item(node_nuc_list, i, node_nuc)
      nuc => nuclides(i)

      ! Check for density
      if (check_for_node(node_nuc, "density")) then
        call get_node_value(node_nuc, "density", density)
        call nuc % set_density(density)
      else
        message = "Missing density for nuclide."
        call fatal_error()
      end if

      ! Check for scattering xs
      if (check_for_node(node_nuc, "xs_file")) then
        call get_node_value(node_nuc, "xs_file", xs_file)
        call nuc % set_xs_file(xs_file)
      else
        message = "Missing cross section file in nuclide."
        call fatal_error()
      end if
    end do

    ! Get tally information
    if (check_for_node(doc, "tally")) then
      call get_node_ptr(doc, "tally", node_tal)

      if (check_for_node(node_tal, "type")) then
        call get_node_value(node_tal, "type", tally_type)
        select case (tally_type)
        case ('equal-lethargy')
          call flux_tal % set_type(EQUAL_LETHARGY)
          call abs_tal % set_type(EQUAL_LETHARGY)
          call scat_tal % set_type(EQUAL_LETHARGY)
          call r2c_tal % set_type(EQUAL_LETHARGY)
          call outscatc_tal % set_type(EQUAL_LETHARGY)
          call winscatc_tal % set_type(EQUAL_LETHARGY)
          call wc_tal % set_type(EQUAL_LETHARGY)
          call p1_scat_tal % set_type(EQUAL_LETHARGY)
          if (check_for_node(node_tal, "nbins")) then
            call get_node_value(node_tal, "nbins", tally_nbins)
            call flux_tal % set_nbins(tally_nbins)
            call abs_tal % set_nbins(tally_nbins)
            call scat_tal % set_nbins(tally_nbins)
            call r2c_tal % set_nbins(tally_nbins)
            call outscatc_tal % set_nbins(tally_nbins)
            call winscatc_tal % set_nbins(tally_nbins)
            call wc_tal % set_nbins(tally_nbins)
            call p1_scat_tal % set_nbins(tally_nbins)
          else
            message = 'Must specify tally bins.'
            call fatal_error()
          end if

        case default
          message = 'Only equal-lethargy allowed.'
          call fatal_error()
        end select
      else
        message = 'Need to specify tally type.'
        call fatal_error() 
      end if

    else
      message = 'Must specify tally information.'
      call fatal_error()
    endif

    ! Initialize tallies
    call flux_tal % initialize()
    call abs_tal % initialize()
    call scat_tal % initialize()
    call r2c_tal % initialize()
    call outscatc_tal % initialize()
    call winscatc_tal % initialize()
    call wc_tal % initialize()
    call p1_scat_tal % initialize()

    ! Close input XML file
    call close_xmldoc(doc)

    ! Read in cross section data
    call read_xs()

  end subroutine read_input_xml

!===============================================================================
! READ_XS
!===============================================================================

  subroutine read_xs()

    character(len=MAX_LINE_LEN) :: filename
    character(len=MAX_WORD_LEN) :: name
    integer :: i
    integer :: npts
    logical :: file_exists
    real(8) :: mass
    real(8), allocatable :: E(:)
    real(8), allocatable :: xs_a(:)
    real(8), allocatable :: xs_s(:)
    type(Node), pointer :: doc => null()
    type(Nuclide), pointer :: nuc => null()

    ! Loop around nuclides
    do i = 1, n_nuclides

      nuc => nuclides(i)

      ! Check if xs file exists 
      filename = nuc % get_xs_file()
      inquire(FILE=filename, EXIST=file_exists)
      if (.not. file_exists) then
        message = "XS file does not exist."
        call fatal_error()
      end if

      ! Parse xs file
      call open_xmldoc(doc, filename)

      ! Read in name field
      call get_node_value(doc, "name", name)
      call nuc % set_name(name)

      ! Read in atomic weight 
      call get_node_value(doc, "mass", mass)
      call nuc % set_mass(mass)

      ! Get size of grid
      npts = get_arraysize_double(doc, "E")
      allocate(E(npts))
      allocate(xs_s(npts))
      allocate(xs_a(npts))
      call get_node_array(doc, "E", E)
      call get_node_array(doc, "xs_s", xs_s)
      xs_a = micro_absxs

      ! Set data in nuclide
      call nuc % set_energy(E)
      call nuc % set_xs_s(xs_s)
      call nuc % set_xs_a(xs_a)

      ! Close xs file
      call close_xmldoc(doc)

      ! Deallocate arrays
      deallocate(E)
      deallocate(xs_s)
      deallocate(xs_a)

    end do

  end subroutine read_xs

end module input_xml
