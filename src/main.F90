program main

  use constants,  only: MAX_COLLISIONS
  use cross_section,  only: calculate_xs
  use finalize,  only: finalize_run
  use global,  only: n_particles
  use initialize,  only: initialize_run
  use particle_class,  only: Particle
  use physics, only: collide
  use tracking, only: transport

  implicit none

  integer :: i
  integer :: j
  type(Particle) :: p

  ! Inititalize
  call initialize_run()

  ! Begin loop around particles
  PARTICLE_LOOP: do i = 1, n_particles

    ! Get initial particle information
    call p % start()

    ! Start collision loop
    COLLISION_LOOP: do j = 1, MAX_COLLISIONS

      ! Calculate cross sections
      call calculate_xs(p)

      ! Transport particle
      call transport(p)

      ! Perform collision physics
      call collide(p)

      ! Check if particle is still alive
      if (.not. p % get_alive()) then
        call p % set_n_collisions(j)
        exit COLLISION_LOOP
      end if

    end do COLLISION_LOOP

  end do PARTICLE_LOOP

  ! Finalize
  call finalize_run()

end program main
