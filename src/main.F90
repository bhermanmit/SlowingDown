program main

  use constants,  only: MAX_COLLISIONS, MIN_ENERGY
  use cross_section,  only: calculate_xs
  use finalize,  only: finalize_run
  use global,  only: n_particles
  use initialize,  only: initialize_run
  use particle_class,  only: Particle
  use physics, only: run_physics

  implicit none

  integer :: i
  integer :: j
  type(Particle) :: p

  ! Inititalize
  call initialize_run(p)

  ! Begin loop around particles
  PARTICLE_LOOP: do i = 1, n_particles

    ! Print particle id in increments of 1000
    if (mod(i, 100000) == 0) print *, 'At Particle:', i

    ! Get initial particle information
    call p % start()

    ! Start collision loop
    COLLISION_LOOP: do j = 1, MAX_COLLISIONS

      ! Start the collision
      call p % begin_collision()

      ! Calculate cross sections
      call calculate_xs(p)

      ! Perform transport and collision physics
      call run_physics(p)

      ! Tally results
      call p % tally()

      ! Kill neutron if below min energy
      if (p % get_energy() <= MIN_ENERGY) call p % set_alive(.false.)

      ! Check if particle is still alive
      if (.not. p % get_alive()) then
        call p % set_n_collisions(j)
        exit COLLISION_LOOP
      end if

    end do COLLISION_LOOP

  end do PARTICLE_LOOP

  ! Finalize
  call finalize_run()
  call p % clear()

end program main
