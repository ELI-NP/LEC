      MODULE shared_data
      USE mpi
      IMPLICIT NONE

      TYPE particle
         REAL(num),DIMENSION(3) :: part_p
         REAL(num),DIMENSION(c_ndims) :: part_pos
         REAL(num) :: weight
         REAL(num) :: charge
         REAL(num) :: mass
         TYPE(particle), POINTER :: next, prev
         INTEGER :: coll_count
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
         REAL(num) :: particle_energy
#endif
      END TYPE particle

      TYPE particle_list
         TYPE(particle), POINTER :: head
         TYPE(particle), POINTER :: tail
         INTEGER(i8) :: count
         INTEGER :: id_update
         ! Pointer is safe if the particles in it are all unambiguously linked
         LOGICAL :: safe

         ! Does this partlist hold copies of particles rather than originals
         LOGICAL :: holds_copies
         TYPE(particle_list), POINTER :: next, prev
      END TYPE particle_list

      TYPE particle_species
         CHARACTER(string_length) :: name
         TYPE(particle_species), POINTER :: next, prev
         REAL(num) :: charge
         REAL(num) :: mass
         REAL(num) :: weight
         INTEGER(i8) :: count
         TYPE(particle_list) :: attached_list
         ! Secondary list
         TYPE(particle_list),DIMENSION(:,:), POINTER :: secondary_list
         REAL(num) :: npart_per_cell
      END TYPE particle_species
      END MODULE shared_data
!-------------------------------------------
      SUBROUTINE init_particle(new_particle)
!-------------------------------------------
      TYPE(particle), POINTER :: new_particle

      new_particle%part_p = 0.0_num
      new_particle%part_pos = 0.0_num
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
      new_particle%weight = 0.0_num
#endif
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
    ! This assigns an optical depth to newly created particle
      new_particle%particle_energy = 0.0_num
#endif
#ifdef PHOTONS
      new_particle%optical_depth = LOG(1.0_num / (1.0_num - random()))
#endif
      END SUBROUTINE init_particle
!---------------------------------------------
      SUBROUTINE create_particle(new_particle)
!---------------------------------------------
      TYPE(particle), POINTER :: new_particle

      ALLOCATE(new_particle)
      CALL init_particle(new_particle)

      END SUBROUTINE create_particle
!---------------------------------------------------------------
      SUBROUTINE add_particle_to_partlist(partlist,new_particle)
!---------------------------------------------------------------
      TYPE(particle_list), INTENT(INOUT) :: partlist
      TYPE(particle), POINTER :: new_particle

      ! Note that this will work even if you are using an unsafe particle list
      ! BE CAREFUL if doing so, it can cause unexpected behaviour

      ! if (!particle) return;
      IF (.NOT. ASSOCIATED(new_particle)) RETURN
      NULLIFY(new_particle%next, new_particle%prev)

      ! Add particle count
      partlist%count = partlist%count + 1
      partlist%id_update = 1
      IF(.NOT. ASSOCIATED(partlist%tail)) THEN
         ! partlist is empty
        partlist%head => new_particle
        partlist%tail => new_particle
        RETURN
      END IF

      partlist%tail%next => new_particle
      new_particle%prev => partlist%tail
      NULLIFY(new_particle%next)
      partlist%tail => new_particle

      END SUBROUTINE add_particle_to_partlist
!------------------------------------------------------------------
      SUBROUTINE generate_photon(generating_electron, iphoton, eta)
!------------------------------------------------------------------
      ! Generates a photon moving in same direction as electron
      ! (generates entirely new photon)

      TYPE(particle), POINTER :: generating_electron
      INTEGER, INTENT(IN) :: iphoton
      REAL(num), INTENT(IN) :: eta
      REAL(num) :: dir_x, dir_y, dir_z, mag_p, generating_gamma
      REAL(num) :: rand_temp, photon_energy
      TYPE(particle), POINTER :: new_photon

      mag_p = MAX(SQRT(generating_electron%part_p(1)**2 &
              + generating_electron%part_p(2)**2 &
              + generating_electron%part_p(3)**2), c_tiny)

      dir_x = generating_electron%part_p(1) / mag_p
      dir_y = generating_electron%part_p(2) / mag_p
      dir_z = generating_electron%part_p(3) / mag_p

      generating_gamma = SQRT(1.0_num + (mag_p / m0 / c)**2)

      ! Determine photon energy

      rand_temp = random()
      photon_energy = calculate_photon_energy(rand_temp, eta, generating_gamma)

      IF(use_radiation_reaction) THEN
      ! Calculate electron recoil
         mag_p = mag_p - photon_energy / c

         generating_electron%part_p(1) = dir_x * mag_p
         generating_electron%part_p(2) = dir_y * mag_p
         generating_electron%part_p(3) = dir_z * mag_p
      END IF

      ! This will only create photons that have energies above a user specified
      ! cutoff and if photon generation is turned on. E+/- recoil is always
      ! considered
      IF (photon_energy > photon_energy_min .AND. produce_photons) THEN
        IF (photon_energy < c_tiny) photon_energy = c_tiny

        CALL create_particle(new_photon)
        new_photon%part_pos = generating_electron%part_pos

        new_photon%part_p(1) = dir_x * photon_energy / c
        new_photon%part_p(2) = dir_y * photon_energy / c
        new_photon%part_p(3) = dir_z * photon_energy / c

        new_photon%optical_depth = reset_optical_depth()
        new_photon%particle_energy = photon_energy
        new_photon%weight = generating_electron%weight

        CALL add_particle_to_partlist(species_list(iphoton)%attached_list, &
          new_photon)
      END IF
      END SUBROUTINE generate_photon
!--------------------------------------
      SUBROUTINE push_photons(ispecies)
!--------------------------------------
      USE shared_data
      REAL(num) :: delta_x, delta_y
      INTEGER,INTENT(IN) :: ispecies
      TYPE(particle), POINTER :: current
      REAL(num) :: current_energy, dtfac, fac
      dtfac = dt * c**2
      ! set current to point to head of list
      current => species_list(ispecies)%attached_list%head

      ! loop over photons
      DO WHILE(ASSOCIATED(current))

      ! Note that this is the energy of a single REAL particle in the
      ! pseudoparticle, NOT the energy of the pseudoparticle
         current_energy = current%particle_energy

         fac = dtfac/current_energy
         delta_x = current%part_p(1)*fac
         delta_y = current%part_p(2)*fac
         current%part_pos = current%part_pos + (/delta_x, delta_y/)
         current => current%next
      END DO
      END SUBROUTINE push_photons
