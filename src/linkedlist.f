      MODULE partlist
      USE P_common
      IMPLICIT NONE
      SAVE
      CONTAINS
!---------------------------------------------
      SUBROUTINE create(new_particle)
!---------------------------------------------
      USE P_common
      TYPE(particle),POINTER :: new_particle

      ALLOCATE(new_particle)
      new_particle%momentum = 0.d0
      new_particle%position = 0.d0
      new_particle%weight = 0.d0

      END SUBROUTINE create
!--------------------------------------------------
      SUBROUTINE add_to_list(partlist,new_particle)
!--------------------------------------------------
      USE P_common
      TYPE(particle_list), INTENT(INOUT) :: partlist
      TYPE(particle), POINTER :: new_particle

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

      END SUBROUTINE add_to_list
      END MODULE partlist
!--------------------------------------------------
      SUBROUTINE create_photon(iphoton)
!--------------------------------------------------
      USE sim_common
      USE mpi_common
      USE R_common
      USE P_common
      USE partlist
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: iphoton
      TYPE(particle),POINTER :: new_photon
      TYPE(particle_species),DIMENSION(:),POINTER :: species_list
      CALL create(new_photon)
      new_photon%position(1) = Xe
      new_photon%position(2) = Ye
      new_photon%position(3) = Ze
      new_photon%momentum(1) = Vx0*ENN
      new_photon%momentum(2) = Vy0*ENN
      new_photon%momentum(3) = Vz0*ENN
      new_photon%energy = ENN*ENE0
      new_photon%weight = wight(i)
c      CALL add_to_list(species_list(iphoton)%attached_list
c     &                      ,new_photon)
      IF(myrank.EQ.0) WRITE(*,*) new_photon%energy
      END SUBROUTINE create_photon
!--------------------------------------
c      SUBROUTINE push_photons(ispecies)
!--------------------------------------
c      USE P_common
c      IMPLICIT NONE
c      REAL(kind=8) :: delta_x, delta_y, delta_z
c      INTEGER,INTENT(IN) :: ispecies
c      TYPE(particle), POINTER :: current
c      REAL(kind=8) :: current_energy, dtfac, fac
c      dtfac = dt * c**2
c      ! set current to point to head of list
c      current => species_list(ispecies)%attached_list%head

      ! loop over photons
c      DO WHILE(ASSOCIATED(current))

      ! Note that this is the energy of a single REAL particle in the
      ! pseudoparticle, NOT the energy of the pseudoparticle
c         current_energy = current%particle_energy

c         fac = dtfac/current_energy
c         delta_x = current%momentum(1)*fac
c         delta_y = current%momentum(2)*fac
c         delta_z = current%momentum(3)*fac
c         current%position(1) = current%position(1) + delta_x
c         current%position(2) = current%position(2) + delta_y
c         current%position(3) = current%position(3) + delta_z
c         current => current%next
c      END DO
c      END SUBROUTINE push_photons
