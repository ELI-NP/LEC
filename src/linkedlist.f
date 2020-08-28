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
      SUBROUTINE create_photon
!--------------------------------------------------
      USE sim_common
      USE mpi_common
      USE R_common
      USE P_common
      USE partlist
      IMPLICIT NONE
      TYPE(particle),POINTER :: new_photon

      CALL create(new_photon)
      new_photon%position(1) = Xe
      new_photon%position(2) = Ye
      new_photon%position(3) = Ze
      new_photon%momentum(1) = Vx0*ENN
      new_photon%momentum(2) = Vy0*ENN
      new_photon%momentum(3) = Vz0*ENN
      new_photon%energy = ENN*ENE0
      new_photon%weight = wight(i)
      CALL add_to_list(species_list(1)%attached_list
     &                      ,new_photon)
c      IF(myrank.EQ.0) WRITE(*,*) species_list(1)%attached_list
      END SUBROUTINE create_photon
!--------------------------------------
      SUBROUTINE push_photons
!--------------------------------------
      USE P_common
      USE sim_common
      USE out_common
      IMPLICIT NONE

      REAL(kind=8) :: delta_x, delta_y, delta_z,photon_x,photon_y
      TYPE(particle), POINTER :: current
      REAL(kind=8) :: current_energy, dtfac, fac

      WRITE(fo_name2,444) TRIM(data_file)//'photon',jobno
      OPEN(34,file=fo_name2,form='formatted')

c      ! set current to point to head of list
      current => species_list(1)%attached_list%head

      ! loop over photons
      DO WHILE(ASSOCIATED(current))

      ! Note that this is the energy of a single REAL particle in the
      ! pseudoparticle, NOT the energy of the pseudoparticle
         current_energy = current%energy

         fac = dt
         TT = dsqrt(Vx**2+Vy**2+Vz**2)
         delta_x = Vx/TT*fac
         delta_y = Vy/TT*fac
         delta_z = Vz/TT*fac
         current%position(1) = current%position(1) + delta_x
         current%position(2) = current%position(2) + delta_y
         current%position(3) = current%position(3) + delta_z

         photon_x = current%position(1)
         photon_y = current%position(2)

         IF(MOD(kstep,ksout).EQ.0) THEN
           WRITE(34,666) photon_x*Rx, photon_y*Rx
         END IF
         current => current%next
      END DO



444   format(A,I3.3,'.dat')
666   format(11(E16.6,1X))
      END SUBROUTINE push_photons
