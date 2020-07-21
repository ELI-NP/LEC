! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
!  Overall period>2^123;
!
      FUNCTION random()
      IMPLICIT NONE
      INTEGER(kind=4) :: kiss,x,y,z,w
      REAL(kind=8) :: random
	COMMON /kisscom/x,y,z,w
C
      x = 69069*x + 1327217885
      y = IEOR(y,ISHFT(y,13))
      y = IEOR(y,ISHFT(y,-17))
      y = IEOR(y,ISHFT(y,5))
      z = 18000*IAND(z,65535) + ISHFT(z,-16)
      w = 30903*IAND(w,65535) + ISHFT(w,-16)
      kiss = x + y + ISHFT(z,16) + w
      random = (DBLE(kiss) + 2147483648.d0)/4294967296.d0
      RETURN
	END
!-------------------------------
      SUBROUTINE rand_init(seed)
!-------------------------------
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: seed
      INTEGER,PARAMETER :: init_x = 123456789
      INTEGER,PARAMETER :: init_y = 362436069
      INTEGER,PARAMETER :: init_z = 521288629
      INTEGER,PARAMETER :: init_w = 916191069
      REAL(kind=4) :: x,y,z,w
	COMMON /kisscom/x,y,z,w
c
      x = init_x + seed
      y = init_y + seed
      z = init_z + seed
      w = init_w + seed

      END SUBROUTINE
c
c     Uncomment the following section to test the random number generator
c     gfortran random.f ; for mpi processes, use mpif90 random.f
c     ./a.out ; for mpi process, use mpirun -np 4 ./a.out
c      PROGRAM main
c	IMPLICIT NONE
c      INCLUDE "mpif.h"
c	INTEGER :: i,seed
c      INTEGER :: myrank,ierr,nprocs
c	REAL(kind=8) :: rand,rand2,random
c      CHARACTER(LEN=100):: fo_name2
c
c      CALL MPI_INIT(IERR)
c      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
c      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)
c	seed = 123456789*(1 + myrank)
c 	CALL rand_init(seed)
c	WRITE(*,*) myrank, seed
c      WRITE(fo_name2,444) 'random',myrank
c      OPEN (1,file=fo_name2,form='formatted',status='unknown')
c
c	DO i = 1,700
c         WRITE(1,*) random(), random()
c	END DO
c
c      CLOSE(1)
c      CALL MPI_FINALIZE(IERR)
c444   FORMAT(A,I3.3,'.dat')
c	STOP
c      END
