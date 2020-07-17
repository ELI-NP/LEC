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
      y = ieor(y,ishft(y,13))
      y = ieor(y,ishft(y,-17))
      y = ieor(y,ishft(y,5))
      z = 18000*iand(z,65535) + ishft(z,-16)
      w = 30903*iand(w,65535) + ishft(w,-16)
      kiss = x + y + ishft (z, 16) + w
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
c
      x = init_x + seed
      y = init_y + seed
      z = init_z + seed
      w = init_w + seed

      END SUBROUTINE
c
      ! Uncomment the following section to test the random number generator
      ! gfortran random.f
      ! ./a.out
c     PROGRAM main
c	IMPLICIT NONE
c	INTEGER :: i
c	REAL(kind=8) :: rand,rand2,random
c	CALL rand_init(94561)
c
c	DO i = 1,10
c        rand = random()
c	   WRITE(*,*) rand
c	END DO
c	STOP
c     END
