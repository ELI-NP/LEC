      MODULE random_common
      INTEGER,PARAMETER :: icpu = 4
      INTEGER,PARAMETER      :: LE0 = 1000   ,LE1=1000
      REAL(kind=8),PARAMETER :: Zcom = 79.d0   ! Z component for nucl
      REAL(kind=8),PARAMETER :: Zcm3 = 4.3d0   ! zcm3 = zcom**(1/3)
c     REAL(kind=8),PARAMETER :: Zcom = 13.d0  ! Z component for nucl
c     REAL(kind=8),PARAMETER :: Zcm3 = 2.35d0 ! zcm3 = zcom**(1/3)
      INTEGER,PARAMETER      :: LPx = 200
      REAL(kind=8),PARAMETER :: BPx = 200.d0
      REAL(kind=8),PARAMETER :: Emax = 10000.d0
c
      REAL(kind=8) :: dp1
      REAL(kind=8),SAVE,DIMENSION(LE0,LPx) :: resultL
      REAL(kind=8),SAVE,DIMENSION(LPx) :: totalH,totalH2
      END MODULE random_common
c
      MODULE sim_common
      USE random_common
      INTEGER      :: i,j,k,ksmax,ksout,kk,kstep,ii,ipl
     &               ,itotal,ksoutP,iconR,SKL,LL,polar,shape
     &               ,photon,species,shot,itotal0,L,OutRad,OutPairs
     &		   ,QED,sampled,sampled2,sampled3,sampled4,load
      LOGICAL      :: emmits,exists,load_particle
      REAL(kind=8) :: Xe,Ye,Ze,T,VX,VY,VZ
      REAL(kind=8) :: AVEX,AVEY,AVEZ,AVBX,AVBY,AVBZ
      REAL(kind=8) :: PXS,PYS,PZS,BXT,BYT,BZT
      REAL(kind=8) :: Rx,Rm,Rt,Rd1,Rd2,Rt2
      REAL(kind=8) :: pin,rin,div,div2,xinit,bin,enum,rmass
      REAL(kind=8) :: Wx0,Wy0,Wz0,Vx0,Vy0,Vz0,Wx,Wy,Wz
      REAL(kind=8) :: gg,ff,dF0,dVx,dVy,dVz
      REAL(kind=8) :: fvx,fvy,fvz
      REAL(kind=8) :: gammi,TTT,ENE00,RRR,RAD,TT1,GAMM
      REAL(kind=8) :: alf,palf,alf2,alpha,comp,sigmax,sigmay,sigmaz
      REAL(kind=8) :: SL,EL,BL,Rd,VL,Em,EV,inc_ang
      REAL(kind=8) :: DX,DT,DTI,DT1,phaseX,phaseY,phaseZ,TT
      REAL(kind=8) :: ww,wp,wk,pp,pw,pi,pi2,ppt,sp,ws,lb,wb
      REAL(kind=8) :: we0i,wpi,wsi,Pksout,wmin,www,Tm,Um,FF1,ENEd
      REAL(kind=8) :: ENEh,ENEv,ENE0,ENE,EmaxV,we0,XI,ENN,p_x,ENEA
      REAL(kind=8) :: x0,y0,z0,Tx,Ty,Tz,percentage
      REAL(kind=8) :: PQM       ! charge/mass ; rmass=-1.d0 for electron
      REAL(kind=8),DIMENSION(:,:),ALLOCATABLE   :: RE,RH
      REAL(kind=8) :: ERD,EKE,EKK,ERK
      REAL(kind=8),DIMENSION(0:3000 + 1,200):: diffC,diffQ,diffD,diffR
      REAL(kind=8),DIMENSION(200) :: totalR,totalC,totalP,totalRC
     &                              ,totalS,totalT
      REAL(kind=8) :: rand,rand1
      INTEGER,SAVE :: tmove,tmove0,tmove1
     &		   ,trdct,trdct0,trdct1
     &		   ,tmall,tmall0,t_max
     &		   ,tinit,tinit0,t_rate
     &		   ,tcurr,tcurr0,tcurr1
      END MODULE sim_common
c
      MODULE mpi_common
      INTEGER :: myrank,ierr,nprocs,thREADs
      END MODULE mpi_common
c
      MODULE out_common
      INTEGER :: jobno
      CHARACTER(LEN=100):: fo_name2,cwd
      CHARACTER(LEN=10) :: filename = 'input.dat'
      CHARACTER(LEN=100):: data_dir,data_file,input_file
      CHARACTER(LEN=*), PARAMETER :: data_dir_file =
     &					'USE_DATA_DIRECTORY'
      END MODULE out_common

      MODULE R_common
      INTEGER      :: Lall,LP,IPTSS,IPTFF,jj
      INTEGER,DIMENSION(7) :: Ne7
      REAL(kind=8) :: wight00,const
      REAL(kind=8),DIMENSION(:),ALLOCATABLE :: wight0
      REAL(kind=8),DIMENSION(:),ALLOCATABLE :: wight
      REAL(kind=8) :: TmY,TmZ
      REAL(kind=8),DIMENSION(:,:,:), ALLOCATABLE:: phtn
      REAL(kind=8) :: wmitT3,vmitT3,wmitTT,vmitTT,hh
      REAL(kind=8),DIMENSION(:),ALLOCATABLE :: wmit3,vmit3
      REAL(kind=8),DIMENSION(:,:),ALLOCATABLE :: emit3,fmit3,qmit3
      REAL(kind=8),DIMENSION(:),ALLOCATABLE :: emitT3,fmitT3,qmitT3
      REAL(kind=8),DIMENSION(:),ALLOCATABLE :: emitTT,fmitTT,qmitTT
      END MODULE R_common
c
      program main
      USE random_common
      USE sim_common
      USE mpi_common
      USE R_common
      USE out_common
      USE omp_lib
      IMPLICIT NONE
      INCLUDE "mpif.h"
c
      NAMELIST /PARAM1/ jobno,ksmax,div,div2			       !time & grid PARAMETERs
      NAMELIST /PARAM2/ SL,Ev,pw,pp,sp				       !laser PARAMETERs
      NAMELIST /PARAM3/ alpha,enum,bin,shot,inc_ang		       !electron PARAMETERs
      NAMELIST /PARAM4/ xinit,rmass,sigmax,sigmay,sigmaz           !configurations
      NAMELIST /PARAM5/ iconR,QED,ipl,shape,load,OutRad,OutPairs   !physical processes
c
      CALL system_clock(tmall0)
      CALL MPI_INIT(IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)	!initialize MPI
c
      tmove = 0 ; tcurr = 0 ; trdct = 0
      tmall = 0 ; tinit = 0
c
      PI  = 4.D0*DATAN(1.D0)
      PI2 = PI*2.D0
c
      CALL getcwd(cwd)
      IF (myrank.EQ.0) THEN
         OPEN(unit=41, status='OLD', file = TRIM(data_dir_file)
     &	, iostat=ierr)
         IF (ierr.EQ.0) THEN
      	    READ(41,'(A)') data_dir
            CLOSE(41)
      	    PRINT*, 'Using data directory "'
     &	          // TRIM(data_dir) // '"'
         ELSE
     	    PRINT*, 'SpecIFy output directory'
      	    READ(*,'(A)') data_dir
         END IF
      END IF
      CALL MPI_BCAST(data_dir,100, MPI_CHARACTER, 0,MPI_COMM_WORLD,ierr)

      input_file = TRIM(ADJUSTL(data_dir))//'/'
     &	    // TRIM(ADJUSTL(filename))

      data_file =  TRIM(ADJUSTL(data_dir))//'/'
c
      INQUIRE(file=input_file, exist=exists)
      IF (.NOT. exists .AND. myrank == 0) THEN
        PRINT *, '*** ERROR ***'
        PRINT *, 'Input deck file "' // TRIM(input_file)
     &            // '" DOes not exist.'
        PRINT *, 'Create the file and rerun the code.'
      END IF
c
      IF(myrank.EQ.0) THEN
         WRITE(*,*) 'Output directory:',data_file
      END IF

      OPEN(8,status="old",file=TRIM(input_file),form='formatted')	!READ input data file
      READ(8,PARAM1)
c
      jobno = jobno + myrank
      WRITE(fo_name2,444) TRIM(ADJUSTL(data_file))//'output',jobno
      OPEN(9,file=fo_name2,form='formatted')
c
      WRITE(9,PARAM1)
      READ(8,PARAM2) ; WRITE(9,PARAM2)
      READ(8,PARAM3) ; WRITE(9,PARAM3)
      READ(8,PARAM4) ; WRITE(9,PARAM4)
      READ(8,PARAM5) ; WRITE(9,PARAM5)
      CLOSE(8)
c
      WRITE(9,*) "Check OPENMP"
!$omp parallel
      thREADs = omp_get_num_thREADs()
!$omp END parallel
      WRITE(9,*) "total number of thREADs=",thREADs
c
      CALL setprm  						!PARAMETER setup
      CALL welcome
c
      CALL system_clock(tinit0)
      CALL setbeam
c
      WRITE(fo_name2,444) TRIM(data_file)//'orbt1q',jobno
      OPEN (10,file=fo_name2,form='formatted',status='REPLACE')
      WRITE(fo_name2,444) TRIM(data_file)//'orbt2q',jobno
      OPEN (11,file=fo_name2,form='formatted',status='REPLACE')
      WRITE(fo_name2,444) TRIM(data_file)//'orbt3q',jobno
      OPEN (12,file=fo_name2,form='formatted',status='REPLACE')
      WRITE(fo_name2,444) TRIM(data_file)//'orbt4q',jobno
      OPEN (13,file=fo_name2,form='formatted',status='REPLACE')
      WRITE(fo_name2,444) TRIM(data_file)//'orbt5q',jobno
      OPEN (14,file=fo_name2,form='formatted',status='REPLACE')
      WRITE(fo_name2,444) TRIM(data_file)//'orbt6q',jobno
      OPEN (15,file=fo_name2,form='formatted',status='REPLACE')
      WRITE(fo_name2,444) TRIM(data_file)//'orbt7q',jobno
      OPEN (16,file=fo_name2,form='formatted',status='REPLACE')
      CALL system_clock(tinit)
      tinit = tinit - tinit0
c
c     orbit calculation
c
      T      = -wp*xinit
      RAD    = 0.d0
      TT1    = 0.d0
      kstep  = 0
c
1000  CONTINUE              			  ! time loop
	KSTEP = KSTEP + 1
        T = T + dt
        CALL system_clock(tmove0)
        IF(iconR.EQ.0) CALL Lorentz
        IF(iconR.EQ.1) CALL Sokolov
        IF(iconR.EQ.2) CALL Landau_Lifshitz
        IF(iconR.EQ.3) CALL qedemmit	  !particle loop
        CALL out_orbit 				  !output extracted electron orbit
        CALL system_clock(tmove1)
        tmove = tmove + (tmove1 - tmove0)
c
        IF(OutRad.EQ.1) CALL store_orbit    !calculate emission
        IF((MOD(kstep,ksout).EQ.0)
     &     .AND.(alpha.GT.0.d0)  ) CALL average_energy !calculate average energy
c
        IF((MOD(kstep,ksmax).EQ.0)
     &     .AND.(alpha.GT.0.d0)  ) CALL histogram
c
        IF((MOD(kstep,ksmax).EQ.0)
     &     .AND.(alpha.GT.0.d0)  ) CALL histogram2d
c
        IF((MOD(kstep,1000).EQ.0)
     &     .AND.(myrank.EQ.0))
     &	WRITE(*,*) "Iteration =", kstep,"; Time step =", Rt*kstep
c---------------------------------------------------------
      IF(kstep.GE.ksmax) GO TO 2500
      GO TO 1000					   ! END time loop
c
      CLOSE(10) ; CLOSE(11) ; CLOSE(12)
      CLOSE(13) ; CLOSE(14) ; CLOSE(15)
      CLOSE(16) ;	CLOSE(20)
c
2500  CONTINUE
      DEALLOCATE(RE,RH)
      IF(OutRad.EQ.1) THEN
         IF(iconR.EQ.3) THEN
            CALL photon_his
         ELSE
            CALL radiation
         END IF
      END IF
c
      IF(OutRad.EQ.1) DEALLOCATE(phtn)
      DEALLOCATE(wight0,wight)
c
333   FORMAT(A,I2.2,'_',I2.2,'.dat')
444   FORMAT(A,I3.3,'.dat')
555   FORMAT(A,I2.2,'_',I2.2,'_',I2.2,'.dat')
666   FORMAT(11(E16.6,1X))
c
9999  CONTINUE
c
      CALL system_clock(tmall,t_rate,t_max)
      tmall = tmall - tmall0
      WRITE(9,*) "all      :", DBLE(tmall)/DBLE(t_rate)
      WRITE(9,*) "initc    :", DBLE(tinit)/DBLE(t_rate)
      WRITE(9,*) "pmove    :", DBLE(tmove)/DBLE(t_rate)
      WRITE(9,*) "rad      :", DBLE(tcurr)/DBLE(t_rate)
      WRITE(9,*) "reduction(rad):", DBLE(trdct)/DBLE(t_rate)
      WRITE(9,*) "others   :", DBLE(tmall - tinit - tmove - tcurr
     &                        -trdct)/DBLE(t_rate)
      CLOSE(9)
      CALL MPI_FINALIZE(IERR)
      IF(myrank.EQ.0)
     &   WRITE(*,*)
     &   "Final runtime =",tmall/t_rate,"seconds"
c
      STOP
      END
!--------------------------------------
      SUBROUTINE setprm
!--------------------------------------
      USE random_common
      USE sim_common
      USE mpi_common
      USE R_common
      USE out_common
      USE omp_lib
      IMPLICIT NONE
      CHARACTER(LEN=8)  :: qed_file = 'vac4.dat'
      CHARACTER(LEN=100):: table_location

      VL = 3.d0*1.0d8              ! light speed           [m/s]
      EL = 2.74d0*1.0d3 *dsqrt(SL) ! peak electric field   [V/m]
      BL = 9.28d0*1.0d - 2*dsqrt(SL) ! peak magnetic field   [Gauss]
      Rm = 1.7d0*1.0d1/BL          ! typical Larmor radius [m]
      Rd1 = Rm/div                 ! grid size             [m]
      Rd2 = pw/div2
      Rd = MIN(Rd1,Rd2)
      Rx = Rm                      ! unit scale            [m]
      Rt = 0.5d0*Rd/VL             ! time step             [s]
      Rt2 = Rt*div                 ! unit time             [s]
      dx = DBLE(1.d0/div)*Rd/Rd1   ! normalized grid size
      dt = 0.5d0*dx                ! normalized time step
      Em = Ev/(0.511d6*rmass)      ! nomalized beam energy
      Vx = dsqrt((Em + 1.d0)**2 - 1.d0)! normalized 4 - velocity (+x)
      Wx = Vx/(Em + 1.d0)          ! normalized 3 - velocity (+x)
      Vy = 0.d0                    ! normalized 4 - velocity (+y)
      Vz = 0.d0                    ! normalized 4 - velocity (+z)
      ppt = pp/VL
      ww = pw/Rd*dx                ! normalized wavelength
      wk = pi2/ww                  ! normalized wave number
      wp = pp/Rd*dx                ! normalized pulse length
      ws = sp/Rd*dx                ! normalized waist radius
      wpi = 1.d0/wp
      wsi = 1.d0/ws
      pin = bin/(0.511d6*rmass)     ! normalized bin size
      rin = Em/6000.d0              ! reference bin size
      wb = alpha/Rd*dx              ! normalized transverse beam size
      PQM =-1.d0/rmass
c
      WRITE(9,*) " "
      WRITE(9,*) "Parameters for pulse laser"
      WRITE(9,*) " "
      IF(ipl.EQ.0) THEN
         WRITE(9,*) "Laser polarization: linear"
         polar = 0.d0
      ELSE
         WRITE(9,*) "Laser polarization: circular"
         polar = 1.d0
      END IF
      WRITE(9,*) " "
      IF(shape.EQ.0) WRITE(9,*) "0th order Gaussian beam"
      IF(shape.EQ.1) WRITE(9,*) "5th order Paraxial approximation"
      WRITE(9,*) " "
      WRITE(9,*) "Laser Intensity           [W/cm2]", SL
      WRITE(9,*) "Peak electric field         [V/m]", EL
      WRITE(9,*) "Peak magnetic field       [Gauss]", BL
      WRITE(9,*) "Larmor radius for light speed [m]", Rm
      WRITE(9,*) "laser wavelength              [m]", pw
      WRITE(9,*) "pulse length                  [m]", pp
      WRITE(9,*) "pulse duration                [s]", ppt
      WRITE(9,*) "pulse duration (FWHM)         [s]", ppt*1.1774d0
      WRITE(9,*) "waist radius (1/e2)           [m]", sp
c
      WRITE(9,*) " "
      WRITE(9,*) "PARAMETERs for electron beam"
      WRITE(9,*) " "
      IF(load.EQ.1) THEN
        load_particle=.TRUE.
        WRITE(9,*) "Load particles from: f_E_smoothed_new_final.txt"
      ELSE
        load_particle=.FALSE.
      END IF
      WRITE(9,*) " "
      WRITE(9,*) "PARAMETERs for electron beam"
      WRITE(9,*) "initial kinetic energy [eV]", Ev
      WRITE(9,*) "initital momentum    [p/mc]", Vx
      WRITE(9,*) "initital 3 - velocity [v/c]", Wx
      IF(alpha.GT.1) THEN
        WRITE(9,*) "momentum spREAD x [%]", sigmax*100
        WRITE(9,*) "momentum spREAD y [%]", sigmay*100
        WRITE(9,*) "momentum spREAD z [%]", sigmaz*100
      END IF
      WRITE(9,*) "number of shot",shot
c
      WRITE(9,*) " "
      WRITE(9,*) "PARAMETERs for computation"
      WRITE(9,*) " "
      WRITE(9,*) "Larmor radius/grid separation      ", div
      WRITE(9,*) "Laser wavelength/grid separation   ", div2
      WRITE(9,*) "grid separation,Rd1,Rd2         [m]", Rd1, Rd2
      WRITE(9,*) "grid separation,Rd=MIN(Rd1,Rd2) [m]", Rd
      WRITE(9,*) "physical time step INTerval     [s]", Rt
      WRITE(9,*) "physical unit time              [s]", Rt2
      WRITE(9,*) "energy bin size of spectrum    [eV]", bin
c
      palf = 6.2d - 24
      alf  = palf*dt/Rt
      alf2 = alf*137.d0*9.d0/4.d0
      comp  = palf*137.d0*VL*3.d0*PI
      WRITE(9,*) "coefficient for radiation [s]", palf
      WRITE(9,*) "Compton wavelength        [m]", comp
      IF(iconR.EQ.0) THEN
        WRITE(9,*) " "
        WRITE(9,*) "Particle pusher: Lorentz"
        WRITE(9,*) " "
      ELSE IF(iconR.EQ.1) THEN
        WRITE(9,*) " "
        WRITE(9,*) "Particle pusher: Sokolov"
        WRITE(9,*) " "
      ELSE IF(iconR.EQ.2) THEN
        WRITE(9,*) " "
        WRITE(9,*) "Particle pusher: reduced Landau Lifshitz"
        WRITE(9,*) " "
      ELSE IF(iconR.EQ.3) THEN
        WRITE(9,*) " "
        WRITE(9,*) "Particle pusher: Lorentz"
        WRITE(9,*) "        USE QED: Stochastic"
        WRITE(9,*) " "
      END IF
c
      IF(QED.EQ.0) THEN
        WRITE(9,*) "USE QED: No, Classical"
        WRITE(9,*) " "
      ELSE
        WRITE(9,*) "USE QED: QED_assisted"
        WRITE(9,*) " "
      END IF
c
      WRITE(9,*) "max time step",ksmax
      ksout = 1
      ksoutP = 1
      IF(ksmax.GE.10000) ksout = ksmax/10000
      IF(ksmax.GE.1000000) ksoutP = ksmax/1000000
      Pksout= DBLE(ksoutP)
c setup for theoretical cross sections
c
      EmaxV= 1000.d0
      we0 =(log(EmaxV) - log(0.001d0))/199.0 !!
      we0i = 1.d0/we0
c
c radiation & pair production in strong - field
c
c     totalR : quantum   total cross section (x energy) ; fcem
c     totalC : classical total cross section (x energy) ; gcem
c     totalP : pair production total cross section      ; dcpr
c     totalS : quantum   total cross section            ; dcem
c     totalT : classical total cross section            ; ecem
c     diffC ; classical differential cross section (x energy)
c     diffQ ;   quantum differential cross section (x energy)
c     diffD ; classical differential cross section
c     diffR ;   quantum differential cross section
c
      CALL getcwd(cwd)
      table_location = TRIM(cwd)//'/TABLE/'//TRIM(ADJUSTL(qed_file))
      IF (myrank.EQ.0) THEN
      INQUIRE(file=TRIM(table_location), exist=exists)
      IF (.NOT.exists) THEN
         WRITE(*,*) '*** ERROR ***'
         WRITE(*,*) 'Unable to find QED tables in the ',
     &   'directory "' // TRIM(table_location) // '"'
      END IF
      END IF

      OPEN(96,file=table_location,form='unformatted')
      READ(96) totalR,totalC,totalP,totalS,totalT
      READ(96) diffC,diffQ,diffD,diffR
      CLOSE(96)

      IF(QED.EQ.1) THEN
        totalRC = totalR/totalC
      ELSE
        totalRC=1.d0
      END IF
c
      RETURN
      END
!--------------------------------------
      SUBROUTINE setbeam
!--------------------------------------
      USE random_common
      USE sim_common
      USE mpi_common
      USE R_common
      USE omp_lib
      USE out_common
      IMPLICIT NONE
      REAL(kind=8) :: aa
c generate initial electron conditions for incident beam
c---------------------------
      IF(alpha.NE.0.d0) THEN
c---------------------------
        sampled = 19 ; sampled2 = sampled/2
        sampled3 = (sampled - 1)**3

        ALLOCATE(Re(11,sampled3))
        ALLOCATE(wight0(sampled3))

        WRITE(9,*) "transverse beam size   [m]", alpha
        WRITE(9,*) "electron number per shot  ", enum
        WRITE(9,*) "incident angle    [degree]", inc_ang
        inc_ang = inc_ang/180.d0*pi

        DO k = 1,sampled - 1
        DO j = 1,sampled - 1
        DO i = 1,sampled - 1
           kk = (k - 1)*(sampled - 1)*(sampled - 1)
     &          + (j - 1)*(sampled - 1) + i
           CALL random_number(rand)
           CALL random_number(rand1)
           phaseX  = (DBLE(k - sampled2))/(sampled2 - 1)*0.707d0
           phaseY  = (DBLE(i - sampled2))/(sampled2 - 1)*0.707d0
           phaseZ  = (DBLE(j - sampled2))/(sampled2 - 1)*0.707d0
           Re(1,kk) = wp*xinit + phaseX*wb
           Re(2,kk) = phaseY*wb
           Re(3,kk) = phaseZ*wb
           wight0(kk) = dexp(-phaseX**2 - phaseY**2 - phaseZ**2)
c
           IF(load_particle) CALL manual_load
c
           Re(4,kk) = Vx*(-1.d0) + sigmax*dcos(2.d0*pi*rand1)
     &		            *Vx*dsqrt(-2.d0*log(rand))
           Re(5,kk) = sigmay*Vx*dcos(2.d0*pi*rand1)
     &		            *dsqrt(-2.d0*log(rand))
           Re(6,kk) = sigmaz*Vx*dsin(2.d0*pi*rand1)
     &	          	  *dsqrt(-2.d0*log(rand))
        END DO
        END DO
        END DO
        wight00 = 0.d0
        DO i = 1,sampled3
             wight00 = wight00 + wight0(i)
        END DO
        wight0 = wight0/wight00
c
        DO i = 1,sampled3
           Xe = Re(1,i)
           Ye = Re(2,i)
           Vx0 = Xe*dcos(inc_ang) - Ye*dsin(inc_ang)
           Vy0 = Xe*dsin(inc_ang) + Ye*dcos(inc_ang)
           Re(1,i) = Vx0
           Re(2,i) = Vy0
        END DO
c
        DO i = 1,sampled3
           Vx0 = Re(4,i)
           Vy0 = Re(5,i)
           Re(4,i) = Vx0*dcos(inc_ang)-Vy0*dsin(inc_ang)
           Re(5,i) = Vx0*dsin(inc_ang)+Vy0*dcos(inc_ang)
        END DO

        itotal = INT(sampled3/DBLE(nprocs))
        ALLOCATE(Rh(11,itotal))
        ALLOCATE(wight(itotal))
        jj = myrank*itotal
        WRITE(9,*) "myrank                            ", myrank
        WRITE(9,*) "total process number, nprocs      ", nprocs
        WRITE(9,*) "calculated particle number, itotal", itotal
        WRITE(9,*) "calculated particle number from",
     &              jj + 1, "to", jj + itotal
        DO i = 1,itotal
           Rh(1,i) = Re(1,i + jj)
           Rh(2,i) = Re(2,i + jj)
           Rh(3,i) = Re(3,i + jj)
           Rh(4,i) = Re(4,i + jj)
           Rh(5,i) = Re(5,i + jj)
           Rh(6,i) = Re(6,i + jj)
           wight(i)= wight0(i + jj)
        END DO
        Re = Rh
c
        CALL histogram
        CALL histogram2d
c----------
      ELSE
c----------
        WRITE(9,*) "single electron"
        WRITE(9,*) "electron number per shot = 1  "
        ALLOCATE(Re(11,1),Rh(11,1))
        ALLOCATE(wight0(1),wight(1))
        itotal = 1
        Xe = wp*xinit
        Ye = 0.d0
        Ze = 0.d0
        inc_ang = inc_ang/180.d0*pi
        Re(1,1)  = Xe*dcos(inc_ang) - Ye*dsin(inc_ang)
        Re(2,1)  = Xe*dsin(inc_ang) + Ye*dcos(inc_ang)
        Re(3,1)  = Ze
        Re(4,1)  = Vx*(-1.d0)*dcos(inc_ang)
        Re(5,1)  = Vx*(-1.d0)*dsin(inc_ang)
        Re(6,1)  = 0.d0
        wight(1) = 1.d0
      END IF
c
      DO i = 1,7
         Ne7(i)=(itotal/8)*i + 1
         WRITE(9,*) "sampling electron number", i, Ne7(i)
      END DO
c
      if(OutRad.EQ.1) THEN
         ALLOCATE(phtn(6,ksmax,itotal))
         phtn = 0.d0
      END IF

444   FORMAT(A,I4.4,'.dat')
666   FORMAT(3(E14.4,1X))
      RETURN
      END
!--------------------------------------
      SUBROUTINE store_orbit
!-------------------------------------
      USE random_common
      USE sim_common
      USE mpi_common
      USE R_common
      USE omp_lib
      IMPLICIT NONE
      REAL(kind=8) :: TTY,TTZ
c
      IF(MOD(kstep,ksoutP).EQ.0) THEN
         j = kstep/ksoutP
         IF(j.LE.1000000) THEN
            DO i=1,itotal
               Vx = Re(4,i)
               Vy = Re(5,i)
               Vz = Re(6,i)
               ff = Re(7,i)
               TTT = Re(9,i)
               Xi = Re(10,i)
               ENE = Re(11,i)
               TTY = dsqrt(Vx**2 + Vy**2)
               TmY = acos(Vx*(-1.d0)/TTY)
               IF(Vy.LT.0.d0) TmY = (-1.d0)*TmY
               TTZ = dsqrt(Vx**2 + Vz**2)
               TmZ = acos(Vx*(-1.d0)/TTZ)
               IF(Vz.LT.0.d0) TmZ=(-1.d0)*TmZ   ! direction angleZ
               phtn(1,j,i) = Xi
               phtn(2,j,i) = ENE
               phtn(3,j,i) = TmY
               phtn(4,j,i) = TTT
               phtn(5,j,i) = ff
               phtn(6,j,i) = TmZ
            END DO
        END IF
      END IF
      RETURN
      END
!--------------------------------------
      SUBROUTINE out_orbit
!--------------------------------------
      USE random_common
      USE sim_common
      USE mpi_common
      USE R_common
      USE omp_lib
      IMPLICIT NONE

      IF(MOD(kstep,ksout).EQ.0) THEN
      DO i = 1,7
         k = Ne7(i)
         Vx = Re(4,k)
         Vy = Re(5,k)
         Vz = Re(6,k)
         ENE = dsqrt(1.d0 + Vx**2 + Vy**2 + Vz**2)
         WRITE(9 + i,666) t*Rt2, Re(1,k)*Rx, Re(2,k)*Rx	!t, x, y
     &             ,Re(4,k), Re(5,k), (ENE)*0.511d6	      !px, py, K.E
     &	       ,abs(TTT)*0.511d6, RAD*0.511d6, Xi	      !Wem,Wrad,chi_e
      END DO
      END IF
 666  FORMAT(11(E16.6,1X))
      RETURN
      END
!------------------------------
      SUBROUTINE average_energy
!------------------------------
      USE random_common
      USE sim_common
      USE mpi_common
      USE R_common
      USE out_common
      USE omp_lib
      IMPLICIT NONE
      REAL(kind=8) :: SIG_neg,SIG_pos
     &		   ,SIGE_pos,SIGE_neg,sss
      INTEGER      :: num_p,num_n,num_pp,num_nn
      INTEGER,DIMENSION(icpu)      :: num_pos, num_neg
      REAL(kind=8),DIMENSION(icpu) :: SE,SE_pos,SE_neg
      INCLUDE "mpif.h"

      SE = 0.d0
      itotal0 = INT(itotal/icpu + 1)
      DO LP = 1,icpu
         IPTSS  = itotal0*(LP - 1) + 1
         IPTFF  = MIN(itotal0* LP,itotal)
         DO L  = IPTSS,IPTFF
            Vx  = Re(4,L)
            Vy  = Re(5,L)
            Vz  = Re(6,L)
            GAMM = sqrt(1.d0 + Vx*Vx + Vy*Vy + Vz*Vz)
            SE(LP) = SE(LP) + GAMM
         END DO
      END DO
c
      GAMM = 0.d0
      DO LP=1,icpu,4
         GAMM = GAMM + SE(LP) + SE(LP + 1)
     &          + SE(LP + 2) + SE(LP + 3)
      END DO

      EKK = 0.d0
c---------------------------------------------------------
      CALL mpi_allreduce(GAMM,EKK,1,MPI_REAL8,MPI_SUM
     &                   ,mpi_comm_world ,ierr)
c
      EKK = EKK/(itotal*nprocs)
      SE_neg = 0.d0 ; SE_pos = 0.d0
      num_pos=0 ; num_neg = 0
      DO LP = 1, icpu
         IPTSS = itotal0*(LP - 1) + 1
         IPTFF = MIN(itotal0*LP, itotal)
         DO L = IPTSS,IPTFF
            Vx  = Re(4,L)
            Vy  = Re(5,L)
            Vz  = Re(6,L)
            GAMM = dsqrt(1.d0 + Vx*Vx + Vy*Vy + Vz*Vz)
            sss = GAMM - EKK
            IF(sss.LT.0.d0) THEN
              SE_neg(LP) = SE_neg(LP) + abs(sss)
              num_neg(LP) = num_neg(LP) + 1
            ELSE
              SE_pos(LP) = SE_pos(LP) + abs(sss)
              num_pos(LP) = num_pos(LP) + 1
            END IF
         END DO
      END DO
c
      SIG_pos = 0.d0 ; SIG_neg = 0.d0
      num_p = 0 ; num_n = 0
      DO LP = 1,icpu,4
         SIG_pos = SIG_pos + SE_pos(LP)  + SE_pos(LP + 1)
     &		 + SE_pos(LP + 2) + SE_pos(LP + 3)
         SIG_neg = SIG_neg + SE_neg(LP)  + SE_neg(LP + 1)
     &		 + SE_neg(LP + 2) + SE_neg(LP + 3)
          num_p = num_p   + num_pos(LP) + num_pos(LP + 1)
     &            + num_pos(LP + 2) + num_pos(LP + 3)
          num_n = num_n  + num_neg(LP) + num_neg(LP + 1)
     &            + num_neg(LP + 2) + num_neg(LP + 3)
      END DO
c
      SIGE_pos = 0.d0 ; SIGE_neg = 0.d0
      num_pp = 0 ; num_nn = 0
      CALL mpi_allreduce(SIG_pos,SIGE_pos,1,MPI_REAL8,MPI_SUM
     &                   ,mpi_comm_world,ierr)
      CALL mpi_allreduce(SIG_neg,SIGE_neg,1,MPI_REAL8,MPI_SUM
     &                   ,mpi_comm_world,ierr)
      CALL mpi_allreduce(num_p,num_pp,1,MPI_INTEGER,MPI_SUM
     &                   ,mpi_comm_world,ierr)
      CALL mpi_allreduce(num_n,num_nn,1,MPI_INTEGER,MPI_SUM
     &                   ,mpi_comm_world,ierr)
c
      SIGE_pos = SIGE_pos/num_pp
      SIGE_neg = SIGE_neg/num_nn
c
      IF(myrank.EQ.0) THEN
        WRITE(fo_name2,444) TRIM(data_file)//'AveEne',jobno
        OPEN (20,file=fo_name2,form='formatted',status='unknown')
        WRITE(20,666) t*Rt2, EKK*0.511d6
     &		 , (EKK + SIGE_pos)*0.511d6
     &		 , (EKK - SIGE_neg)*0.511d6
      END IF
444   FORMAT(A,I3.3,'.dat')
666   FORMAT(11(E16.6,1X))
      RETURN
      END
!-------------------------
      SUBROUTINE radiation
!-------------------------
      USE random_common
      USE sim_common
      USE mpi_common
      USE R_common
      USE out_common
      USE omp_lib
      IMPLICIT NONE
      INCLUDE "mpif.h"
      REAL(kind=8),DIMENSION(:),ALLOCATABLE :: total
      REAL(kind=8),DIMENSION(:,:),ALLOCATABLE :: diff1,diff2
c
      IF(myrank.EQ.0) WRITE(*,*) "calculating radiation"
      ALLOCATE(wmit3(icpu),vmit3(icpu))
      ALLOCATE(emit3(6000,icpu)) ; ALLOCATE(fmit3(6000,icpu))
      ALLOCATE(emitT3(6000)) ; ALLOCATE(fmitT3(6000))
      ALLOCATE(emitTT(6000)) ; ALLOCATE(fmitTT(6000))
      ALLOCATE(total(6000))
      ALLOCATE(diff1(0:3000 + 1,200))
      ALLOCATE(diff2(0:3000 + 1,200))
c
      IF(QED.EQ.1) THEN
        total = totalR
        diff1 = diffR
        diff2 = diffQ
      ELSE IF(QED.EQ.0) THEN
        total = totalC
        diff1 = diffD
        diff2 = diffC
      END IF
c
c distribute electron energy reduction onto emission energy spectum
      CALL system_clock(tcurr0)
!$omp workshare
      ENEh = 1.0d0*Em
      ENEd = ENEh/6000.d0
      Lall = int(ksmax/icpu)
      emit3 = 0.d0 ; fmit3 = 0.d0
      wmit3 = 0.d0 ; vmit3 = 0.d0
!$omp end workshare
c
c     diffC ; classical differential cross section (x energy)
c     diffQ ;   quantum differential cross section (x energy)
c     diffD ; classical differential cross section
c     diffR ;   quantum differential cross section
!$omp parallel do private(LP,IPTSS,IPTFF,j,L,Xi
!$omp&  ,ENE,Um,hh,kk,ff,FF1,gg,k,ENE0,Xe,TTT,PXS)
!$omp&shared(we0i,phtn,PKsout,total,wight
!$omp&  ,emit3,fmit3,diffR,diffQ,wmit3,vmit3)
      DO LP = 1,icpu                 ! parallelization loop
         IPTSS = (LP - 1)*Lall + 1
         IPTFF = LP*Lall
      DO j = IPTSS,IPTFF           ! timestep loop
      DO L = 1,itotal              ! particle loop
         Xi = phtn(1,j,L)        ! quantum PARAMETER
         ENE = phtn(2,j,L)        ! energy
         Um = phtn(4,j,L)*Pksout ! energy defference
         hh = phtn(5,j,L)*Pksout ! energy defference
c
         IF(XI.GT.0.001) THEN
           kk = MIN(IDNINT(log(XI*1000.d0)*we0i + 1.5d0),200)
c INTegration of total cross section
           ff = total(kk)*3000.d0
c coefficient
           FF1 = ENE/3000.d0
           gg = Um/ff*ENEd/FF1*wight(L)
c
           DO k = 1,6000
              ENE0 = (DBLE(k) - 0.5d0)*ENEd
              ii = IDNINT(ENE0/ENE*3000.d0 + 0.5d0)
              Xe = ENE0/ENE*3000.d0 + 0.5d0 - DBLE(ii)
              IF(ii.GT.2999) EXIT
              TTT = diff1(ii,kk)+Xe*(diff1(ii + 1,kk) - diff1(ii,kk))
              PXS = diff2(ii,kk)+Xe*(diff2(ii + 1,kk) - diff2(ii,kk))
              emit3(k,LP) = emit3(k,LP) + TTT*gg/ENE
              fmit3(k,LP) = fmit3(k,LP) + PXS*gg
           END DO
        END IF
        wmit3(LP) = wmit3(LP) + Um*wight(L)
        vmit3(LP) = vmit3(LP) + hh*wight(L)
      END DO
      END DO
      END DO
!$omp end parallel do
c
!$omp parallel do private(i,k,Tm,ff)
!$omp&      shared(emit3,fmit3,emitT3,fmitT3)
      DO i = 1,6000
         Tm = 0.d0
         ff = 0.d0
         DO k = 1,icpu,4
            Tm = Tm + emit3(i,k) + emit3(i,k + 1)
     &	     + emit3(i,k + 2) + emit3(i,k + 3)
            ff = ff + fmit3(i,k) + fmit3(i,k + 1)
     &	     + fmit3(i,k + 2) + fmit3(i,k + 3)
         END DO
         emitT3(i) = Tm
         fmitT3(i) = ff
      END DO
!$omp end parallel do
c
      wmitT3 = 0.d0
      vmitT3 = 0.d0
!$omp parallel do private(k) shared(wmit3,vmit3)
!$omp&         reduction(+:wmitT3,vmitT3)
      DO k=1,icpu,4
         wmitT3 = wmitT3 + wmit3(k) + wmit3(k + 1)
     &		+ wmit3(k + 2) + wmit3(k + 3)
         vmitT3 = vmitT3 + vmit3(k) + vmit3(k + 1)
     &	      + vmit3(k + 2) + vmit3(k + 3)
      END DO
!$omp end parallel do
      CALL system_clock(tcurr1)
      tcurr = tcurr + tcurr1 - tcurr0
c
c summation in MPI processes
c
      IF(kstep.NE.ksmax) RETURN
c
      CALL system_clock(trdct0)
      vmitT3  = ABS(vmitT3)
      CALL mpi_allreduce(emitT3,emitTT,6000,mpi_real8
     &                         ,mpi_sum,mpi_comm_world,ierr)
      CALL mpi_allreduce(fmitT3,fmitTT,6000,mpi_real8
     &                         ,mpi_sum,mpi_comm_world,ierr)
      CALL mpi_allreduce(wmitT3,wmitTT,1,mpi_real8
     &                         ,mpi_sum,mpi_comm_world,ierr)
      CALL mpi_allreduce(vmitT3,vmitTT,1,mpi_real8
     &                         ,mpi_sum,mpi_comm_world,ierr)
      CALL system_clock(trdct1)
      trdct = trdct + trdct1 - trdct0
c
c END distribute onto emission energy spectum
c
      WRITE(9,*) "energy reduction(EM work) [J]",
     &		vmitTT*enum*0.511d6*1.6d-19
      WRITE(9,*) "energy reduction (radiation) [J]",
     &		wmitTT*enum*0.511d6*1.6d-19
c
c output photon spectrum data
c
      const = pin/rin*shot*enum
      IF(myrank.EQ.0) THEN
        WRITE(fo_name2,444) TRIM(data_file)//'phtne',jobno
        OPEN (19,file=fo_name2,form='formatted',status='unknown')
c
        DO i = 1,6000
          ENE0 = (DBLE(i)-0.5d0)*ENEd
          WRITE(19,*) ENE0*0.511d6, emitTT(i)*const, fmitTT(i)*const
        END DO
        CLOSE(19)
      END IF
c
      ff = 0.d0
!$omp parallel do private(i) shared(fmit3) reduction(+:ff)
      DO i = 1,6000
         ff = ff + fmitTT(i)*enum
      END DO
!$omp end parallel do
c
      WRITE(9,*) "total reduction energy [J]", ff*0.511d6*1.6d-19

      IF(OutPairs.EQ.0) RETURN

c     setup for theoretical cross sections
      CALL pair_init
c
      IF(myrank.EQ.0) WRITE(*,*) "Convert Pairs"
c
      ALLOCATE(qmit3(6000,30))
      ALLOCATE(qmitT3(6000) )
      ALLOCATE(qmitTT(6000) )
c
      ENEh = 1.d0*Em
      ENEd = ENEh/6000.d0
      Lall = 6000/30
      qmit3 = 0.d0
      DO LP=1,30
         IPTSS = (LP - 1)*Lall + 1
         IPTFF = LP*Lall
      DO i = IPTSS,IPTFF
         ENE0 = (DBLE(i)-0.5d0)*ENEd
         IF(ENE0.GE.2.d0) THEN
           kk   = min(IDNINT(log(ENE0*0.5d0)/dp1 + 1.d0),LPx)
           ff   = totalH(kk)
           ENE  = (1.d0 - dexp(-3.d - 5*ff))*emitT3(i)
c INTegration of total cross section
           ff = totalH2(kk)*1000.d0
c coefficient
           FF1 = ENE0/1000.d0
           gg = ENE/ff*ENEd/FF1
c
           DO k = 1,6000
              ENEA = (DBLE(k)-0.5d0)*ENEd
              ii = IDNINT(ENEA/ENE0*1000.d0 + 0.5d0)
              Xe = ENEA/ENE0*1000.d0 + 0.5d0 - DBLE(ii)
              IF(ii.GE.1000) EXIT
              IF(ii.GE.1) THEN
              TTT = resultL(ii,kk)
     &              + Xe*(resultL(ii + 1,kk) - resultL(ii,kk))
              qmit3(k,LP) = qmit3(k,LP) + TTT*gg
              END IF
           END DO
         END IF
      END DO
      END DO

      DO i = 1,6000
         ff = 0.d0
         DO k = 1,30
            ff = ff + qmit3(i,k)
         END DO
         qmitT3(i)=ff
      END DO
c
c END convert to the positron track data using material filter
c
c summation in MPI processes
c
      CALL mpi_allreduce(qmitT3,qmitTT,6000,mpi_real8
     &                   ,mpi_sum,mpi_comm_world,ierr)
c
c output positron spectrum data
c
      const = pin/rin*shot*enum
      IF(myrank.EQ.0) THEN
        WRITE(fo_name2,444) TRIM(data_file)//'pairTT', jobno
        OPEN (34,file=fo_name2,form='formatted',status='unknown')
c
        DO i = 1,6000
           ENE0 = (DBLE(i)-0.5d0)*ENEd
           WRITE(34,*) ENE0*0.511d6, qmitTT(i)*const
        END DO
      END IF
c
      CLOSE(34)
c

      DEALLOCATE(wmit3,vmit3)
      DEALLOCATE(emit3,fmit3)
      DEALLOCATE(emitT3,fmitT3)
      DEALLOCATE(emitTT,fmitTT)
      DEALLOCATE(total,diff1,diff2)
      DEALLOCATE(qmit3,qmitT3,qmitTT)
444   FORMAT(A,I3.3,'.dat')
      RETURN
      END
!----------------------
      SUBROUTINE angdis
!----------------------
      USE random_common
      USE sim_common
      USE mpi_common
      USE R_common
      USE out_common
      USE omp_lib
      IMPLICIT NONE
      INCLUDE "mpif.h"
      INTEGER,PARAMETER :: ang = 1024
      REAL(kind=8) :: thetay,thetaz,ang2
      REAL(kind=8),DIMENSION(:,:,:),ALLOCATABLE:: emit2
      REAL(kind=8),DIMENSION(:,:),ALLOCATABLE:: emitT2,emitTT2
c
      ang2 = DBLE(ang/2)
      ALLOCATE(emit2(ang,ang,icpu))
      ALLOCATE(emitT2(ang,ang))
      ALLOCATE(emitTT2(ang,ang))
c
      IF(myrank.EQ.0) WRITE(*,*) "Calculating angular distribution..."
c
!$omp workshare
      Lall = int(ksmax/icpu)
      emit2 = 0.d0
!$omp end workshare
c
!$omp parallel do private(LP,IPTSS,IPTFF,j,L,Xi
!$omp&  ,Um,TmY,TmZ,ii,jj)
!$omp&shared(phtn,PKsout,wight,emit2)
      DO LP = 1,icpu
         IPTSS = (LP-1)*Lall + 1
         IPTFF = LP*Lall
      DO j = IPTSS,IPTFF
      DO L = 1,itotal
         Xi  = phtn(1,j,L)        ! quantum PARAMETER
         TmY = phtn(3,j,L)
         TmZ = phtn(6,j,L)
         Um  = phtn(4,j,L)*Pksout ! energy defference
c
         IF(XI.GT.0.001) THEN
           ii = MIN(IDNINT(TmZ/(pi)*ang2 + 0.5d0 + ang2),ang)
           jj = MIN(IDNINT(TmY/(pi)*ang2 + 0.5d0 + ang2),ang)  ! angle in radian
           emit2(ii,jj,LP) = emit2(ii,jj,LP) + Um*wight(L)
         END IF
      END DO
      END DO
      END DO
!$omp end parallel do

!$omp parallel do private(j,i,Tm,k) shared(emit2,emitT2)
	DO j = 1,ang
	DO i = 1,ang
         Tm = 0.d0
         DO k = 1,icpu,4
            Tm = Tm + emit2(i,j,k ) + emit2(i,j,k + 1)
     &		  + emit2(i,j,k + 2) + emit2(i,j,k + 3)
         END DO
         emitT2(i,j) = Tm
      END DO
	END DO
!$omp end parallel do
c
c summation in MPI processes
c
      IF(kstep.NE.ksmax) RETURN
c
      CALL mpi_allreduce(emitT2,emitTT2,ang*ang,mpi_real8
     &                   ,mpi_sum,mpi_comm_world,ierr)
c END distribute onto emission energy spectum

c output photon spectrum data
      const = pin/rin*shot*enum

      IF(myrank.EQ.0) THEN
        WRITE(fo_name2,444) TRIM(data_file)//'phtnTe', jobno
        OPEN (21,file=fo_name2,form='formatted',status='unknown')
c
        DO j = 1,ang
        DO i = 1,ang
           thetaz=(i - 0.5 - INT(ang2))/INT(ang2)*(pi)
           thetay=(j - 0.5 - INT(ang2))/INT(ang2)*(pi)
           WRITE(21,666) thetaz, thetay, emitTT2(i,j)
        END DO
        END DO
        CLOSE(21)
      END IF
c
      DEALLOCATE(emit2,emitT2,emitTT2)
c
444   FORMAT(A,I3.3,'.dat')
666   FORMAT(3(E14.4,1X))
      RETURN
      END
!-----------------------
      SUBROUTINE lorentz
!-----------------------
      USE random_common
      USE sim_common
      USE mpi_common
      USE R_common
      USE omp_lib
      IMPLICIT NONE
c
      itotal0 = INT(itotal/icpu + 1)
!$omp do
      DO LP = 1,icpu
         IPTSS = itotal0*(LP - 1) + 1
         IPTFF = MIN(itotal0*LP,itotal)
         DO i = IPTSS,IPTFF      ! particle loop
         DT1 = PQM*DT*0.5d0
         Xe = Re(1,i)
         Ye = Re(2,i)
         Ze = Re(3,i)
         Vx0 = Re(4,i)
         Vy0 = Re(5,i)
         Vz0 = Re(6,i)
         ENE00 = dsqrt(1.d0 + Vx0**2 + Vy0**2 + Vz0**2)
c
         IF(shape.EQ.1)THEN
            CALL gaussian
         ELSE
            CALL parax
         END IF
c
         VX = VX0 + AVEX*DT1
         VY = VY0 + AVEY*DT1
         VZ = VZ0 + AVEZ*DT1
         gg = DT1/SQRT(1.0 + VX*VX + VY*VY + VZ*VZ)
c
         BXT = AVBX*gg
         BYT = AVBY*gg
         BZT = AVBZ*gg
c
         PXS = VX + (VY*BZT - VZ*BYT)
         PYS = VY + (VZ*BXT - VX*BZT)
         PZS = VZ + (VX*BYT - VY*BXT)
c
         ff = 2.0/(1.0+(BXT*BXT + BYT*BYT + BZT*BZT))
         VX = VX + AVEX*DT1 + ff*(PYS*BZT - PZS*BYT)
         VY = VY + AVEY*DT1 + ff*(PZS*BXT - PXS*BZT)
         VZ = VZ + AVEZ*DT1 + ff*(PXS*BYT - PYS*BXT)
c
         ENE0 = dsqrt(1.d0 + VX*VX  + VY*VY + VZ*VZ)
         GAMMI = 1.d0/ENE0
         Wx0 = Vx*GAMMI
         Wy0 = Vy*GAMMI
         Wz0 = Vz*GAMMI
c
         DTI = 1.0d0/DT
         BXT = (Vx - Vx0)*DTI ! fVx
         BYT = (Vy - Vy0)*DTI ! fVy
         BZT = (Vz - Vz0)*DTI ! fVz
c
         Vx0 = Vx
         Vy0 = Vy
         Vz0 = Vz
c
c------------------------------------------------------------------
c
         ff = BXT*BXT + BYT*BYT + BZT*BZT
         gg = BXT*Wx0 + BYT*Wy0 + BZT*Wz0
c
         XI = Alf2*ENE0*sqrt(abs(ff - gg**2))*0.66667d0
         PXS = (BXT - Wx0*gg)*Alf
         PYS = (BYT - Wy0*gg)*Alf
         PZS = (BZT - Wz0*gg)*Alf
         RRR = ENE0**2*(BXT*PXS + BYT*PYS + BZT*PZS)*dt
         RAD = RAD + RRR
         TTT = AVEX*(Vx*GAMMI)*dt
     &         + AVEY*(Vy*GAMMI)*dt
     &         + AVEZ*(Vz*GAMMI)*dt
         Re(1,i) = Xe + Wx0*dt
         Re(2,i) = Ye + Wy0*dt
         Re(3,i) = Ze + Wz0*dt
         Re(4,i) = Vx
         Re(5,i) = Vy
         Re(6,i) = Vz
         Re(7,i) = TTT
         Re(8,i) = (ENE00 - ENE)-TTT
         Re(9,i) = RRR !ENE0 - ENE - TTT
c
         Re(10,i) = Xi
         Re(11,i) = ENE0
         END DO
      END DO     ! END particle loop
!$omp end do
      RETURN
      END
!-----------------------
      SUBROUTINE Sokolov
!-----------------------
      USE random_common
      USE sim_common
      USE mpi_common
      USE R_common
      USE omp_lib
      IMPLICIT NONE
c
      itotal0 = INT(itotal/icpu + 1)
!$omp do
      DO LP = 1,icpu
         IPTSS = itotal0*(LP - 1) + 1
         IPTFF = MIN(itotal0*LP,itotal)
         DO i = IPTSS,IPTFF      ! particle loop
         DT1 = PQM*DT*0.5d0
         Xe = Re(1,i)
         Ye = Re(2,i)
         Ze = Re(3,i)
         Vx0 = Re(4,i)
         Vy0 = Re(5,i)
         Vz0 = Re(6,i)
         ENE00 = dsqrt(1.d0 + Vx0**2 + Vy0**2 + Vz0**2)
c
         IF(shape.EQ.1) THEN
           CALL gaussian
         ELSE
           CALL parax
         END IF
c
         VX = VX0 + AVEX*DT1
         VY = VY0 + AVEY*DT1
         VZ = VZ0 + AVEZ*DT1
         gg = DT1/SQRT(1.0 + VX*VX + VY*VY + VZ*VZ)
c
         BXT = AVBX*gg
         BYT = AVBY*gg
         BZT = AVBZ*gg
c
         PXS = VX + (VY*BZT - VZ*BYT)
         PYS = VY + (VZ*BXT - VX*BZT)
         PZS = VZ + (VX*BYT - VY*BXT)
c
         ff  = 2.0/(1.0 + (BXT*BXT + BYT*BYT + BZT*BZT))
         VX  = VX + AVEX*DT1 + ff*(PYS*BZT - PZS*BYT)
         VY  = VY + AVEY*DT1 + ff*(PZS*BXT - PXS*BZT)
         VZ  = VZ + AVEZ*DT1 + ff*(PXS*BYT - PYS*BXT)
c
         ENE0 = dsqrt(1.d0 + VX*VX  + VY*VY + VZ*VZ)
         GAMMI = 1.d0/ENE0
         Wx0 = Vx*GAMMI
         Wy0 = Vy*GAMMI
         Wz0 = Vz*GAMMI
c
         DTI = 1.0d0/DT
         BXT = (Vx - Vx0)*DTI ! fVx
         BYT = (Vy - Vy0)*DTI ! fVy
         BZT = (Vz - Vz0)*DTI ! fVz
c
         Vx0 = Vx
         Vy0 = Vy
         Vz0 = Vz
c------------------------------------------------------------------
         ff = BXT*BXT + BYT*BYT + BZT*BZT
         gg = BXT*Wx0 + BYT*Wy0 + BZT*Wz0
c
         XI = Alf2*ENE0*sqrt(abs(ff - gg**2))*0.66667d0
         IF(XI.GT.0.001d0) THEN
            kk = MIN(IDNINT(log(XI*1000.d0)*we0i + 1.5d0),200)
            ENEh = Alf*totalRC(kk)
         ELSE
            ENEh = Alf
         END IF
c
         GAMMI = ENEh/(1.0 + ENEh*gg)
         PXS = (BXT - Wx0*gg)*GAMMI
         PYS = (BYT - Wy0*gg)*GAMMI
         PZS = (BZT - Wz0*gg)*GAMMI
c
         GAMMI = (BXT*PXS + BYT*PYS + BZT*PZS)*(ENE0**2)
         Vx = Vx0 - ((PYS*AVBZ - PZS*AVBY) + Wx0*GAMMI)*dt
         Vy = Vy0 - ((PZS*AVBX - PXS*AVBZ) + Wy0*GAMMI)*dt
         Vz = Vz0 - ((PXS*AVBY - PYS*AVBX) + Wz0*GAMMI)*dt
c
         ENE = SQRT(1.d0 + VX*VX + VY*VY + VZ*VZ)
c
         GAMMI = 1.d0/ENE
c
         TTT = AVEX*(Vx*GAMMI + PXS)*dt
     &         + AVEY*(Vy*GAMMI + PYS)*dt
     &         + AVEZ*(Vz*GAMMI + PZS)*dt
         RRR = ENE0**2*(BXT*PXS + BYT*PYS + BZT*PZS)*dt
         TT1 = TT1 + TTT
         RAD = RAD + RRR
         Re(1,i) = Xe + Vx*GAMMI*dt + PXS*dt
         Re(2,i) = Ye + Vy*GAMMI*dt + PYS*dt
         Re(3,i) = Ze + Vz*GAMMI*dt + PZS*dt
         Re(4,i) = Vx
         Re(5,i) = Vy
         Re(6,i) = Vz
         Re(7,i) = TTT
         Re(8,i) = (ENE00 - ENE) - TTT
         Re(9,i) = RRR
c
         Re(10,i) = Xi
         Re(11,i) = ENE0
         END DO
      END DO
!$omp end do
      RETURN
      END
!--------------------------------------
      SUBROUTINE Landau_Lifshitz
!--------------------------------------
      USE random_common
      USE sim_common
      USE mpi_common
      USE R_common
      USE omp_lib
      IMPLICIT NONE
      REAL(kind=8) :: LLT2X,LLT2Y,LLT2Z,LL3TX,LL3TY,LL3TZ
      REAL(kind=8) :: BXK,BYK,BZK
c
      itotal0 = INT(itotal/icpu + 1)
!$omp do
      DO LP = 1,icpu
         IPTSS = itotal0*(LP - 1) + 1
         IPTFF = MIN(itotal0*LP,itotal)
         DO i = IPTSS,IPTFF      ! particle loop
         DT1 = PQM*DT*0.5d0
         Xe = Re(1,i)
         Ye = Re(2,i)
         Ze = Re(3,i)
         Vx0 = Re(4,i)
         Vy0 = Re(5,i)
         Vz0 = Re(6,i)
         ENE00 = dsqrt(1.d0 + Vx0**2 + Vy0**2 + Vz0**2)
c
         IF(shape.EQ.0)THEN
           CALL gaussian
         ELSE
           CALL parax
         END IF
c
         VX = VX0 + AVEX*DT1
         VY = VY0 + AVEY*DT1
         VZ = VZ0 + AVEZ*DT1
         gg = DT1/SQRT(1.0 + VX*VX + VY*VY + VZ*VZ)
c
         BXT = AVBX*gg
         BYT = AVBY*gg
         BZT = AVBZ*gg
c
         PXS = VX +(VY*BZT - VZ*BYT)
         PYS = VY +(VZ*BXT - VX*BZT)
         PZS = VZ +(VX*BYT - VY*BXT)
c
         ff  = 2.0/(1.0+(BXT*BXT + BYT*BYT + BZT*BZT))
         VX  = VX + AVEX*DT1 + ff*(PYS*BZT - PZS*BYT)
         VY  = VY + AVEY*DT1 + ff*(PZS*BXT - PXS*BZT)
         VZ  = VZ + AVEZ*DT1 + ff*(PXS*BYT - PYS*BXT)
c
         ENE0 = d sqrt(1.d0 + VX*VX  + VY*VY + VZ*VZ)
         GAMMI = 1.d0/ENE0
         Wx0 = Vx*GAMMI
         Wy0 = Vy*GAMMI
         Wz0 = Vz*GAMMI
c
         DTI = 1.0d0/DT
         BXT = (Vx - Vx0)*DTI ! fVx
         BYT = (Vy - Vy0)*DTI ! fVy
         BZT = (Vz - Vz0)*DTI ! fVz
c
         Vx0 = Vx
         Vy0 = Vy
         Vz0 = Vz
c------------------------------------------------------------------
         ff = BXT*BXT + BYT*BYT + BZT*BZT
         gg = BXT*Wx0 + BYT*Wy0 + BZT*Wz0
c
         XI = Alf2*ENE0*sqrt(abs(ff - gg**2))*0.66667d0
         IF(XI.GT.0.001d0) THEN
            kk    = MIN(IDNINT(log(XI*1000.d0)*we0i + 1.5d0),200)
            ENEh  = Alf*totalRC(kk)
         ELSE
            ENEh  = Alf
         END IF
         GAMMI = ENEh
         PXS   = (BXT - Wx0*gg)*GAMMI
         PYS   = (BYT - Wy0*gg)*GAMMI
         PZS   = (BZT - Wz0*gg)*GAMMI
c
         BXK = BYT*AVBZ - BZT*AVBY
         BYK = - (BXT*AVBZ - BZT*AVBX)
         BZK = BXT*AVBY - BYT*AVBX
c
         LLT2X = ENEh*(BXK + gg*AVEX)
         LLT2Y = ENEh*(BYK + gg*AVEY)
         LLT2Z = ENEh*(BZK + gg*AVEZ)
c
         GAMMI = (BXT*PXS + BYT*PYS + BZT*PZS)*(ENE0**2)
         LL3TX = Wx0*GAMMI
         LL3TY = Wy0*GAMMI
         LL3TZ = Wz0*GAMMI
c
         Vx = Vx0 + dt*(LLT2X - LL3TX)
         Vy = Vy0 + dt*(LLT2Y - LL3TY)
         Vz = Vz0 + dt*(LLT2Z - LL3TZ)
c
         ENE = SQRT(1.d0 + VX*VX + VY*VY + VZ*VZ)
         GAMMI = 1.d0/ENE
         Wx = Vx*GAMMI
         Wy = Vy*GAMMI
         Wz = Vz*GAMMI
         TTT = AVEX*Wx + AVEY*Wy + AVEZ*Wz
         RRR = ENE0**2*(BXT*PXS + BYT*PYS + BZT*PZS)*dt
         TT1 = TT1 + TTT
         RAD = RAD + RRR
         Re(1,i) = Xe + Wx*dt
         Re(2,i) = Ye + Wy*dt
         Re(3,i) = Ze + Wz*dt
         Re(4,i) = Vx
         Re(5,i) = Vy
         Re(6,i) = Vz
         Re(7,i) = TTT
         Re(8,i) = (ENE00 - ENE) - TTT
         Re(9,i) = RRR
         Re(10,i) = Xi
         Re(11,i) = ENE0
      END DO
      END DO
!$omp end do
      RETURN
      END
!-------------------------------------
      SUBROUTINE random
!-------------------------------------
      USE random_common
      USE sim_common
      IMPLICIT NONE
      INTEGER, PARAMETER :: values = 245
      INTEGER :: d
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
c
c-----Dynamic RanDOm Number Generator------
c      CALL date_and_time(values=values)
      CALL random_seed(size = d)
      ALLOCATE(seed(1:d))
      seed = values
      CALL random_seed(put=seed)
c
      RETURN
      END
!-------------------------------------
      SUBROUTINE phemit
!-------------------------------------
      USE random_common
      USE sim_common
      IMPLICIT NONE
      REAL(kind=8) :: det1,det2,redcomp,ss
c
      redcomp = comp/(2.d0*PI)
      det1 = 1.d0/(PI*dsqrt(3.d0))
      det2 = 1.d0/(137.d0*redcomp) ! probability coefficient
c
      ff = totalS(kk)*3000
      ss = totalS(kk)*det1*det2*Rt*VL*1.d0/ENE0 !
c
      IF(ss.GE.1.d - 3) THEN
        WRITE(*,*) " W*dt > 1 --- Probability of emission > 1"
        WRITE(*,*) " Please reduce time step, "
        WRITE(*,*) " i.e. increase div in SLAC.dat"
        photon = -1
        stop
      END IF
c      WRITE(*,*) ss
c
      CALL random_number(rand)
c      WRITE(*,*) rand
      IF(rand.LT.ss) THEN
	  emmits = .TRUE.
      ELSE
	  emmits = .FALSE.
      END IF
c
      RETURN
      END
!-------------------------------------
      SUBROUTINE qmemit
!-------------------------------------
      USE random_common
      USE sim_common
      IMPLICIT NONE
      REAL(kind=8),DIMENSION(6000) :: W
      REAL(kind=8) :: gg1,ss
c
      CALL random_number(rand)
c      WRITE(*,*) rand
c
      DO k = 1,6000
         ii = IDNINT((k - 0.5d0)*0.5d0 + 0.5d0)
         IF(ii.GT.2999) exit
         ss = (k - 0.5d0)*0.5d0 + 0.5d0 - DBLE(ii)
         TTT= diffR(ii,kk) + ss*(diffR(ii + 1,kk) - diffR(ii,kk))
         W(k + 1) = W(k) + TTT*0.5d0/ff
         gg1 = W(k + 1)
         IF(gg1.GE.rand) THEN
           ENN = k + (rand - W(k))/(W(k + 1)-W(k))
         exit
         END IF
      END DO
      RETURN
      END
!-------------------------------------
      SUBROUTINE qedemmit
!-------------------------------------
      USE sim_common
      USE mpi_common
      USE R_common
      USE omp_lib
      IMPLICIT NONE
c
      IF(QED.EQ.0) THEN
        WRITE(*,*) "QED must be turned on"
        WRITE(*,*) "set QED = 1 in SLAC.dat"
        stop
	END IF
c
      itotal0 = INT(itotal/icpu + 1)
!$omp do
      DO LP = 1,icpu
         IPTSS = itotal0*(LP - 1)+1
         IPTFF = MIN(itotal0*LP,itotal)
         DO i = IPTSS,IPTFF      ! particle loop
         DT1 = PQM*DT*0.5d0
         Xe = Re(1,i)
         Ye = Re(2,i)
         Ze = Re(3,i)
         Vx0 = Re(4,i)
         Vy0 = Re(5,i)
         Vz0 = Re(6,i)
         ENE00 = dsqrt(1.d0 + Vx0**2 + Vy0**2 + Vz0**2)
c
         IF(shape.EQ.1)THEN
           CALL gaussian
         ELSE
           CALL parax
         END IF
c
         VX = VX0 + AVEX*DT1
         VY = VY0 + AVEY*DT1
         VZ = VZ0 + AVEZ*DT1
         gg = DT1/SQRT(1.0 + VX*VX + VY*VY + VZ*VZ)
c
         BXT = AVBX*gg
         BYT = AVBY*gg
         BZT = AVBZ*gg
c
         PXS = VX +(VY*BZT - VZ*BYT)
         PYS = VY +(VZ*BXT - VX*BZT)
         PZS = VZ +(VX*BYT - VY*BXT)
c
         ff = 2.0/(1.0+(BXT*BXT + BYT*BYT + BZT*BZT))
         VX = VX + AVEX*DT1 + ff*(PYS*BZT - PZS*BYT)
         VY = VY + AVEY*DT1 + ff*(PZS*BXT - PXS*BZT)
         VZ = VZ + AVEZ*DT1 + ff*(PXS*BYT - PYS*BXT)
c
         ENE0 = dsqrt(1.d0 + VX*VX  + VY*VY + VZ*VZ)
         GAMMI = 1.d0/ENE0
         Wx0 = Vx*GAMMI
         Wy0 = Vy*GAMMI
         Wz0 = Vz*GAMMI
c
	   DTI = 1.0d0/DT
         BXT = (Vx - Vx0)*DTI ! fVx
         BYT = (Vy - Vy0)*DTI ! fVy
         BZT = (Vz - Vz0)*DTI ! fVz
c
         Vx0 = Vx
         Vy0 = Vy
         Vz0 = Vz
c------------------------------------------------------------------
         ff = BXT*BXT + BYT*BYT + BZT*BZT
         gg = BXT*Wx0 + BYT*Wy0 + BZT*Wz0
c
         XI = Alf2*ENE0*sqrt(abs(ff - gg**2))*0.66667d0
         IF(XI.GT.0.001d0) THEN
            kk = MIN(IDNINT(log(XI*1000.d0)*we0i + 1.5d0),200)
            ENEh = Alf*totalRC(kk)
c        calculates whether emission should occur or not
            CALL phemit
         IF(emmits) THEN
            CALL qmemit
            ENN = ENN/6000
            WRITE(*,*) "Photon Energy [MeV] =", ENN*ENE0*0.511
c  ******* Update electron momentum due to recoil *******
            Vx = Vx0*(1.d0 - ENN)
            Vy = Vy0*(1.d0 - ENN)
            Vz = Vz0*(1.d0 - ENN)
            ELSE
            ENEh = Alf
            Vx = Vx0
            Vy = Vy0
            Vz = Vz0
            END IF
         ELSE
            ENEh = Alf
            Vx = Vx0
            Vy = Vy0
            Vz = Vz0
         END IF
c
         GAMMI = ENEh
         PXS = (BXT - Wx0*gg)*GAMMI
         PYS = (BYT - Wy0*gg)*GAMMI
         PZS = (BZT - Wz0*gg)*GAMMI
         ENE = SQRT(1.d0 + VX*VX + VY*VY + VZ*VZ)
         GAMMI = 1.d0/ENE
c
         TTT = AVEX*(Vx*GAMMI)*dt
     &          + AVEY*(Vy*GAMMI)*dt
     &          + AVEZ*(Vz*GAMMI)*dt
         RRR = ENE0**2*(BXT*PXS + BYT*PYS + BZT*PZS)*dt
         TT1 = TT1 + TTT
         RAD = RAD + RRR
         Re(1,i) = Xe + Vx*GAMMI*dt
         Re(2,i) = Ye + Vy*GAMMI*dt
         Re(3,i) = Ze + Vz*GAMMI*dt
         Re(4,i) = Vx
         Re(5,i) = Vy
         Re(6,i) = Vz
         Re(7,i) = TTT
         Re(8,i) = (ENE00 - ENE) - TTT
         Re(9,i) = ENN*ENE0
         Re(10,i) = Xi
         Re(11,i) = ENE0
         END DO
      END DO
!$omp end do
      RETURN
      END
!-------------------------------------
      SUBROUTINE histogram
!-------------------------------------
      USE sim_common
      USE mpi_common
      USE R_common
      USE out_common
      USE omp_lib
      IMPLICIT NONE
      INCLUDE "mpif.h"
      INTEGER :: enegrid
      REAL(kind = 8) :: b4,b5,b6,enediv,aa,sigma,ffsum
      REAL(kind=8),DIMENSION(:), ALLOCATABLE:: his,his_sum
c
      enegrid = 512
      ALLOCATE(his(enegrid))
      ALLOCATE(his_sum(enegrid))

      sigma = MAX(sigmax,sigmay,sigmaz)
!$omp workshare
      enediv  = 1/(Em*(1 + 5.d0*sigma))
      his = 0.d0
      his_sum = 0.d0
!$omp end workshare
c
      itotal0 = INT(itotal/icpu) + 1
!$omp do
      DO LP = 1,icpu
         IPTSS = itotal0*(LP - 1) + 1
         IPTFF = MIN(itotal0*LP,itotal)
      DO k = IPTSS,IPTFF      ! particle loop
         b4 = Re(4,k) ; b5 = Re(5,k) ; b6 = Re(6,k)
         ENE = dsqrt(1.d0 + b4**2 + b5**2 + b6**2) - 1.d0
         i = max(IDNINT(ENE*enediv*enegrid),1)
         i = MIN(enegrid,i)
         his(i) = his(i) + wight(k)
      END DO
      END DO
!$omp end do
c
      ff = 0.d0
      DO i = 1,enegrid
         ff = ff + his(i)
      END DO

      CALL mpi_allreduce(his,his_sum,enegrid,mpi_real8
     &                   ,mpi_sum,mpi_comm_world,ierr)
      CALL mpi_allreduce(ff,ffsum,1,mpi_real8
     &                   ,mpi_sum,mpi_comm_world,ierr)
c
      const = (bin/0.511d6)/(Em*(1 + 5.d0*sigma)/enegrid)*enum
      IF(myrank.EQ.0) THEN
c	  WRITE(*,*) const
        WRITE(fo_name2,444) TRIM(data_file)//'dist_fn'
     &			   ,kstep/100000,jobno
        OPEN (31,file = fo_name2,form='formatted',status='unknown')
        DO i = 1,enegrid
           WRITE(31,666) i/(enegrid*enediv)*0.511d6
     &		      , his_sum(i)
        END DO
        CLOSE(31)
      END IF
c
      IF(kstep.EQ.0) WRITE(9,*) 'Total init. particle', ffsum
      IF(kstep.EQ.ksmax) WRITE(9,*) 'Total final. particle', ffsum
      DEALLOCATE(his,his_sum)

444   FORMAT(A,I4.4,I4.4,'.dat')
666   FORMAT(2(E14.4,1X))
      RETURN
      END
!-------------------------------------
      SUBROUTINE photon_his
!-------------------------------------
      USE sim_common
      USE mpi_common
      USE R_common
      USE out_common
      USE omp_lib
      IMPLICIT NONE
      INCLUDE "mpif.h"
      INTEGER :: enegrid
      REAL(kind=8) :: b4,b5,b6,enediv,aa,sigma,ffsum
      REAL(kind=8),DIMENSION(:), ALLOCATABLE:: his,his_sum
c
      enegrid = 512
      ALLOCATE(his(enegrid))
      ALLOCATE(his_sum(enegrid))

!$omp workshare
      enediv  = 1/Em
      Lall = int(ksmax/icpu)
      his = 0.d0
      his_sum = 0.d0
!$omp end workshare
c
!$omp do
      DO LP = 1,icpu                 ! parallelization loop
         IPTSS = (LP-1)*Lall + 1
         IPTFF = LP*Lall
      DO j=IPTSS,IPTFF           ! timestep loop
      DO L = 1,itotal              ! particle loop
         ENE = phtn(4,j,L)*Pksout*wight(L)
         i = MAX(IDNINT(ENE*enediv*enegrid),1)
         i = MIN(enegrid,i)
         his(i) = his(i) + 1
      END DO
      END DO
      END DO
!$omp end do
c
      ff = 0.d0
      DO i = 1,enegrid
         ff = ff + his(i)
      END DO

      CALL mpi_allreduce(his,his_sum,enegrid,mpi_real8
     &                   ,mpi_sum,mpi_comm_world,ierr)
      CALL mpi_allreduce(ff,ffsum,1,mpi_real8
     &                   ,mpi_sum,mpi_comm_world,ierr)
      const = pin/(Em/enegrid)*enum
      IF(myrank.EQ.0) THEN
        WRITE(*,*) const
        WRITE(fo_name2,444) TRIM(data_file)//'dist_ph',jobno
        OPEN (34,file = fo_name2,form='formatted',status='unknown')
        DO i=1,enegrid
           ENE0 = i*Em/enegrid
           WRITE(34,666) ENE0*0.511d6, his_sum(i)*const
        END DO
        CLOSE(34)
      END IF
c
      IF(kstep.EQ.ksmax) WRITE(9,*) 'Total photon', ffsum
      DEALLOCATE(his,his_sum)

444   FORMAT(A,I4.4,'.dat')
666   FORMAT(2(E14.4,1X))
      RETURN
      END
!-------------------------------------
      SUBROUTINE histogram2d
!-------------------------------------
      USE sim_common
      USE mpi_common
      USE R_common
      USE out_common
      USE omp_lib
      IMPLICIT NONE
      INCLUDE "mpif.h"
      INTEGER :: grid1, grid2
      REAL(kind=8) :: b1,b4,b5,b6,momdiv1,momdiv2,aa
      REAL(kind=8),DIMENSION(:,:), ALLOCATABLE:: his,his_sum
c
      grid1 = 64
      grid2 = 64
      ALLOCATE(his(grid1,grid2))
      ALLOCATE(his_sum(grid1,grid2))

      Vx = dsqrt((Em + 1.d0)**2-1.d0)
!$omp workshare
      momdiv1 = 1.d0/(0.15*Vx)
      momdiv2 = momdiv1
      his = 0.d0
!$omp end workshare
c
      itotal0 = INT(itotal/icpu + 1)
!$omp do
      DO LP = 1,icpu
         IPTSS = itotal0*(LP - 1)+1
         IPTFF = MIN(itotal0*LP,itotal)
      DO k = IPTSS,IPTFF      ! particle loop
         b4 = Re(4,k) ; b5=Re(5,k) ; b6 = Re(6,k)
         b5 = -1.d0*b4*dsin(inc_ang) + b5*dcos(inc_ang)
         i = max(IDNINT(b5*momdiv1*32.d0)+32,1)
         j = max(IDNINT(b6*momdiv2*32.d0)+32,1)
         i = MIN(grid1,i)
         j = MIN(grid2,j)
         his(i,j) = his(i,j) + wight(k)
      END DO
      END DO
!$omp end do
c
      his_sum = 0.d0
      CALL mpi_allreduce(his,his_sum,grid1*grid2,mpi_real8
     &                   ,mpi_sum,mpi_comm_world,ierr)
c
      IF(myrank.EQ.0) THEN
        WRITE(fo_name2,444) TRIM(data_file)//'dist_fn2d'
     &			   ,kstep/100000,jobno
        OPEN (32,file = fo_name2,form='formatted',status='unknown')
c
        DO i = 1,grid1
        DO j = 1,grid2
           WRITE(32,666) (i - 32)/(32*momdiv1)   !py
     &			,(j-32)/(32*momdiv2)   !pz
     &			,his_sum(i,j)
        END DO
        END DO

        CLOSE(32)
	END IF
c
      DEALLOCATE(his,his_sum)
444   FORMAT(A,I4.4,I4.4,'.dat')
666   FORMAT(3(E14.4,1X))
      RETURN
      END
!--------------------------------------
      SUBROUTINE manual_load
!--------------------------------------
      USE random_common
      USE sim_common
      IMPLICIT NONE
      REAL(kind=8),DIMENSION(6000) :: W
      REAL(kind=8) :: gg1,ss,energy
      INTEGER, PARAMETER :: np_local = 116
      INTEGER :: ip
      REAL(kind=8), DIMENSION(np_local) :: Ex_axis, dist_fn

      OPEN(33,file='f_E_smoothed_new_final.txt',status='OLD')

      DO ip = 1,np_local

        READ(33,*) Ex_axis(ip), dist_fn(ip)

        Ex_axis(ip) = Ex_axis(ip)/0.511
        dist_fn(ip) = dist_fn(ip)*5

      ENDDO

      CLOSE(33)
c
      CALL random_number(rand)
c
      W = 0.d0
      DO k = 1,3000
         ii = IDNINT((k - 0.5d0)*116.d0/3000.d0 + 0.5d0)
         IF(ii.GT.115) exit
         ss = (k - 0.5d0)*116.d0/3000.d0 + 0.5d0 - DBLE(ii)
         TTT = dist_fn(ii) + ss*(dist_fn(ii + 1) - dist_fn(ii))
         W(k + 1) = W(k) + TTT*116.d0/3000.d0
         gg1 = W(k + 1)
         IF(rand.LT.gg1) THEN
           ENN = k + (rand-W(k))/(W(k + 1)-W(k))
         EXIT
         END IF
      END DO

      energy = ENN*2300/0.511/3000
      p_x = dsqrt((energy + 1.d0)**2 - 1.d0)
      Vx  = p_x
      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE pair_init
c cross section of Bremsstrahlung & pair production (Bethe-Heitler)
C------------------------------------------------------------------
      USE random_common
      USE mpi_common
      USE out_common
      IMPLICIT NONE
      INTEGER :: i,j
      REAL(kind=8) :: pi,pi2
      REAL(kind=8) :: totalmin,totalmin2,corrct
      REAL(kind=8) :: kx,k0,dk,kk,k,E0,p0,E,p,eL,ec,ee,ppi,p3i,p3j
      REAL(kind=8) :: D1,D2,D3,D4,D5,Fd,Fd0,Fd1,v0,v1,cor0,cor1,cor2
      REAL(kind=8) :: ele,red,totL
c
      REAL(kind=8),DIMENSION(LE0,LPx) :: resultH
      REAL(kind=8),DIMENSION(LE0,LPx) :: resultI
      REAL(kind=8),DIMENSION(LE0,LPx) :: resultJ
      REAL(kind=8),DIMENSION(LE0,LPx) :: resultK
c
      PI  = 4.D0*DATAN(1.D0)
      PI2 = PI*2.D0
c
      IF(myrank.EQ.0) WRITE(*,*) "Bethe - Heitler"
c
      totalmin = log(183.d0/Zcm3)
      totalmin2 = 0.925d0*(Zcom/137.d0)**2 ! for high Z
c     totalmin2 = 1.200d0*(Zcom/137.d0)**2 ! for  low Z
      dp1 = (log(Emax)-log(2.0d0))/BPx  !!
      corrct = 0.6d0/82.d0*Zcom

      DO j = 1,LPx
         kx = dexp(dp1*DBLE(j-1))*2.d0  ! initial photonn energy in x
         k0 = dsqrt(kx**2)
         dk =(k0 - 2.d0)/DBLE(LE0)       ! minus 2*mc^2
         kk = 1.d0/(k0*k0*k0)
         IF(k0.LE.2.d0) CYCLE
         DO i = 1,LE0
            k = dk*(DBLE(i)-0.5d0) ! ratio of emission positron energy
            E0 = k + 1.d0     ! positron energy
            p0 = sqrt(E0**2-1.d0) ! positron momentum

            E = k0 - E0
            p = sqrt(E**2-1.d0)

            eL = 2.d0*log((E*E0 + p*p0 + 1.d0)/k0) ! L
            ec = 2.d0*log(E0 + P0)               ! E+
            ee = 2.d0*log(E + p)               ! E-

            ppi = 1.d0/(p0*p)
            p3i = 1.d0/(p0**3)
            p3j = 1.d0/(p**3)

            D1 = -4.d0/3.d0 - 2.d0*E*E0*(p**2 + p0**2)*ppi**2
            D2 = ec*E*p3i + ee*E0*p3j - ee*ec*ppi
            D3 = (8.d0/3.d0*E*E0*ppi)
            D4 = (k0**2*(ppi**3))*((E0*E)**2 + (p0*p)**2)
            D5 = (0.5d0*k0*ppi)*((E0*E - p0**2)*ec*p3i
     &           +(E0*E - p**2)*ee*p3j + 2.d0*k0*E0*E*(ppi**2))
            Fd = (p0*p*kk*(D1 + D2 + eL*(D4 - D3 - D5)))*(k0 - 2.d0)

            resultH(i,j) = Fd
c
c     relativistic & screening
c
            eL = (E**2 + E0**2 + (2.d0/3.d0)*E*E0)*totalmin
            ec = E*E0/9.d0
            Fd0 = 4.d0*kk*(eL - ec)*(k0 - 2.d0)
c
            resultI(i,j) = Fd0
c
c     Coulomb correction for relativistic case
c
            eL = (E**2 + E0**2 + (2.d0/3.d0)*E*E0)*totalmin2
            Fd1 = -2.d0*kk*eL*(k0 - 2.d0)

            resultJ(i,j) = Fd1 + Fd0
c
c     Coulomb correction for non relativistic case
c
            v0 = p0/dsqrt(1.d0 + p0**2)
            v1 = p/dsqrt(1.d0 + p**2)
            cor0 = corrct/v0
            cor1 = corrct/v1
            cor2 = pi2*pi2*cor0*cor1/(dexp(pi2*cor0) - 1.d0)
     &                          /(1.d0 - dexp(-pi2*cor1))
            resultK(i,j) = Fd*cor2
        END DO
      END DO
c
      resultL = 0.d0
      DO j = 1,LPx
         kx = dexp(dp1*DBLE(j - 1))*2.d0  ! initial photonn energy in x
         k0 = dsqrt(kx**2)
         dk = (k0 - 2.d0)/DBLE(LE0)       ! minus 2*mc^2
         kk = 1.d0/(k0*k0*k0)
c
         DO i=1,LE0
c
            k = dk*(DBLE(i)-0.5d0) ! ratio of emission positron energy
            E0 = k + 1.d0     ! positron energy
            p0 = sqrt(E0**2 - 1.d0) ! positron momentum

            E = k0 - E0
            p = sqrt(E**2 - 1.d0)

            eL = MIN(E,E0)
            ee = E*E0/k0
            IF((ee.GT.8.d0).and.(ee.LT.24.d0)) THEN
              ele =(ee - 8.d0)/16.d0
              red = (1.d0 - ele)*resultH(i,j) + ele*resultJ(i,j)
              resultL(i,j) = red
            ELSE IF(ee.GE.24.d0) THEN
              resultL(i,j) = resultJ(i,j)
            ELSE IF((eL.GT.0.1d0).and.(eL.LT.0.6d0)) THEN
              ele =( eL - 0.1d0)/0.5d0
              red = (1.d0 - ele)*resultK(i,j) + ele*resultH(i,j)
              resultL(i,j) = red
            ELSE IF(eL.LE.0.1d0) THEN
              resultL(i,j) = resultK(i,j)
            ELSE
              resultL(i,j) = resultH(i,j)
            END IF
         END DO
      END DO
c
      DO j=1,LPx
         kx = dexp(dp1*DBLE(j - 1))*2.d0  ! initial photonn energy in x
         k0 = dsqrt(kx**2)
         dk = (k0 - 2.d0)/DBLE(LE0)       ! MINus 2*mc^2
         totL = 0.d0
         DO i = 1,LE0
            IF(k0.GT.2.d0) THEN
              k = dk/(k0 - 2.d0)
              totL = totL + resultL(i,j)*k
            END IF
         END DO
         totalH(j)=totL
      END DO
c
      DO j=1,LPx
         totL = 0.d0
         DO i = 1,LE0
            totL = totL + resultL(i,j)/DBLE(LE0)
         END DO
         totalH2(j)=totL
      END DO

      IF(myrank.EQ.0) THEN
        WRITE(fo_name2,444) TRIM(data_file)//'Bethe-Heitler', jobno
        OPEN (33,file = fo_name2,form='formatted',status='unknown')
c
        DO i=1,LPx
           WRITE(33,*) i, totalH(i), totalH2(i)
        END DO
        CLOSE(33)
      END IF
444   FORMAT(A,I4.4,'.dat')
      RETURN
      END
C-----------------------------------------------------------------
      SUBROUTINE welcome
C-----------------------------------------------------------------
      USE mpi_common
      USE sim_common
      IF(myrank.NE.0) RETURN
      WRITE(*,*) ' '//achar(27)//'[34m'
      WRITE(*,*) ' '//achar(27)//'[1m'
      WRITE(*,*) '	###L      ########E    #######C      '
      WRITE(*,*) '	###L      ########E   ##########C  	 '
      WRITE(*,*) '	###L      ###E       ###C    ###C    '
      WRITE(*,*) '	###L      ########E  ###C            '
      WRITE(*,*) '	###L      ########E  ###C            '
      WRITE(*,*) '	###L      ###E       ###C    ###C    '
      WRITE(*,*) '	#######L  ########E   ##########C    '
      WRITE(*,*) '	#LASER##  ELECTRON#    COLLISION     '
      WRITE(*,*) ' '//achar(27)//'[32m'
      WRITE(*,*) 'Welcome to Laser Electron Collision code (v - 1.3.0)'
      WRITE(*,*) ' '//achar(27)//'[0m'
      WRITE(*,*) '*****************************************************'
      WRITE(*,*) "The code is running on", NPROCS, " processors"
      WRITE(*,*) "                      ", threads, " threads"
      WRITE(*,*) '*****************************************************'
      IF(load.EQ.1)
     &	WRITE(*,*) "Load particles from: f_E_smoothed_new_final.txt"
      IF(OutRad.EQ.1 ) WRITE(*,*)  'Produce radiation'
      IF(OutPairs.EQ.1) WRITE(*,*) 'Produce pairs'
      WRITE(*,*)
      RETURN
      END
