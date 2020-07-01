      module random_common
      integer, parameter :: icpu=4
      integer,parameter      :: LE0 =1000   ,LE1=1000
      real(kind=8),parameter :: Zcom=   79.d0   ! Z component for nucl
      real(kind=8),parameter :: Zcm3=   4.3d0   ! zcm3 = zcom**(1/3)
c     real(kind=8),parameter :: Zcom=   13.d0  ! Z component for nucl
c     real(kind=8),parameter :: Zcm3=   2.35d0 ! zcm3 = zcom**(1/3)
      integer,parameter      :: LPx =200   
      real(kind=8),parameter :: BPx =200.d0
      real(kind=8),parameter :: Emax=10000.d0
c
      real(kind=8) :: dp1
      real(kind=8),save,dimension(LE0,LPx) :: resultL
      real(kind=8),save,dimension(    LPx) :: totalH,totalH2
      end module random_common
c
      module sim_common
      use random_common
      integer      :: i,j,k,ksmax,ksout,kk,kstep,ii,ipl
     &               ,itotal,ksoutP,iconR,SKL,LL,polar,shape
     &               ,photon,species,shot,itotal0,L,OutRad,OutPairs
     &		     ,QED,sampled,sampled2,sampled3,sampled4,load		   
      logical      :: emmits,exists,load_particle
      real(kind=8) :: Xe,Ye,Ze,T,VX,VY,VZ
      real(kind=8) :: AVEX,AVEY,AVEZ,AVBX,AVBY,AVBZ
      real(kind=8) :: PXS,PYS,PZS,BXT,BYT,BZT
      real(kind=8) :: Rx,Rm,Rt,Rd1,Rd2,Rt2
      real(kind=8) :: pin,rin,div,div2,xinit,bin,enum,rmass
      real(kind=8) :: Wx0,Wy0,Wz0,Vx0,Vy0,Vz0,Wx,Wy,Wz
      real(kind=8) :: gg,ff,dF0,dVx,dVy,dVz
      real(kind=8) :: fvx,fvy,fvz
      real(kind=8) :: gammi,TTT,ENE00,RRR,RAD,TT1,GAMM
      real(kind=8) :: alf,palf,alf2,alpha,comp,sigmax,sigmay,sigmaz
      real(kind=8) :: SL,EL,BL,Rd,VL,Em,EV,inc_ang
      real(kind=8) :: DX,DT,DTI,DT1,phaseX,phaseY,phaseZ,TT
      real(kind=8) :: ww,wp,wk,pp,pw,pi,pi2,ppt,sp,ws,lb,wb
      real(kind=8) :: we0i,wpi,wsi,Pksout,wmin,www,Tm,Um,FF1,ENEd
      real(kind=8) :: ENEh,ENEv,ENE0,ENE,EmaxV,we0,XI,ENN,p_x,ENEA
      real(kind=8) :: x0,y0,z0,Tx,Ty,Tz,percentage
      real(kind=8) :: PQM       ! charge/mass ; rmass=-1.d0 for electron
      real(kind=8),dimension(:,:), allocatable   :: RE,RH
      real(kind=8),dimension(:,:), allocatable   :: phn,cphn
      real(kind=8) :: ERD,EKE,EKK,ERK
      real(kind=8),dimension(0:3000+1,200):: diffC,diffQ,diffD,diffR
      real(kind=8),dimension(200) :: totalR,totalC,totalP,totalRC
     &                              ,totalS,totalT
      real(kind=8) :: rand,rand1
      integer,save :: tmove,tmove0,tmove1
     &		   ,trdct,trdct0,trdct1
     &		   ,tmall,tmall0, t_max
     &		   ,tinit,tinit0,t_rate
     &		   ,tcurr,tcurr0,tcurr1
      end module sim_common
c
      module mpi_common
      integer :: myrank,ierr,nprocs,threads
      end module mpi_common
c
      module out_common
      integer :: jobno
      CHARACTER(LEN=100):: fo_name2,cwd
      CHARACTER(LEN=10) :: filename = 'input.dat'
      CHARACTER(LEN=100):: data_dir,data_file,input_file
      CHARACTER(LEN=*), PARAMETER :: data_dir_file = 
     &					'USE_DATA_DIRECTORY'
      end module out_common

      module R_common
      integer      :: Lall,LP,IPTSS,IPTFF,jj
      integer,dimension(7) :: Ne7
      real(kind=8) :: wmitT3,vmitT3,wmitTT,vmitTT,hh
      real(kind=8),dimension(:,:,:), allocatable:: phtn
      real(kind=8),dimension(:  ), allocatable:: wmit3,vmit3
      real(kind=8),dimension(:,:), allocatable:: emit3,fmit3,qmit3
      real(kind=8),dimension(:  ), allocatable:: emitT3,fmitT3,qmitT3
      real(kind=8),dimension(:  ), allocatable:: emitTT,fmitTT,qmitTT
      real(kind=8) :: wight00,const
      real(kind=8),dimension(:), allocatable :: wight0
      real(kind=8),dimension(:), allocatable :: wight  
      real(kind=8) :: TmY,TmZ
      end module R_common
c
      program main
      use random_common
      use sim_common
      use mpi_common
      use R_common
      use out_common
      use omp_lib
      implicit none
      include "mpif.h"
c
      NAMELIST /PARAM1/ jobno,ksmax,div,div2			       !time & grid parameters
      NAMELIST /PARAM2/ SL,Ev,pw,pp,sp				       !laser parameters
      NAMELIST /PARAM3/ alpha,enum,bin,shot,inc_ang		       !electron parameters
      NAMELIST /PARAM4/ xinit,rmass,sigmax,sigmay,sigmaz           !configurations
      NAMELIST /PARAM5/ iconR,QED,ipl,shape,load,OutRad,OutPairs   !physical processes
c
      call system_clock(tmall0)
      call MPI_INIT(IERR)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)	!initialize MPI
c
      tmove=0 ; tcurr=0 ; trdct=0
      tmall=0 ; tinit=0 
c
      PI  = 4.D0*DATAN( 1.D0 )
      PI2 = PI*2.D0
c
      call getcwd(cwd)
      IF (myrank == 0) THEN
         OPEN(unit=41, status='OLD', file=TRIM(data_dir_file)
     &	, iostat=ierr)
         IF (ierr == 0) THEN
      	    READ(41,'(A)') data_dir
            CLOSE(41)
      	    PRINT*, 'Using data directory "' 
     &	          // TRIM(data_dir) // '"'
         ELSE
     	    PRINT*, 'Specify output directory'
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
     &            // '" does not exist.'
        PRINT *, 'Create the file and rerun the code.'
      END IF
c
      if(myrank==0) then 
         WRITE(*,*) 'Output directory:',data_file
      end if

      open ( 8,status="old",file=TRIM(input_file),form='formatted')	!read input data file
      read ( 8,PARAM1) 
c
      jobno = jobno + myrank
      write(fo_name2,444) TRIM(ADJUSTL(data_file))//'output',jobno
      open ( 9,file=fo_name2,form='formatted')
c
      write( 9,PARAM1)
      read ( 8,PARAM2) ; write(9,PARAM2) 
      read ( 8,PARAM3) ; write(9,PARAM3)
      read ( 8,PARAM4) ; write(9,PARAM4) 
      read ( 8,PARAM5) ; write(9,PARAM5)
      close( 8)
c
      write(9,*) "Check openMP"
!$omp parallel
      threads = omp_get_num_threads()
!$omp end parallel
      write(9,*) "total number of threads=",threads
c 
      call setprm  						!parameter setup
      call welcome
c
      call system_clock(tinit0)
      call setbeam
c
      write(fo_name2,444) TRIM(data_file)//'orbt1q',jobno
      open (10,file=fo_name2,form='formatted',status='REPLACE')
      write(fo_name2,444) TRIM(data_file)//'orbt2q',jobno
      open (11,file=fo_name2,form='formatted',status='REPLACE')
      write(fo_name2,444) TRIM(data_file)//'orbt3q',jobno 
      open (12,file=fo_name2,form='formatted',status='REPLACE')
      write(fo_name2,444) TRIM(data_file)//'orbt4q',jobno
      open (13,file=fo_name2,form='formatted',status='REPLACE')
      write(fo_name2,444) TRIM(data_file)//'orbt5q',jobno
      open (14,file=fo_name2,form='formatted',status='REPLACE')
      write(fo_name2,444) TRIM(data_file)//'orbt6q',jobno
      open (15,file=fo_name2,form='formatted',status='REPLACE')
      write(fo_name2,444) TRIM(data_file)//'orbt7q',jobno
      open (17,file=fo_name2,form='formatted',status='REPLACE')
      call system_clock(tinit )
      tinit = tinit-tinit0
c
c     orbit calculation 
c
      T      = -wp*xinit
      RAD    = 0.d0
      TT1    = 0.d0 
      photon = 0
      kstep  = 0 
c
c      call random_seed()
1000	CONTINUE              					! time loop
	   KSTEP = KSTEP + 1
           T=T+dt
	call system_clock(tmove0)        
         if(iconR.eq.0) call Lorentz          
         if(iconR.eq.1) call Sokolov            
         if(iconR.eq.2) call Landau_Lifshitz            
         if(iconR.eq.3) call qedemmit		      !particle loop
         call outorbit2 				      !output extracted electron orbit
 	call system_clock(tmove1)
        tmove = tmove + (tmove1-tmove0)
c
         if(OutRad.eq.1)             call outorbit1 	!store orbit data for emission calculation    
	   if((mod(kstep,ksout).eq.0)
     &     .and.(alpha.gt.0.d0)    ) call avgene	!calculate average energy
c
           if((mod(kstep,ksmax).eq.0)
     &     .and.(alpha.gt.0.d0)    ) call histogram
c
 	   if((mod(kstep,ksmax).eq.0)
     &     .and.(alpha.gt.0.d0)    ) call histogram2d
c
	   if((mod(kstep,1000).eq.0)
     &     .and.(myrank.eq.0))
     &	write(*,*) "Iteration =",kstep,"; Time step =",Rt*kstep
c---------------------------------------------------------
      IF(KSTEP.GE.KSMAX) GO TO 2500
      GO TO 1000						       ! end time loop
c
c     end  orbit calculation
c
      close(10) ; close(11) ; close(12)
      close(13) ; close(14) ; close(15)
      close(16) ;	close(20)
c
2500  CONTINUE
	deallocate(RE,RH)
c
      if(OutRad.eq.1) then
      	call angdis 
	call radiation
      end if
c      if(OutRad.eq.1) call photon_his
c
 333  format(A,I2.2,'_',I2.2,'.dat') 
 444  format(A,I3.3,'.dat')
 555  format(A,I2.2,'_',I2.2,'_',I2.2,'.dat')
 666  format(11(E16.6,1X))
c
      if(OutRad.eq.1) deallocate(phtn)
	deallocate(wight0,wight)
c
9999  continue
c
      call system_clock(tmall,t_rate,t_max)
      tmall = tmall-tmall0
      write(9,*) "all      :",dble(tmall)/dble(t_rate)
      write(9,*) "initc    :",dble(tinit)/dble(t_rate)
      write(9,*) "pmove    :",dble(tmove)/dble(t_rate)
      write(9,*) "rad      :",dble(tcurr)/dble(t_rate)
      write(9,*) "reduction(rad):",dble(trdct)/dble(t_rate)
      write(9,*) "others   :",dble(tmall-tinit-tmove-tcurr
     &				-trdct)/dble(t_rate)
      close(9)
      call MPI_FINALIZE(IERR)
      if(myrank.eq.0) 
     &   write(*,*) 
     &   "Final runtime =",tmall/t_rate,"seconds"
c
      STOP
      END
!--------------------------------------
	subroutine setprm
!--------------------------------------
      use random_common
      use sim_common
      use mpi_common
      use R_common
      use out_common
      use omp_lib
      implicit none
      CHARACTER(LEN=8)  :: qed_file = 'vac4.dat'
      CHARACTER(LEN=100):: table_location

      VL = 3.d0*1.0d8              ! light speed           [m/s]
      EL = 2.74d0*1.0d3 *dsqrt(SL) ! peak electric field   [V/m]
      BL = 9.28d0*1.0d-2*dsqrt(SL) ! peak magnetic field   [Gauss]
      Rm = 1.7d0*1.0d1/BL          ! typical Larmor radius [m]
      Rd1 = Rm/div                 ! grid size             [m]
      Rd2 = pw/div2
      Rd = min(Rd1,Rd2)
      Rx = Rm                      ! unit scale            [m]
      Rt = 0.5d0*Rd/VL             ! time step             [s]
      Rt2= Rt*div                  ! unit time             [s]
      dx = dble(1.d0/div)*Rd/Rd1   ! normalized grid size
      dt = 0.5d0*dx                ! normalized time step
      Em = Ev/(0.511d6*rmass)      ! nomalized beam energy 
      Vx = dsqrt((Em+1.d0)**2-1.d0)! normalized 4-velocity (+x)
      Wx = Vx/(Em+1.d0)            ! normalized 3-velocity (+x)
      Vy = 0.d0                    ! normalized 4-velocity (+y)
      Vz = 0.d0                    ! normalized 4-velocity (+z)
      ppt = pp/VL
      ww = pw/Rd*dx                ! normalized wavelength
      wk = pi2/ww                  ! normalized wave number
      wp = pp/Rd*dx                ! normalized pulse length
      ws = sp/Rd*dx                ! normalized waist radius
      wpi= 1.d0/wp
      wsi= 1.d0/ws
      pin= bin/(0.511d6*rmass)     ! normalized bin size
      rin= Em/6000.d0              ! reference bin size 
      wb = alpha/Rd*dx             ! normalized transverse beam size
      PQM =-1.d0/rmass
c
      write(9,*) " "
      write(9,*) "Parameters for pulse laser"
      write(9,*) " "
      if(ipl.eq.0) then
         write(9,*) "Laser polarization: linear"
         polar = 0.d0
      else
         write(9,*) "Laser polarization: circular"
         polar = 1.d0
      end if
      write(9,*) " "
	if(shape.eq.0) write(9,*) "0th order Gaussian beam"
	if(shape.eq.1) write(9,*) "5th order Paraxial approximation"
	write(9,*) " "
      write(9,*) "Laser Intensity               [W/cm2]",SL
      write(9,*) "Peak electric field             [V/m]",EL
      write(9,*) "Peak magnetic field           [Gauss]",BL
      write(9,*) "Larmor radius for light speed     [m]",Rm
      write(9,*) "laser wavelength                  [m]",pw
      write(9,*) "pulse length                      [m]",pp
      write(9,*) "pulse duration                    [s]",ppt
      write(9,*) "pulse duration (FWHM)             [s]",ppt*1.1774d0
      write(9,*) "waist radius (1/e2)               [m]",sp
c
      write(9,*) " "
      write(9,*) "parameters for electron beam"
      write(9,*) " "
      if(load.eq.1) then
		load_particle=.TRUE.
		write(9,*) "Load particles from: f_E_smoothed_new_final.txt"
      else
		load_particle=.FALSE.
      end if
      write(9,*) " "
      write(9,*) "parameters for electron beam"
      write(9,*) "initial kinetic energy           [eV]",Ev
      write(9,*) "initital momentum              [p/mc]",Vx
      write(9,*) "initital 3-velocity             [v/c]",Wx
      if(alpha.gt.1) then
	     write(9,*) "momentum spread x                 [%]",sigmax*100
	     write(9,*) "momentum spread y                 [%]",sigmay*100
	     write(9,*) "momentum spread z                 [%]",sigmaz*100
      end if
      write(9,*) "number of shot                       ",shot
c
      write(9,*) " "
      write(9,*) "parameters for computation"
      write(9,*) " "
      write(9,*) "Larmor radius / grid separation      ",div
      write(9,*) "Laser wavelength / grid separation   ",div2
      write(9,*) "grid separation,Rd1,Rd2           [m]",Rd1,Rd2
      write(9,*) "grid separation,Rd=min(Rd1,Rd2)   [m]",Rd
      write(9,*) "physical time step interval       [s]",Rt
      write(9,*) "physical unit time                [s]",Rt2
      write(9,*) "energy bin size of spectrum      [eV]",bin
c
      palf = 6.2d-24
      alf  = palf*dt/Rt
      alf2 = alf*137.d0*9.d0/4.d0
      comp  = palf*137.d0*VL*3.d0*PI 
      write(9,*) "coefficient for radiation         [s]",palf 
      write(9,*) "Compton wavelength                [m]",comp
      if(iconR.eq.0) then
         write(9,*) " "
         write(9,*) "Particle pusher: Lorentz"
         write(9,*) " "
      else if(iconR.eq.1) then
         write(9,*) " "
         write(9,*) "Particle pusher: Sokolov"                                                          
         write(9,*) " "
      else if(iconR.eq.2) then
         write(9,*) " "
         write(9,*) "Particle pusher: reduced Landau Lifshitz"
         write(9,*) " "
	else if(iconR.eq.3) then
	 write(9,*) " "
         write(9,*) "Particle pusher: Lorentz"
         write(9,*) "        Use QED: Stochastic"
	 write(9,*) " "
      end if
c	
      if(QED.eq.0) then
         write(9,*) "Use QED: No, Classical"
         write(9,*) " "
      else
         write(9,*) "Use QED: QED_assisted"
         write(9,*) " "
      end if
c
      write(9,*) "max time step",ksmax
      ksout = 1
      ksoutP= 1
      if(ksmax.ge.  10000) ksout = ksmax/10000
      if(ksmax.ge.1000000) ksoutP= ksmax/1000000
      Pksout= dble(ksoutP)
c setup for theoretical cross sections
c
      EmaxV= 1000.d0 
      we0 =(log(EmaxV)-log(0.001d0))/199.0 !! 
      we0i= 1.d0/we0
c
c radiation & pair production in strong-field
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
      call getcwd(cwd)
      table_location = TRIM(cwd)//'/TABLE/'//TRIM(ADJUSTL(qed_file))
      IF (myrank == 0 ) THEN
      INQUIRE(file=TRIM(table_location), exist=exists)
      IF (.NOT.exists) THEN
          WRITE(*,*) '*** ERROR ***'
          WRITE(*,*) 'Unable to find QED tables in the ', 
     &    'directory "' // TRIM(table_location) // '"'
      END IF
      END IF
      
      open(96,file=table_location,form='unformatted')
      read(96) totalR,totalC,totalP,totalS,totalT
      read(96) diffC,diffQ,diffD,diffR
      close(96)

      if(QED.eq.1) then 
         totalRC=totalR/totalC
      else
         totalRC=1.d0
      end if

      return
      end
!--------------------------------------
	subroutine setbeam
!--------------------------------------
      use random_common
      use sim_common
      use mpi_common
      use R_common
      use omp_lib
      use out_common
      implicit none
      real(kind=8) :: aa
c generate initial electron conditions for incident beam  
c---------------------------
      if(alpha.ne.0.d0) then
c---------------------------
        sampled = 12 ; sampled2 = sampled/2 
	sampled3 = (sampled-1)**3

	allocate(Re(11,sampled3))
        allocate(wight0(sampled3))

	write(9,*) "transverse beam size   [m]",alpha
 	write(9,*) "electron number per shot  ",enum
	write(9,*) "incident angle    [degree]",inc_ang
        inc_ang = inc_ang/180.d0*pi

	write(fo_name2,444) TRIM(data_file)//'dist_sp',jobno
        open (10,file=fo_name2,form='formatted',status='unknown')

	do k=1,sampled-1
        do j=1,sampled-1
        do i=1,sampled-1
         kk = (k-1)*(sampled-1)*(sampled-1)+(j-1)*(sampled-1)+i
	   call random_number(rand)
	   call random_number(rand1)
	   phaseX  = (dble(k-sampled2))/(sampled2-1)*0.707d0
         phaseY  = (dble(i-sampled2))/(sampled2-1)*0.707d0
         phaseZ  = (dble(j-sampled2))/(sampled2-1)*0.707d0
         Re(1,kk)= wp*xinit + phaseX*wb                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
         Re(2,kk)= phaseY*wb
         Re(3,kk)= phaseZ*wb
	   write(10,666) Re(1,kk),Re(2,kk),Re(3,kk)
         wight0(kk)= dexp(-phaseX**2-phaseY**2-phaseZ**2)
c
	   if(load_particle) call manual_load
c
	   Re(4,kk) = Vx*(-1.d0) + sigmax*dcos(2.d0*pi*rand1)
     &		  *Vx*dsqrt(-2.d0*log(rand))
	   Re(5,kk) = sigmay*Vx*dcos(2.d0*pi*rand1)
     &		  *dsqrt(-2.d0*log(rand))
	   Re(6,kk) = sigmaz*Vx*dsin(2.d0*pi*rand1)
     &		  *dsqrt(-2.d0*log(rand))
        end do
	end do
	end do
        wight00=0.d0
        do i=1,sampled3
          wight00 = wight00 + wight0(i)
        end do
        wight0=wight0/wight00
c
        do i=1,sampled3
         Xe = Re(1,i)
         Ye = Re(2,i)
         Vx0 = Xe*dcos(inc_ang) - Ye*dsin(inc_ang)
         Vy0 = Xe*dsin(inc_ang) + Ye*dcos(inc_ang)
	   Re(1,i) = Vx0 
	   Re(2,i) = Vy0
        end do
c
        do i=1,sampled3
         Vx0 = Re(4,i)
         Vy0 = Re(5,i)
         Re(4,i) = Vx0*dcos(inc_ang)-Vy0*dsin(inc_ang)
         Re(5,i) = Vx0*dsin(inc_ang)+Vy0*dcos(inc_ang)
        end do

        itotal = int(sampled3/dble(nprocs))
	allocate(Rh(11,itotal))
        allocate(wight(itotal))
        jj = myrank*itotal
        write(9,*) "myrank                            ",myrank
        write(9,*) "total process number, nprocs        ",nprocs
        write(9,*) "calculated particle number, itotal",itotal
        write(9,*) "calculated particle number from",jj+1,"to",jj+itotal
        do i=1,itotal
         Rh(1,i) = Re(1,i+jj)
         Rh(2,i) = Re(2,i+jj)
         Rh(3,i) = Re(3,i+jj)
         Rh(4,i) = Re(4,i+jj)
         Rh(5,i) = Re(5,i+jj)
         Rh(6,i) = Re(6,i+jj)
         wight(i)= wight0(i+jj)
        end do
	Re = Rh
c
      	if(OutRad.eq.1) then
		   allocate(phtn(6,ksmax,itotal))
		   phtn = 0.d0
		end if
c	
	call histogram
 	call histogram2d
c----------
      else
c----------
	   write(9,*) "single electron"
	   write(9,*) "electron number per shot = 1  "
	   allocate(Re(11,1),Rh(11,1))
	   allocate(wight0(1),wight(1))
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
		if(OutRad.eq.1) then
			allocate(phtn(6,ksmax,1)) 
			phtn = 0.d0
		end if
      end if
c
      do i=1,7
         Ne7(i)=(itotal/8)*i+1
         write(9,*) "sampling electron number",i,Ne7(i)
      end do
c     end incident beam generation 
c
444   format(A,I4.4,'.dat')
666   format(3(E14.4,1X))
      return
      end
!--------------------------------------
      subroutine outorbit1
!-------------------------------------
      use random_common
      use sim_common
      use mpi_common
      use R_common
      use omp_lib
      implicit none
      real(kind=8)    :: TTY,TTZ
	if(mod(kstep,ksoutP).eq.0) then 
        j=kstep/ksoutP
        if(j.le.1000000) then 
          do i=1,itotal
             Vx =Re( 4,i)
             Vy =Re( 5,i)
             Vz =Re( 6,i)
	       ff =Re( 7,i)
             TTT=Re( 9,i)
             Xi =Re(10,i)
             ENE=Re(11,i)
             TTY = dsqrt(Vx**2+Vy**2)
             TmY = acos(Vx*(-1.d0)/TTY)
             if(Vy.lt.0.d0) TmY=(-1.d0)*TmY
 		 TTZ = dsqrt(Vx**2+Vz**2)
             TmZ = acos(Vx*(-1.d0)/TTZ)
             if(Vz.lt.0.d0) TmZ=(-1.d0)*TmZ   ! direction angleZ
             phtn(1,j,i)=Xi 
             phtn(2,j,i)=ENE
             phtn(3,j,i)=TmY 
             phtn(4,j,i)=TTT
             phtn(5,j,i)=ff
    		 phtn(6,j,i)=TmZ
          end do
        end if
	end if
      return
      end
!--------------------------------------
	subroutine outorbit2
!--------------------------------------
      use random_common
      use sim_common
      use mpi_common
	use R_common
      use omp_lib
      implicit none

	if(mod(kstep,ksout).eq.0 ) then
      do i=1,7
         k =Ne7(i)
         Vx=Re(4,k)
         Vy=Re(5,k)
         Vz=Re(6,k)
         ENE=dsqrt(1.d0+Vx**2+Vy**2+Vz**2)
         write(9+i,666) t*Rt2,Re(1,k)*Rx,Re(2,k)*Rx	!t, x, y 
     &             ,Re(4,k),Re(5,k),(ENE)*0.511d6	!px, py, K.E
     &	       ,abs(TTT)*0.511d6,RAD*0.511d6,Xi	!Wem,Wrad,chi_e
      end do
	end if
 666  format(11(E16.6,1X))
      return
      end
!--------------------------------------
	subroutine avgene
!--------------------------------------
      use random_common
      use sim_common
      use mpi_common
	use R_common
	use out_common
      use omp_lib
      implicit none
	real(kind=8) :: SIG_neg,SIG_pos
     &		   ,SIGE_pos,SIGE_neg,sss
      integer      :: num_p,num_n,num_pp,num_nn
      integer,dimension(icpu)      :: num_pos, num_neg
	real(kind=8),dimension(icpu) :: SE,SE_pos,SE_neg,SR
	include "mpif.h"

      SE = 0.d0 ; SR = 0.d0
	itotal0 = int(itotal/icpu+1)
      do LP = 1,icpu
         IPTSS  = itotal0*(LP-1)+1
         IPTFF  = min(itotal0* LP,itotal)
         do L  = IPTSS,IPTFF
            Vx  = Re(4,L)
            Vy  = Re(5,L)
            Vz  = Re(6,L)
 		RRR = Re(9,L)
            GAMM = sqrt(1.d0+Vx*Vx+Vy*Vy+Vz*Vz)
            SE(LP)=SE(LP)+GAMM
		SR(LP)=SR(LP)+RRR
         end do
      end do
c
      GAMM=0.d0
      RRR =0.d0
      do LP=1,icpu,4
         GAMM = GAMM + SE(LP)  + SE(LP+1)
     &		   + SE(LP+2)+ SE(LP+3)
	   RRR  = RRR  + SR(LP)  + SR(LP+1)
     &		   + SR(LP+2)+ SR(LP+3)
      end do

      EKK=0.d0
      ERK=0.d0
c---------------------------------------------------------
	  call mpi_allreduce(GAMM,EKK,1,MPI_REAL8,MPI_SUM
     &                                ,mpi_comm_world ,ierr)
	  call mpi_allreduce(RRR,ERK,1,MPI_REAL8,MPI_SUM
     &                                ,mpi_comm_world ,ierr)
c
	EKK = EKK /(itotal*nprocs) ; ERK = ERK /(itotal*nprocs)
	SE_neg = 0.d0 ; SE_pos = 0.d0 
	num_pos=0 ; num_neg=0 
      do LP = 1,icpu
         IPTSS  = itotal0*(LP-1)+1
         IPTFF  = min(itotal0* LP,itotal)
         do L  = IPTSS,IPTFF
            Vx  = Re(4,L)
            Vy  = Re(5,L)
            Vz  = Re(6,L)
            GAMM = sqrt(1.d0+Vx*Vx+Vy*Vy+Vz*Vz)
	      sss = GAMM-EKK
            if(sss.lt.0.d0) then
		  SE_neg(LP) = SE_neg(LP) + abs(sss)
              num_neg(LP) = num_neg(LP) + 1
            else 
      	  SE_pos(LP) = SE_pos(LP) + abs(sss)
              num_pos(LP) = num_pos(LP) + 1 
		end if
         end do
      end do
c
	SIG_pos=0.d0 ; SIG_neg =0.d0
      num_p = 0 ; num_n = 0
      do LP=1,icpu,4
         SIG_pos= SIG_pos+ SE_pos(LP)  + SE_pos( LP+1)
     &		       + SE_pos(LP+2)+ SE_pos( LP+3)
         SIG_neg= SIG_neg+ SE_neg(LP)  + SE_neg( LP+1)
     &		       + SE_neg(LP+2)+ SE_neg( LP+3)
          num_p = num_p  + num_pos(LP) + num_pos(LP+1)
     &                   + num_pos(LP+2)+num_pos(LP+3)
          num_n = num_n  + num_neg(LP) + num_neg(LP+1)
     &                   + num_neg(LP+2)+num_neg(LP+3)
      end do

c     
	SIGE_pos=0.d0 ; SIGE_neg=0.d0 
	num_pp=0 ; num_nn=0
 	call mpi_allreduce(SIG_pos,SIGE_pos,1,MPI_REAL8,MPI_SUM
     &                                ,mpi_comm_world ,ierr)
 	call mpi_allreduce(SIG_neg,SIGE_neg,1,MPI_REAL8,MPI_SUM
     &                                ,mpi_comm_world ,ierr)
 	call mpi_allreduce(num_p,num_pp,1,MPI_INTEGER,MPI_SUM
     &                                ,mpi_comm_world ,ierr)
      call mpi_allreduce(num_n,num_nn,1,MPI_INTEGER,MPI_SUM
     &                                ,mpi_comm_world ,ierr)
c	
	SIGE_pos=SIGE_pos/num_pp
	SIGE_neg=SIGE_neg/num_nn
c
	if(myrank==0) then
        write(fo_name2,444) TRIM(data_file)//'AveEne',jobno
        open (20,file=fo_name2,form='formatted',status='unknown')
        write(20,666) t*Rt2,EKK*0.511d6,ERK*0.511d6
     &		 ,(EKK+SIGE_pos)*0.511d6
     &		 ,(EKK-SIGE_neg)*0.511d6
	end if
 444  format(A,I3.3,'.dat')
 666  format(11(E16.6,1X))
      return
      end
!--------------------------------------
	subroutine radiation
!--------------------------------------
      use random_common
      use sim_common
      use mpi_common
	use R_common
	use out_common
      use omp_lib
      implicit none
	include "mpif.h"
	real(kind=8),dimension(:), allocatable :: total
	real(kind=8),dimension(:,:), allocatable :: diff1,diff2
      if(myrank.eq.0) write(*,*) "calculating radiation"
c
c distribute electron energy reduction onto emission energy spectum

      allocate(wmit3(icpu),vmit3(icpu))
	allocate(emit3(6000,icpu)) ; allocate(fmit3(6000,icpu))
	allocate(emitT3(6000)  ) ; allocate(fmitT3(6000)  ) 
      allocate(emitTT(6000)  ) ; allocate(fmitTT(6000)  )
	allocate(total(6000)) 
      allocate(diff1(0:3000+1,200))
      allocate(diff2(0:3000+1,200))
c
      call system_clock(tcurr0)
      if(QED.eq.1) then
	  total=totalR
        diff1=diffR
	  diff2=diffQ
	else if(QED.eq.0) then
	  total=totalC
        diff1=diffD
	  diff2=diffC
      end if
!$omp workshare
      ENEh = 1.0d0*Em
      ENEd = ENEh/6000.d0
      Lall = int(ksmax/icpu) 
      emit3= 0.d0
      fmit3= 0.d0
      wmit3= 0.d0
      vmit3= 0.d0
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
      do LP=1,icpu                 ! parallelization loop
         IPTSS =    (LP-1)*Lall+1
         IPTFF =     LP   *Lall
      do j=IPTSS,IPTFF           ! timestep loop			
      do L=1,itotal              ! particle loop
         Xi = phtn(1,j,L)        ! quantum parameter
         ENE= phtn(2,j,L)        ! energy
         Um = phtn(4,j,L)*Pksout ! energy defference
         hh = phtn(5,j,L)*Pksout ! energy defference
c
         if(XI.gt.0.001) then
           kk = min(idnint(log(XI*1000.d0)*we0i+1.5d0),200)
c integration of total cross section
           ff = total(kk)*3000.d0
c coefficient
           FF1= ENE/3000.d0
           gg = Um/ff*ENEd/FF1*wight(L)
c
           do k=1,6000
              ENE0 = (dble(k)-0.5d0)*ENEd
              ii   = idnint(ENE0/ENE*3000.d0+0.5d0)
              Xe   = ENE0/ENE*3000.d0+0.5d0 - dble(ii)
              if(ii.gt.2999) exit
              TTT  = diff1(ii,kk)+Xe*(diff1(ii+1,kk)-diff1(ii,kk))
              PXS  = diff2(ii,kk)+Xe*(diff2(ii+1,kk)-diff2(ii,kk))
              emit3(k  ,LP)=emit3(k  ,LP)+TTT*gg/ENE
              fmit3(k  ,LP)=fmit3(k  ,LP)+PXS*gg
           end do
        end if
        wmit3(LP)=wmit3(LP)+Um*wight(L)
        vmit3(LP)=vmit3(LP)+hh*wight(L)
      end do
      end do
      end do
!$omp end parallel do 

c
!$omp parallel do private(i,k,Tm,ff) 
!$omp&      shared(emit3,fmit3,emitT3,fmitT3) 
      do i=1,6000
         Tm = 0.d0
         ff = 0.d0
         do k=1,icpu,4
            Tm = Tm + emit3(i,k)  + emit3(i,k+1)
     &		  + emit3(i,k+2)+ emit3(i,k+3)
            ff = ff + fmit3(i,k)  + fmit3(i,k+1)
     &		  + fmit3(i,k+2)+ fmit3(i,k+3)
         end do
         emitT3(i)=Tm
         fmitT3(i)=ff
      end do
!$omp end parallel do
c
      wmitT3=0.d0
      vmitT3=0.d0
!$omp parallel do private(k) shared(wmit3,vmit3) 
!$omp&         reduction(+:wmitT3,vmitT3)
      do k=1,icpu,4
         wmitT3=wmitT3+wmit3(k)  +wmit3(k+1)
     &		    +wmit3(k+2)+wmit3(k+3)
         vmitT3=vmitT3+vmit3(k)  +vmit3(k+1)
     &		    +vmit3(k+2)+vmit3(k+3)
      end do
!$omp end parallel do
      call system_clock(tcurr1)
      tcurr=tcurr+tcurr1-tcurr0
c
c summation in MPI processes
c
      if(kstep.ne.ksmax) return  
c
c     goto 9999
c
	call system_clock(trdct0)
      vmitT3  = abs(vmitT3)
      call mpi_allreduce(emitT3,emitTT ,6000,mpi_real8
     &                                 ,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(fmitT3,fmitTT ,6000,mpi_real8
     &                                 ,mpi_sum,mpi_comm_world,ierr)
c
      call mpi_allreduce(wmitT3,wmitTT ,   1,mpi_real8
     &                                 ,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(vmitT3,vmitTT ,   1,mpi_real8
     &                                 ,mpi_sum,mpi_comm_world,ierr)
      call system_clock(trdct1)
      trdct=trdct+trdct1-trdct0
c
c end distribute onto emission energy spectum
c
      write(9,*) "energy reduction (  EM work) [J]"
     &						,vmitTT*enum*0.511d6*1.6d-19
      write(9,*) "energy reduction (radiation) [J]"
     &						,wmitTT*enum*0.511d6*1.6d-19
c
c output photon spectrum data
c
	const = pin/rin*shot*enum
	if(myrank==0) then
        write(fo_name2,444) TRIM(data_file)//'phtne',jobno
        open (19,file=fo_name2,form='formatted',status='unknown')
c      
        do i=1,6000
          ENE0 = (dble(i)-0.5d0)*ENEd
          write(19,*) ENE0*0.511d6,emitTT(i)*const,fmitTT(i)*const
        end do 
        close(19)
      end if
c
      ff=0.d0
!$omp parallel do private(k) shared(fmit3) reduction(+:ff)
      do i=1,6000
         ff=ff+fmitTT(i)*enum
      end do
!$omp end parallel do
c
      write(9,*) "total reduction energy [J]",ff*0.511d6*1.6d-19

	if(OutPairs.eq.0) return
c     setup for theoretical cross sections
	call pair_init
c
	if(myrank.eq.0) write(*,*) "Convert Pairs"
c
	allocate(qmit3(6000,30)) 
	allocate(qmitT3(6000)  )
      allocate(qmitTT(6000)  )
c
      ENEh = 1.d0*Em
      ENEd = ENEh/6000.d0
      Lall = 6000/30
      qmit3=0.d0
      do LP=1,30
         IPTSS = (LP-1)*Lall+1
         IPTFF =  LP   *Lall
      do i=IPTSS,IPTFF
         ENE0 = (dble(i)-0.5d0)*ENEd
	   if(ENE0.ge.2.d0) then 
         kk   = min(idnint(log(ENE0*0.5d0)/dp1+1.d0),LPx)
         ff   = totalH(kk)
         ENE  = (1.d0-dexp(-3.d-5*ff))*emitT3(i) 
c integration of total cross section 
         ff = totalH2(kk)*1000.d0
c coefficient 
         FF1= ENE0/1000.d0
         gg = ENE/ff*ENEd/FF1
c
         do k=1,6000
            ENEA = (dble(k)-0.5d0)*ENEd
            ii   = idnint(ENEA/ENE0*1000.d0+0.5d0)
            Xe   = ENEA/ENE0*1000.d0+0.5d0 - dble(ii)
            if(ii.ge.1000) exit
            if(ii.ge.1) then 
              TTT  = resultL(ii,kk)+Xe*(resultL(ii+1,kk)-resultL(ii,kk))
              qmit3(k   ,LP)=qmit3(k   ,LP)+TTT*gg
            end if
         end do
         end if
      end do
      end do

      do i=1,6000
         ff = 0.d0
         do k=1,30
            ff = ff + qmit3(i,k)
         end do
         qmitT3(i)=ff
      end do
c
c end convert to the positron track data using materical filter
c
c
c summation in MPI processes
c
      call mpi_allreduce(qmitT3,qmitTT,6000,mpi_real8
     &                                 ,mpi_sum,mpi_comm_world,ierr)
c
c output positron spectrum data
c
	if(myrank.eq.0) then
         write(fo_name2,444) TRIM(data_file)//'pairTT',jobno
         open (34,file=fo_name2,form='formatted',status='unknown')
c
         do i=1,6000
            ENE0 = (dble(i)-0.5d0)*ENEd
            write(34,*) ENE0*0.511d6,qmitTT(i)*const
         end do
	end if
c
      close(34)

      deallocate(wmit3 , vmit3)
	deallocate(emit3 , fmit3) 
	deallocate(emitT3,fmitT3) 
      deallocate(emitTT,fmitTT) 
      deallocate(total,diff1,diff2)
      deallocate(qmit3,qmitT3,qmitTT)
 444  format(A,I3.3,'.dat')
      return
      end

C=============================================================
      subroutine angdis
C=============================================================
      use random_common
      use sim_common
      use mpi_common
	use R_common
	use out_common
      use omp_lib
      implicit none
      include "mpif.h"
      integer,PARAMETER :: ang = 1024
	real(kind=8) :: thetay,thetaz,ang2
      real(kind=8),dimension(:,:,:), allocatable:: emit2
      real(kind=8),dimension(:,:  ), allocatable:: emitT2,emitTT2
c
	ang2=dble(ang/2)
      allocate(  emit2(ang,ang,icpu))
      allocate( emitT2(ang,ang     ))
	allocate(emitTT2(ang,ang     ))
c
	if(myrank.eq.0) write(*,*) "Calculating angular distribution..."
	
!$omp workshare
	Lall = 1000000/icpu + 1
      emit2= 0.d0
!$omp end workshare
c
!$omp parallel do private(LP,IPTSS,IPTFF,j,L,Xi
!$omp&  ,Um,TmY,TmZ,ii,jj)
!$omp&shared(phtn,PKsout,wight,emit2)
      do LP=1,icpu
         IPTSS = (LP-1)*Lall+1
         IPTFF = min(LP*Lall,1000000)
      do j=IPTSS,IPTFF
      do L=1,itotal
         Xi  = phtn(1,j,L)        ! quantum parameter
         TmY = phtn(3,j,L)
         TmZ = phtn(6,j,L)
         Um  = phtn(4,j,L)*Pksout ! energy defference
c
         if(XI.gt.0.001) then
           ii = min(idnint(TmZ/(pi)*ang2+0.5d0+ang2),ang)
           jj = min(idnint(TmY/(pi)*ang2+0.5d0+ang2),ang)  ! angle in radian
           emit2(ii,jj,LP)=emit2(ii,jj,LP) + Um*wight(L)
         end if
      end do
      end do
      end do
!$omp end parallel do

!$omp parallel do private(j,i,Tm,k) shared(emit2,emitT2) 
	do j=1,ang
	do i=1,ang
         Tm = 0.d0
         do k=1,icpu,4
            Tm = Tm + emit2(i,j,k  ) + emit2(i,j,k+1)
     &		  + emit2(i,j,k+2) + emit2(i,j,k+3)
         end do
         emitT2(i,j)=Tm
      end do
	end do
!$omp end parallel do
c
c summation in MPI processes
c
      if(kstep.ne.ksmax) return
c
      call mpi_allreduce(emitT2,emitTT2,ang*ang,mpi_real8
     &                                 ,mpi_sum,mpi_comm_world,ierr)
c end distribute onto emission energy spectum

c output photon spectrum data
      const = pin/rin*shot*enum

	if(myrank==0) then
        write(fo_name2,444) TRIM(data_file)//'phtnTe',jobno
        open (21,file=fo_name2,form='formatted',status='unknown')
c
        do j=1,ang
	  do i=1,ang
           thetaz=(i-0.5-int(ang2))/int(ang2)*(pi)
	     thetay=(j-0.5-int(ang2))/int(ang2)*(pi)
           write(21,666) thetaz,thetay,emitTT2(i,j)
        end do
        end do
        close(21)
	end if
c
      deallocate(emit2,emitT2,emitTT2)

444   format(A,I3.3,'.dat')
666   format(3(E14.4,1X))
      RETURN
      END
!--------------------------------------
      subroutine Lorentz
!--------------------------------------
      use random_common
      use sim_common
      use mpi_common
	use R_common
      use omp_lib
      implicit none
c
	itotal0 = int(itotal/icpu+1) 
!$omp do            
      do LP = 1,icpu
         IPTSS  = itotal0*(LP-1)+1
         IPTFF  = min(itotal0* LP,itotal)
         do i  = IPTSS,IPTFF      ! particle loop
	   DT1    = PQM*DT*0.5d0
         Xe =Re( 1,i)
         Ye =Re( 2,i)
         Ze =Re( 3,i)
         Vx0=Re( 4,i)
         Vy0=Re( 5,i)
         Vz0=Re( 6,i)
         ENE00=dsqrt( 1.d0 + Vx0**2+Vy0**2+Vz0**2 )
c
	   IF(shape.eq.1)then
	     CALL gaussian
	   else
	     CALL parax
	   END IF
c     
         VX = VX0+ AVEX*DT1
         VY = VY0+ AVEY*DT1
         VZ = VZ0+ AVEZ*DT1
         gg = DT1 / SQRT( 1.0 + VX*VX + VY*VY + VZ*VZ )
c
         BXT = AVBX*gg
         BYT = AVBY*gg
         BZT = AVBZ*gg
c
         PXS = VX + ( VY*BZT - VZ*BYT )
         PYS = VY + ( VZ*BXT - VX*BZT )
         PZS = VZ + ( VX*BYT - VY*BXT )
c     
         ff  = 2.0/( 1.0+ ( BXT*BXT + BYT*BYT + BZT*BZT ) )
         VX  = VX + AVEX*DT1 + ff*( PYS*BZT - PZS*BYT )
         VY  = VY + AVEY*DT1 + ff*( PZS*BXT - PXS*BZT )
         VZ  = VZ + AVEZ*DT1 + ff*( PXS*BYT - PYS*BXT )
c
         ENE0 =dsqrt( 1.d0 + VX*VX  + VY*VY + VZ*VZ )
         GAMMI= 1.d0/ENE0
         Wx0= Vx*GAMMI
         Wy0= Vy*GAMMI
         Wz0= Vz*GAMMI
c
	   DTI = 1.0d0/DT
         BXT = (Vx-Vx0)*DTI ! fVx 
         BYT = (Vy-Vy0)*DTI ! fVy
         BZT = (Vz-Vz0)*DTI ! fVz
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
         PXS   = (BXT - Wx0*gg)*Alf
         PYS   = (BYT - Wy0*gg)*Alf
         PZS   = (BZT - Wz0*gg)*Alf
         RRR = ENE0**2*(BXT*PXS + BYT*PYS + BZT*PZS)*dt
         RAD = RAD + RRR
         TTT  = AVEX*(Vx*GAMMI)*dt
     &            +AVEY*(Vy*GAMMI)*dt
     &            +AVEZ*(Vz*GAMMI)*dt
c
         Re(1,i) = Xe + Wx0*dt
         Re(2,i) = Ye + Wy0*dt
         Re(3,i) = Ze + Wz0*dt
         Re(4,i) = Vx
         Re(5,i) = Vy
         Re(6,i) = Vz
         Re(7,i) = TTT
         Re(8,i) = (ENE00-ENE)-TTT
         Re(9,i) = RRR !ENE0-ENE-TTT
c
         Re(10,i) = Xi
         Re(11,i) = ENE0
c      
	   end do
      end do     ! end particle loop
!$omp end do

      return
      end
!--------------------------------------
      subroutine Sokolov
!--------------------------------------
      use random_common
      use sim_common
      use mpi_common
	use R_common
      use omp_lib
      implicit none
c
	itotal0 = int(itotal/icpu+1) 
!$omp do            
      do LP = 1,icpu
         IPTSS  = itotal0*(LP-1)+1
         IPTFF  = min(itotal0* LP,itotal)
         do i  = IPTSS,IPTFF      ! particle loop
	   DT1    = PQM*DT*0.5d0
         Xe =Re( 1,i)
         Ye =Re( 2,i)
         Ze =Re( 3,i)
         Vx0=Re( 4,i)
         Vy0=Re( 5,i)
         Vz0=Re( 6,i)
         ENE00=dsqrt( 1.d0 + Vx0**2+Vy0**2+Vz0**2 )
c
	   IF(shape.eq.1)then
	     CALL gaussian
	   else
	     CALL parax
	   END IF
c     
         VX = VX0+ AVEX*DT1
         VY = VY0+ AVEY*DT1
         VZ = VZ0+ AVEZ*DT1
         gg = DT1 / SQRT( 1.0 + VX*VX + VY*VY + VZ*VZ )
c
         BXT = AVBX*gg
         BYT = AVBY*gg
         BZT = AVBZ*gg
c
         PXS = VX + ( VY*BZT - VZ*BYT )
         PYS = VY + ( VZ*BXT - VX*BZT )
         PZS = VZ + ( VX*BYT - VY*BXT )
c     
         ff  = 2.0/( 1.0+ ( BXT*BXT + BYT*BYT + BZT*BZT ) )
         VX  = VX + AVEX*DT1 + ff*( PYS*BZT - PZS*BYT )
         VY  = VY + AVEY*DT1 + ff*( PZS*BXT - PXS*BZT )
         VZ  = VZ + AVEZ*DT1 + ff*( PXS*BYT - PYS*BXT )
c
         ENE0 =dsqrt( 1.d0 + VX*VX  + VY*VY + VZ*VZ )
         GAMMI= 1.d0/ENE0
         Wx0= Vx*GAMMI
         Wy0= Vy*GAMMI
         Wz0= Vz*GAMMI
c
	   DTI = 1.0d0/DT
         BXT = (Vx-Vx0)*DTI ! fVx 
         BYT = (Vy-Vy0)*DTI ! fVy
         BZT = (Vz-Vz0)*DTI ! fVz
c
         Vx0 = Vx
         Vy0 = Vy
         Vz0 = Vz
c------------------------------------------------------------------
         ff = BXT*BXT + BYT*BYT + BZT*BZT 
         gg = BXT*Wx0 + BYT*Wy0 + BZT*Wz0
c
         XI = Alf2*ENE0*sqrt(abs(ff - gg**2))*0.66667d0
         if(XI.gt.0.001d0) then
            kk    = min(idnint(log(XI*1000.d0)*we0i+1.5d0),200)
            ENEh  = Alf*totalRC(kk)
         else
            ENEh  = Alf
         end if
c
         GAMMI = ENEh/(1.0+ENEh*gg)
         PXS   = (BXT - Wx0*gg)*GAMMI
         PYS   = (BYT - Wy0*gg)*GAMMI
         PZS   = (BZT - Wz0*gg)*GAMMI
c
         GAMMI = (BXT*PXS + BYT*PYS + BZT*PZS)*(ENE0**2)
         Vx    = Vx0 - (( PYS*AVBZ - PZS*AVBY) + Wx0*GAMMI)*dt
         Vy    = Vy0 - (( PZS*AVBX - PXS*AVBZ) + Wy0*GAMMI)*dt
         Vz    = Vz0 - (( PXS*AVBY - PYS*AVBX) + Wz0*GAMMI)*dt
c
         ENE   = SQRT( 1.d0 + VX*VX + VY*VY + VZ*VZ )
c
         GAMMI = 1.d0/ENE
c
         TTT   = AVEX*(Vx*GAMMI+PXS)*dt
     &          +AVEY*(Vy*GAMMI+PYS)*dt
     &          +AVEZ*(Vz*GAMMI+PZS)*dt
         RRR   = ENE0**2*(BXT*PXS + BYT*PYS + BZT*PZS)*dt
	   TT1   = TT1 + TTT
         RAD     = RAD + RRR
         Re(1,i) = Xe + Vx*GAMMI*dt + PXS*dt
         Re(2,i) = Ye + Vy*GAMMI*dt + PYS*dt
         Re(3,i) = Ze + Vz*GAMMI*dt + PZS*dt
         Re(4,i) = Vx
         Re(5,i) = Vy
         Re(6,i) = Vz
         Re(7,i) = TTT
         Re(8,i) = (ENE00-ENE)-TTT
         Re(9,i) = RRR !ENE0-ENE-TTT
c
         Re(10,i) = Xi
         Re(11,i) = ENE0
	   end do
	end do
!$omp end do
      return
      end
c
!--------------------------------------
      subroutine Landau_Lifshitz
!--------------------------------------
      use random_common
      use sim_common
      use mpi_common
	use R_common
      use omp_lib
      implicit none
      real(kind=8) :: LLT2X,LLT2Y,LLT2Z,LL3TX,LL3TY,LL3TZ
      real(kind=8) :: BXK,BYK,BZK
c
	itotal0 = int(itotal/icpu+1) 
!$omp do            
      do LP = 1,icpu
         IPTSS  = itotal0*(LP-1)+1
         IPTFF  = min(itotal0* LP,itotal)
         do i  = IPTSS,IPTFF      ! particle loop
	   DT1    = PQM*DT*0.5d0
         Xe =Re( 1,i)
         Ye =Re( 2,i)
         Ze =Re( 3,i)
         Vx0=Re( 4,i)
         Vy0=Re( 5,i)
         Vz0=Re( 6,i)
         ENE00=dsqrt( 1.d0 + Vx0**2+Vy0**2+Vz0**2 )
c
	   IF(shape.eq.0)then
	     CALL gaussian
	   else
	     CALL parax
	   END IF
c     
         VX = VX0+ AVEX*DT1
         VY = VY0+ AVEY*DT1
         VZ = VZ0+ AVEZ*DT1
         gg = DT1 / SQRT( 1.0 + VX*VX + VY*VY + VZ*VZ )
c
         BXT = AVBX*gg
         BYT = AVBY*gg
         BZT = AVBZ*gg
c
         PXS = VX + ( VY*BZT - VZ*BYT )
         PYS = VY + ( VZ*BXT - VX*BZT )
         PZS = VZ + ( VX*BYT - VY*BXT )
c     
         ff  = 2.0/( 1.0+ ( BXT*BXT + BYT*BYT + BZT*BZT ) )
         VX  = VX + AVEX*DT1 + ff*( PYS*BZT - PZS*BYT )
         VY  = VY + AVEY*DT1 + ff*( PZS*BXT - PXS*BZT )
         VZ  = VZ + AVEZ*DT1 + ff*( PXS*BYT - PYS*BXT )
c
         ENE0 =dsqrt( 1.d0 + VX*VX  + VY*VY + VZ*VZ )
         GAMMI= 1.d0/ENE0
         Wx0= Vx*GAMMI
         Wy0= Vy*GAMMI
         Wz0= Vz*GAMMI
c
	   DTI = 1.0d0/DT
         BXT = (Vx-Vx0)*DTI ! fVx 
         BYT = (Vy-Vy0)*DTI ! fVy
         BZT = (Vz-Vz0)*DTI ! fVz
c
         Vx0 = Vx
         Vy0 = Vy
         Vz0 = Vz
c------------------------------------------------------------------
         ff = BXT*BXT + BYT*BYT + BZT*BZT 
         gg = BXT*Wx0 + BYT*Wy0 + BZT*Wz0
c
         XI = Alf2*ENE0*sqrt(abs(ff - gg**2))*0.66667d0
         if(XI.gt.0.001d0) then
            kk    = min(idnint(log(XI*1000.d0)*we0i+1.5d0),200)
            ENEh  = Alf*totalRC(kk)
         else
            ENEh  = Alf
         end if
         GAMMI = ENEh
         PXS   = (BXT - Wx0*gg)*GAMMI
         PYS   = (BYT - Wy0*gg)*GAMMI
         PZS   = (BZT - Wz0*gg)*GAMMI
c
         BXK = BYT*AVBZ-BZT*AVBY
         BYK = - (BXT*AVBZ-BZT*AVBX)
         BZK = BXT*AVBY-BYT*AVBX
c
         LLT2X = ENEh*(BXK+gg*AVEX)
         LLT2Y = ENEh*(BYK+gg*AVEY)
         LLT2Z = ENEh*(BZK+gg*AVEZ)
c
	   GAMMI = (BXT*PXS + BYT*PYS + BZT*PZS)*(ENE0**2)
         LL3TX = Wx0*GAMMI
         LL3TY = Wy0*GAMMI
         LL3TZ = Wz0*GAMMI
c
         Vx= Vx0 + dt*(LLT2X-LL3TX)
         Vy= Vy0 + dt*(LLT2Y-LL3TY)
         Vz= Vz0 + dt*(LLT2Z-LL3TZ)
c
         ENE= SQRT( 1.d0 + VX*VX + VY*VY + VZ*VZ )
         GAMMI= 1.d0/ENE
         Wx= Vx*GAMMI
         Wy= Vy*GAMMI
         Wz= Vz*GAMMI
         TTT = AVEX*Wx + AVEY*Wy + AVEZ*Wz
         RRR     = ENE0**2*(BXT*PXS + BYT*PYS + BZT*PZS)*dt
         TT1	  = TT1 + TTT
         RAD     = RAD + RRR
         Re(1,i) = Xe + Wx*dt
         Re(2,i) = Ye + Wy*dt
         Re(3,i) = Ze + Wz*dt
         Re(4,i) = Vx
         Re(5,i) = Vy
         Re(6,i) = Vz
         Re(7,i) = TTT
         Re(8,i) = (ENE00-ENE)-TTT
         Re(9,i) = RRR
         Re(10,i) = Xi
         Re(11,i) = ENE0
	   end do
	end do
!$omp end do
      return
      end
c
!-------------------------------------
      subroutine random
!-------------------------------------
      use random_common
      use sim_common
      implicit none
      integer, parameter :: values=245
      integer :: d
      integer, dimension(:), allocatable :: seed
c
c-----Dynamic Random Number Generator------
c      call date_and_time(values=values)
      call random_seed(size=d)
      allocate(seed(1:d))
      seed = values
      call random_seed(put=seed)
c
      return
      end
c
!-------------------------------------
      subroutine phemit
!-------------------------------------
      use random_common
      use sim_common
      implicit none
      real(kind=8) :: det1,det2,redcomp,ss
c
      redcomp = comp/(2.d0*PI)
      det1 = 1.d0/(PI*dsqrt(3.d0))
      det2 = 1.d0/(137.d0*redcomp) ! probability coefficient
c
      ff = totalS(kk)*3000
      ss = totalS(kk)*det1*det2*Rt*VL*1.d0/ENE0 !
c
      if(ss.ge.1.d-3) then
        write(*,*) " W*dt > 1 --- Probability of emission > 1"
        write(*,*) " Please reduce time step, "
        write(*,*) " i.e. increase div in SLAC.dat"
        photon = -1
        stop
      end if
c      write(*,*) ss
c
      call random_number(rand)
c      write(*,*) rand
      if(rand.lt.ss) then
	  emmits=.TRUE.
      else
	  emmits=.FALSE.
      end if
c
      return
      end
c
!-------------------------------------
      subroutine qmemit
!-------------------------------------
      use random_common
      use sim_common
      implicit none
      real(kind=8),dimension(6000) :: W
      real(kind=8) :: gg1,ss
c
      call random_number(rand)
c      write(*,*) rand
c
      do k=1,6000
         ii   = idnint((k-0.5d0)*0.5d0 + 0.5d0)
         if(ii.gt.2999) exit
         ss   = (k - 0.5d0)*0.5d0 + 0.5d0 - dble(ii)
         TTT  = diffR(ii,kk) + ss*(diffR(ii+1,kk)-diffR(ii,kk))
         W(k+1)   = W(k) + TTT*0.5d0/ff
         gg1 = W(k+1)
         if(gg1.GE.rand) then
           ENN = k + (rand-W(k))/(W(k+1)-W(k))
         exit
         end if
      end do
      return
      end
c
!-------------------------------------
      subroutine qedemmit
!-------------------------------------
      use sim_common
      use mpi_common
	use R_common
      use omp_lib
      implicit none
c
      if(QED.eq.0) then
        write(*,*) "QED must be turned on"
	  write(*,*) "set QED=1 in SLAC.dat"
        stop
	end if
c
	itotal0 = int(itotal/icpu+1) 
!$omp do            
      do LP = 1,icpu
         IPTSS  = itotal0*(LP-1)+1
         IPTFF  = min(itotal0* LP,itotal)
         do i  = IPTSS,IPTFF      ! particle loop
	   DT1    = PQM*DT*0.5d0
         Xe =Re( 1,i)
         Ye =Re( 2,i)
         Ze =Re( 3,i)
         Vx0=Re( 4,i)
         Vy0=Re( 5,i)
         Vz0=Re( 6,i)
         ENE00=dsqrt( 1.d0 + Vx0**2+Vy0**2+Vz0**2 )
c
	   IF(shape.eq.1)then
	     CALL gaussian
	   else
	     CALL parax
	   END IF
c     
         VX = VX0+ AVEX*DT1
         VY = VY0+ AVEY*DT1
         VZ = VZ0+ AVEZ*DT1
         gg = DT1 / SQRT( 1.0 + VX*VX + VY*VY + VZ*VZ )
c
         BXT = AVBX*gg
         BYT = AVBY*gg
         BZT = AVBZ*gg
c
         PXS = VX + ( VY*BZT - VZ*BYT )
         PYS = VY + ( VZ*BXT - VX*BZT )
         PZS = VZ + ( VX*BYT - VY*BXT )
c     
         ff  = 2.0/( 1.0+ ( BXT*BXT + BYT*BYT + BZT*BZT ) )
         VX  = VX + AVEX*DT1 + ff*( PYS*BZT - PZS*BYT )
         VY  = VY + AVEY*DT1 + ff*( PZS*BXT - PXS*BZT )
         VZ  = VZ + AVEZ*DT1 + ff*( PXS*BYT - PYS*BXT )
c
         ENE0 =dsqrt( 1.d0 + VX*VX  + VY*VY + VZ*VZ )
         GAMMI= 1.d0/ENE0
         Wx0= Vx*GAMMI
         Wy0= Vy*GAMMI
         Wz0= Vz*GAMMI
c
	   DTI = 1.0d0/DT
         BXT = (Vx-Vx0)*DTI ! fVx 
         BYT = (Vy-Vy0)*DTI ! fVy
         BZT = (Vz-Vz0)*DTI ! fVz
c
         Vx0 = Vx
         Vy0 = Vy
         Vz0 = Vz
c------------------------------------------------------------------
         ff = BXT*BXT + BYT*BYT + BZT*BZT 
         gg = BXT*Wx0 + BYT*Wy0 + BZT*Wz0
c
         XI = Alf2*ENE0*sqrt(abs(ff - gg**2))*0.66667d0
         if(XI.gt.0.001d0) then
            kk    = min(idnint(log(XI*1000.d0)*we0i+1.5d0),200)
            ENEh  = Alf*totalRC(kk)
c        calculates whether emission should occur or not
            call phemit
         if(emmits) then
            call qmemit
            ENN = ENN/6000
            write(*,*) "Photon Energy [MeV] =",ENN*ENE0*0.511
c
c  ******* Update electron momentum due to recoil *******
            Vx    = Vx0*(1.d0 - ENN)
            Vy    = Vy0*(1.d0 - ENN)
            Vz    = Vz0*(1.d0 - ENN)
            else
            ENEh  = Alf
            Vx    = Vx0
            Vy    = Vy0
            Vz    = Vz0
            end if
         else
            ENEh  = Alf
            Vx    = Vx0
            Vy    = Vy0
            Vz    = Vz0
         end if
c
         GAMMI = ENEh
         PXS   = (BXT - Wx0*gg)*GAMMI
         PYS   = (BYT - Wy0*gg)*GAMMI
         PZS   = (BZT - Wz0*gg)*GAMMI
	   ENE= SQRT( 1.d0 + VX*VX + VY*VY + VZ*VZ )
         GAMMI= 1.d0/ENE
c
         TTT  = AVEX*(Vx*GAMMI)*dt
     &          +AVEY*(Vy*GAMMI)*dt
     &          +AVEZ*(Vz*GAMMI)*dt
         RRR     = ENE0**2*(BXT*PXS + BYT*PYS + BZT*PZS)*dt
         TT1	  = TT1 + TTT
         RAD     = RAD + RRR
         Re(1,i) = Xe + Vx*GAMMI*dt
         Re(2,i) = Ye + Vy*GAMMI*dt
         Re(3,i) = Ze + Vz*GAMMI*dt
         Re(4,i) = Vx
         Re(5,i) = Vy
         Re(6,i) = Vz
         Re(7,i) = TTT
         Re(8,i) = (ENE00-ENE)-TTT
         Re(9,i) = ENN*ENE0
         Re(10,i) = Xi
         Re(11,i) = ENE0
	   end do
	end do
      return
      end
c
!-------------------------------------
      subroutine histogram
!-------------------------------------
      use sim_common
      use mpi_common
      use R_common
	use out_common
      use omp_lib
      implicit none
	include "mpif.h"
      integer :: enegrid
	real(kind=8) :: b4,b5,b6,enediv,aa,sigma,ffsum
      real(kind=8),dimension(:), allocatable:: his,his_sum
c	
      enegrid = 512
      allocate(his(enegrid))
	allocate(his_sum(enegrid)) 

      sigma=max(sigmax,sigmay,sigmaz)
!$omp workshare
      enediv  = 1/(Em*(1+5.d0*sigma))
	his = 0.d0
	his_sum = 0.d0
!$omp end workshare
c
	itotal0 = int(itotal/icpu)+1
!$omp do            
      do LP = 1,icpu
         IPTSS  = itotal0*(LP-1)+1
         IPTFF  = min(itotal0* LP,itotal)
      do k  = IPTSS,IPTFF      ! particle loop
         b4=Re(4,k) ; b5=Re(5,k) ; b6=Re(6,k)
	   ENE = dsqrt(1.d0 + b4**2 + b5**2 + b6**2)-1.d0
         i = max(idnint(ENE*enediv*enegrid),1)
	   i = min(enegrid,i)
         his(i) = his(i) + wight(k)
	end do
	end do
!$omp end do
c     
      ff = 0.d0 
      do i=1,enegrid
         ff = ff + his(i)
      end do

      call mpi_allreduce(his,his_sum,enegrid,mpi_real8
     &                                 ,mpi_sum,mpi_comm_world,ierr)
	call mpi_allreduce(ff,ffsum,1,mpi_real8
     &                                 ,mpi_sum,mpi_comm_world,ierr)
c
      const = (bin/0.511d6)/(Em*(1+5.d0*sigma)/enegrid)*enum
	if(myrank.eq.0) then
c	  write(*,*) const
        write(fo_name2,444) TRIM(data_file)//'dist_fn'
     &			   ,kstep/100000,jobno
        open (31,file=fo_name2,form='formatted',status='unknown')
	  do i=1,enegrid
           write(31,666) i/(enegrid*enediv)*0.511d6
     &		      , his_sum(i)
        end do
        close(31)
	end if
c
	if(kstep.eq.0) write(9,*) 'Total init. particle', ffsum
	if(kstep.eq.ksmax) write(9,*) 'Total final. particle', ffsum
	deallocate(his,his_sum)

444   format(A,I4.4,I4.4,'.dat')
666   format(2(E14.4,1X))
	return
      end
!-------------------------------------
      subroutine photon_his
!-------------------------------------
      use sim_common
      use mpi_common
      use R_common
	use out_common
      use omp_lib
      implicit none
	include "mpif.h"
      integer :: enegrid
	real(kind=8) :: b4,b5,b6,enediv,aa,sigma,ffsum
      real(kind=8),dimension(:), allocatable:: his,his_sum
c	
      enegrid =512
      allocate(his(enegrid))
	allocate(his_sum(enegrid))

!$omp workshare
      enediv  = 1/Em
      Lall = int(ksmax/icpu)
	his = 0.d0
	his_sum = 0.d0
!$omp end workshare
c
!$omp do            
      do LP=1,icpu                 ! parallelization loop
         IPTSS =    (LP-1)*Lall+1
         IPTFF =     LP   *Lall
      do j=IPTSS,IPTFF           ! timestep loop			
      do L=1,itotal              ! particle loop        
         ENE = phtn(4,j,L)*Pksout
         i = max(idnint(ENE*enediv*enegrid),1)
	   i = min(enegrid,i)
         his(i) = his(i) + wight(L)
	end do
	end do
	end do
!$omp end do
c
      ff = 0.d0 
      do i=1,enegrid
         ff = ff + his(i)
      end do

      call mpi_allreduce(his,his_sum,enegrid,mpi_real8
     &                                 ,mpi_sum,mpi_comm_world,ierr)
c
	call mpi_allreduce(ff,ffsum,1,mpi_real8
     &                                 ,mpi_sum,mpi_comm_world,ierr)
      const = pin/(Em/enegrid)
	if(myrank.eq.0) then
	  write(*,*) const
        write(fo_name2,444) TRIM(data_file)//'dist_ph',jobno
        open (34,file=fo_name2,form='formatted',status='unknown')
	  do i=1,enegrid
	     ENE0 = i*Em/enegrid
           write(34,666) ENE0*0.511d6
     &		      , his_sum(i)*const
        end do
        close(34)
	end if
c
	if(kstep.eq.ksmax) write(9,*) 'Total photon', ffsum
	deallocate(his,his_sum)

444   format(A,I4.4,'.dat')
666   format(2(E14.4,1X))
	return
      end
!-------------------------------------
      subroutine histogram2d
!-------------------------------------
      use sim_common
      use mpi_common
      use R_common
	use out_common
      use omp_lib
      implicit none
	include "mpif.h"
      integer :: grid1, grid2
	real(kind=8) :: b1,b4,b5,b6,momdiv1,momdiv2,aa
      real(kind=8),dimension(:,:), allocatable:: his,his_sum
c	
      grid1 = 64
	grid2 = 64
      allocate(his(grid1,grid2))
	allocate(his_sum(grid1,grid2)) 

	Vx = dsqrt((Em+1.d0)**2-1.d0)
!$omp workshare
      momdiv1 = 1.d0/(0.15*Vx)
      momdiv2 = momdiv1
	his = 0.d0
!$omp end workshare
c
	itotal0 = int(itotal/icpu+1) 
!$omp do            
      do LP = 1,icpu
         IPTSS  = itotal0*(LP-1)+1
         IPTFF  = min(itotal0* LP,itotal)
      do k  = IPTSS,IPTFF      ! particle loop
         b4=Re(4,k) ; b5=Re(5,k) ; b6=Re(6,k)
         b5=-1.d0*b4*dsin(inc_ang) + b5*dcos(inc_ang)	   
         i = max(idnint(b5*momdiv1*32.d0)+32,1)
	   j = max(idnint(b6*momdiv2*32.d0)+32,1)
	   i = min(grid1,i)
	   j = min(grid2,j)
         his(i,j) = his(i,j) + wight(k)
	end do
	end do
!$omp end do
c
	his_sum = 0.d0
      call mpi_allreduce(his,his_sum,grid1*grid2,mpi_real8
     &                                 ,mpi_sum,mpi_comm_world,ierr)
c
      if(myrank==0) then
        write(fo_name2,444) TRIM(data_file)//'dist_fn2d'
     &			   ,kstep/100000,jobno
        open (32,file=fo_name2,form='formatted',status='unknown')
c
	  do i=1,grid1
	  do j=1,grid2
           write(32,666) (i-32)/(32*momdiv1)   !py
     &			,(j-32)/(32*momdiv2)   !pz
     &			,his_sum(i,j)
        end do
	  end do

	  close(32)
	end if
c
	deallocate(his,his_sum)

444   format(A,I4.4,I4.4,'.dat')
666   format(3(E14.4,1X))
	return
      end
!--------------------------------------
	subroutine manual_load
!--------------------------------------
      use random_common
      use sim_common
      implicit none
      real(kind=8),dimension(6000) :: W
      real(kind=8) :: gg1,ss,energy
      INTEGER, PARAMETER :: np_local = 116
      INTEGER :: ip
      real(kind=8), DIMENSION(np_local) :: Ex_axis, dist_fn

      OPEN(33,file='f_E_smoothed_new_final.txt',status='OLD')

      DO ip = 1, np_local

        READ(33,*) Ex_axis(ip), dist_fn(ip)

        Ex_axis(ip) = Ex_axis(ip)/0.511
        dist_fn(ip) = dist_fn(ip)*5

      ENDDO

      CLOSE(33)
c
	call random_number(rand)
c
	W = 0.d0
      do k=1,3000
         ii   = idnint((k-0.5d0)*116.d0/3000.d0 + 0.5d0)
         if(ii.gt.115) exit
         ss   = (k - 0.5d0)*116.d0/3000.d0 + 0.5d0 - dble(ii)
         TTT  = dist_fn(ii) + ss*(dist_fn(ii+1)-dist_fn(ii))
         W(k+1)   = W(k) + TTT*116.d0/3000.d0
         gg1 = W(k+1)
         if(rand.lt.gg1) then
           ENN = k + (rand-W(k))/(W(k+1)-W(k))
         exit
         end if
      end do

	energy = ENN*2300/0.511/3000
	p_x = dsqrt((energy+1.d0)**2-1.d0) 
	Vx  = p_x
	return
      end

C-----------------------------------------------------------------
      subroutine pair_init
c
c cross section of Bremmstrahlung & pair production (Bethe-Heitler)
c
C-----------------------------------------------------------------
      use random_common
	use mpi_common
	use out_common
      implicit none
      integer :: i,j
      real(kind=8) :: pi,pi2
      real(kind=8) :: totalmin,totalmin2,corrct
      real(kind=8) :: kx,k0,dk,kk,k,E0,p0,E,p,eL,ec,ee,ppi,p3i,p3j
      real(kind=8) :: D1,D2,D3,D4,D5,Fd,Fd0,Fd1,v0,v1,cor0,cor1,cor2
      real(kind=8) :: ele,red,totL
c
      real(kind=8),dimension(LE0,LPx) :: resultH
      real(kind=8),dimension(LE0,LPx) :: resultI
      real(kind=8),dimension(LE0,LPx) :: resultJ
      real(kind=8),dimension(LE0,LPx) :: resultK
c
      PI  = 4.D0*DATAN( 1.D0 )
      PI2 = PI*2.D0
c
	if(myrank.eq.0) write(*,*) "Bethe-Heitler"
c
      totalmin = log(183.d0/Zcm3)
      totalmin2= 0.925d0*(Zcom/137.d0)**2 ! for high Z
c     totalmin2= 1.200d0*(Zcom/137.d0)**2 ! for  low Z
      dp1=(log(Emax)-log(2.0d0))/BPx  !! 
      corrct  =0.6d0/82.d0*Zcom

      do j=1,LPx
         kx = dexp(dp1*dble(j-1))*2.d0  ! initial photonn energy in x
         k0 = dsqrt(kx**2)
         dk =(k0-2.d0)/dble(LE0)       ! minus 2*mc^2
         kk =1.d0/(k0*k0*k0)
         if(k0.le.2.d0) cycle
         do i=1,LE0
            k=dk*(dble(i)-0.5d0) ! ratio of emission positron energy
            E0=k+1.d0     ! positron energy
            p0=sqrt(E0**2-1.d0) ! positron momentum

            E =k0-E0
            p =sqrt( E**2-1.d0)

            eL =2.d0*log((E*E0+p*p0+1.d0)/k0) ! L
            ec =2.d0*log(E0+P0)               ! E+
            ee =2.d0*log(E +p )               ! E-

            ppi= 1.d0/(p0*p) 
            p3i= 1.d0/(p0**3)
            p3j= 1.d0/(p**3)

            D1 =-4.d0/3.d0 - 2.d0*E*E0*(p**2+p0**2)*ppi**2
            D2 = ec*E*p3i+ee*E0*p3j-ee*ec*ppi
            D3 = (8.d0/3.d0*E*E0*ppi)
            D4 = (k0**2*(ppi**3))*((E0*E)**2+(p0*p)**2)
            D5 = (0.5d0*k0*ppi)*((E0*E-p0**2)*ec*p3i
     &                       +(E0*E- p**2)*ee*p3j+2.d0*k0*E0*E*(ppi**2))

            Fd = (p0*p*kk*(D1+D2+eL*(D4-D3-D5)))*(k0-2.d0)

            resultH(i,j)=Fd
c 
c     relativistic & screening
c
            eL =(E**2+E0**2+(2.d0/3.d0)*E*E0)*totalmin
            ec = E*E0/9.d0
            Fd0= 4.d0*kk*(eL-ec)*(k0-2.d0)
c
            resultI(i,j)=Fd0
c
c     Coulomb correction for relativistic case
c
            eL =(E**2+E0**2+(2.d0/3.d0)*E*E0)*totalmin2
            Fd1= -2.d0*kk*eL*(k0-2.d0)

            resultJ(i,j)=Fd1+Fd0
c         
c     Coulomb correction for non relativistic case
c
            v0 = p0/dsqrt(1.d0+p0**2)
            v1 =  p/dsqrt(1.d0+ p**2)
            cor0 = corrct/v0
            cor1 = corrct/v1
            cor2 =pi2*pi2*cor0*cor1/(dexp( pi2*cor0)-1.d0)
     &                          /(1.d0-dexp(-pi2*cor1))
c
            resultK(i,j)=Fd*cor2
        end do
      end do
c
      resultL=0.d0
      do j=1,LPx
         kx = dexp(dp1*dble(j-1))*2.d0  ! initial photonn energy in x
         k0 = dsqrt(kx**2)
         dk =(k0-2.d0)/dble(LE0)       ! minus 2*mc^2
         kk =1.d0/(k0*k0*k0)
c
         do i=1,LE0
c
            k=dk*(dble(i)-0.5d0) ! ratio of emission positron energy
            E0=k+1.d0     ! positron energy
            p0=sqrt(E0**2-1.d0) ! positron momentum

            E =k0-E0
            p =sqrt( E**2-1.d0)

            eL = min(E,E0)
            ee = E*E0/k0
            if((ee.gt.8.d0).and.(ee.lt.24.d0)) then 
              ele =(ee-8.d0)/16.d0
              red = (1.d0-ele)*resultH(i,j) + ele*resultJ(i,j)
              resultL(i,j)=red
            else if(ee.ge.24.d0) then 
              resultL(i,j)=resultJ(i,j)
            else if((eL.gt.0.1d0).and.(eL.lt.0.6d0)) then 
              ele =(eL-0.1d0)/0.5d0
              red = (1.d0-ele)*resultK(i,j) + ele*resultH(i,j)
              resultL(i,j)=red
            else if(eL.le.0.1d0) then 
              resultL(i,j)=resultK(i,j)
            else  
              resultL(i,j)=resultH(i,j)
            end if    
         end do
      end do
c
      do j=1,LPx
         kx = dexp(dp1*dble(j-1))*2.d0  ! initial photonn energy in x
         k0 = dsqrt(kx**2)
         dk =(k0-2.d0)/dble(LE0)       ! minus 2*mc^2
         totL=0.d0
         do i=1,LE0
            if(k0.gt.2.d0) then 
              k=dk/(k0-2.d0)
              totL = totL+resultL(i,j)*k
            end if
         end do
         totalH(j)=totL
      end do
c
      do j=1,LPx
         totL=0.d0
         do i=1,LE0 
            totL = totL+resultL(i,j)/dble(LE0)
         end do
         totalH2(j)=totL
      end do
	
	if(myrank==0) then
        write(fo_name2,444) TRIM(data_file)//'Bethe-Heitler',jobno
        open (33,file=fo_name2,form='formatted',status='unknown')
c
	  do i=1,LPx
           write(33,*) i, totalH(i), totalH2(i)
	  end do
	  close(33)
	end if 
444   format(A,I4.4,'.dat')
      return 
      end
C-----------------------------------------------------------------
      subroutine welcome
C-----------------------------------------------------------------
	use mpi_common
	use sim_common
      if(myrank.ne.0) return
	write(*,*) ' '//achar(27)//'[34m'
	write(*,*) ' '//achar(27)//'[1m'
	write(*,*) '	###L      ########E    #######C      '
	write(*,*) '	###L      ########E   ##########C  	 '   
	write(*,*) '	###L      ###E       ###C    ###C    '
	write(*,*) '	###L      ########E  ###C            '
	write(*,*) '	###L      ########E  ###C            '
	write(*,*) '	###L      ###E       ###C    ###C    '
	write(*,*) '	#######L  ########E   ##########C    '
	write(*,*) '	LASER###  ELECTRON#    COLLISION     '
	write(*,*) ' '//achar(27)//'[32m'
	write(*,*) 'Welcome to Laser Electron Collision code (v-1.3.0)'
	write(*,*) ' '//achar(27)//'[0m'
      write(*,*) '*****************************************************'
      write(*,*) "The code is running on",NPROCS," processors"
      write(*,*) "		       ",threads," threads"
      write(*,*) '*****************************************************'
	if(load.eq.1) 
     &	write(*,*) "Load particles from: f_E_smoothed_new_final.txt"
      if(OutRad.eq.1  ) write(*,*) 'Produce radiation'
	if(OutPairs.eq.1) write(*,*) 'Produce pairs'
      write(*,*)
      return 
      end
