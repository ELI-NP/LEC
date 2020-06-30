C----------------------------------------------------------
      SUBROUTINE parax
C----------------------------------------------------------
c
      ! Laser field beyond paraxial approximation up to fifth order
c
      use mpi_common
      use random_common
      use sim_common
      implicit none
	real(kind=8):: xr,eps,r,xxi,nu,zeta,eta,rho,w,g,eta0
      real(kind=8):: PsiP,PsiG,PsiR,Psi0,Psi
      real(kind=8):: S0,S2,S3,S4,S5,S6,C1,C2,C3,C4,C5,C6,C7,C8
	real(kind=8):: TermE1,TermE2,TermE3,TermE4,TermE5
	real(kind=8):: TermB1,TermB2,TermB3,TermB4
	real(kind=8):: E1temp,E2temp,E3temp,B1temp,B2temp,B3temp
c
      xr=wk*ws*ws/2.d0
      eps=ws/xr
c
      r=dsqrt(Ye**2+Ze**2)
      xxi=Ye*wsi
      nu=Ze*wsi
      zeta=Xe/xr
c
      rho=r/ws
c
      eta0=0.d0

	Wx0  = (T-Xe)
c
      w = ws*sqrt(1.d0+Xe*Xe/(xr*xr))
c
	 g = 1.d0/cosh(Wx0*wpi)
c
      PsiP = (T-Xe)*Wk;
      PsiG = atan(zeta);
      PsiR = 0.5d0*wk*Xe*r*r/(Xe*Xe+xr*xr)
      Psi0 = 0.0d0;
      Psi = Psi0 + PsiP - PsiR + PsiG;
c
      BXT=ws/w*g*exp(-r*r/(w*w))
c
      S0=sin(Psi)
      S2=(ws/w)**2.d0*sin(Psi+2.d0*PsiG)
      S3=(ws/w)**3.d0*sin(Psi+3.d0*PsiG)
      S4=(ws/w)**4.d0*sin(Psi+4.d0*PsiG)
      S5=(ws/w)**5.d0*sin(Psi+5.d0*PsiG)
      S6=(ws/w)**6.d0*sin(Psi+6.d0*PsiG)
c
      C1=(ws/w)*cos(Psi+PsiG)
      C2=(ws/w)**2.d0*cos(Psi+2.d0*PsiG)
      C3=(ws/w)**3.d0*cos(Psi+3.d0*PsiG)
      C4=(ws/w)**4.d0*cos(Psi+4.d0*PsiG)
      C5=(ws/w)**5.d0*cos(Psi+5.d0*PsiG)
      C6=(ws/w)**6.d0*cos(Psi+6.d0*PsiG)
      C7=(ws/w)**7.d0*cos(Psi+7.d0*PsiG)
c
	TermE1 = xxi*xxi*S2-rho**4.d0*S3/4.d0
      TermE2 = S2/8.d0-rho*rho*S3/4.d0
     &         -rho*rho*(rho*rho-16.d0*xxi*xxi)*S4/16.d0
     &         -rho**4.d0*(rho*rho+2.d0*xxi*xxi)*S5/8.d0
     &         +rho**8.d0*S6/32.d0
	TermE3 = rho**2.d0*S4-rho**4.d0*S5/4.d0
	TermE4 = -C2/2.d0+rho*rho*C3-rho**4.d0*C4/4.d0
	TermE5 = -3.d0*C3/8.d0-3.d0*rho*rho*C4/8.d0
     &	   + 17.d0*rho**4.d0*C5/16.d0-3.d0*rho**6.d0*C6/8.d0
     &         +rho**8.d0*C7/32.d0
	TermB1 = rho*rho*S2/2.d0-rho**4.d0*S3/4.d0
	TermB2 = -S2/8.d0+rho*rho*S3/4.d0+5.d0*rho**4.d0*S4/16.d0
     &	   -rho**6.d0*S5/4.d0+rho**8.d0*S6/32.d0
	TermB3 = C2/2.d0+rho**2.d0*C3/2.d0-rho**4.d0*C4/4.d0
	TermB4 = 3.d0*C3/8.d0+3.d0*rho**2.d0*C4/8.d0
     &	   + 3*rho**4.d0*C5/16.d0-rho**6.d0*C6/4.d0
     &	   + rho**8.d0*C7/32.d0
c
	E1temp= xxi*(eps*C1+eps**3.d0*(TermE4)+eps**5.d0*(TermE5))
c
      E2temp= S0+eps**2.d0*(TermE1)+eps**4.d0*(TermE2)
c
      E3temp= xxi*nu*(eps**2.d0*S2+eps**4.d0*(TermE3))
c
	B1temp=nu*(eps*C1+eps**3.d0*(TermB3)+eps**5.d0*(TermB4))
c
      B2temp=0.d0
c
      B3temp=S0+eps*eps*(TermB1)+eps**4.d0*(TermB2)
c
      AVEX = BXT*E1temp
      AVEY = BXT*E2temp
      AVEZ = BXT*E3temp
      AVBX = BXT*B1temp
      AVBY = BXT*B2temp
      AVBZ = BXT*B3temp
c
      RETURN
      END

C----------------------------------------------------------
	SUBROUTINE gaussian
C----------------------------------------------------------

      ! Spatial and temporal Gaussian shape laser pulse
c
      use random_common
      use sim_common
      implicit none
c---------------------------------------------------------------------
c for integer time step, CEx,y,z,  PEx,y,z and PBx,y,z
c---------------------------------------------------------------------
c
	Wx0  = (T-Xe)		      ! wk = normalized wave number
c
	ff   = (Ye*wsi)**2+(Ze*wsi)**2
	BXT  = dexp(-((Wx0*wpi)**2+ff)) ! initial amplitude of pulse laser
c
	Wx0  = Wx0*wk
c
      BYT  = dsin( Wx0)*BXT
      BZT  = dcos( Wx0)*BXT
c
      AVEX = 0.d0
      AVEY = BYT 
      AVEZ = polar*BZT
      AVBX = 0.d0
      AVBY = -1.d0*AVEZ
      AVBZ = AVEY 
c
      RETURN
      END

