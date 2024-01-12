***********************************************

	function staustau2QQ_M2Vertex()
	implicit none

#include "staustau2QQ_Kinematics.h"
#include "staustau2QQ_Model.h"
#include "staustau2QQ_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

* functions & subroutines
	integer a,b,e,d,i,j,m,n,p,q,x,y
	double precision colfac,cfactor
	double precision findValW
	complex*16 staustau2QQ_M2Vertex,test,Zgluino1l,Zgluino2l,Zgluino3l,higgsgluinocorr(2),higgsgluoncorr(2),Zgluon1l,Zgluon2l,Zgluon3l,DenV(3),DenSH(2)
	complex*16 Zgluino1r,Zgluino2r,Zgluino3r,Zgluon1r,Zgluon2r,Zgluon3r
      complex*16 photongluino1r,photongluino2r,photongluino3r,photongluon1r,photongluon2r,photongluon3r
      complex*16 photongluino1l,photongluino2l,photongluino3l,photongluon1l,photongluon2l,photongluon3l
	double precision g0R,g0L,g1L,g2L,g1R,g2R,g0,M0,M1,M2,k0,k0L,k0R,k1L,k1R,k2L,k2R

	external myCget, myBget

	staustau2QQ_M2Vertex = 0d0

      DenV(1)=0d0
	DenV(2) = sqrtS**2 
	DenV(3) = (sqrtS**2-MZ**2)
	
      do x=1,2
	  DenSH(x) = (sqrtS**2-MHk(x)**2)	
      enddo

!      myCget(55248.057920899111195,59337.89087323519925,434232.05987013,0d0,59337.89087323519925,0d0)
!      

      call myCget((sqrtS*muu3)**2,(sqrtS)**2,(sqrtS*muu4)**2,0d0,(sqrtS*muu3)**2,(sqrtS*muu3)**2)
      call myBget((sqrtS)**2,(sqrtS*muu3)**2,(sqrtS*muu3)**2)

 !************************************************
  !!-----------------Higgs with gluon exch.------ 
  !********************************************** 
  
      higgsgluoncorr(1)=0d0
      higgsgluoncorr(2)=0d0
      
	cfactor = 4d0/3d0
	
      do x=1,2
      
	higgsgluoncorr(x) = cfactor*GS**2*hff(x)/(16d0*pi**2)*(4*xb0 + (4*(sqrtS*muu3)**2-2*(sqrtS)**2)*(xc1+xc2)-2*sqrtS**2*(1-muu3**2-muu4**2)*xc0)
	
      enddo

      
 !************************************************
!!-----------------Higgs with gluino exch.---------  (IR (massive gluino) and UV convergent)
!*************************************************

      higgsgluinocorr(1) = 0d0
      higgsgluinocorr(2) = 0d0
     
      Do i=1,3
      Do j=1,3
      Do a=1,2
      Do b=3,4
      Do e=1,2
      Do d=3,4

      cfactor = 4d0/3d0
      M0= MGl
      M1= MSf(a,b,i)
      M2= MSf(e,d,j)
      
      k0 =GLR(x,a,e,b,i)
      k1L=GluL(a,b,3,i,3)
      k1R=GluR(a,b,3,i,3)
      k2L=GluR(e,d,3,j,3)
      k2R=GluL(e,d,3,j,3)
      
      call myCget((sqrtS*muu3)**2,(sqrtS)**2,(sqrtS*muu4)**2,M0**2,M1**2,M2**2)

      higgsgluinocorr(x) = higgsgluinocorr(x) + cfactor/(16d0*pi**2)*(-k0)*(M0*xc0*k1L*k2L-muu4*sqrtS*xc1*k1R*k2L-muu3*sqrtS*xc2*k1L*k2R)
      
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
  
! !    ************************************************
! !   -----------------Z/gamma with gluon exch.------- 
! !   ********************************************** 
! !   
! !     we treat Z and gamma exchange separately for more checks

      
      photongluon1l=0d0
      photongluon2l=0d0
      photongluon3l=0d0
      photongluon1r=0d0
      photongluon2r=0d0
      photongluon3r=0d0
! 
      call myBget(sqrtS**2,(sqrtS*muu3)**2,(sqrtS*muu4)**2) 
      call myCget((sqrtS*muu3)**2,(sqrtS)**2,(sqrtS*muu4)**2,0d0,(sqrtS*muu4)**2,(sqrtS*muu3)**2)
      
 
      cfactor = 4d0/3d0
      M0= 0d0
      M1= (sqrtS*muu3)
      M2= (sqrtS*muu4)
      
      k0L=BL(2,3,3,3,3)
      k0R=BR(2,3,3,3,3)
      
      
      photongluon1l=photongluon1l + cfactor*GS**2/(16d0*pi**2)*(-2d0*(-xb0+2*xc00-
     -   (sqrtS*muu4)**2*xc1-(sqrtS*muu3)**2*xc2)*k0L-sqrtS**2*(1d0-muu3**2-muu4**2)*(xc0+xc1+xc2)*k0L+2*(xc1+xc2)*k0R*(sqrtS*muu3)*(sqrtS*muu4))
      
      photongluon2l=photongluon2l + cfactor*GS**2/(16d0*pi**2)*(2*M1*k0R*(xc1-xc2-xc0)+2*M2*(xc0-xc2+xc1)*k0L+
     -   (sqrtS*muu4)*(-2*(xc11-xc12+xc1)+2*(xc0+xc1+xc2))*k0R+(sqrtS*muu3)*(-2*(xc12-xc22-xc2)-2*(xc0+xc1+xc2))*k0L) 
      
      photongluon3l=photongluon3l + cfactor/(16d0*pi**2)*(2*M1*(xc1+xc0+xc2)*k0R+2*M2*k0L*(xc0+xc1+xc2)+
     -   (sqrtS*muu4)*(-2*(xc11+xc12+xc1)-2*(xc0+xc1+xc2))*k0R+(sqrtS*muu3)*(-2*(xc12+xc22+xc2)-2*(xc0+xc1+xc2))*k0L)      
     

      photongluon1r=photongluon1r + cfactor*GS**2/(16d0*pi**2)*(-2d0*(-xb0+2*xc00-
     -   (sqrtS*muu4)**2*xc1-(sqrtS*muu3)**2*xc2)*k0R-sqrtS**2*(1-muu3**2-muu4**2)*(xc0+xc1+xc2)*k0R+2*(xc1+xc2)*k0L*(sqrtS*muu3)*(sqrtS*muu4))
      
      photongluon2r=photongluon2r + cfactor*GS**2/(16d0*pi**2)*(2*M1*k0L*(xc1-xc2-xc0)+2*M2*(xc0-xc2+xc1)*k0R+
     -   (sqrtS*muu4)*(-2*(xc11-xc12+xc1)+2*(xc0+xc1+xc2))*k0L+(sqrtS*muu3)*(-2*(xc12-xc22-xc2)-2*(xc0+xc1+xc2))*k0R) 
      
      photongluon3r=photongluon3r + cfactor/(16d0*pi**2)*(2*M1*(xc1+xc0+xc2)*k0L+2*M2*k0R*(xc0+xc1+xc2)+
     -   (sqrtS*muu4)*(-2*(xc11+xc12+xc1)-2*(xc0+xc1+xc2))*k0L+(sqrtS*muu3)*(-2*(xc12+xc22+xc2)-2*(xc0+xc1+xc2))*k0R)     

      Zgluon1l=0d0
      Zgluon2l=0d0
      Zgluon3l=0d0
      Zgluon1r=0d0
      Zgluon2r=0d0
      Zgluon3r=0d0

      cfactor = 4d0/3d0
      M0= 0d0
      M1= (sqrtS*muu3)
      M2= (sqrtS*muu4)
      
      k0L=BL(3,3,3,3,3)
      k0R=BR(3,3,3,3,3)

      
      Zgluon1l=Zgluon1l + cfactor*GS**2/(16d0*pi**2)*(-2d0*M1*M2*xc0*k0R-2d0*(-xb0+2*xc00-
     -   (sqrtS*muu4)**2*xc1-(sqrtS*muu3)**2*xc2)*k0L-sqrtS**2*(1-muu3**2-muu4**2)*(xc0+xc1+xc2)*k0L+2*(xc0+xc1+xc2)*k0R*(sqrtS*muu3)*(sqrtS*muu4))
      
      Zgluon2l=Zgluon2l + cfactor*GS**2/(16d0*pi**2)*(2*M1*k0R*(xc1-xc2-xc0)+2*M2*(xc0-xc2+xc1)*k0L+
     -   (sqrtS*muu4)*(-2*(xc11-xc12+xc1)+2*(xc0+xc1+xc2))*k0R+(sqrtS*muu3)*(-2*(xc12-xc22-xc2)-2*(xc0+xc1+xc2))*k0L) 
      
      Zgluon3l=Zgluon3l + cfactor/(16d0*pi**2)*(2*M1*(xc1+xc0+xc2)*k0R+2*M2*k0L*(xc0+xc1+xc2)+
     -   (sqrtS*muu4)*(-2*(xc11+xc12+xc1)-2*(xc0+xc1+xc2))*k0R+(sqrtS*muu3)*(-2*(xc12+xc22+xc2)-2*(xc0+xc1+xc2))*k0L)      
     
     
      Zgluon1r=Zgluon1r + cfactor*GS**2/(16d0*pi**2)*(-2d0*M1*M2*xc0*k0L-2d0*(-xb0+2*xc00-
     -   (sqrtS*muu4)**2*xc1-(sqrtS*muu3)**2*xc2)*k0R-sqrtS**2*(1-muu3**2-muu4**2)*(xc0+xc1+xc2)*k0R+2*(xc0+xc1+xc2)*k0L*(sqrtS*muu3)*(sqrtS*muu4))
      
      Zgluon2r=Zgluon2r + cfactor*GS**2/(16d0*pi**2)*(2*M1*k0L*(xc1-xc2-xc0)+2*M2*(xc0-xc2+xc1)*k0R+
     -   (sqrtS*muu4)*(-2*(xc11-xc12+xc1)+2*(xc0+xc1+xc2))*k0L+(sqrtS*muu3)*(-2*(xc12-xc22-xc2)-2*(xc0+xc1+xc2))*k0R) 
      
      Zgluon3r=Zgluon3r + cfactor/(16d0*pi**2)*(2*M1*(xc1+xc0+xc2)*k0L+2*M2*k0R*(xc0+xc1+xc2)+
     -   (sqrtS*muu4)*(-2*(xc11+xc12+xc1)-2*(xc0+xc1+xc2))*k0L+(sqrtS*muu3)*(-2*(xc12+xc22+xc2)-2*(xc0+xc1+xc2))*k0R)     
  
!         
!    !************************************************
!   !!-----------------Z/gamma with gluino exch.-------  (IR convergent)
!   !********************************************** 
!   
!     !gamma propagator

      photongluino1l=0d0
      photongluino2l=0d0
      photongluino3l=0d0
      photongluino1r=0d0
      photongluino2r=0d0
      photongluino3r=0d0
        
      Do i=1,3
      Do j=1,3
      Do a=1,2
      Do b=3,4
      Do e=1,2
      Do d=3,4
        
      cfactor = 4d0/3d0
      M0= MGl
      M1= MSf(a,b,i)
      M2= MSf(e,d,j)
      
      k0=VSfSf(2,a,e,b,d,i,j)
      k1L=GluL(a,b,3,i,3)
      k1R=GluR(a,b,3,i,3)
      k2L=GluR(e,d,3,j,3)
      k2R=GluL(e,d,3,j,3)
      
      call myCget((sqrtS*muu3)**2,(sqrtS)**2,(sqrtS*muu4)**2,M0**2,M1**2,M2**2)

      photongluino1l=photongluino1l - cfactor*k0/(16d0*pi**2)*(-2*xc00*k2R*k1L)
      
      photongluino2l=photongluino2l - cfactor*k0/(16d0*pi**2)*(M0*(xc1-xc2)*k1L*k2L-(sqrtS*muu4)*(xc11-xc12)*k2L*k1R+(sqrtS*muu3)*(xc12-xc22)*k2R*k1L)      
      
      photongluino3l=photongluino3l - cfactor*k0/(16d0*pi**2)*(M0*(xc0+xc1+xc2)*k1L*k2L-(sqrtS*muu4)*k2L*k1R*(xc1+xc11+xc12)+(sqrtS*muu3)*k2R*k1L*(xc2+xc12+xc22))
      
      photongluino1r=photongluino1r - cfactor*k0/(16d0*pi**2)*(-2*xc00*k2l*k1r)
      
      photongluino2r=photongluino2r-cfactor*k0/(16d0*pi**2)*(M0*(xc1-xc2)*k1R*k2R-(sqrtS*muu4)*(xc11-xc12)*k2R*k1L+(sqrtS*muu3)*(xc12-xc22)*k2L*k1R)      
      
      photongluino3r=photongluino3r - cfactor*k0/(16d0*pi**2)*(M0*(xc0+xc1+xc2)*k1R*k2R-(sqrtS*muu4)*k2R*k1L*(xc1+xc11+xc12)+(sqrtS*muu3)*k2L*k1R*(xc2+xc12+xc22))       

      enddo  
      enddo
      enddo
      enddo
      enddo
      enddo
      
!       !Z propagator

      Zgluino1l=0d0
      Zgluino2l=0d0
      Zgluino3l=0d0
      Zgluino1r=0d0
      Zgluino2r=0d0
      Zgluino3r=0d0
! !       
      Do i=1,3
      Do j=1,3
      Do a=1,2
      Do b=3,4
      Do e=1,2
      Do d=3,4
        
      cfactor = 4d0/3d0
      M0= MGl
      M1= MSf(a,b,i)
      M2= MSf(e,d,j)
      
      k0=VSfSf(3,a,e,b,d,i,j)
      k1L=GluL(a,b,3,i,3)
      k1R=GluR(a,b,3,i,3)
      k2L=GluR(e,d,3,j,3)
      k2R=GluL(e,d,3,j,3)
      
      call myCget((sqrtS*muu3)**2,(sqrtS)**2,(sqrtS*muu4)**2,M0**2,M1**2,M2**2)

      Zgluino1l=Zgluino1l - cfactor*k0/(16d0*pi**2)*(-2*xc00*k2R*k1L)
      
      Zgluino2l=Zgluino2l - cfactor*k0/(16d0*pi**2)*(M0*(xc1-xc2)*k1L*k2L-(sqrtS*muu4)*(xc11-xc12)*k2L*k1R+(sqrtS*muu3)*(xc12-xc22)*k2R*k1L)      
      
      Zgluino3l=Zgluino3l - cfactor*k0/(16d0*pi**2)*(M0*(xc0+xc1+xc2)*k1L*k2L-(sqrtS*muu4)*k2L*k1R*(xc1+xc11+xc12)+(sqrtS*muu3)*k2R*k1L*(xc2+xc12+xc22))
      
      Zgluino1r=Zgluino1r - cfactor*k0/(16d0*pi**2)*(-2*xc00*k2l*k1r)
      
      Zgluino2r=Zgluino2r - cfactor*k0/(16d0*pi**2)*(M0*(xc1-xc2)*k1R*k2R-(sqrtS*muu4)*(xc11-xc12)*k2R*k1L+(sqrtS*muu3)*(xc12-xc22)*k2L*k1R)      
      
      Zgluino3r=Zgluino3r - cfactor*k0/(16d0*pi**2)*(M0*(xc0+xc1+xc2)*k1R*k2R-(sqrtS*muu4)*k2R*k1L*(xc1+xc11+xc12)+(sqrtS*muu3)*k2L*k1R*(xc2+xc12+xc22))       

      enddo  
      enddo
      enddo
      enddo
      enddo
      enddo


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  !************************************************
  !!-----------------INTERFERENCES------- 
  !********************************************** 
  
    !!---------------------------------Higgs-2-quark vertex-------------------------------- 
  
    !contributions only proportional to tree level
   
      colfac = Nc
  
      Do x=1,2
      Do y=1,2
 
      staustau2QQ_M2Vertex =  staustau2QQ_M2Vertex + 2d0*colfac*hSf(x)*hSf(y)*(higgsgluinocorr(x)+higgsgluoncorr(x))*hff(y)/(DenSH(x)*DenSH(y))*
     -    (-2*(-1 + muu3**2 + 2*muu3*muu4 + muu4**2)*sqrtS**2) 
     
       enddo
      enddo  
!   
!    !!---------------------------------Vector-2-quark vertex-------------------------------- 
! 
! 
       !Tree level proportional:
     
       !There are further kinematical structures that do not factorise with the Born CS
       !the structure are named like in the HitchHiker nomenclature
       !structures B1...,C1...,D1... vanish due to on-shell condition + Dirac equation
! !  
	Do y=2,3
	
	colfac=Nc
	g1L=photongluon1l+photongluino1l
	g2L=BL(y,3,3,3,3)
	g1R=photongluon1r+photongluino1r
	g2R=BR(y,3,3,3,3)
	
      staustau2QQ_M2Vertex =  staustau2QQ_M2Vertex +2d0*colfac*(sqrtS**4*(g1L*(-2*g2R*(-1 + 2*muu1**2 + 2*muu2**2)*muu3*muu4 +
     -   g2L*(2 - 3*muu3**2 + muu3**4 - 3*muu4**2 +
     -   2*muu3**2*muu4**2 + muu4**4 +
     -   4*muu2**2*(-1 + muu3**2 + muu4**2 - tred) + 
     -   4*muu1**2*(-1 + muu2**2 + muu3**2 + muu4**2 - tred) +
     -   4*tred - 4*muu3**2*tred - 4*muu4**2*tred + 4*tred**2)) +
     -    g1R*(-2*g2L*(-1 + 2*muu1**2 + 2*muu2**2)*muu3*muu4 +
     -   g2R*(2 - 3*muu3**2 + muu3**4 - 3*muu4**2 +
     -   2*muu3**2*muu4**2 + muu4**4 +
     -   4*muu2**2*(-1 + muu3**2 + muu4**2 - tred) +
     -   4*muu1**2*(-1 + muu2**2 + muu3**2 + muu4**2 - tred) +
     -   4*tred - 4*muu3**2*tred - 4*muu4**2*tred + 4*tred**2))))*
     -   VSfSf(2,isf1,isf2,itt1,itt2,igen1,igen2)*VSfSf(y,isf1,isf2,itt1,itt2,igen1,igen2)/(DenV(2)*DenV(y))
	
	enddo
! 
	Do y=2,3
	
	colfac=Nc
	g1L=Zgluon1l+Zgluino1l
	g2L=BL(y,3,3,3,3)
	g1R=Zgluon1r+Zgluino1r
	g2R=BR(y,3,3,3,3)
	
      staustau2QQ_M2Vertex =  staustau2QQ_M2Vertex +2d0*colfac*(sqrtS**4*(g1L*(-2*g2R*(-1 + 2*muu1**2 + 2*muu2**2)*muu3*muu4 +
     -   g2L*(2 - 3*muu3**2 + muu3**4 - 3*muu4**2 +
     -   2*muu3**2*muu4**2 + muu4**4 +
     -   4*muu2**2*(-1 + muu3**2 + muu4**2 - tred) + 
     -   4*muu1**2*(-1 + muu2**2 + muu3**2 + muu4**2 - tred) +
     -   4*tred - 4*muu3**2*tred - 4*muu4**2*tred + 4*tred**2)) +
     -    g1R*(-2*g2L*(-1 + 2*muu1**2 + 2*muu2**2)*muu3*muu4 +
     -   g2R*(2 - 3*muu3**2 + muu3**4 - 3*muu4**2 +
     -   2*muu3**2*muu4**2 + muu4**4 +
     -   4*muu2**2*(-1 + muu3**2 + muu4**2 - tred) +
     -   4*muu1**2*(-1 + muu2**2 + muu3**2 + muu4**2 - tred) +
     -   4*tred - 4*muu3**2*tred - 4*muu4**2*tred + 4*tred**2))))*
     -   VSfSf(3,isf1,isf2,itt1,itt2,igen1,igen2)*VSfSf(y,isf1,isf2,itt1,itt2,igen1,igen2)/(DenV(3)*DenV(y))
  
	enddo

!        !structure A2 for photons, Z:
    
	Do y=2,3
	
	colfac=Nc
	g1L=photongluon2l+photongluino2l
	g2L=BL(y,3,3,3,3)
	g1R=photongluon2r+photongluino2r
	g2R=BR(y,3,3,3,3)

      staustau2QQ_M2Vertex =  staustau2QQ_M2Vertex -2d0*sqrtS**5*(g1R*(-(g2R*muu4*
     - (-1 + 2*muu1**2 + muu3**2 + muu4**2 - 2*tred)) +
     - g2L*muu3*(-1 + 2*muu2**2 + muu3**2 + muu4**2 - 2*tred)) +
     - g1L*(-(g2L*muu4*(-1 + 2*muu1**2 + muu3**2 + muu4**2 -
     - 2*tred)) + g2R*muu3*
     - (-1 + 2*muu2**2 + muu3**2 + muu4**2 - 2*tred)))*
     - (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 2*tred)*
     - VSfSf(2,isf1,isf2,itt1,itt2,igen1,igen2)*VSfSf(y,isf1,isf2,itt1,itt2,igen1,igen2)/(2d0*DenV(2)*DenV(y))
     
	enddo
! ! 
	Do y=2,3
	
	colfac=Nc
	g1L=Zgluon2l+Zgluino2l
	g2L=BL(y,3,3,3,3)
	g1R=Zgluon2r+Zgluino2r
	g2R=BR(y,3,3,3,3)

      staustau2QQ_M2Vertex =  staustau2QQ_M2Vertex -2d0*sqrtS**5*(g1R*(-(g2R*muu4*
     - (-1 + 2*muu1**2 + muu3**2 + muu4**2 - 2*tred)) +
     - g2L*muu3*(-1 + 2*muu2**2 + muu3**2 + muu4**2 - 2*tred)) +
     - g1L*(-(g2L*muu4*(-1 + 2*muu1**2 + muu3**2 + muu4**2 -
     - 2*tred)) + g2R*muu3*
     - (-1 + 2*muu2**2 + muu3**2 + muu4**2 - 2*tred)))*
     - (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 2*tred)*
     -     VSfSf(3,isf1,isf2,itt1,itt2,igen1,igen2)*VSfSf(y,isf1,isf2,itt1,itt2,igen1,igen2)/(2d0*DenV(3)*DenV(y))
     
	enddo
! 
! ! 
!      !structure A3 for photons, Z
! !      
	Do y=2,3
	
	colfac=Nc
	g1L=photongluon3l+photongluino3l
	g2L=BL(y,3,3,3,3)
	g1R=photongluon3r+photongluino3r
	g2R=BR(y,3,3,3,3)

      staustau2QQ_M2Vertex =  staustau2QQ_M2Vertex - 2d0*(muu1**2 - muu2**2)*sqrtS**5*
     - (g1R*(-(g2R*muu4*(-1 + 2*muu1**2 + muu3**2 + muu4**2 -
     -  2*tred)) + g2L*muu3*
     - (-1 + 2*muu2**2 + muu3**2 + muu4**2 - 2*tred)) +
     - g1L*(-(g2L*muu4*(-1 + 2*muu1**2 + muu3**2 + muu4**2 -
     - 2*tred)) + g2R*muu3*
     - (-1 + 2*muu2**2 + muu3**2 + muu4**2 - 2*tred)))*
     -    VSfSf(2,isf1,isf2,itt1,itt2,igen1,igen2)*VSfSf(y,isf1,isf2,itt1,itt2,igen1,igen2)/(2d0*DenV(2)*DenV(y))
     
	enddo
!  
	Do y=2,3
	
	colfac=Nc
	g1L=Zgluon3l+Zgluino3l
	g2L=BL(y,3,3,3,3)
	g1R=Zgluon3r+Zgluino3r
	g2R=BR(y,3,3,3,3)

      staustau2QQ_M2Vertex =  staustau2QQ_M2Vertex -2d0*(muu1**2 - muu2**2)*sqrtS**5*
     - (g1R*(-(g2R*muu4*(-1 + 2*muu1**2 + muu3**2 + muu4**2 -
     -  2*tred)) + g2L*muu3*
     - (-1 + 2*muu2**2 + muu3**2 + muu4**2 - 2*tred)) +
     - g1L*(-(g2L*muu4*(-1 + 2*muu1**2 + muu3**2 + muu4**2 -
     - 2*tred)) + g2R*muu3*
     - (-1 + 2*muu2**2 + muu3**2 + muu4**2 - 2*tred)))*
     -    VSfSf(3,isf1,isf2,itt1,itt2,igen1,igen2)*VSfSf(y,isf1,isf2,itt1,itt2,igen1,igen2)/(2d0*DenV(3)*DenV(y))
     
	enddo

! ! 
! ! 
! ! !!---------------------------------Higgs corrected, interfered with vector contribution -------------------------------- 
! ! ! 
! ! ! 
      Do x=1,2
      Do y=2,3
	
      colfac=Nc
      g1l=BL(y,3,3,3,3)
      g1r=BR(y,3,3,3,3)

      staustau2QQ_M2Vertex =  staustau2QQ_M2Vertex +2d0*colfac*hSf(x)*(higgsgluinocorr(x)+higgsgluoncorr(x))*(g1l+g1r)*(sqrtS**3*(muu3**3 + muu3**2*muu4 +
     -    muu4*(-1 + 2*muu1**2 + muu4**2 - 2*tred) +
     -    muu3*(-1 + 2*muu2**2 + muu4**2 - 2*tred)))*(VSfSf(y,isf1,isf2,itt1,itt2,igen1,igen2)/(DenV(y)*DenSH(x)))
	
      enddo
      enddo  
! ! 
! ! ! !!---------------------------------Photon corrected, interfered with Higgs contribution-------------------------------- 
! ! ! 
! ! ! 
! ! 
!      !Tree level proportional
! 
      Do x=1,2
! 	
      colfac=Nc
      g1l=Photongluon1l+photongluino1l
      g1r=Photongluon1r+photongluino1r

      staustau2QQ_M2Vertex =  staustau2QQ_M2Vertex +2d0*colfac*hSf(x)*hff(x)*(g1l+g1r)*(sqrtS**3*(muu3**3 + muu3**2*muu4 +
     -    muu4*(-1 + 2*muu1**2 + muu4**2 - 2*tred) +
     -    muu3*(-1 + 2*muu2**2 + muu4**2 - 2*tred)))*VSfSf(2,isf1,isf2,itt1,itt2,igen1,igen2)/(DenV(2)*DenSH(x))

      enddo  

      
      !structure A2
    
      Do x=1,2
	
      colfac=Nc
      g1l=Photongluon2l+photongluino2l
      g1r=Photongluon2r+photongluino2r   
     
      staustau2QQ_M2Vertex =  staustau2QQ_M2Vertex+2d0*colfac*sqrtS**4*hSf(x)*hff(x)*(g1L+g1R)*(muu3**2+2*muu3*muu4+muu4**2-1)*(muu1**2+muu2**2+muu3**2+muu4**2-2*tred-1)*
     -   (VSfSf(2,isf1,isf2,itt1,itt2,igen1,igen2)/(DenV(2)*DenSH(x))) 

      enddo
      
      !structure A3  
     
      Do x=1,2

	
      colfac=Nc
      g1l=Photongluon3l+photongluino3l
      g1r=Photongluon3r+photongluino3r
  
      staustau2QQ_M2Vertex =  staustau2QQ_M2Vertex+2d0*colfac*sqrtS**4*hSf(x)*hff(x)*(-(g1L+g1R))*(muu1**2-muu2**2)*(muu3**2+2*muu3*muu4+muu4**2-1)*
     -   (VSfSf(2,isf1,isf2,itt1,itt2,igen1,igen2)/(DenV(2)*DenSH(x)))

      enddo
!       
! ! ! !!---------------------------------Z corrected, interfered with Higgs contribution--------------------------------       
! 
!      !Tree level proportional

      Do x=1,2
	
      colfac=Nc
      g1l=Zgluon1l+Zgluino1l
      g1r=Zgluon1r+Zgluino1r

      staustau2QQ_M2Vertex =  staustau2QQ_M2Vertex +2d0*colfac*hSf(x)*hff(x)*(g1l+g1r)*(sqrtS**3*(muu3**3 + muu3**2*muu4 +
     -    muu4*(-1 + 2*muu1**2 + muu4**2 - 2*tred) +
     -    muu3*(-1 + 2*muu2**2 + muu4**2 - 2*tred)))*(VSfSf(3,isf1,isf2,itt1,itt2,igen1,igen2)/(DenV(3)*DenSH(x)))

      enddo  
!       
! !     
! !      !structure 2
! !      
      Do x=1,2
	
      colfac=Nc
      g1l=Zgluon2l+Zgluino2l
      g1r=Zgluon2r+Zgluino2r   
     
      staustau2QQ_M2Vertex =  staustau2QQ_M2Vertex+2d0*colfac*sqrtS**4*hSf(x)*hff(x)*(g1L+g1R)*(muu3**2+2*muu3*muu4+muu4**2-1)*(muu1**2+muu2**2+muu3**2+muu4**2-2*tred-1)*
     -   (VSfSf(3,isf1,isf2,itt1,itt2,igen1,igen2)/(DenV(3)*DenSH(x))) 

      enddo
!       
! !       
! !      !structure 3  
     
      Do x=1,2

      colfac=Nc
      g1l=Zgluon3l+Zgluino3l
      g1r=Zgluon3r+Zgluino3r
  
      staustau2QQ_M2Vertex =  staustau2QQ_M2Vertex+2d0*colfac*sqrtS**4*hSf(x)*hff(x)*(-(g1L+g1R))*(muu1**2-muu2**2)*(muu3**2+2*muu3*muu4+muu4**2-1)*
     -   (VSfSf(3,isf1,isf2,itt1,itt2,igen1,igen2)/(DenV(3)*DenSH(x)))

      enddo
      

!        
	end
