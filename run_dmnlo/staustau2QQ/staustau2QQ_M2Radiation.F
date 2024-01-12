
      double precision function staustau2QQ_M2Radiation()
      implicit none
#include "staustau2QQ_Kinematics.h"
#include "staustau2QQ_Model.h"
#include "staustau2QQ_GenCouplings.h"

      integer x,y,h,hc
      double precision colfac,g0R,g0L,g1L,g2L,g1R,g2R,g0,dipolecont
      double complex test
      complex*16 DenV(3),DenSH(2)

!Propagator
      DenV(1)=0d0
      DenV(2) = sqrtS**2 
      DenV(3) = (sqrtS**2-MZ**2)

      do x=1,2
            DenSH(x) = (sqrtS**2-MHk(x)**2)	
      enddo

      staustau2QQ_M2Radiation = 0d0
          
!       !!Dipole contribution following 0907.0030 (D312,D321)

      dipolecont= 2*gs**2/(sqrtS**2)*CF*(((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) -   !!factorising with Born cross section
     -  (Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*
     -  (muu3**4 - (1 + 3*muu4**2 - x1 - 2*x2)*(1 + muu4**2 - x2) + muu3**2*(-2*muu4**2 - x1 + x2)))/
     -  ((-1 + muu3**2 + muu4**2)*(-1 + muu3**2 - muu4**2 + x2)*Sqrt(-4*muu4**2 + x2**2)))+((2*(-1 + muu3**2 + muu4**2))/(-2 + x1 + x2) - 
     -  (Sqrt(muu3**4 + (-1 + muu4**2)**2 - 2*muu3**2*(1 + muu4**2))*
     -  (3*muu3**4 - muu4**4 + muu3**2*(4 + 2*muu4**2 - 5*x1 - x2) + muu4**2*(-x1 + x2) + (-1 + x1)*(-1 + 2*x1 + x2)))/
     -   ((-1 + muu3**2 + muu4**2)*(1 + muu3**2 - muu4**2 - x1)*Sqrt(-4*muu3**2 + x1**2)))) 

      
      !!Interference terms
      
      
!       Do x = 1,2
!        Do y = 2,3
!       
!       g1l=BL(y,3,3,3,3)
!       g1r=BR(y,3,3,3,3)
!       
!       staustau2QQ_M2Radiation = staustau2QQ_M2Radiation - CF*Nc*gs**2*(VSfSf(y,isf1,isf2,itt1,itt2,igen1,igen2)*hSf(x)*hff(x)/(DenV(y)*DenSH(x))*4*
!      -(2*g1L*k2k3*muu3*p1k1*sqrtS + 2*g1R*k2k3*muu3*p1k1*sqrtS +
!      - 4*g1L*k1k2*muu4*p1k1*sqrtS + 4*g1R*k1k2*muu4*p1k1*sqrtS +
!      - 2*g1L*k1k3*muu4*p1k1*sqrtS + 2*g1R*k1k3*muu4*p1k1*sqrtS -
!      - 4*g1L*k1k2*muu3*p1k2*sqrtS - 4*g1R*k1k2*muu3*p1k2*sqrtS -
!      - 2*g1L*k2k3*muu3*p1k2*sqrtS - 2*g1R*k2k3*muu3*p1k2*sqrtS -
!      - 2*g1L*k1k3*muu4*p1k2*sqrtS - 2*g1R*k1k3*muu4*p1k2*sqrtS -
!      - 2*g1L*k1k2*muu3*p1k3*sqrtS - 2*g1R*k1k2*muu3*p1k3*sqrtS -
!      - 4*g1L*k2k3*muu3*p1k3*sqrtS - 4*g1R*k2k3*muu3*p1k3*sqrtS +
!      - 2*g1L*k1k2*muu4*p1k3*sqrtS + 2*g1R*k1k2*muu4*p1k3*sqrtS +
!      - 4*g1L*k1k3*muu4*p1k3*sqrtS + 4*g1R*k1k3*muu4*p1k3*sqrtS -
!      - 2*g1L*k2k3*muu3*p2k1*sqrtS - 2*g1R*k2k3*muu3*p2k1*sqrtS -
!      - 4*g1L*k1k2*muu4*p2k1*sqrtS - 4*g1R*k1k2*muu4*p2k1*sqrtS -
!      - 2*g1L*k1k3*muu4*p2k1*sqrtS - 2*g1R*k1k3*muu4*p2k1*sqrtS +
!      - 4*g1L*k1k2*muu3*p2k2*sqrtS + 4*g1R*k1k2*muu3*p2k2*sqrtS +
!      - 2*g1L*k2k3*muu3*p2k2*sqrtS + 2*g1R*k2k3*muu3*p2k2*sqrtS +
!      - 2*g1L*k1k3*muu4*p2k2*sqrtS + 2*g1R*k1k3*muu4*p2k2*sqrtS +
!      - 2*g1L*k1k2*muu3*p2k3*sqrtS + 2*g1R*k1k2*muu3*p2k3*sqrtS +
!      - 4*g1L*k2k3*muu3*p2k3*sqrtS + 4*g1R*k2k3*muu3*p2k3*sqrtS -
!      - 2*g1L*k1k2*muu4*p2k3*sqrtS - 2*g1R*k1k2*muu4*p2k3*sqrtS -
!      - 4*g1L*k1k3*muu4*p2k3*sqrtS - 4*g1R*k1k3*muu4*p2k3*sqrtS -
!      - 2*g1L*muu3**2*muu4*p1k1*sqrtS**3 -
!      - 2*g1R*muu3**2*muu4*p1k1*sqrtS**3 -
!      - 2*g1L*muu4**3*p1k1*sqrtS**3 - 2*g1R*muu4**3*p1k1*sqrtS**3 +
!      - 2*g1L*muu3**3*p1k2*sqrtS**3 + 2*g1R*muu3**3*p1k2*sqrtS**3 +
!      - 2*g1L*muu3*muu4**2*p1k2*sqrtS**3 +
!      - 2*g1R*muu3*muu4**2*p1k2*sqrtS**3 -
!      - 2*g1L*muu3**2*muu4*p1k3*sqrtS**3 -
!      - 2*g1R*muu3**2*muu4*p1k3*sqrtS**3 +
!      - 2*g1L*muu3*muu4**2*p1k3*sqrtS**3 +
!      - 2*g1R*muu3*muu4**2*p1k3*sqrtS**3 +
!      - 2*g1L*muu3**2*muu4*p2k1*sqrtS**3 +
!      - 2*g1R*muu3**2*muu4*p2k1*sqrtS**3 +
!      - 2*g1L*muu4**3*p2k1*sqrtS**3 + 2*g1R*muu4**3*p2k1*sqrtS**3 -
!      - 2*g1L*muu3**3*p2k2*sqrtS**3 - 2*g1R*muu3**3*p2k2*sqrtS**3 -
!      - 2*g1L*muu3*muu4**2*p2k2*sqrtS**3 -
!      - 2*g1R*muu3*muu4**2*p2k2*sqrtS**3 +
!      - 2*g1L*muu3**2*muu4*p2k3*sqrtS**3 +
!      - 2*g1R*muu3**2*muu4*p2k3*sqrtS**3 -
!      - 2*g1L*muu3*muu4**2*p2k3*sqrtS**3 -
!      - 2*g1R*muu3*muu4**2*p2k3*sqrtS**3))/(4d0*k1k3**2)
!       
!        enddo
!       enddo 
! 
!   !! DIPOLE 
!   
!   
      Do x = 1,2
      Do y = 2,3
      
      g1l=BL(y,3,3,3,3)
      g1r=BR(y,3,3,3,3)     

      staustau2QQ_M2Radiation = staustau2QQ_M2Radiation - Nc*VSfSf(y,isf1,isf2,itt1,itt2,igen1,igen2)*hSf(x)*hff(x)/(DenV(y)*DenSH(x))*
     -  (4*(-(g1L*muu4*p1k1*sqrtS)/2. - (g1R*muu4*p1k1*sqrtS)/2. +
     -  (g1L*muu3*p1k2*sqrtS)/2. + (g1R*muu3*p1k2*sqrtS)/2. +
     -  (g1L*muu4*p2k1*sqrtS)/2. + (g1R*muu4*p2k1*sqrtS)/2. -
     -  (g1L*muu3*p2k2*sqrtS)/2. - (g1R*muu3*p2k2*sqrtS)/2.))*dipolecont/(4d0*k1k3**2)

      enddo
      enddo   
      
      
      !Same for Higgs
      
      Do x = 1,2
      Do y = 1,2
      

      staustau2QQ_M2Radiation = staustau2QQ_M2Radiation - CF*gs**2*Nc*hSf(x)*hff(x)*hsf(y)*hff(y)/(DenSH(y)*DenSH(x))*
     -   (4*(-8*k1k2**2 - 8*k1k2*k1k3 - 8*k1k2*k2k3 - 16*k1k3*k2k3 +
     -   4*k1k2*muu3**2*sqrtS**2 + 8*k2k3*muu3**2*sqrtS**2 +
     -   8*k1k2*muu3*muu4*sqrtS**2 + 4*k1k2*muu4**2*sqrtS**2 +
     -   8*k1k3*muu4**2*sqrtS**2 - 4*muu3**3*muu4*sqrtS**4 -
     -   4*muu3*muu4**3*sqrtS**4))/(4d0*k1k3**2) 
    
      enddo
      enddo  

      Do x=2,2
      Do y=2,2
	
	colfac=Nc
	g1L=BL(x,3,3,3,3)
	g2L=BL(y,3,3,3,3)
	g1R=BR(x,3,3,3,3)
	g2R=BR(y,3,3,3,3)
	
      staustau2QQ_M2Radiation = staustau2QQ_M2Radiation -CF*gs**2*VSfSf(x,isf1,isf2,itt1,itt2,igen1,igen2)*VSfSf(y,isf1,isf2,itt1,itt2,igen1,igen2)/(DenV(x)*DenV(y))*
     - (4*(4*g1L*g2L*k2k3*p1k1**2 + 4*g1R*g2R*k2k3*p1k1**2 -
     - 8*g1L*g2L*k1k2*p1k1*p1k2 - 8*g1R*g2R*k1k2*p1k1*p1k2 -
     - 4*g1L*g2L*k1k3*p1k1*p1k2 - 4*g1R*g2R*k1k3*p1k1*p1k2 -
     - 4*g1L*g2L*k2k3*p1k1*p1k2 - 4*g1R*g2R*k2k3*p1k1*p1k2 +
     - 4*g1L*g2L*k1k3*p1k2**2 + 4*g1R*g2R*k1k3*p1k2**2 -
     - 4*g1L*g2L*k1k2*p1k1*p1k3 - 4*g1R*g2R*k1k2*p1k1*p1k3 -
     - 4*g1L*g2L*k2k3*p1k1*p1k3 - 4*g1R*g2R*k2k3*p1k1*p1k3 -
     - 4*g1L*g2L*k1k2*p1k2*p1k3 - 4*g1R*g2R*k1k2*p1k2*p1k3 -
     - 4*g1L*g2L*k1k3*p1k2*p1k3 - 4*g1R*g2R*k1k3*p1k2*p1k3 -
     - 8*g1L*g2L*k1k2**2*p1p2 - 8*g1R*g2R*k1k2**2*p1p2 -
     - 8*g1L*g2L*k1k2*k1k3*p1p2 - 8*g1R*g2R*k1k2*k1k3*p1p2 -
     - 8*g1L*g2L*k1k2*k2k3*p1p2 - 8*g1R*g2R*k1k2*k2k3*p1p2 -
     - 8*g1L*g2L*k1k3*k2k3*p1p2 - 8*g1R*g2R*k1k3*k2k3*p1p2 -
     - 8*g1L*g2L*k2k3*p1k1*p2k1 - 8*g1R*g2R*k2k3*p1k1*p2k1 +
     -  8*g1L*g2L*k1k2*p1k2*p2k1 + 8*g1R*g2R*k1k2*p1k2*p2k1 +
     - 4*g1L*g2L*k1k3*p1k2*p2k1 + 4*g1R*g2R*k1k3*p1k2*p2k1 +
     - 4*g1L*g2L*k2k3*p1k2*p2k1 + 4*g1R*g2R*k2k3*p1k2*p2k1 +
     - 4*g1L*g2L*k1k2*p1k3*p2k1 + 4*g1R*g2R*k1k2*p1k3*p2k1 +
     - 4*g1L*g2L*k2k3*p1k3*p2k1 + 4*g1R*g2R*k2k3*p1k3*p2k1 +
     - 4*g1L*g2L*k2k3*p2k1**2 + 4*g1R*g2R*k2k3*p2k1**2 +
     - 8*g1L*g2L*k1k2*p1k1*p2k2 + 8*g1R*g2R*k1k2*p1k1*p2k2 +
     -  4*g1L*g2L*k1k3*p1k1*p2k2 + 4*g1R*g2R*k1k3*p1k1*p2k2 +
     - 4*g1L*g2L*k2k3*p1k1*p2k2 + 4*g1R*g2R*k2k3*p1k1*p2k2 -
     - 8*g1L*g2L*k1k3*p1k2*p2k2 - 8*g1R*g2R*k1k3*p1k2*p2k2 +
     - 4*g1L*g2L*k1k2*p1k3*p2k2 + 4*g1R*g2R*k1k2*p1k3*p2k2 +
     - 4*g1L*g2L*k1k3*p1k3*p2k2 + 4*g1R*g2R*k1k3*p1k3*p2k2 -
     - 8*g1L*g2L*k1k2*p2k1*p2k2 - 8*g1R*g2R*k1k2*p2k1*p2k2 -
     -  4*g1L*g2L*k1k3*p2k1*p2k2 - 4*g1R*g2R*k1k3*p2k1*p2k2 -
     - 4*g1L*g2L*k2k3*p2k1*p2k2 - 4*g1R*g2R*k2k3*p2k1*p2k2 +
     - 4*g1L*g2L*k1k3*p2k2**2 + 4*g1R*g2R*k1k3*p2k2**2 +
     - 4*g1L*g2L*k1k2*p1k1*p2k3 + 4*g1R*g2R*k1k2*p1k1*p2k3 +
     - 4*g1L*g2L*k2k3*p1k1*p2k3 + 4*g1R*g2R*k2k3*p1k1*p2k3 +
     - 4*g1L*g2L*k1k2*p1k2*p2k3 + 4*g1R*g2R*k1k2*p1k2*p2k3 +
     - 4*g1L*g2L*k1k3*p1k2*p2k3 + 4*g1R*g2R*k1k3*p1k2*p2k3 -
     - 4*g1L*g2L*k1k2*p2k1*p2k3 - 4*g1R*g2R*k1k2*p2k1*p2k3 -
     - 4*g1L*g2L*k2k3*p2k1*p2k3 - 4*g1R*g2R*k2k3*p2k1*p2k3 -
     - 4*g1L*g2L*k1k2*p2k2*p2k3 - 4*g1R*g2R*k1k2*p2k2*p2k3 -
     - 4*g1L*g2L*k1k3*p2k2*p2k3 - 4*g1R*g2R*k1k3*p2k2*p2k3 +
     - 4*g1L*g2L*k1k2**2*muu1**2*sqrtS**2 +
     - 4*g1R*g2R*k1k2**2*muu1**2*sqrtS**2 +
     - 4*g1L*g2L*k1k2*k1k3*muu1**2*sqrtS**2 +
     - 4*g1R*g2R*k1k2*k1k3*muu1**2*sqrtS**2 +
     - 4*g1L*g2L*k1k2*k2k3*muu1**2*sqrtS**2 +
     - 4*g1R*g2R*k1k2*k2k3*muu1**2*sqrtS**2 +
     - 4*g1L*g2L*k1k3*k2k3*muu1**2*sqrtS**2 +
     - 4*g1R*g2R*k1k3*k2k3*muu1**2*sqrtS**2 +
     - 4*g1L*g2L*k1k2**2*muu2**2*sqrtS**2 +
     - 4*g1R*g2R*k1k2**2*muu2**2*sqrtS**2 +
     - 4*g1L*g2L*k1k2*k1k3*muu2**2*sqrtS**2 +
     - 4*g1R*g2R*k1k2*k1k3*muu2**2*sqrtS**2 +
     - 4*g1L*g2L*k1k2*k2k3*muu2**2*sqrtS**2 +
     - 4*g1R*g2R*k1k2*k2k3*muu2**2*sqrtS**2 +
     - 4*g1L*g2L*k1k3*k2k3*muu2**2*sqrtS**2 +
     - 4*g1R*g2R*k1k3*k2k3*muu2**2*sqrtS**2 +
     - 4*g1L*g2L*muu3**2*p1k1*p1k2*sqrtS**2 +
     - 4*g1R*g2R*muu3**2*p1k1*p1k2*sqrtS**2 +
     - 4*g1L*g2L*muu4**2*p1k1*p1k2*sqrtS**2 +
     - 4*g1R*g2R*muu4**2*p1k1*p1k2*sqrtS**2 +
     -  4*g1L*g2L*muu4**2*p1k1*p1k3*sqrtS**2 +
     - 4*g1R*g2R*muu4**2*p1k1*p1k3*sqrtS**2 +
     -  4*g1L*g2L*muu3**2*p1k2*p1k3*sqrtS**2 +
     - 4*g1R*g2R*muu3**2*p1k2*p1k3*sqrtS**2 +
     - 4*g1R*g2L*muu3*muu4*p1k3**2*sqrtS**2 +
     - 4*g1L*g2R*muu3*muu4*p1k3**2*sqrtS**2 +
     - 4*g1L*g2L*k1k2*muu3**2*p1p2*sqrtS**2 +
     - 4*g1R*g2R*k1k2*muu3**2*p1p2*sqrtS**2 +
     - 8*g1L*g2L*k2k3*muu3**2*p1p2*sqrtS**2 +
     - 8*g1R*g2R*k2k3*muu3**2*p1p2*sqrtS**2 -
     - 8*g1R*g2L*k1k2*muu3*muu4*p1p2*sqrtS**2 -
     - 8*g1L*g2R*k1k2*muu3*muu4*p1p2*sqrtS**2 +
     - 4*g1L*g2L*k1k2*muu4**2*p1p2*sqrtS**2 +
     - 4*g1R*g2R*k1k2*muu4**2*p1p2*sqrtS**2 +
     - 8*g1L*g2L*k1k3*muu4**2*p1p2*sqrtS**2 +
     - 8*g1R*g2R*k1k3*muu4**2*p1p2*sqrtS**2 -
     - 4*g1L*g2L*muu3**2*p1k2*p2k1*sqrtS**2 -
     - 4*g1R*g2R*muu3**2*p1k2*p2k1*sqrtS**2 -
     - 4*g1L*g2L*muu4**2*p1k2*p2k1*sqrtS**2 -
     - 4*g1R*g2R*muu4**2*p1k2*p2k1*sqrtS**2 -
     - 4*g1L*g2L*muu4**2*p1k3*p2k1*sqrtS**2 -
     - 4*g1R*g2R*muu4**2*p1k3*p2k1*sqrtS**2 -
     - 4*g1L*g2L*muu3**2*p1k1*p2k2*sqrtS**2 -
     - 4*g1R*g2R*muu3**2*p1k1*p2k2*sqrtS**2 -
     - 4*g1L*g2L*muu4**2*p1k1*p2k2*sqrtS**2 -
     - 4*g1R*g2R*muu4**2*p1k1*p2k2*sqrtS**2 -
     - 4*g1L*g2L*muu3**2*p1k3*p2k2*sqrtS**2 -
     - 4*g1R*g2R*muu3**2*p1k3*p2k2*sqrtS**2 +
     - 4*g1L*g2L*muu3**2*p2k1*p2k2*sqrtS**2 +
     - 4*g1R*g2R*muu3**2*p2k1*p2k2*sqrtS**2 +
     - 4*g1L*g2L*muu4**2*p2k1*p2k2*sqrtS**2 +
     - 4*g1R*g2R*muu4**2*p2k1*p2k2*sqrtS**2 -
     - 4*g1L*g2L*muu4**2*p1k1*p2k3*sqrtS**2 -
     - 4*g1R*g2R*muu4**2*p1k1*p2k3*sqrtS**2 -
     - 4*g1L*g2L*muu3**2*p1k2*p2k3*sqrtS**2 -
     - 4*g1R*g2R*muu3**2*p1k2*p2k3*sqrtS**2 -
     - 8*g1R*g2L*muu3*muu4*p1k3*p2k3*sqrtS**2 -
     - 8*g1L*g2R*muu3*muu4*p1k3*p2k3*sqrtS**2 +
     - 4*g1L*g2L*muu4**2*p2k1*p2k3*sqrtS**2 +
     - 4*g1R*g2R*muu4**2*p2k1*p2k3*sqrtS**2 +
     - 4*g1L*g2L*muu3**2*p2k2*p2k3*sqrtS**2 +
     - 4*g1R*g2R*muu3**2*p2k2*p2k3*sqrtS**2 +
     - 4*g1R*g2L*muu3*muu4*p2k3**2*sqrtS**2 +
     - 4*g1L*g2R*muu3*muu4*p2k3**2*sqrtS**2 -
     - 2*g1L*g2L*k1k2*muu1**2*muu3**2*sqrtS**4 -
     - 2*g1R*g2R*k1k2*muu1**2*muu3**2*sqrtS**4 -
     -    4*g1L*g2L*k2k3*muu1**2*muu3**2*sqrtS**4 -
     - 4*g1R*g2R*k2k3*muu1**2*muu3**2*sqrtS**4 -
     - 2*g1L*g2L*k1k2*muu2**2*muu3**2*sqrtS**4 -
     - 2*g1R*g2R*k1k2*muu2**2*muu3**2*sqrtS**4 -
     - 4*g1L*g2L*k2k3*muu2**2*muu3**2*sqrtS**4 -
     - 4*g1R*g2R*k2k3*muu2**2*muu3**2*sqrtS**4 +
     - 4*g1R*g2L*k1k2*muu1**2*muu3*muu4*sqrtS**4 +
     - 4*g1L*g2R*k1k2*muu1**2*muu3*muu4*sqrtS**4 +
     - 4*g1R*g2L*k1k2*muu2**2*muu3*muu4*sqrtS**4 +
     -  4*g1L*g2R*k1k2*muu2**2*muu3*muu4*sqrtS**4 -
     - 2*g1L*g2L*k1k2*muu1**2*muu4**2*sqrtS**4 -
     - 2*g1R*g2R*k1k2*muu1**2*muu4**2*sqrtS**4 -
     - 4*g1L*g2L*k1k3*muu1**2*muu4**2*sqrtS**4 -
     - 4*g1R*g2R*k1k3*muu1**2*muu4**2*sqrtS**4 -
     - 2*g1L*g2L*k1k2*muu2**2*muu4**2*sqrtS**4 -
     - 2*g1R*g2R*k1k2*muu2**2*muu4**2*sqrtS**4 -
     - 4*g1L*g2L*k1k3*muu2**2*muu4**2*sqrtS**4 -
     - 4*g1R*g2R*k1k3*muu2**2*muu4**2*sqrtS**4 +
     - 4*g1R*g2L*muu3**3*muu4*p1p2*sqrtS**4 +
     - 4*g1L*g2R*muu3**3*muu4*p1p2*sqrtS**4 +
     - 4*g1R*g2L*muu3*muu4**3*p1p2*sqrtS**4 +
     - 4*g1L*g2R*muu3*muu4**3*p1p2*sqrtS**4 -
     - 2*g1R*g2L*muu1**2*muu3**3*muu4*sqrtS**6 -
     - 2*g1L*g2R*muu1**2*muu3**3*muu4*sqrtS**6 -
     - 2*g1R*g2L*muu2**2*muu3**3*muu4*sqrtS**6 -
     - 2*g1L*g2R*muu2**2*muu3**3*muu4*sqrtS**6 -
     - 2*g1R*g2L*muu1**2*muu3*muu4**3*sqrtS**6 -
     - 2*g1L*g2R*muu1**2*muu3*muu4**3*sqrtS**6 -
     - 2*g1R*g2L*muu2**2*muu3*muu4**3*sqrtS**6 -
     - 2*g1L*g2R*muu2**2*muu3*muu4**3*sqrtS**6))/(4d0*k1k3**2) 
	
      enddo
      enddo 
      
      Do x=2,2
      Do y=2,2
	
	colfac=Nc
	g1L=BL(x,3,3,3,3)
	g2L=BL(y,3,3,3,3)
	g1R=BR(x,3,3,3,3)
	g2R=BR(y,3,3,3,3)
	
      staustau2QQ_M2Radiation = staustau2QQ_M2Radiation -VSfSf(x,isf1,isf2,itt1,itt2,igen1,igen2)*VSfSf(y,isf1,isf2,itt1,itt2,igen1,igen2)/(DenV(x)*DenV(y))*
     -   (4*(-2*p1p2 + muu1**2*sqrtS**2 + muu2**2*sqrtS**2)*
     - (g1L*g2L*p1k1*p1k2 + g1R*g2R*p1k1*p1k2 + g1L*g2L*k1k2*p1p2 +
     - g1R*g2R*k1k2*p1p2 - g1L*g2L*p1k2*p2k1 - g1R*g2R*p1k2*p2k1 -
     - g1L*g2L*p1k1*p2k2 - g1R*g2R*p1k1*p2k2 + g1L*g2L*p2k1*p2k2 +
     -  g1R*g2R*p2k1*p2k2 - (g1L*g2L*k1k2*muu1**2*sqrtS**2)/2. -
     - (g1R*g2R*k1k2*muu1**2*sqrtS**2)/2. -
     - (g1L*g2L*k1k2*muu2**2*sqrtS**2)/2. -
     - (g1R*g2R*k1k2*muu2**2*sqrtS**2)/2. +
     - g1R*g2L*muu3*muu4*p1p2*sqrtS**2 +
     - g1L*g2R*muu3*muu4*p1p2*sqrtS**2 -
     - (g1R*g2L*muu1**2*muu3*muu4*sqrtS**4)/2. -
     -  (g1L*g2R*muu1**2*muu3*muu4*sqrtS**4)/2. -
     - (g1R*g2L*muu2**2*muu3*muu4*sqrtS**4)/2. -
     - (g1L*g2R*muu2**2*muu3*muu4*sqrtS**4)/2.))*dipolecont/(4d0*k1k3**2) 
	
      enddo
      enddo 
  
      
      
      
      end
