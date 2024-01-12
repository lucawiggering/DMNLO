
	function staustau2QQ_M2Counter()
	implicit none
	
#include "staustau2QQ_Kinematics.h"
#include "staustau2QQ_Model.h"
#include "staustau2QQ_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"
#include "../util/DMNLO_RenScheme.h"

	complex*16 staustau2QQ_M2Counter
	complex*16 Den(4)
	
      integer i,j,k,l,m,n,x,y,h
	double precision colfac,g0R,g0L,g1L,g2L,g1R,g2R
	double complex test
	complex*16 DenV(3),DenSH(2)

! ! Propagator
      DenV(1)=0d0
	DenV(2) = sqrtS**2 
	DenV(3) = (sqrtS**2-MZ**2)
	
	do x=1,2
	  DenSH(x) = (sqrtS**2-MHk(x)**2)	
	enddo
	
	! ! --------------------- Vector boson squared -----------
	
      staustau2QQ_M2Counter = 0d0

	Do x=2,3
	Do y=2,3
	
	colfac=Nc
	g1L=(dMf(3,3)/(Mf(3,3))+dZfL(3,3)/(2d0)+dZfR(3,3)/(2d0))*BL(x,3,3,3,3)
	g2L=BL(y,3,3,3,3)
	g1R=(dMf(3,3)/(Mf(3,3))+dZfL(3,3)/(2d0)+dZfR(3,3)/(2d0))*BR(x,3,3,3,3)
	g2R=BR(y,3,3,3,3)
 
      staustau2QQ_M2Counter = staustau2QQ_M2Counter +2d0*colfac*(sqrtS**4*(g1L*(-2*g2R*(-1 + 2*muu1**2 + 2*muu2**2)*muu3*muu4 +
     - g2L*(2 - 3*muu3**2 + muu3**4 - 3*muu4**2 +
     - 2*muu3**2*muu4**2 + muu4**4 +
     - 4*muu2**2*(-1 + muu3**2 + muu4**2 - tred) + 
     - 4*muu1**2*(-1 + muu2**2 + muu3**2 + muu4**2 - tred) +
     - 4*tred - 4*muu3**2*tred - 4*muu4**2*tred + 4*tred**2)) +
     - g1R*(-2*g2L*(-1 + 2*muu1**2 + 2*muu2**2)*muu3*muu4 +
     - g2R*(2 - 3*muu3**2 + muu3**4 - 3*muu4**2 +
     - 2*muu3**2*muu4**2 + muu4**4 +
     - 4*muu2**2*(-1 + muu3**2 + muu4**2 - tred) +
     - 4*muu1**2*(-1 + muu2**2 + muu3**2 + muu4**2 - tred) +
     - 4*tred - 4*muu3**2*tred - 4*muu4**2*tred + 4*tred**2))))*
     - VSfSf(x,isf1,isf2,itt1,itt2,igen1,igen2)*VSfSf(y,isf1,isf2,itt1,itt2,igen1,igen2)/(2d0*DenV(x)*DenV(y))

	enddo
	enddo
! ! 
! ! 	
! 	! ! --------------------- Higgs boson squared -----------	
! 
	Do x=1,2
	Do y=1,2
	
      colfac=0d0
      colfac=colfac+Nc
     
      staustau2QQ_M2Counter =  staustau2QQ_M2Counter + 2d0*colfac*hSf(x)*hsf(y)*
     -    (dMf(3,3)/(Mf(3,3))+dZfL(3,3)/(2d0)+dZfR(3,3)/(2d0))*hff(x)*hff(y)/(DenSH(x)*DenSH(y))*
     -    (-2*(-1 + muu3**2 + 2*muu3*muu4 + muu4**2)*sqrtS**2)
	
	enddo
	enddo


    ! ! --------------------- Higgs vector interference, Higgs corrected -----------	
!     
	Do x=1,2
	Do y=2,3
	
	colfac=Nc
	g1l=BL(x,3,3,3,3)
	g1r=BR(x,3,3,3,3)

      staustau2QQ_M2Counter = staustau2QQ_M2Counter +2d0*colfac*hSf(x)*hff(x)* (dMf(3,3)/(Mf(3,3))+dZfL(3,3)/(2d0)+dZfR(3,3)/(2d0))*(g1l+g1r)*(sqrtS**3*(muu3**3 + muu3**2*muu4 +
     -    muu4*(-1 + 2*muu1**2 + muu4**2 - 2*tred) +
     -    muu3*(-1 + 2*muu2**2 + muu4**2 - 2*tred)))*(VSfSf(y,isf1,isf2,itt1,itt2,igen1,igen2)/(DenV(y)*DenSH(x)))
	
	enddo
	enddo  
! 
! !     ! ! --------------------- Higgs vector interference, vector corrected -----------	
!     
	Do x=1,2
	Do y=2,3
	
	colfac=Nc
	g1l=(dMf(3,3)/(Mf(3,3))+dZfL(3,3)/(2d0)+dZfR(3,3)/(2d0))*BL(x,3,3,3,3)
	g1r=(dMf(3,3)/(Mf(3,3))+dZfL(3,3)/(2d0)+dZfR(3,3)/(2d0))*BR(x,3,3,3,3)

      staustau2QQ_M2Counter = staustau2QQ_M2Counter +2d0*colfac*hSf(x)*hff(x)*(g1l+g1r)*(sqrtS**3*(muu3**3 + muu3**2*muu4 +
     -    muu4*(-1 + 2*muu1**2 + muu4**2 - 2*tred) +
     -    muu3*(-1 + 2*muu2**2 + muu4**2 - 2*tred)))*(VSfSf(y,isf1,isf2,itt1,itt2,igen1,igen2)/(DenV(y)*DenSH(x)))
	
	enddo
	enddo  

      end
