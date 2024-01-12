      double precision function stsT2QQbar_M2VertexCounter2sqG()
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"
#include "../common/stsT2xx_GenCounterterms.h"

      double precision colorf,Nf
      double complex DenSV


!Propagator
      DenSV = sman

! Number of effectively massless quarks
      Nf = 4d0

***************************************************
* Counterterms for Sfermion-Sfermion-Gluon vertex *
***************************************************
	
c-------------------------------------s-Kanal Gluon---------------------------------------------------------
      colorf = 1/4d0*(Nc**2-1)

      stsT2QQbar_M2VertexCounter2sqG =  -Nf*2*colorf*dble((sqrtS**4*((-GS)*
     &        ((muu3**4 - muu4**2 + muu4**4 + muu3**2*(-1 + 2*muu4**2 - 4*tred) + 
     &              4*muu1**2*(muu2**2 - tred) + 4*tred - 4*muu2**2*tred - 
     &              4*muu4**2*tred + 4*tred**2)*(-GS) +  2*(-1 + 2*muu1**2 + 2*muu2**2)*muu3*muu4*(-GS)
     &           ) + (-GS)*(2*(-1 + 2*muu1**2 + 2*muu2**2)*muu3*muu4*(-GS) + 
     &           (muu3**4 - muu4**2 + muu4**4 + muu3**2*(-1 + 2*muu4**2 - 4*tred) + 
     &              4*muu1**2*(muu2**2 - tred) + 4*tred - 4*muu2**2*tred - 
     &              4*muu4**2*tred + 4*tred**2)*(-GS)))*
     &      dconjg(SfSfg(isf1,itt1,igen1,isf2,itt2,igen2))*
     &      dSfSfg(isf1,itt1,igen1,isf2,itt2,igen2))/(dconjg(DenSV)*DenSV))

      end

