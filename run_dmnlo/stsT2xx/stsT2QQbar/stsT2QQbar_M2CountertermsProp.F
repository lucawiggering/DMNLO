      double precision function stsT2QQbar_M2CountertermsProp()
      implicit none

#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"
#include "../../util/DMNLO_RenScheme.h"


      double complex dPropSG,dPropFG
      double precision colorf,Nf
      double complex DenSV

!Propagator
      DenSV = sman

c------only for debugging:
! 	   dZGlR = 0d0
! 	   dZGlL = 0d0
!        dZG = 0d0
	
! !CountertermKop Gluon: 
! !k^2 in Amplitude
      dPropFG = dZG

      dPropSG = 0d0

! Number of effectively massless quarks
      Nf = 4d0
 
! c======================================================================================================================================
! c===============================================Gluon-Propagator Correction==========================================================
! c======================================================================================================================================

! c==========================================================Amplitude Squares=============================================================
! c-----------------------Amplitude square for S-channel with Gluon -----------

      colorf = (Nc**2-1)/4d0
      stst2qqbar_M2CountertermsProp = Nf*2*colorf*dble( -((sqrtS**4*(dPropSG - dPropFG*sman)*
     &      ((-GS)* ((muu3**4 - muu4**2 + muu4**4 + muu3**2*(-1 + 2*muu4**2- 4*tred) + 
     &              4*muu1**2*(muu2**2 - tred) + 4*tred - 4*muu2**2*tred - 
     &              4*muu4**2*tred + 4*tred**2)*(-GS) + 
     &           2*(-1 + 2*muu1**2 + 2*muu2**2)*muu3*muu4*
     &            (-GS)) +  (-GS)*(2*(-1 + 2*muu1**2 + 2*muu2**2)*muu3*muu4*
     &            (-GS) + (muu3**4 - muu4**2 + muu4**4 + muu3**2*(-1 + 2*muu4**2 - 4*tred) + 
     &              4*muu1**2*(muu2**2 - tred) + 4*tred - 4*muu2**2*tred - 
     &              4*muu4**2*tred + 4*tred**2)*(-GS)))*
     &      dconjg(SfSfg(isf1,itt1,igen1,isf2,itt2,igen2))*
     &      SfSfg(isf1,itt1,igen1,isf2,itt2,igen2))/(dconjg(DenSV)*DenSV**2)))
       end
