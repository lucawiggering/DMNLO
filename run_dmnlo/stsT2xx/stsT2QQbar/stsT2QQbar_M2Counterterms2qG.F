      double precision function stsT2QQbar_M2VertexCounter2qG()
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"
#include "../common/stsT2xx_GenCounterterms.h"

      double precision colfac
      double complex DenSV
      integer fType,fGen
	
!Propagator
      DenSV = sman 

! initilize with zero
      stsT2QQbar_M2VertexCounter2qG = 0d0
	     
! Loop over massless quarks
      do fType = 3,4
      do fGen = 1,2

*************************************************
* Counterterms for fermion-fermion-Gluon vertex *
*************************************************

c---------------------Amplitudensquare for s-channel with Gluon-----------
      colfac=(Nc**2-1)/4d0
      
      stsT2QQbar_M2VertexCounter2qG =  stsT2QQbar_M2VertexCounter2qG - 2*colfac*dreal((sqrtS**4
     &     *(BL(fType,fType,fGen,fGen)*((muu3**4 - muu4**2 + muu4**4 + 
     &              muu3**2*(-1 + 2*muu4**2 - 4*tred) + 
     &              4*muu1**2*(muu2**2 - tred) + 4*tred - 4*muu2**2*tred - 
     &              4*muu4**2*tred + 4*tred**2)*dBL(fType,fType,fGen,fGen) + 
     &           2*(-1 + 2*muu1**2 + 2*muu2**2)*muu3*muu4*
     &            dBR(fType,fType,fGen,fGen)) + 
     &        BR(fType,fType,fGen,fGen)*
     &         (2*(-1 + 2*muu1**2 + 2*muu2**2)*muu3*muu4*
     &            dBL(fType,fType,fGen,fGen) + 
     &           (muu3**4 - muu4**2 + muu4**4 + 
     &              muu3**2*(-1 + 2*muu4**2 - 4*tred) + 
     &              4*muu1**2*(muu2**2 - tred) + 4*tred - 4*muu2**2*tred - 
     &              4*muu4**2*tred + 4*tred**2)*dBR(fType,fType,fGen,fGen)))*
     &      dconjg(SfSfg(isf1,itt1,igen1,isf2,itt2,igen2))*
     &      SfSfg(isf1,itt1,igen1,isf2,itt2,igen2))/
     &    (dconjg(DenSV)*DenSV))

      enddo
      enddo

      end

