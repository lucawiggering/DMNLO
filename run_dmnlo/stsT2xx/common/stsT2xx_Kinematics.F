!**********************************
!**** subroutine to set the kinematics of the ingoing 4-momenta and magnitude of
!**** in/out-going momenta
!***********************************

      subroutine stsT2xx_SetKinematics()

      implicit none

#include "stsT2xx_Kinematics.h"
#include "stsT2xx_Model.h"
	
      muu1 = MSf(iSf1,itt1,igen1)/sqrtS
      muu2 = MSf(iSf2,itt2,igen2)/sqrtS		
      muu3 = 0d0
      muu4 = 0d0
      muu5 = 0d0

      m1 = MSf(iSf1,itt1,igen1)

      end
	
