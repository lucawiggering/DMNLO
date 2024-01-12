      double precision function stsT2xx_M2Real3g()
      implicit none
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"
#include "common/stsT2xx_GenCouplings.h"

      double precision M2 

      M2=MSf2(isf1,itt1,igen1)

!************************************************   |Ms4|^2   ************************************************

      stsT2xx_M2Real3g=((-27*GGGG**2*Nc**2*(-1 + Nc**2)*(M2 - p1p2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(4*(M2 + p1p2)**2))

!************************************************   |Mtt4k1|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-((-2 + Nc**2)*(-1 + Nc**2)**2*(k2k3 + M2 - p1k2 - p1k3 - p1p2 
     &+ p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(2*Nc**2*(-k2k3 + p1k2 + p1k3)**2))

!************************************************   |Mtt4k2|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-((-2 + Nc**2)*(-1 + Nc**2)**2*(k1k3 + M2 - p1k1 - p1k3 - p1p2 
     &+ p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(2*Nc**2*(-k1k3 + p1k1 + p1k3)**2))

!************************************************   |Mtt4k3|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-((-2 + Nc**2)*(-1 + Nc**2)**2*(k1k2 + M2 - p1k1 - p1k2 - p1p2 
     &+ p2k1 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(2*Nc**2*(-k1k2 + p1k1 + p1k2)**2))

!************************************************   |Mbt4k1|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-(((-2 + Nc**2)*(-1 + Nc**2)**2*(M2 - p1k1)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(Nc**2*p1k1**2)))

!************************************************   |Mbt4k2|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-(((-2 + Nc**2)*(-1 + Nc**2)**2*(M2 - p1k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(Nc**2*p1k2**2)))

!************************************************   |Mbt4k3|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-(((-2 + Nc**2)*(-1 + Nc**2)**2*(M2 - p1k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(Nc**2*p1k3**2)))

!************************************************   |Mttk1|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+((GGG**2*(-1 + Nc**2)**2*(k2k3**2 - 2*k2k3*(8*M2 + p1k2 + p1k3) 
     &+ 4*(p1k2**2 + 7*p1k2*p1k3 + p1k3**2))*(k2k3 + M2 - p1k2 - p1k3 - p1p2 + p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(16*k2k3**2*(-k2k3 + p1k2 + p1k3)**2))

!************************************************   |Mttk2|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+((GGG**2*(-1 + Nc**2)**2*(k1k3**2 - 2*k1k3*(8*M2 + p1k1 + p1k3) 
     &+ 4*(p1k1**2 + 7*p1k1*p1k3 + p1k3**2))*(k1k3 + M2 - p1k1 - p1k3 - p1p2 + p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(16*k1k3**2*(-k1k3 + p1k1 + p1k3)**2))

!************************************************   |Mttk3|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+((GGG**2*(-1 + Nc**2)**2*(k1k2**2 - 2*k1k2*(8*M2 + p1k1 + p1k2) 
     &+ 4*(p1k1**2 + 7*p1k1*p1k2 + p1k2**2))*(k1k2 + M2 - p1k1 - p1k2 - p1p2 + p2k1 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(16*k1k2**2*(-k1k2 + p1k1 + p1k2)**2))

!************************************************   |Mbtk1|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-(GGG**2*(-1 + Nc**2)**2*(M2 - p1k1)*(-k1k2**2 - k1k3**2 + 8
     &*k2k3*M2 - 8*k2k3*p1k1 - p1k2**2 - 7*p1k2*p1k3 - p1k3**2 - 8*k2k3*p1p2 + 8*k2k3*p2k1 + 2*p1k2
     &*p2k2 + 7*p1k3*p2k2 - p2k2**2 + k1k2*(-7*k1k3 + 2*p1k2 + 7*p1k3 - 2*p2k2 - 7*p2k3) + k1k3*(7
     &*p1k2 + 2*p1k3 - 7*p2k2 - 2*p2k3) + 7*p1k2*p2k3 + 2*p1k3*p2k3 - 7*p2k2*p2k3 - p2k3**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(8*k2k3**2*p1k1**2))

!************************************************   |Mbtk2|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-(GGG**2*(-1 + Nc**2)**2*(M2 - p1k2)*(-k1k2**2 - k2k3**2 + 8
     &*k1k3*M2 - p1k1**2 - 8*k1k3*p1k2 - 7*p1k1*p1k3 - p1k3**2 - 8*k1k3*p1p2 + 2*p1k1*p2k1 + 7*p1k3
     &*p2k1 - p2k1**2 + 8*k1k3*p2k2 + k1k2*(-7*k2k3 + 2*p1k1 + 7*p1k3 - 2*p2k1 - 7*p2k3) + k2k3*(7
     &*p1k1 + 2*p1k3 - 7*p2k1 - 2*p2k3) + 7*p1k1*p2k3 + 2*p1k3*p2k3 - 7*p2k1*p2k3 - p2k3**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(8*k1k3**2*p1k2**2))

!************************************************   |Mbtk3|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-(GGG**2*(-1 + Nc**2)**2*(M2 - p1k3)*(-k1k3**2 - k2k3**2 + 8
     &*k1k2*M2 - p1k1**2 - 7*p1k1*p1k2 - p1k2**2 - 8*k1k2*p1k3 - 8*k1k2*p1p2 + 2*p1k1*p2k1 + 7*p1k2
     &*p2k1 - p2k1**2 + k1k3*(-7*k2k3 + 2*p1k1 + 7*p1k2 - 2*p2k1 - 7*p2k2) + k2k3*(7*p1k1 + 2*p1k2 
     &- 7*p2k1 - 2*p2k2) + 7*p1k1*p2k2 + 2*p1k2*p2k2 - 7*p2k1*p2k2 - p2k2**2 + 8*k1k2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(8*k1k2**2*p1k3**2))

!************************************************   |Msk1|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-(GGG**4*Nc**2*(-1 + Nc**2)*(16*k2k3**2*M2 + 32*k2k3*M2**2 
     &+ 16*k2k3*M2*p1k1 + k2k3*p1k1**2 + 16*k2k3*M2*p1k2 - 19*k2k3*p1k1*p1k2 + 8*k2k3*p1k2**2 - 4*M2
     &*p1k2**2 + 16*k2k3*M2*p1k3 - 19*k2k3*p1k1*p1k3 + 6*k2k3*p1k2*p1k3 - 28*M2*p1k2*p1k3 + 8*k2k3
     &*p1k3**2 - 4*M2*p1k3**2 - 2*k1k2**2*(M2 - p1p2) - 2*k1k3**2*(M2 - p1p2) - 16*k2k3**2*p1p2 - 16
     &*k2k3*p1k1*p1p2 - 16*k2k3*p1k2*p1p2 - 16*k2k3*p1k3*p1p2 - 32*k2k3*p1p2**2 + 16*k2k3*M2*p2k1 
     &- 2*k2k3*p1k1*p2k1 + 19*k2k3*p1k2*p2k1 + 19*k2k3*p1k3*p2k1 - 16*k2k3*p1p2*p2k1 + k2k3*p2k1**2 
     &+ 16*k2k3*M2*p2k2 + 19*k2k3*p1k1*p2k2 - 16*k2k3*p1k2*p2k2 - 6*k2k3*p1k3*p2k2 - 10*p1k2*p1k3
     &*p2k2 + 10*p1k3**2*p2k2 - 16*k2k3*p1p2*p2k2 + 8*p1k2*p1p2*p2k2 + 28*p1k3*p1p2*p2k2 - 19*k2k3*p2k1
     &*p2k2 + 8*k2k3*p2k2**2 - 4*M2*p2k2**2 + 10*p1k3*p2k2**2 + 16*k2k3*M2*p2k3 + 19*k2k3*p1k1*p2k3 
     &- 6*k2k3*p1k2*p2k3 + 10*p1k2**2*p2k3 - 16*k2k3*p1k3*p2k3 - 10*p1k2*p1k3*p2k3 - 16*k2k3*p1p2
     &*p2k3 + 28*p1k2*p1p2*p2k3 + 8*p1k3*p1p2*p2k3 - 19*k2k3*p2k1*p2k3+ 6*k2k3*p2k2*p2k3 - 28*M2
     &*p2k2*p2k3 - 10*p1k2*p2k2*p2k3 - 10*p1k3*p2k2*p2k3 + 8*k2k3*p2k3**2 - 4*M2*p2k3**2 + 10*p1k2
     &*p2k3**2 + k1k3*(14*p1k2*p1p2 + 4*p1k3*p1p2 - 7*p1k2*p2k1 - 2*p1k3*p2k1 + 14*p1p2*p2k2 + 7*p2k1
     &*p2k2 + p1k1*(7*p1k2 + 2*p1k3 - 7*p2k2 - 2*p2k3) + 4*p1p2*p2k3 + 2*p2k1*p2k3 - 2*M2*(7*p1k2 
     &+ 2*p1k3 + 7*p2k2 + 2*p2k3)) + k1k2*(2*p1k1*p1k2 + 7*p1k1*p1k3 - 14*k1k3*(M2 - p1p2) + 4*p1k2
     &*p1p2 + 14*p1k3*p1p2 - 2*p1k2*p2k1 - 7*p1k3*p2k1 - 2*p1k1*p2k2 + 4*p1p2*p2k2 + 2*p2k1*p2k2 - 7
     &*p1k1*p2k3 + 14*p1p2*p2k3 + 7*p2k1*p2k3 - 2*M2*(2*p1k2 + 7*p1k3 + 2*p2k2 + 7*p2k3)))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(16*k2k3**2*(M2 + p1p2)**2))

!************************************************   |Msk2|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-(GGG**4*Nc**2*(-1 + Nc**2)*(-2*k2k3**2*M2 - 14*k2k3*M2*p1k1 
     &- 4*M2*p1k1**2 + 7*k2k3*p1k1*p1k2 - 4*k2k3*M2*p1k3 - 28*M2*p1k1*p1k3 + 2*k2k3*p1k2*p1k3 - 4
     &*M2*p1k3**2 - 2*k1k2**2*(M2 - p1p2) + 16*k1k3**2*(M2 - p1p2) + 2*k2k3**2*p1p2 + 14*k2k3*p1k1
     &*p1p2 + 4*k2k3*p1k3*p1p2 - 14*k2k3*M2*p2k1 - 7*k2k3*p1k2*p2k1 - 10*p1k1*p1k3*p2k1 + 10*p1k3**2
     &*p2k1 + 14*k2k3*p1p2*p2k1 + 8*p1k1*p1p2*p2k1 + 28*p1k3*p1p2*p2k1- 4*M2*p2k1**2 + 10*p1k3*p2k1
     &**2 - 7*k2k3*p1k1*p2k2 - 2*k2k3*p1k3*p2k2 + 7*k2k3*p2k1*p2k2 - 4*k2k3*M2*p2k3 + 10*p1k1**2*p2k3 
     &- 2*k2k3*p1k2*p2k3 - 10*p1k1*p1k3*p2k3 + 4*k2k3*p1p2*p2k3 + 28*p1k1*p1p2*p2k3 + 8*p1k3*p1p2
     &*p2k3 - 28*M2*p2k1*p2k3 - 10*p1k1*p2k1*p2k3 - 10*p1k3*p2k1*p2k3 + 2*k2k3*p2k2*p2k3 - 4*M2*p2k3
     &**2 + 10*p1k1*p2k3**2 + k1k3*(32*M2**2 + 8*p1k1**2 + p1k2**2 - 19*p1k2*p1k3 + 8*p1k3**2 - 16
     &*p1k2*p1p2 - 16*p1k3*p1p2 - 32*p1p2**2 + 19*p1k2*p2k1 - 6*p1k3*p2k1 - 16*p1p2*p2k1 + 8*p2k1
     &**2 - 2*p1k2*p2k2 + 19*p1k3*p2k2 - 16*p1p2*p2k2 - 19*p2k1*p2k2 + p2k2**2 + 19*p1k2*p2k3 - 16
     &*p1k3*p2k3 - 16*p1p2*p2k3 + 6*p2k1*p2k3 - 19*p2k2*p2k3 + 8*p2k3**2+ 16*M2*(p1k1 + p1k2 + p1k3 
     &+ p2k1 + p2k2 + p2k3) - p1k1*(19*p1k2 - 6*p1k3 + 16*p1p2 + 16*p2k1- 19*p2k2 + 6*p2k3)) + k1k2
     &*(2*p1k1*p1k2 + 7*p1k2*p1k3 - 14*k2k3*(M2 - p1p2) + 4*p1k1*p1p2 + 14*p1k3*p1p2 - 2*p1k2*p2k1 
     &+ 4*p1p2*p2k1 - 2*p1k1*p2k2 - 7*p1k3*p2k2 + 2*p2k1*p2k2 - 7*p1k2*p2k3 + 14*p1p2*p2k3 + 7*p2k2
     &*p2k3 - 2*M2*(2*p1k1 + 7*p1k3 +2*p2k1 + 7*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)
     &/(16*k1k3**2*(M2 + p1p2)**2))

!************************************************   |Msk3|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-(GGG**4*Nc**2*(-1 + Nc**2)*(-2*k2k3**2*M2 - 14*k2k3*M2*p1k1 
     &- 4*M2*p1k1**2 - 4*k2k3*M2*p1k2 - 28*M2*p1k1*p1k2 - 4*M2*p1k2**2 +7*k2k3*p1k1*p1k3 + 2*k2k3
     &*p1k2*p1k3 + 16*k1k2**2*(M2 - p1p2) - 2*k1k3**2*(M2 - p1p2) + 2*k2k3**2*p1p2 + 14*k2k3*p1k1
     &*p1p2 + 4*k2k3*p1k2*p1p2 - 14*k2k3*M2*p2k1 - 10*p1k1*p1k2*p2k1 + 10*p1k2**2*p2k1 - 7*k2k3*p1k3
     &*p2k1 + 14*k2k3*p1p2*p2k1 + 8*p1k1*p1p2*p2k1 + 28*p1k2*p1p2*p2k1 -4*M2*p2k1**2 + 10*p1k2*p2k1
     &**2 - 4*k2k3*M2*p2k2 + 10*p1k1**2*p2k2 - 10*p1k1*p1k2*p2k2 - 2*k2k3*p1k3*p2k2 + 4*k2k3*p1p2
     &*p2k2 + 28*p1k1*p1p2*p2k2 + 8*p1k2*p1p2*p2k2 - 28*M2*p2k1*p2k2 - 10*p1k1*p2k1*p2k2 - 10*p1k2*p2k1
     &*p2k2 - 4*M2*p2k2**2 + 10*p1k1*p2k2**2 - 7*k2k3*p1k1*p2k3 - 2*k2k3*p1k2*p2k3 + 7*k2k3*p2k1
     &*p2k3 + 2*k2k3*p2k2*p2k3 + k1k3*(2*p1k1*p1k3 + 7*p1k2*p1k3 - 14*k2k3*(M2 - p1p2) + 4*p1k1*p1p2 
     &+ 14*p1k2*p1p2 - 2*p1k3*p2k1 + 4*p1p2*p2k1 - 7*p1k3*p2k2 + 14*p1p2*p2k2 - 2*M2*(2*p1k1 + 7
     &*p1k2 + 2*p2k1 + 7*p2k2) - 2*p1k1*p2k3 - 7*p1k2*p2k3 + 2*p2k1*p2k3+7*p2k2*p2k3) + k1k2*(32*M2
     &**2 + 8*p1k1**2 + 8*p1k2**2 - 19*p1k2*p1k3 + p1k3**2 - 16*p1k2*p1p2 - 16*p1k3*p1p2 - 32*p1p2
     &**2 - 6*p1k2*p2k1 + 19*p1k3*p2k1 - 16*p1p2*p2k1 + 8*p2k1**2 - 16*p1k2*p2k2 + 19*p1k3*p2k2 - 16
     &*p1p2*p2k2 + 6*p2k1*p2k2 + 8*p2k2**2 + 19*p1k2*p2k3 - 2*p1k3*p2k3 - 16*p1p2*p2k3 - 19*p2k1*p2k3 
     &- 19*p2k2*p2k3 + p2k3**2 + 16*M2*(p1k1 + p1k2 + p1k3 + p2k1 + p2k2+ p2k3) + p1k1*(6*p1k2 
     &- 19*p1k3 - 16*p1p2 - 16*p2k1 - 6*p2k2 + 19*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)
     &/(16*k1k2**2*(M2 + p1p2)**2))

!************************************************   |Mvk1|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+((-9*GGG**2*(2 - 3*Nc**2 + Nc**4)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(4*k2k3))

!************************************************   |Mvk2|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+((-9*GGG**2*(2 - 3*Nc**2 + Nc**4)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(4*k1k3))

!************************************************   |Mvk3|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+((-9*GGG**2*(2 - 3*Nc**2 + Nc**4)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(4*k1k2))

!************************************************   |MtP123|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-((-1 + Nc**2)**3*(M2 - p1k1)*(k1k2 + M2 - 2*p1k1 - p1k2)*(k1k2 
     &+ M2 - p1k1 - p1k2 - p1p2 + p2k1 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(4.
     &*Nc**2*p1k1**2*(-k1k2 + p1k1 + p1k2)**2))

!************************************************   |MtP132|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-((-1 + Nc**2)**3*(M2 - p1k1)*(k1k3 + M2 - 2*p1k1 - p1k3)*(k1k3 
     &+ M2 - p1k1 - p1k3 - p1p2 + p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(4.
     &*Nc**2*p1k1**2*(-k1k3 + p1k1 + p1k3)**2))

!************************************************   |MtP321|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-((-1 + Nc**2)**3*(k2k3 + M2 - p1k2 - 2*p1k3)*(M2 - p1k3)*(k2k3 
     &+ M2 - p1k2 - p1k3 - p1p2 + p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(4.
     &*Nc**2*p1k3**2*(-k2k3 + p1k2 + p1k3)**2))

!************************************************   |MtP231|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-((-1 + Nc**2)**3*(M2 - p1k2)*(k2k3 + M2 - 2*p1k2 - p1k3)*(k2k3 
     &+ M2 - p1k2 - p1k3 - p1p2 + p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(4.
     &*Nc**2*p1k2**2*(-k2k3 + p1k2 + p1k3)**2))

!************************************************   |MtP213|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-((-1 + Nc**2)**3*(k1k2 + M2 - p1k1 - 2*p1k2)*(M2 - p1k2)*(k1k2 
     &+ M2 - p1k1 - p1k2 - p1p2 + p2k1 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(4.
     &*Nc**2*p1k2**2*(-k1k2 + p1k1 + p1k2)**2))

!************************************************   |MtP321|^2   ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+(-((-1 + Nc**2)**3*(k1k3 + M2 - p1k1 - 2*p1k3)*(M2 - p1k3)*(k1k3 
     &+ M2 - p1k1 - p1k3 - p1p2 + p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(4.
     &*Nc**2*p1k3**2*(-k1k3 + p1k1 + p1k3)**2))

!************************************************   2 Re(Ms4 Mtt4k1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGGG*Nc**2*(-1 + Nc**2)*(-2*M2 + p1k2 + p1k3 + 2*p1p2 
     &- p2k2 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(16.
     &*(k2k3 - p1k2 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Ms4 Mtt4k2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGGG*Nc**2*(-1 + Nc**2)*(-2*M2 + p1k1 + p1k3 + 2*p1p2 
     &- p2k1 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(16.
     &*(k1k3 - p1k1 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Ms4 Mtt4k3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGGG*Nc**2*(-1 + Nc**2)*(-2*M2 + p1k1 + p1k2 + 2*p1p2 
     &- p2k1 - p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(16.
     &*(k1k2 - p1k1 - p1k2)*(M2 + p1p2)))

!************************************************   2 Re(Ms4 Mbt4k1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGGG*Nc**2*(-1 + Nc**2)*(2*M2 - p1k1 - 2*p1p2 + p2k1)
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(16*p1k1*(M2 + p1p2)))

!************************************************   2 Re(Ms4 Mbt4k2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGGG*Nc**2*(-1 + Nc**2)*(2*M2 - p1k2 - 2*p1p2 + p2k2)
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(16*p1k2*(M2 + p1p2)))

!************************************************   2 Re(Ms4 Mbt4k3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGGG*Nc**2*(-1 + Nc**2)*(2*M2 - p1k3 - 2*p1p2 + p2k3)
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(16*p1k3*(M2 + p1p2)))

!************************************************   2 Re(Ms4 Mttk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((-9*GGG*GGGG*Nc**2*(-1 + Nc**2)*(2*k2k3*(M2 - p1p2) 
     &- (p1k2 + p1k3)*(p1k2 + p1k3 - p2k2 - p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(32*k2k3
     &*(k2k3 - p1k2 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Ms4 Mttk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((-9*GGG*GGGG*Nc**2*(-1 + Nc**2)*(2*k1k3*(M2 - p1p2) 
     &- (p1k1 + p1k3)*(p1k1 + p1k3 - p2k1 - p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(32*k1k3
     &*(k1k3 - p1k1 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Ms4 Mttk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((-9*GGG*GGGG*Nc**2*(-1 + Nc**2)*(2*k1k2*(M2 - p1p2) 
     &- (p1k1 + p1k2)*(p1k1 + p1k2 - p2k1 - p2k2))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(32*k1k2
     &*(k1k2 - p1k1 - p1k2)*(M2 + p1p2)))

!************************************************   2 Re(Ms4 Mbtk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((9*GGG*GGGG*Nc**2*(-1 + Nc**2)*(2*M2*p1k2 + p1k1*p1k2 
     &+ 2*M2*p1k3 + p1k1*p1k3 - 2*p1k2*p1p2 - 2*p1k3*p1p2 + k1k2*(-2*M2 + p1k1 + 2*p1p2 - p2k1) + k1k3
     &*(-2*M2 + p1k1 + 2*p1p2 - p2k1) + p1k2*p2k1 + p1k3*p2k1 + 2*M2*p2k2 - 3*p1k1*p2k2 - 2*p1p2
     &*p2k2 + p2k1*p2k2 + 2*M2*p2k3 - 3*p1k1*p2k3 - 2*p1p2*p2k3 + p2k1*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3)/(64*k2k3*p1k1*(M2 + p1p2)))

!************************************************   2 Re(Ms4 Mbtk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((9*GGG*GGGG*Nc**2*(-1 + Nc**2)*(2*M2*p1k1 + p1k1*p1k2 
     &+ 2*M2*p1k3 + p1k2*p1k3 - 2*p1k1*p1p2 - 2*p1k3*p1p2 + 2*M2*p2k1 - 3*p1k2*p2k1 - 2*p1p2*p2k1 
     &+ k1k2*(-2*M2 + p1k2 + 2*p1p2 - p2k2) + k2k3*(-2*M2 + p1k2 + 2*p1p2 - p2k2) + p1k1*p2k2 + p1k3
     &*p2k2 + p2k1*p2k2 + 2*M2*p2k3 - 3*p1k2*p2k3 - 2*p1p2*p2k3 + p2k2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3)/(64*k1k3*p1k2*(M2 + p1p2)))

!************************************************   2 Re(Ms4 Mbtk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((9*GGG*GGGG*Nc**2*(-1 + Nc**2)*(2*M2*p1k1 + 2*M2*p1k2 
     &+ p1k1*p1k3 + p1k2*p1k3 - 2*p1k1*p1p2 - 2*p1k2*p1p2 + 2*M2*p2k1 - 3*p1k3*p2k1 - 2*p1p2*p2k1 
     &+ 2*M2*p2k2 - 3*p1k3*p2k2 - 2*p1p2*p2k2 + k1k3*(-2*M2 + p1k3 + 2*p1p2 - p2k3) + k2k3*(-2*M2 
     &+ p1k3 + 2*p1p2 - p2k3) + p1k1*p2k3 + p1k2*p2k3 + p2k1*p2k3 + p2k2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3)/(64*k1k2*p1k3*(M2 + p1p2)))

!************************************************   2 Re(Ms4 Msk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((9*GGG**2*GGGG*Nc**2*(-1 + Nc**2)*(-3*p1k1*p1k2 + 2*p1k2
     &**2 - 3*p1k1*p1k3 + 4*p1k2*p1k3 + 2*p1k3**2 + 4*k2k3*(M2 - p1p2) - 2*p1k2*p1p2 - 2*p1k3*p1p2 
     &+ 3*p1k2*p2k1 + 3*p1k3*p2k1 + 3*p1k1*p2k2 - 4*p1k2*p2k2 - 4*p1k3*p2k2 - 2*p1p2*p2k2 - 3*p2k1
     &*p2k2 + 2*p2k2**2 + 3*p1k1*p2k3 - 4*p1k2*p2k3 - 4*p1k3*p2k3 - 2*p1p2*p2k3 - 3*p2k1*p2k3 + 4
     &*p2k2*p2k3 + 2*p2k3**2 + 2*M2*(p1k2 + p1k3 + p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(32*k2k3*(M2 + p1p2)**2))

!************************************************   2 Re(Ms4 Msk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((9*GGG**2*GGGG*Nc**2*(-1 + Nc**2)*(2*p1k1**2 - 3*p1k1
     &*p1k2 + 4*p1k1*p1k3 - 3*p1k2*p1k3 + 2*p1k3**2 + 4*k1k3*(M2 - p1p2) - 2*p1k1*p1p2 - 2*p1k3*p1p2 
     &- 4*p1k1*p2k1 + 3*p1k2*p2k1 - 4*p1k3*p2k1 - 2*p1p2*p2k1 + 2*p2k1**2 + 3*p1k1*p2k2 + 3*p1k3
     &*p2k2 - 3*p2k1*p2k2 - 4*p1k1*p2k3 + 3*p1k2*p2k3 - 4*p1k3*p2k3 - 2*p1p2*p2k3 + 4*p2k1*p2k3 - 3
     &*p2k2*p2k3 + 2*p2k3**2 + 2*M2*(p1k1 + p1k3 + p2k1 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(32*k1k3*(M2 + p1p2)**2))

!************************************************   2 Re(Ms4 Msk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((9*GGG**2*GGGG*Nc**2*(-1 + Nc**2)*(2*p1k1**2 + 4*p1k1
     &*p1k2 + 2*p1k2**2 - 3*p1k1*p1k3 - 3*p1k2*p1k3 + 4*k1k2*(M2 - p1p2) - 2*p1k1*p1p2 - 2*p1k2*p1p2 
     &- 4*p1k1*p2k1 - 4*p1k2*p2k1 + 3*p1k3*p2k1 - 2*p1p2*p2k1 + 2*p2k1**2 - 4*p1k1*p2k2 - 4*p1k2
     &*p2k2 + 3*p1k3*p2k2 - 2*p1p2*p2k2 + 4*p2k1*p2k2 + 2*p2k2**2 + 2*M2*(p1k1 + p1k2 + p2k1 + p2k2) 
     &+ 3*p1k1*p2k3 + 3*p1k2*p2k3 - 3*p2k1*p2k3 - 3*p2k2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(32*k1k2*(M2 + p1p2)**2))

!************************************************   2 Re(Ms4 Mvk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Ms4 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Ms4 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Ms4 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGGG*Nc**2*(-1 + Nc**2)*(-4*M2**2 + 2*p1k1**2 + 3*p1k1
     &*p1k2 + 2*p1k2**2 - 2*p1k1*p1p2 + 2*p1k2*p1p2 + 4*p1p2**2 - 4*p1k1*p2k1 - 5*p1k2*p2k1 - 6
     &*p1p2*p2k1 + 2*p2k1**2 + 2*M2*(3*p1k1 + p1k2 + p2k1 - p2k2) + p1k1*p2k2 - 2*p1k2*p2k2 - 2*p1p2
     &*p2k2 + p2k1*p2k2 + k1k2*(-4*M2 - p1k2 + 4*p1p2 + p2k2))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(64*p1k1*(-k1k2 + p1k1 + p1k2)*(M2 + p1p2)))

!************************************************   2 Re(Ms4 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGGG*Nc**2*(-1 + Nc**2)*(-4*M2**2 + 2*p1k1**2 + 3*p1k1
     &*p1k3 + 2*p1k3**2 - 2*p1k1*p1p2 + 2*p1k3*p1p2 + 4*p1p2**2 - 4*p1k1*p2k1 - 5*p1k3*p2k1 - 6
     &*p1p2*p2k1 + 2*p2k1**2 + 2*M2*(3*p1k1 + p1k3 + p2k1 - p2k3) + p1k1*p2k3 - 2*p1k3*p2k3 - 2*p1p2
     &*p2k3 + p2k1*p2k3 + k1k3*(-4*M2 - p1k3 + 4*p1p2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(64*p1k1*(-k1k3 + p1k1 + p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Ms4 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGGG*Nc**2*(-1 + Nc**2)*(-4*M2**2 + 2*p1k2**2 + 3*p1k2
     &*p1k3 + 2*p1k3**2 + 2*p1k2*p1p2 - 2*p1k3*p1p2 + 4*p1p2**2 - 2*p1k2*p2k2 + p1k3*p2k2 - 2*p1p2
     &*p2k2 + k2k3*(-4*M2 - p1k2 + 4*p1p2 + p2k2) - 5*p1k2*p2k3 - 4*p1k3*p2k3 - 6*p1p2*p2k3 + p2k2
     &*p2k3 + 2*p2k3**2 + 2*M2*(p1k2 + 3*p1k3 - p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(64*p1k3*(-k2k3 + p1k2 + p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Ms4 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGGG*Nc**2*(-1 + Nc**2)*(-4*M2**2 + 2*p1k2**2 + 3*p1k2
     &*p1k3 + 2*p1k3**2 - 2*p1k2*p1p2 + 2*p1k3*p1p2 + 4*p1p2**2 - 4*p1k2*p2k2 - 5*p1k3*p2k2 - 6
     &*p1p2*p2k2 + 2*p2k2**2 + 2*M2*(3*p1k2 + p1k3 + p2k2 - p2k3) + p1k2*p2k3 - 2*p1k3*p2k3 - 2*p1p2
     &*p2k3 + p2k2*p2k3 + k2k3*(-4*M2 - p1k3 + 4*p1p2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(64*p1k2*(-k2k3 + p1k2 + p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Ms4 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGGG*Nc**2*(-1 + Nc**2)*(-4*M2**2 + 2*p1k1**2 + 3*p1k1
     &*p1k2 + 2*p1k2**2 + 2*p1k1*p1p2 - 2*p1k2*p1p2 + 4*p1p2**2 - 2*p1k1*p2k1 + p1k2*p2k1 - 2*p1p2
     &*p2k1 + k1k2*(-4*M2 - p1k1 + 4*p1p2 + p2k1) - 5*p1k1*p2k2 - 4*p1k2*p2k2 - 6*p1p2*p2k2 + p2k1
     &*p2k2 + 2*p2k2**2 + 2*M2*(p1k1 + 3*p1k2 - p2k1 + p2k2))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(64*p1k2*(-k1k2 + p1k1 + p1k2)*(M2 + p1p2)))

!************************************************   2 Re(Ms4 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGGG*Nc**2*(-1 + Nc**2)*(-4*M2**2 + 2*p1k1**2 + 3*p1k1
     &*p1k3 + 2*p1k3**2 + 2*p1k1*p1p2 - 2*p1k3*p1p2 + 4*p1p2**2 - 2*p1k1*p2k1 + p1k3*p2k1 - 2*p1p2
     &*p2k1 + k1k3*(-4*M2 - p1k1 + 4*p1p2 + p2k1) - 5*p1k1*p2k3 - 4*p1k3*p2k3 - 6*p1p2*p2k3 + p2k1
     &*p2k3 + 2*p2k3**2 + 2*M2*(p1k1 + 3*p1k3 - p2k1 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(64*p1k3*(-k1k3 + p1k1 + p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Mtt4k1 Mtt4k2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-1 + Nc**2)*(k1k2 + k1k3 + k2k3 + 2*M2 - p1k1 - p1k2 
     &- 2*p1k3 - 2*p1p2 + p2k1 + p2k2 + 2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(8*Nc**2*(k2k3 - p1k2 - p1k3)*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mtt4k1 Mtt4k3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-1 + Nc**2)*(k1k2 + k1k3 + k2k3 + 2*M2 - p1k1 - 2*p1k2 
     &- p1k3 - 2*p1p2 + p2k1 + 2*p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(8*Nc**2*(k1k2 - p1k1 - p1k2)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mtt4k1 Mbt4k1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-1 + Nc**2)*(2 + Nc**2)*(k1k2 + k1k3 + 2*M2 - p1k1 
     &- 2*p1k2 - 2*p1k3 - 2*p1p2 + p2k1)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(4*Nc**2*p1k1*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(Mtt4k1 Mbt4k2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2)**2*(-1 + Nc**2)*(k2k3 + 2*M2 - 3*p1k2 
     &- 2*p1k3 - 2*p1p2 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(32*Nc**2*p1k2*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mtt4k1 Mbt4k3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-2 + Nc**2)**2*(-1 + Nc**2)*(k2k3 + 2*M2 - 2*p1k2 
     &- 3*p1k3 - 2*p1p2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(32*Nc**2*(k2k3 - p1k2 - p1k3)*p1k3))

!************************************************   2 Re(Mtt4k1 Mttk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mtt4k1 Mttk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k3**2 - 4*k2k3*M2 + 6*k2k3*p1k1 
     &+ 2*M2*p1k1 - 3*p1k1**2 + 2*p1k1*p1k2 + 3*k2k3*p1k3 - 2*M2*p1k3 - 4*p1k1*p1k3 + 4*p1k2*p1k3 
     &+ p1k3**2 + 4*k2k3*p1p2 - 2*p1k1*p1p2 + 2*p1k3*p1p2 + 2*M2*p2k1 + 2*p1k1*p2k1 - 4*p1k2*p2k1 - 4
     &*p1k3*p2k1 - 2*p1p2*p2k1 + p2k1**2 + k1k2*(k1k3 - 2*M2 + 3*p1k1 + 2*p1p2 + p2k1) + 2*p1k1*p2k2 
     &- 2*p1k3*p2k2 - k2k3*p2k3 - 2*M2*p2k3 + 8*p1k1*p2k3 - 2*p1k2*p2k3 + 2*p1p2*p2k3 - p2k3**2 
     &- k1k3*(k2k3 - 4*M2 - 2*p1k1 + 6*p1k2 + 5*p1k3 + 4*p1p2 - 2*p2k1 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k1k3*(k1k3 - p1k1 - p1k3)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mtt4k1 Mttk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2**2 - 4*k2k3*M2 + 6*k2k3*p1k1 
     &+ 2*M2*p1k1 - 3*p1k1**2 + 3*k2k3*p1k2 - 2*M2*p1k2 - 4*p1k1*p1k2 + p1k2**2 + 2*p1k1*p1k3 + 4*p1k2
     &*p1k3 + 4*k2k3*p1p2 - 2*p1k1*p1p2 + 2*p1k2*p1p2 + 2*M2*p2k1 + 2*p1k1*p2k1 - 4*p1k2*p2k1 - 4
     &*p1k3*p2k1 - 2*p1p2*p2k1 + p2k1**2 + k1k3*(-2*M2 + 3*p1k1 + 2*p1p2 + p2k1) + k1k2*(k1k3 - k2k3 
     &+ 4*M2 + 2*p1k1 - 5*p1k2 - 6*p1k3 - 4*p1p2 + 2*p2k1 - p2k2) - k2k3*p2k2 - 2*M2*p2k2 + 8*p1k1
     &*p2k2 - 2*p1k3*p2k2 + 2*p1p2*p2k2 - p2k2**2 + 2*p1k1*p2k3 - 2*p1k2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k1k2*(k1k2 - p1k1 - p1k2)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mtt4k1 Mbtk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mtt4k1 Mbtk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(2 - 3*Nc**2 + Nc**4)*(2*k2k3**2 - 2*k1k3*M2 - 4
     &*M2*p1k1 + 3*k1k3*p1k2 + 2*p1k1*p1k2 - 2*M2*p1k3 + 2*p1k1*p1k3 + p1k2*p1k3 + 4*p1k3**2 + 2*k1k3
     &*p1p2 + 4*p1k1*p1p2 + 2*p1k3*p1p2 - 4*M2*p2k1 + 6*p1k2*p2k1 + 2*p1k3*p2k1 + 4*p1p2*p2k1 - k1k3
     &*p2k2 - 6*p1k1*p2k2 - 3*p1k3*p2k2 - 2*p2k1*p2k2 - 2*M2*p2k3 - 4*p1k1*p2k3 + 3*p1k2*p2k3 - 4
     &*p1k3*p2k3 + 2*p1p2*p2k3 - p2k2*p2k3 + k2k3*(-3*p1k1 + p1k2 - 6*p1k3 - p2k1 + p2k2 + 2*p2k3) 
     &+ k1k2*(k2k3 + 2*(p1k2 + p2k2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*k1k3*p1k2*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mtt4k1 Mbtk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(2 - 3*Nc**2 + Nc**4)*(2*k2k3**2 - 2*k1k2*M2 - 4
     &*M2*p1k1 - 2*M2*p1k2 + 2*p1k1*p1k2 + 4*p1k2**2 + 3*k1k2*p1k3 + 2*p1k1*p1k3 + p1k2*p1k3 + 2*k1k2
     &*p1p2 + 4*p1k1*p1p2 + 2*p1k2*p1p2 - 4*M2*p2k1 + 2*p1k2*p2k1 + 6*p1k3*p2k1 + 4*p1p2*p2k1 - 2
     &*M2*p2k2 - 4*p1k1*p2k2 - 4*p1k2*p2k2 + 3*p1k3*p2k2 + 2*p1p2*p2k2 - k1k2*p2k3 - 6*p1k1*p2k3 - 3
     &*p1k2*p2k3 - 2*p2k1*p2k3 - p2k2*p2k3 + k2k3*(-3*p1k1 - 6*p1k2 + p1k3 - p2k1 + 2*p2k2 + p2k3) 
     &+ k1k3*(k2k3 + 2*(p1k3 + p2k2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*k1k2*p1k3*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mtt4k1 Msk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mtt4k1 Msk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(-2*k2k3*M2 - k2k3*p1k1 + 6
     &*M2*p1k1 + p1k1**2 + 5*k2k3*p1k2 - 4*p1k1*p1k2 - 3*k2k3*p1k3 + 6*p1k1*p1k3 - 5*p1k2*p1k3 + 5
     &*p1k3**2 + 2*k2k3*p1p2 - 6*p1k1*p1p2 + k2k3*p2k1 + 6*M2*p2k1 - 2*p1k1*p2k1 - 2*p1k2*p2k1 - 12
     &*p1k3*p2k1 - 6*p1p2*p2k1 + p2k1**2 - 5*k2k3*p2k2 + 10*p1k1*p2k2 + 5*p1k3*p2k2 - 4*p2k1*p2k2 
     &+ 3*k2k3*p2k3 + 5*p1k2*p2k3 - 10*p1k3*p2k3 + 6*p2k1*p2k3 - 5*p2k2*p2k3 + 5*p2k3**2 + k1k2*(2*M2 
     &- p1k1 + p1k2 - 4*p1k3 - 2*p1p2 + p2k1 - p2k2 + 4*p2k3) + k1k3*(6*M2 - p1k1 - p1k2 - 6*p1k3 
     &- 6*p1p2 + p2k1 + p2k2 + 6*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*k1k3*(k2k3 - p1k2 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Mtt4k1 Msk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(-2*k2k3*M2 - k2k3*p1k1 + 6
     &*M2*p1k1 + p1k1**2 - 3*k2k3*p1k2 + 6*p1k1*p1k2 + 5*p1k2**2 + 5*k2k3*p1k3 - 4*p1k1*p1k3 - 5*p1k2
     &*p1k3 + 2*k2k3*p1p2 - 6*p1k1*p1p2 + k2k3*p2k1 + 6*M2*p2k1 - 2*p1k1*p2k1 - 12*p1k2*p2k1 - 2
     &*p1k3*p2k1 - 6*p1p2*p2k1 + p2k1**2 + 3*k2k3*p2k2 - 10*p1k2*p2k2 + 5*p1k3*p2k2 + 6*p2k1*p2k2 
     &+ 5*p2k2**2 + k1k3*(2*M2 - p1k1 - 4*p1k2 + p1k3 - 2*p1p2 + p2k1 + 4*p2k2 - p2k3) - 5*k2k3*p2k3 
     &+ 10*p1k1*p2k3 + 5*p1k2*p2k3 - 4*p2k1*p2k3 - 5*p2k2*p2k3 + k1k2*(6*M2 - p1k1 - 6*p1k2 - p1k3 
     &- 6*p1p2 + p2k1 + 6*p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*k1k2*(k2k3 - p1k2 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Mtt4k1 Mvk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mtt4k1 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGG*(4 - 5*Nc**2 + Nc**4)*(k1k2 + k1k3 + 2*k2k3 - p1k1 
     &- 2*p1k3 + p2k1 + 2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(32*k1k3*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(Mtt4k1 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGG*(4 - 5*Nc**2 + Nc**4)*(k1k2 + k1k3 + 2*k2k3 - p1k1 
     &- 2*p1k2 + p2k1 + 2*p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(32*k1k2*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(Mtt4k1 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2 + Nc**4)*(k1k2 + k1k3 + 2*M2 - p1k1 - 2
     &*p1k2 - 2*p1k3 - 2*p1p2 + p2k1)*(3*k1k2 + 2*M2 - 4*p1k1 - 3*p1k2 - 2*p1p2 + 2*p2k1 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*Nc**2*p1k1
     &*(-k1k2 + p1k1 + p1k2)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mtt4k1 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2 + Nc**4)*(k1k2 + k1k3 + 2*M2 - p1k1 - 2
     &*p1k2 - 2*p1k3 - 2*p1p2 + p2k1)*(3*k1k3 + 2*M2 - 4*p1k1 - 3*p1k3 - 2*p1p2 + 2*p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*Nc**2*p1k1
     &*(-k1k3 + p1k1 + p1k3)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mtt4k1 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2)*(-1 + Nc**2)**2*(k2k3 + 4*M2 - 2*(p1k2 
     &+ 3*p1k3))*(k2k3 + M2 - p1k2 - p1k3 - p1p2 + p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*Nc**2*p1k3*(-k2k3 + p1k2 + p1k3)**2))

!************************************************   2 Re(Mtt4k1 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2)*(-1 + Nc**2)**2*(k2k3 + 4*M2 - 2*(3*p1k2 
     &+ p1k3))*(k2k3 + M2 - p1k2 - p1k3 - p1p2 + p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*Nc**2*p1k2*(-k2k3 + p1k2 + p1k3)**2))

!************************************************   2 Re(Mtt4k1 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((2 - 3*Nc**2 + Nc**4)*(k1k2 + 2*M2 - 2*p1k1 - 3*p1k2 
     &- 2*p1p2 + p2k2)*(k1k2 + k1k3 + 2*k2k3 + 2*M2 - p1k1 - 4*p1k2 - 2*p1k3 - 2*p1p2 + p2k1 + 2*p2k2)
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*Nc**2*p1k2*(-k1k2 + p1k1 + p1k2)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mtt4k1 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((2 - 3*Nc**2 + Nc**4)*(k1k3 + 2*M2 - 2*p1k1 - 3*p1k3 
     &- 2*p1p2 + p2k3)*(k1k2 + k1k3 + 2*k2k3 + 2*M2 - p1k1 - 2*p1k2 - 4*p1k3 - 2*p1p2 + p2k1 + 2*p2k3)
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*Nc**2*p1k3*(-k1k3 + p1k1 + p1k3)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mtt4k2 Mtt4k3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-1 + Nc**2)*(k1k2 + k1k3 + k2k3 + 2*M2 - 2*p1k1 - p1k2 
     &- p1k3 - 2*p1p2 + 2*p2k1 + p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(8*Nc**2*(k1k2 - p1k1 - p1k2)*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mtt4k2 Mbt4k1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2)**2*(-1 + Nc**2)*(k1k3 + 2*M2 - 3*p1k1 
     &- 2*p1k3 - 2*p1p2 + p2k1)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(32*Nc**2*p1k1*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mtt4k2 Mbt4k2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-1 + Nc**2)*(2 + Nc**2)*(k1k2 + k2k3 + 2*M2 - 2*p1k1 
     &- p1k2 - 2*p1k3 - 2*p1p2 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(4*Nc**2*p1k2*(k1k3 - p1k1 - p1k3)))

!************************************************   2 Re(Mtt4k2 Mbt4k3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-2 + Nc**2)**2*(-1 + Nc**2)*(k1k3 + 2*M2 - 2*p1k1 
     &- 3*p1k3 - 2*p1p2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(32*Nc**2*(k1k3 - p1k1 - p1k3)*p1k3))

!************************************************   2 Re(Mtt4k2 Mttk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k2k3**2 + 4*k2k3*M2 - 6*k2k3*p1k1 
     &+ 2*k2k3*p1k2 + 2*M2*p1k2 + 2*p1k1*p1k2 - 3*p1k2**2 - 5*k2k3*p1k3 - 2*M2*p1k3 + 4*p1k1*p1k3 
     &- 4*p1k2*p1k3 + p1k3**2 - 4*k2k3*p1p2 - 2*p1k2*p1p2 + 2*p1k3*p1p2 + 2*p1k2*p2k1 - 2*p1k3*p2k1 
     &+ 2*k2k3*p2k2 + 2*M2*p2k2 - 4*p1k1*p2k2 + 2*p1k2*p2k2 - 4*p1k3*p2k2 - 2*p1p2*p2k2 + p2k2**2 
     &+ k1k2*(k2k3 - 2*M2 + 3*p1k2 + 2*p1p2 + p2k2) - k2k3*p2k3 - 2*M2*p2k3 - 2*p1k1*p2k3 + 8*p1k2*p2k3 
     &+ 2*p1p2*p2k3 - p2k3**2 - k1k3*(k2k3 + 4*M2 - 6*p1k2 - 3*p1k3 - 4*p1p2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k2k3*(k2k3 - p1k2 - p1k3)*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mtt4k2 Mttk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mtt4k2 Mttk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2**2 - 2*k2k3*M2 - 2*M2*p1k1 + p1k1
     &**2 + 3*k2k3*p1k2 + 2*M2*p1k2 - 4*p1k1*p1k2 - 3*p1k2**2 + 4*p1k1*p1k3 + 2*p1k2*p1k3 + 2*k2k3
     &*p1p2 + 2*p1k1*p1p2 - 2*p1k2*p1p2 + k1k3*(-4*M2 + 3*p1k1 + 6*p1k2 + 4*p1p2 - p2k1) - 2*M2
     &*p2k1 + 8*p1k2*p2k1 - 2*p1k3*p2k1 + 2*p1p2*p2k1 - p2k1**2 - k1k2*(k1k3 - k2k3 - 4*M2 + 5*p1k1 
     &- 2*p1k2 + 6*p1k3 + 4*p1p2 + p2k1 - 2*p2k2) + k2k3*p2k2 + 2*M2*p2k2 - 4*p1k1*p2k2 + 2*p1k2*p2k2 
     &- 4*p1k3*p2k2 - 2*p1p2*p2k2 + p2k2**2 - 2*p1k1*p2k3 + 2*p1k2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k1k2*(k1k2 - p1k1 - p1k2)*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mtt4k2 Mbtk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(2 - 3*Nc**2 + Nc**4)*(2*k1k3**2 - 2*k2k3*M2 + 3
     &*k2k3*p1k1 - 4*M2*p1k2 + 2*p1k1*p1k2 - 2*M2*p1k3 + p1k1*p1k3 + 2*p1k2*p1k3 + 4*p1k3**2 + 2*k2k3
     &*p1p2 + 4*p1k2*p1p2 + 2*p1k3*p1p2 - k2k3*p2k1 - 6*p1k2*p2k1 - 3*p1k3*p2k1 - 4*M2*p2k2 + 6*p1k1
     &*p2k2 + 2*p1k3*p2k2 + 4*p1p2*p2k2 - 2*p2k1*p2k2 - 2*M2*p2k3 + 3*p1k1*p2k3 - 4*p1k2*p2k3 - 4
     &*p1k3*p2k3 + 2*p1p2*p2k3 - p2k1*p2k3 + k1k3*(p1k1 - 3*p1k2 - 6*p1k3 + p2k1 - p2k2 + 2*p2k3) 
     &+ k1k2*(k1k3 + 2*(p1k1 + p2k1 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*k2k3*p1k1*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mtt4k2 Mbtk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mtt4k2 Mbtk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(2 - 3*Nc**2 + Nc**4)*(2*k1k3**2 - 2*M2*p1k1 + 4
     &*p1k1**2 - 4*M2*p1k2 + 2*p1k1*p1k2 + 2*k2k3*p1k3 + p1k1*p1k3 + 2*p1k2*p1k3 + 2*p1k1*p1p2 + 4
     &*p1k2*p1p2 + 2*k2k3*p2k1 - 2*M2*p2k1 - 4*p1k1*p2k1 - 4*p1k2*p2k1 + 3*p1k3*p2k1 + 2*p1p2*p2k1 
     &- 4*M2*p2k2 + 2*p1k1*p2k2 + 6*p1k3*p2k2 + 4*p1p2*p2k2 + 2*k2k3*p2k3 - 3*p1k1*p2k3 - 6*p1k2*p2k3 
     &- p2k1*p2k3 - 2*p2k2*p2k3 - k1k2*(2*M2 - 3*p1k3 - 2*p1p2 + p2k3) + k1k3*(k2k3 - 6*p1k1 - 3
     &*p1k2 + p1k3 + 2*p2k1 - p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*k1k2*p1k3*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mtt4k2 Msk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(6*k2k3*M2 - k2k3*p1k1 - k2k3
     &*p1k2 + 6*M2*p1k2 - 4*p1k1*p1k2 + p1k2**2 - 6*k2k3*p1k3 - 5*p1k1*p1k3 + 6*p1k2*p1k3 + 5*p1k3
     &**2 - 6*k2k3*p1p2 - 6*p1k2*p1p2 + k2k3*p2k1 + 10*p1k2*p2k1 + 5*p1k3*p2k1 + k2k3*p2k2 + 6*M2
     &*p2k2 - 2*p1k1*p2k2 - 2*p1k2*p2k2 - 12*p1k3*p2k2 - 6*p1p2*p2k2 - 4*p2k1*p2k2 + p2k2**2 + 6*k2k3
     &*p2k3 + 5*p1k1*p2k3 - 10*p1k3*p2k3 - 5*p2k1*p2k3 + 6*p2k2*p2k3 + 5*p2k3**2 + k1k3*(-2*M2 + 5
     &*p1k1 - p1k2 - 3*p1k3 + 2*p1p2 - 5*p2k1 + p2k2 + 3*p2k3) + k1k2*(2*M2 + p1k1 - p1k2 - 4*p1k3 
     &- 2*p1p2 - p2k1 + p2k2 + 4*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*k2k3*(k1k3 - p1k1 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Mtt4k2 Msk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mtt4k2 Msk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(2*k2k3*M2 - 4*k2k3*p1k1 
     &+ 5*p1k1**2 - k2k3*p1k2 + 6*M2*p1k2 + 6*p1k1*p1k2 + p1k2**2 + k2k3*p1k3 - 5*p1k1*p1k3 - 4*p1k2
     &*p1k3 - 2*k2k3*p1p2 - 6*p1k2*p1p2 + 4*k2k3*p2k1 - 10*p1k1*p2k1 + 5*p1k3*p2k1 + 5*p2k1**2 + k2k3
     &*p2k2 + 6*M2*p2k2 - 12*p1k1*p2k2 - 2*p1k2*p2k2 - 2*p1k3*p2k2 - 6*p1p2*p2k2 + 6*p2k1*p2k2 + p2k2
     &**2 + k1k3*(-2*M2 - 3*p1k1 - p1k2 + 5*p1k3 + 2*p1p2 + 3*p2k1 + p2k2 - 5*p2k3) - k2k3*p2k3 
     &+ 5*p1k1*p2k3 + 10*p1k2*p2k3 - 5*p2k1*p2k3 - 4*p2k2*p2k3 + k1k2*(6*M2 - 6*p1k1 - p1k2 - p1k3 
     &- 6*p1p2 + 6*p2k1 + p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*k1k2*(k1k3 - p1k1 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Mtt4k2 Mvk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGG*(4 - 5*Nc**2 + Nc**4)*(k1k2 + 2*k1k3 + k2k3 - p1k2 
     &- 2*p1k3 + p2k2 + 2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(32*k2k3*(k1k3 - p1k1 - p1k3)))

!************************************************   2 Re(Mtt4k2 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mtt4k2 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGG*(4 - 5*Nc**2 + Nc**4)*(k1k2 + 2*k1k3 + k2k3 - 2
     &*p1k1 - p1k2 + 2*p2k1 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(32*k1k2*(k1k3 - p1k1 - p1k3)))

!************************************************   2 Re(Mtt4k2 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((2 - 3*Nc**2 + Nc**4)*(k1k2 + 2*M2 - 3*p1k1 - 2*p1k2 
     &- 2*p1p2 + p2k1)*(k1k2 + 2*k1k3 + k2k3 + 2*M2 - 4*p1k1 - p1k2 - 2*p1k3 - 2*p1p2 + 2*p2k1 + p2k2)
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*Nc**2*p1k1*(-k1k2 + p1k1 + p1k2)*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mtt4k2 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2)*(-1 + Nc**2)**2*(k1k3 + 4*M2 - 2*(3*p1k1 
     &+ p1k3))*(k1k3 + M2 - p1k1 - p1k3 - p1p2 + p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*Nc**2*p1k1*(-k1k3 + p1k1 + p1k3)**2))

!************************************************   2 Re(Mtt4k2 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((2 - 3*Nc**2 + Nc**4)*(k2k3 + 2*M2 - 2*p1k2 - 3*p1k3 
     &- 2*p1p2 + p2k3)*(k1k2 + 2*k1k3 + k2k3 + 2*M2 - 2*p1k1 - p1k2 - 4*p1k3 - 2*p1p2 + p2k2 + 2*p2k3)
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*Nc**2*p1k3*(-k1k3 + p1k1 + p1k3)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mtt4k2 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-2 + Nc**2 + Nc**4)*(k1k2 + k2k3 + 2*M2 - 2*p1k1 - p1k2 
     &- 2*p1k3 - 2*p1p2 + p2k2)*(3*k2k3 + 2*M2 - 4*p1k2 - 3*p1k3 - 2*p1p2 + 2*p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*Nc**2*p1k2
     &*(k1k3 - p1k1 - p1k3)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mtt4k2 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2 + Nc**4)*(k1k2 + k2k3 + 2*M2 - 2*p1k1 
     &- p1k2 - 2*p1k3 - 2*p1p2 + p2k2)*(3*k1k2 + 2*M2 - 3*p1k1 - 4*p1k2 - 2*p1p2 + p2k1 + 2*p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*Nc**2*p1k2
     &*(-k1k2 + p1k1 + p1k2)*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mtt4k2 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2)*(-1 + Nc**2)**2*(k1k3 + 4*M2 - 2*(p1k1 
     &+ 3*p1k3))*(k1k3 + M2 - p1k1 - p1k3 - p1p2 + p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*Nc**2*p1k3*(-k1k3 + p1k1 + p1k3)**2))

!************************************************   2 Re(Mtt4k3 Mbt4k1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2)**2*(-1 + Nc**2)*(k1k2 + 2*M2 - 3*p1k1 
     &- 2*p1k2 - 2*p1p2 + p2k1)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(32*Nc**2*p1k1*(-k1k2 + p1k1 + p1k2)))

!************************************************   2 Re(Mtt4k3 Mbt4k2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-2 + Nc**2)**2*(-1 + Nc**2)*(k1k2 + 2*M2 - 2*p1k1 
     &- 3*p1k2 - 2*p1p2 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(32*Nc**2*(k1k2 - p1k1 - p1k2)*p1k2))

!************************************************   2 Re(Mtt4k3 Mbt4k3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-1 + Nc**2)*(2 + Nc**2)*(k1k3 + k2k3 + 2*M2 - 2*p1k1 
     &- 2*p1k2 - p1k3 - 2*p1p2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(4*Nc**2*(k1k2 - p1k1 - p1k2)*p1k3))

!************************************************   2 Re(Mtt4k3 Mttk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(-1 + Nc**2)*(k2k3**2 + 4*k2k3*M2 - 6*k2k3*p1k1 
     &- 5*k2k3*p1k2 - 2*M2*p1k2 + 4*p1k1*p1k2 + p1k2**2 + 2*k2k3*p1k3 + 2*M2*p1k3 + 2*p1k1*p1k3 - 4
     &*p1k2*p1k3 - 3*p1k3**2 - 4*k2k3*p1p2 + 2*p1k2*p1p2 - 2*p1k3*p1p2 - 2*p1k2*p2k1 + 2*p1k3*p2k1 
     &- k2k3*p2k2 - 2*M2*p2k2 - 2*p1k1*p2k2 + 8*p1k3*p2k2 + 2*p1p2*p2k2 - p2k2**2 - k1k2*(k2k3 + 4
     &*M2 - 3*p1k2 - 6*p1k3 - 4*p1p2 + p2k2) + 2*k2k3*p2k3 + 2*M2*p2k3 - 4*p1k1*p2k3 - 4*p1k2*p2k3 
     &+ 2*p1k3*p2k3 - 2*p1p2*p2k3 + p2k3**2 + k1k3*(k2k3 - 2*M2 + 3*p1k3 + 2*p1p2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k2k3*(k1k2 - p1k1 - p1k2)*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(Mtt4k3 Mttk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(-1 + Nc**2)*(k1k3**2 - 2*k2k3*M2 - 2*M2*p1k1 
     &+ p1k1**2 + 4*p1k1*p1k2 + 3*k2k3*p1k3 + 2*M2*p1k3 - 4*p1k1*p1k3 + 2*p1k2*p1k3 - 3*p1k3**2 + 2
     &*k2k3*p1p2 + 2*p1k1*p1p2 - 2*p1k3*p1p2 - 2*M2*p2k1 - 2*p1k2*p2k1 + 8*p1k3*p2k1 + 2*p1p2*p2k1 
     &- p2k1**2 - k1k2*(k1k3 + 4*M2 - 3*p1k1 - 6*p1k3 - 4*p1p2 + p2k1) - 2*p1k1*p2k2 + 2*p1k3*p2k2 
     &+ k2k3*p2k3 + 2*M2*p2k3 - 4*p1k1*p2k3 - 4*p1k2*p2k3 + 2*p1k3*p2k3 - 2*p1p2*p2k3 + p2k3**2 + k1k3
     &*(k2k3 + 4*M2 - 5*p1k1 - 6*p1k2 + 2*p1k3 - 4*p1p2 - p2k1 + 2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k1k3*(k1k2 - p1k1 - p1k2)*(k1k3 - p1k1 - p1k3)))

!************************************************   2 Re(Mtt4k3 Mttk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mtt4k3 Mbtk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(2 - 3*Nc**2 + Nc**4)*(2*k1k2**2 + 2*k1k3*p1k1 
     &- 2*M2*p1k2 + p1k1*p1k2 + 4*p1k2**2 - 4*M2*p1k3 + 2*p1k1*p1k3 + 2*p1k2*p1k3 + 2*p1k2*p1p2 + 4
     &*p1k3*p1p2 + 2*k1k3*p2k1 - 3*p1k2*p2k1 - 6*p1k3*p2k1 - k2k3*(2*M2 - 3*p1k1 - 2*p1p2 + p2k1) + 2
     &*k1k3*p2k2 - 2*M2*p2k2 + 3*p1k1*p2k2 - 4*p1k2*p2k2 - 4*p1k3*p2k2 + 2*p1p2*p2k2 - p2k1*p2k2 
     &+ k1k2*(k1k3 + p1k1 - 6*p1k2 - 3*p1k3 + p2k1 + 2*p2k2 - p2k3) - 4*M2*p2k3 + 6*p1k1*p2k3 + 2*p1k2
     &*p2k3 + 4*p1p2*p2k3 - 2*p2k1*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*k2k3*p1k1*(-k1k2 + p1k1 + p1k2)))

!************************************************   2 Re(Mtt4k3 Mbtk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(2 - 3*Nc**2 + Nc**4)*(2*k1k2**2 - 2*M2*p1k1 + 4
     &*p1k1**2 + 2*k2k3*p1k2 + p1k1*p1k2 - 4*M2*p1k3 + 2*p1k1*p1k3 + 2*p1k2*p1k3 + 2*p1k1*p1p2 + 4
     &*p1k3*p1p2 + 2*k2k3*p2k1 - 2*M2*p2k1 - 4*p1k1*p2k1 + 3*p1k2*p2k1 - 4*p1k3*p2k1 + 2*p1p2*p2k1 
     &+ 2*k2k3*p2k2 - 3*p1k1*p2k2 - 6*p1k3*p2k2 - p2k1*p2k2 - k1k3*(2*M2 - 3*p1k2 - 2*p1p2 + p2k2) 
     &+ k1k2*(k2k3 - 6*p1k1 + p1k2 - 3*p1k3 + 2*p2k1 + p2k2 - p2k3) - 4*M2*p2k3 + 2*p1k1*p2k3 + 6*p1k2
     &*p2k3 + 4*p1p2*p2k3 - 2*p2k2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*k1k3*p1k2*(-k1k2 + p1k1 + p1k2)))

!************************************************   2 Re(Mtt4k3 Mbtk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mtt4k3 Msk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*Nc**2*(-1 + Nc**2)*(-6*k2k3*M2 + k2k3*p1k1 
     &+ 6*k2k3*p1k2 + 5*p1k1*p1k2 - 5*p1k2**2 + k2k3*p1k3 - 6*M2*p1k3 + 4*p1k1*p1k3 - 6*p1k2*p1k3 - p1k3
     &**2 + 6*k2k3*p1p2 + 6*p1k3*p1p2 - k2k3*p2k1 - 5*p1k2*p2k1 - 10*p1k3*p2k1 - 6*k2k3*p2k2 - 5
     &*p1k1*p2k2 + 10*p1k2*p2k2 + 5*p2k1*p2k2 - 5*p2k2**2 + k1k3*(-2*M2 - p1k1 + 4*p1k2 + p1k3 + 2
     &*p1p2 + p2k1 - 4*p2k2 - p2k3) + k1k2*(2*M2 - 5*p1k1 + 3*p1k2 + p1k3 - 2*p1p2 + 5*p2k1 - 3*p2k2 
     &- p2k3) - k2k3*p2k3 - 6*M2*p2k3 + 2*p1k1*p2k3 + 12*p1k2*p2k3 + 2*p1k3*p2k3 + 6*p1p2*p2k3 + 4
     &*p2k1*p2k3 - 6*p2k2*p2k3 - p2k3**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*k2k3*(k1k2 - p1k1 - p1k2)*(M2 + p1p2)))

!************************************************   2 Re(Mtt4k3 Msk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*Nc**2*(-1 + Nc**2)*(-2*k2k3*M2 + 4*k2k3*p1k1 
     &- 5*p1k1**2 - k2k3*p1k2 + 5*p1k1*p1k2 + k2k3*p1k3 - 6*M2*p1k3 - 6*p1k1*p1k3 + 4*p1k2*p1k3 - p1k3
     &**2 + 2*k2k3*p1p2 + 6*p1k3*p1p2 - 4*k2k3*p2k1 + 10*p1k1*p2k1 - 5*p1k2*p2k1 - 5*p2k1**2 + k2k3
     &*p2k2 - 5*p1k1*p2k2 - 10*p1k3*p2k2 + 5*p2k1*p2k2 + k1k3*(-6*M2 + 6*p1k1 + p1k2 + p1k3 + 6
     &*p1p2 - 6*p2k1 - p2k2 - p2k3) + k1k2*(2*M2 + 3*p1k1 - 5*p1k2 + p1k3 - 2*p1p2 - 3*p2k1 + 5*p2k2 
     &- p2k3) - k2k3*p2k3 - 6*M2*p2k3 + 12*p1k1*p2k3 + 2*p1k2*p2k3 + 2*p1k3*p2k3 + 6*p1p2*p2k3 - 6
     &*p2k1*p2k3 + 4*p2k2*p2k3 - p2k3**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*k1k3*(k1k2 - p1k1 - p1k2)*(M2 + p1p2)))

!************************************************   2 Re(Mtt4k3 Msk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mtt4k3 Mvk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGG*(4 - 5*Nc**2 + Nc**4)*(2*k1k2 + k1k3 + k2k3 - 2
     &*p1k2 - p1k3 + 2*p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(32*k2k3*(k1k2 - p1k1 - p1k2)))

!************************************************   2 Re(Mtt4k3 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGG*(4 - 5*Nc**2 + Nc**4)*(2*k1k2 + k1k3 + k2k3 - 2
     &*p1k1 - p1k3 + 2*p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(32*k1k3*(k1k2 - p1k1 - p1k2)))

!************************************************   2 Re(Mtt4k3 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mtt4k3 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2)*(-1 + Nc**2)**2*(k1k2 + 4*M2 - 2*(3*p1k1 
     &+ p1k2))*(k1k2 + M2 - p1k1 - p1k2 - p1p2 + p2k1 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*Nc**2*p1k1*(-k1k2 + p1k1 + p1k2)**2))

!************************************************   2 Re(Mtt4k3 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((2 - 3*Nc**2 + Nc**4)*(k1k3 + 2*M2 - 3*p1k1 - 2*p1k3 
     &- 2*p1p2 + p2k1)*(2*k1k2 + k1k3 + k2k3 + 2*M2 - 4*p1k1 - 2*p1k2 - p1k3 - 2*p1p2 + 2*p2k1 + p2k3)
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*Nc**2*p1k1*(-k1k2 + p1k1 + p1k2)*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mtt4k3 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-2 + Nc**2 + Nc**4)*(k1k3 + k2k3 + 2*M2 - 2*p1k1 - 2
     &*p1k2 - p1k3 - 2*p1p2 + p2k3)*(3*k2k3 + 2*M2 - 3*p1k2 - 4*p1k3 - 2*p1p2 + p2k2 + 2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*Nc**2*(k1k2 
     &- p1k1 - p1k2)*p1k3*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mtt4k3 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((2 - 3*Nc**2 + Nc**4)*(k2k3 + 2*M2 - 3*p1k2 - 2*p1k3 
     &- 2*p1p2 + p2k2)*(2*k1k2 + k1k3 + k2k3 + 2*M2 - 2*p1k1 - 4*p1k2 - p1k3 - 2*p1p2 + 2*p2k2 + p2k3)
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*Nc**2*p1k2*(-k1k2 + p1k1 + p1k2)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mtt4k3 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2)*(-1 + Nc**2)**2*(k1k2 + 4*M2 - 2*(p1k1 
     &+ 3*p1k2))*(k1k2 + M2 - p1k1 - p1k2 - p1p2 + p2k1 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*Nc**2*p1k2*(-k1k2 + p1k1 + p1k2)**2))

!************************************************   2 Re(Mtt4k3 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-2 + Nc**2 + Nc**4)*(k1k3 + k2k3 + 2*M2 - 2*p1k1 - 2
     &*p1k2 - p1k3 - 2*p1p2 + p2k3)*(3*k1k3 + 2*M2 - 3*p1k1 - 4*p1k3 - 2*p1p2 + p2k1 + 2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*Nc**2*(k1k2 
     &- p1k1 - p1k2)*p1k3*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mbt4k1 Mbt4k2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-1 + Nc**2)*(k1k2 + 4*M2 - 2*(p1k1 + p1k2))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(8*Nc**2*p1k1*p1k2))

!************************************************   2 Re(Mbt4k1 Mbt4k3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-1 + Nc**2)*(k1k3 + 4*M2 - 2*(p1k1 + p1k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(8*Nc**2*p1k1*p1k3))

!************************************************   2 Re(Mbt4k1 Mttk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mbt4k1 Mttk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(2 - 3*Nc**2 + Nc**4)*(k1k3**2 + 4*p1k1**2 + 6
     &*p1k3**2 + 4*p1k3*p1p2 - 8*M2*p2k1 + 2*p1k3*p2k1 - 4*M2*p2k3 - 2*p1k3*p2k3 + k1k3*(-8*M2 + 2*p1k1 
     &- 5*p1k3 - 4*p1p2 + p2k3) + 2*p1k1*(4*p1k3 + 4*p1p2 + 2*p2k1 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k3*p1k1*(k1k3 - p1k1 - p1k3)))

!************************************************   2 Re(Mbt4k1 Mttk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(2 - 3*Nc**2 + Nc**4)*(k1k2**2 + 4*p1k1**2 + 6
     &*p1k2**2 + 4*p1k2*p1p2 - 8*M2*p2k1 + 2*p1k2*p2k1 - 4*M2*p2k2 - 2*p1k2*p2k2 + k1k2*(-8*M2 + 2*p1k1 
     &- 5*p1k2 - 4*p1p2 + p2k2) + 2*p1k1*(4*p1k2 + 4*p1p2 + 2*p2k1 + p2k2))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k2*p1k1*(k1k2 - p1k1 - p1k2)))

!************************************************   2 Re(Mbt4k1 Mbtk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mbt4k1 Mbtk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2**2 + 6*k2k3*M2 - 3*k2k3*p1k1 
     &+ 2*p1k1**2 - 4*k2k3*p1k2 + 4*p1k1*p1k3 - 6*p1k2*p1k3 - 2*k2k3*p1p2 + 4*p1k1*p1p2 - 4*p1k3*p1p2 
     &+ k2k3*p2k1 - 4*M2*p2k1 - 2*p1k1*p2k1 + 2*p1k2*p2k1 - 2*p1k3*p2k1 + 2*p1k1*p2k2 + 4*p1k3*p2k2 
     &- 2*k1k3*(2*M2 - 3*p1k2 - 2*p1p2 + p2k2) + 4*M2*p2k3 - 2*p1k1*p2k3 - 2*p1k2*p2k3 + k1k2*(2
     &*k2k3 - 3*p1k1 - 2*p1k2 - 3*p1k3 - 4*p1p2 + p2k1 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k1k3*p1k1*p1k2))

!************************************************   2 Re(Mbt4k1 Mbtk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k3**2 + 6*k2k3*M2 - 3*k2k3*p1k1 
     &+ 2*p1k1**2 + 4*p1k1*p1k2 - 4*k2k3*p1k3 - 6*p1k2*p1k3 - 2*k2k3*p1p2 + 4*p1k1*p1p2 - 4*p1k2*p1p2 
     &+ k2k3*p2k1 - 4*M2*p2k1 - 2*p1k1*p2k1 - 2*p1k2*p2k1 + 2*p1k3*p2k1 + 4*M2*p2k2 - 2*p1k1*p2k2 
     &- 2*p1k3*p2k2 + k1k3*(2*k2k3 - 3*p1k1 - 3*p1k2 - 2*p1k3 - 4*p1p2 + p2k1 + p2k2) + 2*p1k1*p2k3 
     &+ 4*p1k2*p2k3 - 2*k1k2*(2*M2 - 3*p1k3 - 2*p1p2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k1k2*p1k1*p1k3))

!************************************************   2 Re(Mbt4k1 Msk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mbt4k1 Msk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(2*k2k3*M2 - k2k3*p1k1 - 2
     &*p1k1**2 + 8*p1k1*p1k2 - 12*p1k1*p1k3 + 10*p1k2*p1k3 - 10*p1k3**2 - 2*k2k3*p1p2 + 12*p1k1*p1p2 
     &+ k2k3*p2k1 - 12*M2*p2k1 + 2*p1k1*p2k1 - 4*p1k2*p2k1 + 6*p1k3*p2k1 - 4*p1k1*p2k2 - 8*p1k3*p2k2 
     &+ k1k3*(-6*M2 + 6*p1k1 - 4*p1k2 + 5*p1k3 + 6*p1p2 - 6*p2k1 + 4*p2k2 - 5*p2k3) + 6*p1k1*p2k3 
     &- 2*p1k2*p2k3 + 10*p1k3*p2k3 + k1k2*(-2*M2 - p1k1 - p1k3 + 2*p1p2 + p2k1 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k3*p1k1*(M2 + p1p2)))

!************************************************   2 Re(Mbt4k1 Msk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(2*k2k3*M2 - k2k3*p1k1 - 2
     &*p1k1**2 - 12*p1k1*p1k2 - 10*p1k2**2 + 8*p1k1*p1k3 + 10*p1k2*p1k3 - 2*k2k3*p1p2 + 12*p1k1*p1p2 
     &+ k2k3*p2k1 - 12*M2*p2k1 + 2*p1k1*p2k1 + 6*p1k2*p2k1 - 4*p1k3*p2k1 + 6*p1k1*p2k2 + 10*p1k2
     &*p2k2 - 2*p1k3*p2k2 + k1k3*(-2*M2 - p1k1 - p1k2 + 2*p1p2 + p2k1 + p2k2) - 4*p1k1*p2k3 - 8*p1k2
     &*p2k3 + k1k2*(-6*M2 + 6*p1k1 + 5*p1k2 - 4*p1k3 + 6*p1p2 - 6*p2k1 - 5*p2k2 + 4*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k2*p1k1*(M2 + p1p2)))

!************************************************   2 Re(Mbt4k1 Mvk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mbt4k1 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGG*(4 - 5*Nc**2 + Nc**4)*(k1k3 - p1k1 - 2*p1k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(16*k1k3*p1k1))

!************************************************   2 Re(Mbt4k1 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGG*(4 - 5*Nc**2 + Nc**4)*(k1k2 - p1k1 - 2*p1k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(16*k1k2*p1k1))

!************************************************   2 Re(Mbt4k1 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2)*(-1 + Nc**2)**2*(M2 - p1k1)*(3*k1k2 + 2
     &*M2 - 4*p1k1 - 3*p1k2 - 2*p1p2 + 2*p2k1 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(16.
     &*Nc**2*p1k1**2*(-k1k2 + p1k1 + p1k2)))

!************************************************   2 Re(Mbt4k1 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2)*(-1 + Nc**2)**2*(M2 - p1k1)*(3*k1k3 + 2
     &*M2 - 4*p1k1 - 3*p1k3 - 2*p1p2 + 2*p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(16.
     &*Nc**2*p1k1**2*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mbt4k1 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-2 + Nc**2 + Nc**4)*(k2k3 + 4*M2 - 2*(p1k2 + 3*p1k3))
     &*(k1k2 + k1k3 + 2*M2 - p1k1 - 2*p1k2 - 2*p1k3 - 2*p1p2 + p2k1)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*Nc**2*p1k1*(k2k3 - p1k2 - p1k3)*p1k3))

!************************************************   2 Re(Mbt4k1 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2 + Nc**4)*(k2k3 + 4*M2 - 2*(3*p1k2 + p1k3))
     &*(k1k2 + k1k3 + 2*M2 - p1k1 - 2*p1k2 - 2*p1k3 - 2*p1p2 + p2k1)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*Nc**2*p1k1*p1k2*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mbt4k1 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((2 - 3*Nc**2 + Nc**4)*(k1k2 + 2*M2 - 2*p1k1 - 2*p1k2)
     &*(k1k2 + 2*M2 - 2*p1k1 - 3*p1k2 - 2*p1p2 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32.
     &*Nc**2*p1k1*p1k2*(-k1k2 + p1k1 + p1k2)))

!************************************************   2 Re(Mbt4k1 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((2 - 3*Nc**2 + Nc**4)*(k1k3 + 2*M2 - 2*p1k1 - 2*p1k3)
     &*(k1k3 + 2*M2 - 2*p1k1 - 3*p1k3 - 2*p1p2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32.
     &*Nc**2*p1k1*p1k3*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mbt4k2 Mbt4k3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-1 + Nc**2)*(k2k3 + 4*M2 - 2*(p1k2 + p1k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(8*Nc**2*p1k2*p1k3))

!************************************************   2 Re(Mbt4k2 Mttk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(2 - 3*Nc**2 + Nc**4)*(k2k3**2 + 4*p1k2**2 + 6
     &*p1k3**2 + 4*p1k3*p1p2 - 8*M2*p2k2 + 2*p1k3*p2k2 - 4*M2*p2k3 - 2*p1k3*p2k3 + k2k3*(-8*M2 + 2*p1k2 
     &- 5*p1k3 - 4*p1p2 + p2k3) + 2*p1k2*(4*p1k3 + 4*p1p2 + 2*p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k2k3*p1k2*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(Mbt4k2 Mttk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mbt4k2 Mttk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(2 - 3*Nc**2 + Nc**4)*(k1k2**2 + 6*p1k1**2 + 4
     &*p1k2**2 + 8*p1k2*p1p2 - 4*M2*p2k1 + 2*p1k2*p2k1 + k1k2*(-8*M2 - 5*p1k1 + 2*p1k2 - 4*p1p2 + p2k1) 
     &- 8*M2*p2k2 + 4*p1k2*p2k2 + 2*p1k1*(4*p1k2 + 2*p1p2 - p2k1 + p2k2))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k2*(k1k2 - p1k1 - p1k2)*p1k2))

!************************************************   2 Re(Mbt4k2 Mbtk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2**2 + k1k3*(6*M2 - 4*p1k1 - 3
     &*p1k2 - 2*p1p2 + p2k2) + k1k2*(2*k1k3 - 2*p1k1 - 3*p1k2 - 3*p1k3 - 4*p1p2 + p2k2 + p2k3) - 2*(
     &-p1k2**2 + 3*p1k1*p1k3 + 2*p1k3*p1p2 - 2*p1k3*p2k1 + k2k3*(2*M2 - 3*p1k1 - 2*p1p2 + p2k1) + 2
     &*M2*p2k2 - p1k1*p2k2 + p1k3*p2k2 - 2*M2*p2k3 + p1k1*p2k3 + p1k2*(-2*p1k3 - 2*p1p2 - p2k1 + p2k2 
     &+ p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k2k3*p1k1*p1k2))

!************************************************   2 Re(Mbt4k2 Mbtk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mbt4k2 Mbtk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k2k3**2 + k1k3*(2*k2k3 + 6*M2 - 3
     &*p1k2 - 4*p1k3 - 2*p1p2 + p2k2) + k2k3*(-3*p1k1 - 3*p1k2 - 2*p1k3 - 4*p1p2 + p2k1 + p2k2) - 2
     &*(-p1k2**2 - 2*p1k2*p1p2 - 2*M2*p2k1 + p1k2*p2k1 + p1k3*p2k1 + 2*M2*p2k2 + p1k2*p2k2 - p1k3*p2k2 
     &+ p1k1*(-2*p1k2 + 3*p1k3 + 2*p1p2 + p2k2 - 2*p2k3) - p1k2*p2k3 + k1k2*(2*M2 - 3*p1k3 - 2*p1p2 
     &+ p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k1k2*p1k2*p1k3))

!************************************************   2 Re(Mbt4k2 Msk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(-6*k2k3*M2 - 4*k2k3*p1k1 
     &+ 6*k2k3*p1k2 + 8*p1k1*p1k2 - 2*p1k2**2 + 5*k2k3*p1k3 + 10*p1k1*p1k3 - 12*p1k2*p1k3 - 10*p1k3
     &**2 + 6*k2k3*p1p2 + 12*p1k2*p1p2 + 4*k2k3*p2k1 - 4*p1k2*p2k1 - 8*p1k3*p2k1 - 6*k2k3*p2k2 - 12
     &*M2*p2k2 - 4*p1k1*p2k2 + 2*p1k2*p2k2 + 6*p1k3*p2k2 + k1k3*(2*M2 - p1k2 - 2*p1p2 + p2k2) - 5*k2k3
     &*p2k3 - 2*p1k1*p2k3 + 6*p1k2*p2k3 + 10*p1k3*p2k3 + k1k2*(-2*M2 - p1k2 - p1k3 + 2*p1p2 + p2k2 
     &+ p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k2k3*p1k2*(M2 + p1p2)))

!************************************************   2 Re(Mbt4k2 Msk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mbt4k2 Msk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(-2*k2k3*M2 - k2k3*p1k1 - 10
     &*p1k1**2 - k2k3*p1k2 - 12*p1k1*p1k2 - 2*p1k2**2 + 10*p1k1*p1k3 + 8*p1k2*p1k3 + 2*k2k3*p1p2 
     &+ 12*p1k2*p1p2 + k2k3*p2k1 + 10*p1k1*p2k1 + 6*p1k2*p2k1 - 2*p1k3*p2k1 + k2k3*p2k2 - 12*M2*p2k2 
     &+ 6*p1k1*p2k2 + 2*p1k2*p2k2 - 4*p1k3*p2k2 + k1k3*(2*M2 - p1k2 - 2*p1p2 + p2k2) - 8*p1k1*p2k3 
     &- 4*p1k2*p2k3 + k1k2*(-6*M2 + 5*p1k1 + 6*p1k2 - 4*p1k3 + 6*p1p2 - 5*p2k1 - 6*p2k2 + 4*p2k3))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k2*p1k2*(M2 + p1p2)))

!************************************************   2 Re(Mbt4k2 Mvk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGG*(4 - 5*Nc**2 + Nc**4)*(k2k3 - p1k2 - 2*p1k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(16*k2k3*p1k2))

!************************************************   2 Re(Mbt4k2 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mbt4k2 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGG*(4 - 5*Nc**2 + Nc**4)*(k1k2 - 2*p1k1 - p1k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(16*k1k2*p1k2))

!************************************************   2 Re(Mbt4k2 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((2 - 3*Nc**2 + Nc**4)*(k1k2 + 2*M2 - 2*p1k1 - 2*p1k2)
     &*(k1k2 + 2*M2 - 3*p1k1 - 2*p1k2 - 2*p1p2 + p2k1)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32.
     &*Nc**2*p1k1*p1k2*(-k1k2 + p1k1 + p1k2)))

!************************************************   2 Re(Mbt4k2 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2 + Nc**4)*(k1k3 + 4*M2 - 2*(3*p1k1 + p1k3))
     &*(k1k2 + k2k3 + 2*M2 - 2*p1k1 - p1k2 - 2*p1k3 - 2*p1p2 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*Nc**2*p1k1*p1k2*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mbt4k2 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((2 - 3*Nc**2 + Nc**4)*(k2k3 + 2*M2 - 2*p1k2 - 2*p1k3)
     &*(k2k3 + 2*M2 - 2*p1k2 - 3*p1k3 - 2*p1p2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32.
     &*Nc**2*p1k2*p1k3*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mbt4k2 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2)*(-1 + Nc**2)**2*(M2 - p1k2)*(3*k2k3 + 2
     &*M2 - 4*p1k2 - 3*p1k3 - 2*p1p2 + 2*p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(16.
     &*Nc**2*p1k2**2*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mbt4k2 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2)*(-1 + Nc**2)**2*(M2 - p1k2)*(3*k1k2 + 2
     &*M2 - 3*p1k1 - 4*p1k2 - 2*p1p2 + p2k1 + 2*p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(16.
     &*Nc**2*p1k2**2*(-k1k2 + p1k1 + p1k2)))

!************************************************   2 Re(Mbt4k2 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-2 + Nc**2 + Nc**4)*(k1k3 + 4*M2 - 2*(p1k1 + 3*p1k3))
     &*(k1k2 + k2k3 + 2*M2 - 2*p1k1 - p1k2 - 2*p1k3 - 2*p1p2 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*Nc**2*p1k2*(k1k3 - p1k1 - p1k3)*p1k3))

!************************************************   2 Re(Mbt4k3 Mttk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(2 - 3*Nc**2 + Nc**4)*(k2k3**2 + 6*p1k2**2 + 4
     &*p1k3**2 + 8*p1k3*p1p2 - 4*M2*p2k2 + 2*p1k3*p2k2 + k2k3*(-8*M2 - 5*p1k2 + 2*p1k3 - 4*p1p2 + p2k2) 
     &- 8*M2*p2k3 + 4*p1k3*p2k3 + 2*p1k2*(4*p1k3 + 2*p1p2 - p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k2k3*(k2k3 - p1k2 - p1k3)*p1k3))

!************************************************   2 Re(Mbt4k3 Mttk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(2 - 3*Nc**2 + Nc**4)*(k1k3**2 + 6*p1k1**2 + 4
     &*p1k3**2 + 8*p1k3*p1p2 - 4*M2*p2k1 + 2*p1k3*p2k1 + k1k3*(-8*M2 - 5*p1k1 + 2*p1k3 - 4*p1p2 + p2k1) 
     &- 8*M2*p2k3 + 4*p1k3*p2k3 + 2*p1k1*(4*p1k3 + 2*p1p2 - p2k1 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k3*(k1k3 - p1k1 - p1k3)*p1k3))

!************************************************   2 Re(Mbt4k3 Mttk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mbt4k3 Mbtk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k3**2 + k1k2*(2*k1k3 + 6*M2 - 4
     &*p1k1 - 3*p1k3 - 2*p1p2 + p2k3) + k1k3*(-2*p1k1 - 3*p1k2 - 3*p1k3 - 4*p1p2 + p2k2 + p2k3) - 2
     &*(-2*p1k2*p1k3 - p1k3**2 + 2*p1k2*p1p2 - 2*p1k3*p1p2 - 2*p1k2*p2k1 - p1k3*p2k1 + k2k3*(2*M2 - 3
     &*p1k1 - 2*p1p2 + p2k1) - 2*M2*p2k2 + p1k3*p2k2 + p1k1*(3*p1k2 + p2k2 - p2k3) + 2*M2*p2k3 + p1k2
     &*p2k3 + p1k3*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k2k3*p1k1*p1k3))

!************************************************   2 Re(Mbt4k3 Mbtk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k2k3**2 + k1k2*(2*k2k3 + 6*M2 - 4
     &*p1k2 - 3*p1k3 - 2*p1p2 + p2k3) + k2k3*(-3*p1k1 - 2*p1k2 - 3*p1k3 - 4*p1p2 + p2k1 + p2k3) - 2
     &*(-p1k3**2 - 2*p1k3*p1p2 - 2*M2*p2k1 + p1k2*p2k1 + p1k3*p2k1 - p1k3*p2k2 + k1k3*(2*M2 - 3*p1k2 
     &- 2*p1p2 + p2k2) + 2*M2*p2k3 - p1k2*p2k3 + p1k3*p2k3 + p1k1*(3*p1k2 - 2*p1k3 + 2*p1p2 - 2*p2k2 
     &+ p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k1k3*p1k2*p1k3))

!************************************************   2 Re(Mbt4k3 Mbtk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mbt4k3 Msk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(-6*k2k3*M2 - 4*k2k3*p1k1 
     &+ 5*k2k3*p1k2 + 10*p1k1*p1k2 - 10*p1k2**2 + 6*k2k3*p1k3 + 8*p1k1*p1k3 - 12*p1k2*p1k3 - 2*p1k3
     &**2 + 6*k2k3*p1p2 + 12*p1k3*p1p2 + 4*k2k3*p2k1 - 8*p1k2*p2k1 - 4*p1k3*p2k1 - 5*k2k3*p2k2 - 2
     &*p1k1*p2k2 + 10*p1k2*p2k2 + 6*p1k3*p2k2 - 6*k2k3*p2k3 - 12*M2*p2k3 - 4*p1k1*p2k3 + 6*p1k2*p2k3 
     &+ 2*p1k3*p2k3 + k1k2*(2*M2 - p1k3 - 2*p1p2 + p2k3) + k1k3*(-2*M2 - p1k2 - p1k3 + 2*p1p2 + p2k2 
     &+ p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k2k3*p1k3*(M2 + p1p2)))

!************************************************   2 Re(Mbt4k3 Msk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(-2*k2k3*M2 - k2k3*p1k1 - 10
     &*p1k1**2 + 10*p1k1*p1k2 - k2k3*p1k3 - 12*p1k1*p1k3 + 8*p1k2*p1k3 - 2*p1k3**2 + 2*k2k3*p1p2 
     &+ 12*p1k3*p1p2 + k2k3*p2k1 + 10*p1k1*p2k1 - 2*p1k2*p2k1 + 6*p1k3*p2k1 - 8*p1k1*p2k2 - 4*p1k3
     &*p2k2 + k1k3*(-6*M2 + 5*p1k1 - 4*p1k2 + 6*p1k3 + 6*p1p2 - 5*p2k1 + 4*p2k2 - 6*p2k3) + k2k3*p2k3 
     &- 12*M2*p2k3 + 6*p1k1*p2k3 - 4*p1k2*p2k3 + 2*p1k3*p2k3 + k1k2*(2*M2 - p1k3 - 2*p1p2 + p2k3))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k3*p1k3*(M2 + p1p2)))

!************************************************   2 Re(Mbt4k3 Msk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mbt4k3 Mvk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGG*(4 - 5*Nc**2 + Nc**4)*(k2k3 - 2*p1k2 - p1k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(16*k2k3*p1k3))

!************************************************   2 Re(Mbt4k3 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((3*GGG*(4 - 5*Nc**2 + Nc**4)*(k1k3 - 2*p1k1 - p1k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(16*k1k3*p1k3))

!************************************************   2 Re(Mbt4k3 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mbt4k3 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2 + Nc**4)*(k1k2 + 4*M2 - 2*(3*p1k1 + p1k2))
     &*(k1k3 + k2k3 + 2*M2 - 2*p1k1 - 2*p1k2 - p1k3 - 2*p1p2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*Nc**2*p1k1*(-k1k2 + p1k1 + p1k2)*p1k3))

!************************************************   2 Re(Mbt4k3 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((2 - 3*Nc**2 + Nc**4)*(k1k3 + 2*M2 - 2*p1k1 - 2*p1k3)
     &*(k1k3 + 2*M2 - 3*p1k1 - 2*p1k3 - 2*p1p2 + p2k1)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32.
     &*Nc**2*p1k1*p1k3*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mbt4k3 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2)*(-1 + Nc**2)**2*(M2 - p1k3)*(3*k2k3 + 2
     &*M2 - 3*p1k2 - 4*p1k3 - 2*p1p2 + p2k2 + 2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(16.
     &*Nc**2*p1k3**2*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mbt4k3 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((2 - 3*Nc**2 + Nc**4)*(k2k3 + 2*M2 - 2*p1k2 - 2*p1k3)
     &*(k2k3 + 2*M2 - 3*p1k2 - 2*p1k3 - 2*p1p2 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32.
     &*Nc**2*p1k2*p1k3*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mbt4k3 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-2 + Nc**2 + Nc**4)*(k1k2 + 4*M2 - 2*(p1k1 + 3*p1k2))
     &*(k1k3 + k2k3 + 2*M2 - 2*p1k1 - 2*p1k2 - p1k3 - 2*p1p2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*Nc**2*(k1k2 - p1k1 - p1k2)*p1k2*p1k3))

!************************************************   2 Re(Mbt4k3 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-2 + Nc**2)*(-1 + Nc**2)**2*(M2 - p1k3)*(3*k1k3 + 2
     &*M2 - 3*p1k1 - 4*p1k3 - 2*p1p2 + p2k1 + 2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(16.
     &*Nc**2*p1k3**2*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mttk1 Mttk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mttk1 Mttk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mttk1 Mbtk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*(-1 + Nc**2)*(k1k2 + k1k3 + 2*M2 - p1k1 - 2
     &*p1k2 - 2*p1k3 - 2*p1p2 + p2k1)*(16*k2k3*M2 - 16*k2k3*p1k1 + k2k3*p1k2 - 4*p1k2**2 + k2k3*p1k3 
     &- 28*p1k2*p1k3 - 4*p1k3**2 + k1k3*(-k2k3 + 14*p1k2 + 4*p1k3) + k1k2*(-k2k3 + 4*p1k2 + 14*p1k3) 
     &- 16*k2k3*p1p2 - k2k3*p2k2 + 4*p1k2*p2k2 + 14*p1k3*p2k2 - k2k3*p2k3 + 14*p1k2*p2k3 + 4*p1k3
     &*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(64*k2k3**2*p1k1*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(Mttk1 Mbtk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*Nc**2*(-1 + Nc**2)*(k2k3**3 + 16*k2k3**2*M2 
     &+ 16*k2k3*M2**2 + 6*k2k3**2*p1k1 + 24*k2k3*M2*p1k1 - 10*k2k3**2*p1k2 - 24*k2k3*M2*p1k2 - 32*k2k3
     &*p1k1*p1k2 - 32*M2*p1k1*p1k2 + 12*k2k3*p1k2**2 + 24*p1k1*p1k2**2 - 12*k2k3*M2*p1k3 - 8*k2k3
     &*p1k1*p1k3 - 16*M2*p1k1*p1k3 - 10*k2k3*p1k2*p1k3 - 16*M2*p1k2*p1k3 + 24*p1k1*p1k2*p1k3 + 12
     &*p1k2**2*p1k3 - 7*k2k3*p1k3**2 - 8*M2*p1k3**2 + 24*p1k2*p1k3**2 + 6*p1k3**3 - 4*k2k3**2*p1p2 
     &- 16*k2k3*M2*p1p2 - 8*k2k3*p1k1*p1p2 + 16*k2k3*p1k2*p1p2 + 16*p1k1*p1k2*p1p2 + 4*k2k3*p1k3*p1p2 
     &+ 8*p1k1*p1k3*p1p2 + 8*p1k2*p1k3*p1p2 + 4*p1k3**2*p1p2 + 2*k2k3**2*p2k1 - 8*k2k3*M2*p2k1 + 8
     &*p1k2**2*p2k1 - 8*k2k3*p1k3*p2k1 + 8*p1k2*p1k3*p2k1 + 8*p1k3**2*p2k1 - 8*k2k3*p1p2*p2k1 + 16
     &*p1k2*p1p2*p2k1 + 8*p1k3*p1p2*p2k1 + 8*k2k3*M2*p2k2 + 16*k2k3*p1k1*p2k2 + 16*M2*p1k1*p2k2 - 4
     &*k2k3*p1k2*p2k2 - 40*p1k1*p1k2*p2k2 + 6*k2k3*p1k3*p2k2 + 8*M2*p1k3*p2k2 - 20*p1k1*p1k3*p2k2 
     &- 20*p1k2*p1k3*p2k2 - 10*p1k3**2*p2k2 - 16*M2*p2k1*p2k2 + 8*p1k2*p2k1*p2k2 + 4*p1k3*p2k1*p2k2 
     &+ 2*k2k3**2*p2k3 + 12*k2k3*M2*p2k3 + 6*k2k3*p1k1*p2k3 + 8*M2*p1k1*p2k3 - 8*k2k3*p1k2*p2k3 - 20
     &*p1k1*p1k2*p2k3 + 4*p1k2**2*p2k3 - 2*k2k3*p1k3*p2k3 + 4*M2*p1k3*p2k3 - 4*p1k1*p1k3*p2k3 - 18
     &*p1k2*p1k3*p2k3 - 4*p1k3**2*p2k3 - 4*k2k3*p1p2*p2k3 + 8*p1k2*p1p2*p2k3 + 4*p1k3*p1p2*p2k3 + 2
     &*k2k3*p2k1*p2k3 - 8*M2*p2k1*p2k3 + 4*p1k2*p2k1*p2k3 - 4*p1k3*p2k1*p2k3 - 8*M2*p2k2*p2k3 + 4*p1k2
     &*p2k2*p2k3 + 2*p1k3*p2k2*p2k3 + k2k3*p2k3**2 - 4*M2*p2k3**2 + 2*p1k2*p2k3**2 - 2*p1k3*p2k3
     &**2 + 2*k1k2*(k2k3**2 + 2*(8*M2**2 + 6*p1k2**2 + p1k3*(3*p1k3 + 4*p1p2 - p2k2 - 2*p2k3) - 2*M2
     &*(6*p1k2 + 3*p1k3 + 4*p1p2 - 2*p2k2 - p2k3) + p1k2*(6*p1k3 + 8*p1p2 - 2*p2k2 - p2k3)) + k2k3
     &*(4*M2 - 4*p1k2 - 6*p1k3 - 4*p1p2 + p2k3)) - k1k3*(k2k3**2 + k2k3*(-4*M2 + 2*p1k2 - 3*p1k3 + 4
     &*p2k2 + p2k3) - 2*(8*M2**2 + 6*p1k2**2 - p1k2*(2*p1k3 - 8*p1p2 + 2*p2k2 + p2k3) + 2*M2*(-6*p1k2 
     &- 4*p1p2 + 2*p2k2 + p2k3) + p1k3*(-p1k3 - 2*p1p2 + 5*p2k2 + p2k3))))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)
     &/(128*k1k3*k2k3*p1k2*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(Mttk1 Mbtk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(-4*(k2k3 - 2*p1k2 - p1k3)
     &*(-2*p1k1*(-M2 + p1k2 + p1k3 + p1p2) + (-M2 + p1k2 + p1k3 + p1p2)*(-k1k3 + p1k1 - p2k1) + (M2 
     &- p1k3 - p1p2)*(k1k2 + k1k3 - p1k1 + p2k1) + 2*M2*(k2k3 - p1k2 + p2k2) - 2*(p1k3 + p1p2)*(k2k3 
     &- p1k2 + p2k2) + (-M2 + p1k2 + p1k3 + p1p2)*(k2k3 - p1k2 + p2k2) + p1k2*(k2k3 + 2*M2 - p1k2 
     &- 2*p1k3 - 2*p1p2 + p2k2 + 2*p2k3) + 2*p1k1*(k2k3 + M2 - p1k3 - p1p2 + p2k2 + 2*p2k3)) - (k2k3 
     &- 4*M2 - 2*p1k2 + 2*p1k3)*((k2k3 - p1k2 + p2k2)*(2*k1k3 + k2k3 - 2*p1k1 - p1k2 + 2*p2k1 
     &+ p2k2) - k1k2*(k2k3 + 4*M2 - p1k2 - 4*p1k3 - 4*p1p2 + p2k2 + 4*p2k3)) + 2*(p1k2 - p1k3)*(2*k2k3
     &**2 - 2*k1k2*M2 - 4*M2*p1k1 - 2*M2*p1k2 + 2*p1k1*p1k2 + 4*p1k2**2 + 3*k1k2*p1k3 + 2*p1k1*p1k3 
     &+ p1k2*p1k3 + 2*k1k2*p1p2 + 4*p1k1*p1p2 + 2*p1k2*p1p2 - 4*M2*p2k1 + 2*p1k2*p2k1 + 6*p1k3*p2k1 
     &+ 4*p1p2*p2k1 - 2*M2*p2k2 - 4*p1k1*p2k2 - 4*p1k2*p2k2 + 3*p1k3*p2k2 + 2*p1p2*p2k2 - k1k2
     &*p2k3 - 6*p1k1*p2k3 - 3*p1k2*p2k3 - 2*p2k1*p2k3 - p2k2*p2k3 + k2k3*(-3*p1k1 - 6*p1k2 + p1k3 - p2k1 
     &+ 2*p2k2 + p2k3) + k1k3*(k2k3 + 2*(p1k3 + p2k2 + p2k3))) - 2*(4*M2 - 3*p1k3)*(2*k2k3**2 
     &+ (k1k2 - 2*p1k1 - p1k2 + 2*p2k1 + p2k2)*(p1k3 - p2k3) + k2k3*(2*M2 + p1k1 - 2*p1k2 - p1k3 - 2
     &*p1p2 - p2k1 + 2*p2k2 + p2k3) + k1k3*(k2k3 + 2*(2*M2 - p1k2 - p1k3 - 2*p1p2 + p2k2 + p2k3))))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(128*k1k2*k2k3*(k2k3 - p1k2 - p1k3)*p1k3))

!************************************************   2 Re(Mttk1 Msk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**3*Nc**2*(-1 + Nc**2)*(32*k2k3**2*M2 + 32*k2k3*M2
     &**2 + 2*k2k3**2*p1k1 + 16*k2k3*M2*p1k1 + 2*k2k3*M2*p1k2 - 19*k2k3*p1k1*p1k2 - 6*k2k3*p1k2**2 
     &- 8*M2*p1k2**2 + 4*p1k1*p1k2**2 + 2*k2k3*M2*p1k3 - 19*k2k3*p1k1*p1k3 - 52*k2k3*p1k2*p1k3 - 56
     &*M2*p1k2*p1k3 + 28*p1k1*p1k2*p1k3 - 6*k2k3*p1k3**2 - 8*M2*p1k3**2 + 4*p1k1*p1k3**2 - 32*k2k3
     &**2*p1p2 - 16*k2k3*p1k1*p1p2 - 34*k2k3*p1k2*p1p2 + 8*p1k2**2*p1p2 - 34*k2k3*p1k3*p1p2 + 56*p1k2
     &*p1k3*p1p2 + 8*p1k3**2*p1p2 - 32*k2k3*p1p2**2 - 2*k2k3**2*p2k1 + 16*k2k3*M2*p2k1 + 3*k2k3*p1k2
     &*p2k1 - 4*p1k2**2*p2k1 + 3*k2k3*p1k3*p2k1 - 28*p1k2*p1k3*p2k1 - 4*p1k3**2*p2k1 - 16*k2k3*p1p2
     &*p2k1 + 34*k2k3*M2*p2k2 + 17*k2k3*p1k1*p2k2 + 6*k2k3*p1k2*p2k2 - 8*M2*p1k2*p2k2 - 4*p1k1*p1k2
     &*p2k2 + 26*k2k3*p1k3*p2k2 - 28*M2*p1k3*p2k2 - 14*p1k1*p1k3*p2k2 - 20*p1k2*p1k3*p2k2 + 20*p1k3
     &**2*p2k2 - 2*k2k3*p1p2*p2k2 + 8*p1k2*p1p2*p2k2 + 28*p1k3*p1p2*p2k2 - k2k3*p2k1*p2k2 + 4*p1k2
     &*p2k1*p2k2 + 14*p1k3*p2k1*p2k2 + 34*k2k3*M2*p2k3 + 17*k2k3*p1k1*p2k3 + 26*k2k3*p1k2*p2k3 - 28
     &*M2*p1k2*p2k3 - 14*p1k1*p1k2*p2k3 + 20*p1k2**2*p2k3 + 6*k2k3*p1k3*p2k3 - 8*M2*p1k3*p2k3 - 4
     &*p1k1*p1k3*p2k3 - 20*p1k2*p1k3*p2k3 - 2*k2k3*p1p2*p2k3 + 28*p1k2*p1p2*p2k3 + 8*p1k3*p1p2*p2k3 
     &- k2k3*p2k1*p2k3 + 14*p1k2*p2k1*p2k3 + 4*p1k3*p2k1*p2k3 + k1k2*(k2k3 - 4*p1k2 - 14*p1k3)*(2
     &*M2 - p1k2 - p1k3 - 2*p1p2 + p2k2 + p2k3) + k1k3*(k2k3 - 14*p1k2 - 4*p1k3)*(2*M2 - p1k2 - p1k3 
     &- 2*p1p2 + p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(64*k2k3**2*(k2k3 - p1k2 
     &- p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Mttk1 Msk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**3*Nc**2*(-1 + Nc**2)*(-8*k2k3*M2**2 + 3*k2k3**2
     &*p1k1 - 4*k2k3*M2*p1k1 + 12*k2k3*p1k1**2 + k2k3**2*p1k2 + 2*k2k3*M2*p1k2 - 24*k2k3*p1k1*p1k2 
     &+ 20*M2*p1k1*p1k2 + 6*k2k3*p1k2**2 + k2k3**2*p1k3 - 6*k2k3*M2*p1k3 + 19*k2k3*p1k1*p1k3 + 4*M2
     &*p1k1*p1k3 - 6*p1k1**2*p1k3 - 19*k2k3*p1k2*p1k3 + 4*M2*p1k2*p1k3 + 24*p1k1*p1k2*p1k3 - 6*p1k2
     &**2*p1k3 + 9*k2k3*p1k3**2 + 8*M2*p1k3**2 - 12*p1k1*p1k3**2 + 24*p1k2*p1k3**2 - 6*p1k3**3 - 8
     &*k2k3*p1k1*p1p2 - 18*k2k3*p1k2*p1p2 + 2*k2k3*p1k3*p1p2 + 12*p1k1*p1k3*p1p2 + 12*p1k2*p1k3*p1p2 
     &+ 12*p1k3**2*p1p2 + 8*k2k3*p1p2**2 - 3*k2k3**2*p2k1 + 16*k2k3*M2*p2k1 - 12*k2k3*p1k1*p2k1 + 8
     &*M2*p1k2*p2k1 - 10*p1k2**2*p2k1 + 11*k2k3*p1k3*p2k1 - 8*M2*p1k3*p2k1 + 6*p1k1*p1k3*p2k1 - 38
     &*p1k2*p1k3*p2k1 + 6*p1k3**2*p2k1 - 4*k2k3*p1p2*p2k1 - 28*p1k2*p1p2*p2k1 - 8*p1k3*p1p2*p2k1 - k2k3
     &**2*p2k2 + 14*k2k3*M2*p2k2 + 28*k2k3*p1k1*p2k2 - 8*M2*p1k1*p2k2 - 6*k2k3*p1k2*p2k2 + 10*p1k1
     &*p1k2*p2k2 + 5*k2k3*p1k3*p2k2 - 16*M2*p1k3*p2k2 - 14*p1k1*p1k3*p2k2 + 8*p1k2*p1k3*p2k2 - 10
     &*p1k3**2*p2k2 + 2*k2k3*p1p2*p2k2 - 12*p1k1*p1p2*p2k2 - 4*k2k3*p2k1*p2k2 + 20*M2*p2k1*p2k2 - 6
     &*p1k2*p2k1*p2k2 + 16*p1k3*p2k1*p2k2 + 6*p1k1*p2k2**2 - 6*p1k3*p2k2**2 + k1k3**2*(2*p1k2 + p1k3 
     &- 2*p2k2 - p2k3) - k2k3**2*p2k3 + 10*k2k3*M2*p2k3 - 29*k2k3*p1k1*p2k3 - 4*M2*p1k1*p2k3 + 6
     &*p1k1**2*p2k3 + 17*k2k3*p1k2*p2k3 + 4*M2*p1k2*p2k3 + 4*p1k1*p1k2*p2k3 + 4*p1k2**2*p2k3 - 8*k2k3
     &*p1k3*p2k3 - 12*M2*p1k3*p2k3 + 18*p1k1*p1k3*p2k3 - 38*p1k2*p1k3*p2k3 + 12*p1k3**2*p2k3 - 6
     &*k2k3*p1p2*p2k3 - 12*p1k1*p1p2*p2k3 - 20*p1k2*p1p2*p2k3 - 28*p1k3*p1p2*p2k3 - k2k3*p2k1*p2k3 
     &+ 16*M2*p2k1*p2k3 - 6*p1k1*p2k1*p2k3 + 10*p1k2*p2k1*p2k3 - 2*p1k3*p2k1*p2k3 - 3*k2k3*p2k2*p2k3 
     &+ 16*M2*p2k2*p2k3 - 2*p1k1*p2k2*p2k3 + 10*p1k3*p2k2*p2k3 - k2k3*p2k3**2 + 20*M2*p2k3**2 - 10
     &*p1k1*p2k3**2 + 14*p1k2*p2k3**2 - 6*p1k3*p2k3**2 + k1k2**2*(6*M2 - p1k2 - 2*p1k3 - 6*p1p2 + p2k2 
     &+ 2*p2k3) + k1k2*(-16*M2**2 - 2*M2*p1k1 + 4*M2*p1k2 - 11*p1k1*p1k2 + 6*p1k2**2 - 6*M2*p1k3 
     &- 6*p1k1*p1k3 - 18*p1k2*p1k3 - 7*p1k3**2 + 2*p1k1*p1p2 + 10*p1k3*p1p2 + 16*p1p2**2 - 6*M2*p2k1 
     &+ 11*p1k2*p2k1 + 8*p1k3*p2k1 + 6*p1p2*p2k1 + 4*M2*p2k2 + 3*p1k1*p2k2 - 6*p1k2*p2k2 - 12*p1k3
     &*p2k2 - 8*p1p2*p2k2 - 3*p2k1*p2k2 - 14*M2*p2k3 - 2*p1k1*p2k3 + 32*p1k2*p2k3 + 8*p1k3*p2k3 
     &+ 10*p1p2*p2k3 - 2*p2k2*p2k3 - p2k3**2 + k2k3*(12*M2 + 2*p1k1 - p1k3 - 12*p1p2 - 2*p2k1 + p2k3) 
     &+ k1k3*(6*M2 + p1k2 - p1k3 - 6*p1p2 - p2k2 + p2k3)) + k1k3*(-8*M2**2 - 12*p1k1*p1k2 + 2*p1k2
     &**2 - 7*p1k1*p1k3 - 17*p1k2*p1k3 - 2*p1k3**2 - 8*p1k1*p1p2 - 12*p1k2*p1p2 + 8*p1k3*p1p2 + 8
     &*p1p2**2 + 4*p1k2*p2k1 + p1k3*p2k1 + 8*p1k1*p2k2 - 2*p1k2*p2k2 - 13*p1k3*p2k2 - 4*p1p2*p2k2 
     &+ 4*M2*(2*p1k1 + 2*p1k2 - 3*p1k3 + 2*p2k2 - p2k3) + 3*p1k1*p2k3 + 29*p1k2*p2k3 + 8*p1p2*p2k3 
     &+ 3*p2k1*p2k3 + p2k2*p2k3 + 2*p2k3**2 + k2k3*(8*M2 - 2*p1k1 + 2*p1k2 - 3*p1k3 - 8*p1p2 + 2*p2k1 - 2*p2k2 
     &+ 3*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k3*k2k3*(k2k3 - p1k2 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Mttk1 Msk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**3*Nc**2*(-1 + Nc**2)*(-8*k2k3*M2**2 + 3*k2k3**2
     &*p1k1 - 4*k2k3*M2*p1k1 + 12*k2k3*p1k1**2 + k2k3**2*p1k2 - 6*k2k3*M2*p1k2 + 19*k2k3*p1k1*p1k2 
     &+ 4*M2*p1k1*p1k2 - 6*p1k1**2*p1k2 + 9*k2k3*p1k2**2 + 8*M2*p1k2**2 - 12*p1k1*p1k2**2 - 6*p1k2
     &**3 + k2k3**2*p1k3 + 2*k2k3*M2*p1k3 - 24*k2k3*p1k1*p1k3 + 20*M2*p1k1*p1k3 - 19*k2k3*p1k2*p1k3 
     &+ 4*M2*p1k2*p1k3 + 24*p1k1*p1k2*p1k3 + 24*p1k2**2*p1k3 + 6*k2k3*p1k3**2 - 6*p1k2*p1k3**2 - 8
     &*k2k3*p1k1*p1p2 + 2*k2k3*p1k2*p1p2 + 12*p1k1*p1k2*p1p2 + 12*p1k2**2*p1p2 - 18*k2k3*p1k3*p1p2 
     &+ 12*p1k2*p1k3*p1p2 + 8*k2k3*p1p2**2 - 3*k2k3**2*p2k1 + 16*k2k3*M2*p2k1 - 12*k2k3*p1k1*p2k1 + 11
     &*k2k3*p1k2*p2k1 - 8*M2*p1k2*p2k1 + 6*p1k1*p1k2*p2k1 + 6*p1k2**2*p2k1 + 8*M2*p1k3*p2k1 - 38
     &*p1k2*p1k3*p2k1 - 10*p1k3**2*p2k1 - 4*k2k3*p1p2*p2k1 - 8*p1k2*p1p2*p2k1 - 28*p1k3*p1p2*p2k1 - k2k3
     &**2*p2k2 + 10*k2k3*M2*p2k2 - 29*k2k3*p1k1*p2k2 - 4*M2*p1k1*p2k2 + 6*p1k1**2*p2k2 - 8*k2k3
     &*p1k2*p2k2 - 12*M2*p1k2*p2k2 + 18*p1k1*p1k2*p2k2 + 12*p1k2**2*p2k2 + 17*k2k3*p1k3*p2k2 + 4*M2
     &*p1k3*p2k2 + 4*p1k1*p1k3*p2k2 - 38*p1k2*p1k3*p2k2 + 4*p1k3**2*p2k2 - 6*k2k3*p1p2*p2k2 - 12*p1k1
     &*p1p2*p2k2 - 28*p1k2*p1p2*p2k2 - 20*p1k3*p1p2*p2k2 - k2k3*p2k1*p2k2 + 16*M2*p2k1*p2k2 - 6*p1k1
     &*p2k1*p2k2 - 2*p1k2*p2k1*p2k2 + 10*p1k3*p2k1*p2k2 - k2k3*p2k2**2 + 20*M2*p2k2**2 - 10*p1k1
     &*p2k2**2 - 6*p1k2*p2k2**2 + 14*p1k3*p2k2**2 + k1k2**2*(p1k2 + 2*p1k3 - p2k2 - 2*p2k3) - k2k3
     &**2*p2k3 + 14*k2k3*M2*p2k3 + 28*k2k3*p1k1*p2k3 - 8*M2*p1k1*p2k3 + 5*k2k3*p1k2*p2k3 - 16*M2*p1k2
     &*p2k3 - 14*p1k1*p1k2*p2k3 - 10*p1k2**2*p2k3 - 6*k2k3*p1k3*p2k3 + 10*p1k1*p1k3*p2k3 + 8*p1k2
     &*p1k3*p2k3 + 2*k2k3*p1p2*p2k3 - 12*p1k1*p1p2*p2k3 - 4*k2k3*p2k1*p2k3 + 20*M2*p2k1*p2k3 + 16*p1k2
     &*p2k1*p2k3 - 6*p1k3*p2k1*p2k3 - 3*k2k3*p2k2*p2k3 + 16*M2*p2k2*p2k3 - 2*p1k1*p2k2*p2k3 + 10
     &*p1k2*p2k2*p2k3 + 6*p1k1*p2k3**2 - 6*p1k2*p2k3**2 + k1k3**2*(6*M2 - 2*p1k2 - p1k3 - 6*p1p2 + 2
     &*p2k2 + p2k3) + k1k3*(-16*M2**2 - 6*p1k1*p1k2 - 7*p1k2**2 - 11*p1k1*p1k3 - 18*p1k2*p1k3 + 6*p1k3
     &**2 + 2*p1k1*p1p2 + 10*p1k2*p1p2 + 16*p1p2**2 + 8*p1k2*p2k1 + 11*p1k3*p2k1 + 6*p1p2*p2k1 
     &- 2*p1k1*p2k2 + 8*p1k2*p2k2 + 32*p1k3*p2k2 + 10*p1p2*p2k2 - p2k2**2 + k2k3*(12*M2 + 2*p1k1 - p1k2 
     &- 12*p1p2 - 2*p2k1 + p2k2) - 2*M2*(p1k1 + 3*p1k2 - 2*p1k3 + 3*p2k1 + 7*p2k2 - 2*p2k3) + 3
     &*p1k1*p2k3 - 12*p1k2*p2k3 - 6*p1k3*p2k3 - 8*p1p2*p2k3 - 3*p2k1*p2k3 - 2*p2k2*p2k3) + k1k2*(-8
     &*M2**2 + 8*M2*p1k1 - 12*M2*p1k2 - 7*p1k1*p1k2 - 2*p1k2**2 + 8*M2*p1k3 - 12*p1k1*p1k3 - 17*p1k2
     &*p1k3 + 2*p1k3**2 - 8*p1k1*p1p2 + 8*p1k2*p1p2 - 12*p1k3*p1p2 + 8*p1p2**2 + p1k2*p2k1 + 4*p1k3
     &*p2k1 - 4*M2*p2k2 + 3*p1k1*p2k2 + 29*p1k3*p2k2 + 8*p1p2*p2k2 + 3*p2k1*p2k2 + 2*p2k2**2 + k2k3
     &*(8*M2 - 2*p1k1 - 3*p1k2 + 2*p1k3 - 8*p1p2 + 2*p2k1 + 3*p2k2 - 2*p2k3) + k1k3*(6*M2 - p1k2 + p1k3 
     &- 6*p1p2 + p2k2 - p2k3) + 8*M2*p2k3 + 8*p1k1*p2k3 - 13*p1k2*p2k3 - 2*p1k3*p2k3 - 4*p1p2*p2k3 + p2k2*p2k3))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2*k2k3*(k2k3 - p1k2 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Mttk1 Mvk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*(2 - 3*Nc**2 + Nc**4)*(2*k2k3**2 + 4*p1k2**2 
     &+ 2*p1k2*(14*p1k3 - 2*p2k2 - 7*p2k3) + 2*p1k3*(2*p1k3 - 7*p2k2 - 2*p2k3) + k2k3*(-16*M2 - 3
     &*p1k2 - 3*p1k3 + 16*p1p2 + p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32.
     &*k2k3**2*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(Mttk1 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*Nc**2*(-1 + Nc**2)*(k1k2*(k2k3 + 10*M2 - 15
     &*p1k2 - 6*p1k3 - 10*p1p2 + p2k2) - k1k3*(k2k3 - 8*M2 + 12*p1k2 + 3*p1k3 + 8*p1p2 + p2k3) + 2*(k2k3
     &**2 + 5*p1k2*p1k3 - 5*p1k3**2 - 2*p1k2*p2k1 + 2*p1k3*p2k1 - p1k3*p2k2 + k2k3*(M2 + 9*p1k1 
     &- 6*p1k2 + 3*p1k3 - p1p2 + p2k2) - 4*p1k2*p2k3 + 5*p1k3*p2k3 + p1k1*(-2*p1k2 - 7*p1k3 + 4*p2k2 
     &+ 5*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k3*k2k3*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(Mttk1 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(k1k2*(k2k3 - 8*M2 + 3*p1k2 
     &+ 12*p1k3 + 8*p1p2 + p2k2) - k1k3*(k2k3 + 10*M2 - 6*p1k2 - 15*p1k3 - 10*p1p2 + p2k3) - 2*(k2k3
     &**2 - 5*p1k2**2 + 5*p1k2*p1k3 + 2*p1k2*p2k1 - 2*p1k3*p2k1 + 5*p1k2*p2k2 - 4*p1k3*p2k2 - p1k2
     &*p2k3 + k2k3*(M2 + 9*p1k1 + 3*p1k2 - 6*p1k3 - p1p2 + p2k3) + p1k1*(-7*p1k2 - 2*p1k3 + 5*p2k2 
     &+ 4*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k2*k2k3*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(Mttk1 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(-1 + Nc**2)*(-k1k2 - k1k3 - 2*M2 + p1k1 + 2*p1k2 
     &+ 2*p1k3 + 2*p1p2 - p2k1)*(-2*k1k2**2 + 2*k1k3**2 + k2k3**2 - 2*k2k3*p1k2 + 4*p1k1*p1k2 + 4
     &*p1k2**2 - 3*k2k3*p1k3 - 4*p1k1*p1k3 + 2*p1k3**2 - 4*k2k3*p1p2 + 8*p1k2*p1p2 + 4*p1k3*p1p2 
     &- 4*p1k2*p2k1 + 4*p1k3*p2k1 - 8*M2*p2k2 + 8*p1k1*p2k2 + 4*p1k2*p2k2 + 2*p1k3*p2k2 - 2*k1k2*(2
     &*M2 - 2*p1k1 + p1k2 - 3*p1k3 + 2*p1p2 + p2k2) + k2k3*p2k3 - 4*M2*p2k3 + 4*p1k1*p2k3 + 2*p1k2*p2k3 
     &- 2*p1k3*p2k3 + k1k3*(3*k2k3 + 4*M2 - 4*p1k1 - 6*p1k2 - 4*p1k3 - 8*p1p2 + 2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**5)/(128*k2k3*p1k1*(-k1k2 + p1k1 + p1k2)*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(Mttk1 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2 + k1k3 + 2*M2 - p1k1 - 2*p1k2 
     &- 2*p1k3 - 2*p1p2 + p2k1)*(2*k1k2**2 - 2*k1k3**2 + k2k3**2 - 3*k2k3*p1k2 - 4*p1k1*p1k2 + 2*p1k2
     &**2 - 2*k2k3*p1k3 + 4*p1k1*p1k3 + 4*p1k3**2 - 4*k2k3*p1p2 + 4*p1k2*p1p2 + 8*p1k3*p1p2 + 4
     &*p1k2*p2k1 - 4*p1k3*p2k1 + k2k3*p2k2 - 4*M2*p2k2 + 4*p1k1*p2k2 - 2*p1k2*p2k2 + 2*p1k3*p2k2 + k1k2
     &*(3*k2k3 + 4*M2 - 4*p1k1 - 4*p1k2 - 6*p1k3 - 8*p1p2 + 2*p2k2) - 8*M2*p2k3 + 8*p1k1*p2k3 + 2
     &*p1k2*p2k3 + 4*p1k3*p2k3 - 2*k1k3*(2*M2 - 2*p1k1 - 3*p1k2 + p1k3 + 2*p1p2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**5)/(128*k2k3*p1k1*(k2k3 - p1k2 - p1k3)*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mttk1 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)**2*(k2k3**2 - 4*k2k3*(M2 + p1k2) + 4
     &*(p1k2**2 + p1k2*p1k3 + p1k3**2))*(k2k3 + M2 - p1k2 - p1k3 - p1p2 + p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**5)/(32*k2k3*p1k3*(-k2k3 + p1k2 + p1k3)**2))

!************************************************   2 Re(Mttk1 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)**2*(k2k3**2 - 4*k2k3*(M2 + p1k3) + 4
     &*(p1k2**2 + p1k2*p1k3 + p1k3**2))*(k2k3 + M2 - p1k2 - p1k3 - p1p2 + p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**5)/(32*k2k3*p1k2*(-k2k3 + p1k2 + p1k3)**2))

!************************************************   2 Re(Mttk1 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2 + k1k3 + 2*k2k3 + 2*M2 - p1k1 
     &- 4*p1k2 - 2*p1k3 - 2*p1p2 + p2k1 + 2*p2k2)*(k2k3**2 + k1k3*(k2k3 - 4*M2 + 2*p1k2 - 2*p1k3) 
     &+ k2k3*(-4*p1k1 - 2*p1k2 - 3*p1k3 - 4*p1p2 + p2k3) + 2*(4*p1k1*p1k2 + 2*p1k2**2 + 2*p1k1*p1k3 
     &+ p1k3**2 + k1k2*(-4*M2 + 2*p1k2 + p1k3) + 4*p1k2*p1p2 + 2*p1k3*p1p2 - 4*M2*p2k2 + 2*p1k2*p2k2 
     &+ p1k3*p2k2 - 2*M2*p2k3 + p1k2*p2k3 - p1k3*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*
     &*5)/(128*k2k3*p1k2*(-k1k2 + p1k1 + p1k2)*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(Mttk1 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2 + k1k3 + 2*k2k3 + 2*M2 - p1k1 
     &- 2*p1k2 - 4*p1k3 - 2*p1p2 + p2k1 + 2*p2k3)*(k2k3**2 + k1k2*(k2k3 - 4*M2 - 2*p1k2 + 2*p1k3) 
     &+ k2k3*(-4*p1k1 - 3*p1k2 - 2*p1k3 - 4*p1p2 + p2k2) + 2*(2*p1k1*p1k2 + p1k2**2 + 4*p1k1*p1k3 
     &+ 2*p1k3**2 + k1k3*(-4*M2 + p1k2 + 2*p1k3) + 2*p1k2*p1p2 + 4*p1k3*p1p2 - 2*M2*p2k2 - p1k2*p2k2 
     &+ p1k3*p2k2 - 4*M2*p2k3 + p1k2*p2k3 + 2*p1k3*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*
     &*5)/(128*k2k3*(k2k3 - p1k2 - p1k3)*p1k3*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mttk2 Mttk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mttk2 Mbtk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*Nc**2*(-1 + Nc**2)*(k1k3**3 + k1k3*(16*M2**2 
     &+ 12*p1k1**2 - 32*p1k1*p1k2 - 10*p1k1*p1k3 - 8*p1k2*p1k3 - 7*p1k3**2 + 16*p1k1*p1p2 - 8*p1k2
     &*p1p2 + 4*p1k3*p1p2 - 4*p1k1*p2k1 + 16*p1k2*p2k1 + 6*p1k3*p2k1 - 8*p1k3*p2k2 - 8*p1p2*p2k2 - 4
     &*M2*(6*p1k1 - 6*p1k2 + 3*p1k3 + 4*p1p2 - 2*p2k1 + 2*p2k2 - 3*p2k3) + k2k3*(4*M2 - 2*p1k1 + 3
     &*p1k3 - 4*p2k1 - p2k3) - 8*p1k1*p2k3 + 6*p1k2*p2k3 - 2*p1k3*p2k3 - 4*p1p2*p2k3 + 2*p2k2*p2k3 
     &+ p2k3**2) + 2*k1k2*(k1k3**2 + 2*(8*M2**2 + 6*p1k1**2 + p1k3*(3*p1k3 + 4*p1p2 - p2k1 - 2*p2k3) 
     &- 2*M2*(6*p1k1 + 3*p1k3 + 4*p1p2 - 2*p2k1 - p2k3) + p1k1*(6*p1k3 + 8*p1p2 - 2*p2k1 - p2k3)) 
     &+ k1k3*(4*M2 - 4*p1k1 - 6*p1k3 - 4*p1p2 + p2k3)) + k1k3**2*(-k2k3 + 2*(8*M2 - 5*p1k1 + 3*p1k2 
     &- 2*p1p2 + p2k2 + p2k3)) + 2*(12*p1k1**2*p1k2 + 6*p1k1**2*p1k3 + 12*p1k1*p1k2*p1k3 + 12*p1k1
     &*p1k3**2 + 3*p1k3**3 + 8*p1k1*p1k2*p1p2 + 4*p1k1*p1k3*p1p2 + 4*p1k2*p1k3*p1p2 + 2*p1k3**2*p1p2 
     &- 20*p1k1*p1k2*p2k1 - 10*p1k1*p1k3*p2k1 - 10*p1k2*p1k3*p2k1 - 5*p1k3**2*p2k1 + 4*p1k1**2*p2k2 
     &+ 4*p1k1*p1k3*p2k2 + 4*p1k3**2*p2k2 + 8*p1k1*p1p2*p2k2 + 4*p1k3*p1p2*p2k2 + 4*p1k1*p2k1*p2k2 
     &+ 2*p1k3*p2k1*p2k2 + 2*p1k1**2*p2k3 - 10*p1k1*p1k2*p2k3 - 9*p1k1*p1k3*p2k3 - 2*p1k2*p1k3
     &*p2k3 - 2*p1k3**2*p2k3 + 4*p1k1*p1p2*p2k3 + 2*p1k3*p1p2*p2k3 + 2*p1k1*p2k1*p2k3 + p1k3*p2k1*p2k3 
     &+ 2*p1k1*p2k2*p2k3 - 2*p1k3*p2k2*p2k3 + p1k1*p2k3**2 - p1k3*p2k3**2 - 2*M2*(2*p1k3**2 + 4
     &*p1k1*(2*p1k2 + p1k3) - 2*p1k3*p2k1 + 4*p2k1*p2k2 + p1k2*(4*p1k3 - 4*p2k1 - 2*p2k3) - p1k3*p2k3 
     &+ 2*p2k1*p2k3 + 2*p2k2*p2k3 + p2k3**2) + k2k3*(8*M2**2 + 6*p1k1**2 - p1k1*(2*p1k3 - 8*p1p2 
     &+ 2*p2k1 + p2k3) + 2*M2*(-6*p1k1 - 4*p1p2 + 2*p2k1 + p2k3) + p1k3*(-p1k3 - 2*p1p2 + 5*p2k1 + p2k3))))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(128*k1k3*k2k3*p1k1*(k1k3 - p1k1 - p1k3)))

!************************************************   2 Re(Mttk2 Mbtk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*(-1 + Nc**2)*(k1k2 + k2k3 + 2*M2 - 2*p1k1 - p1k2 
     &- 2*p1k3 - 2*p1p2 + p2k2)*(k1k2*(-k1k3 + 4*p1k1 + 14*p1k3) - k1k3*(k2k3 - 16*M2 - p1k1 
     &+ 16*p1k2 - p1k3 + 16*p1p2 + p2k1 + p2k3) + 2*(7*k2k3*p1k1 - 2*p1k1**2 + 2*k2k3*p1k3 - 14*p1k1
     &*p1k3 - 2*p1k3**2 + 2*p1k1*p2k1 + 7*p1k3*p2k1 + 7*p1k1*p2k3 + 2*p1k3*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(64*k1k3**2*p1k2*(k1k3 - p1k1 - p1k3)))

!************************************************   2 Re(Mttk2 Mbtk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(-k1k3**3 - 2*k1k3**2*(k2k3 
     &+ 8*M2 + 3*p1k2 - 5*p1k3 - 2*p1p2 + p2k1 + p2k2) + k1k3*(-16*M2**2 + 7*p1k1**2 + 8*p1k1*p1k2 
     &+ 10*p1k1*p1k3 + 32*p1k2*p1k3 - 12*p1k3**2 - 4*p1k1*p1p2 + 8*p1k2*p1p2 - 16*p1k3*p1p2 + 2*p1k1
     &*p2k1 - 6*p1k2*p2k1 + 8*p1k3*p2k1 + 4*p1p2*p2k1 - p2k1**2 - 2*k2k3*(4*M2 - 6*p1k1 - 4*p1k3 
     &- 4*p1p2 + p2k1) + 8*p1k1*p2k2 + 8*p1p2*p2k2 - 2*p2k1*p2k2 + 4*M2*(3*p1k1 - 6*p1k2 + 6*p1k3 
     &+ 4*p1p2 - 3*p2k1 + 2*p2k2 - 2*p2k3) - 6*p1k1*p2k3 - 16*p1k2*p2k3 + 4*p1k3*p2k3) - 2*(3*p1k1
     &**3 + 12*p1k1**2*p1k3 + 12*p1k1*p1k2*p1k3 + 6*p1k1*p1k3**2 + 12*p1k2*p1k3**2 + 2*p1k1**2*p1p2 
     &+ 4*p1k1*p1k2*p1p2 + 4*p1k1*p1k3*p1p2 + 8*p1k2*p1k3*p1p2 - 2*p1k1**2*p2k1 - 2*p1k1*p1k2*p2k1 
     &- 9*p1k1*p1k3*p2k1 - 10*p1k2*p1k3*p2k1 + 2*p1k3**2*p2k1 + 2*p1k1*p1p2*p2k1 + 4*p1k3*p1p2*p2k1 
     &- p1k1*p2k1**2 + p1k3*p2k1**2 + 4*p1k1**2*p2k2 + 4*p1k1*p1k3*p2k2 + 4*p1k3**2*p2k2 + 4*p1k1
     &*p1p2*p2k2 + 8*p1k3*p1p2*p2k2 - 2*p1k1*p2k1*p2k2 + 2*p1k3*p2k1*p2k2 + 2*k2k3*(8*M2**2 + 3*p1k1
     &**2 - 2*M2*(3*p1k1 + 6*p1k3 + 4*p1p2 - p2k1 - 2*p2k3) + p1k3*(6*p1k3 + 8*p1p2 - p2k1 - 2*p2k3) 
     &+ p1k1*(6*p1k3 + 4*p1p2 - 2*p2k1 - p2k3)) - 5*p1k1**2*p2k3 - 10*p1k1*p1k2*p2k3 - 10*p1k1*p1k3
     &*p2k3 - 20*p1k2*p1k3*p2k3 + p1k1*p2k1*p2k3 + 2*p1k3*p2k1*p2k3 + 2*p1k1*p2k2*p2k3 + 4*p1k3*p2k2
     &*p2k3 - 2*M2*(2*p1k1**2 + p1k2*(8*p1k3 - 2*p2k1 - 4*p2k3) + p1k1*(4*p1k2 + 4*p1k3 - p2k1 
     &- 2*p2k3) + (p2k1 + 2*p2k2)*(p2k1 + 2*p2k3))) + k1k2*(k1k3**2 + k1k3*(-4*M2 - 3*p1k1 + 2*p1k3 
     &+ p2k1 + 4*p2k3) - 2*(8*M2**2 - p1k1**2 + p1k3*(6*p1k3 + 8*p1p2 - p2k1 - 2*p2k3) + 2*M2*(-6*p1k3 
     &- 4*p1p2 + p2k1 + 2*p2k3) + p1k1*(-2*p1k3 - 2*p1p2 + p2k1 + 5*p2k3))))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)
     &/(128*k1k2*k1k3*(k1k3 - p1k1 - p1k3)*p1k3))

!************************************************   2 Re(Mttk2 Msk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**3*Nc**2*(-1 + Nc**2)*(-8*k2k3*M2**2 + 2*k2k3**2
     &*p1k1 + 8*k2k3*M2*p1k1 + 2*k2k3*p1k1**2 + 8*k2k3*M2*p1k2 - 12*k2k3*p1k1*p1k2 + 20*M2*p1k1*p1k2 
     &+ k2k3**2*p1k3 - 12*k2k3*M2*p1k3 - 17*k2k3*p1k1*p1k3 + 4*M2*p1k1*p1k3 - 6*p1k1**2*p1k3 - 7
     &*k2k3*p1k2*p1k3 + 4*M2*p1k2*p1k3 + 24*p1k1*p1k2*p1k3 - 6*p1k2**2*p1k3 - 2*k2k3*p1k3**2 + 8*M2
     &*p1k3**2 + 24*p1k1*p1k3**2 - 12*p1k2*p1k3**2 - 6*p1k3**3 - 12*k2k3*p1k1*p1p2 - 8*k2k3*p1k2*p1p2 
     &+ 8*k2k3*p1k3*p1p2 + 12*p1k1*p1k3*p1p2 + 12*p1k2*p1k3*p1p2 + 12*p1k3**2*p1p2 + 8*k2k3*p1p2
     &**2 - 2*k2k3**2*p2k1 + 8*k2k3*M2*p2k1 - 2*k2k3*p1k1*p2k1 + 8*k2k3*p1k2*p2k1 - 8*M2*p1k2*p2k1 
     &+ 10*p1k1*p1k2*p2k1 - 13*k2k3*p1k3*p2k1 - 16*M2*p1k3*p2k1 + 8*p1k1*p1k3*p2k1 - 14*p1k2*p1k3*p2k1 
     &- 10*p1k3**2*p2k1 - 4*k2k3*p1p2*p2k1 - 12*p1k2*p1p2*p2k1 + 6*p1k2*p2k1**2 - 6*p1k3*p2k1**2 
     &+ 4*k2k3*p1k1*p2k2 + 8*M2*p1k1*p2k2 - 10*p1k1**2*p2k2 + k2k3*p1k3*p2k2 - 8*M2*p1k3*p2k2 - 38
     &*p1k1*p1k3*p2k2 + 6*p1k2*p1k3*p2k2 + 6*p1k3**2*p2k2 - 28*p1k1*p1p2*p2k2 - 8*p1k3*p1p2*p2k2 + 20
     &*M2*p2k1*p2k2 - 6*p1k1*p2k1*p2k2 + 16*p1k3*p2k1*p2k2 + k1k3**2*(p1k1 + 3*p1k2 + p1k3 - p2k1 
     &- 3*p2k2 - p2k3) - k2k3**2*p2k3 - 4*k2k3*M2*p2k3 + 29*k2k3*p1k1*p2k3 + 4*M2*p1k1*p2k3 + 4*p1k1
     &**2*p2k3 + 3*k2k3*p1k2*p2k3 - 4*M2*p1k2*p2k3 + 4*p1k1*p1k2*p2k3 + 6*p1k2**2*p2k3 - 12*M2*p1k3
     &*p2k3 - 38*p1k1*p1k3*p2k3 + 18*p1k2*p1k3*p2k3 + 12*p1k3**2*p2k3 + 8*k2k3*p1p2*p2k3 - 20*p1k1
     &*p1p2*p2k3 - 12*p1k2*p1p2*p2k3 - 28*p1k3*p1p2*p2k3 + k2k3*p2k1*p2k3 + 16*M2*p2k1*p2k3 - 2*p1k2
     &*p2k1*p2k3 + 10*p1k3*p2k1*p2k3 + 3*k2k3*p2k2*p2k3 + 16*M2*p2k2*p2k3 + 10*p1k1*p2k2*p2k3 - 6
     &*p1k2*p2k2*p2k3 - 2*p1k3*p2k2*p2k3 + 2*k2k3*p2k3**2 + 20*M2*p2k3**2 + 14*p1k1*p2k3**2 - 10*p1k2
     &*p2k3**2 - 6*p1k3*p2k3**2 + k1k2**2*(6*M2 - p1k1 - 2*p1k3 - 6*p1p2 + p2k1 + 2*p2k3) + k1k2
     &*(-16*M2**2 + 4*M2*p1k1 + 6*p1k1**2 - 2*M2*p1k2 - 11*p1k1*p1k2 - 6*M2*p1k3 - 18*p1k1*p1k3 - 6
     &*p1k2*p1k3 - 7*p1k3**2 + 2*p1k2*p1p2 + 10*p1k3*p1p2 + 16*p1p2**2 + 4*M2*p2k1 - 6*p1k1*p2k1 + 3
     &*p1k2*p2k1 - 12*p1k3*p2k1 - 8*p1p2*p2k1 - 6*M2*p2k2 + 11*p1k1*p2k2 + 8*p1k3*p2k2 + 6*p1p2*p2k2 
     &- 3*p2k1*p2k2 - 14*M2*p2k3 + 32*p1k1*p2k3 - 2*p1k2*p2k3 + 8*p1k3*p2k3 + 10*p1p2*p2k3 - 2*p2k1
     &*p2k3 - p2k3**2 + k2k3*(6*M2 + p1k1 - p1k3 - 6*p1p2 - p2k1 + p2k3) + k1k3*(12*M2 + 2*p1k2 
     &- p1k3 - 12*p1p2 - 2*p2k2 + p2k3)) + k1k3*(-8*M2**2 + 6*p1k1**2 - 24*p1k1*p1k2 + 12*p1k2**2 - 19
     &*p1k1*p1k3 + 19*p1k2*p1k3 + 9*p1k3**2 - 18*p1k1*p1p2 - 8*p1k2*p1p2 + 2*p1k3*p1p2 + 8*p1p2*
     &*2 - 6*p1k1*p2k1 + 28*p1k2*p2k1 + 5*p1k3*p2k1 + 2*p1p2*p2k1 - 12*p1k2*p2k2 + 11*p1k3*p2k2 - 4
     &*p1p2*p2k2 - 4*p2k1*p2k2 + 17*p1k1*p2k3 - 29*p1k2*p2k3 - 8*p1k3*p2k3 - 6*p1p2*p2k3 - 3*p2k1*p2k3 
     &- p2k2*p2k3 - p2k3**2 + k2k3*(8*M2 + 2*p1k1 - 2*p1k2 - 3*p1k3 - 8*p1p2 - 2*p2k1 + 2*p2k2 + 3*p2k3) 
     &+ 2*M2*(p1k1 - 2*p1k2 - 3*p1k3 + 7*p2k1 + 8*p2k2 + 5*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)
     &/(128*k1k3*k2k3*(k1k3 - p1k1 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Mttk2 Msk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**3*Nc**2*(-1 + Nc**2)*(2*k1k3**2*(16*M2 + p1k2 
     &- 16*p1p2 - p2k2) + k1k2*(k1k3 - 4*p1k1 - 14*p1k3)*(2*M2 - p1k1 - p1k3 - 2*p1p2 + p2k1 + p2k3) 
     &+ k1k3*(32*M2**2 - 6*p1k1**2 - 19*p1k1*p1k2 - 52*p1k1*p1k3 - 19*p1k2*p1k3 - 6*p1k3**2 - 34*p1k1
     &*p1p2 - 16*p1k2*p1p2 - 34*p1k3*p1p2 - 32*p1p2**2 + 6*p1k1*p2k1 + 17*p1k2*p2k1 + 26*p1k3*p2k1 
     &- 2*p1p2*p2k1 + 3*p1k1*p2k2 + 3*p1k3*p2k2 - 16*p1p2*p2k2 - p2k1*p2k2 + 26*p1k1*p2k3 + 17*p1k2
     &*p2k3 + 6*p1k3*p2k3 - 2*p1p2*p2k3 - p2k2*p2k3 + k2k3*(2*M2 - p1k1 - p1k3 - 2*p1p2 + p2k1 
     &+ p2k3) + 2*M2*(p1k1 + 8*p1k2 + p1k3 + 17*p2k1 + 8*p2k2 + 17*p2k3)) + 2*(2*p1k1**2*p1k2 + 14*p1k1
     &*p1k2*p1k3 + 2*p1k2*p1k3**2 + 4*p1k1**2*p1p2 + 28*p1k1*p1k3*p1p2 + 4*p1k3**2*p1p2 - 2*p1k1
     &*p1k2*p2k1 - 10*p1k1*p1k3*p2k1 - 7*p1k2*p1k3*p2k1 + 10*p1k3**2*p2k1 + 4*p1k1*p1p2*p2k1 + 14*p1k3
     &*p1p2*p2k1 - 2*p1k1**2*p2k2 - 14*p1k1*p1k3*p2k2 - 2*p1k3**2*p2k2 + 2*p1k1*p2k1*p2k2 + 7*p1k3
     &*p2k1*p2k2 + k2k3*(7*p1k1 + 2*p1k3)*(-2*M2 + p1k1 + p1k3 + 2*p1p2 - p2k1 - p2k3) + 10*p1k1
     &**2*p2k3 - 7*p1k1*p1k2*p2k3 - 10*p1k1*p1k3*p2k3 - 2*p1k2*p1k3*p2k3 + 14*p1k1*p1p2*p2k3 + 4*p1k3
     &*p1p2*p2k3 + 7*p1k1*p2k2*p2k3 + 2*p1k3*p2k2*p2k3 - 2*M2*(2*p1k1**2 + p1k3*(2*p1k3 + 7*p2k1 
     &+ 2*p2k3) + p1k1*(14*p1k3 + 2*p2k1 + 7*p2k3))))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)
     &/(64*k1k3**2*(k1k3 - p1k1 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Mttk2 Msk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**3*Nc**2*(-1 + Nc**2)*(6*k2k3**2*M2 - 16*k2k3*M2
     &**2 - 2*k2k3**2*p1k1 - 6*k2k3*M2*p1k1 - 7*k2k3*p1k1**2 + 8*M2*p1k1**2 - 6*p1k1**3 - 2*k2k3*M2
     &*p1k2 - 6*k2k3*p1k1*p1k2 + 4*M2*p1k1*p1k2 - 12*p1k1**2*p1k2 - 6*p1k1*p1k2**2 - k2k3**2*p1k3 
     &+ 4*k2k3*M2*p1k3 - 18*k2k3*p1k1*p1k3 + 4*M2*p1k1*p1k3 + 24*p1k1**2*p1k3 - 11*k2k3*p1k2*p1k3 + 20
     &*M2*p1k2*p1k3 + 24*p1k1*p1k2*p1k3 + 6*k2k3*p1k3**2 - 6*p1k1*p1k3**2 - 6*k2k3**2*p1p2 + 10*k2k3
     &*p1k1*p1p2 + 12*p1k1**2*p1p2 + 2*k2k3*p1k2*p1p2 + 12*p1k1*p1k2*p1p2 + 12*p1k1*p1k3*p1p2 + 16
     &*k2k3*p1p2**2 + 2*k2k3**2*p2k1 - 14*k2k3*M2*p2k1 + 8*k2k3*p1k1*p2k1 - 12*M2*p1k1*p2k1 + 12
     &*p1k1**2*p2k1 - 2*k2k3*p1k2*p2k1 - 4*M2*p1k2*p2k1 + 18*p1k1*p1k2*p2k1 + 6*p1k2**2*p2k1 + 32*k2k3
     &*p1k3*p2k1 + 4*M2*p1k3*p2k1 - 38*p1k1*p1k3*p2k1 + 4*p1k2*p1k3*p2k1 + 4*p1k3**2*p2k1 + 10*k2k3
     &*p1p2*p2k1 - 28*p1k1*p1p2*p2k1 - 12*p1k2*p1p2*p2k1 - 20*p1k3*p1p2*p2k1 - k2k3*p2k1**2 + 20
     &*M2*p2k1**2 - 6*p1k1*p2k1**2 - 10*p1k2*p2k1**2 + 14*p1k3*p2k1**2 - 6*k2k3*M2*p2k2 + 8*k2k3*p1k1
     &*p2k2 - 8*M2*p1k1*p2k2 + 6*p1k1**2*p2k2 + 6*p1k1*p1k2*p2k2 + 11*k2k3*p1k3*p2k2 + 8*M2*p1k3*p2k2 
     &- 38*p1k1*p1k3*p2k2 - 10*p1k3**2*p2k2 + 6*k2k3*p1p2*p2k2 - 8*p1k1*p1p2*p2k2 - 28*p1k3*p1p2
     &*p2k2 + 16*M2*p2k1*p2k2 - 2*p1k1*p2k1*p2k2 - 6*p1k2*p2k1*p2k2 + 10*p1k3*p2k1*p2k2 + k1k2**2
     &*(p1k1 + 2*p1k3 - p2k1 - 2*p2k3) + k1k3**2*(p1k1 + 3*p1k2 + p1k3 - p2k1 - 3*p2k2 - p2k3) + k2k3
     &**2*p2k3 + 4*k2k3*M2*p2k3 - 12*k2k3*p1k1*p2k3 - 16*M2*p1k1*p2k3 - 10*p1k1**2*p2k3 + 3*k2k3*p1k2
     &*p2k3 - 8*M2*p1k2*p2k3 - 14*p1k1*p1k2*p2k3 - 6*k2k3*p1k3*p2k3 + 8*p1k1*p1k3*p2k3 + 10*p1k2
     &*p1k3*p2k3 - 8*k2k3*p1p2*p2k3 - 12*p1k2*p1p2*p2k3 - 2*k2k3*p2k1*p2k3 + 16*M2*p2k1*p2k3 + 10*p1k1
     &*p2k1*p2k3 - 2*p1k2*p2k1*p2k3 - 3*k2k3*p2k2*p2k3 + 20*M2*p2k2*p2k3 + 16*p1k1*p2k2*p2k3 - 6
     &*p1k3*p2k2*p2k3 - 6*p1k1*p2k3**2 + 6*p1k2*p2k3**2 + k1k2*(-8*M2**2 - 12*M2*p1k1 - 2*p1k1**2 
     &+ 8*M2*p1k2 - 7*p1k1*p1k2 + 8*M2*p1k3 - 17*p1k1*p1k3 - 12*p1k2*p1k3 + 2*p1k3**2 + 8*p1k1*p1p2 
     &- 8*p1k2*p1p2 - 12*p1k3*p1p2 + 8*p1p2**2 - 4*M2*p2k1 + 3*p1k2*p2k1 + 29*p1k3*p2k1 + 8*p1p2*p2k1 
     &+ 2*p2k1**2 + p1k1*p2k2 + 4*p1k3*p2k2 + 3*p2k1*p2k2 + k1k3*(8*M2 - 3*p1k1 - 2*p1k2 + 2*p1k3 
     &- 8*p1p2 + 3*p2k1 + 2*p2k2 - 2*p2k3) + k2k3*(6*M2 - p1k1 + p1k3 - 6*p1p2 + p2k1 - p2k3) + 8
     &*M2*p2k3 - 13*p1k1*p2k3 + 8*p1k2*p2k3 - 2*p1k3*p2k3 - 4*p1p2*p2k3 + p2k1*p2k3) + k1k3*(-8*M2
     &**2 + 9*p1k1**2 + 19*p1k1*p1k2 + 12*p1k2**2 - 19*p1k1*p1k3 - 24*p1k2*p1k3 + 6*p1k3**2 + 2*p1k1
     &*p1p2 - 8*p1k2*p1p2 - 18*p1k3*p1p2 + 8*p1p2**2 - 8*p1k1*p2k1 - 29*p1k2*p2k1 + 17*p1k3*p2k1 - 6
     &*p1p2*p2k1 - p2k1**2 + k2k3*(12*M2 - p1k1 + 2*p1k2 - 12*p1p2 + p2k1 - 2*p2k2) + 11*p1k1*p2k2 
     &- 12*p1k2*p2k2 - 4*p1p2*p2k2 - p2k1*p2k2 + 5*p1k1*p2k3 + 28*p1k2*p2k3 - 6*p1k3*p2k3 + 2*p1p2*p2k3 
     &- 3*p2k1*p2k3 - 4*p2k2*p2k3 + 2*M2*(-3*p1k1 - 2*p1k2 + p1k3 + 5*p2k1 + 8*p2k2 + 7*p2k3)))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2*k1k3*(k1k3 - p1k1 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Mttk2 Mvk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*Nc**2*(-1 + Nc**2)*(2*k1k3**2 - k1k3*k2k3 + 8
     &*k2k3*M2 - 12*k2k3*p1k1 - 4*p1k1*p1k2 - 3*k2k3*p1k3 + 10*p1k1*p1k3 - 14*p1k2*p1k3 - 10*p1k3
     &**2 - 8*k2k3*p1p2 + 8*p1k2*p2k1 - 2*p1k3*p2k1 + k1k2*(k1k3 + 10*M2 - 15*p1k1 - 6*p1k3 - 10*p1p2 
     &+ p2k1) + 2*k1k3*(M2 - 6*p1k1 + 9*p1k2 + 3*p1k3 - p1p2 + p2k1) - 4*p1k1*p2k2 + 4*p1k3*p2k2 
     &- k2k3*p2k3 - 8*p1k1*p2k3 + 10*p1k2*p2k3 + 10*p1k3*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k3*k2k3*(k1k3 - p1k1 - p1k3)))

!************************************************   2 Re(Mttk2 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*(2 - 3*Nc**2 + Nc**4)*(2*k1k3**2 + 4*p1k1**2 
     &+ 2*p1k1*(14*p1k3 - 2*p2k1 - 7*p2k3) + 2*p1k3*(2*p1k3 - 7*p2k1 - 2*p2k3) + k1k3*(-16*M2 - 3
     &*p1k1 - 3*p1k3 + 16*p1p2 + p2k1 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32.
     &*k1k3**2*(k1k3 - p1k1 - p1k3)))

!************************************************   2 Re(Mttk2 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(-2*k1k3**2 - 10*k2k3*M2 
     &+ 6*k2k3*p1k1 + 10*p1k1**2 + 14*p1k1*p1k2 + 15*k2k3*p1k3 - 10*p1k1*p1k3 + 4*p1k2*p1k3 + 10*k2k3
     &*p1p2 - 10*p1k1*p2k1 - 10*p1k2*p2k1 + 8*p1k3*p2k1 + k1k2*(k1k3 - 8*M2 + 3*p1k1 + 12*p1k3 + 8
     &*p1p2 + p2k1) - 4*p1k1*p2k2 + 4*p1k3*p2k2 - k2k3*p2k3 + 2*p1k1*p2k3 - 8*p1k2*p2k3 - k1k3*(k2k3 
     &+ 2*(M2 + 3*p1k1 + 9*p1k2 - 6*p1k3 - p1p2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k2*k1k3*(k1k3 - p1k1 - p1k3)))

!************************************************   2 Re(Mttk2 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2 + 2*k1k3 + k2k3 + 2*M2 - 4*p1k1 
     &- p1k2 - 2*p1k3 - 2*p1p2 + 2*p2k1 + p2k2)*(k1k3**2 + k1k3*(k2k3 - 2*p1k1 - 4*p1k2 - 3*p1k3 
     &- 4*p1p2 + p2k3) + 2*(2*p1k1**2 + 4*p1k1*p1k2 + k2k3*(-2*M2 + p1k1 - p1k3) + 2*p1k2*p1k3 + p1k3
     &**2 + k1k2*(-4*M2 + 2*p1k1 + p1k3) + 4*p1k1*p1p2 + 2*p1k3*p1p2 - 4*M2*p2k1 + 2*p1k1*p2k1 
     &+ p1k3*p2k1 - 2*M2*p2k3 + p1k1*p2k3 - p1k3*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(128.
     &*k1k3*p1k1*(-k1k2 + p1k1 + p1k2)*(k1k3 - p1k1 - p1k3)))

!************************************************   2 Re(Mttk2 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)**2*(k1k3**2 - 4*k1k3*(M2 + p1k3) + 4
     &*(p1k1**2 + p1k1*p1k3 + p1k3**2))*(k1k3 + M2 - p1k1 - p1k3 - p1p2 + p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**5)/(32*k1k3*p1k1*(-k1k3 + p1k1 + p1k3)**2))

!************************************************   2 Re(Mttk2 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2 + 2*k1k3 + k2k3 + 2*M2 - 2*p1k1 
     &- p1k2 - 4*p1k3 - 2*p1p2 + p2k2 + 2*p2k3)*(k1k3**2 + k1k2*(k1k3 - 4*M2 - 2*p1k1 + 2*p1k3) 
     &+ k1k3*(-3*p1k1 - 4*p1k2 - 2*p1k3 - 4*p1p2 + p2k1) + 2*(p1k1**2 + 4*p1k2*p1k3 + 2*p1k3**2 + k2k3
     &*(-4*M2 + p1k1 + 2*p1k3) + 4*p1k3*p1p2 - 2*M2*p2k1 + p1k3*p2k1 - 4*M2*p2k3 + 2*p1k3*p2k3 
     &+ p1k1*(2*p1k2 + 2*p1p2 - p2k1 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(128*k1k3
     &*(k1k3 - p1k1 - p1k3)*p1k3*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mttk2 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2 + k2k3 + 2*M2 - 2*p1k1 - p1k2 
     &- 2*p1k3 - 2*p1p2 + p2k2)*(2*k1k2**2 + k1k3**2 + k1k3*(-3*p1k1 - 2*p1k3 - 4*p1p2 + p2k1) + k1k2
     &*(3*k1k3 + 4*M2 - 4*p1k1 - 4*p1k2 - 6*p1k3 - 8*p1p2 + 2*p2k1) + 2*(-k2k3**2 + p1k1**2 + 2
     &*p1k2*p1k3 + 2*p1k3**2 + 4*p1k3*p1p2 - 2*M2*p2k1 + 2*p1k2*p2k1 + p1k3*p2k1 - 2*p1k3*p2k2 - 4*M2
     &*p2k3 + 4*p1k2*p2k3 + 2*p1k3*p2k3 - k2k3*(2*M2 - 3*p1k1 - 2*p1k2 + p1k3 + 2*p1p2 + p2k3) + p1k1
     &*(-2*p1k2 + 2*p1p2 - p2k1 + 2*p2k2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(128.
     &*k1k3*p1k2*(k1k3 - p1k1 - p1k3)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mttk2 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(-1 + Nc**2)*(-k1k2 - k2k3 - 2*M2 + 2*p1k1 + p1k2 
     &+ 2*p1k3 + 2*p1p2 - p2k2)*(-2*k1k2**2 + k1k3**2 - 2*k1k2*(2*M2 + p1k1 - 2*p1k2 - 3*p1k3 + 2
     &*p1p2 + p2k1) + k1k3*(3*k2k3 - 2*p1k1 - 3*p1k3 - 4*p1p2 + p2k3) + 2*(k2k3**2 + 2*p1k1**2 - 2
     &*p1k2*p1k3 + p1k3**2 + 2*p1k3*p1p2 - 4*M2*p2k1 + 4*p1k2*p2k1 + p1k3*p2k1 + 2*p1k3*p2k2 - 2*M2
     &*p2k3 + 2*p1k2*p2k3 - p1k3*p2k3 + k2k3*(2*M2 - 3*p1k1 - 2*p1k2 - 2*p1k3 - 4*p1p2 + p2k3) + p1k1
     &*(2*p1k2 + 4*p1p2 + 2*p2k1 - 2*p2k2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(128.
     &*k1k3*p1k2*(-k1k2 + p1k1 + p1k2)*(k1k3 - p1k1 - p1k3)))

!************************************************   2 Re(Mttk2 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)**2*(k1k3**2 - 4*k1k3*(M2 + p1k1) + 4
     &*(p1k1**2 + p1k1*p1k3 + p1k3**2))*(k1k3 + M2 - p1k1 - p1k3 - p1p2 + p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**5)/(32*k1k3*p1k3*(-k1k3 + p1k1 + p1k3)**2))

!************************************************   2 Re(Mttk3 Mbtk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*Nc**2*(-1 + Nc**2)*(k1k2**3 + k1k2*(16*M2**2 
     &- 24*M2*p1k1 + 12*p1k1**2 - 12*M2*p1k2 - 10*p1k1*p1k2 - 7*p1k2**2 + 24*M2*p1k3 - 32*p1k1*p1k3 
     &- 8*p1k2*p1k3 - 16*M2*p1p2 + 16*p1k1*p1p2 + 4*p1k2*p1p2 - 8*p1k3*p1p2 + 8*M2*p2k1 - 4*p1k1
     &*p2k1 + 6*p1k2*p2k1 + 16*p1k3*p2k1 + k2k3*(4*M2 - 2*p1k1 + 3*p1k2 - 4*p2k1 - p2k2) + 12*M2*p2k2 
     &- 8*p1k1*p2k2 - 2*p1k2*p2k2 + 6*p1k3*p2k2 - 4*p1p2*p2k2 + p2k2**2 + 2*k1k3*(4*M2 - 4*p1k1 
     &- 6*p1k2 - 4*p1p2 + p2k2) - 8*M2*p2k3 - 8*p1k2*p2k3 - 8*p1p2*p2k3 + 2*p2k2*p2k3) + 2*(-8*M2*p1k1
     &*p1k2 + 6*p1k1**2*p1k2 - 4*M2*p1k2**2 + 12*p1k1*p1k2**2 + 3*p1k2**3 - 16*M2*p1k1*p1k3 + 12
     &*p1k1**2*p1k3 - 8*M2*p1k2*p1k3 + 12*p1k1*p1k2*p1k3 + 4*p1k1*p1k2*p1p2 + 2*p1k2**2*p1p2 + 8*p1k1
     &*p1k3*p1p2 + 4*p1k2*p1k3*p1p2 + 4*M2*p1k2*p2k1 - 10*p1k1*p1k2*p2k1 - 5*p1k2**2*p2k1 + 8*M2*p1k3
     &*p2k1 - 20*p1k1*p1k3*p2k1 - 10*p1k2*p1k3*p2k1 + 2*k1k3*(8*M2**2 + 6*p1k1**2 + p1k2*(3*p1k2 
     &+ 4*p1p2 - p2k1 - 2*p2k2) - 2*M2*(6*p1k1 + 3*p1k2 + 4*p1p2 - 2*p2k1 - p2k2) + p1k1*(6*p1k2 
     &+ 8*p1p2 - 2*p2k1 - p2k2)) + 2*p1k1**2*p2k2 + 2*M2*p1k2*p2k2 - 9*p1k1*p1k2*p2k2 - 2*p1k2**2*p2k2 
     &+ 4*M2*p1k3*p2k2 - 10*p1k1*p1k3*p2k2 - 2*p1k2*p1k3*p2k2 + 4*p1k1*p1p2*p2k2 + 2*p1k2*p1p2*p2k2 
     &- 4*M2*p2k1*p2k2 + 2*p1k1*p2k1*p2k2 + p1k2*p2k1*p2k2 - 2*M2*p2k2**2 + p1k1*p2k2**2 - p1k2
     &*p2k2**2 + k2k3*(8*M2**2 + 6*p1k1**2 - p1k1*(2*p1k2 - 8*p1p2 + 2*p2k1 + p2k2) + 2*M2*(-6*p1k1 
     &- 4*p1p2 + 2*p2k1 + p2k2) + p1k2*(-p1k2 - 2*p1p2 + 5*p2k1 + p2k2)) + 4*p1k1**2*p2k3 + 4*p1k1
     &*p1k2*p2k3 + 4*p1k2**2*p2k3 + 8*p1k1*p1p2*p2k3 + 4*p1k2*p1p2*p2k3 - 8*M2*p2k1*p2k3 + 4*p1k1*p2k1
     &*p2k3 + 2*p1k2*p2k1*p2k3 - 4*M2*p2k2*p2k3 + 2*p1k1*p2k2*p2k3 - 2*p1k2*p2k2*p2k3) + k1k2**2
     &*(2*k1k3 - k2k3 + 2*(8*M2 - 5*p1k1 + 3*p1k3 - 2*p1p2 + p2k2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)
     &/(128*k1k2*k2k3*p1k1*(k1k2 - p1k1 - p1k2)))

!************************************************   2 Re(Mttk3 Mbtk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*(Nc**2 - Nc**4)*(-4*(k1k2 - 2*p1k1 - p1k2)*(
     &-2*k2k3*M2 + k2k3*p1k1 - M2*p1k1 + 2*p1k1**2 + 2*k2k3*p1k2 + p1k1*p1k2 - 2*M2*p1k3 + p1k1*p1k3 
     &+ 2*p1k2*p1k3 + 2*k2k3*p1p2 + p1k1*p1p2 + 2*p1k3*p1p2 + k1k3*(-M2 + p1k2 + p1p2) + k1k2*(-M2 
     &- 2*p1k1 + p1k2 - 2*p1k3 + p1p2) - M2*p2k1 - 2*p1k1*p2k1 + p1k2*p2k1 - 2*p1k3*p2k1 + p1p2*p2k1 
     &- 2*p1k1*p2k2 - 4*p1k3*p2k2 - 2*M2*p2k3 + p1k1*p2k3 + 2*p1k2*p2k3 + 2*p1p2*p2k3) - 2*(p1k1 
     &- p1k2)*(2*k1k2**2 - 2*M2*p1k1 + 4*p1k1**2 + 2*k2k3*p1k2 + p1k1*p1k2 - 4*M2*p1k3 + 2*p1k1*p1k3 
     &+ 2*p1k2*p1k3 + 2*p1k1*p1p2 + 4*p1k3*p1p2 + 2*k2k3*p2k1 - 2*M2*p2k1 - 4*p1k1*p2k1 + 3*p1k2
     &*p2k1 - 4*p1k3*p2k1 + 2*p1p2*p2k1 + 2*k2k3*p2k2 - 3*p1k1*p2k2 - 6*p1k3*p2k2 - p2k1*p2k2 - k1k3
     &*(2*M2 - 3*p1k2 - 2*p1p2 + p2k2) + k1k2*(k2k3 - 6*p1k1 + p1k2 - 3*p1k3 + 2*p2k1 + p2k2 - p2k3) 
     &- 4*M2*p2k3 + 2*p1k1*p2k3 + 6*p1k2*p2k3 + 4*p1p2*p2k3 - 2*p2k2*p2k3) + (k1k2 - 4*M2 - 2*p1k1 
     &+ 2*p1k2)*(k1k2**2 - k1k2*k1k3 + k1k3*(-4*M2 + p1k1 + 4*p1k2 + 4*p1p2 - p2k1 - 4*p2k2) + (p1k1 
     &- p2k1)*(-2*k2k3 + p1k1 + 2*p1k3 - p2k1 - 2*p2k3) + 2*k1k2*(k2k3 - p1k1 - p1k3 + p2k1 + p2k3)) 
     &+ 2*(4*M2 - 3*p1k2)*(2*k1k2**2 + 2*k2k3*(2*M2 - p1k1 - p1k2 - 2*p1p2 + p2k1 + p2k2) + k1k2
     &*(k2k3 + 2*M2 - 2*p1k1 - p1k2 + p1k3 - 2*p1p2 + 2*p2k1 + p2k2 - p2k3) + (p1k2 - p2k2)*(k1k3 
     &- p1k1 - 2*p1k3 + p2k1 + 2*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(128*k1k2*k1k3*(k1k2 - p1k1 - p1k2)*p1k2))

!************************************************   2 Re(Mttk3 Mbtk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*(-1 + Nc**2)*(k1k2*(k1k3 + k2k3 - 16*M2 - p1k1 
     &- p1k2 + 16*p1k3 + 16*p1p2 + p2k1 + p2k2) - 2*(2*k1k3*p1k1 + 7*k2k3*p1k1 - 2*p1k1**2 + 7*k1k3
     &*p1k2 + 2*k2k3*p1k2 - 14*p1k1*p1k2 - 2*p1k2**2 + 2*p1k1*p2k1 + 7*p1k2*p2k1 + 7*p1k1*p2k2 + 2
     &*p1k2*p2k2))*(k1k3 + k2k3 + 2*M2 - 2*p1k1 - 2*p1k2 - p1k3 - 2*p1p2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(64*k1k2**2*(k1k2 - p1k1 - p1k2)*p1k3))

!************************************************   2 Re(Mttk3 Msk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**3*Nc**2*(-1 + Nc**2)*(-8*k2k3*M2**2 + 2*k2k3**2
     &*p1k1 + 8*k2k3*M2*p1k1 + 2*k2k3*p1k1**2 + k2k3**2*p1k2 - 12*k2k3*M2*p1k2 - 17*k2k3*p1k1*p1k2 
     &+ 4*M2*p1k1*p1k2 - 6*p1k1**2*p1k2 - 2*k2k3*p1k2**2 + 8*M2*p1k2**2 + 24*p1k1*p1k2**2 - 6*p1k2
     &**3 + 8*k2k3*M2*p1k3 - 12*k2k3*p1k1*p1k3 + 20*M2*p1k1*p1k3 - 7*k2k3*p1k2*p1k3 + 4*M2*p1k2*p1k3 
     &+ 24*p1k1*p1k2*p1k3 - 12*p1k2**2*p1k3 - 6*p1k2*p1k3**2 - 12*k2k3*p1k1*p1p2 + 8*k2k3*p1k2*p1p2 
     &+ 12*p1k1*p1k2*p1p2 + 12*p1k2**2*p1p2 - 8*k2k3*p1k3*p1p2 + 12*p1k2*p1k3*p1p2 + 8*k2k3*p1p2
     &**2 - 2*k2k3**2*p2k1 + 8*k2k3*M2*p2k1 - 2*k2k3*p1k1*p2k1 - 13*k2k3*p1k2*p2k1 - 16*M2*p1k2*p2k1 
     &+ 8*p1k1*p1k2*p2k1 - 10*p1k2**2*p2k1 + 8*k2k3*p1k3*p2k1 - 8*M2*p1k3*p2k1 + 10*p1k1*p1k3*p2k1 
     &- 14*p1k2*p1k3*p2k1 - 4*k2k3*p1p2*p2k1 - 12*p1k3*p1p2*p2k1 - 6*p1k2*p2k1**2 + 6*p1k3*p2k1**2 
     &- k2k3**2*p2k2 - 4*k2k3*M2*p2k2 + 29*k2k3*p1k1*p2k2 + 4*M2*p1k1*p2k2 + 4*p1k1**2*p2k2 - 12*M2
     &*p1k2*p2k2 - 38*p1k1*p1k2*p2k2 + 12*p1k2**2*p2k2 + 3*k2k3*p1k3*p2k2 - 4*M2*p1k3*p2k2 + 4*p1k1
     &*p1k3*p2k2 + 18*p1k2*p1k3*p2k2 + 6*p1k3**2*p2k2 + 8*k2k3*p1p2*p2k2 - 20*p1k1*p1p2*p2k2 - 28
     &*p1k2*p1p2*p2k2 - 12*p1k3*p1p2*p2k2 + k2k3*p2k1*p2k2 + 16*M2*p2k1*p2k2 + 10*p1k2*p2k1*p2k2 - 2
     &*p1k3*p2k1*p2k2 + 2*k2k3*p2k2**2 + 20*M2*p2k2**2 + 14*p1k1*p2k2**2 - 6*p1k2*p2k2**2 - 10*p1k3
     &*p2k2**2 + k1k3**2*(6*M2 - p1k1 - 2*p1k2 - 6*p1p2 + p2k1 + 2*p2k2) + k1k2**2*(p1k1 + p1k2 + 3
     &*p1k3 - p2k1 - p2k2 - 3*p2k3) + 4*k2k3*p1k1*p2k3 + 8*M2*p1k1*p2k3 - 10*p1k1**2*p2k3 + k2k3*p1k2
     &*p2k3 - 8*M2*p1k2*p2k3 - 38*p1k1*p1k2*p2k3 + 6*p1k2**2*p2k3 + 6*p1k2*p1k3*p2k3 - 28*p1k1*p1p2
     &*p2k3 - 8*p1k2*p1p2*p2k3 + 20*M2*p2k1*p2k3 - 6*p1k1*p2k1*p2k3 + 16*p1k2*p2k1*p2k3 + 3*k2k3
     &*p2k2*p2k3 + 16*M2*p2k2*p2k3 + 10*p1k1*p2k2*p2k3 - 2*p1k2*p2k2*p2k3 - 6*p1k3*p2k2*p2k3 + k1k3
     &*(-16*M2**2 + 6*p1k1**2 - 18*p1k1*p1k2 - 7*p1k2**2 - 11*p1k1*p1k3 - 6*p1k2*p1k3 + 10*p1k2*p1p2 
     &+ 2*p1k3*p1p2 + 16*p1p2**2 - 6*p1k1*p2k1 - 12*p1k2*p2k1 + 3*p1k3*p2k1 - 8*p1p2*p2k1 + 32*p1k1
     &*p2k2 + 8*p1k2*p2k2 - 2*p1k3*p2k2 + 10*p1p2*p2k2 - 2*p2k1*p2k2 - p2k2**2 + k2k3*(6*M2 + p1k1 
     &- p1k2 - 6*p1p2 - p2k1 + p2k2) + 2*M2*(2*p1k1 - 3*p1k2 - p1k3 + 2*p2k1 - 7*p2k2 - 3*p2k3) + 11
     &*p1k1*p2k3 + 8*p1k2*p2k3 + 6*p1p2*p2k3 - 3*p2k1*p2k3) + k1k2*(-8*M2**2 + 2*M2*p1k1 + 6*p1k1
     &**2 - 6*M2*p1k2 - 19*p1k1*p1k2 + 9*p1k2**2 - 4*M2*p1k3 - 24*p1k1*p1k3 + 19*p1k2*p1k3 + 12*p1k3
     &**2 - 18*p1k1*p1p2 + 2*p1k2*p1p2 - 8*p1k3*p1p2 + 8*p1p2**2 + 14*M2*p2k1 - 6*p1k1*p2k1 + 5*p1k2
     &*p2k1 + 28*p1k3*p2k1 + 2*p1p2*p2k1 + 10*M2*p2k2 + 17*p1k1*p2k2 - 8*p1k2*p2k2 - 29*p1k3*p2k2 
     &- 6*p1p2*p2k2 - 3*p2k1*p2k2 - p2k2**2 + k1k3*(12*M2 - p1k2 + 2*p1k3 - 12*p1p2 + p2k2 - 2*p2k3) 
     &+ 16*M2*p2k3 + 11*p1k2*p2k3 - 12*p1k3*p2k3 - 4*p1p2*p2k3 - 4*p2k1*p2k3 - p2k2*p2k3 + k2k3*(8*M2 
     &+ 2*p1k1 - 3*p1k2 - 2*p1k3 - 8*p1p2 - 2*p2k1 + 3*p2k2 + 2*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)
     &/(128*k1k2*k2k3*(k1k2 - p1k1 - p1k2)*(M2 + p1p2)))

!************************************************   2 Re(Mttk3 Msk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**3*Nc**2*(-1 + Nc**2)*(6*k2k3**2*M2 - 16*k2k3*M2
     &**2 - 2*k2k3**2*p1k1 - 6*k2k3*M2*p1k1 - 7*k2k3*p1k1**2 + 8*M2*p1k1**2 - 6*p1k1**3 - k2k3**2*p1k2 
     &+ 4*k2k3*M2*p1k2 - 18*k2k3*p1k1*p1k2 + 4*M2*p1k1*p1k2 + 24*p1k1**2*p1k2 + 6*k2k3*p1k2**2 
     &- 6*p1k1*p1k2**2 - 2*k2k3*M2*p1k3 - 6*k2k3*p1k1*p1k3 + 4*M2*p1k1*p1k3 - 12*p1k1**2*p1k3 - 11
     &*k2k3*p1k2*p1k3 + 20*M2*p1k2*p1k3 + 24*p1k1*p1k2*p1k3 - 6*p1k1*p1k3**2 - 6*k2k3**2*p1p2 + 10*k2k3
     &*p1k1*p1p2 + 12*p1k1**2*p1p2 + 12*p1k1*p1k2*p1p2 + 2*k2k3*p1k3*p1p2 + 12*p1k1*p1k3*p1p2 + 16
     &*k2k3*p1p2**2 + 2*k2k3**2*p2k1 - 14*k2k3*M2*p2k1 + 8*k2k3*p1k1*p2k1 - 12*M2*p1k1*p2k1 + 12
     &*p1k1**2*p2k1 + 32*k2k3*p1k2*p2k1 + 4*M2*p1k2*p2k1 - 38*p1k1*p1k2*p2k1 + 4*p1k2**2*p2k1 - 2*k2k3
     &*p1k3*p2k1 - 4*M2*p1k3*p2k1 + 18*p1k1*p1k3*p2k1 + 4*p1k2*p1k3*p2k1 + 6*p1k3**2*p2k1 + 10*k2k3
     &*p1p2*p2k1 - 28*p1k1*p1p2*p2k1 - 20*p1k2*p1p2*p2k1 - 12*p1k3*p1p2*p2k1 - k2k3*p2k1**2 + 20
     &*M2*p2k1**2 - 6*p1k1*p2k1**2 + 14*p1k2*p2k1**2 - 10*p1k3*p2k1**2 + k1k3**2*(p1k1 + 2*p1k2 - p2k1 
     &- 2*p2k2) + k2k3**2*p2k2 + 4*k2k3*M2*p2k2 - 12*k2k3*p1k1*p2k2 - 16*M2*p1k1*p2k2 - 10*p1k1
     &**2*p2k2 - 6*k2k3*p1k2*p2k2 + 8*p1k1*p1k2*p2k2 + 3*k2k3*p1k3*p2k2 - 8*M2*p1k3*p2k2 - 14*p1k1*p1k3
     &*p2k2 + 10*p1k2*p1k3*p2k2 - 8*k2k3*p1p2*p2k2 - 12*p1k3*p1p2*p2k2 - 2*k2k3*p2k1*p2k2 + 16*M2
     &*p2k1*p2k2 + 10*p1k1*p2k1*p2k2 - 2*p1k3*p2k1*p2k2 - 6*p1k1*p2k2**2 + 6*p1k3*p2k2**2 + k1k2*
     &*2*(p1k1 + p1k2 + 3*p1k3 - p2k1 - p2k2 - 3*p2k3) - 6*k2k3*M2*p2k3 + 8*k2k3*p1k1*p2k3 - 8*M2*p1k1
     &*p2k3 + 6*p1k1**2*p2k3 + 11*k2k3*p1k2*p2k3 + 8*M2*p1k2*p2k3 - 38*p1k1*p1k2*p2k3 - 10*p1k2*
     &*2*p2k3 + 6*p1k1*p1k3*p2k3 + 6*k2k3*p1p2*p2k3 - 8*p1k1*p1p2*p2k3 - 28*p1k2*p1p2*p2k3 + 16*M2*p2k1
     &*p2k3 - 2*p1k1*p2k1*p2k3 + 10*p1k2*p2k1*p2k3 - 6*p1k3*p2k1*p2k3 - 3*k2k3*p2k2*p2k3 + 20*M2
     &*p2k2*p2k3 + 16*p1k1*p2k2*p2k3 - 6*p1k2*p2k2*p2k3 + k1k3*(-8*M2**2 - 2*p1k1**2 - 17*p1k1*p1k2 
     &+ 2*p1k2**2 - 7*p1k1*p1k3 - 12*p1k2*p1k3 + 8*p1k1*p1p2 - 12*p1k2*p1p2 - 8*p1k3*p1p2 + 8*p1p2
     &**2 + 29*p1k2*p2k1 + 3*p1k3*p2k1 + 8*p1p2*p2k1 + 2*p2k1**2 - 4*M2*(3*p1k1 - 2*p1k2 - 2*p1k3 
     &+ p2k1 - 2*p2k2) + k2k3*(6*M2 - p1k1 + p1k2 - 6*p1p2 + p2k1 - p2k2) - 13*p1k1*p2k2 - 2*p1k2*p2k2 
     &+ 8*p1k3*p2k2 - 4*p1p2*p2k2 + p2k1*p2k2 + p1k1*p2k3 + 4*p1k2*p2k3 + 3*p2k1*p2k3) + k1k2*(
     &-8*M2**2 - 6*M2*p1k1 + 9*p1k1**2 + 2*M2*p1k2 - 19*p1k1*p1k2 + 6*p1k2**2 - 4*M2*p1k3 + 19*p1k1
     &*p1k3 - 24*p1k2*p1k3 + 12*p1k3**2 + 2*p1k1*p1p2 - 18*p1k2*p1p2 - 8*p1k3*p1p2 + 8*p1p2**2 + 10
     &*M2*p2k1 - 8*p1k1*p2k1 + 17*p1k2*p2k1 - 29*p1k3*p2k1 - 6*p1p2*p2k1 - p2k1**2 + 14*M2*p2k2 + 5
     &*p1k1*p2k2 - 6*p1k2*p2k2 + 28*p1k3*p2k2 + 2*p1p2*p2k2 - 3*p2k1*p2k2 + k2k3*(12*M2 - p1k1 + 2*p1k3 
     &- 12*p1p2 + p2k1 - 2*p2k3) + 16*M2*p2k3 + 11*p1k1*p2k3 - 12*p1k3*p2k3 - 4*p1p2*p2k3 - p2k1*p2k3 
     &- 4*p2k2*p2k3 + k1k3*(8*M2 - 3*p1k1 + 2*p1k2 - 2*p1k3 - 8*p1p2 + 3*p2k1 - 2*p2k2 + 2*p2k3)))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2*k1k3*(k1k2 - p1k1 - p1k2)*(M2 + p1p2)))

!************************************************   2 Re(Mttk3 Msk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**3*Nc**2*(-1 + Nc**2)*(2*k1k2**2*(16*M2 + p1k3 
     &- 16*p1p2 - p2k3) + k1k2*(32*M2**2 + 2*M2*p1k1 - 6*p1k1**2 + 2*M2*p1k2 - 52*p1k1*p1k2 - 6*p1k2
     &**2 + 16*M2*p1k3 - 19*p1k1*p1k3 - 19*p1k2*p1k3 - 34*p1k1*p1p2 - 34*p1k2*p1p2 - 16*p1k3*p1p2 
     &- 32*p1p2**2 + 34*M2*p2k1 + 6*p1k1*p2k1 + 26*p1k2*p2k1 + 17*p1k3*p2k1 - 2*p1p2*p2k1 + 34*M2*p2k2 
     &+ 26*p1k1*p2k2 + 6*p1k2*p2k2 + 17*p1k3*p2k2 - 2*p1p2*p2k2 + k1k3*(2*M2 - p1k1 - p1k2 - 2*p1p2 
     &+ p2k1 + p2k2) + k2k3*(2*M2 - p1k1 - p1k2 - 2*p1p2 + p2k1 + p2k2) + 16*M2*p2k3 + 3*p1k1*p2k3 
     &+ 3*p1k2*p2k3 - 16*p1p2*p2k3 - p2k1*p2k3 - p2k2*p2k3) - 2*(4*M2*p1k1**2 + 28*M2*p1k1*p1k2 
     &+ 4*M2*p1k2**2 - 2*p1k1**2*p1k3 - 14*p1k1*p1k2*p1k3 - 2*p1k2**2*p1k3 - 4*p1k1**2*p1p2 - 28*p1k1
     &*p1k2*p1p2 - 4*p1k2**2*p1p2 + 4*M2*p1k1*p2k1 + 14*M2*p1k2*p2k1 + 10*p1k1*p1k2*p2k1 - 10*p1k2
     &**2*p2k1 + 2*p1k1*p1k3*p2k1 + 7*p1k2*p1k3*p2k1 - 4*p1k1*p1p2*p2k1 - 14*p1k2*p1p2*p2k1 + 14
     &*M2*p1k1*p2k2 - 10*p1k1**2*p2k2 + 4*M2*p1k2*p2k2 + 10*p1k1*p1k2*p2k2 + 7*p1k1*p1k3*p2k2 + 2*p1k2
     &*p1k3*p2k2 - 14*p1k1*p1p2*p2k2 - 4*p1k2*p1p2*p2k2 + k2k3*(7*p1k1 + 2*p1k2)*(2*M2 - p1k1 - p1k2 
     &- 2*p1p2 + p2k1 + p2k2) + k1k3*(2*p1k1 + 7*p1k2)*(2*M2 - p1k1 - p1k2 - 2*p1p2 + p2k1 + p2k2) 
     &+ 2*p1k1**2*p2k3 + 14*p1k1*p1k2*p2k3 + 2*p1k2**2*p2k3 - 2*p1k1*p2k1*p2k3 - 7*p1k2*p2k1*p2k3 
     &- 7*p1k1*p2k2*p2k3 - 2*p1k2*p2k2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(64*k1k2**2*(k1k2 - p1k1 - p1k2)*(M2 + p1p2)))

!************************************************   2 Re(Mttk3 Mvk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*Nc**2*(-1 + Nc**2)*(2*k1k2**2 + 8*k2k3*M2 - 12
     &*k2k3*p1k1 - 3*k2k3*p1k2 + 10*p1k1*p1k2 - 10*p1k2**2 - 4*p1k1*p1k3 - 14*p1k2*p1k3 - 8*k2k3
     &*p1p2 - 2*p1k2*p2k1 + 8*p1k3*p2k1 + k1k3*(10*M2 - 15*p1k1 - 6*p1k2 - 10*p1p2 + p2k1) + k1k2*(k1k3 
     &- k2k3 + 2*(M2 - 6*p1k1 + 3*p1k2 + 9*p1k3 - p1p2 + p2k1)) - k2k3*p2k2 - 8*p1k1*p2k2 + 10
     &*p1k2*p2k2 + 10*p1k3*p2k2 - 4*p1k1*p2k3 + 4*p1k2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k2*k2k3*(k1k2 - p1k1 - p1k2)))

!************************************************   2 Re(Mttk3 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*Nc**2*(-1 + Nc**2)*(2*k1k2**2 + 10*k2k3*M2 
     &- 6*k2k3*p1k1 - 10*p1k1**2 - 15*k2k3*p1k2 + 10*p1k1*p1k2 - 14*p1k1*p1k3 - 4*p1k2*p1k3 - 10*k2k3
     &*p1p2 + k1k3*(8*M2 - 3*p1k1 - 12*p1k2 - 8*p1p2 - p2k1) + 10*p1k1*p2k1 - 8*p1k2*p2k1 + 10*p1k3
     &*p2k1 + k2k3*p2k2 - 2*p1k1*p2k2 + 8*p1k3*p2k2 + k1k2*(-k1k3 + k2k3 + 2*(M2 + 3*p1k1 - 6*p1k2 
     &+ 9*p1k3 - p1p2 + p2k2)) + 4*p1k1*p2k3 - 4*p1k2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k2*k1k3*(k1k2 - p1k1 - p1k2)))

!************************************************   2 Re(Mttk3 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*(2 - 3*Nc**2 + Nc**4)*(2*k1k2**2 + 4*p1k1**2 
     &+ 2*p1k1*(14*p1k2 - 2*p2k1 - 7*p2k2) + 2*p1k2*(2*p1k2 - 7*p2k1 - 2*p2k2) + k1k2*(-16*M2 - 3
     &*p1k1 - 3*p1k2 + 16*p1p2 + p2k1 + p2k2))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32.
     &*k1k2**2*(k1k2 - p1k1 - p1k2)))

!************************************************   2 Re(Mttk3 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)**2*(k1k2**2 - 4*k1k2*(M2 + p1k2) + 4
     &*(p1k1**2 + p1k1*p1k2 + p1k2**2))*(k1k2 + M2 - p1k1 - p1k2 - p1p2 + p2k1 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**5)/(32*k1k2*p1k1*(-k1k2 + p1k1 + p1k2)**2))

!************************************************   2 Re(Mttk3 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2**2 + k1k2*(k2k3 - 2*p1k1 - 3
     &*p1k2 - 4*p1k3 - 4*p1p2 + p2k2) + 2*(2*p1k1**2 + k2k3*(-2*M2 + p1k1 - p1k2) + p1k2**2 + k1k3*(
     &-4*M2 + 2*p1k1 + p1k2) + 4*p1k1*p1k3 + 2*p1k2*p1k3 + 4*p1k1*p1p2 + 2*p1k2*p1p2 - 4*M2*p2k1 + 2
     &*p1k1*p2k1 + p1k2*p2k1 - 2*M2*p2k2 + p1k1*p2k2 - p1k2*p2k2))*(2*k1k2 + k1k3 + k2k3 + 2*M2 - 4
     &*p1k1 - 2*p1k2 - p1k3 - 2*p1p2 + 2*p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(128.
     &*k1k2*p1k1*(k1k2 - p1k1 - p1k2)*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mttk3 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k3 + k2k3 + 2*M2 - 2*p1k1 - 2*p1k2 
     &- p1k3 - 2*p1p2 + p2k3)*(k1k2**2 + k1k2*(3*k1k3 - 3*p1k1 - 2*p1k2 - 4*p1p2 + p2k1) + 2*(k1k3
     &**2 - k2k3**2 + p1k1**2 + 2*p1k2**2 - 2*p1k1*p1k3 + 2*p1k2*p1k3 + 2*p1k1*p1p2 + 4*p1k2*p1p2 
     &- 2*M2*p2k1 - p1k1*p2k1 + p1k2*p2k1 + 2*p1k3*p2k1 + k1k3*(2*M2 - 2*p1k1 - 3*p1k2 - 2*p1k3 - 4
     &*p1p2 + p2k1) - 4*M2*p2k2 + p1k1*p2k2 + 2*p1k2*p2k2 + 4*p1k3*p2k2 - k2k3*(2*M2 - 3*p1k1 + p1k2 
     &- 2*p1k3 + 2*p1p2 + p2k2) + 2*p1k1*p2k3 - 2*p1k2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**5)/(128*k1k2*(k1k2 - p1k1 - p1k2)*p1k3*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mttk3 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2**2 + k1k2*(k1k3 - 3*p1k1 - 2
     &*p1k2 - 4*p1k3 - 4*p1p2 + p2k1) + 2*(p1k1**2 + 2*p1k2**2 + k1k3*(-2*M2 - p1k1 + p1k2) + k2k3*(
     &-4*M2 + p1k1 + 2*p1k2) + 2*p1k1*p1k3 + 4*p1k2*p1k3 + 2*p1k1*p1p2 + 4*p1k2*p1p2 - 2*M2*p2k1 - p1k1
     &*p2k1 + p1k2*p2k1 - 4*M2*p2k2 + p1k1*p2k2 + 2*p1k2*p2k2))*(2*k1k2 + k1k3 + k2k3 + 2*M2 - 2
     &*p1k1 - 4*p1k2 - p1k3 - 2*p1p2 + 2*p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(128.
     &*k1k2*(k1k2 - p1k1 - p1k2)*p1k2*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mttk3 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)**2*(k1k2**2 - 4*k1k2*(M2 + p1k1) + 4
     &*(p1k1**2 + p1k1*p1k2 + p1k2**2))*(k1k2 + M2 - p1k1 - p1k2 - p1p2 + p2k1 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**5)/(32*k1k2*p1k2*(-k1k2 + p1k1 + p1k2)**2))

!************************************************   2 Re(Mttk3 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(-1 + Nc**2)*(-k1k3 - k2k3 - 2*M2 + 2*p1k1 + 2
     &*p1k2 + p1k3 + 2*p1p2 - p2k3)*(k1k2**2 + k1k2*(3*k2k3 - 2*p1k1 - 3*p1k2 - 4*p1p2 + p2k2) + 2*(
     &-k1k3**2 + k2k3**2 + 2*p1k1**2 + p1k2**2 + 2*p1k1*p1k3 - 2*p1k2*p1k3 + 4*p1k1*p1p2 + 2*p1k2*p1p2 
     &- 4*M2*p2k1 + 2*p1k1*p2k1 + p1k2*p2k1 + 4*p1k3*p2k1 - k1k3*(2*M2 + p1k1 - 3*p1k2 - 2*p1k3 
     &+ 2*p1p2 + p2k1) - 2*M2*p2k2 + p1k1*p2k2 - p1k2*p2k2 + 2*p1k3*p2k2 + k2k3*(2*M2 - 3*p1k1 - 2
     &*p1k2 - 2*p1k3 - 4*p1p2 + p2k2) - 2*p1k1*p2k3 + 2*p1k2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**5)/(128*k1k2*(k1k2 - p1k1 - p1k2)*p1k3*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mbtk1 Mbtk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mbtk1 Mbtk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mbtk1 Msk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**3*Nc**2*(-1 + Nc**2)*(k1k3**2*(4*M2 + 5*p1k2 - 4
     &*p1p2 - 5*p2k2) + k1k2**2*(4*M2 + 5*p1k3 - 4*p1p2 - 5*p2k3) - k1k3*(14*p1k1*p1k2 - 5*p1k2**2 
     &+ 4*p1k1*p1k3 + 5*p1k2*p1k3 + 14*p1k2*p1p2 + 4*p1k3*p1p2 - 28*p1k2*p2k1 - 8*p1k3*p2k1 - 8*k2k3
     &*(2*M2 - p1k1 - 2*p1p2 + p2k1) + 14*p1k1*p2k2 + 10*p1k2*p2k2 - 15*p1k3*p2k2 + 14*p1p2*p2k2 
     &- 5*p2k2**2 + 2*M2*(7*p1k2 + 2*p1k3 - 21*p2k2 - 6*p2k3) + 4*p1k1*p2k3 + 5*p1k2*p2k3 + 4*p1p2
     &*p2k3 + 5*p2k2*p2k3) + k1k2*(-4*M2*p1k2 - 4*p1k1*p1k2 - 14*M2*p1k3 - 14*p1k1*p1k3 - 5*p1k2*p1k3 
     &+ 5*p1k3**2 - 4*p1k2*p1p2 - 14*p1k3*p1p2 + 8*p1k2*p2k1 + 28*p1k3*p2k1 + 8*k2k3*(2*M2 - p1k1 
     &- 2*p1p2 + p2k1) + 12*M2*p2k2 - 4*p1k1*p2k2 - 5*p1k3*p2k2 - 4*p1p2*p2k2 + 42*M2*p2k3 - 14*p1k1
     &*p2k3 + 15*p1k2*p2k3 - 10*p1k3*p2k3 - 14*p1p2*p2k3 - 5*p2k2*p2k3 + 5*p2k3**2 + k1k3*(28*M2 
     &- 5*p1k2 - 5*p1k3 - 28*p1p2 + 5*p2k2 + 5*p2k3)) - 2*(-2*p1k1*p1k2**2 - 14*p1k1*p1k2*p1k3 - 2
     &*p1k1*p1k3**2 - 4*p1k2**2*p1p2 - 28*p1k2*p1k3*p1p2 - 4*p1k3**2*p1p2 + 4*p1k2**2*p2k1 + 28*p1k2
     &*p1k3*p2k1 + 4*p1k3**2*p2k1 + 4*M2*p1k2*p2k2 + 14*M2*p1k3*p2k2 - 5*p1k2*p1k3*p2k2 + 5*p1k3**2
     &*p2k2 + 4*p1k2*p1p2*p2k2 + 14*p1k3*p1p2*p2k2 - 4*p1k2*p2k1*p2k2 - 14*p1k3*p2k1*p2k2 - 4*M2*p2k2
     &**2 + 2*p1k1*p2k2**2 + 5*p1k3*p2k2**2 + 14*M2*p1k2*p2k3 + 5*p1k2**2*p2k3 + 4*M2*p1k3*p2k3 
     &- 5*p1k2*p1k3*p2k3 + 14*p1k2*p1p2*p2k3 + 4*p1k3*p1p2*p2k3 - 14*p1k2*p2k1*p2k3 - 4*p1k3*p2k1*p2k3 
     &- 28*M2*p2k2*p2k3 + 14*p1k1*p2k2*p2k3 - 5*p1k2*p2k2*p2k3 - 5*p1k3*p2k2*p2k3 - 4*M2*p2k3**2 
     &+ 2*p1k1*p2k3**2 + 5*p1k2*p2k3**2 + 4*k2k3*(4*M2**2 - 2*p1k1**2 + p1k1*(p1k2 + p1k3 - 2*p1p2 
     &+ 2*p2k1 - 3*p2k2 - 3*p2k3) - (2*p1p2 - p2k1)*(p1k2 + p1k3 + 2*p1p2 + p2k2 + p2k3) + 2*M2*(-p1k1 
     &+ p1k2 + p1k3 + p2k1 + p2k2 + p2k3))))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(64*k2k3**2*p1k1*(M2 + p1p2)))

!************************************************   2 Re(Mbtk1 Msk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**3*Nc**2*(-1 + Nc**2)*(-8*k2k3*M2**2 + 12*k2k3
     &*M2*p1k1 - 4*k2k3*p1k1**2 - 10*k2k3*M2*p1k2 + k2k3*p1k1*p1k2 + 20*M2*p1k1*p1k2 + 12*k2k3*M2*p1k3 
     &- 6*k2k3*p1k1*p1k3 + 16*M2*p1k1*p1k3 + 6*p1k1**2*p1k3 + 16*M2*p1k2*p1k3 - 24*p1k1*p1k2*p1k3 
     &+ 6*p1k2**2*p1k3 + 20*M2*p1k3**2 + 12*p1k1*p1k3**2 - 24*p1k2*p1k3**2 + 6*p1k3**3 - 4*k2k3*p1k1
     &*p1p2 + 10*k2k3*p1k2*p1p2 - 12*k2k3*p1k3*p1p2 - 12*p1k1*p1k3*p1p2 - 12*p1k2*p1k3*p1p2 - 12
     &*p1k3**2*p1p2 + 8*k2k3*p1p2**2 - 4*k2k3*M2*p2k1 + 4*k2k3*p1k1*p2k1 - 5*k2k3*p1k2*p2k1 - 8*M2
     &*p1k2*p2k1 + 10*p1k1*p1k2*p2k1 + 6*k2k3*p1k3*p2k1 + 8*M2*p1k3*p2k1 + 2*p1k1*p1k3*p2k1 + 24*p1k2
     &*p1k3*p2k1 + 6*p1k3**2*p2k1 - 4*k2k3*p1p2*p2k1 - 12*p1k2*p1p2*p2k1 - 12*p1k3*p1p2*p2k1 + 6*p1k2
     &*p2k1**2 + 12*p1k3*p2k1**2 - 6*k2k3*M2*p2k2 + 7*k2k3*p1k1*p2k2 + 8*M2*p1k1*p2k2 - 10*p1k1
     &**2*p2k2 + 16*M2*p1k3*p2k2 + 18*p1k1*p1k3*p2k2 - 6*p1k2*p1k3*p2k2 + 28*p1k3**2*p2k2 + 6*k2k3*p1p2
     &*p2k2 - 28*p1k1*p1p2*p2k2 - 20*p1k3*p1p2*p2k2 - 3*k2k3*p2k1*p2k2 + 20*M2*p2k1*p2k2 - 6*p1k1
     &*p2k1*p2k2 - 22*p1k3*p2k1*p2k2 + k1k3**2*(-14*M2 + 8*p1k1 - 7*p1k2 + 3*p1k3 + 14*p1p2 - 8*p2k1 
     &+ 7*p2k2 - 3*p2k3) + k1k2**2*(6*M2 - p1k1 + p1k3 - 6*p1p2 + p2k1 - p2k3) + 4*M2*p1k1*p2k3 
     &- 14*p1k1**2*p2k3 - 4*M2*p1k2*p2k3 + 6*p1k1*p1k2*p2k3 - 6*p1k2**2*p2k3 + 12*M2*p1k3*p2k3 - 30
     &*p1k1*p1k3*p2k3 + 20*p1k2*p1k3*p2k3 - 12*p1k3**2*p2k3 - 8*p1k1*p1p2*p2k3 - 28*p1k3*p1p2*p2k3 
     &+ 4*M2*p2k1*p2k3 - 6*p1k1*p2k1*p2k3 + 14*p1k2*p2k1*p2k3 + 14*p1k3*p2k1*p2k3 + 4*M2*p2k2*p2k3 
     &- 16*p1k1*p2k2*p2k3 + 6*p1k2*p2k2*p2k3 - 28*p1k3*p2k2*p2k3 + 8*M2*p2k3**2 - 2*p1k1*p2k3**2 + 4
     &*p1k2*p2k3**2 + 6*p1k3*p2k3**2 + k1k3*(-8*M2**2 + 21*p1k1*p1k2 - 4*p1k2**2 - 2*p1k1*p1k3 + 24
     &*p1k2*p1k3 - 9*p1k3**2 + 18*p1k1*p1p2 + 14*p1k2*p1p2 - 8*p1k3*p1p2 + 8*p1p2**2 - 33*p1k2*p2k1 
     &+ 6*p1k3*p2k1 - 10*p1p2*p2k1 - 3*k2k3*(2*M2 - p1k1 - 2*p1p2 + p2k1) + 3*p1k1*p2k2 + 4*p1k2
     &*p2k2 - 28*p1k3*p2k2 + 14*p1p2*p2k2 + 9*p2k1*p2k2 - 2*p1k1*p2k3 - 10*p1k2*p2k3 + 12*p1k3*p2k3 
     &+ 4*p1p2*p2k3 - 2*p2k1*p2k3 + 14*p2k2*p2k3 - 3*p2k3**2 + 2*M2*(p1k1 - 9*p1k2 - 5*p2k1 - 5*p2k2 
     &+ 2*p2k3)) + k1k2*(-16*M2**2 + 4*M2*p1k1 + 6*p1k1**2 - 2*M2*p1k2 - 11*p1k1*p1k2 + 8*M2*p1k3 
     &+ 19*p1k1*p1k3 - 5*p1k2*p1k3 + 12*p1k3**2 + 2*p1k2*p1p2 - 8*p1k3*p1p2 + 16*p1p2**2 + 4*M2*p2k1 
     &- 6*p1k1*p2k1 + 3*p1k2*p2k1 + 9*p1k3*p2k1 - 8*p1p2*p2k1 + 3*k2k3*(2*M2 - p1k1 - 2*p1p2 + p2k1) 
     &- 6*M2*p2k2 + 11*p1k1*p2k2 + 3*p1k3*p2k2 + 6*p1p2*p2k2 - 3*p2k1*p2k2 + 12*M2*p2k3 - 27*p1k1
     &*p2k3 + 5*p1k2*p2k3 - 10*p1k3*p2k3 - 12*p1p2*p2k3 - p2k1*p2k3 - 3*p2k2*p2k3 - 2*p2k3**2 
     &+ k1k3*(-4*M2 - 4*p1k1 + 2*p1k2 - 5*p1k3 + 4*p1p2 + 4*p2k1 - 2*p2k2 + 5*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)
     &/(128*k1k3*k2k3*p1k1*(M2 + p1p2)))

!************************************************   2 Re(Mbtk1 Msk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**3*Nc**2*(-1 + Nc**2)*(-8*k2k3*M2**2 + 12*k2k3
     &*M2*p1k1 - 4*k2k3*p1k1**2 + 12*k2k3*M2*p1k2 - 6*k2k3*p1k1*p1k2 + 16*M2*p1k1*p1k2 + 6*p1k1**2*p1k2 
     &+ 20*M2*p1k2**2 + 12*p1k1*p1k2**2 + 6*p1k2**3 - 10*k2k3*M2*p1k3 + k2k3*p1k1*p1k3 + 20*M2
     &*p1k1*p1k3 + 16*M2*p1k2*p1k3 - 24*p1k1*p1k2*p1k3 - 24*p1k2**2*p1k3 + 6*p1k2*p1k3**2 - 4*k2k3*p1k1
     &*p1p2 - 12*k2k3*p1k2*p1p2 - 12*p1k1*p1k2*p1p2 - 12*p1k2**2*p1p2 + 10*k2k3*p1k3*p1p2 - 12*p1k2
     &*p1k3*p1p2 + 8*k2k3*p1p2**2 - 4*k2k3*M2*p2k1 + 4*k2k3*p1k1*p2k1 + 6*k2k3*p1k2*p2k1 + 8*M2
     &*p1k2*p2k1 + 2*p1k1*p1k2*p2k1 + 6*p1k2**2*p2k1 - 5*k2k3*p1k3*p2k1 - 8*M2*p1k3*p2k1 + 10*p1k1*p1k3
     &*p2k1 + 24*p1k2*p1k3*p2k1 - 4*k2k3*p1p2*p2k1 - 12*p1k2*p1p2*p2k1 - 12*p1k3*p1p2*p2k1 + 12
     &*p1k2*p2k1**2 + 6*p1k3*p2k1**2 + k1k3**2*(6*M2 - p1k1 + p1k2 - 6*p1p2 + p2k1 - p2k2) + 4*M2*p1k1
     &*p2k2 - 14*p1k1**2*p2k2 + 12*M2*p1k2*p2k2 - 30*p1k1*p1k2*p2k2 - 12*p1k2**2*p2k2 - 4*M2*p1k3
     &*p2k2 + 6*p1k1*p1k3*p2k2 + 20*p1k2*p1k3*p2k2 - 6*p1k3**2*p2k2 - 8*p1k1*p1p2*p2k2 - 28*p1k2*p1p2
     &*p2k2 + 4*M2*p2k1*p2k2 - 6*p1k1*p2k1*p2k2 + 14*p1k2*p2k1*p2k2 + 14*p1k3*p2k1*p2k2 + 8*M2*p2k2
     &**2 - 2*p1k1*p2k2**2 + 6*p1k2*p2k2**2 + 4*p1k3*p2k2**2 - 6*k2k3*M2*p2k3 + 7*k2k3*p1k1*p2k3 
     &+ 8*M2*p1k1*p2k3 - 10*p1k1**2*p2k3 + 16*M2*p1k2*p2k3 + 18*p1k1*p1k2*p2k3 + 28*p1k2**2*p2k3 - 6
     &*p1k2*p1k3*p2k3 + 6*k2k3*p1p2*p2k3 - 28*p1k1*p1p2*p2k3 - 20*p1k2*p1p2*p2k3 - 3*k2k3*p2k1*p2k3 
     &+ 20*M2*p2k1*p2k3 - 6*p1k1*p2k1*p2k3 - 22*p1k2*p2k1*p2k3 + 4*M2*p2k2*p2k3 - 16*p1k1*p2k2*p2k3 
     &- 28*p1k2*p2k2*p2k3 + 6*p1k3*p2k2*p2k3 + k1k2**2*(-14*M2 + 8*p1k1 + 3*p1k2 - 7*p1k3 + 14*p1p2 
     &- 8*p2k1 - 3*p2k2 + 7*p2k3) + k1k3*(-16*M2**2 + 6*p1k1**2 + 19*p1k1*p1k2 + 12*p1k2**2 - 11
     &*p1k1*p1k3 - 5*p1k2*p1k3 - 8*p1k2*p1p2 + 2*p1k3*p1p2 + 16*p1p2**2 - 6*p1k1*p2k1 + 9*p1k2*p2k1 
     &+ 3*p1k3*p2k1 - 8*p1p2*p2k1 + 3*k2k3*(2*M2 - p1k1 - 2*p1p2 + p2k1) - 27*p1k1*p2k2 - 10*p1k2
     &*p2k2 + 5*p1k3*p2k2 - 12*p1p2*p2k2 - p2k1*p2k2 - 2*p2k2**2 + 2*M2*(2*p1k1 + 4*p1k2 - p1k3 + 2
     &*p2k1 + 6*p2k2 - 3*p2k3) + 11*p1k1*p2k3 + 3*p1k2*p2k3 + 6*p1p2*p2k3 - 3*p2k1*p2k3 - 3*p2k2*p2k3) 
     &+ k1k2*(-8*M2**2 + 2*M2*p1k1 - 2*p1k1*p1k2 - 9*p1k2**2 - 18*M2*p1k3 + 21*p1k1*p1k3 + 24*p1k2
     &*p1k3 - 4*p1k3**2 + 18*p1k1*p1p2 - 8*p1k2*p1p2 + 14*p1k3*p1p2 + 8*p1p2**2 - 10*M2*p2k1 + 6
     &*p1k2*p2k1 - 33*p1k3*p2k1 - 10*p1p2*p2k1 - 3*k2k3*(2*M2 - p1k1 - 2*p1p2 + p2k1) + 4*M2*p2k2 
     &- 2*p1k1*p2k2 + 12*p1k2*p2k2 - 10*p1k3*p2k2 + 4*p1p2*p2k2 - 2*p2k1*p2k2 - 3*p2k2**2 + k1k3*(-4
     &*M2 - 4*p1k1 - 5*p1k2 + 2*p1k3 + 4*p1p2 + 4*p2k1 + 5*p2k2 - 2*p2k3) - 10*M2*p2k3 + 3*p1k1*p2k3 
     &- 28*p1k2*p2k3 + 4*p1k3*p2k3 + 14*p1p2*p2k3 + 9*p2k1*p2k3 + 14*p2k2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)
     &/(128*k1k2*k2k3*p1k1*(M2 + p1p2)))

!************************************************   2 Re(Mbtk1 Mvk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*(2 - 3*Nc**2 + Nc**4)*(2*k1k2**2 + 2*k1k3**2 
     &+ k1k3*(-21*p1k2 - 6*p1k3 + 7*p2k2 + 2*p2k3) + k1k2*(14*k1k3 - 6*p1k2 - 21*p1k3 + 2*p2k2 + 7
     &*p2k3) - 2*(-2*p1k2**2 + 4*k2k3*(2*M2 - 3*p1k1 - 2*p1p2 + p2k1) + p1k3*(-2*p1k3 + 7*p2k2 + 2*p2k3) 
     &+ p1k2*(-14*p1k3 + 2*p2k2 + 7*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32.
     &*k2k3**2*p1k1))

!************************************************   2 Re(Mbtk1 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(-5*k1k3**2 + 2*k2k3*M2 - 3
     &*k2k3*p1k1 - 4*p1k1*p1k2 - 14*p1k1*p1k3 + 10*p1k2*p1k3 - 10*p1k3**2 - 2*k2k3*p1p2 + k2k3*p2k1 
     &+ 8*p1k2*p2k1 + 10*p1k3*p2k1 - 4*p1k1*p2k2 - 8*p1k3*p2k2 + k1k3*(8*M2 - 3*p1k1 - 6*p1k2 + 15
     &*p1k3 - 8*p1p2 - p2k1 + 4*p2k2 - 5*p2k3) + 4*p1k1*p2k3 - 2*p1k2*p2k3 + 10*p1k3*p2k3 + k1k2*(5
     &*k1k3 + 10*M2 - 15*p1k1 - 9*p1k3 - 10*p1p2 + p2k1 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k3*k2k3*p1k1))

!************************************************   2 Re(Mbtk1 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*Nc**2*(-1 + Nc**2)*(5*k1k2**2 - 2*k2k3*M2 + 3
     &*k2k3*p1k1 + 14*p1k1*p1k2 + 10*p1k2**2 + 4*p1k1*p1k3 - 10*p1k2*p1k3 + 2*k2k3*p1p2 - k2k3*p2k1 
     &- 10*p1k2*p2k1 - 8*p1k3*p2k1 - 4*p1k1*p2k2 - 10*p1k2*p2k2 + 2*p1k3*p2k2 - k1k3*(10*M2 - 15
     &*p1k1 - 9*p1k2 - 10*p1p2 + p2k1 + p2k2) + k1k2*(-5*k1k3 - 8*M2 + 3*p1k1 - 15*p1k2 + 6*p1k3 + 8
     &*p1p2 + p2k1 + 5*p2k2 - 4*p2k3) + 4*p1k1*p2k3 + 8*p1k2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k2*k2k3*p1k1))

!************************************************   2 Re(Mbtk1 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)**2*(M2 - p1k1)*(3*k1k2**2 + 2*k2k3
     &*M2 - k2k3*p1k2 - 2*M2*p1k2 + 3*p1k2**2 - 4*M2*p1k3 - 2*k2k3*p1p2 + 2*p1k2*p1p2 + 4*p1k3*p1p2 
     &+ 2*k2k3*p2k1 - 2*p1k2*p2k1 - 4*p1k3*p2k1 + k2k3*p2k2 - 2*M2*p2k2 + 4*p1k1*p2k2 - 2*p1k2*p2k2 
     &- 2*p1k3*p2k2 + 2*p1p2*p2k2 - 2*p2k1*p2k2 - p2k2**2 + 2*k1k3*(2*M2 - 2*p1p2 + 2*p2k1 + p2k2) 
     &+ k1k2*(k2k3 + 2*(M2 - 3*p1k2 - p1p2 + p2k1 + p2k2 - 2*p2k3)) - 4*M2*p2k3 + 8*p1k1*p2k3 + 4
     &*p1k2*p2k3 + 4*p1p2*p2k3 - 4*p2k1*p2k3 - 2*p2k2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(32*k2k3*p1k1**2*(-k1k2 + p1k1 + p1k2)))

!************************************************   2 Re(Mbtk1 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)**2*(M2 - p1k1)*(3*k1k3**2 + 2*k2k3
     &*M2 - 4*M2*p1k2 - k2k3*p1k3 - 2*M2*p1k3 + 3*p1k3**2 - 2*k2k3*p1p2 + 4*p1k2*p1p2 + 2*p1k3*p1p2 
     &+ 2*k2k3*p2k1 - 4*p1k2*p2k1 - 2*p1k3*p2k1 - 4*M2*p2k2 + 8*p1k1*p2k2 + 4*p1k3*p2k2 + 4*p1p2*p2k2 
     &- 4*p2k1*p2k2 + k2k3*p2k3 - 2*M2*p2k3 + 4*p1k1*p2k3 - 2*p1k2*p2k3 - 2*p1k3*p2k3 + 2*p1p2*p2k3 
     &- 2*p2k1*p2k3 - 2*p2k2*p2k3 - p2k3**2 + 2*k1k2*(2*M2 - 2*p1p2 + 2*p2k1 + p2k3) + k1k3*(k2k3 
     &+ 2*(M2 - 3*p1k3 - p1p2 + p2k1 - 2*p2k2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(32*k2k3*p1k1**2*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mbtk1 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(-1 + Nc**2)*(k1k2 + k1k3 + 2*M2 - p1k1 - 2*p1k2 
     &- 2*p1k3 - 2*p1p2 + p2k1)*(4*k2k3*M2 - 4*k2k3*p1k1 + k2k3*p1k2 - 4*p1k1*p1k2 - 2*p1k2**2 - k1k3
     &*(k2k3 + 4*M2 - 8*p1k2 - 6*p1k3) + k2k3*p1k3 + 4*p1k1*p1k3 - 4*p1k2*p1k3 - 6*p1k3**2 - k1k2
     &*(k2k3 - 4*M2 - 2*p1k2 + 4*p1k3) - 4*k2k3*p1p2 - 4*p1k2*p1p2 + 4*p1k3*p1p2 - k2k3*p2k2 + 4
     &*M2*p2k2 + 2*p1k2*p2k2 - 4*p1k3*p2k2 - k2k3*p2k3 - 4*M2*p2k3 + 8*p1k2*p2k3 + 6*p1k3*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**5)/(128*k2k3*p1k1*(k2k3 - p1k2 - p1k3)*p1k3))

!************************************************   2 Re(Mbtk1 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2 + k1k3 + 2*M2 - p1k1 - 2*p1k2 
     &- 2*p1k3 - 2*p1p2 + p2k1)*(-4*k2k3*M2 + 4*k2k3*p1k1 - k2k3*p1k2 - 4*p1k1*p1k2 + 6*p1k2**2 + k1k2
     &*(k2k3 + 4*M2 - 6*p1k2 - 8*p1k3) + k1k3*(k2k3 - 4*M2 + 4*p1k2 - 2*p1k3) - k2k3*p1k3 + 4*p1k1
     &*p1k3 + 4*p1k2*p1k3 + 2*p1k3**2 + 4*k2k3*p1p2 - 4*p1k2*p1p2 + 4*p1k3*p1p2 + k2k3*p2k2 + 4
     &*M2*p2k2 - 6*p1k2*p2k2 - 8*p1k3*p2k2 + k2k3*p2k3 - 4*M2*p2k3 + 4*p1k2*p2k3 - 2*p1k3*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**5)/(128*k2k3*p1k1*p1k2*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(Mbtk1 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2 + 2*M2 - 2*p1k1 - 2*p1k2)*(k1k2
     &**2 - 2*k2k3*M2 + 2*k2k3*p1k1 + k2k3*p1k2 + 2*M2*p1k2 - 2*p1k1*p1k2 - 3*p1k2**2 + 4*M2*p1k3 
     &- 4*p1k1*p1k3 + 2*k2k3*p1p2 - 2*p1k2*p1p2 - 4*p1k3*p1p2 - 4*k2k3*p2k1 + 4*p1k2*p2k1 + 8*p1k3
     &*p2k1 - k2k3*p2k2 + 2*M2*p2k2 - 2*p1k1*p2k2 + 2*p1k2*p2k2 + 2*p1k3*p2k2 - 2*p1p2*p2k2 + p2k2
     &**2 + 2*k1k3*(2*M2 - 2*p1k1 - 2*p1k2 - 2*p1p2 + p2k2) + 4*M2*p2k3 - 4*p1k1*p2k3 - 4*p1k2*p2k3 
     &- 4*p1p2*p2k3 + 2*p2k2*p2k3 + k1k2*(2*k1k3 - k2k3 + 2*(M2 - p1k1 + p1k2 + p1k3 - p1p2 + p2k2 
     &+ p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(64*k2k3*p1k1*p1k2*(-k1k2 + p1k1 + p1k2)))

!************************************************   2 Re(Mbtk1 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k3 + 2*M2 - 2*p1k1 - 2*p1k3)*(k1k3
     &**2 - k1k3*k2k3 - 2*k2k3*M2 + 2*k2k3*p1k1 + 4*M2*p1k2 - 4*p1k1*p1k2 + k2k3*p1k3 + 2*M2*p1k3 
     &- 2*p1k1*p1k3 - 3*p1k3**2 + 2*k2k3*p1p2 - 4*p1k2*p1p2 - 2*p1k3*p1p2 - 4*k2k3*p2k1 + 8*p1k2*p2k1 
     &+ 4*p1k3*p2k1 + 4*M2*p2k2 - 4*p1k1*p2k2 - 4*p1k3*p2k2 - 4*p1p2*p2k2 - k2k3*p2k3 + 2*M2*p2k3 
     &- 2*p1k1*p2k3 + 2*p1k2*p2k3 + 2*p1k3*p2k3 - 2*p1p2*p2k3 + 2*p2k2*p2k3 + p2k3**2 + 2*k1k2*(k1k3 
     &+ 2*M2 - 2*p1k1 - 2*p1k3 - 2*p1p2 + p2k3) + 2*k1k3*(M2 - p1k1 + p1k2 + p1k3 - p1p2 + p2k2 
     &+ p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(64*k2k3*p1k1*p1k3*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mbtk2 Mbtk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Mbtk2 Msk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**3*Nc**2*(-1 + Nc**2)*(-14*k2k3**2*M2 - 8*k2k3
     &*M2**2 - 7*k2k3**2*p1k1 - 18*k2k3*M2*p1k1 - 4*k2k3*p1k1**2 + 8*k2k3**2*p1k2 + 2*k2k3*M2*p1k2 
     &+ 21*k2k3*p1k1*p1k2 + 20*M2*p1k1*p1k2 + 3*k2k3**2*p1k3 + 24*k2k3*p1k1*p1k3 + 16*M2*p1k1*p1k3 
     &+ 6*p1k1**2*p1k3 - 2*k2k3*p1k2*p1k3 + 16*M2*p1k2*p1k3 - 24*p1k1*p1k2*p1k3 + 6*p1k2**2*p1k3 - 9
     &*k2k3*p1k3**2 + 20*M2*p1k3**2 - 24*p1k1*p1k3**2 + 12*p1k2*p1k3**2 + 6*p1k3**3 + 14*k2k3**2*p1p2 
     &+ 14*k2k3*p1k1*p1p2 + 18*k2k3*p1k2*p1p2 - 8*k2k3*p1k3*p1p2 - 12*p1k1*p1k3*p1p2 - 12*p1k2*p1k3
     &*p1p2 - 12*p1k3**2*p1p2 + 8*k2k3*p1p2**2 + 7*k2k3**2*p2k1 - 10*k2k3*M2*p2k1 + 4*k2k3*p1k1
     &*p2k1 + 3*k2k3*p1k2*p2k1 + 8*M2*p1k2*p2k1 - 10*p1k2**2*p2k1 - 28*k2k3*p1k3*p2k1 + 16*M2*p1k3*p2k1 
     &- 6*p1k1*p1k3*p2k1 + 18*p1k2*p1k3*p2k1 + 28*p1k3**2*p2k1 + 14*k2k3*p1p2*p2k1 - 28*p1k2*p1p2
     &*p2k1 - 20*p1k3*p1p2*p2k1 - 8*k2k3**2*p2k2 - 10*k2k3*M2*p2k2 - 33*k2k3*p1k1*p2k2 - 8*M2*p1k1
     &*p2k2 + 10*p1k1*p1k2*p2k2 + 6*k2k3*p1k3*p2k2 + 8*M2*p1k3*p2k2 + 24*p1k1*p1k3*p2k2 + 2*p1k2*p1k3
     &*p2k2 + 6*p1k3**2*p2k2 - 10*k2k3*p1p2*p2k2 - 12*p1k1*p1p2*p2k2 - 12*p1k3*p1p2*p2k2 + 9*k2k3
     &*p2k1*p2k2 + 20*M2*p2k1*p2k2 - 6*p1k2*p2k1*p2k2 - 22*p1k3*p2k1*p2k2 + 6*p1k1*p2k2**2 + 12*p1k3
     &*p2k2**2 + k1k3*(-8*M2**2 + p1k1*p1k2 - 4*p1k2**2 - 6*p1k2*p1k3 + 10*p1k1*p1p2 - 4*p1k2*p1p2 
     &- 12*p1k3*p1p2 + 8*p1p2**2 + 7*p1k2*p2k1 + 6*p1p2*p2k1 - 5*p1k1*p2k2 + 4*p1k2*p2k2 + 6*p1k3
     &*p2k2 - 4*p1p2*p2k2 - 3*p2k1*p2k2 - 3*k2k3*(2*M2 - p1k2 - 2*p1p2 + p2k2) - 2*M2*(5*p1k1 - 6*p1k2 
     &- 6*p1k3 + 3*p2k1 + 2*p2k2)) + k1k2**2*(6*M2 - p1k2 + p1k3 - 6*p1p2 + p2k2 - p2k3) - 3*k2k3
     &**2*p2k3 + 4*k2k3*M2*p2k3 - 10*k2k3*p1k1*p2k3 - 4*M2*p1k1*p2k3 - 6*p1k1**2*p2k3 - 2*k2k3*p1k2
     &*p2k3 + 4*M2*p1k2*p2k3 + 6*p1k1*p1k2*p2k3 - 14*p1k2**2*p2k3 + 12*k2k3*p1k3*p2k3 + 12*M2*p1k3
     &*p2k3 + 20*p1k1*p1k3*p2k3 - 30*p1k2*p1k3*p2k3 - 12*p1k3**2*p2k3 + 4*k2k3*p1p2*p2k3 - 8*p1k2
     &*p1p2*p2k3 - 28*p1k3*p1p2*p2k3 + 14*k2k3*p2k1*p2k3 + 4*M2*p2k1*p2k3 + 6*p1k1*p2k1*p2k3 - 16*p1k2
     &*p2k1*p2k3 - 28*p1k3*p2k1*p2k3 - 2*k2k3*p2k2*p2k3 + 4*M2*p2k2*p2k3 + 14*p1k1*p2k2*p2k3 - 6
     &*p1k2*p2k2*p2k3 + 14*p1k3*p2k2*p2k3 - 3*k2k3*p2k3**2 + 8*M2*p2k3**2 + 4*p1k1*p2k3**2 - 2*p1k2
     &*p2k3**2 + 6*p1k3*p2k3**2 + k1k2*(-16*M2**2 - 2*M2*p1k1 + 4*M2*p1k2 - 11*p1k1*p1k2 + 6*p1k2**2 
     &+ 8*M2*p1k3 - 5*p1k1*p1k3 + 19*p1k2*p1k3 + 12*p1k3**2 + 2*p1k1*p1p2 - 8*p1k3*p1p2 + 16*p1p2
     &**2 - 6*M2*p2k1 + 11*p1k2*p2k1 + 3*p1k3*p2k1 + 6*p1p2*p2k1 + 4*M2*p2k2 + 3*p1k1*p2k2 - 6*p1k2
     &*p2k2 + 9*p1k3*p2k2 - 8*p1p2*p2k2 - 3*p2k1*p2k2 + 3*k1k3*(2*M2 - p1k2 - 2*p1p2 + p2k2) + 12*M2
     &*p2k3 + 5*p1k1*p2k3 - 27*p1k2*p2k3 - 10*p1k3*p2k3 - 12*p1p2*p2k3 - 3*p2k1*p2k3 - p2k2*p2k3 
     &- 2*p2k3**2 + k2k3*(-4*M2 + 2*p1k1 - 4*p1k2 - 5*p1k3 + 4*p1p2 - 2*p2k1 + 4*p2k2 + 5*p2k3)))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k3*k2k3*p1k2*(M2 + p1p2)))

!************************************************   2 Re(Mbtk2 Msk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**3*Nc**2*(-1 + Nc**2)*(4*k2k3**2*M2 + 5*k2k3**2
     &*p1k1 - 14*k2k3*M2*p1k1 + 5*k2k3*p1k1**2 - 14*k2k3*p1k1*p1k2 + 4*p1k1**2*p1k2 - 4*k2k3*M2*p1k3 
     &- 5*k2k3*p1k1*p1k3 - 4*k2k3*p1k2*p1k3 + 28*p1k1*p1k2*p1k3 + 4*p1k2*p1k3**2 - 4*k2k3**2*p1p2 
     &- 14*k2k3*p1k1*p1p2 + 8*p1k1**2*p1p2 - 4*k2k3*p1k3*p1p2 + 56*p1k1*p1k3*p1p2 + 8*p1k3**2*p1p2 
     &- 5*k2k3**2*p2k1 + 42*k2k3*M2*p2k1 - 10*k2k3*p1k1*p2k1 - 8*M2*p1k1*p2k1 - 14*k2k3*p1k2*p2k1 
     &+ 15*k2k3*p1k3*p2k1 - 28*M2*p1k3*p2k1 + 10*p1k1*p1k3*p2k1 - 10*p1k3**2*p2k1 - 14*k2k3*p1p2*p2k1 
     &- 8*p1k1*p1p2*p2k1 - 28*p1k3*p1p2*p2k1 + 5*k2k3*p2k1**2 + 8*M2*p2k1**2 - 4*p1k2*p2k1**2 - 10
     &*p1k3*p2k1**2 + 28*k2k3*p1k1*p2k2 - 8*p1k1**2*p2k2 + 8*k2k3*p1k3*p2k2 - 56*p1k1*p1k3*p2k2 - 8
     &*p1k3**2*p2k2 + 8*p1k1*p2k1*p2k2 + 28*p1k3*p2k1*p2k2 + k1k2**2*(4*M2 + 5*p1k3 - 4*p1p2 - 5*p2k3) 
     &+ 12*k2k3*M2*p2k3 - 5*k2k3*p1k1*p2k3 - 28*M2*p1k1*p2k3 - 10*p1k1**2*p2k3 - 4*k2k3*p1k2*p2k3 
     &- 8*M2*p1k3*p2k3 + 10*p1k1*p1k3*p2k3 - 4*k2k3*p1p2*p2k3 - 28*p1k1*p1p2*p2k3 - 8*p1k3*p1p2
     &*p2k3 - 5*k2k3*p2k1*p2k3 + 56*M2*p2k1*p2k3 + 10*p1k1*p2k1*p2k3 - 28*p1k2*p2k1*p2k3 + 10*p1k3
     &*p2k1*p2k3 + 28*p1k1*p2k2*p2k3 + 8*p1k3*p2k2*p2k3 + 8*M2*p2k3**2 - 10*p1k1*p2k3**2 - 4*p1k2*p2k3
     &**2 - 8*k1k3*(4*M2**2 + p1k1*p1k2 - 2*p1k2**2 + p1k2*p1k3 - 2*p1k1*p1p2 - 2*p1k2*p1p2 - 2*p1k3
     &*p1p2 - 4*p1p2**2 - 3*p1k2*p2k1 - 2*p1p2*p2k1 + k2k3*(-2*M2 + p1k2 + 2*p1p2 - p2k2) + p1k1
     &*p2k2 + 2*p1k2*p2k2 + p1k3*p2k2 + 2*p1p2*p2k2 + p2k1*p2k2 - 3*p1k2*p2k3 - 2*p1p2*p2k3 + p2k2
     &*p2k3 + 2*M2*(p1k1 - p1k2 + p1k3 + p2k1 + p2k2 + p2k3)) + k1k2*(-4*M2*p1k1 - 4*p1k1*p1k2 - 14
     &*M2*p1k3 - 5*p1k1*p1k3 - 14*p1k2*p1k3 + 5*p1k3**2 - 4*p1k1*p1p2 - 14*p1k3*p1p2 + 12*M2*p2k1 - 4
     &*p1k2*p2k1 - 5*p1k3*p2k1 - 4*p1p2*p2k1 + 8*p1k1*p2k2 + 28*p1k3*p2k2 + 8*k1k3*(2*M2 - p1k2 - 2
     &*p1p2 + p2k2) + 42*M2*p2k3 + 15*p1k1*p2k3 - 14*p1k2*p2k3 - 10*p1k3*p2k3 - 14*p1p2*p2k3 - 5*p2k1
     &*p2k3 + 5*p2k3**2 + k2k3*(28*M2 - 5*p1k1 - 5*p1k3 - 28*p1p2 + 5*p2k1 + 5*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)
     &/(64*k1k3**2*p1k2*(M2 + p1p2)))

!************************************************   2 Re(Mbtk2 Msk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**3*Nc**2*(-1 + Nc**2)*(6*k2k3**2*M2 - 16*k2k3*M2
     &**2 + k2k3**2*p1k1 + 8*k2k3*M2*p1k1 + 12*k2k3*p1k1**2 + 20*M2*p1k1**2 + 6*p1k1**3 - k2k3**2
     &*p1k2 + 4*k2k3*M2*p1k2 + 19*k2k3*p1k1*p1k2 + 16*M2*p1k1*p1k2 + 12*p1k1**2*p1k2 + 6*k2k3*p1k2*
     &*2 + 6*p1k1*p1k2**2 - 2*k2k3*M2*p1k3 - 5*k2k3*p1k1*p1k3 + 16*M2*p1k1*p1k3 - 24*p1k1**2*p1k3 - 11
     &*k2k3*p1k2*p1k3 + 20*M2*p1k2*p1k3 - 24*p1k1*p1k2*p1k3 + 6*p1k1*p1k3**2 - 6*k2k3**2*p1p2 - 8
     &*k2k3*p1k1*p1p2 - 12*p1k1**2*p1p2 - 12*p1k1*p1k2*p1p2 + 2*k2k3*p1k3*p1p2 - 12*p1k1*p1k3*p1p2 
     &+ 16*k2k3*p1p2**2 - k2k3**2*p2k1 + 12*k2k3*M2*p2k1 - 10*k2k3*p1k1*p2k1 + 12*M2*p1k1*p2k1 - 12
     &*p1k1**2*p2k1 - 27*k2k3*p1k2*p2k1 + 4*M2*p1k2*p2k1 - 30*p1k1*p1k2*p2k1 - 14*p1k2**2*p2k1 + 5
     &*k2k3*p1k3*p2k1 - 4*M2*p1k3*p2k1 + 20*p1k1*p1k3*p2k1 + 6*p1k2*p1k3*p2k1 - 6*p1k3**2*p2k1 - 12
     &*k2k3*p1p2*p2k1 - 28*p1k1*p1p2*p2k1 - 8*p1k2*p1p2*p2k1 - 2*k2k3*p2k1**2 + 8*M2*p2k1**2 + 6*p1k1
     &*p2k1**2 - 2*p1k2*p2k1**2 + 4*p1k3*p2k1**2 + k2k3**2*p2k2 + 4*k2k3*M2*p2k2 + 9*k2k3*p1k1*p2k2 
     &+ 8*M2*p1k1*p2k2 + 6*p1k1**2*p2k2 - 6*k2k3*p1k2*p2k2 + 2*p1k1*p1k2*p2k2 + 3*k2k3*p1k3*p2k2 
     &- 8*M2*p1k3*p2k2 + 24*p1k1*p1k3*p2k2 + 10*p1k2*p1k3*p2k2 - 8*k2k3*p1p2*p2k2 - 12*p1k1*p1p2*p2k2 
     &- 12*p1k3*p1p2*p2k2 - k2k3*p2k1*p2k2 + 4*M2*p2k1*p2k2 + 14*p1k1*p2k1*p2k2 - 6*p1k2*p2k1*p2k2 
     &+ 14*p1k3*p2k1*p2k2 + 12*p1k1*p2k2**2 + 6*p1k3*p2k2**2 - 6*k2k3*M2*p2k3 + 3*k2k3*p1k1*p2k3 
     &+ 16*M2*p1k1*p2k3 + 28*p1k1**2*p2k3 + 11*k2k3*p1k2*p2k3 + 8*M2*p1k2*p2k3 + 18*p1k1*p1k2*p2k3 
     &- 10*p1k2**2*p2k3 - 6*p1k1*p1k3*p2k3 + 6*k2k3*p1p2*p2k3 - 20*p1k1*p1p2*p2k3 - 28*p1k2*p1p2*p2k3 
     &- 3*k2k3*p2k1*p2k3 + 4*M2*p2k1*p2k3 - 28*p1k1*p2k1*p2k3 - 16*p1k2*p2k1*p2k3 + 6*p1k3*p2k1
     &*p2k3 - 3*k2k3*p2k2*p2k3 + 20*M2*p2k2*p2k3 - 22*p1k1*p2k2*p2k3 - 6*p1k2*p2k2*p2k3 + k1k2**2*(
     &-14*M2 + 3*p1k1 + 8*p1k2 - 7*p1k3 + 14*p1p2 - 3*p2k1 - 8*p2k2 + 7*p2k3) + k1k3*(-8*M2**2 - 6
     &*p1k1*p1k2 - 4*p1k2**2 + p1k2*p1k3 - 12*p1k1*p1p2 - 4*p1k2*p1p2 + 10*p1k3*p1p2 + 8*p1p2**2 + 6
     &*p1k1*p2k2 + 4*p1k2*p2k2 - 5*p1k3*p2k2 - 4*p1p2*p2k2 + 3*k2k3*(2*M2 - p1k2 - 2*p1p2 + p2k2) 
     &+ 2*M2*(6*p1k1 + 6*p1k2 - 5*p1k3 - 2*p2k2 - 3*p2k3) + 7*p1k2*p2k3 + 6*p1p2*p2k3 - 3*p2k2*p2k3) 
     &+ k1k2*(-8*M2**2 - 9*p1k1**2 + 2*M2*p1k2 - 2*p1k1*p1k2 - 18*M2*p1k3 + 24*p1k1*p1k3 + 21*p1k2
     &*p1k3 - 4*p1k3**2 - 8*p1k1*p1p2 + 18*p1k2*p1p2 + 14*p1k3*p1p2 + 8*p1p2**2 + 4*M2*p2k1 + 12*p1k1
     &*p2k1 - 2*p1k2*p2k1 - 10*p1k3*p2k1 + 4*p1p2*p2k1 - 3*p2k1**2 - 10*M2*p2k2 + 6*p1k1*p2k2 - 33
     &*p1k3*p2k2 - 10*p1p2*p2k2 - 2*p2k1*p2k2 - 3*k1k3*(2*M2 - p1k2 - 2*p1p2 + p2k2) + k2k3*(-4*M2 
     &- 5*p1k1 - 4*p1k2 + 2*p1k3 + 4*p1p2 + 5*p2k1 + 4*p2k2 - 2*p2k3) - 10*M2*p2k3 - 28*p1k1*p2k3 
     &+ 3*p1k2*p2k3 + 4*p1k3*p2k3 + 14*p1p2*p2k3 + 14*p2k1*p2k3 + 9*p2k2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)
     &/(128*k1k2*k1k3*p1k2*(M2 + p1p2)))

!************************************************   2 Re(Mbtk2 Mvk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(-5*k2k3**2 + 2*k1k3*M2 - 3
     &*k1k3*p1k2 - 4*p1k1*p1k2 + 10*p1k1*p1k3 - 14*p1k2*p1k3 - 10*p1k3**2 - 2*k1k3*p1p2 - 4*p1k2*p2k1 
     &- 8*p1k3*p2k1 + k1k3*p2k2 + 8*p1k1*p2k2 + 10*p1k3*p2k2 + k2k3*(8*M2 - 6*p1k1 - 3*p1k2 + 15
     &*p1k3 - 8*p1p2 + 4*p2k1 - p2k2 - 5*p2k3) - 2*p1k1*p2k3 + 4*p1k2*p2k3 + 10*p1k3*p2k3 + k1k2*(5
     &*k2k3 + 10*M2 - 15*p1k2 - 9*p1k3 - 10*p1p2 + p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k3*k2k3*p1k2))

!************************************************   2 Re(Mbtk2 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*(2 - 3*Nc**2 + Nc**4)*(2*k1k2**2 + 2*k2k3**2 
     &+ k2k3*(-21*p1k1 - 6*p1k3 + 7*p2k1 + 2*p2k3) + k1k2*(14*k2k3 - 6*p1k1 - 21*p1k3 + 2*p2k1 + 7
     &*p2k3) - 2*(-2*p1k1**2 + 4*k1k3*(2*M2 - 3*p1k2 - 2*p1p2 + p2k2) + p1k3*(-2*p1k3 + 7*p2k1 + 2*p2k3) 
     &+ p1k1*(-14*p1k3 + 2*p2k1 + 7*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32.
     &*k1k3**2*p1k2))

!************************************************   2 Re(Mbtk2 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*Nc**2*(-1 + Nc**2)*(5*k1k2**2 - 10*k2k3*M2 
     &+ 9*k2k3*p1k1 + 10*p1k1**2 + 15*k2k3*p1k2 + 14*p1k1*p1k2 - 10*p1k1*p1k3 + 4*p1k2*p1k3 + 10*k2k3
     &*p1p2 - k2k3*p2k1 - 10*p1k1*p2k1 - 4*p1k2*p2k1 + 2*p1k3*p2k1 - k2k3*p2k2 - 10*p1k1*p2k2 - 8
     &*p1k3*p2k2 - k1k3*(2*M2 - 3*p1k2 - 2*p1p2 + p2k2) + k1k2*(-5*k2k3 - 8*M2 - 15*p1k1 + 3*p1k2 + 6
     &*p1k3 + 8*p1p2 + 5*p2k1 + p2k2 - 4*p2k3) + 8*p1k1*p2k3 + 4*p1k2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k2*k1k3*p1k2))

!************************************************   2 Re(Mbtk2 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2 + 2*M2 - 2*p1k1 - 2*p1k2)*(k1k2
     &**2 - k1k2*k1k3 + 4*k2k3*M2 - 4*k2k3*p1k1 + 2*M2*p1k1 - 3*p1k1**2 - 4*k2k3*p1k2 - 2*p1k1*p1k2 
     &+ 4*M2*p1k3 - 4*p1k2*p1k3 - 4*k2k3*p1p2 - 2*p1k1*p1p2 - 4*p1k3*p1p2 + 2*k2k3*p2k1 + 2*M2*p2k1 
     &+ 2*p1k1*p2k1 - 2*p1k2*p2k1 + 2*p1k3*p2k1 - 2*p1p2*p2k1 + p2k1**2 + k1k3*(-2*M2 + p1k1 + 2
     &*p1k2 + 2*p1p2 - p2k1 - 4*p2k2) + 4*p1k1*p2k2 + 8*p1k3*p2k2 + 4*M2*p2k3 - 4*p1k1*p2k3 - 4*p1k2
     &*p2k3 - 4*p1p2*p2k3 + 2*p2k1*p2k3 + 2*k1k2*(k2k3 + M2 + p1k1 - p1k2 + p1k3 - p1p2 + p2k1 + p2k3))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(64*k1k3*p1k1*p1k2*(-k1k2 + p1k1 + p1k2)))

!************************************************   2 Re(Mbtk2 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2 + k2k3 + 2*M2 - 2*p1k1 - p1k2 
     &- 2*p1k3 - 2*p1p2 + p2k2)*(-4*k2k3*M2 + 4*k2k3*p1k1 + 6*p1k1**2 - 4*p1k1*p1k2 + k1k2*(k1k3 
     &+ 4*M2 - 6*p1k1 - 8*p1k3) - 2*k2k3*p1k3 + 4*p1k1*p1k3 + 4*p1k2*p1k3 + 2*p1k3**2 - 4*p1k1*p1p2 
     &+ 4*p1k3*p1p2 + 4*M2*p2k1 - 6*p1k1*p2k1 - 8*p1k3*p2k1 - 4*M2*p2k3 + 4*p1k1*p2k3 - 2*p1k3*p2k3 
     &+ k1k3*(k2k3 - 4*M2 - p1k1 + 4*p1k2 - p1k3 + 4*p1p2 + p2k1 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**5)/(128*k1k3*p1k1*p1k2*(k1k3 - p1k1 - p1k3)))

!************************************************   2 Re(Mbtk2 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k2k3 + 2*M2 - 2*p1k2 - 2*p1k3)*(k2k3
     &**2 + 2*k2k3*M2 + 2*k2k3*p1k1 + 4*M2*p1k1 - 2*k2k3*p1k2 - 4*p1k1*p1k2 + 2*k2k3*p1k3 + 2*M2
     &*p1k3 - 2*p1k2*p1k3 - 3*p1k3**2 - 2*k2k3*p1p2 - 4*p1k1*p1p2 - 2*p1k3*p1p2 + 2*k2k3*p2k1 + 4*M2
     &*p2k1 - 4*p1k2*p2k1 - 4*p1k3*p2k1 - 4*p1p2*p2k1 + 8*p1k1*p2k2 + 4*p1k3*p2k2 + 2*k2k3*p2k3 + 2
     &*M2*p2k3 + 2*p1k1*p2k3 - 2*p1k2*p2k3 + 2*p1k3*p2k3 - 2*p1p2*p2k3 + 2*p2k1*p2k3 + p2k3**2 + 2
     &*k1k2*(k2k3 + 2*M2 - 2*p1k2 - 2*p1k3 - 2*p1p2 + p2k3) - k1k3*(k2k3 + 2*M2 - 2*p1k2 - p1k3 - 2
     &*p1p2 + 4*p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(64*k1k3*p1k2*p1k3*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mbtk2 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)**2*(M2 - p1k2)*(3*k2k3**2 + 4*k1k2
     &*M2 - 4*M2*p1k1 - 2*M2*p1k3 + 3*p1k3**2 - 4*k1k2*p1p2 + 4*p1k1*p1p2 + 2*p1k3*p1p2 - 4*M2*p2k1 
     &+ 8*p1k2*p2k1 + 4*p1k3*p2k1 + 4*p1p2*p2k1 + 4*k1k2*p2k2 - 4*p1k1*p2k2 - 2*p1k3*p2k2 - 4*p2k1
     &*p2k2 + 2*k1k2*p2k3 - 2*M2*p2k3 - 2*p1k1*p2k3 + 4*p1k2*p2k3 - 2*p1k3*p2k3 + 2*p1p2*p2k3 - 2*p2k1
     &*p2k3 - 2*p2k2*p2k3 - p2k3**2 + 2*k2k3*(M2 - 3*p1k3 - p1p2 - 2*p2k1 + p2k2 + p2k3) + k1k3*(k2k3 
     &+ 2*M2 - p1k3 - 2*p1p2 + 2*p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(32*k1k3*p1k2**2*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mbtk2 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)**2*(M2 - p1k2)*(3*k1k2**2 + 4*k2k3
     &*M2 - 2*M2*p1k1 + 3*p1k1**2 - 4*M2*p1k3 - 4*k2k3*p1p2 + 2*p1k1*p1p2 + 4*p1k3*p1p2 + 2*k2k3*p2k1 
     &- 2*M2*p2k1 - 2*p1k1*p2k1 + 4*p1k2*p2k1 - 2*p1k3*p2k1 + 2*p1p2*p2k1 - p2k1**2 + 4*k2k3*p2k2 
     &- 2*p1k1*p2k2 - 4*p1k3*p2k2 - 2*p2k1*p2k2 + k1k3*(2*M2 - p1k1 - 2*p1p2 + p2k1 + 2*p2k2) + k1k2
     &*(k1k3 + 2*(M2 - 3*p1k1 - p1p2 + p2k1 + p2k2 - 2*p2k3)) - 4*M2*p2k3 + 4*p1k1*p2k3 + 8*p1k2
     &*p2k3 + 4*p1p2*p2k3 - 2*p2k1*p2k3 - 4*p2k2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(32*k1k3*p1k2**2*(-k1k2 + p1k1 + p1k2)))

!************************************************   2 Re(Mbtk2 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(-1 + Nc**2)*(k1k2 + k2k3 + 2*M2 - 2*p1k1 - p1k2 
     &- 2*p1k3 - 2*p1p2 + p2k2)*(-4*k2k3*M2 + 8*k2k3*p1k1 - 2*p1k1**2 - 4*p1k1*p1k2 + 6*k2k3*p1k3 
     &- 4*p1k1*p1k3 + 4*p1k2*p1k3 - 6*p1k3**2 - k1k2*(k1k3 - 4*M2 - 2*p1k1 + 4*p1k3) - 4*p1k1*p1p2 
     &+ 4*p1k3*p1p2 + 4*M2*p2k1 + 2*p1k1*p2k1 - 4*p1k3*p2k1 - 4*M2*p2k3 + 8*p1k1*p2k3 + 6*p1k3*p2k3 
     &- k1k3*(k2k3 - 4*M2 - p1k1 + 4*p1k2 - p1k3 + 4*p1p2 + p2k1 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**5)/(128*k1k3*p1k2*(k1k3 - p1k1 - p1k3)*p1k3))

!************************************************   2 Re(Mbtk3 Msk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**3*Nc**2*(-1 + Nc**2)*(-14*k2k3**2*M2 - 8*k2k3
     &*M2**2 - 7*k2k3**2*p1k1 - 18*k2k3*M2*p1k1 - 4*k2k3*p1k1**2 + 3*k2k3**2*p1k2 + 24*k2k3*p1k1*p1k2 
     &+ 16*M2*p1k1*p1k2 + 6*p1k1**2*p1k2 - 9*k2k3*p1k2**2 + 20*M2*p1k2**2 - 24*p1k1*p1k2**2 + 6*p1k2
     &**3 + 8*k2k3**2*p1k3 + 2*k2k3*M2*p1k3 + 21*k2k3*p1k1*p1k3 + 20*M2*p1k1*p1k3 - 2*k2k3*p1k2
     &*p1k3 + 16*M2*p1k2*p1k3 - 24*p1k1*p1k2*p1k3 + 12*p1k2**2*p1k3 + 6*p1k2*p1k3**2 + 14*k2k3**2*p1p2 
     &+ 14*k2k3*p1k1*p1p2 - 8*k2k3*p1k2*p1p2 - 12*p1k1*p1k2*p1p2 - 12*p1k2**2*p1p2 + 18*k2k3*p1k3
     &*p1p2 - 12*p1k2*p1k3*p1p2 + 8*k2k3*p1p2**2 + 7*k2k3**2*p2k1 - 10*k2k3*M2*p2k1 + 4*k2k3*p1k1
     &*p2k1 - 28*k2k3*p1k2*p2k1 + 16*M2*p1k2*p2k1 - 6*p1k1*p1k2*p2k1 + 28*p1k2**2*p2k1 + 3*k2k3*p1k3
     &*p2k1 + 8*M2*p1k3*p2k1 + 18*p1k2*p1k3*p2k1 - 10*p1k3**2*p2k1 + 14*k2k3*p1p2*p2k1 - 20*p1k2*p1p2
     &*p2k1 - 28*p1k3*p1p2*p2k1 - 3*k2k3**2*p2k2 + 4*k2k3*M2*p2k2 - 10*k2k3*p1k1*p2k2 - 4*M2*p1k1
     &*p2k2 - 6*p1k1**2*p2k2 + 12*k2k3*p1k2*p2k2 + 12*M2*p1k2*p2k2 + 20*p1k1*p1k2*p2k2 - 12*p1k2**2
     &*p2k2 - 2*k2k3*p1k3*p2k2 + 4*M2*p1k3*p2k2 + 6*p1k1*p1k3*p2k2 - 30*p1k2*p1k3*p2k2 - 14*p1k3**2
     &*p2k2 + 4*k2k3*p1p2*p2k2 - 28*p1k2*p1p2*p2k2 - 8*p1k3*p1p2*p2k2 + 14*k2k3*p2k1*p2k2 + 4*M2*p2k1
     &*p2k2 + 6*p1k1*p2k1*p2k2 - 28*p1k2*p2k1*p2k2 - 16*p1k3*p2k1*p2k2 - 3*k2k3*p2k2**2 + 8*M2*p2k2
     &**2 + 4*p1k1*p2k2**2 + 6*p1k2*p2k2**2 - 2*p1k3*p2k2**2 - 8*k2k3**2*p2k3 - 10*k2k3*M2*p2k3 
     &- 33*k2k3*p1k1*p2k3 - 8*M2*p1k1*p2k3 + 6*k2k3*p1k2*p2k3 + 8*M2*p1k2*p2k3 + 24*p1k1*p1k2*p2k3 
     &+ 6*p1k2**2*p2k3 + 10*p1k1*p1k3*p2k3 + 2*p1k2*p1k3*p2k3 - 10*k2k3*p1p2*p2k3 - 12*p1k1*p1p2*p2k3 
     &- 12*p1k2*p1p2*p2k3 + 9*k2k3*p2k1*p2k3 + 20*M2*p2k1*p2k3 - 22*p1k2*p2k1*p2k3 - 6*p1k3*p2k1
     &*p2k3 - 2*k2k3*p2k2*p2k3 + 4*M2*p2k2*p2k3 + 14*p1k1*p2k2*p2k3 + 14*p1k2*p2k2*p2k3 - 6*p1k3*p2k2
     &*p2k3 + 6*p1k1*p2k3**2 + 12*p1k2*p2k3**2 + k1k3**2*(6*M2 + p1k2 - p1k3 - 6*p1p2 - p2k2 + p2k3) 
     &+ k1k2*(-8*M2**2 - 10*M2*p1k1 + 12*M2*p1k2 + 12*M2*p1k3 + p1k1*p1k3 - 6*p1k2*p1k3 - 4*p1k3
     &**2 + 10*p1k1*p1p2 - 12*p1k2*p1p2 - 4*p1k3*p1p2 + 8*p1p2**2 - 6*M2*p2k1 + 7*p1k3*p2k1 + 6*p1p2
     &*p2k1 - 4*M2*p2k3 - 5*p1k1*p2k3 + 6*p1k2*p2k3 + 4*p1k3*p2k3 - 4*p1p2*p2k3 - 3*p2k1*p2k3 + 3
     &*k1k3*(2*M2 - p1k3 - 2*p1p2 + p2k3) - 3*k2k3*(2*M2 - p1k3 - 2*p1p2 + p2k3)) + k1k3*(-16*M2**2 
     &- 5*p1k1*p1k2 + 12*p1k2**2 - 11*p1k1*p1k3 + 19*p1k2*p1k3 + 6*p1k3**2 + 2*p1k1*p1p2 - 8*p1k2*p1p2 
     &+ 16*p1p2**2 + 3*p1k2*p2k1 + 11*p1k3*p2k1 + 6*p1p2*p2k1 + 5*p1k1*p2k2 - 10*p1k2*p2k2 - 27
     &*p1k3*p2k2 - 12*p1p2*p2k2 - 3*p2k1*p2k2 - 2*p2k2**2 + 3*p1k1*p2k3 + 9*p1k2*p2k3 - 6*p1k3*p2k3 
     &- 8*p1p2*p2k3 - 3*p2k1*p2k3 - p2k2*p2k3 + k2k3*(-4*M2 + 2*p1k1 - 5*p1k2 - 4*p1k3 + 4*p1p2 - 2*p2k1 + 5*p2k2 + 4*p2k3) 
     &+ M2*(-2*p1k1 + 8*p1k2 + 4*p1k3 - 6*p2k1 + 12*p2k2 + 4*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2*k2k3*p1k3*(M2 + p1p2)))

!************************************************   2 Re(Mbtk3 Msk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**3*Nc**2*(-1 + Nc**2)*(6*k2k3**2*M2 - 16*k2k3*M2
     &**2 + k2k3**2*p1k1 + 8*k2k3*M2*p1k1 + 12*k2k3*p1k1**2 + 20*M2*p1k1**2 + 6*p1k1**3 - 2*k2k3*M2
     &*p1k2 - 5*k2k3*p1k1*p1k2 + 16*M2*p1k1*p1k2 - 24*p1k1**2*p1k2 + 6*p1k1*p1k2**2 - k2k3**2*p1k3 
     &+ 4*k2k3*M2*p1k3 + 19*k2k3*p1k1*p1k3 + 16*M2*p1k1*p1k3 + 12*p1k1**2*p1k3 - 11*k2k3*p1k2*p1k3 
     &+ 20*M2*p1k2*p1k3 - 24*p1k1*p1k2*p1k3 + 6*k2k3*p1k3**2 + 6*p1k1*p1k3**2 - 6*k2k3**2*p1p2 - 8
     &*k2k3*p1k1*p1p2 - 12*p1k1**2*p1p2 + 2*k2k3*p1k2*p1p2 - 12*p1k1*p1k2*p1p2 - 12*p1k1*p1k3*p1p2 
     &+ 16*k2k3*p1p2**2 - k2k3**2*p2k1 + 12*k2k3*M2*p2k1 - 10*k2k3*p1k1*p2k1 + 12*M2*p1k1*p2k1 - 12
     &*p1k1**2*p2k1 + 5*k2k3*p1k2*p2k1 - 4*M2*p1k2*p2k1 + 20*p1k1*p1k2*p2k1 - 6*p1k2**2*p2k1 - 27*k2k3
     &*p1k3*p2k1 + 4*M2*p1k3*p2k1 - 30*p1k1*p1k3*p2k1 + 6*p1k2*p1k3*p2k1 - 14*p1k3**2*p2k1 - 12
     &*k2k3*p1p2*p2k1 - 28*p1k1*p1p2*p2k1 - 8*p1k3*p1p2*p2k1 - 2*k2k3*p2k1**2 + 8*M2*p2k1**2 + 6*p1k1
     &*p2k1**2 + 4*p1k2*p2k1**2 - 2*p1k3*p2k1**2 - 6*k2k3*M2*p2k2 + 3*k2k3*p1k1*p2k2 + 16*M2*p1k1
     &*p2k2 + 28*p1k1**2*p2k2 - 6*p1k1*p1k2*p2k2 + 11*k2k3*p1k3*p2k2 + 8*M2*p1k3*p2k2 + 18*p1k1*p1k3
     &*p2k2 - 10*p1k3**2*p2k2 + 6*k2k3*p1p2*p2k2 - 20*p1k1*p1p2*p2k2 - 28*p1k3*p1p2*p2k2 - 3*k2k3*p2k1
     &*p2k2 + 4*M2*p2k1*p2k2 - 28*p1k1*p2k1*p2k2 + 6*p1k2*p2k1*p2k2 - 16*p1k3*p2k1*p2k2 + k1k3*
     &*2*(-14*M2 + 3*p1k1 - 7*p1k2 + 8*p1k3 + 14*p1p2 - 3*p2k1 + 7*p2k2 - 8*p2k3) + k2k3**2*p2k3 + 4
     &*k2k3*M2*p2k3 + 9*k2k3*p1k1*p2k3 + 8*M2*p1k1*p2k3 + 6*p1k1**2*p2k3 + 3*k2k3*p1k2*p2k3 - 8*M2
     &*p1k2*p2k3 + 24*p1k1*p1k2*p2k3 - 6*k2k3*p1k3*p2k3 + 2*p1k1*p1k3*p2k3 + 10*p1k2*p1k3*p2k3 - 8*k2k3
     &*p1p2*p2k3 - 12*p1k1*p1p2*p2k3 - 12*p1k2*p1p2*p2k3 - k2k3*p2k1*p2k3 + 4*M2*p2k1*p2k3 + 14
     &*p1k1*p2k1*p2k3 + 14*p1k2*p2k1*p2k3 - 6*p1k3*p2k1*p2k3 - 3*k2k3*p2k2*p2k3 + 20*M2*p2k2*p2k3 - 22
     &*p1k1*p2k2*p2k3 - 6*p1k3*p2k2*p2k3 + 12*p1k1*p2k3**2 + 6*p1k2*p2k3**2 + k1k2*(-8*M2**2 + 12
     &*M2*p1k1 - 10*M2*p1k2 + 12*M2*p1k3 - 6*p1k1*p1k3 + p1k2*p1k3 - 4*p1k3**2 - 12*p1k1*p1p2 + 10
     &*p1k2*p1p2 - 4*p1k3*p1p2 + 8*p1p2**2 - 6*M2*p2k2 + 7*p1k3*p2k2 + 6*p1p2*p2k2 - 4*M2*p2k3 + 6*p1k1
     &*p2k3 - 5*p1k2*p2k3 + 4*p1k3*p2k3 - 4*p1p2*p2k3 - 3*p2k2*p2k3 - 3*k1k3*(2*M2 - p1k3 - 2*p1p2 
     &+ p2k3) + 3*k2k3*(2*M2 - p1k3 - 2*p1p2 + p2k3)) + k1k3*(-8*M2**2 - 9*p1k1**2 + 24*p1k1*p1k2 
     &- 4*p1k2**2 - 2*p1k1*p1k3 + 21*p1k2*p1k3 - 8*p1k1*p1p2 + 14*p1k2*p1p2 + 18*p1k3*p1p2 + 8*p1p2
     &**2 + 12*p1k1*p2k1 - 10*p1k2*p2k1 - 2*p1k3*p2k1 + 4*p1p2*p2k1 - 3*p2k1**2 - 28*p1k1*p2k2 + 4
     &*p1k2*p2k2 + 3*p1k3*p2k2 + 14*p1p2*p2k2 + 14*p2k1*p2k2 + 6*p1k1*p2k3 - 33*p1k2*p2k3 - 10*p1p2
     &*p2k3 - 2*p2k1*p2k3 + 9*p2k2*p2k3 + k2k3*(-4*M2 - 5*p1k1 + 2*p1k2 - 4*p1k3 + 4*p1p2 + 5*p2k1 - 2*p2k2 + 4*p2k3) 
     &- 2*M2*(9*p1k2 - p1k3 - 2*p2k1 + 5*p2k2 + 5*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2*k1k3*p1k3*(M2 + p1p2)))

!************************************************   2 Re(Mbtk3 Msk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**3*Nc**2*(-1 + Nc**2)*(4*k2k3**2*M2 + 5*k2k3**2
     &*p1k1 - 14*k2k3*M2*p1k1 + 5*k2k3*p1k1**2 - 4*k2k3*M2*p1k2 - 5*k2k3*p1k1*p1k2 - 14*k2k3*p1k1*p1k3 
     &+ 4*p1k1**2*p1k3 - 4*k2k3*p1k2*p1k3 + 28*p1k1*p1k2*p1k3 + 4*p1k2**2*p1k3 - 4*k2k3**2*p1p2 
     &- 14*k2k3*p1k1*p1p2 + 8*p1k1**2*p1p2 - 4*k2k3*p1k2*p1p2 + 56*p1k1*p1k2*p1p2 + 8*p1k2**2*p1p2 
     &- 5*k2k3**2*p2k1 + 42*k2k3*M2*p2k1 - 10*k2k3*p1k1*p2k1 - 8*M2*p1k1*p2k1 + 15*k2k3*p1k2*p2k1 
     &- 28*M2*p1k2*p2k1 + 10*p1k1*p1k2*p2k1 - 10*p1k2**2*p2k1 - 14*k2k3*p1k3*p2k1 - 14*k2k3*p1p2*p2k1 
     &- 8*p1k1*p1p2*p2k1 - 28*p1k2*p1p2*p2k1 + 5*k2k3*p2k1**2 + 8*M2*p2k1**2 - 10*p1k2*p2k1**2 - 4
     &*p1k3*p2k1**2 + k1k3**2*(4*M2 + 5*p1k2 - 4*p1p2 - 5*p2k2) + 12*k2k3*M2*p2k2 - 5*k2k3*p1k1*p2k2 
     &- 28*M2*p1k1*p2k2 - 10*p1k1**2*p2k2 - 8*M2*p1k2*p2k2 + 10*p1k1*p1k2*p2k2 - 4*k2k3*p1k3*p2k2 
     &- 4*k2k3*p1p2*p2k2 - 28*p1k1*p1p2*p2k2 - 8*p1k2*p1p2*p2k2 - 5*k2k3*p2k1*p2k2 + 56*M2*p2k1*p2k2 
     &+ 10*p1k1*p2k1*p2k2 + 10*p1k2*p2k1*p2k2 - 28*p1k3*p2k1*p2k2 + 8*M2*p2k2**2 - 10*p1k1*p2k2
     &**2 - 4*p1k3*p2k2**2 + 28*k2k3*p1k1*p2k3 - 8*p1k1**2*p2k3 + 8*k2k3*p1k2*p2k3 - 56*p1k1*p1k2*p2k3 
     &- 8*p1k2**2*p2k3 + 8*p1k1*p2k1*p2k3 + 28*p1k2*p2k1*p2k3 + 28*p1k1*p2k2*p2k3 + 8*p1k2*p2k2
     &*p2k3 - k1k3*(5*p1k1*p1k2 - 5*p1k2**2 + 4*p1k1*p1k3 + 14*p1k2*p1k3 + 4*p1k1*p1p2 + 14*p1k2*p1p2 
     &+ 5*p1k2*p2k1 + 4*p1k3*p2k1 + 4*p1p2*p2k1 + 2*M2*(2*p1k1 + 7*p1k2 - 6*p2k1 - 21*p2k2) + k2k3
     &*(-28*M2 + 5*p1k1 + 5*p1k2 + 28*p1p2 - 5*p2k1 - 5*p2k2) - 15*p1k1*p2k2 + 10*p1k2*p2k2 + 14
     &*p1k3*p2k2 + 14*p1p2*p2k2 + 5*p2k1*p2k2 - 5*p2k2**2 - 8*p1k1*p2k3 - 28*p1k2*p2k3) - 8*k1k2*(4
     &*M2**2 + 2*M2*p1k1 + 2*M2*p1k2 - 2*M2*p1k3 + p1k1*p1k3 + p1k2*p1k3 - 2*p1k3**2 - 2*p1k1*p1p2 
     &- 2*p1k2*p1p2 - 2*p1k3*p1p2 - 4*p1p2**2 + 2*M2*p2k1 - 3*p1k3*p2k1 - 2*p1p2*p2k1 + 2*M2*p2k2 - 3
     &*p1k3*p2k2 - 2*p1p2*p2k2 + k1k3*(-2*M2 + p1k3 + 2*p1p2 - p2k3) + k2k3*(-2*M2 + p1k3 + 2*p1p2 
     &- p2k3) + 2*M2*p2k3 + p1k1*p2k3 + p1k2*p2k3 + 2*p1k3*p2k3 + 2*p1p2*p2k3 + p2k1*p2k3 + p2k2*p2k3))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(64*k1k2**2*p1k3*(M2 + p1p2)))

!************************************************   2 Re(Mbtk3 Mvk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(-5*k2k3**2 + 2*k1k2*M2 + 10
     &*p1k1*p1k2 - 10*p1k2**2 - 3*k1k2*p1k3 - 4*p1k1*p1k3 - 14*p1k2*p1k3 - 2*k1k2*p1p2 - 8*p1k2*p2k1 
     &- 4*p1k3*p2k1 - 2*p1k1*p2k2 + 10*p1k2*p2k2 + 4*p1k3*p2k2 + k2k3*(8*M2 - 6*p1k1 + 15*p1k2 
     &- 3*p1k3 - 8*p1p2 + 4*p2k1 - 5*p2k2 - p2k3) + k1k2*p2k3 + 8*p1k1*p2k3 + 10*p1k2*p2k3 + k1k3*(5
     &*k2k3 + 10*M2 - 9*p1k2 - 15*p1k3 - 10*p1p2 + p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k2*k2k3*p1k3))

!************************************************   2 Re(Mbtk3 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*Nc**2*(-1 + Nc**2)*(5*k1k3**2 - 10*k2k3*M2 
     &+ 9*k2k3*p1k1 + 10*p1k1**2 - 10*p1k1*p1k2 + 15*k2k3*p1k3 + 14*p1k1*p1k3 + 4*p1k2*p1k3 + 10*k2k3
     &*p1p2 - k2k3*p2k1 - 10*p1k1*p2k1 + 2*p1k2*p2k1 - 4*p1k3*p2k1 + 8*p1k1*p2k2 + 4*p1k3*p2k2 - k2k3
     &*p2k3 - 10*p1k1*p2k3 - 8*p1k2*p2k3 - k1k2*(2*M2 - 3*p1k3 - 2*p1p2 + p2k3) + k1k3*(-5*k2k3 
     &- 8*M2 - 15*p1k1 + 6*p1k2 + 3*p1k3 + 8*p1p2 + 5*p2k1 - 4*p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k2*k1k3*p1k3))

!************************************************   2 Re(Mbtk3 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*(2 - 3*Nc**2 + Nc**4)*(2*k1k3**2 + 2*k2k3**2 
     &+ k2k3*(-21*p1k1 - 6*p1k2 + 7*p2k1 + 2*p2k2) + k1k3*(14*k2k3 - 6*p1k1 - 21*p1k2 + 2*p2k1 + 7
     &*p2k2) - 2*(-2*p1k1**2 + p1k2*(-2*p1k2 + 7*p2k1 + 2*p2k2) + p1k1*(-14*p1k2 + 2*p2k1 + 7*p2k2) 
     &+ 4*k1k2*(2*M2 - 3*p1k3 - 2*p1p2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32.
     &*k1k2**2*p1k3))

!************************************************   2 Re(Mbtk3 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2*(k1k3 + k2k3 - 4*M2 - p1k1 - p1k2 
     &+ 4*p1k3 + 4*p1p2 + p2k1 + p2k2) - 2*(-3*p1k1**2 - 2*p1k1*p1k2 - p1k2**2 + k2k3*(2*M2 - 2
     &*p1k1 + p1k2) + k1k3*(-2*M2 + 3*p1k1 + 4*p1k2) + 2*p1k1*p1k3 - 2*p1k2*p1k3 + 2*p1k1*p1p2 - 2
     &*p1k2*p1p2 - 2*M2*p2k1 + 3*p1k1*p2k1 + 4*p1k2*p2k1 + 2*M2*p2k2 - 2*p1k1*p2k2 + p1k2*p2k2))*(k1k3 
     &+ k2k3 + 2*M2 - 2*p1k1 - 2*p1k2 - p1k3 - 2*p1p2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**5)/(128*k1k2*p1k1*(k1k2 - p1k1 - p1k2)*p1k3))

!************************************************   2 Re(Mbtk3 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k3 + 2*M2 - 2*p1k1 - 2*p1k3)*(k1k3
     &**2 + 4*k2k3*M2 - 4*k2k3*p1k1 + 2*M2*p1k1 - 3*p1k1**2 + 4*M2*p1k2 - 4*k2k3*p1k3 - 2*p1k1*p1k3 
     &- 4*p1k2*p1k3 - 4*k2k3*p1p2 - 2*p1k1*p1p2 - 4*p1k2*p1p2 + 2*k2k3*p2k1 + 2*M2*p2k1 + 2*p1k1
     &*p2k1 + 2*p1k2*p2k1 - 2*p1k3*p2k1 - 2*p1p2*p2k1 + p2k1**2 + 4*M2*p2k2 - 4*p1k1*p2k2 - 4*p1k3
     &*p2k2 - 4*p1p2*p2k2 + 2*p2k1*p2k2 + 2*k1k3*(k2k3 + M2 + p1k1 + p1k2 - p1k3 - p1p2 + p2k1 + p2k2) 
     &+ 4*p1k1*p2k3 + 8*p1k2*p2k3 - k1k2*(k1k3 + 2*M2 - p1k1 - 2*p1k3 - 2*p1p2 + p2k1 + 4*p2k3))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(64*k1k2*p1k1*p1k3*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mbtk3 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)**2*(M2 - p1k3)*(3*k2k3**2 + 4*k1k3
     &*M2 - 4*M2*p1k1 - 2*M2*p1k2 + 3*p1k2**2 - 4*k1k3*p1p2 + 4*p1k1*p1p2 + 2*p1k2*p1p2 - 4*M2*p2k1 
     &+ 4*p1k2*p2k1 + 8*p1k3*p2k1 + 4*p1p2*p2k1 + 2*k1k3*p2k2 - 2*M2*p2k2 - 2*p1k1*p2k2 - 2*p1k2*p2k2 
     &+ 4*p1k3*p2k2 + 2*p1p2*p2k2 - 2*p2k1*p2k2 - p2k2**2 + 4*k1k3*p2k3 - 4*p1k1*p2k3 - 2*p1k2*p2k3 
     &- 4*p2k1*p2k3 - 2*p2k2*p2k3 + 2*k2k3*(M2 - 3*p1k2 - p1p2 - 2*p2k1 + p2k2 + p2k3) + k1k2*(k2k3 
     &+ 2*M2 - p1k2 - 2*p1p2 + p2k2 + 2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(32*k1k2*p1k3**2*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mbtk3 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k2k3 + 2*M2 - 2*p1k2 - 2*p1k3)*(k2k3
     &**2 + 2*k2k3*M2 + 2*k2k3*p1k1 + 4*M2*p1k1 + 2*k2k3*p1k2 + 2*M2*p1k2 - 3*p1k2**2 - 2*k2k3*p1k3 
     &- 4*p1k1*p1k3 - 2*p1k2*p1k3 - 2*k2k3*p1p2 - 4*p1k1*p1p2 - 2*p1k2*p1p2 + 2*k2k3*p2k1 + 4*M2
     &*p2k1 - 4*p1k2*p2k1 - 4*p1k3*p2k1 - 4*p1p2*p2k1 + 2*k2k3*p2k2 + 2*M2*p2k2 + 2*p1k1*p2k2 + 2*p1k2
     &*p2k2 - 2*p1k3*p2k2 - 2*p1p2*p2k2 + 2*p2k1*p2k2 + p2k2**2 + 2*k1k3*(k2k3 + 2*M2 - 2*p1k2 
     &- 2*p1k3 - 2*p1p2 + p2k2) + 8*p1k1*p2k3 + 4*p1k2*p2k3 - k1k2*(k2k3 + 2*M2 - p1k2 - 2*p1k3 - 2
     &*p1p2 + p2k2 + 4*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(64*k1k2*p1k2*p1k3*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mbtk3 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(4*k2k3*M2 - 8*k2k3*p1k1 + 2*p1k1*
     &*2 - 2*k1k3*(2*M2 + p1k1 - 2*p1k2) - 6*k2k3*p1k2 + 4*p1k1*p1k2 + 6*p1k2**2 + 4*p1k1*p1k3 - 4*p1k2
     &*p1k3 + 4*p1k1*p1p2 - 4*p1k2*p1p2 - 4*M2*p2k1 - 2*p1k1*p2k1 + 4*p1k2*p2k1 + 4*M2*p2k2 - 8
     &*p1k1*p2k2 - 6*p1k2*p2k2 + k1k2*(k1k3 + k2k3 - 4*M2 - p1k1 - p1k2 + 4*p1k3 + 4*p1p2 + p2k1 + p2k2))
     &*(k1k3 + k2k3 + 2*M2 - 2*p1k1 - 2*p1k2 - p1k3 - 2*p1p2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**5)/(128*k1k2*(k1k2 - p1k1 - p1k2)*p1k2*p1k3))

!************************************************   2 Re(Mbtk3 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)**2*(M2 - p1k3)*(3*k1k3**2 + 4*k2k3
     &*M2 - 2*M2*p1k1 + 3*p1k1**2 - 4*M2*p1k2 - 4*k2k3*p1p2 + 2*p1k1*p1p2 + 4*p1k2*p1p2 + 2*k2k3*p2k1 
     &- 2*M2*p2k1 - 2*p1k1*p2k1 - 2*p1k2*p2k1 + 4*p1k3*p2k1 + 2*p1p2*p2k1 - p2k1**2 - 4*M2*p2k2 
     &+ 4*p1k1*p2k2 + 8*p1k3*p2k2 + 4*p1p2*p2k2 - 2*p2k1*p2k2 + 4*k2k3*p2k3 - 2*p1k1*p2k3 - 4*p1k2*p2k3 
     &- 2*p2k1*p2k3 - 4*p2k2*p2k3 + 2*k1k3*(M2 - 3*p1k1 - p1p2 + p2k1 - 2*p2k2 + p2k3) + k1k2*(k1k3 
     &+ 2*M2 - p1k1 - 2*p1p2 + p2k1 + 2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**5)/(32*k1k2*p1k3**2*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Msk1 Msk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**4*Nc**2*(-1 + Nc**2)*(-2*k2k3**2*M2 - 16*k2k3
     &*M2**2 + 6*k2k3*M2*p1k1 - 2*k2k3*M2*p1k2 + 3*k2k3*p1k1*p1k2 + 40*M2*p1k1*p1k2 - 4*k2k3*p1k2**2 
     &+ 8*k2k3*M2*p1k3 + 9*k2k3*p1k1*p1k3 + 20*M2*p1k1*p1k3 - 13*k2k3*p1k2*p1k3 + 20*M2*p1k2*p1k3 
     &+ 7*k2k3*p1k3**2 + 28*M2*p1k3**2 + 4*k1k2**2*(M2 - p1p2) - 2*k1k3**2*(M2 - p1p2) + 2*k2k3**2
     &*p1p2 - 6*k2k3*p1k1*p1p2 + 2*k2k3*p1k2*p1p2 - 8*k2k3*p1k3*p1p2 + 16*k2k3*p1p2**2 + 6*k2k3*M2*p2k1 
     &- 3*k2k3*p1k2*p2k1 - 8*p1k1*p1k2*p2k1 + 8*p1k2**2*p2k1 - 9*k2k3*p1k3*p2k1 - 4*p1k1*p1k3*p2k1 
     &+ 4*p1k2*p1k3*p2k1 + 12*p1k3**2*p2k1 - 6*k2k3*p1p2*p2k1 - 40*p1k2*p1p2*p2k1 - 20*p1k3*p1p2
     &*p2k1 + 8*p1k2*p2k1**2 + 4*p1k3*p2k1**2 - 2*k2k3*M2*p2k2 - 3*k2k3*p1k1*p2k2 + 8*p1k1**2*p2k2 
     &+ 8*k2k3*p1k2*p2k2 - 8*p1k1*p1k2*p2k2 + 13*k2k3*p1k3*p2k2 + 4*p1k1*p1k3*p2k2 - 4*p1k2*p1k3*p2k2 
     &+ 12*p1k3**2*p2k2 + 2*k2k3*p1p2*p2k2 - 40*p1k1*p1p2*p2k2 - 20*p1k3*p1p2*p2k2 + 3*k2k3*p2k1
     &*p2k2 + 40*M2*p2k1*p2k2 - 8*p1k1*p2k1*p2k2 - 8*p1k2*p2k1*p2k2 - 8*p1k3*p2k1*p2k2 - 4*k2k3*p2k2
     &**2 + 8*p1k1*p2k2**2 + 4*p1k3*p2k2**2 + 8*k2k3*M2*p2k3 - 9*k2k3*p1k1*p2k3 + 4*p1k1**2*p2k3 
     &+ 13*k2k3*p1k2*p2k3 - 8*p1k1*p1k2*p2k3 + 4*p1k2**2*p2k3 - 14*k2k3*p1k3*p2k3 - 12*p1k1*p1k3*p2k3 
     &- 12*p1k2*p1k3*p2k3 - 8*k2k3*p1p2*p2k3 - 20*p1k1*p1p2*p2k3 - 20*p1k2*p1p2*p2k3 - 56*p1k3*p1p2
     &*p2k3 + 9*k2k3*p2k1*p2k3 + 20*M2*p2k1*p2k3 - 4*p1k1*p2k1*p2k3 + 4*p1k2*p2k1*p2k3 - 12*p1k3
     &*p2k1*p2k3 - 13*k2k3*p2k2*p2k3 + 20*M2*p2k2*p2k3 + 4*p1k1*p2k2*p2k3 - 4*p1k2*p2k2*p2k3 - 12*p1k3
     &*p2k2*p2k3 + 7*k2k3*p2k3**2 + 28*M2*p2k3**2 + 12*p1k1*p2k3**2 + 12*p1k2*p2k3**2 + k1k3*(-16
     &*M2**2 - 4*p1k1**2 + 3*p1k1*p1k2 - 13*p1k1*p1k3 + 9*p1k2*p1k3 + 7*p1k3**2 + 8*k2k3*(M2 - p1p2) 
     &+ 2*p1k1*p1p2 - 6*p1k2*p1p2 - 8*p1k3*p1p2 + 16*p1p2**2 + 8*p1k1*p2k1 - 3*p1k2*p2k1 + 13*p1k3
     &*p2k1 + 2*p1p2*p2k1 - 4*p2k1**2 - 3*p1k1*p2k2 - 9*p1k3*p2k2 - 6*p1p2*p2k2 + 3*p2k1*p2k2 - 2
     &*M2*(p1k1 - 3*p1k2 - 4*p1k3 + p2k1 - 3*p2k2 - 4*p2k3) + 13*p1k1*p2k3 - 9*p1k2*p2k3 - 14*p1k3
     &*p2k3 - 8*p1p2*p2k3 - 13*p2k1*p2k3 + 9*p2k2*p2k3 + 7*p2k3**2) + k1k2*(-32*M2**2 - 4*M2*p1k1 - 8
     &*p1k1**2 - 4*M2*p1k2 + 14*p1k1*p1k2 - 8*p1k2**2 - 4*M2*p1k3 - p1k1*p1k3 - p1k2*p1k3 + 23*p1k3
     &**2 + 4*k1k3*(M2 - p1p2) + 4*k2k3*(M2 - p1p2) + 4*p1k1*p1p2 + 4*p1k2*p1p2 + 4*p1k3*p1p2 + 32
     &*p1p2**2 - 4*M2*p2k1 + 16*p1k1*p2k1 - 14*p1k2*p2k1 + p1k3*p2k1 + 4*p1p2*p2k1 - 8*p2k1**2 - 4
     &*M2*p2k2 - 14*p1k1*p2k2 + 16*p1k2*p2k2 + p1k3*p2k2 + 4*p1p2*p2k2 + 14*p2k1*p2k2 - 8*p2k2**2 
     &- 4*M2*p2k3 + p1k1*p2k3 + p1k2*p2k3 - 46*p1k3*p2k3 + 4*p1p2*p2k3 - p2k1*p2k3 - p2k2*p2k3 + 23*p2k3**2))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k3*k2k3*(M2 + p1p2)**2))

!************************************************   2 Re(Msk1 Msk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**4*Nc**2*(-1 + Nc**2)*(2*k2k3**2*M2 + 16*k2k3*M2
     &**2 - 6*k2k3*M2*p1k1 - 8*k2k3*M2*p1k2 - 9*k2k3*p1k1*p1k2 - 20*M2*p1k1*p1k2 - 7*k2k3*p1k2**2 
     &- 28*M2*p1k2**2 + 2*k2k3*M2*p1k3 - 3*k2k3*p1k1*p1k3 - 40*M2*p1k1*p1k3 + 13*k2k3*p1k2*p1k3 - 20
     &*M2*p1k2*p1k3 + 4*k2k3*p1k3**2 + 2*k1k2**2*(M2 - p1p2) - 4*k1k3**2*(M2 - p1p2) - 2*k2k3**2*p1p2 
     &+ 6*k2k3*p1k1*p1p2 + 8*k2k3*p1k2*p1p2 - 2*k2k3*p1k3*p1p2 - 16*k2k3*p1p2**2 - 6*k2k3*M2*p2k1 
     &+ 9*k2k3*p1k2*p2k1 + 4*p1k1*p1k2*p2k1 - 12*p1k2**2*p2k1 + 3*k2k3*p1k3*p2k1 + 8*p1k1*p1k3*p2k1 
     &- 4*p1k2*p1k3*p2k1 - 8*p1k3**2*p2k1 + 6*k2k3*p1p2*p2k1 + 20*p1k2*p1p2*p2k1 + 40*p1k3*p1p2
     &*p2k1 - 4*p1k2*p2k1**2 - 8*p1k3*p2k1**2 - 8*k2k3*M2*p2k2 + 9*k2k3*p1k1*p2k2 - 4*p1k1**2*p2k2 
     &+ 14*k2k3*p1k2*p2k2 + 12*p1k1*p1k2*p2k2 - 13*k2k3*p1k3*p2k2 + 8*p1k1*p1k3*p2k2 + 12*p1k2*p1k3
     &*p2k2 - 4*p1k3**2*p2k2 + 8*k2k3*p1p2*p2k2 + 20*p1k1*p1p2*p2k2 + 56*p1k2*p1p2*p2k2 + 20*p1k3*p1p2
     &*p2k2 - 9*k2k3*p2k1*p2k2 - 20*M2*p2k1*p2k2 + 4*p1k1*p2k1*p2k2 + 12*p1k2*p2k1*p2k2 - 4*p1k3
     &*p2k1*p2k2 - 7*k2k3*p2k2**2 - 28*M2*p2k2**2 - 12*p1k1*p2k2**2 - 12*p1k3*p2k2**2 + 2*k2k3*M2*p2k3 
     &+ 3*k2k3*p1k1*p2k3 - 8*p1k1**2*p2k3 - 13*k2k3*p1k2*p2k3 - 4*p1k1*p1k2*p2k3 - 12*p1k2**2*p2k3 
     &- 8*k2k3*p1k3*p2k3 + 8*p1k1*p1k3*p2k3 + 4*p1k2*p1k3*p2k3 - 2*k2k3*p1p2*p2k3 + 40*p1k1*p1p2
     &*p2k3 + 20*p1k2*p1p2*p2k3 - 3*k2k3*p2k1*p2k3 - 40*M2*p2k1*p2k3 + 8*p1k1*p2k1*p2k3 + 8*p1k2*p2k1
     &*p2k3 + 8*p1k3*p2k1*p2k3 + 13*k2k3*p2k2*p2k3 - 20*M2*p2k2*p2k3 - 4*p1k1*p2k2*p2k3 + 12*p1k2
     &*p2k2*p2k3 + 4*p1k3*p2k2*p2k3 + 4*k2k3*p2k3**2 - 8*p1k1*p2k3**2 - 4*p1k2*p2k3**2 + k1k2*(-4*k1k3
     &*M2 - 8*k2k3*M2 + 16*M2**2 + 2*M2*p1k1 + 4*p1k1**2 - 8*M2*p1k2 + 13*p1k1*p1k2 - 7*p1k2**2 
     &- 6*M2*p1k3 - 3*p1k1*p1k3 - 9*p1k2*p1k3 + 4*k1k3*p1p2 + 8*k2k3*p1p2 - 2*p1k1*p1p2 + 8*p1k2*p1p2 
     &+ 6*p1k3*p1p2 - 16*p1p2**2 + 2*M2*p2k1 - 8*p1k1*p2k1 - 13*p1k2*p2k1 + 3*p1k3*p2k1 - 2*p1p2
     &*p2k1 + 4*p2k1**2 - 8*M2*p2k2 - 13*p1k1*p2k2 + 14*p1k2*p2k2 + 9*p1k3*p2k2 + 8*p1p2*p2k2 + 13
     &*p2k1*p2k2 - 7*p2k2**2 - 6*M2*p2k3 + 3*p1k1*p2k3 + 9*p1k2*p2k3 + 6*p1p2*p2k3 - 3*p2k1*p2k3 - 9
     &*p2k2*p2k3) + k1k3*(32*M2**2 + 8*p1k1**2 + p1k1*p1k2 - 23*p1k2**2 - 14*p1k1*p1k3 + p1k2*p1k3 
     &+ 8*p1k3**2 - 4*k2k3*(M2 - p1p2) - 4*p1k1*p1p2 - 4*p1k2*p1p2 - 4*p1k3*p1p2 - 32*p1p2**2 - 16
     &*p1k1*p2k1 - p1k2*p2k1 + 14*p1k3*p2k1 - 4*p1p2*p2k1 + 8*p2k1**2 - p1k1*p2k2 + 46*p1k2*p2k2 - p1k3
     &*p2k2 - 4*p1p2*p2k2 + p2k1*p2k2 - 23*p2k2**2 + 14*p1k1*p2k3 - p1k2*p2k3 - 16*p1k3*p2k3 - 4*p1p2*p2k3 
     &- 14*p2k1*p2k3 + p2k2*p2k3 + 8*p2k3**2 + 4*M2*(p1k1 + p1k2 + p1k3 + p2k1 + p2k2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)
     &/(64*k1k2*k2k3*(M2 + p1p2)**2))

!************************************************   2 Re(Msk1 Mvk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Msk1 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Msk1 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Msk1 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(k1k3**2*(-6*M2 + 2*p1k1 
     &+ 3*p1k2 + 6*p1p2 - 2*p2k1 - 3*p2k2) + 2*k1k2**2*(3*M2 - p1k1 - 3*p1k3 - 3*p1p2 + p2k1 + 3*p2k3) 
     &+ k1k3*(12*M2**2 - 4*p1k1**2 - 10*p1k1*p1k2 - 6*p1k2**2 - 3*p1k2*p1k3 - 16*p1k1*p1p2 - 26
     &*p1k2*p1p2 - 10*p1k3*p1p2 - 20*p1p2**2 + 4*p1k1*p2k1 + 14*p1k2*p2k1 + 2*p1k3*p2k1 + 12*p1p2*p2k1 
     &+ k2k3*(-2*M2 + p1k2 + 2*p1p2 - p2k2) - 8*p1k1*p2k2 + 6*p1k2*p2k2 + 9*p1k3*p2k2 - 6*p1p2*p2k2 
     &+ 4*p2k1*p2k2 - 8*p1k1*p2k3 - 7*p1k2*p2k3 - 10*p1p2*p2k3 + 6*p2k1*p2k3 + p2k2*p2k3 + 2*M2
     &*(9*p1k2 + 7*p1k3 + 4*p1p2 + 2*p2k1 + 7*p2k2 + 3*p2k3)) - k1k2*(-12*M2**2 - 4*p1k1**2 + 6*p1k1
     &*p1k2 + k2k3*p1k3 - 10*p1k1*p1k3 - 6*p1k2*p1k3 - 3*p1k3**2 - 4*p1k1*p1p2 + 2*p1k2*p1p2 - 20
     &*p1k3*p1p2 + 4*p1p2**2 + 4*p1k1*p2k1 - 4*p1k2*p2k1 + 20*p1k3*p2k1 + 4*p1k1*p2k2 - 12*p1k3*p2k2 
     &+ 6*p1p2*p2k2 - 6*p2k1*p2k2 + 3*k1k3*(2*M2 - 2*p1k2 + p1k3 - 2*p1p2 + 2*p2k2 - p2k3) - k2k3
     &*p2k3 - 2*p1k1*p2k3 + 18*p1k2*p2k3 + 2*p1k3*p2k3 + 8*p1p2*p2k3 - 8*p2k1*p2k3 + p2k3**2 + 2*M2
     &*(6*p1k1 + p1k2 + 4*p1p2 - 4*p2k1 - 5*p2k2 + 6*p2k3)) + 2*(2*p1k1**2*p1k2 + 4*p1k1*p1k2**2 - 2
     &*p1k1**2*p1k3 - p1k1*p1k3**2 + 12*p1k1*p1k2*p1p2 + 4*p1k2**2*p1p2 - 6*p1k1*p1k3*p1p2 + 2*p1k3
     &**2*p1p2 + 8*p1k2*p1p2**2 + 4*p1k3*p1p2**2 - 8*p1k1*p1k2*p2k1 - 3*p1k2**2*p2k1 + 8*p1k1*p1k3
     &*p2k1 + 3*p1k2*p1k3*p2k1 - 10*p1k2*p1p2*p2k1 + 4*p1k3*p1p2*p2k1 + 4*p1k2*p2k1**2 - 4*p1k3*p2k1
     &**2 + 8*p1k1**2*p2k2 - p1k1*p1k2*p2k2 - 3*p1k1*p1k3*p2k2 - 6*p1k2*p1k3*p2k2 - 3*p1k3**2*p2k2 
     &+ 2*p1k1*p1p2*p2k2 + 8*p1k2*p1p2*p2k2 + 6*p1k3*p1p2*p2k2 - 6*p1k1*p2k1*p2k2 - 7*p1k2*p2k1*p2k2 
     &- 4*p1k3*p2k1*p2k2 + 7*p1k1*p2k2**2 - 2*p1k1**2*p2k3 + 2*p1k1*p1k2*p2k3 + 6*p1k2**2*p2k3 + 2
     &*p1k1*p1k3*p2k3 + 3*p1k2*p1k3*p2k3 - 2*p1k1*p1p2*p2k3 + 10*p1k2*p1p2*p2k3 + 6*p1k3*p1p2*p2k3 
     &- 9*p1k2*p2k1*p2k3 - 4*p1k3*p2k1*p2k3 + 11*p1k1*p2k2*p2k3 - p1k3*p2k2*p2k3 + 3*p1k1*p2k3**2 
     &+ p1k2*p2k3**2 - 4*M2**2*(p1k2 + 2*p1k3 - p2k2 + p2k3) + k2k3*(4*M2**2 - 2*p1k2*p1p2 - 4*p1p2
     &**2 + p1k2*p2k1 - p1k3*p2k1 + 2*p1p2*p2k1 + p1k3*p2k2 - 2*p1p2*p2k2 + p2k1*p2k2 - p1k1*(p1k2 
     &+ 2*p1p2 + p2k2) + 2*M2*(-p1k1 + p1k2 + p1k3 + p2k1 + p2k2) - p1k2*p2k3 - 2*p1p2*p2k3 + p2k1
     &*p2k3) - 2*M2*(p1k2**2 + 2*p1k3**2 - 2*p1k3*p1p2 + 4*p1k3*p2k1 + 2*p1k3*p2k2 + 2*p1p2*p2k2 - 2
     &*p2k1*p2k2 + 3*p2k2**2 + p1k1*(-3*p1k3 + 5*p2k2 - 5*p2k3) + p1k2*(3*p1k3 + 2*p1p2 - p2k1 
     &+ 2*p2k2 - 2*p2k3) + p1k3*p2k3 - 2*p1p2*p2k3 + 2*p2k1*p2k3 + 5*p2k2*p2k3 + p2k3**2)))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(128*k2k3*p1k1*(-k1k2 + p1k1 + p1k2)*(M2 + p1p2)))

!************************************************   2 Re(Msk1 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*Nc**2*(-1 + Nc**2)*(-2*k1k3**2*(3*M2 - p1k1 
     &- 3*p1k2 - 3*p1p2 + p2k1 + 3*p2k2) + k1k2**2*(6*M2 - 2*p1k1 - 3*p1k3 - 6*p1p2 + 2*p2k1 + 3*p2k3) 
     &+ k1k3*(-12*M2**2 - 4*p1k1**2 + k2k3*p1k2 - 3*p1k2**2 - 6*p1k2*p1k3 - 20*p1k2*p1p2 + 2*p1k3
     &*p1p2 + 4*p1p2**2 + 20*p1k2*p2k1 - 4*p1k3*p2k1 - k2k3*p2k2 + 2*p1k2*p2k2 + 18*p1k3*p2k2 + 8
     &*p1p2*p2k2 - 8*p2k1*p2k2 + p2k2**2 + 2*M2*(6*p1k1 + p1k3 + 4*p1p2 - 4*p2k1 + 6*p2k2 - 5*p2k3) 
     &- 2*p1k1*(5*p1k2 - 3*p1k3 + 2*p1p2 - 2*p2k1 + p2k2 - 2*p2k3) - 12*p1k2*p2k3 + 6*p1p2*p2k3 - 6
     &*p2k1*p2k3) + k1k2*(-12*M2**2 + 4*p1k1**2 - 14*M2*p1k2 - 18*M2*p1k3 + 10*p1k1*p1k3 + 3*p1k2
     &*p1k3 + 6*p1k3**2 - 8*M2*p1p2 + 16*p1k1*p1p2 + 10*p1k2*p1p2 + 26*p1k3*p1p2 + 20*p1p2**2 - 4*M2
     &*p2k1 - 4*p1k1*p2k1 - 2*p1k2*p2k1 - 14*p1k3*p2k1 - 12*p1p2*p2k1 - 6*M2*p2k2 + 8*p1k1*p2k2 + 7
     &*p1k3*p2k2 + 10*p1p2*p2k2 - 6*p2k1*p2k2 - 14*M2*p2k3 + 8*p1k1*p2k3 - 9*p1k2*p2k3 - 6*p1k3*p2k3 
     &+ 6*p1p2*p2k3 - 4*p2k1*p2k3 - p2k2*p2k3 + k2k3*(2*M2 - p1k3 - 2*p1p2 + p2k3) + 3*k1k3*(2*M2 
     &+ p1k2 - 2*p1k3 - 2*p1p2 - p2k2 + 2*p2k3)) - 2*(-2*p1k1**2*p1k2 - p1k1*p1k2**2 + 2*p1k1**2*p1k3 
     &+ 4*p1k1*p1k3**2 - 6*p1k1*p1k2*p1p2 + 2*p1k2**2*p1p2 + 12*p1k1*p1k3*p1p2 + 4*p1k3**2*p1p2 
     &+ 4*p1k2*p1p2**2 + 8*p1k3*p1p2**2 + 8*p1k1*p1k2*p2k1 - 8*p1k1*p1k3*p2k1 + 3*p1k2*p1k3*p2k1 
     &- 3*p1k3**2*p2k1 + 4*p1k2*p1p2*p2k1 - 10*p1k3*p1p2*p2k1 - 4*p1k2*p2k1**2 + 4*p1k3*p2k1**2 - 2
     &*p1k1**2*p2k2 + 2*p1k1*p1k2*p2k2 + 2*p1k1*p1k3*p2k2 + 3*p1k2*p1k3*p2k2 + 6*p1k3**2*p2k2 - 2*p1k1
     &*p1p2*p2k2 + 6*p1k2*p1p2*p2k2 + 10*p1k3*p1p2*p2k2 - 4*p1k2*p2k1*p2k2 - 9*p1k3*p2k1*p2k2 + 3
     &*p1k1*p2k2**2 + p1k3*p2k2**2 - 4*M2**2*(2*p1k2 + p1k3 + p2k2 - p2k3) + 8*p1k1**2*p2k3 - 3*p1k1
     &*p1k2*p2k3 - 3*p1k2**2*p2k3 - p1k1*p1k3*p2k3 - 6*p1k2*p1k3*p2k3 + 2*p1k1*p1p2*p2k3 + 6*p1k2
     &*p1p2*p2k3 + 8*p1k3*p1p2*p2k3 - 6*p1k1*p2k1*p2k3 - 4*p1k2*p2k1*p2k3 - 7*p1k3*p2k1*p2k3 + 11*p1k1
     &*p2k2*p2k3 - p1k2*p2k2*p2k3 + 7*p1k1*p2k3**2 + k2k3*(4*M2**2 - 2*p1k3*p1p2 - 4*p1p2**2 - p1k2
     &*p2k1 + p1k3*p2k1 + 2*p1p2*p2k1 - p1k3*p2k2 - 2*p1p2*p2k2 + p2k1*p2k2 + p1k2*p2k3 - 2*p1p2
     &*p2k3 + p2k1*p2k3 - p1k1*(p1k3 + 2*p1p2 + p2k3) + 2*M2*(-p1k1 + p1k2 + p1k3 + p2k1 + p2k3)) - 2
     &*M2*(2*p1k2**2 + p1k3**2 + 2*p1k3*p1p2 - p1k3*p2k1 - 2*p1k3*p2k2 - 2*p1p2*p2k2 + 2*p2k1*p2k2 
     &+ p2k2**2 + 2*p1k3*p2k3 + 2*p1p2*p2k3 - 2*p2k1*p2k3 + 5*p2k2*p2k3 + 3*p2k3**2 + p1k2*(3*p1k3 
     &- 2*p1p2 + 4*p2k1 + p2k2 + 2*p2k3) + p1k1*(-3*p1k2 - 5*p2k2 + 5*p2k3))))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)
     &/(128*k2k3*p1k1*(-k1k3 + p1k1 + p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Msk1 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(8*k2k3**2*M2 + 8*k2k3*M2
     &**2 + 2*k2k3**2*p1k1 + 4*k2k3*M2*p1k1 + 2*k2k3*M2*p1k2 - 11*k2k3*p1k1*p1k2 + 8*M2*p1k1*p1k2 + 2
     &*k2k3*p1k2**2 - 4*M2*p1k2**2 + 2*p1k1*p1k2**2 - 4*p1k2**3 + 2*k2k3*M2*p1k3 - 3*k2k3*p1k1*p1k3 
     &- 8*M2*p1k1*p1k3 - 4*k2k3*p1k2*p1k3 - 8*M2*p1k2*p1k3 + 4*p1k1*p1k2*p1k3 - 4*p1k2**2*p1k3 - 14
     &*k2k3*p1k3**2 - 12*M2*p1k3**2 + 6*p1k1*p1k3**2 + 4*p1k2*p1k3**2 + 4*p1k3**3 - 8*k2k3**2*p1p2 
     &- 4*k2k3*p1k1*p1p2 - 18*k2k3*p1k2*p1p2 + 8*M2*p1k2*p1p2 - 4*p1k1*p1k2*p1p2 - 4*p1k2**2*p1p2 
     &- 2*k2k3*p1k3*p1p2 - 8*M2*p1k3*p1p2 + 4*p1k1*p1k3*p1p2 + 8*p1k2*p1k3*p1p2 + 20*p1k3**2*p1p2 
     &- 8*k2k3*p1p2**2 - 8*p1k2*p1p2**2 + 8*p1k3*p1p2**2 - 2*k2k3**2*p2k1 + 4*k2k3*M2*p2k1 + 7*k2k3
     &*p1k2*p2k1 - 6*p1k2**2*p2k1 - k2k3*p1k3*p2k1 - 4*p1k2*p1k3*p2k1 - 2*p1k3**2*p2k1 - 4*k2k3*p1p2
     &*p2k1 - 4*p1k2*p1p2*p2k1 + 4*p1k3*p1p2*p2k1 + 18*k2k3*M2*p2k2 - 8*M2**2*p2k2 + 5*k2k3*p1k1*p2k2 
     &- 4*M2*p1k1*p2k2 - 2*k2k3*p1k2*p2k2 + 4*M2*p1k2*p2k2 + 2*p1k1*p1k2*p2k2 + 4*p1k2**2*p2k2 
     &- 6*k2k3*p1k3*p2k2 + 16*M2*p1k3*p2k2 - 12*p1k2*p1k3*p2k2 - 24*p1k3**2*p2k2 - 2*k2k3*p1p2*p2k2 
     &+ 8*M2*p1p2*p2k2 + 4*p1k2*p1p2*p2k2 - 8*p1k3*p1p2*p2k2 - k2k3*p2k1*p2k2 + 4*M2*p2k1*p2k2 + 2
     &*p1k2*p2k1*p2k2 - 4*p1k3*p2k1*p2k2 + 2*k2k3*M2*p2k3 + 8*M2**2*p2k3 + 5*k2k3*p1k1*p2k3 + 4*M2
     &*p1k1*p2k3 + 10*k2k3*p1k2*p2k3 - 24*M2*p1k2*p2k3 - 4*p1k1*p1k2*p2k3 + 16*p1k2**2*p2k3 + 14*k2k3
     &*p1k3*p2k3 - 20*M2*p1k3*p2k3 - 10*p1k1*p1k3*p2k3 + 20*p1k2*p1k3*p2k3 - 4*p1k3**2*p2k3 - 2*k2k3
     &*p1p2*p2k3 - 8*M2*p1p2*p2k3 + 16*p1k2*p1p2*p2k3 + 12*p1k3*p1p2*p2k3 - k2k3*p2k1*p2k3 - 4*M2
     &*p2k1*p2k3 + 8*p1k2*p2k1*p2k3 + 6*p1k3*p2k1*p2k3 + k1k3*(k2k3 + 4*M2 - 8*p1k2 - 6*p1k3)*(2*M2 
     &- p1k2 - p1k3 - 2*p1p2 + p2k2 + p2k3) + k1k2*(k2k3 - 4*M2 - 2*p1k2 + 4*p1k3)*(2*M2 - p1k2 - p1k3 
     &- 2*p1p2 + p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(128*k2k3*(k2k3 - p1k2 - p1k3)*p1k3*(M2 + p1p2)))

!************************************************   2 Re(Msk1 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(8*k2k3**2*M2 + 8*k2k3*M2
     &**2 + 2*k2k3**2*p1k1 + 4*k2k3*M2*p1k1 + 2*k2k3*M2*p1k2 - 3*k2k3*p1k1*p1k2 - 8*M2*p1k1*p1k2 - 14
     &*k2k3*p1k2**2 - 12*M2*p1k2**2 + 6*p1k1*p1k2**2 + 4*p1k2**3 + 2*k2k3*M2*p1k3 - 11*k2k3*p1k1*p1k3 
     &+ 8*M2*p1k1*p1k3 - 4*k2k3*p1k2*p1k3 - 8*M2*p1k2*p1k3 + 4*p1k1*p1k2*p1k3 + 4*p1k2**2*p1k3 
     &+ 2*k2k3*p1k3**2 - 4*M2*p1k3**2 + 2*p1k1*p1k3**2 - 4*p1k2*p1k3**2 - 4*p1k3**3 - 8*k2k3**2*p1p2 
     &- 4*k2k3*p1k1*p1p2 - 2*k2k3*p1k2*p1p2 - 8*M2*p1k2*p1p2 + 4*p1k1*p1k2*p1p2 + 20*p1k2**2*p1p2 
     &- 18*k2k3*p1k3*p1p2 + 8*M2*p1k3*p1p2 - 4*p1k1*p1k3*p1p2 + 8*p1k2*p1k3*p1p2 - 4*p1k3**2*p1p2 
     &- 8*k2k3*p1p2**2 + 8*p1k2*p1p2**2 - 8*p1k3*p1p2**2 - 2*k2k3**2*p2k1 + 4*k2k3*M2*p2k1 - k2k3*p1k2
     &*p2k1 - 2*p1k2**2*p2k1 + 7*k2k3*p1k3*p2k1 - 4*p1k2*p1k3*p2k1 - 6*p1k3**2*p2k1 - 4*k2k3*p1p2
     &*p2k1 + 4*p1k2*p1p2*p2k1 - 4*p1k3*p1p2*p2k1 + 2*k2k3*M2*p2k2 + 8*M2**2*p2k2 + 5*k2k3*p1k1*p2k2 
     &+ 4*M2*p1k1*p2k2 + 14*k2k3*p1k2*p2k2 - 20*M2*p1k2*p2k2 - 10*p1k1*p1k2*p2k2 - 4*p1k2**2*p2k2 
     &+ 10*k2k3*p1k3*p2k2 - 24*M2*p1k3*p2k2 - 4*p1k1*p1k3*p2k2 + 20*p1k2*p1k3*p2k2 + 16*p1k3**2*p2k2 
     &- 2*k2k3*p1p2*p2k2 - 8*M2*p1p2*p2k2 + 12*p1k2*p1p2*p2k2 + 16*p1k3*p1p2*p2k2 - k2k3*p2k1*p2k2 
     &- 4*M2*p2k1*p2k2 + 6*p1k2*p2k1*p2k2 + 8*p1k3*p2k1*p2k2 + 18*k2k3*M2*p2k3 - 8*M2**2*p2k3 
     &+ 5*k2k3*p1k1*p2k3 - 4*M2*p1k1*p2k3 - 6*k2k3*p1k2*p2k3 + 16*M2*p1k2*p2k3 - 24*p1k2**2*p2k3 - 2
     &*k2k3*p1k3*p2k3 + 4*M2*p1k3*p2k3 + 2*p1k1*p1k3*p2k3 - 12*p1k2*p1k3*p2k3 + 4*p1k3**2*p2k3 - 2
     &*k2k3*p1p2*p2k3 + 8*M2*p1p2*p2k3 - 8*p1k2*p1p2*p2k3 + 4*p1k3*p1p2*p2k3 - k2k3*p2k1*p2k3 + 4*M2
     &*p2k1*p2k3 - 4*p1k2*p2k1*p2k3 + 2*p1k3*p2k1*p2k3 + k1k2*(k2k3 + 4*M2 - 6*p1k2 - 8*p1k3)*(2*M2 
     &- p1k2 - p1k3 - 2*p1p2 + p2k2 + p2k3) + k1k3*(k2k3 - 4*M2 + 4*p1k2 - 2*p1k3)*(2*M2 - p1k2 - p1k3 
     &- 2*p1p2 + p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(128*k2k3*p1k2*(k2k3 - p1k2 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Msk1 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Msk1 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Msk2 Msk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**4*Nc**2*(-1 + Nc**2)*(-4*k2k3**2*M2 + 32*k2k3*M2
     &**2 + 4*k2k3*M2*p1k1 - 23*k2k3*p1k1**2 - 28*M2*p1k1**2 + 4*k2k3*M2*p1k2 + k2k3*p1k1*p1k2 - 20
     &*M2*p1k1*p1k2 + 8*k2k3*p1k2**2 + 4*k2k3*M2*p1k3 + k2k3*p1k1*p1k3 - 20*M2*p1k1*p1k3 - 14*k2k3
     &*p1k2*p1k3 - 40*M2*p1k2*p1k3 + 8*k2k3*p1k3**2 + 2*k1k2**2*(M2 - p1p2) + 2*k1k3**2*(M2 - p1p2) 
     &+ 4*k2k3**2*p1p2 - 4*k2k3*p1k1*p1p2 - 4*k2k3*p1k2*p1p2 - 4*k2k3*p1k3*p1p2 - 32*k2k3*p1p2**2 
     &+ 4*k2k3*M2*p2k1 + 46*k2k3*p1k1*p2k1 - k2k3*p1k2*p2k1 + 12*p1k1*p1k2*p2k1 - 4*p1k2**2*p2k1 - k2k3
     &*p1k3*p2k1 + 12*p1k1*p1k3*p2k1 + 8*p1k2*p1k3*p2k1 - 4*p1k3**2*p2k1 - 4*k2k3*p1p2*p2k1 + 56
     &*p1k1*p1p2*p2k1 + 20*p1k2*p1p2*p2k1 + 20*p1k3*p1p2*p2k1 - 23*k2k3*p2k1**2 - 28*M2*p2k1**2 - 12
     &*p1k2*p2k1**2 - 12*p1k3*p2k1**2 + 4*k2k3*M2*p2k2 - k2k3*p1k1*p2k2 - 12*p1k1**2*p2k2 - 16*k2k3
     &*p1k2*p2k2 + 4*p1k1*p1k2*p2k2 + 14*k2k3*p1k3*p2k2 - 4*p1k1*p1k3*p2k2 + 8*p1k2*p1k3*p2k2 - 8
     &*p1k3**2*p2k2 - 4*k2k3*p1p2*p2k2 + 20*p1k1*p1p2*p2k2 + 40*p1k3*p1p2*p2k2 + k2k3*p2k1*p2k2 - 20
     &*M2*p2k1*p2k2 + 12*p1k1*p2k1*p2k2 + 4*p1k2*p2k1*p2k2 - 4*p1k3*p2k1*p2k2 + 8*k2k3*p2k2**2 - 4
     &*p1k1*p2k2**2 - 8*p1k3*p2k2**2 + 4*k2k3*M2*p2k3 - k2k3*p1k1*p2k3 - 12*p1k1**2*p2k3 + 14*k2k3
     &*p1k2*p2k3 - 4*p1k1*p1k2*p2k3 - 8*p1k2**2*p2k3 - 16*k2k3*p1k3*p2k3 + 4*p1k1*p1k3*p2k3 + 8*p1k2
     &*p1k3*p2k3 - 4*k2k3*p1p2*p2k3 + 20*p1k1*p1p2*p2k3 + 40*p1k2*p1p2*p2k3 + k2k3*p2k1*p2k3 - 20*M2
     &*p2k1*p2k3 + 12*p1k1*p2k1*p2k3 - 4*p1k2*p2k1*p2k3 + 4*p1k3*p2k1*p2k3 - 14*k2k3*p2k2*p2k3 - 40
     &*M2*p2k2*p2k3 + 8*p1k1*p2k2*p2k3 + 8*p1k2*p2k2*p2k3 + 8*p1k3*p2k2*p2k3 + 8*k2k3*p2k3**2 - 4
     &*p1k1*p2k3**2 - 8*p1k2*p2k3**2 + k1k2*(-8*k1k3*M2 - 4*k2k3*M2 + 16*M2**2 - 8*M2*p1k1 - 7*p1k1
     &**2 + 2*M2*p1k2 + 13*p1k1*p1k2 + 4*p1k2**2 - 6*M2*p1k3 - 9*p1k1*p1k3 - 3*p1k2*p1k3 + 8*k1k3*p1p2 
     &+ 4*k2k3*p1p2 + 8*p1k1*p1p2 - 2*p1k2*p1p2 + 6*p1k3*p1p2 - 16*p1p2**2 - 8*M2*p2k1 + 14*p1k1
     &*p2k1 - 13*p1k2*p2k1 + 9*p1k3*p2k1 + 8*p1p2*p2k1 - 7*p2k1**2 + 2*M2*p2k2 - 13*p1k1*p2k2 - 8*p1k2
     &*p2k2 + 3*p1k3*p2k2 - 2*p1p2*p2k2 + 13*p2k1*p2k2 + 4*p2k2**2 - 6*M2*p2k3 + 9*p1k1*p2k3 + 3
     &*p1k2*p2k3 + 6*p1p2*p2k3 - 9*p2k1*p2k3 - 3*p2k2*p2k3) + k1k3*(16*M2**2 - 7*p1k1**2 - 9*p1k1*p1k2 
     &+ 13*p1k1*p1k3 - 3*p1k2*p1k3 + 4*p1k3**2 - 4*k2k3*(M2 - p1p2) + 8*p1k1*p1p2 + 6*p1k2*p1p2 
     &- 2*p1k3*p1p2 - 16*p1p2**2 + 14*p1k1*p2k1 + 9*p1k2*p2k1 - 13*p1k3*p2k1 + 8*p1p2*p2k1 - 7*p2k1
     &**2 + 9*p1k1*p2k2 + 3*p1k3*p2k2 + 6*p1p2*p2k2 - 9*p2k1*p2k2 - 2*M2*(4*p1k1 + 3*p1k2 - p1k3 
     &+ 4*p2k1 + 3*p2k2 - p2k3) - 13*p1k1*p2k3 + 3*p1k2*p2k3 - 8*p1k3*p2k3 - 2*p1p2*p2k3 + 13*p2k1*p2k3 
     &- 3*p2k2*p2k3 + 4*p2k3**2))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k2*k1k3*(M2 + p1p2)**2))

!************************************************   2 Re(Msk2 Mvk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Msk2 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Msk2 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Msk2 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Msk2 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(2*k1k3**2*(4*M2 + p1k2 - 4
     &*p1p2 - p2k2) + k1k2*(k1k3 + 4*M2 - 6*p1k1 - 8*p1k3)*(2*M2 - p1k1 - p1k3 - 2*p1p2 + p2k1 + p2k3) 
     &+ k1k3*(8*M2**2 - 14*p1k1**2 - 3*p1k1*p1k2 - 4*p1k1*p1k3 - 11*p1k2*p1k3 + 2*p1k3**2 - 2
     &*p1k1*p1p2 - 4*p1k2*p1p2 - 18*p1k3*p1p2 - 8*p1p2**2 + 14*p1k1*p2k1 + 5*p1k2*p2k1 + 10*p1k3*p2k1 
     &- 2*p1p2*p2k1 - p1k1*p2k2 + 7*p1k3*p2k2 - 4*p1p2*p2k2 - p2k1*p2k2 - 6*p1k1*p2k3 + 5*p1k2*p2k3 
     &- 2*p1k3*p2k3 - 2*p1p2*p2k3 - p2k2*p2k3 + k2k3*(2*M2 - p1k1 - p1k3 - 2*p1p2 + p2k1 + p2k3) 
     &+ 2*M2*(p1k1 + 2*p1k2 + p1k3 + p2k1 + 2*p2k2 + 9*p2k3)) - 2*(-2*p1k1**3 - 3*p1k1**2*p1k2 - 2
     &*p1k1**2*p1k3 - 2*p1k1*p1k2*p1k3 + 2*p1k1*p1k3**2 - p1k2*p1k3**2 + 2*p1k3**3 - 10*p1k1**2*p1p2 
     &- 2*p1k1*p1k2*p1p2 - 4*p1k1*p1k3*p1p2 + 2*p1k2*p1k3*p1p2 + 2*p1k3**2*p1p2 - 4*p1k1*p1p2**2 
     &+ 4*p1k3*p1p2**2 + 2*p1k1**2*p2k1 + 5*p1k1*p1k2*p2k1 - 10*p1k1*p1k3*p2k1 + 2*p1k2*p1k3*p2k1 
     &- 8*p1k3**2*p2k1 - 6*p1k1*p1p2*p2k1 - 8*p1k3*p1p2*p2k1 + p1k1**2*p2k2 + 2*p1k1*p1k3*p2k2 + 3*p1k3
     &**2*p2k2 - 2*p1k1*p1p2*p2k2 + 2*p1k3*p1p2*p2k2 - 3*p1k1*p2k1*p2k2 - 4*p1k3*p2k1*p2k2 - 4*M2
     &**2*(p2k1 - p2k3) + 12*p1k1**2*p2k3 + 6*p1k1*p1k3*p2k3 - p1k2*p1k3*p2k3 - 2*p1k3**2*p2k3 + 4
     &*p1k1*p1p2*p2k3 - 2*p1k3*p1p2*p2k3 + 2*p1k1*p2k2*p2k3 - p1k3*p2k2*p2k3 + k2k3*(2*M2 - 2*p1k1 
     &+ p1k3)*(2*M2 - p1k1 - p1k3 - 2*p1p2 + p2k1 + p2k3) + 2*M2*(3*p1k1**2 + p1k3**2 - 2*p1k3*p1p2 
     &+ 6*p1k3*p2k1 + 2*p1p2*p2k1 + p2k1*p2k2 + p1k1*(2*p1k2 + 2*p1k3 + 2*p1p2 + 5*p2k1 - 4*p2k3) 
     &- p1k3*p2k3 - 2*p1p2*p2k3 - p2k2*p2k3 + p1k2*(-2*p1k3 - p2k1 + p2k3))))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)
     &/(128*k1k3*p1k1*(k1k3 - p1k1 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Msk2 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Msk2 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*Nc**2*(-1 + Nc**2)*(-2*k2k3**2*(3*M2 - 3*p1k1 
     &- p1k2 - 3*p1p2 + 3*p2k1 + p2k2) + k1k2**2*(6*M2 - 2*p1k2 - 3*p1k3 - 6*p1p2 + 2*p2k2 + 3*p2k3) 
     &+ k2k3*(-12*M2**2 - 3*p1k1**2 - 10*p1k1*p1k2 - 4*p1k2**2 - 6*p1k1*p1k3 + 6*p1k2*p1k3 - 20
     &*p1k1*p1p2 - 4*p1k2*p1p2 + 2*p1k3*p1p2 + 4*p1p2**2 + k1k3*(p1k1 - p2k1) + 2*p1k1*p2k1 - 2*p1k2
     &*p2k1 + 18*p1k3*p2k1 + 8*p1p2*p2k1 + p2k1**2 + 20*p1k1*p2k2 + 4*p1k2*p2k2 - 4*p1k3*p2k2 - 8
     &*p2k1*p2k2 + 2*M2*(6*p1k2 + p1k3 + 4*p1p2 + 6*p2k1 - 4*p2k2 - 5*p2k3) - 12*p1k1*p2k3 + 4*p1k2
     &*p2k3 + 6*p1p2*p2k3 - 6*p2k2*p2k3) + k1k2*(-12*M2**2 - 14*M2*p1k1 + 4*p1k2**2 - 18*M2*p1k3 + 3
     &*p1k1*p1k3 + 10*p1k2*p1k3 + 6*p1k3**2 - 8*M2*p1p2 + 10*p1k1*p1p2 + 16*p1k2*p1p2 + 26*p1k3*p1p2 
     &+ 20*p1p2**2 - 6*M2*p2k1 + 8*p1k2*p2k1 + 7*p1k3*p2k1 + 10*p1p2*p2k1 - 4*M2*p2k2 - 2*p1k1*p2k2 
     &- 4*p1k2*p2k2 - 14*p1k3*p2k2 - 12*p1p2*p2k2 - 6*p2k1*p2k2 - 14*M2*p2k3 - 9*p1k1*p2k3 + 8*p1k2
     &*p2k3 - 6*p1k3*p2k3 + 6*p1p2*p2k3 - p2k1*p2k3 - 4*p2k2*p2k3 + k1k3*(2*M2 - p1k3 - 2*p1p2 
     &+ p2k3) + 3*k2k3*(2*M2 + p1k1 - 2*p1k3 - 2*p1p2 - p2k1 + 2*p2k3)) - 2*(-(p1k1**2*p1k2) - 2*p1k1
     &*p1k2**2 + 2*p1k2**2*p1k3 + 4*p1k2*p1k3**2 + 2*p1k1**2*p1p2 - 6*p1k1*p1k2*p1p2 + 12*p1k2*p1k3
     &*p1p2 + 4*p1k3**2*p1p2 + 4*p1k1*p1p2**2 + 8*p1k3*p1p2**2 + 2*p1k1*p1k2*p2k1 - 2*p1k2**2*p2k1 
     &+ 3*p1k1*p1k3*p2k1 + 2*p1k2*p1k3*p2k1 + 6*p1k3**2*p2k1 + 6*p1k1*p1p2*p2k1 - 2*p1k2*p1p2*p2k1 
     &+ 10*p1k3*p1p2*p2k1 + 3*p1k2*p2k1**2 + p1k3*p2k1**2 + 8*p1k1*p1k2*p2k2 + 3*p1k1*p1k3*p2k2 - 8
     &*p1k2*p1k3*p2k2 - 3*p1k3**2*p2k2 + 4*p1k1*p1p2*p2k2 - 10*p1k3*p1p2*p2k2 - 4*p1k1*p2k1*p2k2 
     &- 9*p1k3*p2k1*p2k2 - 4*p1k1*p2k2**2 + 4*p1k3*p2k2**2 - 4*M2**2*(2*p1k1 + p1k3 + p2k1 - p2k3) 
     &- 3*p1k1**2*p2k3 - 3*p1k1*p1k2*p2k3 + 8*p1k2**2*p2k3 - 6*p1k1*p1k3*p2k3 - p1k2*p1k3*p2k3 + 6*p1k1
     &*p1p2*p2k3 + 2*p1k2*p1p2*p2k3 + 8*p1k3*p1p2*p2k3 - p1k1*p2k1*p2k3 + 11*p1k2*p2k1*p2k3 - 4
     &*p1k1*p2k2*p2k3 - 6*p1k2*p2k2*p2k3 - 7*p1k3*p2k2*p2k3 + 7*p1k2*p2k3**2 + k1k3*(4*M2**2 - 2*p1k3
     &*p1p2 - 4*p1p2**2 - p1k3*p2k1 - 2*p1p2*p2k1 - p1k1*p2k2 + p1k3*p2k2 + 2*p1p2*p2k2 + p2k1*p2k2 
     &+ p1k1*p2k3 - 2*p1p2*p2k3 + p2k2*p2k3 - p1k2*(p1k3 + 2*p1p2 + p2k3) + 2*M2*(p1k1 - p1k2 + p1k3 
     &+ p2k2 + p2k3)) - 2*M2*(2*p1k1**2 + p1k3**2 - 5*p1k2*p2k1 - 2*p1p2*p2k1 + p2k1**2 + 2*p2k1
     &*p2k2 + 5*p1k2*p2k3 + 2*p1p2*p2k3 + 5*p2k1*p2k3 - 2*p2k2*p2k3 + 3*p2k3**2 + p1k3*(2*p1p2 - 2*p2k1 
     &- p2k2 + 2*p2k3) + p1k1*(-3*p1k2 + 3*p1k3 - 2*p1p2 + p2k1 + 4*p2k2 + 2*p2k3))))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(128*k1k3*p1k2*(-k2k3 + p1k2 + p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Msk2 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(-6*k2k3**2*M2 + 12*k2k3*M2
     &**2 + 3*k2k3**2*p1k1 + 18*k2k3*M2*p1k1 - 8*M2**2*p1k1 - 6*k2k3*p1k1**2 - 4*M2*p1k1**2 + 2*k2k3
     &**2*p1k2 - 10*k2k3*p1k1*p1k2 + 8*p1k1**2*p1k2 - 4*k2k3*p1k2**2 + 4*p1k1*p1k2**2 + 14*k2k3*M2
     &*p1k3 - 16*M2**2*p1k3 - 3*k2k3*p1k1*p1k3 - 12*M2*p1k1*p1k3 + 12*M2*p1k2*p1k3 - 4*p1k2**2*p1k3 
     &- 8*M2*p1k3**2 - 2*p1k2*p1k3**2 + 6*k2k3**2*p1p2 + 8*k2k3*M2*p1p2 - 26*k2k3*p1k1*p1p2 - 8*M2
     &*p1k1*p1p2 + 8*p1k1**2*p1p2 - 16*k2k3*p1k2*p1p2 + 24*p1k1*p1k2*p1p2 - 10*k2k3*p1k3*p1p2 + 8
     &*M2*p1k3*p1p2 - 12*p1k2*p1k3*p1p2 + 4*p1k3**2*p1p2 - 20*k2k3*p1p2**2 + 16*p1k1*p1p2**2 + 8*p1k3
     &*p1p2**2 + k1k3*k2k3*(-2*M2 + p1k1 + 2*p1p2 - p2k1) - 3*k2k3**2*p2k1 + 14*k2k3*M2*p2k1 + 8*M2
     &**2*p2k1 + 6*k2k3*p1k1*p2k1 - 8*M2*p1k1*p2k1 - 8*k2k3*p1k2*p2k1 - 20*M2*p1k2*p2k1 - 2*p1k1*p1k2
     &*p2k1 + 16*p1k2**2*p2k1 + 9*k2k3*p1k3*p2k1 - 8*M2*p1k3*p2k1 - 12*p1k1*p1k3*p2k1 - 6*p1k2*p1k3
     &*p2k1 - 6*p1k3**2*p2k1 - 6*k2k3*p1p2*p2k1 - 8*M2*p1p2*p2k1 + 16*p1k1*p1p2*p2k1 + 4*p1k2*p1p2
     &*p2k1 + 12*p1k3*p1p2*p2k1 - 12*M2*p2k1**2 + 14*p1k2*p2k1**2 - 2*k2k3**2*p2k2 + 4*k2k3*M2*p2k2 
     &+ 14*k2k3*p1k1*p2k2 + 4*M2*p1k1*p2k2 - 6*p1k1**2*p2k2 + 4*k2k3*p1k2*p2k2 - 16*p1k1*p1k2*p2k2 
     &+ 2*k2k3*p1k3*p2k2 - 16*M2*p1k3*p2k2 + 6*p1k1*p1k3*p2k2 + 16*p1k2*p1k3*p2k2 + 12*k2k3*p1p2
     &*p2k2 - 20*p1k1*p1p2*p2k2 + 8*p1k3*p1p2*p2k2 + 4*k2k3*p2k1*p2k2 + 8*M2*p2k1*p2k2 - 14*p1k1*p2k1
     &*p2k2 - 12*p1k2*p2k1*p2k2 - 8*p1k3*p2k1*p2k2 + 8*p1k1*p2k2**2 - 8*p1k3*p2k2**2 + 6*k2k3*M2
     &*p2k3 - 8*M2**2*p2k3 - 7*k2k3*p1k1*p2k3 + 8*M2*p1k1*p2k3 + 12*p1k1**2*p2k3 - 8*k2k3*p1k2*p2k3 
     &+ 20*M2*p1k2*p2k3 + 4*p1k1*p1k2*p2k3 - 4*p1k2**2*p2k3 - 4*M2*p1k3*p2k3 + 6*p1k1*p1k3*p2k3 + 4
     &*p1k2*p1k3*p2k3 - 10*k2k3*p1p2*p2k3 + 8*M2*p1p2*p2k3 + 20*p1k1*p1p2*p2k3 - 4*p1k2*p1p2*p2k3 
     &+ 12*p1k3*p1p2*p2k3 + k2k3*p2k1*p2k3 - 20*M2*p2k1*p2k3 + 22*p1k2*p2k1*p2k3 - 2*p1k3*p2k1*p2k3 
     &+ 6*k2k3*p2k2*p2k3 - 8*M2*p2k2*p2k3 - 18*p1k1*p2k2*p2k3 - 8*p1k3*p2k2*p2k3 - 4*M2*p2k3**2 + 2
     &*p1k1*p2k3**2 + 6*p1k2*p2k3**2 + 2*k1k2**2*(3*M2 - p1k2 - 3*p1k3 - 3*p1p2 + p2k2 + 3*p2k3) + 2
     &*k1k3*(4*M2**2 - 2*p1k2*p1p2 - 4*p1p2**2 - p1k2*p2k1 + p1k3*p2k1 - 2*p1p2*p2k1 - p1k3*p2k2 
     &+ 2*p1p2*p2k2 + p2k1*p2k2 + 2*M2*(p1k1 - p1k2 + p1k3 + p2k1 + p2k2) - 2*p1p2*p2k3 + p2k2*p2k3 
     &- p1k1*(p1k2 + 2*p1p2 - p2k2 + p2k3)) - k1k2*(-12*M2**2 + 6*p1k1*p1k2 - 4*p1k2**2 + k1k3*p1k3 
     &- 6*p1k1*p1k3 - 10*p1k2*p1k3 - 3*p1k3**2 + 2*p1k1*p1p2 - 4*p1k2*p1p2 - 20*p1k3*p1p2 + 4*p1p2
     &**2 + 4*p1k2*p2k1 - 12*p1k3*p2k1 + 6*p1p2*p2k1 - 4*p1k1*p2k2 + 4*p1k2*p2k2 + 20*p1k3*p2k2 - 6
     &*p2k1*p2k2 + 3*k2k3*(2*M2 - 2*p1k1 + p1k3 - 2*p1p2 + 2*p2k1 - p2k3) - k1k3*p2k3 + 18*p1k1*p2k3 
     &- 2*p1k2*p2k3 + 2*p1k3*p2k3 + 8*p1p2*p2k3 - 8*p2k2*p2k3 + p2k3**2 + 2*M2*(p1k1 + 6*p1k2 + 4*p1p2 
     &- 5*p2k1 - 4*p2k2 + 6*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(128*k1k3*p1k2*(-k1k2 
     &+ p1k1 + p1k2)*(M2 + p1p2)))

!************************************************   2 Re(Msk2 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(2*k1k3**2*(4*M2 + p1k2 - 4
     &*p1p2 - p2k2) + k1k2*(k1k3 - 4*M2 - 2*p1k1 + 4*p1k3)*(2*M2 - p1k1 - p1k3 - 2*p1p2 + p2k1 + p2k3) 
     &+ k1k3*(8*M2**2 + 2*p1k1**2 - 11*p1k1*p1k2 - 4*p1k1*p1k3 - 3*p1k2*p1k3 - 14*p1k3**2 - 18
     &*p1k1*p1p2 - 4*p1k2*p1p2 - 2*p1k3*p1p2 - 8*p1p2**2 - 2*p1k1*p2k1 + 5*p1k2*p2k1 - 6*p1k3*p2k1 
     &- 2*p1p2*p2k1 + 7*p1k1*p2k2 - p1k3*p2k2 - 4*p1p2*p2k2 - p2k1*p2k2 + 10*p1k1*p2k3 + 5*p1k2*p2k3 
     &+ 14*p1k3*p2k3 - 2*p1p2*p2k3 - p2k2*p2k3 + k2k3*(2*M2 - p1k1 - p1k3 - 2*p1p2 + p2k1 + p2k3) 
     &+ 2*M2*(p1k1 + 2*p1k2 + p1k3 + 9*p2k1 + 2*p2k2 + p2k3)) + 2*(-2*p1k1**3 + p1k1**2*p1k2 - 2*p1k1
     &**2*p1k3 + 2*p1k1*p1k2*p1k3 + 2*p1k1*p1k3**2 + 3*p1k2*p1k3**2 + 2*p1k3**3 - 2*p1k1**2*p1p2 
     &- 2*p1k1*p1k2*p1p2 + 4*p1k1*p1k3*p1p2 + 2*p1k2*p1k3*p1p2 + 10*p1k3**2*p1p2 - 4*p1k1*p1p2**2 
     &+ 4*p1k3*p1p2**2 + 2*p1k1**2*p2k1 + p1k1*p1k2*p2k1 - 6*p1k1*p1k3*p2k1 - 12*p1k3**2*p2k1 + 2*p1k1
     &*p1p2*p2k1 - 4*p1k3*p1p2*p2k1 - 3*p1k1**2*p2k2 - 2*p1k1*p1k3*p2k2 - p1k3**2*p2k2 - 2*p1k1
     &*p1p2*p2k2 + 2*p1k3*p1p2*p2k2 + p1k1*p2k1*p2k2 - 2*p1k3*p2k1*p2k2 - 4*M2**2*(p2k1 - p2k3) + 8
     &*p1k1**2*p2k3 - 2*p1k1*p1k2*p2k3 + 10*p1k1*p1k3*p2k3 - 5*p1k2*p1k3*p2k3 - 2*p1k3**2*p2k3 + 8*p1k1
     &*p1p2*p2k3 + 6*p1k3*p1p2*p2k3 + 4*p1k1*p2k2*p2k3 + 3*p1k3*p2k2*p2k3 + k2k3*(2*M2 - 4*p1k1 
     &- 3*p1k3)*(2*M2 - p1k1 - p1k3 - 2*p1p2 + p2k1 + p2k3) - 2*M2*(p1k1**2 + 3*p1k3**2 + 2*p1k3*p1p2 
     &- 4*p1k3*p2k1 - 2*p1p2*p2k1 - p2k1*p2k2 - p1k1*(2*p1k2 - 2*p1k3 + 2*p1p2 + p2k1 - 6*p2k3) 
     &+ p1k2*(2*p1k3 + p2k1 - p2k3) + 5*p1k3*p2k3 + 2*p1p2*p2k3 + p2k2*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)
     &/(128*k1k3*(k1k3 - p1k1 - p1k3)*p1k3*(M2 + p1p2)))

!************************************************   2 Re(Msk3 Mvk1^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Msk3 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Msk3 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Msk3 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(2*k1k2**2*(4*M2 + p1k3 - 4
     &*p1p2 - p2k3) + k1k2*(8*M2**2 + 2*M2*p1k1 - 14*p1k1**2 + 2*M2*p1k2 - 4*p1k1*p1k2 + 2*p1k2**2 
     &+ 4*M2*p1k3 - 3*p1k1*p1k3 - 11*p1k2*p1k3 - 2*p1k1*p1p2 - 18*p1k2*p1p2 - 4*p1k3*p1p2 - 8*p1p2
     &**2 + 2*M2*p2k1 + 14*p1k1*p2k1 + 10*p1k2*p2k1 + 5*p1k3*p2k1 - 2*p1p2*p2k1 + 18*M2*p2k2 - 6*p1k1
     &*p2k2 - 2*p1k2*p2k2 + 5*p1k3*p2k2 - 2*p1p2*p2k2 + k1k3*(2*M2 - p1k1 - p1k2 - 2*p1p2 + p2k1 
     &+ p2k2) + k2k3*(2*M2 - p1k1 - p1k2 - 2*p1p2 + p2k1 + p2k2) + 4*M2*p2k3 - p1k1*p2k3 + 7*p1k2*p2k3 
     &- 4*p1p2*p2k3 - p2k1*p2k3 - p2k2*p2k3) + 2*(-6*M2*p1k1**2 + 2*p1k1**3 - 4*M2*p1k1*p1k2 + 2
     &*p1k1**2*p1k2 - 2*M2*p1k2**2 - 2*p1k1*p1k2**2 - 2*p1k2**3 - 4*M2*p1k1*p1k3 + 3*p1k1**2*p1k3 
     &+ 4*M2*p1k2*p1k3 + 2*p1k1*p1k2*p1k3 + p1k2**2*p1k3 - 4*M2*p1k1*p1p2 + 10*p1k1**2*p1p2 + 4*M2
     &*p1k2*p1p2 + 4*p1k1*p1k2*p1p2 - 2*p1k2**2*p1p2 + 2*p1k1*p1k3*p1p2 - 2*p1k2*p1k3*p1p2 + 4*p1k1
     &*p1p2**2 - 4*p1k2*p1p2**2 + 4*M2**2*p2k1 - 10*M2*p1k1*p2k1 - 2*p1k1**2*p2k1 - 12*M2*p1k2*p2k1 
     &+ 10*p1k1*p1k2*p2k1 + 8*p1k2**2*p2k1 + 2*M2*p1k3*p2k1 - 5*p1k1*p1k3*p2k1 - 2*p1k2*p1k3*p2k1 
     &- 4*M2*p1p2*p2k1 + 6*p1k1*p1p2*p2k1 + 8*p1k2*p1p2*p2k1 - 4*M2**2*p2k2 + 8*M2*p1k1*p2k2 - 12*p1k1
     &**2*p2k2 + 2*M2*p1k2*p2k2 - 6*p1k1*p1k2*p2k2 + 2*p1k2**2*p2k2 - 2*M2*p1k3*p2k2 + p1k2*p1k3
     &*p2k2 + 4*M2*p1p2*p2k2 - 4*p1k1*p1p2*p2k2 + 2*p1k2*p1p2*p2k2 + k1k3*(2*M2 - 3*p1k1 - 4*p1k2)*(2
     &*M2 - p1k1 - p1k2 - 2*p1p2 + p2k1 + p2k2) - k2k3*(2*M2 - 2*p1k1 + p1k2)*(2*M2 - p1k1 - p1k2 
     &- 2*p1p2 + p2k1 + p2k2) - p1k1**2*p2k3 - 2*p1k1*p1k2*p2k3 - 3*p1k2**2*p2k3 + 2*p1k1*p1p2*p2k3 
     &- 2*p1k2*p1p2*p2k3 - 2*M2*p2k1*p2k3 + 3*p1k1*p2k1*p2k3 + 4*p1k2*p2k1*p2k3 + 2*M2*p2k2*p2k3 - 2*p1k1*p2k2*p2k3 
     &+ p1k2*p2k2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(128*k1k2*p1k1*(k1k2 - p1k1 - p1k2)*(M2 + p1p2)))

!************************************************   2 Re(Msk3 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Msk3 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG**2*Nc**2*(-1 + Nc**2)*(-6*k2k3**2*M2 - 12*k2k3
     &*M2**2 + 6*k2k3**2*p1k1 + 16*M2**2*p1k1 - 3*k2k3*p1k1**2 + 8*M2*p1k1**2 + 2*k2k3*M2*p1k2 + 8*M2
     &**2*p1k2 - 6*k2k3*p1k1*p1k2 + 12*M2*p1k1*p1k2 + 4*M2*p1k2**2 + 2*k2k3**2*p1k3 + 12*k2k3*M2*p1k3 
     &- 10*k2k3*p1k1*p1k3 - 12*M2*p1k1*p1k3 + 2*p1k1**2*p1k3 + 6*k2k3*p1k2*p1k3 - 8*p1k2**2*p1k3 
     &- 4*k2k3*p1k3**2 + 4*p1k1*p1k3**2 - 4*p1k2*p1k3**2 + 6*k2k3**2*p1p2 + 8*k2k3*M2*p1p2 - 20*k2k3
     &*p1k1*p1p2 - 8*M2*p1k1*p1p2 - 4*p1k1**2*p1p2 + 2*k2k3*p1k2*p1p2 + 8*M2*p1k2*p1p2 - 8*p1k2
     &**2*p1p2 - 4*k2k3*p1k3*p1p2 + 12*p1k1*p1k3*p1p2 - 24*p1k2*p1k3*p1p2 + 4*k2k3*p1p2**2 - 8*p1k1
     &*p1p2**2 - 16*p1k2*p1p2**2 - 6*k2k3**2*p2k1 + 12*k2k3*M2*p2k1 + 8*M2**2*p2k1 + 2*k2k3*p1k1*p2k1 
     &+ 4*M2*p1k1*p2k1 + 18*k2k3*p1k2*p2k1 - 8*M2*p1k2*p2k1 - 6*p1k1*p1k2*p2k1 - 12*p1k2**2*p2k1 
     &- 2*k2k3*p1k3*p2k1 - 20*M2*p1k3*p2k1 - 4*p1k1*p1k3*p2k1 - 4*p1k2*p1k3*p2k1 + 4*p1k3**2*p2k1 
     &+ 8*k2k3*p1p2*p2k1 - 8*M2*p1p2*p2k1 - 12*p1k1*p1p2*p2k1 - 20*p1k2*p1p2*p2k1 + 4*p1k3*p1p2*p2k1 
     &+ k2k3*p2k1**2 + 4*M2*p2k1**2 - 2*p1k2*p2k1**2 - 6*p1k3*p2k1**2 - 10*k2k3*M2*p2k2 - 8*M2**2
     &*p2k2 - 12*k2k3*p1k1*p2k2 + 8*M2*p1k1*p2k2 + 6*p1k1**2*p2k2 + 8*M2*p1k2*p2k2 + 12*p1k1*p1k2*p2k2 
     &+ 4*k2k3*p1k3*p2k2 + 20*M2*p1k3*p2k2 + 6*p1k1*p1k3*p2k2 + 2*p1k2*p1k3*p2k2 - 16*p1k3**2*p2k2 
     &+ 6*k2k3*p1p2*p2k2 + 8*M2*p1p2*p2k2 - 12*p1k1*p1p2*p2k2 - 16*p1k2*p1p2*p2k2 - 4*p1k3*p1p2
     &*p2k2 + 20*M2*p2k1*p2k2 + 2*p1k1*p2k1*p2k2 - 22*p1k3*p2k1*p2k2 + 12*M2*p2k2**2 - 14*p1k3*p2k2
     &**2 - 2*k2k3**2*p2k3 - 8*k2k3*M2*p2k3 + 20*k2k3*p1k1*p2k3 + 16*M2*p1k1*p2k3 - 4*k2k3*p1k2*p2k3 
     &- 4*M2*p1k2*p2k3 - 6*p1k1*p1k2*p2k3 + 6*p1k2**2*p2k3 + 4*k2k3*p1k3*p2k3 - 16*p1k1*p1k3*p2k3 
     &+ 16*p1k2*p1k3*p2k3 - 8*p1k1*p1p2*p2k3 + 20*p1k2*p1p2*p2k3 - 8*k2k3*p2k1*p2k3 + 8*M2*p2k1*p2k3 
     &+ 8*p1k1*p2k1*p2k3 + 18*p1k2*p2k1*p2k3 - 6*k2k3*p2k2*p2k3 - 8*M2*p2k2*p2k3 + 8*p1k1*p2k2*p2k3 
     &+ 14*p1k2*p2k2*p2k3 + 12*p1k3*p2k2*p2k3 + 8*p1k1*p2k3**2 - 8*p1k2*p2k3**2 + k1k3**2*(6*M2 
     &- 3*p1k2 - 2*p1k3 - 6*p1p2 + 3*p2k2 + 2*p2k3) + k1k2*(-8*M2**2 + k2k3*p1k1 + 2*p1k2*p1k3 + 4
     &*p1k2*p1p2 + 4*p1k3*p1p2 + 8*p1p2**2 - k2k3*p2k1 + 2*p1k2*p2k1 + 4*p1p2*p2k1 - 2*p1k1*p2k2 + 2
     &*p1k3*p2k2 + 4*p1p2*p2k2 + k1k3*(2*M2 - p1k2 - 2*p1p2 + p2k2) + 2*p1k1*p2k3 - 2*p1k2*p2k3 - 4
     &*p1p2*p2k3 - 2*p2k1*p2k3 - 2*p2k2*p2k3 - 4*M2*(p1k1 + p1k2 - p1k3 + p2k2 + p2k3)) + k1k3*(-12
     &*M2**2 + 3*p1k1*p1k2 + 6*p1k2**2 + 10*p1k2*p1k3 + 4*p1k3**2 + 10*p1k1*p1p2 + 26*p1k2*p1p2 + 16
     &*p1k3*p1p2 + 20*p1p2**2 + 7*p1k2*p2k1 + 8*p1k3*p2k1 + 10*p1p2*p2k1 - 9*p1k1*p2k2 - 6*p1k2*p2k2 
     &+ 8*p1k3*p2k2 + 6*p1p2*p2k2 - p2k1*p2k2 + 3*k2k3*(2*M2 + p1k1 - 2*p1k2 - 2*p1p2 - p2k1 + 2
     &*p2k2) - 2*p1k1*p2k3 - 14*p1k2*p2k3 - 4*p1k3*p2k3 - 12*p1p2*p2k3 - 6*p2k1*p2k3 - 4*p2k2*p2k3 
     &- 2*M2*(7*p1k1 + 9*p1k2 + 4*p1p2 + 3*p2k1 + 7*p2k2 + 2*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)
     &/(128*k1k2*p1k3*(-k2k3 + p1k2 + p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Msk3 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(0)

!************************************************   2 Re(Msk3 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(2*k1k2**2*(4*M2 + p1k3 - 4
     &*p1p2 - p2k3) + k1k2*(8*M2**2 + 2*M2*p1k1 + 2*p1k1**2 + 2*M2*p1k2 - 4*p1k1*p1k2 - 14*p1k2**2 
     &+ 4*M2*p1k3 - 11*p1k1*p1k3 - 3*p1k2*p1k3 - 18*p1k1*p1p2 - 2*p1k2*p1p2 - 4*p1k3*p1p2 - 8*p1p2
     &**2 + 18*M2*p2k1 - 2*p1k1*p2k1 - 6*p1k2*p2k1 + 5*p1k3*p2k1 - 2*p1p2*p2k1 + 2*M2*p2k2 + 10*p1k1
     &*p2k2 + 14*p1k2*p2k2 + 5*p1k3*p2k2 - 2*p1p2*p2k2 + k1k3*(2*M2 - p1k1 - p1k2 - 2*p1p2 + p2k1 
     &+ p2k2) + k2k3*(2*M2 - p1k1 - p1k2 - 2*p1p2 + p2k1 + p2k2) + 4*M2*p2k3 + 7*p1k1*p2k3 - p1k2*p2k3 
     &- 4*p1p2*p2k3 - p2k1*p2k3 - p2k2*p2k3) - 2*(2*M2*p1k1**2 + 2*p1k1**3 + 4*M2*p1k1*p1k2 + 2
     &*p1k1**2*p1k2 + 6*M2*p1k2**2 - 2*p1k1*p1k2**2 - 2*p1k2**3 - 4*M2*p1k1*p1k3 - p1k1**2*p1k3 + 4
     &*M2*p1k2*p1k3 - 2*p1k1*p1k2*p1k3 - 3*p1k2**2*p1k3 - 4*M2*p1k1*p1p2 + 2*p1k1**2*p1p2 + 4*M2*p1k2
     &*p1p2 - 4*p1k1*p1k2*p1p2 - 10*p1k2**2*p1p2 + 2*p1k1*p1k3*p1p2 - 2*p1k2*p1k3*p1p2 + 4*p1k1*p1p2
     &**2 - 4*p1k2*p1p2**2 + 4*M2**2*p2k1 - 2*M2*p1k1*p2k1 - 2*p1k1**2*p2k1 - 8*M2*p1k2*p2k1 + 6
     &*p1k1*p1k2*p2k1 + 12*p1k2**2*p2k1 + 2*M2*p1k3*p2k1 - p1k1*p1k3*p2k1 - 4*M2*p1p2*p2k1 - 2*p1k1
     &*p1p2*p2k1 + 4*p1k2*p1p2*p2k1 - 4*M2**2*p2k2 + 12*M2*p1k1*p2k2 - 8*p1k1**2*p2k2 + 10*M2*p1k2
     &*p2k2 - 10*p1k1*p1k2*p2k2 + 2*p1k2**2*p2k2 - 2*M2*p1k3*p2k2 + 2*p1k1*p1k3*p2k2 + 5*p1k2*p1k3*p2k2 
     &+ 4*M2*p1p2*p2k2 - 8*p1k1*p1p2*p2k2 - 6*p1k2*p1p2*p2k2 - k2k3*(2*M2 - 4*p1k1 - 3*p1k2)*(2
     &*M2 - p1k1 - p1k2 - 2*p1p2 + p2k1 + p2k2) + k1k3*(2*M2 + p1k1 - 2*p1k2)*(2*M2 - p1k1 - p1k2 
     &- 2*p1p2 + p2k1 + p2k2) + 3*p1k1**2*p2k3 + 2*p1k1*p1k2*p2k3 + p1k2**2*p2k3 + 2*p1k1*p1p2*p2k3 
     &- 2*p1k2*p1p2*p2k3 - 2*M2*p2k1*p2k3 - p1k1*p2k1*p2k3 + 2*p1k2*p2k1*p2k3 + 2*M2*p2k2*p2k3 
     &- 4*p1k1*p2k2*p2k3 - 3*p1k2*p2k2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(128*k1k2*(k1k2 - p1k1 - p1k2)*p1k2*(M2 + p1p2)))

!************************************************   2 Re(Msk3 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG**2*Nc**2*(-1 + Nc**2)*(-6*k2k3**2*M2 + 12*k2k3*M2
     &**2 + 3*k2k3**2*p1k1 + 18*k2k3*M2*p1k1 - 8*M2**2*p1k1 - 6*k2k3*p1k1**2 - 4*M2*p1k1**2 + 14*k2k3
     &*M2*p1k2 - 16*M2**2*p1k2 - 3*k2k3*p1k1*p1k2 - 12*M2*p1k1*p1k2 - 8*M2*p1k2**2 + 2*k2k3**2*p1k3 
     &- 10*k2k3*p1k1*p1k3 + 8*p1k1**2*p1k3 + 12*M2*p1k2*p1k3 - 2*p1k2**2*p1k3 - 4*k2k3*p1k3**2 
     &+ 4*p1k1*p1k3**2 - 4*p1k2*p1k3**2 + 6*k2k3**2*p1p2 + 8*k2k3*M2*p1p2 - 26*k2k3*p1k1*p1p2 - 8*M2
     &*p1k1*p1p2 + 8*p1k1**2*p1p2 - 10*k2k3*p1k2*p1p2 + 8*M2*p1k2*p1p2 + 4*p1k2**2*p1p2 - 16*k2k3
     &*p1k3*p1p2 + 24*p1k1*p1k3*p1p2 - 12*p1k2*p1k3*p1p2 - 20*k2k3*p1p2**2 + 16*p1k1*p1p2**2 + 8*p1k2
     &*p1p2**2 + k1k2*k2k3*(-2*M2 + p1k1 + 2*p1p2 - p2k1) - 3*k2k3**2*p2k1 + 14*k2k3*M2*p2k1 + 8*M2
     &**2*p2k1 + 6*k2k3*p1k1*p2k1 - 8*M2*p1k1*p2k1 + 9*k2k3*p1k2*p2k1 - 8*M2*p1k2*p2k1 - 12*p1k1*p1k2
     &*p2k1 - 6*p1k2**2*p2k1 - 8*k2k3*p1k3*p2k1 - 20*M2*p1k3*p2k1 - 2*p1k1*p1k3*p2k1 - 6*p1k2*p1k3
     &*p2k1 + 16*p1k3**2*p2k1 - 6*k2k3*p1p2*p2k1 - 8*M2*p1p2*p2k1 + 16*p1k1*p1p2*p2k1 + 12*p1k2*p1p2
     &*p2k1 + 4*p1k3*p1p2*p2k1 - 12*M2*p2k1**2 + 14*p1k3*p2k1**2 + 6*k2k3*M2*p2k2 - 8*M2**2*p2k2 
     &- 7*k2k3*p1k1*p2k2 + 8*M2*p1k1*p2k2 + 12*p1k1**2*p2k2 - 4*M2*p1k2*p2k2 + 6*p1k1*p1k2*p2k2 - 8
     &*k2k3*p1k3*p2k2 + 20*M2*p1k3*p2k2 + 4*p1k1*p1k3*p2k2 + 4*p1k2*p1k3*p2k2 - 4*p1k3**2*p2k2 - 10
     &*k2k3*p1p2*p2k2 + 8*M2*p1p2*p2k2 + 20*p1k1*p1p2*p2k2 + 12*p1k2*p1p2*p2k2 - 4*p1k3*p1p2*p2k2 
     &+ k2k3*p2k1*p2k2 - 20*M2*p2k1*p2k2 - 2*p1k2*p2k1*p2k2 + 22*p1k3*p2k1*p2k2 - 4*M2*p2k2**2 + 2
     &*p1k1*p2k2**2 + 6*p1k3*p2k2**2 - 2*k2k3**2*p2k3 + 4*k2k3*M2*p2k3 + 14*k2k3*p1k1*p2k3 + 4*M2*p1k1
     &*p2k3 - 6*p1k1**2*p2k3 + 2*k2k3*p1k2*p2k3 - 16*M2*p1k2*p2k3 + 6*p1k1*p1k2*p2k3 + 4*k2k3*p1k3
     &*p2k3 - 16*p1k1*p1k3*p2k3 + 16*p1k2*p1k3*p2k3 + 12*k2k3*p1p2*p2k3 - 20*p1k1*p1p2*p2k3 + 8*p1k2
     &*p1p2*p2k3 + 4*k2k3*p2k1*p2k3 + 8*M2*p2k1*p2k3 - 14*p1k1*p2k1*p2k3 - 8*p1k2*p2k1*p2k3 - 12
     &*p1k3*p2k1*p2k3 + 6*k2k3*p2k2*p2k3 - 8*M2*p2k2*p2k3 - 18*p1k1*p2k2*p2k3 - 8*p1k2*p2k2*p2k3 + 8
     &*p1k1*p2k3**2 - 8*p1k2*p2k3**2 + 2*k1k3**2*(3*M2 - 3*p1k2 - p1k3 - 3*p1p2 + 3*p2k2 + p2k3) - k1k3
     &*(-12*M2**2 + k1k2*p1k2 - 6*p1k1*p1k2 - 3*p1k2**2 + 6*p1k1*p1k3 - 10*p1k2*p1k3 - 4*p1k3*
     &*2 + 2*p1k1*p1p2 - 20*p1k2*p1p2 - 4*p1k3*p1p2 + 4*p1p2**2 - 12*p1k2*p2k1 + 4*p1k3*p2k1 + 6*p1p2
     &*p2k1 + 3*k2k3*(2*M2 - 2*p1k1 + p1k2 - 2*p1p2 + 2*p2k1 - p2k2) - k1k2*p2k2 + 18*p1k1*p2k2 + 2
     &*p1k2*p2k2 - 2*p1k3*p2k2 + 8*p1p2*p2k2 + p2k2**2 + 2*M2*(p1k1 + 6*p1k3 + 4*p1p2 - 5*p2k1 + 6
     &*p2k2 - 4*p2k3) - 4*p1k1*p2k3 + 20*p1k2*p2k3 + 4*p1k3*p2k3 - 6*p2k1*p2k3 - 8*p2k2*p2k3) + 2*k1k2
     &*(4*M2**2 - 2*p1k3*p1p2 - 4*p1p2**2 + p1k2*p2k1 - p1k3*p2k1 - 2*p1p2*p2k1 - 2*p1p2*p2k2 - p1k1
     &*(p1k3 + 2*p1p2 + p2k2 - p2k3) - p1k2*p2k3 + 2*p1p2*p2k3 + p2k1*p2k3 + p2k2*p2k3 + 2*M2*(p1k1 + p1k2 - p1k3 + p2k1 + p2k3)))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(128*k1k2*p1k3*(-k1k3 + p1k1 + p1k3)*(M2 + p1p2)))

!************************************************   2 Re(Mvk1 Mvk2^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((9*GGG**2*(2*k1k2 + k1k3 + k2k3)*Nc**2*(-1 + Nc**2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(16*k1k3*k2k3))

!************************************************   2 Re(Mvk1 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((9*GGG**2*(k1k2 + 2*k1k3 + k2k3)*Nc**2*(-1 + Nc**2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(16*k1k2*k2k3))

!************************************************   2 Re(Mvk1 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(2 - 3*Nc**2 + Nc**4)*(3*k1k2**2 + k1k3*(2*M2 - 2
     &*p1k1 - 3*p1k2 - 6*p1p2 + 2*p2k1 - p2k2) + k1k2*(k2k3 - 2*M2 + 2*p1k1 - 9*p1k2 + 3*p1k3 - 6
     &*p1p2 + 4*p2k1 - p2k2 - p2k3) + 2*(3*p1k2**2 + 4*p1k2*p1p2 + 2*p1k3*p1p2 - 3*p1k2*p2k1 + k2k3
     &*(-p1k2 - 2*p1p2 + p2k1) - 4*M2*p2k2 + p1k2*p2k2 + p1k3*p2k2 - 2*M2*p2k3 + p1k2*p2k3 + p1k1*(p1k2 
     &- p1k3 + 6*p2k2 + 3*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*k2k3*p1k1*(-k1k2 + p1k1 + p1k2)))

!************************************************   2 Re(Mvk1 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(2 - 3*Nc**2 + Nc**4)*(3*k1k3**2 + k1k2*(2*M2 - 2
     &*p1k1 - 3*p1k3 - 6*p1p2 + 2*p2k1 - p2k3) + k1k3*(k2k3 - 2*M2 + 2*p1k1 + 3*p1k2 - 9*p1k3 - 6
     &*p1p2 + 4*p2k1 - p2k2 - p2k3) + 2*(3*p1k3**2 + 2*p1k2*p1p2 + 4*p1k3*p1p2 - 3*p1k3*p2k1 + k2k3
     &*(-p1k3 - 2*p1p2 + p2k1) - 2*M2*p2k2 + p1k3*p2k2 - 4*M2*p2k3 + p1k2*p2k3 + p1k3*p2k3 + p1k1*(
     &-p1k2 + p1k3 + 3*p2k2 + 6*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*k2k3*p1k1*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mvk1 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(2 - 3*Nc**2 + Nc**4)*(2*k2k3**2 + k2k3*(-4*M2 
     &- 7*p1k2 + p1k3 + 4*p1p2 + p2k2 + p2k3) + 2*(3*p1k2**2 + p1k3**2 - 2*p1k3*p1p2 - 2*M2*p2k2 + 2
     &*p1k3*p2k2 + p1k2*(2*p1k3 + 2*p1p2 - p2k2 - 4*p2k3) + 2*M2*p2k3 - 3*p1k3*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k2k3*(k2k3 - p1k2 - p1k3)*p1k3))

!************************************************   2 Re(Mvk1 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(2 - 3*Nc**2 + Nc**4)*(2*k2k3**2 + k2k3*(-4*M2 
     &+ p1k2 - 7*p1k3 + 4*p1p2 + p2k2 + p2k3) + 2*(p1k2**2 + 3*p1k3**2 + 2*p1k3*p1p2 + 2*M2*p2k2 - 4
     &*p1k3*p2k2 - 2*M2*p2k3 - p1k3*p2k3 + p1k2*(2*p1k3 - 2*p1p2 - 3*p2k2 + 2*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k2k3*p1k2*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(Mvk1 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2**2 + 4*k2k3*M2 - 4*k2k3*p1k2 
     &- 6*p1k1*p1k2 - 2*p1k2**2 - 6*p1k1*p1k3 - 4*p1k2*p1k3 - 8*p1k2*p1p2 - 4*p1k3*p1p2 - 2*k2k3*p2k1 
     &+ 2*p1k2*p2k1 + 4*p1k3*p2k1 - 2*k2k3*p2k2 + 8*M2*p2k2 - 4*p1k1*p2k2 - 6*p1k2*p2k2 + 6*p1k3
     &*p2k2 + k1k3*(6*M2 - 4*p1k1 - 9*p1k2 - 2*p1p2 + p2k2) + 4*M2*p2k3 - 2*p1k1*p2k3 - 6*p1k2*p2k3 
     &+ k1k2*(2*k1k3 - 3*k2k3 + 6*M2 - 2*p1k1 - 3*p1k2 + 9*p1k3 + 2*p1p2 + p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k2k3*p1k2*(-k1k2 + p1k1 + p1k2)))

!************************************************   2 Re(Mvk1 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k3**2 + 4*k2k3*M2 - 6*p1k1*p1k2 
     &- 4*k2k3*p1k3 - 6*p1k1*p1k3 - 4*p1k2*p1k3 - 2*p1k3**2 - 4*p1k2*p1p2 - 8*p1k3*p1p2 - 2*k2k3*p2k1 
     &+ 4*p1k2*p2k1 + 2*p1k3*p2k1 + 4*M2*p2k2 - 2*p1k1*p2k2 - 6*p1k3*p2k2 - 2*k2k3*p2k3 + 8*M2*p2k3 
     &- 4*p1k1*p2k3 + 6*p1k2*p2k3 - 6*p1k3*p2k3 + k1k2*(2*k1k3 + 6*M2 - 4*p1k1 - 9*p1k3 - 2*p1p2 
     &+ p2k3) + k1k3*(-3*k2k3 + 6*M2 - 2*p1k1 + 9*p1k2 - 3*p1k3 + 2*p1p2 + p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k2k3*p1k3*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mvk2 Mvk3^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((9*GGG**2*(k1k2 + k1k3 + 2*k2k3)*Nc**2*(-1 + Nc**2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(16*k1k2*k1k3))

!************************************************   2 Re(Mvk2 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(k1k2**2 + 4*k1k3*M2 + 6*k2k3*M2 - 9
     &*k2k3*p1k1 - 2*p1k1**2 - 4*k2k3*p1k2 - 6*p1k1*p1k2 - 4*p1k1*p1k3 - 6*p1k2*p1k3 - 2*k2k3*p1p2 
     &- 8*p1k1*p1p2 - 4*p1k3*p1p2 + k2k3*p2k1 + 8*M2*p2k1 - 6*p1k1*p2k1 - 4*p1k2*p2k1 + 6*p1k3*p2k1 
     &+ 2*p1k1*p2k2 + 4*p1k3*p2k2 - 2*k1k3*(2*p1k1 + p2k1 + p2k2) + 4*M2*p2k3 - 6*p1k1*p2k3 - 2*p1k2
     &*p2k3 + k1k2*(-3*k1k3 + 2*k2k3 + 6*M2 - 3*p1k1 - 2*p1k2 + 9*p1k3 + 2*p1p2 + p2k1 + p2k3))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k1k3*p1k1*(-k1k2 + p1k1 + p1k2)))

!************************************************   2 Re(Mvk2 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(2 - 3*Nc**2 + Nc**4)*(2*k1k3**2 + k1k3*(-4*M2 
     &+ p1k1 - 7*p1k3 + 4*p1p2 + p2k1 + p2k3) + 2*(p1k1**2 + 3*p1k3**2 + 2*p1k3*p1p2 + 2*M2*p2k1 - 4
     &*p1k3*p2k1 - 2*M2*p2k3 - p1k3*p2k3 + p1k1*(2*p1k3 - 2*p1p2 - 3*p2k1 + 2*p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k3*p1k1*(k1k3 - p1k1 - p1k3)))

!************************************************   2 Re(Mvk2 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(-1 + Nc**2)*(-3*k1k3*k2k3 + k2k3**2 + 6*k2k3*M2 
     &+ 9*k2k3*p1k1 - 2*k2k3*p1k2 - 6*p1k1*p1k2 - 3*k2k3*p1k3 - 4*p1k1*p1k3 - 6*p1k2*p1k3 - 2*p1k3
     &**2 + 2*k2k3*p1p2 - 4*p1k1*p1p2 - 8*p1k3*p1p2 + k2k3*p2k1 + 4*M2*p2k1 - 2*p1k2*p2k1 - 6*p1k3
     &*p2k1 + 4*p1k1*p2k2 + 2*p1k3*p2k2 + k2k3*p2k3 + 8*M2*p2k3 + 6*p1k1*p2k3 - 4*p1k2*p2k3 - 6*p1k3
     &*p2k3 + k1k2*(2*k2k3 + 6*M2 - 4*p1k2 - 9*p1k3 - 2*p1p2 + p2k3) - 2*k1k3*(-2*M2 + 2*p1k3 + p2k2 
     &+ p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k1k3*p1k3*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mvk2 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(2 - 3*Nc**2 + Nc**4)*(3*k2k3**2 + 2*k1k2*M2 - 2
     &*k1k2*p1k2 - 2*p1k1*p1k2 - 3*k1k2*p1k3 + 2*p1k2*p1k3 + 6*p1k3**2 - 6*k1k2*p1p2 + 4*p1k1*p1p2 
     &+ 8*p1k3*p1p2 - 4*M2*p2k1 + 6*p1k2*p2k1 + 2*p1k3*p2k1 + 2*k1k2*p2k2 - 6*p1k3*p2k2 + k1k3*(k2k3 
     &- 2*p1k3 - 4*p1p2 + 2*p2k2) - k1k2*p2k3 - 8*M2*p2k3 + 2*p1k1*p2k3 + 12*p1k2*p2k3 + 2*p1k3*p2k3 
     &- k2k3*(2*M2 - 3*p1k1 - 2*p1k2 + 9*p1k3 + 6*p1p2 + p2k1 - 4*p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k3*p1k2*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mvk2 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(2 - 3*Nc**2 + Nc**4)*(3*k1k2**2 + k2k3*(2*M2 - 3
     &*p1k1 - 2*p1k2 - 6*p1p2 - p2k1 + 2*p2k2) + k1k2*(k1k3 - 2*M2 - 9*p1k1 + 2*p1k2 + 3*p1k3 - 6
     &*p1p2 - p2k1 + 4*p2k2 - p2k3) + 2*(3*p1k1**2 - p1k2*p1k3 + 2*p1k3*p1p2 - 4*M2*p2k1 + 6*p1k2*p2k1 
     &+ p1k3*p2k1 + k1k3*(-p1k1 - 2*p1p2 + p2k2) - 2*M2*p2k3 + 3*p1k2*p2k3 + p1k1*(p1k2 + 4*p1p2 
     &+ p2k1 - 3*p2k2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64.
     &*k1k3*p1k2*(-k1k2 + p1k1 + p1k2)))

!************************************************   2 Re(Mvk2 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(2 - 3*Nc**2 + Nc**4)*(2*k1k3**2 + k1k3*(-4*M2 
     &- 7*p1k1 + p1k3 + 4*p1p2 + p2k1 + p2k3) + 2*(3*p1k1**2 + p1k3**2 - 2*p1k3*p1p2 - 2*M2*p2k1 + 2
     &*p1k3*p2k1 + p1k1*(2*p1k3 + 2*p1p2 - p2k1 - 4*p2k3) + 2*M2*p2k3 - 3*p1k3*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k3*(k1k3 - p1k1 - p1k3)*p1k3))

!************************************************   2 Re(Mvk3 MtP123^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(2 - 3*Nc**2 + Nc**4)*(2*k1k2**2 + k1k2*(-4*M2 
     &+ p1k1 - 7*p1k2 + 4*p1p2 + p2k1 + p2k2) + 2*(p1k1**2 + 3*p1k2**2 + 2*p1k2*p1p2 + 2*M2*p2k1 - 4
     &*p1k2*p2k1 - 2*M2*p2k2 - p1k2*p2k2 + p1k1*(2*p1k2 - 2*p1p2 - 3*p2k1 + 2*p2k2)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k2*p1k1*(k1k2 - p1k1 - p1k2)))

!************************************************   2 Re(Mvk3 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(-1 + Nc**2)*(3*k1k2*k1k3 - k1k3**2 - 6*k2k3*M2 
     &+ 9*k2k3*p1k1 + 2*p1k1**2 + 4*p1k1*p1k2 + 4*k2k3*p1k3 + 6*p1k1*p1k3 + 6*p1k2*p1k3 + 2*k2k3*p1p2 
     &+ 8*p1k1*p1p2 + 4*p1k2*p1p2 - k2k3*p2k1 - 8*M2*p2k1 + 6*p1k1*p2k1 - 6*p1k2*p2k1 + 4*p1k3
     &*p2k1 - 4*M2*p2k2 + 6*p1k1*p2k2 + 2*p1k3*p2k2 - k1k3*(2*k2k3 + 6*M2 - 3*p1k1 + 9*p1k2 - 2*p1k3 
     &+ 2*p1p2 + p2k1 + p2k2) - 2*p1k1*p2k3 - 4*p1k2*p2k3 + 2*k1k2*(-2*M2 + 2*p1k1 + p2k1 + p2k3))
     &*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k1k2*p1k1*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(Mvk3 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(2 - 3*Nc**2 + Nc**4)*(3*k2k3**2 + 2*k1k3*M2 - 3
     &*k1k3*p1k2 + 6*p1k2**2 - 2*k1k3*p1k3 - 2*p1k1*p1k3 + 2*p1k2*p1k3 - 6*k1k3*p1p2 + 4*p1k1*p1p2 
     &+ 8*p1k2*p1p2 - 4*M2*p2k1 + 2*p1k2*p2k1 + 6*p1k3*p2k1 - k1k3*p2k2 - 8*M2*p2k2 + 2*p1k1*p2k2 
     &+ 2*p1k2*p2k2 + 12*p1k3*p2k2 - k2k3*(2*M2 - 3*p1k1 + 9*p1k2 - 2*p1k3 + 6*p1p2 + p2k1 + p2k2 - 4
     &*p2k3) + 2*k1k3*p2k3 - 6*p1k2*p2k3 + k1k2*(k2k3 - 2*p1k2 - 4*p1p2 + 2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k2*p1k3*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mvk3 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(-1 + Nc**2)*(3*k1k2*k2k3 - k2k3**2 - 6*k2k3*M2 
     &- 9*k2k3*p1k1 + 3*k2k3*p1k2 + 4*p1k1*p1k2 + 2*p1k2**2 + 2*k2k3*p1k3 + 6*p1k1*p1k3 + 6*p1k2*p1k3 
     &- 2*k2k3*p1p2 + 4*p1k1*p1p2 + 8*p1k2*p1p2 - k2k3*p2k1 - 4*M2*p2k1 + 6*p1k2*p2k1 + 2*p1k3
     &*p2k1 - k2k3*p2k2 - 8*M2*p2k2 - 6*p1k1*p2k2 + 6*p1k2*p2k2 + 4*p1k3*p2k2 - k1k3*(2*k2k3 + 6*M2 
     &- 9*p1k2 - 4*p1k3 - 2*p1p2 + p2k2) - 4*p1k1*p2k3 - 2*p1k2*p2k3 + 2*k1k2*(-2*M2 + 2*p1k2 + p2k2 
     &+ p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k1k2*p1k2*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(Mvk3 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-(GGG*(2 - 3*Nc**2 + Nc**4)*(2*k1k2**2 + k1k2*(-4*M2 
     &- 7*p1k1 + p1k2 + 4*p1p2 + p2k1 + p2k2) + 2*(3*p1k1**2 + p1k2**2 - 2*p1k2*p1p2 - 2*M2*p2k1 + 2
     &*p1k2*p2k1 + p1k1*(2*p1k2 + 2*p1p2 - p2k1 - 4*p2k2) + 2*M2*p2k2 - 3*p1k2*p2k2))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k2*(k1k2 - p1k1 - p1k2)*p1k2))

!************************************************   2 Re(Mvk3 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble((GGG*(2 - 3*Nc**2 + Nc**4)*(3*k1k3**2 + 2*k2k3*M2 - 3
     &*k2k3*p1k1 + 6*p1k1**2 - 2*k2k3*p1k3 + 2*p1k1*p1k3 - 2*p1k2*p1k3 - 6*k2k3*p1p2 + 8*p1k1*p1p2 
     &+ 4*p1k2*p1p2 - k2k3*p2k1 - 8*M2*p2k1 + 2*p1k1*p2k1 + 2*p1k2*p2k1 + 12*p1k3*p2k1 - 4*M2*p2k2 
     &+ 2*p1k1*p2k2 + 6*p1k3*p2k2 - k1k3*(2*M2 + 9*p1k1 - 3*p1k2 - 2*p1k3 + 6*p1p2 + p2k1 + p2k2 - 4
     &*p2k3) + 2*k2k3*p2k3 - 6*p1k1*p2k3 + k1k2*(k1k3 - 2*p1k1 - 4*p1p2 + 2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(64*k1k2*p1k3*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(MtP123 MtP132^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-1 + Nc**2)**2*(M2 - p1k1)*(k1k2 + 2*k1k3 + k2k3 + 2
     &*M2 - 4*p1k1 - p1k2 - 2*p1k3 - 2*p1p2 + 2*p2k1 + p2k2)*(2*k1k2 + k1k3 + k2k3 + 2*M2 - 4*p1k1 
     &- 2*p1k2 - p1k3 - 2*p1p2 + 2*p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(32*Nc
     &**2*p1k1**2*(-k1k2 + p1k1 + p1k2)*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(MtP123 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-1 + Nc**4)*(k1k2 + 2*k1k3 + k2k3 + 2*M2 - 2*p1k1 
     &- 2*p1k2 - 2*p1k3)*(k1k2 + k1k3 + 2*M2 - p1k1 - 2*p1k2 - 2*p1k3 - 2*p1p2 + p2k1)*(k1k3 + k2k3 
     &+ 2*M2 - 2*p1k1 - 2*p1k2 - p1k3 - 2*p1p2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(64.
     &*Nc**2*p1k1*(-k1k2 + p1k1 + p1k2)*p1k3*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(MtP123 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-1 + Nc**2)*(k1k2 + 2*M2 - 2*p1k1 - 2*p1k2)*(k1k2 
     &+ k1k3 + 2*M2 - p1k1 - 2*p1k2 - 2*p1k3 - 2*p1p2 + p2k1)*(2*k1k2 + k1k3 + k2k3 + 2*M2 - 2*p1k1 
     &- 4*p1k2 - p1k3 - 2*p1p2 + 2*p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(64*Nc
     &**2*p1k1*p1k2*(-k1k2 + p1k1 + p1k2)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(MtP123 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-1 + Nc**2)**2*(k1k2 + 2*M2 - 2*p1k1 - 2*p1k2)**2*(k1k2 
     &+ M2 - p1k1 - p1k2 - p1p2 + p2k1 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(16.
     &*Nc**2*p1k1*p1k2*(-k1k2 + p1k1 + p1k2)**2))

!************************************************   2 Re(MtP123 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-1 + Nc**2)*(k1k3 + 2*M2 - 2*p1k1 - 2*p1k3)*(k1k2 
     &+ 2*k1k3 + k2k3 + 2*M2 - 4*p1k1 - p1k2 - 2*p1k3 - 2*p1p2 + 2*p2k1 + p2k2)*(k1k3 + k2k3 + 2*M2 
     &- 2*p1k1 - 2*p1k2 - p1k3 - 2*p1p2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(64*Nc
     &**2*p1k1*(-k1k2 + p1k1 + p1k2)*p1k3*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(MtP132 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-1 + Nc**2)*(k1k3 + 2*M2 - 2*p1k1 - 2*p1k3)*(k1k2 
     &+ k1k3 + 2*M2 - p1k1 - 2*p1k2 - 2*p1k3 - 2*p1p2 + p2k1)*(k1k2 + 2*k1k3 + k2k3 + 2*M2 - 2*p1k1 
     &- p1k2 - 4*p1k3 - 2*p1p2 + p2k2 + 2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(64*Nc
     &**2*p1k1*p1k3*(-k1k3 + p1k1 + p1k3)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(MtP132 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-1 + Nc**4)*(2*k1k2 + k1k3 + k2k3 + 2*M2 - 2*p1k1 
     &- 2*p1k2 - 2*p1k3)*(k1k2 + k1k3 + 2*M2 - p1k1 - 2*p1k2 - 2*p1k3 - 2*p1p2 + p2k1)*(k1k2 + k2k3 
     &+ 2*M2 - 2*p1k1 - p1k2 - 2*p1k3 - 2*p1p2 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(64.
     &*Nc**2*p1k1*p1k2*(-k1k3 + p1k1 + p1k3)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(MtP132 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-1 + Nc**2)*(k1k2 + 2*M2 - 2*p1k1 - 2*p1k2)*(k1k2 
     &+ k2k3 + 2*M2 - 2*p1k1 - p1k2 - 2*p1k3 - 2*p1p2 + p2k2)*(2*k1k2 + k1k3 + k2k3 + 2*M2 - 4*p1k1 
     &- 2*p1k2 - p1k3 - 2*p1p2 + 2*p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(64*Nc
     &**2*p1k1*p1k2*(-k1k2 + p1k1 + p1k2)*(-k1k3 + p1k1 + p1k3)))

!************************************************   2 Re(MtP132 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-1 + Nc**2)**2*(k1k3 + 2*M2 - 2*p1k1 - 2*p1k3)**2*(k1k3 
     &+ M2 - p1k1 - p1k3 - p1p2 + p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(16.
     &*Nc**2*p1k1*p1k3*(-k1k3 + p1k1 + p1k3)**2))

!************************************************   2 Re(MtP321 MtP231^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-1 + Nc**2)**2*(k2k3 + 2*M2 - 2*p1k2 - 2*p1k3)**2*(k2k3 
     &+ M2 - p1k2 - p1k3 - p1p2 + p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(16.
     &*Nc**2*p1k2*p1k3*(-k2k3 + p1k2 + p1k3)**2))

!************************************************   2 Re(MtP321 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-1 + Nc**2)*(k2k3 + 2*M2 - 2*p1k2 - 2*p1k3)*(k1k2 
     &+ k1k3 + 2*k2k3 + 2*M2 - p1k1 - 4*p1k2 - 2*p1k3 - 2*p1p2 + p2k1 + 2*p2k2)*(k1k3 + k2k3 + 2*M2 
     &- 2*p1k1 - 2*p1k2 - p1k3 - 2*p1p2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(64*Nc
     &**2*p1k2*(-k1k2 + p1k1 + p1k2)*p1k3*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(MtP321 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-1 + Nc**2)**2*(M2 - p1k3)*(k1k2 + k1k3 + 2*k2k3 + 2
     &*M2 - p1k1 - 2*p1k2 - 4*p1k3 - 2*p1p2 + p2k1 + 2*p2k3)*(k1k2 + 2*k1k3 + k2k3 + 2*M2 - 2*p1k1 
     &- p1k2 - 4*p1k3 - 2*p1p2 + p2k2 + 2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(32*Nc
     &**2*p1k3**2*(-k1k3 + p1k1 + p1k3)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(MtP231 MtP213^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(((-1 + Nc**2)**2*(M2 - p1k2)*(k1k2 + k1k3 + 2*k2k3 + 2
     &*M2 - p1k1 - 4*p1k2 - 2*p1k3 - 2*p1p2 + p2k1 + 2*p2k2)*(2*k1k2 + k1k3 + k2k3 + 2*M2 - 2*p1k1 
     &- 4*p1k2 - p1k3 - 2*p1p2 + 2*p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(32*Nc
     &**2*p1k2**2*(-k1k2 + p1k1 + p1k2)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(MtP231 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-1 + Nc**2)*(k2k3 + 2*M2 - 2*p1k2 - 2*p1k3)*(k1k2 
     &+ k2k3 + 2*M2 - 2*p1k1 - p1k2 - 2*p1k3 - 2*p1p2 + p2k2)*(k1k2 + k1k3 + 2*k2k3 + 2*M2 - p1k1 
     &- 2*p1k2 - 4*p1k3 - 2*p1p2 + p2k1 + 2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(64*Nc
     &**2*p1k2*p1k3*(-k1k3 + p1k1 + p1k3)*(-k2k3 + p1k2 + p1k3)))

!************************************************   2 Re(MtP213 MtP321^\ast)     ************************************************

      stsT2xx_M2Real3g=stsT2xx_M2Real3g+2*dble(-((-1 + Nc**4)*(k1k2 + k1k3 + 2*(k2k3 + M2 - p1k1 - p1k2 
     &- p1k3))*(k1k2 + k2k3 + 2*M2 - 2*p1k1 - p1k2 - 2*p1k3 - 2*p1p2 + p2k2)*(k1k3 + k2k3 + 2*M2 
     &- 2*p1k1 - 2*p1k2 - p1k3 - 2*p1p2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**6)/(64*Nc
     &**2*p1k2*(-k1k2 + p1k1 + p1k2)*p1k3*(-k1k3 + p1k1 + p1k3)))

     
      end
