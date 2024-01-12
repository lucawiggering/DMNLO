      double precision function stsT2xx_M2RealQ()
      implicit none
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"
#include "common/stsT2xx_GenCouplings.h"

      double precision M2 
      
      M2=MSf2(isf1,itt1,igen1)

!************************************************   |MvQ|^2   ************************************************

      stsT2xx_M2RealQ=((Gs**2*(2 - 3*Nc**2 + Nc**4)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(2*k1k2*Nc))

!************************************************   |MtQ|^2   ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+((Gs**2*(-1 + Nc**2)**2*(k1k2*M2 - 2*p1k1*p1k2)*(M2 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(2*k1k2**2*Nc*p2k3**2))

!************************************************   |MbQ|^2   ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+((Gs**2*(-1 + Nc**2)**2*(M2 - p1k3)*(-(k1k3*(k2k3 - p1k2 + p2k2)) 
     &+ (p1k1 - p2k1)*(k2k3 - p1k2 + p2k2) + k1k2*(M2 - p1k3 - p1p2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(4*k1k2**2*Nc*p1k3**2))

!************************************************   |Ms3Q|^2   ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+((GGG**2*Gs**2*Nc*(-1 + Nc**2)*(-2*k2k3*M2*p1k1 - 4*M2*p1k1*p1k2 
     &+ k2k3*p1k1*p1k3 + 2*k1k2**2*(M2 - p1p2) + 2*k2k3*p1k1*p1p2 - 2*k2k3*M2*p2k1 - 2*p1k1*p1k2
     &*p2k1 + 2*p1k2**2*p2k1 - k2k3*p1k3*p2k1 + 2*k2k3*p1p2*p2k1 + 4*p1k2*p1p2*p2k1 + 2*p1k2*p2k1*
     &*2 + 2*p1k1**2*p2k2 - 2*p1k1*p1k2*p2k2 + 4*p1k1*p1p2*p2k2 - 4*M2*p2k1*p2k2 - 2*p1k1*p2k1*p2k2 
     &- 2*p1k2*p2k1*p2k2 + 2*p1k1*p2k2**2 - k2k3*p1k1*p2k3 + k2k3*p2k1*p2k3 + k1k3*(p1k2*p1k3 - 2*k2k3
     &*(M2 - p1p2) + 2*p1k2*p1p2 - p1k3*p2k2 + 2*p1p2*p2k2 - 2*M2*(p1k2 + p2k2) - p1k2*p2k3 + p2k2
     &*p2k3) + k1k2*(4*M2**2 + p1k1**2 + p1k2**2 - 2*p1k2*p1k3 - 2*p1k2*p1p2 - 2*p1k3*p1p2 - 4*p1p2
     &**2 + 2*p1k3*p2k1 - 2*p1p2*p2k1 + p2k1**2 - 2*p1k2*p2k2 + 2*p1k3*p2k2 - 2*p1p2*p2k2 + p2k2
     &**2 - 2*p1k1*(p1k3 + p1p2 + p2k1 - p2k3) + 2*p1k2*p2k3 - 2*p1p2*p2k3 - 2*p2k1*p2k3 - 2*p2k2*p2k3 
     &+ 2*M2*(p1k1 + p1k2 + p1k3 + p2k1 + p2k2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(8*k1k2**2*(M2 + p1p2)**2))

!************************************************   |Ms1Q|^2   ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+((Gs**4*(-1 + Nc**2)**2*(k1k2*(M2**2 + p1k1**2 - p1p2**2 - 2*p1k1
     &*p2k1 + p2k1**2) - 2*(p1k2*p2k1*(-p1p2 + p2k1) + p1k1**2*p2k2 + M2*(p1k1*p1k2 + p2k1*p2k2) 
     &- p1k1*(p1k2*p2k1 + (p1p2 + p2k1)*p2k2)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(4*k2k3
     &**2*Nc*(M2 + p1p2)**2))

!************************************************   |Ms2Q|^2   ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+((Gs**4*(-1 + Nc**2)**2*(k1k2*(M2**2 + p1k2**2 - p1p2**2 - 2*p1k2
     &*p2k2 + p2k2**2) - 2*(p1k2**2*p2k1 + p1k1*p2k2*(-p1p2 + p2k2) + M2*(p1k1*p1k2 + p2k1*p2k2) 
     &- p1k2*(p1p2*p2k1 + (p1k1 + p2k1)*p2k2)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(4*k1k3
     &**2*Nc*(M2 + p1p2)**2))

!************************************************   2 Re(MvQ MtQ^ast)     ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+2*dble((Gs**2*(-2 + Nc**2)*(-1 + Nc**2)*(-(k2k3*p1k1) - k1k3*p1k2 
     &+ k1k2*p1k3 - 2*k1k2*p1p2 + 2*p1k2*p2k1 + 2*p1k1*p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(8*k1k2**2*Nc*p2k3))

!************************************************   2 Re(MvQ MbQ^ast)     ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+2*dble(-(Gs**2*(-2 + Nc**2)*(-1 + Nc**2)*(-3*k2k3*p1k1 + 4*p1k1
     &*p1k2 + k2k3*p2k1 - 2*p1k2*p2k1 - 2*p1k1*p2k2 + k1k3*(2*k2k3 - 3*p1k2 + p2k2) - k1k2*(2*M2 
     &- 3*p1k3 - 2*p1p2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(16
     &*k1k2**2*Nc*p1k3))

!************************************************   2 Re(MvQ Ms3Q^ast)     ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+2*dble(0)

!************************************************   2 Re(MvQ Ms1Q^ast)     ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+2*dble((Gs**3*(4 - 5*Nc**2 + Nc**4)*(k1k2*(p1k1 - p2k1) + p1k2
     &*p2k1 - p1k1*p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(4
     &*k1k2*k2k3*Nc*(M2 + p1p2)))

!************************************************   2 Re(MvQ Ms2Q^ast)     ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+2*dble((Gs**3*(4 - 5*Nc**2 + Nc**4)*(p1k2*p2k1 - p1k1*p2k2 + k1k2
     &*(-p1k2 + p2k2))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(4
     &*k1k2*k1k3*Nc*(M2 + p1p2)))

!************************************************   2 Re(MtQ MbQ^ast)     ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+2*dble(-(Gs**2*(-1 + Nc**2)*(k2k3*p1k1 + k1k3*p1k2 - 2*p1k1*p1k2 
     &+ k1k2*(M2 - p1k3 - p1p2) + p1k2*p2k1 + p1k1*p2k2)*(p1k3 - 2*p1p2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(8*k1k2**2*Nc*p1k3*p2k3))

!************************************************   2 Re(MtQ Ms3Q^ast)     ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+2*dble(-(GGG*Gs**2*Nc*(-1 + Nc**2)*(2*k1k3*M2*p1k2 - k1k3*p1k1
     &*p1k2 + 8*M2*p1k1*p1k2 + k1k3*p1k2**2 - 2*k1k3*p1k2*p1p2 + 2*p1k1*p1k2*p2k1 - 2*p1k2**2*p2k1 
     &+ 4*p1k2*p1k3*p2k1 - 4*p1k2*p1p2*p2k1 + k2k3*(2*M2*p1k1 + p1k1**2 + p1k2*p2k1 - p1k1*(p1k2 + 2
     &*p1p2 + p2k1)) + k1k3*p1k1*p2k2 - 2*p1k1**2*p2k2 - k1k3*p1k2*p2k2 + 2*p1k1*p1k2*p2k2 + 4*p1k1
     &*p1k3*p2k2 - 4*p1k1*p1p2*p2k2 - 4*p1k1*p1k2*p2k3 - 2*p1k2*p2k1*p2k3 - 2*p1k1*p2k2*p2k3 + k1k2
     &*(-4*M2**2 - 2*M2*p1k3 - p1k1*p1k3 - p1k2*p1k3 + k1k3*(M2 - p1p2) + k2k3*(M2 - p1p2) + 2*p1k1
     &*p1p2 + 2*p1k2*p1p2 - 2*p1k3*p1p2 + 4*p1p2**2 - 2*M2*p2k1 + p1k3*p2k1 - 2*M2*p2k2 + p1k3*p2k2 
     &+ 2*M2*p2k3 + 2*p1p2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(16*k1k2**2*(M2 + p1p2)*p2k3))

!************************************************   2 Re(MtQ Ms1Q^ast)     ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+2*dble((Gs**3*(2 - 3*Nc**2 + Nc**4)*(-(k1k3*p1k1*p1k2) - 4*M2
     &*p1k1*p1k2 - k1k3*p1k2*p1p2 + k1k3*p1k2*p2k1 + 2*p1k1*p1k2*p2k1 + 2*p1k2*p1p2*p2k1 - 2*p1k2*p2k1
     &**2 + k2k3*(-p1k1**2 - M2*p2k1 + p1k1*(p1p2 + p2k1)) + k1k3*M2*p2k2 + 2*p1k1**2*p2k2 - 2*p1k1
     &*p1k3*p2k2 + 2*p1k1*p1p2*p2k2 - 2*p1k1*p2k1*p2k2 + 2*p1k1*p1k2*p2k3 + k1k2*(2*M2**2 + (p1k3 
     &- 2*p1p2)*(p1k1 + p1p2 - p2k1) - M2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(8*k1k2
     &*k2k3*Nc*(M2 + p1p2)*p2k3))

!************************************************   2 Re(MtQ Ms2Q^ast)     ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+2*dble((Gs**3*(-1 + Nc**2)*(-4*M2*p1k1*p1k2 - k1k3*p1k2**2 + k1k3
     &*p1k2*p1p2 + k2k3*M2*p2k1 + 2*p1k2**2*p2k1 - 2*p1k2*p1k3*p2k1 + 2*p1k2*p1p2*p2k1 - k2k3*p1k1
     &*(p1k2 + p1p2 - p2k2) - k1k3*M2*p2k2 + k1k3*p1k2*p2k2 + 2*p1k1*p1k2*p2k2 + 2*p1k1*p1p2*p2k2 
     &- 2*p1k2*p2k1*p2k2 - 2*p1k1*p2k2**2 + 2*p1k1*p1k2*p2k3 + k1k2*(2*M2**2 + (p1k3 - 2*p1p2)*(p1k2 
     &+ p1p2 - p2k2) - M2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(4*k1k2*k1k3*Nc*(M2 + p1p2)*p2k3))

!************************************************   2 Re(MbQ Ms3Q^ast)     ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+2*dble((GGG*Gs**2*Nc*(-8*(2*M2 - p1k3 - 2*p1p2 + p2k3)*(p1k1*p1k2 
     &- k2k3*p2k1 - p2k1*p2k2 - k1k3*(k2k3 + p2k2) + k1k2*p2k3) + 8*Nc**2*(2*M2 - p1k3 - 2*p1p2 
     &+ p2k3)*(p1k1*p1k2 - k2k3*p2k1 - p2k1*p2k2 - k1k3*(k2k3 + p2k2) + k1k2*p2k3) - 4*(p1k1 + p1k2 
     &- p1k3 - p2k1 - p2k2 + p2k3)*(-3*k2k3*p1k1 + 4*p1k1*p1k2 + k2k3*p2k1 - 2*p1k2*p2k1 - 2*p1k1
     &*p2k2 + k1k3*(2*k2k3 - 3*p1k2 + p2k2) - k1k2*(2*M2 - 3*p1k3 - 2*p1p2 + p2k3)) + 4*Nc**2*(p1k1 
     &+ p1k2 - p1k3 - p2k1 - p2k2 + p2k3)*(-3*k2k3*p1k1 + 4*p1k1*p1k2 + k2k3*p2k1 - 2*p1k2*p2k1 - 2
     &*p1k1*p2k2 + k1k3*(2*k2k3 - 3*p1k2 + p2k2) - k1k2*(2*M2 - 3*p1k3 - 2*p1p2 + p2k3)) + 4*(k1k3 
     &+ k2k3 - 2*M2 - 2*p1k1 - 2*p1k2 + p1k3 - 2*p1p2 + p2k3)*(k2k3*(p1k1 - p2k1) + (k1k3 - 2*p1k1 
     &+ 2*p2k1)*(p1k2 - p2k2) + k1k2*(2*M2 - p1k3 - 2*p1p2 + p2k3)) - 4*Nc**2*(k1k3 + k2k3 - 2*M2 
     &- 2*p1k1 - 2*p1k2 + p1k3 - 2*p1p2 + p2k3)*(k2k3*(p1k1 - p2k1) + (k1k3 - 2*p1k1 + 2*p2k1)*(p1k2 
     &- p2k2) + k1k2*(2*M2 - p1k3 - 2*p1p2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2**2*p1k3*(M2 + p1p2)))

!************************************************   2 Re(MbQ Ms1Q^ast)     ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+2*dble((Gs**3*(-1 + Nc**2)*(-(k2k3*M2*p1k1) - 3*k2k3*p1k1**2 
     &+ 4*p1k1**2*p1k2 + 3*k2k3*p1k1*p1p2 - 4*p1k1*p1k2*p1p2 - 3*k2k3*M2*p2k1 + 4*k2k3*p1k1*p2k1 + 2
     &*M2*p1k2*p2k1 - 6*p1k1*p1k2*p2k1 + 2*k2k3*p1k3*p2k1 + k2k3*p1p2*p2k1 + 2*p1k2*p1p2*p2k1 - k2k3
     &*p2k1**2 + 2*p1k2*p2k1**2 + 2*M2*p1k1*p2k2 - 2*p1k1**2*p2k2 - 2*p1k1*p1k3*p2k2 + 2*p1k1*p1p2
     &*p2k2 - 4*M2*p2k1*p2k2 + 2*p1k1*p2k1*p2k2 + 2*p1k3*p2k1*p2k2 + k1k3*(2*k2k3*p1k1 + M2*p1k2 - 3
     &*p1k1*p1k2 + p1k2*p1p2 - 2*k2k3*p2k1 + 3*p1k2*p2k1 - M2*p2k2 + p1k1*p2k2 - p1p2*p2k2 - p2k1
     &*p2k2) - 2*k2k3*p1k1*p2k3 + 2*p1k1*p1k2*p2k3 - 2*p1k2*p2k1*p2k3 + k1k2*(2*M2**2 + 3*p1k1*p1k3 
     &+ 2*p1k1*p1p2 - p1k3*p1p2 - 2*p1p2**2 - 3*p1k3*p2k1 - 2*p1p2*p2k1 - p1k1*p2k3 + p1p2*p2k3 + p2k1
     &*p2k3 + M2*(-2*p1k1 - p1k3 + 2*p2k1 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(8*k1k2*k2k3*Nc*p1k3*(M2 + p1p2)))

!************************************************   2 Re(MbQ Ms2Q^ast)     ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+2*dble((Gs**3*(2 - 3*Nc**2 + Nc**4)*(-(k1k3*M2*p1k2) - 3*k1k3
     &*p1k2**2 + 4*p1k1*p1k2**2 + 3*k1k3*p1k2*p1p2 - 4*p1k1*p1k2*p1p2 + 2*M2*p1k2*p2k1 - 2*p1k2**2*p2k1 
     &- 2*p1k2*p1k3*p2k1 + 2*p1k2*p1p2*p2k1 - 3*k1k3*M2*p2k2 + 2*M2*p1k1*p2k2 + 4*k1k3*p1k2*p2k2 
     &- 6*p1k1*p1k2*p2k2 + 2*k1k3*p1k3*p2k2 + k1k3*p1p2*p2k2 + 2*p1k1*p1p2*p2k2 - 4*M2*p2k1*p2k2 
     &+ 2*p1k2*p2k1*p2k2 + 2*p1k3*p2k1*p2k2 - k1k3*p2k2**2 + 2*p1k1*p2k2**2 + k2k3*(M2*p1k1 + 2*k1k3
     &*p1k2 - 3*p1k1*p1k2 + p1k1*p1p2 - M2*p2k1 + p1k2*p2k1 - p1p2*p2k1 - 2*k1k3*p2k2 + 3*p1k1*p2k2 
     &- p2k1*p2k2) - 2*k1k3*p1k2*p2k3 + 2*p1k1*p1k2*p2k3 - 2*p1k1*p2k2*p2k3 + k1k2*(2*M2**2 + 3*p1k2
     &*p1k3 + 2*p1k2*p1p2 - p1k3*p1p2 - 2*p1p2**2 - 3*p1k3*p2k2 - 2*p1p2*p2k2 - p1k2*p2k3 + p1p2
     &*p2k3 + p2k2*p2k3 + M2*(-2*p1k2 - p1k3 + 2*p2k2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(16*k1k2*k1k3*Nc*p1k3*(M2 + p1p2)))

!************************************************   2 Re(Ms3Q Ms1Q^ast)     ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+2*dble((GGG*Gs**3*Nc*(-1 + Nc**2)*(-(k1k3*M2*p1k2) - k1k3*p1k1
     &*p1k2 - 4*M2*p1k1*p1k2 + k1k3*p1k2*p1p2 + k1k3*p1k2*p2k1 + 2*p1k2**2*p2k1 - 2*p1k2*p1k3*p2k1 
     &+ 4*p1k2*p1p2*p2k1 - k2k3*(p1k1**2 + M2*(p1k1 + p2k1) + p2k1*(-p1p2 + p2k1) - p1k1*(p1p2 + 2
     &*p2k1)) - k1k3*M2*p2k2 + k1k3*p1k1*p2k2 - 2*p1k1*p1k2*p2k2 + 2*p1k1*p1k3*p2k2 + k1k3*p1p2*p2k2 
     &+ 4*p1k1*p1p2*p2k2 - k1k3*p2k1*p2k2 - 4*M2*p2k1*p2k2 - 2*p1k2*p2k1*p2k2 + 2*p1k1*p2k2**2 + 2
     &*p1k2*p2k1*p2k3 - 2*p1k1*p2k2*p2k3 + k1k2*(4*M2**2 + 2*p1k1**2 - p1k3*p1p2 - 4*p1p2**2 - 2*p1k2
     &*p2k1 + p1k3*p2k1 - 2*p1p2*p2k1 + 2*p2k1**2 + 2*p2k1*p2k2 - p1p2*p2k3 - p2k1*p2k3 + M2*(2*p1k1 
     &+ p1k3 + 2*p2k1 + p2k3) + p1k1*(2*p1k2 - p1k3 - 2*p1p2 - 4*p2k1 - 2*p2k2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(16*k1k2*k2k3*(M2 + p1p2)**2))

!************************************************   2 Re(Ms3Q Ms2Q^ast)     ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+2*dble((GGG*Gs**3*Nc*(-1 + Nc**2)*(-(k1k3*M2*p1k2) - 4*M2*p1k1
     &*p1k2 - k1k3*p1k2**2 + k1k3*p1k2*p1p2 - 2*p1k1*p1k2*p2k1 + 2*p1k2*p1k3*p2k1 + 4*p1k2*p1p2*p2k1 
     &+ 2*p1k2*p2k1**2 - k1k3*M2*p2k2 + 2*p1k1**2*p2k2 + 2*k1k3*p1k2*p2k2 - 2*p1k1*p1k3*p2k2 + k1k3
     &*p1p2*p2k2 + 4*p1k1*p1p2*p2k2 - 4*M2*p2k1*p2k2 - 2*p1k1*p2k1*p2k2 - k1k3*p2k2**2 + k2k3*(-(M2
     &*(p1k1 + p2k1)) + p2k1*(p1k2 + p1p2 - p2k2) + p1k1*(-p1k2 + p1p2 + p2k2)) - 2*p1k2*p2k1*p2k3 
     &+ 2*p1k1*p2k2*p2k3 + k1k2*(4*M2**2 + 2*p1k2**2 - p1k2*p1k3 - 2*p1k2*p1p2 - p1k3*p1p2 - 4*p1p2
     &**2 - 2*p1k2*p2k1 + 2*p1k1*(p1k2 - p2k2) - 4*p1k2*p2k2 + p1k3*p2k2 - 2*p1p2*p2k2 + 2*p2k1*p2k2 
     &+ 2*p2k2**2 + p1k2*p2k3 - p1p2*p2k3 - p2k2*p2k3 + M2*(2*p1k2 + p1k3 + 2*p2k2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(16*k1k2*k1k3*(M2 + p1p2)**2))

!************************************************   2 Re(Ms1Q Ms2Q^ast)     ************************************************

      stsT2xx_M2RealQ=stsT2xx_M2RealQ+2*dble((Gs**4*(-1 + Nc**2)*(p1k2*p2k1*(p1k2 + p2k1 - p2k2) + p1k1
     &**2*p2k2 + k1k2*(M2**2 + p1k1*p1k2 - p1p2**2 - p1k2*p2k1 - p1k1*p2k2 + p2k1*p2k2) - p1k1*((p2k1 
     &- p2k2)*p2k2 + p1k2*(p2k1 + p2k2)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(4*k1k3
     &*k2k3*Nc*(M2 + p1p2)**2))
     
      end
