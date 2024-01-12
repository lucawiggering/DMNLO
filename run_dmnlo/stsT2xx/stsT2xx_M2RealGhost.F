      double precision function stsT2xx_M2RealGhost()
      implicit none
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"
#include "common/stsT2xx_GenCouplings.h"
 
      double precision M2 

      M2=MSf2(isf1,itt1,igen1)

!##############################################################################################################################
!#################################################    |S1|^2   ################################################################   
!##############################################################################################################################

!************************************************   |MvC123|^2   ************************************************

      stsT2xx_M2RealGhost=0d0

!************************************************   2 Re(MvC123 MtC123^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**2*(2 - 3*Nc**2 + Nc**4)*(k1k2 - 2*p1k1)*(k1k2 
     &- p1k1 + p2k1)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32
     &*k1k2*(k1k2 - p1k1 - p1k2)*(M2 - p1k3 + p1p2 - p2k3)))

!************************************************   2 Re(MvC123 MbC123^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**2*(2 - 3*Nc**2 + Nc**4)*(k1k3 - 2*p1k1)*(k1k3 
     &- p1k1 + p2k1)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32
     &*k1k2*p1k3*(-M2 + p1k3 - p1p2 + p2k3)))

!************************************************   2 Re(MvC123 Ms1C123^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   2 Re(MvC123 Ms2C123^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   2 Re(MvC123 Ms3C123^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   |MtC123|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**2*(-1 + Nc**2)**2*(k1k2 - 2*p1k1)**2*(k1k2 + M2 - p1k1 
     &- p1k2 - p1p2 + p2k1 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(32*k1k2**2*(-k1k2 + p1k1 + p1k2)**2))

!************************************************   2 Re(MtC123 MbC123^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**2*(-1 + Nc**2)*(k1k2 - 2*p1k1)*(-k1k3 + p1k1 
     &- p2k1)*(-k1k3 - k2k3 - 2*M2 + 2*p1k1 + 2*p1k2 + p1k3 + 2*p1p2 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(64*k1k2**2*(k1k2 - p1k1 - p1k2)*p1k3))

!************************************************   2 Re(MtC123 Ms1C123^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k1k2 - 2*p1k1)*(p1k1 - p2k1)
     &*(k1k2 + k1k3 + k2k3 - p1k2 - p1k3 + p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3)/(128*k1k2*k2k3*(k1k2 - p1k1 - p1k2)*(M2 + p1p2)))

!************************************************   2 Re(MtC123 Ms2C123^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k1k2 - 2*p1k1)*(k1k2 - p1k1 
     &+ p2k1)*(p1k1 + p1k3 - p2k1 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2
     &*k1k3*(k1k2 - p1k1 - p1k2)*(M2 + p1p2)))

!************************************************   2 Re(MtC123 Ms3C123^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GGG*GhGhG**2*Nc**2*(-1 + Nc**2)*(-k1k2 + 2*p1k1)*(2*(p1k1 
     &- p2k1)*(k1k2 + p2k1 + p2k2) + (k1k3 + p1k1 + p2k1)*(2*M2 - p1k1 - p1k2 - 2*p1p2 + p2k1 
     &+ p2k2) - (k1k2 - p1k1 + p2k1)*(p1k1 + p1k2 - p1k3 - p2k1 - p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3)/(64*k1k2**2*(k1k2 - p1k1 - p1k2)*(M2 + p1p2)))

!************************************************   |MbC123|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**2*(-1 + Nc**2)**2*(M2 - p1k3)*(k1k3 - p1k1 + p2k1)*
     &*2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(16*k1k2**2*p1k3**2))

!************************************************   2 Re(MbC123 Ms1C123^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**3*Nc**2*(-1 + Nc**2)*(k2k3 - 2*(p1k2 + p1k3))*(p1k1 
     &- p2k1)*(-k1k3 + p1k1 - p2k1)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2*k2k3*p1k3*(M2 + p1p2)))

!************************************************   2 Re(MbC123 Ms2C123^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k1k3 - 2*p1k1)*(k1k3 - p1k1 
     &+ p2k1)*(p1k1 + p1k3 - p2k1 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2*k1k3*p1k3*(M2 + p1p2)))

!************************************************   2 Re(MbC123 Ms3C123^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GGG*GhGhG**2*Nc**2*(-1 + Nc**2)*(k1k3 - p1k1 + p2k1)*(
     &-2*(2*M2 + p1k2 - p1k3)*p2k1 + k2k3*(-p1k1 + p2k1) + k1k3*(-2*M2 + p1k2 + 2*p1p2 - p2k2) + 2
     &*p1k1*(p1k3 + 2*p1p2 + p2k2 - 2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(64*k1k2**2*p1k3*(M2 + p1p2)))

!************************************************   |Ms1C123|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**4*Nc**2*(-1 + Nc**2)*(p1k1 - p2k1)**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(16*k2k3*(M2 + p1p2)**2))

!************************************************   2 Re(Ms1C123 Ms2C123^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**4*(k1k2 + k1k3)*Nc**2*(-1 + Nc**2)*(p1k1 - p2k1)
     &*(p1k1 + p1k3 - p2k1 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k3*k2k3*(M2 + p1p2)**2))

!************************************************   2 Re(Ms1C123 Ms3C123^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GGG*GhGhG**3*Nc**2*(-1 + Nc**2)*(p1k1 - p2k1)*(k2k3*(
     &-p1k1 + p2k1) + k1k2*(p1k2 - p1k3 - p2k2 + p2k3) + 2*(k1k3*p1k2 + p1k2*p2k1 + p1k3*p2k1 - k1k3
     &*p2k2 - p1k1*p2k2 - p1k1*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k2*k2k3*(M2 + p1p2)**2))

!************************************************   |Ms2C123|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(0)

!************************************************   2 Re(Ms2C123 Ms3C123^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GGG*GhGhG**3*(k1k2 - k1k3)*Nc**2*(-1 + Nc**2)*(p1k1 
     &- p2k1)*(p1k1 + p1k3 - p2k1 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k2*k1k3*(M2 + p1p2)**2))

!************************************************   |Ms3C123|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+((GGG**2*GhGhG**2*Nc**2*(-1 + Nc**2)*(k1k2*p1k1*p1k2 - k1k2*p1k1
     &*p1k3 + k1k3**2*(-M2 + p1p2) - k1k2*p1k2*p2k1 + 2*p1k1*p1k2*p2k1 + k1k2*p1k3*p2k1 + 4*p1k1
     &*p1p2*p2k1 - 2*p1k2*p2k1**2 - 2*M2*(p1k1**2 + p2k1**2) - k1k2*p1k1*p2k2 - 2*p1k1**2*p2k2 + k1k2
     &*p2k1*p2k2 + 2*p1k1*p2k1*p2k2 + k1k2*p1k1*p2k3 - k1k2*p2k1*p2k3 + k1k3*(-2*M2*(p1k1 + p2k1) 
     &+ p1k1*(p1k3 + 2*p1p2 - p2k3) + p2k1*(-p1k3 + 2*p1p2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(16*k1k2**2*(M2 + p1p2)**2))

!##############################################################################################################################
!#################################################    |S2|^2   ################################################################   
!##############################################################################################################################

!************************************************   |MvC213|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(0)

!************************************************   2 Re(MvC213 MtC213^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**2*(2 - 3*Nc**2 + Nc**4)*(k1k2 - 2*p1k2)*(k1k2 
     &- p1k2 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32
     &*k1k2*(k1k2 - p1k1 - p1k2)*(M2 - p1k3 + p1p2 - p2k3)))

!************************************************   2 Re(MvC213 MbC213^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**2*(2 - 3*Nc**2 + Nc**4)*(k2k3 - 2*p1k2)*(k2k3 
     &- p1k2 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32
     &*k1k2*p1k3*(-M2 + p1k3 - p1p2 + p2k3)))

!************************************************   2 Re(MvC213 Ms1C213^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   2 Re(MvC213 Ms2C213^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   2 Re(MvC213 Ms3C213^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   |MtC213|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**2*(-1 + Nc**2)**2*(k1k2 - 2*p1k2)**2*(k1k2 + M2 - p1k1 
     &- p1k2 - p1p2 + p2k1 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(32*k1k2**2*(-k1k2 + p1k1 + p1k2)**2))

!************************************************   2 Re(MtC213 MbC213^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**2*(-1 + Nc**2)*(k1k2 - 2*p1k2)*(-k2k3 + p1k2 
     &- p2k2)*(-k1k3 - k2k3 - 2*M2 + 2*p1k1 + 2*p1k2 + p1k3 + 2*p1p2 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(64*k1k2**2*(k1k2 - p1k1 - p1k2)*p1k3))

!************************************************   2 Re(MtC213 Ms1C213^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k1k2 - 2*p1k2)*(k1k2 - p1k2 
     &+ p2k2)*(p1k2 + p1k3 - p2k2 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2
     &*k2k3*(k1k2 - p1k1 - p1k2)*(M2 + p1p2)))

!************************************************   2 Re(MtC213 Ms2C213^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k1k2 - 2*p1k2)*(p1k2 - p2k2)
     &*(k1k2 + k1k3 + k2k3 - p1k1 - p1k3 + p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3)/(128*k1k2*k1k3*(k1k2 - p1k1 - p1k2)*(M2 + p1p2)))

!************************************************   2 Re(MtC213 Ms3C213^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GGG*GhGhG**2*Nc**2*(-1 + Nc**2)*(-k1k2 + 2*p1k2)*(2*(p1k2 
     &- p2k2)*(k1k2 + p2k1 + p2k2) + (k2k3 + p1k2 + p2k2)*(2*M2 - p1k1 - p1k2 - 2*p1p2 + p2k1 
     &+ p2k2) - (k1k2 - p1k2 + p2k2)*(p1k1 + p1k2 - p1k3 - p2k1 - p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3)/(64*k1k2**2*(k1k2 - p1k1 - p1k2)*(M2 + p1p2)))

!************************************************   |MbC213|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**2*(-1 + Nc**2)**2*(M2 - p1k3)*(k2k3 - p1k2 + p2k2)*
     &*2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(16*k1k2**2*p1k3**2))

!************************************************   2 Re(MbC213 Ms1C213^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k2k3 - 2*p1k2)*(k2k3 - p1k2 
     &+ p2k2)*(p1k2 + p1k3 - p2k2 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2*k2k3*p1k3*(M2 + p1p2)))

!************************************************   2 Re(MbC213 Ms2C213^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**3*Nc**2*(-1 + Nc**2)*(k1k3 - 2*(p1k1 + p1k3))*(p1k2 
     &- p2k2)*(-k2k3 + p1k2 - p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2*k1k3*p1k3*(M2 + p1p2)))

!************************************************   2 Re(MbC213 Ms3C213^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GGG*GhGhG**2*Nc**2*(-1 + Nc**2)*(k2k3 - p1k2 + p2k2)*(k2k3
     &*(-2*M2 + p1k1 + 2*p1p2 - p2k1) - 2*(2*M2 + p1k1 - p1k3)*p2k2 + k1k3*(-p1k2 + p2k2) + 2*p1k2
     &*(p1k3 + 2*p1p2 + p2k1 - 2*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(64*k1k2**2*p1k3*(M2 + p1p2)))

!************************************************   |Ms1C213|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(0)

!************************************************   2 Re(Ms1C213 Ms2C213^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**4*(k1k2 + k2k3)*Nc**2*(-1 + Nc**2)*(p1k2 - p2k2)
     &*(p1k2 + p1k3 - p2k2 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k3*k2k3*(M2 + p1p2)**2))

!************************************************   2 Re(Ms1C213 Ms3C213^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GGG*GhGhG**3*(k1k2 - k2k3)*Nc**2*(-1 + Nc**2)*(p1k2 
     &- p2k2)*(p1k2 + p1k3 - p2k2 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k2*k2k3*(M2 + p1p2)**2))

!************************************************   |Ms2C213|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**4*Nc**2*(-1 + Nc**2)*(p1k2 - p2k2)**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(16*k1k3*(M2 + p1p2)**2))

!************************************************   2 Re(Ms2C213 Ms3C213^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GGG*GhGhG**3*Nc**2*(-1 + Nc**2)*(p1k2 - p2k2)*(-(k1k3
     &*p1k2) + 2*k2k3*(p1k1 - p2k1) - 2*p1k2*p2k1 + k1k3*p2k2 + 2*p1k1*p2k2 + 2*p1k3*p2k2 - 2*p1k2*p2k3 
     &+ k1k2*(p1k1 - p1k3 - p2k1 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k2*k1k3*(M2 + p1p2)**2))

!************************************************   |Ms3C213|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GGG**2*GhGhG**2*Nc**2*(-1 + Nc**2)*(k2k3**2*(M2 - p1p2) + 2
     &*(p1k2**2*p2k1 - p1k2*(p1k1 + 2*p1p2 + p2k1)*p2k2 + p1k1*p2k2**2 + M2*(p1k2**2 + p2k2**2)) - k1k2
     &*(p1k2 - p2k2)*(p1k1 - p1k3 - p2k1 + p2k3) + k2k3*(2*M2*(p1k2 + p2k2) + p2k2*(p1k3 - 2*p1p2 
     &- p2k3) + p1k2*(-p1k3 - 2*p1p2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(16*k1k2**2*(M2 + p1p2)**2))


!##############################################################################################################################
!#################################################    |S3|^2   ################################################################   
!##############################################################################################################################

!************************************************   |MvC132|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(0)

!************************************************   2 Re(MvC132 MtC132^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**2*(2 - 3*Nc**2 + Nc**4)*(k1k3 - 2*p1k1)*(k1k3 
     &- p1k1 + p2k1)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32
     &*k1k3*(k1k3 - p1k1 - p1k3)*(M2 - p1k2 + p1p2 - p2k2)))

!************************************************   2 Re(MvC132 MbC132^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**2*(2 - 3*Nc**2 + Nc**4)*(k1k2 - 2*p1k1)*(k1k2 
     &- p1k1 + p2k1)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32
     &*k1k3*p1k2*(-M2 + p1k2 - p1p2 + p2k2)))

!************************************************   2 Re(MvC132 Ms1C132^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   2 Re(MvC132 Ms2C132^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   2 Re(MvC132 Ms3C132^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   |MtC132|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**2*(-1 + Nc**2)**2*(k1k3 - 2*p1k1)**2*(k1k3 + M2 - p1k1 
     &- p1k3 - p1p2 + p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(32*k1k3**2*(-k1k3 + p1k1 + p1k3)**2))

!************************************************   2 Re(MtC132 MbC132^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**2*(-1 + Nc**2)*(k1k3 - 2*p1k1)*(-k1k2 + p1k1 
     &- p2k1)*(-k1k2 - k2k3 - 2*M2 + 2*p1k1 + p1k2 + 2*p1k3 + 2*p1p2 - p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(64*k1k3**2*p1k2*(k1k3 - p1k1 - p1k3)))

!************************************************   2 Re(MtC132 Ms1C132^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k1k3 - 2*p1k1)*(p1k1 - p2k1)
     &*(k1k2 + k1k3 + k2k3 - p1k2 - p1k3 + p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3)/(128*k1k3*k2k3*(k1k3 - p1k1 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(MtC132 Ms2C132^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k1k3 - 2*p1k1)*(k1k3 - p1k1 
     &+ p2k1)*(p1k1 + p1k2 - p2k1 - p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2
     &*k1k3*(k1k3 - p1k1 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(MtC132 Ms3C132^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GGG*GhGhG**2*Nc**2*(-1 + Nc**2)*(-k1k3 + 2*p1k1)*(-((k1k3 
     &- p1k1 + p2k1)*(p1k1 - p1k2 + p1k3 - p2k1 + p2k2 - p2k3)) + 2*(p1k1 - p2k1)*(k1k3 + p2k1 
     &+ p2k3) + (k1k2 + p1k1 + p2k1)*(2*M2 - p1k1 - p1k3 - 2*p1p2 + p2k1 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3)/(64*k1k3**2*(k1k3 - p1k1 - p1k3)*(M2 + p1p2)))

!************************************************   |MbC132|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**2*(-1 + Nc**2)**2*(M2 - p1k2)*(k1k2 - p1k1 + p2k1)*
     &*2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(16*k1k3**2*p1k2**2))

!************************************************   2 Re(MbC132 Ms1C132^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**3*Nc**2*(-1 + Nc**2)*(k2k3 - 2*(p1k2 + p1k3))*(p1k1 
     &- p2k1)*(-k1k2 + p1k1 - p2k1)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k3*k2k3*p1k2*(M2 + p1p2)))

!************************************************   2 Re(MbC132 Ms2C132^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k1k2 - 2*p1k1)*(k1k2 - p1k1 
     &+ p2k1)*(p1k1 + p1k2 - p2k1 - p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2*k1k3*p1k2*(M2 + p1p2)))

!************************************************   2 Re(MbC132 Ms3C132^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GGG*GhGhG**2*Nc**2*(-1 + Nc**2)*(k1k2 - p1k1 + p2k1)*(
     &-2*(2*M2 - p1k2 + p1k3)*p2k1 + k2k3*(-p1k1 + p2k1) + k1k2*(-2*M2 + p1k3 + 2*p1p2 - p2k3) + 2
     &*p1k1*(p1k2 + 2*p1p2 - 2*p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(64*k1k3**2*p1k2*(M2 + p1p2)))

!************************************************   |Ms1C132|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**4*Nc**2*(-1 + Nc**2)*(p1k1 - p2k1)**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(16*k2k3*(M2 + p1p2)**2))

!************************************************   2 Re(Ms1C132 Ms2C132^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**4*(k1k2 + k1k3)*Nc**2*(-1 + Nc**2)*(p1k1 - p2k1)
     &*(p1k1 + p1k2 - p2k1 - p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k2*k2k3*(M2 + p1p2)**2))

!************************************************   2 Re(Ms1C132 Ms3C132^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GGG*GhGhG**3*Nc**2*(-1 + Nc**2)*(p1k1 - p2k1)*(k2k3*(p1k1 
     &- p2k1) + k1k3*(p1k2 - p1k3 - p2k2 + p2k3) + 2*(-(k1k2*p1k3) - p1k2*p2k1 - p1k3*p2k1 + p1k1
     &*p2k2 + k1k2*p2k3 + p1k1*p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k3*k2k3*(M2 + p1p2)**2))

!************************************************   |Ms2C132|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(0)

!************************************************   2 Re(Ms2C132 Ms3C132^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GGG*GhGhG**3*(k1k2 - k1k3)*Nc**2*(-1 + Nc**2)*(p1k1 - p2k1)
     &*(p1k1 + p1k2 - p2k1 - p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k2*k1k3*(M2 + p1p2)**2))

!************************************************   |Ms3C132|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GGG**2*GhGhG**2*Nc**2*(-1 + Nc**2)*(k1k3*p1k1*p1k2 - k1k3*p1k1
     &*p1k3 + k1k2**2*(M2 - p1p2) - k1k3*p1k2*p2k1 + k1k3*p1k3*p2k1 - 2*p1k1*p1k3*p2k1 - 4*p1k1
     &*p1p2*p2k1 + 2*p1k3*p2k1**2 + 2*M2*(p1k1**2 + p2k1**2) - k1k3*p1k1*p2k2 + k1k3*p2k1*p2k2 + k1k2
     &*(2*M2*(p1k1 + p2k1) + p2k1*(p1k2 - 2*p1p2 - p2k2) + p1k1*(-p1k2 - 2*p1p2 + p2k2)) + k1k3*p1k1
     &*p2k3 + 2*p1k1**2*p2k3 - k1k3*p2k1*p2k3 - 2*p1k1*p2k1*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(16*k1k3**2*(M2 + p1p2)**2))


!##############################################################################################################################
!#################################################    |S4|^2   ################################################################   
!##############################################################################################################################

!************************************************   |MvC312|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(0)

!************************************************   2 Re(MvC312 MtC312^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**2*(2 - 3*Nc**2 + Nc**4)*(k1k3 - 2*p1k3)*(k1k3 
     &- p1k3 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32
     &*k1k3*(k1k3 - p1k1 - p1k3)*(M2 - p1k2 + p1p2 - p2k2)))

!************************************************   2 Re(MvC312 MbC312^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**2*(2 - 3*Nc**2 + Nc**4)*(k2k3 - 2*p1k3)*(k2k3 
     &- p1k3 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32
     &*k1k3*p1k2*(-M2 + p1k2 - p1p2 + p2k2)))

!************************************************   2 Re(MvC312 Ms1C312^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   2 Re(MvC312 Ms2C312^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   2 Re(MvC312 Ms3C312^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   |MtC312|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**2*(-1 + Nc**2)**2*(k1k3 - 2*p1k3)**2*(k1k3 + M2 - p1k1 
     &- p1k3 - p1p2 + p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(32*k1k3**2*(-k1k3 + p1k1 + p1k3)**2))

!************************************************   2 Re(MtC312 MbC312^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**2*(-1 + Nc**2)*(k1k3 - 2*p1k3)*(-k1k2 - k2k3 
     &- 2*M2 + 2*p1k1 + p1k2 + 2*p1k3 + 2*p1p2 - p2k2)*(-k2k3 + p1k3 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(64*k1k3**2*p1k2*(k1k3 - p1k1 - p1k3)))

!************************************************   2 Re(MtC312 Ms1C312^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k1k3 - 2*p1k3)*(p1k2 + p1k3 
     &- p2k2 - p2k3)*(k1k3 - p1k3 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k3
     &*k2k3*(k1k3 - p1k1 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(MtC312 Ms2C312^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k1k3 - 2*p1k3)*(k1k2 + k1k3 
     &+ k2k3 - p1k1 - p1k2 + p2k1 + p2k2)*(p1k3 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3)/(128*k1k2*k1k3*(k1k3 - p1k1 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(MtC312 Ms3C312^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GGG*GhGhG**2*Nc**2*(-1 + Nc**2)*(-k1k3 + 2*p1k3)*(-((p1k1 
     &- p1k2 + p1k3 - p2k1 + p2k2 - p2k3)*(k1k3 - p1k3 + p2k3)) + 2*(p1k3 - p2k3)*(k1k3 + p2k1 
     &+ p2k3) + (k2k3 + p1k3 + p2k3)*(2*M2 - p1k1 - p1k3 - 2*p1p2 + p2k1 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3)/(64*k1k3**2*(k1k3 - p1k1 - p1k3)*(M2 + p1p2)))

!************************************************   |MbC312|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**2*(-1 + Nc**2)**2*(M2 - p1k2)*(k2k3 - p1k3 + p2k3)*
     &*2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(16*k1k3**2*p1k2**2))

!************************************************   2 Re(MbC312 Ms1C312^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k2k3 - 2*p1k3)*(p1k2 + p1k3 
     &- p2k2 - p2k3)*(k2k3 - p1k3 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k3*k2k3*p1k2*(M2 + p1p2)))

!************************************************   2 Re(MbC312 Ms2C312^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**3*Nc**2*(-1 + Nc**2)*(k1k2 - 2*(p1k1 + p1k2))*(p1k3 
     &- p2k3)*(-k2k3 + p1k3 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2*k1k3*p1k2*(M2 + p1p2)))

!************************************************   2 Re(MbC312 Ms3C312^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GGG*GhGhG**2*Nc**2*(-1 + Nc**2)*(k2k3 - p1k3 + p2k3)*(k2k3
     &*(-2*M2 + p1k1 + 2*p1p2 - p2k1) + 2*p1k3*(2*p1p2 + p2k1 - 2*p2k2) - 2*(2*M2 + p1k1)*p2k3 
     &+ k1k2*(-p1k3 + p2k3) + 2*p1k2*(p1k3 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(64*k1k3**2*p1k2*(M2 + p1p2)))

!************************************************   |Ms1C312|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(0)

!************************************************   2 Re(Ms1C312 Ms2C312^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**4*(k1k3 + k2k3)*Nc**2*(-1 + Nc**2)*(p1k3 - p2k3)
     &*(p1k2 + p1k3 - p2k2 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k2*k2k3*(M2 + p1p2)**2))

!************************************************   2 Re(Ms1C312 Ms3C312^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GGG*GhGhG**3*(k1k3 - k2k3)*Nc**2*(-1 + Nc**2)*(p1k3 
     &- p2k3)*(p1k2 + p1k3 - p2k2 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k3*k2k3*(M2 + p1p2)**2))

!************************************************   |Ms2C312|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**4*Nc**2*(-1 + Nc**2)*(p1k3 - p2k3)**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(16*k1k2*(M2 + p1p2)**2))

!************************************************   2 Re(Ms2C312 Ms3C312^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GGG*GhGhG**3*Nc**2*(-1 + Nc**2)*(p1k3 - p2k3)*(-2*k2k3
     &*p1k1 + k1k2*p1k3 + 2*k2k3*p2k1 + 2*p1k3*p2k1 + k1k3*(-p1k1 + p1k2 + p2k1 - p2k2) + 2*p1k3*p2k2 
     &- k1k2*p2k3 - 2*p1k1*p2k3 - 2*p1k2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k2*k1k3*(M2 + p1p2)**2))

!************************************************   |Ms3C312|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GGG**2*GhGhG**2*Nc**2*(-1 + Nc**2)*(k2k3**2*(M2 - p1p2) - k1k3
     &*(p1k1 - p1k2 - p2k1 + p2k2)*(p1k3 - p2k3) + k2k3*(-2*p1k3*p1p2 + p1k3*p2k2 - 2*p1p2*p2k3 
     &- p2k2*p2k3 + p1k2*(-p1k3 + p2k3) + 2*M2*(p1k3 + p2k3)) + 2*(p1k3**2*p2k1 - p1k3*(p1k1 + 2*p1p2 
     &+ p2k1)*p2k3 + p1k1*p2k3**2 + M2*(p1k3**2 + p2k3**2)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(16*k1k3**2*(M2 + p1p2)**2))


!##############################################################################################################################
!#################################################    |S5|^2   ################################################################   
!##############################################################################################################################

!************************************************   |MvC321|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(0)

!************************************************   2 Re(MvC321 MtC321^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**2*(2 - 3*Nc**2 + Nc**4)*(k2k3 - 2*p1k3)*(k2k3 
     &- p1k3 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32
     &*k2k3**2*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(MvC321 MbC321^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**2*(2 - 3*Nc**2 + Nc**4)*(k1k3 - 2*p1k3)*(k1k3 
     &- p1k3 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k2k3**2*p1k1))

!************************************************   2 Re(MvC321 Ms1C321^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   2 Re(MvC321 Ms2C321^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   2 Re(MvC321 Ms3C321^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   |MtC321|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**2*(-1 + Nc**2)**2*(k2k3 - 2*p1k3)**2*(k2k3 + M2 - p1k2 
     &- p1k3 - p1p2 + p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(32*k2k3**2*(-k2k3 + p1k2 + p1k3)**2))

!************************************************   2 Re(MtC321 MbC321^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**2*(-1 + Nc**2)*(k2k3 - 2*p1k3)*(-k1k2 - k1k3 
     &- 2*M2 + p1k1 + 2*p1k2 + 2*p1k3 + 2*p1p2 - p2k1)*(-k1k3 + p1k3 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(64*k2k3**2*p1k1*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(MtC321 Ms1C321^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k2k3 - 2*p1k3)*(p1k1 + p1k3 
     &- p2k1 - p2k3)*(k2k3 - p1k3 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k3
     &*k2k3*(k2k3 - p1k2 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(MtC321 Ms2C321^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k2k3 - 2*p1k3)*(k1k2 + k1k3 
     &+ k2k3 - p1k1 - p1k2 + p2k1 + p2k2)*(p1k3 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3)/(128*k1k2*k2k3*(k2k3 - p1k2 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(MtC321 Ms3C321^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GGG*GhGhG**2*Nc**2*(-1 + Nc**2)*(k2k3 - 2*p1k3)*(-2*M2
     &*p1k3 + p1k1*p1k3 + 2*p1k3*p1p2 - p1k3*p2k1 - 2*p1k3*p2k2 + k1k3*(-2*M2 + p1k2 + p1k3 + 2*p1p2 
     &- p2k2 - p2k3) - 2*M2*p2k3 - p1k1*p2k3 + 2*p1k2*p2k3 + 2*p1p2*p2k3 + p2k1*p2k3 + k2k3*(-p1k1 
     &+ p1k2 - p1k3 + p2k1 - p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(64*k2k3**2*(k2k3 - p1k2 - p1k3)*(M2 + p1p2)))

!************************************************   |MbC321|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**2*(-1 + Nc**2)**2*(M2 - p1k1)*(k1k3 - p1k3 + p2k3)*
     &*2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(16*k2k3**2*p1k1**2))

!************************************************   2 Re(MbC321 Ms1C321^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k1k3 - 2*p1k3)*(p1k1 + p1k3 
     &- p2k1 - p2k3)*(k1k3 - p1k3 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k3*k2k3*p1k1*(M2 + p1p2)))

!************************************************   2 Re(MbC321 Ms2C321^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**3*Nc**2*(-1 + Nc**2)*(k1k2 - 2*(p1k1 + p1k2))*(p1k3 
     &- p2k3)*(-k1k3 + p1k3 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2*k2k3*p1k1*(M2 + p1p2)))

!************************************************   2 Re(MbC321 Ms3C321^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GGG*GhGhG**2*Nc**2*(-1 + Nc**2)*(k1k3 - p1k3 + p2k3)*(k1k3
     &*(-2*M2 + p1k2 + 2*p1p2 - p2k2) + 2*p1k3*(2*p1p2 - 2*p2k1 + p2k2) - 2*(2*M2 + p1k2)*p2k3 
     &+ k1k2*(-p1k3 + p2k3) + 2*p1k1*(p1k3 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(64*k2k3**2*p1k1*(M2 + p1p2)))

!************************************************   |Ms1C321|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(0)

!************************************************   2 Re(Ms1C321 Ms2C321^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**4*(k1k3 + k2k3)*Nc**2*(-1 + Nc**2)*(p1k3 - p2k3)
     &*(p1k1 + p1k3 - p2k1 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k2*k1k3*(M2 + p1p2)**2))

!************************************************   2 Re(Ms1C321 Ms3C321^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GGG*GhGhG**3*(k1k3 - k2k3)*Nc**2*(-1 + Nc**2)*(p1k3 - p2k3)
     &*(p1k1 + p1k3 - p2k1 - p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k3*k2k3*(M2 + p1p2)**2))

!************************************************   |Ms2C321|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**4*Nc**2*(-1 + Nc**2)*(p1k3 - p2k3)**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(16*k1k2*(M2 + p1p2)**2))

!************************************************   2 Re(Ms2C321 Ms3C321^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GGG*GhGhG**3*Nc**2*(-1 + Nc**2)*(p1k3 - p2k3)*(-2*k1k3
     &*p1k2 + k1k2*p1k3 + 2*p1k3*p2k1 + 2*k1k3*p2k2 + 2*p1k3*p2k2 + k2k3*(p1k1 - p1k2 - p2k1 + p2k2) 
     &- k1k2*p2k3 - 2*p1k1*p2k3 - 2*p1k2*p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k2*k2k3*(M2 + p1p2)**2))

!************************************************   |Ms3C321|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GGG**2*GhGhG**2*Nc**2*(-1 + Nc**2)*(k1k3**2*(M2 - p1p2) + k2k3
     &*(p1k1 - p1k2 - p2k1 + p2k2)*(p1k3 - p2k3) + k1k3*(-2*p1k3*p1p2 + p1k3*p2k1 - 2*p1p2*p2k3 
     &- p2k1*p2k3 + p1k1*(-p1k3 + p2k3) + 2*M2*(p1k3 + p2k3)) + 2*(p1k3**2*p2k2 - p1k3*(p1k2 + 2*p1p2 
     &+ p2k2)*p2k3 + p1k2*p2k3**2 + M2*(p1k3**2 + p2k3**2)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(16*k2k3**2*(M2 + p1p2)**2))


!##############################################################################################################################
!#################################################    |S6|^2   ################################################################   
!##############################################################################################################################


!************************************************   |MvC231|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(0)

!************************************************   2 Re(MvC231 MtC231^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**2*(2 - 3*Nc**2 + Nc**4)*(k2k3 - 2*p1k2)*(k2k3 
     &- p1k2 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32
     &*k2k3**2*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(MvC231 MbC231^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**2*(2 - 3*Nc**2 + Nc**4)*(k1k2 - 2*p1k2)*(k1k2 
     &- p1k2 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(32*k2k3**2*p1k1))

!************************************************   2 Re(MvC231 Ms1C231^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   2 Re(MvC231 Ms2C231^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   2 Re(MvC231 Ms3C231^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(0)

!************************************************   |MtC231|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**2*(-1 + Nc**2)**2*(k2k3 - 2*p1k2)**2*(k2k3 + M2 - p1k2 
     &- p1k3 - p1p2 + p2k2 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(32*k2k3**2*(-k2k3 + p1k2 + p1k3)**2))

!************************************************   2 Re(MtC231 MbC231^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**2*(-1 + Nc**2)*(k2k3 - 2*p1k2)*(-k1k2 - k1k3 
     &- 2*M2 + p1k1 + 2*p1k2 + 2*p1k3 + 2*p1p2 - p2k1)*(-k1k2 + p1k2 - p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**4)/(64*k2k3**2*p1k1*(k2k3 - p1k2 - p1k3)))

!************************************************   2 Re(MtC231 Ms1C231^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k2k3 - 2*p1k2)*(p1k2 - p2k2)
     &*(k1k2 + k1k3 + k2k3 - p1k1 - p1k3 + p2k1 + p2k3)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**3)/(128*k1k3*k2k3*(k2k3 - p1k2 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(MtC231 Ms2C231^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k2k3 - 2*p1k2)*(p1k1 + p1k2 
     &- p2k1 - p2k2)*(k2k3 - p1k2 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2
     &*k2k3*(k2k3 - p1k2 - p1k3)*(M2 + p1p2)))

!************************************************   2 Re(MtC231 Ms3C231^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GGG*GhGhG**2*Nc**2*(-1 + Nc**2)*(k2k3 - 2*p1k2)*(2*M2
     &*p1k2 - p1k1*p1k2 - 2*p1k2*p1p2 + p1k2*p2k1 + 2*M2*p2k2 + p1k1*p2k2 - 2*p1k3*p2k2 - 2*p1p2*p2k2 
     &- p2k1*p2k2 + 2*p1k2*p2k3 + k2k3*(p1k1 + p1k2 - p1k3 - p2k1 - p2k2 + p2k3) + k1k2*(2*M2 - p1k2 
     &- p1k3 - 2*p1p2 + p2k2 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(64*k2k3**2*(k2k3 - p1k2 - p1k3)*(M2 + p1p2)))

!************************************************   |MbC231|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**2*(-1 + Nc**2)**2*(M2 - p1k1)*(k1k2 - p1k2 + p2k2)*
     &*2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(16*k2k3**2*p1k1**2))

!************************************************   2 Re(MbC231 Ms1C231^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**3*Nc**2*(-1 + Nc**2)*(k1k3 - 2*(p1k1 + p1k3))*(p1k2 
     &- p2k2)*(-k1k2 + p1k2 - p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k3*k2k3*p1k1*(M2 + p1p2)))

!************************************************   2 Re(MbC231 Ms2C231^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GhGhG**3*Nc**2*(-1 + Nc**2)*(k1k2 - 2*p1k2)*(p1k1 + p1k2 
     &- p2k1 - p2k2)*(k1k2 - p1k2 + p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(128*k1k2*k2k3*p1k1*(M2 + p1p2)))

!************************************************   2 Re(MbC231 Ms3C231^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GGG*GhGhG**2*Nc**2*(-1 + Nc**2)*(k1k2 - p1k2 + p2k2)*(
     &-2*(2*M2 + p1k3)*p2k2 + k1k3*(-p1k2 + p2k2) + 2*p1k1*(p1k2 + p2k2) + k1k2*(-2*M2 + p1k3 + 2*p1p2 
     &- p2k3) + 2*p1k2*(2*p1p2 - 2*p2k1 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(64*k2k3**2*p1k1*(M2 + p1p2)))

!************************************************   |Ms1C231|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GhGhG**4*Nc**2*(-1 + Nc**2)*(p1k2 - p2k2)**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &**2)/(16*k1k3*(M2 + p1p2)**2))

!************************************************   2 Re(Ms1C231 Ms2C231^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GhGhG**4*(k1k2 + k2k3)*Nc**2*(-1 + Nc**2)*(p1k2 - p2k2)
     &*(p1k1 + p1k2 - p2k1 - p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k2*k1k3*(M2 + p1p2)**2))

!************************************************   2 Re(Ms1C231 Ms3C231^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble(-(GGG*GhGhG**3*Nc**2*(-1 + Nc**2)*(p1k2 - p2k2)*(k1k3*(p1k2 
     &- p2k2) - 2*(p1k1 + p1k3)*p2k2 + 2*k1k2*(-p1k3 + p2k3) + k2k3*(p1k1 - p1k3 - p2k1 + p2k3) 
     &+ 2*p1k2*(p2k1 + p2k3))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k3*k2k3*(M2 + p1p2)**2))

!************************************************   |Ms2C231|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(0)

!************************************************   2 Re(Ms2C231 Ms3C231^ast)     ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+2*dble((GGG*GhGhG**3*(k1k2 - k2k3)*Nc**2*(-1 + Nc**2)*(p1k2 - p2k2)
     &*(p1k1 + p1k2 - p2k1 - p2k2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(64*k1k2*k2k3*(M2 + p1p2)**2))

!************************************************   |Ms3C231|^2   ************************************************

      stsT2xx_M2RealGhost=stsT2xx_M2RealGhost+(-(GGG**2*GhGhG**2*Nc**2*(-1 + Nc**2)*(k1k2**2*(M2 - p1p2) + k1k2
     &*(-2*p1k2*p1p2 + p1k2*p2k1 - 2*p1p2*p2k2 - p2k1*p2k2 + p1k1*(-p1k2 + p2k2) + 2*M2*(p1k2 + p2k2)) 
     &+ k2k3*(p1k2 - p2k2)*(p1k1 - p1k3 - p2k1 + p2k3) + 2*(p1k3*p2k2**2 + M2*(p1k2**2 + p2k2
     &**2) + p1k2**2*p2k3 - p1k2*p2k2*(p1k3 + 2*p1p2 + p2k3)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(16*k2k3**2*(M2 + p1p2)**2))


      end


