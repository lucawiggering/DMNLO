      double precision function stsT2QQbar_M2Prop()
      implicit none

#include "../common/stsT2xx_Kinematics.h"

      double complex stsT2xx_PiGluonLoopT
      double precision stst2qqbar_M2Tree,Nf

! Number of effectively massless quarks
      Nf = 4d0

      stsT2QQbar_M2Prop = Nf*2*dble(-stsT2xx_PiGluonLoopT(sman)/sman)*stst2qqbar_M2Tree()

      end