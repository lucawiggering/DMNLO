      double precision function stsT2QQbar_M2Vertex2sqG()
      implicit none

#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      double complex DenSV
      double complex VSfSf1m2,VSfSf1p2
      double precision colorf,Nf

!Propagator
      DenSV = sman

! Number of massless quarks
      Nf = 4d0

      call stsT2xx_calc_partial_amplitudes_ASqSqG(VSfSf1m2, VSfSf1p2, isf1,itt1,igen1,isf2,itt2,igen2,m1**2,m1**2,sman)


c----------------------------------------------s-channel Gluon-----------------------------------------------------------
      colorf = 1d0/4d0*(Nc**2-1)

      stsT2QQbar_M2Vertex2sqG =  - Nf*2*colorf*dble(((sqrtS**4*((-GS)*
     &         ((-(muu4**2*VSfSf1m2) + muu4**4*VSfSf1m2 + 4*tred*VSfSf1m2 -
     &              4*muu2**2*tred*VSfSf1m2 - 4*muu4**2*tred*VSfSf1m2 +
     &              4*tred**2*VSfSf1m2 + muu4**2*VSfSf1p2 - muu4**4*VSfSf1p2 +
     &              2*muu4**2*tred*VSfSf1p2 + muu3**4*(VSfSf1m2 + VSfSf1p2) +
     &              muu1**2*(4*muu2**2*VSfSf1m2 - 4*tred*VSfSf1m2 -
     &                 2*muu4**2*VSfSf1p2) +
     &              muu3**2*((-1 + 2*muu4**2 - 4*tred)*VSfSf1m2 +
     &                 (-1 + 2*muu2**2 - 2*tred)*VSfSf1p2))*
     &            (-GS) + 2*muu3*muu4*((-1 + 2*muu1**2 + 2*muu2**2)*VSfSf1m2 +
     &              (muu1**2 - muu2**2)*VSfSf1p2)*(-GS)) +
     &        (-GS)*(2*muu3*muu4*((-1 + 2*muu1**2 + 2*muu2**2)*VSfSf1m2 +
     &              (muu1**2 - muu2**2)*VSfSf1p2)*(-GS) +
     &           (-(muu4**2*VSfSf1m2) + muu4**4*VSfSf1m2 + 4*tred*VSfSf1m2 -
     &              4*muu2**2*tred*VSfSf1m2 - 4*muu4**2*tred*VSfSf1m2 +
     &              4*tred**2*VSfSf1m2 + muu4**2*VSfSf1p2 - muu4**4*VSfSf1p2 +
     &              2*muu4**2*tred*VSfSf1p2 + muu3**4*(VSfSf1m2 + VSfSf1p2) +
     &              muu1**2*(4*muu2**2*VSfSf1m2 - 4*tred*VSfSf1m2 -
     &                 2*muu4**2*VSfSf1p2) +
     &              muu3**2*((-1 + 2*muu4**2 - 4*tred)*VSfSf1m2 +
     &                 (-1 + 2*muu2**2 - 2*tred)*VSfSf1p2))*(-GS)))*
     &      dconjg(SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)))/
     &    (dconjg(DenSV)*DenSV)))

      end


! Below are the results for the squark-gluon-vertex correction obtained with FormCalc 9.9. 
! They are in perfect agreement with the results above.

!      double precision function FormCalc_SquarkGluon_1() ! Gluon-squark loop at p1
!      implicit none
!#include "../common/stsT2xx_Kinematics.h"
!#include "../common/stsT2xx_Model.h"
!#include "../../util/DMNLO_LoopIntegrals.h"
!
!      call myBget(m1**2,0d0,m1**2)
!
!      FormCalc_SquarkGluon_1 = 2*dble((7*gs**6*(m1**4 - tman*uman)*(2*xb0 + xb1))/(12*Pi**2*sman**2))
!
!      end
!
!      double precision function FormCalc_SquarkGluon_2() !Gluon-squark loop at p2
!      implicit none
!#include "../common/stsT2xx_Kinematics.h"
!#include "../common/stsT2xx_Model.h"
!#include "../../util/DMNLO_LoopIntegrals.h"
!
!      call myBget(m1**2,0d0,m1**2)
!
!      FormCalc_SquarkGluon_2 = 2*dble((7*gs**6*(m1**4 - tman*uman)*(2*xb0 +xb1))/(12*Pi**2*sman**2))
!
!      end
!
!
!      double precision function FormCalc_SquarkGluon_5() ! Gluon-exchange (Gluon-Squark-Squark)
!      implicit none
!#include "../common/stsT2xx_Kinematics.h"
!#include "../common/stsT2xx_Model.h"
!#include "../../util/DMNLO_LoopIntegrals.h"
!
!      call myCget(m1**2,sman,m1**2,0d0,m1**2,m1**2)
!
!      FormCalc_SquarkGluon_5 = 2*dble((gs**6*(m1**4 - tman*uman)*(4*xc00+ 4*m1**2*(xc0 + 2*xc1 + xc11 + 2*xc12 + 2*xc2 + xc22) - 
!     &      sman*(2*xc0 + 3*xc1 + xc11 + 2*xc12 + 3*xc2 + xc22)))/(6*Pi**2*sman**2))
!
!      end
!
!      double precision function FormCalc_SquarkGluon_6() ! Squark-exchange (Squark-Gluon-Gluon)
!      implicit none
!#include "../common/stsT2xx_Kinematics.h"
!#include "../common/stsT2xx_Model.h"
!#include "../../util/DMNLO_LoopIntegrals.h"
!
!      call myCget(sman,m1**2,m1**2,0d0,0d0,m1**2)
!      call myBget(sman,0d0,0d0)
!
!      FormCalc_SquarkGluon_6 = 2*dble((-3*gs**6*(m1**4 - tman*uman)*(4*xb0 + sman*(4*xc0 + 4*xc2 + xc22) - 4*(xc00 + m1**2*(2*xc2 + xc22))))/(4*Pi**2*sman**2))
!
!      end
!
!      double precision function FormCalc_SquarkGluon_7() ! Gluino-exchange (Gluino-Quark-Quark)
!      implicit none
!#include "../common/stsT2xx_Kinematics.h"
!#include "../common/stsT2xx_Model.h"
!#include "../../util/DMNLO_LoopIntegrals.h"
!
!      call myCget(m1**2,sman,m1**2,Mgl**2,Mf2(itt1,igen1),Mf2(itt1,igen1))
!      call myBget(sman,Mf2(itt1,igen1),Mf2(itt1,igen1))
!
!      FormCalc_SquarkGluon_7 = 2*dble((gs**6*(m1**4 - tman*uman)*(USf(isf1,1,itt1,igen1)*
!     -      ((xb0 + 2*Mgl**2*xc0 + m1**2*xc1 + Mgl**2*xc1 + m1**2*xc2 + Mgl**2*xc2 + (xc1 + xc2)*Mf2(itt1,igen1))*
!     -          USfC(isf1,1,itt1,igen1) - 2*Mgl*(xc0 + xc1 + xc2)*Mf(itt1,igen1)*USfC(isf1,2,itt1,igen1)) + 
!     -      USf(isf1,2,itt1,igen1)*(-2*Mgl*(xc0 + xc1 + xc2)*Mf(itt1,igen1)*USfC(isf1,1,itt1,igen1) + 
!     -         (xb0 + 2*Mgl**2*xc0 + m1**2*xc1 + Mgl**2*xc1 + m1**2*xc2 + Mgl**2*xc2 + (xc1 + xc2)*Mf2(itt1,igen1))*
!     -          USfC(isf1,2,itt1,igen1))))/(3*Pi**2*sman**2))
!
!      end
!
!      double precision function FormCalc_SquarkGluon_8() ! Quark-exchange (Squark-Gluon-Gluon)
!      implicit none
!#include "../common/stsT2xx_Kinematics.h"
!#include "../common/stsT2xx_Model.h"
!#include "../../util/DMNLO_LoopIntegrals.h"
!
!      call myCget(sman,m1**2,m1**2,Mgl2,Mgl2,Mf2(itt1,igen1))
!      call myBget(sman,Mgl2,Mgl2)
!
!      FormCalc_SquarkGluon_8 = 2*dble((3*gs**6*(m1**4 - tman*uman)*(USf(isf1,1,itt1,igen1)*
!     -       ((-xb0 + (m1**2 + Mgl**2)*(xc0 + xc2) + (-xc0 + xc2)*Mf2(itt1,igen1))*USfC(isf1,1,itt1,igen1) - 
!     -         2*Mgl*xc2*Mf(itt1,igen1)*USfC(isf1,2,itt1,igen1)) + 
!     -      USf(isf1,2,itt1,igen1)*(-2*Mgl*xc2*Mf(itt1,igen1)*USfC(isf1,1,itt1,igen1) + 
!     -         (-xb0 + (m1**2 + Mgl**2)*(xc0 + xc2) + (-xc0 + xc2)*Mf2(itt1,igen1))*USfC(isf1,2,itt1,igen1))))/
!     -  (Pi**2*sman**2))
!
!      end