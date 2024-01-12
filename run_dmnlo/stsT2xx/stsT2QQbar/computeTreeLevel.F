!       double precision function computeTreeLevel2QG(vcAL1,vcAR1,vcBL1,vcBR1,vcCL1,vcCR1,vcDL1,vcDR1,vcAL2,vcAR2,vcBL2,vcBR2,vcCL2,vcCR2,vcAL3,vcAR3,vcBL3,vcBR3,vcCL3,vcCR3)      
      double precision function computeTreeLevel2QG(pvcAL1,pvcAR1 ,pvcBL1,pvcBR1,pvcCL1,pvcCR1,pvcDL1,pvcDR1,pvcAL2,pvcAR2,pvcBL2,pvcBR2,pvcCL2,pvcCR2,pvcAL3,pvcAR3,pvcBL3,pvcBR3,pvcCL3,pvcCR3)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      double precision colorf
      double complex DenSV
      double complex vcAL1,vcAR1,vcBL1,vcBR1,vcCL1,vcCR1,vcDL1,vcDR1
      double complex vcAL2,vcAR2,vcBL2,vcBR2,vcCL2,vcCR2
      double complex vcAL3,vcAR3,vcBL3,vcBR3,vcCL3,vcCR3
      double complex pvcAL1,pvcAR1,pvcBL1,pvcBR1,pvcCL1,pvcCR1,pvcDL1,pvcDR1
      double complex pvcAL2,pvcAR2,pvcBL2,pvcBR2,pvcCL2,pvcCR2
      double complex pvcAL3,pvcAR3,pvcBL3,pvcBR3,pvcCL3,pvcCR3

      DenSV = sman


      vcAR1 = pvcAR1/(4*Pi)**2
      vcAL1 = pvcAL1/(4*Pi)**2

      vcAR2 = pvcAR2/(4*Pi)**2
      vcAL2 = pvcAL2/(4*Pi)**2

      vcAR3 = pvcAR3/(4*Pi)**2
      vcAL3 = pvcAL3/(4*Pi)**2

      vcBL1 = pvcBL1/(4*Pi)**2
      vcBR1 = pvcBR1/(4*Pi)**2

      vcBL2 = pvcBL2/(4*Pi)**2
      vcBR2 = pvcBR2/(4*Pi)**2

      vcBL3 = pvcBL3/(4*Pi)**2
      vcBR3 = pvcBR3/(4*Pi)**2

      vcCL1 = pvcCL1/(4*Pi)**2
      vcCR1 = pvcCL1/(4*Pi)**2

      vcCL2 = pvcCL2/(4*Pi)**2
      vcCR2 = pvcCR2/(4*Pi)**2

      vcCL3 = pvcCL3/(4*Pi)**2
      vcCR3 = pvcCL3/(4*Pi)**2

      vcDL1 = pvcDL1/(4*Pi)**2
      vcDR1 = pvcDR1/(4*Pi)**2


      colorf=(Nc**2-1)/4d0

      computeTreeLevel2QG =  - 2d0*colorf*dreal(((sqrtS**4*((-(muu4**2*vcAL1) + muu4**4*vcAL1 + 4*tred*vcAL1 - 
     -           4*muu2**2*tred*vcAL1 - 4*muu4**2*tred*vcAL1 + 
     -           4*tred**2*vcAL1 + muu3**5*sqrtS*vcAL2 + 
     -           muu3**3*sqrtS*((-2 + 3*muu2**2 + 2*muu4**2 - 4*tred)*vcAL2 - 
     -              muu2**2*vcAL3) + 
     -           muu3*(sqrtS*(-1 + 2*muu2**2 + muu4**2 - 2*tred)*
     -               ((-1 + muu2**2 + muu4**2 - 2*tred)*vcAL2 - muu2**2*vcAL3)
     -               + 2*(-1 + 2*muu2**2)*muu4*vcAR1) + muu4*sqrtS*vcAR2 - 
     -           muu2**2*muu4*sqrtS*vcAR2 - 2*muu4**3*sqrtS*vcAR2 + 
     -           muu2**2*muu4**3*sqrtS*vcAR2 + muu4**5*sqrtS*vcAR2 + 
     -           4*muu4*sqrtS*tred*vcAR2 - 2*muu2**2*muu4*sqrtS*tred*vcAR2 - 
     -           4*muu4**3*sqrtS*tred*vcAR2 + 4*muu4*sqrtS*tred**2*vcAR2 + 
     -           muu3**4*(vcAL1 + muu4*sqrtS*vcAR2) + 
     -           muu2**2*muu4*sqrtS*vcAR3 - muu2**2*muu4**3*sqrtS*vcAR3 + 
     -           2*muu2**2*muu4*sqrtS*tred*vcAR3 + 
     -           2*muu1**4*muu4*sqrtS*(vcAR2 + vcAR3) + 
     -           muu3**2*((-1 + 2*muu4**2 - 4*tred)*vcAL1 + 
     -              muu4*sqrtS*((-2 + muu2**2 + 2*muu4**2 - 4*tred)*vcAR2 - 
     -                 muu2**2*vcAR3)) + 
     -           muu1**2*(-(muu3*sqrtS*vcAL2) + muu3**3*sqrtS*vcAL2 + 
     -              muu3*muu4**2*sqrtS*vcAL2 - muu3*sqrtS*vcAL3 + 
     -              muu3**3*sqrtS*vcAL3 + muu3*muu4**2*sqrtS*vcAL3 + 
     -              4*muu3*muu4*vcAR1 - 3*muu4*sqrtS*vcAR2 + 
     -              3*muu3**2*muu4*sqrtS*vcAR2 + 3*muu4**3*sqrtS*vcAR2 + 
     -              2*muu2**2*(2*vcAL1 + 
     -                 sqrtS*(muu3*(vcAL2 + vcAL3) + muu4*(vcAR2 - vcAR3))) - 
     -              muu4*sqrtS*vcAR3 + muu3**2*muu4*sqrtS*vcAR3 + 
     -              muu4**3*sqrtS*vcAR3 - 
     -              2*tred*(2*vcAL1 + 
     -                 sqrtS*(muu3*(vcAL2 + vcAL3) + muu4*(3*vcAR2 + vcAR3)))))*
     -         (-GS) + 
     -        (muu4**5*sqrtS*vcAL2 + 
     -           muu4**3*sqrtS*((-2 + 3*muu1**2 + muu2**2 - 4*tred)*vcAL2 + 
     -              (muu1**2 - muu2**2)*vcAL3) + 
     -           muu4*sqrtS*(-1 + 2*muu1**2 - 2*tred)*
     -            ((-1 + muu1**2 + muu2**2 - 2*tred)*vcAL2 + 
     -              (muu1**2 - muu2**2)*vcAL3) + muu4**4*vcAR1 + 
     -           4*(muu1**2*(muu2**2 - tred) + tred*(1 - muu2**2 + tred))*
     -            vcAR1 + muu3**4*(muu4*sqrtS*vcAL2 + vcAR1) - 
     -           muu4**2*(vcAR1 + 4*tred*vcAR1) + 
     -           muu3**2*(2*muu4**3*sqrtS*vcAL2 + 
     -              muu4*sqrtS*((-2 + 3*muu1**2 + muu2**2 - 4*tred)*vcAL2 + 
     -                 (muu1**2 - muu2**2)*vcAL3) + 2*muu4**2*vcAR1 - 
     -              (1 + 4*tred)*vcAR1) + muu3**5*sqrtS*vcAR2 + 
     -           muu3**3*sqrtS*((-2 + muu1**2 + 3*muu2**2 + 2*muu4**2 - 4*tred)*
     -               vcAR2 + (muu1**2 - muu2**2)*vcAR3) + 
     -           muu3*(2*(-1 + 2*muu1**2 + 2*muu2**2)*muu4*vcAL1 + 
     -              muu4**4*sqrtS*vcAR2 + 
     -              muu4**2*sqrtS*
     -               ((-2 + muu1**2 + 3*muu2**2 - 4*tred)*vcAR2 + 
     -                 (muu1**2 - muu2**2)*vcAR3) + 
     -              sqrtS*(-1 + 2*muu2**2 - 2*tred)*
     -               ((-1 + muu1**2 + muu2**2 - 2*tred)*vcAR2 + 
     -                (muu1**2 - muu2**2)*vcAR3)))*(-GS))
     -       *dconjg(SfSfg(isf1,itt1,igen1,isf2,itt2,igen2))*
     -      SfSfg(isf1,itt1,igen1,isf2,itt2,igen2))/
     -    (dconjg(DenSV)*DenSV)))     

      end




      double precision function FormCalc_QGV_Triangle1() ! Quark-Quark-Gluon
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      call myCget(0d0,sman,0d0,0d0,0d0,0d0)
      call myBget(sman,0d0,0d0)

      FormCalc_QGV_Triangle1 = 2*dble(-1d0/3d0*(gs**6*(m1**4- tman*uman)*(-xb0 + 2*xc00 + sman*(xc0 + xc1 + xc2)))/(Pi**2*sman**2))

      end

      double precision function FormCalc_QGV_Triangle2(fType,fGen) ! Gluino-Squark-Squark
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      integer fType,fGen,k

      FormCalc_QGV_Triangle2 = 0d0

      do k=1,2

      call myCget(0d0,sman,0d0,Mgl2,MSf2(k,fType,fGen),MSf2(k,fType,fGen))

      FormCalc_QGV_Triangle2 = FormCalc_QGV_Triangle2 + 2*dble((gs**6*(m1**4 - tman*uman)*xc00*(USf(k,1,fType,fGen)*USfC(k,1,fType,fGen) + USf(k,2,fType,fGen)*USfC(k,2,fType,fGen)))/(3*Pi**2*sman**2))

      enddo

      end

      double precision function FormCalc_QGV_Triangle3() ! Quark-Gluon-Gluon
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      call myCget(sman,0d0,0d0,0d0,0d0,0d0)
      call myBget(sman,0d0,0d0)

      FormCalc_QGV_Triangle3 = 2*dble((-3*gs**6*(m1**4- tman*uman)*(xb0 + 2*xc00 + sman*(xc0 + xc2)))/(Pi**2*sman**2))

      end

      double precision function FormCalc_QGV_Triangle4(fType,fGen) ! Squark-Gluino-Gluino
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      integer fType,fGen,k

      FormCalc_QGV_Triangle4 = 0d0

      do k=1,2

      call myCget(sman,0d0,0d0,Mgl2,Mgl2,MSf2(k,fType,fGen))
      call myBget(sman,Mgl2,Mgl2)

      FormCalc_QGV_Triangle4 = FormCalc_QGV_Triangle4 + 2*dble((-3*gs**6*(m1**4 - tman*uman)*(xb0 - Mgl**2*xc0 - 2*xc00 + xc0*MSf2(k,fType,fGen))*
     &    (USf(k,1,fType,fGen)*USfC(k,1,fType,fGen) + USf(k,2,fType,fGen)*USfC(k,2,fType,fGen)))/(2*Pi**2*sman**2))

      enddo

      end


! Results for the boxes produced with FormCalc 9.9. These are in perfect agreement with the results below.

      double precision function FormCalcBox1()
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"


      call myDget(sman,0d0,tman,m1**2,0d0,m1**2,0d0,0d0,0d0,m1**2) 
      call myCget(sman,0d0,0d0,0d0,0d0,0d0)

      FormCalcBox1 = 2*dble((-7*gs**6*(m1**4 - tman*uman)*(4*xc0 - 4*tman*xd0 - 2*xd00 + tman*xd23 - 4*tman*xd3 + m1**2*(4*xd0 - xd23 - xd33) - tman*xd33))/(12.*Pi**2*sman))     

      end      

      double precision function FormCalcBox1u()
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      call myDget(sman,0d0,uman,m1**2,0d0,m1**2,0d0,0d0,0d0,m1**2)
      call myCget(sman,0d0,0d0,0d0,0d0,0d0)

      FormCalcBox1u = 2*dble((gs**6*(-m1**4 + tman*uman)*((-2*(-m1**2 + uman)*xd0)/3. - xd00/3. + (2*(xc0 - uman*xd3))/3. + ((-m1**2 + uman)*xd23 + (-m1**2 - uman)*xd33)/6.))/(Pi**2*sman))     

      end    

      double precision function FormCalcBox2(fType,fGen)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      integer fType,fGen,k

      FormCalcBox2 = 0d0

      do k=1,2

      call myDget(sman,m1**2,tman,0d0,m1**2,0d0,Mgl2,Mgl2,Mf2(itt1,igen1),MSf2(k,fType,fGen))
      call myCget(sman,0d0,0d0,Mgl2,Mgl2,MSf2(k,fType,fGen))

      FormCalcBox2 = FormCalcBox2 + 2*dble( (7*gs**6*(-m1**4 + tman*uman)*(xc0*(USf(isf1,2,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,1,fType,fGen) + 
     &         USf(isf1,1,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,2,fType,fGen)) + 
     &      (-m1**2 + tman)*xd3*(USf(isf1,2,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,1,fType,fGen) + 
     &         USf(isf1,1,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,2,fType,fGen)) - 
     &      Mgl**2*xd0*(USf(isf1,1,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,1,fType,fGen) + 
     &         USf(isf1,2,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,2,fType,fGen) + 
     &         ((m1**2 - Mf2(itt1,igen1))*(USf(isf1,2,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,1,fType,fGen) + 
     &             USf(isf1,1,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,2,fType,fGen)))/Mgl**2) + 
     &      xd2*(Mgl*Mf(itt1,igen1)*USf(isf1,1,itt1,igen1)*USfC(isf1,2,itt1,igen1)*(USf(k,1,fType,fGen)*USfC(k,1,fType,fGen) + USf(k,2,fType,fGen)*USfC(k,2,fType,fGen)) + 
     &         Mgl*(Mf(itt1,igen1)*USf(isf1,2,itt1,igen1)*USfC(isf1,1,itt1,igen1)*(USf(k,1,fType,fGen)*USfC(k,1,fType,fGen) + USf(k,2,fType,fGen)*USfC(k,2,fType,fGen)) + 
     &            Mgl*(-(USf(isf1,1,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,1,fType,fGen)) - 
     &               USf(isf1,2,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,2,fType,fGen) - 
     &               ((m1**2 + Mf2(itt1,igen1))*(USf(isf1,2,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,1,fType,fGen) + 
     &                    USf(isf1,1,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,2,fType,fGen)))/Mgl**2)))))/(6.*Pi**2*sman))     

      enddo

      end    

      double precision function FormCalcBox2u(fType,fGen)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      integer fType,fGen,k

      FormCalcBox2u = 0d0

      do k=1,2

      call myDget(sman,m1**2,uman,0d0,m1**2,0d0,Mgl2,Mgl2,Mf2(itt1,igen1),MSf2(k,fType,fGen))
      call myCget(sman,0d0,0d0,Mgl2,Mgl2,MSf2(k,fType,fGen))

      FormCalcBox2u = FormCalcBox2u + 2*dble((gs**6*(-m1**4 + tman*uman)*(xc0*(USf(isf1,1,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,1,fType,fGen) + 
     -         USf(isf1,2,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,2,fType,fGen)) + 
     -      (-m1**2 + uman)*xd3*(USf(isf1,1,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,1,fType,fGen) + 
     -         USf(isf1,2,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,2,fType,fGen)) - 
     -      Mgl**2*xd0*(USf(isf1,2,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,1,fType,fGen) + 
     -         USf(isf1,1,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,2,fType,fGen) + 
     -         ((m1**2 - Mf2(itt1,igen1))*(USf(isf1,1,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,1,fType,fGen) + 
     -             USf(isf1,2,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,2,fType,fGen)))/Mgl**2) + 
     -      xd2*(Mgl*Mf(itt1,igen1)*USf(isf1,1,itt1,igen1)*USfC(isf1,2,itt1,igen1)*(USf(k,1,fType,fGen)*USfC(k,1,fType,fGen) + USf(k,2,fType,fGen)*USfC(k,2,fType,fGen)) + 
     -         Mgl*(Mf(itt1,igen1)*USf(isf1,2,itt1,igen1)*USfC(isf1,1,itt1,igen1)*(USf(k,1,fType,fGen)*USfC(k,1,fType,fGen) + USf(k,2,fType,fGen)*USfC(k,2,fType,fGen)) + 
     -            Mgl*(-(USf(isf1,2,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,1,fType,fGen)) - 
     -               USf(isf1,1,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,2,fType,fGen) - 
     -               ((m1**2 + Mf2(itt1,igen1))*(USf(isf1,1,itt1,igen1)*USf(k,1,fType,fGen)*USfC(isf1,1,itt1,igen1)*USfC(k,1,fType,fGen) + 
     -                    USf(isf1,2,itt1,igen1)*USf(k,2,fType,fGen)*USfC(isf1,2,itt1,igen1)*USfC(k,2,fType,fGen)))/Mgl**2)))))/(3.*Pi**2*sman))     

      enddo

      end 
     

      double precision function FormCalc_SquarkGluon_1() ! Gluon-squark loop at p1
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      call myBget(m1**2,0d0,m1**2)

      FormCalc_SquarkGluon_1 = 2*dble((7*gs**6*(m1**4 - tman*uman)*(2*xb0 + xb1))/(12*Pi**2*sman**2))

      end

      double precision function FormCalc_SquarkGluon_2() !Gluon-squark loop at p2
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      call myBget(m1**2,0d0,m1**2)

      FormCalc_SquarkGluon_2 = 2*dble((7*gs**6*(m1**4 - tman*uman)*(2*xb0 +xb1))/(12*Pi**2*sman**2))

      end


      double precision function FormCalc_SquarkGluon_5() ! Gluon-exchange (Gluon-Squark-Squark)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      call myCget(m1**2,sman,m1**2,0d0,m1**2,m1**2)

      FormCalc_SquarkGluon_5 = 2*dble((gs**6*(m1**4 - tman*uman)*(4*xc00+ 4*m1**2*(xc0 + 2*xc1 + xc11 + 2*xc12 + 2*xc2 + xc22) - 
     &      sman*(2*xc0 + 3*xc1 + xc11 + 2*xc12 + 3*xc2 + xc22)))/(6*Pi**2*sman**2))

      end

      double precision function FormCalc_SquarkGluon_6() ! Squark-exchange (Squark-Gluon-Gluon)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      call myCget(sman,m1**2,m1**2,0d0,0d0,m1**2)
      call myBget(sman,0d0,0d0)

      FormCalc_SquarkGluon_6 = 2*dble((-3*gs**6*(m1**4 - tman*uman)*(4*xb0 + sman*(4*xc0 + 4*xc2 + xc22) - 4*(xc00 + m1**2*(2*xc2 + xc22))))/(4*Pi**2*sman**2))

      end

      double precision function FormCalc_SquarkGluon_7() ! Gluino-exchange (Gluino-Quark-Quark)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      call myCget(m1**2,sman,m1**2,Mgl**2,Mf2(itt1,igen1),Mf2(itt1,igen1))
      call myBget(sman,Mf2(itt1,igen1),Mf2(itt1,igen1))

      FormCalc_SquarkGluon_7 = 2*dble((gs**6*(m1**4 - tman*uman)*(USf(isf1,1,itt1,igen1)*
     -      ((xb0 + 2*Mgl**2*xc0 + m1**2*xc1 + Mgl**2*xc1 + m1**2*xc2 + Mgl**2*xc2 + (xc1 + xc2)*Mf2(itt1,igen1))*
     -          USfC(isf1,1,itt1,igen1) - 2*Mgl*(xc0 + xc1 + xc2)*Mf(itt1,igen1)*USfC(isf1,2,itt1,igen1)) + 
     -      USf(isf1,2,itt1,igen1)*(-2*Mgl*(xc0 + xc1 + xc2)*Mf(itt1,igen1)*USfC(isf1,1,itt1,igen1) + 
     -         (xb0 + 2*Mgl**2*xc0 + m1**2*xc1 + Mgl**2*xc1 + m1**2*xc2 + Mgl**2*xc2 + (xc1 + xc2)*Mf2(itt1,igen1))*
     -          USfC(isf1,2,itt1,igen1))))/(3*Pi**2*sman**2))

      end

      double precision function FormCalc_SquarkGluon_8() ! Quark-exchange (Squark-Gluon-Gluon)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      call myCget(sman,m1**2,m1**2,Mgl2,Mgl2,Mf2(itt1,igen1))
      call myBget(sman,Mgl2,Mgl2)

      FormCalc_SquarkGluon_8 = 2*dble((3*gs**6*(m1**4 - tman*uman)*(USf(isf1,1,itt1,igen1)*
     -       ((-xb0 + (m1**2 + Mgl**2)*(xc0 + xc2) + (-xc0 + xc2)*Mf2(itt1,igen1))*USfC(isf1,1,itt1,igen1) - 
     -         2*Mgl*xc2*Mf(itt1,igen1)*USfC(isf1,2,itt1,igen1)) + 
     -      USf(isf1,2,itt1,igen1)*(-2*Mgl*xc2*Mf(itt1,igen1)*USfC(isf1,1,itt1,igen1) + 
     -         (-xb0 + (m1**2 + Mgl**2)*(xc0 + xc2) + (-xc0 + xc2)*Mf2(itt1,igen1))*USfC(isf1,2,itt1,igen1))))/
     -  (Pi**2*sman**2))

      end