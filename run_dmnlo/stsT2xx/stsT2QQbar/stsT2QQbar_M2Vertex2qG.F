      double precision function stsT2QQbar_M2Vertex2qG()
      implicit none

#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      integer k,fType,fGen
      double complex g0,g0L,g0R,g1,g1L,g1R,g2,g2L,g2R,g3L,g3R
      double precision M0l,M1l,M2l
      double complex DenSV
      double complex vcAL1,vcAR1,vcBL1,vcBR1,vcCL1,vcCR1,vcDL1,vcDR1
      double complex vcAL2,vcAR2,vcBL2,vcBR2,vcCL2,vcCR2
      double complex vcAL3,vcAR3,vcBL3,vcBR3,vcCL3,vcCR3
	double precision colorf, colori

!Gluon propagator
      DenSV = sman

! initilize with zero
      stsT2QQbar_M2Vertex2qG = 0d0

! Loop over massless quarks
      do fType = 3,4
      do fGen = 1,2

! ****************************************************
! ************ Loop external quark (fType,fGen), quark (fType,fGen), gluon vertex ***********************
! ****************************************************	
 	
      vcAR1 = 0d0
      vcAL1 = 0d0

      vcAR2 = 0d0
      vcAL2 = 0d0

      vcAR3 = 0d0
      vcAL3 = 0d0

      vcBL1 = 0d0
      vcBR1 = 0d0

      vcBL2 = 0d0
      vcBR2 = 0d0

      vcBL3 = 0d0
      vcBR3 = 0d0

      vcCL1 = 0d0
      vcCR1 = 0d0

      vcCL2 = 0d0
      vcCR2 = 0d0

      vcCL3 = 0d0
      vcCR3 = 0d0

      vcDL1 = 0d0
      vcDR1 = 0d0

! Triangle 1: Quark-Quark-Gluon
      
      g0L = BL(fType,fType,fGen,fGen)
      g0R = BR(fType,fType,fGen,fGen)
      g1L = BL(fType,fType,fGen,fGen)
      g1R = BR(fType,fType,fGen,fGen)
      g2L = BL(fType,fType,fGen,fGen)
      g2R = BR(fType,fType,fGen,fGen)

	call myCget(0d0,sman,0d0,0d0,0d0,0d0)
	call myBget(sman,0d0,0d0)	
	
	colori = -1d0/(2*Nc)

      vcAL1 = vcAL1 + colori*(g0L*g1L*g2L*2*(xb0-2*xc00-sqrtS**2*(xc0+xc1+xc2)))
    
      vcAR1 = vcAR1 + colori*(g0R*g1R*g2R*2*(xb0-2*xc00-sqrtS**2*(xc0+xc1+xc2)))

      vcBL2 = vcBL2 + colori*(-2*(xc11+xc12+2*xc1+xc0+xc2)*g2R*g0R*g1R)

      vcBR2 = vcBR2 + colori*(-2*(xc11+xc12+2*xc1+xc0+xc2)*g2L*g0L*g1L)

      vcBL3 = vcBL3 + colori*(2*(-xc11+xc12+xc0+xc2)*g2R*g0R*g1R)

      vcBR3 = vcBR3 + colori*(2*(-xc11+xc12+xc0+xc2)*g2L*g0L*g1L)

      vcCL2 = vcCL2 + colori*(-2*(xc22+xc12+2*xc2+xc0+xc1)*g2L*g1L*g0L)

      vcCR2 = vcCR2 + colori*(-2*(xc22+xc12+2*xc2+xc0+xc1)*g2R*g1R*g0R)

      vcCL3 = vcCL3 + colori*(2*(-xc22+xc12+xc0+xc1)*g2L*g1L*g0L)

      vcCR3 = vcCR3 + colori*(2*(-xc22+xc12+xc0+xc1)*g2R*g1R*g0R)

      vcDL1 = vcDL1 + colori*(-2*(xc2+xc0+xc1)*g2R*g1R*g0R)

      vcDR1 = vcDR1 + colori*(-2*(xc2+xc0+xc1)*g2L*g1L*g0L)

! Triangle 2: Gluino-Squark-Squark

      do k=1,2

      g0 = -GS
      g1L = GluL(k,fType,fType,fGen,fGen)
      g1R = GluR(k,fType,fType,fGen,fGen)
      g2L = GluR(k,fType,fType,fGen,fGen)
      g2R = GluL(k,fType,fType,fGen,fGen)

      M0l = MGl
      M1l = MSf(k,fType,fGen)
      M2l = MSf(k,fType,fGen)
	
      call myCget(0d0,sman,0d0,M0l**2,M1l**2,M2l**2)
      call myBget(sman,M1l**2,M2l**2)	

      colori = -1d0/(2*Nc)

      vcAL1 = vcAL1 + colori*(g0*g1L*g2R*2*(xc00))

      vcAR1 = vcAR1 + colori*g0*g1R*g2L*2*(xc00)

      vcAL2 = vcAL2 + colori*g0*(M0l*(xc2+xc0+xc1)*g1L*g2R)

      vcAR2 = vcAR2 + colori*g0*(M0l*(xc2+xc0+xc1)*g1R*g2L)

      vcAL3 = vcAL3 + colori*g0*(g1L*g2L*M0l*(xc2-xc1))

      vcAR3 = vcAR3 + colori*g0*(g1L*g2L*M0l*(xc2-xc1))

      vcBL2 = vcBL2 + colori*g0*g1R*g2L*(xc1+xc11+xc12)

      vcBR2 = vcBR2 + colori*g0*g1L*g2R*(xc1+xc11+xc12)

      vcBL3 = vcBL3 + colori*g0*g1L*g2R*(xc12-xc11)

      vcBR3 = vcBR3 + colori*g0*g1R*g2L*(xc12-xc11)

      vcCL2 = vcCL2 + colori*g0*g1R*g2L*(xc12-xc22)

      vcCR2 = vcCR2 + colori*g0*g1L*g2R*(xc12-xc22)

      vcCL3 = vcCL3 - colori*g0*g1L*g2R*(xc22+xc12+xc2)

      vcCR3 = vcCR3 - colori*g0*g1R*g2L*(xc22+xc12+xc2)

      enddo
	  	  
! Triangle 3: Quark-Gluon-Gluon

      g0 = GS
      g1L = BL(fType,fType,fGen,fGen)
      g1R = BR(fType,fType,fGen,fGen)
      g2L = BL(fType,fType,fGen,fGen)
      g2R = BR(fType,fType,fGen,fGen)

      call myCget(0d0,sman,0d0,0d0,0d0,0d0)
      call myBget(sman,0d0,0d0)	

      colori = Nc/2d0

      vcAR1 = vcAR1 - colori*g0*(g1R*g2L*(2*(xb0+2*xc00)-2*sqrtS**2*(xc1+xc2)))

      vcAL1 = vcAL1 - colori*g0*(g1L*g2R*(2*(xb0+2*xc00)-2*sqrtS**2*(xc1+xc2)))

      vcBL2 = vcBL2 - colori*g0*(g1R*g2R*(-2*(xc11+xc12+xc1)+2*xc1+xc2))

      vcBR2 = vcBR2 - colori*g0*(g1L*g2L*(-2*(xc11+xc12+xc1)+2*xc1+xc2))

      vcBL3 = vcBL3 - colori*g0*(g1R*g2R*(-2*(xc12-xc11)+2*xc1+xc2))

      vcBR3 = vcBR3 - colori*g0*(g1L*g2L*(-2*(xc12-xc11)+2*xc1+xc2))

      vcCL2 = vcCL2 - colori*g0*g1L*g2L*(2*xc2+xc1+2*(xc12+xc22+xc2))

      vcCR2 = vcCR2 - colori*g0*g1R*g2R*(2*xc2+xc1+2*(xc12+xc22+xc2))

      vcCL3 = vcCL3 - colori*g0*g1L*g2L*(-2*xc2-xc1-2*(xc12-xc22))

      vcCR3 = vcCR3 - colori*g0*g1R*g2R*(-2*xc2-xc1-2*(xc12-xc22))

      vcDL1 = vcDL1 - colori*3*g0*g1R*g2R*(xc2+xc1)

      vcDR1 = vcDR1 - colori*3*g0*g1L*g2L*(xc2+xc1)

! Triangle 4: Squark-Gluino-Gluino

      do k=1,2

      g0 = GS
      g1L = GluL(k,fType,fType,fGen,fGen)
      g1R = GluR(k,fType,fType,fGen,fGen)
      g2L = GluR(k,fType,fType,fGen,fGen)
      g2R = GluL(k,fType,fType,fGen,fGen)

      M0l = MSf(k,fType,fGen)
      M1l = MGl
      M2l = MGl

      call myCget(0d0,sman,0d0,M0l**2,M1l**2,M2l**2)
      call myBget(sman,M1l**2,M2l**2)	
	
      colori = -Nc/2d0

      vcAL1 = vcAL1 + colori*g0*(g1L*g2R*(xb0-2*xc00+M0l**2*xc0-M1l*M2l*xc0))

      vcAR1 = vcAR1 + colori*g0*(g1R*g2L*(xb0-2*xc00+M0l**2*xc0-M1l*M2l*xc0))

      vcAL2 = vcAL2 - colori*g0*(-g1L*g2L*(M2l*xc2-M1l*xc1))

      vcAR2 = vcAR2 - colori*g0*(-g1R*g2R*(M2l*xc2-M1l*xc1))

      vcAL3 = vcAL3 - colori*g0*(-g1L*g2L*(M2l*xc2+M1l*xc1))

      vcAR3 = vcAR3 - colori*g0*(-g1R*g2R*(M2l*xc2+M1l*xc1)) 

      vcBL1 = vcBL1 - colori*g0*(g1R*g2R*(M1l*xc1-M2l*(xc1+xc0)))

      vcBR1 = vcBR1 - colori*g0*(g1L*g2L*(M1l*xc1-M2l*(xc1+xc0)))

      vcBL2 = vcBL2 + colori*g0*g1R*g2L*(xc12+xc11+xc1)

      vcBR2 = vcBR2 + colori*g0*g1L*g2R*(xc12+xc11+xc1)

      vcBL3 = vcBL3 - colori*g0*g1R*g2L*(-xc12+xc11+xc1)

      vcBR3 = vcBR3 - colori*g0*g1L*g2R*(-xc12+xc11+xc1)

      vcCL1 = vcCL1 - colori*g0*(g1L*g2L*(-M2l*xc2+M1l*(xc2+xc0)))

      vcCR1 = vcCR1 - colori*g0*(g1R*g2R*(-M2l*xc2+M1l*(xc2+xc0)))

      vcCL2 = vcCL2 - colori*g0*g1L*g2R*(xc12+xc22+xc2) 

      vcCR2 = vcCR2 - colori*g0*g1R*g2L*(xc12+xc22+xc2)

      vcCL3 = vcCL3 - colori*g0*g1L*g2R*(-xc12+xc22+xc2)

      vcCR3 = vcCR3 - colori*g0*g1R*g2L*(-xc12+xc22+xc2)

      vcDL1 = vcDL1 + colori*g0*g1R*g2L*(xc1+xc0+xc2)

      vcDR1 = vcDR1 + colori*g0*g1L*g2R*(xc1+xc0+xc2)

      enddo

! Take into account the factor 1/(4 Pi)^2

      vcAR1 = vcAR1/(4*Pi)**2
      vcAL1 = vcAL1/(4*Pi)**2

      vcAR2 = vcAR2/(4*Pi)**2
      vcAL2 = vcAL2/(4*Pi)**2

      vcAR3 = vcAR3/(4*Pi)**2
      vcAL3 = vcAL3/(4*Pi)**2

      vcBL1 = vcBL1/(4*Pi)**2
      vcBR1 = vcBR1/(4*Pi)**2

      vcBL2 = vcBL2/(4*Pi)**2
      vcBR2 = vcBR2/(4*Pi)**2

      vcBL3 = vcBL3/(4*Pi)**2
      vcBR3 = vcBR3/(4*Pi)**2

      vcCL1 = vcCL1/(4*Pi)**2
      vcCR1 = vcCL1/(4*Pi)**2

      vcCL2 = vcCL2/(4*Pi)**2
      vcCR2 = vcCR2/(4*Pi)**2

      vcCL3 = vcCL3/(4*Pi)**2
      vcCR3 = vcCL3/(4*Pi)**2

      vcDL1 = vcDL1/(4*Pi)**2
      vcDR1 = vcDR1/(4*Pi)**2

c---------------------Amplitudensquare for s-channel with Gluon-----------

      colorf=(Nc**2-1)/4d0

      stsT2QQbar_M2Vertex2qG =  stsT2QQbar_M2Vertex2qG - 2d0*colorf*dreal(((sqrtS**4*((-(muu4**2*vcAL1) + muu4**4*vcAL1 + 4*tred*vcAL1 - 
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
     -         BL(fType,fType,fGen,fGen) + 
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
     -                (muu1**2 - muu2**2)*vcAR3)))*BR(fType,fType,fGen,fGen))
     -       *dconjg(SfSfg(isf1,itt1,igen1,isf2,itt2,igen2))*
     -      SfSfg(isf1,itt1,igen1,isf2,itt2,igen2))/
     -    (dconjg(DenSV)*DenSV)))     


      enddo
      enddo

	end


! Results produced with FormCalc 9.9. These are in perfect agreement with the expressions for the quark-gluon vertex above.
!      double precision function FormCalc_QGV_Triangle1() ! Quark-Quark-Gluon
!      implicit none
!#include "../common/stsT2xx_Kinematics.h"
!#include "../common/stsT2xx_Model.h"
!#include "../../util/DMNLO_LoopIntegrals.h"
!
!      call myCget(0d0,sman,0d0,0d0,0d0,0d0)
!      call myBget(sman,0d0,0d0)
!
!      FormCalc_QGV_Triangle1 = 2*dble(-1d0/3d0*(gs**6*(m1**4- tman*uman)*(-xb0 + 2*xc00 + sman*(xc0 + xc1 + xc2)))/(Pi**2*sman**2))
!
!      end
!
!      double precision function FormCalc_QGV_Triangle2(fType,fGen) ! Gluino-Squark-Squark
!      implicit none
!#include "../common/stsT2xx_Kinematics.h"
!#include "../common/stsT2xx_Model.h"
!#include "../../util/DMNLO_LoopIntegrals.h"
!
!      integer fType,fGen,k
!
!      FormCalc_QGV_Triangle2 = 0d0
!
!      do k=1,2
!
!      call myCget(0d0,sman,0d0,Mgl2,MSf2(k,fType,fGen),MSf2(k,fType,fGen))
!
!      FormCalc_QGV_Triangle2 = FormCalc_QGV_Triangle2 + 2*dble((gs**6*(m1**4 - tman*uman)*xc00*(USf(k,1,fType,fGen)*USfC(k,1,fType,fGen) + USf(k,2,fType,fGen)*USfC(k,2,fType,fGen)))/(3*Pi**2*sman**2))
!
!      enddo
!
!      end
!
!      double precision function FormCalc_QGV_Triangle3() ! Quark-Gluon-Gluon
!      implicit none
!#include "../common/stsT2xx_Kinematics.h"
!#include "../common/stsT2xx_Model.h"
!#include "../../util/DMNLO_LoopIntegrals.h"
!
!      call myCget(sman,0d0,0d0,0d0,0d0,0d0)
!      call myBget(sman,0d0,0d0)
!
!      FormCalc_QGV_Triangle3 = 2*dble((-3*gs**6*(m1**4- tman*uman)*(xb0 + 2*xc00 + sman*(xc0 + xc2)))/(Pi**2*sman**2))
!
!      end
!
!      double precision function FormCalc_QGV_Triangle4(fType,fGen) ! Squark-Gluino-Gluino
!      implicit none
!#include "../common/stsT2xx_Kinematics.h"
!#include "../common/stsT2xx_Model.h"
!#include "../../util/DMNLO_LoopIntegrals.h"
!
!      integer fType,fGen,k
!
!      FormCalc_QGV_Triangle4 = 0d0
!
!      do k=1,2
!
!      call myCget(sman,0d0,0d0,Mgl2,Mgl2,MSf2(k,fType,fGen))
!      call myBget(sman,Mgl2,Mgl2)
!
!      FormCalc_QGV_Triangle4 = FormCalc_QGV_Triangle4 + 2*dble((-3*gs**6*(m1**4 - tman*uman)*(xb0 - Mgl**2*xc0 - 2*xc00 + xc0*MSf2(k,fType,fGen))*
!     &    (USf(k,1,fType,fGen)*USfC(k,1,fType,fGen) + USf(k,2,fType,fGen)*USfC(k,2,fType,fGen)))/(2*Pi**2*sman**2))
!
!      enddo
!
!      end