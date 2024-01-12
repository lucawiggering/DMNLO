      subroutine stst2gg_calc_A3g(PA_A3g,k1k1,k2k2,k3k3)
      implicit none 
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"
! The first, second and third index for the values 1,2,3 denote the partial amplitude A_{mno} connected to k_m^\rho k_n^\sigma k_o^\tau. 
! For example PA_A3g(1,2,3)-> k_1^\rho k_2^\sigma k_3^\tau.
! The fourth entry of the first index corresponds to the metric tensors. In this case, the second index corresponds to A_{g1}, A_{g2} and A_{g3}
! and the third index corresponds to the index of the momentum.
! For example PA_A3g(4,1,2)-> g^{\rho\sigma} k_2^\tau and PA_A3g(4,1,3)-> g^{\rho\tau} k_3^\sigma.
! Compare this statement to thesis for better understanding.

! All gluons are treated as off-shell . 
! Therefore, the amplitude still depends on the scalar products k1k1,k2k2,k3k3.

      double complex PA_A3g(4,3,3),tempPAkkk(3,3,3),tempPAg(3,3),c0
      integer k,gammaType,gammaGen
      double precision Pi2,GS3,color,k1k1,k2k2,k3k3,k1k3,M2

      c0=(0d0,0d0)
      Pi2=Pi**2
      GS3=GS**3
      k1k3=(k3k3+k1k1-k2k2)/2d0

      PA_A3g=0d0

!     Gluon loop k1
      call myBget(k1k1,0d0,0d0)

      color=Nc/2d0
      tempPAg=reshape((/(-9*GS3*xb0)/(32d0*Pi2),(9*GS3*xb0)/(32d0*Pi2),c0,c0,c0,c0,c0,c0,c0/),(/3,3/))

      PA_A3g(4,:,:) = PA_A3g(4,:,:) + color*tempPAg

!     Gluon loop k2
      call myBget(k2k2,0d0,0d0)

      color=Nc/2d0
      tempPAg=reshape((/c0,c0,c0,(9*GS3*xb0)/(32d0*Pi2),c0,(-9*GS3*xb0)/(32d0*Pi2),c0,c0,c0/),(/3,3/))

      PA_A3g(4,:,:) = PA_A3g(4,:,:) + color*tempPAg

!     Gluon loop k3
      call myBget(k3k3,0d0,0d0)

      color=Nc/2d0
      tempPAg=reshape((/c0,c0,c0,c0,c0,c0,c0,(9*GS3*xb0)/(32d0*Pi2),(-9*GS3*xb0)/(32d0*Pi2)/),(/3,3/))

      PA_A3g(4,:,:) = PA_A3g(4,:,:) + color*tempPAg
     
!     Quark triangles (both directions)
      do gammaType=3,4
      do gammaGen=1,3

! The first two generations are considered massless to avoid numerical instabilities arising from large logarithms.
      M2=Mf2(gammaType,gammaGen)
      if (gammaGen.le.2) then
      M2=0d0
      endif

      call myCget(k1k1,k2k2,k3k3,M2,M2,M2)
      call myBget(k2k2,M2,M2) 
      
      color=1d0/4d0
      tempPAg=reshape((/(GS3*(2*xb0+xb1-2*xc00+4*xc001-2*k1k1*xc11-2*k1k3*xc12+k3k3*xc2))/(2d0*Pi2),
     &(GS3*(xb1+2*xc00+4*xc001-k3k3*xc2))/(2d0*Pi2),
     &(GS3*(xb1+2*xc00+4*xc001-2*k1k3*(xc1+xc11)-k3k3*(2*xc12+xc2)))/(2d0*Pi2),c0,c0,c0,
     &-(GS3*(xb0+xb1-2*xc00-4*xc002+k1k1*(xc1+2*xc12)+2*k1k3*(xc2+xc22)))/(2d0*Pi2),
     &-(GS3*(xb0+xb1-2*xc00-4*xc002+k1k1*xc1))/(2d0*Pi2),
     &(GS3*(xb0-xb1-2*xc00+4*xc002+k1k1*xc1-2*k1k3*xc12-2*k3k3*xc22))/(2d0*Pi2)/),(/3,3/))
      tempPAkkk=reshape((/(2*GS3*(xc11+xc111))/Pi2,c0,
     &(GS3*(2*xc112+xc12))/Pi2,c0,c0,c0,
     &(GS3*(xc1+xc11+2*xc112+xc12))/Pi2,c0,
     &(GS3*(xc12+2*xc122))/Pi2,c0,c0,c0,c0,c0,c0,c0,c0,c0,
     &(GS3*(xc1+xc11+2*(xc112+xc12)))/Pi2,c0,
     &(GS3*(xc12+2*xc122+xc2+xc22))/Pi2,c0,c0,c0,
     &(GS3*(2*xc12+2*xc122+xc2+xc22))/Pi2,c0,(2*GS3*(xc22+xc222))/Pi2/),(/3,3,3/))

      PA_A3g(4,:,:) = PA_A3g(4,:,:) + color*tempPAg
      PA_A3g(1:3,:,:) = PA_A3g(1:3,:,:) + color*tempPAkkk

      enddo
      enddo

!     Gluino triangle
! The coefficiants A_{gi,m} and A{mno} correspond to one half of the quark triangle's coefficiants.
      call myCget(k1k1,k2k2,k3k3,MGl2,MGl2,MGl2)
      call myBget(k2k2,MGl2,MGl2) 

      color=Nc/2d0
      tempPAg=reshape((/(GS3*(2*xb0+xb1-2*xc00+4*xc001-2*k1k1*xc11-2*k1k3*xc12+k3k3*xc2))/(2d0*Pi2),
     &(GS3*(xb1+2*xc00+4*xc001-k3k3*xc2))/(2d0*Pi2),
     &(GS3*(xb1+2*xc00+4*xc001-2*k1k3*(xc1+xc11)-k3k3*(2*xc12+xc2)))/(2d0*Pi2),c0,c0,c0,
     &-(GS3*(xb0+xb1-2*xc00-4*xc002+k1k1*(xc1+2*xc12)+2*k1k3*(xc2+xc22)))/(2d0*Pi2),
     &-(GS3*(xb0+xb1-2*xc00-4*xc002+k1k1*xc1))/(2d0*Pi2),
     &(GS3*(xb0-xb1-2*xc00+4*xc002+k1k1*xc1-2*k1k3*xc12-2*k3k3*xc22))/(2d0*Pi2)/),(/3,3/))/2d0
      tempPAkkk=reshape((/(2*GS3*(xc11+xc111))/Pi2,c0,
     &(GS3*(2*xc112+xc12))/Pi2,c0,c0,c0,
     &(GS3*(xc1+xc11+2*xc112+xc12))/Pi2,c0,
     &(GS3*(xc12+2*xc122))/Pi2,c0,c0,c0,c0,c0,c0,c0,c0,c0,
     &(GS3*(xc1+xc11+2*(xc112+xc12)))/Pi2,c0,
     &(GS3*(xc12+2*xc122+xc2+xc22))/Pi2,c0,c0,c0,
     &(GS3*(2*xc12+2*xc122+xc2+xc22))/Pi2,c0,(2*GS3*(xc22+xc222))/Pi2/),(/3,3,3/))/2d0

      PA_A3g(4,:,:) = PA_A3g(4,:,:) + color*tempPAg
      PA_A3g(1:3,:,:) = PA_A3g(1:3,:,:) + color*tempPAkkk

!     Squark triangle
      do gammaType=3,4
      do gammaGen=1,3
      do k=1,2

      M2=MSf2(k,gammaType,gammaGen)
      call myCget(k1k1,k2k2,k3k3,M2,M2,M2)

      color=1d0/4d0
      tempPAg=reshape((/-((GS3*xc001)/Pi2),-(GS3*(xc00+2*xc001))/(2d0*Pi2),
     &-(GS3*(xc00+2*xc001))/(2d0*Pi2),c0,c0,c0,-(GS3*(xc00+2*xc002))/(2d0*Pi2),
     &-(GS3*(xc00+2*xc002))/(2d0*Pi2),-((GS3*xc002)/Pi2)/),(/3,3/))
      tempPAkkk=reshape((/-(GS3*(xc1+4*(xc11+xc111)))/(4d0*Pi2),c0,
     &-(GS3*(2*xc112+xc12))/(2d0*Pi2),c0,c0,c0,
     &-(GS3*(xc1+2*(xc11+2*xc112+xc12)))/(4d0*Pi2),c0,
     &-(GS3*(xc12+2*xc122))/(2d0*Pi2),c0,c0,c0,c0,c0,c0,c0,c0,c0,
     &-(GS3*(xc0+4*xc1+4*xc11+8*xc112+8*xc12+2*xc2))/(8d0*Pi2),c0,
     &-(GS3*(2*xc12+4*xc122+xc2+2*xc22))/(4d0*Pi2),c0,c0,c0,
     &-(GS3*(xc0+2*xc1+4*(2*xc12+2*xc122+xc2+xc22)))/(8d0*Pi2),c0,
     &-(GS3*(xc2+4*(xc22+xc222)))/(4d0*Pi2)/),(/3,3,3/))
      
      PA_A3g(4,:,:) = PA_A3g(4,:,:) + color*tempPAg
      PA_A3g(1:3,:,:) = PA_A3g(1:3,:,:) + color*tempPAkkk

      enddo
      enddo
      enddo

!     Gluon triangle
      call myCget(k1k1,k2k2,k3k3,0d0,0d0,0d0)
      call myBget(k2k2,0d0,0d0) 

      color=Nc/2d0
      tempPAg=reshape((/(GS3*(6*xb0+2*xb1+2*k3k3*xc0-18*xc001-2*k1k1*xc1+k3k3*xc1-2*k3k3*xc12+2*k1k3*(2*xc0+5*xc1-xc11+4*xc2)))/(16d0*Pi2),
     &-(GS3*(xb0-2*xb1+k3k3*xc0+11*xc00+18*xc001-k3k3*xc1+2*k1k1*(xc1+xc11)+2*k3k3*xc12-k3k3*xc2+k1k3*(3*xc0+13*xc1+2*xc11+2*xc12+3*xc2)))/(16d0*Pi2),
     &-(GS3*(xb0-2*xb1+5*k3k3*xc0+11*xc00+18*xc001+10*k3k3*xc1+2*k1k1*(xc1+xc11)+4*k3k3*xc2+k1k3*(-6*xc1+2*xc12+3*xc2)))/(16d0*Pi2),c0,c0,c0,
     &-(GS3*(3*xb0+2*xb1+11*xc00+18*xc002+k1k3*(3*xc1+2*xc12-6*xc2)+2*k3k3*xc2+k1k1*(5*xc0+4*xc1+10*xc2)+2*k3k3*xc22))/(16d0*Pi2),
     &-(GS3*(3*xb0+2*xb1+11*xc00+18*xc002+k1k1*(xc0-xc1+2*xc12-xc2)+2*k3k3*xc2+2*k3k3*xc22+k1k3*(3*xc0+3*xc1+2*xc12+13*xc2+2*xc22)))/(16d0*Pi2),
     &(GS3*(4*xb0-2*xb1-18*xc002-2*k3k3*xc2+k1k1*(2*xc0-2*xc12+xc2)+2*k1k3*(2*xc0+4*xc1+5*xc2-xc22)))/(16d0*Pi2)/),(/3,3/))
      tempPAkkk=reshape((/(-9*GS3*(xc11+xc111))/(8d0*Pi2),c0,
     &(GS3*(3*xc1-18*xc112-9*xc12-3*xc2))/(16d0*Pi2),c0,c0,c0,
     &-(GS3*(3*xc0+15*xc1+9*xc11+18*xc112+9*xc12+3*xc2))/(16d0*Pi2),c0,
     &-(GS3*(3*xc1+9*xc12+18*xc122-3*xc2))/(16d0*Pi2),c0,c0,c0,c0,c0,c0,c0,c0,c0,
     &-(GS3*(-3*xc0+6*xc1+9*xc11+18*xc112+18*xc12-6*xc2))/(16d0*Pi2),c0,
     &-(GS3*(3*xc0+3*xc1+9*xc12+18*xc122+15*xc2+9*xc22))/(16d0*Pi2),c0,c0,c0,
     &-(GS3*(-3*xc0-6*xc1+18*xc12+18*xc122+6*xc2+9*xc22))/(16d0*Pi2),c0,(-9*GS3*(xc22+xc222))/(8d0*Pi2)/),(/3,3,3/))

      PA_A3g(4,:,:) = PA_A3g(4,:,:) + color*tempPAg
      PA_A3g(1:3,:,:) = PA_A3g(1:3,:,:) + color*tempPAkkk

!     Ghost triangles
      call myCget(k1k1,k2k2,k3k3,0d0,0d0,0d0)

      color=Nc/2d0
      tempPAg=reshape((/(GS3*xc001)/(8d0*Pi2),(GS3*(xc00+2*xc001))/(16d0*Pi2),
     &(GS3*(xc00+2*xc001))/(16d0*Pi2),c0,c0,c0,(GS3*(xc00+2*xc002))/(16d0*Pi2),
     &(GS3*(xc00+2*xc002))/(16d0*Pi2),(GS3*xc002)/(8d0*Pi2)/),(/3,3/))
      tempPAkkk=reshape((/(GS3*(xc11+xc111))/(8d0*Pi2),c0,(GS3*(2*xc112+xc12))/(16d0*Pi2),c0,c0,c0,
     &(GS3*(xc1+xc11+2*xc112+xc12))/(16d0*Pi2),c0,(GS3*(xc12+2*xc122))/(16d0*Pi2),c0,c0,c0,c0,c0,c0,c0,c0,c0,
     &(GS3*(xc1+xc11+2*(xc112+xc12)))/(16d0*Pi2),c0,(GS3*(xc12+2*xc122+xc2+xc22))/(16d0*Pi2),c0,c0,c0,
     &(GS3*(2*xc12+2*xc122+xc2+xc22))/(16d0*Pi2),c0,(GS3*(xc22+xc222))/(8d0*Pi2)/),(/3,3,3/))

      PA_A3g(4,:,:) = PA_A3g(4,:,:) + color*tempPAg
      PA_A3g(1:3,:,:) = PA_A3g(1:3,:,:) + color*tempPAkkk

      end


      subroutine stst2gg_calc_M2Vertex3g(stst2gg_M2Vertex3g)
	     implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"

	
      double complex PA_A3g(4,3,3)
      double complex DenS, DenSft(2,4,3), DenSfu(2,4,3)
      integer i,j,k,m,n,l
      double precision k1k1,k2k2,k3k3,k1p1,k1p2,M2,Krondelta,stst2gg_M2Vertex3g(3)

      M2=MSf2(iSf1,itt1,igen1)
!Propagator
      DenS = sqrtS**2 

      do i=1,2
      do j=1,4
      do k=1,3
      DenSft(i,j,k)=(tman - MSf2(i,j,k)) !+ cI*MSf(i,j,k)*WSf(i,j,k)
      DenSfu(i,j,k)=(uman - MSf2(i,j,k)) !+ cI*MSf(i,j,k)*WSf(i,j,k)
      enddo
      enddo
      enddo

! Values for scalar products for the 3g-vertex correction (k1-> -k2, k2-> p1+p2, k3-> k1)
      k1k1=0d0
      k2k2=sman
      k3k3=0d0

      call stst2gg_calc_A3g(PA_A3g,k1k1,k2k2,k3k3)
      PA_A3g=dconjg(PA_A3g)

! Values for scalar products for the 2->2 process
      k1k2=sman/2d0
      p1p2=(sman-2*M2)/2d0
      k1p1=(M2-tman)/2d0
      k1p2=(M2-uman)/2d0

      stst2gg_M2Vertex3g = 0d0

! only antisymmetric octet contributes
      stst2gg_M2Vertex3g(3) = 2*dble((CF*CV*gs**2*Nc*dconjg(SfSfg(isf1,itt1,igen1,isf2,itt2,igen2))*
     &    ((-1 + FPSk1)*(-1 + FPSk2)*k1k2**4*(k1p1 - k1p2)**2*sman*(PA_A3g(3,1,1) + PA_A3g(3,3,1)) - 
     &      4*FPSk1*FPSk2*(k1p1 - k1p2)**2*(k1p1 + k1p2)**3*sman*
     &       (k1p1*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + k1p2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) - PA_A3g(4,1,2) + PA_A3g(4,3,2)) + 
     &      k1k2**3*((-1 + FPSk1)*(-1 + FPSk2)*k1p1**3*(sman*(PA_A3g(2,1,1) + PA_A3g(2,3,1) - 4*PA_A3g(3,1,1) - PA_A3g(3,1,2) - 4*PA_A3g(3,3,1) - PA_A3g(3,3,2)) + 
     &            4*k1p2*(PA_A3g(1,1,1) + 2*PA_A3g(1,1,3) + PA_A3g(1,3,1) + 2*PA_A3g(1,3,3) - 4*PA_A3g(3,1,1) + PA_A3g(3,1,3) - 4*PA_A3g(3,3,1) + PA_A3g(3,3,3))) + 
     &         k1p2*sman*((-1 + FPSk1)*(-1 + FPSk2)*k1p2**2*(PA_A3g(2,1,1) + PA_A3g(2,3,1) - 4*PA_A3g(3,1,1) - PA_A3g(3,1,2) - 4*PA_A3g(3,3,1) - PA_A3g(3,3,2)) + 
     &            2*(M2 - p1p2)*((-1 + FPSk1)*PA_A3g(4,1,1) + (-1 + FPSk2)*PA_A3g(4,3,3)) + 
     &            k1p2*(2*M2*((-1 + FPSk1)*PA_A3g(2,1,1) + (-1 + FPSk1)*PA_A3g(2,3,1) - (-1 + FPSk2)*(PA_A3g(3,1,2) + PA_A3g(3,3,2))) + 
     &               2*p1p2*((-1 + FPSk1)*PA_A3g(2,1,1) + (-1 + FPSk1)*PA_A3g(2,3,1) - (-1 + FPSk2)*(PA_A3g(3,1,2) + PA_A3g(3,3,2))) - 
     &               (-1 + FPSk1)*(-1 + FPSk2)*(PA_A3g(4,1,1) + PA_A3g(4,2,1) + PA_A3g(4,2,3) + PA_A3g(4,3,3)))) + 
     &         k1p1**2*(-((-1 + FPSk1)*(-1 + FPSk2)*k1p2*sman*(PA_A3g(2,1,1) + PA_A3g(2,3,1) - 4*PA_A3g(3,1,1) - PA_A3g(3,1,2) - 4*PA_A3g(3,3,1) - PA_A3g(3,3,2))) - 
     &            8*(-1 + FPSk1)*(-1 + FPSk2)*k1p2**2*(PA_A3g(1,1,1) + 2*PA_A3g(1,1,3) + PA_A3g(1,3,1) + 2*PA_A3g(1,3,3) - 4*PA_A3g(3,1,1) + PA_A3g(3,1,3) - 
     &               4*PA_A3g(3,3,1) + PA_A3g(3,3,3)) + sman*(2*M2*
     &                ((-1 + FPSk1)*PA_A3g(2,1,1) + (-1 + FPSk1)*PA_A3g(2,3,1) - (-1 + FPSk2)*(PA_A3g(3,1,2) + PA_A3g(3,3,2))) + 
     &               2*p1p2*((-1 + FPSk1)*PA_A3g(2,1,1) + (-1 + FPSk1)*PA_A3g(2,3,1) - (-1 + FPSk2)*(PA_A3g(3,1,2) + PA_A3g(3,3,2))) - 
     &               (-1 + FPSk1)*(-1 + FPSk2)*(PA_A3g(4,1,1) + PA_A3g(4,2,1) + PA_A3g(4,2,3) + PA_A3g(4,3,3)))) + 
     &         k1p1*(-((-1 + FPSk1)*(-1 + FPSk2)*k1p2**2*sman*(PA_A3g(2,1,1) + PA_A3g(2,3,1) - 4*PA_A3g(3,1,1) - PA_A3g(3,1,2) - 4*PA_A3g(3,3,1) - PA_A3g(3,3,2))) + 
     &            4*(-1 + FPSk1)*(-1 + FPSk2)*k1p2**3*(PA_A3g(1,1,1) + 2*PA_A3g(1,1,3) + PA_A3g(1,3,1) + 2*PA_A3g(1,3,3) - 4*PA_A3g(3,1,1) + PA_A3g(3,1,3) - 
     &               4*PA_A3g(3,3,1) + PA_A3g(3,3,3)) + 2*(M2 - p1p2)*sman*((-1 + FPSk1)*PA_A3g(4,1,1) + (-1 + FPSk2)*PA_A3g(4,3,3)) - 
     &            2*k1p2*(2*p1p2*sman*((-1 + FPSk1)*PA_A3g(2,1,1) + (-1 + FPSk1)*PA_A3g(2,3,1) - (-1 + FPSk2)*(PA_A3g(3,1,2) + PA_A3g(3,3,2))) + 
     &               4*p1p2*((-1 + FPSk1)*PA_A3g(4,1,1) - 2*(-1 + FPSk1)*PA_A3g(4,1,3) - (-1 + FPSk2)*(2*PA_A3g(4,3,1) - PA_A3g(4,3,3))) + 
     &               2*M2*(sman*((-1 + FPSk1)*PA_A3g(2,1,1) + (-1 + FPSk1)*PA_A3g(2,3,1) - (-1 + FPSk2)*(PA_A3g(3,1,2) + PA_A3g(3,3,2))) - 
     &                  2*(-1 + FPSk1)*PA_A3g(4,1,1) + 4*(-1 + FPSk1)*PA_A3g(4,1,3) + 2*(-1 + FPSk2)*(2*PA_A3g(4,3,1) - PA_A3g(4,3,3))) - 
     &               (-1 + FPSk1)*(-1 + FPSk2)*sman*(PA_A3g(4,1,1) + PA_A3g(4,2,1) + PA_A3g(4,2,3) + PA_A3g(4,3,3))))) + 
     &      2*k1k2*(-(FPSk1*(k1p1**5*(8*k1p2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + sman*(PA_A3g(2,1,2) + PA_A3g(2,3,2) - 2*(PA_A3g(3,1,2) + PA_A3g(3,3,2)))) - 
     &              2*k1p1**3*(8*k1p2**3*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 
     &                 k1p2**2*(sman*(PA_A3g(2,1,2) + PA_A3g(2,3,2) - 2*(PA_A3g(3,1,2) + PA_A3g(3,3,2))) + 6*PA_A3g(4,3,2)) + 
     &                 sman*(p1p2*(-PA_A3g(4,1,2) + PA_A3g(4,3,2)) + M2*(PA_A3g(4,1,2) + PA_A3g(4,3,2))) - 
     &                 2*k1p2*sman*(PA_A3g(4,2,1) + PA_A3g(4,2,3) + PA_A3g(4,3,1) + PA_A3g(4,3,3))) + 
     &              k1p1*k1p2**2*(8*k1p2**3*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 
     &                 k1p2**2*(sman*(PA_A3g(2,1,2) + PA_A3g(2,3,2) - 2*(PA_A3g(3,1,2) + PA_A3g(3,3,2))) + 12*PA_A3g(4,3,2)) + 
     &                 2*sman*(M2*(-3*PA_A3g(4,1,2) + PA_A3g(4,3,2)) + p1p2*(3*PA_A3g(4,1,2) + PA_A3g(4,3,2))) + 
     &                 4*k1p2*sman*(PA_A3g(4,2,1) + PA_A3g(4,2,3) + PA_A3g(4,3,1) + PA_A3g(4,3,3))) + 
     &              k1p2**3*sman*(k1p2**2*(PA_A3g(2,1,2) + PA_A3g(2,3,2) - 2*(PA_A3g(3,1,2) + PA_A3g(3,3,2))) - 
     &                 2*(p1p2*(-PA_A3g(4,1,2) + PA_A3g(4,3,2)) + M2*(PA_A3g(4,1,2) + PA_A3g(4,3,2))) - 
     &                 k1p2*(2*M2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 2*p1p2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + PA_A3g(4,1,2) + 2*PA_A3g(4,2,1) + 2*PA_A3g(4,2,3) - 
     &                    PA_A3g(4,3,2) + 2*PA_A3g(4,3,3))) + k1p1**4*
     &               (k1p2*(sman*(PA_A3g(2,1,2) + PA_A3g(2,3,2) - 2*(PA_A3g(3,1,2) + PA_A3g(3,3,2))) + 12*PA_A3g(4,3,2)) - 
     &                 sman*(2*M2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 2*p1p2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + PA_A3g(4,1,2) + 2*PA_A3g(4,2,1) + 2*PA_A3g(4,2,3) - 
     &                    PA_A3g(4,3,2) + 2*PA_A3g(4,3,3))) - 2*k1p1**2*k1p2*
     &               (k1p2**2*(sman*(PA_A3g(2,1,2) + PA_A3g(2,3,2) - 2*(PA_A3g(3,1,2) + PA_A3g(3,3,2))) + 6*PA_A3g(4,3,2)) + 
     &                 sman*(M2*(3*PA_A3g(4,1,2) - PA_A3g(4,3,2)) - p1p2*(3*PA_A3g(4,1,2) + PA_A3g(4,3,2))) + 
     &                 k1p2*sman*(-2*M2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) - 2*p1p2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) - PA_A3g(4,1,2) + 2*PA_A3g(4,2,1) + 2*PA_A3g(4,2,3) + 
     &                    4*PA_A3g(4,3,1) + PA_A3g(4,3,2) + 2*PA_A3g(4,3,3))))) + 
     &         FPSk2*(k1p1**5*(8*(-1 + FPSk1)*k1p2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 
     &               sman*(2*(-1 + FPSk1)*PA_A3g(2,1,1) + (-1 + 2*FPSk1)*PA_A3g(2,1,2) - 2*PA_A3g(2,3,1) + 2*FPSk1*PA_A3g(2,3,1) - PA_A3g(2,3,2) + 
     &                  2*FPSk1*PA_A3g(2,3,2) - 2*FPSk1*PA_A3g(3,1,2) - 2*FPSk1*PA_A3g(3,3,2))) + 
     &            k1p2**3*sman*(k1p2**2*(2*(-1 + FPSk1)*PA_A3g(2,1,1) + (-1 + 2*FPSk1)*PA_A3g(2,1,2) - 2*PA_A3g(2,3,1) + 2*FPSk1*PA_A3g(2,3,1) - PA_A3g(2,3,2) + 
     &                  2*FPSk1*PA_A3g(2,3,2) - 2*FPSk1*PA_A3g(3,1,2) - 2*FPSk1*PA_A3g(3,3,2)) - 
     &               2*(p1p2*(PA_A3g(4,1,2) - PA_A3g(4,3,2)) + M2*(PA_A3g(4,1,2) + PA_A3g(4,3,2))) + 
     &               k1p2*(2*M2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 2*p1p2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 2*PA_A3g(4,1,1) - 2*FPSk1*PA_A3g(4,1,1) + PA_A3g(4,1,2) - 
     &                  2*FPSk1*PA_A3g(4,1,2) + 2*PA_A3g(4,2,1) - 2*FPSk1*PA_A3g(4,2,1) + 2*PA_A3g(4,2,3) - 2*FPSk1*PA_A3g(4,2,3) - PA_A3g(4,3,2) + 
     &                  2*FPSk1*PA_A3g(4,3,2) - 2*FPSk1*PA_A3g(4,3,3))) + 
     &            k1p1**4*(k1p2*(sman*(2*(-1 + FPSk1)*PA_A3g(2,1,1) + (-1 + 2*FPSk1)*PA_A3g(2,1,2) - 2*PA_A3g(2,3,1) + 2*FPSk1*PA_A3g(2,3,1) - PA_A3g(2,3,2) + 
     &                     2*FPSk1*PA_A3g(2,3,2) - 2*FPSk1*PA_A3g(3,1,2) - 2*FPSk1*PA_A3g(3,3,2)) + 12*(-((-1 + FPSk1)*PA_A3g(4,1,2)) + FPSk1*PA_A3g(4,3,2))) + 
     &               sman*(2*M2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 2*p1p2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 2*PA_A3g(4,1,1) - 2*FPSk1*PA_A3g(4,1,1) + PA_A3g(4,1,2) - 
     &                  2*FPSk1*PA_A3g(4,1,2) + 2*PA_A3g(4,2,1) - 2*FPSk1*PA_A3g(4,2,1) + 2*PA_A3g(4,2,3) - 2*FPSk1*PA_A3g(4,2,3) - PA_A3g(4,3,2) + 
     &                  2*FPSk1*PA_A3g(4,3,2) - 2*FPSk1*PA_A3g(4,3,3))) - 
     &            2*k1p1**3*(8*(-1 + FPSk1)*k1p2**3*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 
     &               k1p2**2*(sman*(2*(-1 + FPSk1)*PA_A3g(2,1,1) + (-1 + 2*FPSk1)*PA_A3g(2,1,2) - 2*PA_A3g(2,3,1) + 2*FPSk1*PA_A3g(2,3,1) - PA_A3g(2,3,2) + 
     &                     2*FPSk1*PA_A3g(2,3,2) - 2*FPSk1*PA_A3g(3,1,2) - 2*FPSk1*PA_A3g(3,3,2)) - 6*(-1 + FPSk1)*PA_A3g(4,1,2) + 6*FPSk1*PA_A3g(4,3,2)) + 
     &               sman*(p1p2*(PA_A3g(4,1,2) - PA_A3g(4,3,2)) + M2*(PA_A3g(4,1,2) + PA_A3g(4,3,2))) - 
     &               2*k1p2*sman*((-1 + FPSk1)*PA_A3g(4,1,1) + (-1 + FPSk1)*PA_A3g(4,1,3) - PA_A3g(4,2,1) + FPSk1*PA_A3g(4,2,1) - PA_A3g(4,2,3) + 
     &                  FPSk1*PA_A3g(4,2,3) + FPSk1*PA_A3g(4,3,1) + FPSk1*PA_A3g(4,3,3))) + 
     &            k1p1*k1p2**2*(8*(-1 + FPSk1)*k1p2**3*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 
     &               2*sman*(M2*(PA_A3g(4,1,2) - 3*PA_A3g(4,3,2)) + p1p2*(PA_A3g(4,1,2) + 3*PA_A3g(4,3,2))) + 
     &               k1p2**2*(sman*(2*(-1 + FPSk1)*PA_A3g(2,1,1) + (-1 + 2*FPSk1)*PA_A3g(2,1,2) - 2*PA_A3g(2,3,1) + 2*FPSk1*PA_A3g(2,3,1) - PA_A3g(2,3,2) + 
     &                     2*FPSk1*PA_A3g(2,3,2) - 2*FPSk1*PA_A3g(3,1,2) - 2*FPSk1*PA_A3g(3,3,2)) + 12*(-((-1 + FPSk1)*PA_A3g(4,1,2)) + FPSk1*PA_A3g(4,3,2))) + 
     &               4*k1p2*sman*((-1 + FPSk1)*PA_A3g(4,1,1) + (-1 + FPSk1)*PA_A3g(4,1,3) - PA_A3g(4,2,1) + FPSk1*PA_A3g(4,2,1) - PA_A3g(4,2,3) + 
     &                  FPSk1*PA_A3g(4,2,3) + FPSk1*PA_A3g(4,3,1) + FPSk1*PA_A3g(4,3,3))) - 
     &            2*k1p1**2*k1p2*(k1p2**2*(sman*(2*(-1 + FPSk1)*PA_A3g(2,1,1) + (-1 + 2*FPSk1)*PA_A3g(2,1,2) - 2*PA_A3g(2,3,1) + 2*FPSk1*PA_A3g(2,3,1) - 
     &                     PA_A3g(2,3,2) + 2*FPSk1*PA_A3g(2,3,2) - 2*FPSk1*PA_A3g(3,1,2) - 2*FPSk1*PA_A3g(3,3,2)) - 6*(-1 + FPSk1)*PA_A3g(4,1,2) + 
     &                  6*FPSk1*PA_A3g(4,3,2)) - sman*(M2*(PA_A3g(4,1,2) - 3*PA_A3g(4,3,2)) + p1p2*(PA_A3g(4,1,2) + 3*PA_A3g(4,3,2))) + 
     &               k1p2*sman*(2*M2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 2*p1p2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) - 2*PA_A3g(4,1,1) + 2*FPSk1*PA_A3g(4,1,1) + 
     &                  PA_A3g(4,1,2) - 2*FPSk1*PA_A3g(4,1,2) - 4*PA_A3g(4,1,3) + 4*FPSk1*PA_A3g(4,1,3) - 2*PA_A3g(4,2,1) + 2*FPSk1*PA_A3g(4,2,1) - 2*PA_A3g(4,2,3) + 
     &                  2*FPSk1*PA_A3g(4,2,3) + 4*FPSk1*PA_A3g(4,3,1) - PA_A3g(4,3,2) + 2*FPSk1*PA_A3g(4,3,2) + 2*FPSk1*PA_A3g(4,3,3))))) + 
     &      k1k2**2*(k1p1**4*(4*(-1 + FPSk1)*(-1 + FPSk2)*k1p2*(PA_A3g(1,1,2) + PA_A3g(1,3,2) - 5*PA_A3g(2,1,1) - PA_A3g(2,1,3) - 5*PA_A3g(2,3,1) - PA_A3g(2,3,3) + 
     &               5*PA_A3g(3,1,2) + 5*PA_A3g(3,3,2)) - sman*(2*(-1 + FPSk1)*(-1 + 2*FPSk2)*PA_A3g(2,1,1) + (-1 + FPSk1)*(-1 + FPSk2)*PA_A3g(2,1,2) + 
     &               2*PA_A3g(2,3,1) - 2*FPSk1*PA_A3g(2,3,1) - 4*FPSk2*PA_A3g(2,3,1) + 4*FPSk1*FPSk2*PA_A3g(2,3,1) + PA_A3g(2,3,2) - FPSk1*PA_A3g(2,3,2) - 
     &               FPSk2*PA_A3g(2,3,2) + FPSk1*FPSk2*PA_A3g(2,3,2) - 4*PA_A3g(3,1,1) + 4*FPSk1*PA_A3g(3,1,1) + 4*FPSk2*PA_A3g(3,1,1) - 
     &               4*FPSk1*FPSk2*PA_A3g(3,1,1) - 2*PA_A3g(3,1,2) + 4*FPSk1*PA_A3g(3,1,2) + 2*FPSk2*PA_A3g(3,1,2) - 4*FPSk1*FPSk2*PA_A3g(3,1,2) - 4*PA_A3g(3,3,1) + 
     &               4*FPSk1*PA_A3g(3,3,1) + 4*FPSk2*PA_A3g(3,3,1) - 4*FPSk1*FPSk2*PA_A3g(3,3,1) - 2*PA_A3g(3,3,2) + 4*FPSk1*PA_A3g(3,3,2) + 2*FPSk2*PA_A3g(3,3,2) - 
     &               4*FPSk1*FPSk2*PA_A3g(3,3,2))) + k1p1**3*(-4*(-1 + FPSk1)*(-1 + FPSk2)*k1p2**2*
     &             (PA_A3g(1,1,2) + PA_A3g(1,3,2) - 5*PA_A3g(2,1,1) - PA_A3g(2,1,3) - 5*PA_A3g(2,3,1) - PA_A3g(2,3,3) + 5*PA_A3g(3,1,2) + 5*PA_A3g(3,3,2)) - 
     &            sman*(2*M2*(2*(-1 + FPSk1)*PA_A3g(2,1,1) + (-2 + FPSk1 + FPSk2)*PA_A3g(2,1,2) - 2*PA_A3g(2,3,1) + 2*FPSk1*PA_A3g(2,3,1) - 2*PA_A3g(2,3,2) + 
     &                  FPSk1*PA_A3g(2,3,2) + FPSk2*PA_A3g(2,3,2) + 2*PA_A3g(3,1,2) - 2*FPSk2*PA_A3g(3,1,2) + 2*PA_A3g(3,3,2) - 2*FPSk2*PA_A3g(3,3,2)) + 
     &               2*p1p2*(2*(-1 + FPSk1)*PA_A3g(2,1,1) + (-2 + FPSk1 + FPSk2)*PA_A3g(2,1,2) - 2*PA_A3g(2,3,1) + 2*FPSk1*PA_A3g(2,3,1) - 2*PA_A3g(2,3,2) + 
     &                  FPSk1*PA_A3g(2,3,2) + FPSk2*PA_A3g(2,3,2) + 2*PA_A3g(3,1,2) - 2*FPSk2*PA_A3g(3,1,2) + 2*PA_A3g(3,3,2) - 2*FPSk2*PA_A3g(3,3,2)) - 
     &               2*PA_A3g(4,1,1) + 2*FPSk1*PA_A3g(4,1,1) + 4*FPSk2*PA_A3g(4,1,1) - 4*FPSk1*FPSk2*PA_A3g(4,1,1) - PA_A3g(4,1,2) + FPSk1*PA_A3g(4,1,2) + 
     &               FPSk2*PA_A3g(4,1,2) - FPSk1*FPSk2*PA_A3g(4,1,2) - 4*PA_A3g(4,2,1) + 4*FPSk1*PA_A3g(4,2,1) + 4*FPSk2*PA_A3g(4,2,1) - 
     &               4*FPSk1*FPSk2*PA_A3g(4,2,1) - 4*PA_A3g(4,2,3) + 4*FPSk1*PA_A3g(4,2,3) + 4*FPSk2*PA_A3g(4,2,3) - 4*FPSk1*FPSk2*PA_A3g(4,2,3) + PA_A3g(4,3,2) - 
     &               FPSk1*PA_A3g(4,3,2) - FPSk2*PA_A3g(4,3,2) + FPSk1*FPSk2*PA_A3g(4,3,2) - 2*PA_A3g(4,3,3) + 4*FPSk1*PA_A3g(4,3,3) + 2*FPSk2*PA_A3g(4,3,3) - 
     &               4*FPSk1*FPSk2*PA_A3g(4,3,3)) - 4*k1p2*(-4*p1p2*PA_A3g(2,1,2) - 4*p1p2*PA_A3g(2,3,2) - 4*M2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 
     &               (-1 + FPSk1)*(-1 + FPSk2)*sman*(PA_A3g(1,1,1) + PA_A3g(1,1,3) + PA_A3g(1,3,1) + PA_A3g(1,3,3) + PA_A3g(3,1,1) + PA_A3g(3,1,3) + PA_A3g(3,3,1) + 
     &                  PA_A3g(3,3,3)) - PA_A3g(4,1,1) + FPSk1*PA_A3g(4,1,1) + 3*FPSk2*PA_A3g(4,1,1) - 3*FPSk1*FPSk2*PA_A3g(4,1,1) - PA_A3g(4,1,3) + 
     &               FPSk1*PA_A3g(4,1,3) - 3*FPSk2*PA_A3g(4,1,3) + 3*FPSk1*FPSk2*PA_A3g(4,1,3) - 6*PA_A3g(4,2,1) + 2*FPSk1*PA_A3g(4,2,1) + 2*FPSk2*PA_A3g(4,2,1) - 
     &               2*FPSk1*FPSk2*PA_A3g(4,2,1) - 6*PA_A3g(4,2,3) + 2*FPSk1*PA_A3g(4,2,3) + 2*FPSk2*PA_A3g(4,2,3) - 2*FPSk1*FPSk2*PA_A3g(4,2,3) - PA_A3g(4,3,1) - 
     &               3*FPSk1*PA_A3g(4,3,1) + FPSk2*PA_A3g(4,3,1) + 3*FPSk1*FPSk2*PA_A3g(4,3,1) - PA_A3g(4,3,3) + 3*FPSk1*PA_A3g(4,3,3) + FPSk2*PA_A3g(4,3,3) - 
     &               3*FPSk1*FPSk2*PA_A3g(4,3,3))) + k1p1*(4*(-1 + FPSk1)*(-1 + FPSk2)*k1p2**4*
     &             (PA_A3g(1,1,2) + PA_A3g(1,3,2) - 5*PA_A3g(2,1,1) - PA_A3g(2,1,3) - 5*PA_A3g(2,3,1) - PA_A3g(2,3,3) + 5*PA_A3g(3,1,2) + 5*PA_A3g(3,3,2)) - 
     &            4*(M2**2 - p1p2**2)*sman*(PA_A3g(4,1,2) - PA_A3g(4,3,2)) + 
     &            4*k1p2*sman*(2*M2**2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 
     &               M2*(4*p1p2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) - (-2 + FPSk1 + FPSk2)*PA_A3g(4,1,2) - 2*PA_A3g(4,1,3) + 2*FPSk1*PA_A3g(4,1,3) - 2*PA_A3g(4,3,1) + 
     &                  2*FPSk2*PA_A3g(4,3,1) - 2*PA_A3g(4,3,2) + FPSk1*PA_A3g(4,3,2) + FPSk2*PA_A3g(4,3,2)) + 
     &               p1p2*(2*p1p2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) - FPSk2*PA_A3g(4,1,2) + 2*PA_A3g(4,1,3) + 2*PA_A3g(4,2,1) + 2*PA_A3g(4,2,3) + 2*PA_A3g(4,3,1) - 
     &                  2*FPSk2*PA_A3g(4,3,1) - FPSk2*PA_A3g(4,3,2) + FPSk1*(PA_A3g(4,1,2) - 2*PA_A3g(4,1,3) + PA_A3g(4,3,2)))) - 
     &            4*k1p2**3*(-4*p1p2*PA_A3g(2,1,2) - 4*p1p2*PA_A3g(2,3,2) - 4*M2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 
     &               (-1 + FPSk1)*(-1 + FPSk2)*sman*(PA_A3g(1,1,1) + PA_A3g(1,1,3) + PA_A3g(1,3,1) + PA_A3g(1,3,3) + PA_A3g(3,1,1) + PA_A3g(3,1,3) + PA_A3g(3,3,1) + 
     &                  PA_A3g(3,3,3)) - PA_A3g(4,1,1) + FPSk1*PA_A3g(4,1,1) + 3*FPSk2*PA_A3g(4,1,1) - 3*FPSk1*FPSk2*PA_A3g(4,1,1) - PA_A3g(4,1,3) + 
     &               FPSk1*PA_A3g(4,1,3) - 3*FPSk2*PA_A3g(4,1,3) + 3*FPSk1*FPSk2*PA_A3g(4,1,3) - 6*PA_A3g(4,2,1) + 2*FPSk1*PA_A3g(4,2,1) + 2*FPSk2*PA_A3g(4,2,1) - 
     &               2*FPSk1*FPSk2*PA_A3g(4,2,1) - 6*PA_A3g(4,2,3) + 2*FPSk1*PA_A3g(4,2,3) + 2*FPSk2*PA_A3g(4,2,3) - 2*FPSk1*FPSk2*PA_A3g(4,2,3) - PA_A3g(4,3,1) - 
     &               3*FPSk1*PA_A3g(4,3,1) + FPSk2*PA_A3g(4,3,1) + 3*FPSk1*FPSk2*PA_A3g(4,3,1) - PA_A3g(4,3,3) + 3*FPSk1*PA_A3g(4,3,3) + FPSk2*PA_A3g(4,3,3) - 
     &               3*FPSk1*FPSk2*PA_A3g(4,3,3)) + k1p2**2*(2*p1p2*sman*
     &                (2*(-1 + FPSk1)*PA_A3g(2,1,1) + (-2 + FPSk1 + FPSk2)*PA_A3g(2,1,2) - 2*PA_A3g(2,3,1) + 2*FPSk1*PA_A3g(2,3,1) - 2*PA_A3g(2,3,2) + 
     &                  FPSk1*PA_A3g(2,3,2) + FPSk2*PA_A3g(2,3,2) + 2*PA_A3g(3,1,2) - 2*FPSk2*PA_A3g(3,1,2) + 2*PA_A3g(3,3,2) - 2*FPSk2*PA_A3g(3,3,2)) - 
     &               2*(-1 + FPSk1)*(-1 + 2*FPSk2)*sman*PA_A3g(4,1,1) + 24*p1p2*((-1 + FPSk1)*PA_A3g(4,1,2) - (-1 + FPSk2)*PA_A3g(4,3,2)) + 
     &               2*M2*(sman*(2*(-1 + FPSk1)*PA_A3g(2,1,1) + (-2 + FPSk1 + FPSk2)*PA_A3g(2,1,2) - 2*PA_A3g(2,3,1) + 2*FPSk1*PA_A3g(2,3,1) - 2*PA_A3g(2,3,2) + 
     &                     FPSk1*PA_A3g(2,3,2) + FPSk2*PA_A3g(2,3,2) + 2*PA_A3g(3,1,2) - 2*FPSk2*PA_A3g(3,1,2) + 2*PA_A3g(3,3,2) - 2*FPSk2*PA_A3g(3,3,2)) - 
     &                  12*(-1 + FPSk1)*PA_A3g(4,1,2) + 12*(-1 + FPSk2)*PA_A3g(4,3,2)) - 
     &               (-1 + FPSk2)*sman*((-1 + FPSk1)*PA_A3g(4,1,2) + 4*(-1 + FPSk1)*PA_A3g(4,2,1) - 4*PA_A3g(4,2,3) + 4*FPSk1*PA_A3g(4,2,3) + PA_A3g(4,3,2) - 
     &                  FPSk1*PA_A3g(4,3,2) - 2*PA_A3g(4,3,3) + 4*FPSk1*PA_A3g(4,3,3)))) - 
     &         k1p2*sman*(k1p2**3*(2*(-1 + FPSk1)*(-1 + 2*FPSk2)*PA_A3g(2,1,1) + (-1 + FPSk1)*(-1 + FPSk2)*PA_A3g(2,1,2) + 2*PA_A3g(2,3,1) - 
     &               2*FPSk1*PA_A3g(2,3,1) - 4*FPSk2*PA_A3g(2,3,1) + 4*FPSk1*FPSk2*PA_A3g(2,3,1) + PA_A3g(2,3,2) - FPSk1*PA_A3g(2,3,2) - FPSk2*PA_A3g(2,3,2) + 
     &               FPSk1*FPSk2*PA_A3g(2,3,2) - 4*PA_A3g(3,1,1) + 4*FPSk1*PA_A3g(3,1,1) + 4*FPSk2*PA_A3g(3,1,1) - 4*FPSk1*FPSk2*PA_A3g(3,1,1) - 2*PA_A3g(3,1,2) + 
     &               4*FPSk1*PA_A3g(3,1,2) + 2*FPSk2*PA_A3g(3,1,2) - 4*FPSk1*FPSk2*PA_A3g(3,1,2) - 4*PA_A3g(3,3,1) + 4*FPSk1*PA_A3g(3,3,1) + 4*FPSk2*PA_A3g(3,3,1) - 
     &               4*FPSk1*FPSk2*PA_A3g(3,3,1) - 2*PA_A3g(3,3,2) + 4*FPSk1*PA_A3g(3,3,2) + 2*FPSk2*PA_A3g(3,3,2) - 4*FPSk1*FPSk2*PA_A3g(3,3,2)) + 
     &            4*(M2**2 - p1p2**2)*(PA_A3g(4,1,2) - PA_A3g(4,3,2)) + 
     &            k1p2**2*(2*M2*(2*(-1 + FPSk1)*PA_A3g(2,1,1) + (-2 + FPSk1 + FPSk2)*PA_A3g(2,1,2) - 2*PA_A3g(2,3,1) + 2*FPSk1*PA_A3g(2,3,1) - 2*PA_A3g(2,3,2) + 
     &                  FPSk1*PA_A3g(2,3,2) + FPSk2*PA_A3g(2,3,2) + 2*PA_A3g(3,1,2) - 2*FPSk2*PA_A3g(3,1,2) + 2*PA_A3g(3,3,2) - 2*FPSk2*PA_A3g(3,3,2)) + 
     &               2*p1p2*(2*(-1 + FPSk1)*PA_A3g(2,1,1) + (-2 + FPSk1 + FPSk2)*PA_A3g(2,1,2) - 2*PA_A3g(2,3,1) + 2*FPSk1*PA_A3g(2,3,1) - 2*PA_A3g(2,3,2) + 
     &                  FPSk1*PA_A3g(2,3,2) + FPSk2*PA_A3g(2,3,2) + 2*PA_A3g(3,1,2) - 2*FPSk2*PA_A3g(3,1,2) + 2*PA_A3g(3,3,2) - 2*FPSk2*PA_A3g(3,3,2)) - 
     &               2*PA_A3g(4,1,1) + 2*FPSk1*PA_A3g(4,1,1) + 4*FPSk2*PA_A3g(4,1,1) - 4*FPSk1*FPSk2*PA_A3g(4,1,1) - PA_A3g(4,1,2) + FPSk1*PA_A3g(4,1,2) + 
     &               FPSk2*PA_A3g(4,1,2) - FPSk1*FPSk2*PA_A3g(4,1,2) - 4*PA_A3g(4,2,1) + 4*FPSk1*PA_A3g(4,2,1) + 4*FPSk2*PA_A3g(4,2,1) - 
     &               4*FPSk1*FPSk2*PA_A3g(4,2,1) - 4*PA_A3g(4,2,3) + 4*FPSk1*PA_A3g(4,2,3) + 4*FPSk2*PA_A3g(4,2,3) - 4*FPSk1*FPSk2*PA_A3g(4,2,3) + PA_A3g(4,3,2) - 
     &               FPSk1*PA_A3g(4,3,2) - FPSk2*PA_A3g(4,3,2) + FPSk1*FPSk2*PA_A3g(4,3,2) - 2*PA_A3g(4,3,3) + 4*FPSk1*PA_A3g(4,3,3) + 2*FPSk2*PA_A3g(4,3,3) - 
     &               4*FPSk1*FPSk2*PA_A3g(4,3,3)) + 2*k1p2*(2*M2**2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 
     &               p1p2*(2*p1p2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) - 2*(-1 + FPSk1)*PA_A3g(4,1,1) + 2*PA_A3g(4,1,2) - FPSk1*PA_A3g(4,1,2) - FPSk2*PA_A3g(4,1,2) + 
     &                  2*PA_A3g(4,2,1) + 2*PA_A3g(4,2,3) - 2*PA_A3g(4,3,2) + FPSk1*PA_A3g(4,3,2) + FPSk2*PA_A3g(4,3,2) + 2*PA_A3g(4,3,3) - 2*FPSk2*PA_A3g(4,3,3)) + 
     &               M2*(4*p1p2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 2*(-1 + FPSk1)*PA_A3g(4,1,1) + FPSk1*PA_A3g(4,1,2) - FPSk2*PA_A3g(4,1,2) + FPSk1*PA_A3g(4,3,2) - 
     &                  FPSk2*PA_A3g(4,3,2) - 2*PA_A3g(4,3,3) + 2*FPSk2*PA_A3g(4,3,3)))) - 
     &         k1p1**2*(4*(-1 + FPSk1)*(-1 + FPSk2)*k1p2**3*(PA_A3g(1,1,2) + PA_A3g(1,3,2) - 5*PA_A3g(2,1,1) - PA_A3g(2,1,3) - 5*PA_A3g(2,3,1) - PA_A3g(2,3,3) + 
     &               5*PA_A3g(3,1,2) + 5*PA_A3g(3,3,2)) - k1p2*(2*p1p2*sman*
     &                (2*(-1 + FPSk1)*PA_A3g(2,1,1) + (-2 + FPSk1 + FPSk2)*PA_A3g(2,1,2) - 2*PA_A3g(2,3,1) + 2*FPSk1*PA_A3g(2,3,1) - 2*PA_A3g(2,3,2) + 
     &                  FPSk1*PA_A3g(2,3,2) + FPSk2*PA_A3g(2,3,2) + 2*PA_A3g(3,1,2) - 2*FPSk2*PA_A3g(3,1,2) + 2*PA_A3g(3,3,2) - 2*FPSk2*PA_A3g(3,3,2)) - 
     &               2*(-1 + FPSk1)*(-1 + 2*FPSk2)*sman*PA_A3g(4,1,1) + 24*p1p2*((-1 + FPSk1)*PA_A3g(4,1,2) - (-1 + FPSk2)*PA_A3g(4,3,2)) + 
     &               2*M2*(sman*(2*(-1 + FPSk1)*PA_A3g(2,1,1) + (-2 + FPSk1 + FPSk2)*PA_A3g(2,1,2) - 2*PA_A3g(2,3,1) + 2*FPSk1*PA_A3g(2,3,1) - 2*PA_A3g(2,3,2) + 
     &                     FPSk1*PA_A3g(2,3,2) + FPSk2*PA_A3g(2,3,2) + 2*PA_A3g(3,1,2) - 2*FPSk2*PA_A3g(3,1,2) + 2*PA_A3g(3,3,2) - 2*FPSk2*PA_A3g(3,3,2)) - 
     &                  12*(-1 + FPSk1)*PA_A3g(4,1,2) + 12*(-1 + FPSk2)*PA_A3g(4,3,2)) - 
     &               (-1 + FPSk2)*sman*((-1 + FPSk1)*PA_A3g(4,1,2) + 4*(-1 + FPSk1)*PA_A3g(4,2,1) - 4*PA_A3g(4,2,3) + 4*FPSk1*PA_A3g(4,2,3) + PA_A3g(4,3,2) - 
     &                  FPSk1*PA_A3g(4,3,2) - 2*PA_A3g(4,3,3) + 4*FPSk1*PA_A3g(4,3,3))) + 
     &            2*sman*(2*M2**2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 
     &               p1p2*(2*p1p2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) - 2*(-1 + FPSk1)*PA_A3g(4,1,1) + 2*PA_A3g(4,1,2) - FPSk1*PA_A3g(4,1,2) - FPSk2*PA_A3g(4,1,2) + 
     &                  2*PA_A3g(4,2,1) + 2*PA_A3g(4,2,3) - 2*PA_A3g(4,3,2) + FPSk1*PA_A3g(4,3,2) + FPSk2*PA_A3g(4,3,2) + 2*PA_A3g(4,3,3) - 2*FPSk2*PA_A3g(4,3,3)) + 
     &               M2*(4*p1p2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 2*(-1 + FPSk1)*PA_A3g(4,1,1) + FPSk1*PA_A3g(4,1,2) - FPSk2*PA_A3g(4,1,2) + FPSk1*PA_A3g(4,3,2) - 
     &                  FPSk2*PA_A3g(4,3,2) - 2*PA_A3g(4,3,3) + 2*FPSk2*PA_A3g(4,3,3))) - 
     &            2*k1p2**2*(sman*(4*(-1 + FPSk1)*(-1 + FPSk2)*PA_A3g(1,1,1) + 4*(-1 + FPSk1)*(-1 + FPSk2)*PA_A3g(1,1,3) + 4*PA_A3g(1,3,1) - 
     &                  4*FPSk1*PA_A3g(1,3,1) - 4*FPSk2*PA_A3g(1,3,1) + 4*FPSk1*FPSk2*PA_A3g(1,3,1) + 4*PA_A3g(1,3,3) - 4*FPSk1*PA_A3g(1,3,3) - 
     &                  4*FPSk2*PA_A3g(1,3,3) + 4*FPSk1*FPSk2*PA_A3g(1,3,3) + 2*PA_A3g(2,1,1) - 2*FPSk1*PA_A3g(2,1,1) - 4*FPSk2*PA_A3g(2,1,1) + 
     &                  4*FPSk1*FPSk2*PA_A3g(2,1,1) + PA_A3g(2,1,2) - FPSk1*PA_A3g(2,1,2) - FPSk2*PA_A3g(2,1,2) + FPSk1*FPSk2*PA_A3g(2,1,2) + 2*PA_A3g(2,3,1) - 
     &                  2*FPSk1*PA_A3g(2,3,1) - 4*FPSk2*PA_A3g(2,3,1) + 4*FPSk1*FPSk2*PA_A3g(2,3,1) + PA_A3g(2,3,2) - FPSk1*PA_A3g(2,3,2) - FPSk2*PA_A3g(2,3,2) + 
     &                  FPSk1*FPSk2*PA_A3g(2,3,2) - 2*PA_A3g(3,1,2) + 4*FPSk1*PA_A3g(3,1,2) + 2*FPSk2*PA_A3g(3,1,2) - 4*FPSk1*FPSk2*PA_A3g(3,1,2) + 4*PA_A3g(3,1,3) - 
     &                  4*FPSk1*PA_A3g(3,1,3) - 4*FPSk2*PA_A3g(3,1,3) + 4*FPSk1*FPSk2*PA_A3g(3,1,3) - 2*PA_A3g(3,3,2) + 4*FPSk1*PA_A3g(3,3,2) + 
     &                  2*FPSk2*PA_A3g(3,3,2) - 4*FPSk1*FPSk2*PA_A3g(3,3,2) + 4*PA_A3g(3,3,3) - 4*FPSk1*PA_A3g(3,3,3) - 4*FPSk2*PA_A3g(3,3,3) + 
     &                  4*FPSk1*FPSk2*PA_A3g(3,3,3)) - 4*(4*M2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + 4*p1p2*(PA_A3g(2,1,2) + PA_A3g(2,3,2)) + PA_A3g(4,1,1) - 
     &                  FPSk1*PA_A3g(4,1,1) - 3*FPSk2*PA_A3g(4,1,1) + 3*FPSk1*FPSk2*PA_A3g(4,1,1) + PA_A3g(4,1,3) - FPSk1*PA_A3g(4,1,3) + 3*FPSk2*PA_A3g(4,1,3) - 
     &                  3*FPSk1*FPSk2*PA_A3g(4,1,3) + 6*PA_A3g(4,2,1) - 2*FPSk1*PA_A3g(4,2,1) - 2*FPSk2*PA_A3g(4,2,1) + 2*FPSk1*FPSk2*PA_A3g(4,2,1) + 
     &                  6*PA_A3g(4,2,3) - 2*FPSk1*PA_A3g(4,2,3) - 2*FPSk2*PA_A3g(4,2,3) + 2*FPSk1*FPSk2*PA_A3g(4,2,3) + PA_A3g(4,3,1) + 3*FPSk1*PA_A3g(4,3,1) - 
     &                  FPSk2*PA_A3g(4,3,1) - 3*FPSk1*FPSk2*PA_A3g(4,3,1) + PA_A3g(4,3,3) - 3*FPSk1*PA_A3g(4,3,3) - FPSk2*PA_A3g(4,3,3) + 3*FPSk1*FPSk2*PA_A3g(4,3,3))
     &               )))))/(4*k1k2**2*k1p1*k1p2*sman*dconjg(DenS)))



      end
