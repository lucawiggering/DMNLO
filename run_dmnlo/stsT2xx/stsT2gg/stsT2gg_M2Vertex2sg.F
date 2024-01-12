
      subroutine stst2gg_calc_M2Vertex2sg(stst2gg_M2Vertex2sg)
      implicit none

#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"

      double precision M2
      double complex VSfSf1m2,VSfSf1p2
      double complex DenS, DenSft(2,4,3), DenSfu(2,4,3)
      integer i,j,k
      double precision p1p1,p2p2,p1Pp22,k1p1,k1p2,Krondelta,stst2gg_M2Vertex2sg(3)

      stst2gg_M2Vertex2sg = 0d0

      !Propagator
      DenS = sman

      do i=1,2
      do j=1,4
      do k=1,3
      DenSft(i,j,k)=(tman - MSf2(i,j,k)) !+ cI*MSf(i,j,k)*WSf(i,j,k)
      DenSfu(i,j,k)=(uman - MSf2(i,j,k)) !+ cI*MSf(i,j,k)*WSf(i,j,k)
      enddo
      enddo
      enddo

      M2=MSf2(iSf1,itt1,igen1)

! Values for scalar products for the 2->2 process
      k1k2=sman/2d0
      p1p2=(sman-2*M2)/2d0
      k1p1=(M2-tman)/2d0
      k1p2=(M2-uman)/2d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  s-channel correction !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
! Scalar products for s-channel correction 
      p1p1=M2
      p2p2=M2
      p1Pp22=sman

      call stsT2xx_calc_partial_amplitudes_ASqSqG(VSfSf1m2,VSfSf1p2,isf1,itt1,igen1,isf2,itt2,igen2,p1p1,p2p2,p1Pp22)

      VSfSf1m2=dconjg(VSfSf1m2)
      VSfSf1p2=dconjg(VSfSf1p2)

! antisymmetric octet
      stst2gg_M2Vertex2sg(3) = 2*dble(-(CF*CV*gs**2*GS*Nc*(4*(k1p1 - k1p2)*sman*(FPSk1*(FPSk2*(k1p1 - k1p2)*(k1p1 + k1p2)**2*VSfSf1m2 + k1p1*k1p2**2*(VSfSf1m2 - 3*VSfSf1p2) - 
     -             k1p1**3*VSfSf1p2 - k1p2**3*VSfSf1p2 - k1p1**2*k1p2*(VSfSf1m2 + 3*VSfSf1p2)) + 
     -          FPSk2*(-(k1p1**2*k1p2*(VSfSf1m2 - 3*VSfSf1p2)) + k1p1**3*VSfSf1p2 + k1p2**3*VSfSf1p2 + k1p1*k1p2**2*(VSfSf1m2 + 3*VSfSf1p2))) - 
     -       2*k1k2*(k1p1**3*(2*(-5 - 3*FPSk2 + FPSk1*(-3 + 7*FPSk2))*k1p2*VSfSf1m2 + 
     -             sman*(FPSk1*(-1 + 2*FPSk2)*VSfSf1m2 - FPSk1*VSfSf1p2 + FPSk2*(-VSfSf1m2 + VSfSf1p2))) + 
     -          k1p1**2*(-4*(-5 - 3*FPSk2 + FPSk1*(-3 + 7*FPSk2))*k1p2**2*VSfSf1m2 - 
     -             2*sman*((-2 + FPSk1 + FPSk2)*M2*VSfSf1m2 - (-1 + FPSk1 + FPSk2)*p1p2*VSfSf1m2 + (FPSk1 - FPSk2)*M2*VSfSf1p2 + 
     -                (FPSk1 - FPSk2)*p1p2*VSfSf1p2) + k1p2*sman*(FPSk1*(VSfSf1m2 - 2*FPSk2*VSfSf1m2 - VSfSf1p2) + FPSk2*(VSfSf1m2 + VSfSf1p2))) + 
     -          k1p2**2*sman*(4*M2*VSfSf1m2 - 2*p1p2*VSfSf1m2 - FPSk2*
     -              (2*p1p2*(-VSfSf1m2 + VSfSf1p2) + k1p2*(VSfSf1m2 + VSfSf1p2) + 2*M2*(VSfSf1m2 + VSfSf1p2)) + 
     -             FPSk1*((-1 + 2*FPSk2)*k1p2*VSfSf1m2 + k1p2*VSfSf1p2+ 2*M2*(-VSfSf1m2 + VSfSf1p2) + 2*p1p2*(VSfSf1m2 + VSfSf1p2))) + 
     -          k1p1*k1p2*(2*(-5 - 3*FPSk2 + FPSk1*(-3 + 7*FPSk2))*k1p2**2*VSfSf1m2 + 
     -             2*(-((-2 + FPSk1 + FPSk2)*M2) + (-4 + FPSk1 + FPSk2)*p1p2)*sman*VSfSf1m2 + 
     -             k1p2*sman*(FPSk2*(VSfSf1m2 - VSfSf1p2) + FPSk1*(VSfSf1m2 - 2*FPSk2*VSfSf1m2 + VSfSf1p2)))) + 
     -       k1k2**2*((-1 + FPSk1)*(-1 + FPSk2)*k1p1**2*sman*VSfSf1m2 + 
     -          k1p2*sman*((-1 + FPSk1)*(-1 + FPSk2)*k1p2*VSfSf1m2 + 2*p1p2*((-2 + FPSk1 + FPSk2)*VSfSf1m2 + (FPSk1 - FPSk2)*VSfSf1p2) - 
     -             2*M2*((-2 + FPSk1 + FPSk2)*VSfSf1m2 + (-FPSk1 + FPSk2)*VSfSf1p2)) - 
     -          2*k1p1*(k1p2*(8*(-2 + FPSk1 + FPSk2)*M2 - 8*(-2 + FPSk1 + FPSk2)*p1p2 + (-1 + FPSk1)*(-1 + FPSk2)*sman)*VSfSf1m2 + 
     -             sman*(M2*((-2 + FPSk1 + FPSk2)*VSfSf1m2 + (FPSk1 - FPSk2)*VSfSf1p2) - p1p2*((-2 + FPSk1 + FPSk2)*VSfSf1m2 + (-FPSk1 + FPSk2)*VSfSf1p2))))))/
     -  (2*k1k2*k1p1*k1p2*sman*dconjg(DenS)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Ghosts   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Ghost 1
      stst2gg_M2Vertex2sg(3)=stst2gg_M2Vertex2sg(3)-(1d0-FPSk1)*(1d0-FPSk2)*2*dble(-(CF*CV*GS**3*Nc*(tman - uman)*(-2*M2*VSfSf1p2 + uman*(-VSfSf1m2 + VSfSf1p2) + tman*(VSfSf1m2 + VSfSf1p2))*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(4*DenS*dconjg(DenS)))

! Ghost 2
      stst2gg_M2Vertex2sg(3)=stst2gg_M2Vertex2sg(3)-(1d0-FPSk1)*(1d0-FPSk2)*2*dble(-(CF*CV*GS**3*Nc*(tman - uman)*(tman*(VSfSf1m2 - VSfSf1p2) + 2*M2*VSfSf1p2 - uman*(VSfSf1m2 + VSfSf1p2))*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(4*DenS*dconjg(DenS)))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  t-channel correction upper vertex !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sclar products for t-channel correction, upper vertex 
      p1p1=M2
      p2p2=tman
      p1Pp22=0d0
      
      call stsT2xx_calc_partial_amplitudes_ASqSqG(VSfSf1m2,VSfSf1p2,isf1,itt1,igen1,isf2,itt2,igen2,p1p1,p2p2,p1Pp22)

      VSfSf1m2=dconjg(VSfSf1m2)
      VSfSf1p2=dconjg(VSfSf1p2)
 
! singlet 
      stst2gg_M2Vertex2sg(1) = 2*dble( -(gs**2*GS*(-1 + Nc**2)*(-4*FPSk1*FPSk2*k1p1*(k1p1**4 + 2*k1p1**2*k1p2**2 + 4*k1p1*k1p2**3 + k1p2**4)*VSfSf1m2 + 
     -       (-1 + FPSk1)*(-1 + FPSk2)*k1k2**3*k1p1*k1p2*(VSfSf1m2 - VSfSf1p2) + 
     -       2*k1k2*k1p1*(FPSk1*((-1 + 2*FPSk2)*k1p1**3*VSfSf1m2 + k1p1*k1p2**2*(5*(-1 + 2*FPSk2)*VSfSf1m2 - 2*FPSk2*VSfSf1p2) + 
     -             k1p1**2*((-3 + 4*FPSk2)*k1p2*VSfSf1m2 + 2*p1p2*VSfSf1m2 - FPSk2*k1p2*VSfSf1p2) + 
     -             k1p2**2*((-3 + 4*FPSk2)*k1p2*VSfSf1m2 + 2*(2*M2 + p1p2)*VSfSf1m2 - FPSk2*k1p2*VSfSf1p2)) + 
     -          FPSk2*(-(k1p1**3*VSfSf1m2) + k1p1*k1p2**2*(-5*VSfSf1m2+ 2*VSfSf1p2) + k1p1**2*(2*p1p2*VSfSf1m2 + k1p2*(-3*VSfSf1m2 + VSfSf1p2)) + 
     -             k1p2**2*(2*(2*M2 + p1p2)*VSfSf1m2 + k1p2*(-3*VSfSf1m2+ VSfSf1p2)))) - 
     -       k1k2**2*((-1 + FPSk1)*(-1 + FPSk2)*k1p1**3*VSfSf1m2 + 4*k1p2*M2**2*VSfSf1m2 + 
     -          k1p1**2*(2*(-2 + FPSk1 + FPSk2)*p1p2*VSfSf1m2 + 3*(-1 + FPSk1)*(-1 + FPSk2)*k1p2*(2*VSfSf1m2 - VSfSf1p2)) + 
     -          k1p1*(4*p1p2**2*VSfSf1m2 + (-1 + FPSk1)*(-1 + FPSk2)*k1p2**2*(5*VSfSf1m2 - 3*VSfSf1p2) + 
     -             2*k1p2*(p1p2*((-4 + FPSk1 + FPSk2)*VSfSf1m2 + VSfSf1p2 - FPSk1*VSfSf1p2) + M2*(2*(-2 + FPSk1 + FPSk2)*VSfSf1m2 + VSfSf1p2 - FPSk1*VSfSf1p2)))
     -          )))/(6*k1k2**2*k1p1*k1p2*dconjg(DenSft(isf2,itt2,igen2))))

! symmetric octet
      stst2gg_M2Vertex2sg(2) = 2*dble(-(CF*gs**2*GS*(-4 + Nc**2)*(-4*FPSk1*FPSk2*k1p1*(k1p1**4 + 2*k1p1**2*k1p2**2 + 4*k1p1*k1p2**3 + k1p2**4)*VSfSf1m2 + 
     -       (-1 + FPSk1)*(-1 + FPSk2)*k1k2**3*k1p1*k1p2*(VSfSf1m2 - VSfSf1p2) + 
     -       2*k1k2*k1p1*(FPSk1*((-1 + 2*FPSk2)*k1p1**3*VSfSf1m2 + k1p1*k1p2**2*(5*(-1 + 2*FPSk2)*VSfSf1m2 - 2*FPSk2*VSfSf1p2) + 
     -             k1p1**2*((-3 + 4*FPSk2)*k1p2*VSfSf1m2 + 2*p1p2*VSfSf1m2 - FPSk2*k1p2*VSfSf1p2) + 
     -             k1p2**2*((-3 + 4*FPSk2)*k1p2*VSfSf1m2 + 2*(2*M2 + p1p2)*VSfSf1m2 - FPSk2*k1p2*VSfSf1p2)) + 
     -          FPSk2*(-(k1p1**3*VSfSf1m2) + k1p1*k1p2**2*(-5*VSfSf1m2+ 2*VSfSf1p2) + k1p1**2*(2*p1p2*VSfSf1m2 + k1p2*(-3*VSfSf1m2 + VSfSf1p2)) + 
     -             k1p2**2*(2*(2*M2 + p1p2)*VSfSf1m2 + k1p2*(-3*VSfSf1m2+ VSfSf1p2)))) - 
     -       k1k2**2*((-1 + FPSk1)*(-1 + FPSk2)*k1p1**3*VSfSf1m2 + 4*k1p2*M2**2*VSfSf1m2 + 
     -          k1p1**2*(2*(-2 + FPSk1 + FPSk2)*p1p2*VSfSf1m2 + 3*(-1 + FPSk1)*(-1 + FPSk2)*k1p2*(2*VSfSf1m2 - VSfSf1p2)) + 
     -          k1p1*(4*p1p2**2*VSfSf1m2 + (-1 + FPSk1)*(-1 + FPSk2)*k1p2**2*(5*VSfSf1m2 - 3*VSfSf1p2) + 
     -             2*k1p2*(p1p2*((-4 + FPSk1 + FPSk2)*VSfSf1m2 + VSfSf1p2 - FPSk1*VSfSf1p2) + M2*(2*(-2 + FPSk1 + FPSk2)*VSfSf1m2 + VSfSf1p2 - FPSk1*VSfSf1p2)))
     -          )))/(2*k1k2**2*k1p1*k1p2*dconjg(DenSft(isf2,itt2,igen2))))

! anti symmetric octet
      stst2gg_M2Vertex2sg(3) =  stst2gg_M2Vertex2sg(3) + 2*dble( -(CF*CV*gs**2*GS*Nc*(4*FPSk1*FPSk2*k1p1*(k1p1**4 + 2*k1p1**2*k1p2**2 - 4*k1p1*k1p2**3 + k1p2**4)*sman*VSfSf1m2 + 
     -       2*k1k2**3*k1p1*k1p2*((-1 + FPSk1)*(-1 + FPSk2)*k1p1*(VSfSf1m2 - VSfSf1p2) - (-1 + FPSk1)*(-1 + FPSk2)*k1p2*(VSfSf1m2 - VSfSf1p2) + 
     -          2*(M2 - p1p2)*((-2 + FPSk1 + FPSk2)*VSfSf1m2 + VSfSf1p2 - FPSk1*VSfSf1p2)) - 
     -       2*k1k2*k1p1*(FPSk2*(-(k1p1**3*(6*k1p2 + sman)*VSfSf1m2) + k1p1*k1p2**2*(2*k1p2*VSfSf1m2 + 3*sman*VSfSf1m2 - 2*sman*VSfSf1p2) + 
     -             k1p2**2*sman*(2*(-2*M2 + p1p2)*VSfSf1m2 + k1p2*(-VSfSf1m2 + VSfSf1p2)) + 
     -             k1p1**2*(4*k1p2**2*VSfSf1m2 + 2*p1p2*sman*VSfSf1m2 + k1p2*sman*(-VSfSf1m2 + VSfSf1p2))) + 
     -          FPSk1*(k1p1**3*((-6 + 8*FPSk2)*k1p2 + (-1 + 2*FPSk2)*sman)*VSfSf1m2 + 
     -             k1p1*k1p2**2*((2 - 8*FPSk2)*k1p2*VSfSf1m2 + (3 - 6*FPSk2)*sman*VSfSf1m2 + 2*FPSk2*sman*VSfSf1p2) - 
     -             k1p1**2*(4*(-1 + FPSk2)*k1p2**2*VSfSf1m2 - 2*p1p2*sman*VSfSf1m2 + k1p2*sman*(VSfSf1m2 - 2*FPSk2*VSfSf1m2 + FPSk2*VSfSf1p2)) + 
     -             k1p2**2*(-((k1p2 + 4*M2 - 2*p1p2)*sman*VSfSf1m2) + FPSk2*k1p2*(4*k1p2*VSfSf1m2 + 2*sman*VSfSf1m2 - sman*VSfSf1p2)))) + 
     -       k1k2**2*(-4*k1p2*M2**2*sman*VSfSf1m2 + k1p1**3*(4*(-1 + FPSk1*FPSk2)*k1p2*VSfSf1m2 + (-1 + FPSk1)*(-1 + FPSk2)*sman*VSfSf1m2 - 
     -             2*(-1 + FPSk1)*(1 + FPSk2)*k1p2*VSfSf1p2) + k1p1**2*
     -           (2*(-2 + FPSk1 + FPSk2)*p1p2*sman*VSfSf1m2 + k1p2*
     -              (-4*(-2 + FPSk1 + FPSk2)*M2*VSfSf1m2 + 4*(-4 + FPSk1+FPSk2)*p1p2*VSfSf1m2 - (-1 + FPSk1)*(-1 + FPSk2)*sman*(2*VSfSf1m2 - VSfSf1p2)) + 
     -             k1p2**2*(4*(1 + FPSk1 + FPSk2 - 3*FPSk1*FPSk2)*VSfSf1m2 + 2*(-1 + FPSk1)*(1 + 3*FPSk2)*VSfSf1p2)) + 
     -          k1p1*(16*k1p2**2*M2*VSfSf1m2 - 8*k1p2**2*p1p2*VSfSf1m2 +k1p2**2*sman*VSfSf1m2 + 8*k1p2*M2*sman*VSfSf1m2 - 4*k1p2*p1p2*sman*VSfSf1m2 + 
     -             4*p1p2**2*sman*VSfSf1m2 - k1p2**2*sman*VSfSf1p2 - 2*k1p2*M2*sman*VSfSf1p2 + 2*k1p2*p1p2*sman*VSfSf1p2 + 
     -             FPSk2*k1p2*(2*(-2*M2 + p1p2)*sman*VSfSf1m2 - 4*k1p2**2*(VSfSf1m2 - VSfSf1p2) + 
     -                k1p2*(-8*M2*VSfSf1m2 + 8*p1p2*VSfSf1m2 - sman*VSfSf1m2 + sman*VSfSf1p2)) + 
     -             FPSk1*k1p2*(k1p2*(-8*M2*VSfSf1m2 + 8*p1p2*VSfSf1m2 + (-1 + FPSk2)*sman*(VSfSf1m2 - VSfSf1p2)) + 
     -                k1p2**2*((-4 + 8*FPSk2)*VSfSf1m2 - 4*FPSk2*VSfSf1p2) + 2*sman*(p1p2*(VSfSf1m2 - VSfSf1p2) + M2*(-2*VSfSf1m2 + VSfSf1p2)))))))/
     -  (2*k1k2**2*k1p1*k1p2*sman*dconjg(DenSft(isf2,itt2,igen2))))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  t-channel correction lower vertex !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sclar products for t-channel correction, lower vertex 
      p1p1=tman
      p2p2=M2
      p1Pp22=0d0
      
      call stsT2xx_calc_partial_amplitudes_ASqSqG(VSfSf1m2,VSfSf1p2,isf1,itt1,igen1,isf2,itt2,igen2,p1p1,p2p2,p1Pp22)

      VSfSf1m2=dconjg(VSfSf1m2)
      VSfSf1p2=dconjg(VSfSf1p2)
  
! Contraction of t-channel correction (lower vertex) with ...
! singlet
      stst2gg_M2Vertex2sg(1) = stst2gg_M2Vertex2sg(1) + 2*dble( -(gs**2*GS*(-1 + Nc**2)*(-4*FPSk1*FPSk2*k1p1*(k1p1**4 + 2*k1p1**2*k1p2**2 + 4*k1p1*k1p2**3 + k1p2**4)*VSfSf1m2 + 
     -       (-1 + FPSk1)*(-1 + FPSk2)*k1k2**3*k1p1*k1p2*(VSfSf1m2 + VSfSf1p2) + 
     -       2*k1k2*k1p1*(-(FPSk2*(k1p1**3 + 5*k1p1*k1p2**2 + k1p1**2*(3*k1p2 - 2*p1p2) + k1p2**2*(3*k1p2 - 4*M2 - 2*p1p2))*VSfSf1m2) + 
     -          FPSk1*((-1 + 2*FPSk2)*k1p1**3*VSfSf1m2 + k1p1*k1p2**2*(5*(-1 + 2*FPSk2)*VSfSf1m2 + 2*(-1 + FPSk2)*VSfSf1p2) + 
     -             k1p1**2*((-3 + 4*FPSk2)*k1p2*VSfSf1m2 + 2*p1p2*VSfSf1m2 + (-1 + FPSk2)*k1p2*VSfSf1p2) + 
     -             k1p2**2*((-3 + 4*FPSk2)*k1p2*VSfSf1m2 + 2*(2*M2 + p1p2)*VSfSf1m2 + (-1 + FPSk2)*k1p2*VSfSf1p2))) - 
     -       k1k2**2*((-1 + FPSk1)*(-1 + FPSk2)*k1p1**3*VSfSf1m2 + 4*k1p2*M2**2*VSfSf1m2 + 
     -          k1p1**2*(2*(-2 + FPSk1 + FPSk2)*p1p2*VSfSf1m2 + 3*(-1 + FPSk1)*(-1 + FPSk2)*k1p2*(2*VSfSf1m2 + VSfSf1p2)) + 
     -          k1p1*(4*p1p2**2*VSfSf1m2 + (-1 + FPSk1)*(-1 + FPSk2)*k1p2**2*(5*VSfSf1m2 + 3*VSfSf1p2) + 
     -             2*k1p2*(2*(-2 + FPSk1 + FPSk2)*M2*VSfSf1m2 + (-4 + FPSk1 + FPSk2)*p1p2*VSfSf1m2 + (-1 + FPSk2)*M2*VSfSf1p2 + (-1 + FPSk2)*p1p2*VSfSf1p2)))))/
     -  (6*k1k2**2*k1p1*k1p2*dconjg(DenSft(isf1,itt1,igen1))))

! symmetric octet
      stst2gg_M2Vertex2sg(2) =  stst2gg_M2Vertex2sg(2) + 2*dble(  -(CF*gs**2*GS*(-4 + Nc**2)*(-4*FPSk1*FPSk2*k1p1*(k1p1**4 + 2*k1p1**2*k1p2**2 + 4*k1p1*k1p2**3 + k1p2**4)*VSfSf1m2 + 
     -       (-1 + FPSk1)*(-1 + FPSk2)*k1k2**3*k1p1*k1p2*(VSfSf1m2 + VSfSf1p2) + 
     -       2*k1k2*k1p1*(-(FPSk2*(k1p1**3 + 5*k1p1*k1p2**2 + k1p1**2*(3*k1p2 - 2*p1p2) + k1p2**2*(3*k1p2 - 4*M2 - 2*p1p2))*VSfSf1m2) + 
     -          FPSk1*((-1 + 2*FPSk2)*k1p1**3*VSfSf1m2 + k1p1*k1p2**2*(5*(-1 + 2*FPSk2)*VSfSf1m2 + 2*(-1 + FPSk2)*VSfSf1p2) + 
     -             k1p1**2*((-3 + 4*FPSk2)*k1p2*VSfSf1m2 + 2*p1p2*VSfSf1m2 + (-1 + FPSk2)*k1p2*VSfSf1p2) + 
     -             k1p2**2*((-3 + 4*FPSk2)*k1p2*VSfSf1m2 + 2*(2*M2 + p1p2)*VSfSf1m2 + (-1 + FPSk2)*k1p2*VSfSf1p2))) - 
     -       k1k2**2*((-1 + FPSk1)*(-1 + FPSk2)*k1p1**3*VSfSf1m2 + 4*k1p2*M2**2*VSfSf1m2 + 
     -          k1p1**2*(2*(-2 + FPSk1 + FPSk2)*p1p2*VSfSf1m2 + 3*(-1 + FPSk1)*(-1 + FPSk2)*k1p2*(2*VSfSf1m2 + VSfSf1p2)) + 
     -          k1p1*(4*p1p2**2*VSfSf1m2 + (-1 + FPSk1)*(-1 + FPSk2)*k1p2**2*(5*VSfSf1m2 + 3*VSfSf1p2) + 
     -             2*k1p2*(2*(-2 + FPSk1 + FPSk2)*M2*VSfSf1m2 + (-4 + FPSk1 + FPSk2)*p1p2*VSfSf1m2 + (-1 + FPSk2)*M2*VSfSf1p2 + (-1 + FPSk2)*p1p2*VSfSf1p2)))))/
     -  (2*k1k2**2*k1p1*k1p2*dconjg(DenSft(isf1,itt1,igen1))))

! antisymmetric octet
      stst2gg_M2Vertex2sg(3) =  stst2gg_M2Vertex2sg(3) + 2*dble((CF*CV*gs**2*GS*Nc*(-2*sman*(-(sman*(M2 + tman)) + FPSk1*(M2 - tman)*(2*M2 - sman - 2*uman))*
     -       (-(sman*(2*M2*VSfSf1m2 + 2*tman*VSfSf1m2 - M2*VSfSf1p2 + tman*VSfSf1p2)) + 
     -         FPSk2*(M2 - tman)*(4*M2*VSfSf1m2 - 4*uman*VSfSf1m2 - sman*(2*VSfSf1m2 + VSfSf1p2)))*DenSfu(isf1,itt1,igen1) + 
     -      DenSft(isf1,itt1,igen1)*(sman*(-(sman*(-6*M2 + 2*sman + tman + uman)) + 
     -            FPSk1*(4*M2**2 + sman*(tman + uman) + 2*(tman**2 + uman**2) - 2*M2*(sman + 2*(tman + uman))))*
     -          (-(sman*(2*sman*VSfSf1m2 + tman*VSfSf1m2 + uman*VSfSf1m2 + uman*VSfSf1p2 - M2*(6*VSfSf1m2 + VSfSf1p2))) + 
     -            FPSk2*(4*M2**2*VSfSf1m2 + sman*tman*VSfSf1m2 + 2*(tman**2 + uman**2)*VSfSf1m2 + sman*uman*(VSfSf1m2 + VSfSf1p2) - 
     -               M2*(4*(tman + uman)*VSfSf1m2 + sman*(2*VSfSf1m2 + VSfSf1p2)))) + 
     -         sman*(8*M2**2*(3*(-2 + FPSk1 + FPSk2)*sman - (FPSk1*(-1 + FPSk2) - FPSk2)*(tman - uman))*VSfSf1m2 - 
     -            2*(tman - uman)*(-(FPSk2*tman*(3*tman + uman)) + FPSk1*((-3 + 4*FPSk2)*tman**2 + (-1 + 2*FPSk2)*tman*uman - 2*FPSk2*uman**2))*VSfSf1m2 + 
     -            2*sman**3*((-2 + FPSk1 + FPSk2)*VSfSf1m2 + (-1 + FPSk2)*VSfSf1p2) + 
     -            sman**2*(-(uman*((5 + FPSk1*(-5 + FPSk2) - 5*FPSk2)*VSfSf1m2 + (-1 + FPSk1)*(-1 + FPSk2)*VSfSf1p2)) + 
     -               tman*((-7 + FPSk1 + FPSk2 + FPSk1*FPSk2)*VSfSf1m2 + (-1 + FPSk1)*(-1 + FPSk2)*VSfSf1p2)) - 
     -            sman*(tman - uman)*(2*(-1 + FPSk1*FPSk2)*tman*VSfSf1m2 + 2*FPSk2*uman*VSfSf1m2 + (1 + FPSk1)*(-1 + FPSk2)*tman*VSfSf1p2 + 
     -               2*FPSk1*uman*(VSfSf1m2 - 2*FPSk2*VSfSf1m2 + VSfSf1p2 - FPSk2*VSfSf1p2)) + 
     -            M2*(2*(tman - uman)*(FPSk1*(-7 + 10*FPSk2)*tman - FPSk1*(1 + 2*FPSk2)*uman - FPSk2*(7*tman + uman))*VSfSf1m2 - 
     -               2*sman**2*(7*(-2 + FPSk1 + FPSk2)*VSfSf1m2 + 4*(-1 + FPSk2)*VSfSf1p2) - 
     -               sman*(2*(-11 + 3*FPSk2 + FPSk1*(3 + FPSk2))*tman*VSfSf1m2 - 2*(13 + FPSk1*(-9 + FPSk2) - 9*FPSk2)*uman*VSfSf1m2 + 
     -                  (-1 + FPSk1)*(-1 + FPSk2)*tman*VSfSf1p2 + (-1 + FPSk1 + FPSk2 - FPSk1*FPSk2)*uman*VSfSf1p2)))*DenSfu(isf1,itt1,igen1))))/
     -  (4*sman**3*dconjg(DenSft(isf1,itt1,igen1))*DenSft(isf1,itt1,igen1)*DenSfu(isf1,itt1,igen1)))


  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  u-channel correction upper vertex !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sclar products for u-channel correction, upper vertex 
      p1p1=M2
      p2p2=uman
      p1Pp22=0d0
      
      call stsT2xx_calc_partial_amplitudes_ASqSqG(VSfSf1m2,VSfSf1p2,isf1,itt1,igen1,isf2,itt2,igen2,p1p1,p2p2,p1Pp22)

      VSfSf1m2=dconjg(VSfSf1m2)
      VSfSf1p2=dconjg(VSfSf1p2)

! Contraction of u-channel correction (upper vertex) with ...
! s-channel
      stst2gg_M2Vertex2sg(1) = stst2gg_M2Vertex2sg(1) + 2*dble( -(gs**2*GS*(-1 + Nc**2)*(-4*FPSk1*FPSk2*k1p2*(k1p1**4 + 4*k1p1**3*k1p2 + 2*k1p1**2*k1p2**2 + k1p2**4)*VSfSf1m2 + 
     -       (-1 + FPSk1)*(-1 + FPSk2)*k1k2**3*k1p1*k1p2*(VSfSf1m2 - VSfSf1p2) + 
     -       2*k1k2*k1p2*(-(FPSk2*(3*k1p1**3 + 3*k1p1*k1p2**2 + k1p2**2*(k1p2 - 2*p1p2) + k1p1**2*(5*k1p2 - 2*(2*M2 + p1p2)))*VSfSf1m2) + 
     -          FPSk1*(k1p2**2*((-1 + 2*FPSk2)*k1p2 + 2*p1p2)*VSfSf1m2+ k1p1**3*((-3 + 4*FPSk2)*VSfSf1m2 + VSfSf1p2 - FPSk2*VSfSf1p2) + 
     -             k1p1*k1p2**2*((-3 + 4*FPSk2)*VSfSf1m2 + VSfSf1p2 - FPSk2*VSfSf1p2) + 
     -             k1p1**2*(5*(-1 + 2*FPSk2)*k1p2*VSfSf1m2 + 2*(2*M2 + p1p2)*VSfSf1m2 - 2*(-1 + FPSk2)*k1p2*VSfSf1p2))) - 
     -       k1k2**2*(k1p2*((-1 + FPSk1)*k1p2 + 2*p1p2)*((-1 + FPSk2)*k1p2 + 2*p1p2)*VSfSf1m2 + 
     -          (-1 + FPSk1)*(-1 + FPSk2)*k1p1**2*k1p2*(5*VSfSf1m2 - 3*VSfSf1p2) + 
     -          k1p1*(4*M2**2*VSfSf1m2 + 3*(-1 + FPSk1)*(-1 + FPSk2)*k1p2**2*(2*VSfSf1m2 - VSfSf1p2) + 
     -             2*k1p2*(p1p2*((-4 + FPSk1 + FPSk2)*VSfSf1m2 + VSfSf1p2 - FPSk2*VSfSf1p2) + M2*(2*(-2 + FPSk1 + FPSk2)*VSfSf1m2 + VSfSf1p2 - FPSk2*VSfSf1p2)))
     -          )))/(6*k1k2**2*k1p1*k1p2*dconjg(DenSfu(isf2,itt2,igen2))))


! symmetric octet
      stst2gg_M2Vertex2sg(2) =  stst2gg_M2Vertex2sg(2) + 2*dble( -(CF*gs**2*GS*(-4 + Nc**2)*(-4*FPSk1*FPSk2*k1p2*(k1p1**4 + 4*k1p1**3*k1p2 + 2*k1p1**2*k1p2**2 + k1p2**4)*VSfSf1m2 + 
     -       (-1 + FPSk1)*(-1 + FPSk2)*k1k2**3*k1p1*k1p2*(VSfSf1m2 - VSfSf1p2) + 
     -       2*k1k2*k1p2*(-(FPSk2*(3*k1p1**3 + 3*k1p1*k1p2**2 + k1p2**2*(k1p2 - 2*p1p2) + k1p1**2*(5*k1p2 - 2*(2*M2 + p1p2)))*VSfSf1m2) + 
     -          FPSk1*(k1p2**2*((-1 + 2*FPSk2)*k1p2 + 2*p1p2)*VSfSf1m2+ k1p1**3*((-3 + 4*FPSk2)*VSfSf1m2 + VSfSf1p2 - FPSk2*VSfSf1p2) + 
     -             k1p1*k1p2**2*((-3 + 4*FPSk2)*VSfSf1m2 + VSfSf1p2 - FPSk2*VSfSf1p2) + 
     -             k1p1**2*(5*(-1 + 2*FPSk2)*k1p2*VSfSf1m2 + 2*(2*M2 + p1p2)*VSfSf1m2 - 2*(-1 + FPSk2)*k1p2*VSfSf1p2))) - 
     -       k1k2**2*(k1p2*((-1 + FPSk1)*k1p2 + 2*p1p2)*((-1 + FPSk2)*k1p2 + 2*p1p2)*VSfSf1m2 + 
     -          (-1 + FPSk1)*(-1 + FPSk2)*k1p1**2*k1p2*(5*VSfSf1m2 - 3*VSfSf1p2) + 
     -          k1p1*(4*M2**2*VSfSf1m2 + 3*(-1 + FPSk1)*(-1 + FPSk2)*k1p2**2*(2*VSfSf1m2 - VSfSf1p2) + 
     -             2*k1p2*(p1p2*((-4 + FPSk1 + FPSk2)*VSfSf1m2 + VSfSf1p2 - FPSk2*VSfSf1p2) + M2*(2*(-2 + FPSk1 + FPSk2)*VSfSf1m2 + VSfSf1p2 - FPSk2*VSfSf1p2)))
     -          )))/(2*k1k2**2*k1p1*k1p2*dconjg(DenSfu(isf2,itt2,igen2))))

! antisymmetric octet
      stst2gg_M2Vertex2sg(3) =  stst2gg_M2Vertex2sg(3) +2*dble((CF*CV*gs**2*GS*Nc*(-4*FPSk1*FPSk2*k1p2*(k1p1**4 - 4*k1p1**3*k1p2 + 2*k1p1**2*k1p2**2 + k1p2**4)*sman*VSfSf1m2 + 
     -      2*k1k2**3*k1p1*k1p2*((-1 + FPSk1)*(-1 + FPSk2)*k1p1*(VSfSf1m2 - VSfSf1p2) - (-1 + FPSk1)*(-1 + FPSk2)*k1p2*(VSfSf1m2 - VSfSf1p2) - 
     -         2*(M2 - p1p2)*((-2 + FPSk1 + FPSk2)*VSfSf1m2 + VSfSf1p2 - FPSk2*VSfSf1p2)) + 
     -      2*k1k2*k1p2*(FPSk2*(k1p1**3*(2*k1p2 - sman) - k1p2**2*(k1p2 - 2*p1p2)*sman - k1p1*k1p2**2*(6*k1p2 + sman) + 
     -            k1p1**2*(4*k1p2**2 + 3*k1p2*sman + 2*(-2*M2 + p1p2)*sman))*VSfSf1m2 + 
     -         FPSk1*(-(k1p2**2*(k1p2 - 2*p1p2)*sman*VSfSf1m2) + 
     -            k1p1**2*(4*k1p2**2*VSfSf1m2 + 2*(-2*M2 + p1p2)*sman*VSfSf1m2 + k1p2*sman*(3*VSfSf1m2 - 2*VSfSf1p2)) + 
     -            k1p1*k1p2**2*(-6*k1p2*VSfSf1m2 + sman*(-VSfSf1m2 + VSfSf1p2)) + k1p1**3*(2*k1p2*VSfSf1m2 + sman*(-VSfSf1m2 + VSfSf1p2)) + 
     -            FPSk2*(k1p1 - k1p2)*(4*k1p1**3*VSfSf1m2 - 2*k1p2**2*sman*VSfSf1m2 + k1p1*k1p2*(-8*k1p2*VSfSf1m2 - 4*sman*VSfSf1m2 + sman*VSfSf1p2) - 
     -               k1p1**2*(4*k1p2*VSfSf1m2 - 2*sman*VSfSf1m2 + sman*VSfSf1p2)))) - 
     -      k1k2**2*(4*k1p1**2*k1p2**2*VSfSf1m2 - 4*k1p1*k1p2**3*VSfSf1m2 + 16*k1p1**2*k1p2*M2*VSfSf1m2 + 8*k1p1*k1p2**2*M2*VSfSf1m2 - 
     -         8*k1p1**2*k1p2*p1p2*VSfSf1m2 - 16*k1p1*k1p2**2*p1p2*VSfSf1m2 + k1p1**2*k1p2*sman*VSfSf1m2 - 2*k1p1*k1p2**2*sman*VSfSf1m2 + 
     -         k1p2**3*sman*VSfSf1m2 + 8*k1p1*k1p2*M2*sman*VSfSf1m2 - 4*k1p1*M2**2*sman*VSfSf1m2 - 4*k1p1*k1p2*p1p2*sman*VSfSf1m2 - 
     -         4*k1p2**2*p1p2*sman*VSfSf1m2 + 4*k1p2*p1p2**2*sman*VSfSf1m2 - 2*k1p1**2*k1p2**2*VSfSf1p2 + 2*k1p1*k1p2**3*VSfSf1p2 - 
     -         k1p1**2*k1p2*sman*VSfSf1p2 + k1p1*k1p2**2*sman*VSfSf1p2 - 2*k1p1*k1p2*M2*sman*VSfSf1p2 + 2*k1p1*k1p2*p1p2*sman*VSfSf1p2 + 
     -         FPSk1*k1p2*(k1p2*((-1 + FPSk2)*k1p2 + 2*p1p2)*sman*VSfSf1m2 + 4*k1p1**3*((-1 + 2*FPSk2)*VSfSf1m2 + VSfSf1p2 - FPSk2*VSfSf1p2) + 
     -            k1p1**2*((4 - 12*FPSk2)*k1p2*VSfSf1m2 - 8*M2*VSfSf1m2 + 8*p1p2*VSfSf1m2 - sman*VSfSf1m2 + FPSk2*sman*VSfSf1m2 + 
     -               6*(-1 + FPSk2)*k1p2*VSfSf1p2 + sman*VSfSf1p2 - FPSk2*sman*VSfSf1p2) + 
     -            k1p1*(2*(-2*M2 + p1p2)*sman*VSfSf1m2 + FPSk2*k1p2*(2*k1p2 - sman)*(2*VSfSf1m2 - VSfSf1p2) + 2*k1p2**2*VSfSf1p2 + 
     -               k1p2*(-4*M2*VSfSf1m2 + 4*p1p2*VSfSf1m2 + 2*sman*VSfSf1m2 - sman*VSfSf1p2))) - 
     -         FPSk2*k1p2*(4*k1p1**3*VSfSf1m2 + k1p2*(k1p2 - 2*p1p2)*sman*VSfSf1m2 + 
     -            k1p1**2*(8*M2*VSfSf1m2 - 8*p1p2*VSfSf1m2 + sman*VSfSf1m2 - sman*VSfSf1p2 - 2*k1p2*(2*VSfSf1m2 + VSfSf1p2)) + 
     -            k1p1*(2*k1p2**2*VSfSf1p2 + 2*sman*(2*M2*VSfSf1m2 - p1p2*VSfSf1m2 - M2*VSfSf1p2 + p1p2*VSfSf1p2) + 
     -               k1p2*(4*M2*VSfSf1m2 - 4*p1p2*VSfSf1m2 - 2*sman*VSfSf1m2 + sman*VSfSf1p2))))))/(2.*k1k2**2*k1p1*k1p2*sman*dconjg(DenSfu(isf2,itt2,igen2))))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  u-channel correction lower vertex !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sclar products for u-channel correction, lower vertex 
      p1p1=uman
      p2p2=M2
      p1Pp22=0d0
      
      call stsT2xx_calc_partial_amplitudes_ASqSqG(VSfSf1m2,VSfSf1p2,isf1,itt1,igen1,isf2,itt2,igen2,p1p1,p2p2,p1Pp22)

      VSfSf1m2=dconjg(VSfSf1m2)
      VSfSf1p2=dconjg(VSfSf1p2)
  
! singlet
      stst2gg_M2Vertex2sg(1) = stst2gg_M2Vertex2sg(1) + 2*dble( -(gs**2*GS*(-1 + Nc**2)*(-4*FPSk1*FPSk2*k1p2*(k1p1**4 + 4*k1p1**3*k1p2 + 2*k1p1**2*k1p2**2 + k1p2**4)*VSfSf1m2 + 
     -       (-1 + FPSk1)*(-1 + FPSk2)*k1k2**3*k1p1*k1p2*(VSfSf1m2 + VSfSf1p2) - 
     -       k1k2**2*(k1p2*((-1 + FPSk1)*k1p2 + 2*p1p2)*((-1 + FPSk2)*k1p2 + 2*p1p2)*VSfSf1m2 + 
     -          (-1 + FPSk1)*(-1 + FPSk2)*k1p1**2*k1p2*(5*VSfSf1m2 + 3*VSfSf1p2) + 
     -          k1p1*(4*M2**2*VSfSf1m2 + 3*(-1 + FPSk1)*(-1 + FPSk2)*k1p2**2*(2*VSfSf1m2 + VSfSf1p2) + 
     -             2*k1p2*(2*(-2 + FPSk1 + FPSk2)*M2*VSfSf1m2 + (-4 + FPSk1 + FPSk2)*p1p2*VSfSf1m2 + (-1 + FPSk1)*M2*VSfSf1p2 + (-1 + FPSk1)*p1p2*VSfSf1p2))) + 
     -       2*k1k2*k1p2*(FPSk1*(k1p2**2*((-1 + 2*FPSk2)*k1p2 + 2*p1p2)*VSfSf1m2 + k1p1**3*((-3 + 4*FPSk2)*VSfSf1m2 + FPSk2*VSfSf1p2) + 
     -             k1p1*k1p2**2*((-3 + 4*FPSk2)*VSfSf1m2 + FPSk2*VSfSf1p2) + 
     -             k1p1**2*(5*(-1 + 2*FPSk2)*k1p2*VSfSf1m2 + 2*(2*M2 + p1p2)*VSfSf1m2 + 2*FPSk2*k1p2*VSfSf1p2)) - 
     -          FPSk2*(k1p2**2*(k1p2 - 2*p1p2)*VSfSf1m2 + k1p1**3*(3*VSfSf1m2 + VSfSf1p2) + k1p1*k1p2**2*(3*VSfSf1m2 + VSfSf1p2) + 
     -             k1p1**2*(-2*(2*M2 + p1p2)*VSfSf1m2 + k1p2*(5*VSfSf1m2+ 2*VSfSf1p2))))))/(6*k1k2**2*k1p1*k1p2*dconjg(DenSfu(isf1,itt1,igen1))))

! symmetric octet
      stst2gg_M2Vertex2sg(2) =  stst2gg_M2Vertex2sg(2) + 2*dble( -(CF*gs**2*GS*(-4 + Nc**2)*(-4*FPSk1*FPSk2*k1p2*(k1p1**4 + 4*k1p1**3*k1p2 + 2*k1p1**2*k1p2**2 + k1p2**4)*VSfSf1m2 + 
     -       (-1 + FPSk1)*(-1 + FPSk2)*k1k2**3*k1p1*k1p2*(VSfSf1m2 + VSfSf1p2) - 
     -       k1k2**2*(k1p2*((-1 + FPSk1)*k1p2 + 2*p1p2)*((-1 + FPSk2)*k1p2 + 2*p1p2)*VSfSf1m2 + 
     -          (-1 + FPSk1)*(-1 + FPSk2)*k1p1**2*k1p2*(5*VSfSf1m2 + 3*VSfSf1p2) + 
     -          k1p1*(4*M2**2*VSfSf1m2 + 3*(-1 + FPSk1)*(-1 + FPSk2)*k1p2**2*(2*VSfSf1m2 + VSfSf1p2) + 
     -             2*k1p2*(2*(-2 + FPSk1 + FPSk2)*M2*VSfSf1m2 + (-4 + FPSk1 + FPSk2)*p1p2*VSfSf1m2 + (-1 + FPSk1)*M2*VSfSf1p2 + (-1 + FPSk1)*p1p2*VSfSf1p2))) + 
     -       2*k1k2*k1p2*(FPSk1*(k1p2**2*((-1 + 2*FPSk2)*k1p2 + 2*p1p2)*VSfSf1m2 + k1p1**3*((-3 + 4*FPSk2)*VSfSf1m2 + FPSk2*VSfSf1p2) + 
     -             k1p1*k1p2**2*((-3 + 4*FPSk2)*VSfSf1m2 + FPSk2*VSfSf1p2) + 
     -             k1p1**2*(5*(-1 + 2*FPSk2)*k1p2*VSfSf1m2 + 2*(2*M2 + p1p2)*VSfSf1m2 + 2*FPSk2*k1p2*VSfSf1p2)) - 
     -          FPSk2*(k1p2**2*(k1p2 - 2*p1p2)*VSfSf1m2 + k1p1**3*(3*VSfSf1m2 + VSfSf1p2) + k1p1*k1p2**2*(3*VSfSf1m2 + VSfSf1p2) + 
     -             k1p1**2*(-2*(2*M2 + p1p2)*VSfSf1m2 + k1p2*(5*VSfSf1m2+ 2*VSfSf1p2))))))/(2*k1k2**2*k1p1*k1p2*dconjg(DenSfu(isf1,itt1,igen1))))

! antisymmetric octet
      stst2gg_M2Vertex2sg(3) =  stst2gg_M2Vertex2sg(3) +2*dble(  (CF*CV*gs**2*GS*Nc*(sman*(-(sman*(-6*M2 + 2*sman + tman + uman)) + 
     -         FPSk2*(4*M2**2 + sman*(tman + uman) + 2*(tman**2 + uman**2) - 2*M2*(sman + 2*(tman + uman))))*
     -       (-(sman*(2*sman*VSfSf1m2 + tman*VSfSf1m2 + uman*VSfSf1m2 + tman*VSfSf1p2 - M2*(6*VSfSf1m2 + VSfSf1p2))) + 
     -         FPSk1*(4*M2**2*VSfSf1m2 + sman*uman*VSfSf1m2 + 2*(tman**2 + uman**2)*VSfSf1m2 + sman*tman*(VSfSf1m2 + VSfSf1p2) - 
     -            M2*(4*(tman + uman)*VSfSf1m2 + sman*(2*VSfSf1m2 + VSfSf1p2))))*DenSfu(isf1,itt1,igen1) + 
     -      DenSft(isf1,itt1,igen1)*(-2*sman*(FPSk2*(2*M2 - sman - 2*tman)*(M2 - uman) - sman*(M2 + uman))*
     -          (-(sman*(2*M2*VSfSf1m2 + 2*uman*VSfSf1m2 - M2*VSfSf1p2 + uman*VSfSf1p2)) + 
     -            FPSk1*(M2 - uman)*(4*M2*VSfSf1m2 - 4*tman*VSfSf1m2 - sman*(2*VSfSf1m2 + VSfSf1p2))) + 
     -         sman*(8*M2**2*(3*(-2 + FPSk1 + FPSk2)*sman + (FPSk1*(-1 + FPSk2) - FPSk2)*(tman - uman))*VSfSf1m2 - 
     -            2*(tman - uman)*(FPSk1*uman*(tman + 3*uman) + FPSk2*uman*(tman + 3*uman) + 2*FPSk1*FPSk2*(tman**2 - tman*uman - 2*uman**2))*VSfSf1m2 + 
     -            2*sman**3*((-2 + FPSk1 + FPSk2)*VSfSf1m2 + (-1 + FPSk1)*VSfSf1p2) + 
     -            sman**2*(-(tman*((5 + FPSk1*(-5 + FPSk2) - 5*FPSk2)*VSfSf1m2 + (-1 + FPSk1)*(-1 + FPSk2)*VSfSf1p2)) + 
     -               uman*((-7 + FPSk1 + FPSk2 + FPSk1*FPSk2)*VSfSf1m2 + (-1 + FPSk1)*(-1 + FPSk2)*VSfSf1p2)) + 
     -            M2*(2*(tman - uman)*(FPSk2*(tman + 7*uman) + FPSk1*(tman + 2*FPSk2*tman + 7*uman - 10*FPSk2*uman))*VSfSf1m2 - 
     -               2*sman**2*(7*(-2 + FPSk1 + FPSk2)*VSfSf1m2 + 4*(-1 + FPSk1)*VSfSf1p2) + 
     -               sman*(2*(13 + FPSk1*(-9 + FPSk2) - 9*FPSk2)*tman*VSfSf1m2 - 2*(-11 + 3*FPSk2 + FPSk1*(3 + FPSk2))*uman*VSfSf1m2 + 
     -                  (-1 + FPSk1)*(-1 + FPSk2)*tman*VSfSf1p2 - (-1 + FPSk1)*(-1 + FPSk2)*uman*VSfSf1p2)) - 
     -            sman*(tman - uman)*(-2*FPSk2*tman*VSfSf1m2 + 2*uman*VSfSf1m2 - 2*FPSk2*tman*VSfSf1p2 + uman*VSfSf1p2 + FPSk2*uman*VSfSf1p2 + 
     -               2*FPSk1*tman*((-1 + 2*FPSk2)*VSfSf1m2 + FPSk2*VSfSf1p2) - FPSk1*uman*(VSfSf1p2 + FPSk2*(2*VSfSf1m2 + VSfSf1p2))))*DenSfu(isf1,itt1,igen1)))
     -    )/(4*sman**3*dconjg(DenSfu(isf1,itt1,igen1))*DenSft(isf1,itt1,igen1)*DenSfu(isf1,itt1,igen1)))

      end 
      