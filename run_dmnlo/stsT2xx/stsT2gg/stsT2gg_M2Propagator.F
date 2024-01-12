      double complex function stsT2xx_PiSquarkLoop(myP2,i,j,type1,type2,gen1,gen2)
!i,type1,gen1 are sfermion index, type and geneartion of incoming squark
!j,type2,gen2 are sfermion index, type and geneartion of outgoing squark
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"
      
      integer i,j,k,type1,type2,gen1,gen2
      double precision ML0,ML02,A01,myP2,Krondelta

      if (type1.eq.type2 .and. gen1.eq.gen2) then

      stsT2xx_PiSquarkLoop=0d0
! ===== Quark-gluino loop =====  I
      call myBget(myP2,MGl2,Mf2(type1,gen1))
      call myAget(Mf2(type1,gen1))
     
      stsT2xx_PiSquarkLoop=stsT2xx_PiSquarkLoop-GS**2*CF/(8d0*Pi**2)*(Krondelta(i,j)*(A01(MGl)+xa0)+     
     &      xb0*(Krondelta(i,j)*(Mf2(type1,gen1)+Mgl2-myP2)-2d0*Mgl*Mf(type1,gen1)*(USf(i,1,type1,gen1)*USf(j,2,type1,gen1)+USf(j,1,type1,gen1)*USf(i,2,type1,gen1))))   

! ===== Squark-gluon loop ===== II
      call myBget(myP2,0d0,MSf2(i,type1,gen1))
      call myAget(MSf2(i,type1,gen1))

      stsT2xx_PiSquarkLoop=stsT2xx_PiSquarkLoop+GS**2*CF*Krondelta(i,j)/(16d0*Pi**2)*(xa0-2d0*xb0*(MSf2(i,type1,gen1)+myP2))

! ===== Squark one-point loop =====  III
      do k=1,2
        call myAget(MSf2(k,type1,gen1))

        stsT2xx_PiSquarkLoop=stsT2xx_PiSquarkLoop+GS**2*CF/(16d0*Pi**2)*(USf(i,1,type1,gen1)*USf(k,1,type1,gen1)-USf(i,2,type1,gen1)*USf(k,2,type1,gen1))*
     &   (USf(j,1,type1,gen1)*USf(k,1,type1,gen1)-USf(j,2,type1,gen1)*USf(k,2,type1,gen1))*xa0
      
      enddo

! ===== Gluon one-point loop vanishes =====  III
      
      endif
      end 



      subroutine  stst2gg_calc_M2Propagator(stst2gg_M2Propagator)
      implicit none

#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      integer i,j,k
      double precision M2,Krondelta,stst2gg_M2Propagator(3)
      double complex DenS, DenSft(2,4,3), DenSfu(2,4,3)
      double complex PiGT,stsT2xx_PiSquarkLoop,stsT2xx_PiGluonLoopT
      double complex PiSquarkUman,PiSquarkTman

      M2=MSf2(iSf1,itt1,igen1)
    
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

      PiGT=stsT2xx_PiGluonLoopT(sman)

      PiSquarkTman=stsT2xx_PiSquarkLoop(tman,isf1,isf2,itt1,itt2,igen1,igen2) 
      PiSquarkUman=stsT2xx_PiSquarkLoop(uman,isf1,isf2,itt1,itt2,igen1,igen2)  

!********************************************
! Gluon propagator corrections
!*************************************+****

! ! !---------------------2 Re(Ms^{NLO} Ms^{Tree\ast})-----------
! 
      stst2gg_M2Propagator(3) = 2*dble((-PiGT/DenS)* (CF*CV*gs**4*Nc*(16*(-2 + FPSk1 + FPSk2)*M2**3*sman - 2*(-2 + FPSk1 + FPSk2)*sman**3*(tman + uman) + 
     -      2*(-FPSk2 + FPSk1*(-1 + 2*FPSk2))*sman*(tman - uman)**2*(tman + uman) + 
     -      sman**2*((5 + FPSk1*(-5 + FPSk2) - 5*FPSk2)*tman**2 - 2*(1 + FPSk1*(-3 + FPSk2) - 3*FPSk2)*tman*uman + 
     -         (5 + FPSk1*(-5 + FPSk2) - 5*FPSk2)*uman**2) - 2*M2**2*
     -       (10*(-2 + FPSk1 + FPSk2)*sman**2 - (5 + FPSk1 + FPSk2 + FPSk1*FPSk2)*(tman - uman)**2 + 8*(-2 + FPSk1 + FPSk2)*sman*(tman + uman)) + 
     -      2*(tman - uman)**2*(FPSk1*tman*uman + (5 + FPSk2)*tman*uman + FPSk1*FPSk2*(2*tman**2 - 3*tman*uman + 2*uman**2)) + 
     -      2*M2*(2*(-2 + FPSk1 + FPSk2)*sman**3 + 6*(-2 + FPSk1 + FPSk2)*sman**2*(tman + uman) - 
     -         (5 + FPSk1 + FPSk2 + FPSk1*FPSk2)*(tman - uman)**2*(tman + uman) - 
     -         2*sman*((6 - 5*FPSk2 + FPSk1*(-5 + 2*FPSk2))*tman**2 + 2*(-2 + FPSk1*(3 - 2*FPSk2) + 3*FPSk2)*tman*uman + 
     -            (6 - 5*FPSk2 + FPSk1*(-5 + 2*FPSk2))*uman**2)))*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/
     -  (4*DenS*sman*(M2 - tman)*(M2 - uman)))

!! --------------------   -2 Re(Mgh1^{NLO} MGh1^{Tree\ast})  -------------------------------
      stst2gg_M2Propagator(3) = stst2gg_M2Propagator(3) -(1d0-FPSk1)*(1d0-FPSk2)*2*dble((-PiGT/DenS)*(CF*CV*GS**4*Nc*(-tman + uman)**2*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(4*DenS*dconjg(DenS)))

!! --------------------   -2 Re(Mgh1^{NLO} MGh1^{Tree\ast})    -------------------------------
      stst2gg_M2Propagator(3) = stst2gg_M2Propagator(3) -(1d0-FPSk1)*(1d0-FPSk2)*2*dble((-PiGT/DenS)*(CF*CV*GS**4*Nc*(tman- uman)**2*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(4*DenS*dconjg(DenS)))

!********************************************
! Squark propagator corrections 
!********************************************

! singlet --------------    2 Re(Mt^((1)prop) MSinglet^\ast)
      stst2gg_M2Propagator(1) = 2*dble((-PiSquarkTman/DenSft(isf1,itt1,igen1))*(gs**4*(-1 + Nc**2)*(sman*(sman*(52*M2**3 - M2**2*(29*sman + 52*tman + 32*uman) - 
     -            tman*(4*sman**2 + 4*sman*tman + tman**2 + 9*sman*uman + 6*tman*uman + 5*uman**2) + 
     -            M2*(4*sman**2 + 33*sman*tman + 17*tman**2 + 9*sman*uman + 22*tman*uman + 5*uman**2)) + 
     -         FPSk2*(M2 - tman)*(24*M2**3 - sman**2*(2*tman + 3*uman) - 4*M2**2*(5*sman + 11*tman + 7*uman) - sman*(5*tman**2 + 6*tman*uman + 9*uman**2) - 
     -            2*(tman**3 + 3*tman**2*uman + 5*tman*uman**2 + 3*uman**3) + M2*(5*sman**2 + 20*sman*(tman + uman) + 4*(5*tman**2 + 8*tman*uman + 5*uman**2))))
     -        + FPSk1*(M2 - tman)*(sman*(24*M2**3 - sman**2*(2*tman + 3*uman) - 4*M2**2*(5*sman + 11*tman + 7*uman) - 
     -            sman*(5*tman**2 + 6*tman*uman + 9*uman**2) - 2*(tman**3 + 3*tman**2*uman + 5*tman*uman**2 + 3*uman**3) + 
     -            M2*(5*sman**2 + 20*sman*(tman + uman) + 4*(5*tman**2 + 8*tman*uman + 5*uman**2))) + 
     -         FPSk2*(32*M2**4 + sman**3*uman - 8*M2**3*(5*sman + 6*tman + 10*uman) + sman**2*(tman**2 + 6*tman*uman + 5*uman**2) + 
     -            4*M2**2*(3*sman**2 + 12*sman*tman + 8*tman**2 + 18*sman*uman + 20*tman*uman + 20*uman**2) + 
     -            4*sman*(tman**3 + 2*tman**2*uman + 5*tman*uman**2 + 2*uman**3) + 4*(tman**4 + 2*tman**2*uman**2 + 4*tman*uman**3 + uman**4) - 
     -            M2*(sman**3 + 8*sman**2*(tman + 2*uman) + 4*sman*(5*tman**2 + 14*tman*uman + 11*uman**2) + 
     -               16*(tman**3 + tman**2*uman + 4*tman*uman**2 + 2*uman**3)))))*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/
     -  (12*sman**2*(M2 - tman)**2*(M2 - uman)))

! singlet -------------- 2 Re(Mu^((1)prop) MSinglet^\ast)
      stst2gg_M2Propagator(1) = stst2gg_M2Propagator(1) + 2*dble((-PiSquarkUman/DenSfu(isf1,itt1,igen1))*(gs**4*(-1 + Nc**2)*(sman*(sman*(52*M2**3 - M2**2*(29*sman + 32*tman + 52*uman) - 
     -            uman*(4*sman**2 + 9*sman*tman + 5*tman**2 + 4*sman*uman + 6*tman*uman + uman**2) + 
     -            M2*(4*sman**2 + 9*sman*tman + 5*tman**2 + 33*sman*uman + 22*tman*uman + 17*uman**2)) + 
     -         FPSk2*(M2 - uman)*(24*M2**3 - sman**2*(3*tman + 2*uman) - 4*M2**2*(5*sman + 7*tman + 11*uman) - sman*(9*tman**2 + 6*tman*uman + 5*uman**2) - 
     -            2*(3*tman**3 + 5*tman**2*uman + 3*tman*uman**2 + uman**3) + M2*(5*sman**2 + 20*sman*(tman + uman) + 4*(5*tman**2 + 8*tman*uman + 5*uman**2))))
     -        + FPSk1*(M2 - uman)*(sman*(24*M2**3 - sman**2*(3*tman + 2*uman) - 4*M2**2*(5*sman + 7*tman + 11*uman) - 
     -            sman*(9*tman**2 + 6*tman*uman + 5*uman**2) - 2*(3*tman**3 + 5*tman**2*uman + 3*tman*uman**2 + uman**3) + 
     -            M2*(5*sman**2 + 20*sman*(tman + uman) + 4*(5*tman**2 + 8*tman*uman + 5*uman**2))) + 
     -         FPSk2*(32*M2**4 + sman**3*tman - 8*M2**3*(5*sman + 10*tman + 6*uman) + sman**2*(5*tman**2 + 6*tman*uman + uman**2) + 
     -            4*M2**2*(3*sman**2 + 18*sman*tman + 20*tman**2 + 12*sman*uman + 20*tman*uman + 8*uman**2) + 
     -            4*sman*(2*tman**3 + 5*tman**2*uman + 2*tman*uman**2 + uman**3) + 4*(tman**4 + 4*tman**3*uman + 2*tman**2*uman**2 + uman**4) - 
     -            M2*(sman**3 + 8*sman**2*(2*tman + uman) + 4*sman*(11*tman**2 + 14*tman*uman + 5*uman**2) + 
     -               16*(2*tman**3 + 4*tman**2*uman + tman*uman**2 + uman**3)))))*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/
     -  (12*sman**2*(M2 - tman)*(M2 - uman)**2))
     

! symmetric octet --------------    2 Re(Mt^((s))prop) MSOctet^\ast)
      stst2gg_M2Propagator(2) = 2*dble((-PiSquarkTman/DenSft(isf1,itt1,igen1))* (CF*gs**4*(-4 + Nc**2)*(sman*(sman*(52*M2**3 - M2**2*(29*sman + 52*tman + 32*uman) - 
     -            tman*(4*sman**2 + 4*sman*tman + tman**2 + 9*sman*uman + 6*tman*uman + 5*uman**2) + 
     -            M2*(4*sman**2 + 33*sman*tman + 17*tman**2 + 9*sman*uman + 22*tman*uman + 5*uman**2)) + 
     -         FPSk2*(M2 - tman)*(24*M2**3 - sman**2*(2*tman + 3*uman) - 4*M2**2*(5*sman + 11*tman + 7*uman) - sman*(5*tman**2 + 6*tman*uman + 9*uman**2) - 
     -            2*(tman**3 + 3*tman**2*uman + 5*tman*uman**2 + 3*uman**3) + M2*(5*sman**2 + 20*sman*(tman + uman) + 4*(5*tman**2 + 8*tman*uman + 5*uman**2))))
     -        + FPSk1*(M2 - tman)*(sman*(24*M2**3 - sman**2*(2*tman + 3*uman) - 4*M2**2*(5*sman + 11*tman + 7*uman) - 
     -            sman*(5*tman**2 + 6*tman*uman + 9*uman**2) - 2*(tman**3 + 3*tman**2*uman + 5*tman*uman**2 + 3*uman**3) + 
     -            M2*(5*sman**2 + 20*sman*(tman + uman) + 4*(5*tman**2 + 8*tman*uman + 5*uman**2))) + 
     -         FPSk2*(32*M2**4 + sman**3*uman - 8*M2**3*(5*sman + 6*tman + 10*uman) + sman**2*(tman**2 + 6*tman*uman + 5*uman**2) + 
     -            4*M2**2*(3*sman**2 + 12*sman*tman + 8*tman**2 + 18*sman*uman + 20*tman*uman + 20*uman**2) + 
     -            4*sman*(tman**3 + 2*tman**2*uman + 5*tman*uman**2 + 2*uman**3) + 4*(tman**4 + 2*tman**2*uman**2 + 4*tman*uman**3 + uman**4) - 
     -            M2*(sman**3 + 8*sman**2*(tman + 2*uman) + 4*sman*(5*tman**2 + 14*tman*uman + 11*uman**2) + 
     -               16*(tman**3 + tman**2*uman + 4*tman*uman**2 + 2*uman**3)))))*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/
     -  (4*sman**2*(M2 - tman)**2*(M2 - uman)))

! symmetric octet -------------- 2 Re(Mu^(((8s))prop) MSOctet^\ast)
      stst2gg_M2Propagator(2) = stst2gg_M2Propagator(2) + 2*dble((-PiSquarkUman/DenSfu(isf1,itt1,igen1))* (CF*gs**4*(-4 + Nc**2)*(sman*(sman*(52*M2**3 - M2**2*(29*sman + 32*tman + 52*uman) - 
     -            uman*(4*sman**2 + 9*sman*tman + 5*tman**2 + 4*sman*uman + 6*tman*uman + uman**2) + 
     -            M2*(4*sman**2 + 9*sman*tman + 5*tman**2 + 33*sman*uman + 22*tman*uman + 17*uman**2)) + 
     -         FPSk2*(M2 - uman)*(24*M2**3 - sman**2*(3*tman + 2*uman) - 4*M2**2*(5*sman + 7*tman + 11*uman) - sman*(9*tman**2 + 6*tman*uman + 5*uman**2) - 
     -            2*(3*tman**3 + 5*tman**2*uman + 3*tman*uman**2 + uman**3) + M2*(5*sman**2 + 20*sman*(tman + uman) + 4*(5*tman**2 + 8*tman*uman + 5*uman**2))))
     -        + FPSk1*(M2 - uman)*(sman*(24*M2**3 - sman**2*(3*tman + 2*uman) - 4*M2**2*(5*sman + 7*tman + 11*uman) - 
     -            sman*(9*tman**2 + 6*tman*uman + 5*uman**2) - 2*(3*tman**3 + 5*tman**2*uman + 3*tman*uman**2 + uman**3) + 
     -            M2*(5*sman**2 + 20*sman*(tman + uman) + 4*(5*tman**2 + 8*tman*uman + 5*uman**2))) + 
     -         FPSk2*(32*M2**4 + sman**3*tman - 8*M2**3*(5*sman + 10*tman + 6*uman) + sman**2*(5*tman**2 + 6*tman*uman + uman**2) + 
     -            4*M2**2*(3*sman**2 + 18*sman*tman + 20*tman**2 + 12*sman*uman + 20*tman*uman + 8*uman**2) + 
     -            4*sman*(2*tman**3 + 5*tman**2*uman + 2*tman*uman**2 + uman**3) + 4*(tman**4 + 4*tman**3*uman + 2*tman**2*uman**2 + uman**4) - 
     -            M2*(sman**3 + 8*sman**2*(2*tman + uman) + 4*sman*(11*tman**2 + 14*tman*uman + 5*uman**2) + 
     -               16*(2*tman**3 + 4*tman**2*uman + tman*uman**2 + uman**3)))))*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/
     -  (4*sman**2*(M2 - tman)*(M2 - uman)**2))

! antisymmetric octet -------------- 2 Re(Mt^((8A)prop) MAOctet^\ast)
      stst2gg_M2Propagator(3) = stst2gg_M2Propagator(3) +  2*dble((-PiSquarkTman/DenSft(isf1,itt1,igen1))*( -(CF*CV*gs**4*Nc*(8*M2**4*((6 + FPSk1 + FPSk2)*sman - (FPSk1 + FPSk2 + FPSk1*FPSk2)*(tman - uman)) - 
     -       2*M2**3*((-2 + 5*FPSk1 + 5*FPSk2)*sman**2 + sman*((35 + FPSk1*(11 - 5*FPSk2) + 11*FPSk2)*tman + (37 + FPSk1 + FPSk2 + 5*FPSk1*FPSk2)*uman) - 
     -          (tman - uman)*(FPSk1*(11 + 10*FPSk2)*tman + FPSk1*(5 + 6*FPSk2)*uman + FPSk2*(11*tman + 5*uman))) + 
     -       tman*(-4*sman**4 + 2*(-2 + FPSk1 + FPSk2)*sman**3*tman + 
     -          sman**2*((-1 - FPSk1*(-5 + FPSk2) + 5*FPSk2)*tman**2 + (-3 + FPSk1)*(-3 + FPSk2)*tman*uman + 4*uman**2) - 
     -          2*sman*tman*(tman - uman)*(FPSk1*(-1 + 2*FPSk2)*tman + uman + FPSk1*(-2 + 3*FPSk2)*uman - FPSk2*(tman + 2*uman)) - 
     -          2*tman*(tman - uman)*(FPSk1*uman*(3*tman + uman) + FPSk2*uman*(3*tman + uman) + 2*FPSk1*FPSk2*(tman**2 - tman*uman + 2*uman**2))) + 
     -       M2**2*(2*(-10 + FPSk1 + FPSk2)*sman**3 + sman**2*(-(FPSk1*(-21 + FPSk2)*tman) + 21*(-1 + FPSk2)*tman + (25 + FPSk1 + FPSk2 + FPSk1*FPSk2)*uman) + 
     -          2*sman*((10 + 15*FPSk1 + 15*FPSk2 - 12*FPSk1*FPSk2)*tman**2 + (49 - 5*FPSk2 + FPSk1*(-5 + 9*FPSk2))*tman*uman + 
     -             (13 + 2*FPSk2 + FPSk1*(2 + 3*FPSk2))*uman**2) - 
     -          2*(tman - uman)*(FPSk2*(10*tman**2 + 13*tman*uman + uman**2) + 
     -             FPSk1*(10*(1 + FPSk2)*tman**2 + (13 + 10*FPSk2)*tman*uman + (1 + 4*FPSk2)*uman**2))) + 
     -       2*M2*(2*sman**4 - 2*(-6 + FPSk1 + FPSk2)*sman**3*tman + 
     -          sman**2*((1 + FPSk1*(-8 + FPSk2) - 8*FPSk2)*tman**2 + (-9 + FPSk1 + FPSk2 - FPSk1*FPSk2)*tman*uman - 2*uman**2) + 
     -          sman*tman*((1 + 9*FPSk1*(-1 + FPSk2) - 9*FPSk2)*tman**2 + (-11 + FPSk1*(5 - 3*FPSk2) + 5*FPSk2)*tman*uman - 2*(7 + 3*FPSk1*FPSk2)*uman**2) + 
     -          tman*(tman - uman)*(FPSk2*(3*tman**2 + 11*tman*uman + 2*uman**2) + 
     -             FPSk1*((3 + 6*FPSk2)*tman**2 + (11 + 2*FPSk2)*tman*uman + 2*(1 + 4*FPSk2)*uman**2))))*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*
     -     Krondelta(itt1,itt2))/(4*sman**2*(M2 - tman)**2*(M2 - uman))))

! antisymmetric octet -------------- 2 Re(Mu^((8A)prop) MAOctet^\ast)
      stst2gg_M2Propagator(3) = stst2gg_M2Propagator(3) + 2*dble((-PiSquarkUman/DenSfu(isf1,itt1,igen1))*( -(CF*CV*gs**4*Nc*(8*M2**4*((6 + FPSk1 + FPSk2)*sman + (FPSk1 + FPSk2 + FPSk1*FPSk2)*(tman - uman)) - 
     -       2*M2**3*((-2 + 5*FPSk1 + 5*FPSk2)*sman**2 + sman*((37 + FPSk1 + FPSk2 + 5*FPSk1*FPSk2)*tman + (35 + FPSk1*(11 - 5*FPSk2) + 11*FPSk2)*uman) + 
     -          (tman - uman)*(FPSk1*(5 + 6*FPSk2)*tman + FPSk1*(11 + 10*FPSk2)*uman + FPSk2*(5*tman + 11*uman))) + 
     -       uman*(-4*sman**4 + 2*(-2 + FPSk1 + FPSk2)*sman**3*uman + 
     -          2*sman*(tman - uman)*uman*(tman - 2*FPSk2*tman + FPSk1*(-2 + 3*FPSk2)*tman - FPSk2*uman + FPSk1*(-1 + 2*FPSk2)*uman) + 
     -          sman**2*(4*tman**2 + (-3 + FPSk1)*(-3 + FPSk2)*tman*uman + (-1 - FPSk1*(-5 + FPSk2) + 5*FPSk2)*uman**2) + 
     -          2*(tman - uman)*uman*(FPSk2*tman*(tman + 3*uman) + FPSk1*((1 + 4*FPSk2)*tman**2 + (3 - 2*FPSk2)*tman*uman + 2*FPSk2*uman**2))) + 
     -       M2**2*(2*(-10 + FPSk1 + FPSk2)*sman**3 + sman**2*((25 + FPSk1 + FPSk2 + FPSk1*FPSk2)*tman - FPSk1*(-21 + FPSk2)*uman + 21*(-1 + FPSk2)*uman) + 
     -          2*sman*((13 + 2*FPSk2 + FPSk1*(2 + 3*FPSk2))*tman**2 + (49 - 5*FPSk2 + FPSk1*(-5 + 9*FPSk2))*tman*uman + 
     -             (10 + 15*FPSk1 + 15*FPSk2 - 12*FPSk1*FPSk2)*uman**2) + 
     -          2*(tman - uman)*(FPSk2*(tman**2 + 13*tman*uman + 10*uman**2) + 
     -             FPSk1*((1 + 4*FPSk2)*tman**2 + (13 + 10*FPSk2)*tman*uman + 10*(1 + FPSk2)*uman**2))) + 
     -       2*M2*(2*sman**4 - 2*(-6 + FPSk1 + FPSk2)*sman**3*uman + 
     -          sman*uman*(-2*(7 + 3*FPSk1*FPSk2)*tman**2 + (-11 + FPSk1*(5 - 3*FPSk2) + 5*FPSk2)*tman*uman + (1 + 9*FPSk1*(-1 + FPSk2) - 9*FPSk2)*uman**2) + 
     -          sman**2*(-2*tman**2 + (-9 + FPSk1 + FPSk2 - FPSk1*FPSk2)*tman*uman + (1 + FPSk1*(-8 + FPSk2) - 8*FPSk2)*uman**2) - 
     -          (tman - uman)*uman*(FPSk2*(2*tman**2 + 11*tman*uman + 3*uman**2) + 
     -             FPSk1*((2 + 8*FPSk2)*tman**2 + (11 + 2*FPSk2)*tman*uman + 3*(1 + 2*FPSk2)*uman**2))))*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*
     -     Krondelta(itt1,itt2))/(4*sman**2*(M2 - tman)*(M2 - uman)**2)))
 


      end