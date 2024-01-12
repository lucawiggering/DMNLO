      subroutine  stst2gg_calc_M2PropagatorCounter(stst2gg_M2PropagatorCounter) 
      implicit none 
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCouplings.h"
#include "../../util/DMNLO_RenScheme.h"

      double complex DenS, DenSft(2,4,3), DenSfu(2,4,3)
      double precision M2,Krondelta,stst2gg_M2PropagatorCounter(3)
      integer i,j,k

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


!********************************************
! Gluon propagator corrections
!********************************************

      stst2gg_M2PropagatorCounter(3) = 2*dble((-dZG)*(CF*CV*gs**4*Nc*(16*(-2 + FPSk1 + FPSk2)*M2**3*sman - 2*(-2 + FPSk1 + FPSk2)*sman**3*(tman + uman) + 
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
      stst2gg_M2PropagatorCounter(3) = stst2gg_M2PropagatorCounter(3) -(1d0-FPSk1)*(1d0-FPSk2)*2*dble((-dZG)*(CF*CV*GS**4*Nc*(-tman + uman)**2*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(4d0*DenS*dconjg(DenS)))

!! --------------------   -2 Re(Mgh1^{NLO} MGh1^{Tree\ast})    -------------------------------
      stst2gg_M2PropagatorCounter(3) = stst2gg_M2PropagatorCounter(3) -(1d0-FPSk1)*(1d0-FPSk2)*2*dble((-dZG)*(CF*CV*GS**4*Nc*(tman- uman)**2*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(4d0*DenS*dconjg(DenS)))
      
!********************************************
! Squark propagator corrections 
!********************************************
      
! singlet --------------    2 Re(Mt^((1)prop) MSinglet^\ast)

      stst2gg_M2PropagatorCounter(1) = 2*dble((-(dZSf(isf1,isf2,itt1,igen1)-dMSf2(isf1,itt1,igen1)/DenSft(isf1,itt1,igen1)))*
     -     (gs**4*(-1 + Nc**2)*(sman*(sman*(52*M2**3 - M2**2*(29*sman + 52*tman + 32*uman) - 
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
      stst2gg_M2PropagatorCounter(1) = stst2gg_M2PropagatorCounter(1) + 2*dble((-(dZSf(isf1,isf2,itt1,igen1)-dMSf2(isf1,itt1,igen1)/DenSfu(isf1,itt1,igen1)))*
     -     (gs**4*(-1 + Nc**2)*(sman*(sman*(52*M2**3 - M2**2*(29*sman + 32*tman + 52*uman) - 
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


! symmetric octet --------------    2 Re(Mt^((1)prop) MSOctet^\ast)
      stst2gg_M2PropagatorCounter(2) = 2*dble((-(dZSf(isf1,isf2,itt1,igen1)-dMSf2(isf1,itt1,igen1)/DenSft(isf1,itt1,igen1)))*
     -     (CF*gs**4*(-4 + Nc**2)*(sman*(sman*(52*M2**3 - M2**2*(29*sman + 52*tman + 32*uman) - 
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

! symmetric octet -------------- 2 Re(Mu^((1)prop) MSOctet^\ast)
      stst2gg_M2PropagatorCounter(2) =  stst2gg_M2PropagatorCounter(2) + 2*dble((-(dZSf(isf1,isf2,itt1,igen1)-dMSf2(isf1,itt1,igen1)/DenSfu(isf1,itt1,igen1))
     -         *(CF*gs**4*(-4 + Nc**2)*(sman*(sman*(52*M2**3 - M2**2*(29*sman + 32*tman + 52*uman) - 
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
     -  (4*sman**2*(M2 - tman)*(M2 - uman)**2)))

! antisymmetric octet -------------- 2 Re(Mt^((8A)prop) MAOctet^\ast)
      stst2gg_M2PropagatorCounter(3) = stst2gg_M2PropagatorCounter(3) +  2*dble((-(dZSf(isf1,isf2,itt1,igen1)-dMSf2(isf1,itt1,igen1)/DenSft(isf1,itt1,igen1)))*
     -    ( -(CF*CV*gs**4*Nc*(8*M2**4*((6 + FPSk1 + FPSk2)*sman - (FPSk1 + FPSk2 + FPSk1*FPSk2)*(tman - uman)) - 
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
      stst2gg_M2PropagatorCounter(3) = stst2gg_M2PropagatorCounter(3) +2*dble(  (-(dZSf(isf1,isf2,itt1,igen1)-dMSf2(isf1,itt1,igen1)/DenSfu(isf1,itt1,igen1)))*
     -    ( -(CF*CV*gs**4*Nc*(8*M2**4*((6 + FPSk1 + FPSk2)*sman + (FPSk1 + FPSk2 + FPSk1*FPSk2)*(tman - uman)) - 
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
      

      subroutine stst2gg_calc_M2Vertex4Counter(stst2gg_M2Vertex4Counter)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCounterterms.h"
#include "../common/stsT2xx_GenCouplings.h"

      double precision M2,Krondelta,stst2gg_M2Vertex4Counter(3)

      M2=MSf2(iSf1,itt1,igen1)
      
! singlet 
      stst2gg_M2Vertex4Counter(1) = 2*dble(-(gs**2*(-1 + Nc**2)*(FPSk1*(-1 + FPSk2)*(2*M2 - tman - uman)*
     -        (8*M2**2 + sman**2 + 4*sman*(tman + uman) - 8*M2*(sman + tman + uman) + 4*(tman**2 + uman**2)) - 
     -       FPSk2*(2*M2 - tman - uman)*(8*M2**2 + sman**2 + 4*sman*(tman + uman) - 8*M2*(sman + tman + uman) + 4*(tman**2 + uman**2)) + 
     -       sman*(-40*M2**2 + 10*M2*sman + 32*M2*(tman + uman) - 5*sman*(tman + uman) - 4*(tman**2 + 4*tman*uman + uman**2)))*
     -     dconjg(dSfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)))/(12*sman*(M2 - tman)*(M2 - uman)))


! symmetric octet
      stst2gg_M2Vertex4Counter(2) = 2*dble(-(CF*gs**2*(-4 + Nc**2)*(FPSk1*(-1 + FPSk2)*(2*M2 - tman - uman)*
     -        (8*M2**2 + sman**2 + 4*sman*(tman + uman) - 8*M2*(sman + tman + uman) + 4*(tman**2 + uman**2)) - 
     -       FPSk2*(2*M2 - tman - uman)*(8*M2**2 + sman**2 + 4*sman*(tman + uman) - 8*M2*(sman + tman + uman) + 4*(tman**2 + uman**2)) + 
     -       sman*(-40*M2**2 + 10*M2*sman + 32*M2*(tman + uman) - 5*sman*(tman + uman) - 4*(tman**2 + 4*tman*uman + uman**2)))*
     -     dconjg(dSfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)))/(4*sman*(M2 - tman)*(M2 - uman)))

! antisymmetric octet
      stst2gg_M2Vertex4Counter(3) = 0d0
           
      end


      subroutine stst2gg_calc_M2Vertex3gCounter(stst2gg_M2Vertex3gCounter)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCounterterms.h"
#include "../common/stsT2xx_GenCouplings.h"

      double precision M2,Krondelta,stst2gg_M2Vertex3gCounter(3)

      M2=MSf2(iSf1,itt1,igen1)
  
! only antisymmetric part contributes
      stst2gg_M2Vertex3gCounter = 0d0
          
      stst2gg_M2Vertex3gCounter(3) = 2*dble(-(CF*CV*gs**2*Nc*(16*(-2 + FPSk1 + FPSk2)*M2**3*sman - 2*(-2 + FPSk1 + FPSk2)*sman**3*(tman + uman) + 
     -       2*(-FPSk2 + FPSk1*(-1 + 2*FPSk2))*sman*(tman - uman)**2*(tman + uman) + 
     -       sman**2*((5 + FPSk1*(-5 + FPSk2) - 5*FPSk2)*tman**2 - 2*(1 + FPSk1*(-3 + FPSk2) - 3*FPSk2)*tman*uman + 
     -          (5 + FPSk1*(-5 + FPSk2) - 5*FPSk2)*uman**2) - 2*M2**2*
     -        (10*(-2 + FPSk1 + FPSk2)*sman**2 - (5 + FPSk1 + FPSk2 + FPSk1*FPSk2)*(tman - uman)**2 + 8*(-2 + FPSk1 + FPSk2)*sman*(tman + uman)) + 
     -       2*(tman - uman)**2*(FPSk1*tman*uman + (5 + FPSk2)*tman*uman + FPSk1*FPSk2*(2*tman**2 - 3*tman*uman + 2*uman**2)) + 
     -       2*M2*(2*(-2 + FPSk1 + FPSk2)*sman**3 + 6*(-2 + FPSk1 + FPSk2)*sman**2*(tman + uman) - 
     -          (5 + FPSk1 + FPSk2 + FPSk1*FPSk2)*(tman - uman)**2*(tman + uman) - 
     -          2*sman*((6 - 5*FPSk2 + FPSk1*(-5 + 2*FPSk2))*tman**2 + 2*(-2 + FPSk1*(3 - 2*FPSk2) + 3*FPSk2)*tman*uman + 
     -             (6 - 5*FPSk2 + FPSk1*(-5 + 2*FPSk2))*uman**2)))*dconjg(dggg)*dconjg(SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)))/
     -  (4*sman*(M2 - tman)*(M2 - uman)*sman))
     
      end


      subroutine stst2gg_calc_M2Vertex2sgCounter(stst2gg_M2Vertex2sgCounter)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCounterterms.h"
#include "../common/stsT2xx_GenCouplings.h"

      double complex DenS,DenSft(2,4,3), DenSfu(2,4,3)
      double precision M2,Krondelta,stst2gg_M2Vertex2sgCounter(3)
      integer i,j,k

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  s-channel correction !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
! antisymmetric part
      stst2gg_M2Vertex2sgCounter(3) = 2*dble(-(CF*CV*gs**2*GS*Nc*(16*(-2 + FPSk1 + FPSk2)*M2**3*sman - 2*(-2 + FPSk1 + FPSk2)*sman**3*(tman + uman) + 
     -       2*(-FPSk2 + FPSk1*(-1 + 2*FPSk2))*sman*(tman - uman)**2*(tman + uman) + 
     -       sman**2*((5 + FPSk1*(-5 + FPSk2) - 5*FPSk2)*tman**2 - 2*(1 + FPSk1*(-3 + FPSk2) - 3*FPSk2)*tman*uman + 
     -          (5 + FPSk1*(-5 + FPSk2) - 5*FPSk2)*uman**2) - 2*M2**2*
     -        (10*(-2 + FPSk1 + FPSk2)*sman**2 - (5 + FPSk1 + FPSk2 + FPSk1*FPSk2)*(tman - uman)**2 + 8*(-2 + FPSk1 + FPSk2)*sman*(tman + uman)) + 
     -       2*(tman - uman)**2*(FPSk1*tman*uman + (5 + FPSk2)*tman*uman + FPSk1*FPSk2*(2*tman**2 - 3*tman*uman + 2*uman**2)) + 
     -       2*M2*(2*(-2 + FPSk1 + FPSk2)*sman**3 + 6*(-2 + FPSk1 + FPSk2)*sman**2*(tman + uman) - 
     -          (5 + FPSk1 + FPSk2 + FPSk1*FPSk2)*(tman - uman)**2*(tman + uman) - 
     -          2*sman*((6 - 5*FPSk2 + FPSk1*(-5 + 2*FPSk2))*tman**2 + 2*(-2 + FPSk1*(3 - 2*FPSk2) + 3*FPSk2)*tman*uman + 
     -             (6 - 5*FPSk2 + FPSk1*(-5 + 2*FPSk2))*uman**2)))*dconjg(dSfSfg(isf2,itt2,igen2,isf1,itt1,igen1)))/
     -  (4*sman*(M2 - tman)*(M2 - uman)*dconjg(DenS)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Ghosts !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Ghost 1
      stst2gg_M2Vertex2sgCounter(3) = stst2gg_M2Vertex2sgCounter(3)-(1d0-FPSk1)*(1d0-FPSk2)*2*dble((-CF*CV*GS**3*Nc*(tman - uman)**2*dconjg(dSfSfg(isf1,itt1,igen1,isf2,itt2,igen2))*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/
     -  (4*DenS*dconjg(DenS)))

! Ghost 2
      stst2gg_M2Vertex2sgCounter(3) = stst2gg_M2Vertex2sgCounter(3)-(1d0-FPSk1)*(1d0-FPSk2)*2*dble((-CF*CV*GS**3*Nc*(tman - uman)**2*dconjg(dSfSfg(isf1,itt1,igen1,isf2,itt2,igen2))*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/
     -  (4*DenS*dconjg(DenS)))
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  t-channel correction upper+lower vertex !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Singlet      
      stst2gg_M2Vertex2sgCounter(1) = 4*dble((gs**2*GS*(-1 + Nc**2)*(12*M2 - 5*sman + (4*FPSk1*(M2 - tman)**2)/sman + (4*FPSk2*(M2 - tman)**2)/sman - (4*FPSk1*FPSk2*(M2 - tman)**2)/sman + 
     -      (4*(M2 + tman)**2)/(M2 - tman) + (FPSk1*FPSk2*(M2 - tman)**2*(-2*M2 + sman + 2*tman)**2)/(sman**2*(M2 - uman)) - 4*uman + 
     -      (FPSk1*(M2 - tman)*(-2*M2 + sman + 2*tman)*(-6*M2 + 2*sman + tman + uman))/(sman*(M2 - uman)) + 
     -      (FPSk2*(M2 - tman)*(-2*M2 + sman + 2*tman)*(-6*M2 + 2*sman + tman + uman))/(sman*(M2 - uman)) + (-6*M2 + 2*sman + tman + uman)**2/(M2 - uman) + 
     -      (4*FPSk1*(M2 + tman)*(-2*M2 + sman + 2*uman))/sman + (4*FPSk2*(M2 + tman)*(-2*M2 + sman + 2*uman))/sman + 
     -      (2*FPSk1*FPSk2*(M2 - tman)*(-2*M2 + sman + 2*tman)*(-2*M2 + sman + 2*uman))/sman**2 + 
     -      (FPSk1*(-6*M2 + 2*sman + tman + uman)*(-2*M2 + sman + 2*uman))/sman + (FPSk2*(-6*M2 + 2*sman + tman + uman)*(-2*M2 + sman + 2*uman))/sman + 
     -      (FPSk1*(-2*M2 + sman + 2*uman)**2)/sman + (FPSk2*(-2*M2 + sman + 2*uman)**2)/sman - (FPSk1*FPSk2*(-2*M2 + sman + 2*uman)**2)/sman + 
     -      (4*FPSk1*FPSk2*(M2 - tman)*(-2*M2 + sman + 2*uman)**2)/sman**2 + (FPSk1*FPSk2*(M2 - uman)*(-2*M2 + sman + 2*uman)**2)/sman**2)*
     -    dconjg(dSfSfg(isf2,itt2,igen2,isf1,itt1,igen1)))/(12*dconjg(DenSft(isf2,itt2,igen2))))


! symmetric octet     
      stst2gg_M2Vertex2sgCounter(2) = 4*dble((CF*gs**2*GS*(-4 + Nc**2)*(12*M2 - 5*sman + (4*FPSk1*(M2 - tman)**2)/sman + (4*FPSk2*(M2 - tman)**2)/sman - (4*FPSk1*FPSk2*(M2 - tman)**2)/sman + 
     -      (4*(M2 + tman)**2)/(M2 - tman) + (FPSk1*FPSk2*(M2 - tman)**2*(-2*M2 + sman + 2*tman)**2)/(sman**2*(M2 - uman)) - 4*uman + 
     -      (FPSk1*(M2 - tman)*(-2*M2 + sman + 2*tman)*(-6*M2 + 2*sman + tman + uman))/(sman*(M2 - uman)) + 
     -      (FPSk2*(M2 - tman)*(-2*M2 + sman + 2*tman)*(-6*M2 + 2*sman + tman + uman))/(sman*(M2 - uman)) + (-6*M2 + 2*sman + tman + uman)**2/(M2 - uman) + 
     -      (4*FPSk1*(M2 + tman)*(-2*M2 + sman + 2*uman))/sman + (4*FPSk2*(M2 + tman)*(-2*M2 + sman + 2*uman))/sman + 
     -      (2*FPSk1*FPSk2*(M2 - tman)*(-2*M2 + sman + 2*tman)*(-2*M2 + sman + 2*uman))/sman**2 + 
     -      (FPSk1*(-6*M2 + 2*sman + tman + uman)*(-2*M2 + sman + 2*uman))/sman + (FPSk2*(-6*M2 + 2*sman + tman + uman)*(-2*M2 + sman + 2*uman))/sman + 
     -      (FPSk1*(-2*M2 + sman + 2*uman)**2)/sman + (FPSk2*(-2*M2 + sman + 2*uman)**2)/sman - (FPSk1*FPSk2*(-2*M2 + sman + 2*uman)**2)/sman + 
     -      (4*FPSk1*FPSk2*(M2 - tman)*(-2*M2 + sman + 2*uman)**2)/sman**2 + (FPSk1*FPSk2*(M2 - uman)*(-2*M2 + sman + 2*uman)**2)/sman**2)*
     -    dconjg(dSfSfg(isf2,itt2,igen2,isf1,itt1,igen1)))/(4*dconjg(DenSft(isf2,itt2,igen2))))

! antisymmetric octet
      stst2gg_M2Vertex2sgCounter(3) = stst2gg_M2Vertex2sgCounter(3) + 4*dble(-(CF*CV*gs**2*GS*Nc*(8*M2**4*((6 + FPSk1 + FPSk2)*sman - (FPSk1 + FPSk2 + FPSk1*FPSk2)*(tman - uman)) - 
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
     -             FPSk1*((3 + 6*FPSk2)*tman**2 + (11 + 2*FPSk2)*tman*uman + 2*(1 + 4*FPSk2)*uman**2))))*dconjg(dSfSfg(isf2,itt2,igen2,isf1,itt1,igen1)))/
     -  (4*sman**2*(M2 - tman)*(M2 - uman)*dconjg(DenSft(isf2,itt2,igen2))))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   u-channel correction upper+lower vertex !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
! Singlet      
      stst2gg_M2Vertex2sgCounter(1) =  stst2gg_M2Vertex2sgCounter(1) + 4*dble((gs**3*(-1 + Nc**2)*(12*M2 - 5*sman - 4*tman + (FPSk1*(-2*M2 + sman + 2*tman)**2)/sman + (FPSk2*(-2*M2 + sman + 2*tman)**2)/sman - 
     -      (FPSk1*FPSk2*(-2*M2 + sman + 2*tman)**2)/sman + (FPSk1*FPSk2*(M2 - tman)*(-2*M2 + sman + 2*tman)**2)/sman**2 + 
     -      (4*FPSk1*FPSk2*(-2*M2 + sman + 2*tman)**2*(M2 - uman))/sman**2 + (4*FPSk1*(M2 - uman)**2)/sman + (4*FPSk2*(M2 - uman)**2)/sman - 
     -      (4*FPSk1*FPSk2*(M2 - uman)**2)/sman + (4*FPSk1*(-2*M2 + sman + 2*tman)*(M2 + uman))/sman + (4*FPSk2*(-2*M2 + sman + 2*tman)*(M2 + uman))/sman + 
     -      (4*(M2 + uman)**2)/(M2 - uman) + (FPSk1*(-2*M2 + sman + 2*tman)*(-6*M2 + 2*sman + tman + uman))/sman + 
     -      (FPSk2*(-2*M2 + sman + 2*tman)*(-6*M2 + 2*sman + tman + uman))/sman + (-6*M2 + 2*sman + tman + uman)**2/(M2 - tman) + 
     -      (2*FPSk1*FPSk2*(-2*M2 + sman + 2*tman)*(M2 - uman)*(-2*M2 + sman + 2*uman))/sman**2 + 
     -      (FPSk1*(M2 - uman)*(-6*M2 + 2*sman + tman + uman)*(-2*M2 + sman + 2*uman))/(sman*(M2 - tman)) + 
     -      (FPSk2*(M2 - uman)*(-6*M2 + 2*sman + tman + uman)*(-2*M2 + sman + 2*uman))/(sman*(M2 - tman)) + 
     -      (FPSk1*FPSk2*(M2 - uman)**2*(-2*M2 + sman + 2*uman)**2)/(sman**2*(M2 - tman)))*dconjg(dSfSfg(isf1,itt1,igen1,isf2,itt2,igen2)))/
     -  (12*dconjg(DenSfu(isf2,itt2,igen2))))

! symmetric octet     
      stst2gg_M2Vertex2sgCounter(2) =  stst2gg_M2Vertex2sgCounter(2) + 4*dble((CF*gs**3*(-4 + Nc**2)*(12*M2 - 5*sman - 4*tman + (FPSk1*(-2*M2 + sman + 2*tman)**2)/sman + (FPSk2*(-2*M2 + sman + 2*tman)**2)/sman - 
     -      (FPSk1*FPSk2*(-2*M2 + sman + 2*tman)**2)/sman + (FPSk1*FPSk2*(M2 - tman)*(-2*M2 + sman + 2*tman)**2)/sman**2 + 
     -      (4*FPSk1*FPSk2*(-2*M2 + sman + 2*tman)**2*(M2 - uman))/sman**2 + (4*FPSk1*(M2 - uman)**2)/sman + (4*FPSk2*(M2 - uman)**2)/sman - 
     -      (4*FPSk1*FPSk2*(M2 - uman)**2)/sman + (4*FPSk1*(-2*M2 + sman + 2*tman)*(M2 + uman))/sman + (4*FPSk2*(-2*M2 + sman + 2*tman)*(M2 + uman))/sman + 
     -      (4*(M2 + uman)**2)/(M2 - uman) + (FPSk1*(-2*M2 + sman + 2*tman)*(-6*M2 + 2*sman + tman + uman))/sman + 
     -      (FPSk2*(-2*M2 + sman + 2*tman)*(-6*M2 + 2*sman + tman + uman))/sman + (-6*M2 + 2*sman + tman + uman)**2/(M2 - tman) + 
     -      (2*FPSk1*FPSk2*(-2*M2 + sman + 2*tman)*(M2 - uman)*(-2*M2 + sman + 2*uman))/sman**2 + 
     -      (FPSk1*(M2 - uman)*(-6*M2 + 2*sman + tman + uman)*(-2*M2 + sman + 2*uman))/(sman*(M2 - tman)) + 
     -      (FPSk2*(M2 - uman)*(-6*M2 + 2*sman + tman + uman)*(-2*M2 + sman + 2*uman))/(sman*(M2 - tman)) + 
     -      (FPSk1*FPSk2*(M2 - uman)**2*(-2*M2 + sman + 2*uman)**2)/(sman**2*(M2 - tman)))*dconjg(dSfSfg(isf1,itt1,igen1,isf2,itt2,igen2)))/
     -  (4*dconjg(DenSfu(isf2,itt2,igen2))))

! antisymmetric octet
      stst2gg_M2Vertex2sgCounter(3) = stst2gg_M2Vertex2sgCounter(3) + 4*dble(-(CF*CV*gs**3*Nc*(8*M2**4*((6 + FPSk1 + FPSk2)*sman + (FPSk1 + FPSk2 + FPSk1*FPSk2)*(tman - uman)) - 
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
     -             FPSk1*((2 + 8*FPSk2)*tman**2 + (11 + 2*FPSk2)*tman*uman + 3*(1 + 2*FPSk2)*uman**2))))*dconjg(dSfSfg(isf1,itt1,igen1,isf2,itt2,igen2)))/
     -  (4*sman**2*(M2 - tman)*(M2 - uman)*dconjg(DenSfu(isf2,itt2,igen2))))

      end

      subroutine stst2gg_calc_M2VertexGhCounter(M2Color)
      implicit none
#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../common/stsT2xx_GenCounterterms.h"
#include "../../util/DMNLO_RenScheme.h"

      double precision Krondelta,M2Color(3)

      M2Color = 0d0

! Interfernce of MGh1 Counter with MGh1Tree
      M2Color(3) = 2*dble((CF*CV*GS**3*Nc*(-tman + uman)**2*dconjg(dGhGhG)*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(4*sman**2))
  
! Interfernce of MGh2 Counter with MGh2Tree
      M2Color(3) = M2Color(3) + 2*dble((CF*CV*GS**3*Nc*(tman - uman)**2*dconjg(dGhGhG)*Krondelta(igen1,igen2)*Krondelta(isf1,isf2)*Krondelta(itt1,itt2))/(4*sman**2))

      end 