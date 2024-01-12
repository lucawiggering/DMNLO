************************************************

	  function stst2QQ_M2CountertermsProp()
	  implicit none

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
#include "stst2QQ_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"
#include "../util/DMNLO_RenScheme.h"

	  integer i,j,k,l,n,m,tcount,itypeU,itypeT,att,antitt,b,neu
	  complex*16 g0,g0L,g0R,g1,g1L,g1R,g2,g2L,g2R,g3L,g3R,g4L,g4R
	  complex*16 g1Lc, g1Rc, g2Lc, g2Rc, g3Rc, g3Lc, g4Lc, g4Rc
	  double precision Ml0,Ml1,tempCounter, tempH, tempZ, tempt, tempu
	  double precision tempH11, tempZ11, tempt11, tempu11, tempH12, tempZ12, tempt12, tempu12
	  double precision tempH21, tempZ21, tempt21, tempu21, tempH22, tempZ22, tempt22, tempu22
	  complex*16 dPropSGlL,dPropSGlR,dPropFGlL,dPropFGlR,dPropSG,dPropFG
	  double precision swap
	  double precision Mloop0,Mloop1,Mloop2

	  double precision Krondelta, ef, I3f, hf

	  double complex test, Gscalar

	  complex*16 BL1g,BL2g,BL3g,BR1g,BR2g,BR3g

	  double complex nc0,nc1,nc2,nc00,nc11,nc12,nc22

	  complex*16 DenSmV, DenSmH(6)
	  complex*16 DenSmV1, DenSmH1(6)
	  complex*16 tempGluonL, tempGluonR

* functions & subroutines
	  double precision colfac
	  complex*16 stst2QQ_M2CountertermsProp
	  complex*16 B02,C03f
	  external myCget, myBget

	  stst2QQ_M2CountertermsProp = 0d0
	  test=0d0


c 	Propagator--------------------------------------------------------------------------------------

	  do i=1,4 !hier gluino und neutralino unterscheiden
		DenGauT(i)=(tman - MNeu(i)**2)! + cI*MNeu(i)*wNeu(i)
	  enddo
	  DenGauT(5)=(tman - MCha(1)**2)! + cI*MCha(1)*wCha(1)
	  DenGauT(6)=(tman - MCha(2)**2)! + cI*MCha(2)*wCha(2)
	
	  do i=1,4 !hier gluino und neutralino unterscheiden
		DenGauU(i)=(uman - MNeu(i)**2)! + cI*MNeu(i)*wNeu(i)
	  enddo
	  DenGauU(5)=(uman - MCha(1)**2)! + cI*MCha(1)*wCha(1)
	  DenGauU(6)=(uman - MCha(2)**2)! + cI*MCha(2)*wCha(2)


	  DenGlT = (tman-MGl**2)
	  DenGlU = (uman-MGl**2)
	
! 	  only for debugging:
! 	   dZGlR = 0d0
! 	   dZGlL = 0d0
	
!CountertermKop Gluino:
	  dPropFGlL = 1/2d0*(dZGlL+dconjg(dZGlL))
	
	  dPropFGlR = 1/2d0*(dZGlR+dconjg(dZGlR))
	
	  dPropSGlL = MGl/2d0*(dZGlL+dconjg(dZGlR))+dMGl
	
	  dPropSGlR = MGl/2d0*(dZGlR+dconjg(dZGlL))+dMGl

	
	
! !======================================================================================================================================
! ! ===============================================Gluino-Propagator Correction==========================================================
! !======================================================================================================================================

! ! 
! ==========================================================t-channel counter=============================================================
! -----------------------Amplitude square for T-channel with Gluino -----------
	  colfac = Nc/2d0*CF
	
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(g1L)
	  g1Rc = dconjg(g1R)

	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(g2L)
	  g2Rc = dconjg(g2R)
	  
	  colfac = Nc/2d0*CF
	
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(GluLt1(isf1,itt1,ftt1,igen1,fgen1))
	  g1Rc = dconjg(GluRt1(isf1,itt1,ftt1,igen1,fgen1))
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(g2L)
	  g2Rc = dconjg(g2R)

	  stst2QQ_M2CountertermsProp = stst2QQ_M2CountertermsProp + 2d0*colfac*dreal((sqrtS**5*(-(dPropSGlL*
     -         (g1R*(g1Rc*g2R*tred*
     -               (g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -                 g2Lc*muu4*(-muu1**2 + muu3**2 + tred)) + 
     -              g1Lc*g2R*muu3*tred*
     -               (2*g2Lc*muu4*muugl + 
     -                 g2Rc*(-muu2**2 + muu4**2 + tred)) + 
     -              g1Lc*g2L*muu3*muugl*
     -               (2*g2Rc*muu4*tred + 
     -                 g2Lc*muugl*(-muu2**2 + muu4**2 + tred)) + 
     -              g1Rc*g2L*muugl*
     -               (g2Rc*muu4*muugl*(-muu1**2 + muu3**2 + tred) + 
     -                 g2Lc*((muu1 - muu3)*(muu1 + muu3)*
     -                     (muu2 - muu4)*(muu2 + muu4) + tred - 
     -                    (muu1**2 + muu2**2)*tred + tred**2))) + 
     -           g1L*muugl*(g1Lc*g2L*muugl*
     -               (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -                 g2Rc*muu4*(-muu1**2 + muu3**2 + tred)) + 
     -              g1Rc*g2L*muu3*muugl*
     -               (2*g2Rc*muu4*muugl + 
     -                 g2Lc*(-muu2**2 + muu4**2 + tred)) + 
     -              g1Rc*g2R*muu3*
     -               (2*g2Lc*muu4*tred + 
     -                 g2Rc*muugl*(-muu2**2 + muu4**2 + tred)) + 
     -              g1Lc*g2R*
     -               (g2Lc*muu4*muugl*(-muu1**2 + muu3**2 + tred) + 
     -                 g2Rc*((muu1 - muu3)*(muu1 + muu3)*
     -                     (muu2 - muu4)*(muu2 + muu4) + tred - 
     -                    (muu1**2 + muu2**2)*tred + tred**2))))) - 
     -      dPropSGlR*(g1R*muugl*
     -          (g1Rc*g2R*muugl*
     -             (g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2Lc*muu4*(-muu1**2 + muu3**2 + tred)) + 
     -            g1Lc*g2R*muu3*muugl*
     -             (2*g2Lc*muu4*muugl + 
     -               g2Rc*(-muu2**2 + muu4**2 + tred)) + 
     -            g1Lc*g2L*muu3*
     -             (2*g2Rc*muu4*tred + 
     -               g2Lc*muugl*(-muu2**2 + muu4**2 + tred)) + 
     -            g1Rc*g2L*(g2Rc*muu4*muugl*
     -                (-muu1**2 + muu3**2 + tred) + 
     -               g2Lc*((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                   (muu2 + muu4) + tred - 
     -                  (muu1**2 + muu2**2)*tred + tred**2))) + 
     -         g1L*(g1Lc*g2L*tred*
     -             (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2Rc*muu4*(-muu1**2 + muu3**2 + tred)) + 
     -            g1Rc*g2L*muu3*tred*
     -             (2*g2Rc*muu4*muugl + 
     -               g2Lc*(-muu2**2 + muu4**2 + tred)) + 
     -            g1Rc*g2R*muu3*muugl*
     -             (2*g2Lc*muu4*tred + 
     -               g2Rc*muugl*(-muu2**2 + muu4**2 + tred)) + 
     -            g1Lc*g2R*muugl*
     -             (g2Lc*muu4*muugl*(-muu1**2 + muu3**2 + tred) + 
     -               g2Rc*((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                   (muu2 + muu4) + tred - 
     -                  (muu1**2 + muu2**2)*tred + tred**2)))) + 
     -      sqrtS*(dPropFGlR*g1R*tred*
     -          (g1Rc*g2R*muugl*
     -             (g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2Lc*muu4*(-muu1**2 + muu3**2 + tred)) + 
     -            g1Lc*g2R*muu3*muugl*
     -             (2*g2Lc*muu4*muugl + 
     -               g2Rc*(-muu2**2 + muu4**2 + tred)) + 
     -            g1Lc*g2L*muu3*
     -             (2*g2Rc*muu4*tred + 
     -               g2Lc*muugl*(-muu2**2 + muu4**2 + tred)) + 
     -            g1Rc*g2L*(g2Rc*muu4*muugl*
     -                (-muu1**2 + muu3**2 + tred) + 
     -               g2Lc*((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                   (muu2 + muu4) + tred - 
     -                  (muu1**2 + muu2**2)*tred + tred**2))) + 
     -         dPropFGlL*g1R*muugl*
     -          (g1Rc*g2R*tred*
     -             (g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2Lc*muu4*(-muu1**2 + muu3**2 + tred)) + 
     -            g1Lc*g2R*muu3*tred*
     -             (2*g2Lc*muu4*muugl + 
     -               g2Rc*(-muu2**2 + muu4**2 + tred)) + 
     -            g1Lc*g2L*muu3*muugl*
     -             (2*g2Rc*muu4*tred + 
     -               g2Lc*muugl*(-muu2**2 + muu4**2 + tred)) + 
     -            g1Rc*g2L*muugl*
     -             (g2Rc*muu4*muugl*(-muu1**2 + muu3**2 + tred) + 
     -               g2Lc*((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                   (muu2 + muu4) + tred - 
     -                  (muu1**2 + muu2**2)*tred + tred**2))) + 
     -         dPropFGlL*g1L*tred*
     -          (g1Lc*g2L*muugl*
     -             (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2Rc*muu4*(-muu1**2 + muu3**2 + tred)) + 
     -            g1Rc*g2L*muu3*muugl*
     -             (2*g2Rc*muu4*muugl + 
     -               g2Lc*(-muu2**2 + muu4**2 + tred)) + 
     -            g1Rc*g2R*muu3*
     -             (2*g2Lc*muu4*tred + 
     -               g2Rc*muugl*(-muu2**2 + muu4**2 + tred)) + 
     -            g1Lc*g2R*(g2Lc*muu4*muugl*
     -                (-muu1**2 + muu3**2 + tred) + 
     -               g2Rc*((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                   (muu2 + muu4) + tred - 
     -                  (muu1**2 + muu2**2)*tred + tred**2))) + 
     -         dPropFGlR*g1L*muugl*
     -          (g1Lc*g2L*tred*
     -             (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2Rc*muu4*(-muu1**2 + muu3**2 + tred)) + 
     -            g1Rc*g2L*muu3*tred*
     -             (2*g2Rc*muu4*muugl + 
     -               g2Lc*(-muu2**2 + muu4**2 + tred)) + 
     -            g1Rc*g2R*muu3*muugl*
     -             (2*g2Lc*muu4*tred + 
     -               g2Rc*muugl*(-muu2**2 + muu4**2 + tred)) + 
     -            g1Lc*g2R*muugl*
     -             (g2Lc*muu4*muugl*(-muu1**2 + muu3**2 + tred) + 
     -               g2Rc*((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                   (muu2 + muu4) + tred - 
     -                  (muu1**2 + muu2**2)*tred + tred**2))))))/
     -  (DenGlT**2*dconjg(DenGlT)))

! !      
! ! ! ---------------------t-channel Count and u-channel Gluino squared------------------------------
! ! 
! ! 
	  colfac = -1d0/2d0*CF
	
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(GluLt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Rc = dconjg(GluRt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(GluLt1(isf1,itt1,ftt2,igen1,fgen2))
	
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Rc = dconjg(GluRt1(isf1,itt1,ftt2,igen1,fgen2))
	  
	  stst2QQ_M2CountertermsProp = stst2QQ_M2CountertermsProp + 2d0*colfac*
     -  	  dreal((sqrtS**5*(-(dPropSGlR*g1R*muugl*
     -         (g1Rc*g2R*muugl*
     -            (g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -              g2Lc*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)) + 
     -           g1Lc*g2R*muu3*muugl*
     -            (2*g2Lc*muu4*muugl + 
     -              g2Rc*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) +
     -            g1Lc*g2L*muu3*
     -            (g2Rc*muu4*
     -               (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -              g2Lc*muugl*(-muu2**2 + muu4**2 + tred)) + 
     -           g1Rc*g2L*(g2Rc*muu4*muugl*
     -               (-muu1**2 + muu3**2 + tred) + 
     -              g2Lc*((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -                 muu1**2*(-muu2**2 + tred) + 
     -                 muu3**2*(muu4**2 + tred))))) + 
     -      dPropFGlR*g1R*sqrtS*tred*
     -       (g1Rc*g2R*muugl*
     -          (g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2Lc*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -               tred)) + 
     -         g1Lc*g2R*muu3*muugl*
     -          (2*g2Lc*muu4*muugl + 
     -            g2Rc*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -         g1Lc*g2L*muu3*
     -          (g2Rc*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -            g2Lc*muugl*(-muu2**2 + muu4**2 + tred)) + 
     -         g1Rc*g2L*(g2Rc*muu4*muugl*
     -             (-muu1**2 + muu3**2 + tred) + 
     -            g2Lc*((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred)))) + 
     -      dPropFGlL*g1R*muugl*sqrtS*
     -       (g1Rc*g2R*(g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2Lc*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -               tred))*tred + 
     -         g1Lc*g2R*muu3*
     -          (2*g2Lc*muu4*muugl + 
     -            g2Rc*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -          tred + g1Lc*g2L*muu3*muugl*
     -          (g2Rc*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -            g2Lc*muugl*(-muu2**2 + muu4**2 + tred)) + 
     -         g1Rc*g2L*muugl*
     -          (g2Rc*muu4*muugl*(-muu1**2 + muu3**2 + tred) + 
     -            g2Lc*((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred)))) - 
     -      dPropSGlL*g1L*muugl*
     -       (g1Lc*g2L*muugl*
     -          (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2Rc*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -               tred)) + 
     -         g1Rc*g2Lc*muu3*
     -          (g2R*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -            g2L*muugl*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -               tred)) + 
     -         g1Rc*g2Rc*muu3*muugl*
     -          (2*g2L*muu4*muugl + g2R*(-muu2**2 + muu4**2 + tred)) +
     -           g1Lc*g2R*(g2Lc*muu4*muugl*
     -             (-muu1**2 + muu3**2 + tred) + 
     -            g2Rc*((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred)))) + 
     -      dPropFGlL*g1L*sqrtS*tred*
     -       (g1Lc*g2L*muugl*
     -          (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2Rc*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -               tred)) + 
     -         g1Rc*g2Lc*muu3*
     -          (g2R*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -            g2L*muugl*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -               tred)) + 
     -         g1Rc*g2Rc*muu3*muugl*
     -          (2*g2L*muu4*muugl + g2R*(-muu2**2 + muu4**2 + tred)) +
     -          g1Lc*g2R*(g2Lc*muu4*muugl*
     -             (-muu1**2 + muu3**2 + tred) + 
     -            g2Rc*((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred)))) + 
     -      dPropFGlR*g1L*muugl*sqrtS*
     -       (g1Lc*g2L*(g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2Rc*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -               tred))*tred + 
     -         g1Rc*g2Lc*muu3*
     -          (g2R*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2)*
     -             muugl + g2L*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*tred) +
     -          g1Rc*g2Rc*muu3*muugl*
     -          (2*g2L*muu4*tred + 
     -            g2R*muugl*(-muu2**2 + muu4**2 + tred)) + 
     -         g1Lc*g2R*muugl*
     -          (g2Lc*muu4*muugl*(-muu1**2 + muu3**2 + tred) + 
     -            g2Rc*((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred)))) + 
     -      dPropSGlL*g1R*(-(g1Lc*muu3*
     -            (g2R*(2*g2Lc*muu4*muugl + 
     -                 g2Rc*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                    tred))*tred + 
     -              g2L*muugl*
     -               (g2Rc*muu4*
     -                  (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -                 g2Lc*muugl*(-muu2**2 + muu4**2 + tred)))) + 
     -         g1Rc*(g2R*(-(g2Rc*(-1 + muu3**2 + muu4**2)*muugl) - 
     -               g2Lc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred))*
     -             tred + g2L*muugl*
     -             (g2Rc*muu4*muugl*(muu1**2 - muu3**2 - tred) + 
     -               g2Lc*(muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -                  (muu1**2 + muu2**2 + muu3**2 + muu4**2)*
     -                   tred + tred**2)))) + 
     -      dPropSGlR*g1L*(-(g1Rc*muu3*
     -            (g2Lc*(g2R*muu4*
     -                  (-muu1**2 - muu2**2 + muu3**2 + muu4**2)*muugl
     -                   + g2L*
     -                  (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)*
     -                  tred) + 
     -              g2Rc*muugl*
     -               (2*g2L*muu4*tred + 
     -                 g2R*muugl*(-muu2**2 + muu4**2 + tred)))) + 
     -         g1Lc*(g2L*(-(g2Lc*(-1 + muu3**2 + muu4**2)*muugl) - 
     -               g2Rc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred))*
     -             tred + g2R*muugl*
     -             (g2Lc*muu4*muugl*(muu1**2 - muu3**2 - tred) + 
     -               g2Rc*(muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -                  (muu1**2 + muu2**2 + muu3**2 + muu4**2)*
     -                   tred + tred**2))))))/
     -  (DenGlT**2*dconjg(DenGlU)))	  
! 	  
! 	  
! ! ! ! !---------------------t-channel Gluino t-channel Neutralino--------------------------------------------
      do n=1,6
	  colfac = 0d0
	  
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  
	  g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  
	  stst2QQ_M2CountertermsProp = stst2QQ_M2CountertermsProp + 2d0*colfac*
     -  	  dreal((sqrtS**5*(sqrtS*((dPropFGlL*g1R*g2L*g3Rc*g4Lc + 
     -            dPropFGlR*g1L*g2R*g3Lc*g4Rc)*(muu1 - muu3)*
     -          (muu1 + muu3)*(muu2 - muu4)*(muu2 + muu4)*muugl**2 + 
     -         ((dPropFGlR*g1R*g2L*g3Rc*g4Lc + 
     -               dPropFGlL*g1L*g2R*g3Lc*g4Rc)*(muu1 - muu3)*
     -             (muu1 + muu3)*(muu2 - muu4)*(muu2 + muu4) - 
     -            (dPropFGlL + dPropFGlR)*
     -             ((g1L*g2L*g3Rc*g4Lc + g1R*g2R*g3Lc*g4Rc)*muu2**2*
     -                muu3 + 
     -               (g1R*g2R*g3Rc*g4Lc + g1L*g2L*g3Lc*g4Rc)*
     -                (muu1 - muu3)*(muu1 + muu3)*muu4 - 
     -               (g1L*g2L*g3Rc*g4Lc + g1R*g2R*g3Lc*g4Rc)*muu3*
     -                muu4**2)*muugl - 
     -            ((dPropFGlL*g1R*g2L*g3Rc*g4Lc + 
     -                  dPropFGlR*g1L*g2R*g3Lc*g4Rc)*
     -                (-1 + muu1**2 + muu2**2) - 
     -               2*(dPropFGlR*g1L*g2R*g3Rc*g4Lc + 
     -                  dPropFGlL*g1R*g2L*g3Lc*g4Rc)*muu3*muu4)*
     -             muugl**2)*tred + 
     -         (dPropFGlR*g1R*g2L*g3Rc*g4Lc + 
     -            dPropFGlL*g1L*g2R*g3Lc*g4Rc - 
     -            dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2 - 
     -            dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2 - 
     -            dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**2 - 
     -            dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**2 + 
     -            2*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3*muu4 + 
     -            2*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3*muu4 + 
     -            (dPropFGlL + dPropFGlR)*
     -             (g1R*g2R*(g3Lc*g4Rc*muu3 + g3Rc*g4Lc*muu4) + 
     -               g1L*g2L*(g3Rc*g4Lc*muu3 + g3Lc*g4Rc*muu4))*muugl +
     -              (dPropFGlL*g1R*g2L*g3Rc*g4Lc + 
     -               dPropFGlR*g1L*g2R*g3Lc*g4Rc)*muugl**2)*tred**2 + 
     -         (dPropFGlR*g1R*g2L*g3Rc*g4Lc + 
     -            dPropFGlL*g1L*g2R*g3Lc*g4Rc)*tred**3) - 
     -      dPropSGlL*(g1L*muugl*
     -          (2*g2R*g3Rc*g4Lc*muu3*muu4*tred + 
     -            g2L*g3Lc*g4Rc*muu4*muugl*
     -             (-muu1**2 + muu3**2 + tred) + 
     -            g2L*g3Rc*g4Lc*muu3*muugl*
     -             (-muu2**2 + muu4**2 + tred) + 
     -            g2R*g3Lc*g4Rc*
     -             ((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                (muu2 + muu4) + tred - 
     -               (muu1**2 + muu2**2)*tred + tred**2)) + 
     -         g1R*(g2R*tred*
     -             (g3Rc*g4Lc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -               g3Lc*g4Rc*muu3*(-muu2**2 + muu4**2 + tred)) + 
     -            g2L*muugl*(2*g3Lc*g4Rc*muu3*muu4*tred + 
     -               g3Rc*g4Lc*
     -                ((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                   (muu2 + muu4) + tred - 
     -                  (muu1**2 + muu2**2)*tred + tred**2)))) - 
     -      dPropSGlR*(g1R*muugl*
     -          (2*g2L*g3Lc*g4Rc*muu3*muu4*tred + 
     -            g2R*g3Rc*g4Lc*muu4*muugl*
     -             (-muu1**2 + muu3**2 + tred) + 
     -            g2R*g3Lc*g4Rc*muu3*muugl*
     -             (-muu2**2 + muu4**2 + tred) + 
     -            g2L*g3Rc*g4Lc*
     -             ((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                (muu2 + muu4) + tred - 
     -               (muu1**2 + muu2**2)*tred + tred**2)) + 
     -         g1L*(g2L*tred*
     -             (g3Lc*g4Rc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -               g3Rc*g4Lc*muu3*(-muu2**2 + muu4**2 + tred)) + 
     -            g2R*muugl*(2*g3Rc*g4Lc*muu3*muu4*tred + 
     -               g3Lc*g4Rc*
     -                ((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                   (muu2 + muu4) + tred - 
     -                  (muu1**2 + muu2**2)*tred + tred**2)))) + 
     -      (dPropSGlR*muugl*
     -          (g1R*g2L*(g3Rc*g4Rc*(muu1 - muu3)*(muu1 + muu3)*
     -                muu4 + 
     -               g3Lc*g4Lc*muu3*(muu2 - muu4)*(muu2 + muu4)) + 
     -            g1L*g2R*(g3Lc*g4Lc*(muu1 - muu3)*(muu1 + muu3)*
     -                muu4 + 
     -               g3Rc*g4Rc*muu3*(muu2 - muu4)*(muu2 + muu4)) - 
     -            g1R*g2R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -               g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*muugl) + 
     -         (-((dPropFGlL*g1R*g2L*g3Lc*g4Lc + 
     -                 dPropFGlR*g1L*g2R*g3Rc*g4Rc)*muu2**2*muu3) + 
     -            (dPropFGlR*g1L*g2R*g3Lc*g4Lc + 
     -               dPropFGlL*g1R*g2L*g3Rc*g4Rc)*(-muu1 + muu3)*
     -             (muu1 + muu3)*muu4 + 
     -            (dPropFGlL*g1R*g2L*g3Lc*g4Lc + 
     -               dPropFGlR*g1L*g2R*g3Rc*g4Rc)*muu3*muu4**2)*
     -          muugl**2*sqrtS - 
     -         dPropSGlR*(g1L*g2L*
     -             (2*g3Rc*g4Rc*muu3*muu4 + 
     -               g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)) + 
     -            g1L*g2R*(g3Rc*g4Rc*muu3 + g3Lc*g4Lc*muu4)*muugl + 
     -            g1R*g2L*(g3Lc*g4Lc*muu3 + g3Rc*g4Rc*muu4)*muugl)*
     -          tred + (-(dPropFGlR*g1R*g2L*g3Lc*g4Lc*muu2**2*muu3) - 
     -            dPropFGlL*g1L*g2R*g3Rc*g4Rc*muu2**2*muu3 - 
     -            dPropFGlL*g1L*g2R*g3Lc*g4Lc*muu1**2*muu4 - 
     -            dPropFGlR*g1R*g2L*g3Rc*g4Rc*muu1**2*muu4 + 
     -            dPropFGlL*g1L*g2R*g3Lc*g4Lc*muu3**2*muu4 + 
     -            dPropFGlR*g1R*g2L*g3Rc*g4Rc*muu3**2*muu4 + 
     -            dPropFGlR*g1R*g2L*g3Lc*g4Lc*muu3*muu4**2 + 
     -            dPropFGlL*g1L*g2R*g3Rc*g4Rc*muu3*muu4**2 + 
     -            (dPropFGlL + dPropFGlR)*
     -             (g1L*g2L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -                  g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)) + 
     -               g1R*g2R*
     -                (2*g3Lc*g4Lc*muu3*muu4 + 
     -                  g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)))*muugl + 
     -            (dPropFGlR*g1L*g2R*
     -                (g3Rc*g4Rc*muu3 + g3Lc*g4Lc*muu4) + 
     -               dPropFGlL*g1R*g2L*
     -                (g3Lc*g4Lc*muu3 + g3Rc*g4Rc*muu4))*muugl**2)*
     -          sqrtS*tred + 
     -         (dPropFGlL*g1L*g2R*(g3Rc*g4Rc*muu3 + g3Lc*g4Lc*muu4) + 
     -            dPropFGlR*g1R*g2L*(g3Lc*g4Lc*muu3 + g3Rc*g4Rc*muu4))
     -           *sqrtS*tred**2 + 
     -         dPropSGlL*(g1R*g2L*
     -             (g3Rc*g4Rc*(muu1 - muu3)*(muu1 + muu3)*muu4 + 
     -               g3Lc*g4Lc*muu3*(muu2 - muu4)*(muu2 + muu4))*muugl
     -              - g1R*(2*g2R*g3Lc*g4Lc*muu3*muu4 + 
     -               g2R*g3Rc*g4Rc*(-1 + muu3**2 + muu4**2) + 
     -               g2L*(g3Lc*g4Lc*muu3 + g3Rc*g4Rc*muu4)*muugl)*tred
     -              - g1L*muugl*
     -             (2*g2L*g3Rc*g4Rc*muu3*muu4*muugl + 
     -               g2L*g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2R*g3Lc*g4Lc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -               g2R*g3Rc*g4Rc*muu3*(-muu2**2 + muu4**2 + tred))))
     -        *muupt(n)))/(DenGlT**2*dconjg(DenGauT(n))))	  
	  
	  enddo
! !---------------------t-channel Gluino u-channel Neutralino--------------------------------------------------------
      do n=1,6
	  colfac = CF*Nc
	  
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
! 	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))	  
	  
	  stst2QQ_M2CountertermsProp = stst2QQ_M2CountertermsProp + 2d0*colfac*
     -  	  dreal((sqrtS**5*(sqrtS*(-(dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu1**2*
     -            muu2**2*muugl**2) - 
     -         dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu1**2*muu2**2*muugl**2 - 
     -         dPropFGlR*g1L*g2R*g3Rc*g4Lc*muu1**2*muu3*muu4*
     -          muugl**2 - dPropFGlL*g1R*g2L*g3Lc*g4Rc*muu1**2*muu3*
     -          muu4*muugl**2 - 
     -         dPropFGlR*g1L*g2R*g3Rc*g4Lc*muu2**2*muu3*muu4*
     -          muugl**2 - dPropFGlL*g1R*g2L*g3Lc*g4Rc*muu2**2*muu3*
     -          muu4*muugl**2 + 
     -         dPropFGlR*g1L*g2R*g3Rc*g4Lc*muu3**3*muu4*muugl**2 + 
     -         dPropFGlL*g1R*g2L*g3Lc*g4Rc*muu3**3*muu4*muugl**2 + 
     -         dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu3**2*muu4**2*muugl**2 + 
     -         dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu3**2*muu4**2*muugl**2 + 
     -         dPropFGlR*g1L*g2R*g3Rc*g4Lc*muu3*muu4**3*muugl**2 + 
     -         dPropFGlL*g1R*g2L*g3Lc*g4Rc*muu3*muu4**3*muugl**2 + 
     -         (-(dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu2**2) - 
     -            dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu2**2 - 
     -            dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu1**2*muu3*muu4 - 
     -            dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu1**2*muu3*muu4 - 
     -            dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu2**2*muu3*muu4 - 
     -            dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu2**2*muu3*muu4 + 
     -            dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3**3*muu4 + 
     -            dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3**3*muu4 + 
     -            dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**2*muu4**2 + 
     -            dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*muu4**2 + 
     -            dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3*muu4**3 + 
     -            dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3*muu4**3 + 
     -            (dPropFGlL + dPropFGlR)*
     -             ((g1L*g2L*g3Rc*g4Lc + g1R*g2R*g3Lc*g4Rc)*muu3*
     -                (-1 + muu2**2 + muu3**2) + 
     -               (g1R*g2R*g3Rc*g4Lc + g1L*g2L*g3Lc*g4Rc)*
     -                (-1 + muu1**2 + 2*muu3**2)*muu4 + 
     -               2*(g1L*g2L*g3Rc*g4Lc + g1R*g2R*g3Lc*g4Rc)*muu3*
     -                muu4**2 + 
     -               (g1R*g2R*g3Rc*g4Lc + g1L*g2L*g3Lc*g4Rc)*muu4**3)*
     -             muugl + (dPropFGlL*g1R*g2L*g3Rc*g4Lc + 
     -               dPropFGlR*g1L*g2R*g3Lc*g4Rc)*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2)*
     -             muugl**2)*tred - 
     -         (-((dPropFGlR*g1R*g2L*g3Rc*g4Lc + 
     -                 dPropFGlL*g1L*g2R*g3Lc*g4Rc)*
     -               (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2)) + 
     -            (dPropFGlL + dPropFGlR)*
     -             (g1R*g2R*(g3Lc*g4Rc*muu3 + g3Rc*g4Lc*muu4) + 
     -               g1L*g2L*(g3Rc*g4Lc*muu3 + g3Lc*g4Rc*muu4))*muugl +
     -              (dPropFGlL*g1R*g2L*g3Rc*g4Lc + 
     -               dPropFGlR*g1L*g2R*g3Lc*g4Rc)*muugl**2)*tred**2 - 
     -         (dPropFGlR*g1R*g2L*g3Rc*g4Lc + 
     -            dPropFGlL*g1L*g2R*g3Lc*g4Rc)*tred**3) + 
     -      dPropSGlL*(-(g1L*muugl*
     -            (g2R*g3Rc*g4Lc*muu3*muu4*
     -               (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -              g2L*g3Lc*g4Rc*muu4*muugl*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -              g2L*g3Rc*g4Lc*muu3*muugl*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred) + 
     -              g2R*g3Lc*g4Rc*
     -               ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -                 muu1**2*(-muu2**2 + tred) + 
     -                 muu3**2*(muu4**2 + tred)))) + 
     -         g1R*(g2R*(-(g3Rc*g4Lc*muu4*
     -                  (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)) -
     -                 g3Lc*g4Rc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             tred + g2L*muugl*
     -             (g3Lc*g4Rc*muu3*muu4*
     -                (muu1**2 + muu2**2 - muu3**2 - muu4**2) + 
     -               g3Rc*g4Lc*
     -                (muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -                  (muu1**2 + muu2**2 + muu3**2 + muu4**2)*
     -                   tred + tred**2)))) + 
     -      dPropSGlR*(-(g1R*muugl*
     -            (g2L*g3Lc*g4Rc*muu3*muu4*
     -               (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -              g2R*g3Rc*g4Lc*muu4*muugl*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -              g2R*g3Lc*g4Rc*muu3*muugl*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred) + 
     -              g2L*g3Rc*g4Lc*
     -               ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -                 muu1**2*(-muu2**2 + tred) + 
     -                 muu3**2*(muu4**2 + tred)))) + 
     -         g1L*(g2L*(-(g3Lc*g4Rc*muu4*
     -                  (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)) -
     -                 g3Rc*g4Lc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             tred + g2R*muugl*
     -             (g3Rc*g4Lc*muu3*muu4*
     -                (muu1**2 + muu2**2 - muu3**2 - muu4**2) + 
     -               g3Lc*g4Rc*
     -                (muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -                  (muu1**2 + muu2**2 + muu3**2 + muu4**2)*
     -                   tred + tred**2)))) + 
     -      (dPropSGlR*muugl*
     -          (g1R*g2L*(g3Rc*g4Rc*(muu1 - muu3)*(muu1 + muu3)*
     -                muu4 + 
     -               g3Lc*g4Lc*muu3*(muu2 - muu4)*(muu2 + muu4)) + 
     -            g1L*g2R*(g3Lc*g4Lc*(muu1 - muu3)*(muu1 + muu3)*
     -                muu4 + 
     -               g3Rc*g4Rc*muu3*(muu2 - muu4)*(muu2 + muu4)) - 
     -            g1R*g2R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -               g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*muugl) + 
     -         (-((dPropFGlL*g1R*g2L*g3Lc*g4Lc + 
     -                 dPropFGlR*g1L*g2R*g3Rc*g4Rc)*muu2**2*muu3) + 
     -            (dPropFGlR*g1L*g2R*g3Lc*g4Lc + 
     -               dPropFGlL*g1R*g2L*g3Rc*g4Rc)*(-muu1 + muu3)*
     -             (muu1 + muu3)*muu4 + 
     -            (dPropFGlL*g1R*g2L*g3Lc*g4Lc + 
     -               dPropFGlR*g1L*g2R*g3Rc*g4Rc)*muu3*muu4**2)*
     -          muugl**2*sqrtS - 
     -         dPropSGlR*(g1L*g2L*
     -             (2*g3Rc*g4Rc*muu3*muu4 + 
     -               g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)) + 
     -            g1L*g2R*(g3Rc*g4Rc*muu3 + g3Lc*g4Lc*muu4)*muugl + 
     -            g1R*g2L*(g3Lc*g4Lc*muu3 + g3Rc*g4Rc*muu4)*muugl)*
     -          tred + (-(dPropFGlR*g1R*g2L*g3Lc*g4Lc*muu2**2*muu3) - 
     -            dPropFGlL*g1L*g2R*g3Rc*g4Rc*muu2**2*muu3 - 
     -            dPropFGlL*g1L*g2R*g3Lc*g4Lc*muu1**2*muu4 - 
     -            dPropFGlR*g1R*g2L*g3Rc*g4Rc*muu1**2*muu4 + 
     -            dPropFGlL*g1L*g2R*g3Lc*g4Lc*muu3**2*muu4 + 
     -            dPropFGlR*g1R*g2L*g3Rc*g4Rc*muu3**2*muu4 + 
     -            dPropFGlR*g1R*g2L*g3Lc*g4Lc*muu3*muu4**2 + 
     -            dPropFGlL*g1L*g2R*g3Rc*g4Rc*muu3*muu4**2 + 
     -            (dPropFGlL + dPropFGlR)*
     -             (g1L*g2L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -                  g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)) + 
     -               g1R*g2R*
     -                (2*g3Lc*g4Lc*muu3*muu4 + 
     -                  g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)))*muugl + 
     -            (dPropFGlR*g1L*g2R*
     -                (g3Rc*g4Rc*muu3 + g3Lc*g4Lc*muu4) + 
     -               dPropFGlL*g1R*g2L*
     -                (g3Lc*g4Lc*muu3 + g3Rc*g4Rc*muu4))*muugl**2)*
     -          sqrtS*tred + 
     -         (dPropFGlL*g1L*g2R*(g3Rc*g4Rc*muu3 + g3Lc*g4Lc*muu4) + 
     -            dPropFGlR*g1R*g2L*(g3Lc*g4Lc*muu3 + g3Rc*g4Rc*muu4))
     -           *sqrtS*tred**2 + 
     -         dPropSGlL*(g1R*g2L*
     -             (g3Rc*g4Rc*(muu1 - muu3)*(muu1 + muu3)*muu4 + 
     -               g3Lc*g4Lc*muu3*(muu2 - muu4)*(muu2 + muu4))*muugl
     -              - g1R*(2*g2R*g3Lc*g4Lc*muu3*muu4 + 
     -               g2R*g3Rc*g4Rc*(-1 + muu3**2 + muu4**2) + 
     -               g2L*(g3Lc*g4Lc*muu3 + g3Rc*g4Rc*muu4)*muugl)*tred
     -              - g1L*muugl*
     -             (2*g2L*g3Rc*g4Rc*muu3*muu4*muugl + 
     -               g2L*g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2R*g3Lc*g4Lc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -               g2R*g3Rc*g4Rc*muu3*(-muu2**2 + muu4**2 + tred))))
     -        *muupu(n)))/(DenGlT**2*dconjg(DenGauU(n))))	 
	  
	  enddo
	  
! ! 	  
! ! ! ==========================================================u-channel correction=============================================================
! ! 	  
! ! 
! ! !---------------------u-channel Gluino squared------------------------------
! ! ! ! 
	  colfac = Nc/2d0*CF
    
	
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1Lc = dconjg(g1L)
	  g1Rc = dconjg(g1R)
	
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2Lc = dconjg(g2L)
	  g2Rc = dconjg(g2R)	  
	  
	  stst2QQ_M2CountertermsProp = stst2QQ_M2CountertermsProp + 2d0*colfac*
     -  	  dreal((sqrtS**5*(-(dPropSGlR*
     -         (g1R*muugl*(g1Rc*g2R*muugl*
     -               (g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -                 g2Lc*muu4*
     -                  (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)) +
     -                g1Lc*g2R*muu3*muugl*
     -               (2*g2Lc*muu4*muugl + 
     -                 g2Rc*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                    tred)) + 
     -              g1Lc*g2L*muu3*
     -               (2*g2Rc*muu4*
     -                  (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                    tred) + 
     -                 g2Lc*muugl*
     -                  (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) +
     -                g1Rc*g2L*
     -               (g2Rc*muu4*muugl*
     -                  (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -                 g2Lc*(muu3**4 - muu4**2 + 
     -                    muu3**2*
     -                     (-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -                    (muu4**2 - tred)**2 + 
     -                    muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -                    muu2**2*tred))) + 
     -           g1L*(g1Rc*g2R*muu3*muugl*
     -               (2*g2Lc*muu4*
     -                  (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                    tred) + 
     -                 g2Rc*muugl*
     -                  (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) +
     -                g1Lc*g2L*
     -               (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -                 g2Rc*muu4*
     -                  (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred))*
     -               (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                 tred) + 
     -              g1Rc*g2L*muu3*
     -               (2*g2Rc*muu4*muugl + 
     -                 g2Lc*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -                    tred))*
     -               (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                 tred) + 
     -              g1Lc*g2R*muugl*
     -               (g2Lc*muu4*muugl*
     -                  (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -                 g2Rc*(muu3**4 - muu4**2 + 
     -                    muu3**2*
     -                     (-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -                    (muu4**2 - tred)**2 + 
     -                    muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -                    muu2**2*tred))))) - 
     -      dPropSGlL*(g1L*g1Rc*muu3*muugl*
     -          (g2Rc*muugl*(2*g2L*muu4*muugl + 
     -               g2R*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) +
     -             g2Lc*(2*g2R*muu4*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred) + 
     -               g2L*muugl*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))) + 
     -         g1Lc*g1R*muu3*
     -          (g2L*muugl*(2*g2Rc*muu4*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred) + 
     -               g2Lc*muugl*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -            g2R*(2*g2Lc*muu4*muugl + 
     -               g2Rc*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))
     -              *(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)) + 
     -         g1L*g1Lc*muugl*
     -          (g2L*muugl*(g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2Rc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)) + 
     -            g2R*(g2Lc*muu4*muugl*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g2Rc*(muu3**4 - muu4**2 + 
     -                  muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -                  (muu4**2 - tred)**2 + 
     -                  muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -                  muu2**2*tred))) + 
     -         g1R*g1Rc*(g2Rc*muugl*
     -             (g2R*(-1 + muu3**2 + muu4**2)*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred) + 
     -               g2L*muu4*muugl*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)) + 
     -            g2Lc*(g2R*muu4*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g2L*muugl*
     -                (muu3**4 - muu4**2 + 
     -                  muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -                  (muu4**2 - tred)**2 + 
     -                  muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -                  muu2**2*tred)))) + 
     -      sqrtS*(dPropFGlR*g1R*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          (g1Rc*g2R*muugl*
     -             (g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2Lc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)) + 
     -            g1Lc*g2R*muu3*muugl*
     -             (2*g2Lc*muu4*muugl + 
     -               g2Rc*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))
     -              + g1Lc*g2L*muu3*
     -             (2*g2Rc*muu4*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred) + 
     -               g2Lc*muugl*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -            g1Rc*g2L*(g2Rc*muu4*muugl*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g2Lc*(muu3**4 - muu4**2 + 
     -                  muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -                  (muu4**2 - tred)**2 + 
     -                  muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -                  muu2**2*tred))) + 
     -         dPropFGlL*g1L*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          (g1Lc*g2L*muugl*
     -             (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2Rc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)) + 
     -            g1Rc*g2Rc*muu3*muugl*
     -             (2*g2L*muu4*muugl + 
     -               g2R*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) +
     -              g1Rc*g2Lc*muu3*
     -             (2*g2R*muu4*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred) + 
     -               g2L*muugl*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -            g1Lc*g2R*(g2Lc*muu4*muugl*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g2Rc*(muu3**4 - muu4**2 + 
     -                  muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -                  (muu4**2 - tred)**2 + 
     -                  muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -                  muu2**2*tred))) + 
     -         dPropFGlR*g1L*muugl*
     -          (g1Rc*g2R*muu3*muugl*
     -             (2*g2Lc*muu4*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred) + 
     -               g2Rc*muugl*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -            g1Lc*g2L*(g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2Rc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred))*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)
     -              + g1Rc*g2L*muu3*
     -             (2*g2Rc*muu4*muugl + 
     -               g2Lc*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))
     -              *(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred) + 
     -            g1Lc*g2R*muugl*
     -             (g2Lc*muu4*muugl*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g2Rc*(muu3**4 - muu4**2 + 
     -                  muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -                  (muu4**2 - tred)**2 + 
     -                  muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -                  muu2**2*tred))) + 
     -         dPropFGlL*g1R*muugl*
     -          (g1Rc*g2Rc*muugl*
     -             (g2R*(-1 + muu3**2 + muu4**2)*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred) + 
     -               g2L*muu4*muugl*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)) + 
     -            g1Lc*g2L*muu3*muugl*
     -             (2*g2Rc*muu4*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred) + 
     -               g2Lc*muugl*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -            g1Lc*g2R*muu3*
     -             (2*g2Lc*muu4*muugl + 
     -               g2Rc*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))
     -              *(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred) + 
     -            g1Rc*g2Lc*(g2R*muu4*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g2L*muugl*
     -                (muu3**4 - muu4**2 + 
     -                  muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -                  (muu4**2 - tred)**2 + 
     -                  muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -                  muu2**2*tred))))))/(DenGlU**2*dconjg(DenGlU)))	 
! 	  
! !      
! ! ! ---------------------t-channel and u-channel Gluino squared------------------------------
! ! 
! ! 
	  colfac = -1d0/2d0*CF
	
	  g1Lc = dconjg(GluLt1(isf1,itt1,ftt1,igen1,fgen1))
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	
	  g1Rc = dconjg(GluRt1(isf1,itt1,ftt1,igen1,fgen1))
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	
	  g2Lc = dconjg(GluLt1(isf2,itt2,ftt2,igen2,fgen2))
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	
	  g2Rc = dconjg(GluRt1(isf2,itt2,ftt2,igen2,fgen2))
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  
! 	  	  
	  stst2QQ_M2CountertermsProp = stst2QQ_M2CountertermsProp + 2d0*colfac*
     -  	  dreal((sqrtS**5*(dPropSGlL*
     -       (g1Lc*g1R*muu3*(g2L*muugl*
     -             (g2Rc*muu4*
     -                (muu1**2 + muu2**2 - muu3**2 - muu4**2) - 
     -               g2Lc*muugl*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) - 
     -            g2R*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*(2*g2Lc*muu4*muugl + 
     -               g2Rc*(-muu2**2 + muu4**2 + tred))) - 
     -         g1L*g1Rc*muu3*muugl*
     -          (g2Rc*muugl*(2*g2L*muu4*muugl + 
     -               g2R*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) +
     -              g2Lc*(g2R*muu4*
     -                (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -               g2L*muugl*(-muu2**2 + muu4**2 + tred))) + 
     -         g1R*g1Rc*(-(g2Rc*muugl*
     -               (g2R*(-1 + muu3**2 + muu4**2)*
     -                  (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                    tred) + 
     -                 g2L*muu4*muugl*
     -                  (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred))) +
     -              g2Lc*(g2R*muu4*(muu1**2 - muu3**2 - tred)*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred) + 
     -               g2L*muugl*
     -                (muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -                  (muu1**2 + muu2**2 + muu3**2 + muu4**2)*
     -                   tred + tred**2))) - 
     -         g1L*g1Lc*muugl*
     -          (g2L*muugl*(g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2Rc*muu4*(-muu1**2 + muu3**2 + tred)) + 
     -            g2R*(g2Lc*muu4*muugl*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g2Rc*((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -                  muu1**2*(-muu2**2 + tred) + 
     -                  muu3**2*(muu4**2 + tred))))) + 
     -      dPropSGlR*(-(g1R*muugl*
     -            (g1Lc*g2L*muu3*
     -               (g2Rc*muu4*
     -                  (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -                 g2Lc*muugl*
     -                  (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) +
     -                g1Rc*g2R*muugl*
     -               (g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -                 g2Lc*muu4*(-muu1**2 + muu3**2 + tred)) + 
     -              g1Lc*g2R*muu3*muugl*
     -               (2*g2Lc*muu4*muugl + 
     -                 g2Rc*(-muu2**2 + muu4**2 + tred)) + 
     -              g1Rc*g2L*
     -               (g2Rc*muu4*muugl*
     -                  (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -                 g2Lc*((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -                    muu1**2*(-muu2**2 + tred) + 
     -                    muu3**2*(muu4**2 + tred))))) + 
     -         g1L*(g1Rc*g2R*muu3*muugl*
     -             (g2Lc*muu4*
     -                (muu1**2 + muu2**2 - muu3**2 - muu4**2) - 
     -               g2Rc*muugl*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -            g1Lc*g2L*(-(g2Lc*(-1 + muu3**2 + muu4**2)*muugl) + 
     -               g2Rc*muu4*(muu1**2 - muu3**2 - tred))*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)
     -              - g1Rc*g2L*muu3*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*(2*g2Rc*muu4*muugl + 
     -               g2Lc*(-muu2**2 + muu4**2 + tred)) - 
     -            g1Lc*g2R*muugl*
     -             (g2Lc*muu4*muugl*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g2Rc*((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -                  muu1**2*(-muu2**2 + tred) + 
     -                  muu3**2*(muu4**2 + tred))))) + 
     -      sqrtS*(dPropFGlR*g1L*muugl*
     -          (g1Rc*g2R*muu3*muugl*
     -             (g2Lc*muu4*
     -                (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -               g2Rc*muugl*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) - 
     -            g1Lc*g2L*(-(g2Lc*(-1 + muu3**2 + muu4**2)*muugl) + 
     -               g2Rc*muu4*(muu1**2 - muu3**2 - tred))*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)
     -              + g1Rc*g2L*muu3*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*(2*g2Rc*muu4*muugl + 
     -               g2Lc*(-muu2**2 + muu4**2 + tred)) + 
     -            g1Lc*g2R*muugl*
     -             (g2Lc*muu4*muugl*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g2Rc*((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -                  muu1**2*(-muu2**2 + tred) + 
     -                  muu3**2*(muu4**2 + tred)))) + 
     -         dPropFGlL*g1R*muugl*
     -          (g1Rc*g2Rc*muugl*
     -             (g2R*(-1 + muu3**2 + muu4**2)*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred) + 
     -               g2L*muu4*muugl*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)) + 
     -            g1Lc*g2L*muu3*muugl*
     -             (g2Rc*muu4*
     -                (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -               g2Lc*muugl*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -            g1Lc*g2R*muu3*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred)*(2*g2Lc*muu4*muugl + 
     -               g2Rc*(-muu2**2 + muu4**2 + tred)) + 
     -            g1Rc*g2Lc*(-(g2R*muu4*(muu1**2 - muu3**2 - tred)*
     -                  (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                    tred)) + 
     -               g2L*muugl*
     -                ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -                  muu1**2*(-muu2**2 + tred) + 
     -                  muu3**2*(muu4**2 + tred)))) - 
     -         dPropFGlR*g1R*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          (-(g1Lc*muu3*
     -               (g2L*(g2Rc*muu4*
     -                     (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -                    g2Lc*muugl*
     -                     (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)
     -                    ) + 
     -                 g2R*muugl*
     -                  (2*g2Lc*muu4*muugl + 
     -                    g2Rc*(-muu2**2 + muu4**2 + tred)))) + 
     -            g1Rc*(-(g2R*muugl*
     -                  (g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -                    g2Lc*muu4*(-muu1**2 + muu3**2 + tred))) + 
     -               g2L*(-(g2Rc*muu4*muugl*
     -                     (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)
     -                     ) + 
     -                  g2Lc*
     -                   (muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -                     (muu1**2 + muu2**2 + muu3**2 + muu4**2)*
     -                      tred + tred**2)))) - 
     -         dPropFGlL*g1L*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)*
     -          (-(g1Rc*muu3*
     -               (g2Rc*muugl*
     -                  (2*g2L*muu4*muugl + 
     -                    g2R*
     -                     (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)
     -                    ) + 
     -                 g2Lc*(g2R*muu4*
     -                     (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -                    g2L*muugl*(-muu2**2 + muu4**2 + tred)))) + 
     -            g1Lc*(-(g2L*muugl*
     -                  (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -                    g2Rc*muu4*(-muu1**2 + muu3**2 + tred))) + 
     -               g2R*(-(g2Lc*muu4*muugl*
     -                     (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)
     -                     ) + 
     -                  g2Rc*
     -                   (muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -                     (muu1**2 + muu2**2 + muu3**2 + muu4**2)*
     -                      tred + tred**2)))))))/
     -  (DenGlU**2*dconjg(DenGlT)))	 
	          
! !---------------------u-channel Gluino t-channel Neutralino---------------------------------------
      do n=1,6
	  colfac = CF*Nc
	  
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)

	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)

	  g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))

	  g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2)) 
	  
	  	  
	  stst2QQ_M2CountertermsProp = stst2QQ_M2CountertermsProp + 2d0*colfac*
     -  	  dreal((sqrtS**5*(-(dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3) + 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu1**2*muu2**2*muu3 + 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu2**4*muu3 + 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3**3 - 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4 + 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu1**4*muu4 + 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu2**2*muu4 + 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu3**2*muu4 - 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu2**2*muu3**2*muu4 - 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu3**4*muu4 + 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu3*muu4**2 - 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3*muu4**2 - 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu3**3*muu4**2 + 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4**3 - 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu3**2*muu4**3 - 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu3*muu4**4 + 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu1**2*muu2**2*muugl + 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu2**2*muugl + 
     -      dPropSGlL*g1L*g2R*g3Rc*g4Lc*muu1**2*muu3*muu4*muugl + 
     -      dPropSGlL*g1R*g2L*g3Lc*g4Rc*muu1**2*muu3*muu4*muugl + 
     -      dPropSGlL*g1L*g2R*g3Rc*g4Lc*muu2**2*muu3*muu4*muugl + 
     -      dPropSGlL*g1R*g2L*g3Lc*g4Rc*muu2**2*muu3*muu4*muugl - 
     -      dPropSGlL*g1L*g2R*g3Rc*g4Lc*muu3**3*muu4*muugl - 
     -      dPropSGlL*g1R*g2L*g3Lc*g4Rc*muu3**3*muu4*muugl - 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu3**2*muu4**2*muugl - 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*muu4**2*muugl - 
     -      dPropSGlL*g1L*g2R*g3Rc*g4Lc*muu3*muu4**3*muugl - 
     -      dPropSGlL*g1R*g2L*g3Lc*g4Rc*muu3*muu4**3*muugl + 
     -      dPropSGlL*g1L*g2L*g3Rc*g4Lc*muu2**2*muu3*muugl**2 + 
     -      dPropSGlL*g1L*g2L*g3Lc*g4Rc*muu1**2*muu4*muugl**2 - 
     -      dPropSGlL*g1L*g2L*g3Lc*g4Rc*muu3**2*muu4*muugl**2 - 
     -      dPropSGlL*g1L*g2L*g3Rc*g4Lc*muu3*muu4**2*muugl**2 + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu2**2*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu2**2*sqrtS - 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**4*muu2**2*sqrtS - 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**4*muu2**2*sqrtS - 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu2**4*sqrtS - 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu2**4*sqrtS - 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu2**2*muu3**2*
     -       sqrtS - dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu2**2*
     -       muu3**2*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu1**2*muu3*muu4*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu1**2*muu3*muu4*sqrtS - 
     -      dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu1**4*muu3*muu4*sqrtS - 
     -      dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu1**4*muu3*muu4*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu2**2*muu3*muu4*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu2**2*muu3*muu4*sqrtS - 
     -      2*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu1**2*muu2**2*muu3*muu4*
     -       sqrtS - 2*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu1**2*muu2**2*
     -       muu3*muu4*sqrtS - 
     -      dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu2**4*muu3*muu4*sqrtS - 
     -      dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu2**4*muu3*muu4*sqrtS - 
     -      dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3**3*muu4*sqrtS - 
     -      dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3**3*muu4*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3**5*muu4*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3**5*muu4*sqrtS - 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu2**2*muu4**2*
     -       sqrtS - dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu2**2*
     -       muu4**2*sqrtS - 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**2*muu4**2*sqrtS - 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*muu4**2*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu3**2*muu4**2*
     -       sqrtS + dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu3**2*
     -       muu4**2*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**2*muu3**2*muu4**2*
     -       sqrtS + dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*muu3**2*
     -       muu4**2*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**4*muu4**2*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**4*muu4**2*sqrtS - 
     -      dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3*muu4**3*sqrtS - 
     -      dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3*muu4**3*sqrtS + 
     -      2*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3**3*muu4**3*sqrtS + 
     -      2*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3**3*muu4**3*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**2*muu4**4*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*muu4**4*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3*muu4**5*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3*muu4**5*sqrtS + 
     -      dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu2**2*muu3*muugl*sqrtS + 
     -      dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu2**2*muu3*muugl*sqrtS + 
     -      dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3*muugl*sqrtS + 
     -      dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3*muugl*sqrtS - 
     -      dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu1**2*muu2**2*muu3*muugl*
     -       sqrtS - dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu1**2*muu2**2*muu3*
     -       muugl*sqrtS - dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu1**2*
     -       muu2**2*muu3*muugl*sqrtS - 
     -      dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu1**2*muu2**2*muu3*muugl*
     -       sqrtS - dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu2**4*muu3*muugl*
     -       sqrtS - dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu2**4*muu3*muugl*
     -       sqrtS - dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu2**4*muu3*muugl*
     -       sqrtS - dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu2**4*muu3*muugl*
     -       sqrtS - dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu2**2*muu3**3*
     -       muugl*sqrtS - dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu2**2*
     -       muu3**3*muugl*sqrtS - 
     -      dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3**3*muugl*sqrtS - 
     -      dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3**3*muugl*sqrtS + 
     -      dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4*muugl*sqrtS + 
     -      dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4*muugl*sqrtS + 
     -      dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu1**2*muu4*muugl*sqrtS + 
     -      dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu1**2*muu4*muugl*sqrtS - 
     -      dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu1**4*muu4*muugl*sqrtS - 
     -      dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu1**4*muu4*muugl*sqrtS - 
     -      dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu1**4*muu4*muugl*sqrtS - 
     -      dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu1**4*muu4*muugl*sqrtS - 
     -      dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu2**2*muu4*muugl*
     -       sqrtS - dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu1**2*muu2**2*muu4*
     -       muugl*sqrtS - dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu1**2*
     -       muu2**2*muu4*muugl*sqrtS - 
     -      dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu1**2*muu2**2*muu4*muugl*
     -       sqrtS - dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu3**2*muu4*muugl*
     -       sqrtS - dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu3**2*muu4*muugl*
     -       sqrtS - dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu3**2*muu4*muugl*
     -       sqrtS - dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu3**2*muu4*muugl*
     -       sqrtS + dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu2**2*muu3**2*muu4*
     -       muugl*sqrtS + dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu2**2*
     -       muu3**2*muu4*muugl*sqrtS + 
     -      dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu2**2*muu3**2*muu4*muugl*
     -       sqrtS + dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu2**2*muu3**2*muu4*
     -       muugl*sqrtS + dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu3**4*muu4*
     -       muugl*sqrtS + dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu3**4*muu4*
     -       muugl*sqrtS + dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu3**4*muu4*
     -       muugl*sqrtS + dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu3**4*muu4*
     -       muugl*sqrtS - dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu3*muu4**2*
     -       muugl*sqrtS - dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu3*muu4**2*
     -       muugl*sqrtS - dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu3*muu4**2*
     -       muugl*sqrtS - dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu3*muu4**2*
     -       muugl*sqrtS + dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu1**2*muu3*
     -       muu4**2*muugl*sqrtS + 
     -      dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu1**2*muu3*muu4**2*muugl*
     -       sqrtS + dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3*muu4**2*
     -       muugl*sqrtS + dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3*
     -       muu4**2*muugl*sqrtS + 
     -      dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu3**3*muu4**2*muugl*sqrtS + 
     -      dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu3**3*muu4**2*muugl*sqrtS + 
     -      dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu3**3*muu4**2*muugl*sqrtS + 
     -      dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu3**3*muu4**2*muugl*sqrtS - 
     -      dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4**3*muugl*sqrtS - 
     -      dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4**3*muugl*sqrtS - 
     -      dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu1**2*muu4**3*muugl*sqrtS - 
     -      dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu1**2*muu4**3*muugl*sqrtS + 
     -      dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu3**2*muu4**3*muugl*sqrtS + 
     -      dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu3**2*muu4**3*muugl*sqrtS + 
     -      dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu3**2*muu4**3*muugl*sqrtS + 
     -      dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu3**2*muu4**3*muugl*sqrtS + 
     -      dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu3*muu4**4*muugl*sqrtS + 
     -      dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu3*muu4**4*muugl*sqrtS + 
     -      dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu3*muu4**4*muugl*sqrtS + 
     -      dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu3*muu4**4*muugl*sqrtS - 
     -      dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu1**2*muu2**2*muugl**2*
     -       sqrtS - dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu1**2*muu2**2*
     -       muugl**2*sqrtS - 
     -      dPropFGlR*g1L*g2R*g3Rc*g4Lc*muu1**2*muu3*muu4*muugl**2*
     -       sqrtS - dPropFGlL*g1R*g2L*g3Lc*g4Rc*muu1**2*muu3*muu4*
     -       muugl**2*sqrtS - 
     -      dPropFGlR*g1L*g2R*g3Rc*g4Lc*muu2**2*muu3*muu4*muugl**2*
     -       sqrtS - dPropFGlL*g1R*g2L*g3Lc*g4Rc*muu2**2*muu3*muu4*
     -       muugl**2*sqrtS + 
     -      dPropFGlR*g1L*g2R*g3Rc*g4Lc*muu3**3*muu4*muugl**2*sqrtS + 
     -      dPropFGlL*g1R*g2L*g3Lc*g4Rc*muu3**3*muu4*muugl**2*sqrtS + 
     -      dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu3**2*muu4**2*muugl**2*
     -       sqrtS + dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu3**2*muu4**2*
     -       muugl**2*sqrtS + 
     -      dPropFGlR*g1L*g2R*g3Rc*g4Lc*muu3*muu4**3*muugl**2*sqrtS + 
     -      dPropFGlL*g1R*g2L*g3Lc*g4Rc*muu3*muu4**3*muugl**2*sqrtS + 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu3*tred - 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3*tred - 
     -      2*dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3*tred - 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu3**3*tred + 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu4*tred - 
     -      2*dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4*tred - 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu2**2*muu4*tred - 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu4**3*tred + 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muugl*tred + 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muugl*tred - 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu1**2*muugl*tred - 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muugl*tred - 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu2**2*muugl*tred - 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*muugl*tred - 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu3**2*muugl*tred - 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*muugl*tred - 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu4**2*muugl*tred - 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu4**2*muugl*tred - 
     -      dPropSGlL*g1L*g2L*g3Rc*g4Lc*muu3*muugl**2*tred - 
     -      dPropSGlL*g1L*g2L*g3Lc*g4Rc*muu4*muugl**2*tred + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*sqrtS*tred + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*sqrtS*tred - 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*sqrtS*tred - 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*sqrtS*tred + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**4*sqrtS*tred + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**4*sqrtS*tred - 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**2*sqrtS*tred - 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*sqrtS*tred + 
     -      3*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu2**2*sqrtS*
     -       tred + 3*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu2**2*
     -       sqrtS*tred + dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**4*sqrtS*
     -       tred + dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**4*sqrtS*tred - 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**2*sqrtS*tred - 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*sqrtS*tred + 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu3**2*sqrtS*
     -       tred + 2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu3**2*
     -       sqrtS*tred + 2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**2*
     -       muu3**2*sqrtS*tred + 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*muu3**2*sqrtS*
     -       tred + dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**4*sqrtS*tred + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**4*sqrtS*tred + 
     -      dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu1**2*muu3*muu4*sqrtS*
     -       tred + dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu1**2*muu3*muu4*
     -       sqrtS*tred + dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu2**2*muu3*
     -       muu4*sqrtS*tred + 
     -      dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu2**2*muu3*muu4*sqrtS*
     -       tred - dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3**3*muu4*sqrtS*
     -       tred - dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3**3*muu4*sqrtS*
     -       tred - 2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu4**2*sqrtS*
     -       tred - 2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu4**2*sqrtS*
     -       tred + 2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu4**2*
     -       sqrtS*tred + 2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*
     -       muu4**2*sqrtS*tred + 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**2*muu4**2*sqrtS*
     -       tred + 2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*muu4**2*
     -       sqrtS*tred + dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**2*muu4**2*
     -       sqrtS*tred + dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*muu4**2*
     -       sqrtS*tred - dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3*muu4**3*
     -       sqrtS*tred - dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3*muu4**3*
     -       sqrtS*tred + dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu4**4*sqrtS*
     -       tred + dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu4**4*sqrtS*tred - 
     -      dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu3*muugl*sqrtS*tred - 
     -      dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu3*muugl*sqrtS*tred - 
     -      dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu3*muugl*sqrtS*tred - 
     -      dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu3*muugl*sqrtS*tred + 
     -      dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu1**2*muu3*muugl*sqrtS*
     -       tred + dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu1**2*muu3*muugl*
     -       sqrtS*tred + dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3*
     -       muugl*sqrtS*tred + 
     -      dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3*muugl*sqrtS*
     -       tred + 2*dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu2**2*muu3*muugl*
     -       sqrtS*tred + 2*dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu2**2*muu3*
     -       muugl*sqrtS*tred + 
     -      2*dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3*muugl*sqrtS*
     -       tred + 2*dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3*muugl*
     -       sqrtS*tred + dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu3**3*muugl*
     -       sqrtS*tred + dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu3**3*muugl*
     -       sqrtS*tred + dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu3**3*muugl*
     -       sqrtS*tred + dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu3**3*muugl*
     -       sqrtS*tred - dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu4*muugl*
     -       sqrtS*tred - dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu4*muugl*
     -       sqrtS*tred - dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu4*muugl*
     -       sqrtS*tred - dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu4*muugl*
     -       sqrtS*tred + 2*dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4*
     -       muugl*sqrtS*tred + 
     -      2*dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4*muugl*sqrtS*
     -       tred + 2*dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu1**2*muu4*muugl*
     -       sqrtS*tred + 2*dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu1**2*muu4*
     -       muugl*sqrtS*tred + 
     -      dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu2**2*muu4*muugl*sqrtS*
     -       tred + dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu2**2*muu4*muugl*
     -       sqrtS*tred + dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu2**2*muu4*
     -       muugl*sqrtS*tred + 
     -      dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu2**2*muu4*muugl*sqrtS*
     -       tred + dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu4**3*muugl*sqrtS*
     -       tred + dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu4**3*muugl*sqrtS*
     -       tred + dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu4**3*muugl*sqrtS*
     -       tred + dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu4**3*muugl*sqrtS*
     -       tred - dPropFGlL*g1R*g2L*g3Rc*g4Lc*muugl**2*sqrtS*tred - 
     -      dPropFGlR*g1L*g2R*g3Lc*g4Rc*muugl**2*sqrtS*tred + 
     -      dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu1**2*muugl**2*sqrtS*tred + 
     -      dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu1**2*muugl**2*sqrtS*tred + 
     -      dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu2**2*muugl**2*sqrtS*tred + 
     -      dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu2**2*muugl**2*sqrtS*tred + 
     -      dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu3**2*muugl**2*sqrtS*tred + 
     -      dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu3**2*muugl**2*sqrtS*tred + 
     -      dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu4**2*muugl**2*sqrtS*tred + 
     -      dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu4**2*muugl**2*sqrtS*tred + 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu3*tred**2 + 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu4*tred**2 + 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muugl*tred**2 + 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muugl*tred**2 + 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*sqrtS*tred**2 + 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*sqrtS*tred**2 - 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*sqrtS*tred**2 - 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*sqrtS*tred**2 - 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**2*sqrtS*tred**2 - 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*sqrtS*tred**2 - 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**2*sqrtS*tred**2 - 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*sqrtS*tred**2 - 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu4**2*sqrtS*tred**2 - 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu4**2*sqrtS*tred**2 - 
     -      dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu3*muugl*sqrtS*tred**2 - 
     -      dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu3*muugl*sqrtS*tred**2 - 
     -      dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu3*muugl*sqrtS*tred**2 - 
     -      dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu3*muugl*sqrtS*tred**2 - 
     -      dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu4*muugl*sqrtS*tred**2 - 
     -      dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu4*muugl*sqrtS*tred**2 - 
     -      dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu4*muugl*sqrtS*tred**2 - 
     -      dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu4*muugl*sqrtS*tred**2 - 
     -      dPropFGlL*g1R*g2L*g3Rc*g4Lc*muugl**2*sqrtS*tred**2 - 
     -      dPropFGlR*g1L*g2R*g3Lc*g4Rc*muugl**2*sqrtS*tred**2 + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*sqrtS*tred**3 + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*sqrtS*tred**3 + 
     -      dPropSGlR*g1R*muugl*
     -       (g2L*g3Lc*g4Rc*muu3*muu4*
     -          (muu1**2 + muu2**2 - muu3**2 - muu4**2) + 
     -         g2R*g3Rc*g4Lc*muu4*muugl*(muu1**2 - muu3**2 - tred) + 
     -         g2R*g3Lc*g4Rc*muu3*muugl*(muu2**2 - muu4**2 - tred) + 
     -         g2L*g3Rc*g4Lc*
     -          (muu1**2*(muu2**2 - tred) + 
     -            tred*(1 - muu2**2 - muu4**2 + tred) - 
     -            muu3**2*(muu4**2 + tred))) + 
     -      dPropSGlR*g1L*(g2L*
     -          (g3Lc*g4Rc*muu4*(muu1**2 - muu3**2 - tred) + 
     -            g3Rc*g4Lc*muu3*(muu2**2 - muu4**2 - tred))*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred) + 
     -         g2R*muugl*(g3Rc*g4Lc*muu3*muu4*
     -             (muu1**2 + muu2**2 - muu3**2 - muu4**2) + 
     -            g3Lc*g4Rc*(muu1**2*(muu2**2 - tred) + 
     -               tred*(1 - muu2**2 - muu4**2 + tred) - 
     -               muu3**2*(muu4**2 + tred)))) + 
     -      (-(dPropSGlR*(g1L*
     -               (g2R*muugl*
     -                  (g3Lc*g4Lc*muu4*
     -                     (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)
     -                      + g3Rc*g4Rc*muu3*
     -                     (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)
     -                    ) + 
     -                 g2L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -                    g3Lc*g4Lc*(-1 + muu3**2 + muu4**2))*
     -                  (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                    tred)) + 
     -              g1R*muugl*
     -               (2*g2R*g3Lc*g4Lc*muu3*muu4*muugl + 
     -                 g2R*g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -                 g2L*g3Rc*g4Rc*muu4*
     -                  (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -                 g2L*g3Lc*g4Lc*muu3*
     -                  (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))))
     -           - dPropSGlL*
     -          (g1R*(g2L*muugl*
     -                (g3Rc*g4Rc*muu4*
     -                   (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) +
     -                    g3Lc*g4Lc*muu3*
     -                   (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) +
     -                 g2R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -                  g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)) + 
     -            g1L*muugl*(2*g2L*g3Rc*g4Rc*muu3*muu4*muugl + 
     -               g2L*g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2R*g3Lc*g4Lc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g2R*g3Rc*g4Rc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))) + 
     -         sqrtS*(dPropFGlR*g1L*muugl*
     -             (g2R*muugl*
     -                (g3Lc*g4Lc*muu4*
     -                   (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) +
     -                    g3Rc*g4Rc*muu3*
     -                   (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) +
     -                 g2L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -                  g3Lc*g4Lc*(-1 + muu3**2 + muu4**2))*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)) + 
     -            dPropFGlL*g1R*muugl*
     -             (g2L*muugl*
     -                (g3Rc*g4Rc*muu4*
     -                   (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) +
     -                    g3Lc*g4Lc*muu3*
     -                   (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) +
     -                 g2R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -                  g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)) + 
     -            dPropFGlR*g1R*
     -             (2*g2R*g3Lc*g4Lc*muu3*muu4*muugl + 
     -               g2R*g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2L*g3Rc*g4Rc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g2L*g3Lc*g4Lc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)
     -              + dPropFGlL*g1L*
     -             (2*g2L*g3Rc*g4Rc*muu3*muu4*muugl + 
     -               g2L*g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2R*g3Lc*g4Lc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g2R*g3Rc*g4Rc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)
     -            ))*muupt(n)))/(DenGlU**2*dconjg(DenGauT(n))))	 
 	  enddo
! !---------------------u-channel Gluino u-channel Neutralino--------------------------------------------------
      do n=1,6
	  colfac = 0d0
	  
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  
	  
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))	  
	  	  
	  stst2QQ_M2CountertermsProp = stst2QQ_M2CountertermsProp + 2d0*colfac*
     -  	  dreal((sqrtS**5*(-(dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu3) + 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3 + 
     -      2*dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3 - 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu1**2*muu2**2*muu3 - 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu2**4*muu3 + 
     -      2*dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu3**3 - 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3**3 - 
     -      2*dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3**3 - 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu3**5 - 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu4 + 
     -      2*dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4 - 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu1**4*muu4 + 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu2**2*muu4 - 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu2**2*muu4 + 
     -      3*dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu3**2*muu4 - 
     -      3*dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu3**2*muu4 - 
     -      2*dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu2**2*muu3**2*muu4 - 
     -      2*dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu3**4*muu4 + 
     -      3*dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu3*muu4**2 - 
     -      2*dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3*muu4**2 - 
     -      3*dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3*muu4**2 - 
     -      3*dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu3**3*muu4**2 + 
     -      2*dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu4**3 - 
     -      2*dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4**3 - 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu2**2*muu4**3 - 
     -      3*dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu3**2*muu4**3 - 
     -      2*dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu3*muu4**4 - 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu4**5 - 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu1**2*muu2**2*muugl - 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu2**2*muugl + 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu3**2*muugl + 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*muugl - 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu2**2*muu3**2*muugl - 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*muu3**2*muugl - 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu3**4*muugl - 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu3**4*muugl + 
     -      2*dPropSGlL*g1L*g2R*g3Rc*g4Lc*muu3*muu4*muugl + 
     -      2*dPropSGlL*g1R*g2L*g3Lc*g4Rc*muu3*muu4*muugl - 
     -      2*dPropSGlL*g1L*g2R*g3Rc*g4Lc*muu1**2*muu3*muu4*muugl - 
     -      2*dPropSGlL*g1R*g2L*g3Lc*g4Rc*muu1**2*muu3*muu4*muugl - 
     -      2*dPropSGlL*g1L*g2R*g3Rc*g4Lc*muu2**2*muu3*muu4*muugl - 
     -      2*dPropSGlL*g1R*g2L*g3Lc*g4Rc*muu2**2*muu3*muu4*muugl - 
     -      2*dPropSGlL*g1L*g2R*g3Rc*g4Lc*muu3**3*muu4*muugl - 
     -      2*dPropSGlL*g1R*g2L*g3Lc*g4Rc*muu3**3*muu4*muugl + 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu4**2*muugl + 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu4**2*muugl - 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu1**2*muu4**2*muugl - 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu4**2*muugl - 
     -      3*dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu3**2*muu4**2*muugl - 
     -      3*dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*muu4**2*muugl - 
     -      2*dPropSGlL*g1L*g2R*g3Rc*g4Lc*muu3*muu4**3*muugl - 
     -      2*dPropSGlL*g1R*g2L*g3Lc*g4Rc*muu3*muu4**3*muugl - 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu4**4*muugl - 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu4**4*muugl + 
     -      dPropSGlL*g1L*g2L*g3Rc*g4Lc*muu3*muugl**2 - 
     -      dPropSGlL*g1L*g2L*g3Rc*g4Lc*muu2**2*muu3*muugl**2 - 
     -      dPropSGlL*g1L*g2L*g3Rc*g4Lc*muu3**3*muugl**2 + 
     -      dPropSGlL*g1L*g2L*g3Lc*g4Rc*muu4*muugl**2 - 
     -      dPropSGlL*g1L*g2L*g3Lc*g4Rc*muu1**2*muu4*muugl**2 - 
     -      2*dPropSGlL*g1L*g2L*g3Lc*g4Rc*muu3**2*muu4*muugl**2 - 
     -      2*dPropSGlL*g1L*g2L*g3Rc*g4Lc*muu3*muu4**2*muugl**2 - 
     -      dPropSGlL*g1L*g2L*g3Lc*g4Rc*muu4**3*muugl**2 - 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu2**2*sqrtS - 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu2**2*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**4*muu2**2*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**4*muu2**2*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu2**4*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu2**4*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**2*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*sqrtS - 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu3**2*sqrtS - 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu3**2*sqrtS - 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**2*muu3**2*sqrtS - 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*muu3**2*sqrtS + 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu2**2*muu3**2*
     -       sqrtS + 2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu2**2*
     -       muu3**2*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**4*muu3**2*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**4*muu3**2*sqrtS - 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**4*sqrtS - 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**4*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu3**4*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu3**4*sqrtS + 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**2*muu3**4*sqrtS + 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*muu3**4*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**6*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**6*sqrtS + 
     -      2*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3*muu4*sqrtS + 
     -      2*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3*muu4*sqrtS - 
     -      4*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu1**2*muu3*muu4*sqrtS - 
     -      4*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu1**2*muu3*muu4*sqrtS + 
     -      2*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu1**4*muu3*muu4*sqrtS + 
     -      2*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu1**4*muu3*muu4*sqrtS - 
     -      4*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu2**2*muu3*muu4*sqrtS - 
     -      4*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu2**2*muu3*muu4*sqrtS + 
     -      4*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu1**2*muu2**2*muu3*muu4*
     -       sqrtS + 4*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu1**2*muu2**2*
     -       muu3*muu4*sqrtS + 
     -      2*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu2**4*muu3*muu4*sqrtS + 
     -      2*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu2**4*muu3*muu4*sqrtS - 
     -      4*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3**3*muu4*sqrtS - 
     -      4*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3**3*muu4*sqrtS + 
     -      4*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu1**2*muu3**3*muu4*
     -       sqrtS + 4*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu1**2*muu3**3*
     -       muu4*sqrtS + 4*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu2**2*
     -       muu3**3*muu4*sqrtS + 
     -      4*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu2**2*muu3**3*muu4*
     -       sqrtS + 2*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3**5*muu4*
     -       sqrtS + 2*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3**5*muu4*
     -       sqrtS + dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu4**2*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu4**2*sqrtS - 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu4**2*sqrtS - 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu4**2*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**4*muu4**2*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**4*muu4**2*sqrtS - 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**2*muu4**2*sqrtS - 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*muu4**2*sqrtS + 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu2**2*muu4**2*
     -       sqrtS + 2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu2**2*
     -       muu4**2*sqrtS - 
     -      5*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**2*muu4**2*sqrtS - 
     -      5*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*muu4**2*sqrtS + 
     -      4*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu3**2*muu4**2*
     -       sqrtS + 4*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu3**2*
     -       muu4**2*sqrtS + 
     -      4*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**2*muu3**2*muu4**2*
     -       sqrtS + 4*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*muu3**2*
     -       muu4**2*sqrtS + 
     -      4*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**4*muu4**2*sqrtS + 
     -      4*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**4*muu4**2*sqrtS - 
     -      4*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3*muu4**3*sqrtS - 
     -      4*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3*muu4**3*sqrtS + 
     -      4*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu1**2*muu3*muu4**3*
     -       sqrtS + 4*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu1**2*muu3*
     -       muu4**3*sqrtS + 
     -      4*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu2**2*muu3*muu4**3*
     -       sqrtS + 4*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu2**2*muu3*
     -       muu4**3*sqrtS + 
     -      4*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3**3*muu4**3*sqrtS + 
     -      4*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3**3*muu4**3*sqrtS - 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu4**4*sqrtS - 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu4**4*sqrtS + 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu4**4*sqrtS + 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu4**4*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**2*muu4**4*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*muu4**4*sqrtS + 
     -      4*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**2*muu4**4*sqrtS + 
     -      4*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*muu4**4*sqrtS + 
     -      2*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3*muu4**5*sqrtS + 
     -      2*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3*muu4**5*sqrtS + 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu4**6*sqrtS + 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu4**6*sqrtS + 
     -      dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu3*muugl*sqrtS + 
     -      dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu3*muugl*sqrtS + 
     -      dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu3*muugl*sqrtS + 
     -      dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu3*muugl*sqrtS - 
     -      dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu1**2*muu3*muugl*sqrtS - 
     -      dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu1**2*muu3*muugl*sqrtS - 
     -      dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3*muugl*sqrtS - 
     -      dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3*muugl*sqrtS - 
     -      2*dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu2**2*muu3*muugl*sqrtS - 
     -      2*dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu2**2*muu3*muugl*sqrtS - 
     -      2*dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3*muugl*sqrtS - 
     -      2*dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3*muugl*sqrtS + 
     -      dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu1**2*muu2**2*muu3*muugl*
     -       sqrtS + dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu1**2*muu2**2*muu3*
     -       muugl*sqrtS + dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu1**2*
     -       muu2**2*muu3*muugl*sqrtS + 
     -      dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu1**2*muu2**2*muu3*muugl*
     -       sqrtS + dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu2**4*muu3*muugl*
     -       sqrtS + dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu2**4*muu3*muugl*
     -       sqrtS + dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu2**4*muu3*muugl*
     -       sqrtS + dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu2**4*muu3*muugl*
     -       sqrtS - 2*dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu3**3*muugl*
     -       sqrtS - 2*dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu3**3*muugl*
     -       sqrtS - 2*dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu3**3*muugl*
     -       sqrtS - 2*dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu3**3*muugl*
     -       sqrtS + dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu1**2*muu3**3*
     -       muugl*sqrtS + dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu1**2*
     -       muu3**3*muugl*sqrtS + 
     -      dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3**3*muugl*sqrtS + 
     -      dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3**3*muugl*sqrtS + 
     -      2*dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu2**2*muu3**3*muugl*
     -       sqrtS + 2*dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu2**2*muu3**3*
     -       muugl*sqrtS + 2*dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu2**2*
     -       muu3**3*muugl*sqrtS + 
     -      2*dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3**3*muugl*
     -       sqrtS + dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu3**5*muugl*
     -       sqrtS + dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu3**5*muugl*
     -       sqrtS + dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu3**5*muugl*
     -       sqrtS + dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu3**5*muugl*
     -       sqrtS + dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu4*muugl*sqrtS + 
     -      dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu4*muugl*sqrtS + 
     -      dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu4*muugl*sqrtS + 
     -      dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu4*muugl*sqrtS - 
     -      2*dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4*muugl*sqrtS - 
     -      2*dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4*muugl*sqrtS - 
     -      2*dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu1**2*muu4*muugl*sqrtS - 
     -      2*dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu1**2*muu4*muugl*sqrtS + 
     -      dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu1**4*muu4*muugl*sqrtS + 
     -      dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu1**4*muu4*muugl*sqrtS + 
     -      dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu1**4*muu4*muugl*sqrtS + 
     -      dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu1**4*muu4*muugl*sqrtS - 
     -      dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu2**2*muu4*muugl*sqrtS - 
     -      dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu2**2*muu4*muugl*sqrtS - 
     -      dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu2**2*muu4*muugl*sqrtS - 
     -      dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu2**2*muu4*muugl*sqrtS + 
     -      dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu2**2*muu4*muugl*
     -       sqrtS + dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu1**2*muu2**2*muu4*
     -       muugl*sqrtS + dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu1**2*
     -       muu2**2*muu4*muugl*sqrtS + 
     -      dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu1**2*muu2**2*muu4*muugl*
     -       sqrtS - 3*dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu3**2*muu4*muugl*
     -       sqrtS - 3*dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu3**2*muu4*muugl*
     -       sqrtS - 3*dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu3**2*muu4*muugl*
     -       sqrtS - 3*dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu3**2*muu4*muugl*
     -       sqrtS + 3*dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu3**2*
     -       muu4*muugl*sqrtS + 
     -      3*dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu1**2*muu3**2*muu4*muugl*
     -       sqrtS + 3*dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu1**2*muu3**2*
     -       muu4*muugl*sqrtS + 
     -      3*dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu1**2*muu3**2*muu4*muugl*
     -       sqrtS + 2*dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu2**2*muu3**2*
     -       muu4*muugl*sqrtS + 
     -      2*dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu2**2*muu3**2*muu4*muugl*
     -       sqrtS + 2*dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu2**2*muu3**2*
     -       muu4*muugl*sqrtS + 
     -      2*dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu2**2*muu3**2*muu4*muugl*
     -       sqrtS + 2*dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu3**4*muu4*muugl*
     -       sqrtS + 2*dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu3**4*muu4*muugl*
     -       sqrtS + 2*dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu3**4*muu4*muugl*
     -       sqrtS + 2*dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu3**4*muu4*muugl*
     -       sqrtS - 3*dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu3*muu4**2*muugl*
     -       sqrtS - 3*dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu3*muu4**2*muugl*
     -       sqrtS - 3*dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu3*muu4**2*muugl*
     -       sqrtS - 3*dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu3*muu4**2*muugl*
     -       sqrtS + 2*dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu1**2*muu3*
     -       muu4**2*muugl*sqrtS + 
     -      2*dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu1**2*muu3*muu4**2*muugl*
     -       sqrtS + 2*dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3*
     -       muu4**2*muugl*sqrtS + 
     -      2*dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3*muu4**2*muugl*
     -       sqrtS + 3*dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu2**2*muu3*
     -       muu4**2*muugl*sqrtS + 
     -      3*dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu2**2*muu3*muu4**2*muugl*
     -       sqrtS + 3*dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3*
     -       muu4**2*muugl*sqrtS + 
     -      3*dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3*muu4**2*muugl*
     -       sqrtS + 3*dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu3**3*muu4**2*
     -       muugl*sqrtS + 3*dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu3**3*
     -       muu4**2*muugl*sqrtS + 
     -      3*dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu3**3*muu4**2*muugl*
     -       sqrtS + 3*dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu3**3*muu4**2*
     -       muugl*sqrtS - 2*dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu4**3*
     -       muugl*sqrtS - 2*dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu4**3*
     -       muugl*sqrtS - 2*dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu4**3*
     -       muugl*sqrtS - 2*dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu4**3*
     -       muugl*sqrtS + 2*dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*
     -       muu4**3*muugl*sqrtS + 
     -      2*dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4**3*muugl*
     -       sqrtS + 2*dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu1**2*muu4**3*
     -       muugl*sqrtS + 2*dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu1**2*
     -       muu4**3*muugl*sqrtS + 
     -      dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu2**2*muu4**3*muugl*sqrtS + 
     -      dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu2**2*muu4**3*muugl*sqrtS + 
     -      dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu2**2*muu4**3*muugl*sqrtS + 
     -      dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu2**2*muu4**3*muugl*sqrtS + 
     -      3*dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu3**2*muu4**3*muugl*
     -       sqrtS + 3*dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu3**2*muu4**3*
     -       muugl*sqrtS + 3*dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu3**2*
     -       muu4**3*muugl*sqrtS + 
     -      3*dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu3**2*muu4**3*muugl*
     -       sqrtS + 2*dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu3*muu4**4*muugl*
     -       sqrtS + 2*dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu3*muu4**4*muugl*
     -       sqrtS + 2*dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu3*muu4**4*muugl*
     -       sqrtS + 2*dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu3*muu4**4*muugl*
     -       sqrtS + dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu4**5*muugl*
     -       sqrtS + dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu4**5*muugl*
     -       sqrtS + dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu4**5*muugl*
     -       sqrtS + dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu4**5*muugl*
     -       sqrtS + dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu1**2*muu2**2*
     -       muugl**2*sqrtS + 
     -      dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu1**2*muu2**2*muugl**2*
     -       sqrtS - dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu3**2*muugl**2*
     -       sqrtS - dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu3**2*muugl**2*
     -       sqrtS + dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu2**2*muu3**2*
     -       muugl**2*sqrtS + 
     -      dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu2**2*muu3**2*muugl**2*
     -       sqrtS + dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu3**4*muugl**2*
     -       sqrtS + dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu3**4*muugl**2*
     -       sqrtS - 2*dPropFGlR*g1L*g2R*g3Rc*g4Lc*muu3*muu4*muugl**2*
     -       sqrtS - 2*dPropFGlL*g1R*g2L*g3Lc*g4Rc*muu3*muu4*muugl**2*
     -       sqrtS + 2*dPropFGlR*g1L*g2R*g3Rc*g4Lc*muu1**2*muu3*muu4*
     -       muugl**2*sqrtS + 
     -      2*dPropFGlL*g1R*g2L*g3Lc*g4Rc*muu1**2*muu3*muu4*muugl**2*
     -       sqrtS + 2*dPropFGlR*g1L*g2R*g3Rc*g4Lc*muu2**2*muu3*muu4*
     -       muugl**2*sqrtS + 
     -      2*dPropFGlL*g1R*g2L*g3Lc*g4Rc*muu2**2*muu3*muu4*muugl**2*
     -       sqrtS + 2*dPropFGlR*g1L*g2R*g3Rc*g4Lc*muu3**3*muu4*
     -       muugl**2*sqrtS + 
     -      2*dPropFGlL*g1R*g2L*g3Lc*g4Rc*muu3**3*muu4*muugl**2*
     -       sqrtS - dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu4**2*muugl**2*
     -       sqrtS - dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu4**2*muugl**2*
     -       sqrtS + dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu1**2*muu4**2*
     -       muugl**2*sqrtS + 
     -      dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu1**2*muu4**2*muugl**2*
     -       sqrtS + 3*dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu3**2*muu4**2*
     -       muugl**2*sqrtS + 
     -      3*dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu3**2*muu4**2*muugl**2*
     -       sqrtS + 2*dPropFGlR*g1L*g2R*g3Rc*g4Lc*muu3*muu4**3*
     -       muugl**2*sqrtS + 
     -      2*dPropFGlL*g1R*g2L*g3Lc*g4Rc*muu3*muu4**3*muugl**2*
     -       sqrtS + dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu4**4*muugl**2*
     -       sqrtS + dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu4**4*muugl**2*
     -       sqrtS - 2*dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu3*tred + 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3*tred + 
     -      2*dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3*tred + 
     -      2*dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu3**3*tred - 
     -      2*dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu4*tred + 
     -      2*dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4*tred + 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu2**2*muu4*tred + 
     -      3*dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu3**2*muu4*tred + 
     -      3*dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu3*muu4**2*tred + 
     -      2*dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu4**3*tred - 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muugl*tred - 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muugl*tred + 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu1**2*muugl*tred + 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muugl*tred + 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu2**2*muugl*tred + 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*muugl*tred + 
     -      2*dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu3**2*muugl*tred + 
     -      2*dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*muugl*tred + 
     -      2*dPropSGlL*g1L*g2R*g3Rc*g4Lc*muu3*muu4*muugl*tred + 
     -      2*dPropSGlL*g1R*g2L*g3Lc*g4Rc*muu3*muu4*muugl*tred + 
     -      2*dPropSGlL*g1R*g2L*g3Rc*g4Lc*muu4**2*muugl*tred + 
     -      2*dPropSGlL*g1L*g2R*g3Lc*g4Rc*muu4**2*muugl*tred + 
     -      dPropSGlL*g1L*g2L*g3Rc*g4Lc*muu3*muugl**2*tred + 
     -      dPropSGlL*g1L*g2L*g3Lc*g4Rc*muu4*muugl**2*tred - 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*sqrtS*tred - 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*sqrtS*tred + 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*sqrtS*tred + 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*sqrtS*tred - 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**4*sqrtS*tred - 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**4*sqrtS*tred + 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**2*sqrtS*tred + 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*sqrtS*tred - 
     -      3*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu2**2*sqrtS*
     -       tred - 3*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu2**2*
     -       sqrtS*tred - dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**4*sqrtS*
     -       tred - dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**4*sqrtS*tred + 
     -      4*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**2*sqrtS*tred + 
     -      4*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*sqrtS*tred - 
     -      3*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu3**2*sqrtS*
     -       tred - 3*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*muu3**2*
     -       sqrtS*tred - 4*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**2*
     -       muu3**2*sqrtS*tred - 
     -      4*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*muu3**2*sqrtS*
     -       tred - 3*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**4*sqrtS*
     -       tred - 3*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**4*sqrtS*
     -       tred + 4*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3*muu4*sqrtS*
     -       tred + 4*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3*muu4*sqrtS*
     -       tred - 4*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu1**2*muu3*muu4*
     -       sqrtS*tred - 4*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu1**2*muu3*
     -       muu4*sqrtS*tred - 
     -      4*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu2**2*muu3*muu4*sqrtS*
     -       tred - 4*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu2**2*muu3*muu4*
     -       sqrtS*tred - 4*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3**3*muu4*
     -       sqrtS*tred - 4*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3**3*muu4*
     -       sqrtS*tred + 4*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu4**2*sqrtS*
     -       tred + 4*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu4**2*sqrtS*
     -       tred - 4*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*muu4**2*
     -       sqrtS*tred - 4*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*
     -       muu4**2*sqrtS*tred - 
     -      3*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**2*muu4**2*sqrtS*
     -       tred - 3*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*muu4**2*
     -       sqrtS*tred - 7*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**2*
     -       muu4**2*sqrtS*tred - 
     -      7*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*muu4**2*sqrtS*
     -       tred - 4*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3*muu4**3*sqrtS*
     -       tred - 4*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3*muu4**3*sqrtS*
     -       tred - 3*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu4**4*sqrtS*
     -       tred - 3*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu4**4*sqrtS*
     -       tred + 2*dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu3*muugl*sqrtS*
     -       tred + 2*dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu3*muugl*sqrtS*
     -       tred + 2*dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu3*muugl*sqrtS*
     -       tred + 2*dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu3*muugl*sqrtS*
     -       tred - dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu1**2*muu3*muugl*
     -       sqrtS*tred - dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu1**2*muu3*
     -       muugl*sqrtS*tred - 
     -      dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3*muugl*sqrtS*
     -       tred - dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu1**2*muu3*muugl*
     -       sqrtS*tred - 2*dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu2**2*muu3*
     -       muugl*sqrtS*tred - 
     -      2*dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu2**2*muu3*muugl*sqrtS*
     -       tred - 2*dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3*muugl*
     -       sqrtS*tred - 2*dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu2**2*muu3*
     -       muugl*sqrtS*tred - 
     -      2*dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu3**3*muugl*sqrtS*tred - 
     -      2*dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu3**3*muugl*sqrtS*tred - 
     -      2*dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu3**3*muugl*sqrtS*tred - 
     -      2*dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu3**3*muugl*sqrtS*tred + 
     -      2*dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu4*muugl*sqrtS*tred + 
     -      2*dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu4*muugl*sqrtS*tred + 
     -      2*dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu4*muugl*sqrtS*tred + 
     -      2*dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu4*muugl*sqrtS*tred - 
     -      2*dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4*muugl*sqrtS*
     -       tred - 2*dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu1**2*muu4*muugl*
     -       sqrtS*tred - 2*dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu1**2*muu4*
     -       muugl*sqrtS*tred - 
     -      2*dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu1**2*muu4*muugl*sqrtS*
     -       tred - dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu2**2*muu4*muugl*
     -       sqrtS*tred - dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu2**2*muu4*
     -       muugl*sqrtS*tred - 
     -      dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu2**2*muu4*muugl*sqrtS*
     -       tred - dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu2**2*muu4*muugl*
     -       sqrtS*tred - 3*dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu3**2*muu4*
     -       muugl*sqrtS*tred - 
     -      3*dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu3**2*muu4*muugl*sqrtS*
     -       tred - 3*dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu3**2*muu4*muugl*
     -       sqrtS*tred - 3*dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu3**2*muu4*
     -       muugl*sqrtS*tred - 
     -      3*dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu3*muu4**2*muugl*sqrtS*
     -       tred - 3*dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu3*muu4**2*muugl*
     -       sqrtS*tred - 3*dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu3*muu4**2*
     -       muugl*sqrtS*tred - 
     -      3*dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu3*muu4**2*muugl*sqrtS*
     -       tred - 2*dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu4**3*muugl*sqrtS*
     -       tred - 2*dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu4**3*muugl*sqrtS*
     -       tred - 2*dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu4**3*muugl*sqrtS*
     -       tred - 2*dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu4**3*muugl*sqrtS*
     -       tred + dPropFGlL*g1R*g2L*g3Rc*g4Lc*muugl**2*sqrtS*tred + 
     -      dPropFGlR*g1L*g2R*g3Lc*g4Rc*muugl**2*sqrtS*tred - 
     -      dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu1**2*muugl**2*sqrtS*tred - 
     -      dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu1**2*muugl**2*sqrtS*tred - 
     -      dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu2**2*muugl**2*sqrtS*tred - 
     -      dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu2**2*muugl**2*sqrtS*tred - 
     -      2*dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu3**2*muugl**2*sqrtS*
     -       tred - 2*dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu3**2*muugl**2*
     -       sqrtS*tred - 2*dPropFGlR*g1L*g2R*g3Rc*g4Lc*muu3*muu4*
     -       muugl**2*sqrtS*tred - 
     -      2*dPropFGlL*g1R*g2L*g3Lc*g4Rc*muu3*muu4*muugl**2*sqrtS*
     -       tred - 2*dPropFGlL*g1R*g2L*g3Rc*g4Lc*muu4**2*muugl**2*
     -       sqrtS*tred - 2*dPropFGlR*g1L*g2R*g3Lc*g4Rc*muu4**2*
     -       muugl**2*sqrtS*tred - 
     -      dPropSGlL*g1R*g2R*g3Lc*g4Rc*muu3*tred**2 - 
     -      dPropSGlL*g1R*g2R*g3Rc*g4Lc*muu4*tred**2 - 
     -      dPropSGlL*g1R*g2L*g3Rc*g4Lc*muugl*tred**2 - 
     -      dPropSGlL*g1L*g2R*g3Lc*g4Rc*muugl*tred**2 - 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*sqrtS*tred**2 - 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*sqrtS*tred**2 + 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu1**2*sqrtS*tred**2 + 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu1**2*sqrtS*tred**2 + 
     -      2*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu2**2*sqrtS*tred**2 + 
     -      2*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu2**2*sqrtS*tred**2 + 
     -      3*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu3**2*sqrtS*tred**2 + 
     -      3*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu3**2*sqrtS*tred**2 + 
     -      2*dPropFGlL*g1L*g2R*g3Rc*g4Lc*muu3*muu4*sqrtS*tred**2 + 
     -      2*dPropFGlR*g1R*g2L*g3Lc*g4Rc*muu3*muu4*sqrtS*tred**2 + 
     -      3*dPropFGlR*g1R*g2L*g3Rc*g4Lc*muu4**2*sqrtS*tred**2 + 
     -      3*dPropFGlL*g1L*g2R*g3Lc*g4Rc*muu4**2*sqrtS*tred**2 + 
     -      dPropFGlL*g1L*g2L*g3Rc*g4Lc*muu3*muugl*sqrtS*tred**2 + 
     -      dPropFGlR*g1L*g2L*g3Rc*g4Lc*muu3*muugl*sqrtS*tred**2 + 
     -      dPropFGlL*g1R*g2R*g3Lc*g4Rc*muu3*muugl*sqrtS*tred**2 + 
     -      dPropFGlR*g1R*g2R*g3Lc*g4Rc*muu3*muugl*sqrtS*tred**2 + 
     -      dPropFGlL*g1R*g2R*g3Rc*g4Lc*muu4*muugl*sqrtS*tred**2 + 
     -      dPropFGlR*g1R*g2R*g3Rc*g4Lc*muu4*muugl*sqrtS*tred**2 + 
     -      dPropFGlL*g1L*g2L*g3Lc*g4Rc*muu4*muugl*sqrtS*tred**2 + 
     -      dPropFGlR*g1L*g2L*g3Lc*g4Rc*muu4*muugl*sqrtS*tred**2 + 
     -      dPropFGlL*g1R*g2L*g3Rc*g4Lc*muugl**2*sqrtS*tred**2 + 
     -      dPropFGlR*g1L*g2R*g3Lc*g4Rc*muugl**2*sqrtS*tred**2 - 
     -      dPropFGlR*g1R*g2L*g3Rc*g4Lc*sqrtS*tred**3 - 
     -      dPropFGlL*g1L*g2R*g3Lc*g4Rc*sqrtS*tred**3 - 
     -      dPropSGlR*(g1R*muugl*
     -          (2*g2L*g3Lc*g4Rc*muu3*muu4*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)
     -              + g2R*g3Rc*g4Lc*muu4*muugl*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g2R*g3Lc*g4Rc*muu3*muugl*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred) + 
     -            g2L*g3Rc*g4Lc*
     -             (muu3**4 - muu4**2 + muu4**4 + 
     -               muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -               muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -               muu2**2*tred - 2*muu4**2*tred + tred**2)) + 
     -         g1L*(g2L*(g3Lc*g4Rc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g3Rc*g4Lc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)
     -              + g2R*muugl*
     -             (2*g3Rc*g4Lc*muu3*muu4*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred) + 
     -               g3Lc*g4Rc*
     -                (muu3**4 - muu4**2 + muu4**4 + 
     -                  muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -                  muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -                  muu2**2*tred - 2*muu4**2*tred + tred**2)))) + 
     -      (-(dPropSGlR*(g1L*
     -               (g2R*muugl*
     -                  (g3Lc*g4Lc*muu4*
     -                     (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)
     -                      + g3Rc*g4Rc*muu3*
     -                     (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)
     -                    ) + 
     -                 g2L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -                    g3Lc*g4Lc*(-1 + muu3**2 + muu4**2))*
     -                  (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                    tred)) + 
     -              g1R*muugl*
     -               (2*g2R*g3Lc*g4Lc*muu3*muu4*muugl + 
     -                 g2R*g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -                 g2L*g3Rc*g4Rc*muu4*
     -                  (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -                 g2L*g3Lc*g4Lc*muu3*
     -                  (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))))
     -           - dPropSGlL*
     -          (g1R*(g2L*muugl*
     -                (g3Rc*g4Rc*muu4*
     -                   (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) +
     -                   g3Lc*g4Lc*muu3*
     -                   (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) +
     -                 g2R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -                  g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)) + 
     -            g1L*muugl*(2*g2L*g3Rc*g4Rc*muu3*muu4*muugl + 
     -               g2L*g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2R*g3Lc*g4Lc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g2R*g3Rc*g4Rc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))) + 
     -         sqrtS*(dPropFGlR*g1L*muugl*
     -             (g2R*muugl*
     -                (g3Lc*g4Lc*muu4*
     -                   (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) +
     -                    g3Rc*g4Rc*muu3*
     -                   (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) +
     -                 g2L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -                  g3Lc*g4Lc*(-1 + muu3**2 + muu4**2))*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)) + 
     -            dPropFGlL*g1R*muugl*
     -             (g2L*muugl*
     -                (g3Rc*g4Rc*muu4*
     -                   (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) +
     -                     g3Lc*g4Lc*muu3*
     -                   (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) +
     -                 g2R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -                  g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*
     -                (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                  tred)) + 
     -            dPropFGlR*g1R*
     -             (2*g2R*g3Lc*g4Lc*muu3*muu4*muugl + 
     -               g2R*g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2L*g3Rc*g4Rc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g2L*g3Lc*g4Lc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)
     -              + dPropFGlL*g1L*
     -             (2*g2L*g3Rc*g4Rc*muu3*muu4*muugl + 
     -               g2L*g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -               g2R*g3Lc*g4Lc*muu4*
     -                (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -               g2R*g3Rc*g4Rc*muu3*
     -                (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)
     -            ))*muupu(n)))/(DenGlU**2*dconjg(DenGauU(n))))	 
! 	  
      enddo
      
	  end
