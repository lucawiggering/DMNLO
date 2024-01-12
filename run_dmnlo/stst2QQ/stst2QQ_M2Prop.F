************************************************

	  function stst2QQ_M2Prop()
	  implicit none

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
#include "stst2QQ_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	  integer i,j,k,l,n,m,tcount,itypeU,itypeT,att,antitt,b,neu,genloop,typeloop
	  complex*16 g0,g1
	  complex*16 g1Lc, g1Rc, g2Lc, g2Rc, g3Rc, g3Lc, g4Lc, g4Rc
  	  complex*16 g1L, g1R, g2L, g2R, g3R, g3L, g4L, g4R
	  double precision Ml0,Ml1,tempCounter, tempH, tempZ, tempt, tempu
	  double precision tempH11, tempZ11, tempt11, tempu11, tempH12, tempZ12, tempt12, tempu12
	  double precision tempH21, tempZ21, tempt21, tempu21, tempH22, tempZ22, tempt22, tempu22
	  complex*16 PiSL, PiSR, PiFL, PiFR, PiGT, PiGB, testB, testT,tffAT,tffAL
	  double precision swap, symfac
	  double precision Mloop0,Mloop1,Mloop2

	  double precision Krondelta, ef, I3f, hf

	  double complex test, Gscalar

	  complex*16 BL1g,BL2g,BL3g,BR1g,BR2g,BR3g

	  double complex nc0,nc1,nc2,nc00,nc11,nc12,nc22

	  complex*16 DenSmV, DenSmH(6)
	  complex*16 DenSmV1, DenSmH1(6)
	  complex*16 tempGluonL, tempGluonR

* functions & subroutines
	  double precision colfac,colorf
	
	  complex*16 stst2qq_M2Prop
	
	  complex*16 B02,C03f
	  external myCget, myBget

	  stst2qq_M2Prop = 0d0
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

! 	
! !======================================================================================================================================
! ! ===============================================Gluino-Propagator Correction==========================================================
! !======================================================================================================================================

!--------------------------------------------Gluon-Gluino- Loop------------------------------------------------------------------------

	  PiSL = 0d0
	  PiSR = 0d0
	
	  PiFL = 0d0
	  PiFR = 0d0
	
	  Ml0 = 0d0
	  Ml1 = MGl
	
	  call myBget(tred*sqrtS**2,ML0**2,ML1**2)
	
	  g0 = GS
	  g1 = GS
	
	  colorf = Nc
	
	  PiSL = PiSL - g0*g1*colorf*Ml1/MGl*(4*xb0)
	  PiSR = PiSR - g0*g1*colorf*Ml1/MGl*(4*xb0)

	  PiFL = PiFL + g0*g1*colorf*(2*xb0+2*xb1)
	  PiFR = PiFR + g0*g1*colorf*(2*xb0+2*xb1)
	
! 	
! !--------------------------------------------Squark-Quark- Loop------------------------------------------------------------------------
	  do m=1,2
	  do typeloop=3,4
	    do genloop = 1,3

	    Ml0 = MSf(m,typeloop,genloop)
	    Ml1 = Mf(typeloop, genloop)
	
	    call myBget(tred*sqrtS**2,ML0**2,ML1**2)
	
	    g1R = -sqrt2*GS*USf(m,1,typeloop,genloop)
	    g1L = sqrt2*GS*USf(m,2,typeloop,genloop)
	    g2R = sqrt2*GS*USf(m,2,typeloop,genloop)
	    g2L = -sqrt2*GS*USf(m,1,typeloop,genloop)
	
	    colorf = 1d0/2d0
	
	    PiSL = PiSL + g1L*g2L*colorf*Ml1/MGl*xb0
	    PiSR = PiSR + g1R*g2R*colorf*Ml1/MGl*xb0
	  
	    PiFL = PiFL + g1L*g2R*colorf*(xb0+xb1)
	    PiFR = PiFR + g1R*g2L*colorf*(xb0+xb1)

	    
	    enddo
	  enddo
	  enddo
	
! 	
! ! !--------------------------------------------Squark-Antiquark- Loop------------------------------------------------------------------------
	  do m=1,2
	  do typeloop=3,4
	    do genloop = 1,3
	    Ml0 = MSf(m,typeloop,genloop)
	    Ml1 = Mf(typeloop, genloop)
	
	    call myBget(tred*sqrtS**2,ML0**2,ML1**2)
	
	    g1L = -sqrt2*GS*USf(m,1,typeloop,genloop)
	    g1R = sqrt2*GS*USf(m,2,typeloop,genloop)
	    g2L = sqrt2*GS*USf(m,2,typeloop,genloop)
	    g2R = -sqrt2*GS*USf(m,1,typeloop,genloop)
	
	    colorf = 1d0/2d0
	
	    PiSL = PiSL + g1R*g2R*colorf*Ml1/MGl*xb0
	    PiSR = PiSR + g1L*g2L*colorf*Ml1/MGl*xb0
	  
	    PiFL = PiFL - g1R*g2L*colorf*(xb0+xb1)
	    PiFR = PiFR - g1L*g2R*colorf*(xb0+xb1)
	    
        
	    enddo
	  enddo
	  enddo
! ! 
	  PiSL = PiSL/(16d0*pi**2)
	  PiSR = PiSR/(16d0*pi**2)
	
	  PiFL = PiFL/(16d0*pi**2)
	  PiFR = PiFR/(16d0*pi**2)


! ==========================================================t-channel correction=============================================================
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
	  
	  stst2qq_M2Prop = stst2qq_M2Prop + colfac*2d0*dreal((sqrtS**6*(g1R*(-(g1Lc*g2L*muu3*
     -            (g2Lc*muugl*(muu2**2 - muu4**2 - tred) - 
     -              2*g2Rc*muu4*tred)*
     -            (muugl**2*(PiFL + PiSL + PiSR) + PiFR*tred)) + 
     -         g1Rc*g2R*muugl*
     -          (muugl**2*PiSR + (PiFL + PiFR + PiSL)*tred)*
     -          (g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2Lc*muu4*(-muu1**2 + muu3**2 + tred)) + 
     -         g1Lc*g2R*muu3*muugl*
     -          (muugl**2*PiSR + (PiFL + PiFR + PiSL)*tred)*
     -          (2*g2Lc*muu4*muugl + g2Rc*(-muu2**2 + muu4**2 + tred))
     -           + g1Rc*g2L*(muugl**2*(PiFL + PiSL + PiSR) + 
     -            PiFR*tred)*
     -          (g2Rc*muu4*muugl*(-muu1**2 + muu3**2 + tred) + 
     -            g2Lc*((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                (muu2 + muu4) + tred - 
     -               (muu1**2 + muu2**2)*tred + tred**2))) + 
     -      g1L*(-(g1Rc*g2R*muu3*
     -            (g2Rc*muugl*(muu2**2 - muu4**2 - tred) - 
     -              2*g2Lc*muu4*tred)*
     -            (muugl**2*(PiFR + PiSL + PiSR) + PiFL*tred)) + 
     -         g1Lc*g2L*muugl*
     -          (muugl**2*PiSL + (PiFL + PiFR + PiSR)*tred)*
     -          (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2Rc*muu4*(-muu1**2 + muu3**2 + tred)) + 
     -         g1Rc*g2L*muu3*muugl*
     -          (muugl**2*PiSL + (PiFL + PiFR + PiSR)*tred)*
     -          (2*g2Rc*muu4*muugl + g2Lc*(-muu2**2 + muu4**2 + tred))
     -           + g1Lc*g2R*(muugl**2*(PiFR + PiSL + PiSR) + 
     -            PiFL*tred)*
     -          (g2Lc*muu4*muugl*(-muu1**2 + muu3**2 + tred) + 
     -            g2Rc*((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -                (muu2 + muu4) + tred - 
     -               (muu1**2 + muu2**2)*tred + tred**2)))))/
     -  (DenGlT**2*dconjg(DenGlT)))
     
! ! ---------------------t-channel and u-channel Gluino squared------------------------------
! 
! 
	  colfac = -1d0/2d0*CF
	
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(GluLt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Rc = dconjg(GluRt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(GluLt1(isf1,itt1,ftt2,igen1,fgen2))
	
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Rc = dconjg(GluRt1(isf1,itt1,ftt2,igen1,fgen2))
	  
	  
	  stst2qq_M2Prop = stst2qq_M2Prop + colfac*2d0*dreal((sqrtS**6*(g1R*(-(g1Lc*g2L*muu3*
     -            (g2Rc*muu4*
     -               (muu1**2 + muu2**2 - muu3**2 - muu4**2) + 
     -              g2Lc*muugl*(muu2**2 - muu4**2 - tred))*
     -            (muugl**2*(PiFL + PiSL + PiSR) + PiFR*tred)) + 
     -         g1Rc*g2R*muugl*
     -          (g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2Lc*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -               tred))*(muugl**2*PiSR + 
     -            (PiFL + PiFR + PiSL)*tred) + 
     -         g1Lc*g2R*muu3*muugl*
     -          (2*g2Lc*muu4*muugl + 
     -            g2Rc*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -          (muugl**2*PiSR + (PiFL + PiFR + PiSL)*tred) + 
     -         g1Rc*g2L*(muugl**2*(PiFL + PiSL + PiSR) + PiFR*tred)*
     -          (g2Rc*muu4*muugl*(-muu1**2 + muu3**2 + tred) + 
     -            g2Lc*((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred)))) + 
     -      g1L*(g1Lc*g2L*muugl*
     -          (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2Rc*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -               tred))*(muugl**2*PiSL + 
     -            (PiFL + PiFR + PiSR)*tred) + 
     -         g1Rc*g2Rc*muu3*muugl*
     -          (-(g2R*(muu2**2 - muu4**2 - tred)*
     -               (muugl**2*(PiFR + PiSL + PiSR) + PiFL*tred)) + 
     -            2*g2L*muu4*muugl*
     -             (muugl**2*PiSL + (PiFL + PiFR + PiSR)*tred)) + 
     -         g1Rc*g2Lc*muu3*
     -          (g2R*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2)*
     -             (muugl**2*(PiFR + PiSL + PiSR) + PiFL*tred) + 
     -            g2L*muugl*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -               tred)*(muugl**2*PiSL + (PiFL + PiFR + PiSR)*tred)
     -            ) + g1Lc*g2R*
     -          (muugl**2*(PiFR + PiSL + PiSR) + PiFL*tred)*
     -          (g2Lc*muu4*muugl*(-muu1**2 + muu3**2 + tred) + 
     -            g2Rc*((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -               muu1**2*(-muu2**2 + tred) + 
     -               muu3**2*(muu4**2 + tred))))))/
     -  (DenGlT**2*dconjg(DenGlU)))	  

! ! !---------------------t-channel Gluino t-channel Neutralino--------------------------------------------
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

	  stst2qq_M2Prop = stst2qq_M2Prop + colfac*2d0*dreal((sqrtS**6*(-(g1R*g2R*muugl*
     -         (g3Rc*g4Lc*muu4*(muu1**2 - muu3**2 - tred) + 
     -           g3Lc*g4Rc*muu3*(muu2**2 - muu4**2 - tred))*
     -         (muugl**2*PiSR + (PiFL + PiFR + PiSL)*tred)) - 
     -      g1L*g2L*muugl*(g3Lc*g4Rc*muu4*
     -          (muu1**2 - muu3**2 - tred) + 
     -         g3Rc*g4Lc*muu3*(muu2**2 - muu4**2 - tred))*
     -       (muugl**2*PiSL + (PiFL + PiFR + PiSR)*tred) + 
     -      g1R*g2L*(muugl**2*(PiFL + PiSL + PiSR) + PiFR*tred)*
     -       (2*g3Lc*g4Rc*muu3*muu4*tred + 
     -         g3Rc*g4Lc*((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -             (muu2 + muu4) + tred - (muu1**2 + muu2**2)*tred + 
     -            tred**2)) + 
     -      g1L*g2R*(muugl**2*(PiFR + PiSL + PiSR) + PiFL*tred)*
     -       (2*g3Rc*g4Lc*muu3*muu4*tred + 
     -         g3Lc*g4Rc*((muu1 - muu3)*(muu1 + muu3)*(muu2 - muu4)*
     -             (muu2 + muu4) + tred - (muu1**2 + muu2**2)*tred + 
     -            tred**2)) + 
     -      (-(g1L*g2R*(g3Lc*g4Lc*muu4*(muu1**2 - muu3**2 - tred) + 
     -              g3Rc*g4Rc*muu3*(muu2**2 - muu4**2 - tred))*
     -            (muugl**2*(PiFR + PiSL + PiSR) + PiFL*tred)) - 
     -         g1R*g2L*(g3Rc*g4Rc*muu4*(muu1**2 - muu3**2 - tred) + 
     -            g3Lc*g4Lc*muu3*(muu2**2 - muu4**2 - tred))*
     -          (muugl**2*(PiFL + PiSL + PiSR) + PiFR*tred) + 
     -         g1R*g2R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -            g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*muugl*
     -          (muugl**2*PiSR + (PiFL + PiFR + PiSL)*tred) + 
     -         g1L*g2L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -            g3Lc*g4Lc*(-1 + muu3**2 + muu4**2))*muugl*
     -          (muugl**2*PiSL + (PiFL + PiFR + PiSR)*tred))*muupt(n))
     -    )/(DenGlT**2*dconjg(DenGauT(n))))	  
	  
	  enddo
	  
	  
!---------------------t-channel Gluino u-channel Neutralino--------------------------------------------------------
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
	  
	  stst2qq_M2Prop = stst2qq_M2Prop + colfac*2d0*dreal( (sqrtS**6*(g1R*g2R*muugl*
     -       (g3Rc*g4Lc*muu4*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -         g3Lc*g4Rc*muu3*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -       (muugl**2*PiSR + (PiFL + PiFR + PiSL)*tred) + 
     -      g1L*g2L*muugl*(g3Lc*g4Rc*muu4*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -         g3Rc*g4Lc*muu3*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -       (muugl**2*PiSL + (PiFL + PiFR + PiSR)*tred) - 
     -      g1R*g2L*(muugl**2*(PiFL + PiSL + PiSR) + PiFR*tred)*
     -       (g3Lc*g4Rc*muu3*muu4*
     -          (muu1**2 + muu2**2 - muu3**2 - muu4**2) + 
     -         g3Rc*g4Lc*(muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -            (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred + 
     -            tred**2)) - 
     -      g1L*g2R*(muugl**2*(PiFR + PiSL + PiSR) + PiFL*tred)*
     -       (g3Rc*g4Lc*muu3*muu4*
     -          (muu1**2 + muu2**2 - muu3**2 - muu4**2) + 
     -         g3Lc*g4Rc*(muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -            (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred + 
     -            tred**2)) + 
     -      (-(g1L*g2R*(g3Lc*g4Lc*muu4*(muu1**2 - muu3**2 - tred) + 
     -              g3Rc*g4Rc*muu3*(muu2**2 - muu4**2 - tred))*
     -            (muugl**2*(PiFR + PiSL + PiSR) + PiFL*tred)) - 
     -         g1R*g2L*(g3Rc*g4Rc*muu4*(muu1**2 - muu3**2 - tred) + 
     -            g3Lc*g4Lc*muu3*(muu2**2 - muu4**2 - tred))*
     -          (muugl**2*(PiFL + PiSL + PiSR) + PiFR*tred) + 
     -         g1R*g2R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -            g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*muugl*
     -          (muugl**2*PiSR + (PiFL + PiFR + PiSL)*tred) + 
     -         g1L*g2L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -            g3Lc*g4Lc*(-1 + muu3**2 + muu4**2))*muugl*
     -          (muugl**2*PiSL + (PiFL + PiFR + PiSR)*tred))*muupu(n))
     -    )/(DenGlT**2*dconjg(DenGauU(n))))	  	  
	  
	  enddo
	  
! 	  
! 	  ! ==========================================================u-channel correction=============================================================
! 	  
!--------------------------------------------Gluon-Gluino- Loop------------------------------------------------------------------------

	  PiSL = 0d0
	  PiSR = 0d0
	
	  PiFL = 0d0
	  PiFR = 0d0
	
	  Ml0 = 0d0
	  Ml1 = MGl
	
	  call myBget(ured*sqrtS**2,ML0**2,ML1**2)
	
	  g0 = GS
	  g1 = GS
	
	  colorf = Nc
	
	  PiSL = PiSL - g0*g1*colorf*Ml1/MGl*(4*xb0)
	  PiSR = PiSR - g0*g1*colorf*Ml1/MGl*(4*xb0)

	  PiFL = PiFL + g0*g1*colorf*(2*xb0+2*xb1)
	  PiFR = PiFR + g0*g1*colorf*(2*xb0+2*xb1)
	
! 	
! !--------------------------------------------Squark-Quark- Loop------------------------------------------------------------------------
	  do m=1,2
	  do typeloop=3,4
	    do genloop = 1,3

	    Ml0 = MSf(m,typeloop,genloop)
	    Ml1 = Mf(typeloop, genloop)
	
	    call myBget(ured*sqrtS**2,ML0**2,ML1**2)
	
	    g1R = -sqrt2*GS*USf(m,1,typeloop,genloop)
	    g1L = sqrt2*GS*USf(m,2,typeloop,genloop)
	    g2R = sqrt2*GS*USf(m,2,typeloop,genloop)
	    g2L = -sqrt2*GS*USf(m,1,typeloop,genloop)
	
	    colorf = 1d0/2d0
	
	    PiSL = PiSL + g1L*g2L*colorf*Ml1/MGl*xb0
	    PiSR = PiSR + g1R*g2R*colorf*Ml1/MGl*xb0
	  
	    PiFL = PiFL + g1L*g2R*colorf*(xb0+xb1)
	    PiFR = PiFR + g1R*g2L*colorf*(xb0+xb1)

	    
	    enddo
	  enddo
	  enddo
	
! 	
! ! !--------------------------------------------Squark-Antiquark- Loop------------------------------------------------------------------------
	  do m=1,2
	  do typeloop=3,4
	    do genloop = 1,3
	    Ml0 = MSf(m,typeloop,genloop)
	    Ml1 = Mf(typeloop, genloop)
	
	    call myBget(ured*sqrtS**2,ML0**2,ML1**2)
	
	    g1L = -sqrt2*GS*USf(m,1,typeloop,genloop)
	    g1R = sqrt2*GS*USf(m,2,typeloop,genloop)
	    g2L = sqrt2*GS*USf(m,2,typeloop,genloop)
	    g2R = -sqrt2*GS*USf(m,1,typeloop,genloop)
	
	    colorf = 1d0/2d0
	
	    PiSL = PiSL + g1R*g2R*colorf*Ml1/MGl*xb0
	    PiSR = PiSR + g1L*g2L*colorf*Ml1/MGl*xb0
	  
	    PiFL = PiFL - g1R*g2L*colorf*(xb0+xb1)
	    PiFR = PiFR - g1L*g2R*colorf*(xb0+xb1)
	    
        
	    enddo
	  enddo
	  enddo
! ! 
	  PiSL = PiSL/(16d0*pi**2)
	  PiSR = PiSR/(16d0*pi**2)
	
	  PiFL = PiFL/(16d0*pi**2)
	  PiFR = PiFR/(16d0*pi**2)


! 
! ! !---------------------u-channel Gluino squared------------------------------
! ! ! 
	  colfac = Nc/2d0*CF
    
	
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1Lc = dconjg(g1L)
	  g1Rc = dconjg(g1R)
	
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2Lc = dconjg(g2L)
	  g2Rc = dconjg(g2R)

	  stst2qq_M2Prop = stst2qq_M2Prop + colfac*2d0*dreal((sqrtS**6*(g1L*(g1Lc*g2L*muugl*
     -          (-PiFR + muugl**2*PiSL - PiSR + 
     -            PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred) + 
     -            (PiFR + PiSR)*
     -             (muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred))*
     -          (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2Rc*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -               tred)) + 
     -         g1Rc*g2Rc*muu3*muugl*
     -          (2*g2L*muu4*muugl*
     -             (-PiFR + muugl**2*PiSL - PiSR + 
     -               PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 - tred) + 
     -               (PiFR + PiSR)*
     -                (muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)) +
     -             g2R*(muugl**2*(PiFR + PiSL + PiSR) + 
     -               PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 - tred))*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -         g1Rc*g2Lc*muu3*
     -          (2*g2R*muu4*(muugl**2*(PiFR + PiSL + PiSR) + 
     -               PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 - tred))*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)
     -              + g2L*muugl*
     -             (-PiFR + muugl**2*PiSL - PiSR + 
     -               PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 - tred) + 
     -               (PiFR + PiSR)*
     -                (muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred))*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -         g1Lc*g2R*(muugl**2*(PiFR + PiSL + PiSR) + 
     -            PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred))*(g2Lc*muu4*muugl*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g2Rc*(muu3**4 - muu4**2 + 
     -               muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -               (muu4**2 - tred)**2 + 
     -               muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -               muu2**2*tred))) + 
     -      g1R*(g1Rc*g2Rc*muugl*
     -          (g2R*(-1 + muu3**2 + muu4**2)*muugl*
     -             (-PiFR - PiSL + muugl**2*PiSR + 
     -               PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 - tred) + 
     -               (PiFR + PiSL)*
     -                (muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)) + 
     -              g2L*muu4*
     -             (muugl**2*(PiFL + PiSL + PiSR) + 
     -               PiFR*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 - tred))*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)) + 
     -         g1Lc*g2R*muu3*muugl*
     -          (-PiFR - PiSL + muugl**2*PiSR + 
     -            PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred) + 
     -            (PiFR + PiSL)*
     -             (muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred))*
     -          (2*g2Lc*muu4*muugl + 
     -            g2Rc*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -         g1Lc*g2L*muu3*
     -          (muugl**2*(PiFL + PiSL + PiSR) + 
     -            PiFR*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred))*(2*g2Rc*muu4*
     -             (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)
     -              + g2Lc*muugl*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -         g1Rc*g2Lc*(g2R*muu4*muugl*
     -             (-PiFR - PiSL + muugl**2*PiSR + 
     -               PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 - tred) + 
     -               (PiFR + PiSL)*
     -                (muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred))*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g2L*(muugl**2*(PiFL + PiSL + PiSR) + 
     -               PiFR*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 - tred))*
     -             (muu3**4 - muu4**2 + 
     -               muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -               (muu4**2 - tred)**2 + 
     -               muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -               muu2**2*tred)))))/(DenGlU**2*dconjg(DenGlU)))	  	  
       
!      
! ! ---------------------t-channel and u-channel Gluino squared------------------------------
! 
! 
	  colfac = -1d0/2d0*CF
	
	  g1Lc = dconjg(GluLt1(isf1,itt1,ftt1,igen1,fgen1))
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	
	  g1Rc = dconjg(GluRt1(isf1,itt1,ftt1,igen1,fgen1))
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	
	  g2Lc = dconjg(GluLt1(isf2,itt2,ftt2,igen2,fgen2))
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	
	  g2Rc = dconjg(GluRt1(isf2,itt2,ftt2,igen2,fgen2))
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  
	  stst2qq_M2Prop = stst2qq_M2Prop + colfac*2d0*dreal((sqrtS**6*(g1L*(g1Rc*g2Lc*muu3*
     -          (g2R*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2)*
     -             (muugl**2*(PiFR + PiSL + PiSR) + 
     -               PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 - tred)) - 
     -            g2L*muugl*(-PiFR + muugl**2*PiSL - PiSR + 
     -               PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 - tred) + 
     -               (PiFR + PiSR)*
     -                (muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred))*
     -             (muu2**2 - muu4**2 - tred)) + 
     -         g1Rc*g2Rc*muu3*muugl*
     -          (2*g2L*muu4*muugl*
     -             (-PiFR + muugl**2*PiSL - PiSR + 
     -               PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 - tred) + 
     -               (PiFR + PiSR)*
     -                (muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)) +
     -              g2R*(muugl**2*(PiFR + PiSL + PiSR) + 
     -               PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 - tred))*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -         g1Lc*g2L*muugl*
     -          (-PiFR + muugl**2*PiSL - PiSR + 
     -            PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred) + 
     -            (PiFR + PiSR)*
     -             (muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred))*
     -          (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2Rc*muu4*(-muu1**2 + muu3**2 + tred)) - 
     -         g1Lc*g2R*(muugl**2*(PiFR + PiSL + PiSR) + 
     -            PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred))*(-(g2Lc*muu4*muugl*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)) + 
     -            g2Rc*(muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -               (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred + 
     -               tred**2))) + 
     -      g1R*(g1Rc*g2Rc*muugl*
     -          (g2R*(-1 + muu3**2 + muu4**2)*muugl*
     -             (-PiFR - PiSL + muugl**2*PiSR + 
     -               PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 - tred) + 
     -               (PiFR + PiSL)*
     -                (muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)) +
     -              g2L*muu4*
     -             (muugl**2*(PiFL + PiSL + PiSR) + 
     -               PiFR*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 - tred))*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)) - 
     -         g1Lc*g2L*muu3*
     -          (muugl**2*(PiFL + PiSL + PiSR) + 
     -            PiFR*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred))*(g2Rc*muu4*
     -             (muu1**2 + muu2**2 - muu3**2 - muu4**2) - 
     -            g2Lc*muugl*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -         g1Lc*g2R*muu3*muugl*
     -          (-PiFR - PiSL + muugl**2*PiSR + 
     -            PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred) + 
     -            (PiFR + PiSL)*
     -             (muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred))*
     -          (2*g2Lc*muu4*muugl + g2Rc*(-muu2**2 + muu4**2 + tred))
     -           + g1Rc*g2Lc*
     -          (-(g2R*muu4*muugl*
     -               (-PiFR - PiSL + muugl**2*PiSR + 
     -                 PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                    muu4**2 - tred) + 
     -                 (PiFR + PiSL)*
     -                  (muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred)
     -                 )*(muu1**2 - muu3**2 - tred)) - 
     -            g2L*(muugl**2*(PiFL + PiSL + PiSR) + 
     -               PiFR*(-1 + muu1**2 + muu2**2 + muu3**2 + 
     -                  muu4**2 - tred))*
     -             (muu1**2*muu2**2 - muu3**2*muu4**2 + tred - 
     -               (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred + 
     -               tred**2)))))/(DenGlU**2*dconjg(DenGlT)))	  
         
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
	  
	  
	  stst2qq_M2Prop = stst2qq_M2Prop + colfac*2d0*dreal( (sqrtS**6*(-(g1R*g2R*muugl*
     -         (g3Rc*g4Lc*muu4*(muu1**2 - muu3**2 - tred) + 
     -           g3Lc*g4Rc*muu3*(muu2**2 - muu4**2 - tred))*
     -         (-PiSL + muu1**2*PiSL + muu2**2*PiSL + muu3**2*PiSL + 
     -           muu4**2*PiSL + muugl**2*PiSR + 
     -           PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -              tred) + PiFR*
     -            (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred) -
     -            PiSL*tred)) - 
     -      g1L*g2L*muugl*(g3Lc*g4Rc*muu4*
     -          (muu1**2 - muu3**2 - tred) + 
     -         g3Rc*g4Lc*muu3*(muu2**2 - muu4**2 - tred))*
     -       (muugl**2*PiSL - PiSR + muu1**2*PiSR + muu2**2*PiSR + 
     -         muu3**2*PiSR + muu4**2*PiSR + 
     -         PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -            tred) + PiFR*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred) - 
     -         PiSR*tred) - g1R*g2L*
     -       (muugl**2*(PiFL + PiSL + PiSR) + 
     -         PiFR*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -            tred))*(g3Lc*g4Rc*muu3*muu4*
     -          (muu1**2 + muu2**2 - muu3**2 - muu4**2) + 
     -         g3Rc*g4Lc*(muu1**2*(muu2**2 - tred) + 
     -            tred*(1 - muu2**2 - muu4**2 + tred) - 
     -            muu3**2*(muu4**2 + tred))) - 
     -      g1L*g2R*(muugl**2*(PiFR + PiSL + PiSR) + 
     -         PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -            tred))*(g3Rc*g4Lc*muu3*muu4*
     -          (muu1**2 + muu2**2 - muu3**2 - muu4**2) + 
     -         g3Lc*g4Rc*(muu1**2*(muu2**2 - tred) + 
     -            tred*(1 - muu2**2 - muu4**2 + tred) - 
     -            muu3**2*(muu4**2 + tred))) + 
     -      (g1R*g2L*(muugl**2*(PiFL + PiSL + PiSR) + 
     -            PiFR*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred))*(g3Rc*g4Rc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g3Lc*g4Lc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -         g1L*g2R*(muugl**2*(PiFR + PiSL + PiSR) + 
     -            PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred))*(g3Lc*g4Lc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g3Rc*g4Rc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -         g1R*g2R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -            g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*muugl*
     -          (-PiSL + muu1**2*PiSL + muu2**2*PiSL + muu3**2*PiSL + 
     -            muu4**2*PiSL + muugl**2*PiSR + 
     -            PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred) + 
     -            PiFR*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred) - PiSL*tred) + 
     -         g1L*g2L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -            g3Lc*g4Lc*(-1 + muu3**2 + muu4**2))*muugl*
     -          (muugl**2*PiSL - PiSR + muu1**2*PiSR + muu2**2*PiSR + 
     -            muu3**2*PiSR + muu4**2*PiSR + 
     -            PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred) + 
     -            PiFR*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred) - PiSR*tred))*muupt(n)))/
     -  (DenGlU**2*dconjg(DenGauT(n))))	  
	  
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
	  
	  stst2qq_M2Prop = stst2qq_M2Prop + colfac*2d0*dreal((sqrtS**6*(g1R*g2R*muugl*
     -       (g3Rc*g4Lc*muu4*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -         g3Lc*g4Rc*muu3*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -       (-PiSL + muu1**2*PiSL + muu2**2*PiSL + muu3**2*PiSL + 
     -         muu4**2*PiSL + muugl**2*PiSR + 
     -         PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -            tred) + PiFR*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred) - 
     -         PiSL*tred) + g1L*g2L*muugl*
     -       (g3Lc*g4Rc*muu4*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -         g3Rc*g4Lc*muu3*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred))*
     -       (muugl**2*PiSL - PiSR + muu1**2*PiSR + muu2**2*PiSR + 
     -         muu3**2*PiSR + muu4**2*PiSR + 
     -         PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -            tred) + PiFR*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred) - 
     -         PiSR*tred) + g1R*g2L*
     -       (muugl**2*(PiFL + PiSL + PiSR) + 
     -         PiFR*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -            tred))*(2*g3Lc*g4Rc*muu3*muu4*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred) + 
     -         g3Rc*g4Lc*(muu3**4 - muu4**2 + muu4**4 + 
     -            muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -            muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -            muu2**2*tred - 2*muu4**2*tred + tred**2)) + 
     -      g1L*g2R*(muugl**2*(PiFR + PiSL + PiSR) + 
     -         PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -            tred))*(2*g3Rc*g4Lc*muu3*muu4*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred) + 
     -         g3Lc*g4Rc*(muu3**4 - muu4**2 + muu4**4 + 
     -            muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -            muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -            muu2**2*tred - 2*muu4**2*tred + tred**2)) + 
     -      (g1R*g2L*(muugl**2*(PiFL + PiSL + PiSR) + 
     -            PiFR*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred))*(g3Rc*g4Rc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g3Lc*g4Lc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -         g1L*g2R*(muugl**2*(PiFR + PiSL + PiSR) + 
     -            PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred))*(g3Lc*g4Lc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g3Rc*g4Rc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -         g1R*g2R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -            g3Rc*g4Rc*(-1 + muu3**2 + muu4**2))*muugl*
     -          (-PiSL + muu1**2*PiSL + muu2**2*PiSL + muu3**2*PiSL + 
     -            muu4**2*PiSL + muugl**2*PiSR + 
     -            PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred) + 
     -            PiFR*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred) - PiSL*tred) + 
     -         g1L*g2L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -            g3Lc*g4Lc*(-1 + muu3**2 + muu4**2))*muugl*
     -          (muugl**2*PiSL - PiSR + muu1**2*PiSR + muu2**2*PiSR + 
     -            muu3**2*PiSR + muu4**2*PiSR + 
     -            PiFL*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred) + 
     -            PiFR*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -               tred) - PiSR*tred))*muupu(n)))/
     -  (DenGlU**2*dconjg(DenGauU(n))))	  
	  
	  enddo
	 
	  end
