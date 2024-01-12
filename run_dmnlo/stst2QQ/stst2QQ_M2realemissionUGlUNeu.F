************************************************

	  function stst2QQ_M2realemissionUGLUNeu()
	  IMPLICIT NONE

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
#include "stst2QQ_GenCouplings.h"

 	  double precision findValW
      integer i,j,k,l,m,n,neuc

	  double precision colfac
  	  double precision muuq,muuqc
  	  double precision MGaugino(6), wGaugino(6)
	  double complex g1L, g1Lc, g2L, g2Lc, g3L, g3Lc, g4L, g4Lc
	  double complex g1R, g1Rc, g2R, g2Rc, g3R, g3Rc, g4R, g4Rc
	  double complex gA, gB, gD, gE, gF
	  double complex gAc, gBc, gDc, gEc, gFc
	  double complex test
	  
	  complex*16 Denlegp1,Denlegp1c,Denlegp2,Denlegp2c,Denlegk1,Denlegk1c,Denlegk2,Denlegk2c
	  complex*16 DenU1, DenU2,DenU1c, DenU2c


	  complex*16 stst2QQ_M2realemissionUGLUNeu
	
	  stst2QQ_M2realemissionUGLUNeu = 0d0
	  test = 0d0
      
c 	Propagator--------------------------------------------------------------------------------------

	  Denlegp1 = -2d0*p1k3
	  Denlegp1c = dconjg(Denlegp1)
	  Denlegp2 = -2d0*p2k3
	  Denlegp2c = dconjg(Denlegp2)
	  Denlegk1 = 2d0*k1k3
	  Denlegk1c = dconjg(Denlegk1)
	  Denlegk2 = 2d0*k2k3
	  Denlegk2c = dconjg(Denlegk2)
	  
	  
      do i=1,4
	  MGaugino(i) = MNeu(i)
	  wGaugino(i) = wNeu(i)
	  enddo
	  MGaugino(5) = MCha(1)
	  wGaugino(5) = wCha(1)
	  MGaugino(6) = MCha(2)
	  wGaugino(6) = wCha(2) 
	  
c-------------------------------------------------------------------------------------------------------
!     print*,'--------------------------Tree-level-----------------------------------'
	  gA = -GS
	  gAc = -GS
	
	  gB = -GS
	  gBc = -GS
	
	  gD = -GS
	  gDc = -GS
	
	  gE = -GS
	  gEc = -GS
	
	  gF = GS
	  gFc = GS
	  
	  do neuc = 1,6
! ! ---------------------U-channel Gluino squared------------------------------
! 
      DenU1 = muu4**2*sqrtS**2+muu1**2*sqrtS**2-2*p1k2-MGl**2
	  DenU1c = muu4**2*sqrtS**2+muu1**2*sqrtS**2-2*p1k2-MGaugino(neuc)**2!- cI*MGaugino(neu)*wGaugino(neu)
	  DenU2 = muu3**2*sqrtS**2+muu2**2*sqrtS**2-2*p2k1-MGl**2
	  DenU2c = muu3**2*sqrtS**2+muu2**2*sqrtS**2-2*p2k1-MGaugino(neuc)**2!- cI*MGaugino(neu)*wGaugino(neu)
	  
!-------------Masses------------------------------------------------------------------------------------------------------------------------------------
	  muuq = MGl/sqrtS
	  muuqc = MGaugino(neuc)/sqrtS
	  
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1Lc = dconjg(ALt1(neuc,isf2,itt2,ftt1,igen2,fgen1))
	  g1Rc = dconjg(ARt1(neuc,isf2,itt2,ftt1,igen2,fgen1))
	
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2Lc = dconjg(ALt1(neuc,isf1,itt1,ftt2,igen1,fgen2))
	  g2Rc = dconjg(ARt1(neuc,isf1,itt1,ftt2,igen1,fgen2))

	  
! ! ---------------------p1-p1------------------------------
! 
	  colfac = 0d0
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(-8*gA*gAc*(p1k3 - muu1**2*sqrtS**2)*
     -    (-2*(g1R*g1Rc*g2L*g2Lc + g1L*g1Lc*g2R*g2Rc)*p2k1*p2k2 + 
     -      (g1Lc*g1R*muu3*(-2*g2L*g2Rc*muu4*p2k1 + g2R*g2Rc*muuq*(-k1k2 + p2k2) + 
     -            g2L*g2Lc*muuqc*(-k1k2 + p2k2)) + 
     -         g1L*g1Rc*muu3*(-2*g2Lc*g2R*muu4*p2k1 + g2L*g2Lc*muuq*(-k1k2 + p2k2) + 
     -            g2R*g2Rc*muuqc*(-k1k2 + p2k2)) + 
     -         g1R*g1Rc*(-(g2R*muuq*(g2Rc*k1k2*muuqc + g2Lc*muu4*p2k1)) + 
     -            g2L*(g2Lc*k1k2*(muu2 - muu3)*(muu2 + muu3) - g2Rc*muu4*muuqc*p2k1 + 
     -               2*g2Lc*muu3**2*p2k2)) + 
     -         g1L*g1Lc*(-(g2L*muuq*(g2Lc*k1k2*muuqc + g2Rc*muu4*p2k1)) + 
     -            g2R*(g2Rc*k1k2*(muu2 - muu3)*(muu2 + muu3) - g2Lc*muu4*muuqc*p2k1 + 
     -               2*g2Rc*muu3**2*p2k2)))*sqrtS**2 + 
     -      muu3*muu4*(g1L*(g1Rc*g2Lc*g2R*(muu2**2 + muu3**2) + g1Rc*g2L*g2Rc*muuq*muuqc + 
     -            g1Lc*muu3*(g2L*g2Rc*muuq + g2Lc*g2R*muuqc)) + 
     -         g1R*(g1Lc*g2L*g2Rc*(muu2**2 + muu3**2) + g1Lc*g2Lc*g2R*muuq*muuqc + 
     -            g1Rc*muu3*(g2Lc*g2R*muuq + g2L*g2Rc*muuqc)))*sqrtS**4))/
     -  (Denlegp1*Denlegp1c*DenU1*DenU1c)
     
! ---------------------p1-p2------------------------------

	  colfac = Nc/2d0*CF
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(4*gA*gBc*(p1k3 - 2*p1p2 + p2k3)*
     -    ((g1R*g1Rc*g2L*g2Lc + g1L*g1Lc*g2R*g2Rc)*(k1k2*p1p2 - p1k2*p2k1 - p1k1*p2k2) + 
     -      (g1L*g1Lc*(g2L*muuq*(g2Lc*k1k2*muuqc + g2Rc*muu4*(k1k2 - p1k1)) + 
     -            g2R*g2Rc*muu3**2*p1k2 + g2R*muu4*(g2Rc*muu4 + g2Lc*muuqc)*p2k1) + 
     -         g1R*g1Rc*(g2L*g2Lc*muu3**2*p1k2 + g2Rc*muuqc*(g2R*k1k2*muuq + g2L*muu4*p2k1) + 
     -            g2Lc*muu4*(g2R*muuq*(k1k2 - p1k1) + g2L*muu4*p2k1)) + 
     -         g1L*g1Rc*muu3*(g2L*g2Lc*muuq*p1k2 + g2R*g2Rc*muuqc*(k1k2 - p2k2) + 
     -            g2Lc*g2R*muu4*(k1k2 - p1k1 + p1p2 - p2k2)) + 
     -         g1Lc*g1R*muu3*(g2R*g2Rc*muuq*p1k2 + g2L*g2Lc*muuqc*(k1k2 - p2k2) + 
     -            g2L*g2Rc*muu4*(k1k2 - p1k1 + p1p2 - p2k2)))*sqrtS**2 - 
     -      muu3*muu4*(g1L*g1Lc*g2R*muu3*(g2Rc*muu4 + g2Lc*muuqc) + 
     -         g1Lc*g1R*g2R*muuq*(g2Rc*muu4 + g2Lc*muuqc) + 
     -         g1R*g1Rc*g2L*muu3*(g2Lc*muu4 + g2Rc*muuqc) + 
     -         g1L*g1Rc*g2L*muuq*(g2Lc*muu4 + g2Rc*muuqc))*sqrtS**4))/
     -  (Denlegp1*Denlegp2c*DenU1*DenU2c)
     
! ---------------------p2-p1-------------------------------

	  colfac = Nc/2d0*CF
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(4*gAc*gB*(p1k3 - 2*p1p2 + p2k3)*
     -    ((g1R*g1Rc*g2L*g2Lc + g1L*g1Lc*g2R*g2Rc)*(k1k2*p1p2 - p1k2*p2k1 - p1k1*p2k2) + 
     -      (g1L*g1Lc*(g2Lc*muuqc*(g2L*k1k2*muuq + g2R*muu4*(k1k2 - p1k1)) + 
     -            g2R*g2Rc*muu3**2*p1k2 + g2Rc*muu4*(g2R*muu4 + g2L*muuq)*p2k1) + 
     -         g1R*g1Rc*(g2L*g2Lc*muu3**2*p1k2 + g2R*muuq*(g2Rc*k1k2*muuqc + g2Lc*muu4*p2k1) + 
     -            g2L*muu4*(g2Rc*muuqc*(k1k2 - p1k1) + g2Lc*muu4*p2k1)) + 
     -         g1L*g1Rc*muu3*(g2R*g2Rc*muuqc*p1k2 + g2L*g2Lc*muuq*(k1k2 - p2k2) + 
     -            g2Lc*g2R*muu4*(k1k2 - p1k1 + p1p2 - p2k2)) + 
     -         g1Lc*g1R*muu3*(g2L*g2Lc*muuqc*p1k2 + g2R*g2Rc*muuq*(k1k2 - p2k2) + 
     -            g2L*g2Rc*muu4*(k1k2 - p1k1 + p1p2 - p2k2)))*sqrtS**2 - 
     -      muu3*muu4*(g1R*g2Lc*(g2L*muu4 + g2R*muuq)*(g1Rc*muu3 + g1Lc*muuqc) + 
     -         g1L*g2Rc*(g2R*muu4 + g2L*muuq)*(g1Lc*muu3 + g1Rc*muuqc))*sqrtS**4))/
     -  (Denlegp1c*Denlegp2*DenU1c*DenU2)     
	  
! ! ---------------------p1-k1------------------------------

	  colfac = Nc/2d0*CF
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(4*gA*gDc*(g1L*g1Rc*muu3*sqrtS**2*
     -       (g2R*g2Rc*muuqc*(k1k3*p1k2 + k1k2*(k1k3 - 2*p1k1 - p1k3) + k2k3*(-p1k1 + p1p2) - 
     -            k1k3*p2k2 + 2*p1k1*p2k2 + p1k3*p2k2 - p1k2*p2k3) + 
     -         g2L*g2Rc*muu4*muuq*muuqc*(-k1k3 + 2*p1k1 + p1k3)*sqrtS**2 + 
     -         g2L*g2Lc*muuq*(p1k2*(k1k3 - 2*(p1k1 + p1k3)) + 
     -            (k2k3*muu1**2 + muu4**2*(-k1k3 + 2*p1k1 + p1k3))*sqrtS**2) + 
     -         g2Lc*g2R*muu4*(k1k2*(k1k3 - 2*p1k1 - p1k3) - 
     -            (k2k3 - 2*(p1k1 + p1k3))*(p1k1 - p1p2) + (2*p1k1 + p1k3)*p2k2 - p1k2*p2k3 + 
     -            muu1**2*p2k3*sqrtS**2 - k1k3*(p1k1 - p1k2 - p1p2 + p2k2 + muu1**2*sqrtS**2)))+
     -       g1Lc*g1R*muu3*sqrtS**2*(g2L*g2Lc*muuqc*
     -          (k1k3*p1k2 + k1k2*(k1k3 - 2*p1k1 - p1k3) + k2k3*(-p1k1 + p1p2) - k1k3*p2k2 + 
     -            2*p1k1*p2k2 + p1k3*p2k2 - p1k2*p2k3) + 
     -         g2Lc*g2R*muu4*muuq*muuqc*(-k1k3 + 2*p1k1 + p1k3)*sqrtS**2 + 
     -         g2R*g2Rc*muuq*(p1k2*(k1k3 - 2*(p1k1 + p1k3)) + 
     -            (k2k3*muu1**2 + muu4**2*(-k1k3 + 2*p1k1 + p1k3))*sqrtS**2) + 
     -         g2L*g2Rc*muu4*(k1k2*(k1k3 - 2*p1k1 - p1k3) - 
     -            (k2k3 - 2*(p1k1 + p1k3))*(p1k1 - p1p2) + (2*p1k1 + p1k3)*p2k2 - p1k2*p2k3 + 
     -            muu1**2*p2k3*sqrtS**2 - k1k3*(p1k1 - p1k2 - p1p2 + p2k2 + muu1**2*sqrtS**2)))+
     -       g1R*g1Rc*(g2L*(-(g2Lc*(k1k3 - 2*(p1k1 + p1k3))*
     -               (-(k1k2*p1p2) + p1k2*p2k1 + p1k1*p2k2)) + 
     -            (g2Rc*muu4*muuqc*(-(p1k3*p2k1) + k1k3*(p1p2 + p2k1) - p1k1*(2*p2k1 + p2k3)) + 
     -               g2Lc*(-2*muu3**2*p1k2*(p1k1 + p1k3) + 
     -                  k1k3*(muu3**2*p1k2 + muu4**2*(p1p2 + p2k1) - muu1**2*p2k2) + 
     -                  muu1**2*(-(k2k3*p2k1) + k1k2*p2k3) - 
     -                  muu4**2*(p1k3*p2k1 + p1k1*(2*p2k1 + p2k3))))*sqrtS**2 + 
     -            muu3**2*(g2Lc*k2k3*muu1**2 + g2Lc*muu4**2*(-k1k3 + 2*p1k1 + p1k3) + 
     -               g2Rc*muu4*muuqc*(-k1k3 + 2*p1k1 + p1k3))*sqrtS**4) - 
     -         g2R*muuq*sqrtS**2*(g2Rc*muuqc*
     -             (k2k3*p1k1 - k1k3*p1k2 + k1k2*(-k1k3 + 2*p1k1 + p1k3)) + 
     -            g2Lc*muu4*(k1k2*(-k1k3 + 2*p1k1 + p1k3) + p1k1*(k2k3 - 2*(p1k1 + p1k3)) + 
     -               k1k3*(p1k1 - p1k2 + muu1**2*sqrtS**2)))) + 
     -      g1L*g1Lc*(g2R*(-(g2Rc*(k1k3 - 2*(p1k1 + p1k3))*
     -               (-(k1k2*p1p2) + p1k2*p2k1 + p1k1*p2k2)) + 
     -            (g2Lc*muu4*muuqc*(-(p1k3*p2k1) + k1k3*(p1p2 + p2k1) - p1k1*(2*p2k1 + p2k3)) + 
     -               g2Rc*(-2*muu3**2*p1k2*(p1k1 + p1k3) + 
     -                  k1k3*(muu3**2*p1k2 + muu4**2*(p1p2 + p2k1) - muu1**2*p2k2) + 
     -                  muu1**2*(-(k2k3*p2k1) + k1k2*p2k3) - 
     -                  muu4**2*(p1k3*p2k1 + p1k1*(2*p2k1 + p2k3))))*sqrtS**2 + 
     -            muu3**2*(g2Rc*k2k3*muu1**2 + g2Rc*muu4**2*(-k1k3 + 2*p1k1 + p1k3) + 
     -               g2Lc*muu4*muuqc*(-k1k3 + 2*p1k1 + p1k3))*sqrtS**4) - 
     -         g2L*muuq*sqrtS**2*(g2Lc*muuqc*
     -             (k2k3*p1k1 - k1k3*p1k2 + k1k2*(-k1k3 + 2*p1k1 + p1k3)) + 
     -            g2Rc*muu4*(k1k2*(-k1k3 + 2*p1k1 + p1k3) + p1k1*(k2k3 - 2*(p1k1 + p1k3)) + 
     -               k1k3*(p1k1 - p1k2 + muu1**2*sqrtS**2))))))/
     -  (Denlegk1c*Denlegp1*DenU1*DenU1c)
     
     
! ! ---------------------k1-p1-------------------------------

	  colfac = Nc/2d0*CF
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(4*gAc*gD*(g1L*g1Rc*muu3*sqrtS**2*
     -       (g2L*g2Lc*muuq*(k1k3*p1k2 + k1k2*(k1k3 - 2*p1k1 - p1k3) + k2k3*(-p1k1 + p1p2) - 
     -            k1k3*p2k2 + 2*p1k1*p2k2 + p1k3*p2k2 - p1k2*p2k3) + 
     -         g2L*g2Rc*muu4*muuq*muuqc*(-k1k3 + 2*p1k1 + p1k3)*sqrtS**2 + 
     -         g2R*g2Rc*muuqc*(p1k2*(k1k3 - 2*(p1k1 + p1k3)) + 
     -            (k2k3*muu1**2 + muu4**2*(-k1k3 + 2*p1k1 + p1k3))*sqrtS**2) + 
     -         g2Lc*g2R*muu4*(k1k2*(k1k3 - 2*p1k1 - p1k3) - 
     -            (k2k3 - 2*(p1k1 + p1k3))*(p1k1 - p1p2) + (2*p1k1 + p1k3)*p2k2 - p1k2*p2k3 + 
     -            muu1**2*p2k3*sqrtS**2 - k1k3*(p1k1 - p1k2 - p1p2 + p2k2 + muu1**2*sqrtS**2))) +
     -        g1Lc*g1R*muu3*sqrtS**2*(g2R*g2Rc*muuq*
     -          (k1k3*p1k2 + k1k2*(k1k3 - 2*p1k1 - p1k3) + k2k3*(-p1k1 + p1p2) - k1k3*p2k2 + 
     -            2*p1k1*p2k2 + p1k3*p2k2 - p1k2*p2k3) + 
     -         g2Lc*g2R*muu4*muuq*muuqc*(-k1k3 + 2*p1k1 + p1k3)*sqrtS**2 + 
     -         g2L*g2Lc*muuqc*(p1k2*(k1k3 - 2*(p1k1 + p1k3)) + 
     -            (k2k3*muu1**2 + muu4**2*(-k1k3 + 2*p1k1 + p1k3))*sqrtS**2) + 
     -         g2L*g2Rc*muu4*(k1k2*(k1k3 - 2*p1k1 - p1k3) - 
     -            (k2k3 - 2*(p1k1 + p1k3))*(p1k1 - p1p2) + (2*p1k1 + p1k3)*p2k2 - p1k2*p2k3 + 
     -            muu1**2*p2k3*sqrtS**2 - k1k3*(p1k1 - p1k2 - p1p2 + p2k2 + muu1**2*sqrtS**2)))+
     -        g1R*g1Rc*(g2L*(-(g2Lc*(k1k3 - 2*(p1k1 + p1k3))*
     -               (-(k1k2*p1p2) + p1k2*p2k1 + p1k1*p2k2)) - 
     -            (g2Rc*muu4*muuqc*(k1k3*(p1k1 - p1k2) + k1k2*(-k1k3 + 2*p1k1 + p1k3) + 
     -                  p1k1*(k2k3 - 2*(p1k1 + p1k3))) + 
     -               g2Lc*(2*muu3**2*p1k2*(p1k1 + p1k3) - 
     -                  k1k3*(muu3**2*p1k2 + muu4**2*(p1p2 + p2k1) - muu1**2*p2k2) + 
     -                  muu1**2*(k2k3*p2k1 - k1k2*p2k3) + 
     -                  muu4**2*(p1k3*p2k1 + p1k1*(2*p2k1 + p2k3))))*sqrtS**2 + 
     -            (-(g2Rc*k1k3*muu1**2*muu4*muuqc) + 
     -               g2Lc*muu3**2*(k2k3*muu1**2 + muu4**2*(-k1k3 + 2*p1k1 + p1k3)))*sqrtS**4) + 
     -         g2R*muuq*sqrtS**2*(g2Rc*muuqc*
     -             (-(k2k3*p1k1) + k1k3*p1k2 + k1k2*(k1k3 - 2*p1k1 - p1k3)) + 
     -            g2Lc*muu4*(k1k3*p1p2 + k1k3*p2k1 - 2*p1k1*p2k1 - p1k3*p2k1 - p1k1*p2k3 + 
     -               muu3**2*(-k1k3 + 2*p1k1 + p1k3)*sqrtS**2))) + 
     -      g1L*g1Lc*(g2R*(-(g2Rc*(k1k3 - 2*(p1k1 + p1k3))*
     -               (-(k1k2*p1p2) + p1k2*p2k1 + p1k1*p2k2)) - 
     -            (g2Lc*muu4*muuqc*(k1k3*(p1k1 - p1k2) + k1k2*(-k1k3 + 2*p1k1 + p1k3) + 
     -                  p1k1*(k2k3 - 2*(p1k1 + p1k3))) + 
     -               g2Rc*(2*muu3**2*p1k2*(p1k1 + p1k3) - 
     -                  k1k3*(muu3**2*p1k2 + muu4**2*(p1p2 + p2k1) - muu1**2*p2k2) + 
     -                  muu1**2*(k2k3*p2k1 - k1k2*p2k3) + 
     -                  muu4**2*(p1k3*p2k1 + p1k1*(2*p2k1 + p2k3))))*sqrtS**2 + 
     -            (-(g2Lc*k1k3*muu1**2*muu4*muuqc) + 
     -               g2Rc*muu3**2*(k2k3*muu1**2 + muu4**2*(-k1k3 + 2*p1k1 + p1k3)))*sqrtS**4) + 
     -         g2L*muuq*sqrtS**2*(g2Lc*muuqc*
     -             (-(k2k3*p1k1) + k1k3*p1k2 + k1k2*(k1k3 - 2*p1k1 - p1k3)) + 
     -            g2Rc*muu4*(k1k3*p1p2 + k1k3*p2k1 - 2*p1k1*p2k1 - p1k3*p2k1 - p1k1*p2k3 + 
     -               muu3**2*(-k1k3 + 2*p1k1 + p1k3)*sqrtS**2)))))/
     -  (Denlegk1*Denlegp1c*DenU1*DenU1c)
	  ! 	  
! ! ---------------------p1-k2------------------------------

	  colfac = 0d0
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(4*gA*gEc*(g1L*(2*g1Lc*g2R*g2Rc*p2k1*
     -          (-(p1k3*p2k2) + k2k3*(p1p2 + p2k2) - p1k2*(2*p2k2 + p2k3)) - 
     -         (g1Lc*(g2L*g2Lc*muuq*muuqc*
     -                (-(k2k3*p1k1) + k1k3*p1k2 + k1k2*(-k2k3 + 2*p1k2 + p1k3)) + 
     -               g2Lc*g2R*muu4*muuqc*
     -                (k1k3*p1p2 + (-k2k3 + 2*p1k2 + p1k3)*p2k1 - p1k1*p2k3) + 
     -               g2L*g2Rc*muu4*muuq*
     -                (-(k1k3*p1p2) + (-k2k3 + 2*p1k2 + p1k3)*p2k1 + p1k1*p2k3) + 
     -               g2R*g2Rc*(k2k3*muu2**2*p1k1 - k1k3*muu2**2*p1k2 + k1k3*muu3**2*p1k2 + 
     -                  k1k2*(muu2 - muu3)*(muu2 + muu3)*(k2k3 - 2*p1k2 - p1k3) - 
     -                  4*muu3**2*p1k2*p2k2 - 2*muu3**2*p1k3*p2k2 + 
     -                  k2k3*muu3**2*(-p1k1 + 2*(p1p2 + p2k2)) - 2*muu3**2*p1k2*p2k3)) + 
     -            g1Rc*muu3*(g2L*g2Lc*muuq*
     -                (k1k3*p1k2 + k1k2*(-k2k3 + 2*p1k2 + p1k3) - (2*p1k2 + p1k3)*p2k2 + 
     -                  k2k3*(-p1k1 + p1p2 + p2k2) - p1k2*p2k3) + 
     -               g2R*(2*g2Lc*muu4*(-k2k3 + 2*p1k2 + p1k3)*p2k1 + 
     -                  g2Rc*muuqc*(k1k3*p1k2 + k1k2*(-k2k3 + 2*p1k2 + p1k3) - 
     -                     (2*p1k2 + p1k3)*p2k2 + k2k3*(-p1k1 + p1p2 + p2k2) - p1k2*p2k3))))*
     -          sqrtS**2 - muu3*muu4*(g1Rc*g2Lc*g2R*(muu2**2 + muu3**2) + 
     -            g1Rc*g2L*g2Rc*muuq*muuqc + g1Lc*muu3*(g2L*g2Rc*muuq + g2Lc*g2R*muuqc))*
     -          (k2k3 - 2*p1k2 - p1k3)*sqrtS**4) + 
     -      g1R*(2*g1Rc*g2L*g2Lc*p2k1*(-(p1k3*p2k2) + k2k3*(p1p2 + p2k2) - 
     -            p1k2*(2*p2k2 + p2k3)) - 
     -         (g1Rc*(g2R*g2Rc*muuq*muuqc*
     -                (-(k2k3*p1k1) + k1k3*p1k2 + k1k2*(-k2k3 + 2*p1k2 + p1k3)) + 
     -               g2L*g2Rc*muu4*muuqc*
     -                (k1k3*p1p2 + (-k2k3 + 2*p1k2 + p1k3)*p2k1 - p1k1*p2k3) + 
     -               g2Lc*g2R*muu4*muuq*
     -                (-(k1k3*p1p2) + (-k2k3 + 2*p1k2 + p1k3)*p2k1 + p1k1*p2k3) + 
     -               g2L*g2Lc*(k2k3*muu2**2*p1k1 - k1k3*muu2**2*p1k2 + k1k3*muu3**2*p1k2 + 
     -                  k1k2*(muu2 - muu3)*(muu2 + muu3)*(k2k3 - 2*p1k2 - p1k3) - 
     -                  4*muu3**2*p1k2*p2k2 - 2*muu3**2*p1k3*p2k2 + 
     -                  k2k3*muu3**2*(-p1k1 + 2*(p1p2 + p2k2)) - 2*muu3**2*p1k2*p2k3)) + 
     -            g1Lc*muu3*(g2R*g2Rc*muuq*
     -                (k1k3*p1k2 + k1k2*(-k2k3 + 2*p1k2 + p1k3) - (2*p1k2 + p1k3)*p2k2 + 
     -                  k2k3*(-p1k1 + p1p2 + p2k2) - p1k2*p2k3) + 
     -               g2L*(2*g2Rc*muu4*(-k2k3 + 2*p1k2 + p1k3)*p2k1 + 
     -                  g2Lc*muuqc*(k1k3*p1k2 + k1k2*(-k2k3 + 2*p1k2 + p1k3) - 
     -                     (2*p1k2 + p1k3)*p2k2 + k2k3*(-p1k1 + p1p2 + p2k2) - p1k2*p2k3))))*
     -          sqrtS**2 - muu3*muu4*(g1Lc*g2L*g2Rc*(muu2**2 + muu3**2) + 
     -            g1Lc*g2Lc*g2R*muuq*muuqc + g1Rc*muu3*(g2Lc*g2R*muuq + g2L*g2Rc*muuqc))*
     -          (k2k3 - 2*p1k2 - p1k3)*sqrtS**4)))/(Denlegk2c*Denlegp1*DenU1*DenU2c)
     
! ! ---------------------k2-p1------------------------------

	  colfac = 0d0
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(4*gAc*gE*(g1L*(2*g1Lc*g2R*g2Rc*p2k1*
     -          (-(p1k3*p2k2) + k2k3*(p1p2 + p2k2) - p1k2*(2*p2k2 + p2k3)) - 
     -         (g1Lc*(g2L*g2Lc*muuq*muuqc*
     -                (-(k2k3*p1k1) + k1k3*p1k2 + k1k2*(-k2k3 + 2*p1k2 + p1k3)) + 
     -               g2L*g2Rc*muu4*muuq*
     -                (k1k3*p1p2 + (-k2k3 + 2*p1k2 + p1k3)*p2k1 - p1k1*p2k3) + 
     -               g2Lc*g2R*muu4*muuqc*
     -                (-(k1k3*p1p2) + (-k2k3 + 2*p1k2 + p1k3)*p2k1 + p1k1*p2k3) + 
     -               g2R*g2Rc*(k2k3*muu2**2*p1k1 - k1k3*muu2**2*p1k2 + k1k3*muu3**2*p1k2 + 
     -                  k1k2*(muu2 - muu3)*(muu2 + muu3)*(k2k3 - 2*p1k2 - p1k3) - 
     -                  4*muu3**2*p1k2*p2k2 - 2*muu3**2*p1k3*p2k2 + 
     -                  k2k3*muu3**2*(-p1k1 + 2*(p1p2 + p2k2)) - 2*muu3**2*p1k2*p2k3)) + 
     -            g1Rc*muu3*(g2L*g2Lc*muuq*
     -                (k1k3*p1k2 + k1k2*(-k2k3 + 2*p1k2 + p1k3) - (2*p1k2 + p1k3)*p2k2 + 
     -                  k2k3*(-p1k1 + p1p2 + p2k2) - p1k2*p2k3) + 
     -               g2R*(2*g2Lc*muu4*(-k2k3 + 2*p1k2 + p1k3)*p2k1 + 
     -                  g2Rc*muuqc*(k1k3*p1k2 + k1k2*(-k2k3 + 2*p1k2 + p1k3) - 
     -                     (2*p1k2 + p1k3)*p2k2 + k2k3*(-p1k1 + p1p2 + p2k2) - p1k2*p2k3))))*
     -          sqrtS**2 - muu3*muu4*(g1Rc*g2Lc*g2R*(muu2**2 + muu3**2) + 
     -            g1Rc*g2L*g2Rc*muuq*muuqc + g1Lc*muu3*(g2L*g2Rc*muuq + g2Lc*g2R*muuqc))*
     -          (k2k3 - 2*p1k2 - p1k3)*sqrtS**4) + 
     -      g1R*(2*g1Rc*g2L*g2Lc*p2k1*(-(p1k3*p2k2) + k2k3*(p1p2 + p2k2) - 
     -            p1k2*(2*p2k2 + p2k3)) - 
     -         (g1Rc*(g2R*g2Rc*muuq*muuqc*
     -                (-(k2k3*p1k1) + k1k3*p1k2 + k1k2*(-k2k3 + 2*p1k2 + p1k3)) + 
     -               g2Lc*g2R*muu4*muuq*
     -                (k1k3*p1p2 + (-k2k3 + 2*p1k2 + p1k3)*p2k1 - p1k1*p2k3) + 
     -               g2L*g2Rc*muu4*muuqc*
     -                (-(k1k3*p1p2) + (-k2k3 + 2*p1k2 + p1k3)*p2k1 + p1k1*p2k3) + 
     -               g2L*g2Lc*(k2k3*muu2**2*p1k1 - k1k3*muu2**2*p1k2 + k1k3*muu3**2*p1k2 + 
     -                  k1k2*(muu2 - muu3)*(muu2 + muu3)*(k2k3 - 2*p1k2 - p1k3) - 
     -                  4*muu3**2*p1k2*p2k2 - 2*muu3**2*p1k3*p2k2 + 
     -                  k2k3*muu3**2*(-p1k1 + 2*(p1p2 + p2k2)) - 2*muu3**2*p1k2*p2k3)) + 
     -            g1Lc*muu3*(g2R*g2Rc*muuq*
     -                (k1k3*p1k2 + k1k2*(-k2k3 + 2*p1k2 + p1k3) - (2*p1k2 + p1k3)*p2k2 + 
     -                  k2k3*(-p1k1 + p1p2 + p2k2) - p1k2*p2k3) + 
     -               g2L*(2*g2Rc*muu4*(-k2k3 + 2*p1k2 + p1k3)*p2k1 + 
     -                  g2Lc*muuqc*(k1k3*p1k2 + k1k2*(-k2k3 + 2*p1k2 + p1k3) - 
     -                     (2*p1k2 + p1k3)*p2k2 + k2k3*(-p1k1 + p1p2 + p2k2) - p1k2*p2k3))))*
     -          sqrtS**2 - muu3*muu4*(g1Lc*g2L*g2Rc*(muu2**2 + muu3**2) + 
     -            g1Lc*g2Lc*g2R*muuq*muuqc + g1Rc*muu3*(g2Lc*g2R*muuq + g2L*g2Rc*muuqc))*
     -          (k2k3 - 2*p1k2 - p1k3)*sqrtS**4)))/(Denlegk2*Denlegp1c*DenU1c*DenU2)
! 	  
! ! ---------------------p2-p2------------------------------
! 
	  colfac = 0d0
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(-8*gB*gBc*(p2k3 - muu2**2*sqrtS**2)*
     -    (-2*(g1R*g1Rc*g2L*g2Lc + g1L*g1Lc*g2R*g2Rc)*p1k1*p1k2 + 
     -      (g1L*g1Lc*(g2Rc*k1k2*(g2R*(muu1 - muu4)*(muu1 + muu4) - g2L*muu4*muuq) - 
     -            g2Lc*k1k2*(g2R*muu4 + g2L*muuq)*muuqc + 
     -            muu4*(2*g2R*g2Rc*muu4 + g2L*g2Rc*muuq + g2Lc*g2R*muuqc)*p1k1) + 
     -         g1R*g1Rc*(g2Lc*k1k2*(g2L*(muu1 - muu4)*(muu1 + muu4) - g2R*muu4*muuq) - 
     -            g2Rc*k1k2*(g2L*muu4 + g2R*muuq)*muuqc + 
     -            muu4*(2*g2L*g2Lc*muu4 + g2Lc*g2R*muuq + g2L*g2Rc*muuqc)*p1k1) - 
     -         g1Lc*g1R*muu3*(2*g2L*g2Rc*muu4 + g2R*g2Rc*muuq + g2L*g2Lc*muuqc)*p1k2 - 
     -         g1L*g1Rc*muu3*(2*g2Lc*g2R*muu4 + g2L*g2Lc*muuq + g2R*g2Rc*muuqc)*p1k2)*sqrtS**2 +
     -       muu3*muu4*(g1L*g1Rc*(g2Lc*g2R*(muu1**2 + muu4**2) + g2L*g2Lc*muu4*muuq + 
     -            g2Rc*(g2R*muu4 + g2L*muuq)*muuqc) + 
     -         g1Lc*g1R*(g2L*g2Rc*(muu1**2 + muu4**2) + g2L*g2Lc*muu4*muuqc + 
     -            g2R*muuq*(g2Rc*muu4 + g2Lc*muuqc)))*sqrtS**4))/
     -  (Denlegp2*Denlegp2c*DenU2*DenU2c)

! ! ---------------------p2-k1------------------------------
! 
	  colfac = 0d0

	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(4*gB*gDc*(g1R*g1Rc*(g2R*muuq*(g2Rc*muuqc*
     -             (-(k2k3*p2k1) + k1k3*p2k2 + k1k2*(k1k3 - 2*p2k1 - p2k3)) + 
     -            g2Lc*muu4*(-(k2k3*p2k1) + 2*p1k1*p2k1 + p1k3*p2k1 - 
     -               k1k3*(p1k1 + p1p2 - p2k2) + k1k2*(k1k3 - 2*p2k1 - p2k3) + p1k1*p2k3))*
     -          sqrtS**2 + g2L*(2*g2Lc*p1k2*
     -             (k1k3*(p1k1 + p1p2) - p1k3*p2k1 - p1k1*(2*p2k1 + p2k3)) + 
     -            g2Rc*muu4*muuqc*(-(k2k3*p2k1) + 2*p1k1*p2k1 + p1k3*p2k1 - 
     -               k1k3*(p1k1 + p1p2 - p2k2) + k1k2*(k1k3 - 2*p2k1 - p2k3) + p1k1*p2k3)*
     -             sqrtS**2 + g2Lc*(-2*k1k3*muu4**2*(p1k1 + p1p2) + k2k3*muu1**2*p2k1 + 
     -               k1k3*(-muu1 + muu4)*(muu1 + muu4)*p2k2 - 
     -               k1k2*(muu1 - muu4)*(muu1 + muu4)*(k1k3 - 2*p2k1 - p2k3) + 
     -               muu4**2*(-(k2k3*p2k1) + 4*p1k1*p2k1 + 2*p1k3*p2k1 + 2*p1k1*p2k3))*sqrtS**2)
     -         ) + g1L*g1Lc*(g2L*muuq*(g2Lc*muuqc*
     -             (-(k2k3*p2k1) + k1k3*p2k2 + k1k2*(k1k3 - 2*p2k1 - p2k3)) + 
     -            g2Rc*muu4*(-(k2k3*p2k1) + 2*p1k1*p2k1 + p1k3*p2k1 - 
     -               k1k3*(p1k1 + p1p2 - p2k2) + k1k2*(k1k3 - 2*p2k1 - p2k3) + p1k1*p2k3))*
     -          sqrtS**2 + g2R*(2*g2Rc*p1k2*
     -             (k1k3*(p1k1 + p1p2) - p1k3*p2k1 - p1k1*(2*p2k1 + p2k3)) + 
     -            g2Lc*muu4*muuqc*(-(k2k3*p2k1) + 2*p1k1*p2k1 + p1k3*p2k1 - 
     -               k1k3*(p1k1 + p1p2 - p2k2) + k1k2*(k1k3 - 2*p2k1 - p2k3) + p1k1*p2k3)*
     -             sqrtS**2 + g2Rc*(-2*k1k3*muu4**2*(p1k1 + p1p2) + k2k3*muu1**2*p2k1 + 
     -               k1k3*(-muu1 + muu4)*(muu1 + muu4)*p2k2 - 
     -               k1k2*(muu1 - muu4)*(muu1 + muu4)*(k1k3 - 2*p2k1 - p2k3) + 
     -               muu4**2*(-(k2k3*p2k1) + 4*p1k1*p2k1 + 2*p1k3*p2k1 + 2*p1k1*p2k3))*sqrtS**2)
     -         ) + g1L*g1Rc*muu3*sqrtS**2*
     -       (g2R*g2Rc*muuqc*(k1k3*p1k2 - k2k3*p1p2 + p1k3*p2k2 - p1k2*(2*p2k1 + p2k3)) - 
     -         g2Rc*muu4*(g2R*muu4 + g2L*muuq)*muuqc*(k1k3 - 2*p2k1 - p2k3)*sqrtS**2 - 
     -         g2Lc*g2R*muu4*(k1k3 - 2*p2k1 - p2k3)*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2) + 
     -         g2L*g2Lc*muuq*(k2k3*p1p2 - p1k3*p2k2 + 
     -            (k1k3 - 2*p2k1 - p2k3)*(p1k2 - muu4**2*sqrtS**2))) + 
     -      g1Lc*g1R*muu3*sqrtS**2*(g2R*g2Rc*muuq*
     -          (k1k3*p1k2 + k2k3*p1p2 - p1k3*p2k2 - p1k2*(2*p2k1 + p2k3)) - 
     -         g2R*muu4*muuq*(g2Rc*muu4 + g2Lc*muuqc)*(k1k3 - 2*p2k1 - p2k3)*sqrtS**2 - 
     -         g2L*g2Rc*muu4*(k1k3 - 2*p2k1 - p2k3)*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2) + 
     -         g2L*g2Lc*muuqc*(-(k2k3*p1p2) + p1k3*p2k2 + 
     -            (k1k3 - 2*p2k1 - p2k3)*(p1k2 - muu4**2*sqrtS**2)))))/
     -  (Denlegk1c*Denlegp2*DenU1c*DenU2)
     
     
! ! ---------------------k1-p2-------------------------------
! 
	  colfac = 0d0

	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(4*gBc*gD*(g1R*g1Rc*(g2R*muuq*(g2Rc*muuqc*
     -             (-(k2k3*p2k1) + k1k3*p2k2 + k1k2*(k1k3 - 2*p2k1 - p2k3)) + 
     -            g2Lc*muu4*(-(k2k3*p2k1) + 2*p1k1*p2k1 + p1k3*p2k1 - 
     -               k1k3*(p1k1 + p1p2 - p2k2) + k1k2*(k1k3 - 2*p2k1 - p2k3) + p1k1*p2k3))*
     -          sqrtS**2 + g2L*(2*g2Lc*p1k2*
     -             (k1k3*(p1k1 + p1p2) - p1k3*p2k1 - p1k1*(2*p2k1 + p2k3)) + 
     -            g2Rc*muu4*muuqc*(-(k2k3*p2k1) + 2*p1k1*p2k1 + p1k3*p2k1 - 
     -               k1k3*(p1k1 + p1p2 - p2k2) + k1k2*(k1k3 - 2*p2k1 - p2k3) + p1k1*p2k3)*
     -             sqrtS**2 + g2Lc*(-2*k1k3*muu4**2*(p1k1 + p1p2) + k2k3*muu1**2*p2k1 + 
     -               k1k3*(-muu1 + muu4)*(muu1 + muu4)*p2k2 - 
     -               k1k2*(muu1 - muu4)*(muu1 + muu4)*(k1k3 - 2*p2k1 - p2k3) + 
     -               muu4**2*(-(k2k3*p2k1) + 4*p1k1*p2k1 + 2*p1k3*p2k1 + 2*p1k1*p2k3))*sqrtS**2)
     -         ) + g1L*g1Lc*(g2L*muuq*(g2Lc*muuqc*
     -             (-(k2k3*p2k1) + k1k3*p2k2 + k1k2*(k1k3 - 2*p2k1 - p2k3)) + 
     -            g2Rc*muu4*(-(k2k3*p2k1) + 2*p1k1*p2k1 + p1k3*p2k1 - 
     -               k1k3*(p1k1 + p1p2 - p2k2) + k1k2*(k1k3 - 2*p2k1 - p2k3) + p1k1*p2k3))*
     -          sqrtS**2 + g2R*(2*g2Rc*p1k2*
     -             (k1k3*(p1k1 + p1p2) - p1k3*p2k1 - p1k1*(2*p2k1 + p2k3)) + 
     -            g2Lc*muu4*muuqc*(-(k2k3*p2k1) + 2*p1k1*p2k1 + p1k3*p2k1 - 
     -               k1k3*(p1k1 + p1p2 - p2k2) + k1k2*(k1k3 - 2*p2k1 - p2k3) + p1k1*p2k3)*
     -             sqrtS**2 + g2Rc*(-2*k1k3*muu4**2*(p1k1 + p1p2) + k2k3*muu1**2*p2k1 + 
     -               k1k3*(-muu1 + muu4)*(muu1 + muu4)*p2k2 - 
     -               k1k2*(muu1 - muu4)*(muu1 + muu4)*(k1k3 - 2*p2k1 - p2k3) + 
     -               muu4**2*(-(k2k3*p2k1) + 4*p1k1*p2k1 + 2*p1k3*p2k1 + 2*p1k1*p2k3))*sqrtS**2)
     -         ) + g1Lc*g1R*muu3*sqrtS**2*
     -       (g2R*g2Rc*muuq*(k1k3*p1k2 - k2k3*p1p2 + p1k3*p2k2 - p1k2*(2*p2k1 + p2k3)) - 
     -         g2R*muu4*muuq*(g2Rc*muu4 + g2Lc*muuqc)*(k1k3 - 2*p2k1 - p2k3)*sqrtS**2 - 
     -         g2L*g2Rc*muu4*(k1k3 - 2*p2k1 - p2k3)*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2) + 
     -         g2L*g2Lc*muuqc*(k2k3*p1p2 - p1k3*p2k2 + 
     -            (k1k3 - 2*p2k1 - p2k3)*(p1k2 - muu4**2*sqrtS**2))) + 
     -      g1L*g1Rc*muu3*sqrtS**2*(g2R*g2Rc*muuqc*
     -          (k1k3*p1k2 + k2k3*p1p2 - p1k3*p2k2 - p1k2*(2*p2k1 + p2k3)) - 
     -         g2Rc*muu4*(g2R*muu4 + g2L*muuq)*muuqc*(k1k3 - 2*p2k1 - p2k3)*sqrtS**2 - 
     -         g2Lc*g2R*muu4*(k1k3 - 2*p2k1 - p2k3)*(-2*p1k2 + (muu1**2 + muu4**2)*sqrtS**2) + 
     -         g2L*g2Lc*muuq*(-(k2k3*p1p2) + p1k3*p2k2 + 
     -            (k1k3 - 2*p2k1 - p2k3)*(p1k2 - muu4**2*sqrtS**2)))))/
     -  (Denlegk1*Denlegp2c*DenU1*DenU2c)
	  
	  
	  
! ! ---------------------p2-k2------------------------------
! 
	  colfac = Nc/2d0*CF
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(4*gB*gEc*(g1L*g1Lc*(g2R*g2Rc*(k1k2*p1p2 - p1k2*p2k1 - p1k1*p2k2)*
     -          (k2k3 - 2*(p2k2 + p2k3)) + 
     -         (g2Lc*muuqc*(g2L*muuq*(k2k3*p2k1 - k1k3*p2k2 + k1k2*(k2k3 - 2*p2k2 - p2k3)) + 
     -               g2R*muu4*(-(k2k3*p1k1) + k1k3*p1p2 + k2k3*p2k1 - p1k3*p2k1 - k1k3*p2k2 + 
     -                  2*p1k1*p2k2 + k1k2*(k2k3 - 2*p2k2 - p2k3) + p1k1*p2k3)) - 
     -            g2Rc*(g2R*(k1k3*muu2**2*p1k2 - k1k2*muu2**2*p1k3 + 
     -                  k2k3*(muu2**2*p1k1 - muu3**2*(p1k2 + p1p2) - muu4**2*p2k1) + 
     -                  2*muu3**2*p1k2*p2k2 + muu3**2*p1k3*p2k2 + 2*muu4**2*p2k1*p2k2 + 
     -                  muu3**2*p1k2*p2k3 + 2*muu4**2*p2k1*p2k3) + 
     -               g2L*muu4*muuq*p2k1*(-k2k3 + 2*(p2k2 + p2k3))))*sqrtS**2 + 
     -         g2Rc*muu4*(g2R*muu4 + g2L*muuq)*(k1k3*muu2**2 + muu3**2*(-k2k3 + 2*p2k2 + p2k3))*
     -          sqrtS**4) + g1R*g1Rc*(g2L*g2Lc*(k1k2*p1p2 - p1k2*p2k1 - p1k1*p2k2)*
     -          (k2k3 - 2*(p2k2 + p2k3)) + 
     -         (g2R*g2Rc*muuq*muuqc*(k2k3*p2k1 - k1k3*p2k2 + k1k2*(k2k3 - 2*p2k2 - p2k3)) + 
     -            g2L*g2Rc*muu4*muuqc*(k1k3*p1p2 - p1k3*p2k1 + k2k3*(-p1k1 + p2k1) - 
     -               k1k3*p2k2 + 2*p1k1*p2k2 + k1k2*(k2k3 - 2*p2k2 - p2k3) + p1k1*p2k3) - 
     -            g2L*g2Lc*(k1k3*muu2**2*p1k2 - k1k2*muu2**2*p1k3 + 
     -               k2k3*(muu2**2*p1k1 - muu3**2*(p1k2 + p1p2) - muu4**2*p2k1) + 
     -               2*muu3**2*p1k2*p2k2 + muu3**2*p1k3*p2k2 + 2*muu4**2*p2k1*p2k2 + 
     -               muu3**2*p1k2*p2k3 + 2*muu4**2*p2k1*p2k3) + 
     -            g2Lc*g2R*muu4*muuq*p2k1*(k2k3 - 2*(p2k2 + p2k3)))*sqrtS**2 + 
     -         g2Lc*muu4*(g2L*muu4 + g2R*muuq)*(k1k3*muu2**2 + muu3**2*(-k2k3 + 2*p2k2 + p2k3))*
     -          sqrtS**4) - g1L*g1Rc*muu3*sqrtS**2*
     -       (g2R*g2Rc*muuqc*(-(k2k3*(p1k2 + p1p2)) + p1k3*p2k2 + p1k2*(2*p2k2 + p2k3)) + 
     -         g2Rc*muu4*(g2R*muu4 + g2L*muuq)*muuqc*(k2k3 - 2*p2k2 - p2k3)*sqrtS**2 + 
     -         g2L*g2Lc*muuq*(k1k2*(-k2k3 + 2*p2k2 + p2k3) + p2k2*(k1k3 - 2*(p2k2 + p2k3)) + 
     -            k2k3*(-p2k1 + p2k2 + muu2**2*sqrtS**2)) + 
     -         g2Lc*g2R*muu4*(k1k3*(-p1p2 + p2k2) - 2*p2k2*(p1k1 - p1p2 + p2k2) - 
     -            (p1k1 - 2*p1p2 + 2*p2k2)*p2k3 + k1k2*(-k2k3 + 2*p2k2 + p2k3) + 
     -            p1k3*(p2k1 - muu2**2*sqrtS**2) + 
     -            k2k3*(p1k1 - p1p2 - p2k1 + p2k2 + muu2**2*sqrtS**2))) - 
     -      g1Lc*g1R*muu3*sqrtS**2*(g2Lc*g2R*muu4*muuq*muuqc*(k2k3 - 2*p2k2 - p2k3)*sqrtS**2 + 
     -         g2L*g2Lc*muuqc*(-(k2k3*p1k2) - k2k3*p1p2 + 2*p1k2*p2k2 + p1k3*p2k2 + p1k2*p2k3 + 
     -            muu4**2*(k2k3 - 2*p2k2 - p2k3)*sqrtS**2) + 
     -         g2R*g2Rc*muuq*(k1k2*(-k2k3 + 2*p2k2 + p2k3) + p2k2*(k1k3 - 2*(p2k2 + p2k3)) + 
     -            k2k3*(-p2k1 + p2k2 + muu2**2*sqrtS**2)) + 
     -         g2L*g2Rc*muu4*(k1k3*(-p1p2 + p2k2) - 2*p2k2*(p1k1 - p1p2 + p2k2) - 
     -            (p1k1 - 2*p1p2 + 2*p2k2)*p2k3 + k1k2*(-k2k3 + 2*p2k2 + p2k3) + 
     -            p1k3*(p2k1 - muu2**2*sqrtS**2) + 
     -            k2k3*(p1k1 - p1p2 - p2k1 + p2k2 + muu2**2*sqrtS**2)))))/
     -  (Denlegk2c*Denlegp2*DenU2*DenU2c)
     
! ! ---------------------k2-p2-------------------------------
! 
	  colfac = Nc/2d0*CF
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(4*gBc*gE*(g1L*g1Lc*(-(g2R*g2Rc*(-(k1k2*p1p2) + p1k2*p2k1 + p1k1*p2k2)*
     -            (k2k3 - 2*(p2k2 + p2k3))) + 
     -         (g2L*muuq*(g2Lc*muuqc*(k2k3*p2k1 - k1k3*p2k2 + k1k2*(k2k3 - 2*p2k2 - p2k3)) + 
     -               g2Rc*muu4*(-(k2k3*p1k1) + k1k3*p1p2 + k2k3*p2k1 - p1k3*p2k1 - k1k3*p2k2 + 
     -                  2*p1k1*p2k2 + k1k2*(k2k3 - 2*p2k2 - p2k3) + p1k1*p2k3)) - 
     -            g2R*(g2Rc*(k1k3*muu2**2*p1k2 - k1k2*muu2**2*p1k3 + 
     -                  k2k3*(muu2**2*p1k1 - muu3**2*(p1k2 + p1p2) - muu4**2*p2k1) + 
     -                  2*muu3**2*p1k2*p2k2 + muu3**2*p1k3*p2k2 + 2*muu4**2*p2k1*p2k2 + 
     -                  muu3**2*p1k2*p2k3 + 2*muu4**2*p2k1*p2k3) + 
     -               g2Lc*muu4*muuqc*p2k1*(-k2k3 + 2*(p2k2 + p2k3))))*sqrtS**2 + 
     -         g2R*muu4*(g2Rc*muu4 + g2Lc*muuqc)*
     -          (k1k3*muu2**2 + muu3**2*(-k2k3 + 2*p2k2 + p2k3))*sqrtS**4) + 
     -      g1R*g1Rc*(-(g2L*g2Lc*(-(k1k2*p1p2) + p1k2*p2k1 + p1k1*p2k2)*
     -            (k2k3 - 2*(p2k2 + p2k3))) + 
     -         (g2R*g2Rc*muuq*muuqc*(k2k3*p2k1 - k1k3*p2k2 + k1k2*(k2k3 - 2*p2k2 - p2k3)) + 
     -            g2Lc*g2R*muu4*muuq*(k1k3*p1p2 - p1k3*p2k1 + k2k3*(-p1k1 + p2k1) - k1k3*p2k2 + 
     -               2*p1k1*p2k2 + k1k2*(k2k3 - 2*p2k2 - p2k3) + p1k1*p2k3) - 
     -            g2L*g2Lc*(k1k3*muu2**2*p1k2 - k1k2*muu2**2*p1k3 + 
     -               k2k3*(muu2**2*p1k1 - muu3**2*(p1k2 + p1p2) - muu4**2*p2k1) + 
     -               2*muu3**2*p1k2*p2k2 + muu3**2*p1k3*p2k2 + 2*muu4**2*p2k1*p2k2 + 
     -               muu3**2*p1k2*p2k3 + 2*muu4**2*p2k1*p2k3) + 
     -            g2L*g2Rc*muu4*muuqc*p2k1*(k2k3 - 2*(p2k2 + p2k3)))*sqrtS**2 + 
     -         g2L*muu4*(g2Lc*muu4 + g2Rc*muuqc)*
     -          (k1k3*muu2**2 + muu3**2*(-k2k3 + 2*p2k2 + p2k3))*sqrtS**4) - 
     -      g1L*g1Rc*muu3*sqrtS**2*(g2L*g2Rc*muu4*muuq*muuqc*(k2k3 - 2*p2k2 - p2k3)*sqrtS**2 + 
     -         g2L*g2Lc*muuq*(-(k2k3*p1k2) - k2k3*p1p2 + 2*p1k2*p2k2 + p1k3*p2k2 + p1k2*p2k3 + 
     -            muu4**2*(k2k3 - 2*p2k2 - p2k3)*sqrtS**2) + 
     -         g2R*g2Rc*muuqc*(k1k2*(-k2k3 + 2*p2k2 + p2k3) + p2k2*(k1k3 - 2*(p2k2 + p2k3)) + 
     -            k2k3*(-p2k1 + p2k2 + muu2**2*sqrtS**2)) + 
     -         g2Lc*g2R*muu4*(k1k3*(-p1p2 + p2k2) - 2*p2k2*(p1k1 - p1p2 + p2k2) - 
     -            (p1k1 - 2*p1p2 + 2*p2k2)*p2k3 + k1k2*(-k2k3 + 2*p2k2 + p2k3) + 
     -            p1k3*(p2k1 - muu2**2*sqrtS**2) + 
     -            k2k3*(p1k1 - p1p2 - p2k1 + p2k2 + muu2**2*sqrtS**2))) - 
     -      g1Lc*g1R*muu3*sqrtS**2*(g2R*g2Rc*muuq*
     -          (-(k2k3*(p1k2 + p1p2)) + p1k3*p2k2 + p1k2*(2*p2k2 + p2k3)) + 
     -         g2R*muu4*muuq*(g2Rc*muu4 + g2Lc*muuqc)*(k2k3 - 2*p2k2 - p2k3)*sqrtS**2 + 
     -         g2L*g2Lc*muuqc*(k1k2*(-k2k3 + 2*p2k2 + p2k3) + p2k2*(k1k3 - 2*(p2k2 + p2k3)) + 
     -            k2k3*(-p2k1 + p2k2 + muu2**2*sqrtS**2)) + 
     -         g2L*g2Rc*muu4*(k1k3*(-p1p2 + p2k2) - 2*p2k2*(p1k1 - p1p2 + p2k2) - 
     -            (p1k1 - 2*p1p2 + 2*p2k2)*p2k3 + k1k2*(-k2k3 + 2*p2k2 + p2k3) + 
     -            p1k3*(p2k1 - muu2**2*sqrtS**2) + 
     -            k2k3*(p1k1 - p1p2 - p2k1 + p2k2 + muu2**2*sqrtS**2)))))/
     -  (Denlegk2*Denlegp2c*DenU2*DenU2c)     

! ! ---------------------k1-k1------------------------------
! 
	  colfac = 0d0
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(8*gD*gDc*(g1L*g1Rc*muu3*sqrtS**2*(k1k3 + muu3**2*sqrtS**2)*
     -       (-((2*g2Lc*g2R*muu4 + g2L*g2Lc*muuq + g2R*g2Rc*muuqc)*p1k2) + 
     -         muu4*(g2Lc*g2R*(muu1**2 + muu4**2) + g2L*g2Lc*muu4*muuq + 
     -            g2Rc*(g2R*muu4 + g2L*muuq)*muuqc)*sqrtS**2) + 
     -      g1Lc*g1R*muu3*sqrtS**2*(k1k3 + muu3**2*sqrtS**2)*
     -       (-((2*g2L*g2Rc*muu4 + g2R*g2Rc*muuq + g2L*g2Lc*muuqc)*p1k2) + 
     -         muu4*(g2L*g2Rc*(muu1**2 + muu4**2) + g2L*g2Lc*muu4*muuqc + 
     -            g2R*muuq*(g2Rc*muu4 + g2Lc*muuqc))*sqrtS**2) + 
     -      g1R*g1Rc*(-(g2R*muuq*sqrtS**2*
     -            (g2Rc*muuqc*(-(k1k3*k2k3) + (k1k2 + k2k3)*muu3**2*sqrtS**2) + 
     -              g2Lc*muu4*(k1k3*(-k2k3 + p1k3) + 
     -                 muu3**2*(k1k2 + k2k3 - p1k1 - p1k3)*sqrtS**2))) + 
     -         g2L*(g2Rc*muu4*muuqc*sqrtS**2*
     -             (k1k3*(k2k3 - p1k3) + muu3**2*(-k1k2 - k2k3 + p1k1 + p1k3)*sqrtS**2) + 
     -            g2Lc*(2*k1k3*p1k2*p1k3 + 
     -               k1k3*(-(k2k3*muu1**2) + k2k3*muu4**2 - 2*muu4**2*p1k3)*sqrtS**2 + 
     -               muu3**2*sqrtS**2*((k1k2 + k2k3)*(muu1 - muu4)*(muu1 + muu4)*sqrtS**2 - 
     -                  2*(p1k1 + p1k3)*(p1k2 - muu4**2*sqrtS**2))))) + 
     -      g1L*g1Lc*(-(g2L*muuq*sqrtS**2*
     -            (g2Lc*muuqc*(-(k1k3*k2k3) + (k1k2 + k2k3)*muu3**2*sqrtS**2) + 
     -              g2Rc*muu4*(k1k3*(-k2k3 + p1k3) + 
     -                 muu3**2*(k1k2 + k2k3 - p1k1 - p1k3)*sqrtS**2))) + 
     -         g2R*(g2Lc*muu4*muuqc*sqrtS**2*
     -             (k1k3*(k2k3 - p1k3) + muu3**2*(-k1k2 - k2k3 + p1k1 + p1k3)*sqrtS**2) + 
     -            g2Rc*(2*k1k3*p1k2*p1k3 + 
     -               k1k3*(-(k2k3*muu1**2) + k2k3*muu4**2 - 2*muu4**2*p1k3)*sqrtS**2 + 
     -               muu3**2*sqrtS**2*((k1k2 + k2k3)*(muu1 - muu4)*(muu1 + muu4)*sqrtS**2 - 
     -                  2*(p1k1 + p1k3)*(p1k2 - muu4**2*sqrtS**2)))))))/
     -  (Denlegk1*Denlegk1c*DenU1*DenU1c)
! ! ---------------------k1-k2------------------------------
! 
	  colfac = Nc/2d0*CF
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(-4*gD*gEc*(-(g1L*g1Rc*muu3*sqrtS**2*
     -         (g2R*g2Rc*muuqc*(k2k3*p1k1 - k1k3*p1k2 - k1k2*(2*p1k2 + p1k3)) + 
     -           g2Rc*muu4*muuqc*((2*k1k2 + k1k3 + k2k3)*(g2R*muu4 + g2L*muuq) - g2R*muu4*p1k3)*
     -            sqrtS**2 + g2L*g2Lc*muuq*
     -            (-2*k1k2**2 - 2*k1k3*k2k3 - k2k3*p2k1 + k1k3*p2k2 + 2*k2k3*p2k2 + 
     -              2*k2k3*p2k3 + k1k2*(-2*(k1k3 + k2k3 - p2k2) + p2k3) + 
     -              (k2k3*muu3**2 + muu4**2*(k1k3 - p2k3))*sqrtS**2) + 
     -           g2Lc*g2R*muu4*(-2*k1k2**2 + k2k3*p1k1 - k2k3*p1p2 - k2k3*p2k1 + p1k3*p2k1 + 
     -              2*k2k3*p2k2 - p1k3*p2k2 + 2*k2k3*p2k3 - p1k1*p2k3 + p1k2*p2k3 - 
     -              2*p1k3*p2k3 + k1k2*
     -               (-2*k1k3 - 2*k2k3 + 2*p1k1 + p1k3 - 2*p1p2 + 2*p2k2 + p2k3) + 
     -              (muu3**2*(k2k3 - p1k3) - muu4**2*p2k3)*sqrtS**2 + 
     -              k1k3*(-2*k2k3 + 2*p1k1 - p1k2 + 2*p1k3 - p1p2 + p2k2 + muu4**2*sqrtS**2))))-
     -       g1Lc*g1R*muu3*sqrtS**2*(g2R*g2Rc*muuq*
     -          (-2*k1k2**2 - k2k3*p2k1 + k1k3*(-2*k2k3 + p2k2) + 
     -            k1k2*(-2*(k1k3 + k2k3 - p2k2) + p2k3) + 2*k2k3*(p2k2 + p2k3)) + 
     -         g2R*muuq*(g2Rc*k2k3*muu3**2 + g2Lc*(2*k1k2 + k1k3 + k2k3)*muu4*muuqc + 
     -            g2Rc*muu4**2*(k1k3 - p2k3))*sqrtS**2 + 
     -         g2L*g2Lc*muuqc*(k2k3*p1k1 - 2*k1k2*p1k2 - k1k3*p1k2 - k1k2*p1k3 + 
     -            muu4**2*(2*k1k2 + k1k3 + k2k3 - p1k3)*sqrtS**2) + 
     -         g2L*g2Rc*muu4*(-2*k1k2**2 + k2k3*p1k1 - k2k3*p1p2 - k2k3*p2k1 + p1k3*p2k1 + 
     -            2*k2k3*p2k2 - p1k3*p2k2 + 2*k2k3*p2k3 - p1k1*p2k3 + p1k2*p2k3 - 2*p1k3*p2k3 + 
     -            k1k2*(-2*k1k3 - 2*k2k3 + 2*p1k1 + p1k3 - 2*p1p2 + 2*p2k2 + p2k3) + 
     -            (muu3**2*(k2k3 - p1k3) - muu4**2*p2k3)*sqrtS**2 + 
     -            k1k3*(-2*k2k3 + 2*p1k1 - p1k2 + 2*p1k3 - p1p2 + p2k2 + muu4**2*sqrtS**2))) + 
     -      g1R*g1Rc*(g2L*(2*g2Lc*(k1k2**2*p1p2 + k2k3*p1k1*(p2k1 - p2k2) + 
     -               k1k3*(-p1k1 + p1k2)*p2k2 + 
     -               k1k2*((k1k3 + k2k3)*p1p2 - p1k3*p2k1 - p1k1*p2k2 - p1k2*(p2k1 + p2k3))) + 
     -            (g2Rc*muu4*muuqc*(2*k1k2**2 + 2*k1k2*(k1k3 + k2k3 - p1k1) - k2k3*p1k1 + 
     -                  k1k3*(2*k2k3 - 2*p1k1 + p1k2 - 2*p1k3) - k1k2*p1k3) + 
     -               g2Lc*(k1k3*muu3**2*p1k2 + k1k2*muu3**2*(2*p1k2 + p1k3) - 
     -                  k1k3*muu4**2*p1p2 - k2k3*muu3**2*(p1k1 + p1p2) + k2k3*muu4**2*p2k1 - 
     -                  muu4**2*p1k3*p2k1 - k1k3*muu4**2*p2k2 + muu3**2*p1k3*p2k2 + 
     -                  muu4**2*p1k1*p2k3 - muu3**2*p1k2*p2k3 + k1k2*muu4**2*(2*p2k1 + p2k3)))*
     -             sqrtS**2 - muu4*(g2Rc*muuqc*(k1k3*muu4**2 + muu3**2*(k2k3 - p1k3)) + 
     -               g2Lc*muu3**2*muu4*(2*k1k2 + k1k3 + k2k3 - p1k3 - p2k3))*sqrtS**4) + 
     -         g2R*muuq*sqrtS**2*(g2Rc*muuqc*
     -             (2*(k1k2 + k1k3)*(k1k2 + k2k3) - (k2k3*muu3**2 + k1k3*muu4**2)*sqrtS**2) + 
     -            g2Lc*muu4*(2*k1k2*p2k1 + k2k3*p2k1 - k1k3*p2k2 + k1k2*p2k3 - 
     -               muu3**2*(2*k1k2 + k1k3 + k2k3 - p2k3)*sqrtS**2))) + 
     -      g1L*g1Lc*(g2R*(2*g2Rc*(k1k2**2*p1p2 + k2k3*p1k1*(p2k1 - p2k2) + 
     -               k1k3*(-p1k1 + p1k2)*p2k2 + 
     -               k1k2*((k1k3 + k2k3)*p1p2 - p1k3*p2k1 - p1k1*p2k2 - p1k2*(p2k1 + p2k3))) + 
     -            (g2Lc*muu4*muuqc*(2*k1k2**2 + 2*k1k2*(k1k3 + k2k3 - p1k1) - k2k3*p1k1 + 
     -                  k1k3*(2*k2k3 - 2*p1k1 + p1k2 - 2*p1k3) - k1k2*p1k3) + 
     -               g2Rc*(k1k3*muu3**2*p1k2 + k1k2*muu3**2*(2*p1k2 + p1k3) - 
     -                  k1k3*muu4**2*p1p2 - k2k3*muu3**2*(p1k1 + p1p2) + k2k3*muu4**2*p2k1 - 
     -                  muu4**2*p1k3*p2k1 - k1k3*muu4**2*p2k2 + muu3**2*p1k3*p2k2 + 
     -                  muu4**2*p1k1*p2k3 - muu3**2*p1k2*p2k3 + k1k2*muu4**2*(2*p2k1 + p2k3)))*
     -             sqrtS**2 - muu4*(g2Lc*muuqc*(k1k3*muu4**2 + muu3**2*(k2k3 - p1k3)) + 
     -               g2Rc*muu3**2*muu4*(2*k1k2 + k1k3 + k2k3 - p1k3 - p2k3))*sqrtS**4) + 
     -         g2L*muuq*sqrtS**2*(g2Lc*muuqc*
     -             (2*(k1k2 + k1k3)*(k1k2 + k2k3) - (k2k3*muu3**2 + k1k3*muu4**2)*sqrtS**2) + 
     -            g2Rc*muu4*(2*k1k2*p2k1 + k2k3*p2k1 - k1k3*p2k2 + k1k2*p2k3 - 
     -               muu3**2*(2*k1k2 + k1k3 + k2k3 - p2k3)*sqrtS**2)))))/
     -  (Denlegk1*Denlegk2c*DenU1*DenU2c)
     
     
! ! ---------------------k2-k1-------------------------------
! 
	  colfac = Nc/2d0*CF
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(-4*gDc*gE*(-(g1L*g1Rc*muu3*sqrtS**2*
     -         (g2R*g2Rc*muuqc*(-2*k1k2**2 - k2k3*p2k1 + k1k3*(-2*k2k3 + p2k2) + 
     -              k1k2*(-2*(k1k3 + k2k3 - p2k2) + p2k3) + 2*k2k3*(p2k2 + p2k3)) + 
     -           g2Rc*muuqc*(g2R*k2k3*muu3**2 + g2L*(2*k1k2 + k1k3 + k2k3)*muu4*muuq + 
     -              g2R*muu4**2*(k1k3 - p2k3))*sqrtS**2 + 
     -           g2L*g2Lc*muuq*(k2k3*p1k1 - 2*k1k2*p1k2 - k1k3*p1k2 - k1k2*p1k3 + 
     -              muu4**2*(2*k1k2 + k1k3 + k2k3 - p1k3)*sqrtS**2) + 
     -           g2Lc*g2R*muu4*(-2*k1k2**2 + k2k3*p1k1 - k2k3*p1p2 - k2k3*p2k1 + p1k3*p2k1 + 
     -              2*k2k3*p2k2 - p1k3*p2k2 + 2*k2k3*p2k3 - p1k1*p2k3 + p1k2*p2k3 - 
     -              2*p1k3*p2k3 + k1k2*
     -               (-2*k1k3 - 2*k2k3 + 2*p1k1 + p1k3 - 2*p1p2 + 2*p2k2 + p2k3) + 
     -              (muu3**2*(k2k3 - p1k3) - muu4**2*p2k3)*sqrtS**2 + 
     -              k1k3*(-2*k2k3 + 2*p1k1 - p1k2 + 2*p1k3 - p1p2 + p2k2 + muu4**2*sqrtS**2)))) -
     -        g1Lc*g1R*muu3*sqrtS**2*(g2R*g2Rc*muuq*
     -          (k2k3*p1k1 - k1k3*p1k2 - k1k2*(2*p1k2 + p1k3)) + 
     -         g2R*muu4*muuq*((2*k1k2 + k1k3 + k2k3)*(g2Rc*muu4 + g2Lc*muuqc) - g2Rc*muu4*p1k3)*
     -          sqrtS**2 + g2L*g2Lc*muuqc*
     -          (-2*k1k2**2 - 2*k1k3*k2k3 - k2k3*p2k1 + k1k3*p2k2 + 2*k2k3*p2k2 + 2*k2k3*p2k3 + 
     -            k1k2*(-2*(k1k3 + k2k3 - p2k2) + p2k3) + 
     -            (k2k3*muu3**2 + muu4**2*(k1k3 - p2k3))*sqrtS**2) + 
     -         g2L*g2Rc*muu4*(-2*k1k2**2 + k2k3*p1k1 - k2k3*p1p2 - k2k3*p2k1 + p1k3*p2k1 + 
     -            2*k2k3*p2k2 - p1k3*p2k2 + 2*k2k3*p2k3 - p1k1*p2k3 + p1k2*p2k3 - 2*p1k3*p2k3 + 
     -            k1k2*(-2*k1k3 - 2*k2k3 + 2*p1k1 + p1k3 - 2*p1p2 + 2*p2k2 + p2k3) + 
     -            (muu3**2*(k2k3 - p1k3) - muu4**2*p2k3)*sqrtS**2 + 
     -            k1k3*(-2*k2k3 + 2*p1k1 - p1k2 + 2*p1k3 - p1p2 + p2k2 + muu4**2*sqrtS**2))) + 
     -      g1R*g1Rc*(g2L*(2*g2Lc*(k1k2**2*p1p2 + k2k3*p1k1*(p2k1 - p2k2) + 
     -               k1k3*(-p1k1 + p1k2)*p2k2 + 
     -               k1k2*((k1k3 + k2k3)*p1p2 - p1k3*p2k1 - p1k1*p2k2 - p1k2*(p2k1 + p2k3))) + 
     -            (g2Rc*muu4*muuqc*(k2k3*p2k1 - k1k3*p2k2 + k1k2*(2*p2k1 + p2k3)) + 
     -               g2Lc*(k1k3*muu3**2*p1k2 + k1k2*muu3**2*(2*p1k2 + p1k3) - 
     -                  k1k3*muu4**2*p1p2 - k2k3*muu3**2*(p1k1 + p1p2) + k2k3*muu4**2*p2k1 - 
     -                  muu4**2*p1k3*p2k1 - k1k3*muu4**2*p2k2 + muu3**2*p1k3*p2k2 + 
     -                  muu4**2*p1k1*p2k3 - muu3**2*p1k2*p2k3 + k1k2*muu4**2*(2*p2k1 + p2k3)))*
     -             sqrtS**2 - muu3**2*muu4*
     -             (g2Rc*muuqc*(2*k1k2 + k1k3 + k2k3 - p2k3) + 
     -               g2Lc*muu4*(2*k1k2 + k1k3 + k2k3 - p1k3 - p2k3))*sqrtS**4) - 
     -         g2R*muuq*sqrtS**2*(g2Rc*muuqc*
     -             (-2*(k1k2 + k1k3)*(k1k2 + k2k3) + (k2k3*muu3**2 + k1k3*muu4**2)*sqrtS**2) + 
     -            g2Lc*muu4*(-2*k1k2**2 + k2k3*p1k1 + k1k2*(-2*(k1k3 + k2k3 - p1k1) + p1k3) + 
     -               muu3**2*(k2k3 - p1k3)*sqrtS**2 + 
     -               k1k3*(-2*k2k3 + 2*p1k1 - p1k2 + 2*p1k3 + muu4**2*sqrtS**2)))) + 
     -      g1L*g1Lc*(g2R*(2*g2Rc*(k1k2**2*p1p2 + k2k3*p1k1*(p2k1 - p2k2) + 
     -               k1k3*(-p1k1 + p1k2)*p2k2 + 
     -               k1k2*((k1k3 + k2k3)*p1p2 - p1k3*p2k1 - p1k1*p2k2 - p1k2*(p2k1 + p2k3))) + 
     -            (g2Lc*muu4*muuqc*(k2k3*p2k1 - k1k3*p2k2 + k1k2*(2*p2k1 + p2k3)) + 
     -               g2Rc*(k1k3*muu3**2*p1k2 + k1k2*muu3**2*(2*p1k2 + p1k3) - 
     -                  k1k3*muu4**2*p1p2 - k2k3*muu3**2*(p1k1 + p1p2) + k2k3*muu4**2*p2k1 - 
     -                  muu4**2*p1k3*p2k1 - k1k3*muu4**2*p2k2 + muu3**2*p1k3*p2k2 + 
     -                  muu4**2*p1k1*p2k3 - muu3**2*p1k2*p2k3 + k1k2*muu4**2*(2*p2k1 + p2k3)))*
     -             sqrtS**2 - muu3**2*muu4*
     -             (g2Lc*muuqc*(2*k1k2 + k1k3 + k2k3 - p2k3) + 
     -               g2Rc*muu4*(2*k1k2 + k1k3 + k2k3 - p1k3 - p2k3))*sqrtS**4) - 
     -         g2L*muuq*sqrtS**2*(g2Lc*muuqc*
     -             (-2*(k1k2 + k1k3)*(k1k2 + k2k3) + (k2k3*muu3**2 + k1k3*muu4**2)*sqrtS**2) + 
     -            g2Rc*muu4*(-2*k1k2**2 + k2k3*p1k1 + k1k2*(-2*(k1k3 + k2k3 - p1k1) + p1k3) + 
     -               muu3**2*(k2k3 - p1k3)*sqrtS**2 + 
     -               k1k3*(-2*k2k3 + 2*p1k1 - p1k2 + 2*p1k3 + muu4**2*sqrtS**2))))))/
     -  (Denlegk1c*Denlegk2*DenU1c*DenU2)

! ! ---------------------k2-k2------------------------------
! 
	  colfac = 0d0
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(8*gE*gEc*(2*(g1R*g1Rc*g2L*g2Lc + g1L*g1Lc*g2R*g2Rc)*k2k3*p2k1*p2k3 - 
     -      (g1Lc*g1R*k2k3*muu3*(2*g2L*g2Rc*muu4*p2k1 + g2R*g2Rc*muuq*(-k1k3 + p2k3) + 
     -            g2L*g2Lc*muuqc*(-k1k3 + p2k3)) + 
     -         g1L*g1Rc*k2k3*muu3*(2*g2Lc*g2R*muu4*p2k1 + g2L*g2Lc*muuq*(-k1k3 + p2k3) + 
     -            g2R*g2Rc*muuqc*(-k1k3 + p2k3)) + 
     -         g1R*g1Rc*(g2L*g2Rc*k2k3*muu4*muuqc*p2k1 + 
     -            g2R*k2k3*muuq*(-(g2Rc*k1k3*muuqc) + g2Lc*muu4*p2k1) + 
     -            g2L*g2Lc*(k1k3*k2k3*(muu2 - muu3)*(muu2 + muu3) + 2*k2k3*muu3**2*p2k3 + 
     -               2*muu4**2*p2k1*(p2k2 + p2k3))) + 
     -         g1L*g1Lc*(g2Lc*g2R*k2k3*muu4*muuqc*p2k1 + 
     -            g2L*k2k3*muuq*(-(g2Lc*k1k3*muuqc) + g2Rc*muu4*p2k1) + 
     -            g2R*g2Rc*(k1k3*k2k3*(muu2 - muu3)*(muu2 + muu3) + 2*k2k3*muu3**2*p2k3 + 
     -               2*muu4**2*p2k1*(p2k2 + p2k3))))*sqrtS**2 + 
     -      muu4*(g1R*g1Rc*(-(g2R*g2Rc*(k1k2 + k1k3)*muu4*muuq*muuqc) + 
     -            g2Lc*g2R*muuq*(k2k3*muu3**2 - muu4**2*p2k1) + 
     -            g2L*g2Rc*muuqc*(k2k3*muu3**2 - muu4**2*p2k1) + 
     -            g2L*g2Lc*muu4*(k1k2*(muu2 - muu3)*(muu2 + muu3) + 
     -               k1k3*(muu2 - muu3)*(muu2 + muu3) + 2*muu3**2*(p2k2 + p2k3))) + 
     -         g1L*g1Lc*(-(g2L*g2Lc*(k1k2 + k1k3)*muu4*muuq*muuqc) + 
     -            g2L*g2Rc*muuq*(k2k3*muu3**2 - muu4**2*p2k1) + 
     -            g2Lc*g2R*muuqc*(k2k3*muu3**2 - muu4**2*p2k1) + 
     -            g2R*g2Rc*muu4*(k1k2*(muu2 - muu3)*(muu2 + muu3) + 
     -               k1k3*(muu2 - muu3)*(muu2 + muu3) + 2*muu3**2*(p2k2 + p2k3))) + 
     -         g1L*g1Rc*muu3*(g2Rc*muuqc*
     -             (g2L*k2k3*muuq + g2R*muu4*(-k1k2 - k1k3 + p2k2 + p2k3)) + 
     -            g2Lc*(g2R*k2k3*(muu2**2 + muu3**2) - 2*g2R*muu4**2*p2k1 + 
     -               g2L*muu4*muuq*(-k1k2 - k1k3 + p2k2 + p2k3))) + 
     -         g1Lc*g1R*muu3*(g2R*muuq*
     -             (g2Lc*k2k3*muuqc + g2Rc*muu4*(-k1k2 - k1k3 + p2k2 + p2k3)) + 
     -            g2L*(g2Rc*k2k3*(muu2**2 + muu3**2) - 2*g2Rc*muu4**2*p2k1 + 
     -               g2Lc*muu4*muuqc*(-k1k2 - k1k3 + p2k2 + p2k3))))*sqrtS**4 + 
     -      muu3*muu4**3*(g1L*(g1Rc*g2Lc*g2R*(muu2**2 + muu3**2) + g1Rc*g2L*g2Rc*muuq*muuqc + 
     -            g1Lc*muu3*(g2L*g2Rc*muuq + g2Lc*g2R*muuqc)) + 
     -         g1R*(g1Lc*g2L*g2Rc*(muu2**2 + muu3**2) + g1Lc*g2Lc*g2R*muuq*muuqc + 
     -            g1Rc*muu3*(g2Lc*g2R*muuq + g2L*g2Rc*muuqc)))*sqrtS**6))/
     -  (Denlegk2*Denlegk2c*DenU2*DenU2c)

! ! ! ---------------------prop-p1------------------------------
! 
	  colfac = 0d0
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(2*gA*gFc*(-2*(g1R*g1Rc*g2L*g2Lc + g1L*g1Lc*g2R*g2Rc)*p2k1*
     -       (k2k3*p1p2 - p1k3*p2k2 - p1k2*p2k3) + 
     -      (g1L*g1Rc*muu3*(g2R*(2*g2Lc*muu4*(-k2k3 + 2*p1k2 + p1k3)*p2k1 + 
     -               g2Rc*muuqc*(k2k3*(-p1k1 + p1p2 - 2*p2k1) + p1k3*(k1k2 - p2k2) + 
     -                  p1k2*(k1k3 + 4*p2k1 - p2k3))) + 
     -            g2L*g2Lc*muuq*(k1k3*p1k2 + k1k2*p1k3 + k2k3*(-p1k1 + p1p2) - p1k3*p2k2 - 
     -               p1k2*p2k3)) + g1Lc*g1R*muu3*
     -          (g2L*(2*g2Rc*muu4*(-k2k3 + 2*p1k2 + p1k3)*p2k1 + 
     -               g2Lc*muuqc*(k2k3*(-p1k1 + p1p2 - 2*p2k1) + p1k3*(k1k2 - p2k2) + 
     -                  p1k2*(k1k3 + 4*p2k1 - p2k3))) + 
     -            g2R*g2Rc*muuq*(k1k3*p1k2 + k1k2*p1k3 + k2k3*(-p1k1 + p1p2) - p1k3*p2k2 - 
     -               p1k2*p2k3)) + g1L*g1Lc*
     -          (g2L*g2Lc*muuq*muuqc*(2*p1k2*p2k1 - k2k3*(p1k1 + p2k1) - 2*p1k1*p2k2 + 
     -               k1k3*(p1k2 + p2k2) + k1k2*(p1k3 + 2*p1p2 - p2k3)) + 
     -            g2L*g2Rc*muu4*muuq*(2*k1k2*p1p2 - k1k3*p1p2 + (-k2k3 + 2*p1k2 + p1k3)*p2k1 + 
     -               k1k3*p2k2 - 2*p1k1*p2k2 - k1k2*p2k3 + p1k1*p2k3) + 
     -            g2R*g2Rc*(-(k1k2*muu2**2*p1k3) + k1k2*muu3**2*p1k3 - 2*k1k2*muuqc**2*p1p2 + 
     -               2*muuqc**2*p1k2*p2k1 + 4*muu4**2*p1p2*p2k1 + 
     -               k2k3*(muu2**2*p1k1 - muu3**2*(p1k1 - 2*p1p2) - muuqc**2*p2k1) + 
     -               2*muuqc**2*p1k1*p2k2 - 2*muu3**2*p1k3*p2k2 - 4*muu1**2*p2k1*p2k2 - 
     -               k1k3*((muu2 - muu3)*(muu2 + muu3)*p1k2 + muuqc**2*p2k2) + 
     -               k1k2*muuqc**2*p2k3 - 2*muu3**2*p1k2*p2k3 - 2*muu4**2*p2k1*p2k3) + 
     -            g2Lc*g2R*muu4*muuqc*(k1k3*(p1p2 - p2k2) + 
     -               p2k1*(-k2k3 + 2*p1k2 + p1k3 + 4*p1p2 - 2*p2k3) + p1k1*(2*p2k2 - p2k3) + 
     -               k1k2*(-2*p1p2 + p2k3))) + 
     -         g1R*g1Rc*(g2R*g2Rc*muuq*muuqc*
     -             (2*p1k2*p2k1 - k2k3*(p1k1 + p2k1) - 2*p1k1*p2k2 + k1k3*(p1k2 + p2k2) + 
     -               k1k2*(p1k3 + 2*p1p2 - p2k3)) + 
     -            g2Lc*g2R*muu4*muuq*(2*k1k2*p1p2 - k1k3*p1p2 + (-k2k3 + 2*p1k2 + p1k3)*p2k1 + 
     -               k1k3*p2k2 - 2*p1k1*p2k2 - k1k2*p2k3 + p1k1*p2k3) + 
     -            g2L*g2Lc*(-(k1k2*muu2**2*p1k3) + k1k2*muu3**2*p1k3 - 2*k1k2*muuqc**2*p1p2 + 
     -               2*muuqc**2*p1k2*p2k1 + 4*muu4**2*p1p2*p2k1 + 
     -               k2k3*(muu2**2*p1k1 - muu3**2*(p1k1 - 2*p1p2) - muuqc**2*p2k1) + 
     -               2*muuqc**2*p1k1*p2k2 - 2*muu3**2*p1k3*p2k2 - 4*muu1**2*p2k1*p2k2 - 
     -               k1k3*((muu2 - muu3)*(muu2 + muu3)*p1k2 + muuqc**2*p2k2) + 
     -               k1k2*muuqc**2*p2k3 - 2*muu3**2*p1k2*p2k3 - 2*muu4**2*p2k1*p2k3) + 
     -            g2L*g2Rc*muu4*muuqc*(k1k3*(p1p2 - p2k2) + 
     -               p2k1*(-k2k3 + 2*p1k2 + p1k3 + 4*p1p2 - 2*p2k3) + p1k1*(2*p2k2 - p2k3) + 
     -               k1k2*(-2*p1p2 + p2k3))))*sqrtS**2 + 
     -      (g1L*g1Rc*muu3*(g2L*g2Rc*muu4*muuq*muuqc*
     -             (-k1k3 + k2k3 + 2*p1k1 - 2*p1k2 - p1k3 - 2*p1p2 + p2k3) + 
     -            g2R*g2Rc*muuqc*((muu2**2 + muu3**2)*(k2k3 - 2*p1k2) + 
     -               2*muu1**2*(-k1k2 + p2k2) + muu4**2*(-k1k3 + 2*p1k1 - 2*p1p2 + p2k3)) + 
     -            g2L*g2Lc*muuq*(muuqc**2*(k2k3 - 2*p1k2) + 2*muu1**2*(-k1k2 + p2k2) + 
     -               muu4**2*(-k1k3 + 2*p1k1 - 2*p1p2 + p2k3)) + 
     -            g2Lc*g2R*muu4*((muu2**2 + muu3**2)*(k2k3 - 2*p1k2 - p1k3) - 4*muu1**2*p2k1 + 
     -               muuqc**2*(-k1k3 + 2*p1k1 - 2*p1p2 + p2k3))) + 
     -         g1Lc*g1R*muu3*(g2Lc*g2R*muu4*muuq*muuqc*
     -             (-k1k3 + k2k3 + 2*p1k1 - 2*p1k2 - p1k3 - 2*p1p2 + p2k3) + 
     -            g2L*g2Lc*muuqc*((muu2**2 + muu3**2)*(k2k3 - 2*p1k2) + 
     -               2*muu1**2*(-k1k2 + p2k2) + muu4**2*(-k1k3 + 2*p1k1 - 2*p1p2 + p2k3)) + 
     -            g2R*g2Rc*muuq*(muuqc**2*(k2k3 - 2*p1k2) + 2*muu1**2*(-k1k2 + p2k2) + 
     -               muu4**2*(-k1k3 + 2*p1k1 - 2*p1p2 + p2k3)) + 
     -            g2L*g2Rc*muu4*((muu2**2 + muu3**2)*(k2k3 - 2*p1k2 - p1k3) - 4*muu1**2*p2k1 + 
     -               muuqc**2*(-k1k3 + 2*p1k1 - 2*p1p2 + p2k3))) + 
     -         g1R*g1Rc*(g2R*g2Rc*muuq*muuqc*
     -             (-2*k1k2*muu1**2 - muu4**2*(k1k3 - 2*p1k1) + muu3**2*(k2k3 - 2*p1k2)) + 
     -            g2Lc*g2R*muu4*muuq*(-(muuqc**2*(k1k3 - 2*p1k1)) + 
     -               muu3**2*(k2k3 - 2*p1k2 - p1k3) - 2*muu1**2*p2k1) + 
     -            g2L*g2Rc*muu4*muuqc*(k1k3*(muu2 - muu3)*(muu2 + muu3) - 
     -               2*(muu2**2*p1k1 + muu1**2*p2k1) + 
     -               muu3**2*(k2k3 + 2*p1k1 - 2*p1k2 - p1k3 - 4*p1p2 + 2*p2k3)) + 
     -            g2L*g2Lc*(2*k1k2*muu1**2*(muu2 - muu3)*(muu2 + muu3) + 
     -               k1k3*(muu2 - muu3)*(muu2 + muu3)*muu4**2 - 2*muu2**2*muu4**2*p1k1 + 
     -               muu3**2*(muuqc**2*(k2k3 - 2*p1k2) + 4*muu1**2*p2k2 + 
     -                  2*muu4**2*(p1k1 - 2*p1p2 + p2k3)))) + 
     -         g1L*g1Lc*(g2L*g2Lc*muuq*muuqc*
     -             (-2*k1k2*muu1**2 - muu4**2*(k1k3 - 2*p1k1) + muu3**2*(k2k3 - 2*p1k2)) + 
     -            g2L*g2Rc*muu4*muuq*(-(muuqc**2*(k1k3 - 2*p1k1)) + 
     -               muu3**2*(k2k3 - 2*p1k2 - p1k3) - 2*muu1**2*p2k1) + 
     -            g2Lc*g2R*muu4*muuqc*(k1k3*(muu2 - muu3)*(muu2 + muu3) - 
     -               2*(muu2**2*p1k1 + muu1**2*p2k1) + 
     -               muu3**2*(k2k3 + 2*p1k1 - 2*p1k2 - p1k3 - 4*p1p2 + 2*p2k3)) + 
     -            g2R*g2Rc*(2*k1k2*muu1**2*(muu2 - muu3)*(muu2 + muu3) + 
     -               k1k3*(muu2 - muu3)*(muu2 + muu3)*muu4**2 - 2*muu2**2*muu4**2*p1k1 + 
     -               muu3**2*(muuqc**2*(k2k3 - 2*p1k2) + 4*muu1**2*p2k2 + 
     -                  2*muu4**2*(p1k1 - 2*p1p2 + p2k3)))))*sqrtS**4 + 
     -      2*muu1**2*muu3*muu4*(g1L*(g1Rc*g2Lc*g2R*(muu2**2 + muu3**2) + 
     -            g1Rc*g2L*g2Rc*muuq*muuqc + g1Lc*muu3*(g2L*g2Rc*muuq + g2Lc*g2R*muuqc)) + 
     -         g1R*(g1Lc*g2L*g2Rc*(muu2**2 + muu3**2) + g1Lc*g2Lc*g2R*muuq*muuqc + 
     -            g1Rc*muu3*(g2Lc*g2R*muuq + g2L*g2Rc*muuqc)))*sqrtS**6))/
     -  (Denlegp1*DenU1*DenU1c*DenU2c)

! ! ---------------------prop-p2------------------------------
! 
	  colfac = 0d0
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu + colfac*2d0*(2*gB*gFc*(g1L*g1Rc*muu3*sqrtS**2*
     -       (g2L*muuq*(g2Lc*(k1k2*p1k3 - 2*k1k2*p1p2 + k2k3*(-p1k1 + p1p2) + 2*p1k1*p2k2 - 
     -               p1k3*p2k2 + p1k2*(k1k3 - 2*p2k1 - p2k3)) + 
     -            (g2Rc*muu4*muuqc*(-k1k3 + k2k3 - p1k3 + 2*(p1p2 + p2k1 - p2k2) + p2k3) + 
     -               g2Lc*(2*muu2**2*p1k2 + muuqc**2*(k2k3 - 2*p2k2) + 
     -                  muu4**2*(-k1k3 + 2*p2k1 + p2k3)))*sqrtS**2 - 
     -            2*muu2**2*muu4*(g2Lc*muu4 + g2Rc*muuqc)*sqrtS**4) + 
     -         g2R*(2*g2Lc*k1k3*muu4*p1k2 + 
     -            g2Rc*muuqc*(k2k3*p1k1 + k1k3*p1k2 - k1k2*p1k3 + 2*p1k2*p1k3 + 2*k1k2*p1p2 - 
     -               k2k3*p1p2 - 4*p1k2*p1p2 - 2*p1k2*p2k1 - 2*p1k1*p2k2 + p1k3*p2k2 - 
     -               p1k2*p2k3 + (k2k3*(-muu1**2 + muu4**2) + 2*(muu2**2*p1k2 + muu1**2*p2k2) + 
     -                  muu4**2*(-k1k3 - 2*p1k3 + 4*p1p2 + 2*p2k1 - 2*p2k2 + p2k3))*sqrtS**2 - 
     -               2*muu2**2*muu4**2*sqrtS**4) + 
     -            g2Lc*muu4*(-2*p1k2*(2*p2k1 + p2k3 - 2*muu2**2*sqrtS**2) + 
     -               sqrtS**2*(muuqc**2*(k2k3 - p1k3 + 2*p1p2 - 2*p2k2) - 
     -                  (muu1**2 + muu4**2)*(k1k3 - 2*p2k1 - p2k3 + 2*muu2**2*sqrtS**2))))) + 
     -      g1Lc*g1R*muu3*sqrtS**2*(g2R*muuq*
     -          (g2Rc*(k1k2*p1k3 - 2*k1k2*p1p2 + k2k3*(-p1k1 + p1p2) + 2*p1k1*p2k2 - 
     -               p1k3*p2k2 + p1k2*(k1k3 - 2*p2k1 - p2k3)) + 
     -            (g2Lc*muu4*muuqc*(-k1k3 + k2k3 - p1k3 + 2*(p1p2 + p2k1 - p2k2) + p2k3) + 
     -               g2Rc*(2*muu2**2*p1k2 + muuqc**2*(k2k3 - 2*p2k2) + 
     -                  muu4**2*(-k1k3 + 2*p2k1 + p2k3)))*sqrtS**2 - 
     -            2*muu2**2*muu4*(g2Rc*muu4 + g2Lc*muuqc)*sqrtS**4) + 
     -         g2L*(2*g2Rc*k1k3*muu4*p1k2 + 
     -            g2Lc*muuqc*(k2k3*p1k1 + k1k3*p1k2 - k1k2*p1k3 + 2*p1k2*p1k3 + 2*k1k2*p1p2 - 
     -               k2k3*p1p2 - 4*p1k2*p1p2 - 2*p1k2*p2k1 - 2*p1k1*p2k2 + p1k3*p2k2 - 
     -               p1k2*p2k3 + (k2k3*(-muu1**2 + muu4**2) + 2*(muu2**2*p1k2 + muu1**2*p2k2) + 
     -                  muu4**2*(-k1k3 - 2*p1k3 + 4*p1p2 + 2*p2k1 - 2*p2k2 + p2k3))*sqrtS**2 - 
     -               2*muu2**2*muu4**2*sqrtS**4) + 
     -            g2Rc*muu4*(-2*p1k2*(2*p2k1 + p2k3 - 2*muu2**2*sqrtS**2) + 
     -               sqrtS**2*(muuqc**2*(k2k3 - p1k3 + 2*p1p2 - 2*p2k2) - 
     -                  (muu1**2 + muu4**2)*(k1k3 - 2*p2k1 - p2k3 + 2*muu2**2*sqrtS**2))))) + 
     -      g1R*g1Rc*(g2R*muuq*sqrtS**2*
     -          (g2Rc*muuqc*(-(k2k3*p1k1) + k1k3*p1k2 + k1k2*p1k3 - 2*k1k2*p1p2 - k2k3*p2k1 - 
     -               2*p1k2*p2k1 + k1k3*p2k2 + 2*p1k1*p2k2 - k1k2*p2k3 + 
     -               (2*k1k2*muu2**2 - muu4**2*(k1k3 - 2*p2k1) + muu3**2*(k2k3 - 2*p2k2))*
     -                sqrtS**2) + g2Lc*muu4*
     -             (-(k1k3*p1p2) - k2k3*p2k1 + p1k3*p2k1 + k1k3*p2k2 - k1k2*p2k3 + p1k1*p2k3 + 
     -               (2*muu2**2*(k1k2 - p1k1) - muuqc**2*(k1k3 - 2*p2k1) + 
     -                  muu3**2*(k2k3 - p1k3 + 2*p1p2 - 2*p2k2))*sqrtS**2)) + 
     -         g2L*(g2Rc*muu4*muuqc*sqrtS**2*
     -             (2*k1k3*p1k2 - k1k3*p1p2 - k2k3*p2k1 - 4*p1k2*p2k1 + p1k3*p2k1 + k1k3*p2k2 - 
     -               k1k2*p2k3 + p1k1*p2k3 - 
     -               (2*muu2**2*(-k1k2 + p1k1) + (muu1**2 + muu4**2)*(k1k3 - 2*p2k1) + 
     -                  muu3**2*(-k2k3 + p1k3 - 2*p1p2 + 2*p2k2))*sqrtS**2) + 
     -            g2Lc*(-2*p1k2*(p1k3*p2k1 + p1k1*p2k3) + 
     -               2*p1k2*(2*muu2**2*p1k1 + muu3**2*(p1k3 - 2*p1p2) - muuqc**2*p2k1)*
     -                sqrtS**2 + sqrtS**2*
     -                (-(k1k2*muuqc**2*(p1k3 - 2*p1p2)) - 2*muuqc**2*p1k1*p2k2 + 
     -                  2*muu1**2*muu3**2*p2k2*sqrtS**2 + 
     -                  k1k2*(muu1 - muu4)*(muu1 + muu4)*(p2k3 - 2*muu2**2*sqrtS**2) + 
     -                  2*muu4**2*(p1k3*p2k1 + p1k1*p2k3 - 
     -                     (2*muu2**2*p1k1 - muuqc**2*p2k1 + muu3**2*(p1k3 - 2*p1p2 + p2k2))*
     -                      sqrtS**2) + 
     -                  k2k3*(muuqc**2*p1k1 + 
     -                     (muu1 - muu4)*(muu1 + muu4)*(p2k1 - muu3**2*sqrtS**2))) + 
     -               k1k3*(p1k2*(2*p1p2 + muuqc**2*sqrtS**2) - 
     -                  sqrtS**2*(muu1**2*p2k2 + muu4**2*(2*p1p2 - p2k2 + muuqc**2*sqrtS**2)))))
     -         ) + g1L*g1Lc*(g2L*muuq*sqrtS**2*
     -          (g2Lc*muuqc*(-(k2k3*p1k1) + k1k3*p1k2 + k1k2*p1k3 - 2*k1k2*p1p2 - k2k3*p2k1 - 
     -               2*p1k2*p2k1 + k1k3*p2k2 + 2*p1k1*p2k2 - k1k2*p2k3 + 
     -               (2*k1k2*muu2**2 - muu4**2*(k1k3 - 2*p2k1) + muu3**2*(k2k3 - 2*p2k2))*
     -                sqrtS**2) + g2Rc*muu4*
     -             (-(k1k3*p1p2) - k2k3*p2k1 + p1k3*p2k1 + k1k3*p2k2 - k1k2*p2k3 + p1k1*p2k3 + 
     -               (2*muu2**2*(k1k2 - p1k1) - muuqc**2*(k1k3 - 2*p2k1) + 
     -                  muu3**2*(k2k3 - p1k3 + 2*p1p2 - 2*p2k2))*sqrtS**2)) + 
     -         g2R*(g2Lc*muu4*muuqc*sqrtS**2*
     -             (2*k1k3*p1k2 - k1k3*p1p2 - k2k3*p2k1 - 4*p1k2*p2k1 + p1k3*p2k1 + k1k3*p2k2 - 
     -               k1k2*p2k3 + p1k1*p2k3 - 
     -               (2*muu2**2*(-k1k2 + p1k1) + (muu1**2 + muu4**2)*(k1k3 - 2*p2k1) + 
     -                  muu3**2*(-k2k3 + p1k3 - 2*p1p2 + 2*p2k2))*sqrtS**2) + 
     -            g2Rc*(-2*p1k2*(p1k3*p2k1 + p1k1*p2k3) + 
     -               2*p1k2*(2*muu2**2*p1k1 + muu3**2*(p1k3 - 2*p1p2) - muuqc**2*p2k1)*
     -                sqrtS**2 + sqrtS**2*
     -                (-(k1k2*muuqc**2*(p1k3 - 2*p1p2)) - 2*muuqc**2*p1k1*p2k2 + 
     -                  2*muu1**2*muu3**2*p2k2*sqrtS**2 + 
     -                  k1k2*(muu1 - muu4)*(muu1 + muu4)*(p2k3 - 2*muu2**2*sqrtS**2) + 
     -                  2*muu4**2*(p1k3*p2k1 + p1k1*p2k3 - 
     -                     (2*muu2**2*p1k1 - muuqc**2*p2k1 + muu3**2*(p1k3 - 2*p1p2 + p2k2))*
     -                      sqrtS**2) + 
     -                  k2k3*(muuqc**2*p1k1 + 
     -                     (muu1 - muu4)*(muu1 + muu4)*(p2k1 - muu3**2*sqrtS**2))) + 
     -               k1k3*(p1k2*(2*p1p2 + muuqc**2*sqrtS**2) - 
     -                  sqrtS**2*(muu1**2*p2k2 + muu4**2*(2*p1p2 - p2k2 + muuqc**2*sqrtS**2)))))
     -         )))/(Denlegp2*DenU1c*DenU2*DenU2c)
!      
! ! ---------------------prop-k1------------------------------
! 
	  colfac = 0d0
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu+ colfac*2d0*(4*gD*gFc*(4*(g1R*g1Rc*g2L*g2Lc + g1L*g1Lc*g2R*g2Rc)*p1k2*(p1k1 + p1k3)*p2k1 + 
     -      (g1R*g1Rc*(2*g2Lc*g2R*muu4*muuq*(k1k2 + k2k3 - p1k1 - p1k3)*p2k1 + 
     -            g2R*g2Rc*muuq*muuqc*(k2k3*p1k1 + k1k3*p1k2 - k1k2*p1k3 + 
     -               2*(k1k2 + k2k3)*p2k1) + 
     -            2*g2L*g2Rc*muu4*muuqc*(k1k3*p1k2 + (k1k2 + k2k3 - p1k1 - p1k3)*p2k1) + 
     -            g2L*g2Lc*(k1k3*muuqc**2*p1k2 - 2*muu3**2*p1k2*(p1k1 + 2*p1k3 + p1p2) - 
     -               4*muu4**2*(p1k1 + p1k3)*p2k1 + 
     -               k1k2*(muuqc**2*p1k3 - 2*muu1**2*p2k1 + 2*muu4**2*p2k1) - 
     -               k2k3*(muuqc**2*p1k1 + 2*(muu1 - muu4)*(muu1 + muu4)*p2k1))) + 
     -         g1L*g1Lc*(2*g2L*g2Rc*muu4*muuq*(k1k2 + k2k3 - p1k1 - p1k3)*p2k1 + 
     -            g2L*g2Lc*muuq*muuqc*(k2k3*p1k1 + k1k3*p1k2 - k1k2*p1k3 + 
     -               2*(k1k2 + k2k3)*p2k1) + 
     -            2*g2Lc*g2R*muu4*muuqc*(k1k3*p1k2 + (k1k2 + k2k3 - p1k1 - p1k3)*p2k1) + 
     -            g2R*g2Rc*(k1k3*muuqc**2*p1k2 - 2*muu3**2*p1k2*(p1k1 + 2*p1k3 + p1p2) - 
     -               4*muu4**2*(p1k1 + p1k3)*p2k1 + 
     -               k1k2*(muuqc**2*p1k3 - 2*muu1**2*p2k1 + 2*muu4**2*p2k1) - 
     -               k2k3*(muuqc**2*p1k1 + 2*(muu1 - muu4)*(muu1 + muu4)*p2k1))) + 
     -         g1Lc*g1R*muu3*(g2R*g2Rc*muuq*
     -             (k2k3*(p1k1 - p1p2) - k1k2*(p1k3 + p1p2) + (p1k1 + p1k3)*p2k2 + 
     -               p1k2*(k1k3 + p2k1 - p2k3)) + 
     -            g2L*(2*g2Rc*muu4*p1k2*(k1k3 + p2k1 - p2k3) + 
     -               g2Lc*muuqc*(k1k3*p1k2 - 4*p1k2*p1k3 + k2k3*(-p1k1 + p1p2) + 
     -                  k1k2*(p1k3 + p1p2) + p1k2*p2k1 - p1k3*p2k2 - p1k1*(2*p1k2 + p2k2) - 
     -                  p1k2*p2k3))) + 
     -         g1L*g1Rc*muu3*(g2L*g2Lc*muuq*
     -             (k2k3*(p1k1 - p1p2) - k1k2*(p1k3 + p1p2) + (p1k1 + p1k3)*p2k2 + 
     -               p1k2*(k1k3 + p2k1 - p2k3)) + 
     -            g2R*(2*g2Lc*muu4*p1k2*(k1k3 + p2k1 - p2k3) + 
     -               g2Rc*muuqc*(k1k3*p1k2 - 4*p1k2*p1k3 + k2k3*(-p1k1 + p1p2) + 
     -                  k1k2*(p1k3 + p1p2) + p1k2*p2k1 - p1k3*p2k2 - p1k1*(2*p1k2 + p2k2) - 
     -                  p1k2*p2k3))))*sqrtS**2 - 
     -      (g1L*g1Lc*(g2R*g2Rc*(k1k2*muu3**2*(-muu1**2 + muu4**2) + 
     -               2*k2k3*muu3**2*(-muu1**2 + muu4**2) + k1k3*muu4**2*muuqc**2 - 
     -               2*muu3**2*muu4**2*p1k1 + muu3**2*muuqc**2*p1k2 - 4*muu3**2*muu4**2*p1k3 - 
     -               2*muu3**2*muu4**2*p1p2 + muu3**2*(-muu1 + muu4)*(muu1 + muu4)*p2k2) + 
     -            g2L*g2Lc*muuq*muuqc*(k1k3*muu4**2 + muu3**2*(k1k2 + 2*k2k3 + p1k2 + p2k2)) + 
     -            g2L*g2Rc*muu4*muuq*(k1k3*muuqc**2 + 
     -               muu3**2*(k1k2 + 2*k2k3 - p1k1 - 2*p1k3 - p1p2 + p2k2)) + 
     -            g2Lc*g2R*muu4*muuqc*(k1k3*(muu1**2 + muu4**2) + 
     -               muu3**2*(k1k2 + 2*k2k3 - p1k1 + 2*p1k2 - 2*p1k3 - p1p2 + p2k2))) + 
     -         g1R*g1Rc*(g2L*g2Lc*(k1k2*muu3**2*(-muu1**2 + muu4**2) + 
     -               2*k2k3*muu3**2*(-muu1**2 + muu4**2) + k1k3*muu4**2*muuqc**2 - 
     -               2*muu3**2*muu4**2*p1k1 + muu3**2*muuqc**2*p1k2 - 4*muu3**2*muu4**2*p1k3 - 
     -               2*muu3**2*muu4**2*p1p2 + muu3**2*(-muu1 + muu4)*(muu1 + muu4)*p2k2) + 
     -            g2R*g2Rc*muuq*muuqc*(k1k3*muu4**2 + muu3**2*(k1k2 + 2*k2k3 + p1k2 + p2k2)) + 
     -            g2Lc*g2R*muu4*muuq*(k1k3*muuqc**2 + 
     -               muu3**2*(k1k2 + 2*k2k3 - p1k1 - 2*p1k3 - p1p2 + p2k2)) + 
     -            g2L*g2Rc*muu4*muuqc*(k1k3*(muu1**2 + muu4**2) + 
     -               muu3**2*(k1k2 + 2*k2k3 - p1k1 + 2*p1k2 - 2*p1k3 - p1p2 + p2k2))) + 
     -         g1Lc*g1R*muu3*(g2R*g2Rc*muuq*
     -             ((k1k2 + 2*k2k3)*muuqc**2 + muu3**2*p1k2 + muu4**2*(k1k3 + p2k1 - p2k3)) + 
     -            g2L*g2Rc*muu4*(2*muu3**2*p1k2 + muuqc**2*(k1k2 + 2*k2k3 - p1k1 - 2*p1k3) + 
     -               (muu1**2 + muu4**2)*(k1k3 + p2k1 - p2k3)) + 
     -            g2L*g2Lc*muuqc*(k1k2*(-muu1**2 + muu4**2) + 2*k2k3*(-muu1**2 + muu4**2) + 
     -               muu3**2*p1k2 + muu4**2*(k1k3 - 2*p1k1 - 4*p1k3 + p2k1 - p2k3)) + 
     -            g2Lc*g2R*muu4*muuq*muuqc*(k1k2 + k1k3 + 2*k2k3 - p1k1 - 2*p1k3 + p2k1 - p2k3))
     -           + g1L*g1Rc*muu3*(g2L*g2Lc*muuq*
     -             ((k1k2 + 2*k2k3)*muuqc**2 + muu3**2*p1k2 + muu4**2*(k1k3 + p2k1 - p2k3)) + 
     -            g2Lc*g2R*muu4*(2*muu3**2*p1k2 + muuqc**2*(k1k2 + 2*k2k3 - p1k1 - 2*p1k3) + 
     -               (muu1**2 + muu4**2)*(k1k3 + p2k1 - p2k3)) + 
     -            g2R*g2Rc*muuqc*(k1k2*(-muu1**2 + muu4**2) + 2*k2k3*(-muu1**2 + muu4**2) + 
     -               muu3**2*p1k2 + muu4**2*(k1k3 - 2*p1k1 - 4*p1k3 + p2k1 - p2k3)) + 
     -            g2L*g2Rc*muu4*muuq*muuqc*(k1k2 + k1k3 + 2*k2k3 - p1k1 - 2*p1k3 + p2k1 - p2k3))
     -         )*sqrtS**4 + muu3**2*muu4*
     -       (g1L*(g1Rc*muu3 + g1Lc*muuqc)*
     -          (g2Lc*g2R*(muu1**2 + muu4**2) + g2L*g2Lc*muu4*muuq + 
     -            g2Rc*(g2R*muu4 + g2L*muuq)*muuqc) + 
     -         g1R*(g1Lc*muu3 + g1Rc*muuqc)*
     -          (g2L*g2Rc*(muu1**2 + muu4**2) + g2L*g2Lc*muu4*muuqc + 
     -            g2R*muuq*(g2Rc*muu4 + g2Lc*muuqc)))*sqrtS**6))/(Denlegk1*DenU1*DenU1c*DenU2c)

! ! ---------------------prop-k2------------------------------
! 
	  colfac = 0d0
! 
	  stst2QQ_M2realemissionUGLUNeu = stst2QQ_M2realemissionUGLUNeu+colfac*2d0*(4*gE*gFc*(-4*(g1R*g1Rc*g2L*g2Lc + g1L*g1Lc*g2R*g2Rc)*p1k2*p2k1*(p2k2 + p2k3) + 
     -      (-2*g1Lc*g1R*muu3*(g2L*g2Rc*muu4*(k2k3 + p1k2 - p1k3)*p2k1 + 
     -            g2R*g2Rc*muuq*p1k2*(k1k2 + k1k3 - p2k2 - p2k3) + 
     -            g2L*g2Lc*muuqc*(k1k2*p1k2 + k1k3*p1k2 + k2k3*p2k1 - p1k2*(p2k2 + p2k3))) - 
     -         2*g1L*g1Rc*muu3*(g2Lc*g2R*muu4*(k2k3 + p1k2 - p1k3)*p2k1 + 
     -            g2L*g2Lc*muuq*p1k2*(k1k2 + k1k3 - p2k2 - p2k3) + 
     -            g2R*g2Rc*muuqc*(k1k2*p1k2 + k1k3*p1k2 + k2k3*p2k1 - p1k2*(p2k2 + p2k3))) + 
     -         g1L*g1Lc*(-(g2L*muuq*(g2Lc*muuqc*
     -                  (2*k1k2*p1k2 + 2*k1k3*p1k2 + k2k3*p2k1 + k1k3*p2k2 - k1k2*p2k3) + 
     -                 g2Rc*muu4*((k2k3 + p1k2 - p1k3)*p2k1 + k1k3*(-p1p2 + p2k2) - 
     -                    k1k2*(p1p2 + p2k3) + p1k1*(p2k2 + p2k3)))) + 
     -            g2R*(g2Rc*(2*k1k2*(muu2 - muu3)*(muu2 + muu3)*p1k2 + 
     -                  2*k1k3*(muu2 - muu3)*(muu2 + muu3)*p1k2 - k2k3*muuqc**2*p2k1 + 
     -                  k1k3*muuqc**2*p2k2 - k1k2*muuqc**2*p2k3 + 
     -                  4*muu3**2*p1k2*(p2k2 + p2k3) + 2*muu4**2*p2k1*(p1p2 + p2k2 + 2*p2k3)) + 
     -               g2Lc*muu4*muuqc*(k1k3*(-p1p2 + p2k2) - k1k2*(p1p2 + p2k3) + 
     -                  p1k1*(p2k2 + p2k3) + p2k1*(-k2k3 - p1k2 + p1k3 + 2*p2k2 + 4*p2k3)))) + 
     -         g1R*g1Rc*(-(g2R*muuq*(g2Rc*muuqc*
     -                  (2*k1k2*p1k2 + 2*k1k3*p1k2 + k2k3*p2k1 + k1k3*p2k2 - k1k2*p2k3) + 
     -                 g2Lc*muu4*((k2k3 + p1k2 - p1k3)*p2k1 + k1k3*(-p1p2 + p2k2) - 
     -                    k1k2*(p1p2 + p2k3) + p1k1*(p2k2 + p2k3)))) + 
     -            g2L*(g2Lc*(2*k1k2*(muu2 - muu3)*(muu2 + muu3)*p1k2 + 
     -                  2*k1k3*(muu2 - muu3)*(muu2 + muu3)*p1k2 - k2k3*muuqc**2*p2k1 + 
     -                  k1k3*muuqc**2*p2k2 - k1k2*muuqc**2*p2k3 + 
     -                  4*muu3**2*p1k2*(p2k2 + p2k3) + 2*muu4**2*p2k1*(p1p2 + p2k2 + 2*p2k3)) + 
     -               g2Rc*muu4*muuqc*(k1k3*(-p1p2 + p2k2) - k1k2*(p1p2 + p2k3) + 
     -                  p1k1*(p2k2 + p2k3) + p2k1*(-k2k3 - p1k2 + p1k3 + 2*p2k2 + 4*p2k3)))))*
     -       sqrtS**2 + (g1Lc*g1R*muu3*
     -          (g2L*(g2Rc*muu4*((muu2**2 + muu3**2)*(k2k3 + p1k2 - p1k3) + 2*muu4**2*p2k1 + 
     -                  muuqc**2*(k1k2 + 2*k1k3 - p2k2 - 2*p2k3)) + 
     -               g2Lc*muuqc*(k2k3*(muu2**2 + muu3**2) + 
     -                  muu4**2*(k1k2 + 2*k1k3 + p1k1 - p1p2 + 2*p2k1 - p2k2 - 2*p2k3))) + 
     -            g2R*muuq*(g2Rc*k2k3*muuqc**2 + 
     -               g2Lc*muu4*muuqc*(k1k2 + 2*k1k3 + k2k3 + p1k2 - p1k3 - p2k2 - 2*p2k3) + 
     -               g2Rc*muu4**2*(k1k2 + 2*k1k3 + p1k1 - p1p2 - p2k2 - 2*p2k3))) + 
     -         g1L*g1Rc*muu3*(g2Lc*(g2R*muu4*
     -                ((muu2**2 + muu3**2)*(k2k3 + p1k2 - p1k3) + 2*muu4**2*p2k1 + 
     -                  muuqc**2*(k1k2 + 2*k1k3 - p2k2 - 2*p2k3)) + 
     -               g2L*muuq*(k2k3*muuqc**2 + 
     -                  muu4**2*(k1k2 + 2*k1k3 + p1k1 - p1p2 - p2k2 - 2*p2k3))) + 
     -            g2Rc*muuqc*(g2R*k2k3*(muu2**2 + muu3**2) + 
     -               g2L*muu4*muuq*(k1k2 + 2*k1k3 + k2k3 + p1k2 - p1k3 - p2k2 - 2*p2k3) + 
     -               g2R*muu4**2*(k1k2 + 2*k1k3 + p1k1 - p1p2 + 2*p2k1 - p2k2 - 2*p2k3))) + 
     -         g1R*g1Rc*(g2Lc*g2R*muu4*muuq*
     -             ((k1k2 + 2*k1k3)*muuqc**2 + muu3**2*(k2k3 + p1k2 - p1k3) + muu4**2*p2k1) + 
     -            g2R*g2Rc*muuq*muuqc*(k2k3*muu3**2 + muu4**2*(k1k2 + 2*k1k3 + p1k1 + p2k1)) + 
     -            g2L*g2Rc*muu4*muuqc*(k1k2*(-muu2**2 + muu3**2) + 
     -               2*k1k3*(-muu2**2 + muu3**2) + muu4**2*p2k1 + 
     -               muu3**2*(k2k3 + p1k2 - p1k3 - 2*p2k2 - 4*p2k3)) + 
     -            g2L*g2Lc*(k1k2*(-muu2**2 + muu3**2)*muu4**2 + 
     -               2*k1k3*(-muu2**2 + muu3**2)*muu4**2 + k2k3*muu3**2*muuqc**2 - 
     -               muu2**2*muu4**2*p1k1 + muu3**2*muu4**2*p1k1 - 2*muu3**2*muu4**2*p1p2 + 
     -               muu4**2*muuqc**2*p2k1 - 2*muu3**2*muu4**2*p2k2 - 4*muu3**2*muu4**2*p2k3)) +
     -          g1L*g1Lc*(g2L*g2Rc*muu4*muuq*
     -             ((k1k2 + 2*k1k3)*muuqc**2 + muu3**2*(k2k3 + p1k2 - p1k3) + muu4**2*p2k1) + 
     -            g2L*g2Lc*muuq*muuqc*(k2k3*muu3**2 + muu4**2*(k1k2 + 2*k1k3 + p1k1 + p2k1)) + 
     -            g2Lc*g2R*muu4*muuqc*(k1k2*(-muu2**2 + muu3**2) + 
     -               2*k1k3*(-muu2**2 + muu3**2) + muu4**2*p2k1 + 
     -               muu3**2*(k2k3 + p1k2 - p1k3 - 2*p2k2 - 4*p2k3)) + 
     -            g2R*g2Rc*(k1k2*(-muu2**2 + muu3**2)*muu4**2 + 
     -               2*k1k3*(-muu2**2 + muu3**2)*muu4**2 + k2k3*muu3**2*muuqc**2 - 
     -               muu2**2*muu4**2*p1k1 + muu3**2*muu4**2*p1k1 - 2*muu3**2*muu4**2*p1p2 + 
     -               muu4**2*muuqc**2*p2k1 - 2*muu3**2*muu4**2*p2k2 - 4*muu3**2*muu4**2*p2k3)))*
     -       sqrtS**4 - muu3*muu4**2*(g1Lc*g1R*g2L*(muu2**2 + muu3**2)*
     -          (g2Rc*muu4 + g2Lc*muuqc) + g1Lc*g1R*g2R*muuq*muuqc*(g2Lc*muu4 + g2Rc*muuqc) + 
     -         g1R*g1Rc*muu3*(g2Lc*g2R*muu4*muuq + g2L*g2Rc*muu4*muuqc + g2R*g2Rc*muuq*muuqc + 
     -            g2L*g2Lc*muuqc**2) + 
     -         g1L*g1Rc*(g2Lc*g2R*(muu2**2 + muu3**2)*muu4 + 
     -            g2Rc*(g2R*(muu2**2 + muu3**2) + g2L*muu4*muuq)*muuqc + g2L*g2Lc*muuq*muuqc**2)
     -           + g1L*g1Lc*muu3*(g2L*muuq*(g2Rc*muu4 + g2Lc*muuqc) + 
     -            g2R*muuqc*(g2Lc*muu4 + g2Rc*muuqc)))*sqrtS**6))/(Denlegk2*DenU1c*DenU2*DenU2c)

!      


        enddo ! End of Neutralino neuc loop
	end
