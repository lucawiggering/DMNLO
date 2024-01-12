************************************************

	function stst2QQ_M2Boxcancel3()
	implicit none

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
#include "stst2QQ_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	  integer i,j,k,l,n,m,b,d,h,neu,neuc
	  complex*16 g0,g0L,g0R,g1,g1L,g1R,g2,g2L,g2R,g3,g3L,g3R
      complex*16 g1Lc, g2Lc, g3Lc, g4L, g4Lc
	  complex*16 g1Rc, g2Rc, g3Rc, g4R, g4Rc
	  double precision Ml0,Ml1,Ml3
	  complex*16 vbA, vbBL, vbBR, vbCLp1, vbCRp1, vbCLp2, vbCRp2, vbDp1, vbDk1, vbDk2
	  double precision Krondelta, ef, I3f, hf
	  complex*16 test,x12,x23, x13
	  complex*16 cspen,LLi3
      double precision d0modt,d0modu,mtreett,mtreetu,mtreeuu
      double precision mtreeinttt(6),mtreeintuu(6),mtreeinttu(6),mtreeintut(6)
      double precision mtreeneutt(6,6),mtreeneuuu(6,6),mtreeneutu(6,6)
	
	double complex KappaFun,testx13,ieps
	
	double precision testtreebox, boxfactor,boxfactor1
* functions & subroutines
	double precision colorf, colori,symfac
	double precision colboxglt,colboxglu,colboxneut,colboxneuu,colampglt,colampglu,colampneut,colampneuu 
	double precision colbox,colampgltu,colampglint,colampglneutu,colampglneuint,colampneutu,colampneuint
	double complex stst2QQ_M2Boxcancel3
	external myCget, myBget, myDget

	
	
	stst2QQ_M2Boxcancel3 = 0d0
	test = 0d0
	ieps = dcmplx(0d0,1d-17)

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

	  
! !!!!!Remove negative crosssection in resummation area
! !Define Treelevel
! 
! ! ! ---------------------t-channel Gluino squared------------------------------
! 	
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(g1L)
	  g1Rc = dconjg(g1R)

	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(g2L)
	  g2Rc = dconjg(g2R)
	

	  mtreett = (-((sqrtS**4*(g1R*
     -            (g1Rc*g2R*muugl*(g2Rc*(-1 + muu3**2 + muu4**2)* 
     -              muugl + g2Lc*muu4*(-muu1**2 + muu3**2 + tred)) + 
     -           g1Lc*g2R*muu3*muugl*
     -            (2*g2Lc*muu4*muugl + 
     -              g2Rc*(-muu2**2 + muu4**2 + tred)) + 
     -           g1Lc*g2L*muu3*
     -            (2*g2Rc*muu4*tred + 
     -              g2Lc*muugl*(-muu2**2 + muu4**2 + tred)) + 
     -           g1Rc*g2L*(g2Rc*muu4*muugl*
     -               (-muu1**2 + muu3**2 + tred) + 
     -              g2Lc*(muu3**2*muu4**2 + 
     -                 muu1**2*(muu2**2 - muu4**2 - tred) + tred + 
     -                 tred**2 - muu2**2*(muu3**2 + tred)))) + 
     -        g1L*(g1Lc*g2L*muugl*
     -            (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -              g2Rc*muu4*(-muu1**2 + muu3**2 + tred)) + 
     -           g1Rc*g2L*muu3*muugl*
     -            (2*g2Rc*muu4*muugl + 
     -              g2Lc*(-muu2**2 + muu4**2 + tred)) + 
     -           g1Rc*g2R*muu3*
     -            (2*g2Lc*muu4*tred + 
     -              g2Rc*muugl*(-muu2**2 + muu4**2 + tred)) + 
     -           g1Lc*g2R*(g2Lc*muu4*muugl*
     -               (-muu1**2 + muu3**2 + tred) + 
     -              g2Rc*(muu3**2*muu4**2 + 
     -                 muu1**2*(muu2**2 - muu4**2 - tred) + tred + 
     -                 tred**2 - muu2**2*(muu3**2 + tred))))))/
     -    (DenGlT*dconjg(DenGlT))))
     

! !---------------------u-channel Gluino squared------------------------------
! ! ! 
 
	
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1Lc = dconjg(g1L)
	  g1Rc = dconjg(g1R)
	
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2Lc = dconjg(g2L)
	  g2Rc = dconjg(g2R)
	

	  mtreeuu = (-((sqrtS**4*(g1R*(g1Rc*g2R*muugl*
     -            (g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -              g2Lc*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)) + 
     -           g1Lc*g2R*muu3*muugl*
     -            (2*g2Lc*muu4*muugl + 
     -              g2Rc*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) +
     -             g1Lc*g2L*muu3*
     -            (2*g2Rc*muu4*
     -               (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                 tred) + 
     -              g2Lc*muugl*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -           g1Rc*g2L*(g2Rc*muu4*muugl*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -              g2Lc*(muu3**4 - muu4**2 + muu4**4 + 
     -                 muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -                 muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -                 muu2**2*tred - 2*muu4**2*tred + tred**2))) + 
     -        g1L*(g1Lc*g2L*muugl*
     -            (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -              g2Rc*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred)) + 
     -           g1Rc*g2Rc*muu3*muugl*
     -            (2*g2L*muu4*muugl + 
     -              g2R*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) +
     -             g1Rc*g2Lc*muu3*
     -            (2*g2R*muu4*
     -               (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                 tred) + 
     -              g2L*muugl*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -           g1Lc*g2R*(g2Lc*muu4*muugl*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -              g2Rc*(muu3**4 - muu4**2 + muu4**4 + 
     -                 muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -                 muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -                 muu2**2*tred - 2*muu4**2*tred + tred**2)))))/
     -    (DenGlU*dconjg(DenGlU))))
! 	!      
! ---------------------t-channel and u-channel Gluino squared------------------------------	
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(GluLt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Rc = dconjg(GluRt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(GluLt1(isf1,itt1,ftt2,igen1,fgen2))
	
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Rc = dconjg(GluRt1(isf1,itt1,ftt2,igen1,fgen2))

	   mtreetu=((-sqrtS**4*(g1R*(g1Rc*g2R*muugl*
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
     -      g1L*(g1Lc*g2L*muugl*
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
     -               muu3**2*(muu4**2 + tred))))))/
     -  (DenGlT*dconjg(DenGlU)))
     
     
     
     
      do n=1,6
      do k=1,6
     

! ! !---------------------t-channel Neutralino squared------------------------------
! ! !         

	  g3L = ALt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3R = ARt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))
	
	  g4L = ALt1(k,isf2,itt2,ftt2,igen2,fgen2)
      g4R = ARt1(k,isf2,itt2,ftt2,igen2,fgen2)
	  g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  	

	  Mtreeneutt(k,n) = (-((sqrtS**4*
     -        (g3R*g3Rc*g4L*g4Lc*muu1**2*muu2**2 + 
     -        g3L*g3Lc*g4R*g4Rc*muu1**2*muu2**2 - 
     -        g3R*g3Rc*g4L*g4Lc*muu2**2*muu3**2 - 
     -        g3L*g3Lc*g4R*g4Rc*muu2**2*muu3**2 - 
     -        g3R*g3Rc*g4L*g4Lc*muu1**2*muu4**2 - 
     -        g3L*g3Lc*g4R*g4Rc*muu1**2*muu4**2 + 
     -        g3R*g3Rc*g4L*g4Lc*muu3**2*muu4**2 + 
     -        g3L*g3Lc*g4R*g4Rc*muu3**2*muu4**2 + 
     -        g3R*g3Rc*g4L*g4Lc*tred + g3L*g3Lc*g4R*g4Rc*tred - 
     -        g3R*g3Rc*g4L*g4Lc*muu1**2*tred - 
     -        g3L*g3Lc*g4R*g4Rc*muu1**2*tred - 
     -        g3R*g3Rc*g4L*g4Lc*muu2**2*tred - 
     -        g3L*g3Lc*g4R*g4Rc*muu2**2*tred + 
     -        2*g3L*g3Rc*g4Lc*g4R*muu3*muu4*tred + 
     -        2*g3Lc*g3R*g4L*g4Rc*muu3*muu4*tred + 
     -        g3R*g3Rc*g4L*g4Lc*tred**2 + g3L*g3Lc*g4R*g4Rc*tred**2 + 
     -        (g3L*g4L*(g3Lc*g4Rc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -              g3Rc*g4Lc*muu3*(-muu2**2 + muu4**2 + tred)) + 
     -           g3R*g4R*(g3Rc*g4Lc*muu4*
     -               (-muu1**2 + muu3**2 + tred) + 
     -              g3Lc*g4Rc*muu3*(-muu2**2 + muu4**2 + tred)))*
     -         muuqt(k) + muupt(n)*
     -         (g3Lc*g4Lc*(g3L*g4R*muu4*(-muu1**2 + muu3**2 + tred) + 
     -              g3R*g4L*muu3*(-muu2**2 + muu4**2 + tred)) + 
     -           g3Rc*g4Rc*(g3R*g4L*muu4*
     -               (-muu1**2 + muu3**2 + tred) + 
     -              g3L*g4R*muu3*(-muu2**2 + muu4**2 + tred)) + 
     -           (g3L*g4L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -                 g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)) + 
     -              g3R*g4R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -                 g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)))*muuqt(k)))
     -      )/(dconjg(DenGauT(n))*DenGauT(k))))
!      
!       print*,'tt',stst2qq_M2Tree-test
! 	  test = stst2qq_M2Tree

     
!---------------------u-channel Neutralino squared------------------------------

	  
	
	  g3L = ALt1(k,isf2,itt2,ftt1,igen2,fgen1)
	  g3R = ARt1(k,isf2,itt2,ftt1,igen2,fgen1)
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4L = ALt1(k,isf1,itt1,ftt2,igen1,fgen2)
	  g4R = ARt1(k,isf1,itt1,ftt2,igen1,fgen2)
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	

	  Mtreeneuuu(k,n) =(-((sqrtS**4*(g3R*g3Rc*
     -        g4L*g4Lc*muu1**2*muu2**2 + 
     -        g3L*g3Lc*g4R*g4Rc*muu1**2*muu2**2 - 
     -        g3R*g3Rc*g4L*g4Lc*muu3**2 - g3L*g3Lc*g4R*g4Rc*muu3**2 + 
     -        g3R*g3Rc*g4L*g4Lc*muu2**2*muu3**2 + 
     -        g3L*g3Lc*g4R*g4Rc*muu2**2*muu3**2 + 
     -        g3R*g3Rc*g4L*g4Lc*muu3**4 + g3L*g3Lc*g4R*g4Rc*muu3**4 - 
     -        2*g3L*g3Rc*g4Lc*g4R*muu3*muu4 - 
     -        2*g3Lc*g3R*g4L*g4Rc*muu3*muu4 + 
     -        2*g3L*g3Rc*g4Lc*g4R*muu1**2*muu3*muu4 + 
     -        2*g3Lc*g3R*g4L*g4Rc*muu1**2*muu3*muu4 + 
     -        2*g3L*g3Rc*g4Lc*g4R*muu2**2*muu3*muu4 + 
     -        2*g3Lc*g3R*g4L*g4Rc*muu2**2*muu3*muu4 + 
     -        2*g3L*g3Rc*g4Lc*g4R*muu3**3*muu4 + 
     -        2*g3Lc*g3R*g4L*g4Rc*muu3**3*muu4 - 
     -        g3R*g3Rc*g4L*g4Lc*muu4**2 - g3L*g3Lc*g4R*g4Rc*muu4**2 + 
     -        g3R*g3Rc*g4L*g4Lc*muu1**2*muu4**2 + 
     -        g3L*g3Lc*g4R*g4Rc*muu1**2*muu4**2 + 
     -        3*g3R*g3Rc*g4L*g4Lc*muu3**2*muu4**2 + 
     -        3*g3L*g3Lc*g4R*g4Rc*muu3**2*muu4**2 + 
     -        2*g3L*g3Rc*g4Lc*g4R*muu3*muu4**3 + 
     -        2*g3Lc*g3R*g4L*g4Rc*muu3*muu4**3 + 
     -        g3R*g3Rc*g4L*g4Lc*muu4**4 + g3L*g3Lc*g4R*g4Rc*muu4**4 + 
     -        g3R*g3Rc*g4L*g4Lc*tred + g3L*g3Lc*g4R*g4Rc*tred - 
     -        g3R*g3Rc*g4L*g4Lc*muu1**2*tred - 
     -        g3L*g3Lc*g4R*g4Rc*muu1**2*tred - 
     -        g3R*g3Rc*g4L*g4Lc*muu2**2*tred - 
     -        g3L*g3Lc*g4R*g4Rc*muu2**2*tred - 
     -        2*g3R*g3Rc*g4L*g4Lc*muu3**2*tred - 
     -        2*g3L*g3Lc*g4R*g4Rc*muu3**2*tred - 
     -        2*g3L*g3Rc*g4Lc*g4R*muu3*muu4*tred - 
     -        2*g3Lc*g3R*g4L*g4Rc*muu3*muu4*tred - 
     -        2*g3R*g3Rc*g4L*g4Lc*muu4**2*tred - 
     -        2*g3L*g3Lc*g4R*g4Rc*muu4**2*tred + 
     -        g3R*g3Rc*g4L*g4Lc*tred**2 + g3L*g3Lc*g4R*g4Rc*tred**2 + 
     -        (g3L*g4L*(g3Lc*g4Rc*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -              g3Rc*g4Lc*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -           g3R*g4R*(g3Rc*g4Lc*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -              g3Lc*g4Rc*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)))*
     -         muuqu(k) + muupu(n)*
     -         (g3Lc*g4Lc*(g3L*g4R*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -              g3R*g4L*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -           g3Rc*g4Rc*(g3R*g4L*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -              g3L*g4R*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -           (g3L*g4L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -                 g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)) + 
     -              g3R*g4R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -                 g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)))*muuqu(k)))
     -      )/(dconjg(DenGauU(n))*DenGauU(k))))

!---------------------t-channel and u-channel Neutralino squared------------------------------

	  g3L = ALt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g3R = ARt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4L = ALt1(k,isf2,itt2,ftt2,igen2,fgen2)
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
	  g4R = ARt1(k,isf2,itt2,ftt2,igen2,fgen2)
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))

	  Mtreeneutu(k,n) =((-sqrtS**4*
     -      (-(g3R*g3Rc*g4L*g4Lc*muu1**2*muu2**2) - 
     -      g3L*g3Lc*g4R*g4Rc*muu1**2*muu2**2 - 
     -      g3L*g3Rc*g4Lc*g4R*muu1**2*muu3*muu4 - 
     -      g3Lc*g3R*g4L*g4Rc*muu1**2*muu3*muu4 - 
     -      g3L*g3Rc*g4Lc*g4R*muu2**2*muu3*muu4 - 
     -      g3Lc*g3R*g4L*g4Rc*muu2**2*muu3*muu4 + 
     -      g3L*g3Rc*g4Lc*g4R*muu3**3*muu4 + 
     -      g3Lc*g3R*g4L*g4Rc*muu3**3*muu4 + 
     -      g3R*g3Rc*g4L*g4Lc*muu3**2*muu4**2 + 
     -      g3L*g3Lc*g4R*g4Rc*muu3**2*muu4**2 + 
     -      g3L*g3Rc*g4Lc*g4R*muu3*muu4**3 + 
     -      g3Lc*g3R*g4L*g4Rc*muu3*muu4**3 - g3R*g3Rc*g4L*g4Lc*tred - 
     -      g3L*g3Lc*g4R*g4Rc*tred + g3R*g3Rc*g4L*g4Lc*muu1**2*tred + 
     -      g3L*g3Lc*g4R*g4Rc*muu1**2*tred + 
     -      g3R*g3Rc*g4L*g4Lc*muu2**2*tred + 
     -      g3L*g3Lc*g4R*g4Rc*muu2**2*tred + 
     -      g3R*g3Rc*g4L*g4Lc*muu3**2*tred + 
     -      g3L*g3Lc*g4R*g4Rc*muu3**2*tred + 
     -      g3R*g3Rc*g4L*g4Lc*muu4**2*tred + 
     -      g3L*g3Lc*g4R*g4Rc*muu4**2*tred - 
     -      g3R*g3Rc*g4L*g4Lc*tred**2 - g3L*g3Lc*g4R*g4Rc*tred**2 + 
     -      (g3L*g4L*(g3Lc*g4Rc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g3Rc*g4Lc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -         g3R*g4R*(g3Rc*g4Lc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g3Lc*g4Rc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)))*
     -       muuqt(k) + muupu(n)*
     -       (g3Lc*g4Lc*(g3L*g4R*muu4*(-muu1**2 + muu3**2 + tred) + 
     -            g3R*g4L*muu3*(-muu2**2 + muu4**2 + tred)) + 
     -         g3Rc*g4Rc*(g3R*g4L*muu4*(-muu1**2 + muu3**2 + tred) + 
     -            g3L*g4R*muu3*(-muu2**2 + muu4**2 + tred)) + 
     -         (g3L*g4L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -               g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)) + 
     -            g3R*g4R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -               g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)))*muuqt(k))))/
     -  (dconjg(DenGauU(n))*DenGauT(k)))
     
        enddo
        
        
        
        
! ! !---------------------t-channel Gluino t-channel Neutralino--------------------------------------------	  
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  
	  g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))

	  mtreeinttt(n) = (-sqrtS**4*(g1R*(2*g2L*g3Lc*g4Rc*muu3*muu4*tred + 
     -         g2R*g3Rc*g4Lc*muu4*muugl*(-muu1**2 + muu3**2 + tred) + 
     -         g2R*g3Lc*g4Rc*muu3*muugl*(-muu2**2 + muu4**2 + tred) + 
     -         g2L*g3Rc*g4Lc*
     -          (muu3**2*muu4**2 + 
     -            muu1**2*(muu2**2 - muu4**2 - tred) + tred + 
     -            tred**2 - muu2**2*(muu3**2 + tred))) + 
     -      g1L*(2*g2R*g3Rc*g4Lc*muu3*muu4*tred + 
     -         g2L*g3Lc*g4Rc*muu4*muugl*(-muu1**2 + muu3**2 + tred) + 
     -         g2L*g3Rc*g4Lc*muu3*muugl*(-muu2**2 + muu4**2 + tred) + 
     -         g2R*g3Lc*g4Rc*
     -          (muu3**2*muu4**2 + 
     -            muu1**2*(muu2**2 - muu4**2 - tred) + tred + 
     -            tred**2 - muu2**2*(muu3**2 + tred))) + 
     -      (g1R*(2*g2R*g3Lc*g4Lc*muu3*muu4*muugl + 
     -            g2R*g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2L*g3Rc*g4Rc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -            g2L*g3Lc*g4Lc*muu3*(-muu2**2 + muu4**2 + tred)) + 
     -         g1L*(2*g2L*g3Rc*g4Rc*muu3*muu4*muugl + 
     -            g2L*g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2R*g3Lc*g4Lc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -            g2R*g3Rc*g4Rc*muu3*(-muu2**2 + muu4**2 + tred)))*
     -       muupt(n)))/(DenGlT*dconjg(DenGauT(n)))        
     
    
! 	  
!---------------------u-channel Gluino u-channel Neutralino--------------------------------------------------	  
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  
	  
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
      mtreeintuu(n) = (-sqrtS**4*(g1R*(2*g2L*g3Lc*g4Rc*muu3*muu4*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred) + 
     -         g2R*g3Rc*g4Lc*muu4*muugl*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -         g2R*g3Lc*g4Rc*muu3*muugl*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred) + 
     -         g2L*g3Rc*g4Lc*
     -          (muu3**4 - muu4**2 + muu4**4 + 
     -            muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -            muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -            muu2**2*tred - 2*muu4**2*tred + tred**2)) + 
     -      g1L*(2*g2R*g3Rc*g4Lc*muu3*muu4*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tred) + 
     -         g2L*g3Lc*g4Rc*muu4*muugl*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -         g2L*g3Rc*g4Lc*muu3*muugl*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred) + 
     -         g2R*g3Lc*g4Rc*
     -          (muu3**4 - muu4**2 + muu4**4 + 
     -            muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tred) + 
     -            muu1**2*(muu2**2 + muu4**2 - tred) + tred - 
     -            muu2**2*tred - 2*muu4**2*tred + tred**2)) + 
     -      (g1R*(2*g2R*g3Lc*g4Lc*muu3*muu4*muugl + 
     -            g2R*g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2L*g3Rc*g4Rc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g2L*g3Lc*g4Lc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -         g1L*(2*g2L*g3Rc*g4Rc*muu3*muu4*muugl + 
     -            g2L*g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2R*g3Lc*g4Lc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g2R*g3Rc*g4Rc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)))*
     -       muupu(n)))/(DenGlU*dconjg(DenGauU(n)))        
        
!         
!---------------------t-channel Gluino u-channel Neutralino--------------------------------------------------------	  
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
! 	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
	
      mtreeinttu(n) = (-sqrtS**4*(g1R*(g2L*g3Lc*g4Rc*muu3*muu4*
     -          (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -         g2R*g3Rc*g4Lc*muu4*muugl*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -         g2R*g3Lc*g4Rc*muu3*muugl*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred) + 
     -         g2L*g3Rc*g4Lc*
     -          ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -            muu1**2*(-muu2**2 + tred) + muu3**2*(muu4**2 + tred)
     -            )) + g1L*(g2R*g3Rc*g4Lc*muu3*muu4*
     -          (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -         g2L*g3Lc*g4Rc*muu4*muugl*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -         g2L*g3Rc*g4Lc*muu3*muugl*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred) + 
     -         g2R*g3Lc*g4Rc*
     -          ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -            muu1**2*(-muu2**2 + tred) + muu3**2*(muu4**2 + tred)
     -            )) + (g1R*(2*g2R*g3Lc*g4Lc*muu3*muu4*muugl + 
     -            g2R*g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2L*g3Rc*g4Rc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -            g2L*g3Lc*g4Lc*muu3*(-muu2**2 + muu4**2 + tred)) + 
     -         g1L*(2*g2L*g3Rc*g4Rc*muu3*muu4*muugl + 
     -            g2L*g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2R*g3Lc*g4Lc*muu4*(-muu1**2 + muu3**2 + tred) + 
     -            g2R*g3Rc*g4Rc*muu3*(-muu2**2 + muu4**2 + tred)))*
     -       muupu(n)))/(DenGlT*dconjg(DenGauU(n)))   
     
    
!         
!---------------------u-channel Gluino t-channel Neutralino---------------------------------------	  
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)

	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)

	  g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))

	  g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  
	  mtreeintut(n) = (-sqrtS**4*(g1R*(g2L*g3Lc*g4Rc*muu3*muu4*
     -          (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -         g2R*g3Rc*g4Lc*muu4*muugl*(-muu1**2 + muu3**2 + tred) + 
     -         g2R*g3Lc*g4Rc*muu3*muugl*(-muu2**2 + muu4**2 + tred) + 
     -         g2L*g3Rc*g4Lc*
     -          ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -            muu1**2*(-muu2**2 + tred) + muu3**2*(muu4**2 + tred)
     -            )) + g1L*(g2R*g3Rc*g4Lc*muu3*muu4*
     -          (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -         g2L*g3Lc*g4Rc*muu4*muugl*(-muu1**2 + muu3**2 + tred) + 
     -         g2L*g3Rc*g4Lc*muu3*muugl*(-muu2**2 + muu4**2 + tred) + 
     -         g2R*g3Lc*g4Rc*
     -          ((-1 + muu2**2 + muu4**2 - tred)*tred + 
     -            muu1**2*(-muu2**2 + tred) + muu3**2*(muu4**2 + tred)
     -            )) + (g1R*(2*g2R*g3Lc*g4Lc*muu3*muu4*muugl + 
     -            g2R*g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2L*g3Rc*g4Rc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g2L*g3Lc*g4Lc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)) + 
     -         g1L*(2*g2L*g3Rc*g4Rc*muu3*muu4*muugl + 
     -            g2L*g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2R*g3Lc*g4Lc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tred) + 
     -            g2R*g3Rc*g4Rc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred)))*
     -       muupt(n)))/(DenGlU*dconjg(DenGauT(n)))
     
     

! ! 
        enddo
 

! !      Defining color for Box and amplitude

        colboxglt = CF/2d0*(Nc+1d0)
        colboxglu = -CF/2d0*(Nc+1d0)
        colboxneut = -CF*Nc
        colboxneuu = CF*Nc
        
        colampglt = -(Nc**2-1d0)/(2d0*(2d0*Nc*(Nc-1d0)))
        colampglu = (Nc**2-1d0)/(2d0*(2d0*Nc*(Nc-1d0)))
        colampneut = Nc*(Nc-1d0)/(2d0*(Nc*(Nc-1d0)))
        colampneuu = -Nc*(Nc-1d0)/(2d0*(Nc*(Nc-1d0)))
        
! c---------------------Calculation of triplet part of Tree-level-------------------------------------------
! 
!         stst2qq_M2Tree3 = (Nc**2-1d0)**2/(8*Nc*(Nc-1d0))*(GltGlt + GluGlu - GltGlu)+
!      -      Nc*(Nc-1d0)/(2d0)*(NeutNeut + NeuuNeuu - NeutNeuu) -
!      -      (Nc**2-1d0)/(4d0)*(GltNeut + GluNeuu - GltNeuu - GluNeut)
            
!             colbox = (-(Nc+1d0)/(2*Nc))
!             
!             colampgltu = (Nc**2-1d0)**2/(8*Nc*(Nc-1d0))
!             colampglint = -colampgltu 
!             colampneutu = Nc*(Nc-1d0)/(2d0)
!             colampneuint = -colampneutu
!             colampglneutu = -(Nc**2-1d0)/(4d0)
!             colampglneuint = -colampglneutu


c-------------------------------------------------------------------------------------------------------------         
c---------------------------------FOR GLUINO BOX1 CALCULATION-------------------------------------------------

      x13 = (-(1-vrel)/cdsqrt(dcmplx(1d0-vrel**2)))
      x12 = (cdsqrt(dcmplx(1d0-4*muu1*muugl/(muu3**2-(muu1-muugl)**2)))-1d0)/(cdsqrt(dcmplx(1d0-4*muu1*muugl/(muu3**2-(muu1-muugl)**2)))+1d0)		
      x23 = (cdsqrt(dcmplx(1d0-4*muu2*muugl/(muu4**2-(muu2-muugl)**2)))-1d0)/(cdsqrt(dcmplx(1d0-4*muu2*muugl/(muu4**2-(muu2-muugl)**2)))+1d0)  

      
!       d0modt = 1d0/(sqrtS**4*muu1*muu2*(tred-muugl**2))*(-dsqrt(1-vrel**2)/(2*vrel))*(cdlog(dcmplx(x12))**2 +
!      -       cdlog(dcmplx(x23))**2+cspen(dcmplx(1d0-x13/(x12*x23))) +cspen(dcmplx(1d0-x13*x12*x23)) +cspen(dcmplx(1d0-x13*x23/x12))+cspen(dcmplx(1d0-x13*x12/x23)) )

      d0modt = 1d0/(sqrtS**4*muu1*muu2*(tred-muugl**2))*(-1d0/(2*vrel))*(cdlog(dcmplx(x12))**2 +
     -       cdlog(dcmplx(x23))**2+LLi3(x13,x12,x23)+LLi3(x13,1d0/x12,x23)+LLi3(x13,x12,1d0/x23)+LLi3(x13,1d0/x12,1d0/x23))
     
!             print*,'x12',x12,'x13',x13,'x23',x23
!             print*,'d0modt gluino []',(cdlog(dcmplx(x12))**2 +
!      -       cdlog(dcmplx(x23))**2+LLi3(x13,x12,x23)+LLi3(x13,1d0/x12,x23)+LLi3(x13,x12,1d0/x23)+LLi3(x13,1d0/x12,1d0/x23)),Pi**2
!      
! !             print*,'check 2*Li2(2)',LLi3(x13,1d0/x12,x23)+LLi3(x13,x12,1d0/x23),Pi**2/2d0,2d0*cspen(dcmplx(1d0-x13*x23/x12))
!      
!             print*,'check A.29',LLi3(-1d0,1d0/x12,1d0/x12),-LLi3(-1d0,x12,x12)-1d0/2d0*cdlog(-x12**2)!+Pi**2/2d0+cI*Pi*cdlog(x12**2/(x12**2+1d0)**2),cspen(dcmplx((1d0+x12**2)/x12**2))
! !             print*,'check A.27',LLi3(-1d0,x12,x12),

! 
!       d0modt = 1d0/(sqrtS**4*muu1*muu2*(tred-muugl**2))*(-1d0/(2*vrel))*(cdlog(dcmplx(x12))**2 +
!      -       cdlog(dcmplx(x23))**2+LLi3(x13,x12,x23)+LLi3(x13,1d0/x12,x23)+LLi3(x13,x12,1d0/x23)+LLi3(x13,1d0/x12,1d0/x23))

      x13 = (-(1-vrel)/cdsqrt(dcmplx(1d0-vrel**2)))
      x12 = (cdsqrt(dcmplx(1d0-4*muu2*muugl/(muu3**2-(muu2-muugl)**2)))-1d0)/(cdsqrt(dcmplx(1d0-4*muu2*muugl/(muu3**2-(muu2-muugl)**2)))+1d0)		
      x23 = (cdsqrt(dcmplx(1d0-4*muu1*muugl/(muu4**2-(muu1-muugl)**2)))-1d0)/(cdsqrt(dcmplx(1d0-4*muu1*muugl/(muu4**2-(muu1-muugl)**2)))+1d0)		


!       d0modu = 1d0/(sqrtS**4*muu1*muu2*(ured-muugl**2))*(-cdsqrt(dcmplx(1d0-vrel**2))/(2*vrel))*(cdlog(dcmplx(x12))**2 +
!      -       cdlog(dcmplx(x23))**2+LLi3(x13,x12,x23)+LLi3(x13,1d0/x12,x23)+LLi3(x13,x12,1d0/x23)+LLi3(x13,1d0/x12,1d0/x23))
      d0modu = 1d0/(sqrtS**4*muu1*muu2*(ured-muugl**2))*(-1d0/(2*vrel))*(cdlog(dcmplx(x12))**2 +
     -       cdlog(dcmplx(x23))**2+LLi3(x13,x12,x23)+LLi3(x13,1d0/x12,x23)+LLi3(x13,x12,1d0/x23)+LLi3(x13,1d0/x12,1d0/x23))


! c----------------------Boxcorrection tbox-t-------------------------------
	  colorf= colboxglt*colampglt
!       colorf = 
	  
	  boxfactor = GS**2*sqrtS**2*2d0*(1d0-muu1**2-muu2**2)*d0modt/(4d0*Pi)**2*DenGlT
	  
! 	  	  print*,'boxfac gluino',2d0*boxfactor,-GS**2/(4d0*Pi)*Pi/vrel,2d0*boxfactor*(-(Nc+1d0)/(2*Nc))
	  	  
! 	  	  print*,'compare factors box',vrel,GS**2
! 
! 	  print*,'boxfac',boxfactor*DenGlT
	  

	  stst2QQ_M2Boxcancel3 = stst2QQ_M2Boxcancel3 + 2d0*colorf*boxfactor*mtreett
! 	  	  print*,'Abzug tt',- 2d0*colorf*boxfactor*mtreett*DenGlT

c----------------------Boxcorrection ubox-u-------------------------------
	  colorf= colboxglu*colampglu
	  
	  boxfactor = GS**2*sqrtS**2*2d0*(1d0-muu1**2-muu2**2)*d0modu/(4d0*Pi)**2*DenGlU
	  
! 	  	  	  print*,'boxfac gluino',2d0*boxfactor,-GS**2/(4d0*Pi)*Pi/vrel,2d0*boxfactor*(-(Nc+1d0)/(2*Nc))

	        
	  stst2QQ_M2Boxcancel3 = stst2QQ_M2Boxcancel3 + 2d0*colorf*boxfactor*mtreeuu
! 	  print*,'Abzug uu',- 2d0*colorf*boxfactor*mtreeuu*DenGlU
	  
c----------------------Boxcorrection tbox-u-------------------------------
	  colorf= colboxglt*colampglu
	    
	  boxfactor = GS**2*sqrtS**2*2d0*(1d0-muu1**2-muu2**2)*d0modt/(4d0*Pi)**2*DenGlT
	

	  stst2QQ_M2Boxcancel3 = stst2QQ_M2Boxcancel3 + 2d0*colorf*boxfactor*mtreetu
! 	  print*,'Abzug tu',-2d0*colorf*boxfactor*mtreetu*DenGlT
	  
c----------------------Boxcorrection ubox-t-------------------------------
	  colorf= colboxglu*colampglt
	    
	  boxfactor = GS**2*sqrtS**2*2d0*(1d0-muu1**2-muu2**2)*d0modu/(4d0*Pi)**2*DenGlu
	

	  stst2QQ_M2Boxcancel3 = stst2QQ_M2Boxcancel3 + 2d0*colorf*boxfactor*mtreetu
! 	  	  print*,'Abzug ut',- 2d0*colorf*boxfactor*mtreetu*DenGlu
	  
	  
!         print*,'Nach gluino mit gluino',stst2QQ_M2Boxcancel3
        
! 	  
	  do neu = 1,6
!         
c----------------------Boxcorrection tbox-tNeu-------------------------------
	  colorf= colboxglt*colampneut
	  
	  boxfactor = GS**2*sqrtS**2*2d0*(1d0-muu1**2-muu2**2)*d0modt/(4d0*Pi)**2*DenGlT
	

	  stst2QQ_M2Boxcancel3 = stst2QQ_M2Boxcancel3 + 2d0*colorf*boxfactor*mtreeinttt(neu)
	  
! 	  print*,'Boxcorrectiontt',- 2d0*colorf*boxfactor*mtreeinttt(neu)*DenGlT
	  
c----------------------Boxcorrection ubox-u-------------------------------
	  colorf= colboxglu*colampneuu
	  
	  boxfactor = GS**2*sqrtS**2*2d0*(1d0-muu1**2-muu2**2)*d0modu/(4d0*Pi)**2*DenGlU
	

	  stst2QQ_M2Boxcancel3 = stst2QQ_M2Boxcancel3 + 2d0*colorf*boxfactor*mtreeintuu(neu)
! 	  
! c----------------------Boxcorrection tbox-u-------------------------------
	  colorf= colboxglt*colampneuu
! 	  print*,'color box tu gl',colorf
	    
	  boxfactor = GS**2*sqrtS**2*2d0*(1d0-muu1**2-muu2**2)*d0modt/(4d0*Pi)**2*DenGlT
	

	  stst2QQ_M2Boxcancel3 = stst2QQ_M2Boxcancel3 + 2d0*colorf*boxfactor*mtreeinttu(neu)
	  
	  
c----------------------Boxcorrection ubox-t-------------------------------
	  colorf= colboxglu*colampneut
! 	  print*,'color box ut gl',colorf
	  boxfactor = GS**2*sqrtS**2*2d0*(1d0-muu1**2-muu2**2)*d0modu/(4d0*Pi)**2*DenGlu
	

	  stst2QQ_M2Boxcancel3 = stst2QQ_M2Boxcancel3 + 2d0*colorf*boxfactor*mtreeintut(neu)

! 	  	  print*,'boxfac',2d0*boxfactor*colboxglu,-GS**2/(4d0*Pi)*Pi/vrel*colboxglu
! 	  
! ! ! 	  print*,'Boxcorrectiontt',- 2d0*colorf*boxfactor*mtreeintut(neu)*DenGlu
! 	  
	  enddo
! ! ! 	  	  
! ! 	  print*,'Box',stst2QQ_M2Boxcancel3	  
! 
! 
!         print*,'Am nach Gluino korrektur ganz',stst2QQ_M2Boxcancel3
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 


! 	  	  
! ! ! ! !         FOR NEUTRALINO BOX1 CALCULATION
! 
      do neu = 1,6
      x13 = (-(1-vrel)/cdsqrt(dcmplx(1d0-vrel**2)))
      
      x12 = (cdsqrt(dcmplx(1d0-4*muu1*muuqt(neu)/(muu3**2-(muu1-muuqt(neu))**2)))-1d0)/
     -      (cdsqrt(dcmplx(1d0-4*muu1*muuqt(neu)/(muu3**2-(muu1-muuqt(neu))**2)))+1d0)		
     
     
      x23 = (cdsqrt(dcmplx(1d0-4*muu2*muuqt(neu)/(muu4**2-(muu2-muuqt(neu))**2)))-1d0)/
     -      (cdsqrt(dcmplx(1d0-4*muu2*muuqt(neu)/(muu4**2-(muu2-muuqt(neu))**2)))+1d0)  
!       

!                     print*,'x12',x12,'x13',x13,'x23',x23
!       d0modt = 1d0/(sqrtS**4*muu1*muu2*(tred-muuqt(neu)**2))*(-dsqrt(1-vrel**2)/(2*vrel))*(cdlog(dcmplx(x12))**2 +
!      -       cdlog(dcmplx(x23))**2+cspen(dcmplx(1d0-x13/(x12*x23))) +cspen(dcmplx(1d0-x13*x12*x23)) +cspen(dcmplx(1d0-x13*x23/x12))+cspen(dcmplx(1d0-x13*x12/x23)) )

!       d0modt = 1d0/(sqrtS**4*muu1*muu2*(tred-muuqt(neu)**2))*(-1/(2*vrel))*(cdlog(dcmplx(x12))**2 +
!      -       cdlog(dcmplx(x23))**2+cspen(dcmplx(1d0-x13/(x12*x23))) +cspen(dcmplx(1d0-x13*x12*x23)) +cspen(dcmplx(1d0-x13*x23/x12))+cspen(dcmplx(1d0-x13*x12/x23)) )
     
!       d0modt = 1d0/(sqrtS**4*muu1*muu2*(tred-muuqt(neu)**2))*(-cdsqrt(dcmplx(1d0-vrel**2))/(2*vrel))*(cdlog(dcmplx(x12))**2 +
!      -       cdlog(dcmplx(x23))**2+LLi3(x13,x12,x23)+LLi3(x13,1d0/x12,x23)+LLi3(x13,x12,1d0/x23)+LLi3(x13,1d0/x12,1d0/x23))

      d0modt = 1d0/(sqrtS**4*muu1*muu2*(tred-muuqt(neu)**2))*(-1d0/(2*vrel))*(cdlog(dcmplx(x12))**2 +
     -       cdlog(dcmplx(x23))**2+LLi3(x13,x12,x23)+LLi3(x13,1d0/x12,x23)+LLi3(x13,x12,1d0/x23)+LLi3(x13,1d0/x12,1d0/x23))

!         print*,'pi^2 factor',(cdlog(dcmplx(x12))**2 +
!      -       cdlog(dcmplx(x23))**2+LLi3(x13,x12,x23)+LLi3(x13,1d0/x12,x23)+LLi3(x13,x12,1d0/x23)+LLi3(x13,1d0/x12,1d0/x23)),Pi**2
        
!         print*,'x13',x13,'x12',x12,'x23',x23
!         print*,'LLi3(x13,x12,x23)',LLi3(x13,x12,x23),'LLi3(x13,1d0/x12,x23)',LLi3(x13,1d0/x12,x23)
!         print*,'LLi3(x13,x12,1d0/x23)',LLi3(x13,x12,1d0/x23)
!         print*,'LLi3(x13,1d0/x12,1d0/x23)',LLi3(x13,1d0/x12,1d0/x23),cspen(dcmplx(1d0-x13/(x12*x23)))
!         print*,'topmasse',Mf(3,3)
      
!       print*,'x12,x13,x23',x12,x13,x23
!       
!       print*,'test',neu,(cdsqrt(dcmplx(1d0-4*muu1*muuqt(neu)/(muu3**2-(muu1-muuqt(neu))**2)))-1d0),(cdsqrt(dcmplx(1d0-4*muu1*muuqt(neu)/(muu3**2-(muu1-muuqt(neu))**2)))+1d0)	
!       print*,'test1',neu,cdsqrt(dcmplx(1d0-4*muu1*muuqt(neu)/(muu3**2-(muu1-muuqt(neu))**2))),cdsqrt(dcmplx(1d0-4*muu1*muuqt(neu)/(muu3**2-(muu1-muuqt(neu))**2)))
!       print*,'test2',neu,-4*muu1*muuqt(neu),(muu3**2-(muu1-muuqt(neu))**2),dcmplx(1d0-4*muu1*muuqt(neu)/(muu3**2-(muu1-muuqt(neu))**2))
! !       print*,'MNeu3',muuqt(neu), MNeu(3)
!         print*,'Neutralinomassen',MNeu(1), MNeu(2),MSf(1,3,3),MGl, Mf(3,3)
!       print*,'d0mod neu',d0modt
!       print*,'x12*x13*x23',x12*x13*x23, LLi3(x13,x12,x23),cspen(dcmplx(1d0-x13*x12*x23))
!       print*,'LLi3 function',cspen(dcmplx(1d0-x13*x12*x23)),cspen(dcmplx(1d0-x13*x12/x23)),cspen(dcmplx(1d0-x13*x23/x12)) ,cspen(dcmplx(1d0-x13/(x12*x23)))
! ! 
      do neuc = 1,6
! ! !!!!!!!!!!Boxcorrection tbox-t-------------------------------
	  colorf= colboxneut*colampneut
	    
	  boxfactor = GS**2*sqrtS**2*2d0*(1d0-muu1**2-muu2**2)*d0modt/(4d0*Pi)**2
	  
! 	  	  	  print*,'boxfac neutralino',2d0*boxfactor,-GS**2/(4d0*Pi)*Pi/vrel,2d0*boxfactor*(-(Nc+1d0)/(2*Nc))
	
	  stst2QQ_M2Boxcancel3 = stst2QQ_M2Boxcancel3 + 2d0*colorf*boxfactor*mtreeneutt(neu,neuc)*DenGauT(neu)
	  
! 	  
! !!!!!!!!!!Boxcorrection tbox-u-------------------------------
	  colorf= colboxneut*colampneuu
	    
	  boxfactor = GS**2*sqrtS**2*2d0*(1d0-muu1**2-muu2**2)*d0modt/(4d0*Pi)**2
	
	  stst2QQ_M2Boxcancel3 = stst2QQ_M2Boxcancel3 + 2d0*colorf*boxfactor*mtreeneutu(neu,neuc)*DenGauT(neu)
!         print*,'boxcorr',neu,neuc,- 2d0*colorf*boxfactor*mtreeneutu(neu,neuc)*DenGauT(neu),boxfactor*DenGauT(neu),mtreeneutu(neu,neuc)
!     
      enddo
      
      

	  
!!!!!!!!!!Boxcorrection tbox-tGl-------------------------------
	  colorf= colboxneut*colampglt
	  
	  boxfactor = GS**2*sqrtS**2*2d0*(1d0-muu1**2-muu2**2)*d0modt/(4d0*Pi)**2
	

	  stst2QQ_M2Boxcancel3 = stst2QQ_M2Boxcancel3 + 2d0*colorf*boxfactor*mtreeinttt(neu)*DenGauT(neu)
	  
! ! 	  print*,'Boxcorrectiontt',- 2d0*colorf*boxfactor*mtreeinttt(neu)*DenGlT

! ! !!!!!!!!!!Boxcorrection tbox-u-------------------------------
	  colorf= colboxneut*colampglu
! 	  print*,'color box tu neu',colorf
	    
	  boxfactor = GS**2*sqrtS**2*2d0*(1d0-muu1**2-muu2**2)*d0modt/(4d0*Pi)**2
	
	  stst2QQ_M2Boxcancel3 = stst2QQ_M2Boxcancel3 + 2d0*colorf*boxfactor*mtreeintut(neu)*DenGauT(neu)

	  
	  
	  
	  
	  
	  x13 = (-(1-vrel)/cdsqrt(dcmplx(1d0-vrel**2)))
      x12 = (cdsqrt(dcmplx(1d0-4*muu2*muuqt(neu)/(muu3**2-(muu2-muuqt(neu))**2)))-1d0)/
     -      (cdsqrt(dcmplx(1d0-4*muu2*muuqt(neu)/(muu3**2-(muu2-muuqt(neu))**2)))+1d0)		
      x23 = (cdsqrt(dcmplx(1d0-4*muu1*muuqt(neu)/(muu4**2-(muu1-muuqt(neu))**2)))-1d0)/
     -      (cdsqrt(dcmplx(1d0-4*muu1*muuqt(neu)/(muu4**2-(muu1-muuqt(neu))**2)))+1d0)		


!       d0modu = 1d0/(sqrtS**4*muu1*muu2*(ured-muuqt(neu)**2))*(-cdsqrt(dcmplx(1d0-vrel**2))/(2*vrel))*(cdlog(dcmplx(x12))**2 +
!      -       cdlog(dcmplx(x23))**2+LLi3(x13,x12,x23)+LLi3(x13,1d0/x12,x23)+LLi3(x13,x12,1d0/x23)+LLi3(x13,1d0/x12,1d0/x23))

      d0modu = 1d0/(sqrtS**4*muu1*muu2*(ured-muuqt(neu)**2))*(-1d0/(2*vrel))*(cdlog(dcmplx(x12))**2 +
     -       cdlog(dcmplx(x23))**2+LLi3(x13,x12,x23)+LLi3(x13,1d0/x12,x23)+LLi3(x13,x12,1d0/x23)+LLi3(x13,1d0/x12,1d0/x23))


      do neuc = 1,6
!!!!!!!!!!Boxcorrection ubox-u-------------------------------
	  colorf= colboxneuu*colampneuu
	    
	  boxfactor = GS**2*sqrtS**2*2d0*(1d0-muu1**2-muu2**2)*d0modu/(4d0*Pi)**2
	
	  stst2QQ_M2Boxcancel3 = stst2QQ_M2Boxcancel3 + 2d0*colorf*boxfactor*mtreeneuuu(neu,neuc)*DenGauU(neu)
	  
! 	  
! !!!!!!!!!!Boxcorrection ubox-t-------------------------------
	  colorf= colboxneuu*colampneut
	    
	  boxfactor = GS**2*sqrtS**2*2d0*(1d0-muu1**2-muu2**2)*d0modu/(4d0*Pi)**2
	
	  stst2QQ_M2Boxcancel3 = stst2QQ_M2Boxcancel3 + 2d0*colorf*boxfactor*mtreeneutu(neuc,neu)*DenGauU(neu)
!     
      enddo
! ! 	  
! ! ! ! !!!!!!!!!!Boxcorrection ubox-uGl-------------------------------
	  colorf= colboxneuu*colampglu
	  
	  boxfactor = GS**2*sqrtS**2*2d0*(1d0-muu1**2-muu2**2)*d0modu/(4d0*Pi)**2
	

	  stst2QQ_M2Boxcancel3 = stst2QQ_M2Boxcancel3 + 2d0*colorf*boxfactor*mtreeintuu(neu)*DenGauU(neu)
	  
! 	  print*,'Boxcorrectiontt',- 2d0*colorf*boxfactor*mtreeinttt(neu)*DenGlT

! !!!!!!!!!!Boxcorrection ubox-t-------------------------------
	  colorf= colboxneuu*colampglt
! 	  print*,'color box ut neu',colorf
	    
	  boxfactor = GS**2*sqrtS**2*2d0*(1d0-muu1**2-muu2**2)*d0modu/(4d0*Pi)**2
	
	  stst2QQ_M2Boxcancel3 = stst2QQ_M2Boxcancel3 + 2d0*colorf*boxfactor*mtreeinttu(neu)*DenGauU(neu)
        
      enddo
        
!         print*,'Am ENDE',stst2QQ_M2Boxcancel3
      
	end
                                   
