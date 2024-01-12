************************************************

	  function stst2qq_M2Counterterm()
	  implicit none

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
#include "stst2QQ_GenCouplings.h"
#include "stst2QQ_GenCounterterms.h"


 	  double precision findValW
      integer i,j,k,l,m,n

	  double precision colfac
	  double complex g1L, g1Lc, g2L, g2Lc, g3L, g3Lc, g4L, g4Lc
	  double complex g1R, g1Rc, g2R, g2Rc, g3R, g3Rc, g4R, g4Rc
	  double complex test

	  complex*16 stst2qq_M2Counterterm

	
	  stst2qq_M2Counterterm = 0d0
	  test = 0d0
      
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
! !  
c-------------------------------------------------------------------------------------------------------

! ! ---------------------t-channel Gluino squared------------------------------
! 
	  colfac = Nc/2d0*CF
	
	  g1L = dGluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = dGluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(GluLt1(isf1,itt1,ftt1,igen1,fgen1))
	  g1Rc = dconjg(GluRt1(isf1,itt1,ftt1,igen1,fgen1))
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(g2L)
	  g2Rc = dconjg(g2R)
	

	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*dreal(-((sqrtS**4*(g1R*
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
!      
!      
! ! ! ! !---------------------u-channel Gluino squared------------------------------
! ! ! ! ! ! 
	  colfac = Nc/2d0*CF
    
	
	  g1L = dGluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = dGluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1Lc = dconjg(GluLt1(isf2,itt2,ftt1,igen2,fgen1))
	  g1Rc = dconjg(GluRt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2Lc = dconjg(g2L)
	  g2Rc = dconjg(g2R)
	

	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-((sqrtS**4*(g1R*(g1Rc*g2R*muugl*
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
! ! 	!      
! 
!      
! ! ! ---------------------t-channel Counter and u-channel Gluino squared------------------------------
! ! 
! 
	  colfac = -1d0/2d0*CF
	
	  g1L = dGluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(GluLt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g1R = dGluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Rc = dconjg(GluRt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(GluLt1(isf1,itt1,ftt2,igen1,fgen2))
	
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Rc = dconjg(GluRt1(isf1,itt1,ftt2,igen1,fgen2))

	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*((-sqrtS**4*(g1R*(g1Rc*g2R*muugl*
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
!      
! !      
! ! ---------------------t-channel and u-channel Counter Gluino squared------------------------------
! 

	  colfac = -1d0/2d0*CF
	
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(dGluLt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Rc = dconjg(dGluRt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(GluLt1(isf1,itt1,ftt2,igen1,fgen2))
	
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Rc = dconjg(GluRt1(isf1,itt1,ftt2,igen1,fgen2))

	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*((-sqrtS**4*(g1R*(g1Rc*g2R*muugl*
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
! ! 
! 
! ! ! ! !-------------------------------------------------------------------------------
! ! ! ! !-------------------------------Neutralino squared------------------------------
! ! ! ! !-------------------------------------------------------------------------------
! ! ! ! 
! ! ! ! !---------------------t-channel Neutralino squared------------------------------
! ! ! !         
      do n=1,6
      do k=1,6
! ! 
! !         
!         
	  colfac = Nc**2
	  g3L = dALt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3R = dARt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))
	
	  g4L = ALt1(k,isf2,itt2,ftt2,igen2,fgen2)
      g4R = ARt1(k,isf2,itt2,ftt2,igen2,fgen2)
	  g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  	


	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*( -((sqrtS**4*
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
! ! 
! !      
! ! ! ! !---------------------u-channel Neutralino squared------------------------------
! ! ! ! ! 
	  colfac = Nc**2
	
	  g3L = dALt1(k,isf2,itt2,ftt1,igen2,fgen1)
	  g3R = dARt1(k,isf2,itt2,ftt1,igen2,fgen1)
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4L = ALt1(k,isf1,itt1,ftt2,igen1,fgen2)
	  g4R = ARt1(k,isf1,itt1,ftt2,igen1,fgen2)
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  
	  
	

	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-((sqrtS**4*(g3R*g3Rc*
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
     
     
! !      
! ! !---------------------t-channel Counter and u-channel Neutralino squared------------------------------
! 
	  colfac = Nc
	
	  g3L = dALt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
! 	
	  g3R = dARt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4L = ALt1(k,isf2,itt2,ftt2,igen2,fgen2)
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
	  g4R = ARt1(k,isf2,itt2,ftt2,igen2,fgen2)
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))

	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*((-sqrtS**4*
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
     
! ! ! !---------------------t-channel and u-channel Counter Neutralino squared------------------------------
! 
	  colfac = Nc
	
	  g3L = ALt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3Lc = dconjg(dALt1(n,isf2,itt2,ftt1,igen2,fgen1))
! 	
	  g3R = ARt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3Rc = dconjg(dARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4L = ALt1(k,isf2,itt2,ftt2,igen2,fgen2)
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
	  g4R = ARt1(k,isf2,itt2,ftt2,igen2,fgen2)
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))

	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0* colfac*((-sqrtS**4*
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
! 
        enddo
! ! ! !        
! ! ! !------------------------------------------------------------------------------------------------------
! ! ! !----------------------------------------INTEFERENZ-PARTS---------------------------------------------- 
! ! ! !------------------------------------------------------------------------------------------------------     
! ! ! 
!---------------------t-channel Counter Gluino t-channel Neutralino--------------------------------------------
	  colfac = 0d0
	  
	  g1L = dGluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = dGluRt1(isf1,itt1,ftt1,igen1,fgen1)
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  
	  g3Rc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Lc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))
	
	  g4Rc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Lc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))

	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-sqrtS**4*(g1R*(2*g2L*g3Lc*g4Rc*muu3*muu4*tred + 
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
! !---------------------t-channel Gluino t-channel Counter Neutralino--------------------------------------------
	  colfac = 0d0
	  
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  
	  g3Rc = dconjg(dALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Lc = dconjg(dARt1(n,isf1,itt1,ftt1,igen1,fgen1))
	
	  g4Rc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Lc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))

	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-sqrtS**4*(g1R*(2*g2L*g3Lc*g4Rc*muu3*muu4*tred + 
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
! ! 	  
! ! !---------------------u-channel Gluino Counter u-channel Neutralino--------------------------------------------------
	  colfac = 0d0
	  
	  g1L = dGluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = dGluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  
	  
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-sqrtS**4*(g1R*(2*g2L*g3Lc*g4Rc*muu3*muu4*
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
        
! ! !---------------------u-channel Gluino u-channel Counter Neutralino--------------------------------------------------
	  colfac = 0d0
	  
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  
	  
	  g3Lc = dconjg(dALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(dARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-sqrtS**4*(g1R*(2*g2L*g3Lc*g4Rc*muu3*muu4*
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
! !---------------------t-channel Counter Gluino u-channel Neutralino--------------------------------------------------------
	  colfac = CF*Nc
	  
	  g1L = dGluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = dGluRt1(isf1,itt1,ftt1,igen1,fgen1)
! 	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
	
	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-sqrtS**4*(g1R*(g2L*g3Lc*g4Rc*muu3*muu4*
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
! !      
! ! !---------------------t-channel Gluino u-channel Counter Neutralino--------------------------------------------------------
	  colfac = CF*Nc
	  
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
! 	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  
	  g3Lc = dconjg(dALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(dARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
	
	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-sqrtS**4*(g1R*(g2L*g3Lc*g4Rc*muu3*muu4*
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
     
! !     
! !         
! !---------------------u-channel Counter Gluino t-channel Neutralino---------------------------------------
	  colfac = CF*Nc
	  
	  g1L = dGluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = dGluRt1(isf2,itt2,ftt1,igen2,fgen1)

	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)

	  g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))

	  g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  
	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-sqrtS**4*(g1R*(g2L*g3Lc*g4Rc*muu3*muu4*
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


! ! !---------------------u-channel Gluino t-channel Counter Neutralino---------------------------------------
	  colfac = CF*Nc
	  
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)

	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)

	  g3Lc = dconjg(dALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(dARt1(n,isf1,itt1,ftt1,igen1,fgen1))

	  g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  
	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-sqrtS**4*(g1R*(g2L*g3Lc*g4Rc*muu3*muu4*
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
! 
        enddo
! 
!         
! c-------------------------------------------------------------------------------------------------------
! c-------------------------------------------------------------------------------------------------------
! c-------------------------------------------------------------------------------------------------------
! c-------------------------------------------------------------------------------------------------------
! c-------------------------------------------------------------------------------------------------------
! c-------------------------------------------------------------------------------------------------------
! c-------------------------------------------------------------------------------------------------------
! c--------------------------------------lower Counterterm------------------------------------------------
! c-------------------------------------------------------------------------------------------------------
! c-------------------------------------------------------------------------------------------------------
! c-------------------------------------------------------------------------------------------------------
! c-------------------------------------------------------------------------------------------------------
! c-------------------------------------------------------------------------------------------------------
! c-------------------------------------------------------------------------------------------------------
! 
! 	          test = stst2qq_M2Counterterm
! 
! 
! ! ! ---------------------t-channel Gluino squared------------------------------
! ! ! ! 
	  colfac = Nc/2d0*CF
	
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(g1L)
	  g1Rc = dconjg(g1R)

	
! 	  g2L = dGluR(isf2,itt2,ftt2,igen2,fgen2)
! 	  g2R = dGluL(isf2,itt2,ftt2,igen2,fgen2)
	  g2L = dGluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = dGluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(GluLt1(isf2,itt2,ftt2,igen2,fgen2))
	  g2Rc = dconjg(GluRt1(isf2,itt2,ftt2,igen2,fgen2))
	

	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*dreal(-((sqrtS**4*(g1R*
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
!      

! !      
! ! ! ! !---------------------u-channel Gluino squared------------------------------
! ! ! ! ! 
	  colfac = Nc/2d0*CF
    
	
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1Lc = dconjg(g1L)
	  g1Rc = dconjg(g1R)
	
	  g2L = dGluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = dGluRt1(isf1,itt1,ftt2,igen1,fgen2)
! 	  g2L = dGluR(isf1,itt1,ftt2,igen1,fgen2)
! 	  g2R = dGluL(isf1,itt1,ftt2,igen1,fgen2)
	  g2Lc = dconjg(GluLt1(isf1,itt1,ftt2,igen1,fgen2))
	  g2Rc = dconjg(GluRt1(isf1,itt1,ftt2,igen1,fgen2))
	

	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-((sqrtS**4*(g1R*(g1Rc*g2R*muugl*
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
     
!      
! ! ! ! ---------------------t-channel Counter and u-channel Gluino squared------------------------------
! ! 
! 
	  colfac = -1d0/2d0*CF
	
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(GluLt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Rc = dconjg(GluRt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g2L = dGluLt1(isf2,itt2,ftt2,igen2,fgen2)
! 	  g2L = dGluR(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(GluLt1(isf1,itt1,ftt2,igen1,fgen2))
	
	  g2R = dGluRt1(isf2,itt2,ftt2,igen2,fgen2)
! 	  g2R = dGluL(isf2,itt2,ftt2,igen2,fgen2)
	  g2Rc = dconjg(GluRt1(isf1,itt1,ftt2,igen1,fgen2))

	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*((-sqrtS**4*(g1R*(g1Rc*g2R*muugl*
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
! ! 
! ! 
! ! ! ---------------------t-channel and u-channel Counter Gluino squared------------------------------
! ! 
! ! ! 
	  colfac = -1d0/2d0*CF
	
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(GluLt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Rc = dconjg(GluRt1(isf2,itt2,ftt1,igen2,fgen1))
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(dGluLt1(isf1,itt1,ftt2,igen1,fgen2))
!       g2Lc = dconjg(dGluR(isf1,itt1,ftt2,igen1,fgen2))

	
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Rc = dconjg(dGluRt1(isf1,itt1,ftt2,igen1,fgen2))
! 	  g2Rc = dconjg(dGluL(isf1,itt1,ftt2,igen1,fgen2))


	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*((-sqrtS**4*(g1R*(g1Rc*g2R*muugl*
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
! 

! ! ! ! 
! ! ! ! ! ! !-------------------------------------------------------------------------------
! ! ! ! ! ! !-------------------------------Neutralino squared------------------------------
! ! ! ! ! ! !-------------------------------------------------------------------------------
! ! ! ! ! ! 
! ! ! ! ! ! !---------------------t-channel Neutralino squared------------------------------
! ! ! ! ! ! !         
      do n=1,6
      do k=1,6
! ! ! 
! !         
! !         
	  colfac = Nc**2
	  g3L = ALt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3R = ARt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))
	
	  g4L = dALt1(k,isf2,itt2,ftt2,igen2,fgen2)
      g4R = dARt1(k,isf2,itt2,ftt2,igen2,fgen2)
	  g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))


	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*( -((sqrtS**4*
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

 
     
! !---------------------u-channel Neutralino squared------------------------------

	  colfac = Nc**2
	
	  g3L = ALt1(k,isf2,itt2,ftt1,igen2,fgen1)
	  g3R = ARt1(k,isf2,itt2,ftt1,igen2,fgen1)
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4L = dALt1(k,isf1,itt1,ftt2,igen1,fgen2)
	  g4R = dARt1(k,isf1,itt1,ftt2,igen1,fgen2)
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
        
	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-((sqrtS**4*(g3R*g3Rc*
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
     
     
! !      
! !      
! !---------------------t-channel Counter and u-channel Neutralino squared------------------------------

	  colfac = Nc
	
	  g3L = ALt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g3R = ARt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4L = dALt1(k,isf2,itt2,ftt2,igen2,fgen2)
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
	  g4R = dARt1(k,isf2,itt2,ftt2,igen2,fgen2)
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))

	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*((-sqrtS**4*
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
! 
! ! 
! !---------------------t-channel and u-channel Counter Neutralino squared------------------------------
! 
	  colfac = Nc
	
	  g3L = ALt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
! 	
	  g3R = ARt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4L = ALt1(k,isf2,itt2,ftt2,igen2,fgen2)
	  g4Lc = dconjg(dALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
	  g4R = ARt1(k,isf2,itt2,ftt2,igen2,fgen2)
	  g4Rc = dconjg(dARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  

	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*((-sqrtS**4*
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
! ! 
        enddo
! ! ! ! !        
! ! ! ! ! !------------------------------------------------------------------------------------------------------
! ! ! ! ! !----------------------------------------INTEFERENZ-PARTS---------------------------------------------- 
! ! ! ! ! !------------------------------------------------------------------------------------------------------     
! ! ! ! ! 
! ! ! ! ! !---------------------t-channel Counter Gluino t-channel Neutralino--------------------------------------------
	  colfac = 0d0
	  
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	
	  g2L = dGluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = dGluRt1(isf2,itt2,ftt2,igen2,fgen2)
! 	  g2L = dGluR(isf2,itt2,ftt2,igen2,fgen2)
! 	  g2R = dGluL(isf2,itt2,ftt2,igen2,fgen2)
! 	  
	  g3Rc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Lc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))
	
	  g4Rc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Lc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))

	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-sqrtS**4*(g1R*(2*g2L*g3Lc*g4Rc*muu3*muu4*tred + 
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
	  
! !---------------------t-channel Gluino t-channel Counter Neutralino--------------------------------------------
	  colfac = 0d0
	  
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  
	  g3Rc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Lc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))
	
	  g4Rc = dconjg(dALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Lc = dconjg(dARt1(n,isf2,itt2,ftt2,igen2,fgen2))

	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-sqrtS**4*(g1R*(2*g2L*g3Lc*g4Rc*muu3*muu4*tred + 
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
! ! 	  
! ! ! !---------------------u-channel Counter Gluino u-channel Neutralino--------------------------------------------------
	  colfac = 0d0
	  
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  
	  g2L = dGluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = dGluRt1(isf1,itt1,ftt2,igen1,fgen2)
! 	  g2L = dGluR(isf1,itt1,ftt2,igen1,fgen2)
! 	  g2R = dGluL(isf1,itt1,ftt2,igen1,fgen2)
! 	  
	  
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-sqrtS**4*(g1R*(2*g2L*g3Lc*g4Rc*muu3*muu4*
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
        
!         
! ! !---------------------u-channel Gluino u-channel Counter Neutralino--------------------------------------------------
	  colfac = 0d0
	  
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  
	  
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(dALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(dARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-sqrtS**4*(g1R*(2*g2L*g3Lc*g4Rc*muu3*muu4*
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
        
! !         
! ! !---------------------t-channel Counter Gluino u-channel Neutralino--------------------------------------------------------
	  colfac = CF*Nc
	  
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
! 	
	  g2L = dGluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = dGluRt1(isf2,itt2,ftt2,igen2,fgen2)
! 	  g2L = dGluR(isf2,itt2,ftt2,igen2,fgen2)
! 	  g2R = dGluL(isf2,itt2,ftt2,igen2,fgen2)
! 	  
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
	
	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0* colfac*(-sqrtS**4*(g1R*(g2L*g3Lc*g4Rc*muu3*muu4*
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
! !     
! !     
! ! !---------------------t-channel Gluino u-channel Counter Neutralino--------------------------------------------------------
	  colfac = CF*Nc
	  
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
! 	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4Lc = dconjg(dALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(dARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	
	
	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-sqrtS**4*(g1R*(g2L*g3Lc*g4Rc*muu3*muu4*
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
! !      
! ! !         
! ! !---------------------u-channel Counter Gluino t-channel Neutralino---------------------------------------
	  colfac = CF*Nc
	  
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)

	  g2L = dGluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = dGluRt1(isf1,itt1,ftt2,igen1,fgen2)
!       g2L = dGluR(isf1,itt1,ftt2,igen1,fgen2)
! 	  g2R = dGluL(isf1,itt1,ftt2,igen1,fgen2)

	  g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))

	  g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  
	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-sqrtS**4*(g1R*(g2L*g3Lc*g4Rc*muu3*muu4*
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
! ! !---------------------u-channel Gluino t-channel Counter Neutralino---------------------------------------
	  colfac = CF*Nc
	  
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)

	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)

	  g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))

	  g4Lc = dconjg(dALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(dARt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  
	  stst2qq_M2Counterterm = stst2qq_M2Counterterm + 2d0*colfac*(-sqrtS**4*(g1R*(g2L*g3Lc*g4Rc*muu3*muu4*
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
! 
! ! ! 
        enddo
! ! 


	end
