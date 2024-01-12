************************************************

	function stst2QQ_M2Boxsextet()
	implicit none

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
#include "stst2QQ_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer i,j,k,l,n,m,b,d,h
	complex*16 g0,g0L,g0R,g1,g1L,g1R,g2,g2L,g2R,g3,g3L,g3R
	double complex g1Lc, g2Lc, g3Lc, g4L, g4Lc
	double complex g1Rc, g2Rc, g3Rc, g4R, g4Rc
	double precision Ml0,Ml1,Ml3
	complex*16 vbA, vbBL, vbBR, vbCLp1, vbCRp1, vbCLp2, vbCRp2, vbDp1, vbDk1, vbDk2
	double precision Krondelta, ef, I3f, hf
	complex*16 test
	
	double complex KappaFun,testx13,ieps
	
	double precision boxfactor,boxfactor1,testtreebox
* functions & subroutines
	double precision colorf, colori,symfac
	double complex stst2QQ_M2Boxsextet
	external myCget, myBget, myDget

	stst2QQ_M2Boxsextet = 0d0
	test = 0d0
	ieps = dcmplx(0d0,1d-17)

! 
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
! !Propagator
! 
! 	do i=1,4 !here distinguish between gluino und neutralino
! 		DenGauT(i)=(tman - MNeu(i)**2) !+ cI*MNeu(i)*wNeu(i)
! 	enddo
! 	DenGauT(5)=(tman - MCha(1)**2) !+ cI*MCha(1)*wCha(1)
! 	DenGauT(6)=(tman - MCha(2)**2) !+ cI*MCha(2)*wCha(2)
! 
! 	DenGlT = (tman-MGl**2)

! 	print*,'--------------------------Box-correction squark identical color and quarks identical color-----------------------------------'
!-------------------------------------------------------------------------------------------------------------------------------------------
!IMPORTANT: Boxen taken from stsT2QQbar_M2Boxijts and modified couplings 0 and 1 

!Splitting boxendiagrams in colorstructure Delta(i,j) Delta(t,s)	
	
	vbA = 0d0
	vbBL = 0d0
	vbBR = 0d0
	vbCLp1 = 0d0
	vbCRp1 = 0d0
	vbCLp2 = 0d0
	vbCRp2 = 0d0
	vbDp1 = 0d0
	vbDk1 = 0d0
	vbDk2 = 0d0
! 	
! 	
!-------------------------------------------------------------------------------------------------
!---------------------------------here QCD boxdiagramms-------------------------------------------
!-------------------------------------------------------------------------------------------------
! ! 
! ! !-------------------------------------------------------------------------------------------------
! ! !--------------------------------------t-channels-------------------------------------------------
! ! !-------------------------------------------------------------------------------------------------
! ! ! 
! ! ! ! ! ****************************************************
! ! ! ! ! ************ Box 1 Gluon-Squark (i,k,m)- Gluino- Squark (j,l,n) ***********************
! ! ! ! ! ****************************************************
! ! ! !Coupling 0 is L R exchanged (GluLt2->GluLt1), Coupling 1 gains minus sign due to opposite charge flow
	  do i = 1,2
	    do j = 1,2
	      do k = 3,4
	        do l = 3,4
              do m = 3,3
                do n = 3,3

          
		  
		  Ml0 = 0d0
		  Ml1 = MSf(i,k,m)
		  Ml2 = MGl
		  Ml3 = MSf(j,l,n)
		  
		  g0L = GluLt1(j,l,ftt2,n,fgen2)
		  g0R = GluRt1(j,l,ftt2,n,fgen2)
		  g1 = VSfSf(1,isf2,j,itt2,l,igen2,n)
		  g2 = VSfSf(1,i,isf1,k,itt1,m,igen1)
		  g3L = GluLt1(i,k,ftt1,m,fgen1) 
		  g3R = GluRt1(i,k,ftt1,m,fgen1)
		  
		  call myCget(muu3**2*sqrtS**2,muu4**2*sqrtS**2,sqrtS**2,Ml1**2,Ml2**2,Ml3**2)
		  call myDget(muu1**2*sqrtS**2,muu3**2*sqrtS**2,muu4**2*sqrtS**2,muu2**2*sqrtS**2,tred*sqrtS**2,sqrtS**2, Ml0**2,Ml1**2,Ml2**2,Ml3**2)

   
	  colori = CF/2d0*(Nc-1d0)
	
	  vbA = vbA + colori*0d0
	  
	  vbBR = vbBR + colori*g1*g2*(g0R*g3R*Ml2*(xc0 + 
     -       sqrtS**2*(2*(-1 + muu1**2 + muu2**2)*xd0 + 
     -          (-1 + 3*muu1**2 + muu2**2)*xd1 + muu1**2*xd2 + 
     -          muu2**2*xd2 - muu3**2*xd2 - muu4**2*xd2 + 2*tred*xd2 - 
     -          xd3 + muu1**2*xd3 + 3*muu2**2*xd3)) + 
     -    g0R*g3L*muu3*sqrtS*
     -     (xc0 + xc1 + xc2 - 2*sqrtS**2*xd0 + 
     -       2*muu1**2*sqrtS**2*xd0 + 2*muu2**2*sqrtS**2*xd0 + 
     -       2*xd00 - sqrtS**2*xd1 + 3*muu1**2*sqrtS**2*xd1 + 
     -       muu2**2*sqrtS**2*xd1 - sqrtS**2*xd12 + 
     -       3*muu1**2*sqrtS**2*xd12 + muu2**2*sqrtS**2*xd12 - 
     -       sqrtS**2*xd13 + 3*muu1**2*sqrtS**2*xd13 + 
     -       muu2**2*sqrtS**2*xd13 - 2*sqrtS**2*xd2 + 
     -       3*muu1**2*sqrtS**2*xd2 + 3*muu2**2*sqrtS**2*xd2 - 
     -       muu3**2*sqrtS**2*xd2 - muu4**2*sqrtS**2*xd2 + 
     -       2*sqrtS**2*tred*xd2 + muu1**2*sqrtS**2*xd22 + 
     -       muu2**2*sqrtS**2*xd22 - muu3**2*sqrtS**2*xd22 - 
     -       muu4**2*sqrtS**2*xd22 + 2*sqrtS**2*tred*xd22 - 
     -       sqrtS**2*xd23 + 2*muu1**2*sqrtS**2*xd23 + 
     -       4*muu2**2*sqrtS**2*xd23 - muu3**2*sqrtS**2*xd23 - 
     -       muu4**2*sqrtS**2*xd23 + 2*sqrtS**2*tred*xd23 - 
     -       3*sqrtS**2*xd3 + 3*muu1**2*sqrtS**2*xd3 + 
     -       5*muu2**2*sqrtS**2*xd3 - sqrtS**2*xd33 + 
     -       muu1**2*sqrtS**2*xd33 + 3*muu2**2*sqrtS**2*xd33) + 
     -    g0L*g3R*muu4*sqrtS*
     -     (-xc2 - 2*xd00 - sqrtS**2*
     -        ((-1 + 3*muu1**2 + muu2**2)*xd13 + muu2**2*xd23 - 
     -          muu3**2*xd23 - muu4**2*xd23 + 2*tred*xd23 - 2*xd3 + 
     -          2*muu2**2*xd3 - xd33 + 3*muu2**2*xd33 + 
     -          muu1**2*(xd23 + 2*xd3 + xd33))))
	  
	  vbBL = vbBL + colori*g1*g2*(g0L*g3L*Ml2*(xc0 + 
     -       sqrtS**2*(2*(-1 + muu1**2 + muu2**2)*xd0 + 
     -          (-1 + 3*muu1**2 + muu2**2)*xd1 + muu1**2*xd2 + 
     -          muu2**2*xd2 - muu3**2*xd2 - muu4**2*xd2 + 2*tred*xd2 - 
     -          xd3 + muu1**2*xd3 + 3*muu2**2*xd3)) + 
     -    g0L*g3R*muu3*sqrtS*
     -     (xc0 + xc1 + xc2 - 2*sqrtS**2*xd0 + 
     -       2*muu1**2*sqrtS**2*xd0 + 2*muu2**2*sqrtS**2*xd0 + 
     -       2*xd00 - sqrtS**2*xd1 + 3*muu1**2*sqrtS**2*xd1 + 
     -       muu2**2*sqrtS**2*xd1 - sqrtS**2*xd12 + 
     -       3*muu1**2*sqrtS**2*xd12 + muu2**2*sqrtS**2*xd12 - 
     -       sqrtS**2*xd13 + 3*muu1**2*sqrtS**2*xd13 + 
     -       muu2**2*sqrtS**2*xd13 - 2*sqrtS**2*xd2 + 
     -       3*muu1**2*sqrtS**2*xd2 + 3*muu2**2*sqrtS**2*xd2 - 
     -       muu3**2*sqrtS**2*xd2 - muu4**2*sqrtS**2*xd2 + 
     -       2*sqrtS**2*tred*xd2 + muu1**2*sqrtS**2*xd22 + 
     -       muu2**2*sqrtS**2*xd22 - muu3**2*sqrtS**2*xd22 - 
     -       muu4**2*sqrtS**2*xd22 + 2*sqrtS**2*tred*xd22 - 
     -       sqrtS**2*xd23 + 2*muu1**2*sqrtS**2*xd23 + 
     -       4*muu2**2*sqrtS**2*xd23 - muu3**2*sqrtS**2*xd23 - 
     -       muu4**2*sqrtS**2*xd23 + 2*sqrtS**2*tred*xd23 - 
     -       3*sqrtS**2*xd3 + 3*muu1**2*sqrtS**2*xd3 + 
     -       5*muu2**2*sqrtS**2*xd3 - sqrtS**2*xd33 + 
     -       muu1**2*sqrtS**2*xd33 + 3*muu2**2*sqrtS**2*xd33) + 
     -    g0R*g3L*muu4*sqrtS*
     -     (-xc2 - 2*xd00 - sqrtS**2*
     -        ((-1 + 3*muu1**2 + muu2**2)*xd13 + muu2**2*xd23 - 
     -          muu3**2*xd23 - muu4**2*xd23 + 2*tred*xd23 - 2*xd3 + 
     -          2*muu2**2*xd3 - xd33 + 3*muu2**2*xd33 + 
     -          muu1**2*(xd23 + 2*xd3 + xd33))))
	  
	  vbCLp1 = vbCLp1 + colori*(-(g0L*g1*g2*g3R*(xc0 + 4*xd00 + 
     -      sqrtS**2*(2*(-1 + muu1**2 + muu2**2)*xd0 + 
     -         (-3 + 5*muu1**2 + 3*muu2**2)*xd1 - xd11 + 
     -         3*muu1**2*xd11 + muu2**2*xd11 - xd12 + 4*muu1**2*xd12 + 
     -         2*muu2**2*xd12 - muu3**2*xd12 - muu4**2*xd12 + 
     -         2*tred*xd12 - 2*xd13 + 4*muu1**2*xd13 + 
     -         4*muu2**2*xd13 - 2*xd2 + 3*muu1**2*xd2 + 
     -         3*muu2**2*xd2 - muu3**2*xd2 - muu4**2*xd2 + 
     -         2*tred*xd2 + muu1**2*xd22 + muu2**2*xd22 - 
     -         muu3**2*xd22 - muu4**2*xd22 + 2*tred*xd22 - xd23 + 
     -         2*muu1**2*xd23 + 4*muu2**2*xd23 - muu3**2*xd23 - 
     -         muu4**2*xd23 + 2*tred*xd23 - 3*xd3 + 3*muu1**2*xd3 + 
     -         5*muu2**2*xd3 - xd33 + muu1**2*xd33 + 3*muu2**2*xd33))))
	  
	  vbCRp1 = vbCRp1 +colori*(-(g0R*g1*g2*g3L*(xc0 + 4*xd00 + 
     -      sqrtS**2*(2*(-1 + muu1**2 + muu2**2)*xd0 + 
     -         (-3 + 5*muu1**2 + 3*muu2**2)*xd1 - xd11 + 
     -         3*muu1**2*xd11 + muu2**2*xd11 - xd12 + 4*muu1**2*xd12 + 
     -         2*muu2**2*xd12 - muu3**2*xd12 - muu4**2*xd12 + 
     -         2*tred*xd12 - 2*xd13 + 4*muu1**2*xd13 + 
     -         4*muu2**2*xd13 - 2*xd2 + 3*muu1**2*xd2 + 
     -         3*muu2**2*xd2 - muu3**2*xd2 - muu4**2*xd2 + 
     -         2*tred*xd2 + muu1**2*xd22 + muu2**2*xd22 - 
     -         muu3**2*xd22 - muu4**2*xd22 + 2*tred*xd22 - xd23 + 
     -         2*muu1**2*xd23 + 4*muu2**2*xd23 - muu3**2*xd23 - 
     -         muu4**2*xd23 + 2*tred*xd23 - 3*xd3 + 3*muu1**2*xd3 + 
     -         5*muu2**2*xd3 - xd33 + muu1**2*xd33 + 3*muu2**2*xd33))))
     
	  vbCLp2 = vbCLp2 + colori*0d0
	  vbCRp2 = vbCRp2 + colori*0d0
	  vbDp1 = vbDp1 + colori*0d0
	  vbDk1 = vBDk1 + colori*0d0
	  vbDk2 = vbDk2 + colori*0d0
	
		        enddo
		      enddo
	        enddo
	      enddo
	    enddo
	  enddo

! ! ! ! ! ! ! ****************************************************
! ! ! ! ! ! ! ************ Box 1.2 Gluon-Squark (i,k,m)- Gluino- Quark (l,n) ***********************
! ! ! ! ! ! ! ****************************************************
! ! ! !Coupling 1 is L R exchanged (GluLt2->GluLt1), Coupling 0 gains a minus sign (here minus times minus results in positiv coupling)
	  do i = 1,2
       do k = 3,4
	    do l = 3,4
	     do m = 3,3
		  do n = 3,3

		  Ml0 = 0d0
		  Ml1 = MSf(i,k,m)
		  Ml2 = MGl
		  Ml3 = Mf(l,n)
! 		  
		  g0L = GS*Krondelta(l,ftt2)*Krondelta(n,fgen2)
		  g0R = GS*Krondelta(l,ftt2)*Krondelta(n,fgen2)
		  g1L = GluLt1(isf2,itt2,l,igen2,n) 
		  g1R = GluRt1(isf2,itt2,l,igen2,n)
		  g2 = VSfSf(1,i,isf1,k,itt1,m,igen1)
		  g3L = GluLt1(i,k,ftt1,m,fgen1) 
		  g3R = GluRt1(i,k,ftt1,m,fgen1)
		  
! 		
		  call myCget(muu3**2*sqrtS**2,ured*sqrtS**2,muu2*sqrtS**2,Ml1**2,Ml2**2,Ml3**2)
		  call myDget(muu1**2*sqrtS**2,muu3**2*sqrtS**2,muu2**2*sqrtS**2,muu4**2*sqrtS**2,tred*sqrtS**2,ured*sqrtS**2, Ml0**2,Ml1**2,Ml2**2,Ml3**2)
		  		  
	colori = CF*Nc*(1d0/2d0+CF)
	
	vbA = vbA + colori*0d0
	
	vbBL = vbBL + colori*g2*(g0L*(g1L*(g3R*muu3*sqrtS*(xc0 + xc1 + 2*xc2 + Ml0**2*xd0 + 2*sqrtS**2*xd0 - 2*muu2**2*sqrtS**2*xd0 - 2*muu3**2*sqrtS**2*xd0 + 2*sqrtS**2*tred*xd0 + sqrtS**2*xd1 + 
     -             2*muu1**2*sqrtS**2*xd1 - muu2**2*sqrtS**2*xd1 - muu3**2*sqrtS**2*xd1 + sqrtS**2*tred*xd1 + sqrtS**2*xd12 - muu2**2*sqrtS**2*xd12 - muu3**2*sqrtS**2*xd12 + 
     -             sqrtS**2*tred*xd12 + Ml0**2*xd2 + 2*sqrtS**2*xd2 - 3*muu2**2*sqrtS**2*xd2 - 2*muu3**2*sqrtS**2*xd2 + muu4**2*sqrtS**2*xd2 + 3*sqrtS**2*tred*xd2 - 
     -             muu2**2*sqrtS**2*xd22 + muu4**2*sqrtS**2*xd22 + sqrtS**2*tred*xd22 + 2*muu4**2*sqrtS**2*xd23 + 2*sqrtS**2*xd3 - 2*muu2**2*sqrtS**2*xd3 - 2*muu3**2*sqrtS**2*xd3 + 
     -             muu4**2*sqrtS**2*xd3 + 2*sqrtS**2*tred*xd3) + g3L*Ml2*(xc0 + Ml0**2*xd0 + 
     -             sqrtS**2*(-2*(-1 + muu2**2 + muu3**2 - tred)*xd0 + (1 + 2*muu1**2 - muu2**2 - muu3**2 + tred)*xd1 + 2*muu1**2*xd2 - muu2**2*xd2 + muu4**2*xd2 + tred*xd2 + 2*xd3 - 
     -                2*muu2**2*xd3 - 2*muu3**2*xd3 + muu4**2*xd3 + 2*tred*xd3))) + 
     -       g1R*Ml3*(g3R*Ml2*muu3*sqrtS*xd2 + g3L*(xc0 + Ml0**2*xd0 + sqrtS**2*(muu1**2*(2*xd0 + 3*xd1 + 2*xd2) - 2*(-1 + muu2**2 + muu3**2)*xd3 + tred*(xd2 + 2*xd3))))) + 
     -    g0R*muu4*sqrtS*(g1L*Ml3*(g3L*Ml2 + g3R*muu3*sqrtS)*xd3 - g1R*(g3R*Ml2*muu3*sqrtS*xd2 + 
     -          g3L*(xc0 + xc2 + Ml0**2*xd0 + 2*muu1**2*sqrtS**2*xd0 - 2*xd00 + 3*muu1**2*sqrtS**2*xd1 - sqrtS**2*xd13 + muu2**2*sqrtS**2*xd13 + muu3**2*sqrtS**2*xd13 - 
     -             sqrtS**2*tred*xd13 + 2*muu1**2*sqrtS**2*xd2 + sqrtS**2*tred*xd2 + muu2**2*sqrtS**2*xd23 - muu4**2*sqrtS**2*xd23 - sqrtS**2*tred*xd23 - Ml0**2*xd3 + 
     -             2*muu1**2*sqrtS**2*xd3 - 2*muu4**2*sqrtS**2*xd33))))
     
     
	vbBR = vbBR + colori*(g2*(g0R*(-(g1L*Ml3*(-(g3L*Ml2*muu3*sqrtS*xd2) + g3R*(xc0 + Ml0**2*xd0 + sqrtS**2*(muu1**2*(-2*xd0 + xd1 + 2*xd2) - tred*(xd2 - 2*xd3) - 2*(-1 + muu2**2 + muu3**2)*xd3)))) + 
     -         g1R*(g3L*muu3*sqrtS*(xc0 - xc1 - 2*xc2 + Ml0**2*xd0 + 2*sqrtS**2*xd0 - 2*muu2**2*sqrtS**2*xd0 - 2*muu3**2*sqrtS**2*xd0 + 2*sqrtS**2*tred*xd0 + sqrtS**2*xd1 + 
     -               2*muu1**2*sqrtS**2*xd1 - muu2**2*sqrtS**2*xd1 - muu3**2*sqrtS**2*xd1 + sqrtS**2*tred*xd1 - sqrtS**2*xd12 + muu2**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 - 
     -               sqrtS**2*tred*xd12 - Ml0**2*xd2 - 2*sqrtS**2*xd2 + muu2**2*sqrtS**2*xd2 + 2*muu3**2*sqrtS**2*xd2 + muu4**2*sqrtS**2*xd2 - sqrtS**2*tred*xd2 + muu2**2*sqrtS**2*xd22 - 
     -               muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 - 2*muu4**2*sqrtS**2*xd23 + 2*sqrtS**2*xd3 - 2*muu2**2*sqrtS**2*xd3 - 2*muu3**2*sqrtS**2*xd3 + muu4**2*sqrtS**2*xd3 + 
     -               2*sqrtS**2*tred*xd3) + g3R*Ml2*(xc0 + Ml0**2*xd0 + sqrtS**2*
     -                (-2*(-1 + muu2**2 + muu3**2 - tred)*xd0 + (1 + 2*muu1**2 - muu2**2 - muu3**2 + tred)*xd1 + 2*muu1**2*xd2 - muu2**2*xd2 + muu4**2*xd2 + tred*xd2 + 2*xd3 - 
     -                  2*muu2**2*xd3 - 2*muu3**2*xd3 + muu4**2*xd3 + 2*tred*xd3)))) + 
     -      g0L*muu4*sqrtS*(g1R*Ml3*(g3R*Ml2 + g3L*muu3*sqrtS)*xd3 + g1L*(-(g3L*Ml2*muu3*sqrtS*xd2) + 
     -            g3R*(xc0 + xc2 + Ml0**2*xd0 - 2*muu1**2*sqrtS**2*xd0 - 2*xd00 + muu1**2*sqrtS**2*xd1 - sqrtS**2*xd13 + muu2**2*sqrtS**2*xd13 + muu3**2*sqrtS**2*xd13 - 
     -               sqrtS**2*tred*xd13 + 2*muu1**2*sqrtS**2*xd2 - sqrtS**2*tred*xd2 + muu2**2*sqrtS**2*xd23 - muu4**2*sqrtS**2*xd23 - sqrtS**2*tred*xd23 - Ml0**2*xd3 - 
     -               2*muu1**2*sqrtS**2*xd3 - 2*muu4**2*sqrtS**2*xd33)))))
     
	vbCLp1 = vbCLp1 + colori*(g2*(g0R*muu4*sqrtS*(-3d0*g1L*g3R*Ml3*xd3 + g1R*(g3L*muu3*sqrtS*(2*xd0 + xd1 - 2*xd2 + 2*xd3) + g3R*Ml2*(2*xd0 + xd1 + xd2 + 2*xd3))) + 
     -      g0L*(-(g1R*Ml3*(g3L*muu3*sqrtS*(2*xd0 + xd1 - 2*xd2) + g3R*Ml2*(2*xd0 + xd1 + xd2))) + 
     -         g1L*(-2*g3L*Ml2*muu3*sqrtS*xd2 + g3R*(xc2 + Ml0**2*(xd0 + xd1 + xd2) + 
     -               sqrtS**2*(2*(-1 + muu2**2 + muu3**2 - tred)*xd0 - (-1 + 2*muu1**2 + muu2**2 + muu3**2 - tred)*xd1 + xd11 - muu2**2*xd11 - muu3**2*xd11 + tred*xd11 + xd12 - 
     -                  2*muu2**2*xd12 - muu3**2*xd12 + muu4**2*xd12 + 2*tred*xd12 + 2*muu4**2*xd13 + 2*xd2 - muu2**2*xd2 - 2*muu3**2*xd2 - muu4**2*xd2 - tred*xd2 - muu2**2*xd22 + 
     -                  muu4**2*xd22 + tred*xd22 + 2*muu4**2*xd23 - 2*xd3 + 2*muu2**2*xd3 + 2*muu3**2*xd3 + muu4**2*xd3 - 2*tred*xd3))))))
     
	vbCRp1 = vbCRp1 + colori*(g2*(g0L*muu4*sqrtS*(-3d0*g1R*g3L*Ml3*xd3 + g1L*(g3R*muu3*sqrtS*(2*xd0 + xd1 - 2*xd2 + 2*xd3) + g3L*Ml2*(2*xd0 + xd1 + xd2 + 2*xd3))) + 
     -      g0R*(-(g1L*Ml3*(g3R*muu3*sqrtS*(2*xd0 + xd1 - 2*xd2) + g3L*Ml2*(2*xd0 + xd1 + xd2))) + 
     -         g1R*(-2*g3R*Ml2*muu3*sqrtS*xd2 + g3L*(xc2 + Ml0**2*(xd0 + xd1 + xd2) + 
     -               sqrtS**2*(2*(-1 + muu2**2 + muu3**2 - tred)*xd0 - (-1 + 2*muu1**2 + muu2**2 + muu3**2 - tred)*xd1 + xd11 - muu2**2*xd11 - muu3**2*xd11 + tred*xd11 + xd12 - 
     -                  2*muu2**2*xd12 - muu3**2*xd12 + muu4**2*xd12 + 2*tred*xd12 + 2*muu4**2*xd13 + 2*xd2 - muu2**2*xd2 - 2*muu3**2*xd2 - muu4**2*xd2 - tred*xd2 - muu2**2*xd22 + 
     -                  muu4**2*xd22 + tred*xd22 + 2*muu4**2*xd23 - 2*xd3 + 2*muu2**2*xd3 + 2*muu3**2*xd3 + muu4**2*xd3 - 2*tred*xd3))))))
!      
	  vbCLp2 = vbCLp2 + colori*0d0
	  vbCRp2 = vbCRp2 + colori*0d0
	  vbDp1 = vbDp1 + colori*0d0
	  vbDk1 = vBDk1 + colori*0d0
	  vbDk2 = vbDk2 + colori*0d0
	
	      enddo
	     enddo
	    enddo
	   enddo
	  enddo
! 
! ! ! ! ****************************************************
! ! ! ! ************ Box 1.3 Gluino-Quark (k,m)- Gluon- Squark (j,l,n) ***********************
! ! ! ! ****************************************************
! !Coupling 0 is L R exchanged (GluLt2->GluLt1),Coupling 1 gains minus sign due to opposite charge flow
	do j = 1,2
	    do k = 3,4
	      do l = 3,4
		do m = 3,3
		  do n = 3,3
		  
		  Ml0 = 0d0
		  Ml1 = MSf(j,l,n)
		  Ml2 = MGl
		  Ml3 = Mf(k,m)
! 		  
		  
		  g0L = GluLt1(j,l,ftt2,n,fgen2) 
		  g0R = GluRt1(j,l,ftt2,n,fgen2)
		  g1 = -VSfSf(1,isf2,j,itt2,l,igen2,n)
		  g2L = GluLt1(isf1,itt1,k,igen1,m) 
		  g2R = GluRt1(isf1,itt1,k,igen1,m)
		  g3L = -GS*Krondelta(k,ftt1)*Krondelta(m,fgen1)
		  g3R = -GS*Krondelta(k,ftt1)*Krondelta(m,fgen1)
! 		
		  call myCget(muu3**2*sqrtS**2,ured*sqrtS**2,muu2*sqrtS**2,Ml1**2,Ml2**2,Ml3**2)
		  call myDget(muu2**2*sqrtS**2,muu4**2*sqrtS**2,muu1**2*sqrtS**2,muu3**2*sqrtS**2,tred*sqrtS**2,ured*sqrtS**2, Ml0**2,Ml1**2,Ml2**2,Ml3**2)
		  		  
	colori = CF*Nc*(1d0/2d0+CF)

	vbA = vbA + colori*0d0
	
	vbBL = vbBL + colori*g1*(g0R*muu4*sqrtS*(g2R*Ml2*(-(g3R*Ml3*(2*xd0 + xd1)) + g3L*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3))) - 
     -       g2L*(g3R*(xc0 - xc1 - xc2 + 2*Ml0**2*xd0 + 2*muu3**2*sqrtS**2*xd0 + Ml0**2*xd1 + 2*sqrtS**2*xd1 - 2*muu1**2*sqrtS**2*xd1 + muu3**2*sqrtS**2*xd1 - 2*muu4**2*sqrtS**2*xd1 + 
     -             2*sqrtS**2*tred*xd1 + sqrtS**2*xd11 - muu1**2*sqrtS**2*xd11 - muu4**2*sqrtS**2*xd11 + sqrtS**2*tred*xd11 - muu1**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 + 
     -             sqrtS**2*tred*xd12 + 2*muu3**2*sqrtS**2*xd13 + 2*muu3**2*sqrtS**2*xd2 + 2*sqrtS**2*tred*xd2 + 4*muu3**2*sqrtS**2*xd3) - g3L*Ml3*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3)))
     -       ) + g0L*(g2L*Ml2*(g3L*Ml3*muu3*sqrtS*(2*xd0 + xd1 + xd2 + xd3) + 
     -          g3R*(xc0 + Ml0**2*xd0 - sqrtS**2*(2*(-1 + muu1**2 + muu3**2 + muu4**2 - tred)*xd0 + (-1 + muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 - tred)*xd1 + muu1**2*xd2 - 
     -                2*muu2**2*xd2 + 2*muu4**2*xd2 - tred*xd2 - 2*xd3 + 2*muu1**2*xd3 + muu3**2*xd3 + 2*muu4**2*xd3 - 2*tred*xd3))) + 
     -       g2R*(g3R*Ml3*(xc0 + Ml0**2*xd0 + sqrtS**2*(tred*xd2 + muu2**2*(2*xd0 + 3*xd1 + 2*xd2) + 2*xd3 - 2*muu1**2*xd3 - muu3**2*xd3 + 2*tred*xd3 - 
     -                muu4**2*(2*xd0 + xd1 + 2*(xd2 + xd3)))) + g3L*muu3*sqrtS*
     -           (xc0 + 2*sqrtS**2*xd0 - 2*muu1**2*sqrtS**2*xd0 - 2*muu2**2*sqrtS**2*xd0 + 2*sqrtS**2*tred*xd0 + 2*xd00 + 3*sqrtS**2*xd1 - 3*muu1**2*sqrtS**2*xd1 - muu2**2*sqrtS**2*xd1 - 
     -             2*muu4**2*sqrtS**2*xd1 + 3*sqrtS**2*tred*xd1 + sqrtS**2*xd11 - muu1**2*sqrtS**2*xd11 - muu4**2*sqrtS**2*xd11 + sqrtS**2*tred*xd11 + sqrtS**2*xd12 - 
     -             2*muu1**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 - muu4**2*sqrtS**2*xd12 + 2*sqrtS**2*tred*xd12 + sqrtS**2*xd13 - muu1**2*sqrtS**2*xd13 + 2*muu3**2*sqrtS**2*xd13 - 
     -             muu4**2*sqrtS**2*xd13 + sqrtS**2*tred*xd13 + 2*sqrtS**2*xd2 - 3*muu1**2*sqrtS**2*xd2 - 2*muu2**2*sqrtS**2*xd2 + muu3**2*sqrtS**2*xd2 + 4*sqrtS**2*tred*xd2 - 
     -             muu1**2*sqrtS**2*xd22 + muu3**2*sqrtS**2*xd22 + sqrtS**2*tred*xd22 - muu1**2*sqrtS**2*xd23 + 3*muu3**2*sqrtS**2*xd23 + sqrtS**2*tred*xd23 + 2*sqrtS**2*xd3 - 
     -             2*muu1**2*sqrtS**2*xd3 - 2*muu2**2*sqrtS**2*xd3 + 3*muu3**2*sqrtS**2*xd3 + 2*sqrtS**2*tred*xd3 + Ml0**2*(2*xd0 + xd1 + xd2 + xd3) + 2*muu3**2*sqrtS**2*xd33))))
	
	
     
	vbBR = vbBR + colori*g1*(g0L*muu4*sqrtS*(g2L*Ml2*(-(g3L*Ml3*(2*xd0 + xd1)) + g3R*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3))) - 
     -       g2R*(g3L*(xc0 - xc1 - xc2 + 2*Ml0**2*xd0 + 2*muu3**2*sqrtS**2*xd0 + Ml0**2*xd1 + 2*sqrtS**2*xd1 - 2*muu1**2*sqrtS**2*xd1 + muu3**2*sqrtS**2*xd1 - 2*muu4**2*sqrtS**2*xd1 + 
     -             2*sqrtS**2*tred*xd1 + sqrtS**2*xd11 - muu1**2*sqrtS**2*xd11 - muu4**2*sqrtS**2*xd11 + sqrtS**2*tred*xd11 - muu1**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 + 
     -             sqrtS**2*tred*xd12 + 2*muu3**2*sqrtS**2*xd13 + 2*muu3**2*sqrtS**2*xd2 + 2*sqrtS**2*tred*xd2 + 4*muu3**2*sqrtS**2*xd3) - g3R*Ml3*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3)))
     -       ) + g0R*(g2R*Ml2*(g3R*Ml3*muu3*sqrtS*(2*xd0 + xd1 + xd2 + xd3) + 
     -          g3L*(xc0 + Ml0**2*xd0 - sqrtS**2*(2*(-1 + muu1**2 + muu3**2 + muu4**2 - tred)*xd0 + (-1 + muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 - tred)*xd1 + muu1**2*xd2 - 
     -                2*muu2**2*xd2 + 2*muu4**2*xd2 - tred*xd2 - 2*xd3 + 2*muu1**2*xd3 + muu3**2*xd3 + 2*muu4**2*xd3 - 2*tred*xd3))) + 
     -       g2L*(g3L*Ml3*(xc0 + Ml0**2*xd0 + sqrtS**2*(tred*xd2 + muu2**2*(2*xd0 + 3*xd1 + 2*xd2) + 2*xd3 - 2*muu1**2*xd3 - muu3**2*xd3 + 2*tred*xd3 - 
     -                muu4**2*(2*xd0 + xd1 + 2*(xd2 + xd3)))) + g3R*muu3*sqrtS*
     -           (xc0 + 2*sqrtS**2*xd0 - 2*muu1**2*sqrtS**2*xd0 - 2*muu2**2*sqrtS**2*xd0 + 2*sqrtS**2*tred*xd0 + 2*xd00 + 3*sqrtS**2*xd1 - 3*muu1**2*sqrtS**2*xd1 - muu2**2*sqrtS**2*xd1 - 
     -             2*muu4**2*sqrtS**2*xd1 + 3*sqrtS**2*tred*xd1 + sqrtS**2*xd11 - muu1**2*sqrtS**2*xd11 - muu4**2*sqrtS**2*xd11 + sqrtS**2*tred*xd11 + sqrtS**2*xd12 - 
     -             2*muu1**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 - muu4**2*sqrtS**2*xd12 + 2*sqrtS**2*tred*xd12 + sqrtS**2*xd13 - muu1**2*sqrtS**2*xd13 + 2*muu3**2*sqrtS**2*xd13 - 
     -             muu4**2*sqrtS**2*xd13 + sqrtS**2*tred*xd13 + 2*sqrtS**2*xd2 - 3*muu1**2*sqrtS**2*xd2 - 2*muu2**2*sqrtS**2*xd2 + muu3**2*sqrtS**2*xd2 + 4*sqrtS**2*tred*xd2 - 
     -             muu1**2*sqrtS**2*xd22 + muu3**2*sqrtS**2*xd22 + sqrtS**2*tred*xd22 - muu1**2*sqrtS**2*xd23 + 3*muu3**2*sqrtS**2*xd23 + sqrtS**2*tred*xd23 + 2*sqrtS**2*xd3 - 
     -             2*muu1**2*sqrtS**2*xd3 - 2*muu2**2*sqrtS**2*xd3 + 3*muu3**2*sqrtS**2*xd3 + 2*sqrtS**2*tred*xd3 + Ml0**2*(2*xd0 + xd1 + xd2 + xd3) + 2*muu3**2*sqrtS**2*xd33))))


     
	vbCLp1 = vbCLp1 - colori*(g1*(g0L*(g2L*Ml2*(g3L*Ml3*(2*xd0 + xd1 + xd2) - g3R*muu3*sqrtS*(2*xd0 + xd1 + xd2 + 2*xd3)) + 
     -         g2R*(-(g3R*Ml3*muu3*sqrtS*xd3) + g3L*(2*xc0 + xc2 + 3*Ml0**2*xd0 + 2*sqrtS**2*xd0 - 2*muu1**2*sqrtS**2*xd0 - 2*muu4**2*sqrtS**2*xd0 + 2*sqrtS**2*tred*xd0 + Ml0**2*xd1 + 
     -               3*sqrtS**2*xd1 - 3*muu1**2*sqrtS**2*xd1 + 2*muu2**2*sqrtS**2*xd1 - 3*muu4**2*sqrtS**2*xd1 + 3*sqrtS**2*tred*xd1 + sqrtS**2*xd11 - muu1**2*sqrtS**2*xd11 - 
     -               muu4**2*sqrtS**2*xd11 + sqrtS**2*tred*xd11 + sqrtS**2*xd12 - 2*muu1**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 - muu4**2*sqrtS**2*xd12 + 2*sqrtS**2*tred*xd12 + 
     -               2*muu3**2*sqrtS**2*xd13 + Ml0**2*xd2 + 2*sqrtS**2*xd2 - 3*muu1**2*sqrtS**2*xd2 + muu3**2*sqrtS**2*xd2 - 2*muu4**2*sqrtS**2*xd2 + 5*sqrtS**2*tred*xd2 - 
     -               muu1**2*sqrtS**2*xd22 + muu3**2*sqrtS**2*xd22 + sqrtS**2*tred*xd22 + 2*muu3**2*sqrtS**2*xd23 + 2*sqrtS**2*xd3 - 2*muu1**2*sqrtS**2*xd3 + 3*muu3**2*sqrtS**2*xd3 - 
     -               2*muu4**2*sqrtS**2*xd3 + 2*sqrtS**2*tred*xd3))) + g0R*muu4*sqrtS*
     -       (2*g2R*g3L*Ml2*xd2 + g2L*(g3L*Ml3*(2*xd0 + xd1 + 2*xd2) - g3R*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3))))))
     
	vbCRp1 = vbCRp1 - colori*(g1*(g0R*(g2R*Ml2*(g3R*Ml3*(2*xd0 + xd1 + xd2) - g3L*muu3*sqrtS*(2*xd0 + xd1 + xd2 + 2*xd3)) + 
     -         g2L*(-(g3L*Ml3*muu3*sqrtS*xd3) + g3R*(2*xc0 + xc2 + 3*Ml0**2*xd0 + 2*sqrtS**2*xd0 - 2*muu1**2*sqrtS**2*xd0 - 2*muu4**2*sqrtS**2*xd0 + 2*sqrtS**2*tred*xd0 + Ml0**2*xd1 + 
     -               3*sqrtS**2*xd1 - 3*muu1**2*sqrtS**2*xd1 + 2*muu2**2*sqrtS**2*xd1 - 3*muu4**2*sqrtS**2*xd1 + 3*sqrtS**2*tred*xd1 + sqrtS**2*xd11 - muu1**2*sqrtS**2*xd11 - 
     -               muu4**2*sqrtS**2*xd11 + sqrtS**2*tred*xd11 + sqrtS**2*xd12 - 2*muu1**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 - muu4**2*sqrtS**2*xd12 + 2*sqrtS**2*tred*xd12 + 
     -               2*muu3**2*sqrtS**2*xd13 + Ml0**2*xd2 + 2*sqrtS**2*xd2 - 3*muu1**2*sqrtS**2*xd2 + muu3**2*sqrtS**2*xd2 - 2*muu4**2*sqrtS**2*xd2 + 5*sqrtS**2*tred*xd2 - 
     -               muu1**2*sqrtS**2*xd22 + muu3**2*sqrtS**2*xd22 + sqrtS**2*tred*xd22 + 2*muu3**2*sqrtS**2*xd23 + 2*sqrtS**2*xd3 - 2*muu1**2*sqrtS**2*xd3 + 3*muu3**2*sqrtS**2*xd3 - 
     -               2*muu4**2*sqrtS**2*xd3 + 2*sqrtS**2*tred*xd3))) + g0L*muu4*sqrtS*
     -       (2*g2L*g3R*Ml2*xd2 + g2R*(g3R*Ml3*(2*xd0 + xd1 + 2*xd2) - g3L*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3))))))
     
!      
	  vbCLp2 = vbCLp2 + colori*0d0
	  vbCRp2 = vbCRp2 + colori*0d0
	  vbDp1 = vbDp1 + colori*0d0
	  vbDk1 = vBDk1 + colori*0d0
	  vbDk2 = vbDk2 + colori*0d0
	
		  enddo
		enddo
	      enddo
	    enddo
	enddo
! 	
! 
! ! ****************************************************
! ! ************ Box 4 Gluino-Quark(k,m)-Gluon- Quark(l,n) ***********************
! ! ****************************************************
! ! !Coupling 1 is L R exchanged (GluLt2->GluLt1), Coupling 0 gains a minus sign (here minus times minus results in positiv coupling)
! 
	  do k = 3,4
       do l = 3,4
        do m = 3,3
		 do n = 3,3
		
		  Ml0 = MGl
		  Ml1 = Mf(k,m)
		  Ml2 = 0
		  Ml3 = Mf(l,n)
		  
		  g0 = GS*Krondelta(l,ftt2)*Krondelta(n,fgen2)
		  g1L = GluLt1(isf2,itt2,l,igen2,n)
		  g1R = GluRt1(isf2,itt2,l,igen2,n)
		  g2L = GluLt1(isf1,itt1,k,igen1,m)
		  g2R = GluRt1(isf1,itt1,k,igen1,m)
		  g3 = -GS*Krondelta(k,ftt1)*Krondelta(m,fgen1)
		  
		  call myCget(muu3**2*sqrtS**2,muu4**2*sqrtS**2,sqrtS**2,Ml1**2,Ml2**2,Ml3**2)
		  call myDget(muu1**2*sqrtS**2,muu3**2*sqrtS**2,muu4**2*sqrtS**2,muu2**2*sqrtS**2,tred*sqrtS**2,sqrtS**2, Ml0**2,Ml1**2,Ml2**2,Ml3**2)
		  		  
	colori = CF/2d0*(Nc-1d0)

	vbA = vbA + colori*0d0
	
	vbBL = vbBL - colori*2*g0*g3*((g1R*(g2R*Ml0*(Ml3*muu3*sqrtS*(xd2 + xd3) + Ml1*(2*Ml3*xd0 - muu4*sqrtS*(xd0 + xd3))) + 
     -       g2L*(-(muu4*sqrtS*(xc0 + xc2 + Ml0**2*xd0 + muu1**2*sqrtS**2*xd1 + muu1**2*sqrtS**2*xd2 + 
     -               Ml0**2*xd3)) + Ml3*(2*xc0 + 2*Ml0**2*xd0 + 
     -             sqrtS*(Ml1*muu3*(xd2 + xd3) + muu1**2*sqrtS*(2*xd1 + xd2 + xd3) + 
     -                sqrtS*(-(muu3**2*xd2) + tred*xd2 + (-1 + muu2**2)*xd3))))) + 
     -    g1L*(g2L*Ml0*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*
     -           (-(Ml3*muu4*xd3) + Ml1*muu3*(xd0 + xd2 + xd3) + 
     -             sqrtS*((-1 + muu1**2 + muu2**2)*xd0 + (-1 + 3*muu1**2 + muu2**2)*xd1 + muu1**2*xd2 + 
     -                muu2**2*xd2 - muu3**2*xd2 - muu4**2*xd2 + 2*tred*xd2 - xd3 + muu1**2*xd3 + 3*muu2**2*xd3)))
     -        + g2R*(muu3*sqrtS*(xc0 + xc1 + xc2 + Ml0**2*xd0 + muu1**2*sqrtS**2*xd1 + Ml0**2*xd2 + sqrtS**2*xd2 - 
     -             muu2**2*sqrtS**2*xd2 - muu3**2*sqrtS**2*xd2 + sqrtS**2*tred*xd2 + Ml0**2*xd3) + 
     -          Ml1*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*
     -              (-(Ml3*muu4*xd3) + sqrtS*((-1 + muu1**2 + muu2**2)*xd1 + (-muu4**2 + tred)*xd2 + 
     -                   muu2**2*(xd2 + 2*xd3))))))))
          
	vbBR = vbBR - colori*2*g0*g3*(g1L*(g2L*Ml0*(Ml3*muu3*sqrtS*(xd2 + xd3) + Ml1*(2*Ml3*xd0 - muu4*sqrtS*(xd0 + xd3))) + 
     -       g2R*(-(muu4*sqrtS*(xc0 + xc2 + Ml0**2*xd0 + muu1**2*sqrtS**2*xd1 + muu1**2*sqrtS**2*xd2 + 
     -               Ml0**2*xd3)) + Ml3*(2*xc0 + 2*Ml0**2*xd0 + 
     -             sqrtS*(Ml1*muu3*(xd2 + xd3) + muu1**2*sqrtS*(2*xd1 + xd2 + xd3) + 
     -                sqrtS*(-(muu3**2*xd2) + tred*xd2 + (-1 + muu2**2)*xd3))))) + 
     -    g1R*(g2R*Ml0*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*
     -           (-(Ml3*muu4*xd3) + Ml1*muu3*(xd0 + xd2 + xd3) + 
     -             sqrtS*((-1 + muu1**2 + muu2**2)*xd0 + (-1 + 3*muu1**2 + muu2**2)*xd1 + muu1**2*xd2 + 
     -                muu2**2*xd2 - muu3**2*xd2 - muu4**2*xd2 + 2*tred*xd2 - xd3 + muu1**2*xd3 + 3*muu2**2*xd3)))
     -        + g2L*(muu3*sqrtS*(xc0 + xc1 + xc2 + Ml0**2*xd0 + muu1**2*sqrtS**2*xd1 + Ml0**2*xd2 + sqrtS**2*xd2 - 
     -             muu2**2*sqrtS**2*xd2 - muu3**2*sqrtS**2*xd2 + sqrtS**2*tred*xd2 + Ml0**2*xd3) + 
     -          Ml1*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*
     -              (-(Ml3*muu4*xd3) + sqrtS*((-1 + muu1**2 + muu2**2)*xd1 + (-muu4**2 + tred)*xd2 + 
     -                   muu2**2*(xd2 + 2*xd3)))))))
	
	vbCLp1 = vbCLp1 + colori* 2*g0*g3*(g1R*(g2R*Ml0*Ml3*(xd0 + xd1 + xd2 + xd3) + 
     -       g2L*(-(muu3*muu4*sqrtS**2*xd2) + Ml1*Ml3*(xd1 + xd2 + xd3))) + 
     -    g1L*(g2L*Ml0*Ml1*(xd0 + xd1 + xd2 + xd3) + 
     -       g2R*(xc0 + Ml0**2*(2*xd0 + xd1 + xd2 + xd3) + 
     -          sqrtS**2*(muu1**2*xd1 - (-1 + muu3**2 + muu4**2 - tred)*xd2 + muu2**2*xd3))))
	
	vbCRp1 = vbCRp1 + colori* 2*g0*g3*(g1L*(g2L*Ml0*Ml3*(xd0 + xd1 + xd2 + xd3) + 
     -       g2R*(-(muu3*muu4*sqrtS**2*xd2) + Ml1*Ml3*(xd1 + xd2 + xd3))) + 
     -    g1R*(g2R*Ml0*Ml1*(xd0 + xd1 + xd2 + xd3) + 
     -       g2L*(xc0 + Ml0**2*(2*xd0 + xd1 + xd2 + xd3) + 
     -          sqrtS**2*(muu1**2*xd1 - (-1 + muu3**2 + muu4**2 - tred)*xd2 + muu2**2*xd3))))

     
      vbCLp2 = vbCLp2 + colori*0d0
	  vbCRp2 = vbCRp2 + colori*0d0
	  vbDp1 = vbDp1 + colori*0d0
	  vbDk1 = vBDk1 + colori*0d0
	  vbDk2 = vbDk2 + colori*0d0
		enddo
	      enddo
	    enddo
	  enddo
! ! 	  
! ! ! 
! ! ! 
! ! ! ! !-------------------------------------------------------------------------------------------------
! ! ! ! ! !--------------------------------------u-channels-------------------------------------------------
! ! ! ! ! !-------------------------------------------------------------------------------------------------
! ! ! ! ! !Couplings has to be exchanged to t-channels with exchanging (isf1<->isf2,itt1<->itt2,igen1<->igen2)
! ! ! ! ! ! ! ****************************************************
! ! ! ! ! ! ! ************ Box 1 Gluon-Squark (i,k,m)- Gluino- Squark (j,l,n) ***********************
! ! ! ! ! ! ! ****************************************************
! ! ! ! ! !Coupling 0 is L R exchanged (GluLt2->GluLt1) to sqsq*->qqbar case, calculation in mathematicafile, Minussign of charge flow in g1
! ! ! ! ! !g1(2->1) and g2(1->2) exchange incoming sfermion index
! ! ! ! ! !myDget muu1->muu2, muu2->muu1, tred->ured
	  do i = 1,2
	   do j = 1,2
	    do k = 3,4
	     do l = 3,4
          do m = 3,3
           do n = 3,3
	  
		  
		  Ml0 = 0d0
		  Ml1 = MSf(i,k,m)
		  Ml2 = MGl
		  Ml3 = MSf(j,l,n)
		  
		  g0L = GluLt1(j,l,ftt2,n,fgen2)
		  g0R = GluRt1(j,l,ftt2,n,fgen2)
		  g1 = VSfSf(1,isf1,j,itt1,l,igen1,n)
		  g2 = VSfSf(1,i,isf2,k,itt2,m,igen2)
		  g3L = GluLt1(i,k,ftt1,m,fgen1) 
		  g3R = GluRt1(i,k,ftt1,m,fgen1)
		  
		  call myCget(muu3**2*sqrtS**2,muu4**2*sqrtS**2,sqrtS**2,Ml1**2,Ml2**2,Ml3**2)
		  call myDget(muu2**2*sqrtS**2,muu3**2*sqrtS**2,muu4**2*sqrtS**2,muu1**2*sqrtS**2,ured*sqrtS**2,sqrtS**2, Ml0**2,Ml1**2,Ml2**2,Ml3**2)

  
      colori = CF/2d0*(Nc-1d0)
	
	  vbA = vbA + colori*0d0
	  
	    vbBL = vbBL + colori*(g1*g2*(g0L*g3L*Ml2*(xc0 + 
     -       sqrtS**2*(2*(-1 + muu1**2 + muu2**2)*xd0 + 
     -          (-1 + muu1**2 + 3*muu2**2)*xd1 - 2*xd2 + 
     -          3*muu1**2*xd2 + 3*muu2**2*xd2 + muu3**2*xd2 + 
     -          muu4**2*xd2 - 2*tred*xd2 - xd3 + 3*muu1**2*xd3 + 
     -          muu2**2*xd3)) + 
     -    g0L*g3R*muu3*sqrtS*
     -     (xc0 + xc1 + xc2 - 2*sqrtS**2*xd0 + 
     -       2*muu1**2*sqrtS**2*xd0 + 2*muu2**2*sqrtS**2*xd0 + 
     -       2*xd00 - sqrtS**2*xd1 + muu1**2*sqrtS**2*xd1 + 
     -       3*muu2**2*sqrtS**2*xd1 - sqrtS**2*xd12 + 
     -       muu1**2*sqrtS**2*xd12 + 3*muu2**2*sqrtS**2*xd12 - 
     -       sqrtS**2*xd13 + muu1**2*sqrtS**2*xd13 + 
     -       3*muu2**2*sqrtS**2*xd13 - 4*sqrtS**2*xd2 + 
     -       5*muu1**2*sqrtS**2*xd2 + 5*muu2**2*sqrtS**2*xd2 + 
     -       muu3**2*sqrtS**2*xd2 + muu4**2*sqrtS**2*xd2 - 
     -       2*sqrtS**2*tred*xd2 - 2*sqrtS**2*xd22 + 
     -       3*muu1**2*sqrtS**2*xd22 + 3*muu2**2*sqrtS**2*xd22 + 
     -       muu3**2*sqrtS**2*xd22 + muu4**2*sqrtS**2*xd22 - 
     -       2*sqrtS**2*tred*xd22 - 3*sqrtS**2*xd23 + 
     -       6*muu1**2*sqrtS**2*xd23 + 4*muu2**2*sqrtS**2*xd23 + 
     -       muu3**2*sqrtS**2*xd23 + muu4**2*sqrtS**2*xd23 - 
     -       2*sqrtS**2*tred*xd23 - 3*sqrtS**2*xd3 + 
     -       5*muu1**2*sqrtS**2*xd3 + 3*muu2**2*sqrtS**2*xd3 - 
     -       sqrtS**2*xd33 + 3*muu1**2*sqrtS**2*xd33 + 
     -       muu2**2*sqrtS**2*xd33) - 
     -    g0R*g3L*muu4*sqrtS*
     -     (xc2 + 2*xd00 + sqrtS**2*
     -        ((-1 + muu1**2 + 3*muu2**2)*xd13 + 
     -          (-2 + 3*muu1**2 + 3*muu2**2 + muu3**2 + muu4**2 - 
     -             2*tred)*xd23 - 2*xd3 + 2*muu1**2*xd3 + 
     -          2*muu2**2*xd3 - xd33 + 3*muu1**2*xd33 + muu2**2*xd33))))
     
     
	    vbBR = vbBR + colori*(g1*g2*(g0R*g3R*Ml2*(xc0 + 
     -       sqrtS**2*(2*(-1 + muu1**2 + muu2**2)*xd0 + 
     -          (-1 + muu1**2 + 3*muu2**2)*xd1 - 2*xd2 + 
     -          3*muu1**2*xd2 + 3*muu2**2*xd2 + muu3**2*xd2 + 
     -          muu4**2*xd2 - 2*tred*xd2 - xd3 + 3*muu1**2*xd3 + 
     -          muu2**2*xd3)) + 
     -    g0R*g3L*muu3*sqrtS*
     -     (xc0 + xc1 + xc2 - 2*sqrtS**2*xd0 + 
     -       2*muu1**2*sqrtS**2*xd0 + 2*muu2**2*sqrtS**2*xd0 + 
     -       2*xd00 - sqrtS**2*xd1 + muu1**2*sqrtS**2*xd1 + 
     -       3*muu2**2*sqrtS**2*xd1 - sqrtS**2*xd12 + 
     -       muu1**2*sqrtS**2*xd12 + 3*muu2**2*sqrtS**2*xd12 - 
     -       sqrtS**2*xd13 + muu1**2*sqrtS**2*xd13 + 
     -       3*muu2**2*sqrtS**2*xd13 - 4*sqrtS**2*xd2 + 
     -       5*muu1**2*sqrtS**2*xd2 + 5*muu2**2*sqrtS**2*xd2 + 
     -       muu3**2*sqrtS**2*xd2 + muu4**2*sqrtS**2*xd2 - 
     -       2*sqrtS**2*tred*xd2 - 2*sqrtS**2*xd22 + 
     -       3*muu1**2*sqrtS**2*xd22 + 3*muu2**2*sqrtS**2*xd22 + 
     -       muu3**2*sqrtS**2*xd22 + muu4**2*sqrtS**2*xd22 - 
     -       2*sqrtS**2*tred*xd22 - 3*sqrtS**2*xd23 + 
     -       6*muu1**2*sqrtS**2*xd23 + 4*muu2**2*sqrtS**2*xd23 + 
     -       muu3**2*sqrtS**2*xd23 + muu4**2*sqrtS**2*xd23 - 
     -       2*sqrtS**2*tred*xd23 - 3*sqrtS**2*xd3 + 
     -       5*muu1**2*sqrtS**2*xd3 + 3*muu2**2*sqrtS**2*xd3 - 
     -       sqrtS**2*xd33 + 3*muu1**2*sqrtS**2*xd33 + 
     -       muu2**2*sqrtS**2*xd33) - 
     -    g0L*g3R*muu4*sqrtS*
     -     (xc2 + 2*xd00 + sqrtS**2*
     -        ((-1 + muu1**2 + 3*muu2**2)*xd13 + 
     -          (-2 + 3*muu1**2 + 3*muu2**2 + muu3**2 + muu4**2 - 
     -             2*tred)*xd23 - 2*xd3 + 2*muu1**2*xd3 + 
     -          2*muu2**2*xd3 - xd33 + 3*muu1**2*xd33 + muu2**2*xd33))))
      
	    vbCLp2 = vbCLp2 + colori*( -(g0L*g1*g2*g3R*(xc0 + 4*xd00 + 
     -      sqrtS**2*(2*(-1 + muu1**2 + muu2**2)*xd0 + 
     -         (-3 + 3*muu1**2 + 5*muu2**2)*xd1 - xd11 + 
     -         muu1**2*xd11 + 3*muu2**2*xd11 - 3*xd12 + 
     -         4*muu1**2*xd12 + 6*muu2**2*xd12 + muu3**2*xd12 + 
     -         muu4**2*xd12 - 2*tred*xd12 - 2*xd13 + 4*muu1**2*xd13 + 
     -         4*muu2**2*xd13 - 4*xd2 + 5*muu1**2*xd2 + 
     -         5*muu2**2*xd2 + muu3**2*xd2 + muu4**2*xd2 - 
     -         2*tred*xd2 - 2*xd22 + 3*muu1**2*xd22 + 
     -         3*muu2**2*xd22 + muu3**2*xd22 + muu4**2*xd22 - 
     -         2*tred*xd22 - 3*xd23 + 6*muu1**2*xd23 + 
     -         4*muu2**2*xd23 + muu3**2*xd23 + muu4**2*xd23 - 
     -         2*tred*xd23 - 3*xd3 + 5*muu1**2*xd3 + 3*muu2**2*xd3 - 
     -         xd33 + 3*muu1**2*xd33 + muu2**2*xd33))))
   
	    vbCRp2 = vbCRp2 + colori*(-(g0R*g1*g2*g3L*(xc0 + 4*xd00 + 
     -      sqrtS**2*(2*(-1 + muu1**2 + muu2**2)*xd0 + 
     -         (-3 + 3*muu1**2 + 5*muu2**2)*xd1 - xd11 + 
     -         muu1**2*xd11 + 3*muu2**2*xd11 - 3*xd12 + 
     -         4*muu1**2*xd12 + 6*muu2**2*xd12 + muu3**2*xd12 + 
     -         muu4**2*xd12 - 2*tred*xd12 - 2*xd13 + 4*muu1**2*xd13 + 
     -         4*muu2**2*xd13 - 4*xd2 + 5*muu1**2*xd2 + 
     -         5*muu2**2*xd2 + muu3**2*xd2 + muu4**2*xd2 - 
     -         2*tred*xd2 - 2*xd22 + 3*muu1**2*xd22 + 
     -         3*muu2**2*xd22 + muu3**2*xd22 + muu4**2*xd22 - 
     -         2*tred*xd22 - 3*xd23 + 6*muu1**2*xd23 + 
     -         4*muu2**2*xd23 + muu3**2*xd23 + muu4**2*xd23 - 
     -         2*tred*xd23 - 3*xd3 + 5*muu1**2*xd3 + 3*muu2**2*xd3 - 
     -         xd33 + 3*muu1**2*xd33 + muu2**2*xd33))))
     
      vbCLp1 = vbCLp1 + colori*0d0
	  vbCRp1 = vbCRp1 + colori*0d0
	  vbDp1 = vbDp1 + colori*0d0
	  vbDk1 = vBDk1 + colori*0d0
	  vbDk2 = vbDk2 + colori*0d0
! 	
		   enddo
		  enddo
	     enddo
	    enddo
	   enddo
	  enddo
	  
! 	  
! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ****************************************************
! ! ! ! ! ! ! ! ! ! ************ Box 1.2 Gluon-Squark (i,k,m)- Gluino- Quark (l,n) ***********************
! ! ! ! ! ! ! ! ! ! ****************************************************
! ! ! ! ! ! !Coupling 1 is L R exchanged (GluLt2->GluLt1), Coupling 0 gains a minus sign (here minus times minus results in positiv coupling)
! ! ! ! ! ! !g1(2->1) and g2(1->2) exchange incoming sfermion index
! ! ! ! !myCget,myDget muu1->muu2, muu2->muu1, tred->ured, ured->tred
	  do i = 1,2
       do k = 3,4
	    do l = 3,4
	     do m = 3,3
		  do n = 3,3

		  Ml0 = 0d0
		  Ml1 = MSf(i,k,m)
		  Ml2 = MGl
		  Ml3 = Mf(l,n)
! 		  
		  g0L = GS*Krondelta(l,ftt2)*Krondelta(n,fgen2)
		  g0R = GS*Krondelta(l,ftt2)*Krondelta(n,fgen2)
		  g1L = GluLt1(isf1,itt1,l,igen1,n) 
		  g1R = GluRt1(isf1,itt1,l,igen1,n)
		  g2 = VSfSf(1,i,isf2,k,itt2,m,igen2)
		  g3L = GluLt1(i,k,ftt1,m,fgen1) 
		  g3R = GluRt1(i,k,ftt1,m,fgen1)
		  
! 		
		  call myCget(muu3**2*sqrtS**2,tred*sqrtS**2,muu1*sqrtS**2,Ml1**2,Ml2**2,Ml3**2)
		  call myDget(muu2**2*sqrtS**2,muu3**2*sqrtS**2,muu1**2*sqrtS**2,muu4**2*sqrtS**2,ured*sqrtS**2,tred*sqrtS**2, Ml0**2,Ml1**2,Ml2**2,Ml3**2)
		  		  
	colori = Nc*CF*(CF+1d0/2d0)

	vbA = vbA + colori*0d0
	
	vbBL = vbBL + colori*g2*(-(g0L*(g1R*g3R*Ml2*Ml3*muu3*sqrtS*(2*xd0 + xd1) + 
     -         g1L*g3R*muu3*sqrtS*(xc0 - xc1 - 2*xc2 + 2*Ml0**2*xd0 + 2*muu4**2*sqrtS**2*xd0 + Ml0**2*xd1 + 
     -            2*muu2**2*sqrtS**2*xd1 + 3*muu4**2*sqrtS**2*xd1 - 2*sqrtS**2*tred*xd1 + muu2**2*sqrtS**2*xd11 + 
     -            muu4**2*sqrtS**2*xd11 - sqrtS**2*tred*xd11 - sqrtS**2*xd12 + muu2**2*sqrtS**2*xd12 + 
     -            muu3**2*sqrtS**2*xd12 + 2*muu4**2*sqrtS**2*xd12 - sqrtS**2*tred*xd12 + 2*muu4**2*sqrtS**2*xd13 - 
     -            2*sqrtS**2*xd2 + 2*muu1**2*sqrtS**2*xd2 + 2*muu2**2*sqrtS**2*xd2 + 2*muu3**2*sqrtS**2*xd2 + 
     -            4*muu4**2*sqrtS**2*xd2 - 2*sqrtS**2*tred*xd2 + 4*muu4**2*sqrtS**2*xd3) - 
     -         g1R*g3L*Ml3*(xc0 + Ml0**2*xd0 + 
     -            sqrtS**2*(-xd2 + muu1**2*xd2 + muu4**2*xd2 - tred*xd2 - muu3**2*(2*xd0 + xd1 + xd2) + 
     -               muu4**2*xd3 - 2*tred*xd3 + muu2**2*(2*xd0 + 3*xd1 + 3*xd2 + 2*xd3))) - 
     -         g1L*g3L*Ml2*(xc0 + Ml0**2*xd0 + 
     -            sqrtS**2*(-xd2 - muu3**2*xd2 + muu4**2*xd2 + muu4**2*xd3 - tred*(2*xd0 + xd1 + xd2 + 2*xd3) + 
     -               muu2**2*(2*xd0 + 3*xd1 + 3*xd2 + 2*xd3))))) + 
     -    g0R*muu4*sqrtS*(g1L*g3L*Ml2*Ml3*(2*xd0 + xd1 + xd2 + xd3) + 
     -       g1R*g3R*Ml2*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3)) + 
     -       g1L*g3R*Ml3*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3)) + 
     -       g1R*g3L*(xc0 - xc2 + 2*Ml0**2*xd0 + 2*muu3**2*sqrtS**2*xd0 + 2*muu4**2*sqrtS**2*xd0 - 
     -          2*sqrtS**2*tred*xd0 + 2*xd00 + Ml0**2*xd1 + 2*muu2**2*sqrtS**2*xd1 + muu3**2*sqrtS**2*xd1 + 
     -          3*muu4**2*sqrtS**2*xd1 - 3*sqrtS**2*tred*xd1 + muu2**2*sqrtS**2*xd11 + muu4**2*sqrtS**2*xd11 - 
     -          sqrtS**2*tred*xd11 - sqrtS**2*xd12 + 2*muu2**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 + 
     -          3*muu4**2*sqrtS**2*xd12 - 2*sqrtS**2*tred*xd12 + muu2**2*sqrtS**2*xd13 + 3*muu4**2*sqrtS**2*xd13 - 
     -          sqrtS**2*tred*xd13 + Ml0**2*xd2 - 2*sqrtS**2*xd2 + muu1**2*sqrtS**2*xd2 + 2*muu2**2*sqrtS**2*xd2 + 
     -          4*muu3**2*sqrtS**2*xd2 + 5*muu4**2*sqrtS**2*xd2 - 4*sqrtS**2*tred*xd2 - sqrtS**2*xd22 + 
     -          muu2**2*sqrtS**2*xd22 + muu3**2*sqrtS**2*xd22 + 2*muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 - 
     -          sqrtS**2*xd23 + muu2**2*sqrtS**2*xd23 + muu3**2*sqrtS**2*xd23 + 4*muu4**2*sqrtS**2*xd23 - 
     -          sqrtS**2*tred*xd23 + Ml0**2*xd3 + 2*muu3**2*sqrtS**2*xd3 + 5*muu4**2*sqrtS**2*xd3 - 
     -          2*sqrtS**2*tred*xd3 + 2*muu4**2*sqrtS**2*xd33)))
     
     
	vbBR = vbBR + colori*g2*(-(g0R*(g1L*g3L*Ml2*Ml3*muu3*sqrtS*(2*xd0 + xd1) + 
     -         g1R*g3L*muu3*sqrtS*(xc0 - xc1 - 2*xc2 + 2*Ml0**2*xd0 + 2*muu4**2*sqrtS**2*xd0 + Ml0**2*xd1 + 
     -            2*muu2**2*sqrtS**2*xd1 + 3*muu4**2*sqrtS**2*xd1 - 2*sqrtS**2*tred*xd1 + muu2**2*sqrtS**2*xd11 + 
     -            muu4**2*sqrtS**2*xd11 - sqrtS**2*tred*xd11 - sqrtS**2*xd12 + muu2**2*sqrtS**2*xd12 + 
     -            muu3**2*sqrtS**2*xd12 + 2*muu4**2*sqrtS**2*xd12 - sqrtS**2*tred*xd12 + 2*muu4**2*sqrtS**2*xd13 - 
     -            2*sqrtS**2*xd2 + 2*muu1**2*sqrtS**2*xd2 + 2*muu2**2*sqrtS**2*xd2 + 2*muu3**2*sqrtS**2*xd2 + 
     -            4*muu4**2*sqrtS**2*xd2 - 2*sqrtS**2*tred*xd2 + 4*muu4**2*sqrtS**2*xd3) - 
     -         g1L*g3R*Ml3*(xc0 + Ml0**2*xd0 + 
     -            sqrtS**2*(-xd2 + muu1**2*xd2 + muu4**2*xd2 - tred*xd2 - muu3**2*(2*xd0 + xd1 + xd2) + 
     -               muu4**2*xd3 - 2*tred*xd3 + muu2**2*(2*xd0 + 3*xd1 + 3*xd2 + 2*xd3))) - 
     -         g1R*g3R*Ml2*(xc0 + Ml0**2*xd0 + 
     -            sqrtS**2*(-xd2 - muu3**2*xd2 + muu4**2*xd2 + muu4**2*xd3 - tred*(2*xd0 + xd1 + xd2 + 2*xd3) + 
     -               muu2**2*(2*xd0 + 3*xd1 + 3*xd2 + 2*xd3))))) + 
     -    g0L*muu4*sqrtS*(g1R*g3R*Ml2*Ml3*(2*xd0 + xd1 + xd2 + xd3) + 
     -       g1L*g3L*Ml2*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3)) + 
     -       g1R*g3L*Ml3*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3)) + 
     -       g1L*g3R*(xc0 - xc2 + 2*Ml0**2*xd0 + 2*muu3**2*sqrtS**2*xd0 + 2*muu4**2*sqrtS**2*xd0 - 
     -          2*sqrtS**2*tred*xd0 + 2*xd00 + Ml0**2*xd1 + 2*muu2**2*sqrtS**2*xd1 + muu3**2*sqrtS**2*xd1 + 
     -          3*muu4**2*sqrtS**2*xd1 - 3*sqrtS**2*tred*xd1 + muu2**2*sqrtS**2*xd11 + muu4**2*sqrtS**2*xd11 - 
     -          sqrtS**2*tred*xd11 - sqrtS**2*xd12 + 2*muu2**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 + 
     -          3*muu4**2*sqrtS**2*xd12 - 2*sqrtS**2*tred*xd12 + muu2**2*sqrtS**2*xd13 + 3*muu4**2*sqrtS**2*xd13 - 
     -          sqrtS**2*tred*xd13 + Ml0**2*xd2 - 2*sqrtS**2*xd2 + muu1**2*sqrtS**2*xd2 + 2*muu2**2*sqrtS**2*xd2 + 
     -          4*muu3**2*sqrtS**2*xd2 + 5*muu4**2*sqrtS**2*xd2 - 4*sqrtS**2*tred*xd2 - sqrtS**2*xd22 + 
     -          muu2**2*sqrtS**2*xd22 + muu3**2*sqrtS**2*xd22 + 2*muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 - 
     -          sqrtS**2*xd23 + muu2**2*sqrtS**2*xd23 + muu3**2*sqrtS**2*xd23 + 4*muu4**2*sqrtS**2*xd23 - 
     -          sqrtS**2*tred*xd23 + Ml0**2*xd3 + 2*muu3**2*sqrtS**2*xd3 + 5*muu4**2*sqrtS**2*xd3 - 
     -          2*sqrtS**2*tred*xd3 + 2*muu4**2*sqrtS**2*xd33)))
	
	
	vbCLp1 = vbCLp1 + colori*g2*(g0L*(2*g1L*g3L*Ml2*muu3*sqrtS*xd2 + g1R*g3R*Ml2*Ml3*(2*xd0 + xd1 + xd2) + 
     -       g1R*g3L*Ml3*muu3*sqrtS*(2*xd0 + xd1 + 2*xd2) + 
     -       g1L*g3R*(2*xc0 - xc2 + 3*Ml0**2*xd0 + 2*muu2**2*sqrtS**2*xd0 + 2*muu4**2*sqrtS**2*xd0 - 
     -          2*sqrtS**2*tred*xd0 + Ml0**2*xd1 + 5*muu2**2*sqrtS**2*xd1 + 3*muu4**2*sqrtS**2*xd1 - 
     -          3*sqrtS**2*tred*xd1 + muu2**2*sqrtS**2*xd11 + muu4**2*sqrtS**2*xd11 - sqrtS**2*tred*xd11 - 
     -          sqrtS**2*xd12 + 2*muu2**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 + 3*muu4**2*sqrtS**2*xd12 - 
     -          2*sqrtS**2*tred*xd12 + 2*muu4**2*sqrtS**2*xd13 + Ml0**2*xd2 - 3*sqrtS**2*xd2 + 
     -          2*muu1**2*sqrtS**2*xd2 + 5*muu2**2*sqrtS**2*xd2 + 3*muu3**2*sqrtS**2*xd2 + 6*muu4**2*sqrtS**2*xd2 - 
     -          5*sqrtS**2*tred*xd2 - sqrtS**2*xd22 + muu2**2*sqrtS**2*xd22 + muu3**2*sqrtS**2*xd22 + 
     -          2*muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 + 2*muu4**2*sqrtS**2*xd23 + 2*muu2**2*sqrtS**2*xd3 + 
     -          5*muu4**2*sqrtS**2*xd3 - 2*sqrtS**2*tred*xd3)) - 
     -    g0R*muu4*sqrtS*(g1L*g3R*Ml3*xd3 + g1R*g3R*Ml2*(2*xd0 + xd1 + xd2 + 2*xd3) + 
     -       g1R*g3L*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3))))
	
	
	vbCRp1 = vbCRp1 + colori*g2*(g0R*(2*g1R*g3R*Ml2*muu3*sqrtS*xd2 + g1L*g3L*Ml2*Ml3*(2*xd0 + xd1 + xd2) + 
     -       g1L*g3R*Ml3*muu3*sqrtS*(2*xd0 + xd1 + 2*xd2) + 
     -       g1R*g3L*(2*xc0 - xc2 + 3*Ml0**2*xd0 + 2*muu2**2*sqrtS**2*xd0 + 2*muu4**2*sqrtS**2*xd0 - 
     -          2*sqrtS**2*tred*xd0 + Ml0**2*xd1 + 5*muu2**2*sqrtS**2*xd1 + 3*muu4**2*sqrtS**2*xd1 - 
     -          3*sqrtS**2*tred*xd1 + muu2**2*sqrtS**2*xd11 + muu4**2*sqrtS**2*xd11 - sqrtS**2*tred*xd11 - 
     -          sqrtS**2*xd12 + 2*muu2**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 + 3*muu4**2*sqrtS**2*xd12 - 
     -          2*sqrtS**2*tred*xd12 + 2*muu4**2*sqrtS**2*xd13 + Ml0**2*xd2 - 3*sqrtS**2*xd2 + 
     -          2*muu1**2*sqrtS**2*xd2 + 5*muu2**2*sqrtS**2*xd2 + 3*muu3**2*sqrtS**2*xd2 + 6*muu4**2*sqrtS**2*xd2 - 
     -          5*sqrtS**2*tred*xd2 - sqrtS**2*xd22 + muu2**2*sqrtS**2*xd22 + muu3**2*sqrtS**2*xd22 + 
     -          2*muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 + 2*muu4**2*sqrtS**2*xd23 + 2*muu2**2*sqrtS**2*xd3 + 
     -          5*muu4**2*sqrtS**2*xd3 - 2*sqrtS**2*tred*xd3)) - 
     -    g0L*muu4*sqrtS*(g1R*g3L*Ml3*xd3 + g1L*g3L*Ml2*(2*xd0 + xd1 + xd2 + 2*xd3) + 
     -       g1L*g3R*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3))))
	
	
      vbCLp2 = vbCLp2 + colori*0d0
	  vbCRp2 = vbCRp2 + colori*0d0
	  vbDp1 = vbDp1 + colori*0d0
	  vbDk1 = vBDk1 + colori*0d0
	  vbDk2 = vbDk2 + colori*0d0
	
		  enddo
		 enddo
	    enddo
	   enddo
	  enddo

! ! ! ! ****************************************************
! ! ! ! ! ************ Box 1.3 Gluino-Quark (k,m)- Gluon- Squark (j,l,n) ***********************
! ! ! ! ! ****************************************************
! ! !Coupling 0 is L R exchanged (GluLt2->GluLt1), minus sign in coupling 1 from opposite charge flow
! ! !g1(2->1) and g2(1->2) exchange incoming sfermion index
! ! !myCget,myDget muu1->muu2, muu2->muu1, tred->ured, ured->tred
	do j = 1,2
	    do k = 3,4
	      do l = 3,4
		do m = 3,3
		  do n = 3,3
		  
		  Ml0 = 0d0
		  Ml1 = MSf(j,l,n)
		  Ml2 = MGl
		  Ml3 = Mf(k,m)
! 		  
		  
		  g0L = GluLt1(j,l,ftt2,n,fgen2) 
		  g0R = GluRt1(j,l,ftt2,n,fgen2)
		  g1 = -VSfSf(1,isf1,j,itt1,l,igen1,n)
		  g2L = GluLt1(isf2,itt2,k,igen2,m) 
		  g2R = GluRt1(isf2,itt2,k,igen2,m)
		  g3L = -GS*Krondelta(k,ftt1)*Krondelta(m,fgen1)
		  g3R = -GS*Krondelta(k,ftt1)*Krondelta(m,fgen1)
! 		
		  call myCget(muu3**2*sqrtS**2,tred*sqrtS**2,muu1*sqrtS**2,Ml1**2,Ml2**2,Ml3**2)
		  call myDget(muu1**2*sqrtS**2,muu4**2*sqrtS**2,muu2**2*sqrtS**2,muu3**2*sqrtS**2,ured*sqrtS**2,tred*sqrtS**2, Ml0**2,Ml1**2,Ml2**2,Ml3**2)
		  		  
	colori = Nc*CF*(CF+1d0/2d0)

	vbA = vbA + colori*0d0
	
	vbBL = vbBL + colori* g1*(g0R*muu4*sqrtS*(g2R*Ml2*(g3R*Ml3 - g3L*muu3*sqrtS)*xd2 + g2L*g3L*Ml3*muu3*sqrtS*xd3 + 
     -       g2L*g3R*(2*xc0 + xc1 + xc2 + Ml0**2*xd0 + 2*muu1**2*sqrtS**2*xd0 + 2*muu3**2*sqrtS**2*xd0 - 
     -          2*sqrtS**2*tred*xd0 - Ml0**2*xd1 + 3*muu1**2*sqrtS**2*xd1 + muu3**2*sqrtS**2*xd1 - 
     -          sqrtS**2*tred*xd1 + muu1**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 - sqrtS**2*tred*xd12 - 
     -          sqrtS**2*xd2 + 3*muu1**2*sqrtS**2*xd2 + 4*muu3**2*sqrtS**2*xd2 + muu4**2*sqrtS**2*xd2 - 
     -          3*sqrtS**2*tred*xd2 - sqrtS**2*xd22 + muu1**2*sqrtS**2*xd22 + 2*muu3**2*sqrtS**2*xd22 + 
     -          muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 + 2*muu3**2*sqrtS**2*xd23 + 2*muu1**2*sqrtS**2*xd3 + 
     -          3*muu3**2*sqrtS**2*xd3 - 2*sqrtS**2*tred*xd3)) + 
     -    g0L*(g2L*g3L*Ml2*Ml3*muu3*sqrtS*xd3 + 
     -       g2R*g3R*Ml3*(xc0 + Ml0**2*xd0 + sqrtS**2*
     -           ((-1 + muu2**2 + muu3**2 + muu4**2 - tred)*xd2 + 2*(muu3**2 - tred)*xd3 + 
     -             muu1**2*(2*xd0 + 3*xd1 + 3*xd2 + 2*xd3))) + 
     -       g2L*g3R*Ml2*(xc0 + Ml0**2*xd0 + sqrtS**2*
     -           (-2*tred*xd0 - tred*xd1 - xd2 + muu4**2*xd2 - tred*xd2 - 2*tred*xd3 + 
     -             muu1**2*(2*xd0 + 3*xd1 + 3*xd2 + 2*xd3) + muu3**2*(2*xd0 + xd1 + 2*xd2 + 3*xd3))) + 
     -       g2R*g3L*muu3*sqrtS*(-2*xc0 - 2*muu1**2*sqrtS**2*xd0 + 2*xd00 - 3*muu1**2*sqrtS**2*xd1 + 
     -          muu1**2*sqrtS**2*xd13 + muu3**2*sqrtS**2*xd13 - sqrtS**2*tred*xd13 + sqrtS**2*xd2 - 
     -          3*muu1**2*sqrtS**2*xd2 - muu2**2*sqrtS**2*xd2 - muu3**2*sqrtS**2*xd2 - muu4**2*sqrtS**2*xd2 + 
     -          sqrtS**2*tred*xd2 - sqrtS**2*xd23 + muu1**2*sqrtS**2*xd23 + 2*muu3**2*sqrtS**2*xd23 + 
     -          muu4**2*sqrtS**2*xd23 - sqrtS**2*tred*xd23 - 2*muu1**2*sqrtS**2*xd3 + 
     -          Ml0**2*(-xd0 + xd1 + xd2 + xd3) + 2*muu3**2*sqrtS**2*xd33)))
	
     
	vbBR = vbBR + colori*g1*(g0L*muu4*sqrtS*(g2L*Ml2*(g3L*Ml3 - g3R*muu3*sqrtS)*xd2 + g2R*g3R*Ml3*muu3*sqrtS*xd3 + 
     -       g2R*g3L*(2*xc0 + xc1 + xc2 + Ml0**2*xd0 + 2*muu1**2*sqrtS**2*xd0 + 2*muu3**2*sqrtS**2*xd0 - 
     -          2*sqrtS**2*tred*xd0 - Ml0**2*xd1 + 3*muu1**2*sqrtS**2*xd1 + muu3**2*sqrtS**2*xd1 - 
     -          sqrtS**2*tred*xd1 + muu1**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 - sqrtS**2*tred*xd12 - 
     -          sqrtS**2*xd2 + 3*muu1**2*sqrtS**2*xd2 + 4*muu3**2*sqrtS**2*xd2 + muu4**2*sqrtS**2*xd2 - 
     -          3*sqrtS**2*tred*xd2 - sqrtS**2*xd22 + muu1**2*sqrtS**2*xd22 + 2*muu3**2*sqrtS**2*xd22 + 
     -          muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 + 2*muu3**2*sqrtS**2*xd23 + 2*muu1**2*sqrtS**2*xd3 + 
     -          3*muu3**2*sqrtS**2*xd3 - 2*sqrtS**2*tred*xd3)) + 
     -    g0R*(g2R*g3R*Ml2*Ml3*muu3*sqrtS*xd3 + 
     -       g2L*g3L*Ml3*(xc0 + Ml0**2*xd0 + sqrtS**2*
     -           ((-1 + muu2**2 + muu3**2 + muu4**2 - tred)*xd2 + 2*(muu3**2 - tred)*xd3 + 
     -             muu1**2*(2*xd0 + 3*xd1 + 3*xd2 + 2*xd3))) + 
     -       g2R*g3L*Ml2*(xc0 + Ml0**2*xd0 + sqrtS**2*
     -           (-2*tred*xd0 - tred*xd1 - xd2 + muu4**2*xd2 - tred*xd2 - 2*tred*xd3 + 
     -             muu1**2*(2*xd0 + 3*xd1 + 3*xd2 + 2*xd3) + muu3**2*(2*xd0 + xd1 + 2*xd2 + 3*xd3))) + 
     -       g2L*g3R*muu3*sqrtS*(-2*xc0 - 2*muu1**2*sqrtS**2*xd0 + 2*xd00 - 3*muu1**2*sqrtS**2*xd1 + 
     -          muu1**2*sqrtS**2*xd13 + muu3**2*sqrtS**2*xd13 - sqrtS**2*tred*xd13 + sqrtS**2*xd2 - 
     -          3*muu1**2*sqrtS**2*xd2 - muu2**2*sqrtS**2*xd2 - muu3**2*sqrtS**2*xd2 - muu4**2*sqrtS**2*xd2 + 
     -          sqrtS**2*tred*xd2 - sqrtS**2*xd23 + muu1**2*sqrtS**2*xd23 + 2*muu3**2*sqrtS**2*xd23 + 
     -          muu4**2*sqrtS**2*xd23 - sqrtS**2*tred*xd23 - 2*muu1**2*sqrtS**2*xd3 + 
     -          Ml0**2*(-xd0 + xd1 + xd2 + xd3) + 2*muu3**2*sqrtS**2*xd33)))
     
	vbCLp1 = vbCLp1 + colori* g1*(g0L*(g2L*g3L*Ml2*Ml3*(2*xd0 + xd1 + xd2) - g2R*g3R*Ml3*muu3*sqrtS*xd3 - 
     -       g2L*g3R*Ml2*muu3*sqrtS*(2*xd0 + xd1 + xd2 + 2*xd3) + 
     -       g2R*g3L*(4*xc0 - xc2 + 3*Ml0**2*xd0 + 2*muu1**2*sqrtS**2*xd0 + 2*muu3**2*sqrtS**2*xd0 - 
     -          2*sqrtS**2*tred*xd0 - Ml0**2*xd1 + 5*muu1**2*sqrtS**2*xd1 + 3*muu3**2*sqrtS**2*xd1 - 
     -          3*sqrtS**2*tred*xd1 + muu1**2*sqrtS**2*xd11 + muu3**2*sqrtS**2*xd11 - sqrtS**2*tred*xd11 - 
     -          sqrtS**2*xd12 + 2*muu1**2*sqrtS**2*xd12 + 3*muu3**2*sqrtS**2*xd12 + muu4**2*sqrtS**2*xd12 - 
     -          2*sqrtS**2*tred*xd12 + 2*muu3**2*sqrtS**2*xd13 - Ml0**2*xd2 - 3*sqrtS**2*xd2 + 
     -          5*muu1**2*sqrtS**2*xd2 + 2*muu2**2*sqrtS**2*xd2 + 6*muu3**2*sqrtS**2*xd2 + 3*muu4**2*sqrtS**2*xd2 - 
     -          5*sqrtS**2*tred*xd2 - sqrtS**2*xd22 + muu1**2*sqrtS**2*xd22 + 2*muu3**2*sqrtS**2*xd22 + 
     -          muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 + 2*muu3**2*sqrtS**2*xd23 + 2*muu1**2*sqrtS**2*xd3 + 
     -          5*muu3**2*sqrtS**2*xd3 - 2*sqrtS**2*tred*xd3)) + 
     -    g0R*muu4*sqrtS*(2*g2R*g3L*Ml2*xd2 + g2L*g3L*Ml3*(2*xd0 + xd1 + 2*xd2) - 
     -       g2L*g3R*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3))))
     
	vbCRp1 = vbCRp1 + colori*g1*(g0R*(g2R*g3R*Ml2*Ml3*(2*xd0 + xd1 + xd2) - g2L*g3L*Ml3*muu3*sqrtS*xd3 - 
     -       g2R*g3L*Ml2*muu3*sqrtS*(2*xd0 + xd1 + xd2 + 2*xd3) + 
     -       g2L*g3R*(4*xc0 - xc2 + 3*Ml0**2*xd0 + 2*muu1**2*sqrtS**2*xd0 + 2*muu3**2*sqrtS**2*xd0 - 
     -          2*sqrtS**2*tred*xd0 - Ml0**2*xd1 + 5*muu1**2*sqrtS**2*xd1 + 3*muu3**2*sqrtS**2*xd1 - 
     -          3*sqrtS**2*tred*xd1 + muu1**2*sqrtS**2*xd11 + muu3**2*sqrtS**2*xd11 - sqrtS**2*tred*xd11 - 
     -          sqrtS**2*xd12 + 2*muu1**2*sqrtS**2*xd12 + 3*muu3**2*sqrtS**2*xd12 + muu4**2*sqrtS**2*xd12 - 
     -          2*sqrtS**2*tred*xd12 + 2*muu3**2*sqrtS**2*xd13 - Ml0**2*xd2 - 3*sqrtS**2*xd2 + 
     -          5*muu1**2*sqrtS**2*xd2 + 2*muu2**2*sqrtS**2*xd2 + 6*muu3**2*sqrtS**2*xd2 + 3*muu4**2*sqrtS**2*xd2 - 
     -          5*sqrtS**2*tred*xd2 - sqrtS**2*xd22 + muu1**2*sqrtS**2*xd22 + 2*muu3**2*sqrtS**2*xd22 + 
     -          muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 + 2*muu3**2*sqrtS**2*xd23 + 2*muu1**2*sqrtS**2*xd3 + 
     -          5*muu3**2*sqrtS**2*xd3 - 2*sqrtS**2*tred*xd3)) + 
     -    g0L*muu4*sqrtS*(2*g2L*g3R*Ml2*xd2 + g2R*g3R*Ml3*(2*xd0 + xd1 + 2*xd2) - 
     -       g2R*g3L*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3))))
     
!      
      vbCLp2 = vbCLp2 + colori*0d0
	  vbCRp2 = vbCRp2 + colori*0d0
	  vbDp1 = vbDp1 + colori*0d0
	  vbDk1 = vBDk1 + colori*0d0
	  vbDk2 = vbDk2 + colori*0d0
	
		  enddo
		enddo
	      enddo
	    enddo
	enddo
! 	
! ! 
! ! ! ****************************************************
! ! ! ************ Box 4 Gluino-Quark(k,m)-Gluon- Quark(l,n) ***********************
! ! ! ****************************************************
! ! !Coupling 1 is L R exchanged (GluLt2->GluLt1), Coupling 0 gains a minus sign (here minus times minus results in positiv coupling)
! ! !g1(2->1) and g2(1->2) exchange incoming sfermion index
! ! !myCget,myDget muu1->muu2, muu2->muu1, tred->ured, ured->tred
! ! ! 
	  do k = 3,4
       do l = 3,4
        do m = 3,3
		 do n = 3,3
		
		  Ml0 = MGl
		  Ml1 = Mf(k,m)
		  Ml2 = 0
		  Ml3 = Mf(l,n)
		  
		  g0 = GS*Krondelta(l,ftt2)*Krondelta(n,fgen2)
		  g1L = GluLt1(isf1,itt1,l,igen1,n)
		  g1R = GluRt1(isf1,itt1,l,igen1,n)
		  g2L = GluLt1(isf2,itt2,k,igen2,m)
		  g2R = GluRt1(isf2,itt2,k,igen2,m)
		  g3 = -GS*Krondelta(k,ftt1)*Krondelta(m,fgen1)
		  
		  call myCget(muu3**2*sqrtS**2,muu4**2*sqrtS**2,sqrtS**2,Ml1**2,Ml2**2,Ml3**2)
		  call myDget(muu2**2*sqrtS**2,muu3**2*sqrtS**2,muu4**2*sqrtS**2,muu1**2*sqrtS**2,ured*sqrtS**2,sqrtS**2, Ml0**2,Ml1**2,Ml2**2,Ml3**2)
		  		  
	colori = CF/2d0*(Nc-1)

	vbA = vbA + colori*0d0
	
	vbBL = vbBL + colori*(-2*g0*g3*(g1R*g2R*Ml0*(2*Ml1*Ml3*xd0 - Ml3*muu3*sqrtS*(xd0 + xd1) + Ml1*muu4*sqrtS*(xd1 + xd2)) - 
     -    g1L*g2R*muu3*sqrtS*(-xc1 - xc2 + Ml0**2*xd0 + Ml0**2*xd1 + muu1**2*sqrtS**2*xd2 + muu1**2*sqrtS**2*xd3) + 
     -    g1L*g2R*Ml1*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*
     -        (Ml3*muu4*(xd1 + xd2) + sqrtS*((-1 + muu1**2 + muu2**2)*xd1 + 
     -             (-1 + 2*muu1**2 + muu2**2 + muu3**2 - tred)*xd2 + 2*muu1**2*xd3))) + 
     -    g1L*g2L*Ml0*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*
     -        (-(Ml1*muu3*xd1) + Ml3*muu4*(xd0 + xd1 + xd2) + 
     -          sqrtS*((-1 + muu1**2 + muu2**2)*xd0 + (-1 + muu1**2 + 3*muu2**2)*xd1 - 2*xd2 + 3*muu1**2*xd2 + 
     -             3*muu2**2*xd2 + muu3**2*xd2 + muu4**2*xd2 - 2*tred*xd2 - xd3 + 3*muu1**2*xd3 + muu2**2*xd3))) + 
     -    g1R*g2L*(muu4*sqrtS*(-xc2 + Ml0**2*(xd0 + xd1 + xd2) + 
     -          sqrtS**2*((muu3**2 - tred)*xd2 + muu1**2*(xd2 + xd3))) + 
     -       Ml3*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*(-(Ml1*muu3*xd1) + 
     -          sqrtS*((-1 + muu1**2 + muu4**2 - tred)*xd2 + (-1 + muu1**2)*xd3 + muu2**2*(2*xd1 + 2*xd2 + xd3)))))))
          
	vbBR = vbBR + colori*(-2*g0*g3*(g1L*g2L*Ml0*(2*Ml1*Ml3*xd0 - Ml3*muu3*sqrtS*(xd0 + xd1) + Ml1*muu4*sqrtS*(xd1 + xd2)) - 
     -    g1R*g2L*muu3*sqrtS*(-xc1 - xc2 + Ml0**2*xd0 + Ml0**2*xd1 + muu1**2*sqrtS**2*xd2 + muu1**2*sqrtS**2*xd3) + 
     -    g1R*g2L*Ml1*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*
     -        (Ml3*muu4*(xd1 + xd2) + sqrtS*((-1 + muu1**2 + muu2**2)*xd1 + 
     -             (-1 + 2*muu1**2 + muu2**2 + muu3**2 - tred)*xd2 + 2*muu1**2*xd3))) + 
     -    g1R*g2R*Ml0*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*
     -        (-(Ml1*muu3*xd1) + Ml3*muu4*(xd0 + xd1 + xd2) + 
     -          sqrtS*((-1 + muu1**2 + muu2**2)*xd0 + (-1 + muu1**2 + 3*muu2**2)*xd1 - 2*xd2 + 3*muu1**2*xd2 + 
     -             3*muu2**2*xd2 + muu3**2*xd2 + muu4**2*xd2 - 2*tred*xd2 - xd3 + 3*muu1**2*xd3 + muu2**2*xd3))) + 
     -    g1L*g2R*(muu4*sqrtS*(-xc2 + Ml0**2*(xd0 + xd1 + xd2) + 
     -          sqrtS**2*((muu3**2 - tred)*xd2 + muu1**2*(xd2 + xd3))) + 
     -       Ml3*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*(-(Ml1*muu3*xd1) + 
     -             sqrtS*((-1 + muu1**2 + muu4**2 - tred)*xd2 + (-1 + muu1**2)*xd3 + muu2**2*(2*xd1 + 2*xd2 + xd3)))))))
	
	vbCLp1 = vbCLp1 + colori*(-2*g0*g3*(-(g1R*g2L*muu3*muu4*sqrtS**2*xd2) + g1R*g2L*Ml1*Ml3*(xd1 + xd2 + xd3) + 
     -    g1L*g2L*Ml0*Ml1*(xd0 + xd1 + xd2 + xd3) + g1R*g2R*Ml0*Ml3*(xd0 + xd1 + xd2 + xd3) + 
     -    g1L*g2R*(xc0 + Ml0**2*(2*xd0 + xd1 + xd2 + xd3) + 
     -       sqrtS**2*(-(tred*xd2) + muu2**2*(xd1 + xd2) + muu1**2*(xd2 + xd3)))))
	
	vbCRp1 = vbCRp1 + colori*(-2*g0*g3*(-(g1L*g2R*muu3*muu4*sqrtS**2*xd2) + g1L*g2R*Ml1*Ml3*(xd1 + xd2 + xd3) + 
     -    g1R*g2R*Ml0*Ml1*(xd0 + xd1 + xd2 + xd3) + g1L*g2L*Ml0*Ml3*(xd0 + xd1 + xd2 + xd3) + 
     -    g1R*g2L*(xc0 + Ml0**2*(2*xd0 + xd1 + xd2 + xd3) + 
     -       sqrtS**2*(-(tred*xd2) + muu2**2*(xd1 + xd2) + muu1**2*(xd2 + xd3))))) 

     
	  vbCLp2 = vbCLp2 + colori*0d0
	  vbCRp2 = vbCRp2 + colori*0d0
	  vbDp1 = vbDp1 + colori*0d0
	  vbDk1 = vBDk1 + colori*0d0
	  vbDk2 = vbDk2 + colori*0d0
		enddo
	      enddo
	    enddo
	  enddo
! ! ! !-------------------------------------------------------------------------------------------------
! ! ! !-------------------------------------------------------------------------------------------------
! ! ! !-------------------------------------------------------------------------------------------------
! ! ! !-------------------------------------------------------------------------------------------------
! ! ! !-------------------------------------------------------------------------------------------------
! ! ! !-------------------------------------------------------------------------------------------------
! ! ! !-----------------here boxdiagramms with neutralino-------------------------------------------
! ! ! !-------------------------------------------------------------------------------------------------
! ! ! !-------------------------------------------------------------------------------------------------
! ! ! !-------------------------------------------------------------------------------------------------
! ! ! !-------------------------------------------------------------------------------------------------
! ! ! !-------------------------------------------------------------------------------------------------
! ! ! !-------------------------------------------------------------------------------------------------
! ! 
! ! 
! ! 
! ! ! -------------------------------------------------------------------------------------------------
! ! ! --------------------------------------t-channels-------------------------------------------------
! ! ! -------------------------------------------------------------------------------------------------
! ! 
	do b = 1,6 	!Start Neutralino loop
! ****************************************************
! ************ Box 1 Gluon-Squark (i,k,m)- Neutralino- Squark (j,l,n) ***********************
! ****************************************************
! Coupling 0 L<->R , Coupling 1 gains minus sign due to charge flow change

	do i = 1,2
	  do j = 1,2
	    do k = 3,4
	      do l = 3,4
		do m = 3,3
		  do n = 3,3
		  	
		  
		  Ml0 = 0d0
		  Ml1 = MSf(i,k,m)
		 
		  if(b.LE.4)then
		    Ml2 = MNeu(b)
		  else
		    if(b.EQ.5) then
		      Ml2 = MCha(1)
		    else
		      Ml2 = MCha(2)
		    endif
		  endif
		  
		  Ml3 = MSf(j,l,n)
		  
		  
		  g0L = AL(b,j,l,ftt2,n,fgen2)
		  g0R = AR(b,j,l,ftt2,n,fgen2)
		  g1 = VSfSf(1,isf2,j,itt2,l,igen2,n)
		  g2 = VSfSf(1,i,isf1,k,itt1,m,igen1)
		  g3L = AL(b,i,k,ftt1,m,fgen1) 
		  g3R = AR(b,i,k,ftt1,m,fgen1)
		  
		  call myCget(muu3**2*sqrtS**2,muu4**2*sqrtS**2,sqrtS**2,Ml1**2,Ml2**2,Ml3**2)
		  call myDget(muu1**2*sqrtS**2,muu3**2*sqrtS**2,muu4**2*sqrtS**2,muu2**2*sqrtS**2,tred*sqrtS**2,sqrtS**2, Ml0**2,Ml1**2,Ml2**2,Ml3**2)

	  colori = Nc*CF
 	
	  vbA = vbA + colori*0d0
	  
	  vbBR = vbBR + colori*g1*g2*(g0R*g3R*Ml2*(xc0 + 
     -       sqrtS**2*(2*(-1 + muu1**2 + muu2**2)*xd0 + 
     -          (-1 + 3*muu1**2 + muu2**2)*xd1 + muu1**2*xd2 + 
     -          muu2**2*xd2 - muu3**2*xd2 - muu4**2*xd2 + 2*tred*xd2 - 
     -          xd3 + muu1**2*xd3 + 3*muu2**2*xd3)) + 
     -    g0R*g3L*muu3*sqrtS*
     -     (xc0 + xc1 + xc2 - 2*sqrtS**2*xd0 + 
     -       2*muu1**2*sqrtS**2*xd0 + 2*muu2**2*sqrtS**2*xd0 + 
     -       2*xd00 - sqrtS**2*xd1 + 3*muu1**2*sqrtS**2*xd1 + 
     -       muu2**2*sqrtS**2*xd1 - sqrtS**2*xd12 + 
     -       3*muu1**2*sqrtS**2*xd12 + muu2**2*sqrtS**2*xd12 - 
     -       sqrtS**2*xd13 + 3*muu1**2*sqrtS**2*xd13 + 
     -       muu2**2*sqrtS**2*xd13 - 2*sqrtS**2*xd2 + 
     -       3*muu1**2*sqrtS**2*xd2 + 3*muu2**2*sqrtS**2*xd2 - 
     -       muu3**2*sqrtS**2*xd2 - muu4**2*sqrtS**2*xd2 + 
     -       2*sqrtS**2*tred*xd2 + muu1**2*sqrtS**2*xd22 + 
     -       muu2**2*sqrtS**2*xd22 - muu3**2*sqrtS**2*xd22 - 
     -       muu4**2*sqrtS**2*xd22 + 2*sqrtS**2*tred*xd22 - 
     -       sqrtS**2*xd23 + 2*muu1**2*sqrtS**2*xd23 + 
     -       4*muu2**2*sqrtS**2*xd23 - muu3**2*sqrtS**2*xd23 - 
     -       muu4**2*sqrtS**2*xd23 + 2*sqrtS**2*tred*xd23 - 
     -       3*sqrtS**2*xd3 + 3*muu1**2*sqrtS**2*xd3 + 
     -       5*muu2**2*sqrtS**2*xd3 - sqrtS**2*xd33 + 
     -       muu1**2*sqrtS**2*xd33 + 3*muu2**2*sqrtS**2*xd33) + 
     -    g0L*g3R*muu4*sqrtS*
     -     (-xc2 - 2*xd00 - sqrtS**2*
     -        ((-1 + 3*muu1**2 + muu2**2)*xd13 + muu2**2*xd23 - 
     -          muu3**2*xd23 - muu4**2*xd23 + 2*tred*xd23 - 2*xd3 + 
     -          2*muu2**2*xd3 - xd33 + 3*muu2**2*xd33 + 
     -          muu1**2*(xd23 + 2*xd3 + xd33))))
	  
	  vbBL = vbBL + colori*g1*g2*(g0L*g3L*Ml2*(xc0 + 
     -       sqrtS**2*(2*(-1 + muu1**2 + muu2**2)*xd0 + 
     -          (-1 + 3*muu1**2 + muu2**2)*xd1 + muu1**2*xd2 + 
     -          muu2**2*xd2 - muu3**2*xd2 - muu4**2*xd2 + 2*tred*xd2 - 
     -          xd3 + muu1**2*xd3 + 3*muu2**2*xd3)) + 
     -    g0L*g3R*muu3*sqrtS*
     -     (xc0 + xc1 + xc2 - 2*sqrtS**2*xd0 + 
     -       2*muu1**2*sqrtS**2*xd0 + 2*muu2**2*sqrtS**2*xd0 + 
     -       2*xd00 - sqrtS**2*xd1 + 3*muu1**2*sqrtS**2*xd1 + 
     -       muu2**2*sqrtS**2*xd1 - sqrtS**2*xd12 + 
     -       3*muu1**2*sqrtS**2*xd12 + muu2**2*sqrtS**2*xd12 - 
     -       sqrtS**2*xd13 + 3*muu1**2*sqrtS**2*xd13 + 
     -       muu2**2*sqrtS**2*xd13 - 2*sqrtS**2*xd2 + 
     -       3*muu1**2*sqrtS**2*xd2 + 3*muu2**2*sqrtS**2*xd2 - 
     -       muu3**2*sqrtS**2*xd2 - muu4**2*sqrtS**2*xd2 + 
     -       2*sqrtS**2*tred*xd2 + muu1**2*sqrtS**2*xd22 + 
     -       muu2**2*sqrtS**2*xd22 - muu3**2*sqrtS**2*xd22 - 
     -       muu4**2*sqrtS**2*xd22 + 2*sqrtS**2*tred*xd22 - 
     -       sqrtS**2*xd23 + 2*muu1**2*sqrtS**2*xd23 + 
     -       4*muu2**2*sqrtS**2*xd23 - muu3**2*sqrtS**2*xd23 - 
     -       muu4**2*sqrtS**2*xd23 + 2*sqrtS**2*tred*xd23 - 
     -       3*sqrtS**2*xd3 + 3*muu1**2*sqrtS**2*xd3 + 
     -       5*muu2**2*sqrtS**2*xd3 - sqrtS**2*xd33 + 
     -       muu1**2*sqrtS**2*xd33 + 3*muu2**2*sqrtS**2*xd33) + 
     -    g0R*g3L*muu4*sqrtS*
     -     (-xc2 - 2*xd00 - sqrtS**2*
     -        ((-1 + 3*muu1**2 + muu2**2)*xd13 + muu2**2*xd23 - 
     -          muu3**2*xd23 - muu4**2*xd23 + 2*tred*xd23 - 2*xd3 + 
     -          2*muu2**2*xd3 - xd33 + 3*muu2**2*xd33 + 
     -          muu1**2*(xd23 + 2*xd3 + xd33))))
	  
	  vbCLp1 = vbCLp1 + colori*(-(g0L*g1*g2*g3R*(xc0 + 4*xd00 + 
     -      sqrtS**2*(2*(-1 + muu1**2 + muu2**2)*xd0 + 
     -         (-3 + 5*muu1**2 + 3*muu2**2)*xd1 - xd11 + 
     -         3*muu1**2*xd11 + muu2**2*xd11 - xd12 + 4*muu1**2*xd12 + 
     -         2*muu2**2*xd12 - muu3**2*xd12 - muu4**2*xd12 + 
     -         2*tred*xd12 - 2*xd13 + 4*muu1**2*xd13 + 
     -         4*muu2**2*xd13 - 2*xd2 + 3*muu1**2*xd2 + 
     -         3*muu2**2*xd2 - muu3**2*xd2 - muu4**2*xd2 + 
     -         2*tred*xd2 + muu1**2*xd22 + muu2**2*xd22 - 
     -         muu3**2*xd22 - muu4**2*xd22 + 2*tred*xd22 - xd23 + 
     -         2*muu1**2*xd23 + 4*muu2**2*xd23 - muu3**2*xd23 - 
     -         muu4**2*xd23 + 2*tred*xd23 - 3*xd3 + 3*muu1**2*xd3 + 
     -         5*muu2**2*xd3 - xd33 + muu1**2*xd33 + 3*muu2**2*xd33))))
	  
	  vbCRp1 = vbCRp1 +colori*(-(g0R*g1*g2*g3L*(xc0 + 4*xd00 + 
     -      sqrtS**2*(2*(-1 + muu1**2 + muu2**2)*xd0 + 
     -         (-3 + 5*muu1**2 + 3*muu2**2)*xd1 - xd11 + 
     -         3*muu1**2*xd11 + muu2**2*xd11 - xd12 + 4*muu1**2*xd12 + 
     -         2*muu2**2*xd12 - muu3**2*xd12 - muu4**2*xd12 + 
     -         2*tred*xd12 - 2*xd13 + 4*muu1**2*xd13 + 
     -         4*muu2**2*xd13 - 2*xd2 + 3*muu1**2*xd2 + 
     -         3*muu2**2*xd2 - muu3**2*xd2 - muu4**2*xd2 + 
     -         2*tred*xd2 + muu1**2*xd22 + muu2**2*xd22 - 
     -         muu3**2*xd22 - muu4**2*xd22 + 2*tred*xd22 - xd23 + 
     -         2*muu1**2*xd23 + 4*muu2**2*xd23 - muu3**2*xd23 - 
     -         muu4**2*xd23 + 2*tred*xd23 - 3*xd3 + 3*muu1**2*xd3 + 
     -         5*muu2**2*xd3 - xd33 + muu1**2*xd33 + 3*muu2**2*xd33))))
     
	
     
	vbCLp2 = vbCLp2 + colori*0d0
	vbCRp2 = vbCRp2 + colori*0d0
	vbDp1 = vbDp1 + colori*0d0
	vbDk1 = vBDk1 + colori*0d0
	vbDk2 = vbDk2 + colori*0d0

		    
		  enddo
		enddo
	      enddo
	    enddo
	  enddo
	enddo
! ! ! ! 
! ! ! ! ! ! ! ****************************************************
! ! ! ! ! ! ! ************ Box 1.2 Gluon-Squark (i,k,m)- Neutralino- Quark (l,n) ***********************
! ! ! ! ! ! ! ****************************************************
! ! ! ! ! ! ! ! ! Coupling 1 L<->R , Coupling 0 gains minus sign due fermion flow change
! ! 
	   do i = 1,2
	    do k = 3,4
	      do l = 3,4
		do m = 3,3
		  do n = 3,3
 
		  Ml0 = 0d0
		  Ml1 = MSf(i,k,m)

		  if(b.LE.4)then
		    Ml2 = MNeu(b)
		  else
		    if(b.EQ.5) then
		      Ml2 = MCha(1)
		    else
		      Ml2 = MCha(2)
		    endif
		  endif

		  Ml3 = Mf(l,n)
! 		  
		  g0L = GS*Krondelta(l,ftt2)*Krondelta(n,fgen2)
		  g0R = GS*Krondelta(l,ftt2)*Krondelta(n,fgen2)
		  g1L = AL(b,isf2,itt2,l,igen2,n) 
		  g1R = AR(b,isf2,itt2,l,igen2,n)
		  g2 = VSfSf(1,i,isf1,k,itt1,m,igen1)
		  g3L = AL(b,i,k,ftt1,m,fgen1) 
		  g3R = AR(b,i,k,ftt1,m,fgen1)
		  
! 		
		  call myCget(muu3**2*sqrtS**2,ured*sqrtS**2,muu2*sqrtS**2,Ml1**2,Ml2**2,Ml3**2)
		  call myDget(muu1**2*sqrtS**2,muu3**2*sqrtS**2,muu2**2*sqrtS**2,muu4**2*sqrtS**2,tred*sqrtS**2,ured*sqrtS**2, Ml0**2,Ml1**2,Ml2**2,Ml3**2)
		  		  
	colori = Nc*CF

	vbA = vbA + colori*0d0
	
	vbBL = vbBL + colori*g2*(g0L*(g1L*(g3R*muu3*sqrtS*(xc0 + xc1 + 2*xc2 + Ml0**2*xd0 + 2*sqrtS**2*xd0 - 2*muu2**2*sqrtS**2*xd0 - 2*muu3**2*sqrtS**2*xd0 + 2*sqrtS**2*tred*xd0 + sqrtS**2*xd1 + 
     -             2*muu1**2*sqrtS**2*xd1 - muu2**2*sqrtS**2*xd1 - muu3**2*sqrtS**2*xd1 + sqrtS**2*tred*xd1 + sqrtS**2*xd12 - muu2**2*sqrtS**2*xd12 - muu3**2*sqrtS**2*xd12 + 
     -             sqrtS**2*tred*xd12 + Ml0**2*xd2 + 2*sqrtS**2*xd2 - 3*muu2**2*sqrtS**2*xd2 - 2*muu3**2*sqrtS**2*xd2 + muu4**2*sqrtS**2*xd2 + 3*sqrtS**2*tred*xd2 - 
     -             muu2**2*sqrtS**2*xd22 + muu4**2*sqrtS**2*xd22 + sqrtS**2*tred*xd22 + 2*muu4**2*sqrtS**2*xd23 + 2*sqrtS**2*xd3 - 2*muu2**2*sqrtS**2*xd3 - 2*muu3**2*sqrtS**2*xd3 + 
     -             muu4**2*sqrtS**2*xd3 + 2*sqrtS**2*tred*xd3) + g3L*Ml2*(xc0 + Ml0**2*xd0 + 
     -             sqrtS**2*(-2*(-1 + muu2**2 + muu3**2 - tred)*xd0 + (1 + 2*muu1**2 - muu2**2 - muu3**2 + tred)*xd1 + 2*muu1**2*xd2 - muu2**2*xd2 + muu4**2*xd2 + tred*xd2 + 2*xd3 - 
     -                2*muu2**2*xd3 - 2*muu3**2*xd3 + muu4**2*xd3 + 2*tred*xd3))) + 
     -       g1R*Ml3*(g3R*Ml2*muu3*sqrtS*xd2 + g3L*(xc0 + Ml0**2*xd0 + sqrtS**2*(muu1**2*(2*xd0 + 3*xd1 + 2*xd2) - 2*(-1 + muu2**2 + muu3**2)*xd3 + tred*(xd2 + 2*xd3))))) + 
     -    g0R*muu4*sqrtS*(g1L*Ml3*(g3L*Ml2 + g3R*muu3*sqrtS)*xd3 - g1R*(g3R*Ml2*muu3*sqrtS*xd2 + 
     -          g3L*(xc0 + xc2 + Ml0**2*xd0 + 2*muu1**2*sqrtS**2*xd0 - 2*xd00 + 3*muu1**2*sqrtS**2*xd1 - sqrtS**2*xd13 + muu2**2*sqrtS**2*xd13 + muu3**2*sqrtS**2*xd13 - 
     -             sqrtS**2*tred*xd13 + 2*muu1**2*sqrtS**2*xd2 + sqrtS**2*tred*xd2 + muu2**2*sqrtS**2*xd23 - muu4**2*sqrtS**2*xd23 - sqrtS**2*tred*xd23 - Ml0**2*xd3 + 
     -             2*muu1**2*sqrtS**2*xd3 - 2*muu4**2*sqrtS**2*xd33))))
          
	vbBR = vbBR + colori*(g2*(g0R*(-(g1L*Ml3*(-(g3L*Ml2*muu3*sqrtS*xd2) + g3R*(xc0 + Ml0**2*xd0 + sqrtS**2*(muu1**2*(-2*xd0 + xd1 + 2*xd2) - tred*(xd2 - 2*xd3) - 2*(-1 + muu2**2 + muu3**2)*xd3)))) + 
     -         g1R*(g3L*muu3*sqrtS*(xc0 - xc1 - 2*xc2 + Ml0**2*xd0 + 2*sqrtS**2*xd0 - 2*muu2**2*sqrtS**2*xd0 - 2*muu3**2*sqrtS**2*xd0 + 2*sqrtS**2*tred*xd0 + sqrtS**2*xd1 + 
     -               2*muu1**2*sqrtS**2*xd1 - muu2**2*sqrtS**2*xd1 - muu3**2*sqrtS**2*xd1 + sqrtS**2*tred*xd1 - sqrtS**2*xd12 + muu2**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 - 
     -               sqrtS**2*tred*xd12 - Ml0**2*xd2 - 2*sqrtS**2*xd2 + muu2**2*sqrtS**2*xd2 + 2*muu3**2*sqrtS**2*xd2 + muu4**2*sqrtS**2*xd2 - sqrtS**2*tred*xd2 + muu2**2*sqrtS**2*xd22 - 
     -               muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 - 2*muu4**2*sqrtS**2*xd23 + 2*sqrtS**2*xd3 - 2*muu2**2*sqrtS**2*xd3 - 2*muu3**2*sqrtS**2*xd3 + muu4**2*sqrtS**2*xd3 + 
     -               2*sqrtS**2*tred*xd3) + g3R*Ml2*(xc0 + Ml0**2*xd0 + sqrtS**2*
     -                (-2*(-1 + muu2**2 + muu3**2 - tred)*xd0 + (1 + 2*muu1**2 - muu2**2 - muu3**2 + tred)*xd1 + 2*muu1**2*xd2 - muu2**2*xd2 + muu4**2*xd2 + tred*xd2 + 2*xd3 - 
     -                  2*muu2**2*xd3 - 2*muu3**2*xd3 + muu4**2*xd3 + 2*tred*xd3)))) + 
     -      g0L*muu4*sqrtS*(g1R*Ml3*(g3R*Ml2 + g3L*muu3*sqrtS)*xd3 + g1L*(-(g3L*Ml2*muu3*sqrtS*xd2) + 
     -            g3R*(xc0 + xc2 + Ml0**2*xd0 - 2*muu1**2*sqrtS**2*xd0 - 2*xd00 + muu1**2*sqrtS**2*xd1 - sqrtS**2*xd13 + muu2**2*sqrtS**2*xd13 + muu3**2*sqrtS**2*xd13 - 
     -               sqrtS**2*tred*xd13 + 2*muu1**2*sqrtS**2*xd2 - sqrtS**2*tred*xd2 + muu2**2*sqrtS**2*xd23 - muu4**2*sqrtS**2*xd23 - sqrtS**2*tred*xd23 - Ml0**2*xd3 - 
     -               2*muu1**2*sqrtS**2*xd3 - 2*muu4**2*sqrtS**2*xd33)))))
     
	vbCLp1 = vbCLp1 + colori*(g2*(g0R*muu4*sqrtS*(-3d0*g1L*g3R*Ml3*xd3 + g1R*(g3L*muu3*sqrtS*(2*xd0 + xd1 - 2*xd2 + 2*xd3) + g3R*Ml2*(2*xd0 + xd1 + xd2 + 2*xd3))) + 
     -      g0L*(-(g1R*Ml3*(g3L*muu3*sqrtS*(2*xd0 + xd1 - 2*xd2) + g3R*Ml2*(2*xd0 + xd1 + xd2))) + 
     -         g1L*(-2*g3L*Ml2*muu3*sqrtS*xd2 + g3R*(xc2 + Ml0**2*(xd0 + xd1 + xd2) + 
     -               sqrtS**2*(2*(-1 + muu2**2 + muu3**2 - tred)*xd0 - (-1 + 2*muu1**2 + muu2**2 + muu3**2 - tred)*xd1 + xd11 - muu2**2*xd11 - muu3**2*xd11 + tred*xd11 + xd12 - 
     -                  2*muu2**2*xd12 - muu3**2*xd12 + muu4**2*xd12 + 2*tred*xd12 + 2*muu4**2*xd13 + 2*xd2 - muu2**2*xd2 - 2*muu3**2*xd2 - muu4**2*xd2 - tred*xd2 - muu2**2*xd22 + 
     -                  muu4**2*xd22 + tred*xd22 + 2*muu4**2*xd23 - 2*xd3 + 2*muu2**2*xd3 + 2*muu3**2*xd3 + muu4**2*xd3 - 2*tred*xd3))))))
     
	vbCRp1 = vbCRp1 + colori*(g2*(g0L*muu4*sqrtS*(-3d0*g1R*g3L*Ml3*xd3 + g1L*(g3R*muu3*sqrtS*(2*xd0 + xd1 - 2*xd2 + 2*xd3) + g3L*Ml2*(2*xd0 + xd1 + xd2 + 2*xd3))) + 
     -      g0R*(-(g1L*Ml3*(g3R*muu3*sqrtS*(2*xd0 + xd1 - 2*xd2) + g3L*Ml2*(2*xd0 + xd1 + xd2))) + 
     -         g1R*(-2*g3R*Ml2*muu3*sqrtS*xd2 + g3L*(xc2 + Ml0**2*(xd0 + xd1 + xd2) + 
     -               sqrtS**2*(2*(-1 + muu2**2 + muu3**2 - tred)*xd0 - (-1 + 2*muu1**2 + muu2**2 + muu3**2 - tred)*xd1 + xd11 - muu2**2*xd11 - muu3**2*xd11 + tred*xd11 + xd12 - 
     -                  2*muu2**2*xd12 - muu3**2*xd12 + muu4**2*xd12 + 2*tred*xd12 + 2*muu4**2*xd13 + 2*xd2 - muu2**2*xd2 - 2*muu3**2*xd2 - muu4**2*xd2 - tred*xd2 - muu2**2*xd22 + 
     -                  muu4**2*xd22 + tred*xd22 + 2*muu4**2*xd23 - 2*xd3 + 2*muu2**2*xd3 + 2*muu3**2*xd3 + muu4**2*xd3 - 2*tred*xd3))))))
     
      vbCLp2 = vbCLp2 + colori*0d0
	  vbCRp2 = vbCRp2 + colori*0d0
	  vbDp1 = vbDp1 + colori*0d0
	  vbDk1 = vBDk1 + colori*0d0
	  vbDk2 = vbDk2 + colori*0d0
	
		   
		  enddo
		enddo
	      enddo
	    enddo
	enddo



! ! ! 
! ! ! ! ! ****************************************************
! ! ! ! ! ************ Box 1.3 Neutralino-Squark (i,k,m)- Gluon- Quark (l,n) ***********************
! ! ! ! ! ****************************************************
! ! ! ! ! Coupling 0 L<->R , Coupling 1 gains minus sign due to charge flow change

	do j = 1,2
	    do k = 3,4
	      do l = 3,4
		do m = 3,3
		  do n = 3,3
		  
		  Ml0 = 0d0
		  Ml1 = MSf(j,l,n)
		  
		  if(b.LE.4)then
		    Ml2 = MNeu(b)
		  else
		    if(b.EQ.5) then
		      Ml2 = MCha(1)
		    else
		      Ml2 = MCha(2)
		    endif
		  endif
		  
		  Ml3 = Mf(k,m)
		  
		  g0L = AL(b,j,l,ftt2,n,fgen2) 
		  g0R = AR(b,j,l,ftt2,n,fgen2)
		  g1 = -VSfSf(1,isf2,j,itt2,l,igen2,n)
		  g2L = AL(b,isf1,itt1,k,igen1,m) 
		  g2R = AR(b,isf1,itt1,k,igen1,m)
		  g3L = -GS*Krondelta(k,ftt1)*Krondelta(m,fgen1)
		  g3R = -GS*Krondelta(k,ftt1)*Krondelta(m,fgen1)
! 		
		  call myCget(muu3**2*sqrtS**2,ured*sqrtS**2,muu2*sqrtS**2,Ml1**2,Ml2**2,Ml3**2)
		  call myDget(muu2**2*sqrtS**2,muu4**2*sqrtS**2,muu1**2*sqrtS**2,muu3**2*sqrtS**2,tred*sqrtS**2,ured*sqrtS**2, Ml0**2,Ml1**2,Ml2**2,Ml3**2)
		 		  
	colori = Nc*CF

	vbA = vbA + colori*0d0
	
	vbBL = vbBL + colori*g1*(g0R*muu4*sqrtS*(g2R*Ml2*(-(g3R*Ml3*(2*xd0 + xd1)) + g3L*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3))) - 
     -       g2L*(g3R*(xc0 - xc1 - xc2 + 2*Ml0**2*xd0 + 2*muu3**2*sqrtS**2*xd0 + Ml0**2*xd1 + 2*sqrtS**2*xd1 - 2*muu1**2*sqrtS**2*xd1 + muu3**2*sqrtS**2*xd1 - 2*muu4**2*sqrtS**2*xd1 + 
     -             2*sqrtS**2*tred*xd1 + sqrtS**2*xd11 - muu1**2*sqrtS**2*xd11 - muu4**2*sqrtS**2*xd11 + sqrtS**2*tred*xd11 - muu1**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 + 
     -             sqrtS**2*tred*xd12 + 2*muu3**2*sqrtS**2*xd13 + 2*muu3**2*sqrtS**2*xd2 + 2*sqrtS**2*tred*xd2 + 4*muu3**2*sqrtS**2*xd3) - g3L*Ml3*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3)))
     -       ) + g0L*(g2L*Ml2*(g3L*Ml3*muu3*sqrtS*(2*xd0 + xd1 + xd2 + xd3) + 
     -          g3R*(xc0 + Ml0**2*xd0 - sqrtS**2*(2*(-1 + muu1**2 + muu3**2 + muu4**2 - tred)*xd0 + (-1 + muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 - tred)*xd1 + muu1**2*xd2 - 
     -                2*muu2**2*xd2 + 2*muu4**2*xd2 - tred*xd2 - 2*xd3 + 2*muu1**2*xd3 + muu3**2*xd3 + 2*muu4**2*xd3 - 2*tred*xd3))) + 
     -       g2R*(g3R*Ml3*(xc0 + Ml0**2*xd0 + sqrtS**2*(tred*xd2 + muu2**2*(2*xd0 + 3*xd1 + 2*xd2) + 2*xd3 - 2*muu1**2*xd3 - muu3**2*xd3 + 2*tred*xd3 - 
     -                muu4**2*(2*xd0 + xd1 + 2*(xd2 + xd3)))) + g3L*muu3*sqrtS*
     -           (xc0 + 2*sqrtS**2*xd0 - 2*muu1**2*sqrtS**2*xd0 - 2*muu2**2*sqrtS**2*xd0 + 2*sqrtS**2*tred*xd0 + 2*xd00 + 3*sqrtS**2*xd1 - 3*muu1**2*sqrtS**2*xd1 - muu2**2*sqrtS**2*xd1 - 
     -             2*muu4**2*sqrtS**2*xd1 + 3*sqrtS**2*tred*xd1 + sqrtS**2*xd11 - muu1**2*sqrtS**2*xd11 - muu4**2*sqrtS**2*xd11 + sqrtS**2*tred*xd11 + sqrtS**2*xd12 - 
     -             2*muu1**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 - muu4**2*sqrtS**2*xd12 + 2*sqrtS**2*tred*xd12 + sqrtS**2*xd13 - muu1**2*sqrtS**2*xd13 + 2*muu3**2*sqrtS**2*xd13 - 
     -             muu4**2*sqrtS**2*xd13 + sqrtS**2*tred*xd13 + 2*sqrtS**2*xd2 - 3*muu1**2*sqrtS**2*xd2 - 2*muu2**2*sqrtS**2*xd2 + muu3**2*sqrtS**2*xd2 + 4*sqrtS**2*tred*xd2 - 
     -             muu1**2*sqrtS**2*xd22 + muu3**2*sqrtS**2*xd22 + sqrtS**2*tred*xd22 - muu1**2*sqrtS**2*xd23 + 3*muu3**2*sqrtS**2*xd23 + sqrtS**2*tred*xd23 + 2*sqrtS**2*xd3 - 
     -             2*muu1**2*sqrtS**2*xd3 - 2*muu2**2*sqrtS**2*xd3 + 3*muu3**2*sqrtS**2*xd3 + 2*sqrtS**2*tred*xd3 + Ml0**2*(2*xd0 + xd1 + xd2 + xd3) + 2*muu3**2*sqrtS**2*xd33))))
     
	vbBR = vbBR + colori*g1*(g0L*muu4*sqrtS*(g2L*Ml2*(-(g3L*Ml3*(2*xd0 + xd1)) + g3R*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3))) - 
     -       g2R*(g3L*(xc0 - xc1 - xc2 + 2*Ml0**2*xd0 + 2*muu3**2*sqrtS**2*xd0 + Ml0**2*xd1 + 2*sqrtS**2*xd1 - 2*muu1**2*sqrtS**2*xd1 + muu3**2*sqrtS**2*xd1 - 2*muu4**2*sqrtS**2*xd1 + 
     -             2*sqrtS**2*tred*xd1 + sqrtS**2*xd11 - muu1**2*sqrtS**2*xd11 - muu4**2*sqrtS**2*xd11 + sqrtS**2*tred*xd11 - muu1**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 + 
     -             sqrtS**2*tred*xd12 + 2*muu3**2*sqrtS**2*xd13 + 2*muu3**2*sqrtS**2*xd2 + 2*sqrtS**2*tred*xd2 + 4*muu3**2*sqrtS**2*xd3) - g3R*Ml3*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3)))
     -       ) + g0R*(g2R*Ml2*(g3R*Ml3*muu3*sqrtS*(2*xd0 + xd1 + xd2 + xd3) + 
     -          g3L*(xc0 + Ml0**2*xd0 - sqrtS**2*(2*(-1 + muu1**2 + muu3**2 + muu4**2 - tred)*xd0 + (-1 + muu1**2 - 2*muu2**2 + muu3**2 + muu4**2 - tred)*xd1 + muu1**2*xd2 - 
     -                2*muu2**2*xd2 + 2*muu4**2*xd2 - tred*xd2 - 2*xd3 + 2*muu1**2*xd3 + muu3**2*xd3 + 2*muu4**2*xd3 - 2*tred*xd3))) + 
     -       g2L*(g3L*Ml3*(xc0 + Ml0**2*xd0 + sqrtS**2*(tred*xd2 + muu2**2*(2*xd0 + 3*xd1 + 2*xd2) + 2*xd3 - 2*muu1**2*xd3 - muu3**2*xd3 + 2*tred*xd3 - 
     -                muu4**2*(2*xd0 + xd1 + 2*(xd2 + xd3)))) + g3R*muu3*sqrtS*
     -           (xc0 + 2*sqrtS**2*xd0 - 2*muu1**2*sqrtS**2*xd0 - 2*muu2**2*sqrtS**2*xd0 + 2*sqrtS**2*tred*xd0 + 2*xd00 + 3*sqrtS**2*xd1 - 3*muu1**2*sqrtS**2*xd1 - muu2**2*sqrtS**2*xd1 - 
     -             2*muu4**2*sqrtS**2*xd1 + 3*sqrtS**2*tred*xd1 + sqrtS**2*xd11 - muu1**2*sqrtS**2*xd11 - muu4**2*sqrtS**2*xd11 + sqrtS**2*tred*xd11 + sqrtS**2*xd12 - 
     -             2*muu1**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 - muu4**2*sqrtS**2*xd12 + 2*sqrtS**2*tred*xd12 + sqrtS**2*xd13 - muu1**2*sqrtS**2*xd13 + 2*muu3**2*sqrtS**2*xd13 - 
     -             muu4**2*sqrtS**2*xd13 + sqrtS**2*tred*xd13 + 2*sqrtS**2*xd2 - 3*muu1**2*sqrtS**2*xd2 - 2*muu2**2*sqrtS**2*xd2 + muu3**2*sqrtS**2*xd2 + 4*sqrtS**2*tred*xd2 - 
     -             muu1**2*sqrtS**2*xd22 + muu3**2*sqrtS**2*xd22 + sqrtS**2*tred*xd22 - muu1**2*sqrtS**2*xd23 + 3*muu3**2*sqrtS**2*xd23 + sqrtS**2*tred*xd23 + 2*sqrtS**2*xd3 - 
     -             2*muu1**2*sqrtS**2*xd3 - 2*muu2**2*sqrtS**2*xd3 + 3*muu3**2*sqrtS**2*xd3 + 2*sqrtS**2*tred*xd3 + Ml0**2*(2*xd0 + xd1 + xd2 + xd3) + 2*muu3**2*sqrtS**2*xd33))))
     
	vbCLp1 = vbCLp1 - colori*(g1*(g0L*(g2L*Ml2*(g3L*Ml3*(2*xd0 + xd1 + xd2) - g3R*muu3*sqrtS*(2*xd0 + xd1 + xd2 + 2*xd3)) + 
     -         g2R*(-(g3R*Ml3*muu3*sqrtS*xd3) + g3L*(2*xc0 + xc2 + 3*Ml0**2*xd0 + 2*sqrtS**2*xd0 - 2*muu1**2*sqrtS**2*xd0 - 2*muu4**2*sqrtS**2*xd0 + 2*sqrtS**2*tred*xd0 + Ml0**2*xd1 + 
     -               3*sqrtS**2*xd1 - 3*muu1**2*sqrtS**2*xd1 + 2*muu2**2*sqrtS**2*xd1 - 3*muu4**2*sqrtS**2*xd1 + 3*sqrtS**2*tred*xd1 + sqrtS**2*xd11 - muu1**2*sqrtS**2*xd11 - 
     -               muu4**2*sqrtS**2*xd11 + sqrtS**2*tred*xd11 + sqrtS**2*xd12 - 2*muu1**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 - muu4**2*sqrtS**2*xd12 + 2*sqrtS**2*tred*xd12 + 
     -               2*muu3**2*sqrtS**2*xd13 + Ml0**2*xd2 + 2*sqrtS**2*xd2 - 3*muu1**2*sqrtS**2*xd2 + muu3**2*sqrtS**2*xd2 - 2*muu4**2*sqrtS**2*xd2 + 5*sqrtS**2*tred*xd2 - 
     -               muu1**2*sqrtS**2*xd22 + muu3**2*sqrtS**2*xd22 + sqrtS**2*tred*xd22 + 2*muu3**2*sqrtS**2*xd23 + 2*sqrtS**2*xd3 - 2*muu1**2*sqrtS**2*xd3 + 3*muu3**2*sqrtS**2*xd3 - 
     -               2*muu4**2*sqrtS**2*xd3 + 2*sqrtS**2*tred*xd3))) + g0R*muu4*sqrtS*
     -       (2*g2R*g3L*Ml2*xd2 + g2L*(g3L*Ml3*(2*xd0 + xd1 + 2*xd2) - g3R*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3))))))
     
	vbCRp1 = vbCRp1 - colori*(g1*(g0R*(g2R*Ml2*(g3R*Ml3*(2*xd0 + xd1 + xd2) - g3L*muu3*sqrtS*(2*xd0 + xd1 + xd2 + 2*xd3)) + 
     -         g2L*(-(g3L*Ml3*muu3*sqrtS*xd3) + g3R*(2*xc0 + xc2 + 3*Ml0**2*xd0 + 2*sqrtS**2*xd0 - 2*muu1**2*sqrtS**2*xd0 - 2*muu4**2*sqrtS**2*xd0 + 2*sqrtS**2*tred*xd0 + Ml0**2*xd1 + 
     -               3*sqrtS**2*xd1 - 3*muu1**2*sqrtS**2*xd1 + 2*muu2**2*sqrtS**2*xd1 - 3*muu4**2*sqrtS**2*xd1 + 3*sqrtS**2*tred*xd1 + sqrtS**2*xd11 - muu1**2*sqrtS**2*xd11 - 
     -               muu4**2*sqrtS**2*xd11 + sqrtS**2*tred*xd11 + sqrtS**2*xd12 - 2*muu1**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 - muu4**2*sqrtS**2*xd12 + 2*sqrtS**2*tred*xd12 + 
     -               2*muu3**2*sqrtS**2*xd13 + Ml0**2*xd2 + 2*sqrtS**2*xd2 - 3*muu1**2*sqrtS**2*xd2 + muu3**2*sqrtS**2*xd2 - 2*muu4**2*sqrtS**2*xd2 + 5*sqrtS**2*tred*xd2 - 
     -               muu1**2*sqrtS**2*xd22 + muu3**2*sqrtS**2*xd22 + sqrtS**2*tred*xd22 + 2*muu3**2*sqrtS**2*xd23 + 2*sqrtS**2*xd3 - 2*muu1**2*sqrtS**2*xd3 + 3*muu3**2*sqrtS**2*xd3 - 
     -               2*muu4**2*sqrtS**2*xd3 + 2*sqrtS**2*tred*xd3))) + g0L*muu4*sqrtS*
     -       (2*g2L*g3R*Ml2*xd2 + g2R*(g3R*Ml3*(2*xd0 + xd1 + 2*xd2) - g3L*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3))))))
     

      vbCLp2 = vbCLp2 + colori*0d0
	  vbCRp2 = vbCRp2 + colori*0d0
	  vbDp1 = vbDp1 + colori*0d0
	  vbDk1 = vBDk1 + colori*0d0
	  vbDk2 = vbDk2 + colori*0d0
	
		  enddo
		enddo
	      enddo
	    enddo
	enddo
 
! 
! ! ! ! ! 
! ! ! ! ! ! ****************************************************
! ! ! ! ! ! ! ************ Box 4 Neutralino-Quark(k,m)-Gluon- Quark(l,n) ***********************
! ! ! ! ! ! ! ****************************************************
! ! ! !! ! ! Coupling 1 L<->R , Coupling 0 gains minus sign due to fermion flow change
! !  	
	  do k = 3,4
	    do l = 3,4
	      do m = 3,3
		do n = 3,3
		  
		
		  if(b.LE.4)then
		    Ml0 = MNeu(b)
		  else
		    if(b.EQ.5) then
		      Ml0 = MCha(1)
		    else
		      Ml0 = MCha(2)
		    endif
		  endif
		  Ml1 = Mf(k,m)
		  Ml2 = 0
		  Ml3 = Mf(l,n)
		  
		  g0 = GS*Krondelta(l,ftt2)*Krondelta(n,fgen2)
		  g1L = AL(b,isf2,itt2,l,igen2,n)
		  g1R = AR(b,isf2,itt2,l,igen2,n)
		  g2L = AL(b,isf1,itt1,k,igen1,m)
		  g2R = AR(b,isf1,itt1,k,igen1,m)
		  g3 = -GS*Krondelta(k,ftt1)*Krondelta(m,fgen1)
		  
		  call myCget(muu3**2*sqrtS**2,muu4**2*sqrtS**2,sqrtS**2,Ml1**2,Ml2**2,Ml3**2)
		  call myDget(muu1**2*sqrtS**2,muu3**2*sqrtS**2,muu4**2*sqrtS**2,muu2**2*sqrtS**2,tred*sqrtS**2,sqrtS**2, Ml0**2,Ml1**2,Ml2**2,Ml3**2)
		
	colori = Nc*CF

	vbA = vbA + colori*0d0
	
	vbBL = vbBL - colori*2*g0*g3*(g1R*(g2R*Ml0*(Ml3*muu3*sqrtS*(xd2 + xd3) + Ml1*(2*Ml3*xd0 - muu4*sqrtS*(xd0 + xd3))) + 
     -       g2L*(-(muu4*sqrtS*(xc0 + xc2 + Ml0**2*xd0 + muu1**2*sqrtS**2*xd1 + muu1**2*sqrtS**2*xd2 + 
     -               Ml0**2*xd3)) + Ml3*(2*xc0 + 2*Ml0**2*xd0 + 
     -             sqrtS*(Ml1*muu3*(xd2 + xd3) + muu1**2*sqrtS*(2*xd1 + xd2 + xd3) + 
     -                sqrtS*(-(muu3**2*xd2) + tred*xd2 + (-1 + muu2**2)*xd3))))) + 
     -    g1L*(g2L*Ml0*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*
     -           (-(Ml3*muu4*xd3) + Ml1*muu3*(xd0 + xd2 + xd3) + 
     -             sqrtS*((-1 + muu1**2 + muu2**2)*xd0 + (-1 + 3*muu1**2 + muu2**2)*xd1 + muu1**2*xd2 + 
     -                muu2**2*xd2 - muu3**2*xd2 - muu4**2*xd2 + 2*tred*xd2 - xd3 + muu1**2*xd3 + 3*muu2**2*xd3)))
     -        + g2R*(muu3*sqrtS*(xc0 + xc1 + xc2 + Ml0**2*xd0 + muu1**2*sqrtS**2*xd1 + Ml0**2*xd2 + sqrtS**2*xd2 - 
     -             muu2**2*sqrtS**2*xd2 - muu3**2*sqrtS**2*xd2 + sqrtS**2*tred*xd2 + Ml0**2*xd3) + 
     -          Ml1*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*
     -              (-(Ml3*muu4*xd3) + sqrtS*((-1 + muu1**2 + muu2**2)*xd1 + (-muu4**2 + tred)*xd2 + 
     -                   muu2**2*(xd2 + 2*xd3)))))))
     
	vbBR = vbBR - colori*2*g0*g3*(g1L*(g2L*Ml0*(Ml3*muu3*sqrtS*(xd2 + xd3) + Ml1*(2*Ml3*xd0 - muu4*sqrtS*(xd0 + xd3))) + 
     -       g2R*(-(muu4*sqrtS*(xc0 + xc2 + Ml0**2*xd0 + muu1**2*sqrtS**2*xd1 + muu1**2*sqrtS**2*xd2 + 
     -               Ml0**2*xd3)) + Ml3*(2*xc0 + 2*Ml0**2*xd0 + 
     -             sqrtS*(Ml1*muu3*(xd2 + xd3) + muu1**2*sqrtS*(2*xd1 + xd2 + xd3) + 
     -                sqrtS*(-(muu3**2*xd2) + tred*xd2 + (-1 + muu2**2)*xd3))))) + 
     -    g1R*(g2R*Ml0*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*
     -           (-(Ml3*muu4*xd3) + Ml1*muu3*(xd0 + xd2 + xd3) + 
     -             sqrtS*((-1 + muu1**2 + muu2**2)*xd0 + (-1 + 3*muu1**2 + muu2**2)*xd1 + muu1**2*xd2 + 
     -                muu2**2*xd2 - muu3**2*xd2 - muu4**2*xd2 + 2*tred*xd2 - xd3 + muu1**2*xd3 + 3*muu2**2*xd3)))
     -        + g2L*(muu3*sqrtS*(xc0 + xc1 + xc2 + Ml0**2*xd0 + muu1**2*sqrtS**2*xd1 + Ml0**2*xd2 + sqrtS**2*xd2 - 
     -             muu2**2*sqrtS**2*xd2 - muu3**2*sqrtS**2*xd2 + sqrtS**2*tred*xd2 + Ml0**2*xd3) + 
     -          Ml1*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*
     -              (-(Ml3*muu4*xd3) + sqrtS*((-1 + muu1**2 + muu2**2)*xd1 + (-muu4**2 + tred)*xd2 + 
     -                   muu2**2*(xd2 + 2*xd3)))))))
	
	vbCLp1 = vbCLp1 + colori*2*g0*g3*(g1R*(g2R*Ml0*Ml3*(xd0 + xd1 + xd2 + xd3) + 
     -       g2L*(-(muu3*muu4*sqrtS**2*xd2) + Ml1*Ml3*(xd1 + xd2 + xd3))) + 
     -    g1L*(g2L*Ml0*Ml1*(xd0 + xd1 + xd2 + xd3) + 
     -       g2R*(xc0 + Ml0**2*(2*xd0 + xd1 + xd2 + xd3) + 
     -          sqrtS**2*(muu1**2*xd1 - (-1 + muu3**2 + muu4**2 - tred)*xd2 + muu2**2*xd3))))
	
	vbCRp1 = vbCRp1 + colori*2*g0*g3*(g1L*(g2L*Ml0*Ml3*(xd0 + xd1 + xd2 + xd3) + 
     -       g2R*(-(muu3*muu4*sqrtS**2*xd2) + Ml1*Ml3*(xd1 + xd2 + xd3))) + 
     -    g1R*(g2R*Ml0*Ml1*(xd0 + xd1 + xd2 + xd3) + 
     -       g2L*(xc0 + Ml0**2*(2*xd0 + xd1 + xd2 + xd3) + 
     -          sqrtS**2*(muu1**2*xd1 - (-1 + muu3**2 + muu4**2 - tred)*xd2 + muu2**2*xd3))))
	

      vbCLp2 = vbCLp2 + colori*0d0
	  vbCRp2 = vbCRp2 + colori*0d0
	  vbDp1 = vbDp1 + colori*0d0
	  vbDk1 = vBDk1 + colori*0d0
	  vbDk2 = vbDk2 + colori*0d0
	
		  enddo
		enddo
	      enddo
	    enddo
! ! 	    
! ! ! ! ! !-------------------------------------------------------------------------------------------------
! ! ! ! ! !--------------------------------------u-channels-------------------------------------------------
! ! ! ! ! !-------------------------------------------------------------------------------------------------
! ! ! ! 
! ! ! ! ! ! ! ****************************************************
! ! ! ! ! ! ! ************ Box 1 Gluon-Squark (i,k,m)- Neutralino- Squark (j,l,n) ***********************
! ! ! ! ! ! ! ****************************************************
! ! ! ! ! ! Coupling 0 L<->R , Coupling 1 gains minus sign due to charge flow change
! ! 
	  do i = 1,2
	   do j = 1,2
        do k = 3,4
	     do l = 3,4
		  do m = 3,3
		   do n = 3,3
		  
		  Ml0 = 0d0
		  Ml1 = MSf(i,k,m)
		 
		  if(b.LE.4)then
		    Ml2 = MNeu(b)
		  else
		    if(b.EQ.5) then
		      Ml2 = MCha(1)
		    else
		      Ml2 = MCha(2)
		    endif
		  endif
		  
		  Ml3 = MSf(j,l,n)
		  
		  
		  g0L = AL(b,j,l,ftt2,n,fgen2)
		  g0R = AR(b,j,l,ftt2,n,fgen2)
		  g1 = VSfSf(1,isf1,j,itt1,l,igen1,n)
		  g2 = VSfSf(1,i,isf2,k,itt2,m,igen2)
		  g3L = AL(b,i,k,ftt1,m,fgen1) 
		  g3R = AR(b,i,k,ftt1,m,fgen1)
		  
		  call myCget(muu3**2*sqrtS**2,muu4**2*sqrtS**2,sqrtS**2,Ml1**2,Ml2**2,Ml3**2)
		  call myDget(muu2**2*sqrtS**2,muu3**2*sqrtS**2,muu4**2*sqrtS**2,muu1**2*sqrtS**2,ured*sqrtS**2,sqrtS**2, Ml0**2,Ml1**2,Ml2**2,Ml3**2)
		  		  
	  colori = Nc*CF
	  
	  
	    vbBL = vbBL + colori*(g1*g2*(g0L*g3L*Ml2*(xc0 + 
     -       sqrtS**2*(2*(-1 + muu1**2 + muu2**2)*xd0 + 
     -          (-1 + muu1**2 + 3*muu2**2)*xd1 - 2*xd2 + 
     -          3*muu1**2*xd2 + 3*muu2**2*xd2 + muu3**2*xd2 + 
     -          muu4**2*xd2 - 2*tred*xd2 - xd3 + 3*muu1**2*xd3 + 
     -          muu2**2*xd3)) + 
     -    g0L*g3R*muu3*sqrtS*
     -     (xc0 + xc1 + xc2 - 2*sqrtS**2*xd0 + 
     -       2*muu1**2*sqrtS**2*xd0 + 2*muu2**2*sqrtS**2*xd0 + 
     -       2*xd00 - sqrtS**2*xd1 + muu1**2*sqrtS**2*xd1 + 
     -       3*muu2**2*sqrtS**2*xd1 - sqrtS**2*xd12 + 
     -       muu1**2*sqrtS**2*xd12 + 3*muu2**2*sqrtS**2*xd12 - 
     -       sqrtS**2*xd13 + muu1**2*sqrtS**2*xd13 + 
     -       3*muu2**2*sqrtS**2*xd13 - 4*sqrtS**2*xd2 + 
     -       5*muu1**2*sqrtS**2*xd2 + 5*muu2**2*sqrtS**2*xd2 + 
     -       muu3**2*sqrtS**2*xd2 + muu4**2*sqrtS**2*xd2 - 
     -       2*sqrtS**2*tred*xd2 - 2*sqrtS**2*xd22 + 
     -       3*muu1**2*sqrtS**2*xd22 + 3*muu2**2*sqrtS**2*xd22 + 
     -       muu3**2*sqrtS**2*xd22 + muu4**2*sqrtS**2*xd22 - 
     -       2*sqrtS**2*tred*xd22 - 3*sqrtS**2*xd23 + 
     -       6*muu1**2*sqrtS**2*xd23 + 4*muu2**2*sqrtS**2*xd23 + 
     -       muu3**2*sqrtS**2*xd23 + muu4**2*sqrtS**2*xd23 - 
     -       2*sqrtS**2*tred*xd23 - 3*sqrtS**2*xd3 + 
     -       5*muu1**2*sqrtS**2*xd3 + 3*muu2**2*sqrtS**2*xd3 - 
     -       sqrtS**2*xd33 + 3*muu1**2*sqrtS**2*xd33 + 
     -       muu2**2*sqrtS**2*xd33) - 
     -    g0R*g3L*muu4*sqrtS*
     -     (xc2 + 2*xd00 + sqrtS**2*
     -        ((-1 + muu1**2 + 3*muu2**2)*xd13 + 
     -          (-2 + 3*muu1**2 + 3*muu2**2 + muu3**2 + muu4**2 - 
     -             2*tred)*xd23 - 2*xd3 + 2*muu1**2*xd3 + 
     -          2*muu2**2*xd3 - xd33 + 3*muu1**2*xd33 + muu2**2*xd33))))
     
     
	    vbBR = vbBR + colori*(g1*g2*(g0R*g3R*Ml2*(xc0 + 
     -       sqrtS**2*(2*(-1 + muu1**2 + muu2**2)*xd0 + 
     -          (-1 + muu1**2 + 3*muu2**2)*xd1 - 2*xd2 + 
     -          3*muu1**2*xd2 + 3*muu2**2*xd2 + muu3**2*xd2 + 
     -          muu4**2*xd2 - 2*tred*xd2 - xd3 + 3*muu1**2*xd3 + 
     -          muu2**2*xd3)) + 
     -    g0R*g3L*muu3*sqrtS*
     -     (xc0 + xc1 + xc2 - 2*sqrtS**2*xd0 + 
     -       2*muu1**2*sqrtS**2*xd0 + 2*muu2**2*sqrtS**2*xd0 + 
     -       2*xd00 - sqrtS**2*xd1 + muu1**2*sqrtS**2*xd1 + 
     -       3*muu2**2*sqrtS**2*xd1 - sqrtS**2*xd12 + 
     -       muu1**2*sqrtS**2*xd12 + 3*muu2**2*sqrtS**2*xd12 - 
     -       sqrtS**2*xd13 + muu1**2*sqrtS**2*xd13 + 
     -       3*muu2**2*sqrtS**2*xd13 - 4*sqrtS**2*xd2 + 
     -       5*muu1**2*sqrtS**2*xd2 + 5*muu2**2*sqrtS**2*xd2 + 
     -       muu3**2*sqrtS**2*xd2 + muu4**2*sqrtS**2*xd2 - 
     -       2*sqrtS**2*tred*xd2 - 2*sqrtS**2*xd22 + 
     -       3*muu1**2*sqrtS**2*xd22 + 3*muu2**2*sqrtS**2*xd22 + 
     -       muu3**2*sqrtS**2*xd22 + muu4**2*sqrtS**2*xd22 - 
     -       2*sqrtS**2*tred*xd22 - 3*sqrtS**2*xd23 + 
     -       6*muu1**2*sqrtS**2*xd23 + 4*muu2**2*sqrtS**2*xd23 + 
     -       muu3**2*sqrtS**2*xd23 + muu4**2*sqrtS**2*xd23 - 
     -       2*sqrtS**2*tred*xd23 - 3*sqrtS**2*xd3 + 
     -       5*muu1**2*sqrtS**2*xd3 + 3*muu2**2*sqrtS**2*xd3 - 
     -       sqrtS**2*xd33 + 3*muu1**2*sqrtS**2*xd33 + 
     -       muu2**2*sqrtS**2*xd33) - 
     -    g0L*g3R*muu4*sqrtS*
     -     (xc2 + 2*xd00 + sqrtS**2*
     -        ((-1 + muu1**2 + 3*muu2**2)*xd13 + 
     -          (-2 + 3*muu1**2 + 3*muu2**2 + muu3**2 + muu4**2 - 
     -             2*tred)*xd23 - 2*xd3 + 2*muu1**2*xd3 + 
     -          2*muu2**2*xd3 - xd33 + 3*muu1**2*xd33 + muu2**2*xd33))))
      
	    vbCLp2 = vbCLp2 + colori*( -(g0L*g1*g2*g3R*(xc0 + 4*xd00 + 
     -      sqrtS**2*(2*(-1 + muu1**2 + muu2**2)*xd0 + 
     -         (-3 + 3*muu1**2 + 5*muu2**2)*xd1 - xd11 + 
     -         muu1**2*xd11 + 3*muu2**2*xd11 - 3*xd12 + 
     -         4*muu1**2*xd12 + 6*muu2**2*xd12 + muu3**2*xd12 + 
     -         muu4**2*xd12 - 2*tred*xd12 - 2*xd13 + 4*muu1**2*xd13 + 
     -         4*muu2**2*xd13 - 4*xd2 + 5*muu1**2*xd2 + 
     -         5*muu2**2*xd2 + muu3**2*xd2 + muu4**2*xd2 - 
     -         2*tred*xd2 - 2*xd22 + 3*muu1**2*xd22 + 
     -         3*muu2**2*xd22 + muu3**2*xd22 + muu4**2*xd22 - 
     -         2*tred*xd22 - 3*xd23 + 6*muu1**2*xd23 + 
     -         4*muu2**2*xd23 + muu3**2*xd23 + muu4**2*xd23 - 
     -         2*tred*xd23 - 3*xd3 + 5*muu1**2*xd3 + 3*muu2**2*xd3 - 
     -         xd33 + 3*muu1**2*xd33 + muu2**2*xd33))))
   
	    vbCRp2 = vbCRp2 + colori*(-(g0R*g1*g2*g3L*(xc0 + 4*xd00 + 
     -      sqrtS**2*(2*(-1 + muu1**2 + muu2**2)*xd0 + 
     -         (-3 + 3*muu1**2 + 5*muu2**2)*xd1 - xd11 + 
     -         muu1**2*xd11 + 3*muu2**2*xd11 - 3*xd12 + 
     -         4*muu1**2*xd12 + 6*muu2**2*xd12 + muu3**2*xd12 + 
     -         muu4**2*xd12 - 2*tred*xd12 - 2*xd13 + 4*muu1**2*xd13 + 
     -         4*muu2**2*xd13 - 4*xd2 + 5*muu1**2*xd2 + 
     -         5*muu2**2*xd2 + muu3**2*xd2 + muu4**2*xd2 - 
     -         2*tred*xd2 - 2*xd22 + 3*muu1**2*xd22 + 
     -         3*muu2**2*xd22 + muu3**2*xd22 + muu4**2*xd22 - 
     -         2*tred*xd22 - 3*xd23 + 6*muu1**2*xd23 + 
     -         4*muu2**2*xd23 + muu3**2*xd23 + muu4**2*xd23 - 
     -         2*tred*xd23 - 3*xd3 + 5*muu1**2*xd3 + 3*muu2**2*xd3 - 
     -         xd33 + 3*muu1**2*xd33 + muu2**2*xd33))))
      vbCLp1 = vbCLp1 + colori*0d0
	  vbCRp1 = vbCRp1 + colori*0d0
	  vbDp1 = vbDp1 + colori*0d0
	  vbDk1 = vBDk1 + colori*0d0
	  vbDk2 = vbDk2 + colori*0d0

		    
		   enddo
		  enddo
	     enddo
	    enddo
	   enddo
	  enddo
! ! 
! ! ! ! ! ! ! ****************************************************
! ! ! ! ! ! ! ************ Box 1.2 Gluon-Squark (i,k,m)- Neutralino- Quark (l,n) ***********************
! ! ! ! ! ! ! ****************************************************
! ! ! ! ! ! ! ! ! Coupling 1 L<->R , Coupling 0 gains minus sign due fermion flow change
! ! ! 
	   do i = 1,2
	    do k = 3,4
	      do l = 3,4
		do m = 3,3
		  do n = 3,3
 
		  Ml0 = 0d0
		  Ml1 = MSf(i,k,m)

		  if(b.LE.4)then
		    Ml2 = MNeu(b)
		  else
		    if(b.EQ.5) then
		      Ml2 = MCha(1)
		    else
		      Ml2 = MCha(2)
		    endif
		  endif

		  Ml3 = Mf(l,n)
! 		  
		  g0L = GS*Krondelta(l,ftt2)*Krondelta(n,fgen2)
		  g0R = GS*Krondelta(l,ftt2)*Krondelta(n,fgen2)
		  g1L = AL(b,isf1,itt1,l,igen1,n) 
		  g1R = AR(b,isf1,itt1,l,igen1,n)
		  g2 = VSfSf(1,i,isf2,k,itt2,m,igen2)
		  g3L = AL(b,i,k,ftt1,m,fgen1) 
		  g3R = AR(b,i,k,ftt1,m,fgen1)
		  
! 		
		  call myCget(muu3**2*sqrtS**2,ured*sqrtS**2,muu2*sqrtS**2,Ml1**2,Ml2**2,Ml3**2)
		  call myDget(muu2**2*sqrtS**2,muu3**2*sqrtS**2,muu1**2*sqrtS**2,muu4**2*sqrtS**2,ured*sqrtS**2,tred*sqrtS**2, Ml0**2,Ml1**2,Ml2**2,Ml3**2)
		  		  
	colori = Nc*CF

	vbA = vbA + colori*0d0
	
	vbBL = vbBL + colori*( g2*(-(g0L*(g1R*g3R*Ml2*Ml3*muu3*sqrtS*(2*xd0 + xd1) + 
     -         g1L*g3R*muu3*sqrtS*(xc0 - xc1 - 2*xc2 + 2*Ml0**2*xd0 + 2*muu4**2*sqrtS**2*xd0 + Ml0**2*xd1 + 
     -            2*muu2**2*sqrtS**2*xd1 + 3*muu4**2*sqrtS**2*xd1 - 2*sqrtS**2*tred*xd1 + muu2**2*sqrtS**2*xd11 + 
     -            muu4**2*sqrtS**2*xd11 - sqrtS**2*tred*xd11 - sqrtS**2*xd12 + muu2**2*sqrtS**2*xd12 + 
     -            muu3**2*sqrtS**2*xd12 + 2*muu4**2*sqrtS**2*xd12 - sqrtS**2*tred*xd12 + 2*muu4**2*sqrtS**2*xd13 - 
     -            2*sqrtS**2*xd2 + 2*muu1**2*sqrtS**2*xd2 + 2*muu2**2*sqrtS**2*xd2 + 2*muu3**2*sqrtS**2*xd2 + 
     -            4*muu4**2*sqrtS**2*xd2 - 2*sqrtS**2*tred*xd2 + 4*muu4**2*sqrtS**2*xd3) - 
     -         g1R*g3L*Ml3*(xc0 + Ml0**2*xd0 + 
     -            sqrtS**2*(-xd2 + muu1**2*xd2 + muu4**2*xd2 - tred*xd2 - muu3**2*(2*xd0 + xd1 + xd2) + 
     -               muu4**2*xd3 - 2*tred*xd3 + muu2**2*(2*xd0 + 3*xd1 + 3*xd2 + 2*xd3))) - 
     -         g1L*g3L*Ml2*(xc0 + Ml0**2*xd0 + 
     -            sqrtS**2*(-xd2 - muu3**2*xd2 + muu4**2*xd2 + muu4**2*xd3 - tred*(2*xd0 + xd1 + xd2 + 2*xd3) + 
     -               muu2**2*(2*xd0 + 3*xd1 + 3*xd2 + 2*xd3))))) + 
     -    g0R*muu4*sqrtS*(g1L*g3L*Ml2*Ml3*(2*xd0 + xd1 + xd2 + xd3) + 
     -       g1R*g3R*Ml2*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3)) + 
     -       g1L*g3R*Ml3*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3)) + 
     -       g1R*g3L*(xc0 - xc2 + 2*Ml0**2*xd0 + 2*muu3**2*sqrtS**2*xd0 + 2*muu4**2*sqrtS**2*xd0 - 
     -          2*sqrtS**2*tred*xd0 + 2*xd00 + Ml0**2*xd1 + 2*muu2**2*sqrtS**2*xd1 + muu3**2*sqrtS**2*xd1 + 
     -          3*muu4**2*sqrtS**2*xd1 - 3*sqrtS**2*tred*xd1 + muu2**2*sqrtS**2*xd11 + muu4**2*sqrtS**2*xd11 - 
     -          sqrtS**2*tred*xd11 - sqrtS**2*xd12 + 2*muu2**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 + 
     -          3*muu4**2*sqrtS**2*xd12 - 2*sqrtS**2*tred*xd12 + muu2**2*sqrtS**2*xd13 + 3*muu4**2*sqrtS**2*xd13 - 
     -          sqrtS**2*tred*xd13 + Ml0**2*xd2 - 2*sqrtS**2*xd2 + muu1**2*sqrtS**2*xd2 + 2*muu2**2*sqrtS**2*xd2 + 
     -          4*muu3**2*sqrtS**2*xd2 + 5*muu4**2*sqrtS**2*xd2 - 4*sqrtS**2*tred*xd2 - sqrtS**2*xd22 + 
     -          muu2**2*sqrtS**2*xd22 + muu3**2*sqrtS**2*xd22 + 2*muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 - 
     -          sqrtS**2*xd23 + muu2**2*sqrtS**2*xd23 + muu3**2*sqrtS**2*xd23 + 4*muu4**2*sqrtS**2*xd23 - 
     -          sqrtS**2*tred*xd23 + Ml0**2*xd3 + 2*muu3**2*sqrtS**2*xd3 + 5*muu4**2*sqrtS**2*xd3 - 
     -          2*sqrtS**2*tred*xd3 + 2*muu4**2*sqrtS**2*xd33))))
          
	vbBR = vbBR + colori*(g2*(-(g0R*(g1L*g3L*Ml2*Ml3*muu3*sqrtS*(2*xd0 + xd1) + 
     -         g1R*g3L*muu3*sqrtS*(xc0 - xc1 - 2*xc2 + 2*Ml0**2*xd0 + 2*muu4**2*sqrtS**2*xd0 + Ml0**2*xd1 + 
     -            2*muu2**2*sqrtS**2*xd1 + 3*muu4**2*sqrtS**2*xd1 - 2*sqrtS**2*tred*xd1 + muu2**2*sqrtS**2*xd11 + 
     -            muu4**2*sqrtS**2*xd11 - sqrtS**2*tred*xd11 - sqrtS**2*xd12 + muu2**2*sqrtS**2*xd12 + 
     -            muu3**2*sqrtS**2*xd12 + 2*muu4**2*sqrtS**2*xd12 - sqrtS**2*tred*xd12 + 2*muu4**2*sqrtS**2*xd13 - 
     -            2*sqrtS**2*xd2 + 2*muu1**2*sqrtS**2*xd2 + 2*muu2**2*sqrtS**2*xd2 + 2*muu3**2*sqrtS**2*xd2 + 
     -            4*muu4**2*sqrtS**2*xd2 - 2*sqrtS**2*tred*xd2 + 4*muu4**2*sqrtS**2*xd3) - 
     -         g1L*g3R*Ml3*(xc0 + Ml0**2*xd0 + 
     -            sqrtS**2*(-xd2 + muu1**2*xd2 + muu4**2*xd2 - tred*xd2 - muu3**2*(2*xd0 + xd1 + xd2) + 
     -               muu4**2*xd3 - 2*tred*xd3 + muu2**2*(2*xd0 + 3*xd1 + 3*xd2 + 2*xd3))) - 
     -         g1R*g3R*Ml2*(xc0 + Ml0**2*xd0 + 
     -            sqrtS**2*(-xd2 - muu3**2*xd2 + muu4**2*xd2 + muu4**2*xd3 - tred*(2*xd0 + xd1 + xd2 + 2*xd3) + 
     -               muu2**2*(2*xd0 + 3*xd1 + 3*xd2 + 2*xd3))))) + 
     -    g0L*muu4*sqrtS*(g1R*g3R*Ml2*Ml3*(2*xd0 + xd1 + xd2 + xd3) + 
     -       g1L*g3L*Ml2*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3)) + 
     -       g1R*g3L*Ml3*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3)) + 
     -       g1L*g3R*(xc0 - xc2 + 2*Ml0**2*xd0 + 2*muu3**2*sqrtS**2*xd0 + 2*muu4**2*sqrtS**2*xd0 - 
     -          2*sqrtS**2*tred*xd0 + 2*xd00 + Ml0**2*xd1 + 2*muu2**2*sqrtS**2*xd1 + muu3**2*sqrtS**2*xd1 + 
     -          3*muu4**2*sqrtS**2*xd1 - 3*sqrtS**2*tred*xd1 + muu2**2*sqrtS**2*xd11 + muu4**2*sqrtS**2*xd11 - 
     -          sqrtS**2*tred*xd11 - sqrtS**2*xd12 + 2*muu2**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 + 
     -          3*muu4**2*sqrtS**2*xd12 - 2*sqrtS**2*tred*xd12 + muu2**2*sqrtS**2*xd13 + 3*muu4**2*sqrtS**2*xd13 - 
     -          sqrtS**2*tred*xd13 + Ml0**2*xd2 - 2*sqrtS**2*xd2 + muu1**2*sqrtS**2*xd2 + 2*muu2**2*sqrtS**2*xd2 + 
     -          4*muu3**2*sqrtS**2*xd2 + 5*muu4**2*sqrtS**2*xd2 - 4*sqrtS**2*tred*xd2 - sqrtS**2*xd22 + 
     -          muu2**2*sqrtS**2*xd22 + muu3**2*sqrtS**2*xd22 + 2*muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 - 
     -          sqrtS**2*xd23 + muu2**2*sqrtS**2*xd23 + muu3**2*sqrtS**2*xd23 + 4*muu4**2*sqrtS**2*xd23 - 
     -          sqrtS**2*tred*xd23 + Ml0**2*xd3 + 2*muu3**2*sqrtS**2*xd3 + 5*muu4**2*sqrtS**2*xd3 - 
     -          2*sqrtS**2*tred*xd3 + 2*muu4**2*sqrtS**2*xd33))))
     
	vbCLp1 = vbCLp1 + colori*(g2*(g0L*(2*g1L*g3L*Ml2*muu3*sqrtS*xd2 + g1R*g3R*Ml2*Ml3*(2*xd0 + xd1 + xd2) + 
     -       g1R*g3L*Ml3*muu3*sqrtS*(2*xd0 + xd1 + 2*xd2) + 
     -       g1L*g3R*(2*xc0 - xc2 + 3*Ml0**2*xd0 + 2*muu2**2*sqrtS**2*xd0 + 2*muu4**2*sqrtS**2*xd0 - 
     -          2*sqrtS**2*tred*xd0 + Ml0**2*xd1 + 5*muu2**2*sqrtS**2*xd1 + 3*muu4**2*sqrtS**2*xd1 - 
     -          3*sqrtS**2*tred*xd1 + muu2**2*sqrtS**2*xd11 + muu4**2*sqrtS**2*xd11 - sqrtS**2*tred*xd11 - 
     -          sqrtS**2*xd12 + 2*muu2**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 + 3*muu4**2*sqrtS**2*xd12 - 
     -          2*sqrtS**2*tred*xd12 + 2*muu4**2*sqrtS**2*xd13 + Ml0**2*xd2 - 3*sqrtS**2*xd2 + 
     -          2*muu1**2*sqrtS**2*xd2 + 5*muu2**2*sqrtS**2*xd2 + 3*muu3**2*sqrtS**2*xd2 + 6*muu4**2*sqrtS**2*xd2 - 
     -          5*sqrtS**2*tred*xd2 - sqrtS**2*xd22 + muu2**2*sqrtS**2*xd22 + muu3**2*sqrtS**2*xd22 + 
     -          2*muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 + 2*muu4**2*sqrtS**2*xd23 + 2*muu2**2*sqrtS**2*xd3 + 
     -          5*muu4**2*sqrtS**2*xd3 - 2*sqrtS**2*tred*xd3)) - 
     -    g0R*muu4*sqrtS*(g1L*g3R*Ml3*xd3 + g1R*g3R*Ml2*(2*xd0 + xd1 + xd2 + 2*xd3) + 
     -       g1R*g3L*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3)))))
     
	vbCRp1 = vbCRp1 + colori*(g2*(g0R*(2*g1R*g3R*Ml2*muu3*sqrtS*xd2 + g1L*g3L*Ml2*Ml3*(2*xd0 + xd1 + xd2) + 
     -       g1L*g3R*Ml3*muu3*sqrtS*(2*xd0 + xd1 + 2*xd2) + 
     -       g1R*g3L*(2*xc0 - xc2 + 3*Ml0**2*xd0 + 2*muu2**2*sqrtS**2*xd0 + 2*muu4**2*sqrtS**2*xd0 - 
     -          2*sqrtS**2*tred*xd0 + Ml0**2*xd1 + 5*muu2**2*sqrtS**2*xd1 + 3*muu4**2*sqrtS**2*xd1 - 
     -          3*sqrtS**2*tred*xd1 + muu2**2*sqrtS**2*xd11 + muu4**2*sqrtS**2*xd11 - sqrtS**2*tred*xd11 - 
     -          sqrtS**2*xd12 + 2*muu2**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 + 3*muu4**2*sqrtS**2*xd12 - 
     -          2*sqrtS**2*tred*xd12 + 2*muu4**2*sqrtS**2*xd13 + Ml0**2*xd2 - 3*sqrtS**2*xd2 + 
     -          2*muu1**2*sqrtS**2*xd2 + 5*muu2**2*sqrtS**2*xd2 + 3*muu3**2*sqrtS**2*xd2 + 6*muu4**2*sqrtS**2*xd2 - 
     -          5*sqrtS**2*tred*xd2 - sqrtS**2*xd22 + muu2**2*sqrtS**2*xd22 + muu3**2*sqrtS**2*xd22 + 
     -          2*muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 + 2*muu4**2*sqrtS**2*xd23 + 2*muu2**2*sqrtS**2*xd3 + 
     -          5*muu4**2*sqrtS**2*xd3 - 2*sqrtS**2*tred*xd3)) - 
     -    g0L*muu4*sqrtS*(g1R*g3L*Ml3*xd3 + g1L*g3L*Ml2*(2*xd0 + xd1 + xd2 + 2*xd3) + 
     -       g1L*g3R*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3)))))
     

      vbCLp2 = vbCLp2 + colori*0d0
	  vbCRp2 = vbCRp2 + colori*0d0
	  vbDp1 = vbDp1 + colori*0d0
	  vbDk1 = vBDk1 + colori*0d0
	  vbDk2 = vbDk2 + colori*0d0
	
		   
		  enddo
		enddo
	      enddo
	    enddo
	enddo



! ! 
! ! ! ****************************************************
! ! ! ************ Box 1.3 Neutralino-Squark (i,k,m)- Gluon- Quark (l,n) ***********************
! ! ! ****************************************************
! ! ! Coupling 0 L<->R , Coupling 1 gains minus sign due to charge flow change

	do j = 1,2
	    do k = 3,4
	      do l = 3,4
		do m = 3,3
		  do n = 3,3
		  
		  Ml0 = 0d0
		  Ml1 = MSf(j,l,n)
		  
		  if(b.LE.4)then
		    Ml2 = MNeu(b)
		  else
		    if(b.EQ.5) then
		      Ml2 = MCha(1)
		    else
		      Ml2 = MCha(2)
		    endif
		  endif
		  
		  Ml3 = Mf(k,m)
		  
		  g0L = AL(b,j,l,ftt2,n,fgen2) 
		  g0R = AR(b,j,l,ftt2,n,fgen2)
		  g1 = -VSfSf(1,isf1,j,itt1,l,igen1,n)
		  g2L = AL(b,isf2,itt2,k,igen2,m) 
		  g2R = AR(b,isf2,itt2,k,igen2,m)
		  g3L = -GS*Krondelta(k,ftt1)*Krondelta(m,fgen1)
		  g3R = -GS*Krondelta(k,ftt1)*Krondelta(m,fgen1)
! 		
		  call myCget(muu3**2*sqrtS**2,ured*sqrtS**2,muu2*sqrtS**2,Ml1**2,Ml2**2,Ml3**2)
		  call myDget(muu1**2*sqrtS**2,muu4**2*sqrtS**2,muu2**2*sqrtS**2,muu3**2*sqrtS**2,ured*sqrtS**2,tred*sqrtS**2, Ml0**2,Ml1**2,Ml2**2,Ml3**2)
		 		  
	colori = Nc*CF

	vbA = vbA + colori*0d0
	
	vbBL = vbBL + colori*(g1*(g0R*muu4*sqrtS*(g2R*Ml2*(g3R*Ml3 - g3L*muu3*sqrtS)*xd2 + g2L*g3L*Ml3*muu3*sqrtS*xd3 + 
     -       g2L*g3R*(2*xc0 + xc1 + xc2 + Ml0**2*xd0 + 2*muu1**2*sqrtS**2*xd0 + 2*muu3**2*sqrtS**2*xd0 - 
     -          2*sqrtS**2*tred*xd0 - Ml0**2*xd1 + 3*muu1**2*sqrtS**2*xd1 + muu3**2*sqrtS**2*xd1 - 
     -          sqrtS**2*tred*xd1 + muu1**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 - sqrtS**2*tred*xd12 - 
     -          sqrtS**2*xd2 + 3*muu1**2*sqrtS**2*xd2 + 4*muu3**2*sqrtS**2*xd2 + muu4**2*sqrtS**2*xd2 - 
     -          3*sqrtS**2*tred*xd2 - sqrtS**2*xd22 + muu1**2*sqrtS**2*xd22 + 2*muu3**2*sqrtS**2*xd22 + 
     -          muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 + 2*muu3**2*sqrtS**2*xd23 + 2*muu1**2*sqrtS**2*xd3 + 
     -          3*muu3**2*sqrtS**2*xd3 - 2*sqrtS**2*tred*xd3)) + 
     -    g0L*(g2L*g3L*Ml2*Ml3*muu3*sqrtS*xd3 + 
     -       g2R*g3R*Ml3*(xc0 + Ml0**2*xd0 + sqrtS**2*
     -           ((-1 + muu2**2 + muu3**2 + muu4**2 - tred)*xd2 + 2*(muu3**2 - tred)*xd3 + 
     -             muu1**2*(2*xd0 + 3*xd1 + 3*xd2 + 2*xd3))) + 
     -       g2L*g3R*Ml2*(xc0 + Ml0**2*xd0 + sqrtS**2*
     -           (-2*tred*xd0 - tred*xd1 - xd2 + muu4**2*xd2 - tred*xd2 - 2*tred*xd3 + 
     -             muu1**2*(2*xd0 + 3*xd1 + 3*xd2 + 2*xd3) + muu3**2*(2*xd0 + xd1 + 2*xd2 + 3*xd3))) + 
     -       g2R*g3L*muu3*sqrtS*(-2*xc0 - 2*muu1**2*sqrtS**2*xd0 + 2*xd00 - 3*muu1**2*sqrtS**2*xd1 + 
     -          muu1**2*sqrtS**2*xd13 + muu3**2*sqrtS**2*xd13 - sqrtS**2*tred*xd13 + sqrtS**2*xd2 - 
     -          3*muu1**2*sqrtS**2*xd2 - muu2**2*sqrtS**2*xd2 - muu3**2*sqrtS**2*xd2 - muu4**2*sqrtS**2*xd2 + 
     -          sqrtS**2*tred*xd2 - sqrtS**2*xd23 + muu1**2*sqrtS**2*xd23 + 2*muu3**2*sqrtS**2*xd23 + 
     -          muu4**2*sqrtS**2*xd23 - sqrtS**2*tred*xd23 - 2*muu1**2*sqrtS**2*xd3 + 
     -          Ml0**2*(-xd0 + xd1 + xd2 + xd3) + 2*muu3**2*sqrtS**2*xd33))))
     
	vbBR = vbBR + colori*(g1*(g0L*muu4*sqrtS*(g2L*Ml2*(g3L*Ml3 - g3R*muu3*sqrtS)*xd2 + g2R*g3R*Ml3*muu3*sqrtS*xd3 + 
     -       g2R*g3L*(2*xc0 + xc1 + xc2 + Ml0**2*xd0 + 2*muu1**2*sqrtS**2*xd0 + 2*muu3**2*sqrtS**2*xd0 - 
     -          2*sqrtS**2*tred*xd0 - Ml0**2*xd1 + 3*muu1**2*sqrtS**2*xd1 + muu3**2*sqrtS**2*xd1 - 
     -          sqrtS**2*tred*xd1 + muu1**2*sqrtS**2*xd12 + muu3**2*sqrtS**2*xd12 - sqrtS**2*tred*xd12 - 
     -          sqrtS**2*xd2 + 3*muu1**2*sqrtS**2*xd2 + 4*muu3**2*sqrtS**2*xd2 + muu4**2*sqrtS**2*xd2 - 
     -          3*sqrtS**2*tred*xd2 - sqrtS**2*xd22 + muu1**2*sqrtS**2*xd22 + 2*muu3**2*sqrtS**2*xd22 + 
     -          muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 + 2*muu3**2*sqrtS**2*xd23 + 2*muu1**2*sqrtS**2*xd3 + 
     -          3*muu3**2*sqrtS**2*xd3 - 2*sqrtS**2*tred*xd3)) + 
     -    g0R*(g2R*g3R*Ml2*Ml3*muu3*sqrtS*xd3 + 
     -       g2L*g3L*Ml3*(xc0 + Ml0**2*xd0 + sqrtS**2*
     -           ((-1 + muu2**2 + muu3**2 + muu4**2 - tred)*xd2 + 2*(muu3**2 - tred)*xd3 + 
     -             muu1**2*(2*xd0 + 3*xd1 + 3*xd2 + 2*xd3))) + 
     -       g2R*g3L*Ml2*(xc0 + Ml0**2*xd0 + sqrtS**2*
     -           (-2*tred*xd0 - tred*xd1 - xd2 + muu4**2*xd2 - tred*xd2 - 2*tred*xd3 + 
     -             muu1**2*(2*xd0 + 3*xd1 + 3*xd2 + 2*xd3) + muu3**2*(2*xd0 + xd1 + 2*xd2 + 3*xd3))) + 
     -       g2L*g3R*muu3*sqrtS*(-2*xc0 - 2*muu1**2*sqrtS**2*xd0 + 2*xd00 - 3*muu1**2*sqrtS**2*xd1 + 
     -          muu1**2*sqrtS**2*xd13 + muu3**2*sqrtS**2*xd13 - sqrtS**2*tred*xd13 + sqrtS**2*xd2 - 
     -          3*muu1**2*sqrtS**2*xd2 - muu2**2*sqrtS**2*xd2 - muu3**2*sqrtS**2*xd2 - muu4**2*sqrtS**2*xd2 + 
     -          sqrtS**2*tred*xd2 - sqrtS**2*xd23 + muu1**2*sqrtS**2*xd23 + 2*muu3**2*sqrtS**2*xd23 + 
     -          muu4**2*sqrtS**2*xd23 - sqrtS**2*tred*xd23 - 2*muu1**2*sqrtS**2*xd3 + 
     -          Ml0**2*(-xd0 + xd1 + xd2 + xd3) + 2*muu3**2*sqrtS**2*xd33))))
     
	vbCLp1 = vbCLp1 + colori*( g1*(g0L*(g2L*g3L*Ml2*Ml3*(2*xd0 + xd1 + xd2) - g2R*g3R*Ml3*muu3*sqrtS*xd3 - 
     -       g2L*g3R*Ml2*muu3*sqrtS*(2*xd0 + xd1 + xd2 + 2*xd3) + 
     -       g2R*g3L*(4*xc0 - xc2 + 3*Ml0**2*xd0 + 2*muu1**2*sqrtS**2*xd0 + 2*muu3**2*sqrtS**2*xd0 - 
     -          2*sqrtS**2*tred*xd0 - Ml0**2*xd1 + 5*muu1**2*sqrtS**2*xd1 + 3*muu3**2*sqrtS**2*xd1 - 
     -          3*sqrtS**2*tred*xd1 + muu1**2*sqrtS**2*xd11 + muu3**2*sqrtS**2*xd11 - sqrtS**2*tred*xd11 - 
     -          sqrtS**2*xd12 + 2*muu1**2*sqrtS**2*xd12 + 3*muu3**2*sqrtS**2*xd12 + muu4**2*sqrtS**2*xd12 - 
     -          2*sqrtS**2*tred*xd12 + 2*muu3**2*sqrtS**2*xd13 - Ml0**2*xd2 - 3*sqrtS**2*xd2 + 
     -          5*muu1**2*sqrtS**2*xd2 + 2*muu2**2*sqrtS**2*xd2 + 6*muu3**2*sqrtS**2*xd2 + 3*muu4**2*sqrtS**2*xd2 - 
     -          5*sqrtS**2*tred*xd2 - sqrtS**2*xd22 + muu1**2*sqrtS**2*xd22 + 2*muu3**2*sqrtS**2*xd22 + 
     -          muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 + 2*muu3**2*sqrtS**2*xd23 + 2*muu1**2*sqrtS**2*xd3 + 
     -          5*muu3**2*sqrtS**2*xd3 - 2*sqrtS**2*tred*xd3)) + 
     -    g0R*muu4*sqrtS*(2*g2R*g3L*Ml2*xd2 + g2L*g3L*Ml3*(2*xd0 + xd1 + 2*xd2) - 
     -       g2L*g3R*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3)))))
     
	vbCRp1 = vbCRp1 + colori*(g1*(g0R*(g2R*g3R*Ml2*Ml3*(2*xd0 + xd1 + xd2) - g2L*g3L*Ml3*muu3*sqrtS*xd3 - 
     -       g2R*g3L*Ml2*muu3*sqrtS*(2*xd0 + xd1 + xd2 + 2*xd3) + 
     -       g2L*g3R*(4*xc0 - xc2 + 3*Ml0**2*xd0 + 2*muu1**2*sqrtS**2*xd0 + 2*muu3**2*sqrtS**2*xd0 - 
     -          2*sqrtS**2*tred*xd0 - Ml0**2*xd1 + 5*muu1**2*sqrtS**2*xd1 + 3*muu3**2*sqrtS**2*xd1 - 
     -          3*sqrtS**2*tred*xd1 + muu1**2*sqrtS**2*xd11 + muu3**2*sqrtS**2*xd11 - sqrtS**2*tred*xd11 - 
     -          sqrtS**2*xd12 + 2*muu1**2*sqrtS**2*xd12 + 3*muu3**2*sqrtS**2*xd12 + muu4**2*sqrtS**2*xd12 - 
     -          2*sqrtS**2*tred*xd12 + 2*muu3**2*sqrtS**2*xd13 - Ml0**2*xd2 - 3*sqrtS**2*xd2 + 
     -          5*muu1**2*sqrtS**2*xd2 + 2*muu2**2*sqrtS**2*xd2 + 6*muu3**2*sqrtS**2*xd2 + 3*muu4**2*sqrtS**2*xd2 - 
     -          5*sqrtS**2*tred*xd2 - sqrtS**2*xd22 + muu1**2*sqrtS**2*xd22 + 2*muu3**2*sqrtS**2*xd22 + 
     -          muu4**2*sqrtS**2*xd22 - sqrtS**2*tred*xd22 + 2*muu3**2*sqrtS**2*xd23 + 2*muu1**2*sqrtS**2*xd3 + 
     -          5*muu3**2*sqrtS**2*xd3 - 2*sqrtS**2*tred*xd3)) + 
     -    g0L*muu4*sqrtS*(2*g2L*g3R*Ml2*xd2 + g2R*g3R*Ml3*(2*xd0 + xd1 + 2*xd2) - 
     -       g2R*g3L*muu3*sqrtS*(2*xd0 + xd1 + 2*(xd2 + xd3)))))
     
!      

      vbCLp2 = vbCLp2 + colori*0d0
	  vbCRp2 = vbCRp2 + colori*0d0
	  vbDp1 = vbDp1 + colori*0d0
	  vbDk1 = vBDk1 + colori*0d0
	  vbDk2 = vbDk2 + colori*0d0
	
		  enddo
		enddo
	      enddo
	    enddo
	enddo
 

! ! ! ! 
! ! ! ! ****************************************************
! ! ! ! ! ************ Box 4 Neutralino-Quark(k,m)-Gluon- Quark(l,n) ***********************
! ! ! ! ! ****************************************************
! !! ! ! Coupling 1 L<->R , Coupling 0 gains minus sign due to fermion flow change
 	
	  do k = 3,4
	    do l = 3,4
	      do m = 3,3
		do n = 3,3
		  
		
		  if(b.LE.4)then
		    Ml0 = MNeu(b)
		  else
		    if(b.EQ.5) then
		      Ml0 = MCha(1)
		    else
		      Ml0 = MCha(2)
		    endif
		  endif
		  Ml1 = Mf(k,m)
		  Ml2 = 0
		  Ml3 = Mf(l,n)
		  
		  g0 = GS*Krondelta(l,ftt2)*Krondelta(n,fgen2)
		  g1L = AL(b,isf1,itt1,l,igen1,n)
		  g1R = AR(b,isf1,itt1,l,igen1,n)
		  g2L = AL(b,isf2,itt2,k,igen2,m)
		  g2R = AR(b,isf2,itt2,k,igen2,m)
		  g3 = -GS*Krondelta(k,ftt1)*Krondelta(m,fgen1)
		  
		  call myCget(muu3**2*sqrtS**2,muu4**2*sqrtS**2,sqrtS**2,Ml1**2,Ml2**2,Ml3**2)
		  call myDget(muu2**2*sqrtS**2,muu3**2*sqrtS**2,muu4**2*sqrtS**2,muu1**2*sqrtS**2,ured*sqrtS**2,sqrtS**2, Ml0**2,Ml1**2,Ml2**2,Ml3**2)
		
	colori = Nc*CF

	vbA = vbA + colori*0d0
	
	vbBL = vbBL + colori*(-2*g0*g3*(g1R*g2R*Ml0*(2*Ml1*Ml3*xd0 - Ml3*muu3*sqrtS*(xd0 + xd1) + Ml1*muu4*sqrtS*(xd1 + xd2)) - 
     -    g1L*g2R*muu3*sqrtS*(-xc1 - xc2 + Ml0**2*xd0 + Ml0**2*xd1 + muu1**2*sqrtS**2*xd2 + muu1**2*sqrtS**2*xd3) + 
     -    g1L*g2R*Ml1*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*
     -        (Ml3*muu4*(xd1 + xd2) + sqrtS*((-1 + muu1**2 + muu2**2)*xd1 + 
     -             (-1 + 2*muu1**2 + muu2**2 + muu3**2 - tred)*xd2 + 2*muu1**2*xd3))) + 
     -    g1L*g2L*Ml0*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*
     -        (-(Ml1*muu3*xd1) + Ml3*muu4*(xd0 + xd1 + xd2) + 
     -          sqrtS*((-1 + muu1**2 + muu2**2)*xd0 + (-1 + muu1**2 + 3*muu2**2)*xd1 - 2*xd2 + 3*muu1**2*xd2 + 
     -             3*muu2**2*xd2 + muu3**2*xd2 + muu4**2*xd2 - 2*tred*xd2 - xd3 + 3*muu1**2*xd3 + muu2**2*xd3))) + 
     -    g1R*g2L*(muu4*sqrtS*(-xc2 + Ml0**2*(xd0 + xd1 + xd2) + 
     -          sqrtS**2*((muu3**2 - tred)*xd2 + muu1**2*(xd2 + xd3))) + 
     -       Ml3*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*(-(Ml1*muu3*xd1) + 
     -             sqrtS*((-1 + muu1**2 + muu4**2 - tred)*xd2 + (-1 + muu1**2)*xd3 + muu2**2*(2*xd1 + 2*xd2 + xd3)))))))
     
	vbBR = vbBR + colori*(-2*g0*g3*(g1L*g2L*Ml0*(2*Ml1*Ml3*xd0 - Ml3*muu3*sqrtS*(xd0 + xd1) + Ml1*muu4*sqrtS*(xd1 + xd2)) - 
     -    g1R*g2L*muu3*sqrtS*(-xc1 - xc2 + Ml0**2*xd0 + Ml0**2*xd1 + muu1**2*sqrtS**2*xd2 + muu1**2*sqrtS**2*xd3) + 
     -    g1R*g2L*Ml1*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*
     -        (Ml3*muu4*(xd1 + xd2) + sqrtS*((-1 + muu1**2 + muu2**2)*xd1 + 
     -             (-1 + 2*muu1**2 + muu2**2 + muu3**2 - tred)*xd2 + 2*muu1**2*xd3))) + 
     -    g1R*g2R*Ml0*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*
     -        (-(Ml1*muu3*xd1) + Ml3*muu4*(xd0 + xd1 + xd2) + 
     -          sqrtS*((-1 + muu1**2 + muu2**2)*xd0 + (-1 + muu1**2 + 3*muu2**2)*xd1 - 2*xd2 + 3*muu1**2*xd2 + 
     -             3*muu2**2*xd2 + muu3**2*xd2 + muu4**2*xd2 - 2*tred*xd2 - xd3 + 3*muu1**2*xd3 + muu2**2*xd3))) + 
     -    g1L*g2R*(muu4*sqrtS*(-xc2 + Ml0**2*(xd0 + xd1 + xd2) + 
     -          sqrtS**2*((muu3**2 - tred)*xd2 + muu1**2*(xd2 + xd3))) + 
     -       Ml3*(2*xc0 + 2*Ml0**2*xd0 + sqrtS*(-(Ml1*muu3*xd1) + 
     -             sqrtS*((-1 + muu1**2 + muu4**2 - tred)*xd2 + (-1 + muu1**2)*xd3 + muu2**2*(2*xd1 + 2*xd2 + xd3)))))))
	
	vbCLp1 = vbCLp1 + colori*( -2*g0*g3*(-(g1R*g2L*muu3*muu4*sqrtS**2*xd2) + g1R*g2L*Ml1*Ml3*(xd1 + xd2 + xd3) + 
     -    g1L*g2L*Ml0*Ml1*(xd0 + xd1 + xd2 + xd3) + g1R*g2R*Ml0*Ml3*(xd0 + xd1 + xd2 + xd3) + 
     -    g1L*g2R*(xc0 + Ml0**2*(2*xd0 + xd1 + xd2 + xd3) + 
     -       sqrtS**2*(-(tred*xd2) + muu2**2*(xd1 + xd2) + muu1**2*(xd2 + xd3)))))
	
	vbCRp1 = vbCRp1 + colori*( -2*g0*g3*(-(g1L*g2R*muu3*muu4*sqrtS**2*xd2) + g1L*g2R*Ml1*Ml3*(xd1 + xd2 + xd3) + 
     -    g1R*g2R*Ml0*Ml1*(xd0 + xd1 + xd2 + xd3) + g1L*g2L*Ml0*Ml3*(xd0 + xd1 + xd2 + xd3) + 
     -    g1R*g2L*(xc0 + Ml0**2*(2*xd0 + xd1 + xd2 + xd3) + 
     -       sqrtS**2*(-(tred*xd2) + muu2**2*(xd1 + xd2) + muu1**2*(xd2 + xd3)))))
	

      vbCLp2 = vbCLp2 + colori*0d0
	  vbCRp2 = vbCRp2 + colori*0d0
	  vbDp1 = vbDp1 + colori*0d0
	  vbDk1 = vBDk1 + colori*0d0
	  vbDk2 = vbDk2 + colori*0d0
	
		  enddo
		enddo
	      enddo
	    enddo
! 	    
! ! ! ! ! 	    
	  enddo		!End of neutralino-loop	
! ! ! ! ! 

! ! ! !-------------------------------------------------------------------------------------------------
! ! ! !-----------------Hier /(4PI^2) mitnehmen---------------------------------------------------------
! ! ! !-------------------------------------------------------------------------------------------------
! ! ! 	  
      vbA = vbA/(4d0*Pi)**2
	  vbBL = vbBL/(4d0*Pi)**2
	  vbBR = vbBR/(4d0*Pi)**2
	  vbCLp1 = vbCLp1/(4d0*Pi)**2
	  vbCRp1 = vbCRp1/(4d0*Pi)**2
	  vbCLp2 = vbCLp2/(4d0*Pi)**2
	  vbCRp2 = vbCRp2/(4d0*Pi)**2
	  vbDp1 = vbDp1/(4d0*Pi)**2
	  vbDk1 = vbDk1/(4d0*Pi)**2
	  vbDk2 = vbDk2/(4d0*Pi)**2

! ! 
! ! 	
! ! c--------------------------------------------------------------------------------------------------------------------------
! ! ! c----------------------------------------------Amplitudenquadrate----------------------------------------------------------
! ! ! c--------------------------------------------------------------------------------------------------------------------------
!
! c----------------------------------------------Box-T-Kanal Gluino----------------------------------------------------------
! ! 
	  g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
	  g1Lc = dconjg(g1L)
	  g1Rc = dconjg(g1R)

	
	  g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
	  g2Lc = dconjg(g2L)
	  g2Rc = dconjg(g2R)
	  
	  
	  colorf = (Nc**2-1d0)/2d0/(2*Nc*(Nc+1))
	
	  stst2QQ_M2Boxsextet = stst2QQ_M2Boxsextet + 2d0*colorf*dreal((sqrtS**3*(g1Lc*g2Rc*
     -       (muu3**2*muu4*vbA + muu3*muu4**2*vbA + muu3*tred*vbA + 
     -         muu4*tred*vbA + muu3**2*muu4*vbBL + muu4*tred*vbBL + 
     -         muu3*muu4**2*vbBR + muu3*tred*vbBR + 
     -         muu3**3*muu4*sqrtS*vbCLp1 - 
     -         muu3*muu4*sqrtS*tred*vbCLp1 - 
     -         muu3*muu4**3*sqrtS*vbCLp2 + 
     -         muu3*muu4*sqrtS*tred*vbCLp2 - sqrtS*tred*vbCRp1 + 
     -         muu3**2*sqrtS*tred*vbCRp1 - sqrtS*tred**2*vbCRp1 + 
     -         muu1**2*(-(muu4*(vbA + vbBL + muu3*sqrtS*vbCLp1)) + 
     -            muu4**2*sqrtS*vbCRp1 + sqrtS*tred*(vbCRp1 - vbCRp2))
     -           + sqrtS*tred*(1 - muu4**2 + tred)*vbCRp2 - 
     -         muu2**2*(muu3*(vbA + vbBR - muu4*sqrtS*vbCLp2) + 
     -            sqrtS*(muu1**2 - tred)*(vbCRp1 - vbCRp2) + 
     -            muu3**2*sqrtS*vbCRp2)) + 
     -      g1Rc*g2Rc*muugl*((-1 + muu3 + muu4)*(1 + muu3 + muu4)*
     -          vbA - vbBR + muu3**2*(vbBR + muu4*sqrtS*vbCLp1) + 
     -         muu4*(muu4*vbBR - muu4**2*sqrtS*vbCLp2 + 
     -            sqrtS*((muu1**2 - tred)*(vbCLp1 - vbCLp2) + vbCLp2))
     -           + muu3**3*sqrtS*vbCRp1 + 
     -         muu3*(2*muu4*vbBL + 
     -            sqrtS*(-1 + muu2**2 - tred)*vbCRp1 - 
     -            muu4**2*sqrtS*vbCRp2 + 
     -            sqrtS*(-muu2**2 + tred)*vbCRp2)) + 
     -      g1Lc*g2Lc*muugl*((-1 + muu3 + muu4)*(1 + muu3 + muu4)*
     -          vbA - vbBL + muu3**3*sqrtS*vbCLp1 + 
     -         muu3*(2*muu4*vbBR + 
     -            sqrtS*(-1 + muu2**2 - tred)*vbCLp1 - 
     -            muu4**2*sqrtS*vbCLp2 + 
     -            sqrtS*(-muu2**2 + tred)*vbCLp2) + 
     -         muu3**2*(vbBL + muu4*sqrtS*vbCRp1) + 
     -         muu4*(muu4*vbBL - muu4**2*sqrtS*vbCRp2 + 
     -            sqrtS*((muu1**2 - tred)*(vbCRp1 - vbCRp2) + vbCRp2))
     -         ) + g1Rc*g2Lc*
     -       (muu3**2*muu4*vbA + muu3*muu4**2*vbA + muu3*tred*vbA + 
     -         muu4*tred*vbA + muu3*muu4**2*vbBL + muu3*tred*vbBL + 
     -         muu3**2*muu4*vbBR + muu4*tred*vbBR - 
     -         sqrtS*tred*vbCLp1 + muu3**2*sqrtS*tred*vbCLp1 - 
     -         sqrtS*tred**2*vbCLp1 + sqrtS*tred*vbCLp2 - 
     -         muu4**2*sqrtS*tred*vbCLp2 + sqrtS*tred**2*vbCLp2 + 
     -         muu3**3*muu4*sqrtS*vbCRp1 - 
     -         muu3*muu4*sqrtS*tred*vbCRp1 + 
     -         muu1**2*(muu4**2*sqrtS*vbCLp1 + 
     -            sqrtS*tred*(vbCLp1 - vbCLp2) - 
     -            muu4*(vbA + vbBR + muu3*sqrtS*vbCRp1)) - 
     -         muu3*muu4**3*sqrtS*vbCRp2 + 
     -         muu3*muu4*sqrtS*tred*vbCRp2 - 
     -         muu2**2*(sqrtS*(muu1**2 - tred)*(vbCLp1 - vbCLp2) + 
     -            muu3**2*sqrtS*vbCLp2 + 
     -            muu3*(vbA + vbBL - muu4*sqrtS*vbCRp2)))))/
     -  dconjg(DenGlT))
     
! 
! ! ! c----------------------------------------------Box-U-Kanal Gluino----------------------------------------------------------
	  g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
	  g1Lc = dconjg(g1L)
	  g1Rc = dconjg(g1R)
	
	  g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
	  g2Lc = dconjg(g2L)
	  g2Rc = dconjg(g2R)
	  
	  colorf = (Nc**2-1d0)/2d0/(2*Nc*(Nc+1))
	
	  stst2QQ_M2Boxsextet = stst2QQ_M2Boxsextet + 2d0*colorf*dreal((sqrtS**3*(g1Lc*g2Rc*
     -       (muu4**3*(vbA + vbBL) + 
     -         muu4*(-1 + muu1**2 - tred)*(vbA + vbBL) + 
     -         muu3**3*(vbA + vbBR + muu4*sqrtS*vbCLp1) + 
     -         muu3*((-1 + muu2**2 + 2*muu4**2 - tred)*vbA + 
     -            (-1 + muu2**2 + 2*muu4**2 - tred)*vbBR + 
     -            muu4*sqrtS*
     -             ((-1 + 2*muu1**2 + muu2**2 - tred)*vbCLp1 + 
     -               vbCLp2 - 
     -               (muu1**2 + 2*muu2**2 + muu4**2 - tred)*vbCLp2)) +
     -           muu3**4*sqrtS*vbCRp1 + 
     -         sqrtS*((muu1**2 - tred)*(muu2**2 - tred) + tred)*
     -          (vbCRp1 - vbCRp2) - muu4**4*sqrtS*vbCRp2 + 
     -         muu4**2*sqrtS*
     -          (-(tred*vbCRp1) + vbCRp2 - muu1**2*vbCRp2 + 
     -            2*tred*vbCRp2) + 
     -         muu3**2*(2*muu4*(vbA + vbBL) + 
     -            sqrtS*(-1 + muu2**2 - 2*tred)*vbCRp1 + 
     -            muu4**2*sqrtS*(vbCRp1 - vbCRp2) + sqrtS*tred*vbCRp2)
     -         ) + g1Rc*g2Rc*muugl*
     -       ((-1 + muu3 + muu4)*(1 + muu3 + muu4)*vbA - vbBR + 
     -         muu3**2*(vbBR + muu4*sqrtS*vbCLp1) + 
     -         muu4*(muu4*vbBR - muu4**2*sqrtS*vbCLp2 + 
     -            sqrtS*((muu1**2 - tred)*(vbCLp1 - vbCLp2) + vbCLp2))
     -           + muu3**3*sqrtS*vbCRp1 + 
     -         muu3*(2*muu4*vbBL + 
     -            sqrtS*(-1 + muu2**2 - tred)*vbCRp1 - 
     -            muu4**2*sqrtS*vbCRp2 + 
     -            sqrtS*(-muu2**2 + tred)*vbCRp2)) + 
     -      g1Lc*g2Lc*muugl*((-1 + muu3 + muu4)*(1 + muu3 + muu4)*
     -          vbA - vbBL + muu3**3*sqrtS*vbCLp1 + 
     -         muu3*(2*muu4*vbBR + 
     -            sqrtS*(-1 + muu2**2 - tred)*vbCLp1 - 
     -            muu4**2*sqrtS*vbCLp2 + 
     -            sqrtS*(-muu2**2 + tred)*vbCLp2) + 
     -         muu3**2*(vbBL + muu4*sqrtS*vbCRp1) + 
     -         muu4*(muu4*vbBL - muu4**2*sqrtS*vbCRp2 + 
     -            sqrtS*((muu1**2 - tred)*(vbCRp1 - vbCRp2) + vbCRp2))
     -         ) + g1Rc*g2Lc*
     -       (muu4**3*(vbA + vbBR) + 
     -         muu4*(-1 + muu1**2 - tred)*(vbA + vbBR) + 
     -         muu3**4*sqrtS*vbCLp1 + 
     -         sqrtS*((muu1**2 - tred)*(muu2**2 - tred) + tred)*
     -          (vbCLp1 - vbCLp2) - muu4**4*sqrtS*vbCLp2 + 
     -         muu4**2*sqrtS*
     -          (-(tred*vbCLp1) + vbCLp2 - muu1**2*vbCLp2 + 
     -            2*tred*vbCLp2) + 
     -         muu3**2*(2*muu4*(vbA + vbBR) + 
     -            sqrtS*(-1 + muu2**2 - 2*tred)*vbCLp1 + 
     -            muu4**2*sqrtS*(vbCLp1 - vbCLp2) + sqrtS*tred*vbCLp2)
     -           + muu3**3*(vbA + vbBL + muu4*sqrtS*vbCRp1) + 
     -         muu3*((-1 + muu2**2 + 2*muu4**2 - tred)*vbA + 
     -            (-1 + muu2**2 + 2*muu4**2 - tred)*vbBL + 
     -            muu4*sqrtS*
     -             ((-1 + 2*muu1**2 + muu2**2 - tred)*vbCRp1 + 
     -               vbCRp2 - 
     -               (muu1**2 + 2*muu2**2 + muu4**2 - tred)*vbCRp2))))
     -    )/dconjg(DenGlU))
! ! 

! ! ! ! c----------------------------------------------Box-T-Kanal Neutralino----------------------------------------------------------
      do n=1,6
      
	  g3L = ALt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3R = ARt1(k,isf1,itt1,ftt1,igen1,fgen1)
	  g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
	  g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))
	
	  g4L = ALt1(k,isf2,itt2,ftt2,igen2,fgen2)
      g4R = ARt1(k,isf2,itt2,ftt2,igen2,fgen2)
	  g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))
	  
	colorf = Nc*(Nc+1d0)/(2*Nc*(Nc+1))
	
	stst2QQ_M2Boxsextet = stst2QQ_M2Boxsextet + 2d0*colorf*dreal((sqrtS**3*(g3Lc*g4Rc*
     -       (muu3**2*muu4*vbA + muu3*muu4**2*vbA + muu3*tred*vbA + 
     -         muu4*tred*vbA + muu3**2*muu4*vbBL + muu4*tred*vbBL + 
     -         muu3*muu4**2*vbBR + muu3*tred*vbBR + 
     -         muu3**3*muu4*sqrtS*vbCLp1 - 
     -         muu3*muu4*sqrtS*tred*vbCLp1 - 
     -         muu3*muu4**3*sqrtS*vbCLp2 + 
     -         muu3*muu4*sqrtS*tred*vbCLp2 - sqrtS*tred*vbCRp1 + 
     -         muu3**2*sqrtS*tred*vbCRp1 - sqrtS*tred**2*vbCRp1 + 
     -         muu1**2*(-(muu4*(vbA + vbBL + muu3*sqrtS*vbCLp1)) + 
     -            muu4**2*sqrtS*vbCRp1 + sqrtS*tred*(vbCRp1 - vbCRp2))
     -           + sqrtS*tred*(1 - muu4**2 + tred)*vbCRp2 - 
     -         muu2**2*(muu3*(vbA + vbBR - muu4*sqrtS*vbCLp2) + 
     -            sqrtS*(muu1**2 - tred)*(vbCRp1 - vbCRp2) + 
     -            muu3**2*sqrtS*vbCRp2)) + 
     -      g3Rc*g4Lc*(muu3**2*muu4*vbA + muu3*muu4**2*vbA + 
     -         muu3*tred*vbA + muu4*tred*vbA + muu3*muu4**2*vbBL + 
     -         muu3*tred*vbBL + muu3**2*muu4*vbBR + muu4*tred*vbBR - 
     -         sqrtS*tred*vbCLp1 + muu3**2*sqrtS*tred*vbCLp1 - 
     -         sqrtS*tred**2*vbCLp1 + sqrtS*tred*vbCLp2 - 
     -         muu4**2*sqrtS*tred*vbCLp2 + sqrtS*tred**2*vbCLp2 + 
     -         muu3**3*muu4*sqrtS*vbCRp1 - 
     -         muu3*muu4*sqrtS*tred*vbCRp1 + 
     -         muu1**2*(muu4**2*sqrtS*vbCLp1 + 
     -            sqrtS*tred*(vbCLp1 - vbCLp2) - 
     -            muu4*(vbA + vbBR + muu3*sqrtS*vbCRp1)) - 
     -         muu3*muu4**3*sqrtS*vbCRp2 + 
     -         muu3*muu4*sqrtS*tred*vbCRp2 - 
     -         muu2**2*(sqrtS*(muu1**2 - tred)*(vbCLp1 - vbCLp2) + 
     -            muu3**2*sqrtS*vbCLp2 + 
     -            muu3*(vbA + vbBL - muu4*sqrtS*vbCRp2))) + 
     -      (g3Rc*g4Rc*((-1 + muu3 + muu4)*(1 + muu3 + muu4)*vbA - 
     -            vbBR + muu3**2*(vbBR + muu4*sqrtS*vbCLp1) + 
     -            muu4*(muu4*vbBR - muu4**2*sqrtS*vbCLp2 + 
     -               sqrtS*((muu1**2 - tred)*(vbCLp1 - vbCLp2) + 
     -                  vbCLp2)) + muu3**3*sqrtS*vbCRp1 + 
     -            muu3*(2*muu4*vbBL + 
     -               sqrtS*(-1 + muu2**2 - tred)*vbCRp1 - 
     -               muu4**2*sqrtS*vbCRp2 + 
     -               sqrtS*(-muu2**2 + tred)*vbCRp2)) + 
     -         g3Lc*g4Lc*((-1 + muu3 + muu4)*(1 + muu3 + muu4)*vbA - 
     -            vbBL + muu3**3*sqrtS*vbCLp1 + 
     -            muu3*(2*muu4*vbBR + 
     -               sqrtS*(-1 + muu2**2 - tred)*vbCLp1 - 
     -               muu4**2*sqrtS*vbCLp2 + 
     -               sqrtS*(-muu2**2 + tred)*vbCLp2) + 
     -            muu3**2*(vbBL + muu4*sqrtS*vbCRp1) + 
     -            muu4*(muu4*vbBL - muu4**2*sqrtS*vbCRp2 + 
     -               sqrtS*((muu1**2 - tred)*(vbCRp1 - vbCRp2) + 
     -                  vbCRp2))))*muupt(n)))/dconjg(DenGauT(n)))

! ! ! c----------------------------------------------Box-U-Kanal Neutralino----------------------------------------------------------
	  g3L = ALt1(k,isf2,itt2,ftt1,igen2,fgen1)
	  g3R = ARt1(k,isf2,itt2,ftt1,igen2,fgen1)
	  g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
	  g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))
	
	  g4L = ALt1(k,isf1,itt1,ftt2,igen1,fgen2)
	  g4R = ARt1(k,isf1,itt1,ftt2,igen1,fgen2)
	  g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))
	  
	  colorf = Nc*(Nc+1d0)/(2*Nc*(Nc+1))
	
	  stst2QQ_M2Boxsextet = stst2QQ_M2Boxsextet + 2d0*colorf*dreal((sqrtS**3*(g3Lc*g4Rc*
     -       (muu4**3*(vbA + vbBL) + 
     -         muu4*(-1 + muu1**2 - tred)*(vbA + vbBL) + 
     -         muu3**3*(vbA + vbBR + muu4*sqrtS*vbCLp1) + 
     -         muu3*((-1 + muu2**2 + 2*muu4**2 - tred)*vbA + 
     -            (-1 + muu2**2 + 2*muu4**2 - tred)*vbBR + 
     -            muu4*sqrtS*
     -             ((-1 + 2*muu1**2 + muu2**2 - tred)*vbCLp1 + 
     -               vbCLp2 - 
     -               (muu1**2 + 2*muu2**2 + muu4**2 - tred)*vbCLp2)) +
     -           muu3**4*sqrtS*vbCRp1 + 
     -         sqrtS*((muu1**2 - tred)*(muu2**2 - tred) + tred)*
     -          (vbCRp1 - vbCRp2) - muu4**4*sqrtS*vbCRp2 + 
     -         muu4**2*sqrtS*
     -          (-(tred*vbCRp1) + vbCRp2 - muu1**2*vbCRp2 + 
     -            2*tred*vbCRp2) + 
     -         muu3**2*(2*muu4*(vbA + vbBL) + 
     -            sqrtS*(-1 + muu2**2 - 2*tred)*vbCRp1 + 
     -            muu4**2*sqrtS*(vbCRp1 - vbCRp2) + sqrtS*tred*vbCRp2)
     -         ) + g3Rc*g4Lc*
     -       (muu4**3*(vbA + vbBR) + 
     -         muu4*(-1 + muu1**2 - tred)*(vbA + vbBR) + 
     -         muu3**4*sqrtS*vbCLp1 + 
     -         sqrtS*((muu1**2 - tred)*(muu2**2 - tred) + tred)*
     -          (vbCLp1 - vbCLp2) - muu4**4*sqrtS*vbCLp2 + 
     -         muu4**2*sqrtS*
     -          (-(tred*vbCLp1) + vbCLp2 - muu1**2*vbCLp2 + 
     -            2*tred*vbCLp2) + 
     -         muu3**2*(2*muu4*(vbA + vbBR) + 
     -            sqrtS*(-1 + muu2**2 - 2*tred)*vbCLp1 + 
     -            muu4**2*sqrtS*(vbCLp1 - vbCLp2) + sqrtS*tred*vbCLp2)
     -           + muu3**3*(vbA + vbBL + muu4*sqrtS*vbCRp1) + 
     -         muu3*((-1 + muu2**2 + 2*muu4**2 - tred)*vbA + 
     -            (-1 + muu2**2 + 2*muu4**2 - tred)*vbBL + 
     -            muu4*sqrtS*
     -             ((-1 + 2*muu1**2 + muu2**2 - tred)*vbCRp1 + 
     -               vbCRp2 - 
     -               (muu1**2 + 2*muu2**2 + muu4**2 - tred)*vbCRp2))) +
     -       (g3Rc*g4Rc*((-1 + muu3 + muu4)*(1 + muu3 + muu4)*vbA - 
     -            vbBR + muu3**2*(vbBR + muu4*sqrtS*vbCLp1) + 
     -            muu4*(muu4*vbBR - muu4**2*sqrtS*vbCLp2 + 
     -               sqrtS*((muu1**2 - tred)*(vbCLp1 - vbCLp2) + 
     -                  vbCLp2)) + muu3**3*sqrtS*vbCRp1 + 
     -            muu3*(2*muu4*vbBL + 
     -               sqrtS*(-1 + muu2**2 - tred)*vbCRp1 - 
     -               muu4**2*sqrtS*vbCRp2 + 
     -               sqrtS*(-muu2**2 + tred)*vbCRp2)) + 
     -         g3Lc*g4Lc*((-1 + muu3 + muu4)*(1 + muu3 + muu4)*vbA - 
     -            vbBL + muu3**3*sqrtS*vbCLp1 + 
     -            muu3*(2*muu4*vbBR + 
     -               sqrtS*(-1 + muu2**2 - tred)*vbCLp1 - 
     -               muu4**2*sqrtS*vbCLp2 + 
     -               sqrtS*(-muu2**2 + tred)*vbCLp2) + 
     -            muu3**2*(vbBL + muu4*sqrtS*vbCRp1) + 
     -            muu4*(muu4*vbBL - muu4**2*sqrtS*vbCRp2 + 
     -               sqrtS*((muu1**2 - tred)*(vbCRp1 - vbCRp2) + 
     -                  vbCRp2))))*muupu(n)))/dconjg(DenGauU(n)))
! ! ! 
! ! ! ! 
      enddo

      
    
	
! 	
	end
