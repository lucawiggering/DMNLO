! =========================================================================
	  double precision function NeuQ2qx_M2BoxH()

	  implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "NeuQ2qx_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"


	  double complex B1vL,B1vR,B1sL,B1sR,B2vL,B2vR,B2sL,B2sR,B3vL,B3vR,B3sL,B3sR,B4vL,B4vR,B4sL,B4sR

	  double complex A1Lc,A1Rc,A2Lc,A2Rc,A3Lc,A3Rc
	  double complex FLc,FRc
	  double complex Gc
	  double complex HLc,HRc

	  double precision mq1c,mq2c,mq3c

	  double complex g0,g0L,g0R,g1,g1L,g1R,g2,g2L,g2R,g3,g3L,g3R
	  double precision Mloop0,Mloop1,Mloop2,Mloop3

	  double precision boxGluonT,boxGluinoT,boxGluonS,boxGluinoS
	  double precision MbsH, MbtH, MbuH

	  integer i,k,l
	  double complex xc0k,xc1k,xc2k

	  integer sumU, sumUc, sumT, sumTc
	  integer sumUend, sumUcend, sumTend, sumTcend

  	  external myCget
  	  external myDget

	  sumU = 0
	  sumUc = 0
	  sumT = 0
	  sumTc = 0

	  if ((xtype.eq.5).or.(xtype.eq.6)) then

	  sumTend = 1
	  sumTcend = 1
	  sumUend = 3
	  sumUcend = 3

	  elseif (xtype.eq.7) then

	  sumTend = 0
	  sumTcend = 0
	  sumUend = 3
	  sumUcend = 3

	  elseif (xtype.eq.8) then

	  sumTend = 1
	  sumTcend = 1
	  sumUend = 1
	  sumUcend = 1

	  endif


	  NeuQ2qx_M2BoxH = 0d0

! ------------------------- Box Gluon s-channel -------------------------
!
! 	  print*,'------------------------------------------------------------- Box Gluon s-channel --------------------------------'

	    B4sL = 0d0
	    B4sR = 0d0
	    B4vL = 0d0
	    B4vR = 0d0

	  boxGluonS = 0d0

	  g0  = -GS 
	  g1L =  AL(chimeigen,sfmeigen,sftype,sfgen)
	  g1R =  AR(chimeigen,sfmeigen,sftype,sfgen)
	  g2L =  FL(xtype,pstype,psgen)
	  g2R =  FR(xtype,pstype,psgen)
	  g3  = -GS

	  Mloop0 = 0d0
	  Mloop1 = m2
	  Mloop2 = Mf(pstype,psgen)
	  Mloop3 = m3

	  call myDget(m2**2,m1**2,m4**2,m3**2,sman,uman,Mloop0**2,Mloop1**2,Mloop2**2,Mloop3**2)


	  call myCget(m2**2,sman,m1**2,Mloop1**2,Mloop0**2,Mloop2**2)
	  xc0k = xc0
	  xc1k = xc1
	  xc2k = xc2

	  call myCget(m1**2,m4**2,uman,Mloop1**2,Mloop2**2,Mloop3**2)




	  B4sL =  - (g0*g3*(g1L*(-(g2R*(m3 - Mloop3)*(xc0 + Mloop0**2*xd0 + sman*xd2 + m2**2*(2*xd0 + 3*xd1 + 2*xd2))) + 
     &          g2R*(-m3**3 + 2*m3**2*Mloop3 + m3*Mloop3**2 + 2*m2**2*(-m3 + Mloop3) - 2*Mloop3*uman)*xd3 + 
     &          g2L*Mloop2*(xc0 + Mloop0**2*xd0 - 2*uman*xd0 - uman*xd1 + m1**2*xd2 - tman*xd2 - uman*xd2 + m3*Mloop3*xd3 - 2*uman*xd3 + 
     &             m2**2*(2*xd0 + 3*(xd1 + xd2) + 2*xd3) + m3**2*(2*xd0 + xd1 + 2*xd2 + 3*xd3))) + 
     &       g1R*m1*(g2R*Mloop2*(-m3 + Mloop3)*xd2 + g2L*
     &           (xc0 + xc2k + Mloop0**2*xd0 - 2*uman*xd0 - uman*xd1 + (m1**2 + Mloop3**2 - tman)*xd2 - 3*uman*xd2 + m3*Mloop3*xd3 - 
     &             2*uman*xd3 + m2**2*(2*xd0 + 3*(xd1 + xd2) + 2*xd3) + m3**2*(2*xd0 + xd1 + 3*(xd2 + xd3))))))/(16.*Pi**2)

	  B4sR =  - (g0*g3*(g1R*(-(g2L*(m3 - Mloop3)*(xc0 + Mloop0**2*xd0 + sman*xd2 + m2**2*(2*xd0 + 3*xd1 + 2*xd2))) + 
     &          g2L*(-m3**3 + 2*m3**2*Mloop3 + m3*Mloop3**2 + 2*m2**2*(-m3 + Mloop3) - 2*Mloop3*uman)*xd3 + 
     &          g2R*Mloop2*(xc0 + Mloop0**2*xd0 - 2*uman*xd0 - uman*xd1 + m1**2*xd2 - tman*xd2 - uman*xd2 + m3*Mloop3*xd3 - 2*uman*xd3 + 
     &             m2**2*(2*xd0 + 3*(xd1 + xd2) + 2*xd3) + m3**2*(2*xd0 + xd1 + 2*xd2 + 3*xd3))) + 
     &       g1L*m1*(g2L*Mloop2*(-m3 + Mloop3)*xd2 + g2R*
     &           (xc0 + xc2k + Mloop0**2*xd0 - 2*uman*xd0 - uman*xd1 + (m1**2 + Mloop3**2 - tman)*xd2 - 3*uman*xd2 + m3*Mloop3*xd3 - 
     &             2*uman*xd3 + m2**2*(2*xd0 + 3*(xd1 + xd2) + 2*xd3) + m3**2*(2*xd0 + xd1 + 3*(xd2 + xd3))))))/(16.*Pi**2)

	  B4vR =  - (g0*g3*(g1R*m1*(2*g2R*Mloop2*xd2 + g2L*(-m3 + Mloop3)*(2*xd0 + xd1 + 2*xd2) - 2*g2L*m3*xd3) + 
     &       g1L*(g2L*Mloop2*((-m3 + Mloop3)*(2*xd0 + xd1 + xd2) - 2*m3*xd3) + 
     &          g2R*(3*xc0 - xc0k - xc1k + 3*Mloop0**2*xd0 - 2*uman*xd0 + Mloop3**2*xd1 - 3*uman*xd1 + 
     &             (m1**2 + Mloop3**2 + 2*sman - tman)*xd2 - 3*uman*xd2 - m3*Mloop3*xd3 - 2*uman*xd3 + 
     &             m2**2*(2*xd0 + 5*xd1 + 3*xd2 + 2*xd3) + m3**2*(2*xd0 + 2*xd1 + 3*xd2 + 5*xd3)))))/(16.*Pi**2)

	  B4vL =  - (g0*g3*(g1L*m1*(2*g2L*Mloop2*xd2 + g2R*(-m3 + Mloop3)*(2*xd0 + xd1 + 2*xd2) - 2*g2R*m3*xd3) + 
     &       g1R*(g2R*Mloop2*((-m3 + Mloop3)*(2*xd0 + xd1 + xd2) - 2*m3*xd3) + 
     &          g2L*(3*xc0 - xc0k - xc1k + 3*Mloop0**2*xd0 - 2*uman*xd0 + Mloop3**2*xd1 - 3*uman*xd1 + 
     &             (m1**2 + Mloop3**2 + 2*sman - tman)*xd2 - 3*uman*xd2 - m3*Mloop3*xd3 - 2*uman*xd3 + 
     &             m2**2*(2*xd0 + 5*xd1 + 3*xd2 + 2*xd3) + m3**2*(2*xd0 + 2*xd1 + 3*xd2 + 5*xd3)))))/(16.*Pi**2)




! --- tree-level - s-channel ---

	    A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
	    A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
	    FLc = dconjg(FL(xtype,pstype,psgen))
	    FRc = dconjg(FR(xtype,pstype,psgen))

	    mq1c = Mf(pstype,psgen)

	    boxGluonS = boxGluonS + CF*NC*2d0*MbsH(B4vL,B4vR,B4sL,B4sR,A1Lc,A1Rc,FLc,FRc,mq1c)
	    



! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    boxGluonS = boxGluonS + CF*NC*2d0*MbtH(B4vL,B4vR,B4sL,B4sR,A2Lc,A2Rc,Gc,mq2c)

	  enddo


! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    boxGluonS = boxGluonS + CF*NC*2d0*MbuH(B4vL,B4vR,B4sL,B4sR,A3Lc,A3Rc,HLc,HRc,mq3c)

	  enddo

	  NeuQ2qx_M2BoxH = NeuQ2qx_M2BoxH + boxGluonS
!  	  print *,"boxGluonS_RESULT",boxGluonS
! 	  print*,'------------------------------------------------------------- End of Box Gluon s-channel --------------------------------'
!-----------------------------------------------------------------------------------------------------

! ------------------------- Box Gluon t-channel -------------------------
!
! 	print*,'------------------------------------------------------------- Box Gluon t-channel --------------------------------'

	    B1sL = 0d0
	    B1sR = 0d0
	    B1vL = 0d0
	    B1vR = 0d0


            boxGluonT = 0d0

	  do sumT = 0, sumTend



	  g0 = -GS
	  g1 = -GS
	  g2L = AL(chimeigen,ptmeigen + sumT,pttype,ptgen)
	  g2R = AR(chimeigen,ptmeigen + sumT,pttype,ptgen)
	  g3 = G(xtype,sfmeigen,ptmeigen + sumT,sftype,sfgen)

	  Mloop0 = 0d0
	  Mloop1 = m3
	  Mloop2 = MSf(ptmeigen + sumT,pttype,ptgen)
	  Mloop3 = m2



 	    call myDget(m3**2,m1**2,m4**2,m2**2,tman,uman,Mloop0**2,Mloop1**2,Mloop2**2,Mloop3**2)


	    call myCget(m1**2,m4**2,uman,Mloop1**2,Mloop2**2,Mloop3**2)


	  B1sL = B1sL - (g0*g1*g3*(2*g2L*(m2**2 + m3**2 - uman)*xd0 + g2L*(xc0 + Mloop0**2*xd0) + g2L*(2*m2**2 + m3*(3*m3 + Mloop1) - 2*uman)*xd1 + 
     &       (g2R*m1*(m3 - Mloop1) + g2L*(-m1**2 + 2*m2**2 + m3*(2*m3 + Mloop1) + tman - 2*uman))*xd2 + g2L*(3*m2**2 + m3**2 - uman)*xd3))/(16.*Pi**2)

	  B1sR = B1sR - (g0*g1*g3*(2*g2R*(m2**2 + m3**2 - uman)*xd0 + g2R*(xc0 + Mloop0**2*xd0) + g2R*(2*m2**2 + m3*(3*m3 + Mloop1) - 2*uman)*xd1 + 
     &       (g2L*m1*(m3 - Mloop1) + g2R*(-m1**2 + 2*m2**2 + m3*(2*m3 + Mloop1) + tman - 2*uman))*xd2 + g2R*(3*m2**2 + m3**2 - uman)*xd3))/(16.*Pi**2)

	  B1vR = B1vR - (g0*g1*g3*(2*g2L*(-m3 + Mloop1)*xd0 - 2*g2L*m3*xd1 - 2*(g2R*m1 + g2L*m3)*xd2 + g2L*(-m3 + Mloop1)*xd3))/(16.*Pi**2)

	  B1vL = B1vL - (g0*g1*g3*(2*g2R*(-m3 + Mloop1)*xd0 - 2*g2R*m3*xd1 - 2*(g2L*m1 + g2R*m3)*xd2 + g2R*(-m3 + Mloop1)*xd3))/(16.*Pi**2)



	  enddo


! --- tree-level - s-channel ---

	    A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
	    A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
	    FLc = dconjg(FL(xtype,pstype,psgen))
	    FRc = dconjg(FR(xtype,pstype,psgen))

	    mq1c = Mf(pstype,psgen)

	    boxGluonT = boxGluonT + CF*NC*2d0*MbsH(B1vL,B1vR,B1sL,B1sR,A1Lc,A1Rc,FLc,FRc,mq1c)
!	    print *,"boxGluonT",boxGluonT



! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    boxGluonT = boxGluonT + CF*NC*2d0*MbtH(B1vL,B1vR,B1sL,B1sR,A2Lc,A2Rc,Gc,mq2c)
!	    print *,"boxGluonT",boxGluonT

	  enddo


! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    boxGluonT = boxGluonT + CF*NC*2d0*MbuH(B1vL,B1vR,B1sL,B1sR,A3Lc,A3Rc,HLc,HRc,mq3c)
!	    print *,"boxGluonT",boxGluonT
	  enddo

!  	  print *,"boxGluonT_RESULT",boxGluonT
	  NeuQ2qx_M2BoxH = NeuQ2qx_M2BoxH + boxGluonT
!	  print *,"NeuQ2qx_M2BoxH",NeuQ2qx_M2BoxH

! 	print*,'------------------------------------------------------------- End of Box Gluon t-channel --------------------------------'
!-----------------------------------------------------------------------------------------------------



! ------------------------- Box Gluino s-channel -------------------------
!
! 	print*,'------------------------------------------------------------- Box Gluino s-channel --------------------------------'


	    B2sL = 0d0
	    B2sR = 0d0
	    B2vL = 0d0
	    B2vR = 0d0

            boxGluinoS = 0d0


	  do i = 1,2
		do k = 1,2



	  g0L =  sqrt2*GS*USf(sfmeigen,2,sftype,sfgen)
	  g0R = -sqrt2*GS*USf(sfmeigen,1,sftype,sfgen)
	  g1L = dconjg(AR(chimeigen,i,pstype,psgen))
	  g1R = dconjg(AL(chimeigen,i,pstype,psgen))
	  g2 = G(xtype,i,k,pstype,psgen)
	  g3L =  sqrt2*GS*USf(k,2,pttype,ptgen)
	  g3R = -sqrt2*GS*USf(k,1,pttype,ptgen)

	  Mloop0 = MGl
	  Mloop1 = Mf(pstype,psgen)
	  Mloop2 = MSf(i,pstype,psgen)
	  Mloop3 = MSf(k,pttype,ptgen)


	  call myDget(m2**2,m1**2,m4**2,m3**2,sman,uman,Mloop0**2,Mloop1**2,Mloop2**2,Mloop3**2)
	  call myCget(m1**2,m4**2,uman,Mloop1**2,Mloop2**2,Mloop3**2)



	  B2sL = B2sL + (g2*(g0L*g1L*g3L*Mloop0*Mloop1*xd0 + g0R*g1L*g3L*(xc0 + Mloop0**2*xd0) + g0R*g1L*g3L*m2**2*xd1 - 
     &       g3L*(m1*(g0R*g1L*m1 + g0L*g1R*Mloop0 + g0R*g1R*Mloop1) - g0R*g1L*sman)*xd2 - 
     &       g1L*g3R*m3*(g0R*Mloop0 + g0L*Mloop1)*xd3))/(16.*Pi**2)

	  B2sR = B2sR + (g2*(g0R*g1R*g3R*Mloop0*Mloop1*xd0 + g0L*g1R*g3R*(xc0 + Mloop0**2*xd0) + g0L*g1R*g3R*m2**2*xd1 - 
     &       g3R*(m1*(g0L*g1R*m1 + g0R*g1L*Mloop0 + g0L*g1L*Mloop1) - g0L*g1R*sman)*xd2 - 
     &       g1R*g3L*m3*(g0L*Mloop0 + g0R*Mloop1)*xd3))/(16.*Pi**2)

	  B2vR = B2vR + (g2*(-(g0R*g1L*g3R*Mloop0*xd0) - g1L*g3R*(g0R*Mloop0 + g0L*Mloop1)*xd1 - 
     &       g3R*(g0L*g1R*m1 + g0R*g1L*Mloop0 + g0L*g1L*Mloop1)*xd2 + g0R*g1L*g3L*m3*xd3))/(16.*Pi**2)

	  B2vL = B2vL + (g2*(-(g0L*g1R*g3L*Mloop0*xd0) - g1R*g3L*(g0L*Mloop0 + g0R*Mloop1)*xd1 - 
     &       g3L*(g0R*g1L*m1 + g0L*g1R*Mloop0 + g0R*g1R*Mloop1)*xd2 + g0L*g1R*g3R*m3*xd3))/(16.*Pi**2)



		enddo
	  enddo


! --- tree-level - s-channel ---

	    A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
	    A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
	    FLc = dconjg(FL(xtype,pstype,psgen))
	    FRc = dconjg(FR(xtype,pstype,psgen))

	    mq1c = Mf(pstype,psgen)

	    boxGluinoS = boxGluinoS + CF*NC*2d0*MbsH(B2vL,B2vR,B2sL,B2sR,A1Lc,A1Rc,FLc,FRc,mq1c)



! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    boxGluinoS = boxGluinoS + CF*NC*2d0*MbtH(B2vL,B2vR,B2sL,B2sR,A2Lc,A2Rc,Gc,mq2c)

	  enddo


! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    boxGluinoS = boxGluinoS + CF*NC*2d0*MbuH(B2vL,B2vR,B2sL,B2sR,A3Lc,A3Rc,HLc,HRc,mq3c)

	  enddo

	  NeuQ2qx_M2BoxH = NeuQ2qx_M2BoxH + boxGluinoS
! 	  print *,"boxGluinoS_RESULT",boxGluinoS

!	  print *,"NeuQ2qx_M2BoxH",NeuQ2qx_M2BoxH


! 	print*,'------------------------------------------------------------- End of Box Gluino s-channel --------------------------------'
!-----------------------------------------------------------------------------------------------------

!------------------------- Box Gluino t-channel -------------------------
!
! 	print*,'------------------------------------------------------------- Box Gluino t-channel --------------------------------'


	    B3sL = 0d0
	    B3sR = 0d0
	    B3vL = 0d0
	    B3vR = 0d0

	  boxGluinoT = 0d0

	  do sumT = 0, sumTend

	  g0L =  sqrt2*GS*USf(sfmeigen,2,sftype,sfgen)
	  g0R = -sqrt2*GS*USf(sfmeigen,1,sftype,sfgen)
	  g1L =  sqrt2*GS*USf(ptmeigen+sumT,2,pttype,ptgen)
	  g1R = -sqrt2*GS*USf(ptmeigen+sumT,1,pttype,ptgen)
	  g2L =  dconjg(AR(chimeigen,ptmeigen+sumT,pttype,ptgen))
	  g2R =  dconjg(AL(chimeigen,ptmeigen+sumT,pttype,ptgen))
	  g3L =  FL(xtype,pstype,psgen) 
	  g3R =  FR(xtype,pstype,psgen) 

	  Mloop0 = MGl
	  Mloop1 = MSf(ptmeigen+sumT,pttype,ptgen)
	  Mloop2 = Mf(ftype,fgen)
	  Mloop3 = Mf(pstype,psgen)


	  call myDget(m3**2,m1**2,m4**2,m2**2,tman,uman,Mloop0**2,Mloop1**2,Mloop2**2,Mloop3**2)



	  call myCget(m3**2,tman,m1**2,Mloop1**2,Mloop0**2,Mloop2**2)
	  xc0k = xc0
	  xc1k = xc1
	  xc2k = xc2



	  call myCget(m1**2,m4**2,uman,Mloop1**2,Mloop2**2,Mloop3**2)



	  B3sL = B3sL + (g0R*(-(g1R*m3*Mloop0*(g2R*g3L*m1*xd1 + g2L*(g3L*Mloop2*(xd1 + xd2) + g3R*Mloop3*(xd0 + xd1 + xd2)))) + 
     &  g1L*(g2R*m1*(g3R*Mloop2*Mloop3*xd2 + g3L*(xc0 + xc2k + Mloop0**2*xd0 + (m1**2 + m2**2 + m3**2 + Mloop3**2 - sman - uman)*xd2 + 
     &  m2**2*xd3)) + g2L*(g3R*Mloop3*(xc0 + Mloop0**2*xd0 + tman*xd2 + m2**2*xd3 - uman*xd3 + m3**2*(xd1 + xd3)) + 
     &  g3L*Mloop2*(xc0 + Mloop0**2*xd0 + m1**2*xd2 - sman*xd2 + m2**2*(xd2 + xd3))))) + 
     &  g0L*(-(g1R*m3*(g2R*g3L*m1*Mloop3*xd1 + g2L*(g3L*Mloop2*Mloop3*(xd1 + xd2) + 
     &  g3R*(xc0 - xc0k - xc1k + Mloop0**2*xd0 - uman*xd1 + m1**2*xd2 - (sman + uman)*xd2 + m3**2*(xd1 + xd2) + 
     &  Mloop3**2*(xd1 + xd2) + m2**2*(xd2 + xd3))))) + 
     &  g1L*Mloop0*(g2R*m1*(g3R*Mloop2*xd2 + g3L*Mloop3*(xd0 + xd2)) + 
     &  g2L*(g3L*Mloop2*Mloop3*xd0 + g3R*(xc0 + Mloop0**2*xd0 - uman*xd0 - uman*xd1 + tman*xd2 - uman*xd2 - uman*xd3 + 
     &  m3**2*(xd0 + 2*xd1 + xd2 + xd3) + m2**2*(xd0 + xd1 + xd2 + 2*xd3))))))/(16.*Pi**2)

	  B3sR = B3sR + (g0L*(-(g1L*m3*Mloop0*(g2L*g3R*m1*xd1 + g2R*(g3R*Mloop2*(xd1 + xd2) + g3L*Mloop3*(xd0 + xd1 + xd2)))) + 
     & g1R*(g2L*m1*(g3L*Mloop2*Mloop3*xd2 + g3R*(xc0 + xc2k + Mloop0**2*xd0 + (m1**2 + m2**2 + m3**2 + Mloop3**2 - sman - uman)*xd2 + 
     & m2**2*xd3)) + g2R*(g3L*Mloop3*(xc0 + Mloop0**2*xd0 + tman*xd2 + m2**2*xd3 - uman*xd3 + m3**2*(xd1 + xd3)) + 
     & g3R*Mloop2*(xc0 + Mloop0**2*xd0 + m1**2*xd2 - sman*xd2 + m2**2*(xd2 + xd3))))) + 
     & g0R*(-(g1L*m3*(g2L*g3R*m1*Mloop3*xd1 + g2R*(g3R*Mloop2*Mloop3*(xd1 + xd2) + 
     & g3L*(xc0 - xc0k - xc1k + Mloop0**2*xd0 - uman*xd1 + m1**2*xd2 - (sman + uman)*xd2 + m3**2*(xd1 + xd2) + 
     & Mloop3**2*(xd1 + xd2) + m2**2*(xd2 + xd3))))) + 
     & g1R*Mloop0*(g2L*m1*(g3L*Mloop2*xd2 + g3R*Mloop3*(xd0 + xd2)) + 
     & g2R*(g3R*Mloop2*Mloop3*xd0 + g3L*(xc0 + Mloop0**2*xd0 - uman*xd0 - uman*xd1 + tman*xd2 - uman*xd2 - uman*xd3 + 
     & m3**2*(xd0 + 2*xd1 + xd2 + xd3) + m2**2*(xd0 + xd1 + xd2 + 2*xd3))))))/(16.*Pi**2)

	  B3vR = B3vR + (-(g0R*(g1L*m3*(-(g3L*(g2R*m1*xd1 + g2L*Mloop2*(xd1 + xd2))) + g2L*g3R*Mloop3*xd3) + 
     & g1R*Mloop0*(g2R*g3L*m1*(xd0 + xd2 + xd3) + g2L*(g3R*Mloop3*xd3 + g3L*Mloop2*(xd0 + xd3))))) + 
     & g0L*(-(g1L*g2L*g3R*m3*Mloop0*(xd0 + xd1 + xd2 + xd3)) + 
     & g1R*(g2R*m1*(g3R*Mloop2*xd2 - g3L*Mloop3*xd3) + 
     & g2L*(-(g3L*Mloop2*Mloop3*xd3) + g3R*(xc0 + Mloop0**2*xd0 + m3**2*xd1 + tman*xd2 + (m2 - Mloop3)*(m2 + Mloop3)*xd3)))))/(16.*Pi**2)

	  B3vL = B3vL + (-(g0L*(g1R*m3*(-(g3R*(g2L*m1*xd1 + g2R*Mloop2*(xd1 + xd2))) + g2R*g3L*Mloop3*xd3) + 
     & g1L*Mloop0*(g2L*g3R*m1*(xd0 + xd2 + xd3) + g2R*(g3L*Mloop3*xd3 + g3R*Mloop2*(xd0 + xd3))))) + 
     & g0R*(-(g1R*g2R*g3L*m3*Mloop0*(xd0 + xd1 + xd2 + xd3)) + 
     & g1L*(g2L*m1*(g3L*Mloop2*xd2 - g3R*Mloop3*xd3) + 
     & g2R*(-(g3R*Mloop2*Mloop3*xd3) + g3L*(xc0 + Mloop0**2*xd0 + m3**2*xd1 + tman*xd2 + (m2 - Mloop3)*(m2 + Mloop3)*xd3)))))/(16.*Pi**2)




	    enddo

! --- tree-level - s-channel ---

	    A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
	    A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
	    FLc = dconjg(FL(xtype,pstype,psgen))
	    FRc = dconjg(FR(xtype,pstype,psgen))

	    mq1c = Mf(pstype,psgen)

	    boxGluinoT = boxGluinoT + CF*NC*2d0*MbsH(B3vL,B3vR,B3sL,B3sR,A1Lc,A1Rc,FLc,FRc,mq1c)
!	    print *,"CF*NC*2d0*MbsH(B3vL,B3vR,B3sL,B3sR,A1Lc,A1Rc,FLc,FRc,mq1c)",CF*NC*2d0*MbsH(B3vL,B3vR,B3sL,B3sR,A1Lc,A1Rc,FLc,FRc,mq1c)


! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    boxGluinoT = boxGluinoT + CF*NC*2d0*MbtH(B3vL,B3vR,B3sL,B3sR,A2Lc,A2Rc,Gc,mq2c)
!	    print *,"CF*NC*2d0*MbtH(B3vL,B3vR,B3sL,B3sR,A2Lc,A2Rc,Gc,mq2c)",CF*NC*2d0*MbtH(B3vL,B3vR,B3sL,B3sR,A2Lc,A2Rc,Gc,mq2c)
	  enddo


! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    boxGluinoT = boxGluinoT + CF*NC*2d0*MbuH(B3vL,B3vR,B3sL,B3sR,A3Lc,A3Rc,HLc,HRc,mq3c)
!	    print *,"boxGluinoT",boxGluinoT
	  enddo

!	    print *,"CF*NC*2d0*MbuH(B3vL,B3vR,B3sL,B3sR,A3Lc,A3Rc,HLc,HRc,mq3c)",CF*NC*2d0*MbuH(B3vL,B3vR,B3sL,B3sR,A3Lc,A3Rc,HLc,HRc,mq3c)
	  NeuQ2qx_M2BoxH = NeuQ2qx_M2BoxH + boxGluinoT
!  	  print *,"boxGluinoT_RESULT",boxGluinoT
! 	  print*,'------------------------------------------------------------- End of Box Gluino t-channel --------------------------------'
! -----------------------------------------------------------------------------------------------------





	   end


