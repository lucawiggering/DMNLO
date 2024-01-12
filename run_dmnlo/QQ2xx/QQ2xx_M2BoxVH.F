***********************************************************************************************************************************

      subroutine QQ2xx_M2GluonBoxVH(argx,argyh,VH_TX,VH_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"
#include "../util/DMNLO_LoopIntegrals.h"

      integer i,j,argx,argyh
      complex*16 g0t,g1t,g2t,g3t,g0u,g1u,g2u,g3u,g0Lt,g0Rt,g1Lt,g1Rt,g2Lt,g2Rt,g3Lt,g3Rt,g0Lu,g0Ru,g1Lu,g1Ru,g2Lu,g2Ru,g3Lu,g3Ru,
     &           xb01t,xb02t,xc01t,xc02t,xc03t,xc04t,xc11t,xc12t,xc13t,xc21t,xc22t,xc23t,xd0t,xd1t,xd2t,xd3t,xb01u,xb02u,xc01u,xc02u,xc03u,xc04u,xc11u,xc12u,xc13u,xc21u,xc22u,xc23u,xd0u,xd1u,xd2u,xd3u
      double precision M0Lpt,M1Lpt,M2Lpt,M3Lpt,M0Lpu,M1Lpu,M2Lpu,M3Lpu,VH_TX(6),VH_UX(6),Box_TQ,Box_TS,Box_TV,Box_TT,Box_TU,Box_UQ,Box_US,Box_UV,Box_UT,Box_UU
        
      external myBget,myCget,myDget
  
! =================================================================================================================================	
!                             ---                     Vector-Higgs-final-states                    ---
! =================================================================================================================================

      Box_TT = 0d0
      Box_TU = 0d0
      Box_TQ = 0d0
      Box_TS = 0d0
      Box_TV = 0d0 

      Box_UT = 0d0
      Box_UU = 0d0   
      Box_UQ = 0d0
      Box_US = 0d0 
      Box_UV = 0d0   

      do i= 1, 2

!=====================================================!Box(t)_Gluon================================================================

        M0Lpt=MV(Looptype0)
        M1Lpt=MSf(Loopmeigen1,Looptype1,Loopgen1)
        M2Lpt=MSf(i,Looptypet,Loopgent)
        M3Lpt=MSf(Loopmeigen2,Looptype2,Loopgen2)

        g0t=VSfSf(Looptype0,Loopmeigen1,sfmeigen1,Looptype1,sftype1,Loopgen1,sfgen1)
        g1t=VSfSft(argx,i,Loopmeigen1,Looptypet,Looptype1,Loopgent,Loopgen1)
        g2t=HSfSfu(argyh,i,Loopmeigen2,Looptypet,Looptype2,Loopgent,Loopgen2)
        g3t=VSfSf(Looptype0,sfmeigen2,Loopmeigen2,sftype2,Looptype2,sfgen2,Loopgen2)

!=====================================================!Box(u)_Gluon================================================================
       
        M0Lpu=MV(Looptype0)
        M1Lpu=MSf(Loopmeigen1,Looptype1,Loopgen1)
        M2Lpu=MSf(i,Looptypeu,Loopgenu)
        M3Lpu=MSf(Loopmeigen2,Looptype2,Loopgen2)

        g0u=VSfSf(Looptype0,Loopmeigen1,sfmeigen1,Looptype1,sftype1,Loopgen1,sfgen1)
        g1u=HSfSfu(argyh,Loopmeigen1,i,Looptype1,Looptypeu,Loopgen1,Loopgenu)
        g2u=VSfSft(argx,Loopmeigen2,i,Looptype2,Looptypeu,Loopgen2,Loopgenu)
        g3u=VSfSf(Looptype0,sfmeigen2,Loopmeigen2,sftype2,Looptype2,sfgen2,Loopgen2)

!=====================================================!call Loopfunctions==========================================================

!t
        call myCget(m1**2,m3**2,tman,M0Lpt**2,M1Lpt**2,M2Lpt**2)
       
        xc01t = xc0    
        xc12t = xc1
        xc21t = xc2

        call myCget(tman,m4**2,m2**2,M0Lpt**2,M2Lpt**2,M3Lpt**2)

        xc02t = xc0
        xc11t = xc1
        xc22t = xc2

        call myCget(m3**2,m4**2,sqrtS**2,M1Lpt**2,M2Lpt**2,M3Lpt**2)

        xc03t = xc0
        xc13t = xc1
        xc23t = xc2

        call myDget(m1**2, m3**2, m4**2, m2**2, tman, sqrtS**2, M0Lpt**2, M1Lpt**2, M2Lpt**2, M3Lpt**2)

        xd0t = xd0
        xd1t = xd1
        xd2t = xd2
        xd3t = xd3
!u

        call myCget(m1**2,m4**2,uman,M0Lpu**2,M1Lpu**2,M2Lpu**2)
       
        xc01u = xc0    
        xc12u = xc1
        xc21u = xc2

        call myCget(uman,m3**2,m2**2,M0Lpu**2,M2Lpu**2,M3Lpu**2)

        xc02u = xc0
        xc11u = xc1
        xc22u = xc2

        call myCget(m4**2,m3**2,sqrtS**2,M1Lpu**2,M2Lpu**2,M3Lpu**2)

        xc03u = xc0
        xc13u = xc1
        xc23u = xc2

        call myDget(m1**2, m4**2, m3**2, m2**2, uman, sqrtS**2, M0Lpu**2, M1Lpu**2, M2Lpu**2, M3Lpu**2)

        xd0u = xd0
        xd1u = xd1
        xd2u = xd2
        xd3u = xd3

!===================================================================================================================================

        do j=1,  2   

!Box(t)*t + Box(u)*t
      
        Coupl_Ac = dconjg(VSfSft(argx,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))       
        Coupl_Bc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
        sfmc = MSf(j,proptypetc,propgentc)
        
!t
	Box_TT = Box_TT + 2*colorfactorTU*dreal(  (Coupl_Ac*Coupl_Bc*g0t*g1t*g2t*g3t*((m1**2 - tman)*(xc01t + xc02t - xc03t + 2*xc11t - 2*xc13t + 2*xc21t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 
     &           2*sqrtS**2*xd0t - 2*M0Lpt**2*xd2t + 2*m1**2*xd2t + 2*M1Lpt**2*xd2t + 2*m2**2*xd2t + 2*M3Lpt**2*xd2t - 4*sqrtS**2*xd2t) + 
     &        2*(-((3*m1**2 - m3**2 + tman)*(xc01t + xc02t + xc11t + xc12t + xc21t + xc23t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 
     &                 M0Lpt**2*xd1t + m1**2*xd1t + M1Lpt**2*xd1t + m2**2*xd1t + M3Lpt**2*xd1t - 2*sqrtS**2*xd1t - M0Lpt**2*xd2t + m1**2*xd2t + M1Lpt**2*xd2t + m2**2*xd2t + M3Lpt**2*xd2t - 
     &                 2*sqrtS**2*xd2t))/2. - ((2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(xc22t - xc23t + (-M0Lpt**2 + m1**2 + M1Lpt**2 + m2**2 + M3Lpt**2 - 2*sqrtS**2)*xd3t))/2.)))/
     &    (16.*pi**2*(-sfmc**2 + tman))  )
!u
        Box_UT = Box_UT + 2*colorfactorTU*dreal(  (Coupl_Ac*Coupl_Bc*g0u*g1u*g2u*g3u*((2*m1**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*uman)*
     &         (xc01u + xc02u - xc03u + 2*xc11u - 2*xc13u + 2*xc21u - M0Lpu**2*xd0u + m1**2*xd0u + M1Lpu**2*xd0u + m2**2*xd0u + M3Lpu**2*xd0u - 2*sqrtS**2*xd0u - 2*M0Lpu**2*xd2u + 
     &           2*m1**2*xd2u + 2*M1Lpu**2*xd2u + 2*m2**2*xd2u + 2*M3Lpu**2*xd2u - 4*sqrtS**2*xd2u) - 
     &        (3*m1**2 - m3**2 + tman)*(xc01u + xc02u + xc03u + 2*xc11u + 2*xc12u + 2*xc21u + 2*xc23u - M0Lpu**2*xd0u + m1**2*xd0u + M1Lpu**2*xd0u + m2**2*xd0u + M3Lpu**2*xd0u - 
     &           2*sqrtS**2*xd0u - 2*M0Lpu**2*xd1u + 2*m1**2*xd1u + 2*M1Lpu**2*xd1u + 2*m2**2*xd1u + 2*M3Lpu**2*xd1u - 4*sqrtS**2*xd1u - 2*M0Lpu**2*xd2u + 2*m1**2*xd2u + 2*M1Lpu**2*xd2u + 
     &           2*m2**2*xd2u + 2*M3Lpu**2*xd2u - 4*sqrtS**2*xd2u) - (2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*
     &         (xc01u + xc02u - xc03u + 2*xc22u - 2*xc23u - M0Lpu**2*xd0u + m1**2*xd0u + M1Lpu**2*xd0u + m2**2*xd0u + M3Lpu**2*xd0u - 2*sqrtS**2*xd0u - 2*M0Lpu**2*xd3u + 2*m1**2*xd3u + 
     &           2*M1Lpu**2*xd3u + 2*m2**2*xd3u + 2*M3Lpu**2*xd3u - 4*sqrtS**2*xd3u)))/(32.*pi**2*(-sfmc**2 + tman))  )

!Box(t)*u + Box(u)*u
 
        Coupl_Ac = dconjg(VSfSft(argx,sfmeigen2,j,sftype2,proptypeu,sfgen2,propgenu))
        Coupl_Bc = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
        sfmc = MSf(j,proptypeuc,propgenuc)
        
!t
	Box_TU = Box_TU + 2*colorfactorTU*dreal(  (Coupl_Ac*Coupl_Bc*g0t*g1t*g2t*g3t*((-m2**2 + uman)*(xc01t + xc02t - xc03t + 2*xc11t - 2*xc13t + 2*xc21t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 
     &           2*sqrtS**2*xd0t - 2*M0Lpt**2*xd2t + 2*m1**2*xd2t + 2*M1Lpt**2*xd2t + 2*m2**2*xd2t + 2*M3Lpt**2*xd2t - 4*sqrtS**2*xd2t) + 
     &        2*(-((3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(xc01t + xc02t + xc11t + xc12t + xc21t + xc23t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 
     &                 2*sqrtS**2*xd0t - M0Lpt**2*xd1t + m1**2*xd1t + M1Lpt**2*xd1t + m2**2*xd1t + M3Lpt**2*xd1t - 2*sqrtS**2*xd1t - M0Lpt**2*xd2t + m1**2*xd2t + M1Lpt**2*xd2t + m2**2*xd2t + 
     &                 M3Lpt**2*xd2t - 2*sqrtS**2*xd2t))/2. - ((3*m2**2 - m3**2 + uman)*(xc22t - xc23t + (-M0Lpt**2 + m1**2 + M1Lpt**2 + m2**2 + M3Lpt**2 - 2*sqrtS**2)*xd3t))/2.)))/
     &    (16.*pi**2*(-sfmc**2 + uman))  )
!u
	Box_UU = Box_UU + 2*colorfactorTU*dreal(   (Coupl_Ac*Coupl_Bc*g0u*g1u*g2u*g3u*(-((2*m2**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*tman)*
     &           (xc01u + xc02u - xc03u + 2*xc11u - 2*xc13u + 2*xc21u - M0Lpu**2*xd0u + m1**2*xd0u + M1Lpu**2*xd0u + m2**2*xd0u + M3Lpu**2*xd0u - 2*sqrtS**2*xd0u - 2*M0Lpu**2*xd2u + 
     &             2*m1**2*xd2u + 2*M1Lpu**2*xd2u + 2*m2**2*xd2u + 2*M3Lpu**2*xd2u - 4*sqrtS**2*xd2u)) - 
     &        (3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(xc01u + xc02u + xc03u + 2*xc11u + 2*xc12u + 2*xc21u + 2*xc23u - M0Lpu**2*xd0u + m1**2*xd0u + M1Lpu**2*xd0u + m2**2*xd0u + 
     &           M3Lpu**2*xd0u - 2*sqrtS**2*xd0u - 2*M0Lpu**2*xd1u + 2*m1**2*xd1u + 2*M1Lpu**2*xd1u + 2*m2**2*xd1u + 2*M3Lpu**2*xd1u - 4*sqrtS**2*xd1u - 2*M0Lpu**2*xd2u + 2*m1**2*xd2u + 
     &           2*M1Lpu**2*xd2u + 2*m2**2*xd2u + 2*M3Lpu**2*xd2u - 4*sqrtS**2*xd2u) - 
     &        (3*m2**2 - m3**2 + uman)*(xc01u + xc02u - xc03u + 2*xc22u - 2*xc23u - M0Lpu**2*xd0u + m1**2*xd0u + M1Lpu**2*xd0u + m2**2*xd0u + M3Lpu**2*xd0u - 2*sqrtS**2*xd0u - 
     &           2*M0Lpu**2*xd3u + 2*m1**2*xd3u + 2*M1Lpu**2*xd3u + 2*m2**2*xd3u + 2*M3Lpu**2*xd3u - 4*sqrtS**2*xd3u)))/(32.*pi**2*(-sfmc**2 + uman))  )


        enddo

!==================================================================================================================================     
!Box(t)*s(Higgs)+Box(u)*s(Higgs)       

        do j=sumstartH, sumendH

        Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(JLR(argx,j,argyh))
        propac = dconjg(DenSmH(j))

!t        
	Box_TS = Box_TS + 2*colorfactorSH*dreal(  (Coupl_Ac*Coupl_Bc*g0t*g1t*g2t*g3t*(((m1**2 + m2**2 + m3**2 - m4**2 + sqrtS**2 - tman - uman)*
     &           (xc01t + xc02t - xc03t + 2*xc11t - 2*xc13t + 2*xc21t - (M0Lpt**2 - m1**2 - M1Lpt**2 - m2**2 - M3Lpt**2 + 2*sqrtS**2)*(xd0t + 2*xd2t)))/2. + 
     &        2*(-((2*m1**2 - m2**2 + m4**2 + sqrtS**2 - uman)*(xc01t + xc02t + xc11t + xc12t + xc21t + xc23t - 
     &                 (M0Lpt**2 - m1**2 - M1Lpt**2 - m2**2 - M3Lpt**2 + 2*sqrtS**2)*(xd0t + xd1t + xd2t)))/2. - 
     &           ((m1**2 - 2*m2**2 - m4**2 - sqrtS**2 + tman)*(xc22t - xc23t + (-M0Lpt**2 + m1**2 + M1Lpt**2 + m2**2 + M3Lpt**2 - 2*sqrtS**2)*xd3t))/2.)))/(16.*pi**2*propac)  )
!u
	Box_US = Box_US + 2*colorfactorSH*dreal(  (Coupl_Ac*Coupl_Bc*g0u*g1u*g2u*g3u*((m1**2 + m2**2 + 4*m4**2 - tman - uman)*
     &         (xc01u + xc02u - xc03u + 2*xc11u - 2*xc13u + 2*xc21u - (M0Lpu**2 - m1**2 - M1Lpu**2 - m2**2 - M3Lpu**2 + 2*sqrtS**2)*(xd0u + 2*xd2u)) - 
     &        (2*m1**2 - m2**2 + m4**2 + sqrtS**2 - uman)*(xc01u + xc02u + xc03u + 2*xc11u + 2*xc12u + 2*xc21u + 2*xc23u - 
     &           (M0Lpu**2 - m1**2 - M1Lpu**2 - m2**2 - M3Lpu**2 + 2*sqrtS**2)*(xd0u + 2*(xd1u + xd2u))) + 
     &        (-m1**2 + 2*m2**2 + m4**2 + sqrtS**2 - tman)*(xc01u + xc02u - xc03u + 2*xc22u - 2*xc23u - (M0Lpu**2 - m1**2 - M1Lpu**2 - m2**2 - M3Lpu**2 + 2*sqrtS**2)*(xd0u + 2*xd3u))))/
     &    (32.*pi**2*propac)  )
  

        enddo

!==================================================================================================================================
!Box(t)*Vector + Box(u)*Vector

        do j=sumstartV, sumendV

        Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = dconjg(LLR(j,argx,argyh))
        propac = dconjg(DenSmV(j))

!t        
        Box_TV = Box_TV + 2*colorfactorSV*dreal(  -(Coupl_Ac*Coupl_Bc*g0t*g1t*g2t*g3t*(((m1**2 + m3**2 - tman)*(xc01t + xc02t - xc03t + 2*xc11t - 2*xc13t + 2*xc21t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + 
     &              M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 2*M0Lpt**2*xd2t + 2*m1**2*xd2t + 2*M1Lpt**2*xd2t + 2*m2**2*xd2t + 2*M3Lpt**2*xd2t - 4*sqrtS**2*xd2t))/2. - 
     &         ((m2**2 + m3**2 - uman)*(xc01t + xc02t - xc03t + 2*xc11t - 2*xc13t + 2*xc21t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 
     &              2*M0Lpt**2*xd2t + 2*m1**2*xd2t + 2*M1Lpt**2*xd2t + 2*m2**2*xd2t + 2*M3Lpt**2*xd2t - 4*sqrtS**2*xd2t))/2. - 
     &         2*m1**2*(xc01t + xc02t + xc11t + xc12t + xc21t + xc23t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - M0Lpt**2*xd1t + m1**2*xd1t + 
     &            M1Lpt**2*xd1t + m2**2*xd1t + M3Lpt**2*xd1t - 2*sqrtS**2*xd1t - M0Lpt**2*xd2t + m1**2*xd2t + M1Lpt**2*xd2t + m2**2*xd2t + M3Lpt**2*xd2t - 2*sqrtS**2*xd2t) + 
     &         (-m1**2 - m2**2 + sqrtS**2)*(xc22t - xc23t + (-M0Lpt**2 + m1**2 + M1Lpt**2 + m2**2 + M3Lpt**2 - 2*sqrtS**2)*xd3t) - 
     &         2*(((m1**2 + m2**2 - sqrtS**2)*(xc01t + xc02t + xc11t + xc12t + xc21t + xc23t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 
     &                 M0Lpt**2*xd1t + m1**2*xd1t + M1Lpt**2*xd1t + m2**2*xd1t + M3Lpt**2*xd1t - 2*sqrtS**2*xd1t - M0Lpt**2*xd2t + m1**2*xd2t + M1Lpt**2*xd2t + m2**2*xd2t + M3Lpt**2*xd2t - 
     &                 2*sqrtS**2*xd2t))/2. + m2**2*(xc22t - xc23t + (-M0Lpt**2 + m1**2 + M1Lpt**2 + m2**2 + M3Lpt**2 - 2*sqrtS**2)*xd3t))))/(16.*pi**2*propac)  )
!u        
        Box_UV = Box_UV + 2*colorfactorSV*dreal(  -(Coupl_Ac*Coupl_Bc*g0u*g1u*g2u*g3u*(-((m2**2 + m4**2 - tman)*(xc01u + xc02u - xc03u + 2*xc11u - 2*xc13u + 2*xc21u - M0Lpu**2*xd0u + m1**2*xd0u + M1Lpu**2*xd0u + m2**2*xd0u + 
     &              M3Lpu**2*xd0u - 2*sqrtS**2*xd0u - 2*M0Lpu**2*xd2u + 2*m1**2*xd2u + 2*M1Lpu**2*xd2u + 2*m2**2*xd2u + 2*M3Lpu**2*xd2u - 4*sqrtS**2*xd2u)) + 
     &         (m1**2 + m4**2 - uman)*(xc01u + xc02u - xc03u + 2*xc11u - 2*xc13u + 2*xc21u - M0Lpu**2*xd0u + m1**2*xd0u + M1Lpu**2*xd0u + m2**2*xd0u + M3Lpu**2*xd0u - 2*sqrtS**2*xd0u - 
     &            2*M0Lpu**2*xd2u + 2*m1**2*xd2u + 2*M1Lpu**2*xd2u + 2*m2**2*xd2u + 2*M3Lpu**2*xd2u - 4*sqrtS**2*xd2u) - 
     &         2*m1**2*(xc01u + xc02u + xc03u + 2*xc11u + 2*xc12u + 2*xc21u + 2*xc23u - M0Lpu**2*xd0u + m1**2*xd0u + M1Lpu**2*xd0u + m2**2*xd0u + M3Lpu**2*xd0u - 2*sqrtS**2*xd0u - 
     &            2*M0Lpu**2*xd1u + 2*m1**2*xd1u + 2*M1Lpu**2*xd1u + 2*m2**2*xd1u + 2*M3Lpu**2*xd1u - 4*sqrtS**2*xd1u - 2*M0Lpu**2*xd2u + 2*m1**2*xd2u + 2*M1Lpu**2*xd2u + 2*m2**2*xd2u + 
     &            2*M3Lpu**2*xd2u - 4*sqrtS**2*xd2u) - (m1**2 + m2**2 - sqrtS**2)*
     &          (xc01u + xc02u + xc03u + 2*xc11u + 2*xc12u + 2*xc21u + 2*xc23u - M0Lpu**2*xd0u + m1**2*xd0u + M1Lpu**2*xd0u + m2**2*xd0u + M3Lpu**2*xd0u - 2*sqrtS**2*xd0u - 2*M0Lpu**2*xd1u + 
     &            2*m1**2*xd1u + 2*M1Lpu**2*xd1u + 2*m2**2*xd1u + 2*M3Lpu**2*xd1u - 4*sqrtS**2*xd1u - 2*M0Lpu**2*xd2u + 2*m1**2*xd2u + 2*M1Lpu**2*xd2u + 2*m2**2*xd2u + 2*M3Lpu**2*xd2u - 
     &            4*sqrtS**2*xd2u) - 2*m2**2*(xc01u + xc02u - xc03u + 2*xc22u - 2*xc23u - M0Lpu**2*xd0u + m1**2*xd0u + M1Lpu**2*xd0u + m2**2*xd0u + M3Lpu**2*xd0u - 2*sqrtS**2*xd0u - 
     &            2*M0Lpu**2*xd3u + 2*m1**2*xd3u + 2*M1Lpu**2*xd3u + 2*m2**2*xd3u + 2*M3Lpu**2*xd3u - 4*sqrtS**2*xd3u) - 
     &         (m1**2 + m2**2 - sqrtS**2)*(xc01u + xc02u - xc03u + 2*xc22u - 2*xc23u - M0Lpu**2*xd0u + m1**2*xd0u + M1Lpu**2*xd0u + m2**2*xd0u + M3Lpu**2*xd0u - 2*sqrtS**2*xd0u - 
     &            2*M0Lpu**2*xd3u + 2*m1**2*xd3u + 2*M1Lpu**2*xd3u + 2*m2**2*xd3u + 2*M3Lpu**2*xd3u - 4*sqrtS**2*xd3u)))/(32.*pi**2*propac)  )

        enddo           

!================================================================================================================================== 

      enddo

      VH_TX(1) = Box_TQ
      VH_TX(2) = Box_TS
      VH_TX(3) = Box_TV
      VH_TX(4) = Box_TT
      VH_TX(5) = Box_TU

      VH_TX(6) = VH_TX(1) + VH_TX(2) + VH_TX(3) + VH_TX(4) + VH_TX(5)

      VH_UX(1) = Box_UQ
      VH_UX(2) = Box_US
      VH_UX(3) = Box_UV
      VH_UX(4) = Box_UT
      VH_UX(5) = Box_UU

      VH_UX(6) = VH_UX(1) + VH_UX(2) + VH_UX(3) + VH_UX(4) + VH_UX(5)

!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine QQ2xx_M2GluonTriangleBoxVH(argx,argyh,VH_TX,VH_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer j,argx,argyh
      complex*16 CA1,CA3,CB1,CB3
      double precision VH_TX(6),VH_UX(6),Box_TQ,Box_TS,Box_TV,Box_TT,Box_TU,Box_UQ,Box_US,Box_UV,Box_UT,Box_UU
  
! =================================================================================================================================	
!                             ---                     Vector-Higgs-final-states                    ---
! =================================================================================================================================

      Box_TT = 0d0
      Box_TU = 0d0
      Box_TQ = 0d0
      Box_TS = 0d0
      Box_TV = 0d0 

      Box_UT = 0d0
      Box_UU = 0d0   
      Box_UQ = 0d0
      Box_US = 0d0 
      Box_UV = 0d0   
       
!=====================================================!Box-Triangle_Gluon(Typ2t)=====================================================
       
      call QQ2xx_BT2VH(0,m2**2,m4**2,tman,Looptype0,argyh,argx,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,Loopmeigen2,Looptype2,Loopgen2,Looptypet,Loopgent,CA1,CA3) 

      CA1 = CA1/(16d0*pi**2)
      CA3 = CA3/(16d0*pi**2)
!=====================================================!Box-Triangle_Gluon(Typ1u)====================================================

      call QQ2xx_BT1HV(1,uman,m4**2,m1**2,Looptype0,argyh,argx,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,Looptypeu,Loopgenu,Loopmeigen1,Looptype1,Loopgen1,CB1,CB3) 

      CB1 = CB1/(16d0*pi**2)
      CB3 = CB3/(16d0*pi**2)
!====================================================================================================================================

      do j=1,  2   

!BoxTri(t)*t + BoxTri(u)*t
      
        Coupl_Ac = dconjg(VSfSft(argx,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))       
        Coupl_Bc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
        sfmc = MSf(j,proptypetc,propgentc)
        
!t
	Box_TT = Box_TT + 2*colorfactorTU*dreal(  (Coupl_Ac*Coupl_Bc*(CA3*(m1**2 - m3**2 + 3*tman) + CA1*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)))/(2.*(sfmc**2 - tman))  )
!u
        Box_UT = Box_UT + 2*colorfactorTU*dreal(  (Coupl_Ac*Coupl_Bc*(CB3*(3*m1**2 - m3**2 + tman) + CB1*(2*m1**2 + m2**2 - m3**2 - 2*m4**2 + uman)))/(2.*(sfmc**2 - tman))  )

!BoxTri(t)*u + BoxTri(u)*u
 
        Coupl_Ac = dconjg(VSfSft(argx,sfmeigen2,j,sftype2,proptypeu,sfgen2,propgenu))
        Coupl_Bc = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
        sfmc = MSf(j,proptypeuc,propgenuc)
        
!t
	Box_TU = Box_TU + 2*colorfactorTU*dreal(  (Coupl_Ac*Coupl_Bc*(CA3*(m1**2 + 2*m2**2 - m3**2 - 2*m4**2 + tman) + CA1*(3*m2**2 - m3**2 + uman)))/(2.*(sfmc**2 - uman))  )
!u
	Box_UU = Box_UU + 2*colorfactorTU*dreal(  (Coupl_Ac*Coupl_Bc*(CB3*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman) + CB1*(m2**2 - m3**2 + 3*uman)))/(2.*(sfmc**2 - uman))  )


      enddo

!==================================================================================================================================     
       
      do j=sumstartH, sumendH

!BoxTri(t)*s(Higgs)+BoxTri(u)*s(Higgs)

        Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(JLR(argx,j,argyh))
        propac = dconjg(DenSmH(j))

!t        
	Box_TS = Box_TS + 2*colorfactorSH*dreal(  -(Coupl_Ac*Coupl_Bc*(CA3*(m1**2 - 2*m2**2 - m3**2 + 2*m4**2 + tman) + CA1*(m1**2 - 2*m2**2 - m4**2 - sqrtS**2 + tman)))/(2.*propac)  )
!u
	Box_US = Box_US + 2*colorfactorSH*dreal(   (Coupl_Ac*Coupl_Bc*(CB1*(-2*m1**2 + m2**2 - m3**2 + 2*m4**2 + uman) + CB3*(-2*m1**2 + m2**2 - m4**2 - sqrtS**2 + uman)))/(2.*propac)  )
  

      enddo

!==================================================================================================================================

      do j=sumstartV, sumendV

!BoxTri(t)*Vector + BoxTri(u)*Vector

        Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = dconjg(LLR(j,argx,argyh))
        propac = dconjg(DenSmV(j))

!t        
        Box_TV = Box_TV + 2*colorfactorSV*dreal( (Coupl_Ac*Coupl_Bc*(CA1*(m1**2 + 3*m2**2 - sqrtS**2) + CA3*(2*m1**2 + 2*m2**2 - sqrtS**2 + tman - uman)))/(2.*propac)  )
!u        
        Box_UV = Box_UV + 2*colorfactorSV*dreal( (Coupl_Ac*Coupl_Bc*(CB3*(3*m1**2 + m2**2 - sqrtS**2) + CB1*(2*m1**2 + 2*m2**2 - sqrtS**2 - tman + uman)))/(2.*propac)  )

      enddo           
    
!==================================================================================================================================

      VH_TX(1) = Box_TQ
      VH_TX(2) = Box_TS
      VH_TX(3) = Box_TV
      VH_TX(4) = Box_TT
      VH_TX(5) = Box_TU

      VH_TX(6) = VH_TX(1) + VH_TX(2) + VH_TX(3) + VH_TX(4) + VH_TX(5)

      VH_UX(1) = Box_UQ
      VH_UX(2) = Box_US
      VH_UX(3) = Box_UV
      VH_UX(4) = Box_UT
      VH_UX(5) = Box_UU

      VH_UX(6) = VH_UX(1) + VH_UX(2) + VH_UX(3) + VH_UX(4) + VH_UX(5)

!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine QQ2xx_M2sQuarkTriangleBoxVH(argx,argyh,VH_TX,VH_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer j,argx,argyh
      complex*16 CA1,CA3,CB1,CB3
      double precision VH_TX(6),VH_UX(6),Box_TQ,Box_TS,Box_TV,Box_TT,Box_TU,Box_UQ,Box_US,Box_UV,Box_UT,Box_UU

! =================================================================================================================================
!                             ---                     Vector-Higgs-final-states                    ---
! =================================================================================================================================

      Box_TT = 0d0
      Box_TU = 0d0
      Box_TQ = 0d0
      Box_TS = 0d0
      Box_TV = 0d0 

      Box_UT = 0d0
      Box_UU = 0d0   
      Box_UQ = 0d0
      Box_US = 0d0 
      Box_UV = 0d0  

!=====================================================!Box-Triangle_Gluon(Typ2t)====================================================

      call  QQ2xx_BTsf4VH2(0,m4**2,m3**2,sqrtS**2,argx,argyh,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,CA1,CA3) 
       
      CA1 = CA1/(16d0*pi**2)
      CA3 = CA3/(16d0*pi**2)
!=====================================================!Box-Triangle_Gluon(Typ1u)====================================================
       
      call  QQ2xx_BTsf4VH1(1,sqrtS**2,m3**2,m4**2,argx,argyh,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,CB1,CB3) 

      CB1 = CB1/(16d0*pi**2)
      CB3 = CB3/(16d0*pi**2)
!===================================================================================================================================

      do j=1,  2   

!BoxTri(t)*t + BoxTri(u)*t
      
        Coupl_Ac = dconjg(VSfSft(argx,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))       
        Coupl_Bc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
        sfmc = MSf(j,proptypetc,propgentc)
        
!t
	Box_TT = Box_TT + 2*dreal(  (Coupl_Ac*Coupl_Bc*(CA1*(2*m1**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*uman) + CA3*(m1**2 - 3*m2**2 - 2*m3**2 + 2*sqrtS**2 + tman + uman)))/(2.*(sfmc**2 - tman))  )
!u
        Box_UT = Box_UT + 2*dreal( -(Coupl_Ac*Coupl_Bc*(CB3*(2*m1**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*uman) + CB1*(m1**2 - 3*m2**2 - 2*m3**2 + 2*sqrtS**2 + tman + uman)))/(2.*(sfmc**2 - tman))  )

!BoxTri(t)*u + BoxTri(u)*u
 
        Coupl_Ac = dconjg(VSfSft(argx,sfmeigen2,j,sftype2,proptypeu,sfgen2,propgenu))
        Coupl_Bc = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
        sfmc = MSf(j,proptypeuc,propgenuc)
        
!t
	Box_TU = Box_TU + 2*dreal(   (Coupl_Ac*Coupl_Bc*(CA1*(m1**2 - m2**2 - 2*m4**2 + tman - uman) + CA3*(m1**2 - 3*m2**2 - 2*m4**2 + tman + uman)))/(2.*(sfmc**2 - uman))  )
!u
	Box_UU = Box_UU + 2*dreal(  -(Coupl_Ac*Coupl_Bc*(CB3*(m1**2 - m2**2 - 2*m4**2 + tman - uman) + CB1*(m1**2 - 3*m2**2 - 2*m4**2 + tman + uman)))/(2.*(sfmc**2 - uman))  )


      enddo

!==================================================================================================================================     
!BoxTri(t)*s(Higgs)+BoxTri(u)*s(Higgs)
       
      do j=sumstartH, sumendH

        Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(JLR(argx,j,argyh))
        propac = dconjg(DenSmH(j))

!t        
	Box_TS = Box_TS + 2*dreal(  (Coupl_Ac*Coupl_Bc*(CA1*m3**2 + CA3*m3**2 - 3*CA1*m4**2 - CA3*m4**2 - (CA1 + 3*CA3)*sqrtS**2))/(2.*propac)  )
!u
	Box_US = Box_US + 2*dreal(  (Coupl_Ac*Coupl_Bc*(CB3*(-m3**2 + 3*m4**2 + sqrtS**2) + CB1*(-m3**2 + m4**2 + 3*sqrtS**2)))/(2.*propac)  )
  

      enddo

!==================================================================================================================================
!BoxTri(t)*Vector + BoxTri(u)*Vector

      do j=sumstartV, sumendV

        Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = dconjg(LLR(j,argx,argyh))
        propac = dconjg(DenSmV(j))

!t        
        Box_TV = Box_TV + 2*dreal(   (Coupl_Ac*Coupl_Bc*(2*CA3*(m1**2 - m2**2) + CA1*(m1**2 - m2**2 + tman - uman)))/(2.*propac)  )
!u        
        Box_UV = Box_UV + 2*dreal(  -(Coupl_Ac*Coupl_Bc*(2*CB1*(m1**2 - m2**2) + CB3*(m1**2 - m2**2 + tman - uman)))/(2.*propac)  )

      enddo           
    
!==================================================================================================================================

      VH_TX(1) = Box_TQ
      VH_TX(2) = Box_TS
      VH_TX(3) = Box_TV
      VH_TX(4) = Box_TT
      VH_TX(5) = Box_TU

      VH_TX(6) = VH_TX(1) + VH_TX(2) + VH_TX(3) + VH_TX(4) + VH_TX(5)

      VH_UX(1) = Box_UQ
      VH_UX(2) = Box_US
      VH_UX(3) = Box_UV
      VH_UX(4) = Box_UT
      VH_UX(5) = Box_UU

      VH_UX(6) = VH_UX(1) + VH_UX(2) + VH_UX(3) + VH_UX(4) + VH_UX(5)

!===================================================================================================================================

      end




!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine QQ2xx_M2GluinoBoxVH(argx,argyh,VH_TX,VH_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"
#include "../util/DMNLO_LoopIntegrals.h"

      integer j,argx,argyh
      complex*16 g0t,g1t,g2t,g3t,g0u,g1u,g2u,g3u,g0Lt,g0Rt,g1Lt,g1Rt,g2Lt,g2Rt,g3Lt,g3Rt,g0Lu,g0Ru,g1Lu,g1Ru,g2Lu,g2Ru,g3Lu,g3Ru,
     &           xb01t,xb02t,xc01t,xc02t,xc03t,xc04t,xc11t,xc12t,xc13t,xc21t,xc22t,xc23t,xd0t,xd1t,xd2t,xd3t,xb01u,xb02u,xc01u,xc02u,xc03u,xc04u,xc11u,xc12u,xc13u,xc21u,xc22u,xc23u,xd0u,xd1u,xd2u,xd3u
      double precision M0Lpt,M1Lpt,M2Lpt,M3Lpt,M0Lpu,M1Lpu,M2Lpu,M3Lpu,VH_TX(6),VH_UX(6),Box_TQ,Box_TS,Box_TV,Box_TT,Box_TU,Box_UQ,Box_US,Box_UV,Box_UT,Box_UU
        
      external myBget,myCget,myDget
  
! =================================================================================================================================	
!                             ---                     Vector-Higgs-final-states                    ---
! =================================================================================================================================

      Box_TT = 0d0
      Box_TU = 0d0
      Box_TQ = 0d0
      Box_TS = 0d0
      Box_TV = 0d0 

      Box_UT = 0d0
      Box_UU = 0d0   
      Box_UQ = 0d0
      Box_US = 0d0 
      Box_UV = 0d0    

!=====================================================!Box(t)_Gluino===============================================================


      M0Lpt=MGl
      M1Lpt=Mf(Looptype1,Loopgen1)
      M2Lpt=Mf(Looptypet,Loopgent)
      M3Lpt=Mf(Looptype2,Loopgen2)

      g0Lt=GluR(sfmeigen1,sftype1,Looptype1,sfgen1,Loopgen1)
      g0Rt=GluL(sfmeigen1,sftype1,Looptype1,sfgen1,Loopgen1)
      g1Lt=BLt(argx,Looptype1,Looptypet,Loopgen1,Loopgent)
      g1Rt=BRt(argx,Looptype1,Looptypet,Loopgen1,Loopgent)
      g2Lt=DLu(argyh,Looptypet,Looptype2,Loopgent,Loopgen2)
      g2Rt=DRu(argyh,Looptypet,Looptype2,Loopgent,Loopgen2)
      g3Lt=GluL(sfmeigen2,sftype2,Looptype2,sfgen2,Loopgen2)
      g3Rt=GluR(sfmeigen2,sftype2,Looptype2,sfgen2,Loopgen2)

!=====================================================!Box(u)_Gluino================================================================


      M0Lpu=MGl
      M1Lpu=Mf(Looptype1,Loopgen1)
      M2Lpu=Mf(Looptypeu,Loopgenu)
      M3Lpu=Mf(Looptype2,Loopgen2)

      g0Lu=GluR(sfmeigen1,sftype1,Looptype1,sfgen1,Loopgen1)
      g0Ru=GluL(sfmeigen1,sftype1,Looptype1,sfgen1,Loopgen1)
      g1Lu=DLu(argyh,Looptype1,Looptypeu,Loopgen1,Loopgenu)
      g1Ru=DRu(argyh,Looptype1,Looptypeu,Loopgen1,Loopgenu)
      g2Lu=BLt(argx,Looptypeu,Looptype2,Loopgenu,Loopgen2)
      g2Ru=BRt(argx,Looptypeu,Looptype2,Loopgenu,Loopgen2)
      g3Lu=GluL(sfmeigen2,sftype2,Looptype2,sfgen2,Loopgen2)
      g3Ru=GluR(sfmeigen2,sftype2,Looptype2,sfgen2,Loopgen2)

!=====================================================!call Loopfunctions==========================================================

!t
      call myBget(m4**2,M2Lpt**2,M3Lpt**2)

      xb01t = xb0

      call myBget(tman,M0Lpt**2,M2Lpt**2)

      xb02t = xb0

      call myCget(m1**2,sqrtS**2,m2**2,M0Lpt**2,M1Lpt**2,M3Lpt**2)
       
      xc01t = xc0    

      call myCget(m1**2,m3**2,tman,M0Lpt**2,M1Lpt**2,M2Lpt**2)
       
      xc02t = xc0    
      xc12t = xc1
      xc22t = xc2

      call myCget(tman,m4**2,m2**2,M0Lpt**2,M2Lpt**2,M3Lpt**2)

      xc03t = xc0

      call myCget(m3**2,m4**2,sqrtS**2,M1Lpt**2,M2Lpt**2,M3Lpt**2)

      xc04t = xc0
      xc11t = xc1
      xc21t = xc2

      call myDget(m1**2, m3**2, m4**2, m2**2, tman, sqrtS**2, M0Lpt**2, M1Lpt**2, M2Lpt**2, M3Lpt**2)

      xd0t = xd0
      xd1t = xd1
      xd2t = xd2
      xd3t = xd3

!u
 
      call myBget(m3**2,M2Lpu**2,M3Lpu**2)

      xb01u = xb0

      call myBget(uman,M0Lpu**2,M2Lpu**2)

      xb02u = xb0

      call myCget(m1**2,sqrtS**2,m2**2,M0Lpu**2,M1Lpu**2,M3Lpu**2)
       
      xc01u = xc0    

      call myCget(m1**2,m4**2,uman,M0Lpu**2,M1Lpu**2,M2Lpu**2)
       
      xc02u = xc0    

      call myCget(uman,m3**2,m2**2,M0Lpu**2,M2Lpu**2,M3Lpu**2)

      xc03u = xc0
      xc12u = xc1
      xc22u = xc2

      call myCget(m4**2,m3**2,sqrtS**2,M1Lpu**2,M2Lpu**2,M3Lpu**2)

      xc04u = xc0
      xc11u = xc1
      xc21u = xc2

      call myDget(m1**2, m4**2, m3**2, m2**2, uman, sqrtS**2, M0Lpu**2, M1Lpu**2, M2Lpu**2, M3Lpu**2)

      xd0u = xd0
      xd1u = xd1
      xd2u = xd2
      xd3u = xd3

!===================================================================================================================================


      do j=1,  2   

!Box(t)*t + Box(u)*t
      
        Coupl_Ac = dconjg(VSfSft(argx,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))       
        Coupl_Bc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
        sfmc = MSf(j,proptypetc,propgentc)
        
!t
	Box_TT = Box_TT + 2*colorfactorTU*dreal( (Coupl_Ac*Coupl_Bc*(2*(m1**2 - tman)*(g0Rt*(g1Lt*(g2Rt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*xd2t + g3Lt*(2*xc11t + (2*M0Lpt**2 - m1**2 - m2**2 + sqrtS**2)*xd2t)) + 
     &                 g2Lt*(g3Rt*M0Lpt*(xc02t + xc03t + 2*xc22t + M1Lpt**2*xd0t + M3Lpt**2*xd0t - sqrtS**2*xd0t + m3**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t - sqrtS**2*xd2t) + 
     &                    g3Lt*M3Lpt*(xc03t + xc04t + 2*xc11t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t))) + 
     &              g1Rt*M1Lpt*(2*g2Lt*M2Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*xd2t + g2Rt*
     &                  (2*g3Rt*M0Lpt*M3Lpt*(xd0t + xd2t) + g3Lt*(xc02t + xc04t + 2*xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t - m2**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t + tman*xd2t)))
     &              ) + g0Lt*(g1Rt*(g2Lt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*xd2t + g3Rt*(2*xc11t + (2*M0Lpt**2 - m1**2 - m2**2 + sqrtS**2)*xd2t)) + 
     &                 g2Rt*(g3Lt*M0Lpt*(xc02t + xc03t + 2*xc22t + M1Lpt**2*xd0t + M3Lpt**2*xd0t - sqrtS**2*xd0t + m3**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t - sqrtS**2*xd2t) + 
     &                    g3Rt*M3Lpt*(xc03t + xc04t + 2*xc11t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t))) + 
     &              g1Lt*M1Lpt*(2*g2Rt*M2Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*xd2t + g2Lt*
     &                  (2*g3Lt*M0Lpt*M3Lpt*(xd0t + xd2t) + g3Rt*(xc02t + xc04t + 2*xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t - m2**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t + tman*xd2t)))
     &              )) + (3*m1**2 - m3**2 + tman)*(-(g0Rt*(g1Lt*(g2Rt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*(xd0t + xd1t + xd2t) + 
     &                      g3Lt*(xc02t - xc04t - 2*xc21t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 2*M0Lpt**2*xd1t - m1**2*xd1t - m2**2*xd1t + sqrtS**2*xd1t + 2*M0Lpt**2*xd2t - 
     &                         m1**2*xd2t - m2**2*xd2t + sqrtS**2*xd2t)) + g2Lt*
     &                    (g3Rt*M0Lpt*(xc01t + 2*xc02t + xc03t + 2*xc12t + 2*xc22t + M1Lpt**2*xd0t + M2Lpt**2*xd0t + 2*M3Lpt**2*xd0t - m4**2*xd0t - sqrtS**2*xd0t + m3**2*xd1t + 
     &                         2*M3Lpt**2*xd1t - m4**2*xd1t - sqrtS**2*xd1t + m3**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t - sqrtS**2*xd2t) + 
     &                      g3Lt*M3Lpt*(xc01t + xc03t - 2*xc21t + 2*M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t + 2*M0Lpt**2*xd1t - m1**2*xd1t + m3**2*xd1t - 
     &                         tman*xd1t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t))) + 
     &                g1Rt*M1Lpt*(2*g2Lt*M2Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*(xd1t + xd2t) + 
     &                   g2Rt*(2*g3Rt*M0Lpt*M3Lpt*(xd0t + xd1t + xd2t) + g3Lt*(xc02t + xc04t + 2*xc12t + 2*xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t - m2**2*xd1t + 2*M3Lpt**2*xd1t - 
     &                         m4**2*xd1t + tman*xd1t - m2**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t + tman*xd2t))))) - 
     &           g0Lt*(g1Rt*(g2Lt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*(xd0t + xd1t + xd2t) + 
     &                    g3Rt*(xc02t - xc04t - 2*xc21t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 2*M0Lpt**2*xd1t - m1**2*xd1t - m2**2*xd1t + sqrtS**2*xd1t + 2*M0Lpt**2*xd2t - 
     &                       m1**2*xd2t - m2**2*xd2t + sqrtS**2*xd2t)) + g2Rt*(g3Lt*M0Lpt*
     &                     (xc01t + 2*xc02t + xc03t + 2*xc12t + 2*xc22t + M1Lpt**2*xd0t + M2Lpt**2*xd0t + 2*M3Lpt**2*xd0t - m4**2*xd0t - sqrtS**2*xd0t + m3**2*xd1t + 2*M3Lpt**2*xd1t - 
     &                       m4**2*xd1t - sqrtS**2*xd1t + m3**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t - sqrtS**2*xd2t) + 
     &                    g3Rt*M3Lpt*(xc01t + xc03t - 2*xc21t + 2*M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t + 2*M0Lpt**2*xd1t - m1**2*xd1t + m3**2*xd1t - 
     &                       tman*xd1t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t))) + 
     &              g1Lt*M1Lpt*(2*g2Rt*M2Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*(xd1t + xd2t) + 
     &                 g2Lt*(2*g3Lt*M0Lpt*M3Lpt*(xd0t + xd1t + xd2t) + g3Rt*(xc02t + xc04t + 2*xc12t + 2*xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t - m2**2*xd1t + 2*M3Lpt**2*xd1t - 
     &                       m4**2*xd1t + tman*xd1t - m2**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t + tman*xd2t))))) + 
     &        (2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(g0Rt*(g1Rt*M1Lpt*
     &               (-2*g2Lt*M2Lpt*(g3Lt*M3Lpt*xd3t + g3Rt*M0Lpt*(xd0t + xd3t)) + 
     &                 g2Rt*(-2*g3Rt*M0Lpt*M3Lpt*xd3t + g3Lt*(xc01t + xc04t + M0Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t + m2**2*xd3t - 2*M3Lpt**2*xd3t + m4**2*xd3t - tman*xd3t))) + 
     &              g1Lt*(-(g2Rt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*xd3t + g3Lt*(xc03t + xc04t + 2*xc21t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 2*M0Lpt**2*xd3t - m1**2*xd3t - m2**2*xd3t + 
     &                         sqrtS**2*xd3t))) + g2Lt*(g3Rt*M0Lpt*(xc01t + xc03t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - m3**2*xd0t - m3**2*xd3t - 2*M3Lpt**2*xd3t + m4**2*xd3t + sqrtS**2*xd3t) + 
     &                    g3Lt*M3Lpt*(-2*xc21t + (-2*M0Lpt**2 + m1**2 - m3**2 + tman)*xd3t)))) + 
     &           g0Lt*(g1Lt*M1Lpt*(-2*g2Rt*M2Lpt*(g3Rt*M3Lpt*xd3t + g3Lt*M0Lpt*(xd0t + xd3t)) + 
     &                 g2Lt*(-2*g3Lt*M0Lpt*M3Lpt*xd3t + g3Rt*(xc01t + xc04t + M0Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t + m2**2*xd3t - 2*M3Lpt**2*xd3t + m4**2*xd3t - tman*xd3t))) + 
     &              g1Rt*(-(g2Lt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*xd3t + g3Rt*(xc03t + xc04t + 2*xc21t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 2*M0Lpt**2*xd3t - m1**2*xd3t - m2**2*xd3t + 
     &                         sqrtS**2*xd3t))) + g2Rt*(g3Lt*M0Lpt*(xc01t + xc03t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - m3**2*xd0t - m3**2*xd3t - 2*M3Lpt**2*xd3t + m4**2*xd3t + sqrtS**2*xd3t) + 
     &                    g3Rt*M3Lpt*(-2*xc21t + (-2*M0Lpt**2 + m1**2 - m3**2 + tman)*xd3t)))))))/(32.*pi**2*(-sfmc**2 + tman))  )
!u
        Box_UT = Box_UT + 2*colorfactorTU*dreal(  (Coupl_Ac*Coupl_Bc*((3*m1**2 - m3**2 + tman)*(g0Ru*(-(g1Ru*M1Lpu*(2*g2Ru*M2Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*(xd1u + xd2u) + 
     &                   g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(xd0u + xd1u + xd2u) + g3Lu*(xc02u - xc04u - 2*xc21u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd1u - m2**2*xd1u + 
     &                         m3**2*xd1u - uman*xd1u + 2*M0Lpu**2*xd2u - m2**2*xd2u + m3**2*xd2u - uman*xd2u)))) + 
     &              g1Lu*(-(g2Lu*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*(xd0u + xd1u + xd2u) + 
     &                      g3Lu*(xc02u - xc04u - 2*xc21u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd1u - m1**2*xd1u - m2**2*xd1u + sqrtS**2*xd1u + 2*M0Lpu**2*xd2u - 
     &                         m1**2*xd2u - m2**2*xd2u + sqrtS**2*xd2u))) + g2Ru*
     &                  (g3Ru*M0Lpu*(xc01u - xc03u - 2*xc12u - M1Lpu**2*xd0u + M2Lpu**2*xd0u - m3**2*xd0u + sqrtS**2*xd0u - 2*M1Lpu**2*xd1u - m3**2*xd1u + m4**2*xd1u + sqrtS**2*xd1u - 
     &                       2*M1Lpu**2*xd2u - m3**2*xd2u + m4**2*xd2u + sqrtS**2*xd2u) + 
     &                    g3Lu*M3Lpu*(xc01u - xc03u - 2*xc12u + m1**2*xd0u - M1Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u + m1**2*xd1u - 2*M1Lpu**2*xd1u + m4**2*xd1u - uman*xd1u + m1**2*xd2u - 
     &                       2*M1Lpu**2*xd2u + m4**2*xd2u - uman*xd2u)))) + g0Lu*
     &            (-(g1Lu*M1Lpu*(2*g2Lu*M2Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*(xd1u + xd2u) + 
     &                   g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(xd0u + xd1u + xd2u) + g3Ru*(xc02u - xc04u - 2*xc21u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd1u - m2**2*xd1u + 
     &                         m3**2*xd1u - uman*xd1u + 2*M0Lpu**2*xd2u - m2**2*xd2u + m3**2*xd2u - uman*xd2u)))) + 
     &              g1Ru*(-(g2Ru*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*(xd0u + xd1u + xd2u) + 
     &                      g3Ru*(xc02u - xc04u - 2*xc21u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd1u - m1**2*xd1u - m2**2*xd1u + sqrtS**2*xd1u + 2*M0Lpu**2*xd2u - 
     &                         m1**2*xd2u - m2**2*xd2u + sqrtS**2*xd2u))) + g2Lu*
     &                  (g3Lu*M0Lpu*(xc01u - xc03u - 2*xc12u - M1Lpu**2*xd0u + M2Lpu**2*xd0u - m3**2*xd0u + sqrtS**2*xd0u - 2*M1Lpu**2*xd1u - m3**2*xd1u + m4**2*xd1u + sqrtS**2*xd1u - 
     &                       2*M1Lpu**2*xd2u - m3**2*xd2u + m4**2*xd2u + sqrtS**2*xd2u) + 
     &                    g3Ru*M3Lpu*(xc01u - xc03u - 2*xc12u + m1**2*xd0u - M1Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u + m1**2*xd1u - 2*M1Lpu**2*xd1u + m4**2*xd1u - uman*xd1u + m1**2*xd2u - 
     &                       2*M1Lpu**2*xd2u + m4**2*xd2u - uman*xd2u))))) + (2*m1**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*uman)*
     &         (g0Ru*(g1Ru*M1Lpu*(2*g2Ru*M2Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*xd2u + 
     &                 g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(xd0u + xd2u) + g3Lu*(xc02u + xc04u + 2*xc11u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd2u - m2**2*xd2u + m3**2*xd2u - 
     &                       uman*xd2u))) + g1Lu*(g2Lu*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*xd2u + g3Lu*(2*xc11u + (2*M0Lpu**2 - m1**2 - m2**2 + sqrtS**2)*xd2u)) + 
     &                 g2Ru*(g3Ru*M0Lpu*(xc02u + xc03u + 2*xc12u + M1Lpu**2*xd0u + M3Lpu**2*xd0u - sqrtS**2*xd0u + 2*M1Lpu**2*xd2u + m3**2*xd2u - m4**2*xd2u - sqrtS**2*xd2u) + 
     &                    g3Lu*M3Lpu*(xc03u + xc04u + 2*xc12u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u - m1**2*xd2u + 2*M1Lpu**2*xd2u - m4**2*xd2u + uman*xd2u)))) + 
     &           g0Lu*(g1Lu*M1Lpu*(2*g2Lu*M2Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*xd2u + 
     &                 g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(xd0u + xd2u) + g3Ru*(xc02u + xc04u + 2*xc11u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd2u - m2**2*xd2u + m3**2*xd2u - 
     &                       uman*xd2u))) + g1Ru*(g2Ru*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*xd2u + g3Ru*(2*xc11u + (2*M0Lpu**2 - m1**2 - m2**2 + sqrtS**2)*xd2u)) + 
     &                 g2Lu*(g3Lu*M0Lpu*(xc02u + xc03u + 2*xc12u + M1Lpu**2*xd0u + M3Lpu**2*xd0u - sqrtS**2*xd0u + 2*M1Lpu**2*xd2u + m3**2*xd2u - m4**2*xd2u - sqrtS**2*xd2u) + 
     &                    g3Ru*M3Lpu*(xc03u + xc04u + 2*xc12u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u - m1**2*xd2u + 2*M1Lpu**2*xd2u - m4**2*xd2u + uman*xd2u))))) + 
     &        (2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(-(g0Ru*(g1Lu*(g2Lu*M2Lpu*
     &                    (2*g3Ru*M0Lpu*M3Lpu*xd3u + g3Lu*(xc03u + xc04u + 2*xc21u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + 2*M0Lpu**2*xd3u - m1**2*xd3u - m2**2*xd3u + sqrtS**2*xd3u)) + 
     &                   g2Ru*(g3Ru*M0Lpu*(xc01u + xc03u + 2*xc22u + M1Lpu**2*xd0u + M2Lpu**2*xd0u - m4**2*xd0u + 2*M1Lpu**2*xd3u + m3**2*xd3u - m4**2*xd3u - sqrtS**2*xd3u) + 
     &                      g3Lu*M3Lpu*(2*xc22u - (m1**2 - 2*M1Lpu**2 + m4**2 - uman)*xd3u))) + 
     &                g1Ru*M1Lpu*(2*g2Ru*M2Lpu*(g3Lu*M3Lpu*xd3u + g3Ru*M0Lpu*(xd0u + xd3u)) + 
     &                   g2Lu*(2*g3Ru*M0Lpu*M3Lpu*xd3u + g3Lu*(xc01u + xc04u + 2*xc21u + M0Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u + 2*M0Lpu**2*xd3u - m2**2*xd3u + m3**2*xd3u - uman*xd3u))))) 
     &            -   g0Lu*(g1Ru*(g2Ru*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*xd3u + g3Ru*(xc03u + xc04u + 2*xc21u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + 2*M0Lpu**2*xd3u - m1**2*xd3u - m2**2*xd3u + 
     &                       sqrtS**2*xd3u)) + g2Lu*(g3Lu*M0Lpu*(xc01u + xc03u + 2*xc22u + M1Lpu**2*xd0u + M2Lpu**2*xd0u - m4**2*xd0u + 2*M1Lpu**2*xd3u + m3**2*xd3u - m4**2*xd3u - 
     &                       sqrtS**2*xd3u) + g3Ru*M3Lpu*(2*xc22u - (m1**2 - 2*M1Lpu**2 + m4**2 - uman)*xd3u))) + 
     &              g1Lu*M1Lpu*(2*g2Lu*M2Lpu*(g3Ru*M3Lpu*xd3u + g3Lu*M0Lpu*(xd0u + xd3u)) + 
     &                 g2Ru*(2*g3Lu*M0Lpu*M3Lpu*xd3u + g3Ru*(xc01u + xc04u + 2*xc21u + M0Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u + 2*M0Lpu**2*xd3u - m2**2*xd3u + m3**2*xd3u - uman*xd3u)))))))/
     &    (32.*pi**2*(-sfmc**2 + tman))   )

!Box(t)*u + Box(u)*u
 
        Coupl_Ac = dconjg(VSfSft(argx,sfmeigen2,j,sftype2,proptypeu,sfgen2,propgenu))
        Coupl_Bc = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
        sfmc = MSf(j,proptypeuc,propgenuc)
        
!t
	Box_TU = Box_TU + 2*colorfactorTU*dreal( (Coupl_Ac*Coupl_Bc*((-m2**2 + uman)*(g0Rt*(g1Lt*(g2Rt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*xd2t + g3Lt*(2*xc11t + (2*M0Lpt**2 - m1**2 - m2**2 + sqrtS**2)*xd2t)) + 
     &                 g2Lt*(g3Rt*M0Lpt*(xc02t + xc03t + 2*xc22t + M1Lpt**2*xd0t + M3Lpt**2*xd0t - sqrtS**2*xd0t + m3**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t - sqrtS**2*xd2t) + 
     &                    g3Lt*M3Lpt*(xc03t + xc04t + 2*xc11t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t))) + 
     &              g1Rt*M1Lpt*(2*g2Lt*M2Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*xd2t + g2Rt*
     &                  (2*g3Rt*M0Lpt*M3Lpt*(xd0t + xd2t) + g3Lt*(xc02t + xc04t + 2*xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t - m2**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t + tman*xd2t)))
     &              ) + g0Lt*(g1Rt*(g2Lt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*xd2t + g3Rt*(2*xc11t + (2*M0Lpt**2 - m1**2 - m2**2 + sqrtS**2)*xd2t)) + 
     &                 g2Rt*(g3Lt*M0Lpt*(xc02t + xc03t + 2*xc22t + M1Lpt**2*xd0t + M3Lpt**2*xd0t - sqrtS**2*xd0t + m3**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t - sqrtS**2*xd2t) + 
     &                    g3Rt*M3Lpt*(xc03t + xc04t + 2*xc11t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t))) + 
     &              g1Lt*M1Lpt*(2*g2Rt*M2Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*xd2t + g2Lt*
     &                  (2*g3Lt*M0Lpt*M3Lpt*(xd0t + xd2t) + g3Rt*(xc02t + xc04t + 2*xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t - m2**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t + tman*xd2t)))
     &              )) - ((3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(g0Rt*
     &              (g1Lt*(g2Rt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*(xd0t + xd1t + xd2t) + 
     &                      g3Lt*(xc02t - xc04t - 2*xc21t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 2*M0Lpt**2*xd1t - m1**2*xd1t - m2**2*xd1t + sqrtS**2*xd1t + 2*M0Lpt**2*xd2t - 
     &                         m1**2*xd2t - m2**2*xd2t + sqrtS**2*xd2t)) + g2Lt*
     &                    (g3Rt*M0Lpt*(xc01t + 2*xc02t + xc03t + 2*xc12t + 2*xc22t + M1Lpt**2*xd0t + M2Lpt**2*xd0t + 2*M3Lpt**2*xd0t - m4**2*xd0t - sqrtS**2*xd0t + m3**2*xd1t + 
     &                         2*M3Lpt**2*xd1t - m4**2*xd1t - sqrtS**2*xd1t + m3**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t - sqrtS**2*xd2t) + 
     &                      g3Lt*M3Lpt*(xc01t + xc03t - 2*xc21t + 2*M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t + 2*M0Lpt**2*xd1t - m1**2*xd1t + m3**2*xd1t - 
     &                         tman*xd1t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t))) + 
     &                g1Rt*M1Lpt*(2*g2Lt*M2Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*(xd1t + xd2t) + 
     &                   g2Rt*(2*g3Rt*M0Lpt*M3Lpt*(xd0t + xd1t + xd2t) + g3Lt*(xc02t + xc04t + 2*xc12t + 2*xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t - m2**2*xd1t + 2*M3Lpt**2*xd1t - 
     &                         m4**2*xd1t + tman*xd1t - m2**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t + tman*xd2t)))) + 
     &             g0Lt*(g1Rt*(g2Lt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*(xd0t + xd1t + xd2t) + 
     &                      g3Rt*(xc02t - xc04t - 2*xc21t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 2*M0Lpt**2*xd1t - m1**2*xd1t - m2**2*xd1t + sqrtS**2*xd1t + 2*M0Lpt**2*xd2t - 
     &                         m1**2*xd2t - m2**2*xd2t + sqrtS**2*xd2t)) + g2Rt*
     &                    (g3Lt*M0Lpt*(xc01t + 2*xc02t + xc03t + 2*xc12t + 2*xc22t + M1Lpt**2*xd0t + M2Lpt**2*xd0t + 2*M3Lpt**2*xd0t - m4**2*xd0t - sqrtS**2*xd0t + m3**2*xd1t + 
     &                         2*M3Lpt**2*xd1t - m4**2*xd1t - sqrtS**2*xd1t + m3**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t - sqrtS**2*xd2t) + 
     &                      g3Rt*M3Lpt*(xc01t + xc03t - 2*xc21t + 2*M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t + 2*M0Lpt**2*xd1t - m1**2*xd1t + m3**2*xd1t - 
     &                         tman*xd1t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t))) + 
     &                g1Lt*M1Lpt*(2*g2Rt*M2Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*(xd1t + xd2t) + 
     &                   g2Lt*(2*g3Lt*M0Lpt*M3Lpt*(xd0t + xd1t + xd2t) + g3Rt*(xc02t + xc04t + 2*xc12t + 2*xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t - m2**2*xd1t + 2*M3Lpt**2*xd1t - 
     &                         m4**2*xd1t + tman*xd1t - m2**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t + tman*xd2t))))))/2. + 
     &        ((3*m2**2 - m3**2 + uman)*(g0Rt*(g1Rt*M1Lpt*(-2*g2Lt*M2Lpt*(g3Lt*M3Lpt*xd3t + g3Rt*M0Lpt*(xd0t + xd3t)) + 
     &                   g2Rt*(-2*g3Rt*M0Lpt*M3Lpt*xd3t + g3Lt*(xc01t + xc04t + M0Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t + m2**2*xd3t - 2*M3Lpt**2*xd3t + m4**2*xd3t - tman*xd3t))) + 
     &                g1Lt*(-(g2Rt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*xd3t + g3Lt*(xc03t + xc04t + 2*xc21t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 2*M0Lpt**2*xd3t - m1**2*xd3t - m2**2*xd3t + 
     &                           sqrtS**2*xd3t))) + g2Lt*(g3Rt*M0Lpt*(xc01t + xc03t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - m3**2*xd0t - m3**2*xd3t - 2*M3Lpt**2*xd3t + m4**2*xd3t + sqrtS**2*xd3t) + 
     &                      g3Lt*M3Lpt*(-2*xc21t + (-2*M0Lpt**2 + m1**2 - m3**2 + tman)*xd3t)))) + 
     &             g0Lt*(g1Lt*M1Lpt*(-2*g2Rt*M2Lpt*(g3Rt*M3Lpt*xd3t + g3Lt*M0Lpt*(xd0t + xd3t)) + 
     &                   g2Lt*(-2*g3Lt*M0Lpt*M3Lpt*xd3t + g3Rt*(xc01t + xc04t + M0Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t + m2**2*xd3t - 2*M3Lpt**2*xd3t + m4**2*xd3t - tman*xd3t))) + 
     &                g1Rt*(-(g2Lt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*xd3t + g3Rt*(xc03t + xc04t + 2*xc21t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 2*M0Lpt**2*xd3t - m1**2*xd3t - m2**2*xd3t + 
     &                           sqrtS**2*xd3t))) + g2Rt*(g3Lt*M0Lpt*(xc01t + xc03t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - m3**2*xd0t - m3**2*xd3t - 2*M3Lpt**2*xd3t + m4**2*xd3t + sqrtS**2*xd3t) + 
     &                      g3Rt*M3Lpt*(-2*xc21t + (-2*M0Lpt**2 + m1**2 - m3**2 + tman)*xd3t))))))/2.))/(16.*pi**2*(-sfmc**2 + uman))  )
!u
	Box_UU = Box_UU + 2*colorfactorTU*dreal(  (Coupl_Ac*Coupl_Bc*((3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(g0Ru*
     &            (-(g1Ru*M1Lpu*(2*g2Ru*M2Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*(xd1u + xd2u) + 
     &                   g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(xd0u + xd1u + xd2u) + g3Lu*(xc02u - xc04u - 2*xc21u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd1u - m2**2*xd1u + 
     &                         m3**2*xd1u - uman*xd1u + 2*M0Lpu**2*xd2u - m2**2*xd2u + m3**2*xd2u - uman*xd2u)))) + 
     &              g1Lu*(-(g2Lu*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*(xd0u + xd1u + xd2u) + 
     &                      g3Lu*(xc02u - xc04u - 2*xc21u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd1u - m1**2*xd1u - m2**2*xd1u + sqrtS**2*xd1u + 2*M0Lpu**2*xd2u - 
     &                         m1**2*xd2u - m2**2*xd2u + sqrtS**2*xd2u))) + g2Ru*
     &                  (g3Ru*M0Lpu*(xc01u - xc03u - 2*xc12u - M1Lpu**2*xd0u + M2Lpu**2*xd0u - m3**2*xd0u + sqrtS**2*xd0u - 2*M1Lpu**2*xd1u - m3**2*xd1u + m4**2*xd1u + sqrtS**2*xd1u - 
     &                       2*M1Lpu**2*xd2u - m3**2*xd2u + m4**2*xd2u + sqrtS**2*xd2u) + 
     &                    g3Lu*M3Lpu*(xc01u - xc03u - 2*xc12u + m1**2*xd0u - M1Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u + m1**2*xd1u - 2*M1Lpu**2*xd1u + m4**2*xd1u - uman*xd1u + m1**2*xd2u - 
     &                       2*M1Lpu**2*xd2u + m4**2*xd2u - uman*xd2u)))) + g0Lu*
     &            (-(g1Lu*M1Lpu*(2*g2Lu*M2Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*(xd1u + xd2u) + 
     &                   g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(xd0u + xd1u + xd2u) + g3Ru*(xc02u - xc04u - 2*xc21u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd1u - m2**2*xd1u + 
     &                         m3**2*xd1u - uman*xd1u + 2*M0Lpu**2*xd2u - m2**2*xd2u + m3**2*xd2u - uman*xd2u)))) + 
     &              g1Ru*(-(g2Ru*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*(xd0u + xd1u + xd2u) + 
     &                      g3Ru*(xc02u - xc04u - 2*xc21u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd1u - m1**2*xd1u - m2**2*xd1u + sqrtS**2*xd1u + 2*M0Lpu**2*xd2u - 
     &                         m1**2*xd2u - m2**2*xd2u + sqrtS**2*xd2u))) + g2Lu*
     &                  (g3Lu*M0Lpu*(xc01u - xc03u - 2*xc12u - M1Lpu**2*xd0u + M2Lpu**2*xd0u - m3**2*xd0u + sqrtS**2*xd0u - 2*M1Lpu**2*xd1u - m3**2*xd1u + m4**2*xd1u + sqrtS**2*xd1u - 
     &                       2*M1Lpu**2*xd2u - m3**2*xd2u + m4**2*xd2u + sqrtS**2*xd2u) + 
     &                    g3Ru*M3Lpu*(xc01u - xc03u - 2*xc12u + m1**2*xd0u - M1Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u + m1**2*xd1u - 2*M1Lpu**2*xd1u + m4**2*xd1u - uman*xd1u + m1**2*xd2u - 
     &                       2*M1Lpu**2*xd2u + m4**2*xd2u - uman*xd2u))))) + (2*m2**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*tman)*
     &         (-(g0Ru*(g1Ru*M1Lpu*(2*g2Ru*M2Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*xd2u + 
     &                   g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(xd0u + xd2u) + g3Lu*(xc02u + xc04u + 2*xc11u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd2u - m2**2*xd2u + m3**2*xd2u - 
     &                         uman*xd2u))) + g1Lu*(g2Lu*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*xd2u + g3Lu*(2*xc11u + (2*M0Lpu**2 - m1**2 - m2**2 + sqrtS**2)*xd2u)) + 
     &                   g2Ru*(g3Ru*M0Lpu*(xc02u + xc03u + 2*xc12u + M1Lpu**2*xd0u + M3Lpu**2*xd0u - sqrtS**2*xd0u + 2*M1Lpu**2*xd2u + m3**2*xd2u - m4**2*xd2u - sqrtS**2*xd2u) + 
     &                      g3Lu*M3Lpu*(xc03u + xc04u + 2*xc12u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u - m1**2*xd2u + 2*M1Lpu**2*xd2u - m4**2*xd2u + uman*xd2u))))) - 
     &           g0Lu*(g1Lu*M1Lpu*(2*g2Lu*M2Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*xd2u + 
     &                 g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(xd0u + xd2u) + g3Ru*(xc02u + xc04u + 2*xc11u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd2u - m2**2*xd2u + m3**2*xd2u - 
     &                       uman*xd2u))) + g1Ru*(g2Ru*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*xd2u + g3Ru*(2*xc11u + (2*M0Lpu**2 - m1**2 - m2**2 + sqrtS**2)*xd2u)) + 
     &                 g2Lu*(g3Lu*M0Lpu*(xc02u + xc03u + 2*xc12u + M1Lpu**2*xd0u + M3Lpu**2*xd0u - sqrtS**2*xd0u + 2*M1Lpu**2*xd2u + m3**2*xd2u - m4**2*xd2u - sqrtS**2*xd2u) + 
     &                    g3Ru*M3Lpu*(xc03u + xc04u + 2*xc12u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u - m1**2*xd2u + 2*M1Lpu**2*xd2u - m4**2*xd2u + uman*xd2u))))) + 
     &        (3*m2**2 - m3**2 + uman)*(-(g0Ru*(g1Lu*(g2Lu*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*xd3u + 
     &                      g3Lu*(xc03u + xc04u + 2*xc21u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + 2*M0Lpu**2*xd3u - m1**2*xd3u - m2**2*xd3u + sqrtS**2*xd3u)) + 
     &                   g2Ru*(g3Ru*M0Lpu*(xc01u + xc03u + 2*xc22u + M1Lpu**2*xd0u + M2Lpu**2*xd0u - m4**2*xd0u + 2*M1Lpu**2*xd3u + m3**2*xd3u - m4**2*xd3u - sqrtS**2*xd3u) + 
     &                      g3Lu*M3Lpu*(2*xc22u - (m1**2 - 2*M1Lpu**2 + m4**2 - uman)*xd3u))) + 
     &                g1Ru*M1Lpu*(2*g2Ru*M2Lpu*(g3Lu*M3Lpu*xd3u + g3Ru*M0Lpu*(xd0u + xd3u)) + 
     &                   g2Lu*(2*g3Ru*M0Lpu*M3Lpu*xd3u + g3Lu*(xc01u + xc04u + 2*xc21u + M0Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u + 2*M0Lpu**2*xd3u - m2**2*xd3u + m3**2*xd3u - uman*xd3u))))) 
     &            -   g0Lu*(g1Ru*(g2Ru*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*xd3u + g3Ru*(xc03u + xc04u + 2*xc21u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + 2*M0Lpu**2*xd3u - m1**2*xd3u - m2**2*xd3u + 
     &                       sqrtS**2*xd3u)) + g2Lu*(g3Lu*M0Lpu*(xc01u + xc03u + 2*xc22u + M1Lpu**2*xd0u + M2Lpu**2*xd0u - m4**2*xd0u + 2*M1Lpu**2*xd3u + m3**2*xd3u - m4**2*xd3u - 
     &                       sqrtS**2*xd3u) + g3Ru*M3Lpu*(2*xc22u - (m1**2 - 2*M1Lpu**2 + m4**2 - uman)*xd3u))) + 
     &              g1Lu*M1Lpu*(2*g2Lu*M2Lpu*(g3Ru*M3Lpu*xd3u + g3Lu*M0Lpu*(xd0u + xd3u)) + 
     &                 g2Ru*(2*g3Lu*M0Lpu*M3Lpu*xd3u + g3Ru*(xc01u + xc04u + 2*xc21u + M0Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u + 2*M0Lpu**2*xd3u - m2**2*xd3u + m3**2*xd3u - uman*xd3u)))))))/
     &    (32.*pi**2*(-sfmc**2 + uman)) )


      enddo

!==================================================================================================================================     
!Box(t)*s(Higgs)+Box(u)*s(Higgs)
       
      do j=sumstartH, sumendH

        Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(JLR(argx,j,argyh))
        propac = dconjg(DenSmH(j))

!t        
	Box_TS = Box_TS + 2*colorfactorSH*dreal( (Coupl_Ac*Coupl_Bc*((m1**2 + m2**2 + m3**2 - m4**2 + sqrtS**2 - tman - uman)*
     &         (g0Lt*(g1Lt*M1Lpt*(2*g2Rt*M2Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*xd2t + 
     &                 g2Lt*(2*g3Lt*M0Lpt*M3Lpt*(xd0t + xd2t) + g3Rt*(xc02t + xc04t + 2*xc22t + M0Lpt**2*xd0t + M3Lpt**2*xd0t + (2*M3Lpt**2 - m4**2 + tman)*xd2t - m2**2*(xd0t + xd2t)))) + 
     &              g1Rt*(g2Lt*M2Lpt*(2*g3Rt*xc11t + 2*g3Lt*M0Lpt*M3Lpt*xd2t + g3Rt*(2*M0Lpt**2 - m1**2 - m2**2 + sqrtS**2)*xd2t) + 
     &                 g2Rt*(g3Rt*M3Lpt*(xc03t + xc04t + 2*xc11t + M1Lpt**2*xd0t + m3**2*xd2t - tman*xd2t - m1**2*(xd0t + xd2t) + M0Lpt**2*(xd0t + 2*xd2t)) + 
     &                    g3Lt*M0Lpt*(xc02t + xc03t + 2*xc22t + M1Lpt**2*xd0t + (m3 - m4)*(m3 + m4)*xd2t - sqrtS**2*(xd0t + xd2t) + M3Lpt**2*(xd0t + 2*xd2t))))) + 
     &           g0Rt*(g1Rt*M1Lpt*(2*g2Lt*M2Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*xd2t + 
     &                 g2Rt*(2*g3Rt*M0Lpt*M3Lpt*(xd0t + xd2t) + g3Lt*(xc02t + xc04t + 2*xc22t + M0Lpt**2*xd0t + M3Lpt**2*xd0t + (2*M3Lpt**2 - m4**2 + tman)*xd2t - m2**2*(xd0t + xd2t)))) + 
     &              g1Lt*(g2Rt*M2Lpt*(2*g3Lt*xc11t + 2*g3Rt*M0Lpt*M3Lpt*xd2t + g3Lt*(2*M0Lpt**2 - m1**2 - m2**2 + sqrtS**2)*xd2t) + 
     &                 g2Lt*(g3Lt*M3Lpt*(xc03t + xc04t + 2*xc11t + M1Lpt**2*xd0t + m3**2*xd2t - tman*xd2t - m1**2*(xd0t + xd2t) + M0Lpt**2*(xd0t + 2*xd2t)) + 
     &                    g3Rt*M0Lpt*(xc02t + xc03t + 2*xc22t + M1Lpt**2*xd0t + (m3 - m4)*(m3 + m4)*xd2t - sqrtS**2*(xd0t + xd2t) + M3Lpt**2*(xd0t + 2*xd2t)))))) + 
     &        (2*m1**2 - m2**2 + m4**2 + sqrtS**2 - uman)*(-(g0Rt*(g1Rt*M1Lpt*(2*g2Lt*M2Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*(xd1t + xd2t) + 
     &                   g2Rt*(2*g3Rt*M0Lpt*M3Lpt*(xd0t + xd1t + xd2t) + g3Lt*(xc02t + xc04t + 2*xc12t + 2*xc22t + M0Lpt**2*xd0t + M3Lpt**2*xd0t + 2*M3Lpt**2*xd1t - m4**2*xd1t + tman*xd1t + 
     &                         (2*M3Lpt**2 - m4**2 + tman)*xd2t - m2**2*(xd0t + xd1t + xd2t)))) + 
     &                g1Lt*(g2Lt*(g3Lt*M3Lpt*(xc01t + xc03t - 2*xc21t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t - tman*xd1t - tman*xd2t + m3**2*(xd1t + xd2t) + 
     &                         2*M0Lpt**2*(xd0t + xd1t + xd2t) - m1**2*(xd0t + xd1t + xd2t)) + 
     &                      g3Rt*M0Lpt*(xc01t + 2*xc02t + xc03t + 2*xc12t + 2*xc22t + M1Lpt**2*xd0t + M2Lpt**2*xd0t + m3**2*(xd1t + xd2t) + (2*M3Lpt**2 - m4**2 - sqrtS**2)*(xd0t + xd1t + xd2t))
     &                      ) + g2Rt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*(xd0t + xd1t + xd2t) + 
     &                      g3Lt*(xc02t - xc04t - 2*xc21t + M3Lpt**2*xd0t - (m1 - sqrtS)*(m1 + sqrtS)*(xd1t + xd2t) - m2**2*(xd0t + xd1t + xd2t) + M0Lpt**2*(xd0t + 2*(xd1t + xd2t))))))) - 
     &           g0Lt*(g1Lt*M1Lpt*(2*g2Rt*M2Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*(xd1t + xd2t) + 
     &                 g2Lt*(2*g3Lt*M0Lpt*M3Lpt*(xd0t + xd1t + xd2t) + g3Rt*(xc02t + xc04t + 2*xc12t + 2*xc22t + M0Lpt**2*xd0t + M3Lpt**2*xd0t + 2*M3Lpt**2*xd1t - m4**2*xd1t + tman*xd1t + 
     &                       (2*M3Lpt**2 - m4**2 + tman)*xd2t - m2**2*(xd0t + xd1t + xd2t)))) + 
     &              g1Rt*(g2Rt*(g3Rt*M3Lpt*(xc01t + xc03t - 2*xc21t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t - tman*xd1t - tman*xd2t + m3**2*(xd1t + xd2t) + 
     &                       2*M0Lpt**2*(xd0t + xd1t + xd2t) - m1**2*(xd0t + xd1t + xd2t)) + 
     &                    g3Lt*M0Lpt*(xc01t + 2*xc02t + xc03t + 2*xc12t + 2*xc22t + M1Lpt**2*xd0t + M2Lpt**2*xd0t + m3**2*(xd1t + xd2t) + (2*M3Lpt**2 - m4**2 - sqrtS**2)*(xd0t + xd1t + xd2t))) 
     &                  + g2Lt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*(xd0t + xd1t + xd2t) + 
     &                    g3Rt*(xc02t - xc04t - 2*xc21t + M3Lpt**2*xd0t - (m1 - sqrtS)*(m1 + sqrtS)*(xd1t + xd2t) - m2**2*(xd0t + xd1t + xd2t) + M0Lpt**2*(xd0t + 2*(xd1t + xd2t))))))) + 
     &        (m1**2 - 2*m2**2 - m4**2 - sqrtS**2 + tman)*(g0Rt*(g1Rt*M1Lpt*(-2*g2Lt*M2Lpt*(g3Lt*M3Lpt*xd3t + g3Rt*M0Lpt*(xd0t + xd3t)) + 
     &                 g2Rt*(-2*g3Rt*M0Lpt*M3Lpt*xd3t + g3Lt*(xc01t + xc04t + (M0Lpt**2 + M2Lpt**2 - tman)*xd0t + (m2**2 - 2*M3Lpt**2 + m4**2 - tman)*xd3t))) + 
     &              g1Lt*(g2Lt*(g3Rt*M0Lpt*(xc01t + xc03t + (M1Lpt**2 + M2Lpt**2 - m3**2)*xd0t + (-m3**2 - 2*M3Lpt**2 + m4**2 + sqrtS**2)*xd3t) + 
     &                    g3Lt*M3Lpt*(-2*xc21t + (-2*M0Lpt**2 + m1**2 - m3**2 + tman)*xd3t)) - 
     &                 g2Rt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*xd3t + g3Lt*(xc03t + xc04t + 2*xc21t + M1Lpt**2*xd0t - m2**2*xd3t + sqrtS**2*xd3t - m1**2*(xd0t + xd3t) + M0Lpt**2*(xd0t + 2*xd3t))))) + 
     &           g0Lt*(g1Lt*M1Lpt*(-2*g2Rt*M2Lpt*(g3Rt*M3Lpt*xd3t + g3Lt*M0Lpt*(xd0t + xd3t)) + 
     &                 g2Lt*(-2*g3Lt*M0Lpt*M3Lpt*xd3t + g3Rt*(xc01t + xc04t + (M0Lpt**2 + M2Lpt**2 - tman)*xd0t + (m2**2 - 2*M3Lpt**2 + m4**2 - tman)*xd3t))) + 
     &              g1Rt*(g2Rt*(g3Lt*M0Lpt*(xc01t + xc03t + (M1Lpt**2 + M2Lpt**2 - m3**2)*xd0t + (-m3**2 - 2*M3Lpt**2 + m4**2 + sqrtS**2)*xd3t) + 
     &                    g3Rt*M3Lpt*(-2*xc21t + (-2*M0Lpt**2 + m1**2 - m3**2 + tman)*xd3t)) - 
     &                 g2Lt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*xd3t + g3Rt*(xc03t + xc04t + 2*xc21t + M1Lpt**2*xd0t - m2**2*xd3t + sqrtS**2*xd3t - m1**2*(xd0t + xd3t) + M0Lpt**2*(xd0t + 2*xd3t))))))))/
     &    (32.*pi**2*propac)  )
!u
	Box_US = Box_US + 2*colorfactorSH*dreal(  (Coupl_Ac*Coupl_Bc*(-((m1**2 + m2**2 + 4*m4**2 - tman - uman)*(-(g0Lu*(g1Lu*M1Lpu*
     &                   (2*g2Lu*M2Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*xd2u + g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(xd0u + xd2u) + 
     &                        g3Ru*(xc02u + xc04u + 2*xc11u + M3Lpu**2*xd0u + m3**2*xd2u - uman*xd2u - m2**2*(xd0u + xd2u) + M0Lpu**2*(xd0u + 2*xd2u)))) + 
     &                  g1Ru*(g2Ru*M2Lpu*(2*g3Ru*xc11u + 2*g3Lu*M0Lpu*M3Lpu*xd2u + g3Ru*(2*M0Lpu**2 - m1**2 - m2**2 + sqrtS**2)*xd2u) + 
     &                     g2Lu*(g3Ru*M3Lpu*(xc03u + xc04u + 2*xc12u + M0Lpu**2*xd0u + M1Lpu**2*xd0u + (2*M1Lpu**2 - m4**2 + uman)*xd2u - m1**2*(xd0u + xd2u)) + 
     &                        g3Lu*M0Lpu*(xc02u + xc03u + 2*xc12u + M3Lpu**2*xd0u + (m3 - m4)*(m3 + m4)*xd2u - sqrtS**2*(xd0u + xd2u) + M1Lpu**2*(xd0u + 2*xd2u)))))) - 
     &             g0Ru*(g1Ru*M1Lpu*(2*g2Ru*M2Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*xd2u + 
     &                   g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(xd0u + xd2u) + g3Lu*(xc02u + xc04u + 2*xc11u + M3Lpu**2*xd0u + m3**2*xd2u - uman*xd2u - m2**2*(xd0u + xd2u) + M0Lpu**2*(xd0u + 2*xd2u)))) + 
     &                g1Lu*(g2Lu*M2Lpu*(2*g3Lu*xc11u + 2*g3Ru*M0Lpu*M3Lpu*xd2u + g3Lu*(2*M0Lpu**2 - m1**2 - m2**2 + sqrtS**2)*xd2u) + 
     &                   g2Ru*(g3Lu*M3Lpu*(xc03u + xc04u + 2*xc12u + M0Lpu**2*xd0u + M1Lpu**2*xd0u + (2*M1Lpu**2 - m4**2 + uman)*xd2u - m1**2*(xd0u + xd2u)) + 
     &                      g3Ru*M0Lpu*(xc02u + xc03u + 2*xc12u + M3Lpu**2*xd0u + (m3 - m4)*(m3 + m4)*xd2u - sqrtS**2*(xd0u + xd2u) + M1Lpu**2*(xd0u + 2*xd2u))))))) + 
     &        (2*m1**2 - m2**2 + m4**2 + sqrtS**2 - uman)*(g0Ru*(-(g1Ru*M1Lpu*(2*g2Ru*M2Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*(xd1u + xd2u) + 
     &                   g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(xd0u + xd1u + xd2u) + g3Lu*(xc02u - xc04u - 2*xc21u + M3Lpu**2*xd0u + (m3**2 - uman)*(xd1u + xd2u) - m2**2*(xd0u + xd1u + xd2u) + 
     &                         M0Lpu**2*(xd0u + 2*(xd1u + xd2u)))))) + g1Lu*(-(g2Lu*M2Lpu*
     &                    (2*g3Ru*M0Lpu*M3Lpu*(xd0u + xd1u + xd2u) + g3Lu*(xc02u - xc04u - 2*xc21u + M3Lpu**2*xd0u - (m1 - sqrtS)*(m1 + sqrtS)*(xd1u + xd2u) - m2**2*(xd0u + xd1u + xd2u) + 
     &                         M0Lpu**2*(xd0u + 2*(xd1u + xd2u))))) + g2Ru*(g3Ru*M0Lpu*
     &                     (xc01u - xc03u - 2*xc12u + M2Lpu**2*xd0u - m3**2*xd0u + sqrtS**2*xd0u - m3**2*xd1u + m4**2*xd1u + sqrtS**2*xd1u - m3**2*xd2u + m4**2*xd2u + sqrtS**2*xd2u - 
     &                       M1Lpu**2*(xd0u + 2*(xd1u + xd2u))) + g3Lu*M3Lpu*(xc01u - xc03u - 2*xc12u + M2Lpu**2*xd0u - uman*xd0u - uman*xd1u - uman*xd2u + m4**2*(xd1u + xd2u) + 
     &                       m1**2*(xd0u + xd1u + xd2u) - M1Lpu**2*(xd0u + 2*(xd1u + xd2u)))))) + 
     &           g0Lu*(-(g1Lu*M1Lpu*(2*g2Lu*M2Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*(xd1u + xd2u) + 
     &                   g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(xd0u + xd1u + xd2u) + g3Ru*(xc02u - xc04u - 2*xc21u + M3Lpu**2*xd0u + (m3**2 - uman)*(xd1u + xd2u) - m2**2*(xd0u + xd1u + xd2u) + 
     &                         M0Lpu**2*(xd0u + 2*(xd1u + xd2u)))))) + g1Ru*(-(g2Ru*M2Lpu*
     &                    (2*g3Lu*M0Lpu*M3Lpu*(xd0u + xd1u + xd2u) + g3Ru*(xc02u - xc04u - 2*xc21u + M3Lpu**2*xd0u - (m1 - sqrtS)*(m1 + sqrtS)*(xd1u + xd2u) - m2**2*(xd0u + xd1u + xd2u) + 
     &                         M0Lpu**2*(xd0u + 2*(xd1u + xd2u))))) + g2Lu*(g3Lu*M0Lpu*
     &                     (xc01u - xc03u - 2*xc12u + M2Lpu**2*xd0u - m3**2*xd0u + sqrtS**2*xd0u - m3**2*xd1u + m4**2*xd1u + sqrtS**2*xd1u - m3**2*xd2u + m4**2*xd2u + sqrtS**2*xd2u - 
     &                       M1Lpu**2*(xd0u + 2*(xd1u + xd2u))) + g3Ru*M3Lpu*(xc01u - xc03u - 2*xc12u + M2Lpu**2*xd0u - uman*xd0u - uman*xd1u - uman*xd2u + m4**2*(xd1u + xd2u) + 
     &                       m1**2*(xd0u + xd1u + xd2u) - M1Lpu**2*(xd0u + 2*(xd1u + xd2u))))))) + 
     &        (m1**2 - 2*m2**2 - m4**2 - sqrtS**2 + tman)*(-(g0Lu*(g1Lu*M1Lpu*(2*g2Lu*M2Lpu*(g3Ru*M3Lpu*xd3u + g3Lu*M0Lpu*(xd0u + xd3u)) + 
     &                   g2Ru*(2*g3Lu*M0Lpu*M3Lpu*xd3u + g3Ru*(xc01u + xc04u + 2*xc21u + M2Lpu**2*xd0u + (-m2**2 + m3**2)*xd3u - uman*(xd0u + xd3u) + M0Lpu**2*(xd0u + 2*xd3u)))) + 
     &                g1Ru*(g2Ru*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*xd3u + g3Ru*(xc03u + xc04u + 2*xc21u + M1Lpu**2*xd0u - m2**2*xd3u + sqrtS**2*xd3u - m1**2*(xd0u + xd3u) + 
     &                         M0Lpu**2*(xd0u + 2*xd3u))) + g2Lu*(g3Ru*M3Lpu*(2*xc22u - (m1**2 - 2*M1Lpu**2 + m4**2 - uman)*xd3u) + 
     &                      g3Lu*M0Lpu*(xc01u + xc03u + 2*xc22u + (M2Lpu - m4)*(M2Lpu + m4)*xd0u + m3**2*xd3u - (m4**2 + sqrtS**2)*xd3u + M1Lpu**2*(xd0u + 2*xd3u)))))) - 
     &           g0Ru*(g1Ru*M1Lpu*(2*g2Ru*M2Lpu*(g3Lu*M3Lpu*xd3u + g3Ru*M0Lpu*(xd0u + xd3u)) + 
     &                 g2Lu*(2*g3Ru*M0Lpu*M3Lpu*xd3u + g3Lu*(xc01u + xc04u + 2*xc21u + M2Lpu**2*xd0u + (-m2**2 + m3**2)*xd3u - uman*(xd0u + xd3u) + M0Lpu**2*(xd0u + 2*xd3u)))) + 
     &              g1Lu*(g2Lu*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*xd3u + g3Lu*(xc03u + xc04u + 2*xc21u + M1Lpu**2*xd0u - m2**2*xd3u + sqrtS**2*xd3u - m1**2*(xd0u + xd3u) + M0Lpu**2*(xd0u + 2*xd3u))) + 
     &                 g2Ru*(g3Lu*M3Lpu*(2*xc22u - (m1**2 - 2*M1Lpu**2 + m4**2 - uman)*xd3u) + 
     &                    g3Ru*M0Lpu*(xc01u + xc03u + 2*xc22u + (M2Lpu - m4)*(M2Lpu + m4)*xd0u + m3**2*xd3u - (m4**2 + sqrtS**2)*xd3u + M1Lpu**2*(xd0u + 2*xd3u))))))))/(32.*pi**2*propac)  )
  

      enddo

!==================================================================================================================================
!Box(t)*Vector + Box(u)*Vector

      do j=sumstartV, sumendV

        Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = dconjg(LLR(j,argx,argyh))
        propac = dconjg(DenSmV(j))

!t        
        Box_TV = Box_TV + 2*colorfactorSV*dreal(   -(Coupl_Ac*Coupl_Bc*(-((m1**2 - m2**2 - tman + uman)*(-(g0Rt*(g1Lt*(g2Rt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*xd2t + g3Lt*(2*xc11t + (2*M0Lpt**2 - m1**2 - m2**2 + sqrtS**2)*xd2t)) + 
     &                      g2Lt*(g3Rt*M0Lpt*(xc02t + xc03t + 2*xc22t + M1Lpt**2*xd0t + M3Lpt**2*xd0t - sqrtS**2*xd0t + m3**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t - sqrtS**2*xd2t) + 
     &                         g3Lt*M3Lpt*(xc03t + xc04t + 2*xc11t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t))) + 
     &                   g1Rt*M1Lpt*(2*g2Lt*M2Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*xd2t + 
     &                      g2Rt*(2*g3Rt*M0Lpt*M3Lpt*(xd0t + xd2t) + g3Lt*(xc02t + xc04t + 2*xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t - m2**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t + 
     &                            tman*xd2t))))) - g0Lt*(g1Rt*(g2Lt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*xd2t + g3Rt*(2*xc11t + (2*M0Lpt**2 - m1**2 - m2**2 + sqrtS**2)*xd2t)) + 
     &                    g2Rt*(g3Lt*M0Lpt*(xc02t + xc03t + 2*xc22t + M1Lpt**2*xd0t + M3Lpt**2*xd0t - sqrtS**2*xd0t + m3**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t - sqrtS**2*xd2t) + 
     &                       g3Rt*M3Lpt*(xc03t + xc04t + 2*xc11t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t))) + 
     &                 g1Lt*M1Lpt*(2*g2Rt*M2Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*xd2t + 
     &                    g2Lt*(2*g3Lt*M0Lpt*M3Lpt*(xd0t + xd2t) + g3Rt*(xc02t + xc04t + 2*xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t - m2**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t + 
     &                          tman*xd2t)))))) + (3*m1**2 + m2**2 - sqrtS**2)*
     &          (-(g0Rt*(g1Lt*(g2Rt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*(xd0t + xd1t + xd2t) + 
     &                       g3Lt*(xc02t - xc04t - 2*xc21t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 2*M0Lpt**2*xd1t - m1**2*xd1t - m2**2*xd1t + sqrtS**2*xd1t + 2*M0Lpt**2*xd2t - 
     &                          m1**2*xd2t - m2**2*xd2t + sqrtS**2*xd2t)) + g2Lt*
     &                     (g3Rt*M0Lpt*(xc01t + 2*xc02t + xc03t + 2*xc12t + 2*xc22t + M1Lpt**2*xd0t + M2Lpt**2*xd0t + 2*M3Lpt**2*xd0t - m4**2*xd0t - sqrtS**2*xd0t + m3**2*xd1t + 
     &                          2*M3Lpt**2*xd1t - m4**2*xd1t - sqrtS**2*xd1t + m3**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t - sqrtS**2*xd2t) + 
     &                       g3Lt*M3Lpt*(xc01t + xc03t - 2*xc21t + 2*M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t + 2*M0Lpt**2*xd1t - m1**2*xd1t + m3**2*xd1t - 
     &                          tman*xd1t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t))) + 
     &                 g1Rt*M1Lpt*(2*g2Lt*M2Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*(xd1t + xd2t) + 
     &                    g2Rt*(2*g3Rt*M0Lpt*M3Lpt*(xd0t + xd1t + xd2t) + g3Lt*
     &                        (xc02t + xc04t + 2*xc12t + 2*xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t - m2**2*xd1t + 2*M3Lpt**2*xd1t - m4**2*xd1t + tman*xd1t - m2**2*xd2t + 
     &                          2*M3Lpt**2*xd2t - m4**2*xd2t + tman*xd2t))))) - 
     &            g0Lt*(g1Rt*(g2Lt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*(xd0t + xd1t + xd2t) + 
     &                     g3Rt*(xc02t - xc04t - 2*xc21t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 2*M0Lpt**2*xd1t - m1**2*xd1t - m2**2*xd1t + sqrtS**2*xd1t + 2*M0Lpt**2*xd2t - 
     &                        m1**2*xd2t - m2**2*xd2t + sqrtS**2*xd2t)) + g2Rt*
     &                   (g3Lt*M0Lpt*(xc01t + 2*xc02t + xc03t + 2*xc12t + 2*xc22t + M1Lpt**2*xd0t + M2Lpt**2*xd0t + 2*M3Lpt**2*xd0t - m4**2*xd0t - sqrtS**2*xd0t + m3**2*xd1t + 
     &                        2*M3Lpt**2*xd1t - m4**2*xd1t - sqrtS**2*xd1t + m3**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t - sqrtS**2*xd2t) + 
     &                     g3Rt*M3Lpt*(xc01t + xc03t - 2*xc21t + 2*M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t + 2*M0Lpt**2*xd1t - m1**2*xd1t + m3**2*xd1t - 
     &                        tman*xd1t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t))) + 
     &               g1Lt*M1Lpt*(2*g2Rt*M2Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*(xd1t + xd2t) + 
     &                  g2Lt*(2*g3Lt*M0Lpt*M3Lpt*(xd0t + xd1t + xd2t) + g3Rt*(xc02t + xc04t + 2*xc12t + 2*xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t - m2**2*xd1t + 2*M3Lpt**2*xd1t - 
     &                        m4**2*xd1t + tman*xd1t - m2**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t + tman*xd2t))))) + 
     &         (m1**2 + 3*m2**2 - sqrtS**2)*(g0Rt*(g1Rt*M1Lpt*(-2*g2Lt*M2Lpt*(g3Lt*M3Lpt*xd3t + g3Rt*M0Lpt*(xd0t + xd3t)) + 
     &                  g2Rt*(-2*g3Rt*M0Lpt*M3Lpt*xd3t + g3Lt*(xc01t + xc04t + M0Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t + m2**2*xd3t - 2*M3Lpt**2*xd3t + m4**2*xd3t - tman*xd3t))) + 
     &               g1Lt*(-(g2Rt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*xd3t + g3Lt*(xc03t + xc04t + 2*xc21t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 2*M0Lpt**2*xd3t - m1**2*xd3t - m2**2*xd3t + 
     &                          sqrtS**2*xd3t))) + g2Lt*(g3Rt*M0Lpt*(xc01t + xc03t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - m3**2*xd0t - m3**2*xd3t - 2*M3Lpt**2*xd3t + m4**2*xd3t + sqrtS**2*xd3t) + 
     &                     g3Lt*M3Lpt*(-2*xc21t + (-2*M0Lpt**2 + m1**2 - m3**2 + tman)*xd3t)))) + 
     &            g0Lt*(g1Lt*M1Lpt*(-2*g2Rt*M2Lpt*(g3Rt*M3Lpt*xd3t + g3Lt*M0Lpt*(xd0t + xd3t)) + 
     &                  g2Lt*(-2*g3Lt*M0Lpt*M3Lpt*xd3t + g3Rt*(xc01t + xc04t + M0Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t + m2**2*xd3t - 2*M3Lpt**2*xd3t + m4**2*xd3t - tman*xd3t))) + 
     &               g1Rt*(-(g2Lt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*xd3t + g3Rt*(xc03t + xc04t + 2*xc21t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 2*M0Lpt**2*xd3t - m1**2*xd3t - m2**2*xd3t + 
     &                          sqrtS**2*xd3t))) + g2Rt*(g3Lt*M0Lpt*(xc01t + xc03t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - m3**2*xd0t - m3**2*xd3t - 2*M3Lpt**2*xd3t + m4**2*xd3t + sqrtS**2*xd3t) + 
     &                     g3Rt*M3Lpt*(-2*xc21t + (-2*M0Lpt**2 + m1**2 - m3**2 + tman)*xd3t)))))))/(32.*pi**2*propac)   )
!u        
        Box_UV = Box_UV + 2*colorfactorSV*dreal(    -(Coupl_Ac*Coupl_Bc*((3*m1**2 + m2**2 - sqrtS**2)*(g0Ru*(-(g1Ru*M1Lpu*(2*g2Ru*M2Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*(xd1u + xd2u) + 
     &                    g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(xd0u + xd1u + xd2u) + g3Lu*
     &                        (xc02u - xc04u - 2*xc21u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd1u - m2**2*xd1u + m3**2*xd1u - uman*xd1u + 2*M0Lpu**2*xd2u - m2**2*xd2u + 
     &                          m3**2*xd2u - uman*xd2u)))) + g1Lu*(-(g2Lu*M2Lpu*
     &                     (2*g3Ru*M0Lpu*M3Lpu*(xd0u + xd1u + xd2u) + g3Lu*(xc02u - xc04u - 2*xc21u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd1u - m1**2*xd1u - m2**2*xd1u + 
     &                          sqrtS**2*xd1u + 2*M0Lpu**2*xd2u - m1**2*xd2u - m2**2*xd2u + sqrtS**2*xd2u))) + 
     &                  g2Ru*(g3Ru*M0Lpu*(xc01u - xc03u - 2*xc12u - M1Lpu**2*xd0u + M2Lpu**2*xd0u - m3**2*xd0u + sqrtS**2*xd0u - 2*M1Lpu**2*xd1u - m3**2*xd1u + m4**2*xd1u + sqrtS**2*xd1u - 
     &                        2*M1Lpu**2*xd2u - m3**2*xd2u + m4**2*xd2u + sqrtS**2*xd2u) + 
     &                     g3Lu*M3Lpu*(xc01u - xc03u - 2*xc12u + m1**2*xd0u - M1Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u + m1**2*xd1u - 2*M1Lpu**2*xd1u + m4**2*xd1u - uman*xd1u + m1**2*xd2u - 
     &                        2*M1Lpu**2*xd2u + m4**2*xd2u - uman*xd2u)))) + g0Lu*
     &             (-(g1Lu*M1Lpu*(2*g2Lu*M2Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*(xd1u + xd2u) + 
     &                    g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(xd0u + xd1u + xd2u) + g3Ru*
     &                        (xc02u - xc04u - 2*xc21u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd1u - m2**2*xd1u + m3**2*xd1u - uman*xd1u + 2*M0Lpu**2*xd2u - m2**2*xd2u + 
     &                          m3**2*xd2u - uman*xd2u)))) + g1Ru*(-(g2Ru*M2Lpu*
     &                     (2*g3Lu*M0Lpu*M3Lpu*(xd0u + xd1u + xd2u) + g3Ru*(xc02u - xc04u - 2*xc21u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd1u - m1**2*xd1u - m2**2*xd1u + 
     &                          sqrtS**2*xd1u + 2*M0Lpu**2*xd2u - m1**2*xd2u - m2**2*xd2u + sqrtS**2*xd2u))) + 
     &                  g2Lu*(g3Lu*M0Lpu*(xc01u - xc03u - 2*xc12u - M1Lpu**2*xd0u + M2Lpu**2*xd0u - m3**2*xd0u + sqrtS**2*xd0u - 2*M1Lpu**2*xd1u - m3**2*xd1u + m4**2*xd1u + sqrtS**2*xd1u - 
     &                        2*M1Lpu**2*xd2u - m3**2*xd2u + m4**2*xd2u + sqrtS**2*xd2u) + 
     &                     g3Ru*M3Lpu*(xc01u - xc03u - 2*xc12u + m1**2*xd0u - M1Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u + m1**2*xd1u - 2*M1Lpu**2*xd1u + m4**2*xd1u - uman*xd1u + m1**2*xd2u - 
     &                        2*M1Lpu**2*xd2u + m4**2*xd2u - uman*xd2u))))) + (-m1**2 + m2**2 - tman + uman)*
     &          (-(g0Ru*(g1Ru*M1Lpu*(2*g2Ru*M2Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*xd2u + 
     &                    g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(xd0u + xd2u) + g3Lu*(xc02u + xc04u + 2*xc11u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd2u - m2**2*xd2u + m3**2*xd2u - 
     &                          uman*xd2u))) + g1Lu*(g2Lu*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*xd2u + g3Lu*(2*xc11u + (2*M0Lpu**2 - m1**2 - m2**2 + sqrtS**2)*xd2u)) + 
     &                    g2Ru*(g3Ru*M0Lpu*(xc02u + xc03u + 2*xc12u + M1Lpu**2*xd0u + M3Lpu**2*xd0u - sqrtS**2*xd0u + 2*M1Lpu**2*xd2u + m3**2*xd2u - m4**2*xd2u - sqrtS**2*xd2u) + 
     &                       g3Lu*M3Lpu*(xc03u + xc04u + 2*xc12u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u - m1**2*xd2u + 2*M1Lpu**2*xd2u - m4**2*xd2u + uman*xd2u))))) - 
     &            g0Lu*(g1Lu*M1Lpu*(2*g2Lu*M2Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*xd2u + 
     &                  g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(xd0u + xd2u) + g3Ru*(xc02u + xc04u + 2*xc11u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd2u - m2**2*xd2u + m3**2*xd2u - 
     &                        uman*xd2u))) + g1Ru*(g2Ru*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*xd2u + g3Ru*(2*xc11u + (2*M0Lpu**2 - m1**2 - m2**2 + sqrtS**2)*xd2u)) + 
     &                  g2Lu*(g3Lu*M0Lpu*(xc02u + xc03u + 2*xc12u + M1Lpu**2*xd0u + M3Lpu**2*xd0u - sqrtS**2*xd0u + 2*M1Lpu**2*xd2u + m3**2*xd2u - m4**2*xd2u - sqrtS**2*xd2u) + 
     &                     g3Ru*M3Lpu*(xc03u + xc04u + 2*xc12u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u - m1**2*xd2u + 2*M1Lpu**2*xd2u - m4**2*xd2u + uman*xd2u))))) + 
     &         (m1**2 + 3*m2**2 - sqrtS**2)*(-(g0Ru*(g1Lu*(g2Lu*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*xd3u + 
     &                       g3Lu*(xc03u + xc04u + 2*xc21u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + 2*M0Lpu**2*xd3u - m1**2*xd3u - m2**2*xd3u + sqrtS**2*xd3u)) + 
     &                    g2Ru*(g3Ru*M0Lpu*(xc01u + xc03u + 2*xc22u + M1Lpu**2*xd0u + M2Lpu**2*xd0u - m4**2*xd0u + 2*M1Lpu**2*xd3u + m3**2*xd3u - m4**2*xd3u - sqrtS**2*xd3u) + 
     &                       g3Lu*M3Lpu*(2*xc22u - (m1**2 - 2*M1Lpu**2 + m4**2 - uman)*xd3u))) + 
     &                 g1Ru*M1Lpu*(2*g2Ru*M2Lpu*(g3Lu*M3Lpu*xd3u + g3Ru*M0Lpu*(xd0u + xd3u)) + 
     &                    g2Lu*(2*g3Ru*M0Lpu*M3Lpu*xd3u + g3Lu*(xc01u + xc04u + 2*xc21u + M0Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u + 2*M0Lpu**2*xd3u - m2**2*xd3u + m3**2*xd3u - uman*xd3u))))) 
     &             -   g0Lu*(g1Ru*(g2Ru*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*xd3u + g3Ru*(xc03u + xc04u + 2*xc21u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + 2*M0Lpu**2*xd3u - m1**2*xd3u - m2**2*xd3u + 
     &                        sqrtS**2*xd3u)) + g2Lu*(g3Lu*M0Lpu*(xc01u + xc03u + 2*xc22u + M1Lpu**2*xd0u + M2Lpu**2*xd0u - m4**2*xd0u + 2*M1Lpu**2*xd3u + m3**2*xd3u - m4**2*xd3u - 
     &                        sqrtS**2*xd3u) + g3Ru*M3Lpu*(2*xc22u - (m1**2 - 2*M1Lpu**2 + m4**2 - uman)*xd3u))) + 
     &               g1Lu*M1Lpu*(2*g2Lu*M2Lpu*(g3Ru*M3Lpu*xd3u + g3Lu*M0Lpu*(xd0u + xd3u)) + 
     &                  g2Ru*(2*g3Lu*M0Lpu*M3Lpu*xd3u + g3Ru*(xc01u + xc04u + 2*xc21u + M0Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u + 2*M0Lpu**2*xd3u - m2**2*xd3u + m3**2*xd3u - uman*xd3u)))))))/
     &    (32.*pi**2*propac)    )

      enddo           
 
!==================================================================================================================================

      VH_TX(1) = Box_TQ
      VH_TX(2) = Box_TS
      VH_TX(3) = Box_TV
      VH_TX(4) = Box_TT
      VH_TX(5) = Box_TU

      VH_TX(6) = VH_TX(1) + VH_TX(2) + VH_TX(3) + VH_TX(4) + VH_TX(5)

      VH_UX(1) = Box_UQ
      VH_UX(2) = Box_US
      VH_UX(3) = Box_UV
      VH_UX(4) = Box_UT
      VH_UX(5) = Box_UU

      VH_UX(6) = VH_UX(1) + VH_UX(2) + VH_UX(3) + VH_UX(4) + VH_UX(5)

!===================================================================================================================================

      end
