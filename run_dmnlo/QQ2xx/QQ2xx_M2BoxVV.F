!**********************************************************************************************************************************

      subroutine QQ2xx_M2GluonBoxVV(argx,argy,VV_TX,VV_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"
#include "../util/DMNLO_LoopIntegrals.h"

      integer i,j,argx,argy
      complex*16 g0t,g1t,g2t,g3t,g0u,g1u,g2u,g3u,g0Lt,g0Rt,g1Lt,g1Rt,g2Lt,g2Rt,g3Lt,g3Rt,g0Lu,g0Ru,g1Lu,g1Ru,g2Lu,g2Ru,g3Lu,g3Ru,
     &           xb01t,xb02t,xc01t,xc02t,xc03t,xc04t,xc11t,xc12t,xc13t,xc21t,xc22t,xc23t,xc001t,xc002t,xc003t,xc111t,xc112t,xc113t,xc121t,
     &           xc122t,xc123t,xc221t,xc222t,xc223t,xd0t,xd1t,xd2t,xd3t,xd00t,xd11t,xd12t,xd13t,xd22t,xd23t,xd33t,xb01u,xb02u,xc01u,xc02u,
     &           xc03u,xc04u,xc11u,xc12u,xc13u,xc21u,xc22u,xc23u,xc001u,xc002u,xc003u,xc111u,xc112u,xc113u,xc121u,xc122u,xc123u,xc221u,xc222u,
     &           xc223u,xd0u,xd1u,xd2u,xd3u,xd00u,xd11u,xd12u,xd13u,xd22u,xd23u,xd33u
      double precision M0Lpt,M1Lpt,M2Lpt,M3Lpt,M0Lpu,M1Lpu,M2Lpu,M3Lpu,VV_TX(6),VV_UX(6),Box_TQ,Box_TS,Box_TV,Box_TT,Box_TU,Box_UQ,Box_US,Box_UV,Box_UT,Box_UU
        
      external myBget,myCget,myDget
  
! ==================================================================================================================================
!                             ---                 Higgs-Vector-final-states                    ---
! ==================================================================================================================================

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

! =================================================================================================================================	
!                             ---                     Vector-Vector-final-states                    ---
! =================================================================================================================================

       do i= 1, 2

!=====================================================!Box(t)_Gluon================================================================

        M0Lpt=MV(Looptype0)
        M1Lpt=MSf(Loopmeigen1,Looptype1,Loopgen1)
        M2Lpt=MSf(i,Looptypet,Loopgent)
        M3Lpt=MSf(Loopmeigen2,Looptype2,Loopgen2)

        g0t=VSfSf(Looptype0,Loopmeigen1,sfmeigen1,Looptype1,sftype1,Loopgen1,sfgen1)
        g1t=VSfSft(xtype,i,Loopmeigen1,Looptypet,Looptype1,Loopgent,Loopgen1)
        g2t=VSfSfu(ytype,Loopmeigen2,i,Looptype2,Looptypet,Loopgen2,Loopgent)
        g3t=VSfSf(Looptype0,sfmeigen2,Loopmeigen2,sftype2,Looptype2,sfgen2,Loopgen2)

!=====================================================!Box(u)_Gluon================================================================
       
        M0Lpu=MV(Looptype0)
        M1Lpu=MSf(Loopmeigen1,Looptype1,Loopgen1)
        M2Lpu=MSf(i,Looptypeu,Loopgenu)
        M3Lpu=MSf(Loopmeigen2,Looptype2,Loopgen2)

        g0u=VSfSf(Looptype0,Loopmeigen1,sfmeigen1,Looptype1,sftype1,Loopgen1,sfgen1)
        g1u=VSfSfu(ytype,i,Loopmeigen1,Looptypeu,Looptype1,Loopgenu,Loopgen1)
        g2u=VSfSft(xtype,Loopmeigen2,i,Looptype2,Looptypeu,Loopgen2,Loopgenu)
        g3u=VSfSf(Looptype0,sfmeigen2,Loopmeigen2,sftype2,Looptype2,sfgen2,Loopgen2)

!=====================================================!call Loopfunctions==========================================================

!t
        call myCget(m1**2,m3**2,tman,M0Lpt**2,M1Lpt**2,M2Lpt**2)
       
        xc01t = xc0    
        xc12t = xc1
        xc21t = xc2
        xc003t = xc00
        xc113t = xc11
        xc123t = xc12
        xc223t = xc22

        call myCget(tman,m4**2,m2**2,M0Lpt**2,M2Lpt**2,M3Lpt**2)

        xc02t = xc0
        xc11t = xc1
        xc22t = xc2
        xc002t = xc00
        xc111t = xc11
        xc122t = xc12
        xc222t = xc22

        call myCget(m3**2,m4**2,sqrtS**2,M1Lpt**2,M2Lpt**2,M3Lpt**2)

        xc03t = xc0
        xc13t = xc1
        xc23t = xc2
        xc001t = xc00
        xc112t = xc11
        xc121t = xc12
        xc221t = xc22

        call myDget(m1**2, m3**2, m4**2, m2**2, tman, sqrtS**2, M0Lpt**2, M1Lpt**2, M2Lpt**2, M3Lpt**2)

        xd0t = xd0
        xd1t = xd1
        xd2t = xd2
        xd3t = xd3
        xd00t = xd00
        xd11t = xd11
        xd12t = xd12  
        xd13t = xd13     
        xd22t = xd22
        xd23t = xd23 
        xd33t = xd33

!u
        call myCget(m1**2,m4**2,uman,M0Lpu**2,M1Lpu**2,M2Lpu**2)
       
        xc01u = xc0    
        xc12u = xc1
        xc21u = xc2
        xc003u = xc00
        xc113u = xc11
        xc123u = xc12
        xc223u = xc22

        call myCget(uman,m3**2,m2**2,M0Lpu**2,M2Lpu**2,M3Lpu**2)

        xc02u = xc0
        xc11u = xc1
        xc22u = xc2
        xc002u = xc00
        xc111u = xc11
        xc122u = xc12
        xc222u = xc22

        call myCget(m4**2,m3**2,sqrtS**2,M1Lpu**2,M2Lpu**2,M3Lpu**2)

        xc03u = xc0
        xc13u = xc1
        xc23u = xc2
        xc001u = xc00
        xc112u = xc11
        xc121u = xc12
        xc221u = xc22

        call myDget(m1**2, m4**2, m3**2, m2**2, uman, sqrtS**2, M0Lpu**2, M1Lpu**2, M2Lpu**2, M3Lpu**2)

        xd0u = xd0
        xd1u = xd1
        xd2u = xd2
        xd3u = xd3
        xd00u = xd00
        xd11u = xd11
        xd12u = xd12  
        xd13u = xd13     
        xd22u = xd22
        xd23u = xd23 
        xd33u = xd33

!==================================================================================================================================

        do j=1,  2   

!Box(t)*t + Box(u)*t
       
        Coupl_Ac = dconjg(VSfSft(xtype,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))       
        Coupl_Bc = dconjg(VSfSfu(ytype,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
        sfmc = MSf(j,proptypetc,propgentc)
        
!t
	Box_TT = Box_TT + 2*colorfactorTU*dreal( -(Coupl_Ac*Coupl_Bc*g0t*g1t*g2t*g3t*(-4*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*uman)*
     &        (-xc001t + xc002t + xc003t - M0Lpt**2*xd00t + m1**2*xd00t + M1Lpt**2*xd00t + m2**2*xd00t + M3Lpt**2*xd00t - 2*sqrtS**2*xd00t) + 
     &       (m1**2 - tman)*((2*m2**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*uman)*
     &           (xc01t + xc02t - xc03t + 4*xc111t - 4*xc112t + 4*xc11t - 4*xc13t + 4*xc21t + 4*xc223t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 
     &             2*sqrtS**2*xd0t - 4*M0Lpt**2*xd22t + 4*m1**2*xd22t + 4*M1Lpt**2*xd22t + 4*m2**2*xd22t + 4*M3Lpt**2*xd22t - 8*sqrtS**2*xd22t - 4*M0Lpt**2*xd2t + 4*m1**2*xd2t + 
     &             4*M1Lpt**2*xd2t + 4*m2**2*xd2t + 4*M3Lpt**2*xd2t - 8*sqrtS**2*xd2t) + 
     &          (3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*(xc01t + xc02t + xc03t + 4*xc111t + 4*xc11t + 4*xc121t + 4*xc123t + 2*xc12t + 2*xc13t + 4*xc21t + 4*xc223t + 2*xc23t - 
     &             M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 4*M0Lpt**2*xd12t + 4*m1**2*xd12t + 4*M1Lpt**2*xd12t + 4*m2**2*xd12t + 
     &             4*M3Lpt**2*xd12t - 8*sqrtS**2*xd12t - 2*M0Lpt**2*xd1t + 2*m1**2*xd1t + 2*M1Lpt**2*xd1t + 2*m2**2*xd1t + 2*M3Lpt**2*xd1t - 4*sqrtS**2*xd1t - 4*M0Lpt**2*xd22t + 
     &             4*m1**2*xd22t + 4*M1Lpt**2*xd22t + 4*m2**2*xd22t + 4*M3Lpt**2*xd22t - 8*sqrtS**2*xd22t - 4*M0Lpt**2*xd2t + 4*m1**2*xd2t + 4*M1Lpt**2*xd2t + 4*m2**2*xd2t + 
     &             4*M3Lpt**2*xd2t - 8*sqrtS**2*xd2t) + (3*m2**2 - m4**2 + tman)*
     &           (xc01t + xc02t - xc03t + 2*xc11t - 4*xc121t + 4*xc122t - 2*xc13t + 2*xc21t + 2*xc22t - 2*xc23t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 
     &             2*sqrtS**2*xd0t - 4*M0Lpt**2*xd23t + 4*m1**2*xd23t + 4*M1Lpt**2*xd23t + 4*m2**2*xd23t + 4*M3Lpt**2*xd23t - 8*sqrtS**2*xd23t - 2*M0Lpt**2*xd2t + 2*m1**2*xd2t + 
     &             2*M1Lpt**2*xd2t + 2*m2**2*xd2t + 2*M3Lpt**2*xd2t - 4*sqrtS**2*xd2t - 2*M0Lpt**2*xd3t + 2*m1**2*xd3t + 2*M1Lpt**2*xd3t + 2*m2**2*xd3t + 2*M3Lpt**2*xd3t - 4*sqrtS**2*xd3t)) 
     &        + (3*m1**2 - m3**2 + tman)*(-((3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*
     &             (xc01t + xc02t + 2*xc111t + 2*xc113t + 3*xc11t + 4*xc123t + 3*xc12t + 3*xc21t - 2*xc221t + 2*xc223t - xc23t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + 
     &               M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 2*M0Lpt**2*xd11t + 2*m1**2*xd11t + 2*M1Lpt**2*xd11t + 2*m2**2*xd11t + 2*M3Lpt**2*xd11t - 4*sqrtS**2*xd11t - 4*M0Lpt**2*xd12t + 
     &               4*m1**2*xd12t + 4*M1Lpt**2*xd12t + 4*m2**2*xd12t + 4*M3Lpt**2*xd12t - 8*sqrtS**2*xd12t - 3*M0Lpt**2*xd1t + 3*m1**2*xd1t + 3*M1Lpt**2*xd1t + 3*m2**2*xd1t + 
     &               3*M3Lpt**2*xd1t - 6*sqrtS**2*xd1t - 2*M0Lpt**2*xd22t + 2*m1**2*xd22t + 2*M1Lpt**2*xd22t + 2*m2**2*xd22t + 2*M3Lpt**2*xd22t - 4*sqrtS**2*xd22t - 3*M0Lpt**2*xd2t + 
     &               3*m1**2*xd2t + 3*M1Lpt**2*xd2t + 3*m2**2*xd2t + 3*M3Lpt**2*xd2t - 6*sqrtS**2*xd2t)) - 
     &          (2*m2**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*uman)*(xc01t + xc02t + 2*xc111t + 3*xc11t + 2*xc121t + 2*xc123t + xc12t + 3*xc21t + 2*xc223t + xc23t - M0Lpt**2*xd0t + m1**2*xd0t + 
     &             M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 2*M0Lpt**2*xd12t + 2*m1**2*xd12t + 2*M1Lpt**2*xd12t + 2*m2**2*xd12t + 2*M3Lpt**2*xd12t - 4*sqrtS**2*xd12t - 
     &             M0Lpt**2*xd1t + m1**2*xd1t + M1Lpt**2*xd1t + m2**2*xd1t + M3Lpt**2*xd1t - 2*sqrtS**2*xd1t - 2*M0Lpt**2*xd22t + 2*m1**2*xd22t + 2*M1Lpt**2*xd22t + 2*m2**2*xd22t + 
     &             2*M3Lpt**2*xd22t - 4*sqrtS**2*xd22t - 3*M0Lpt**2*xd2t + 3*m1**2*xd2t + 3*M1Lpt**2*xd2t + 3*m2**2*xd2t + 3*M3Lpt**2*xd2t - 6*sqrtS**2*xd2t) - 
     &          (3*m2**2 - m4**2 + tman)*(xc01t + xc02t + xc11t + 2*xc122t + xc12t + xc21t + 2*xc221t + 2*xc22t + xc23t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + 
     &             M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 2*M0Lpt**2*xd13t + 2*m1**2*xd13t + 2*M1Lpt**2*xd13t + 2*m2**2*xd13t + 2*M3Lpt**2*xd13t - 4*sqrtS**2*xd13t - M0Lpt**2*xd1t + m1**2*xd1t + 
     &             M1Lpt**2*xd1t + m2**2*xd1t + M3Lpt**2*xd1t - 2*sqrtS**2*xd1t - 2*M0Lpt**2*xd23t + 2*m1**2*xd23t + 2*M1Lpt**2*xd23t + 2*m2**2*xd23t + 2*M3Lpt**2*xd23t - 4*sqrtS**2*xd23t - 
     &             M0Lpt**2*xd2t + m1**2*xd2t + M1Lpt**2*xd2t + m2**2*xd2t + M3Lpt**2*xd2t - 2*sqrtS**2*xd2t - 2*M0Lpt**2*xd3t + 2*m1**2*xd3t + 2*M1Lpt**2*xd3t + 2*m2**2*xd3t + 
     &             2*M3Lpt**2*xd3t - 4*sqrtS**2*xd3t)) + (-2*m1**2 - 3*m2**2 - m3**2 + 2*sqrtS**2 + uman)*
     &        ((2*m2**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*uman)*(-2*xc121t + 2*xc122t + xc22t - xc23t - 2*M0Lpt**2*xd23t + 2*m1**2*xd23t + 2*M1Lpt**2*xd23t + 2*m2**2*xd23t + 
     &             2*M3Lpt**2*xd23t - 4*sqrtS**2*xd23t - M0Lpt**2*xd3t + m1**2*xd3t + M1Lpt**2*xd3t + m2**2*xd3t + M3Lpt**2*xd3t - 2*sqrtS**2*xd3t) + 
     &          (3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*(2*xc122t + 2*xc221t + xc22t + xc23t - 2*M0Lpt**2*xd13t + 2*m1**2*xd13t + 2*M1Lpt**2*xd13t + 2*m2**2*xd13t + 
     &             2*M3Lpt**2*xd13t - 4*sqrtS**2*xd13t - 2*M0Lpt**2*xd23t + 2*m1**2*xd23t + 2*M1Lpt**2*xd23t + 2*m2**2*xd23t + 2*M3Lpt**2*xd23t - 4*sqrtS**2*xd23t - M0Lpt**2*xd3t + 
     &             m1**2*xd3t + M1Lpt**2*xd3t + m2**2*xd3t + M3Lpt**2*xd3t - 2*sqrtS**2*xd3t) + 
     &          (3*m2**2 - m4**2 + tman)*(-2*xc221t + 2*xc222t + xc22t - xc23t - 2*M0Lpt**2*xd33t + 2*m1**2*xd33t + 2*M1Lpt**2*xd33t + 2*m2**2*xd33t + 2*M3Lpt**2*xd33t - 4*sqrtS**2*xd33t - 
     &             M0Lpt**2*xd3t + m1**2*xd3t + M1Lpt**2*xd3t + m2**2*xd3t + M3Lpt**2*xd3t - 2*sqrtS**2*xd3t))))/(32.*pi**2*(-sfmc**2 + tman)) )

!u
	Box_UT = Box_UT + 2*colorfactorTU*dreal(   -(Coupl_Ac*Coupl_Bc*g0u*g1u*g2u*g3u*(-((m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xc01u) - (m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xc01u - 
     &       (3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc01u - 
     &       (2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc01u - (m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xc02u - 
     &       (m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xc02u - (3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc02u - 
     &       (2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc02u - (m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xc03u + 
     &       (m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xc03u - 4*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xc111u - 
     &       2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc111u - 2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc113u - 
     &       4*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xc11u - 2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xc11u - 
     &       3*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc11u - 
     &       (2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc11u - 4*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xc121u + 
     &       4*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xc121u - 2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xc122u - 
     &       4*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xc122u - 
     &       2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc122u - 4*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xc123u - 
     &       4*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc123u - 2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xc12u - 
     &       3*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc12u - 
     &       (2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc12u - 2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xc13u + 
     &       2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xc13u - 4*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xc21u - 
     &       2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xc21u - 3*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc21u - 
     &       (2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc21u - 2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xc221u + 
     &       2*(3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xc221u + 2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc221u - 
     &       2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc221u - 
     &       2*(3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xc222u - 4*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xc223u - 
     &       2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc223u - (3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xc22u - 
     &       2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xc22u - (3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xc22u - 
     &       2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc22u - 2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xc23u - 
     &       (3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xc23u + 2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xc23u + 
     &       (3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xc23u + (3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc23u - 
     &       (2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xc23u - 
     &       4*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*uman)*
     &        (-xc001u + xc002u + xc003u - M0Lpu**2*xd00u + m1**2*xd00u + M1Lpu**2*xd00u + m2**2*xd00u + M3Lpu**2*xd00u - 2*sqrtS**2*xd00u) + 
     &       M0Lpu**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd0u - m1**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd0u - M1Lpu**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd0u - 
     &       m2**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd0u - M3Lpu**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd0u + 2*sqrtS**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd0u + 
     &       M0Lpu**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd0u - m1**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd0u - 
     &       M1Lpu**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd0u - m2**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd0u - 
     &       M3Lpu**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd0u + M0Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd0u - 
     &       m1**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd0u - M1Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd0u - 
     &       m2**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd0u - M3Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd0u + 
     &       M0Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd0u - 
     &       m1**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd0u - 
     &       M1Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd0u - 
     &       m2**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd0u - 
     &       M3Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd0u - 
     &       2*sqrtS**2*(m2**2 - tman)*(-2*m1**2 - 3*m2**2 - m3**2 + 2*sqrtS**2 + uman)*xd0u - 2*sqrtS**2*(3*m1**2 - m3**2 + tman)*(-3*m1**2 - 2*m2**2 - m4**2 + 2*sqrtS**2 + uman)*xd0u + 
     &       2*sqrtS**2*(-2*m1**2 - 3*m2**2 - m3**2 + 2*sqrtS**2 + uman)*(-3*m1**2 - 2*m2**2 - m4**2 + 2*sqrtS**2 + uman)*xd0u + 
     &       2*M0Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd11u - 
     &       2*m1**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd11u - 
     &       2*M1Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd11u - 
     &       2*m2**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd11u - 
     &       2*M3Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd11u - 
     &       4*sqrtS**2*(3*m1**2 - m3**2 + tman)*(-3*m1**2 - 2*m2**2 - m4**2 + 2*sqrtS**2 + uman)*xd11u + 4*M0Lpu**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd12u - 
     &       4*m1**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd12u - 4*M1Lpu**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd12u - 4*m2**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd12u - 
     &       4*M3Lpu**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd12u + 8*sqrtS**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd12u + 
     &       4*M0Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd12u - 
     &       4*m1**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd12u - 
     &       4*M1Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd12u - 
     &       4*m2**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd12u - 
     &       4*M3Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd12u - 
     &       8*sqrtS**2*(3*m1**2 - m3**2 + tman)*(-3*m1**2 - 2*m2**2 - m4**2 + 2*sqrtS**2 + uman)*xd12u + 2*M0Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd13u - 
     &       2*m1**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd13u - 2*M1Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd13u - 
     &       2*m2**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd13u - 2*M3Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd13u + 
     &       4*sqrtS**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd13u + 
     &       2*M0Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd13u - 
     &       2*m1**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd13u - 
     &       2*M1Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd13u - 
     &       2*m2**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd13u - 
     &       2*M3Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd13u + 
     &       4*sqrtS**2*(-2*m1**2 - 3*m2**2 - m3**2 + 2*sqrtS**2 + uman)*(-3*m1**2 - 2*m2**2 - m4**2 + 2*sqrtS**2 + uman)*xd13u + 2*M0Lpu**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd1u - 
     &       2*m1**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd1u - 2*M1Lpu**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd1u - 2*m2**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd1u - 
     &       2*M3Lpu**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd1u + 4*sqrtS**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd1u + 
     &       3*M0Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd1u - 
     &       3*m1**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd1u - 
     &       3*M1Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd1u - 
     &       3*m2**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd1u - 
     &       3*M3Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd1u + 
     &       M0Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd1u - 
     &       m1**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd1u - 
     &       M1Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd1u - 
     &       m2**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd1u - 
     &       M3Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd1u - 
     &       6*sqrtS**2*(3*m1**2 - m3**2 + tman)*(-3*m1**2 - 2*m2**2 - m4**2 + 2*sqrtS**2 + uman)*xd1u + 
     &       2*sqrtS**2*(-2*m1**2 - 3*m2**2 - m3**2 + 2*sqrtS**2 + uman)*(-3*m1**2 - 2*m2**2 - m4**2 + 2*sqrtS**2 + uman)*xd1u + 4*M0Lpu**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd22u - 
     &       4*m1**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd22u - 4*M1Lpu**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd22u - 4*m2**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd22u - 
     &       4*M3Lpu**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd22u + 8*sqrtS**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd22u + 
     &       2*M0Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd22u - 
     &       2*m1**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd22u - 
     &       2*M1Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd22u - 
     &       2*m2**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd22u - 
     &       2*M3Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd22u - 
     &       4*sqrtS**2*(3*m1**2 - m3**2 + tman)*(-3*m1**2 - 2*m2**2 - m4**2 + 2*sqrtS**2 + uman)*xd22u + 2*M0Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd23u - 
     &       2*m1**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd23u - 2*M1Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd23u - 
     &       2*m2**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd23u - 2*M3Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd23u + 
     &       4*sqrtS**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd23u + 4*M0Lpu**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd23u - 
     &       4*m1**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd23u - 4*M1Lpu**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd23u - 
     &       4*m2**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd23u - 4*M3Lpu**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd23u + 
     &       2*M0Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd23u - 
     &       2*m1**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd23u - 
     &       2*M1Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd23u - 
     &       2*m2**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd23u - 
     &       2*M3Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd23u - 
     &       8*sqrtS**2*(m2**2 - tman)*(-2*m1**2 - 3*m2**2 - m3**2 + 2*sqrtS**2 + uman)*xd23u + 
     &       4*sqrtS**2*(-2*m1**2 - 3*m2**2 - m3**2 + 2*sqrtS**2 + uman)*(-3*m1**2 - 2*m2**2 - m4**2 + 2*sqrtS**2 + uman)*xd23u + 4*M0Lpu**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd2u - 
     &       4*m1**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd2u - 4*M1Lpu**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd2u - 4*m2**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd2u - 
     &       4*M3Lpu**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd2u + 8*sqrtS**2*(m2**2 - tman)*(3*m1**2 - m3**2 + tman)*xd2u + 
     &       2*M0Lpu**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd2u - 2*m1**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd2u - 
     &       2*M1Lpu**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd2u - 2*m2**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd2u - 
     &       2*M3Lpu**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd2u + 3*M0Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd2u - 
     &       3*m1**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd2u - 
     &       3*M1Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd2u - 
     &       3*m2**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd2u - 
     &       3*M3Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd2u + 
     &       M0Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd2u - 
     &       m1**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd2u - 
     &       M1Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd2u - 
     &       m2**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd2u - 
     &       M3Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd2u - 
     &       4*sqrtS**2*(m2**2 - tman)*(-2*m1**2 - 3*m2**2 - m3**2 + 2*sqrtS**2 + uman)*xd2u - 6*sqrtS**2*(3*m1**2 - m3**2 + tman)*(-3*m1**2 - 2*m2**2 - m4**2 + 2*sqrtS**2 + uman)*xd2u + 
     &       2*sqrtS**2*(-2*m1**2 - 3*m2**2 - m3**2 + 2*sqrtS**2 + uman)*(-3*m1**2 - 2*m2**2 - m4**2 + 2*sqrtS**2 + uman)*xd2u + 
     &       2*M0Lpu**2*(3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd33u - 
     &       2*m1**2*(3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd33u - 
     &       2*M1Lpu**2*(3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd33u - 
     &       2*m2**2*(3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd33u - 
     &       2*M3Lpu**2*(3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd33u - 
     &       4*sqrtS**2*(3*m2**2 - m4**2 + tman)*(-2*m1**2 - 3*m2**2 - m3**2 + 2*sqrtS**2 + uman)*xd33u + M0Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd3u - 
     &       m1**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd3u - M1Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd3u - 
     &       m2**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd3u - M3Lpu**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd3u + 
     &       2*sqrtS**2*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman)*xd3u + 2*M0Lpu**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd3u - 
     &       2*m1**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd3u - 2*M1Lpu**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd3u - 
     &       2*m2**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd3u - 2*M3Lpu**2*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd3u + 
     &       M0Lpu**2*(3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd3u - m1**2*(3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd3u - 
     &       M1Lpu**2*(3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd3u - m2**2*(3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd3u - 
     &       M3Lpu**2*(3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*xd3u + 
     &       2*M0Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd3u - 
     &       2*m1**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd3u - 
     &       2*M1Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd3u - 
     &       2*m2**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd3u - 
     &       2*M3Lpu**2*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*xd3u - 
     &       4*sqrtS**2*(m2**2 - tman)*(-2*m1**2 - 3*m2**2 - m3**2 + 2*sqrtS**2 + uman)*xd3u - 2*sqrtS**2*(3*m2**2 - m4**2 + tman)*(-2*m1**2 - 3*m2**2 - m3**2 + 2*sqrtS**2 + uman)*xd3u + 
     &       4*sqrtS**2*(-2*m1**2 - 3*m2**2 - m3**2 + 2*sqrtS**2 + uman)*(-3*m1**2 - 2*m2**2 - m4**2 + 2*sqrtS**2 + uman)*xd3u + 
     &       (2*m1**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*uman)*((m2**2 - tman)*
     &           (xc01u + xc02u - xc03u + 4*xc111u - 4*xc112u + 4*xc11u - 4*xc13u + 4*xc21u + 4*xc223u - M0Lpu**2*xd0u + m1**2*xd0u + M1Lpu**2*xd0u + m2**2*xd0u + M3Lpu**2*xd0u - 
     &             2*sqrtS**2*xd0u - 4*M0Lpu**2*xd22u + 4*m1**2*xd22u + 4*M1Lpu**2*xd22u + 4*m2**2*xd22u + 4*M3Lpu**2*xd22u - 8*sqrtS**2*xd22u - 4*M0Lpu**2*xd2u + 4*m1**2*xd2u + 
     &             4*M1Lpu**2*xd2u + 4*m2**2*xd2u + 4*M3Lpu**2*xd2u - 8*sqrtS**2*xd2u) + 
     &          (3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*(xc01u + xc02u + 2*xc111u + 3*xc11u + 2*xc121u + 2*xc123u + xc12u + 3*xc21u + 2*xc223u + xc23u - M0Lpu**2*xd0u + m1**2*xd0u + 
     &             M1Lpu**2*xd0u + m2**2*xd0u + M3Lpu**2*xd0u - 2*sqrtS**2*xd0u - 2*M0Lpu**2*xd12u + 2*m1**2*xd12u + 2*M1Lpu**2*xd12u + 2*m2**2*xd12u + 2*M3Lpu**2*xd12u - 4*sqrtS**2*xd12u - 
     &             M0Lpu**2*xd1u + m1**2*xd1u + M1Lpu**2*xd1u + m2**2*xd1u + M3Lpu**2*xd1u - 2*sqrtS**2*xd1u - 2*M0Lpu**2*xd22u + 2*m1**2*xd22u + 2*M1Lpu**2*xd22u + 2*m2**2*xd22u + 
     &             2*M3Lpu**2*xd22u - 4*sqrtS**2*xd22u - 3*M0Lpu**2*xd2u + 3*m1**2*xd2u + 3*M1Lpu**2*xd2u + 3*m2**2*xd2u + 3*M3Lpu**2*xd2u - 6*sqrtS**2*xd2u) + 
     &          (3*m2**2 - m4**2 + tman)*(-2*xc121u + 2*xc122u + xc22u - xc23u - 2*M0Lpu**2*xd23u + 2*m1**2*xd23u + 2*M1Lpu**2*xd23u + 2*m2**2*xd23u + 2*M3Lpu**2*xd23u - 4*sqrtS**2*xd23u - 
     &             M0Lpu**2*xd3u + m1**2*xd3u + M1Lpu**2*xd3u + m2**2*xd3u + M3Lpu**2*xd3u - 2*sqrtS**2*xd3u))))/(32.*pi**2*(-sfmc**2 + tman))  )

!Box(t)*u + Box(u)*u
 
        Coupl_Ac = dconjg(VSfSfu(ytype,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
        Coupl_Bc = dconjg(VSfSft(xtype,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
        sfmc = MSf(j,proptypeuc,propgenuc)
        

!t
	Box_TU = Box_TU + 2*colorfactorTU*dreal(   -(Coupl_Ac*Coupl_Bc*g0t*g1t*g2t*g3t*(-4*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*tman)*
     &        (-xc001t + xc002t + xc003t - M0Lpt**2*xd00t + m1**2*xd00t + M1Lpt**2*xd00t + m2**2*xd00t + M3Lpt**2*xd00t - 2*sqrtS**2*xd00t) + 
     &       (-m2**2 + uman)*(-((2*m1**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*tman)*
     &             (xc01t + xc02t - xc03t + 4*xc111t - 4*xc112t + 4*xc11t - 4*xc13t + 4*xc21t + 4*xc223t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 
     &               2*sqrtS**2*xd0t - 4*M0Lpt**2*xd22t + 4*m1**2*xd22t + 4*M1Lpt**2*xd22t + 4*m2**2*xd22t + 4*M3Lpt**2*xd22t - 8*sqrtS**2*xd22t - 4*M0Lpt**2*xd2t + 4*m1**2*xd2t + 
     &               4*M1Lpt**2*xd2t + 4*m2**2*xd2t + 4*M3Lpt**2*xd2t - 8*sqrtS**2*xd2t)) + 
     &          (3*m1**2 - m4**2 + uman)*(xc01t + xc02t + xc03t + 4*xc111t + 4*xc11t + 4*xc121t + 4*xc123t + 2*xc12t + 2*xc13t + 4*xc21t + 4*xc223t + 2*xc23t - M0Lpt**2*xd0t + m1**2*xd0t + 
     &             M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 4*M0Lpt**2*xd12t + 4*m1**2*xd12t + 4*M1Lpt**2*xd12t + 4*m2**2*xd12t + 4*M3Lpt**2*xd12t - 8*sqrtS**2*xd12t - 
     &             2*M0Lpt**2*xd1t + 2*m1**2*xd1t + 2*M1Lpt**2*xd1t + 2*m2**2*xd1t + 2*M3Lpt**2*xd1t - 4*sqrtS**2*xd1t - 4*M0Lpt**2*xd22t + 4*m1**2*xd22t + 4*M1Lpt**2*xd22t + 4*m2**2*xd22t + 
     &             4*M3Lpt**2*xd22t - 8*sqrtS**2*xd22t - 4*M0Lpt**2*xd2t + 4*m1**2*xd2t + 4*M1Lpt**2*xd2t + 4*m2**2*xd2t + 4*M3Lpt**2*xd2t - 8*sqrtS**2*xd2t) + 
     &          (2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(xc01t + xc02t - xc03t + 2*xc11t - 4*xc121t + 4*xc122t - 2*xc13t + 2*xc21t + 2*xc22t - 2*xc23t - M0Lpt**2*xd0t + m1**2*xd0t + 
     &             M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 4*M0Lpt**2*xd23t + 4*m1**2*xd23t + 4*M1Lpt**2*xd23t + 4*m2**2*xd23t + 4*M3Lpt**2*xd23t - 8*sqrtS**2*xd23t - 
     &             2*M0Lpt**2*xd2t + 2*m1**2*xd2t + 2*M1Lpt**2*xd2t + 2*m2**2*xd2t + 2*M3Lpt**2*xd2t - 4*sqrtS**2*xd2t - 2*M0Lpt**2*xd3t + 2*m1**2*xd3t + 2*M1Lpt**2*xd3t + 2*m2**2*xd3t + 
     &             2*M3Lpt**2*xd3t - 4*sqrtS**2*xd3t)) + (3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*
     &        (-((3*m1**2 - m4**2 + uman)*(xc01t + xc02t + 2*xc111t + 2*xc113t + 3*xc11t + 4*xc123t + 3*xc12t + 3*xc21t - 2*xc221t + 2*xc223t - xc23t - M0Lpt**2*xd0t + m1**2*xd0t + 
     &               M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 2*M0Lpt**2*xd11t + 2*m1**2*xd11t + 2*M1Lpt**2*xd11t + 2*m2**2*xd11t + 2*M3Lpt**2*xd11t - 
     &               4*sqrtS**2*xd11t - 4*M0Lpt**2*xd12t + 4*m1**2*xd12t + 4*M1Lpt**2*xd12t + 4*m2**2*xd12t + 4*M3Lpt**2*xd12t - 8*sqrtS**2*xd12t - 3*M0Lpt**2*xd1t + 3*m1**2*xd1t + 
     &               3*M1Lpt**2*xd1t + 3*m2**2*xd1t + 3*M3Lpt**2*xd1t - 6*sqrtS**2*xd1t - 2*M0Lpt**2*xd22t + 2*m1**2*xd22t + 2*M1Lpt**2*xd22t + 2*m2**2*xd22t + 2*M3Lpt**2*xd22t - 
     &               4*sqrtS**2*xd22t - 3*M0Lpt**2*xd2t + 3*m1**2*xd2t + 3*M1Lpt**2*xd2t + 3*m2**2*xd2t + 3*M3Lpt**2*xd2t - 6*sqrtS**2*xd2t)) + 
     &          (2*m1**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*tman)*(xc01t + xc02t + 2*xc111t + 3*xc11t + 2*xc121t + 2*xc123t + xc12t + 3*xc21t + 2*xc223t + xc23t - M0Lpt**2*xd0t + m1**2*xd0t + 
     &             M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 2*M0Lpt**2*xd12t + 2*m1**2*xd12t + 2*M1Lpt**2*xd12t + 2*m2**2*xd12t + 2*M3Lpt**2*xd12t - 4*sqrtS**2*xd12t - 
     &             M0Lpt**2*xd1t + m1**2*xd1t + M1Lpt**2*xd1t + m2**2*xd1t + M3Lpt**2*xd1t - 2*sqrtS**2*xd1t - 2*M0Lpt**2*xd22t + 2*m1**2*xd22t + 2*M1Lpt**2*xd22t + 2*m2**2*xd22t + 
     &             2*M3Lpt**2*xd22t - 4*sqrtS**2*xd22t - 3*M0Lpt**2*xd2t + 3*m1**2*xd2t + 3*M1Lpt**2*xd2t + 3*m2**2*xd2t + 3*M3Lpt**2*xd2t - 6*sqrtS**2*xd2t) - 
     &          (2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(xc01t + xc02t + xc11t + 2*xc122t + xc12t + xc21t + 2*xc221t + 2*xc22t + xc23t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + 
     &             m2**2*xd0t + M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 2*M0Lpt**2*xd13t + 2*m1**2*xd13t + 2*M1Lpt**2*xd13t + 2*m2**2*xd13t + 2*M3Lpt**2*xd13t - 4*sqrtS**2*xd13t - M0Lpt**2*xd1t + 
     &             m1**2*xd1t + M1Lpt**2*xd1t + m2**2*xd1t + M3Lpt**2*xd1t - 2*sqrtS**2*xd1t - 2*M0Lpt**2*xd23t + 2*m1**2*xd23t + 2*M1Lpt**2*xd23t + 2*m2**2*xd23t + 2*M3Lpt**2*xd23t - 
     &             4*sqrtS**2*xd23t - M0Lpt**2*xd2t + m1**2*xd2t + M1Lpt**2*xd2t + m2**2*xd2t + M3Lpt**2*xd2t - 2*sqrtS**2*xd2t - 2*M0Lpt**2*xd3t + 2*m1**2*xd3t + 2*M1Lpt**2*xd3t + 
     &             2*m2**2*xd3t + 2*M3Lpt**2*xd3t - 4*sqrtS**2*xd3t)) + (-3*m2**2 + m3**2 - uman)*
     &        (-((2*m1**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*tman)*(-2*xc121t + 2*xc122t + xc22t - xc23t - 2*M0Lpt**2*xd23t + 2*m1**2*xd23t + 2*M1Lpt**2*xd23t + 2*m2**2*xd23t + 
     &               2*M3Lpt**2*xd23t - 4*sqrtS**2*xd23t - M0Lpt**2*xd3t + m1**2*xd3t + M1Lpt**2*xd3t + m2**2*xd3t + M3Lpt**2*xd3t - 2*sqrtS**2*xd3t)) + 
     &          (3*m1**2 - m4**2 + uman)*(2*xc122t + 2*xc221t + xc22t + xc23t - 2*M0Lpt**2*xd13t + 2*m1**2*xd13t + 2*M1Lpt**2*xd13t + 2*m2**2*xd13t + 2*M3Lpt**2*xd13t - 4*sqrtS**2*xd13t - 
     &             2*M0Lpt**2*xd23t + 2*m1**2*xd23t + 2*M1Lpt**2*xd23t + 2*m2**2*xd23t + 2*M3Lpt**2*xd23t - 4*sqrtS**2*xd23t - M0Lpt**2*xd3t + m1**2*xd3t + M1Lpt**2*xd3t + m2**2*xd3t + 
     &             M3Lpt**2*xd3t - 2*sqrtS**2*xd3t) + (2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*
     &           (-2*xc221t + 2*xc222t + xc22t - xc23t - 2*M0Lpt**2*xd33t + 2*m1**2*xd33t + 2*M1Lpt**2*xd33t + 2*m2**2*xd33t + 2*M3Lpt**2*xd33t - 4*sqrtS**2*xd33t - M0Lpt**2*xd3t + 
     &             m1**2*xd3t + M1Lpt**2*xd3t + m2**2*xd3t + M3Lpt**2*xd3t - 2*sqrtS**2*xd3t))))/(32.*pi**2*(-sfmc**2 + uman))  )

!u
	Box_UU = Box_UU + 2*colorfactorTU*dreal(  -(Coupl_Ac*Coupl_Bc*g0u*g1u*g2u*g3u*((3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xc01u + (m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xc01u - 
     &       (3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xc01u - (3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xc01u + 
     &       (3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xc02u + (m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xc02u - 
     &       (3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xc02u - (3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xc02u + 
     &       (3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xc03u + (-m1**2 + uman)*(3*m2**2 - m3**2 + uman)*xc03u + 
     &       4*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xc111u - 2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xc111u - 
     &       2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xc113u + 4*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xc11u + 
     &       2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xc11u - 3*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xc11u - 
     &       (3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xc11u + 4*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xc121u + 
     &       4*(-m1**2 + uman)*(3*m2**2 - m3**2 + uman)*xc121u - 2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xc122u + 
     &       4*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xc122u - 2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xc122u + 
     &       4*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xc123u - 4*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xc123u + 
     &       2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xc12u - 3*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xc12u - 
     &       (3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xc12u + 2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xc13u + 
     &       2*(-m1**2 + uman)*(3*m2**2 - m3**2 + uman)*xc13u + 4*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xc21u + 2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xc21u - 
     &       3*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xc21u - (3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xc21u - 
     &       2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xc221u + 
     &       2*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman)*xc221u + 2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xc221u - 
     &       2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xc221u - 2*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman)*xc222u + 
     &       4*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xc223u - 2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xc223u - 
     &       (3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xc22u - 
     &       (2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman)*xc22u + 2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xc22u - 
     &       2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xc22u - (3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xc23u + 
     &       2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xc23u + (2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman)*xc23u + 
     &       2*(-m1**2 + uman)*(3*m2**2 - m3**2 + uman)*xc23u + (3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xc23u - 
     &       (3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xc23u - 4*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*tman)*
     &        (-xc001u + xc002u + xc003u - M0Lpu**2*xd00u + m1**2*xd00u + M1Lpu**2*xd00u + m2**2*xd00u + M3Lpu**2*xd00u - 2*sqrtS**2*xd00u) - 
     &       M0Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd0u + m1**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd0u + 
     &       M1Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd0u + m2**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd0u + 
     &       M3Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd0u + 2*sqrtS**2*(-3*m1**2 - 2*m2**2 - m3**2 + 2*sqrtS**2 + tman)*(m1**2 - uman)*xd0u + 
     &       m1**2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xd0u + M1Lpu**2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xd0u + m2**2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xd0u + 
     &       M3Lpu**2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xd0u + M0Lpu**2*(-m1**2 + uman)*(3*m2**2 - m3**2 + uman)*xd0u + 2*sqrtS**2*(-m1**2 + uman)*(3*m2**2 - m3**2 + uman)*xd0u + 
     &       M0Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd0u - m1**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd0u - 
     &       M1Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd0u - m2**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd0u - 
     &       M3Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd0u - 
     &       2*sqrtS**2*(-3*m1**2 - 2*m2**2 - m3**2 + 2*sqrtS**2 + tman)*(3*m1**2 - m4**2 + uman)*xd0u + M0Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd0u - 
     &       m1**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd0u - M1Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd0u - 
     &       m2**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd0u - M3Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd0u + 
     &       2*sqrtS**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd0u + 2*M0Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd11u - 
     &       2*m1**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd11u - 
     &       2*M1Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd11u - 
     &       2*m2**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd11u - 
     &       2*M3Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd11u - 
     &       4*sqrtS**2*(-3*m1**2 - 2*m2**2 - m3**2 + 2*sqrtS**2 + tman)*(3*m1**2 - m4**2 + uman)*xd11u - 4*M0Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd12u + 
     &       4*m1**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd12u + 4*M1Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd12u + 
     &       4*m2**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd12u + 4*M3Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd12u + 
     &       8*sqrtS**2*(-3*m1**2 - 2*m2**2 - m3**2 + 2*sqrtS**2 + tman)*(m1**2 - uman)*xd12u + 4*M0Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd12u - 
     &       4*m1**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd12u - 
     &       4*M1Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd12u - 
     &       4*m2**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd12u - 
     &       4*M3Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd12u - 
     &       8*sqrtS**2*(-3*m1**2 - 2*m2**2 - m3**2 + 2*sqrtS**2 + tman)*(3*m1**2 - m4**2 + uman)*xd12u + 
     &       2*M0Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xd13u - 
     &       2*m1**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xd13u - 
     &       2*M1Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xd13u - 
     &       2*m2**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xd13u - 
     &       2*M3Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xd13u + 
     &       4*sqrtS**2*(-3*m1**2 - 2*m2**2 - m3**2 + 2*sqrtS**2 + tman)*(-2*m1**2 - 3*m2**2 - m4**2 + 2*sqrtS**2 + tman)*xd13u + 
     &       2*M0Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd13u - 2*m1**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd13u - 
     &       2*M1Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd13u - 2*m2**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd13u - 
     &       2*M3Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd13u + 4*sqrtS**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd13u - 
     &       2*M0Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd1u + 2*m1**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd1u + 
     &       2*M1Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd1u + 2*m2**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd1u + 
     &       2*M3Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd1u + 4*sqrtS**2*(-3*m1**2 - 2*m2**2 - m3**2 + 2*sqrtS**2 + tman)*(m1**2 - uman)*xd1u + 
     &       3*M0Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd1u - 
     &       3*m1**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd1u - 
     &       3*M1Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd1u - 
     &       3*m2**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd1u - 
     &       3*M3Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd1u - 
     &       6*sqrtS**2*(-3*m1**2 - 2*m2**2 - m3**2 + 2*sqrtS**2 + tman)*(3*m1**2 - m4**2 + uman)*xd1u + M0Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd1u - 
     &       m1**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd1u - M1Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd1u - 
     &       m2**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd1u - M3Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd1u + 
     &       2*sqrtS**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd1u - 4*M0Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd22u + 
     &       4*m1**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd22u + 4*M1Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd22u + 
     &       4*m2**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd22u + 4*M3Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd22u + 
     &       8*sqrtS**2*(-3*m1**2 - 2*m2**2 - m3**2 + 2*sqrtS**2 + tman)*(m1**2 - uman)*xd22u + 2*M0Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd22u - 
     &       2*m1**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd22u - 
     &       2*M1Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd22u - 
     &       2*m2**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd22u - 
     &       2*M3Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd22u - 
     &       4*sqrtS**2*(-3*m1**2 - 2*m2**2 - m3**2 + 2*sqrtS**2 + tman)*(3*m1**2 - m4**2 + uman)*xd22u + 
     &       2*M0Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xd23u - 
     &       2*m1**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xd23u - 
     &       2*M1Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xd23u - 
     &       2*m2**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xd23u - 
     &       2*M3Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xd23u + 
     &       4*sqrtS**2*(-3*m1**2 - 2*m2**2 - m3**2 + 2*sqrtS**2 + tman)*(-2*m1**2 - 3*m2**2 - m4**2 + 2*sqrtS**2 + tman)*xd23u + 4*m1**2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xd23u + 
     &       4*M1Lpu**2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xd23u + 4*m2**2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xd23u + 4*M3Lpu**2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xd23u + 
     &       4*M0Lpu**2*(-m1**2 + uman)*(3*m2**2 - m3**2 + uman)*xd23u + 8*sqrtS**2*(-m1**2 + uman)*(3*m2**2 - m3**2 + uman)*xd23u + 
     &       2*M0Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd23u - 2*m1**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd23u - 
     &       2*M1Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd23u - 2*m2**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd23u - 
     &       2*M3Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd23u + 4*sqrtS**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd23u - 
     &       4*M0Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd2u + 4*m1**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd2u + 
     &       4*M1Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd2u + 4*m2**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd2u + 
     &       4*M3Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman)*xd2u + 8*sqrtS**2*(-3*m1**2 - 2*m2**2 - m3**2 + 2*sqrtS**2 + tman)*(m1**2 - uman)*xd2u + 
     &       2*m1**2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xd2u + 2*M1Lpu**2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xd2u + 2*m2**2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xd2u + 
     &       2*M3Lpu**2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xd2u + 2*M0Lpu**2*(-m1**2 + uman)*(3*m2**2 - m3**2 + uman)*xd2u + 4*sqrtS**2*(-m1**2 + uman)*(3*m2**2 - m3**2 + uman)*xd2u + 
     &       3*M0Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd2u - 
     &       3*m1**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd2u - 
     &       3*M1Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd2u - 
     &       3*m2**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd2u - 
     &       3*M3Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman)*xd2u - 
     &       6*sqrtS**2*(-3*m1**2 - 2*m2**2 - m3**2 + 2*sqrtS**2 + tman)*(3*m1**2 - m4**2 + uman)*xd2u + M0Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd2u - 
     &       m1**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd2u - M1Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd2u - 
     &       m2**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd2u - M3Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd2u + 
     &       2*sqrtS**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd2u + 2*M0Lpu**2*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman)*xd33u - 
     &       2*m1**2*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman)*xd33u - 
     &       2*M1Lpu**2*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman)*xd33u - 
     &       2*m2**2*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman)*xd33u - 
     &       2*M3Lpu**2*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman)*xd33u - 
     &       4*sqrtS**2*(-2*m1**2 - 3*m2**2 - m4**2 + 2*sqrtS**2 + tman)*(3*m2**2 - m3**2 + uman)*xd33u + 
     &       M0Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xd3u - 
     &       m1**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xd3u - 
     &       M1Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xd3u - 
     &       m2**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xd3u - 
     &       M3Lpu**2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*xd3u + 
     &       2*sqrtS**2*(-3*m1**2 - 2*m2**2 - m3**2 + 2*sqrtS**2 + tman)*(-2*m1**2 - 3*m2**2 - m4**2 + 2*sqrtS**2 + tman)*xd3u + 
     &       M0Lpu**2*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman)*xd3u - m1**2*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman)*xd3u - 
     &       M1Lpu**2*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman)*xd3u - m2**2*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman)*xd3u - 
     &       M3Lpu**2*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman)*xd3u - 
     &       2*sqrtS**2*(-2*m1**2 - 3*m2**2 - m4**2 + 2*sqrtS**2 + tman)*(3*m2**2 - m3**2 + uman)*xd3u + 2*m1**2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xd3u + 
     &       2*M1Lpu**2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xd3u + 2*m2**2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xd3u + 2*M3Lpu**2*(m1**2 - uman)*(3*m2**2 - m3**2 + uman)*xd3u + 
     &       2*M0Lpu**2*(-m1**2 + uman)*(3*m2**2 - m3**2 + uman)*xd3u + 4*sqrtS**2*(-m1**2 + uman)*(3*m2**2 - m3**2 + uman)*xd3u + 
     &       2*M0Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd3u - 2*m1**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd3u - 
     &       2*M1Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd3u - 2*m2**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd3u - 
     &       2*M3Lpu**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd3u + 4*sqrtS**2*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)*xd3u + 
     &       (-2*m2**2 - m3**2 - 3*m4**2 + sqrtS**2 + 2*tman)*(-((m1**2 - uman)*
     &             (xc01u + xc02u - xc03u + 4*xc111u - 4*xc112u + 4*xc11u - 4*xc13u + 4*xc21u + 4*xc223u - M0Lpu**2*xd0u + m1**2*xd0u + M1Lpu**2*xd0u + m2**2*xd0u + M3Lpu**2*xd0u - 
     &               2*sqrtS**2*xd0u - 4*M0Lpu**2*xd22u + 4*m1**2*xd22u + 4*M1Lpu**2*xd22u + 4*m2**2*xd22u + 4*M3Lpu**2*xd22u - 8*sqrtS**2*xd22u - 4*M0Lpu**2*xd2u + 4*m1**2*xd2u + 
     &               4*M1Lpu**2*xd2u + 4*m2**2*xd2u + 4*M3Lpu**2*xd2u - 8*sqrtS**2*xd2u)) + 
     &          (3*m1**2 - m4**2 + uman)*(xc01u + xc02u + 2*xc111u + 3*xc11u + 2*xc121u + 2*xc123u + xc12u + 3*xc21u + 2*xc223u + xc23u - M0Lpu**2*xd0u + m1**2*xd0u + M1Lpu**2*xd0u + 
     &             m2**2*xd0u + M3Lpu**2*xd0u - 2*sqrtS**2*xd0u - 2*M0Lpu**2*xd12u + 2*m1**2*xd12u + 2*M1Lpu**2*xd12u + 2*m2**2*xd12u + 2*M3Lpu**2*xd12u - 4*sqrtS**2*xd12u - M0Lpu**2*xd1u + 
     &             m1**2*xd1u + M1Lpu**2*xd1u + m2**2*xd1u + M3Lpu**2*xd1u - 2*sqrtS**2*xd1u - 2*M0Lpu**2*xd22u + 2*m1**2*xd22u + 2*M1Lpu**2*xd22u + 2*m2**2*xd22u + 2*M3Lpu**2*xd22u - 
     &             4*sqrtS**2*xd22u - 3*M0Lpu**2*xd2u + 3*m1**2*xd2u + 3*M1Lpu**2*xd2u + 3*m2**2*xd2u + 3*M3Lpu**2*xd2u - 6*sqrtS**2*xd2u) + 
     &          (2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(-2*xc121u + 2*xc122u + xc22u - xc23u - 2*M0Lpu**2*xd23u + 2*m1**2*xd23u + 2*M1Lpu**2*xd23u + 2*m2**2*xd23u + 
     &             2*M3Lpu**2*xd23u - 4*sqrtS**2*xd23u - M0Lpu**2*xd3u + m1**2*xd3u + M1Lpu**2*xd3u + m2**2*xd3u + M3Lpu**2*xd3u - 2*sqrtS**2*xd3u))))/(32.*pi**2*(-sfmc**2 + uman))  )

        enddo

!================================================================================================================================== 
!Box(t)*q + Box(u)*q

        Coupl_Ac = dconjg(VVSfSf(xtype,ytype,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))

!t
	Box_TQ = Box_TQ + 2*colorfactorQ*dreal(  (Coupl_Ac*g0t*g1t*g2t*g3t*(32*(xc001t - xc002t - xc003t + M0Lpt**2*xd00t - m1**2*xd00t - M1Lpt**2*xd00t - m2**2*xd00t - M3Lpt**2*xd00t + 2*sqrtS**2*xd00t) - 
     &      2*m3**2*(xc01t + xc02t - xc03t + 4*xc111t - 4*xc112t + 4*xc11t - 4*xc13t + 4*xc21t + 4*xc223t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 
     &         2*sqrtS**2*xd0t - 4*M0Lpt**2*xd22t + 4*m1**2*xd22t + 4*M1Lpt**2*xd22t + 4*m2**2*xd22t + 4*M3Lpt**2*xd22t - 8*sqrtS**2*xd22t - 4*M0Lpt**2*xd2t + 4*m1**2*xd2t + 
     &         4*M1Lpt**2*xd2t + 4*m2**2*xd2t + 4*M3Lpt**2*xd2t - 8*sqrtS**2*xd2t) + 
     &      (m1**2 + m3**2 - tman)*(xc01t + xc02t + xc03t + 4*xc111t + 4*xc11t + 4*xc121t + 4*xc123t + 2*xc12t + 2*xc13t + 4*xc21t + 4*xc223t + 2*xc23t - M0Lpt**2*xd0t + m1**2*xd0t + 
     &         M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 4*M0Lpt**2*xd12t + 4*m1**2*xd12t + 4*M1Lpt**2*xd12t + 4*m2**2*xd12t + 4*M3Lpt**2*xd12t - 8*sqrtS**2*xd12t - 
     &         2*M0Lpt**2*xd1t + 2*m1**2*xd1t + 2*M1Lpt**2*xd1t + 2*m2**2*xd1t + 2*M3Lpt**2*xd1t - 4*sqrtS**2*xd1t - 4*M0Lpt**2*xd22t + 4*m1**2*xd22t + 4*M1Lpt**2*xd22t + 4*m2**2*xd22t + 
     &         4*M3Lpt**2*xd22t - 8*sqrtS**2*xd22t - 4*M0Lpt**2*xd2t + 4*m1**2*xd2t + 4*M1Lpt**2*xd2t + 4*m2**2*xd2t + 4*M3Lpt**2*xd2t - 8*sqrtS**2*xd2t) - 
     &      (m2**2 + m3**2 - uman)*(xc01t + xc02t - xc03t + 2*xc11t - 4*xc121t + 4*xc122t - 2*xc13t + 2*xc21t + 2*xc22t - 2*xc23t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + 
     &         M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 4*M0Lpt**2*xd23t + 4*m1**2*xd23t + 4*M1Lpt**2*xd23t + 4*m2**2*xd23t + 4*M3Lpt**2*xd23t - 8*sqrtS**2*xd23t - 2*M0Lpt**2*xd2t + 2*m1**2*xd2t + 
     &         2*M1Lpt**2*xd2t + 2*m2**2*xd2t + 2*M3Lpt**2*xd2t - 4*sqrtS**2*xd2t - 2*M0Lpt**2*xd3t + 2*m1**2*xd3t + 2*M1Lpt**2*xd3t + 2*m2**2*xd3t + 2*M3Lpt**2*xd3t - 4*sqrtS**2*xd3t) + 
     &      2*(-2*m1**2*(xc01t + xc02t + 2*xc111t + 2*xc113t + 3*xc11t + 4*xc123t + 3*xc12t + 3*xc21t - 2*xc221t + 2*xc223t - xc23t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + 
     &            m2**2*xd0t + M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 2*M0Lpt**2*xd11t + 2*m1**2*xd11t + 2*M1Lpt**2*xd11t + 2*m2**2*xd11t + 2*M3Lpt**2*xd11t - 4*sqrtS**2*xd11t - 
     &            4*M0Lpt**2*xd12t + 4*m1**2*xd12t + 4*M1Lpt**2*xd12t + 4*m2**2*xd12t + 4*M3Lpt**2*xd12t - 8*sqrtS**2*xd12t - 3*M0Lpt**2*xd1t + 3*m1**2*xd1t + 3*M1Lpt**2*xd1t + 
     &            3*m2**2*xd1t + 3*M3Lpt**2*xd1t - 6*sqrtS**2*xd1t - 2*M0Lpt**2*xd22t + 2*m1**2*xd22t + 2*M1Lpt**2*xd22t + 2*m2**2*xd22t + 2*M3Lpt**2*xd22t - 4*sqrtS**2*xd22t - 
     &            3*M0Lpt**2*xd2t + 3*m1**2*xd2t + 3*M1Lpt**2*xd2t + 3*m2**2*xd2t + 3*M3Lpt**2*xd2t - 6*sqrtS**2*xd2t) + 
     &         (m1**2 + m3**2 - tman)*(xc01t + xc02t + 2*xc111t + 3*xc11t + 2*xc121t + 2*xc123t + xc12t + 3*xc21t + 2*xc223t + xc23t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + 
     &            m2**2*xd0t + M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 2*M0Lpt**2*xd12t + 2*m1**2*xd12t + 2*M1Lpt**2*xd12t + 2*m2**2*xd12t + 2*M3Lpt**2*xd12t - 4*sqrtS**2*xd12t - M0Lpt**2*xd1t + 
     &            m1**2*xd1t + M1Lpt**2*xd1t + m2**2*xd1t + M3Lpt**2*xd1t - 2*sqrtS**2*xd1t - 2*M0Lpt**2*xd22t + 2*m1**2*xd22t + 2*M1Lpt**2*xd22t + 2*m2**2*xd22t + 2*M3Lpt**2*xd22t - 
     &            4*sqrtS**2*xd22t - 3*M0Lpt**2*xd2t + 3*m1**2*xd2t + 3*M1Lpt**2*xd2t + 3*m2**2*xd2t + 3*M3Lpt**2*xd2t - 6*sqrtS**2*xd2t) - 
     &         (m1**2 + m2**2 - sqrtS**2)*(xc01t + xc02t + xc11t + 2*xc122t + xc12t + xc21t + 2*xc221t + 2*xc22t + xc23t - M0Lpt**2*xd0t + m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + 
     &            M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 2*M0Lpt**2*xd13t + 2*m1**2*xd13t + 2*M1Lpt**2*xd13t + 2*m2**2*xd13t + 2*M3Lpt**2*xd13t - 4*sqrtS**2*xd13t - M0Lpt**2*xd1t + m1**2*xd1t + 
     &            M1Lpt**2*xd1t + m2**2*xd1t + M3Lpt**2*xd1t - 2*sqrtS**2*xd1t - 2*M0Lpt**2*xd23t + 2*m1**2*xd23t + 2*M1Lpt**2*xd23t + 2*m2**2*xd23t + 2*M3Lpt**2*xd23t - 4*sqrtS**2*xd23t - 
     &            M0Lpt**2*xd2t + m1**2*xd2t + M1Lpt**2*xd2t + m2**2*xd2t + M3Lpt**2*xd2t - 2*sqrtS**2*xd2t - 2*M0Lpt**2*xd3t + 2*m1**2*xd3t + 2*M1Lpt**2*xd3t + 2*m2**2*xd3t + 
     &            2*M3Lpt**2*xd3t - 4*sqrtS**2*xd3t) - (m2**2 + m3**2 - uman)*
     &          (-2*xc121t + 2*xc122t + xc22t - xc23t - 2*M0Lpt**2*xd23t + 2*m1**2*xd23t + 2*M1Lpt**2*xd23t + 2*m2**2*xd23t + 2*M3Lpt**2*xd23t - 4*sqrtS**2*xd23t - M0Lpt**2*xd3t + 
     &            m1**2*xd3t + M1Lpt**2*xd3t + m2**2*xd3t + M3Lpt**2*xd3t - 2*sqrtS**2*xd3t) - 
     &         (m1**2 + m2**2 - sqrtS**2)*(2*xc122t + 2*xc221t + xc22t + xc23t - 2*M0Lpt**2*xd13t + 2*m1**2*xd13t + 2*M1Lpt**2*xd13t + 2*m2**2*xd13t + 2*M3Lpt**2*xd13t - 4*sqrtS**2*xd13t - 
     &            2*M0Lpt**2*xd23t + 2*m1**2*xd23t + 2*M1Lpt**2*xd23t + 2*m2**2*xd23t + 2*M3Lpt**2*xd23t - 4*sqrtS**2*xd23t - M0Lpt**2*xd3t + m1**2*xd3t + M1Lpt**2*xd3t + m2**2*xd3t + 
     &            M3Lpt**2*xd3t - 2*sqrtS**2*xd3t) - 2*m2**2*(-2*xc221t + 2*xc222t + xc22t - xc23t - 2*M0Lpt**2*xd33t + 2*m1**2*xd33t + 2*M1Lpt**2*xd33t + 2*m2**2*xd33t + 2*M3Lpt**2*xd33t - 
     &            4*sqrtS**2*xd33t - M0Lpt**2*xd3t + m1**2*xd3t + M1Lpt**2*xd3t + m2**2*xd3t + M3Lpt**2*xd3t - 2*sqrtS**2*xd3t))))/(32.*pi**2)  )

!u
	Box_UQ = Box_UQ + 2*colorfactorQ*dreal(  -(Coupl_Ac*g0u*g1u*g2u*g3u*(-32*xc001u + 32*xc002u + 32*xc003u + 3*m1**2*xc01u + 3*m2**2*xc01u - 2*sqrtS**2*xc01u - tman*xc01u + 3*uman*xc01u + 3*m1**2*xc02u + 3*m2**2*xc02u - 
     &       2*sqrtS**2*xc02u - tman*xc02u + 3*uman*xc02u - m1**2*xc03u - m2**2*xc03u - 4*m4**2*xc03u + tman*xc03u + uman*xc03u + 8*uman*xc111u - 8*m4**2*xc112u + 8*m1**2*xc113u + 
     &       4*m1**2*xc11u + 4*m2**2*xc11u - 2*sqrtS**2*xc11u - 2*tman*xc11u + 10*uman*xc11u - 8*m1**2*xc121u - 8*m2**2*xc121u - 16*m4**2*xc121u + 8*tman*xc121u + 8*uman*xc121u + 
     &       8*m1**2*xc122u + 16*m2**2*xc122u + 8*m4**2*xc122u - 8*sqrtS**2*xc122u - 8*tman*xc122u + 8*m1**2*xc123u - 8*m4**2*xc123u + 8*uman*xc123u + 10*m1**2*xc12u + 2*m2**2*xc12u - 
     &       4*m4**2*xc12u - 2*sqrtS**2*xc12u + 4*uman*xc12u - 2*m1**2*xc13u - 2*m2**2*xc13u - 12*m4**2*xc13u + 2*tman*xc13u + 2*uman*xc13u + 4*m1**2*xc21u + 4*m2**2*xc21u - 
     &       2*sqrtS**2*xc21u - 2*tman*xc21u + 10*uman*xc21u - 8*sqrtS**2*xc221u + 8*m2**2*xc222u + 8*uman*xc223u + 6*m1**2*xc22u + 14*m2**2*xc22u + 4*m4**2*xc22u - 6*sqrtS**2*xc22u - 
     &       4*tman*xc22u - 4*m1**2*xc23u - 4*m2**2*xc23u - 8*m4**2*xc23u - 4*sqrtS**2*xc23u + 4*tman*xc23u + 4*uman*xc23u - 32*M0Lpu**2*xd00u + 32*m1**2*xd00u + 32*M1Lpu**2*xd00u + 
     &       32*m2**2*xd00u + 32*M3Lpu**2*xd00u - 64*sqrtS**2*xd00u - 3*M0Lpu**2*m1**2*xd0u + 3*m1**4*xd0u + 3*m1**2*M1Lpu**2*xd0u - 3*M0Lpu**2*m2**2*xd0u + 6*m1**2*m2**2*xd0u + 
     &       3*M1Lpu**2*m2**2*xd0u + 3*m2**4*xd0u + 3*m1**2*M3Lpu**2*xd0u + 3*m2**2*M3Lpu**2*xd0u + 2*M0Lpu**2*sqrtS**2*xd0u - 8*m1**2*sqrtS**2*xd0u - 2*M1Lpu**2*sqrtS**2*xd0u - 
     &       8*m2**2*sqrtS**2*xd0u - 2*M3Lpu**2*sqrtS**2*xd0u + 4*sqrtS**4*xd0u + M0Lpu**2*tman*xd0u - m1**2*tman*xd0u - M1Lpu**2*tman*xd0u - m2**2*tman*xd0u - M3Lpu**2*tman*xd0u + 
     &       2*sqrtS**2*tman*xd0u - 3*M0Lpu**2*uman*xd0u + 3*m1**2*uman*xd0u + 3*M1Lpu**2*uman*xd0u + 3*m2**2*uman*xd0u + 3*M3Lpu**2*uman*xd0u - 6*sqrtS**2*uman*xd0u - 
     &       8*M0Lpu**2*m1**2*xd11u + 8*m1**4*xd11u + 8*m1**2*M1Lpu**2*xd11u + 8*m1**2*m2**2*xd11u + 8*m1**2*M3Lpu**2*xd11u - 16*m1**2*sqrtS**2*xd11u - 8*M0Lpu**2*m1**2*xd12u + 
     &       8*m1**4*xd12u + 8*m1**2*M1Lpu**2*xd12u + 8*m1**2*m2**2*xd12u + 8*m1**2*M3Lpu**2*xd12u + 8*M0Lpu**2*m4**2*xd12u - 8*m1**2*m4**2*xd12u - 8*M1Lpu**2*m4**2*xd12u - 
     &       8*m2**2*m4**2*xd12u - 8*M3Lpu**2*m4**2*xd12u - 16*m1**2*sqrtS**2*xd12u + 16*m4**2*sqrtS**2*xd12u - 8*M0Lpu**2*uman*xd12u + 8*m1**2*uman*xd12u + 8*M1Lpu**2*uman*xd12u + 
     &       8*m2**2*uman*xd12u + 8*M3Lpu**2*uman*xd12u - 16*sqrtS**2*uman*xd12u - 8*M0Lpu**2*m1**2*xd13u + 8*m1**4*xd13u + 8*m1**2*M1Lpu**2*xd13u - 8*M0Lpu**2*m2**2*xd13u + 
     &       16*m1**2*m2**2*xd13u + 8*M1Lpu**2*m2**2*xd13u + 8*m2**4*xd13u + 8*m1**2*M3Lpu**2*xd13u + 8*m2**2*M3Lpu**2*xd13u + 8*M0Lpu**2*sqrtS**2*xd13u - 24*m1**2*sqrtS**2*xd13u - 
     &       8*M1Lpu**2*sqrtS**2*xd13u - 24*m2**2*sqrtS**2*xd13u - 8*M3Lpu**2*sqrtS**2*xd13u + 16*sqrtS**4*xd13u - 10*M0Lpu**2*m1**2*xd1u + 10*m1**4*xd1u + 10*m1**2*M1Lpu**2*xd1u - 
     &       2*M0Lpu**2*m2**2*xd1u + 12*m1**2*m2**2*xd1u + 2*M1Lpu**2*m2**2*xd1u + 2*m2**4*xd1u + 10*m1**2*M3Lpu**2*xd1u + 2*m2**2*M3Lpu**2*xd1u + 4*M0Lpu**2*m4**2*xd1u - 
     &       4*m1**2*m4**2*xd1u - 4*M1Lpu**2*m4**2*xd1u - 4*m2**2*m4**2*xd1u - 4*M3Lpu**2*m4**2*xd1u + 2*M0Lpu**2*sqrtS**2*xd1u - 22*m1**2*sqrtS**2*xd1u - 2*M1Lpu**2*sqrtS**2*xd1u - 
     &       6*m2**2*sqrtS**2*xd1u - 2*M3Lpu**2*sqrtS**2*xd1u + 8*m4**2*sqrtS**2*xd1u + 4*sqrtS**4*xd1u - 4*M0Lpu**2*uman*xd1u + 4*m1**2*uman*xd1u + 4*M1Lpu**2*uman*xd1u + 
     &       4*m2**2*uman*xd1u + 4*M3Lpu**2*uman*xd1u - 8*sqrtS**2*uman*xd1u - 8*M0Lpu**2*uman*xd22u + 8*m1**2*uman*xd22u + 8*M1Lpu**2*uman*xd22u + 8*m2**2*uman*xd22u + 
     &       8*M3Lpu**2*uman*xd22u - 16*sqrtS**2*uman*xd22u - 8*M0Lpu**2*m1**2*xd23u + 8*m1**4*xd23u + 8*m1**2*M1Lpu**2*xd23u - 16*M0Lpu**2*m2**2*xd23u + 24*m1**2*m2**2*xd23u + 
     &       16*M1Lpu**2*m2**2*xd23u + 16*m2**4*xd23u + 8*m1**2*M3Lpu**2*xd23u + 16*m2**2*M3Lpu**2*xd23u - 8*M0Lpu**2*m4**2*xd23u + 8*m1**2*m4**2*xd23u + 8*M1Lpu**2*m4**2*xd23u + 
     &       8*m2**2*m4**2*xd23u + 8*M3Lpu**2*m4**2*xd23u + 8*M0Lpu**2*sqrtS**2*xd23u - 24*m1**2*sqrtS**2*xd23u - 8*M1Lpu**2*sqrtS**2*xd23u - 40*m2**2*sqrtS**2*xd23u - 
     &       8*M3Lpu**2*sqrtS**2*xd23u - 16*m4**2*sqrtS**2*xd23u + 16*sqrtS**4*xd23u + 8*M0Lpu**2*tman*xd23u - 8*m1**2*tman*xd23u - 8*M1Lpu**2*tman*xd23u - 8*m2**2*tman*xd23u - 
     &       8*M3Lpu**2*tman*xd23u + 16*sqrtS**2*tman*xd23u - 4*M0Lpu**2*m1**2*xd2u + 4*m1**4*xd2u + 4*m1**2*M1Lpu**2*xd2u - 4*M0Lpu**2*m2**2*xd2u + 8*m1**2*m2**2*xd2u + 
     &       4*M1Lpu**2*m2**2*xd2u + 4*m2**4*xd2u + 4*m1**2*M3Lpu**2*xd2u + 4*m2**2*M3Lpu**2*xd2u + 2*M0Lpu**2*sqrtS**2*xd2u - 10*m1**2*sqrtS**2*xd2u - 2*M1Lpu**2*sqrtS**2*xd2u - 
     &       10*m2**2*sqrtS**2*xd2u - 2*M3Lpu**2*sqrtS**2*xd2u + 4*sqrtS**4*xd2u + 2*M0Lpu**2*tman*xd2u - 2*m1**2*tman*xd2u - 2*M1Lpu**2*tman*xd2u - 2*m2**2*tman*xd2u - 
     &       2*M3Lpu**2*tman*xd2u + 4*sqrtS**2*tman*xd2u - 10*M0Lpu**2*uman*xd2u + 10*m1**2*uman*xd2u + 10*M1Lpu**2*uman*xd2u + 10*m2**2*uman*xd2u + 10*M3Lpu**2*uman*xd2u - 
     &       20*sqrtS**2*uman*xd2u - 8*M0Lpu**2*m2**2*xd33u + 8*m1**2*m2**2*xd33u + 8*M1Lpu**2*m2**2*xd33u + 8*m2**4*xd33u + 8*m2**2*M3Lpu**2*xd33u - 16*m2**2*sqrtS**2*xd33u - 
     &       6*M0Lpu**2*m1**2*xd3u + 6*m1**4*xd3u + 6*m1**2*M1Lpu**2*xd3u - 14*M0Lpu**2*m2**2*xd3u + 20*m1**2*m2**2*xd3u + 14*M1Lpu**2*m2**2*xd3u + 14*m2**4*xd3u + 6*m1**2*M3Lpu**2*xd3u + 
     &       14*m2**2*M3Lpu**2*xd3u - 4*M0Lpu**2*m4**2*xd3u + 4*m1**2*m4**2*xd3u + 4*M1Lpu**2*m4**2*xd3u + 4*m2**2*m4**2*xd3u + 4*M3Lpu**2*m4**2*xd3u + 6*M0Lpu**2*sqrtS**2*xd3u - 
     &       18*m1**2*sqrtS**2*xd3u - 6*M1Lpu**2*sqrtS**2*xd3u - 34*m2**2*sqrtS**2*xd3u - 6*M3Lpu**2*sqrtS**2*xd3u - 8*m4**2*sqrtS**2*xd3u + 12*sqrtS**4*xd3u + 4*M0Lpu**2*tman*xd3u - 
     &       4*m1**2*tman*xd3u - 4*M1Lpu**2*tman*xd3u - 4*m2**2*tman*xd3u - 4*M3Lpu**2*tman*xd3u + 8*sqrtS**2*tman*xd3u))/(32.*pi**2)  )

!==================================================================================================================================   
!Box(t)*s(Higgs) + Box(u)*s(Higgs)  
       
        do j=sumstartH, sumendH

        Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(LLR(xtype,ytype,j))
        propac = dconjg(DenSmH(j))
        
	Box_TS = Box_TS + 2*colorfactorSH*dreal(  (Coupl_Ac*Coupl_Bc*g0t*g1t*g2t*g3t*(-32*xc001t + 32*xc002t + 32*xc003t - 2*sqrtS**2*xc01t + 3*tman*xc01t - uman*xc01t - 2*sqrtS**2*xc02t + 3*tman*xc02t - uman*xc02t + tman*xc03t + 
     &      uman*xc03t + 8*tman*xc111t - 2*sqrtS**2*xc11t + 10*tman*xc11t - 2*uman*xc11t + 8*tman*xc121t + 8*uman*xc121t - 8*sqrtS**2*xc122t - 8*uman*xc122t + 8*tman*xc123t - 
     &      2*sqrtS**2*xc12t + 4*tman*xc12t + 2*tman*xc13t + 2*uman*xc13t - 2*sqrtS**2*xc21t + 10*tman*xc21t - 2*uman*xc21t - 8*sqrtS**2*xc221t + 8*tman*xc223t - 6*sqrtS**2*xc22t - 
     &      4*uman*xc22t - 4*sqrtS**2*xc23t + 4*tman*xc23t + 4*uman*xc23t - 4*m3**2*(xc03t + 2*xc112t + 4*xc121t - 2*xc122t + 2*xc123t + xc12t + 3*xc13t - xc22t + 2*xc23t) + 
     &      32*(-M0Lpt**2 + M1Lpt**2 + M3Lpt**2)*xd00t - 64*sqrtS**2*xd00t + 2*M0Lpt**2*sqrtS**2*xd0t - 2*M1Lpt**2*sqrtS**2*xd0t - 2*M3Lpt**2*sqrtS**2*xd0t + 4*sqrtS**4*xd0t - 
     &      3*M0Lpt**2*tman*xd0t + 3*M1Lpt**2*tman*xd0t + 3*M3Lpt**2*tman*xd0t - 6*sqrtS**2*tman*xd0t + M0Lpt**2*uman*xd0t - M1Lpt**2*uman*xd0t - M3Lpt**2*uman*xd0t + 2*sqrtS**2*uman*xd0t + 
     &      8*M0Lpt**2*m3**2*xd12t - 8*M1Lpt**2*m3**2*xd12t - 8*m3**2*M3Lpt**2*xd12t + 16*m3**2*sqrtS**2*xd12t - 8*M0Lpt**2*tman*xd12t + 8*M1Lpt**2*tman*xd12t + 8*M3Lpt**2*tman*xd12t - 
     &      16*sqrtS**2*tman*xd12t + 8*M0Lpt**2*sqrtS**2*xd13t - 8*M1Lpt**2*sqrtS**2*xd13t - 8*M3Lpt**2*sqrtS**2*xd13t + 16*sqrtS**4*xd13t + 4*M0Lpt**2*m3**2*xd1t - 4*M1Lpt**2*m3**2*xd1t - 
     &      4*m3**2*M3Lpt**2*xd1t + 2*M0Lpt**2*sqrtS**2*xd1t - 2*M1Lpt**2*sqrtS**2*xd1t + 8*m3**2*sqrtS**2*xd1t - 2*M3Lpt**2*sqrtS**2*xd1t + 4*sqrtS**4*xd1t - 4*M0Lpt**2*tman*xd1t + 
     &      4*M1Lpt**2*tman*xd1t + 4*M3Lpt**2*tman*xd1t - 8*sqrtS**2*tman*xd1t - 8*M0Lpt**2*tman*xd22t + 8*M1Lpt**2*tman*xd22t + 8*M3Lpt**2*tman*xd22t - 16*sqrtS**2*tman*xd22t - 
     &      8*M0Lpt**2*m3**2*xd23t + 8*M1Lpt**2*m3**2*xd23t + 8*m3**2*M3Lpt**2*xd23t + 8*M0Lpt**2*sqrtS**2*xd23t - 8*M1Lpt**2*sqrtS**2*xd23t - 16*m3**2*sqrtS**2*xd23t - 
     &      8*M3Lpt**2*sqrtS**2*xd23t + 16*sqrtS**4*xd23t + 8*M0Lpt**2*uman*xd23t - 8*M1Lpt**2*uman*xd23t - 8*M3Lpt**2*uman*xd23t + 16*sqrtS**2*uman*xd23t + 2*M0Lpt**2*sqrtS**2*xd2t - 
     &      2*M1Lpt**2*sqrtS**2*xd2t - 2*M3Lpt**2*sqrtS**2*xd2t + 4*sqrtS**4*xd2t - 10*M0Lpt**2*tman*xd2t + 10*M1Lpt**2*tman*xd2t + 10*M3Lpt**2*tman*xd2t - 20*sqrtS**2*tman*xd2t + 
     &      2*M0Lpt**2*uman*xd2t - 2*M1Lpt**2*uman*xd2t - 2*M3Lpt**2*uman*xd2t + 4*sqrtS**2*uman*xd2t + 
     &      2*(M0Lpt**2 - M1Lpt**2 - M3Lpt**2 + 2*sqrtS**2)*(-2*m3**2 + 3*sqrtS**2 + 2*uman)*xd3t + m1**4*(3*xd0t + 8*xd11t + 8*xd12t + 8*xd13t + 10*xd1t + 8*xd23t + 4*xd2t + 6*xd3t) + 
     &      m1**2*(3*xc01t + 3*xc02t - xc03t + 8*xc113t + 4*xc11t - 8*xc121t + 8*xc122t + 8*xc123t + 10*xc12t - 2*xc13t + 4*xc21t + 6*xc22t - 4*xc23t + 32*xd00t + 6*m2**2*xd0t + 
     &         3*M3Lpt**2*xd0t - 8*sqrtS**2*xd0t + 3*tman*xd0t - uman*xd0t + 8*m2**2*xd11t + 8*M3Lpt**2*xd11t - 16*sqrtS**2*xd11t + 8*m2**2*xd12t - 8*m3**2*xd12t + 8*M3Lpt**2*xd12t - 
     &         16*sqrtS**2*xd12t + 8*tman*xd12t + 16*m2**2*xd13t + 8*M3Lpt**2*xd13t - 24*sqrtS**2*xd13t + 12*m2**2*xd1t - 4*m3**2*xd1t + 10*M3Lpt**2*xd1t - 22*sqrtS**2*xd1t + 4*tman*xd1t + 
     &         8*tman*xd22t + 24*m2**2*xd23t + 8*m3**2*xd23t + 8*M3Lpt**2*xd23t - 24*sqrtS**2*xd23t - 8*uman*xd23t + 8*m2**2*xd2t + 4*M3Lpt**2*xd2t - 10*sqrtS**2*xd2t + 10*tman*xd2t - 
     &         2*uman*xd2t + 8*m2**2*xd33t + 2*(10*m2**2 + 2*m3**2 + 3*M3Lpt**2 - 9*sqrtS**2 - 2*uman)*xd3t - 
     &         M0Lpt**2*(3*xd0t + 8*xd11t + 8*xd12t + 8*xd13t + 10*xd1t + 8*xd23t + 4*xd2t + 6*xd3t) + M1Lpt**2*(3*xd0t + 8*xd11t + 8*xd12t + 8*xd13t + 10*xd1t + 8*xd23t + 4*xd2t + 6*xd3t)) 
     &       + m2**4*(3*xd0t + 2*(4*xd13t + xd1t + 8*xd23t + 2*xd2t + 4*xd33t + 7*xd3t)) + 
     &      m2**2*(3*xc01t + 3*xc02t - xc03t + 4*xc11t - 8*xc121t + 16*xc122t + 2*xc12t - 2*xc13t + 4*xc21t + 8*xc222t + 14*xc22t - 4*xc23t + 32*xd00t + 3*M3Lpt**2*xd0t - 8*sqrtS**2*xd0t + 
     &         3*tman*xd0t - uman*xd0t + 8*tman*xd12t + 8*M3Lpt**2*xd13t - 24*sqrtS**2*xd13t + 2*M3Lpt**2*xd1t - 6*sqrtS**2*xd1t + 4*tman*xd1t + 8*tman*xd22t + 16*M3Lpt**2*xd23t - 
     &         40*sqrtS**2*xd23t - 8*uman*xd23t + 4*M3Lpt**2*xd2t - 10*sqrtS**2*xd2t + 10*tman*xd2t - 2*uman*xd2t + 8*M3Lpt**2*xd33t - 16*sqrtS**2*xd33t + 14*M3Lpt**2*xd3t - 
     &         34*sqrtS**2*xd3t - 4*uman*xd3t + 4*m3**2*(-2*xd12t - xd1t + 2*xd23t + xd3t) - M0Lpt**2*(3*xd0t + 2*(4*xd13t + xd1t + 8*xd23t + 2*xd2t + 4*xd33t + 7*xd3t)) + 
     &         M1Lpt**2*(3*xd0t + 2*(4*xd13t + xd1t + 8*xd23t + 2*xd2t + 4*xd33t + 7*xd3t)))))/(32.*pi**2*propac)  )

	Box_US = Box_US + 2*colorfactorSH*dreal(   (Coupl_Ac*Coupl_Bc*g0u*g1u*g2u*g3u*(-32*xc001u + 32*xc002u + 32*xc003u - 2*sqrtS**2*xc01u - tman*xc01u + 3*uman*xc01u - 2*sqrtS**2*xc02u - tman*xc02u + 3*uman*xc02u + tman*xc03u + 
     &      uman*xc03u + 8*uman*xc111u - 2*sqrtS**2*xc11u - 2*tman*xc11u + 10*uman*xc11u + 8*tman*xc121u + 8*uman*xc121u - 8*sqrtS**2*xc122u - 8*tman*xc122u + 8*uman*xc123u - 
     &      2*sqrtS**2*xc12u + 4*uman*xc12u + 2*tman*xc13u + 2*uman*xc13u - 2*sqrtS**2*xc21u - 2*tman*xc21u + 10*uman*xc21u - 8*sqrtS**2*xc221u + 8*uman*xc223u - 6*sqrtS**2*xc22u - 
     &      4*tman*xc22u - 4*sqrtS**2*xc23u + 4*tman*xc23u + 4*uman*xc23u - 4*m4**2*(xc03u + 2*xc112u + 4*xc121u - 2*xc122u + 2*xc123u + xc12u + 3*xc13u - xc22u + 2*xc23u) + 
     &      32*(-M0Lpu**2 + M1Lpu**2 + M3Lpu**2)*xd00u - 64*sqrtS**2*xd00u + 2*M0Lpu**2*sqrtS**2*xd0u - 2*M1Lpu**2*sqrtS**2*xd0u - 2*M3Lpu**2*sqrtS**2*xd0u + 4*sqrtS**4*xd0u + 
     &      M0Lpu**2*tman*xd0u - M1Lpu**2*tman*xd0u - M3Lpu**2*tman*xd0u + 2*sqrtS**2*tman*xd0u - 3*M0Lpu**2*uman*xd0u + 3*M1Lpu**2*uman*xd0u + 3*M3Lpu**2*uman*xd0u - 6*sqrtS**2*uman*xd0u + 
     &      8*M0Lpu**2*m4**2*xd12u - 8*M1Lpu**2*m4**2*xd12u - 8*M3Lpu**2*m4**2*xd12u + 16*m4**2*sqrtS**2*xd12u - 8*M0Lpu**2*uman*xd12u + 8*M1Lpu**2*uman*xd12u + 8*M3Lpu**2*uman*xd12u - 
     &      16*sqrtS**2*uman*xd12u + 8*M0Lpu**2*sqrtS**2*xd13u - 8*M1Lpu**2*sqrtS**2*xd13u - 8*M3Lpu**2*sqrtS**2*xd13u + 16*sqrtS**4*xd13u + 4*M0Lpu**2*m4**2*xd1u - 4*M1Lpu**2*m4**2*xd1u - 
     &      4*M3Lpu**2*m4**2*xd1u + 2*M0Lpu**2*sqrtS**2*xd1u - 2*M1Lpu**2*sqrtS**2*xd1u - 2*M3Lpu**2*sqrtS**2*xd1u + 8*m4**2*sqrtS**2*xd1u + 4*sqrtS**4*xd1u - 4*M0Lpu**2*uman*xd1u + 
     &      4*M1Lpu**2*uman*xd1u + 4*M3Lpu**2*uman*xd1u - 8*sqrtS**2*uman*xd1u - 8*M0Lpu**2*uman*xd22u + 8*M1Lpu**2*uman*xd22u + 8*M3Lpu**2*uman*xd22u - 16*sqrtS**2*uman*xd22u - 
     &      8*M0Lpu**2*m4**2*xd23u + 8*M1Lpu**2*m4**2*xd23u + 8*M3Lpu**2*m4**2*xd23u + 8*M0Lpu**2*sqrtS**2*xd23u - 8*M1Lpu**2*sqrtS**2*xd23u - 8*M3Lpu**2*sqrtS**2*xd23u - 
     &      16*m4**2*sqrtS**2*xd23u + 16*sqrtS**4*xd23u + 8*M0Lpu**2*tman*xd23u - 8*M1Lpu**2*tman*xd23u - 8*M3Lpu**2*tman*xd23u + 16*sqrtS**2*tman*xd23u + 2*M0Lpu**2*sqrtS**2*xd2u - 
     &      2*M1Lpu**2*sqrtS**2*xd2u - 2*M3Lpu**2*sqrtS**2*xd2u + 4*sqrtS**4*xd2u + 2*M0Lpu**2*tman*xd2u - 2*M1Lpu**2*tman*xd2u - 2*M3Lpu**2*tman*xd2u + 4*sqrtS**2*tman*xd2u - 
     &      10*M0Lpu**2*uman*xd2u + 10*M1Lpu**2*uman*xd2u + 10*M3Lpu**2*uman*xd2u - 20*sqrtS**2*uman*xd2u + 
     &      2*(M0Lpu**2 - M1Lpu**2 - M3Lpu**2 + 2*sqrtS**2)*(-2*m4**2 + 3*sqrtS**2 + 2*tman)*xd3u + m1**4*(3*xd0u + 8*xd11u + 8*xd12u + 8*xd13u + 10*xd1u + 8*xd23u + 4*xd2u + 6*xd3u) + 
     &      m1**2*(3*xc01u + 3*xc02u - xc03u + 8*xc113u + 4*xc11u - 8*xc121u + 8*xc122u + 8*xc123u + 10*xc12u - 2*xc13u + 4*xc21u + 6*xc22u - 4*xc23u + 32*xd00u + 6*m2**2*xd0u + 
     &         3*M3Lpu**2*xd0u - 8*sqrtS**2*xd0u - tman*xd0u + 3*uman*xd0u + 8*m2**2*xd11u + 8*M3Lpu**2*xd11u - 16*sqrtS**2*xd11u + 8*m2**2*xd12u + 8*M3Lpu**2*xd12u - 8*m4**2*xd12u - 
     &         16*sqrtS**2*xd12u + 8*uman*xd12u + 16*m2**2*xd13u + 8*M3Lpu**2*xd13u - 24*sqrtS**2*xd13u + 12*m2**2*xd1u + 10*M3Lpu**2*xd1u - 4*m4**2*xd1u - 22*sqrtS**2*xd1u + 4*uman*xd1u + 
     &         8*uman*xd22u + 24*m2**2*xd23u + 8*M3Lpu**2*xd23u + 8*m4**2*xd23u - 24*sqrtS**2*xd23u - 8*tman*xd23u + 8*m2**2*xd2u + 4*M3Lpu**2*xd2u - 10*sqrtS**2*xd2u - 2*tman*xd2u + 
     &         10*uman*xd2u + 8*m2**2*xd33u + 2*(10*m2**2 + 3*M3Lpu**2 + 2*m4**2 - 9*sqrtS**2 - 2*tman)*xd3u - 
     &         M0Lpu**2*(3*xd0u + 8*xd11u + 8*xd12u + 8*xd13u + 10*xd1u + 8*xd23u + 4*xd2u + 6*xd3u) + M1Lpu**2*(3*xd0u + 8*xd11u + 8*xd12u + 8*xd13u + 10*xd1u + 8*xd23u + 4*xd2u + 6*xd3u)) 
     &       + m2**4*(3*xd0u + 2*(4*xd13u + xd1u + 8*xd23u + 2*xd2u + 4*xd33u + 7*xd3u)) + 
     &      m2**2*(3*xc01u + 3*xc02u - xc03u + 4*xc11u - 8*xc121u + 16*xc122u + 2*xc12u - 2*xc13u + 4*xc21u + 8*xc222u + 14*xc22u - 4*xc23u + 32*xd00u + 3*M3Lpu**2*xd0u - 8*sqrtS**2*xd0u - 
     &         tman*xd0u + 3*uman*xd0u - 8*m4**2*xd12u + 8*uman*xd12u + 8*M3Lpu**2*xd13u - 24*sqrtS**2*xd13u + 2*M3Lpu**2*xd1u - 4*m4**2*xd1u - 6*sqrtS**2*xd1u + 4*uman*xd1u + 8*uman*xd22u + 
     &         16*M3Lpu**2*xd23u + 8*m4**2*xd23u - 40*sqrtS**2*xd23u - 8*tman*xd23u + 4*M3Lpu**2*xd2u - 10*sqrtS**2*xd2u - 2*tman*xd2u + 10*uman*xd2u + 8*M3Lpu**2*xd33u - 16*sqrtS**2*xd33u + 
     &         2*(7*M3Lpu**2 + 2*m4**2 - 17*sqrtS**2 - 2*tman)*xd3u - M0Lpu**2*(3*xd0u + 2*(4*xd13u + xd1u + 8*xd23u + 2*xd2u + 4*xd33u + 7*xd3u)) + 
     &         M1Lpu**2*(3*xd0u + 2*(4*xd13u + xd1u + 8*xd23u + 2*xd2u + 4*xd33u + 7*xd3u)))))/(32.*pi**2*propac)  )
  
        enddo

!==================================================================================================================================  
!Box(t)*Vector + Box(u)*Vector

        do j=sumstartV, sumendV

        Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = dconjg(VVV(j,xtype,ytype))
        propac = dconjg(DenSmV(j))
        
        Box_TV = Box_TV + 2*colorfactorSV*dreal(  -(Coupl_Ac*Coupl_Bc*g0t*g1t*g2t*g3t*(m1**4*xc01t - 20*m1**2*m2**2*xc01t + 3*m2**4*xc01t - 17*m1**2*m3**2*xc01t - 7*m2**2*m3**2*xc01t - 9*m1**2*m4**2*xc01t - 15*m2**2*m4**2*xc01t - 
     &       m1**2*sqrtS**2*xc01t - 3*m2**2*sqrtS**2*xc01t + 6*m3**2*sqrtS**2*xc01t + 6*m4**2*sqrtS**2*xc01t + 2*sqrtS**4*xc01t + m1**4*xc02t - 20*m1**2*m2**2*xc02t + 3*m2**4*xc02t - 
     &       17*m1**2*m3**2*xc02t - 7*m2**2*m3**2*xc02t - 9*m1**2*m4**2*xc02t - 15*m2**2*m4**2*xc02t - m1**2*sqrtS**2*xc02t - 3*m2**2*sqrtS**2*xc02t + 6*m3**2*sqrtS**2*xc02t + 
     &       6*m4**2*sqrtS**2*xc02t + 2*sqrtS**4*xc02t + m1**4*xc03t - m2**4*xc03t - 3*m1**2*m3**2*xc03t + 3*m2**2*m3**2*xc03t - 3*m1**2*m4**2*xc03t + 3*m2**2*m4**2*xc03t + 
     &       m1**2*sqrtS**2*xc03t - m2**2*sqrtS**2*xc03t - 16*m1**2*m3**2*xc111t - 16*m2**2*m3**2*xc111t - 16*m1**2*m4**2*xc111t - 16*m2**2*m4**2*xc111t + 8*m1**2*sqrtS**2*xc111t + 
     &       8*m2**2*sqrtS**2*xc111t + 8*m3**2*sqrtS**2*xc111t + 8*m4**2*sqrtS**2*xc111t - 4*sqrtS**4*xc111t - 12*m1**2*m3**2*xc112t + 12*m2**2*m3**2*xc112t - 4*m1**2*m4**2*xc112t + 
     &       4*m2**2*m4**2*xc112t + 4*m1**2*sqrtS**2*xc112t - 4*m2**2*sqrtS**2*xc112t + 12*m1**2*m3**2*xc113t + 4*m2**2*m3**2*xc113t - 12*m1**2*m4**2*xc113t - 4*m2**2*m4**2*xc113t - 
     &       4*m3**2*sqrtS**2*xc113t + 4*m4**2*sqrtS**2*xc113t - 20*m1**2*m2**2*xc11t + 4*m2**4*xc11t - 36*m1**2*m3**2*xc11t - 20*m2**2*m3**2*xc11t - 24*m1**2*m4**2*xc11t - 
     &       32*m2**2*m4**2*xc11t + 8*m1**2*sqrtS**2*xc11t + 4*m2**2*sqrtS**2*xc11t + 14*m3**2*sqrtS**2*xc11t + 14*m4**2*sqrtS**2*xc11t - 2*sqrtS**4*xc11t - 32*m1**2*m3**2*xc121t + 
     &       32*m2**2*m3**2*xc121t + 8*m1**2*sqrtS**2*xc121t - 8*m2**2*sqrtS**2*xc121t - 16*m1**2*m3**2*xc122t - 32*m2**2*m3**2*xc122t - 16*m2**2*m4**2*xc122t + 4*m1**2*sqrtS**2*xc122t + 
     &       12*m2**2*sqrtS**2*xc122t + 12*m3**2*sqrtS**2*xc122t + 4*m4**2*sqrtS**2*xc122t - 4*sqrtS**4*xc122t - 16*m1**2*m3**2*xc123t - 32*m1**2*m4**2*xc123t - 16*m2**2*m4**2*xc123t + 
     &       12*m1**2*sqrtS**2*xc123t + 4*m2**2*sqrtS**2*xc123t + 4*m3**2*sqrtS**2*xc123t + 12*m4**2*sqrtS**2*xc123t - 4*sqrtS**4*xc123t + 2*m1**4*xc12t - 20*m1**2*m2**2*xc12t + 
     &       2*m2**4*xc12t - 8*m1**2*m3**2*xc12t - 24*m1**2*m4**2*xc12t - 16*m2**2*m4**2*xc12t - 4*m2**2*sqrtS**2*xc12t + 2*m3**2*sqrtS**2*xc12t + 10*m4**2*sqrtS**2*xc12t + 
     &       2*sqrtS**4*xc12t + 2*m1**4*xc13t - 2*m2**4*xc13t - 12*m1**2*m3**2*xc13t + 12*m2**2*m3**2*xc13t - 8*m1**2*m4**2*xc13t + 8*m2**2*m4**2*xc13t + 4*m1**2*sqrtS**2*xc13t - 
     &       4*m2**2*sqrtS**2*xc13t - 20*m1**2*m2**2*xc21t + 4*m2**4*xc21t - 36*m1**2*m3**2*xc21t - 20*m2**2*m3**2*xc21t - 24*m1**2*m4**2*xc21t - 32*m2**2*m4**2*xc21t + 
     &       8*m1**2*sqrtS**2*xc21t + 4*m2**2*sqrtS**2*xc21t + 14*m3**2*sqrtS**2*xc21t + 14*m4**2*sqrtS**2*xc21t - 2*sqrtS**4*xc21t - 16*m1**2*m3**2*xc221t + 16*m2**2*m3**2*xc221t + 
     &       16*m1**2*m4**2*xc221t - 16*m2**2*m4**2*xc221t - 4*m1**2*m3**2*xc222t - 12*m2**2*m3**2*xc222t + 4*m1**2*m4**2*xc222t + 12*m2**2*m4**2*xc222t + 4*m3**2*sqrtS**2*xc222t - 
     &       4*m4**2*sqrtS**2*xc222t - 16*m1**2*m3**2*xc223t - 16*m2**2*m3**2*xc223t - 16*m1**2*m4**2*xc223t - 16*m2**2*m4**2*xc223t + 8*m1**2*sqrtS**2*xc223t + 8*m2**2*sqrtS**2*xc223t + 
     &       8*m3**2*sqrtS**2*xc223t + 8*m4**2*sqrtS**2*xc223t - 4*sqrtS**4*xc223t + 2*m1**4*xc22t - 20*m1**2*m2**2*xc22t + 2*m2**4*xc22t - 16*m1**2*m3**2*xc22t - 24*m2**2*m3**2*xc22t - 
     &       8*m2**2*m4**2*xc22t - 4*m1**2*sqrtS**2*xc22t + 10*m3**2*sqrtS**2*xc22t + 2*m4**2*sqrtS**2*xc22t + 2*sqrtS**4*xc22t - 24*m1**2*m3**2*xc23t + 24*m2**2*m3**2*xc23t + 
     &       8*m1**2*m4**2*xc23t - 8*m2**2*m4**2*xc23t + 4*m1**2*sqrtS**2*xc23t - 4*m2**2*sqrtS**2*xc23t - M0Lpt**2*m1**4*xd0t + m1**6*xd0t + m1**4*M1Lpt**2*xd0t + 
     &       20*M0Lpt**2*m1**2*m2**2*xd0t - 19*m1**4*m2**2*xd0t - 20*m1**2*M1Lpt**2*m2**2*xd0t - 3*M0Lpt**2*m2**4*xd0t - 17*m1**2*m2**4*xd0t + 3*M1Lpt**2*m2**4*xd0t + 3*m2**6*xd0t + 
     &       17*M0Lpt**2*m1**2*m3**2*xd0t - 17*m1**4*m3**2*xd0t - 17*m1**2*M1Lpt**2*m3**2*xd0t + 7*M0Lpt**2*m2**2*m3**2*xd0t - 24*m1**2*m2**2*m3**2*xd0t - 7*M1Lpt**2*m2**2*m3**2*xd0t - 
     &       7*m2**4*m3**2*xd0t + m1**4*M3Lpt**2*xd0t - 20*m1**2*m2**2*M3Lpt**2*xd0t + 3*m2**4*M3Lpt**2*xd0t - 17*m1**2*m3**2*M3Lpt**2*xd0t - 7*m2**2*m3**2*M3Lpt**2*xd0t + 
     &       9*M0Lpt**2*m1**2*m4**2*xd0t - 9*m1**4*m4**2*xd0t - 9*m1**2*M1Lpt**2*m4**2*xd0t + 15*M0Lpt**2*m2**2*m4**2*xd0t - 24*m1**2*m2**2*m4**2*xd0t - 15*M1Lpt**2*m2**2*m4**2*xd0t - 
     &       15*m2**4*m4**2*xd0t - 9*m1**2*M3Lpt**2*m4**2*xd0t - 15*m2**2*M3Lpt**2*m4**2*xd0t + M0Lpt**2*m1**2*sqrtS**2*xd0t - 3*m1**4*sqrtS**2*xd0t - m1**2*M1Lpt**2*sqrtS**2*xd0t + 
     &       3*M0Lpt**2*m2**2*sqrtS**2*xd0t + 36*m1**2*m2**2*sqrtS**2*xd0t - 3*M1Lpt**2*m2**2*sqrtS**2*xd0t - 9*m2**4*sqrtS**2*xd0t - 6*M0Lpt**2*m3**2*sqrtS**2*xd0t + 
     &       40*m1**2*m3**2*sqrtS**2*xd0t + 6*M1Lpt**2*m3**2*sqrtS**2*xd0t + 20*m2**2*m3**2*sqrtS**2*xd0t - m1**2*M3Lpt**2*sqrtS**2*xd0t - 3*m2**2*M3Lpt**2*sqrtS**2*xd0t + 
     &       6*m3**2*M3Lpt**2*sqrtS**2*xd0t - 6*M0Lpt**2*m4**2*sqrtS**2*xd0t + 24*m1**2*m4**2*sqrtS**2*xd0t + 6*M1Lpt**2*m4**2*sqrtS**2*xd0t + 36*m2**2*m4**2*sqrtS**2*xd0t + 
     &       6*M3Lpt**2*m4**2*sqrtS**2*xd0t - 2*M0Lpt**2*sqrtS**4*xd0t + 4*m1**2*sqrtS**4*xd0t + 2*M1Lpt**2*sqrtS**4*xd0t + 8*m2**2*sqrtS**4*xd0t - 12*m3**2*sqrtS**4*xd0t + 
     &       2*M3Lpt**2*sqrtS**4*xd0t - 12*m4**2*sqrtS**4*xd0t - 4*sqrtS**6*xd0t - 12*M0Lpt**2*m1**2*m3**2*xd11t + 12*m1**4*m3**2*xd11t + 12*m1**2*M1Lpt**2*m3**2*xd11t - 
     &       4*M0Lpt**2*m2**2*m3**2*xd11t + 16*m1**2*m2**2*m3**2*xd11t + 4*M1Lpt**2*m2**2*m3**2*xd11t + 4*m2**4*m3**2*xd11t + 12*m1**2*m3**2*M3Lpt**2*xd11t + 4*m2**2*m3**2*M3Lpt**2*xd11t + 
     &       12*M0Lpt**2*m1**2*m4**2*xd11t - 12*m1**4*m4**2*xd11t - 12*m1**2*M1Lpt**2*m4**2*xd11t + 4*M0Lpt**2*m2**2*m4**2*xd11t - 16*m1**2*m2**2*m4**2*xd11t - 4*M1Lpt**2*m2**2*m4**2*xd11t - 
     &       4*m2**4*m4**2*xd11t - 12*m1**2*M3Lpt**2*m4**2*xd11t - 4*m2**2*M3Lpt**2*m4**2*xd11t + 4*M0Lpt**2*m3**2*sqrtS**2*xd11t - 28*m1**2*m3**2*sqrtS**2*xd11t - 
     &       4*M1Lpt**2*m3**2*sqrtS**2*xd11t - 12*m2**2*m3**2*sqrtS**2*xd11t - 4*m3**2*M3Lpt**2*sqrtS**2*xd11t - 4*M0Lpt**2*m4**2*sqrtS**2*xd11t + 28*m1**2*m4**2*sqrtS**2*xd11t + 
     &       4*M1Lpt**2*m4**2*sqrtS**2*xd11t + 12*m2**2*m4**2*sqrtS**2*xd11t + 4*M3Lpt**2*m4**2*sqrtS**2*xd11t + 8*m3**2*sqrtS**4*xd11t - 8*m4**2*sqrtS**4*xd11t + 
     &       16*M0Lpt**2*m1**2*m3**2*xd12t - 16*m1**4*m3**2*xd12t - 16*m1**2*M1Lpt**2*m3**2*xd12t - 16*m1**2*m2**2*m3**2*xd12t - 16*m1**2*m3**2*M3Lpt**2*xd12t + 
     &       32*M0Lpt**2*m1**2*m4**2*xd12t - 32*m1**4*m4**2*xd12t - 32*m1**2*M1Lpt**2*m4**2*xd12t + 16*M0Lpt**2*m2**2*m4**2*xd12t - 48*m1**2*m2**2*m4**2*xd12t - 
     &       16*M1Lpt**2*m2**2*m4**2*xd12t - 16*m2**4*m4**2*xd12t - 32*m1**2*M3Lpt**2*m4**2*xd12t - 16*m2**2*M3Lpt**2*m4**2*xd12t - 12*M0Lpt**2*m1**2*sqrtS**2*xd12t + 
     &       12*m1**4*sqrtS**2*xd12t + 12*m1**2*M1Lpt**2*sqrtS**2*xd12t - 4*M0Lpt**2*m2**2*sqrtS**2*xd12t + 16*m1**2*m2**2*sqrtS**2*xd12t + 4*M1Lpt**2*m2**2*sqrtS**2*xd12t + 
     &       4*m2**4*sqrtS**2*xd12t - 4*M0Lpt**2*m3**2*sqrtS**2*xd12t + 36*m1**2*m3**2*sqrtS**2*xd12t + 4*M1Lpt**2*m3**2*sqrtS**2*xd12t + 4*m2**2*m3**2*sqrtS**2*xd12t + 
     &       12*m1**2*M3Lpt**2*sqrtS**2*xd12t + 4*m2**2*M3Lpt**2*sqrtS**2*xd12t + 4*m3**2*M3Lpt**2*sqrtS**2*xd12t - 12*M0Lpt**2*m4**2*sqrtS**2*xd12t + 76*m1**2*m4**2*sqrtS**2*xd12t + 
     &       12*M1Lpt**2*m4**2*sqrtS**2*xd12t + 44*m2**2*m4**2*sqrtS**2*xd12t + 12*M3Lpt**2*m4**2*sqrtS**2*xd12t + 4*M0Lpt**2*sqrtS**4*xd12t - 28*m1**2*sqrtS**4*xd12t - 
     &       4*M1Lpt**2*sqrtS**4*xd12t - 12*m2**2*sqrtS**4*xd12t - 8*m3**2*sqrtS**4*xd12t - 4*M3Lpt**2*sqrtS**4*xd12t - 24*m4**2*sqrtS**4*xd12t + 8*sqrtS**6*xd12t + 
     &       8*M0Lpt**2*m1**2*m3**2*xd13t - 8*m1**4*m3**2*xd13t - 8*m1**2*M1Lpt**2*m3**2*xd13t - 8*M0Lpt**2*m2**2*m3**2*xd13t + 8*M1Lpt**2*m2**2*m3**2*xd13t + 8*m2**4*m3**2*xd13t - 
     &       8*m1**2*m3**2*M3Lpt**2*xd13t + 8*m2**2*m3**2*M3Lpt**2*xd13t - 8*M0Lpt**2*m1**2*m4**2*xd13t + 8*m1**4*m4**2*xd13t + 8*m1**2*M1Lpt**2*m4**2*xd13t + 8*M0Lpt**2*m2**2*m4**2*xd13t - 
     &       8*M1Lpt**2*m2**2*m4**2*xd13t - 8*m2**4*m4**2*xd13t + 8*m1**2*M3Lpt**2*m4**2*xd13t - 8*m2**2*M3Lpt**2*m4**2*xd13t + 16*m1**2*m3**2*sqrtS**2*xd13t - 
     &       16*m2**2*m3**2*sqrtS**2*xd13t - 16*m1**2*m4**2*sqrtS**2*xd13t + 16*m2**2*m4**2*sqrtS**2*xd13t - 2*M0Lpt**2*m1**4*xd1t + 2*m1**6*xd1t + 2*m1**4*M1Lpt**2*xd1t + 
     &       20*M0Lpt**2*m1**2*m2**2*xd1t - 18*m1**4*m2**2*xd1t - 20*m1**2*M1Lpt**2*m2**2*xd1t - 2*M0Lpt**2*m2**4*xd1t - 18*m1**2*m2**4*xd1t + 2*M1Lpt**2*m2**4*xd1t + 2*m2**6*xd1t + 
     &       8*M0Lpt**2*m1**2*m3**2*xd1t - 8*m1**4*m3**2*xd1t - 8*m1**2*M1Lpt**2*m3**2*xd1t - 8*m1**2*m2**2*m3**2*xd1t + 2*m1**4*M3Lpt**2*xd1t - 20*m1**2*m2**2*M3Lpt**2*xd1t + 
     &       2*m2**4*M3Lpt**2*xd1t - 8*m1**2*m3**2*M3Lpt**2*xd1t + 24*M0Lpt**2*m1**2*m4**2*xd1t - 24*m1**4*m4**2*xd1t - 24*m1**2*M1Lpt**2*m4**2*xd1t + 16*M0Lpt**2*m2**2*m4**2*xd1t - 
     &       40*m1**2*m2**2*m4**2*xd1t - 16*M1Lpt**2*m2**2*m4**2*xd1t - 16*m2**4*m4**2*xd1t - 24*m1**2*M3Lpt**2*m4**2*xd1t - 16*m2**2*M3Lpt**2*m4**2*xd1t - 4*m1**4*sqrtS**2*xd1t + 
     &       4*M0Lpt**2*m2**2*sqrtS**2*xd1t + 36*m1**2*m2**2*sqrtS**2*xd1t - 4*M1Lpt**2*m2**2*sqrtS**2*xd1t - 8*m2**4*sqrtS**2*xd1t - 2*M0Lpt**2*m3**2*sqrtS**2*xd1t + 
     &       18*m1**2*m3**2*sqrtS**2*xd1t + 2*M1Lpt**2*m3**2*sqrtS**2*xd1t + 2*m2**2*m3**2*sqrtS**2*xd1t - 4*m2**2*M3Lpt**2*sqrtS**2*xd1t + 2*m3**2*M3Lpt**2*sqrtS**2*xd1t - 
     &       10*M0Lpt**2*m4**2*sqrtS**2*xd1t + 58*m1**2*m4**2*sqrtS**2*xd1t + 10*M1Lpt**2*m4**2*sqrtS**2*xd1t + 42*m2**2*m4**2*sqrtS**2*xd1t + 10*M3Lpt**2*m4**2*sqrtS**2*xd1t - 
     &       2*M0Lpt**2*sqrtS**4*xd1t + 2*m1**2*sqrtS**4*xd1t + 2*M1Lpt**2*sqrtS**4*xd1t + 10*m2**2*sqrtS**4*xd1t - 4*m3**2*sqrtS**4*xd1t + 2*M3Lpt**2*sqrtS**4*xd1t - 
     &       20*m4**2*sqrtS**4*xd1t - 4*sqrtS**6*xd1t + 16*M0Lpt**2*m1**2*m3**2*xd22t - 16*m1**4*m3**2*xd22t - 16*m1**2*M1Lpt**2*m3**2*xd22t + 16*M0Lpt**2*m2**2*m3**2*xd22t - 
     &       32*m1**2*m2**2*m3**2*xd22t - 16*M1Lpt**2*m2**2*m3**2*xd22t - 16*m2**4*m3**2*xd22t - 16*m1**2*m3**2*M3Lpt**2*xd22t - 16*m2**2*m3**2*M3Lpt**2*xd22t + 
     &       16*M0Lpt**2*m1**2*m4**2*xd22t - 16*m1**4*m4**2*xd22t - 16*m1**2*M1Lpt**2*m4**2*xd22t + 16*M0Lpt**2*m2**2*m4**2*xd22t - 32*m1**2*m2**2*m4**2*xd22t - 
     &       16*M1Lpt**2*m2**2*m4**2*xd22t - 16*m2**4*m4**2*xd22t - 16*m1**2*M3Lpt**2*m4**2*xd22t - 16*m2**2*M3Lpt**2*m4**2*xd22t - 8*M0Lpt**2*m1**2*sqrtS**2*xd22t + 8*m1**4*sqrtS**2*xd22t + 
     &       8*m1**2*M1Lpt**2*sqrtS**2*xd22t - 8*M0Lpt**2*m2**2*sqrtS**2*xd22t + 16*m1**2*m2**2*sqrtS**2*xd22t + 8*M1Lpt**2*m2**2*sqrtS**2*xd22t + 8*m2**4*sqrtS**2*xd22t - 
     &       8*M0Lpt**2*m3**2*sqrtS**2*xd22t + 40*m1**2*m3**2*sqrtS**2*xd22t + 8*M1Lpt**2*m3**2*sqrtS**2*xd22t + 40*m2**2*m3**2*sqrtS**2*xd22t + 8*m1**2*M3Lpt**2*sqrtS**2*xd22t + 
     &       8*m2**2*M3Lpt**2*sqrtS**2*xd22t + 8*m3**2*M3Lpt**2*sqrtS**2*xd22t - 8*M0Lpt**2*m4**2*sqrtS**2*xd22t + 40*m1**2*m4**2*sqrtS**2*xd22t + 8*M1Lpt**2*m4**2*sqrtS**2*xd22t + 
     &       40*m2**2*m4**2*sqrtS**2*xd22t + 8*M3Lpt**2*m4**2*sqrtS**2*xd22t + 4*M0Lpt**2*sqrtS**4*xd22t - 20*m1**2*sqrtS**4*xd22t - 4*M1Lpt**2*sqrtS**4*xd22t - 20*m2**2*sqrtS**4*xd22t - 
     &       16*m3**2*sqrtS**4*xd22t - 4*M3Lpt**2*sqrtS**4*xd22t - 16*m4**2*sqrtS**4*xd22t + 8*sqrtS**6*xd22t + 16*M0Lpt**2*m1**2*m3**2*xd23t - 16*m1**4*m3**2*xd23t - 
     &       16*m1**2*M1Lpt**2*m3**2*xd23t + 32*M0Lpt**2*m2**2*m3**2*xd23t - 48*m1**2*m2**2*m3**2*xd23t - 32*M1Lpt**2*m2**2*m3**2*xd23t - 32*m2**4*m3**2*xd23t - 
     &       16*m1**2*m3**2*M3Lpt**2*xd23t - 32*m2**2*m3**2*M3Lpt**2*xd23t + 16*M0Lpt**2*m2**2*m4**2*xd23t - 16*m1**2*m2**2*m4**2*xd23t - 16*M1Lpt**2*m2**2*m4**2*xd23t - 
     &       16*m2**4*m4**2*xd23t - 16*m2**2*M3Lpt**2*m4**2*xd23t - 4*M0Lpt**2*m1**2*sqrtS**2*xd23t + 4*m1**4*sqrtS**2*xd23t + 4*m1**2*M1Lpt**2*sqrtS**2*xd23t - 
     &       12*M0Lpt**2*m2**2*sqrtS**2*xd23t + 16*m1**2*m2**2*sqrtS**2*xd23t + 12*M1Lpt**2*m2**2*sqrtS**2*xd23t + 12*m2**4*sqrtS**2*xd23t - 12*M0Lpt**2*m3**2*sqrtS**2*xd23t + 
     &       44*m1**2*m3**2*sqrtS**2*xd23t + 12*M1Lpt**2*m3**2*sqrtS**2*xd23t + 76*m2**2*m3**2*sqrtS**2*xd23t + 4*m1**2*M3Lpt**2*sqrtS**2*xd23t + 12*m2**2*M3Lpt**2*sqrtS**2*xd23t + 
     &       12*m3**2*M3Lpt**2*sqrtS**2*xd23t - 4*M0Lpt**2*m4**2*sqrtS**2*xd23t + 4*m1**2*m4**2*sqrtS**2*xd23t + 4*M1Lpt**2*m4**2*sqrtS**2*xd23t + 36*m2**2*m4**2*sqrtS**2*xd23t + 
     &       4*M3Lpt**2*m4**2*sqrtS**2*xd23t + 4*M0Lpt**2*sqrtS**4*xd23t - 12*m1**2*sqrtS**4*xd23t - 4*M1Lpt**2*sqrtS**4*xd23t - 28*m2**2*sqrtS**4*xd23t - 24*m3**2*sqrtS**4*xd23t - 
     &       4*M3Lpt**2*sqrtS**4*xd23t - 8*m4**2*sqrtS**4*xd23t + 8*sqrtS**6*xd23t + 20*M0Lpt**2*m1**2*m2**2*xd2t - 20*m1**4*m2**2*xd2t - 20*m1**2*M1Lpt**2*m2**2*xd2t - 
     &       4*M0Lpt**2*m2**4*xd2t - 16*m1**2*m2**4*xd2t + 4*M1Lpt**2*m2**4*xd2t + 4*m2**6*xd2t + 36*M0Lpt**2*m1**2*m3**2*xd2t - 36*m1**4*m3**2*xd2t - 36*m1**2*M1Lpt**2*m3**2*xd2t + 
     &       20*M0Lpt**2*m2**2*m3**2*xd2t - 56*m1**2*m2**2*m3**2*xd2t - 20*M1Lpt**2*m2**2*m3**2*xd2t - 20*m2**4*m3**2*xd2t - 20*m1**2*m2**2*M3Lpt**2*xd2t + 4*m2**4*M3Lpt**2*xd2t - 
     &       36*m1**2*m3**2*M3Lpt**2*xd2t - 20*m2**2*m3**2*M3Lpt**2*xd2t + 24*M0Lpt**2*m1**2*m4**2*xd2t - 24*m1**4*m4**2*xd2t - 24*m1**2*M1Lpt**2*m4**2*xd2t + 32*M0Lpt**2*m2**2*m4**2*xd2t - 
     &       56*m1**2*m2**2*m4**2*xd2t - 32*M1Lpt**2*m2**2*m4**2*xd2t - 32*m2**4*m4**2*xd2t - 24*m1**2*M3Lpt**2*m4**2*xd2t - 32*m2**2*M3Lpt**2*m4**2*xd2t - 8*M0Lpt**2*m1**2*sqrtS**2*xd2t + 
     &       8*m1**4*sqrtS**2*xd2t + 8*m1**2*M1Lpt**2*sqrtS**2*xd2t - 4*M0Lpt**2*m2**2*sqrtS**2*xd2t + 52*m1**2*m2**2*sqrtS**2*xd2t + 4*M1Lpt**2*m2**2*sqrtS**2*xd2t - 4*m2**4*sqrtS**2*xd2t - 
     &       14*M0Lpt**2*m3**2*sqrtS**2*xd2t + 86*m1**2*m3**2*sqrtS**2*xd2t + 14*M1Lpt**2*m3**2*sqrtS**2*xd2t + 54*m2**2*m3**2*sqrtS**2*xd2t + 8*m1**2*M3Lpt**2*sqrtS**2*xd2t + 
     &       4*m2**2*M3Lpt**2*sqrtS**2*xd2t + 14*m3**2*M3Lpt**2*sqrtS**2*xd2t - 14*M0Lpt**2*m4**2*sqrtS**2*xd2t + 62*m1**2*m4**2*sqrtS**2*xd2t + 14*M1Lpt**2*m4**2*sqrtS**2*xd2t + 
     &       78*m2**2*m4**2*sqrtS**2*xd2t + 14*M3Lpt**2*m4**2*sqrtS**2*xd2t + 2*M0Lpt**2*sqrtS**4*xd2t - 18*m1**2*sqrtS**4*xd2t - 2*M1Lpt**2*sqrtS**4*xd2t - 10*m2**2*sqrtS**4*xd2t - 
     &       28*m3**2*sqrtS**4*xd2t - 2*M3Lpt**2*sqrtS**4*xd2t - 28*m4**2*sqrtS**4*xd2t + 4*sqrtS**6*xd2t + 4*M0Lpt**2*m1**2*m3**2*xd33t - 4*m1**4*m3**2*xd33t - 
     &       4*m1**2*M1Lpt**2*m3**2*xd33t + 12*M0Lpt**2*m2**2*m3**2*xd33t - 16*m1**2*m2**2*m3**2*xd33t - 12*M1Lpt**2*m2**2*m3**2*xd33t - 12*m2**4*m3**2*xd33t - 4*m1**2*m3**2*M3Lpt**2*xd33t - 
     &       12*m2**2*m3**2*M3Lpt**2*xd33t - 4*M0Lpt**2*m1**2*m4**2*xd33t + 4*m1**4*m4**2*xd33t + 4*m1**2*M1Lpt**2*m4**2*xd33t - 12*M0Lpt**2*m2**2*m4**2*xd33t + 16*m1**2*m2**2*m4**2*xd33t + 
     &       12*M1Lpt**2*m2**2*m4**2*xd33t + 12*m2**4*m4**2*xd33t + 4*m1**2*M3Lpt**2*m4**2*xd33t + 12*m2**2*M3Lpt**2*m4**2*xd33t - 4*M0Lpt**2*m3**2*sqrtS**2*xd33t + 
     &       12*m1**2*m3**2*sqrtS**2*xd33t + 4*M1Lpt**2*m3**2*sqrtS**2*xd33t + 28*m2**2*m3**2*sqrtS**2*xd33t + 4*m3**2*M3Lpt**2*sqrtS**2*xd33t + 4*M0Lpt**2*m4**2*sqrtS**2*xd33t - 
     &       12*m1**2*m4**2*sqrtS**2*xd33t - 4*M1Lpt**2*m4**2*sqrtS**2*xd33t - 28*m2**2*m4**2*sqrtS**2*xd33t - 4*M3Lpt**2*m4**2*sqrtS**2*xd33t - 8*m3**2*sqrtS**4*xd33t + 
     &       8*m4**2*sqrtS**4*xd33t - 2*M0Lpt**2*m1**4*xd3t + 2*m1**6*xd3t + 2*m1**4*M1Lpt**2*xd3t + 20*M0Lpt**2*m1**2*m2**2*xd3t - 18*m1**4*m2**2*xd3t - 20*m1**2*M1Lpt**2*m2**2*xd3t - 
     &       2*M0Lpt**2*m2**4*xd3t - 18*m1**2*m2**4*xd3t + 2*M1Lpt**2*m2**4*xd3t + 2*m2**6*xd3t + 16*M0Lpt**2*m1**2*m3**2*xd3t - 16*m1**4*m3**2*xd3t - 16*m1**2*M1Lpt**2*m3**2*xd3t + 
     &       24*M0Lpt**2*m2**2*m3**2*xd3t - 40*m1**2*m2**2*m3**2*xd3t - 24*M1Lpt**2*m2**2*m3**2*xd3t - 24*m2**4*m3**2*xd3t + 2*m1**4*M3Lpt**2*xd3t - 20*m1**2*m2**2*M3Lpt**2*xd3t + 
     &       2*m2**4*M3Lpt**2*xd3t - 16*m1**2*m3**2*M3Lpt**2*xd3t - 24*m2**2*m3**2*M3Lpt**2*xd3t + 8*M0Lpt**2*m2**2*m4**2*xd3t - 8*m1**2*m2**2*m4**2*xd3t - 8*M1Lpt**2*m2**2*m4**2*xd3t - 
     &       8*m2**4*m4**2*xd3t - 8*m2**2*M3Lpt**2*m4**2*xd3t + 4*M0Lpt**2*m1**2*sqrtS**2*xd3t - 8*m1**4*sqrtS**2*xd3t - 4*m1**2*M1Lpt**2*sqrtS**2*xd3t + 36*m1**2*m2**2*sqrtS**2*xd3t - 
     &       4*m2**4*sqrtS**2*xd3t - 10*M0Lpt**2*m3**2*sqrtS**2*xd3t + 42*m1**2*m3**2*sqrtS**2*xd3t + 10*M1Lpt**2*m3**2*sqrtS**2*xd3t + 58*m2**2*m3**2*sqrtS**2*xd3t - 
     &       4*m1**2*M3Lpt**2*sqrtS**2*xd3t + 10*m3**2*M3Lpt**2*sqrtS**2*xd3t - 2*M0Lpt**2*m4**2*sqrtS**2*xd3t + 2*m1**2*m4**2*sqrtS**2*xd3t + 2*M1Lpt**2*m4**2*sqrtS**2*xd3t + 
     &       18*m2**2*m4**2*sqrtS**2*xd3t + 2*M3Lpt**2*m4**2*sqrtS**2*xd3t - 2*M0Lpt**2*sqrtS**4*xd3t + 10*m1**2*sqrtS**4*xd3t + 2*M1Lpt**2*sqrtS**4*xd3t + 2*m2**2*sqrtS**4*xd3t - 
     &       20*m3**2*sqrtS**4*xd3t + 2*M3Lpt**2*sqrtS**4*xd3t - 4*m4**2*sqrtS**4*xd3t - 4*sqrtS**6*xd3t + 
     &       tman**2*(5*xc01t + 5*xc02t + xc03t + 12*xc111t + 16*xc11t + 16*xc121t - 4*xc122t + 12*xc123t + 6*xc12t + 2*xc13t + 16*xc21t + 12*xc223t - 2*xc22t + 8*xc23t - 5*M0Lpt**2*xd0t + 
     &          5*m1**2*xd0t + 5*M1Lpt**2*xd0t + 5*m2**2*xd0t + 5*M3Lpt**2*xd0t - 10*sqrtS**2*xd0t - 12*M0Lpt**2*xd12t + 12*m1**2*xd12t + 12*M1Lpt**2*xd12t + 12*m2**2*xd12t + 
     &          12*M3Lpt**2*xd12t - 24*sqrtS**2*xd12t - 6*M0Lpt**2*xd1t + 6*m1**2*xd1t + 6*M1Lpt**2*xd1t + 6*m2**2*xd1t + 6*M3Lpt**2*xd1t - 12*sqrtS**2*xd1t - 12*M0Lpt**2*xd22t + 
     &          12*m1**2*xd22t + 12*M1Lpt**2*xd22t + 12*m2**2*xd22t + 12*M3Lpt**2*xd22t - 24*sqrtS**2*xd22t + 4*M0Lpt**2*xd23t - 4*m1**2*xd23t - 4*M1Lpt**2*xd23t - 4*m2**2*xd23t - 
     &          4*M3Lpt**2*xd23t + 8*sqrtS**2*xd23t - 16*M0Lpt**2*xd2t + 16*m1**2*xd2t + 16*M1Lpt**2*xd2t + 16*m2**2*xd2t + 16*M3Lpt**2*xd2t - 32*sqrtS**2*xd2t + 2*M0Lpt**2*xd3t - 
     &          2*m1**2*xd3t - 2*M1Lpt**2*xd3t - 2*m2**2*xd3t - 2*M3Lpt**2*xd3t + 4*sqrtS**2*xd3t) - 
     &       uman**2*(xc01t + xc02t + xc03t + 4*xc111t + 4*xc11t + 16*xc121t - 12*xc122t + 4*xc123t + 2*xc12t + 2*xc13t + 4*xc21t + 4*xc223t - 6*xc22t + 8*xc23t - M0Lpt**2*xd0t + 
     &          m1**2*xd0t + M1Lpt**2*xd0t + m2**2*xd0t + M3Lpt**2*xd0t - 2*sqrtS**2*xd0t - 4*M0Lpt**2*xd12t + 4*m1**2*xd12t + 4*M1Lpt**2*xd12t + 4*m2**2*xd12t + 4*M3Lpt**2*xd12t - 
     &          8*sqrtS**2*xd12t - 2*M0Lpt**2*xd1t + 2*m1**2*xd1t + 2*M1Lpt**2*xd1t + 2*m2**2*xd1t + 2*M3Lpt**2*xd1t - 4*sqrtS**2*xd1t - 4*M0Lpt**2*xd22t + 4*m1**2*xd22t + 4*M1Lpt**2*xd22t + 
     &          4*m2**2*xd22t + 4*M3Lpt**2*xd22t - 8*sqrtS**2*xd22t + 12*M0Lpt**2*xd23t - 12*m1**2*xd23t - 12*M1Lpt**2*xd23t - 12*m2**2*xd23t - 12*M3Lpt**2*xd23t + 24*sqrtS**2*xd23t - 
     &          4*M0Lpt**2*xd2t + 4*m1**2*xd2t + 4*M1Lpt**2*xd2t + 4*m2**2*xd2t + 4*M3Lpt**2*xd2t - 8*sqrtS**2*xd2t + 6*M0Lpt**2*xd3t - 6*m1**2*xd3t - 6*M1Lpt**2*xd3t - 6*m2**2*xd3t - 
     &          6*M3Lpt**2*xd3t + 12*sqrtS**2*xd3t) + uman*(48*xc001t - 48*xc002t - 48*xc003t + 8*m1**2*xc01t + 6*m2**2*xc01t + 3*m3**2*xc01t + 3*m4**2*xc01t - sqrtS**2*xc01t + 
     &          8*m1**2*xc02t + 6*m2**2*xc02t + 3*m3**2*xc02t + 3*m4**2*xc02t - sqrtS**2*xc02t + 2*m2**2*xc03t + 3*m3**2*xc03t - m4**2*xc03t - sqrtS**2*xc03t + 8*m1**2*xc111t + 
     &          8*m2**2*xc111t + 8*m3**2*xc111t + 8*m4**2*xc111t - 8*sqrtS**2*xc111t + 4*m3**2*xc112t - 4*m4**2*xc112t + 4*sqrtS**2*xc112t - 4*m1**2*xc113t + 4*m2**2*xc113t - 
     &          4*sqrtS**2*xc113t + 16*m1**2*xc11t + 12*m2**2*xc11t + 10*m3**2*xc11t + 10*m4**2*xc11t - 6*sqrtS**2*xc11t + 16*m1**2*xc121t + 16*m2**2*xc121t + 24*m3**2*xc121t + 
     &          8*m4**2*xc121t - 4*m1**2*xc122t - 12*m2**2*xc122t - 12*m3**2*xc122t - 4*m4**2*xc122t + 8*sqrtS**2*xc122t + 4*m1**2*xc123t + 12*m2**2*xc123t + 12*m3**2*xc123t + 
     &          4*m4**2*xc123t - 8*sqrtS**2*xc123t + 4*m1**2*xc12t + 12*m2**2*xc12t + 6*m3**2*xc12t + 2*m4**2*xc12t - 6*sqrtS**2*xc12t + 4*m2**2*xc13t + 8*m3**2*xc13t - 4*m4**2*xc13t + 
     &          16*m1**2*xc21t + 12*m2**2*xc21t + 10*m3**2*xc21t + 10*m4**2*xc21t - 6*sqrtS**2*xc21t + 16*sqrtS**2*xc221t + 4*m1**2*xc222t - 4*m2**2*xc222t - 4*sqrtS**2*xc222t + 
     &          8*m1**2*xc223t + 8*m2**2*xc223t + 8*m3**2*xc223t + 8*m4**2*xc223t - 8*sqrtS**2*xc223t + 4*m1**2*xc22t - 4*m2**2*xc22t - 6*m3**2*xc22t - 2*m4**2*xc22t + 2*sqrtS**2*xc22t + 
     &          8*m1**2*xc23t + 8*m2**2*xc23t + 12*m3**2*xc23t + 4*m4**2*xc23t + 8*sqrtS**2*xc23t + 48*M0Lpt**2*xd00t - 48*m1**2*xd00t - 48*M1Lpt**2*xd00t - 48*m2**2*xd00t - 
     &          48*M3Lpt**2*xd00t + 96*sqrtS**2*xd00t - 8*M0Lpt**2*m1**2*xd0t + 8*m1**4*xd0t + 8*m1**2*M1Lpt**2*xd0t - 6*M0Lpt**2*m2**2*xd0t + 14*m1**2*m2**2*xd0t + 6*M1Lpt**2*m2**2*xd0t + 
     &          6*m2**4*xd0t - 3*M0Lpt**2*m3**2*xd0t + 3*m1**2*m3**2*xd0t + 3*M1Lpt**2*m3**2*xd0t + 3*m2**2*m3**2*xd0t + 8*m1**2*M3Lpt**2*xd0t + 6*m2**2*M3Lpt**2*xd0t + 
     &          3*m3**2*M3Lpt**2*xd0t - 3*M0Lpt**2*m4**2*xd0t + 3*m1**2*m4**2*xd0t + 3*M1Lpt**2*m4**2*xd0t + 3*m2**2*m4**2*xd0t + 3*M3Lpt**2*m4**2*xd0t + M0Lpt**2*sqrtS**2*xd0t - 
     &          17*m1**2*sqrtS**2*xd0t - M1Lpt**2*sqrtS**2*xd0t - 13*m2**2*sqrtS**2*xd0t - 6*m3**2*sqrtS**2*xd0t - M3Lpt**2*sqrtS**2*xd0t - 6*m4**2*sqrtS**2*xd0t + 2*sqrtS**4*xd0t + 
     &          4*M0Lpt**2*m1**2*xd11t - 4*m1**4*xd11t - 4*m1**2*M1Lpt**2*xd11t - 4*M0Lpt**2*m2**2*xd11t + 4*M1Lpt**2*m2**2*xd11t + 4*m2**4*xd11t - 4*m1**2*M3Lpt**2*xd11t + 
     &          4*m2**2*M3Lpt**2*xd11t + 4*M0Lpt**2*sqrtS**2*xd11t + 4*m1**2*sqrtS**2*xd11t - 4*M1Lpt**2*sqrtS**2*xd11t - 12*m2**2*sqrtS**2*xd11t - 4*M3Lpt**2*sqrtS**2*xd11t + 
     &          8*sqrtS**4*xd11t - 4*M0Lpt**2*m1**2*xd12t + 4*m1**4*xd12t + 4*m1**2*M1Lpt**2*xd12t - 12*M0Lpt**2*m2**2*xd12t + 16*m1**2*m2**2*xd12t + 12*M1Lpt**2*m2**2*xd12t + 
     &          12*m2**4*xd12t - 12*M0Lpt**2*m3**2*xd12t + 12*m1**2*m3**2*xd12t + 12*M1Lpt**2*m3**2*xd12t + 12*m2**2*m3**2*xd12t + 4*m1**2*M3Lpt**2*xd12t + 12*m2**2*M3Lpt**2*xd12t + 
     &          12*m3**2*M3Lpt**2*xd12t - 4*M0Lpt**2*m4**2*xd12t + 4*m1**2*m4**2*xd12t + 4*M1Lpt**2*m4**2*xd12t + 4*m2**2*m4**2*xd12t + 4*M3Lpt**2*m4**2*xd12t + 8*M0Lpt**2*sqrtS**2*xd12t - 
     &          16*m1**2*sqrtS**2*xd12t - 8*M1Lpt**2*sqrtS**2*xd12t - 32*m2**2*sqrtS**2*xd12t - 24*m3**2*sqrtS**2*xd12t - 8*M3Lpt**2*sqrtS**2*xd12t - 8*m4**2*sqrtS**2*xd12t + 
     &          16*sqrtS**4*xd12t - 8*M0Lpt**2*sqrtS**2*xd13t + 8*m1**2*sqrtS**2*xd13t + 8*M1Lpt**2*sqrtS**2*xd13t + 8*m2**2*sqrtS**2*xd13t + 8*M3Lpt**2*sqrtS**2*xd13t - 16*sqrtS**4*xd13t - 
     &          4*M0Lpt**2*m1**2*xd1t + 4*m1**4*xd1t + 4*m1**2*M1Lpt**2*xd1t - 12*M0Lpt**2*m2**2*xd1t + 16*m1**2*m2**2*xd1t + 12*M1Lpt**2*m2**2*xd1t + 12*m2**4*xd1t - 6*M0Lpt**2*m3**2*xd1t + 
     &          6*m1**2*m3**2*xd1t + 6*M1Lpt**2*m3**2*xd1t + 6*m2**2*m3**2*xd1t + 4*m1**2*M3Lpt**2*xd1t + 12*m2**2*M3Lpt**2*xd1t + 6*m3**2*M3Lpt**2*xd1t - 2*M0Lpt**2*m4**2*xd1t + 
     &          2*m1**2*m4**2*xd1t + 2*M1Lpt**2*m4**2*xd1t + 2*m2**2*m4**2*xd1t + 2*M3Lpt**2*m4**2*xd1t + 6*M0Lpt**2*sqrtS**2*xd1t - 14*m1**2*sqrtS**2*xd1t - 6*M1Lpt**2*sqrtS**2*xd1t - 
     &          30*m2**2*sqrtS**2*xd1t - 12*m3**2*sqrtS**2*xd1t - 6*M3Lpt**2*sqrtS**2*xd1t - 4*m4**2*sqrtS**2*xd1t + 12*sqrtS**4*xd1t - 8*M0Lpt**2*m1**2*xd22t + 8*m1**4*xd22t + 
     &          8*m1**2*M1Lpt**2*xd22t - 8*M0Lpt**2*m2**2*xd22t + 16*m1**2*m2**2*xd22t + 8*M1Lpt**2*m2**2*xd22t + 8*m2**4*xd22t - 8*M0Lpt**2*m3**2*xd22t + 8*m1**2*m3**2*xd22t + 
     &          8*M1Lpt**2*m3**2*xd22t + 8*m2**2*m3**2*xd22t + 8*m1**2*M3Lpt**2*xd22t + 8*m2**2*M3Lpt**2*xd22t + 8*m3**2*M3Lpt**2*xd22t - 8*M0Lpt**2*m4**2*xd22t + 8*m1**2*m4**2*xd22t + 
     &          8*M1Lpt**2*m4**2*xd22t + 8*m2**2*m4**2*xd22t + 8*M3Lpt**2*m4**2*xd22t + 8*M0Lpt**2*sqrtS**2*xd22t - 24*m1**2*sqrtS**2*xd22t - 8*M1Lpt**2*sqrtS**2*xd22t - 
     &          24*m2**2*sqrtS**2*xd22t - 16*m3**2*sqrtS**2*xd22t - 8*M3Lpt**2*sqrtS**2*xd22t - 16*m4**2*sqrtS**2*xd22t + 16*sqrtS**4*xd22t + 4*M0Lpt**2*m1**2*xd23t - 4*m1**4*xd23t - 
     &          4*m1**2*M1Lpt**2*xd23t + 12*M0Lpt**2*m2**2*xd23t - 16*m1**2*m2**2*xd23t - 12*M1Lpt**2*m2**2*xd23t - 12*m2**4*xd23t + 12*M0Lpt**2*m3**2*xd23t - 12*m1**2*m3**2*xd23t - 
     &          12*M1Lpt**2*m3**2*xd23t - 12*m2**2*m3**2*xd23t - 4*m1**2*M3Lpt**2*xd23t - 12*m2**2*M3Lpt**2*xd23t - 12*m3**2*M3Lpt**2*xd23t + 4*M0Lpt**2*m4**2*xd23t - 4*m1**2*m4**2*xd23t - 
     &          4*M1Lpt**2*m4**2*xd23t - 4*m2**2*m4**2*xd23t - 4*M3Lpt**2*m4**2*xd23t - 8*M0Lpt**2*sqrtS**2*xd23t + 16*m1**2*sqrtS**2*xd23t + 8*M1Lpt**2*sqrtS**2*xd23t + 
     &          32*m2**2*sqrtS**2*xd23t + 24*m3**2*sqrtS**2*xd23t + 8*M3Lpt**2*sqrtS**2*xd23t + 8*m4**2*sqrtS**2*xd23t - 16*sqrtS**4*xd23t - 16*M0Lpt**2*m1**2*xd2t + 16*m1**4*xd2t + 
     &          16*m1**2*M1Lpt**2*xd2t - 12*M0Lpt**2*m2**2*xd2t + 28*m1**2*m2**2*xd2t + 12*M1Lpt**2*m2**2*xd2t + 12*m2**4*xd2t - 10*M0Lpt**2*m3**2*xd2t + 10*m1**2*m3**2*xd2t + 
     &          10*M1Lpt**2*m3**2*xd2t + 10*m2**2*m3**2*xd2t + 16*m1**2*M3Lpt**2*xd2t + 12*m2**2*M3Lpt**2*xd2t + 10*m3**2*M3Lpt**2*xd2t - 10*M0Lpt**2*m4**2*xd2t + 10*m1**2*m4**2*xd2t + 
     &          10*M1Lpt**2*m4**2*xd2t + 10*m2**2*m4**2*xd2t + 10*M3Lpt**2*m4**2*xd2t + 6*M0Lpt**2*sqrtS**2*xd2t - 38*m1**2*sqrtS**2*xd2t - 6*M1Lpt**2*sqrtS**2*xd2t - 
     &          30*m2**2*sqrtS**2*xd2t - 20*m3**2*sqrtS**2*xd2t - 6*M3Lpt**2*sqrtS**2*xd2t - 20*m4**2*sqrtS**2*xd2t + 12*sqrtS**4*xd2t - 4*M0Lpt**2*m1**2*xd33t + 4*m1**4*xd33t + 
     &          4*m1**2*M1Lpt**2*xd33t + 4*M0Lpt**2*m2**2*xd33t - 4*M1Lpt**2*m2**2*xd33t - 4*m2**4*xd33t + 4*m1**2*M3Lpt**2*xd33t - 4*m2**2*M3Lpt**2*xd33t + 4*M0Lpt**2*sqrtS**2*xd33t - 
     &          12*m1**2*sqrtS**2*xd33t - 4*M1Lpt**2*sqrtS**2*xd33t + 4*m2**2*sqrtS**2*xd33t - 4*M3Lpt**2*sqrtS**2*xd33t + 8*sqrtS**4*xd33t - 4*M0Lpt**2*m1**2*xd3t + 4*m1**4*xd3t + 
     &          4*m1**2*M1Lpt**2*xd3t + 4*M0Lpt**2*m2**2*xd3t - 4*M1Lpt**2*m2**2*xd3t - 4*m2**4*xd3t + 6*M0Lpt**2*m3**2*xd3t - 6*m1**2*m3**2*xd3t - 6*M1Lpt**2*m3**2*xd3t - 
     &          6*m2**2*m3**2*xd3t + 4*m1**2*M3Lpt**2*xd3t - 4*m2**2*M3Lpt**2*xd3t - 6*m3**2*M3Lpt**2*xd3t + 2*M0Lpt**2*m4**2*xd3t - 2*m1**2*m4**2*xd3t - 2*M1Lpt**2*m4**2*xd3t - 
     &          2*m2**2*m4**2*xd3t - 2*M3Lpt**2*m4**2*xd3t - 2*M0Lpt**2*sqrtS**2*xd3t - 6*m1**2*sqrtS**2*xd3t + 2*M1Lpt**2*sqrtS**2*xd3t + 10*m2**2*sqrtS**2*xd3t + 12*m3**2*sqrtS**2*xd3t + 
     &          2*M3Lpt**2*sqrtS**2*xd3t + 4*m4**2*sqrtS**2*xd3t - 4*sqrtS**4*xd3t) - 
     &          tman*(48*xc001t - 48*xc002t - 48*xc003t - 2*m1**2*xc01t + 3*m3**2*xc01t + 3*m4**2*xc01t + 3*sqrtS**2*xc01t + 4*uman*xc01t - 2*m1**2*xc02t + 3*m3**2*xc02t + 3*m4**2*xc02t + 
     &          3*sqrtS**2*xc02t + 4*uman*xc02t + 2*m1**2*xc03t + 3*m3**2*xc03t - m4**2*xc03t - sqrtS**2*xc03t + 8*m1**2*xc111t + 8*m2**2*xc111t + 8*m3**2*xc111t + 8*m4**2*xc111t - 
     &          8*sqrtS**2*xc111t + 8*uman*xc111t + 4*m3**2*xc112t - 4*m4**2*xc112t + 4*sqrtS**2*xc112t - 4*m1**2*xc113t + 4*m2**2*xc113t - 4*sqrtS**2*xc113t + 4*m1**2*xc11t + 
     &          8*m2**2*xc11t + 10*m3**2*xc11t + 10*m4**2*xc11t - 2*sqrtS**2*xc11t + 12*uman*xc11t + 16*m1**2*xc121t + 16*m2**2*xc121t + 24*m3**2*xc121t + 8*m4**2*xc121t - 4*m1**2*xc122t - 
     &          12*m2**2*xc122t - 12*m3**2*xc122t - 4*m4**2*xc122t + 8*sqrtS**2*xc122t + 8*uman*xc122t + 4*m1**2*xc123t + 12*m2**2*xc123t + 12*m3**2*xc123t + 4*m4**2*xc123t - 
     &          8*sqrtS**2*xc123t + 8*uman*xc123t - 4*m1**2*xc12t + 4*m2**2*xc12t + 6*m3**2*xc12t + 2*m4**2*xc12t - 2*sqrtS**2*xc12t + 4*uman*xc12t + 4*m1**2*xc13t + 8*m3**2*xc13t - 
     &          4*m4**2*xc13t + 4*m1**2*xc21t + 8*m2**2*xc21t + 10*m3**2*xc21t + 10*m4**2*xc21t - 2*sqrtS**2*xc21t + 12*uman*xc21t + 16*sqrtS**2*xc221t + 4*m1**2*xc222t - 4*m2**2*xc222t - 
     &          4*sqrtS**2*xc222t + 8*m1**2*xc223t + 8*m2**2*xc223t + 8*m3**2*xc223t + 8*m4**2*xc223t - 8*sqrtS**2*xc223t + 8*uman*xc223t - 4*m1**2*xc22t - 12*m2**2*xc22t - 6*m3**2*xc22t - 
     &          2*m4**2*xc22t + 6*sqrtS**2*xc22t + 4*uman*xc22t + 8*m1**2*xc23t + 8*m2**2*xc23t + 12*m3**2*xc23t + 4*m4**2*xc23t + 8*sqrtS**2*xc23t + 48*M0Lpt**2*xd00t - 48*m1**2*xd00t - 
     &          48*M1Lpt**2*xd00t - 48*m2**2*xd00t - 48*M3Lpt**2*xd00t + 96*sqrtS**2*xd00t + 2*M0Lpt**2*m1**2*xd0t - 2*m1**4*xd0t - 2*m1**2*M1Lpt**2*xd0t - 2*m1**2*m2**2*xd0t - 
     &          3*M0Lpt**2*m3**2*xd0t + 3*m1**2*m3**2*xd0t + 3*M1Lpt**2*m3**2*xd0t + 3*m2**2*m3**2*xd0t - 2*m1**2*M3Lpt**2*xd0t + 3*m3**2*M3Lpt**2*xd0t - 3*M0Lpt**2*m4**2*xd0t + 
     &          3*m1**2*m4**2*xd0t + 3*M1Lpt**2*m4**2*xd0t + 3*m2**2*m4**2*xd0t + 3*M3Lpt**2*m4**2*xd0t - 3*M0Lpt**2*sqrtS**2*xd0t + 7*m1**2*sqrtS**2*xd0t + 3*M1Lpt**2*sqrtS**2*xd0t + 
     &          3*m2**2*sqrtS**2*xd0t - 6*m3**2*sqrtS**2*xd0t + 3*M3Lpt**2*sqrtS**2*xd0t - 6*m4**2*sqrtS**2*xd0t - 6*sqrtS**4*xd0t - 4*M0Lpt**2*uman*xd0t + 4*m1**2*uman*xd0t + 
     &          4*M1Lpt**2*uman*xd0t + 4*m2**2*uman*xd0t + 4*M3Lpt**2*uman*xd0t - 8*sqrtS**2*uman*xd0t + 4*M0Lpt**2*m1**2*xd11t - 4*m1**4*xd11t - 4*m1**2*M1Lpt**2*xd11t - 
     &          4*M0Lpt**2*m2**2*xd11t + 4*M1Lpt**2*m2**2*xd11t + 4*m2**4*xd11t - 4*m1**2*M3Lpt**2*xd11t + 4*m2**2*M3Lpt**2*xd11t + 4*M0Lpt**2*sqrtS**2*xd11t + 4*m1**2*sqrtS**2*xd11t - 
     &          4*M1Lpt**2*sqrtS**2*xd11t - 12*m2**2*sqrtS**2*xd11t - 4*M3Lpt**2*sqrtS**2*xd11t + 8*sqrtS**4*xd11t - 4*M0Lpt**2*m1**2*xd12t + 4*m1**4*xd12t + 4*m1**2*M1Lpt**2*xd12t - 
     &          12*M0Lpt**2*m2**2*xd12t + 16*m1**2*m2**2*xd12t + 12*M1Lpt**2*m2**2*xd12t + 12*m2**4*xd12t - 12*M0Lpt**2*m3**2*xd12t + 12*m1**2*m3**2*xd12t + 12*M1Lpt**2*m3**2*xd12t + 
     &          12*m2**2*m3**2*xd12t + 4*m1**2*M3Lpt**2*xd12t + 12*m2**2*M3Lpt**2*xd12t + 12*m3**2*M3Lpt**2*xd12t - 4*M0Lpt**2*m4**2*xd12t + 4*m1**2*m4**2*xd12t + 4*M1Lpt**2*m4**2*xd12t + 
     &          4*m2**2*m4**2*xd12t + 4*M3Lpt**2*m4**2*xd12t + 8*M0Lpt**2*sqrtS**2*xd12t - 16*m1**2*sqrtS**2*xd12t - 8*M1Lpt**2*sqrtS**2*xd12t - 32*m2**2*sqrtS**2*xd12t - 
     &          24*m3**2*sqrtS**2*xd12t - 8*M3Lpt**2*sqrtS**2*xd12t - 8*m4**2*sqrtS**2*xd12t + 16*sqrtS**4*xd12t - 8*M0Lpt**2*uman*xd12t + 8*m1**2*uman*xd12t + 8*M1Lpt**2*uman*xd12t + 
     &          8*m2**2*uman*xd12t + 8*M3Lpt**2*uman*xd12t - 16*sqrtS**2*uman*xd12t - 8*M0Lpt**2*sqrtS**2*xd13t + 8*m1**2*sqrtS**2*xd13t + 8*M1Lpt**2*sqrtS**2*xd13t + 
     &          8*m2**2*sqrtS**2*xd13t + 8*M3Lpt**2*sqrtS**2*xd13t - 16*sqrtS**4*xd13t + 4*M0Lpt**2*m1**2*xd1t - 4*m1**4*xd1t - 4*m1**2*M1Lpt**2*xd1t - 4*M0Lpt**2*m2**2*xd1t + 
     &          4*M1Lpt**2*m2**2*xd1t + 4*m2**4*xd1t - 6*M0Lpt**2*m3**2*xd1t + 6*m1**2*m3**2*xd1t + 6*M1Lpt**2*m3**2*xd1t + 6*m2**2*m3**2*xd1t - 4*m1**2*M3Lpt**2*xd1t + 
     &          4*m2**2*M3Lpt**2*xd1t + 6*m3**2*M3Lpt**2*xd1t - 2*M0Lpt**2*m4**2*xd1t + 2*m1**2*m4**2*xd1t + 2*M1Lpt**2*m4**2*xd1t + 2*m2**2*m4**2*xd1t + 2*M3Lpt**2*m4**2*xd1t + 
     &          2*M0Lpt**2*sqrtS**2*xd1t + 6*m1**2*sqrtS**2*xd1t - 2*M1Lpt**2*sqrtS**2*xd1t - 10*m2**2*sqrtS**2*xd1t - 12*m3**2*sqrtS**2*xd1t - 2*M3Lpt**2*sqrtS**2*xd1t - 
     &          4*m4**2*sqrtS**2*xd1t + 4*sqrtS**4*xd1t - 4*M0Lpt**2*uman*xd1t + 4*m1**2*uman*xd1t + 4*M1Lpt**2*uman*xd1t + 4*m2**2*uman*xd1t + 4*M3Lpt**2*uman*xd1t - 8*sqrtS**2*uman*xd1t - 
     &          8*M0Lpt**2*m1**2*xd22t + 8*m1**4*xd22t + 8*m1**2*M1Lpt**2*xd22t - 8*M0Lpt**2*m2**2*xd22t + 16*m1**2*m2**2*xd22t + 8*M1Lpt**2*m2**2*xd22t + 8*m2**4*xd22t - 
     &          8*M0Lpt**2*m3**2*xd22t + 8*m1**2*m3**2*xd22t + 8*M1Lpt**2*m3**2*xd22t + 8*m2**2*m3**2*xd22t + 8*m1**2*M3Lpt**2*xd22t + 8*m2**2*M3Lpt**2*xd22t + 8*m3**2*M3Lpt**2*xd22t - 
     &          8*M0Lpt**2*m4**2*xd22t + 8*m1**2*m4**2*xd22t + 8*M1Lpt**2*m4**2*xd22t + 8*m2**2*m4**2*xd22t + 8*M3Lpt**2*m4**2*xd22t + 8*M0Lpt**2*sqrtS**2*xd22t - 24*m1**2*sqrtS**2*xd22t - 
     &          8*M1Lpt**2*sqrtS**2*xd22t - 24*m2**2*sqrtS**2*xd22t - 16*m3**2*sqrtS**2*xd22t - 8*M3Lpt**2*sqrtS**2*xd22t - 16*m4**2*sqrtS**2*xd22t + 16*sqrtS**4*xd22t - 
     &          8*M0Lpt**2*uman*xd22t + 8*m1**2*uman*xd22t + 8*M1Lpt**2*uman*xd22t + 8*m2**2*uman*xd22t + 8*M3Lpt**2*uman*xd22t - 16*sqrtS**2*uman*xd22t + 4*M0Lpt**2*m1**2*xd23t - 
     &          4*m1**4*xd23t - 4*m1**2*M1Lpt**2*xd23t + 12*M0Lpt**2*m2**2*xd23t - 16*m1**2*m2**2*xd23t - 12*M1Lpt**2*m2**2*xd23t - 12*m2**4*xd23t + 12*M0Lpt**2*m3**2*xd23t - 
     &          12*m1**2*m3**2*xd23t - 12*M1Lpt**2*m3**2*xd23t - 12*m2**2*m3**2*xd23t - 4*m1**2*M3Lpt**2*xd23t - 12*m2**2*M3Lpt**2*xd23t - 12*m3**2*M3Lpt**2*xd23t + 4*M0Lpt**2*m4**2*xd23t - 
     &          4*m1**2*m4**2*xd23t - 4*M1Lpt**2*m4**2*xd23t - 4*m2**2*m4**2*xd23t - 4*M3Lpt**2*m4**2*xd23t - 8*M0Lpt**2*sqrtS**2*xd23t + 16*m1**2*sqrtS**2*xd23t + 
     &          8*M1Lpt**2*sqrtS**2*xd23t + 32*m2**2*sqrtS**2*xd23t + 24*m3**2*sqrtS**2*xd23t + 8*M3Lpt**2*sqrtS**2*xd23t + 8*m4**2*sqrtS**2*xd23t - 16*sqrtS**4*xd23t - 
     &          8*M0Lpt**2*uman*xd23t + 8*m1**2*uman*xd23t + 8*M1Lpt**2*uman*xd23t + 8*m2**2*uman*xd23t + 8*M3Lpt**2*uman*xd23t - 16*sqrtS**2*uman*xd23t - 4*M0Lpt**2*m1**2*xd2t + 
     &          4*m1**4*xd2t + 4*m1**2*M1Lpt**2*xd2t - 8*M0Lpt**2*m2**2*xd2t + 12*m1**2*m2**2*xd2t + 8*M1Lpt**2*m2**2*xd2t + 8*m2**4*xd2t - 10*M0Lpt**2*m3**2*xd2t + 10*m1**2*m3**2*xd2t + 
     &          10*M1Lpt**2*m3**2*xd2t + 10*m2**2*m3**2*xd2t + 4*m1**2*M3Lpt**2*xd2t + 8*m2**2*M3Lpt**2*xd2t + 10*m3**2*M3Lpt**2*xd2t - 10*M0Lpt**2*m4**2*xd2t + 10*m1**2*m4**2*xd2t + 
     &          10*M1Lpt**2*m4**2*xd2t + 10*m2**2*m4**2*xd2t + 10*M3Lpt**2*m4**2*xd2t + 2*M0Lpt**2*sqrtS**2*xd2t - 10*m1**2*sqrtS**2*xd2t - 2*M1Lpt**2*sqrtS**2*xd2t - 
     &          18*m2**2*sqrtS**2*xd2t - 20*m3**2*sqrtS**2*xd2t - 2*M3Lpt**2*sqrtS**2*xd2t - 20*m4**2*sqrtS**2*xd2t + 4*sqrtS**4*xd2t - 12*M0Lpt**2*uman*xd2t + 12*m1**2*uman*xd2t + 
     &          12*M1Lpt**2*uman*xd2t + 12*m2**2*uman*xd2t + 12*M3Lpt**2*uman*xd2t - 24*sqrtS**2*uman*xd2t - 4*M0Lpt**2*m1**2*xd33t + 4*m1**4*xd33t + 4*m1**2*M1Lpt**2*xd33t + 
     &          4*M0Lpt**2*m2**2*xd33t - 4*M1Lpt**2*m2**2*xd33t - 4*m2**4*xd33t + 4*m1**2*M3Lpt**2*xd33t - 4*m2**2*M3Lpt**2*xd33t + 4*M0Lpt**2*sqrtS**2*xd33t - 12*m1**2*sqrtS**2*xd33t - 
     &          4*M1Lpt**2*sqrtS**2*xd33t + 4*m2**2*sqrtS**2*xd33t - 4*M3Lpt**2*sqrtS**2*xd33t + 8*sqrtS**4*xd33t + 4*M0Lpt**2*m1**2*xd3t - 4*m1**4*xd3t - 4*m1**2*M1Lpt**2*xd3t + 
     &          12*M0Lpt**2*m2**2*xd3t - 16*m1**2*m2**2*xd3t - 12*M1Lpt**2*m2**2*xd3t - 12*m2**4*xd3t + 6*M0Lpt**2*m3**2*xd3t - 6*m1**2*m3**2*xd3t - 6*M1Lpt**2*m3**2*xd3t - 
     &          6*m2**2*m3**2*xd3t - 4*m1**2*M3Lpt**2*xd3t - 12*m2**2*M3Lpt**2*xd3t - 6*m3**2*M3Lpt**2*xd3t + 2*M0Lpt**2*m4**2*xd3t - 2*m1**2*m4**2*xd3t - 2*M1Lpt**2*m4**2*xd3t - 
     &          2*m2**2*m4**2*xd3t - 2*M3Lpt**2*m4**2*xd3t - 6*M0Lpt**2*sqrtS**2*xd3t + 14*m1**2*sqrtS**2*xd3t + 6*M1Lpt**2*sqrtS**2*xd3t + 30*m2**2*sqrtS**2*xd3t + 12*m3**2*sqrtS**2*xd3t + 
     &          6*M3Lpt**2*sqrtS**2*xd3t + 4*m4**2*sqrtS**2*xd3t - 12*sqrtS**4*xd3t - 4*M0Lpt**2*uman*xd3t + 4*m1**2*uman*xd3t + 4*M1Lpt**2*uman*xd3t + 4*m2**2*uman*xd3t + 
     &          4*M3Lpt**2*uman*xd3t - 8*sqrtS**2*uman*xd3t)))/(64.*pi**2*propac)  )

	Box_UV = Box_UV + 2*colorfactorSV*dreal(  -(Coupl_Ac*Coupl_Bc*g0u*g1u*g2u*g3u*(-(m1**4*xc01u) + 20*m1**2*m2**2*xc01u - 3*m2**4*xc01u + 9*m1**2*m3**2*xc01u + 15*m2**2*m3**2*xc01u + 17*m1**2*m4**2*xc01u + 7*m2**2*m4**2*xc01u + 
     &       m1**2*sqrtS**2*xc01u + 3*m2**2*sqrtS**2*xc01u - 6*m3**2*sqrtS**2*xc01u - 6*m4**2*sqrtS**2*xc01u - 2*sqrtS**4*xc01u - m1**4*xc02u + 20*m1**2*m2**2*xc02u - 3*m2**4*xc02u + 
     &       9*m1**2*m3**2*xc02u + 15*m2**2*m3**2*xc02u + 17*m1**2*m4**2*xc02u + 7*m2**2*m4**2*xc02u + m1**2*sqrtS**2*xc02u + 3*m2**2*sqrtS**2*xc02u - 6*m3**2*sqrtS**2*xc02u - 
     &       6*m4**2*sqrtS**2*xc02u - 2*sqrtS**4*xc02u - m1**4*xc03u + m2**4*xc03u + 3*m1**2*m3**2*xc03u - 3*m2**2*m3**2*xc03u + 3*m1**2*m4**2*xc03u - 3*m2**2*m4**2*xc03u - 
     &       m1**2*sqrtS**2*xc03u + m2**2*sqrtS**2*xc03u + 16*m1**2*m3**2*xc111u + 16*m2**2*m3**2*xc111u + 16*m1**2*m4**2*xc111u + 16*m2**2*m4**2*xc111u - 8*m1**2*sqrtS**2*xc111u - 
     &       8*m2**2*sqrtS**2*xc111u - 8*m3**2*sqrtS**2*xc111u - 8*m4**2*sqrtS**2*xc111u + 4*sqrtS**4*xc111u + 4*m1**2*m3**2*xc112u - 4*m2**2*m3**2*xc112u + 12*m1**2*m4**2*xc112u - 
     &       12*m2**2*m4**2*xc112u - 4*m1**2*sqrtS**2*xc112u + 4*m2**2*sqrtS**2*xc112u + 12*m1**2*m3**2*xc113u + 4*m2**2*m3**2*xc113u - 12*m1**2*m4**2*xc113u - 4*m2**2*m4**2*xc113u - 
     &       4*m3**2*sqrtS**2*xc113u + 4*m4**2*sqrtS**2*xc113u + 20*m1**2*m2**2*xc11u - 4*m2**4*xc11u + 24*m1**2*m3**2*xc11u + 32*m2**2*m3**2*xc11u + 36*m1**2*m4**2*xc11u + 
     &       20*m2**2*m4**2*xc11u - 8*m1**2*sqrtS**2*xc11u - 4*m2**2*sqrtS**2*xc11u - 14*m3**2*sqrtS**2*xc11u - 14*m4**2*sqrtS**2*xc11u + 2*sqrtS**4*xc11u + 32*m1**2*m4**2*xc121u - 
     &       32*m2**2*m4**2*xc121u - 8*m1**2*sqrtS**2*xc121u + 8*m2**2*sqrtS**2*xc121u + 16*m2**2*m3**2*xc122u + 16*m1**2*m4**2*xc122u + 32*m2**2*m4**2*xc122u - 4*m1**2*sqrtS**2*xc122u - 
     &       12*m2**2*sqrtS**2*xc122u - 4*m3**2*sqrtS**2*xc122u - 12*m4**2*sqrtS**2*xc122u + 4*sqrtS**4*xc122u + 32*m1**2*m3**2*xc123u + 16*m2**2*m3**2*xc123u + 16*m1**2*m4**2*xc123u - 
     &       12*m1**2*sqrtS**2*xc123u - 4*m2**2*sqrtS**2*xc123u - 12*m3**2*sqrtS**2*xc123u - 4*m4**2*sqrtS**2*xc123u + 4*sqrtS**4*xc123u - 2*m1**4*xc12u + 20*m1**2*m2**2*xc12u - 
     &       2*m2**4*xc12u + 24*m1**2*m3**2*xc12u + 16*m2**2*m3**2*xc12u + 8*m1**2*m4**2*xc12u + 4*m2**2*sqrtS**2*xc12u - 10*m3**2*sqrtS**2*xc12u - 2*m4**2*sqrtS**2*xc12u - 
     &       2*sqrtS**4*xc12u - 2*m1**4*xc13u + 2*m2**4*xc13u + 8*m1**2*m3**2*xc13u - 8*m2**2*m3**2*xc13u + 12*m1**2*m4**2*xc13u - 12*m2**2*m4**2*xc13u - 4*m1**2*sqrtS**2*xc13u + 
     &       4*m2**2*sqrtS**2*xc13u + 20*m1**2*m2**2*xc21u - 4*m2**4*xc21u + 24*m1**2*m3**2*xc21u + 32*m2**2*m3**2*xc21u + 36*m1**2*m4**2*xc21u + 20*m2**2*m4**2*xc21u - 
     &       8*m1**2*sqrtS**2*xc21u - 4*m2**2*sqrtS**2*xc21u - 14*m3**2*sqrtS**2*xc21u - 14*m4**2*sqrtS**2*xc21u + 2*sqrtS**4*xc21u - 16*m1**2*m3**2*xc221u + 16*m2**2*m3**2*xc221u + 
     &       16*m1**2*m4**2*xc221u - 16*m2**2*m4**2*xc221u - 4*m1**2*m3**2*xc222u - 12*m2**2*m3**2*xc222u + 4*m1**2*m4**2*xc222u + 12*m2**2*m4**2*xc222u + 4*m3**2*sqrtS**2*xc222u - 
     &       4*m4**2*sqrtS**2*xc222u + 16*m1**2*m3**2*xc223u + 16*m2**2*m3**2*xc223u + 16*m1**2*m4**2*xc223u + 16*m2**2*m4**2*xc223u - 8*m1**2*sqrtS**2*xc223u - 8*m2**2*sqrtS**2*xc223u - 
     &       8*m3**2*sqrtS**2*xc223u - 8*m4**2*sqrtS**2*xc223u + 4*sqrtS**4*xc223u - 2*m1**4*xc22u + 20*m1**2*m2**2*xc22u - 2*m2**4*xc22u + 8*m2**2*m3**2*xc22u + 16*m1**2*m4**2*xc22u + 
     &       24*m2**2*m4**2*xc22u + 4*m1**2*sqrtS**2*xc22u - 2*m3**2*sqrtS**2*xc22u - 10*m4**2*sqrtS**2*xc22u - 2*sqrtS**4*xc22u - 8*m1**2*m3**2*xc23u + 8*m2**2*m3**2*xc23u + 
     &       24*m1**2*m4**2*xc23u - 24*m2**2*m4**2*xc23u - 4*m1**2*sqrtS**2*xc23u + 4*m2**2*sqrtS**2*xc23u + M0Lpu**2*m1**4*xd0u - m1**6*xd0u - m1**4*M1Lpu**2*xd0u - 
     &       20*M0Lpu**2*m1**2*m2**2*xd0u + 19*m1**4*m2**2*xd0u + 20*m1**2*M1Lpu**2*m2**2*xd0u + 3*M0Lpu**2*m2**4*xd0u + 17*m1**2*m2**4*xd0u - 3*M1Lpu**2*m2**4*xd0u - 3*m2**6*xd0u - 
     &       9*M0Lpu**2*m1**2*m3**2*xd0u + 9*m1**4*m3**2*xd0u + 9*m1**2*M1Lpu**2*m3**2*xd0u - 15*M0Lpu**2*m2**2*m3**2*xd0u + 24*m1**2*m2**2*m3**2*xd0u + 15*M1Lpu**2*m2**2*m3**2*xd0u + 
     &       15*m2**4*m3**2*xd0u - m1**4*M3Lpu**2*xd0u + 20*m1**2*m2**2*M3Lpu**2*xd0u - 3*m2**4*M3Lpu**2*xd0u + 9*m1**2*m3**2*M3Lpu**2*xd0u + 15*m2**2*m3**2*M3Lpu**2*xd0u - 
     &       17*M0Lpu**2*m1**2*m4**2*xd0u + 17*m1**4*m4**2*xd0u + 17*m1**2*M1Lpu**2*m4**2*xd0u - 7*M0Lpu**2*m2**2*m4**2*xd0u + 24*m1**2*m2**2*m4**2*xd0u + 7*M1Lpu**2*m2**2*m4**2*xd0u + 
     &       7*m2**4*m4**2*xd0u + 17*m1**2*M3Lpu**2*m4**2*xd0u + 7*m2**2*M3Lpu**2*m4**2*xd0u - M0Lpu**2*m1**2*sqrtS**2*xd0u + 3*m1**4*sqrtS**2*xd0u + m1**2*M1Lpu**2*sqrtS**2*xd0u - 
     &       3*M0Lpu**2*m2**2*sqrtS**2*xd0u - 36*m1**2*m2**2*sqrtS**2*xd0u + 3*M1Lpu**2*m2**2*sqrtS**2*xd0u + 9*m2**4*sqrtS**2*xd0u + 6*M0Lpu**2*m3**2*sqrtS**2*xd0u - 
     &       24*m1**2*m3**2*sqrtS**2*xd0u - 6*M1Lpu**2*m3**2*sqrtS**2*xd0u - 36*m2**2*m3**2*sqrtS**2*xd0u + m1**2*M3Lpu**2*sqrtS**2*xd0u + 3*m2**2*M3Lpu**2*sqrtS**2*xd0u - 
     &       6*m3**2*M3Lpu**2*sqrtS**2*xd0u + 6*M0Lpu**2*m4**2*sqrtS**2*xd0u - 40*m1**2*m4**2*sqrtS**2*xd0u - 6*M1Lpu**2*m4**2*sqrtS**2*xd0u - 20*m2**2*m4**2*sqrtS**2*xd0u - 
     &       6*M3Lpu**2*m4**2*sqrtS**2*xd0u + 2*M0Lpu**2*sqrtS**4*xd0u - 4*m1**2*sqrtS**4*xd0u - 2*M1Lpu**2*sqrtS**4*xd0u - 8*m2**2*sqrtS**4*xd0u + 12*m3**2*sqrtS**4*xd0u - 
     &       2*M3Lpu**2*sqrtS**4*xd0u + 12*m4**2*sqrtS**4*xd0u + 4*sqrtS**6*xd0u - 12*M0Lpu**2*m1**2*m3**2*xd11u + 12*m1**4*m3**2*xd11u + 12*m1**2*M1Lpu**2*m3**2*xd11u - 
     &       4*M0Lpu**2*m2**2*m3**2*xd11u + 16*m1**2*m2**2*m3**2*xd11u + 4*M1Lpu**2*m2**2*m3**2*xd11u + 4*m2**4*m3**2*xd11u + 12*m1**2*m3**2*M3Lpu**2*xd11u + 4*m2**2*m3**2*M3Lpu**2*xd11u + 
     &       12*M0Lpu**2*m1**2*m4**2*xd11u - 12*m1**4*m4**2*xd11u - 12*m1**2*M1Lpu**2*m4**2*xd11u + 4*M0Lpu**2*m2**2*m4**2*xd11u - 16*m1**2*m2**2*m4**2*xd11u - 4*M1Lpu**2*m2**2*m4**2*xd11u - 
     &       4*m2**4*m4**2*xd11u - 12*m1**2*M3Lpu**2*m4**2*xd11u - 4*m2**2*M3Lpu**2*m4**2*xd11u + 4*M0Lpu**2*m3**2*sqrtS**2*xd11u - 28*m1**2*m3**2*sqrtS**2*xd11u - 
     &       4*M1Lpu**2*m3**2*sqrtS**2*xd11u - 12*m2**2*m3**2*sqrtS**2*xd11u - 4*m3**2*M3Lpu**2*sqrtS**2*xd11u - 4*M0Lpu**2*m4**2*sqrtS**2*xd11u + 28*m1**2*m4**2*sqrtS**2*xd11u + 
     &       4*M1Lpu**2*m4**2*sqrtS**2*xd11u + 12*m2**2*m4**2*sqrtS**2*xd11u + 4*M3Lpu**2*m4**2*sqrtS**2*xd11u + 8*m3**2*sqrtS**4*xd11u - 8*m4**2*sqrtS**4*xd11u - 
     &       32*M0Lpu**2*m1**2*m3**2*xd12u + 32*m1**4*m3**2*xd12u + 32*m1**2*M1Lpu**2*m3**2*xd12u - 16*M0Lpu**2*m2**2*m3**2*xd12u + 48*m1**2*m2**2*m3**2*xd12u + 
     &       16*M1Lpu**2*m2**2*m3**2*xd12u + 16*m2**4*m3**2*xd12u + 32*m1**2*m3**2*M3Lpu**2*xd12u + 16*m2**2*m3**2*M3Lpu**2*xd12u - 16*M0Lpu**2*m1**2*m4**2*xd12u + 16*m1**4*m4**2*xd12u + 
     &       16*m1**2*M1Lpu**2*m4**2*xd12u + 16*m1**2*m2**2*m4**2*xd12u + 16*m1**2*M3Lpu**2*m4**2*xd12u + 12*M0Lpu**2*m1**2*sqrtS**2*xd12u - 12*m1**4*sqrtS**2*xd12u - 
     &       12*m1**2*M1Lpu**2*sqrtS**2*xd12u + 4*M0Lpu**2*m2**2*sqrtS**2*xd12u - 16*m1**2*m2**2*sqrtS**2*xd12u - 4*M1Lpu**2*m2**2*sqrtS**2*xd12u - 4*m2**4*sqrtS**2*xd12u + 
     &       12*M0Lpu**2*m3**2*sqrtS**2*xd12u - 76*m1**2*m3**2*sqrtS**2*xd12u - 12*M1Lpu**2*m3**2*sqrtS**2*xd12u - 44*m2**2*m3**2*sqrtS**2*xd12u - 12*m1**2*M3Lpu**2*sqrtS**2*xd12u - 
     &       4*m2**2*M3Lpu**2*sqrtS**2*xd12u - 12*m3**2*M3Lpu**2*sqrtS**2*xd12u + 4*M0Lpu**2*m4**2*sqrtS**2*xd12u - 36*m1**2*m4**2*sqrtS**2*xd12u - 4*M1Lpu**2*m4**2*sqrtS**2*xd12u - 
     &       4*m2**2*m4**2*sqrtS**2*xd12u - 4*M3Lpu**2*m4**2*sqrtS**2*xd12u - 4*M0Lpu**2*sqrtS**4*xd12u + 28*m1**2*sqrtS**4*xd12u + 4*M1Lpu**2*sqrtS**4*xd12u + 12*m2**2*sqrtS**4*xd12u + 
     &       24*m3**2*sqrtS**4*xd12u + 4*M3Lpu**2*sqrtS**4*xd12u + 8*m4**2*sqrtS**4*xd12u - 8*sqrtS**6*xd12u + 8*M0Lpu**2*m1**2*m3**2*xd13u - 8*m1**4*m3**2*xd13u - 
     &       8*m1**2*M1Lpu**2*m3**2*xd13u - 8*M0Lpu**2*m2**2*m3**2*xd13u + 8*M1Lpu**2*m2**2*m3**2*xd13u + 8*m2**4*m3**2*xd13u - 8*m1**2*m3**2*M3Lpu**2*xd13u + 8*m2**2*m3**2*M3Lpu**2*xd13u - 
     &       8*M0Lpu**2*m1**2*m4**2*xd13u + 8*m1**4*m4**2*xd13u + 8*m1**2*M1Lpu**2*m4**2*xd13u + 8*M0Lpu**2*m2**2*m4**2*xd13u - 8*M1Lpu**2*m2**2*m4**2*xd13u - 8*m2**4*m4**2*xd13u + 
     &       8*m1**2*M3Lpu**2*m4**2*xd13u - 8*m2**2*M3Lpu**2*m4**2*xd13u + 16*m1**2*m3**2*sqrtS**2*xd13u - 16*m2**2*m3**2*sqrtS**2*xd13u - 16*m1**2*m4**2*sqrtS**2*xd13u + 
     &       16*m2**2*m4**2*sqrtS**2*xd13u + 2*M0Lpu**2*m1**4*xd1u - 2*m1**6*xd1u - 2*m1**4*M1Lpu**2*xd1u - 20*M0Lpu**2*m1**2*m2**2*xd1u + 18*m1**4*m2**2*xd1u + 
     &       20*m1**2*M1Lpu**2*m2**2*xd1u + 2*M0Lpu**2*m2**4*xd1u + 18*m1**2*m2**4*xd1u - 2*M1Lpu**2*m2**4*xd1u - 2*m2**6*xd1u - 24*M0Lpu**2*m1**2*m3**2*xd1u + 24*m1**4*m3**2*xd1u + 
     &       24*m1**2*M1Lpu**2*m3**2*xd1u - 16*M0Lpu**2*m2**2*m3**2*xd1u + 40*m1**2*m2**2*m3**2*xd1u + 16*M1Lpu**2*m2**2*m3**2*xd1u + 16*m2**4*m3**2*xd1u - 2*m1**4*M3Lpu**2*xd1u + 
     &       20*m1**2*m2**2*M3Lpu**2*xd1u - 2*m2**4*M3Lpu**2*xd1u + 24*m1**2*m3**2*M3Lpu**2*xd1u + 16*m2**2*m3**2*M3Lpu**2*xd1u - 8*M0Lpu**2*m1**2*m4**2*xd1u + 8*m1**4*m4**2*xd1u + 
     &       8*m1**2*M1Lpu**2*m4**2*xd1u + 8*m1**2*m2**2*m4**2*xd1u + 8*m1**2*M3Lpu**2*m4**2*xd1u + 4*m1**4*sqrtS**2*xd1u - 4*M0Lpu**2*m2**2*sqrtS**2*xd1u - 36*m1**2*m2**2*sqrtS**2*xd1u + 
     &       4*M1Lpu**2*m2**2*sqrtS**2*xd1u + 8*m2**4*sqrtS**2*xd1u + 10*M0Lpu**2*m3**2*sqrtS**2*xd1u - 58*m1**2*m3**2*sqrtS**2*xd1u - 10*M1Lpu**2*m3**2*sqrtS**2*xd1u - 
     &       42*m2**2*m3**2*sqrtS**2*xd1u + 4*m2**2*M3Lpu**2*sqrtS**2*xd1u - 10*m3**2*M3Lpu**2*sqrtS**2*xd1u + 2*M0Lpu**2*m4**2*sqrtS**2*xd1u - 18*m1**2*m4**2*sqrtS**2*xd1u - 
     &       2*M1Lpu**2*m4**2*sqrtS**2*xd1u - 2*m2**2*m4**2*sqrtS**2*xd1u - 2*M3Lpu**2*m4**2*sqrtS**2*xd1u + 2*M0Lpu**2*sqrtS**4*xd1u - 2*m1**2*sqrtS**4*xd1u - 2*M1Lpu**2*sqrtS**4*xd1u - 
     &       10*m2**2*sqrtS**4*xd1u + 20*m3**2*sqrtS**4*xd1u - 2*M3Lpu**2*sqrtS**4*xd1u + 4*m4**2*sqrtS**4*xd1u + 4*sqrtS**6*xd1u - 16*M0Lpu**2*m1**2*m3**2*xd22u + 16*m1**4*m3**2*xd22u + 
     &       16*m1**2*M1Lpu**2*m3**2*xd22u - 16*M0Lpu**2*m2**2*m3**2*xd22u + 32*m1**2*m2**2*m3**2*xd22u + 16*M1Lpu**2*m2**2*m3**2*xd22u + 16*m2**4*m3**2*xd22u + 
     &       16*m1**2*m3**2*M3Lpu**2*xd22u + 16*m2**2*m3**2*M3Lpu**2*xd22u - 16*M0Lpu**2*m1**2*m4**2*xd22u + 16*m1**4*m4**2*xd22u + 16*m1**2*M1Lpu**2*m4**2*xd22u - 
     &       16*M0Lpu**2*m2**2*m4**2*xd22u + 32*m1**2*m2**2*m4**2*xd22u + 16*M1Lpu**2*m2**2*m4**2*xd22u + 16*m2**4*m4**2*xd22u + 16*m1**2*M3Lpu**2*m4**2*xd22u + 
     &       16*m2**2*M3Lpu**2*m4**2*xd22u + 8*M0Lpu**2*m1**2*sqrtS**2*xd22u - 8*m1**4*sqrtS**2*xd22u - 8*m1**2*M1Lpu**2*sqrtS**2*xd22u + 8*M0Lpu**2*m2**2*sqrtS**2*xd22u - 
     &       16*m1**2*m2**2*sqrtS**2*xd22u - 8*M1Lpu**2*m2**2*sqrtS**2*xd22u - 8*m2**4*sqrtS**2*xd22u + 8*M0Lpu**2*m3**2*sqrtS**2*xd22u - 40*m1**2*m3**2*sqrtS**2*xd22u - 
     &       8*M1Lpu**2*m3**2*sqrtS**2*xd22u - 40*m2**2*m3**2*sqrtS**2*xd22u - 8*m1**2*M3Lpu**2*sqrtS**2*xd22u - 8*m2**2*M3Lpu**2*sqrtS**2*xd22u - 8*m3**2*M3Lpu**2*sqrtS**2*xd22u + 
     &       8*M0Lpu**2*m4**2*sqrtS**2*xd22u - 40*m1**2*m4**2*sqrtS**2*xd22u - 8*M1Lpu**2*m4**2*sqrtS**2*xd22u - 40*m2**2*m4**2*sqrtS**2*xd22u - 8*M3Lpu**2*m4**2*sqrtS**2*xd22u - 
     &       4*M0Lpu**2*sqrtS**4*xd22u + 20*m1**2*sqrtS**4*xd22u + 4*M1Lpu**2*sqrtS**4*xd22u + 20*m2**2*sqrtS**4*xd22u + 16*m3**2*sqrtS**4*xd22u + 4*M3Lpu**2*sqrtS**4*xd22u + 
     &       16*m4**2*sqrtS**4*xd22u - 8*sqrtS**6*xd22u - 16*M0Lpu**2*m2**2*m3**2*xd23u + 16*m1**2*m2**2*m3**2*xd23u + 16*M1Lpu**2*m2**2*m3**2*xd23u + 16*m2**4*m3**2*xd23u + 
     &       16*m2**2*m3**2*M3Lpu**2*xd23u - 16*M0Lpu**2*m1**2*m4**2*xd23u + 16*m1**4*m4**2*xd23u + 16*m1**2*M1Lpu**2*m4**2*xd23u - 32*M0Lpu**2*m2**2*m4**2*xd23u + 
     &       48*m1**2*m2**2*m4**2*xd23u + 32*M1Lpu**2*m2**2*m4**2*xd23u + 32*m2**4*m4**2*xd23u + 16*m1**2*M3Lpu**2*m4**2*xd23u + 32*m2**2*M3Lpu**2*m4**2*xd23u + 
     &       4*M0Lpu**2*m1**2*sqrtS**2*xd23u - 4*m1**4*sqrtS**2*xd23u - 4*m1**2*M1Lpu**2*sqrtS**2*xd23u + 12*M0Lpu**2*m2**2*sqrtS**2*xd23u - 16*m1**2*m2**2*sqrtS**2*xd23u - 
     &       12*M1Lpu**2*m2**2*sqrtS**2*xd23u - 12*m2**4*sqrtS**2*xd23u + 4*M0Lpu**2*m3**2*sqrtS**2*xd23u - 4*m1**2*m3**2*sqrtS**2*xd23u - 4*M1Lpu**2*m3**2*sqrtS**2*xd23u - 
     &       36*m2**2*m3**2*sqrtS**2*xd23u - 4*m1**2*M3Lpu**2*sqrtS**2*xd23u - 12*m2**2*M3Lpu**2*sqrtS**2*xd23u - 4*m3**2*M3Lpu**2*sqrtS**2*xd23u + 12*M0Lpu**2*m4**2*sqrtS**2*xd23u - 
     &       44*m1**2*m4**2*sqrtS**2*xd23u - 12*M1Lpu**2*m4**2*sqrtS**2*xd23u - 76*m2**2*m4**2*sqrtS**2*xd23u - 12*M3Lpu**2*m4**2*sqrtS**2*xd23u - 4*M0Lpu**2*sqrtS**4*xd23u + 
     &       12*m1**2*sqrtS**4*xd23u + 4*M1Lpu**2*sqrtS**4*xd23u + 28*m2**2*sqrtS**4*xd23u + 8*m3**2*sqrtS**4*xd23u + 4*M3Lpu**2*sqrtS**4*xd23u + 24*m4**2*sqrtS**4*xd23u - 8*sqrtS**6*xd23u - 
     &       20*M0Lpu**2*m1**2*m2**2*xd2u + 20*m1**4*m2**2*xd2u + 20*m1**2*M1Lpu**2*m2**2*xd2u + 4*M0Lpu**2*m2**4*xd2u + 16*m1**2*m2**4*xd2u - 4*M1Lpu**2*m2**4*xd2u - 4*m2**6*xd2u - 
     &       24*M0Lpu**2*m1**2*m3**2*xd2u + 24*m1**4*m3**2*xd2u + 24*m1**2*M1Lpu**2*m3**2*xd2u - 32*M0Lpu**2*m2**2*m3**2*xd2u + 56*m1**2*m2**2*m3**2*xd2u + 32*M1Lpu**2*m2**2*m3**2*xd2u + 
     &       32*m2**4*m3**2*xd2u + 20*m1**2*m2**2*M3Lpu**2*xd2u - 4*m2**4*M3Lpu**2*xd2u + 24*m1**2*m3**2*M3Lpu**2*xd2u + 32*m2**2*m3**2*M3Lpu**2*xd2u - 36*M0Lpu**2*m1**2*m4**2*xd2u + 
     &       36*m1**4*m4**2*xd2u + 36*m1**2*M1Lpu**2*m4**2*xd2u - 20*M0Lpu**2*m2**2*m4**2*xd2u + 56*m1**2*m2**2*m4**2*xd2u + 20*M1Lpu**2*m2**2*m4**2*xd2u + 20*m2**4*m4**2*xd2u + 
     &       36*m1**2*M3Lpu**2*m4**2*xd2u + 20*m2**2*M3Lpu**2*m4**2*xd2u + 8*M0Lpu**2*m1**2*sqrtS**2*xd2u - 8*m1**4*sqrtS**2*xd2u - 8*m1**2*M1Lpu**2*sqrtS**2*xd2u + 
     &       4*M0Lpu**2*m2**2*sqrtS**2*xd2u - 52*m1**2*m2**2*sqrtS**2*xd2u - 4*M1Lpu**2*m2**2*sqrtS**2*xd2u + 4*m2**4*sqrtS**2*xd2u + 14*M0Lpu**2*m3**2*sqrtS**2*xd2u - 
     &       62*m1**2*m3**2*sqrtS**2*xd2u - 14*M1Lpu**2*m3**2*sqrtS**2*xd2u - 78*m2**2*m3**2*sqrtS**2*xd2u - 8*m1**2*M3Lpu**2*sqrtS**2*xd2u - 4*m2**2*M3Lpu**2*sqrtS**2*xd2u - 
     &       14*m3**2*M3Lpu**2*sqrtS**2*xd2u + 14*M0Lpu**2*m4**2*sqrtS**2*xd2u - 86*m1**2*m4**2*sqrtS**2*xd2u - 14*M1Lpu**2*m4**2*sqrtS**2*xd2u - 54*m2**2*m4**2*sqrtS**2*xd2u - 
     &       14*M3Lpu**2*m4**2*sqrtS**2*xd2u - 2*M0Lpu**2*sqrtS**4*xd2u + 18*m1**2*sqrtS**4*xd2u + 2*M1Lpu**2*sqrtS**4*xd2u + 10*m2**2*sqrtS**4*xd2u + 28*m3**2*sqrtS**4*xd2u + 
     &       2*M3Lpu**2*sqrtS**4*xd2u + 28*m4**2*sqrtS**4*xd2u - 4*sqrtS**6*xd2u + 4*M0Lpu**2*m1**2*m3**2*xd33u - 4*m1**4*m3**2*xd33u - 4*m1**2*M1Lpu**2*m3**2*xd33u + 
     &       12*M0Lpu**2*m2**2*m3**2*xd33u - 16*m1**2*m2**2*m3**2*xd33u - 12*M1Lpu**2*m2**2*m3**2*xd33u - 12*m2**4*m3**2*xd33u - 4*m1**2*m3**2*M3Lpu**2*xd33u - 
     &       12*m2**2*m3**2*M3Lpu**2*xd33u - 4*M0Lpu**2*m1**2*m4**2*xd33u + 4*m1**4*m4**2*xd33u + 4*m1**2*M1Lpu**2*m4**2*xd33u - 12*M0Lpu**2*m2**2*m4**2*xd33u + 16*m1**2*m2**2*m4**2*xd33u + 
     &       12*M1Lpu**2*m2**2*m4**2*xd33u + 12*m2**4*m4**2*xd33u + 4*m1**2*M3Lpu**2*m4**2*xd33u + 12*m2**2*M3Lpu**2*m4**2*xd33u - 4*M0Lpu**2*m3**2*sqrtS**2*xd33u + 
     &       12*m1**2*m3**2*sqrtS**2*xd33u + 4*M1Lpu**2*m3**2*sqrtS**2*xd33u + 28*m2**2*m3**2*sqrtS**2*xd33u + 4*m3**2*M3Lpu**2*sqrtS**2*xd33u + 4*M0Lpu**2*m4**2*sqrtS**2*xd33u - 
     &       12*m1**2*m4**2*sqrtS**2*xd33u - 4*M1Lpu**2*m4**2*sqrtS**2*xd33u - 28*m2**2*m4**2*sqrtS**2*xd33u - 4*M3Lpu**2*m4**2*sqrtS**2*xd33u - 8*m3**2*sqrtS**4*xd33u + 
     &       8*m4**2*sqrtS**4*xd33u + 2*M0Lpu**2*m1**4*xd3u - 2*m1**6*xd3u - 2*m1**4*M1Lpu**2*xd3u - 20*M0Lpu**2*m1**2*m2**2*xd3u + 18*m1**4*m2**2*xd3u + 20*m1**2*M1Lpu**2*m2**2*xd3u + 
     &       2*M0Lpu**2*m2**4*xd3u + 18*m1**2*m2**4*xd3u - 2*M1Lpu**2*m2**4*xd3u - 2*m2**6*xd3u - 8*M0Lpu**2*m2**2*m3**2*xd3u + 8*m1**2*m2**2*m3**2*xd3u + 8*M1Lpu**2*m2**2*m3**2*xd3u + 
     &       8*m2**4*m3**2*xd3u - 2*m1**4*M3Lpu**2*xd3u + 20*m1**2*m2**2*M3Lpu**2*xd3u - 2*m2**4*M3Lpu**2*xd3u + 8*m2**2*m3**2*M3Lpu**2*xd3u - 16*M0Lpu**2*m1**2*m4**2*xd3u + 
     &       16*m1**4*m4**2*xd3u + 16*m1**2*M1Lpu**2*m4**2*xd3u - 24*M0Lpu**2*m2**2*m4**2*xd3u + 40*m1**2*m2**2*m4**2*xd3u + 24*M1Lpu**2*m2**2*m4**2*xd3u + 24*m2**4*m4**2*xd3u + 
     &       16*m1**2*M3Lpu**2*m4**2*xd3u + 24*m2**2*M3Lpu**2*m4**2*xd3u - 4*M0Lpu**2*m1**2*sqrtS**2*xd3u + 8*m1**4*sqrtS**2*xd3u + 4*m1**2*M1Lpu**2*sqrtS**2*xd3u - 
     &       36*m1**2*m2**2*sqrtS**2*xd3u + 4*m2**4*sqrtS**2*xd3u + 2*M0Lpu**2*m3**2*sqrtS**2*xd3u - 2*m1**2*m3**2*sqrtS**2*xd3u - 2*M1Lpu**2*m3**2*sqrtS**2*xd3u - 
     &       18*m2**2*m3**2*sqrtS**2*xd3u + 4*m1**2*M3Lpu**2*sqrtS**2*xd3u - 2*m3**2*M3Lpu**2*sqrtS**2*xd3u + 10*M0Lpu**2*m4**2*sqrtS**2*xd3u - 42*m1**2*m4**2*sqrtS**2*xd3u - 
     &       10*M1Lpu**2*m4**2*sqrtS**2*xd3u - 58*m2**2*m4**2*sqrtS**2*xd3u - 10*M3Lpu**2*m4**2*sqrtS**2*xd3u + 2*M0Lpu**2*sqrtS**4*xd3u - 10*m1**2*sqrtS**4*xd3u - 2*M1Lpu**2*sqrtS**4*xd3u - 
     &       2*m2**2*sqrtS**4*xd3u + 4*m3**2*sqrtS**4*xd3u - 2*M3Lpu**2*sqrtS**4*xd3u + 20*m4**2*sqrtS**4*xd3u + 4*sqrtS**6*xd3u - 
     &       uman**2*(5*xc01u + 5*xc02u + xc03u + 12*xc111u + 16*xc11u + 16*xc121u - 4*xc122u + 12*xc123u + 6*xc12u + 2*xc13u + 16*xc21u + 12*xc223u - 2*xc22u + 8*xc23u - 5*M0Lpu**2*xd0u + 
     &          5*m1**2*xd0u + 5*M1Lpu**2*xd0u + 5*m2**2*xd0u + 5*M3Lpu**2*xd0u - 10*sqrtS**2*xd0u - 12*M0Lpu**2*xd12u + 12*m1**2*xd12u + 12*M1Lpu**2*xd12u + 12*m2**2*xd12u + 
     &          12*M3Lpu**2*xd12u - 24*sqrtS**2*xd12u - 6*M0Lpu**2*xd1u + 6*m1**2*xd1u + 6*M1Lpu**2*xd1u + 6*m2**2*xd1u + 6*M3Lpu**2*xd1u - 12*sqrtS**2*xd1u - 12*M0Lpu**2*xd22u + 
     &          12*m1**2*xd22u + 12*M1Lpu**2*xd22u + 12*m2**2*xd22u + 12*M3Lpu**2*xd22u - 24*sqrtS**2*xd22u + 4*M0Lpu**2*xd23u - 4*m1**2*xd23u - 4*M1Lpu**2*xd23u - 4*m2**2*xd23u - 
     &          4*M3Lpu**2*xd23u + 8*sqrtS**2*xd23u - 16*M0Lpu**2*xd2u + 16*m1**2*xd2u + 16*M1Lpu**2*xd2u + 16*m2**2*xd2u + 16*M3Lpu**2*xd2u - 32*sqrtS**2*xd2u + 2*M0Lpu**2*xd3u - 
     &          2*m1**2*xd3u - 2*M1Lpu**2*xd3u - 2*m2**2*xd3u - 2*M3Lpu**2*xd3u + 4*sqrtS**2*xd3u) + 
     &       tman**2*(xc01u + xc02u + xc03u + 4*xc111u + 4*xc11u + 16*xc121u - 12*xc122u + 4*xc123u + 2*xc12u + 2*xc13u + 4*xc21u + 4*xc223u - 6*xc22u + 8*xc23u - M0Lpu**2*xd0u + 
     &          m1**2*xd0u + M1Lpu**2*xd0u + m2**2*xd0u + M3Lpu**2*xd0u - 2*sqrtS**2*xd0u - 4*M0Lpu**2*xd12u + 4*m1**2*xd12u + 4*M1Lpu**2*xd12u + 4*m2**2*xd12u + 4*M3Lpu**2*xd12u - 
     &          8*sqrtS**2*xd12u - 2*M0Lpu**2*xd1u + 2*m1**2*xd1u + 2*M1Lpu**2*xd1u + 2*m2**2*xd1u + 2*M3Lpu**2*xd1u - 4*sqrtS**2*xd1u - 4*M0Lpu**2*xd22u + 4*m1**2*xd22u + 4*M1Lpu**2*xd22u + 
     &          4*m2**2*xd22u + 4*M3Lpu**2*xd22u - 8*sqrtS**2*xd22u + 12*M0Lpu**2*xd23u - 12*m1**2*xd23u - 12*M1Lpu**2*xd23u - 12*m2**2*xd23u - 12*M3Lpu**2*xd23u + 24*sqrtS**2*xd23u - 
     &          4*M0Lpu**2*xd2u + 4*m1**2*xd2u + 4*M1Lpu**2*xd2u + 4*m2**2*xd2u + 4*M3Lpu**2*xd2u - 8*sqrtS**2*xd2u + 6*M0Lpu**2*xd3u - 6*m1**2*xd3u - 6*M1Lpu**2*xd3u - 6*m2**2*xd3u - 
     &          6*M3Lpu**2*xd3u + 12*sqrtS**2*xd3u) + uman*(48*xc001u - 48*xc002u - 48*xc003u - 2*m1**2*xc01u + 3*m3**2*xc01u + 3*m4**2*xc01u + 3*sqrtS**2*xc01u - 2*m1**2*xc02u + 
     &          3*m3**2*xc02u + 3*m4**2*xc02u + 3*sqrtS**2*xc02u + 2*m1**2*xc03u - m3**2*xc03u + 3*m4**2*xc03u - sqrtS**2*xc03u + 8*m1**2*xc111u + 8*m2**2*xc111u + 8*m3**2*xc111u + 
     &          8*m4**2*xc111u - 8*sqrtS**2*xc111u - 4*m3**2*xc112u + 4*m4**2*xc112u + 4*sqrtS**2*xc112u - 4*m1**2*xc113u + 4*m2**2*xc113u - 4*sqrtS**2*xc113u + 4*m1**2*xc11u + 
     &          8*m2**2*xc11u + 10*m3**2*xc11u + 10*m4**2*xc11u - 2*sqrtS**2*xc11u + 16*m1**2*xc121u + 16*m2**2*xc121u + 8*m3**2*xc121u + 24*m4**2*xc121u - 4*m1**2*xc122u - 12*m2**2*xc122u - 
     &          4*m3**2*xc122u - 12*m4**2*xc122u + 8*sqrtS**2*xc122u + 4*m1**2*xc123u + 12*m2**2*xc123u + 4*m3**2*xc123u + 12*m4**2*xc123u - 8*sqrtS**2*xc123u - 4*m1**2*xc12u + 
     &          4*m2**2*xc12u + 2*m3**2*xc12u + 6*m4**2*xc12u - 2*sqrtS**2*xc12u + 4*m1**2*xc13u - 4*m3**2*xc13u + 8*m4**2*xc13u + 4*m1**2*xc21u + 8*m2**2*xc21u + 10*m3**2*xc21u + 
     &          10*m4**2*xc21u - 2*sqrtS**2*xc21u + 16*sqrtS**2*xc221u + 4*m1**2*xc222u - 4*m2**2*xc222u - 4*sqrtS**2*xc222u + 8*m1**2*xc223u + 8*m2**2*xc223u + 8*m3**2*xc223u + 
     &          8*m4**2*xc223u - 8*sqrtS**2*xc223u - 4*m1**2*xc22u - 12*m2**2*xc22u - 2*m3**2*xc22u - 6*m4**2*xc22u + 6*sqrtS**2*xc22u + 8*m1**2*xc23u + 8*m2**2*xc23u + 4*m3**2*xc23u + 
     &          12*m4**2*xc23u + 8*sqrtS**2*xc23u + 48*M0Lpu**2*xd00u - 48*m1**2*xd00u - 48*M1Lpu**2*xd00u - 48*m2**2*xd00u - 48*M3Lpu**2*xd00u + 96*sqrtS**2*xd00u + 2*M0Lpu**2*m1**2*xd0u - 
     &          2*m1**4*xd0u - 2*m1**2*M1Lpu**2*xd0u - 2*m1**2*m2**2*xd0u - 3*M0Lpu**2*m3**2*xd0u + 3*m1**2*m3**2*xd0u + 3*M1Lpu**2*m3**2*xd0u + 3*m2**2*m3**2*xd0u - 2*m1**2*M3Lpu**2*xd0u + 
     &          3*m3**2*M3Lpu**2*xd0u - 3*M0Lpu**2*m4**2*xd0u + 3*m1**2*m4**2*xd0u + 3*M1Lpu**2*m4**2*xd0u + 3*m2**2*m4**2*xd0u + 3*M3Lpu**2*m4**2*xd0u - 3*M0Lpu**2*sqrtS**2*xd0u + 
     &          7*m1**2*sqrtS**2*xd0u + 3*M1Lpu**2*sqrtS**2*xd0u + 3*m2**2*sqrtS**2*xd0u - 6*m3**2*sqrtS**2*xd0u + 3*M3Lpu**2*sqrtS**2*xd0u - 6*m4**2*sqrtS**2*xd0u - 6*sqrtS**4*xd0u + 
     &          4*M0Lpu**2*m1**2*xd11u - 4*m1**4*xd11u - 4*m1**2*M1Lpu**2*xd11u - 4*M0Lpu**2*m2**2*xd11u + 4*M1Lpu**2*m2**2*xd11u + 4*m2**4*xd11u - 4*m1**2*M3Lpu**2*xd11u + 
     &          4*m2**2*M3Lpu**2*xd11u + 4*M0Lpu**2*sqrtS**2*xd11u + 4*m1**2*sqrtS**2*xd11u - 4*M1Lpu**2*sqrtS**2*xd11u - 12*m2**2*sqrtS**2*xd11u - 4*M3Lpu**2*sqrtS**2*xd11u + 
     &          8*sqrtS**4*xd11u - 4*M0Lpu**2*m1**2*xd12u + 4*m1**4*xd12u + 4*m1**2*M1Lpu**2*xd12u - 12*M0Lpu**2*m2**2*xd12u + 16*m1**2*m2**2*xd12u + 12*M1Lpu**2*m2**2*xd12u + 
     &          12*m2**4*xd12u - 4*M0Lpu**2*m3**2*xd12u + 4*m1**2*m3**2*xd12u + 4*M1Lpu**2*m3**2*xd12u + 4*m2**2*m3**2*xd12u + 4*m1**2*M3Lpu**2*xd12u + 12*m2**2*M3Lpu**2*xd12u + 
     &          4*m3**2*M3Lpu**2*xd12u - 12*M0Lpu**2*m4**2*xd12u + 12*m1**2*m4**2*xd12u + 12*M1Lpu**2*m4**2*xd12u + 12*m2**2*m4**2*xd12u + 12*M3Lpu**2*m4**2*xd12u + 
     &          8*M0Lpu**2*sqrtS**2*xd12u - 16*m1**2*sqrtS**2*xd12u - 8*M1Lpu**2*sqrtS**2*xd12u - 32*m2**2*sqrtS**2*xd12u - 8*m3**2*sqrtS**2*xd12u - 8*M3Lpu**2*sqrtS**2*xd12u - 
     &          24*m4**2*sqrtS**2*xd12u + 16*sqrtS**4*xd12u - 8*M0Lpu**2*sqrtS**2*xd13u + 8*m1**2*sqrtS**2*xd13u + 8*M1Lpu**2*sqrtS**2*xd13u + 8*m2**2*sqrtS**2*xd13u + 
     &          8*M3Lpu**2*sqrtS**2*xd13u - 16*sqrtS**4*xd13u + 4*M0Lpu**2*m1**2*xd1u - 4*m1**4*xd1u - 4*m1**2*M1Lpu**2*xd1u - 4*M0Lpu**2*m2**2*xd1u + 4*M1Lpu**2*m2**2*xd1u + 4*m2**4*xd1u - 
     &          2*M0Lpu**2*m3**2*xd1u + 2*m1**2*m3**2*xd1u + 2*M1Lpu**2*m3**2*xd1u + 2*m2**2*m3**2*xd1u - 4*m1**2*M3Lpu**2*xd1u + 4*m2**2*M3Lpu**2*xd1u + 2*m3**2*M3Lpu**2*xd1u - 
     &          6*M0Lpu**2*m4**2*xd1u + 6*m1**2*m4**2*xd1u + 6*M1Lpu**2*m4**2*xd1u + 6*m2**2*m4**2*xd1u + 6*M3Lpu**2*m4**2*xd1u + 2*M0Lpu**2*sqrtS**2*xd1u + 6*m1**2*sqrtS**2*xd1u - 
     &          2*M1Lpu**2*sqrtS**2*xd1u - 10*m2**2*sqrtS**2*xd1u - 4*m3**2*sqrtS**2*xd1u - 2*M3Lpu**2*sqrtS**2*xd1u - 12*m4**2*sqrtS**2*xd1u + 4*sqrtS**4*xd1u - 8*M0Lpu**2*m1**2*xd22u + 
     &          8*m1**4*xd22u + 8*m1**2*M1Lpu**2*xd22u - 8*M0Lpu**2*m2**2*xd22u + 16*m1**2*m2**2*xd22u + 8*M1Lpu**2*m2**2*xd22u + 8*m2**4*xd22u - 8*M0Lpu**2*m3**2*xd22u + 
     &          8*m1**2*m3**2*xd22u + 8*M1Lpu**2*m3**2*xd22u + 8*m2**2*m3**2*xd22u + 8*m1**2*M3Lpu**2*xd22u + 8*m2**2*M3Lpu**2*xd22u + 8*m3**2*M3Lpu**2*xd22u - 8*M0Lpu**2*m4**2*xd22u + 
     &          8*m1**2*m4**2*xd22u + 8*M1Lpu**2*m4**2*xd22u + 8*m2**2*m4**2*xd22u + 8*M3Lpu**2*m4**2*xd22u + 8*M0Lpu**2*sqrtS**2*xd22u - 24*m1**2*sqrtS**2*xd22u - 
     &          8*M1Lpu**2*sqrtS**2*xd22u - 24*m2**2*sqrtS**2*xd22u - 16*m3**2*sqrtS**2*xd22u - 8*M3Lpu**2*sqrtS**2*xd22u - 16*m4**2*sqrtS**2*xd22u + 16*sqrtS**4*xd22u + 
     &          4*M0Lpu**2*m1**2*xd23u - 4*m1**4*xd23u - 4*m1**2*M1Lpu**2*xd23u + 12*M0Lpu**2*m2**2*xd23u - 16*m1**2*m2**2*xd23u - 12*M1Lpu**2*m2**2*xd23u - 12*m2**4*xd23u + 
     &          4*M0Lpu**2*m3**2*xd23u - 4*m1**2*m3**2*xd23u - 4*M1Lpu**2*m3**2*xd23u - 4*m2**2*m3**2*xd23u - 4*m1**2*M3Lpu**2*xd23u - 12*m2**2*M3Lpu**2*xd23u - 4*m3**2*M3Lpu**2*xd23u + 
     &          12*M0Lpu**2*m4**2*xd23u - 12*m1**2*m4**2*xd23u - 12*M1Lpu**2*m4**2*xd23u - 12*m2**2*m4**2*xd23u - 12*M3Lpu**2*m4**2*xd23u - 8*M0Lpu**2*sqrtS**2*xd23u + 
     &          16*m1**2*sqrtS**2*xd23u + 8*M1Lpu**2*sqrtS**2*xd23u + 32*m2**2*sqrtS**2*xd23u + 8*m3**2*sqrtS**2*xd23u + 8*M3Lpu**2*sqrtS**2*xd23u + 24*m4**2*sqrtS**2*xd23u - 
     &          16*sqrtS**4*xd23u - 4*M0Lpu**2*m1**2*xd2u + 4*m1**4*xd2u + 4*m1**2*M1Lpu**2*xd2u - 8*M0Lpu**2*m2**2*xd2u + 12*m1**2*m2**2*xd2u + 8*M1Lpu**2*m2**2*xd2u + 8*m2**4*xd2u - 
     &          10*M0Lpu**2*m3**2*xd2u + 10*m1**2*m3**2*xd2u + 10*M1Lpu**2*m3**2*xd2u + 10*m2**2*m3**2*xd2u + 4*m1**2*M3Lpu**2*xd2u + 8*m2**2*M3Lpu**2*xd2u + 10*m3**2*M3Lpu**2*xd2u - 
     &          10*M0Lpu**2*m4**2*xd2u + 10*m1**2*m4**2*xd2u + 10*M1Lpu**2*m4**2*xd2u + 10*m2**2*m4**2*xd2u + 10*M3Lpu**2*m4**2*xd2u + 2*M0Lpu**2*sqrtS**2*xd2u - 10*m1**2*sqrtS**2*xd2u - 
     &          2*M1Lpu**2*sqrtS**2*xd2u - 18*m2**2*sqrtS**2*xd2u - 20*m3**2*sqrtS**2*xd2u - 2*M3Lpu**2*sqrtS**2*xd2u - 20*m4**2*sqrtS**2*xd2u + 4*sqrtS**4*xd2u - 4*M0Lpu**2*m1**2*xd33u + 
     &          4*m1**4*xd33u + 4*m1**2*M1Lpu**2*xd33u + 4*M0Lpu**2*m2**2*xd33u - 4*M1Lpu**2*m2**2*xd33u - 4*m2**4*xd33u + 4*m1**2*M3Lpu**2*xd33u - 4*m2**2*M3Lpu**2*xd33u + 
     &          4*M0Lpu**2*sqrtS**2*xd33u - 12*m1**2*sqrtS**2*xd33u - 4*M1Lpu**2*sqrtS**2*xd33u + 4*m2**2*sqrtS**2*xd33u - 4*M3Lpu**2*sqrtS**2*xd33u + 8*sqrtS**4*xd33u + 
     &          4*M0Lpu**2*m1**2*xd3u - 4*m1**4*xd3u - 4*m1**2*M1Lpu**2*xd3u + 12*M0Lpu**2*m2**2*xd3u - 16*m1**2*m2**2*xd3u - 12*M1Lpu**2*m2**2*xd3u - 12*m2**4*xd3u + 2*M0Lpu**2*m3**2*xd3u - 
     &          2*m1**2*m3**2*xd3u - 2*M1Lpu**2*m3**2*xd3u - 2*m2**2*m3**2*xd3u - 4*m1**2*M3Lpu**2*xd3u - 12*m2**2*M3Lpu**2*xd3u - 2*m3**2*M3Lpu**2*xd3u + 6*M0Lpu**2*m4**2*xd3u - 
     &          6*m1**2*m4**2*xd3u - 6*M1Lpu**2*m4**2*xd3u - 6*m2**2*m4**2*xd3u - 6*M3Lpu**2*m4**2*xd3u - 6*M0Lpu**2*sqrtS**2*xd3u + 14*m1**2*sqrtS**2*xd3u + 6*M1Lpu**2*sqrtS**2*xd3u + 
     &          30*m2**2*sqrtS**2*xd3u + 4*m3**2*sqrtS**2*xd3u + 6*M3Lpu**2*sqrtS**2*xd3u + 12*m4**2*sqrtS**2*xd3u - 12*sqrtS**4*xd3u) - 
     &          tman*(48*xc001u - 48*xc002u - 48*xc003u + 8*m1**2*xc01u + 6*m2**2*xc01u + 3*m3**2*xc01u + 3*m4**2*xc01u - sqrtS**2*xc01u - 4*uman*xc01u + 8*m1**2*xc02u + 6*m2**2*xc02u + 
     &          3*m3**2*xc02u + 3*m4**2*xc02u - sqrtS**2*xc02u - 4*uman*xc02u + 2*m2**2*xc03u - m3**2*xc03u + 3*m4**2*xc03u - sqrtS**2*xc03u + 8*m1**2*xc111u + 8*m2**2*xc111u + 
     &          8*m3**2*xc111u + 8*m4**2*xc111u - 8*sqrtS**2*xc111u - 8*uman*xc111u - 4*m3**2*xc112u + 4*m4**2*xc112u + 4*sqrtS**2*xc112u - 4*m1**2*xc113u + 4*m2**2*xc113u - 
     &          4*sqrtS**2*xc113u + 16*m1**2*xc11u + 12*m2**2*xc11u + 10*m3**2*xc11u + 10*m4**2*xc11u - 6*sqrtS**2*xc11u - 12*uman*xc11u + 16*m1**2*xc121u + 16*m2**2*xc121u + 
     &          8*m3**2*xc121u + 24*m4**2*xc121u - 4*m1**2*xc122u - 12*m2**2*xc122u - 4*m3**2*xc122u - 12*m4**2*xc122u + 8*sqrtS**2*xc122u - 8*uman*xc122u + 4*m1**2*xc123u + 
     &          12*m2**2*xc123u + 4*m3**2*xc123u + 12*m4**2*xc123u - 8*sqrtS**2*xc123u - 8*uman*xc123u + 4*m1**2*xc12u + 12*m2**2*xc12u + 2*m3**2*xc12u + 6*m4**2*xc12u - 6*sqrtS**2*xc12u - 
     &          4*uman*xc12u + 4*m2**2*xc13u - 4*m3**2*xc13u + 8*m4**2*xc13u + 16*m1**2*xc21u + 12*m2**2*xc21u + 10*m3**2*xc21u + 10*m4**2*xc21u - 6*sqrtS**2*xc21u - 12*uman*xc21u + 
     &          16*sqrtS**2*xc221u + 4*m1**2*xc222u - 4*m2**2*xc222u - 4*sqrtS**2*xc222u + 8*m1**2*xc223u + 8*m2**2*xc223u + 8*m3**2*xc223u + 8*m4**2*xc223u - 8*sqrtS**2*xc223u - 
     &          8*uman*xc223u + 4*m1**2*xc22u - 4*m2**2*xc22u - 2*m3**2*xc22u - 6*m4**2*xc22u + 2*sqrtS**2*xc22u - 4*uman*xc22u + 8*m1**2*xc23u + 8*m2**2*xc23u + 4*m3**2*xc23u + 
     &          12*m4**2*xc23u + 8*sqrtS**2*xc23u + 48*M0Lpu**2*xd00u - 48*m1**2*xd00u - 48*M1Lpu**2*xd00u - 48*m2**2*xd00u - 48*M3Lpu**2*xd00u + 96*sqrtS**2*xd00u - 8*M0Lpu**2*m1**2*xd0u + 
     &          8*m1**4*xd0u + 8*m1**2*M1Lpu**2*xd0u - 6*M0Lpu**2*m2**2*xd0u + 14*m1**2*m2**2*xd0u + 6*M1Lpu**2*m2**2*xd0u + 6*m2**4*xd0u - 3*M0Lpu**2*m3**2*xd0u + 3*m1**2*m3**2*xd0u + 
     &          3*M1Lpu**2*m3**2*xd0u + 3*m2**2*m3**2*xd0u + 8*m1**2*M3Lpu**2*xd0u + 6*m2**2*M3Lpu**2*xd0u + 3*m3**2*M3Lpu**2*xd0u - 3*M0Lpu**2*m4**2*xd0u + 3*m1**2*m4**2*xd0u + 
     &          3*M1Lpu**2*m4**2*xd0u + 3*m2**2*m4**2*xd0u + 3*M3Lpu**2*m4**2*xd0u + M0Lpu**2*sqrtS**2*xd0u - 17*m1**2*sqrtS**2*xd0u - M1Lpu**2*sqrtS**2*xd0u - 13*m2**2*sqrtS**2*xd0u - 
     &          6*m3**2*sqrtS**2*xd0u - M3Lpu**2*sqrtS**2*xd0u - 6*m4**2*sqrtS**2*xd0u + 2*sqrtS**4*xd0u + 4*M0Lpu**2*uman*xd0u - 4*m1**2*uman*xd0u - 4*M1Lpu**2*uman*xd0u - 
     &          4*m2**2*uman*xd0u - 4*M3Lpu**2*uman*xd0u + 8*sqrtS**2*uman*xd0u + 4*M0Lpu**2*m1**2*xd11u - 4*m1**4*xd11u - 4*m1**2*M1Lpu**2*xd11u - 4*M0Lpu**2*m2**2*xd11u + 
     &          4*M1Lpu**2*m2**2*xd11u + 4*m2**4*xd11u - 4*m1**2*M3Lpu**2*xd11u + 4*m2**2*M3Lpu**2*xd11u + 4*M0Lpu**2*sqrtS**2*xd11u + 4*m1**2*sqrtS**2*xd11u - 4*M1Lpu**2*sqrtS**2*xd11u - 
     &          12*m2**2*sqrtS**2*xd11u - 4*M3Lpu**2*sqrtS**2*xd11u + 8*sqrtS**4*xd11u - 4*M0Lpu**2*m1**2*xd12u + 4*m1**4*xd12u + 4*m1**2*M1Lpu**2*xd12u - 12*M0Lpu**2*m2**2*xd12u + 
     &          16*m1**2*m2**2*xd12u + 12*M1Lpu**2*m2**2*xd12u + 12*m2**4*xd12u - 4*M0Lpu**2*m3**2*xd12u + 4*m1**2*m3**2*xd12u + 4*M1Lpu**2*m3**2*xd12u + 4*m2**2*m3**2*xd12u + 
     &          4*m1**2*M3Lpu**2*xd12u + 12*m2**2*M3Lpu**2*xd12u + 4*m3**2*M3Lpu**2*xd12u - 12*M0Lpu**2*m4**2*xd12u + 12*m1**2*m4**2*xd12u + 12*M1Lpu**2*m4**2*xd12u + 12*m2**2*m4**2*xd12u + 
     &          12*M3Lpu**2*m4**2*xd12u + 8*M0Lpu**2*sqrtS**2*xd12u - 16*m1**2*sqrtS**2*xd12u - 8*M1Lpu**2*sqrtS**2*xd12u - 32*m2**2*sqrtS**2*xd12u - 8*m3**2*sqrtS**2*xd12u - 
     &          8*M3Lpu**2*sqrtS**2*xd12u - 24*m4**2*sqrtS**2*xd12u + 16*sqrtS**4*xd12u + 8*M0Lpu**2*uman*xd12u - 8*m1**2*uman*xd12u - 8*M1Lpu**2*uman*xd12u - 8*m2**2*uman*xd12u - 
     &          8*M3Lpu**2*uman*xd12u + 16*sqrtS**2*uman*xd12u - 8*M0Lpu**2*sqrtS**2*xd13u + 8*m1**2*sqrtS**2*xd13u + 8*M1Lpu**2*sqrtS**2*xd13u + 8*m2**2*sqrtS**2*xd13u + 
     &          8*M3Lpu**2*sqrtS**2*xd13u - 16*sqrtS**4*xd13u - 4*M0Lpu**2*m1**2*xd1u + 4*m1**4*xd1u + 4*m1**2*M1Lpu**2*xd1u - 12*M0Lpu**2*m2**2*xd1u + 16*m1**2*m2**2*xd1u + 
     &          12*M1Lpu**2*m2**2*xd1u + 12*m2**4*xd1u - 2*M0Lpu**2*m3**2*xd1u + 2*m1**2*m3**2*xd1u + 2*M1Lpu**2*m3**2*xd1u + 2*m2**2*m3**2*xd1u + 4*m1**2*M3Lpu**2*xd1u + 
     &          12*m2**2*M3Lpu**2*xd1u + 2*m3**2*M3Lpu**2*xd1u - 6*M0Lpu**2*m4**2*xd1u + 6*m1**2*m4**2*xd1u + 6*M1Lpu**2*m4**2*xd1u + 6*m2**2*m4**2*xd1u + 6*M3Lpu**2*m4**2*xd1u + 
     &          6*M0Lpu**2*sqrtS**2*xd1u - 14*m1**2*sqrtS**2*xd1u - 6*M1Lpu**2*sqrtS**2*xd1u - 30*m2**2*sqrtS**2*xd1u - 4*m3**2*sqrtS**2*xd1u - 6*M3Lpu**2*sqrtS**2*xd1u - 
     &          12*m4**2*sqrtS**2*xd1u + 12*sqrtS**4*xd1u + 4*M0Lpu**2*uman*xd1u - 4*m1**2*uman*xd1u - 4*M1Lpu**2*uman*xd1u - 4*m2**2*uman*xd1u - 4*M3Lpu**2*uman*xd1u + 
     &          8*sqrtS**2*uman*xd1u - 8*M0Lpu**2*m1**2*xd22u + 8*m1**4*xd22u + 8*m1**2*M1Lpu**2*xd22u - 8*M0Lpu**2*m2**2*xd22u + 16*m1**2*m2**2*xd22u + 8*M1Lpu**2*m2**2*xd22u + 
     &          8*m2**4*xd22u - 8*M0Lpu**2*m3**2*xd22u + 8*m1**2*m3**2*xd22u + 8*M1Lpu**2*m3**2*xd22u + 8*m2**2*m3**2*xd22u + 8*m1**2*M3Lpu**2*xd22u + 8*m2**2*M3Lpu**2*xd22u + 
     &          8*m3**2*M3Lpu**2*xd22u - 8*M0Lpu**2*m4**2*xd22u + 8*m1**2*m4**2*xd22u + 8*M1Lpu**2*m4**2*xd22u + 8*m2**2*m4**2*xd22u + 8*M3Lpu**2*m4**2*xd22u + 8*M0Lpu**2*sqrtS**2*xd22u - 
     &          24*m1**2*sqrtS**2*xd22u - 8*M1Lpu**2*sqrtS**2*xd22u - 24*m2**2*sqrtS**2*xd22u - 16*m3**2*sqrtS**2*xd22u - 8*M3Lpu**2*sqrtS**2*xd22u - 16*m4**2*sqrtS**2*xd22u + 
     &          16*sqrtS**4*xd22u + 8*M0Lpu**2*uman*xd22u - 8*m1**2*uman*xd22u - 8*M1Lpu**2*uman*xd22u - 8*m2**2*uman*xd22u - 8*M3Lpu**2*uman*xd22u + 16*sqrtS**2*uman*xd22u + 
     &          4*M0Lpu**2*m1**2*xd23u - 4*m1**4*xd23u - 4*m1**2*M1Lpu**2*xd23u + 12*M0Lpu**2*m2**2*xd23u - 16*m1**2*m2**2*xd23u - 12*M1Lpu**2*m2**2*xd23u - 12*m2**4*xd23u + 
     &          4*M0Lpu**2*m3**2*xd23u - 4*m1**2*m3**2*xd23u - 4*M1Lpu**2*m3**2*xd23u - 4*m2**2*m3**2*xd23u - 4*m1**2*M3Lpu**2*xd23u - 12*m2**2*M3Lpu**2*xd23u - 4*m3**2*M3Lpu**2*xd23u + 
     &          12*M0Lpu**2*m4**2*xd23u - 12*m1**2*m4**2*xd23u - 12*M1Lpu**2*m4**2*xd23u - 12*m2**2*m4**2*xd23u - 12*M3Lpu**2*m4**2*xd23u - 8*M0Lpu**2*sqrtS**2*xd23u + 
     &          16*m1**2*sqrtS**2*xd23u + 8*M1Lpu**2*sqrtS**2*xd23u + 32*m2**2*sqrtS**2*xd23u + 8*m3**2*sqrtS**2*xd23u + 8*M3Lpu**2*sqrtS**2*xd23u + 24*m4**2*sqrtS**2*xd23u - 
     &          16*sqrtS**4*xd23u + 8*M0Lpu**2*uman*xd23u - 8*m1**2*uman*xd23u - 8*M1Lpu**2*uman*xd23u - 8*m2**2*uman*xd23u - 8*M3Lpu**2*uman*xd23u + 16*sqrtS**2*uman*xd23u - 
     &          16*M0Lpu**2*m1**2*xd2u + 16*m1**4*xd2u + 16*m1**2*M1Lpu**2*xd2u - 12*M0Lpu**2*m2**2*xd2u + 28*m1**2*m2**2*xd2u + 12*M1Lpu**2*m2**2*xd2u + 12*m2**4*xd2u - 
     &          10*M0Lpu**2*m3**2*xd2u + 10*m1**2*m3**2*xd2u + 10*M1Lpu**2*m3**2*xd2u + 10*m2**2*m3**2*xd2u + 16*m1**2*M3Lpu**2*xd2u + 12*m2**2*M3Lpu**2*xd2u + 10*m3**2*M3Lpu**2*xd2u - 
     &          10*M0Lpu**2*m4**2*xd2u + 10*m1**2*m4**2*xd2u + 10*M1Lpu**2*m4**2*xd2u + 10*m2**2*m4**2*xd2u + 10*M3Lpu**2*m4**2*xd2u + 6*M0Lpu**2*sqrtS**2*xd2u - 38*m1**2*sqrtS**2*xd2u - 
     &          6*M1Lpu**2*sqrtS**2*xd2u - 30*m2**2*sqrtS**2*xd2u - 20*m3**2*sqrtS**2*xd2u - 6*M3Lpu**2*sqrtS**2*xd2u - 20*m4**2*sqrtS**2*xd2u + 12*sqrtS**4*xd2u + 12*M0Lpu**2*uman*xd2u - 
     &          12*m1**2*uman*xd2u - 12*M1Lpu**2*uman*xd2u - 12*m2**2*uman*xd2u - 12*M3Lpu**2*uman*xd2u + 24*sqrtS**2*uman*xd2u - 4*M0Lpu**2*m1**2*xd33u + 4*m1**4*xd33u + 
     &          4*m1**2*M1Lpu**2*xd33u + 4*M0Lpu**2*m2**2*xd33u - 4*M1Lpu**2*m2**2*xd33u - 4*m2**4*xd33u + 4*m1**2*M3Lpu**2*xd33u - 4*m2**2*M3Lpu**2*xd33u + 4*M0Lpu**2*sqrtS**2*xd33u - 
     &          12*m1**2*sqrtS**2*xd33u - 4*M1Lpu**2*sqrtS**2*xd33u + 4*m2**2*sqrtS**2*xd33u - 4*M3Lpu**2*sqrtS**2*xd33u + 8*sqrtS**4*xd33u - 4*M0Lpu**2*m1**2*xd3u + 4*m1**4*xd3u + 
     &          4*m1**2*M1Lpu**2*xd3u + 4*M0Lpu**2*m2**2*xd3u - 4*M1Lpu**2*m2**2*xd3u - 4*m2**4*xd3u + 2*M0Lpu**2*m3**2*xd3u - 2*m1**2*m3**2*xd3u - 2*M1Lpu**2*m3**2*xd3u - 
     &          2*m2**2*m3**2*xd3u + 4*m1**2*M3Lpu**2*xd3u - 4*m2**2*M3Lpu**2*xd3u - 2*m3**2*M3Lpu**2*xd3u + 6*M0Lpu**2*m4**2*xd3u - 6*m1**2*m4**2*xd3u - 6*M1Lpu**2*m4**2*xd3u - 
     &          6*m2**2*m4**2*xd3u - 6*M3Lpu**2*m4**2*xd3u - 2*M0Lpu**2*sqrtS**2*xd3u - 6*m1**2*sqrtS**2*xd3u + 2*M1Lpu**2*sqrtS**2*xd3u + 10*m2**2*sqrtS**2*xd3u + 4*m3**2*sqrtS**2*xd3u + 
     &          2*M3Lpu**2*sqrtS**2*xd3u + 12*m4**2*sqrtS**2*xd3u - 4*sqrtS**4*xd3u + 4*M0Lpu**2*uman*xd3u - 4*m1**2*uman*xd3u - 4*M1Lpu**2*uman*xd3u - 4*m2**2*uman*xd3u - 
     &          4*M3Lpu**2*uman*xd3u + 8*sqrtS**2*uman*xd3u)))/(64.*pi**2*propac)  )

        enddo           

!================================================================================================================================== 

      enddo

      VV_TX(1) = Box_TQ
      VV_TX(2) = Box_TS
      VV_TX(3) = Box_TV
      VV_TX(4) = Box_TT
      VV_TX(5) = Box_TU

      VV_TX(6) = VV_TX(1) + VV_TX(2) + VV_TX(3) + VV_TX(4) + VV_TX(5)

      VV_UX(1) = Box_UQ
      VV_UX(2) = Box_US
      VV_UX(3) = Box_UV
      VV_UX(4) = Box_UT
      VV_UX(5) = Box_UU

      VV_UX(6) = VV_UX(1) + VV_UX(2) + VV_UX(3) + VV_UX(4) + VV_UX(5)

!==================================================================================================================================

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine QQ2xx_M2GluonTriangleBoxVV(argx,argy,VV_TX,VV_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer j,argx,argy
      complex*16 CA0t,CA1t,CA2t,CA3t,CA4t,CA0u,CA1u,CA2u,CA3u,CA4u,CB0t,CB1t,CB2t,CB3t,CB4t,CB0u,CB1u,CB2u,CB3u,CB4u
      double precision VV_TX(6),VV_UX(6),Box_TQ,Box_TS,Box_TV,Box_TT,Box_TU,Box_UQ,Box_US,Box_UV,Box_UT,Box_UU

      external QQ2xx_BT1VV, QQ2xx_BT2VV

! =================================================================================================================================
!                             ---                     Higgs-Vector-final-states                    ---
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

!=====================================================!Box-Triangle(t,Typ1)_Gluon==================================================

      call QQ2xx_BT1VV(0,tman,m3**2,m1**2,Looptype0,xtype,ytype,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,Looptypet,Loopgent,Loopmeigen1,Looptype1,Loopgen1,CA0t,CA1t,CA2t,CA3t,CA4t) 

      CA0t = CA0t/(16d0*pi**2)
      CA1t = CA1t/(16d0*pi**2)
      CA2t = CA2t/(16d0*pi**2)
      CA3t = CA3t/(16d0*pi**2)
      CA4t = CA4t/(16d0*pi**2)
!=====================================================!Box-Triangle(t,Typ2)_Gluon==================================================

      call QQ2xx_BT2VV(0,m2**2,m4**2,tman,Looptype0,ytype,xtype,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,Loopmeigen2,Looptype2,Loopgen2,Looptypet,Loopgent,CB0t,CB1t,CB2t,CB3t,CB4t)

      CB0t = CB0t/(16d0*pi**2)
      CB1t = CB1t/(16d0*pi**2)
      CB2t = CB2t/(16d0*pi**2)
      CB3t = CB3t/(16d0*pi**2)
      CB4t = CB4t/(16d0*pi**2) 
!=====================================================!Box-Triangle(u,Typ1)_Gluon==================================================

      call QQ2xx_BT1VV(1,uman,m4**2,m1**2,Looptype0,ytype,xtype,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,Looptypeu,Loopgenu,Loopmeigen1,Looptype1,Loopgen1,CA0u,CA1u,CA2u,CA3u,CA4u) 

      CA0u = CA0u/(16d0*pi**2)
      CA1u = CA1u/(16d0*pi**2)
      CA2u = CA2u/(16d0*pi**2)
      CA3u = CA3u/(16d0*pi**2)
      CA4u = CA4u/(16d0*pi**2)
!=====================================================!Box-Triangle(u,Typ2)_Gluon==================================================

      call QQ2xx_BT2VV(1,m2**2,m3**2,uman,Looptype0,xtype,ytype,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,Loopmeigen2,Looptype2,Loopgen2,Looptypeu,Loopgenu,CB0u,CB1u,CB2u,CB3u,CB4u) 

      CB0u = CB0u/(16d0*pi**2)
      CB1u = CB1u/(16d0*pi**2)
      CB2u = CB2u/(16d0*pi**2)
      CB3u = CB3u/(16d0*pi**2)
      CB4u = CB4u/(16d0*pi**2)
!==================================================================================================================================

      do j=1,  2   

!Box(t)*t + Box(u)*t
       
        Coupl_Ac = dconjg(VSfSft(xtype,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))       
        Coupl_Bc = dconjg(VSfSfu(ytype,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
        sfmc = MSf(j,proptypetc,propgentc)
        
!t
	Box_TT = Box_TT + 2*colorfactorTU*dreal(   -(Coupl_Ac*Coupl_Bc*(-2*CA3t*(m2**2 - tman)*(3*m1**2 - m3**2 + tman) + 2*CB2t*(-m1**2 + tman)*(3*m2**2 - m4**2 + tman) + CA3t*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman) + 
     &       CB2t*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman) + 2*CA0t*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*uman) + 
     &       2*CB0t*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*uman) - 2*CB4t*(m1**2 - tman)*(2*m2**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*uman) + 
     &       CB4t*(3*m1**2 - m3**2 + tman)*(2*m2**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*uman) - 2*CA1t*(m2**2 - tman)*(2*m1**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*uman) + 
     &       CA1t*(3*m2**2 - m4**2 + tman)*(2*m1**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*uman) - 2*CA1t*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman) + 
     &       CA1t*(3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman) + CB1t*(3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman) + 
     &       CB3t*(2*m2**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*uman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman) - 2*CB4t*(m1**2 - tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman) + 
     &       CA4t*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman) + CB4t*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman) + 
     &       CA2t*(2*m1**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman) + 
     &       CA2t*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman) + 
     &       CB3t*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)))/(4.*(-sfmc**2 + tman))  )

!u
	Box_UT = Box_UT + 2*colorfactorTU*dreal(  -(Coupl_Ac*Coupl_Bc*(2*CB4u*(m2**2 - tman)*(3*m1**2 - m3**2 + tman) + 2*CA1u*(m1**2 - tman)*(3*m2**2 - m4**2 + tman) + CA2u*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman) + 
     &       CB3u*(3*m1**2 - m3**2 + tman)*(3*m2**2 - m4**2 + tman) + 2*CA0u*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*uman) + 
     &       2*CB0u*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*uman) - 2*CA1u*(m1**2 - tman)*(2*m2**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*uman) - 
     &       CA2u*(3*m1**2 - m3**2 + tman)*(2*m2**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*uman) - 2*CB4u*(m2**2 - tman)*(2*m1**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*uman) - 
     &       CB3u*(3*m2**2 - m4**2 + tman)*(2*m1**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*uman) + 2*CB2u*(m2**2 - tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman) + 
     &       CA1u*(3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman) + CB1u*(3*m2**2 - m4**2 + tman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman) - 
     &       CA1u*(2*m2**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*uman)*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman) + 2*CA3u*(m1**2 - tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman) + 
     &       CA4u*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman) + CB4u*(3*m1**2 - m3**2 + tman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman) - 
     &       CB4u*(2*m1**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman) + 
     &       CA3u*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman) + 
     &       CB2u*(2*m1**2 + 3*m2**2 + m3**2 - 2*sqrtS**2 - uman)*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)))/(4.*(-sfmc**2 + tman))  )

!Box(t)*u + Box(u)*u
 
        Coupl_Ac = dconjg(VSfSfu(ytype,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
        Coupl_Bc = dconjg(VSfSft(xtype,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
        sfmc = MSf(j,proptypeuc,propgenuc)
        

!t
        Box_TU = Box_TU + 2*colorfactorTU*dreal(   -(Coupl_Ac*Coupl_Bc*(2*CA0t*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*tman) + 2*CB0t*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*tman) - 
     &       CB4t*(2*m1**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*tman)*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman) - 
     &       CA1t*(2*m2**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman) + 
     &       CA3t*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman) + 
     &       CB2t*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman) - 2*CA1t*(2*m2**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*tman)*(m1**2 - uman) + 
     &       2*CA3t*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman) - 2*CB4t*(2*m1**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*tman)*(m2**2 - uman) + 
     &       2*CB2t*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(m2**2 - uman) + CB3t*(2*m1**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*tman)*(-3*m2**2 + m3**2 - uman) + 
     &       CA2t*(2*m2**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*tman)*(-3*m1**2 + m4**2 - uman) + CA1t*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman) + 
     &       CB1t*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman) + 2*CA1t*(m1**2 - uman)*(3*m2**2 - m3**2 + uman) + 
     &       CA4t*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman) + CB4t*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman) + 
     &       2*CB4t*(m2**2 - uman)*(3*m1**2 - m4**2 + uman) + CA2t*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman) + CB3t*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)))/
     &       (4.*(-sfmc**2 + uman))  )

!u
	Box_UU = Box_UU + 2*colorfactorTU*dreal(  -(Coupl_Ac*Coupl_Bc*(2*CA0u*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*tman) + 2*CB0u*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*tman) + 
     &       CA2u*(2*m1**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*tman)*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman) + 
     &       CB3u*(2*m2**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman) + 
     &       CA2u*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman) + 
     &       CB3u*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman) - 2*CB4u*(2*m2**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*tman)*(m1**2 - uman) - 
     &       2*CB4u*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(m1**2 - uman) - 2*CA1u*(2*m1**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*tman)*(m2**2 - uman) - 
     &       2*CA1u*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(m2**2 - uman) - CA1u*(2*m1**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*tman)*(-3*m2**2 + m3**2 - uman) - 
     &       CB4u*(2*m2**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*tman)*(-3*m1**2 + m4**2 - uman) + CA1u*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman) + 
     &       CB1u*(2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(3*m2**2 - m3**2 + uman) + 2*CB2u*(-m1**2 + uman)*(3*m2**2 - m3**2 + uman) + 
     &       CA4u*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman) + CB4u*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(3*m1**2 - m4**2 + uman) - 
     &       2*CA3u*(m2**2 - uman)*(3*m1**2 - m4**2 + uman) + CA3u*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman) + CB2u*(3*m2**2 - m3**2 + uman)*(3*m1**2 - m4**2 + uman)))/
     &       (4.*(-sfmc**2 + uman))  )

       enddo

!================================================================================================================================== 
!Box(t)*q + Box(u)*q

        Coupl_Ac = dconjg(VVSfSf(xtype,ytype,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))

!t
	Box_TQ = Box_TQ + 2*colorfactorQ*dreal(  (Coupl_Ac*(8*CA0t + 8*CB0t + 2*CA2t*m1**2 + 2*CA3t*m1**2 + 2*CA4t*m1**2 + CB2t*m1**2 + CB3t*m1**2 + CA2t*m2**2 + CA3t*m2**2 + 2*CB1t*m2**2 + 2*CB2t*m2**2 + 2*CB3t*m2**2 + CB2t*m3**2 + 
     &      CB3t*m3**2 + CA2t*m4**2 + CA3t*m4**2 - CA2t*sqrtS**2 - CA3t*sqrtS**2 - CB2t*sqrtS**2 - CB3t*sqrtS**2 + 2*CA1t*tman + 2*CB4t*tman - CA2t*uman - CA3t*uman - CB2t*uman - CB3t*uman))/
     &      2.  )

!u
	Box_UQ = Box_UQ + 2*colorfactorQ*dreal(  (Coupl_Ac*(8*CA0u + 8*CB0u + 2*CA2u*m1**2 + 2*CA3u*m1**2 + 2*CA4u*m1**2 + CB2u*m1**2 + CB3u*m1**2 + CA2u*m2**2 + CA3u*m2**2 + 2*CB1u*m2**2 + 2*CB2u*m2**2 + 2*CB3u*m2**2 + CA2u*m3**2 + 
     &      CA3u*m3**2 + CB2u*m4**2 + CB3u*m4**2 - CA2u*sqrtS**2 - CA3u*sqrtS**2 - CB2u*sqrtS**2 - CB3u*sqrtS**2 - CA2u*tman - CA3u*tman - CB2u*tman - CB3u*tman + 2*CA1u*uman + 2*CB4u*uman))/
     &      2.  )

!==================================================================================================================================     
!Box(t)*s(Higgs) + Box(u)*s(Higgs)
       
      do j=sumstartH, sumendH

        Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(LLR(xtype,ytype,j))
        propac = dconjg(DenSmH(j))
        
	Box_TS = Box_TS + 2*colorfactorSH*dreal(  -(Coupl_Ac*Coupl_Bc*(8*CA0t + 8*CB0t + 2*CA2t*m1**2 + 2*CA3t*m1**2 + 2*CA4t*m1**2 + CB2t*m1**2 + CB3t*m1**2 + CA2t*m2**2 + CA3t*m2**2 + 2*CB1t*m2**2 + 2*CB2t*m2**2 + 2*CB3t*m2**2 + 
     &       CB2t*m3**2 + CB3t*m3**2 + CA2t*m4**2 + CA3t*m4**2 - CA2t*sqrtS**2 - CA3t*sqrtS**2 - CB2t*sqrtS**2 - CB3t*sqrtS**2 + 2*CA1t*tman + 2*CB4t*tman - CA2t*uman - CA3t*uman - 
     &       CB2t*uman - CB3t*uman))/(2.*propac)  )

	Box_US = Box_US + 2*colorfactorSH*dreal(   -(Coupl_Ac*Coupl_Bc*(8*CA0u + 8*CB0u + 2*CA2u*m1**2 + 2*CA3u*m1**2 + 2*CA4u*m1**2 + CB2u*m1**2 + CB3u*m1**2 + CA2u*m2**2 + CA3u*m2**2 + 2*CB1u*m2**2 + 2*CB2u*m2**2 + 2*CB3u*m2**2 + 
     &       CA2u*m3**2 + CA3u*m3**2 + CB2u*m4**2 + CB3u*m4**2 - CA2u*sqrtS**2 - CA3u*sqrtS**2 - CB2u*sqrtS**2 - CB3u*sqrtS**2 - CA2u*tman - CA3u*tman - CB2u*tman - CB3u*tman + 2*CA1u*uman + 
     &       2*CB4u*uman))/(2.*propac)  )
  

      enddo

!================================================================================================================================== 
!Box(t)*Vector + Box(u)*Vector

      do j=sumstartV, sumendV

        Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = dconjg(VVV(j,xtype,ytype))
        propac = dconjg(DenSmV(j))
        
        Box_TV = Box_TV + 2*colorfactorSV*dreal(  -(Coupl_Ac*Coupl_Bc*(-(CB2t*m1**4) + CB3t*m1**4 + 3*CB2t*m1**2*m2**2 - 3*CB3t*m1**2*m2**2 - 2*CB2t*m2**4 + 2*CB3t*m2**4 + 4*CA1t*m1**2*m3**2 - 3*CA4t*m1**2*m3**2 + CB1t*m1**2*m3**2 + 
     &       CB2t*m1**2*m3**2 + 3*CB3t*m1**2*m3**2 + 4*CB4t*m1**2*m3**2 + 4*CA1t*m2**2*m3**2 - CA4t*m2**2*m3**2 + 3*CB1t*m2**2*m3**2 - CB2t*m2**2*m3**2 + 9*CB3t*m2**2*m3**2 + 
     &       4*CB4t*m2**2*m3**2 + 4*CA1t*m1**2*m4**2 + 3*CA4t*m1**2*m4**2 - CB1t*m1**2*m4**2 + 2*CB2t*m1**2*m4**2 - 2*CB3t*m1**2*m4**2 + 4*CB4t*m1**2*m4**2 + 4*CA1t*m2**2*m4**2 + 
     &       CA4t*m2**2*m4**2 - 3*CB1t*m2**2*m4**2 + 6*CB2t*m2**2*m4**2 - 2*CB3t*m2**2*m4**2 + 4*CB4t*m2**2*m4**2 - 2*CA1t*m1**2*sqrtS**2 + 2*CB2t*m1**2*sqrtS**2 - 3*CB3t*m1**2*sqrtS**2 - 
     &       2*CB4t*m1**2*sqrtS**2 - 2*CA1t*m2**2*sqrtS**2 + 2*CB2t*m2**2*sqrtS**2 - 5*CB3t*m2**2*sqrtS**2 - 2*CB4t*m2**2*sqrtS**2 - 2*CA1t*m3**2*sqrtS**2 + CA4t*m3**2*sqrtS**2 - 
     &       CB1t*m3**2*sqrtS**2 - 3*CB3t*m3**2*sqrtS**2 - 2*CB4t*m3**2*sqrtS**2 - 2*CA1t*m4**2*sqrtS**2 - CA4t*m4**2*sqrtS**2 + CB1t*m4**2*sqrtS**2 - 2*CB2t*m4**2*sqrtS**2 + 
     &       CB3t*m4**2*sqrtS**2 - 2*CB4t*m4**2*sqrtS**2 + CA1t*sqrtS**4 - CB2t*sqrtS**4 + 2*CB3t*sqrtS**4 + CB4t*sqrtS**4 - 12*CA0t*tman - 12*CB0t*tman + 2*CA1t*m1**2*tman - 
     &       CA4t*m1**2*tman + CB1t*m1**2*tman - 2*CB2t*m1**2*tman + CB3t*m1**2*tman + 2*CB4t*m1**2*tman + 2*CA1t*m2**2*tman + CA4t*m2**2*tman - CB1t*m2**2*tman + CB2t*m2**2*tman - 
     &       4*CB3t*m2**2*tman + 2*CB4t*m2**2*tman + 2*CA1t*m3**2*tman - CB2t*m3**2*tman - 2*CB3t*m3**2*tman + 2*CB4t*m3**2*tman + 2*CA1t*m4**2*tman - CB3t*m4**2*tman + 2*CB4t*m4**2*tman - 
     &       2*CA1t*sqrtS**2*tman - CA4t*sqrtS**2*tman - CB1t*sqrtS**2*tman + 2*CB2t*sqrtS**2*tman - 2*CB4t*sqrtS**2*tman - 3*CA1t*tman**2 + CB3t*tman**2 - 3*CB4t*tman**2 + 12*CA0t*uman + 
     &       12*CB0t*uman - 2*CA1t*m1**2*uman + CA4t*m1**2*uman - CB1t*m1**2*uman + CB2t*m1**2*uman - 2*CB4t*m1**2*uman - 2*CA1t*m2**2*uman - CA4t*m2**2*uman + CB1t*m2**2*uman + 
     &       3*CB3t*m2**2*uman - 2*CB4t*m2**2*uman - 2*CA1t*m3**2*uman + CB2t*m3**2*uman + 2*CB3t*m3**2*uman - 2*CB4t*m3**2*uman - 2*CA1t*m4**2*uman + CB3t*m4**2*uman - 2*CB4t*m4**2*uman + 
     &       2*CA1t*sqrtS**2*uman + CA4t*sqrtS**2*uman + CB1t*sqrtS**2*uman - 2*CB2t*sqrtS**2*uman + 2*CB4t*sqrtS**2*uman + 2*CA1t*tman*uman + CB2t*tman*uman + CB3t*tman*uman + 
     &       2*CB4t*tman*uman + CA1t*uman**2 - CB2t*uman**2 - 2*CB3t*uman**2 + CB4t*uman**2 + 
     &       CA2t*(2*m1**4 + m2**4 + m3**2*sqrtS**2 - 3*m4**2*sqrtS**2 + 2*sqrtS**4 - m3**2*tman - 2*m4**2*tman + tman**2 + m2**2*(-2*m3**2 + 3*m4**2 - 3*sqrtS**2 + tman) - 
     &       m1**2*(3*m2**2 + 2*m3**2 - 9*m4**2 + 5*sqrtS**2 + 4*tman - 3*uman) + m3**2*uman + 2*m4**2*uman + tman*uman - 2*uman**2) + 
     &       CA3t*(-2*m1**4 - m2**4 - 2*m3**2*sqrtS**2 - sqrtS**4 - m4**2*tman + 2*sqrtS**2*tman + m1**2*(3*m2**2 + 6*m3**2 - m4**2 + 2*sqrtS**2 + tman) + m4**2*uman - 2*sqrtS**2*uman + 
     &       tman*uman - uman**2 + m2**2*(2*m3**2 + m4**2 + 2*sqrtS**2 - 2*tman + uman))))/(4.*propac)  )

        Box_UV = Box_UV + 2*colorfactorSV*dreal(  (Coupl_Ac*Coupl_Bc*(-(CB2u*m1**4) + CB3u*m1**4 + 3*CB2u*m1**2*m2**2 - 3*CB3u*m1**2*m2**2 - 2*CB2u*m2**4 + 2*CB3u*m2**4 + 4*CA1u*m1**2*m3**2 + 3*CA4u*m1**2*m3**2 - CB1u*m1**2*m3**2 + 
     &      2*CB2u*m1**2*m3**2 - 2*CB3u*m1**2*m3**2 + 4*CB4u*m1**2*m3**2 + 4*CA1u*m2**2*m3**2 + CA4u*m2**2*m3**2 - 3*CB1u*m2**2*m3**2 + 6*CB2u*m2**2*m3**2 - 2*CB3u*m2**2*m3**2 + 
     &      4*CB4u*m2**2*m3**2 + 4*CA1u*m1**2*m4**2 - 3*CA4u*m1**2*m4**2 + CB1u*m1**2*m4**2 + CB2u*m1**2*m4**2 + 3*CB3u*m1**2*m4**2 + 4*CB4u*m1**2*m4**2 + 4*CA1u*m2**2*m4**2 - 
     &      CA4u*m2**2*m4**2 + 3*CB1u*m2**2*m4**2 - CB2u*m2**2*m4**2 + 9*CB3u*m2**2*m4**2 + 4*CB4u*m2**2*m4**2 - 2*CA1u*m1**2*sqrtS**2 + 2*CB2u*m1**2*sqrtS**2 - 3*CB3u*m1**2*sqrtS**2 - 
     &      2*CB4u*m1**2*sqrtS**2 - 2*CA1u*m2**2*sqrtS**2 + 2*CB2u*m2**2*sqrtS**2 - 5*CB3u*m2**2*sqrtS**2 - 2*CB4u*m2**2*sqrtS**2 - 2*CA1u*m3**2*sqrtS**2 - CA4u*m3**2*sqrtS**2 + 
     &      CB1u*m3**2*sqrtS**2 - 2*CB2u*m3**2*sqrtS**2 + CB3u*m3**2*sqrtS**2 - 2*CB4u*m3**2*sqrtS**2 - 2*CA1u*m4**2*sqrtS**2 + CA4u*m4**2*sqrtS**2 - CB1u*m4**2*sqrtS**2 - 
     &      3*CB3u*m4**2*sqrtS**2 - 2*CB4u*m4**2*sqrtS**2 + CA1u*sqrtS**4 - CB2u*sqrtS**4 + 2*CB3u*sqrtS**4 + CB4u*sqrtS**4 + 12*CA0u*tman + 12*CB0u*tman - 2*CA1u*m1**2*tman + 
     &      CA4u*m1**2*tman - CB1u*m1**2*tman + CB2u*m1**2*tman - 2*CB4u*m1**2*tman - 2*CA1u*m2**2*tman - CA4u*m2**2*tman + CB1u*m2**2*tman + 3*CB3u*m2**2*tman - 2*CB4u*m2**2*tman - 
     &      2*CA1u*m3**2*tman + CB3u*m3**2*tman - 2*CB4u*m3**2*tman - 2*CA1u*m4**2*tman + CB2u*m4**2*tman + 2*CB3u*m4**2*tman - 2*CB4u*m4**2*tman + 2*CA1u*sqrtS**2*tman + 
     &      CA4u*sqrtS**2*tman + CB1u*sqrtS**2*tman - 2*CB2u*sqrtS**2*tman + 2*CB4u*sqrtS**2*tman + CA1u*tman**2 - CB2u*tman**2 - 2*CB3u*tman**2 + CB4u*tman**2 - 12*CA0u*uman - 
     &      12*CB0u*uman + 2*CA1u*m1**2*uman - CA4u*m1**2*uman + CB1u*m1**2*uman - 2*CB2u*m1**2*uman + CB3u*m1**2*uman + 2*CB4u*m1**2*uman + 2*CA1u*m2**2*uman + CA4u*m2**2*uman - 
     &      CB1u*m2**2*uman + CB2u*m2**2*uman - 4*CB3u*m2**2*uman + 2*CB4u*m2**2*uman + 2*CA1u*m3**2*uman - CB3u*m3**2*uman + 2*CB4u*m3**2*uman + 2*CA1u*m4**2*uman - CB2u*m4**2*uman - 
     &      2*CB3u*m4**2*uman + 2*CB4u*m4**2*uman - 2*CA1u*sqrtS**2*uman - CA4u*sqrtS**2*uman - CB1u*sqrtS**2*uman + 2*CB2u*sqrtS**2*uman - 2*CB4u*sqrtS**2*uman + 2*CA1u*tman*uman + 
     &      CB2u*tman*uman + CB3u*tman*uman + 2*CB4u*tman*uman - 3*CA1u*uman**2 + CB3u*uman**2 - 3*CB4u*uman**2 + 
     &      CA3u*(-2*m1**4 - m2**4 - 2*m4**2*sqrtS**2 - sqrtS**4 + m3**2*tman - 2*sqrtS**2*tman - tman**2 + m2**2*(m3**2 + 2*m4**2 + 2*sqrtS**2 + tman - 2*uman) - m3**2*uman + 
     &      2*sqrtS**2*uman + tman*uman + m1**2*(3*m2**2 - m3**2 + 6*m4**2 + 2*sqrtS**2 + uman)) + 
     &      CA2u*(2*m1**4 + m2**4 - 3*m3**2*sqrtS**2 + m4**2*sqrtS**2 + 2*sqrtS**4 + 2*m3**2*tman + m4**2*tman - 2*tman**2 - 2*m3**2*uman - m4**2*uman + tman*uman + uman**2 + 
     &      m2**2*(3*m3**2 - 2*m4**2 - 3*sqrtS**2 + uman) - m1**2*(3*m2**2 - 9*m3**2 + 2*m4**2 + 5*sqrtS**2 - 3*tman + 4*uman))))/(4.*propac)  )

      enddo        
   
!================================================================================================================================== 

      VV_TX(1) = Box_TQ
      VV_TX(2) = Box_TS
      VV_TX(3) = Box_TV
      VV_TX(4) = Box_TT
      VV_TX(5) = Box_TU

      VV_TX(6) = VV_TX(1) + VV_TX(2) + VV_TX(3) + VV_TX(4) + VV_TX(5)

      VV_UX(1) = Box_UQ
      VV_UX(2) = Box_US
      VV_UX(3) = Box_UV
      VV_UX(4) = Box_UT
      VV_UX(5) = Box_UU

      VV_UX(6) = VV_UX(1) + VV_UX(2) + VV_UX(3) + VV_UX(4) + VV_UX(5)

!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine QQ2xx_M2sQuarkTriangleBoxVV(argx,argy,VV_TX,VV_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer j,argx,argy
      complex*16 CA0t,CA1t,CA2t,CA3t,CA4t,CB0u,CB1u,CB2u,CB3u,CB4u
      double precision VV_TX(6),VV_UX(6),Box_TQ,Box_TS,Box_TV,Box_TT,Box_TU,Box_UQ,Box_US,Box_UV,Box_UT,Box_UU

      external QQ2xx_BTsf4VV

! =================================================================================================================================
!                             ---                     Higgs-Vector-final-states                    ---
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

!====================================================!Box-Triangel(t) (4-sQuark-Vertex)============================================

      call QQ2xx_BTsf4VV(0,sqrtS**2,m4**2,m3**2,ytype,xtype,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,CA0t,CA1t,CA2t,CA3t,CA4t) 

      CA0t = CA0t/(16d0*pi**2)
      CA1t = CA1t/(16d0*pi**2)
      CA2t = CA2t/(16d0*pi**2)
      CA3t = CA3t/(16d0*pi**2)
      CA4t = CA4t/(16d0*pi**2)
!====================================================!Box-Triangel(u) (4-sQuark-Vertex)============================================

      call QQ2xx_BTsf4VV(1,sqrtS**2,m3**2,m4**2,xtype,ytype,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,CB0u,CB1u,CB2u,CB3u,CB4u) 

      CB0u = CB0u/(16d0*pi**2)
      CB1u = CB1u/(16d0*pi**2)
      CB2u = CB2u/(16d0*pi**2)
      CB3u = CB3u/(16d0*pi**2)
      CB4u = CB4u/(16d0*pi**2)
!==================================================================================================================================

      do j=1, 2   

!Box(t)*t + Box(u)*t
       
        Coupl_Ac = dconjg(VSfSft(xtype,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))       
        Coupl_Bc = dconjg(VSfSfu(ytype,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
        sfmc = MSf(j,proptypetc,propgentc)

!t        
        Box_TT = Box_TT + 2*dreal(   (Coupl_Ac*Coupl_Bc*(6*CA2t*m1**4 - 2*CA2t*m1**2*m2**2 - 2*CA3t*m1**2*m2**2 - 4*CA4t*m1**2*m2**2 + 6*CA3t*m2**4 - 3*CA3t*m1**2*m3**2 - 6*CA4t*m1**2*m3**2 + 13*CA3t*m2**2*m3**2 + 
     &      6*CA3t*m3**4 + 4*CA2t*m1**2*m4**2 - CA3t*m1**2*m4**2 - 2*CA4t*m1**2*m4**2 + 3*CA3t*m2**2*m4**2 + 2*CA3t*m3**2*m4**2 - 4*CA2t*m1**2*sqrtS**2 + CA3t*m1**2*sqrtS**2 + 
     &      2*CA4t*m1**2*sqrtS**2 - 7*CA3t*m2**2*sqrtS**2 - 8*CA3t*m3**2*sqrtS**2 - 2*CA3t*m4**2*sqrtS**2 + 2*CA3t*sqrtS**4 - 8*CA2t*m1**2*tman + 2*CA2t*m2**2*tman - 2*CA3t*m2**2*tman + 
     &      4*CA4t*m2**2*tman - 3*CA3t*m3**2*tman + 6*CA4t*m3**2*tman - 4*CA2t*m4**2*tman - CA3t*m4**2*tman + 2*CA4t*m4**2*tman + 4*CA2t*sqrtS**2*tman + CA3t*sqrtS**2*tman - 
     &      2*CA4t*sqrtS**2*tman + 2*CA2t*tman**2 + 2*CA0t*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*uman) - 2*CA2t*m1**2*uman + 2*CA3t*m1**2*uman + 4*CA4t*m1**2*uman - 
     &      8*CA3t*m2**2*uman - 7*CA3t*m3**2*uman - CA3t*m4**2*uman + 5*CA3t*sqrtS**2*uman + 2*CA2t*tman*uman + 2*CA3t*tman*uman - 4*CA4t*tman*uman + 2*CA3t*uman**2 + 
     &      CA1t*(3*m1**2 - m2**2 + 2*m4**2 - 2*sqrtS**2 - tman - uman)*(m1**2 - 3*m2**2 - 2*m3**2 + 2*sqrtS**2 + tman + uman)))/(4.*(sfmc**2 - tman))   )
!u
        Box_UT = Box_UT + 2*dreal(   (Coupl_Ac*Coupl_Bc*(6*CB3u*m1**4 - 2*CB2u*m1**2*m2**2 - 2*CB3u*m1**2*m2**2 - 4*CB4u*m1**2*m2**2 + 6*CB2u*m2**4 + 3*CB3u*m1**2*m3**2 + 4*CB2u*m2**2*m3**2 - CB3u*m2**2*m3**2 - 
     &      2*CB4u*m2**2*m3**2 + 13*CB3u*m1**2*m4**2 - 3*CB3u*m2**2*m4**2 - 6*CB4u*m2**2*m4**2 + 2*CB3u*m3**2*m4**2 + 6*CB3u*m4**4 - 7*CB3u*m1**2*sqrtS**2 - 4*CB2u*m2**2*sqrtS**2 + 
     &      CB3u*m2**2*sqrtS**2 + 2*CB4u*m2**2*sqrtS**2 - 2*CB3u*m3**2*sqrtS**2 - 8*CB3u*m4**2*sqrtS**2 + 2*CB3u*sqrtS**4 + 2*CB2u*m1**2*tman - 2*CB3u*m1**2*tman + 4*CB4u*m1**2*tman - 
     &      8*CB2u*m2**2*tman - 4*CB2u*m3**2*tman - CB3u*m3**2*tman + 2*CB4u*m3**2*tman - 3*CB3u*m4**2*tman + 6*CB4u*m4**2*tman + 4*CB2u*sqrtS**2*tman + CB3u*sqrtS**2*tman - 
     &      2*CB4u*sqrtS**2*tman + 2*CB2u*tman**2 + 2*CB0u*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*uman) - 8*CB3u*m1**2*uman - 2*CB2u*m2**2*uman + 2*CB3u*m2**2*uman + 
     &      4*CB4u*m2**2*uman - CB3u*m3**2*uman - 7*CB3u*m4**2*uman + 5*CB3u*sqrtS**2*uman + 2*CB2u*tman*uman + 2*CB3u*tman*uman - 4*CB4u*tman*uman + 2*CB3u*uman**2 + 
     &      CB1u*(3*m1**2 - m2**2 + 2*m4**2 - 2*sqrtS**2 - tman - uman)*(m1**2 - 3*m2**2 - 2*m3**2 + 2*sqrtS**2 + tman + uman)))/(4.*(sfmc**2 - tman))  )

!Box(t)*u + Box(u)*u
 
        Coupl_Ac = dconjg(VSfSfu(ytype,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
        Coupl_Bc = dconjg(VSfSft(xtype,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
        sfmc = MSf(j,proptypeuc,propgenuc)
        
!t        
	Box_TU = Box_TU + 2*dreal(   (Coupl_Ac*Coupl_Bc*(6*CA3t*m1**4 - 2*CA2t*m1**2*m2**2 - 2*CA3t*m1**2*m2**2 - 4*CA4t*m1**2*m2**2 + 6*CA2t*m2**4 + 13*CA3t*m1**2*m3**2 - 3*CA3t*m2**2*m3**2 - 6*CA4t*m2**2*m3**2 + 
     &      6*CA3t*m3**4 + 3*CA3t*m1**2*m4**2 + 4*CA2t*m2**2*m4**2 - CA3t*m2**2*m4**2 - 2*CA4t*m2**2*m4**2 + 2*CA3t*m3**2*m4**2 - 7*CA3t*m1**2*sqrtS**2 - 4*CA2t*m2**2*sqrtS**2 + 
     &      CA3t*m2**2*sqrtS**2 + 2*CA4t*m2**2*sqrtS**2 - 8*CA3t*m3**2*sqrtS**2 - 2*CA3t*m4**2*sqrtS**2 + 2*CA3t*sqrtS**4 + 
     &      2*CA0t*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*tman) - 8*CA3t*m1**2*tman - 2*CA2t*m2**2*tman + 2*CA3t*m2**2*tman + 4*CA4t*m2**2*tman - 7*CA3t*m3**2*tman - 
     &      CA3t*m4**2*tman + 5*CA3t*sqrtS**2*tman + 2*CA3t*tman**2 + 2*CA2t*m1**2*uman - 2*CA3t*m1**2*uman + 4*CA4t*m1**2*uman - 8*CA2t*m2**2*uman - 3*CA3t*m3**2*uman + 6*CA4t*m3**2*uman - 
     &      4*CA2t*m4**2*uman - CA3t*m4**2*uman + 2*CA4t*m4**2*uman + 4*CA2t*sqrtS**2*uman + CA3t*sqrtS**2*uman - 2*CA4t*sqrtS**2*uman + 2*CA2t*tman*uman + 2*CA3t*tman*uman - 
     &      4*CA4t*tman*uman + 2*CA2t*uman**2 + CA1t*(3*m1**2 - m2**2 + 2*m3**2 - 2*sqrtS**2 - tman - uman)*(m1**2 - 3*m2**2 - 2*m4**2 + 2*sqrtS**2 + tman + uman)))/(4.*(sfmc**2 - uman))  )
!u
        Box_UU = Box_UU + 2*dreal(   (Coupl_Ac*Coupl_Bc*(6*CB2u*m1**4 - 2*CB2u*m1**2*m2**2 - 2*CB3u*m1**2*m2**2 - 4*CB4u*m1**2*m2**2 + 6*CB3u*m2**4 + 4*CB2u*m1**2*m3**2 - CB3u*m1**2*m3**2 - 2*CB4u*m1**2*m3**2 + 
     &      3*CB3u*m2**2*m3**2 - 3*CB3u*m1**2*m4**2 - 6*CB4u*m1**2*m4**2 + 13*CB3u*m2**2*m4**2 + 2*CB3u*m3**2*m4**2 + 6*CB3u*m4**4 - 4*CB2u*m1**2*sqrtS**2 + CB3u*m1**2*sqrtS**2 + 
     &      2*CB4u*m1**2*sqrtS**2 - 7*CB3u*m2**2*sqrtS**2 - 2*CB3u*m3**2*sqrtS**2 - 8*CB3u*m4**2*sqrtS**2 + 2*CB3u*sqrtS**4 + 
     &      2*CB0u*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*tman) - 2*CB2u*m1**2*tman + 2*CB3u*m1**2*tman + 4*CB4u*m1**2*tman - 8*CB3u*m2**2*tman - CB3u*m3**2*tman - 
     &      7*CB3u*m4**2*tman + 5*CB3u*sqrtS**2*tman + 2*CB3u*tman**2 - 8*CB2u*m1**2*uman + 2*CB2u*m2**2*uman - 2*CB3u*m2**2*uman + 4*CB4u*m2**2*uman - 4*CB2u*m3**2*uman - CB3u*m3**2*uman + 
     &      2*CB4u*m3**2*uman - 3*CB3u*m4**2*uman + 6*CB4u*m4**2*uman + 4*CB2u*sqrtS**2*uman + CB3u*sqrtS**2*uman - 2*CB4u*sqrtS**2*uman + 2*CB2u*tman*uman + 2*CB3u*tman*uman - 
     &      4*CB4u*tman*uman + 2*CB2u*uman**2 + CB1u*(3*m1**2 - m2**2 + 2*m3**2 - 2*sqrtS**2 - tman - uman)*(m1**2 - 3*m2**2 - 2*m4**2 + 2*sqrtS**2 + tman + uman)))/(4.*(sfmc**2 - uman))  )


      enddo

!================================================================================================================================== 
!Box(t)*q + Box(u)*q

        Coupl_Ac = dconjg(VVSfSf(xtype,ytype,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))

!t
        Box_TQ = Box_TQ + 2*dreal(   (Coupl_Ac*(8*CA0t + CA3t*m1**2 + CA3t*m2**2 + 2*CA3t*m3**2 + 2*CA4t*m3**2 + 2*CA1t*sqrtS**2 - CA3t*tman + CA2t*(m1**2 + m2**2 + 2*m3**2 - tman - uman) - CA3t*uman))/2.  )
!u
        Box_UQ = Box_UQ + 2*dreal(   (Coupl_Ac*(8*CB0u + CB3u*m1**2 + CB3u*m2**2 + 2*CB3u*m4**2 + 2*CB4u*m4**2 + 2*CB1u*sqrtS**2 - CB3u*tman + CB2u*(m1**2 + m2**2 + 2*m4**2 - tman - uman) - CB3u*uman))/2.  )

!==================================================================================================================================  
!Box(t)*s(Higgs) + Box(u)*s(Higgs)  
       
      do j=sumstartH, sumendH

        Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(LLR(xtype,ytype,j))
        propac = dconjg(DenSmH(j))

!t        
        Box_TS = Box_TS + 2*dreal(  -(Coupl_Ac*Coupl_Bc*(8*CA0t + CA3t*m1**2 + CA3t*m2**2 + 2*CA3t*m3**2 + 2*CA4t*m3**2 + 2*CA1t*sqrtS**2 - CA3t*tman + CA2t*(m1**2 + m2**2 + 2*m3**2 - tman - uman) - CA3t*uman))/(2.*propac)   )
!u
        Box_US = Box_US + 2*dreal(   -(Coupl_Ac*Coupl_Bc*(8*CB0u + CB3u*m1**2 + CB3u*m2**2 + 2*CB3u*m4**2 + 2*CB4u*m4**2 + 2*CB1u*sqrtS**2 - CB3u*tman + CB2u*(m1**2 + m2**2 + 2*m4**2 - tman - uman) - CB3u*uman))/(2.*propac)  )

      enddo

!================================================================================================================================== 
!Box(t)*Vector + Box(u)*Vector

      do j=sumstartV, sumendV

        Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = dconjg(VVV(j,xtype,ytype))
        propac = dconjg(DenSmV(j))
        
!t        
        Box_TV = Box_TV + 2*dreal(    (Coupl_Ac*Coupl_Bc*(4*CA1t*m1**2*m3**2 + 3*CA4t*m1**2*m3**2 - 4*CA1t*m2**2*m3**2 - 3*CA4t*m2**2*m3**2 - 4*CA1t*m1**2*m4**2 + CA4t*m1**2*m4**2 + 4*CA1t*m2**2*m4**2 - CA4t*m2**2*m4**2 - 
     &      CA4t*m1**2*sqrtS**2 + CA4t*m2**2*sqrtS**2 + 12*CA0t*tman + CA4t*m3**2*tman - CA4t*m4**2*tman + 4*CA1t*sqrtS**2*tman + CA4t*sqrtS**2*tman + 
     &      CA3t*(m1**4 - m2**4 + 2*m1**2*(4*m3**2 - m4**2 - sqrtS**2 + tman - 2*uman) + (4*m3**2 + 2*m4**2 + 2*sqrtS**2 - 3*tman - 3*uman)*(tman - uman) + 
     &         2*m2**2*(-4*m3**2 + m4**2 + sqrtS**2 + 2*tman - uman)) - 12*CA0t*uman - CA4t*m3**2*uman + CA4t*m4**2*uman - 4*CA1t*sqrtS**2*uman - CA4t*sqrtS**2*uman + 
     &      CA2t*(-m1**4 + m2**4 + 2*m1**2*(m4**2 + tman) + (2*m3**2 - 2*sqrtS**2 - tman - uman)*(tman - uman) - 2*m2**2*(m4**2 + uman))))/(4.*propac)  )
!u
        Box_UV = Box_UV + 2*dreal(   (Coupl_Ac*Coupl_Bc*(4*CB1u*m1**2*m3**2 - CB4u*m1**2*m3**2 - 4*CB1u*m2**2*m3**2 + CB4u*m2**2*m3**2 - 4*CB1u*m1**2*m4**2 - 3*CB4u*m1**2*m4**2 + 4*CB1u*m2**2*m4**2 + 3*CB4u*m2**2*m4**2 + 
     &      CB4u*m1**2*sqrtS**2 - CB4u*m2**2*sqrtS**2 + 12*CB0u*tman - CB4u*m3**2*tman + CB4u*m4**2*tman + 4*CB1u*sqrtS**2*tman + CB4u*sqrtS**2*tman - 12*CB0u*uman + CB4u*m3**2*uman - 
     &      CB4u*m4**2*uman - 4*CB1u*sqrtS**2*uman - CB4u*sqrtS**2*uman + CB2u*
     &       (m1**4 - m2**4 + 2*m2**2*(m3**2 + tman) + (2*m4**2 - 2*sqrtS**2 - tman - uman)*(tman - uman) - 2*m1**2*(m3**2 + uman)) + 
     &      CB3u*(-m1**4 + m2**4 + (2*m3**2 + 4*m4**2 + 2*sqrtS**2 - 3*tman - 3*uman)*(tman - uman) + 2*m1**2*(m3**2 - 4*m4**2 + sqrtS**2 + 2*tman - uman) - 
     &         2*m2**2*(m3**2 - 4*m4**2 + sqrtS**2 - tman + 2*uman))))/(4.*propac)  )

      enddo           

!================================================================================================================================== 

      VV_TX(1) = Box_TQ
      VV_TX(2) = Box_TS
      VV_TX(3) = Box_TV
      VV_TX(4) = Box_TT
      VV_TX(5) = Box_TU

      VV_TX(6) = VV_TX(1) + VV_TX(2) + VV_TX(3) + VV_TX(4) + VV_TX(5)

      VV_UX(1) = Box_UQ
      VV_UX(2) = Box_US
      VV_UX(3) = Box_UV
      VV_UX(4) = Box_UT
      VV_UX(5) = Box_UU

      VV_UX(6) = VV_UX(1) + VV_UX(2) + VV_UX(3) + VV_UX(4) + VV_UX(5)

!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************




      subroutine QQ2xx_M2GluonsQuarkLoopBoxVV(argx,argy,VV_TX,VV_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer j,argx,argy
      complex*16 CAt,CAu,QQ2xx_BSVV
      double precision VV_TX(6),VV_UX(6),Box_TQ,Box_TS,Box_TV,Box_TT,Box_TU,Box_UQ,Box_US,Box_UV,Box_UT,Box_UU

! =================================================================================================================================
!                             ---                     Higgs-Vector-final-states                    ---
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

!========================================!Box-Selfenergy(t)_Gluon+Box-Selfenergy(u)_Gluon==========================================

!t- and u-channel both proportional to the metric=> no kinematical difference between t- and u-channel

      CAt= 1/(16d0*pi**2)*QQ2xx_BSVV(0,tman,Looptype0,xtype,ytype,Looptypet,Loopgent,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2)

      CAu= 1/(16d0*pi**2)*QQ2xx_BSVV(1,uman,Looptype0,ytype,xtype,Looptypeu,Loopgenu,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2)
  !==================================================================================================================================

      do j=1,  2   

!Box(t)*t + Box(u)*t
       
        Coupl_Ac = dconjg(VSfSft(xtype,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))       
        Coupl_Bc = dconjg(VSfSfu(ytype,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
        sfmc = MSf(j,proptypetc,propgentc)

!t        
	Box_TT = Box_TT + 2*colorfactorTU*dreal( (CAt*Coupl_Ac*Coupl_Bc*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*uman))/(2.*(sfmc**2 - tman))   )
!u
	Box_UT = Box_UT + 2*colorfactorTU*dreal( (CAu*Coupl_Ac*Coupl_Bc*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*uman))/(2.*(sfmc**2 - tman))   )

!Box(t)*u + Box(u)*u
 
        Coupl_Ac = dconjg(VSfSfu(ytype,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
        Coupl_Bc = dconjg(VSfSft(xtype,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
        sfmc = MSf(j,proptypeuc,propgenuc)

!t        
	Box_TU = Box_TU + 2*colorfactorTU*dreal(  (CAt*Coupl_Ac*Coupl_Bc*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*tman))/(2.*(sfmc**2 - uman))  )
!u
	Box_UU = Box_UU + 2*colorfactorTU*dreal(  (CAu*Coupl_Ac*Coupl_Bc*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*tman))/(2.*(sfmc**2 - uman))  )

      enddo

!================================================================================================================================== 
!Box(t)*q + Box(u)*q

        Coupl_Ac = dconjg(VVSfSf(xtype,ytype,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))

!t
	Box_TQ = Box_TQ + 2*colorfactorQ*dreal(  4*CAt*Coupl_Ac  )
!u
	Box_UQ = Box_UQ + 2*colorfactorQ*dreal(  4*CAu*Coupl_Ac  )

!==================================================================================================================================     
!Box(t)*s(Higgs)
       
      do j=sumstartH, sumendH

        Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(LLR(xtype,ytype,j))
        propac = dconjg(DenSmH(j))

!t        
	Box_TS = Box_TS + 2*colorfactorSH*dreal(  (-4*CAt*Coupl_Ac*Coupl_Bc)/propac  )
!u
	Box_US = Box_US + 2*colorfactorSH*dreal(  (-4*CAu*Coupl_Ac*Coupl_Bc)/propac  )

      enddo

!==================================================================================================================================
!Box(t)*Vector

      do j=sumstartV, sumendV

        Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = dconjg(VVV(j,xtype,ytype))
        propac = dconjg(DenSmV(j))

!t        
        Box_TV = Box_TV + 2*colorfactorSV*dreal(  (3*CAt*Coupl_Ac*Coupl_Bc*(tman - uman))/propac  )
!u
        Box_UV = Box_UV + 2*colorfactorSV*dreal(  (3*CAu*Coupl_Ac*Coupl_Bc*(tman - uman))/propac  )

      enddo           
    
!================================================================================================================================== 

      VV_TX(1) = Box_TQ
      VV_TX(2) = Box_TS
      VV_TX(3) = Box_TV
      VV_TX(4) = Box_TT
      VV_TX(5) = Box_TU

      VV_TX(6) = VV_TX(1) + VV_TX(2) + VV_TX(3) + VV_TX(4) + VV_TX(5)

      VV_UX(1) = Box_UQ
      VV_UX(2) = Box_US
      VV_UX(3) = Box_UV
      VV_UX(4) = Box_UT
      VV_UX(5) = Box_UU

      VV_UX(6) = VV_UX(1) + VV_UX(2) + VV_UX(3) + VV_UX(4) + VV_UX(5)

!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine QQ2xx_M2GluinoBoxVV(argx,argy,VV_TX,VV_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"
#include "../util/DMNLO_LoopIntegrals.h"

      integer j,argx,argy
      complex*16 g0t,g1t,g2t,g3t,g0u,g1u,g2u,g3u,g0Lt,g0Rt,g1Lt,g1Rt,g2Lt,g2Rt,g3Lt,g3Rt,g0Lu,g0Ru,g1Lu,g1Ru,g2Lu,g2Ru,g3Lu,g3Ru,
     &           xb01t,xb02t,xc01t,xc02t,xc03t,xc04t,xc11t,xc12t,xc13t,xc21t,xc22t,xc23t,xc001t,xc002t,xc003t,xc111t,xc112t,xc113t,xc121t,
     &           xc122t,xc123t,xc221t,xc222t,xc223t,xd0t,xd1t,xd2t,xd3t,xd00t,xd11t,xd12t,xd13t,xd22t,xd23t,xd33t,xb01u,xb02u,xc01u,xc02u,
     &           xc03u,xc04u,xc11u,xc12u,xc13u,xc21u,xc22u,xc23u,xc001u,xc002u,xc003u,xc111u,xc112u,xc113u,xc121u,xc122u,xc123u,xc221u,xc222u,
     &           xc223u,xd0u,xd1u,xd2u,xd3u,xd00u,xd11u,xd12u,xd13u,xd22u,xd23u,xd33u
      double precision M0Lpt,M1Lpt,M2Lpt,M3Lpt,M0Lpu,M1Lpu,M2Lpu,M3Lpu,VV_TX(6),VV_UX(6),Box_TQ,Box_TS,Box_TV,Box_TT,Box_TU,Box_UQ,Box_US,Box_UV,Box_UT,Box_UU
        
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
      g1Lt=BLt(xtype,Looptype1,Looptypet,Loopgen1,Loopgent)
      g1Rt=BRt(xtype,Looptype1,Looptypet,Loopgen1,Loopgent)
      g2Lt=BLu(ytype,Looptypet,Looptype2,Loopgent,Loopgen2)
      g2Rt=BRu(ytype,Looptypet,Looptype2,Loopgent,Loopgen2)
      g3Lt=GluL(sfmeigen2,sftype2,Looptype2,sfgen2,Loopgen2)
      g3Rt=GluR(sfmeigen2,sftype2,Looptype2,sfgen2,Loopgen2)

!=====================================================!Box(u)_Gluino================================================================


      M0Lpu=MGl
      M1Lpu=Mf(Looptype1,Loopgen1)
      M2Lpu=Mf(Looptypeu,Loopgenu)
      M3Lpu=Mf(Looptype2,Loopgen2)

      g0Lu=GluR(sfmeigen1,sftype1,Looptype1,sfgen1,Loopgen1)
      g0Ru=GluL(sfmeigen1,sftype1,Looptype1,sfgen1,Loopgen1)
      g1Lu=BLu(ytype,Looptype1,Looptypeu,Loopgen1,Loopgenu)
      g1Ru=BRu(ytype,Looptype1,Looptypeu,Loopgen1,Loopgenu)
      g2Lu=BLt(xtype,Looptypeu,Looptype2,Loopgenu,Loopgen2)
      g2Ru=BRt(xtype,Looptypeu,Looptype2,Loopgenu,Loopgen2)
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
      xc13t = xc1
      xc23t = xc2

      call myCget(m3**2,m4**2,sqrtS**2,M1Lpt**2,M2Lpt**2,M3Lpt**2)

      xc04t = xc0
      xc11t = xc1
      xc21t = xc2
      xc001t = xc00
      xc111t = xc11
      xc121t = xc12
      xc221t = xc22

      call myDget(m1**2, m3**2, m4**2, m2**2, tman, sqrtS**2, M0Lpt**2, M1Lpt**2, M2Lpt**2, M3Lpt**2)

      xd0t = xd0
      xd1t = xd1
      xd2t = xd2
      xd3t = xd3
      xd00t = xd00
      xd11t = xd11
      xd12t = xd12  
      xd13t = xd13     
      xd22t = xd22
      xd23t = xd23 
      xd33t = xd33

!u
      call myBget(m3**2,M2Lpu**2,M3Lpu**2)

      xb01u = xb0

      call myBget(uman,M0Lpu**2,M2Lpu**2)

      xb02u = xb0

      call myCget(m1**2,sqrtS**2,m2**2,M0Lpu**2,M1Lpu**2,M3Lpu**2)
       
      xc01u = xc0    

      call myCget(m1**2,m4**2,uman,M0Lpu**2,M1Lpu**2,M2Lpu**2)
       
      xc02u = xc0    
      xc12u = xc1
      xc22u = xc2

      call myCget(uman,m3**2,m2**2,M0Lpu**2,M2Lpu**2,M3Lpu**2)

      xc03u = xc0
      xc13u = xc1
      xc23u = xc2

      call myCget(m4**2,m3**2,sqrtS**2,M1Lpu**2,M2Lpu**2,M3Lpu**2)

      xc04u = xc0
      xc11u = xc1
      xc21u = xc2
      xc001u = xc00
      xc111u = xc11
      xc121u = xc12
      xc221u = xc22

      call myDget(m1**2, m4**2, m3**2, m2**2, uman, sqrtS**2, M0Lpu**2, M1Lpu**2, M2Lpu**2, M3Lpu**2)

      xd0u = xd0
      xd1u = xd1
      xd2u = xd2
      xd3u = xd3
      xd00u = xd00
      xd11u = xd11
      xd12u = xd12  
      xd13u = xd13     
      xd22u = xd22
      xd23u = xd23 
      xd33u = xd33

!==================================================================================================================================

      do j=1,  2   

!Box(t)*t + Box(u)*t
       
      Coupl_Ac = dconjg(VSfSft(xtype,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))       
      Coupl_Bc = dconjg(VSfSfu(ytype,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
      sfmc = MSf(j,proptypetc,propgentc)
        
!t
      Box_TT = Box_TT + 2*colorfactorTU*dreal(  -(Coupl_Ac*Coupl_Bc*(((6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*uman)*
     &          (g0Rt*(2*g1Rt*M1Lpt*(-(g2Lt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*xd0t + g3Lt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t))) + 
     &                  g2Rt*(g3Rt*M0Lpt*(xc01t + xc02t - 4*xd00t + M2Lpt**2*xd0t + M3Lpt**2*xd0t - m4**2*xd0t) + 
     &                     g3Lt*M3Lpt*(xc01t + xc04t - 4*xd00t + M0Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t))) + 
     &               g1Lt*(-2*g2Rt*M2Lpt*(g3Lt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Rt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &                  g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(xc01t + xc03t - 4*xd00t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - m3**2*xd0t) + 
     &                     g3Lt*(2*xb01t + 2*xb02t - 8*xc001t + 2*M0Lpt**2*xc01t - m1**2*xc01t - m2**2*xc01t + sqrtS**2*xc01t - m1**2*xc02t + 2*M1Lpt**2*xc02t - m3**2*xc02t + tman*xc02t - 
     &                        m2**2*xc03t + 2*M3Lpt**2*xc03t - m4**2*xc03t + tman*xc03t + 2*M1Lpt**2*xc04t + m3**2*xc04t - m4**2*xc04t + sqrtS**2*xc04t + 4*m3**2*xc11t + 2*m3**2*xc21t - 
     &                        2*m4**2*xc21t + 2*sqrtS**2*xc21t - 8*M0Lpt**2*xd00t + 4*m1**2*xd00t + 4*m2**2*xd00t - 4*sqrtS**2*xd00t - M1Lpt**2*m2**2*xd0t + 2*M0Lpt**2*M2Lpt**2*xd0t - 
     &                        m1**2*M2Lpt**2*xd0t - m2**2*M2Lpt**2*xd0t - M0Lpt**2*m3**2*xd0t + m2**2*m3**2*xd0t - m1**2*M3Lpt**2*xd0t + 2*M1Lpt**2*M3Lpt**2*xd0t - m3**2*M3Lpt**2*xd0t - 
     &                        M0Lpt**2*m4**2*xd0t + m1**2*m4**2*xd0t - M1Lpt**2*m4**2*xd0t + M0Lpt**2*sqrtS**2*xd0t + M2Lpt**2*sqrtS**2*xd0t + M1Lpt**2*tman*xd0t + M3Lpt**2*tman*xd0t - 
     &                        sqrtS**2*tman*xd0t)))) + g0Lt*(2*g1Lt*M1Lpt*(-(g2Rt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*xd0t + g3Rt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t))) + 
     &                  g2Lt*(g3Lt*M0Lpt*(xc01t + xc02t - 4*xd00t + M2Lpt**2*xd0t + M3Lpt**2*xd0t - m4**2*xd0t) + 
     &                     g3Rt*M3Lpt*(xc01t + xc04t - 4*xd00t + M0Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t))) + 
     &               g1Rt*(-2*g2Lt*M2Lpt*(g3Rt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Lt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &                  g2Rt*(2*g3Lt*M0Lpt*M3Lpt*(xc01t + xc03t - 4*xd00t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - m3**2*xd0t) + 
     &                     g3Rt*(2*xb01t + 2*xb02t - 8*xc001t + 2*M0Lpt**2*xc01t - m1**2*xc01t - m2**2*xc01t + sqrtS**2*xc01t - m1**2*xc02t + 2*M1Lpt**2*xc02t - m3**2*xc02t + tman*xc02t - 
     &                        m2**2*xc03t + 2*M3Lpt**2*xc03t - m4**2*xc03t + tman*xc03t + 2*M1Lpt**2*xc04t + m3**2*xc04t - m4**2*xc04t + sqrtS**2*xc04t + 4*m3**2*xc11t + 2*m3**2*xc21t - 
     &                        2*m4**2*xc21t + 2*sqrtS**2*xc21t - 8*M0Lpt**2*xd00t + 4*m1**2*xd00t + 4*m2**2*xd00t - 4*sqrtS**2*xd00t - M1Lpt**2*m2**2*xd0t + 2*M0Lpt**2*M2Lpt**2*xd0t - 
     &                        m1**2*M2Lpt**2*xd0t - m2**2*M2Lpt**2*xd0t - M0Lpt**2*m3**2*xd0t + m2**2*m3**2*xd0t - m1**2*M3Lpt**2*xd0t + 2*M1Lpt**2*M3Lpt**2*xd0t - m3**2*M3Lpt**2*xd0t - 
     &                        M0Lpt**2*m4**2*xd0t + m1**2*m4**2*xd0t - M1Lpt**2*m4**2*xd0t + M0Lpt**2*sqrtS**2*xd0t + M2Lpt**2*sqrtS**2*xd0t + M1Lpt**2*tman*xd0t + M3Lpt**2*tman*xd0t - 
     &                        sqrtS**2*tman*xd0t))))))/2. - 2*((m1**2 - tman)*
     &           ((-2*m2**2 - 3*m3**2 - m4**2 + sqrtS**2 + 2*uman)*(g0Rt*(2*g1Rt*g2Rt*M1Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*(xd22t + xd2t) + 
     &                   g1Lt*g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(xd22t + xd2t) + g3Lt*(2*xc111t + 2*xc11t + (2*M0Lpt**2 - m1**2 - m2**2 + sqrtS**2)*(xd22t + xd2t)))) + 
     &                g0Lt*(2*g1Lt*g2Lt*M1Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*(xd22t + xd2t) + 
     &                   g1Rt*g2Rt*(2*g3Lt*M0Lpt*M3Lpt*(xd22t + xd2t) + g3Rt*(2*xc111t + 2*xc11t + (2*M0Lpt**2 - m1**2 - m2**2 + sqrtS**2)*(xd22t + xd2t))))) - 
     &             ((3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*(g0Rt*(2*g1Rt*g2Rt*M1Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*(2*xd12t + xd1t + 2*(xd22t + xd2t)) + 
     &                     g1Lt*(-2*g2Rt*M2Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*xd2t + 
     &                        g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(xd0t + 2*xd12t + xd1t + 2*xd22t + 3*xd2t) + 
     &                           g3Lt*(xc02t - xc04t - 2*xc11t - 4*xc121t - 2*xc21t + 2*xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 4*M0Lpt**2*xd12t - 2*m1**2*xd12t - 
     &                              2*m2**2*xd12t + 2*sqrtS**2*xd12t + 2*M0Lpt**2*xd1t - m1**2*xd1t - m2**2*xd1t + sqrtS**2*xd1t + 4*M0Lpt**2*xd22t - 2*m1**2*xd22t - 2*m2**2*xd22t + 
     &                              2*sqrtS**2*xd22t + 4*M0Lpt**2*xd2t - 2*m1**2*xd2t - 3*m2**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t + 2*sqrtS**2*xd2t + tman*xd2t)))) + 
     &                  g0Lt*(2*g1Lt*g2Lt*M1Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*(2*xd12t + xd1t + 2*(xd22t + xd2t)) + 
     &                     g1Rt*(-2*g2Lt*M2Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*xd2t + 
     &                        g2Rt*(2*g3Lt*M0Lpt*M3Lpt*(xd0t + 2*xd12t + xd1t + 2*xd22t + 3*xd2t) + 
     &                           g3Rt*(xc02t - xc04t - 2*xc11t - 4*xc121t - 2*xc21t + 2*xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 4*M0Lpt**2*xd12t - 2*m1**2*xd12t - 
     &                              2*m2**2*xd12t + 2*sqrtS**2*xd12t + 2*M0Lpt**2*xd1t - m1**2*xd1t - m2**2*xd1t + sqrtS**2*xd1t + 4*M0Lpt**2*xd22t - 2*m1**2*xd22t - 2*m2**2*xd22t + 
     &                              2*sqrtS**2*xd22t + 4*M0Lpt**2*xd2t - 2*m1**2*xd2t - 3*m2**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t + 2*sqrtS**2*xd2t + tman*xd2t))))))/2. + 
     &             ((3*m2**2 - m4**2 + tman)*(-(g0Rt*(2*g1Rt*M1Lpt*(g2Lt*g3Lt*M2Lpt*xd2t + g2Rt*(g3Lt*M3Lpt*(2*xd23t + xd3t) + g3Rt*M0Lpt*(xd0t + 2*xd23t + xd2t + xd3t))) + 
     &                       g1Lt*(2*g2Rt*g3Rt*M0Lpt*M2Lpt*xd2t + g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(2*xd23t + xd3t) + 
     &                             g3Lt*(xc03t + xc04t + 2*xc11t + 4*xc121t + 2*xc21t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 4*M0Lpt**2*xd23t - 2*m1**2*xd23t - 2*m2**2*xd23t + 
     &                                2*sqrtS**2*xd23t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t + 2*M0Lpt**2*xd3t - m1**2*xd3t - m2**2*xd3t + sqrtS**2*xd3t))))) - 
     &                  g0Lt*(2*g1Lt*M1Lpt*(g2Rt*g3Rt*M2Lpt*xd2t + g2Lt*(g3Rt*M3Lpt*(2*xd23t + xd3t) + g3Lt*M0Lpt*(xd0t + 2*xd23t + xd2t + xd3t))) + 
     &                     g1Rt*(2*g2Lt*g3Lt*M0Lpt*M2Lpt*xd2t + g2Rt*(2*g3Lt*M0Lpt*M3Lpt*(2*xd23t + xd3t) + 
     &                           g3Rt*(xc03t + xc04t + 2*xc11t + 4*xc121t + 2*xc21t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 4*M0Lpt**2*xd23t - 2*m1**2*xd23t - 2*m2**2*xd23t + 
     &                              2*sqrtS**2*xd23t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t + 2*M0Lpt**2*xd3t - m1**2*xd3t - m2**2*xd3t + sqrtS**2*xd3t))))))/2.) + 
     &          ((-2*m1**2 - 3*m2**2 - m3**2 + 2*sqrtS**2 + uman)*(2*(-3*m2**2 + m4**2 - tman)*
     &                (g0Rt*(2*g1Rt*g2Rt*M1Lpt*(g3Lt*M3Lpt*xd33t + g3Rt*M0Lpt*(xd33t + xd3t)) + 
     &                     g1Lt*g2Lt*(2*g3Rt*M0Lpt*M3Lpt*xd33t + g3Lt*(xc21t + 2*xc221t + xc23t + 2*M0Lpt**2*xd33t - m1**2*xd33t - m2**2*xd33t + sqrtS**2*xd33t + M0Lpt**2*xd3t - m1**2*xd3t + 
     &                           M1Lpt**2*xd3t))) + g0Lt*(2*g1Lt*g2Lt*M1Lpt*(g3Rt*M3Lpt*xd33t + g3Lt*M0Lpt*(xd33t + xd3t)) + 
     &                     g1Rt*g2Rt*(2*g3Lt*M0Lpt*M3Lpt*xd33t + g3Rt*(xc21t + 2*xc221t + xc23t + 2*M0Lpt**2*xd33t - m1**2*xd33t - m2**2*xd33t + sqrtS**2*xd33t + M0Lpt**2*xd3t - m1**2*xd3t + 
     &                           M1Lpt**2*xd3t)))) + (2*m2**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*uman)*
     &                (-(g0Rt*(2*g1Rt*M1Lpt*(-(g2Lt*g3Lt*M2Lpt*xd2t) + g2Rt*(g3Lt*M3Lpt*(2*xd23t + xd3t) + g3Rt*M0Lpt*(xd0t + 2*xd23t + xd2t + xd3t))) + 
     &                       g1Lt*(-2*g2Rt*g3Rt*M0Lpt*M2Lpt*xd2t + g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(2*xd23t + xd3t) + 
     &                             g3Lt*(xc03t + xc04t + 4*xc121t + 2*xc13t + 2*xc21t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 4*M0Lpt**2*xd23t - 2*m1**2*xd23t - 2*m2**2*xd23t + 
     &                                2*sqrtS**2*xd23t - m1**2*xd2t + 2*M1Lpt**2*xd2t - m3**2*xd2t + tman*xd2t + 2*M0Lpt**2*xd3t - m1**2*xd3t - m2**2*xd3t + sqrtS**2*xd3t))))) - 
     &                  g0Lt*(2*g1Lt*M1Lpt*(-(g2Rt*g3Rt*M2Lpt*xd2t) + g2Lt*(g3Rt*M3Lpt*(2*xd23t + xd3t) + g3Lt*M0Lpt*(xd0t + 2*xd23t + xd2t + xd3t))) + 
     &                     g1Rt*(-2*g2Lt*g3Lt*M0Lpt*M2Lpt*xd2t + g2Rt*(2*g3Lt*M0Lpt*M3Lpt*(2*xd23t + xd3t) + 
     &                           g3Rt*(xc03t + xc04t + 4*xc121t + 2*xc13t + 2*xc21t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 4*M0Lpt**2*xd23t - 2*m1**2*xd23t - 2*m2**2*xd23t + 
     &                              2*sqrtS**2*xd23t - m1**2*xd2t + 2*M1Lpt**2*xd2t - m3**2*xd2t + tman*xd2t + 2*M0Lpt**2*xd3t - m1**2*xd3t - m2**2*xd3t + sqrtS**2*xd3t))))) + 
     &               (3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*(g0Rt*(-2*g1Rt*M1Lpt*
     &                      (-(g2Lt*g3Lt*M2Lpt*(xd1t + xd2t)) + g2Rt*(g3Lt*M3Lpt*(2*xd13t + 2*xd23t + xd3t) + g3Rt*M0Lpt*(xd0t + 2*xd13t + xd1t + 2*xd23t + xd2t + xd3t))) + 
     &                     g1Lt*(2*g2Rt*M2Lpt*(g3Lt*M3Lpt*xd3t + g3Rt*M0Lpt*(xd0t + xd1t + xd2t + xd3t)) + 
     &                        g2Lt*(-4*g3Rt*M0Lpt*M3Lpt*(xd13t + xd23t + xd3t) + 
     &                           g3Lt*(xc01t - xc03t - 2*xc13t + 2*xc21t + 4*xc221t + m1**2*xd0t - M1Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t - 4*M0Lpt**2*xd13t + 2*m1**2*xd13t + 
     &                              2*m2**2*xd13t - 2*sqrtS**2*xd13t + m1**2*xd1t - 2*M1Lpt**2*xd1t + m3**2*xd1t - tman*xd1t - 4*M0Lpt**2*xd23t + 2*m1**2*xd23t + 2*m2**2*xd23t - 
     &                              2*sqrtS**2*xd23t + m1**2*xd2t - 2*M1Lpt**2*xd2t + m3**2*xd2t - tman*xd2t - 2*M0Lpt**2*xd3t + m1**2*xd3t + 2*m2**2*xd3t - 2*M3Lpt**2*xd3t + m4**2*xd3t - 
     &                              sqrtS**2*xd3t - tman*xd3t)))) + g0Lt*(-2*g1Lt*M1Lpt*
     &                      (-(g2Rt*g3Rt*M2Lpt*(xd1t + xd2t)) + g2Lt*(g3Rt*M3Lpt*(2*xd13t + 2*xd23t + xd3t) + g3Lt*M0Lpt*(xd0t + 2*xd13t + xd1t + 2*xd23t + xd2t + xd3t))) + 
     &                     g1Rt*(2*g2Lt*M2Lpt*(g3Rt*M3Lpt*xd3t + g3Lt*M0Lpt*(xd0t + xd1t + xd2t + xd3t)) + 
     &                        g2Rt*(-4*g3Lt*M0Lpt*M3Lpt*(xd13t + xd23t + xd3t) + 
     &                           g3Rt*(xc01t - xc03t - 2*xc13t + 2*xc21t + 4*xc221t + m1**2*xd0t - M1Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t - 4*M0Lpt**2*xd13t + 2*m1**2*xd13t + 
     &                              2*m2**2*xd13t - 2*sqrtS**2*xd13t + m1**2*xd1t - 2*M1Lpt**2*xd1t + m3**2*xd1t - tman*xd1t - 4*M0Lpt**2*xd23t + 2*m1**2*xd23t + 2*m2**2*xd23t - 
     &                              2*sqrtS**2*xd23t + m1**2*xd2t - 2*M1Lpt**2*xd2t + m3**2*xd2t - tman*xd2t - 2*M0Lpt**2*xd3t + m1**2*xd3t + 2*m2**2*xd3t - 2*M3Lpt**2*xd3t + m4**2*xd3t - 
     &                              sqrtS**2*xd3t - tman*xd3t)))))))/4. - ((3*m1**2 - m3**2 + tman)*
     &             (-2*(3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*(g0Rt*
     &                   (2*g1Rt*g2Rt*M1Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*(xd11t + 2*xd12t + xd1t + xd22t + xd2t) + 
     &                     g1Lt*g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(xd0t + xd11t + 2*xd12t + 2*xd1t + xd22t + 2*xd2t) + 
     &                        g3Lt*(xc02t + xc12t + xc21t + 2*xc221t + xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 2*M0Lpt**2*xd11t - m1**2*xd11t - m2**2*xd11t + sqrtS**2*xd11t + 
     &                           4*M0Lpt**2*xd12t - 2*m1**2*xd12t - 2*m2**2*xd12t + 2*sqrtS**2*xd12t + 3*M0Lpt**2*xd1t - m1**2*xd1t - 2*m2**2*xd1t + M3Lpt**2*xd1t + sqrtS**2*xd1t + 
     &                           2*M0Lpt**2*xd22t - m1**2*xd22t - m2**2*xd22t + sqrtS**2*xd22t + 3*M0Lpt**2*xd2t - m1**2*xd2t - 2*m2**2*xd2t + M3Lpt**2*xd2t + sqrtS**2*xd2t))) + 
     &                  g0Lt*(2*g1Lt*g2Lt*M1Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*(xd11t + 2*xd12t + xd1t + xd22t + xd2t) + 
     &                     g1Rt*g2Rt*(2*g3Lt*M0Lpt*M3Lpt*(xd0t + xd11t + 2*xd12t + 2*xd1t + xd22t + 2*xd2t) + 
     &                        g3Rt*(xc02t + xc12t + xc21t + 2*xc221t + xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 2*M0Lpt**2*xd11t - m1**2*xd11t - m2**2*xd11t + sqrtS**2*xd11t + 
     &                           4*M0Lpt**2*xd12t - 2*m1**2*xd12t - 2*m2**2*xd12t + 2*sqrtS**2*xd12t + 3*M0Lpt**2*xd1t - m1**2*xd1t - 2*m2**2*xd1t + M3Lpt**2*xd1t + sqrtS**2*xd1t + 
     &                           2*M0Lpt**2*xd22t - m1**2*xd22t - m2**2*xd22t + sqrtS**2*xd22t + 3*M0Lpt**2*xd2t - m1**2*xd2t - 2*m2**2*xd2t + M3Lpt**2*xd2t + sqrtS**2*xd2t)))) + 
     &               (2*m2**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*uman)*(-(g0Rt*
     &                     (2*g1Rt*g2Rt*M1Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*(2*xd12t + xd1t + 2*(xd22t + xd2t)) + 
     &                       g1Lt*(2*g2Rt*M2Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*xd2t + 
     &                          g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(xd0t + 2*xd12t + xd1t + 2*xd22t + 3*xd2t) + 
     &                             g3Lt*(xc02t - xc04t - 4*xc121t - 2*xc21t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 4*M0Lpt**2*xd12t - 2*m1**2*xd12t - 2*m2**2*xd12t + 
     &                                2*sqrtS**2*xd12t + 2*M0Lpt**2*xd1t - m1**2*xd1t - m2**2*xd1t + sqrtS**2*xd1t + 4*M0Lpt**2*xd22t - 2*m1**2*xd22t - 2*m2**2*xd22t + 2*sqrtS**2*xd22t + 
     &                                6*M0Lpt**2*xd2t - 2*m1**2*xd2t - 3*m2**2*xd2t + m4**2*xd2t + 2*sqrtS**2*xd2t - tman*xd2t))))) - 
     &                  g0Lt*(2*g1Lt*g2Lt*M1Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*(2*xd12t + xd1t + 2*(xd22t + xd2t)) + 
     &                     g1Rt*(2*g2Lt*M2Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*xd2t + 
     &                        g2Rt*(2*g3Lt*M0Lpt*M3Lpt*(xd0t + 2*xd12t + xd1t + 2*xd22t + 3*xd2t) + 
     &                           g3Rt*(xc02t - xc04t - 4*xc121t - 2*xc21t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 4*M0Lpt**2*xd12t - 2*m1**2*xd12t - 2*m2**2*xd12t + 2*sqrtS**2*xd12t + 
     &                              2*M0Lpt**2*xd1t - m1**2*xd1t - m2**2*xd1t + sqrtS**2*xd1t + 4*M0Lpt**2*xd22t - 2*m1**2*xd22t - 2*m2**2*xd22t + 2*sqrtS**2*xd22t + 6*M0Lpt**2*xd2t - 
     &                              2*m1**2*xd2t - 3*m2**2*xd2t + m4**2*xd2t + 2*sqrtS**2*xd2t - tman*xd2t))))) + 
     &               (3*m2**2 - m4**2 + tman)*(-(g0Rt*(2*g1Rt*M1Lpt*(g2Lt*g3Lt*M2Lpt*(xd1t + xd2t) + 
     &                          g2Rt*(g3Lt*M3Lpt*(2*xd13t + 2*xd23t + xd3t) + g3Rt*M0Lpt*(xd0t + 2*xd13t + xd1t + 2*xd23t + xd2t + xd3t))) + 
     &                       g1Lt*(2*g2Rt*M2Lpt*(g3Lt*M3Lpt*xd3t + g3Rt*M0Lpt*(xd0t + xd1t + xd2t + xd3t)) + 
     &                          g2Lt*(4*g3Rt*M0Lpt*M3Lpt*(xd13t + xd23t + xd3t) + 
     &                             g3Lt*(xc01t + xc03t - 2*xc21t - 4*xc221t + 2*M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t + 4*M0Lpt**2*xd13t - 2*m1**2*xd13t - 
     &                                2*m2**2*xd13t + 2*sqrtS**2*xd13t + 2*M0Lpt**2*xd1t - m1**2*xd1t + m3**2*xd1t - tman*xd1t + 4*M0Lpt**2*xd23t - 2*m1**2*xd23t - 2*m2**2*xd23t + 
     &                                2*sqrtS**2*xd23t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t + 4*M0Lpt**2*xd3t - m1**2*xd3t - 2*m2**2*xd3t + m4**2*xd3t + sqrtS**2*xd3t - 
     &                                tman*xd3t))))) - g0Lt*(2*g1Lt*M1Lpt*(g2Rt*g3Rt*M2Lpt*(xd1t + xd2t) + 
     &                        g2Lt*(g3Rt*M3Lpt*(2*xd13t + 2*xd23t + xd3t) + g3Lt*M0Lpt*(xd0t + 2*xd13t + xd1t + 2*xd23t + xd2t + xd3t))) + 
     &                     g1Rt*(2*g2Lt*M2Lpt*(g3Rt*M3Lpt*xd3t + g3Lt*M0Lpt*(xd0t + xd1t + xd2t + xd3t)) + 
     &                        g2Rt*(4*g3Lt*M0Lpt*M3Lpt*(xd13t + xd23t + xd3t) + 
     &                           g3Rt*(xc01t + xc03t - 2*xc21t - 4*xc221t + 2*M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t + 4*M0Lpt**2*xd13t - 2*m1**2*xd13t - 
     &                              2*m2**2*xd13t + 2*sqrtS**2*xd13t + 2*M0Lpt**2*xd1t - m1**2*xd1t + m3**2*xd1t - tman*xd1t + 4*M0Lpt**2*xd23t - 2*m1**2*xd23t - 2*m2**2*xd23t + 
     &                              2*sqrtS**2*xd23t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t + 4*M0Lpt**2*xd3t - m1**2*xd3t - 2*m2**2*xd3t + m4**2*xd3t + sqrtS**2*xd3t - 
     &                              tman*xd3t)))))))/4.)))/(32.*pi**2*(-sfmc**2 + tman))  )

!u
      Box_UT = Box_UT + 2*colorfactorTU*dreal(  -(Coupl_Ac*Coupl_Bc*((6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*uman)*
     &        (g0Ru*(2*g1Ru*M1Lpu*(-(g2Lu*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*xd0u + g3Lu*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u))) + 
     &                g2Ru*(g3Ru*M0Lpu*(xc01u + xc02u - 4*xd00u + M2Lpu**2*xd0u - m3**2*xd0u + M3Lpu**2*xd0u) + 
     &                   g3Lu*M3Lpu*(xc01u + xc04u - 4*xd00u + M0Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u))) + 
     &             g1Lu*(-2*g2Ru*M2Lpu*(g3Lu*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Ru*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &                g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(xc01u + xc03u - 4*xd00u + M1Lpu**2*xd0u + M2Lpu**2*xd0u - m4**2*xd0u) + 
     &                   g3Lu*(2*xb01u + 2*xb02u - 8*xc001u + 2*M0Lpu**2*xc01u - m1**2*xc01u - m2**2*xc01u + sqrtS**2*xc01u - m1**2*xc02u + 2*M1Lpu**2*xc02u - m4**2*xc02u + uman*xc02u - 
     &                      m2**2*xc03u - m3**2*xc03u + 2*M3Lpu**2*xc03u + uman*xc03u + 2*M1Lpu**2*xc04u - m3**2*xc04u + m4**2*xc04u + sqrtS**2*xc04u + 4*m4**2*xc11u - 2*m3**2*xc21u + 
     &                      2*m4**2*xc21u + 2*sqrtS**2*xc21u - 8*M0Lpu**2*xd00u + 4*m1**2*xd00u + 4*m2**2*xd00u - 4*sqrtS**2*xd00u - M1Lpu**2*m2**2*xd0u + 2*M0Lpu**2*M2Lpu**2*xd0u - 
     &                      m1**2*M2Lpu**2*xd0u - m2**2*M2Lpu**2*xd0u - M0Lpu**2*m3**2*xd0u + m1**2*m3**2*xd0u - M1Lpu**2*m3**2*xd0u - m1**2*M3Lpu**2*xd0u + 2*M1Lpu**2*M3Lpu**2*xd0u - 
     &                      M0Lpu**2*m4**2*xd0u + m2**2*m4**2*xd0u - M3Lpu**2*m4**2*xd0u + M0Lpu**2*sqrtS**2*xd0u + M2Lpu**2*sqrtS**2*xd0u + M1Lpu**2*uman*xd0u + M3Lpu**2*uman*xd0u - 
     &                      sqrtS**2*uman*xd0u)))) + g0Lu*(2*g1Lu*M1Lpu*(-(g2Ru*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*xd0u + g3Ru*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u))) + 
     &                g2Lu*(g3Lu*M0Lpu*(xc01u + xc02u - 4*xd00u + M2Lpu**2*xd0u - m3**2*xd0u + M3Lpu**2*xd0u) + 
     &                   g3Ru*M3Lpu*(xc01u + xc04u - 4*xd00u + M0Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u))) + 
     &             g1Ru*(-2*g2Lu*M2Lpu*(g3Ru*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Lu*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &                g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(xc01u + xc03u - 4*xd00u + M1Lpu**2*xd0u + M2Lpu**2*xd0u - m4**2*xd0u) + 
     &                   g3Ru*(2*xb01u + 2*xb02u - 8*xc001u + 2*M0Lpu**2*xc01u - m1**2*xc01u - m2**2*xc01u + sqrtS**2*xc01u - m1**2*xc02u + 2*M1Lpu**2*xc02u - m4**2*xc02u + uman*xc02u - 
     &                      m2**2*xc03u - m3**2*xc03u + 2*M3Lpu**2*xc03u + uman*xc03u + 2*M1Lpu**2*xc04u - m3**2*xc04u + m4**2*xc04u + sqrtS**2*xc04u + 4*m4**2*xc11u - 2*m3**2*xc21u + 
     &                      2*m4**2*xc21u + 2*sqrtS**2*xc21u - 8*M0Lpu**2*xd00u + 4*m1**2*xd00u + 4*m2**2*xd00u - 4*sqrtS**2*xd00u - M1Lpu**2*m2**2*xd0u + 2*M0Lpu**2*M2Lpu**2*xd0u - 
     &                      m1**2*M2Lpu**2*xd0u - m2**2*M2Lpu**2*xd0u - M0Lpu**2*m3**2*xd0u + m1**2*m3**2*xd0u - M1Lpu**2*m3**2*xd0u - m1**2*M3Lpu**2*xd0u + 2*M1Lpu**2*M3Lpu**2*xd0u - 
     &                      M0Lpu**2*m4**2*xd0u + m2**2*m4**2*xd0u - M3Lpu**2*m4**2*xd0u + M0Lpu**2*sqrtS**2*xd0u + M2Lpu**2*sqrtS**2*xd0u + M1Lpu**2*uman*xd0u + M3Lpu**2*uman*xd0u - 
     &                      sqrtS**2*uman*xd0u))))) - (2*m1**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*uman)*
     &        (4*(-m2**2 + tman)*(g0Ru*(2*g1Ru*g2Ru*M1Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*(xd22u + xd2u) + 
     &                g1Lu*g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(xd22u + xd2u) + g3Lu*(2*xc111u + 2*xc11u + (2*M0Lpu**2 - m1**2 - m2**2 + sqrtS**2)*(xd22u + xd2u)))) + 
     &             g0Lu*(2*g1Lu*g2Lu*M1Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*(xd22u + xd2u) + 
     &                g1Ru*g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(xd22u + xd2u) + g3Ru*(2*xc111u + 2*xc11u + (2*M0Lpu**2 - m1**2 - m2**2 + sqrtS**2)*(xd22u + xd2u))))) - 
     &          (3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*(g0Ru*(2*g1Ru*g2Ru*M1Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*(2*xd12u + xd1u + 2*(xd22u + xd2u)) + 
     &                g1Lu*(2*g2Ru*M2Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*xd2u + g2Lu*
     &                    (2*g3Ru*M0Lpu*M3Lpu*(xd0u + 2*xd12u + xd1u + 2*xd22u + 3*xd2u) + 
     &                      g3Lu*(xc02u - xc04u - 4*xc121u - 2*xc21u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 4*M0Lpu**2*xd12u - 2*m1**2*xd12u - 2*m2**2*xd12u + 2*sqrtS**2*xd12u + 
     &                         2*M0Lpu**2*xd1u - m1**2*xd1u - m2**2*xd1u + sqrtS**2*xd1u + 4*M0Lpu**2*xd22u - 2*m1**2*xd22u - 2*m2**2*xd22u + 2*sqrtS**2*xd22u + 6*M0Lpu**2*xd2u - 
     &                         2*m1**2*xd2u - 3*m2**2*xd2u + m3**2*xd2u + 2*sqrtS**2*xd2u - uman*xd2u)))) + 
     &             g0Lu*(2*g1Lu*g2Lu*M1Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*(2*xd12u + xd1u + 2*(xd22u + xd2u)) + 
     &                g1Ru*(2*g2Lu*M2Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*xd2u + g2Ru*
     &                    (2*g3Lu*M0Lpu*M3Lpu*(xd0u + 2*xd12u + xd1u + 2*xd22u + 3*xd2u) + 
     &                      g3Ru*(xc02u - xc04u - 4*xc121u - 2*xc21u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 4*M0Lpu**2*xd12u - 2*m1**2*xd12u - 2*m2**2*xd12u + 2*sqrtS**2*xd12u + 
     &                         2*M0Lpu**2*xd1u - m1**2*xd1u - m2**2*xd1u + sqrtS**2*xd1u + 4*M0Lpu**2*xd22u - 2*m1**2*xd22u - 2*m2**2*xd22u + 2*sqrtS**2*xd22u + 6*M0Lpu**2*xd2u - 
     &                         2*m1**2*xd2u - 3*m2**2*xd2u + m3**2*xd2u + 2*sqrtS**2*xd2u - uman*xd2u))))) + 
     &          (3*m2**2 - m4**2 + tman)*(-(g0Ru*(2*g1Ru*M1Lpu*(-(g2Lu*g3Lu*M2Lpu*xd2u) + g2Ru*(g3Lu*M3Lpu*(2*xd23u + xd3u) + g3Ru*M0Lpu*(xd0u + 2*xd23u + xd2u + xd3u))) + 
     &                  g1Lu*(-2*g2Ru*g3Ru*M0Lpu*M2Lpu*xd2u + g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(2*xd23u + xd3u) + 
     &                        g3Lu*(xc03u + xc04u + 4*xc121u + 2*xc13u + 2*xc21u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + 4*M0Lpu**2*xd23u - 2*m1**2*xd23u - 2*m2**2*xd23u + 
     &                           2*sqrtS**2*xd23u - m1**2*xd2u + 2*M1Lpu**2*xd2u - m4**2*xd2u + uman*xd2u + 2*M0Lpu**2*xd3u - m1**2*xd3u - m2**2*xd3u + sqrtS**2*xd3u))))) - 
     &             g0Lu*(2*g1Lu*M1Lpu*(-(g2Ru*g3Ru*M2Lpu*xd2u) + g2Lu*(g3Ru*M3Lpu*(2*xd23u + xd3u) + g3Lu*M0Lpu*(xd0u + 2*xd23u + xd2u + xd3u))) + 
     &                g1Ru*(-2*g2Lu*g3Lu*M0Lpu*M2Lpu*xd2u + g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(2*xd23u + xd3u) + 
     &                      g3Ru*(xc03u + xc04u + 4*xc121u + 2*xc13u + 2*xc21u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + 4*M0Lpu**2*xd23u - 2*m1**2*xd23u - 2*m2**2*xd23u + 
     &                         2*sqrtS**2*xd23u - m1**2*xd2u + 2*M1Lpu**2*xd2u - m4**2*xd2u + uman*xd2u + 2*M0Lpu**2*xd3u - m1**2*xd3u - m2**2*xd3u + sqrtS**2*xd3u)))))) - 
     &       2*(3*m1**2 - m3**2 + tman)*((3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*
     &           (g0Ru*(2*g1Ru*g2Ru*M1Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*(xd11u + 2*xd12u + xd1u + xd22u + xd2u) + 
     &                g1Lu*g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(xd0u + xd11u + 2*xd12u + 2*xd1u + xd22u + 2*xd2u) + 
     &                   g3Lu*(xc02u + xc12u + xc21u + 2*xc221u + xc22u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd11u - m1**2*xd11u - m2**2*xd11u + sqrtS**2*xd11u + 
     &                      4*M0Lpu**2*xd12u - 2*m1**2*xd12u - 2*m2**2*xd12u + 2*sqrtS**2*xd12u + 3*M0Lpu**2*xd1u - m1**2*xd1u - 2*m2**2*xd1u + M3Lpu**2*xd1u + sqrtS**2*xd1u + 
     &                      2*M0Lpu**2*xd22u - m1**2*xd22u - m2**2*xd22u + sqrtS**2*xd22u + 3*M0Lpu**2*xd2u - m1**2*xd2u - 2*m2**2*xd2u + M3Lpu**2*xd2u + sqrtS**2*xd2u))) + 
     &             g0Lu*(2*g1Lu*g2Lu*M1Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*(xd11u + 2*xd12u + xd1u + xd22u + xd2u) + 
     &                g1Ru*g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(xd0u + xd11u + 2*xd12u + 2*xd1u + xd22u + 2*xd2u) + 
     &                   g3Ru*(xc02u + xc12u + xc21u + 2*xc221u + xc22u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd11u - m1**2*xd11u - m2**2*xd11u + sqrtS**2*xd11u + 
     &                      4*M0Lpu**2*xd12u - 2*m1**2*xd12u - 2*m2**2*xd12u + 2*sqrtS**2*xd12u + 3*M0Lpu**2*xd1u - m1**2*xd1u - 2*m2**2*xd1u + M3Lpu**2*xd1u + sqrtS**2*xd1u + 
     &                      2*M0Lpu**2*xd22u - m1**2*xd22u - m2**2*xd22u + sqrtS**2*xd22u + 3*M0Lpu**2*xd2u - m1**2*xd2u - 2*m2**2*xd2u + M3Lpu**2*xd2u + sqrtS**2*xd2u)))) - 
     &          (-m2**2 + tman)*(g0Ru*(2*g1Ru*g2Ru*M1Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*(2*xd12u + xd1u + 2*(xd22u + xd2u)) + 
     &                g1Lu*(-2*g2Ru*M2Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*xd2u + g2Lu*
     &                    (2*g3Ru*M0Lpu*M3Lpu*(xd0u + 2*xd12u + xd1u + 2*xd22u + 3*xd2u) + 
     &                      g3Lu*(xc02u - xc04u - 2*xc11u - 4*xc121u - 2*xc21u + 2*xc22u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 4*M0Lpu**2*xd12u - 2*m1**2*xd12u - 2*m2**2*xd12u + 
     &                         2*sqrtS**2*xd12u + 2*M0Lpu**2*xd1u - m1**2*xd1u - m2**2*xd1u + sqrtS**2*xd1u + 4*M0Lpu**2*xd22u - 2*m1**2*xd22u - 2*m2**2*xd22u + 2*sqrtS**2*xd22u + 
     &                         4*M0Lpu**2*xd2u - 2*m1**2*xd2u - 3*m2**2*xd2u - m3**2*xd2u + 2*M3Lpu**2*xd2u + 2*sqrtS**2*xd2u + uman*xd2u)))) + 
     &             g0Lu*(2*g1Lu*g2Lu*M1Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*(2*xd12u + xd1u + 2*(xd22u + xd2u)) + 
     &                g1Ru*(-2*g2Lu*M2Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*xd2u + g2Ru*
     &                    (2*g3Lu*M0Lpu*M3Lpu*(xd0u + 2*xd12u + xd1u + 2*xd22u + 3*xd2u) + 
     &                      g3Ru*(xc02u - xc04u - 2*xc11u - 4*xc121u - 2*xc21u + 2*xc22u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 4*M0Lpu**2*xd12u - 2*m1**2*xd12u - 2*m2**2*xd12u + 
     &                         2*sqrtS**2*xd12u + 2*M0Lpu**2*xd1u - m1**2*xd1u - m2**2*xd1u + sqrtS**2*xd1u + 4*M0Lpu**2*xd22u - 2*m1**2*xd22u - 2*m2**2*xd22u + 2*sqrtS**2*xd22u + 
     &                         4*M0Lpu**2*xd2u - 2*m1**2*xd2u - 3*m2**2*xd2u - m3**2*xd2u + 2*M3Lpu**2*xd2u + 2*sqrtS**2*xd2u + uman*xd2u))))) - 
     &          ((3*m2**2 - m4**2 + tman)*(g0Ru*(-2*g1Ru*M1Lpu*(-(g2Lu*g3Lu*M2Lpu*(xd1u + xd2u)) + 
     &                     g2Ru*(g3Lu*M3Lpu*(2*xd13u + 2*xd23u + xd3u) + g3Ru*M0Lpu*(xd0u + 2*xd13u + xd1u + 2*xd23u + xd2u + xd3u))) + 
     &                  g1Lu*(2*g2Ru*M2Lpu*(g3Lu*M3Lpu*xd3u + g3Ru*M0Lpu*(xd0u + xd1u + xd2u + xd3u)) + 
     &                     g2Lu*(-4*g3Ru*M0Lpu*M3Lpu*(xd13u + xd23u + xd3u) + 
     &                        g3Lu*(xc01u - xc03u - 2*xc13u + 2*xc21u + 4*xc221u + m1**2*xd0u - M1Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u - 4*M0Lpu**2*xd13u + 2*m1**2*xd13u + 
     &                           2*m2**2*xd13u - 2*sqrtS**2*xd13u + m1**2*xd1u - 2*M1Lpu**2*xd1u + m4**2*xd1u - uman*xd1u - 4*M0Lpu**2*xd23u + 2*m1**2*xd23u + 2*m2**2*xd23u - 
     &                           2*sqrtS**2*xd23u + m1**2*xd2u - 2*M1Lpu**2*xd2u + m4**2*xd2u - uman*xd2u - 2*M0Lpu**2*xd3u + m1**2*xd3u + 2*m2**2*xd3u + m3**2*xd3u - 2*M3Lpu**2*xd3u - 
     &                           sqrtS**2*xd3u - uman*xd3u)))) + g0Lu*(-2*g1Lu*M1Lpu*
     &                   (-(g2Ru*g3Ru*M2Lpu*(xd1u + xd2u)) + g2Lu*(g3Ru*M3Lpu*(2*xd13u + 2*xd23u + xd3u) + g3Lu*M0Lpu*(xd0u + 2*xd13u + xd1u + 2*xd23u + xd2u + xd3u))) + 
     &                  g1Ru*(2*g2Lu*M2Lpu*(g3Ru*M3Lpu*xd3u + g3Lu*M0Lpu*(xd0u + xd1u + xd2u + xd3u)) + 
     &                     g2Ru*(-4*g3Lu*M0Lpu*M3Lpu*(xd13u + xd23u + xd3u) + 
     &                        g3Ru*(xc01u - xc03u - 2*xc13u + 2*xc21u + 4*xc221u + m1**2*xd0u - M1Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u - 4*M0Lpu**2*xd13u + 2*m1**2*xd13u + 
     &                           2*m2**2*xd13u - 2*sqrtS**2*xd13u + m1**2*xd1u - 2*M1Lpu**2*xd1u + m4**2*xd1u - uman*xd1u - 4*M0Lpu**2*xd23u + 2*m1**2*xd23u + 2*m2**2*xd23u - 
     &                           2*sqrtS**2*xd23u + m1**2*xd2u - 2*M1Lpu**2*xd2u + m4**2*xd2u - uman*xd2u - 2*M0Lpu**2*xd3u + m1**2*xd3u + 2*m2**2*xd3u + m3**2*xd3u - 2*M3Lpu**2*xd3u - 
     &                           sqrtS**2*xd3u - uman*xd3u))))))/2.) - 2*(-2*m1**2 - 3*m2**2 - m3**2 + 2*sqrtS**2 + uman)*
     &        ((-3*m2**2 + m4**2 - tman)*(g0Ru*(2*g1Ru*g2Ru*M1Lpu*(g3Lu*M3Lpu*xd33u + g3Ru*M0Lpu*(xd33u + xd3u)) + 
     &                g1Lu*g2Lu*(2*g3Ru*M0Lpu*M3Lpu*xd33u + g3Lu*(xc21u + 2*xc221u + xc23u + 2*M0Lpu**2*xd33u - m1**2*xd33u - m2**2*xd33u + sqrtS**2*xd33u + M0Lpu**2*xd3u - m1**2*xd3u + 
     &                      M1Lpu**2*xd3u))) + g0Lu*(2*g1Lu*g2Lu*M1Lpu*(g3Ru*M3Lpu*xd33u + g3Lu*M0Lpu*(xd33u + xd3u)) + 
     &                g1Ru*g2Ru*(2*g3Lu*M0Lpu*M3Lpu*xd33u + g3Ru*(xc21u + 2*xc221u + xc23u + 2*M0Lpu**2*xd33u - m1**2*xd33u - m2**2*xd33u + sqrtS**2*xd33u + M0Lpu**2*xd3u - m1**2*xd3u + 
     &                      M1Lpu**2*xd3u)))) + (-m2**2 + tman)*(g0Ru*(2*g1Ru*M1Lpu*(g2Lu*g3Lu*M2Lpu*xd2u + g2Ru*(g3Lu*M3Lpu*(2*xd23u + xd3u) + g3Ru*M0Lpu*(xd0u + 2*xd23u + xd2u + xd3u))) + 
     &                g1Lu*(2*g2Ru*g3Ru*M0Lpu*M2Lpu*xd2u + g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(2*xd23u + xd3u) + 
     &                      g3Lu*(xc03u + xc04u + 2*xc11u + 4*xc121u + 2*xc21u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + 4*M0Lpu**2*xd23u - 2*m1**2*xd23u - 2*m2**2*xd23u + 
     &                         2*sqrtS**2*xd23u + 2*M0Lpu**2*xd2u - m1**2*xd2u + m4**2*xd2u - uman*xd2u + 2*M0Lpu**2*xd3u - m1**2*xd3u - m2**2*xd3u + sqrtS**2*xd3u)))) + 
     &             g0Lu*(2*g1Lu*M1Lpu*(g2Ru*g3Ru*M2Lpu*xd2u + g2Lu*(g3Ru*M3Lpu*(2*xd23u + xd3u) + g3Lu*M0Lpu*(xd0u + 2*xd23u + xd2u + xd3u))) + 
     &                g1Ru*(2*g2Lu*g3Lu*M0Lpu*M2Lpu*xd2u + g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(2*xd23u + xd3u) + 
     &                      g3Ru*(xc03u + xc04u + 2*xc11u + 4*xc121u + 2*xc21u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + 4*M0Lpu**2*xd23u - 2*m1**2*xd23u - 2*m2**2*xd23u + 
     &                         2*sqrtS**2*xd23u + 2*M0Lpu**2*xd2u - m1**2*xd2u + m4**2*xd2u - uman*xd2u + 2*M0Lpu**2*xd3u - m1**2*xd3u - m2**2*xd3u + sqrtS**2*xd3u))))) - 
     &          ((3*m1**2 + 2*m2**2 + m4**2 - 2*sqrtS**2 - uman)*(g0Ru*(2*g1Ru*M1Lpu*
     &                   (g2Lu*g3Lu*M2Lpu*(xd1u + xd2u) + g2Ru*(g3Lu*M3Lpu*(2*xd13u + 2*xd23u + xd3u) + g3Ru*M0Lpu*(xd0u + 2*xd13u + xd1u + 2*xd23u + xd2u + xd3u))) + 
     &                  g1Lu*(2*g2Ru*M2Lpu*(g3Lu*M3Lpu*xd3u + g3Ru*M0Lpu*(xd0u + xd1u + xd2u + xd3u)) + 
     &                     g2Lu*(4*g3Ru*M0Lpu*M3Lpu*(xd13u + xd23u + xd3u) + 
     &                        g3Lu*(xc01u + xc03u - 2*xc21u - 4*xc221u + 2*M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u + 4*M0Lpu**2*xd13u - 2*m1**2*xd13u - 
     &                           2*m2**2*xd13u + 2*sqrtS**2*xd13u + 2*M0Lpu**2*xd1u - m1**2*xd1u + m4**2*xd1u - uman*xd1u + 4*M0Lpu**2*xd23u - 2*m1**2*xd23u - 2*m2**2*xd23u + 
     &                           2*sqrtS**2*xd23u + 2*M0Lpu**2*xd2u - m1**2*xd2u + m4**2*xd2u - uman*xd2u + 4*M0Lpu**2*xd3u - m1**2*xd3u - 2*m2**2*xd3u + m3**2*xd3u + sqrtS**2*xd3u - 
     &                           uman*xd3u)))) + g0Lu*(2*g1Lu*M1Lpu*(g2Ru*g3Ru*M2Lpu*(xd1u + xd2u) + 
     &                     g2Lu*(g3Ru*M3Lpu*(2*xd13u + 2*xd23u + xd3u) + g3Lu*M0Lpu*(xd0u + 2*xd13u + xd1u + 2*xd23u + xd2u + xd3u))) + 
     &                  g1Ru*(2*g2Lu*M2Lpu*(g3Ru*M3Lpu*xd3u + g3Lu*M0Lpu*(xd0u + xd1u + xd2u + xd3u)) + 
     &                     g2Ru*(4*g3Lu*M0Lpu*M3Lpu*(xd13u + xd23u + xd3u) + 
     &                        g3Ru*(xc01u + xc03u - 2*xc21u - 4*xc221u + 2*M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u + 4*M0Lpu**2*xd13u - 2*m1**2*xd13u - 
     &                           2*m2**2*xd13u + 2*sqrtS**2*xd13u + 2*M0Lpu**2*xd1u - m1**2*xd1u + m4**2*xd1u - uman*xd1u + 4*M0Lpu**2*xd23u - 2*m1**2*xd23u - 2*m2**2*xd23u + 
     &                           2*sqrtS**2*xd23u + 2*M0Lpu**2*xd2u - m1**2*xd2u + m4**2*xd2u - uman*xd2u + 4*M0Lpu**2*xd3u - m1**2*xd3u - 2*m2**2*xd3u + m3**2*xd3u + sqrtS**2*xd3u - 
     &                           uman*xd3u))))))/2.)))/(64.*pi**2*(-sfmc**2 + tman))  )


!Box(t)*u + Box(u)*u
 
      Coupl_Ac = dconjg(VSfSfu(ytype,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
      Coupl_Bc = dconjg(VSfSft(xtype,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
      sfmc = MSf(j,proptypeuc,propgenuc)
        

!t
      Box_TU = Box_TU + 2*colorfactorTU*dreal(  -(Coupl_Ac*Coupl_Bc*((6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*tman)*
     &        (g0Rt*(2*g1Rt*M1Lpt*(-(g2Lt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*xd0t + g3Lt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t))) + 
     &                g2Rt*(g3Rt*M0Lpt*(xc01t + xc02t - 4*xd00t + M2Lpt**2*xd0t + M3Lpt**2*xd0t - m4**2*xd0t) + 
     &                   g3Lt*M3Lpt*(xc01t + xc04t - 4*xd00t + M0Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t))) + 
     &             g1Lt*(-2*g2Rt*M2Lpt*(g3Lt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Rt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &                g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(xc01t + xc03t - 4*xd00t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - m3**2*xd0t) + 
     &                   g3Lt*(2*xb01t + 2*xb02t - 8*xc001t + 2*M0Lpt**2*xc01t - m1**2*xc01t - m2**2*xc01t + sqrtS**2*xc01t - m1**2*xc02t + 2*M1Lpt**2*xc02t - m3**2*xc02t + tman*xc02t - 
     &                      m2**2*xc03t + 2*M3Lpt**2*xc03t - m4**2*xc03t + tman*xc03t + 2*M1Lpt**2*xc04t + m3**2*xc04t - m4**2*xc04t + sqrtS**2*xc04t + 4*m3**2*xc11t + 2*m3**2*xc21t - 
     &                      2*m4**2*xc21t + 2*sqrtS**2*xc21t - 8*M0Lpt**2*xd00t + 4*m1**2*xd00t + 4*m2**2*xd00t - 4*sqrtS**2*xd00t - M1Lpt**2*m2**2*xd0t + 2*M0Lpt**2*M2Lpt**2*xd0t - 
     &                      m1**2*M2Lpt**2*xd0t - m2**2*M2Lpt**2*xd0t - M0Lpt**2*m3**2*xd0t + m2**2*m3**2*xd0t - m1**2*M3Lpt**2*xd0t + 2*M1Lpt**2*M3Lpt**2*xd0t - m3**2*M3Lpt**2*xd0t - 
     &                      M0Lpt**2*m4**2*xd0t + m1**2*m4**2*xd0t - M1Lpt**2*m4**2*xd0t + M0Lpt**2*sqrtS**2*xd0t + M2Lpt**2*sqrtS**2*xd0t + M1Lpt**2*tman*xd0t + M3Lpt**2*tman*xd0t - 
     &                      sqrtS**2*tman*xd0t)))) + g0Lt*(2*g1Lt*M1Lpt*(-(g2Rt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*xd0t + g3Rt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t))) + 
     &                g2Lt*(g3Lt*M0Lpt*(xc01t + xc02t - 4*xd00t + M2Lpt**2*xd0t + M3Lpt**2*xd0t - m4**2*xd0t) + 
     &                   g3Rt*M3Lpt*(xc01t + xc04t - 4*xd00t + M0Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t))) + 
     &             g1Rt*(-2*g2Lt*M2Lpt*(g3Rt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Lt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &                g2Rt*(2*g3Lt*M0Lpt*M3Lpt*(xc01t + xc03t - 4*xd00t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - m3**2*xd0t) + 
     &                   g3Rt*(2*xb01t + 2*xb02t - 8*xc001t + 2*M0Lpt**2*xc01t - m1**2*xc01t - m2**2*xc01t + sqrtS**2*xc01t - m1**2*xc02t + 2*M1Lpt**2*xc02t - m3**2*xc02t + tman*xc02t - 
     &                      m2**2*xc03t + 2*M3Lpt**2*xc03t - m4**2*xc03t + tman*xc03t + 2*M1Lpt**2*xc04t + m3**2*xc04t - m4**2*xc04t + sqrtS**2*xc04t + 4*m3**2*xc11t + 2*m3**2*xc21t - 
     &                      2*m4**2*xc21t + 2*sqrtS**2*xc21t - 8*M0Lpt**2*xd00t + 4*m1**2*xd00t + 4*m2**2*xd00t - 4*sqrtS**2*xd00t - M1Lpt**2*m2**2*xd0t + 2*M0Lpt**2*M2Lpt**2*xd0t - 
     &                      m1**2*M2Lpt**2*xd0t - m2**2*M2Lpt**2*xd0t - M0Lpt**2*m3**2*xd0t + m2**2*m3**2*xd0t - m1**2*M3Lpt**2*xd0t + 2*M1Lpt**2*M3Lpt**2*xd0t - m3**2*M3Lpt**2*xd0t - 
     &                      M0Lpt**2*m4**2*xd0t + m1**2*m4**2*xd0t - M1Lpt**2*m4**2*xd0t + M0Lpt**2*sqrtS**2*xd0t + M2Lpt**2*sqrtS**2*xd0t + M1Lpt**2*tman*xd0t + M3Lpt**2*tman*xd0t - 
     &                      sqrtS**2*tman*xd0t))))) - 2*(-m2**2 + uman)*(2*(2*m1**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*tman)*
     &           (g0Rt*(2*g1Rt*g2Rt*M1Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*(xd22t + xd2t) + 
     &                g1Lt*g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(xd22t + xd2t) + g3Lt*(2*xc111t + 2*xc11t + (2*M0Lpt**2 - m1**2 - m2**2 + sqrtS**2)*(xd22t + xd2t)))) + 
     &             g0Lt*(2*g1Lt*g2Lt*M1Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*(xd22t + xd2t) + 
     &                g1Rt*g2Rt*(2*g3Lt*M0Lpt*M3Lpt*(xd22t + xd2t) + g3Rt*(2*xc111t + 2*xc11t + (2*M0Lpt**2 - m1**2 - m2**2 + sqrtS**2)*(xd22t + xd2t))))) + 
     &          (3*m1**2 - m4**2 + uman)*(-(g0Rt*(2*g1Rt*g2Rt*M1Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*(2*xd12t + xd1t + 2*(xd22t + xd2t)) + 
     &                  g1Lt*(-2*g2Rt*M2Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*xd2t + 
     &                     g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(xd0t + 2*xd12t + xd1t + 2*xd22t + 3*xd2t) + 
     &                        g3Lt*(xc02t - xc04t - 2*xc11t - 4*xc121t - 2*xc21t + 2*xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 4*M0Lpt**2*xd12t - 2*m1**2*xd12t - 2*m2**2*xd12t + 
     &                           2*sqrtS**2*xd12t + 2*M0Lpt**2*xd1t - m1**2*xd1t - m2**2*xd1t + sqrtS**2*xd1t + 4*M0Lpt**2*xd22t - 2*m1**2*xd22t - 2*m2**2*xd22t + 2*sqrtS**2*xd22t + 
     &                           4*M0Lpt**2*xd2t - 2*m1**2*xd2t - 3*m2**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t + 2*sqrtS**2*xd2t + tman*xd2t))))) - 
     &             g0Lt*(2*g1Lt*g2Lt*M1Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*(2*xd12t + xd1t + 2*(xd22t + xd2t)) + 
     &                g1Rt*(-2*g2Lt*M2Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*xd2t + g2Rt*
     &                    (2*g3Lt*M0Lpt*M3Lpt*(xd0t + 2*xd12t + xd1t + 2*xd22t + 3*xd2t) + 
     &                      g3Rt*(xc02t - xc04t - 2*xc11t - 4*xc121t - 2*xc21t + 2*xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 4*M0Lpt**2*xd12t - 2*m1**2*xd12t - 2*m2**2*xd12t + 
     &                         2*sqrtS**2*xd12t + 2*M0Lpt**2*xd1t - m1**2*xd1t - m2**2*xd1t + sqrtS**2*xd1t + 4*M0Lpt**2*xd22t - 2*m1**2*xd22t - 2*m2**2*xd22t + 2*sqrtS**2*xd22t + 
     &                         4*M0Lpt**2*xd2t - 2*m1**2*xd2t - 3*m2**2*xd2t + 2*M3Lpt**2*xd2t - m4**2*xd2t + 2*sqrtS**2*xd2t + tman*xd2t))))) + 
     &          (2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(-(g0Rt*(2*g1Rt*M1Lpt*
     &                   (g2Lt*g3Lt*M2Lpt*xd2t + g2Rt*(g3Lt*M3Lpt*(2*xd23t + xd3t) + g3Rt*M0Lpt*(xd0t + 2*xd23t + xd2t + xd3t))) + 
     &                  g1Lt*(2*g2Rt*g3Rt*M0Lpt*M2Lpt*xd2t + g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(2*xd23t + xd3t) + 
     &                        g3Lt*(xc03t + xc04t + 2*xc11t + 4*xc121t + 2*xc21t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 4*M0Lpt**2*xd23t - 2*m1**2*xd23t - 2*m2**2*xd23t + 
     &                           2*sqrtS**2*xd23t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t + 2*M0Lpt**2*xd3t - m1**2*xd3t - m2**2*xd3t + sqrtS**2*xd3t))))) - 
     &             g0Lt*(2*g1Lt*M1Lpt*(g2Rt*g3Rt*M2Lpt*xd2t + g2Lt*(g3Rt*M3Lpt*(2*xd23t + xd3t) + g3Lt*M0Lpt*(xd0t + 2*xd23t + xd2t + xd3t))) + 
     &                g1Rt*(2*g2Lt*g3Lt*M0Lpt*M2Lpt*xd2t + g2Rt*(2*g3Lt*M0Lpt*M3Lpt*(2*xd23t + xd3t) + 
     &                      g3Rt*(xc03t + xc04t + 2*xc11t + 4*xc121t + 2*xc21t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 4*M0Lpt**2*xd23t - 2*m1**2*xd23t - 2*m2**2*xd23t + 
     &                         2*sqrtS**2*xd23t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t + 2*M0Lpt**2*xd3t - m1**2*xd3t - m2**2*xd3t + sqrtS**2*xd3t)))))) - 
     &       (-3*m2**2 + m3**2 - uman)*(2*(-2*m1**2 - 3*m2**2 - m4**2 + 2*sqrtS**2 + tman)*
     &           (g0Rt*(2*g1Rt*g2Rt*M1Lpt*(g3Lt*M3Lpt*xd33t + g3Rt*M0Lpt*(xd33t + xd3t)) + 
     &                g1Lt*g2Lt*(2*g3Rt*M0Lpt*M3Lpt*xd33t + g3Lt*(xc21t + 2*xc221t + xc23t + 2*M0Lpt**2*xd33t - m1**2*xd33t - m2**2*xd33t + sqrtS**2*xd33t + M0Lpt**2*xd3t - m1**2*xd3t + 
     &                      M1Lpt**2*xd3t))) + g0Lt*(2*g1Lt*g2Lt*M1Lpt*(g3Rt*M3Lpt*xd33t + g3Lt*M0Lpt*(xd33t + xd3t)) + 
     &                g1Rt*g2Rt*(2*g3Lt*M0Lpt*M3Lpt*xd33t + g3Rt*(xc21t + 2*xc221t + xc23t + 2*M0Lpt**2*xd33t - m1**2*xd33t - m2**2*xd33t + sqrtS**2*xd33t + M0Lpt**2*xd3t - m1**2*xd3t + 
     &                      M1Lpt**2*xd3t)))) + (2*m1**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*tman)*
     &           (g0Rt*(2*g1Rt*M1Lpt*(-(g2Lt*g3Lt*M2Lpt*xd2t) + g2Rt*(g3Lt*M3Lpt*(2*xd23t + xd3t) + g3Rt*M0Lpt*(xd0t + 2*xd23t + xd2t + xd3t))) + 
     &                g1Lt*(-2*g2Rt*g3Rt*M0Lpt*M2Lpt*xd2t + g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(2*xd23t + xd3t) + 
     &                      g3Lt*(xc03t + xc04t + 4*xc121t + 2*xc13t + 2*xc21t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 4*M0Lpt**2*xd23t - 2*m1**2*xd23t - 2*m2**2*xd23t + 
     &                         2*sqrtS**2*xd23t - m1**2*xd2t + 2*M1Lpt**2*xd2t - m3**2*xd2t + tman*xd2t + 2*M0Lpt**2*xd3t - m1**2*xd3t - m2**2*xd3t + sqrtS**2*xd3t)))) + 
     &             g0Lt*(2*g1Lt*M1Lpt*(-(g2Rt*g3Rt*M2Lpt*xd2t) + g2Lt*(g3Rt*M3Lpt*(2*xd23t + xd3t) + g3Lt*M0Lpt*(xd0t + 2*xd23t + xd2t + xd3t))) + 
     &                g1Rt*(-2*g2Lt*g3Lt*M0Lpt*M2Lpt*xd2t + g2Rt*(2*g3Lt*M0Lpt*M3Lpt*(2*xd23t + xd3t) + 
     &                      g3Rt*(xc03t + xc04t + 4*xc121t + 2*xc13t + 2*xc21t + M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + 4*M0Lpt**2*xd23t - 2*m1**2*xd23t - 2*m2**2*xd23t + 
     &                         2*sqrtS**2*xd23t - m1**2*xd2t + 2*M1Lpt**2*xd2t - m3**2*xd2t + tman*xd2t + 2*M0Lpt**2*xd3t - m1**2*xd3t - m2**2*xd3t + sqrtS**2*xd3t))))) + 
     &          (3*m1**2 - m4**2 + uman)*(g0Rt*(-2*g1Rt*M1Lpt*(-(g2Lt*g3Lt*M2Lpt*(xd1t + xd2t)) + 
     &                   g2Rt*(g3Lt*M3Lpt*(2*xd13t + 2*xd23t + xd3t) + g3Rt*M0Lpt*(xd0t + 2*xd13t + xd1t + 2*xd23t + xd2t + xd3t))) + 
     &                g1Lt*(2*g2Rt*M2Lpt*(g3Lt*M3Lpt*xd3t + g3Rt*M0Lpt*(xd0t + xd1t + xd2t + xd3t)) + 
     &                   g2Lt*(-4*g3Rt*M0Lpt*M3Lpt*(xd13t + xd23t + xd3t) + g3Lt*
     &                       (xc01t - xc03t - 2*xc13t + 2*xc21t + 4*xc221t + m1**2*xd0t - M1Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t - 4*M0Lpt**2*xd13t + 2*m1**2*xd13t + 2*m2**2*xd13t - 
     &                         2*sqrtS**2*xd13t + m1**2*xd1t - 2*M1Lpt**2*xd1t + m3**2*xd1t - tman*xd1t - 4*M0Lpt**2*xd23t + 2*m1**2*xd23t + 2*m2**2*xd23t - 2*sqrtS**2*xd23t + m1**2*xd2t - 
     &                         2*M1Lpt**2*xd2t + m3**2*xd2t - tman*xd2t - 2*M0Lpt**2*xd3t + m1**2*xd3t + 2*m2**2*xd3t - 2*M3Lpt**2*xd3t + m4**2*xd3t - sqrtS**2*xd3t - tman*xd3t)))) + 
     &             g0Lt*(-2*g1Lt*M1Lpt*(-(g2Rt*g3Rt*M2Lpt*(xd1t + xd2t)) + g2Lt*(g3Rt*M3Lpt*(2*xd13t + 2*xd23t + xd3t) + g3Lt*M0Lpt*(xd0t + 2*xd13t + xd1t + 2*xd23t + xd2t + xd3t))) + 
     &                g1Rt*(2*g2Lt*M2Lpt*(g3Rt*M3Lpt*xd3t + g3Lt*M0Lpt*(xd0t + xd1t + xd2t + xd3t)) + 
     &                   g2Rt*(-4*g3Lt*M0Lpt*M3Lpt*(xd13t + xd23t + xd3t) + g3Rt*
     &                       (xc01t - xc03t - 2*xc13t + 2*xc21t + 4*xc221t + m1**2*xd0t - M1Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t - 4*M0Lpt**2*xd13t + 2*m1**2*xd13t + 2*m2**2*xd13t - 
     &                         2*sqrtS**2*xd13t + m1**2*xd1t - 2*M1Lpt**2*xd1t + m3**2*xd1t - tman*xd1t - 4*M0Lpt**2*xd23t + 2*m1**2*xd23t + 2*m2**2*xd23t - 2*sqrtS**2*xd23t + m1**2*xd2t - 
     &                         2*M1Lpt**2*xd2t + m3**2*xd2t - tman*xd2t - 2*M0Lpt**2*xd3t + m1**2*xd3t + 2*m2**2*xd3t - 2*M3Lpt**2*xd3t + m4**2*xd3t - sqrtS**2*xd3t - tman*xd3t)))))) + 
     &       (3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*(-2*(3*m1**2 - m4**2 + uman)*
     &           (g0Rt*(2*g1Rt*g2Rt*M1Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*(xd11t + 2*xd12t + xd1t + xd22t + xd2t) + 
     &                g1Lt*g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(xd0t + xd11t + 2*xd12t + 2*xd1t + xd22t + 2*xd2t) + 
     &                   g3Lt*(xc02t + xc12t + xc21t + 2*xc221t + xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 2*M0Lpt**2*xd11t - m1**2*xd11t - m2**2*xd11t + sqrtS**2*xd11t + 
     &                      4*M0Lpt**2*xd12t - 2*m1**2*xd12t - 2*m2**2*xd12t + 2*sqrtS**2*xd12t + 3*M0Lpt**2*xd1t - m1**2*xd1t - 2*m2**2*xd1t + M3Lpt**2*xd1t + sqrtS**2*xd1t + 
     &                      2*M0Lpt**2*xd22t - m1**2*xd22t - m2**2*xd22t + sqrtS**2*xd22t + 3*M0Lpt**2*xd2t - m1**2*xd2t - 2*m2**2*xd2t + M3Lpt**2*xd2t + sqrtS**2*xd2t))) + 
     &             g0Lt*(2*g1Lt*g2Lt*M1Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*(xd11t + 2*xd12t + xd1t + xd22t + xd2t) + 
     &                g1Rt*g2Rt*(2*g3Lt*M0Lpt*M3Lpt*(xd0t + xd11t + 2*xd12t + 2*xd1t + xd22t + 2*xd2t) + 
     &                   g3Rt*(xc02t + xc12t + xc21t + 2*xc221t + xc22t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 2*M0Lpt**2*xd11t - m1**2*xd11t - m2**2*xd11t + sqrtS**2*xd11t + 
     &                      4*M0Lpt**2*xd12t - 2*m1**2*xd12t - 2*m2**2*xd12t + 2*sqrtS**2*xd12t + 3*M0Lpt**2*xd1t - m1**2*xd1t - 2*m2**2*xd1t + M3Lpt**2*xd1t + sqrtS**2*xd1t + 
     &                      2*M0Lpt**2*xd22t - m1**2*xd22t - m2**2*xd22t + sqrtS**2*xd22t + 3*M0Lpt**2*xd2t - m1**2*xd2t - 2*m2**2*xd2t + M3Lpt**2*xd2t + sqrtS**2*xd2t)))) + 
     &          (2*m1**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*tman)*(g0Rt*(2*g1Rt*g2Rt*M1Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*(2*xd12t + xd1t + 2*(xd22t + xd2t)) + 
     &                g1Lt*(2*g2Rt*M2Lpt*(g3Rt*M0Lpt + g3Lt*M3Lpt)*xd2t + g2Lt*
     &                    (2*g3Rt*M0Lpt*M3Lpt*(xd0t + 2*xd12t + xd1t + 2*xd22t + 3*xd2t) + 
     &                      g3Lt*(xc02t - xc04t - 4*xc121t - 2*xc21t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 4*M0Lpt**2*xd12t - 2*m1**2*xd12t - 2*m2**2*xd12t + 2*sqrtS**2*xd12t + 
     &                         2*M0Lpt**2*xd1t - m1**2*xd1t - m2**2*xd1t + sqrtS**2*xd1t + 4*M0Lpt**2*xd22t - 2*m1**2*xd22t - 2*m2**2*xd22t + 2*sqrtS**2*xd22t + 6*M0Lpt**2*xd2t - 
     &                         2*m1**2*xd2t - 3*m2**2*xd2t + m4**2*xd2t + 2*sqrtS**2*xd2t - tman*xd2t)))) + 
     &             g0Lt*(2*g1Lt*g2Lt*M1Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*(2*xd12t + xd1t + 2*(xd22t + xd2t)) + 
     &                g1Rt*(2*g2Lt*M2Lpt*(g3Lt*M0Lpt + g3Rt*M3Lpt)*xd2t + g2Rt*
     &                    (2*g3Lt*M0Lpt*M3Lpt*(xd0t + 2*xd12t + xd1t + 2*xd22t + 3*xd2t) + 
     &                      g3Rt*(xc02t - xc04t - 4*xc121t - 2*xc21t + M0Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 4*M0Lpt**2*xd12t - 2*m1**2*xd12t - 2*m2**2*xd12t + 2*sqrtS**2*xd12t + 
     &                         2*M0Lpt**2*xd1t - m1**2*xd1t - m2**2*xd1t + sqrtS**2*xd1t + 4*M0Lpt**2*xd22t - 2*m1**2*xd22t - 2*m2**2*xd22t + 2*sqrtS**2*xd22t + 6*M0Lpt**2*xd2t - 
     &                         2*m1**2*xd2t - 3*m2**2*xd2t + m4**2*xd2t + 2*sqrtS**2*xd2t - tman*xd2t))))) + 
     &          (2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(-(g0Rt*(2*g1Rt*M1Lpt*
     &                   (g2Lt*g3Lt*M2Lpt*(xd1t + xd2t) + g2Rt*(g3Lt*M3Lpt*(2*xd13t + 2*xd23t + xd3t) + g3Rt*M0Lpt*(xd0t + 2*xd13t + xd1t + 2*xd23t + xd2t + xd3t))) + 
     &                  g1Lt*(2*g2Rt*M2Lpt*(g3Lt*M3Lpt*xd3t + g3Rt*M0Lpt*(xd0t + xd1t + xd2t + xd3t)) + 
     &                     g2Lt*(4*g3Rt*M0Lpt*M3Lpt*(xd13t + xd23t + xd3t) + 
     &                        g3Lt*(xc01t + xc03t - 2*xc21t - 4*xc221t + 2*M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t + 4*M0Lpt**2*xd13t - 2*m1**2*xd13t - 
     &                           2*m2**2*xd13t + 2*sqrtS**2*xd13t + 2*M0Lpt**2*xd1t - m1**2*xd1t + m3**2*xd1t - tman*xd1t + 4*M0Lpt**2*xd23t - 2*m1**2*xd23t - 2*m2**2*xd23t + 
     &                           2*sqrtS**2*xd23t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t + 4*M0Lpt**2*xd3t - m1**2*xd3t - 2*m2**2*xd3t + m4**2*xd3t + sqrtS**2*xd3t - 
     &                           tman*xd3t))))) - g0Lt*(2*g1Lt*M1Lpt*(g2Rt*g3Rt*M2Lpt*(xd1t + xd2t) + 
     &                   g2Lt*(g3Rt*M3Lpt*(2*xd13t + 2*xd23t + xd3t) + g3Lt*M0Lpt*(xd0t + 2*xd13t + xd1t + 2*xd23t + xd2t + xd3t))) + 
     &                g1Rt*(2*g2Lt*M2Lpt*(g3Rt*M3Lpt*xd3t + g3Lt*M0Lpt*(xd0t + xd1t + xd2t + xd3t)) + 
     &                   g2Rt*(4*g3Lt*M0Lpt*M3Lpt*(xd13t + xd23t + xd3t) + g3Rt*
     &                       (xc01t + xc03t - 2*xc21t - 4*xc221t + 2*M0Lpt**2*xd0t - m1**2*xd0t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t + 4*M0Lpt**2*xd13t - 2*m1**2*xd13t - 
     &                         2*m2**2*xd13t + 2*sqrtS**2*xd13t + 2*M0Lpt**2*xd1t - m1**2*xd1t + m3**2*xd1t - tman*xd1t + 4*M0Lpt**2*xd23t - 2*m1**2*xd23t - 2*m2**2*xd23t + 
     &                         2*sqrtS**2*xd23t + 2*M0Lpt**2*xd2t - m1**2*xd2t + m3**2*xd2t - tman*xd2t + 4*M0Lpt**2*xd3t - m1**2*xd3t - 2*m2**2*xd3t + m4**2*xd3t + sqrtS**2*xd3t - tman*xd3t)
     &                      )))))))/(64.*pi**2*(-sfmc**2 + uman))  )

!u
      Box_UU = Box_UU + 2*colorfactorTU*dreal(  -(Coupl_Ac*Coupl_Bc*((6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*tman)*
     &        (g0Ru*(2*g1Ru*M1Lpu*(-(g2Lu*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*xd0u + g3Lu*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u))) + 
     &                g2Ru*(g3Ru*M0Lpu*(xc01u + xc02u - 4*xd00u + M2Lpu**2*xd0u - m3**2*xd0u + M3Lpu**2*xd0u) + 
     &                   g3Lu*M3Lpu*(xc01u + xc04u - 4*xd00u + M0Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u))) + 
     &             g1Lu*(-2*g2Ru*M2Lpu*(g3Lu*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Ru*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &                g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(xc01u + xc03u - 4*xd00u + M1Lpu**2*xd0u + M2Lpu**2*xd0u - m4**2*xd0u) + 
     &                   g3Lu*(2*xb01u + 2*xb02u - 8*xc001u + 2*M0Lpu**2*xc01u - m1**2*xc01u - m2**2*xc01u + sqrtS**2*xc01u - m1**2*xc02u + 2*M1Lpu**2*xc02u - m4**2*xc02u + uman*xc02u - 
     &                      m2**2*xc03u - m3**2*xc03u + 2*M3Lpu**2*xc03u + uman*xc03u + 2*M1Lpu**2*xc04u - m3**2*xc04u + m4**2*xc04u + sqrtS**2*xc04u + 4*m4**2*xc11u - 2*m3**2*xc21u + 
     &                      2*m4**2*xc21u + 2*sqrtS**2*xc21u - 8*M0Lpu**2*xd00u + 4*m1**2*xd00u + 4*m2**2*xd00u - 4*sqrtS**2*xd00u - M1Lpu**2*m2**2*xd0u + 2*M0Lpu**2*M2Lpu**2*xd0u - 
     &                      m1**2*M2Lpu**2*xd0u - m2**2*M2Lpu**2*xd0u - M0Lpu**2*m3**2*xd0u + m1**2*m3**2*xd0u - M1Lpu**2*m3**2*xd0u - m1**2*M3Lpu**2*xd0u + 2*M1Lpu**2*M3Lpu**2*xd0u - 
     &                      M0Lpu**2*m4**2*xd0u + m2**2*m4**2*xd0u - M3Lpu**2*m4**2*xd0u + M0Lpu**2*sqrtS**2*xd0u + M2Lpu**2*sqrtS**2*xd0u + M1Lpu**2*uman*xd0u + M3Lpu**2*uman*xd0u - 
     &                      sqrtS**2*uman*xd0u)))) + g0Lu*(2*g1Lu*M1Lpu*(-(g2Ru*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*xd0u + g3Ru*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u))) + 
     &                g2Lu*(g3Lu*M0Lpu*(xc01u + xc02u - 4*xd00u + M2Lpu**2*xd0u - m3**2*xd0u + M3Lpu**2*xd0u) + 
     &                   g3Ru*M3Lpu*(xc01u + xc04u - 4*xd00u + M0Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u))) + 
     &             g1Ru*(-2*g2Lu*M2Lpu*(g3Ru*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Lu*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &                g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(xc01u + xc03u - 4*xd00u + M1Lpu**2*xd0u + M2Lpu**2*xd0u - m4**2*xd0u) + 
     &                   g3Ru*(2*xb01u + 2*xb02u - 8*xc001u + 2*M0Lpu**2*xc01u - m1**2*xc01u - m2**2*xc01u + sqrtS**2*xc01u - m1**2*xc02u + 2*M1Lpu**2*xc02u - m4**2*xc02u + uman*xc02u - 
     &                      m2**2*xc03u - m3**2*xc03u + 2*M3Lpu**2*xc03u + uman*xc03u + 2*M1Lpu**2*xc04u - m3**2*xc04u + m4**2*xc04u + sqrtS**2*xc04u + 4*m4**2*xc11u - 2*m3**2*xc21u + 
     &                      2*m4**2*xc21u + 2*sqrtS**2*xc21u - 8*M0Lpu**2*xd00u + 4*m1**2*xd00u + 4*m2**2*xd00u - 4*sqrtS**2*xd00u - M1Lpu**2*m2**2*xd0u + 2*M0Lpu**2*M2Lpu**2*xd0u - 
     &                      m1**2*M2Lpu**2*xd0u - m2**2*M2Lpu**2*xd0u - M0Lpu**2*m3**2*xd0u + m1**2*m3**2*xd0u - M1Lpu**2*m3**2*xd0u - m1**2*M3Lpu**2*xd0u + 2*M1Lpu**2*M3Lpu**2*xd0u - 
     &                      M0Lpu**2*m4**2*xd0u + m2**2*m4**2*xd0u - M3Lpu**2*m4**2*xd0u + M0Lpu**2*sqrtS**2*xd0u + M2Lpu**2*sqrtS**2*xd0u + M1Lpu**2*uman*xd0u + M3Lpu**2*uman*xd0u - 
     &                      sqrtS**2*uman*xd0u))))) - (-2*m2**2 - m3**2 - 3*m4**2 + sqrtS**2 + 2*tman)*
     &        (4*(m1**2 - uman)*(g0Ru*(2*g1Ru*g2Ru*M1Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*(xd22u + xd2u) + 
     &                g1Lu*g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(xd22u + xd2u) + g3Lu*(2*xc111u + 2*xc11u + (2*M0Lpu**2 - m1**2 - m2**2 + sqrtS**2)*(xd22u + xd2u)))) + 
     &             g0Lu*(2*g1Lu*g2Lu*M1Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*(xd22u + xd2u) + 
     &                g1Ru*g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(xd22u + xd2u) + g3Ru*(2*xc111u + 2*xc11u + (2*M0Lpu**2 - m1**2 - m2**2 + sqrtS**2)*(xd22u + xd2u))))) + 
     &          (3*m1**2 - m4**2 + uman)*(-(g0Ru*(2*g1Ru*g2Ru*M1Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*(2*xd12u + xd1u + 2*(xd22u + xd2u)) + 
     &                  g1Lu*(2*g2Ru*M2Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*xd2u + g2Lu*
     &                      (2*g3Ru*M0Lpu*M3Lpu*(xd0u + 2*xd12u + xd1u + 2*xd22u + 3*xd2u) + 
     &                        g3Lu*(xc02u - xc04u - 4*xc121u - 2*xc21u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 4*M0Lpu**2*xd12u - 2*m1**2*xd12u - 2*m2**2*xd12u + 2*sqrtS**2*xd12u + 
     &                           2*M0Lpu**2*xd1u - m1**2*xd1u - m2**2*xd1u + sqrtS**2*xd1u + 4*M0Lpu**2*xd22u - 2*m1**2*xd22u - 2*m2**2*xd22u + 2*sqrtS**2*xd22u + 6*M0Lpu**2*xd2u - 
     &                           2*m1**2*xd2u - 3*m2**2*xd2u + m3**2*xd2u + 2*sqrtS**2*xd2u - uman*xd2u))))) - 
     &             g0Lu*(2*g1Lu*g2Lu*M1Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*(2*xd12u + xd1u + 2*(xd22u + xd2u)) + 
     &                g1Ru*(2*g2Lu*M2Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*xd2u + g2Ru*
     &                    (2*g3Lu*M0Lpu*M3Lpu*(xd0u + 2*xd12u + xd1u + 2*xd22u + 3*xd2u) + 
     &                      g3Ru*(xc02u - xc04u - 4*xc121u - 2*xc21u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 4*M0Lpu**2*xd12u - 2*m1**2*xd12u - 2*m2**2*xd12u + 2*sqrtS**2*xd12u + 
     &                         2*M0Lpu**2*xd1u - m1**2*xd1u - m2**2*xd1u + sqrtS**2*xd1u + 4*M0Lpu**2*xd22u - 2*m1**2*xd22u - 2*m2**2*xd22u + 2*sqrtS**2*xd22u + 6*M0Lpu**2*xd2u - 
     &                         2*m1**2*xd2u - 3*m2**2*xd2u + m3**2*xd2u + 2*sqrtS**2*xd2u - uman*xd2u))))) + 
     &          (2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(-(g0Ru*(2*g1Ru*M1Lpu*
     &                   (-(g2Lu*g3Lu*M2Lpu*xd2u) + g2Ru*(g3Lu*M3Lpu*(2*xd23u + xd3u) + g3Ru*M0Lpu*(xd0u + 2*xd23u + xd2u + xd3u))) + 
     &                  g1Lu*(-2*g2Ru*g3Ru*M0Lpu*M2Lpu*xd2u + g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(2*xd23u + xd3u) + 
     &                        g3Lu*(xc03u + xc04u + 4*xc121u + 2*xc13u + 2*xc21u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + 4*M0Lpu**2*xd23u - 2*m1**2*xd23u - 2*m2**2*xd23u + 
     &                           2*sqrtS**2*xd23u - m1**2*xd2u + 2*M1Lpu**2*xd2u - m4**2*xd2u + uman*xd2u + 2*M0Lpu**2*xd3u - m1**2*xd3u - m2**2*xd3u + sqrtS**2*xd3u))))) - 
     &             g0Lu*(2*g1Lu*M1Lpu*(-(g2Ru*g3Ru*M2Lpu*xd2u) + g2Lu*(g3Ru*M3Lpu*(2*xd23u + xd3u) + g3Lu*M0Lpu*(xd0u + 2*xd23u + xd2u + xd3u))) + 
     &                g1Ru*(-2*g2Lu*g3Lu*M0Lpu*M2Lpu*xd2u + g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(2*xd23u + xd3u) + 
     &                      g3Ru*(xc03u + xc04u + 4*xc121u + 2*xc13u + 2*xc21u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + 4*M0Lpu**2*xd23u - 2*m1**2*xd23u - 2*m2**2*xd23u + 
     &                         2*sqrtS**2*xd23u - m1**2*xd2u + 2*M1Lpu**2*xd2u - m4**2*xd2u + uman*xd2u + 2*M0Lpu**2*xd3u - m1**2*xd3u - m2**2*xd3u + sqrtS**2*xd3u)))))) - 
     &       2*(3*m1**2 + 2*m2**2 + m3**2 - 2*sqrtS**2 - tman)*((3*m1**2 - m4**2 + uman)*
     &           (g0Ru*(2*g1Ru*g2Ru*M1Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*(xd11u + 2*xd12u + xd1u + xd22u + xd2u) + 
     &                g1Lu*g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(xd0u + xd11u + 2*xd12u + 2*xd1u + xd22u + 2*xd2u) + 
     &                   g3Lu*(xc02u + xc12u + xc21u + 2*xc221u + xc22u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd11u - m1**2*xd11u - m2**2*xd11u + sqrtS**2*xd11u + 
     &                      4*M0Lpu**2*xd12u - 2*m1**2*xd12u - 2*m2**2*xd12u + 2*sqrtS**2*xd12u + 3*M0Lpu**2*xd1u - m1**2*xd1u - 2*m2**2*xd1u + M3Lpu**2*xd1u + sqrtS**2*xd1u + 
     &                      2*M0Lpu**2*xd22u - m1**2*xd22u - m2**2*xd22u + sqrtS**2*xd22u + 3*M0Lpu**2*xd2u - m1**2*xd2u - 2*m2**2*xd2u + M3Lpu**2*xd2u + sqrtS**2*xd2u))) + 
     &             g0Lu*(2*g1Lu*g2Lu*M1Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*(xd11u + 2*xd12u + xd1u + xd22u + xd2u) + 
     &                g1Ru*g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(xd0u + xd11u + 2*xd12u + 2*xd1u + xd22u + 2*xd2u) + 
     &                   g3Ru*(xc02u + xc12u + xc21u + 2*xc221u + xc22u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 2*M0Lpu**2*xd11u - m1**2*xd11u - m2**2*xd11u + sqrtS**2*xd11u + 
     &                      4*M0Lpu**2*xd12u - 2*m1**2*xd12u - 2*m2**2*xd12u + 2*sqrtS**2*xd12u + 3*M0Lpu**2*xd1u - m1**2*xd1u - 2*m2**2*xd1u + M3Lpu**2*xd1u + sqrtS**2*xd1u + 
     &                      2*M0Lpu**2*xd22u - m1**2*xd22u - m2**2*xd22u + sqrtS**2*xd22u + 3*M0Lpu**2*xd2u - m1**2*xd2u - 2*m2**2*xd2u + M3Lpu**2*xd2u + sqrtS**2*xd2u)))) - 
     &          (m1**2 - uman)*(g0Ru*(2*g1Ru*g2Ru*M1Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*(2*xd12u + xd1u + 2*(xd22u + xd2u)) + 
     &                g1Lu*(-2*g2Ru*M2Lpu*(g3Ru*M0Lpu + g3Lu*M3Lpu)*xd2u + g2Lu*
     &                    (2*g3Ru*M0Lpu*M3Lpu*(xd0u + 2*xd12u + xd1u + 2*xd22u + 3*xd2u) + 
     &                      g3Lu*(xc02u - xc04u - 2*xc11u - 4*xc121u - 2*xc21u + 2*xc22u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 4*M0Lpu**2*xd12u - 2*m1**2*xd12u - 2*m2**2*xd12u + 
     &                         2*sqrtS**2*xd12u + 2*M0Lpu**2*xd1u - m1**2*xd1u - m2**2*xd1u + sqrtS**2*xd1u + 4*M0Lpu**2*xd22u - 2*m1**2*xd22u - 2*m2**2*xd22u + 2*sqrtS**2*xd22u + 
     &                         4*M0Lpu**2*xd2u - 2*m1**2*xd2u - 3*m2**2*xd2u - m3**2*xd2u + 2*M3Lpu**2*xd2u + 2*sqrtS**2*xd2u + uman*xd2u)))) + 
     &             g0Lu*(2*g1Lu*g2Lu*M1Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*(2*xd12u + xd1u + 2*(xd22u + xd2u)) + 
     &                g1Ru*(-2*g2Lu*M2Lpu*(g3Lu*M0Lpu + g3Ru*M3Lpu)*xd2u + g2Ru*
     &                    (2*g3Lu*M0Lpu*M3Lpu*(xd0u + 2*xd12u + xd1u + 2*xd22u + 3*xd2u) + 
     &                      g3Ru*(xc02u - xc04u - 2*xc11u - 4*xc121u - 2*xc21u + 2*xc22u + M0Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 4*M0Lpu**2*xd12u - 2*m1**2*xd12u - 2*m2**2*xd12u + 
     &                         2*sqrtS**2*xd12u + 2*M0Lpu**2*xd1u - m1**2*xd1u - m2**2*xd1u + sqrtS**2*xd1u + 4*M0Lpu**2*xd22u - 2*m1**2*xd22u - 2*m2**2*xd22u + 2*sqrtS**2*xd22u + 
     &                         4*M0Lpu**2*xd2u - 2*m1**2*xd2u - 3*m2**2*xd2u - m3**2*xd2u + 2*M3Lpu**2*xd2u + 2*sqrtS**2*xd2u + uman*xd2u))))) - 
     &          ((2*m1**2 + 3*m2**2 + m4**2 - 2*sqrtS**2 - tman)*(g0Ru*(-2*g1Ru*M1Lpu*
     &                   (-(g2Lu*g3Lu*M2Lpu*(xd1u + xd2u)) + g2Ru*(g3Lu*M3Lpu*(2*xd13u + 2*xd23u + xd3u) + g3Ru*M0Lpu*(xd0u + 2*xd13u + xd1u + 2*xd23u + xd2u + xd3u))) + 
     &                  g1Lu*(2*g2Ru*M2Lpu*(g3Lu*M3Lpu*xd3u + g3Ru*M0Lpu*(xd0u + xd1u + xd2u + xd3u)) + 
     &                     g2Lu*(-4*g3Ru*M0Lpu*M3Lpu*(xd13u + xd23u + xd3u) + 
     &                        g3Lu*(xc01u - xc03u - 2*xc13u + 2*xc21u + 4*xc221u + m1**2*xd0u - M1Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u - 4*M0Lpu**2*xd13u + 2*m1**2*xd13u + 2*m2**2*xd13u - 
     &                           2*sqrtS**2*xd13u + m1**2*xd1u - 2*M1Lpu**2*xd1u + m4**2*xd1u - uman*xd1u - 4*M0Lpu**2*xd23u + 2*m1**2*xd23u + 2*m2**2*xd23u - 2*sqrtS**2*xd23u + m1**2*xd2u - 
     &                           2*M1Lpu**2*xd2u + m4**2*xd2u - uman*xd2u - 2*M0Lpu**2*xd3u + m1**2*xd3u + 2*m2**2*xd3u + m3**2*xd3u - 2*M3Lpu**2*xd3u - sqrtS**2*xd3u - uman*xd3u)))) + 
     &               g0Lu*(-2*g1Lu*M1Lpu*(-(g2Ru*g3Ru*M2Lpu*(xd1u + xd2u)) + g2Lu*(g3Ru*M3Lpu*(2*xd13u + 2*xd23u + xd3u) + g3Lu*M0Lpu*(xd0u + 2*xd13u + xd1u + 2*xd23u + xd2u + xd3u))) + 
     &                  g1Ru*(2*g2Lu*M2Lpu*(g3Ru*M3Lpu*xd3u + g3Lu*M0Lpu*(xd0u + xd1u + xd2u + xd3u)) + 
     &                     g2Ru*(-4*g3Lu*M0Lpu*M3Lpu*(xd13u + xd23u + xd3u) + 
     &                        g3Ru*(xc01u - xc03u - 2*xc13u + 2*xc21u + 4*xc221u + m1**2*xd0u - M1Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u - 4*M0Lpu**2*xd13u + 2*m1**2*xd13u + 2*m2**2*xd13u - 
     &                           2*sqrtS**2*xd13u + m1**2*xd1u - 2*M1Lpu**2*xd1u + m4**2*xd1u - uman*xd1u - 4*M0Lpu**2*xd23u + 2*m1**2*xd23u + 2*m2**2*xd23u - 2*sqrtS**2*xd23u + m1**2*xd2u - 
     &                           2*M1Lpu**2*xd2u + m4**2*xd2u - uman*xd2u - 2*M0Lpu**2*xd3u + m1**2*xd3u + 2*m2**2*xd3u + m3**2*xd3u - 2*M3Lpu**2*xd3u - sqrtS**2*xd3u - uman*xd3u))))))/2.) - 
     &       2*(-3*m2**2 + m3**2 - uman)*((-2*m1**2 - 3*m2**2 - m4**2 + 2*sqrtS**2 + tman)*
     &           (g0Ru*(2*g1Ru*g2Ru*M1Lpu*(g3Lu*M3Lpu*xd33u + g3Ru*M0Lpu*(xd33u + xd3u)) + 
     &                g1Lu*g2Lu*(2*g3Ru*M0Lpu*M3Lpu*xd33u + g3Lu*(xc21u + 2*xc221u + xc23u + 2*M0Lpu**2*xd33u - m1**2*xd33u - m2**2*xd33u + sqrtS**2*xd33u + M0Lpu**2*xd3u - m1**2*xd3u + 
     &                      M1Lpu**2*xd3u))) + g0Lu*(2*g1Lu*g2Lu*M1Lpu*(g3Ru*M3Lpu*xd33u + g3Lu*M0Lpu*(xd33u + xd3u)) + 
     &                g1Ru*g2Ru*(2*g3Lu*M0Lpu*M3Lpu*xd33u + g3Ru*(xc21u + 2*xc221u + xc23u + 2*M0Lpu**2*xd33u - m1**2*xd33u - m2**2*xd33u + sqrtS**2*xd33u + M0Lpu**2*xd3u - m1**2*xd3u + 
     &                      M1Lpu**2*xd3u)))) + (m1**2 - uman)*(g0Ru*(2*g1Ru*M1Lpu*(g2Lu*g3Lu*M2Lpu*xd2u + g2Ru*(g3Lu*M3Lpu*(2*xd23u + xd3u) + g3Ru*M0Lpu*(xd0u + 2*xd23u + xd2u + xd3u))) + 
     &                g1Lu*(2*g2Ru*g3Ru*M0Lpu*M2Lpu*xd2u + g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(2*xd23u + xd3u) + 
     &                      g3Lu*(xc03u + xc04u + 2*xc11u + 4*xc121u + 2*xc21u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + 4*M0Lpu**2*xd23u - 2*m1**2*xd23u - 2*m2**2*xd23u + 
     &                         2*sqrtS**2*xd23u + 2*M0Lpu**2*xd2u - m1**2*xd2u + m4**2*xd2u - uman*xd2u + 2*M0Lpu**2*xd3u - m1**2*xd3u - m2**2*xd3u + sqrtS**2*xd3u)))) + 
     &             g0Lu*(2*g1Lu*M1Lpu*(g2Ru*g3Ru*M2Lpu*xd2u + g2Lu*(g3Ru*M3Lpu*(2*xd23u + xd3u) + g3Lu*M0Lpu*(xd0u + 2*xd23u + xd2u + xd3u))) + 
     &                g1Ru*(2*g2Lu*g3Lu*M0Lpu*M2Lpu*xd2u + g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(2*xd23u + xd3u) + 
     &                      g3Ru*(xc03u + xc04u + 2*xc11u + 4*xc121u + 2*xc21u + M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + 4*M0Lpu**2*xd23u - 2*m1**2*xd23u - 2*m2**2*xd23u + 
     &                         2*sqrtS**2*xd23u + 2*M0Lpu**2*xd2u - m1**2*xd2u + m4**2*xd2u - uman*xd2u + 2*M0Lpu**2*xd3u - m1**2*xd3u - m2**2*xd3u + sqrtS**2*xd3u))))) + 
     &          ((3*m1**2 - m4**2 + uman)*(-(g0Ru*(2*g1Ru*M1Lpu*(g2Lu*g3Lu*M2Lpu*(xd1u + xd2u) + 
     &                       g2Ru*(g3Lu*M3Lpu*(2*xd13u + 2*xd23u + xd3u) + g3Ru*M0Lpu*(xd0u + 2*xd13u + xd1u + 2*xd23u + xd2u + xd3u))) + 
     &                    g1Lu*(2*g2Ru*M2Lpu*(g3Lu*M3Lpu*xd3u + g3Ru*M0Lpu*(xd0u + xd1u + xd2u + xd3u)) + 
     &                       g2Lu*(4*g3Ru*M0Lpu*M3Lpu*(xd13u + xd23u + xd3u) + 
     &                          g3Lu*(xc01u + xc03u - 2*xc21u - 4*xc221u + 2*M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u + 4*M0Lpu**2*xd13u - 2*m1**2*xd13u - 
     &                             2*m2**2*xd13u + 2*sqrtS**2*xd13u + 2*M0Lpu**2*xd1u - m1**2*xd1u + m4**2*xd1u - uman*xd1u + 4*M0Lpu**2*xd23u - 2*m1**2*xd23u - 2*m2**2*xd23u + 
     &                             2*sqrtS**2*xd23u + 2*M0Lpu**2*xd2u - m1**2*xd2u + m4**2*xd2u - uman*xd2u + 4*M0Lpu**2*xd3u - m1**2*xd3u - 2*m2**2*xd3u + m3**2*xd3u + sqrtS**2*xd3u - 
     &                             uman*xd3u))))) - g0Lu*(2*g1Lu*M1Lpu*(g2Ru*g3Ru*M2Lpu*(xd1u + xd2u) + 
     &                     g2Lu*(g3Ru*M3Lpu*(2*xd13u + 2*xd23u + xd3u) + g3Lu*M0Lpu*(xd0u + 2*xd13u + xd1u + 2*xd23u + xd2u + xd3u))) + 
     &                  g1Ru*(2*g2Lu*M2Lpu*(g3Ru*M3Lpu*xd3u + g3Lu*M0Lpu*(xd0u + xd1u + xd2u + xd3u)) + 
     &                     g2Ru*(4*g3Lu*M0Lpu*M3Lpu*(xd13u + xd23u + xd3u) + 
     &                        g3Ru*(xc01u + xc03u - 2*xc21u - 4*xc221u + 2*M0Lpu**2*xd0u - m1**2*xd0u + M1Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u + 4*M0Lpu**2*xd13u - 2*m1**2*xd13u - 
     &                           2*m2**2*xd13u + 2*sqrtS**2*xd13u + 2*M0Lpu**2*xd1u - m1**2*xd1u + m4**2*xd1u - uman*xd1u + 4*M0Lpu**2*xd23u - 2*m1**2*xd23u - 2*m2**2*xd23u + 
     &                           2*sqrtS**2*xd23u + 2*M0Lpu**2*xd2u - m1**2*xd2u + m4**2*xd2u - uman*xd2u + 4*M0Lpu**2*xd3u - m1**2*xd3u - 2*m2**2*xd3u + m3**2*xd3u + sqrtS**2*xd3u - 
     &                           uman*xd3u))))))/2.)))/(64.*pi**2*(-sfmc**2 + uman))  )

      enddo

!================================================================================================================================== 
!Box(t)*q + Box(u)*q

      Coupl_Ac = dconjg(VVSfSf(xtype,ytype,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))

!t
      Box_TQ = Box_TQ + 2*colorfactorQ*dreal(  (Coupl_Ac*(4*(g0Rt*(2*g1Rt*M1Lpt*(-(g2Lt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*xd0t + g3Lt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t))) + 
     &               g2Rt*(g3Rt*M0Lpt*(xc01t + xc02t - 4*xd00t + M2Lpt**2*xd0t + M3Lpt**2*xd0t - m4**2*xd0t) + 
     &                  g3Lt*M3Lpt*(xc01t + xc04t - 4*xd00t + M0Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t))) + 
     &            g1Lt*(-2*g2Rt*M2Lpt*(g3Lt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Rt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &               g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(xc01t + xc03t - 4*xd00t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - m3**2*xd0t) + 
     &                  g3Lt*(2*xb01t + 2*xb02t - 8*xc001t + 2*M0Lpt**2*xc01t - m1**2*xc01t - m2**2*xc01t + sqrtS**2*xc01t - m1**2*xc02t + 2*M1Lpt**2*xc02t - m3**2*xc02t + tman*xc02t - 
     &                     m2**2*xc03t + 2*M3Lpt**2*xc03t - m4**2*xc03t + tman*xc03t + 2*M1Lpt**2*xc04t + m3**2*xc04t - m4**2*xc04t + sqrtS**2*xc04t + 4*m3**2*xc11t + 2*m3**2*xc21t - 
     &                     2*m4**2*xc21t + 2*sqrtS**2*xc21t - 8*M0Lpt**2*xd00t + 4*m1**2*xd00t + 4*m2**2*xd00t - 4*sqrtS**2*xd00t - M1Lpt**2*m2**2*xd0t + 2*M0Lpt**2*M2Lpt**2*xd0t - 
     &                     m1**2*M2Lpt**2*xd0t - m2**2*M2Lpt**2*xd0t - M0Lpt**2*m3**2*xd0t + m2**2*m3**2*xd0t - m1**2*M3Lpt**2*xd0t + 2*M1Lpt**2*M3Lpt**2*xd0t - m3**2*M3Lpt**2*xd0t - 
     &                     M0Lpt**2*m4**2*xd0t + m1**2*m4**2*xd0t - M1Lpt**2*m4**2*xd0t + M0Lpt**2*sqrtS**2*xd0t + M2Lpt**2*sqrtS**2*xd0t + M1Lpt**2*tman*xd0t + M3Lpt**2*tman*xd0t - 
     &                     sqrtS**2*tman*xd0t)))) + g0Lt*(2*g1Lt*M1Lpt*(-(g2Rt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*xd0t + g3Rt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t))) + 
     &               g2Lt*(g3Lt*M0Lpt*(xc01t + xc02t - 4*xd00t + M2Lpt**2*xd0t + M3Lpt**2*xd0t - m4**2*xd0t) + 
     &                  g3Rt*M3Lpt*(xc01t + xc04t - 4*xd00t + M0Lpt**2*xd0t + M2Lpt**2*xd0t - tman*xd0t))) + 
     &            g1Rt*(-2*g2Lt*M2Lpt*(g3Rt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Lt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &               g2Rt*(2*g3Lt*M0Lpt*M3Lpt*(xc01t + xc03t - 4*xd00t + M1Lpt**2*xd0t + M2Lpt**2*xd0t - m3**2*xd0t) + 
     &                  g3Rt*(2*xb01t + 2*xb02t - 8*xc001t + 2*M0Lpt**2*xc01t - m1**2*xc01t - m2**2*xc01t + sqrtS**2*xc01t - m1**2*xc02t + 2*M1Lpt**2*xc02t - m3**2*xc02t + tman*xc02t - 
     &                     m2**2*xc03t + 2*M3Lpt**2*xc03t - m4**2*xc03t + tman*xc03t + 2*M1Lpt**2*xc04t + m3**2*xc04t - m4**2*xc04t + sqrtS**2*xc04t + 4*m3**2*xc11t + 2*m3**2*xc21t - 
     &                     2*m4**2*xc21t + 2*sqrtS**2*xc21t - 8*M0Lpt**2*xd00t + 4*m1**2*xd00t + 4*m2**2*xd00t - 4*sqrtS**2*xd00t - M1Lpt**2*m2**2*xd0t + 2*M0Lpt**2*M2Lpt**2*xd0t - 
     &                     m1**2*M2Lpt**2*xd0t - m2**2*M2Lpt**2*xd0t - M0Lpt**2*m3**2*xd0t + m2**2*m3**2*xd0t - m1**2*M3Lpt**2*xd0t + 2*M1Lpt**2*M3Lpt**2*xd0t - m3**2*M3Lpt**2*xd0t - 
     &                     M0Lpt**2*m4**2*xd0t + m1**2*m4**2*xd0t - M1Lpt**2*m4**2*xd0t + M0Lpt**2*sqrtS**2*xd0t + M2Lpt**2*sqrtS**2*xd0t + M1Lpt**2*tman*xd0t + M3Lpt**2*tman*xd0t - 
     &                     sqrtS**2*tman*xd0t))))) - 2*g0Rt*(2*g1Rt*g2Rt*M1Lpt*
     &          (g3Lt*M3Lpt*(2*tman*xd12t + 2*m2**2*xd13t - 2*sqrtS**2*xd13t + tman*xd1t + 2*tman*xd22t + 4*m2**2*xd23t - 2*sqrtS**2*xd23t - 2*uman*xd23t + 2*tman*xd2t + 2*m2**2*xd33t + 
     &               2*m2**2*xd3t - sqrtS**2*xd3t - uman*xd3t + m3**2*(-2*xd12t - xd1t + 2*xd23t + xd3t) + m1**2*(2*xd11t + 2*xd12t + 2*xd13t + xd1t + 2*xd23t + xd3t)) + 
     &            g3Rt*M0Lpt*(m3**2*xd0t - sqrtS**2*xd0t - uman*xd0t - 2*m3**2*xd12t + 2*tman*xd12t - 2*sqrtS**2*xd13t - m3**2*xd1t - sqrtS**2*xd1t + tman*xd1t + 2*tman*xd22t + 
     &               2*m3**2*xd23t - 2*sqrtS**2*xd23t - 2*uman*xd23t + m3**2*xd2t - sqrtS**2*xd2t + 2*tman*xd2t - uman*xd2t + m3**2*xd3t - sqrtS**2*xd3t - uman*xd3t + 
     &               m1**2*(xd0t + 2*xd11t + 2*xd12t + 2*xd13t + 2*xd1t + 2*xd23t + xd2t + xd3t) + m2**2*(2*xd0t + 2*xd13t + xd1t + 4*xd23t + 2*xd2t + 2*xd33t + 4*xd3t))) + 
     &         g1Lt*g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(tman*xd0t + 2*tman*xd12t + 2*m2**2*xd13t - 2*sqrtS**2*xd13t + tman*xd1t + 2*tman*xd22t + 4*m2**2*xd23t - 2*sqrtS**2*xd23t - 2*uman*xd23t + 
     &               3*tman*xd2t + 2*m2**2*xd33t - m3**2*(xd0t + 2*xd12t + xd1t - 2*xd23t + xd2t - xd3t) + 3*m2**2*xd3t - 2*sqrtS**2*xd3t - uman*xd3t + 
     &               m1**2*(xd0t + 2*xd11t + 2*xd12t + 2*xd13t + 3*xd1t + 2*xd23t + xd2t + 2*xd3t)) + 
     &            g3Lt*(tman*xc02t + 2*m2**2*xc03t - sqrtS**2*xc03t - uman*xc03t + m2**2*xc04t - tman*xc04t - uman*xc04t + m2**2*xc11t - tman*xc11t - uman*xc11t + 4*m2**2*xc121t - 
     &               4*tman*xc121t - 4*uman*xc121t + 2*m2**2*xc13t - sqrtS**2*xc13t - uman*xc13t + 2*m2**2*xc21t + 2*sqrtS**2*xc21t - 2*tman*xc21t - 2*uman*xc21t + 4*sqrtS**2*xc221t + 
     &               tman*xc22t + 2*m2**2*xc23t + 2*M0Lpt**2*m2**2*xd0t + 2*M1Lpt**2*m2**2*xd0t - M0Lpt**2*sqrtS**2*xd0t - M1Lpt**2*sqrtS**2*xd0t + M0Lpt**2*tman*xd0t - m2**2*tman*xd0t + 
     &               M3Lpt**2*tman*xd0t - M0Lpt**2*uman*xd0t - M1Lpt**2*uman*xd0t + 4*M0Lpt**2*tman*xd12t - 2*m2**2*tman*xd12t + 2*sqrtS**2*tman*xd12t + 4*M0Lpt**2*m2**2*xd13t - 
     &               2*m2**4*xd13t - 4*M0Lpt**2*sqrtS**2*xd13t + 4*m2**2*sqrtS**2*xd13t - 2*sqrtS**4*xd13t + M0Lpt**2*m2**2*xd1t + M1Lpt**2*m2**2*xd1t - M0Lpt**2*sqrtS**2*xd1t - 
     &               M1Lpt**2*sqrtS**2*xd1t + 2*M0Lpt**2*tman*xd1t - m2**2*tman*xd1t + sqrtS**2*tman*xd1t + 4*M0Lpt**2*tman*xd22t - 2*m2**2*tman*xd22t + 2*sqrtS**2*tman*xd22t + 
     &               8*M0Lpt**2*m2**2*xd23t - 4*m2**4*xd23t - 4*M0Lpt**2*sqrtS**2*xd23t + 6*m2**2*sqrtS**2*xd23t - 2*sqrtS**4*xd23t - 4*M0Lpt**2*uman*xd23t + 2*m2**2*uman*xd23t - 
     &               2*sqrtS**2*uman*xd23t + 2*M0Lpt**2*m2**2*xd2t + 2*M1Lpt**2*m2**2*xd2t - M0Lpt**2*sqrtS**2*xd2t - M1Lpt**2*sqrtS**2*xd2t + 5*M0Lpt**2*tman*xd2t - 3*m2**2*tman*xd2t + 
     &               M3Lpt**2*tman*xd2t + 2*sqrtS**2*tman*xd2t - M0Lpt**2*uman*xd2t - M1Lpt**2*uman*xd2t + 4*M0Lpt**2*m2**2*xd33t - 2*m2**4*xd33t + 2*m2**2*sqrtS**2*xd33t + 
     &               7*M0Lpt**2*m2**2*xd3t + 2*M1Lpt**2*m2**2*xd3t - 3*m2**4*xd3t + m2**2*M3Lpt**2*xd3t - 3*M0Lpt**2*sqrtS**2*xd3t + 4*m2**2*sqrtS**2*xd3t - M3Lpt**2*sqrtS**2*xd3t - 
     &               sqrtS**4*xd3t - 2*M0Lpt**2*uman*xd3t + m2**2*uman*xd3t - sqrtS**2*uman*xd3t - m1**4*(xd0t + 2*xd11t + 2*xd12t + 2*xd13t + 2*xd1t + 2*xd23t + xd2t + xd3t) + 
     &               m3**2*(-xc02t + xc03t + 2*xc04t + 4*xc111t + 6*xc11t + 8*xc121t + xc13t + 4*xc21t - xc22t + M1Lpt**2*xd0t + m2**2*xd0t - M3Lpt**2*xd0t - 4*M0Lpt**2*xd12t + 
     &                  2*m2**2*xd12t - 2*sqrtS**2*xd12t - 2*M0Lpt**2*xd1t + m2**2*xd1t - sqrtS**2*xd1t + 4*M0Lpt**2*xd23t - 2*m2**2*xd23t + 2*sqrtS**2*xd23t + M1Lpt**2*xd2t + m2**2*xd2t - 
     &                  M3Lpt**2*xd2t + 2*M0Lpt**2*xd3t - m2**2*xd3t + sqrtS**2*xd3t) + 
     &               m1**2*(xc02t + xc03t + xc04t + xc11t + 4*xc121t + 2*xc12t + xc13t + 2*xc21t + xc22t + 2*M0Lpt**2*xd0t + M1Lpt**2*xd0t - 3*m2**2*xd0t - m3**2*xd0t + M3Lpt**2*xd0t + 
     &                  sqrtS**2*xd0t + uman*xd0t + 4*M0Lpt**2*xd11t - 2*m2**2*xd11t + 2*sqrtS**2*xd11t + 4*M0Lpt**2*xd12t - 2*m2**2*xd12t + 2*m3**2*xd12t + 2*sqrtS**2*xd12t - 2*tman*xd12t + 
     &                  4*M0Lpt**2*xd13t - 4*m2**2*xd13t + 4*sqrtS**2*xd13t + 5*M0Lpt**2*xd1t + M1Lpt**2*xd1t - 4*m2**2*xd1t + m3**2*xd1t + 2*M3Lpt**2*xd1t + 2*sqrtS**2*xd1t - tman*xd1t - 
     &                  2*tman*xd22t + 4*M0Lpt**2*xd23t - 6*m2**2*xd23t - 2*m3**2*xd23t + 4*sqrtS**2*xd23t + 2*uman*xd23t + 2*M0Lpt**2*xd2t + M1Lpt**2*xd2t - 3*m2**2*xd2t - m3**2*xd2t + 
     &                  M3Lpt**2*xd2t + sqrtS**2*xd2t - 2*tman*xd2t + uman*xd2t - 2*m2**2*xd33t + 3*M0Lpt**2*xd3t - 6*m2**2*xd3t - m3**2*xd3t + M3Lpt**2*xd3t + 2*sqrtS**2*xd3t + uman*xd3t))))
     &        - 2*g0Lt*(2*g1Lt*g2Lt*M1Lpt*(g3Rt*M3Lpt*(2*tman*xd12t + 2*m2**2*xd13t - 2*sqrtS**2*xd13t + tman*xd1t + 2*tman*xd22t + 4*m2**2*xd23t - 2*sqrtS**2*xd23t - 2*uman*xd23t + 
     &               2*tman*xd2t + 2*m2**2*xd33t + 2*m2**2*xd3t - sqrtS**2*xd3t - uman*xd3t + m3**2*(-2*xd12t - xd1t + 2*xd23t + xd3t) + 
     &               m1**2*(2*xd11t + 2*xd12t + 2*xd13t + xd1t + 2*xd23t + xd3t)) + 
     &            g3Lt*M0Lpt*(m3**2*xd0t - sqrtS**2*xd0t - uman*xd0t - 2*m3**2*xd12t + 2*tman*xd12t - 2*sqrtS**2*xd13t - m3**2*xd1t - sqrtS**2*xd1t + tman*xd1t + 2*tman*xd22t + 
     &               2*m3**2*xd23t - 2*sqrtS**2*xd23t - 2*uman*xd23t + m3**2*xd2t - sqrtS**2*xd2t + 2*tman*xd2t - uman*xd2t + m3**2*xd3t - sqrtS**2*xd3t - uman*xd3t + 
     &               m1**2*(xd0t + 2*xd11t + 2*xd12t + 2*xd13t + 2*xd1t + 2*xd23t + xd2t + xd3t) + m2**2*(2*xd0t + 2*xd13t + xd1t + 4*xd23t + 2*xd2t + 2*xd33t + 4*xd3t))) + 
     &         g1Rt*g2Rt*(2*g3Lt*M0Lpt*M3Lpt*(tman*xd0t + 2*tman*xd12t + 2*m2**2*xd13t - 2*sqrtS**2*xd13t + tman*xd1t + 2*tman*xd22t + 4*m2**2*xd23t - 2*sqrtS**2*xd23t - 2*uman*xd23t + 
     &               3*tman*xd2t + 2*m2**2*xd33t - m3**2*(xd0t + 2*xd12t + xd1t - 2*xd23t + xd2t - xd3t) + 3*m2**2*xd3t - 2*sqrtS**2*xd3t - uman*xd3t + 
     &               m1**2*(xd0t + 2*xd11t + 2*xd12t + 2*xd13t + 3*xd1t + 2*xd23t + xd2t + 2*xd3t)) + 
     &            g3Rt*(tman*xc02t + 2*m2**2*xc03t - sqrtS**2*xc03t - uman*xc03t + m2**2*xc04t - tman*xc04t - uman*xc04t + m2**2*xc11t - tman*xc11t - uman*xc11t + 4*m2**2*xc121t - 
     &               4*tman*xc121t - 4*uman*xc121t + 2*m2**2*xc13t - sqrtS**2*xc13t - uman*xc13t + 2*m2**2*xc21t + 2*sqrtS**2*xc21t - 2*tman*xc21t - 2*uman*xc21t + 4*sqrtS**2*xc221t + 
     &               tman*xc22t + 2*m2**2*xc23t + 2*M0Lpt**2*m2**2*xd0t + 2*M1Lpt**2*m2**2*xd0t - M0Lpt**2*sqrtS**2*xd0t - M1Lpt**2*sqrtS**2*xd0t + M0Lpt**2*tman*xd0t - m2**2*tman*xd0t + 
     &               M3Lpt**2*tman*xd0t - M0Lpt**2*uman*xd0t - M1Lpt**2*uman*xd0t + 4*M0Lpt**2*tman*xd12t - 2*m2**2*tman*xd12t + 2*sqrtS**2*tman*xd12t + 4*M0Lpt**2*m2**2*xd13t - 
     &               2*m2**4*xd13t - 4*M0Lpt**2*sqrtS**2*xd13t + 4*m2**2*sqrtS**2*xd13t - 2*sqrtS**4*xd13t + M0Lpt**2*m2**2*xd1t + M1Lpt**2*m2**2*xd1t - M0Lpt**2*sqrtS**2*xd1t - 
     &               M1Lpt**2*sqrtS**2*xd1t + 2*M0Lpt**2*tman*xd1t - m2**2*tman*xd1t + sqrtS**2*tman*xd1t + 4*M0Lpt**2*tman*xd22t - 2*m2**2*tman*xd22t + 2*sqrtS**2*tman*xd22t + 
     &               8*M0Lpt**2*m2**2*xd23t - 4*m2**4*xd23t - 4*M0Lpt**2*sqrtS**2*xd23t + 6*m2**2*sqrtS**2*xd23t - 2*sqrtS**4*xd23t - 4*M0Lpt**2*uman*xd23t + 2*m2**2*uman*xd23t - 
     &               2*sqrtS**2*uman*xd23t + 2*M0Lpt**2*m2**2*xd2t + 2*M1Lpt**2*m2**2*xd2t - M0Lpt**2*sqrtS**2*xd2t - M1Lpt**2*sqrtS**2*xd2t + 5*M0Lpt**2*tman*xd2t - 3*m2**2*tman*xd2t + 
     &               M3Lpt**2*tman*xd2t + 2*sqrtS**2*tman*xd2t - M0Lpt**2*uman*xd2t - M1Lpt**2*uman*xd2t + 4*M0Lpt**2*m2**2*xd33t - 2*m2**4*xd33t + 2*m2**2*sqrtS**2*xd33t + 
     &               7*M0Lpt**2*m2**2*xd3t + 2*M1Lpt**2*m2**2*xd3t - 3*m2**4*xd3t + m2**2*M3Lpt**2*xd3t - 3*M0Lpt**2*sqrtS**2*xd3t + 4*m2**2*sqrtS**2*xd3t - M3Lpt**2*sqrtS**2*xd3t - 
     &               sqrtS**4*xd3t - 2*M0Lpt**2*uman*xd3t + m2**2*uman*xd3t - sqrtS**2*uman*xd3t - m1**4*(xd0t + 2*xd11t + 2*xd12t + 2*xd13t + 2*xd1t + 2*xd23t + xd2t + xd3t) + 
     &               m3**2*(-xc02t + xc03t + 2*xc04t + 4*xc111t + 6*xc11t + 8*xc121t + xc13t + 4*xc21t - xc22t + M1Lpt**2*xd0t + m2**2*xd0t - M3Lpt**2*xd0t - 4*M0Lpt**2*xd12t + 
     &                  2*m2**2*xd12t - 2*sqrtS**2*xd12t - 2*M0Lpt**2*xd1t + m2**2*xd1t - sqrtS**2*xd1t + 4*M0Lpt**2*xd23t - 2*m2**2*xd23t + 2*sqrtS**2*xd23t + M1Lpt**2*xd2t + m2**2*xd2t - 
     &                  M3Lpt**2*xd2t + 2*M0Lpt**2*xd3t - m2**2*xd3t + sqrtS**2*xd3t) + 
     &               m1**2*(xc02t + xc03t + xc04t + xc11t + 4*xc121t + 2*xc12t + xc13t + 2*xc21t + xc22t + 2*M0Lpt**2*xd0t + M1Lpt**2*xd0t - 3*m2**2*xd0t - m3**2*xd0t + M3Lpt**2*xd0t + 
     &                  sqrtS**2*xd0t + uman*xd0t + 4*M0Lpt**2*xd11t - 2*m2**2*xd11t + 2*sqrtS**2*xd11t + 4*M0Lpt**2*xd12t - 2*m2**2*xd12t + 2*m3**2*xd12t + 2*sqrtS**2*xd12t - 2*tman*xd12t + 
     &                  4*M0Lpt**2*xd13t - 4*m2**2*xd13t + 4*sqrtS**2*xd13t + 5*M0Lpt**2*xd1t + M1Lpt**2*xd1t - 4*m2**2*xd1t + m3**2*xd1t + 2*M3Lpt**2*xd1t + 2*sqrtS**2*xd1t - tman*xd1t - 
     &                  2*tman*xd22t + 4*M0Lpt**2*xd23t - 6*m2**2*xd23t - 2*m3**2*xd23t + 4*sqrtS**2*xd23t + 2*uman*xd23t + 2*M0Lpt**2*xd2t + M1Lpt**2*xd2t - 3*m2**2*xd2t - m3**2*xd2t + 
     &                  M3Lpt**2*xd2t + sqrtS**2*xd2t - 2*tman*xd2t + uman*xd2t - 2*m2**2*xd33t + 3*M0Lpt**2*xd3t - 6*m2**2*xd3t - m3**2*xd3t + M3Lpt**2*xd3t + 2*sqrtS**2*xd3t + uman*xd3t))))
     &      ))/(32.*pi**2)  )

!u
      Box_UQ = Box_UQ + 2*colorfactorQ*dreal(  (Coupl_Ac*(4*(g0Ru*(2*g1Ru*M1Lpu*(-(g2Lu*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*xd0u + g3Lu*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u))) + 
     &               g2Ru*(g3Ru*M0Lpu*(xc01u + xc02u - 4*xd00u + M2Lpu**2*xd0u - m3**2*xd0u + M3Lpu**2*xd0u) + 
     &                  g3Lu*M3Lpu*(xc01u + xc04u - 4*xd00u + M0Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u))) + 
     &            g1Lu*(-2*g2Ru*M2Lpu*(g3Lu*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Ru*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &               g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(xc01u + xc03u - 4*xd00u + M1Lpu**2*xd0u + M2Lpu**2*xd0u - m4**2*xd0u) + 
     &                  g3Lu*(2*xb01u + 2*xb02u - 8*xc001u + 2*M0Lpu**2*xc01u - m1**2*xc01u - m2**2*xc01u + sqrtS**2*xc01u - m1**2*xc02u + 2*M1Lpu**2*xc02u - m4**2*xc02u + uman*xc02u - 
     &                     m2**2*xc03u - m3**2*xc03u + 2*M3Lpu**2*xc03u + uman*xc03u + 2*M1Lpu**2*xc04u - m3**2*xc04u + m4**2*xc04u + sqrtS**2*xc04u + 4*m4**2*xc11u - 2*m3**2*xc21u + 
     &                     2*m4**2*xc21u + 2*sqrtS**2*xc21u - 8*M0Lpu**2*xd00u + 4*m1**2*xd00u + 4*m2**2*xd00u - 4*sqrtS**2*xd00u - M1Lpu**2*m2**2*xd0u + 2*M0Lpu**2*M2Lpu**2*xd0u - 
     &                     m1**2*M2Lpu**2*xd0u - m2**2*M2Lpu**2*xd0u - M0Lpu**2*m3**2*xd0u + m1**2*m3**2*xd0u - M1Lpu**2*m3**2*xd0u - m1**2*M3Lpu**2*xd0u + 2*M1Lpu**2*M3Lpu**2*xd0u - 
     &                     M0Lpu**2*m4**2*xd0u + m2**2*m4**2*xd0u - M3Lpu**2*m4**2*xd0u + M0Lpu**2*sqrtS**2*xd0u + M2Lpu**2*sqrtS**2*xd0u + M1Lpu**2*uman*xd0u + M3Lpu**2*uman*xd0u - 
     &                     sqrtS**2*uman*xd0u)))) + g0Lu*(2*g1Lu*M1Lpu*(-(g2Ru*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*xd0u + g3Ru*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u))) + 
     &               g2Lu*(g3Lu*M0Lpu*(xc01u + xc02u - 4*xd00u + M2Lpu**2*xd0u - m3**2*xd0u + M3Lpu**2*xd0u) + 
     &                  g3Ru*M3Lpu*(xc01u + xc04u - 4*xd00u + M0Lpu**2*xd0u + M2Lpu**2*xd0u - uman*xd0u))) + 
     &            g1Ru*(-2*g2Lu*M2Lpu*(g3Ru*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Lu*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &               g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(xc01u + xc03u - 4*xd00u + M1Lpu**2*xd0u + M2Lpu**2*xd0u - m4**2*xd0u) + 
     &                  g3Ru*(2*xb01u + 2*xb02u - 8*xc001u + 2*M0Lpu**2*xc01u - m1**2*xc01u - m2**2*xc01u + sqrtS**2*xc01u - m1**2*xc02u + 2*M1Lpu**2*xc02u - m4**2*xc02u + uman*xc02u - 
     &                     m2**2*xc03u - m3**2*xc03u + 2*M3Lpu**2*xc03u + uman*xc03u + 2*M1Lpu**2*xc04u - m3**2*xc04u + m4**2*xc04u + sqrtS**2*xc04u + 4*m4**2*xc11u - 2*m3**2*xc21u + 
     &                     2*m4**2*xc21u + 2*sqrtS**2*xc21u - 8*M0Lpu**2*xd00u + 4*m1**2*xd00u + 4*m2**2*xd00u - 4*sqrtS**2*xd00u - M1Lpu**2*m2**2*xd0u + 2*M0Lpu**2*M2Lpu**2*xd0u - 
     &                     m1**2*M2Lpu**2*xd0u - m2**2*M2Lpu**2*xd0u - M0Lpu**2*m3**2*xd0u + m1**2*m3**2*xd0u - M1Lpu**2*m3**2*xd0u - m1**2*M3Lpu**2*xd0u + 2*M1Lpu**2*M3Lpu**2*xd0u - 
     &                     M0Lpu**2*m4**2*xd0u + m2**2*m4**2*xd0u - M3Lpu**2*m4**2*xd0u + M0Lpu**2*sqrtS**2*xd0u + M2Lpu**2*sqrtS**2*xd0u + M1Lpu**2*uman*xd0u + M3Lpu**2*uman*xd0u - 
     &                     sqrtS**2*uman*xd0u))))) - 2*g0Ru*(2*g1Ru*g2Ru*M1Lpu*
     &          (g3Lu*M3Lpu*(2*uman*xd12u + 2*m2**2*xd13u - 2*sqrtS**2*xd13u + uman*xd1u + 2*uman*xd22u + 4*m2**2*xd23u - 2*sqrtS**2*xd23u - 2*tman*xd23u + 2*uman*xd2u + 2*m2**2*xd33u + 
     &               2*m2**2*xd3u - sqrtS**2*xd3u - tman*xd3u + m4**2*(-2*xd12u - xd1u + 2*xd23u + xd3u) + m1**2*(2*xd11u + 2*xd12u + 2*xd13u + xd1u + 2*xd23u + xd3u)) + 
     &            g3Ru*M0Lpu*(m4**2*xd0u - sqrtS**2*xd0u - tman*xd0u - 2*m4**2*xd12u + 2*uman*xd12u - 2*sqrtS**2*xd13u - m4**2*xd1u - sqrtS**2*xd1u + uman*xd1u + 2*uman*xd22u + 
     &               2*m4**2*xd23u - 2*sqrtS**2*xd23u - 2*tman*xd23u + m4**2*xd2u - sqrtS**2*xd2u - tman*xd2u + 2*uman*xd2u + m4**2*xd3u - sqrtS**2*xd3u - tman*xd3u + 
     &               m1**2*(xd0u + 2*xd11u + 2*xd12u + 2*xd13u + 2*xd1u + 2*xd23u + xd2u + xd3u) + m2**2*(2*xd0u + 2*xd13u + xd1u + 4*xd23u + 2*xd2u + 2*xd33u + 4*xd3u))) + 
     &         g1Lu*g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(uman*xd0u + 2*uman*xd12u + 2*m2**2*xd13u - 2*sqrtS**2*xd13u + uman*xd1u + 2*uman*xd22u + 4*m2**2*xd23u - 2*sqrtS**2*xd23u - 2*tman*xd23u + 
     &               3*uman*xd2u + 2*m2**2*xd33u - m4**2*(xd0u + 2*xd12u + xd1u - 2*xd23u + xd2u - xd3u) + 3*m2**2*xd3u - 2*sqrtS**2*xd3u - tman*xd3u + 
     &               m1**2*(xd0u + 2*xd11u + 2*xd12u + 2*xd13u + 3*xd1u + 2*xd23u + xd2u + 2*xd3u)) + 
     &            g3Lu*(uman*xc02u + 2*m2**2*xc03u - sqrtS**2*xc03u - tman*xc03u + m2**2*xc04u - tman*xc04u - uman*xc04u + m2**2*xc11u - tman*xc11u - uman*xc11u + 4*m2**2*xc121u - 
     &               4*tman*xc121u - 4*uman*xc121u + 2*m2**2*xc13u - sqrtS**2*xc13u - tman*xc13u + 2*m2**2*xc21u + 2*sqrtS**2*xc21u - 2*tman*xc21u - 2*uman*xc21u + 4*sqrtS**2*xc221u + 
     &               uman*xc22u + 2*m2**2*xc23u + 2*M0Lpu**2*m2**2*xd0u + 2*M1Lpu**2*m2**2*xd0u - M0Lpu**2*sqrtS**2*xd0u - M1Lpu**2*sqrtS**2*xd0u - M0Lpu**2*tman*xd0u - M1Lpu**2*tman*xd0u + 
     &               M0Lpu**2*uman*xd0u - m2**2*uman*xd0u + M3Lpu**2*uman*xd0u + 4*M0Lpu**2*uman*xd12u - 2*m2**2*uman*xd12u + 2*sqrtS**2*uman*xd12u + 4*M0Lpu**2*m2**2*xd13u - 2*m2**4*xd13u - 
     &               4*M0Lpu**2*sqrtS**2*xd13u + 4*m2**2*sqrtS**2*xd13u - 2*sqrtS**4*xd13u + M0Lpu**2*m2**2*xd1u + M1Lpu**2*m2**2*xd1u - M0Lpu**2*sqrtS**2*xd1u - M1Lpu**2*sqrtS**2*xd1u + 
     &               2*M0Lpu**2*uman*xd1u - m2**2*uman*xd1u + sqrtS**2*uman*xd1u + 4*M0Lpu**2*uman*xd22u - 2*m2**2*uman*xd22u + 2*sqrtS**2*uman*xd22u + 8*M0Lpu**2*m2**2*xd23u - 
     &               4*m2**4*xd23u - 4*M0Lpu**2*sqrtS**2*xd23u + 6*m2**2*sqrtS**2*xd23u - 2*sqrtS**4*xd23u - 4*M0Lpu**2*tman*xd23u + 2*m2**2*tman*xd23u - 2*sqrtS**2*tman*xd23u + 
     &               2*M0Lpu**2*m2**2*xd2u + 2*M1Lpu**2*m2**2*xd2u - M0Lpu**2*sqrtS**2*xd2u - M1Lpu**2*sqrtS**2*xd2u - M0Lpu**2*tman*xd2u - M1Lpu**2*tman*xd2u + 5*M0Lpu**2*uman*xd2u - 
     &               3*m2**2*uman*xd2u + M3Lpu**2*uman*xd2u + 2*sqrtS**2*uman*xd2u + 4*M0Lpu**2*m2**2*xd33u - 2*m2**4*xd33u + 2*m2**2*sqrtS**2*xd33u + 7*M0Lpu**2*m2**2*xd3u + 
     &               2*M1Lpu**2*m2**2*xd3u - 3*m2**4*xd3u + m2**2*M3Lpu**2*xd3u - 3*M0Lpu**2*sqrtS**2*xd3u + 4*m2**2*sqrtS**2*xd3u - M3Lpu**2*sqrtS**2*xd3u - sqrtS**4*xd3u - 
     &               2*M0Lpu**2*tman*xd3u + m2**2*tman*xd3u - sqrtS**2*tman*xd3u - m1**4*(xd0u + 2*xd11u + 2*xd12u + 2*xd13u + 2*xd1u + 2*xd23u + xd2u + xd3u) + 
     &               m4**2*(-xc02u + xc03u + 2*xc04u + 4*xc111u + 6*xc11u + 8*xc121u + xc13u + 4*xc21u - xc22u + M1Lpu**2*xd0u + m2**2*xd0u - M3Lpu**2*xd0u - 4*M0Lpu**2*xd12u + 
     &                  2*m2**2*xd12u - 2*sqrtS**2*xd12u - 2*M0Lpu**2*xd1u + m2**2*xd1u - sqrtS**2*xd1u + 4*M0Lpu**2*xd23u - 2*m2**2*xd23u + 2*sqrtS**2*xd23u + M1Lpu**2*xd2u + m2**2*xd2u - 
     &                  M3Lpu**2*xd2u + 2*M0Lpu**2*xd3u - m2**2*xd3u + sqrtS**2*xd3u) + 
     &               m1**2*(xc02u + xc03u + xc04u + xc11u + 4*xc121u + 2*xc12u + xc13u + 2*xc21u + xc22u + 2*M0Lpu**2*xd0u + M1Lpu**2*xd0u - 3*m2**2*xd0u + M3Lpu**2*xd0u - m4**2*xd0u + 
     &                  sqrtS**2*xd0u + tman*xd0u + 4*M0Lpu**2*xd11u - 2*m2**2*xd11u + 2*sqrtS**2*xd11u + 4*M0Lpu**2*xd12u - 2*m2**2*xd12u + 2*m4**2*xd12u + 2*sqrtS**2*xd12u - 2*uman*xd12u + 
     &                  4*M0Lpu**2*xd13u - 4*m2**2*xd13u + 4*sqrtS**2*xd13u + 5*M0Lpu**2*xd1u + M1Lpu**2*xd1u - 4*m2**2*xd1u + 2*M3Lpu**2*xd1u + m4**2*xd1u + 2*sqrtS**2*xd1u - uman*xd1u - 
     &                  2*uman*xd22u + 4*M0Lpu**2*xd23u - 6*m2**2*xd23u - 2*m4**2*xd23u + 4*sqrtS**2*xd23u + 2*tman*xd23u + 2*M0Lpu**2*xd2u + M1Lpu**2*xd2u - 3*m2**2*xd2u + M3Lpu**2*xd2u - 
     &                  m4**2*xd2u + sqrtS**2*xd2u + tman*xd2u - 2*uman*xd2u - 2*m2**2*xd33u + 3*M0Lpu**2*xd3u - 6*m2**2*xd3u + M3Lpu**2*xd3u - m4**2*xd3u + 2*sqrtS**2*xd3u + tman*xd3u)))) - 
     &      2*g0Lu*(2*g1Lu*g2Lu*M1Lpu*(g3Ru*M3Lpu*(2*uman*xd12u + 2*m2**2*xd13u - 2*sqrtS**2*xd13u + uman*xd1u + 2*uman*xd22u + 4*m2**2*xd23u - 2*sqrtS**2*xd23u - 2*tman*xd23u + 
     &               2*uman*xd2u + 2*m2**2*xd33u + 2*m2**2*xd3u - sqrtS**2*xd3u - tman*xd3u + m4**2*(-2*xd12u - xd1u + 2*xd23u + xd3u) + 
     &               m1**2*(2*xd11u + 2*xd12u + 2*xd13u + xd1u + 2*xd23u + xd3u)) + 
     &            g3Lu*M0Lpu*(m4**2*xd0u - sqrtS**2*xd0u - tman*xd0u - 2*m4**2*xd12u + 2*uman*xd12u - 2*sqrtS**2*xd13u - m4**2*xd1u - sqrtS**2*xd1u + uman*xd1u + 2*uman*xd22u + 
     &               2*m4**2*xd23u - 2*sqrtS**2*xd23u - 2*tman*xd23u + m4**2*xd2u - sqrtS**2*xd2u - tman*xd2u + 2*uman*xd2u + m4**2*xd3u - sqrtS**2*xd3u - tman*xd3u + 
     &               m1**2*(xd0u + 2*xd11u + 2*xd12u + 2*xd13u + 2*xd1u + 2*xd23u + xd2u + xd3u) + m2**2*(2*xd0u + 2*xd13u + xd1u + 4*xd23u + 2*xd2u + 2*xd33u + 4*xd3u))) + 
     &         g1Ru*g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(uman*xd0u + 2*uman*xd12u + 2*m2**2*xd13u - 2*sqrtS**2*xd13u + uman*xd1u + 2*uman*xd22u + 4*m2**2*xd23u - 2*sqrtS**2*xd23u - 2*tman*xd23u + 
     &               3*uman*xd2u + 2*m2**2*xd33u - m4**2*(xd0u + 2*xd12u + xd1u - 2*xd23u + xd2u - xd3u) + 3*m2**2*xd3u - 2*sqrtS**2*xd3u - tman*xd3u + 
     &               m1**2*(xd0u + 2*xd11u + 2*xd12u + 2*xd13u + 3*xd1u + 2*xd23u + xd2u + 2*xd3u)) + 
     &            g3Ru*(uman*xc02u + 2*m2**2*xc03u - sqrtS**2*xc03u - tman*xc03u + m2**2*xc04u - tman*xc04u - uman*xc04u + m2**2*xc11u - tman*xc11u - uman*xc11u + 4*m2**2*xc121u - 
     &               4*tman*xc121u - 4*uman*xc121u + 2*m2**2*xc13u - sqrtS**2*xc13u - tman*xc13u + 2*m2**2*xc21u + 2*sqrtS**2*xc21u - 2*tman*xc21u - 2*uman*xc21u + 4*sqrtS**2*xc221u + 
     &               uman*xc22u + 2*m2**2*xc23u + 2*M0Lpu**2*m2**2*xd0u + 2*M1Lpu**2*m2**2*xd0u - M0Lpu**2*sqrtS**2*xd0u - M1Lpu**2*sqrtS**2*xd0u - M0Lpu**2*tman*xd0u - M1Lpu**2*tman*xd0u + 
     &               M0Lpu**2*uman*xd0u - m2**2*uman*xd0u + M3Lpu**2*uman*xd0u + 4*M0Lpu**2*uman*xd12u - 2*m2**2*uman*xd12u + 2*sqrtS**2*uman*xd12u + 4*M0Lpu**2*m2**2*xd13u - 2*m2**4*xd13u - 
     &               4*M0Lpu**2*sqrtS**2*xd13u + 4*m2**2*sqrtS**2*xd13u - 2*sqrtS**4*xd13u + M0Lpu**2*m2**2*xd1u + M1Lpu**2*m2**2*xd1u - M0Lpu**2*sqrtS**2*xd1u - M1Lpu**2*sqrtS**2*xd1u + 
     &               2*M0Lpu**2*uman*xd1u - m2**2*uman*xd1u + sqrtS**2*uman*xd1u + 4*M0Lpu**2*uman*xd22u - 2*m2**2*uman*xd22u + 2*sqrtS**2*uman*xd22u + 8*M0Lpu**2*m2**2*xd23u - 
     &               4*m2**4*xd23u - 4*M0Lpu**2*sqrtS**2*xd23u + 6*m2**2*sqrtS**2*xd23u - 2*sqrtS**4*xd23u - 4*M0Lpu**2*tman*xd23u + 2*m2**2*tman*xd23u - 2*sqrtS**2*tman*xd23u + 
     &               2*M0Lpu**2*m2**2*xd2u + 2*M1Lpu**2*m2**2*xd2u - M0Lpu**2*sqrtS**2*xd2u - M1Lpu**2*sqrtS**2*xd2u - M0Lpu**2*tman*xd2u - M1Lpu**2*tman*xd2u + 5*M0Lpu**2*uman*xd2u - 
     &               3*m2**2*uman*xd2u + M3Lpu**2*uman*xd2u + 2*sqrtS**2*uman*xd2u + 4*M0Lpu**2*m2**2*xd33u - 2*m2**4*xd33u + 2*m2**2*sqrtS**2*xd33u + 7*M0Lpu**2*m2**2*xd3u + 
     &               2*M1Lpu**2*m2**2*xd3u - 3*m2**4*xd3u + m2**2*M3Lpu**2*xd3u - 3*M0Lpu**2*sqrtS**2*xd3u + 4*m2**2*sqrtS**2*xd3u - M3Lpu**2*sqrtS**2*xd3u - sqrtS**4*xd3u - 
     &               2*M0Lpu**2*tman*xd3u + m2**2*tman*xd3u - sqrtS**2*tman*xd3u - m1**4*(xd0u + 2*xd11u + 2*xd12u + 2*xd13u + 2*xd1u + 2*xd23u + xd2u + xd3u) + 
     &               m4**2*(-xc02u + xc03u + 2*xc04u + 4*xc111u + 6*xc11u + 8*xc121u + xc13u + 4*xc21u - xc22u + M1Lpu**2*xd0u + m2**2*xd0u - M3Lpu**2*xd0u - 4*M0Lpu**2*xd12u + 
     &                  2*m2**2*xd12u - 2*sqrtS**2*xd12u - 2*M0Lpu**2*xd1u + m2**2*xd1u - sqrtS**2*xd1u + 4*M0Lpu**2*xd23u - 2*m2**2*xd23u + 2*sqrtS**2*xd23u + M1Lpu**2*xd2u + m2**2*xd2u - 
     &                  M3Lpu**2*xd2u + 2*M0Lpu**2*xd3u - m2**2*xd3u + sqrtS**2*xd3u) + 
     &               m1**2*(xc02u + xc03u + xc04u + xc11u + 4*xc121u + 2*xc12u + xc13u + 2*xc21u + xc22u + 2*M0Lpu**2*xd0u + M1Lpu**2*xd0u - 3*m2**2*xd0u + M3Lpu**2*xd0u - m4**2*xd0u + 
     &                  sqrtS**2*xd0u + tman*xd0u + 4*M0Lpu**2*xd11u - 2*m2**2*xd11u + 2*sqrtS**2*xd11u + 4*M0Lpu**2*xd12u - 2*m2**2*xd12u + 2*m4**2*xd12u + 2*sqrtS**2*xd12u - 2*uman*xd12u + 
     &                  4*M0Lpu**2*xd13u - 4*m2**2*xd13u + 4*sqrtS**2*xd13u + 5*M0Lpu**2*xd1u + M1Lpu**2*xd1u - 4*m2**2*xd1u + 2*M3Lpu**2*xd1u + m4**2*xd1u + 2*sqrtS**2*xd1u - uman*xd1u - 
     &                  2*uman*xd22u + 4*M0Lpu**2*xd23u - 6*m2**2*xd23u - 2*m4**2*xd23u + 4*sqrtS**2*xd23u + 2*tman*xd23u + 2*M0Lpu**2*xd2u + M1Lpu**2*xd2u - 3*m2**2*xd2u + M3Lpu**2*xd2u - 
     &                  m4**2*xd2u + sqrtS**2*xd2u + tman*xd2u - 2*uman*xd2u - 2*m2**2*xd33u + 3*M0Lpu**2*xd3u - 6*m2**2*xd3u + M3Lpu**2*xd3u - m4**2*xd3u + 2*sqrtS**2*xd3u + tman*xd3u))))))/
     &  (32.*pi**2)  )

!==================================================================================================================================  
!Box(t)*s(Higgs) + Box(u)*s(Higgs)  
       
      do j=sumstartH, sumendH

      Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
      Coupl_Bc = dconjg(LLR(xtype,ytype,j))
      propac = dconjg(DenSmH(j))
        
      Box_TS = Box_TS + 2*colorfactorSH*dreal(  (Coupl_Ac*Coupl_Bc*(g0Rt*(2*g1Rt*M1Lpt*(2*g2Lt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*xd0t + g3Lt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t)) + 
     &            g2Rt*(g3Lt*M3Lpt*(-2*xc01t - 2*xc04t + 8*xd00t - 2*M0Lpt**2*xd0t - 2*M2Lpt**2*xd0t + tman*(2*xd0t + 2*xd12t + xd1t + 2*(xd22t + xd2t)) - 
     &                  (sqrtS**2 + uman)*(2*xd23t + xd3t) + m3**2*(-2*xd12t - xd1t + 2*xd23t + xd3t) + m1**2*(2*xd11t + 2*xd12t + 2*xd13t + xd1t + 2*xd23t + xd3t) + 
     &                  2*(-(sqrtS**2*xd13t) + m2**2*(xd13t + 2*xd23t + xd33t + xd3t))) + 
     &               g3Rt*M0Lpt*(-2*xc01t - 2*xc02t + 8*xd00t - 2*M2Lpt**2*xd0t + m3**2*xd0t - 2*M3Lpt**2*xd0t + 2*m4**2*xd0t - sqrtS**2*xd0t - uman*xd0t - 2*m3**2*xd12t + 2*tman*xd12t - 
     &                  2*sqrtS**2*xd13t - m3**2*xd1t - sqrtS**2*xd1t + tman*xd1t + 2*tman*xd22t + 2*m3**2*xd23t - 2*sqrtS**2*xd23t - 2*uman*xd23t + m3**2*xd2t - sqrtS**2*xd2t + 
     &                  2*tman*xd2t - uman*xd2t + (m3**2 - sqrtS**2 - uman)*xd3t + m1**2*(xd0t + 2*xd11t + 2*xd12t + 2*xd13t + 2*xd1t + 2*xd23t + xd2t + xd3t) + 
     &                  m2**2*(2*xd0t + 2*xd13t + xd1t + 2*(2*xd23t + xd2t + xd33t + 2*xd3t))))) + 
     &         g1Lt*(4*g2Rt*M2Lpt*(g3Lt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Rt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &            g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(-2*xc01t - 2*xc03t + 8*xd00t - 2*M1Lpt**2*xd0t - 2*M2Lpt**2*xd0t + m3**2*xd0t + tman*xd0t - 2*m3**2*xd12t + 2*tman*xd12t - 2*sqrtS**2*xd13t - 
     &                  m3**2*xd1t + tman*xd1t + 2*tman*xd22t + 2*m3**2*xd23t - 2*sqrtS**2*xd23t - 2*uman*xd23t - m3**2*xd2t + 3*tman*xd2t + (m3**2 - 2*sqrtS**2 - uman)*xd3t + 
     &                  m1**2*(xd0t + 2*xd11t + 2*xd12t + 2*xd13t + 3*xd1t + 2*xd23t + xd2t + 2*xd3t) + m2**2*(2*(xd13t + 2*xd23t + xd33t) + 3*xd3t)) + 
     &               g3Lt*(-4*xb01t - 4*xb02t + 16*xc001t + 2*m2**2*xc01t - 2*sqrtS**2*xc01t - 4*M1Lpt**2*xc02t + m3**2*xc02t - tman*xc02t + 4*m2**2*xc03t + m3**2*xc03t - 4*M3Lpt**2*xc03t + 
     &                  2*m4**2*xc03t - sqrtS**2*xc03t - 2*tman*xc03t - uman*xc03t - 4*M1Lpt**2*xc04t + m2**2*xc04t + 2*m4**2*xc04t - 2*sqrtS**2*xc04t - tman*xc04t - uman*xc04t + 
     &                  4*m3**2*xc111t + m2**2*xc11t - 2*m3**2*xc11t - tman*xc11t - uman*xc11t + 4*m2**2*xc121t + 8*m3**2*xc121t - 4*tman*xc121t - 4*uman*xc121t + 2*m2**2*xc13t + 
     &                  m3**2*xc13t - sqrtS**2*xc13t - uman*xc13t + 2*m2**2*xc21t + 4*m4**2*xc21t - 2*sqrtS**2*xc21t - 2*tman*xc21t - 2*uman*xc21t + 4*sqrtS**2*xc221t - m3**2*xc22t + 
     &                  tman*xc22t + 2*m2**2*xc23t - 8*m2**2*xd00t + 8*sqrtS**2*xd00t + 4*M1Lpt**2*m2**2*xd0t + 2*m2**2*M2Lpt**2*xd0t + M1Lpt**2*m3**2*xd0t - m2**2*m3**2*xd0t - 
     &                  4*M1Lpt**2*M3Lpt**2*xd0t + m3**2*M3Lpt**2*xd0t + 2*M1Lpt**2*m4**2*xd0t - M1Lpt**2*sqrtS**2*xd0t - 2*M2Lpt**2*sqrtS**2*xd0t - 2*M1Lpt**2*tman*xd0t - m2**2*tman*xd0t - 
     &                  M3Lpt**2*tman*xd0t + 2*sqrtS**2*tman*xd0t - M1Lpt**2*uman*xd0t + 2*m2**2*m3**2*xd12t - 2*m3**2*sqrtS**2*xd12t - 2*m2**2*tman*xd12t + 2*sqrtS**2*tman*xd12t - 
     &                  2*m2**4*xd13t + 4*m2**2*sqrtS**2*xd13t - 2*sqrtS**4*xd13t + M1Lpt**2*m2**2*xd1t + m2**2*m3**2*xd1t - M1Lpt**2*sqrtS**2*xd1t - m3**2*sqrtS**2*xd1t - m2**2*tman*xd1t + 
     &                  sqrtS**2*tman*xd1t - 2*m2**2*tman*xd22t + 2*sqrtS**2*tman*xd22t - 4*m2**4*xd23t - 2*m2**2*m3**2*xd23t + 6*m2**2*sqrtS**2*xd23t + 2*m3**2*sqrtS**2*xd23t - 
     &                  2*sqrtS**4*xd23t + 2*m2**2*uman*xd23t - 2*sqrtS**2*uman*xd23t + 2*M1Lpt**2*m2**2*xd2t + M1Lpt**2*m3**2*xd2t + m2**2*m3**2*xd2t - m3**2*M3Lpt**2*xd2t - 
     &                  M1Lpt**2*sqrtS**2*xd2t - 3*m2**2*tman*xd2t + M3Lpt**2*tman*xd2t + 2*sqrtS**2*tman*xd2t - M1Lpt**2*uman*xd2t - 2*m2**4*xd33t + 2*m2**2*sqrtS**2*xd33t + 
     &                  (2*M1Lpt**2*m2**2 - (m2 - sqrtS)*(m2 + sqrtS)*(3*m2**2 + m3**2 - M3Lpt**2 - sqrtS**2 - uman))*xd3t - 
     &                  m1**4*(xd0t + 2*xd11t + 2*xd12t + 2*xd13t + 2*xd1t + 2*xd23t + xd2t + xd3t) + 
     &                  m1**2*(2*xc01t + 3*xc02t + xc03t + xc04t + xc11t + 4*xc121t + 2*xc12t + xc13t + 2*xc21t + xc22t - 8*xd00t + 2*M2Lpt**2*xd0t - m3**2*xd0t + 3*M3Lpt**2*xd0t - 
     &                     2*m4**2*xd0t + sqrtS**2*xd0t + uman*xd0t + 2*sqrtS**2*xd11t + 2*m3**2*xd12t + 2*sqrtS**2*xd12t - 2*tman*xd12t + 4*sqrtS**2*xd13t + m3**2*xd1t + 2*M3Lpt**2*xd1t + 
     &                     2*sqrtS**2*xd1t - tman*xd1t - 2*tman*xd22t - 2*m3**2*xd23t + 4*sqrtS**2*xd23t + 2*uman*xd23t - m3**2*xd2t + M3Lpt**2*xd2t + sqrtS**2*xd2t - 2*tman*xd2t + 
     &                     uman*xd2t + M1Lpt**2*(xd0t + xd1t + xd2t) - m3**2*xd3t + M3Lpt**2*xd3t + 2*sqrtS**2*xd3t + uman*xd3t - 
     &                     m2**2*(3*xd0t + 2*xd11t + 2*xd12t + 4*xd13t + 4*xd1t + 6*xd23t + 3*xd2t + 2*xd33t + 6*xd3t)) + 
     &                  M0Lpt**2*(-4*xc01t + 16*xd00t - 4*M2Lpt**2*xd0t + 2*m3**2*xd0t + 2*m4**2*xd0t - 3*sqrtS**2*xd0t + tman*xd0t - uman*xd0t - 4*m3**2*xd12t + 4*tman*xd12t - 
     &                     4*sqrtS**2*xd13t - 2*m3**2*xd1t - sqrtS**2*xd1t + 2*tman*xd1t + 4*tman*xd22t + 4*m3**2*xd23t - 4*sqrtS**2*xd23t - 4*uman*xd23t - sqrtS**2*xd2t + 5*tman*xd2t - 
     &                     uman*xd2t + 2*m3**2*xd3t - 3*sqrtS**2*xd3t - 2*uman*xd3t + m1**2*(2*xd0t + 4*xd11t + 4*xd12t + 4*xd13t + 5*xd1t + 4*xd23t + 2*xd2t + 3*xd3t) + 
     &                     m2**2*(2*xd0t + 4*xd13t + xd1t + 8*xd23t + 2*xd2t + 4*xd33t + 7*xd3t)))))) + 
     &      g0Lt*(2*g1Lt*M1Lpt*(2*g2Rt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*xd0t + g3Rt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t)) + 
     &            g2Lt*(g3Rt*M3Lpt*(-2*xc01t - 2*xc04t + 8*xd00t - 2*M0Lpt**2*xd0t - 2*M2Lpt**2*xd0t + tman*(2*xd0t + 2*xd12t + xd1t + 2*(xd22t + xd2t)) - 
     &                  (sqrtS**2 + uman)*(2*xd23t + xd3t) + m3**2*(-2*xd12t - xd1t + 2*xd23t + xd3t) + m1**2*(2*xd11t + 2*xd12t + 2*xd13t + xd1t + 2*xd23t + xd3t) + 
     &                  2*(-(sqrtS**2*xd13t) + m2**2*(xd13t + 2*xd23t + xd33t + xd3t))) + 
     &               g3Lt*M0Lpt*(-2*xc01t - 2*xc02t + 8*xd00t - 2*M2Lpt**2*xd0t + m3**2*xd0t - 2*M3Lpt**2*xd0t + 2*m4**2*xd0t - sqrtS**2*xd0t - uman*xd0t - 2*m3**2*xd12t + 2*tman*xd12t - 
     &                  2*sqrtS**2*xd13t - m3**2*xd1t - sqrtS**2*xd1t + tman*xd1t + 2*tman*xd22t + 2*m3**2*xd23t - 2*sqrtS**2*xd23t - 2*uman*xd23t + m3**2*xd2t - sqrtS**2*xd2t + 
     &                  2*tman*xd2t - uman*xd2t + (m3**2 - sqrtS**2 - uman)*xd3t + m1**2*(xd0t + 2*xd11t + 2*xd12t + 2*xd13t + 2*xd1t + 2*xd23t + xd2t + xd3t) + 
     &                  m2**2*(2*xd0t + 2*xd13t + xd1t + 2*(2*xd23t + xd2t + xd33t + 2*xd3t))))) + 
     &         g1Rt*(4*g2Lt*M2Lpt*(g3Rt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Lt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &            g2Rt*(2*g3Lt*M0Lpt*M3Lpt*(-2*xc01t - 2*xc03t + 8*xd00t - 2*M1Lpt**2*xd0t - 2*M2Lpt**2*xd0t + m3**2*xd0t + tman*xd0t - 2*m3**2*xd12t + 2*tman*xd12t - 2*sqrtS**2*xd13t - 
     &                  m3**2*xd1t + tman*xd1t + 2*tman*xd22t + 2*m3**2*xd23t - 2*sqrtS**2*xd23t - 2*uman*xd23t - m3**2*xd2t + 3*tman*xd2t + (m3**2 - 2*sqrtS**2 - uman)*xd3t + 
     &                  m1**2*(xd0t + 2*xd11t + 2*xd12t + 2*xd13t + 3*xd1t + 2*xd23t + xd2t + 2*xd3t) + m2**2*(2*(xd13t + 2*xd23t + xd33t) + 3*xd3t)) + 
     &               g3Rt*(-4*xb01t - 4*xb02t + 16*xc001t + 2*m2**2*xc01t - 2*sqrtS**2*xc01t - 4*M1Lpt**2*xc02t + m3**2*xc02t - tman*xc02t + 4*m2**2*xc03t + m3**2*xc03t - 4*M3Lpt**2*xc03t + 
     &                  2*m4**2*xc03t - sqrtS**2*xc03t - 2*tman*xc03t - uman*xc03t - 4*M1Lpt**2*xc04t + m2**2*xc04t + 2*m4**2*xc04t - 2*sqrtS**2*xc04t - tman*xc04t - uman*xc04t + 
     &                  4*m3**2*xc111t + m2**2*xc11t - 2*m3**2*xc11t - tman*xc11t - uman*xc11t + 4*m2**2*xc121t + 8*m3**2*xc121t - 4*tman*xc121t - 4*uman*xc121t + 2*m2**2*xc13t + 
     &                  m3**2*xc13t - sqrtS**2*xc13t - uman*xc13t + 2*m2**2*xc21t + 4*m4**2*xc21t - 2*sqrtS**2*xc21t - 2*tman*xc21t - 2*uman*xc21t + 4*sqrtS**2*xc221t - m3**2*xc22t + 
     &                  tman*xc22t + 2*m2**2*xc23t - 8*m2**2*xd00t + 8*sqrtS**2*xd00t + 4*M1Lpt**2*m2**2*xd0t + 2*m2**2*M2Lpt**2*xd0t + M1Lpt**2*m3**2*xd0t - m2**2*m3**2*xd0t - 
     &                  4*M1Lpt**2*M3Lpt**2*xd0t + m3**2*M3Lpt**2*xd0t + 2*M1Lpt**2*m4**2*xd0t - M1Lpt**2*sqrtS**2*xd0t - 2*M2Lpt**2*sqrtS**2*xd0t - 2*M1Lpt**2*tman*xd0t - m2**2*tman*xd0t - 
     &                  M3Lpt**2*tman*xd0t + 2*sqrtS**2*tman*xd0t - M1Lpt**2*uman*xd0t + 2*m2**2*m3**2*xd12t - 2*m3**2*sqrtS**2*xd12t - 2*m2**2*tman*xd12t + 2*sqrtS**2*tman*xd12t - 
     &                  2*m2**4*xd13t + 4*m2**2*sqrtS**2*xd13t - 2*sqrtS**4*xd13t + M1Lpt**2*m2**2*xd1t + m2**2*m3**2*xd1t - M1Lpt**2*sqrtS**2*xd1t - m3**2*sqrtS**2*xd1t - m2**2*tman*xd1t + 
     &                  sqrtS**2*tman*xd1t - 2*m2**2*tman*xd22t + 2*sqrtS**2*tman*xd22t - 4*m2**4*xd23t - 2*m2**2*m3**2*xd23t + 6*m2**2*sqrtS**2*xd23t + 2*m3**2*sqrtS**2*xd23t - 
     &                  2*sqrtS**4*xd23t + 2*m2**2*uman*xd23t - 2*sqrtS**2*uman*xd23t + 2*M1Lpt**2*m2**2*xd2t + M1Lpt**2*m3**2*xd2t + m2**2*m3**2*xd2t - m3**2*M3Lpt**2*xd2t - 
     &                  M1Lpt**2*sqrtS**2*xd2t - 3*m2**2*tman*xd2t + M3Lpt**2*tman*xd2t + 2*sqrtS**2*tman*xd2t - M1Lpt**2*uman*xd2t - 2*m2**4*xd33t + 2*m2**2*sqrtS**2*xd33t + 
     &                  (2*M1Lpt**2*m2**2 - (m2 - sqrtS)*(m2 + sqrtS)*(3*m2**2 + m3**2 - M3Lpt**2 - sqrtS**2 - uman))*xd3t - 
     &                  m1**4*(xd0t + 2*xd11t + 2*xd12t + 2*xd13t + 2*xd1t + 2*xd23t + xd2t + xd3t) + 
     &                  m1**2*(2*xc01t + 3*xc02t + xc03t + xc04t + xc11t + 4*xc121t + 2*xc12t + xc13t + 2*xc21t + xc22t - 8*xd00t + 2*M2Lpt**2*xd0t - m3**2*xd0t + 3*M3Lpt**2*xd0t - 
     &                     2*m4**2*xd0t + sqrtS**2*xd0t + uman*xd0t + 2*sqrtS**2*xd11t + 2*m3**2*xd12t + 2*sqrtS**2*xd12t - 2*tman*xd12t + 4*sqrtS**2*xd13t + m3**2*xd1t + 2*M3Lpt**2*xd1t + 
     &                     2*sqrtS**2*xd1t - tman*xd1t - 2*tman*xd22t - 2*m3**2*xd23t + 4*sqrtS**2*xd23t + 2*uman*xd23t - m3**2*xd2t + M3Lpt**2*xd2t + sqrtS**2*xd2t - 2*tman*xd2t + 
     &                     uman*xd2t + M1Lpt**2*(xd0t + xd1t + xd2t) - m3**2*xd3t + M3Lpt**2*xd3t + 2*sqrtS**2*xd3t + uman*xd3t - 
     &                     m2**2*(3*xd0t + 2*xd11t + 2*xd12t + 4*xd13t + 4*xd1t + 6*xd23t + 3*xd2t + 2*xd33t + 6*xd3t)) + 
     &                  M0Lpt**2*(-4*xc01t + 16*xd00t - 4*M2Lpt**2*xd0t + 2*m3**2*xd0t + 2*m4**2*xd0t - 3*sqrtS**2*xd0t + tman*xd0t - uman*xd0t - 4*m3**2*xd12t + 4*tman*xd12t - 
     &                     4*sqrtS**2*xd13t - 2*m3**2*xd1t - sqrtS**2*xd1t + 2*tman*xd1t + 4*tman*xd22t + 4*m3**2*xd23t - 4*sqrtS**2*xd23t - 4*uman*xd23t - sqrtS**2*xd2t + 5*tman*xd2t - 
     &                     uman*xd2t + 2*m3**2*xd3t - 3*sqrtS**2*xd3t - 2*uman*xd3t + m1**2*(2*xd0t + 4*xd11t + 4*xd12t + 4*xd13t + 5*xd1t + 4*xd23t + 2*xd2t + 3*xd3t) + 
     &                     m2**2*(2*xd0t + 4*xd13t + xd1t + 8*xd23t + 2*xd2t + 4*xd33t + 7*xd3t))))))))/(16.*pi**2*propac)  )

      Box_US = Box_US + 2*colorfactorSH*dreal(   (Coupl_Ac*Coupl_Bc*(g0Ru*(2*g1Ru*M1Lpu*(2*g2Lu*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*xd0u + g3Lu*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u)) + 
     &            g2Ru*(g3Lu*M3Lpu*(-2*xc01u - 2*xc04u + 8*xd00u - 2*M0Lpu**2*xd0u - 2*M2Lpu**2*xd0u + uman*(2*xd0u + 2*xd12u + xd1u + 2*(xd22u + xd2u)) - 
     &                  (sqrtS**2 + tman)*(2*xd23u + xd3u) + m4**2*(-2*xd12u - xd1u + 2*xd23u + xd3u) + m1**2*(2*xd11u + 2*xd12u + 2*xd13u + xd1u + 2*xd23u + xd3u) + 
     &                  2*(-(sqrtS**2*xd13u) + m2**2*(xd13u + 2*xd23u + xd33u + xd3u))) + 
     &               g3Ru*M0Lpu*(-2*xc01u - 2*xc02u + 8*xd00u - 2*M2Lpu**2*xd0u + 2*m3**2*xd0u - 2*M3Lpu**2*xd0u + m4**2*xd0u - sqrtS**2*xd0u - tman*xd0u - 2*m4**2*xd12u + 2*uman*xd12u - 
     &                  2*sqrtS**2*xd13u - m4**2*xd1u - sqrtS**2*xd1u + uman*xd1u + 2*uman*xd22u + 2*m4**2*xd23u - 2*sqrtS**2*xd23u - 2*tman*xd23u + m4**2*xd2u - sqrtS**2*xd2u - tman*xd2u + 
     &                  2*uman*xd2u + (m4**2 - sqrtS**2 - tman)*xd3u + m1**2*(xd0u + 2*xd11u + 2*xd12u + 2*xd13u + 2*xd1u + 2*xd23u + xd2u + xd3u) + 
     &                  m2**2*(2*xd0u + 2*xd13u + xd1u + 2*(2*xd23u + xd2u + xd33u + 2*xd3u))))) + 
     &         g1Lu*(4*g2Ru*M2Lpu*(g3Lu*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Ru*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &            g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(-2*xc01u - 2*xc03u + 8*xd00u - 2*M1Lpu**2*xd0u - 2*M2Lpu**2*xd0u + m4**2*xd0u + uman*xd0u - 2*m4**2*xd12u + 2*uman*xd12u - 2*sqrtS**2*xd13u - 
     &                  m4**2*xd1u + uman*xd1u + 2*uman*xd22u + 2*m4**2*xd23u - 2*sqrtS**2*xd23u - 2*tman*xd23u - m4**2*xd2u + 3*uman*xd2u + (m4**2 - 2*sqrtS**2 - tman)*xd3u + 
     &                  m1**2*(xd0u + 2*xd11u + 2*xd12u + 2*xd13u + 3*xd1u + 2*xd23u + xd2u + 2*xd3u) + m2**2*(2*(xd13u + 2*xd23u + xd33u) + 3*xd3u)) + 
     &               g3Lu*(-4*xb01u - 4*xb02u + 16*xc001u + 2*m2**2*xc01u - 2*sqrtS**2*xc01u - 4*M1Lpu**2*xc02u + m4**2*xc02u - uman*xc02u + 4*m2**2*xc03u + 2*m3**2*xc03u - 
     &                  4*M3Lpu**2*xc03u + m4**2*xc03u - sqrtS**2*xc03u - tman*xc03u - 2*uman*xc03u - 4*M1Lpu**2*xc04u + m2**2*xc04u + 2*m3**2*xc04u - 2*sqrtS**2*xc04u - tman*xc04u - 
     &                  uman*xc04u + 4*m4**2*xc111u + m2**2*xc11u - 2*m4**2*xc11u - tman*xc11u - uman*xc11u + 4*m2**2*xc121u + 8*m4**2*xc121u - 4*tman*xc121u - 4*uman*xc121u + 
     &                  2*m2**2*xc13u + m4**2*xc13u - sqrtS**2*xc13u - tman*xc13u + 2*m2**2*xc21u + 4*m3**2*xc21u - 2*sqrtS**2*xc21u - 2*tman*xc21u - 2*uman*xc21u + 4*sqrtS**2*xc221u - 
     &                  m4**2*xc22u + uman*xc22u + 2*m2**2*xc23u - 8*m2**2*xd00u + 8*sqrtS**2*xd00u + 4*M1Lpu**2*m2**2*xd0u + 2*m2**2*M2Lpu**2*xd0u + 2*M1Lpu**2*m3**2*xd0u - 
     &                  4*M1Lpu**2*M3Lpu**2*xd0u + M1Lpu**2*m4**2*xd0u - m2**2*m4**2*xd0u + M3Lpu**2*m4**2*xd0u - M1Lpu**2*sqrtS**2*xd0u - 2*M2Lpu**2*sqrtS**2*xd0u - M1Lpu**2*tman*xd0u - 
     &                  2*M1Lpu**2*uman*xd0u - m2**2*uman*xd0u - M3Lpu**2*uman*xd0u + 2*sqrtS**2*uman*xd0u + 2*m2**2*m4**2*xd12u - 2*m4**2*sqrtS**2*xd12u - 2*m2**2*uman*xd12u + 
     &                  2*sqrtS**2*uman*xd12u - 2*m2**4*xd13u + 4*m2**2*sqrtS**2*xd13u - 2*sqrtS**4*xd13u + M1Lpu**2*m2**2*xd1u + m2**2*m4**2*xd1u - M1Lpu**2*sqrtS**2*xd1u - 
     &                  m4**2*sqrtS**2*xd1u - m2**2*uman*xd1u + sqrtS**2*uman*xd1u - 2*m2**2*uman*xd22u + 2*sqrtS**2*uman*xd22u - 4*m2**4*xd23u - 2*m2**2*m4**2*xd23u + 
     &                  6*m2**2*sqrtS**2*xd23u + 2*m4**2*sqrtS**2*xd23u - 2*sqrtS**4*xd23u + 2*m2**2*tman*xd23u - 2*sqrtS**2*tman*xd23u + 2*M1Lpu**2*m2**2*xd2u + M1Lpu**2*m4**2*xd2u + 
     &                  m2**2*m4**2*xd2u - M3Lpu**2*m4**2*xd2u - M1Lpu**2*sqrtS**2*xd2u - M1Lpu**2*tman*xd2u - 3*m2**2*uman*xd2u + M3Lpu**2*uman*xd2u + 2*sqrtS**2*uman*xd2u - 2*m2**4*xd33u + 
     &                  2*m2**2*sqrtS**2*xd33u + (2*M1Lpu**2*m2**2 - (m2 - sqrtS)*(m2 + sqrtS)*(3*m2**2 - M3Lpu**2 + m4**2 - sqrtS**2 - tman))*xd3u - 
     &                  m1**4*(xd0u + 2*xd11u + 2*xd12u + 2*xd13u + 2*xd1u + 2*xd23u + xd2u + xd3u) + 
     &                  m1**2*(2*xc01u + 3*xc02u + xc03u + xc04u + xc11u + 4*xc121u + 2*xc12u + xc13u + 2*xc21u + xc22u - 8*xd00u + 2*M2Lpu**2*xd0u - 2*m3**2*xd0u + 3*M3Lpu**2*xd0u - 
     &                     m4**2*xd0u + sqrtS**2*xd0u + tman*xd0u + 2*sqrtS**2*xd11u + 2*m4**2*xd12u + 2*sqrtS**2*xd12u - 2*uman*xd12u + 4*sqrtS**2*xd13u + 2*M3Lpu**2*xd1u + m4**2*xd1u + 
     &                     2*sqrtS**2*xd1u - uman*xd1u - 2*uman*xd22u - 2*m4**2*xd23u + 4*sqrtS**2*xd23u + 2*tman*xd23u + M3Lpu**2*xd2u - m4**2*xd2u + sqrtS**2*xd2u + tman*xd2u - 
     &                     2*uman*xd2u + M1Lpu**2*(xd0u + xd1u + xd2u) + (M3Lpu**2 - m4**2 + 2*sqrtS**2 + tman)*xd3u - 
     &                     m2**2*(3*xd0u + 2*xd11u + 2*xd12u + 4*xd13u + 4*xd1u + 6*xd23u + 3*xd2u + 2*xd33u + 6*xd3u)) + 
     &                  M0Lpu**2*(-4*xc01u + 16*xd00u - 4*M2Lpu**2*xd0u + 2*m3**2*xd0u + 2*m4**2*xd0u - 3*sqrtS**2*xd0u - tman*xd0u + uman*xd0u - 4*m4**2*xd12u + 4*uman*xd12u - 
     &                     4*sqrtS**2*xd13u - 2*m4**2*xd1u - sqrtS**2*xd1u + 2*uman*xd1u + 4*uman*xd22u + 4*m4**2*xd23u - 4*sqrtS**2*xd23u - 4*tman*xd23u - sqrtS**2*xd2u - tman*xd2u + 
     &                     5*uman*xd2u + 2*m4**2*xd3u - 3*sqrtS**2*xd3u - 2*tman*xd3u + m1**2*(2*xd0u + 4*xd11u + 4*xd12u + 4*xd13u + 5*xd1u + 4*xd23u + 2*xd2u + 3*xd3u) + 
     &                     m2**2*(2*xd0u + 4*xd13u + xd1u + 8*xd23u + 2*xd2u + 4*xd33u + 7*xd3u)))))) + 
     &      g0Lu*(2*g1Lu*M1Lpu*(2*g2Ru*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*xd0u + g3Ru*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u)) + 
     &            g2Lu*(g3Ru*M3Lpu*(-2*xc01u - 2*xc04u + 8*xd00u - 2*M0Lpu**2*xd0u - 2*M2Lpu**2*xd0u + uman*(2*xd0u + 2*xd12u + xd1u + 2*(xd22u + xd2u)) - 
     &                  (sqrtS**2 + tman)*(2*xd23u + xd3u) + m4**2*(-2*xd12u - xd1u + 2*xd23u + xd3u) + m1**2*(2*xd11u + 2*xd12u + 2*xd13u + xd1u + 2*xd23u + xd3u) + 
     &                  2*(-(sqrtS**2*xd13u) + m2**2*(xd13u + 2*xd23u + xd33u + xd3u))) + 
     &               g3Lu*M0Lpu*(-2*xc01u - 2*xc02u + 8*xd00u - 2*M2Lpu**2*xd0u + 2*m3**2*xd0u - 2*M3Lpu**2*xd0u + m4**2*xd0u - sqrtS**2*xd0u - tman*xd0u - 2*m4**2*xd12u + 2*uman*xd12u - 
     &                  2*sqrtS**2*xd13u - m4**2*xd1u - sqrtS**2*xd1u + uman*xd1u + 2*uman*xd22u + 2*m4**2*xd23u - 2*sqrtS**2*xd23u - 2*tman*xd23u + m4**2*xd2u - sqrtS**2*xd2u - tman*xd2u + 
     &                  2*uman*xd2u + (m4**2 - sqrtS**2 - tman)*xd3u + m1**2*(xd0u + 2*xd11u + 2*xd12u + 2*xd13u + 2*xd1u + 2*xd23u + xd2u + xd3u) + 
     &                  m2**2*(2*xd0u + 2*xd13u + xd1u + 2*(2*xd23u + xd2u + xd33u + 2*xd3u))))) + 
     &         g1Ru*(4*g2Lu*M2Lpu*(g3Ru*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Lu*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &            g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(-2*xc01u - 2*xc03u + 8*xd00u - 2*M1Lpu**2*xd0u - 2*M2Lpu**2*xd0u + m4**2*xd0u + uman*xd0u - 2*m4**2*xd12u + 2*uman*xd12u - 2*sqrtS**2*xd13u - 
     &                  m4**2*xd1u + uman*xd1u + 2*uman*xd22u + 2*m4**2*xd23u - 2*sqrtS**2*xd23u - 2*tman*xd23u - m4**2*xd2u + 3*uman*xd2u + (m4**2 - 2*sqrtS**2 - tman)*xd3u + 
     &                  m1**2*(xd0u + 2*xd11u + 2*xd12u + 2*xd13u + 3*xd1u + 2*xd23u + xd2u + 2*xd3u) + m2**2*(2*(xd13u + 2*xd23u + xd33u) + 3*xd3u)) + 
     &               g3Ru*(-4*xb01u - 4*xb02u + 16*xc001u + 2*m2**2*xc01u - 2*sqrtS**2*xc01u - 4*M1Lpu**2*xc02u + m4**2*xc02u - uman*xc02u + 4*m2**2*xc03u + 2*m3**2*xc03u - 
     &                  4*M3Lpu**2*xc03u + m4**2*xc03u - sqrtS**2*xc03u - tman*xc03u - 2*uman*xc03u - 4*M1Lpu**2*xc04u + m2**2*xc04u + 2*m3**2*xc04u - 2*sqrtS**2*xc04u - tman*xc04u - 
     &                  uman*xc04u + 4*m4**2*xc111u + m2**2*xc11u - 2*m4**2*xc11u - tman*xc11u - uman*xc11u + 4*m2**2*xc121u + 8*m4**2*xc121u - 4*tman*xc121u - 4*uman*xc121u + 
     &                  2*m2**2*xc13u + m4**2*xc13u - sqrtS**2*xc13u - tman*xc13u + 2*m2**2*xc21u + 4*m3**2*xc21u - 2*sqrtS**2*xc21u - 2*tman*xc21u - 2*uman*xc21u + 4*sqrtS**2*xc221u - 
     &                  m4**2*xc22u + uman*xc22u + 2*m2**2*xc23u - 8*m2**2*xd00u + 8*sqrtS**2*xd00u + 4*M1Lpu**2*m2**2*xd0u + 2*m2**2*M2Lpu**2*xd0u + 2*M1Lpu**2*m3**2*xd0u - 
     &                  4*M1Lpu**2*M3Lpu**2*xd0u + M1Lpu**2*m4**2*xd0u - m2**2*m4**2*xd0u + M3Lpu**2*m4**2*xd0u - M1Lpu**2*sqrtS**2*xd0u - 2*M2Lpu**2*sqrtS**2*xd0u - M1Lpu**2*tman*xd0u - 
     &                  2*M1Lpu**2*uman*xd0u - m2**2*uman*xd0u - M3Lpu**2*uman*xd0u + 2*sqrtS**2*uman*xd0u + 2*m2**2*m4**2*xd12u - 2*m4**2*sqrtS**2*xd12u - 2*m2**2*uman*xd12u + 
     &                  2*sqrtS**2*uman*xd12u - 2*m2**4*xd13u + 4*m2**2*sqrtS**2*xd13u - 2*sqrtS**4*xd13u + M1Lpu**2*m2**2*xd1u + m2**2*m4**2*xd1u - M1Lpu**2*sqrtS**2*xd1u - 
     &                  m4**2*sqrtS**2*xd1u - m2**2*uman*xd1u + sqrtS**2*uman*xd1u - 2*m2**2*uman*xd22u + 2*sqrtS**2*uman*xd22u - 4*m2**4*xd23u - 2*m2**2*m4**2*xd23u + 
     &                  6*m2**2*sqrtS**2*xd23u + 2*m4**2*sqrtS**2*xd23u - 2*sqrtS**4*xd23u + 2*m2**2*tman*xd23u - 2*sqrtS**2*tman*xd23u + 2*M1Lpu**2*m2**2*xd2u + M1Lpu**2*m4**2*xd2u + 
     &                  m2**2*m4**2*xd2u - M3Lpu**2*m4**2*xd2u - M1Lpu**2*sqrtS**2*xd2u - M1Lpu**2*tman*xd2u - 3*m2**2*uman*xd2u + M3Lpu**2*uman*xd2u + 2*sqrtS**2*uman*xd2u - 2*m2**4*xd33u + 
     &                  2*m2**2*sqrtS**2*xd33u + (2*M1Lpu**2*m2**2 - (m2 - sqrtS)*(m2 + sqrtS)*(3*m2**2 - M3Lpu**2 + m4**2 - sqrtS**2 - tman))*xd3u - 
     &                  m1**4*(xd0u + 2*xd11u + 2*xd12u + 2*xd13u + 2*xd1u + 2*xd23u + xd2u + xd3u) + 
     &                  m1**2*(2*xc01u + 3*xc02u + xc03u + xc04u + xc11u + 4*xc121u + 2*xc12u + xc13u + 2*xc21u + xc22u - 8*xd00u + 2*M2Lpu**2*xd0u - 2*m3**2*xd0u + 3*M3Lpu**2*xd0u - 
     &                     m4**2*xd0u + sqrtS**2*xd0u + tman*xd0u + 2*sqrtS**2*xd11u + 2*m4**2*xd12u + 2*sqrtS**2*xd12u - 2*uman*xd12u + 4*sqrtS**2*xd13u + 2*M3Lpu**2*xd1u + m4**2*xd1u + 
     &                     2*sqrtS**2*xd1u - uman*xd1u - 2*uman*xd22u - 2*m4**2*xd23u + 4*sqrtS**2*xd23u + 2*tman*xd23u + M3Lpu**2*xd2u - m4**2*xd2u + sqrtS**2*xd2u + tman*xd2u - 
     &                     2*uman*xd2u + M1Lpu**2*(xd0u + xd1u + xd2u) + (M3Lpu**2 - m4**2 + 2*sqrtS**2 + tman)*xd3u - 
     &                     m2**2*(3*xd0u + 2*xd11u + 2*xd12u + 4*xd13u + 4*xd1u + 6*xd23u + 3*xd2u + 2*xd33u + 6*xd3u)) + 
     &                  M0Lpu**2*(-4*xc01u + 16*xd00u - 4*M2Lpu**2*xd0u + 2*m3**2*xd0u + 2*m4**2*xd0u - 3*sqrtS**2*xd0u - tman*xd0u + uman*xd0u - 4*m4**2*xd12u + 4*uman*xd12u - 
     &                     4*sqrtS**2*xd13u - 2*m4**2*xd1u - sqrtS**2*xd1u + 2*uman*xd1u + 4*uman*xd22u + 4*m4**2*xd23u - 4*sqrtS**2*xd23u - 4*tman*xd23u - sqrtS**2*xd2u - tman*xd2u + 
     &                     5*uman*xd2u + 2*m4**2*xd3u - 3*sqrtS**2*xd3u - 2*tman*xd3u + m1**2*(2*xd0u + 4*xd11u + 4*xd12u + 4*xd13u + 5*xd1u + 4*xd23u + 2*xd2u + 3*xd3u) + 
     &                     m2**2*(2*xd0u + 4*xd13u + xd1u + 8*xd23u + 2*xd2u + 4*xd33u + 7*xd3u))))))))/(16.*pi**2*propac)  )
  

      enddo


!==================================================================================================================================  
!Box(t)*Vector

      do j=sumstartV, sumendV

      Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
      Coupl_Bc = dconjg(VVV(j,xtype,ytype))
      propac = dconjg(DenSmV(j))
        
      Box_TV = Box_TV + 2*colorfactorSV*dreal(  (Coupl_Ac*Coupl_Bc*(g0Rt*(2*g1Rt*M1Lpt*(g2Lt*M2Lpt*(12*g3Rt*M0Lpt*M3Lpt*(-tman + uman)*xd0t + 
     &               g3Lt*(-2*m1**4*xd1t + 20*m1**2*m2**2*xd1t - 2*m2**4*xd1t + 4*m1**2*m3**2*xd1t + 4*m2**2*m3**2*xd1t + 4*m1**2*m4**2*xd1t + 4*m2**2*m4**2*xd1t + 6*m1**2*sqrtS**2*xd1t + 
     &                  6*m2**2*sqrtS**2*xd1t - 2*m3**2*sqrtS**2*xd1t - 2*m4**2*sqrtS**2*xd1t - 4*sqrtS**4*xd1t - 2*m1**4*xd2t + 6*m1**2*m2**2*xd2t - 4*m2**4*xd2t - 2*m1**2*m3**2*xd2t - 
     &                  10*m2**2*m3**2*xd2t + 4*m1**2*m4**2*xd2t + 8*m2**2*m4**2*xd2t + 5*m1**2*sqrtS**2*xd2t + 7*m2**2*sqrtS**2*xd2t + 3*m3**2*sqrtS**2*xd2t - 3*m4**2*sqrtS**2*xd2t - 
     &                  3*sqrtS**4*xd2t - tman**2*xd2t + uman**2*xd2t + uman*
     &                   (6*xc02t + 6*xc04t + 6*M0Lpt**2*xd0t - 6*m2**2*xd0t + 6*M3Lpt**2*xd0t - 4*m1**2*xd1t - 4*m2**2*xd1t + 2*sqrtS**2*xd1t + m1**2*xd2t - 3*m2**2*xd2t - m3**2*xd2t - 
     &                     m4**2*xd2t - 2*sqrtS**2*xd2t) + tman*(-6*xc02t - 6*xc04t - 6*M0Lpt**2*xd0t + 6*m2**2*xd0t - 6*M3Lpt**2*xd0t - 4*m1**2*xd1t - 4*m2**2*xd1t + 2*sqrtS**2*xd1t - 
     &                     3*m1**2*xd2t + 5*m2**2*xd2t + m3**2*xd2t + m4**2*xd2t + 2*sqrtS**2*xd2t))) + 
     &            g2Rt*(g3Rt*M0Lpt*(4*m1**2*m3**2*xd0t + 8*m2**2*m3**2*xd0t + 4*m2**2*m4**2*xd0t - m1**2*sqrtS**2*xd0t - 3*m2**2*sqrtS**2*xd0t - 3*m3**2*sqrtS**2*xd0t - m4**2*sqrtS**2*xd0t + 
     &                  sqrtS**4*xd0t - 6*m1**2*m3**2*xd11t - 2*m2**2*m3**2*xd11t + 6*m1**2*m4**2*xd11t + 2*m2**2*m4**2*xd11t + 2*m3**2*sqrtS**2*xd11t - 2*m4**2*sqrtS**2*xd11t + 
     &                  8*m1**2*m3**2*xd12t + 16*m1**2*m4**2*xd12t + 8*m2**2*m4**2*xd12t - 6*m1**2*sqrtS**2*xd12t - 2*m2**2*sqrtS**2*xd12t - 2*m3**2*sqrtS**2*xd12t - 6*m4**2*sqrtS**2*xd12t + 
     &                  2*sqrtS**4*xd12t + 4*m1**2*m3**2*xd13t - 4*m2**2*m3**2*xd13t - 4*m1**2*m4**2*xd13t + 4*m2**2*m4**2*xd13t + 6*m1**2*m3**2*xd1t - 2*m2**2*m3**2*xd1t + 
     &                  6*m1**2*m4**2*xd1t + 6*m2**2*m4**2*xd1t - 3*m1**2*sqrtS**2*xd1t - m2**2*sqrtS**2*xd1t - m3**2*sqrtS**2*xd1t - 3*m4**2*sqrtS**2*xd1t + sqrtS**4*xd1t + 
     &                  8*m1**2*m3**2*xd22t + 8*m2**2*m3**2*xd22t + 8*m1**2*m4**2*xd22t + 8*m2**2*m4**2*xd22t - 4*m1**2*sqrtS**2*xd22t - 4*m2**2*sqrtS**2*xd22t - 4*m3**2*sqrtS**2*xd22t - 
     &                  4*m4**2*sqrtS**2*xd22t + 2*sqrtS**4*xd22t + 8*m1**2*m3**2*xd23t + 16*m2**2*m3**2*xd23t + 8*m2**2*m4**2*xd23t - 2*m1**2*sqrtS**2*xd23t - 6*m2**2*sqrtS**2*xd23t - 
     &                  6*m3**2*sqrtS**2*xd23t - 2*m4**2*sqrtS**2*xd23t + 2*sqrtS**4*xd23t + 12*m1**2*m3**2*xd2t + 16*m2**2*m3**2*xd2t + 8*m1**2*m4**2*xd2t + 12*m2**2*m4**2*xd2t - 
     &                  5*m1**2*sqrtS**2*xd2t - 7*m2**2*sqrtS**2*xd2t - 7*m3**2*sqrtS**2*xd2t - 5*m4**2*sqrtS**2*xd2t + 3*sqrtS**4*xd2t + 2*m1**2*m3**2*xd33t + 6*m2**2*m3**2*xd33t - 
     &                  2*m1**2*m4**2*xd33t - 6*m2**2*m4**2*xd33t - 2*m3**2*sqrtS**2*xd33t + 2*m4**2*sqrtS**2*xd33t + uman**2*(-3*xd0t + 2*xd12t + xd1t + 2*xd22t - 6*xd23t - xd2t - 3*xd3t) + 
     &                  6*m1**2*m3**2*xd3t + 14*m2**2*m3**2*xd3t - 2*m1**2*m4**2*xd3t - 2*m2**2*m4**2*xd3t - m1**2*sqrtS**2*xd3t - 3*m2**2*sqrtS**2*xd3t - 5*m3**2*sqrtS**2*xd3t + 
     &                  m4**2*sqrtS**2*xd3t + sqrtS**4*xd3t + tman**2*(xd0t - 6*xd12t - 3*xd1t - 6*xd22t + 2*xd23t - 5*xd2t + xd3t) - 
     &                  uman*(6*xc01t + 6*xc02t - 24*xd00t - m1**2*xd0t - 3*m2**2*xd0t + 6*M2Lpt**2*xd0t - 3*m3**2*xd0t + 6*M3Lpt**2*xd0t - 7*m4**2*xd0t + 2*sqrtS**2*xd0t - 2*m1**2*xd11t + 
     &                     2*m2**2*xd11t - 2*sqrtS**2*xd11t + 2*m1**2*xd12t + 6*m2**2*xd12t + 6*m3**2*xd12t + 2*m4**2*xd12t - 4*sqrtS**2*xd12t + 4*sqrtS**2*xd13t + m1**2*xd1t + 
     &                     3*m2**2*xd1t + 3*m3**2*xd1t + m4**2*xd1t + 4*m1**2*xd22t + 4*m2**2*xd22t + 4*m3**2*xd22t + 4*m4**2*xd22t - 4*sqrtS**2*xd22t - 2*m1**2*xd23t - 6*m2**2*xd23t - 
     &                     6*m3**2*xd23t - 2*m4**2*xd23t + 4*sqrtS**2*xd23t + 3*m1**2*xd2t + m2**2*xd2t + m3**2*xd2t + 3*m4**2*xd2t - 2*sqrtS**2*xd2t + 2*m1**2*xd33t - 2*m2**2*xd33t - 
     &                     2*sqrtS**2*xd33t + m1**2*xd3t - 5*m2**2*xd3t - 3*m3**2*xd3t - m4**2*xd3t) + 
     &                  tman*(6*xc01t + 6*xc02t - 24*xd00t - m1**2*xd0t - 3*m2**2*xd0t + 6*M2Lpt**2*xd0t - 3*m3**2*xd0t + 6*M3Lpt**2*xd0t - 7*m4**2*xd0t + 2*sqrtS**2*xd0t + 2*uman*xd0t - 
     &                     2*m1**2*xd11t + 2*m2**2*xd11t - 2*sqrtS**2*xd11t + 2*m1**2*xd12t + 6*m2**2*xd12t + 6*m3**2*xd12t + 2*m4**2*xd12t - 4*sqrtS**2*xd12t + 4*uman*xd12t + 
     &                     4*sqrtS**2*xd13t + m1**2*xd1t + 3*m2**2*xd1t + 3*m3**2*xd1t + m4**2*xd1t + 2*uman*xd1t + 4*m1**2*xd22t + 4*m2**2*xd22t + 4*m3**2*xd22t + 4*m4**2*xd22t - 
     &                     4*sqrtS**2*xd22t + 4*uman*xd22t - 2*m1**2*xd23t - 6*m2**2*xd23t - 6*m3**2*xd23t - 2*m4**2*xd23t + 4*sqrtS**2*xd23t + 4*uman*xd23t + 3*m1**2*xd2t + m2**2*xd2t + 
     &                     m3**2*xd2t + 3*m4**2*xd2t - 2*sqrtS**2*xd2t + 6*uman*xd2t + 2*m1**2*xd33t - 2*m2**2*xd33t - 2*sqrtS**2*xd33t + m1**2*xd3t - 5*m2**2*xd3t - 3*m3**2*xd3t - 
     &                     m4**2*xd3t + 2*uman*xd3t)) + g3Lt*M3Lpt*(-6*m1**2*m3**2*xd11t - 2*m2**2*m3**2*xd11t + 6*m1**2*m4**2*xd11t + 2*m2**2*m4**2*xd11t + 2*m3**2*sqrtS**2*xd11t - 
     &                  2*m4**2*sqrtS**2*xd11t + 8*m1**2*m3**2*xd12t + 16*m1**2*m4**2*xd12t + 8*m2**2*m4**2*xd12t - 6*m1**2*sqrtS**2*xd12t - 2*m2**2*sqrtS**2*xd12t - 2*m3**2*sqrtS**2*xd12t - 
     &                  6*m4**2*sqrtS**2*xd12t + 2*sqrtS**4*xd12t + 4*m1**2*m3**2*xd13t - 4*m2**2*m3**2*xd13t - 4*m1**2*m4**2*xd13t + 4*m2**2*m4**2*xd13t + 4*m1**2*m3**2*xd1t + 
     &                  8*m1**2*m4**2*xd1t + 4*m2**2*m4**2*xd1t - 3*m1**2*sqrtS**2*xd1t - m2**2*sqrtS**2*xd1t - m3**2*sqrtS**2*xd1t - 3*m4**2*sqrtS**2*xd1t + sqrtS**4*xd1t + 
     &                  8*m1**2*m3**2*xd22t + 8*m2**2*m3**2*xd22t + 8*m1**2*m4**2*xd22t + 8*m2**2*m4**2*xd22t - 4*m1**2*sqrtS**2*xd22t - 4*m2**2*sqrtS**2*xd22t - 4*m3**2*sqrtS**2*xd22t - 
     &                  4*m4**2*sqrtS**2*xd22t + 2*sqrtS**4*xd22t + 8*m1**2*m3**2*xd23t + 16*m2**2*m3**2*xd23t + 8*m2**2*m4**2*xd23t - 2*m1**2*sqrtS**2*xd23t - 6*m2**2*sqrtS**2*xd23t - 
     &                  6*m3**2*sqrtS**2*xd23t - 2*m4**2*sqrtS**2*xd23t + 2*sqrtS**4*xd23t + 8*m1**2*m3**2*xd2t + 8*m2**2*m3**2*xd2t + 8*m1**2*m4**2*xd2t + 8*m2**2*m4**2*xd2t - 
     &                  4*m1**2*sqrtS**2*xd2t - 4*m2**2*sqrtS**2*xd2t - 4*m3**2*sqrtS**2*xd2t - 4*m4**2*sqrtS**2*xd2t + 2*sqrtS**4*xd2t + 2*m1**2*m3**2*xd33t + 6*m2**2*m3**2*xd33t - 
     &                  2*m1**2*m4**2*xd33t - 6*m2**2*m4**2*xd33t - 2*m3**2*sqrtS**2*xd33t + 2*m4**2*sqrtS**2*xd33t + uman**2*(2*xd12t + xd1t + 2*xd22t - 6*xd23t + 2*xd2t - 3*xd3t) + 
     &                  4*m1**2*m3**2*xd3t + 8*m2**2*m3**2*xd3t + 4*m2**2*m4**2*xd3t - m1**2*sqrtS**2*xd3t - 3*m2**2*sqrtS**2*xd3t - 3*m3**2*sqrtS**2*xd3t - m4**2*sqrtS**2*xd3t + 
     &                  sqrtS**4*xd3t + tman**2*(-6*xd0t - 6*xd12t - 3*xd1t - 6*xd22t + 2*xd23t - 6*xd2t + xd3t) - 
     &                  uman*(6*xc01t + 6*xc04t - 24*xd00t + 6*M0Lpt**2*xd0t + 6*M2Lpt**2*xd0t - 2*m1**2*xd11t + 2*m2**2*xd11t - 2*sqrtS**2*xd11t + 2*m1**2*xd12t + 6*m2**2*xd12t + 
     &                     6*m3**2*xd12t + 2*m4**2*xd12t - 4*sqrtS**2*xd12t + 4*sqrtS**2*xd13t + m1**2*xd1t + 3*m2**2*xd1t + 3*m3**2*xd1t + m4**2*xd1t - 2*sqrtS**2*xd1t + 4*m1**2*xd22t + 
     &                     4*m2**2*xd22t + 4*m3**2*xd22t + 4*m4**2*xd22t - 4*sqrtS**2*xd22t - 2*m1**2*xd23t - 6*m2**2*xd23t - 6*m3**2*xd23t - 2*m4**2*xd23t + 4*sqrtS**2*xd23t + 
     &                     4*m1**2*xd2t + 4*m2**2*xd2t + 4*m3**2*xd2t + 4*m4**2*xd2t - 4*sqrtS**2*xd2t + 2*m1**2*xd33t - 2*m2**2*xd33t - 2*sqrtS**2*xd33t - m1**2*xd3t - 3*m2**2*xd3t - 
     &                     3*m3**2*xd3t - m4**2*xd3t + 2*sqrtS**2*xd3t) + tman*
     &                   (6*xc01t + 6*xc04t - 24*xd00t + 6*M0Lpt**2*xd0t + 6*M2Lpt**2*xd0t + 6*uman*xd0t - 2*m1**2*xd11t + 2*m2**2*xd11t - 2*sqrtS**2*xd11t + 2*m1**2*xd12t + 6*m2**2*xd12t + 
     &                     6*m3**2*xd12t + 2*m4**2*xd12t - 4*sqrtS**2*xd12t + 4*uman*xd12t + 4*sqrtS**2*xd13t + m1**2*xd1t + 3*m2**2*xd1t + 3*m3**2*xd1t + m4**2*xd1t - 2*sqrtS**2*xd1t + 
     &                     2*uman*xd1t + 4*m1**2*xd22t + 4*m2**2*xd22t + 4*m3**2*xd22t + 4*m4**2*xd22t - 4*sqrtS**2*xd22t + 4*uman*xd22t - 2*m1**2*xd23t - 6*m2**2*xd23t - 6*m3**2*xd23t - 
     &                     2*m4**2*xd23t + 4*sqrtS**2*xd23t + 4*uman*xd23t + 4*m1**2*xd2t + 4*m2**2*xd2t + 4*m3**2*xd2t + 4*m4**2*xd2t - 4*sqrtS**2*xd2t + 4*uman*xd2t + 2*m1**2*xd33t - 
     &                     2*m2**2*xd33t - 2*sqrtS**2*xd33t - m1**2*xd3t - 3*m2**2*xd3t - 3*m3**2*xd3t - m4**2*xd3t + 2*sqrtS**2*xd3t + 2*uman*xd3t)))) + 
     &         g1Lt*(-2*g2Rt*M2Lpt*(g3Lt*M3Lpt*(-2*m1**4*xd2t - 14*m1**2*m2**2*xd2t - 14*m1**2*m3**2*xd2t - 6*m2**2*m3**2*xd2t + 4*m1**2*m4**2*xd2t + m1**2*sqrtS**2*xd2t - 
     &                  m2**2*sqrtS**2*xd2t + 5*m3**2*sqrtS**2*xd2t - m4**2*sqrtS**2*xd2t + sqrtS**4*xd2t + tman**2*xd2t - uman**2*xd2t + 2*m1**4*xd3t - 20*m1**2*m2**2*xd3t + 2*m2**4*xd3t - 
     &                  4*m1**2*m3**2*xd3t - 4*m2**2*m3**2*xd3t - 4*m1**2*m4**2*xd3t - 4*m2**2*m4**2*xd3t - 6*m1**2*sqrtS**2*xd3t - 6*m2**2*sqrtS**2*xd3t + 2*m3**2*sqrtS**2*xd3t + 
     &                  2*m4**2*sqrtS**2*xd3t + 4*sqrtS**4*xd3t + uman*(-6*xc03t - 6*xc04t - 6*M0Lpt**2*xd0t + 6*m1**2*xd0t - 6*M1Lpt**2*xd0t + 9*m1**2*xd2t + m2**2*xd2t + m3**2*xd2t + 
     &                     m4**2*xd2t + 4*m1**2*xd3t + 4*m2**2*xd3t - 2*sqrtS**2*xd3t) + 
     &                  tman*(6*xc03t + 6*xc04t + 6*M0Lpt**2*xd0t - 6*m1**2*xd0t + 6*M1Lpt**2*xd0t + m1**2*xd2t + 5*m2**2*xd2t - m3**2*xd2t - m4**2*xd2t - 4*sqrtS**2*xd2t + 4*m1**2*xd3t + 
     &                     4*m2**2*xd3t - 2*sqrtS**2*xd3t)) + 2*g3Rt*M0Lpt*(m1**4*xd0t - 10*m1**2*m2**2*xd0t + m2**4*xd0t - 2*m1**2*m3**2*xd0t - 2*m2**2*m3**2*xd0t - 2*m1**2*m4**2*xd0t - 
     &                  2*m2**2*m4**2*xd0t - 3*m1**2*sqrtS**2*xd0t - 3*m2**2*sqrtS**2*xd0t + m3**2*sqrtS**2*xd0t + m4**2*sqrtS**2*xd0t + 2*sqrtS**4*xd0t + m1**4*xd1t - 10*m1**2*m2**2*xd1t + 
     &                  m2**4*xd1t - 2*m1**2*m3**2*xd1t - 2*m2**2*m3**2*xd1t - 2*m1**2*m4**2*xd1t - 2*m2**2*m4**2*xd1t - 3*m1**2*sqrtS**2*xd1t - 3*m2**2*sqrtS**2*xd1t + m3**2*sqrtS**2*xd1t + 
     &                  m4**2*sqrtS**2*xd1t + 2*sqrtS**4*xd1t - 10*m1**2*m2**2*xd2t + 2*m2**4*xd2t - 6*m1**2*m3**2*xd2t + 2*m2**2*m3**2*xd2t - 4*m2**2*m4**2*xd2t - 2*m1**2*sqrtS**2*xd2t - 
     &                  4*m2**2*sqrtS**2*xd2t + m3**2*sqrtS**2*xd2t + m4**2*sqrtS**2*xd2t + 2*sqrtS**4*xd2t + tman**2*xd2t - uman**2*xd2t + m1**4*xd3t - 10*m1**2*m2**2*xd3t + m2**4*xd3t - 
     &                  2*m1**2*m3**2*xd3t - 2*m2**2*m3**2*xd3t - 2*m1**2*m4**2*xd3t - 2*m2**2*m4**2*xd3t - 3*m1**2*sqrtS**2*xd3t - 3*m2**2*sqrtS**2*xd3t + m3**2*sqrtS**2*xd3t + 
     &                  m4**2*sqrtS**2*xd3t + 2*sqrtS**4*xd3t + tman*(3*xc02t + 3*xc03t + 2*m1**2*xd0t + 3*M1Lpt**2*xd0t + 2*m2**2*xd0t + 3*M3Lpt**2*xd0t - 4*sqrtS**2*xd0t + 2*m1**2*xd1t + 
     &                     2*m2**2*xd1t - sqrtS**2*xd1t + 2*m1**2*xd2t - m3**2*xd2t - m4**2*xd2t - 3*sqrtS**2*xd2t + 2*m1**2*xd3t + 2*m2**2*xd3t - sqrtS**2*xd3t) + 
     &                  uman*(-3*xc02t - 3*xc03t + 2*m1**2*xd0t - 3*M1Lpt**2*xd0t + 2*m2**2*xd0t - 3*M3Lpt**2*xd0t + 2*sqrtS**2*xd0t + 2*m1**2*xd1t + 2*m2**2*xd1t - sqrtS**2*xd1t + 
     &                     4*m1**2*xd2t + 2*m2**2*xd2t + m3**2*xd2t + m4**2*xd2t + sqrtS**2*xd2t + 2*m1**2*xd3t + 2*m2**2*xd3t - sqrtS**2*xd3t))) + 
     &            g2Lt*(2*g3Rt*M0Lpt*M3Lpt*(4*m1**2*m3**2*xd0t + 8*m1**2*m4**2*xd0t + 4*m2**2*m4**2*xd0t - 3*m1**2*sqrtS**2*xd0t - m2**2*sqrtS**2*xd0t - m3**2*sqrtS**2*xd0t - 
     &                  3*m4**2*sqrtS**2*xd0t + sqrtS**4*xd0t - 6*m1**2*m3**2*xd11t - 2*m2**2*m3**2*xd11t + 6*m1**2*m4**2*xd11t + 2*m2**2*m4**2*xd11t + 2*m3**2*sqrtS**2*xd11t - 
     &                  2*m4**2*sqrtS**2*xd11t + 8*m1**2*m3**2*xd12t + 16*m1**2*m4**2*xd12t + 8*m2**2*m4**2*xd12t - 6*m1**2*sqrtS**2*xd12t - 2*m2**2*sqrtS**2*xd12t - 2*m3**2*sqrtS**2*xd12t - 
     &                  6*m4**2*sqrtS**2*xd12t + 2*sqrtS**4*xd12t + 4*m1**2*m3**2*xd13t - 4*m2**2*m3**2*xd13t - 4*m1**2*m4**2*xd13t + 4*m2**2*m4**2*xd13t - 2*m1**2*m3**2*xd1t - 
     &                  2*m2**2*m3**2*xd1t + 14*m1**2*m4**2*xd1t + 6*m2**2*m4**2*xd1t - 3*m1**2*sqrtS**2*xd1t - m2**2*sqrtS**2*xd1t + m3**2*sqrtS**2*xd1t - 5*m4**2*sqrtS**2*xd1t + 
     &                  sqrtS**4*xd1t + 8*m1**2*m3**2*xd22t + 8*m2**2*m3**2*xd22t + 8*m1**2*m4**2*xd22t + 8*m2**2*m4**2*xd22t - 4*m1**2*sqrtS**2*xd22t - 4*m2**2*sqrtS**2*xd22t - 
     &                  4*m3**2*sqrtS**2*xd22t - 4*m4**2*sqrtS**2*xd22t + 2*sqrtS**4*xd22t + 8*m1**2*m3**2*xd23t + 16*m2**2*m3**2*xd23t + 8*m2**2*m4**2*xd23t - 2*m1**2*sqrtS**2*xd23t - 
     &                  6*m2**2*sqrtS**2*xd23t - 6*m3**2*sqrtS**2*xd23t - 2*m4**2*sqrtS**2*xd23t + 2*sqrtS**4*xd23t + 12*m1**2*m3**2*xd2t + 8*m2**2*m3**2*xd2t + 16*m1**2*m4**2*xd2t + 
     &                  12*m2**2*m4**2*xd2t - 7*m1**2*sqrtS**2*xd2t - 5*m2**2*sqrtS**2*xd2t - 5*m3**2*sqrtS**2*xd2t - 7*m4**2*sqrtS**2*xd2t + 3*sqrtS**4*xd2t + 2*m1**2*m3**2*xd33t + 
     &                  6*m2**2*m3**2*xd33t - 2*m1**2*m4**2*xd33t - 6*m2**2*m4**2*xd33t - 2*m3**2*sqrtS**2*xd33t + 2*m4**2*sqrtS**2*xd33t + 
     &                  uman**2*(xd0t + 2*xd12t + xd1t + 2*xd22t - 6*xd23t + 3*xd2t - 3*xd3t) + 6*m1**2*m3**2*xd3t + 6*m2**2*m3**2*xd3t - 2*m1**2*m4**2*xd3t + 6*m2**2*m4**2*xd3t - 
     &                  m1**2*sqrtS**2*xd3t - 3*m2**2*sqrtS**2*xd3t - 3*m3**2*sqrtS**2*xd3t - m4**2*sqrtS**2*xd3t + sqrtS**4*xd3t + 
     &                  tman**2*(-3*xd0t - 6*xd12t - 3*xd1t - 6*xd22t + 2*xd23t - 9*xd2t + xd3t) - 
     &                  uman*(6*xc01t + 6*xc03t - 24*xd00t + m1**2*xd0t + 6*M1Lpt**2*xd0t + 3*m2**2*xd0t + 6*M2Lpt**2*xd0t - 3*m3**2*xd0t + m4**2*xd0t - 2*sqrtS**2*xd0t - 2*m1**2*xd11t + 
     &                     2*m2**2*xd11t - 2*sqrtS**2*xd11t + 2*m1**2*xd12t + 6*m2**2*xd12t + 6*m3**2*xd12t + 2*m4**2*xd12t - 4*sqrtS**2*xd12t + 4*sqrtS**2*xd13t - m1**2*xd1t + 
     &                     5*m2**2*xd1t + 3*m3**2*xd1t + m4**2*xd1t - 4*sqrtS**2*xd1t + 4*m1**2*xd22t + 4*m2**2*xd22t + 4*m3**2*xd22t + 4*m4**2*xd22t - 4*sqrtS**2*xd22t - 2*m1**2*xd23t - 
     &                     6*m2**2*xd23t - 6*m3**2*xd23t - 2*m4**2*xd23t + 4*sqrtS**2*xd23t + 5*m1**2*xd2t + 7*m2**2*xd2t + 7*m3**2*xd2t + 5*m4**2*xd2t - 6*sqrtS**2*xd2t + 2*m1**2*xd33t - 
     &                     2*m2**2*xd33t - 2*sqrtS**2*xd33t - m1**2*xd3t - 3*m2**2*xd3t - 3*m3**2*xd3t - m4**2*xd3t + 4*sqrtS**2*xd3t) + 
     &                  tman*(6*xc01t + 6*xc03t - 24*xd00t + m1**2*xd0t + 6*M1Lpt**2*xd0t + 3*m2**2*xd0t + 6*M2Lpt**2*xd0t - 3*m3**2*xd0t + m4**2*xd0t - 2*sqrtS**2*xd0t + 2*uman*xd0t - 
     &                     2*m1**2*xd11t + 2*m2**2*xd11t - 2*sqrtS**2*xd11t + 2*m1**2*xd12t + 6*m2**2*xd12t + 6*m3**2*xd12t + 2*m4**2*xd12t - 4*sqrtS**2*xd12t + 4*uman*xd12t + 
     &                     4*sqrtS**2*xd13t - m1**2*xd1t + 5*m2**2*xd1t + 3*m3**2*xd1t + m4**2*xd1t - 4*sqrtS**2*xd1t + 2*uman*xd1t + 4*m1**2*xd22t + 4*m2**2*xd22t + 4*m3**2*xd22t + 
     &                     4*m4**2*xd22t - 4*sqrtS**2*xd22t + 4*uman*xd22t - 2*m1**2*xd23t - 6*m2**2*xd23t - 6*m3**2*xd23t - 2*m4**2*xd23t + 4*sqrtS**2*xd23t + 4*uman*xd23t + 5*m1**2*xd2t + 
     &                     7*m2**2*xd2t + 7*m3**2*xd2t + 5*m4**2*xd2t - 6*sqrtS**2*xd2t + 6*uman*xd2t + 2*m1**2*xd33t - 2*m2**2*xd33t - 2*sqrtS**2*xd33t - m1**2*xd3t - 3*m2**2*xd3t - 
     &                     3*m3**2*xd3t - m4**2*xd3t + 4*sqrtS**2*xd3t + 2*uman*xd3t)) + 
     &               g3Lt*(-2*m1**4*xc01t + 20*m1**2*m2**2*xc01t - 2*m2**4*xc01t + 4*m1**2*m3**2*xc01t + 4*m2**2*m3**2*xc01t + 4*m1**2*m4**2*xc01t + 4*m2**2*m4**2*xc01t + 
     &                  6*m1**2*sqrtS**2*xc01t + 6*m2**2*sqrtS**2*xc01t - 2*m3**2*sqrtS**2*xc01t - 2*m4**2*sqrtS**2*xc01t - 4*sqrtS**4*xc01t + 4*m1**2*m3**2*xc02t + 8*m1**2*m4**2*xc02t + 
     &                  4*m2**2*m4**2*xc02t - 3*m1**2*sqrtS**2*xc02t - m2**2*sqrtS**2*xc02t - m3**2*sqrtS**2*xc02t - 3*m4**2*sqrtS**2*xc02t + sqrtS**4*xc02t + 4*m1**2*m3**2*xc03t + 
     &                  8*m2**2*m3**2*xc03t + 4*m2**2*m4**2*xc03t - m1**2*sqrtS**2*xc03t - 3*m2**2*sqrtS**2*xc03t - 3*m3**2*sqrtS**2*xc03t - m4**2*sqrtS**2*xc03t + sqrtS**4*xc03t - 
     &                  8*m1**2*m3**2*xc04t + 8*m2**2*m3**2*xc04t + 2*m1**2*sqrtS**2*xc04t - 2*m2**2*sqrtS**2*xc04t - 12*m1**2*m3**2*xc111t + 12*m2**2*m3**2*xc111t - 4*m1**2*m4**2*xc111t + 
     &                  4*m2**2*m4**2*xc111t + 4*m1**2*sqrtS**2*xc111t - 4*m2**2*sqrtS**2*xc111t + 2*m1**4*xc11t - 2*m2**4*xc11t - 12*m1**2*m3**2*xc11t + 12*m2**2*m3**2*xc11t - 
     &                  8*m1**2*m4**2*xc11t + 8*m2**2*m4**2*xc11t + 4*m1**2*sqrtS**2*xc11t - 4*m2**2*sqrtS**2*xc11t - 32*m1**2*m3**2*xc121t + 32*m2**2*m3**2*xc121t + 
     &                  8*m1**2*sqrtS**2*xc121t - 8*m2**2*sqrtS**2*xc121t - 6*m1**2*m3**2*xc12t - 2*m2**2*m3**2*xc12t + 6*m1**2*m4**2*xc12t + 2*m2**2*m4**2*xc12t + 2*m3**2*sqrtS**2*xc12t - 
     &                  2*m4**2*sqrtS**2*xc12t + 2*m1**4*xc13t - 6*m1**2*m2**2*xc13t + 4*m2**4*xc13t + 6*m1**2*m3**2*xc13t + 18*m2**2*m3**2*xc13t - 4*m1**2*m4**2*xc13t - 
     &                  4*m2**2*m4**2*xc13t - 6*m1**2*sqrtS**2*xc13t - 10*m2**2*sqrtS**2*xc13t - 6*m3**2*sqrtS**2*xc13t + 2*m4**2*sqrtS**2*xc13t + 4*sqrtS**4*xc13t - 24*m1**2*m3**2*xc21t + 
     &                  24*m2**2*m3**2*xc21t + 8*m1**2*m4**2*xc21t - 8*m2**2*m4**2*xc21t + 4*m1**2*sqrtS**2*xc21t - 4*m2**2*sqrtS**2*xc21t - 16*m1**2*m3**2*xc221t + 16*m2**2*m3**2*xc221t + 
     &                  16*m1**2*m4**2*xc221t - 16*m2**2*m4**2*xc221t - 2*m1**4*xc22t - 14*m1**2*m2**2*xc22t - 10*m1**2*m3**2*xc22t - 6*m2**2*m3**2*xc22t + 12*m1**2*m4**2*xc22t + 
     &                  4*m2**2*m4**2*xc22t - 2*m1**2*sqrtS**2*xc22t - 2*m2**2*sqrtS**2*xc22t + 4*m3**2*sqrtS**2*xc22t - 4*m4**2*sqrtS**2*xc22t + 2*sqrtS**4*xc22t + 2*m1**2*m3**2*xc23t + 
     &                  6*m2**2*m3**2*xc23t - 2*m1**2*m4**2*xc23t - 6*m2**2*m4**2*xc23t - 2*m3**2*sqrtS**2*xc23t + 2*m4**2*sqrtS**2*xc23t - 2*M0Lpt**2*m1**4*xd0t + 
     &                  20*M0Lpt**2*m1**2*m2**2*xd0t - 2*M0Lpt**2*m2**4*xd0t - 2*m1**4*M2Lpt**2*xd0t + 20*m1**2*m2**2*M2Lpt**2*xd0t - 2*m2**4*M2Lpt**2*xd0t + 12*M0Lpt**2*m1**2*m3**2*xd0t - 
     &                  4*m1**4*m3**2*xd0t + 4*m1**2*M1Lpt**2*m3**2*xd0t + 12*M0Lpt**2*m2**2*m3**2*xd0t - 12*m1**2*m2**2*m3**2*xd0t + 8*M1Lpt**2*m2**2*m3**2*xd0t + 
     &                  4*m1**2*M2Lpt**2*m3**2*xd0t + 4*m2**2*M2Lpt**2*m3**2*xd0t + 4*m1**2*m3**2*M3Lpt**2*xd0t + 12*M0Lpt**2*m1**2*m4**2*xd0t + 12*M0Lpt**2*m2**2*m4**2*xd0t - 
     &                  12*m1**2*m2**2*m4**2*xd0t + 4*M1Lpt**2*m2**2*m4**2*xd0t - 4*m2**4*m4**2*xd0t + 4*m1**2*M2Lpt**2*m4**2*xd0t + 4*m2**2*M2Lpt**2*m4**2*xd0t + 
     &                  8*m1**2*M3Lpt**2*m4**2*xd0t + 4*m2**2*M3Lpt**2*m4**2*xd0t + 2*M0Lpt**2*m1**2*sqrtS**2*xd0t + m1**4*sqrtS**2*xd0t - m1**2*M1Lpt**2*sqrtS**2*xd0t + 
     &                  2*M0Lpt**2*m2**2*sqrtS**2*xd0t + 6*m1**2*m2**2*sqrtS**2*xd0t - 3*M1Lpt**2*m2**2*sqrtS**2*xd0t + m2**4*sqrtS**2*xd0t + 6*m1**2*M2Lpt**2*sqrtS**2*xd0t + 
     &                  6*m2**2*M2Lpt**2*sqrtS**2*xd0t - 6*M0Lpt**2*m3**2*sqrtS**2*xd0t + 3*m1**2*m3**2*sqrtS**2*xd0t - 3*M1Lpt**2*m3**2*sqrtS**2*xd0t + m2**2*m3**2*sqrtS**2*xd0t - 
     &                  2*M2Lpt**2*m3**2*sqrtS**2*xd0t - 3*m1**2*M3Lpt**2*sqrtS**2*xd0t - m2**2*M3Lpt**2*sqrtS**2*xd0t - m3**2*M3Lpt**2*sqrtS**2*xd0t - 6*M0Lpt**2*m4**2*sqrtS**2*xd0t + 
     &                  m1**2*m4**2*sqrtS**2*xd0t - M1Lpt**2*m4**2*sqrtS**2*xd0t + 3*m2**2*m4**2*sqrtS**2*xd0t - 2*M2Lpt**2*m4**2*sqrtS**2*xd0t - 3*M3Lpt**2*m4**2*sqrtS**2*xd0t - 
     &                  2*M0Lpt**2*sqrtS**4*xd0t - m1**2*sqrtS**4*xd0t + M1Lpt**2*sqrtS**4*xd0t - m2**2*sqrtS**4*xd0t - 4*M2Lpt**2*sqrtS**4*xd0t + M3Lpt**2*sqrtS**4*xd0t - 
     &                  12*M0Lpt**2*m1**2*m3**2*xd11t + 6*m1**4*m3**2*xd11t - 4*M0Lpt**2*m2**2*m3**2*xd11t + 8*m1**2*m2**2*m3**2*xd11t + 2*m2**4*m3**2*xd11t + 12*M0Lpt**2*m1**2*m4**2*xd11t - 
     &                  6*m1**4*m4**2*xd11t + 4*M0Lpt**2*m2**2*m4**2*xd11t - 8*m1**2*m2**2*m4**2*xd11t - 2*m2**4*m4**2*xd11t + 4*M0Lpt**2*m3**2*sqrtS**2*xd11t - 
     &                  8*m1**2*m3**2*sqrtS**2*xd11t - 4*m2**2*m3**2*sqrtS**2*xd11t - 4*M0Lpt**2*m4**2*sqrtS**2*xd11t + 8*m1**2*m4**2*sqrtS**2*xd11t + 4*m2**2*m4**2*sqrtS**2*xd11t + 
     &                  2*m3**2*sqrtS**4*xd11t - 2*m4**2*sqrtS**4*xd11t + 16*M0Lpt**2*m1**2*m3**2*xd12t - 8*m1**4*m3**2*xd12t - 8*m1**2*m2**2*m3**2*xd12t + 32*M0Lpt**2*m1**2*m4**2*xd12t - 
     &                  16*m1**4*m4**2*xd12t + 16*M0Lpt**2*m2**2*m4**2*xd12t - 24*m1**2*m2**2*m4**2*xd12t - 8*m2**4*m4**2*xd12t - 12*M0Lpt**2*m1**2*sqrtS**2*xd12t + 6*m1**4*sqrtS**2*xd12t - 
     &                  4*M0Lpt**2*m2**2*sqrtS**2*xd12t + 8*m1**2*m2**2*sqrtS**2*xd12t + 2*m2**4*sqrtS**2*xd12t - 4*M0Lpt**2*m3**2*sqrtS**2*xd12t + 10*m1**2*m3**2*sqrtS**2*xd12t + 
     &                  2*m2**2*m3**2*sqrtS**2*xd12t - 12*M0Lpt**2*m4**2*sqrtS**2*xd12t + 22*m1**2*m4**2*sqrtS**2*xd12t + 14*m2**2*m4**2*sqrtS**2*xd12t + 4*M0Lpt**2*sqrtS**4*xd12t - 
     &                  8*m1**2*sqrtS**4*xd12t - 4*m2**2*sqrtS**4*xd12t - 2*m3**2*sqrtS**4*xd12t - 6*m4**2*sqrtS**4*xd12t + 2*sqrtS**6*xd12t + 8*M0Lpt**2*m1**2*m3**2*xd13t - 
     &                  4*m1**4*m3**2*xd13t - 8*M0Lpt**2*m2**2*m3**2*xd13t + 4*m2**4*m3**2*xd13t - 8*M0Lpt**2*m1**2*m4**2*xd13t + 4*m1**4*m4**2*xd13t + 8*M0Lpt**2*m2**2*m4**2*xd13t - 
     &                  4*m2**4*m4**2*xd13t + 4*m1**2*m3**2*sqrtS**2*xd13t - 4*m2**2*m3**2*sqrtS**2*xd13t - 4*m1**2*m4**2*sqrtS**2*xd13t + 4*m2**2*m4**2*sqrtS**2*xd13t - 
     &                  2*M0Lpt**2*m1**4*xd1t + 2*m1**4*M1Lpt**2*xd1t + 20*M0Lpt**2*m1**2*m2**2*xd1t - 20*m1**2*M1Lpt**2*m2**2*xd1t - 2*M0Lpt**2*m2**4*xd1t + 2*M1Lpt**2*m2**4*xd1t + 
     &                  8*M0Lpt**2*m1**2*m3**2*xd1t - 8*m1**4*m3**2*xd1t - 2*m1**2*M1Lpt**2*m3**2*xd1t + 24*m1**2*m2**2*m3**2*xd1t - 6*M1Lpt**2*m2**2*m3**2*xd1t + 4*m1**2*m3**4*xd1t + 
     &                  4*m2**2*m3**4*xd1t - 6*m1**2*m3**2*M3Lpt**2*xd1t - 2*m2**2*m3**2*M3Lpt**2*xd1t + 24*M0Lpt**2*m1**2*m4**2*xd1t - 6*m1**4*m4**2*xd1t - 6*m1**2*M1Lpt**2*m4**2*xd1t + 
     &                  16*M0Lpt**2*m2**2*m4**2*xd1t - 20*m1**2*m2**2*m4**2*xd1t - 2*M1Lpt**2*m2**2*m4**2*xd1t - 6*m2**4*m4**2*xd1t + 4*m1**2*m3**2*m4**2*xd1t + 4*m2**2*m3**2*m4**2*xd1t + 
     &                  6*m1**2*M3Lpt**2*m4**2*xd1t + 2*m2**2*M3Lpt**2*m4**2*xd1t + 3*m1**4*sqrtS**2*xd1t - 6*m1**2*M1Lpt**2*sqrtS**2*xd1t + 4*M0Lpt**2*m2**2*sqrtS**2*xd1t + 
     &                  4*m1**2*m2**2*sqrtS**2*xd1t - 6*M1Lpt**2*m2**2*sqrtS**2*xd1t + m2**4*sqrtS**2*xd1t - 2*M0Lpt**2*m3**2*sqrtS**2*xd1t + 11*m1**2*m3**2*sqrtS**2*xd1t + 
     &                  2*M1Lpt**2*m3**2*sqrtS**2*xd1t + 5*m2**2*m3**2*sqrtS**2*xd1t - 2*m3**4*sqrtS**2*xd1t + 2*m3**2*M3Lpt**2*sqrtS**2*xd1t - 10*M0Lpt**2*m4**2*sqrtS**2*xd1t + 
     &                  11*m1**2*m4**2*sqrtS**2*xd1t + 2*M1Lpt**2*m4**2*sqrtS**2*xd1t + 9*m2**2*m4**2*sqrtS**2*xd1t - 2*m3**2*m4**2*sqrtS**2*xd1t - 2*M3Lpt**2*m4**2*sqrtS**2*xd1t - 
     &                  2*M0Lpt**2*sqrtS**4*xd1t - 4*m1**2*sqrtS**4*xd1t + 4*M1Lpt**2*sqrtS**4*xd1t - 2*m2**2*sqrtS**4*xd1t - 5*m3**2*sqrtS**4*xd1t - 3*m4**2*sqrtS**4*xd1t + sqrtS**6*xd1t + 
     &                  16*M0Lpt**2*m1**2*m3**2*xd22t - 8*m1**4*m3**2*xd22t + 16*M0Lpt**2*m2**2*m3**2*xd22t - 16*m1**2*m2**2*m3**2*xd22t - 8*m2**4*m3**2*xd22t + 
     &                  16*M0Lpt**2*m1**2*m4**2*xd22t - 8*m1**4*m4**2*xd22t + 16*M0Lpt**2*m2**2*m4**2*xd22t - 16*m1**2*m2**2*m4**2*xd22t - 8*m2**4*m4**2*xd22t - 
     &                  8*M0Lpt**2*m1**2*sqrtS**2*xd22t + 4*m1**4*sqrtS**2*xd22t - 8*M0Lpt**2*m2**2*sqrtS**2*xd22t + 8*m1**2*m2**2*sqrtS**2*xd22t + 4*m2**4*sqrtS**2*xd22t - 
     &                  8*M0Lpt**2*m3**2*sqrtS**2*xd22t + 12*m1**2*m3**2*sqrtS**2*xd22t + 12*m2**2*m3**2*sqrtS**2*xd22t - 8*M0Lpt**2*m4**2*sqrtS**2*xd22t + 12*m1**2*m4**2*sqrtS**2*xd22t + 
     &                  12*m2**2*m4**2*sqrtS**2*xd22t + 4*M0Lpt**2*sqrtS**4*xd22t - 6*m1**2*sqrtS**4*xd22t - 6*m2**2*sqrtS**4*xd22t - 4*m3**2*sqrtS**4*xd22t - 4*m4**2*sqrtS**4*xd22t + 
     &                  2*sqrtS**6*xd22t + 16*M0Lpt**2*m1**2*m3**2*xd23t - 8*m1**4*m3**2*xd23t + 32*M0Lpt**2*m2**2*m3**2*xd23t - 24*m1**2*m2**2*m3**2*xd23t - 16*m2**4*m3**2*xd23t + 
     &                  16*M0Lpt**2*m2**2*m4**2*xd23t - 8*m1**2*m2**2*m4**2*xd23t - 8*m2**4*m4**2*xd23t - 4*M0Lpt**2*m1**2*sqrtS**2*xd23t + 2*m1**4*sqrtS**2*xd23t - 
     &                  12*M0Lpt**2*m2**2*sqrtS**2*xd23t + 8*m1**2*m2**2*sqrtS**2*xd23t + 6*m2**4*sqrtS**2*xd23t - 12*M0Lpt**2*m3**2*sqrtS**2*xd23t + 14*m1**2*m3**2*sqrtS**2*xd23t + 
     &                  22*m2**2*m3**2*sqrtS**2*xd23t - 4*M0Lpt**2*m4**2*sqrtS**2*xd23t + 2*m1**2*m4**2*sqrtS**2*xd23t + 10*m2**2*m4**2*sqrtS**2*xd23t + 4*M0Lpt**2*sqrtS**4*xd23t - 
     &                  4*m1**2*sqrtS**4*xd23t - 8*m2**2*sqrtS**4*xd23t - 6*m3**2*sqrtS**4*xd23t - 2*m4**2*sqrtS**4*xd23t + 2*sqrtS**6*xd23t + 2*m1**4*M1Lpt**2*xd2t + 
     &                  20*M0Lpt**2*m1**2*m2**2*xd2t - 6*m1**2*M1Lpt**2*m2**2*xd2t - 4*M0Lpt**2*m2**4*xd2t + 4*M1Lpt**2*m2**4*xd2t + 36*M0Lpt**2*m1**2*m3**2*xd2t - 14*m1**4*m3**2*xd2t + 
     &                  6*m1**2*M1Lpt**2*m3**2*xd2t + 20*M0Lpt**2*m2**2*m3**2*xd2t - 22*m1**2*m2**2*m3**2*xd2t + 18*M1Lpt**2*m2**2*m3**2*xd2t - 12*m2**4*m3**2*xd2t - 2*m1**2*m3**4*xd2t - 
     &                  10*m2**2*m3**4*xd2t - 2*m1**4*M3Lpt**2*xd2t - 14*m1**2*m2**2*M3Lpt**2*xd2t - 10*m1**2*m3**2*M3Lpt**2*xd2t - 6*m2**2*m3**2*M3Lpt**2*xd2t + 
     &                  24*M0Lpt**2*m1**2*m4**2*xd2t - 6*m1**4*m4**2*xd2t - 4*m1**2*M1Lpt**2*m4**2*xd2t + 32*M0Lpt**2*m2**2*m4**2*xd2t - 14*m1**2*m2**2*m4**2*xd2t - 
     &                  4*M1Lpt**2*m2**2*m4**2*xd2t - 12*m2**4*m4**2*xd2t + 18*m1**2*m3**2*m4**2*xd2t + 14*m2**2*m3**2*m4**2*xd2t + 12*m1**2*M3Lpt**2*m4**2*xd2t + 
     &                  4*m2**2*M3Lpt**2*m4**2*xd2t - 4*m1**2*m4**4*xd2t - 8*M0Lpt**2*m1**2*sqrtS**2*xd2t + 5*m1**4*sqrtS**2*xd2t - 6*m1**2*M1Lpt**2*sqrtS**2*xd2t - 
     &                  4*M0Lpt**2*m2**2*sqrtS**2*xd2t + 14*m1**2*m2**2*sqrtS**2*xd2t - 10*M1Lpt**2*m2**2*sqrtS**2*xd2t + 5*m2**4*sqrtS**2*xd2t - 14*M0Lpt**2*m3**2*sqrtS**2*xd2t + 
     &                  20*m1**2*m3**2*sqrtS**2*xd2t - 6*M1Lpt**2*m3**2*sqrtS**2*xd2t + 20*m2**2*m3**2*sqrtS**2*xd2t + 3*m3**4*sqrtS**2*xd2t - 2*m1**2*M3Lpt**2*sqrtS**2*xd2t - 
     &                  2*m2**2*M3Lpt**2*sqrtS**2*xd2t + 4*m3**2*M3Lpt**2*sqrtS**2*xd2t - 14*M0Lpt**2*m4**2*sqrtS**2*xd2t + 12*m1**2*m4**2*sqrtS**2*xd2t + 2*M1Lpt**2*m4**2*sqrtS**2*xd2t + 
     &                  16*m2**2*m4**2*sqrtS**2*xd2t - 8*m3**2*m4**2*sqrtS**2*xd2t - 4*M3Lpt**2*m4**2*sqrtS**2*xd2t + m4**4*sqrtS**2*xd2t + 2*M0Lpt**2*sqrtS**4*xd2t - 7*m1**2*sqrtS**4*xd2t + 
     &                  4*M1Lpt**2*sqrtS**4*xd2t - 7*m2**2*sqrtS**4*xd2t - 7*m3**2*sqrtS**4*xd2t + 2*M3Lpt**2*sqrtS**4*xd2t - 5*m4**2*sqrtS**4*xd2t + 2*sqrtS**6*xd2t + 2*tman**3*xd2t + 
     &                  4*M0Lpt**2*m1**2*m3**2*xd33t - 2*m1**4*m3**2*xd33t + 12*M0Lpt**2*m2**2*m3**2*xd33t - 8*m1**2*m2**2*m3**2*xd33t - 6*m2**4*m3**2*xd33t - 4*M0Lpt**2*m1**2*m4**2*xd33t + 
     &                  2*m1**4*m4**2*xd33t - 12*M0Lpt**2*m2**2*m4**2*xd33t + 8*m1**2*m2**2*m4**2*xd33t + 6*m2**4*m4**2*xd33t - 4*M0Lpt**2*m3**2*sqrtS**2*xd33t + 
     &                  4*m1**2*m3**2*sqrtS**2*xd33t + 8*m2**2*m3**2*sqrtS**2*xd33t + 4*M0Lpt**2*m4**2*sqrtS**2*xd33t - 4*m1**2*m4**2*sqrtS**2*xd33t - 8*m2**2*m4**2*sqrtS**2*xd33t - 
     &                  2*m3**2*sqrtS**4*xd33t + 2*m4**2*sqrtS**4*xd33t - 2*M0Lpt**2*m1**4*xd3t + 20*M0Lpt**2*m1**2*m2**2*xd3t - 2*M0Lpt**2*m2**4*xd3t + 16*M0Lpt**2*m1**2*m3**2*xd3t - 
     &                  6*m1**4*m3**2*xd3t + 2*m1**2*M1Lpt**2*m3**2*xd3t + 24*M0Lpt**2*m2**2*m3**2*xd3t - 20*m1**2*m2**2*m3**2*xd3t + 6*M1Lpt**2*m2**2*m3**2*xd3t - 6*m2**4*m3**2*xd3t + 
     &                  2*m1**4*M3Lpt**2*xd3t - 20*m1**2*m2**2*M3Lpt**2*xd3t + 2*m2**4*M3Lpt**2*xd3t - 2*m1**2*m3**2*M3Lpt**2*xd3t - 6*m2**2*m3**2*M3Lpt**2*xd3t - 
     &                  2*m1**2*M1Lpt**2*m4**2*xd3t + 8*M0Lpt**2*m2**2*m4**2*xd3t + 24*m1**2*m2**2*m4**2*xd3t - 6*M1Lpt**2*m2**2*m4**2*xd3t - 8*m2**4*m4**2*xd3t + 4*m1**2*m3**2*m4**2*xd3t + 
     &                  4*m2**2*m3**2*m4**2*xd3t - 6*m1**2*M3Lpt**2*m4**2*xd3t - 2*m2**2*M3Lpt**2*m4**2*xd3t + 4*m1**2*m4**4*xd3t + 4*m2**2*m4**4*xd3t + 4*M0Lpt**2*m1**2*sqrtS**2*xd3t + 
     &                  m1**4*sqrtS**2*xd3t + 4*m1**2*m2**2*sqrtS**2*xd3t + 3*m2**4*sqrtS**2*xd3t - 10*M0Lpt**2*m3**2*sqrtS**2*xd3t + 9*m1**2*m3**2*sqrtS**2*xd3t - 
     &                  2*M1Lpt**2*m3**2*sqrtS**2*xd3t + 11*m2**2*m3**2*sqrtS**2*xd3t - 6*m1**2*M3Lpt**2*sqrtS**2*xd3t - 6*m2**2*M3Lpt**2*sqrtS**2*xd3t + 2*m3**2*M3Lpt**2*sqrtS**2*xd3t - 
     &                  2*M0Lpt**2*m4**2*sqrtS**2*xd3t + 5*m1**2*m4**2*sqrtS**2*xd3t + 2*M1Lpt**2*m4**2*sqrtS**2*xd3t + 11*m2**2*m4**2*sqrtS**2*xd3t - 2*m3**2*m4**2*sqrtS**2*xd3t + 
     &                  2*M3Lpt**2*m4**2*sqrtS**2*xd3t - 2*m4**4*sqrtS**2*xd3t - 2*M0Lpt**2*sqrtS**4*xd3t - 2*m1**2*sqrtS**4*xd3t - 4*m2**2*sqrtS**4*xd3t - 3*m3**2*sqrtS**4*xd3t + 
     &                  4*M3Lpt**2*sqrtS**4*xd3t - 5*m4**2*sqrtS**4*xd3t + sqrtS**6*xd3t + 
     &                  uman**2*(xc02t - 3*xc03t - 4*xc04t - 2*xc11t - 16*xc121t - 4*xc13t - 8*xc21t - 2*M0Lpt**2*xd0t + 3*m1**2*xd0t - 3*M1Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 
     &                     4*M0Lpt**2*xd12t - 2*m1**2*xd12t - 2*m2**2*xd12t + 2*sqrtS**2*xd12t + 2*M0Lpt**2*xd1t - m1**2*xd1t - m2**2*xd1t + sqrtS**2*xd1t + 4*M0Lpt**2*xd22t - 
     &                     2*m1**2*xd22t - 2*m2**2*xd22t + 2*sqrtS**2*xd22t - 12*M0Lpt**2*xd23t + 6*m1**2*xd23t + 6*m2**2*xd23t - 6*sqrtS**2*xd23t + 4*M0Lpt**2*xd2t + m1**2*xd2t - 
     &                     4*M1Lpt**2*xd2t - 3*m2**2*xd2t + m3**2*xd2t + m4**2*xd2t + 2*sqrtS**2*xd2t - 6*M0Lpt**2*xd3t + 3*m1**2*xd3t + 3*m2**2*xd3t - 3*sqrtS**2*xd3t) + 
     &                  tman**2*(3*xc02t + 7*xc03t + 4*xc04t + 2*xc11t + 16*xc121t + 2*xc13t + 8*xc21t - 2*xc22t - 2*M0Lpt**2*xd0t + 3*m1**2*xd0t + 7*M1Lpt**2*xd0t + 7*m2**2*xd0t + 
     &                     3*M3Lpt**2*xd0t - 8*sqrtS**2*xd0t - 12*M0Lpt**2*xd12t + 6*m1**2*xd12t + 6*m2**2*xd12t - 6*sqrtS**2*xd12t - 6*M0Lpt**2*xd1t + 7*m1**2*xd1t + 7*m2**2*xd1t - 
     &                     5*sqrtS**2*xd1t - 12*M0Lpt**2*xd22t + 6*m1**2*xd22t + 6*m2**2*xd22t - 6*sqrtS**2*xd22t + 4*M0Lpt**2*xd23t - 2*m1**2*xd23t - 2*m2**2*xd23t + 2*sqrtS**2*xd23t - 
     &                     16*M0Lpt**2*xd2t + 9*m1**2*xd2t + 2*M1Lpt**2*xd2t + 9*m2**2*xd2t - 3*m3**2*xd2t - 2*M3Lpt**2*xd2t - 3*m4**2*xd2t - 12*sqrtS**2*xd2t + 2*M0Lpt**2*xd3t + 
     &                     3*m1**2*xd3t + 3*m2**2*xd3t - sqrtS**2*xd3t) + uman*
     &                   (-12*xb01t - 12*xb02t + 48*xc001t - 12*M0Lpt**2*xc01t + 2*m1**2*xc01t + 2*m2**2*xc01t - 4*sqrtS**2*xc01t + 5*m1**2*xc02t - 12*M1Lpt**2*xc02t - 3*m2**2*xc02t + 
     &                     3*m3**2*xc02t - m4**2*xc02t + 2*sqrtS**2*xc02t + m1**2*xc03t + 9*m2**2*xc03t + 3*m3**2*xc03t - 12*M3Lpt**2*xc03t + 7*m4**2*xc03t - 2*sqrtS**2*xc03t + 
     &                     4*m1**2*xc04t - 12*M1Lpt**2*xc04t + 4*m2**2*xc04t + 8*m4**2*xc04t - 6*sqrtS**2*xc04t + 4*m3**2*xc111t - 4*m4**2*xc111t + 4*sqrtS**2*xc111t + 4*m2**2*xc11t - 
     &                     16*m3**2*xc11t - 4*m4**2*xc11t + 16*m1**2*xc121t + 16*m2**2*xc121t + 24*m3**2*xc121t + 8*m4**2*xc121t + 2*m1**2*xc12t - 2*m2**2*xc12t + 2*sqrtS**2*xc12t + 
     &                     6*m2**2*xc13t + 4*m3**2*xc13t + 2*m4**2*xc13t + 8*m1**2*xc21t + 8*m2**2*xc21t + 16*m4**2*xc21t - 4*sqrtS**2*xc21t + 16*sqrtS**2*xc221t + 8*m1**2*xc22t - 
     &                     2*m2**2*xc22t - 2*m3**2*xc22t + 2*sqrtS**2*xc22t - 2*m1**2*xc23t + 2*m2**2*xc23t + 2*sqrtS**2*xc23t + 48*M0Lpt**2*xd00t - 24*m1**2*xd00t - 24*m2**2*xd00t + 
     &                     24*sqrtS**2*xd00t - 4*M0Lpt**2*m1**2*xd0t - m1**4*xd0t + m1**2*M1Lpt**2*xd0t - 4*M0Lpt**2*m2**2*xd0t - 2*m1**2*m2**2*xd0t + 9*M1Lpt**2*m2**2*xd0t + 3*m2**4*xd0t - 
     &                     12*M0Lpt**2*M2Lpt**2*xd0t + 2*m1**2*M2Lpt**2*xd0t + 2*m2**2*M2Lpt**2*xd0t + 6*M0Lpt**2*m3**2*xd0t - 3*m1**2*m3**2*xd0t + 3*M1Lpt**2*m3**2*xd0t - 
     &                     3*m2**2*m3**2*xd0t + 5*m1**2*M3Lpt**2*xd0t - 12*M1Lpt**2*M3Lpt**2*xd0t - 3*m2**2*M3Lpt**2*xd0t + 3*m3**2*M3Lpt**2*xd0t + 6*M0Lpt**2*m4**2*xd0t - 
     &                     7*m1**2*m4**2*xd0t + 7*M1Lpt**2*m4**2*xd0t + m2**2*m4**2*xd0t - M3Lpt**2*m4**2*xd0t - 4*M0Lpt**2*sqrtS**2*xd0t + 2*m1**2*sqrtS**2*xd0t - 2*M1Lpt**2*sqrtS**2*xd0t - 
     &                     2*m2**2*sqrtS**2*xd0t - 4*M2Lpt**2*sqrtS**2*xd0t + 2*M3Lpt**2*sqrtS**2*xd0t + 4*M0Lpt**2*m1**2*xd11t - 2*m1**4*xd11t - 4*M0Lpt**2*m2**2*xd11t + 2*m2**4*xd11t + 
     &                     4*M0Lpt**2*sqrtS**2*xd11t - 4*m2**2*sqrtS**2*xd11t + 2*sqrtS**4*xd11t - 4*M0Lpt**2*m1**2*xd12t + 2*m1**4*xd12t - 12*M0Lpt**2*m2**2*xd12t + 8*m1**2*m2**2*xd12t + 
     &                     6*m2**4*xd12t - 12*M0Lpt**2*m3**2*xd12t + 6*m1**2*m3**2*xd12t + 6*m2**2*m3**2*xd12t - 4*M0Lpt**2*m4**2*xd12t + 2*m1**2*m4**2*xd12t + 2*m2**2*m4**2*xd12t + 
     &                     8*M0Lpt**2*sqrtS**2*xd12t - 6*m1**2*sqrtS**2*xd12t - 10*m2**2*sqrtS**2*xd12t - 6*m3**2*sqrtS**2*xd12t - 2*m4**2*sqrtS**2*xd12t + 4*sqrtS**4*xd12t - 
     &                     8*M0Lpt**2*sqrtS**2*xd13t + 4*m1**2*sqrtS**2*xd13t + 4*m2**2*sqrtS**2*xd13t - 4*sqrtS**4*xd13t - 4*M0Lpt**2*m1**2*xd1t + m1**4*xd1t + 4*m1**2*M1Lpt**2*xd1t - 
     &                     12*M0Lpt**2*m2**2*xd1t + 2*m1**2*m2**2*xd1t + 4*M1Lpt**2*m2**2*xd1t + 5*m2**4*xd1t - 6*M0Lpt**2*m3**2*xd1t - m1**2*m3**2*xd1t - m2**2*m3**2*xd1t + 
     &                     2*m1**2*M3Lpt**2*xd1t - 2*m2**2*M3Lpt**2*xd1t - 2*M0Lpt**2*m4**2*xd1t + m1**2*m4**2*xd1t + m2**2*m4**2*xd1t + 6*M0Lpt**2*sqrtS**2*xd1t - m1**2*sqrtS**2*xd1t - 
     &                     4*M1Lpt**2*sqrtS**2*xd1t - 7*m2**2*sqrtS**2*xd1t - m3**2*sqrtS**2*xd1t + 2*M3Lpt**2*sqrtS**2*xd1t - m4**2*sqrtS**2*xd1t + 2*sqrtS**4*xd1t - 
     &                     8*M0Lpt**2*m1**2*xd22t + 4*m1**4*xd22t - 8*M0Lpt**2*m2**2*xd22t + 8*m1**2*m2**2*xd22t + 4*m2**4*xd22t - 8*M0Lpt**2*m3**2*xd22t + 4*m1**2*m3**2*xd22t + 
     &                     4*m2**2*m3**2*xd22t - 8*M0Lpt**2*m4**2*xd22t + 4*m1**2*m4**2*xd22t + 4*m2**2*m4**2*xd22t + 8*M0Lpt**2*sqrtS**2*xd22t - 8*m1**2*sqrtS**2*xd22t - 
     &                     8*m2**2*sqrtS**2*xd22t - 4*m3**2*sqrtS**2*xd22t - 4*m4**2*sqrtS**2*xd22t + 4*sqrtS**4*xd22t + 4*M0Lpt**2*m1**2*xd23t - 2*m1**4*xd23t + 12*M0Lpt**2*m2**2*xd23t - 
     &                     8*m1**2*m2**2*xd23t - 6*m2**4*xd23t + 12*M0Lpt**2*m3**2*xd23t - 6*m1**2*m3**2*xd23t - 6*m2**2*m3**2*xd23t + 4*M0Lpt**2*m4**2*xd23t - 2*m1**2*m4**2*xd23t - 
     &                     2*m2**2*m4**2*xd23t - 8*M0Lpt**2*sqrtS**2*xd23t + 6*m1**2*sqrtS**2*xd23t + 10*m2**2*sqrtS**2*xd23t + 6*m3**2*sqrtS**2*xd23t + 2*m4**2*sqrtS**2*xd23t - 
     &                     4*sqrtS**4*xd23t - 16*M0Lpt**2*m1**2*xd2t + 3*m1**4*xd2t - 12*M0Lpt**2*m2**2*xd2t + 6*m1**2*m2**2*xd2t + 6*M1Lpt**2*m2**2*xd2t + 7*m2**4*xd2t - 
     &                     10*M0Lpt**2*m3**2*xd2t + 2*m1**2*m3**2*xd2t + 4*M1Lpt**2*m3**2*xd2t + 4*m2**2*m3**2*xd2t - m3**4*xd2t + 8*m1**2*M3Lpt**2*xd2t - 2*m2**2*M3Lpt**2*xd2t - 
     &                     2*m3**2*M3Lpt**2*xd2t - 10*M0Lpt**2*m4**2*xd2t - 6*m1**2*m4**2*xd2t + 2*M1Lpt**2*m4**2*xd2t + 4*m2**2*m4**2*xd2t - 2*m3**2*m4**2*xd2t - m4**4*xd2t + 
     &                     6*M0Lpt**2*sqrtS**2*xd2t - 6*m1**2*sqrtS**2*xd2t - 10*m2**2*sqrtS**2*xd2t - 6*m3**2*sqrtS**2*xd2t + 2*M3Lpt**2*sqrtS**2*xd2t - 4*m4**2*sqrtS**2*xd2t + 
     &                     4*sqrtS**4*xd2t - 4*M0Lpt**2*m1**2*xd33t + 2*m1**4*xd33t + 4*M0Lpt**2*m2**2*xd33t - 2*m2**4*xd33t + 4*M0Lpt**2*sqrtS**2*xd33t - 4*m1**2*sqrtS**2*xd33t + 
     &                     2*sqrtS**4*xd33t - 4*M0Lpt**2*m1**2*xd3t + m1**4*xd3t - 2*m1**2*M1Lpt**2*xd3t + 4*M0Lpt**2*m2**2*xd3t - 6*m1**2*m2**2*xd3t + 2*M1Lpt**2*m2**2*xd3t - 3*m2**4*xd3t + 
     &                     6*M0Lpt**2*m3**2*xd3t - 3*m1**2*m3**2*xd3t - 3*m2**2*m3**2*xd3t + 4*m1**2*M3Lpt**2*xd3t + 4*m2**2*M3Lpt**2*xd3t + 2*M0Lpt**2*m4**2*xd3t - 5*m1**2*m4**2*xd3t - 
     &                     5*m2**2*m4**2*xd3t - 2*M0Lpt**2*sqrtS**2*xd3t + m1**2*sqrtS**2*xd3t + 2*M1Lpt**2*sqrtS**2*xd3t + 7*m2**2*sqrtS**2*xd3t + 3*m3**2*sqrtS**2*xd3t - 
     &                     4*M3Lpt**2*sqrtS**2*xd3t + 3*m4**2*sqrtS**2*xd3t - 2*sqrtS**4*xd3t) + 
     &                  tman*(12*xb01t + 12*xb02t - 48*xc001t + 12*M0Lpt**2*xc01t - 10*m1**2*xc01t - 10*m2**2*xc01t + 8*sqrtS**2*xc01t - 5*m1**2*xc02t + 12*M1Lpt**2*xc02t + 3*m2**2*xc02t - 
     &                     3*m3**2*xc02t + m4**2*xc02t - 2*sqrtS**2*xc02t - 4*uman*xc02t - m1**2*xc03t - 9*m2**2*xc03t - 3*m3**2*xc03t + 12*M3Lpt**2*xc03t - 7*m4**2*xc03t + 
     &                     2*sqrtS**2*xc03t - 4*uman*xc03t - 4*m1**2*xc04t + 12*M1Lpt**2*xc04t - 4*m2**2*xc04t - 8*m4**2*xc04t + 6*sqrtS**2*xc04t - 4*m3**2*xc111t + 4*m4**2*xc111t - 
     &                     4*sqrtS**2*xc111t - 4*m1**2*xc11t + 16*m3**2*xc11t + 4*m4**2*xc11t - 16*m1**2*xc121t - 16*m2**2*xc121t - 24*m3**2*xc121t - 8*m4**2*xc121t - 2*m1**2*xc12t + 
     &                     2*m2**2*xc12t - 2*sqrtS**2*xc12t + 2*m1**2*xc13t - 8*m2**2*xc13t - 4*m3**2*xc13t - 2*m4**2*xc13t + 2*uman*xc13t - 8*m1**2*xc21t - 8*m2**2*xc21t - 16*m4**2*xc21t + 
     &                     4*sqrtS**2*xc21t - 16*sqrtS**2*xc221t + 2*m1**2*xc22t + 8*m2**2*xc22t + 2*m3**2*xc22t - 6*sqrtS**2*xc22t + 2*uman*xc22t + 2*m1**2*xc23t - 2*m2**2*xc23t - 
     &                     2*sqrtS**2*xc23t - 48*M0Lpt**2*xd00t + 24*m1**2*xd00t + 24*m2**2*xd00t - 24*sqrtS**2*xd00t - 4*M0Lpt**2*m1**2*xd0t + 3*m1**4*xd0t - m1**2*M1Lpt**2*xd0t - 
     &                     4*M0Lpt**2*m2**2*xd0t - 18*m1**2*m2**2*xd0t - 9*M1Lpt**2*m2**2*xd0t - m2**4*xd0t + 12*M0Lpt**2*M2Lpt**2*xd0t - 10*m1**2*M2Lpt**2*xd0t - 10*m2**2*M2Lpt**2*xd0t - 
     &                     6*M0Lpt**2*m3**2*xd0t - m1**2*m3**2*xd0t - 3*M1Lpt**2*m3**2*xd0t - m2**2*m3**2*xd0t - 5*m1**2*M3Lpt**2*xd0t + 12*M1Lpt**2*M3Lpt**2*xd0t + 3*m2**2*M3Lpt**2*xd0t - 
     &                     3*m3**2*M3Lpt**2*xd0t - 6*M0Lpt**2*m4**2*xd0t + 3*m1**2*m4**2*xd0t - 7*M1Lpt**2*m4**2*xd0t - 5*m2**2*m4**2*xd0t + M3Lpt**2*m4**2*xd0t + 8*M0Lpt**2*sqrtS**2*xd0t - 
     &                     8*m1**2*sqrtS**2*xd0t + 2*M1Lpt**2*sqrtS**2*xd0t - 4*m2**2*sqrtS**2*xd0t + 8*M2Lpt**2*sqrtS**2*xd0t + 2*m3**2*sqrtS**2*xd0t - 2*M3Lpt**2*sqrtS**2*xd0t + 
     &                     2*m4**2*sqrtS**2*xd0t + 4*sqrtS**4*xd0t + 4*M0Lpt**2*uman*xd0t + 2*m1**2*uman*xd0t - 4*M1Lpt**2*uman*xd0t + 2*m2**2*uman*xd0t - 4*M3Lpt**2*uman*xd0t + 
     &                     4*sqrtS**2*uman*xd0t - 4*M0Lpt**2*m1**2*xd11t + 2*m1**4*xd11t + 4*M0Lpt**2*m2**2*xd11t - 2*m2**4*xd11t - 4*M0Lpt**2*sqrtS**2*xd11t + 4*m2**2*sqrtS**2*xd11t - 
     &                     2*sqrtS**4*xd11t + 4*M0Lpt**2*m1**2*xd12t - 2*m1**4*xd12t + 12*M0Lpt**2*m2**2*xd12t - 8*m1**2*m2**2*xd12t - 6*m2**4*xd12t + 12*M0Lpt**2*m3**2*xd12t - 
     &                     6*m1**2*m3**2*xd12t - 6*m2**2*m3**2*xd12t + 4*M0Lpt**2*m4**2*xd12t - 2*m1**2*m4**2*xd12t - 2*m2**2*m4**2*xd12t - 8*M0Lpt**2*sqrtS**2*xd12t + 
     &                     6*m1**2*sqrtS**2*xd12t + 10*m2**2*sqrtS**2*xd12t + 6*m3**2*sqrtS**2*xd12t + 2*m4**2*sqrtS**2*xd12t - 4*sqrtS**4*xd12t + 8*M0Lpt**2*uman*xd12t - 
     &                     4*m1**2*uman*xd12t - 4*m2**2*uman*xd12t + 4*sqrtS**2*uman*xd12t + 8*M0Lpt**2*sqrtS**2*xd13t - 4*m1**2*sqrtS**2*xd13t - 4*m2**2*sqrtS**2*xd13t + 4*sqrtS**4*xd13t - 
     &                     4*M0Lpt**2*m1**2*xd1t + m1**4*xd1t + 4*m1**2*M1Lpt**2*xd1t + 4*M0Lpt**2*m2**2*xd1t - 22*m1**2*m2**2*xd1t + 4*M1Lpt**2*m2**2*xd1t - 3*m2**4*xd1t + 
     &                     6*M0Lpt**2*m3**2*xd1t - 11*m1**2*m3**2*xd1t - 11*m2**2*m3**2*xd1t - 2*m1**2*M3Lpt**2*xd1t + 2*m2**2*M3Lpt**2*xd1t + 2*M0Lpt**2*m4**2*xd1t - 5*m1**2*m4**2*xd1t - 
     &                     5*m2**2*m4**2*xd1t - 2*M0Lpt**2*sqrtS**2*xd1t - 5*m1**2*sqrtS**2*xd1t + m2**2*sqrtS**2*xd1t + 7*m3**2*sqrtS**2*xd1t - 2*M3Lpt**2*sqrtS**2*xd1t + 
     &                     3*m4**2*sqrtS**2*xd1t + 2*sqrtS**4*xd1t + 4*M0Lpt**2*uman*xd1t + 2*m1**2*uman*xd1t + 2*m2**2*uman*xd1t + 8*M0Lpt**2*m1**2*xd22t - 4*m1**4*xd22t + 
     &                     8*M0Lpt**2*m2**2*xd22t - 8*m1**2*m2**2*xd22t - 4*m2**4*xd22t + 8*M0Lpt**2*m3**2*xd22t - 4*m1**2*m3**2*xd22t - 4*m2**2*m3**2*xd22t + 8*M0Lpt**2*m4**2*xd22t - 
     &                     4*m1**2*m4**2*xd22t - 4*m2**2*m4**2*xd22t - 8*M0Lpt**2*sqrtS**2*xd22t + 8*m1**2*sqrtS**2*xd22t + 8*m2**2*sqrtS**2*xd22t + 4*m3**2*sqrtS**2*xd22t + 
     &                     4*m4**2*sqrtS**2*xd22t - 4*sqrtS**4*xd22t + 8*M0Lpt**2*uman*xd22t - 4*m1**2*uman*xd22t - 4*m2**2*uman*xd22t + 4*sqrtS**2*uman*xd22t - 4*M0Lpt**2*m1**2*xd23t + 
     &                     2*m1**4*xd23t - 12*M0Lpt**2*m2**2*xd23t + 8*m1**2*m2**2*xd23t + 6*m2**4*xd23t - 12*M0Lpt**2*m3**2*xd23t + 6*m1**2*m3**2*xd23t + 6*m2**2*m3**2*xd23t - 
     &                     4*M0Lpt**2*m4**2*xd23t + 2*m1**2*m4**2*xd23t + 2*m2**2*m4**2*xd23t + 8*M0Lpt**2*sqrtS**2*xd23t - 6*m1**2*sqrtS**2*xd23t - 10*m2**2*sqrtS**2*xd23t - 
     &                     6*m3**2*sqrtS**2*xd23t - 2*m4**2*sqrtS**2*xd23t + 4*sqrtS**4*xd23t + 8*M0Lpt**2*uman*xd23t - 4*m1**2*uman*xd23t - 4*m2**2*uman*xd23t + 4*sqrtS**2*uman*xd23t + 
     &                     4*M0Lpt**2*m1**2*xd2t - 3*m1**4*xd2t + 2*m1**2*M1Lpt**2*xd2t + 8*M0Lpt**2*m2**2*xd2t - 26*m1**2*m2**2*xd2t - 8*M1Lpt**2*m2**2*xd2t - 3*m2**4*xd2t + 
     &                     10*M0Lpt**2*m3**2*xd2t - 16*m1**2*m3**2*xd2t - 4*M1Lpt**2*m3**2*xd2t + 2*m2**2*m3**2*xd2t + m3**4*xd2t + 2*m1**2*M3Lpt**2*xd2t + 8*m2**2*M3Lpt**2*xd2t + 
     &                     2*m3**2*M3Lpt**2*xd2t + 10*M0Lpt**2*m4**2*xd2t - 4*m1**2*m4**2*xd2t - 2*M1Lpt**2*m4**2*xd2t - 18*m2**2*m4**2*xd2t + 2*m3**2*m4**2*xd2t + m4**4*xd2t - 
     &                     2*M0Lpt**2*sqrtS**2*xd2t + 2*m1**2*sqrtS**2*xd2t + 2*m2**2*sqrtS**2*xd2t + 8*m3**2*sqrtS**2*xd2t - 6*M3Lpt**2*sqrtS**2*xd2t + 10*m4**2*sqrtS**2*xd2t + 
     &                     12*M0Lpt**2*uman*xd2t + 2*m1**2*uman*xd2t + 2*M1Lpt**2*uman*xd2t - 2*m2**2*uman*xd2t + 2*m3**2*uman*xd2t + 2*M3Lpt**2*uman*xd2t + 2*m4**2*uman*xd2t + 
     &                     6*sqrtS**2*uman*xd2t - 2*uman**2*xd2t + 4*M0Lpt**2*m1**2*xd33t - 2*m1**4*xd33t - 4*M0Lpt**2*m2**2*xd33t + 2*m2**4*xd33t - 4*M0Lpt**2*sqrtS**2*xd33t + 
     &                     4*m1**2*sqrtS**2*xd33t - 2*sqrtS**4*xd33t - 4*M0Lpt**2*m1**2*xd3t + m1**4*xd3t + 2*m1**2*M1Lpt**2*xd3t - 12*M0Lpt**2*m2**2*xd3t - 14*m1**2*m2**2*xd3t - 
     &                     2*M1Lpt**2*m2**2*xd3t + 5*m2**4*xd3t - 6*M0Lpt**2*m3**2*xd3t - m1**2*m3**2*xd3t - m2**2*m3**2*xd3t + 4*m1**2*M3Lpt**2*xd3t + 4*m2**2*M3Lpt**2*xd3t - 
     &                     2*M0Lpt**2*m4**2*xd3t - 7*m1**2*m4**2*xd3t - 7*m2**2*m4**2*xd3t + 6*M0Lpt**2*sqrtS**2*xd3t - 7*m1**2*sqrtS**2*xd3t - 2*M1Lpt**2*sqrtS**2*xd3t - 
     &                     13*m2**2*sqrtS**2*xd3t - m3**2*sqrtS**2*xd3t + 3*m4**2*sqrtS**2*xd3t + 6*sqrtS**4*xd3t + 4*M0Lpt**2*uman*xd3t + 2*m1**2*uman*xd3t + 2*m2**2*uman*xd3t))))) + 
     &      g0Lt*(2*g1Lt*M1Lpt*(g2Rt*M2Lpt*(12*g3Lt*M0Lpt*M3Lpt*(-tman + uman)*xd0t + 
     &               g3Rt*(-2*m1**4*xd1t + 20*m1**2*m2**2*xd1t - 2*m2**4*xd1t + 4*m1**2*m3**2*xd1t + 4*m2**2*m3**2*xd1t + 4*m1**2*m4**2*xd1t + 4*m2**2*m4**2*xd1t + 6*m1**2*sqrtS**2*xd1t + 
     &                  6*m2**2*sqrtS**2*xd1t - 2*m3**2*sqrtS**2*xd1t - 2*m4**2*sqrtS**2*xd1t - 4*sqrtS**4*xd1t - 2*m1**4*xd2t + 6*m1**2*m2**2*xd2t - 4*m2**4*xd2t - 2*m1**2*m3**2*xd2t - 
     &                  10*m2**2*m3**2*xd2t + 4*m1**2*m4**2*xd2t + 8*m2**2*m4**2*xd2t + 5*m1**2*sqrtS**2*xd2t + 7*m2**2*sqrtS**2*xd2t + 3*m3**2*sqrtS**2*xd2t - 3*m4**2*sqrtS**2*xd2t - 
     &                  3*sqrtS**4*xd2t - tman**2*xd2t + uman**2*xd2t + uman*
     &                   (6*xc02t + 6*xc04t + 6*M0Lpt**2*xd0t - 6*m2**2*xd0t + 6*M3Lpt**2*xd0t - 4*m1**2*xd1t - 4*m2**2*xd1t + 2*sqrtS**2*xd1t + m1**2*xd2t - 3*m2**2*xd2t - m3**2*xd2t - 
     &                     m4**2*xd2t - 2*sqrtS**2*xd2t) + tman*(-6*xc02t - 6*xc04t - 6*M0Lpt**2*xd0t + 6*m2**2*xd0t - 6*M3Lpt**2*xd0t - 4*m1**2*xd1t - 4*m2**2*xd1t + 2*sqrtS**2*xd1t - 
     &                     3*m1**2*xd2t + 5*m2**2*xd2t + m3**2*xd2t + m4**2*xd2t + 2*sqrtS**2*xd2t))) + 
     &            g2Lt*(g3Lt*M0Lpt*(4*m1**2*m3**2*xd0t + 8*m2**2*m3**2*xd0t + 4*m2**2*m4**2*xd0t - m1**2*sqrtS**2*xd0t - 3*m2**2*sqrtS**2*xd0t - 3*m3**2*sqrtS**2*xd0t - m4**2*sqrtS**2*xd0t + 
     &                  sqrtS**4*xd0t - 6*m1**2*m3**2*xd11t - 2*m2**2*m3**2*xd11t + 6*m1**2*m4**2*xd11t + 2*m2**2*m4**2*xd11t + 2*m3**2*sqrtS**2*xd11t - 2*m4**2*sqrtS**2*xd11t + 
     &                  8*m1**2*m3**2*xd12t + 16*m1**2*m4**2*xd12t + 8*m2**2*m4**2*xd12t - 6*m1**2*sqrtS**2*xd12t - 2*m2**2*sqrtS**2*xd12t - 2*m3**2*sqrtS**2*xd12t - 6*m4**2*sqrtS**2*xd12t + 
     &                  2*sqrtS**4*xd12t + 4*m1**2*m3**2*xd13t - 4*m2**2*m3**2*xd13t - 4*m1**2*m4**2*xd13t + 4*m2**2*m4**2*xd13t + 6*m1**2*m3**2*xd1t - 2*m2**2*m3**2*xd1t + 
     &                  6*m1**2*m4**2*xd1t + 6*m2**2*m4**2*xd1t - 3*m1**2*sqrtS**2*xd1t - m2**2*sqrtS**2*xd1t - m3**2*sqrtS**2*xd1t - 3*m4**2*sqrtS**2*xd1t + sqrtS**4*xd1t + 
     &                  8*m1**2*m3**2*xd22t + 8*m2**2*m3**2*xd22t + 8*m1**2*m4**2*xd22t + 8*m2**2*m4**2*xd22t - 4*m1**2*sqrtS**2*xd22t - 4*m2**2*sqrtS**2*xd22t - 4*m3**2*sqrtS**2*xd22t - 
     &                  4*m4**2*sqrtS**2*xd22t + 2*sqrtS**4*xd22t + 8*m1**2*m3**2*xd23t + 16*m2**2*m3**2*xd23t + 8*m2**2*m4**2*xd23t - 2*m1**2*sqrtS**2*xd23t - 6*m2**2*sqrtS**2*xd23t - 
     &                  6*m3**2*sqrtS**2*xd23t - 2*m4**2*sqrtS**2*xd23t + 2*sqrtS**4*xd23t + 12*m1**2*m3**2*xd2t + 16*m2**2*m3**2*xd2t + 8*m1**2*m4**2*xd2t + 12*m2**2*m4**2*xd2t - 
     &                  5*m1**2*sqrtS**2*xd2t - 7*m2**2*sqrtS**2*xd2t - 7*m3**2*sqrtS**2*xd2t - 5*m4**2*sqrtS**2*xd2t + 3*sqrtS**4*xd2t + 2*m1**2*m3**2*xd33t + 6*m2**2*m3**2*xd33t - 
     &                  2*m1**2*m4**2*xd33t - 6*m2**2*m4**2*xd33t - 2*m3**2*sqrtS**2*xd33t + 2*m4**2*sqrtS**2*xd33t + uman**2*(-3*xd0t + 2*xd12t + xd1t + 2*xd22t - 6*xd23t - xd2t - 3*xd3t) + 
     &                  6*m1**2*m3**2*xd3t + 14*m2**2*m3**2*xd3t - 2*m1**2*m4**2*xd3t - 2*m2**2*m4**2*xd3t - m1**2*sqrtS**2*xd3t - 3*m2**2*sqrtS**2*xd3t - 5*m3**2*sqrtS**2*xd3t + 
     &                  m4**2*sqrtS**2*xd3t + sqrtS**4*xd3t + tman**2*(xd0t - 6*xd12t - 3*xd1t - 6*xd22t + 2*xd23t - 5*xd2t + xd3t) - 
     &                  uman*(6*xc01t + 6*xc02t - 24*xd00t - m1**2*xd0t - 3*m2**2*xd0t + 6*M2Lpt**2*xd0t - 3*m3**2*xd0t + 6*M3Lpt**2*xd0t - 7*m4**2*xd0t + 2*sqrtS**2*xd0t - 2*m1**2*xd11t + 
     &                     2*m2**2*xd11t - 2*sqrtS**2*xd11t + 2*m1**2*xd12t + 6*m2**2*xd12t + 6*m3**2*xd12t + 2*m4**2*xd12t - 4*sqrtS**2*xd12t + 4*sqrtS**2*xd13t + m1**2*xd1t + 
     &                     3*m2**2*xd1t + 3*m3**2*xd1t + m4**2*xd1t + 4*m1**2*xd22t + 4*m2**2*xd22t + 4*m3**2*xd22t + 4*m4**2*xd22t - 4*sqrtS**2*xd22t - 2*m1**2*xd23t - 6*m2**2*xd23t - 
     &                     6*m3**2*xd23t - 2*m4**2*xd23t + 4*sqrtS**2*xd23t + 3*m1**2*xd2t + m2**2*xd2t + m3**2*xd2t + 3*m4**2*xd2t - 2*sqrtS**2*xd2t + 2*m1**2*xd33t - 2*m2**2*xd33t - 
     &                     2*sqrtS**2*xd33t + m1**2*xd3t - 5*m2**2*xd3t - 3*m3**2*xd3t - m4**2*xd3t) + 
     &                  tman*(6*xc01t + 6*xc02t - 24*xd00t - m1**2*xd0t - 3*m2**2*xd0t + 6*M2Lpt**2*xd0t - 3*m3**2*xd0t + 6*M3Lpt**2*xd0t - 7*m4**2*xd0t + 2*sqrtS**2*xd0t + 2*uman*xd0t - 
     &                     2*m1**2*xd11t + 2*m2**2*xd11t - 2*sqrtS**2*xd11t + 2*m1**2*xd12t + 6*m2**2*xd12t + 6*m3**2*xd12t + 2*m4**2*xd12t - 4*sqrtS**2*xd12t + 4*uman*xd12t + 
     &                     4*sqrtS**2*xd13t + m1**2*xd1t + 3*m2**2*xd1t + 3*m3**2*xd1t + m4**2*xd1t + 2*uman*xd1t + 4*m1**2*xd22t + 4*m2**2*xd22t + 4*m3**2*xd22t + 4*m4**2*xd22t - 
     &                     4*sqrtS**2*xd22t + 4*uman*xd22t - 2*m1**2*xd23t - 6*m2**2*xd23t - 6*m3**2*xd23t - 2*m4**2*xd23t + 4*sqrtS**2*xd23t + 4*uman*xd23t + 3*m1**2*xd2t + m2**2*xd2t + 
     &                     m3**2*xd2t + 3*m4**2*xd2t - 2*sqrtS**2*xd2t + 6*uman*xd2t + 2*m1**2*xd33t - 2*m2**2*xd33t - 2*sqrtS**2*xd33t + m1**2*xd3t - 5*m2**2*xd3t - 3*m3**2*xd3t - 
     &                     m4**2*xd3t + 2*uman*xd3t)) + g3Rt*M3Lpt*(-6*m1**2*m3**2*xd11t - 2*m2**2*m3**2*xd11t + 6*m1**2*m4**2*xd11t + 2*m2**2*m4**2*xd11t + 2*m3**2*sqrtS**2*xd11t - 
     &                  2*m4**2*sqrtS**2*xd11t + 8*m1**2*m3**2*xd12t + 16*m1**2*m4**2*xd12t + 8*m2**2*m4**2*xd12t - 6*m1**2*sqrtS**2*xd12t - 2*m2**2*sqrtS**2*xd12t - 2*m3**2*sqrtS**2*xd12t - 
     &                  6*m4**2*sqrtS**2*xd12t + 2*sqrtS**4*xd12t + 4*m1**2*m3**2*xd13t - 4*m2**2*m3**2*xd13t - 4*m1**2*m4**2*xd13t + 4*m2**2*m4**2*xd13t + 4*m1**2*m3**2*xd1t + 
     &                  8*m1**2*m4**2*xd1t + 4*m2**2*m4**2*xd1t - 3*m1**2*sqrtS**2*xd1t - m2**2*sqrtS**2*xd1t - m3**2*sqrtS**2*xd1t - 3*m4**2*sqrtS**2*xd1t + sqrtS**4*xd1t + 
     &                  8*m1**2*m3**2*xd22t + 8*m2**2*m3**2*xd22t + 8*m1**2*m4**2*xd22t + 8*m2**2*m4**2*xd22t - 4*m1**2*sqrtS**2*xd22t - 4*m2**2*sqrtS**2*xd22t - 4*m3**2*sqrtS**2*xd22t - 
     &                  4*m4**2*sqrtS**2*xd22t + 2*sqrtS**4*xd22t + 8*m1**2*m3**2*xd23t + 16*m2**2*m3**2*xd23t + 8*m2**2*m4**2*xd23t - 2*m1**2*sqrtS**2*xd23t - 6*m2**2*sqrtS**2*xd23t - 
     &                  6*m3**2*sqrtS**2*xd23t - 2*m4**2*sqrtS**2*xd23t + 2*sqrtS**4*xd23t + 8*m1**2*m3**2*xd2t + 8*m2**2*m3**2*xd2t + 8*m1**2*m4**2*xd2t + 8*m2**2*m4**2*xd2t - 
     &                  4*m1**2*sqrtS**2*xd2t - 4*m2**2*sqrtS**2*xd2t - 4*m3**2*sqrtS**2*xd2t - 4*m4**2*sqrtS**2*xd2t + 2*sqrtS**4*xd2t + 2*m1**2*m3**2*xd33t + 6*m2**2*m3**2*xd33t - 
     &                  2*m1**2*m4**2*xd33t - 6*m2**2*m4**2*xd33t - 2*m3**2*sqrtS**2*xd33t + 2*m4**2*sqrtS**2*xd33t + uman**2*(2*xd12t + xd1t + 2*xd22t - 6*xd23t + 2*xd2t - 3*xd3t) + 
     &                  4*m1**2*m3**2*xd3t + 8*m2**2*m3**2*xd3t + 4*m2**2*m4**2*xd3t - m1**2*sqrtS**2*xd3t - 3*m2**2*sqrtS**2*xd3t - 3*m3**2*sqrtS**2*xd3t - m4**2*sqrtS**2*xd3t + 
     &                  sqrtS**4*xd3t + tman**2*(-6*xd0t - 6*xd12t - 3*xd1t - 6*xd22t + 2*xd23t - 6*xd2t + xd3t) - 
     &                  uman*(6*xc01t + 6*xc04t - 24*xd00t + 6*M0Lpt**2*xd0t + 6*M2Lpt**2*xd0t - 2*m1**2*xd11t + 2*m2**2*xd11t - 2*sqrtS**2*xd11t + 2*m1**2*xd12t + 6*m2**2*xd12t + 
     &                     6*m3**2*xd12t + 2*m4**2*xd12t - 4*sqrtS**2*xd12t + 4*sqrtS**2*xd13t + m1**2*xd1t + 3*m2**2*xd1t + 3*m3**2*xd1t + m4**2*xd1t - 2*sqrtS**2*xd1t + 4*m1**2*xd22t + 
     &                     4*m2**2*xd22t + 4*m3**2*xd22t + 4*m4**2*xd22t - 4*sqrtS**2*xd22t - 2*m1**2*xd23t - 6*m2**2*xd23t - 6*m3**2*xd23t - 2*m4**2*xd23t + 4*sqrtS**2*xd23t + 
     &                     4*m1**2*xd2t + 4*m2**2*xd2t + 4*m3**2*xd2t + 4*m4**2*xd2t - 4*sqrtS**2*xd2t + 2*m1**2*xd33t - 2*m2**2*xd33t - 2*sqrtS**2*xd33t - m1**2*xd3t - 3*m2**2*xd3t - 
     &                     3*m3**2*xd3t - m4**2*xd3t + 2*sqrtS**2*xd3t) + tman*
     &                   (6*xc01t + 6*xc04t - 24*xd00t + 6*M0Lpt**2*xd0t + 6*M2Lpt**2*xd0t + 6*uman*xd0t - 2*m1**2*xd11t + 2*m2**2*xd11t - 2*sqrtS**2*xd11t + 2*m1**2*xd12t + 6*m2**2*xd12t + 
     &                     6*m3**2*xd12t + 2*m4**2*xd12t - 4*sqrtS**2*xd12t + 4*uman*xd12t + 4*sqrtS**2*xd13t + m1**2*xd1t + 3*m2**2*xd1t + 3*m3**2*xd1t + m4**2*xd1t - 2*sqrtS**2*xd1t + 
     &                     2*uman*xd1t + 4*m1**2*xd22t + 4*m2**2*xd22t + 4*m3**2*xd22t + 4*m4**2*xd22t - 4*sqrtS**2*xd22t + 4*uman*xd22t - 2*m1**2*xd23t - 6*m2**2*xd23t - 6*m3**2*xd23t - 
     &                     2*m4**2*xd23t + 4*sqrtS**2*xd23t + 4*uman*xd23t + 4*m1**2*xd2t + 4*m2**2*xd2t + 4*m3**2*xd2t + 4*m4**2*xd2t - 4*sqrtS**2*xd2t + 4*uman*xd2t + 2*m1**2*xd33t - 
     &                     2*m2**2*xd33t - 2*sqrtS**2*xd33t - m1**2*xd3t - 3*m2**2*xd3t - 3*m3**2*xd3t - m4**2*xd3t + 2*sqrtS**2*xd3t + 2*uman*xd3t)))) + 
     &         g1Rt*(-2*g2Lt*M2Lpt*(g3Rt*M3Lpt*(-2*m1**4*xd2t - 14*m1**2*m2**2*xd2t - 14*m1**2*m3**2*xd2t - 6*m2**2*m3**2*xd2t + 4*m1**2*m4**2*xd2t + m1**2*sqrtS**2*xd2t - 
     &                  m2**2*sqrtS**2*xd2t + 5*m3**2*sqrtS**2*xd2t - m4**2*sqrtS**2*xd2t + sqrtS**4*xd2t + tman**2*xd2t - uman**2*xd2t + 2*m1**4*xd3t - 20*m1**2*m2**2*xd3t + 2*m2**4*xd3t - 
     &                  4*m1**2*m3**2*xd3t - 4*m2**2*m3**2*xd3t - 4*m1**2*m4**2*xd3t - 4*m2**2*m4**2*xd3t - 6*m1**2*sqrtS**2*xd3t - 6*m2**2*sqrtS**2*xd3t + 2*m3**2*sqrtS**2*xd3t + 
     &                  2*m4**2*sqrtS**2*xd3t + 4*sqrtS**4*xd3t + uman*(-6*xc03t - 6*xc04t - 6*M0Lpt**2*xd0t + 6*m1**2*xd0t - 6*M1Lpt**2*xd0t + 9*m1**2*xd2t + m2**2*xd2t + m3**2*xd2t + 
     &                     m4**2*xd2t + 4*m1**2*xd3t + 4*m2**2*xd3t - 2*sqrtS**2*xd3t) + 
     &                  tman*(6*xc03t + 6*xc04t + 6*M0Lpt**2*xd0t - 6*m1**2*xd0t + 6*M1Lpt**2*xd0t + m1**2*xd2t + 5*m2**2*xd2t - m3**2*xd2t - m4**2*xd2t - 4*sqrtS**2*xd2t + 4*m1**2*xd3t + 
     &                     4*m2**2*xd3t - 2*sqrtS**2*xd3t)) + 2*g3Lt*M0Lpt*(m1**4*xd0t - 10*m1**2*m2**2*xd0t + m2**4*xd0t - 2*m1**2*m3**2*xd0t - 2*m2**2*m3**2*xd0t - 2*m1**2*m4**2*xd0t - 
     &                  2*m2**2*m4**2*xd0t - 3*m1**2*sqrtS**2*xd0t - 3*m2**2*sqrtS**2*xd0t + m3**2*sqrtS**2*xd0t + m4**2*sqrtS**2*xd0t + 2*sqrtS**4*xd0t + m1**4*xd1t - 10*m1**2*m2**2*xd1t + 
     &                  m2**4*xd1t - 2*m1**2*m3**2*xd1t - 2*m2**2*m3**2*xd1t - 2*m1**2*m4**2*xd1t - 2*m2**2*m4**2*xd1t - 3*m1**2*sqrtS**2*xd1t - 3*m2**2*sqrtS**2*xd1t + m3**2*sqrtS**2*xd1t + 
     &                  m4**2*sqrtS**2*xd1t + 2*sqrtS**4*xd1t - 10*m1**2*m2**2*xd2t + 2*m2**4*xd2t - 6*m1**2*m3**2*xd2t + 2*m2**2*m3**2*xd2t - 4*m2**2*m4**2*xd2t - 2*m1**2*sqrtS**2*xd2t - 
     &                  4*m2**2*sqrtS**2*xd2t + m3**2*sqrtS**2*xd2t + m4**2*sqrtS**2*xd2t + 2*sqrtS**4*xd2t + tman**2*xd2t - uman**2*xd2t + m1**4*xd3t - 10*m1**2*m2**2*xd3t + m2**4*xd3t - 
     &                  2*m1**2*m3**2*xd3t - 2*m2**2*m3**2*xd3t - 2*m1**2*m4**2*xd3t - 2*m2**2*m4**2*xd3t - 3*m1**2*sqrtS**2*xd3t - 3*m2**2*sqrtS**2*xd3t + m3**2*sqrtS**2*xd3t + 
     &                  m4**2*sqrtS**2*xd3t + 2*sqrtS**4*xd3t + tman*(3*xc02t + 3*xc03t + 2*m1**2*xd0t + 3*M1Lpt**2*xd0t + 2*m2**2*xd0t + 3*M3Lpt**2*xd0t - 4*sqrtS**2*xd0t + 2*m1**2*xd1t + 
     &                     2*m2**2*xd1t - sqrtS**2*xd1t + 2*m1**2*xd2t - m3**2*xd2t - m4**2*xd2t - 3*sqrtS**2*xd2t + 2*m1**2*xd3t + 2*m2**2*xd3t - sqrtS**2*xd3t) + 
     &                  uman*(-3*xc02t - 3*xc03t + 2*m1**2*xd0t - 3*M1Lpt**2*xd0t + 2*m2**2*xd0t - 3*M3Lpt**2*xd0t + 2*sqrtS**2*xd0t + 2*m1**2*xd1t + 2*m2**2*xd1t - sqrtS**2*xd1t + 
     &                     4*m1**2*xd2t + 2*m2**2*xd2t + m3**2*xd2t + m4**2*xd2t + sqrtS**2*xd2t + 2*m1**2*xd3t + 2*m2**2*xd3t - sqrtS**2*xd3t))) + 
     &            g2Rt*(2*g3Lt*M0Lpt*M3Lpt*(4*m1**2*m3**2*xd0t + 8*m1**2*m4**2*xd0t + 4*m2**2*m4**2*xd0t - 3*m1**2*sqrtS**2*xd0t - m2**2*sqrtS**2*xd0t - m3**2*sqrtS**2*xd0t - 
     &                  3*m4**2*sqrtS**2*xd0t + sqrtS**4*xd0t - 6*m1**2*m3**2*xd11t - 2*m2**2*m3**2*xd11t + 6*m1**2*m4**2*xd11t + 2*m2**2*m4**2*xd11t + 2*m3**2*sqrtS**2*xd11t - 
     &                  2*m4**2*sqrtS**2*xd11t + 8*m1**2*m3**2*xd12t + 16*m1**2*m4**2*xd12t + 8*m2**2*m4**2*xd12t - 6*m1**2*sqrtS**2*xd12t - 2*m2**2*sqrtS**2*xd12t - 2*m3**2*sqrtS**2*xd12t - 
     &                  6*m4**2*sqrtS**2*xd12t + 2*sqrtS**4*xd12t + 4*m1**2*m3**2*xd13t - 4*m2**2*m3**2*xd13t - 4*m1**2*m4**2*xd13t + 4*m2**2*m4**2*xd13t - 2*m1**2*m3**2*xd1t - 
     &                  2*m2**2*m3**2*xd1t + 14*m1**2*m4**2*xd1t + 6*m2**2*m4**2*xd1t - 3*m1**2*sqrtS**2*xd1t - m2**2*sqrtS**2*xd1t + m3**2*sqrtS**2*xd1t - 5*m4**2*sqrtS**2*xd1t + 
     &                  sqrtS**4*xd1t + 8*m1**2*m3**2*xd22t + 8*m2**2*m3**2*xd22t + 8*m1**2*m4**2*xd22t + 8*m2**2*m4**2*xd22t - 4*m1**2*sqrtS**2*xd22t - 4*m2**2*sqrtS**2*xd22t - 
     &                  4*m3**2*sqrtS**2*xd22t - 4*m4**2*sqrtS**2*xd22t + 2*sqrtS**4*xd22t + 8*m1**2*m3**2*xd23t + 16*m2**2*m3**2*xd23t + 8*m2**2*m4**2*xd23t - 2*m1**2*sqrtS**2*xd23t - 
     &                  6*m2**2*sqrtS**2*xd23t - 6*m3**2*sqrtS**2*xd23t - 2*m4**2*sqrtS**2*xd23t + 2*sqrtS**4*xd23t + 12*m1**2*m3**2*xd2t + 8*m2**2*m3**2*xd2t + 16*m1**2*m4**2*xd2t + 
     &                  12*m2**2*m4**2*xd2t - 7*m1**2*sqrtS**2*xd2t - 5*m2**2*sqrtS**2*xd2t - 5*m3**2*sqrtS**2*xd2t - 7*m4**2*sqrtS**2*xd2t + 3*sqrtS**4*xd2t + 2*m1**2*m3**2*xd33t + 
     &                  6*m2**2*m3**2*xd33t - 2*m1**2*m4**2*xd33t - 6*m2**2*m4**2*xd33t - 2*m3**2*sqrtS**2*xd33t + 2*m4**2*sqrtS**2*xd33t + 
     &                  uman**2*(xd0t + 2*xd12t + xd1t + 2*xd22t - 6*xd23t + 3*xd2t - 3*xd3t) + 6*m1**2*m3**2*xd3t + 6*m2**2*m3**2*xd3t - 2*m1**2*m4**2*xd3t + 6*m2**2*m4**2*xd3t - 
     &                  m1**2*sqrtS**2*xd3t - 3*m2**2*sqrtS**2*xd3t - 3*m3**2*sqrtS**2*xd3t - m4**2*sqrtS**2*xd3t + sqrtS**4*xd3t + 
     &                  tman**2*(-3*xd0t - 6*xd12t - 3*xd1t - 6*xd22t + 2*xd23t - 9*xd2t + xd3t) - 
     &                  uman*(6*xc01t + 6*xc03t - 24*xd00t + m1**2*xd0t + 6*M1Lpt**2*xd0t + 3*m2**2*xd0t + 6*M2Lpt**2*xd0t - 3*m3**2*xd0t + m4**2*xd0t - 2*sqrtS**2*xd0t - 2*m1**2*xd11t + 
     &                     2*m2**2*xd11t - 2*sqrtS**2*xd11t + 2*m1**2*xd12t + 6*m2**2*xd12t + 6*m3**2*xd12t + 2*m4**2*xd12t - 4*sqrtS**2*xd12t + 4*sqrtS**2*xd13t - m1**2*xd1t + 
     &                     5*m2**2*xd1t + 3*m3**2*xd1t + m4**2*xd1t - 4*sqrtS**2*xd1t + 4*m1**2*xd22t + 4*m2**2*xd22t + 4*m3**2*xd22t + 4*m4**2*xd22t - 4*sqrtS**2*xd22t - 2*m1**2*xd23t - 
     &                     6*m2**2*xd23t - 6*m3**2*xd23t - 2*m4**2*xd23t + 4*sqrtS**2*xd23t + 5*m1**2*xd2t + 7*m2**2*xd2t + 7*m3**2*xd2t + 5*m4**2*xd2t - 6*sqrtS**2*xd2t + 2*m1**2*xd33t - 
     &                     2*m2**2*xd33t - 2*sqrtS**2*xd33t - m1**2*xd3t - 3*m2**2*xd3t - 3*m3**2*xd3t - m4**2*xd3t + 4*sqrtS**2*xd3t) + 
     &                  tman*(6*xc01t + 6*xc03t - 24*xd00t + m1**2*xd0t + 6*M1Lpt**2*xd0t + 3*m2**2*xd0t + 6*M2Lpt**2*xd0t - 3*m3**2*xd0t + m4**2*xd0t - 2*sqrtS**2*xd0t + 2*uman*xd0t - 
     &                     2*m1**2*xd11t + 2*m2**2*xd11t - 2*sqrtS**2*xd11t + 2*m1**2*xd12t + 6*m2**2*xd12t + 6*m3**2*xd12t + 2*m4**2*xd12t - 4*sqrtS**2*xd12t + 4*uman*xd12t + 
     &                     4*sqrtS**2*xd13t - m1**2*xd1t + 5*m2**2*xd1t + 3*m3**2*xd1t + m4**2*xd1t - 4*sqrtS**2*xd1t + 2*uman*xd1t + 4*m1**2*xd22t + 4*m2**2*xd22t + 4*m3**2*xd22t + 
     &                     4*m4**2*xd22t - 4*sqrtS**2*xd22t + 4*uman*xd22t - 2*m1**2*xd23t - 6*m2**2*xd23t - 6*m3**2*xd23t - 2*m4**2*xd23t + 4*sqrtS**2*xd23t + 4*uman*xd23t + 5*m1**2*xd2t + 
     &                     7*m2**2*xd2t + 7*m3**2*xd2t + 5*m4**2*xd2t - 6*sqrtS**2*xd2t + 6*uman*xd2t + 2*m1**2*xd33t - 2*m2**2*xd33t - 2*sqrtS**2*xd33t - m1**2*xd3t - 3*m2**2*xd3t - 
     &                     3*m3**2*xd3t - m4**2*xd3t + 4*sqrtS**2*xd3t + 2*uman*xd3t)) + 
     &               g3Rt*(-2*m1**4*xc01t + 20*m1**2*m2**2*xc01t - 2*m2**4*xc01t + 4*m1**2*m3**2*xc01t + 4*m2**2*m3**2*xc01t + 4*m1**2*m4**2*xc01t + 4*m2**2*m4**2*xc01t + 
     &                  6*m1**2*sqrtS**2*xc01t + 6*m2**2*sqrtS**2*xc01t - 2*m3**2*sqrtS**2*xc01t - 2*m4**2*sqrtS**2*xc01t - 4*sqrtS**4*xc01t + 4*m1**2*m3**2*xc02t + 8*m1**2*m4**2*xc02t + 
     &                  4*m2**2*m4**2*xc02t - 3*m1**2*sqrtS**2*xc02t - m2**2*sqrtS**2*xc02t - m3**2*sqrtS**2*xc02t - 3*m4**2*sqrtS**2*xc02t + sqrtS**4*xc02t + 4*m1**2*m3**2*xc03t + 
     &                  8*m2**2*m3**2*xc03t + 4*m2**2*m4**2*xc03t - m1**2*sqrtS**2*xc03t - 3*m2**2*sqrtS**2*xc03t - 3*m3**2*sqrtS**2*xc03t - m4**2*sqrtS**2*xc03t + sqrtS**4*xc03t - 
     &                  8*m1**2*m3**2*xc04t + 8*m2**2*m3**2*xc04t + 2*m1**2*sqrtS**2*xc04t - 2*m2**2*sqrtS**2*xc04t - 12*m1**2*m3**2*xc111t + 12*m2**2*m3**2*xc111t - 4*m1**2*m4**2*xc111t + 
     &                  4*m2**2*m4**2*xc111t + 4*m1**2*sqrtS**2*xc111t - 4*m2**2*sqrtS**2*xc111t + 2*m1**4*xc11t - 2*m2**4*xc11t - 12*m1**2*m3**2*xc11t + 12*m2**2*m3**2*xc11t - 
     &                  8*m1**2*m4**2*xc11t + 8*m2**2*m4**2*xc11t + 4*m1**2*sqrtS**2*xc11t - 4*m2**2*sqrtS**2*xc11t - 32*m1**2*m3**2*xc121t + 32*m2**2*m3**2*xc121t + 
     &                  8*m1**2*sqrtS**2*xc121t - 8*m2**2*sqrtS**2*xc121t - 6*m1**2*m3**2*xc12t - 2*m2**2*m3**2*xc12t + 6*m1**2*m4**2*xc12t + 2*m2**2*m4**2*xc12t + 2*m3**2*sqrtS**2*xc12t - 
     &                  2*m4**2*sqrtS**2*xc12t + 2*m1**4*xc13t - 6*m1**2*m2**2*xc13t + 4*m2**4*xc13t + 6*m1**2*m3**2*xc13t + 18*m2**2*m3**2*xc13t - 4*m1**2*m4**2*xc13t - 
     &                  4*m2**2*m4**2*xc13t - 6*m1**2*sqrtS**2*xc13t - 10*m2**2*sqrtS**2*xc13t - 6*m3**2*sqrtS**2*xc13t + 2*m4**2*sqrtS**2*xc13t + 4*sqrtS**4*xc13t - 24*m1**2*m3**2*xc21t + 
     &                  24*m2**2*m3**2*xc21t + 8*m1**2*m4**2*xc21t - 8*m2**2*m4**2*xc21t + 4*m1**2*sqrtS**2*xc21t - 4*m2**2*sqrtS**2*xc21t - 16*m1**2*m3**2*xc221t + 16*m2**2*m3**2*xc221t + 
     &                  16*m1**2*m4**2*xc221t - 16*m2**2*m4**2*xc221t - 2*m1**4*xc22t - 14*m1**2*m2**2*xc22t - 10*m1**2*m3**2*xc22t - 6*m2**2*m3**2*xc22t + 12*m1**2*m4**2*xc22t + 
     &                  4*m2**2*m4**2*xc22t - 2*m1**2*sqrtS**2*xc22t - 2*m2**2*sqrtS**2*xc22t + 4*m3**2*sqrtS**2*xc22t - 4*m4**2*sqrtS**2*xc22t + 2*sqrtS**4*xc22t + 2*m1**2*m3**2*xc23t + 
     &                  6*m2**2*m3**2*xc23t - 2*m1**2*m4**2*xc23t - 6*m2**2*m4**2*xc23t - 2*m3**2*sqrtS**2*xc23t + 2*m4**2*sqrtS**2*xc23t - 2*M0Lpt**2*m1**4*xd0t + 
     &                  20*M0Lpt**2*m1**2*m2**2*xd0t - 2*M0Lpt**2*m2**4*xd0t - 2*m1**4*M2Lpt**2*xd0t + 20*m1**2*m2**2*M2Lpt**2*xd0t - 2*m2**4*M2Lpt**2*xd0t + 12*M0Lpt**2*m1**2*m3**2*xd0t - 
     &                  4*m1**4*m3**2*xd0t + 4*m1**2*M1Lpt**2*m3**2*xd0t + 12*M0Lpt**2*m2**2*m3**2*xd0t - 12*m1**2*m2**2*m3**2*xd0t + 8*M1Lpt**2*m2**2*m3**2*xd0t + 
     &                  4*m1**2*M2Lpt**2*m3**2*xd0t + 4*m2**2*M2Lpt**2*m3**2*xd0t + 4*m1**2*m3**2*M3Lpt**2*xd0t + 12*M0Lpt**2*m1**2*m4**2*xd0t + 12*M0Lpt**2*m2**2*m4**2*xd0t - 
     &                  12*m1**2*m2**2*m4**2*xd0t + 4*M1Lpt**2*m2**2*m4**2*xd0t - 4*m2**4*m4**2*xd0t + 4*m1**2*M2Lpt**2*m4**2*xd0t + 4*m2**2*M2Lpt**2*m4**2*xd0t + 
     &                  8*m1**2*M3Lpt**2*m4**2*xd0t + 4*m2**2*M3Lpt**2*m4**2*xd0t + 2*M0Lpt**2*m1**2*sqrtS**2*xd0t + m1**4*sqrtS**2*xd0t - m1**2*M1Lpt**2*sqrtS**2*xd0t + 
     &                  2*M0Lpt**2*m2**2*sqrtS**2*xd0t + 6*m1**2*m2**2*sqrtS**2*xd0t - 3*M1Lpt**2*m2**2*sqrtS**2*xd0t + m2**4*sqrtS**2*xd0t + 6*m1**2*M2Lpt**2*sqrtS**2*xd0t + 
     &                  6*m2**2*M2Lpt**2*sqrtS**2*xd0t - 6*M0Lpt**2*m3**2*sqrtS**2*xd0t + 3*m1**2*m3**2*sqrtS**2*xd0t - 3*M1Lpt**2*m3**2*sqrtS**2*xd0t + m2**2*m3**2*sqrtS**2*xd0t - 
     &                  2*M2Lpt**2*m3**2*sqrtS**2*xd0t - 3*m1**2*M3Lpt**2*sqrtS**2*xd0t - m2**2*M3Lpt**2*sqrtS**2*xd0t - m3**2*M3Lpt**2*sqrtS**2*xd0t - 6*M0Lpt**2*m4**2*sqrtS**2*xd0t + 
     &                  m1**2*m4**2*sqrtS**2*xd0t - M1Lpt**2*m4**2*sqrtS**2*xd0t + 3*m2**2*m4**2*sqrtS**2*xd0t - 2*M2Lpt**2*m4**2*sqrtS**2*xd0t - 3*M3Lpt**2*m4**2*sqrtS**2*xd0t - 
     &                  2*M0Lpt**2*sqrtS**4*xd0t - m1**2*sqrtS**4*xd0t + M1Lpt**2*sqrtS**4*xd0t - m2**2*sqrtS**4*xd0t - 4*M2Lpt**2*sqrtS**4*xd0t + M3Lpt**2*sqrtS**4*xd0t - 
     &                  12*M0Lpt**2*m1**2*m3**2*xd11t + 6*m1**4*m3**2*xd11t - 4*M0Lpt**2*m2**2*m3**2*xd11t + 8*m1**2*m2**2*m3**2*xd11t + 2*m2**4*m3**2*xd11t + 12*M0Lpt**2*m1**2*m4**2*xd11t - 
     &                  6*m1**4*m4**2*xd11t + 4*M0Lpt**2*m2**2*m4**2*xd11t - 8*m1**2*m2**2*m4**2*xd11t - 2*m2**4*m4**2*xd11t + 4*M0Lpt**2*m3**2*sqrtS**2*xd11t - 
     &                  8*m1**2*m3**2*sqrtS**2*xd11t - 4*m2**2*m3**2*sqrtS**2*xd11t - 4*M0Lpt**2*m4**2*sqrtS**2*xd11t + 8*m1**2*m4**2*sqrtS**2*xd11t + 4*m2**2*m4**2*sqrtS**2*xd11t + 
     &                  2*m3**2*sqrtS**4*xd11t - 2*m4**2*sqrtS**4*xd11t + 16*M0Lpt**2*m1**2*m3**2*xd12t - 8*m1**4*m3**2*xd12t - 8*m1**2*m2**2*m3**2*xd12t + 32*M0Lpt**2*m1**2*m4**2*xd12t - 
     &                  16*m1**4*m4**2*xd12t + 16*M0Lpt**2*m2**2*m4**2*xd12t - 24*m1**2*m2**2*m4**2*xd12t - 8*m2**4*m4**2*xd12t - 12*M0Lpt**2*m1**2*sqrtS**2*xd12t + 6*m1**4*sqrtS**2*xd12t - 
     &                  4*M0Lpt**2*m2**2*sqrtS**2*xd12t + 8*m1**2*m2**2*sqrtS**2*xd12t + 2*m2**4*sqrtS**2*xd12t - 4*M0Lpt**2*m3**2*sqrtS**2*xd12t + 10*m1**2*m3**2*sqrtS**2*xd12t + 
     &                  2*m2**2*m3**2*sqrtS**2*xd12t - 12*M0Lpt**2*m4**2*sqrtS**2*xd12t + 22*m1**2*m4**2*sqrtS**2*xd12t + 14*m2**2*m4**2*sqrtS**2*xd12t + 4*M0Lpt**2*sqrtS**4*xd12t - 
     &                  8*m1**2*sqrtS**4*xd12t - 4*m2**2*sqrtS**4*xd12t - 2*m3**2*sqrtS**4*xd12t - 6*m4**2*sqrtS**4*xd12t + 2*sqrtS**6*xd12t + 8*M0Lpt**2*m1**2*m3**2*xd13t - 
     &                  4*m1**4*m3**2*xd13t - 8*M0Lpt**2*m2**2*m3**2*xd13t + 4*m2**4*m3**2*xd13t - 8*M0Lpt**2*m1**2*m4**2*xd13t + 4*m1**4*m4**2*xd13t + 8*M0Lpt**2*m2**2*m4**2*xd13t - 
     &                  4*m2**4*m4**2*xd13t + 4*m1**2*m3**2*sqrtS**2*xd13t - 4*m2**2*m3**2*sqrtS**2*xd13t - 4*m1**2*m4**2*sqrtS**2*xd13t + 4*m2**2*m4**2*sqrtS**2*xd13t - 
     &                  2*M0Lpt**2*m1**4*xd1t + 2*m1**4*M1Lpt**2*xd1t + 20*M0Lpt**2*m1**2*m2**2*xd1t - 20*m1**2*M1Lpt**2*m2**2*xd1t - 2*M0Lpt**2*m2**4*xd1t + 2*M1Lpt**2*m2**4*xd1t + 
     &                  8*M0Lpt**2*m1**2*m3**2*xd1t - 8*m1**4*m3**2*xd1t - 2*m1**2*M1Lpt**2*m3**2*xd1t + 24*m1**2*m2**2*m3**2*xd1t - 6*M1Lpt**2*m2**2*m3**2*xd1t + 4*m1**2*m3**4*xd1t + 
     &                  4*m2**2*m3**4*xd1t - 6*m1**2*m3**2*M3Lpt**2*xd1t - 2*m2**2*m3**2*M3Lpt**2*xd1t + 24*M0Lpt**2*m1**2*m4**2*xd1t - 6*m1**4*m4**2*xd1t - 6*m1**2*M1Lpt**2*m4**2*xd1t + 
     &                  16*M0Lpt**2*m2**2*m4**2*xd1t - 20*m1**2*m2**2*m4**2*xd1t - 2*M1Lpt**2*m2**2*m4**2*xd1t - 6*m2**4*m4**2*xd1t + 4*m1**2*m3**2*m4**2*xd1t + 4*m2**2*m3**2*m4**2*xd1t + 
     &                  6*m1**2*M3Lpt**2*m4**2*xd1t + 2*m2**2*M3Lpt**2*m4**2*xd1t + 3*m1**4*sqrtS**2*xd1t - 6*m1**2*M1Lpt**2*sqrtS**2*xd1t + 4*M0Lpt**2*m2**2*sqrtS**2*xd1t + 
     &                  4*m1**2*m2**2*sqrtS**2*xd1t - 6*M1Lpt**2*m2**2*sqrtS**2*xd1t + m2**4*sqrtS**2*xd1t - 2*M0Lpt**2*m3**2*sqrtS**2*xd1t + 11*m1**2*m3**2*sqrtS**2*xd1t + 
     &                  2*M1Lpt**2*m3**2*sqrtS**2*xd1t + 5*m2**2*m3**2*sqrtS**2*xd1t - 2*m3**4*sqrtS**2*xd1t + 2*m3**2*M3Lpt**2*sqrtS**2*xd1t - 10*M0Lpt**2*m4**2*sqrtS**2*xd1t + 
     &                  11*m1**2*m4**2*sqrtS**2*xd1t + 2*M1Lpt**2*m4**2*sqrtS**2*xd1t + 9*m2**2*m4**2*sqrtS**2*xd1t - 2*m3**2*m4**2*sqrtS**2*xd1t - 2*M3Lpt**2*m4**2*sqrtS**2*xd1t - 
     &                  2*M0Lpt**2*sqrtS**4*xd1t - 4*m1**2*sqrtS**4*xd1t + 4*M1Lpt**2*sqrtS**4*xd1t - 2*m2**2*sqrtS**4*xd1t - 5*m3**2*sqrtS**4*xd1t - 3*m4**2*sqrtS**4*xd1t + sqrtS**6*xd1t + 
     &                  16*M0Lpt**2*m1**2*m3**2*xd22t - 8*m1**4*m3**2*xd22t + 16*M0Lpt**2*m2**2*m3**2*xd22t - 16*m1**2*m2**2*m3**2*xd22t - 8*m2**4*m3**2*xd22t + 
     &                  16*M0Lpt**2*m1**2*m4**2*xd22t - 8*m1**4*m4**2*xd22t + 16*M0Lpt**2*m2**2*m4**2*xd22t - 16*m1**2*m2**2*m4**2*xd22t - 8*m2**4*m4**2*xd22t - 
     &                  8*M0Lpt**2*m1**2*sqrtS**2*xd22t + 4*m1**4*sqrtS**2*xd22t - 8*M0Lpt**2*m2**2*sqrtS**2*xd22t + 8*m1**2*m2**2*sqrtS**2*xd22t + 4*m2**4*sqrtS**2*xd22t - 
     &                  8*M0Lpt**2*m3**2*sqrtS**2*xd22t + 12*m1**2*m3**2*sqrtS**2*xd22t + 12*m2**2*m3**2*sqrtS**2*xd22t - 8*M0Lpt**2*m4**2*sqrtS**2*xd22t + 12*m1**2*m4**2*sqrtS**2*xd22t + 
     &                  12*m2**2*m4**2*sqrtS**2*xd22t + 4*M0Lpt**2*sqrtS**4*xd22t - 6*m1**2*sqrtS**4*xd22t - 6*m2**2*sqrtS**4*xd22t - 4*m3**2*sqrtS**4*xd22t - 4*m4**2*sqrtS**4*xd22t + 
     &                  2*sqrtS**6*xd22t + 16*M0Lpt**2*m1**2*m3**2*xd23t - 8*m1**4*m3**2*xd23t + 32*M0Lpt**2*m2**2*m3**2*xd23t - 24*m1**2*m2**2*m3**2*xd23t - 16*m2**4*m3**2*xd23t + 
     &                  16*M0Lpt**2*m2**2*m4**2*xd23t - 8*m1**2*m2**2*m4**2*xd23t - 8*m2**4*m4**2*xd23t - 4*M0Lpt**2*m1**2*sqrtS**2*xd23t + 2*m1**4*sqrtS**2*xd23t - 
     &                  12*M0Lpt**2*m2**2*sqrtS**2*xd23t + 8*m1**2*m2**2*sqrtS**2*xd23t + 6*m2**4*sqrtS**2*xd23t - 12*M0Lpt**2*m3**2*sqrtS**2*xd23t + 14*m1**2*m3**2*sqrtS**2*xd23t + 
     &                  22*m2**2*m3**2*sqrtS**2*xd23t - 4*M0Lpt**2*m4**2*sqrtS**2*xd23t + 2*m1**2*m4**2*sqrtS**2*xd23t + 10*m2**2*m4**2*sqrtS**2*xd23t + 4*M0Lpt**2*sqrtS**4*xd23t - 
     &                  4*m1**2*sqrtS**4*xd23t - 8*m2**2*sqrtS**4*xd23t - 6*m3**2*sqrtS**4*xd23t - 2*m4**2*sqrtS**4*xd23t + 2*sqrtS**6*xd23t + 2*m1**4*M1Lpt**2*xd2t + 
     &                  20*M0Lpt**2*m1**2*m2**2*xd2t - 6*m1**2*M1Lpt**2*m2**2*xd2t - 4*M0Lpt**2*m2**4*xd2t + 4*M1Lpt**2*m2**4*xd2t + 36*M0Lpt**2*m1**2*m3**2*xd2t - 14*m1**4*m3**2*xd2t + 
     &                  6*m1**2*M1Lpt**2*m3**2*xd2t + 20*M0Lpt**2*m2**2*m3**2*xd2t - 22*m1**2*m2**2*m3**2*xd2t + 18*M1Lpt**2*m2**2*m3**2*xd2t - 12*m2**4*m3**2*xd2t - 2*m1**2*m3**4*xd2t - 
     &                  10*m2**2*m3**4*xd2t - 2*m1**4*M3Lpt**2*xd2t - 14*m1**2*m2**2*M3Lpt**2*xd2t - 10*m1**2*m3**2*M3Lpt**2*xd2t - 6*m2**2*m3**2*M3Lpt**2*xd2t + 
     &                  24*M0Lpt**2*m1**2*m4**2*xd2t - 6*m1**4*m4**2*xd2t - 4*m1**2*M1Lpt**2*m4**2*xd2t + 32*M0Lpt**2*m2**2*m4**2*xd2t - 14*m1**2*m2**2*m4**2*xd2t - 
     &                  4*M1Lpt**2*m2**2*m4**2*xd2t - 12*m2**4*m4**2*xd2t + 18*m1**2*m3**2*m4**2*xd2t + 14*m2**2*m3**2*m4**2*xd2t + 12*m1**2*M3Lpt**2*m4**2*xd2t + 
     &                  4*m2**2*M3Lpt**2*m4**2*xd2t - 4*m1**2*m4**4*xd2t - 8*M0Lpt**2*m1**2*sqrtS**2*xd2t + 5*m1**4*sqrtS**2*xd2t - 6*m1**2*M1Lpt**2*sqrtS**2*xd2t - 
     &                  4*M0Lpt**2*m2**2*sqrtS**2*xd2t + 14*m1**2*m2**2*sqrtS**2*xd2t - 10*M1Lpt**2*m2**2*sqrtS**2*xd2t + 5*m2**4*sqrtS**2*xd2t - 14*M0Lpt**2*m3**2*sqrtS**2*xd2t + 
     &                  20*m1**2*m3**2*sqrtS**2*xd2t - 6*M1Lpt**2*m3**2*sqrtS**2*xd2t + 20*m2**2*m3**2*sqrtS**2*xd2t + 3*m3**4*sqrtS**2*xd2t - 2*m1**2*M3Lpt**2*sqrtS**2*xd2t - 
     &                  2*m2**2*M3Lpt**2*sqrtS**2*xd2t + 4*m3**2*M3Lpt**2*sqrtS**2*xd2t - 14*M0Lpt**2*m4**2*sqrtS**2*xd2t + 12*m1**2*m4**2*sqrtS**2*xd2t + 2*M1Lpt**2*m4**2*sqrtS**2*xd2t + 
     &                  16*m2**2*m4**2*sqrtS**2*xd2t - 8*m3**2*m4**2*sqrtS**2*xd2t - 4*M3Lpt**2*m4**2*sqrtS**2*xd2t + m4**4*sqrtS**2*xd2t + 2*M0Lpt**2*sqrtS**4*xd2t - 7*m1**2*sqrtS**4*xd2t + 
     &                  4*M1Lpt**2*sqrtS**4*xd2t - 7*m2**2*sqrtS**4*xd2t - 7*m3**2*sqrtS**4*xd2t + 2*M3Lpt**2*sqrtS**4*xd2t - 5*m4**2*sqrtS**4*xd2t + 2*sqrtS**6*xd2t + 2*tman**3*xd2t + 
     &                  4*M0Lpt**2*m1**2*m3**2*xd33t - 2*m1**4*m3**2*xd33t + 12*M0Lpt**2*m2**2*m3**2*xd33t - 8*m1**2*m2**2*m3**2*xd33t - 6*m2**4*m3**2*xd33t - 4*M0Lpt**2*m1**2*m4**2*xd33t + 
     &                  2*m1**4*m4**2*xd33t - 12*M0Lpt**2*m2**2*m4**2*xd33t + 8*m1**2*m2**2*m4**2*xd33t + 6*m2**4*m4**2*xd33t - 4*M0Lpt**2*m3**2*sqrtS**2*xd33t + 
     &                  4*m1**2*m3**2*sqrtS**2*xd33t + 8*m2**2*m3**2*sqrtS**2*xd33t + 4*M0Lpt**2*m4**2*sqrtS**2*xd33t - 4*m1**2*m4**2*sqrtS**2*xd33t - 8*m2**2*m4**2*sqrtS**2*xd33t - 
     &                  2*m3**2*sqrtS**4*xd33t + 2*m4**2*sqrtS**4*xd33t - 2*M0Lpt**2*m1**4*xd3t + 20*M0Lpt**2*m1**2*m2**2*xd3t - 2*M0Lpt**2*m2**4*xd3t + 16*M0Lpt**2*m1**2*m3**2*xd3t - 
     &                  6*m1**4*m3**2*xd3t + 2*m1**2*M1Lpt**2*m3**2*xd3t + 24*M0Lpt**2*m2**2*m3**2*xd3t - 20*m1**2*m2**2*m3**2*xd3t + 6*M1Lpt**2*m2**2*m3**2*xd3t - 6*m2**4*m3**2*xd3t + 
     &                  2*m1**4*M3Lpt**2*xd3t - 20*m1**2*m2**2*M3Lpt**2*xd3t + 2*m2**4*M3Lpt**2*xd3t - 2*m1**2*m3**2*M3Lpt**2*xd3t - 6*m2**2*m3**2*M3Lpt**2*xd3t - 
     &                  2*m1**2*M1Lpt**2*m4**2*xd3t + 8*M0Lpt**2*m2**2*m4**2*xd3t + 24*m1**2*m2**2*m4**2*xd3t - 6*M1Lpt**2*m2**2*m4**2*xd3t - 8*m2**4*m4**2*xd3t + 4*m1**2*m3**2*m4**2*xd3t + 
     &                  4*m2**2*m3**2*m4**2*xd3t - 6*m1**2*M3Lpt**2*m4**2*xd3t - 2*m2**2*M3Lpt**2*m4**2*xd3t + 4*m1**2*m4**4*xd3t + 4*m2**2*m4**4*xd3t + 4*M0Lpt**2*m1**2*sqrtS**2*xd3t + 
     &                  m1**4*sqrtS**2*xd3t + 4*m1**2*m2**2*sqrtS**2*xd3t + 3*m2**4*sqrtS**2*xd3t - 10*M0Lpt**2*m3**2*sqrtS**2*xd3t + 9*m1**2*m3**2*sqrtS**2*xd3t - 
     &                  2*M1Lpt**2*m3**2*sqrtS**2*xd3t + 11*m2**2*m3**2*sqrtS**2*xd3t - 6*m1**2*M3Lpt**2*sqrtS**2*xd3t - 6*m2**2*M3Lpt**2*sqrtS**2*xd3t + 2*m3**2*M3Lpt**2*sqrtS**2*xd3t - 
     &                  2*M0Lpt**2*m4**2*sqrtS**2*xd3t + 5*m1**2*m4**2*sqrtS**2*xd3t + 2*M1Lpt**2*m4**2*sqrtS**2*xd3t + 11*m2**2*m4**2*sqrtS**2*xd3t - 2*m3**2*m4**2*sqrtS**2*xd3t + 
     &                  2*M3Lpt**2*m4**2*sqrtS**2*xd3t - 2*m4**4*sqrtS**2*xd3t - 2*M0Lpt**2*sqrtS**4*xd3t - 2*m1**2*sqrtS**4*xd3t - 4*m2**2*sqrtS**4*xd3t - 3*m3**2*sqrtS**4*xd3t + 
     &                  4*M3Lpt**2*sqrtS**4*xd3t - 5*m4**2*sqrtS**4*xd3t + sqrtS**6*xd3t + 
     &                  uman**2*(xc02t - 3*xc03t - 4*xc04t - 2*xc11t - 16*xc121t - 4*xc13t - 8*xc21t - 2*M0Lpt**2*xd0t + 3*m1**2*xd0t - 3*M1Lpt**2*xd0t - m2**2*xd0t + M3Lpt**2*xd0t + 
     &                     4*M0Lpt**2*xd12t - 2*m1**2*xd12t - 2*m2**2*xd12t + 2*sqrtS**2*xd12t + 2*M0Lpt**2*xd1t - m1**2*xd1t - m2**2*xd1t + sqrtS**2*xd1t + 4*M0Lpt**2*xd22t - 
     &                     2*m1**2*xd22t - 2*m2**2*xd22t + 2*sqrtS**2*xd22t - 12*M0Lpt**2*xd23t + 6*m1**2*xd23t + 6*m2**2*xd23t - 6*sqrtS**2*xd23t + 4*M0Lpt**2*xd2t + m1**2*xd2t - 
     &                     4*M1Lpt**2*xd2t - 3*m2**2*xd2t + m3**2*xd2t + m4**2*xd2t + 2*sqrtS**2*xd2t - 6*M0Lpt**2*xd3t + 3*m1**2*xd3t + 3*m2**2*xd3t - 3*sqrtS**2*xd3t) + 
     &                  tman**2*(3*xc02t + 7*xc03t + 4*xc04t + 2*xc11t + 16*xc121t + 2*xc13t + 8*xc21t - 2*xc22t - 2*M0Lpt**2*xd0t + 3*m1**2*xd0t + 7*M1Lpt**2*xd0t + 7*m2**2*xd0t + 
     &                     3*M3Lpt**2*xd0t - 8*sqrtS**2*xd0t - 12*M0Lpt**2*xd12t + 6*m1**2*xd12t + 6*m2**2*xd12t - 6*sqrtS**2*xd12t - 6*M0Lpt**2*xd1t + 7*m1**2*xd1t + 7*m2**2*xd1t - 
     &                     5*sqrtS**2*xd1t - 12*M0Lpt**2*xd22t + 6*m1**2*xd22t + 6*m2**2*xd22t - 6*sqrtS**2*xd22t + 4*M0Lpt**2*xd23t - 2*m1**2*xd23t - 2*m2**2*xd23t + 2*sqrtS**2*xd23t - 
     &                     16*M0Lpt**2*xd2t + 9*m1**2*xd2t + 2*M1Lpt**2*xd2t + 9*m2**2*xd2t - 3*m3**2*xd2t - 2*M3Lpt**2*xd2t - 3*m4**2*xd2t - 12*sqrtS**2*xd2t + 2*M0Lpt**2*xd3t + 
     &                     3*m1**2*xd3t + 3*m2**2*xd3t - sqrtS**2*xd3t) + uman*
     &                   (-12*xb01t - 12*xb02t + 48*xc001t - 12*M0Lpt**2*xc01t + 2*m1**2*xc01t + 2*m2**2*xc01t - 4*sqrtS**2*xc01t + 5*m1**2*xc02t - 12*M1Lpt**2*xc02t - 3*m2**2*xc02t + 
     &                     3*m3**2*xc02t - m4**2*xc02t + 2*sqrtS**2*xc02t + m1**2*xc03t + 9*m2**2*xc03t + 3*m3**2*xc03t - 12*M3Lpt**2*xc03t + 7*m4**2*xc03t - 2*sqrtS**2*xc03t + 
     &                     4*m1**2*xc04t - 12*M1Lpt**2*xc04t + 4*m2**2*xc04t + 8*m4**2*xc04t - 6*sqrtS**2*xc04t + 4*m3**2*xc111t - 4*m4**2*xc111t + 4*sqrtS**2*xc111t + 4*m2**2*xc11t - 
     &                     16*m3**2*xc11t - 4*m4**2*xc11t + 16*m1**2*xc121t + 16*m2**2*xc121t + 24*m3**2*xc121t + 8*m4**2*xc121t + 2*m1**2*xc12t - 2*m2**2*xc12t + 2*sqrtS**2*xc12t + 
     &                     6*m2**2*xc13t + 4*m3**2*xc13t + 2*m4**2*xc13t + 8*m1**2*xc21t + 8*m2**2*xc21t + 16*m4**2*xc21t - 4*sqrtS**2*xc21t + 16*sqrtS**2*xc221t + 8*m1**2*xc22t - 
     &                     2*m2**2*xc22t - 2*m3**2*xc22t + 2*sqrtS**2*xc22t - 2*m1**2*xc23t + 2*m2**2*xc23t + 2*sqrtS**2*xc23t + 48*M0Lpt**2*xd00t - 24*m1**2*xd00t - 24*m2**2*xd00t + 
     &                     24*sqrtS**2*xd00t - 4*M0Lpt**2*m1**2*xd0t - m1**4*xd0t + m1**2*M1Lpt**2*xd0t - 4*M0Lpt**2*m2**2*xd0t - 2*m1**2*m2**2*xd0t + 9*M1Lpt**2*m2**2*xd0t + 3*m2**4*xd0t - 
     &                     12*M0Lpt**2*M2Lpt**2*xd0t + 2*m1**2*M2Lpt**2*xd0t + 2*m2**2*M2Lpt**2*xd0t + 6*M0Lpt**2*m3**2*xd0t - 3*m1**2*m3**2*xd0t + 3*M1Lpt**2*m3**2*xd0t - 
     &                     3*m2**2*m3**2*xd0t + 5*m1**2*M3Lpt**2*xd0t - 12*M1Lpt**2*M3Lpt**2*xd0t - 3*m2**2*M3Lpt**2*xd0t + 3*m3**2*M3Lpt**2*xd0t + 6*M0Lpt**2*m4**2*xd0t - 
     &                     7*m1**2*m4**2*xd0t + 7*M1Lpt**2*m4**2*xd0t + m2**2*m4**2*xd0t - M3Lpt**2*m4**2*xd0t - 4*M0Lpt**2*sqrtS**2*xd0t + 2*m1**2*sqrtS**2*xd0t - 2*M1Lpt**2*sqrtS**2*xd0t - 
     &                     2*m2**2*sqrtS**2*xd0t - 4*M2Lpt**2*sqrtS**2*xd0t + 2*M3Lpt**2*sqrtS**2*xd0t + 4*M0Lpt**2*m1**2*xd11t - 2*m1**4*xd11t - 4*M0Lpt**2*m2**2*xd11t + 2*m2**4*xd11t + 
     &                     4*M0Lpt**2*sqrtS**2*xd11t - 4*m2**2*sqrtS**2*xd11t + 2*sqrtS**4*xd11t - 4*M0Lpt**2*m1**2*xd12t + 2*m1**4*xd12t - 12*M0Lpt**2*m2**2*xd12t + 8*m1**2*m2**2*xd12t + 
     &                     6*m2**4*xd12t - 12*M0Lpt**2*m3**2*xd12t + 6*m1**2*m3**2*xd12t + 6*m2**2*m3**2*xd12t - 4*M0Lpt**2*m4**2*xd12t + 2*m1**2*m4**2*xd12t + 2*m2**2*m4**2*xd12t + 
     &                     8*M0Lpt**2*sqrtS**2*xd12t - 6*m1**2*sqrtS**2*xd12t - 10*m2**2*sqrtS**2*xd12t - 6*m3**2*sqrtS**2*xd12t - 2*m4**2*sqrtS**2*xd12t + 4*sqrtS**4*xd12t - 
     &                     8*M0Lpt**2*sqrtS**2*xd13t + 4*m1**2*sqrtS**2*xd13t + 4*m2**2*sqrtS**2*xd13t - 4*sqrtS**4*xd13t - 4*M0Lpt**2*m1**2*xd1t + m1**4*xd1t + 4*m1**2*M1Lpt**2*xd1t - 
     &                     12*M0Lpt**2*m2**2*xd1t + 2*m1**2*m2**2*xd1t + 4*M1Lpt**2*m2**2*xd1t + 5*m2**4*xd1t - 6*M0Lpt**2*m3**2*xd1t - m1**2*m3**2*xd1t - m2**2*m3**2*xd1t + 
     &                     2*m1**2*M3Lpt**2*xd1t - 2*m2**2*M3Lpt**2*xd1t - 2*M0Lpt**2*m4**2*xd1t + m1**2*m4**2*xd1t + m2**2*m4**2*xd1t + 6*M0Lpt**2*sqrtS**2*xd1t - m1**2*sqrtS**2*xd1t - 
     &                     4*M1Lpt**2*sqrtS**2*xd1t - 7*m2**2*sqrtS**2*xd1t - m3**2*sqrtS**2*xd1t + 2*M3Lpt**2*sqrtS**2*xd1t - m4**2*sqrtS**2*xd1t + 2*sqrtS**4*xd1t - 
     &                     8*M0Lpt**2*m1**2*xd22t + 4*m1**4*xd22t - 8*M0Lpt**2*m2**2*xd22t + 8*m1**2*m2**2*xd22t + 4*m2**4*xd22t - 8*M0Lpt**2*m3**2*xd22t + 4*m1**2*m3**2*xd22t + 
     &                     4*m2**2*m3**2*xd22t - 8*M0Lpt**2*m4**2*xd22t + 4*m1**2*m4**2*xd22t + 4*m2**2*m4**2*xd22t + 8*M0Lpt**2*sqrtS**2*xd22t - 8*m1**2*sqrtS**2*xd22t - 
     &                     8*m2**2*sqrtS**2*xd22t - 4*m3**2*sqrtS**2*xd22t - 4*m4**2*sqrtS**2*xd22t + 4*sqrtS**4*xd22t + 4*M0Lpt**2*m1**2*xd23t - 2*m1**4*xd23t + 12*M0Lpt**2*m2**2*xd23t - 
     &                     8*m1**2*m2**2*xd23t - 6*m2**4*xd23t + 12*M0Lpt**2*m3**2*xd23t - 6*m1**2*m3**2*xd23t - 6*m2**2*m3**2*xd23t + 4*M0Lpt**2*m4**2*xd23t - 2*m1**2*m4**2*xd23t - 
     &                     2*m2**2*m4**2*xd23t - 8*M0Lpt**2*sqrtS**2*xd23t + 6*m1**2*sqrtS**2*xd23t + 10*m2**2*sqrtS**2*xd23t + 6*m3**2*sqrtS**2*xd23t + 2*m4**2*sqrtS**2*xd23t - 
     &                     4*sqrtS**4*xd23t - 16*M0Lpt**2*m1**2*xd2t + 3*m1**4*xd2t - 12*M0Lpt**2*m2**2*xd2t + 6*m1**2*m2**2*xd2t + 6*M1Lpt**2*m2**2*xd2t + 7*m2**4*xd2t - 
     &                     10*M0Lpt**2*m3**2*xd2t + 2*m1**2*m3**2*xd2t + 4*M1Lpt**2*m3**2*xd2t + 4*m2**2*m3**2*xd2t - m3**4*xd2t + 8*m1**2*M3Lpt**2*xd2t - 2*m2**2*M3Lpt**2*xd2t - 
     &                     2*m3**2*M3Lpt**2*xd2t - 10*M0Lpt**2*m4**2*xd2t - 6*m1**2*m4**2*xd2t + 2*M1Lpt**2*m4**2*xd2t + 4*m2**2*m4**2*xd2t - 2*m3**2*m4**2*xd2t - m4**4*xd2t + 
     &                     6*M0Lpt**2*sqrtS**2*xd2t - 6*m1**2*sqrtS**2*xd2t - 10*m2**2*sqrtS**2*xd2t - 6*m3**2*sqrtS**2*xd2t + 2*M3Lpt**2*sqrtS**2*xd2t - 4*m4**2*sqrtS**2*xd2t + 
     &                     4*sqrtS**4*xd2t - 4*M0Lpt**2*m1**2*xd33t + 2*m1**4*xd33t + 4*M0Lpt**2*m2**2*xd33t - 2*m2**4*xd33t + 4*M0Lpt**2*sqrtS**2*xd33t - 4*m1**2*sqrtS**2*xd33t + 
     &                     2*sqrtS**4*xd33t - 4*M0Lpt**2*m1**2*xd3t + m1**4*xd3t - 2*m1**2*M1Lpt**2*xd3t + 4*M0Lpt**2*m2**2*xd3t - 6*m1**2*m2**2*xd3t + 2*M1Lpt**2*m2**2*xd3t - 3*m2**4*xd3t + 
     &                     6*M0Lpt**2*m3**2*xd3t - 3*m1**2*m3**2*xd3t - 3*m2**2*m3**2*xd3t + 4*m1**2*M3Lpt**2*xd3t + 4*m2**2*M3Lpt**2*xd3t + 2*M0Lpt**2*m4**2*xd3t - 5*m1**2*m4**2*xd3t - 
     &                     5*m2**2*m4**2*xd3t - 2*M0Lpt**2*sqrtS**2*xd3t + m1**2*sqrtS**2*xd3t + 2*M1Lpt**2*sqrtS**2*xd3t + 7*m2**2*sqrtS**2*xd3t + 3*m3**2*sqrtS**2*xd3t - 
     &                     4*M3Lpt**2*sqrtS**2*xd3t + 3*m4**2*sqrtS**2*xd3t - 2*sqrtS**4*xd3t) + 
     &                  tman*(12*xb01t + 12*xb02t - 48*xc001t + 12*M0Lpt**2*xc01t - 10*m1**2*xc01t - 10*m2**2*xc01t + 8*sqrtS**2*xc01t - 5*m1**2*xc02t + 12*M1Lpt**2*xc02t + 3*m2**2*xc02t - 
     &                     3*m3**2*xc02t + m4**2*xc02t - 2*sqrtS**2*xc02t - 4*uman*xc02t - m1**2*xc03t - 9*m2**2*xc03t - 3*m3**2*xc03t + 12*M3Lpt**2*xc03t - 7*m4**2*xc03t + 
     &                     2*sqrtS**2*xc03t - 4*uman*xc03t - 4*m1**2*xc04t + 12*M1Lpt**2*xc04t - 4*m2**2*xc04t - 8*m4**2*xc04t + 6*sqrtS**2*xc04t - 4*m3**2*xc111t + 4*m4**2*xc111t - 
     &                     4*sqrtS**2*xc111t - 4*m1**2*xc11t + 16*m3**2*xc11t + 4*m4**2*xc11t - 16*m1**2*xc121t - 16*m2**2*xc121t - 24*m3**2*xc121t - 8*m4**2*xc121t - 2*m1**2*xc12t + 
     &                     2*m2**2*xc12t - 2*sqrtS**2*xc12t + 2*m1**2*xc13t - 8*m2**2*xc13t - 4*m3**2*xc13t - 2*m4**2*xc13t + 2*uman*xc13t - 8*m1**2*xc21t - 8*m2**2*xc21t - 16*m4**2*xc21t + 
     &                     4*sqrtS**2*xc21t - 16*sqrtS**2*xc221t + 2*m1**2*xc22t + 8*m2**2*xc22t + 2*m3**2*xc22t - 6*sqrtS**2*xc22t + 2*uman*xc22t + 2*m1**2*xc23t - 2*m2**2*xc23t - 
     &                     2*sqrtS**2*xc23t - 48*M0Lpt**2*xd00t + 24*m1**2*xd00t + 24*m2**2*xd00t - 24*sqrtS**2*xd00t - 4*M0Lpt**2*m1**2*xd0t + 3*m1**4*xd0t - m1**2*M1Lpt**2*xd0t - 
     &                     4*M0Lpt**2*m2**2*xd0t - 18*m1**2*m2**2*xd0t - 9*M1Lpt**2*m2**2*xd0t - m2**4*xd0t + 12*M0Lpt**2*M2Lpt**2*xd0t - 10*m1**2*M2Lpt**2*xd0t - 10*m2**2*M2Lpt**2*xd0t - 
     &                     6*M0Lpt**2*m3**2*xd0t - m1**2*m3**2*xd0t - 3*M1Lpt**2*m3**2*xd0t - m2**2*m3**2*xd0t - 5*m1**2*M3Lpt**2*xd0t + 12*M1Lpt**2*M3Lpt**2*xd0t + 3*m2**2*M3Lpt**2*xd0t - 
     &                     3*m3**2*M3Lpt**2*xd0t - 6*M0Lpt**2*m4**2*xd0t + 3*m1**2*m4**2*xd0t - 7*M1Lpt**2*m4**2*xd0t - 5*m2**2*m4**2*xd0t + M3Lpt**2*m4**2*xd0t + 8*M0Lpt**2*sqrtS**2*xd0t - 
     &                     8*m1**2*sqrtS**2*xd0t + 2*M1Lpt**2*sqrtS**2*xd0t - 4*m2**2*sqrtS**2*xd0t + 8*M2Lpt**2*sqrtS**2*xd0t + 2*m3**2*sqrtS**2*xd0t - 2*M3Lpt**2*sqrtS**2*xd0t + 
     &                     2*m4**2*sqrtS**2*xd0t + 4*sqrtS**4*xd0t + 4*M0Lpt**2*uman*xd0t + 2*m1**2*uman*xd0t - 4*M1Lpt**2*uman*xd0t + 2*m2**2*uman*xd0t - 4*M3Lpt**2*uman*xd0t + 
     &                     4*sqrtS**2*uman*xd0t - 4*M0Lpt**2*m1**2*xd11t + 2*m1**4*xd11t + 4*M0Lpt**2*m2**2*xd11t - 2*m2**4*xd11t - 4*M0Lpt**2*sqrtS**2*xd11t + 4*m2**2*sqrtS**2*xd11t - 
     &                     2*sqrtS**4*xd11t + 4*M0Lpt**2*m1**2*xd12t - 2*m1**4*xd12t + 12*M0Lpt**2*m2**2*xd12t - 8*m1**2*m2**2*xd12t - 6*m2**4*xd12t + 12*M0Lpt**2*m3**2*xd12t - 
     &                     6*m1**2*m3**2*xd12t - 6*m2**2*m3**2*xd12t + 4*M0Lpt**2*m4**2*xd12t - 2*m1**2*m4**2*xd12t - 2*m2**2*m4**2*xd12t - 8*M0Lpt**2*sqrtS**2*xd12t + 
     &                     6*m1**2*sqrtS**2*xd12t + 10*m2**2*sqrtS**2*xd12t + 6*m3**2*sqrtS**2*xd12t + 2*m4**2*sqrtS**2*xd12t - 4*sqrtS**4*xd12t + 8*M0Lpt**2*uman*xd12t - 
     &                     4*m1**2*uman*xd12t - 4*m2**2*uman*xd12t + 4*sqrtS**2*uman*xd12t + 8*M0Lpt**2*sqrtS**2*xd13t - 4*m1**2*sqrtS**2*xd13t - 4*m2**2*sqrtS**2*xd13t + 4*sqrtS**4*xd13t - 
     &                     4*M0Lpt**2*m1**2*xd1t + m1**4*xd1t + 4*m1**2*M1Lpt**2*xd1t + 4*M0Lpt**2*m2**2*xd1t - 22*m1**2*m2**2*xd1t + 4*M1Lpt**2*m2**2*xd1t - 3*m2**4*xd1t + 
     &                     6*M0Lpt**2*m3**2*xd1t - 11*m1**2*m3**2*xd1t - 11*m2**2*m3**2*xd1t - 2*m1**2*M3Lpt**2*xd1t + 2*m2**2*M3Lpt**2*xd1t + 2*M0Lpt**2*m4**2*xd1t - 5*m1**2*m4**2*xd1t - 
     &                     5*m2**2*m4**2*xd1t - 2*M0Lpt**2*sqrtS**2*xd1t - 5*m1**2*sqrtS**2*xd1t + m2**2*sqrtS**2*xd1t + 7*m3**2*sqrtS**2*xd1t - 2*M3Lpt**2*sqrtS**2*xd1t + 
     &                     3*m4**2*sqrtS**2*xd1t + 2*sqrtS**4*xd1t + 4*M0Lpt**2*uman*xd1t + 2*m1**2*uman*xd1t + 2*m2**2*uman*xd1t + 8*M0Lpt**2*m1**2*xd22t - 4*m1**4*xd22t + 
     &                     8*M0Lpt**2*m2**2*xd22t - 8*m1**2*m2**2*xd22t - 4*m2**4*xd22t + 8*M0Lpt**2*m3**2*xd22t - 4*m1**2*m3**2*xd22t - 4*m2**2*m3**2*xd22t + 8*M0Lpt**2*m4**2*xd22t - 
     &                     4*m1**2*m4**2*xd22t - 4*m2**2*m4**2*xd22t - 8*M0Lpt**2*sqrtS**2*xd22t + 8*m1**2*sqrtS**2*xd22t + 8*m2**2*sqrtS**2*xd22t + 4*m3**2*sqrtS**2*xd22t + 
     &                     4*m4**2*sqrtS**2*xd22t - 4*sqrtS**4*xd22t + 8*M0Lpt**2*uman*xd22t - 4*m1**2*uman*xd22t - 4*m2**2*uman*xd22t + 4*sqrtS**2*uman*xd22t - 4*M0Lpt**2*m1**2*xd23t + 
     &                     2*m1**4*xd23t - 12*M0Lpt**2*m2**2*xd23t + 8*m1**2*m2**2*xd23t + 6*m2**4*xd23t - 12*M0Lpt**2*m3**2*xd23t + 6*m1**2*m3**2*xd23t + 6*m2**2*m3**2*xd23t - 
     &                     4*M0Lpt**2*m4**2*xd23t + 2*m1**2*m4**2*xd23t + 2*m2**2*m4**2*xd23t + 8*M0Lpt**2*sqrtS**2*xd23t - 6*m1**2*sqrtS**2*xd23t - 10*m2**2*sqrtS**2*xd23t - 
     &                     6*m3**2*sqrtS**2*xd23t - 2*m4**2*sqrtS**2*xd23t + 4*sqrtS**4*xd23t + 8*M0Lpt**2*uman*xd23t - 4*m1**2*uman*xd23t - 4*m2**2*uman*xd23t + 4*sqrtS**2*uman*xd23t + 
     &                     4*M0Lpt**2*m1**2*xd2t - 3*m1**4*xd2t + 2*m1**2*M1Lpt**2*xd2t + 8*M0Lpt**2*m2**2*xd2t - 26*m1**2*m2**2*xd2t - 8*M1Lpt**2*m2**2*xd2t - 3*m2**4*xd2t + 
     &                     10*M0Lpt**2*m3**2*xd2t - 16*m1**2*m3**2*xd2t - 4*M1Lpt**2*m3**2*xd2t + 2*m2**2*m3**2*xd2t + m3**4*xd2t + 2*m1**2*M3Lpt**2*xd2t + 8*m2**2*M3Lpt**2*xd2t + 
     &                     2*m3**2*M3Lpt**2*xd2t + 10*M0Lpt**2*m4**2*xd2t - 4*m1**2*m4**2*xd2t - 2*M1Lpt**2*m4**2*xd2t - 18*m2**2*m4**2*xd2t + 2*m3**2*m4**2*xd2t + m4**4*xd2t - 
     &                     2*M0Lpt**2*sqrtS**2*xd2t + 2*m1**2*sqrtS**2*xd2t + 2*m2**2*sqrtS**2*xd2t + 8*m3**2*sqrtS**2*xd2t - 6*M3Lpt**2*sqrtS**2*xd2t + 10*m4**2*sqrtS**2*xd2t + 
     &                     12*M0Lpt**2*uman*xd2t + 2*m1**2*uman*xd2t + 2*M1Lpt**2*uman*xd2t - 2*m2**2*uman*xd2t + 2*m3**2*uman*xd2t + 2*M3Lpt**2*uman*xd2t + 2*m4**2*uman*xd2t + 
     &                     6*sqrtS**2*uman*xd2t - 2*uman**2*xd2t + 4*M0Lpt**2*m1**2*xd33t - 2*m1**4*xd33t - 4*M0Lpt**2*m2**2*xd33t + 2*m2**4*xd33t - 4*M0Lpt**2*sqrtS**2*xd33t + 
     &                     4*m1**2*sqrtS**2*xd33t - 2*sqrtS**4*xd33t - 4*M0Lpt**2*m1**2*xd3t + m1**4*xd3t + 2*m1**2*M1Lpt**2*xd3t - 12*M0Lpt**2*m2**2*xd3t - 14*m1**2*m2**2*xd3t - 
     &                     2*M1Lpt**2*m2**2*xd3t + 5*m2**4*xd3t - 6*M0Lpt**2*m3**2*xd3t - m1**2*m3**2*xd3t - m2**2*m3**2*xd3t + 4*m1**2*M3Lpt**2*xd3t + 4*m2**2*M3Lpt**2*xd3t - 
     &                     2*M0Lpt**2*m4**2*xd3t - 7*m1**2*m4**2*xd3t - 7*m2**2*m4**2*xd3t + 6*M0Lpt**2*sqrtS**2*xd3t - 7*m1**2*sqrtS**2*xd3t - 2*M1Lpt**2*sqrtS**2*xd3t - 
     &                     13*m2**2*sqrtS**2*xd3t - m3**2*sqrtS**2*xd3t + 3*m4**2*sqrtS**2*xd3t + 6*sqrtS**4*xd3t + 4*M0Lpt**2*uman*xd3t + 2*m1**2*uman*xd3t + 2*m2**2*uman*xd3t)))))))/
     &  (64.*pi**2*propac)  )

      Box_UV = Box_UV + 2*colorfactorSV*dreal( -(Coupl_Ac*Coupl_Bc*(g0Ru*(2*g1Ru*M1Lpu*(g2Lu*M2Lpu*(12*g3Ru*M0Lpu*M3Lpu*(tman - uman)*xd0u + 
     &                g3Lu*(-2*m1**4*xd1u + 20*m1**2*m2**2*xd1u - 2*m2**4*xd1u + 4*m1**2*m3**2*xd1u + 4*m2**2*m3**2*xd1u + 4*m1**2*m4**2*xd1u + 4*m2**2*m4**2*xd1u + 6*m1**2*sqrtS**2*xd1u + 
     &                   6*m2**2*sqrtS**2*xd1u - 2*m3**2*sqrtS**2*xd1u - 2*m4**2*sqrtS**2*xd1u - 4*sqrtS**4*xd1u - 2*m1**4*xd2u + 6*m1**2*m2**2*xd2u - 4*m2**4*xd2u + 4*m1**2*m3**2*xd2u + 
     &                   8*m2**2*m3**2*xd2u - 2*m1**2*m4**2*xd2u - 10*m2**2*m4**2*xd2u + 5*m1**2*sqrtS**2*xd2u + 7*m2**2*sqrtS**2*xd2u - 3*m3**2*sqrtS**2*xd2u + 3*m4**2*sqrtS**2*xd2u - 
     &                   3*sqrtS**4*xd2u + tman**2*xd2u - uman**2*xd2u + tman*
     &                    (6*xc02u + 6*xc04u + 6*M0Lpu**2*xd0u - 6*m2**2*xd0u + 6*M3Lpu**2*xd0u - 4*m1**2*xd1u - 4*m2**2*xd1u + 2*sqrtS**2*xd1u + m1**2*xd2u - 3*m2**2*xd2u - m3**2*xd2u - 
     &                      m4**2*xd2u - 2*sqrtS**2*xd2u) + uman*(-6*xc02u - 6*xc04u - 6*M0Lpu**2*xd0u + 6*m2**2*xd0u - 6*M3Lpu**2*xd0u - 4*m1**2*xd1u - 4*m2**2*xd1u + 2*sqrtS**2*xd1u - 
     &                      3*m1**2*xd2u + 5*m2**2*xd2u + m3**2*xd2u + m4**2*xd2u + 2*sqrtS**2*xd2u))) + 
     &             g2Ru*(g3Ru*M0Lpu*(4*m2**2*m3**2*xd0u + 4*m1**2*m4**2*xd0u + 8*m2**2*m4**2*xd0u - m1**2*sqrtS**2*xd0u - 3*m2**2*sqrtS**2*xd0u - m3**2*sqrtS**2*xd0u - 
     &                   3*m4**2*sqrtS**2*xd0u + sqrtS**4*xd0u + 6*m1**2*m3**2*xd11u + 2*m2**2*m3**2*xd11u - 6*m1**2*m4**2*xd11u - 2*m2**2*m4**2*xd11u - 2*m3**2*sqrtS**2*xd11u + 
     &                   2*m4**2*sqrtS**2*xd11u + 16*m1**2*m3**2*xd12u + 8*m2**2*m3**2*xd12u + 8*m1**2*m4**2*xd12u - 6*m1**2*sqrtS**2*xd12u - 2*m2**2*sqrtS**2*xd12u - 
     &                   6*m3**2*sqrtS**2*xd12u - 2*m4**2*sqrtS**2*xd12u + 2*sqrtS**4*xd12u - 4*m1**2*m3**2*xd13u + 4*m2**2*m3**2*xd13u + 4*m1**2*m4**2*xd13u - 4*m2**2*m4**2*xd13u + 
     &                   6*m1**2*m3**2*xd1u + 6*m2**2*m3**2*xd1u + 6*m1**2*m4**2*xd1u - 2*m2**2*m4**2*xd1u - 3*m1**2*sqrtS**2*xd1u - m2**2*sqrtS**2*xd1u - 3*m3**2*sqrtS**2*xd1u - 
     &                   m4**2*sqrtS**2*xd1u + sqrtS**4*xd1u + 8*m1**2*m3**2*xd22u + 8*m2**2*m3**2*xd22u + 8*m1**2*m4**2*xd22u + 8*m2**2*m4**2*xd22u - 4*m1**2*sqrtS**2*xd22u - 
     &                   4*m2**2*sqrtS**2*xd22u - 4*m3**2*sqrtS**2*xd22u - 4*m4**2*sqrtS**2*xd22u + 2*sqrtS**4*xd22u + 8*m2**2*m3**2*xd23u + 8*m1**2*m4**2*xd23u + 16*m2**2*m4**2*xd23u - 
     &                   2*m1**2*sqrtS**2*xd23u - 6*m2**2*sqrtS**2*xd23u - 2*m3**2*sqrtS**2*xd23u - 6*m4**2*sqrtS**2*xd23u + 2*sqrtS**4*xd23u + 8*m1**2*m3**2*xd2u + 12*m2**2*m3**2*xd2u + 
     &                   12*m1**2*m4**2*xd2u + 16*m2**2*m4**2*xd2u - 5*m1**2*sqrtS**2*xd2u - 7*m2**2*sqrtS**2*xd2u - 5*m3**2*sqrtS**2*xd2u - 7*m4**2*sqrtS**2*xd2u + 3*sqrtS**4*xd2u - 
     &                   2*m1**2*m3**2*xd33u - 6*m2**2*m3**2*xd33u + 2*m1**2*m4**2*xd33u + 6*m2**2*m4**2*xd33u + 2*m3**2*sqrtS**2*xd33u - 2*m4**2*sqrtS**2*xd33u + 
     &                   tman**2*(-3*xd0u + 2*xd12u + xd1u + 2*xd22u - 6*xd23u - xd2u - 3*xd3u) - 2*m1**2*m3**2*xd3u - 2*m2**2*m3**2*xd3u + 6*m1**2*m4**2*xd3u + 14*m2**2*m4**2*xd3u - 
     &                   m1**2*sqrtS**2*xd3u - 3*m2**2*sqrtS**2*xd3u + m3**2*sqrtS**2*xd3u - 5*m4**2*sqrtS**2*xd3u + sqrtS**4*xd3u + 
     &                   uman**2*(xd0u - 6*xd12u - 3*xd1u - 6*xd22u + 2*xd23u - 5*xd2u + xd3u) + 
     &                   uman*(6*xc01u + 6*xc02u - 24*xd00u - m1**2*xd0u - 3*m2**2*xd0u + 6*M2Lpu**2*xd0u - 7*m3**2*xd0u + 6*M3Lpu**2*xd0u - 3*m4**2*xd0u + 2*sqrtS**2*xd0u - 2*m1**2*xd11u + 
     &                      2*m2**2*xd11u - 2*sqrtS**2*xd11u + 2*m1**2*xd12u + 6*m2**2*xd12u + 2*m3**2*xd12u + 6*m4**2*xd12u - 4*sqrtS**2*xd12u + 4*sqrtS**2*xd13u + m1**2*xd1u + 
     &                      3*m2**2*xd1u + m3**2*xd1u + 3*m4**2*xd1u + 4*m1**2*xd22u + 4*m2**2*xd22u + 4*m3**2*xd22u + 4*m4**2*xd22u - 4*sqrtS**2*xd22u - 2*m1**2*xd23u - 6*m2**2*xd23u - 
     &                      2*m3**2*xd23u - 6*m4**2*xd23u + 4*sqrtS**2*xd23u + 3*m1**2*xd2u + m2**2*xd2u + 3*m3**2*xd2u + m4**2*xd2u - 2*sqrtS**2*xd2u + 2*m1**2*xd33u - 2*m2**2*xd33u - 
     &                      2*sqrtS**2*xd33u + m1**2*xd3u - 5*m2**2*xd3u - m3**2*xd3u - 3*m4**2*xd3u) - 
     &                   tman*(6*xc01u + 6*xc02u - 24*xd00u - m1**2*xd0u - 3*m2**2*xd0u + 6*M2Lpu**2*xd0u - 7*m3**2*xd0u + 6*M3Lpu**2*xd0u - 3*m4**2*xd0u + 2*sqrtS**2*xd0u - 2*uman*xd0u - 
     &                      2*m1**2*xd11u + 2*m2**2*xd11u - 2*sqrtS**2*xd11u + 2*m1**2*xd12u + 6*m2**2*xd12u + 2*m3**2*xd12u + 6*m4**2*xd12u - 4*sqrtS**2*xd12u - 4*uman*xd12u + 
     &                      4*sqrtS**2*xd13u + m1**2*xd1u + 3*m2**2*xd1u + m3**2*xd1u + 3*m4**2*xd1u - 2*uman*xd1u + 4*m1**2*xd22u + 4*m2**2*xd22u + 4*m3**2*xd22u + 4*m4**2*xd22u - 
     &                      4*sqrtS**2*xd22u - 4*uman*xd22u - 2*m1**2*xd23u - 6*m2**2*xd23u - 2*m3**2*xd23u - 6*m4**2*xd23u + 4*sqrtS**2*xd23u - 4*uman*xd23u + 3*m1**2*xd2u + m2**2*xd2u + 
     &                      3*m3**2*xd2u + m4**2*xd2u - 2*sqrtS**2*xd2u - 6*uman*xd2u + 2*m1**2*xd33u - 2*m2**2*xd33u - 2*sqrtS**2*xd33u + m1**2*xd3u - 5*m2**2*xd3u - m3**2*xd3u - 
     &                      3*m4**2*xd3u - 2*uman*xd3u)) + g3Lu*M3Lpu*(6*m1**2*m3**2*xd11u + 2*m2**2*m3**2*xd11u - 6*m1**2*m4**2*xd11u - 2*m2**2*m4**2*xd11u - 2*m3**2*sqrtS**2*xd11u + 
     &                   2*m4**2*sqrtS**2*xd11u + 16*m1**2*m3**2*xd12u + 8*m2**2*m3**2*xd12u + 8*m1**2*m4**2*xd12u - 6*m1**2*sqrtS**2*xd12u - 2*m2**2*sqrtS**2*xd12u - 
     &                   6*m3**2*sqrtS**2*xd12u - 2*m4**2*sqrtS**2*xd12u + 2*sqrtS**4*xd12u - 4*m1**2*m3**2*xd13u + 4*m2**2*m3**2*xd13u + 4*m1**2*m4**2*xd13u - 4*m2**2*m4**2*xd13u + 
     &                   8*m1**2*m3**2*xd1u + 4*m2**2*m3**2*xd1u + 4*m1**2*m4**2*xd1u - 3*m1**2*sqrtS**2*xd1u - m2**2*sqrtS**2*xd1u - 3*m3**2*sqrtS**2*xd1u - m4**2*sqrtS**2*xd1u + 
     &                   sqrtS**4*xd1u + 8*m1**2*m3**2*xd22u + 8*m2**2*m3**2*xd22u + 8*m1**2*m4**2*xd22u + 8*m2**2*m4**2*xd22u - 4*m1**2*sqrtS**2*xd22u - 4*m2**2*sqrtS**2*xd22u - 
     &                   4*m3**2*sqrtS**2*xd22u - 4*m4**2*sqrtS**2*xd22u + 2*sqrtS**4*xd22u + 8*m2**2*m3**2*xd23u + 8*m1**2*m4**2*xd23u + 16*m2**2*m4**2*xd23u - 2*m1**2*sqrtS**2*xd23u - 
     &                   6*m2**2*sqrtS**2*xd23u - 2*m3**2*sqrtS**2*xd23u - 6*m4**2*sqrtS**2*xd23u + 2*sqrtS**4*xd23u + 8*m1**2*m3**2*xd2u + 8*m2**2*m3**2*xd2u + 8*m1**2*m4**2*xd2u + 
     &                   8*m2**2*m4**2*xd2u - 4*m1**2*sqrtS**2*xd2u - 4*m2**2*sqrtS**2*xd2u - 4*m3**2*sqrtS**2*xd2u - 4*m4**2*sqrtS**2*xd2u + 2*sqrtS**4*xd2u - 2*m1**2*m3**2*xd33u - 
     &                   6*m2**2*m3**2*xd33u + 2*m1**2*m4**2*xd33u + 6*m2**2*m4**2*xd33u + 2*m3**2*sqrtS**2*xd33u - 2*m4**2*sqrtS**2*xd33u + 
     &                   tman**2*(2*xd12u + xd1u + 2*xd22u - 6*xd23u + 2*xd2u - 3*xd3u) + 4*m2**2*m3**2*xd3u + 4*m1**2*m4**2*xd3u + 8*m2**2*m4**2*xd3u - m1**2*sqrtS**2*xd3u - 
     &                   3*m2**2*sqrtS**2*xd3u - m3**2*sqrtS**2*xd3u - 3*m4**2*sqrtS**2*xd3u + sqrtS**4*xd3u + uman**2*(-6*xd0u - 6*xd12u - 3*xd1u - 6*xd22u + 2*xd23u - 6*xd2u + xd3u) + 
     &                   uman*(6*xc01u + 6*xc04u - 24*xd00u + 6*M0Lpu**2*xd0u + 6*M2Lpu**2*xd0u - 2*m1**2*xd11u + 2*m2**2*xd11u - 2*sqrtS**2*xd11u + 2*m1**2*xd12u + 6*m2**2*xd12u + 
     &                      2*m3**2*xd12u + 6*m4**2*xd12u - 4*sqrtS**2*xd12u + 4*sqrtS**2*xd13u + m1**2*xd1u + 3*m2**2*xd1u + m3**2*xd1u + 3*m4**2*xd1u - 2*sqrtS**2*xd1u + 4*m1**2*xd22u + 
     &                      4*m2**2*xd22u + 4*m3**2*xd22u + 4*m4**2*xd22u - 4*sqrtS**2*xd22u - 2*m1**2*xd23u - 6*m2**2*xd23u - 2*m3**2*xd23u - 6*m4**2*xd23u + 4*sqrtS**2*xd23u + 
     &                      4*m1**2*xd2u + 4*m2**2*xd2u + 4*m3**2*xd2u + 4*m4**2*xd2u - 4*sqrtS**2*xd2u + 2*m1**2*xd33u - 2*m2**2*xd33u - 2*sqrtS**2*xd33u - m1**2*xd3u - 3*m2**2*xd3u - 
     &                      m3**2*xd3u - 3*m4**2*xd3u + 2*sqrtS**2*xd3u) + tman*
     &                    (-6*xc01u - 6*xc04u + 24*xd00u - 6*M0Lpu**2*xd0u - 6*M2Lpu**2*xd0u + 6*uman*xd0u + 2*m1**2*xd11u - 2*m2**2*xd11u + 2*sqrtS**2*xd11u - 2*m1**2*xd12u - 
     &                      6*m2**2*xd12u - 2*m3**2*xd12u - 6*m4**2*xd12u + 4*sqrtS**2*xd12u + 4*uman*xd12u - 4*sqrtS**2*xd13u - m1**2*xd1u - 3*m2**2*xd1u - m3**2*xd1u - 3*m4**2*xd1u + 
     &                      2*sqrtS**2*xd1u + 2*uman*xd1u - 4*m1**2*xd22u - 4*m2**2*xd22u - 4*m3**2*xd22u - 4*m4**2*xd22u + 4*sqrtS**2*xd22u + 4*uman*xd22u + 2*m1**2*xd23u + 6*m2**2*xd23u + 
     &                      2*m3**2*xd23u + 6*m4**2*xd23u - 4*sqrtS**2*xd23u + 4*uman*xd23u - 4*m1**2*xd2u - 4*m2**2*xd2u - 4*m3**2*xd2u - 4*m4**2*xd2u + 4*sqrtS**2*xd2u + 4*uman*xd2u - 
     &                      2*m1**2*xd33u + 2*m2**2*xd33u + 2*sqrtS**2*xd33u + m1**2*xd3u + 3*m2**2*xd3u + m3**2*xd3u + 3*m4**2*xd3u - 2*sqrtS**2*xd3u + 2*uman*xd3u)))) + 
     &          g1Lu*(2*g2Ru*M2Lpu*(-2*g3Ru*M0Lpu*(m1**4*xd0u - 10*m1**2*m2**2*xd0u + m2**4*xd0u - 2*m1**2*m3**2*xd0u - 2*m2**2*m3**2*xd0u - 2*m1**2*m4**2*xd0u - 2*m2**2*m4**2*xd0u - 
     &                   3*m1**2*sqrtS**2*xd0u - 3*m2**2*sqrtS**2*xd0u + m3**2*sqrtS**2*xd0u + m4**2*sqrtS**2*xd0u + 2*sqrtS**4*xd0u + m1**4*xd1u - 10*m1**2*m2**2*xd1u + m2**4*xd1u - 
     &                   2*m1**2*m3**2*xd1u - 2*m2**2*m3**2*xd1u - 2*m1**2*m4**2*xd1u - 2*m2**2*m4**2*xd1u - 3*m1**2*sqrtS**2*xd1u - 3*m2**2*sqrtS**2*xd1u + m3**2*sqrtS**2*xd1u + 
     &                   m4**2*sqrtS**2*xd1u + 2*sqrtS**4*xd1u - 10*m1**2*m2**2*xd2u + 2*m2**4*xd2u - 4*m2**2*m3**2*xd2u - 6*m1**2*m4**2*xd2u + 2*m2**2*m4**2*xd2u - 2*m1**2*sqrtS**2*xd2u - 
     &                   4*m2**2*sqrtS**2*xd2u + m3**2*sqrtS**2*xd2u + m4**2*sqrtS**2*xd2u + 2*sqrtS**4*xd2u - tman**2*xd2u + uman**2*xd2u + m1**4*xd3u - 10*m1**2*m2**2*xd3u + m2**4*xd3u - 
     &                   2*m1**2*m3**2*xd3u - 2*m2**2*m3**2*xd3u - 2*m1**2*m4**2*xd3u - 2*m2**2*m4**2*xd3u - 3*m1**2*sqrtS**2*xd3u - 3*m2**2*sqrtS**2*xd3u + m3**2*sqrtS**2*xd3u + 
     &                   m4**2*sqrtS**2*xd3u + 2*sqrtS**4*xd3u + uman*(3*xc02u + 3*xc03u + 2*m1**2*xd0u + 3*M1Lpu**2*xd0u + 2*m2**2*xd0u + 3*M3Lpu**2*xd0u - 4*sqrtS**2*xd0u + 2*m1**2*xd1u + 
     &                      2*m2**2*xd1u - sqrtS**2*xd1u + 2*m1**2*xd2u - m3**2*xd2u - m4**2*xd2u - 3*sqrtS**2*xd2u + 2*m1**2*xd3u + 2*m2**2*xd3u - sqrtS**2*xd3u) + 
     &                   tman*(-3*xc02u - 3*xc03u + 2*m1**2*xd0u - 3*M1Lpu**2*xd0u + 2*m2**2*xd0u - 3*M3Lpu**2*xd0u + 2*sqrtS**2*xd0u + 2*m1**2*xd1u + 2*m2**2*xd1u - sqrtS**2*xd1u + 
     &                      4*m1**2*xd2u + 2*m2**2*xd2u + m3**2*xd2u + m4**2*xd2u + sqrtS**2*xd2u + 2*m1**2*xd3u + 2*m2**2*xd3u - sqrtS**2*xd3u)) + 
     &                g3Lu*M3Lpu*(2*m1**4*xd2u + 14*m1**2*m2**2*xd2u - 4*m1**2*m3**2*xd2u + 14*m1**2*m4**2*xd2u + 6*m2**2*m4**2*xd2u - m1**2*sqrtS**2*xd2u + m2**2*sqrtS**2*xd2u + 
     &                   m3**2*sqrtS**2*xd2u - 5*m4**2*sqrtS**2*xd2u - sqrtS**4*xd2u + tman**2*xd2u - uman**2*xd2u - 2*m1**4*xd3u + 20*m1**2*m2**2*xd3u - 2*m2**4*xd3u + 4*m1**2*m3**2*xd3u + 
     &                   4*m2**2*m3**2*xd3u + 4*m1**2*m4**2*xd3u + 4*m2**2*m4**2*xd3u + 6*m1**2*sqrtS**2*xd3u + 6*m2**2*sqrtS**2*xd3u - 2*m3**2*sqrtS**2*xd3u - 2*m4**2*sqrtS**2*xd3u - 
     &                   4*sqrtS**4*xd3u - uman*(6*xc03u + 6*xc04u + 6*M0Lpu**2*xd0u - 6*m1**2*xd0u + 6*M1Lpu**2*xd0u + m1**2*xd2u + 5*m2**2*xd2u - m3**2*xd2u - m4**2*xd2u - 
     &                      4*sqrtS**2*xd2u + 4*m1**2*xd3u + 4*m2**2*xd3u - 2*sqrtS**2*xd3u) + 
     &                   tman*(6*xc03u + 6*xc04u + 6*M0Lpu**2*xd0u - 6*m1**2*xd0u + 6*M1Lpu**2*xd0u - 9*m1**2*xd2u - m2**2*xd2u - m3**2*xd2u - m4**2*xd2u - 4*m1**2*xd3u - 4*m2**2*xd3u + 
     &                      2*sqrtS**2*xd3u))) + g2Lu*(2*g3Ru*M0Lpu*M3Lpu*(8*m1**2*m3**2*xd0u + 4*m2**2*m3**2*xd0u + 4*m1**2*m4**2*xd0u - 3*m1**2*sqrtS**2*xd0u - m2**2*sqrtS**2*xd0u - 
     &                   3*m3**2*sqrtS**2*xd0u - m4**2*sqrtS**2*xd0u + sqrtS**4*xd0u + 6*m1**2*m3**2*xd11u + 2*m2**2*m3**2*xd11u - 6*m1**2*m4**2*xd11u - 2*m2**2*m4**2*xd11u - 
     &                   2*m3**2*sqrtS**2*xd11u + 2*m4**2*sqrtS**2*xd11u + 16*m1**2*m3**2*xd12u + 8*m2**2*m3**2*xd12u + 8*m1**2*m4**2*xd12u - 6*m1**2*sqrtS**2*xd12u - 
     &                   2*m2**2*sqrtS**2*xd12u - 6*m3**2*sqrtS**2*xd12u - 2*m4**2*sqrtS**2*xd12u + 2*sqrtS**4*xd12u - 4*m1**2*m3**2*xd13u + 4*m2**2*m3**2*xd13u + 4*m1**2*m4**2*xd13u - 
     &                   4*m2**2*m4**2*xd13u + 14*m1**2*m3**2*xd1u + 6*m2**2*m3**2*xd1u - 2*m1**2*m4**2*xd1u - 2*m2**2*m4**2*xd1u - 3*m1**2*sqrtS**2*xd1u - m2**2*sqrtS**2*xd1u - 
     &                   5*m3**2*sqrtS**2*xd1u + m4**2*sqrtS**2*xd1u + sqrtS**4*xd1u + 8*m1**2*m3**2*xd22u + 8*m2**2*m3**2*xd22u + 8*m1**2*m4**2*xd22u + 8*m2**2*m4**2*xd22u - 
     &                   4*m1**2*sqrtS**2*xd22u - 4*m2**2*sqrtS**2*xd22u - 4*m3**2*sqrtS**2*xd22u - 4*m4**2*sqrtS**2*xd22u + 2*sqrtS**4*xd22u + 8*m2**2*m3**2*xd23u + 8*m1**2*m4**2*xd23u + 
     &                   16*m2**2*m4**2*xd23u - 2*m1**2*sqrtS**2*xd23u - 6*m2**2*sqrtS**2*xd23u - 2*m3**2*sqrtS**2*xd23u - 6*m4**2*sqrtS**2*xd23u + 2*sqrtS**4*xd23u + 16*m1**2*m3**2*xd2u + 
     &                   12*m2**2*m3**2*xd2u + 12*m1**2*m4**2*xd2u + 8*m2**2*m4**2*xd2u - 7*m1**2*sqrtS**2*xd2u - 5*m2**2*sqrtS**2*xd2u - 7*m3**2*sqrtS**2*xd2u - 5*m4**2*sqrtS**2*xd2u + 
     &                   3*sqrtS**4*xd2u - 2*m1**2*m3**2*xd33u - 6*m2**2*m3**2*xd33u + 2*m1**2*m4**2*xd33u + 6*m2**2*m4**2*xd33u + 2*m3**2*sqrtS**2*xd33u - 2*m4**2*sqrtS**2*xd33u + 
     &                   tman**2*(xd0u + 2*xd12u + xd1u + 2*xd22u - 6*xd23u + 3*xd2u - 3*xd3u) - 2*m1**2*m3**2*xd3u + 6*m2**2*m3**2*xd3u + 6*m1**2*m4**2*xd3u + 6*m2**2*m4**2*xd3u - 
     &                   m1**2*sqrtS**2*xd3u - 3*m2**2*sqrtS**2*xd3u - m3**2*sqrtS**2*xd3u - 3*m4**2*sqrtS**2*xd3u + sqrtS**4*xd3u + 
     &                   uman**2*(-3*xd0u - 6*xd12u - 3*xd1u - 6*xd22u + 2*xd23u - 9*xd2u + xd3u) + 
     &                   uman*(6*xc01u + 6*xc03u - 24*xd00u + m1**2*xd0u + 6*M1Lpu**2*xd0u + 3*m2**2*xd0u + 6*M2Lpu**2*xd0u + m3**2*xd0u - 3*m4**2*xd0u - 2*sqrtS**2*xd0u - 2*m1**2*xd11u + 
     &                      2*m2**2*xd11u - 2*sqrtS**2*xd11u + 2*m1**2*xd12u + 6*m2**2*xd12u + 2*m3**2*xd12u + 6*m4**2*xd12u - 4*sqrtS**2*xd12u + 4*sqrtS**2*xd13u - m1**2*xd1u + 
     &                      5*m2**2*xd1u + m3**2*xd1u + 3*m4**2*xd1u - 4*sqrtS**2*xd1u + 4*m1**2*xd22u + 4*m2**2*xd22u + 4*m3**2*xd22u + 4*m4**2*xd22u - 4*sqrtS**2*xd22u - 2*m1**2*xd23u - 
     &                      6*m2**2*xd23u - 2*m3**2*xd23u - 6*m4**2*xd23u + 4*sqrtS**2*xd23u + 5*m1**2*xd2u + 7*m2**2*xd2u + 5*m3**2*xd2u + 7*m4**2*xd2u - 6*sqrtS**2*xd2u + 2*m1**2*xd33u - 
     &                      2*m2**2*xd33u - 2*sqrtS**2*xd33u - m1**2*xd3u - 3*m2**2*xd3u - m3**2*xd3u - 3*m4**2*xd3u + 4*sqrtS**2*xd3u) + 
     &                   tman*(-6*xc01u - 6*xc03u + 24*xd00u - m1**2*xd0u - 6*M1Lpu**2*xd0u - 3*m2**2*xd0u - 6*M2Lpu**2*xd0u - m3**2*xd0u + 3*m4**2*xd0u + 2*sqrtS**2*xd0u + 2*uman*xd0u + 
     &                      2*m1**2*xd11u - 2*m2**2*xd11u + 2*sqrtS**2*xd11u - 2*m1**2*xd12u - 6*m2**2*xd12u - 2*m3**2*xd12u - 6*m4**2*xd12u + 4*sqrtS**2*xd12u + 4*uman*xd12u - 
     &                      4*sqrtS**2*xd13u + m1**2*xd1u - 5*m2**2*xd1u - m3**2*xd1u - 3*m4**2*xd1u + 4*sqrtS**2*xd1u + 2*uman*xd1u - 4*m1**2*xd22u - 4*m2**2*xd22u - 4*m3**2*xd22u - 
     &                      4*m4**2*xd22u + 4*sqrtS**2*xd22u + 4*uman*xd22u + 2*m1**2*xd23u + 6*m2**2*xd23u + 2*m3**2*xd23u + 6*m4**2*xd23u - 4*sqrtS**2*xd23u + 4*uman*xd23u - 5*m1**2*xd2u - 
     &                      7*m2**2*xd2u - 5*m3**2*xd2u - 7*m4**2*xd2u + 6*sqrtS**2*xd2u + 6*uman*xd2u - 2*m1**2*xd33u + 2*m2**2*xd33u + 2*sqrtS**2*xd33u + m1**2*xd3u + 3*m2**2*xd3u + 
     &                      m3**2*xd3u + 3*m4**2*xd3u - 4*sqrtS**2*xd3u + 2*uman*xd3u)) + 
     &                g3Lu*(-2*m1**4*xc01u + 20*m1**2*m2**2*xc01u - 2*m2**4*xc01u + 4*m1**2*m3**2*xc01u + 4*m2**2*m3**2*xc01u + 4*m1**2*m4**2*xc01u + 4*m2**2*m4**2*xc01u + 
     &                   6*m1**2*sqrtS**2*xc01u + 6*m2**2*sqrtS**2*xc01u - 2*m3**2*sqrtS**2*xc01u - 2*m4**2*sqrtS**2*xc01u - 4*sqrtS**4*xc01u + 8*m1**2*m3**2*xc02u + 4*m2**2*m3**2*xc02u + 
     &                   4*m1**2*m4**2*xc02u - 3*m1**2*sqrtS**2*xc02u - m2**2*sqrtS**2*xc02u - 3*m3**2*sqrtS**2*xc02u - m4**2*sqrtS**2*xc02u + sqrtS**4*xc02u + 4*m2**2*m3**2*xc03u + 
     &                   4*m1**2*m4**2*xc03u + 8*m2**2*m4**2*xc03u - m1**2*sqrtS**2*xc03u - 3*m2**2*sqrtS**2*xc03u - m3**2*sqrtS**2*xc03u - 3*m4**2*sqrtS**2*xc03u + sqrtS**4*xc03u - 
     &                   8*m1**2*m4**2*xc04u + 8*m2**2*m4**2*xc04u + 2*m1**2*sqrtS**2*xc04u - 2*m2**2*sqrtS**2*xc04u - 4*m1**2*m3**2*xc111u + 4*m2**2*m3**2*xc111u - 12*m1**2*m4**2*xc111u + 
     &                   12*m2**2*m4**2*xc111u + 4*m1**2*sqrtS**2*xc111u - 4*m2**2*sqrtS**2*xc111u + 2*m1**4*xc11u - 2*m2**4*xc11u - 8*m1**2*m3**2*xc11u + 8*m2**2*m3**2*xc11u - 
     &                   12*m1**2*m4**2*xc11u + 12*m2**2*m4**2*xc11u + 4*m1**2*sqrtS**2*xc11u - 4*m2**2*sqrtS**2*xc11u - 32*m1**2*m4**2*xc121u + 32*m2**2*m4**2*xc121u + 
     &                   8*m1**2*sqrtS**2*xc121u - 8*m2**2*sqrtS**2*xc121u + 6*m1**2*m3**2*xc12u + 2*m2**2*m3**2*xc12u - 6*m1**2*m4**2*xc12u - 2*m2**2*m4**2*xc12u - 2*m3**2*sqrtS**2*xc12u + 
     &                   2*m4**2*sqrtS**2*xc12u + 2*m1**4*xc13u - 6*m1**2*m2**2*xc13u + 4*m2**4*xc13u - 4*m1**2*m3**2*xc13u - 4*m2**2*m3**2*xc13u + 6*m1**2*m4**2*xc13u + 
     &                   18*m2**2*m4**2*xc13u - 6*m1**2*sqrtS**2*xc13u - 10*m2**2*sqrtS**2*xc13u + 2*m3**2*sqrtS**2*xc13u - 6*m4**2*sqrtS**2*xc13u + 4*sqrtS**4*xc13u + 8*m1**2*m3**2*xc21u - 
     &                   8*m2**2*m3**2*xc21u - 24*m1**2*m4**2*xc21u + 24*m2**2*m4**2*xc21u + 4*m1**2*sqrtS**2*xc21u - 4*m2**2*sqrtS**2*xc21u + 16*m1**2*m3**2*xc221u - 16*m2**2*m3**2*xc221u - 
     &                   16*m1**2*m4**2*xc221u + 16*m2**2*m4**2*xc221u - 2*m1**4*xc22u - 14*m1**2*m2**2*xc22u + 12*m1**2*m3**2*xc22u + 4*m2**2*m3**2*xc22u - 10*m1**2*m4**2*xc22u - 
     &                   6*m2**2*m4**2*xc22u - 2*m1**2*sqrtS**2*xc22u - 2*m2**2*sqrtS**2*xc22u - 4*m3**2*sqrtS**2*xc22u + 4*m4**2*sqrtS**2*xc22u + 2*sqrtS**4*xc22u - 2*m1**2*m3**2*xc23u - 
     &                   6*m2**2*m3**2*xc23u + 2*m1**2*m4**2*xc23u + 6*m2**2*m4**2*xc23u + 2*m3**2*sqrtS**2*xc23u - 2*m4**2*sqrtS**2*xc23u - 2*M0Lpu**2*m1**4*xd0u + 
     &                   20*M0Lpu**2*m1**2*m2**2*xd0u - 2*M0Lpu**2*m2**4*xd0u - 2*m1**4*M2Lpu**2*xd0u + 20*m1**2*m2**2*M2Lpu**2*xd0u - 2*m2**4*M2Lpu**2*xd0u + 12*M0Lpu**2*m1**2*m3**2*xd0u + 
     &                   12*M0Lpu**2*m2**2*m3**2*xd0u - 12*m1**2*m2**2*m3**2*xd0u + 4*M1Lpu**2*m2**2*m3**2*xd0u - 4*m2**4*m3**2*xd0u + 4*m1**2*M2Lpu**2*m3**2*xd0u + 
     &                   4*m2**2*M2Lpu**2*m3**2*xd0u + 8*m1**2*m3**2*M3Lpu**2*xd0u + 4*m2**2*m3**2*M3Lpu**2*xd0u + 12*M0Lpu**2*m1**2*m4**2*xd0u - 4*m1**4*m4**2*xd0u + 
     &                   4*m1**2*M1Lpu**2*m4**2*xd0u + 12*M0Lpu**2*m2**2*m4**2*xd0u - 12*m1**2*m2**2*m4**2*xd0u + 8*M1Lpu**2*m2**2*m4**2*xd0u + 4*m1**2*M2Lpu**2*m4**2*xd0u + 
     &                   4*m2**2*M2Lpu**2*m4**2*xd0u + 4*m1**2*M3Lpu**2*m4**2*xd0u + 2*M0Lpu**2*m1**2*sqrtS**2*xd0u + m1**4*sqrtS**2*xd0u - m1**2*M1Lpu**2*sqrtS**2*xd0u + 
     &                   2*M0Lpu**2*m2**2*sqrtS**2*xd0u + 6*m1**2*m2**2*sqrtS**2*xd0u - 3*M1Lpu**2*m2**2*sqrtS**2*xd0u + m2**4*sqrtS**2*xd0u + 6*m1**2*M2Lpu**2*sqrtS**2*xd0u + 
     &                   6*m2**2*M2Lpu**2*sqrtS**2*xd0u - 6*M0Lpu**2*m3**2*sqrtS**2*xd0u + m1**2*m3**2*sqrtS**2*xd0u - M1Lpu**2*m3**2*sqrtS**2*xd0u + 3*m2**2*m3**2*sqrtS**2*xd0u - 
     &                   2*M2Lpu**2*m3**2*sqrtS**2*xd0u - 3*m1**2*M3Lpu**2*sqrtS**2*xd0u - m2**2*M3Lpu**2*sqrtS**2*xd0u - 3*m3**2*M3Lpu**2*sqrtS**2*xd0u - 6*M0Lpu**2*m4**2*sqrtS**2*xd0u + 
     &                   3*m1**2*m4**2*sqrtS**2*xd0u - 3*M1Lpu**2*m4**2*sqrtS**2*xd0u + m2**2*m4**2*sqrtS**2*xd0u - 2*M2Lpu**2*m4**2*sqrtS**2*xd0u - M3Lpu**2*m4**2*sqrtS**2*xd0u - 
     &                   2*M0Lpu**2*sqrtS**4*xd0u - m1**2*sqrtS**4*xd0u + M1Lpu**2*sqrtS**4*xd0u - m2**2*sqrtS**4*xd0u - 4*M2Lpu**2*sqrtS**4*xd0u + M3Lpu**2*sqrtS**4*xd0u + 
     &                   12*M0Lpu**2*m1**2*m3**2*xd11u - 6*m1**4*m3**2*xd11u + 4*M0Lpu**2*m2**2*m3**2*xd11u - 8*m1**2*m2**2*m3**2*xd11u - 2*m2**4*m3**2*xd11u - 
     &                   12*M0Lpu**2*m1**2*m4**2*xd11u + 6*m1**4*m4**2*xd11u - 4*M0Lpu**2*m2**2*m4**2*xd11u + 8*m1**2*m2**2*m4**2*xd11u + 2*m2**4*m4**2*xd11u - 
     &                   4*M0Lpu**2*m3**2*sqrtS**2*xd11u + 8*m1**2*m3**2*sqrtS**2*xd11u + 4*m2**2*m3**2*sqrtS**2*xd11u + 4*M0Lpu**2*m4**2*sqrtS**2*xd11u - 8*m1**2*m4**2*sqrtS**2*xd11u - 
     &                   4*m2**2*m4**2*sqrtS**2*xd11u - 2*m3**2*sqrtS**4*xd11u + 2*m4**2*sqrtS**4*xd11u + 32*M0Lpu**2*m1**2*m3**2*xd12u - 16*m1**4*m3**2*xd12u + 
     &                   16*M0Lpu**2*m2**2*m3**2*xd12u - 24*m1**2*m2**2*m3**2*xd12u - 8*m2**4*m3**2*xd12u + 16*M0Lpu**2*m1**2*m4**2*xd12u - 8*m1**4*m4**2*xd12u - 8*m1**2*m2**2*m4**2*xd12u - 
     &                   12*M0Lpu**2*m1**2*sqrtS**2*xd12u + 6*m1**4*sqrtS**2*xd12u - 4*M0Lpu**2*m2**2*sqrtS**2*xd12u + 8*m1**2*m2**2*sqrtS**2*xd12u + 2*m2**4*sqrtS**2*xd12u - 
     &                   12*M0Lpu**2*m3**2*sqrtS**2*xd12u + 22*m1**2*m3**2*sqrtS**2*xd12u + 14*m2**2*m3**2*sqrtS**2*xd12u - 4*M0Lpu**2*m4**2*sqrtS**2*xd12u + 10*m1**2*m4**2*sqrtS**2*xd12u + 
     &                   2*m2**2*m4**2*sqrtS**2*xd12u + 4*M0Lpu**2*sqrtS**4*xd12u - 8*m1**2*sqrtS**4*xd12u - 4*m2**2*sqrtS**4*xd12u - 6*m3**2*sqrtS**4*xd12u - 2*m4**2*sqrtS**4*xd12u + 
     &                   2*sqrtS**6*xd12u - 8*M0Lpu**2*m1**2*m3**2*xd13u + 4*m1**4*m3**2*xd13u + 8*M0Lpu**2*m2**2*m3**2*xd13u - 4*m2**4*m3**2*xd13u + 8*M0Lpu**2*m1**2*m4**2*xd13u - 
     &                   4*m1**4*m4**2*xd13u - 8*M0Lpu**2*m2**2*m4**2*xd13u + 4*m2**4*m4**2*xd13u - 4*m1**2*m3**2*sqrtS**2*xd13u + 4*m2**2*m3**2*sqrtS**2*xd13u + 
     &                   4*m1**2*m4**2*sqrtS**2*xd13u - 4*m2**2*m4**2*sqrtS**2*xd13u - 2*M0Lpu**2*m1**4*xd1u + 2*m1**4*M1Lpu**2*xd1u + 20*M0Lpu**2*m1**2*m2**2*xd1u - 
     &                   20*m1**2*M1Lpu**2*m2**2*xd1u - 2*M0Lpu**2*m2**4*xd1u + 2*M1Lpu**2*m2**4*xd1u + 24*M0Lpu**2*m1**2*m3**2*xd1u - 6*m1**4*m3**2*xd1u - 6*m1**2*M1Lpu**2*m3**2*xd1u + 
     &                   16*M0Lpu**2*m2**2*m3**2*xd1u - 20*m1**2*m2**2*m3**2*xd1u - 2*M1Lpu**2*m2**2*m3**2*xd1u - 6*m2**4*m3**2*xd1u + 6*m1**2*m3**2*M3Lpu**2*xd1u + 
     &                   2*m2**2*m3**2*M3Lpu**2*xd1u + 8*M0Lpu**2*m1**2*m4**2*xd1u - 8*m1**4*m4**2*xd1u - 2*m1**2*M1Lpu**2*m4**2*xd1u + 24*m1**2*m2**2*m4**2*xd1u - 
     &                   6*M1Lpu**2*m2**2*m4**2*xd1u + 4*m1**2*m3**2*m4**2*xd1u + 4*m2**2*m3**2*m4**2*xd1u - 6*m1**2*M3Lpu**2*m4**2*xd1u - 2*m2**2*M3Lpu**2*m4**2*xd1u + 4*m1**2*m4**4*xd1u + 
     &                   4*m2**2*m4**4*xd1u + 3*m1**4*sqrtS**2*xd1u - 6*m1**2*M1Lpu**2*sqrtS**2*xd1u + 4*M0Lpu**2*m2**2*sqrtS**2*xd1u + 4*m1**2*m2**2*sqrtS**2*xd1u - 
     &                   6*M1Lpu**2*m2**2*sqrtS**2*xd1u + m2**4*sqrtS**2*xd1u - 10*M0Lpu**2*m3**2*sqrtS**2*xd1u + 11*m1**2*m3**2*sqrtS**2*xd1u + 2*M1Lpu**2*m3**2*sqrtS**2*xd1u + 
     &                   9*m2**2*m3**2*sqrtS**2*xd1u - 2*m3**2*M3Lpu**2*sqrtS**2*xd1u - 2*M0Lpu**2*m4**2*sqrtS**2*xd1u + 11*m1**2*m4**2*sqrtS**2*xd1u + 2*M1Lpu**2*m4**2*sqrtS**2*xd1u + 
     &                   5*m2**2*m4**2*sqrtS**2*xd1u - 2*m3**2*m4**2*sqrtS**2*xd1u + 2*M3Lpu**2*m4**2*sqrtS**2*xd1u - 2*m4**4*sqrtS**2*xd1u - 2*M0Lpu**2*sqrtS**4*xd1u - 
     &                   4*m1**2*sqrtS**4*xd1u + 4*M1Lpu**2*sqrtS**4*xd1u - 2*m2**2*sqrtS**4*xd1u - 3*m3**2*sqrtS**4*xd1u - 5*m4**2*sqrtS**4*xd1u + sqrtS**6*xd1u + 
     &                   16*M0Lpu**2*m1**2*m3**2*xd22u - 8*m1**4*m3**2*xd22u + 16*M0Lpu**2*m2**2*m3**2*xd22u - 16*m1**2*m2**2*m3**2*xd22u - 8*m2**4*m3**2*xd22u + 
     &                   16*M0Lpu**2*m1**2*m4**2*xd22u - 8*m1**4*m4**2*xd22u + 16*M0Lpu**2*m2**2*m4**2*xd22u - 16*m1**2*m2**2*m4**2*xd22u - 8*m2**4*m4**2*xd22u - 
     &                   8*M0Lpu**2*m1**2*sqrtS**2*xd22u + 4*m1**4*sqrtS**2*xd22u - 8*M0Lpu**2*m2**2*sqrtS**2*xd22u + 8*m1**2*m2**2*sqrtS**2*xd22u + 4*m2**4*sqrtS**2*xd22u - 
     &                   8*M0Lpu**2*m3**2*sqrtS**2*xd22u + 12*m1**2*m3**2*sqrtS**2*xd22u + 12*m2**2*m3**2*sqrtS**2*xd22u - 8*M0Lpu**2*m4**2*sqrtS**2*xd22u + 12*m1**2*m4**2*sqrtS**2*xd22u + 
     &                   12*m2**2*m4**2*sqrtS**2*xd22u + 4*M0Lpu**2*sqrtS**4*xd22u - 6*m1**2*sqrtS**4*xd22u - 6*m2**2*sqrtS**4*xd22u - 4*m3**2*sqrtS**4*xd22u - 4*m4**2*sqrtS**4*xd22u + 
     &                   2*sqrtS**6*xd22u + 16*M0Lpu**2*m2**2*m3**2*xd23u - 8*m1**2*m2**2*m3**2*xd23u - 8*m2**4*m3**2*xd23u + 16*M0Lpu**2*m1**2*m4**2*xd23u - 8*m1**4*m4**2*xd23u + 
     &                   32*M0Lpu**2*m2**2*m4**2*xd23u - 24*m1**2*m2**2*m4**2*xd23u - 16*m2**4*m4**2*xd23u - 4*M0Lpu**2*m1**2*sqrtS**2*xd23u + 2*m1**4*sqrtS**2*xd23u - 
     &                   12*M0Lpu**2*m2**2*sqrtS**2*xd23u + 8*m1**2*m2**2*sqrtS**2*xd23u + 6*m2**4*sqrtS**2*xd23u - 4*M0Lpu**2*m3**2*sqrtS**2*xd23u + 2*m1**2*m3**2*sqrtS**2*xd23u + 
     &                   10*m2**2*m3**2*sqrtS**2*xd23u - 12*M0Lpu**2*m4**2*sqrtS**2*xd23u + 14*m1**2*m4**2*sqrtS**2*xd23u + 22*m2**2*m4**2*sqrtS**2*xd23u + 4*M0Lpu**2*sqrtS**4*xd23u - 
     &                   4*m1**2*sqrtS**4*xd23u - 8*m2**2*sqrtS**4*xd23u - 2*m3**2*sqrtS**4*xd23u - 6*m4**2*sqrtS**4*xd23u + 2*sqrtS**6*xd23u + 2*m1**4*M1Lpu**2*xd2u + 
     &                   20*M0Lpu**2*m1**2*m2**2*xd2u - 6*m1**2*M1Lpu**2*m2**2*xd2u - 4*M0Lpu**2*m2**4*xd2u + 4*M1Lpu**2*m2**4*xd2u + 24*M0Lpu**2*m1**2*m3**2*xd2u - 6*m1**4*m3**2*xd2u - 
     &                   4*m1**2*M1Lpu**2*m3**2*xd2u + 32*M0Lpu**2*m2**2*m3**2*xd2u - 14*m1**2*m2**2*m3**2*xd2u - 4*M1Lpu**2*m2**2*m3**2*xd2u - 12*m2**4*m3**2*xd2u - 4*m1**2*m3**4*xd2u - 
     &                   2*m1**4*M3Lpu**2*xd2u - 14*m1**2*m2**2*M3Lpu**2*xd2u + 12*m1**2*m3**2*M3Lpu**2*xd2u + 4*m2**2*m3**2*M3Lpu**2*xd2u + 36*M0Lpu**2*m1**2*m4**2*xd2u - 
     &                   14*m1**4*m4**2*xd2u + 6*m1**2*M1Lpu**2*m4**2*xd2u + 20*M0Lpu**2*m2**2*m4**2*xd2u - 22*m1**2*m2**2*m4**2*xd2u + 18*M1Lpu**2*m2**2*m4**2*xd2u - 12*m2**4*m4**2*xd2u + 
     &                   18*m1**2*m3**2*m4**2*xd2u + 14*m2**2*m3**2*m4**2*xd2u - 10*m1**2*M3Lpu**2*m4**2*xd2u - 6*m2**2*M3Lpu**2*m4**2*xd2u - 2*m1**2*m4**4*xd2u - 10*m2**2*m4**4*xd2u - 
     &                   8*M0Lpu**2*m1**2*sqrtS**2*xd2u + 5*m1**4*sqrtS**2*xd2u - 6*m1**2*M1Lpu**2*sqrtS**2*xd2u - 4*M0Lpu**2*m2**2*sqrtS**2*xd2u + 14*m1**2*m2**2*sqrtS**2*xd2u - 
     &                   10*M1Lpu**2*m2**2*sqrtS**2*xd2u + 5*m2**4*sqrtS**2*xd2u - 14*M0Lpu**2*m3**2*sqrtS**2*xd2u + 12*m1**2*m3**2*sqrtS**2*xd2u + 2*M1Lpu**2*m3**2*sqrtS**2*xd2u + 
     &                   16*m2**2*m3**2*sqrtS**2*xd2u + m3**4*sqrtS**2*xd2u - 2*m1**2*M3Lpu**2*sqrtS**2*xd2u - 2*m2**2*M3Lpu**2*sqrtS**2*xd2u - 4*m3**2*M3Lpu**2*sqrtS**2*xd2u - 
     &                   14*M0Lpu**2*m4**2*sqrtS**2*xd2u + 20*m1**2*m4**2*sqrtS**2*xd2u - 6*M1Lpu**2*m4**2*sqrtS**2*xd2u + 20*m2**2*m4**2*sqrtS**2*xd2u - 8*m3**2*m4**2*sqrtS**2*xd2u + 
     &                   4*M3Lpu**2*m4**2*sqrtS**2*xd2u + 3*m4**4*sqrtS**2*xd2u + 2*M0Lpu**2*sqrtS**4*xd2u - 7*m1**2*sqrtS**4*xd2u + 4*M1Lpu**2*sqrtS**4*xd2u - 7*m2**2*sqrtS**4*xd2u - 
     &                   5*m3**2*sqrtS**4*xd2u + 2*M3Lpu**2*sqrtS**4*xd2u - 7*m4**2*sqrtS**4*xd2u + 2*sqrtS**6*xd2u + 2*uman**3*xd2u - 4*M0Lpu**2*m1**2*m3**2*xd33u + 2*m1**4*m3**2*xd33u - 
     &                   12*M0Lpu**2*m2**2*m3**2*xd33u + 8*m1**2*m2**2*m3**2*xd33u + 6*m2**4*m3**2*xd33u + 4*M0Lpu**2*m1**2*m4**2*xd33u - 2*m1**4*m4**2*xd33u + 
     &                   12*M0Lpu**2*m2**2*m4**2*xd33u - 8*m1**2*m2**2*m4**2*xd33u - 6*m2**4*m4**2*xd33u + 4*M0Lpu**2*m3**2*sqrtS**2*xd33u - 4*m1**2*m3**2*sqrtS**2*xd33u - 
     &                   8*m2**2*m3**2*sqrtS**2*xd33u - 4*M0Lpu**2*m4**2*sqrtS**2*xd33u + 4*m1**2*m4**2*sqrtS**2*xd33u + 8*m2**2*m4**2*sqrtS**2*xd33u + 2*m3**2*sqrtS**4*xd33u - 
     &                   2*m4**2*sqrtS**4*xd33u - 2*M0Lpu**2*m1**4*xd3u + 20*M0Lpu**2*m1**2*m2**2*xd3u - 2*M0Lpu**2*m2**4*xd3u - 2*m1**2*M1Lpu**2*m3**2*xd3u + 8*M0Lpu**2*m2**2*m3**2*xd3u + 
     &                   24*m1**2*m2**2*m3**2*xd3u - 6*M1Lpu**2*m2**2*m3**2*xd3u - 8*m2**4*m3**2*xd3u + 4*m1**2*m3**4*xd3u + 4*m2**2*m3**4*xd3u + 2*m1**4*M3Lpu**2*xd3u - 
     &                   20*m1**2*m2**2*M3Lpu**2*xd3u + 2*m2**4*M3Lpu**2*xd3u - 6*m1**2*m3**2*M3Lpu**2*xd3u - 2*m2**2*m3**2*M3Lpu**2*xd3u + 16*M0Lpu**2*m1**2*m4**2*xd3u - 
     &                   6*m1**4*m4**2*xd3u + 2*m1**2*M1Lpu**2*m4**2*xd3u + 24*M0Lpu**2*m2**2*m4**2*xd3u - 20*m1**2*m2**2*m4**2*xd3u + 6*M1Lpu**2*m2**2*m4**2*xd3u - 6*m2**4*m4**2*xd3u + 
     &                   4*m1**2*m3**2*m4**2*xd3u + 4*m2**2*m3**2*m4**2*xd3u - 2*m1**2*M3Lpu**2*m4**2*xd3u - 6*m2**2*M3Lpu**2*m4**2*xd3u + 4*M0Lpu**2*m1**2*sqrtS**2*xd3u + 
     &                   m1**4*sqrtS**2*xd3u + 4*m1**2*m2**2*sqrtS**2*xd3u + 3*m2**4*sqrtS**2*xd3u - 2*M0Lpu**2*m3**2*sqrtS**2*xd3u + 5*m1**2*m3**2*sqrtS**2*xd3u + 
     &                   2*M1Lpu**2*m3**2*sqrtS**2*xd3u + 11*m2**2*m3**2*sqrtS**2*xd3u - 2*m3**4*sqrtS**2*xd3u - 6*m1**2*M3Lpu**2*sqrtS**2*xd3u - 6*m2**2*M3Lpu**2*sqrtS**2*xd3u + 
     &                   2*m3**2*M3Lpu**2*sqrtS**2*xd3u - 10*M0Lpu**2*m4**2*sqrtS**2*xd3u + 9*m1**2*m4**2*sqrtS**2*xd3u - 2*M1Lpu**2*m4**2*sqrtS**2*xd3u + 11*m2**2*m4**2*sqrtS**2*xd3u - 
     &                   2*m3**2*m4**2*sqrtS**2*xd3u + 2*M3Lpu**2*m4**2*sqrtS**2*xd3u - 2*M0Lpu**2*sqrtS**4*xd3u - 2*m1**2*sqrtS**4*xd3u - 4*m2**2*sqrtS**4*xd3u - 5*m3**2*sqrtS**4*xd3u + 
     &                   4*M3Lpu**2*sqrtS**4*xd3u - 3*m4**2*sqrtS**4*xd3u + sqrtS**6*xd3u + 
     &                   tman**2*(xc02u - 3*xc03u - 4*xc04u - 2*xc11u - 16*xc121u - 4*xc13u - 8*xc21u - 2*M0Lpu**2*xd0u + 3*m1**2*xd0u - 3*M1Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 
     &                      4*M0Lpu**2*xd12u - 2*m1**2*xd12u - 2*m2**2*xd12u + 2*sqrtS**2*xd12u + 2*M0Lpu**2*xd1u - m1**2*xd1u - m2**2*xd1u + sqrtS**2*xd1u + 4*M0Lpu**2*xd22u - 
     &                      2*m1**2*xd22u - 2*m2**2*xd22u + 2*sqrtS**2*xd22u - 12*M0Lpu**2*xd23u + 6*m1**2*xd23u + 6*m2**2*xd23u - 6*sqrtS**2*xd23u + 4*M0Lpu**2*xd2u + m1**2*xd2u - 
     &                      4*M1Lpu**2*xd2u - 3*m2**2*xd2u + m3**2*xd2u + m4**2*xd2u + 2*sqrtS**2*xd2u - 2*uman*xd2u - 6*M0Lpu**2*xd3u + 3*m1**2*xd3u + 3*m2**2*xd3u - 3*sqrtS**2*xd3u) + 
     &                   uman**2*(3*xc02u + 7*xc03u + 4*xc04u + 2*xc11u + 16*xc121u + 2*xc13u + 8*xc21u - 2*xc22u - 2*M0Lpu**2*xd0u + 3*m1**2*xd0u + 7*M1Lpu**2*xd0u + 7*m2**2*xd0u + 
     &                      3*M3Lpu**2*xd0u - 8*sqrtS**2*xd0u - 12*M0Lpu**2*xd12u + 6*m1**2*xd12u + 6*m2**2*xd12u - 6*sqrtS**2*xd12u - 6*M0Lpu**2*xd1u + 7*m1**2*xd1u + 7*m2**2*xd1u - 
     &                      5*sqrtS**2*xd1u - 12*M0Lpu**2*xd22u + 6*m1**2*xd22u + 6*m2**2*xd22u - 6*sqrtS**2*xd22u + 4*M0Lpu**2*xd23u - 2*m1**2*xd23u - 2*m2**2*xd23u + 2*sqrtS**2*xd23u - 
     &                      16*M0Lpu**2*xd2u + 9*m1**2*xd2u + 2*M1Lpu**2*xd2u + 9*m2**2*xd2u - 3*m3**2*xd2u - 2*M3Lpu**2*xd2u - 3*m4**2*xd2u - 12*sqrtS**2*xd2u + 2*M0Lpu**2*xd3u + 
     &                      3*m1**2*xd3u + 3*m2**2*xd3u - sqrtS**2*xd3u) + uman*
     &                    (12*xb01u + 12*xb02u - 48*xc001u + 12*M0Lpu**2*xc01u - 10*m1**2*xc01u - 10*m2**2*xc01u + 8*sqrtS**2*xc01u - 5*m1**2*xc02u + 12*M1Lpu**2*xc02u + 3*m2**2*xc02u + 
     &                      m3**2*xc02u - 3*m4**2*xc02u - 2*sqrtS**2*xc02u - m1**2*xc03u - 9*m2**2*xc03u - 7*m3**2*xc03u + 12*M3Lpu**2*xc03u - 3*m4**2*xc03u + 2*sqrtS**2*xc03u - 
     &                      4*m1**2*xc04u + 12*M1Lpu**2*xc04u - 4*m2**2*xc04u - 8*m3**2*xc04u + 6*sqrtS**2*xc04u + 4*m3**2*xc111u - 4*m4**2*xc111u - 4*sqrtS**2*xc111u - 4*m1**2*xc11u + 
     &                      4*m3**2*xc11u + 16*m4**2*xc11u - 16*m1**2*xc121u - 16*m2**2*xc121u - 8*m3**2*xc121u - 24*m4**2*xc121u - 2*m1**2*xc12u + 2*m2**2*xc12u - 2*sqrtS**2*xc12u + 
     &                      2*m1**2*xc13u - 8*m2**2*xc13u - 2*m3**2*xc13u - 4*m4**2*xc13u - 8*m1**2*xc21u - 8*m2**2*xc21u - 16*m3**2*xc21u + 4*sqrtS**2*xc21u - 16*sqrtS**2*xc221u + 
     &                      2*m1**2*xc22u + 8*m2**2*xc22u + 2*m4**2*xc22u - 6*sqrtS**2*xc22u + 2*m1**2*xc23u - 2*m2**2*xc23u - 2*sqrtS**2*xc23u - 48*M0Lpu**2*xd00u + 24*m1**2*xd00u + 
     &                      24*m2**2*xd00u - 24*sqrtS**2*xd00u - 4*M0Lpu**2*m1**2*xd0u + 3*m1**4*xd0u - m1**2*M1Lpu**2*xd0u - 4*M0Lpu**2*m2**2*xd0u - 18*m1**2*m2**2*xd0u - 
     &                      9*M1Lpu**2*m2**2*xd0u - m2**4*xd0u + 12*M0Lpu**2*M2Lpu**2*xd0u - 10*m1**2*M2Lpu**2*xd0u - 10*m2**2*M2Lpu**2*xd0u - 6*M0Lpu**2*m3**2*xd0u + 3*m1**2*m3**2*xd0u - 
     &                      7*M1Lpu**2*m3**2*xd0u - 5*m2**2*m3**2*xd0u - 5*m1**2*M3Lpu**2*xd0u + 12*M1Lpu**2*M3Lpu**2*xd0u + 3*m2**2*M3Lpu**2*xd0u + m3**2*M3Lpu**2*xd0u - 
     &                      6*M0Lpu**2*m4**2*xd0u - m1**2*m4**2*xd0u - 3*M1Lpu**2*m4**2*xd0u - m2**2*m4**2*xd0u - 3*M3Lpu**2*m4**2*xd0u + 8*M0Lpu**2*sqrtS**2*xd0u - 8*m1**2*sqrtS**2*xd0u + 
     &                      2*M1Lpu**2*sqrtS**2*xd0u - 4*m2**2*sqrtS**2*xd0u + 8*M2Lpu**2*sqrtS**2*xd0u + 2*m3**2*sqrtS**2*xd0u - 2*M3Lpu**2*sqrtS**2*xd0u + 2*m4**2*sqrtS**2*xd0u + 
     &                      4*sqrtS**4*xd0u - 4*M0Lpu**2*m1**2*xd11u + 2*m1**4*xd11u + 4*M0Lpu**2*m2**2*xd11u - 2*m2**4*xd11u - 4*M0Lpu**2*sqrtS**2*xd11u + 4*m2**2*sqrtS**2*xd11u - 
     &                      2*sqrtS**4*xd11u + 4*M0Lpu**2*m1**2*xd12u - 2*m1**4*xd12u + 12*M0Lpu**2*m2**2*xd12u - 8*m1**2*m2**2*xd12u - 6*m2**4*xd12u + 4*M0Lpu**2*m3**2*xd12u - 
     &                      2*m1**2*m3**2*xd12u - 2*m2**2*m3**2*xd12u + 12*M0Lpu**2*m4**2*xd12u - 6*m1**2*m4**2*xd12u - 6*m2**2*m4**2*xd12u - 8*M0Lpu**2*sqrtS**2*xd12u + 
     &                      6*m1**2*sqrtS**2*xd12u + 10*m2**2*sqrtS**2*xd12u + 2*m3**2*sqrtS**2*xd12u + 6*m4**2*sqrtS**2*xd12u - 4*sqrtS**4*xd12u + 8*M0Lpu**2*sqrtS**2*xd13u - 
     &                      4*m1**2*sqrtS**2*xd13u - 4*m2**2*sqrtS**2*xd13u + 4*sqrtS**4*xd13u - 4*M0Lpu**2*m1**2*xd1u + m1**4*xd1u + 4*m1**2*M1Lpu**2*xd1u + 4*M0Lpu**2*m2**2*xd1u - 
     &                      22*m1**2*m2**2*xd1u + 4*M1Lpu**2*m2**2*xd1u - 3*m2**4*xd1u + 2*M0Lpu**2*m3**2*xd1u - 5*m1**2*m3**2*xd1u - 5*m2**2*m3**2*xd1u - 2*m1**2*M3Lpu**2*xd1u + 
     &                      2*m2**2*M3Lpu**2*xd1u + 6*M0Lpu**2*m4**2*xd1u - 11*m1**2*m4**2*xd1u - 11*m2**2*m4**2*xd1u - 2*M0Lpu**2*sqrtS**2*xd1u - 5*m1**2*sqrtS**2*xd1u + 
     &                      m2**2*sqrtS**2*xd1u + 3*m3**2*sqrtS**2*xd1u - 2*M3Lpu**2*sqrtS**2*xd1u + 7*m4**2*sqrtS**2*xd1u + 2*sqrtS**4*xd1u + 8*M0Lpu**2*m1**2*xd22u - 4*m1**4*xd22u + 
     &                      8*M0Lpu**2*m2**2*xd22u - 8*m1**2*m2**2*xd22u - 4*m2**4*xd22u + 8*M0Lpu**2*m3**2*xd22u - 4*m1**2*m3**2*xd22u - 4*m2**2*m3**2*xd22u + 8*M0Lpu**2*m4**2*xd22u - 
     &                      4*m1**2*m4**2*xd22u - 4*m2**2*m4**2*xd22u - 8*M0Lpu**2*sqrtS**2*xd22u + 8*m1**2*sqrtS**2*xd22u + 8*m2**2*sqrtS**2*xd22u + 4*m3**2*sqrtS**2*xd22u + 
     &                      4*m4**2*sqrtS**2*xd22u - 4*sqrtS**4*xd22u - 4*M0Lpu**2*m1**2*xd23u + 2*m1**4*xd23u - 12*M0Lpu**2*m2**2*xd23u + 8*m1**2*m2**2*xd23u + 6*m2**4*xd23u - 
     &                      4*M0Lpu**2*m3**2*xd23u + 2*m1**2*m3**2*xd23u + 2*m2**2*m3**2*xd23u - 12*M0Lpu**2*m4**2*xd23u + 6*m1**2*m4**2*xd23u + 6*m2**2*m4**2*xd23u + 
     &                      8*M0Lpu**2*sqrtS**2*xd23u - 6*m1**2*sqrtS**2*xd23u - 10*m2**2*sqrtS**2*xd23u - 2*m3**2*sqrtS**2*xd23u - 6*m4**2*sqrtS**2*xd23u + 4*sqrtS**4*xd23u + 
     &                      4*M0Lpu**2*m1**2*xd2u - 3*m1**4*xd2u + 2*m1**2*M1Lpu**2*xd2u + 8*M0Lpu**2*m2**2*xd2u - 26*m1**2*m2**2*xd2u - 8*M1Lpu**2*m2**2*xd2u - 3*m2**4*xd2u + 
     &                      10*M0Lpu**2*m3**2*xd2u - 4*m1**2*m3**2*xd2u - 2*M1Lpu**2*m3**2*xd2u - 18*m2**2*m3**2*xd2u + m3**4*xd2u + 2*m1**2*M3Lpu**2*xd2u + 8*m2**2*M3Lpu**2*xd2u + 
     &                      10*M0Lpu**2*m4**2*xd2u - 16*m1**2*m4**2*xd2u - 4*M1Lpu**2*m4**2*xd2u + 2*m2**2*m4**2*xd2u + 2*m3**2*m4**2*xd2u + 2*M3Lpu**2*m4**2*xd2u + m4**4*xd2u - 
     &                      2*M0Lpu**2*sqrtS**2*xd2u + 2*m1**2*sqrtS**2*xd2u + 2*m2**2*sqrtS**2*xd2u + 10*m3**2*sqrtS**2*xd2u - 6*M3Lpu**2*sqrtS**2*xd2u + 8*m4**2*sqrtS**2*xd2u + 
     &                      4*M0Lpu**2*m1**2*xd33u - 2*m1**4*xd33u - 4*M0Lpu**2*m2**2*xd33u + 2*m2**4*xd33u - 4*M0Lpu**2*sqrtS**2*xd33u + 4*m1**2*sqrtS**2*xd33u - 2*sqrtS**4*xd33u - 
     &                      4*M0Lpu**2*m1**2*xd3u + m1**4*xd3u + 2*m1**2*M1Lpu**2*xd3u - 12*M0Lpu**2*m2**2*xd3u - 14*m1**2*m2**2*xd3u - 2*M1Lpu**2*m2**2*xd3u + 5*m2**4*xd3u - 
     &                      2*M0Lpu**2*m3**2*xd3u - 7*m1**2*m3**2*xd3u - 7*m2**2*m3**2*xd3u + 4*m1**2*M3Lpu**2*xd3u + 4*m2**2*M3Lpu**2*xd3u - 6*M0Lpu**2*m4**2*xd3u - m1**2*m4**2*xd3u - 
     &                      m2**2*m4**2*xd3u + 6*M0Lpu**2*sqrtS**2*xd3u - 7*m1**2*sqrtS**2*xd3u - 2*M1Lpu**2*sqrtS**2*xd3u - 13*m2**2*sqrtS**2*xd3u + 3*m3**2*sqrtS**2*xd3u - 
     &                      m4**2*sqrtS**2*xd3u + 6*sqrtS**4*xd3u) + tman*(-12*xb01u - 12*xb02u + 48*xc001u - 12*M0Lpu**2*xc01u + 2*m1**2*xc01u + 2*m2**2*xc01u - 4*sqrtS**2*xc01u + 
     &                      5*m1**2*xc02u - 12*M1Lpu**2*xc02u - 3*m2**2*xc02u - m3**2*xc02u + 3*m4**2*xc02u + 2*sqrtS**2*xc02u - 4*uman*xc02u + m1**2*xc03u + 9*m2**2*xc03u + 7*m3**2*xc03u - 
     &                      12*M3Lpu**2*xc03u + 3*m4**2*xc03u - 2*sqrtS**2*xc03u - 4*uman*xc03u + 4*m1**2*xc04u - 12*M1Lpu**2*xc04u + 4*m2**2*xc04u + 8*m3**2*xc04u - 6*sqrtS**2*xc04u - 
     &                      4*m3**2*xc111u + 4*m4**2*xc111u + 4*sqrtS**2*xc111u + 4*m2**2*xc11u - 4*m3**2*xc11u - 16*m4**2*xc11u + 16*m1**2*xc121u + 16*m2**2*xc121u + 8*m3**2*xc121u + 
     &                      24*m4**2*xc121u + 2*m1**2*xc12u - 2*m2**2*xc12u + 2*sqrtS**2*xc12u + 6*m2**2*xc13u + 2*m3**2*xc13u + 4*m4**2*xc13u + 2*uman*xc13u + 8*m1**2*xc21u + 
     &                      8*m2**2*xc21u + 16*m3**2*xc21u - 4*sqrtS**2*xc21u + 16*sqrtS**2*xc221u + 8*m1**2*xc22u - 2*m2**2*xc22u - 2*m4**2*xc22u + 2*sqrtS**2*xc22u + 2*uman*xc22u - 
     &                      2*m1**2*xc23u + 2*m2**2*xc23u + 2*sqrtS**2*xc23u + 48*M0Lpu**2*xd00u - 24*m1**2*xd00u - 24*m2**2*xd00u + 24*sqrtS**2*xd00u - 4*M0Lpu**2*m1**2*xd0u - m1**4*xd0u + 
     &                      m1**2*M1Lpu**2*xd0u - 4*M0Lpu**2*m2**2*xd0u - 2*m1**2*m2**2*xd0u + 9*M1Lpu**2*m2**2*xd0u + 3*m2**4*xd0u - 12*M0Lpu**2*M2Lpu**2*xd0u + 2*m1**2*M2Lpu**2*xd0u + 
     &                      2*m2**2*M2Lpu**2*xd0u + 6*M0Lpu**2*m3**2*xd0u - 7*m1**2*m3**2*xd0u + 7*M1Lpu**2*m3**2*xd0u + m2**2*m3**2*xd0u + 5*m1**2*M3Lpu**2*xd0u - 
     &                      12*M1Lpu**2*M3Lpu**2*xd0u - 3*m2**2*M3Lpu**2*xd0u - m3**2*M3Lpu**2*xd0u + 6*M0Lpu**2*m4**2*xd0u - 3*m1**2*m4**2*xd0u + 3*M1Lpu**2*m4**2*xd0u - 
     &                      3*m2**2*m4**2*xd0u + 3*M3Lpu**2*m4**2*xd0u - 4*M0Lpu**2*sqrtS**2*xd0u + 2*m1**2*sqrtS**2*xd0u - 2*M1Lpu**2*sqrtS**2*xd0u - 2*m2**2*sqrtS**2*xd0u - 
     &                      4*M2Lpu**2*sqrtS**2*xd0u + 2*M3Lpu**2*sqrtS**2*xd0u + 4*M0Lpu**2*uman*xd0u + 2*m1**2*uman*xd0u - 4*M1Lpu**2*uman*xd0u + 2*m2**2*uman*xd0u - 4*M3Lpu**2*uman*xd0u + 
     &                      4*sqrtS**2*uman*xd0u + 4*M0Lpu**2*m1**2*xd11u - 2*m1**4*xd11u - 4*M0Lpu**2*m2**2*xd11u + 2*m2**4*xd11u + 4*M0Lpu**2*sqrtS**2*xd11u - 4*m2**2*sqrtS**2*xd11u + 
     &                      2*sqrtS**4*xd11u - 4*M0Lpu**2*m1**2*xd12u + 2*m1**4*xd12u - 12*M0Lpu**2*m2**2*xd12u + 8*m1**2*m2**2*xd12u + 6*m2**4*xd12u - 4*M0Lpu**2*m3**2*xd12u + 
     &                      2*m1**2*m3**2*xd12u + 2*m2**2*m3**2*xd12u - 12*M0Lpu**2*m4**2*xd12u + 6*m1**2*m4**2*xd12u + 6*m2**2*m4**2*xd12u + 8*M0Lpu**2*sqrtS**2*xd12u - 
     &                      6*m1**2*sqrtS**2*xd12u - 10*m2**2*sqrtS**2*xd12u - 2*m3**2*sqrtS**2*xd12u - 6*m4**2*sqrtS**2*xd12u + 4*sqrtS**4*xd12u + 8*M0Lpu**2*uman*xd12u - 
     &                      4*m1**2*uman*xd12u - 4*m2**2*uman*xd12u + 4*sqrtS**2*uman*xd12u - 8*M0Lpu**2*sqrtS**2*xd13u + 4*m1**2*sqrtS**2*xd13u + 4*m2**2*sqrtS**2*xd13u - 4*sqrtS**4*xd13u - 
     &                      4*M0Lpu**2*m1**2*xd1u + m1**4*xd1u + 4*m1**2*M1Lpu**2*xd1u - 12*M0Lpu**2*m2**2*xd1u + 2*m1**2*m2**2*xd1u + 4*M1Lpu**2*m2**2*xd1u + 5*m2**4*xd1u - 
     &                      2*M0Lpu**2*m3**2*xd1u + m1**2*m3**2*xd1u + m2**2*m3**2*xd1u + 2*m1**2*M3Lpu**2*xd1u - 2*m2**2*M3Lpu**2*xd1u - 6*M0Lpu**2*m4**2*xd1u - m1**2*m4**2*xd1u - 
     &                      m2**2*m4**2*xd1u + 6*M0Lpu**2*sqrtS**2*xd1u - m1**2*sqrtS**2*xd1u - 4*M1Lpu**2*sqrtS**2*xd1u - 7*m2**2*sqrtS**2*xd1u - m3**2*sqrtS**2*xd1u + 
     &                      2*M3Lpu**2*sqrtS**2*xd1u - m4**2*sqrtS**2*xd1u + 2*sqrtS**4*xd1u + 4*M0Lpu**2*uman*xd1u + 2*m1**2*uman*xd1u + 2*m2**2*uman*xd1u - 8*M0Lpu**2*m1**2*xd22u + 
     &                      4*m1**4*xd22u - 8*M0Lpu**2*m2**2*xd22u + 8*m1**2*m2**2*xd22u + 4*m2**4*xd22u - 8*M0Lpu**2*m3**2*xd22u + 4*m1**2*m3**2*xd22u + 4*m2**2*m3**2*xd22u - 
     &                      8*M0Lpu**2*m4**2*xd22u + 4*m1**2*m4**2*xd22u + 4*m2**2*m4**2*xd22u + 8*M0Lpu**2*sqrtS**2*xd22u - 8*m1**2*sqrtS**2*xd22u - 8*m2**2*sqrtS**2*xd22u - 
     &                      4*m3**2*sqrtS**2*xd22u - 4*m4**2*sqrtS**2*xd22u + 4*sqrtS**4*xd22u + 8*M0Lpu**2*uman*xd22u - 4*m1**2*uman*xd22u - 4*m2**2*uman*xd22u + 4*sqrtS**2*uman*xd22u + 
     &                      4*M0Lpu**2*m1**2*xd23u - 2*m1**4*xd23u + 12*M0Lpu**2*m2**2*xd23u - 8*m1**2*m2**2*xd23u - 6*m2**4*xd23u + 4*M0Lpu**2*m3**2*xd23u - 2*m1**2*m3**2*xd23u - 
     &                      2*m2**2*m3**2*xd23u + 12*M0Lpu**2*m4**2*xd23u - 6*m1**2*m4**2*xd23u - 6*m2**2*m4**2*xd23u - 8*M0Lpu**2*sqrtS**2*xd23u + 6*m1**2*sqrtS**2*xd23u + 
     &                      10*m2**2*sqrtS**2*xd23u + 2*m3**2*sqrtS**2*xd23u + 6*m4**2*sqrtS**2*xd23u - 4*sqrtS**4*xd23u + 8*M0Lpu**2*uman*xd23u - 4*m1**2*uman*xd23u - 4*m2**2*uman*xd23u + 
     &                      4*sqrtS**2*uman*xd23u - 16*M0Lpu**2*m1**2*xd2u + 3*m1**4*xd2u - 12*M0Lpu**2*m2**2*xd2u + 6*m1**2*m2**2*xd2u + 6*M1Lpu**2*m2**2*xd2u + 7*m2**4*xd2u - 
     &                      10*M0Lpu**2*m3**2*xd2u - 6*m1**2*m3**2*xd2u + 2*M1Lpu**2*m3**2*xd2u + 4*m2**2*m3**2*xd2u - m3**4*xd2u + 8*m1**2*M3Lpu**2*xd2u - 2*m2**2*M3Lpu**2*xd2u - 
     &                      10*M0Lpu**2*m4**2*xd2u + 2*m1**2*m4**2*xd2u + 4*M1Lpu**2*m4**2*xd2u + 4*m2**2*m4**2*xd2u - 2*m3**2*m4**2*xd2u - 2*M3Lpu**2*m4**2*xd2u - m4**4*xd2u + 
     &                      6*M0Lpu**2*sqrtS**2*xd2u - 6*m1**2*sqrtS**2*xd2u - 10*m2**2*sqrtS**2*xd2u - 4*m3**2*sqrtS**2*xd2u + 2*M3Lpu**2*sqrtS**2*xd2u - 6*m4**2*sqrtS**2*xd2u + 
     &                      4*sqrtS**4*xd2u + 12*M0Lpu**2*uman*xd2u + 2*m1**2*uman*xd2u + 2*M1Lpu**2*uman*xd2u - 2*m2**2*uman*xd2u + 2*m3**2*uman*xd2u + 2*M3Lpu**2*uman*xd2u + 
     &                      2*m4**2*uman*xd2u + 6*sqrtS**2*uman*xd2u - 4*M0Lpu**2*m1**2*xd33u + 2*m1**4*xd33u + 4*M0Lpu**2*m2**2*xd33u - 2*m2**4*xd33u + 4*M0Lpu**2*sqrtS**2*xd33u - 
     &                      4*m1**2*sqrtS**2*xd33u + 2*sqrtS**4*xd33u - 4*M0Lpu**2*m1**2*xd3u + m1**4*xd3u - 2*m1**2*M1Lpu**2*xd3u + 4*M0Lpu**2*m2**2*xd3u - 6*m1**2*m2**2*xd3u + 
     &                      2*M1Lpu**2*m2**2*xd3u - 3*m2**4*xd3u + 2*M0Lpu**2*m3**2*xd3u - 5*m1**2*m3**2*xd3u - 5*m2**2*m3**2*xd3u + 4*m1**2*M3Lpu**2*xd3u + 4*m2**2*M3Lpu**2*xd3u + 
     &                      6*M0Lpu**2*m4**2*xd3u - 3*m1**2*m4**2*xd3u - 3*m2**2*m4**2*xd3u - 2*M0Lpu**2*sqrtS**2*xd3u + m1**2*sqrtS**2*xd3u + 2*M1Lpu**2*sqrtS**2*xd3u + 
     &                      7*m2**2*sqrtS**2*xd3u + 3*m3**2*sqrtS**2*xd3u - 4*M3Lpu**2*sqrtS**2*xd3u + 3*m4**2*sqrtS**2*xd3u - 2*sqrtS**4*xd3u + 4*M0Lpu**2*uman*xd3u + 2*m1**2*uman*xd3u + 
     &                      2*m2**2*uman*xd3u))))) + g0Lu*(2*g1Lu*M1Lpu*(g2Ru*M2Lpu*
     &              (12*g3Lu*M0Lpu*M3Lpu*(tman - uman)*xd0u + g3Ru*(-2*m1**4*xd1u + 20*m1**2*m2**2*xd1u - 2*m2**4*xd1u + 4*m1**2*m3**2*xd1u + 4*m2**2*m3**2*xd1u + 4*m1**2*m4**2*xd1u + 
     &                   4*m2**2*m4**2*xd1u + 6*m1**2*sqrtS**2*xd1u + 6*m2**2*sqrtS**2*xd1u - 2*m3**2*sqrtS**2*xd1u - 2*m4**2*sqrtS**2*xd1u - 4*sqrtS**4*xd1u - 2*m1**4*xd2u + 
     &                   6*m1**2*m2**2*xd2u - 4*m2**4*xd2u + 4*m1**2*m3**2*xd2u + 8*m2**2*m3**2*xd2u - 2*m1**2*m4**2*xd2u - 10*m2**2*m4**2*xd2u + 5*m1**2*sqrtS**2*xd2u + 
     &                   7*m2**2*sqrtS**2*xd2u - 3*m3**2*sqrtS**2*xd2u + 3*m4**2*sqrtS**2*xd2u - 3*sqrtS**4*xd2u + tman**2*xd2u - uman**2*xd2u + 
     &                   tman*(6*xc02u + 6*xc04u + 6*M0Lpu**2*xd0u - 6*m2**2*xd0u + 6*M3Lpu**2*xd0u - 4*m1**2*xd1u - 4*m2**2*xd1u + 2*sqrtS**2*xd1u + m1**2*xd2u - 3*m2**2*xd2u - m3**2*xd2u - 
     &                      m4**2*xd2u - 2*sqrtS**2*xd2u) + uman*(-6*xc02u - 6*xc04u - 6*M0Lpu**2*xd0u + 6*m2**2*xd0u - 6*M3Lpu**2*xd0u - 4*m1**2*xd1u - 4*m2**2*xd1u + 2*sqrtS**2*xd1u - 
     &                      3*m1**2*xd2u + 5*m2**2*xd2u + m3**2*xd2u + m4**2*xd2u + 2*sqrtS**2*xd2u))) + 
     &             g2Lu*(g3Lu*M0Lpu*(4*m2**2*m3**2*xd0u + 4*m1**2*m4**2*xd0u + 8*m2**2*m4**2*xd0u - m1**2*sqrtS**2*xd0u - 3*m2**2*sqrtS**2*xd0u - m3**2*sqrtS**2*xd0u - 
     &                   3*m4**2*sqrtS**2*xd0u + sqrtS**4*xd0u + 6*m1**2*m3**2*xd11u + 2*m2**2*m3**2*xd11u - 6*m1**2*m4**2*xd11u - 2*m2**2*m4**2*xd11u - 2*m3**2*sqrtS**2*xd11u + 
     &                   2*m4**2*sqrtS**2*xd11u + 16*m1**2*m3**2*xd12u + 8*m2**2*m3**2*xd12u + 8*m1**2*m4**2*xd12u - 6*m1**2*sqrtS**2*xd12u - 2*m2**2*sqrtS**2*xd12u - 
     &                   6*m3**2*sqrtS**2*xd12u - 2*m4**2*sqrtS**2*xd12u + 2*sqrtS**4*xd12u - 4*m1**2*m3**2*xd13u + 4*m2**2*m3**2*xd13u + 4*m1**2*m4**2*xd13u - 4*m2**2*m4**2*xd13u + 
     &                   6*m1**2*m3**2*xd1u + 6*m2**2*m3**2*xd1u + 6*m1**2*m4**2*xd1u - 2*m2**2*m4**2*xd1u - 3*m1**2*sqrtS**2*xd1u - m2**2*sqrtS**2*xd1u - 3*m3**2*sqrtS**2*xd1u - 
     &                   m4**2*sqrtS**2*xd1u + sqrtS**4*xd1u + 8*m1**2*m3**2*xd22u + 8*m2**2*m3**2*xd22u + 8*m1**2*m4**2*xd22u + 8*m2**2*m4**2*xd22u - 4*m1**2*sqrtS**2*xd22u - 
     &                   4*m2**2*sqrtS**2*xd22u - 4*m3**2*sqrtS**2*xd22u - 4*m4**2*sqrtS**2*xd22u + 2*sqrtS**4*xd22u + 8*m2**2*m3**2*xd23u + 8*m1**2*m4**2*xd23u + 16*m2**2*m4**2*xd23u - 
     &                   2*m1**2*sqrtS**2*xd23u - 6*m2**2*sqrtS**2*xd23u - 2*m3**2*sqrtS**2*xd23u - 6*m4**2*sqrtS**2*xd23u + 2*sqrtS**4*xd23u + 8*m1**2*m3**2*xd2u + 12*m2**2*m3**2*xd2u + 
     &                   12*m1**2*m4**2*xd2u + 16*m2**2*m4**2*xd2u - 5*m1**2*sqrtS**2*xd2u - 7*m2**2*sqrtS**2*xd2u - 5*m3**2*sqrtS**2*xd2u - 7*m4**2*sqrtS**2*xd2u + 3*sqrtS**4*xd2u - 
     &                   2*m1**2*m3**2*xd33u - 6*m2**2*m3**2*xd33u + 2*m1**2*m4**2*xd33u + 6*m2**2*m4**2*xd33u + 2*m3**2*sqrtS**2*xd33u - 2*m4**2*sqrtS**2*xd33u + 
     &                   tman**2*(-3*xd0u + 2*xd12u + xd1u + 2*xd22u - 6*xd23u - xd2u - 3*xd3u) - 2*m1**2*m3**2*xd3u - 2*m2**2*m3**2*xd3u + 6*m1**2*m4**2*xd3u + 14*m2**2*m4**2*xd3u - 
     &                   m1**2*sqrtS**2*xd3u - 3*m2**2*sqrtS**2*xd3u + m3**2*sqrtS**2*xd3u - 5*m4**2*sqrtS**2*xd3u + sqrtS**4*xd3u + 
     &                   uman**2*(xd0u - 6*xd12u - 3*xd1u - 6*xd22u + 2*xd23u - 5*xd2u + xd3u) + 
     &                   uman*(6*xc01u + 6*xc02u - 24*xd00u - m1**2*xd0u - 3*m2**2*xd0u + 6*M2Lpu**2*xd0u - 7*m3**2*xd0u + 6*M3Lpu**2*xd0u - 3*m4**2*xd0u + 2*sqrtS**2*xd0u - 2*m1**2*xd11u + 
     &                      2*m2**2*xd11u - 2*sqrtS**2*xd11u + 2*m1**2*xd12u + 6*m2**2*xd12u + 2*m3**2*xd12u + 6*m4**2*xd12u - 4*sqrtS**2*xd12u + 4*sqrtS**2*xd13u + m1**2*xd1u + 
     &                      3*m2**2*xd1u + m3**2*xd1u + 3*m4**2*xd1u + 4*m1**2*xd22u + 4*m2**2*xd22u + 4*m3**2*xd22u + 4*m4**2*xd22u - 4*sqrtS**2*xd22u - 2*m1**2*xd23u - 6*m2**2*xd23u - 
     &                      2*m3**2*xd23u - 6*m4**2*xd23u + 4*sqrtS**2*xd23u + 3*m1**2*xd2u + m2**2*xd2u + 3*m3**2*xd2u + m4**2*xd2u - 2*sqrtS**2*xd2u + 2*m1**2*xd33u - 2*m2**2*xd33u - 
     &                      2*sqrtS**2*xd33u + m1**2*xd3u - 5*m2**2*xd3u - m3**2*xd3u - 3*m4**2*xd3u) - 
     &                   tman*(6*xc01u + 6*xc02u - 24*xd00u - m1**2*xd0u - 3*m2**2*xd0u + 6*M2Lpu**2*xd0u - 7*m3**2*xd0u + 6*M3Lpu**2*xd0u - 3*m4**2*xd0u + 2*sqrtS**2*xd0u - 2*uman*xd0u - 
     &                      2*m1**2*xd11u + 2*m2**2*xd11u - 2*sqrtS**2*xd11u + 2*m1**2*xd12u + 6*m2**2*xd12u + 2*m3**2*xd12u + 6*m4**2*xd12u - 4*sqrtS**2*xd12u - 4*uman*xd12u + 
     &                      4*sqrtS**2*xd13u + m1**2*xd1u + 3*m2**2*xd1u + m3**2*xd1u + 3*m4**2*xd1u - 2*uman*xd1u + 4*m1**2*xd22u + 4*m2**2*xd22u + 4*m3**2*xd22u + 4*m4**2*xd22u - 
     &                      4*sqrtS**2*xd22u - 4*uman*xd22u - 2*m1**2*xd23u - 6*m2**2*xd23u - 2*m3**2*xd23u - 6*m4**2*xd23u + 4*sqrtS**2*xd23u - 4*uman*xd23u + 3*m1**2*xd2u + m2**2*xd2u + 
     &                      3*m3**2*xd2u + m4**2*xd2u - 2*sqrtS**2*xd2u - 6*uman*xd2u + 2*m1**2*xd33u - 2*m2**2*xd33u - 2*sqrtS**2*xd33u + m1**2*xd3u - 5*m2**2*xd3u - m3**2*xd3u - 
     &                      3*m4**2*xd3u - 2*uman*xd3u)) + g3Ru*M3Lpu*(6*m1**2*m3**2*xd11u + 2*m2**2*m3**2*xd11u - 6*m1**2*m4**2*xd11u - 2*m2**2*m4**2*xd11u - 2*m3**2*sqrtS**2*xd11u + 
     &                   2*m4**2*sqrtS**2*xd11u + 16*m1**2*m3**2*xd12u + 8*m2**2*m3**2*xd12u + 8*m1**2*m4**2*xd12u - 6*m1**2*sqrtS**2*xd12u - 2*m2**2*sqrtS**2*xd12u - 
     &                   6*m3**2*sqrtS**2*xd12u - 2*m4**2*sqrtS**2*xd12u + 2*sqrtS**4*xd12u - 4*m1**2*m3**2*xd13u + 4*m2**2*m3**2*xd13u + 4*m1**2*m4**2*xd13u - 4*m2**2*m4**2*xd13u + 
     &                   8*m1**2*m3**2*xd1u + 4*m2**2*m3**2*xd1u + 4*m1**2*m4**2*xd1u - 3*m1**2*sqrtS**2*xd1u - m2**2*sqrtS**2*xd1u - 3*m3**2*sqrtS**2*xd1u - m4**2*sqrtS**2*xd1u + 
     &                   sqrtS**4*xd1u + 8*m1**2*m3**2*xd22u + 8*m2**2*m3**2*xd22u + 8*m1**2*m4**2*xd22u + 8*m2**2*m4**2*xd22u - 4*m1**2*sqrtS**2*xd22u - 4*m2**2*sqrtS**2*xd22u - 
     &                   4*m3**2*sqrtS**2*xd22u - 4*m4**2*sqrtS**2*xd22u + 2*sqrtS**4*xd22u + 8*m2**2*m3**2*xd23u + 8*m1**2*m4**2*xd23u + 16*m2**2*m4**2*xd23u - 2*m1**2*sqrtS**2*xd23u - 
     &                   6*m2**2*sqrtS**2*xd23u - 2*m3**2*sqrtS**2*xd23u - 6*m4**2*sqrtS**2*xd23u + 2*sqrtS**4*xd23u + 8*m1**2*m3**2*xd2u + 8*m2**2*m3**2*xd2u + 8*m1**2*m4**2*xd2u + 
     &                   8*m2**2*m4**2*xd2u - 4*m1**2*sqrtS**2*xd2u - 4*m2**2*sqrtS**2*xd2u - 4*m3**2*sqrtS**2*xd2u - 4*m4**2*sqrtS**2*xd2u + 2*sqrtS**4*xd2u - 2*m1**2*m3**2*xd33u - 
     &                   6*m2**2*m3**2*xd33u + 2*m1**2*m4**2*xd33u + 6*m2**2*m4**2*xd33u + 2*m3**2*sqrtS**2*xd33u - 2*m4**2*sqrtS**2*xd33u + 
     &                   tman**2*(2*xd12u + xd1u + 2*xd22u - 6*xd23u + 2*xd2u - 3*xd3u) + 4*m2**2*m3**2*xd3u + 4*m1**2*m4**2*xd3u + 8*m2**2*m4**2*xd3u - m1**2*sqrtS**2*xd3u - 
     &                   3*m2**2*sqrtS**2*xd3u - m3**2*sqrtS**2*xd3u - 3*m4**2*sqrtS**2*xd3u + sqrtS**4*xd3u + uman**2*(-6*xd0u - 6*xd12u - 3*xd1u - 6*xd22u + 2*xd23u - 6*xd2u + xd3u) + 
     &                   uman*(6*xc01u + 6*xc04u - 24*xd00u + 6*M0Lpu**2*xd0u + 6*M2Lpu**2*xd0u - 2*m1**2*xd11u + 2*m2**2*xd11u - 2*sqrtS**2*xd11u + 2*m1**2*xd12u + 6*m2**2*xd12u + 
     &                      2*m3**2*xd12u + 6*m4**2*xd12u - 4*sqrtS**2*xd12u + 4*sqrtS**2*xd13u + m1**2*xd1u + 3*m2**2*xd1u + m3**2*xd1u + 3*m4**2*xd1u - 2*sqrtS**2*xd1u + 4*m1**2*xd22u + 
     &                      4*m2**2*xd22u + 4*m3**2*xd22u + 4*m4**2*xd22u - 4*sqrtS**2*xd22u - 2*m1**2*xd23u - 6*m2**2*xd23u - 2*m3**2*xd23u - 6*m4**2*xd23u + 4*sqrtS**2*xd23u + 
     &                      4*m1**2*xd2u + 4*m2**2*xd2u + 4*m3**2*xd2u + 4*m4**2*xd2u - 4*sqrtS**2*xd2u + 2*m1**2*xd33u - 2*m2**2*xd33u - 2*sqrtS**2*xd33u - m1**2*xd3u - 3*m2**2*xd3u - 
     &                      m3**2*xd3u - 3*m4**2*xd3u + 2*sqrtS**2*xd3u) + tman*
     &                    (-6*xc01u - 6*xc04u + 24*xd00u - 6*M0Lpu**2*xd0u - 6*M2Lpu**2*xd0u + 6*uman*xd0u + 2*m1**2*xd11u - 2*m2**2*xd11u + 2*sqrtS**2*xd11u - 2*m1**2*xd12u - 
     &                      6*m2**2*xd12u - 2*m3**2*xd12u - 6*m4**2*xd12u + 4*sqrtS**2*xd12u + 4*uman*xd12u - 4*sqrtS**2*xd13u - m1**2*xd1u - 3*m2**2*xd1u - m3**2*xd1u - 3*m4**2*xd1u + 
     &                      2*sqrtS**2*xd1u + 2*uman*xd1u - 4*m1**2*xd22u - 4*m2**2*xd22u - 4*m3**2*xd22u - 4*m4**2*xd22u + 4*sqrtS**2*xd22u + 4*uman*xd22u + 2*m1**2*xd23u + 6*m2**2*xd23u + 
     &                      2*m3**2*xd23u + 6*m4**2*xd23u - 4*sqrtS**2*xd23u + 4*uman*xd23u - 4*m1**2*xd2u - 4*m2**2*xd2u - 4*m3**2*xd2u - 4*m4**2*xd2u + 4*sqrtS**2*xd2u + 4*uman*xd2u - 
     &                      2*m1**2*xd33u + 2*m2**2*xd33u + 2*sqrtS**2*xd33u + m1**2*xd3u + 3*m2**2*xd3u + m3**2*xd3u + 3*m4**2*xd3u - 2*sqrtS**2*xd3u + 2*uman*xd3u)))) + 
     &          g1Ru*(-2*g2Lu*M2Lpu*(g3Ru*M3Lpu*(-2*m1**4*xd2u - 14*m1**2*m2**2*xd2u + 4*m1**2*m3**2*xd2u - 14*m1**2*m4**2*xd2u - 6*m2**2*m4**2*xd2u + m1**2*sqrtS**2*xd2u - 
     &                   m2**2*sqrtS**2*xd2u - m3**2*sqrtS**2*xd2u + 5*m4**2*sqrtS**2*xd2u + sqrtS**4*xd2u - tman**2*xd2u + uman**2*xd2u + 2*m1**4*xd3u - 20*m1**2*m2**2*xd3u + 2*m2**4*xd3u - 
     &                   4*m1**2*m3**2*xd3u - 4*m2**2*m3**2*xd3u - 4*m1**2*m4**2*xd3u - 4*m2**2*m4**2*xd3u - 6*m1**2*sqrtS**2*xd3u - 6*m2**2*sqrtS**2*xd3u + 2*m3**2*sqrtS**2*xd3u + 
     &                   2*m4**2*sqrtS**2*xd3u + 4*sqrtS**4*xd3u + tman*(-6*xc03u - 6*xc04u - 6*M0Lpu**2*xd0u + 6*m1**2*xd0u - 6*M1Lpu**2*xd0u + 9*m1**2*xd2u + m2**2*xd2u + m3**2*xd2u + 
     &                      m4**2*xd2u + 4*m1**2*xd3u + 4*m2**2*xd3u - 2*sqrtS**2*xd3u) + 
     &                   uman*(6*xc03u + 6*xc04u + 6*M0Lpu**2*xd0u - 6*m1**2*xd0u + 6*M1Lpu**2*xd0u + m1**2*xd2u + 5*m2**2*xd2u - m3**2*xd2u - m4**2*xd2u - 4*sqrtS**2*xd2u + 4*m1**2*xd3u + 
     &                      4*m2**2*xd3u - 2*sqrtS**2*xd3u)) + 2*g3Lu*M0Lpu*(m1**4*xd0u - 10*m1**2*m2**2*xd0u + m2**4*xd0u - 2*m1**2*m3**2*xd0u - 2*m2**2*m3**2*xd0u - 2*m1**2*m4**2*xd0u - 
     &                   2*m2**2*m4**2*xd0u - 3*m1**2*sqrtS**2*xd0u - 3*m2**2*sqrtS**2*xd0u + m3**2*sqrtS**2*xd0u + m4**2*sqrtS**2*xd0u + 2*sqrtS**4*xd0u + m1**4*xd1u - 10*m1**2*m2**2*xd1u + 
     &                   m2**4*xd1u - 2*m1**2*m3**2*xd1u - 2*m2**2*m3**2*xd1u - 2*m1**2*m4**2*xd1u - 2*m2**2*m4**2*xd1u - 3*m1**2*sqrtS**2*xd1u - 3*m2**2*sqrtS**2*xd1u + 
     &                   m3**2*sqrtS**2*xd1u + m4**2*sqrtS**2*xd1u + 2*sqrtS**4*xd1u - 10*m1**2*m2**2*xd2u + 2*m2**4*xd2u - 4*m2**2*m3**2*xd2u - 6*m1**2*m4**2*xd2u + 2*m2**2*m4**2*xd2u - 
     &                   2*m1**2*sqrtS**2*xd2u - 4*m2**2*sqrtS**2*xd2u + m3**2*sqrtS**2*xd2u + m4**2*sqrtS**2*xd2u + 2*sqrtS**4*xd2u - tman**2*xd2u + uman**2*xd2u + m1**4*xd3u - 
     &                   10*m1**2*m2**2*xd3u + m2**4*xd3u - 2*m1**2*m3**2*xd3u - 2*m2**2*m3**2*xd3u - 2*m1**2*m4**2*xd3u - 2*m2**2*m4**2*xd3u - 3*m1**2*sqrtS**2*xd3u - 
     &                   3*m2**2*sqrtS**2*xd3u + m3**2*sqrtS**2*xd3u + m4**2*sqrtS**2*xd3u + 2*sqrtS**4*xd3u + 
     &                   uman*(3*xc02u + 3*xc03u + 2*m1**2*xd0u + 3*M1Lpu**2*xd0u + 2*m2**2*xd0u + 3*M3Lpu**2*xd0u - 4*sqrtS**2*xd0u + 2*m1**2*xd1u + 2*m2**2*xd1u - sqrtS**2*xd1u + 
     &                      2*m1**2*xd2u - m3**2*xd2u - m4**2*xd2u - 3*sqrtS**2*xd2u + 2*m1**2*xd3u + 2*m2**2*xd3u - sqrtS**2*xd3u) + 
     &                   tman*(-3*xc02u - 3*xc03u + 2*m1**2*xd0u - 3*M1Lpu**2*xd0u + 2*m2**2*xd0u - 3*M3Lpu**2*xd0u + 2*sqrtS**2*xd0u + 2*m1**2*xd1u + 2*m2**2*xd1u - sqrtS**2*xd1u + 
     &                      4*m1**2*xd2u + 2*m2**2*xd2u + m3**2*xd2u + m4**2*xd2u + sqrtS**2*xd2u + 2*m1**2*xd3u + 2*m2**2*xd3u - sqrtS**2*xd3u))) + 
     &             g2Ru*(2*g3Lu*M0Lpu*M3Lpu*(8*m1**2*m3**2*xd0u + 4*m2**2*m3**2*xd0u + 4*m1**2*m4**2*xd0u - 3*m1**2*sqrtS**2*xd0u - m2**2*sqrtS**2*xd0u - 3*m3**2*sqrtS**2*xd0u - 
     &                   m4**2*sqrtS**2*xd0u + sqrtS**4*xd0u + 6*m1**2*m3**2*xd11u + 2*m2**2*m3**2*xd11u - 6*m1**2*m4**2*xd11u - 2*m2**2*m4**2*xd11u - 2*m3**2*sqrtS**2*xd11u + 
     &                   2*m4**2*sqrtS**2*xd11u + 16*m1**2*m3**2*xd12u + 8*m2**2*m3**2*xd12u + 8*m1**2*m4**2*xd12u - 6*m1**2*sqrtS**2*xd12u - 2*m2**2*sqrtS**2*xd12u - 
     &                   6*m3**2*sqrtS**2*xd12u - 2*m4**2*sqrtS**2*xd12u + 2*sqrtS**4*xd12u - 4*m1**2*m3**2*xd13u + 4*m2**2*m3**2*xd13u + 4*m1**2*m4**2*xd13u - 4*m2**2*m4**2*xd13u + 
     &                   14*m1**2*m3**2*xd1u + 6*m2**2*m3**2*xd1u - 2*m1**2*m4**2*xd1u - 2*m2**2*m4**2*xd1u - 3*m1**2*sqrtS**2*xd1u - m2**2*sqrtS**2*xd1u - 5*m3**2*sqrtS**2*xd1u + 
     &                   m4**2*sqrtS**2*xd1u + sqrtS**4*xd1u + 8*m1**2*m3**2*xd22u + 8*m2**2*m3**2*xd22u + 8*m1**2*m4**2*xd22u + 8*m2**2*m4**2*xd22u - 4*m1**2*sqrtS**2*xd22u - 
     &                   4*m2**2*sqrtS**2*xd22u - 4*m3**2*sqrtS**2*xd22u - 4*m4**2*sqrtS**2*xd22u + 2*sqrtS**4*xd22u + 8*m2**2*m3**2*xd23u + 8*m1**2*m4**2*xd23u + 16*m2**2*m4**2*xd23u - 
     &                   2*m1**2*sqrtS**2*xd23u - 6*m2**2*sqrtS**2*xd23u - 2*m3**2*sqrtS**2*xd23u - 6*m4**2*sqrtS**2*xd23u + 2*sqrtS**4*xd23u + 16*m1**2*m3**2*xd2u + 12*m2**2*m3**2*xd2u + 
     &                   12*m1**2*m4**2*xd2u + 8*m2**2*m4**2*xd2u - 7*m1**2*sqrtS**2*xd2u - 5*m2**2*sqrtS**2*xd2u - 7*m3**2*sqrtS**2*xd2u - 5*m4**2*sqrtS**2*xd2u + 3*sqrtS**4*xd2u - 
     &                   2*m1**2*m3**2*xd33u - 6*m2**2*m3**2*xd33u + 2*m1**2*m4**2*xd33u + 6*m2**2*m4**2*xd33u + 2*m3**2*sqrtS**2*xd33u - 2*m4**2*sqrtS**2*xd33u + 
     &                   tman**2*(xd0u + 2*xd12u + xd1u + 2*xd22u - 6*xd23u + 3*xd2u - 3*xd3u) - 2*m1**2*m3**2*xd3u + 6*m2**2*m3**2*xd3u + 6*m1**2*m4**2*xd3u + 6*m2**2*m4**2*xd3u - 
     &                   m1**2*sqrtS**2*xd3u - 3*m2**2*sqrtS**2*xd3u - m3**2*sqrtS**2*xd3u - 3*m4**2*sqrtS**2*xd3u + sqrtS**4*xd3u + 
     &                   uman**2*(-3*xd0u - 6*xd12u - 3*xd1u - 6*xd22u + 2*xd23u - 9*xd2u + xd3u) + 
     &                   uman*(6*xc01u + 6*xc03u - 24*xd00u + m1**2*xd0u + 6*M1Lpu**2*xd0u + 3*m2**2*xd0u + 6*M2Lpu**2*xd0u + m3**2*xd0u - 3*m4**2*xd0u - 2*sqrtS**2*xd0u - 2*m1**2*xd11u + 
     &                      2*m2**2*xd11u - 2*sqrtS**2*xd11u + 2*m1**2*xd12u + 6*m2**2*xd12u + 2*m3**2*xd12u + 6*m4**2*xd12u - 4*sqrtS**2*xd12u + 4*sqrtS**2*xd13u - m1**2*xd1u + 
     &                      5*m2**2*xd1u + m3**2*xd1u + 3*m4**2*xd1u - 4*sqrtS**2*xd1u + 4*m1**2*xd22u + 4*m2**2*xd22u + 4*m3**2*xd22u + 4*m4**2*xd22u - 4*sqrtS**2*xd22u - 2*m1**2*xd23u - 
     &                      6*m2**2*xd23u - 2*m3**2*xd23u - 6*m4**2*xd23u + 4*sqrtS**2*xd23u + 5*m1**2*xd2u + 7*m2**2*xd2u + 5*m3**2*xd2u + 7*m4**2*xd2u - 6*sqrtS**2*xd2u + 2*m1**2*xd33u - 
     &                      2*m2**2*xd33u - 2*sqrtS**2*xd33u - m1**2*xd3u - 3*m2**2*xd3u - m3**2*xd3u - 3*m4**2*xd3u + 4*sqrtS**2*xd3u) + 
     &                   tman*(-6*xc01u - 6*xc03u + 24*xd00u - m1**2*xd0u - 6*M1Lpu**2*xd0u - 3*m2**2*xd0u - 6*M2Lpu**2*xd0u - m3**2*xd0u + 3*m4**2*xd0u + 2*sqrtS**2*xd0u + 2*uman*xd0u + 
     &                      2*m1**2*xd11u - 2*m2**2*xd11u + 2*sqrtS**2*xd11u - 2*m1**2*xd12u - 6*m2**2*xd12u - 2*m3**2*xd12u - 6*m4**2*xd12u + 4*sqrtS**2*xd12u + 4*uman*xd12u - 
     &                      4*sqrtS**2*xd13u + m1**2*xd1u - 5*m2**2*xd1u - m3**2*xd1u - 3*m4**2*xd1u + 4*sqrtS**2*xd1u + 2*uman*xd1u - 4*m1**2*xd22u - 4*m2**2*xd22u - 4*m3**2*xd22u - 
     &                      4*m4**2*xd22u + 4*sqrtS**2*xd22u + 4*uman*xd22u + 2*m1**2*xd23u + 6*m2**2*xd23u + 2*m3**2*xd23u + 6*m4**2*xd23u - 4*sqrtS**2*xd23u + 4*uman*xd23u - 5*m1**2*xd2u - 
     &                      7*m2**2*xd2u - 5*m3**2*xd2u - 7*m4**2*xd2u + 6*sqrtS**2*xd2u + 6*uman*xd2u - 2*m1**2*xd33u + 2*m2**2*xd33u + 2*sqrtS**2*xd33u + m1**2*xd3u + 3*m2**2*xd3u + 
     &                      m3**2*xd3u + 3*m4**2*xd3u - 4*sqrtS**2*xd3u + 2*uman*xd3u)) + 
     &                g3Ru*(-2*m1**4*xc01u + 20*m1**2*m2**2*xc01u - 2*m2**4*xc01u + 4*m1**2*m3**2*xc01u + 4*m2**2*m3**2*xc01u + 4*m1**2*m4**2*xc01u + 4*m2**2*m4**2*xc01u + 
     &                   6*m1**2*sqrtS**2*xc01u + 6*m2**2*sqrtS**2*xc01u - 2*m3**2*sqrtS**2*xc01u - 2*m4**2*sqrtS**2*xc01u - 4*sqrtS**4*xc01u + 8*m1**2*m3**2*xc02u + 4*m2**2*m3**2*xc02u + 
     &                   4*m1**2*m4**2*xc02u - 3*m1**2*sqrtS**2*xc02u - m2**2*sqrtS**2*xc02u - 3*m3**2*sqrtS**2*xc02u - m4**2*sqrtS**2*xc02u + sqrtS**4*xc02u + 4*m2**2*m3**2*xc03u + 
     &                   4*m1**2*m4**2*xc03u + 8*m2**2*m4**2*xc03u - m1**2*sqrtS**2*xc03u - 3*m2**2*sqrtS**2*xc03u - m3**2*sqrtS**2*xc03u - 3*m4**2*sqrtS**2*xc03u + sqrtS**4*xc03u - 
     &                   8*m1**2*m4**2*xc04u + 8*m2**2*m4**2*xc04u + 2*m1**2*sqrtS**2*xc04u - 2*m2**2*sqrtS**2*xc04u - 4*m1**2*m3**2*xc111u + 4*m2**2*m3**2*xc111u - 12*m1**2*m4**2*xc111u + 
     &                   12*m2**2*m4**2*xc111u + 4*m1**2*sqrtS**2*xc111u - 4*m2**2*sqrtS**2*xc111u + 2*m1**4*xc11u - 2*m2**4*xc11u - 8*m1**2*m3**2*xc11u + 8*m2**2*m3**2*xc11u - 
     &                   12*m1**2*m4**2*xc11u + 12*m2**2*m4**2*xc11u + 4*m1**2*sqrtS**2*xc11u - 4*m2**2*sqrtS**2*xc11u - 32*m1**2*m4**2*xc121u + 32*m2**2*m4**2*xc121u + 
     &                   8*m1**2*sqrtS**2*xc121u - 8*m2**2*sqrtS**2*xc121u + 6*m1**2*m3**2*xc12u + 2*m2**2*m3**2*xc12u - 6*m1**2*m4**2*xc12u - 2*m2**2*m4**2*xc12u - 2*m3**2*sqrtS**2*xc12u + 
     &                   2*m4**2*sqrtS**2*xc12u + 2*m1**4*xc13u - 6*m1**2*m2**2*xc13u + 4*m2**4*xc13u - 4*m1**2*m3**2*xc13u - 4*m2**2*m3**2*xc13u + 6*m1**2*m4**2*xc13u + 
     &                   18*m2**2*m4**2*xc13u - 6*m1**2*sqrtS**2*xc13u - 10*m2**2*sqrtS**2*xc13u + 2*m3**2*sqrtS**2*xc13u - 6*m4**2*sqrtS**2*xc13u + 4*sqrtS**4*xc13u + 8*m1**2*m3**2*xc21u - 
     &                   8*m2**2*m3**2*xc21u - 24*m1**2*m4**2*xc21u + 24*m2**2*m4**2*xc21u + 4*m1**2*sqrtS**2*xc21u - 4*m2**2*sqrtS**2*xc21u + 16*m1**2*m3**2*xc221u - 16*m2**2*m3**2*xc221u - 
     &                   16*m1**2*m4**2*xc221u + 16*m2**2*m4**2*xc221u - 2*m1**4*xc22u - 14*m1**2*m2**2*xc22u + 12*m1**2*m3**2*xc22u + 4*m2**2*m3**2*xc22u - 10*m1**2*m4**2*xc22u - 
     &                   6*m2**2*m4**2*xc22u - 2*m1**2*sqrtS**2*xc22u - 2*m2**2*sqrtS**2*xc22u - 4*m3**2*sqrtS**2*xc22u + 4*m4**2*sqrtS**2*xc22u + 2*sqrtS**4*xc22u - 2*m1**2*m3**2*xc23u - 
     &                   6*m2**2*m3**2*xc23u + 2*m1**2*m4**2*xc23u + 6*m2**2*m4**2*xc23u + 2*m3**2*sqrtS**2*xc23u - 2*m4**2*sqrtS**2*xc23u - 2*M0Lpu**2*m1**4*xd0u + 
     &                   20*M0Lpu**2*m1**2*m2**2*xd0u - 2*M0Lpu**2*m2**4*xd0u - 2*m1**4*M2Lpu**2*xd0u + 20*m1**2*m2**2*M2Lpu**2*xd0u - 2*m2**4*M2Lpu**2*xd0u + 12*M0Lpu**2*m1**2*m3**2*xd0u + 
     &                   12*M0Lpu**2*m2**2*m3**2*xd0u - 12*m1**2*m2**2*m3**2*xd0u + 4*M1Lpu**2*m2**2*m3**2*xd0u - 4*m2**4*m3**2*xd0u + 4*m1**2*M2Lpu**2*m3**2*xd0u + 
     &                   4*m2**2*M2Lpu**2*m3**2*xd0u + 8*m1**2*m3**2*M3Lpu**2*xd0u + 4*m2**2*m3**2*M3Lpu**2*xd0u + 12*M0Lpu**2*m1**2*m4**2*xd0u - 4*m1**4*m4**2*xd0u + 
     &                   4*m1**2*M1Lpu**2*m4**2*xd0u + 12*M0Lpu**2*m2**2*m4**2*xd0u - 12*m1**2*m2**2*m4**2*xd0u + 8*M1Lpu**2*m2**2*m4**2*xd0u + 4*m1**2*M2Lpu**2*m4**2*xd0u + 
     &                   4*m2**2*M2Lpu**2*m4**2*xd0u + 4*m1**2*M3Lpu**2*m4**2*xd0u + 2*M0Lpu**2*m1**2*sqrtS**2*xd0u + m1**4*sqrtS**2*xd0u - m1**2*M1Lpu**2*sqrtS**2*xd0u + 
     &                   2*M0Lpu**2*m2**2*sqrtS**2*xd0u + 6*m1**2*m2**2*sqrtS**2*xd0u - 3*M1Lpu**2*m2**2*sqrtS**2*xd0u + m2**4*sqrtS**2*xd0u + 6*m1**2*M2Lpu**2*sqrtS**2*xd0u + 
     &                   6*m2**2*M2Lpu**2*sqrtS**2*xd0u - 6*M0Lpu**2*m3**2*sqrtS**2*xd0u + m1**2*m3**2*sqrtS**2*xd0u - M1Lpu**2*m3**2*sqrtS**2*xd0u + 3*m2**2*m3**2*sqrtS**2*xd0u - 
     &                   2*M2Lpu**2*m3**2*sqrtS**2*xd0u - 3*m1**2*M3Lpu**2*sqrtS**2*xd0u - m2**2*M3Lpu**2*sqrtS**2*xd0u - 3*m3**2*M3Lpu**2*sqrtS**2*xd0u - 6*M0Lpu**2*m4**2*sqrtS**2*xd0u + 
     &                   3*m1**2*m4**2*sqrtS**2*xd0u - 3*M1Lpu**2*m4**2*sqrtS**2*xd0u + m2**2*m4**2*sqrtS**2*xd0u - 2*M2Lpu**2*m4**2*sqrtS**2*xd0u - M3Lpu**2*m4**2*sqrtS**2*xd0u - 
     &                   2*M0Lpu**2*sqrtS**4*xd0u - m1**2*sqrtS**4*xd0u + M1Lpu**2*sqrtS**4*xd0u - m2**2*sqrtS**4*xd0u - 4*M2Lpu**2*sqrtS**4*xd0u + M3Lpu**2*sqrtS**4*xd0u + 
     &                   12*M0Lpu**2*m1**2*m3**2*xd11u - 6*m1**4*m3**2*xd11u + 4*M0Lpu**2*m2**2*m3**2*xd11u - 8*m1**2*m2**2*m3**2*xd11u - 2*m2**4*m3**2*xd11u - 
     &                   12*M0Lpu**2*m1**2*m4**2*xd11u + 6*m1**4*m4**2*xd11u - 4*M0Lpu**2*m2**2*m4**2*xd11u + 8*m1**2*m2**2*m4**2*xd11u + 2*m2**4*m4**2*xd11u - 
     &                   4*M0Lpu**2*m3**2*sqrtS**2*xd11u + 8*m1**2*m3**2*sqrtS**2*xd11u + 4*m2**2*m3**2*sqrtS**2*xd11u + 4*M0Lpu**2*m4**2*sqrtS**2*xd11u - 8*m1**2*m4**2*sqrtS**2*xd11u - 
     &                   4*m2**2*m4**2*sqrtS**2*xd11u - 2*m3**2*sqrtS**4*xd11u + 2*m4**2*sqrtS**4*xd11u + 32*M0Lpu**2*m1**2*m3**2*xd12u - 16*m1**4*m3**2*xd12u + 
     &                   16*M0Lpu**2*m2**2*m3**2*xd12u - 24*m1**2*m2**2*m3**2*xd12u - 8*m2**4*m3**2*xd12u + 16*M0Lpu**2*m1**2*m4**2*xd12u - 8*m1**4*m4**2*xd12u - 8*m1**2*m2**2*m4**2*xd12u - 
     &                   12*M0Lpu**2*m1**2*sqrtS**2*xd12u + 6*m1**4*sqrtS**2*xd12u - 4*M0Lpu**2*m2**2*sqrtS**2*xd12u + 8*m1**2*m2**2*sqrtS**2*xd12u + 2*m2**4*sqrtS**2*xd12u - 
     &                   12*M0Lpu**2*m3**2*sqrtS**2*xd12u + 22*m1**2*m3**2*sqrtS**2*xd12u + 14*m2**2*m3**2*sqrtS**2*xd12u - 4*M0Lpu**2*m4**2*sqrtS**2*xd12u + 10*m1**2*m4**2*sqrtS**2*xd12u + 
     &                   2*m2**2*m4**2*sqrtS**2*xd12u + 4*M0Lpu**2*sqrtS**4*xd12u - 8*m1**2*sqrtS**4*xd12u - 4*m2**2*sqrtS**4*xd12u - 6*m3**2*sqrtS**4*xd12u - 2*m4**2*sqrtS**4*xd12u + 
     &                   2*sqrtS**6*xd12u - 8*M0Lpu**2*m1**2*m3**2*xd13u + 4*m1**4*m3**2*xd13u + 8*M0Lpu**2*m2**2*m3**2*xd13u - 4*m2**4*m3**2*xd13u + 8*M0Lpu**2*m1**2*m4**2*xd13u - 
     &                   4*m1**4*m4**2*xd13u - 8*M0Lpu**2*m2**2*m4**2*xd13u + 4*m2**4*m4**2*xd13u - 4*m1**2*m3**2*sqrtS**2*xd13u + 4*m2**2*m3**2*sqrtS**2*xd13u + 
     &                   4*m1**2*m4**2*sqrtS**2*xd13u - 4*m2**2*m4**2*sqrtS**2*xd13u - 2*M0Lpu**2*m1**4*xd1u + 2*m1**4*M1Lpu**2*xd1u + 20*M0Lpu**2*m1**2*m2**2*xd1u - 
     &                   20*m1**2*M1Lpu**2*m2**2*xd1u - 2*M0Lpu**2*m2**4*xd1u + 2*M1Lpu**2*m2**4*xd1u + 24*M0Lpu**2*m1**2*m3**2*xd1u - 6*m1**4*m3**2*xd1u - 6*m1**2*M1Lpu**2*m3**2*xd1u + 
     &                   16*M0Lpu**2*m2**2*m3**2*xd1u - 20*m1**2*m2**2*m3**2*xd1u - 2*M1Lpu**2*m2**2*m3**2*xd1u - 6*m2**4*m3**2*xd1u + 6*m1**2*m3**2*M3Lpu**2*xd1u + 
     &                   2*m2**2*m3**2*M3Lpu**2*xd1u + 8*M0Lpu**2*m1**2*m4**2*xd1u - 8*m1**4*m4**2*xd1u - 2*m1**2*M1Lpu**2*m4**2*xd1u + 24*m1**2*m2**2*m4**2*xd1u - 
     &                   6*M1Lpu**2*m2**2*m4**2*xd1u + 4*m1**2*m3**2*m4**2*xd1u + 4*m2**2*m3**2*m4**2*xd1u - 6*m1**2*M3Lpu**2*m4**2*xd1u - 2*m2**2*M3Lpu**2*m4**2*xd1u + 4*m1**2*m4**4*xd1u + 
     &                   4*m2**2*m4**4*xd1u + 3*m1**4*sqrtS**2*xd1u - 6*m1**2*M1Lpu**2*sqrtS**2*xd1u + 4*M0Lpu**2*m2**2*sqrtS**2*xd1u + 4*m1**2*m2**2*sqrtS**2*xd1u - 
     &                   6*M1Lpu**2*m2**2*sqrtS**2*xd1u + m2**4*sqrtS**2*xd1u - 10*M0Lpu**2*m3**2*sqrtS**2*xd1u + 11*m1**2*m3**2*sqrtS**2*xd1u + 2*M1Lpu**2*m3**2*sqrtS**2*xd1u + 
     &                   9*m2**2*m3**2*sqrtS**2*xd1u - 2*m3**2*M3Lpu**2*sqrtS**2*xd1u - 2*M0Lpu**2*m4**2*sqrtS**2*xd1u + 11*m1**2*m4**2*sqrtS**2*xd1u + 2*M1Lpu**2*m4**2*sqrtS**2*xd1u + 
     &                   5*m2**2*m4**2*sqrtS**2*xd1u - 2*m3**2*m4**2*sqrtS**2*xd1u + 2*M3Lpu**2*m4**2*sqrtS**2*xd1u - 2*m4**4*sqrtS**2*xd1u - 2*M0Lpu**2*sqrtS**4*xd1u - 
     &                   4*m1**2*sqrtS**4*xd1u + 4*M1Lpu**2*sqrtS**4*xd1u - 2*m2**2*sqrtS**4*xd1u - 3*m3**2*sqrtS**4*xd1u - 5*m4**2*sqrtS**4*xd1u + sqrtS**6*xd1u + 
     &                   16*M0Lpu**2*m1**2*m3**2*xd22u - 8*m1**4*m3**2*xd22u + 16*M0Lpu**2*m2**2*m3**2*xd22u - 16*m1**2*m2**2*m3**2*xd22u - 8*m2**4*m3**2*xd22u + 
     &                   16*M0Lpu**2*m1**2*m4**2*xd22u - 8*m1**4*m4**2*xd22u + 16*M0Lpu**2*m2**2*m4**2*xd22u - 16*m1**2*m2**2*m4**2*xd22u - 8*m2**4*m4**2*xd22u - 
     &                   8*M0Lpu**2*m1**2*sqrtS**2*xd22u + 4*m1**4*sqrtS**2*xd22u - 8*M0Lpu**2*m2**2*sqrtS**2*xd22u + 8*m1**2*m2**2*sqrtS**2*xd22u + 4*m2**4*sqrtS**2*xd22u - 
     &                   8*M0Lpu**2*m3**2*sqrtS**2*xd22u + 12*m1**2*m3**2*sqrtS**2*xd22u + 12*m2**2*m3**2*sqrtS**2*xd22u - 8*M0Lpu**2*m4**2*sqrtS**2*xd22u + 12*m1**2*m4**2*sqrtS**2*xd22u + 
     &                   12*m2**2*m4**2*sqrtS**2*xd22u + 4*M0Lpu**2*sqrtS**4*xd22u - 6*m1**2*sqrtS**4*xd22u - 6*m2**2*sqrtS**4*xd22u - 4*m3**2*sqrtS**4*xd22u - 4*m4**2*sqrtS**4*xd22u + 
     &                   2*sqrtS**6*xd22u + 16*M0Lpu**2*m2**2*m3**2*xd23u - 8*m1**2*m2**2*m3**2*xd23u - 8*m2**4*m3**2*xd23u + 16*M0Lpu**2*m1**2*m4**2*xd23u - 8*m1**4*m4**2*xd23u + 
     &                   32*M0Lpu**2*m2**2*m4**2*xd23u - 24*m1**2*m2**2*m4**2*xd23u - 16*m2**4*m4**2*xd23u - 4*M0Lpu**2*m1**2*sqrtS**2*xd23u + 2*m1**4*sqrtS**2*xd23u - 
     &                   12*M0Lpu**2*m2**2*sqrtS**2*xd23u + 8*m1**2*m2**2*sqrtS**2*xd23u + 6*m2**4*sqrtS**2*xd23u - 4*M0Lpu**2*m3**2*sqrtS**2*xd23u + 2*m1**2*m3**2*sqrtS**2*xd23u + 
     &                   10*m2**2*m3**2*sqrtS**2*xd23u - 12*M0Lpu**2*m4**2*sqrtS**2*xd23u + 14*m1**2*m4**2*sqrtS**2*xd23u + 22*m2**2*m4**2*sqrtS**2*xd23u + 4*M0Lpu**2*sqrtS**4*xd23u - 
     &                   4*m1**2*sqrtS**4*xd23u - 8*m2**2*sqrtS**4*xd23u - 2*m3**2*sqrtS**4*xd23u - 6*m4**2*sqrtS**4*xd23u + 2*sqrtS**6*xd23u + 2*m1**4*M1Lpu**2*xd2u + 
     &                   20*M0Lpu**2*m1**2*m2**2*xd2u - 6*m1**2*M1Lpu**2*m2**2*xd2u - 4*M0Lpu**2*m2**4*xd2u + 4*M1Lpu**2*m2**4*xd2u + 24*M0Lpu**2*m1**2*m3**2*xd2u - 6*m1**4*m3**2*xd2u - 
     &                   4*m1**2*M1Lpu**2*m3**2*xd2u + 32*M0Lpu**2*m2**2*m3**2*xd2u - 14*m1**2*m2**2*m3**2*xd2u - 4*M1Lpu**2*m2**2*m3**2*xd2u - 12*m2**4*m3**2*xd2u - 4*m1**2*m3**4*xd2u - 
     &                   2*m1**4*M3Lpu**2*xd2u - 14*m1**2*m2**2*M3Lpu**2*xd2u + 12*m1**2*m3**2*M3Lpu**2*xd2u + 4*m2**2*m3**2*M3Lpu**2*xd2u + 36*M0Lpu**2*m1**2*m4**2*xd2u - 
     &                   14*m1**4*m4**2*xd2u + 6*m1**2*M1Lpu**2*m4**2*xd2u + 20*M0Lpu**2*m2**2*m4**2*xd2u - 22*m1**2*m2**2*m4**2*xd2u + 18*M1Lpu**2*m2**2*m4**2*xd2u - 12*m2**4*m4**2*xd2u + 
     &                   18*m1**2*m3**2*m4**2*xd2u + 14*m2**2*m3**2*m4**2*xd2u - 10*m1**2*M3Lpu**2*m4**2*xd2u - 6*m2**2*M3Lpu**2*m4**2*xd2u - 2*m1**2*m4**4*xd2u - 10*m2**2*m4**4*xd2u - 
     &                   8*M0Lpu**2*m1**2*sqrtS**2*xd2u + 5*m1**4*sqrtS**2*xd2u - 6*m1**2*M1Lpu**2*sqrtS**2*xd2u - 4*M0Lpu**2*m2**2*sqrtS**2*xd2u + 14*m1**2*m2**2*sqrtS**2*xd2u - 
     &                   10*M1Lpu**2*m2**2*sqrtS**2*xd2u + 5*m2**4*sqrtS**2*xd2u - 14*M0Lpu**2*m3**2*sqrtS**2*xd2u + 12*m1**2*m3**2*sqrtS**2*xd2u + 2*M1Lpu**2*m3**2*sqrtS**2*xd2u + 
     &                   16*m2**2*m3**2*sqrtS**2*xd2u + m3**4*sqrtS**2*xd2u - 2*m1**2*M3Lpu**2*sqrtS**2*xd2u - 2*m2**2*M3Lpu**2*sqrtS**2*xd2u - 4*m3**2*M3Lpu**2*sqrtS**2*xd2u - 
     &                   14*M0Lpu**2*m4**2*sqrtS**2*xd2u + 20*m1**2*m4**2*sqrtS**2*xd2u - 6*M1Lpu**2*m4**2*sqrtS**2*xd2u + 20*m2**2*m4**2*sqrtS**2*xd2u - 8*m3**2*m4**2*sqrtS**2*xd2u + 
     &                   4*M3Lpu**2*m4**2*sqrtS**2*xd2u + 3*m4**4*sqrtS**2*xd2u + 2*M0Lpu**2*sqrtS**4*xd2u - 7*m1**2*sqrtS**4*xd2u + 4*M1Lpu**2*sqrtS**4*xd2u - 7*m2**2*sqrtS**4*xd2u - 
     &                   5*m3**2*sqrtS**4*xd2u + 2*M3Lpu**2*sqrtS**4*xd2u - 7*m4**2*sqrtS**4*xd2u + 2*sqrtS**6*xd2u + 2*uman**3*xd2u - 4*M0Lpu**2*m1**2*m3**2*xd33u + 2*m1**4*m3**2*xd33u - 
     &                   12*M0Lpu**2*m2**2*m3**2*xd33u + 8*m1**2*m2**2*m3**2*xd33u + 6*m2**4*m3**2*xd33u + 4*M0Lpu**2*m1**2*m4**2*xd33u - 2*m1**4*m4**2*xd33u + 
     &                   12*M0Lpu**2*m2**2*m4**2*xd33u - 8*m1**2*m2**2*m4**2*xd33u - 6*m2**4*m4**2*xd33u + 4*M0Lpu**2*m3**2*sqrtS**2*xd33u - 4*m1**2*m3**2*sqrtS**2*xd33u - 
     &                   8*m2**2*m3**2*sqrtS**2*xd33u - 4*M0Lpu**2*m4**2*sqrtS**2*xd33u + 4*m1**2*m4**2*sqrtS**2*xd33u + 8*m2**2*m4**2*sqrtS**2*xd33u + 2*m3**2*sqrtS**4*xd33u - 
     &                   2*m4**2*sqrtS**4*xd33u - 2*M0Lpu**2*m1**4*xd3u + 20*M0Lpu**2*m1**2*m2**2*xd3u - 2*M0Lpu**2*m2**4*xd3u - 2*m1**2*M1Lpu**2*m3**2*xd3u + 8*M0Lpu**2*m2**2*m3**2*xd3u + 
     &                   24*m1**2*m2**2*m3**2*xd3u - 6*M1Lpu**2*m2**2*m3**2*xd3u - 8*m2**4*m3**2*xd3u + 4*m1**2*m3**4*xd3u + 4*m2**2*m3**4*xd3u + 2*m1**4*M3Lpu**2*xd3u - 
     &                   20*m1**2*m2**2*M3Lpu**2*xd3u + 2*m2**4*M3Lpu**2*xd3u - 6*m1**2*m3**2*M3Lpu**2*xd3u - 2*m2**2*m3**2*M3Lpu**2*xd3u + 16*M0Lpu**2*m1**2*m4**2*xd3u - 
     &                   6*m1**4*m4**2*xd3u + 2*m1**2*M1Lpu**2*m4**2*xd3u + 24*M0Lpu**2*m2**2*m4**2*xd3u - 20*m1**2*m2**2*m4**2*xd3u + 6*M1Lpu**2*m2**2*m4**2*xd3u - 6*m2**4*m4**2*xd3u + 
     &                   4*m1**2*m3**2*m4**2*xd3u + 4*m2**2*m3**2*m4**2*xd3u - 2*m1**2*M3Lpu**2*m4**2*xd3u - 6*m2**2*M3Lpu**2*m4**2*xd3u + 4*M0Lpu**2*m1**2*sqrtS**2*xd3u + 
     &                   m1**4*sqrtS**2*xd3u + 4*m1**2*m2**2*sqrtS**2*xd3u + 3*m2**4*sqrtS**2*xd3u - 2*M0Lpu**2*m3**2*sqrtS**2*xd3u + 5*m1**2*m3**2*sqrtS**2*xd3u + 
     &                   2*M1Lpu**2*m3**2*sqrtS**2*xd3u + 11*m2**2*m3**2*sqrtS**2*xd3u - 2*m3**4*sqrtS**2*xd3u - 6*m1**2*M3Lpu**2*sqrtS**2*xd3u - 6*m2**2*M3Lpu**2*sqrtS**2*xd3u + 
     &                   2*m3**2*M3Lpu**2*sqrtS**2*xd3u - 10*M0Lpu**2*m4**2*sqrtS**2*xd3u + 9*m1**2*m4**2*sqrtS**2*xd3u - 2*M1Lpu**2*m4**2*sqrtS**2*xd3u + 11*m2**2*m4**2*sqrtS**2*xd3u - 
     &                   2*m3**2*m4**2*sqrtS**2*xd3u + 2*M3Lpu**2*m4**2*sqrtS**2*xd3u - 2*M0Lpu**2*sqrtS**4*xd3u - 2*m1**2*sqrtS**4*xd3u - 4*m2**2*sqrtS**4*xd3u - 5*m3**2*sqrtS**4*xd3u + 
     &                   4*M3Lpu**2*sqrtS**4*xd3u - 3*m4**2*sqrtS**4*xd3u + sqrtS**6*xd3u + 
     &                   tman**2*(xc02u - 3*xc03u - 4*xc04u - 2*xc11u - 16*xc121u - 4*xc13u - 8*xc21u - 2*M0Lpu**2*xd0u + 3*m1**2*xd0u - 3*M1Lpu**2*xd0u - m2**2*xd0u + M3Lpu**2*xd0u + 
     &                      4*M0Lpu**2*xd12u - 2*m1**2*xd12u - 2*m2**2*xd12u + 2*sqrtS**2*xd12u + 2*M0Lpu**2*xd1u - m1**2*xd1u - m2**2*xd1u + sqrtS**2*xd1u + 4*M0Lpu**2*xd22u - 
     &                      2*m1**2*xd22u - 2*m2**2*xd22u + 2*sqrtS**2*xd22u - 12*M0Lpu**2*xd23u + 6*m1**2*xd23u + 6*m2**2*xd23u - 6*sqrtS**2*xd23u + 4*M0Lpu**2*xd2u + m1**2*xd2u - 
     &                      4*M1Lpu**2*xd2u - 3*m2**2*xd2u + m3**2*xd2u + m4**2*xd2u + 2*sqrtS**2*xd2u - 2*uman*xd2u - 6*M0Lpu**2*xd3u + 3*m1**2*xd3u + 3*m2**2*xd3u - 3*sqrtS**2*xd3u) + 
     &                   uman**2*(3*xc02u + 7*xc03u + 4*xc04u + 2*xc11u + 16*xc121u + 2*xc13u + 8*xc21u - 2*xc22u - 2*M0Lpu**2*xd0u + 3*m1**2*xd0u + 7*M1Lpu**2*xd0u + 7*m2**2*xd0u + 
     &                      3*M3Lpu**2*xd0u - 8*sqrtS**2*xd0u - 12*M0Lpu**2*xd12u + 6*m1**2*xd12u + 6*m2**2*xd12u - 6*sqrtS**2*xd12u - 6*M0Lpu**2*xd1u + 7*m1**2*xd1u + 7*m2**2*xd1u - 
     &                      5*sqrtS**2*xd1u - 12*M0Lpu**2*xd22u + 6*m1**2*xd22u + 6*m2**2*xd22u - 6*sqrtS**2*xd22u + 4*M0Lpu**2*xd23u - 2*m1**2*xd23u - 2*m2**2*xd23u + 2*sqrtS**2*xd23u - 
     &                      16*M0Lpu**2*xd2u + 9*m1**2*xd2u + 2*M1Lpu**2*xd2u + 9*m2**2*xd2u - 3*m3**2*xd2u - 2*M3Lpu**2*xd2u - 3*m4**2*xd2u - 12*sqrtS**2*xd2u + 2*M0Lpu**2*xd3u + 
     &                      3*m1**2*xd3u + 3*m2**2*xd3u - sqrtS**2*xd3u) + uman*
     &                    (12*xb01u + 12*xb02u - 48*xc001u + 12*M0Lpu**2*xc01u - 10*m1**2*xc01u - 10*m2**2*xc01u + 8*sqrtS**2*xc01u - 5*m1**2*xc02u + 12*M1Lpu**2*xc02u + 3*m2**2*xc02u + 
     &                      m3**2*xc02u - 3*m4**2*xc02u - 2*sqrtS**2*xc02u - m1**2*xc03u - 9*m2**2*xc03u - 7*m3**2*xc03u + 12*M3Lpu**2*xc03u - 3*m4**2*xc03u + 2*sqrtS**2*xc03u - 
     &                      4*m1**2*xc04u + 12*M1Lpu**2*xc04u - 4*m2**2*xc04u - 8*m3**2*xc04u + 6*sqrtS**2*xc04u + 4*m3**2*xc111u - 4*m4**2*xc111u - 4*sqrtS**2*xc111u - 4*m1**2*xc11u + 
     &                      4*m3**2*xc11u + 16*m4**2*xc11u - 16*m1**2*xc121u - 16*m2**2*xc121u - 8*m3**2*xc121u - 24*m4**2*xc121u - 2*m1**2*xc12u + 2*m2**2*xc12u - 2*sqrtS**2*xc12u + 
     &                      2*m1**2*xc13u - 8*m2**2*xc13u - 2*m3**2*xc13u - 4*m4**2*xc13u - 8*m1**2*xc21u - 8*m2**2*xc21u - 16*m3**2*xc21u + 4*sqrtS**2*xc21u - 16*sqrtS**2*xc221u + 
     &                      2*m1**2*xc22u + 8*m2**2*xc22u + 2*m4**2*xc22u - 6*sqrtS**2*xc22u + 2*m1**2*xc23u - 2*m2**2*xc23u - 2*sqrtS**2*xc23u - 48*M0Lpu**2*xd00u + 24*m1**2*xd00u + 
     &                      24*m2**2*xd00u - 24*sqrtS**2*xd00u - 4*M0Lpu**2*m1**2*xd0u + 3*m1**4*xd0u - m1**2*M1Lpu**2*xd0u - 4*M0Lpu**2*m2**2*xd0u - 18*m1**2*m2**2*xd0u - 
     &                      9*M1Lpu**2*m2**2*xd0u - m2**4*xd0u + 12*M0Lpu**2*M2Lpu**2*xd0u - 10*m1**2*M2Lpu**2*xd0u - 10*m2**2*M2Lpu**2*xd0u - 6*M0Lpu**2*m3**2*xd0u + 3*m1**2*m3**2*xd0u - 
     &                      7*M1Lpu**2*m3**2*xd0u - 5*m2**2*m3**2*xd0u - 5*m1**2*M3Lpu**2*xd0u + 12*M1Lpu**2*M3Lpu**2*xd0u + 3*m2**2*M3Lpu**2*xd0u + m3**2*M3Lpu**2*xd0u - 
     &                      6*M0Lpu**2*m4**2*xd0u - m1**2*m4**2*xd0u - 3*M1Lpu**2*m4**2*xd0u - m2**2*m4**2*xd0u - 3*M3Lpu**2*m4**2*xd0u + 8*M0Lpu**2*sqrtS**2*xd0u - 8*m1**2*sqrtS**2*xd0u + 
     &                      2*M1Lpu**2*sqrtS**2*xd0u - 4*m2**2*sqrtS**2*xd0u + 8*M2Lpu**2*sqrtS**2*xd0u + 2*m3**2*sqrtS**2*xd0u - 2*M3Lpu**2*sqrtS**2*xd0u + 2*m4**2*sqrtS**2*xd0u + 
     &                      4*sqrtS**4*xd0u - 4*M0Lpu**2*m1**2*xd11u + 2*m1**4*xd11u + 4*M0Lpu**2*m2**2*xd11u - 2*m2**4*xd11u - 4*M0Lpu**2*sqrtS**2*xd11u + 4*m2**2*sqrtS**2*xd11u - 
     &                      2*sqrtS**4*xd11u + 4*M0Lpu**2*m1**2*xd12u - 2*m1**4*xd12u + 12*M0Lpu**2*m2**2*xd12u - 8*m1**2*m2**2*xd12u - 6*m2**4*xd12u + 4*M0Lpu**2*m3**2*xd12u - 
     &                      2*m1**2*m3**2*xd12u - 2*m2**2*m3**2*xd12u + 12*M0Lpu**2*m4**2*xd12u - 6*m1**2*m4**2*xd12u - 6*m2**2*m4**2*xd12u - 8*M0Lpu**2*sqrtS**2*xd12u + 
     &                      6*m1**2*sqrtS**2*xd12u + 10*m2**2*sqrtS**2*xd12u + 2*m3**2*sqrtS**2*xd12u + 6*m4**2*sqrtS**2*xd12u - 4*sqrtS**4*xd12u + 8*M0Lpu**2*sqrtS**2*xd13u - 
     &                      4*m1**2*sqrtS**2*xd13u - 4*m2**2*sqrtS**2*xd13u + 4*sqrtS**4*xd13u - 4*M0Lpu**2*m1**2*xd1u + m1**4*xd1u + 4*m1**2*M1Lpu**2*xd1u + 4*M0Lpu**2*m2**2*xd1u - 
     &                      22*m1**2*m2**2*xd1u + 4*M1Lpu**2*m2**2*xd1u - 3*m2**4*xd1u + 2*M0Lpu**2*m3**2*xd1u - 5*m1**2*m3**2*xd1u - 5*m2**2*m3**2*xd1u - 2*m1**2*M3Lpu**2*xd1u + 
     &                      2*m2**2*M3Lpu**2*xd1u + 6*M0Lpu**2*m4**2*xd1u - 11*m1**2*m4**2*xd1u - 11*m2**2*m4**2*xd1u - 2*M0Lpu**2*sqrtS**2*xd1u - 5*m1**2*sqrtS**2*xd1u + 
     &                      m2**2*sqrtS**2*xd1u + 3*m3**2*sqrtS**2*xd1u - 2*M3Lpu**2*sqrtS**2*xd1u + 7*m4**2*sqrtS**2*xd1u + 2*sqrtS**4*xd1u + 8*M0Lpu**2*m1**2*xd22u - 4*m1**4*xd22u + 
     &                      8*M0Lpu**2*m2**2*xd22u - 8*m1**2*m2**2*xd22u - 4*m2**4*xd22u + 8*M0Lpu**2*m3**2*xd22u - 4*m1**2*m3**2*xd22u - 4*m2**2*m3**2*xd22u + 8*M0Lpu**2*m4**2*xd22u - 
     &                      4*m1**2*m4**2*xd22u - 4*m2**2*m4**2*xd22u - 8*M0Lpu**2*sqrtS**2*xd22u + 8*m1**2*sqrtS**2*xd22u + 8*m2**2*sqrtS**2*xd22u + 4*m3**2*sqrtS**2*xd22u + 
     &                      4*m4**2*sqrtS**2*xd22u - 4*sqrtS**4*xd22u - 4*M0Lpu**2*m1**2*xd23u + 2*m1**4*xd23u - 12*M0Lpu**2*m2**2*xd23u + 8*m1**2*m2**2*xd23u + 6*m2**4*xd23u - 
     &                      4*M0Lpu**2*m3**2*xd23u + 2*m1**2*m3**2*xd23u + 2*m2**2*m3**2*xd23u - 12*M0Lpu**2*m4**2*xd23u + 6*m1**2*m4**2*xd23u + 6*m2**2*m4**2*xd23u + 
     &                      8*M0Lpu**2*sqrtS**2*xd23u - 6*m1**2*sqrtS**2*xd23u - 10*m2**2*sqrtS**2*xd23u - 2*m3**2*sqrtS**2*xd23u - 6*m4**2*sqrtS**2*xd23u + 4*sqrtS**4*xd23u + 
     &                      4*M0Lpu**2*m1**2*xd2u - 3*m1**4*xd2u + 2*m1**2*M1Lpu**2*xd2u + 8*M0Lpu**2*m2**2*xd2u - 26*m1**2*m2**2*xd2u - 8*M1Lpu**2*m2**2*xd2u - 3*m2**4*xd2u + 
     &                      10*M0Lpu**2*m3**2*xd2u - 4*m1**2*m3**2*xd2u - 2*M1Lpu**2*m3**2*xd2u - 18*m2**2*m3**2*xd2u + m3**4*xd2u + 2*m1**2*M3Lpu**2*xd2u + 8*m2**2*M3Lpu**2*xd2u + 
     &                      10*M0Lpu**2*m4**2*xd2u - 16*m1**2*m4**2*xd2u - 4*M1Lpu**2*m4**2*xd2u + 2*m2**2*m4**2*xd2u + 2*m3**2*m4**2*xd2u + 2*M3Lpu**2*m4**2*xd2u + m4**4*xd2u - 
     &                      2*M0Lpu**2*sqrtS**2*xd2u + 2*m1**2*sqrtS**2*xd2u + 2*m2**2*sqrtS**2*xd2u + 10*m3**2*sqrtS**2*xd2u - 6*M3Lpu**2*sqrtS**2*xd2u + 8*m4**2*sqrtS**2*xd2u + 
     &                      4*M0Lpu**2*m1**2*xd33u - 2*m1**4*xd33u - 4*M0Lpu**2*m2**2*xd33u + 2*m2**4*xd33u - 4*M0Lpu**2*sqrtS**2*xd33u + 4*m1**2*sqrtS**2*xd33u - 2*sqrtS**4*xd33u - 
     &                      4*M0Lpu**2*m1**2*xd3u + m1**4*xd3u + 2*m1**2*M1Lpu**2*xd3u - 12*M0Lpu**2*m2**2*xd3u - 14*m1**2*m2**2*xd3u - 2*M1Lpu**2*m2**2*xd3u + 5*m2**4*xd3u - 
     &                      2*M0Lpu**2*m3**2*xd3u - 7*m1**2*m3**2*xd3u - 7*m2**2*m3**2*xd3u + 4*m1**2*M3Lpu**2*xd3u + 4*m2**2*M3Lpu**2*xd3u - 6*M0Lpu**2*m4**2*xd3u - m1**2*m4**2*xd3u - 
     &                      m2**2*m4**2*xd3u + 6*M0Lpu**2*sqrtS**2*xd3u - 7*m1**2*sqrtS**2*xd3u - 2*M1Lpu**2*sqrtS**2*xd3u - 13*m2**2*sqrtS**2*xd3u + 3*m3**2*sqrtS**2*xd3u - 
     &                      m4**2*sqrtS**2*xd3u + 6*sqrtS**4*xd3u) + tman*(-12*xb01u - 12*xb02u + 48*xc001u - 12*M0Lpu**2*xc01u + 2*m1**2*xc01u + 2*m2**2*xc01u - 4*sqrtS**2*xc01u + 
     &                      5*m1**2*xc02u - 12*M1Lpu**2*xc02u - 3*m2**2*xc02u - m3**2*xc02u + 3*m4**2*xc02u + 2*sqrtS**2*xc02u - 4*uman*xc02u + m1**2*xc03u + 9*m2**2*xc03u + 7*m3**2*xc03u - 
     &                      12*M3Lpu**2*xc03u + 3*m4**2*xc03u - 2*sqrtS**2*xc03u - 4*uman*xc03u + 4*m1**2*xc04u - 12*M1Lpu**2*xc04u + 4*m2**2*xc04u + 8*m3**2*xc04u - 6*sqrtS**2*xc04u - 
     &                      4*m3**2*xc111u + 4*m4**2*xc111u + 4*sqrtS**2*xc111u + 4*m2**2*xc11u - 4*m3**2*xc11u - 16*m4**2*xc11u + 16*m1**2*xc121u + 16*m2**2*xc121u + 8*m3**2*xc121u + 
     &                      24*m4**2*xc121u + 2*m1**2*xc12u - 2*m2**2*xc12u + 2*sqrtS**2*xc12u + 6*m2**2*xc13u + 2*m3**2*xc13u + 4*m4**2*xc13u + 2*uman*xc13u + 8*m1**2*xc21u + 
     &                      8*m2**2*xc21u + 16*m3**2*xc21u - 4*sqrtS**2*xc21u + 16*sqrtS**2*xc221u + 8*m1**2*xc22u - 2*m2**2*xc22u - 2*m4**2*xc22u + 2*sqrtS**2*xc22u + 2*uman*xc22u - 
     &                      2*m1**2*xc23u + 2*m2**2*xc23u + 2*sqrtS**2*xc23u + 48*M0Lpu**2*xd00u - 24*m1**2*xd00u - 24*m2**2*xd00u + 24*sqrtS**2*xd00u - 4*M0Lpu**2*m1**2*xd0u - m1**4*xd0u + 
     &                      m1**2*M1Lpu**2*xd0u - 4*M0Lpu**2*m2**2*xd0u - 2*m1**2*m2**2*xd0u + 9*M1Lpu**2*m2**2*xd0u + 3*m2**4*xd0u - 12*M0Lpu**2*M2Lpu**2*xd0u + 2*m1**2*M2Lpu**2*xd0u + 
     &                      2*m2**2*M2Lpu**2*xd0u + 6*M0Lpu**2*m3**2*xd0u - 7*m1**2*m3**2*xd0u + 7*M1Lpu**2*m3**2*xd0u + m2**2*m3**2*xd0u + 5*m1**2*M3Lpu**2*xd0u - 
     &                      12*M1Lpu**2*M3Lpu**2*xd0u - 3*m2**2*M3Lpu**2*xd0u - m3**2*M3Lpu**2*xd0u + 6*M0Lpu**2*m4**2*xd0u - 3*m1**2*m4**2*xd0u + 3*M1Lpu**2*m4**2*xd0u - 
     &                      3*m2**2*m4**2*xd0u + 3*M3Lpu**2*m4**2*xd0u - 4*M0Lpu**2*sqrtS**2*xd0u + 2*m1**2*sqrtS**2*xd0u - 2*M1Lpu**2*sqrtS**2*xd0u - 2*m2**2*sqrtS**2*xd0u - 
     &                      4*M2Lpu**2*sqrtS**2*xd0u + 2*M3Lpu**2*sqrtS**2*xd0u + 4*M0Lpu**2*uman*xd0u + 2*m1**2*uman*xd0u - 4*M1Lpu**2*uman*xd0u + 2*m2**2*uman*xd0u - 4*M3Lpu**2*uman*xd0u + 
     &                      4*sqrtS**2*uman*xd0u + 4*M0Lpu**2*m1**2*xd11u - 2*m1**4*xd11u - 4*M0Lpu**2*m2**2*xd11u + 2*m2**4*xd11u + 4*M0Lpu**2*sqrtS**2*xd11u - 4*m2**2*sqrtS**2*xd11u + 
     &                      2*sqrtS**4*xd11u - 4*M0Lpu**2*m1**2*xd12u + 2*m1**4*xd12u - 12*M0Lpu**2*m2**2*xd12u + 8*m1**2*m2**2*xd12u + 6*m2**4*xd12u - 4*M0Lpu**2*m3**2*xd12u + 
     &                      2*m1**2*m3**2*xd12u + 2*m2**2*m3**2*xd12u - 12*M0Lpu**2*m4**2*xd12u + 6*m1**2*m4**2*xd12u + 6*m2**2*m4**2*xd12u + 8*M0Lpu**2*sqrtS**2*xd12u - 
     &                      6*m1**2*sqrtS**2*xd12u - 10*m2**2*sqrtS**2*xd12u - 2*m3**2*sqrtS**2*xd12u - 6*m4**2*sqrtS**2*xd12u + 4*sqrtS**4*xd12u + 8*M0Lpu**2*uman*xd12u - 
     &                      4*m1**2*uman*xd12u - 4*m2**2*uman*xd12u + 4*sqrtS**2*uman*xd12u - 8*M0Lpu**2*sqrtS**2*xd13u + 4*m1**2*sqrtS**2*xd13u + 4*m2**2*sqrtS**2*xd13u - 4*sqrtS**4*xd13u - 
     &                      4*M0Lpu**2*m1**2*xd1u + m1**4*xd1u + 4*m1**2*M1Lpu**2*xd1u - 12*M0Lpu**2*m2**2*xd1u + 2*m1**2*m2**2*xd1u + 4*M1Lpu**2*m2**2*xd1u + 5*m2**4*xd1u - 
     &                      2*M0Lpu**2*m3**2*xd1u + m1**2*m3**2*xd1u + m2**2*m3**2*xd1u + 2*m1**2*M3Lpu**2*xd1u - 2*m2**2*M3Lpu**2*xd1u - 6*M0Lpu**2*m4**2*xd1u - m1**2*m4**2*xd1u - 
     &                      m2**2*m4**2*xd1u + 6*M0Lpu**2*sqrtS**2*xd1u - m1**2*sqrtS**2*xd1u - 4*M1Lpu**2*sqrtS**2*xd1u - 7*m2**2*sqrtS**2*xd1u - m3**2*sqrtS**2*xd1u + 
     &                      2*M3Lpu**2*sqrtS**2*xd1u - m4**2*sqrtS**2*xd1u + 2*sqrtS**4*xd1u + 4*M0Lpu**2*uman*xd1u + 2*m1**2*uman*xd1u + 2*m2**2*uman*xd1u - 8*M0Lpu**2*m1**2*xd22u + 
     &                      4*m1**4*xd22u - 8*M0Lpu**2*m2**2*xd22u + 8*m1**2*m2**2*xd22u + 4*m2**4*xd22u - 8*M0Lpu**2*m3**2*xd22u + 4*m1**2*m3**2*xd22u + 4*m2**2*m3**2*xd22u - 
     &                      8*M0Lpu**2*m4**2*xd22u + 4*m1**2*m4**2*xd22u + 4*m2**2*m4**2*xd22u + 8*M0Lpu**2*sqrtS**2*xd22u - 8*m1**2*sqrtS**2*xd22u - 8*m2**2*sqrtS**2*xd22u - 
     &                      4*m3**2*sqrtS**2*xd22u - 4*m4**2*sqrtS**2*xd22u + 4*sqrtS**4*xd22u + 8*M0Lpu**2*uman*xd22u - 4*m1**2*uman*xd22u - 4*m2**2*uman*xd22u + 4*sqrtS**2*uman*xd22u + 
     &                      4*M0Lpu**2*m1**2*xd23u - 2*m1**4*xd23u + 12*M0Lpu**2*m2**2*xd23u - 8*m1**2*m2**2*xd23u - 6*m2**4*xd23u + 4*M0Lpu**2*m3**2*xd23u - 2*m1**2*m3**2*xd23u - 
     &                      2*m2**2*m3**2*xd23u + 12*M0Lpu**2*m4**2*xd23u - 6*m1**2*m4**2*xd23u - 6*m2**2*m4**2*xd23u - 8*M0Lpu**2*sqrtS**2*xd23u + 6*m1**2*sqrtS**2*xd23u + 
     &                      10*m2**2*sqrtS**2*xd23u + 2*m3**2*sqrtS**2*xd23u + 6*m4**2*sqrtS**2*xd23u - 4*sqrtS**4*xd23u + 8*M0Lpu**2*uman*xd23u - 4*m1**2*uman*xd23u - 4*m2**2*uman*xd23u + 
     &                      4*sqrtS**2*uman*xd23u - 16*M0Lpu**2*m1**2*xd2u + 3*m1**4*xd2u - 12*M0Lpu**2*m2**2*xd2u + 6*m1**2*m2**2*xd2u + 6*M1Lpu**2*m2**2*xd2u + 7*m2**4*xd2u - 
     &                      10*M0Lpu**2*m3**2*xd2u - 6*m1**2*m3**2*xd2u + 2*M1Lpu**2*m3**2*xd2u + 4*m2**2*m3**2*xd2u - m3**4*xd2u + 8*m1**2*M3Lpu**2*xd2u - 2*m2**2*M3Lpu**2*xd2u - 
     &                      10*M0Lpu**2*m4**2*xd2u + 2*m1**2*m4**2*xd2u + 4*M1Lpu**2*m4**2*xd2u + 4*m2**2*m4**2*xd2u - 2*m3**2*m4**2*xd2u - 2*M3Lpu**2*m4**2*xd2u - m4**4*xd2u + 
     &                      6*M0Lpu**2*sqrtS**2*xd2u - 6*m1**2*sqrtS**2*xd2u - 10*m2**2*sqrtS**2*xd2u - 4*m3**2*sqrtS**2*xd2u + 2*M3Lpu**2*sqrtS**2*xd2u - 6*m4**2*sqrtS**2*xd2u + 
     &                      4*sqrtS**4*xd2u + 12*M0Lpu**2*uman*xd2u + 2*m1**2*uman*xd2u + 2*M1Lpu**2*uman*xd2u - 2*m2**2*uman*xd2u + 2*m3**2*uman*xd2u + 2*M3Lpu**2*uman*xd2u + 
     &                      2*m4**2*uman*xd2u + 6*sqrtS**2*uman*xd2u - 4*M0Lpu**2*m1**2*xd33u + 2*m1**4*xd33u + 4*M0Lpu**2*m2**2*xd33u - 2*m2**4*xd33u + 4*M0Lpu**2*sqrtS**2*xd33u - 
     &                      4*m1**2*sqrtS**2*xd33u + 2*sqrtS**4*xd33u - 4*M0Lpu**2*m1**2*xd3u + m1**4*xd3u - 2*m1**2*M1Lpu**2*xd3u + 4*M0Lpu**2*m2**2*xd3u - 6*m1**2*m2**2*xd3u + 
     &                      2*M1Lpu**2*m2**2*xd3u - 3*m2**4*xd3u + 2*M0Lpu**2*m3**2*xd3u - 5*m1**2*m3**2*xd3u - 5*m2**2*m3**2*xd3u + 4*m1**2*M3Lpu**2*xd3u + 4*m2**2*M3Lpu**2*xd3u + 
     &                      6*M0Lpu**2*m4**2*xd3u - 3*m1**2*m4**2*xd3u - 3*m2**2*m4**2*xd3u - 2*M0Lpu**2*sqrtS**2*xd3u + m1**2*sqrtS**2*xd3u + 2*M1Lpu**2*sqrtS**2*xd3u + 
     &                      7*m2**2*sqrtS**2*xd3u + 3*m3**2*sqrtS**2*xd3u - 4*M3Lpu**2*sqrtS**2*xd3u + 3*m4**2*sqrtS**2*xd3u - 2*sqrtS**4*xd3u + 4*M0Lpu**2*uman*xd3u + 2*m1**2*uman*xd3u + 
     &                      2*m2**2*uman*xd3u)))))))/(64.*pi**2*propac)  )

      enddo           
    
!================================================================================================================================== 

      VV_TX(1) = Box_TQ
      VV_TX(2) = Box_TS
      VV_TX(3) = Box_TV
      VV_TX(4) = Box_TT
      VV_TX(5) = Box_TU

      VV_TX(6) = VV_TX(1) + VV_TX(2) + VV_TX(3) + VV_TX(4) + VV_TX(5)

      VV_UX(1) = Box_UQ
      VV_UX(2) = Box_US
      VV_UX(3) = Box_UV
      VV_UX(4) = Box_UT
      VV_UX(5) = Box_UU

      VV_UX(6) = VV_UX(1) + VV_UX(2) + VV_UX(3) + VV_UX(4) + VV_UX(5)

!==================================================================================================================================

      end
