***********************************************************************************************************************************

      subroutine QQ2xx_M2GluonBoxHH(argxh,argyh,HH_TX,HH_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"
#include "../util/DMNLO_LoopIntegrals.h"

      integer i,j,argxh,argyh
      complex*16 g0t,g1t,g2t,g3t,g0u,g1u,g2u,g3u,g0Lt,g0Rt,g1Lt,g1Rt,g2Lt,g2Rt,g3Lt,g3Rt,g0Lu,g0Ru,g1Lu,g1Ru,g2Lu,g2Ru,g3Lu,g3Ru,
     &           xb01t,xb02t,xc01t,xc02t,xc03t,xc04t,xc11t,xc21t,xd0t,xb01u,xb02u,xc01u,xc02u,xc03u,xc04u,xc11u,xc21u,xd0u,C0t,C0u
      double precision M0Lpt,M1Lpt,M2Lpt,M3Lpt,M0Lpu,M1Lpu,M2Lpu,M3Lpu,HH_TX(6),HH_UX(6),Box_TQ,Box_TS,Box_TV,Box_TT,Box_TU,Box_UQ,Box_US,Box_UV,Box_UT,Box_UU
        
      external myBget,myCget,myDget
  
! =================================================================================================================================	
!                             ---                     Higgs-Higgs-final-states                    ---
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
        g1t=HSfSft(argxh,Loopmeigen1,i,Looptype1,Looptypet,Loopgen1,Loopgent)
        g2t=HSfSfu(argyh,i,Loopmeigen2,Looptypet,Looptype2,Loopgent,Loopgen2)
        g3t=VSfSf(Looptype0,sfmeigen2,Loopmeigen2,sftype2,Looptype2,sfgen2,Loopgen2)

!=====================================================!Box(u)_Gluon================================================================
       
        M0Lpu=MV(Looptype0)
        M1Lpu=MSf(Loopmeigen1,Looptype1,Loopgen1)
        M2Lpu=MSf(i,Looptypeu,Loopgenu)
        M3Lpu=MSf(Loopmeigen2,Looptype2,Loopgen2)

        g0u=VSfSf(Looptype0,Loopmeigen1,sfmeigen1,Looptype1,sftype1,Loopgen1,sfgen1)
        g1u=HSfSfu(argyh,Loopmeigen1,i,Looptype1,Looptypeu,Loopgen1,Loopgenu)
        g2u=HSfSft(argxh,i,Loopmeigen2,Looptypeu,Looptype2,Loopgenu,Loopgen2)
        g3u=VSfSf(Looptype0,sfmeigen2,Loopmeigen2,sftype2,Looptype2,sfgen2,Loopgen2)

!=====================================================!call Loopfunctions==========================================================

!t
        call myCget(m1**2,m3**2,tman,M0Lpt**2,M1Lpt**2,M2Lpt**2)
       
        xc01t = xc0    

        call myCget(tman,m4**2,m2**2,M0Lpt**2,M2Lpt**2,M3Lpt**2)

        xc02t = xc0

        call myCget(m3**2,m4**2,sqrtS**2,M1Lpt**2,M2Lpt**2,M3Lpt**2)

        xc03t = xc0
        xc11t = xc1
        xc21t = xc2

        call myDget(m1**2, m3**2, m4**2, m2**2, tman, sqrtS**2, M0Lpt**2, M1Lpt**2, M2Lpt**2, M3Lpt**2)

        xd0t = xd0

!u
        call myCget(m1**2,m4**2,uman,M0Lpu**2,M1Lpu**2,M2Lpu**2)
       
        xc01u = xc0    

        call myCget(uman,m3**2,m2**2,M0Lpu**2,M2Lpu**2,M3Lpu**2)

        xc02u = xc0

        call myCget(m4**2,m3**2,sqrtS**2,M1Lpu**2,M2Lpu**2,M3Lpu**2)

        xc03u = xc0
        xc11u = xc1
        xc21u = xc2

        call myDget(m1**2, m4**2, m3**2, m2**2, uman, sqrtS**2, M0Lpu**2, M1Lpu**2, M2Lpu**2, M3Lpu**2)

        xd0u = xd0

!==================================================================================================================================

        do j=1,  2   

!Box(t)*t + Box(u)*t
       
        Coupl_Ac = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))        
        Coupl_Bc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
        sfmc = MSf(j,proptypetc,propgentc)
        
!t
	Box_TT = Box_TT + 2*colorfactorTU*dreal(  (Coupl_Ac*Coupl_Bc*g0t*g1t*g2t*g3t*(xc01t + xc02t - xc03t + (-M0Lpt**2 + m1**2 + M1Lpt**2 + m2**2 + M3Lpt**2 - 2*sqrtS**2)*xd0t))/(16.*pi**2*(-sfmc**2 + tman ))  )
!u
	Box_UT = Box_UT + 2*colorfactorTU*dreal(  (Coupl_Ac*Coupl_Bc*g0u*g1u*g2u*g3u*(xc01u + xc02u - xc03u + (-M0Lpu**2 + m1**2 + M1Lpu**2 + m2**2 + M3Lpu**2 - 2*sqrtS**2)*xd0u))/(16.*pi**2*(-sfmc**2 + tman ))  )

!Box(t)*u + Box(u)*u
 
        Coupl_Ac = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
        Coupl_Bc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
        sfmc = MSf(j,proptypeuc,propgenuc)
        
!t
	Box_TU = Box_TU + 2*colorfactorTU*dreal(  (Coupl_Ac*Coupl_Bc*g0t*g1t*g2t*g3t*(xc01t + xc02t - xc03t + (-M0Lpt**2 + m1**2 + M1Lpt**2 + m2**2 + M3Lpt**2 - 2*sqrtS**2)*xd0t))/(16.*pi**2*(-sfmc**2 + uman ))  )
!u
	Box_UU = Box_UU + 2*colorfactorTU*dreal(  (Coupl_Ac*Coupl_Bc*g0u*g1u*g2u*g3u*(xc01u + xc02u - xc03u + (-M0Lpu**2 + m1**2 + M1Lpu**2 + m2**2 + M3Lpu**2 - 2*sqrtS**2)*xd0u))/(16.*pi**2*(-sfmc**2 + uman ))  )


        enddo

!================================================================================================================================== 
!Box(t)*q + Box(u)*q

        Coupl_Ac = dconjg(HHSfSf(argxh,argyh,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))

!t
	Box_TQ = Box_TQ + 2*colorfactorQ*dreal(  -(Coupl_Ac*g0t*g1t*g2t*g3t*(xc01t + xc02t - xc03t + (-M0Lpt**2 + m1**2 + M1Lpt**2 + m2**2 + M3Lpt**2 - 2*sqrtS**2)*xd0t))/(16.*pi**2)   )
!u
	Box_UQ = Box_UQ + 2*colorfactorQ*dreal(  -(Coupl_Ac*g0u*g1u*g2u*g3u*(xc01u + xc02u - xc03u + (-M0Lpu**2 + m1**2 + M1Lpu**2 + m2**2 + M3Lpu**2 - 2*sqrtS**2)*xd0u))/(16.*pi**2)   )

!==================================================================================================================================     
!Box(t)*s(Higgs) + Box(u)*s(Higgs)
       
        do j=sumstartH, sumendH

        Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(ILR(j,argxh,argyh))
        propac = dconjg(DenSmH(j))

!t
	Box_TS = Box_TS + 2*colorfactorSH*dreal(  (Coupl_Ac*Coupl_Bc*g0t*g1t*g2t*g3t*(xc01t + xc02t - xc03t + (-M0Lpt**2 + m1**2 + M1Lpt**2 + m2**2 + M3Lpt**2 - 2*sqrtS**2)*xd0t))/(16.*pi**2*propac)  )
!u
	Box_US = Box_US + 2*colorfactorSH*dreal(  (Coupl_Ac*Coupl_Bc*g0u*g1u*g2u*g3u*(xc01u + xc02u - xc03u + (-M0Lpu**2 + m1**2 + M1Lpu**2 + m2**2 + M3Lpu**2 - 2*sqrtS**2)*xd0u))/(16.*pi**2*propac)  )
  

        enddo

!==================================================================================================================================
!Box(t)*Vector + Box(u)*Vector

        do j=sumstartV,sumendV

        Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = -dconjg(JLR(j,argyh,argxh))
        propac = dconjg(DenSmV(j))
        
!t
	Box_TV = Box_TV + 2*colorfactorSV*dreal(  -(Coupl_Ac*Coupl_Bc*g0t*g1t*g2t*g3t*(tman - uman)*(xc01t + xc02t - xc03t + (-M0Lpt**2 + m1**2 + M1Lpt**2 + m2**2 + M3Lpt**2 - 2*sqrtS**2)*xd0t))/(16.*pi**2*propac)  )
!u
	Box_UV = Box_UV + 2*colorfactorSV*dreal(  -(Coupl_Ac*Coupl_Bc*g0u*g1u*g2u*g3u*(tman - uman)*(xc01u + xc02u - xc03u + (-M0Lpu**2 + m1**2 + M1Lpu**2 + m2**2 + M3Lpu**2 - 2*sqrtS**2)*xd0u))/(16.*pi**2*propac)  )

        enddo           
   
!===================================================================================================================================

      enddo

      HH_TX(1) = Box_TQ
      HH_TX(2) = Box_TS
      HH_TX(3) = Box_TV
      HH_TX(4) = Box_TT
      HH_TX(5) = Box_TU

      HH_TX(6) = HH_TX(1) + HH_TX(2) + HH_TX(3) + HH_TX(4) + HH_TX(5)

      HH_UX(1) = Box_UQ
      HH_UX(2) = Box_US
      HH_UX(3) = Box_UV
      HH_UX(4) = Box_UT
      HH_UX(5) = Box_UU

      HH_UX(6) = HH_UX(1) + HH_UX(2) + HH_UX(3) + HH_UX(4) + HH_UX(5)

!==================================================================================================================================

      end




***********************************************************************************************************************************
***********************************************************************************************************************************





      subroutine QQ2xx_M2TriangleBoxHH(argxh,argyh,HH_TX,HH_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer j,argxh,argyh
      complex*16 C0t,C0u,QQ2xx_BTsf4SS
      double precision   HH_TX(6),HH_UX(6),Box_TQ,Box_TS,Box_TV,Box_TT,
     &                   Box_TU,Box_UQ,Box_US,Box_UV,Box_UT,Box_UU
  
! =================================================================================================================================	
!                             ---                     Higgs-Higgs-final-states                    ---
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

!====================================================!Box-Triangle(t) (one box-propagator replaced by 4-sQuark-Vertex)============================================

      C0t = 1d0/(16d0*pi**2)*QQ2xx_BTsf4SS(0,sqrtS**2,m4**2,m3**2,argyh,argxh,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2) 

!====================================================!Box-Triangle(u) (one box-propagator replaced by 4-sQuark-Vertex)============================================

      C0u = 1d0/(16d0*pi**2)*QQ2xx_BTsf4SS(1,sqrtS**2,m3**2,m4**2,argxh,argyh,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2) 

!=================================================================================================================================================================

      do j=1,  2   

!Box(t)*t + Box(u)*t
       
        Coupl_Ac = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))        
        Coupl_Bc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
        sfmc = MSf(j,proptypetc,propgentc)

!t
	Box_TT = Box_TT + 2*dreal(  (C0t*Coupl_Ac*Coupl_Bc)/(sfmc**2 - tman )   )
!u
	Box_UT = Box_UT + 2*dreal(  (C0u*Coupl_Ac*Coupl_Bc)/(sfmc**2 - tman )   )

!Box(t)*u + Box(u)*u
 
        Coupl_Ac = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
        Coupl_Bc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
        sfmc = MSf(j,proptypeuc,propgenuc)
        
!t
	Box_TU = Box_TU + 2*dreal(  (C0t*Coupl_Ac*Coupl_Bc)/(sfmc**2 - uman )   )
!u
	Box_UU = Box_UU + 2*dreal(  (C0u*Coupl_Ac*Coupl_Bc)/(sfmc**2 - uman )   )


      enddo

!================================================================================================================================== 
!Box(t)*q + Box(u)*q

        Coupl_Ac = dconjg(HHSfSf(argxh,argyh,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))

!t
        Box_TQ = Box_TQ + 2*dreal(  C0t*Coupl_Ac  )
!u
        Box_UQ = Box_UQ + 2*dreal(  C0u*Coupl_Ac  )

!==================================================================================================================================     
!Box(t)*s(Higgs) + Box(u)*s(Higgs)
       
      do j=sumstartH, sumendH

        Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(ILR(j,argxh,argyh))
        propac = dconjg(DenSmH(j))

!t
	Box_TS = Box_TS + 2*dreal( -((C0t*Coupl_Ac*Coupl_Bc)/propac)   )
!u
	Box_US = Box_US + 2*dreal( -((C0u*Coupl_Ac*Coupl_Bc)/propac)   )
  
      enddo

!==================================================================================================================================
!Box(t)*s(Vector) + Box(u)*s(Vector)

      do j=sumstartV, sumendV

        Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = -dconjg(JLR(j,argyh,argxh))
        propac = dconjg(DenSmV(j))
        
!t
	Box_TV = Box_TV + 2*dreal(  (C0t*Coupl_Ac*Coupl_Bc*(tman - uman))/propac  )
!u
	Box_UV = Box_UV + 2*dreal(  (C0u*Coupl_Ac*Coupl_Bc*(tman - uman))/propac  )

      enddo           

!==================================================================================================================================

      HH_TX(1) = Box_TQ
      HH_TX(2) = Box_TS
      HH_TX(3) = Box_TV
      HH_TX(4) = Box_TT
      HH_TX(5) = Box_TU

      HH_TX(6) = HH_TX(1) + HH_TX(2) + HH_TX(3) + HH_TX(4) + HH_TX(5)

      HH_UX(1) = Box_UQ
      HH_UX(2) = Box_US
      HH_UX(3) = Box_UV
      HH_UX(4) = Box_UT
      HH_UX(5) = Box_UU

      HH_UX(6) = HH_UX(1) + HH_UX(2) + HH_UX(3) + HH_UX(4) + HH_UX(5)

!===================================================================================================================================

      end



***********************************************************************************************************************************
***********************************************************************************************************************************



      subroutine QQ2xx_M2GluinoBoxHH(argxh,argyh,HH_TX,HH_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"
#include "../util/DMNLO_LoopIntegrals.h"

      integer j,argxh,argyh
      complex*16 g0t,g1t,g2t,g3t,g0u,g1u,g2u,g3u,g0Lt,g0Rt,g1Lt,g1Rt,g2Lt,g2Rt,g3Lt,g3Rt,g0Lu,g0Ru,g1Lu,g1Ru,g2Lu,g2Ru,g3Lu,g3Ru,
     &           xb01t,xb02t,xc01t,xc02t,xc03t,xc04t,xc11t,xc21t,xd0t,xb01u,xb02u,xc01u,xc02u,xc03u,xc04u,xc11u,xc21u,xd0u,C0t,C0u
      double precision M0Lpt,M1Lpt,M2Lpt,M3Lpt,M0Lpu,M1Lpu,M2Lpu,M3Lpu,HH_TX(6),HH_UX(6),Box_TQ,Box_TS,Box_TV,Box_TT,Box_TU,Box_UQ,Box_US,Box_UV,Box_UT,Box_UU
        
      external myBget,myCget,myDget
  
! =================================================================================================================================	
!                             ---                     Higgs-Higgs-final-states                    ---
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
      g1Lt=DLt(argxh,Looptype1,Looptypet,Loopgen1,Loopgent)
      g1Rt=DRt(argxh,Looptype1,Looptypet,Loopgen1,Loopgent)
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
      g2Lu=DLt(argxh,Looptypeu,Looptype2,Loopgenu,Loopgen2)
      g2Ru=DRt(argxh,Looptypeu,Looptype2,Loopgenu,Loopgen2)
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

      call myCget(tman,m4**2,m2**2,M0Lpt**2,M2Lpt**2,M3Lpt**2)

      xc03t = xc0

      call myCget(m3**2,m4**2,sqrtS**2,M1Lpt**2,M2Lpt**2,M3Lpt**2)

      xc04t = xc0
      xc11t = xc1
      xc21t = xc2

      call myDget(m1**2, m3**2, m4**2, m2**2, tman, sqrtS**2, M0Lpt**2, M1Lpt**2, M2Lpt**2, M3Lpt**2)

      xd0t = xd0

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

      call myCget(m4**2,m3**2,sqrtS**2,M1Lpu**2,M2Lpu**2,M3Lpu**2)

      xc04u = xc0
      xc11u = xc1
      xc21u = xc2

      call myDget(m1**2, m4**2, m3**2, m2**2, uman, sqrtS**2, M0Lpu**2, M1Lpu**2, M2Lpu**2, M3Lpu**2)

      xd0u = xd0

!==================================================================================================================================

      do j=1,  2   

!Box(t)*t + Box(u)*t
       
        Coupl_Ac = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))        
        Coupl_Bc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
        sfmc = MSf(j,proptypetc,propgentc)
        
!t
	Box_TT = Box_TT + 2*colorfactorTU*dreal(  -(Coupl_Ac*Coupl_Bc*(g0Rt*(2*g1Rt*M1Lpt*(g2Rt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*xd0t + g3Lt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t)) + 
     &             g2Lt*(g3Rt*M0Lpt*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + g3Lt*M3Lpt*(xc01t + xc04t + (M0Lpt**2 + M2Lpt**2 - tman)*xd0t))) + 
     &          g1Lt*(2*g2Lt*M2Lpt*(g3Lt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Rt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &             g2Rt*(2*g3Rt*M0Lpt*M3Lpt*(xc01t + xc03t + (M1Lpt**2 + M2Lpt**2 - m3**2)*xd0t) + 
     &                g3Lt*(2*xb01t + 2*xb02t - m2**2*xc01t + sqrtS**2*xc01t + 2*M1Lpt**2*xc02t - m3**2*xc02t + tman*xc02t - m2**2*xc03t + 2*M3Lpt**2*xc03t - m4**2*xc03t + tman*xc03t + 
     &                   2*M1Lpt**2*xc04t + m3**2*xc04t - m4**2*xc04t + sqrtS**2*xc04t + 4*m3**2*xc11t + 2*m3**2*xc21t - 2*m4**2*xc21t + 2*sqrtS**2*xc21t + 
     &                   (m2**2*(-M2Lpt**2 + m3**2) - m3**2*M3Lpt**2 + M2Lpt**2*sqrtS**2 - M1Lpt**2*(m2**2 - 2*M3Lpt**2 + m4**2 - tman) + M3Lpt**2*tman - sqrtS**2*tman)*xd0t - 
     &                   m1**2*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + M0Lpt**2*(2*xc01t + (2*M2Lpt**2 - m3**2 - m4**2 + sqrtS**2)*xd0t))))) + 
     &       g0Lt*(2*g1Lt*M1Lpt*(g2Lt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*xd0t + g3Rt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t)) + 
     &             g2Rt*(g3Lt*M0Lpt*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + g3Rt*M3Lpt*(xc01t + xc04t + (M0Lpt**2 + M2Lpt**2 - tman)*xd0t))) + 
     &          g1Rt*(2*g2Rt*M2Lpt*(g3Rt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Lt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &             g2Lt*(2*g3Lt*M0Lpt*M3Lpt*(xc01t + xc03t + (M1Lpt**2 + M2Lpt**2 - m3**2)*xd0t) + 
     &                g3Rt*(2*xb01t + 2*xb02t - m2**2*xc01t + sqrtS**2*xc01t + 2*M1Lpt**2*xc02t - m3**2*xc02t + tman*xc02t - m2**2*xc03t + 2*M3Lpt**2*xc03t - m4**2*xc03t + tman*xc03t + 
     &                   2*M1Lpt**2*xc04t + m3**2*xc04t - m4**2*xc04t + sqrtS**2*xc04t + 4*m3**2*xc11t + 2*m3**2*xc21t - 2*m4**2*xc21t + 2*sqrtS**2*xc21t + 
     &                   (m2**2*(-M2Lpt**2 + m3**2) - m3**2*M3Lpt**2 + M2Lpt**2*sqrtS**2 - M1Lpt**2*(m2**2 - 2*M3Lpt**2 + m4**2 - tman) + M3Lpt**2*tman - sqrtS**2*tman)*xd0t - 
     &                   m1**2*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + M0Lpt**2*(2*xc01t + (2*M2Lpt**2 - m3**2 - m4**2 + sqrtS**2)*xd0t)))))))/(32.*pi**2*(sfmc**2 - tman ))  )
!u
	Box_UT = Box_UT + 2*colorfactorTU*dreal(   -(Coupl_Ac*Coupl_Bc*(g0Ru*(2*g1Ru*M1Lpu*(g2Ru*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*xd0u + g3Lu*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u)) + 
     &             g2Lu*(g3Ru*M0Lpu*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + g3Lu*M3Lpu*(xc01u + xc04u + (M0Lpu**2 + M2Lpu**2 - uman)*xd0u))) + 
     &          g1Lu*(2*g2Lu*M2Lpu*(g3Lu*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Ru*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &             g2Ru*(2*g3Ru*M0Lpu*M3Lpu*(xc01u + xc03u + (M1Lpu**2 + M2Lpu**2 - m4**2)*xd0u) + 
     &                g3Lu*(2*xb01u + 2*xb02u - m2**2*xc01u + sqrtS**2*xc01u + (2*M1Lpu**2 - m4**2 + uman)*xc02u - m2**2*xc03u - m3**2*xc03u + 2*M3Lpu**2*xc03u + uman*xc03u + 
     &                   2*M1Lpu**2*xc04u - m3**2*xc04u + m4**2*xc04u + sqrtS**2*xc04u + 4*m4**2*xc11u - 2*m3**2*xc21u + 2*m4**2*xc21u + 2*sqrtS**2*xc21u + 
     &                   (-(M3Lpu**2*m4**2) + m2**2*(-M2Lpu**2 + m4**2) + M2Lpu**2*sqrtS**2 - M1Lpu**2*(m2**2 + m3**2 - 2*M3Lpu**2 - uman) + M3Lpu**2*uman - sqrtS**2*uman)*xd0u - 
     &                   m1**2*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + M0Lpu**2*(2*xc01u + (2*M2Lpu**2 - m3**2 - m4**2 + sqrtS**2)*xd0u))))) + 
     &       g0Lu*(2*g1Lu*M1Lpu*(g2Lu*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*xd0u + g3Ru*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u)) + 
     &             g2Ru*(g3Lu*M0Lpu*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + g3Ru*M3Lpu*(xc01u + xc04u + (M0Lpu**2 + M2Lpu**2 - uman)*xd0u))) + 
     &          g1Ru*(2*g2Ru*M2Lpu*(g3Ru*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Lu*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &             g2Lu*(2*g3Lu*M0Lpu*M3Lpu*(xc01u + xc03u + (M1Lpu**2 + M2Lpu**2 - m4**2)*xd0u) + 
     &                g3Ru*(2*xb01u + 2*xb02u - m2**2*xc01u + sqrtS**2*xc01u + (2*M1Lpu**2 - m4**2 + uman)*xc02u - m2**2*xc03u - m3**2*xc03u + 2*M3Lpu**2*xc03u + uman*xc03u + 
     &                   2*M1Lpu**2*xc04u - m3**2*xc04u + m4**2*xc04u + sqrtS**2*xc04u + 4*m4**2*xc11u - 2*m3**2*xc21u + 2*m4**2*xc21u + 2*sqrtS**2*xc21u + 
     &                   (-(M3Lpu**2*m4**2) + m2**2*(-M2Lpu**2 + m4**2) + M2Lpu**2*sqrtS**2 - M1Lpu**2*(m2**2 + m3**2 - 2*M3Lpu**2 - uman) + M3Lpu**2*uman - sqrtS**2*uman)*xd0u - 
     &                   m1**2*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + M0Lpu**2*(2*xc01u + (2*M2Lpu**2 - m3**2 - m4**2 + sqrtS**2)*xd0u)))))))/(32.*pi**2*(sfmc**2 - tman ))  )

!Box(t)*u + Box(u)*u
 
        Coupl_Ac = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
        Coupl_Bc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
        sfmc = MSf(j,proptypeuc,propgenuc)
        
!t
	Box_TU = Box_TU + 2*colorfactorTU*dreal(  -(Coupl_Ac*Coupl_Bc*(g0Rt*(2*g1Rt*M1Lpt*(g2Rt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*xd0t + g3Lt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t)) + 
     &             g2Lt*(g3Rt*M0Lpt*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + g3Lt*M3Lpt*(xc01t + xc04t + (M0Lpt**2 + M2Lpt**2 - tman)*xd0t))) + 
     &          g1Lt*(2*g2Lt*M2Lpt*(g3Lt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Rt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &             g2Rt*(2*g3Rt*M0Lpt*M3Lpt*(xc01t + xc03t + (M1Lpt**2 + M2Lpt**2 - m3**2)*xd0t) + 
     &                g3Lt*(2*xb01t + 2*xb02t - m2**2*xc01t + sqrtS**2*xc01t + 2*M1Lpt**2*xc02t - m3**2*xc02t + tman*xc02t - m2**2*xc03t + 2*M3Lpt**2*xc03t - m4**2*xc03t + tman*xc03t + 
     &                   2*M1Lpt**2*xc04t + m3**2*xc04t - m4**2*xc04t + sqrtS**2*xc04t + 4*m3**2*xc11t + 2*m3**2*xc21t - 2*m4**2*xc21t + 2*sqrtS**2*xc21t + 
     &                   (m2**2*(-M2Lpt**2 + m3**2) - m3**2*M3Lpt**2 + M2Lpt**2*sqrtS**2 - M1Lpt**2*(m2**2 - 2*M3Lpt**2 + m4**2 - tman) + M3Lpt**2*tman - sqrtS**2*tman)*xd0t - 
     &                   m1**2*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + M0Lpt**2*(2*xc01t + (2*M2Lpt**2 - m3**2 - m4**2 + sqrtS**2)*xd0t))))) + 
     &       g0Lt*(2*g1Lt*M1Lpt*(g2Lt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*xd0t + g3Rt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t)) + 
     &             g2Rt*(g3Lt*M0Lpt*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + g3Rt*M3Lpt*(xc01t + xc04t + (M0Lpt**2 + M2Lpt**2 - tman)*xd0t))) + 
     &          g1Rt*(2*g2Rt*M2Lpt*(g3Rt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Lt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &             g2Lt*(2*g3Lt*M0Lpt*M3Lpt*(xc01t + xc03t + (M1Lpt**2 + M2Lpt**2 - m3**2)*xd0t) + 
     &                g3Rt*(2*xb01t + 2*xb02t - m2**2*xc01t + sqrtS**2*xc01t + 2*M1Lpt**2*xc02t - m3**2*xc02t + tman*xc02t - m2**2*xc03t + 2*M3Lpt**2*xc03t - m4**2*xc03t + tman*xc03t + 
     &                   2*M1Lpt**2*xc04t + m3**2*xc04t - m4**2*xc04t + sqrtS**2*xc04t + 4*m3**2*xc11t + 2*m3**2*xc21t - 2*m4**2*xc21t + 2*sqrtS**2*xc21t + 
     &                   (m2**2*(-M2Lpt**2 + m3**2) - m3**2*M3Lpt**2 + M2Lpt**2*sqrtS**2 - M1Lpt**2*(m2**2 - 2*M3Lpt**2 + m4**2 - tman) + M3Lpt**2*tman - sqrtS**2*tman)*xd0t - 
     &                   m1**2*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + M0Lpt**2*(2*xc01t + (2*M2Lpt**2 - m3**2 - m4**2 + sqrtS**2)*xd0t)))))))/(32.*pi**2*(sfmc**2 - uman ))  )
!u
	Box_UU = Box_UU + 2*colorfactorTU*dreal(  -(Coupl_Ac*Coupl_Bc*(g0Ru*(2*g1Ru*M1Lpu*(g2Ru*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*xd0u + g3Lu*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u)) + 
     &             g2Lu*(g3Ru*M0Lpu*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + g3Lu*M3Lpu*(xc01u + xc04u + (M0Lpu**2 + M2Lpu**2 - uman)*xd0u))) + 
     &          g1Lu*(2*g2Lu*M2Lpu*(g3Lu*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Ru*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &             g2Ru*(2*g3Ru*M0Lpu*M3Lpu*(xc01u + xc03u + (M1Lpu**2 + M2Lpu**2 - m4**2)*xd0u) + 
     &                g3Lu*(2*xb01u + 2*xb02u - m2**2*xc01u + sqrtS**2*xc01u + (2*M1Lpu**2 - m4**2 + uman)*xc02u - m2**2*xc03u - m3**2*xc03u + 2*M3Lpu**2*xc03u + uman*xc03u + 
     &                   2*M1Lpu**2*xc04u - m3**2*xc04u + m4**2*xc04u + sqrtS**2*xc04u + 4*m4**2*xc11u - 2*m3**2*xc21u + 2*m4**2*xc21u + 2*sqrtS**2*xc21u + 
     &                   (-(M3Lpu**2*m4**2) + m2**2*(-M2Lpu**2 + m4**2) + M2Lpu**2*sqrtS**2 - M1Lpu**2*(m2**2 + m3**2 - 2*M3Lpu**2 - uman) + M3Lpu**2*uman - sqrtS**2*uman)*xd0u - 
     &                   m1**2*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + M0Lpu**2*(2*xc01u + (2*M2Lpu**2 - m3**2 - m4**2 + sqrtS**2)*xd0u))))) + 
     &       g0Lu*(2*g1Lu*M1Lpu*(g2Lu*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*xd0u + g3Ru*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u)) + 
     &             g2Ru*(g3Lu*M0Lpu*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + g3Ru*M3Lpu*(xc01u + xc04u + (M0Lpu**2 + M2Lpu**2 - uman)*xd0u))) + 
     &          g1Ru*(2*g2Ru*M2Lpu*(g3Ru*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Lu*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &             g2Lu*(2*g3Lu*M0Lpu*M3Lpu*(xc01u + xc03u + (M1Lpu**2 + M2Lpu**2 - m4**2)*xd0u) + 
     &                g3Ru*(2*xb01u + 2*xb02u - m2**2*xc01u + sqrtS**2*xc01u + (2*M1Lpu**2 - m4**2 + uman)*xc02u - m2**2*xc03u - m3**2*xc03u + 2*M3Lpu**2*xc03u + uman*xc03u + 
     &                   2*M1Lpu**2*xc04u - m3**2*xc04u + m4**2*xc04u + sqrtS**2*xc04u + 4*m4**2*xc11u - 2*m3**2*xc21u + 2*m4**2*xc21u + 2*sqrtS**2*xc21u + 
     &                   (-(M3Lpu**2*m4**2) + m2**2*(-M2Lpu**2 + m4**2) + M2Lpu**2*sqrtS**2 - M1Lpu**2*(m2**2 + m3**2 - 2*M3Lpu**2 - uman) + M3Lpu**2*uman - sqrtS**2*uman)*xd0u - 
     &                   m1**2*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + M0Lpu**2*(2*xc01u + (2*M2Lpu**2 - m3**2 - m4**2 + sqrtS**2)*xd0u)))))))/(32.*pi**2*(sfmc**2 - uman ))  )


        enddo

!================================================================================================================================== 
!Box(t)*q + Box(u)*q

        Coupl_Ac = dconjg(HHSfSf(argxh,argyh,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))

!t
	Box_TQ = Box_TQ + 2*colorfactorQ*dreal(  -(Coupl_Ac*(g0Rt*(2*g1Rt*M1Lpt*(g2Rt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*xd0t + g3Lt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t)) + 
     &             g2Lt*(g3Rt*M0Lpt*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + g3Lt*M3Lpt*(xc01t + xc04t + (M0Lpt**2 + M2Lpt**2 - tman)*xd0t))) + 
     &          g1Lt*(2*g2Lt*M2Lpt*(g3Lt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Rt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &             g2Rt*(2*g3Rt*M0Lpt*M3Lpt*(xc01t + xc03t + (M1Lpt**2 + M2Lpt**2 - m3**2)*xd0t) + 
     &                g3Lt*(2*xb01t + 2*xb02t - m2**2*xc01t + sqrtS**2*xc01t + 2*M1Lpt**2*xc02t - m3**2*xc02t + tman*xc02t - m2**2*xc03t + 2*M3Lpt**2*xc03t - m4**2*xc03t + tman*xc03t + 
     &                   2*M1Lpt**2*xc04t + m3**2*xc04t - m4**2*xc04t + sqrtS**2*xc04t + 4*m3**2*xc11t + 2*m3**2*xc21t - 2*m4**2*xc21t + 2*sqrtS**2*xc21t + 
     &                   (m2**2*(-M2Lpt**2 + m3**2) - m3**2*M3Lpt**2 + M2Lpt**2*sqrtS**2 - M1Lpt**2*(m2**2 - 2*M3Lpt**2 + m4**2 - tman) + M3Lpt**2*tman - sqrtS**2*tman)*xd0t - 
     &                   m1**2*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + M0Lpt**2*(2*xc01t + (2*M2Lpt**2 - m3**2 - m4**2 + sqrtS**2)*xd0t))))) + 
     &       g0Lt*(2*g1Lt*M1Lpt*(g2Lt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*xd0t + g3Rt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t)) + 
     &             g2Rt*(g3Lt*M0Lpt*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + g3Rt*M3Lpt*(xc01t + xc04t + (M0Lpt**2 + M2Lpt**2 - tman)*xd0t))) + 
     &          g1Rt*(2*g2Rt*M2Lpt*(g3Rt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Lt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &             g2Lt*(2*g3Lt*M0Lpt*M3Lpt*(xc01t + xc03t + (M1Lpt**2 + M2Lpt**2 - m3**2)*xd0t) + 
     &                g3Rt*(2*xb01t + 2*xb02t - m2**2*xc01t + sqrtS**2*xc01t + 2*M1Lpt**2*xc02t - m3**2*xc02t + tman*xc02t - m2**2*xc03t + 2*M3Lpt**2*xc03t - m4**2*xc03t + tman*xc03t + 
     &                   2*M1Lpt**2*xc04t + m3**2*xc04t - m4**2*xc04t + sqrtS**2*xc04t + 4*m3**2*xc11t + 2*m3**2*xc21t - 2*m4**2*xc21t + 2*sqrtS**2*xc21t + 
     &                   (m2**2*(-M2Lpt**2 + m3**2) - m3**2*M3Lpt**2 + M2Lpt**2*sqrtS**2 - M1Lpt**2*(m2**2 - 2*M3Lpt**2 + m4**2 - tman) + M3Lpt**2*tman - sqrtS**2*tman)*xd0t - 
     &                   m1**2*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + M0Lpt**2*(2*xc01t + (2*M2Lpt**2 - m3**2 - m4**2 + sqrtS**2)*xd0t)))))))/(32.*pi**2)  )
!u
	Box_UQ = Box_UQ + 2*colorfactorQ*dreal(  -(Coupl_Ac*(g0Ru*(2*g1Ru*M1Lpu*(g2Ru*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*xd0u + g3Lu*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u)) + 
     &             g2Lu*(g3Ru*M0Lpu*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + g3Lu*M3Lpu*(xc01u + xc04u + (M0Lpu**2 + M2Lpu**2 - uman)*xd0u))) + 
     &          g1Lu*(2*g2Lu*M2Lpu*(g3Lu*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Ru*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &             g2Ru*(2*g3Ru*M0Lpu*M3Lpu*(xc01u + xc03u + (M1Lpu**2 + M2Lpu**2 - m4**2)*xd0u) + 
     &                g3Lu*(2*xb01u + 2*xb02u - m2**2*xc01u + sqrtS**2*xc01u + (2*M1Lpu**2 - m4**2 + uman)*xc02u - m2**2*xc03u - m3**2*xc03u + 2*M3Lpu**2*xc03u + uman*xc03u + 
     &                   2*M1Lpu**2*xc04u - m3**2*xc04u + m4**2*xc04u + sqrtS**2*xc04u + 4*m4**2*xc11u - 2*m3**2*xc21u + 2*m4**2*xc21u + 2*sqrtS**2*xc21u + 
     &                   (-(M3Lpu**2*m4**2) + m2**2*(-M2Lpu**2 + m4**2) + M2Lpu**2*sqrtS**2 - M1Lpu**2*(m2**2 + m3**2 - 2*M3Lpu**2 - uman) + M3Lpu**2*uman - sqrtS**2*uman)*xd0u - 
     &                   m1**2*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + M0Lpu**2*(2*xc01u + (2*M2Lpu**2 - m3**2 - m4**2 + sqrtS**2)*xd0u))))) + 
     &       g0Lu*(2*g1Lu*M1Lpu*(g2Lu*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*xd0u + g3Ru*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u)) + 
     &             g2Ru*(g3Lu*M0Lpu*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + g3Ru*M3Lpu*(xc01u + xc04u + (M0Lpu**2 + M2Lpu**2 - uman)*xd0u))) + 
     &          g1Ru*(2*g2Ru*M2Lpu*(g3Ru*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Lu*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &             g2Lu*(2*g3Lu*M0Lpu*M3Lpu*(xc01u + xc03u + (M1Lpu**2 + M2Lpu**2 - m4**2)*xd0u) + 
     &                g3Ru*(2*xb01u + 2*xb02u - m2**2*xc01u + sqrtS**2*xc01u + (2*M1Lpu**2 - m4**2 + uman)*xc02u - m2**2*xc03u - m3**2*xc03u + 2*M3Lpu**2*xc03u + uman*xc03u + 
     &                   2*M1Lpu**2*xc04u - m3**2*xc04u + m4**2*xc04u + sqrtS**2*xc04u + 4*m4**2*xc11u - 2*m3**2*xc21u + 2*m4**2*xc21u + 2*sqrtS**2*xc21u + 
     &                   (-(M3Lpu**2*m4**2) + m2**2*(-M2Lpu**2 + m4**2) + M2Lpu**2*sqrtS**2 - M1Lpu**2*(m2**2 + m3**2 - 2*M3Lpu**2 - uman) + M3Lpu**2*uman - sqrtS**2*uman)*xd0u - 
     &                   m1**2*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + M0Lpu**2*(2*xc01u + (2*M2Lpu**2 - m3**2 - m4**2 + sqrtS**2)*xd0u)))))))/(32.*pi**2)  )

!==================================================================================================================================   
!Box(t)*s(Higgs) + Box(u)*s(Higgs)
       
        do j= sumstartH, sumendH

        Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(ILR(j,argxh,argyh))
        propac = dconjg(DenSmH(j))

!t
	Box_TS = Box_TS + 2*colorfactorSH*dreal(  (Coupl_Ac*Coupl_Bc*(g0Rt*(2*g1Rt*M1Lpt*(g2Rt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*xd0t + g3Lt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t)) + 
     &            g2Lt*(g3Rt*M0Lpt*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + g3Lt*M3Lpt*(xc01t + xc04t + (M0Lpt**2 + M2Lpt**2 - tman)*xd0t))) + 
     &         g1Lt*(2*g2Lt*M2Lpt*(g3Lt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Rt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &            g2Rt*(2*g3Rt*M0Lpt*M3Lpt*(xc01t + xc03t + (M1Lpt**2 + M2Lpt**2 - m3**2)*xd0t) + 
     &               g3Lt*(2*xb01t + 2*xb02t - m2**2*xc01t + sqrtS**2*xc01t + 2*M1Lpt**2*xc02t - m3**2*xc02t + tman*xc02t - m2**2*xc03t + 2*M3Lpt**2*xc03t - m4**2*xc03t + tman*xc03t + 
     &                  2*M1Lpt**2*xc04t + m3**2*xc04t - m4**2*xc04t + sqrtS**2*xc04t + 4*m3**2*xc11t + 2*m3**2*xc21t - 2*m4**2*xc21t + 2*sqrtS**2*xc21t + 
     &                  (m2**2*(-M2Lpt**2 + m3**2) - m3**2*M3Lpt**2 + M2Lpt**2*sqrtS**2 - M1Lpt**2*(m2**2 - 2*M3Lpt**2 + m4**2 - tman) + M3Lpt**2*tman - sqrtS**2*tman)*xd0t - 
     &                  m1**2*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + M0Lpt**2*(2*xc01t + (2*M2Lpt**2 - m3**2 - m4**2 + sqrtS**2)*xd0t))))) + 
     &      g0Lt*(2*g1Lt*M1Lpt*(g2Lt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*xd0t + g3Rt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t)) + 
     &            g2Rt*(g3Lt*M0Lpt*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + g3Rt*M3Lpt*(xc01t + xc04t + (M0Lpt**2 + M2Lpt**2 - tman)*xd0t))) + 
     &         g1Rt*(2*g2Rt*M2Lpt*(g3Rt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Lt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &            g2Lt*(2*g3Lt*M0Lpt*M3Lpt*(xc01t + xc03t + (M1Lpt**2 + M2Lpt**2 - m3**2)*xd0t) + 
     &               g3Rt*(2*xb01t + 2*xb02t - m2**2*xc01t + sqrtS**2*xc01t + 2*M1Lpt**2*xc02t - m3**2*xc02t + tman*xc02t - m2**2*xc03t + 2*M3Lpt**2*xc03t - m4**2*xc03t + tman*xc03t + 
     &                  2*M1Lpt**2*xc04t + m3**2*xc04t - m4**2*xc04t + sqrtS**2*xc04t + 4*m3**2*xc11t + 2*m3**2*xc21t - 2*m4**2*xc21t + 2*sqrtS**2*xc21t + 
     &                  (m2**2*(-M2Lpt**2 + m3**2) - m3**2*M3Lpt**2 + M2Lpt**2*sqrtS**2 - M1Lpt**2*(m2**2 - 2*M3Lpt**2 + m4**2 - tman) + M3Lpt**2*tman - sqrtS**2*tman)*xd0t - 
     &                  m1**2*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + M0Lpt**2*(2*xc01t + (2*M2Lpt**2 - m3**2 - m4**2 + sqrtS**2)*xd0t)))))))/(32.*pi**2*propac)  )
!u
	Box_US = Box_US + 2*colorfactorSH*dreal(     (Coupl_Ac*Coupl_Bc*(g0Ru*(2*g1Ru*M1Lpu*(g2Ru*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*xd0u + g3Lu*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u)) + 
     &            g2Lu*(g3Ru*M0Lpu*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + g3Lu*M3Lpu*(xc01u + xc04u + (M0Lpu**2 + M2Lpu**2 - uman)*xd0u))) + 
     &         g1Lu*(2*g2Lu*M2Lpu*(g3Lu*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Ru*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &            g2Ru*(2*g3Ru*M0Lpu*M3Lpu*(xc01u + xc03u + (M1Lpu**2 + M2Lpu**2 - m4**2)*xd0u) + 
     &               g3Lu*(2*xb01u + 2*xb02u - m2**2*xc01u + sqrtS**2*xc01u + (2*M1Lpu**2 - m4**2 + uman)*xc02u - m2**2*xc03u - m3**2*xc03u + 2*M3Lpu**2*xc03u + uman*xc03u + 
     &                  2*M1Lpu**2*xc04u - m3**2*xc04u + m4**2*xc04u + sqrtS**2*xc04u + 4*m4**2*xc11u - 2*m3**2*xc21u + 2*m4**2*xc21u + 2*sqrtS**2*xc21u + 
     &                  (-(M3Lpu**2*m4**2) + m2**2*(-M2Lpu**2 + m4**2) + M2Lpu**2*sqrtS**2 - M1Lpu**2*(m2**2 + m3**2 - 2*M3Lpu**2 - uman) + M3Lpu**2*uman - sqrtS**2*uman)*xd0u - 
     &                  m1**2*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + M0Lpu**2*(2*xc01u + (2*M2Lpu**2 - m3**2 - m4**2 + sqrtS**2)*xd0u))))) + 
     &      g0Lu*(2*g1Lu*M1Lpu*(g2Lu*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*xd0u + g3Ru*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u)) + 
     &            g2Ru*(g3Lu*M0Lpu*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + g3Ru*M3Lpu*(xc01u + xc04u + (M0Lpu**2 + M2Lpu**2 - uman)*xd0u))) + 
     &         g1Ru*(2*g2Ru*M2Lpu*(g3Ru*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Lu*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &            g2Lu*(2*g3Lu*M0Lpu*M3Lpu*(xc01u + xc03u + (M1Lpu**2 + M2Lpu**2 - m4**2)*xd0u) + 
     &               g3Ru*(2*xb01u + 2*xb02u - m2**2*xc01u + sqrtS**2*xc01u + (2*M1Lpu**2 - m4**2 + uman)*xc02u - m2**2*xc03u - m3**2*xc03u + 2*M3Lpu**2*xc03u + uman*xc03u + 
     &                  2*M1Lpu**2*xc04u - m3**2*xc04u + m4**2*xc04u + sqrtS**2*xc04u + 4*m4**2*xc11u - 2*m3**2*xc21u + 2*m4**2*xc21u + 2*sqrtS**2*xc21u + 
     &                  (-(M3Lpu**2*m4**2) + m2**2*(-M2Lpu**2 + m4**2) + M2Lpu**2*sqrtS**2 - M1Lpu**2*(m2**2 + m3**2 - 2*M3Lpu**2 - uman) + M3Lpu**2*uman - sqrtS**2*uman)*xd0u - 
     &                  m1**2*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + M0Lpu**2*(2*xc01u + (2*M2Lpu**2 - m3**2 - m4**2 + sqrtS**2)*xd0u)))))))/(32.*pi**2*propac)  )
  

        enddo

!==================================================================================================================================
!Box(t)*s(Vector) + Box(u)*s(Vector)

        do j=sumstartV, sumendV

        Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = -dconjg(JLR(j,argyh,argxh))
        propac = dconjg(DenSmV(j))
        
!t
	Box_TV = Box_TV + 2*colorfactorSV*dreal(   -(Coupl_Ac*Coupl_Bc*(tman - uman)*(g0Rt*(2*g1Rt*M1Lpt*(g2Rt*M2Lpt*(2*g3Rt*M0Lpt*M3Lpt*xd0t + g3Lt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t)) + 
     &             g2Lt*(g3Rt*M0Lpt*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + g3Lt*M3Lpt*(xc01t + xc04t + (M0Lpt**2 + M2Lpt**2 - tman)*xd0t))) + 
     &          g1Lt*(2*g2Lt*M2Lpt*(g3Lt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Rt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &             g2Rt*(2*g3Rt*M0Lpt*M3Lpt*(xc01t + xc03t + (M1Lpt**2 + M2Lpt**2 - m3**2)*xd0t) + 
     &                g3Lt*(2*xb01t + 2*xb02t - m2**2*xc01t + sqrtS**2*xc01t + 2*M1Lpt**2*xc02t - m3**2*xc02t + tman*xc02t - m2**2*xc03t + 2*M3Lpt**2*xc03t - m4**2*xc03t + tman*xc03t + 
     &                   2*M1Lpt**2*xc04t + m3**2*xc04t - m4**2*xc04t + sqrtS**2*xc04t + 4*m3**2*xc11t + 2*m3**2*xc21t - 2*m4**2*xc21t + 2*sqrtS**2*xc21t + 
     &                   (m2**2*(-M2Lpt**2 + m3**2) - m3**2*M3Lpt**2 + M2Lpt**2*sqrtS**2 - M1Lpt**2*(m2**2 - 2*M3Lpt**2 + m4**2 - tman) + M3Lpt**2*tman - sqrtS**2*tman)*xd0t - 
     &                   m1**2*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + M0Lpt**2*(2*xc01t + (2*M2Lpt**2 - m3**2 - m4**2 + sqrtS**2)*xd0t))))) + 
     &       g0Lt*(2*g1Lt*M1Lpt*(g2Lt*M2Lpt*(2*g3Lt*M0Lpt*M3Lpt*xd0t + g3Rt*(xc02t + xc04t + (M0Lpt**2 - m2**2 + M3Lpt**2)*xd0t)) + 
     &             g2Rt*(g3Lt*M0Lpt*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + g3Rt*M3Lpt*(xc01t + xc04t + (M0Lpt**2 + M2Lpt**2 - tman)*xd0t))) + 
     &          g1Rt*(2*g2Rt*M2Lpt*(g3Rt*M3Lpt*(xc03t + xc04t + (M0Lpt**2 - m1**2 + M1Lpt**2)*xd0t) + g3Lt*M0Lpt*(xc02t + xc03t + (M1Lpt**2 + M3Lpt**2 - sqrtS**2)*xd0t)) + 
     &             g2Lt*(2*g3Lt*M0Lpt*M3Lpt*(xc01t + xc03t + (M1Lpt**2 + M2Lpt**2 - m3**2)*xd0t) + 
     &                g3Rt*(2*xb01t + 2*xb02t - m2**2*xc01t + sqrtS**2*xc01t + 2*M1Lpt**2*xc02t - m3**2*xc02t + tman*xc02t - m2**2*xc03t + 2*M3Lpt**2*xc03t - m4**2*xc03t + tman*xc03t + 
     &                   2*M1Lpt**2*xc04t + m3**2*xc04t - m4**2*xc04t + sqrtS**2*xc04t + 4*m3**2*xc11t + 2*m3**2*xc21t - 2*m4**2*xc21t + 2*sqrtS**2*xc21t + 
     &                   (m2**2*(-M2Lpt**2 + m3**2) - m3**2*M3Lpt**2 + M2Lpt**2*sqrtS**2 - M1Lpt**2*(m2**2 - 2*M3Lpt**2 + m4**2 - tman) + M3Lpt**2*tman - sqrtS**2*tman)*xd0t - 
     &                   m1**2*(xc01t + xc02t + (M2Lpt**2 + M3Lpt**2 - m4**2)*xd0t) + M0Lpt**2*(2*xc01t + (2*M2Lpt**2 - m3**2 - m4**2 + sqrtS**2)*xd0t)))))))/(32.*pi**2*propac)   )
!u
	Box_UV = Box_UV + 2*colorfactorSV*dreal(  -(Coupl_Ac*Coupl_Bc*(tman - uman)*(g0Ru*(2*g1Ru*M1Lpu*(g2Ru*M2Lpu*(2*g3Ru*M0Lpu*M3Lpu*xd0u + g3Lu*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u)) + 
     &             g2Lu*(g3Ru*M0Lpu*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + g3Lu*M3Lpu*(xc01u + xc04u + (M0Lpu**2 + M2Lpu**2 - uman)*xd0u))) + 
     &          g1Lu*(2*g2Lu*M2Lpu*(g3Lu*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Ru*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &             g2Ru*(2*g3Ru*M0Lpu*M3Lpu*(xc01u + xc03u + (M1Lpu**2 + M2Lpu**2 - m4**2)*xd0u) + 
     &                g3Lu*(2*xb01u + 2*xb02u - m2**2*xc01u + sqrtS**2*xc01u + (2*M1Lpu**2 - m4**2 + uman)*xc02u - m2**2*xc03u - m3**2*xc03u + 2*M3Lpu**2*xc03u + uman*xc03u + 
     &                   2*M1Lpu**2*xc04u - m3**2*xc04u + m4**2*xc04u + sqrtS**2*xc04u + 4*m4**2*xc11u - 2*m3**2*xc21u + 2*m4**2*xc21u + 2*sqrtS**2*xc21u + 
     &                   (-(M3Lpu**2*m4**2) + m2**2*(-M2Lpu**2 + m4**2) + M2Lpu**2*sqrtS**2 - M1Lpu**2*(m2**2 + m3**2 - 2*M3Lpu**2 - uman) + M3Lpu**2*uman - sqrtS**2*uman)*xd0u - 
     &                   m1**2*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + M0Lpu**2*(2*xc01u + (2*M2Lpu**2 - m3**2 - m4**2 + sqrtS**2)*xd0u))))) + 
     &       g0Lu*(2*g1Lu*M1Lpu*(g2Lu*M2Lpu*(2*g3Lu*M0Lpu*M3Lpu*xd0u + g3Ru*(xc02u + xc04u + (M0Lpu**2 - m2**2 + M3Lpu**2)*xd0u)) + 
     &             g2Ru*(g3Lu*M0Lpu*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + g3Ru*M3Lpu*(xc01u + xc04u + (M0Lpu**2 + M2Lpu**2 - uman)*xd0u))) + 
     &          g1Ru*(2*g2Ru*M2Lpu*(g3Ru*M3Lpu*(xc03u + xc04u + (M0Lpu**2 - m1**2 + M1Lpu**2)*xd0u) + g3Lu*M0Lpu*(xc02u + xc03u + (M1Lpu**2 + M3Lpu**2 - sqrtS**2)*xd0u)) + 
     &             g2Lu*(2*g3Lu*M0Lpu*M3Lpu*(xc01u + xc03u + (M1Lpu**2 + M2Lpu**2 - m4**2)*xd0u) + 
     &                g3Ru*(2*xb01u + 2*xb02u - m2**2*xc01u + sqrtS**2*xc01u + (2*M1Lpu**2 - m4**2 + uman)*xc02u - m2**2*xc03u - m3**2*xc03u + 2*M3Lpu**2*xc03u + uman*xc03u + 
     &                   2*M1Lpu**2*xc04u - m3**2*xc04u + m4**2*xc04u + sqrtS**2*xc04u + 4*m4**2*xc11u - 2*m3**2*xc21u + 2*m4**2*xc21u + 2*sqrtS**2*xc21u + 
     &                   (-(M3Lpu**2*m4**2) + m2**2*(-M2Lpu**2 + m4**2) + M2Lpu**2*sqrtS**2 - M1Lpu**2*(m2**2 + m3**2 - 2*M3Lpu**2 - uman) + M3Lpu**2*uman - sqrtS**2*uman)*xd0u - 
     &                   m1**2*(xc01u + xc02u + (M2Lpu**2 - m3**2 + M3Lpu**2)*xd0u) + M0Lpu**2*(2*xc01u + (2*M2Lpu**2 - m3**2 - m4**2 + sqrtS**2)*xd0u)))))))/(32.*pi**2*propac)   )

      enddo           

!===================================================================================================================================

      HH_TX(1) = Box_TQ
      HH_TX(2) = Box_TS
      HH_TX(3) = Box_TV
      HH_TX(4) = Box_TT
      HH_TX(5) = Box_TU

      HH_TX(6) = HH_TX(1) + HH_TX(2) + HH_TX(3) + HH_TX(4) + HH_TX(5)

      HH_UX(1) = Box_UQ
      HH_UX(2) = Box_US
      HH_UX(3) = Box_UV
      HH_UX(4) = Box_UT
      HH_UX(5) = Box_UU

      HH_UX(6) = HH_UX(1) + HH_UX(2) + HH_UX(3) + HH_UX(4) + HH_UX(5)

!===================================================================================================================================

      end
