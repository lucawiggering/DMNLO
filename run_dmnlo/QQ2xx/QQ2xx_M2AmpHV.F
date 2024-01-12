!**********************************************************************************************************************************

      subroutine QQ2xx_M2AmpHV(mytree,channel_,corrProp_,corrVert_,corrBox_,argxh,argy,AmpTree,AmpProp,AmpVert,AmpBox)

      implicit none

#include "QQ2xx_Model.h"
	

      integer i,mytree,channel_,corrProp_,corrVert_,corrVert_S,corrVert_V,corrVert_T,corrVert_U,corrBox_,argxh,argy
      double precision HV_SX(6),HV_VX(6),HV_TX(6),HV_UX(6),HV_TXa(6),HV_TXb(6),HV_UXa(6),HV_UXb(6),AmpTree,AmpTreeCh(25),AmpProp,AmpVert,AmpBox  

      external  QQ2xx_CorrTXVA,QQ2xx_CorrTXVB, QQ2xx_CorrTXSA, QQ2xx_CorrTXSB, QQ2xx_M2GluonBoxHV, QQ2xx_M2GluinoBoxHV,  QQ2xx_M2GluonTriangleBoxHV, QQ2xx_M2TriangleBoxHV    

!initialize further argxh/argy dependent O(alpha_s) corrections

      call QQ2xx_CorrTXSA(mytree,argxh)
      call QQ2xx_CorrUXSB(mytree,argxh) 

      if((argxh.ne.4).and.(argxh.ne.6)) then                !For Goldstone bosons these corrections are already initialized
      call QQ2xx_CorrTXVB(mytree,argy)
      call QQ2xx_CorrUXVA(mytree,argy)  
      endif

!tree level:

      call HV_SX_(1,argxh,argy,HV_SX)                       !channel_=1,2,3,4,5,6 -> HV_SQ,HV_SS,HV_SV,HV_ST,HV_SU,HV_SX
      call HV_VX_(1,argxh,argy,HV_VX)                       !channel_=1,2,3,4,5,6 -> HV_VQ,HV_VS,HV_HH,HV_VT,HV_VU,HV_VX
      call HV_TX_(1,1,1,argxh,argy,HV_TX)                   !channel_=1,2,3,4,5,6 -> HV_TQ,HV_TS,HV_TV,HV_TT,HV_TU,HV_TX
      call HV_UX_(1,1,1,argxh,argy,HV_UX)                   !channel_=1,2,3,4,5,6 -> HV_UQ,HV_US,HV_UV,HV_UT,HV_UU,HV_UX


      if(channel_.eq.99) then

        AmpTree = HV_SX(6) + HV_VX(6) + HV_TX(6) + HV_UX(6)

      else

        do i = 1, 5

          AmpTreeCh(i) = 0d0
          AmpTreeCh(i+5) = HV_SX(i)
          AmpTreeCh(i+10) = HV_VX(i)
          AmpTreeCh(i+15) = HV_TX(i)
          AmpTreeCh(i+20) = HV_UX(i)

        enddo

        AmpTree = AmpTreeCh(channel_)

      endif


      if(mytree.eq.0) then
                                     
!propagator corrections:

        if(corrProp_.eq.99) then                            !full propagator corrections

          call HV_TX_(1,1,7,argxh,argy,HV_TX)
          call HV_UX_(1,1,7,argxh,argy,HV_UX) 

          AmpProp = HV_TX(6) + HV_UX(6)  

        else                                              

          AmpProp = 0d0

        endif                      

!vertex corrections:

        if(corrVert_.eq.99) then                            !full Vertex corrections

          corrVert_S = 6
          corrVert_V = 8
          corrVert_T = 6
          corrVert_U = 8

          call HV_SX_(corrVert_S,argxh,argy,HV_SX)
          call HV_VX_(corrVert_V,argxh,argy,HV_VX)
          call HV_TX_(corrVert_T,1,1,argxh,argy,HV_TXa)
          call HV_TX_(1,corrVert_U,1,argxh,argy,HV_TXb)
          call HV_UX_(corrVert_U,1,1,argxh,argy,HV_UXa)
          call HV_UX_(1,corrVert_T,1,argxh,argy,HV_UXb)
             

          AmpVert = HV_SX(6) + HV_VX(6) + HV_TXa(6) + HV_TXb(6) + HV_UXa(6) + HV_UXb(6)  


          call QQ2xx_M2GluonTriangleBoxHV(argxh,argy,HV_TX,HV_UX)    !1,2,3,4,5,6->HV_TQ,HV_TS,HV_TV,HV_TT,HV_TU,HV_resTX/HV_UQ,HV_US,HV_UV,HV_UT,HV_UU,HV_resUX

          AmpVert = AmpVert + HV_TX(6) + HV_UX(6)

          if(Looptype0.eq.1) then                           !pure sQuark triangle (O(gs^2)) only added in the QCD-correction case

            call QQ2xx_M2sQuarkTriangleBoxHV(argxh,argy,HV_TX,HV_UX)

            AmpVert = AmpVert + HV_TX(6) + HV_UX(6)

          endif
          
        else if(corrVert_.eq.999) then                      !Sommerfeld corrections
        
          corrVert_S = 2
          corrVert_V = 2

          call HV_SX_(corrVert_S,argxh,argy,HV_SX)
          call HV_VX_(corrVert_V,argxh,argy,HV_VX)            

          AmpVert = HV_SX(6) + HV_VX(6) 

        else

          AmpVert = 0d0   

        endif   
 
!box corrections:

        if(corrBox_.eq.99) then                             !full Box corrections

          call QQ2xx_M2GluonBoxHV(argxh,argy,HV_TX,HV_UX)   !1,2,3,4,5,6->HV_TQ,HV_TS,HV_TV,HV_TT,HV_TU,HV_resTX/HV_UQ,HV_US,HV_UV,HV_UT,HV_UU,HV_resUX

          AmpBox = HV_TX(6) + HV_UX(6)

          call QQ2xx_M2GluinoBoxHV(argxh,argy,HV_TX,HV_UX)

          AmpBox = AmpBox + HV_TX(6) + HV_UX(6)
         
        else if(corrBox_.eq.999) then                       !Sommerfeld corrections
        
          call QQ2xx_M2GluonBoxHV(argxh,argy,HV_TX,HV_UX)   !1,2,3,4,5,6->HV_TQ,HV_TS,HV_TV,HV_TT,HV_TU,HV_resTX/HV_UQ,HV_US,HV_UV,HV_UT,HV_UU,HV_resUX

          AmpBox = HV_TX(6) + HV_UX(6)

        else

          AmpBox = 0d0

        endif

      else

          AmpProp = 0d0

          AmpVert = 0d0

          AmpBox = 0d0
 
      endif                                                                                                                           

!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************


	
      subroutine HV_SX_(corr,argxh,argy,HV_SX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corr,argxh,argy
      double precision Amp_SQ,Amp_SS,MssHSV,Amp_SV,MssHVSV,Amp_ST,Amp_SU,MstHSV,MsuHSV,HV_SX(6)

      Amp_SQ = 0d0
      Amp_SS = 0d0
      Amp_SV = 0d0
      Amp_ST = 0d0
      Amp_SU = 0d0

      do i= sumstartH,sumendH                                   

!==================================================================================================================================

        SS = (1d0,0d0)
        Coupl_A = Coupl_AS(corr,i)
        Coupl_B = JLR(argy,argxh,i)
        propa = DenSmH(i)

!=======================================================Scalar*q===================================================================

        Amp_SQ = 0d0

!=====================================================Scalar*Scalar================================================================

        do j= sumstartH,sumendH

        Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(JLR(argy,argxh,j))
        propac = dconjg(DenSmH(j))
        
	Amp_SS = Amp_SS + MssHSV(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac) 

        enddo  

!=====================================================Scalar*Vector================================================================

        do j=sumstartV,sumendV 

        Coupl_Ac =  dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = dconjg(LLR(j,argy,argxh))
        propac = dconjg(DenSmV(j))
  
        Amp_SV = Amp_SV + MssHVSV(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)   

        enddo

!=======================================================Scalar*t===================================================================

        do j=1,2

        Coupl_Ac = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))
        Coupl_Bc = dconjg(VSfSfu(argy,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
        sfmc = MSf(j,proptypetc,propgentc)         

	Amp_ST = Amp_ST + MstHSV(SS,Coupl_A,Coupl_Bc,Coupl_B,Coupl_Ac,propa,sfmc) 

        enddo

!=======================================================Scalar*u===================================================================

        do j=1,2

        Coupl_Ac = dconjg(VSfSfu(argy,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
        Coupl_Bc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
        sfmc = MSf(j,proptypeuc,propgenuc) 
        
        Amp_SU = Amp_SU + MsuHSV(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,sfmc)

        enddo

!==================================================================================================================================

      enddo

      HV_SX(1) = Amp_SQ
      HV_SX(2) = Amp_SS
      HV_SX(3) = Amp_SV
      HV_SX(4) = Amp_ST
      HV_SX(5) = Amp_SU

      HV_SX(6) = HV_SX(1) + HV_SX(2) + HV_SX(3) + HV_SX(4) + HV_SX(5)

!==================================================================================================================================


      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine HV_VX_(corr,argxh,argy,HV_VX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corr,argxh,argy
      double precision Amp_VQ,Amp_VS,MssVHSV,Amp_VV,MssVSV,Amp_VT,Amp_VU,MstVSV,MsuVSV,HV_VX(6)

      Amp_VQ = 0d0
      Amp_VS = 0d0
      Amp_VV = 0d0
      Amp_VT = 0d0
      Amp_VU = 0d0

      do i=sumstartV,sumendV 

!==================================================================================================================================

        SS = (1d0,0d0)
        Coupl_AP = Coupl_APV(corr,i)
        Coupl_AM = Coupl_AMV(corr,i)
        Coupl_B = LLR(i,argy,argxh)
        propa = DenSmV(i)

!=========================================================Vector*q=================================================================

        Amp_VQ = 0d0

!=======================================================Vector*Scalar==============================================================

        do j=sumstartH,sumendH 

        Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(JLR(argy,argxh,j))
        propac = dconjg(DenSmH(j))
 
        Amp_VS = Amp_VS + MssVHSV(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)                 

        enddo

!=======================================================Vector*Vector==============================================================
    
        do j=sumstartV,sumendV 

        Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = dconjg(LLR(j,argy,argxh))
        propac = dconjg(DenSmV(j))
        
	Amp_VV = Amp_VV + MssVSV(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

        enddo

!========================================================Vector*t==================================================================

        do j=1,2 

        Coupl_Ac = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))
        Coupl_Bc = dconjg(VSfSfu(argy,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
        sfmc =  MSf(j,proptypetc,propgentc)
 
        Amp_VT = Amp_VT + MstVSV(SS,Coupl_AP,Coupl_AM,Coupl_Bc,Coupl_B,Coupl_Ac,propa,sfmc)    

        enddo

!========================================================Vector*u==================================================================

        do j=1,2  

        Coupl_Ac = dconjg(VSfSfu(argy,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
        Coupl_Bc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
        sfmc =  MSf(j,proptypeuc,propgenuc)
  
        Amp_VU = Amp_VU + MsuVSV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,sfmc) 

        enddo

!==================================================================================================================================

      enddo

      HV_VX(1) = Amp_VQ
      HV_VX(2) = Amp_VS
      HV_VX(3) = Amp_VV
      HV_VX(4) = Amp_VT
      HV_VX(5) = Amp_VU

      HV_VX(6) = HV_VX(1) + HV_VX(2) + HV_VX(3) + HV_VX(4) + HV_VX(5)

!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine HV_TX_(corrA,corrB,corrProp_,argxh,argy,HV_TX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corrA,corrB,corrProp_,argxh,argy
      double precision Amp_TQ,Amp_TS,MtsHSV,Amp_TV,MtsVSV,Amp_TT,Amp_TU,MtHV,MtuSV,HV_TX(6)

      Amp_TQ = 0d0
      Amp_TS = 0d0
      Amp_TV = 0d0
      Amp_TT = 0d0
      Amp_TU = 0d0

      do i=1,2 

!==================================================================================================================================

        sfm = MSf(i,proptypet,propgent)
        SS = SST(corrProp_,i)
        Coupl_A =  Coupl_AT(corrA,i)
        Coupl_BP = Coupl_BPT(corrB,i)
        Coupl_BM = Coupl_BMT(corrB,i)

!==========================================================t*q=====================================================================

        Amp_TQ = 0d0

!========================================================t*Scalar==================================================================

        do j=sumstartH,sumendH 

        Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(JLR(argy,argxh,j))
        propac = dconjg(DenSmH(j))

	Amp_TS = Amp_TS + MtsHSV(SS,Coupl_BP,Coupl_BM,Coupl_Ac,Coupl_A,Coupl_Bc,propac,sfm)         

        enddo

!========================================================t*Vector==================================================================
    
        do j=sumstartV,sumendV 

        Coupl_Ac =  dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = dconjg(LLR(j,argy,argxh))
        propac = dconjg(DenSmV(j)) 

	Amp_TV = Amp_TV + MtsVSV(SS,Coupl_BP,Coupl_BM,Coupl_Ac,Coupl_A,Coupl_Bc,propac,sfm)

        enddo

!==========================================================t*t=====================================================================

        do j=1,2  

        Coupl_Ac = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))
        Coupl_Bc = dconjg(VSfSfu(argy,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
        sfmc = MSf(j,proptypetc,propgentc)
        
        Amp_TT = Amp_TT + colorfactorVH_HV*MtHV(SS,Coupl_BP,Coupl_BM,Coupl_Bc,Coupl_A,Coupl_Ac,sfm,sfmc) 

        enddo

!==========================================================t*u=====================================================================

        do j=1,2  

        Coupl_Ac = dconjg(VSfSfu(argy,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
        Coupl_Bc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
        sfmc = MSf(j,proptypeuc,propgenuc)

	Amp_TU = Amp_TU + colorfactorVH_HV*MtuSV(SS,Coupl_BP,Coupl_BM,Coupl_Ac,Coupl_A,Coupl_Bc,sfm,sfmc)

        enddo

!==================================================================================================================================

      enddo

      HV_TX(1) = Amp_TQ
      HV_TX(2) = Amp_TS
      HV_TX(3) = Amp_TV
      HV_TX(4) = Amp_TT
      HV_TX(5) = Amp_TU

      HV_TX(6) = HV_TX(1) + HV_TX(2) + HV_TX(3) + HV_TX(4) + HV_TX(5)

!==================================================================================================================================

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine HV_UX_(corrA,corrB,corrProp_,argxh,argy,HV_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corrA,corrB,corrProp_,argxh,argy
      double precision Amp_UQ,Amp_US,MusHSV,Amp_UV,MusVSV,Amp_UT,Amp_UU,MutSV,MuHV,HV_UX(6)

      external QQ2xx_SSVV,QQ2xx_SSVGl,QQ2xx_VQQQQ,QQ2xx_QQVVp1,QQ2xx_QQVVp2

      Amp_UQ = 0d0
      Amp_US = 0d0
      Amp_UV = 0d0
      Amp_UT = 0d0
      Amp_UU = 0d0

      do i=1,2

!==================================================================================================================================

        sfm = MSf(i,proptypeu,propgenu)
        SS = SSU(corrProp_,i)
        Coupl_AP = Coupl_APU(corrA,i)
        Coupl_AM = Coupl_AMU(corrA,i)
        Coupl_B = Coupl_BU(corrB,i)

!==========================================================u*q=====================================================================

        Amp_UQ = 0d0

!========================================================u*Scalar===================================================================

        do j=sumstartH,sumendH 

        Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(JLR(argy,argxh,j))
        propac = dconjg(DenSmH(j))

	Amp_US = Amp_US + MusHSV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propac,sfm)            

        enddo

!=======================================================u*Vector===================================================================
    
        do j=sumstartV,sumendV 

        Coupl_Ac =  dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = dconjg(LLR(j,argy,argxh))
        propac = dconjg(DenSmV(j))

        Amp_UV = Amp_UV + MusVSV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propac,sfm)

        enddo

!==========================================================u*t=====================================================================

        do j=1,2  

        Coupl_Ac = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))
        Coupl_Bc = dconjg(VSfSfu(argy,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
        sfmc = MSf(j,proptypetc,propgentc)

	Amp_UT = Amp_UT + colorfactorVH_HV*MutSV(SS,Coupl_AP,Coupl_AM,Coupl_Bc,Coupl_B,Coupl_Ac,sfm,sfmc)

        enddo

!==========================================================u*u=====================================================================

        do j=1,2

        Coupl_Ac = dconjg(VSfSfu(argy,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
        Coupl_Bc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
        sfmc = MSf(j,proptypeuc,propgenuc)

	Amp_UU = Amp_UU + colorfactorVH_HV*MuHV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,sfm,sfmc)

        enddo

!==================================================================================================================================

      enddo

      HV_UX(1) = Amp_UQ
      HV_UX(2) = Amp_US
      HV_UX(3) = Amp_UV
      HV_UX(4) = Amp_UT
      HV_UX(5) = Amp_UU

      HV_UX(6) = HV_UX(1) + HV_UX(2) + HV_UX(3) + HV_UX(4) + HV_UX(5)

!==================================================================================================================================

      end


