!**********************************************************************************************************************************

      subroutine QQ2xx_M2AmpVH(mytree,channel_,corrProp_,corrVert_,corrBox_,argx,argyh,AmpTree,AmpProp,AmpVert,AmpBox)

      implicit none

#include "QQ2xx_Model.h"
	

      integer i,mytree,channel_,corrProp_,corrVert_,corrVert_S,corrVert_V,corrVert_T,corrVert_U,corrBox_,argx,argyh
      double precision VH_SX(6),VH_VX(6),VH_TX(6),VH_UX(6),VH_TXa(6),VH_TXb(6),VH_UXa(6),VH_UXb(6),AmpTree,AmpTreeCh(25),AmpProp,AmpVert,AmpBox  

      external  QQ2xx_CorrTXVA,QQ2xx_CorrTXVB, QQ2xx_CorrTXSA, QQ2xx_CorrTXSB, QQ2xx_M2GluonBoxVH, QQ2xx_M2GluinoBoxVH, QQ2xx_M2GluonTriangleBoxVH, QQ2xx_M2TriangleBoxVH    

!initialize further argx/argyh dependent O(alpha_s) corrections

      call QQ2xx_CorrTXSB(mytree,argyh)
      call QQ2xx_CorrUXSA(mytree,argyh) 

      if((argyh.ne.4).and.(argyh.ne.6)) then                !For Goldstone bosons these corrections are already initialized
      call QQ2xx_CorrTXVA(mytree,argx)
      call QQ2xx_CorrUXVB(mytree,argx)  
      endif

!tree level:

      call VH_SX_(1,argx,argyh,VH_SX)                       !channel_=1,2,3,4,5,6 -> VH_SQ,VH_SS,VH_SV,VH_ST,VH_SU,VH_SX
      call VH_VX_(1,argx,argyh,VH_VX)                       !channel_=1,2,3,4,5,6 -> VH_VQ,VH_VS,VH_HH,VH_VT,VH_VU,VH_VX
      call VH_TX_(1,1,1,argx,argyh,VH_TX)                   !channel_=1,2,3,4,5,6 -> VH_TQ,VH_TS,VH_TV,VH_TT,VH_TU,VH_TX
      call VH_UX_(1,1,1,argx,argyh,VH_UX)                   !channel_=1,2,3,4,5,6 -> VH_UQ,VH_US,VH_UV,VH_UT,VH_UU,VH_UX


      if(channel_.eq.99) then

        AmpTree = VH_SX(6) + VH_VX(6) + VH_TX(6) + VH_UX(6)

      else

        do i = 1, 5

          AmpTreeCh(i) = 0d0
          AmpTreeCh(i+5) = VH_SX(i)
          AmpTreeCh(i+10) = VH_VX(i)
          AmpTreeCh(i+15) = VH_TX(i)
          AmpTreeCh(i+20) = VH_UX(i)

        enddo

        AmpTree = AmpTreeCh(channel_)

      endif


      if(mytree.eq.0) then
                                     
!propagator corrections:

        if(corrProp_.eq.99) then                            !full propagator corrections

          call VH_TX_(1,1,7,argx,argyh,VH_TX)
          call VH_UX_(1,1,7,argx,argyh,VH_UX) 

          AmpProp = VH_TX(6) + VH_UX(6)

        else                                              

          AmpProp = 0d0

        endif                      

!vertex corrections:

        if(corrVert_.eq.99) then                            !full Vertex corrections

          corrVert_S = 6
          corrVert_V = 8
          corrVert_T = 8
          corrVert_U = 6

          call VH_SX_(corrVert_S,argx,argyh,VH_SX)
          call VH_VX_(corrVert_V,argx,argyh,VH_VX)
          call VH_TX_(corrVert_T,1,1,argx,argyh,VH_TXa)
          call VH_TX_(1,corrVert_U,1,argx,argyh,VH_TXb)
          call VH_UX_(corrVert_U,1,1,argx,argyh,VH_UXa)
          call VH_UX_(1,corrVert_T,1,argx,argyh,VH_UXb)
             

          AmpVert = VH_SX(6) + VH_VX(6) + VH_TXa(6) + VH_TXb(6) + VH_UXa(6) + VH_UXb(6)  


          call QQ2xx_M2GluonTriangleBoxVH(argx,argyh,VH_TX,VH_UX)   !1,2,3,4,5,6->VH_TQ,VH_TS,VH_TV,VH_TT,VH_TU,VH_TX/VH_UQ,VH_US,VH_UV,VH_UT,VH_UU,VH_UX

          AmpVert = AmpVert + VH_TX(6) + VH_UX(6)

          if(Looptype0.eq.1) then                           !pure sQuark triangle (O(gs^2)) only added in the QCD-correction case

            call QQ2xx_M2sQuarkTriangleBoxVH(argx,argyh,VH_TX,VH_UX)

            AmpVert = AmpVert + VH_TX(6) + VH_UX(6)

          endif
          
        else if(corrVert_.eq.999) then                      !Sommerfeld corrections
        
          corrVert_S = 2
          corrVert_V = 2

          call VH_SX_(corrVert_S,argx,argyh,VH_SX)
          call VH_VX_(corrVert_V,argx,argyh,VH_VX)
           
          AmpVert = VH_SX(6) + VH_VX(6)

        else

          AmpVert = 0d0   

        endif   
 
!box corrections:

        if(corrBox_.eq.99) then                             !full Box corrections

          call QQ2xx_M2GluonBoxVH(argx,argyh,VH_TX,VH_UX)   !1,2,3,4,5,6->VH_TQ,VH_TS,VH_TV,VH_TT,VH_TU,VH_TX/VH_UQ,VH_US,VH_UV,VH_UT,VH_UU,VH_UX

          AmpBox = VH_TX(6) + VH_UX(6)

          call QQ2xx_M2GluinoBoxVH(argx,argyh,VH_TX,VH_UX)

          AmpBox = AmpBox + VH_TX(6) + VH_UX(6)
          
        else if(corrBox_.eq.999) then                       !Sommerfeld corrections
        
          call QQ2xx_M2GluonBoxVH(argx,argyh,VH_TX,VH_UX)   !1,2,3,4,5,6->VH_TQ,VH_TS,VH_TV,VH_TT,VH_TU,VH_TX/VH_UQ,VH_US,VH_UV,VH_UT,VH_UU,VH_UX

          AmpBox = VH_TX(6) + VH_UX(6)

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


	
      subroutine VH_SX_(corr,argx,argyh,VH_SX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corr,argx,argyh
      double precision Amp_SQ,Amp_SS,MssHVS,Amp_SV,MssHVVS,Amp_ST,Amp_SU,MstHVS,MsuHVS,VH_SX(6)

      Amp_SQ = 0d0
      Amp_SS = 0d0
      Amp_SV = 0d0
      Amp_ST = 0d0
      Amp_SU = 0d0

      do i= sumstartH,sumendH                                   

!==================================================================================================================================

        SS = (1d0,0d0)
        Coupl_A = Coupl_AS(corr,i)
        Coupl_B = JLR(argx,i,argyh)
        propa = DenSmH(i)

!=======================================================Scalar*q===================================================================

        Amp_SQ = 0d0

!=====================================================Scalar*Scalar================================================================

        do j= sumstartH,sumendH

          Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl_Bc = dconjg(JLR(argx,j,argyh))
          propac = dconjg(DenSmH(j))
        
	  Amp_SS = Amp_SS + MssHVS(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)  

        enddo  

!=====================================================Scalar*Vector================================================================

        do j=sumstartV,sumendV 

          Coupl_Ac =  dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl_Bc = dconjg(LLR(j,argx,argyh))
          propac = dconjg(DenSmV(j))

          Amp_SV = Amp_SV + MssHVVS(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)          

        enddo

!=======================================================Scalar*t===================================================================

        do j=1,2

          Coupl_Ac = dconjg(VSfSft(argx,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))
          Coupl_Bc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
          sfmc = MSf(j,proptypetc,propgentc) 
         
	  Amp_ST= Amp_ST + MstHVS(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,sfmc) 

        enddo

!=======================================================Scalar*u===================================================================

        do j=1,2

          Coupl_Ac = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
          Coupl_Bc = dconjg(VSfSft(argx,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
          sfmc = MSf(j,proptypeuc,propgenuc) 

	  Amp_SU = Amp_SU + MsuHVS(SS,Coupl_A,Coupl_Bc,Coupl_B,Coupl_Ac,propa,sfmc) 

        enddo

!==================================================================================================================================

      enddo

      VH_SX(1) = Amp_SQ
      VH_SX(2) = Amp_SS
      VH_SX(3) = Amp_SV
      VH_SX(4) = Amp_ST
      VH_SX(5) = Amp_SU

      VH_SX(6) = VH_SX(1) + VH_SX(2) + VH_SX(3) + VH_SX(4) + VH_SX(5)

!==================================================================================================================================


      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine VH_VX_(corr,argx,argyh,VH_VX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corr,argx,argyh
      double precision Amp_VQ,Amp_VS,MssVVS,Amp_VV,MssVHVS,Amp_VT,Amp_VU,MstVVS,MsuVVS,VH_VX(6)

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
        Coupl_B = LLR(i,argx,argyh)
        propa = DenSmV(i)

!=========================================================Vector*q=================================================================

        Amp_VQ = 0d0

!=======================================================Vector*Scalar==============================================================

        do j=sumstartH,sumendH 

          Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl_Bc = dconjg(JLR(argx,j,argyh))
          propac = dconjg(DenSmH(j)) 
 
          Amp_VS = Amp_VS + MssVHVS(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)               

        enddo

!=======================================================Vector*Vector==============================================================
    
        do j=sumstartV,sumendV 

          Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl_Bc = dconjg(LLR(j,argx,argyh))
          propac = dconjg(DenSmV(j))
        
	  Amp_VV = Amp_VV + MssVVS(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac) 

        enddo

!========================================================Vector*t==================================================================

        do j=1,2  

          Coupl_Ac = dconjg(VSfSft(argx,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))
          Coupl_Bc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
          sfmc =  MSf(j,proptypetc,propgentc)

          Amp_VT = Amp_VT + MstVVS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,sfmc) 

        enddo

!========================================================Vector*u==================================================================

        do j=1,2  

          Coupl_Ac = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
          Coupl_Bc = dconjg(VSfSft(argx,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
          sfmc =  MSf(j,proptypeuc,propgenuc)  

          Amp_VU = Amp_VU + MsuVVS(SS,Coupl_AP,Coupl_AM,Coupl_Bc,Coupl_B,Coupl_Ac,propa,sfmc)  

        enddo

!==================================================================================================================================

      enddo

      VH_VX(1) = Amp_VQ
      VH_VX(2) = Amp_VS
      VH_VX(3) = Amp_VV
      VH_VX(4) = Amp_VT
      VH_VX(5) = Amp_VU

      VH_VX(6) = VH_VX(1) + VH_VX(2) + VH_VX(3) + VH_VX(4) + VH_VX(5)

!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine VH_TX_(corrA,corrB,corrProp_,argx,argyh,VH_TX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corrA,corrB,corrProp_,argx,argyh
      double precision Amp_TQ,Amp_TS,MtsHVS,Amp_TV,MtsVVS,Amp_TT,Amp_TU,MtVS,MtuVS,VH_TX(6)

      external QQ2xx_SSVV,QQ2xx_SSVGl,QQ2xx_VQQQQ,QQ2xx_QQVVp1,QQ2xx_QQVVp2

      Amp_TQ = 0d0
      Amp_TS = 0d0
      Amp_TV = 0d0
      Amp_TT = 0d0
      Amp_TU = 0d0

      do i=1,2 

!==================================================================================================================================

        sfm = MSf(i,proptypet,propgent)
        SS = SST(corrProp_,i)
        Coupl_AP = Coupl_APT(corrA,i)
        Coupl_AM = Coupl_AMT(corrA,i)
        Coupl_B = Coupl_BT(corrB,i)

!==========================================================t*q=====================================================================

        Amp_TQ = 0d0

!========================================================t*Scalar==================================================================

        do j=sumstartH,sumendH 

          Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl_Bc = dconjg(JLR(argx,j,argyh))
          propac = dconjg(DenSmH(j))

          Amp_TS = Amp_TS + MtsHVS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propac,sfm)            

        enddo

!========================================================t*Vector==================================================================
    
        do j=sumstartV,sumendV 

          Coupl_Ac =  dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl_Bc = dconjg(LLR(j,argx,argyh))
          propac = dconjg(DenSmV(j))

          Amp_TV = Amp_TV + MtsVVS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propac,sfm)  

        enddo

!==========================================================t*t=====================================================================

        do j=1,2  

          Coupl_Ac = dconjg(VSfSft(argx,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))
          Coupl_Bc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
          sfmc = MSf(j,proptypetc,propgentc)

	  Amp_TT = Amp_TT + colorfactorVH_HV*MtVS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,sfm,sfmc) 

        enddo

!==========================================================t*u=====================================================================

        do j=1,2  

          Coupl_Ac = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
          Coupl_Bc = dconjg(VSfSft(argx,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
          sfmc = MSf(j,proptypeuc,propgenuc)

	  Amp_TU = Amp_TU + colorfactorVH_HV*MtuVS(SS,Coupl_AP,Coupl_AM,Coupl_Bc,Coupl_B,Coupl_Ac,sfm,sfmc)

        enddo

!==================================================================================================================================

      enddo

      VH_TX(1) = Amp_TQ
      VH_TX(2) = Amp_TS
      VH_TX(3) = Amp_TV
      VH_TX(4) = Amp_TT
      VH_TX(5) = Amp_TU

      VH_TX(6) = VH_TX(1) + VH_TX(2) + VH_TX(3) + VH_TX(4) + VH_TX(5)

!==================================================================================================================================

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine VH_UX_(corrA,corrB,corrProp_,argx,argyh,VH_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corrA,corrB,corrProp_,argx,argyh
      double precision Amp_UQ,Amp_US,MusHVS,Amp_UV,MusVVS,Amp_UT,Amp_UU,MutVS,MuVS,VH_UX(6)

      Amp_UQ = 0d0
      Amp_US = 0d0
      Amp_UV = 0d0
      Amp_UT = 0d0
      Amp_UU = 0d0

      do i=1,2

!==================================================================================================================================

        sfm = MSf(i,proptypeu,propgenu)
        SS = SSU(corrProp_,i)
        Coupl_A = Coupl_AU(corrA,i)
        Coupl_BP = Coupl_BPU(corrB,i)
        Coupl_BM = Coupl_BMU(corrB,i)

!==========================================================u*q=====================================================================

        Amp_UQ = 0d0

!========================================================u*Scalar===================================================================

        do j=sumstartH,sumendH 

          Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl_Bc = dconjg(JLR(argx,j,argyh))
          propac = dconjg(DenSmH(j))

          Amp_US = Amp_US + MusHVS(SS,Coupl_BP,Coupl_BM,Coupl_Ac,Coupl_A,Coupl_Bc,propac,sfm)              

        enddo

!=======================================================u*Vector===================================================================
    
        do j=sumstartV,sumendV 

          Coupl_Ac =  dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl_Bc = dconjg(LLR(j,argx,argyh))
          propac = dconjg(DenSmV(j))

	  Amp_UV = Amp_UV + MusVVS(SS,Coupl_BP,Coupl_BM,Coupl_Ac,Coupl_A,Coupl_Bc,propac,sfm)  

        enddo

!==========================================================u*t=====================================================================

        do j=1,2  

          Coupl_Ac = dconjg(VSfSft(argx,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))
          Coupl_Bc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
          sfmc = MSf(j,proptypetc,propgentc)

          Amp_UT = Amp_UT + colorfactorVH_HV*MutVS(SS,Coupl_BP,Coupl_BM,Coupl_Ac,Coupl_A,Coupl_Bc,sfm,sfmc)

        enddo

!==========================================================u*u=====================================================================

        do j=1,2

          Coupl_Ac = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
          Coupl_Bc = dconjg(VSfSft(argx,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
          sfmc = MSf(j,proptypeuc,propgenuc) 

	  Amp_UU = Amp_UU + colorfactorVH_HV*MuVS(SS,Coupl_BP,Coupl_BM,Coupl_Bc,Coupl_A,Coupl_Ac,sfm,sfmc) 

        enddo

!==================================================================================================================================

      enddo

      VH_UX(1) = Amp_UQ
      VH_UX(2) = Amp_US
      VH_UX(3) = Amp_UV
      VH_UX(4) = Amp_UT
      VH_UX(5) = Amp_UU

      VH_UX(6) = VH_UX(1) + VH_UX(2) + VH_UX(3) + VH_UX(4) + VH_UX(5)

!==================================================================================================================================

      end


