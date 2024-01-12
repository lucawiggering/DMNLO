!**********************************************************************************************************************************

      subroutine QQ2xx_M2AmpHH(mytree,channel_,corrProp_,corrVert_,corrBox_,argxh,argyh,AmpTree,AmpProp,AmpVert,AmpBox)

      implicit none

#include "QQ2xx_Model.h"
	

      integer i,mytree,channel_,corrProp_,corrVert_,corrVert_Q,corrVert_S,corrVert_V,corrVert_T,corrVert_U,corrBox_,argxh,argyh
      double precision HH_QX(6),HH_SX(6),HH_VX(6),HH_TX(6),HH_UX(6),HH_TXa(6),HH_TXb(6),HH_UXa(6),HH_UXb(6),AmpTree,AmpTreeCh(25),AmpProp,AmpVert,AmpBox 

      external  QQ2xx_QXS, QQ2xx_CorrTXSA, QQ2xx_CorrTXSB, QQ2xx_CorrUXSA, QQ2xx_CorrUXSB, QQ2xx_M2GluonBoxHH, QQ2xx_M2GluinoBoxHH, QQ2xx_M2TriangleBoxHH    

!initialize further argxh/argyh dependent O(alpha_s) corrections

      call QQ2xx_CorrQXS(mytree,argxh,argyh)
      call QQ2xx_CorrTXSA(mytree,argxh)
      call QQ2xx_CorrUXSB(mytree,argxh)                  
      call QQ2xx_CorrTXSB(mytree,argyh)            
      call QQ2xx_CorrUXSA(mytree,argyh)

!tree level:

      call HH_QX_(1,argxh,argyh,HH_QX)                      !channel_=1,2,3,4,5,6 -> HH_QQ,HH_QS,HH_QV,HH_QT,HH_QU,HH_QX
      call HH_SX_(1,argxh,argyh,HH_SX)                      !channel_=1,2,3,4,5,6 -> HH_SQ,HH_SS,HH_SV,HH_ST,HH_SU,HH_SX
      call HH_VX_(1,argxh,argyh,HH_VX)                      !channel_=1,2,3,4,5,6 -> HH_VQ,HH_VS,HH_HH,HH_VT,HH_VU,HH_VX
      call HH_TX_(1,1,1,argxh,argyh,HH_TX)                  !channel_=1,2,3,4,5,6 -> HH_TQ,HH_TS,HH_TV,HH_TT,HH_TU,HH_TX
      call HH_UX_(1,1,1,argxh,argyh,HH_UX)                  !channel_=1,2,3,4,5,6 -> HH_UQ,HH_US,HH_UV,HH_UT,HH_UU,HH_UX

      if(channel_.eq.99) then

        AmpTree = HH_QX(6) + HH_SX(6) + HH_VX(6) + HH_TX(6) + HH_UX(6)

      else

        do i = 1, 5

          AmpTreeCh(i) = HH_QX(i)
          AmpTreeCh(i+5) = HH_SX(i)
          AmpTreeCh(i+10) = HH_VX(i)
          AmpTreeCh(i+15) = HH_TX(i)
          AmpTreeCh(i+20) = HH_UX(i)

        enddo

        AmpTree = AmpTreeCh(channel_)

      endif


      if(mytree.eq.0) then
                                     
!propagator corrections:

        if(corrProp_.eq.99) then                            !full propagator corrections

          call HH_TX_(1,1,7,argxh,argyh,HH_TX)
          call HH_UX_(1,1,7,argxh,argyh,HH_UX) 

          AmpProp = HH_TX(6) + HH_UX(6)              

        else                                              

          AmpProp = 0d0

        endif                      

!vertex corrections:

        if(corrVert_.eq.99) then                            !full Vertex corrections

          corrVert_Q = 5
          corrVert_S = 6
          corrVert_V = 8
          corrVert_T = 6
          corrVert_U = 6

          call HH_QX_(corrVert_Q,argxh,argyh,HH_QX)
          call HH_SX_(corrVert_S,argxh,argyh,HH_SX)
          call HH_VX_(corrVert_V,argxh,argyh,HH_VX)
          call HH_TX_(corrVert_T,1,1,argxh,argyh,HH_TXa)
          call HH_TX_(1,corrVert_T,1,argxh,argyh,HH_TXb)
          call HH_UX_(corrVert_U,1,1,argxh,argyh,HH_UXa)
          call HH_UX_(1,corrVert_U,1,argxh,argyh,HH_UXb)
             

          AmpVert = HH_QX(6) + HH_SX(6) + HH_VX(6) +  HH_TXa(6) + HH_TXb(6) + HH_UXa(6) + HH_UXb(6)  


          if(Looptype0.eq.1) then                           !pure sQuark triangle (O(gs^2)) only added in the QCD-correction case (correction to the squark-squark-Higgs-Higgs-Vertex)

            call QQ2xx_M2TriangleBoxHH(argxh,argyh,HH_TX,HH_UX)   !1,2,3,4,5,6->HH_TQ,HH_TS,HH_TV,HH_TT,HH_TU,HH_TX/HH_UQ,HH_US,HH_UV,HH_UT,HH_UU,HH_UX

            AmpVert = AmpVert + HH_TX(6) + HH_UX(6)

          endif
          
        else if(corrVert_.eq.999) then                      !Sommerfeld corrections
        
          corrVert_Q = 2
          corrVert_S = 2
          corrVert_V = 2

          call HH_QX_(corrVert_Q,argxh,argyh,HH_QX)
          call HH_SX_(corrVert_S,argxh,argyh,HH_SX)
          call HH_VX_(corrVert_V,argxh,argyh,HH_VX)
          
          AmpVert = HH_QX(6) + HH_SX(6) + HH_VX(6) 

        else

          AmpVert = 0d0   

        endif   
 
!box corrections:

        if(corrBox_.eq.99) then                             !full Box corrections

          call QQ2xx_M2GluonBoxHH(argxh,argyh,HH_TX,HH_UX)  !1,2,3,4,5,6->HH_TQ,HH_TS,HH_TV,HH_TT,HH_TU,HH_TX/HH_UQ,HH_US,HH_UV,HH_UT,HH_UU,HH_UX

          AmpBox =  HH_TX(6) + HH_UX(6)

          call QQ2xx_M2GluinoBoxHH(argxh,argyh,HH_TX,HH_UX)

          AmpBox = AmpBox + HH_TX(6) + HH_UX(6)
          
        else if(corrBox_.eq.999) then                       !Sommerfeld corrections
        
          call QQ2xx_M2GluonBoxHH(argxh,argyh,HH_TX,HH_UX) 
          
          AmpBox =  HH_TX(6) + HH_UX(6)

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



      subroutine HH_QX_(corr,argxh,argyh,HH_QX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer j,corr,argxh,argyh
      double precision Amp_QQ,MqqSS,Amp_QS,MqsHSS,Amp_QV,MqsVSS,Amp_QT,Amp_QU,MtuqSS,HH_QX(6)

      Amp_QS = 0d0
      Amp_QV = 0d0
      Amp_QT = 0d0
      Amp_QU = 0d0

!==================================================================================================================================

      SS = (1d0,0d0)
      Coupl_A = Coupl_AQS(corr) 

!==========================================================q*q=====================================================================

      Coupl_Ac = dconjg(HHSfSf(argxh,argyh,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))     
 
      Amp_QQ = MqqSS(Coupl_A,Coupl_Ac)

!========================================================q*Scalar==================================================================

      do j= sumstartH,sumendH

        Coupl_Ac =  dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(ILR(j,argxh,argyh))
        propac = dconjg(DenSmH(j))        

 	Amp_QS = Amp_QS + MqsHSS(Coupl_A,Coupl_Ac,Coupl_Bc,propac)  

      enddo

!========================================================q*Vector==================================================================

      do j= sumstartV,sumendV

        Coupl_Ac =  dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = dconjg(JLR(j,argyh,argxh))
        propac = dconjg(DenSmV(j))        

	Amp_QV = Amp_QV + MqsVSS(Coupl_A,Coupl_Ac,Coupl_Bc,propac)   

      enddo

!==========================================================q*t=====================================================================      

      do j= 1,2

        Coupl_Ac =  dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))
        Coupl_Bc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
        sfmc = MSf(j,proptypetc,propgentc)       

	Amp_QT = Amp_QT + MtuqSS(SS,Coupl_A,Coupl_Ac,Coupl_Bc,tman,sfmc)

      enddo 

!==========================================================q*u=====================================================================       

      do j= 1,2

        Coupl_Ac = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
        Coupl_Bc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
        sfmc = MSf(j,proptypeuc,propgenuc)     

        Amp_QU = Amp_QU + MtuqSS(SS,Coupl_A,Coupl_Ac,Coupl_Bc,uman,sfmc) 

      enddo 

!==================================================================================================================================

      HH_QX(1) = Amp_QQ 
      HH_QX(2) = Amp_QS
      HH_QX(3) = Amp_QV
      HH_QX(4) = Amp_QT
      HH_QX(5) = Amp_QU

      HH_QX(6) = HH_QX(1) + HH_QX(2) + HH_QX(3) + HH_QX(4) + HH_QX(5)

!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************


	
      subroutine HH_SX_(corr,argxh,argyh,HH_SX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corr,argxh,argyh
      double precision Amp_SQ,MsqHSS,Amp_SS,MssHSS,Amp_SV,MssHVSS,Amp_ST,Amp_SU,MstuHSS,HH_SX(6)

      Amp_SQ = 0d0
      Amp_SS = 0d0
      Amp_SV = 0d0
      Amp_ST = 0d0
      Amp_SU = 0d0

      do i= sumstartH,sumendH                                   

!==================================================================================================================================

      SS = (1d0,0d0)
      Coupl_A = Coupl_AS(corr,i)
      Coupl_B = ILR(i,argxh,argyh)
      propa = DenSmH(i)

!=========================================================Scalar*q=================================================================

        Coupl_Ac = dconjg(HHSfSf(argxh,argyh,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))   

	Amp_SQ = Amp_SQ + MsqHSS(Coupl_A,Coupl_Ac,Coupl_B,propa) 

!=======================================================Scalar*Scalar==============================================================

        do j= sumstartH,sumendH

          Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl_Bc = dconjg(ILR(j,argxh,argyh))
          propac = dconjg(DenSmH(j))        

	  Amp_SS = Amp_SS + MssHSS(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)   

        enddo  

!=======================================================Scalar*Vector==============================================================

        do j=sumstartV,sumendV 

          Coupl_Ac =  dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl_Bc = dconjg(JLR(j,argyh,argxh))
          propac = dconjg(DenSmV(j)) 
 
          Amp_SV = Amp_SV + MssHVSS(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)          

        enddo

!=========================================================Scalar*t==================================================================

        do j=1,2

          Coupl_Ac = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))
          Coupl_Bc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
          sfmc = MSf(j,proptypetc,propgentc)  

	  Amp_ST = Amp_ST + MstuHSS(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,tman,sfmc) 

        enddo

!=========================================================Scalar*u=================================================================

        do j=1,2

          Coupl_Ac = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
          Coupl_Bc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
          sfmc = MSf(j,proptypeuc,propgenuc)  

	  Amp_SU = Amp_SU + MstuHSS(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,uman,sfmc) 

        enddo

!==================================================================================================================================

      enddo

      HH_SX(1) = Amp_SQ
      HH_SX(2) = Amp_SS
      HH_SX(3) = Amp_SV
      HH_SX(4) = Amp_ST
      HH_SX(5) = Amp_SU

      HH_SX(6) = HH_SX(1) + HH_SX(2) + HH_SX(3) + HH_SX(4) + HH_SX(5)

!==================================================================================================================================


      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine HH_VX_(corr,argxh,argyh,HH_VX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corr,argxh,argyh
      double precision Amp_VQ,MsqVSS,Amp_VS,MssVHSS,Amp_VV,MssVSS,Amp_VT,Amp_VU,MstuVSS,HH_VX(6)

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
        Coupl_B = JLR(i,argyh,argxh)
        propa = DenSmV(i)

!==========================================================Vector*q================================================================

        Coupl_Ac = dconjg(HHSfSf(argxh,argyh,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1)) 

 
        Amp_VQ = Amp_VQ + MsqVSS(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,propa) 

!========================================================Vector*Scalar=============================================================

        do j=sumstartH,sumendH 

          Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl_Bc = dconjg(ILR(j,argxh,argyh))
          propac = dconjg(DenSmH(j)) 
 
          Amp_VS = Amp_VS + MssVHSS(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)               

        enddo

!========================================================Vector-Vector=============================================================
    
        do j=sumstartV,sumendV 
 
          Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl_Bc = dconjg(JLR(j,argyh,argxh))
          propac = dconjg(DenSmV(j))        

	  Amp_VV = Amp_VV + MssVSS(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)   

        enddo

!========================================================Vector*t==================================================================

        do j=1,2  

          Coupl_Ac = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))
          Coupl_Bc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
          sfmc =  MSf(j,proptypetc,propgentc) 
 
          Amp_VT = Amp_VT + MstuVSS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,tman,sfmc) 

        enddo

!========================================================Vector*u==================================================================

        do j=1,2  

          Coupl_Ac = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
          Coupl_Bc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
          sfmc =  MSf(j,proptypeuc,propgenuc) 
 
          Amp_VU = Amp_VU + MstuVSS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,uman,sfmc)  

        enddo

!==================================================================================================================================

      enddo

      HH_VX(1) = Amp_VQ
      HH_VX(2) = Amp_VS
      HH_VX(3) = Amp_VV
      HH_VX(4) = Amp_VT
      HH_VX(5) = Amp_VU

      HH_VX(6) = HH_VX(1) + HH_VX(2) + HH_VX(3) + HH_VX(4) + HH_VX(5)

!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine HH_TX_(corrA,corrB,corrProp_,argxh,argyh,HH_TX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corrA,corrB,corrProp_,argxh,argyh
      double precision Amp_TQ,MtuqSS,Amp_TS,MstuHSS,Amp_TV,MtusVSS,Amp_TT,Amp_TU,MtuSS,HH_TX(6)

      Amp_TQ = 0d0
      Amp_TS = 0d0
      Amp_TV = 0d0
      Amp_TT = 0d0
      Amp_TU = 0d0

      do i=1,2 

!==================================================================================================================================

        sfm = MSf(i,proptypet,propgent)
        SS = SST(corrProp_,i)
        Coupl_A = Coupl_AT(corrA,i)
        Coupl_B = Coupl_BT(corrB,i)

!==========================================================t*q=====================================================================

        Coupl_Ac = dconjg(HHSfSf(argxh,argyh,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))

	Amp_TQ = Amp_TQ + MtuqSS(SS,Coupl_A,Coupl_Ac,Coupl_B,tman,sfm)     

!========================================================t*Scalar==================================================================

        do j=sumstartH,sumendH 

          Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl_Bc = dconjg(ILR(j,argxh,argyh))
          propac = dconjg(DenSmH(j)) 

	  Amp_TS = Amp_TS + MstuHSS(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propac,tman,sfm)        

        enddo

!========================================================t*Vector==================================================================
    
        do j=sumstartV,sumendV 

          Coupl_Ac =  dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))   
          Coupl_Bc = dconjg(JLR(j,argyh,argxh))
          propac = dconjg(DenSmV(j))   

          Amp_TV = Amp_TV + MtusVSS(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propac,tman,sfm)    

        enddo

!========================================================t*t=======================================================================

        do j=1,2  

          Coupl_Ac = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))
          Coupl_Bc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
          sfmc = MSf(j,proptypetc,propgentc)        

	  Amp_TT = Amp_TT + MtuSS(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,tman,tman,sfm,sfmc)  

        enddo

!========================================================t*u=======================================================================

        do j=1,2  

          Coupl_Ac = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
          Coupl_Bc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
          sfmc = MSf(j,proptypeuc,propgenuc)

	  Amp_TU = Amp_TU + MtuSS(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,tman,uman,sfm,sfmc) 

        enddo

!==================================================================================================================================

      enddo

      HH_TX(1) = Amp_TQ
      HH_TX(2) = Amp_TS
      HH_TX(3) = Amp_TV
      HH_TX(4) = Amp_TT
      HH_TX(5) = Amp_TU

      HH_TX(6) = HH_TX(1) + HH_TX(2) + HH_TX(3) + HH_TX(4) + HH_TX(5)

!==================================================================================================================================

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine HH_UX_(corrA,corrB,corrProp_,argxh,argyh,HH_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corr,corrA,corrB,corrProp_,argxh,argyh
      double precision Amp_UQ,MtuqSS,Amp_US,MstuHSS,Amp_UV,MtusVSS,Amp_UT,Amp_UU,MtuSS,HH_UX(6)

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
        Coupl_B = Coupl_BU(corrB,i)

!==========================================================u*q=====================================================================

        Coupl_Ac = dconjg(HHSfSf(argxh,argyh,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))    

	Amp_UQ = Amp_UQ + MtuqSS(SS,Coupl_A,Coupl_Ac,Coupl_B,uman,sfm)   

!========================================================u*Scalar==================================================================

        do j=sumstartH,sumendH 

          Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl_Bc = dconjg(ILR(j,argxh,argyh))
          propac = dconjg(DenSmH(j))
        
          Amp_US = Amp_US + MstuHSS(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propac,uman,sfm)              

        enddo

!========================================================u*Vector==================================================================
    
        do j=sumstartV,sumendV 

          Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl_Bc = dconjg(JLR(j,argyh,argxh))
          propac = dconjg(DenSmV(j))     

	  Amp_UV = Amp_UV + MtusVSS(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propac,uman,sfm)    

        enddo

!========================================================u*t=======================================================================

        do j=1,2  

          Coupl_Ac = dconjg(HSfSft(argxh,sfmeigen1,j,sftype1,proptypetc,sfgen1,propgentc))
          Coupl_Bc = dconjg(HSfSfu(argyh,j,sfmeigen2,proptypetc,sftype2,propgentc,sfgen2))
          sfmc = MSf(j,proptypetc,propgentc)        

	  Amp_UT = Amp_UT + MtuSS(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,uman,tman,sfm,sfmc)   

        enddo

!========================================================u*u=======================================================================

        do j=1,2  

          Coupl_Ac = dconjg(HSfSfu(argyh,sfmeigen1,j,sftype1,proptypeuc,sfgen1,propgenuc))
          Coupl_Bc = dconjg(HSfSft(argxh,j,sfmeigen2,proptypeuc,sftype2,propgenuc,sfgen2))
          sfmc = MSf(j,proptypeuc,propgenuc)        

	  Amp_UU = Amp_UU + MtuSS(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,uman,uman,sfm,sfmc) 

        enddo

!==================================================================================================================================

      enddo

      HH_UX(1) = Amp_UQ
      HH_UX(2) = Amp_US
      HH_UX(3) = Amp_UV
      HH_UX(4) = Amp_UT
      HH_UX(5) = Amp_UU

      HH_UX(6) = HH_UX(1) + HH_UX(2) + HH_UX(3) + HH_UX(4) + HH_UX(5)

!==================================================================================================================================

      end


