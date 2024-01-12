!**********************************************************************************************************************************

      subroutine QQ2xx_M2AmpVV(mytree,channel_,corrProp_,corrVert_,corrBox_,argx,argy,AmpTree,AmpProp,AmpVert,AmpBox)

      implicit none

#include "QQ2xx_Model.h"
	

      integer i,mytree,channel_,corrProp_,corrVert_,corrVert_Q,corrVert_S,corrVert_V,corrVert_T,corrVert_U,corrBox_,argx,argy
      double precision VV_QX(6),VV_SX(6),VV_VX(6),VV_TX(6),VV_UX(6),VV_TXa(6),VV_TXb(6),VV_UXa(6),VV_UXb(6),AmpTree,AmpTreeCh(25),AmpProp,AmpVert,AmpBox  

      external  QQ2xx_CorrQXV,QQ2xx_CorrTXSA, QQ2xx_CorrTXSB, QQ2xx_CorrUXSA, QQ2xx_CorrUXSB, QQ2xx_M2GluonBoxVV, QQ2xx_M2GluinoBoxVV, QQ2xx_M2TriangleBoxVV    

! calculate further argx/argy dependent O(alpha_s) corrections (effective vertices)

      call QQ2xx_CorrQXV(mytree,argx,argy)
      call QQ2xx_CorrTXVA(mytree,argx)
      call QQ2xx_CorrUXVB(mytree,argx)                  
      call QQ2xx_CorrTXVB(mytree,argy)            
      call QQ2xx_CorrUXVA(mytree,argy)

! tree-level amplitudes:

      call VV_QX_(1,argx,argy,VV_QX)                        !channel_=1,2,3,4,5,6 -> VV_QQ,VV_QS,VV_QV,VV_QT,VV_QU,VV_QX
      call VV_SX_(1,argx,argy,VV_SX)                        !channel_=1,2,3,4,5,6 -> VV_SQ,VV_SS,VV_SV,VV_ST,VV_SU,VV_SX
      call VV_VX_(1,argx,argy,VV_VX)                        !channel_=1,2,3,4,5,6 -> VV_VQ,VV_VS,VV_VV,VV_VT,VV_VU,VV_VX
      call VV_TX_(1,1,1,argx,argy,VV_TX)                    !channel_=1,2,3,4,5,6 -> VV_TQ,VV_TS,VV_TV,VV_TT,VV_TU,VV_TX
      call VV_UX_(1,1,1,argx,argy,VV_UX)                    !channel_=1,2,3,4,5,6 -> VV_UQ,VV_US,VV_UV,VV_UT,VV_UU,VV_UX


      if(channel_.eq.99) then

        AmpTree = VV_QX(6) + VV_SX(6) + VV_VX(6) + VV_TX(6) + VV_UX(6)

      else

        do i = 1, 5

          AmpTreeCh(i) = VV_QX(i)
          AmpTreeCh(i+5) = VV_SX(i)
          AmpTreeCh(i+10) = VV_VX(i)
          AmpTreeCh(i+15) = VV_TX(i)
          AmpTreeCh(i+20) = VV_UX(i)

        enddo

        AmpTree = AmpTreeCh(channel_)

      endif


      if(mytree.eq.0) then
                                     
!propagator corrections:

        if(corrProp_.eq.99) then                            !full propagator corrections

          call VV_TX_(1,1,7,argx,argy,VV_TX)
          call VV_UX_(1,1,7,argx,argy,VV_UX) 

          AmpProp = VV_TX(6) + VV_UX(6)

        else                                              

          AmpProp = 0d0

        endif                      

!vertex corrections:

        if(corrVert_.eq.99) then                            !full Vertex corrections

          corrVert_Q = 5
          corrVert_S = 6
          corrVert_V = 8
          corrVert_T = 8
          corrVert_U = 8

          call VV_QX_(corrVert_Q,argx,argy,VV_QX)         
          call VV_SX_(corrVert_S,argx,argy,VV_SX)
          call VV_VX_(corrVert_V,argx,argy,VV_VX)
          call VV_TX_(corrVert_T,1,1,argx,argy,VV_TXa)
          call VV_TX_(1,corrVert_T,1,argx,argy,VV_TXb)
          call VV_UX_(corrVert_U,1,1,argx,argy,VV_UXa)
          call VV_UX_(1,corrVert_U,1,argx,argy,VV_UXb)
             

          AmpVert = VV_QX(6) + VV_SX(6) + VV_VX(6) + VV_TXa(6) + VV_TXb(6) + VV_UXa(6) + VV_UXb(6)  


          call QQ2xx_M2GluonTriangleBoxVV(argx,argy,VV_TX,VV_UX)    !1,2,3,4,5,6->VV_TQ,VV_TS,VV_TV,VV_TT,VV_TU,VV_TX/VV_UQ,VV_US,VV_UV,VV_UT,VV_UU,VV_UX

          AmpVert = AmpVert + VV_TX(6) + VV_UX(6)

          call QQ2xx_M2GluonsQuarkLoopBoxVV(argx,argy,VV_TX,VV_UX)

          AmpVert = AmpVert + VV_TX(6) + VV_UX(6)

          if(Looptype0.eq.1) then                           !pure sQuark triangle (O(gs^2)) only added in the QCD-correction case

            call QQ2xx_M2sQuarkTriangleBoxVV(argx,argy,VV_TX,VV_UX)

            AmpVert = AmpVert + VV_TX(6) + VV_UX(6)

          endif
          
        else if(corrVert_.eq.999) then                      !Sommerfeld corrections

          corrVert_Q = 2
          corrVert_S = 2
          corrVert_V = 2

          call VV_QX_(corrVert_Q,argx,argy,VV_QX)         
          call VV_SX_(corrVert_S,argx,argy,VV_SX)
          call VV_VX_(corrVert_V,argx,argy,VV_VX)           

          AmpVert = VV_QX(6) + VV_SX(6) + VV_VX(6)   

        else

          AmpVert = 0d0   

        endif   
 
!box corrections:

        if(corrBox_.eq.99) then                             !full Box corrections

          call QQ2xx_M2GluonBoxVV(argx,argy,VV_TX,VV_UX)    !1,2,3,4,5,6->VV_TQ,VV_TS,VV_TV,VV_TT,VV_TU,VV_TX/VV_UQ,VV_US,VV_UV,VV_UT,VV_UU,VV_UX

          AmpBox = VV_TX(6) + VV_UX(6)

          call QQ2xx_M2GluinoBoxVV(argx,argy,VV_TX,VV_UX)

          AmpBox = AmpBox + VV_TX(6) + VV_UX(6)
          
        else if(corrBox_.eq.999) then                       !Sommerfeld corrections

          call QQ2xx_M2GluonBoxVV(argx,argy,VV_TX,VV_UX)    !1,2,3,4,5,6->VV_TQ,VV_TS,VV_TV,VV_TT,VV_TU,VV_TX/VV_UQ,VV_US,VV_UV,VV_UT,VV_UU,VV_UX

          AmpBox = VV_TX(6) + VV_UX(6)  

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


      subroutine VV_QX_(corr,argx,argy,VV_QX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer j,corr,argx,argy
      double precision Amp_QQ,MqqVV,Amp_QS,MqsHVV,Amp_QV,MqsVVV,Amp_QT,Amp_QU,MqtVV,MquVV,VV_QX(6)

      Amp_QS = 0d0
      Amp_QV = 0d0
      Amp_QT = 0d0
      Amp_QU = 0d0

!==================================================================================================================================

      SS = (1d0,0d0)
      Coupl_A = Coupl_AQV(corr) 

!============================================================q*q===================================================================

      Coupl_Ac = dconjg(VVSfSf(argx,argy,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))   
 
      Amp_QQ = colorfactorVV_QQ*MqqVV(Coupl_A,Coupl_Ac)

!==========================================================q*Scalar================================================================

      do j= sumstartH,sumendH

        Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
        Coupl_Bc = dconjg(LLR(argx,argy,j))
        propac = dconjg(DenSmH(j))      

 	Amp_QS = Amp_QS + MqsHVV(Coupl_A,Coupl_Ac,Coupl_Bc,propac)  

      enddo

!==========================================================q*Vector================================================================

      do j= sumstartVV,sumendV

        Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
        Coupl_Bc = dconjg(VVV(j,argx,argy))
        propac = dconjg(DenSmV(j))    

	Amp_QV = Amp_QV + colorfactorVV_QS*MqsVVV(Coupl_A,Coupl_Ac,Coupl_Bc,propac)   

      enddo

!=============================================================q*t==================================================================       

      do j= 1,2
      
        Coupl_Ac = dconjg(VSfSft(argx,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))
        Coupl_Bc = dconjg(VSfSfu(argy,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
        sfmc = MSf(j,proptypetc,propgentc)

	Amp_QT = Amp_QT + colorfactorVV_QTU*MqtVV(SS,Coupl_A,Coupl_Ac,Coupl_Bc,tman,sfmc)

      enddo 

!=============================================================q*u==================================================================        

      do j= 1,2

        Coupl_Ac = dconjg(VSfSfu(argy,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
        Coupl_Bc = dconjg(VSfSft(argx,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
        sfmc = MSf(j,proptypeuc,propgenuc)  

        Amp_QU = Amp_QU + colorfactorVV_QTU*MquVV(SS,Coupl_A,Coupl_Ac,Coupl_Bc,uman,sfmc)

      enddo 

!==================================================================================================================================

      VV_QX(1) = Amp_QQ 
      VV_QX(2) = Amp_QS
      VV_QX(3) = Amp_QV
      VV_QX(4) = Amp_QT
      VV_QX(5) = Amp_QU

      VV_QX(6) = VV_QX(1) + VV_QX(2) + VV_QX(3) + VV_QX(4) + VV_QX(5)

!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************


	
      subroutine VV_SX_(corr,argx,argy,VV_SX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corr,argx,argy
      double precision Amp_SQ,MsqHVV,Amp_SS,MssHVV,VV_SS,Amp_SV,MssHVVV,Amp_ST,Amp_SU,MstHVV,MsuHVV,VV_SX(6),MssHGhGh,MssVHGhGh,MssHVGhGh
      complex*16 Coupl_AGh,Coupl_BGh,Coupl_APGh,Coupl_AMGh,propa_Gh

      Amp_SQ = 0d0
      Amp_SS = 0d0
      Amp_SV = 0d0
      Amp_ST = 0d0
      Amp_SU = 0d0

      do i= sumstartH,sumendH                                   

!==================================================================================================================================

        SS = (1d0,0d0)
        Coupl_A = Coupl_AS(corr,i)
        Coupl_B = LLR(argx,argy,i)
        propa = DenSmH(i)

!=========================================================Scalar*q=================================================================
 
        Coupl_Ac = dconjg(VVSfSf(argx,argy,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))

	Amp_SQ = Amp_SQ + MsqHVV(Coupl_A,Coupl_Ac,Coupl_B,propa)

!=======================================================Scalar*Scalar==============================================================

        do j= sumstartH,sumendH

          Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl_Bc = dconjg(LLR(argx,argy,j))
          propac = dconjg(DenSmH(j))

	  Amp_SS = Amp_SS + MssHVV(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac) 

!Ghost-contribution (general*dconjg(amp))
          Coupl_AGh = Coupl_AS(corr,i)/2d0	
          Coupl_BGh = GGH(argx,argy,i)
          propa_Gh = DenSmH(i)
          Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl_Bc = GGH(argy,argx,j)
          propac = dconjg(DenSmH(j))        

	  Amp_SS = Amp_SS + MssHGhGh(Coupl_AGh,Coupl_Ac,Coupl_BGh,Coupl_Bc,propa_Gh,propac) 


!Ghost-contribution (amp*dconjg(general))
          Coupl_AGh = HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)
          Coupl_BGh = GGH(argx,argy,j)
          propa_Gh = DenSmH(j)
          Coupl_Ac = dconjg(Coupl_AS(corr,i))/2d0
          Coupl_Bc = GGH(argy,argx,i)
          propac = dconjg(DenSmH(i))        

	  Amp_SS = Amp_SS + MssHGhGh(Coupl_AGh,Coupl_Ac,Coupl_BGh,Coupl_Bc,propa_Gh,propac) 

        enddo  

!=====================================================Scalar*Vector================================================================

        do j=sumstartVV,sumendV 

          Coupl_Ac =  dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl_Bc = dconjg(VVV(j,argx,argy))
          propac = dconjg(DenSmV(j)) 
 
          Amp_SV = Amp_SV + MssHVVV(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)   

!Ghost-contribution (general*dconjg(Amp))
          Coupl_AGh = Coupl_AS(corr,i)/2d0	
          Coupl_BGh = GGH(argx,argy,i)
          propa_Gh = DenSmH(i)
          Coupl_APc = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl_AMc = (0d0,0d0)
          Coupl_Bc = dconjg(GGV(argx,argy,j))
          propac = dconjg(DenSmV(j))    

          Amp_SV = Amp_SV + MssHVGhGh(Coupl_AGh,Coupl_APc,Coupl_AMc,Coupl_BGh,Coupl_Bc,propa_Gh,propac)

!Ghost-contribution (amp*dconjg(general))
          Coupl_APGh = VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)
          Coupl_AMGh = (0d0,0d0)
          Coupl_BGh = GGV(argx,argy,j)
          propa_Gh = DenSmV(j)
          Coupl_Ac = dconjg(Coupl_AS(corr,i))/2d0
          Coupl_Bc = GGH(argy,argx,i)
          propac = dconjg(DenSmH(i))
 
          Amp_SV = Amp_SV + MssVHGhGh(Coupl_APGh,Coupl_AMGh,Coupl_Ac,Coupl_BGh,Coupl_Bc,propa_Gh,propac) 

        enddo

!=======================================================Scalar*t===================================================================

        do j=1,2

          Coupl_Ac = dconjg(VSfSft(argx,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))
          Coupl_Bc = dconjg(VSfSfu(argy,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
          sfmc = MSf(j,proptypetc,propgentc) 

	  Amp_ST = Amp_ST + MstHVV(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,tman,sfmc)   

        enddo

!=======================================================Scalar*u===================================================================

        do j=1,2

          Coupl_Ac = dconjg(VSfSfu(argy,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
          Coupl_Bc = dconjg(VSfSft(argx,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
          sfmc = MSf(j,proptypeuc,propgenuc) 

	  Amp_SU = Amp_SU + MsuHVV(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,uman,sfmc)  

        enddo

!==================================================================================================================================

      enddo

      VV_SX(1) = Amp_SQ
      VV_SX(2) = Amp_SS
      VV_SX(3) = Amp_SV
      VV_SX(4) = Amp_ST
      VV_SX(5) = Amp_SU

      VV_SX(6) = VV_SX(1) + VV_SX(2) + VV_SX(3) + VV_SX(4) + VV_SX(5)

!==================================================================================================================================


      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine VV_VX_(corr,argx,argy,VV_VX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corr,argx,argy
      double precision Amp_VQ,MsqVVV,Amp_VS,MssVHVV,Amp_VV,MssVVV,Amp_VT,Amp_VU,VV_VX(6),MssVHGhGh,MssHVGhGh,MssVGhGha,MssVGhGhb
      complex*16 MstVVV,MsuVVV,Coupl_AGh,Coupl_BGh,Coupl_APGh,Coupl_AMGh,propa_Gh

      Amp_VQ = 0d0
      Amp_VS = 0d0
      Amp_VV = 0d0
      Amp_VT = 0d0
      Amp_VU = 0d0

      do i=sumstartVV,sumendV 

!==================================================================================================================================

        SS = (1d0,0d0)
        Coupl_AP = Coupl_APV(corr,i)
        Coupl_AM = Coupl_AMV(corr,i)
        Coupl_B = VVV(i,argx,argy)
        propa = DenSmV(i)

!==========================================================Vector*q================================================================

        Coupl_Ac = dconjg(VVSfSf(argx,argy,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))
 
        Amp_VQ = Amp_VQ + colorfactorVV_QS*MsqVVV(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,propa)

!========================================================Vector*Scalar=============================================================

        do j=sumstartH,sumendH 

          Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl_Bc = dconjg(LLR(argx,argy,j))
          propac = dconjg(DenSmH(j))
 
          Amp_VS = Amp_VS + MssVHVV(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)     

!Ghost-contribution (general*dconjg(amp))
          Coupl_APGh = Coupl_APV(corr,i)/2d0
          Coupl_AMGh = Coupl_AMV(corr,i)/2d0
          Coupl_BGh = GGV(argx,argy,i)
          propa_Gh = DenSmV(i)
          Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl_Bc = GGH(argy,argx,j)
          propac = dconjg(DenSmH(j))
  
          Amp_VS = Amp_VS + MssVHGhGh(Coupl_APGh,Coupl_AMGh,Coupl_Ac,Coupl_BGh,Coupl_Bc,propa_Gh,propac)  

!Ghost-contribution (amp*dconjg(general)) 
          Coupl_AGh = HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)
          Coupl_BGh = GGH(argx,argy,j)
          propa_Gh = DenSmH(j)
          Coupl_APc = dconjg(Coupl_APV(corr,i))/2d0
          Coupl_AMc = dconjg(Coupl_AMV(corr,i))/2d0
          Coupl_Bc = dconjg(GGV(argx,argy,i))
          propac = dconjg(DenSmV(i))
 
          Amp_VS = Amp_VS + MssHVGhGh(Coupl_AGh,Coupl_APc,Coupl_AMc,Coupl_BGh,Coupl_Bc,propa_Gh,propac)      

        enddo

!========================================================Vector-Vector=============================================================
    
        do j=sumstartVV,sumendV 

          Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl_Bc = dconjg(VVV(j,argx,argy))
          propac = dconjg(DenSmV(j))

	  Amp_VV = Amp_VV + colorfactorVV_SS*MssVVV(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

!Ghost-contribution  (general*dconjg(amp))
          Coupl_APGh = Coupl_APV(corr,i)/2d0
          Coupl_AMGh = Coupl_AMV(corr,i)/2d0
          Coupl_BGh = GGV(argx,argy,i)
          propa_Gh = DenSmV(i)
          Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl_Bc = dconjg(GGV(argx,argy,j))
          propac = dconjg(DenSmV(j))
        
          Amp_VV = Amp_VV + colorfactorVV_SS*MssVGhGha(Coupl_APGh,Coupl_AMGh,Coupl_Ac,Coupl_BGh,Coupl_Bc,propa_Gh,propac)

!Ghost-contribution (amp*dconjg(general))
          Coupl_AGh = VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)
          Coupl_BGh = GGV(argx,argy,j)
          propa_Gh = DenSmV(j)
          Coupl_APc = dconjg(Coupl_APV(corr,i))/2d0
          Coupl_AMc = dconjg(Coupl_AMV(corr,i))/2d0
          Coupl_Bc = dconjg(GGV(argx,argy,i))
          propac = dconjg(DenSmV(i))
        
          Amp_VV = Amp_VV + colorfactorVV_SS*MssVGhGhb(Coupl_AGh,Coupl_APc,Coupl_AMc,Coupl_BGh,Coupl_Bc,propa_Gh,propac) 

        enddo

!==========================================================Vector*t================================================================

        do j=1,2  

          Coupl_Ac = dconjg(VSfSft(argx,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))
          Coupl_Bc = dconjg(VSfSfu(argy,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
          sfmc =  MSf(j,proptypetc,propgentc) 
 
          Amp_VT = Amp_VT + dreal(colorfactorVV_STU*MstVVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,tman,sfmc))  

        enddo

!========================================================Vector*u==================================================================

        do j=1,2  

          Coupl_Ac = dconjg(VSfSfu(argy,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
          Coupl_Bc = dconjg(VSfSft(argx,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
          sfmc =  MSf(j,proptypeuc,propgenuc)

          Amp_VU = Amp_VU + dreal(dconjg(colorfactorVV_STU)*MsuVVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,uman,sfmc)) 

        enddo

!==================================================================================================================================

      enddo

      VV_VX(1) = Amp_VQ
      VV_VX(2) = Amp_VS
      VV_VX(3) = Amp_VV
      VV_VX(4) = Amp_VT
      VV_VX(5) = Amp_VU

      VV_VX(6) = VV_VX(1) + VV_VX(2) + VV_VX(3) + VV_VX(4) + VV_VX(5)

!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine VV_TX_(corrA,corrB,corrProp_,argx,argy,VV_TX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corrA,corrB,corrProp_,argx,argy
      double precision Amp_TQ,MtqVV,Amp_TS,MtsHVV,Amp_TV,Amp_TT,Amp_TU,MttVV,MtuVV,VV_TX(6)
      complex*16 MtsVVV

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
        Coupl_BP = Coupl_BPT(corrB,i)
        Coupl_BM = Coupl_BMT(corrB,i)

!==========================================================t*q=====================================================================

        Coupl_Ac = dconjg(VVSfSf(argx,argy,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))

	Amp_TQ = Amp_TQ + colorfactorVV_QTU*MtqVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,tman,sfm)     

!========================================================t*Scalar==================================================================

        do j=sumstartH,sumendH 

          Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl_Bc = dconjg(LLR(argx,argy,j))
          propac = dconjg(DenSmH(j))    

	  Amp_TS = Amp_TS + MtsHVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,propac,tman,sfm)      

        enddo

!========================================================t*Vector==================================================================
    
        do j=sumstartVV,sumendV 

          Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl_Bc = dconjg(VVV(j,argx,argy))
          propac = dconjg(DenSmV(j)) 

          Amp_TV = Amp_TV + dreal(dconjg(colorfactorVV_STU)*MtsVVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,propac,tman,sfm) )  

        enddo

!========================================================t*t=======================================================================

        do j=1,2  

          Coupl_Ac = dconjg(VSfSft(argx,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))
          Coupl_Bc = dconjg(VSfSfu(argy,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
          sfmc = MSf(j,proptypetc,propgentc)  

	  Amp_TT = Amp_TT + colorfactorVV_TT*MttVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,tman,tman,sfm,sfmc) 

        enddo

!========================================================t*u=======================================================================

        do j=1,2  

          Coupl_Ac = dconjg(VSfSfu(argy,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
          Coupl_Bc = dconjg(VSfSft(argx,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
          sfmc = MSf(j,proptypeuc,propgenuc)   

	  Amp_TU = Amp_TU + colorfactorVV_TU*MtuVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,tman,uman,sfm,sfmc)  

        enddo

!==================================================================================================================================

      enddo

      VV_TX(1) = Amp_TQ
      VV_TX(2) = Amp_TS
      VV_TX(3) = Amp_TV
      VV_TX(4) = Amp_TT
      VV_TX(5) = Amp_TU

      VV_TX(6) = VV_TX(1) + VV_TX(2) + VV_TX(3) + VV_TX(4) + VV_TX(5)

!==================================================================================================================================

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine VV_UX_(corrA,corrB,corrProp_,argx,argy,VV_UX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corr,corrA,corrB,corrProp_,argx,argy
      double precision Amp_UQ,MuqVV,Amp_US,MusHVV,VV_US,Amp_UV,Amp_UT,Amp_UU,MutVV,MuuVV,VV_UX(6)
      complex*16 MusVVV

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
        Coupl_BP = Coupl_BPU(corrB,i)
        Coupl_BM = Coupl_BMU(corrB,i)

!==========================================================u*q=====================================================================

        Coupl_Ac = dconjg(VVSfSf(argx,argy,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1))  

	Amp_UQ = Amp_UQ + colorfactorVV_QTU*MuqVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,uman,sfm)

!========================================================u*Scalar==================================================================

        do j=sumstartH,sumendH 

          Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl_Bc = dconjg(LLR(argx,argy,j))
          propac = dconjg(DenSmH(j))   
        
          Amp_US = Amp_US + MusHVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,propac,uman,sfm)               

        enddo

!========================================================u*Vector==================================================================
    
        do j=sumstartVV,sumendV 

          Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl_Bc = dconjg(VVV(j,argx,argy))
          propac = dconjg(DenSmV(j))     

	  Amp_UV = Amp_UV + dreal(colorfactorVV_STU*MusVVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,propac,uman,sfm))

        enddo

!========================================================u*t=======================================================================

        do j=1,2  

          Coupl_Ac = dconjg(VSfSft(argx,j,sfmeigen1,proptypetc,sftype1,propgentc,sfgen1))
          Coupl_Bc = dconjg(VSfSfu(argy,sfmeigen2,j,sftype2,proptypetc,sfgen2,propgentc))
          sfmc = MSf(j,proptypetc,propgentc) 

	  Amp_UT = Amp_UT +  colorfactorVV_TU*MutVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,uman,tman,sfm,sfmc)

        enddo

!========================================================u*u=======================================================================

        do j=1,2  

          Coupl_Ac = dconjg(VSfSfu(argy,j,sfmeigen1,proptypeuc,sftype1,propgenuc,sfgen1))
          Coupl_Bc = dconjg(VSfSft(argx,sfmeigen2,j,sftype2,proptypeuc,sfgen2,propgenuc))
          sfmc = MSf(j,proptypeuc,propgenuc)        

	  Amp_UU = Amp_UU + colorfactorVV_TT*MuuVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,uman,uman,sfm,sfmc) 

        enddo

!==================================================================================================================================

      enddo

      VV_UX(1) = Amp_UQ
      VV_UX(2) = Amp_US
      VV_UX(3) = Amp_UV
      VV_UX(4) = Amp_UT
      VV_UX(5) = Amp_UU

      VV_UX(6) = VV_UX(1) + VV_UX(2) + VV_UX(3) + VV_UX(4) + VV_UX(5)

!==================================================================================================================================

      end


