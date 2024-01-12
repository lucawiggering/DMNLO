!**********************************************************************************************************************************

      subroutine QQ2xx_M2AmpFFbar(mytree,channel_,corrVert_,f3_type,f3_gen,f4_type,f4_gen,AmpTree,AmpProp,AmpVert,AmpBox)

      implicit none

#include "QQ2xx_Model.h"
	

      integer i,mytree,channel_,corrVert_,corrVert_S,corrVert_V,f3_type,f3_gen,f4_type,f4_gen
      double precision FFbar_SX(6),FFbar_VX(6),FFbar_TX(6),AmpTree,AmpTreeCh(25),AmpProp,AmpVert,AmpBox  

!full tree level

      call FFbar_SX_(1,f3_type,f3_gen,f4_type,f4_gen,FFbar_SX)
      call FFbar_VX_(1,f3_type,f3_gen,f4_type,f4_gen,FFbar_VX)
      call FFbar_TX_(f3_type,f3_gen,f4_type,f4_gen,FFbar_TX)

      if(channel_.eq.99) then

        AmpTree = FFbar_SX(6) + FFbar_VX(6) + FFbar_TX(6)

      else

        do i = 1, 5

          AmpTreeCh(i) = 0d0
          AmpTreeCh(i+5) = FFbar_SX(i)
          AmpTreeCh(i+10) = FFbar_VX(i)
          AmpTreeCh(i+15) = FFbar_TX(i)
          AmpTreeCh(i+20) = 0d0

        enddo

        AmpTree = AmpTreeCh(channel_)

      endif

      if(mytree.eq.0) then
                                     
!propagator corrections
 
        AmpProp = 0d0

!vertex corrections

        if(corrVert_.eq.99) then                                !full Vertex correction

          corrVert_S = 6
          corrVert_V = 8

          call FFbar_SX_(corrVert_S,f3_type,f3_gen,f4_type,f4_gen,FFbar_SX)
          call FFbar_VX_(corrVert_V,f3_type,f3_gen,f4_type,f4_gen,FFbar_VX)

          AmpVert =  FFbar_SX(6) + FFbar_VX(6) 
          
        else if(corrVert_.eq.999) then                          !Sommerfeld corrections  
        
          corrVert_S = 2
          corrVert_V = 2

          call FFbar_SX_(corrVert_S,f3_type,f3_gen,f4_type,f4_gen,FFbar_SX)
          call FFbar_VX_(corrVert_V,f3_type,f3_gen,f4_type,f4_gen,FFbar_VX)

          AmpVert =  FFbar_SX(6) + FFbar_VX(6) 

        else

         AmpVert = 0d0

        endif

!box corrections

        AmpBox = 0d0

      else

       AmpProp = 0d0
       AmpVert = 0d0
       AmpBox = 0d0
 
      endif
                                                                                                                              

!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************


	
      subroutine FFbar_SX_(corr,f3_type,f3_gen,f4_type,f4_gen,FFbar_SX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corr,f3_type,f3_gen,f4_type,f4_gen
      double precision Amp_SQ,Amp_SS,MssHFFbar,Amp_SV,MssHVFFbar,Amp_ST,Amp_SU,FFbar_SX(6)

      Amp_SQ = 0d0
      Amp_SS = 0d0
      Amp_SV = 0d0
      Amp_ST = 0d0
      Amp_SU = 0d0

      do i= sumstartH,sumendH                                   

!==================================================================================================================================

        Coupl_A = Coupl_AS(corr,i)
        Coupl_BP = DL(i,f4_type,f3_type,f4_gen,f3_gen)
        Coupl_BM = DR(i,f4_type,f3_type,f4_gen,f3_gen)
        propa = DenSmH(i)

!=========================================================Scalar*q=================================================================

        Amp_SQ = 0d0

!=======================================================Scalar*Scalar==============================================================

        do j= sumstartH,sumendH

          Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl_BPc = dconjg(DL(j,f4_type,f3_type,f4_gen,f3_gen))
          Coupl_BMc = dconjg(DR(j,f4_type,f3_type,f4_gen,f3_gen))
          propac = dconjg(DenSmH(j))    

	  Amp_SS = Amp_SS + NC*MssHFFbar(Coupl_A,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_BPc,Coupl_BMc,propa,propac)  

        enddo  

!=======================================================Scalar*Vector==============================================================

        do j=sumstartV,sumendV 

          Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl_BPc = dconjg(BL(j,f4_type,f3_type,f4_gen,f3_gen))
          Coupl_BMc = dconjg(BR(j,f4_type,f3_type,f4_gen,f3_gen))
          propac = dconjg(DenSmV(j))
 
          Amp_SV = Amp_SV + NC*MssHVFFbar(Coupl_A,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_BPc,Coupl_BMc,propa,propac)         

        enddo

!=========================================================Scalar*t==================================================================

        Amp_ST = 0d0  !at tree level completely covered by FFbar_TX

!=========================================================Scalar*u=================================================================

        Amp_SU = 0d0 

!==================================================================================================================================

      enddo

      FFbar_SX(1) = Amp_SQ
      FFbar_SX(2) = Amp_SS
      FFbar_SX(3) = Amp_SV
      FFbar_SX(4) = Amp_ST
      FFbar_SX(5) = Amp_SU

      FFbar_SX(6) = FFbar_SX(1) + FFbar_SX(2) + FFbar_SX(3) + FFbar_SX(4) + FFbar_SX(5)

!==================================================================================================================================


      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine FFbar_VX_(corr,f3_type,f3_gen,f4_type,f4_gen,FFbar_VX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,corr,f3_type,f3_gen,f4_type,f4_gen
      double precision Amp_VQ,Amp_VS,MssVHFFbar,Amp_VV,MssVFFbar,Amp_VT,Amp_VU,FFbar_VU,FFbar_VX(6)

      Amp_VQ = 0d0
      Amp_VS = 0d0
      Amp_VV = 0d0
      Amp_VT = 0d0
      Amp_VU = 0d0

!==========================================================Vector*q================================================================
 
      Amp_VQ = 0d0

!========================================================Vector*Scalar=============================================================

      do i=sumstartV,sumendV 
        do j=sumstartH,sumendH 

          Coupl_AP = Coupl_APV(corr,i)
          Coupl_AM = Coupl_AMV(corr,i)
          Coupl_BP = BL(i,f4_type,f3_type,f4_gen,f3_gen)
          Coupl_BM = BR(i,f4_type,f3_type,f4_gen,f3_gen)
          propa = DenSmV(i)
          Coupl_Ac = dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))
          Coupl_BPc = dconjg(DL(j,f4_type,f3_type,f4_gen,f3_gen))
          Coupl_BMc = dconjg(DR(j,f4_type,f3_type,f4_gen,f3_gen))
          propac = dconjg(DenSmH(j))
 
          Amp_VS = Amp_VS + NC*MssVHFFbar(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_BPc,Coupl_BMc,propa,propac)              

        enddo
       enddo

!========================================================Vector-Vector=============================================================
!Gluon exchange

      Coupl_AM = (0d0,0d0)
      Coupl_AP = VSfSf(1,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)
      Coupl_Ac = dconjg(VSfSf(1,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
      Coupl_BP = BL(1,f4_type,f3_type,f4_gen,f3_gen)
      Coupl_BM = BR(1,f4_type,f3_type,f4_gen,f3_gen)
      Coupl_BPc = dconjg(BL(1,f4_type,f3_type,f4_gen,f3_gen))
      Coupl_BMc = dconjg(BR(1,f4_type,f3_type,f4_gen,f3_gen))
      propa = DenSmV(1)
      propac = dconjg(DenSmV(1))

      Amp_VV = Amp_VV + 2d0/3d0*MssVFFbar(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_BPc,Coupl_BMc,propa,propac) 

!elw. Vector exchange

      do i=sumstartV,sumendV 
        do j=sumstartV,sumendV 

          Coupl_AP = Coupl_APV(corr,i)
          Coupl_AM = Coupl_AMV(corr,i)
          Coupl_BP = BL(i,f4_type,f3_type,f4_gen,f3_gen)
          Coupl_BM = BR(i,f4_type,f3_type,f4_gen,f3_gen)
          propa = DenSmV(i)
          Coupl_Ac = dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))
          Coupl_BPc = dconjg(BL(j,f4_type,f3_type,f4_gen,f3_gen))
          Coupl_BMc = dconjg(BR(j,f4_type,f3_type,f4_gen,f3_gen))
          propac = dconjg(DenSmV(j))

	  Amp_VV = Amp_VV + NC*MssVFFbar(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_BPc,Coupl_BMc,propa,propac)

        enddo
      enddo

!========================================================Vector*t==================================================================

      Amp_VT = 0d0   !at tree level completely covered by FFbar_TX

!========================================================Vector*u==================================================================

      Amp_VU = 0d0   

!==================================================================================================================================

      FFbar_VX(1) = Amp_VQ
      FFbar_VX(2) = Amp_VS
      FFbar_VX(3) = Amp_VV
      FFbar_VX(4) = Amp_VT
      FFbar_VX(5) = Amp_VU

      FFbar_VX(6) = FFbar_VX(1) + FFbar_VX(2) + FFbar_VX(3) + FFbar_VX(4) + FFbar_VX(5)

!==================================================================================================================================


      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************

      subroutine FFbar_TX_(f3_type,f3_gen,f4_type,f4_gen,FFbar_TX)

      implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_Init2to2.h"

      integer i,j,f3_type,f3_gen,f4_type,f4_gen
      double precision colorfacST,Amp_TQ,Amp_TS,Amp_TV,Amp_TT,Amp_TU,FFbar_TX(6)

      Amp_TQ = 0d0
      Amp_TS = 0d0
      Amp_TV = 0d0
      Amp_TT = 0d0
      Amp_TU = 0d0

!==========================================================t*q=====================================================================

      Amp_TQ = 0d0    

!========================================================t*Scalar==================================================================

      do j=sumstartH,sumendH 

!t(Gluino)*Higgs

	Amp_TS = Amp_TS + 2d0*CF*dreal( (    dconjg(DL(j,f4_type,f3_type,f4_gen,f3_gen))*(m3*GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*((m2**2 - m4**2 - tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) 
     &                        -2*m4*MGl*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) + GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(MGl*(m1**2 + m2**2 - tman - uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + m4*(m1**2 - m3**2 - tman)
     &                        *GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + dconjg(DR(j,f4_type,f3_type,f4_gen,f3_gen))*(-(m3*GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(2*m4*MGl*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + (-m2**2 + m4**2 + tman)
     &                        *GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) +  GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(m4*(m1**2 - m3**2 - tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + MGl*(m1**2 + m2**2 - tman - uman)
     &                        *GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) )*dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2))/( (tman-MGl2)*dconjg(DenSmH(j)) )      )       

        do i= 1, 6

!t(Neutralino/Chargino)*Higgs

	  Amp_TS =  Amp_TS + 2d0*dreal( (   
     &                        AL(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(dconjg(AL(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*(m4*(m1**2 - m3**2 - tman)*dconjg(DL(j,f4_type,f3_type,f4_gen,f3_gen)) + m3*(m2**2 - m4**2 - tman)*dconjg(DR(j,f4_type,f3_type,f4_gen,f3_gen))) + 
     &                        dconjg(AR(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*((m1**2 + m2**2 - tman - uman)*dconjg(DL(j,f4_type,f3_type,f4_gen,f3_gen)) - 2*m3*m4*dconjg(DR(j,f4_type,f3_type,f4_gen,f3_gen)))*MNeuCha(i)) + 
     &                        AR(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(dconjg(AR(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*(m3*(m2**2 - m4**2 - tman)*dconjg(DL(j,f4_type,f3_type,f4_gen,f3_gen)) + m4*(m1**2 - m3**2 - tman)*dconjg(DR(j,f4_type,f3_type,f4_gen,f3_gen))) + 
     &                        dconjg(AL(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*(-2*m3*m4*dconjg(DL(j,f4_type,f3_type,f4_gen,f3_gen)) + (m1**2 + m2**2 - tman - uman)*dconjg(DR(j,f4_type,f3_type,f4_gen,f3_gen)))*MNeuCha(i))   )
     &                       *dconjg(HSfSf(j,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)) /( (tman-MNeuCha2(i))*dconjg(DenSmH(j)) )    )

        enddo
      enddo 

!========================================================t*Vector==================================================================
    
      do j=sumstartVV,sumendV 

        if(j.eq.1) then 
          colorfacST=-2d0/9d0    !Gluino-Gluon
        else
          colorfacST=CF          !Gluino-(A,Z)
        endif

        Amp_TV = Amp_TV - 2d0*colorfacST*dreal( (   dconjg(BL(j,f4_type,f3_type,f4_gen,f3_gen))*(m3*GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(MGl*(m1**2 - m2**2 + tman - uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &                        m4*(m1**2 + m2**2 - m3**2 - m4**2 + 2*tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) + 
     &                        GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*((-(m2**2*m3**2) + m1**2*(2*m2**2 - m4**2) + m3**2*tman + m4**2*tman - 2*tman*uman )*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &                        m4*MGl*(-m1**2 + m2**2 + tman - uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + 
     &                        dconjg(BR(j,f4_type,f3_type,f4_gen,f3_gen))*(m3*GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(m4*(m1**2 + m2**2 - m3**2 - m4**2 + 2*tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &                        MGl*(m1**2 - m2**2 + tman - uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) + 
     &                        GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(m4*MGl*(-m1**2 + m2**2 + tman - uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &                        (-(m2**2*m3**2) + m1**2*(2*m2**2 - m4**2) + m3**2*tman + m4**2*tman - 2*tman*uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)))   )
     &                       *dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))   /( (tman-MGl2)*dconjg(DenSmV(j)) )   )  

        if(j.eq.1) then 
          colorfacST=CF    !Neutralino/Chargino-Gluon
        else
          colorfacST=1d0   !Neutralino/Chargino-(A,Z)/(W+-)
        endif


        do i=1, 6 

!t(Neutralino)*Vector          

	  Amp_TV = Amp_TV - 2d0*colorfacST*dreal(    ( AR(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(dconjg(AR(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &                        (m3*m4*(m1**2 + m2**2 - m3**2 - m4**2 + 2*tman)*dconjg(BR(j,f4_type,f3_type,f4_gen,f3_gen)) + 
     &                        dconjg(BL(j,f4_type,f3_type,f4_gen,f3_gen))*(-(m2**2*m3**2) + m1**2*(2*m2**2 - m4**2) + m3**2*tman + m4**2*tman - 2*tman*uman )) + 
     &                        dconjg(AL(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*(m4*(-m1**2 + m2**2 + tman - uman)*dconjg(BL(j,f4_type,f3_type,f4_gen,f3_gen)) + m3*(m1**2 - m2**2 + tman - uman)*dconjg(BR(j,f4_type,f3_type,f4_gen,f3_gen)))*MNeuCha(i)) + 
     &                        AL(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(dconjg(AL(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &                        (m3*m4*(m1**2 + m2**2 - m3**2 - m4**2 + 2*tman)*dconjg(BL(j,f4_type,f3_type,f4_gen,f3_gen)) + 
     &                        dconjg(BR(j,f4_type,f3_type,f4_gen,f3_gen))*(-(m2**2*m3**2) + m1**2*(2*m2**2 - m4**2) + m3**2*tman + m4**2*tman - 2*tman*uman)) + 
     &                        dconjg(AR(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*(m3*(m1**2 - m2**2 + tman - uman)*dconjg(BL(j,f4_type,f3_type,f4_gen,f3_gen)) + m4*(-m1**2 + m2**2 + tman - uman)*dconjg(BR(j,f4_type,f3_type,f4_gen,f3_gen)))*MNeuCha(i)))
     &                       *dconjg(VSfSf(j,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1))  /( (tman-MNeuCha2(i))*dconjg(DenSmV(j)) )   )  
                     
        enddo
      enddo      

!========================================================t*t=======================================================================  
!Gluino exchange 

      Amp_TT = Amp_TT +   (2d0/3d0)*dreal( (   -(dconjg(GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*(dconjg(GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &                    (m3*GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(MGl*(-m2**2 + m4**2 + tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 2*m4*tman*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) + 
     &                    GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*((-(m2**2*m3**2) + m3**2*m4**2 + m1**2*(m2**2 - m4**2) + m3**2*tman + m4**2*tman - tman*uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &                    m4*MGl*(-m1**2 + m3**2 + tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) + 
     &                    MGl*dconjg(GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*(m3*GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*
     &                    (2*m4*MGl*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + (-m2**2 + m4**2 + tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) + 
     &                    GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(m4*(-m1**2 + m3**2 + tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + MGl*(-m1**2 - m2**2 + tman + uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))))) + 
     &                    dconjg(GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*(MGl*dconjg(GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &                    (m3*GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*((m2**2 - m4**2 - tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) - 2*m4*MGl*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) + 
     &                    GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(MGl*(m1**2 + m2**2 - tman - uman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + m4*(m1**2 - m3**2 - tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))) - 
     &                    dconjg(GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*(m3*GluL(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(2*m4*tman*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + MGl*(-m2**2 + m4**2 + tman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) + 
     &                    GluR(sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(m4*MGl*(-m1**2 + m3**2 + tman)*GluL(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen) + 
     &                    (-(m2**2*m3**2) + m3**2*m4**2 + m1**2*(m2**2 - m4**2) + m3**2*tman + m4**2*tman - tman*uman)*GluR(sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))))   )/((tman-MGl2)*(tman-MGl2)) )


!Neutralino/Chargino-exchange

      do i= 1, 6
        do j=1, 6        
 
          Amp_TT = Amp_TT + NC*dreal( ( -(AR(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(AL(j,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)*(m3*dconjg(AL(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*
     &                   (2*m4*tman*dconjg(AR(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) + (-m2**2 + m4**2 + tman)*dconjg(AL(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*MNeuCha(i)) + 
     &                   dconjg(AR(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*(m4*(-m1**2 + m3**2 + tman)*dconjg(AR(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)) + 
     &                   (-m1**2 - m2**2 + tman + uman)*dconjg(AL(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*MNeuCha(i))*MNeuCha(j)) + 
     &                   AR(j,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)*(m4*dconjg(AL(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*MNeuCha(i)*
     &                   ((-m1**2 + m3**2 + tman)*dconjg(AR(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)) + 2*m3*dconjg(AL(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*MNeuCha(j)) + 
     &                   dconjg(AR(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*((-(m2**2*m3**2) + m3**2*m4**2 + m1**2*(m2**2 - m4**2) + m3**2*tman + m4**2*tman - tman*uman)*dconjg(AR(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)) + 
     &                   m3*(-m2**2 + m4**2 + tman)*dconjg(AL(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*MNeuCha(j))))) - 
     &                   AL(i,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)*(AR(j,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)*(m4*dconjg(AL(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*
     &                   (2*m3*tman*dconjg(AR(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)) + (-m1**2 + m3**2 + tman)*dconjg(AL(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*MNeuCha(j)) + 
     &                   dconjg(AR(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*MNeuCha(i)*(m3*(-m2**2 + m4**2 + tman)*dconjg(AR(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)) + 
     &                   (-m1**2 - m2**2 + tman + uman)*dconjg(AL(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*MNeuCha(j))) + 
     &                   AL(j,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen)*(m4*dconjg(AR(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*MNeuCha(i)*
     &                   ((-m1**2 + m3**2 + tman)*dconjg(AL(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)) + 2*m3*dconjg(AR(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*MNeuCha(j)) + 
     &                   dconjg(AL(i,sfmeigen2,sftype2,f4_type,sfgen2,f4_gen))*((-(m2**2*m3**2) + m3**2*m4**2 + m1**2*(m2**2 - m4**2) + m3**2*tman + m4**2*tman - tman*uman)*dconjg(AL(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen)) + 
     &                   m3*(-m2**2 + m4**2 + tman)*dconjg(AR(j,sfmeigen1,sftype1,f3_type,sfgen1,f3_gen))*MNeuCha(j)))) )/((tman-MNeuCha2(i))*(tman-MNeuCha2(j))) )

        enddo
      enddo

!========================================================t*u=======================================================================

      Amp_TU = 0d0 

!==================================================================================================================================

      FFbar_TX(1) = Amp_TQ
      FFbar_TX(2) = Amp_TS
      FFbar_TX(3) = Amp_TV
      FFbar_TX(4) = Amp_TT
      FFbar_TX(5) = Amp_TU

      FFbar_TX(6) = FFbar_TX(1) + FFbar_TX(2) + FFbar_TX(3) + FFbar_TX(4) + FFbar_TX(5)

!==================================================================================================================================

      end

