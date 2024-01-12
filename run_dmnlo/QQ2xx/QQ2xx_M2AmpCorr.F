!**********************************************************************************************************************************

      subroutine QQ2xx_CorrQXS(mytree,x_arg,y_arg)

      implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"
#include "QQ2xx_Init2to2.h"	

      integer mytree,x_arg,y_arg
      complex*16 QQ2xx_SSSSV,QQ2xx_SSQQQQ

!==================================================================================================================================

!Corrections to Q*X (Scalar final states)

!==================================================================================================================================

!usual tree level
      Coupl_AQS(1) = HHSfSf(x_arg,y_arg,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1)


      if(mytree.eq.0) then

!Scalar-Scalar-Scalar-Scalar-Externals with Gluon Loop
      Coupl_AQS(2) = 1d0/(8d0*pi**2)*colorfactorSSSSV*QQ2xx_SSSSV(m1**2,m2**2,sqrtS**2,(sqrtS**2+m1**2-m2**2)/2d0,x_arg,y_arg,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,Looptype0,Loopmeigen1,Looptype1,Loopgen1,Loopmeigen2,Looptype2,Loopgen2)

!*Scalar-Scalar-Scalar-Scalar-Externals with Scalar Loop
      Coupl_AQS(3) = 1d0/(8d0*pi**2)*colorfactorQQQQ*QQ2xx_SSQQQQ(sqrtS**2,x_arg,y_arg,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2) 

!Vertex Counterterm  
      Coupl_AQS(4) = 2d0*dHHSfSf(x_arg,y_arg,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1)

!full Vertex correction
      Coupl_AQS(5) = Coupl_AQS(2) + Coupl_AQS(3) + Coupl_AQS(4)

      endif

!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************




      subroutine QQ2xx_CorrQXV(mytree,x_arg,y_arg)

      implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"
#include "QQ2xx_Init2to2.h"	

      integer mytree,x_arg,y_arg
      complex*16 QQ2xx_VVQQVV,QQ2xx_VVQQQQ

!==================================================================================================================================

Corrections to Q*X (Vector final states)

!==================================================================================================================================

!usual tree level
      Coupl_AQV(1) = VVSfSf(x_arg,y_arg,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1)


      if(mytree.eq.0) then

!Scalar-Scalar-Scalar-Scalar-Externals with Gluon Loop
      Coupl_AQV(2) = 1d0/(8d0*pi**2)*colorfactorSSVVV*QQ2xx_VVQQVV(m2**2,sqrtS**2,m1**2,-(sqrtS**2-m1**2-m2**2)/2d0,Looptype0,x_arg,y_arg,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,Loopmeigen1,Looptype1,Loopgen1,Loopmeigen2,Looptype2,Loopgen2) 

!*Scalar-Scalar-Scalar-Scalar-Externals with Scalar Loop
      Coupl_AQV(3) = 1d0/(8d0*pi**2)*colorfactorQQQQ*QQ2xx_VVQQQQ(sqrtS**2,x_arg,y_arg,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2)

!Vertex Counterterm  
      Coupl_AQV(4) = 2d0*dVVSfSf(x_arg,y_arg,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1)

!full Vertex correction
      Coupl_AQV(5) = Coupl_AQV(2) + Coupl_AQV(3) + Coupl_AQV(4) 

      endif

!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************




      subroutine QQ2xx_CorrSX(mytree)

      implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"
#include "QQ2xx_Init2to2.h"	

      integer i,mytree
      complex*16 QQ2xx_SSSV,QQ2xx_SSSSV,QQ2xx_SSSGl,QQ2xx_SQQQQ

!==================================================================================================================================

!Corrections to S*X

!==================================================================================================================================

      do i=sumstartH,sumendH 

!usual tree level
      Coupl_AS(1,i) = HSfSf(i,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)


      if(mytree.eq.0) then

!Scalar-Scalar-Scalar-Externals with Gluon Loop
      Coupl_AS(2,i) = 1d0/(8d0*pi**2)*colorfactorSSSV*QQ2xx_SSSV(m1**2,m2**2,sqrtS**2,(sqrtS**2+m1**2-m2**2)/2d0,i,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,Looptype0,Loopmeigen1,Looptype1,Loopgen1,Loopmeigen2,Looptype2,Loopgen2) 

!*Scalar-Scalar-Scalar-Externals with Gluino Loop
      Coupl_AS(3,i) = 1d0/(8d0*pi**2)*colorfactorSSSGl*QQ2xx_SSSGl(m2**2,sqrtS**2,m1**2,-(sqrtS**2-m1**2-m2**2)/2d0,i,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2)

!*Scalar-Scalar-Scalar-Externals with 4-Scalar coupling
      Coupl_AS(4,i) = 1d0/(8d0*pi**2)*colorfactorQQQQ*QQ2xx_SQQQQ(sqrtS**2,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,i) 

!Vertex Counterterm  
      Coupl_AS(5,i) = 2d0*dHSfSf(i,sfmeigen1,sfmeigen2,sftype1,sftype2,sfgen1,sfgen2)

!full Vertex correction
      Coupl_AS(6,i) = Coupl_AS(2,i) + Coupl_AS(3,i) + Coupl_AS(4,i) + Coupl_AS(5,i) 

      endif


      enddo
!==================================================================================================================================

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine QQ2xx_CorrVX(mytree)

      implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"
#include "QQ2xx_Init2to2.h"	

      integer i,mytree
      complex*16 AP,AM

      external QQ2xx_SSVV, QQ2xx_SSVGl, QQ2xx_VQQQQ, QQ2xx_QQVVp1, QQ2xx_QQVVp2
!==================================================================================================================================

!Corrections to V*X

!==================================================================================================================================

      do i=sumstartVV,sumendV

!usual tree level
      Coupl_APV(1,i) = VSfSf(i,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)
      Coupl_AMV(1,i) = (0d0,0d0)


      if(mytree.eq.0) then

!Scalar-Scalar-Vector-Externals with Gluon Loop
      call QQ2xx_SSVV(m2**2,sqrtS**2,m1**2,-(sqrtS**2-m1**2-m2**2)/2d0,i,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,Looptype0,Loopmeigen1,Looptype1,Loopgen1,Loopmeigen2,Looptype2,Loopgen2,AP,AM) 

      Coupl_APV(2,i) = - 1d0/(8*pi**2)*colorfactorSSVV*AP        !- from kinematics AP(p2-p1)->(-AP)*(p1-p2)
      Coupl_AMV(2,i) =   1d0/(8*pi**2)*colorfactorSSVV*AM

!Scalar-Scalar-Vector-Externals with Gluino Loop
      call QQ2xx_SSVGl(m2**2,sqrtS**2,m1**2,-(sqrtS**2-m1**2-m2**2)/2d0,i,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,AP,AM)

      Coupl_APV(3,i) = - 1d0/(8*pi**2)*colorfactorSSVGL*AP       !- from kinematics AP(p2-p1)->(-AP)*(p1-p2)
      Coupl_AMV(3,i) =   1d0/(8*pi**2)*colorfactorSSVGl*AM

!Scalar-Scalar-Vector-Externals with 4-Scalar coupling
      call QQ2xx_VQQQQ(sqrtS**2,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,i,AP,AM) 

      Coupl_APV(4,i) = 1d0/(8*pi**2)*colorfactorQQQQ*AP
      Coupl_AMV(4,i) = 1d0/(8*pi**2)*colorfactorQQQQ*AM

!Scalar-Scalar-Vector-Externals with Scalar-Scalar-Vector-Vector-coupling
      call QQ2xx_QQVVp1(m1**2,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,i,Looptype0,Loopmeigen1,Looptype1,Loopgen1,AP,AM) 

      Coupl_APV(5,i) = 1d0/(8*pi**2)*colorfactorQQVV*AP
      Coupl_AMV(5,i) = 1d0/(8*pi**2)*colorfactorQQVV*AM

!Scalar-Scalar-Vector-Externals with Scalar-Scalar-Vector-Vector-coupling
      call QQ2xx_QQVVp2(m2**2,sfmeigen1,sftype1,sfgen1,sfmeigen2,sftype2,sfgen2,i,Looptype0,Loopmeigen2,Looptype2,Loopgen2,AP,AM) 

      Coupl_APV(6,i) =  1d0/(8*pi**2)*colorfactorQQVV*AP
      Coupl_AMV(6,i) = -1d0/(8*pi**2)*colorfactorQQVV*AM       !- from kinematics and changing flow-direction in loop relativ between p1-corr. and p2-corr.

!Vertex Counterterm  
      Coupl_APV(7,i) = 2d0*dVSfSf(i,sfmeigen2,sfmeigen1,sftype2,sftype1,sfgen2,sfgen1)  
      Coupl_AMV(7,i) = (0d0,0d0)

!full Vertex correction
      Coupl_APV(8,i) = Coupl_APV(2,i) + Coupl_APV(3,i) + Coupl_APV(4,i) + Coupl_APV(5,i) + Coupl_APV(6,i) + Coupl_APV(7,i)
      Coupl_AMV(8,i) = Coupl_AMV(2,i) + Coupl_AMV(3,i) + Coupl_AMV(4,i) + Coupl_AMV(5,i) + Coupl_AMV(6,i) + Coupl_AMV(7,i)

      endif


      enddo

!==================================================================================================================================

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************


      subroutine QQ2xx_CorrTXProp(mytree)

      implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"
#include "QQ2xx_Init2to2.h"	

      integer i,mytree
      double precision sfm2
      complex*16 QQ2xx_VS, QQ2xx_FGl, QQ2xx_Sfer, QQ2xx_V

!==================================================================================================================================

!Corrections to T*X (Propagator)

!==================================================================================================================================

      do i=1,2

!usual tree level
      SST(1,i) = (1d0,0d0)

      enddo


      if(mytree.eq.0) then

      call QQ2xx_CalcPropSquark(tman,proptypet,propgent)

      do i=1,2

      sfm2 = MSf(i,Looptypeout1,Loopgenout1)

!Propagator Gluon Loop
      SST(2,i) = -1d0/(8d0*pi**2)*colorfactorVS*QQ2xx_VS(tman,Looptype0,i,Looptypet,Loopgent,i,proptypet,propgent,i,Looptypeout1,Loopgenout1)/(tman-sfm2**2)  

!Propagator Gluino Loop
      SST(3,i) = -1d0/(8d0*pi**2)*colorfactorGl*QQ2xx_FGl(tman,Looptype0,i,proptypet,propgent,i,Looptypeout1,Loopgenout1)/(tman-sfm2**2)

!Propagator Scalar Loop
      SST(4,i) = -1d0/(8d0*pi**2)*colorfactorS*QQ2xx_Sfer(i,proptypet,propgent,i,Looptypeout1,Loopgenout1)/(tman-sfm2**2)

!Propagator Vector Loop
      SST(5,i) = -1d0/(8d0*pi**2)*colorfactorV*QQ2xx_V(Looptype0,i,proptypet,propgent,i,Looptypeout1,Loopgenout1)/(tman-sfm2**2)

!Propagator Counterterm
      SST(6,i) = -2d0*Prop(i,i,proptypet,propgent)/(tman - sfm2**2)

!full Propagator correction
      SST(7,i) = SST(2,i) + SST(3,i) + SST(4,i) + SST(5,i) + SST(6,i)

      enddo

      endif

!==================================================================================================================================

      end

!**********************************************************************************************************************************
!**********************************************************************************************************************************

      subroutine QQ2xx_CorrTXSA(mytree,x_arg)

      implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"
#include "QQ2xx_Init2to2.h"	

      integer i,mytree,x_arg
      complex*16 QQ2xx_SSSV,QQ2xx_SSSGl,QQ2xx_SQQQQ

!==================================================================================================================================

!Corrections to T*X (Vertex/Scalar/Atype)

!==================================================================================================================================

      do i=1,2

!usual tree level
      Coupl_AT(1,i) = HSfSft(x_arg,sfmeigen1,i,sftype1,proptypet,sfgen1,propgent)     


      if(mytree.eq.0) then

!*Scalar-Scalar-Scalar-Externals with Gluon Loop 
      Coupl_AT(2,i) = 1d0/(8d0*pi**2)*colorfactorSSSV*QQ2xx_SSSV(m1**2,tman,m3**2,-(tman-m1**2-m3**2)/2d0,x_arg,sfmeigen1,sftype1,sfgen1,i,proptypet,propgent,Looptype0,Loopmeigen1,Looptype1,Loopgen1,i,Looptypet,Loopgent)

!*Scalar-Scalar-Scalar-Externals with Gluino Loop 
      Coupl_AT(3,i) = 1d0/(8d0*pi**2)*colorfactorSSSGl*QQ2xx_SSSGl(tman,m3**2,m1**2,(tman+m1**2-m3**2)/2d0,x_arg,sfmeigen1,sftype1,sfgen1,i,proptypet,propgent)

!Scalar-Scalar-Scalar-Externals with 4-Scalar-Scalar coupling 
      Coupl_AT(4,i) = 1d0/(8d0*pi**2)*colorfactorQQQQ*QQ2xx_SQQQQ(m3**2,sfmeigen1,sftype1,sfgen1,i,proptypet,propgent,x_arg)

!Vertex Counterterm 
      Coupl_AT(5,i) = 2d0*dHSfSft(x_arg,sfmeigen1,i,sftype1,proptypet,sfgen1,propgent)

!full Vertex correction 
      Coupl_AT(6,i) = Coupl_AT(2,i) + Coupl_AT(3,i) + Coupl_AT(4,i) + Coupl_AT(5,i) 

      endif


      enddo

!==================================================================================================================================
  
      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************

      subroutine QQ2xx_CorrTXSB(mytree,y_arg)

      implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"
#include "QQ2xx_Init2to2.h"	

      integer i,mytree,y_arg
      complex*16 QQ2xx_SSSV,QQ2xx_SSSGl,QQ2xx_SQQQQ

!==================================================================================================================================

!Corrections to T*X (Vertex/Scalar/B-type)

!==================================================================================================================================

      do i=1,2

!usual tree level
      Coupl_BT(1,i) = HSfSfu(y_arg,i,sfmeigen2,proptypet,sftype2,propgent,sfgen2)


      if(mytree.eq.0) then

!*Scalar-Scalar-Scalar-Externals with Gluon Loop 
      Coupl_BT(2,i) = 1d0/(8d0*pi**2)*colorfactorSSSV*QQ2xx_SSSV(tman,m2**2,m4**2,(tman-m2**2+m4**2)/2d0,y_arg,i,proptypet,propgent,sfmeigen2,sftype2,sfgen2,Looptype0,i,Looptypet,Loopgent,Loopmeigen2,Looptype2,Loopgen2)

!*Scalar-Scalar-Scalar-Externals with Gluino Loop 
      Coupl_BT(3,i) = 1d0/(8d0*pi**2)*colorfactorSSSGl*QQ2xx_SSSGl(m2**2,m4**2,tman,(tman+m2**2-m4**2)/2d0,y_arg,i,proptypet,propgent,sfmeigen2,sftype2,sfgen2)

!Scalar-Scalar-Scalar-Externals with 4-Scalar-Scalar coupling 
      Coupl_BT(4,i) = 1d0/(8d0*pi**2)*colorfactorQQQQ*QQ2xx_SQQQQ(m4**2,i,proptypet,propgent,sfmeigen2,sftype2,sfgen2,y_arg)

!Vertex Counterterm 
      Coupl_BT(5,i) = 2d0*dHSfSfu(y_arg,i,sfmeigen2,proptypet,sftype2,propgent,sfgen2)

!full Vertex correction 
      Coupl_BT(6,i) = Coupl_BT(2,i) + Coupl_BT(3,i) + Coupl_BT(4,i) + Coupl_BT(5,i)

      endif
 

      enddo

!==================================================================================================================================

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************

      subroutine QQ2xx_CorrTXVA(mytree,x_arg)

      implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"
#include "QQ2xx_Init2to2.h"	

      integer i,mytree,x_arg
      complex*16 AP,AM

      external QQ2xx_SSVV, QQ2xx_SSVGl, QQ2xx_VQQQQ, QQ2xx_QQVVp1, QQ2xx_QQVVp2
!==================================================================================================================================

!Corrections to T*X (Vertex/Vector/A-type)

!==================================================================================================================================


      do i=1,2

!usual tree level
      Coupl_APT(1,i) = VSfSft(x_arg,i,sfmeigen1,proptypet,sftype1,propgent,sfgen1)
      Coupl_AMT(1,i) = (0d0,0d0)


      if(mytree.eq.0) then

!Scalar-Scalar-Vector-Externals with Gluon Loop 
      call QQ2xx_SSVV(tman,m3**2,m1**2,(tman+m1**2-m3**2)/2d0,x_arg,sfmeigen1,sftype1,sfgen1,i,proptypet,propgent,Looptype0,Loopmeigen1,Looptype1,Loopgen1,i,Looptypet,Loopgent,AP,AM)

      Coupl_APT(2,i) = - 1d0/(8d0*pi**2)*colorfactorSSVV*AP
      Coupl_AMT(2,i) =   1d0/(8d0*pi**2)*colorfactorSSVV*AM

!Scalar-Scalar-Vector-Externals with Gluino Loop 
      call QQ2xx_SSVGl(tman,m3**2,m1**2,(tman+m1**2-m3**2)/2d0,x_arg,sfmeigen1,sftype1,sfgen1,i,proptypet,propgent,AP,AM)

      Coupl_APT(3,i) = - 1d0/(8d0*pi**2)*colorfactorSSVGl*AP
      Coupl_AMT(3,i) =   1d0/(8d0*pi**2)*colorfactorSSVGl*AM

!Scalar-Scalar-Vector-Externals with 4-Scalar-Scalar coupling 
      call QQ2xx_VQQQQ(m3**2,sfmeigen1,sftype1,sfgen1,i,proptypet,propgent,x_arg,AP,AM) 

      Coupl_APT(4,i) = 1d0/(8d0*pi**2)*colorfactorQQQQ*AP
      Coupl_AMT(4,i) = 1d0/(8d0*pi**2)*colorfactorQQQQ*AM

!Scalar-Scalar-Vector-Externals with Scalar-Scalar-Vector-Vector coupling 
      call QQ2xx_QQVVp1(m1**2,sfmeigen1,sftype1,sfgen1,i,proptypet,propgent,x_arg,Looptype0,Loopmeigen1,Looptype1,Loopgen1,AP,AM)

      Coupl_APT(5,i) =   1d0/(8d0*pi**2)*colorfactorQQVV*AP
      Coupl_AMT(5,i) = - 1d0/(8d0*pi**2)*colorfactorQQVV*AM

!Scalar-Scalar-Vector-Externals with Scalar-Scalar-Vector-Vector coupling 
      call QQ2xx_QQVVp2(tman,sfmeigen1,sftype1,sfgen1,i,proptypet,propgent,x_arg,Looptype0,i,Looptypet,Loopgent,AP,AM)

      Coupl_APT(6,i) = 1d0/(8d0*pi**2)*colorfactorQQVV*AP
      Coupl_AMT(6,i) = 1d0/(8d0*pi**2)*colorfactorQQVV*AM

!Vertex Counterterm 

      Coupl_APT(7,i) = 2d0*dVSfSft(x_arg,i,sfmeigen1,proptypet,sftype1,propgent,sfgen1)
      Coupl_AMT(7,i) = (0d0,0d0)

!full Vertex correction 
      Coupl_APT(8,i) = Coupl_APT(2,i) + Coupl_APT(3,i) + Coupl_APT(4,i) + Coupl_APT(5,i) + Coupl_APT(6,i) + Coupl_APT(7,i)
      Coupl_AMT(8,i) = Coupl_AMT(2,i) + Coupl_AMT(3,i) + Coupl_AMT(4,i) + Coupl_AMT(5,i) + Coupl_AMT(6,i) + Coupl_AMT(7,i) 

      endif


      enddo

!==================================================================================================================================

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************

      subroutine QQ2xx_CorrTXVB(mytree,y_arg)

      implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"
#include "QQ2xx_Init2to2.h"	

      integer i,mytree,y_arg
      complex*16 AP,AM

      external QQ2xx_SSVV, QQ2xx_SSVGl, QQ2xx_QQVVp1, QQ2xx_QQVVp2, QQ2xx_VQQQQ
!==================================================================================================================================

!Corrections to T*X (Vertex/Vector/B-type)

!==================================================================================================================================


      do i=1,2

!usual tree level
      Coupl_BPT(1,i) = VSfSfu(y_arg,sfmeigen2,i,sftype2,proptypet,sfgen2,propgent)
      Coupl_BMT(1,i) = (0d0,0d0)


      if(mytree.eq.0) then

!Scalar-Scalar-Vector-Externals with Gluon Loop 
      call QQ2xx_SSVV(m2**2,m4**2,tman,(tman+m2**2-m4**2)/2d0,y_arg,i,proptypet,propgent,sfmeigen2,sftype2,sfgen2,Looptype0,i,Looptypet,Loopgent,Loopmeigen2,Looptype2,Loopgen2,AP,AM)

      Coupl_BPT(2,i) = - 1d0/(8d0*pi**2)*colorfactorSSVV*AP
      Coupl_BMT(2,i) =   1d0/(8d0*pi**2)*colorfactorSSVV*AM

!Scalar-Scalar-Vector-Externals with Gluino Loop 
      call QQ2xx_SSVGl(m2**2,m4**2,tman,(tman+m2**2-m4**2)/2d0,y_arg,i,proptypet,propgent,sfmeigen2,sftype2,sfgen2,AP,AM)

      Coupl_BPT(3,i) = - 1d0/(8d0*pi**2)*colorfactorSSVGl*AP
      Coupl_BMT(3,i) = 1d0/(8d0*pi**2)*colorfactorSSVGl*AM

!Scalar-Scalar-Vector-Externals with 4-Scalar-Scalar coupling 
      call QQ2xx_VQQQQ(m4**2,i,proptypet,propgent,sfmeigen2,sftype2,sfgen2,y_arg,AP,AM)

      Coupl_BPT(4,i) = 1d0/(8d0*pi**2)*colorfactorQQQQ*AP
      Coupl_BMT(4,i) = 1d0/(8d0*pi**2)*colorfactorQQQQ*AM 

!Scalar-Scalar-Vector-Externals with Scalar-Scalar-Vector-Vector coupling 
      call QQ2xx_QQVVp1(tman,i,proptypet,propgent,sfmeigen2,sftype2,sfgen2,y_arg,Looptype0,i,Looptypet,Loopgent,AP,AM)

      Coupl_BPT(5,i) =   1d0/(8d0*pi**2)*colorfactorQQVV*AP
      Coupl_BMT(5,i) = - 1d0/(8d0*pi**2)*colorfactorQQVV*AM

!Scalar-Scalar-Vector-Externals with Scalar-Scalar-Vector-Vector coupling 
      call QQ2xx_QQVVp2(m2**2,i,proptypet,propgent,sfmeigen2,sftype2,sfgen2,y_arg,Looptype0,Loopmeigen2,Looptype2,Loopgen2,AP,AM) 

      Coupl_BPT(6,i) = 1d0/(8d0*pi**2)*colorfactorQQVV*AP
      Coupl_BMT(6,i) = 1d0/(8d0*pi**2)*colorfactorQQVV*AM

!Vertex Counterterm 
      Coupl_BPT(7,i) = 2d0*dVSfSfu(y_arg,sfmeigen2,i,sftype2,proptypet,sfgen2,propgent)
      Coupl_BMT(7,i) = (0d0,0d0)

!full Vertex correction 
      Coupl_BPT(8,i) = Coupl_BPT(2,i) + Coupl_BPT(3,i) + Coupl_BPT(4,i) + Coupl_BPT(5,i) + Coupl_BPT(6,i) + Coupl_BPT(7,i)
      Coupl_BMT(8,i) = Coupl_BMT(2,i) + Coupl_BMT(3,i) + Coupl_BMT(4,i) + Coupl_BMT(5,i) + Coupl_BMT(6,i) + Coupl_BMT(7,i)   

      endif


      enddo

!==================================================================================================================================

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************

      subroutine QQ2xx_CorrUXProp(mytree)

      implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"
#include "QQ2xx_Init2to2.h"	

      integer i,mytree
      double precision sfm2
      complex*16 QQ2xx_VS, QQ2xx_FGl, QQ2xx_Sfer, QQ2xx_V

!==================================================================================================================================

!Corrections to U*X (Propagator)

!==================================================================================================================================

      do i=1,2

!usual tree level
      SSU(1,i) = (1d0,0d0)

      enddo


      if(mytree.eq.0) then

      call QQ2xx_CalcPropSquark(uman,proptypeu,propgenu)

      do i=1,2

      sfm2 = MSf(i,Looptypeout2,Loopgenout2)      

!Propagator Gluon Loop
      SSU(2,i) = -1d0/(8d0*pi**2)*colorfactorVS*QQ2xx_VS(uman,Looptype0,i,Looptypeu,Loopgenu,i,proptypeu,propgenu,i,Looptypeout2,Loopgenout2)/(uman-sfm2**2)  

!Propagator Gluino Loop
      SSU(3,i) = -1d0/(8d0*pi**2)*colorfactorGl*QQ2xx_FGl(uman,Looptype0,i,proptypeu,propgenu,i,Looptypeout2,Loopgenout2)/(uman-sfm2**2)

!Propagator Scalar Loop
      SSU(4,i) = -1d0/(8d0*pi**2)*colorfactorS*QQ2xx_Sfer(i,proptypeu,propgenu,i,Looptypeout2,Loopgenout2)/(uman-sfm2**2)

!Propagator Vector Loop
      SSU(5,i) = -1d0/(8d0*pi**2)*colorfactorV*QQ2xx_V(Looptype0,i,proptypeu,propgenu,i,Looptypeout2,Loopgenout2)/(uman-sfm2**2)

!Propagator Counterterm
      SSU(6,i) = -2d0*Prop(i,i,proptypeu,propgenu)/(uman - sfm2**2)

!full Propagator correction
      SSU(7,i) = SSU(2,i) + SSU(3,i) + SSU(4,i) + SSU(5,i) + SSU(6,i)
 
      enddo

      endif

!==================================================================================================================================

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************


      subroutine QQ2xx_CorrUXSA(mytree,y_arg)

      implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"
#include "QQ2xx_Init2to2.h"	

      integer i,mytree,y_arg
      complex*16 QQ2xx_SSSV,QQ2xx_SSSGl,QQ2xx_SQQQQ

!==================================================================================================================================

!Corrections to U*X (Vertex/Scalar/A-type)

!==================================================================================================================================

      do i=1,2

!usual tree level
      Coupl_AU(1,i) = HSfSfu(y_arg,sfmeigen1,i,sftype1,proptypeu,sfgen1,propgenu)  


      if(mytree.eq.0) then   

!Scalar-Scalar-Scalar-Externals with Gluon Loop
      Coupl_AU(2,i) = 1d0/(8d0*pi**2)*colorfactorSSSV*QQ2xx_SSSV(m1**2,uman,m4**2,-(uman-m1**2-m4**2)/2d0,y_arg,sfmeigen1,sftype1,sfgen1,i,proptypeu,propgenu,Looptype0,Loopmeigen1,Looptype1,Loopgen1,i,Looptypeu,Loopgenu)

!Scalar-Scalar-Scalar-Externals with Gluino Loop
      Coupl_AU(3,i) = 1d0/(8d0*pi**2)*colorfactorSSSGl*QQ2xx_SSSGl(uman,m4**2,m1**2,(uman+m1**2-m4**2)/2d0,y_arg,sfmeigen1,sftype1,sfgen1,i,proptypeu,propgenu)

!Scalar-Scalar-Scalar-Externals with 4-Scalar coupling
      Coupl_AU(4,i) = 1d0/(8d0*pi**2)*colorfactorQQQQ*QQ2xx_SQQQQ(m4**2,sfmeigen1,sftype1,sfgen1,i,proptypeu,propgenu,y_arg)

!Vertex Counterterm 
      Coupl_AU(5,i) = 2d0*dHSfSfu(y_arg,sfmeigen1,i,sftype1,proptypeu,sfgen1,propgenu)

!full Vertex correction 
      Coupl_AU(6,i) = Coupl_AU(2,i) + Coupl_AU(3,i) + Coupl_AU(4,i) + Coupl_AU(5,i) 

      endif


      enddo

!==================================================================================================================================

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine QQ2xx_CorrUXSB(mytree,x_arg)

      implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"
#include "QQ2xx_Init2to2.h"	

      integer i,mytree,x_arg
      complex*16 QQ2xx_SSSV,QQ2xx_SSSGl,QQ2xx_SQQQQ

!==================================================================================================================================

!Corrections to U*X (Vertex/Scalar/B-type)

!==================================================================================================================================

      do i=1,2

!usual tree level 
      Coupl_BU(1,i) = HSfSft(x_arg,i,sfmeigen2,proptypeu,sftype2,propgenu,sfgen2)


      if(mytree.eq.0) then

!Scalar-Scalar-Scalar-Externals with Gluon Loop    
      Coupl_BU(2,i) = 1d0/(8d0*pi**2)*colorfactorSSSV*QQ2xx_SSSV(uman,m2**2,m3**2,(uman-m2**2+m3**2)/2d0,x_arg,i,proptypeu,propgenu,sfmeigen2,sftype2,sfgen2,Looptype0,i,Looptypeu,Loopgenu,Loopmeigen2,Looptype2,Loopgen2)

!*Scalar-Scalar-Scalar-Externals with Gluino Loop  
      Coupl_BU(3,i) = 1d0/(8d0*pi**2)*colorfactorSSSGl*QQ2xx_SSSGl(m2**2,m3**2,uman,(uman+m2**2-m3**2)/2d0,x_arg,i,proptypeu,propgenu,sfmeigen2,sftype2,sfgen2)

!Scalar-Scalar-Scalar-Externals with 4-Scalar coupling 
      Coupl_BU(4,i) = 1d0/(8d0*pi**2)*colorfactorQQQQ*QQ2xx_SQQQQ(m3**2,i,proptypeu,propgenu,sfmeigen2,sftype2,sfgen2,x_arg)

!Vertex Counterterm 
      Coupl_BU(5,i) = 2d0*dHSfSft(x_arg,i,sfmeigen2,proptypeu,sftype2,propgenu,sfgen2)

!full Vertex correction 
      Coupl_BU(6,i) = Coupl_BU(2,i) + Coupl_BU(3,i) + Coupl_BU(4,i) + Coupl_BU(5,i) 

      endif


      enddo

!==================================================================================================================================

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************


      subroutine QQ2xx_CorrUXVA(mytree,y_arg)

      implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"
#include "QQ2xx_Init2to2.h"	

      integer i,mytree,y_arg
      complex*16 AP,AM

      external QQ2xx_SSVV, QQ2xx_SSVGl, QQ2xx_QQVVp1, QQ2xx_QQVVp2, QQ2xx_VQQQQ
!==================================================================================================================================

!Corrections to U*X (Vertex/Vector/Atype)

!==================================================================================================================================

      do i=1,2

!usual tree level
      Coupl_APU(1,i) = VSfSfu(y_arg,i,sfmeigen1,proptypeu,sftype1,propgenu,sfgen1)
      Coupl_AMU(1,i) = (0d0,0d0)


      if(mytree.eq.0) then

!Scalar-Scalar-Vector-Externals with Gluon Loop
      call QQ2xx_SSVV(uman,m4**2,m1**2,(uman+m1**2-m4**2)/2d0,y_arg,sfmeigen1,sftype1,sfgen1,i,proptypeu,propgenu,Looptype0,Loopmeigen1,Looptype1,Loopgen1,i,Looptypeu,Loopgenu,AP,AM) 

      Coupl_APU(2,i) = - 1d0/(8d0*pi**2)*colorfactorSSVV*AP
      Coupl_AMU(2,i) =   1d0/(8d0*pi**2)*colorfactorSSVV*AM

!Scalar-Scalar-Vector-Externals with Gluino Loop
      call QQ2xx_SSVGl(uman,m4**2,m1**2,(uman+m1**2-m4**2)/2d0,y_arg,sfmeigen1,sftype1,sfgen1,i,proptypeu,propgenu,AP,AM)

      Coupl_APU(3,i) = - 1d0/(8d0*pi**2)*colorfactorSSVGl*AP
      Coupl_AMU(3,i) =   1d0/(8d0*pi**2)*colorfactorSSVGl*AM

!Scalar-Scalar-Vector-Externals with 4-Scalar coupling
      call QQ2xx_VQQQQ(m4**2,sfmeigen1,sftype1,sfgen1,i,proptypeu,propgenu,y_arg,AP,AM)

      Coupl_APU(4,i) = 1d0/(8d0*pi**2)*colorfactorQQQQ*AP
      Coupl_AMU(4,i) = 1d0/(8d0*pi**2)*colorfactorQQQQ*AM

!Scalar-Scalar-Vector-Externals with Scalar-Scalar-Vector-Vector coupling
      call QQ2xx_QQVVp1(m1**2,sfmeigen1,sftype1,sfgen1,i,proptypeu,propgenu,y_arg,Looptype0,Loopmeigen1,Looptype1,Loopgen1,AP,AM) 

      Coupl_APU(5,i) =   1d0/(8d0*pi**2)*colorfactorQQVV*AP
      Coupl_AMU(5,i) = - 1d0/(8d0*pi**2)*colorfactorQQVV*AM

!Scalar-Scalar-Vector-Externals with Scalar-Scalar-Vector-Vector coupling
      call QQ2xx_QQVVp2(uman,sfmeigen1,sftype1,sfgen1,i,proptypeu,propgenu,y_arg,Looptype0,i,Looptypeu,Loopgenu,AP,AM)

      Coupl_APU(6,i) = 1d0/(8d0*pi**2)*colorfactorQQVV*AP
      Coupl_AMU(6,i) = 1d0/(8d0*pi**2)*colorfactorQQVV*AM

!Vertex Counterterm 
      Coupl_APU(7,i) = 2d0*dVSfSfu(y_arg,i,sfmeigen1,proptypeu,sftype1,propgenu,sfgen1) 
      Coupl_AMU(7,i) = (0d0,0d0)

!full Vertex correction 
      Coupl_APU(8,i) = Coupl_APU(2,i) + Coupl_APU(3,i) + Coupl_APU(4,i) + Coupl_APU(5,i) + Coupl_APU(6,i) + Coupl_APU(7,i)
      Coupl_AMU(8,i) = Coupl_AMU(2,i) + Coupl_AMU(3,i) + Coupl_AMU(4,i) + Coupl_AMU(5,i) + Coupl_AMU(6,i) + Coupl_AMU(7,i)

      endif

      enddo  

!==================================================================================================================================

      end   


!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine QQ2xx_CorrUXVB(mytree,x_arg)

      implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "QQ2xx_GenCouplings.h"
#include "QQ2xx_GenCounterterms.h"
#include "QQ2xx_Init2to2.h"	

      integer i,mytree,x_arg
      complex*16 AP,AM

      external QQ2xx_SSVV, QQ2xx_SSVGl, QQ2xx_VQQQQ, QQ2xx_QQVVp1, QQ2xx_QQVVp2
!==================================================================================================================================

!Corrections to U*X (Vertex/Vector/Btype)

!==================================================================================================================================

      do i=1,2

!usual tree level
      Coupl_BPU(1,i) = VSfSft(x_arg,sfmeigen2,i,sftype2,proptypeu,sfgen2,propgenu)
      Coupl_BMU(1,i) = (0d0,0d0)    


      if(mytree.eq.0) then

!Scalar-Scalar-Vector-Externals with Gluon Loop
      call QQ2xx_SSVV(m2**2,m3**2,uman,(uman+m2**2-m3**2)/2d0,x_arg,i,proptypeu,propgenu,sfmeigen2,sftype2,sfgen2,Looptype0,i,Looptypeu,Loopgenu,Loopmeigen2,Looptype2,Loopgen2,AP,AM)

      Coupl_BPU(2,i) = - 1d0/(8d0*pi**2)*colorfactorSSVV*AP
      Coupl_BMU(2,i) =   1d0/(8d0*pi**2)*colorfactorSSVV*AM

!Scalar-Scalar-Vector-Externals with Gluino Loop
      call QQ2xx_SSVGl(m2**2,m3**2,uman,(uman+m2**2-m3**2)/2d0,x_arg,i,proptypeu,propgenu,sfmeigen2,sftype2,sfgen2,AP,AM) 

      Coupl_BPU(3,i) = - 1d0/(8d0*pi**2)*colorfactorSSVGl*AP
      Coupl_BMU(3,i) =   1d0/(8d0*pi**2)*colorfactorSSVGl*AM

!Scalar-Scalar-Vector-Externals with 4-Scalar coupling
      call QQ2xx_VQQQQ(m3**2,i,proptypeu,propgenu,sfmeigen2,sftype2,sfgen2,x_arg,AP,AM)

      Coupl_BPU(4,i) = 1d0/(8d0*pi**2)*colorfactorQQQQ*AP
      Coupl_BMU(4,i) = 1d0/(8d0*pi**2)*colorfactorQQQQ*AM

!Scalar-Scalar-Vector-Externals with Scalar-Scalar-Vector-Vector coupling
      call QQ2xx_QQVVp1(uman,i,proptypeu,propgenu,sfmeigen2,sftype2,sfgen2,x_arg,Looptype0,i,Looptypeu,Loopgenu,AP,AM)

      Coupl_BPU(5,i) =   1d0/(8d0*pi**2)*colorfactorQQVV*AP
      Coupl_BMU(5,i) = - 1d0/(8d0*pi**2)*colorfactorQQVV*AM

!Scalar-Scalar-Vector-Externals with Scalar-Scalar-Vector-Vector coupling
      call QQ2xx_QQVVp2(m2**2,i,proptypeu,propgenu,sfmeigen2,sftype2,sfgen2,x_arg,Looptype0,Loopmeigen2,Looptype2,Loopgen2,AP,AM)

      Coupl_BPU(6,i) =   1d0/(8d0*pi**2)*colorfactorQQVV*AP
      Coupl_BMU(6,i) =   1d0/(8d0*pi**2)*colorfactorQQVV*AM

!Vertex Counterterm 
      Coupl_BPU(7,i) = 2d0*dVSfSft(x_arg,sfmeigen2,i,sftype2,proptypeu,sfgen2,propgenu)
      Coupl_BMU(7,i) = (0d0,0d0)  

!full Vertex correction 
      Coupl_BPU(8,i) = Coupl_BPU(2,i) + Coupl_BPU(3,i) + Coupl_BPU(4,i) + Coupl_BPU(5,i) + Coupl_BPU(6,i) + Coupl_BPU(7,i)
      Coupl_BMU(8,i) = Coupl_BMU(2,i) + Coupl_BMU(3,i) + Coupl_BMU(4,i) + Coupl_BMU(5,i) + Coupl_BMU(6,i) + Coupl_BMU(7,i) 

      endif


      enddo  

!==================================================================================================================================

      end





!A-type coupling couples to incoming sfermion (momentum p1)
!B-type coupling couples to incoming anti-sfermion (momentum p2)



