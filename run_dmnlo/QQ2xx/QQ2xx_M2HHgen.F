!1->2

	double precision function MHSS_width(s_in,m_a,m_b,Coupl_A,Coupl_Ac) 

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	double precision s_in,m_a,m_b
	complex*16 Coupl_A,Coupl_Ac


!===============================================Scalar-decay===========================================


        MHSS_width = dreal(Coupl_A*Coupl_Ac )
           

!      print*,'MHSS_width= ',MHSS_width

!======================================================================================================

      end


!2->2
!======================================================================================================

	double precision function MqqSS(Coupl_A,Coupl_Ac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac



!===============================================qq=====================================================


      MqqSS = dreal(Coupl_A*Coupl_Ac )
           

!      print*,'MqqSS= ',MqqSS

!======================================================================================================

      end




	double precision function MssHSS(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Higgs-exchange)========================================


      MssHSS = dreal( Coupl_A*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac) )
           

!      print*,'MssHSS= ',MssHSS

!======================================================================================================

      end




	double precision function MssVSS(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Vector-exchange)========================================


      MssVSS = dreal( (Coupl_AM*(m2**2-m1**2)+Coupl_AP*(tman-uman))*(tman-uman)*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac)  )
           

!      print*,'MssVSS= ',MssVSS

!======================================================================================================

      end




	double precision function MssVHSS(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Vector-Higgs-exchange)========================================

! - from Vector-Propagator

      MssVHSS = dreal( (Coupl_AM*(m2**2-m1**2)+Coupl_AP*(tman-uman))*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac)  )  
           

!      print*,'MssVHSS= ',MssVHSS

!======================================================================================================

      end




	double precision function MssHVSS(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Higgs-Vector-exchange)========================================

! - from Vector-Propagator

      MssHVSS = (tman-uman)*dreal( Coupl_A*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac)  )  
           

!      print*,'MssHVSS= ',MssHVSS

!======================================================================================================

      end




	double precision function MtuSS(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,man,manc,sfm,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc
	double precision man,manc,sfm,sfmc



!=============================================t/u-channel==============================================


      MtuSS =  dreal( SS*Coupl_A*Coupl_B*Coupl_Ac*Coupl_Bc/((man-sfm**2)*(manc-sfmc**2))  )
           

!      print*,'MtuSS= ',MtuSS

!======================================================================================================

      end




	double precision function MsqHSS(Coupl_A,Coupl_Ac,Coupl_B,propa)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_B,propa



!====================================s*q-channel (Higgs-exchange)========================================

!- from i-combination

      MsqHSS = -dreal( Coupl_A*Coupl_Ac*Coupl_B/propa )
           
!      print*,'MsqHSS= ',MsqHSS

!======================================================================================================

      end




	double precision function MsqVSS(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,propa)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,propa



!====================================s*q-channel (Vector-exchange)========================================

! - from Vector-Propagator * -1 from i-combination

      MsqVSS = -dreal( (Coupl_AM*(m2**2-m1**2)+Coupl_AP*(tman-uman))*Coupl_Ac*Coupl_B/propa  )  
           

!      print*,'MsqVSS= ',MsqVSS

!======================================================================================================

      end




	double precision function MqsHSS(Coupl_A,Coupl_Ac,Coupl_Bc,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_Bc,propac



!====================================q*s-channel (Higgs-exchange)========================================

!- from i-combination

      MqsHSS = -dreal( Coupl_A*Coupl_Ac*Coupl_Bc/propac )
           
!      print*,'MqsHSS= ',MqsHSS

!======================================================================================================

      end




	double precision function MqsVSS(Coupl_A,Coupl_Ac,Coupl_Bc,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_Bc,propac



!====================================q*s-channel (Vector-exchange)========================================

!- from i-combination * -1 from propagator

      MqsVSS = -(tman-uman)*dreal( Coupl_A*Coupl_Ac*Coupl_Bc/propac )
           
!      print*,'MqsVSS= ',MqsVSS

!======================================================================================================

      end




	double precision function MtuqSS(SS,Coupl_A,Coupl_Ac,Coupl_B,man,sfm)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_A,Coupl_Ac,Coupl_B
	double precision man,sfm



!=============================================t/u*q-channel============================================

!- from i-combination

      MtuqSS = - dreal( SS*Coupl_A*Coupl_Ac*Coupl_B/(man-sfm**2 ) )
           
!      print*,'MtuqSS= ',MtuqSS

!======================================================================================================

      end




	double precision function MstuHSS(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,manc,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa
	double precision manc,sfmc



!====================================s*t/u-channel (Higgs-exchange)========================================


      MstuHSS =  dreal( SS*Coupl_A*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*(manc-sfmc**2 )) )
           
!      print*,'MstuHSS= ',MstuHSS

!======================================================================================================

      end




	double precision function MstuVSS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,manc,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa
	double precision manc,sfmc



!====================================s*t/u-channel (Vector-exchange)========================================

! - from Vector-Propagator 

      MstuVSS = dreal(  (Coupl_AM*(m2**2-m1**2)+Coupl_AP*(tman-uman))*SS*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*(manc-sfmc**2))   )
           

!      print*,'MstuVSS= ',MstuVSS

!======================================================================================================

      end




	double precision function MtusVSS(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propac,man,sfm)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propac
	double precision man,sfm



!====================================t/u*s-channel (Vector-exchange)========================================

! - from Vector-Propagator 

      MtusVSS = dreal(  Coupl_A*(tman-uman)*SS*Coupl_Ac*Coupl_B*Coupl_Bc/(propac*(man-sfm**2)) )
           

!      print*,'MtusVSS= ',MtusVSS

!======================================================================================================

      end
