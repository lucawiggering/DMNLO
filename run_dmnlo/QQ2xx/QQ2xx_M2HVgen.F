!1->2

	double precision function MHSV_width(s_in,m_a,m_b,Coupl_A,Coupl_Ac) 

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	double precision s_in,m_a,m_b
	complex*16 Coupl_A,Coupl_Ac


!===============================================Scalar-decay===========================================

!- from Vector-Vector-contraction

        MHSV_width = -dreal(Coupl_A*Coupl_Ac*(2*m_a**2 - m_b**2 + 2*s_in))
           

!      print*,'MHSV_width= ',MHSV_width

!======================================================================================================

      end


!2->2
!======================================================================================================

	double precision function MssHSV(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Higgs-exchange)========================================

!- from Vector-Vector-contraction

      MssHSV = dreal( -(2*m3**2 - m4**2 + 2*sqrtS**2)*Coupl_A*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac)   )   
           

!      print*,'MssHSV= ',MssHSV

!======================================================================================================

      end




	double precision function MssVSV(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Vector-exchange)========================================

! (- from Vector-Propagator)^2 * -from Vector-Contraction

      MssVSV = - dreal( (Coupl_AM*(m1**2 - m2**2) + Coupl_AP*(2*m1**2 + 2*m2**2 - sqrtS**2))*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac)  )   
           

!      print*,'MssVSV= ',MssVSV

!======================================================================================================

      end




	double precision function MssVHSV(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Vector*Higgs-exchange)=================================

! - from Vector-Propagator* -from Vector-Contraction*- from JLR-defi

      MssVHSV =  -dreal( (Coupl_AM*(m3**2 - m4**2 + 3*sqrtS**2) + Coupl_AP*(3*m1**2 - 3*m2**2 - tman + uman))/2.
     &                  *Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac)  )           
           

!      print*,'MssVHSV= ',MssVHSV

!======================================================================================================

      end




	double precision function MssHVSV(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Higgs*Vector-exchange)=================================

!- from Vector-Vector-contraction * - from Vector-propagator * - from JLR-Defi

      MssHVSV = -(3*m1**2 - 3*m2**2 - tman + uman)/2.*dreal(Coupl_A*Coupl_B*Coupl_Ac*Coupl_Bc/(propa*propac)  )       
           

!      print*,'MssHVSV= ',MssHVSV

!======================================================================================================

      end




	double precision function MtHV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,sfm,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc
	double precision sfm,sfmc



!====================================t-channel (sfermion-exchange)=====================================

!- from Vector-contraction

      MtHV = -dreal( SS*(Coupl_AM*(-m2**2 + tman) + Coupl_AP*(2*m2**2 - m4**2 + 2*tman))*Coupl_B*Coupl_Ac*Coupl_Bc/((tman-sfm**2)*(tman-sfmc**2))  )  
           

!      print*,'MtHV= ',MtHV

!======================================================================================================

      end




	double precision function MuHV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,sfm,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc
	double precision sfm,sfmc



!====================================u-channel (sfermion-exchange)=====================================

!- from Vector-contraction

      MuHV = -dreal( SS*(Coupl_AM*(m1**2 - uman) + Coupl_AP*(2*m1**2 - m4**2 + 2*uman))*Coupl_B*Coupl_Ac*Coupl_Bc/((uman-sfm**2)*(uman-sfmc**2))  )  
           

!      print*,'MuHV= ',MuHV

!======================================================================================================

      end




	double precision function MtuSV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,sfm,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc
	double precision sfm,sfmc



!====================================tu-channel (sfermion-exchange)====================================

!- from Vector-contraction

      MtuSV = -dreal( SS*(Coupl_AP*(2*m1**2 + 2*m2**2 - m3**2 - sqrtS**2) + Coupl_AM*(m1**2 - uman))*Coupl_B*Coupl_Ac*Coupl_Bc/((tman-sfm**2)*(uman-sfmc**2))  )  
           

!      print*,'MtuSV= ',MtuSV

!======================================================================================================

      end




	double precision function MutSV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,sfm,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc
	double precision sfm,sfmc



!====================================ut-channel (sfermion-exchange)====================================

!- from Vector-contraction

      MutSV = -dreal( SS*(Coupl_AP*(2*m1**2 + 2*m2**2 - m3**2 - sqrtS**2) + Coupl_AM*(-m2**2 + tman))*Coupl_B*Coupl_Ac*Coupl_Bc/((uman-sfm**2)*(tman-sfmc**2))  )  
           

!      print*,'MutSV= ',MutSV

!======================================================================================================

      end




	double precision function MstHSV(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa
	double precision sfmc



!====================================s*t-channel (Higgs-exchange)======================================

!- from Vector-Vector-contraction * - from JLR-Defi.

      MstHSV = (m1**2 - 2*m2**2 - 2*m3**2 + uman)*dreal( SS*Coupl_A*Coupl_B*Coupl_Ac*Coupl_Bc/(propa*(tman-sfmc**2))   )  
           
!      print*,'MstHSV= ',MstHSV

!======================================================================================================

      end




	double precision function MsuHSV(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa
	double precision sfmc



!====================================s*u-channel (Higgs-exchange)======================================

!- from Vector-Vector-contraction * - from JLR-defi.

      MsuHSV = (2*m1**2 - m2**2 + 2*m3**2 - tman)*dreal( SS*Coupl_A*Coupl_B*Coupl_Ac*Coupl_Bc/(propa*(uman-sfmc**2))   )   
           
!      print*,'MsuHSV= ',MsuHSV

!======================================================================================================

      end




	double precision function MtsHSV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propac,sfm)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AM,Coupl_AP,Coupl_Ac,Coupl_B,Coupl_Bc,propac
	double precision sfm



!====================================t*s-channel (Higgs-exchange)======================================

!- from Vector-Vector-contraction * - from JLR-Defi.

      MtsHSV = dreal( SS*(Coupl_AM*(-m3**2 + sqrtS**2) + Coupl_AP*(m1**2 - 2*m2**2 - 2*m3**2 + uman))*Coupl_B*Coupl_Ac*Coupl_Bc/(propac*(tman-sfm**2))   )   
           
!      print*,'MtsHSV= ',MtsHSV

!======================================================================================================

      end




	double precision function MusHSV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propac,sfm)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propac
	double precision sfm



!====================================u*s-channel (Higgs-exchange)======================================

!- from Vector-Vector-contraction * - from JLR-Defi.

      MusHSV = dreal( SS*(Coupl_AM*(-m3**2 + sqrtS**2) + Coupl_AP*(2*m1**2 - m2**2 + 2*m3**2 - tman))*Coupl_B*Coupl_Ac*Coupl_Bc/(propac*(uman-sfm**2))   )    
           
!      print*,'MusHSV= ',MusHSV

!======================================================================================================

      end




	double precision function MtsVSV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propac,sfm)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propac
	double precision sfm



!====================================t*s-channel (Vector-exchange)========================================

!- from Vector-contraction - from Vector-Propagator 

      MtsVSV = dreal( (Coupl_AM*(m1**2 - m2**2 + tman - uman) + Coupl_AP*(m1**2 + 3*m2**2 - 2*m3**2 - 2*m4**2 + 3*tman + uman))/2.
     &               *SS*Coupl_Ac*Coupl_B*Coupl_Bc/(propac*(tman-sfm**2))  )   
           

!      print*,'MtsVSV= ',MtsVSV

!======================================================================================================

      end




	double precision function MusVSV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propac,sfm)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propac
	double precision sfm



!====================================u*s-channel (Vector-exchange)========================================

!- from Vector-contraction - from Vector-Propagator 

      MusVSV = dreal( (Coupl_AM*(m1**2 - m2**2 + tman - uman) + Coupl_AP*(3*m1**2 + m2**2 - 2*m3**2 - 2*m4**2 + tman + 3*uman))/2.
     &               *SS*Coupl_Ac*Coupl_B*Coupl_Bc/(propac*(uman-sfm**2))  )   
           

!      print*,'MusVSV= ',MusVSV

!======================================================================================================

      end




	double precision function MstVSV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa
	double precision sfmc



!====================================s*t-channel (Vector-exchange)========================================

!- from Vector-contraction * - from Vector-Propagator 

      MstVSV = dreal( (Coupl_AM*(m1**2 - 3*m2**2 - 2*m3**2 + tman + uman) + Coupl_AP*(m1**2 + 3*m2**2 - 2*m3**2 - 2*m4**2 + 3*tman + uman))/2.
     &               *SS*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*(tman-sfmc**2))  )      
           

!      print*,'MstVSV= ',MstVSV

!======================================================================================================

      end



	double precision function MsuVSV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa
	double precision sfmc



!====================================s*u-channel (Vector-exchange)========================================

!- from Vector contraction * - from Vector-Propagator 

      MsuVSV = dreal( (Coupl_AM*(3*m1**2 - m2**2 + 2*m3**2 - tman - uman) + Coupl_AP*(3*m1**2 + m2**2 - 2*m3**2 - 2*m4**2 + tman + 3*uman))/2.
     &               *SS*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*(uman-sfmc**2))  )     
           

!      print*,'MsuVSV= ',MsuVSV

!======================================================================================================

      end
