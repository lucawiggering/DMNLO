!1->2

	double precision function MHVS_width(s_in,m_a,m_b,Coupl_A,Coupl_Ac) 

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	double precision s_in,m_a,m_b
	complex*16 Coupl_A,Coupl_Ac


!===============================================Scalar-decay===========================================

!- from Vector-Vector-contraction

        MHVS_width = -dreal(Coupl_A*Coupl_Ac*(-m_a**2 + 2*(m_b**2 + s_in)) )
           

!      print*,'MHVS_width= ',MHVS_width

!======================================================================================================

      end


!2->2
!======================================================================================================

	double precision function MssHVS(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Higgs-exchange)========================================

!- from Vector-Vector-contraction

      MssHVS = dreal( -(2*sqrtS**2+2*m4**2-m3**2)*Coupl_A*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac)   )   
           

!      print*,'MssHVS= ',MssHVS

!======================================================================================================

      end




	double precision function MssVVS(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Vector-exchange)========================================

! (- from Vector-Propagator)^2 * -from Vector-Contraction

      MssVVS = - dreal( (Coupl_AM*(m1**2 - m2**2) + Coupl_AP*(2*m1**2 + 2*m2**2 - sqrtS**2))*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac)  )   
           

!      print*,'MssVVS= ',MssVVS

!======================================================================================================

      end




	double precision function MssVHVS(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Vector*Higgs-exchange)=================================

! - from Vector-Propagator* -from Vector-Contraction

      MssVHVS =  dreal( (Coupl_AM*(3*m1**2 + 3*m2**2 + 2*m3**2 + 4*m4**2 - 3*tman - 3*uman) + Coupl_AP*(3*m1**2 - 3*m2**2 + tman - uman))/2.
     &                  *Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac)  )           
           

!      print*,'MssVHVS= ',MssVHVS

!======================================================================================================

      end




	double precision function MssHVVS(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Higgs*Vector-exchange)=================================

!- from Vector-Vector-contraction, - from Vector-propagator 

      MssHVVS = (-uman+tman+3*m1**2- 3*m2**2)/2d0*dreal(Coupl_A*Coupl_B*Coupl_Ac*Coupl_Bc/(propa*propac)  )       
           

!      print*,'MssHVVS= ',MssHVVS

!======================================================================================================

      end




	double precision function MtVS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,sfm,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc
	double precision sfm,sfmc



!====================================t-channel (sfermion-exchange)=====================================

!- from Vector-contraction

      MtVS = -dreal( SS*(Coupl_AM*(m1**2 - tman) + Coupl_AP*(2*m1**2 - m3**2 + 2*tman))*Coupl_B*Coupl_Ac*Coupl_Bc/((tman-sfm**2)*(tman-sfmc**2))  )  
           

!      print*,'MtVS= ',MtVS

!======================================================================================================

      end




	double precision function MuVS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,sfm,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc
	double precision sfm,sfmc



!====================================u-channel (sfermion-exchange)=====================================

!- from Vector-contraction

      MuVS = -dreal( SS*(Coupl_AM*(-m2**2 + uman) + Coupl_AP*(2*m2**2 - m3**2 + 2*uman))*Coupl_B*Coupl_Ac*Coupl_Bc/((uman-sfm**2)*(uman-sfmc**2))  )  
           

!      print*,'MuVS= ',MuVS

!======================================================================================================

      end




	double precision function MtuVS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,sfm,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc
	double precision sfm,sfmc



!====================================tu-channel (sfermion-exchange)====================================

!- from Vector-contraction

      MtuVS = -dreal( SS*(Coupl_AP*(2*m1**2 + 2*m2**2 - m4**2 - sqrtS**2) + Coupl_AM*(-m2**2 + uman))*Coupl_B*Coupl_Ac*Coupl_Bc/((tman-sfm**2)*(uman-sfmc**2))  )  
           

!      print*,'MtuVS= ',MtuVS

!======================================================================================================

      end




	double precision function MutVS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,sfm,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc
	double precision sfm,sfmc



!====================================ut-channel (sfermion-exchange)====================================

!- from Vector-contraction

      MutVS = -dreal( SS*(Coupl_AP*(2*m1**2 + 2*m2**2 - m4**2 - sqrtS**2) - Coupl_AM*(tman-m1**2))*Coupl_B*Coupl_Ac*Coupl_Bc/((uman-sfm**2)*(tman-sfmc**2))  )  
           

!      print*,'MutVS= ',MutVS

!======================================================================================================

      end




	double precision function MstHVS(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa
	double precision sfmc



!====================================s*t-channel (Higgs-exchange)======================================

!- from Vector-Vector-contraction

      MstHVS = - (2*m1**2 - m2**2 + 2*m4**2 - uman)*dreal( SS*Coupl_A*Coupl_B*Coupl_Ac*Coupl_Bc/(propa*(tman-sfmc**2))   )  
           
!      print*,'MstHVS= ',MstHVS

!======================================================================================================

      end




	double precision function MsuHVS(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa
	double precision sfmc



!====================================s*u-channel (Higgs-exchange)======================================

!- from Vector-Vector-contraction

      MsuHVS = - (m1**2 - 2*m2**2 - 2*m4**2 + tman )*dreal( SS*Coupl_A*Coupl_B*Coupl_Ac*Coupl_Bc/(propa*(uman-sfmc**2))   )   
           
!      print*,'MsuHVS= ',MsuHVS

!======================================================================================================

      end




	double precision function MtsHVS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propac,sfm)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AM,Coupl_AP,Coupl_Ac,Coupl_B,Coupl_Bc,propac
	double precision sfm



!====================================t*s-channel (Higgs-exchange)======================================

!- from Vector-Vector-contraction

      MtsHVS = - dreal( SS*(Coupl_AM*(-m4**2 + sqrtS**2) + Coupl_AP*(2*m1**2 - m2**2 + 2*m4**2 - uman))*Coupl_B*Coupl_Ac*Coupl_Bc/(propac*(tman-sfm**2))   ) 
           
!      print*,'MtsHVS= ',MtsHVS

!======================================================================================================

      end




	double precision function MusHVS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa
	double precision sfmc



!====================================u*s-channel (Higgs-exchange)======================================

!- from Vector-Vector-contraction

      MusHVS = -dreal( SS*(Coupl_AM*(-m4**2 + sqrtS**2) + Coupl_AP*(m1**2 - 2*m2**2 - 2*m4**2 + tman))*Coupl_B*Coupl_Ac*Coupl_Bc/(propa*(uman-sfmc**2))   )  
           
!      print*,'MusHVS= ',MusHVS

!======================================================================================================

      end




	double precision function MtsVVS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propac,sfm)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propac
	double precision sfm



!====================================t*s-channel (Vector-exchange)========================================

!- from Vector-contraction - from Vector-Propagator 

      MtsVVS = dreal( ((Coupl_AM*(m1**2 - m2**2 - tman + uman) + Coupl_AP*(3*m1**2 + m2**2 - 2*m3**2 - 2*m4**2 + 3*tman + uman))/2.)
     &               *SS*Coupl_Ac*Coupl_B*Coupl_Bc/(propac*(tman-sfm**2))  )   
           

!      print*,'MtsVVS= ',MtsVVS

!======================================================================================================

      end




	double precision function MusVVS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propac,sfm)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propac
	double precision sfm



!====================================u*s-channel (Vector-exchange)========================================

!- from Vector-contraction - from Vector-Propagator 

      MusVVS = dreal( ((Coupl_AM*(m1**2 - m2**2 - tman + uman) + Coupl_AP*(m1**2 + 3*m2**2 - 2*m3**2 - 2*m4**2 + tman + 3*uman))/2.)
     &               *SS*Coupl_Ac*Coupl_B*Coupl_Bc/(propac*(uman-sfm**2))  )   
           

!      print*,'MusVVS= ',MusVVS

!======================================================================================================

      end




	double precision function MstVVS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa
	double precision sfmc



!====================================s*t-channel (Vector-exchange)========================================

!- from Vector-contraction * - from Vector-Propagator 

      MstVVS = dreal( (Coupl_AM*(3*m1**2 - m2**2 + 2*m4**2 - tman - uman) + Coupl_AP*(3*m1**2 + m2**2 - 2*m3**2 - 2*m4**2 + 3*tman + uman))/2.
     &               *SS*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*(tman-sfmc**2))  )      
           

!      print*,'MstVVS= ',MstVVS

!======================================================================================================

      end




	double precision function MsuVVS(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa
	double precision sfmc



!====================================s*u-channel (Vector-exchange)========================================

!- from Vector contraction * - from Vector-Propagator 

      MsuVVS = dreal( ((Coupl_AM*(m1**2 - 3*m2**2 - 2*m4**2 + tman + uman) + Coupl_AP*(m1**2 + 3*m2**2 - 2*m3**2 - 2*m4**2 + tman + 3*uman))/2.)
     &               *SS*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*(uman-sfmc**2))  )     
           

!      print*,'MsuVVS= ',MsuVVS

!======================================================================================================

      end
