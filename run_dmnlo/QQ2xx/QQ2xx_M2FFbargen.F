!1->2

	double precision function MHFFbar_width(s_in,m_a,m_b,Coupl_AL,Coupl_AR,Coupl_ALc,Coupl_ARc) 

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	double precision s_in,m_a,m_b
	complex*16 Coupl_AL,Coupl_ALc,Coupl_AR,Coupl_ARc


!===============================================Scalar-decay===========================================


        MHFFbar_width = dreal(-(Coupl_AL*(2*Coupl_ARc*m_a*m_b + Coupl_ALc*(m_a**2 + m_b**2 - s_in))) - Coupl_AR*(2*Coupl_ALc*m_a*m_b + Coupl_ARc*(m_a**2 + m_b**2 - s_in)))           

!      print*,'MHFFbar_width= ',MHFFbar_width

!======================================================================================================

      end



	double precision function MVFFbar_width(s_in,m_a,m_b,Coupl_AL,Coupl_AR,Coupl_ALc,Coupl_ARc) 

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	double precision s_in,m_a,m_b
	complex*16 Coupl_AL,Coupl_ALc,Coupl_AR,Coupl_ARc


!===============================================Vector-decay===========================================

!- from Vector contraction

        MVFFbar_width = dreal(-2*(Coupl_AL*(-4*Coupl_ARc*m_a*m_b + Coupl_ALc*(m_a**2 + m_b**2 - s_in)) + Coupl_AR*(-4*Coupl_ALc*m_a*m_b + Coupl_ARc*(m_a**2 + m_b**2 - s_in))))           

!      print*,'MVFFbar_width= ',MVFFbar_width

!======================================================================================================

      end



!2->2
!======================================================================================================

	double precision function MssHFFbar(Coupl_A,Coupl_Ac,Coupl_BL,Coupl_BR,Coupl_BLc,Coupl_BRc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_BL,Coupl_BR,Coupl_BLc,Coupl_BRc,propa,propac



!====================================s-channel (Higgs-exchange)========================================


      MssHFFbar = dreal( ( -(Coupl_BLc*((m3**2 + m4**2 - sqrtS**2)*Coupl_BL + 2*m3*m4*Coupl_BR)) - Coupl_BRc
     &                           *(2*m3*m4*Coupl_BL + (m3**2 + m4**2 - sqrtS**2)*Coupl_BR) )*Coupl_A*Coupl_Ac
     &                           /(propa*propac) )
           

!      print*,'MssHFFbar= ',MssHFFbar

!======================================================================================================

      end




	double precision function MssVFFbar(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BL,Coupl_BR,Coupl_BLc,Coupl_BRc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BL,Coupl_BR,Coupl_BLc,Coupl_BRc,propa,propac



!====================================s-channel (Vector-exchange)========================================


      MssVFFbar= dreal( -(Coupl_Ac*(Coupl_AP*(Coupl_BL*(2*Coupl_BRc*m3*m4*(2*m1**2 + 2*m2**2 - sqrtS**2) + 
     &                       Coupl_BLc*(m3**4 + m4**4 - m4**2*sqrtS**2 + m3**2*(2*m4**2 - sqrtS**2 - 4*tman) + 4*m1**2*(m2**2 - tman) - 4*m2**2*tman - 4*m4**2*tman + 4*sqrtS**2*tman + 4*tman**2)) + 
     &                       Coupl_BR*(2*Coupl_BLc*m3*m4*(2*m1**2 + 2*m2**2 - sqrtS**2) + Coupl_BRc*
     &                       (m3**4 + m4**4 - m4**2*sqrtS**2 + m3**2*(2*m4**2 - sqrtS**2 - 4*tman) + 4*m1**2*(m2**2 - tman) - 4*m2**2*tman - 4*m4**2*tman + 4*sqrtS**2*tman + 4*tman**2))) + 
     &                       Coupl_AM*(Coupl_BL*(2*Coupl_BRc*(m1**2 - m2**2)*m3*m4 + Coupl_BLc*(2*m2**2*m3**2 + m3**4 - m3**2*(sqrtS**2 + 2*tman) + m4**2*(-2*m1**2 - m4**2 + sqrtS**2 + 2*tman))) + 
     &                       Coupl_BR*(2*Coupl_BLc*(m1**2 - m2**2)*m3*m4 + Coupl_BRc*(2*m2**2*m3**2 + m3**4 - m3**2*(sqrtS**2 + 2*tman) + m4**2*(-2*m1**2 - m4**2 + sqrtS**2 + 2*tman))))))/(propa*propac)   )
           

!      print*,'MssVFFbar= ',MssVFFbar

!======================================================================================================

      end




	double precision function MssVHFFbar(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BL,Coupl_BR,Coupl_BLc,Coupl_BRc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BL,Coupl_BR,Coupl_BLc,Coupl_BRc,propa,propac



!====================================s-channel (Vector*Higgs-exchange)==================================
!- from Vector-propagator

      MssVHFFbar= -dreal( (Coupl_Ac*((Coupl_AM - Coupl_AP)*(Coupl_BR*(Coupl_BRc*m3*(m2**2 + m4**2 - tman) - Coupl_BLc*m4*(m2**2 + m3**2 - uman)) + 
     &                        Coupl_BL*(Coupl_BLc*m3*(m2**2 + m4**2 - tman) - Coupl_BRc*m4*(m2**2 + m3**2 - uman))) + 
     &                        (Coupl_AM + Coupl_AP)*(Coupl_BL*(-(Coupl_BRc*m4*(m1**2 + m3**2 - tman)) + Coupl_BLc*m3*(m1**2 + m4**2 - uman)) + 
     &                        Coupl_BR*(-(Coupl_BLc*m4*(m1**2 + m3**2 - tman)) + Coupl_BRc*m3*(m1**2 + m4**2 - uman)))))/(propa*propac)  )
           

!      print*,'MssVHFFbar= ',MssVHFFbar

!======================================================================================================

      end




	double precision function MssHVFFbar(Coupl_A,Coupl_Ac,Coupl_BL,Coupl_BR,Coupl_BLc,Coupl_BRc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_BL,Coupl_BR,Coupl_BLc,Coupl_BRc,propa,propac



!====================================s-channel (Higgs*Vector-exchange)=================================
!- from Vector-propagator

      MssHVFFbar= -dreal(Coupl_A*Coupl_Ac*(Coupl_BRc*(Coupl_BR*m3*(m1**2 - m2**2 + tman - uman) + Coupl_BL*m4*(-m1**2 + m2**2 + tman - uman)) + 
     &                           Coupl_BLc*(Coupl_BL*m3*(m1**2 - m2**2 + tman - uman) + Coupl_BR*m4*(-m1**2 + m2**2 + tman - uman)))/(propa*propac)  )
           

!      print*,'MssHVFFbar= ',MssHVFFbar

!======================================================================================================

      end
