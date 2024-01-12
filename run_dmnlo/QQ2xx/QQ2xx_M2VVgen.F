!1->2

	double precision function MHVV_width(s_in,m_a,m_b,Coupl_A,Coupl_Ac) 

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	double precision s_in,m_a,m_b
	complex*16 Coupl_A,Coupl_Ac


!===============================================Scalar-decay===========================================


        MHVV_width = dreal(4*Coupl_A*Coupl_Ac)
           

!      print*,'MHVV_width= ',MHVV_width

!======================================================================================================

      end



	double precision function MHGhGh_width(s_in,m_a,m_b,Coupl_A,Coupl_Ac) 

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	double precision s_in,m_a,m_b
	complex*16 Coupl_A,Coupl_Ac


!===============================================Scalar-decay===========================================


        MHGhGh_width = -2d0*dreal(Coupl_A*Coupl_Ac)
           

!      print*,'MHGhGh_width= ',MHGhGh_width

!======================================================================================================

      end



!2->2
!======================================================================================================

	double precision function MqqVV(Coupl_A,Coupl_Ac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac



!=================================================qq===================================================


      MqqVV = dreal(4*Coupl_A*Coupl_Ac )  
           

!      print*,'MqqVV= ',MqqVV    

!======================================================================================================

      end




	double precision function MssHVV(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Scalar-exchange)========================================


      MssHVV = dreal( 4*Coupl_A*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac) )    
           

!      print*,'MssHVV= ',MssHVV

!======================================================================================================

      end




	double precision function MssHGhGh(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Scalar-exchange)========================================
!Ghost-Ghost (see CalcHEP-Manual Appendix C + Peskin p. 515), cancels nonphysical polarizations included in pol. sum =-g_(mu,nu) (gauge-dependent)
!Ghost Antighost + Antighost Ghost:change of flux-direction gives two diagramms=> factor of two
!for ghost-masses see Peskin p. 743 eq 21.52 + 21.15. (in the Feynman-gauge (psi=1) they are the Gauge-Boson masses)
!For a calculation see Peskin p. 538 ff espacially eq. 16.40+16.42->Kinematics ~((p1-p2)*k1)*((p1-p2)*(-k2)) 
!different coupling for conjugated diagramm->treat diagramm and conjugated diagramm sepeartly with in general different (not just conjugated to each other) couplings


      MssHGhGh = - 2*dreal( Coupl_A*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac) )    
           

!      print*,'MssHGhGh= ',MssHGhGh

!======================================================================================================

      end




	double precision function MssVVV(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Vector-exchange)========================================

!(- from propagator)^2*(-from Vector-contraction)^2

      MssVVV = dreal( (-(Coupl_AP*(m1**4 + m2**4 + 10*m3**4 + 20*m3**2*m4**2 + 10*m4**4 - 2*m2**2*(m3**2 + m4**2) - 2*m1**2*(17*m2**2 + m3**2 + m4**2) - 18*m3**2*tman - 18*m4**2*tman + 3*tman**2 - 
     &         18*m3**2*uman - 18*m4**2*uman + 26*tman*uman + 3*uman**2)) + Coupl_AM*
     &     (-m1**4 + m2**4 - 5*(m3**2 - m4**2)*(tman - uman) + m1**2*(m3**2 + m4**2 + tman + uman) - m2**2*(m3**2 + m4**2 + tman + uman)))/2.*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac)  )
           

!      print*,'MssVVV= ',MssVVV            

!======================================================================================================

      end




	double precision function MssVGhGha(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Vector-exchange)========================================
!Ghost-Ghost (see CalcHEP-Manual Appendix C + Peskin p. 515), cancels nonphysical polarizations included in pol. sum =-g_(mu,nu) (gauge-dependent)
!Ghost Antighost + Antighost Ghost:change of flux-direction gives two diagramms=> factor of two
!for ghost-masses see Peskin p. 743 eq 21.52 + 21.15. (in the Feynman-gauge (psi=1) they are the Gauge-Boson masses)
!For a calculation see Peskin p. 538 ff espacially eq. 16.40+16.42->Kinematics ~((p1-p2)*k1)*((p1-p2)*(-k2)) 
!different coupling for conjugated diagramm->treat diagramm and conjugated diagramm sepeartly with in general different (not just conjugated to each other) couplings


!(- from propagator)^2*(-from Vector-contraction)^2

!s(general)*dconjg(s(tree))

      MssVGhGha = -2*dreal( (-((m1**2 - m2**2 + tman - uman)*(Coupl_AM*(m1**2 + m2**2 + 2*m3**2 - tman - uman) + Coupl_AP*(m1**2 - m2**2 - tman + uman)))/4.)*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac)  )
           
!      print*,'MssVGhGha= ',MssVGhGha            

!======================================================================================================

      end




	double precision function MssVGhGhb(Coupl_A,Coupl_APc,Coupl_AMc,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_APc,Coupl_AMc,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Vector-exchange)========================================
!Ghost-Ghost (see CalcHEP-Manual Appendix C + Peskin p. 515), cancels nonphysical polarizations included in pol. sum =-g_(mu,nu) (gauge-dependent)
!Ghost Antighost + Antighost Ghost:change of flux-direction gives two diagramms=> factor of two
!for ghost-masses see Peskin p. 743 eq 21.52 + 21.15. (in the Feynman-gauge (psi=1) they are the Gauge-Boson masses)
!For a calculation see Peskin p. 538 ff espacially eq. 16.40+16.42->Kinematics ~((p1-p2)*k1)*((p1-p2)*(-k2)) 
!different coupling for conjugated diagramm->treat diagramm and conjugated diagramm sepeartly with in general different (not just conjugated to each other) couplings


!(- from propagator)^2*(-from Vector-contraction)^2

!s(tree)*dconjg(s(general))

      MssVGhGhb = -2*dreal( (-((Coupl_AMc*(m1**2 + m2**2 + 2*m4**2 - tman - uman) + Coupl_APc*(m1**2 - m2**2 + tman - uman))*(m1**2 - m2**2 - tman + uman))/4.)
     &            *Coupl_A*Coupl_B*Coupl_Bc/(propa*propac)  )
           
!      print*,'MssVGhGhb= ',MssVGhGhb            

!======================================================================================================

      end




	double precision function MssVHVV(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Vector*Scalar)==========================================

!- from Vector-Propagator*(-from Vector-contraction)^2

      MssVHVV = - dreal( 3*(Coupl_AM*(-m3**2 + m4**2) + Coupl_AP*(tman - uman))*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac)  )   
           

!      print*,'MssVHVV= ',MssVHVV

!======================================================================================================

      end




	double precision function MssVHGhGh(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Vector*Scalar)==========================================
!Ghost-Ghost (see CalcHEP-Manual Appendix C + Peskin p. 515), cancels nonphysical polarizations included in pol. sum =-g_(mu,nu) (gauge-dependent)
!Ghost Antighost + Antighost Ghost:change of flux-direction gives two diagramms=> factor of two
!for ghost-masses see Peskin p. 743 eq 21.52 + 21.15. (in the Feynman-gauge (psi=1) they are the Gauge-Boson masses)
!For a calculation see Peskin p. 538 ff espacially eq. 16.40+16.42->Kinematics ~((p1-p2)*k1)*((p1-p2)*(-k2)) 
!different coupling for conjugated diagramm->treat diagramm and conjugated diagramm sepeartly with in general different (not just conjugated to each other) couplings

      if(((sftype1.eq.2).or.(sftype1.eq.4)).and.((sftype2.eq.1).or.(sftype2.eq.3))) then

      MssVHGhGh = -2*dreal( (Coupl_AM*(m1**2 + m2**2 + 2*m4**2 - tman - uman) + Coupl_AP*(m1**2 - m2**2 + tman - uman))/2.
     &                    *Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac)  )  

      else

      MssVHGhGh = -2*dreal( (Coupl_AM*(m1**2 + m2**2 + 2*m3**2 - tman - uman) + Coupl_AP*(m1**2 - m2**2 - tman + uman))/2.
     &                    *Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac)  )   

      endif
           

!      print*,'MssVHGhGh= ',MssVHGhGh

!======================================================================================================

      end




	double precision function MssHVVV(Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Scalar*Vector)==========================================

! - from Vector-Propagator*(-from Vector-contraction)^2

      MssHVVV = - 3*(tman - uman)*dreal( Coupl_A*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*propac)  )    
           

!      print*,'MssHVVV= ',MssHVVV

!======================================================================================================

      end




	double precision function MssHVGhGh(Coupl_A,Coupl_APc,Coupl_AMc,Coupl_B,Coupl_Bc,propa,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_APc,Coupl_AMc,Coupl_B,Coupl_Bc,propa,propac



!====================================s-channel (Scalar*Vector)==========================================
!Ghost-Ghost (see CalcHEP-Manual Appendix C + Peskin p. 515), cancels nonphysical polarizations included in pol. sum =-g_(mu,nu) (gauge-dependent)
!Ghost Antighost + Antighost Ghost:change of flux-direction gives two diagramms=> factor of two
!for ghost-masses see Peskin p. 743 eq 21.52 + 21.15. (in the Feynman-gauge (psi=1) they are the Gauge-Boson masses)
!For a calculation see Peskin p. 538 ff espacially eq. 16.40+16.42->Kinematics ~((p1-p2)*k1)*((p1-p2)*(-k2)) 
!different coupling for conjugated diagramm->treat diagramm and conjugated diagramm sepeartly with in general different (not just conjugated to each other) couplings

      if(((sftype1.eq.2).or.(sftype1.eq.4)).and.((sftype2.eq.1).or.(sftype2.eq.3))) then

      MssHVGhGh = -2*dreal( (-(Coupl_AMc*(m1**2 + m2**2 + 2*m3**2 - tman - uman) + Coupl_APc*(m1**2 - m2**2 - tman + uman))/2.)
     &                    *Coupl_A*Coupl_B*Coupl_Bc/(propa*propac)  ) 

      else 

      MssHVGhGh = -2*dreal( (-(Coupl_AMc*(m1**2 + m2**2 + 2*m4**2 - tman - uman) + Coupl_APc*(m1**2 - m2**2 + tman - uman))/2.)
     &                    *Coupl_A*Coupl_B*Coupl_Bc/(propa*propac)  )  

      endif           

!      print*,'MssHVGhGh= ',MssHVGhGh

!======================================================================================================

      end




	double precision function MttVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,man,manc,sfm,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc
	double precision man,manc,sfm,sfmc



!=================================================t*t-channel==========================================

!(-from Vector-contraction)^2

      MttVV = dreal( SS*(Coupl_AM*(m1**2 - tman) + Coupl_AP*(2*m1**2 - m3**2 + 2*tman))*(Coupl_BM*(-m2**2 + tman) + Coupl_BP*(2*m2**2 - m4**2 + 2*tman))*Coupl_Ac*Coupl_Bc/((man-sfm**2)*(manc-sfmc**2))  )   
           

!      print*,'MttVV= ',MttVV   

!======================================================================================================

      end




	double precision function MtuVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,man,manc,sfm,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc
	double precision man,manc,sfm,sfmc



!=================================================t*u-channel==========================================

!(-from Vector-contraction)^2

      MtuVV = dreal( SS*(Coupl_AM*(-m2**2 + uman) + Coupl_AP*(m1**2 + m2**2 - m3**2 - 2*m4**2 + tman + uman))*(Coupl_BM*(m1**2 - uman) + Coupl_BP*(m1**2 + m2**2 - 2*m3**2 - m4**2 + tman + uman))
     &              *Coupl_Ac*Coupl_Bc/((man-sfm**2)*(manc-sfmc**2))  )
           

!      print*,'MtuVV= ',MtuVV    

!======================================================================================================

      end




	double precision function MutVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,man,manc,sfm,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc
	double precision man,manc,sfm,sfmc



!=================================================u*t-channel==========================================


      MutVV = dreal( SS*(Coupl_BM*(m1**2 - tman) + Coupl_BP*(m1**2 + m2**2 - m3**2 - 2*m4**2 + tman + uman))*(Coupl_AM*(-m2**2 + tman) + Coupl_AP*(m1**2 + m2**2 - 2*m3**2 - m4**2 + tman + uman))
     &              *Coupl_Ac*Coupl_Bc/((man-sfm**2)*(manc-sfmc**2))  )
           

!      print*,'MutVV= ',MutVV    

!======================================================================================================

      end




	double precision function MuuVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,man,manc,sfm,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc
	double precision man,manc,sfm,sfmc



!=================================================u*u-channel==========================================


      MuuVV = dreal( SS*(Coupl_BP*(-2*m2**2 + m3**2 - 2*uman) + Coupl_BM*(m2**2 - uman))*(Coupl_AP*(-2*m1**2 + m4**2 - 2*uman) + Coupl_AM*(-m1**2 + uman))*Coupl_Ac*Coupl_Bc/((man-sfm**2)*(manc-sfmc**2))  )
           

!      print*,'MuuVV= ',MuuVV   

!======================================================================================================

      end




	double precision function MsqHVV(Coupl_A,Coupl_Ac,Coupl_B,propa)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_B,propa



!====================================s*q-channel (Scalar-exchange)========================================

!- from i-combination

      MsqHVV = -4*dreal( Coupl_A*Coupl_Ac*Coupl_B/propa )
           
!      print*,'MsqHVV= ',MsqHVV    

!======================================================================================================

      end




	double precision function MsqVVV(Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,propa)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,propa



!====================================s*q-channel (Vector-exchange)========================================

! - from Vector-Propagator * -1 from i-combination

      MsqVVV = dreal( 3*(Coupl_AM*(-m3**2 + m4**2) + Coupl_AP*(tman - uman))*Coupl_Ac*Coupl_B/propa  )  
           

!      print*,'MsqVVV= ',MsqVVV    

!======================================================================================================

      end




	double precision function MqsHVV(Coupl_A,Coupl_Ac,Coupl_Bc,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_Bc,propac



!====================================q*s-channel (Scalar-exchange)========================================

!- from i-combination

      MqsHVV = -4*dreal( Coupl_A*Coupl_Ac*Coupl_Bc/propac )   
           
!      print*,'MqsHVV= ',MqsHVV    

!======================================================================================================

      end




	double precision function MqsVVV(Coupl_A,Coupl_Ac,Coupl_Bc,propac)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 Coupl_A,Coupl_Ac,Coupl_Bc,propac



!====================================q*s-channel (Vector-exchange)========================================

!- from i-combination * -1 from propagator 

      MqsVVV = 3*(tman - uman)*dreal( Coupl_A*Coupl_Ac*Coupl_Bc/propac )  
           
!      print*,'MqsVVV= ',MqsVVV

!======================================================================================================

      end




	double precision function MtqVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,man,sfm)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM
	double precision man,sfm



!===============================================t*q-channel ===========================================

!- from i-combination

      MtqVV = -dreal( SS*(-(Coupl_AM*(Coupl_BM*(m3**2 + m4**2 - sqrtS**2) + Coupl_BP*(2*m2**2 + 3*m3**2 + m4**2 - sqrtS**2 - 2*uman))) + 
     &                      Coupl_AP*(Coupl_BP*(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*uman) + Coupl_BM*(2*m1**2 + m3**2 + 3*m4**2 - sqrtS**2 - 2*uman)))/2.
     &                  *Coupl_Ac/(man-sfm**2 ) )
           
!      print*,'MtqVV= ',MtqVV    

!======================================================================================================

      end




	double precision function MuqVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,man,sfm)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM
	double precision man,sfm



!==============================================u*q-channel ============================================

!- from i-combination

      MuqVV = -dreal( SS*(Coupl_AM*(-(Coupl_BM*(m3**2 + m4**2 - sqrtS**2)) + Coupl_BP*(2*m1**2 + m3**2 - m4**2 - sqrtS**2 - 2*uman)) + 
     &                    Coupl_AP*(-(Coupl_BM*(2*m2**2 - m3**2 + m4**2 - sqrtS**2 - 2*uman)) + Coupl_BP*(2*m1**2 + 2*m2**2 - m3**2 - m4**2 - sqrtS**2 + 4*uman)))/2.
     &                  *Coupl_Ac/(man-sfm**2 ) )
           
!      print*,'MuqVV= ',MuqVV        

!======================================================================================================

      end




	double precision function MqtVV(SS,Coupl_A,Coupl_Ac,Coupl_Bc,man,sfm)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_A,Coupl_Ac,Coupl_Bc
	double precision man,sfm



!==============================================q*t-channel ============================================

!- from i-combination

      MqtVV = -(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*uman)/2.*dreal( SS*Coupl_A*Coupl_Ac*Coupl_Bc/(man-sfm**2 ) )           
           
!      print*,'MqtVV= ',MqtVV

!======================================================================================================

      end




	double precision function MquVV(SS,Coupl_A,Coupl_Ac,Coupl_Bc,man,sfm)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_A,Coupl_Ac,Coupl_Bc
	double precision man,sfm



!==============================================q*u-channel ============================================

!- from i-combination

      MquVV = -(6*m1**2 + 6*m2**2 + 3*m3**2 + 3*m4**2 - 5*sqrtS**2 - 4*tman)/2.*dreal( SS*Coupl_A*Coupl_Ac*Coupl_Bc/(man-sfm**2 ) )          
           
!      print*,'MquVV= ',MquVV

!======================================================================================================

      end




	double precision function MstHVV(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,manc,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa
	double precision manc,sfmc



!======================================s*t-channel (Scalar-exchange)====================================


      MstHVV = dreal( SS*(m1**2 + m2**2 - 2*m3**2 - 2*m4**2 + 5*tman + uman)/2.*Coupl_A*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*(manc-sfmc**2 )))
           
!      print*,'MstHVV= ',MstHVV    

!======================================================================================================

      end




	double precision function MsuHVV(SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa,manc,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_A,Coupl_Ac,Coupl_B,Coupl_Bc,propa
	double precision manc,sfmc



!======================================s*u-channel (Scalar-exchange)====================================


      MsuHVV = dreal( SS*(m1**2 + m2**2 - 2*m3**2 - 2*m4**2 + tman + 5*uman)/2.*Coupl_A*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*(manc-sfmc**2 )))
           
!      print*,'MsuHVV= ',MsuHVV        

!======================================================================================================

      end




	double precision function MtsHVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,propa,manc,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,propa
	double precision manc,sfmc



!=======================================t*s-channel (Scalar-exchange)===================================


      MtsHVV = dreal( SS*(Coupl_AM*(Coupl_BM*(m1**2 + m2**2 - tman - uman) + Coupl_BP*(m1**2 - m2**2 - 2*m3**2 - tman + uman)) + 
     &    Coupl_AP*(Coupl_BM*(m1**2 - m2**2 + 2*m4**2 + tman - uman) + Coupl_BP*(m1**2 + m2**2 - 2*m3**2 - 2*m4**2 + 5*tman + uman)))/2.*Coupl_Ac*Coupl_Bc/(propa*(manc-sfmc**2 )))
           
!      print*,'MtsHVV= ',MtsHVV    

!======================================================================================================

      end




	double precision function MusHVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,propa,manc,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,propa
	double precision manc,sfmc



!==========================================u*s-channel (Scalar-exchange)================================


      MusHVV = dreal( SS*(Coupl_AM*(Coupl_BM*(m1**2 + m2**2 - tman - uman) + Coupl_BP*(m1**2 - m2**2 - 2*m4**2 + tman - uman)) + 
     &    Coupl_AP*(Coupl_BM*(m1**2 - m2**2 + 2*m3**2 - tman + uman) + Coupl_BP*(m1**2 + m2**2 - 2*m3**2 - 2*m4**2 + tman + 5*uman)))/2.*Coupl_Ac*Coupl_Bc/(propa*(manc-sfmc**2 )))
           
!      print*,'MusHVV= ',MusHVV    

!======================================================================================================

      end




	complex*16 function MstVVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,manc,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa
	double precision manc,sfmc



!=======================================s*t-channel (Vector-exchange)==================================

! - from Vector-Propagator * (- from Vector-contraction)^2

      MstVVV = - ( (Coupl_AP*(m1**4 + m2**4 + 4*m3**4 + 8*m3**2*m4**2 + 4*m4**4 - 8*m3**2*tman - 8*m4**2*tman + 5*tman**2 - 2*m1**2*(5*m2**2 + 2*m3**2 + tman - uman) - 2*m2**2*(2*m4**2 + tman - uman) - 
     &           4*m3**2*uman - 4*m4**2*uman + 2*tman*uman + uman**2) + Coupl_AM*(m1**4 - m2**4 + (m3**2 - m4**2)*(2*m3**2 + 2*m4**2 - 3*tman - uman) - m2**2*(m3**2 - 3*m4**2 - 3*tman + uman) + 
     &           m1**2*(-3*m3**2 + m4**2 - 3*tman + uman)))/2.*SS*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*(manc-sfmc**2))  )  
           

!      print*,'MstVVV= ',MstVVV    

!======================================================================================================

      end




	complex*16 function MsuVVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa,manc,sfmc)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_B,Coupl_Bc,propa
	double precision manc,sfmc



!====================================s*u-channel (Vector-exchange)=====================================

! - from Vector-Propagator * (- from Vector-contraction)^2

      MsuVVV = -( (Coupl_AM*(-m1**4 + m2**4 + (m3**2 - m4**2)*(2*m3**2 + 2*m4**2 - tman - 3*uman) - m1**2*(m3**2 - 3*m4**2 + tman - 3*uman) + m2**2*(-3*m3**2 + m4**2 + tman - 3*uman)) - 
     &          Coupl_AP*(m1**4 + m2**4 + 4*m3**4 + 8*m3**2*m4**2 + 4*m4**4 - 4*m3**2*tman - 4*m4**2*tman + tman**2 - 8*m3**2*uman - 8*m4**2*uman + 2*tman*uman + 5*uman**2 - 
     &          2*m2**2*(2*m3**2 - tman + uman) - 2*m1**2*(5*m2**2 + 2*m4**2 - tman + uman)))/2.*SS*Coupl_Ac*Coupl_B*Coupl_Bc/(propa*(manc-sfmc**2))  )  
           

!      print*,'MsuVVV= ',MsuVVV    

!======================================================================================================

      end




	complex*16 function MtsVVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,propac,man,sfm)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,propac
	double precision man,sfm



!====================================t*s-channel (Vector-exchange)=====================================

! - from Vector-Propagator * (- from Vector-contraction)^2

      MtsVVV = -( SS*(Coupl_AP*(Coupl_BM*(m1**4 - 2*m3**2*m4**2 - 2*m4**4 + m1**2*(7*m2**2 - 2*m3**2 - m4**2 - tman) + 2*m3**2*tman + 5*m4**2*tman - 2*tman**2 + m2**2*(-2*m3**2 + m4**2 + 2*tman - uman) + 
     &          2*m3**2*uman + 3*m4**2*uman - 5*tman*uman - uman**2) + Coupl_BP*(m1**4 + m2**4 + 4*m3**4 + 8*m3**2*m4**2 + 4*m4**4 - 8*m3**2*tman - 8*m4**2*tman + 5*tman**2 - 
     &          2*m1**2*(5*m2**2 + 2*m3**2 + tman - uman) - 2*m2**2*(2*m4**2 + tman - uman) - 4*m3**2*uman - 4*m4**2*uman + 2*tman*uman + uman**2)) + 
     &          Coupl_AM*(Coupl_BP*(-m2**4 + 2*m3**4 + 2*m3**2*m4**2 - 5*m3**2*tman - 2*m4**2*tman + 2*tman**2 + m2**2*(m3**2 + 2*m4**2 + tman) - 3*m3**2*uman - 2*m4**2*uman + 5*tman*uman + uman**2 + 
     &          m1**2*(-7*m2**2 - m3**2 + 2*m4**2 - 2*tman + uman)) + Coupl_BM*(-((m3**2 + m4**2 - tman - uman)*(tman - uman)) + m2**2*(m3**2 - m4**2 - tman + uman) + 
     &          m1**2*(-m3**2 + m4**2 - tman + uman))))/2.*Coupl_Ac*Coupl_Bc/(propac*(man-sfm**2))  )  
           

!      print*,'MtsVVV= ',MtsVVV   

!======================================================================================================

      end




	complex*16 function MusVVV(SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,propac,man,sfm)

	implicit none
	
#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"

	
	complex*16 SS,Coupl_AP,Coupl_AM,Coupl_Ac,Coupl_BP,Coupl_BM,Coupl_Bc,propac
	double precision man,sfm



!====================================u*s-channel (Vector-exchange)=====================================

!- from Vector-Propagator * (- from Vector-contraction)^2

      MusVVV = -( SS*(Coupl_AP*(Coupl_BM*(-m1**4 + 2*m3**4 + 2*m3**2*m4**2 - 3*m3**2*tman - 2*m4**2*tman + tman**2 + m2**2*(-m3**2 + 2*m4**2 + tman - 2*uman) - 5*m3**2*uman - 2*m4**2*uman + 5*tman*uman + 2*uman**2 + 
     &          m1**2*(-7*m2**2 + m3**2 + 2*m4**2 + uman)) - Coupl_BP*(m1**4 + m2**4 + 4*m3**4 + 8*m3**2*m4**2 + 4*m4**4 - 4*m3**2*tman - 4*m4**2*tman + tman**2 - 8*m3**2*uman - 8*m4**2*uman + 
     &          2*tman*uman + 5*uman**2 - 2*m2**2*(2*m3**2 - tman + uman) - 2*m1**2*(5*m2**2 + 2*m4**2 - tman + uman))) + 
     &          Coupl_AM*(Coupl_BM*(-((m3**2 + m4**2 - tman - uman)*(tman - uman)) + m2**2*(m3**2 - m4**2 - tman + uman) + m1**2*(-m3**2 + m4**2 - tman + uman)) + 
     &          Coupl_BP*(m2**4 - 2*m3**2*m4**2 - 2*m4**4 + 2*m3**2*tman + 3*m4**2*tman - tman**2 + 2*m3**2*uman + 5*m4**2*uman - 5*tman*uman - 2*uman**2 - m2**2*(2*m3**2 + m4**2 + uman) + 
     &          m1**2*(7*m2**2 - 2*m3**2 + m4**2 - tman + 2*uman))))/2.*Coupl_Ac*Coupl_Bc/(propac*(man-sfm**2))  )  
           

!      print*,'MusVVV= ',MusVVV

!======================================================================================================

      end
