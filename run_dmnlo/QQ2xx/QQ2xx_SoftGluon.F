!***************************************************************************************************************************
	double precision function QQ2xx_SoftGluon()

        implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

        double precision IntegralIab_

					
        QQ2xx_SoftGluon = -((GS**2*CF)/(16d0*pi**3))*(IntegralIab_(11) + IntegralIab_(22) - 2d0*IntegralIab_(12))	!HitchHiker_v1.3 eq. 9.20 + 9.21 

        end

!****************************************************************************************************************************

	double precision function IntegralIab_(ij)

        implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

        integer ij
        double precision Ei, Ej, ki, kj, mi, mj, mi_, mj_
        double precision kikj, alpha_soft, v_soft
        double complex cspen

        external CMSenergies_
        
	set_masses: select case (ij)
	   case(11) set_masses 
	   	mi = m1
                mj = m1
                mi_= m2
                mj_= m2        
  	   case(22) set_masses
	   	mi = m2
                mj = m2
                mi_= m1
                mj_= m1
	   case(33) set_masses
	   	mi = m3
                mj = m3
                mi_= m4
                mj_= m4
	   case(44) set_masses 
	   	mi = m4
                mj = m4
                mi_= m3
                mj_= m3
	   case(12) set_masses
	   	mi = m1
                mj = m2
                mi_= m2
                mj_= m1
	   case(21) set_masses
	   	mi = m2
                mj = m1
                mi_= m1
                mj_= m2
	   case(13) set_masses
	   	mi = m1
                mj = m3
                mi_= m2
                mj_= m4
	   case(31) set_masses
	   	mi = m3
                mj = m1
                mi_= m4
                mj_= m2
	   case(14) set_masses
	   	mi = m1
                mj = m4
                mi_= m2
                mj_= m3
	   case(41) set_masses
	   	mi = m4
                mj = m1
                mi_= m3
                mj_= m2
	   case(23) set_masses
	   	mi = m2
                mj = m3
                mi_= m1
                mj_= m4
	   case(32) set_masses
	   	mi = m3
                mj = m2
                mi_= m4
                mj_= m1
	   case(24) set_masses
	   	mi = m2
                mj = m4
                mi_= m1
                mj_= m3
	   case(42) set_masses
	   	mi = m4
                mj = m2
                mi_= m3
                mj_= m1
	   case(34) set_masses
	   	mi = m3
                mj = m4
                mi_= m4
                mj_= m3
	   case(43) set_masses
	   	mi = m4
                mj = m3
                mi_= m3
                mj_= m4
	   case DEFAULT
		mi = 0d0
		mj = 0d0
                mi_= 0d0
                mj_= 0d0
	end select set_masses

	
	call CMSenergies_(mi,mj,mi_,mj_,Ei,Ej,ki,kj,ij,alpha_soft,kikj,v_soft)

	if ((ij.eq.11).or.(ij.eq.22).or.(ij.eq.33).or.(ij.eq.44)) then

        IntegralIab_ = 2d0*pi*(-IRdiv + dlog(4d0*cutoff**2/muSc**2) + (Ei)/(ki)*dlog((Ei - ki)/(Ei + ki)))              
        
        !from HitchHiker_v1.3 eq 9.60 + 9.85

        else

        IntegralIab_ = (4d0*pi*alpha_soft*kikj)/(alpha_soft**2*mi**2 - mj**2)*(((-IRdiv + dlog(4d0*cutoff**2/muSc**2))*dlog(alpha_soft**2*mi**2/mj**2))/2d0 + dreal((cdlog(dcmplx((Ei - ki)/(Ei + ki)))**2 - cdlog(dcmplx((Ej - kj)/(Ej + kj)))**2)/4d0 + cspen(dcmplx(1 - alpha_soft*(Ei - ki)/v_soft)) - cspen(dcmplx(1 - (Ej - kj)/v_soft)) + cspen(dcmplx(1 - alpha_soft*(Ei + ki)/v_soft)) - cspen(dcmplx(1 - (Ej + kj)/v_soft))))

        !from HitchHiker_v1.3 eq 9.59 + 9.85 

        !general structure from HitchHiker_v1.3 eq 9.68 + 9.58 + 9.64 + 9.85 

        endif

        end

!***************************************************************************************************************************************

        subroutine CMSenergies_(mi,mj,mi_,mj_,Ei,Ej,ki,kj,ij,alpha_soft,kikj,v_soft)

        implicit none

#include "QQ2xx_Kinematics.h"
#include "QQ2xx_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

        integer ij
        double precision Ei, Ej, ki, kj, mi, mj, mi_, mj_
        double precision kikj, alpha_soft, alpha_soft_minus, alpha_soft_plus, v_soft, kappa, sman

        sman = sqrtS**2

        Ei = (sman + mi**2 - mi_**2)/(2d0*sqrtS)    
        Ej = (sman + mj**2 - mj_**2 )/(2d0*sqrtS) 

        ki = (kappa(sman, mi**2, mi_**2))/(2d0*sqrtS)
        kj = (kappa(sman, mj**2, mj_**2))/(2d0*sqrtS)


	if ((ij.eq.12).or.(ij.eq.21).or.(ij.eq.34).or.(ij.eq.43)) then

        kikj = (sman - mi**2 - mj**2)/2d0

        else if ((ij.eq.13).or.(ij.eq.31).or.(ij.eq.24).or.(ij.eq.42)) then

        kikj = (mi**2 + mj**2 - tman)/2d0

        else if ((ij.eq.23).or.(ij.eq.32).or.(ij.eq.14).or.(ij.eq.41)) then

        kikj = (mi**2 + mj**2 - uman)/2d0

        endif

        alpha_soft_plus = (2d0*kikj + dsqrt(4d0*kikj**2 - 4d0*mi**2*mj**2))/(2*mi**2)             !HitchHiker_v1.3 eq. 9.24 
        alpha_soft_minus = (2d0*kikj - dsqrt(4d0*kikj**2 - 4d0*mi**2*mj**2))/(2*mi**2)


        if (((alpha_soft_plus*Ei - Ej)/(Ej)).gt.0d0) then                                         !HitchHiker_v1.3 eq. 9.25 

          alpha_soft = alpha_soft_plus

	else if (((alpha_soft_minus*Ei - Ej)/(Ej)).gt.0d0) then

          alpha_soft = alpha_soft_minus

        else 
             alpha_soft = 0d0      

        endif

         v_soft = (alpha_soft**2*mi**2 - mj**2)/(2d0*(alpha_soft*Ei - Ej))  


        end
