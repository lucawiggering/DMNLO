!***************************************************************************************************************************
	double precision function ChiChi2QQ_SoftGluon()

        implicit none

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

        double precision IntegralIabChiChi_

					
        ChiChi2QQ_SoftGluon = -((GS**2*CF)/(16d0*pi**3))*(IntegralIabChiChi_(33)  + IntegralIabChiChi_(44) - 2d0*IntegralIabChiChi_(34))	!HitchHiker Glg. 9.20+9.21 (el. Ladung -> Farbe) + only initial-state radiation included

c        print*,"Check5",ChiChi2QQ_SoftGluon,IntegralIabChiChi_(11),IntegralIabChiChi_(22) ,IntegralIabChiChi_(12),cutoff

        end

!****************************************************************************************************************************

	double precision function IntegralIabChiChi_(ij)

        implicit none

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

        integer ij
        double precision Ei, Ej, ki, kj, mi, mj, mi_, mj_
        double precision kikj, alpha_soft, v_soft
        double complex cspen

        external CMSenergiesChiChi_
        
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


!         print *, "cutoff in softgluon", cutoff

	call CMSenergiesChiChi_(mi,mj,mi_,mj_,Ei,Ej,ki,kj,ij,alpha_soft,kikj,v_soft)

	if ((ij.eq.11).or.(ij.eq.22).or.(ij.eq.33).or.(ij.eq.44)) then

          if(EpsPole.eq.-1d0) then

          IntegralIabChiChi_ = 2d0*pi*(-1d0)   !siehe HitchHiker Glg. 9.60 mit Glg. 9.85

          else

          IntegralIabChiChi_ = 2d0*pi*(-IRdiv + dlog(4d0*cutoff**2/muSc**2) + (Ei)/(ki)*dlog((Ei - ki)/(Ei + ki)))   !siehe HitchHiker Glg. 9.60 mit Glg. 9.85

          endif

!        print*,"Check6 ",pi,IRdiv, LEg, Ln4pi, dlog(4d0*cutoff**2/muSc**2) , (Ei)/(ki)*dlog((Ei - ki)/(Ei + ki)),muSc,cutoff

        else
        
          if(EpsPole.eq.-1d0) then

          IntegralIabChiChi_ = (4d0*pi*alpha_soft*kikj)/(alpha_soft**2*mi**2 - mj**2)*(-1d0 *dlog(alpha_soft**2*mi**2/mj**2))/2d0

          else

          IntegralIabChiChi_ = (4d0*pi*alpha_soft*kikj)/(alpha_soft**2*mi**2 - mj**2)*(((-IRdiv + dlog(4d0*cutoff**2/muSc**2))*dlog(alpha_soft**2*mi**2/mj**2))/2d0 + dreal((cdlog(dcmplx((Ei - ki)/(Ei + ki)))**2 - cdlog(dcmplx((Ej - kj)/(Ej + kj)))**2)/4d0 + cspen(dcmplx(1 - alpha_soft*(Ei - ki)/v_soft)) - cspen(dcmplx(1 - (Ej - kj)/v_soft)) + cspen(dcmplx(1 - alpha_soft*(Ei + ki)/v_soft)) - cspen(dcmplx(1 - (Ej + kj)/v_soft))))

          endif

        !siehe HitchHiker Glg. 9.59 mit Glg. 9.85

!        print*,"Check7 ",pi,alpha_soft,IRdiv, LEg, Ln4pi, dlog(4d0*cutoff**2/muSc**2) , (Ei)/(ki)*dlog((Ei - ki)/(Ei + ki)),muSc,cutoff,cspen(dcmplx(1 - alpha_soft*(Ei - ki)/v_soft))

        !allg. Struktur gemäß HitchHiker Glg. 9.68 und 9.58 oder 9.64 und 9.85

        endif

        end

!***************************************************************************************************************************************

        subroutine CMSenergiesChiChi_(mi,mj,mi_,mj_,Ei,Ej,ki,kj,ij,alpha_soft,kikj,v_soft)

        implicit none

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

        integer ij
        double precision Ei, Ej, ki, kj, mi, mj, mi_, mj_
        double precision kikj, alpha_soft, alpha_soft_minus, alpha_soft_plus, v_soft, kappa, sman

        sman = sqrtS**2

        Ei = (sman + mi**2 - mi_**2)/(2d0*sqrtS)    !Verbindung zum HitchHiker-Kapitel zum Phase-Space slicing: Ei=a^0, Ej=b^0, ki=a^i, kj=b^j, a,b sind die 4er-Imp. der abstrahlenden Teilchen 
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

        alpha_soft_plus = (2d0*kikj + dsqrt(4d0*kikj**2 - 4d0*mi**2*mj**2))/(2*mi**2)             !die zwei Lösungen der Bedingung 9.24 im HitchHiker
        alpha_soft_minus = (2d0*kikj - dsqrt(4d0*kikj**2 - 4d0*mi**2*mj**2))/(2*mi**2)


        if (((alpha_soft_plus*Ei - Ej)/(Ej)).gt.0d0) then   !Wahl des alphas gemäß HitchHiker Bedingung 9.25 

          alpha_soft = alpha_soft_plus

	else if (((alpha_soft_minus*Ei - Ej)/(Ej)).gt.0d0) then

          alpha_soft = alpha_soft_minus

        else 
             alpha_soft = 0d0      

        endif

         v_soft = (alpha_soft**2*mi**2 - mj**2)/(2d0*(alpha_soft*Ei - Ej))  


        end
