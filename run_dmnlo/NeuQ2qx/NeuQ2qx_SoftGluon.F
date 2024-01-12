! =================================================================
        double precision function NeuQ2qx_SoftPhoton()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

        double precision IntegralIab

                    
        NeuQ2qx_SoftPhoton = -((EL**2*(2.d0/3.d0)**2)/(16d0*pi**3))*(IntegralIab(22) + IntegralIab(33) - 2d0*IntegralIab(23))       

        end

! =================================================================
	double precision function NeuQ2qx_HardCollGluon()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

        double precision Azero, Aepsilon, nf, cutswitch, switch
        double complex cspen 

        cutswitch = 1d0
        switch = 1d0
        
        Azero = 0d0
        Aepsilon = 0d0
        nf = 4d0
        
        if (xtype.eq.1) then
 
! gluon -> gluon gluon
        Aepsilon = CV*(switch*11d0/6d0 + cutswitch*dlog(4d0*cutoffS**2) - cutswitch*2d0*dlog(dsqrt(sman) - m3**2/dsqrt(sman)))
! ! without region m        
!         Azero = CV*(67d0/18d0 - pi**2/3d0 - (dlog(2d0*cutoffS) - dlog(dsqrt(sman) - m3**2/dsqrt(sman)))**2 - dlog(cutoffC)*(11d0/6d0 + dlog(4d0*cutoffS**2) - 2d0*dlog(dsqrt(sman) - m3**2/dsqrt(sman))))
!corrected, with region m    
        Azero = CV*(67d0/18d0 - pi**2/3d0 - (dlog(2d0*cutoffS) - dlog(dsqrt(sman) - m3**2/dsqrt(sman)))**2  - dlog(cutoffC)*(11d0/6d0 + dlog(4d0*cutoffS**2) - 2d0*dlog(dsqrt(sman) - m3**2/dsqrt(sman)))
     &     + 2d0*cspen(dcmplx(cutoffC/(2d0*cutoffS*dsqrt(sman)))) - 1d0/6d0)

! gluon -> quark quark
        Aepsilon = Aepsilon - (nf/3d0) 
        Azero = Azero + nf/3d0*(dlog(cutoffC) - 5d0/3d0)
        

      
        if (EpsPole.eq.-1d0) then
             NeuQ2qx_HardCollGluon = (GS**2/(8d0*pi**2))*Aepsilon
             
        else if (EpsPole.eq.0d0) then        
             NeuQ2qx_HardCollGluon = (GS**2/(8d0*pi**2))*(IRdiv*Aepsilon + Aepsilon*dlog(muSc**2) + Azero)	
        else 
             NeuQ2qx_HardCollGluon = 0d0
        endif	
        
        else
        
        NeuQ2qx_HardCollGluon = 0d0
        
        endif
        
        end


! =================================================================       
	double precision function NeuQ2qx_SoftGluon()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

        double precision IntegralIab, IntegralIsoft, IntegralIsoftcoll
         
        
! == only gluon final state ==
	if (xtype.eq.1) then

! only soft	
        NeuQ2qx_SoftGluon = 0
	NeuQ2qx_SoftGluon = -((GS**2)/(16d0*pi**3))*(CF*IntegralIsoft(22) + CF*IntegralIsoft(33) - 2d0*(CF-(CV/2d0))*IntegralIsoft(23))
	
! soft and collinear	
 	NeuQ2qx_SoftGluon = NeuQ2qx_SoftGluon - ((GS**2)/(16d0*pi**3))*CV*(1d0/sman)*(- 4d0*(m2**2 - tman)*IntegralIsoftcoll(24) - 2d0*(sman - m3**2)*IntegralIsoftcoll(35) - 2d0*(sman - m3**2)*IntegralIsoftcoll(45))
 	 	
!  	NeuQ2qx_SoftGluon = (((-GS**2*CV)/(16d0*pi**2))*dlog(1d0/m3**2)) + (((GS**2*CV)/(16d0*pi**2))*dlog(muSc**2))  + (((GS**2*CV)/(16d0*pi**2))*dlog(muSc**2))  - (((GS**2*CV)/(16d0*pi**2))*dlog((1d0*(m2**2-tman)**2)/(m2**2*(sman-m3**2)**(-2))))

	else
	
! == all final states except gluon ==
          NeuQ2qx_SoftGluon = -((GS**2*CF)/(16d0*pi**3))*(IntegralIab(22) + IntegralIab(33) - 2d0*IntegralIab(23))


        
        endif

        end

! =================================================================

	double precision function IntegralIsoftcoll(ij)

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"


        integer ij
        double precision IntegralI_11n, IntegralI_01n
        double precision aI, bI, AAI, BBI, BCI, costheta, kappa

	costheta = 2d0*sman*(tman - m1**2 - m3**2 + (sman + m1**2 - m2**2)*(sman + m3**2 - m4**2)/(2d0*sman))/(kappa(sman,m1**2, m2**2)*kappa(sman,m3**2,m4**2))

	
	if (ij.eq.24) then

	! E2
	aI = (sman + m2**2 - m1**2)/sman
	! beta1
	bI = kappa(sman,m1**2,m2**2)/sman

	
	! E4
	AAI = (sman - m3**2)/sman
	! beta2
	BBI = (sman - m3**2)/sman*costheta

	
	elseif (ij.eq.45) then
	
	! aI
	aI = 1d0
	! zero
	bI = 0d0
	
	! E4
	AAI = (sman - m3**2)/sman
	! beta2
	BBI = (sman - m3**2)/sman*costheta

	elseif (ij.eq.35) then
	
	! aI
	aI = 1d0
	! bI
	bI = 1d0
	
	!E3
	AAI = (sman + m3**2)/sman
	! BCI = dsqrt(B**2 + C**2) = +beta_2
	BCI = (sman - m3**2)/sman
	
	else
	
	aI = 1d0
	bI = 1d0
	
	AAI = 1d0
	BBI = 1d0
	BCI = 1d0
	
	endif
	
	if ((ij.eq.24).or.(ij.eq.45)) then

        IntegralIsoftcoll = IntegralI_11n(aI,bI,AAI,BBI)
        
        elseif (ij.eq.35) then

        IntegralIsoftcoll = IntegralI_01n(AAI,BCI)
        
        else

        IntegralIsoftcoll = 0d0

        endif

        end
    
    
! =================================================================
! Integral taken from Harris & Smith, Nucl. Phys. B 452 (1995)

	double precision function IntegralI_01n(AAI,BCI)

        implicit none
        
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"
 
        double precision AAI, BCI
        double precision zI,yI
        double complex cspen 

        yI = dlog((AAI + BCI)/(AAI - BCI))
        zI = 2d0*cspen(dcmplx((2d0*BCI)/(AAI + BCI))) + (1/2d0)*(dlog((AAI + BCI)/(AAI - BCI)))**2
        
        if (EpsPole.eq.(-2d0)) then  
           IntegralI_01n = 0d0
        elseif (EpsPole.eq.(-1d0)) then  
           IntegralI_01n = -pi/(BCI)*(1/2d0)*yI
        else
           IntegralI_01n = pi/(BCI)*(1/2d0)*(-IRdiv*yI + yI*dlog(4d0*cutoffS**2/muSc**2) - zI)        
        endif

        
        end   
     
    
! =================================================================
! Integral taken from Harris & Smith, Nucl. Phys. B 452 (1995)
! not needed for softcollinear cases
! 
! 	double precision function IntegralI_02n(AAI,BCI)
! 
!         implicit none
!         
! #include "NeuQ2qx_Kinematics.h"
! #include "NeuQ2qx_Model.h"
! #include "../util/DMNLO_LoopIntegrals.h"
!  
!         double precision AAI, BCI
!         double precision zI
!         double complex cspen 
! 
!         
!         zI = (AAI/BCI)*dlog((AAI + BCI)/(AAI - BCI))
!    
! 
!         IntegralI_02n = 2d0*pi/(AAI**2 - BCI**2)*(1/2d0)*(-IRdiv - zI + dlog(4d0*cutoffS**2/muSc**2))
!         
!         if (EpsPole.eq.(-2d0)) then  
!            IntegralI_02n = 0d0
!         elseif (EpsPole.eq.(-1d0)) then  
!            IntegralI_02n = -2d0*pi/(AAI**2 - BCI**2)*(1/2d0)
!         endif
! 
!         end   
!  
 
! =================================================================
! Integral recalculated for our case on the basis of the PhD thesis Bojak

	double precision function IntegralI_11n(aI,bI,AAI,BBI)
! Factor -1/(2epsilon) is defined into myX,myY,myZ (compare to eq. 9.60 of PhD Julia)
        implicit none
        
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"        
#include "../util/DMNLO_LoopIntegrals.h"

        double precision aI, bI, AAI, BBI,myX,myY,myZ,arg
        double complex cspen 

        arg = muSc**2/(4*cutoffS**2)
        
        myX = (pi/(AAI*aI - bI*BBI))*(1/2d0)
        myY = -(pi/(AAI*aI - bI*BBI))*dlog((aI*AAI - bI*BBI)**2/(AAI**2*(aI**2 - bI**2)))/2d0
        myZ = (pi/(AAI*aI - bI*BBI))*(-1d0)/(-2d0)*(dlog((aI-bI)*AAI/(AAI*aI - bI*BBI)))**2 - (1/2d0)*(dlog((aI + bI)/(aI - bI)))**2 + 2d0*cspen(dcmplx(bI*(BBI - AAI)/(AAI*(aI - bI)))) - 2d0*cspen(dcmplx(-bI*(AAI + BBI)/(AAI*aI - bI*BBI)))
        
        if (EpsPole.eq.(-2d0)) then  
            IntegralI_11n = myX
        elseif (EpsPole.eq.(-1d0)) then 
            IntegralI_11n = myY + myX*dlog(arg)
        else 
            IntegralI_11n =  myX*IRdiv**2 +(myY + myX*dlog(arg))*IRdiv+ myZ + myY*dlog(arg) + (myX*(-2*Pi**2+ 3*dlog(arg)**2))/6d0
        endif

        end   
        
    
! =================================================================

	double precision function IntegralIsoft(ij)

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

        integer ij
        double precision Ei, Ej, ki, kj, mi, mj, mi_, mj_
        double precision kikj, alpha_soft, v_soft, kappa
        double complex cspen
        

        external CMSenergies,SetMasses 

        call SetMasses(ij,mi,mj,mi_,mj_)
	call CMSenergies(mi,mj,mi_,mj_,Ei,Ej,ki,kj,ij,alpha_soft,kikj,v_soft)


	if ((ij.eq.11).or.(ij.eq.22).or.(ij.eq.33).or.(ij.eq.44)) then

       if (EpsPole.eq.-1d0) then             
        IntegralIsoft = -2d0*pi
       else if (EpsPole.eq.0d0) then             
        IntegralIsoft = 2d0*pi*((-IRdiv) + dlog(4d0*cutoff**2/muSc**2) + (Ei)/(ki)*dlog((Ei - ki)/(Ei + ki)))
       else
        IntegralIsoft = 0d0
       endif 
        
        else

       if (EpsPole.eq.-1d0) then             
        IntegralIsoft = (4d0*pi*alpha_soft*kikj)/(alpha_soft**2*mi**2 - mj**2)*(-1d0)*dlog(alpha_soft**2*mi**2/mj**2)/2d0
       else if (EpsPole.eq.0d0) then             
        IntegralIsoft = (4d0*pi*alpha_soft*kikj)/(alpha_soft**2*mi**2 - mj**2)*(((-IRdiv + 
     &                   dlog(4d0*cutoff**2/muSc**2))*dlog(alpha_soft**2*mi**2/mj**2))/2d0 + 
     &                   dreal((cdlog(dcmplx((Ei - ki)/(Ei + ki)))**2 - cdlog(dcmplx((Ej - kj)/(Ej + kj)))**2)/4d0 + 
     &                   cspen(dcmplx(1 - alpha_soft*(Ei - ki)/v_soft)) - cspen(dcmplx(1 - (Ej - kj)/v_soft)) + 
     &                   cspen(dcmplx(1 - alpha_soft*(Ei + ki)/v_soft)) - cspen(dcmplx(1 - (Ej + kj)/v_soft))))
       else
        IntegralIsoft = 0d0
       endif 

        endif

        end
        
! =================================================================

	double precision function IntegralIab(ij)

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

        integer ij
        double precision Ei, Ej, ki, kj, mi, mj, mi_, mj_
        double precision kikj, alpha_soft, v_soft
        double complex cspen

        external CMSenergies,SetMasses 


        call SetMasses(ij,mi,mj,mi_,mj_)

        
	call CMSenergies(mi,mj,mi_,mj_,Ei,Ej,ki,kj,ij,alpha_soft,kikj,v_soft)
	       

	if ((ij.eq.11).or.(ij.eq.22).or.(ij.eq.33).or.(ij.eq.44)) then

        IntegralIab = 2d0*pi*(-IRdiv + dlog(4d0*cutoff**2/muSc**2) + (Ei)/(ki)*dlog((Ei - ki)/(Ei + ki)))

        else

        IntegralIab = (4d0*pi*alpha_soft*kikj)/(alpha_soft**2*mi**2 - mj**2)*(((-IRdiv + dlog(4d0*cutoff**2/muSc**2))*dlog(alpha_soft**2*mi**2/mj**2))/2d0 + dreal((cdlog(dcmplx((Ei - ki)/(Ei + ki)))**2 - cdlog(dcmplx((Ej - kj)/(Ej + kj)))**2)/4d0 + cspen(dcmplx(1 - alpha_soft*(Ei - ki)/v_soft)) - cspen(dcmplx(1 - (Ej - kj)/v_soft)) + cspen(dcmplx(1 - alpha_soft*(Ei + ki)/v_soft)) - cspen(dcmplx(1 - (Ej + kj)/v_soft))))

        endif

        end
  
  
! =================================================================

        subroutine SetMasses(ij,mi,mj,mi_,mj_)       
        
        implicit none 

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"     

        integer ij
        double precision mi,mj,mi_,mj_

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

        end
        
        
! =================================================================

        subroutine CMSenergies(mi,mj,mi_,mj_,Ei,Ej,ki,kj,ij,alpha_soft,kikj,v_soft)

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

        integer ij
        double precision Ei, Ej, ki, kj, mi, mj, mi_, mj_
        double precision kikj, alpha_soft, alpha_soft_minus, alpha_soft_plus, v_soft, kappa


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

        alpha_soft_plus = (2d0*kikj + dsqrt(4d0*kikj**2 - 4d0*mi**2*mj**2))/(2*mi**2)
        alpha_soft_minus = (2d0*kikj - dsqrt(4d0*kikj**2 - 4d0*mi**2*mj**2))/(2*mi**2)


        if (((alpha_soft_plus*Ei - Ej)/(Ej)).gt.0d0) then

          alpha_soft = alpha_soft_plus

	else if (((alpha_soft_minus*Ei - Ej)/(Ej)).gt.0d0) then

          alpha_soft = alpha_soft_minus

        else 
             alpha_soft = 0d0      

        endif

         v_soft = (alpha_soft**2*mi**2 - mj**2)/(2d0*(alpha_soft*Ei - Ej))  


        end
