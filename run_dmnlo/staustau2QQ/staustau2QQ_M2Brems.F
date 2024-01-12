 !***************************************************************************************************************************
	 function staustau2QQ_M2Brems()

        implicit none

#include "staustau2QQ_Kinematics.h"
#include "staustau2QQ_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"
#include "staustau2QQ_GenCouplings.h"

      double precision Ihl, Ihl1, Ihl2, Ihl11, Ihl22, Ihl12, Ih1l2, Ih2l1, Colfac
      double precision g1l,g1r,g2l,g2r,kallen,beta0,beta1,beta2   
      
       complex*16 staustau2QQ_M2Brems
       complex*16 CSPEN
	complex*16 DenV(3),DenSH(2)    
       integer x,y

! ! Propagator
      DenV(1)=0d0
	DenV(2) = sqrtS**2 
	DenV(3) = (sqrtS**2-MZ**2)
	
	do x=1,2
	  DenSH(x) = (sqrtS**2-MHk(x)**2)	
	enddo       

      staustau2QQ_M2Brems = 0d0

	!!!---------------------------Denner's  integrals for (2 to 3)-process phase space integration-----------------------
	
	!!From Ansgar Denner: Habilitation thesis (1993), App. D: bremsstrahlung integrals
	
	!!  For the used integration subroutine of the DiffCS, the integrals need a prefactor of 1d0/(32d0*pi**3)
	!!  to compensate the different phase space element. Do not forget the correct flux factor.
	
	kallen =0d0
	beta0=0d0
	beta1=0d0
	beta2=0d0
	
	!define Källen function
	
	kallen = sqrt(sqrtS**4+muu3**4*sqrtS**4+muu4**4*sqrtS**4-2*muu3**2*muu4**2*sqrtS**4-2*muu3**2*sqrtS**4-2*muu4**2*sqrtS**4)
	beta0= (sqrtS**2-muu3**2*sqrtS**2-muu4**2*sqrtS**2+kallen)/(2*muu3*muu4*sqrtS**2)
	beta1= (sqrtS**2-muu3**2*sqrtS**2+muu4**2*sqrtS**2-kallen)/(2*muu3*sqrtS**2)
	beta2= (sqrtS**2+muu3**2*sqrtS**2-muu4**2*sqrtS**2-kallen)/(2*muu4*sqrtS**2)
	
	!!  btw: beta0*beta1*beta2=1d0
	
	!!h and l denote the upper and lower indices 0,1,2 arising from momenta p0,p1,p2 in num./denom.
	
	!CONVERGENT INTEGRALS
	
	Ihl=(1d0/(4d0*sqrtS**2))*(kallen/2d0*(sqrtS**2*(1+muu3**2+muu4**2))+2*sqrtS**4*muu3**2*dlog(beta2)+
     -  2*sqrtS**4*muu4**2*dlog(beta1)+2*sqrtS**4*muu3**2*muu4**2*dlog(beta0))
	
	Ihl1=(1d0/(4d0*sqrtS**2))*(-2*sqrtS**2*dlog(beta2)-2*sqrtS**2*muu4**2*dlog(beta0)-kallen) 
	Ihl2=(1d0/(4d0*sqrtS**2))*(-2*sqrtS**2*dlog(beta1)-2*sqrtS**2*muu3**2*dlog(beta0)-kallen) 
	Ih1l2=(1d0/(4d0*sqrtS**2))*(muu3**4*sqrtS**4*dlog(beta0)-sqrtS**2*(2*muu4**2-2*muu3**2-sqrtS**2)*dlog(beta1)-(kallen/4d0)*(muu4**2*sqrtS**2-
     -   3*muu3**2*sqrtS**2+5*sqrtS**2))
	Ih2l1=(1d0/(4d0*sqrtS**2))*(muu4**4*sqrtS**4*dlog(beta0)-sqrtS**2*(2*muu3**2-2*muu4**2-sqrtS**2)*dlog(beta2)-(kallen/4d0)*(muu3**2*sqrtS**2-
     -   3*muu4**2*sqrtS**2+5*sqrtS**2))    
     
     !!INTEGRALS WITH IR-SINGULARITY
     
     !!CSPEN is Spence function respectively the Eulerian dilogarithm Li_2
     
      if (EpsPole.eq.0d0) then      !! REGULAR PART
     
      Ihl11= (1d0)/(4d0*muu3**2*sqrtS**4)*(-kallen-(sqrtS**2-muu4**sqrtS**2)*dlog(beta0/beta2)-muu3**sqrtS**2*dlog(beta1)-
     -  kallen*((-IRdiv+dlog(muSc**2))/2d0+dlog(sqrtS**3*muu3*muu4/(kallen**2))))
     
      Ihl22= (1d0)/(4d0*muu4**2*sqrtS**4)*(-kallen-(sqrtS**2-muu3**sqrtS**2)*dlog(beta0/beta1)-muu3**sqrtS**2*dlog(beta2)-
     -  kallen*((-IRdiv+dlog(muSc**2))/2d0+dlog(sqrtS**3*muu3*muu4/(kallen**2))))
     
      Ihl12= 1d0/(4d0*sqrtS**2)*(-2d0*dlog(beta0)*((-IRdiv+dlog(muSc**2))/2d0+dlog(sqrtS**3*muu3*muu4/(kallen**2)))+
     -  2*(dlog(beta0))**2-(dlog(beta2))**2-(dlog(beta1))**2+2*CSPEN(DCMPLX(1-beta0**2))-CSPEN(DCMPLX(1-beta1**2))-CSPEN(DCMPLX(1-beta2**2)))  
     
      else if (Epspole.eq.-1d0) then  !! ONLY SINGULAR PART COMPENSATING THE INFRARED POLES IN LOOP FUNCTIONS
     
      Ihl11= (-kallen*(-IRdiv+dlog(muSc**2)))/(8d0*muu3**2*sqrtS**4)         !!responsible for IR singularities in counter terms
     
      Ihl22= (-kallen*(-IRdiv+dlog(muSc**2)))/(8d0*muu4**2*sqrtS**4)
     
      Ihl12= (-dlog(beta0)*(-IRdiv+dlog(muSc**2)))/(4d0*sqrtS**2)            !!responsible for IR singularities in vertex corrections
      
      else
      
      Ihl11= 0d0
     
      Ihl22= 0d0
     
      Ihl12= 0d0
      
      end if

! ------------------------------------------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------------------------------------------
   
      Colfac = 4d0   !Universal for this process

      
!***************************************************************************************************************************         
!! !---------Higgs squared for all possible gluon emissions------------------   (analytically checked with S.S.) -------
!***************************************************************************************************************************  

!       Do x=1,2
!       Do y=1,2
!       
!       !prefactor from dPS3 and flux, 1/(32*pi**3) compensates different definitions of PS3
! 
!       staustau2QQ_M2Brems =  staustau2QQ_M2Brems +1d0/(32d0*pi**3*2d0*Sqrt(sqrtS**4-4*muu1**2*sqrtS**4))*colfac*hff(x)*GS**2*HSf(x)*hff(y)*Hsf(y)/(DenSH(x)*DenSH(y))*  
!      -  ((8*muu3**2*sqrtS**4*Ihl11)- (16*muu3**4*sqrtS**4*Ihl11) -(16*muu3**3*muu4*sqrtS**4*Ihl11) -
!      - (16*muu3**2*sqrtS**2*Ihl1) -(16*muu3*muu4*sqrtS**2*Ihl1) -
!      - 16*Ih2l1 + (8*muu4**2*sqrtS**4*Ihl22)- (16*muu3**2*muu4**2*sqrtS**4*Ihl22) - (16*muu3*muu4**3*sqrtS**4*Ihl22) -
!      - (16*muu3*muu4*sqrtS**2*Ihl2)- (16*muu4**2*sqrtS**2*Ihl2) -
!      - (16*Ih1l2) - 2d0*((-4*sqrtS**4+16*muu3**2*sqrtS**4-16*muu3**4*sqrtS**4+8*muu3*muu4*sqrtS**4-16*muu3**3*muu4*sqrtS**4)*Ihl12 +
!      - (8*sqrtS**2-16*muu3**2*sqrtS**2-8*muu3*muu4*sqrtS**2+8*muu4**2*sqrtS**2)*Ihl2 +
!      - (8*sqrtS**2-8*muu3**2*sqrtS**2-8*muu3*muu4*sqrtS**2)*Ihl1 -16*Ihl))
! 
!       enddo
!       enddo

 !******************************************** 
 !--------------Vector squared------------------   
!********************************************** 
! ! 
! ! 
!       Do x=2,3
!       Do y=2,3
! 
!       g1l=BL(x,4,4,3,3)
!       g1r=BR(x,4,4,3,3)
!       g2l=BL(y,4,4,3,3)
!       g2r=BR(y,4,4,3,3)
! 
!        staustau2QQ_M2Brems =  staustau2QQ_M2Brems -1d0/(32d0*pi**3*2d0*Sqrt(sqrtS**4-4*muu1**2*sqrtS**4))*colfac*
!      -  (2*muu1**2+2*muu2**2-1)*sqrtS**2*GS**2*VSfSf(x,isf1,isf2,itt1,itt2,igen1,igen2)*VSfSf(y,isf1,isf2,itt1,itt2,igen1,igen2)/(DenV(x)*DenV(y))*
!      -  ((-8*(g1L*g2L+g1R*g2R)*muu3**2*sqrtS**4*Ihl11) + (16*(g1L*g2L+g1R*g2R)*muu3**4*sqrtS**4*Ihl11) -
!      -  (32*(g1R*g2L+g1L*g2R)*muu3**3*muu4*sqrtS**4*Ihl11)-(8*(g1L*g2L+g1R*g2R)*muu4**2*sqrtS**4*Ihl22) +
!      -  (16*(g1L*g2L+g1R*g2R)*muu3**2*muu4**2*sqrtS**4*Ihl22) -(32*(g1L*g2L+g1R*g2R)*muu3*muu4**3*sqrtS**4*Ihl22)+   !!!First 4 lines responsible for counter terms 
!      -  (8*(g1L*g2L+g1R*g2R)*sqrtS**4*Ihl12)-(32*(g1L*g2L+g1R*g2R)*muu3**2*sqrtS**4*Ihl12) +                          !!!Next 3 lines compensate vertex corrections 
!      -  (32*(g1L*g2L+g1R*g2R)*muu3**4*sqrtS**4*Ihl12)+
!      -  (32*(g1R*g2L+g1L*g2R)*muu3*muu4*sqrtS**4*Ihl12) - (64*(g1R*g2L+g1L*g2R)*muu3**3*muu4*sqrtS**4*Ihl12)  +
!      -  (16*(g1L*g2L+g1R*g2R)*muu3**2*sqrtS**2*Ihl1) -                                                                !!!FROM HERE: CONVERGENT
!      -  32*(g1R*g2L+g1L*g2R)*muu3*muu4*sqrtS**2*Ihl1-(16*(g1L*g2L+g1R*g2R)*sqrtS**2*Ihl1)  +
!      -  (32*(g1L*g2L+g1R*g2R)*muu3**2*sqrtS**2)*Ihl2 -(32*(g1R*g2L+g1L*g2R)*muu3*muu4*sqrtS**2*Ihl2) -
!      -  (16*(g1L*g2L+g1R*g2R)*muu4**2*sqrtS**2*Ihl2) -(32*(g1R*g2L+g1L*g2R)*muu3*muu4*sqrtS**2*Ihl2) +
!      -  (16*(g1L*g2L+g1R*g2R)*muu4**2*sqrtS**2*Ihl2) -(16*(g1L*g2L+g1R*g2R)*sqrtS**2*Ihl1) +
!      -  (16*(g1L*g2L+g1R*g2R)*muu3**2*sqrtS**2*Ihl1) -(32*(g1R*g2L+g1L*g2R)*muu3*muu4*sqrtS**2*Ihl1) +
!      -  (16*(g1L*g2L+g1R*g2R)*Ih2l1)+(16*(g1L*g2L+g1R*g2R)*Ih1l2))
! ! 
!       enddo 
!       enddo

      end
#include "../util/loop/DMNLO_LoopFunctions.F"
      
  
