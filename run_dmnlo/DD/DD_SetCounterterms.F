!==================================================================================================
!=======================	Calculates all needed coupling counterterms	   ========================
!======================= Copy and expansion of ChiChi2QQ_SetCounterterms.F ========================
!==================================================================================================

	subroutine DD_SetCounterTerms
	implicit none

#include "DD_Model.h"
#include "DD_GenCouplings.h"
#include "DD_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"


	integer a,i,j,k,l
	integer icha, indi, indj, tt, gen
	double complex DLred(6),DRred(6)
	double precision dMf_high(4,3)

! functions from DMNLO_Kinematics.F
	double precision Krondelta, ef, I3f, hf
	
! default settings for direct detection
	icha = 0
	indi = 1
	indj = 1

! Loop over tt and gen
	do tt = 3,4
	do gen = 1,3
	
! Check mixing angle renormalization constants for massless case
! Notwendig um 0*NaN aus Renscheme in diesem Fall zu korrigieren
   
	if(Mf(tt,gen).eq.0d0) then
		do k=1,2
			do l=1,2
				dUSf(k, l, tt, gen) = (0.0d0,0.0d0)
			enddo
	    enddo
	endif

! Use finite masses in mass counterterms
	dMf_high(tt,gen) = Mf_high(tt,gen)*GS**2*CF/(16d0*pi**2)*(-2d0*UVdiv)
	
!	write(*,*)'tt,gen',tt,gen
!	write(*,*)'dMf(tt,gen)',dMf(tt,gen)
!	write(*,*)'dMf_high(tt,gen)',dMf_high(tt,gen)
	
! And deactivate this stuff for the top quark which might be on-shell but anyway massive
	dMf_high(3,3) = dMf(3,3)
	
! Here only the SUSY contributions are included in the counterterm. This is needed for the Z-Gluino Vertex!
	call DMNLO_CalcPif_SUSY(Mf(tt,gen)**2,tt,gen)
	call DMNLO_CalcdPif_SUSY(Mf(tt,gen)**2,tt,gen)
	dZfL_SUSY(tt,gen) = -PifVL_SUSY(tt,gen) - Mf2(tt,gen)*(dPifVL_SUSY(tt,gen) + dPifVR_SUSY(tt,gen)) - Mf(tt,gen)*dPifSL_SUSY(tt,gen) - Mf(tt,gen)*dPifSR_SUSY(tt,gen)
	dZfR_SUSY(tt,gen) = -PifVR_SUSY(tt,gen) - Mf2(tt,gen)*(dPifVR_SUSY(tt,gen) + dPifVL_SUSY(tt,gen)) - Mf(tt,gen)*dPifSR_SUSY(tt,gen) - Mf(tt,gen)*dPifSL_SUSY(tt,gen)
	
!************************
! Coupling counterterms *
!************************

        if(icha.eq.1) then

! W+-f-f

	dBL(tt,gen) = BL(tt,gen)*(1/2d0*dZfL(tt+1,gen) + 1/2d0*dconjg(dZfL(tt,gen)))
	dBR(tt,gen) = BR(tt,gen)*(1/2d0*dZfR(tt+1,gen) + 1/2d0*dconjg(dZfR(tt,gen)))

	dBpL(tt,gen) = (0d0,0d0)
	dBpR(tt,gen) = (0d0,0d0)

        else

! Z-f-f

	dBL(tt,gen) = BL(tt,gen)*(1/2d0*dZfL_SUSY(tt,gen) + 1/2d0*dconjg(dZfL_SUSY(tt,gen)))
	dBR(tt,gen) = BR(tt,gen)*(1/2d0*dZfR_SUSY(tt,gen) + 1/2d0*dconjg(dZfR_SUSY(tt,gen)))

	dBpL(tt,gen) = BpL(tt,gen)*(1/2d0*dZfL(tt,gen) + 1/2d0*dconjg(dZfL(tt,gen)))
	dBpR(tt,gen) = BpL(tt,gen)*(1/2d0*dZfR(tt,gen) + 1/2d0*dconjg(dZfR(tt,gen)))

        endif
	
	if (icha.eq.1) then
		
! FerFerHiggs
	 if ((tt.eq.1).or.(tt.eq.3)) then

  	  DLred(1) = -EL*CA/(SW*2d0*MW*SB) 
  	  DRred(1) = -EL*CA/(SW*2d0*MW*SB) 
  	  DLred(2) = -EL*SA/(SW*2d0*MW*SB) 
  	  DRred(2) = -EL*SA/(SW*2d0*MW*SB) 
  	  DLred(3) = -cI*EL*CB/(SW*2d0*MW*SB) 
  	  DRred(3) = cI*EL*CB/(SW*2d0*MW*SB) 
  	  DLred(4) = -cI*EL/(SW*2d0*MW) 
  	  DRred(4) = cI*EL/(SW*2d0*MW) 

	 else

	  DLred(1) = EL*SA/(SW*2d0*MW*CB) 
	  DRred(1) = EL*SA/(SW*2d0*MW*CB) 
	  DLred(2) = -EL*CA/(SW*2d0*MW*CB) 
	  DRred(2) = -EL*CA/(SW*2d0*MW*CB) 
	  DLred(3) = -cI*EL*SB/(SW*2d0*MW*CB)
 	  DRred(3) = cI*EL*SB/(SW*2d0*MW*CB) 
	  DLred(4) = cI*EL/(SW*2d0*MW)
 	  DRred(4) = -cI*EL/(SW*2d0*MW) 

	 endif

! Fer1Fer2Higgs (tt = type of the up-type fermion)
	DLred(5) = EL*CB/(SW*sqrt2*MW*SB)
	DRred(5) = EL*SB/(SW*sqrt2*MW*CB)
	DLred(6) = EL/(SW*sqrt2*MW)
	DRred(6) = -EL/(SW*sqrt2*MW)

	else
! Yukawa coupling of the neutral Higgses to fermions

! FerFerHiggs
	if ((tt.eq.1).or.(tt.eq.3)) then

	DLred(1) = -EL*CA/(SW*2d0*MW*SB) 
	DRred(1) = -EL*CA/(SW*2d0*MW*SB) 
	DLred(2) = -EL*SA/(SW*2d0*MW*SB) 
	DRred(2) = -EL*SA/(SW*2d0*MW*SB) 
	DLred(3) = -cI*EL*CB/(SW*2d0*MW*SB) 
	DRred(3) = cI*EL*CB/(SW*2d0*MW*SB) 
	DLred(4) = -cI*EL/(SW*2d0*MW) 
	DRred(4) = cI*EL/(SW*2d0*MW) 	
	
	else

	DLred(1) = EL*SA/(SW*2d0*MW*CB) 
	DRred(1) = EL*SA/(SW*2d0*MW*CB) 
	DLred(2) = -EL*CA/(SW*2d0*MW*CB) 
	DRred(2) = -EL*CA/(SW*2d0*MW*CB) 
	DLred(3) = -cI*EL*SB/(SW*2d0*MW*CB)
 	DRred(3) = cI*EL*SB/(SW*2d0*MW*CB) 
	DLred(4) = cI*EL/(SW*2d0*MW)												
 	DRred(4) = -cI*EL/(SW*2d0*MW) 

	endif

! Fer1Fer2Higgs
	DLred(5) = 0d0
	DRred(5) = 0d0
	DLred(6) = 0d0
	DRred(6) = 0d0
	
	endif

! H0-f-f
	dDL(1,tt,gen) = DLred(1)*(dMf_high(tt,gen) + 1/2d0*dZfL(tt,gen)*Mf_high(tt,gen) + 1/2d0*dconjg(dZfR(tt,gen))*Mf_high(tt,gen))
	dDR(1,tt,gen) = DRred(1)*(dMf_high(tt,gen) + 1/2d0*dZfR(tt,gen)*Mf_high(tt,gen) + 1/2d0*dconjg(dZfL(tt,gen))*Mf_high(tt,gen))
	dDL(2,tt,gen) = DLred(2)*(dMf_high(tt,gen) + 1/2d0*dZfL(tt,gen)*Mf_high(tt,gen) + 1/2d0*dconjg(dZfR(tt,gen))*Mf_high(tt,gen))
	dDR(2,tt,gen) = DRred(2)*(dMf_high(tt,gen) + 1/2d0*dZfR(tt,gen)*Mf_high(tt,gen) + 1/2d0*dconjg(dZfL(tt,gen))*Mf_high(tt,gen))
	dDL(3,tt,gen) = DLred(3)*(dMf_high(tt,gen) + 1/2d0*dZfL(tt,gen)*Mf_high(tt,gen) + 1/2d0*dconjg(dZfR(tt,gen))*Mf_high(tt,gen))
	dDR(3,tt,gen) = DRred(3)*(dMf_high(tt,gen) + 1/2d0*dZfR(tt,gen)*Mf_high(tt,gen) + 1/2d0*dconjg(dZfL(tt,gen))*Mf_high(tt,gen))
	dDL(4,tt,gen) = DLred(4)*(dMf_high(tt,gen) + 1/2d0*dZfL(tt,gen)*Mf_high(tt,gen) + 1/2d0*dconjg(dZfR(tt,gen))*Mf_high(tt,gen))
	dDR(4,tt,gen) = DRred(4)*(dMf_high(tt,gen) + 1/2d0*dZfR(tt,gen)*Mf_high(tt,gen) + 1/2d0*dconjg(dZfL(tt,gen))*Mf_high(tt,gen)) 
	
	if (icha.eq.1) then
		
	dDL(5,tt,gen) = DLred(5)*(dMf_high(tt,gen) + 1/2d0*dZfL(tt+1,gen)*Mf_high(tt,gen) + 1/2d0*dconjg(dZfR(tt,gen))*Mf_high(tt,gen))
	dDR(5,tt,gen) = DRred(5)*(dMf_high(tt+1,gen) + 1/2d0*dZfR(tt+1,gen)*Mf_high(tt+1,gen) + 1/2d0*dconjg(dZfL(tt,gen))*Mf_high(tt+1,gen)) 
	dDL(6,tt,gen) = DLred(6)*(dMf_high(tt,gen) + 1/2d0*dZfL(tt+1,gen)*Mf_high(tt,gen) + 1/2d0*dconjg(dZfR(tt,gen))*Mf_high(tt,gen))
	dDR(6,tt,gen) = DRred(6)*(dMf_high(tt+1,gen) + 1/2d0*dZfR(tt+1,gen)*Mf_high(tt+1,gen) + 1/2d0*dconjg(dZfL(tt,gen))*Mf_high(tt+1,gen)) 

	endif

! Neu-Sf-f

	if (icha.eq.0) then

	do a=1,4
	 do i =1,2

	if ((tt.eq.1).or.(tt.eq.3)) then

	dELc(a,i,tt,gen) = -EL*(1d0/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*(dMf_high(tt,gen)*USf(i,1,tt,gen) + Mf_high(tt,gen)*dUSf(i,1,tt,gen)) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*(dUSf(i,2,tt,gen))) +
     & 	                 ELc(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) + ELc(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) + ELc(a,i,tt,gen)*1/2d0*dconjg(dZfR(tt,gen))

	dER(a,i,tt,gen) =  -EL*(1d0/(SW*sqrt2*MW*SB)*ZNeu(a,4)*(dMf_high(tt,gen)*USf(i,2,tt,gen) + Mf_high(tt,gen)*dUSf(i,2,tt,gen)) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*(dUSf(i,1,tt,gen))) +
     & 	                  ER(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  ER(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  ER(a,i,tt,gen)*1/2d0*dconjg(dZfL(tt,gen)) 

	dFL(a,i,tt,gen) =  -EL*(1d0/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*(dMf_high(tt,gen)*USf(i,2,tt,gen) + Mf_high(tt,gen)*dUSf(i,2,tt,gen)) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeuC(a,1) + I3f(tt)*ZNeuC(a,2))*(dUSf(i,1,tt,gen))) +
     & 	                  FL(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  FL(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  FL(a,i,tt,gen)*1/2d0*dZfL(tt,gen) 

	dFR(a,i,tt,gen) = -EL*(1d0/(SW*sqrt2*MW*SB)*ZNeu(a,4)*(dMf_high(tt,gen)*USf(i,1,tt,gen) + Mf_high(tt,gen)*dUSf(i,1,tt,gen)) -
     &   			sqrt2*ef(tt)/CW*ZNeu(a,1)*(dUSf(i,2,tt,gen))) +
     & 	                  FR(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  FR(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  FR(a,i,tt,gen)*1/2d0*dZfR(tt,gen)
		
	else

	dELc(a,i,tt,gen) = -EL*(1d0/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*(dMf_high(tt,gen)*USf(i,1,tt,gen) + Mf_high(tt,gen)*dUSf(i,1,tt,gen)) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*(dUSf(i,2,tt,gen))) +
     & 	                 ELc(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) + ELc(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) + ELc(a,i,tt,gen)*1/2d0*dconjg(dZfR(tt,gen))

	dER(a,i,tt,gen) =  -EL*(1d0/(SW*sqrt2*MW*CB)*ZNeu(a,3)*(dMf_high(tt,gen)*USf(i,2,tt,gen) + Mf_high(tt,gen)*dUSf(i,2,tt,gen)) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*(dUSf(i,1,tt,gen))) +
     & 	                  ER(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  ER(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  ER(a,i,tt,gen)*1/2d0*dconjg(dZfL(tt,gen)) 

	dFL(a,i,tt,gen) =  -EL*(1d0/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*(dMf_high(tt,gen)*USf(i,2,tt,gen) + Mf_high(tt,gen)*dUSf(i,2,tt,gen)) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeuC(a,1) + I3f(tt)*ZNeuC(a,2))*(dUSf(i,1,tt,gen))) +
     & 	                  FL(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  FL(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  FL(a,i,tt,gen)*1/2d0*dZfL(tt,gen) 

	dFR(a,i,tt,gen) = -EL*(1d0/(SW*sqrt2*MW*CB)*ZNeu(a,3)*(dMf_high(tt,gen)*USf(i,1,tt,gen) + Mf_high(tt,gen)*dUSf(i,1,tt,gen)) -
     &   			sqrt2*ef(tt)/CW*ZNeu(a,1)*(dUSf(i,2,tt,gen))) +
     & 	                  FR(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  FR(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  FR(a,i,tt,gen)*1/2d0*dZfR(tt,gen)

	endif

	dGL(a,i,tt,gen) =   dFL(a,i,tt,gen)
	dGR(a,i,tt,gen) =   dFR(a,i,tt,gen)
	dHL(a,i,tt,gen) =  dELc(a,i,tt,gen)
	dHR(a,i,tt,gen) =   dER(a,i,tt,gen)
	
	 enddo
	enddo

	 else if (icha.eq.1) then

	 do i= 1,2
	  do a= 1,4
   
	  if ((tt.eq.1).or.(tt.eq.3)) then
   
	dHL(a,i,tt,gen) = -EL*(1d0/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*(USf(i,1,tt,gen)*dMf_high(tt,gen) + Mf_high(tt,gen)*dUSf(i,1,tt,gen)) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*(dUSf(i,2,tt,gen))) +
     & 	                 HL(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) + HL(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) + HL(a,i,tt,gen)*1/2d0*dconjg(dZfR(tt,gen))

	dHR(a,i,tt,gen) =  -EL*(1d0/(SW*sqrt2*MW*SB)*ZNeu(a,4)*(USf(i,2,tt,gen)*dMf_high(tt,gen) + Mf_high(tt,gen)*dUSf(i,2,tt,gen)) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*(dUSf(i,1,tt,gen))) +
     & 	                  HR(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  HR(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  HR(a,i,tt,gen)*1/2d0*dconjg(dZfL(tt,gen)) 

	dFL(a,i,tt,gen) =  -EL*(1d0/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*(USf(i,2,tt+1,gen)*dMf_high(tt+1,gen) + Mf_high(tt+1,gen)*dUSf(i,2,tt+1,gen)) +
     &   			sqrt2/SW*((ef(tt+1) - I3f(tt+1))*SW/CW*ZNeuC(a,1) + I3f(tt+1)*ZNeuC(a,2))*(dUSf(i,1,tt+1,gen))) +
     & 	                  FL(a,1,tt,gen)*1/2d0*dZSf(1,i,tt+1,gen) +  FL(a,2,tt,gen)*1/2d0*dZSf(2,i,tt+1,gen) +  FL(a,i,tt,gen)*1/2d0*dZfL(tt+1,gen) 

	dFR(a,i,tt,gen) = -EL*(1d0/(SW*sqrt2*MW*CB)*ZNeu(a,3)*(USf(i,1,tt+1,gen)*dMf_high(tt+1,gen) + Mf_high(tt+1,gen)*dUSf(i,1,tt+1,gen)) -
     &   			sqrt2*ef(tt+1)/CW*ZNeu(a,1)*(dUSf(i,2,tt+1,gen))) +
     & 	                  FR(a,1,tt,gen)*1/2d0*dZSf(1,i,tt+1,gen) +  FR(a,2,tt,gen)*1/2d0*dZSf(2,i,tt+1,gen) +  FR(a,i,tt,gen)*1/2d0*dZfR(tt+1,gen)
	  
	  else
! this is not used <-> type is up-type in chargino-neutralino          
   
	dELc(a,i,tt,gen) = 0d0
	dER(a,i,tt,gen) =  0d0
	dGL(a,i,tt,gen) =  0d0
	dGR(a,i,tt,gen) =  0d0
   
	  endif
   
	  enddo
	  do a= 1,2

	    if ((tt.eq.1).or.(tt.eq.3)) then
          
	    dGL(a,i,tt,gen) =  -EL/SW*VChaC(a,1)*dUSf(i,1,tt,gen) + 
     &	                  EL/(SW*sqrt2*MW*SB)*VChaC(a,2)*(USf(i,2,tt,gen)*dMf_high(tt,gen) + Mf_high(tt,gen)*dUSf(i,2,tt,gen)) + 
     &				GL(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  GL(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  GL(a,i,tt,gen)*1/2d0*dZfL(tt+1,gen)

	    dGR(a,i,tt,gen) =  EL/(SW*sqrt2*MW*CB)*UCha(a,2)*(USf(i,1,tt,gen)*dMf_high(tt+1,gen) + Mf_high(tt+1,gen)*dUSf(i,1,tt,gen)) +
     &			     GR(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  GR(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  GR(a,i,tt,gen)*1/2d0*dZfR(tt+1,gen)
	
	    dELc(a,i,tt,gen) =  EL/(SW*sqrt2*MW*SB)*VChaC(a,2)*(USf(i,1,tt+1,gen)*dMf_high(tt,gen) + Mf_high(tt,gen)*dUSf(i,1,tt+1,gen)) + 
     &			     ELc(a,1,tt,gen)*1/2d0*dZSf(1,i,tt+1,gen) +  ELc(a,2,tt,gen)*1/2d0*dZSf(2,i,tt+1,gen) +  ELc(a,i,tt,gen)*1/2d0*dconjg(dZfR(tt,gen))

	    dER(a,i,tt,gen) = -EL/SW*UCha(a,1)*dUSf(i,1,tt+1,gen) + 
     &	                 EL/(SW*sqrt2*MW*CB)*UCha(a,2)*(USf(i,2,tt+1,gen)*dMf_high(tt+1,gen) + Mf_high(tt+1,gen)*dUSf(i,2,tt+1,gen)) +  
     &			     ER(a,1,tt,gen)*1/2d0*dZSf(1,i,tt+1,gen) +  ER(a,2,tt,gen)*1/2d0*dZSf(2,i,tt+1,gen) +  ER(a,i,tt,gen)*1/2d0*dconjg(dZfL(tt,gen))       


	    
	    else
		
! this is not used - type is up-type in chargino-neutralino          
	    FL(a,i,tt,gen) =  0d0
	    FR(a,i,tt,gen) =  0d0
	
	    HL(a,i,tt,gen) = 0d0
	    HR(a,i,tt,gen) = 0d0
          
	    endif

	  enddo
  	 enddo 

	 else if (icha.eq.2) then
	 
	 do i= 1,2
	  do a= 1,2
	 
	  if ((tt.eq.1).or.(tt.eq.3)) then

c FL=conj(ER), FR=conj(EL)
	 
	  dELc(a,i,tt,gen) =EL/(SW*sqrt2*MW*SB)*VChaC(a,2)*(Mf_high(tt,gen)*dUSf(i,1,tt+1,gen) + dMf_high(tt,gen)*USf(i,1,tt+1,gen))  +
     &			    ELc(a,1,tt,gen)*1/2d0*dZSf(1,i,tt+1,gen) +  ELc(a,2,tt,gen)*1/2d0*dZSf(2,i,tt+1,gen) +  ELc(a,i,tt,gen)*1/2d0*dconjg(dZfR(tt,gen))

	  dER(a,i,tt,gen) = -EL/SW*UCha(a,1)*dUSf(i,1,tt+1,gen) + EL/(SW*sqrt2*MW*CB)*UCha(a,2)*(dMf_high(tt+1,gen)*USf(i,2,tt+1,gen) + Mf_high(tt+1,gen)*dUSf(i,2,tt+1,gen)) +
     &	                    ER(a,1,tt,gen)*1/2d0*dZSf(1,i,tt+1,gen) +  ER(a,2,tt,gen)*1/2d0*dZSf(2,i,tt+1,gen) +  ER(a,i,tt,gen)*1/2d0*dconjg(dZfL(tt,gen))

	  dFL(a,i,tt,gen) =  -EL/SW*UChaC(a,1)*dUSf(i,1,tt+1,gen) + EL/(SW*sqrt2*MW*CB)*UChaC(a,2)*(dMf_high(tt+1,gen)*USf(i,2,tt+1,gen) + Mf_high(tt+1,gen)*dUSf(i,2,tt+1,gen)) + 
     &			    FL(a,1,tt,gen)*1/2d0*dZSf(1,i,tt+1,gen) +  FL(a,2,tt,gen)*1/2d0*dZSf(2,i,tt+1,gen) +  FL(a,i,tt,gen)*1/2d0*dZfL(tt,gen)

	  dFR(a,i,tt,gen) = EL/(SW*sqrt2*MW*SB)*VCha(a,2)*(dUSf(i,1,tt+1,gen)*Mf_high(tt,gen) + USf(i,1,tt+1,gen)*dMf_high(tt,gen)) +  
     &			    FR(a,1,tt,gen)*1/2d0*dZSf(1,i,tt+1,gen) +  FR(a,2,tt,gen)*1/2d0*dZSf(2,i,tt+1,gen) +  FR(a,i,tt,gen)*1/2d0*dZfR(tt,gen)
	  
	  dGL(a,i,tt,gen) =  0d0
	  dGR(a,i,tt,gen) =  0d0
	  dHL(a,i,tt,gen) =  0d0
	  dHR(a,i,tt,gen) =  0d0
	  
	  else
	 
	  dELc(a,i,tt,gen) = 0d0
	  dER(a,i,tt,gen) =  0d0
	  dFL(a,i,tt,gen) =  0d0
	  dFR(a,i,tt,gen) =  0d0

c HL=conj(GR), HR=conj(GL)
	 
	  dGL(a,i,tt,gen) =  -EL/SW*VChaC(a,1)*dUSf(i,1,tt-1,gen) + 
     &		  	    EL/(SW*sqrt2*MW*SB)*VChaC(a,2)*(USf(i,2,tt-1,gen)*dMf_high(tt-1,gen) + Mf_high(tt-1,gen)*dUSf(i,2,tt-1,gen)) +
     &			    GL(a,1,tt,gen)*1/2d0*dZSf(1,i,tt-1,gen) +  GL(a,2,tt,gen)*1/2d0*dZSf(2,i,tt-1,gen) +  GL(a,i,tt,gen)*1/2d0*dZfL(tt,gen)

	  dGR(a,i,tt,gen) =  EL/(SW*sqrt2*MW*CB)*UCha(a,2)*(USf(i,1,tt-1,gen)*dMf_high(tt,gen) + Mf_high(tt,gen)*dUSf(i,1,tt-1,gen)) +
     &	                GR(a,1,tt,gen)*1/2d0*dZSf(1,i,tt-1,gen) +  GR(a,2,tt,gen)*1/2d0*dZSf(2,i,tt-1,gen) +  GR(a,i,tt,gen)*1/2d0*dZfR(tt,gen)
 
	  dHL(a,i,tt,gen) =  EL/(SW*sqrt2*MW*CB)*UChaC(a,2)*(USf(i,1,tt-1,gen)*dMf_high(tt,gen) + Mf_high(tt,gen)*dUSf(i,1,tt-1,gen)) + 
     &			    HL(a,1,tt,gen)*1/2d0*dZSf(1,i,tt-1,gen) +  HL(a,2,tt,gen)*1/2d0*dZSf(2,i,tt-1,gen) +  HL(a,i,tt,gen)*1/2d0*dconjg(dZfR(tt,gen))
  
	  dHR(a,i,tt,gen) = -EL/SW*VCha(a,1)*dUSf(i,1,tt-1,gen) + 
     &			   EL/(SW*sqrt2*MW*SB)*VCha(a,2)*(USf(i,2,tt-1,gen)*dMf_high(tt-1,gen) + Mf_high(tt-1,gen)*dUSf(i,2,tt-1,gen)) +  
     &			    HR(a,1,tt,gen)*1/2d0*dZSf(1,i,tt-1,gen) +  HR(a,2,tt,gen)*1/2d0*dZSf(2,i,tt-1,gen) +  HR(a,i,tt,gen)*1/2d0*dconjg(dZfL(tt,gen)) 
	  
	  endif
	 
	  enddo
	 enddo 

	endif
	
	
	
! Counterterms for EFT vertex diagrams

	do a=1,4
	 do i =1,2
	
! This guy is constructed analogous to the H-q-q Counterterm
			dEFT_SI_L(a,i,tt,gen) = 1/2d0*dZfL(tt,gen)+ 1/2d0*dconjg(dZfR(tt,gen)) + dMf_high(tt,gen)/Mf_high(tt,gen) 
			dEFT_SI_R(a,i,tt,gen) = 1/2d0*dZfR(tt,gen)+ 1/2d0*dconjg(dZfL(tt,gen)) + dMf_high(tt,gen)/Mf_high(tt,gen)
			
! This guy is constructed analogous to the Z-q-q Counterterm
			dEFT_SD_L(a,i,tt,gen) = 1/2d0*dZfL(tt,gen) + 1/2d0*dconjg(dZfL(tt,gen)) - 1/2d0*dZfL_SUSY(tt,gen) - 1/2d0*dconjg(dZfL_SUSY(tt,gen))
			dEFT_SD_R(a,i,tt,gen) = 1/2d0*dZfR(tt,gen) + 1/2d0*dconjg(dZfR(tt,gen)) - 1/2d0*dZfL_SUSY(tt,gen) - 1/2d0*dconjg(dZfL_SUSY(tt,gen))
			
! Add an additional finite contribution to for the axial anomaly, see Larin paper (9302240v1) eq.(33). The sign is opposite as the counterterms are
! defined in agreement with Buras and Hill (9806471 and 1409.8290), the Z here is Z^-1 in Larin. If we perform IR checks, we deactivate the constant
! term.
			if(EpsPole.ge.0d0) then
				dEFT_SD_L(a,i,tt,gen) = dEFT_SD_L(a,i,tt,gen) + 4d0*CF*Alfas/(4d0*pi)
				dEFT_SD_R(a,i,tt,gen) = dEFT_SD_R(a,i,tt,gen) + 4d0*CF*Alfas/(4d0*pi)
			endif
	 enddo
	enddo
		

! End of loops over tt and gen
	enddo
	enddo
	
	
	end
		
! =======================================
! QCD Self-energy of the quark (SUSY contributions only)
! =======================================

	subroutine DMNLO_CalcPif_SUSY(myP2,tt,gen)

#include "DD_Model.h"
#include "DD_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer tt,gen,m
	double precision myP2,ML0,ML1
	double complex g1L,g1R,g2L,g2R,g1,g2

	PifVR_SUSY(tt,gen) = 0d0
	PifVL_SUSY(tt,gen) = 0d0
	PifSR_SUSY(tt,gen) = 0d0
	PifSL_SUSY(tt,gen) = 0d0

! ===== Gluino Loop =====

	do m=1,2

	 	g1L = -sqrt2*GS*USf(m,1,tt,gen)
	 	g1R = sqrt2*GS*USf(m,2,tt,gen)
	 	g2L = sqrt2*GS*USf(m,2,tt,gen)
	 	g2R = -sqrt2*GS*USf(m,1,tt,gen)

	 	ML0 = MSf(m,tt,gen) 
	 	ML1 = MGl 

	 	call myBget(myP2,ML0**2,ML1**2)
         
		PifVR_SUSY(tt,gen) = PifVR_SUSY(tt,gen) + CF/(16d0*pi**2)*g1R*g2L*(xb0 + xb1)
		PifVL_SUSY(tt,gen) = PifVL_SUSY(tt,gen) + CF/(16d0*pi**2)*g1L*g2R*(xb0 + xb1)
		PifSR_SUSY(tt,gen) = PifSR_SUSY(tt,gen) + CF/(16d0*pi**2)*g1R*g2R*ML1*xb0
		PifSL_SUSY(tt,gen) = PifSL_SUSY(tt,gen) + CF/(16d0*pi**2)*g1L*g2L*ML1*xb0

	enddo

	end	


! ===================================================
! QCD derivative of Self-energy of the quark (SUSY contributions only)  
! ===================================================

	subroutine DMNLO_CalcdPif_SUSY(myP2,tt,gen)

#include "DD_Model.h"
#include "DD_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer tt,gen,m
	double precision myP2,ML0,ML1
	double complex g1L,g1R,g2L,g2R,g1,g2


	dPifVR_SUSY(tt,gen) = 0d0
	dPifVL_SUSY(tt,gen) = 0d0
	dPifSR_SUSY(tt,gen) = 0d0
	dPifSL_SUSY(tt,gen) = 0d0

! ===== Gluino Loop =====

	do m=1,2

	 	g1L = -sqrt2*GS*USf(m,1,tt,gen)
	 	g1R = sqrt2*GS*USf(m,2,tt,gen)
	 	g2L = sqrt2*GS*USf(m,2,tt,gen)
	 	g2R = -sqrt2*GS*USf(m,1,tt,gen)

	 	ML0 = MSf(m,tt,gen) 
	 	ML1 = MGl

	 	call myBPget(myP2,ML1**2,ML0**2)

		dPifVR_SUSY(tt,gen) = dPifVR_SUSY(tt,gen) - CF/(16d0*pi**2)*g1R*g2L*xbp1
		dPifVL_SUSY(tt,gen) = dPifVL_SUSY(tt,gen) - CF/(16d0*pi**2)*g1L*g2R*xbp1
		dPifSR_SUSY(tt,gen) = dPifSR_SUSY(tt,gen) + CF/(16d0*pi**2)*g1R*g2R*ML1*xbp0
		dPifSL_SUSY(tt,gen) = dPifSL_SUSY(tt,gen) + CF/(16d0*pi**2)*g1L*g2L*ML1*xbp0

	enddo

	end	
	

