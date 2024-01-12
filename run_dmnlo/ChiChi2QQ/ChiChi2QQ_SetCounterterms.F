! =====================================================================
	subroutine ChiChi2QQ_SetCounterTerms()
	implicit none

#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_GenCouplings.h"
#include "ChiChi2QQ_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"


	integer a,i,j,k,l,m,t,g
	double precision mbmass,mtmass, mtit
        double complex DLred(6),DRred(6)

! functions from DMNLO_Kinematics.F
	double precision Krondelta, ef, I3f, hf

 !mixing angle renormalization constants

	do t=3,4                                 !Noch von mir eingebaut, vlt sogar überflüssig!?!?!?!?!?!?!?!?!
	  do g=1,3
	    do k=1,2
	      do l=1,2
                
               if(Mf(t,g).eq.0d0) then

	         dUSf(k, l, t, g) = (0D0,0D0)

               endif

	      enddo
	    enddo
	  enddo
	enddo  

!************************
! Coupling counterterms *
!************************

        if(icha.eq.1) then

! W+-f-f

	dBL = BL*(1/2d0*dZfL(tt+1,gen) + 1/2d0*dconjg(dZfL(tt,gen)))
	dBR = BR*(1/2d0*dZfR(tt+1,gen) + 1/2d0*dconjg(dZfR(tt,gen)))

	dBpL = (0d0,0d0)
	dBpR = (0d0,0d0)

        else

! Z-f-f

	dBL = BL*(1/2d0*dZfL(tt,gen) + 1/2d0*dconjg(dZfL(tt,gen)))
	dBR = BR*(1/2d0*dZfR(tt,gen) + 1/2d0*dconjg(dZfR(tt,gen)))

	dBpL = BpL*(1/2d0*dZfL(tt,gen) + 1/2d0*dconjg(dZfL(tt,gen)))
	dBpR = BpL*(1/2d0*dZfR(tt,gen) + 1/2d0*dconjg(dZfR(tt,gen)))

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
	dDL(1) = DLred(1)*(dMf(tt,gen) + 1/2d0*dZfL(tt,gen)*Mf(tt,gen) + 1/2d0*dconjg(dZfR(tt,gen))*Mf(tt,gen))
	dDR(1) = DRred(1)*(dMf(tt,gen) + 1/2d0*dZfR(tt,gen)*Mf(tt,gen) + 1/2d0*dconjg(dZfL(tt,gen))*Mf(tt,gen))
	dDL(2) = DLred(2)*(dMf(tt,gen) + 1/2d0*dZfL(tt,gen)*Mf(tt,gen) + 1/2d0*dconjg(dZfR(tt,gen))*Mf(tt,gen))
	dDR(2) = DRred(2)*(dMf(tt,gen) + 1/2d0*dZfR(tt,gen)*Mf(tt,gen) + 1/2d0*dconjg(dZfL(tt,gen))*Mf(tt,gen))
	dDL(3) = DLred(3)*(dMf(tt,gen) + 1/2d0*dZfL(tt,gen)*Mf(tt,gen) + 1/2d0*dconjg(dZfR(tt,gen))*Mf(tt,gen))
	dDR(3) = DRred(3)*(dMf(tt,gen) + 1/2d0*dZfR(tt,gen)*Mf(tt,gen) + 1/2d0*dconjg(dZfL(tt,gen))*Mf(tt,gen))
	dDL(4) = DLred(4)*(dMf(tt,gen) + 1/2d0*dZfL(tt,gen)*Mf(tt,gen) + 1/2d0*dconjg(dZfR(tt,gen))*Mf(tt,gen))
	dDR(4) = DRred(4)*(dMf(tt,gen) + 1/2d0*dZfR(tt,gen)*Mf(tt,gen) + 1/2d0*dconjg(dZfL(tt,gen))*Mf(tt,gen)) 

	if (icha.eq.1) then
		
	dDL(5) = DLred(5)*(dMf(tt,gen) + 1/2d0*dZfL(tt+1,gen)*Mf(tt,gen) + 1/2d0*dconjg(dZfR(tt,gen))*Mf(tt,gen))
	dDR(5) = DRred(5)*(dMf(tt+1,gen) + 1/2d0*dZfR(tt+1,gen)*Mf(tt+1,gen) + 1/2d0*dconjg(dZfL(tt,gen))*Mf(tt+1,gen)) 
	dDL(6) = DLred(6)*(dMf(tt,gen) + 1/2d0*dZfL(tt+1,gen)*Mf(tt,gen) + 1/2d0*dconjg(dZfR(tt,gen))*Mf(tt,gen))
	dDR(6) = DRred(6)*(dMf(tt+1,gen) + 1/2d0*dZfR(tt+1,gen)*Mf(tt+1,gen) + 1/2d0*dconjg(dZfL(tt,gen))*Mf(tt+1,gen)) 

	endif

! Neu-Sf-f

	if (icha.eq.0) then

	do a=1,4
	 do i =1,2

	if ((tt.eq.1).or.(tt.eq.3)) then

	dELc(a,i,tt,gen) = -EL*(1d0/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*(dMf(tt,gen)*USf(i,1,tt,gen) + Mf(tt,gen)*dUSf(i,1,tt,gen)) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*(dUSf(i,2,tt,gen))) +
     & 	                 ELc(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) + ELc(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) + ELc(a,i,tt,gen)*1/2d0*dconjg(dZfR(tt,gen))

	dER(a,i,tt,gen) =  -EL*(1d0/(SW*sqrt2*MW*SB)*ZNeu(a,4)*(dMf(tt,gen)*USf(i,2,tt,gen) + Mf(tt,gen)*dUSf(i,2,tt,gen)) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*(dUSf(i,1,tt,gen))) +
     & 	                  ER(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  ER(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  ER(a,i,tt,gen)*1/2d0*dconjg(dZfL(tt,gen)) 

	dFL(a,i,tt,gen) =  -EL*(1d0/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*(dMf(tt,gen)*USf(i,2,tt,gen) + Mf(tt,gen)*dUSf(i,2,tt,gen)) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeuC(a,1) + I3f(tt)*ZNeuC(a,2))*(dUSf(i,1,tt,gen))) +
     & 	                  FL(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  FL(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  FL(a,i,tt,gen)*1/2d0*dZfL(tt,gen) 

	dFR(a,i,tt,gen) = -EL*(1d0/(SW*sqrt2*MW*SB)*ZNeu(a,4)*(dMf(tt,gen)*USf(i,1,tt,gen) + Mf(tt,gen)*dUSf(i,1,tt,gen)) -
     &   			sqrt2*ef(tt)/CW*ZNeu(a,1)*(dUSf(i,2,tt,gen))) +
     & 	                  FR(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  FR(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  FR(a,i,tt,gen)*1/2d0*dZfR(tt,gen)
		
	else

	dELc(a,i,tt,gen) = -EL*(1d0/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*(dMf(tt,gen)*USf(i,1,tt,gen) + Mf(tt,gen)*dUSf(i,1,tt,gen)) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*(dUSf(i,2,tt,gen))) +
     & 	                 ELc(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) + ELc(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) + ELc(a,i,tt,gen)*1/2d0*dconjg(dZfR(tt,gen))

	dER(a,i,tt,gen) =  -EL*(1d0/(SW*sqrt2*MW*CB)*ZNeu(a,3)*(dMf(tt,gen)*USf(i,2,tt,gen) + Mf(tt,gen)*dUSf(i,2,tt,gen)) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*(dUSf(i,1,tt,gen))) +
     & 	                  ER(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  ER(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  ER(a,i,tt,gen)*1/2d0*dconjg(dZfL(tt,gen)) 

	dFL(a,i,tt,gen) =  -EL*(1d0/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*(dMf(tt,gen)*USf(i,2,tt,gen) + Mf(tt,gen)*dUSf(i,2,tt,gen)) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeuC(a,1) + I3f(tt)*ZNeuC(a,2))*(dUSf(i,1,tt,gen))) +
     & 	                  FL(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  FL(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  FL(a,i,tt,gen)*1/2d0*dZfL(tt,gen) 

	dFR(a,i,tt,gen) = -EL*(1d0/(SW*sqrt2*MW*CB)*ZNeu(a,3)*(dMf(tt,gen)*USf(i,1,tt,gen) + Mf(tt,gen)*dUSf(i,1,tt,gen)) -
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
   
	dHL(a,i,tt,gen) = -EL*(1d0/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*(USf(i,1,tt,gen)*dMf(tt,gen) + Mf(tt,gen)*dUSf(i,1,tt,gen)) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*(dUSf(i,2,tt,gen))) +
     & 	                 HL(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) + HL(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) + HL(a,i,tt,gen)*1/2d0*dconjg(dZfR(tt,gen))

	dHR(a,i,tt,gen) =  -EL*(1d0/(SW*sqrt2*MW*SB)*ZNeu(a,4)*(USf(i,2,tt,gen)*dMf(tt,gen) + Mf(tt,gen)*dUSf(i,2,tt,gen)) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*(dUSf(i,1,tt,gen))) +
     & 	                  HR(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  HR(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  HR(a,i,tt,gen)*1/2d0*dconjg(dZfL(tt,gen)) 

	dFL(a,i,tt,gen) =  -EL*(1d0/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*(USf(i,2,tt+1,gen)*dMf(tt+1,gen) + Mf(tt+1,gen)*dUSf(i,2,tt+1,gen)) +
     &   			sqrt2/SW*((ef(tt+1) - I3f(tt+1))*SW/CW*ZNeuC(a,1) + I3f(tt+1)*ZNeuC(a,2))*(dUSf(i,1,tt+1,gen))) +
     & 	                  FL(a,1,tt,gen)*1/2d0*dZSf(1,i,tt+1,gen) +  FL(a,2,tt,gen)*1/2d0*dZSf(2,i,tt+1,gen) +  FL(a,i,tt,gen)*1/2d0*dZfL(tt+1,gen) 

	dFR(a,i,tt,gen) = -EL*(1d0/(SW*sqrt2*MW*CB)*ZNeu(a,3)*(USf(i,1,tt+1,gen)*dMf(tt+1,gen) + Mf(tt+1,gen)*dUSf(i,1,tt+1,gen)) -
     &   			sqrt2*ef(tt+1)/CW*ZNeu(a,1)*(dUSf(i,2,tt+1,gen))) +
     & 	                  FR(a,1,tt,gen)*1/2d0*dZSf(1,i,tt+1,gen) +  FR(a,2,tt,gen)*1/2d0*dZSf(2,i,tt+1,gen) +  FR(a,i,tt,gen)*1/2d0*dZfR(tt+1,gen)

!          print*,'Counter-Check2',a,i,tt,gen,dHL(a,i,tt,gen),dHR(a,i,tt,gen),dFL(a,i,tt,gen),dFR(a,i,tt,gen) 
	  
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
     &	                  EL/(SW*sqrt2*MW*SB)*VChaC(a,2)*(USf(i,2,tt,gen)*dMf(tt,gen) + Mf(tt,gen)*dUSf(i,2,tt,gen)) + 
     &				GL(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  GL(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  GL(a,i,tt,gen)*1/2d0*dZfL(tt+1,gen)

	    dGR(a,i,tt,gen) =  EL/(SW*sqrt2*MW*CB)*UCha(a,2)*(USf(i,1,tt,gen)*dMf(tt+1,gen) + Mf(tt+1,gen)*dUSf(i,1,tt,gen)) +
     &			     GR(a,1,tt,gen)*1/2d0*dZSf(1,i,tt,gen) +  GR(a,2,tt,gen)*1/2d0*dZSf(2,i,tt,gen) +  GR(a,i,tt,gen)*1/2d0*dZfR(tt+1,gen)
	
	    dELc(a,i,tt,gen) =  EL/(SW*sqrt2*MW*SB)*VChaC(a,2)*(USf(i,1,tt+1,gen)*dMf(tt,gen) + Mf(tt,gen)*dUSf(i,1,tt+1,gen)) + 
     &			     ELc(a,1,tt,gen)*1/2d0*dZSf(1,i,tt+1,gen) +  ELc(a,2,tt,gen)*1/2d0*dZSf(2,i,tt+1,gen) +  ELc(a,i,tt,gen)*1/2d0*dconjg(dZfR(tt,gen))

	    dER(a,i,tt,gen) = -EL/SW*UCha(a,1)*dUSf(i,1,tt+1,gen) + 
     &	                 EL/(SW*sqrt2*MW*CB)*UCha(a,2)*(USf(i,2,tt+1,gen)*dMf(tt+1,gen) + Mf(tt+1,gen)*dUSf(i,2,tt+1,gen)) +  
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
	 
	  dELc(a,i,tt,gen) =EL/(SW*sqrt2*MW*SB)*VChaC(a,2)*(Mf(tt,gen)*dUSf(i,1,tt+1,gen) + dMf(tt,gen)*USf(i,1,tt+1,gen))  +
     &			    ELc(a,1,tt,gen)*1/2d0*dZSf(1,i,tt+1,gen) +  ELc(a,2,tt,gen)*1/2d0*dZSf(2,i,tt+1,gen) +  ELc(a,i,tt,gen)*1/2d0*dconjg(dZfR(tt,gen))

	  dER(a,i,tt,gen) = -EL/SW*UCha(a,1)*dUSf(i,1,tt+1,gen) + EL/(SW*sqrt2*MW*CB)*UCha(a,2)*(dMf(tt+1,gen)*USf(i,2,tt+1,gen) + Mf(tt+1,gen)*dUSf(i,2,tt+1,gen)) +
     &	                    ER(a,1,tt,gen)*1/2d0*dZSf(1,i,tt+1,gen) +  ER(a,2,tt,gen)*1/2d0*dZSf(2,i,tt+1,gen) +  ER(a,i,tt,gen)*1/2d0*dconjg(dZfL(tt,gen))

	  dFL(a,i,tt,gen) =  -EL/SW*UChaC(a,1)*dUSf(i,1,tt+1,gen) + EL/(SW*sqrt2*MW*CB)*UChaC(a,2)*(dMf(tt+1,gen)*USf(i,2,tt+1,gen) + Mf(tt+1,gen)*dUSf(i,2,tt+1,gen)) + 
     &			    FL(a,1,tt,gen)*1/2d0*dZSf(1,i,tt+1,gen) +  FL(a,2,tt,gen)*1/2d0*dZSf(2,i,tt+1,gen) +  FL(a,i,tt,gen)*1/2d0*dZfL(tt,gen)

	  dFR(a,i,tt,gen) = EL/(SW*sqrt2*MW*SB)*VCha(a,2)*(dUSf(i,1,tt+1,gen)*Mf(tt,gen) + USf(i,1,tt+1,gen)*dMf(tt,gen)) +  
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
     &		  	    EL/(SW*sqrt2*MW*SB)*VChaC(a,2)*(USf(i,2,tt-1,gen)*dMf(tt-1,gen) + Mf(tt-1,gen)*dUSf(i,2,tt-1,gen)) +
     &			    GL(a,1,tt,gen)*1/2d0*dZSf(1,i,tt-1,gen) +  GL(a,2,tt,gen)*1/2d0*dZSf(2,i,tt-1,gen) +  GL(a,i,tt,gen)*1/2d0*dZfL(tt,gen)

	  dGR(a,i,tt,gen) =  EL/(SW*sqrt2*MW*CB)*UCha(a,2)*(USf(i,1,tt-1,gen)*dMf(tt,gen) + Mf(tt,gen)*dUSf(i,1,tt-1,gen)) +
     &	                GR(a,1,tt,gen)*1/2d0*dZSf(1,i,tt-1,gen) +  GR(a,2,tt,gen)*1/2d0*dZSf(2,i,tt-1,gen) +  GR(a,i,tt,gen)*1/2d0*dZfR(tt,gen)
 
	  dHL(a,i,tt,gen) =  EL/(SW*sqrt2*MW*CB)*UChaC(a,2)*(USf(i,1,tt-1,gen)*dMf(tt,gen) + Mf(tt,gen)*dUSf(i,1,tt-1,gen)) + 
     &			    HL(a,1,tt,gen)*1/2d0*dZSf(1,i,tt-1,gen) +  HL(a,2,tt,gen)*1/2d0*dZSf(2,i,tt-1,gen) +  HL(a,i,tt,gen)*1/2d0*dconjg(dZfR(tt,gen))
  
	  dHR(a,i,tt,gen) = -EL/SW*VCha(a,1)*dUSf(i,1,tt-1,gen) + 
     &			   EL/(SW*sqrt2*MW*SB)*VCha(a,2)*(USf(i,2,tt-1,gen)*dMf(tt-1,gen) + Mf(tt-1,gen)*dUSf(i,2,tt-1,gen)) +  
     &			    HR(a,1,tt,gen)*1/2d0*dZSf(1,i,tt-1,gen) +  HR(a,2,tt,gen)*1/2d0*dZSf(2,i,tt-1,gen) +  HR(a,i,tt,gen)*1/2d0*dconjg(dZfL(tt,gen)) 
	  
	  endif
	 
	  enddo
	 enddo 

	endif

c	print*,'dEL1,dER1=',dELc(1,1,tt,gen),dER(1,1,tt,gen)
c	print*,'dEL2,dER2=',dELc(1,2,tt,gen),dER(1,2,tt,gen)

	end

!**************************
! Propagator counterterms *
!**************************

	subroutine ChiChi2QQ_CalcProp(myP2,t,g)

#include "ChiChi2QQ_GenCounterterms.h"
#include "ChiChi2QQ_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer t,g,i,j
	double precision myP2, Krondelta

	 do i = 1, 2
	  do j = 1, 2

	  Prop(i,j,t,g) = 1/2d0*dZSf(i,j,t,g)*(myP2 - MSf2(i,t,g)) + 1/2d0*dZSf(j,i,t,g)*(myP2 - MSf2(j,t,g)) - Krondelta(i,j)*dMSf2(i,t,g)

	  enddo
 	 enddo

	end

	

