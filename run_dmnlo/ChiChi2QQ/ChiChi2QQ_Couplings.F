	subroutine ChiChi2QQ_SetCouplings
	implicit none

#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_GenCouplings.h"
#include "ChiChi2QQ_GenCounterterms.h"

	double complex GLR(6,2,2,4,3),ZLR(2,2,4,3)
	integer i,j,k,l,a  
	double precision DmbResults(5), DeltaTlighth, DeltaTheavyH, DeltaTA
	
! functions from DMNLO_Kinematics.F
	double precision Krondelta, ef, I3f, hf
		
	external DMNLO_CalcDeltaMb

	
! ******************************************************
! Gaugino-Vector boson & Fermion-Vector boson couplings
! ******************************************************
	if (icha.eq.2) then
	
	ApL = -EL*Krondelta(indi,indj)	
	BpL = -EL*ef(tt)

	AL = EL/(SW*CW)*(-VCha(indj,1)*VChaC(indi,1) - 1/2d0*VCha(indj,2)*VChaC(indi,2) + Krondelta(indi,indj)*SW**2)
	AR = EL/(SW*CW)*(-UChaC(indj,1)*UCha(indi,1) - 1/2d0*UChaC(indj,2)*UCha(indi,2) + Krondelta(indi,indj)*SW**2)
	BL = EL/(SW*CW)*(ef(tt)*SW**2-I3f(tt))
	BR = EL/(SW*CW)*(ef(tt)*SW**2)
	
	else if (icha.eq.1) then

	AL = EL/SW*( ZNeu(indj,2)*VChaC(indi,1) - 1/sqrt2*ZNeu(indj,4)*VChaC(indi,2))
	AR = EL/SW*( ZNeuC(indj,2)*UCha(indi,1) + 1/sqrt2*ZNeuC(indj,3)*UCha(indi,2))
	BL = -EL/(SW*sqrt2)
	BR = 0d0

	else if (icha.eq.0) then

	AL = EL/(2*SW*CW)*(- ZNeu(indj,3)*ZNeuC(indi,3) + ZNeu(indj,4)*ZNeuC(indi,4))
	AR = EL/(2*SW*CW)*(  ZNeuC(indj,3)*ZNeu(indi,3) - ZNeuC(indj,4)*ZNeu(indi,4))
	BL = EL/(SW*CW)*(ef(tt)*SW**2-I3f(tt))
	BR = EL/(SW*CW)*(ef(tt)*SW**2)

	endif
	
! ************************	
! Gaugino-Higgs couplings
! ************************
! Higgs index Hk0 = (h0,H0,A0,G0) = (1,2,3,4) or Hk = (h0,H0,A0,G0,Hp,Gp) = (1,2,3,4,5,6)

	if (icha.eq.0) then
! NeuNeuHiggs
	CL(1) = -EL/(2d0*SW)*((-SA)*(ZNeuC(indi,3)*ZNeuC(indj,2) + ZNeuC(indj,3)*ZNeuC(indi,2) - SW/CW*(ZNeuC(indi,3)*ZNeuC(indj,1) + ZNeuC(indj,3)*ZNeuC(indi,1)))
     &	        + (-CA)*(ZNeuC(indi,4)*ZNeuC(indj,2) + ZNeuC(indj,4)*ZNeuC(indi,2) - SW/CW*(ZNeuC(indi,4)*ZNeuC(indj,1) + ZNeuC(indj,4)*ZNeuC(indi,1))) )
	CR(1) = -EL/(2d0*SW)*((-SA)*( ZNeu(indi,3)*ZNeu(indj,2)  +  ZNeu(indj,3)*ZNeu(indi,2)  - SW/CW*( ZNeu(indi,3)*ZNeu(indj,1)  +  ZNeu(indj,3)*ZNeu(indi,1)))
     &	        + (-CA)*( ZNeu(indi,4)*ZNeu(indj,2)  +  ZNeu(indj,4)*ZNeu(indi,2)  - SW/CW*( ZNeu(indi,4)*ZNeu(indj,1)  +  ZNeu(indj,4)*ZNeu(indi,1))) )
	CL(2) = -EL/(2d0*SW)*((+CA)*(ZNeuC(indi,3)*ZNeuC(indj,2) + ZNeuC(indj,3)*ZNeuC(indi,2) - SW/CW*(ZNeuC(indi,3)*ZNeuC(indj,1) + ZNeuC(indj,3)*ZNeuC(indi,1)))
     &	        + (-SA)*(ZNeuC(indi,4)*ZNeuC(indj,2) + ZNeuC(indj,4)*ZNeuC(indi,2) - SW/CW*(ZNeuC(indi,4)*ZNeuC(indj,1) + ZNeuC(indj,4)*ZNeuC(indi,1))) )
	CR(2) = -EL/(2d0*SW)*((+CA)*( ZNeu(indi,3)*ZNeu(indj,2)  +  ZNeu(indj,3)*ZNeu(indi,2)  - SW/CW*( ZNeu(indi,3)*ZNeu(indj,1)  +  ZNeu(indj,3)*ZNeu(indi,1)))
     &	        + (-SA)*( ZNeu(indi,4)*ZNeu(indj,2)  +  ZNeu(indj,4)*ZNeu(indi,2)  - SW/CW*( ZNeu(indi,4)*ZNeu(indj,1)  +  ZNeu(indj,4)*ZNeu(indi,1))) )
	CL(3) = -cI*EL/(2d0*SW)*((-SB)*(ZNeuC(indi,3)*ZNeuC(indj,2) + ZNeuC(indj,3)*ZNeuC(indi,2) - SW/CW*(ZNeuC(indi,3)*ZNeuC(indj,1) + ZNeuC(indj,3)*ZNeuC(indi,1)))
     &	           + (+CB)*(ZNeuC(indi,4)*ZNeuC(indj,2) + ZNeuC(indj,4)*ZNeuC(indi,2) - SW/CW*(ZNeuC(indi,4)*ZNeuC(indj,1) + ZNeuC(indj,4)*ZNeuC(indi,1))) )
	CR(3) = +cI*EL/(2d0*SW)*((-SB)*( ZNeu(indi,3)*ZNeu(indj,2)  +  ZNeu(indj,3)*ZNeu(indi,2)  - SW/CW*( ZNeu(indi,3)*ZNeu(indj,1)  +  ZNeu(indj,3)*ZNeu(indi,1)))
     &	           + (+CB)*( ZNeu(indi,4)*ZNeu(indj,2)  +  ZNeu(indj,4)*ZNeu(indi,2)  - SW/CW*( ZNeu(indi,4)*ZNeu(indj,1)  +  ZNeu(indj,4)*ZNeu(indi,1))) )
	CL(4) = -cI*EL/(2d0*SW)*((+CB)*(ZNeuC(indi,3)*ZNeuC(indj,2) + ZNeuC(indj,3)*ZNeuC(indi,2) - SW/CW*(ZNeuC(indi,3)*ZNeuC(indj,1) + ZNeuC(indj,3)*ZNeuC(indi,1)))
     &	           + (+SB)*(ZNeuC(indi,4)*ZNeuC(indj,2) + ZNeuC(indj,4)*ZNeuC(indi,2) - SW/CW*(ZNeuC(indi,4)*ZNeuC(indj,1) + ZNeuC(indj,4)*ZNeuC(indi,1))) )
	CR(4) = +cI*EL/(2d0*SW)*((+CB)*( ZNeu(indi,3)*ZNeu(indj,2)  +  ZNeu(indj,3)*ZNeu(indi,2)  - SW/CW*( ZNeu(indi,3)*ZNeu(indj,1)  +  ZNeu(indj,3)*ZNeu(indi,1)))
     &	           + (+SB)*( ZNeu(indi,4)*ZNeu(indj,2)  +  ZNeu(indj,4)*ZNeu(indi,2)  - SW/CW*( ZNeu(indi,4)*ZNeu(indj,1)  +  ZNeu(indj,4)*ZNeu(indi,1))) )

! NeuChaHiggs
	CL(5) = 0d0
	CR(5) = 0d0
	CL(6) = 0d0
	CR(6) = 0d0

	else if (icha.eq.1) then
! NeuNeuHiggs
	CL(1) = 0d0
	CR(1) = 0d0
	CL(2) = 0d0	
	CR(2) = 0d0
	CL(3) = 0d0
	CR(3) = 0d0  
	CL(4) = 0d0
	CR(4) = 0d0

! NeuChaHiggs
	CL(5) = -EL/SW*CB*(VChaC(indi,1)*ZNeuC(indj,4) + 1/sqrt2*(ZNeuC(indj,2) + SW/CW*ZNeuC(indj,1))*VChaC(indi,2)) 
	CR(5) = -EL/SW*SB*(UCha(indi,1)*ZNeu(indj,3) - 1/sqrt2*(ZNeu(indj,2) + SW/CW*ZNeu(indj,1))*UCha(indi,2)) 
	CL(6) = -EL/SW*SB*(VChaC(indi,1)*ZNeuC(indj,4) + 1/sqrt2*(ZNeuC(indj,2) + SW/CW*ZNeuC(indj,1))*VChaC(indi,2)) 
	CR(6) =  EL/SW*CB*(UCha(indi,1)*ZNeu(indj,3) - 1/sqrt2*(ZNeu(indj,2) + SW/CW*ZNeu(indj,1))*UCha(indi,2))  
	
	else if (icha.eq.2) then	
! ChaChaHiggs
	CL(1) = -EL/(SW*sqrt2)*(-SA*VChaC(indi,1)*UChaC(indj,2) + CA*VChaC(indi,2)*UChaC(indj,1))
	CR(1) = -EL/(SW*sqrt2)*(-SA*VCha(indj,1)*UCha(indi,2)   + CA*VCha(indj,2)*UCha(indi,1))
	CL(2) = -EL/(SW*sqrt2)*( CA*VChaC(indi,1)*UChaC(indj,2) + SA*VChaC(indi,2)*UChaC(indj,1))
	CR(2) = -EL/(SW*sqrt2)*( CA*VCha(indj,1)*UCha(indi,2) + SA*VCha(indj,2)*UCha(indi,1))
	CL(3) = -cI*EL/(SW*sqrt2)*(-SB*VChaC(indi,1)*UChaC(indj,2) - CB*VChaC(indi,2)*UChaC(indj,1))
	CR(3) = cI*EL/(SW*sqrt2)*(-SB*VCha(indj,1)*UCha(indi,2) - CB*VCha(indj,2)*UCha(indi,1))
	CL(4) = -cI*EL/(SW*sqrt2)*( CB*VChaC(indi,1)*UChaC(indj,2) - SB*VChaC(indi,2)*UChaC(indj,1))
	CR(4) = cI*EL/(SW*sqrt2)*( CB*VCha(indj,1)*UCha(indi,2) - SB*VCha(indj,2)*UCha(indi,1))

! NeuChaHiggs
	CL(5) = 0d0
	CR(5) = 0d0
	CL(6) = 0d0
	CR(6) = 0d0

	endif
	
!**************************	
! Fermion-Higgs couplings *
!**************************

	if (icha.eq.1) then
		
! FerFerHiggs
	 if ((tt.eq.1).or.(tt.eq.3)) then

  	  DL(1) = -EL*Mf(tt,gen)*CA/(SW*2d0*MW*SB) 
  	  DR(1) = -EL*Mf(tt,gen)*CA/(SW*2d0*MW*SB) 
  	  DL(2) = -EL*Mf(tt,gen)*SA/(SW*2d0*MW*SB) 
  	  DR(2) = -EL*Mf(tt,gen)*SA/(SW*2d0*MW*SB) 
  	  DL(3) = -cI*EL*Mf(tt,gen)*CB/(SW*2d0*MW*SB) 
  	  DR(3) = cI*EL*Mf(tt,gen)*CB/(SW*2d0*MW*SB) 
  	  DL(4) = -cI*EL*Mf(tt,gen)/(SW*2d0*MW) 
  	  DR(4) = cI*EL*Mf(tt,gen)/(SW*2d0*MW) 

	 else

	  DL(1) = EL*Mf(tt,gen)*SA/(SW*2d0*MW*CB) 
	  DR(1) = EL*Mf(tt,gen)*SA/(SW*2d0*MW*CB) 
	  DL(2) = -EL*Mf(tt,gen)*CA/(SW*2d0*MW*CB) 
	  DR(2) = -EL*Mf(tt,gen)*CA/(SW*2d0*MW*CB) 
	  DL(3) = -cI*EL*Mf(tt,gen)*SB/(SW*2d0*MW*CB)
 	  DR(3) = cI*EL*Mf(tt,gen)*SB/(SW*2d0*MW*CB) 
	  DL(4) = cI*EL*Mf(tt,gen)/(SW*2d0*MW)
 	  DR(4) = -cI*EL*Mf(tt,gen)/(SW*2d0*MW) 

	 endif

! Fer1Fer2Higgs (tt = type of the up-type fermion)
	DL(5) = EL*Mf(tt,gen)*CB/(SW*sqrt2*MW*SB)
	DR(5) = EL*Mf(tt+1,gen)*SB/(SW*sqrt2*MW*CB)
	DL(6) = EL*Mf(tt,gen)/(SW*sqrt2*MW)
	DR(6) = -EL*Mf(tt+1,gen)/(SW*sqrt2*MW)

	else
! Yukawa coupling of the neutral Higgses to fermions

! FerFerHiggs
	if ((tt.eq.1).or.(tt.eq.3)) then

	DL(1) = -EL*Mf(tt,gen)*CA/(SW*2d0*MW*SB) 
	DR(1) = -EL*Mf(tt,gen)*CA/(SW*2d0*MW*SB) 
	DL(2) = -EL*Mf(tt,gen)*SA/(SW*2d0*MW*SB) 
	DR(2) = -EL*Mf(tt,gen)*SA/(SW*2d0*MW*SB) 
	DL(3) = -cI*EL*Mf(tt,gen)*CB/(SW*2d0*MW*SB) 
	DR(3) = cI*EL*Mf(tt,gen)*CB/(SW*2d0*MW*SB) 
	DL(4) = -cI*EL*Mf(tt,gen)/(SW*2d0*MW) 
	DR(4) = cI*EL*Mf(tt,gen)/(SW*2d0*MW) 	
	
	else

	DL(1) = EL*Mf(tt,gen)*SA/(SW*2d0*MW*CB) 
	DR(1) = EL*Mf(tt,gen)*SA/(SW*2d0*MW*CB) 
	DL(2) = -EL*Mf(tt,gen)*CA/(SW*2d0*MW*CB) 
	DR(2) = -EL*Mf(tt,gen)*CA/(SW*2d0*MW*CB) 
	DL(3) = -cI*EL*Mf(tt,gen)*SB/(SW*2d0*MW*CB)
 	DR(3) = cI*EL*Mf(tt,gen)*SB/(SW*2d0*MW*CB) 
	DL(4) = cI*EL*Mf(tt,gen)/(SW*2d0*MW)												
 	DR(4) = -cI*EL*Mf(tt,gen)/(SW*2d0*MW) 

	endif

! Fer1Fer2Higgs
	DL(5) = 0d0
	DR(5) = 0d0
	DL(6) = 0d0
	DR(6) = 0d0
	
	endif
	
!*************************************
! Gaugino-Fermion-Sfermion couplings *
!*************************************

! NeuFerSfer
	if (icha.eq.0) then

	 do a= 1, 4
	  do i = 1, 2
        
	  if ((tt.eq.1).or.(tt.eq.3)) then
        
	  HL(a,i,tt,gen) = -EL*(Mf(tt,gen)/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*USf(i,2,tt,gen))
	  HR(a,i,tt,gen) =  -EL*(Mf(tt,gen)/(SW*sqrt2*MW*SB)*ZNeu(a,4)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*USf(i,1,tt,gen))
	  GL(a,i,tt,gen) =  -EL*(Mf(tt,gen)/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeuC(a,1) + I3f(tt)*ZNeuC(a,2))*USf(i,1,tt,gen))
	  GR(a,i,tt,gen) =  -EL*(Mf(tt,gen)/(SW*sqrt2*MW*SB)*ZNeu(a,4)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeu(a,1)*USf(i,2,tt,gen))
        
	  FL(a,i,tt,gen) =  -EL*(Mf(tt,gen)/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeuC(a,1) + I3f(tt)*ZNeuC(a,2))*USf(i,1,tt,gen))
	  FR(a,i,tt,gen) =  -EL*(Mf(tt,gen)/(SW*sqrt2*MW*SB)*ZNeu(a,4)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeu(a,1)*USf(i,2,tt,gen))
	  ELc(a,i,tt,gen) =  -EL*(Mf(tt,gen)/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*USf(i,2,tt,gen))
	  ER(a,i,tt,gen) =  -EL*(Mf(tt,gen)/(SW*sqrt2*MW*SB)*ZNeu(a,4)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*USf(i,1,tt,gen))
	  
	  else
        
c Moritz: H und G sind u-Kanal, E und F sind t-Kanal

	  HL(a,i,tt,gen) = -EL*(Mf(tt,gen)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*USf(i,2,tt,gen))
	  HR(a,i,tt,gen) =  -EL*(Mf(tt,gen)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*USf(i,1,tt,gen))
	  GL(a,i,tt,gen) =  -EL*(Mf(tt,gen)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeuC(a,1) + I3f(tt)*ZNeuC(a,2))*USf(i,1,tt,gen))
	  GR(a,i,tt,gen) =  -EL*(Mf(tt,gen)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeu(a,1)*USf(i,2,tt,gen))
        
	  FL(a,i,tt,gen) =  -EL*(Mf(tt,gen)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeuC(a,1) + I3f(tt)*ZNeuC(a,2))*USf(i,1,tt,gen))
	  FR(a,i,tt,gen) =  -EL*(Mf(tt,gen)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeu(a,1)*USf(i,2,tt,gen))
	  ELc(a,i,tt,gen) =  -EL*(Mf(tt,gen)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*USf(i,2,tt,gen))
	  ER(a,i,tt,gen) =  -EL*(Mf(tt,gen)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*USf(i,1,tt,gen))
		
		
	  endif
        
	  enddo
       enddo 

	 else if (icha.eq.1) then
c NOTE: tt is type of outgoing quark (not antiquark), e.g. tt=top=3 for chi+ chi0 -> t bbar 
	 do i = 1, 2

! neutralino couplings
	  do a= 1, 4
   	  if ((tt.eq.1).or.(tt.eq.3)) then
   
	  HL(a,i,tt,gen) = -EL*(Mf(tt,gen)/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*USf(i,2,tt,gen))
	  HR(a,i,tt,gen) =  -EL*(Mf(tt,gen)/(SW*sqrt2*MW*SB)*ZNeu(a,4)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*USf(i,1,tt,gen))
   
	  FL(a,i,tt,gen) =  -EL*(Mf(tt+1,gen)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,2,tt+1,gen) +
     &   			sqrt2/SW*((ef(tt+1) - I3f(tt+1))*SW/CW*ZNeuC(a,1) + I3f(tt+1)*ZNeuC(a,2))*USf(i,1,tt+1,gen))
	  FR(a,i,tt,gen) =  -EL*(Mf(tt+1,gen)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,1,tt+1,gen) -
     &   			sqrt2*ef(tt+1)/CW*ZNeu(a,1)*USf(i,2,tt+1,gen))
	  
	  else
   
	  HL(a,i,tt,gen) = -EL*(Mf(tt,gen)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*USf(i,2,tt,gen))
	  HR(a,i,tt,gen) =  -EL*(Mf(tt,gen)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*USf(i,1,tt,gen))
   
	  FL(a,i,tt,gen) =  -EL*(Mf(tt-1,gen)/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*USf(i,2,tt-1,gen) +
     &   			sqrt2/SW*((ef(tt-1) - I3f(tt-1))*SW/CW*ZNeuC(a,1) + I3f(tt-1)*ZNeuC(a,2))*USf(i,1,tt-1,gen))
	  FR(a,i,tt,gen) =  -EL*(Mf(tt-1,gen)/(SW*sqrt2*MW*SB)*ZNeu(a,4)*USf(i,1,tt-1,gen) -
     &   			sqrt2*ef(tt-1)/CW*ZNeu(a,1)*USf(i,2,tt-1,gen))
   
	  endif
   	  enddo

! chargino couplings
	  do a= 1, 2

	    if ((tt.eq.1).or.(tt.eq.3)) then

	    GL(a,i,tt,gen) =  -EL/SW*VChaC(a,1)*USf(i,1,tt,gen) + EL*Mf(tt,gen)/(SW*sqrt2*MW*SB)*VChaC(a,2)*USf(i,2,tt,gen)
	    GR(a,i,tt,gen) =  EL*Mf(tt+1,gen)/(SW*sqrt2*MW*CB)*UCha(a,2)*USf(i,1,tt,gen)

	    ELc(a,i,tt,gen) = EL*Mf(tt,gen)/(SW*sqrt2*MW*SB)*VChaC(a,2)*USf(i,1,tt+1,gen) 
	    ER(a,i,tt,gen) = -EL/SW*UCha(a,1)*USf(i,1,tt+1,gen) + EL*Mf(tt+1,gen)/(SW*sqrt2*MW*CB)*UCha(a,2)*USf(i,2,tt+1,gen) 
	    
	    else

	    GL(a,i,tt,gen) =  -EL/SW*UChaC(a,1)*USf(i,1,tt,gen) + EL*Mf(tt,gen)/(SW*sqrt2*MW*CB)*UChaC(a,2)*USf(i,2,tt,gen)
	    GR(a,i,tt,gen) =  EL*Mf(tt-1,gen)/(SW*sqrt2*MW*SB)*VCha(a,2)*USf(i,1,tt,gen)
	
	    ELc(a,i,tt,gen) =  EL*Mf(tt,gen)/(SW*sqrt2*MW*CB)*UChaC(a,2)*USf(i,1,tt-1,gen)  
	    ER(a,i,tt,gen) = -EL/SW*VCha(a,1)*USf(i,1,tt-1,gen) + EL*Mf(tt-1,gen)/(SW*sqrt2*MW*SB)*VCha(a,2)*USf(i,2,tt-1,gen) 
          
	    endif

	  enddo

  	 enddo 

! cha-cha case
	 else if (icha.eq.2) then

	 do i = 1, 2
	  do a= 1, 2

	  if ((tt.eq.1).or.(tt.eq.3)) then

c H=G=0 to "switch off" u-channel
	  HL(a,i,tt,gen) =  0d0
	  HR(a,i,tt,gen) =  0d0
	  GL(a,i,tt,gen) =  0d0
	  GR(a,i,tt,gen) =  0d0

c FL=conj(ER), FR=conj(EL)
	  ELc(a,i,tt,gen) =  EL*Mf(tt,gen)/(SW*sqrt2*MW*SB)*VChaC(a,2)*USf(i,1,tt+1,gen) 
	  ER(a,i,tt,gen) =  -EL/SW*UCha(a,1)*USf(i,1,tt+1,gen) + EL*Mf(tt+1,gen)/(SW*sqrt2*MW*CB)*UCha(a,2)*USf(i,2,tt+1,gen) 

	  FL(a,i,tt,gen) =  -EL/SW*UChaC(a,1)*USf(i,1,tt+1,gen) + EL*Mf(tt+1,gen)/(SW*sqrt2*MW*CB)*UChaC(a,2)*USf(i,2,tt+1,gen)
	  FR(a,i,tt,gen)  =  EL*Mf(tt,gen)/(SW*sqrt2*MW*SB)*VCha(a,2)*USf(i,1,tt+1,gen)
	  
	  else

c HL=conj(GR), HR=conj(GL)
	  GL(a,i,tt,gen)  =  -EL/SW*VChaC(a,1)*USf(i,1,tt-1,gen) + EL*Mf(tt-1,gen)/(SW*sqrt2*MW*SB)*VChaC(a,2)*USf(i,2,tt-1,gen)
	  GR(a,i,tt,gen) =  EL*Mf(tt,gen)/(SW*sqrt2*MW*CB)*UCha(a,2)*USf(i,1,tt-1,gen)

	  HL(a,i,tt,gen) =  EL*Mf(tt,gen)/(SW*sqrt2*MW*CB)*UChaC(a,2)*USf(i,1,tt-1,gen)
	  HR(a,i,tt,gen)  =  -EL/SW*VCha(a,1)*USf(i,1,tt-1,gen) + EL*Mf(tt-1,gen)/(SW*sqrt2*MW*SB)*VCha(a,2)*USf(i,2,tt-1,gen)

c E=F=0 to "switch off" t-channel	  
	  FL(a,i,tt,gen) =  0d0
	  FR(a,i,tt,gen) =  0d0
	  ELc(a,i,tt,gen) =  0d0
	  ER(a,i,tt,gen) =  0d0

	  endif

	  enddo
	 enddo 

	endif
	
!*******************************************
! Vector boson-Sfermion-Sfermion couplings *
!*******************************************

	if (icha.eq.1) then
	 do i=1,2
	  do j=1,2

! W-boson	    
	    ZSfSf(i,j,tt,gen) = -EL/(sqrt2*SW)*USf(i,1,tt,gen)*USf(j,1,tt+1,gen)
	
	  enddo
	 enddo

	else

	  ZLR(1,1,tt,gen) = -EL/(SW*CW)*(I3f(tt)-ef(tt)*SW2)
	  ZLR(2,2,tt,gen) = EL/(SW*CW)*ef(tt)*SW2
	  ZLR(1,2,tt,gen) = 0d0 
	  ZLR(2,1,tt,gen) = 0d0
        
	  do i=1,2
	   do j=1,2
        
	    ZSfSf(i,j,tt,gen) = 0d0
        
	    do k=1,2
	     do l=1,2

! Z-boson	  
	      ZSfSf(i,j,tt,gen) = ZSfSf(i,j,tt,gen) + ZLR(k,l,tt,gen)*USf(i,k,tt,gen)*USf(j,l,tt,gen)
        
   	     enddo
	    enddo
! Photon	    
	    GSfSf(i,j,tt,gen) = -EL*ef(tt)*Krondelta(i,j)
	
 	   enddo
	  enddo

	endif

!**************************
! Higgs-Sfermion-Sfermion *
!**************************

	if ((tt.eq.1).or.(tt.eq.3)) then

	GLR(1,1,1,tt,gen) = -sqrt2*hf(tt,gen)*Mf(tt,gen)*CA + EL*MZ*(I3f(tt)-ef(tt)*SW2)*SAB/(SW*CW) 
	GLR(1,2,2,tt,gen) = -sqrt2*hf(tt,gen)*Mf(tt,gen)*CA + EL*MZ*ef(tt)*SW2*SAB/(SW*CW) 
	GLR(1,1,2,tt,gen) = -hf(tt,gen)/sqrt2*(Af(tt,gen)*CA + MUE*SA) 
	GLR(1,2,1,tt,gen) = -hf(tt,gen)/sqrt2*(Af(tt,gen)*CA + MUE*SA) 

	GLR(2,1,1,tt,gen) = -sqrt2*hf(tt,gen)*Mf(tt,gen)*SA - EL*MZ*(I3f(tt)-ef(tt)*SW2)*CAB/(SW*CW) 
	GLR(2,2,2,tt,gen) = -sqrt2*hf(tt,gen)*Mf(tt,gen)*SA - EL*MZ*ef(tt)*SW2*CAB/(SW*CW) 
	GLR(2,1,2,tt,gen) = -hf(tt,gen)/sqrt2*(Af(tt,gen)*SA - MUE*CA) 
	GLR(2,2,1,tt,gen) = -hf(tt,gen)/sqrt2*(Af(tt,gen)*SA - MUE*CA) 

	GLR(3,1,1,tt,gen) = 0d0
	GLR(3,2,2,tt,gen) = 0d0
	GLR(3,1,2,tt,gen) = cI*hf(tt,gen)/sqrt2*(Af(tt,gen)*CB + MUE*SB) 
	GLR(3,2,1,tt,gen) = -cI*hf(tt,gen)/sqrt2*(Af(tt,gen)*CB + MUE*SB) 

	GLR(4,1,1,tt,gen) = 0d0
	GLR(4,2,2,tt,gen) = 0d0
	GLR(4,1,2,tt,gen) = cI*hf(tt,gen)/sqrt2*(Af(tt,gen)*SB - MUE*CB) 
	GLR(4,2,1,tt,gen) = -cI*hf(tt,gen)/sqrt2*(Af(tt,gen)*SB - MUE*CB) 

! G_{LR}^{tb} where tt = up-type
	GLR(5,1,1,tt,gen) = hf(tt+1,gen)*Mf(tt+1,gen)*SB + hf(tt,gen)*Mf(tt,gen)*CB - EL*MW/(sqrt2*SW)*S2B
	GLR(5,2,2,tt,gen) = hf(tt,gen)*Mf(tt+1,gen)*CB + hf(tt+1,gen)*Mf(tt,gen)*SB
	GLR(5,1,2,tt,gen) = hf(tt+1,gen)*(Af(tt+1,gen)*SB + MUE*CB) 
	GLR(5,2,1,tt,gen) = hf(tt,gen)*(Af(tt,gen)*CB + MUE*SB)

	GLR(6,1,1,tt,gen) = -hf(tt+1,gen)*Mf(tt+1,gen)*CB + hf(tt,gen)*Mf(tt,gen)*SB + EL*MW/(sqrt2*SW)*S2B
	GLR(6,2,2,tt,gen) = hf(tt,gen)*Mf(tt+1,gen)*SB - hf(tt+1,gen)*Mf(tt,gen)*CB
	GLR(6,1,2,tt,gen) = hf(tt+1,gen)*(-Af(tt+1,gen)*CB + MUE*SB) 
	GLR(6,2,1,tt,gen) = hf(tt,gen)*(Af(tt,gen)*SB - MUE*CB)

	else

	GLR(1,1,1,tt,gen) = sqrt2*hf(tt,gen)*Mf(tt,gen)*SA + EL*MZ*(I3f(tt)-ef(tt)*SW2)*SAB/(SW*CW) 
	GLR(1,2,2,tt,gen) = sqrt2*hf(tt,gen)*Mf(tt,gen)*SA + EL*MZ*ef(tt)*SW2*SAB/(SW*CW) 
	GLR(1,1,2,tt,gen) = hf(tt,gen)/sqrt2*(Af(tt,gen)*SA + MUE*CA) 
	GLR(1,2,1,tt,gen) = hf(tt,gen)/sqrt2*(Af(tt,gen)*SA + MUE*CA) 

	GLR(2,1,1,tt,gen) = -sqrt2*hf(tt,gen)*Mf(tt,gen)*CA - EL*MZ*(I3f(tt)-ef(tt)*SW2)*CAB/(SW*CW) 
	GLR(2,2,2,tt,gen) = -sqrt2*hf(tt,gen)*Mf(tt,gen)*CA - EL*MZ*ef(tt)*SW2*CAB/(SW*CW) 
	GLR(2,1,2,tt,gen) = -hf(tt,gen)/sqrt2*(Af(tt,gen)*CA - MUE*SA) 
	GLR(2,2,1,tt,gen) = -hf(tt,gen)/sqrt2*(Af(tt,gen)*CA - MUE*SA) 

	GLR(3,1,1,tt,gen) = 0d0
	GLR(3,2,2,tt,gen) = 0d0
	GLR(3,1,2,tt,gen) = cI*hf(tt,gen)/sqrt2*(Af(tt,gen)*SB + MUE*CB) 
	GLR(3,2,1,tt,gen) = -cI*hf(tt,gen)/sqrt2*(Af(tt,gen)*SB + MUE*CB) 

	GLR(4,1,1,tt,gen) = 0d0
	GLR(4,2,2,tt,gen) = 0d0
	GLR(4,1,2,tt,gen) = -cI*hf(tt,gen)/sqrt2*(Af(tt,gen)*CB - MUE*SB) 
	GLR(4,2,1,tt,gen) = cI*hf(tt,gen)/sqrt2*(Af(tt,gen)*CB - MUE*SB) 

! G_{LR}^{bt} where tt = down-type
	GLR(5,1,1,tt,gen) = hf(tt,gen)*Mf(tt,gen)*SB + hf(tt-1,gen)*Mf(tt-1,gen)*CB - EL*MW/(sqrt2*SW)*S2B
	GLR(5,2,2,tt,gen) = hf(tt-1,gen)*Mf(tt,gen)*CB + hf(tt,gen)*Mf(tt-1,gen)*SB
	GLR(5,1,2,tt,gen) = hf(tt-1,gen)*(Af(tt-1,gen)*CB + MUE*SB) 
	GLR(5,2,1,tt,gen) = hf(tt,gen)*(Af(tt,gen)*SB + MUE*CB)

	GLR(6,1,1,tt,gen) = -hf(tt,gen)*Mf(tt,gen)*CB + hf(tt-1,gen)*Mf(tt-1,gen)*SB + EL*MW/(sqrt2*SW)*S2B
	GLR(6,2,2,tt,gen) = hf(tt-1,gen)*Mf(tt,gen)*SB - hf(tt,gen)*Mf(tt-1,gen)*CB
	GLR(6,1,2,tt,gen) = hf(tt-1,gen)*(Af(tt-1,gen)*SB - MUE*CB) 
	GLR(6,2,1,tt,gen) = hf(tt,gen)*(-Af(tt,gen)*CB + MUE*SB)

	endif

! neutral Higgs
	do a=1,4
	do i=1,2
	 do j=1,2

	  HSfSf(a,i,j,tt,gen) = (0d0,0d0) 

	  do k=1,2
	   do l=1,2
	
	    HSfSf(a,i,j,tt,gen) = HSfSf(a,i,j,tt,gen) + GLR(a,k,l,tt,gen)*USf(i,k,tt,gen)*USf(j,l,tt,gen)

   	   enddo
	  enddo
 	 enddo
	enddo
	enddo

! charged Higgs
	do a=5,6
	do i=1,2
	 do j=1,2
	
	  HSfSf(a,i,j,tt,gen) = (0d0,0d0) 
	
	  do k=1,2
	   do l=1,2
	
	    HSfSf(a,i,j,tt,gen) = HSfSf(a,i,j,tt,gen) + GLR(a,k,l,tt,gen)*USf(i,k,tt,gen)*USf(j,l,tt+1,gen)
	
	   enddo
	  enddo
	 enddo
	enddo
	enddo
	
!******** HIGHER ORDER CORRECTIONS & MicrOmegas COUPLINGS *******************************

! higher order corrections													

**********************************************************************************************************************************
*** Note Patrick: Here are the modified Yukawa Couplings. This is in accordance with ArXiv ePrint: 1001.1935. The second lines ***
*** like "+ DL(1)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD + DeltaMbQCD*CA/(SA*TB))" avoid double counting of the LO and  ***
*** NLO terms. This can easily seen by expanding the first lines in DeltaMb. Note that only the QCD part has to be substracted ***
*** as we don't calculate weak corrections and DeltaMb = DeltaMbQCD + DeltaMbELW.                                              ***
**********************************************************************************************************************************

** Delta QCD for Higgs coupling - 2,3,4-loop contributions

	DeltaQCD = 1d0 + (Alfas/pi)*5.67d0 + (Alfas/pi)**2*(35.94d0 - 1.359d0*5) + (Alfas/pi)**3*(164.14d0 - 25.77d0*5 + 0.259d0*25)  
     &            + (Alfas/pi)**4*(39.34d0 - 220.9d0*5 + 9.685d0*25 - 0.0205d0*125)
	 
** Delta T for top induced corrections to Yukawa couplings

	DeltaTlighth = (Alfas/pi)**2*(1.57d0 - 2.0d0/3.0d0*dlog(QScale**2/Mf(3,3)**2) + 1.0d0/9.0d0*(dlog(Mf(4,3)**2/QScale**2)**2))*(CA/(SA*TB))
	DeltaTheavyH = (Alfas/pi)**2*(1.57d0 - 2.0d0/3.0d0*dlog(QScale**2/Mf(3,3)**2) + 1.0d0/9.0d0*(dlog(Mf(4,3)**2/QScale**2)**2))*(SA/(TB*CA))
	DeltaTA		 = (Alfas/pi)**2*(23.0d0/6.0d0 - dlog(QScale**2/Mf(3,3)**2) + 1.0d0/6.0d0*(dlog(Mf(4,3)**2/QScale**2)**2))*(1/TB**2)
	
** Delta_Mb & the resummation
	call DMNLO_CalcDeltaMb(DmbResults,3)

! Default Setting is zero
	do i=1,6
	  	DLho(i) = 0d0
	  	DRho(i) = 0d0
	enddo
	
	do i = 1, 2
		do a = 1,4
		
		ELcho(a,i,tt,gen) = 0d0
		ERho(a,i,tt,gen) =  0d0
		FLho(a,i,tt,gen) =  0d0
		FRho(a,i,tt,gen) =  0d0

		GLho(a,i,tt,gen) =  0d0
		GRho(a,i,tt,gen) =  0d0
		HLho(a,i,tt,gen) =  0d0
		HRho(a,i,tt,gen) =  0d0
	   
		enddo
	enddo
	  	  
! Resummation is only applied for generation three	  
	if (gen.eq.3)	then
	
! Down-Type Quarks
		if (tt.eq.4)	then
	 
		DLho(1) = DL(1)*dsqrt(DeltaQCD + DeltaTlighth)*(1d0/(1d0+DeltaMb)*(1d0-DeltaMb*CA/(SA*TB))) 
     &	  	    + DL(1)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD + DeltaMbQCD*CA/(SA*TB))
		DRho(1) = DR(1)*dsqrt(DeltaQCD + DeltaTlighth)*(1d0/(1d0+DeltaMb)*(1d0-DeltaMb*CA/(SA*TB))) 
     &	  	    + DR(1)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD + DeltaMbQCD*CA/(SA*TB))

		DLho(2) = DL(2)*dsqrt(DeltaQCD + DeltaTheavyH)*(1d0/(1d0+DeltaMb)*(1d0+DeltaMb*SA/(CA*TB))) 
     &	  	    + DL(2)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD - DeltaMbQCD*SA/(CA*TB))
		DRho(2) = DR(2)*dsqrt(DeltaQCD + DeltaTheavyH)*(1d0/(1d0+DeltaMb)*(1d0+DeltaMb*SA/(CA*TB))) 
     &	  	    + DR(2)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD - DeltaMbQCD*SA/(CA*TB))

		DLho(3) = DL(3)*dsqrt(DeltaQCD + DeltaTA)*(1d0/(1d0 + DeltaMb)*(1d0 - DeltaMb/(TB**2))) 
     &	  	    + DL(3)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD + DeltaMbQCD/(TB**2))
		DRho(3) = DR(3)*dsqrt(DeltaQCD + DeltaTA)*(1d0/(1d0 + DeltaMb)*(1d0 - DeltaMb/(TB**2))) 
     &	  	    + DR(3)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD + DeltaMbQCD/(TB**2))
	
		DLho(4) = DL(4)*dsqrt(DeltaQCD) + DL(4)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0)
		DRho(4) = DR(4)*dsqrt(DeltaQCD) + DR(4)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0)

		if (icha.eq.0)	then
		
		do i = 1, 2
			do a = 1, 4
		
			ELcho(a,i,tt,gen) = - EL*Mf(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,1,tt,gen)
			ERho(a,i,tt,gen) =  - EL*Mf(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,2,tt,gen) 
			FLho(a,i,tt,gen) =  - EL*Mf(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,2,tt,gen)
			FRho(a,i,tt,gen) =  - EL*Mf(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,1,tt,gen)

			GLho(a,i,tt,gen) =  - EL*Mf(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,2,tt,gen)
			GRho(a,i,tt,gen) =  - EL*Mf(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,1,tt,gen)
			HLho(a,i,tt,gen) =  - EL*Mf(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,1,tt,gen)
			HRho(a,i,tt,gen) =  - EL*Mf(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,2,tt,gen)
	 
			enddo
		enddo
		
		else if (icha.eq.2)	then
		
		do i = 1, 2
			do a = 1, 2
			
			GRho(a,i,tt,gen) =  EL*Mf(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*UCha(a,2)*USf(i,1,tt-1,gen)
			HLho(a,i,tt,gen) =  EL*Mf(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*UChaC(a,2)*USf(i,1,tt-1,gen)
	  
			enddo
		enddo
		
		endif
		
! Up-Type Quarks
		else if (tt.eq.3) then
		
		do i=1, 4
			
			DLho(i) = DL(i)*dsqrt(DeltaQCD-1d0-(Alfas/pi)*5.67d0)
			DRho(i) = DR(i)*dsqrt(DeltaQCD-1d0-(Alfas/pi)*5.67d0)
			
		enddo
			
		if (icha.eq.1)	then
		
			DLho(5) = DL(5)*dsqrt(DeltaQCD)*(1d0/(1d0+DeltaMb)) 
     &	  	    + DL(5)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD)
			DRho(5) = DR(5)*dsqrt(DeltaQCD)*(1d0/(1d0+DeltaMb)) 
     &	  	    + DR(5)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD)

			DLho(6) = DL(6)*dsqrt(DeltaQCD)*(1d0/(1d0+DeltaMb)) 
     &	  	    + DL(6)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD)
			DRho(6) = DR(6)*dsqrt(DeltaQCD)*(1d0/(1d0+DeltaMb)) 
     &	  	    + DR(6)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD)
	 
		do i = 1, 2
			do a = 1, 4
		
			FLho(a,i,tt,gen) =  -EL*Mf(tt+1,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,2,tt+1,gen)
			FRho(a,i,tt,gen) =  -EL*Mf(tt+1,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,1,tt+1,gen)
	 
			enddo
			
			do a = 1, 2
			
			GRho(a,i,tt,gen) =  EL*Mf(tt+1,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*UCha(a,2)*USf(i,1,tt,gen)
			ERho(a,i,tt,gen) =  EL*Mf(tt+1,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*UCha(a,2)*USf(i,2,tt+1,gen)
			
			enddo
		enddo
		
		else if (icha.eq.2)	then
		
		do i = 1, 2
			do a = 1, 2
			
			ERho(a,i,tt,gen) = EL*Mf(tt+1,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*UCha(a,2)*USf(i,2,tt+1,gen) 
			FLho(a,i,tt,gen) = EL*Mf(tt+1,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*UChaC(a,2)*USf(i,2,tt+1,gen)
			
			enddo
		enddo
	  
		endif
	 endif	
	endif
		
		
! Couplings for comparison with MicrOmegas

! Mfm(i) input from MO (mssm_ini.F)		
! i=1,2,3 => h0, H0, A0

c	if ((tt.eq.4).and.(gen.eq.3).and.(imicr.eq.1)) then	
c	 Mfm(1) = mbQ/(1d0+DeltaMB)*dsqrt(DeltaQCD)*(1d0 - DeltaMB*CA/(TB*SA))
c	 Mfm(2) = mbQ/(1d0+DeltaMB)*dsqrt(DeltaQCD)*(1d0 + DeltaMB*SA/(TB*CA))
c	 Mfm(3) = mbQ/(1d0+DeltaMB)*dsqrt(DeltaQCD)*(1d0 - DeltaMB/(TB**2))
c	 Mfm(4) = mbQ/(1d0+DeltaMB)*dsqrt(DeltaQCD)
c	 Mfm(5) = mbQ/(1d0+DeltaMB)*dsqrt(DeltaQCD)		
c 	endif
c
c** Micromegas effective top mass
c	if ((tt.eq.3).and.(gen.eq.3)) then
c	 Mfm(1) = mtinp
c	 Mfm(2) = mtinp
c	 Mfm(3) = mtinp
c	 Mfm(4) = mtinp
c	 Mfm(5) = mtinp
c	endif

	 Mfm(1) = Mf(tt,gen)
	 Mfm(2) = Mf(tt,gen)
	 Mfm(3) = Mf(tt,gen)
	 Mfm(4) = Mf(tt,gen)
	 Mfm(5) = Mf(tt,gen)

! FerFerHiggs & NeuFerSfer for MicrOmegas including the resummation

	if (((tt.eq.3).or.(tt.eq.4)).and.(gen.eq.3)) then

! bottom quark resummation
  	 if(tt.eq.4) then

	 	DLM(1) = EL*Mfm(1)*SA/(SW*2d0*MW*CB) 
	 	DRM(1) = EL*Mfm(1)*SA/(SW*2d0*MW*CB) 
	 	DLM(2) = -EL*Mfm(2)*CA/(SW*2d0*MW*CB) 
	 	DRM(2) = -EL*Mfm(2)*CA/(SW*2d0*MW*CB) 
	 	DLM(3) = -cI*EL*Mfm(3)*SB/(SW*2d0*MW*CB)
 	 	DRM(3) = cI*EL*Mfm(3)*SB/(SW*2d0*MW*CB) 
	 	DLM(4) = cI*EL*Mfm(4)/(SW*2d0*MW)
 	 	DRM(4) = -cI*EL*Mfm(4)/(SW*2d0*MW) 
! top quark resummation
	 else	 
		DLM(1) = -EL*Mfm(1)*CA/(SW*2d0*MW*SB) 
		DRM(1) = -EL*Mfm(1)*CA/(SW*2d0*MW*SB) 
		DLM(2) = -EL*Mfm(2)*SA/(SW*2d0*MW*SB) 
		DRM(2) = -EL*Mfm(2)*SA/(SW*2d0*MW*SB) 
		DLM(3) = -cI*EL*Mfm(3)*CB/(SW*2d0*MW*SB) 
		DRM(3) = cI*EL*Mfm(3)*CB/(SW*2d0*MW*SB) 
		DLM(4) = -cI*EL*Mfm(4)/(SW*2d0*MW) 
		DRM(4) = cI*EL*Mfm(4)/(SW*2d0*MW) 
	 endif

	 do i = 1, 2

	  if ((tt.eq.1).or.(tt.eq.3)) then

c =================================================
c NO SWITCH E<->H, F<->G BELOW THIS LINE YET !!!!!!
c =================================================
	ELcM(i,tt,gen) = EL*Mfm(5)/(SW*sqrt2*MW*SB)*ZNeuC(1,4)*USf(i,1,tt,gen) -
     - 			sqrt2*EL*ef(tt)/CW*ZNeuC(1,1)*USf(i,2,tt,gen)
	ERM(i,tt,gen) =  EL*Mfm(5)/(SW*sqrt2*MW*SB)*ZNeu(1,4)*USf(i,2,tt,gen) +
     - 			sqrt2*EL/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(1,1) + I3f(tt)*ZNeu(1,2))*USf(i,1,tt,gen)
	FLM(i,tt,gen) =  EL*Mfm(5)/(SW*sqrt2*MW*SB)*ZNeuC(1,4)*USf(i,2,tt,gen) +
     - 			sqrt2*EL/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeuC(1,1) + I3f(tt)*ZNeuC(1,2))*USf(i,1,tt,gen)
	FRM(i,tt,gen) =  EL*Mfm(5)/(SW*sqrt2*MW*SB)*ZNeu(1,4)*USf(i,1,tt,gen) -
     - 			sqrt2*EL*ef(tt)/CW*ZNeu(1,1)*USf(i,2,tt,gen)

	GLM(i,tt,gen) =  EL*Mfm(5)/(SW*sqrt2*MW*SB)*ZNeuC(1,4)*USf(i,2,tt,gen) +
     - 			sqrt2*EL/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeuC(1,1) + I3f(tt)*ZNeuC(1,2))*USf(i,1,tt,gen)
	GRM(i,tt,gen) =  EL*Mfm(5)/(SW*sqrt2*MW*SB)*ZNeu(1,4)*USf(i,1,tt,gen) -
     - 			sqrt2*EL*ef(tt)/CW*ZNeu(1,1)*USf(i,2,tt,gen)
	HLM(i,tt,gen) =  EL*Mfm(5)/(SW*sqrt2*MW*SB)*ZNeuC(1,4)*USf(i,1,tt,gen) -
     - 			sqrt2*EL*ef(tt)/CW*ZNeuC(1,1)*USf(i,2,tt,gen)
	HRM(i,tt,gen) =  EL*Mfm(5)/(SW*sqrt2*MW*SB)*ZNeu(1,4)*USf(i,2,tt,gen) +
     - 			sqrt2*EL/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(1,1) + I3f(tt)*ZNeu(1,2))*USf(i,1,tt,gen)
	
	else
	
	  ELcM(i,tt,gen) = EL*Mfm(5)/(SW*sqrt2*MW*CB)*ZNeuC(1,3)*USf(i,1,tt,gen) -
     - 			sqrt2*EL*ef(tt)/CW*ZNeuC(1,1)*USf(i,2,tt,gen)
	  ERM(i,tt,gen) =  EL*Mfm(5)/(SW*sqrt2*MW*CB)*ZNeu(1,3)*USf(i,2,tt,gen) +
     - 			sqrt2*EL/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(1,1) + I3f(tt)*ZNeu(1,2))*USf(i,1,tt,gen)
	  FLM(i,tt,gen) =  EL*Mfm(5)/(SW*sqrt2*MW*CB)*ZNeuC(1,3)*USf(i,2,tt,gen) +
     - 			sqrt2*EL/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeuC(1,1) + I3f(tt)*ZNeuC(1,2))*USf(i,1,tt,gen)
	  FRM(i,tt,gen) =  EL*Mfm(5)/(SW*sqrt2*MW*CB)*ZNeu(1,3)*USf(i,1,tt,gen) -
     - 			sqrt2*EL*ef(tt)/CW*ZNeu(1,1)*USf(i,2,tt,gen)

	  GLM(i,tt,gen) =  EL*Mfm(5)/(SW*sqrt2*MW*CB)*ZNeuC(1,3)*USf(i,2,tt,gen) +
     - 			sqrt2*EL/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeuC(1,1) + I3f(tt)*ZNeuC(1,2))*USf(i,1,tt,gen)
	  GRM(i,tt,gen) =  EL*Mfm(5)/(SW*sqrt2*MW*CB)*ZNeu(1,3)*USf(i,1,tt,gen) -
     - 			sqrt2*EL*ef(tt)/CW*ZNeu(1,1)*USf(i,2,tt,gen)
	  HLM(i,tt,gen) =  EL*Mfm(5)/(SW*sqrt2*MW*CB)*ZNeuC(1,3)*USf(i,1,tt,gen) -
     - 			sqrt2*EL*ef(tt)/CW*ZNeuC(1,1)*USf(i,2,tt,gen)
	  HRM(i,tt,gen) =  EL*Mfm(5)/(SW*sqrt2*MW*CB)*ZNeu(1,3)*USf(i,2,tt,gen) +
     - 			sqrt2*EL/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(1,1) + I3f(tt)*ZNeu(1,2))*USf(i,1,tt,gen)
	
	 endif
	
	 enddo
	 
	 endif

	
	end
