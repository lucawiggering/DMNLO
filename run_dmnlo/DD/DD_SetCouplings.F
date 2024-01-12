!==================================================================================================
!======================= Calculates all needed elementary couplings	 ==============================
!======================= Copy and expansion of ChiChi2QQ_Couplings.F ==============================
!==================================================================================================

	subroutine DD_SetCouplings
	implicit none

#include "DD_Model.h"
#include "DD_GenCouplings.h"

	double complex GLR(6,2,2,4,3),ZLR(2,2,4,3)
	integer i,j,k,l,a
	integer icha, indi, indj, tt, gen	
	
! functions from DMNLO_Kinematics.F
	double precision Krondelta, ef, I3f, hf
	double precision hfh
	
c	open (unit=20,file="couplings.txt",action="write",status="old",position="append")
	
! default settings for direct detection
	icha = 0
	indi = 1
	indj = 1

! Loop over tt and gen
	do tt = 3,4
	do gen = 1,3
	
! ******************************************************
! Gaugino-Vector boson & Fermion-Vector boson couplings
! ******************************************************
	if (icha.eq.2) then
	
	ApL(tt,gen) = -EL*Krondelta(indi,indj)	
	BpL(tt,gen) = -EL*ef(tt)

	AL(tt,gen) = EL/(SW*CW)*(-VCha(indj,1)*VChaC(indi,1) - 1/2d0*VCha(indj,2)*VChaC(indi,2) + Krondelta(indi,indj)*SW**2)
	AR(tt,gen) = EL/(SW*CW)*(-UChaC(indj,1)*UCha(indi,1) - 1/2d0*UChaC(indj,2)*UCha(indi,2) + Krondelta(indi,indj)*SW**2)
	BL(tt,gen) = EL/(SW*CW)*(ef(tt)*SW**2-I3f(tt))
	BR(tt,gen) = EL/(SW*CW)*(ef(tt)*SW**2)
	
	else if (icha.eq.1) then

	AL(tt,gen) = EL/SW*( ZNeu(indj,2)*VChaC(indi,1) - 1/sqrt2*ZNeu(indj,4)*VChaC(indi,2))
	AR(tt,gen) = EL/SW*( ZNeuC(indj,2)*UCha(indi,1) + 1/sqrt2*ZNeuC(indj,3)*UCha(indi,2))
	BL(tt,gen) = -EL/(SW*sqrt2)
	BR(tt,gen) = 0d0

	else if (icha.eq.0) then

	AL(tt,gen) = EL/(2*SW*CW)*(- ZNeu(indj,3)*ZNeuC(indi,3) + ZNeu(indj,4)*ZNeuC(indi,4))
	AR(tt,gen) = EL/(2*SW*CW)*(  ZNeuC(indj,3)*ZNeu(indi,3) - ZNeuC(indj,4)*ZNeu(indi,4))
	BL(tt,gen) = EL/(SW*CW)*(ef(tt)*SW**2-I3f(tt))
	BR(tt,gen) = EL/(SW*CW)*(ef(tt)*SW**2)

	endif
	
! ************************	
! Gaugino-Higgs couplings
! ************************
! Higgs index Hk0 = (h0,H0,A0,G0) = (1,2,3,4) or Hk = (h0,H0,A0,G0,Hp,Gp) = (1,2,3,4,5,6)

	if (icha.eq.0) then
! NeuNeuHiggs
	CL(1,tt,gen) = -EL/(2d0*SW)*((-SA)*(ZNeuC(indi,3)*ZNeuC(indj,2) + ZNeuC(indj,3)*ZNeuC(indi,2) - SW/CW*(ZNeuC(indi,3)*ZNeuC(indj,1) + ZNeuC(indj,3)*ZNeuC(indi,1)))
     &	        + (-CA)*(ZNeuC(indi,4)*ZNeuC(indj,2) + ZNeuC(indj,4)*ZNeuC(indi,2) - SW/CW*(ZNeuC(indi,4)*ZNeuC(indj,1) + ZNeuC(indj,4)*ZNeuC(indi,1))) )
	CR(1,tt,gen) = -EL/(2d0*SW)*((-SA)*( ZNeu(indi,3)*ZNeu(indj,2)  +  ZNeu(indj,3)*ZNeu(indi,2)  - SW/CW*( ZNeu(indi,3)*ZNeu(indj,1)  +  ZNeu(indj,3)*ZNeu(indi,1)))
     &	        + (-CA)*( ZNeu(indi,4)*ZNeu(indj,2)  +  ZNeu(indj,4)*ZNeu(indi,2)  - SW/CW*( ZNeu(indi,4)*ZNeu(indj,1)  +  ZNeu(indj,4)*ZNeu(indi,1))) )
	CL(2,tt,gen) = -EL/(2d0*SW)*((+CA)*(ZNeuC(indi,3)*ZNeuC(indj,2) + ZNeuC(indj,3)*ZNeuC(indi,2) - SW/CW*(ZNeuC(indi,3)*ZNeuC(indj,1) + ZNeuC(indj,3)*ZNeuC(indi,1)))
     &	        + (-SA)*(ZNeuC(indi,4)*ZNeuC(indj,2) + ZNeuC(indj,4)*ZNeuC(indi,2) - SW/CW*(ZNeuC(indi,4)*ZNeuC(indj,1) + ZNeuC(indj,4)*ZNeuC(indi,1))) )
	CR(2,tt,gen) = -EL/(2d0*SW)*((+CA)*( ZNeu(indi,3)*ZNeu(indj,2)  +  ZNeu(indj,3)*ZNeu(indi,2)  - SW/CW*( ZNeu(indi,3)*ZNeu(indj,1)  +  ZNeu(indj,3)*ZNeu(indi,1)))
     &	        + (-SA)*( ZNeu(indi,4)*ZNeu(indj,2)  +  ZNeu(indj,4)*ZNeu(indi,2)  - SW/CW*( ZNeu(indi,4)*ZNeu(indj,1)  +  ZNeu(indj,4)*ZNeu(indi,1))) )
	CL(3,tt,gen) = -cI*EL/(2d0*SW)*((-SB)*(ZNeuC(indi,3)*ZNeuC(indj,2) + ZNeuC(indj,3)*ZNeuC(indi,2) - SW/CW*(ZNeuC(indi,3)*ZNeuC(indj,1) + ZNeuC(indj,3)*ZNeuC(indi,1)))
     &	           + (+CB)*(ZNeuC(indi,4)*ZNeuC(indj,2) + ZNeuC(indj,4)*ZNeuC(indi,2) - SW/CW*(ZNeuC(indi,4)*ZNeuC(indj,1) + ZNeuC(indj,4)*ZNeuC(indi,1))) )
	CR(3,tt,gen) = +cI*EL/(2d0*SW)*((-SB)*( ZNeu(indi,3)*ZNeu(indj,2)  +  ZNeu(indj,3)*ZNeu(indi,2)  - SW/CW*( ZNeu(indi,3)*ZNeu(indj,1)  +  ZNeu(indj,3)*ZNeu(indi,1)))
     &	           + (+CB)*( ZNeu(indi,4)*ZNeu(indj,2)  +  ZNeu(indj,4)*ZNeu(indi,2)  - SW/CW*( ZNeu(indi,4)*ZNeu(indj,1)  +  ZNeu(indj,4)*ZNeu(indi,1))) )
	CL(4,tt,gen) = -cI*EL/(2d0*SW)*((+CB)*(ZNeuC(indi,3)*ZNeuC(indj,2) + ZNeuC(indj,3)*ZNeuC(indi,2) - SW/CW*(ZNeuC(indi,3)*ZNeuC(indj,1) + ZNeuC(indj,3)*ZNeuC(indi,1)))
     &	           + (+SB)*(ZNeuC(indi,4)*ZNeuC(indj,2) + ZNeuC(indj,4)*ZNeuC(indi,2) - SW/CW*(ZNeuC(indi,4)*ZNeuC(indj,1) + ZNeuC(indj,4)*ZNeuC(indi,1))) )
	CR(4,tt,gen) = +cI*EL/(2d0*SW)*((+CB)*( ZNeu(indi,3)*ZNeu(indj,2)  +  ZNeu(indj,3)*ZNeu(indi,2)  - SW/CW*( ZNeu(indi,3)*ZNeu(indj,1)  +  ZNeu(indj,3)*ZNeu(indi,1)))
     &	           + (+SB)*( ZNeu(indi,4)*ZNeu(indj,2)  +  ZNeu(indj,4)*ZNeu(indi,2)  - SW/CW*( ZNeu(indi,4)*ZNeu(indj,1)  +  ZNeu(indj,4)*ZNeu(indi,1))) )

! NeuChaHiggs
	CL(5,tt,gen) = 0d0
	CR(5,tt,gen) = 0d0
	CL(6,tt,gen) = 0d0
	CR(6,tt,gen) = 0d0

	else if (icha.eq.1) then
! NeuNeuHiggs
	CL(1,tt,gen) = 0d0
	CR(1,tt,gen) = 0d0
	CL(2,tt,gen) = 0d0	
	CR(2,tt,gen) = 0d0
	CL(3,tt,gen) = 0d0
	CR(3,tt,gen) = 0d0  
	CL(4,tt,gen) = 0d0
	CR(4,tt,gen) = 0d0

! NeuChaHiggs
	CL(5,tt,gen) = -EL/SW*CB*(VChaC(indi,1)*ZNeuC(indj,4) + 1/sqrt2*(ZNeuC(indj,2) + SW/CW*ZNeuC(indj,1))*VChaC(indi,2)) 
	CR(5,tt,gen) = -EL/SW*SB*(UCha(indi,1)*ZNeu(indj,3) - 1/sqrt2*(ZNeu(indj,2) + SW/CW*ZNeu(indj,1))*UCha(indi,2)) 
	CL(6,tt,gen) = -EL/SW*SB*(VChaC(indi,1)*ZNeuC(indj,4) + 1/sqrt2*(ZNeuC(indj,2) + SW/CW*ZNeuC(indj,1))*VChaC(indi,2)) 
	CR(6,tt,gen) =  EL/SW*CB*(UCha(indi,1)*ZNeu(indj,3) - 1/sqrt2*(ZNeu(indj,2) + SW/CW*ZNeu(indj,1))*UCha(indi,2))  
	
	else if (icha.eq.2) then	
! ChaChaHiggs
	CL(1,tt,gen) = -EL/(SW*sqrt2)*(-SA*VChaC(indi,1)*UChaC(indj,2) + CA*VChaC(indi,2)*UChaC(indj,1))
	CR(1,tt,gen) = -EL/(SW*sqrt2)*(-SA*VCha(indj,1)*UCha(indi,2)   + CA*VCha(indj,2)*UCha(indi,1))
	CL(2,tt,gen) = -EL/(SW*sqrt2)*( CA*VChaC(indi,1)*UChaC(indj,2) + SA*VChaC(indi,2)*UChaC(indj,1))
	CR(2,tt,gen) = -EL/(SW*sqrt2)*( CA*VCha(indj,1)*UCha(indi,2) + SA*VCha(indj,2)*UCha(indi,1))
	CL(3,tt,gen) = -cI*EL/(SW*sqrt2)*(-SB*VChaC(indi,1)*UChaC(indj,2) - CB*VChaC(indi,2)*UChaC(indj,1))
	CR(3,tt,gen) = cI*EL/(SW*sqrt2)*(-SB*VCha(indj,1)*UCha(indi,2) - CB*VCha(indj,2)*UCha(indi,1))
	CL(4,tt,gen) = -cI*EL/(SW*sqrt2)*( CB*VChaC(indi,1)*UChaC(indj,2) - SB*VChaC(indi,2)*UChaC(indj,1))
	CR(4,tt,gen) = cI*EL/(SW*sqrt2)*( CB*VCha(indj,1)*UCha(indi,2) - SB*VCha(indj,2)*UCha(indi,1))

! NeuChaHiggs
	CL(5,tt,gen) = 0d0
	CR(5,tt,gen) = 0d0
	CL(6,tt,gen) = 0d0
	CR(6,tt,gen) = 0d0

	endif
	
!**************************	
! Fermion-Higgs couplings *
!**************************

	if (icha.eq.1) then
		
! FerFerHiggs
	 if ((tt.eq.1).or.(tt.eq.3)) then

  	  DL(1,tt,gen) = -EL*Mf_high(tt,gen)*CA/(SW*2d0*MW*SB) 
  	  DR(1,tt,gen) = -EL*Mf_high(tt,gen)*CA/(SW*2d0*MW*SB) 
  	  DL(2,tt,gen) = -EL*Mf_high(tt,gen)*SA/(SW*2d0*MW*SB) 
  	  DR(2,tt,gen) = -EL*Mf_high(tt,gen)*SA/(SW*2d0*MW*SB) 
  	  DL(3,tt,gen) = -cI*EL*Mf_high(tt,gen)*CB/(SW*2d0*MW*SB) 
  	  DR(3,tt,gen) = cI*EL*Mf_high(tt,gen)*CB/(SW*2d0*MW*SB) 
  	  DL(4,tt,gen) = -cI*EL*Mf_high(tt,gen)/(SW*2d0*MW) 
  	  DR(4,tt,gen) = cI*EL*Mf_high(tt,gen)/(SW*2d0*MW) 

	 else

	  DL(1,tt,gen) = EL*Mf_high(tt,gen)*SA/(SW*2d0*MW*CB) 
	  DR(1,tt,gen) = EL*Mf_high(tt,gen)*SA/(SW*2d0*MW*CB) 
	  DL(2,tt,gen) = -EL*Mf_high(tt,gen)*CA/(SW*2d0*MW*CB) 
	  DR(2,tt,gen) = -EL*Mf_high(tt,gen)*CA/(SW*2d0*MW*CB) 
	  DL(3,tt,gen) = -cI*EL*Mf_high(tt,gen)*SB/(SW*2d0*MW*CB)
 	  DR(3,tt,gen) = cI*EL*Mf_high(tt,gen)*SB/(SW*2d0*MW*CB) 
	  DL(4,tt,gen) = cI*EL*Mf_high(tt,gen)/(SW*2d0*MW)
 	  DR(4,tt,gen) = -cI*EL*Mf_high(tt,gen)/(SW*2d0*MW) 

	 endif

! Fer1Fer2Higgs (tt = type of the up-type fermion)
	DL(5,tt,gen) = EL*Mf_high(tt,gen)*CB/(SW*sqrt2*MW*SB)
	DR(5,tt,gen) = EL*Mf_high(tt+1,gen)*SB/(SW*sqrt2*MW*CB)
	DL(6,tt,gen) = EL*Mf_high(tt,gen)/(SW*sqrt2*MW)
	DR(6,tt,gen) = -EL*Mf_high(tt+1,gen)/(SW*sqrt2*MW)

	else
! Yukawa coupling of the neutral Higgses to fermions

! FerFerHiggs
	if ((tt.eq.1).or.(tt.eq.3)) then

	DL(1,tt,gen) = -EL*Mf_high(tt,gen)*CA/(SW*2d0*MW*SB) 
	DR(1,tt,gen) = -EL*Mf_high(tt,gen)*CA/(SW*2d0*MW*SB) 
	DL(2,tt,gen) = -EL*Mf_high(tt,gen)*SA/(SW*2d0*MW*SB) 
	DR(2,tt,gen) = -EL*Mf_high(tt,gen)*SA/(SW*2d0*MW*SB) 
	DL(3,tt,gen) = -cI*EL*Mf_high(tt,gen)*CB/(SW*2d0*MW*SB) 
	DR(3,tt,gen) = cI*EL*Mf_high(tt,gen)*CB/(SW*2d0*MW*SB) 
	DL(4,tt,gen) = -cI*EL*Mf_high(tt,gen)/(SW*2d0*MW) 
	DR(4,tt,gen) = cI*EL*Mf_high(tt,gen)/(SW*2d0*MW) 	
	
	else

	DL(1,tt,gen) = EL*Mf_high(tt,gen)*SA/(SW*2d0*MW*CB) 
	DR(1,tt,gen) = EL*Mf_high(tt,gen)*SA/(SW*2d0*MW*CB) 
	DL(2,tt,gen) = -EL*Mf_high(tt,gen)*CA/(SW*2d0*MW*CB) 
	DR(2,tt,gen) = -EL*Mf_high(tt,gen)*CA/(SW*2d0*MW*CB) 
	DL(3,tt,gen) = -cI*EL*Mf_high(tt,gen)*SB/(SW*2d0*MW*CB)
 	DR(3,tt,gen) = cI*EL*Mf_high(tt,gen)*SB/(SW*2d0*MW*CB) 
	DL(4,tt,gen) = cI*EL*Mf_high(tt,gen)/(SW*2d0*MW)												
 	DR(4,tt,gen) = -cI*EL*Mf_high(tt,gen)/(SW*2d0*MW) 

	endif

! Fer1Fer2Higgs
	DL(5,tt,gen) = 0d0
	DR(5,tt,gen) = 0d0
	DL(6,tt,gen) = 0d0
	DR(6,tt,gen) = 0d0
	
	endif
	
!*************************************
! Gaugino-Fermion-Sfermion couplings *
!*************************************

! NeuFerSfer
	if (icha.eq.0) then

	 do a= 1, 4
	  do i = 1, 2
        
	  if ((tt.eq.1).or.(tt.eq.3)) then
        
	  HL(a,i,tt,gen) = -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*USf(i,2,tt,gen))
	  HR(a,i,tt,gen) =  -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*SB)*ZNeu(a,4)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*USf(i,1,tt,gen))
	  GL(a,i,tt,gen) =  -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeuC(a,1) + I3f(tt)*ZNeuC(a,2))*USf(i,1,tt,gen))
	  GR(a,i,tt,gen) =  -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*SB)*ZNeu(a,4)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeu(a,1)*USf(i,2,tt,gen))
        
	  FL(a,i,tt,gen) =  -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeuC(a,1) + I3f(tt)*ZNeuC(a,2))*USf(i,1,tt,gen))
	  FR(a,i,tt,gen) =  -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*SB)*ZNeu(a,4)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeu(a,1)*USf(i,2,tt,gen))
	  ELc(a,i,tt,gen) =  -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*USf(i,2,tt,gen))
	  ER(a,i,tt,gen) =  -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*SB)*ZNeu(a,4)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*USf(i,1,tt,gen))
	  
	  else

	  HL(a,i,tt,gen) = -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*USf(i,2,tt,gen))
	  HR(a,i,tt,gen) =  -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*USf(i,1,tt,gen))
	  GL(a,i,tt,gen) =  -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeuC(a,1) + I3f(tt)*ZNeuC(a,2))*USf(i,1,tt,gen))
	  GR(a,i,tt,gen) =  -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeu(a,1)*USf(i,2,tt,gen))
        
	  FL(a,i,tt,gen) =  -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeuC(a,1) + I3f(tt)*ZNeuC(a,2))*USf(i,1,tt,gen))
	  FR(a,i,tt,gen) =  -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeu(a,1)*USf(i,2,tt,gen))
	  ELc(a,i,tt,gen) =  -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*USf(i,2,tt,gen))
	  ER(a,i,tt,gen) =  -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,2,tt,gen) +
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
   
	  HL(a,i,tt,gen) = -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*USf(i,2,tt,gen))
	  HR(a,i,tt,gen) =  -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*SB)*ZNeu(a,4)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*USf(i,1,tt,gen))
   
	  FL(a,i,tt,gen) =  -EL*(Mf_high(tt+1,gen)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,2,tt+1,gen) +
     &   			sqrt2/SW*((ef(tt+1) - I3f(tt+1))*SW/CW*ZNeuC(a,1) + I3f(tt+1)*ZNeuC(a,2))*USf(i,1,tt+1,gen))
	  FR(a,i,tt,gen) =  -EL*(Mf_high(tt+1,gen)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,1,tt+1,gen) -
     &   			sqrt2*ef(tt+1)/CW*ZNeu(a,1)*USf(i,2,tt+1,gen))
	  
	  else
   
	  HL(a,i,tt,gen) = -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,1,tt,gen) -
     &   			sqrt2*ef(tt)/CW*ZNeuC(a,1)*USf(i,2,tt,gen))
	  HR(a,i,tt,gen) =  -EL*(Mf_high(tt,gen)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,2,tt,gen) +
     &   			sqrt2/SW*((ef(tt) - I3f(tt))*SW/CW*ZNeu(a,1) + I3f(tt)*ZNeu(a,2))*USf(i,1,tt,gen))
   
	  FL(a,i,tt,gen) =  -EL*(Mf_high(tt-1,gen)/(SW*sqrt2*MW*SB)*ZNeuC(a,4)*USf(i,2,tt-1,gen) +
     &   			sqrt2/SW*((ef(tt-1) - I3f(tt-1))*SW/CW*ZNeuC(a,1) + I3f(tt-1)*ZNeuC(a,2))*USf(i,1,tt-1,gen))
	  FR(a,i,tt,gen) =  -EL*(Mf_high(tt-1,gen)/(SW*sqrt2*MW*SB)*ZNeu(a,4)*USf(i,1,tt-1,gen) -
     &   			sqrt2*ef(tt-1)/CW*ZNeu(a,1)*USf(i,2,tt-1,gen))
   
	  endif
   	  enddo

! chargino couplings
	  do a= 1, 2

	    if ((tt.eq.1).or.(tt.eq.3)) then

	    GL(a,i,tt,gen) =  -EL/SW*VChaC(a,1)*USf(i,1,tt,gen) + EL*Mf_high(tt,gen)/(SW*sqrt2*MW*SB)*VChaC(a,2)*USf(i,2,tt,gen)
	    GR(a,i,tt,gen) =  EL*Mf_high(tt+1,gen)/(SW*sqrt2*MW*CB)*UCha(a,2)*USf(i,1,tt,gen)

	    ELc(a,i,tt,gen) = EL*Mf_high(tt,gen)/(SW*sqrt2*MW*SB)*VChaC(a,2)*USf(i,1,tt+1,gen) 
	    ER(a,i,tt,gen) = -EL/SW*UCha(a,1)*USf(i,1,tt+1,gen) + EL*Mf_high(tt+1,gen)/(SW*sqrt2*MW*CB)*UCha(a,2)*USf(i,2,tt+1,gen) 
	    
	    else

	    GL(a,i,tt,gen) =  -EL/SW*UChaC(a,1)*USf(i,1,tt,gen) + EL*Mf_high(tt,gen)/(SW*sqrt2*MW*CB)*UChaC(a,2)*USf(i,2,tt,gen)
	    GR(a,i,tt,gen) =  EL*Mf_high(tt-1,gen)/(SW*sqrt2*MW*SB)*VCha(a,2)*USf(i,1,tt,gen)
	
	    ELc(a,i,tt,gen) =  EL*Mf_high(tt,gen)/(SW*sqrt2*MW*CB)*UChaC(a,2)*USf(i,1,tt-1,gen)  
	    ER(a,i,tt,gen) = -EL/SW*VCha(a,1)*USf(i,1,tt-1,gen) + EL*Mf_high(tt-1,gen)/(SW*sqrt2*MW*SB)*VCha(a,2)*USf(i,2,tt-1,gen) 
          
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
	  ELc(a,i,tt,gen) =  EL*Mf_high(tt,gen)/(SW*sqrt2*MW*SB)*VChaC(a,2)*USf(i,1,tt+1,gen) 
	  ER(a,i,tt,gen) =  -EL/SW*UCha(a,1)*USf(i,1,tt+1,gen) + EL*Mf_high(tt+1,gen)/(SW*sqrt2*MW*CB)*UCha(a,2)*USf(i,2,tt+1,gen) 

	  FL(a,i,tt,gen) =  -EL/SW*UChaC(a,1)*USf(i,1,tt+1,gen) + EL*Mf_high(tt+1,gen)/(SW*sqrt2*MW*CB)*UChaC(a,2)*USf(i,2,tt+1,gen)
	  FR(a,i,tt,gen)  =  EL*Mf_high(tt,gen)/(SW*sqrt2*MW*SB)*VCha(a,2)*USf(i,1,tt+1,gen)
	  
	  else

c HL=conj(GR), HR=conj(GL)
	  GL(a,i,tt,gen)  =  -EL/SW*VChaC(a,1)*USf(i,1,tt-1,gen) + EL*Mf_high(tt-1,gen)/(SW*sqrt2*MW*SB)*VChaC(a,2)*USf(i,2,tt-1,gen)
	  GR(a,i,tt,gen) =  EL*Mf_high(tt,gen)/(SW*sqrt2*MW*CB)*UCha(a,2)*USf(i,1,tt-1,gen)

	  HL(a,i,tt,gen) =  EL*Mf_high(tt,gen)/(SW*sqrt2*MW*CB)*UChaC(a,2)*USf(i,1,tt-1,gen)
	  HR(a,i,tt,gen)  =  -EL/SW*VCha(a,1)*USf(i,1,tt-1,gen) + EL*Mf_high(tt-1,gen)/(SW*sqrt2*MW*SB)*VCha(a,2)*USf(i,2,tt-1,gen)

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

	GLR(1,1,1,tt,gen) = -sqrt2*hfh(tt,gen)*Mf_high(tt,gen)*CA + EL*MZ*(I3f(tt)-ef(tt)*SW2)*SAB/(SW*CW) 
	GLR(1,2,2,tt,gen) = -sqrt2*hfh(tt,gen)*Mf_high(tt,gen)*CA + EL*MZ*ef(tt)*SW2*SAB/(SW*CW) 
	GLR(1,1,2,tt,gen) = -hfh(tt,gen)/sqrt2*(Af(tt,gen)*CA + MUE*SA) 
	GLR(1,2,1,tt,gen) = -hfh(tt,gen)/sqrt2*(Af(tt,gen)*CA + MUE*SA) 

	GLR(2,1,1,tt,gen) = -sqrt2*hfh(tt,gen)*Mf_high(tt,gen)*SA - EL*MZ*(I3f(tt)-ef(tt)*SW2)*CAB/(SW*CW) 
	GLR(2,2,2,tt,gen) = -sqrt2*hfh(tt,gen)*Mf_high(tt,gen)*SA - EL*MZ*ef(tt)*SW2*CAB/(SW*CW) 
	GLR(2,1,2,tt,gen) = -hfh(tt,gen)/sqrt2*(Af(tt,gen)*SA - MUE*CA) 
	GLR(2,2,1,tt,gen) = -hfh(tt,gen)/sqrt2*(Af(tt,gen)*SA - MUE*CA) 

	GLR(3,1,1,tt,gen) = 0d0
	GLR(3,2,2,tt,gen) = 0d0
	GLR(3,1,2,tt,gen) = cI*hfh(tt,gen)/sqrt2*(Af(tt,gen)*CB + MUE*SB) 
	GLR(3,2,1,tt,gen) = -cI*hfh(tt,gen)/sqrt2*(Af(tt,gen)*CB + MUE*SB) 

	GLR(4,1,1,tt,gen) = 0d0
	GLR(4,2,2,tt,gen) = 0d0
	GLR(4,1,2,tt,gen) = cI*hfh(tt,gen)/sqrt2*(Af(tt,gen)*SB - MUE*CB) 
	GLR(4,2,1,tt,gen) = -cI*hfh(tt,gen)/sqrt2*(Af(tt,gen)*SB - MUE*CB) 

! G_{LR}^{tb} where tt = up-type
	GLR(5,1,1,tt,gen) = hfh(tt+1,gen)*Mf_high(tt+1,gen)*SB + hfh(tt,gen)*Mf_high(tt,gen)*CB - EL*MW/(sqrt2*SW)*S2B
	GLR(5,2,2,tt,gen) = hfh(tt,gen)*Mf_high(tt+1,gen)*CB + hfh(tt+1,gen)*Mf_high(tt,gen)*SB
	GLR(5,1,2,tt,gen) = hfh(tt+1,gen)*(Af(tt+1,gen)*SB + MUE*CB) 
	GLR(5,2,1,tt,gen) = hfh(tt,gen)*(Af(tt,gen)*CB + MUE*SB)

	GLR(6,1,1,tt,gen) = -hfh(tt+1,gen)*Mf_high(tt+1,gen)*CB + hfh(tt,gen)*Mf_high(tt,gen)*SB + EL*MW/(sqrt2*SW)*S2B
	GLR(6,2,2,tt,gen) = hfh(tt,gen)*Mf_high(tt+1,gen)*SB - hfh(tt+1,gen)*Mf_high(tt,gen)*CB
	GLR(6,1,2,tt,gen) = hfh(tt+1,gen)*(-Af(tt+1,gen)*CB + MUE*SB) 
	GLR(6,2,1,tt,gen) = hfh(tt,gen)*(Af(tt,gen)*SB - MUE*CB)

	else

	GLR(1,1,1,tt,gen) = sqrt2*hfh(tt,gen)*Mf_high(tt,gen)*SA + EL*MZ*(I3f(tt)-ef(tt)*SW2)*SAB/(SW*CW) 
	GLR(1,2,2,tt,gen) = sqrt2*hfh(tt,gen)*Mf_high(tt,gen)*SA + EL*MZ*ef(tt)*SW2*SAB/(SW*CW) 
	GLR(1,1,2,tt,gen) = hfh(tt,gen)/sqrt2*(Af(tt,gen)*SA + MUE*CA) 
	GLR(1,2,1,tt,gen) = hfh(tt,gen)/sqrt2*(Af(tt,gen)*SA + MUE*CA) 

	GLR(2,1,1,tt,gen) = -sqrt2*hfh(tt,gen)*Mf_high(tt,gen)*CA - EL*MZ*(I3f(tt)-ef(tt)*SW2)*CAB/(SW*CW) 
	GLR(2,2,2,tt,gen) = -sqrt2*hfh(tt,gen)*Mf_high(tt,gen)*CA - EL*MZ*ef(tt)*SW2*CAB/(SW*CW) 
	GLR(2,1,2,tt,gen) = -hfh(tt,gen)/sqrt2*(Af(tt,gen)*CA - MUE*SA) 
	GLR(2,2,1,tt,gen) = -hfh(tt,gen)/sqrt2*(Af(tt,gen)*CA - MUE*SA) 

	GLR(3,1,1,tt,gen) = 0d0
	GLR(3,2,2,tt,gen) = 0d0
	GLR(3,1,2,tt,gen) = cI*hfh(tt,gen)/sqrt2*(Af(tt,gen)*SB + MUE*CB) 
	GLR(3,2,1,tt,gen) = -cI*hfh(tt,gen)/sqrt2*(Af(tt,gen)*SB + MUE*CB) 

	GLR(4,1,1,tt,gen) = 0d0
	GLR(4,2,2,tt,gen) = 0d0
	GLR(4,1,2,tt,gen) = -cI*hfh(tt,gen)/sqrt2*(Af(tt,gen)*CB - MUE*SB) 
	GLR(4,2,1,tt,gen) = cI*hfh(tt,gen)/sqrt2*(Af(tt,gen)*CB - MUE*SB) 

! G_{LR}^{bt} where tt = down-type
	GLR(5,1,1,tt,gen) = hfh(tt,gen)*Mf_high(tt,gen)*SB + hfh(tt-1,gen)*Mf_high(tt-1,gen)*CB - EL*MW/(sqrt2*SW)*S2B
	GLR(5,2,2,tt,gen) = hfh(tt-1,gen)*Mf_high(tt,gen)*CB + hfh(tt,gen)*Mf_high(tt-1,gen)*SB
	GLR(5,1,2,tt,gen) = hfh(tt-1,gen)*(Af(tt-1,gen)*CB + MUE*SB) 
	GLR(5,2,1,tt,gen) = hfh(tt,gen)*(Af(tt,gen)*SB + MUE*CB)

	GLR(6,1,1,tt,gen) = -hfh(tt,gen)*Mf_high(tt,gen)*CB + hfh(tt-1,gen)*Mf_high(tt-1,gen)*SB + EL*MW/(sqrt2*SW)*S2B
	GLR(6,2,2,tt,gen) = hfh(tt-1,gen)*Mf_high(tt,gen)*SB - hfh(tt,gen)*Mf_high(tt-1,gen)*CB
	GLR(6,1,2,tt,gen) = hfh(tt-1,gen)*(Af(tt-1,gen)*SB - MUE*CB) 
	GLR(6,2,1,tt,gen) = hfh(tt,gen)*(-Af(tt,gen)*CB + MUE*SB)

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

!***************************************************************** 
! Vertex Gluino - Squark (outgoing) - Quark (ingoing) *
!*****************************************************************

	do i = 1, 2
		IL(i,tt,gen) = - GS*sqrt2*USf(i,1,tt,gen)
		IR(i,tt,gen) =   GS*sqrt2*USf(i,2,tt,gen)
	enddo
	
!***************************************************************** 
! Vertex Gluino - Squark (ingoing) - Quark (outgoing) *
!*****************************************************************

	do i = 1, 2
		JL(i,tt,gen) = dconjg(IR(i,tt,gen))
		JR(i,tt,gen) = dconjg(IL(i,tt,gen))
	enddo
		  
	
	
	
! End of loops over tt and gen
	enddo
	enddo
	
c	write(20,*),REAL(AR(3,1)*(BL(3,1)-BR(3,1))),REAL(AR(4,1)*(BL(4,1)-BR(4,1))),REAL(CR(1,3,1)*DL(1,3,1)),REAL(CR(1,4,1)*DL(1,4,1)),REAL(CR(2,3,1)*DL(2,3,1)),REAL(CR(2,4,1)*DL(2,4,1)),REAL(FL(1,1,3,1)*dconjg(FR(1,1,3,1))),REAL(FL(1,1,4,1)*dconjg(FR(1,1,4,1))),REAL(FR(1,1,3,1)*dconjg(FR(1,1,3,1)) + FL(1,1,3,1)*dconjg(FL(1,1,3,1))),REAL(FR(1,1,4,1)*dconjg(FR(1,1,4,1)) + FL(1,1,4,1)*dconjg(FL(1,1,4,1)))
c	close (unit=20)
	
	end

***************** Yukawa coupling at high scale ***********************

	double precision function hfh(type, gen)
	implicit none
#include "DD_Model.h"

	integer type, gen
	
	if ((type.eq.1).or.(type.eq.3)) then
	
  	 hfh = EL*Mf_high(type,gen)/(sqrt2*MW*SW*SB) 
  	 
	else
	
  	 hfh = EL*Mf_high(type,gen)/(sqrt2*MW*SW*CB) 
  	 
	endif	
	end
