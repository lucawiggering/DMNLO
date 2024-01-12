!==================================================================================================
!======================= Calculates effective couplings, like Delta m_b, Delta_QCD etc. ===========
!=======================		Copy and expansion of ChiChi2QQ_Couplings.F				===========
!=======================				Author: Patrick Steppeler						===========
!==================================================================================================

	subroutine DD_SetCouplingsHO
	implicit none

#include "DD_Model.h"
#include "DD_GenCouplings.h"
#include "DD_GenCounterterms.h"

	double precision DmbResults(5), DeltaTlighth, DeltaTheavyH, DeltaTA
	integer i,j,k,l,a
	integer icha, indi, indj, tt, gen
	
! functions from DMNLO_Kinematics.F
	double precision Krondelta
	
! default settings for direct detection
	icha = 0
	indi = 1
	indj = 1

! Loop over tt and gen
	do tt = 3,4
	do gen = 1,3	
	
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
	NNLO_Res=1		!(0 = Old DMNLO variant, 1 = Newer Spira variant)
	Res_Style=1		!Switch the NNLO Resummation on and off (1 = on, 0 = off). The use NNLO_Res = 1 sets Res_Style to 1 automatically.
	call DMNLO_CalcDeltaMb(DmbResults,gen)

cc	DeltaMb = 0d0
cc	DeltaMbQCD = 0d0
	
!	write(*,*)'tt,gen',tt,gen
!	write(*,*)'QCD NLO',DmbResults(1)
!	write(*,*)'QCD NNLO',DmbResults(2)
!	write(*,*)'ELW NLO',DmbResults(3)
!	write(*,*)'ELW NNLO',DmbResults(4)
!	write(*,*)'Full Delta m_b',DmbResults(5)
!	write(*,*)'Delta QCD',DeltaQCD
!	write(*,*)'DeltaTlighth',DeltaTlighth
!	write(*,*)'DeltaTheavyH',DeltaTheavyH
!	write(*,*)'-------------------------'

! Default Setting is zero
	do i=1,6
	  	DLho(i,tt,gen) = 0d0
	  	DRho(i,tt,gen) = 0d0
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
	  	  
! Note Patrick: Uncomment to apply resummation only for generation three. This is not recommandable in the case of direct detection as the small 
! masses in the Yukawa couplings of the first generation quarks cancel out in the end. In addition there are no (sizeable) electroweak contributions
! to Delta m_b for the first to generations as these are related to the top Yukawa coupling. The electroweak contributions tend to reduce Delta m_b
! normally. The lines below are more general than needed and include also the case icha = 1,2 which is not needed for direct detection.

! On the other hand using the resummation for the first generation quarks requires a consistent renormalization scheme including RGEs for the light
! quarks. As this is currently not available in DM@NLO, we choose a conservative approach and activate the resummation only for the third generation.
  
	if (gen.eq.3)	then
	
! Down-Type Quarks
		if (tt.eq.4)	then
	 
		DLho(1,tt,gen) = DL(1,tt,gen)*dsqrt(DeltaQCD + DeltaTlighth)*(1d0/(1d0+DeltaMb)*(1d0-DeltaMb*CA/(SA*TB))) 
     &	  	    + DL(1,tt,gen)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD + DeltaMbQCD*CA/(SA*TB))
		DRho(1,tt,gen) = DR(1,tt,gen)*dsqrt(DeltaQCD + DeltaTlighth)*(1d0/(1d0+DeltaMb)*(1d0-DeltaMb*CA/(SA*TB))) 
     &	  	    + DR(1,tt,gen)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD + DeltaMbQCD*CA/(SA*TB))

		DLho(2,tt,gen) = DL(2,tt,gen)*dsqrt(DeltaQCD + DeltaTheavyH)*(1d0/(1d0+DeltaMb)*(1d0+DeltaMb*SA/(CA*TB))) 
     &	  	    + DL(2,tt,gen)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD - DeltaMbQCD*SA/(CA*TB))
		DRho(2,tt,gen) = DR(2,tt,gen)*dsqrt(DeltaQCD + DeltaTheavyH)*(1d0/(1d0+DeltaMb)*(1d0+DeltaMb*SA/(CA*TB))) 
     &	  	    + DR(2,tt,gen)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD - DeltaMbQCD*SA/(CA*TB))

		DLho(3,tt,gen) = DL(3,tt,gen)*dsqrt(DeltaQCD + DeltaTA)*(1d0/(1d0 + DeltaMb)*(1d0 - DeltaMb/(TB**2))) 
     &	  	    + DL(3,tt,gen)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD + DeltaMbQCD/(TB**2))
		DRho(3,tt,gen) = DR(3,tt,gen)*dsqrt(DeltaQCD + DeltaTA)*(1d0/(1d0 + DeltaMb)*(1d0 - DeltaMb/(TB**2))) 
     &	  	    + DR(3,tt,gen)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD + DeltaMbQCD/(TB**2))
	
		DLho(4,tt,gen) = DL(4,tt,gen)*dsqrt(DeltaQCD) + DL(4,tt,gen)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0)
		DRho(4,tt,gen) = DR(4,tt,gen)*dsqrt(DeltaQCD) + DR(4,tt,gen)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0)

		if (icha.eq.0)	then
		
		do i = 1, 2
			do a = 1, 4
		
			ELcho(a,i,tt,gen) = - EL*Mf_high(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,1,tt,gen)
			ERho(a,i,tt,gen) =  - EL*Mf_high(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,2,tt,gen) 
			FLho(a,i,tt,gen) =  - EL*Mf_high(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,2,tt,gen)
			FRho(a,i,tt,gen) =  - EL*Mf_high(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,1,tt,gen)

			GLho(a,i,tt,gen) =  - EL*Mf_high(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,2,tt,gen)
			GRho(a,i,tt,gen) =  - EL*Mf_high(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,1,tt,gen)
			HLho(a,i,tt,gen) =  - EL*Mf_high(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,1,tt,gen)
			HRho(a,i,tt,gen) =  - EL*Mf_high(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,2,tt,gen)
	 
			enddo
		enddo
		
		else if (icha.eq.2)	then
		
		do i = 1, 2
			do a = 1, 2
			
			GRho(a,i,tt,gen) =  EL*Mf_high(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*UCha(a,2)*USf(i,1,tt-1,gen)
			HLho(a,i,tt,gen) =  EL*Mf_high(tt,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*UChaC(a,2)*USf(i,1,tt-1,gen)
	  
			enddo
		enddo
		
		endif
		
! Up-Type Quarks
		else if (tt.eq.3) then
		
		do i=1, 4
			
			DLho(i,tt,gen) = DL(i,tt,gen)*dsqrt(DeltaQCD-1d0-(Alfas/pi)*5.67d0)
			DRho(i,tt,gen) = DR(i,tt,gen)*dsqrt(DeltaQCD-1d0-(Alfas/pi)*5.67d0)
			
		enddo
			
		if (icha.eq.1)	then
		
			DLho(5,tt,gen) = DL(5,tt,gen)*dsqrt(DeltaQCD)*(1d0/(1d0+DeltaMb)) 
     &	  	    + DL(5,tt,gen)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD)
			DRho(5,tt,gen) = DR(5,tt,gen)*dsqrt(DeltaQCD)*(1d0/(1d0+DeltaMb)) 
     &	  	    + DR(5,tt,gen)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD)

			DLho(6,tt,gen) = DL(6,tt,gen)*dsqrt(DeltaQCD)*(1d0/(1d0+DeltaMb)) 
     &	  	    + DL(6,tt,gen)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD)
			DRho(6,tt,gen) = DR(6,tt,gen)*dsqrt(DeltaQCD)*(1d0/(1d0+DeltaMb)) 
     &	  	    + DR(6,tt,gen)*( -1d0 - 0.5d0*(Alfas/pi)*5.67d0 + DeltaMbQCD)
	 
		do i = 1, 2
			do a = 1, 4
		
			FLho(a,i,tt,gen) =  -EL*Mf_high(tt+1,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeuC(a,3)*USf(i,2,tt+1,gen)
			FRho(a,i,tt,gen) =  -EL*Mf_high(tt+1,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*ZNeu(a,3)*USf(i,1,tt+1,gen)
	 
			enddo
			
			do a = 1, 2
			
			GRho(a,i,tt,gen) =  EL*Mf_high(tt+1,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*UCha(a,2)*USf(i,1,tt,gen)
			ERho(a,i,tt,gen) =  EL*Mf_high(tt+1,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*UCha(a,2)*USf(i,2,tt+1,gen)
			
			enddo
		enddo
		
		else if (icha.eq.2)	then
		
		do i = 1, 2
			do a = 1, 2
			
			ERho(a,i,tt,gen) = EL*Mf_high(tt+1,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*UCha(a,2)*USf(i,2,tt+1,gen) 
			FLho(a,i,tt,gen) = EL*Mf_high(tt+1,gen)*(1d0/(1d0 + DeltaMb) - 1d0 + DeltaMbQCD)/(SW*sqrt2*MW*CB)*UChaC(a,2)*USf(i,2,tt+1,gen)
			
			enddo
		enddo
	  
		endif
	 endif
	 	
	endif		  
	
	
!	write(*,*)'tt,gen',tt,gen
!	write(*,*)'--------------'
!	write(*,*)'DLho(1,tt,gen)',DLho(1,tt,gen)
!	write(*,*)'DRho(1,tt,gen)',DRho(1,tt,gen)
!	write(*,*)'DLho(2,tt,gen)',DLho(2,tt,gen)
!	write(*,*)'DRho(2,tt,gen)',DRho(2,tt,gen)
!	write(*,*)'--------------'
!	write(*,*)'DL(1,tt,gen)',DL(1,tt,gen)
!	write(*,*)'DR(1,tt,gen)',DR(1,tt,gen)
!	write(*,*)'DL(2,tt,gen)',DL(2,tt,gen)
!	write(*,*)'DR(2,tt,gen)',DR(2,tt,gen)
!	write(*,*)'--------------'
!	write(*,*)'ELcho(1,1,tt,gen)',ELcho(1,1,tt,gen)
!	write(*,*)'ERho(1,1,tt,gen)',ERho(1,1,tt,gen)
!	write(*,*)'FLho(1,1,tt,gen)',FLho(1,1,tt,gen)
!	write(*,*)'FRho(1,1,tt,gen)',FRho(1,1,tt,gen)
!	write(*,*)'--------------'
!	write(*,*)'ELc(1,1,tt,gen)',ELc(1,1,tt,gen)
!	write(*,*)'ER(1,1,tt,gen)',ER(1,1,tt,gen)
!	write(*,*)'FL(1,1,tt,gen)',FL(1,1,tt,gen)
!	write(*,*)'FR(1,1,tt,gen)',FR(1,1,tt,gen)
!	write(*,*)'--------------'
!	write(*,*)'GLho(1,1,tt,gen)',GLho(1,1,tt,gen)
!	write(*,*)'GRho(1,1,tt,gen)',GRho(1,1,tt,gen)
!	write(*,*)'HLho(1,1,tt,gen)',HLho(1,1,tt,gen)
!	write(*,*)'HRho(1,1,tt,gen)',HRho(1,1,tt,gen)
!	write(*,*)'--------------'
!	write(*,*)'GL(1,1,tt,gen)',GL(1,1,tt,gen)
!	write(*,*)'GR(1,1,tt,gen)',GR(1,1,tt,gen)
!	write(*,*)'HL(1,1,tt,gen)',HL(1,1,tt,gen)
!	write(*,*)'HR(1,1,tt,gen)',HR(1,1,tt,gen)
!	write(*,*)'--------------'
!	write(*,*)''
	
! End of loops over tt and gen
	enddo
	enddo
	
	end
