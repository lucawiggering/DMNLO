!==================================================================================================
!======================= Calculates propagator contributions to four fermion couplings ============
!=======================				Author: Patrick Steppeler	   	               ============
!==================================================================================================

	subroutine DD_Prop
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_Model.h"
#include "DD_GenCouplings.h"
	
!----------------------- Declare variables --------------------------------------------------------
	complex*16 axialfactor(2,3), scalarfactor(2,3)
	integer nucleon, uptype, gen, tt, i, j
	double complex DDPropLoop
	logical DD_udsCheck
	
!======================= ACTUAL CALCULATIONS ======================================================
	
!======================= UNSUPPRESSED CONTRIBUTIONS ===============================================

!----------------------- Unsuppressed scalar contributions ----------------------------------------
	do gen = 1,3
		do uptype = 1,2
			scalarfactor(uptype,gen) = 0.0d0
		enddo
	enddo
	
!----------------------- Squark contribution s- & u-channel
	do gen = 1,3
		do uptype = 1,2
		
			tt = uptype + 2
			
			do i = 1,2
			do j = 1,2
				
				if (Denom_Mode==0) then
				
				scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 0.25d0*CF*FL(1,i,tt,gen)
     &										*DDPropLoop(s_Mandel(tt,gen),i,j,tt,gen)*dconjg(FR(1,j,tt,gen))/
     &										(MSf2(i,tt,gen)*MSf2(j,tt,gen))
	 
				scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 0.25d0*CF*FL(1,i,tt,gen)
     &										*DDPropLoop(u_Mandel(tt,gen),i,j,tt,gen)*dconjg(FR(1,j,tt,gen))/
     &										(MSf2(i,tt,gen)*MSf2(j,tt,gen))
				
				else
				
				scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 0.25d0*CF*FL(1,i,tt,gen)
     &										*DDPropLoop(s_Mandel(tt,gen),i,j,tt,gen)*dconjg(FR(1,j,tt,gen))/
     &										((MSf2(i,tt,gen)-s_Mandel(tt,gen))*(MSf2(j,tt,gen)-s_Mandel(tt,gen)))
	 
				scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 0.25d0*CF*FL(1,i,tt,gen)
     &										*DDPropLoop(u_Mandel(tt,gen),i,j,tt,gen)*dconjg(FR(1,j,tt,gen))/
     &										((MSf2(i,tt,gen)-u_Mandel(tt,gen))*(MSf2(j,tt,gen)-u_Mandel(tt,gen)))
				
				endif
				
			enddo
			enddo
		enddo
	enddo
	
!----------------------- Total unsuppressed scalar contribution	
	do gen = 1,3
!		gen=3
		do uptype = 1,2
			do nucleon = 1,2
				tt = uptype + 2
				g_4Fs(nucleon) = g_4Fs(nucleon) + f_T(nucleon,uptype,gen)*m_Nucleon(nucleon)
     &							 /Mf_low(tt,gen)*scalarfactor(uptype,gen)*Running_SI(tt,gen)
				g_4FsContribs(nucleon,4) = g_4FsContribs(nucleon,4) + f_T(nucleon,uptype,gen)*m_Nucleon(nucleon)
     &							 /Mf_low(tt,gen)*scalarfactor(uptype,gen)*Running_SI(tt,gen)
			enddo
		enddo
	enddo
	
!----------------------- Unsuppressed axial couplings ---------------------------------------------
	do gen = 1,3
		do uptype = 1,2
			axialfactor(uptype,gen) = 0.0d0
		enddo
	enddo

!Note: For simplicity we loop over all quark types in the following, but for spin dependent only u,d and s contribute. This is
!	   taken into account by the uds Check.
	
!----------------------- Squark contribution s- & u-channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
			do j = 1,2
	
			tt = uptype+2
			
			if(DD_udsCheck(tt,gen)) then
			
			if (Denom_Mode==0) then
		
				axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.125d0*CF*DDPropLoop(s_Mandel(tt,gen),i,j,tt,gen)*(FR(1,i,tt,gen)*
     &				dconjg(FR(1,j,tt,gen)) + FL(1,i,tt,gen)*dconjg(FL(1,j,tt,gen)))/(MSf2(i,tt,gen)*MSf2(j,tt,gen))
	 
				axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.125d0*CF*DDPropLoop(u_Mandel(tt,gen),i,j,tt,gen)*(FR(1,i,tt,gen)*
     &				dconjg(FR(1,j,tt,gen)) + FL(1,i,tt,gen)*dconjg(FL(1,j,tt,gen)))/(MSf2(i,tt,gen)*MSf2(j,tt,gen))
		
			else
		
				axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.125d0*CF*DDPropLoop(s_Mandel(tt,gen),i,j,tt,gen)*(FR(1,i,tt,gen)*
     &				dconjg(FR(1,j,tt,gen)) + FL(1,i,tt,gen)*dconjg(FL(1,j,tt,gen)))/((-s_Mandel(tt,gen)+MSf2(i,tt,gen))*(MSf2(j,tt,gen)-s_Mandel(tt,gen)))
	 
				axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.125d0*CF*DDPropLoop(u_Mandel(tt,gen),i,j,tt,gen)*(FR(1,i,tt,gen)*
     &				dconjg(FR(1,j,tt,gen)) + FL(1,i,tt,gen)*dconjg(FL(1,j,tt,gen)))/((-u_Mandel(tt,gen)+MSf2(i,tt,gen))*(MSf2(j,tt,gen)-u_Mandel(tt,gen)))
	 
			endif
			
			endif
			
			enddo
			enddo
		enddo
	enddo
		
!----------------------- Total unsuppressed axial contribution
	do gen = 1,3
!	gen=3
		do uptype = 1,2
			do nucleon = 1,2
				g_4Fa(nucleon) = g_4Fa(nucleon) + Delta_q(nucleon,uptype,gen)
     &							 *axialfactor(uptype,gen)*Running_SD(tt,gen)
				g_4FaContribs(nucleon,4) = g_4FaContribs(nucleon,4) + Delta_q(nucleon,uptype,gen)
     &							 *axialfactor(uptype,gen)*Running_SD(tt,gen)
			enddo
		enddo
	enddo


! ====================== SUPPRESSED CONTRIBUTIONS =================================================

	if(Supp_Mode==1.OR.Supp_Mode==2) then

!----------------------- Suppressed scalar contributions ------------------------------------------

!----------------------- Total suppressed scalar contribution

!----------------------- Suppressed axial contributions -------------------------------------------
		
!----------------------- Total suppressed axial contribution

	endif				! End of Supp_Mode if
	
! ====================== STRONGLY SUPPRESSED CONTRIBUTIONS ========================================

	if(Supp_Mode==2) then

!----------------------- Strongly suppressed scalar contributions ---------------------------------

!----------------------- Total strongly suppressed scalar contribution

!----------------------- Strongly suppressed axial contributions ----------------------------------
		
!----------------------- Total strongly suppressed axial contribution

	endif				! End of Supp_Mode if
	
!======================= END OF PROGRAM ===========================================================
	
	g_4Fsp = g_4Fs(1)
	g_4Fsn = g_4Fs(2)
	g_4Fap = g_4Fa(1)
	g_4Fan = g_4Fa(2)
	
	end



!==================================================================================================
!======================= Calculates the actual propagator corrections to squark propagators =======
!=======================			Copied from ChiChi2QQ_M2Propagator.F					=======
!==================================================================================================

	function DDPropLoop(p2,k,l,tt,gen)

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_Model.h"
#include "DD_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"
	
!----------------------- Declare variables --------------------------------------------------------
	
	double precision p2
	integer tt, gen, m, k, l
	double complex g1L, g1R, g2L, g2R, ML0, ML1, g1, g2
	double precision A01, Krondelta
	double complex B02
	double complex DDPropLoop

!======================= ACTUAL CALCULATIONS ======================================================
	
	DDPropLoop = 0d0

! *** Gluon Loop ***
! Note: gluon one-point diagram is zero due to A0(mg=0)=0

	g1 = -GS*Krondelta(k,l)
	g2 = -GS*Krondelta(k,l)
	
	ML0 = 0d0 
	ML1 = MSf(k,tt,gen)

	DDPropLoop = DDPropLoop - 
     &		(g1*g2*(2*A01(ML0) - A01(ML1) + (-ML0**2 + 2*ML1**2 + 2*p2)*B02(p2,ML0,ML1)))/(16.*Pi**2)
		
! *** Gluino Loop ***
	g1L = sqrt2*GS*USf(k,2,tt,gen)
	g1R = -sqrt2*GS*USf(k,1,tt,gen)
	g2L = -sqrt2*GS*USf(l,1,tt,gen)
	g2R = sqrt2*GS*USf(l,2,tt,gen)
      
	ML0 = Mf(tt,gen) 
	ML1 = MGl
	
	DDPropLoop = DDPropLoop - 
     &		(2*(g1L*g2L + g1R*g2R)*ML0*ML1*B02(p2,ML0,ML1) + 
     &     	(g1R*g2L + g1L*g2R)*(A01(ML0) + A01(ML1) + (ML0**2 + ML1**2 - p2)*B02(p2,ML0,ML1)))/(16.*Pi**2)

! *** 4 Squark Loop ***
	do m=1,2
      
		g1 = -GS**2*(USf(m,1,tt,gen)*USf(l,1,tt,gen)-USf(m,2,tt,gen)*USf(l,2,tt,gen))*
     &              (USf(m,1,tt,gen)*USf(k,1,tt,gen)-USf(m,2,tt,gen)*USf(k,2,tt,gen))
      
		ML0 = MSf(m,tt,gen)
	
		DDPropLoop = DDPropLoop - (g1*A01(ML0))/(16.*Pi**2)	
	
	enddo
	
	end


