!==================================================================================================
!======================= Calculates propagator counterterms to four fermion couplings =============
!=======================				Author: Patrick Steppeler	   	              =============
!==================================================================================================

	subroutine DD_PropCounter	
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_Model.h"
#include "DD_GenCouplings.h"
#include "DD_GenCounterterms.h"
	
!----------------------- Declare variables --------------------------------------------------------
	complex*16 axialfactor(2,3), scalarfactor(2,3)
	integer nucleon, uptype, gen, tt, i, j
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
				
				call DD_CalcProp(s_Mandel(tt,gen),tt,gen)
				scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 0.25d0*FL(1,i,tt,gen)
     &										*DDPropCounter(i,j,tt,gen)*dconjg(FR(1,j,tt,gen))/
     &										(MSf2(i,tt,gen)*MSf2(j,tt,gen))
				
				call DD_CalcProp(u_Mandel(tt,gen),tt,gen)
				scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 0.25d0*FL(1,i,tt,gen)
     &										*DDPropCounter(i,j,tt,gen)*dconjg(FR(1,j,tt,gen))/
     &										(MSf2(i,tt,gen)*MSf2(j,tt,gen))
				else
				
				call DD_CalcProp(s_Mandel(tt,gen),tt,gen)
				scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 0.25d0*FL(1,i,tt,gen)
     &										*DDPropCounter(i,j,tt,gen)*dconjg(FR(1,j,tt,gen))/
     &										((-s_Mandel(tt,gen)+MSf2(i,tt,gen))*(MSf2(j,tt,gen)-s_Mandel(tt,gen)))
				
				call DD_CalcProp(u_Mandel(tt,gen),tt,gen)
				scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 0.25d0*FL(1,i,tt,gen)
     &										*DDPropCounter(i,j,tt,gen)*dconjg(FR(1,j,tt,gen))/
     &										((-u_Mandel(tt,gen)+MSf2(i,tt,gen))*(MSf2(j,tt,gen)-u_Mandel(tt,gen)))
				endif
	 
			enddo
			enddo
		enddo
	enddo
	
!----------------------- Total unsuppressed scalar contribution	
	do gen = 1,3
!	gen = 3
		do uptype = 1,2
			do nucleon = 1,2
				tt = uptype + 2
				g_4Fs(nucleon) = g_4Fs(nucleon) + f_T(nucleon,uptype,gen)*m_Nucleon(nucleon)
     &							 /Mf_low(tt,gen)*scalarfactor(uptype,gen)*Running_SI(tt,gen)
				g_4FsContribs(nucleon,5) = g_4FsContribs(nucleon,5) + f_T(nucleon,uptype,gen)*m_Nucleon(nucleon)
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
			
			if(Denom_Mode==0) then
		
				call DD_CalcProp(s_Mandel(tt,gen),tt,gen)
				axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.125d0*DDPropCounter(i,j,tt,gen)*(FR(1,i,tt,gen)*
     &				dconjg(FR(1,j,tt,gen)) + FL(1,i,tt,gen)*dconjg(FL(1,j,tt,gen)))/(MSf2(i,tt,gen)*MSf2(j,tt,gen))
		
				call DD_CalcProp(u_Mandel(tt,gen),tt,gen)
				axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.125d0*DDPropCounter(i,j,tt,gen)*(FR(1,i,tt,gen)*
     &				dconjg(FR(1,j,tt,gen)) + FL(1,i,tt,gen)*dconjg(FL(1,j,tt,gen)))/(MSf2(i,tt,gen)*MSf2(j,tt,gen))
		
			else
		
				call DD_CalcProp(s_Mandel(tt,gen),tt,gen)
				axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.125d0*DDPropCounter(i,j,tt,gen)*(FR(1,i,tt,gen)*
     &				dconjg(FR(1,j,tt,gen)) + FL(1,i,tt,gen)*dconjg(FL(1,j,tt,gen)))/((-s_Mandel(tt,gen)+MSf2(i,tt,gen))*(MSf2(j,tt,gen)-s_Mandel(tt,gen)))
		
				call DD_CalcProp(u_Mandel(tt,gen),tt,gen)
				axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.125d0*DDPropCounter(i,j,tt,gen)*(FR(1,i,tt,gen)*
     &				dconjg(FR(1,j,tt,gen)) + FL(1,i,tt,gen)*dconjg(FL(1,j,tt,gen)))/((-u_Mandel(tt,gen)+MSf2(i,tt,gen))*(MSf2(j,tt,gen)-u_Mandel(tt,gen)))
	 
			endif
			
			endif
	 
			enddo
			enddo
		enddo
	enddo
	
!----------------------- Total unsuppressed axial contribution
	do gen = 1,3
!	gen = 3
		do uptype = 1,2
			do nucleon = 1,2
				g_4Fa(nucleon) = g_4Fa(nucleon) + Delta_q(nucleon,uptype,gen)
     &							 *axialfactor(uptype,gen)*Running_SD(tt,gen)
				g_4FaContribs(nucleon,5) = g_4FaContribs(nucleon,5) + Delta_q(nucleon,uptype,gen)
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
!======================= Calculates actual propagator counterterms  ===============================
!=======================	Copied from ChiChi2QQ_SetCounterms.F	===============================
!==================================================================================================

	subroutine DD_CalcProp(myP2,t,g)
	implicit none
	
!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_Model.h"
#include "DD_GenCouplings.h"
#include "DD_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"

!----------------------- Declare variables --------------------------------------------------------
	integer t,g,i,j
	double precision myP2, Krondelta

!======================= ACTUAL CALCULATIONS ======================================================

	 do i = 1, 2
	  do j = 1, 2

	  DDPropCounter(i,j,t,g) = 1/2d0*dZSf(i,j,t,g)*(myP2 - MSf2(i,t,g)) + 1/2d0*dZSf(j,i,t,g)*(myP2 - MSf2(j,t,g)) - Krondelta(i,j)*dMSf2(i,t,g)
	  
	  enddo
 	 enddo

	end


