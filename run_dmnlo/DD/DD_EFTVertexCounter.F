!==================================================================================================
!========= Here the Gluin three-point vertex from the effective field theory side gets ============
!========= its counterterm. This cancels its UV divergence and the IR divegences of    ============
!=========		the quark wave function renormalization dZf appearing for V5-V8.	   ============
!=========							Author: Patrick Steppeler						   ============
!==================================================================================================

	subroutine DD_EFTVertexcounter	
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_Model.h"
#include "DD_GenCouplings.h"
#include "DD_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"
	
!----------------------- Declare variables --------------------------------------------------------
	complex*16 axialfactor(2,3), scalarfactor(2,3), TreeCoupling
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

!----------------------- EFT Vertex s-channel
	do gen = 1,3
		do uptype = 1,2
			tt = uptype + 2
			do i = 1,2
						
			TreeCoupling = - 0.25d0*FL(1,i,tt,gen)*dconjg(FR(1,i,tt,gen))/(MSf(i,tt,gen)**2-s_Mandel(tt,gen))
			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 0.5d0*TreeCoupling*(dEFT_SI_L(1,i,tt,gen) + dconjg(dEFT_SI_R(1,i,tt,gen)))

			enddo
		enddo
	enddo
	
!----------------------- EFT Vertex u-channel
	do gen = 1,3
		do uptype = 1,2
			tt = uptype + 2
			do i = 1,2

			TreeCoupling = - 0.25d0*FL(1,i,tt,gen)*dconjg(FR(1,i,tt,gen))/(MSf(i,tt,gen)**2-u_Mandel(tt,gen))
			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 0.5d0*TreeCoupling*(dEFT_SI_L(1,i,tt,gen) + dconjg(dEFT_SI_R(1,i,tt,gen)))
	 
			enddo
		enddo
	enddo
	
!----------------------- Total unsuppressed scalar contribution
	do gen = 1,3
!	gen = 3
		do uptype = 1,2
			do nucleon = 1,2
				tt = uptype + 2
				g_4Fs(nucleon) = g_4Fs(nucleon) + f_T(nucleon,uptype,gen)*m_Nucleon(nucleon)/Mf_low(tt,gen)*scalarfactor(uptype,gen)*Running_SI(tt,gen)
				g_4FsContribs(nucleon,9) = g_4FsContribs(nucleon,9) + f_T(nucleon,uptype,gen)*m_Nucleon(nucleon)/Mf_low(tt,gen)*scalarfactor(uptype,gen)*Running_SI(tt,gen)
			enddo
		enddo
	enddo

!----------------------- Unsuppressed axial contributions -----------------------------------------
	do gen = 1,3
		do uptype = 1,2
			axialfactor(uptype,gen) = 0.0d0
		enddo
	enddo

!Note: For simplicity we loop over all quark types in the following, but for spin dependent only u,d and s contribute. This is
!	   taken into account by the uds Check.

!----------------------- EFT Vertex s-channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
			tt = uptype+2
			
			if(DD_udsCheck(tt,gen)) then
			
			TreeCoupling = 0.125d0*(FR(1,i,tt,gen)*dconjg(FR(1,i,tt,gen)) + FL(1,i,tt,gen)*dconjg(FL(1,i,tt,gen)))/(MSf2(i,tt,gen)-s_Mandel(tt,gen))
			axialfactor(uptype,gen) = axialfactor(uptype,gen) - 0.5d0*TreeCoupling*(dEFT_SD_L(1,i,tt,gen) + dEFT_SD_R(1,i,tt,gen))
			
			endif
			
			enddo
		enddo
	enddo
	
!----------------------- EFT Vertex u-channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
			tt = uptype+2
			
			if(DD_udsCheck(tt,gen)) then
			
			TreeCoupling = 0.125d0*(FR(1,i,tt,gen)*dconjg(FR(1,i,tt,gen)) + FL(1,i,tt,gen)*dconjg(FL(1,i,tt,gen)))/(MSf2(i,tt,gen)-u_Mandel(tt,gen))
			axialfactor(uptype,gen) = axialfactor(uptype,gen) - 0.5d0*TreeCoupling*(dEFT_SD_L(1,i,tt,gen) + dEFT_SD_R(1,i,tt,gen))
				
			endif
			
			enddo
		enddo
	enddo
	
!----------------------- Total unsuppressed axial contribution
	do gen = 1,3
!		gen = 3
		do uptype = 1,2
			do nucleon = 1,2
				g_4Fa(nucleon) = g_4Fa(nucleon) + Delta_q(nucleon,uptype,gen)*axialfactor(uptype,gen)*Running_SD(tt,gen)
				g_4FaContribs(nucleon,9) = g_4FaContribs(nucleon,9) + Delta_q(nucleon,uptype,gen)*axialfactor(uptype,gen)*Running_SD(tt,gen)
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


