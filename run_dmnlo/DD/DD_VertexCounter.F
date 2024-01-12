!==================================================================================================
!======================= Calculates vertex counterterms to four fermion couplings =================
!=======================				Author: Patrick Steppeler	   	          =================
!==================================================================================================

	subroutine DD_VertexCounter
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_Model.h"
#include "DD_GenCouplings.h"
#include "DD_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"

!----------------------- Declare variables --------------------------------------------------------
	complex*16 ALV1, ARV1, ALV2, ARV2, ALV5, ARV5, ALV6, ARV6, ALV7, ARV7, ALV8, ARV8
	complex*16 AL1V3, AR1V3, AL2V3, AR2V3, AL3V3, AR3V3 
	complex*16 AL1V4, AR1V4, AL2V4, AR2V4, AL3V4, AR3V4
	complex*16 axialfactor(2,3), scalarfactor(2,3)
	integer nucleon, uptype, gen, tt, i
	logical DD_udsCheck
	
!======================= ACTUAL CALCULATIONS ======================================================
	
!======================= UNSUPPRESSED CONTRIBUTIONS ===============================================

!----------------------- Unsuppressed scalar contributions ----------------------------------------
	do gen = 1,3
		do uptype = 1,2
			scalarfactor(uptype,gen) = 0.0d0
		enddo
	enddo
	
!----------------------- Higgs contribution (h0 & H0) V1 & V2: Gluon & Gluino Vertex, t-Channel
!	do gen = 1,3
!!	gen = 3
!		do uptype = 1,2
!			tt = uptype + 2
!				
!			do i = 1,2
!				
!				ALV1 = dDL(i,tt,gen)
!				ARV1 = dDR(i,tt,gen)
!				
!				if (Denom_Mode==0) then
!				
!				scalarfactor(uptype,gen) = scalarfactor(uptype,gen) + CR(i,tt,gen)
!     &										*(ALV1 + ARV1)/(2.0d0*MHk(i)**2)
!				
!				else
!				
!				scalarfactor(uptype,gen) = scalarfactor(uptype,gen) + CR(i,tt,gen)
!     &										*(ALV1 + ARV1)/(2.0d0*(MHk(i)**2-t_Mandel(tt,gen)))
!	 
!				endif
!				
!			enddo
!		enddo
!	enddo
		
!----------------------- Squark contribution V5 & V6: Gluon & Gluino Vertex, Outgoing Squark, s-Channel
	do gen = 1,3
		do uptype = 1,2
			tt = uptype + 2
	 
			do i = 1,2
			
			ALV5 = dFL(1,i,tt,gen)
			ARV5 = dFR(1,i,tt,gen)
			
			if (Denom_Mode==0) then
			
			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV5 + ALV5)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) + 
     &			(ARV5 - ALV5)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/MSf(i,tt,gen)**2
			
			else
			
			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV5 + ALV5)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) + 
     &			(ARV5 - ALV5)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/(MSf(i,tt,gen)**2-s_Mandel(tt,gen))
	 
			endif
	 
			enddo
		enddo
	enddo

!----------------------- Squark contribution V7 & V8: Gluon & Gluino Vertex, Ingoing Squark, s-Channel
	do gen = 1,3
		do uptype = 1,2
			tt = uptype + 2
	 
			do i = 1,2
			
			ALV7 = dELc(1,i,tt,gen)
			ARV7 = dER(1,i,tt,gen)
			
			if (Denom_Mode==0) then
			
			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV7 + ALV7)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) + 
     &			(ARV7 - ALV7)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/MSf(i,tt,gen)**2
			
			else
			
			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV7 + ALV7)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) + 
     &			(ARV7 - ALV7)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/(MSf(i,tt,gen)**2-s_Mandel(tt,gen))
	 
			endif
	 
			enddo
		enddo
	enddo
	
!----------------------- Squark contribution V5 & V6: Gluon & Gluino Vertex, Outgoing Squark, u-Channel
	do gen = 1,3
		do uptype = 1,2
			tt = uptype + 2
			
			do i = 1,2
			
			ALV5 = dFL(1,i,tt,gen)
			ARV5 = dFR(1,i,tt,gen)
			
			if (Denom_Mode==0) then
			
			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV5 + ALV5)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) + 
     &			(ARV5 - ALV5)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/MSf(i,tt,gen)**2
			
			else
			
			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV5 + ALV5)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) + 
     &			(ARV5 - ALV5)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/(MSf(i,tt,gen)**2-u_Mandel(tt,gen))
	 
			endif
	 
			enddo
		enddo
	enddo

!----------------------- Squark contribution V7 & V8: Gluon & Gluino Vertex, Ingoing Squark, u-Channel
	do gen = 1,3
		do uptype = 1,2
			tt = uptype + 2
			
			do i = 1,2
			
			ALV7 = dELc(1,i,tt,gen)
			ARV7 = dER(1,i,tt,gen)
			
			if (Denom_Mode==0) then
			
			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV7 + ALV7)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) + 
     &			(ARV7 - ALV7)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/MSf(i,tt,gen)**2
			
			else
			
			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV7 + ALV7)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) + 
     &			(ARV7 - ALV7)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/(MSf(i,tt,gen)**2-u_Mandel(tt,gen))
	 
			endif
	 
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
				g_4FsContribs(nucleon,3) = g_4FsContribs(nucleon,3) + f_T(nucleon,uptype,gen)*m_Nucleon(nucleon)
     &							 /Mf_low(tt,gen)*scalarfactor(uptype,gen)*Running_SI(tt,gen)
	 
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
	
!----------------------- Z contribution V3 & V4: Gluon & Gluino Vertex, t-Channel
	do gen = 1,3
		do uptype = 1,2
		tt = uptype+2
	
		if(DD_udsCheck(tt,gen)) then
		
			AL1V3 = dBL(tt,gen)
			AR1V3 = dBR(tt,gen)
	
		if (Denom_Mode==0) then
	
			axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.5d0*AR(tt,gen)*(AL1V3-AR1V3)/MZ2
	
		else
	
			axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.5d0*AR(tt,gen)*(AL1V3-AR1V3)/(MZ2-t_Mandel(tt,gen))
	
		endif
		
		endif
		
		enddo
	enddo
	
!----------------------- Squark contribution V5 & V6: Gluon & Gluino Vertex, Outgoing Squark, s-Channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
			tt = uptype + 2
			
			if(DD_udsCheck(tt,gen)) then
			
				ALV5 = dFL(1,i,tt,gen)
				ARV5 = dFR(1,i,tt,gen)
		
			if (Denom_Mode==0) then
		
				axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &				(-(ARV5 + ALV5)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) + 
     &				(ARV5 - ALV5)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/MSf(i,tt,gen)**2
		
			else
		
				axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &				(-(ARV5 + ALV5)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) + 
     &				(ARV5 - ALV5)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/(MSf(i,tt,gen)**2-s_Mandel(tt,gen))
	 
			endif
	 
			endif
			
			enddo
		enddo
	enddo
		
!----------------------- Squark contribution V7 & V8: Gluon & Gluino Vertex, Ingoing Squark, s-Channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
			tt = uptype + 2
			
			if(DD_udsCheck(tt,gen)) then

				ALV7 = dELc(1,i,tt,gen)
				ARV7 = dER(1,i,tt,gen)
		
			if (Denom_Mode==0) then
		
				axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &				(-(ARV7 + ALV7)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) + 
     &				(ARV7 - ALV7)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/MSf(i,tt,gen)**2
		
			else
		
				axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &				(-(ARV7 + ALV7)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) + 
     &				(ARV7 - ALV7)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/(MSf(i,tt,gen)**2-s_Mandel(tt,gen))
	 
			endif
			
			endif
			
			enddo
		enddo
	enddo

!----------------------- Squark contribution V5 & V6: Gluon & Gluino Vertex, Outgoing Squark, u-Channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
			tt = uptype + 2
			
			if(DD_udsCheck(tt,gen)) then

				ALV5 = dFL(1,i,tt,gen)
				ARV5 = dFR(1,i,tt,gen)
		
			if (Denom_Mode==0) then
		
				axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &				(-(ARV5 + ALV5)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) + 
     &				(ARV5 - ALV5)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/MSf(i,tt,gen)**2
		
			else
		
				axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &				(-(ARV5 + ALV5)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) + 
     &				(ARV5 - ALV5)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/(MSf(i,tt,gen)**2-u_Mandel(tt,gen))
	 
			endif
		
			endif
			
			enddo
		enddo
	enddo

!----------------------- Squark contribution V7 & V8: Gluon & Gluino Vertex, Ingoing Squark, u-Channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
			tt = uptype + 2
			
			if(DD_udsCheck(tt,gen)) then

				ALV7 = dELc(1,i,tt,gen)
				ARV7 = dER(1,i,tt,gen)
			
			if (Denom_Mode==0) then
		
			axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &			(-(ARV7 + ALV7)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) + 
     &			(ARV7 - ALV7)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/MSf(i,tt,gen)**2
		
			else
		
			axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &			(-(ARV7 + ALV7)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) + 
     &			(ARV7 - ALV7)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/(MSf(i,tt,gen)**2-u_Mandel(tt,gen))
	 
			endif
		
			endif
			
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
				g_4FaContribs(nucleon,3) = g_4FaContribs(nucleon,3) + Delta_q(nucleon,uptype,gen)
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


