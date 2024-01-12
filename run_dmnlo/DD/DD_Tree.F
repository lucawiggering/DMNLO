!==================================================================================================
!======================= Calculates tree level contributions to four fermion couplings ============
!=======================				Author: Patrick Steppeler	   	               ============
!==================================================================================================

	subroutine DD_Tree
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_Model.h"
#include "DD_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	
!----------------------- Declare variables --------------------------------------------------------
	complex*16 axialfactor(2,3), scalarfactor(2,3), Higgsbeitrag(2,3), Zbeitrag(2,3), Squarkbeitrag(2,3)
	integer nucleon, uptype, gen, tt, i
	logical DD_udsCheck
	
!======================= ACTUAL CALCULATIONS ======================================================
! Tree-Level is deactivated for IR Checks

	if(EpsPole.ge.0d0) then
	
!======================= UNSUPPRESSED CONTRIBUTIONS ===============================================

!----------------------- Unsuppressed scalar contributions ----------------------------------------
	do gen = 1,3
		do uptype = 1,2
			scalarfactor(uptype,gen) = 0.0d0
		enddo
	enddo
	
!----------------------- T1 Higgs contribution (h0 & H0) t-Channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
				tt = uptype + 2

				if (Denom_Mode==0) then

					scalarfactor(uptype,gen) = scalarfactor(uptype,gen) + CR(i,tt,gen)
     &										*DL(i,tt,gen)/MHk(i)**2

				else

					scalarfactor(uptype,gen) = scalarfactor(uptype,gen) + CR(i,tt,gen)
     &										*DL(i,tt,gen)/(MHk(i)**2-t_Mandel(tt,gen))
	 	 
				endif
								
			enddo
			
			Higgsbeitrag(uptype,gen) = scalarfactor(uptype,gen)
			
		enddo
	enddo
	
!----------------------- T3 & T4 Squark contribution s- & u-Channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
				tt = uptype + 2

				if (Denom_Mode==0) then

					scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 0.5d0*FL(1,i,tt,gen)
     &										*dconjg(FR(1,i,tt,gen))/MSf(i,tt,gen)**2

				else

					scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 0.25d0*FL(1,i,tt,gen)
     &										*dconjg(FR(1,i,tt,gen))/(MSf(i,tt,gen)**2 - s_Mandel(tt,gen))

					scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 0.25d0*FL(1,i,tt,gen)
     &										*dconjg(FR(1,i,tt,gen))/(MSf(i,tt,gen)**2 - u_Mandel(tt,gen))

				endif

			enddo

			Squarkbeitrag(uptype,gen) = scalarfactor(uptype,gen) - Higgsbeitrag(uptype,gen)

		enddo
	enddo
	
!----------------------- Total unsuppressed scalar contribution	
	do gen = 1,3
!		gen = 3
		do uptype = 1,2
			do nucleon = 1,2
				tt = uptype + 2
				g_4Fs(nucleon) = g_4Fs(nucleon) + f_T(nucleon,uptype,gen)*m_Nucleon(nucleon)
     &							 /Mf_low(tt,gen)*scalarfactor(uptype,gen)*Running_SI(tt,gen)
				g_4FsContribs(nucleon,1) = g_4FsContribs(nucleon,1) + f_T(nucleon,uptype,gen)*m_Nucleon(nucleon)
     &							 /Mf_low(tt,gen)*scalarfactor(uptype,gen)*Running_SI(tt,gen)
			enddo
			
c			write(*,*)'Squark Anteil SI:',abs(Squarkbeitrag(uptype,gen))/(abs(Squarkbeitrag(uptype,gen)) + abs(Higgsbeitrag(uptype,gen)))
						
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
	
!----------------------- T2 Z contribution t-channel
	do gen = 1,3
		do uptype = 1,2
			tt = uptype+2

			if(DD_udsCheck(tt,gen))then

			if (Denom_Mode ==0) then

				axialfactor(uptype,gen) = 0.5d0*AR(tt,gen)*(BL(tt,gen)-BR(tt,gen))/MZ2

			else

				axialfactor(uptype,gen) = 0.5d0*AR(tt,gen)*(BL(tt,gen)-BR(tt,gen))/(MZ2-t_Mandel(tt,gen))

			endif

			endif

			Zbeitrag(uptype,gen) = axialfactor(uptype,gen)

		enddo
	enddo
	
!----------------------- T3 & T4 Squark contribution s- & u-Channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
			tt = uptype+2

			if(DD_udsCheck(tt,gen)) then

			if (Denom_Mode ==0) then

				axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.25d0*(FR(1,i,tt,gen)*
     &				dconjg(FR(1,i,tt,gen)) + FL(1,i,tt,gen)*dconjg(FL(1,i,tt,gen)))/MSf2(i,tt,gen)

			else

				axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.125d0*(FR(1,i,tt,gen)*
     &				dconjg(FR(1,i,tt,gen)) + FL(1,i,tt,gen)*dconjg(FL(1,i,tt,gen)))/(MSf2(i,tt,gen)-s_Mandel(tt,gen))

				axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.125d0*(FR(1,i,tt,gen)*
     &				dconjg(FR(1,i,tt,gen)) + FL(1,i,tt,gen)*dconjg(FL(1,i,tt,gen)))/(MSf2(i,tt,gen)-u_Mandel(tt,gen))
	 
			endif

			endif
						
			enddo
			
			Squarkbeitrag(uptype,gen) = axialfactor(uptype,gen) - Zbeitrag(uptype,gen) 
			
		enddo
	enddo
		
!----------------------- Total unsuppressed axial contribution
	do gen = 1,3
		do uptype = 1,2
			do nucleon = 1,2
				g_4Fa(nucleon) = g_4Fa(nucleon) + Delta_q(nucleon,uptype,gen)
     &							 *axialfactor(uptype,gen)*Running_SD(tt,gen)
				g_4FaContribs(nucleon,1) = g_4FaContribs(nucleon,1) + Delta_q(nucleon,uptype,gen)
     &							 *axialfactor(uptype,gen)*Running_SD(tt,gen)
			enddo
			
c			write(*,*)'Squark Anteil SD:',abs(Squarkbeitrag(uptype,gen))/(abs(Squarkbeitrag(uptype,gen)) + abs(Zbeitrag(uptype,gen)))
c			write(*,*)'Zbeitrag SD:',abs(Zbeitrag(uptype,gen))
			
		enddo
	enddo
	
! ====================== SUPPRESSED CONTRIBUTIONS =================================================

	if(Supp_Mode==1.OR.Supp_Mode==2) then
	
	endif				! End of Supp-Mode if
	
! ====================== STRONGLY SUPPRESSED CONTRIBUTIONS ========================================

	if(Supp_Mode==2) then
	
	endif				! End of Supp-Mode if
	
!======================= END OF PROGRAM ===========================================================
	
	g_4Fsp = g_4Fs(1)
	g_4Fsn = g_4Fs(2)
	g_4Fap = g_4Fa(1)
	g_4Fan = g_4Fa(2)
	
	endif
	
	end


