!==================================================================================================
!======================= Calculates the neutralino nucleon cross section	=======================
!=======================			Author: Patrick Steppeler				=======================
!==================================================================================================

	subroutine DD_SigmaNucleon	
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_Model.h"

!----------------------- Declare variables --------------------------------------------------------
	integer nucleon, contrib
	double precision gev2cm2	! Conversion factor to the output format in cm^2
	
!======================= ACTUAL CALCULATIONS ======================================================
	gev2cm2=(197.327053d-16)**2
	
! Check imaginary part of the four fermion coupling. A too large imaginary part signalises trouble
! in the loop functions.
	if(abs(1000d0*AIMAG(g_4Fs(1))).gt.abs(REAL(g_4Fs(1))).or.abs(1000d0*AIMAG(g_4Fs(2))).gt.abs(REAL(g_4Fs(2))).or.
     & abs(1000d0*AIMAG(g_4Fa(1))).gt.abs(REAL(g_4Fa(1))).or.abs(1000d0*AIMAG(g_4Fa(2))).gt.abs(REAL(g_4Fa(2)))) then
		! write(*,*)'WARNING! IMAGINARY PART OF FOUR FERMION COUPLINGS TO BIG!'
		! write(*,*)'gFs1 = ',g_4Fs(1)
		! write(*,*)'gFs2 = ',g_4Fs(2)
		! write(*,*)'gFa1 = ',g_4Fa(1)
		! write(*,*)'gFa2 = ',g_4Fa(2)
	endif
	
	do nucleon=1,2
		sigmaXs(nucleon) = gev2cm2*mu_Nucleon(nucleon)**2/pi * dconjg(g_4Fs(nucleon))*g_4Fs(nucleon)
		sigmaXa(nucleon) = gev2cm2*3.0d0*mu_Nucleon(nucleon)**2/pi * dconjg(g_4Fa(nucleon))*g_4Fa(nucleon)
		do contrib=1,9
			sigmaXsContribs(nucleon,contrib) = gev2cm2*mu_Nucleon(nucleon)**2/pi * 
     &			g_4FsContribs(nucleon, contrib)*dconjg(g_4FsContribs(nucleon, contrib))
			sigmaXaContribs(nucleon,contrib) = gev2cm2*3.0d0*mu_Nucleon(nucleon)**2/pi * 
     &			g_4FaContribs(nucleon, contrib)*dconjg(g_4FaContribs(nucleon, contrib))
		enddo
	enddo
	
	sigmaXPs = sigmaXs(1)
	sigmaXPa = sigmaXa(1)
	sigmaXNs = sigmaXs(2)
	sigmaXNa = sigmaXa(2)
	

!======================= END OF PROGRAM ===========================================================

	end


