!==================================================================================================
!==================== This routine calculates the running from the high =========================== 
!====================   scale ~ 1 TeV to the low scale ~ 2-5 GeV	    ===========================
!====================			Author: Patrick Steppeler		        ===========================
!==================================================================================================

	subroutine DD_Running
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_Model.h"
#include "DD_GenCounterterms.h"

!----------------------- Declare variables --------------------------------------------------------
	integer nf											! Number of flavours
	integer gen											! Generation
	integer tt											! Fermion type
	double precision alpha_high, alpha_low, beta0		! alpha_s at different scales, first
														! coefficient of the QCD beta function
	double precision ASDR, MTDR
	double complex dMf_Temp								! To avoid that dMf(3,3) gets replaced when
														! calling MTDR
														
!======================= ACTUAL CALCULATIONS ======================================================
	nf = 5d0
	beta0 = 11d0 - 2d0/3d0*nf
	alpha_high = Alfas
	alpha_low = (4d0*pi)*ASDR(Q_low**2)

	do tt = 3,4
		do gen = 1,3
			Running_SI(tt,gen) = Mf_low(tt,gen)/Mf_high(tt,gen)
			Running_SD(tt,gen) = dexp(2*nf/(pi*beta0)*(alpha_high-alpha_low))
c			Running_SD(tt,gen) = 1d0
c			write(*,*)'Running_SD = ',Running_SD(tt,gen)
c			write(*,*)'Running_SI = ',Running_SI(tt,gen)
		enddo
	enddo
	
! In case of top quark we divide by m_t^DRbar(Q_high) and not by m_t^OS
	dMf_Temp = dMf(3,3)
	Running_SI(3,3) = Mf_low(3,3)/MTDR(MTos,QScale)
	dMf(3,3) = dMf_Temp
	
c	write(*,*)'Old top Running = ',Mf_low(3,3)/Mf_high(3,3)
c	write(*,*)'New top Running = ',Running_SI(3,3)
	

	
!======================= END OF PROGRAM ===========================================================

	end


