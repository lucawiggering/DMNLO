!==================================================================================================
!==================== Small function that checks if quark is u,d or s =============================
!====================			Author: Patrick Steppeler		      =============================
!==================================================================================================

	function DD_udsCheck(tt,gen)	
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Declare variables --------------------------------------------------------
	integer gen				! Generation
	integer tt				! Fermion type
	logical DD_udsCheck		
	
!======================= ACTUAL CALCULATIONS ======================================================
	DD_udsCheck=.false.
!	DD_udsCheck=.true.		! Deactivate uds-Check
	
	if((gen==1).and.((tt==3).or.(tt==4))) then
		DD_udsCheck=.true.
	else if((gen==2).and.(tt==4)) then
		DD_udsCheck=.true.
	endif
	
!======================= END OF PROGRAM ===========================================================

	end


