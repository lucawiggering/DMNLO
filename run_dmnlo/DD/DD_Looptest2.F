!==================================================================================================
!============= Tests the loopfunction in the case of vanishing quark masses =======================
!=======================       Author: Patrick Steppeler		  =================================
!==================================================================================================

	subroutine looptest2(Danger, result)	
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_GenCouplings.h"
#include "DD_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"
	
!----------------------- Declare variables --------------------------------------------------------
	double precision Danger
	double precision result(20)
	complex *16 C03
	double precision MGluino_Marcel, t_Marcel
	
!----------------------- Declare external subroutines ---------------------------------------------
	external DMNLO_ModelPara, DMNLO_ModelIni, DMNLO_ModelDigest, Init_RenScheme
	
!======================= ACTUAL CALCULATIONS ======================================================

	MGluino_Marcel = 1516.0d0
	t_Marcel = -123456.0d0
	
	call DMNLO_ModelPara()
	call DMNLO_ModelIni()
	call Init_RenScheme(1,MGluino_Marcel)
!	call DMNLO_ModelDigest
	
!----------------------- Renormalization Scales for Loop integrals --------------------------------
	UVdiv = 0d1
	IRdiv = 0d1
 	EpsPole = 0d0
	
!----------------------- Call Loop Functions ------------------------------------------------------
!	call myCget(Danger**2, -20.0d0, Danger**2, 0.0d0, Danger**2, Danger**2)
!	call myCget(0.0d0, -20.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0)

	

	call myCget(0d0, t_Marcel, MGluino_Marcel**2, 0d0, 0d0, MGluino_Marcel**2)
	call myBget(t_Marcel, 0d0, MGluino_Marcel**2)
	result(1) = dreal(xc0)
	
	write(*,*)'-------- LOOPTEST ---------'
	write(*,*)'t_Mandel = ', t_Marcel
	write(*,*)'M_Gluino = ', MGluino_Marcel
	write(*,*)'C0 = ', xc0
	write(*,*)'B0 = ', xb0
	
!	write(*,*)'Denner C0 = ', C03(Mf(3,3)**2, Danger, Mf(3,3)**2, MGl, MSf(1,3,3), MSf(1,3,3))
	
!======================= END OF PROGRAM ===========================================================

	end

