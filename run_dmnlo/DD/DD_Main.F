!==================================================================================================
!======================= Main file for direct detection in DM@NLO =================================
!=======================       Author: Patrick Steppeler		  =================================
!==================================================================================================

	subroutine directdetection(p_cm, costh, inLoop_Mode, inDenom_Mode, inBox_Mode, inNucl_Mode, Qout, results,iflags,slhafilename)	
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_GenCouplings.h"
#include "DD_GenCounterterms.h"
#include "DD_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"
	
!----------------------- Declare variables --------------------------------------------------------
	double precision p_cm			! Momentum in CMS frame in GeV
	double precision costh			! Scattering angle
	
	double precision Qout			! Renormalization scale
	double precision results(8)		! Result array for neutralino nucleon cross section
	integer inSupp_Mode				! 0: Only unsuppressed, 1: Up to suppressed, 2: Up to strongly suppressed 
	integer inLoop_Mode				! 0: Only tree level, 1: Add loops
	integer inDenom_Mode			! 0: Only masses in denominators of propagators, 1: Complete kinematics
	integer inBox_Mode				! 0: Old B1/B3 amplitudes, 1: New B1/B3 amplitudes which are correct in the limit v -> 0
									! When performing IR/Checks for v != 0 one should use Box_mode = 0, otherwise Box_Mode = 1
	integer inNucl_Mode				! Choose nuclear input (f_T's & Deltas)
									! 1: micrOMEGAs 2.4.1 setup, 2: micrOMEGAs 4.0.3 setup, 3: DarkSUSY setup, 4 DM@NLO setup
	integer counter
	integer i,j,k,l,tt,gen
	integer iflags(3),MOswitch

	character slhafilename*200
	
!----------------------- Declare external subroutines ---------------------------------------------
	external DMNLO_ModelPara, DMNLO_ModelIni, DMNLO_ModelDigest, Init_RenScheme
	external DD_Init, DD_SetKinematics, DD_SetCouplings, DD_Tree, DD_PrintModel
	external DD_Vertex, DD_VertexCounter, DD_Prop, DD_PropCounter, DD_Box, DD_SigmaNucleon
	
!======================= ACTUAL CALCULATIONS ======================================================
	
!----------------------- Loop over div
!	do counter = 0,8

!----------------------- Renormalization Scales for Loop integrals --------------------------------
!	UVdiv = 10.0d0**counter
!	IRdiv = 10.0d0**counter
	UVdiv = 0d3
	IRdiv = 0d3
	EpsPole = 0d0		!1: UV poles, 0: Full result, -1 IR single poles, -2 IR double poles
	
	QScale = Qout
	Qscalealphas = Qscale

	ischeme = iflags(1)
	MOswitch = iflags(2)
	choosesol = iflags(3)
	
	call DMNLO_ModelPara(MOswitch,slhafilename)
	call DMNLO_ModelIni(slhafilename)

	call Init_RenScheme(ischeme)	!0: Pure DR Bar scheme, 1: Hybrid on-shell/DR Bar scheme
!	call DMNLO_ModelDigest
	if(renfail.eq.0) stop
		
	inSupp_Mode = 0 !  0 -> only unsupressed, 1 -> up to supressed (not implemented yet), 2 -> up to strongly supressed (not implemented yet)
	Supp_Mode 	= inSupp_Mode
	Loop_Mode 	= inLoop_Mode
	Denom_Mode 	= inDenom_Mode
	Box_Mode 	= inBox_Mode
	Nucl_Mode 	= inNucl_Mode
	
	call DD_Init
	call DD_SetKinematics(p_cm, costh)
	call DD_SetCouplings
	call DD_Running
	call DD_Tree
	
	if (Loop_Mode==1) then

		call DD_SetCounterterms
		call DD_SetCouplingsHO
		call DD_Vertex
		call DD_VertexCounter
		call DD_Prop
		call DD_PropCounter
		call DD_Box
		call DD_TreeHO
		call DD_EFTVertex
		call DD_EFTVertexCounter
		
	endif
	
!	call DD_Looptest
	call DD_SigmaNucleon
	call DD_PrintModel(results)
	  
!======================= END OF PROGRAM ===========================================================

! End of div loop
!	enddo
	
	end

!----------------------- Include external subroutines ---------------------------------------------
#include "DD_Init.F"
#include "DD_EFTVertex.F"
#include "DD_EFTVertexCounter.F"
#include "DD_udsCheck.F"
#include "DD_Kinematics.F"
#include "DD_SetCouplings.F"
#include "DD_SetCouplingsHO.F"
#include "DD_SetCounterterms.F"
#include "DD_Tree.F"
#include "DD_TreeHO.F"
#include "DD_Vertex.F"
#include "DD_VertexCounter.F"
#include "DD_Prop.F"
#include "DD_PropCounter.F"
#include "DD_Box.F"
#include "DD_PrintModel.F"
#include "DD_SigmaNucleon.F"
#include "DD_Running.F"
#include "DD_Looptest.F"

