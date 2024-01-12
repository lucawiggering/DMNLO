!==================================================================================================
!======================= Prints all relevant DD parameters	=======================================
!=======================     Author: Patrick Steppeler		=======================================
!==================================================================================================

	subroutine DD_PrintModel(results)	
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_GenCouplings.h"
#include "DD_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

!----------------------- Declare variables --------------------------------------------------------
	integer nucleon, contrib
	double precision completeresults(2,2,10)	! Result array for neutralino nucleon cross section
												! First index:	Nucleon type 1=Proton, 2=Neutron
												! Second index:	Spin type 1=Scalar, 2=Axial
												! Third index: Contribution 1=Total, 2=Tree, 3=Vertex
												! 4=VertexC, 5=Prop, 6=PropC, 7=Box, 8=TreeHo
												! 9=EFTVertex, 10=EFTVertexC
	double precision results(8)				! Output array for neutralino nucleon cross section
		
!======================= ACTUAL CALCULATIONS ======================================================

	! write(*,*)'----------------------------------------'
	! write(*,*)'DD results @ UVdiv, IRdiv = ',UVdiv,IRdiv
	! write(*,*)''
	! write(*,*)'g_4Fsp: TOTAL / Tree', g_4Fsp, g_4FsContribs(1,1)
	! write(*,*)'g_4Fsp: Vertex / Vertexcounter', g_4FsContribs(1,2), g_4FsContribs(1,3)
	! write(*,*)'g_4Fsp: Prop / Propcounter', g_4FsContribs(1,4), g_4FsContribs(1,5)
	! write(*,*)'g_4Fsp: Box / TreeHO', g_4FsContribs(1,6), g_4FsContribs(1,7)
	! write(*,*)'g_4Fsp: EFTVertex / EFTVertexCounter', g_4FsContribs(1,8), g_4FsContribs(1,9)
	! write(*,*)''
	! write(*,*)'g_4Fsn: TOTAL / Tree', g_4Fsn, g_4FsContribs(2,1)
	! write(*,*)'g_4Fsn: Vertex / Vertexcounter', g_4FsContribs(2,2), g_4FsContribs(2,3)
	! write(*,*)'g_4Fsn: Prop / Propcounter', g_4FsContribs(2,4), g_4FsContribs(2,5)
	! write(*,*)'g_4Fsn: Box / TreeHO', g_4FsContribs(2,6), g_4FsContribs(2,7)
	! write(*,*)'g_4Fsn: EFTVertex / EFTVertexCounter', g_4FsContribs(2,8), g_4FsContribs(2,9)
	! write(*,*)''
	! write(*,*)'g_4Fap: TOTAL / Tree', g_4Fap, g_4FaContribs(1,1)
	! write(*,*)'g_4Fap: Vertex / Vertexcounter', g_4FaContribs(1,2), g_4FaContribs(1,3)
	! write(*,*)'g_4Fap: Prop / Propcounter', g_4FaContribs(1,4), g_4FaContribs(1,5)
	! write(*,*)'g_4Fap: Box / TreeHO', g_4FaContribs(1,6), g_4FaContribs(1,7)
	! write(*,*)'g_4Fap: EFTVertex / EFTVertexCounter', g_4FaContribs(1,8), g_4FaContribs(1,9)
	! write(*,*)''
	! write(*,*)'g_4Fan: TOTAL / Tree', g_4Fan, g_4FaContribs(2,1)
	! write(*,*)'g_4Fan: Vertex / Vertexcounter', g_4FaContribs(2,2), g_4FaContribs(2,3)
	! write(*,*)'g_4Fan: Prop / Propcounter', g_4FaContribs(2,4), g_4FaContribs(2,5)
	! write(*,*)'g_4Fan: Box / TreeHO', g_4FaContribs(2,6), g_4FaContribs(2,7)
	! write(*,*)'g_4Fan: EFTVertex / EFTVertexCounter', g_4FaContribs(2,8), g_4FaContribs(2,9)
	! write(*,*)''
	! write(*,*)'----------------------------------------'
	
	do nucleon = 1,2
		completeresults(nucleon,1,1) = sigmaXs(nucleon)
		completeresults(nucleon,2,1) = sigmaXa(nucleon)
		do contrib=2,10
			completeresults(nucleon,1,contrib) = sigmaXsContribs(nucleon,(contrib-1))
			completeresults(nucleon,2,contrib) = sigmaXaContribs(nucleon,(contrib-1))
		enddo
	enddo
	
	results(1) = sigmaXs(1)
	results(2) = sigmaXsContribs(1,1)
	results(3) = sigmaXa(1)
	results(4) = sigmaXaContribs(1,1)
	results(5) = sigmaXs(2)
	results(6) = sigmaXsContribs(2,1)
	results(7) = sigmaXa(2)
	results(8) = sigmaXaContribs(2,1)	

!======================= END OF PROGRAM ===========================================================

	end


