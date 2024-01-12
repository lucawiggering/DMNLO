! =========================================================================
	double precision function NeuQ2qx_M2Vertex()

	implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	double precision NeuQ2qx_M2VertexH, NeuQ2qx_M2VertexV, NeuQ2qx_M2VertexG

	external myCget

	 NeuQ2qx_M2Vertex = 0d0
	 
! according to the xtype the correct file with the amplitdues is chosen
! need to think about of implementing the charged Higgs with the charged W or the Higgses
   	
   	
	 if (xtype.eq.1) then
c		print *, "xtype",xtype
	 	NeuQ2qx_M2Vertex = 0d0
	 	NeuQ2qx_M2Vertex = NeuQ2qx_M2VertexG()
	 	
	 else if ((xtype.eq.2).or.(xtype.eq.3).or.(xtype.eq.4)) then
c		print *, "xtype",xtype
	 	NeuQ2qx_M2Vertex = 0d0
	 	NeuQ2qx_M2Vertex = NeuQ2qx_M2VertexV()
!                print*, "NeuQ2qx_M2Vertex", NeuQ2qx_M2Vertex
	 	
	 else if ((xtype.eq.5).or.(xtype.eq.6).or.(xtype.eq.7).or.(xtype.eq.8)) then
		NeuQ2qx_M2Vertex = 0d0
	 	NeuQ2qx_M2Vertex = NeuQ2qx_M2VertexH()
c	 	print *, "xtype",xtype
c	 	print *, "NeuQ2qx_M2VertexH",NeuQ2qx_M2Vertex
	 	
	 else
	 	print *, "xtype",xtype
	        NeuQ2qx_M2Vertex = 0d0
	 endif

c        print*, "NeuQ2qx_M2Vertex", NeuQ2qx_M2Vertex
        


	end


