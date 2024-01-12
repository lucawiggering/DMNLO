! =========================================================================
	double precision function NeuQ2qx_M2Box()

	implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	double precision NeuQ2qx_M2BoxH, NeuQ2qx_M2BoxV, NeuQ2qx_M2BoxG

	external myCget

	 NeuQ2qx_M2Box = 0d0
	 
   	
	 if ((xtype.eq.1)) then
c		print *, "xtype",xtype
!	 	NeuQ2qx_M2Box = 0d0
	 	NeuQ2qx_M2Box = NeuQ2qx_M2BoxG()
	 	
	 else if ((xtype.eq.2).or.(xtype.eq.3).or.(xtype.eq.4)) then
c		print *, "xtype",xtype
	 	NeuQ2qx_M2Box = 0d0
	 	NeuQ2qx_M2Box = NeuQ2qx_M2BoxV()
c        print*, "NeuQ2qx_M2Box", NeuQ2qx_M2Box
c        print*, ""
	 	
	 else if ((xtype.eq.5).or.(xtype.eq.6).or.(xtype.eq.7).or.(xtype.eq.8)) then
		NeuQ2qx_M2Box = 0d0
	 	NeuQ2qx_M2Box = NeuQ2qx_M2BoxH()
c	 	print *, "xtype",xtype
c	 	print *, "NeuQ2qx_M2BoxH",NeuQ2qx_M2Box
	 	
	 else
	 	print *, "xtype",xtype
	        NeuQ2qx_M2Box = 0d0
	 endif

c        print*, "NeuQ2qx_M2Box", NeuQ2qx_M2Box
        


	end


