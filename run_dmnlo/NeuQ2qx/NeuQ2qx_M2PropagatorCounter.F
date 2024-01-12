! =========================================================================
	double precision function NeuQ2qx_M2PropagatorCounter()

	implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	double precision NeuQ2qx_M2TreeGAred, NeuQ2qx_M2TreeV, NeuQ2qx_M2PropagatorCounterH, NeuQ2qx_M2PropagatorCounterV

	external myCget

	 NeuQ2qx_M2PropagatorCounter = 0d0

! according to the xtype the correct file with the amplitdues is chosen
! need to think about of implementing the charged Higgs with the charged W or the Higgses


         if ((xtype.eq.1).or.(xtype.eq.2).or.(xtype.eq.3).or.(xtype.eq.4)) then
	 	NeuQ2qx_M2PropagatorCounter = 0d0
	 	NeuQ2qx_M2PropagatorCounter = NeuQ2qx_M2PropagatorCounterV()
!		print *, "xtype",xtype
!	 	print *, "NeuQ2qx_M2PropagatorCounterV",NeuQ2qx_M2PropagatorCounter


	 else if ((xtype.eq.5).or.(xtype.eq.6).or.(xtype.eq.7).or.(xtype.eq.8)) then
		NeuQ2qx_M2PropagatorCounter = 0d0
	 	NeuQ2qx_M2PropagatorCounter = NeuQ2qx_M2PropagatorCounterH()
!	 	print *, "xtype",xtype
!	 	print *, "NeuQ2qx_M2PropagatorCounterH",NeuQ2qx_M2PropagatorCounter

	 else
	 	print *, "xtype",xtype
	        NeuQ2qx_M2PropagatorCounter = 0d0
	 endif

c        print*, "NeuQ2qx_M2PropagatorCounter", NeuQ2qx_M2PropagatorCounter


	  end


