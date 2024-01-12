! =========================================================================
	double precision function NeuQ2qx_M2Propagator()

	  implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	  double precision NeuQ2qx_M2TreeGAred, NeuQ2qx_M2PropagatorH, NeuQ2qx_M2PropagatorV

	  external myCget

	  NeuQ2qx_M2Propagator = 0d0

! according to the xtype the correct file with the amplitdues is chosen
! need to think about of implementing the charged Higgs with the charged W or the Higgses


	  if ((xtype.eq.1).or.(xtype.eq.2).or.(xtype.eq.3).or.(xtype.eq.4)) then
c	  print *, "xtype",xtype
	  NeuQ2qx_M2Propagator = 0d0
          NeuQ2qx_M2Propagator = NeuQ2qx_M2PropagatorV()

	  else if ((xtype.eq.5).or.(xtype.eq.6).or.(xtype.eq.7).or.(xtype.eq.8)) then
	  NeuQ2qx_M2Propagator = 0d0
	  NeuQ2qx_M2Propagator = NeuQ2qx_M2PropagatorH()

	  else
	  print *, "xtype",xtype
	  NeuQ2qx_M2Propagator = 0d0
	  endif

	end
