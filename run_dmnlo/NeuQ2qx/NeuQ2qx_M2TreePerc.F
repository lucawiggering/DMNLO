************************************************

	double precision function NeuQ2qx_M2TreePerc(ichannel)

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"


	integer i,j,ichannel
	complex*16 DenSmZ, DenSmH(4)
	double precision DenSmZ1, DenSmH1(4)

	NeuQ2qx_M2treePerc = 0d0
	
! ===== One contribution =====
	if (ichannel.eq.1) then

		NeuQ2qx_M2TreePerc = NeuQ2qx_M2TreePerc + 0d0		

! ===== Another contribution =====
	elseif (ichannel.eq.2) then

		NeuQ2qx_M2TreePerc = NeuQ2qx_M2TreePerc + 0d0		

! ===== Yet another contribution =====
	elseif (ichannel.eq.3) then

		NeuQ2qx_M2TreePerc = NeuQ2qx_M2TreePerc + 0d0		

	endif

	end
