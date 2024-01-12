      double precision function NeuQ2qx_M2Tree()

      implicit none
	
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

	double precision NeuQ2qx_M2TreeGA, NeuQ2qx_M2TreeH, NeuQ2qx_M2TreeV, NeuQ2qx_M2TreeGAred, NeuQ2qx_M2TreeHred, NeuQ2qx_M2TreeHgen, NeuQ2qx_M2TreeVgen

	NeuQ2qx_M2Tree	= 0d0
! according to the xtype the correct file with the amplitdues is chosen
! need to think about of implementing the charged Higgs with the charged W or the Higgses
   	
   	
	 if ((xtype.eq.1).or.(xtype.eq.2)) then
	 	NeuQ2qx_M2Tree = 0d0
	 	NeuQ2qx_M2Tree = NeuQ2qx_M2TreeVgen()
	 	
	 else if ((xtype.eq.3).or.(xtype.eq.4)) then
	 	NeuQ2qx_M2Tree = 0d0
	 	NeuQ2qx_M2Tree = NeuQ2qx_M2TreeVgen()
	 	
	 else if ((xtype.eq.5).or.(xtype.eq.6).or.(xtype.eq.7).or.(xtype.eq.8)) then
		NeuQ2qx_M2Tree = 0d0
	 	NeuQ2qx_M2Tree = NeuQ2qx_M2TreeHgen()
	 else
	    NeuQ2qx_M2Tree = 0d0
	 endif
	
	
	
	end

      double precision function NeuQ2qx_M2TreeDipole(sdip,tdip,udip)

      implicit none
	
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

	double precision NeuQ2qx_M2TreeGA_dipole,NeuQ2qx_M2TreeWZ_dipole, NeuQ2qx_M2TreeH_Dipole
      double precision sdip,tdip,udip
	NeuQ2qx_M2TreeDipole	= 0d0
      
      if ((xtype.eq.1).or.(xtype.eq.2)) then ! gluon or photon
            NeuQ2qx_M2TreeDipole =   NeuQ2qx_M2TreeGA_dipole(sdip,tdip,udip)

      else if ((xtype.eq.3).or.(xtype.eq.4)) then ! W or Z
            NeuQ2qx_M2TreeDipole =  NeuQ2qx_M2TreeWZ_dipole(sdip,tdip,udip)	 

      else if ((xtype.eq.5).or.(xtype.eq.6).or.(xtype.eq.7).or.(xtype.eq.8)) then ! Higgses
            NeuQ2qx_M2TreeDipole = NeuQ2qx_M2TreeH_Dipole(sdip,tdip,udip)
      else
            NeuQ2qx_M2TreeDipole = 0d0
      endif	
	
	end

