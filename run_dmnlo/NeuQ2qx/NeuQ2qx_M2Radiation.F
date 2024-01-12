      double precision function NeuQ2qx_M2Radiation()

      implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

      double precision NeuQ2qx_M2RadiationHSS, NeuQ2qx_M2RadiationHTT, NeuQ2qx_M2RadiationHUU, NeuQ2qx_M2RadiationHST, NeuQ2qx_M2RadiationHSU, NeuQ2qx_M2RadiationHTU
      double precision NeuQ2qx_M2RadiationVSS, NeuQ2qx_M2RadiationVTT, NeuQ2qx_M2RadiationVUU, NeuQ2qx_M2RadiationVST, NeuQ2qx_M2RadiationVSU, NeuQ2qx_M2RadiationVTU	


      NeuQ2qx_M2Radiation = 0d0


      if (xtype.eq.1) then
      NeuQ2qx_M2Radiation = 0d0

      else if ((xtype.eq.2).or.(xtype.eq.3).or.(xtype.eq.4)) then

            NeuQ2qx_M2Radiation = NeuQ2qx_M2RadiationVSS() + NeuQ2qx_M2RadiationVST() + NeuQ2qx_M2RadiationVSU() + NeuQ2qx_M2RadiationVTT() + NeuQ2qx_M2RadiationVTU() + NeuQ2qx_M2RadiationVUU()

      else if ((xtype.eq.5).or.(xtype.eq.6).or.(xtype.eq.7).or.(xtype.eq.8)) then

            NeuQ2qx_M2Radiation = NeuQ2qx_M2RadiationHSS() + NeuQ2qx_M2RadiationHTT() + NeuQ2qx_M2RadiationHUU() + NeuQ2qx_M2RadiationHST() + NeuQ2qx_M2RadiationHTU() + NeuQ2qx_M2RadiationHSU()  

      else
      print *, "xtype",xtype
      NeuQ2qx_M2Radiation = 0d0
      endif



      end
