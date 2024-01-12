      double precision function NeuQ2qx_M2DipolePlusG_Gluon()
      implicit none 

#include "NeuQ2qx_DipoleKinematics.h"
#include "NeuQ2qx_Model.h"

      double precision Nf,diptreeX,diptree1,NeuQ2qx_M2TreeDipole
      double precision tdip,udip,tdip1,udip1
      double precision Jgg,Jqqbar,JhatStSt,VgFinite,VststFinite

! number of effectively massless quark flavours
      Nf = 4d0

! Mandelstam variables
      tdip = q2X
      udip = m1**2 + m2**2 + m3**2 + m4**2 - sdip - tdip

      tdip1 = q2
      udip1 = m1**2 + m2**2 + m3**2 + m4**2 - sdip1 - tdip1

! leading order matrix elements squared
      diptreeX = NeuQ2qx_M2TreeDipole(sdip,tdip,udip)
      diptree1 = NeuQ2qx_M2TreeDipole(sdip1,tdip1,udip1)

! The squared tree level matrix element is contained in the J-terms in the following definitions
      Jgg = dQ2*diptreeX/xdip**2*dlog(1+Adipx)*2/(Rdip*(1-xdip)) + yAprimePlus*APlus*dlog(1+APlus)*2*diptree1
     &           + dQ2*diptreeX/(xdip**2*(1-xdip))*(Rdip**2/6d0 - 2) + 11/(6*(1 - xdip))*diptree1

      Jqqbar = dQ2*diptreeX/(xdip**2*(1-xdip))*(1-Rdip**2/3d0) - 2/(3*(1-xdip))*diptree1

      JhatStSt = 2*(dQ2/(Rdip*xdip)*dlog(1+BdipX)/(1-xdip)*diptreeX + yBprimePlus*Bplus*dlog(1+Bplus)*diptree1
     &               - dQ2/(Rdip*xdip)*dlog(1+AdipX)/(xdip*(1-xdip))*diptreeX - yAprimePlus*APlus*dlog(1+APlus)*diptree1
     &                  - dQ2/(Rdip*xdip)*BdipX/(1+BdipX)/(1-xdip)*(1+Rdip)*diptreeX/2d0 - yBprimePlus*BPlus**2/(1+BPlus)*diptree1)

! Definition of flavour functions
      VgFinite = CV*Jgg + TF*Nf*Jqqbar
      VststFinite = CF*JhatStSt

! FE-IS and IE-FS
      NeuQ2qx_M2DipolePlusG_Gluon = gs**2/(16*Pi**2)*(VgFinite + CV/CF*VststFinite)

      end


      double precision function NeuQ2qx_M2DipolePlusG_Quark()
      implicit none 
#include "NeuQ2qx_DipoleKinematics.h"
#include "NeuQ2qx_Model.h"

      double precision Nf,diptreeX,diptree1, NeuQ2qx_M2TreeDipole
      double precision tdip,udip,tdip1,udip1,JgQ,Jstst,VqFinite,VststFinite

! Mandelstam variables
      udip = q2X
      tdip = m1**2 + m2**2 + m3**2 + m4**2 - sdip - udip

      udip1 = q2
      tdip1 = m1**2 + m2**2 + m3**2 + m4**2 - sdip1 - udip1

! leading order matrix elements squared
      diptreeX = NeuQ2qx_M2TreeDipole(sdip,tdip,udip)
      diptree1 = NeuQ2qx_M2TreeDipole(sdip1,tdip1,udip1)
      
! The squared tree level matrix element is contained in the J-terms in the following definitions
      JgQ = dQ2*2/(xdip**2*(1-xdip))*((xdip-1)**2/(4*(xdip*(etaJX-1)+1)**2) -1 + dlog(1+AdipX)/(Rdip*vdipX))*diptreeX
     &          - 2/(1-xdip)*(-1 + dlog(1+Adip1)/vdip1)*diptree1

      Jstst = dQ2*diptreeX/(1-xdip)*2/(vDipX*Rdip*xdip**2)*(xdip*dlog(1+Bdipx) - dlog(1+Adipx) - 4*etaAX*xdip**2*((etaJX-1)*xdip+1)/(1-vdipX**2*Rdip**2)*vdipX*Rdip)
     &          - diptree1/(1-xdip)*2/vDip1*(dlog(1+Bdip1) - dlog(1+Adip1) - 4*etaA1*etaJ1/(1-vDip1**2)*vDip1)   

! Definition of flavour functions
      VqFinite = CF*JgQ
      VststFinite = CF*Jstst
      
! FE-IS and IE-FS
      NeuQ2qx_M2DipolePlusG_Quark = gs**2/(16*Pi**2)*(2-CV/CF)*(VqFinite+ VststFinite)

      end



      double precision function NeuQ2qx_M2DipolePlusEW()
      implicit none 
#include "NeuQ2qx_DipoleKinematics.h"
#include "NeuQ2qx_Model.h"

      double precision Nf,diptreeX,diptree1, NeuQ2qx_M2TreeDipole
      double precision tdip,udip,tdip1,udip1,JgQ,Jstst,VqFinite,VststFinite

! Mandelstam variables
      udip = q2X
      tdip = m1**2 + m2**2 + m3**2 + m4**2 - sdip - udip

      udip1 = q2
      tdip1 = m1**2 + m2**2 + m3**2 + m4**2 - sdip1 - udip1

! leading order matrix elements squared
      diptree1 = NeuQ2qx_M2TreeDipole(sdip1,tdip1,udip1)
      diptreeX = NeuQ2qx_M2TreeDipole(sdip,tdip,udip)

!=============================================================================================================================================
! The squared tree level matrix element is contained in the J-terms in the following definitions
      JgQ = dQ2*2/(xdip**2*(1-xdip))*((xdip-1)**2/(4*(xdip*(etaJX-1)+1)**2) -1 + dlog(1+AdipX)/(Rdip*vdipX))*diptreeX
     &          - 2/(1-xdip)*(-1 + dlog(1+Adip1)/vdip1)*diptree1

      Jstst = dQ2*diptreeX/(1-xdip)*2/(vDipX*Rdip*xdip**2)*(xdip*dlog(1+Bdipx) - dlog(1+Adipx) - 4*etaAX*xdip**2*((etaJX-1)*xdip+1)/(1-vdipX**2*Rdip**2)*vdipX*Rdip)
     &          - diptree1/(1-xdip)*2/vDip1*(dlog(1+Bdip1) - dlog(1+Adip1) - 4*etaA1*etaJ1/(1-vDip1**2)*vDip1)   

! Definition of flavour functions
      VqFinite = CF*JgQ
      VststFinite = CF*Jstst
      
! FE-IS and IE-FS
      NeuQ2qx_M2DipolePlusEW = gs**2/(8*Pi**2)*(VqFinite+ VststFinite)

      end

