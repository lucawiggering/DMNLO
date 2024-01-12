      double precision function stsT2gg_M2DipolePlusFEIS(Q2IsTman)
      implicit none
#include "../common/stsT2xx_DipoleKinematics.h"
#include "../common/stsT2xx_Model.h"

      double precision tdip1,udip1,tdip,udip,colorDiptree1,colorDiptreeX
      double precision stsT2gg_M2TreeDipoleT1T3,stsT2gg_M2TreeDipoleT2T3
      double precision Jgg,Jqqbar,JhatStSt,VgFinite,VststFinite,Nf
      integer Q2IsTman

! number of effectively massless quark flavours
      Nf = 4d0
      
! Q2IsTman should take the value 1 if Q^2 corresponds to the Mandelstam variable t and 
! take 0 if Q^2 corresponds to the Mandelstam variable u
      if (Q2IsTman.eq.1) then
      tdip1 = q2
      udip1 = m2Sum - sdip1 - q2

      tdip = q2X
      udip = m2Sum - sdip - q2X

      colorDiptree1 = stsT2gg_M2TreeDipoleT1T3(sdip1,tdip1,udip1)
      colorDiptreeX = stsT2gg_M2TreeDipoleT1T3(sdip,tdip,udip)

      elseif (Q2IsTman.eq.0) then
      udip1 = q2
      tdip1 = m2Sum - sdip1  - q2

      udip = q2X
      tdip = m2Sum - sdip - q2X

      colorDiptree1 = stsT2gg_M2TreeDipoleT2T3(sdip1,tdip1,udip1)
      colorDiptreeX = stsT2gg_M2TreeDipoleT2T3(sdip,tdip,udip)

      else 
      
      write(*,*) "Q2IsTman has an invalid value. Q2IsTman = ",Q2IsTman

      endif
       
! The squared tree level matrix element is contained in the J-terms in the following definitions
      Jgg = dQ2*colordiptreeX/xdip**2*dlog(1+Adipx)*2/(Rdip*(1-xdip)) + yAprimePlus*APlus*dlog(1+APlus)*2*colordiptree1
     &           + dQ2*colordiptreeX/(xdip**2*(1-xdip))*(Rdip**2/6d0 - 2) + 11/(6*(1 - xdip))*colordiptree1

      Jqqbar = dQ2*colordiptreeX/(xdip**2*(1-xdip))*(1-Rdip**2/3d0) - 2/(3*(1-xdip))*colordiptree1

      JhatStSt = 2*(dQ2/(Rdip*xdip)*dlog(1+BdipX)/(1-xdip)*colordiptreeX+ yBprimePlus*Bplus*dlog(1+Bplus)*colordiptree1
     &               - dQ2/(Rdip*xdip)*dlog(1+AdipX)/(xdip*(1-xdip))*colordiptreeX - yAprimePlus*APlus*dlog(1+APlus)*colordiptree1
     &                  - dQ2/(Rdip*xdip)*BdipX/(1+BdipX)/(1-xdip)*(1+Rdip)*colordiptreeX/2d0 - yBprimePlus*BPlus**2/(1+BPlus)*colordiptree1)

! Definition of flavour functions
      VgFinite = CV*Jgg + TF*Nf*Jqqbar
      VststFinite = CF*JhatStSt

! FE-IS and IE-FS
      stsT2gg_M2DipolePlusFEIS = -gs**2/(4*Pi**2)*(VgFinite/CV + VststFinite/CF)

      end

      double precision function stsT2gg_M2DipolePlusIEIS()
      implicit none
#include "../common/stsT2xx_DipoleKinematics.h"
#include "../common/stsT2xx_Model.h"

      double precision tdip,udip,tdip1,udip1,colorDiptree1,colorDiptreeX
      double precision stsT2gg_M2TreeDipoleT1T2

! Mandelstam variables
      tdip = q2X
      udip = m2Sum - sdip - q2X

      tdip1 = q2
      udip1 = m2Sum - sdip1 - q2

! colour correlated leading order matrix elements squared
      colorDiptreeX = stsT2gg_M2TreeDipoleT1T2(sdip,tdip,udip)/CF
      colorDiptree1 = stsT2gg_M2TreeDipoleT1T2(sdip1,tdip1,udip1)/CF

      stsT2gg_M2DipolePlusIEIS = -gs**2*CF/(2*Pi**2*(1-xdip))*(d2*dlog(1+Cdip)-1)*(colorDiptreeX*dQ2/sqrtlambdaSdip - colorDiptree1/SqrtLambdaAB) 

      end