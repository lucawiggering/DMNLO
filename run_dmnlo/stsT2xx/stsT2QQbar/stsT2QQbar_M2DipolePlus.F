      double precision function stsT2QQbar_M2DipolePlusFEIS(Q2IsTman)
      implicit none
#include "../common/stsT2xx_DipoleKinematics.h"
#include "../common/stsT2xx_Model.h"

      double precision tdip1,udip1,tdip,udip,colorDiptree1,colorDiptreeX
      double precision stsT2QQbar_M2TreeDipole
      double precision JgQ,JhatStSt,VqFinite,VststFinite
      integer Q2IsTman
      
! Q2IsTman should take the value 1 if Q^2 corresponds to the Mandelstam variable t and 
! take 0 if Q^2 corresponds to the Mandelstam variable u
      if (Q2IsTman.eq.1) then
      tdip1 = q2
      udip1 = m2Sum - sdip1 - q2

      tdip = q2X
      udip = m2Sum - sdip - q2X

      colorDiptree1 = (2-Nc**2)/(2*Nc)*stsT2QQbar_M2TreeDipole(sdip1,tdip1,udip1)/CF ! T1 T3
      colorDiptreeX = (2-Nc**2)/(2*Nc)*stsT2QQbar_M2TreeDipole(sdip, tdip,udip)/CF

      elseif (Q2IsTman.eq.0) then
      udip1 = q2
      tdip1 = m2Sum - sdip1  - q2

      udip = q2X
      tdip = m2Sum - sdip - q2X

      colorDiptree1 = stsT2QQbar_M2TreeDipole(sdip1,tdip1,udip1)/(-Nc*CF) !T1 T4
      colorDiptreeX = stsT2QQbar_M2TreeDipole(sdip, tdip,udip)/(-Nc*CF)

      else 
      
      write(*,*) "Q2IsTman has an invalid value. Q2IsTman = ",Q2IsTman

      endif
       
! The squared tree level matrix element is contained in the J-terms in the following definitions
      JgQ = (-6*colorDiptreeX*dQ2)/(4*xdip**2*(1-xdip)) + (6*colorDiptree1)/(4*(1-xdip))
     &    + (2*colorDiptreeX*dlog(1+AdipX)*dQ2)/(xdip**2*(1-xdip)*Rdip) + 2*yAprimePlus*APlus*dlog(1+APlus)*colorDiptree1

      JhatStSt = 2*(dQ2/(Rdip*xdip)*dlog(1+BdipX)/(1-xdip)*colordiptreeX+ yBprimePlus*Bplus*dlog(1+Bplus)*colordiptree1
     &               - dQ2/(Rdip*xdip)*dlog(1+AdipX)/(xdip*(1-xdip))*colordiptreeX - yAprimePlus*APlus*dlog(1+APlus)*colordiptree1
     &                  - dQ2/(Rdip*xdip)*BdipX/(1+BdipX)/(1-xdip)*(1+Rdip)*colordiptreeX/2d0 - yBprimePlus*BPlus**2/(1+BPlus)*colordiptree1)

! Definition of flavour functions
      VqFinite = CF*JgQ
      VststFinite = CF*JhatStSt

! FE-IS and IE-FS
      stsT2QQbar_M2DipolePlusFEIS = -gs**2/(4*Pi**2)*(VqFinite + VststFinite)

      end

      double precision function stsT2QQbar_M2DipolePlusIEIS()
      implicit none
#include "../common/stsT2xx_DipoleKinematics.h"
#include "../common/stsT2xx_Model.h"

      double precision tdip,udip,tdip1,udip1,colorDiptree1,colorDiptreeX
      double precision stsT2QQbar_M2TreeDipole

! Mandelstam variables
      tdip = q2X
      udip = m2Sum - sdip - q2X

      tdip1 = q2
      udip1 = m2Sum - sdip1 - q2

! colour correlated leading order matrix elements squared
      colorDiptree1 = stsT2QQbar_M2TreeDipole(sdip1,tdip1,udip1)/(2*Nc*CF) !T1 T2
      colorDiptreeX = stsT2QQbar_M2TreeDipole(sdip, tdip,udip)/(2*Nc*CF)

      stsT2QQbar_M2DipolePlusIEIS = -gs**2*CF/(2*Pi**2*(1-xdip))*(d2*dlog(1+Cdip)-1)*(colorDiptreeX*dQ2/sqrtlambdaSdip - colorDiptree1/SqrtLambdaAB) 

      end