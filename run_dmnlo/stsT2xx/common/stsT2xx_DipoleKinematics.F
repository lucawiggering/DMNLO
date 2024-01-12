      subroutine stsT2xx_SetDipoleKinematicsFEIS(xi,x0,sman,ma2,mmb2,mk2,mj2,intfac)
      implicit none 
! For the evaluation of the "plus"-distribution one has to distinguish between three different dependencies $f(x,Q^2)$. 
! For this purpose, we introduce two different variables for $Q^2$. These are q2 which is generated from the integration limits
! of $Q^2$ evaluated at x=1 $Q^2_\pm(1)$ and q2X which is generated from the $x$-dependent integration limits $Q^2_\pm(x)$ instead. 
! Variables with the suffix "X" refer to f(x,q2X), variables with the suffix "1" refer to f(1,q2) and those with the suffix "plus" correspond to f(x,q2).
#include "stsT2xx_Model.h"
#include "stsT2xx_DipoleKinematics.h"
      double precision xi(2),ma2,mmb2,mk2,mj2,intfac,kappa,alpha,beta,prefac,x0
      double precision q2max, q2min,q2maxAt1, q2minAt1,sman, zP, zM,rhoPlus, calAPlus,calBPlus

! Sum of squared masses for calculation of tdip/udip
      m2Sum = ma2 + mmb2 + mk2 + mj2

! set squared CM entergy for x=1
      sdip1 = sman

! generate value for x      
      xdip = (1d0-x0)*xi(1) + x0   

! integration limits for Q2
      prefac = 1d0/(2*(xdip*sman + (1d0-xdip)*(mmb2 - xdip*ma2)))

      alpha =  xdip**2*(ma2**2 + 2*ma2*(mmb2+mk2)-(mmb2-sman)**2)+ 2*mmb2*(ma2+mj2)
     &                 -xdip*(ma2**2 + ma2*(4*mmb2+mj2+mk2-sman) - (mmb2-sman)*(mmb2 - mj2 - mk2))

      beta =  xdip*kappa(ma2,mmb2,sman)*dsqrt((ma2-mmb2)**2*(1d0-xdip)**2 + (1d0-xdip)*(2*ma2*(mj2 + mk2*(2*xdip-1d0)-sman*xdip)
     &               -2*mmb2*(mj2+mk2-sman*xdip))+ kappa(xdip*sman,mk2,mj2)**2)

      q2max = prefac*(alpha + beta)
      q2min = prefac*(alpha - beta)

! q2minAt1 and q2maxAt1 correspond to the integration limits of Q2 evaluated at P2=mj2 which is equivalent to x=1
      q2minAt1 = mmb2 + mk2 - ((-ma2 + mmb2 + sman)*(mk2 - mj2 + sman))/(2*sman) - (kappa(sman,ma2,mmb2)*kappa(sman,mk2,mj2))/(2*sman)
      q2maxAt1 = mmb2 + mk2 - ((-ma2 + mmb2 + sman)*(mk2 - mj2 + sman))/(2*sman) + (kappa(sman,ma2,mmb2)*kappa(sman,mk2,mj2))/(2*sman)

! dQ2 is necessary for the substitution from q2 to q2X
      dQ2=(q2max-q2min)/(q2maxAt1-q2minAt1)

! q2 corresponds to Q2 for the integration limits evaluated at x=1. This variable has to be used in the part of the plus-distribution 
! which is proportional to the Dirac delta distribution
      q2 = (q2maxAt1-q2minAt1)*xi(2) + q2minAt1

! q2X is generated from the x-dependent integration limits
      q2X = dQ2*(q2-q2minAt1)+q2min

! jacobi factor 
      intfac = (q2maxAt1-q2minAt1)*(1d0-x0)/(8*Pi*kappa(sman,ma2,mmb2))

      if (q2X.gt.0d0) then
      if (xdip.le.(-q2X+ma2+mj2)/(2*dsqrt(ma2)*(dsqrt(ma2)-dsqrt(q2X)))) then
      write(*,*) "A higher value for x0 has to be chosen!"
      endif
      endif

! Define other variables that occur often and depend on x and Q2 
      q2bar1 = q2 - ma2 -mj2
      q2barX = q2X - ma2 - mj2

      sqrtLambda1 = kappa(q2,mj2,ma2)
      sqrtLambdaX = kappa(q2X,mj2,ma2)

      Rdip = dsqrt((q2barX + 2*ma2*xdip)**2 - 4*ma2*q2X*xdip**2)/sqrtLambdaX

      sdip = ma2 + mmb2 + 1d0/Rdip*(xdip*(sman-ma2-mmb2)+(q2barX+2*ma2*xdip)/(2*q2X)*(mmb2-mk2+q2X))
     &             - (q2X + ma2 -mj2)/(2*q2X)*(mmb2-mk2+q2X)

      etaA1 = ma2/(-q2bar1)
      etaAX = ma2/(-q2barX)

      etaJ1 = mj2/(-q2bar1)
      etaJX = mj2/(-q2barX)

      vDip1 = sqrtLambda1/(-q2bar1)
      vDipX = sqrtLambdaX/(-q2barX)

! Integration limits of z. These are useful for the definition of A(x) and B(x)
      zP = (1d0-xdip)/2d0*(-q2barX+ sqrtLambdaX*Rdip)/(xdip*mj2- q2barX*(1d0-xdip))
      zM = (1d0-xdip)/2d0*(-q2barX- sqrtLambdaX*Rdip)/(xdip*mj2- q2barX*(1d0-xdip))

! The following two variables are also global
      AdipX = (zP-zM)/(1d0-xdip+zM)
      BdipX = (zP-zM)/zM

      if (mj2.eq.0d0) then
      rhoPlus = dsqrt(1+ 4*etaA1*(xdip-1)*xdip)
      yAprimePlus = ((3 - 4*xdip)*etaA1 - 1)/rhoPlus**3
      calAPlus = 2*((1 - etaA1)*xdip-2)/(rhoPlus*(2*xdip - 3 - rhoPlus)) ! \mathcal{A}
      Aplus = 1d0/((1-xdip)*calAPlus)

      yBprimePlus = etaA1*(1 - 2*xdip)/rhoPlus**3
      calBPlus = 2*etaA1*xdip/(rhoPlus**2 + rhoPlus)  ! \mathcal{B}
      Bplus = 1d0/((1-xdip)*calBPlus)

      else
      Adip1 = 2*vdip1/(2*etaJ1+1-vdip1)
      Bdip1 = 2*vdip1/(1-vdip1)
      endif

      end

      subroutine stsT2xx_SetDipoleKinematicsIEIS(xi, x0,sman,ma2,mmb2,m12,m22,intfac)
      implicit none 
#include "stsT2xx_Model.h"
#include "stsT2xx_DipoleKinematics.h"
      double precision xi(2),ma2,mmb2,m12,m22,intfac,kappa,x0,sman
      double precision q2max,q2min,q2maxAt1,q2minAt1

! Sum of squared masses for calculation of tdip/udip
      m2Sum = ma2 + mmb2 + m12 + m22

! set squared CM entergy for x=1
      sdip1 = sman

! Useful variables
      sbar = sman - ma2 -mmb2
      SqrtLambdaAB = kappa(sman,ma2,mmb2)
      d1 = sbar + 2*ma2 - SqrtLambdaAB
      d2 = sbar/SqrtLambdaAB
      Cdip = 2*SqrtLambdaAB/d1

! Determine x
      xdip = (1d0-x0)*xi(1)+x0

! Determine reduced (dipole) CM energy and associated variables
      sdip = sbar*xdip + ma2 + mmb2
      sqrtlambdaSdip = kappa(sdip,ma2,mmb2)

! integration limits for q^2
      q2max = ma2 + m12 - (sdip + ma2 - mmb2)*(sdip + m12- m22)/(2*sdip)+ sqrtlambdaSdip*kappa(sdip,m12,m22)/(2*sdip)
      q2min = ma2 + m12 - (sdip + ma2 - mmb2)*(sdip + m12- m22)/(2*sdip)- sqrtlambdaSdip*kappa(sdip,m12,m22)/(2*sdip)

! q2minAt1 and q2maxAt1 correspond to the integration limits of Q2 evaluated at x=1
      q2maxAt1 = ma2 + m12 - (sman + ma2 - mmb2)*(sman + m12- m22)/(2*sman)+ SqrtLambdaAB*kappa(sman,m12,m22)/(2*sman)
      q2minAt1 = ma2 + m12 - (sman + ma2 - mmb2)*(sman + m12- m22)/(2*sman)- SqrtLambdaAB*kappa(sman,m12,m22)/(2*sman)

! dQ2 is necessary for the substitution from q2 to q2X
      dQ2=(q2max-q2min)/(q2maxAt1-q2minAt1)

! q2 corresponds to q^2 for the integration limits evaluated at x=1. This variable has to be used in the part of the plus-distribution 
! which is proportional to the Dirac delta distribution
      q2 = (q2maxAt1-q2minAt1)*xi(2) + q2minAt1

! q2X is generated from the x-dependent integration limits
      q2X = dQ2*(q2-q2minAt1)+q2min

! jacobi factor (the term 1/lambda is defined here into the "plus"-distribution as one has to distinguish between sdip and sman)
      intfac = (1d0-x0)*(q2max-q2min)/(8*Pi)

      end