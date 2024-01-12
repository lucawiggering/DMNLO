      double precision function stsT2xx_M2DipoleReal3g()
      implicit none
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"

      double precision kappa,M2,sdip,tdip,udip,yijk,zi,zj,M2_S,M2_8S,M2_8A,tree,treeT1T2,colorTree
      double precision P2ab,xiab,SqrtlambdaAB,Rdip,xija
      double precision stsT2gg_M2TreeDipoleT1T2, stsT2gg_M2TreeDipoleT1T3,stsT2gg_M2TreeDipoleT2T3,stsT2gg_M2TreeDipole
      double precision stsT2xx_M2Dipole_SpinTensor_D12_a,stsT2xx_M2Dipole_SpinTensor_D12_3,stsT2xx_M2Dipole_SpinTensor_D23_b
      double precision stsT2xx_M2Dipole_SpinTensor_D13_2,stsT2xx_M2Dipole_SpinTensor_D23_1,stsT2xx_M2Dipole_SpinTensor_D12_b
      double precision stsT2xx_M2Dipole_SpinTensor_D23_a,stsT2xx_M2Dipole_SpinTensor_D13_a,stsT2xx_M2Dipole_SpinTensor_D13_b

      M2 = MSf2(iSf1,itt1,igen1)

! Definition of colour charge operators: T1 = T_Squark, T2 = T_Antisquark, T3 = GluonA, T4 = GluonB (A (a) and B (b) refer to colour index)
! Definition of momenta: squark(p1)+antisquark(p2) --> g^a(k1)+g^b(k2)+g^c(k3)

! *******************************************************************************************************************************************
! **************************************************** FE-FS ********************************************************************************
! *******************************************************************************************************************************************

! ---------- D_{12,3}  --------    g-> g(k1)+g(k2) -----------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------------------------------------------------

      yijk = k1k2/(k1k2 + k1k3 + k2k3)
      zi = k1k3/(k1k3 + k2k3)
      zj = 1-zi

      sdip = sman
      tdip = M2 - 2*(p1k1 + p1k2 + (p1k3*yijk)/(-1 + yijk))  ! tdip=(\tilde{p}_{ij}-p1)^2
      udip = M2 - 2*(p2k1 + p2k2 + (p2k3*yijk)/(-1 + yijk))  ! udip=(\tilde{p}_{ij}-p2)^2

      colorTree = ((CF-CV)*stsT2gg_M2TreeDipole(sdip,tdip,udip)+stsT2gg_M2TreeDipoleT1T2(sdip,tdip,udip))/CV  ! T3 T4/CA

! diagonal part
      stsT2xx_M2DipoleReal3g = -(2*CV*gs**2*colorTree*(-2 + 1/(1 - (1 - yijk)*zi) + 1/(1 - (1 - yijk)*zj)))/k1k2
    
! spin correlation tensor
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g + (-2*CV*gs**2*stsT2xx_M2Dipole_SpinTensor_D12_3(sdip,tdip,udip))/k1k2**2

! ---------- D_{13,2}  --------    g-> g(k1)+g(k3) -----------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------------------------------------------------

      yijk = k1k3/(k1k2 + k1k3 + k2k3)
      zi = k1k2/(k1k2 + k2k3)
      zj = 1-zi

      sdip = sman
      tdip = M2 - 2*(p1k1 + p1k3 + (p1k2*yijk)/(-1 + yijk))  ! tdip=(\tilde{p}_{ij}-p1)^2
      udip = M2 - 2*(p2k1 + p2k3 + (p2k2*yijk)/(-1 + yijk))  ! udip=(\tilde{p}_{ij}-p2)^2

      colorTree = ((CF-CV)*stsT2gg_M2TreeDipole(sdip,tdip,udip)+stsT2gg_M2TreeDipoleT1T2(sdip,tdip,udip))/CV  ! T3 T4/CA

! diagonal part
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (2*CV*gs**2*colorTree*(-2 + 1/(1 - (1 - yijk)*zi) + 1/(1 - (1 - yijk)*zj)))/k1k3

! spin correlation tensor
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g + (-2*CV*gs**2*stsT2xx_M2Dipole_SpinTensor_D13_2(sdip,tdip,udip))/k1k3**2


! ---------- D_{23,1}  --------    g-> g(k2)+g(k3) -----------------------------------------------------------------------------------------

      yijk = k2k3/(k1k2 + k1k3 + k2k3)
      zi = k1k2/(k1k2 + k1k3)
      zj = 1-zi

      sdip = sman
      tdip = M2 - 2*(p1k2 + p1k3 + (p1k1*yijk)/(-1 + yijk))  ! tdip=(\tilde{p}_{ij}-p1)^2
      udip = M2 - 2*(p2k2 + p2k3 + (p2k1*yijk)/(-1 + yijk))  ! udip=(\tilde{p}_{ij}-p2)^2

      colorTree = ((CF-CV)*stsT2gg_M2TreeDipole(sdip,tdip,udip)+stsT2gg_M2TreeDipoleT1T2(sdip,tdip,udip))/CV  ! T3 T4/CA

! diagonal part 
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (2*CV*gs**2*colorTree*(-2 + 1/(1 - (1 - yijk)*zi) + 1/(1 - (1 - yijk)*zj)))/k2k3

! spin correlation tensor
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g + (-2*CV*gs**2*stsT2xx_M2Dipole_SpinTensor_D23_1(sdip,tdip,udip))/k2k3**2

! *******************************************************************************************************************************************
! **************************************************** IE-IS ********************************************************************************
! *******************************************************************************************************************************************
      SqrtlambdaAB = kappa(sman,M2,M2)     

!---------- D^{a1,b}  --------    \tilde{q} -> \tilde{q}(p_a/p_1)+g(k1) ----------------------------------------------------------

      xiab = (-p1k1 + p1p2 - p2k1)/p1p2

      if (xiab.gt.x0IEISg) then
      P2ab = 2*M2 - 2*p1k1 + 2*p1p2 - 2*p2k1

      sdip = P2ab
      tdip = M2 + 2*(k1k2 - p1k2 - p2k2) - 2*p2k2 - ((2*M2 + 2*p1p2 + P2ab - 2*p2k1)*
     &     ((2*k1k2*M2 - P2ab*p2k2 - 2*M2*(p1k2 + p2k2))*SqrtlambdaAB + (-2*M2*p1k2 + 2*p1p2*p2k2)*kappa(P2ab,M2,M2)))/
     &   (P2ab*(3*M2 + p1p2 - p2k1)*SqrtlambdaAB + 2*(M2**2 - M2*p1k1 + p1p2*(-p1p2 + p2k1))*kappa(P2ab,M2,M2))         ! tdip = (p2-\tilde{k}_2)^2 = M2 - 2 p2 Lambda k2
      udip =  M2 + 2*(k1k3 - p1k3 - p2k3) - 2*p2k3 - ((2*M2 + 2*p1p2 + P2ab - 2*p2k1)*
     &     ((2*k1k3*M2 - P2ab*p2k3 - 2*M2*(p1k3 + p2k3))*SqrtlambdaAB + (-2*M2*p1k3 + 2*p1p2*p2k3)*kappa(P2ab,M2,M2)))/
     &   (P2ab*(3*M2 + p1p2 - p2k1)*SqrtlambdaAB + 2*(M2**2 - M2*p1k1 + p1p2*(-p1p2 + p2k1))*kappa(P2ab,M2,M2))         ! udip = (p2-\tilde{k}_3)^2 = M2 - 2 p2 Lambda k3

      colorTree = stsT2gg_M2TreeDipoleT1T2(sdip,tdip,udip)/CF  ! T1 T2/CF

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*gs**2*colorTree*(-2 + 2/(1 - xiab) - (M2*xiab)/p1k1))/(p1k1*xiab)

      endif

!---------- D^{a2,b}  --------    \tilde{q} -> \tilde{q}(p_a/p_1)+g(k2) ----------------------------------------------------------

      xiab = (-p1k2 + p1p2 - p2k2)/p1p2

      if (xiab.gt.x0IEISg) then
      P2ab = 2*M2 - 2*p1k2 + 2*p1p2 - 2*p2k2

      sdip = P2ab
      tdip = M2 + 2*(k1k2 - p1k1 - p2k1) - 2*p2k1 - ((2*M2 + 2*p1p2 + P2ab - 2*p2k2)*
     &     ((2*k1k2*M2 - P2ab*p2k1 - 2*M2*(p1k1 + p2k1))*SqrtlambdaAB + (-2*M2*p1k1 + 2*p1p2*p2k1)*kappa(P2ab,M2,M2)))/
     &   (P2ab*(3*M2 + p1p2 - p2k2)*SqrtlambdaAB + 2*(M2**2 - M2*p1k2 + p1p2*(-p1p2 + p2k2))*kappa(P2ab,M2,M2)) ! tdip = (p2-\tilde{k}_1)^2 = M2 - 2 p2 Lambda k1
      udip = M2 + 2*(k2k3 - p1k3 - p2k3) - 2*p2k3 - ((2*M2 + 2*p1p2 + P2ab - 2*p2k2)*
     &     ((2*k2k3*M2 - P2ab*p2k3 - 2*M2*(p1k3 + p2k3))*SqrtlambdaAB + (-2*M2*p1k3 + 2*p1p2*p2k3)*kappa(P2ab,M2,M2)))/
     &   (P2ab*(3*M2 + p1p2 - p2k2)*SqrtlambdaAB + 2*(M2**2 - M2*p1k2 + p1p2*(-p1p2 + p2k2))*kappa(P2ab,M2,M2)) ! udip = (p2-\tilde{k}_3)^2 = M2 - 2 p2 Lambda k3

      colorTree = stsT2gg_M2TreeDipoleT1T2(sdip,tdip,udip)/CF  ! T1 T2/CF

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 + 2/(1 - xiab) - (M2*xiab)/p1k2))/(p1k2*xiab)

      endif

!---------- D^{a3,b}  --------    \tilde{q} -> \tilde{q}(p_a/p_1)+g(k3) ----------------------------------------------------------

      xiab = (-p1k3 + p1p2 - p2k3)/p1p2

      if (xiab.gt.x0IEISg) then
      P2ab = 2*M2 - 2*p1k3 + 2*p1p2 - 2*p2k3

      sdip = P2ab
      tdip =  M2 + 2*(k1k3 - p1k1 - p2k1) - 2*p2k1 - ((2*M2 + 2*p1p2 + P2ab - 2*p2k3)*
     &     ((2*k1k3*M2 - P2ab*p2k1 - 2*M2*(p1k1 + p2k1))*SqrtlambdaAB + (-2*M2*p1k1 + 2*p1p2*p2k1)*kappa(P2ab,M2,M2)))/
     &   (P2ab*(3*M2 + p1p2 - p2k3)*SqrtlambdaAB + 2*(M2**2 - M2*p1k3 + p1p2*(-p1p2 + p2k3))*kappa(P2ab,M2,M2))  ! tdip = (p2-\tilde{k}_1)^2 = M2 - 2 p2 Lambda k1
      udip = M2 + 2*(k2k3 - p1k2 - p2k2) - 2*p2k2 - ((2*M2 + 2*p1p2 + P2ab - 2*p2k3)*
     &     ((2*k2k3*M2 - P2ab*p2k2 - 2*M2*(p1k2 + p2k2))*SqrtlambdaAB + (-2*M2*p1k2 + 2*p1p2*p2k2)*kappa(P2ab,M2,M2)))/
     &   (P2ab*(3*M2 + p1p2 - p2k3)*SqrtlambdaAB + 2*(M2**2 - M2*p1k3 + p1p2*(-p1p2 + p2k3))*kappa(P2ab,M2,M2))  ! udip = (p2-\tilde{k}_2)^2 = M2 - 2 p2 Lambda k2

      colorTree = stsT2gg_M2TreeDipoleT1T2(sdip,tdip,udip)/CF  ! T1 T2/CF

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 + 2/(1 - xiab) - (M2*xiab)/p1k3))/(p1k3*xiab)

      endif


!---------- D^{b1,a}  --------    \tilde{q}^{\ast} -> \tilde{q}^{\ast}(p_b/p_2)+g(k1) ----------------------------------------------------------

      xiab = (-p1k1 + p1p2 - p2k1)/p1p2

      if (xiab.gt.x0IEISg) then
      P2ab = 2*M2 - 2*p1k1 + 2*p1p2 - 2*p2k1

      sdip = P2ab
      tdip = M2 - 2*p1k2 + 2*(k1k2 - p1k2 - p2k2) - ((2*M2 - 2*p1k1 + 2*p1p2 + P2ab)*
     &     ((2*k1k2*M2 - p1k2*P2ab - 2*M2*(p1k2 + p2k2))*SqrtlambdaAB + 2*(p1k2*p1p2 - M2*p2k2)*kappa(P2ab,M2,M2)))/
     &   ((3*M2 - p1k1 + p1p2)*P2ab*SqrtlambdaAB + 2*(M2**2 + (p1k1 - p1p2)*p1p2 - M2*p2k1)*kappa(P2ab,M2,M2))    ! tdip = (p1-\tilde{k}_2)^2 = M2 - 2 p1 Lambda k2
      udip =  M2 - 2*p1k3 + 2*(k1k3 - p1k3 - p2k3) - ((2*M2 - 2*p1k1 + 2*p1p2 + P2ab)*
     &     ((2*k1k3*M2 - p1k3*P2ab - 2*M2*(p1k3 + p2k3))*SqrtlambdaAB + 2*(p1k3*p1p2 - M2*p2k3)*kappa(P2ab,M2,M2)))/
     &   ((3*M2 - p1k1 + p1p2)*P2ab*SqrtlambdaAB + 2*(M2**2 + (p1k1 - p1p2)*p1p2 - M2*p2k1)*kappa(P2ab,M2,M2))    ! udip = (p1-\tilde{k}_3)^2 = M2 - 2 p1 Lambda k3

      colorTree = stsT2gg_M2TreeDipoleT1T2(sdip,tdip,udip)/CF  ! T1 T2/CF

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 + 2/(1 - xiab) - (M2*xiab)/p2k1))/(p2k1*xiab)

      endif

!---------- D^{b2,a}  --------    \tilde{q}^{\ast} -> \tilde{q}^{\ast}(p_b/p_2)+g(k2) ----------------------------------------------------------

      xiab = (-p1k2 + p1p2 - p2k2)/p1p2

      if (xiab.gt.x0IEISg) then
      P2ab = 2*M2 - 2*p1k2 + 2*p1p2 - 2*p2k2

      sdip = P2ab
      tdip = M2 - 2*p1k1 + 2*(k1k2 - p1k1 - p2k1) - ((2*M2 - 2*p1k2 + 2*p1p2 + P2ab)*
     &     ((2*k1k2*M2 - p1k1*P2ab - 2*M2*(p1k1 + p2k1))*SqrtlambdaAB + 2*(p1k1*p1p2 - M2*p2k1)*kappa(P2ab,M2,M2)))/
     &   ((3*M2 - p1k2 + p1p2)*P2ab*SqrtlambdaAB + 2*(M2**2 + (p1k2 - p1p2)*p1p2 - M2*p2k2)*kappa(P2ab,M2,M2))  ! tdip = (p1-\tilde{k}_1)^2 = M2 - 2 p1 Lambda k1
      udip = M2 - 2*p1k3 + 2*(k2k3 - p1k3 - p2k3) - ((2*M2 - 2*p1k2 + 2*p1p2 + P2ab)*
     &     ((2*k2k3*M2 - p1k3*P2ab - 2*M2*(p1k3 + p2k3))*SqrtlambdaAB + 2*(p1k3*p1p2 - M2*p2k3)*kappa(P2ab,M2,M2)))/
     &   ((3*M2 - p1k2 + p1p2)*P2ab*SqrtlambdaAB + 2*(M2**2 + (p1k2 - p1p2)*p1p2 - M2*p2k2)*kappa(P2ab,M2,M2))  ! udip = (p1-\tilde{k}_3)^2 = M2 - 2 p1 Lambda k3

      colorTree = stsT2gg_M2TreeDipoleT1T2(sdip,tdip,udip)/CF  ! T1 T2/CF

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 + 2/(1 - xiab) - (M2*xiab)/p2k2))/(p2k2*xiab)

      endif

!---------- D^{b3,a}  --------    \tilde{q}^{\ast} -> \tilde{q}^{\ast}(p_b/p_2)+g(k3) ----------------------------------------------------------

      xiab = (-p1k3 + p1p2 - p2k3)/p1p2

      if (xiab.gt.x0IEISg) then
      P2ab = 2*M2 - 2*p1k3 + 2*p1p2 - 2*p2k3

      sdip = P2ab
      tdip = M2 - 2*p1k1 + 2*(k1k3 - p1k1 - p2k1) - ((2*M2 - 2*p1k3 + 2*p1p2 + P2ab)*
     &     ((2*k1k3*M2 - p1k1*P2ab - 2*M2*(p1k1 + p2k1))*SqrtlambdaAB + 2*(p1k1*p1p2 - M2*p2k1)*kappa(P2ab,M2,M2)))/
     &   ((3*M2 - p1k3 + p1p2)*P2ab*SqrtlambdaAB + 2*(M2**2 + (p1k3 - p1p2)*p1p2 - M2*p2k3)*kappa(P2ab,M2,M2))     ! tdip = (p1-\tilde{k}_1)^2 = M2 - 2 p1 Lambda k1
      udip = M2 - 2*p1k2 + 2*(k2k3 - p1k2 - p2k2) - ((2*M2 - 2*p1k3 + 2*p1p2 + P2ab)*
     &     ((2*k2k3*M2 - p1k2*P2ab - 2*M2*(p1k2 + p2k2))*SqrtlambdaAB + 2*(p1k2*p1p2 - M2*p2k2)*kappa(P2ab,M2,M2)))/
     &   ((3*M2 - p1k3 + p1p2)*P2ab*SqrtlambdaAB + 2*(M2**2 + (p1k3 - p1p2)*p1p2 - M2*p2k3)*kappa(P2ab,M2,M2))     ! udip = (p1-\tilde{k}_2)^2 = M2 - 2 p1 Lambda k2

      colorTree = stsT2gg_M2TreeDipoleT1T2(sdip,tdip,udip)/CF  ! T1 T2/CF

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 + 2/(1 - xiab) - (M2*xiab)/p2k3))/(p2k3*xiab)

      endif

! *******************************************************************************************************************************************
! **************************************************** IE-FS ********************************************************************************
! *******************************************************************************************************************************************
!---------- D^{a2}_1  --------    \tilde{q} -> \tilde{q}(p_a/p_1)+g(k2) ----------------------------------------------------------

      xija = (-k1k2 + p1k1 + p1k2)/(p1k1 + p1k2)
      if (xija.gt.x0FEISg) then

      zi = p1k2/(p1k1 + p1k2)

      tdip = M2 - 2*p2k3 !Q^2=tdip=(p2-k3)^2, Q=p2-k3
      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p2k3)*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))
      udip = M2 - 2*((p1k3*xija)/Rdip + p2k3*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))

      colorTree = stsT2gg_M2TreeDipoleT1T3(sdip,tdip,udip)/CF  ! T1 T3/CF

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 - (M2*xija)/p1k2 + 2/(1 - xija + zi)))/(p1k2*xija)

      endif

!---------- D^{a3}_1  --------    \tilde{q} -> \tilde{q}(p_a/p_1)+g(k3) ----------------------------------------------------------

      xija = (-k1k3 + p1k1 + p1k3)/(p1k1 + p1k3)
      if (xija.gt.x0FEISg) then

      zi = p1k3/(p1k1 + p1k3)

      tdip = M2 - 2*p2k2 !Q^2=tdip=(p2-k2)^2, Q=p2-k2
      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p2k2)*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))
      udip = M2 - 2*((p1k2*xija)/Rdip + p2k2*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))

      colorTree = stsT2gg_M2TreeDipoleT1T3(sdip,tdip,udip)/CF   ! T1 T3/CF
 
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 - (M2*xija)/p1k3 + 2/(1 - xija + zi)))/(p1k3*xija)

      endif

!---------- D^{a1}_2  --------    \tilde{q} -> \tilde{q}(p_a/p_1)+g(k1) ----------------------------------------------------------

      xija = (-k1k2 + p1k1 + p1k2)/(p1k1 + p1k2)
      if (xija.gt.x0FEISg) then

      zi = p1k1/(p1k1 + p1k2)

      tdip = M2 - 2*p2k3 !Q^2=tdip=(p2-k3)^2, Q=p2-k3
      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p2k3)*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))
      udip = M2 - 2*((p1k3*xija)/Rdip + p2k3*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))

      colorTree = stsT2gg_M2TreeDipoleT1T3(sdip,tdip,udip)/CF  ! T1 T3/CF

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 - (M2*xija)/p1k1 + 2/(1 - xija + zi)))/(p1k1*xija)

      endif

!---------- D^{a2}_3  --------    \tilde{q} -> \tilde{q}(p_a/p_1)+g(k2) ----------------------------------------------------------

      xija = (-k2k3 + p1k2 + p1k3)/(p1k2 + p1k3)
      if (xija.gt.x0FEISg) then

      zi = p1k2/(p1k2 + p1k3)

      udip = M2 - 2*p2k1 !Q^2=udip=(p2-k1)^2, Q=p2-k1
      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p2k1)*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))
      tdip = M2 - 2*((p1k1*xija)/Rdip + p2k1*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))

      colorTree = stsT2gg_M2TreeDipoleT2T3(sdip,tdip,udip)/CF  ! T1 T4/CF = T2 T3/CF

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 - (M2*xija)/p1k2 + 2/(1 - xija + zi)))/(p1k2*xija)

      endif

!---------- D^{a3}_2  --------    \tilde{q} -> \tilde{q}(p_a/p_1)+g(k3) ----------------------------------------------------------

      xija = (-k2k3 + p1k2 + p1k3)/(p1k2 + p1k3)
      if (xija.gt.x0FEISg) then

      zi = p1k3/(p1k2 + p1k3)

      udip = M2 - 2*p2k1 !Q^2=udip=(p2-k1)^2, Q=p2-k1
      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p2k1)*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))
      tdip = M2 - 2*((p1k1*xija)/Rdip + p2k1*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))

      colorTree = stsT2gg_M2TreeDipoleT2T3(sdip,tdip,udip)/CF  ! T1 T4/CF = T2 T3/CF

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 - (M2*xija)/p1k3 + 2/(1 - xija + zi)))/(p1k3*xija)

      endif

!---------- D^{a1}_3  --------    \tilde{q} -> \tilde{q}(p_a/p_1)+g(k1) ----------------------------------------------------------

      xija = (-k1k3 + p1k1 + p1k3)/(p1k1 + p1k3)
      if (xija.gt.x0FEISg) then

      zi =  p1k1/(p1k1 + p1k3)

      udip = M2 - 2*p2k2 !Q^2=udip=(p2-k2)^2, Q=p2-k2
      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p2k2)*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))
      tdip = M2 - 2*((p1k2*xija)/Rdip + p2k2*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))

      colorTree = stsT2gg_M2TreeDipoleT2T3(sdip,tdip,udip)/CF  ! T1 T4/CF = T2 T3/CF

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 - (M2*xija)/p1k1 + 2/(1 - xija + zi)))/(p1k1*xija)

      endif


!---------- D^{b2}_1  --------    \tilde{q} -> \tilde{q}(p_b/p_2)+g(k2) ----------------------------------------------------------

      xija = (-k1k2 + p2k1 + p2k2)/(p2k1 + p2k2)
      if (xija.gt.x0FEISg) then

      zi =  p2k2/(p2k1 + p2k2)

      tdip = M2 - 2*p1k3 !Q^2=tdip=(p1-k3)^2, Q=p1-k3
      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p1k3)*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))
      udip = M2 - 2*((p2k3*xija)/Rdip + p1k3*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))

      colorTree = stsT2gg_M2TreeDipoleT1T3(sdip,tdip,udip)/CF  ! T2 T4/CF=T1 T3/CF

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 - (M2*xija)/p2k2 + 2/(1 - xija + zi)))/(p2k2*xija)

      endif


!---------- D^{b3}_1  --------    \tilde{q} -> \tilde{q}(p_b/p_2)+g(k3) ----------------------------------------------------------

      xija = (-k1k3 + p2k1 + p2k3)/(p2k1 + p2k3)
      if (xija.gt.x0FEISg) then

      zi =  p2k3/(p2k1 + p2k3)

      tdip = M2 - 2*p1k2 !Q^2=tdip=(p1-k2)^2, Q=p1-k2
      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p1k2)*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))
      udip = M2 - 2*((p2k2*xija)/Rdip + p1k2*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))

      colorTree = stsT2gg_M2TreeDipoleT1T3(sdip,tdip,udip)/CF !  T2 T4/CF = T1 T3/CF  

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 - (M2*xija)/p2k3 + 2/(1 - xija + zi)))/(p2k3*xija)

      endif

!---------- D^{b1}_2  --------    \tilde{q} -> \tilde{q}(p_b/p_2)+g(k1) ----------------------------------------------------------

      xija = (-k1k2 + p2k1 + p2k2)/(p2k1 + p2k2)
      if (xija.gt.x0FEISg) then

      zi =  p2k1/(p2k1 + p2k2)

      tdip = M2 - 2*p1k3 !Q^2=tdip=(p1-k3)^2, Q=p1-k3
      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p1k3)*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))
      udip = M2 - 2*((p2k3*xija)/Rdip + p1k3*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))

      colorTree = stsT2gg_M2TreeDipoleT1T3(sdip,tdip,udip)/CF !  T2 T4/CF = T1 T3/CF  

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 - (M2*xija)/p2k1 + 2/(1 - xija + zi)))/(p2k1*xija)

      endif

!---------- D^{b2}_3  --------    \tilde{q} -> \tilde{q}(p_b/p_2)+g(k2) ----------------------------------------------------------

      xija = (-k2k3 + p2k2 + p2k3)/(p2k2 + p2k3)
      if (xija.gt.x0FEISg) then

      zi =  p2k2/(p2k2 + p2k3)

      udip = M2 - 2*p1k1 !Q^2=udip=(p1-k1)^2, Q=p1-k1
      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p1k1)*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))
      tdip = M2 - 2*((p2k1*xija)/Rdip + p1k1*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))

      colorTree = stsT2gg_M2TreeDipoleT2T3(sdip,tdip,udip)/CF !  T2 T3/CF  

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 - (M2*xija)/p2k2 + 2/(1 - xija + zi)))/(p2k2*xija)

      endif

!---------- D^{b3}_2  --------    \tilde{q} -> \tilde{q}(p_b/p_2)+g(k3) ----------------------------------------------------------

      xija = (-k2k3 + p2k2 + p2k3)/(p2k2 + p2k3)
      if (xija.gt.x0FEISg) then

      zi =  p2k3/(p2k2 + p2k3)

      udip = M2 - 2*p1k1 !Q^2=udip=(p1-k1)^2, Q=p1-k1
      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p1k1)*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))
      tdip = M2 - 2*((p2k1*xija)/Rdip + p1k1*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))

      colorTree = stsT2gg_M2TreeDipoleT2T3(sdip,tdip,udip)/CF !  T2 T3/CF 

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 - (M2*xija)/p2k3 + 2/(1 - xija + zi)))/(p2k3*xija)

      endif

!---------- D^{b1}_3  --------    \tilde{q} -> \tilde{q}(p_b/p_2)+g(k1) ----------------------------------------------------------

      xija = (-k1k3 + p2k1 + p2k3)/(p2k1 + p2k3)
      if (xija.gt.x0FEISg) then

      zi =  p2k1/(p2k1 + p2k3)

      udip = M2 - 2*p1k2 !Q^2=udip=(p1-k2)^2, Q=p1-k2
      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p1k2)*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))
      tdip = M2 - 2*((p2k2*xija)/Rdip + p1k2*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))

      colorTree = stsT2gg_M2TreeDipoleT2T3(sdip,tdip,udip)/CF !  T2 T3/CF   

      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (CF*colorTree*gs**2*(-2 - (M2*xija)/p2k1 + 2/(1 - xija + zi)))/(p2k1*xija)

      endif


! *******************************************************************************************************************************************
! **************************************************** FE-IS ********************************************************************************
! *******************************************************************************************************************************************
!---------- D_{12}^a  --------   g -> g(k1) + g(k2) ----------------------------------------------------------

      xija = (-k1k2 + p1k1 + p1k2)/(p1k1 + p1k2)
      if (xija.gt.x0FEISg) then

      tdip = M2 - 2*p2k3 !Q^2=tdip=(p2-k3)^2, Q=p2-k3
      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p2k3)*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))
      udip = M2 - 2*((p1k3*xija)/Rdip + p2k3*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))  

      zi =  p1k1/(p1k1 + p1k2)

      colorTree = stsT2gg_M2TreeDipoleT1T3(sdip,tdip,udip)/CV  ! T1 T3/CA

! diagonal part
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (2*colorTree*CV*gs**2*(-2 + 1/(2 - xija - zi) + 1/(1 - xija + zi)))/(k1k2*xija)

! spin correlation tensor
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (2*CV*gs**2*stsT2xx_M2Dipole_SpinTensor_D12_a(sdip,tdip,udip))/(k1k2**2*xija)

      endif

!---------- D_{23}^a  --------   g -> g(k2) + g(k3) ----------------------------------------------------------

      xija = (-k2k3 + p1k2 + p1k3)/(p1k2 + p1k3)
      if (xija.gt.x0FEISg) then

      tdip = M2 - 2*p2k1 !Q^2=tdip=(p2-k1)^2, Q=p2-k1
      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p2k1)*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))
      udip = M2 - 2*((p1k1*xija)/Rdip + p2k1*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))

      zi =  p1k2/(p1k2 + p1k3)

      colorTree = stsT2gg_M2TreeDipoleT1T3(sdip,tdip,udip)/CV   ! T1 T3/CA

! diagonal part
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (2*colorTree*CV*gs**2*(-2 + 1/(2 - xija - zi) + 1/(1 - xija + zi)))/(k2k3*xija)

! spin correlation tensor
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (2*CV*gs**2*stsT2xx_M2Dipole_SpinTensor_D23_a(sdip,tdip,udip))/(k2k3**2*xija)

      endif

!---------- D_{13}^a  --------   g -> g(k1) + g(k3) ----------------------------------------------------------

      xija = (-k1k3 + p1k1 + p1k3)/(p1k1 + p1k3)
      if (xija.gt.x0FEISg) then

      tdip = M2 - 2*p2k2 !Q^2=tdip=(p2-k2)^2, Q=p2-k2
      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p2k2)*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))
      udip = M2 - 2*((p1k2*xija)/Rdip + p2k2*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))

      zi =  p1k1/(p1k1 + p1k3)

      colorTree = stsT2gg_M2TreeDipoleT1T3(sdip,tdip,udip)/CV    ! T1 T3/CA

! diagonal part
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (2*colorTree*CV*gs**2*(-2 + 1/(2 - xija - zi) + 1/(1 - xija + zi)))/(k1k3*xija)

! spin correlation tensor
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (2*CV*gs**2*stsT2xx_M2Dipole_SpinTensor_D13_a(sdip,tdip,udip))/(k1k3**2*xija)

      endif

!---------- D_{12}^b  --------   g -> g(k1) + g(k2) ----------------------------------------------------------

      xija = (-k1k2 + p2k1 + p2k2)/(p2k1 + p2k2)
      if (xija.gt.x0FEISg) then

      udip = M2 - 2*p1k3 !Q^2=udip=(p1-k3)^2, Q=p1-k3
      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p1k3)*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))
      tdip = M2 - 2*((p2k3*xija)/Rdip + p1k3*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))

      zi =  p2k1/(p2k1 + p2k2)

      colorTree = stsT2gg_M2TreeDipoleT2T3(sdip,tdip,udip)/CV   ! T2 T3/CA

! diagonal part
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (2*colorTree*CV*gs**2*(-2 + 1/(2 - xija - zi) + 1/(1 - xija + zi)))/(k1k2*xija)

! spin correlation tensor
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (2*CV*gs**2*stsT2xx_M2Dipole_SpinTensor_D12_b(sdip,tdip,udip))/(k1k2**2*xija)
      
      endif

!---------- D_{13}^b  --------   g -> g(k1) + g(k3) ----------------------------------------------------------

      xija = (-k1k3 + p2k1 + p2k3)/(p2k1 + p2k3)
      if (xija.gt.x0FEISg) then

      udip = M2 - 2*p1k2 !Q^2=udip=(p1-k2)^2, Q=p1-k2
      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p1k2)*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))
      tdip = M2 - 2*((p2k2*xija)/Rdip + p1k2*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))

      zi =  p2k1/(p2k1 + p2k3)

      colorTree = stsT2gg_M2TreeDipoleT2T3(sdip,tdip,udip)/CV   ! T2 T3/CA

! diagonal part
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (2*colorTree*CV*gs**2*(-2 + 1/(2 - xija - zi) + 1/(1 - xija + zi)))/(k1k3*xija)

! spin correlation tensor
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (2*CV*gs**2*stsT2xx_M2Dipole_SpinTensor_D13_b(sdip,tdip,udip))/(k1k3**2*xija)
    
      endif

!---------- D_{23}^b  --------   g -> g(k2) + g(k3) ----------------------------------------------------------

      xija = (-k2k3 + p2k2 + p2k3)/(p2k2 + p2k3)
      if (xija.gt.x0FEISg) then

      udip = M2 - 2*p1k1 !Q^2=udip=(p1-k1)^2, Q=p1-k1
      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p1k1)*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))
      tdip = M2 - 2*((p2k1*xija)/Rdip + p1k1*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))

      zi =  p2k2/(p2k2 + p2k3)

      colorTree = stsT2gg_M2TreeDipoleT2T3(sdip,tdip,udip)/CV   ! T2 T3/CA

! diagonal part
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (2*colorTree*CV*gs**2*(-2 + 1/(2 - xija - zi) + 1/(1 - xija + zi)))/(k2k3*xija)

! spin correlation tensor
      stsT2xx_M2DipoleReal3g = stsT2xx_M2DipoleReal3g - (2*CV*gs**2*stsT2xx_M2Dipole_SpinTensor_D23_b(sdip,tdip,udip))/(k2k3**2*xija)

      endif

      end 





      double precision function stsT2xx_M2DipoleRealQ()
      implicit none
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"

      double precision kappa,M2,sdip,tdip,udip,zi, M2_S,M2_8S,M2_8A,tree
      double precision treeT1T2,treeT1T3,treeT1T4,treeT2T3,treeT2T4
      double precision colorTree,Rdip,xija,xiab,P2ab,SqrtlambdaAB,yijk,zj
      double precision stsT2gg_M2TreeDipole, stsT2gg_M2TreeDipoleT1T2
      double precision stsT2gg_M2TreeDipoleT1T3,stsT2QQbar_M2TreeDipole
      double precision stsT2xx_M2Dipole_SpinTensor_D12_3,stsT2xx_M2Dipole_SpinTensor_D12_a
      double precision stsT2xx_M2Dipole_SpinTensor_D12_b,stsT2gg_M2TreeDipoleT2T3

      M2 = MSf2(iSf1,itt1,igen1)

! Definition of colour charge operators:
!       -   For gg tree level: T1 = T_Squark, T2 = T_Antisquark, T3 = GluonA, T4 = GluonB
!       -   For QQbar tree level: T1 = T_Squark, T2 = T_Antisquark, T3 = Quark, T4 = Antiquark
! Definition of momenta: squark(t,p1)+antisquark(s,p2) --> Quark(r,k1)+Antiquark(u,k2)+g(a,k3) (s,t,r,u,a are color indices)

! *******************************************************************************************************************************************
! **************************************************** FE-FS ********************************************************************************
! *******************************************************************************************************************************************
!---------- D_{12,3}  --------  g -> q(k1)+qbar(k2)

      yijk = k1k2/(k1k2 + k1k3 + k2k3)

      sdip = sman
      tdip = M2-2/(1-yijk)*p2k3  ! tdip=(\tilde{p}_k-p2)^2
      udip = M2-2/(1-yijk)*p1k3  ! udip=(\tilde{p}_k-p1)^2

      zi = k1k3/(k1k3 + k2k3)

      colorTree = ((CF-CV)*stsT2gg_M2TreeDipole(sdip,tdip,udip)+stsT2gg_M2TreeDipoleT1T2(sdip,tdip,udip))/CV  ! T3 T4/CA

! diagonal part
      stsT2xx_M2DipoleRealQ = -(colorTree*gs**2*TF)/k1k2

! spin correlation tensor
      stsT2xx_M2DipoleRealQ =  stsT2xx_M2DipoleRealQ + (2*TF*gs**2*stsT2xx_M2Dipole_SpinTensor_D12_3(sdip,tdip,udip))/k1k2**2

! ---------- D_{13,2}  --------    q-> q(k1)+g(k3) -----------------------------------------------------------------------------------------

      yijk = k1k3/(k1k2 + k1k3 + k2k3)

      sdip = sman
      udip = M2-2/(1-yijk)*p1k2 ! udip=(\tilde{p}_k-p1)^2
      tdip = M2-2/(1-yijk)*p2k2 ! tdip=(\tilde{p}_k-p2)^2

      zi = k2k3/(k1k2 + k2k3)
      zj = 1-zi

      colorTree = stsT2QQbar_M2TreeDipole(sdip,tdip,udip)/(2*Nc)/CF  ! T3 T4/CA

      stsT2xx_M2DipoleRealQ = stsT2xx_M2DipoleRealQ - (CF*colorTree*gs**2*(-1 - zj + 2/(1 - (1 - yijk)*zj)))/k1k3

! ---------- D_{23,1}  --------    qbar-> qbar(k2)+g(k3) -----------------------------------------------------------------------------------------

      yijk = k2k3/(k1k2 + k1k3 + k2k3)

      sdip = sman
      udip = M2-2/(1-yijk)*p2k1  ! udip=(\tilde{p}_k-p2)^2
      tdip = M2-2/(1-yijk)*p1k1  ! tdip=(\tilde{p}_k-p1)^2

      zi = k1k3/(k1k2 + k1k3)
      zj = 1-zi
      
      colorTree = stsT2QQbar_M2TreeDipole(sdip,tdip,udip)/(2*Nc)/CF  ! T3 T4/CA

      stsT2xx_M2DipoleRealQ = stsT2xx_M2DipoleRealQ - (CF*colorTree*gs**2*(-1 - zj + 2/(1 - (1 - yijk)*zj)))/k2k3


! *******************************************************************************************************************************************
! **************************************************** IE-IS ********************************************************************************
! *******************************************************************************************************************************************
      SqrtlambdaAB = kappa(sman,M2,M2)     

!---------- D^{a3,b} --------
      xiab = (-p1k3 + p1p2 - p2k3)/p1p2

      if (xiab.gt.x0IEISQ) then
      P2ab = 2*M2 - 2*p1k3 + 2*p1p2 - 2*p2k3

      sdip = P2ab
      tdip =  M2 + 2*(k1k3 - p1k1 - p2k1) - 2*p2k1 - ((2*M2 + 2*p1p2 + P2ab - 2*p2k3)*
     &     ((2*k1k3*M2 - P2ab*p2k1 - 2*M2*(p1k1 + p2k1))*SqrtlambdaAB + (-2*M2*p1k1 + 2*p1p2*p2k1)*kappa(P2ab,M2,M2)))/
     &   (P2ab*(3*M2 + p1p2 - p2k3)*SqrtlambdaAB + 2*(M2**2 - M2*p1k3 + p1p2*(-p1p2 + p2k3))*kappa(P2ab,M2,M2))  ! tdip = (p2-\tilde{k}_1)^2 = M2 - 2 p2 Lambda k1
      udip = M2 + 2*(k2k3 - p1k2 - p2k2) - 2*p2k2 - ((2*M2 + 2*p1p2 + P2ab - 2*p2k3)*
     &     ((2*k2k3*M2 - P2ab*p2k2 - 2*M2*(p1k2 + p2k2))*SqrtlambdaAB + (-2*M2*p1k2 + 2*p1p2*p2k2)*kappa(P2ab,M2,M2)))/
     &   (P2ab*(3*M2 + p1p2 - p2k3)*SqrtlambdaAB + 2*(M2**2 - M2*p1k3 + p1p2*(-p1p2 + p2k3))*kappa(P2ab,M2,M2))  ! udip = (p2-\tilde{k}_2)^2 = M2 - 2 p2 Lambda k2

      colorTree = stsT2QQbar_M2TreeDipole(sdip,tdip,udip)/(2*Nc)/CF   ! T1 T2/CF

      stsT2xx_M2DipoleRealQ = stsT2xx_M2DipoleRealQ - (CF*colorTree*gs**2*(-2 + 2/(1 - xiab) - (M2*xiab)/p1k3))/(p1k3*xiab)

      endif

!---------- D^{b3,a}  --------   

      xiab = (-p1k3 + p1p2 - p2k3)/p1p2

      if (xiab.gt.x0IEISQ) then
      P2ab = 2*M2 - 2*p1k3 + 2*p1p2 - 2*p2k3

      sdip = P2ab
      tdip = M2 - 2*p1k1 + 2*(k1k3 - p1k1 - p2k1) - ((2*M2 - 2*p1k3 + 2*p1p2 + P2ab)*
     &     ((2*k1k3*M2 - p1k1*P2ab - 2*M2*(p1k1 + p2k1))*SqrtlambdaAB + 2*(p1k1*p1p2 - M2*p2k1)*kappa(P2ab,M2,M2)))/
     &   ((3*M2 - p1k3 + p1p2)*P2ab*SqrtlambdaAB + 2*(M2**2 + (p1k3 - p1p2)*p1p2 - M2*p2k3)*kappa(P2ab,M2,M2))     ! tdip = (p1-\tilde{k}_1)^2 = M2 - 2 p1 Lambda k1
      udip = M2 - 2*p1k2 + 2*(k2k3 - p1k2 - p2k2) - ((2*M2 - 2*p1k3 + 2*p1p2 + P2ab)*
     &     ((2*k2k3*M2 - p1k2*P2ab - 2*M2*(p1k2 + p2k2))*SqrtlambdaAB + 2*(p1k2*p1p2 - M2*p2k2)*kappa(P2ab,M2,M2)))/
     &   ((3*M2 - p1k3 + p1p2)*P2ab*SqrtlambdaAB + 2*(M2**2 + (p1k3 - p1p2)*p1p2 - M2*p2k3)*kappa(P2ab,M2,M2))     ! udip = (p1-\tilde{k}_2)^2 = M2 - 2 p1 Lambda k2


      colorTree = stsT2QQbar_M2TreeDipole(sdip,tdip,udip)/(2*Nc)/CF   ! T1 T2/CF

      stsT2xx_M2DipoleRealQ = stsT2xx_M2DipoleRealQ - (CF*colorTree*gs**2*(-2 + 2/(1 - xiab) - (M2*xiab)/p2k3))/(p2k3*xiab)

      endif


! *******************************************************************************************************************************************
! **************************************************** IE-FS ********************************************************************************
! *******************************************************************************************************************************************

!---------- D^{a3}_1  --------   \tilde{q} -> \tilde{q}(p_a/p_1)+g(k3) ----------------------------------------------------------

      xija = (-k1k3 + p1k1 + p1k3)/(p1k1 + p1k3)
      if (xija.gt.x0FEISQ) then

      zi = p1k3/(p1k1 + p1k3)

      tdip = M2 - 2*p2k2 !Q^2=tdip=(p2-k2)^2, Q=p2-k2
      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p2k2)*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))
      udip = M2 - 2*((p1k2*xija)/Rdip + p2k2*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))

      treeT1T3 = (2-Nc**2)/(2*Nc)*stsT2QQbar_M2TreeDipole(sdip,tdip,udip)
      colorTree = treeT1T3/CF   ! T1 T3/CF
 
      stsT2xx_M2DipoleRealQ = stsT2xx_M2DipoleRealQ - (CF*colorTree*gs**2*(-2 - (M2*xija)/p1k3 + 2/(1 - xija + zi)))/(p1k3*xija)

      endif

!---------- D^{a3}_2  --------   \tilde{q} -> \tilde{q}(p_a/p_1)+g(k3) ----------------------------------------------------------

      xija = (-k2k3 + p1k2 + p1k3)/(p1k2 + p1k3)
      if (xija.gt.x0FEISQ) then

      zi = p1k3/(p1k2 + p1k3)

      udip = M2 - 2*p2k1 !Q^2=udip=(p2-k1)^2, Q=p2-k1
      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p2k1)*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))
      tdip = M2 - 2*((p1k1*xija)/Rdip + p2k1*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))

      treeT1T4 = stsT2QQbar_M2TreeDipole(sdip,tdip,udip)/(-Nc)
      colorTree = treeT1T4/CF   ! T1 T4/CF

      stsT2xx_M2DipoleRealQ = stsT2xx_M2DipoleRealQ - (CF*colorTree*gs**2*(-2 - (M2*xija)/p1k3 + 2/(1 - xija + zi)))/(p1k3*xija)

      endif

!---------- D^{b3}_1  --------    \tilde{q} -> \tilde{q}(p_b/p_2)+g(k3) ----------------------------------------------------------

      xija = (-k1k3 + p2k1 + p2k3)/(p2k1 + p2k3)
      if (xija.gt.x0FEISQ) then

      zi =  p2k3/(p2k1 + p2k3)

      udip = M2 - 2*p1k2 !Q^2=udip=(p1-k2)^2, Q=p1-k2
      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p1k2)*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))
      tdip = M2 - 2*((p2k2*xija)/Rdip + p1k2*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))

      treeT2T3 = stsT2QQbar_M2TreeDipole(sdip,tdip,udip)/(-Nc)
      colorTree = treeT2T3/CF ! T2 T3/CF

      stsT2xx_M2DipoleRealQ = stsT2xx_M2DipoleRealQ - (CF*colorTree*gs**2*(-2 - (M2*xija)/p2k3 + 2/(1 - xija + zi)))/(p2k3*xija)

      endif

!---------- D^{b3}_2  --------    \tilde{q} -> \tilde{q}(p_b/p_2)+g(k3) ----------------------------------------------------------

      xija = (-k2k3 + p2k2 + p2k3)/(p2k2 + p2k3)
      if (xija.gt.x0FEISQ) then

      zi =  p2k3/(p2k2 + p2k3)

      tdip = M2 - 2*p1k1 !Q^2=tdip=(p1-k1)^2, Q=p1-k1
      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p1k1)*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))
      udip = M2 - 2*((p2k1*xija)/Rdip + p1k1*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))

      treeT2T4 = (2-Nc**2)/(2*Nc)*stsT2QQbar_M2TreeDipole(sdip,tdip,udip)
      colorTree = treeT2T4/CF ! T2 T4/CF = T1 T3/CF

      stsT2xx_M2DipoleRealQ = stsT2xx_M2DipoleRealQ - (CF*colorTree*gs**2*(-2 - (M2*xija)/p2k3 + 2/(1 - xija + zi)))/(p2k3*xija)

      endif



! *******************************************************************************************************************************************
! **************************************************** FE-IS ********************************************************************************
! *******************************************************************************************************************************************
!---------- D_{12}^a  --------  g -> q(k1)+qbar(k2) ----------------------------------------------------------
      xija = (-k1k2 + p1k1 + p1k2)/(p1k1 + p1k2)
      
      if (xija.gt.x0FEISg) then

      tdip = M2 - 2*p2k3 !Q^2=tdip=(p2-k3)^2, Q=p2-k3
      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p2k3)*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip))) ! sdip=(\tilde{p}_1+p2)^2
      udip = M2 - 2*((p1k3*xija)/Rdip + p2k3*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))  ! udip=(\tilde{p}_1-k3)^2

      zi =  p1k1/(p1k1 + p1k2)

      colorTree = stsT2gg_M2TreeDipoleT1T3(sdip,tdip,udip)/CV  ! T1 T3/CA

! diagonal part
      stsT2xx_M2DipoleRealQ = stsT2xx_M2DipoleRealQ - (colorTree*gs**2*TF)/(k1k2*xija)

! spin correlation part
      stsT2xx_M2DipoleRealQ = stsT2xx_M2DipoleRealQ + (2*stsT2xx_M2Dipole_SpinTensor_D12_a(sdip,tdip,udip)*gs**2*TF)/(k1k2**2*xija)

      endif


!---------- D_{12}^b  -------- g -> q(k1)+qbar(k2) ----------------------------------------------------------
      xija = (-k1k2 + p2k1 + p2k2)/(p2k1 + p2k2)

      if (xija.gt.x0FEISg) then

      udip = M2 - 2*p1k3 !Q^2=udip=(p1-k3)^2, Q=p1-k3
      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p1k3)*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip))) ! sdip=(\tilde{p}_2+p1)^2
      tdip = M2 - 2*((p2k3*xija)/Rdip + p1k3*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))  ! tdip=(\tilde{p}_2-k3)^2

      zi =  p2k1/(p2k1 + p2k2)

      colorTree = stsT2gg_M2TreeDipoleT2T3(sdip,tdip,udip)/CV   ! T2 T3/CA

! diagonal part
      stsT2xx_M2DipoleRealQ = stsT2xx_M2DipoleRealQ - (colorTree*gs**2*TF)/(k1k2*xija)

! spin correlation part
      stsT2xx_M2DipoleRealQ = stsT2xx_M2DipoleRealQ + (2*stsT2xx_M2Dipole_SpinTensor_D12_b(sdip,tdip,udip)*gs**2*TF)/(k1k2**2*xija)
      endif

!---------- D_{31}^a  -------- q -> q(k1)+g(k3) ----------------------------------------------------------

      xija = (-k1k3 + p1k1 + p1k3)/(p1k1 + p1k3)

      if (xija.gt.x0FEISQ) then

      tdip = M2 - 2*p2k2 !Q^2=tdip=(p2-k2)^2, Q=p2-k2
      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p2k2)*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip))) ! sdip=(\tilde{p}_1+p2)^2
      udip = M2 - 2*((p1k2*xija)/Rdip + p2k2*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip))) ! udip=(\tilde{p}_1-k2)^2

      zi = p1k3/(p1k1 + p1k3)
      zj = 1-zi
            
      treeT1T3 = (2-Nc**2)/(2*Nc)*stsT2QQbar_M2TreeDipole(sdip,tdip,udip)
      colorTree = treeT1T3/CF ! T1 T3/CF

      stsT2xx_M2DipoleRealQ = stsT2xx_M2DipoleRealQ -((CF*colorTree*gs**2*(-2 + zi + 2/(2 - xija - zj)))/(k1k3*xija))

      endif

!---------- D_{32}^a  -------- qbar -> qbar(k2)+g(k3) ----------------------------------------------------------

      xija = (-k2k3 + p1k2 + p1k3)/(p1k2 + p1k3)

      if (xija.gt.x0FEISQ) then

      udip = M2 - 2*p2k1 !Q^2=udip=(p2-k1)^2, Q=p2-k1
      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p2k1)*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip))) ! sdip=(\tilde{p}_1+p2)^2
      tdip = M2 - 2*((p1k1*xija)/Rdip + p2k1*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip))) ! tdip=(\tilde{p}_1-k2)^2
      
      zi = p1k3/(p1k2 + p1k3)
      zj = 1-zi
            
      treeT1T4 = stsT2QQbar_M2TreeDipole(sdip,tdip,udip)/(-Nc)
      colorTree = treeT1T4/CF  ! T1 T4/CF 

      stsT2xx_M2DipoleRealQ = stsT2xx_M2DipoleRealQ -((CF*colorTree*gs**2*(-2 + zi + 2/(2 - xija - zj)))/(k2k3*xija))

      endif


!---------- D_{31}^b  -------- q -> q(k1)+g(k3) ----------------------------------------------------------

      xija = (-k1k3 + p2k1 + p2k3)/(p2k1 + p2k3)

      if (xija.gt.x0FEISQ) then

      udip = M2 - 2*p1k2 !Q^2=udip=(p1-k2)^2, Q=p1-k2
      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p1k2)*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip))) ! sdip=(\tilde{p}_2+p1)^2
      tdip = M2 - 2*((p2k2*xija)/Rdip + p1k2*(-(M2 + udip)/(2*udip) + (-M2 + udip + 2*M2*xija)/(2*Rdip*udip)))  ! tdip=(\tilde{p}_1-k2)^2
      
      zi = p2k3/(p2k1 + p2k3)
      zj = 1-zi
            
      treeT2T3 = stsT2QQbar_M2TreeDipole(sdip,tdip,udip)/(-Nc)
      colorTree = treeT2T3/CF ! T2 T3/CF

      stsT2xx_M2DipoleRealQ = stsT2xx_M2DipoleRealQ -((CF*colorTree*gs**2*(-2 + zi + 2/(2 - xija - zj)))/(k1k3*xija))

      endif


!---------- D_{32}^b  -------- qbar -> qbar(k2)+g(k3) ----------------------------------------------------------

      xija = (-k2k3 + p2k2 + p2k3)/(p2k2 + p2k3)

      if (xija.gt.x0FEISQ) then

      tdip = M2 - 2*p1k1 !Q^2=tdip=(p1-k1)^2, Q=p1-k1
      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)
      sdip = 2*M2 + 2*((p1p2*xija)/Rdip + (M2 - p1k1)*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip))) ! sdip=(\tilde{p}_2+p1)^2
      udip = M2 - 2*((p2k1*xija)/Rdip + p1k1*(-(M2 + tdip)/(2*tdip) + (-M2 + tdip + 2*M2*xija)/(2*Rdip*tdip)))  ! udip=(\tilde{p}_1-k1)^2
      
      zi = p2k3/(p2k2 + p2k3)
      zj = 1-zi
            
      treeT1T3 = (2-Nc**2)/(2*Nc)*stsT2QQbar_M2TreeDipole(sdip,tdip,udip)
      colorTree = treeT1T3/CF   ! T2 T4/CF=T1 T3/CF

      stsT2xx_M2DipoleRealQ = stsT2xx_M2DipoleRealQ -((CF*colorTree*gs**2*(-2 + zi + 2/(2 - xija - zj)))/(k2k3*xija))

      endif
      
      end


      double precision function stsT2xx_M2Dipole_SpinTensor_D12_3(sdip,tdip,udip)
      implicit none
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"
#include "common/stsT2xx_GenCouplings.h"

      double precision sdip,tdip,udip,M2
      double precision yijk,zi,zj

      M2=MSf2(iSf1,itt1,igen1)

      yijk = k1k2/(k1k2 + k1k3 + k2k3)
      zi = k1k3/(k1k3 + k2k3)
      zj = 1-zi

!************************************************   Ms^ast VggP Ms   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_3 = ((GGG**2*Nc*(-1 + Nc**2)
     &   *((-4*k1k3**2*M2**2*zi**2)/(-1 + yijk)**2 + (16*k1k3**2*M2*sdip*zi**2)/(-1 + yijk)**2 - (16*k1k3**2
     &   *p1p2*sdip*zi**2)/(-1 + yijk)**2 + (4*k1k3**2*M2*tdip*zi**2)/(-1 + yijk)**2 + (4*k1k3**2*M2*udip*zi**2)/(
     &   -1 + yijk)**2 - (4*k1k3**2*tdip*udip*zi**2)/(-1 + yijk)**2 - (4*k1k3*M2*p1k1*sdip*zi**2)/(-1 + yijk) 
     &   + (4*k1k3*M2*p2k1*sdip*zi**2)/(-1 + yijk) + (8*k1k3*p1k1*sdip*tdip*zi**2)/(-1 + yijk) - (8*k1k3*p2k1
     &   *sdip*tdip*zi**2)/(-1 + yijk) - (4*k1k3*p1k1*sdip*udip*zi**2)/(-1 + yijk) + (4*k1k3*p2k1*sdip*udip
     &   *zi**2)/(-1 + yijk) + (4*k1k3*M2**2*(k1k2*(-1 + yijk) + k1k3*yijk)*zi**2)/(-1 + yijk)**2 - (16*k1k3*M2
     &   *sdip*(k1k2*(-1 + yijk) + k1k3*yijk)*zi**2)/(-1 + yijk)**2 + (16*k1k3*p1p2*sdip*(k1k2*(-1 + yijk) + k1k3
     &   *yijk)*zi**2)/(-1 + yijk)**2 - (2*k1k3*M2*tdip*(k1k2*(-1 + yijk) + k1k3*yijk)*zi**2)/(-1 + yijk)*
     &   *2 - (4*k1k3*tdip**2*(k1k2*(-1 + yijk) + k1k3*yijk)*zi**2)/(-1 + yijk)**2 - (6*k1k3*M2*udip*(k1k2*(-1 
     &   + yijk) + k1k3*yijk)*zi**2)/(-1 + yijk)**2 + (10*k1k3*tdip*udip*(k1k2*(-1 + yijk) + k1k3*yijk)*zi**2)/(
     &   -1 + yijk)**2 - (2*k1k3*udip**2*(k1k2*(-1 + yijk) + k1k3*yijk)*zi**2)/(-1 + yijk)**2 + (16*k1k3**2
     &   *M2*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi**2)/(-1 + yijk)**3 - (8*k1k3**2*tdip*((p2k1 + p2k2)*(-1 
     &   + yijk) + p2k3*yijk)*zi**2)/(-1 + yijk)**3 - (8*k1k3**2*udip*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)
     &   *zi**2)/(-1 + yijk)**3 + (8*k1k3*p1k1*sdip*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi**2)/(-1 + yijk)
     &   **2 - (8*k1k3*p2k1*sdip*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi**2)/(-1 + yijk)**2 - (16*k1k3*M2*(k1k2
     &   *(-1 + yijk) + k1k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi**2)/(-1 + yijk)**3 + (4*k1k3
     &   *tdip*(k1k2*(-1 + yijk) + k1k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi**2)/(-1 + yijk)**3 + (12
     &   *k1k3*udip*(k1k2*(-1 + yijk) + k1k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi**2)/(-1 + yijk)
     &   **3 - (4*p1k1*sdip*(k1k2*(-1 + yijk) + k1k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi**2)/(
     &   -1 + yijk)**2 + (4*p2k1*sdip*(k1k2*(-1 + yijk) + k1k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi
     &   **2)/(-1 + yijk)**2 + (8*tdip*(k1k2*(-1 + yijk) + k1k3*yijk)**2*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)
     &   *zi**2)/(-1 + yijk)**3 - (8*udip*(k1k2*(-1 + yijk) + k1k3*yijk)**2*((p2k1 + p2k2)*(-1 + yijk) + p2k3
     &   *yijk)*zi**2)/(-1 + yijk)**3 - (16*k1k3**2*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)**2*zi**2)/(-1 + yijk)
     &   **4 + (16*k1k3*(k1k2*(-1 + yijk) + k1k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)**2*zi**2)/(
     &   -1 + yijk)**4 + 2*M2*p1k1*sdip*(k1k2 + (k1k3*yijk)/(-1 + yijk))*zi**2 - 2*M2*p2k1*sdip*(k1k2 + (k1k3
     &   *yijk)/(-1 + yijk))*zi**2 - 4*p1k1*sdip*tdip*(k1k2 + (k1k3*yijk)/(-1 + yijk))*zi**2 + 4*p2k1*sdip*tdip
     &   *(k1k2 + (k1k3*yijk)/(-1 + yijk))*zi**2 + 2*p1k1*sdip*udip*(k1k2 + (k1k3*yijk)/(-1 + yijk))*zi**2 - 2
     &   *p2k1*sdip*udip*(k1k2 + (k1k3*yijk)/(-1 + yijk))*zi**2 + 4*M2*sdip*(k1k2 + (k1k3*yijk)/(-1 + yijk))*
     &   *2*zi**2 - 4*p1p2*sdip*(k1k2 + (k1k3*yijk)/(-1 + yijk))**2*zi**2 - 4*M2*tdip*(k1k2 + (k1k3*yijk)/(-1 
     &   + yijk))**2*zi**2 + 6*tdip**2*(k1k2 + (k1k3*yijk)/(-1 + yijk))**2*zi**2 + 4*M2*udip*(k1k2 + (k1k3*yijk)/(
     &   -1 + yijk))**2*zi**2 - 8*tdip*udip*(k1k2 + (k1k3*yijk)/(-1 + yijk))**2*zi**2 + 2*udip**2*(k1k2 + (k1k3
     &   *yijk)/(-1 + yijk))**2*zi**2 - k1k2*M2**2*sdip*zi*zj + 4*k1k2*M2*sdip*tdip*zi*zj - 4*k1k2*sdip*tdip
     &   **2*zi*zj - 2*k1k2*M2*sdip*udip*zi*zj + 4*k1k2*sdip*tdip*udip*zi*zj - k1k2*sdip*udip**2*zi*zj + (8
     &   *k1k3*k2k3*M2**2*zi*zj)/(-1 + yijk)**2 - (32*k1k3*k2k3*M2*sdip*zi*zj)/(-1 + yijk)**2 + (32*k1k3*k2k3*p1p2
     &   *sdip*zi*zj)/(-1 + yijk)**2 - (8*k1k3*k2k3*M2*tdip*zi*zj)/(-1 + yijk)**2 - (8*k1k3*k2k3*M2*udip*zi
     &   *zj)/(-1 + yijk)**2 + (8*k1k3*k2k3*tdip*udip*zi*zj)/(-1 + yijk)**2 + (4*k2k3*M2*p1k1*sdip*zi*zj)/(-1 
     &   + yijk) + (4*k1k3*M2*p1k2*sdip*zi*zj)/(-1 + yijk) - (4*k2k3*M2*p2k1*sdip*zi*zj)/(-1 + yijk) - (4*k1k3
     &   *M2*p2k2*sdip*zi*zj)/(-1 + yijk) - (8*k2k3*p1k1*sdip*tdip*zi*zj)/(-1 + yijk) - (8*k1k3*p1k2*sdip*tdip
     &   *zi*zj)/(-1 + yijk) + (8*k2k3*p2k1*sdip*tdip*zi*zj)/(-1 + yijk) + (8*k1k3*p2k2*sdip*tdip*zi*zj)/(-1 
     &   + yijk) + (4*k2k3*p1k1*sdip*udip*zi*zj)/(-1 + yijk) + (4*k1k3*p1k2*sdip*udip*zi*zj)/(-1 + yijk) - (4*k2k3
     &   *p2k1*sdip*udip*zi*zj)/(-1 + yijk) - (4*k1k3*p2k2*sdip*udip*zi*zj)/(-1 + yijk) - (4*k2k3*M2**2*(k1k2
     &   *(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 + yijk)**2 + (16*k2k3*M2*sdip*(k1k2*(-1 + yijk) + k1k3*yijk)*zi
     &   *zj)/(-1 + yijk)**2 - (16*k2k3*p1p2*sdip*(k1k2*(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 + yijk)**2 + (2*k2k3
     &   *M2*tdip*(k1k2*(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 + yijk)**2 + (4*k2k3*tdip**2*(k1k2*(-1 + yijk) 
     &   + k1k3*yijk)*zi*zj)/(-1 + yijk)**2 + (6*k2k3*M2*udip*(k1k2*(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 + yijk)
     &   **2 - (10*k2k3*tdip*udip*(k1k2*(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 + yijk)**2 + (2*k2k3*udip**2*(k1k2
     &   *(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 + yijk)**2 - (4*k1k3*M2**2*(k1k2*(-1 + yijk) + k2k3*yijk)*zi*zj)/(
     &   -1 + yijk)**2 + (16*k1k3*M2*sdip*(k1k2*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)**2 - (16*k1k3*p1p2
     &   *sdip*(k1k2*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)**2 + (2*k1k3*M2*tdip*(k1k2*(-1 + yijk) + k2k3
     &   *yijk)*zi*zj)/(-1 + yijk)**2 + (4*k1k3*tdip**2*(k1k2*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)*
     &   *2 + (6*k1k3*M2*udip*(k1k2*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)**2 - (10*k1k3*tdip*udip*(k1k2*(
     &   -1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)**2 + (2*k1k3*udip**2*(k1k2*(-1 + yijk) + k2k3*yijk)*zi*zj)/(
     &   -1 + yijk)**2 - (8*M2*sdip*(k1k2*(-1 + yijk) + k1k3*yijk)*(k1k2*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 
     &   + yijk)**2 + (8*p1p2*sdip*(k1k2*(-1 + yijk) + k1k3*yijk)*(k1k2*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 
     &   + yijk)**2 + (8*M2*tdip*(k1k2*(-1 + yijk) + k1k3*yijk)*(k1k2*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)
     &   **2 - (12*tdip**2*(k1k2*(-1 + yijk) + k1k3*yijk)*(k1k2*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)
     &   **2 - (8*M2*udip*(k1k2*(-1 + yijk) + k1k3*yijk)*(k1k2*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)**2 
     &   + (16*tdip*udip*(k1k2*(-1 + yijk) + k1k3*yijk)*(k1k2*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)**2 - (4
     &   *udip**2*(k1k2*(-1 + yijk) + k1k3*yijk)*(k1k2*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)**2 - (32
     &   *k1k3*k2k3*M2*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi*zj)/(-1 + yijk)**3 + (16*k1k3*k2k3*tdip*((p2k1 
     &   + p2k2)*(-1 + yijk) + p2k3*yijk)*zi*zj)/(-1 + yijk)**3 + (16*k1k3*k2k3*udip*((p2k1 + p2k2)*(-1 + yijk) 
     &   + p2k3*yijk)*zi*zj)/(-1 + yijk)**3 - (8*k2k3*p1k1*sdip*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi
     &   *zj)/(-1 + yijk)**2 - (8*k1k3*p1k2*sdip*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi*zj)/(-1 + yijk)**2 
     &   + (8*k2k3*p2k1*sdip*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi*zj)/(-1 + yijk)**2 + (8*k1k3*p2k2*sdip
     &   *((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi*zj)/(-1 + yijk)**2 + (16*k2k3*M2*(k1k2*(-1 + yijk) + k1k3
     &   *yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi*zj)/(-1 + yijk)**3 - (4*k2k3*tdip*(k1k2*(-1 + yijk) 
     &   + k1k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi*zj)/(-1 + yijk)**3 - (12*k2k3*udip*(k1k2*(
     &   -1 + yijk) + k1k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi*zj)/(-1 + yijk)**3 + (4*p1k2*sdip*(k1k2
     &   *(-1 + yijk) + k1k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi*zj)/(-1 + yijk)**2 - (4*p2k2
     &   *sdip*(k1k2*(-1 + yijk) + k1k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi*zj)/(-1 + yijk)**2 + (16
     &   *k1k3*M2*(k1k2*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi*zj)/(-1 + yijk)
     &   **3 - (4*k1k3*tdip*(k1k2*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi*zj)/(-1 
     &   + yijk)**3 - (12*k1k3*udip*(k1k2*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi
     &   *zj)/(-1 + yijk)**3 + (4*p1k1*sdip*(k1k2*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)
     &   *zi*zj)/(-1 + yijk)**2 - (4*p2k1*sdip*(k1k2*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3
     &   *yijk)*zi*zj)/(-1 + yijk)**2 - (16*tdip*(k1k2*(-1 + yijk) + k1k3*yijk)*(k1k2*(-1 + yijk) + k2k3*yijk)
     &   *((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi*zj)/(-1 + yijk)**3 + (16*udip*(k1k2*(-1 + yijk) + k1k3
     &   *yijk)*(k1k2*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zi*zj)/(-1 + yijk)**3 
     &   + (32*k1k3*k2k3*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)**2*zi*zj)/(-1 + yijk)**4 - (16*k2k3*(k1k2*(-1 
     &   + yijk) + k1k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)**2*zi*zj)/(-1 + yijk)**4 - (16*k1k3*(k1k2
     &   *(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)**2*zi*zj)/(-1 + yijk)**4 - 2*M2*p1k2
     &   *sdip*(k1k2 + (k1k3*yijk)/(-1 + yijk))*zi*zj + 2*M2*p2k2*sdip*(k1k2 + (k1k3*yijk)/(-1 + yijk))*zi
     &   *zj + 4*p1k2*sdip*tdip*(k1k2 + (k1k3*yijk)/(-1 + yijk))*zi*zj - 4*p2k2*sdip*tdip*(k1k2 + (k1k3*yijk)/(
     &   -1 + yijk))*zi*zj - 2*p1k2*sdip*udip*(k1k2 + (k1k3*yijk)/(-1 + yijk))*zi*zj + 2*p2k2*sdip*udip*(k1k2 
     &   + (k1k3*yijk)/(-1 + yijk))*zi*zj - 2*M2*p1k1*sdip*(k1k2 + (k2k3*yijk)/(-1 + yijk))*zi*zj + 2*M2*p2k1
     &   *sdip*(k1k2 + (k2k3*yijk)/(-1 + yijk))*zi*zj + 4*p1k1*sdip*tdip*(k1k2 + (k2k3*yijk)/(-1 + yijk))*zi*zj 
     &   - 4*p2k1*sdip*tdip*(k1k2 + (k2k3*yijk)/(-1 + yijk))*zi*zj - 2*p1k1*sdip*udip*(k1k2 + (k2k3*yijk)/(-1 
     &   + yijk))*zi*zj + 2*p2k1*sdip*udip*(k1k2 + (k2k3*yijk)/(-1 + yijk))*zi*zj + 4*k1k2*M2*sdip*(p2k1 + p2k2 
     &   + (p2k3*yijk)/(-1 + yijk))*zi*zj - 8*k1k2*sdip*tdip*(p2k1 + p2k2 + (p2k3*yijk)/(-1 + yijk))*zi*zj 
     &   + 4*k1k2*sdip*udip*(p2k1 + p2k2 + (p2k3*yijk)/(-1 + yijk))*zi*zj - 4*k1k2*sdip*(p2k1 + p2k2 + (p2k3*yijk)/(
     &   -1 + yijk))**2*zi*zj - (4*k2k3**2*M2**2*zj**2)/(-1 + yijk)**2 + (16*k2k3**2*M2*sdip*zj**2)/(-1 
     &   + yijk)**2 - (16*k2k3**2*p1p2*sdip*zj**2)/(-1 + yijk)**2 + (4*k2k3**2*M2*tdip*zj**2)/(-1 + yijk)**2 
     &   + (4*k2k3**2*M2*udip*zj**2)/(-1 + yijk)**2 - (4*k2k3**2*tdip*udip*zj**2)/(-1 + yijk)**2 - (4*k2k3*M2*p1k2
     &   *sdip*zj**2)/(-1 + yijk) + (4*k2k3*M2*p2k2*sdip*zj**2)/(-1 + yijk) + (8*k2k3*p1k2*sdip*tdip*zj**2)/(
     &   -1 + yijk) - (8*k2k3*p2k2*sdip*tdip*zj**2)/(-1 + yijk) - (4*k2k3*p1k2*sdip*udip*zj**2)/(-1 + yijk) 
     &   + (4*k2k3*p2k2*sdip*udip*zj**2)/(-1 + yijk) + (4*k2k3*M2**2*(k1k2*(-1 + yijk) + k2k3*yijk)*zj**2)/(-1 
     &   + yijk)**2 - (16*k2k3*M2*sdip*(k1k2*(-1 + yijk) + k2k3*yijk)*zj**2)/(-1 + yijk)**2 + (16*k2k3*p1p2*sdip
     &   *(k1k2*(-1 + yijk) + k2k3*yijk)*zj**2)/(-1 + yijk)**2 - (2*k2k3*M2*tdip*(k1k2*(-1 + yijk) + k2k3*yijk)
     &   *zj**2)/(-1 + yijk)**2 - (4*k2k3*tdip**2*(k1k2*(-1 + yijk) + k2k3*yijk)*zj**2)/(-1 + yijk)**2 - (6
     &   *k2k3*M2*udip*(k1k2*(-1 + yijk) + k2k3*yijk)*zj**2)/(-1 + yijk)**2 + (10*k2k3*tdip*udip*(k1k2*(-1 + yijk) 
     &   + k2k3*yijk)*zj**2)/(-1 + yijk)**2 - (2*k2k3*udip**2*(k1k2*(-1 + yijk) + k2k3*yijk)*zj**2)/(-1 
     &   + yijk)**2 + (16*k2k3**2*M2*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zj**2)/(-1 + yijk)**3 - (8*k2k3*
     &   *2*tdip*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zj**2)/(-1 + yijk)**3 - (8*k2k3**2*udip*((p2k1 + p2k2)
     &   *(-1 + yijk) + p2k3*yijk)*zj**2)/(-1 + yijk)**3 + (8*k2k3*p1k2*sdip*((p2k1 + p2k2)*(-1 + yijk) + p2k3
     &   *yijk)*zj**2)/(-1 + yijk)**2 - (8*k2k3*p2k2*sdip*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zj**2)/(-1 
     &   + yijk)**2 - (16*k2k3*M2*(k1k2*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)*zj**2)/(
     &   -1 + yijk)**3 + (4*k2k3*tdip*(k1k2*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3*yijk)
     &   *zj**2)/(-1 + yijk)**3 + (12*k2k3*udip*(k1k2*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) + p2k3
     &   *yijk)*zj**2)/(-1 + yijk)**3 - (4*p1k2*sdip*(k1k2*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k2)*(-1 + yijk) 
     &   + p2k3*yijk)*zj**2)/(-1 + yijk)**2 + (4*p2k2*sdip*(k1k2*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k2)*(
     &   -1 + yijk) + p2k3*yijk)*zj**2)/(-1 + yijk)**2 + (8*tdip*(k1k2*(-1 + yijk) + k2k3*yijk)**2*((p2k1 + p2k2)
     &   *(-1 + yijk) + p2k3*yijk)*zj**2)/(-1 + yijk)**3 - (8*udip*(k1k2*(-1 + yijk) + k2k3*yijk)**2*((p2k1 
     &   + p2k2)*(-1 + yijk) + p2k3*yijk)*zj**2)/(-1 + yijk)**3 - (16*k2k3**2*((p2k1 + p2k2)*(-1 + yijk) 
     &   + p2k3*yijk)**2*zj**2)/(-1 + yijk)**4 + (16*k2k3*(k1k2*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k2)*(-1 
     &   + yijk) + p2k3*yijk)**2*zj**2)/(-1 + yijk)**4 + 2*M2*p1k2*sdip*(k1k2 + (k2k3*yijk)/(-1 + yijk))*zj**2 
     &   - 2*M2*p2k2*sdip*(k1k2 + (k2k3*yijk)/(-1 + yijk))*zj**2 - 4*p1k2*sdip*tdip*(k1k2 + (k2k3*yijk)/(-1 
     &   + yijk))*zj**2 + 4*p2k2*sdip*tdip*(k1k2 + (k2k3*yijk)/(-1 + yijk))*zj**2 + 2*p1k2*sdip*udip*(k1k2 
     &   + (k2k3*yijk)/(-1 + yijk))*zj**2 - 2*p2k2*sdip*udip*(k1k2 + (k2k3*yijk)/(-1 + yijk))*zj**2 
     &   + 4*M2*sdip*(k1k2 + (k2k3*yijk)/(-1 + yijk))**2*zj**2 - 4*p1p2*sdip*(k1k2 + (k2k3*yijk)/(-1 
     &   + yijk))**2*zj**2 - 4*M2*tdip*(k1k2 + (k2k3*yijk)/(-1 + yijk))**2*zj**2 + 6*tdip**2*(k1k2 + (k2k3*yijk)/(-1 
     &   + yijk))**2*zj**2 + 4*M2*udip*(k1k2 + (k2k3*yijk)/(-1 + yijk))**2*zj**2 - 8*tdip*udip*(k1k2 + (k2k3*yijk)/(-1 
     &   + yijk))**2*zj**2 + 2*udip**2*(k1k2 + (k2k3*yijk)/(-1 + yijk))**2*zj**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(8*sdip**3))

!************************************************   2 Re(Ms VggP Mt^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_3 = stsT2xx_M2Dipole_SpinTensor_D12_3 + 2*dble(-0.125d0*(GGG*Nc
     &   *(-1 + Nc**2)*(-2*p1k1*(-1 + yijk)*zi + k1k3*yijk*zi + k1k2*(-1 + yijk)*(zi - zj) - 2*p1k2*zj - k2k3*yijk
     &   *zj + 2*p1k2*yijk*zj)*(M2*p2k1*sdip*zi - 2*p2k1**2*sdip*zi - 2*p2k1*p2k2*sdip*zi - 2*p2k1*sdip*tdip
     &   *zi + p2k1*sdip*udip*zi - 3*M2*p2k1*sdip*yijk*zi + 6*p2k1**2*sdip*yijk*zi + 6*p2k1*p2k2*sdip*yijk*zi 
     &   + 2*p2k1*p2k3*sdip*yijk*zi + 6*p2k1*sdip*tdip*yijk*zi - 3*p2k1*sdip*udip*yijk*zi + 3*M2*p2k1*sdip*yijk
     &   **2*zi - 6*p2k1**2*sdip*yijk**2*zi - 6*p2k1*p2k2*sdip*yijk**2*zi - 4*p2k1*p2k3*sdip*yijk**2*zi - 6
     &   *p2k1*sdip*tdip*yijk**2*zi + 3*p2k1*sdip*udip*yijk**2*zi - M2*p2k1*sdip*yijk**3*zi + 2*p2k1**2*sdip*yijk
     &   **3*zi + 2*p2k1*p2k2*sdip*yijk**3*zi + 2*p2k1*p2k3*sdip*yijk**3*zi + 2*p2k1*sdip*tdip*yijk**3*zi - p2k1
     &   *sdip*udip*yijk**3*zi + k1k3*(4*p2k2**2 + 4*p1p2*sdip - 4*p2k2*tdip - 2*tdip**2 + 2*p2k2*udip + 4
     &   *p2k1**2*(-1 + yijk)**2 + 2*M2**2*(-1 + yijk)**3 - 8*p2k2**2*yijk - 8*p2k2*p2k3*yijk - 10*p1p2*sdip*yijk 
     &   + 14*p2k2*tdip*yijk + 4*p2k3*tdip*yijk + 7*tdip**2*yijk - 6*p2k2*udip*yijk - 2*p2k3*udip*yijk - tdip
     &   *udip*yijk + 4*p2k2**2*yijk**2 + 8*p2k2*p2k3*yijk**2 + 4*p2k3**2*yijk**2 + 8*p1p2*sdip*yijk**2 - 16
     &   *p2k2*tdip*yijk**2 - 10*p2k3*tdip*yijk**2 - 8*tdip**2*yijk**2 + 6*p2k2*udip*yijk**2 + 4*p2k3*udip*yijk
     &   **2 + 2*tdip*udip*yijk**2 - 2*p1p2*sdip*yijk**3 + 6*p2k2*tdip*yijk**3 + 6*p2k3*tdip*yijk**3 + 3*tdip
     &   **2*yijk**3 - 2*p2k2*udip*yijk**3 - 2*p2k3*udip*yijk**3 - tdip*udip*yijk**3 + 2*p2k1*(-1 + yijk)*(-udip 
     &   + 4*p2k2*(-1 + yijk) + 4*p2k3*yijk + 2*udip*yijk - udip*yijk**2 + tdip*(2 - 5*yijk + 3*yijk**2)) 
     &   - M2*(-1 + yijk)*(-4*sdip + 4*tdip - 2*p2k3*yijk + 6*sdip*yijk - 9*tdip*yijk + udip*yijk + 4*p2k3*yijk
     &   **2 - 2*sdip*yijk**2 + 5*tdip*yijk**2 - udip*yijk**2 + p2k1*(2 - 6*yijk + 4*yijk**2) + p2k2*(2 - 6*yijk 
     &   + 4*yijk**2)))*zi + k1k2*(-1 + yijk)**2*(2*M2**2*(-1 + yijk) - 2*p1p2*sdip*(-1 + yijk) + (3*tdip 
     &   - udip)*(-tdip + 2*p2k1*(-1 + yijk) + 2*p2k2*(-1 + yijk) + 2*p2k3*yijk + tdip*yijk) - M2*(2*sdip - 5
     &   *tdip + udip + 4*p2k1*(-1 + yijk) + 4*p2k2*(-1 + yijk) + 4*p2k3*yijk - 2*sdip*yijk + 5*tdip*yijk - udip
     &   *yijk))*(zi - zj) + 2*k2k3*M2**2*zj - 2*k2k3*M2*p2k1*zj - 4*k2k3*p2k1**2*zj - 2*k2k3*M2*p2k2*zj - 8
     &   *k2k3*p2k1*p2k2*zj - 4*k2k3*p2k2**2*zj + 4*k2k3*M2*sdip*zj - 4*k2k3*p1p2*sdip*zj - M2*p2k2*sdip*zj + 2
     &   *p2k1*p2k2*sdip*zj + 2*p2k2**2*sdip*zj - 4*k2k3*M2*tdip*zj + 4*k2k3*p2k1*tdip*zj + 4*k2k3*p2k2*tdip*zj 
     &   + 2*p2k2*sdip*tdip*zj + 2*k2k3*tdip**2*zj - 2*k2k3*p2k1*udip*zj - 2*k2k3*p2k2*udip*zj - p2k2*sdip*udip
     &   *zj - 6*k2k3*M2**2*yijk*zj + 8*k2k3*M2*p2k1*yijk*zj + 8*k2k3*p2k1**2*yijk*zj + 8*k2k3*M2*p2k2*yijk
     &   *zj + 16*k2k3*p2k1*p2k2*yijk*zj + 8*k2k3*p2k2**2*yijk*zj + 2*k2k3*M2*p2k3*yijk*zj + 8*k2k3*p2k1*p2k3
     &   *yijk*zj + 8*k2k3*p2k2*p2k3*yijk*zj - 10*k2k3*M2*sdip*yijk*zj + 10*k2k3*p1p2*sdip*yijk*zj + 3*M2*p2k2
     &   *sdip*yijk*zj - 6*p2k1*p2k2*sdip*yijk*zj - 6*p2k2**2*sdip*yijk*zj - 2*p2k2*p2k3*sdip*yijk*zj + 13*k2k3
     &   *M2*tdip*yijk*zj - 14*k2k3*p2k1*tdip*yijk*zj - 14*k2k3*p2k2*tdip*yijk*zj - 4*k2k3*p2k3*tdip*yijk*zj 
     &   - 6*p2k2*sdip*tdip*yijk*zj - 7*k2k3*tdip**2*yijk*zj - k2k3*M2*udip*yijk*zj + 6*k2k3*p2k1*udip*yijk*zj 
     &   + 6*k2k3*p2k2*udip*yijk*zj + 2*k2k3*p2k3*udip*yijk*zj + 3*p2k2*sdip*udip*yijk*zj + k2k3*tdip*udip*yijk
     &   *zj + 6*k2k3*M2**2*yijk**2*zj - 10*k2k3*M2*p2k1*yijk**2*zj - 4*k2k3*p2k1**2*yijk**2*zj - 10*k2k3*M2
     &   *p2k2*yijk**2*zj - 8*k2k3*p2k1*p2k2*yijk**2*zj - 4*k2k3*p2k2**2*yijk**2*zj - 6*k2k3*M2*p2k3*yijk**2*zj 
     &   - 8*k2k3*p2k1*p2k3*yijk**2*zj - 8*k2k3*p2k2*p2k3*yijk**2*zj - 4*k2k3*p2k3**2*yijk**2*zj + 8*k2k3*M2
     &   *sdip*yijk**2*zj - 8*k2k3*p1p2*sdip*yijk**2*zj - 3*M2*p2k2*sdip*yijk**2*zj + 6*p2k1*p2k2*sdip*yijk**2
     &   *zj + 6*p2k2**2*sdip*yijk**2*zj + 4*p2k2*p2k3*sdip*yijk**2*zj - 14*k2k3*M2*tdip*yijk**2*zj + 16*k2k3*p2k1
     &   *tdip*yijk**2*zj + 16*k2k3*p2k2*tdip*yijk**2*zj + 10*k2k3*p2k3*tdip*yijk**2*zj + 6*p2k2*sdip*tdip
     &   *yijk**2*zj + 8*k2k3*tdip**2*yijk**2*zj + 2*k2k3*M2*udip*yijk**2*zj - 6*k2k3*p2k1*udip*yijk**2*zj - 6
     &   *k2k3*p2k2*udip*yijk**2*zj - 4*k2k3*p2k3*udip*yijk**2*zj - 3*p2k2*sdip*udip*yijk**2*zj - 2*k2k3*tdip*udip
     &   *yijk**2*zj - 2*k2k3*M2**2*yijk**3*zj + 4*k2k3*M2*p2k1*yijk**3*zj + 4*k2k3*M2*p2k2*yijk**3*zj + 4
     &   *k2k3*M2*p2k3*yijk**3*zj - 2*k2k3*M2*sdip*yijk**3*zj + 2*k2k3*p1p2*sdip*yijk**3*zj + M2*p2k2*sdip*yijk
     &   **3*zj - 2*p2k1*p2k2*sdip*yijk**3*zj - 2*p2k2**2*sdip*yijk**3*zj - 2*p2k2*p2k3*sdip*yijk**3*zj + 5*k2k3
     &   *M2*tdip*yijk**3*zj - 6*k2k3*p2k1*tdip*yijk**3*zj - 6*k2k3*p2k2*tdip*yijk**3*zj - 6*k2k3*p2k3*tdip*yijk**3*zj 
     &   - 2*p2k2*sdip*tdip*yijk**3*zj - 3*k2k3*tdip**2*yijk**3*zj - k2k3*M2*udip*yijk**3*zj + 2*k2k3*p2k1*udip*yijk**3*zj 
     &   + 2*k2k3*p2k2*udip*yijk**3*zj + 2*k2k3*p2k3*udip*yijk**3*zj + p2k2*sdip*udip*yijk**3*zj + k2k3*tdip
     &   *udip*yijk**3*zj)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(sdip**2*(M2 - tdip)*(-1 + yijk)**4))

!************************************************   2 Re(Ms VggP Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_3 = stsT2xx_M2Dipole_SpinTensor_D12_3 + 2*dble(-0.125d0*(GGG*Nc
     &   *(-1 + Nc**2)*(-2*p2k1*(-1 + yijk)*zi + k1k3*yijk*zi + k1k2*(-1 + yijk)*(zi - zj) - 2*p2k2*zj - k2k3*yijk
     &   *zj + 2*p2k2*yijk*zj)*(M2*p1k1*sdip*zi - 2*p1k1*p2k1*sdip*zi - 2*p1k1*p2k2*sdip*zi - 2*p1k1*sdip*tdip
     &   *zi + p1k1*sdip*udip*zi - 2*M2*p1k1*sdip*yijk*zi + 4*p1k1*p2k1*sdip*yijk*zi + 4*p1k1*p2k2*sdip*yijk
     &   *zi + 2*p1k1*p2k3*sdip*yijk*zi + 4*p1k1*sdip*tdip*yijk*zi - 2*p1k1*sdip*udip*yijk*zi + M2*p1k1*sdip
     &   *yijk**2*zi - 2*p1k1*p2k1*sdip*yijk**2*zi - 2*p1k1*p2k2*sdip*yijk**2*zi - 2*p1k1*p2k3*sdip*yijk**2*zi 
     &   - 2*p1k1*sdip*tdip*yijk**2*zi + p1k1*sdip*udip*yijk**2*zi + k1k3*(-2*p2k1*tdip - 2*p2k2*tdip - 4*p2k1
     &   *udip - 4*p2k2*udip - 3*tdip*udip + 2*p2k1*tdip*yijk + 2*p2k2*tdip*yijk + 2*p2k3*tdip*yijk - tdip**2
     &   *yijk + 8*p2k1*udip*yijk + 8*p2k2*udip*yijk + 4*p2k3*udip*yijk + 7*tdip*udip*yijk - udip**2*yijk + tdip
     &   **2*yijk**2 - 4*p2k1*udip*yijk**2 - 4*p2k2*udip*yijk**2 - 4*p2k3*udip*yijk**2 - 4*tdip*udip*yijk**2 
     &   + udip**2*yijk**2 + M2**2*(-3 + 5*yijk - 2*yijk**2) - 2*p1p2*sdip*(2 - 3*yijk + yijk**2) + M2*(4*sdip 
     &   + 3*tdip + 3*udip - 6*p2k3*yijk - 6*sdip*yijk - 5*tdip*yijk - 5*udip*yijk + 4*p2k3*yijk**2 + 2*sdip
     &   *yijk**2 + 2*tdip*yijk**2 + 2*udip*yijk**2 + 2*p2k1*(3 - 5*yijk + 2*yijk**2) + 2*p2k2*(3 - 5*yijk + 2
     &   *yijk**2)))*zi - k1k2*(-1 + yijk)*(tdip**2 - 4*p2k1*udip - 4*p2k2*udip - 4*tdip*udip + udip**2 + 2*M2
     &   **2*(-1 + yijk) + 2*p1p2*sdip*(-1 + yijk) - tdip**2*yijk + 4*p2k1*udip*yijk + 4*p2k2*udip*yijk + 4*p2k3
     &   *udip*yijk + 4*tdip*udip*yijk - udip**2*yijk - 2*M2*(-sdip - tdip - udip + 2*p2k1*(-1 + yijk) + 2*p2k2
     &   *(-1 + yijk) + 2*p2k3*yijk + sdip*yijk + tdip*yijk + udip*yijk))*(zi - zj) + 3*k2k3*M2**2*zj - 6*k2k3
     &   *M2*p2k1*zj - 6*k2k3*M2*p2k2*zj - 4*k2k3*M2*sdip*zj - M2*p1k2*sdip*zj + 4*k2k3*p1p2*sdip*zj + 2*p1k2
     &   *p2k1*sdip*zj + 2*p1k2*p2k2*sdip*zj - 3*k2k3*M2*tdip*zj + 2*k2k3*p2k1*tdip*zj + 2*k2k3*p2k2*tdip*zj 
     &   + 2*p1k2*sdip*tdip*zj - 3*k2k3*M2*udip*zj + 4*k2k3*p2k1*udip*zj + 4*k2k3*p2k2*udip*zj - p1k2*sdip*udip
     &   *zj + 3*k2k3*tdip*udip*zj - 5*k2k3*M2**2*yijk*zj + 10*k2k3*M2*p2k1*yijk*zj + 10*k2k3*M2*p2k2*yijk*zj 
     &   + 6*k2k3*M2*p2k3*yijk*zj + 6*k2k3*M2*sdip*yijk*zj + 2*M2*p1k2*sdip*yijk*zj - 6*k2k3*p1p2*sdip*yijk*zj 
     &   - 4*p1k2*p2k1*sdip*yijk*zj - 4*p1k2*p2k2*sdip*yijk*zj - 2*p1k2*p2k3*sdip*yijk*zj + 5*k2k3*M2*tdip*yijk
     &   *zj - 2*k2k3*p2k1*tdip*yijk*zj - 2*k2k3*p2k2*tdip*yijk*zj - 2*k2k3*p2k3*tdip*yijk*zj - 4*p1k2*sdip
     &   *tdip*yijk*zj + k2k3*tdip**2*yijk*zj + 5*k2k3*M2*udip*yijk*zj - 8*k2k3*p2k1*udip*yijk*zj - 8*k2k3*p2k2
     &   *udip*yijk*zj - 4*k2k3*p2k3*udip*yijk*zj + 2*p1k2*sdip*udip*yijk*zj - 7*k2k3*tdip*udip*yijk*zj + k2k3
     &   *udip**2*yijk*zj + 2*k2k3*M2**2*yijk**2*zj - 4*k2k3*M2*p2k1*yijk**2*zj - 4*k2k3*M2*p2k2*yijk**2*zj 
     &   - 4*k2k3*M2*p2k3*yijk**2*zj - 2*k2k3*M2*sdip*yijk**2*zj - M2*p1k2*sdip*yijk**2*zj + 2*k2k3*p1p2*sdip*yijk
     &   **2*zj + 2*p1k2*p2k1*sdip*yijk**2*zj + 2*p1k2*p2k2*sdip*yijk**2*zj + 2*p1k2*p2k3*sdip*yijk**2*zj 
     &   - 2*k2k3*M2*tdip*yijk**2*zj + 2*p1k2*sdip*tdip*yijk**2*zj - k2k3*tdip**2*yijk**2*zj - 2*k2k3*M2*udip*yijk**2*zj 
     &  + 4*k2k3*p2k1*udip*yijk**2*zj + 4*k2k3*p2k2*udip*yijk**2*zj + 4*k2k3*p2k3*udip*yijk**2*zj - p1k2*sdip*udip*yijk**2*zj 
     &  + 4*k2k3*tdip*udip*yijk**2*zj - k2k3*udip**2*yijk**2*zj)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(sdip**2*(M2 - udip)*(-1 + yijk)**3))

!************************************************   2 Re(Ms VggP Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_3 = stsT2xx_M2Dipole_SpinTensor_D12_3 + 2*dble(0)

!************************************************   Mt^ast VggP Mt   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_3 = stsT2xx_M2Dipole_SpinTensor_D12_3 +((Nc*(-1 + Nc**2)*(M2*sdip
     &   *(-1 + yijk) - 2*(p2k1 + p2k2)*(M2 - tdip)*(-1 + yijk) + 2*p2k3*(-M2 + tdip)*yijk)*(-2*p1k1*(-1 + yijk)
     &   *zi + k1k3*yijk*zi + k1k2*(-1 + yijk)*(zi - zj) - 2*p1k2*zj - k2k3*yijk*zj + 2*p1k2*yijk*zj)**2
     &   *SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(2*sdip*(M2 - tdip)**2*(-1 + yijk)**3))

!************************************************   2 Re(Mt VggP Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_3 = stsT2xx_M2Dipole_SpinTensor_D12_3 + 2*dble(0)

!************************************************   2 Re(Mt VggP Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_3 = stsT2xx_M2Dipole_SpinTensor_D12_3 + 2*dble((Nc*(-1 + Nc**2)
     &   *(-2*p1k1*(-1 + yijk)*zi + k1k3*yijk*zi + k1k2*(-1 + yijk)*(zi - zj) - 2*p1k2*zj - k2k3*yijk*zj + 2
     &   *p1k2*yijk*zj)*(p2k1*sdip*zi - 2*p2k1*sdip*yijk*zi + p2k1*sdip*yijk**2*zi + k1k3*(2*p2k1*(-1 + yijk) 
     &   + 2*p2k2*(-1 + yijk) + yijk*(M2 + 2*p2k3 + tdip*(-1 + yijk) - M2*yijk))*zi - k1k2*(M2 - tdip)*(-1 + yijk)
     &   **2*(zi - zj) + 2*k2k3*p2k1*zj + 2*k2k3*p2k2*zj - p2k2*sdip*zj - k2k3*M2*yijk*zj - 2*k2k3*p2k1*yijk
     &   *zj - 2*k2k3*p2k2*yijk*zj - 2*k2k3*p2k3*yijk*zj + 2*p2k2*sdip*yijk*zj + k2k3*tdip*yijk*zj + k2k3*M2
     &   *yijk**2*zj - p2k2*sdip*yijk**2*zj - k2k3*tdip*yijk**2*zj)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(4.
     &   *sdip*(M2 - tdip)*(-1 + yijk)**3))

!************************************************   Mu^ast VggP Mu   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_3 = stsT2xx_M2Dipole_SpinTensor_D12_3 +((Nc*(-1 + Nc**2)*(-M2
     &   **2 - tdip*udip + M2*(sdip + tdip + udip))*(-2*p2k1*(-1 + yijk)*zi + k1k3*yijk*zi + k1k2*(-1 + yijk)*(zi 
     &   - zj) - 2*p2k2*zj - k2k3*yijk*zj + 2*p2k2*yijk*zj)**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(2
     &   *sdip*(M2 - udip)**2*(-1 + yijk)**2))

!************************************************   2 Re(Mu VggP Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_3 = stsT2xx_M2Dipole_SpinTensor_D12_3 + 2*dble(-0.25d0*(Nc*(-1 
     &   + Nc**2)*(-2*p2k1*(-1 + yijk)*zi + k1k3*yijk*zi + k1k2*(-1 + yijk)*(zi - zj) - 2*p2k2*zj - k2k3*yijk*zj 
     &   + 2*p2k2*yijk*zj)*(p1k1*sdip*zi - p1k1*sdip*yijk*zi + k1k3*(tdip + M2*(-1 + yijk) - udip*yijk)*zi 
     &   + k1k2*(M2 - udip)*(-1 + yijk)*(zi - zj) + k2k3*M2*zj - p1k2*sdip*zj - k2k3*tdip*zj - k2k3*M2*yijk*zj 
     &   + p1k2*sdip*yijk*zj + k2k3*udip*yijk*zj)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(sdip
     &   *(M2 - udip)*(-1 + yijk)**2))

!************************************************   Mv^ast VggP Mv   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_3 = stsT2xx_M2Dipole_SpinTensor_D12_3 +((Nc*(-1 + Nc**2)*(2*yijk
     &   *(k1k3*zi - k2k3*zj)**2 + k1k2*(-1 + yijk)*(2*k1k3*zi*(zi - zj) + zj*(-2*k2k3*zi + sdip*zi - sdip*yijk
     &   *zi + 2*k2k3*zj)))*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(2*sdip*(-1 + yijk)**2))


      end



      double precision function stsT2xx_M2Dipole_SpinTensor_D13_2(sdip,tdip,udip)
      implicit none
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"
#include "common/stsT2xx_GenCouplings.h"

      double precision sdip,tdip,udip,M2
      double precision yijk,zi,zj

      M2=MSf2(iSf1,itt1,igen1)

      yijk = k1k3/(k1k2 + k1k3 + k2k3)
      zi = k1k2/(k1k2 + k2k3)
      zj = 1-zi

!************************************************   Ms^ast VggP Ms   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_2 = ((GGG**2*Nc*(-1 + Nc**2)
     &   *((-4*k1k2**2*M2**2*zi**2)/(-1 + yijk)**2 + (16*k1k2**2*M2*sdip*zi**2)/(-1 + yijk)**2 - (16*k1k2**2
     &   *p1p2*sdip*zi**2)/(-1 + yijk)**2 + (4*k1k2**2*M2*tdip*zi**2)/(-1 + yijk)**2 + (4*k1k2**2*M2*udip*zi**2)/(
     &   -1 + yijk)**2 - (4*k1k2**2*tdip*udip*zi**2)/(-1 + yijk)**2 - (4*k1k2*M2*p1k1*sdip*zi**2)/(-1 + yijk) 
     &   + (4*k1k2*M2*p2k1*sdip*zi**2)/(-1 + yijk) + (8*k1k2*p1k1*sdip*tdip*zi**2)/(-1 + yijk) - (8*k1k2*p2k1
     &   *sdip*tdip*zi**2)/(-1 + yijk) - (4*k1k2*p1k1*sdip*udip*zi**2)/(-1 + yijk) + (4*k1k2*p2k1*sdip*udip
     &   *zi**2)/(-1 + yijk) + (4*k1k2*M2**2*(k1k3*(-1 + yijk) + k1k2*yijk)*zi**2)/(-1 + yijk)**2 - (16*k1k2*M2
     &   *sdip*(k1k3*(-1 + yijk) + k1k2*yijk)*zi**2)/(-1 + yijk)**2 + (16*k1k2*p1p2*sdip*(k1k3*(-1 + yijk) + k1k2
     &   *yijk)*zi**2)/(-1 + yijk)**2 - (2*k1k2*M2*tdip*(k1k3*(-1 + yijk) + k1k2*yijk)*zi**2)/(-1 + yijk)*
     &   *2 - (4*k1k2*tdip**2*(k1k3*(-1 + yijk) + k1k2*yijk)*zi**2)/(-1 + yijk)**2 - (6*k1k2*M2*udip*(k1k3*(-1 
     &   + yijk) + k1k2*yijk)*zi**2)/(-1 + yijk)**2 + (10*k1k2*tdip*udip*(k1k3*(-1 + yijk) + k1k2*yijk)*zi**2)/(
     &   -1 + yijk)**2 - (2*k1k2*udip**2*(k1k3*(-1 + yijk) + k1k2*yijk)*zi**2)/(-1 + yijk)**2 + (16*k1k2**2
     &   *M2*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi**2)/(-1 + yijk)**3 - (8*k1k2**2*tdip*((p2k1 + p2k3)*(-1 
     &   + yijk) + p2k2*yijk)*zi**2)/(-1 + yijk)**3 - (8*k1k2**2*udip*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)
     &   *zi**2)/(-1 + yijk)**3 + (8*k1k2*p1k1*sdip*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi**2)/(-1 + yijk)
     &   **2 - (8*k1k2*p2k1*sdip*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi**2)/(-1 + yijk)**2 - (16*k1k2*M2*(k1k3
     &   *(-1 + yijk) + k1k2*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi**2)/(-1 + yijk)**3 + (4*k1k2
     &   *tdip*(k1k3*(-1 + yijk) + k1k2*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi**2)/(-1 + yijk)**3 + (12
     &   *k1k2*udip*(k1k3*(-1 + yijk) + k1k2*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi**2)/(-1 + yijk)
     &   **3 - (4*p1k1*sdip*(k1k3*(-1 + yijk) + k1k2*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi**2)/(
     &   -1 + yijk)**2 + (4*p2k1*sdip*(k1k3*(-1 + yijk) + k1k2*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi
     &   **2)/(-1 + yijk)**2 + (8*tdip*(k1k3*(-1 + yijk) + k1k2*yijk)**2*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)
     &   *zi**2)/(-1 + yijk)**3 - (8*udip*(k1k3*(-1 + yijk) + k1k2*yijk)**2*((p2k1 + p2k3)*(-1 + yijk) + p2k2
     &   *yijk)*zi**2)/(-1 + yijk)**3 - (16*k1k2**2*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)**2*zi**2)/(-1 + yijk)
     &   **4 + (16*k1k2*(k1k3*(-1 + yijk) + k1k2*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)**2*zi**2)/(
     &   -1 + yijk)**4 + 2*M2*p1k1*sdip*(k1k3 + (k1k2*yijk)/(-1 + yijk))*zi**2 - 2*M2*p2k1*sdip*(k1k3 + (k1k2
     &   *yijk)/(-1 + yijk))*zi**2 - 4*p1k1*sdip*tdip*(k1k3 + (k1k2*yijk)/(-1 + yijk))*zi**2 + 4*p2k1*sdip*tdip
     &   *(k1k3 + (k1k2*yijk)/(-1 + yijk))*zi**2 + 2*p1k1*sdip*udip*(k1k3 + (k1k2*yijk)/(-1 + yijk))*zi**2 - 2
     &   *p2k1*sdip*udip*(k1k3 + (k1k2*yijk)/(-1 + yijk))*zi**2 + 4*M2*sdip*(k1k3 + (k1k2*yijk)/(-1 + yijk))*
     &   *2*zi**2 - 4*p1p2*sdip*(k1k3 + (k1k2*yijk)/(-1 + yijk))**2*zi**2 - 4*M2*tdip*(k1k3 + (k1k2*yijk)/(-1 
     &   + yijk))**2*zi**2 + 6*tdip**2*(k1k3 + (k1k2*yijk)/(-1 + yijk))**2*zi**2 + 4*M2*udip*(k1k3 + (k1k2*yijk)/(
     &   -1 + yijk))**2*zi**2 - 8*tdip*udip*(k1k3 + (k1k2*yijk)/(-1 + yijk))**2*zi**2 + 2*udip**2*(k1k3 + (k1k2
     &   *yijk)/(-1 + yijk))**2*zi**2 - k1k3*M2**2*sdip*zi*zj + 4*k1k3*M2*sdip*tdip*zi*zj - 4*k1k3*sdip*tdip
     &   **2*zi*zj - 2*k1k3*M2*sdip*udip*zi*zj + 4*k1k3*sdip*tdip*udip*zi*zj - k1k3*sdip*udip**2*zi*zj + (8
     &   *k1k2*k2k3*M2**2*zi*zj)/(-1 + yijk)**2 - (32*k1k2*k2k3*M2*sdip*zi*zj)/(-1 + yijk)**2 + (32*k1k2*k2k3*p1p2
     &   *sdip*zi*zj)/(-1 + yijk)**2 - (8*k1k2*k2k3*M2*tdip*zi*zj)/(-1 + yijk)**2 - (8*k1k2*k2k3*M2*udip*zi
     &   *zj)/(-1 + yijk)**2 + (8*k1k2*k2k3*tdip*udip*zi*zj)/(-1 + yijk)**2 + (4*k2k3*M2*p1k1*sdip*zi*zj)/(-1 
     &   + yijk) + (4*k1k2*M2*p1k3*sdip*zi*zj)/(-1 + yijk) - (4*k2k3*M2*p2k1*sdip*zi*zj)/(-1 + yijk) - (4*k1k2
     &   *M2*p2k3*sdip*zi*zj)/(-1 + yijk) - (8*k2k3*p1k1*sdip*tdip*zi*zj)/(-1 + yijk) - (8*k1k2*p1k3*sdip*tdip
     &   *zi*zj)/(-1 + yijk) + (8*k2k3*p2k1*sdip*tdip*zi*zj)/(-1 + yijk) + (8*k1k2*p2k3*sdip*tdip*zi*zj)/(-1 
     &   + yijk) + (4*k2k3*p1k1*sdip*udip*zi*zj)/(-1 + yijk) + (4*k1k2*p1k3*sdip*udip*zi*zj)/(-1 + yijk) - (4*k2k3
     &   *p2k1*sdip*udip*zi*zj)/(-1 + yijk) - (4*k1k2*p2k3*sdip*udip*zi*zj)/(-1 + yijk) - (4*k2k3*M2**2*(k1k3
     &   *(-1 + yijk) + k1k2*yijk)*zi*zj)/(-1 + yijk)**2 + (16*k2k3*M2*sdip*(k1k3*(-1 + yijk) + k1k2*yijk)*zi
     &   *zj)/(-1 + yijk)**2 - (16*k2k3*p1p2*sdip*(k1k3*(-1 + yijk) + k1k2*yijk)*zi*zj)/(-1 + yijk)**2 + (2*k2k3
     &   *M2*tdip*(k1k3*(-1 + yijk) + k1k2*yijk)*zi*zj)/(-1 + yijk)**2 + (4*k2k3*tdip**2*(k1k3*(-1 + yijk) 
     &   + k1k2*yijk)*zi*zj)/(-1 + yijk)**2 + (6*k2k3*M2*udip*(k1k3*(-1 + yijk) + k1k2*yijk)*zi*zj)/(-1 + yijk)
     &   **2 - (10*k2k3*tdip*udip*(k1k3*(-1 + yijk) + k1k2*yijk)*zi*zj)/(-1 + yijk)**2 + (2*k2k3*udip**2*(k1k3
     &   *(-1 + yijk) + k1k2*yijk)*zi*zj)/(-1 + yijk)**2 - (4*k1k2*M2**2*(k1k3*(-1 + yijk) + k2k3*yijk)*zi*zj)/(
     &   -1 + yijk)**2 + (16*k1k2*M2*sdip*(k1k3*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)**2 - (16*k1k2*p1p2
     &   *sdip*(k1k3*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)**2 + (2*k1k2*M2*tdip*(k1k3*(-1 + yijk) + k2k3
     &   *yijk)*zi*zj)/(-1 + yijk)**2 + (4*k1k2*tdip**2*(k1k3*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)*
     &   *2 + (6*k1k2*M2*udip*(k1k3*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)**2 - (10*k1k2*tdip*udip*(k1k3*(
     &   -1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)**2 + (2*k1k2*udip**2*(k1k3*(-1 + yijk) + k2k3*yijk)*zi*zj)/(
     &   -1 + yijk)**2 - (8*M2*sdip*(k1k3*(-1 + yijk) + k1k2*yijk)*(k1k3*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 
     &   + yijk)**2 + (8*p1p2*sdip*(k1k3*(-1 + yijk) + k1k2*yijk)*(k1k3*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 
     &   + yijk)**2 + (8*M2*tdip*(k1k3*(-1 + yijk) + k1k2*yijk)*(k1k3*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)
     &   **2 - (12*tdip**2*(k1k3*(-1 + yijk) + k1k2*yijk)*(k1k3*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)
     &   **2 - (8*M2*udip*(k1k3*(-1 + yijk) + k1k2*yijk)*(k1k3*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)**2 
     &   + (16*tdip*udip*(k1k3*(-1 + yijk) + k1k2*yijk)*(k1k3*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)**2 - (4
     &   *udip**2*(k1k3*(-1 + yijk) + k1k2*yijk)*(k1k3*(-1 + yijk) + k2k3*yijk)*zi*zj)/(-1 + yijk)**2 - (32
     &   *k1k2*k2k3*M2*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi*zj)/(-1 + yijk)**3 + (16*k1k2*k2k3*tdip*((p2k1 
     &   + p2k3)*(-1 + yijk) + p2k2*yijk)*zi*zj)/(-1 + yijk)**3 + (16*k1k2*k2k3*udip*((p2k1 + p2k3)*(-1 + yijk) 
     &   + p2k2*yijk)*zi*zj)/(-1 + yijk)**3 - (8*k2k3*p1k1*sdip*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi
     &   *zj)/(-1 + yijk)**2 - (8*k1k2*p1k3*sdip*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi*zj)/(-1 + yijk)**2 
     &   + (8*k2k3*p2k1*sdip*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi*zj)/(-1 + yijk)**2 + (8*k1k2*p2k3*sdip
     &   *((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi*zj)/(-1 + yijk)**2 + (16*k2k3*M2*(k1k3*(-1 + yijk) + k1k2
     &   *yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi*zj)/(-1 + yijk)**3 - (4*k2k3*tdip*(k1k3*(-1 + yijk) 
     &   + k1k2*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi*zj)/(-1 + yijk)**3 - (12*k2k3*udip*(k1k3*(
     &   -1 + yijk) + k1k2*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi*zj)/(-1 + yijk)**3 + (4*p1k3*sdip*(k1k3
     &   *(-1 + yijk) + k1k2*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi*zj)/(-1 + yijk)**2 - (4*p2k3
     &   *sdip*(k1k3*(-1 + yijk) + k1k2*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi*zj)/(-1 + yijk)**2 + (16
     &   *k1k2*M2*(k1k3*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi*zj)/(-1 + yijk)
     &   **3 - (4*k1k2*tdip*(k1k3*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi*zj)/(-1 
     &   + yijk)**3 - (12*k1k2*udip*(k1k3*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi
     &   *zj)/(-1 + yijk)**3 + (4*p1k1*sdip*(k1k3*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)
     &   *zi*zj)/(-1 + yijk)**2 - (4*p2k1*sdip*(k1k3*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2
     &   *yijk)*zi*zj)/(-1 + yijk)**2 - (16*tdip*(k1k3*(-1 + yijk) + k1k2*yijk)*(k1k3*(-1 + yijk) + k2k3*yijk)
     &   *((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi*zj)/(-1 + yijk)**3 + (16*udip*(k1k3*(-1 + yijk) + k1k2
     &   *yijk)*(k1k3*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zi*zj)/(-1 + yijk)**3 
     &   + (32*k1k2*k2k3*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)**2*zi*zj)/(-1 + yijk)**4 - (16*k2k3*(k1k3*(-1 
     &   + yijk) + k1k2*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)**2*zi*zj)/(-1 + yijk)**4 - (16*k1k2*(k1k3
     &   *(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)**2*zi*zj)/(-1 + yijk)**4 - 2*M2*p1k3
     &   *sdip*(k1k3 + (k1k2*yijk)/(-1 + yijk))*zi*zj + 2*M2*p2k3*sdip*(k1k3 + (k1k2*yijk)/(-1 + yijk))*zi
     &   *zj + 4*p1k3*sdip*tdip*(k1k3 + (k1k2*yijk)/(-1 + yijk))*zi*zj - 4*p2k3*sdip*tdip*(k1k3 + (k1k2*yijk)/(
     &   -1 + yijk))*zi*zj - 2*p1k3*sdip*udip*(k1k3 + (k1k2*yijk)/(-1 + yijk))*zi*zj + 2*p2k3*sdip*udip*(k1k3 
     &   + (k1k2*yijk)/(-1 + yijk))*zi*zj - 2*M2*p1k1*sdip*(k1k3 + (k2k3*yijk)/(-1 + yijk))*zi*zj + 2*M2*p2k1
     &   *sdip*(k1k3 + (k2k3*yijk)/(-1 + yijk))*zi*zj + 4*p1k1*sdip*tdip*(k1k3 + (k2k3*yijk)/(-1 + yijk))*zi*zj 
     &   - 4*p2k1*sdip*tdip*(k1k3 + (k2k3*yijk)/(-1 + yijk))*zi*zj - 2*p1k1*sdip*udip*(k1k3 + (k2k3*yijk)/(-1 
     &   + yijk))*zi*zj + 2*p2k1*sdip*udip*(k1k3 + (k2k3*yijk)/(-1 + yijk))*zi*zj + 4*k1k3*M2*sdip*(p2k1 + p2k3 
     &   + (p2k2*yijk)/(-1 + yijk))*zi*zj - 8*k1k3*sdip*tdip*(p2k1 + p2k3 + (p2k2*yijk)/(-1 + yijk))*zi*zj 
     &   + 4*k1k3*sdip*udip*(p2k1 + p2k3 + (p2k2*yijk)/(-1 + yijk))*zi*zj - 4*k1k3*sdip*(p2k1 + p2k3 + (p2k2*yijk)/(
     &   -1 + yijk))**2*zi*zj - (4*k2k3**2*M2**2*zj**2)/(-1 + yijk)**2 + (16*k2k3**2*M2*sdip*zj**2)/(-1 
     &   + yijk)**2 - (16*k2k3**2*p1p2*sdip*zj**2)/(-1 + yijk)**2 + (4*k2k3**2*M2*tdip*zj**2)/(-1 + yijk)**2 
     &   + (4*k2k3**2*M2*udip*zj**2)/(-1 + yijk)**2 - (4*k2k3**2*tdip*udip*zj**2)/(-1 + yijk)**2 - (4*k2k3*M2*p1k3
     &   *sdip*zj**2)/(-1 + yijk) + (4*k2k3*M2*p2k3*sdip*zj**2)/(-1 + yijk) + (8*k2k3*p1k3*sdip*tdip*zj**2)/(
     &   -1 + yijk) - (8*k2k3*p2k3*sdip*tdip*zj**2)/(-1 + yijk) - (4*k2k3*p1k3*sdip*udip*zj**2)/(-1 + yijk) 
     &   + (4*k2k3*p2k3*sdip*udip*zj**2)/(-1 + yijk) + (4*k2k3*M2**2*(k1k3*(-1 + yijk) + k2k3*yijk)*zj**2)/(-1 
     &   + yijk)**2 - (16*k2k3*M2*sdip*(k1k3*(-1 + yijk) + k2k3*yijk)*zj**2)/(-1 + yijk)**2 + (16*k2k3*p1p2*sdip
     &   *(k1k3*(-1 + yijk) + k2k3*yijk)*zj**2)/(-1 + yijk)**2 - (2*k2k3*M2*tdip*(k1k3*(-1 + yijk) + k2k3*yijk)
     &   *zj**2)/(-1 + yijk)**2 - (4*k2k3*tdip**2*(k1k3*(-1 + yijk) + k2k3*yijk)*zj**2)/(-1 + yijk)**2 - (6
     &   *k2k3*M2*udip*(k1k3*(-1 + yijk) + k2k3*yijk)*zj**2)/(-1 + yijk)**2 + (10*k2k3*tdip*udip*(k1k3*(-1 + yijk) 
     &   + k2k3*yijk)*zj**2)/(-1 + yijk)**2 - (2*k2k3*udip**2*(k1k3*(-1 + yijk) + k2k3*yijk)*zj**2)/(-1 
     &   + yijk)**2 + (16*k2k3**2*M2*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zj**2)/(-1 + yijk)**3 - (8*k2k3*
     &   *2*tdip*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zj**2)/(-1 + yijk)**3 - (8*k2k3**2*udip*((p2k1 + p2k3)
     &   *(-1 + yijk) + p2k2*yijk)*zj**2)/(-1 + yijk)**3 + (8*k2k3*p1k3*sdip*((p2k1 + p2k3)*(-1 + yijk) + p2k2
     &   *yijk)*zj**2)/(-1 + yijk)**2 - (8*k2k3*p2k3*sdip*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zj**2)/(-1 
     &   + yijk)**2 - (16*k2k3*M2*(k1k3*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)*zj**2)/(
     &   -1 + yijk)**3 + (4*k2k3*tdip*(k1k3*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2*yijk)
     &   *zj**2)/(-1 + yijk)**3 + (12*k2k3*udip*(k1k3*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k3)*(-1 + yijk) + p2k2
     &   *yijk)*zj**2)/(-1 + yijk)**3 - (4*p1k3*sdip*(k1k3*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k3)*(-1 + yijk) 
     &   + p2k2*yijk)*zj**2)/(-1 + yijk)**2 + (4*p2k3*sdip*(k1k3*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k3)*(
     &   -1 + yijk) + p2k2*yijk)*zj**2)/(-1 + yijk)**2 + (8*tdip*(k1k3*(-1 + yijk) + k2k3*yijk)**2*((p2k1 + p2k3)
     &   *(-1 + yijk) + p2k2*yijk)*zj**2)/(-1 + yijk)**3 - (8*udip*(k1k3*(-1 + yijk) + k2k3*yijk)**2*((p2k1 
     &   + p2k3)*(-1 + yijk) + p2k2*yijk)*zj**2)/(-1 + yijk)**3 - (16*k2k3**2*((p2k1 + p2k3)*(-1 + yijk) 
     &   + p2k2*yijk)**2*zj**2)/(-1 + yijk)**4 + (16*k2k3*(k1k3*(-1 + yijk) + k2k3*yijk)*((p2k1 + p2k3)*(-1 + yijk) 
     &   + p2k2*yijk)**2*zj**2)/(-1 + yijk)**4 + 2*M2*p1k3*sdip*(k1k3 + (k2k3*yijk)/(-1 + yijk))*zj**2 - 2*M2
     &   *p2k3*sdip*(k1k3 + (k2k3*yijk)/(-1 + yijk))*zj**2 - 4*p1k3*sdip*tdip*(k1k3 + (k2k3*yijk)/(-1 + yijk))*zj**2 
     &   + 4*p2k3*sdip*tdip*(k1k3 + (k2k3*yijk)/(-1 + yijk))*zj**2 + 2*p1k3*sdip*udip*(k1k3 + (k2k3*yijk)/(-1 + yijk))*zj**2 
     &   - 2*p2k3*sdip*udip*(k1k3 + (k2k3*yijk)/(-1 + yijk))*zj**2 + 4*M2*sdip*(k1k3 + (k2k3*yijk)/(-1 + yijk))**2*zj**2 
     &   - 4*p1p2*sdip*(k1k3 + (k2k3*yijk)/(-1 + yijk))**2*zj**2 - 4*M2*tdip*(k1k3 + (k2k3*yijk)/(-1 + yijk))**2*zj**2 +
     &   6*tdip**2*(k1k3 + (k2k3*yijk)/(-1 + yijk))**2*zj**2 + 4*M2*udip*(k1k3 + (k2k3*yijk)/(-1 + yijk))**2*zj**2 
     &   - 8*tdip*udip*(k1k3 + (k2k3*yijk)/(-1 + yijk))**2*zj**2 + 2*udip**2*(k1k3 + (k2k3*yijk)/(-1 + yijk))**2
     &   *zj**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(8*sdip**3))

!************************************************   2 Re(Ms VggP Mt^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_2 = stsT2xx_M2Dipole_SpinTensor_D13_2 + 2*dble(-0.125d0*(GGG*Nc
     &   *(-1 + Nc**2)*(-2*p1k1*(-1 + yijk)*zi + k1k2*yijk*zi + k1k3*(-1 + yijk)*(zi - zj) - 2*p1k3*zj - k2k3*yijk
     &   *zj + 2*p1k3*yijk*zj)*(M2*p2k1*sdip*zi - 2*p2k1**2*sdip*zi - 2*p2k1*p2k3*sdip*zi - 2*p2k1*sdip*tdip
     &   *zi + p2k1*sdip*udip*zi - 3*M2*p2k1*sdip*yijk*zi + 6*p2k1**2*sdip*yijk*zi + 2*p2k1*p2k2*sdip*yijk*zi 
     &   + 6*p2k1*p2k3*sdip*yijk*zi + 6*p2k1*sdip*tdip*yijk*zi - 3*p2k1*sdip*udip*yijk*zi + 3*M2*p2k1*sdip*yijk
     &   **2*zi - 6*p2k1**2*sdip*yijk**2*zi - 4*p2k1*p2k2*sdip*yijk**2*zi - 6*p2k1*p2k3*sdip*yijk**2*zi - 6
     &   *p2k1*sdip*tdip*yijk**2*zi + 3*p2k1*sdip*udip*yijk**2*zi - M2*p2k1*sdip*yijk**3*zi + 2*p2k1**2*sdip*yijk
     &   **3*zi + 2*p2k1*p2k2*sdip*yijk**3*zi + 2*p2k1*p2k3*sdip*yijk**3*zi + 2*p2k1*sdip*tdip*yijk**3*zi - p2k1
     &   *sdip*udip*yijk**3*zi + k1k2*(4*p2k3**2 + 4*p1p2*sdip - 4*p2k3*tdip - 2*tdip**2 + 2*p2k3*udip + 4
     &   *p2k1**2*(-1 + yijk)**2 + 2*M2**2*(-1 + yijk)**3 - 8*p2k2*p2k3*yijk - 8*p2k3**2*yijk - 10*p1p2*sdip*yijk 
     &   + 4*p2k2*tdip*yijk + 14*p2k3*tdip*yijk + 7*tdip**2*yijk - 2*p2k2*udip*yijk - 6*p2k3*udip*yijk - tdip
     &   *udip*yijk + 4*p2k2**2*yijk**2 + 8*p2k2*p2k3*yijk**2 + 4*p2k3**2*yijk**2 + 8*p1p2*sdip*yijk**2 - 10
     &   *p2k2*tdip*yijk**2 - 16*p2k3*tdip*yijk**2 - 8*tdip**2*yijk**2 + 4*p2k2*udip*yijk**2 + 6*p2k3*udip*yijk
     &   **2 + 2*tdip*udip*yijk**2 - 2*p1p2*sdip*yijk**3 + 6*p2k2*tdip*yijk**3 + 6*p2k3*tdip*yijk**3 + 3*tdip
     &   **2*yijk**3 - 2*p2k2*udip*yijk**3 - 2*p2k3*udip*yijk**3 - tdip*udip*yijk**3 + 2*p2k1*(-1 + yijk)*(-udip 
     &   + 4*p2k3*(-1 + yijk) + 4*p2k2*yijk + 2*udip*yijk - udip*yijk**2 + tdip*(2 - 5*yijk + 3*yijk**2)) 
     &   - M2*(-1 + yijk)*(-4*sdip + 4*tdip - 2*p2k2*yijk + 6*sdip*yijk - 9*tdip*yijk + udip*yijk + 4*p2k2*yijk
     &   **2 - 2*sdip*yijk**2 + 5*tdip*yijk**2 - udip*yijk**2 + p2k1*(2 - 6*yijk + 4*yijk**2) + p2k3*(2 - 6*yijk 
     &   + 4*yijk**2)))*zi + k1k3*(-1 + yijk)**2*(2*M2**2*(-1 + yijk) - 2*p1p2*sdip*(-1 + yijk) + (3*tdip 
     &   - udip)*(-tdip + 2*p2k1*(-1 + yijk) + 2*p2k3*(-1 + yijk) + 2*p2k2*yijk + tdip*yijk) - M2*(2*sdip - 5
     &   *tdip + udip + 4*p2k1*(-1 + yijk) + 4*p2k3*(-1 + yijk) + 4*p2k2*yijk - 2*sdip*yijk + 5*tdip*yijk - udip
     &   *yijk))*(zi - zj) + 2*k2k3*M2**2*zj - 2*k2k3*M2*p2k1*zj - 4*k2k3*p2k1**2*zj - 2*k2k3*M2*p2k3*zj - 8
     &   *k2k3*p2k1*p2k3*zj - 4*k2k3*p2k3**2*zj + 4*k2k3*M2*sdip*zj - 4*k2k3*p1p2*sdip*zj - M2*p2k3*sdip*zj + 2
     &   *p2k1*p2k3*sdip*zj + 2*p2k3**2*sdip*zj - 4*k2k3*M2*tdip*zj + 4*k2k3*p2k1*tdip*zj + 4*k2k3*p2k3*tdip*zj 
     &   + 2*p2k3*sdip*tdip*zj + 2*k2k3*tdip**2*zj - 2*k2k3*p2k1*udip*zj - 2*k2k3*p2k3*udip*zj - p2k3*sdip*udip
     &   *zj - 6*k2k3*M2**2*yijk*zj + 8*k2k3*M2*p2k1*yijk*zj + 8*k2k3*p2k1**2*yijk*zj + 2*k2k3*M2*p2k2*yijk
     &   *zj + 8*k2k3*p2k1*p2k2*yijk*zj + 8*k2k3*M2*p2k3*yijk*zj + 16*k2k3*p2k1*p2k3*yijk*zj + 8*k2k3*p2k2*p2k3
     &   *yijk*zj + 8*k2k3*p2k3**2*yijk*zj - 10*k2k3*M2*sdip*yijk*zj + 10*k2k3*p1p2*sdip*yijk*zj + 3*M2*p2k3
     &   *sdip*yijk*zj - 6*p2k1*p2k3*sdip*yijk*zj - 2*p2k2*p2k3*sdip*yijk*zj - 6*p2k3**2*sdip*yijk*zj + 13*k2k3
     &   *M2*tdip*yijk*zj - 14*k2k3*p2k1*tdip*yijk*zj - 4*k2k3*p2k2*tdip*yijk*zj - 14*k2k3*p2k3*tdip*yijk*zj 
     &   - 6*p2k3*sdip*tdip*yijk*zj - 7*k2k3*tdip**2*yijk*zj - k2k3*M2*udip*yijk*zj + 6*k2k3*p2k1*udip*yijk*zj 
     &   + 2*k2k3*p2k2*udip*yijk*zj + 6*k2k3*p2k3*udip*yijk*zj + 3*p2k3*sdip*udip*yijk*zj + k2k3*tdip*udip*yijk
     &   *zj + 6*k2k3*M2**2*yijk**2*zj - 10*k2k3*M2*p2k1*yijk**2*zj - 4*k2k3*p2k1**2*yijk**2*zj - 6*k2k3*M2*p2k2
     &   *yijk**2*zj - 8*k2k3*p2k1*p2k2*yijk**2*zj - 4*k2k3*p2k2**2*yijk**2*zj - 10*k2k3*M2*p2k3*yijk**2*zj 
     &   - 8*k2k3*p2k1*p2k3*yijk**2*zj - 8*k2k3*p2k2*p2k3*yijk**2*zj - 4*k2k3*p2k3**2*yijk**2*zj + 8*k2k3*M2
     &   *sdip*yijk**2*zj - 8*k2k3*p1p2*sdip*yijk**2*zj - 3*M2*p2k3*sdip*yijk**2*zj + 6*p2k1*p2k3*sdip*yijk**2
     &   *zj + 4*p2k2*p2k3*sdip*yijk**2*zj + 6*p2k3**2*sdip*yijk**2*zj - 14*k2k3*M2*tdip*yijk**2*zj + 16*k2k3*p2k1
     &   *tdip*yijk**2*zj + 10*k2k3*p2k2*tdip*yijk**2*zj + 16*k2k3*p2k3*tdip*yijk**2*zj + 6*p2k3*sdip*tdip
     &   *yijk**2*zj + 8*k2k3*tdip**2*yijk**2*zj + 2*k2k3*M2*udip*yijk**2*zj - 6*k2k3*p2k1*udip*yijk**2*zj - 4
     &   *k2k3*p2k2*udip*yijk**2*zj - 6*k2k3*p2k3*udip*yijk**2*zj - 3*p2k3*sdip*udip*yijk**2*zj - 2*k2k3*tdip*udip
     &   *yijk**2*zj - 2*k2k3*M2**2*yijk**3*zj + 4*k2k3*M2*p2k1*yijk**3*zj + 4*k2k3*M2*p2k2*yijk**3*zj + 4
     &   *k2k3*M2*p2k3*yijk**3*zj - 2*k2k3*M2*sdip*yijk**3*zj + 2*k2k3*p1p2*sdip*yijk**3*zj + M2*p2k3*sdip*yijk
     &   **3*zj - 2*p2k1*p2k3*sdip*yijk**3*zj - 2*p2k2*p2k3*sdip*yijk**3*zj - 2*p2k3**2*sdip*yijk**3*zj + 5*k2k3
     &   *M2*tdip*yijk**3*zj - 6*k2k3*p2k1*tdip*yijk**3*zj - 6*k2k3*p2k2*tdip*yijk**3*zj - 6*k2k3*p2k3*tdip*yijk**3*zj 
     &   - 2*p2k3*sdip*tdip*yijk**3*zj - 3*k2k3*tdip**2*yijk**3*zj - k2k3*M2*udip*yijk**3*zj + 2*k2k3*p2k1*udip*yijk**3*zj 
     &   + 2*k2k3*p2k2*udip*yijk**3*zj + 2*k2k3*p2k3*udip*yijk**3*zj + p2k3*sdip*udip*yijk**3*zj + k2k3*tdip*udip*yijk**3*zj)
     &   *SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(sdip**2*(M2 - tdip)*(-1 + yijk)**4))

!************************************************   2 Re(Ms VggP Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_2 = stsT2xx_M2Dipole_SpinTensor_D13_2 + 2*dble(-0.125d0*(GGG*Nc
     &   *(-1 + Nc**2)*(-2*p2k1*(-1 + yijk)*zi + k1k2*yijk*zi + k1k3*(-1 + yijk)*(zi - zj) - 2*p2k3*zj - k2k3*yijk
     &   *zj + 2*p2k3*yijk*zj)*(M2*p1k1*sdip*zi - 2*p1k1*p2k1*sdip*zi - 2*p1k1*p2k3*sdip*zi - 2*p1k1*sdip*tdip
     &   *zi + p1k1*sdip*udip*zi - 2*M2*p1k1*sdip*yijk*zi + 4*p1k1*p2k1*sdip*yijk*zi + 2*p1k1*p2k2*sdip*yijk
     &   *zi + 4*p1k1*p2k3*sdip*yijk*zi + 4*p1k1*sdip*tdip*yijk*zi - 2*p1k1*sdip*udip*yijk*zi + M2*p1k1*sdip
     &   *yijk**2*zi - 2*p1k1*p2k1*sdip*yijk**2*zi - 2*p1k1*p2k2*sdip*yijk**2*zi - 2*p1k1*p2k3*sdip*yijk**2*zi 
     &   - 2*p1k1*sdip*tdip*yijk**2*zi + p1k1*sdip*udip*yijk**2*zi + k1k2*(-2*p2k1*tdip - 2*p2k3*tdip - 4*p2k1
     &   *udip - 4*p2k3*udip - 3*tdip*udip + 2*p2k1*tdip*yijk + 2*p2k2*tdip*yijk + 2*p2k3*tdip*yijk - tdip**2
     &   *yijk + 8*p2k1*udip*yijk + 4*p2k2*udip*yijk + 8*p2k3*udip*yijk + 7*tdip*udip*yijk - udip**2*yijk + tdip
     &   **2*yijk**2 - 4*p2k1*udip*yijk**2 - 4*p2k2*udip*yijk**2 - 4*p2k3*udip*yijk**2 - 4*tdip*udip*yijk**2 
     &   + udip**2*yijk**2 + M2**2*(-3 + 5*yijk - 2*yijk**2) - 2*p1p2*sdip*(2 - 3*yijk + yijk**2) + M2*(4*sdip 
     &   + 3*tdip + 3*udip - 6*p2k2*yijk - 6*sdip*yijk - 5*tdip*yijk - 5*udip*yijk + 4*p2k2*yijk**2 + 2*sdip
     &   *yijk**2 + 2*tdip*yijk**2 + 2*udip*yijk**2 + 2*p2k1*(3 - 5*yijk + 2*yijk**2) + 2*p2k3*(3 - 5*yijk + 2
     &   *yijk**2)))*zi - k1k3*(-1 + yijk)*(tdip**2 - 4*p2k1*udip - 4*p2k3*udip - 4*tdip*udip + udip**2 + 2*M2
     &   **2*(-1 + yijk) + 2*p1p2*sdip*(-1 + yijk) - tdip**2*yijk + 4*p2k1*udip*yijk + 4*p2k2*udip*yijk + 4*p2k3
     &   *udip*yijk + 4*tdip*udip*yijk - udip**2*yijk - 2*M2*(-sdip - tdip - udip + 2*p2k1*(-1 + yijk) + 2*p2k3
     &   *(-1 + yijk) + 2*p2k2*yijk + sdip*yijk + tdip*yijk + udip*yijk))*(zi - zj) + 3*k2k3*M2**2*zj - 6*k2k3
     &   *M2*p2k1*zj - 6*k2k3*M2*p2k3*zj - 4*k2k3*M2*sdip*zj - M2*p1k3*sdip*zj + 4*k2k3*p1p2*sdip*zj + 2*p1k3
     &   *p2k1*sdip*zj + 2*p1k3*p2k3*sdip*zj - 3*k2k3*M2*tdip*zj + 2*k2k3*p2k1*tdip*zj + 2*k2k3*p2k3*tdip*zj 
     &   + 2*p1k3*sdip*tdip*zj - 3*k2k3*M2*udip*zj + 4*k2k3*p2k1*udip*zj + 4*k2k3*p2k3*udip*zj - p1k3*sdip*udip
     &   *zj + 3*k2k3*tdip*udip*zj - 5*k2k3*M2**2*yijk*zj + 10*k2k3*M2*p2k1*yijk*zj + 6*k2k3*M2*p2k2*yijk*zj 
     &   + 10*k2k3*M2*p2k3*yijk*zj + 6*k2k3*M2*sdip*yijk*zj + 2*M2*p1k3*sdip*yijk*zj - 6*k2k3*p1p2*sdip*yijk*zj 
     &   - 4*p1k3*p2k1*sdip*yijk*zj - 2*p1k3*p2k2*sdip*yijk*zj - 4*p1k3*p2k3*sdip*yijk*zj + 5*k2k3*M2*tdip*yijk
     &   *zj - 2*k2k3*p2k1*tdip*yijk*zj - 2*k2k3*p2k2*tdip*yijk*zj - 2*k2k3*p2k3*tdip*yijk*zj - 4*p1k3*sdip
     &   *tdip*yijk*zj + k2k3*tdip**2*yijk*zj + 5*k2k3*M2*udip*yijk*zj - 8*k2k3*p2k1*udip*yijk*zj - 4*k2k3*p2k2
     &   *udip*yijk*zj - 8*k2k3*p2k3*udip*yijk*zj + 2*p1k3*sdip*udip*yijk*zj - 7*k2k3*tdip*udip*yijk*zj + k2k3
     &   *udip**2*yijk*zj + 2*k2k3*M2**2*yijk**2*zj - 4*k2k3*M2*p2k1*yijk**2*zj - 4*k2k3*M2*p2k2*yijk**2*zj 
     &   - 4*k2k3*M2*p2k3*yijk**2*zj - 2*k2k3*M2*sdip*yijk**2*zj - M2*p1k3*sdip*yijk**2*zj + 2*k2k3*p1p2*sdip*yijk
     &   **2*zj + 2*p1k3*p2k1*sdip*yijk**2*zj + 2*p1k3*p2k2*sdip*yijk**2*zj + 2*p1k3*p2k3*sdip*yijk**2*zj 
     &   - 2*k2k3*M2*tdip*yijk**2*zj + 2*p1k3*sdip*tdip*yijk**2*zj - k2k3*tdip**2*yijk**2*zj - 2*k2k3*M2*udip*yijk**2*zj 
     &   + 4*k2k3*p2k1*udip*yijk**2*zj + 4*k2k3*p2k2*udip*yijk**2*zj + 4*k2k3*p2k3*udip*yijk**2*zj - p1k3*sdip*udip*yijk**2*zj 
     &   + 4*k2k3*tdip*udip*yijk**2*zj - k2k3*udip**2*yijk**2*zj)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(sdip**2*(M2 - udip)*(-1 + yijk)**3))

!************************************************   2 Re(Ms VggP Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_2 = stsT2xx_M2Dipole_SpinTensor_D13_2 + 2*dble(0)

!************************************************   Mt^ast VggP Mt   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_2 = stsT2xx_M2Dipole_SpinTensor_D13_2 +((Nc*(-1 + Nc**2)*(M2*sdip
     &   *(-1 + yijk) - 2*(p2k1 + p2k3)*(M2 - tdip)*(-1 + yijk) + 2*p2k2*(-M2 + tdip)*yijk)*(-2*p1k1*(-1 + yijk)
     &   *zi + k1k2*yijk*zi + k1k3*(-1 + yijk)*(zi - zj) - 2*p1k3*zj - k2k3*yijk*zj + 2*p1k3*yijk*zj)**2
     &   *SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(2*sdip*(M2 - tdip)**2*(-1 + yijk)**3))

!************************************************   2 Re(Mt VggP Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_2 = stsT2xx_M2Dipole_SpinTensor_D13_2 + 2*dble(0)

!************************************************   2 Re(Mt VggP Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_2 = stsT2xx_M2Dipole_SpinTensor_D13_2 + 2*dble((Nc*(-1 + Nc**2)
     &   *(-2*p1k1*(-1 + yijk)*zi + k1k2*yijk*zi + k1k3*(-1 + yijk)*(zi - zj) - 2*p1k3*zj - k2k3*yijk*zj + 2
     &   *p1k3*yijk*zj)*(p2k1*sdip*zi - 2*p2k1*sdip*yijk*zi + p2k1*sdip*yijk**2*zi + k1k2*(2*p2k1*(-1 + yijk) 
     &   + 2*p2k3*(-1 + yijk) + yijk*(M2 + 2*p2k2 + tdip*(-1 + yijk) - M2*yijk))*zi - k1k3*(M2 - tdip)*(-1 + yijk)
     &   **2*(zi - zj) + 2*k2k3*p2k1*zj + 2*k2k3*p2k3*zj - p2k3*sdip*zj - k2k3*M2*yijk*zj - 2*k2k3*p2k1*yijk
     &   *zj - 2*k2k3*p2k2*yijk*zj - 2*k2k3*p2k3*yijk*zj + 2*p2k3*sdip*yijk*zj + k2k3*tdip*yijk*zj + k2k3*M2
     &   *yijk**2*zj - p2k3*sdip*yijk**2*zj - k2k3*tdip*yijk**2*zj)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(4.
     &   *sdip*(M2 - tdip)*(-1 + yijk)**3))

!************************************************   Mu^ast VggP Mu   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_2 = stsT2xx_M2Dipole_SpinTensor_D13_2 +((Nc*(-1 + Nc**2)*(-M2
     &   **2 - tdip*udip + M2*(sdip + tdip + udip))*(-2*p2k1*(-1 + yijk)*zi + k1k2*yijk*zi + k1k3*(-1 + yijk)*(zi 
     &   - zj) - 2*p2k3*zj - k2k3*yijk*zj + 2*p2k3*yijk*zj)**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(2
     &   *sdip*(M2 - udip)**2*(-1 + yijk)**2))

!************************************************   2 Re(Mu VggP Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_2 = stsT2xx_M2Dipole_SpinTensor_D13_2 + 2*dble(-0.25d0*(Nc*(-1 
     &   + Nc**2)*(-2*p2k1*(-1 + yijk)*zi + k1k2*yijk*zi + k1k3*(-1 + yijk)*(zi - zj) - 2*p2k3*zj - k2k3*yijk*zj 
     &   + 2*p2k3*yijk*zj)*(p1k1*sdip*zi - p1k1*sdip*yijk*zi + k1k2*(tdip + M2*(-1 + yijk) - udip*yijk)*zi 
     &   + k1k3*(M2 - udip)*(-1 + yijk)*(zi - zj) + k2k3*M2*zj - p1k3*sdip*zj - k2k3*tdip*zj - k2k3*M2*yijk*zj 
     &   + p1k3*sdip*yijk*zj + k2k3*udip*yijk*zj)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(sdip
     &   *(M2 - udip)*(-1 + yijk)**2))

!************************************************   Mv^ast VggP Mv   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_2 = stsT2xx_M2Dipole_SpinTensor_D13_2 +((Nc*(-1 + Nc**2)*(2*k1k2
     &   **2*yijk*zi**2 + 2*k1k2*zi*(k1k3*(-1 + yijk)*(zi - zj) - 2*k2k3*yijk*zj) + zj*(-(k1k3*(-1 + yijk)*(sdip
     &   *(-1 + yijk)*zi + 2*k2k3*(zi - zj))) + 2*k2k3**2*yijk*zj))*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **2)/(2*sdip*(-1 + yijk)**2))



      end



      double precision function stsT2xx_M2Dipole_SpinTensor_D23_1(sdip,tdip,udip)
      implicit none
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"
#include "common/stsT2xx_GenCouplings.h"

      double precision sdip,tdip,udip,M2
      double precision yijk,zi,zj

      M2=MSf2(iSf1,itt1,igen1)

      yijk = k2k3/(k1k2 + k1k3 + k2k3)
      zi = k1k2/(k1k2 + k1k3)
      zj = 1-zi

!************************************************   Ms^ast VggP Ms   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_1 = ((GGG**2*Nc*(-1 + Nc**2)
     &   *((-4*k1k2**2*M2**2*zi**2)/(-1 + yijk)**2 + (16*k1k2**2*M2*sdip*zi**2)/(-1 + yijk)**2 - (16*k1k2**2
     &   *p1p2*sdip*zi**2)/(-1 + yijk)**2 + (4*k1k2**2*M2*tdip*zi**2)/(-1 + yijk)**2 + (4*k1k2**2*M2*udip*zi**2)/(
     &   -1 + yijk)**2 - (4*k1k2**2*tdip*udip*zi**2)/(-1 + yijk)**2 - (4*k1k2*M2*p1k2*sdip*zi**2)/(-1 + yijk) 
     &   + (4*k1k2*M2*p2k2*sdip*zi**2)/(-1 + yijk) + (8*k1k2*p1k2*sdip*tdip*zi**2)/(-1 + yijk) - (8*k1k2*p2k2
     &   *sdip*tdip*zi**2)/(-1 + yijk) - (4*k1k2*p1k2*sdip*udip*zi**2)/(-1 + yijk) + (4*k1k2*p2k2*sdip*udip
     &   *zi**2)/(-1 + yijk) + (4*k1k2*M2**2*(k2k3*(-1 + yijk) + k1k2*yijk)*zi**2)/(-1 + yijk)**2 - (16*k1k2*M2
     &   *sdip*(k2k3*(-1 + yijk) + k1k2*yijk)*zi**2)/(-1 + yijk)**2 + (16*k1k2*p1p2*sdip*(k2k3*(-1 + yijk) + k1k2
     &   *yijk)*zi**2)/(-1 + yijk)**2 - (2*k1k2*M2*tdip*(k2k3*(-1 + yijk) + k1k2*yijk)*zi**2)/(-1 + yijk)*
     &   *2 - (4*k1k2*tdip**2*(k2k3*(-1 + yijk) + k1k2*yijk)*zi**2)/(-1 + yijk)**2 - (6*k1k2*M2*udip*(k2k3*(-1 
     &   + yijk) + k1k2*yijk)*zi**2)/(-1 + yijk)**2 + (10*k1k2*tdip*udip*(k2k3*(-1 + yijk) + k1k2*yijk)*zi**2)/(
     &   -1 + yijk)**2 - (2*k1k2*udip**2*(k2k3*(-1 + yijk) + k1k2*yijk)*zi**2)/(-1 + yijk)**2 + (16*k1k2**2
     &   *M2*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi**2)/(-1 + yijk)**3 - (8*k1k2**2*tdip*((p2k2 + p2k3)*(-1 
     &   + yijk) + p2k1*yijk)*zi**2)/(-1 + yijk)**3 - (8*k1k2**2*udip*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)
     &   *zi**2)/(-1 + yijk)**3 + (8*k1k2*p1k2*sdip*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi**2)/(-1 + yijk)
     &   **2 - (8*k1k2*p2k2*sdip*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi**2)/(-1 + yijk)**2 - (16*k1k2*M2*(k2k3
     &   *(-1 + yijk) + k1k2*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi**2)/(-1 + yijk)**3 + (4*k1k2
     &   *tdip*(k2k3*(-1 + yijk) + k1k2*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi**2)/(-1 + yijk)**3 + (12
     &   *k1k2*udip*(k2k3*(-1 + yijk) + k1k2*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi**2)/(-1 + yijk)
     &   **3 - (4*p1k2*sdip*(k2k3*(-1 + yijk) + k1k2*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi**2)/(
     &   -1 + yijk)**2 + (4*p2k2*sdip*(k2k3*(-1 + yijk) + k1k2*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi
     &   **2)/(-1 + yijk)**2 + (8*tdip*(k2k3*(-1 + yijk) + k1k2*yijk)**2*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)
     &   *zi**2)/(-1 + yijk)**3 - (8*udip*(k2k3*(-1 + yijk) + k1k2*yijk)**2*((p2k2 + p2k3)*(-1 + yijk) + p2k1
     &   *yijk)*zi**2)/(-1 + yijk)**3 - (16*k1k2**2*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)**2*zi**2)/(-1 + yijk)
     &   **4 + (16*k1k2*(k2k3*(-1 + yijk) + k1k2*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)**2*zi**2)/(
     &   -1 + yijk)**4 + 2*M2*p1k2*sdip*(k2k3 + (k1k2*yijk)/(-1 + yijk))*zi**2 - 2*M2*p2k2*sdip*(k2k3 + (k1k2
     &   *yijk)/(-1 + yijk))*zi**2 - 4*p1k2*sdip*tdip*(k2k3 + (k1k2*yijk)/(-1 + yijk))*zi**2 + 4*p2k2*sdip*tdip
     &   *(k2k3 + (k1k2*yijk)/(-1 + yijk))*zi**2 + 2*p1k2*sdip*udip*(k2k3 + (k1k2*yijk)/(-1 + yijk))*zi**2 - 2
     &   *p2k2*sdip*udip*(k2k3 + (k1k2*yijk)/(-1 + yijk))*zi**2 + 4*M2*sdip*(k2k3 + (k1k2*yijk)/(-1 + yijk))*
     &   *2*zi**2 - 4*p1p2*sdip*(k2k3 + (k1k2*yijk)/(-1 + yijk))**2*zi**2 - 4*M2*tdip*(k2k3 + (k1k2*yijk)/(-1 
     &   + yijk))**2*zi**2 + 6*tdip**2*(k2k3 + (k1k2*yijk)/(-1 + yijk))**2*zi**2 + 4*M2*udip*(k2k3 + (k1k2*yijk)/(
     &   -1 + yijk))**2*zi**2 - 8*tdip*udip*(k2k3 + (k1k2*yijk)/(-1 + yijk))**2*zi**2 + 2*udip**2*(k2k3 + (k1k2
     &   *yijk)/(-1 + yijk))**2*zi**2 - k2k3*M2**2*sdip*zi*zj + 4*k2k3*M2*sdip*tdip*zi*zj - 4*k2k3*sdip*tdip
     &   **2*zi*zj - 2*k2k3*M2*sdip*udip*zi*zj + 4*k2k3*sdip*tdip*udip*zi*zj - k2k3*sdip*udip**2*zi*zj + (8
     &   *k1k2*k1k3*M2**2*zi*zj)/(-1 + yijk)**2 - (32*k1k2*k1k3*M2*sdip*zi*zj)/(-1 + yijk)**2 + (32*k1k2*k1k3*p1p2
     &   *sdip*zi*zj)/(-1 + yijk)**2 - (8*k1k2*k1k3*M2*tdip*zi*zj)/(-1 + yijk)**2 - (8*k1k2*k1k3*M2*udip*zi
     &   *zj)/(-1 + yijk)**2 + (8*k1k2*k1k3*tdip*udip*zi*zj)/(-1 + yijk)**2 + (4*k1k3*M2*p1k2*sdip*zi*zj)/(-1 
     &   + yijk) + (4*k1k2*M2*p1k3*sdip*zi*zj)/(-1 + yijk) - (4*k1k3*M2*p2k2*sdip*zi*zj)/(-1 + yijk) - (4*k1k2
     &   *M2*p2k3*sdip*zi*zj)/(-1 + yijk) - (8*k1k3*p1k2*sdip*tdip*zi*zj)/(-1 + yijk) - (8*k1k2*p1k3*sdip*tdip
     &   *zi*zj)/(-1 + yijk) + (8*k1k3*p2k2*sdip*tdip*zi*zj)/(-1 + yijk) + (8*k1k2*p2k3*sdip*tdip*zi*zj)/(-1 
     &   + yijk) + (4*k1k3*p1k2*sdip*udip*zi*zj)/(-1 + yijk) + (4*k1k2*p1k3*sdip*udip*zi*zj)/(-1 + yijk) - (4*k1k3
     &   *p2k2*sdip*udip*zi*zj)/(-1 + yijk) - (4*k1k2*p2k3*sdip*udip*zi*zj)/(-1 + yijk) - (4*k1k3*M2**2*(k2k3
     &   *(-1 + yijk) + k1k2*yijk)*zi*zj)/(-1 + yijk)**2 + (16*k1k3*M2*sdip*(k2k3*(-1 + yijk) + k1k2*yijk)*zi
     &   *zj)/(-1 + yijk)**2 - (16*k1k3*p1p2*sdip*(k2k3*(-1 + yijk) + k1k2*yijk)*zi*zj)/(-1 + yijk)**2 + (2*k1k3
     &   *M2*tdip*(k2k3*(-1 + yijk) + k1k2*yijk)*zi*zj)/(-1 + yijk)**2 + (4*k1k3*tdip**2*(k2k3*(-1 + yijk) 
     &   + k1k2*yijk)*zi*zj)/(-1 + yijk)**2 + (6*k1k3*M2*udip*(k2k3*(-1 + yijk) + k1k2*yijk)*zi*zj)/(-1 + yijk)
     &   **2 - (10*k1k3*tdip*udip*(k2k3*(-1 + yijk) + k1k2*yijk)*zi*zj)/(-1 + yijk)**2 + (2*k1k3*udip**2*(k2k3
     &   *(-1 + yijk) + k1k2*yijk)*zi*zj)/(-1 + yijk)**2 - (4*k1k2*M2**2*(k2k3*(-1 + yijk) + k1k3*yijk)*zi*zj)/(
     &   -1 + yijk)**2 + (16*k1k2*M2*sdip*(k2k3*(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 + yijk)**2 - (16*k1k2*p1p2
     &   *sdip*(k2k3*(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 + yijk)**2 + (2*k1k2*M2*tdip*(k2k3*(-1 + yijk) + k1k3
     &   *yijk)*zi*zj)/(-1 + yijk)**2 + (4*k1k2*tdip**2*(k2k3*(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 + yijk)*
     &   *2 + (6*k1k2*M2*udip*(k2k3*(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 + yijk)**2 - (10*k1k2*tdip*udip*(k2k3*(
     &   -1 + yijk) + k1k3*yijk)*zi*zj)/(-1 + yijk)**2 + (2*k1k2*udip**2*(k2k3*(-1 + yijk) + k1k3*yijk)*zi*zj)/(
     &   -1 + yijk)**2 - (8*M2*sdip*(k2k3*(-1 + yijk) + k1k2*yijk)*(k2k3*(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 
     &   + yijk)**2 + (8*p1p2*sdip*(k2k3*(-1 + yijk) + k1k2*yijk)*(k2k3*(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 
     &   + yijk)**2 + (8*M2*tdip*(k2k3*(-1 + yijk) + k1k2*yijk)*(k2k3*(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 + yijk)
     &   **2 - (12*tdip**2*(k2k3*(-1 + yijk) + k1k2*yijk)*(k2k3*(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 + yijk)
     &   **2 - (8*M2*udip*(k2k3*(-1 + yijk) + k1k2*yijk)*(k2k3*(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 + yijk)**2 
     &   + (16*tdip*udip*(k2k3*(-1 + yijk) + k1k2*yijk)*(k2k3*(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 + yijk)**2 - (4
     &   *udip**2*(k2k3*(-1 + yijk) + k1k2*yijk)*(k2k3*(-1 + yijk) + k1k3*yijk)*zi*zj)/(-1 + yijk)**2 - (32
     &   *k1k2*k1k3*M2*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi*zj)/(-1 + yijk)**3 + (16*k1k2*k1k3*tdip*((p2k2 
     &   + p2k3)*(-1 + yijk) + p2k1*yijk)*zi*zj)/(-1 + yijk)**3 + (16*k1k2*k1k3*udip*((p2k2 + p2k3)*(-1 + yijk) 
     &   + p2k1*yijk)*zi*zj)/(-1 + yijk)**3 - (8*k1k3*p1k2*sdip*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi
     &   *zj)/(-1 + yijk)**2 - (8*k1k2*p1k3*sdip*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi*zj)/(-1 + yijk)**2 
     &   + (8*k1k3*p2k2*sdip*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi*zj)/(-1 + yijk)**2 + (8*k1k2*p2k3*sdip
     &   *((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi*zj)/(-1 + yijk)**2 + (16*k1k3*M2*(k2k3*(-1 + yijk) + k1k2
     &   *yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi*zj)/(-1 + yijk)**3 - (4*k1k3*tdip*(k2k3*(-1 + yijk) 
     &   + k1k2*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi*zj)/(-1 + yijk)**3 - (12*k1k3*udip*(k2k3*(
     &   -1 + yijk) + k1k2*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi*zj)/(-1 + yijk)**3 + (4*p1k3*sdip*(k2k3
     &   *(-1 + yijk) + k1k2*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi*zj)/(-1 + yijk)**2 - (4*p2k3
     &   *sdip*(k2k3*(-1 + yijk) + k1k2*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi*zj)/(-1 + yijk)**2 + (16
     &   *k1k2*M2*(k2k3*(-1 + yijk) + k1k3*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi*zj)/(-1 + yijk)
     &   **3 - (4*k1k2*tdip*(k2k3*(-1 + yijk) + k1k3*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi*zj)/(-1 
     &   + yijk)**3 - (12*k1k2*udip*(k2k3*(-1 + yijk) + k1k3*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi
     &   *zj)/(-1 + yijk)**3 + (4*p1k2*sdip*(k2k3*(-1 + yijk) + k1k3*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)
     &   *zi*zj)/(-1 + yijk)**2 - (4*p2k2*sdip*(k2k3*(-1 + yijk) + k1k3*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1
     &   *yijk)*zi*zj)/(-1 + yijk)**2 - (16*tdip*(k2k3*(-1 + yijk) + k1k2*yijk)*(k2k3*(-1 + yijk) + k1k3*yijk)
     &   *((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi*zj)/(-1 + yijk)**3 + (16*udip*(k2k3*(-1 + yijk) + k1k2
     &   *yijk)*(k2k3*(-1 + yijk) + k1k3*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zi*zj)/(-1 + yijk)**3 
     &   + (32*k1k2*k1k3*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)**2*zi*zj)/(-1 + yijk)**4 - (16*k1k3*(k2k3*(-1 
     &   + yijk) + k1k2*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)**2*zi*zj)/(-1 + yijk)**4 - (16*k1k2*(k2k3
     &   *(-1 + yijk) + k1k3*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)**2*zi*zj)/(-1 + yijk)**4 - 2*M2*p1k3
     &   *sdip*(k2k3 + (k1k2*yijk)/(-1 + yijk))*zi*zj + 2*M2*p2k3*sdip*(k2k3 + (k1k2*yijk)/(-1 + yijk))*zi
     &   *zj + 4*p1k3*sdip*tdip*(k2k3 + (k1k2*yijk)/(-1 + yijk))*zi*zj - 4*p2k3*sdip*tdip*(k2k3 + (k1k2*yijk)/(
     &   -1 + yijk))*zi*zj - 2*p1k3*sdip*udip*(k2k3 + (k1k2*yijk)/(-1 + yijk))*zi*zj + 2*p2k3*sdip*udip*(k2k3 
     &   + (k1k2*yijk)/(-1 + yijk))*zi*zj - 2*M2*p1k2*sdip*(k2k3 + (k1k3*yijk)/(-1 + yijk))*zi*zj + 2*M2*p2k2
     &   *sdip*(k2k3 + (k1k3*yijk)/(-1 + yijk))*zi*zj + 4*p1k2*sdip*tdip*(k2k3 + (k1k3*yijk)/(-1 + yijk))*zi*zj 
     &   - 4*p2k2*sdip*tdip*(k2k3 + (k1k3*yijk)/(-1 + yijk))*zi*zj - 2*p1k2*sdip*udip*(k2k3 + (k1k3*yijk)/(-1 
     &   + yijk))*zi*zj + 2*p2k2*sdip*udip*(k2k3 + (k1k3*yijk)/(-1 + yijk))*zi*zj + 4*k2k3*M2*sdip*(p2k2 + p2k3 
     &   + (p2k1*yijk)/(-1 + yijk))*zi*zj - 8*k2k3*sdip*tdip*(p2k2 + p2k3 + (p2k1*yijk)/(-1 + yijk))*zi*zj 
     &   + 4*k2k3*sdip*udip*(p2k2 + p2k3 + (p2k1*yijk)/(-1 + yijk))*zi*zj - 4*k2k3*sdip*(p2k2 + p2k3 + (p2k1*yijk)/(
     &   -1 + yijk))**2*zi*zj - (4*k1k3**2*M2**2*zj**2)/(-1 + yijk)**2 + (16*k1k3**2*M2*sdip*zj**2)/(-1 
     &   + yijk)**2 - (16*k1k3**2*p1p2*sdip*zj**2)/(-1 + yijk)**2 + (4*k1k3**2*M2*tdip*zj**2)/(-1 + yijk)**2 
     &   + (4*k1k3**2*M2*udip*zj**2)/(-1 + yijk)**2 - (4*k1k3**2*tdip*udip*zj**2)/(-1 + yijk)**2 - (4*k1k3*M2*p1k3
     &   *sdip*zj**2)/(-1 + yijk) + (4*k1k3*M2*p2k3*sdip*zj**2)/(-1 + yijk) + (8*k1k3*p1k3*sdip*tdip*zj**2)/(
     &   -1 + yijk) - (8*k1k3*p2k3*sdip*tdip*zj**2)/(-1 + yijk) - (4*k1k3*p1k3*sdip*udip*zj**2)/(-1 + yijk) 
     &   + (4*k1k3*p2k3*sdip*udip*zj**2)/(-1 + yijk) + (4*k1k3*M2**2*(k2k3*(-1 + yijk) + k1k3*yijk)*zj**2)/(-1 
     &   + yijk)**2 - (16*k1k3*M2*sdip*(k2k3*(-1 + yijk) + k1k3*yijk)*zj**2)/(-1 + yijk)**2 + (16*k1k3*p1p2*sdip
     &   *(k2k3*(-1 + yijk) + k1k3*yijk)*zj**2)/(-1 + yijk)**2 - (2*k1k3*M2*tdip*(k2k3*(-1 + yijk) + k1k3*yijk)
     &   *zj**2)/(-1 + yijk)**2 - (4*k1k3*tdip**2*(k2k3*(-1 + yijk) + k1k3*yijk)*zj**2)/(-1 + yijk)**2 - (6
     &   *k1k3*M2*udip*(k2k3*(-1 + yijk) + k1k3*yijk)*zj**2)/(-1 + yijk)**2 + (10*k1k3*tdip*udip*(k2k3*(-1 + yijk) 
     &   + k1k3*yijk)*zj**2)/(-1 + yijk)**2 - (2*k1k3*udip**2*(k2k3*(-1 + yijk) + k1k3*yijk)*zj**2)/(-1 
     &   + yijk)**2 + (16*k1k3**2*M2*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zj**2)/(-1 + yijk)**3 - (8*k1k3*
     &   *2*tdip*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zj**2)/(-1 + yijk)**3 - (8*k1k3**2*udip*((p2k2 + p2k3)
     &   *(-1 + yijk) + p2k1*yijk)*zj**2)/(-1 + yijk)**3 + (8*k1k3*p1k3*sdip*((p2k2 + p2k3)*(-1 + yijk) + p2k1
     &   *yijk)*zj**2)/(-1 + yijk)**2 - (8*k1k3*p2k3*sdip*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zj**2)/(-1 
     &   + yijk)**2 - (16*k1k3*M2*(k2k3*(-1 + yijk) + k1k3*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)*zj**2)/(
     &   -1 + yijk)**3 + (4*k1k3*tdip*(k2k3*(-1 + yijk) + k1k3*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)
     &   *zj**2)/(-1 + yijk)**3 + (12*k1k3*udip*(k2k3*(-1 + yijk) + k1k3*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1
     &   *yijk)*zj**2)/(-1 + yijk)**3 - (4*p1k3*sdip*(k2k3*(-1 + yijk) + k1k3*yijk)*((p2k2 + p2k3)*(-1 + yijk) 
     &   + p2k1*yijk)*zj**2)/(-1 + yijk)**2 + (4*p2k3*sdip*(k2k3*(-1 + yijk) + k1k3*yijk)*((p2k2 + p2k3)*(
     &   -1 + yijk) + p2k1*yijk)*zj**2)/(-1 + yijk)**2 + (8*tdip*(k2k3*(-1 + yijk) + k1k3*yijk)**2*((p2k2 + p2k3)
     &   *(-1 + yijk) + p2k1*yijk)*zj**2)/(-1 + yijk)**3 - (8*udip*(k2k3*(-1 + yijk) + k1k3*yijk)**2*((p2k2 + p2k3)
     &   *(-1 + yijk) + p2k1*yijk)*zj**2)/(-1 + yijk)**3 - (16*k1k3**2*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)**2*zj**2)/(-1 
     &   + yijk)**4 + (16*k1k3*(k2k3*(-1 + yijk) + k1k3*yijk)*((p2k2 + p2k3)*(-1 + yijk) + p2k1*yijk)**2*zj**2)/(-1 
     &   + yijk)**4 + 2*M2*p1k3*sdip*(k2k3 + (k1k3*yijk)/(-1 + yijk))*zj**2 - 2*M2*p2k3*sdip*(k2k3 + (k1k3*yijk)/(-1 
     &   + yijk))*zj**2 - 4*p1k3*sdip*tdip*(k2k3 + (k1k3*yijk)/(-1 + yijk))*zj**2 + 4*p2k3*sdip*tdip*(k2k3 + (k1k3*yijk)/(-1 
     &   + yijk))*zj**2 + 2*p1k3*sdip*udip*(k2k3 + (k1k3*yijk)/(-1 + yijk))*zj**2 - 2*p2k3*sdip*udip*(k2k3 + (k1k3*yijk)/(-1 
     &   + yijk))*zj**2 + 4*M2*sdip*(k2k3 + (k1k3*yijk)/(-1 + yijk))**2*zj**2 - 4*p1p2*sdip*(k2k3 + (k1k3*yijk)/(-1 + yijk))**2*zj**2 
     &   - 4*M2*tdip*(k2k3 + (k1k3*yijk)/(-1 + yijk))**2*zj**2 + 6*tdip**2*(k2k3 + (k1k3*yijk)/(-1 + yijk))**2*zj**2 + 4*M2*udip*(k2k3 
     &   + (k1k3*yijk)/(-1 + yijk))**2*zj**2 - 8*tdip*udip*(k2k3 + (k1k3*yijk)/(-1 + yijk))**2*zj**2 + 2*udip**2*(k2k3 + (k1k3*yijk)/(-1 
     &   + yijk))**2*zj**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(8*sdip**3))

!************************************************   2 Re(Ms VggP Mt^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_1 = stsT2xx_M2Dipole_SpinTensor_D23_1 + 2*dble(-0.125d0*(GGG*Nc
     &   *(-1 + Nc**2)*(-2*p1k2*(-1 + yijk)*zi + k1k2*yijk*zi + k2k3*(-1 + yijk)*(zi - zj) - 2*p1k3*zj - k1k3*yijk
     &   *zj + 2*p1k3*yijk*zj)*(M2*p2k2*sdip*zi - 2*p2k2**2*sdip*zi - 2*p2k2*p2k3*sdip*zi - 2*p2k2*sdip*tdip
     &   *zi + p2k2*sdip*udip*zi - 3*M2*p2k2*sdip*yijk*zi + 2*p2k1*p2k2*sdip*yijk*zi + 6*p2k2**2*sdip*yijk*zi 
     &   + 6*p2k2*p2k3*sdip*yijk*zi + 6*p2k2*sdip*tdip*yijk*zi - 3*p2k2*sdip*udip*yijk*zi + 3*M2*p2k2*sdip*yijk
     &   **2*zi - 4*p2k1*p2k2*sdip*yijk**2*zi - 6*p2k2**2*sdip*yijk**2*zi - 6*p2k2*p2k3*sdip*yijk**2*zi - 6
     &   *p2k2*sdip*tdip*yijk**2*zi + 3*p2k2*sdip*udip*yijk**2*zi - M2*p2k2*sdip*yijk**3*zi + 2*p2k1*p2k2*sdip
     &   *yijk**3*zi + 2*p2k2**2*sdip*yijk**3*zi + 2*p2k2*p2k3*sdip*yijk**3*zi + 2*p2k2*sdip*tdip*yijk**3*zi - p2k2
     &   *sdip*udip*yijk**3*zi + k1k2*(4*p2k3**2 + 4*p1p2*sdip - 4*p2k3*tdip - 2*tdip**2 + 2*p2k3*udip + 4
     &   *p2k2**2*(-1 + yijk)**2 + 2*M2**2*(-1 + yijk)**3 - 8*p2k1*p2k3*yijk - 8*p2k3**2*yijk - 10*p1p2*sdip*yijk 
     &   + 4*p2k1*tdip*yijk + 14*p2k3*tdip*yijk + 7*tdip**2*yijk - 2*p2k1*udip*yijk - 6*p2k3*udip*yijk - tdip
     &   *udip*yijk + 4*p2k1**2*yijk**2 + 8*p2k1*p2k3*yijk**2 + 4*p2k3**2*yijk**2 + 8*p1p2*sdip*yijk**2 - 10
     &   *p2k1*tdip*yijk**2 - 16*p2k3*tdip*yijk**2 - 8*tdip**2*yijk**2 + 4*p2k1*udip*yijk**2 + 6*p2k3*udip*yijk
     &   **2 + 2*tdip*udip*yijk**2 - 2*p1p2*sdip*yijk**3 + 6*p2k1*tdip*yijk**3 + 6*p2k3*tdip*yijk**3 + 3*tdip
     &   **2*yijk**3 - 2*p2k1*udip*yijk**3 - 2*p2k3*udip*yijk**3 - tdip*udip*yijk**3 + 2*p2k2*(-1 + yijk)*(-udip 
     &   + 4*p2k3*(-1 + yijk) + 4*p2k1*yijk + 2*udip*yijk - udip*yijk**2 + tdip*(2 - 5*yijk + 3*yijk**2)) 
     &   - M2*(-1 + yijk)*(-4*sdip + 4*tdip - 2*p2k1*yijk + 6*sdip*yijk - 9*tdip*yijk + udip*yijk + 4*p2k1*yijk
     &   **2 - 2*sdip*yijk**2 + 5*tdip*yijk**2 - udip*yijk**2 + p2k2*(2 - 6*yijk + 4*yijk**2) + p2k3*(2 - 6*yijk 
     &   + 4*yijk**2)))*zi + k2k3*(-1 + yijk)**2*(2*M2**2*(-1 + yijk) - 2*p1p2*sdip*(-1 + yijk) + (3*tdip 
     &   - udip)*(-tdip + 2*p2k2*(-1 + yijk) + 2*p2k3*(-1 + yijk) + 2*p2k1*yijk + tdip*yijk) - M2*(2*sdip - 5
     &   *tdip + udip + 4*p2k2*(-1 + yijk) + 4*p2k3*(-1 + yijk) + 4*p2k1*yijk - 2*sdip*yijk + 5*tdip*yijk - udip
     &   *yijk))*(zi - zj) + 2*k1k3*M2**2*zj - 2*k1k3*M2*p2k2*zj - 4*k1k3*p2k2**2*zj - 2*k1k3*M2*p2k3*zj - 8
     &   *k1k3*p2k2*p2k3*zj - 4*k1k3*p2k3**2*zj + 4*k1k3*M2*sdip*zj - 4*k1k3*p1p2*sdip*zj - M2*p2k3*sdip*zj + 2
     &   *p2k2*p2k3*sdip*zj + 2*p2k3**2*sdip*zj - 4*k1k3*M2*tdip*zj + 4*k1k3*p2k2*tdip*zj + 4*k1k3*p2k3*tdip*zj 
     &   + 2*p2k3*sdip*tdip*zj + 2*k1k3*tdip**2*zj - 2*k1k3*p2k2*udip*zj - 2*k1k3*p2k3*udip*zj - p2k3*sdip*udip
     &   *zj - 6*k1k3*M2**2*yijk*zj + 2*k1k3*M2*p2k1*yijk*zj + 8*k1k3*M2*p2k2*yijk*zj + 8*k1k3*p2k1*p2k2*yijk
     &   *zj + 8*k1k3*p2k2**2*yijk*zj + 8*k1k3*M2*p2k3*yijk*zj + 8*k1k3*p2k1*p2k3*yijk*zj + 16*k1k3*p2k2*p2k3
     &   *yijk*zj + 8*k1k3*p2k3**2*yijk*zj - 10*k1k3*M2*sdip*yijk*zj + 10*k1k3*p1p2*sdip*yijk*zj + 3*M2*p2k3
     &   *sdip*yijk*zj - 2*p2k1*p2k3*sdip*yijk*zj - 6*p2k2*p2k3*sdip*yijk*zj - 6*p2k3**2*sdip*yijk*zj + 13*k1k3
     &   *M2*tdip*yijk*zj - 4*k1k3*p2k1*tdip*yijk*zj - 14*k1k3*p2k2*tdip*yijk*zj - 14*k1k3*p2k3*tdip*yijk*zj 
     &   - 6*p2k3*sdip*tdip*yijk*zj - 7*k1k3*tdip**2*yijk*zj - k1k3*M2*udip*yijk*zj + 2*k1k3*p2k1*udip*yijk*zj 
     &   + 6*k1k3*p2k2*udip*yijk*zj + 6*k1k3*p2k3*udip*yijk*zj + 3*p2k3*sdip*udip*yijk*zj + k1k3*tdip*udip*yijk
     &   *zj + 6*k1k3*M2**2*yijk**2*zj - 6*k1k3*M2*p2k1*yijk**2*zj - 4*k1k3*p2k1**2*yijk**2*zj - 10*k1k3*M2*p2k2
     &   *yijk**2*zj - 8*k1k3*p2k1*p2k2*yijk**2*zj - 4*k1k3*p2k2**2*yijk**2*zj - 10*k1k3*M2*p2k3*yijk**2*zj 
     &   - 8*k1k3*p2k1*p2k3*yijk**2*zj - 8*k1k3*p2k2*p2k3*yijk**2*zj - 4*k1k3*p2k3**2*yijk**2*zj + 8*k1k3*M2
     &   *sdip*yijk**2*zj - 8*k1k3*p1p2*sdip*yijk**2*zj - 3*M2*p2k3*sdip*yijk**2*zj + 4*p2k1*p2k3*sdip*yijk**2
     &   *zj + 6*p2k2*p2k3*sdip*yijk**2*zj + 6*p2k3**2*sdip*yijk**2*zj - 14*k1k3*M2*tdip*yijk**2*zj + 10*k1k3*p2k1
     &   *tdip*yijk**2*zj + 16*k1k3*p2k2*tdip*yijk**2*zj + 16*k1k3*p2k3*tdip*yijk**2*zj + 6*p2k3*sdip*tdip
     &   *yijk**2*zj + 8*k1k3*tdip**2*yijk**2*zj + 2*k1k3*M2*udip*yijk**2*zj - 4*k1k3*p2k1*udip*yijk**2*zj - 6
     &   *k1k3*p2k2*udip*yijk**2*zj - 6*k1k3*p2k3*udip*yijk**2*zj - 3*p2k3*sdip*udip*yijk**2*zj - 2*k1k3*tdip*udip
     &   *yijk**2*zj - 2*k1k3*M2**2*yijk**3*zj + 4*k1k3*M2*p2k1*yijk**3*zj + 4*k1k3*M2*p2k2*yijk**3*zj + 4
     &   *k1k3*M2*p2k3*yijk**3*zj - 2*k1k3*M2*sdip*yijk**3*zj + 2*k1k3*p1p2*sdip*yijk**3*zj + M2*p2k3*sdip*yijk
     &   **3*zj - 2*p2k1*p2k3*sdip*yijk**3*zj - 2*p2k2*p2k3*sdip*yijk**3*zj - 2*p2k3**2*sdip*yijk**3*zj + 5*k1k3
     &   *M2*tdip*yijk**3*zj - 6*k1k3*p2k1*tdip*yijk**3*zj - 6*k1k3*p2k2*tdip*yijk**3*zj - 6*k1k3*p2k3*tdip*yijk**3
     &   *zj - 2*p2k3*sdip*tdip*yijk**3*zj - 3*k1k3*tdip**2*yijk**3*zj - k1k3*M2*udip*yijk**3*zj + 2*k1k3*p2k1*udip
     &   *yijk**3*zj + 2*k1k3*p2k2*udip*yijk**3*zj + 2*k1k3*p2k3*udip*yijk**3*zj + p2k3*sdip*udip*yijk**3*zj + k1k3*tdip*udip
     &   *yijk**3*zj)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(sdip**2*(M2 - tdip)*(-1 + yijk)**4))

!************************************************   2 Re(Ms VggP Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_1 = stsT2xx_M2Dipole_SpinTensor_D23_1 + 2*dble(-0.125d0*(GGG*Nc
     &   *(-1 + Nc**2)*(-2*p2k2*(-1 + yijk)*zi + k1k2*yijk*zi + k2k3*(-1 + yijk)*(zi - zj) - 2*p2k3*zj - k1k3*yijk
     &   *zj + 2*p2k3*yijk*zj)*(M2*p1k2*sdip*zi - 2*p1k2*p2k2*sdip*zi - 2*p1k2*p2k3*sdip*zi - 2*p1k2*sdip*tdip
     &   *zi + p1k2*sdip*udip*zi - 2*M2*p1k2*sdip*yijk*zi + 2*p1k2*p2k1*sdip*yijk*zi + 4*p1k2*p2k2*sdip*yijk
     &   *zi + 4*p1k2*p2k3*sdip*yijk*zi + 4*p1k2*sdip*tdip*yijk*zi - 2*p1k2*sdip*udip*yijk*zi + M2*p1k2*sdip
     &   *yijk**2*zi - 2*p1k2*p2k1*sdip*yijk**2*zi - 2*p1k2*p2k2*sdip*yijk**2*zi - 2*p1k2*p2k3*sdip*yijk**2*zi 
     &   - 2*p1k2*sdip*tdip*yijk**2*zi + p1k2*sdip*udip*yijk**2*zi + k1k2*(-2*p2k2*tdip - 2*p2k3*tdip - 4*p2k2
     &   *udip - 4*p2k3*udip - 3*tdip*udip + 2*p2k1*tdip*yijk + 2*p2k2*tdip*yijk + 2*p2k3*tdip*yijk - tdip**2
     &   *yijk + 4*p2k1*udip*yijk + 8*p2k2*udip*yijk + 8*p2k3*udip*yijk + 7*tdip*udip*yijk - udip**2*yijk + tdip
     &   **2*yijk**2 - 4*p2k1*udip*yijk**2 - 4*p2k2*udip*yijk**2 - 4*p2k3*udip*yijk**2 - 4*tdip*udip*yijk**2 
     &   + udip**2*yijk**2 + M2**2*(-3 + 5*yijk - 2*yijk**2) - 2*p1p2*sdip*(2 - 3*yijk + yijk**2) + M2*(4*sdip 
     &   + 3*tdip + 3*udip - 6*p2k1*yijk - 6*sdip*yijk - 5*tdip*yijk - 5*udip*yijk + 4*p2k1*yijk**2 + 2*sdip
     &   *yijk**2 + 2*tdip*yijk**2 + 2*udip*yijk**2 + 2*p2k2*(3 - 5*yijk + 2*yijk**2) + 2*p2k3*(3 - 5*yijk + 2
     &   *yijk**2)))*zi - k2k3*(-1 + yijk)*(tdip**2 - 4*p2k2*udip - 4*p2k3*udip - 4*tdip*udip + udip**2 + 2*M2
     &   **2*(-1 + yijk) + 2*p1p2*sdip*(-1 + yijk) - tdip**2*yijk + 4*p2k1*udip*yijk + 4*p2k2*udip*yijk + 4*p2k3
     &   *udip*yijk + 4*tdip*udip*yijk - udip**2*yijk - 2*M2*(-sdip - tdip - udip + 2*p2k2*(-1 + yijk) + 2*p2k3
     &   *(-1 + yijk) + 2*p2k1*yijk + sdip*yijk + tdip*yijk + udip*yijk))*(zi - zj) + 3*k1k3*M2**2*zj - 6*k1k3
     &   *M2*p2k2*zj - 6*k1k3*M2*p2k3*zj - 4*k1k3*M2*sdip*zj - M2*p1k3*sdip*zj + 4*k1k3*p1p2*sdip*zj + 2*p1k3
     &   *p2k2*sdip*zj + 2*p1k3*p2k3*sdip*zj - 3*k1k3*M2*tdip*zj + 2*k1k3*p2k2*tdip*zj + 2*k1k3*p2k3*tdip*zj 
     &   + 2*p1k3*sdip*tdip*zj - 3*k1k3*M2*udip*zj + 4*k1k3*p2k2*udip*zj + 4*k1k3*p2k3*udip*zj - p1k3*sdip*udip
     &   *zj + 3*k1k3*tdip*udip*zj - 5*k1k3*M2**2*yijk*zj + 6*k1k3*M2*p2k1*yijk*zj + 10*k1k3*M2*p2k2*yijk*zj 
     &   + 10*k1k3*M2*p2k3*yijk*zj + 6*k1k3*M2*sdip*yijk*zj + 2*M2*p1k3*sdip*yijk*zj - 6*k1k3*p1p2*sdip*yijk*zj 
     &   - 2*p1k3*p2k1*sdip*yijk*zj - 4*p1k3*p2k2*sdip*yijk*zj - 4*p1k3*p2k3*sdip*yijk*zj + 5*k1k3*M2*tdip*yijk
     &   *zj - 2*k1k3*p2k1*tdip*yijk*zj - 2*k1k3*p2k2*tdip*yijk*zj - 2*k1k3*p2k3*tdip*yijk*zj - 4*p1k3*sdip
     &   *tdip*yijk*zj + k1k3*tdip**2*yijk*zj + 5*k1k3*M2*udip*yijk*zj - 4*k1k3*p2k1*udip*yijk*zj - 8*k1k3*p2k2
     &   *udip*yijk*zj - 8*k1k3*p2k3*udip*yijk*zj + 2*p1k3*sdip*udip*yijk*zj - 7*k1k3*tdip*udip*yijk*zj + k1k3
     &   *udip**2*yijk*zj + 2*k1k3*M2**2*yijk**2*zj - 4*k1k3*M2*p2k1*yijk**2*zj - 4*k1k3*M2*p2k2*yijk**2*zj 
     &   - 4*k1k3*M2*p2k3*yijk**2*zj - 2*k1k3*M2*sdip*yijk**2*zj - M2*p1k3*sdip*yijk**2*zj + 2*k1k3*p1p2*sdip*yijk
     &   **2*zj + 2*p1k3*p2k1*sdip*yijk**2*zj + 2*p1k3*p2k2*sdip*yijk**2*zj + 2*p1k3*p2k3*sdip*yijk**2*zj 
     &   - 2*k1k3*M2*tdip*yijk**2*zj + 2*p1k3*sdip*tdip*yijk**2*zj - k1k3*tdip**2*yijk**2*zj - 2*k1k3*M2*udip*yijk**2*zj 
     &  + 4*k1k3*p2k1*udip*yijk**2*zj + 4*k1k3*p2k2*udip*yijk**2*zj + 4*k1k3*p2k3*udip*yijk**2*zj - p1k3*sdip*udip*yijk**2*zj 
     &  + 4*k1k3*tdip*udip*yijk**2*zj - k1k3*udip**2*yijk**2*zj)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(sdip**2*(M2 - udip)*(-1 + yijk)**3))

!************************************************   2 Re(Ms VggP Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_1 = stsT2xx_M2Dipole_SpinTensor_D23_1 + 2*dble(0)

!************************************************   Mt^ast VggP Mt   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_1 = stsT2xx_M2Dipole_SpinTensor_D23_1 +((Nc*(-1 + Nc**2)*(M2*sdip
     &   *(-1 + yijk) - 2*(p2k2 + p2k3)*(M2 - tdip)*(-1 + yijk) + 2*p2k1*(-M2 + tdip)*yijk)*(-2*p1k2*(-1 + yijk)
     &   *zi + k1k2*yijk*zi + k2k3*(-1 + yijk)*(zi - zj) - 2*p1k3*zj - k1k3*yijk*zj + 2*p1k3*yijk*zj)**2
     &   *SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(2*sdip*(M2 - tdip)**2*(-1 + yijk)**3))

!************************************************   2 Re(Mt VggP Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_1 = stsT2xx_M2Dipole_SpinTensor_D23_1 + 2*dble(0)

!************************************************   2 Re(Mt VggP Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_1 = stsT2xx_M2Dipole_SpinTensor_D23_1 + 2*dble((Nc*(-1 + Nc**2)
     &   *(-2*p1k2*(-1 + yijk)*zi + k1k2*yijk*zi + k2k3*(-1 + yijk)*(zi - zj) - 2*p1k3*zj - k1k3*yijk*zj + 2
     &   *p1k3*yijk*zj)*(p2k2*sdip*zi - 2*p2k2*sdip*yijk*zi + p2k2*sdip*yijk**2*zi + k1k2*(2*p2k2*(-1 + yijk) 
     &   + 2*p2k3*(-1 + yijk) + yijk*(M2 + 2*p2k1 + tdip*(-1 + yijk) - M2*yijk))*zi - k2k3*(M2 - tdip)*(-1 + yijk)
     &   **2*(zi - zj) + 2*k1k3*p2k2*zj + 2*k1k3*p2k3*zj - p2k3*sdip*zj - k1k3*M2*yijk*zj - 2*k1k3*p2k1*yijk
     &   *zj - 2*k1k3*p2k2*yijk*zj - 2*k1k3*p2k3*yijk*zj + 2*p2k3*sdip*yijk*zj + k1k3*tdip*yijk*zj + k1k3*M2
     &   *yijk**2*zj - p2k3*sdip*yijk**2*zj - k1k3*tdip*yijk**2*zj)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(4.
     &   *sdip*(M2 - tdip)*(-1 + yijk)**3))

!************************************************   Mu^ast VggP Mu   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_1 = stsT2xx_M2Dipole_SpinTensor_D23_1 +((Nc*(-1 + Nc**2)*(-M2
     &   **2 - tdip*udip + M2*(sdip + tdip + udip))*(-2*p2k2*(-1 + yijk)*zi + k1k2*yijk*zi + k2k3*(-1 + yijk)*(zi 
     &   - zj) - 2*p2k3*zj - k1k3*yijk*zj + 2*p2k3*yijk*zj)**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(2
     &   *sdip*(M2 - udip)**2*(-1 + yijk)**2))

!************************************************   2 Re(Mu VggP Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_1 = stsT2xx_M2Dipole_SpinTensor_D23_1 + 2*dble(-0.25d0*(Nc*(-1 
     &   + Nc**2)*(-2*p2k2*(-1 + yijk)*zi + k1k2*yijk*zi + k2k3*(-1 + yijk)*(zi - zj) - 2*p2k3*zj - k1k3*yijk*zj 
     &   + 2*p2k3*yijk*zj)*(p1k2*sdip*zi - p1k2*sdip*yijk*zi + k1k2*(tdip + M2*(-1 + yijk) - udip*yijk)*zi 
     &   + k2k3*(M2 - udip)*(-1 + yijk)*(zi - zj) + k1k3*M2*zj - p1k3*sdip*zj - k1k3*tdip*zj - k1k3*M2*yijk*zj 
     &   + p1k3*sdip*yijk*zj + k1k3*udip*yijk*zj)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(sdip
     &   *(M2 - udip)*(-1 + yijk)**2))

!************************************************   Mv^ast VggP Mv   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_1 = stsT2xx_M2Dipole_SpinTensor_D23_1 +((Nc*(-1 + Nc**2)*(2*k1k2
     &   **2*yijk*zi**2 + 2*k1k2*zi*(k2k3*(-1 + yijk)*(zi - zj) - 2*k1k3*yijk*zj) + zj*(-(k2k3*sdip*(-1 + yijk)
     &   **2*zi) - 2*k1k3*k2k3*(-1 + yijk)*(zi - zj) + 2*k1k3**2*yijk*zj))*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **2)/(2*sdip*(-1 + yijk)**2))



      end


      double precision function stsT2xx_M2Dipole_SpinTensor_D12_a(sdip,tdip,udip)
      implicit none
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"
#include "common/stsT2xx_GenCouplings.h"

      double precision sdip,tdip,udip,M2,zi,zM,zim,zjm
      double precision xija,Rdip,kappa,ab1,ab2

      M2=MSf2(iSf1,itt1,igen1)

      xija = (-k1k2 + p1k1 + p1k2)/(p1k1 + p1k2)

      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)

      zi =  p1k1/(p1k1 + p1k2)

      zM = (1-Rdip)/2d0

      ab1 = 1 - (M2 + tdip)/(2*tdip) + (tdip - M2 + 2*M2*xija)/(2*Rdip*tdip)
      ab2 = xija/Rdip

      zim = zi - zM
      zjm = 1 - zi - zM

!************************************************   Ms^ast Vgg Ms   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_a = ((GGG**2*Nc**2*(-1 + Nc
     &   **2)*(28*k1k3**2*M2*sdip*zim**2 + 8*k1k3*M2*p2k1*sdip*zim**2 + 8*k1k3**2*p2k3*sdip*zim**2 - 16*k1k3*p2k1
     &   *p2k3*sdip*zim**2 - 8*k1k3**2*sdip**2*zim**2 + 4*k1k3**2*M2*tdip*zim**2 - 8*k1k3**2*p2k3*tdip*zim
     &   **2 - 4*k1k3**2*sdip*tdip*zim**2 + 8*k1k3*p2k1*sdip*tdip*zim**2 - 4*k1k3**2*tdip**2*zim**2 - 4*k1k3**2
     &   *M2*udip*zim**2 + 8*k1k3**2*p2k3*udip*zim**2 + 8*k1k3**2*sdip*udip*zim**2 - 16*k1k3*p2k1*sdip*udip*zim
     &   **2 + 4*k1k3**2*tdip*udip*zim**2 - 56*k1k3*k2k3*M2*sdip*zim*zjm - k1k2*M2**2*sdip*zim*zjm - 8*k2k3*M2
     &   *p2k1*sdip*zim*zjm - 8*k1k3*M2*p2k2*sdip*zim*zjm - 16*k1k3*k2k3*p2k3*sdip*zim*zjm + 4*k1k2*M2*p2k3*sdip
     &   *zim*zjm + 16*k2k3*p2k1*p2k3*sdip*zim*zjm + 16*k1k3*p2k2*p2k3*sdip*zim*zjm - 4*k1k2*p2k3**2*sdip*zim
     &   *zjm + 16*k1k3*k2k3*sdip**2*zim*zjm - 8*k1k3*k2k3*M2*tdip*zim*zjm + 16*k1k3*k2k3*p2k3*tdip*zim*zjm 
     &   + 8*k1k3*k2k3*sdip*tdip*zim*zjm - 2*k1k2*M2*sdip*tdip*zim*zjm - 8*k2k3*p2k1*sdip*tdip*zim*zjm - 8*k1k3
     &   *p2k2*sdip*tdip*zim*zjm + 4*k1k2*p2k3*sdip*tdip*zim*zjm + 8*k1k3*k2k3*tdip**2*zim*zjm - k1k2*sdip*tdip
     &   **2*zim*zjm + 8*k1k3*k2k3*M2*udip*zim*zjm - 16*k1k3*k2k3*p2k3*udip*zim*zjm - 16*k1k3*k2k3*sdip*udip
     &   *zim*zjm + 4*k1k2*M2*sdip*udip*zim*zjm + 16*k2k3*p2k1*sdip*udip*zim*zjm + 16*k1k3*p2k2*sdip*udip*zim
     &   *zjm - 8*k1k2*p2k3*sdip*udip*zim*zjm - 8*k1k3*k2k3*tdip*udip*zim*zjm + 4*k1k2*sdip*tdip*udip*zim*zjm 
     &   - 4*k1k2*sdip*udip**2*zim*zjm + 28*k2k3**2*M2*sdip*zjm**2 + 8*k2k3*M2*p2k2*sdip*zjm**2 + 8*k2k3**2*p2k3
     &   *sdip*zjm**2 - 16*k2k3*p2k2*p2k3*sdip*zjm**2 - 8*k2k3**2*sdip**2*zjm**2 + 4*k2k3**2*M2*tdip*zjm**2 
     &   - 8*k2k3**2*p2k3*tdip*zjm**2 - 4*k2k3**2*sdip*tdip*zjm**2 + 8*k2k3*p2k2*sdip*tdip*zjm**2 - 4*k2k3**2*tdip
     &   **2*zjm**2 - 4*k2k3**2*M2*udip*zjm**2 + 8*k2k3**2*p2k3*udip*zjm**2 + 8*k2k3**2*sdip*udip*zjm**2 - 16
     &   *k2k3*p2k2*sdip*udip*zjm**2 + 4*k2k3**2*tdip*udip*zjm**2 + 2*ab2**2*(M2**2 + 4*p2k3**2 - sdip**2 - sdip
     &   *tdip + M2*(-4*p2k3 + 3*sdip + 2*tdip - 4*udip) + 2*sdip*udip - 2*tdip*udip + 3*udip**2 + 2*p2k3
     &   *(sdip - 2*tdip + 4*udip))*(p1k1*zim - p1k2*zjm)**2 + 2*ab1**2*(M2**2 + 4*p2k3**2 - sdip**2 - sdip*tdip 
     &   + M2*(-4*p2k3 + 3*sdip + 2*tdip - 4*udip) + 2*sdip*udip - 2*tdip*udip + 3*udip**2 + 2*p2k3*(sdip - 2
     &   *tdip + 4*udip))*(k1k3*zim - p2k1*zim + (-k2k3 + p2k2)*zjm)**2 - 2*ab1*(k1k3*zim - p2k1*zim + (-k2k3 
     &   + p2k2)*zjm)*(2*M2*p2k1*sdip*zim - 4*p2k1*p2k3*sdip*zim + 2*p2k1*sdip*tdip*zim - 4*p2k1*sdip*udip*zim 
     &   + k1k3*(M2**2 + 4*p2k3**2 - 4*sdip**2 - 3*sdip*tdip - 2*tdip**2 + M2*(-4*p2k3 + 13*sdip + 3*tdip - 5
     &   *udip) + 6*sdip*udip + tdip*udip + 2*udip**2 + 2*p2k3*(3*sdip - 3*tdip + 5*udip))*zim - k2k3*M2**2*zjm 
     &   + 4*k2k3*M2*p2k3*zjm - 4*k2k3*p2k3**2*zjm - 13*k2k3*M2*sdip*zjm - 2*M2*p2k2*sdip*zjm - 6*k2k3*p2k3
     &   *sdip*zjm + 4*p2k2*p2k3*sdip*zjm + 4*k2k3*sdip**2*zjm - 3*k2k3*M2*tdip*zjm + 6*k2k3*p2k3*tdip*zjm + 3
     &   *k2k3*sdip*tdip*zjm - 2*p2k2*sdip*tdip*zjm + 2*k2k3*tdip**2*zjm + 5*k2k3*M2*udip*zjm - 10*k2k3*p2k3*udip
     &   *zjm - 6*k2k3*sdip*udip*zjm + 4*p2k2*sdip*udip*zjm - k2k3*tdip*udip*zjm - 2*k2k3*udip**2*zjm + 2*ab2
     &   *(M2**2 + 4*p2k3**2 - sdip**2 - sdip*tdip + M2*(-4*p2k3 + 3*sdip + 2*tdip - 4*udip) + 2*sdip*udip 
     &   - 2*tdip*udip + 3*udip**2 + 2*p2k3*(sdip - 2*tdip + 4*udip))*(p1k1*zim - p1k2*zjm)) + 2*ab2*(p1k1*zim 
     &   - p1k2*zjm)*(-4*p2k1*p2k3*sdip*zim + 2*p2k1*sdip*tdip*zim - 4*p2k1*sdip*udip*zim + k1k3*(M2**2 + 4*p2k3
     &   **2 + 6*p2k3*sdip - 4*sdip**2 - 6*p2k3*tdip - 3*sdip*tdip - 2*tdip**2 + M2*(-4*p2k3 + 13*sdip + 3*tdip 
     &   - 5*udip) + 10*p2k3*udip + 6*sdip*udip + tdip*udip + 2*udip**2)*zim - k2k3*M2**2*zjm - 4*k2k3*p2k3
     &   **2*zjm - 6*k2k3*p2k3*sdip*zjm + 4*p2k2*p2k3*sdip*zjm + 4*k2k3*sdip**2*zjm + 6*k2k3*p2k3*tdip*zjm 
     &   + 3*k2k3*sdip*tdip*zjm - 2*p2k2*sdip*tdip*zjm + 2*k2k3*tdip**2*zjm - 10*k2k3*p2k3*udip*zjm - 6*k2k3
     &   *sdip*udip*zjm + 4*p2k2*sdip*udip*zjm - k2k3*tdip*udip*zjm - 2*k2k3*udip**2*zjm + M2*(2*p2k1*sdip*zim 
     &   - 2*p2k2*sdip*zjm + k2k3*(4*p2k3 - 13*sdip - 3*tdip + 5*udip)*zjm)))
     &   *SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(16*CV*sdip**3))

!************************************************   2 Re(Ms Vgg Mt^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_a = stsT2xx_M2Dipole_SpinTensor_D12_a + 2*dble((GGG*Nc**2*(-1 
     &   + Nc**2)*(-(ab1*M2**2*p2k1*zim) + 4*ab1*M2*p2k1*p2k3*zim - 4*ab1*p2k1*p2k3**2*zim - M2*p2k1*sdip*zim 
     &   - 4*ab1*M2*p2k1*sdip*zim + 2*p2k1*p2k3*sdip*zim + ab1*p2k1*sdip**2*zim + 4*ab1*p2k1*p2k3*tdip*zim - p2k1
     &   *sdip*tdip*zim + 2*ab1*M2*p2k1*udip*zim - 8*ab1*p2k1*p2k3*udip*zim + 2*p2k1*sdip*udip*zim - ab1*p2k1
     &   *udip**2*zim + k1k3*((-1 + ab1)*M2**2 + 2*sdip**2 + 4*p2k3*tdip - 6*p2k3*udip - tdip*udip + M2*((2 
     &   - 4*ab1)*p2k3 + 4*(-2 + ab1)*sdip + tdip + udip - 2*ab1*udip) + ab1*(4*p2k3**2 - sdip**2 - 4*p2k3*tdip 
     &   + 8*p2k3*udip + udip**2))*zim + k2k3*M2**2*zjm - ab1*k2k3*M2**2*zjm + ab1*M2**2*p2k2*zjm - 2*k2k3*M2
     &   *p2k3*zjm + 4*ab1*k2k3*M2*p2k3*zjm - 4*ab1*M2*p2k2*p2k3*zjm - 4*ab1*k2k3*p2k3**2*zjm + 4*ab1*p2k2*p2k3
     &   **2*zjm + 8*k2k3*M2*sdip*zjm - 4*ab1*k2k3*M2*sdip*zjm + M2*p2k2*sdip*zjm + 4*ab1*M2*p2k2*sdip*zjm 
     &   - 2*p2k2*p2k3*sdip*zjm - 2*k2k3*sdip**2*zjm + ab1*k2k3*sdip**2*zjm - ab1*p2k2*sdip**2*zjm - k2k3*M2*tdip
     &   *zjm - 4*k2k3*p2k3*tdip*zjm + 4*ab1*k2k3*p2k3*tdip*zjm - 4*ab1*p2k2*p2k3*tdip*zjm + p2k2*sdip*tdip
     &   *zjm - k2k3*M2*udip*zjm + 2*ab1*k2k3*M2*udip*zjm - 2*ab1*M2*p2k2*udip*zjm + 6*k2k3*p2k3*udip*zjm - 8*ab1
     &   *k2k3*p2k3*udip*zjm + 8*ab1*p2k2*p2k3*udip*zjm - 2*p2k2*sdip*udip*zjm + k2k3*tdip*udip*zjm - ab1*k2k3
     &   *udip**2*zjm + ab1*p2k2*udip**2*zjm - ab2*(M2**2 + 4*p2k3**2 - sdip**2 - 4*p2k3*tdip + 8*p2k3*udip 
     &   + udip**2 - 2*M2*(2*p2k3 - 2*sdip + udip))*(p1k1*zim - p1k2*zjm))*((-2 + ab1)*k1k3*zim + ab2*(-(p1k1*zim) 
     &   + p1k2*zjm) - (-2 + ab1)*(p2k1*zim + k2k3*zjm - p2k2*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(8*CV*sdip**2*(M2 - tdip)))

!************************************************   2 Re(Ms Vgg Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_a = stsT2xx_M2Dipole_SpinTensor_D12_a + 2*dble(0)

!************************************************   2 Re(Ms Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_a = stsT2xx_M2Dipole_SpinTensor_D12_a + 2*dble((GGG*Nc**2*(-1 
     &   + Nc**2)*(2*k1k3**2*sdip*zim**2 - 4*k1k3*p2k1*sdip*zim**2 + 2*k1k3**2*tdip*zim**2 - 2*k1k3**2*udip*zim
     &   **2 - 4*k1k3*k2k3*sdip*zim*zjm + k1k2*M2*sdip*zim*zjm + 4*k2k3*p2k1*sdip*zim*zjm + 4*k1k3*p2k2*sdip
     &   *zim*zjm - 2*k1k2*p2k3*sdip*zim*zjm - 4*k1k3*k2k3*tdip*zim*zjm + k1k2*sdip*tdip*zim*zjm + 4*k1k3*k2k3
     &   *udip*zim*zjm - 2*k1k2*sdip*udip*zim*zjm + 2*k2k3**2*sdip*zjm**2 - 4*k2k3*p2k2*sdip*zjm**2 + 2*k2k3**2
     &   *tdip*zjm**2 - 2*k2k3**2*udip*zjm**2 - ab2**2*(M2 - 2*p2k3 - sdip - udip)*(p1k1*zim - p1k2*zjm)**2 - ab1
     &   **2*(M2 - 2*p2k3 - sdip - udip)*(k1k3*zim - p2k1*zim + (-k2k3 + p2k2)*zjm)**2 + ab2*(p1k1*zim - p1k2
     &   *zjm)*(-2*p2k1*sdip*zim + 3*k1k3*(sdip + tdip - udip)*zim + 2*p2k2*sdip*zjm - 3*k2k3*(sdip + tdip 
     &   - udip)*zjm) + ab1*(k1k3*zim - p2k1*zim + (-k2k3 + p2k2)*zjm)*(2*p2k1*sdip*zim - 3*k1k3*(sdip + tdip 
     &   - udip)*zim + 3*k2k3*sdip*zjm - 2*p2k2*sdip*zjm + 3*k2k3*tdip*zjm - 3*k2k3*udip*zjm + 2*ab2*(M2 - 2*p2k3 
     &   - sdip - udip)*(p1k1*zim - p1k2*zjm)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(8*CV*sdip**2))

!************************************************   Mt^ast Vgg Mt   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_a = stsT2xx_M2Dipole_SpinTensor_D12_a +(((-1 + Nc**2)**2*(-8*M2
     &   *p2k3 + 4*M2*sdip + 8*p2k3*udip)*(-(ab2*p1k1*zim) + ab2*p1k2*zjm + (-2 + ab1)*(k1k3*zim - p2k1*zim 
     &   + (-k2k3 + p2k2)*zjm))**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(8*CV*sdip*(M2 - tdip)**2))

!************************************************   2 Re(Mt Vgg Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_a = stsT2xx_M2Dipole_SpinTensor_D12_a + 2*dble(((-1 + Nc**2)*(M2
     &   **2 - sdip**2 - 2*p2k3*tdip + 2*M2*(p2k3 + sdip - udip) + udip**2)*(-(ab2*p1k1*zim) + 2*p2k1*zim + ab2
     &   *p1k2*zjm - 2*p2k2*zjm + ab1*(k1k3*zim - p2k1*zim - k2k3*zjm + p2k2*zjm))*(ab2*p1k1*zim - ab2*p1k2
     &   *zjm - (-2 + ab1)*(k1k3*zim - p2k1*zim + (-k2k3 + p2k2)*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **4)/(4*CV*sdip*(M2 - tdip)*(M2 - udip)))

!************************************************   2 Re(Mt Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_a = stsT2xx_M2Dipole_SpinTensor_D12_a + 2*dble(((2 - 3*Nc**2 
     &   + Nc**4)*(2*ab2*p1k1*p2k3*zim + 2*ab1*(-k1k3 + p2k1)*p2k3*zim - p2k1*sdip*zim + k1k3*(M2 - udip)*zim 
     &   - 2*ab2*p1k2*p2k3*zjm + 2*ab1*(k2k3 - p2k2)*p2k3*zjm + p2k2*sdip*zjm + k2k3*(-M2 + udip)*zjm)*(ab2*p1k1
     &   *zim - ab2*p1k2*zjm - (-2 + ab1)*(k1k3*zim - p2k1*zim + (-k2k3 + p2k2)*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(4*CV*sdip*(M2 - tdip)))

!************************************************   Mu^ast Vgg Mu   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_a = stsT2xx_M2Dipole_SpinTensor_D12_a +(((-1 + Nc**2)*(M2**2 
     &   + tdip*udip - M2*(sdip + tdip + udip))*(-(ab2*p1k1*zim) + 2*p2k1*zim + ab2*p1k2*zjm - 2*p2k2*zjm + ab1
     &   *(k1k3*zim - p2k1*zim - k2k3*zjm + p2k2*zjm))**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(2*CV*sdip*(M2 - udip)**2))

!************************************************   2 Re(Mu Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_a = stsT2xx_M2Dipole_SpinTensor_D12_a + 2*dble(((-1 + Nc**2)*(
     &   -(ab1*M2*p2k1*zim) - p2k1*sdip*zim + ab1*p2k1*sdip*zim + ab1*p2k1*udip*zim + k1k3*((-1 + ab1)*M2 + sdip 
     &   - ab1*sdip + tdip - ab1*udip)*zim + k2k3*M2*zjm - ab1*k2k3*M2*zjm + ab1*M2*p2k2*zjm - k2k3*sdip*zjm 
     &   + ab1*k2k3*sdip*zjm + p2k2*sdip*zjm - ab1*p2k2*sdip*zjm - k2k3*tdip*zjm + ab1*k2k3*udip*zjm - ab1*p2k2
     &   *udip*zjm - ab2*(M2 - sdip - udip)*(p1k1*zim - p1k2*zjm))*(-(ab2*p1k1*zim) + 2*p2k1*zim + ab2*p1k2
     &   *zjm - 2*p2k2*zjm + ab1*(k1k3*zim - p2k1*zim - k2k3*zjm + p2k2*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(2*CV*sdip*(M2 - udip)))

!************************************************   Mv^ast Vgg Mv   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_a = stsT2xx_M2Dipole_SpinTensor_D12_a +(((4 - 5*Nc**2 + Nc**4)
     &   *(-(k1k2*sdip*zim*zjm) - 2*ab2*(k1k3*zim - k2k3*zjm)*(p1k1*zim - p1k2*zjm) + 2*ab1*(k1k3*zim - k2k3*zjm)
     &   *(k1k3*zim - p2k1*zim - k2k3*zjm + p2k2*zjm))*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(4*CV*sdip))


      end


      double precision function stsT2xx_M2Dipole_SpinTensor_D23_a(sdip,tdip,udip)
      implicit none
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"
#include "common/stsT2xx_GenCouplings.h"

      double precision sdip,tdip,udip,M2,zi,zM,zim,zjm
      double precision xija,Rdip,kappa,ab1,ab2

      M2=MSf2(iSf1,itt1,igen1)

      xija = (-k2k3 + p1k2 + p1k3)/(p1k2 + p1k3)

      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)

      zi =  p1k2/(p1k2 + p1k3)

      zM = (1-Rdip)/2d0

      ab1 = 1 - (M2 + tdip)/(2*tdip) + (tdip - M2 + 2*M2*xija)/(2*Rdip*tdip)
      ab2 = xija/Rdip

      zim = zi - zM
      zjm = 1 - zi - zM

!************************************************   Ms^ast Vgg Ms   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_a = ((GGG**2*Nc*(-1 + Nc**2)
     &   *(28*k1k2**2*M2*sdip*zim**2 + 8*k1k2**2*p2k1*sdip*zim**2 + 8*k1k2*M2*p2k2*sdip*zim**2 - 16*k1k2*p2k1
     &   *p2k2*sdip*zim**2 - 8*k1k2**2*sdip**2*zim**2 + 4*k1k2**2*M2*tdip*zim**2 - 8*k1k2**2*p2k1*tdip*zim**2 
     &   - 4*k1k2**2*sdip*tdip*zim**2 + 8*k1k2*p2k2*sdip*tdip*zim**2 - 4*k1k2**2*tdip**2*zim**2 - 4*k1k2**2*M2
     &   *udip*zim**2 + 8*k1k2**2*p2k1*udip*zim**2 + 8*k1k2**2*sdip*udip*zim**2 - 16*k1k2*p2k2*sdip*udip*zim*
     &   *2 + 4*k1k2**2*tdip*udip*zim**2 - 56*k1k2*k1k3*M2*sdip*zim*zjm - k2k3*M2**2*sdip*zim*zjm - 16*k1k2*k1k3
     &   *p2k1*sdip*zim*zjm + 4*k2k3*M2*p2k1*sdip*zim*zjm - 4*k2k3*p2k1**2*sdip*zim*zjm - 8*k1k3*M2*p2k2*sdip
     &   *zim*zjm + 16*k1k3*p2k1*p2k2*sdip*zim*zjm - 8*k1k2*M2*p2k3*sdip*zim*zjm + 16*k1k2*p2k1*p2k3*sdip*zim
     &   *zjm + 16*k1k2*k1k3*sdip**2*zim*zjm - 8*k1k2*k1k3*M2*tdip*zim*zjm + 16*k1k2*k1k3*p2k1*tdip*zim*zjm + 8
     &   *k1k2*k1k3*sdip*tdip*zim*zjm - 2*k2k3*M2*sdip*tdip*zim*zjm + 4*k2k3*p2k1*sdip*tdip*zim*zjm - 8*k1k3*p2k2
     &   *sdip*tdip*zim*zjm - 8*k1k2*p2k3*sdip*tdip*zim*zjm + 8*k1k2*k1k3*tdip**2*zim*zjm - k2k3*sdip*tdip
     &   **2*zim*zjm + 8*k1k2*k1k3*M2*udip*zim*zjm - 16*k1k2*k1k3*p2k1*udip*zim*zjm - 16*k1k2*k1k3*sdip*udip*zim
     &   *zjm + 4*k2k3*M2*sdip*udip*zim*zjm - 8*k2k3*p2k1*sdip*udip*zim*zjm + 16*k1k3*p2k2*sdip*udip*zim*zjm 
     &   + 16*k1k2*p2k3*sdip*udip*zim*zjm - 8*k1k2*k1k3*tdip*udip*zim*zjm + 4*k2k3*sdip*tdip*udip*zim*zjm - 4
     &   *k2k3*sdip*udip**2*zim*zjm + 28*k1k3**2*M2*sdip*zjm**2 + 8*k1k3**2*p2k1*sdip*zjm**2 + 8*k1k3*M2*p2k3*sdip
     &   *zjm**2 - 16*k1k3*p2k1*p2k3*sdip*zjm**2 - 8*k1k3**2*sdip**2*zjm**2 + 4*k1k3**2*M2*tdip*zjm**2 - 8
     &   *k1k3**2*p2k1*tdip*zjm**2 - 4*k1k3**2*sdip*tdip*zjm**2 + 8*k1k3*p2k3*sdip*tdip*zjm**2 - 4*k1k3**2*tdip
     &   **2*zjm**2 - 4*k1k3**2*M2*udip*zjm**2 + 8*k1k3**2*p2k1*udip*zjm**2 + 8*k1k3**2*sdip*udip*zjm**2 - 16
     &   *k1k3*p2k3*sdip*udip*zjm**2 + 4*k1k3**2*tdip*udip*zjm**2 + 2*ab2**2*(M2**2 + 4*p2k1**2 - sdip**2 - sdip
     &   *tdip + M2*(-4*p2k1 + 3*sdip + 2*tdip - 4*udip) + 2*sdip*udip - 2*tdip*udip + 3*udip**2 + 2*p2k1*(sdip 
     &   - 2*tdip + 4*udip))*(p1k2*zim - p1k3*zjm)**2 + 2*ab1**2*(M2**2 + 4*p2k1**2 - sdip**2 - sdip*tdip 
     &   + M2*(-4*p2k1 + 3*sdip + 2*tdip - 4*udip) + 2*sdip*udip - 2*tdip*udip + 3*udip**2 + 2*p2k1*(sdip - 2*tdip 
     &   + 4*udip))*(k1k2*zim - p2k2*zim + (-k1k3 + p2k3)*zjm)**2 - 2*ab1*(k1k2*zim - p2k2*zim + (-k1k3 + p2k3)
     &   *zjm)*(2*M2*p2k2*sdip*zim - 4*p2k1*p2k2*sdip*zim + 2*p2k2*sdip*tdip*zim - 4*p2k2*sdip*udip*zim 
     &   + k1k2*(M2**2 + 4*p2k1**2 - 4*sdip**2 - 3*sdip*tdip - 2*tdip**2 + M2*(-4*p2k1 + 13*sdip + 3*tdip - 5*udip) 
     &   + 6*sdip*udip + tdip*udip + 2*udip**2 + 2*p2k1*(3*sdip - 3*tdip + 5*udip))*zim - k1k3*M2**2*zjm 
     &   + 4*k1k3*M2*p2k1*zjm - 4*k1k3*p2k1**2*zjm - 13*k1k3*M2*sdip*zjm - 6*k1k3*p2k1*sdip*zjm - 2*M2*p2k3*sdip
     &   *zjm + 4*p2k1*p2k3*sdip*zjm + 4*k1k3*sdip**2*zjm - 3*k1k3*M2*tdip*zjm + 6*k1k3*p2k1*tdip*zjm + 3*k1k3
     &   *sdip*tdip*zjm - 2*p2k3*sdip*tdip*zjm + 2*k1k3*tdip**2*zjm + 5*k1k3*M2*udip*zjm - 10*k1k3*p2k1*udip
     &   *zjm - 6*k1k3*sdip*udip*zjm + 4*p2k3*sdip*udip*zjm - k1k3*tdip*udip*zjm - 2*k1k3*udip**2*zjm + 2*ab2
     &   *(M2**2 + 4*p2k1**2 - sdip**2 - sdip*tdip + M2*(-4*p2k1 + 3*sdip + 2*tdip - 4*udip) + 2*sdip*udip - 2
     &   *tdip*udip + 3*udip**2 + 2*p2k1*(sdip - 2*tdip + 4*udip))*(p1k2*zim - p1k3*zjm)) + 2*ab2*(p1k2*zim - p1k3
     &   *zjm)*(-4*p2k1*p2k2*sdip*zim + 2*p2k2*sdip*tdip*zim - 4*p2k2*sdip*udip*zim + k1k2*(M2**2 + 4*p2k1
     &   **2 + 6*p2k1*sdip - 4*sdip**2 - 6*p2k1*tdip - 3*sdip*tdip - 2*tdip**2 + M2*(-4*p2k1 + 13*sdip + 3*tdip 
     &   - 5*udip) + 10*p2k1*udip + 6*sdip*udip + tdip*udip + 2*udip**2)*zim - k1k3*M2**2*zjm - 4*k1k3*p2k1*
     &   *2*zjm - 6*k1k3*p2k1*sdip*zjm + 4*p2k1*p2k3*sdip*zjm + 4*k1k3*sdip**2*zjm + 6*k1k3*p2k1*tdip*zjm 
     &   + 3*k1k3*sdip*tdip*zjm - 2*p2k3*sdip*tdip*zjm + 2*k1k3*tdip**2*zjm - 10*k1k3*p2k1*udip*zjm - 6*k1k3*sdip*
     &   udip*zjm + 4*p2k3*sdip*udip*zjm - k1k3*tdip*udip*zjm - 2*k1k3*udip**2*zjm + M2*(2*p2k2*sdip*zim - 2*p2k3*sdip
     &   *zjm + k1k3*(4*p2k1 - 13*sdip - 3*tdip + 5*udip)*zjm)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(16*sdip**3))

!************************************************   2 Re(Ms Vgg Mt^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_a = stsT2xx_M2Dipole_SpinTensor_D23_a + 2*dble((GGG*Nc*(-1 + Nc
     &   **2)*(-(ab1*M2**2*p2k2*zim) + 4*ab1*M2*p2k1*p2k2*zim - 4*ab1*p2k1**2*p2k2*zim - M2*p2k2*sdip*zim - 4
     &   *ab1*M2*p2k2*sdip*zim + 2*p2k1*p2k2*sdip*zim + ab1*p2k2*sdip**2*zim + 4*ab1*p2k1*p2k2*tdip*zim - p2k2
     &   *sdip*tdip*zim + 2*ab1*M2*p2k2*udip*zim - 8*ab1*p2k1*p2k2*udip*zim + 2*p2k2*sdip*udip*zim - ab1*p2k2
     &   *udip**2*zim + k1k2*((-1 + ab1)*M2**2 + 2*sdip**2 + 4*p2k1*tdip - 6*p2k1*udip - tdip*udip + M2*((2 - 4
     &   *ab1)*p2k1 + 4*(-2 + ab1)*sdip + tdip + udip - 2*ab1*udip) + ab1*(4*p2k1**2 - sdip**2 - 4*p2k1*tdip 
     &   + 8*p2k1*udip + udip**2))*zim + k1k3*M2**2*zjm - ab1*k1k3*M2**2*zjm - 2*k1k3*M2*p2k1*zjm + 4*ab1*k1k3
     &   *M2*p2k1*zjm - 4*ab1*k1k3*p2k1**2*zjm + ab1*M2**2*p2k3*zjm - 4*ab1*M2*p2k1*p2k3*zjm + 4*ab1*p2k1**2*p2k3
     &   *zjm + 8*k1k3*M2*sdip*zjm - 4*ab1*k1k3*M2*sdip*zjm + M2*p2k3*sdip*zjm + 4*ab1*M2*p2k3*sdip*zjm - 2
     &   *p2k1*p2k3*sdip*zjm - 2*k1k3*sdip**2*zjm + ab1*k1k3*sdip**2*zjm - ab1*p2k3*sdip**2*zjm - k1k3*M2*tdip
     &   *zjm - 4*k1k3*p2k1*tdip*zjm + 4*ab1*k1k3*p2k1*tdip*zjm - 4*ab1*p2k1*p2k3*tdip*zjm + p2k3*sdip*tdip*zjm 
     &   - k1k3*M2*udip*zjm + 2*ab1*k1k3*M2*udip*zjm + 6*k1k3*p2k1*udip*zjm - 8*ab1*k1k3*p2k1*udip*zjm - 2*ab1
     &   *M2*p2k3*udip*zjm + 8*ab1*p2k1*p2k3*udip*zjm - 2*p2k3*sdip*udip*zjm + k1k3*tdip*udip*zjm - ab1*k1k3
     &   *udip**2*zjm + ab1*p2k3*udip**2*zjm - ab2*(M2**2 + 4*p2k1**2 - sdip**2 - 4*p2k1*tdip + 8*p2k1*udip + udip
     &   **2 - 2*M2*(2*p2k1 - 2*sdip + udip))*(p1k2*zim - p1k3*zjm))*((-2 + ab1)*k1k2*zim + ab2*(-(p1k2*zim) + p1k3*zjm) 
     &   - (-2 + ab1)*(p2k2*zim + k1k3*zjm - p2k3*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(8*sdip**2*(M2 - tdip)))

!************************************************   2 Re(Ms Vgg Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_a = stsT2xx_M2Dipole_SpinTensor_D23_a + 2*dble(0)

!************************************************   2 Re(Ms Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_a = stsT2xx_M2Dipole_SpinTensor_D23_a + 2*dble((GGG*Nc*(-1 + Nc
     &   **2)*(2*k1k2**2*sdip*zim**2 - 4*k1k2*p2k2*sdip*zim**2 + 2*k1k2**2*tdip*zim**2 - 2*k1k2**2*udip*zim*
     &   *2 - 4*k1k2*k1k3*sdip*zim*zjm + k2k3*M2*sdip*zim*zjm - 2*k2k3*p2k1*sdip*zim*zjm + 4*k1k3*p2k2*sdip*zim
     &   *zjm + 4*k1k2*p2k3*sdip*zim*zjm - 4*k1k2*k1k3*tdip*zim*zjm + k2k3*sdip*tdip*zim*zjm + 4*k1k2*k1k3*udip
     &   *zim*zjm - 2*k2k3*sdip*udip*zim*zjm + 2*k1k3**2*sdip*zjm**2 - 4*k1k3*p2k3*sdip*zjm**2 + 2*k1k3**2*tdip
     &   *zjm**2 - 2*k1k3**2*udip*zjm**2 - ab2**2*(M2 - 2*p2k1 - sdip - udip)*(p1k2*zim - p1k3*zjm)**2 - ab1
     &   **2*(M2 - 2*p2k1 - sdip - udip)*(k1k2*zim - p2k2*zim + (-k1k3 + p2k3)*zjm)**2 + ab2*(p1k2*zim - p1k3
     &   *zjm)*(-2*p2k2*sdip*zim + 3*k1k2*(sdip + tdip - udip)*zim + 2*p2k3*sdip*zjm - 3*k1k3*(sdip + tdip - udip)
     &   *zjm) + ab1*(k1k2*zim - p2k2*zim + (-k1k3 + p2k3)*zjm)*(2*p2k2*sdip*zim - 3*k1k2*(sdip + tdip - udip)
     &   *zim + 3*k1k3*sdip*zjm - 2*p2k3*sdip*zjm + 3*k1k3*tdip*zjm - 3*k1k3*udip*zjm + 2*ab2*(M2 - 2*p2k1 
     &   - sdip - udip)*(p1k2*zim - p1k3*zjm)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(8*sdip**2))

!************************************************   Mt^ast Vgg Mt   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_a = stsT2xx_M2Dipole_SpinTensor_D23_a +(((-1 + Nc**2)**2*(M2*(
     &   -2*p2k1 + sdip) + 2*p2k1*udip)*(-(ab2*p1k2*zim) + ab2*p1k3*zjm + (-2 + ab1)*(k1k2*zim - p2k2*zim + (
     &   -k1k3 + p2k3)*zjm))**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(2*Nc*sdip*(M2 - tdip)**2))

!************************************************   2 Re(Mt Vgg Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_a = stsT2xx_M2Dipole_SpinTensor_D23_a + 2*dble(((-1 + Nc**2)*(M2
     &   **2 - sdip**2 - 2*p2k1*tdip + 2*M2*(p2k1 + sdip - udip) + udip**2)*(-(ab2*p1k2*zim) + 2*p2k2*zim + ab2
     &   *p1k3*zjm - 2*p2k3*zjm + ab1*(k1k2*zim - p2k2*zim - k1k3*zjm + p2k3*zjm))*(ab2*p1k2*zim - ab2*p1k3
     &   *zjm - (-2 + ab1)*(k1k2*zim - p2k2*zim + (-k1k3 + p2k3)*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **4)/(4*Nc*sdip*(M2 - tdip)*(M2 - udip)))

!************************************************   2 Re(Mt Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_a = stsT2xx_M2Dipole_SpinTensor_D23_a + 2*dble(((2 - 3*Nc**2 
     &   + Nc**4)*(2*ab2*p1k2*p2k1*zim + 2*ab1*p2k1*(-k1k2 + p2k2)*zim - p2k2*sdip*zim + k1k2*(M2 - udip)*zim 
     &   - 2*ab2*p1k3*p2k1*zjm + 2*ab1*p2k1*(k1k3 - p2k3)*zjm + p2k3*sdip*zjm + k1k3*(-M2 + udip)*zjm)*(ab2*p1k2
     &   *zim - ab2*p1k3*zjm - (-2 + ab1)*(k1k2*zim - p2k2*zim + (-k1k3 + p2k3)*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(4*Nc*sdip*(M2 - tdip)))

!************************************************   Mu^ast Vgg Mu   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_a = stsT2xx_M2Dipole_SpinTensor_D23_a +(((-1 + Nc**2)*(M2**2 
     &   + tdip*udip - M2*(sdip + tdip + udip))*(-(ab2*p1k2*zim) + 2*p2k2*zim + ab2*p1k3*zjm - 2*p2k3*zjm + ab1
     &   *(k1k2*zim - p2k2*zim - k1k3*zjm + p2k3*zjm))**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(2*Nc*sdip*(M2 - udip)**2))

!************************************************   2 Re(Mu Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_a = stsT2xx_M2Dipole_SpinTensor_D23_a + 2*dble(((-1 + Nc**2)*(
     &   -(ab1*M2*p2k2*zim) - p2k2*sdip*zim + ab1*p2k2*sdip*zim + ab1*p2k2*udip*zim + k1k2*((-1 + ab1)*M2 + sdip 
     &   - ab1*sdip + tdip - ab1*udip)*zim + k1k3*M2*zjm - ab1*k1k3*M2*zjm + ab1*M2*p2k3*zjm - k1k3*sdip*zjm 
     &   + ab1*k1k3*sdip*zjm + p2k3*sdip*zjm - ab1*p2k3*sdip*zjm - k1k3*tdip*zjm + ab1*k1k3*udip*zjm - ab1*p2k3
     &   *udip*zjm - ab2*(M2 - sdip - udip)*(p1k2*zim - p1k3*zjm))*(-(ab2*p1k2*zim) + 2*p2k2*zim + ab2*p1k3
     &   *zjm - 2*p2k3*zjm + ab1*(k1k2*zim - p2k2*zim - k1k3*zjm + p2k3*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(2*Nc*sdip*(M2 - udip)))

!************************************************   Mv^ast Vgg Mv   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_a = stsT2xx_M2Dipole_SpinTensor_D23_a +(((4 - 5*Nc**2 + Nc**4)
     &   *(-(k2k3*sdip*zim*zjm) - 2*ab2*(k1k2*zim - k1k3*zjm)*(p1k2*zim - p1k3*zjm) + 2*ab1*(k1k2*zim - k1k3*zjm)
     &   *(k1k2*zim - p2k2*zim - k1k3*zjm + p2k3*zjm))*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(4*Nc*sdip))

      end

      double precision function stsT2xx_M2Dipole_SpinTensor_D13_a(sdip,tdip,udip)
      implicit none
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"
#include "common/stsT2xx_GenCouplings.h"

      double precision sdip,tdip,udip,M2,zi,zM,zim,zjm
      double precision xija,Rdip,kappa,ab1,ab2

      M2=MSf2(iSf1,itt1,igen1)

      xija = (-k1k3 + p1k1 + p1k3)/(p1k1 + p1k3)

      Rdip = dsqrt(-4*M2*tdip*xija**2 + (-M2 + tdip + 2*M2*xija)**2)/kappa(tdip,0d0,M2)

      zi =  p1k1/(p1k1 + p1k3)

      zM = (1-Rdip)/2d0

      ab1 = 1 - (M2 + tdip)/(2*tdip) + (tdip - M2 + 2*M2*xija)/(2*Rdip*tdip)
      ab2 = xija/Rdip

      zim = zi - zM
      zjm = 1 - zi - zM

!************************************************   Ms^ast Vgg Ms   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_a = ((GGG**2*Nc*(-1 + Nc**2)
     &   *(2*k1k2*M2**2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*zim**2 + 3*M2**2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))**2*zim
     &   **2 + 2*k1k2*M2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 
     &   + ab1)*ab1*(M2 - 2*p2k2))*zim**2 - 4*M2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))**2*(ab2**2*M2 + (-1 + 2*ab1)
     &   *ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*zim**2 - 8*k1k2**2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(
     &   -p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))**2*zim**2 - 8*k1k2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2
     &   **2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))**2*zim**2 - 4*k1k2**2*M2*p2k2
     &   *zim**2 - 10*k1k2*M2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*p2k2*zim**2 - 8*M2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))
     &   **2*p2k2*zim**2 + 8*k1k2**2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2
     &   *p2k2))*p2k2*zim**2 + 12*k1k2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) 
     &   + (-1 + ab1)*ab1*(M2 - 2*p2k2))*p2k2*zim**2 + 8*(ab2*p1k1 + ab1*(-k1k2 + p2k1))**2*(ab2**2*M2 + (
     &   -1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*p2k2*zim**2 + 4*k1k2*(ab2*p1k1 + ab1
     &   *(-k1k2 + p2k1))*p2k2**2*zim**2 + 4*(ab2*p1k1 + ab1*(-k1k2 + p2k1))**2*p2k2**2*zim**2 + 4*k1k2*M2*p2k1
     &   *sdip*zim**2 - 4*k1k2*M2*(ab2*p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*sdip*zim**2 + 16*k1k2*M2*(ab2*p1k1 + ab1
     &   *(-k1k2 + p2k1))*sdip*zim**2 + 2*M2*p2k1*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*sdip*zim**2 - 2*M2*(ab2*p1k1 
     &   + (-1 + ab1)*(-k1k2 + p2k1))*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*sdip*zim**2 + 4*M2*(ab2*p1k1 + ab1*(
     &   -k1k2 + p2k1))**2*sdip*zim**2 - 4*k1k2*p2k1*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)
     &   *ab1*(M2 - 2*p2k2))*sdip*zim**2 + 4*k1k2*(ab2*p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*(ab2**2*M2 + (-1 + 2
     &   *ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*sdip*zim**2 - 2*p2k1*(ab2*p1k1 + ab1*(-k1k2 
     &   + p2k1))*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*sdip*zim**2 + 2
     &   *(ab2*p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2**2*M2 + (-1 + 2*ab1)*ab2
     &   *(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*sdip*zim**2 - 4*k1k2*p2k1*p2k2*sdip*zim**2 + 4*k1k2*(ab2
     &   *p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*p2k2*sdip*zim**2 - 2*p2k1*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*p2k2*sdip
     &   *zim**2 + 2*(ab2*p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*p2k2*sdip*zim*
     &   *2 + 4*k1k2**2*(4*M2 - sdip)*sdip*zim**2 - 4*k1k2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*sdip**2*zim**2 - (ab2
     &   *p1k1 + ab1*(-k1k2 + p2k1))**2*sdip**2*zim**2 + 4*k1k2**2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) 
     &   + (-1 + ab1)*ab1*(M2 - 2*p2k2))*(M2 - udip)*zim**2 - 4*k1k2*M2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*udip
     &   *zim**2 - 6*M2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))**2*udip*zim**2 - 2*k1k2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))
     &   *(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*udip*zim**2 + 4*(ab2
     &   *p1k1 + ab1*(-k1k2 + p2k1))**2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))
     &   *udip*zim**2 + 4*k1k2**2*p2k2*udip*zim**2 + 10*k1k2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*p2k2*udip*zim
     &   **2 + 8*(ab2*p1k1 + ab1*(-k1k2 + p2k1))**2*p2k2*udip*zim**2 - 4*k1k2*p2k1*sdip*udip*zim**2 + 4*k1k2
     &   *(ab2*p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*sdip*udip*zim**2 - 2*p2k1*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*sdip
     &   *udip*zim**2 + 2*(ab2*p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*sdip*udip*zim
     &   **2 + 2*k1k2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*udip**2*zim**2 + 3*(ab2*p1k1 + ab1*(-k1k2 + p2k1))**2*udip
     &   **2*zim**2 - 2*k2k3*M2**2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*zim*zjm - 2*k2k3*M2*(ab2*p1k1 + ab1*(-k1k2 
     &   + p2k1))*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*zim*zjm + 16
     &   *k1k2*k2k3*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))**2*zim*zjm + 8
     &   *k2k3*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 
     &   - 2*p2k2))**2*zim*zjm + 8*k1k2*k2k3*M2*p2k2*zim*zjm + 10*k2k3*M2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*p2k2
     &   *zim*zjm - 16*k1k2*k2k3*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))
     &   *p2k2*zim*zjm - 12*k2k3*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) 
     &   + (-1 + ab1)*ab1*(M2 - 2*p2k2))*p2k2*zim*zjm - 4*k2k3*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*p2k2**2*zim*zjm 
     &   - 2*k1k2*M2**2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*zim*zjm - 6*M2**2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2
     &   *p1k3 + ab1*(-k2k3 + p2k3))*zim*zjm - 2*k1k2*M2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 
     &   + ab1)*ab1*(M2 - 2*p2k2))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*zim*zjm + 8*M2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))
     &   *(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*(ab2*p1k3 + ab1*(-k2k3 
     &   + p2k3))*zim*zjm + 8*k1k2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))
     &   **2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*zim*zjm + 10*k1k2*M2*p2k2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*zim
     &   *zjm + 16*M2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*p2k2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*zim*zjm - 12*k1k2
     &   *(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*p2k2*(ab2*p1k3 + ab1*(-k2k3 
     &   + p2k3))*zim*zjm - 16*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) 
     &   + (-1 + ab1)*ab1*(M2 - 2*p2k2))*p2k2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*zim*zjm - 4*k1k2*p2k2**2*(ab2
     &   *p1k3 + ab1*(-k2k3 + p2k3))*zim*zjm - 8*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*p2k2**2*(ab2*p1k3 + ab1*(-k2k3 
     &   + p2k3))*zim*zjm - 2*k1k3*M2**2*sdip*zim*zjm - 4*k2k3*M2*p2k1*sdip*zim*zjm + 4*k2k3*M2*(ab2*p1k1 + (
     &   -1 + ab1)*(-k1k2 + p2k1))*sdip*zim*zjm - 16*k2k3*M2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*sdip*zim*zjm + 4
     &   *k1k3*M2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*sdip*zim*zjm + 4
     &   *k2k3*p2k1*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*sdip*zim*zjm 
     &   - 4*k2k3*(ab2*p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)
     &   *ab1*(M2 - 2*p2k2))*sdip*zim*zjm - 2*k1k3*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)
     &   *ab1*(M2 - 2*p2k2))**2*sdip*zim*zjm + 4*k1k3*M2*p2k2*sdip*zim*zjm + 4*k2k3*p2k1*p2k2*sdip*zim*zjm 
     &   - 4*k2k3*(ab2*p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*p2k2*sdip*zim*zjm - 4*k1k3*(ab2**2*M2 + (-1 + 2*ab1)*ab2
     &   *(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*p2k2*sdip*zim*zjm - 2*k1k3*p2k2**2*sdip*zim*zjm - 4
     &   *k1k2*M2*p2k3*sdip*zim*zjm - 2*M2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*p2k3*sdip*zim*zjm + 4*k1k2*(ab2**2
     &   *M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*p2k3*sdip*zim*zjm + 2*(ab2*p1k1 
     &   + ab1*(-k1k2 + p2k1))*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*p2k3
     &   *sdip*zim*zjm + 4*k1k2*p2k2*p2k3*sdip*zim*zjm + 2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*p2k2*p2k3*sdip*zim
     &   *zjm + 4*k1k2*M2*(ab2*p1k3 + (-1 + ab1)*(-k2k3 + p2k3))*sdip*zim*zjm + 2*M2*(ab2*p1k1 + ab1*(-k1k2 
     &   + p2k1))*(ab2*p1k3 + (-1 + ab1)*(-k2k3 + p2k3))*sdip*zim*zjm - 4*k1k2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(
     &   -p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*(ab2*p1k3 + (-1 + ab1)*(-k2k3 + p2k3))*sdip*zim*zjm - 2
     &   *(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2
     &   *p2k2))*(ab2*p1k3 + (-1 + ab1)*(-k2k3 + p2k3))*sdip*zim*zjm - 4*k1k2*p2k2*(ab2*p1k3 + (-1 + ab1)*(-k2k3 
     &   + p2k3))*sdip*zim*zjm - 2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*p2k2*(ab2*p1k3 + (-1 + ab1)*(-k2k3 + p2k3))
     &   *sdip*zim*zjm - 16*k1k2*M2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip*zim*zjm - 2*M2*p2k1*(ab2*p1k3 + ab1
     &   *(-k2k3 + p2k3))*sdip*zim*zjm + 2*M2*(ab2*p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*(ab2*p1k3 + ab1*(-k2k3 
     &   + p2k3))*sdip*zim*zjm - 8*M2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip*zim
     &   *zjm + 2*p2k1*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*(ab2*p1k3 
     &   + ab1*(-k2k3 + p2k3))*sdip*zim*zjm - 2*(ab2*p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*(ab2**2*M2 + (-1 + 2*ab1)
     &   *ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip*zim*zjm + 2
     &   *p2k1*p2k2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip*zim*zjm - 2*(ab2*p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*p2k2
     &   *(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip*zim*zjm - 8*k1k2*k2k3*(4*M2 - sdip)*sdip*zim*zjm + 4*k2k3*(ab2
     &   *p1k1 + ab1*(-k1k2 + p2k1))*sdip**2*zim*zjm + 4*k1k2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip**2*zim*zjm 
     &   + 2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip**2*zim*zjm - 8*k1k2*k2k3*(ab2
     &   **2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*(M2 - udip)*zim*zjm + 4*k2k3
     &   *M2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*udip*zim*zjm + 2*k2k3*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2**2*M2 
     &   + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*udip*zim*zjm - 8*k1k2*k2k3*p2k2*udip
     &   *zim*zjm - 10*k2k3*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*p2k2*udip*zim*zjm + 4*k1k2*M2*(ab2*p1k3 + ab1*(
     &   -k2k3 + p2k3))*udip*zim*zjm + 12*M2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*udip
     &   *zim*zjm + 2*k1k2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*(ab2
     &   *p1k3 + ab1*(-k2k3 + p2k3))*udip*zim*zjm - 8*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2**2*M2 + (-1 + 2*ab1)
     &   *ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*udip*zim*zjm - 10
     &   *k1k2*p2k2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*udip*zim*zjm - 16*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*p2k2*(ab2
     &   *p1k3 + ab1*(-k2k3 + p2k3))*udip*zim*zjm + 4*k1k3*M2*sdip*udip*zim*zjm + 4*k2k3*p2k1*sdip*udip*zim
     &   *zjm - 4*k2k3*(ab2*p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*sdip*udip*zim*zjm - 4*k1k3*(ab2**2*M2 + (-1 + 2
     &   *ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*sdip*udip*zim*zjm - 4*k1k3*p2k2*sdip*udip*zim
     &   *zjm + 4*k1k2*p2k3*sdip*udip*zim*zjm + 2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*p2k3*sdip*udip*zim*zjm - 4*k1k2
     &   *(ab2*p1k3 + (-1 + ab1)*(-k2k3 + p2k3))*sdip*udip*zim*zjm - 2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2
     &   *p1k3 + (-1 + ab1)*(-k2k3 + p2k3))*sdip*udip*zim*zjm + 2*p2k1*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip*udip
     &   *zim*zjm - 2*(ab2*p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip*udip*zim*zjm 
     &   - 2*k2k3*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*udip**2*zim*zjm - 2*k1k2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))
     &   *udip**2*zim*zjm - 6*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*udip**2*zim*zjm 
     &   - 2*k1k3*sdip*udip**2*zim*zjm - 8*k2k3**2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1
     &   *(M2 - 2*p2k2))**2*zjm**2 - 4*k2k3**2*M2*p2k2*zjm**2 + 8*k2k3**2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 
     &   + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*p2k2*zjm**2 + 2*k2k3*M2**2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))
     &   *zjm**2 + 2*k2k3*M2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*(ab2
     &   *p1k3 + ab1*(-k2k3 + p2k3))*zjm**2 - 8*k2k3*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)
     &   *ab1*(M2 - 2*p2k2))**2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*zjm**2 - 10*k2k3*M2*p2k2*(ab2*p1k3 + ab1*(-k2k3 
     &   + p2k3))*zjm**2 + 12*k2k3*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))
     &   *p2k2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*zjm**2 + 4*k2k3*p2k2**2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*zjm*
     &   *2 + 3*M2**2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))**2*zjm**2 - 4*M2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) 
     &   + (-1 + ab1)*ab1*(M2 - 2*p2k2))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))**2*zjm**2 - 8*M2*p2k2*(ab2*p1k3 
     &   + ab1*(-k2k3 + p2k3))**2*zjm**2 + 8*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 
     &   - 2*p2k2))*p2k2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))**2*zjm**2 + 4*p2k2**2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))
     &   **2*zjm**2 + 4*k2k3*M2*p2k3*sdip*zjm**2 - 4*k2k3*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 
     &   + ab1)*ab1*(M2 - 2*p2k2))*p2k3*sdip*zjm**2 - 4*k2k3*p2k2*p2k3*sdip*zjm**2 - 4*k2k3*M2*(ab2*p1k3 + (-1 
     &   + ab1)*(-k2k3 + p2k3))*sdip*zjm**2 + 4*k2k3*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)
     &   *ab1*(M2 - 2*p2k2))*(ab2*p1k3 + (-1 + ab1)*(-k2k3 + p2k3))*sdip*zjm**2 + 4*k2k3*p2k2*(ab2*p1k3 + (-1 
     &   + ab1)*(-k2k3 + p2k3))*sdip*zjm**2 + 16*k2k3*M2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip*zjm**2 + 2*M2*p2k3
     &   *(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip*zjm**2 - 2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (
     &   -1 + ab1)*ab1*(M2 - 2*p2k2))*p2k3*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip*zjm**2 - 2*p2k2*p2k3*(ab2*p1k3 
     &   + ab1*(-k2k3 + p2k3))*sdip*zjm**2 - 2*M2*(ab2*p1k3 + (-1 + ab1)*(-k2k3 + p2k3))*(ab2*p1k3 + ab1*(-k2k3 
     &   + p2k3))*sdip*zjm**2 + 2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))
     &   *(ab2*p1k3 + (-1 + ab1)*(-k2k3 + p2k3))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip*zjm**2 + 2*p2k2*(ab2*p1k3 
     &   + (-1 + ab1)*(-k2k3 + p2k3))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip*zjm**2 + 4*M2*(ab2*p1k3 + ab1*(-k2k3 
     &   + p2k3))**2*sdip*zjm**2 + 4*k2k3**2*(4*M2 - sdip)*sdip*zjm**2 - 4*k2k3*(ab2*p1k3 + ab1*(-k2k3 + p2k3))
     &   *sdip**2*zjm**2 - (ab2*p1k3 + ab1*(-k2k3 + p2k3))**2*sdip**2*zjm**2 + 4*k2k3**2*(ab2**2*M2 + (-1 + 2*ab1)
     &   *ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*(M2 - udip)*zjm**2 + 4*k2k3**2*p2k2*udip*zjm**2 
     &   - 4*k2k3*M2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*udip*zjm**2 - 2*k2k3*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 
     &   + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*udip*zjm**2 + 10*k2k3*p2k2*(ab2
     &   *p1k3 + ab1*(-k2k3 + p2k3))*udip*zjm**2 - 6*M2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))**2*udip*zjm**2 + 4*
     &   (ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*(ab2*p1k3 + ab1*(-k2k3 
     &   + p2k3))**2*udip*zjm**2 + 8*p2k2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))**2*udip*zjm**2 - 4*k2k3*p2k3*sdip*udip*zjm**2 
     &   + 4*k2k3*(ab2*p1k3 + (-1 + ab1)*(-k2k3 + p2k3))*sdip*udip*zjm**2 - 2*p2k3*(ab2*p1k3 + ab1*(-k2k3 + p2k3))
     &   *sdip*udip*zjm**2 + 2*(ab2*p1k3 + (-1 + ab1)*(-k2k3 + p2k3))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip*udip*zjm**2
     &   + 2*k2k3*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*udip**2*zjm**2 + 3*(ab2*p1k3 + ab1*(-k2k3 + p2k3))**2*udip**2*zjm**2)
     &  *SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(8*sdip**3))

!************************************************   2 Re(Ms Vgg Mt^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_a = stsT2xx_M2Dipole_SpinTensor_D13_a + 2*dble((GGG*Nc*(-1 + Nc
     &   **2)*((-2 + ab1)*k1k2*zim + ab2*(-(p1k1*zim) + p1k3*zjm) - (-2 + ab1)*(p2k1*zim + k2k3*zjm - p2k3*zjm))
     &   *(8*ab2**3*M2*p2k2*(-(p1k1*zim) + p1k3*zjm) + 8*ab1**3*(M2 - 2*p2k2)*p2k2*(k1k2*zim - p2k1*zim + (
     &   -k2k3 + p2k3)*zjm) - 2*ab2**2*(4*p1k1*(p1k2 - p1p2)*p2k2*zim - M2*p2k1*sdip*zim + k1k2*M2*(M2 + 4*p2k2 
     &   - udip)*zim - k2k3*M2**2*zjm - 4*k2k3*M2*p2k2*zjm - 4*p1k2*p1k3*p2k2*zjm + 4*p1k3*p1p2*p2k2*zjm + M2
     &   *p2k3*sdip*zjm + k2k3*M2*udip*zjm) + 2*(p2k1*p2k2*sdip*zim + p2k1*sdip*udip*zim + k1k2*(3*M2*p2k2 - 4
     &   *M2*sdip + sdip**2 - 3*p2k2*udip)*zim - p2k2*p2k3*sdip*zjm - k2k3*sdip**2*zjm + 3*k2k3*p2k2*udip*zjm 
     &   - p2k3*sdip*udip*zjm + M2*(-(p2k1*sdip*zim) - 3*k2k3*p2k2*zjm + 4*k2k3*sdip*zjm + p2k3*sdip*zjm)) - 2
     &   *ab1**2*(k1k2*(M2**2 + M2*(6*p2k2 - udip) + 2*p2k2*(4*ab2*p1k2 - 4*ab2*p1p2 - 8*p2k2 + udip))*zim + 4
     &   *ab2*p2k2*(M2*p1k1*zim - 2*p1k2*p2k1*zim + 2*p1p2*p2k1*zim - 2*p1k1*p2k2*zim - 2*k2k3*p1k2*zjm - M2
     &   *p1k3*zjm + 2*k2k3*p1p2*zjm + 2*p1k3*p2k2*zjm + 2*p1k2*p2k3*zjm - 2*p1p2*p2k3*zjm) - (M2 - 2*p2k2)*(p2k1
     &   *(4*p2k2 + sdip)*zim - p2k3*(4*p2k2 + sdip)*zjm + k2k3*(M2 + 8*p2k2 - udip)*zjm)) + ab2*(-4*p1k1*p2k2
     &   **2*zim + 2*p1k2*p2k1*sdip*zim - 2*p1p2*p2k1*sdip*zim + p1k1*sdip**2*zim - 2*k1k2*(p1k2 - p1p2)*(4
     &   *p2k2 - udip)*zim - 8*p1k1*p2k2*udip*zim - p1k1*udip**2*zim + 8*k2k3*p1k2*p2k2*zjm - 8*k2k3*p1p2*p2k2
     &   *zjm + 4*p1k3*p2k2**2*zjm - 2*p1k2*p2k3*sdip*zjm + 2*p1p2*p2k3*sdip*zjm - p1k3*sdip**2*zjm - 2*k2k3*p1k2
     &   *udip*zjm + 2*k2k3*p1p2*udip*zjm + 8*p1k3*p2k2*udip*zjm + p1k3*udip**2*zjm + M2**2*(-(p1k1*zim) + p1k3
     &   *zjm) - 2*M2*(k1k2*(p1k2 - p1p2)*zim - p1k1*(4*p2k2 - 2*sdip + udip)*zim + k2k3*(-p1k2 + p1p2)*zjm 
     &   + p1k3*(4*p2k2 - 2*sdip + udip)*zjm)) + ab1*(16*ab2**2*p1k1*p1k2*p2k2*zim - 16*ab2**2*p1k1*p1p2*p2k2
     &   *zim - 8*ab2*p1k2*p2k1*p2k2*zim + 8*ab2*p1p2*p2k1*p2k2*zim - 16*ab2*p1k1*p2k2**2*zim - 4*p2k1*p2k2**2
     &   *zim - 4*ab2*p1k2*p2k1*sdip*zim + 4*ab2*p1p2*p2k1*sdip*zim + 4*p2k1*p2k2*sdip*zim + p2k1*sdip**2*zim 
     &   - 8*p2k1*p2k2*udip*zim - p2k1*udip**2*zim + k1k2*(3*M2**2 - 12*p2k2**2 - sdip**2 + 4*ab2*(p1k2 - p1p2)
     &   *(6*p2k2 - udip) + 4*M2*(ab2*(p1k2 - p1p2) - p2k2 + 2*ab2**2*p2k2 + sdip - udip) + 12*p2k2*udip + udip
     &   **2)*zim - 24*ab2*k2k3*p1k2*p2k2*zjm - 16*ab2**2*p1k2*p1k3*p2k2*zjm + 24*ab2*k2k3*p1p2*p2k2*zjm + 16
     &   *ab2**2*p1k3*p1p2*p2k2*zjm + 12*k2k3*p2k2**2*zjm + 16*ab2*p1k3*p2k2**2*zjm + 8*ab2*p1k2*p2k2*p2k3*zjm 
     &   - 8*ab2*p1p2*p2k2*p2k3*zjm + 4*p2k2**2*p2k3*zjm + 4*ab2*p1k2*p2k3*sdip*zjm - 4*ab2*p1p2*p2k3*sdip*zjm 
     &   - 4*p2k2*p2k3*sdip*zjm + k2k3*sdip**2*zjm - p2k3*sdip**2*zjm + 4*ab2*k2k3*p1k2*udip*zjm - 4*ab2*k2k3
     &   *p1p2*udip*zjm - 12*k2k3*p2k2*udip*zjm + 8*p2k2*p2k3*udip*zjm - k2k3*udip**2*zjm + p2k3*udip**2*zjm 
     &   + M2**2*(-(p2k1*zim) + (-3*k2k3 + p2k3)*zjm) + M2*(2*p2k1*(4*p2k2 - 3*sdip + udip)*zim - 2*p2k3*(4*p2k2 
     &   - 3*sdip + udip)*zjm + 4*k2k3*(p2k2 - sdip + udip)*zjm + 4*ab2*(2*p1k1*p2k2*zim + k2k3*(-p1k2 + p1p2)*zjm 
     &   - 2*p1k3*p2k2*zjm) - 8*ab2**2*p2k2*(p2k1*zim + (k2k3 - p2k3)*zjm))))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **3)/(8*sdip**2*(M2 - tdip)))

!************************************************   2 Re(Ms Vgg Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_a = stsT2xx_M2Dipole_SpinTensor_D13_a + 2*dble(0)

!************************************************   2 Re(Ms Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_a = stsT2xx_M2Dipole_SpinTensor_D13_a + 2*dble((GGG*Nc*(-1 + Nc
     &   **2)*(2*k1k2**2*M2*zim**2 + 3*k1k2*M2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*zim**2 - M2*(ab2*p1k1 + ab1*(
     &   -k1k2 + p2k1))**2*zim**2 - 4*k1k2**2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 
     &   - 2*p2k2))*zim**2 - 6*k1k2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) 
     &   + (-1 + ab1)*ab1*(M2 - 2*p2k2))*zim**2 + 2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))**2*p2k2*zim**2 - 2*k1k2
     &   *p2k1*sdip*zim**2 + 2*k1k2*(ab2*p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*sdip*zim**2 - p2k1*(ab2*p1k1 + ab1
     &   *(-k1k2 + p2k1))*sdip*zim**2 + (ab2*p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*(ab2*p1k1 + ab1*(-k1k2 + p2k1))
     &   *sdip*zim**2 - 2*k1k2**2*udip*zim**2 - 3*k1k2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*udip*zim**2 + (ab2*p1k1 
     &   + ab1*(-k1k2 + p2k1))**2*udip*zim**2 - 4*k1k2*k2k3*M2*zim*zjm - 3*k2k3*M2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))
     &   *zim*zjm + 8*k1k2*k2k3*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))
     &   *zim*zjm + 6*k2k3*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (
     &   -1 + ab1)*ab1*(M2 - 2*p2k2))*zim*zjm - 3*k1k2*M2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*zim*zjm + 2*M2*(ab2
     &   *p1k1 + ab1*(-k1k2 + p2k1))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*zim*zjm + 6*k1k2*(ab2**2*M2 + (-1 + 2*ab1)
     &   *ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*zim*zjm - 4*(ab2
     &   *p1k1 + ab1*(-k1k2 + p2k1))*p2k2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*zim*zjm + 2*k1k3*M2*sdip*zim*zjm + 2
     &   *k2k3*p2k1*sdip*zim*zjm - 2*k2k3*(ab2*p1k1 + (-1 + ab1)*(-k1k2 + p2k1))*sdip*zim*zjm - 2*k1k3*(ab2**2
     &   *M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*sdip*zim*zjm - 2*k1k3*p2k2*sdip
     &   *zim*zjm + 2*k1k2*p2k3*sdip*zim*zjm + (ab2*p1k1 + ab1*(-k1k2 + p2k1))*p2k3*sdip*zim*zjm - 2*k1k2*(ab2
     &   *p1k3 + (-1 + ab1)*(-k2k3 + p2k3))*sdip*zim*zjm - (ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2*p1k3 + (-1 + ab1)
     &   *(-k2k3 + p2k3))*sdip*zim*zjm + p2k1*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip*zim*zjm - (ab2*p1k1 + (-1 
     &   + ab1)*(-k1k2 + p2k1))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*sdip*zim*zjm + 4*k1k2*k2k3*udip*zim*zjm + 3*k2k3
     &   *(ab2*p1k1 + ab1*(-k1k2 + p2k1))*udip*zim*zjm + 3*k1k2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*udip*zim*zjm 
     &   - 2*(ab2*p1k1 + ab1*(-k1k2 + p2k1))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*udip*zim*zjm - 2*k1k3*sdip*udip
     &   *zim*zjm + 2*k2k3**2*M2*zjm**2 - 4*k2k3**2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)
     &   *ab1*(M2 - 2*p2k2))*zjm**2 + 3*k2k3*M2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*zjm**2 - 6*k2k3*(ab2**2*M2 + (
     &   -1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*zjm**2 
     &   - M2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))**2*zjm**2 + 2*p2k2*(ab2*p1k3 + ab1*(-k2k3 + p2k3))**2*zjm**2 - 2
     &   *k2k3*p2k3*sdip*zjm**2 + 2*k2k3*(ab2*p1k3 + (-1 + ab1)*(-k2k3 + p2k3))*sdip*zjm**2 - p2k3*(ab2*p1k3 
     &   + ab1*(-k2k3 + p2k3))*sdip*zjm**2 + (ab2*p1k3 + (-1 + ab1)*(-k2k3 + p2k3))*(ab2*p1k3 + ab1*(-k2k3 + p2k3))
     &   *sdip*zjm**2 - 2*k2k3**2*udip*zjm**2 - 3*k2k3*(ab2*p1k3 + ab1*(-k2k3 + p2k3))*udip*zjm**2 + (ab2*p1k3 
     &   + ab1*(-k2k3 + p2k3))**2*udip*zjm**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(8*sdip**2))

!************************************************   Mt^ast Vgg Mt   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_a = stsT2xx_M2Dipole_SpinTensor_D13_a +(((-1 + Nc**2)**2*(M2*(
     &   -2*p2k2 + sdip) + 2*p2k2*udip)*(-(ab2*p1k1*zim) + ab2*p1k3*zjm + (-2 + ab1)*(k1k2*zim - p2k1*zim + (
     &   -k2k3 + p2k3)*zjm))**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(2*Nc*sdip*(M2 - tdip)**2))

!************************************************   2 Re(Mt Vgg Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_a = stsT2xx_M2Dipole_SpinTensor_D13_a + 2*dble(((-1 + Nc**2)*(M2
     &   **2 + 4*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*p2k2 - sdip**2 
     &   + 2*M2*(sdip - udip) + udip**2)*(-(ab2*p1k1*zim) + 2*p2k1*zim + ab2*p1k3*zjm - 2*p2k3*zjm + ab1*(k1k2
     &   *zim - p2k1*zim - k2k3*zjm + p2k3*zjm))*(ab2*p1k1*zim - ab2*p1k3*zjm - (-2 + ab1)*(k1k2*zim - p2k1*zim 
     &   + (-k2k3 + p2k3)*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(4*Nc*sdip*(M2 - tdip)*(M2 - udip)))

!************************************************   2 Re(Mt Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_a = stsT2xx_M2Dipole_SpinTensor_D13_a + 2*dble(((2 - 3*Nc**2 
     &   + Nc**4)*(2*ab2*p1k1*p2k2*zim + 2*ab1*(-k1k2 + p2k1)*p2k2*zim - p2k1*sdip*zim + k1k2*(M2 - udip)*zim 
     &   - 2*ab2*p1k3*p2k2*zjm + 2*ab1*p2k2*(k2k3 - p2k3)*zjm + p2k3*sdip*zjm + k2k3*(-M2 + udip)*zjm)*(ab2*p1k1
     &   *zim - ab2*p1k3*zjm - (-2 + ab1)*(k1k2*zim - p2k1*zim + (-k2k3 + p2k3)*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(4*Nc*sdip*(M2 - tdip)))

!************************************************   Mu^ast Vgg Mu   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_a = stsT2xx_M2Dipole_SpinTensor_D13_a +(((-1 + Nc**2)*(-(M2*sdip) 
     &   + 2*ab2**2*M2*(M2 - udip) + 2*ab2*(p1k2 - p1p2)*(M2 - udip) + 2*ab1**2*(M2 - 2*p2k2)*(M2 - udip) 
     &   - 2*ab1*(M2 + 2*ab2*p1k2 - 2*ab2*p1p2 - 2*p2k2)*(M2 - udip))*(-(ab2*p1k1*zim) + 2*p2k1*zim + ab2*p1k3
     &   *zjm - 2*p2k3*zjm + ab1*(k1k2*zim - p2k1*zim - k2k3*zjm + p2k3*zjm))**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **4)/(2*Nc*sdip*(M2 - udip)**2))

!************************************************   2 Re(Mu Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_a = stsT2xx_M2Dipole_SpinTensor_D13_a + 2*dble(((-1 + Nc**2)*(
     &   -(ab2*M2*p1k1*zim) + ab1*M2*(k1k2 - p2k1)*zim - 2*k1k2*(ab2**2*M2 + (-1 + 2*ab1)*ab2*(-p1k2 + p1p2) 
     &   + (-1 + ab1)*ab1*(M2 - 2*p2k2))*zim + ab2*p1k1*sdip*zim + (k1k2 - p2k1)*sdip*zim + ab1*(-k1k2 + p2k1)
     &   *sdip*zim + ab2*p1k1*udip*zim + ab1*(-k1k2 + p2k1)*udip*zim + ab2*M2*p1k3*zjm + 2*k2k3*(ab2**2*M2 + (
     &   -1 + 2*ab1)*ab2*(-p1k2 + p1p2) + (-1 + ab1)*ab1*(M2 - 2*p2k2))*zjm + ab1*M2*(-k2k3 + p2k3)*zjm - ab2*p1k3
     &   *sdip*zjm + ab1*(k2k3 - p2k3)*sdip*zjm + (-k2k3 + p2k3)*sdip*zjm - ab2*p1k3*udip*zjm + ab1*(k2k3 
     &   - p2k3)*udip*zjm)*(-(ab2*p1k1*zim) + 2*p2k1*zim + ab2*p1k3*zjm - 2*p2k3*zjm + ab1*(k1k2*zim - p2k1*zim 
     &   - k2k3*zjm + p2k3*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(2*Nc*sdip*(M2 - udip)))

!************************************************   Mv^ast Vgg Mv   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_a = stsT2xx_M2Dipole_SpinTensor_D13_a +(((4 - 5*Nc**2 + Nc**4)
     &   *(-(k1k3*sdip*zim*zjm) - 2*ab2*(k1k2*zim - k2k3*zjm)*(p1k1*zim - p1k3*zjm) + 2*ab1*(k1k2*zim - k2k3*zjm)
     &   *(k1k2*zim - p2k1*zim - k2k3*zjm + p2k3*zjm))*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(4*Nc*sdip))


      end

      double precision function stsT2xx_M2Dipole_SpinTensor_D12_b(sdip,tdip,udip)
      implicit none
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"
#include "common/stsT2xx_GenCouplings.h"

      double precision sdip,tdip,udip,M2,zi,zM,zim,zjm
      double precision xija,Rdip,kappa,ab1,ab2

      M2=MSf2(iSf1,itt1,igen1)

      xija = (-k1k2 + p2k1 + p2k2)/(p2k1 + p2k2)

      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)

      zi =  p2k1/(p2k1 + p2k2)

      zM = (1-Rdip)/2d0

      ab1 = 1 - (M2 + udip)/(2*udip) + (udip - M2 + 2*M2*xija)/(2*Rdip*udip)
      ab2 = xija/Rdip

      zim = zi - zM
      zjm = 1 - zi - zM

      
!************************************************   Ms^ast Vgg Ms   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_b = ((GGG**2*Nc*(-1 + Nc**2)
     &   *(2*k1k3*M2**2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*zim**2 + 2*M2**2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)**2*zim
     &   **2 + 8*k1k3*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*((-1 + ab1)*p1k3 + ab2*p2k3)**2*zim**2 + 8*(ab1*(-k1k3 
     &   + p1k1) + ab2*p2k1)**2*((-1 + ab1)*p1k3 + ab2*p2k3)**2*zim**2 + 16*k1k3**2*M2*sdip*zim**2 + 32*k1k3
     &   *M2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*sdip*zim**2 - 2*M2*p1k1*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*sdip*zim*
     &   *2 + 2*M2*((-1 + ab1)*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*sdip*zim**2 + 8*k1k3
     &   *p1k1*((-1 + ab1)*p1k3 + ab2*p2k3)*sdip*zim**2 - 8*k1k3*((-1 + ab1)*(-k1k3 + p1k1) + ab2*p2k1)*((-1 + ab1)
     &   *p1k3 + ab2*p2k3)*sdip*zim**2 + 4*p1k1*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*((-1 + ab1)*p1k3 + ab2*p2k3)
     &   *sdip*zim**2 - 4*((-1 + ab1)*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*((-1 + ab1)
     &   *p1k3 + ab2*p2k3)*sdip*zim**2 - 8*k1k3*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*sdip**2*zim**2 - 2*(ab1*(-k1k3 
     &   + p1k1) + ab2*p2k1)**2*sdip**2*zim**2 - 8*k1k3**2*sdip*(-2*M2 + sdip)*zim**2 + 8*k1k3**2*((-1 + ab1)
     &   *p1k3 + ab2*p2k3)*(M2 - tdip)*zim**2 + 4*k1k3*p1k1*sdip*(M2 - tdip)*zim**2 - 4*k1k3*((-1 + ab1)*(-k1k3 
     &   + p1k1) + ab2*p2k1)*sdip*(M2 - tdip)*zim**2 - 4*k1k3**2*(M2 - tdip)**2*zim**2 + 6*k1k3*M2*(ab1*(-k1k3 
     &   + p1k1) + ab2*p2k1)*tdip*zim**2 - 2*p1k1*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*sdip*tdip*zim**2 + 2*((-1 
     &   + ab1)*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*sdip*tdip*zim**2 - 4*k1k3*(ab1*(-k1k3 
     &   + p1k1) + ab2*p2k1)*tdip**2*zim**2 - 4*k1k3*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*((-1 + ab1)*p1k3 + ab2
     &   *p2k3)*(2*M2 + 3*tdip - 5*udip)*zim**2 - 8*(ab1*(-k1k3 + p1k1) + ab2*p2k1)**2*((-1 + ab1)*p1k3 + ab2
     &   *p2k3)*(M2 + tdip - 2*udip)*zim**2 + 4*M2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)**2*(2*sdip + tdip - 2*udip)
     &   *zim**2 - 8*k1k3**2*((-1 + ab1)*p1k3 + ab2*p2k3)*(M2 - udip)*zim**2 - 8*k1k3*p1k1*sdip*(M2 - udip)*zim
     &   **2 + 8*k1k3*((-1 + ab1)*(-k1k3 + p1k1) + ab2*p2k1)*sdip*(M2 - udip)*zim**2 + 4*k1k3**2*(M2 - tdip)
     &   *(M2 - udip)*zim**2 - 10*k1k3*M2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*udip*zim**2 + 4*p1k1*(ab1*(-k1k3 + p1k1) 
     &   + ab2*p2k1)*sdip*udip*zim**2 - 4*((-1 + ab1)*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k1k3 + p1k1) + ab2
     &   *p2k1)*sdip*udip*zim**2 + 2*k1k3*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*tdip*udip*zim**2 - 4*(ab1*(-k1k3 + p1k1) 
     &   + ab2*p2k1)**2*tdip*udip*zim**2 + 4*k1k3*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*udip**2*zim**2 + 6*(ab1
     &   *(-k1k3 + p1k1) + ab2*p2k1)**2*udip**2*zim**2 - 2*k2k3*M2**2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*zim*zjm 
     &   - 2*k1k3*M2**2*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*zim*zjm - 4*M2**2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*(ab1
     &   *(-k2k3 + p1k2) + ab2*p2k2)*zim*zjm + 16*M2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k2) + ab2
     &   *p2k2)*((-1 + ab1)*p1k3 + ab2*p2k3)*zim*zjm - 8*k2k3*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*((-1 + ab1)*p1k3 
     &   + ab2*p2k3)**2*zim*zjm - 8*k1k3*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*((-1 + ab1)*p1k3 + ab2*p2k3)**2*zim
     &   *zjm - 16*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*((-1 + ab1)*p1k3 + ab2*p2k3)
     &   **2*zim*zjm - 32*k1k3*k2k3*M2*sdip*zim*zjm - 32*k2k3*M2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*sdip*zim*zjm 
     &   + 2*M2*p1k2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*sdip*zim*zjm - 2*M2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*((-1 
     &   + ab1)*(-k2k3 + p1k2) + ab2*p2k2)*sdip*zim*zjm + 2*M2*p1k1*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*sdip*zim*zjm 
     &   - 2*M2*((-1 + ab1)*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*sdip*zim*zjm - 16*M2
     &   *(ab1*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*sdip*zim*zjm - 8*k2k3*p1k1*((-1 + ab1)
     &   *p1k3 + ab2*p2k3)*sdip*zim*zjm - 8*k1k3*p1k2*((-1 + ab1)*p1k3 + ab2*p2k3)*sdip*zim*zjm + 8*k2k3*((
     &   -1 + ab1)*(-k1k3 + p1k1) + ab2*p2k1)*((-1 + ab1)*p1k3 + ab2*p2k3)*sdip*zim*zjm - 4*p1k2*(ab1*(-k1k3 + p1k1) 
     &   + ab2*p2k1)*((-1 + ab1)*p1k3 + ab2*p2k3)*sdip*zim*zjm + 8*k1k3*((-1 + ab1)*(-k2k3 + p1k2) + ab2
     &   *p2k2)*((-1 + ab1)*p1k3 + ab2*p2k3)*sdip*zim*zjm + 4*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*((-1 + ab1)*(-k2k3 
     &   + p1k2) + ab2*p2k2)*((-1 + ab1)*p1k3 + ab2*p2k3)*sdip*zim*zjm - 4*p1k1*(ab1*(-k2k3 + p1k2) + ab2*p2k2)
     &   *((-1 + ab1)*p1k3 + ab2*p2k3)*sdip*zim*zjm + 4*((-1 + ab1)*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 
     &   + p1k2) + ab2*p2k2)*((-1 + ab1)*p1k3 + ab2*p2k3)*sdip*zim*zjm - 4*k1k2*((-1 + ab1)*p1k3 + ab2*p2k3)
     &   **2*sdip*zim*zjm + 8*k2k3*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*sdip**2*zim*zjm + 8*k1k3*(ab1*(-k2k3 + p1k2) 
     &   + ab2*p2k2)*sdip**2*zim*zjm + 4*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*sdip
     &   **2*zim*zjm + 16*k1k3*k2k3*sdip*(-2*M2 + sdip)*zim*zjm - 16*k1k3*k2k3*((-1 + ab1)*p1k3 + ab2*p2k3)*(M2 
     &   - tdip)*zim*zjm - 4*k2k3*p1k1*sdip*(M2 - tdip)*zim*zjm - 4*k1k3*p1k2*sdip*(M2 - tdip)*zim*zjm + 4*k2k3
     &   *((-1 + ab1)*(-k1k3 + p1k1) + ab2*p2k1)*sdip*(M2 - tdip)*zim*zjm + 4*k1k3*((-1 + ab1)*(-k2k3 + p1k2) 
     &   + ab2*p2k2)*sdip*(M2 - tdip)*zim*zjm - 4*k1k2*((-1 + ab1)*p1k3 + ab2*p2k3)*sdip*(M2 - tdip)*zim*zjm 
     &   + 8*k1k3*k2k3*(M2 - tdip)**2*zim*zjm - k1k2*sdip*(M2 - tdip)**2*zim*zjm - 6*k2k3*M2*(ab1*(-k1k3 + p1k1) 
     &   + ab2*p2k1)*tdip*zim*zjm - 8*M2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k2) + ab2*p2k2)
     &   *tdip*zim*zjm + 16*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*((-1 + ab1)*p1k3 + ab2
     &   *p2k3)*tdip*zim*zjm + 2*p1k2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*sdip*tdip*zim*zjm - 2*(ab1*(-k1k3 + p1k1) 
     &   + ab2*p2k1)*((-1 + ab1)*(-k2k3 + p1k2) + ab2*p2k2)*sdip*tdip*zim*zjm + 2*p1k1*(ab1*(-k2k3 + p1k2) 
     &   + ab2*p2k2)*sdip*tdip*zim*zjm - 2*((-1 + ab1)*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k2) + ab2
     &   *p2k2)*sdip*tdip*zim*zjm + 4*k2k3*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*tdip**2*zim*zjm + 4*k1k3*(ab1*(-k2k3 
     &   + p1k2) + ab2*p2k2)*tdip**2*zim*zjm + 4*k2k3*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*((-1 + ab1)*p1k3 + ab2
     &   *p2k3)*(2*M2 + 3*tdip - 5*udip)*zim*zjm + 4*k1k3*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*((-1 + ab1)*p1k3 + ab2
     &   *p2k3)*(2*M2 + 3*tdip - 5*udip)*zim*zjm - 2*k1k3*M2*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*(16*sdip + 3*tdip 
     &   - 5*udip)*zim*zjm + 16*k1k3*k2k3*((-1 + ab1)*p1k3 + ab2*p2k3)*(M2 - udip)*zim*zjm + 8*k2k3*p1k1*sdip
     &   *(M2 - udip)*zim*zjm + 8*k1k3*p1k2*sdip*(M2 - udip)*zim*zjm - 8*k2k3*((-1 + ab1)*(-k1k3 + p1k1) + ab2
     &   *p2k1)*sdip*(M2 - udip)*zim*zjm - 8*k1k3*((-1 + ab1)*(-k2k3 + p1k2) + ab2*p2k2)*sdip*(M2 - udip)*zim
     &   *zjm + 8*k1k2*((-1 + ab1)*p1k3 + ab2*p2k3)*sdip*(M2 - udip)*zim*zjm - 8*k1k3*k2k3*(M2 - tdip)*(M2 - udip)
     &   *zim*zjm + 4*k1k2*sdip*(M2 - tdip)*(M2 - udip)*zim*zjm - 4*k1k2*sdip*(M2 - udip)**2*zim*zjm + 10
     &   *k2k3*M2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*udip*zim*zjm + 16*M2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 
     &   + p1k2) + ab2*p2k2)*udip*zim*zjm - 32*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k2) + ab2*p2k2)
     &   *((-1 + ab1)*p1k3 + ab2*p2k3)*udip*zim*zjm - 4*p1k2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*sdip*udip*zim
     &   *zjm + 4*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*((-1 + ab1)*(-k2k3 + p1k2) + ab2*p2k2)*sdip*udip*zim*zjm - 4
     &   *p1k1*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*sdip*udip*zim*zjm + 4*((-1 + ab1)*(-k1k3 + p1k1) + ab2*p2k1)*(ab1
     &   *(-k2k3 + p1k2) + ab2*p2k2)*sdip*udip*zim*zjm - 2*k2k3*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*tdip*udip*zim
     &   *zjm - 2*k1k3*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*tdip*udip*zim*zjm + 8*(ab1*(-k1k3 + p1k1) + ab2*p2k1)
     &   *(ab1*(-k2k3 + p1k2) + ab2*p2k2)*tdip*udip*zim*zjm - 4*k2k3*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*udip**2*zim
     &   *zjm - 4*k1k3*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*udip**2*zim*zjm - 12*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*(ab1
     &   *(-k2k3 + p1k2) + ab2*p2k2)*udip**2*zim*zjm + 2*k2k3*M2**2*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*zjm**2 
     &   + 2*M2**2*(ab1*(-k2k3 + p1k2) + ab2*p2k2)**2*zjm**2 + 8*k2k3*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*((-1 + ab1)
     &   *p1k3 + ab2*p2k3)**2*zjm**2 + 8*(ab1*(-k2k3 + p1k2) + ab2*p2k2)**2*((-1 + ab1)*p1k3 + ab2*p2k3)**2
     &   *zjm**2 + 16*k2k3**2*M2*sdip*zjm**2 - 2*M2*p1k2*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*sdip*zjm**2 + 2*M2*((
     &   -1 + ab1)*(-k2k3 + p1k2) + ab2*p2k2)*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*sdip*zjm**2 + 8*k2k3*p1k2*((-1 
     &   + ab1)*p1k3 + ab2*p2k3)*sdip*zjm**2 - 8*k2k3*((-1 + ab1)*(-k2k3 + p1k2) + ab2*p2k2)*((-1 + ab1)*p1k3 
     &   + ab2*p2k3)*sdip*zjm**2 + 4*p1k2*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*((-1 + ab1)*p1k3 + ab2*p2k3)*sdip*zjm
     &   **2 - 4*((-1 + ab1)*(-k2k3 + p1k2) + ab2*p2k2)*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*((-1 + ab1)*p1k3 + ab2
     &   *p2k3)*sdip*zjm**2 - 8*k2k3*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*sdip**2*zjm**2 - 2*(ab1*(-k2k3 + p1k2) + ab2
     &   *p2k2)**2*sdip**2*zjm**2 - 8*k2k3**2*sdip*(-2*M2 + sdip)*zjm**2 + 8*k2k3**2*((-1 + ab1)*p1k3 + ab2
     &   *p2k3)*(M2 - tdip)*zjm**2 + 4*k2k3*p1k2*sdip*(M2 - tdip)*zjm**2 - 4*k2k3*((-1 + ab1)*(-k2k3 + p1k2) 
     &   + ab2*p2k2)*sdip*(M2 - tdip)*zjm**2 - 4*k2k3**2*(M2 - tdip)**2*zjm**2 - 2*p1k2*(ab1*(-k2k3 + p1k2) + ab2
     &   *p2k2)*sdip*tdip*zjm**2 + 2*((-1 + ab1)*(-k2k3 + p1k2) + ab2*p2k2)*(ab1*(-k2k3 + p1k2) + ab2*p2k2)
     &   *sdip*tdip*zjm**2 - 4*k2k3*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*tdip**2*zjm**2 - 4*k2k3*(ab1*(-k2k3 + p1k2) 
     &   + ab2*p2k2)*((-1 + ab1)*p1k3 + ab2*p2k3)*(2*M2 + 3*tdip - 5*udip)*zjm**2 + 2*k2k3*M2*(ab1*(-k2k3 + p1k2) 
     &   + ab2*p2k2)*(16*sdip + 3*tdip - 5*udip)*zjm**2 - 8*(ab1*(-k2k3 + p1k2) + ab2*p2k2)**2*((-1 + ab1)*p1k3 
     &   + ab2*p2k3)*(M2 + tdip - 2*udip)*zjm**2 + 4*M2*(ab1*(-k2k3 + p1k2) + ab2*p2k2)**2*(2*sdip + tdip 
     &   - 2*udip)*zjm**2 - 8*k2k3**2*((-1 + ab1)*p1k3 + ab2*p2k3)*(M2 - udip)*zjm**2 - 8*k2k3*p1k2*sdip*(M2 
     &   - udip)*zjm**2 + 8*k2k3*((-1 + ab1)*(-k2k3 + p1k2) + ab2*p2k2)*sdip*(M2 - udip)*zjm**2 + 4*k2k3**2*(M2 
     &   - tdip)*(M2 - udip)*zjm**2 + 4*p1k2*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*sdip*udip*zjm**2 - 4*((-1 + ab1)*(-k2k3 
     &   + p1k2) + ab2*p2k2)*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*sdip*udip*zjm**2 + 2*k2k3*(ab1*(-k2k3 + p1k2) + ab2*p2k2)
     &   *tdip*udip*zjm**2 - 4*(ab1*(-k2k3 + p1k2) + ab2*p2k2)**2*tdip*udip*zjm**2 + 4*k2k3*(ab1*(-k2k3 + p1k2) 
     &   + ab2*p2k2)*udip**2*zjm**2 + 6*(ab1*(-k2k3 + p1k2) + ab2*p2k2)**2*udip**2*zjm**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(16*sdip**3))

!************************************************   2 Re(Ms Vgg Mt^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_b = stsT2xx_M2Dipole_SpinTensor_D12_b + 2*dble(0)

!************************************************   2 Re(Ms Vgg Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_b = stsT2xx_M2Dipole_SpinTensor_D12_b + 2*dble(-0.125d0*(GGG*Nc
     &   *(-1 + Nc**2)*(4*ab1*M2*(k1k3 - p1k1)*sdip*zim + M2*p1k1*sdip*zim + 2*p1k1*p1k3*sdip*zim - 2*ab1*p1k1
     &   *p1k3*sdip*zim - 4*ab2*M2*p2k1*sdip*zim - 2*ab2*p1k1*p2k3*sdip*zim + ab1*(-k1k3 + p1k1)*sdip**2*zim + ab2
     &   *p2k1*sdip**2*zim + 3*ab1*M2*(k1k3 - p1k1)*tdip*zim + 2*ab1*(k1k3 - p1k1)*p1k3*tdip*zim + 2*ab1**2
     &   *(-k1k3 + p1k1)*p1k3*tdip*zim - 3*ab2*M2*p2k1*tdip*zim - 2*ab2*p1k3*p2k1*tdip*zim + 2*ab1*ab2*p1k3*p2k1
     &   *tdip*zim + 2*ab1*ab2*(-k1k3 + p1k1)*p2k3*tdip*zim + 2*ab2**2*p2k1*p2k3*tdip*zim + p1k1*sdip*tdip*zim 
     &   + 3*ab1*M2*(-k1k3 + p1k1)*udip*zim + 2*ab1**2*(k1k3 - p1k1)*p1k3*udip*zim + 2*ab1*(-k1k3 + p1k1)*p1k3
     &   *udip*zim + 3*ab2*M2*p2k1*udip*zim + 2*ab2*p1k3*p2k1*udip*zim - 2*ab1*ab2*p1k3*p2k1*udip*zim + 2*ab1
     &   *ab2*(k1k3 - p1k1)*p2k3*udip*zim - 2*ab2**2*p2k1*p2k3*udip*zim - 2*p1k1*sdip*udip*zim + 3*ab1*(-k1k3 
     &   + p1k1)*tdip*udip*zim + 3*ab2*p2k1*tdip*udip*zim + 3*ab1*(k1k3 - p1k1)*udip**2*zim - 3*ab2*p2k1*udip
     &   **2*zim + k1k3*(M2**2 - 8*M2*sdip + 2*sdip**2 - 2*(-1 + ab1)*p1k3*(M2 - tdip) - 4*M2*tdip + tdip**2 
     &   + 2*ab2*p2k3*(-M2 + tdip) + 2*M2*udip + 2*tdip*udip - 2*udip**2)*zim - M2*p1k2*sdip*zjm + 4*ab1*M2*(
     &   -k2k3 + p1k2)*sdip*zjm - 2*p1k2*p1k3*sdip*zjm + 2*ab1*p1k2*p1k3*sdip*zjm + 4*ab2*M2*p2k2*sdip*zjm + 2
     &   *ab2*p1k2*p2k3*sdip*zjm + ab1*(k2k3 - p1k2)*sdip**2*zjm - ab2*p2k2*sdip**2*zjm + 3*ab1*M2*(-k2k3 + p1k2)
     &   *tdip*zjm + 2*ab1**2*(k2k3 - p1k2)*p1k3*tdip*zjm + 2*ab1*(-k2k3 + p1k2)*p1k3*tdip*zjm + 3*ab2*M2*p2k2
     &   *tdip*zjm + 2*ab2*p1k3*p2k2*tdip*zjm - 2*ab1*ab2*p1k3*p2k2*tdip*zjm + 2*ab1*ab2*(k2k3 - p1k2)*p2k3
     &   *tdip*zjm - 2*ab2**2*p2k2*p2k3*tdip*zjm - p1k2*sdip*tdip*zjm + 3*ab1*M2*(k2k3 - p1k2)*udip*zjm + 2*ab1
     &   *(k2k3 - p1k2)*p1k3*udip*zjm + 2*ab1**2*(-k2k3 + p1k2)*p1k3*udip*zjm - 3*ab2*M2*p2k2*udip*zjm - 2*ab2
     &   *p1k3*p2k2*udip*zjm + 2*ab1*ab2*p1k3*p2k2*udip*zjm + 2*ab1*ab2*(-k2k3 + p1k2)*p2k3*udip*zjm + 2*ab2*
     &   *2*p2k2*p2k3*udip*zjm + 2*p1k2*sdip*udip*zjm + 3*ab1*(k2k3 - p1k2)*tdip*udip*zjm - 3*ab2*p2k2*tdip*udip
     &   *zjm + 3*ab1*(-k2k3 + p1k2)*udip**2*zjm + 3*ab2*p2k2*udip**2*zjm - k2k3*(M2**2 - 8*M2*sdip + 2*sdip
     &   **2 - 2*(-1 + ab1)*p1k3*(M2 - tdip) - 4*M2*tdip + tdip**2 + 2*ab2*p2k3*(-M2 + tdip) + 2*M2*udip 
     &   + 2*tdip*udip - 2*udip**2)*zjm)*(ab2*p2k1*zim - ab2*p2k2*zjm - (-2 + ab1)*(k1k3*zim - p1k1*zim 
     &   + (-k2k3 + p1k2)*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(sdip**2*(M2 - udip)))

!************************************************   2 Re(Ms Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_b = stsT2xx_M2Dipole_SpinTensor_D12_b + 2*dble((GGG*Nc*(-1 + Nc
     &   **2)*(M2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)**2*zim**2 - 2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)**2*((-1 + ab1)
     &   *p1k3 + ab2*p2k3)*zim**2 - 2*k1k3*p1k1*sdip*zim**2 + 2*k1k3*((-1 + ab1)*(-k1k3 + p1k1) + ab2*p2k1)*sdip
     &   *zim**2 - p1k1*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*sdip*zim**2 + ((-1 + ab1)*(-k1k3 + p1k1) + ab2*p2k1)
     &   *(ab1*(-k1k3 + p1k1) + ab2*p2k1)*sdip*zim**2 - 3*k1k3*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*tdip*zim**2 + 3
     &   *k1k3*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*udip*zim**2 - (ab1*(-k1k3 + p1k1) + ab2*p2k1)**2*udip*zim**2 
     &   + 2*k1k3**2*(-tdip + udip)*zim**2 - 2*M2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k2) + ab2*p2k2)
     &   *zim*zjm + 4*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*((-1 + ab1)*p1k3 + ab2
     &   *p2k3)*zim*zjm - k1k2*M2*sdip*zim*zjm + 2*k2k3*p1k1*sdip*zim*zjm + 2*k1k3*p1k2*sdip*zim*zjm - 2*k2k3
     &   *((-1 + ab1)*(-k1k3 + p1k1) + ab2*p2k1)*sdip*zim*zjm + p1k2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*sdip*zim*zjm 
     &   - 2*k1k3*((-1 + ab1)*(-k2k3 + p1k2) + ab2*p2k2)*sdip*zim*zjm - (ab1*(-k1k3 + p1k1) + ab2*p2k1)*((
     &   -1 + ab1)*(-k2k3 + p1k2) + ab2*p2k2)*sdip*zim*zjm + p1k1*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*sdip*zim*zjm 
     &   - ((-1 + ab1)*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*sdip*zim*zjm + 2*k1k2*((-1 
     &   + ab1)*p1k3 + ab2*p2k3)*sdip*zim*zjm + 3*k2k3*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*tdip*zim*zjm + 3*k1k3*(ab1
     &   *(-k2k3 + p1k2) + ab2*p2k2)*tdip*zim*zjm - k1k2*sdip*tdip*zim*zjm + 4*k1k3*k2k3*(tdip - udip)*zim*zjm 
     &   - 3*k2k3*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*udip*zim*zjm - 3*k1k3*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*udip
     &   *zim*zjm + 2*(ab1*(-k1k3 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*udip*zim*zjm + 2*k1k2*sdip
     &   *udip*zim*zjm + M2*(ab1*(-k2k3 + p1k2) + ab2*p2k2)**2*zjm**2 - 2*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*
     &   *2*((-1 + ab1)*p1k3 + ab2*p2k3)*zjm**2 - 2*k2k3*p1k2*sdip*zjm**2 + 2*k2k3*((-1 + ab1)*(-k2k3 + p1k2) 
     &   + ab2*p2k2)*sdip*zjm**2 - p1k2*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*sdip*zjm**2 + ((-1 + ab1)*(-k2k3 + p1k2) 
     &   + ab2*p2k2)*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*sdip*zjm**2 - 3*k2k3*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*tdip
     &   *zjm**2 + 3*k2k3*(ab1*(-k2k3 + p1k2) + ab2*p2k2)*udip*zjm**2 - (ab1*(-k2k3 + p1k2) + ab2*p2k2)**2*udip*zjm**2 
     &   + 2*k2k3**2*(-tdip + udip)*zjm**2)*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(8*sdip**2))

!************************************************   Mt^ast Vgg Mt   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_b = stsT2xx_M2Dipole_SpinTensor_D12_b +(-0.5*((-1 + Nc**2)*(M2
     &   *sdip - 2*(-1 + ab1)*p1k3*(M2 - udip) + 2*ab2*p2k3*(-M2 + udip))*(2*p1k1*zim - ab2*p2k1*zim - 2*p1k2
     &   *zjm + ab2*p2k2*zjm + ab1*(k1k3*zim - p1k1*zim - k2k3*zjm + p1k2*zjm))**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **4)/(Nc*sdip*(M2 - tdip)**2))

!************************************************   2 Re(Mt Vgg Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_b = stsT2xx_M2Dipole_SpinTensor_D12_b + 2*dble(-0.25d0*((-1 + Nc
     &   **2)*(M2**2 + 2*M2*sdip - sdip**2 + 2*(-1 + ab1)*p1k3*(M2 - tdip) + 2*ab2*p2k3*(M2 - tdip) - 2*M2*udip 
     &   + udip**2)*(-2*p1k1*zim + ab2*p2k1*zim + 2*p1k2*zjm - ab2*p2k2*zjm + ab1*(-(k1k3*zim) + p1k1*zim + k2k3
     &   *zjm - p1k2*zjm))*(ab2*p2k1*zim - ab2*p2k2*zjm - (-2 + ab1)*(k1k3*zim - p1k1*zim + (-k2k3 + p1k2)
     &   *zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(Nc*sdip*(M2 - tdip)*(M2 - udip)))

!************************************************   2 Re(Mt Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_b = stsT2xx_M2Dipole_SpinTensor_D12_b + 2*dble(((-1 + Nc**2)*(2
     &   *p1k1*zim - ab2*p2k1*zim - 2*p1k2*zjm + ab2*p2k2*zjm + ab1*(k1k3*zim - p1k1*zim - k2k3*zjm + p1k2*zjm))
     &   *(2*ab2*p1k3*p2k1*zim - 2*ab2**2*p2k1*p2k3*zim - p1k1*sdip*zim + ab2*p2k1*sdip*zim + k1k3*(-M2 + 2
     &   *ab1**2*p1k3 + sdip - ab1*(2*p1k3 - 2*ab2*p2k3 + sdip) + udip)*zim + k2k3*M2*zjm - 2*ab2*p1k3*p2k2*zjm 
     &   + 2*ab2**2*p2k2*p2k3*zjm - k2k3*sdip*zjm + p1k2*sdip*zjm - ab2*p2k2*sdip*zjm - k2k3*udip*zjm - 2*ab1
     &   **2*p1k3*(p1k1*zim + (k2k3 - p1k2)*zjm) + ab1*(p1k1*(2*p1k3 - 2*ab2*p2k3 + sdip)*zim + (k2k3 - p1k2)
     &   *(2*p1k3 + sdip)*zjm - 2*ab2*(p1k3*p2k1*zim - p1k3*p2k2*zjm + k2k3*p2k3*zjm - p1k2*p2k3*zjm)))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(2*Nc*sdip*(M2 - tdip)))

!************************************************   Mu^ast Vgg Mu   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_b = stsT2xx_M2Dipole_SpinTensor_D12_b +(((-1 + Nc**2)**2*(-M2
     &   **2 - tdip*udip + M2*(sdip + tdip + udip))*(-(ab2*p2k1*zim) + ab2*p2k2*zjm + (-2 + ab1)*(k1k3*zim - p1k1
     &   *zim + (-k2k3 + p1k2)*zjm))**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(2*Nc*sdip*(M2 - udip)**2))

!************************************************   2 Re(Mu Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_b = stsT2xx_M2Dipole_SpinTensor_D12_b + 2*dble(((2 - 3*Nc**2 
     &   + Nc**4)*((-2 + ab1)*k1k3*zim - (-2 + ab1)*p1k1*zim - ab2*p2k1*zim + 2*k2k3*zjm - ab1*k2k3*zjm - 2*p1k2
     &   *zjm + ab1*p1k2*zjm + ab2*p2k2*zjm)*(-(ab2*M2*p2k1*zim) + p1k1*sdip*zim + ab2*p2k1*udip*zim + k1k3*((
     &   -1 + ab1)*M2 + tdip - ab1*udip)*zim + k2k3*M2*zjm + ab2*M2*p2k2*zjm - p1k2*sdip*zjm - k2k3*tdip*zjm 
     &   - ab2*p2k2*udip*zjm - ab1*(M2 - udip)*(p1k1*zim + k2k3*zjm - p1k2*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(4*Nc*sdip*(M2 - udip)))

!************************************************   Mv^ast Vgg Mv   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D12_b = stsT2xx_M2Dipole_SpinTensor_D12_b +(((4 - 5*Nc**2 + Nc**4)
     &   *(-(k1k2*sdip*zim*zjm) + 2*ab1*(k1k3*zim - k2k3*zjm)*(k1k3*zim - p1k1*zim - k2k3*zjm + p1k2*zjm) - 2
     &   *ab2*(k1k3*zim - k2k3*zjm)*(p2k1*zim - p2k2*zjm))*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(4*Nc*sdip))



      end


      double precision function stsT2xx_M2Dipole_SpinTensor_D13_b(sdip,tdip,udip)
      implicit none
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"
#include "common/stsT2xx_GenCouplings.h"

      double precision sdip,tdip,udip,M2,zi,zM,zim,zjm
      double precision xija,Rdip,kappa,ab1,ab2

      M2=MSf2(iSf1,itt1,igen1)

      xija = (-k1k3 + p2k1 + p2k3)/(p2k1 + p2k3)

      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)

      zi =  p2k1/(p2k1 + p2k3)

      zM = (1-Rdip)/2d0

      ab1 = 1 - (M2 + udip)/(2*udip) + (udip - M2 + 2*M2*xija)/(2*Rdip*udip)
      ab2 = xija/Rdip

      zim = zi - zM
      zjm = 1 - zi - zM

      
!************************************************   Ms^ast Vgg Ms   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_b = ((GGG**2*Nc*(-1 + Nc**2)
     &   *(2*k1k2*M2**2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*zim**2 - 10*k1k2*M2*p1k2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*zim*
     &   *2 + 4*k1k2*p1k2**2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*zim**2 + 3*M2**2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)**2*zim**2 
     &   - 8*M2*p1k2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)**2*zim**2 + 4*p1k2**2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)**2*zim**2 + 8
     &   *k1k2**2*p1k2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*zim**2 + 12*k1k2*p1k2
     &   *(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))
     &   *zim**2 + 8*p1k2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)**2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2
     &   *(p1p2 - p2k2))*zim**2 - 8*k1k2**2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))
     &   **2*zim**2 - 8*k1k2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)
     &   *ab2*(p1p2 - p2k2))**2*zim**2 + 8*k1k2**2*M2*sdip*zim**2 - 4*k1k2*p1k1*p1k2*sdip*zim**2 + 4*k1k2*p1k2*((-1 + ab1)
     &   *(-k1k2 + p1k1) + ab2*p2k1)*sdip*zim**2 + 16*k1k2*M2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*sdip*zim**2 - 2*p1k1*p1k2
     &   *(ab1*(-k1k2 + p1k1) + ab2*p2k1)*sdip*zim**2 + 2*p1k2*((-1 + ab1)*(-k1k2 + p1k1) + ab2*p2k1)*(ab1*(-k1k2 + p1k1) 
     &   + ab2*p2k1)*sdip*zim**2 + 4*M2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)**2*sdip*zim**2 - 4*k1k2*p1k1*(ab2**2*M2 
     &   + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*sdip*zim**2 + 4*k1k2*((-1 + ab1)*(-k1k2 + p1k1) 
     &   + ab2*p2k1)*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*sdip*zim**2 - 2*p1k1
     &   *(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*sdip
     &   *zim**2 + 2*((-1 + ab1)*(-k1k2 + p1k1) + ab2*p2k1)*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab2**2*M2 + (-1 + ab1)
     &   *ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*sdip*zim**2 - 4*k1k2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*sdip
     &   **2*zim**2 - (ab1*(-k1k2 + p1k1) + ab2*p2k1)**2*sdip**2*zim**2 - 4*k1k2**2*sdip*(-2*M2 + sdip)*zim**2 - 4*k1k2
     &   **2*p1k2*(M2 - tdip)*zim**2 + 4*k1k2**2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))
     &   *(M2 - tdip)*zim**2 + 2*k1k2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (
     &   -1 + 2*ab1)*ab2*(p1p2 - p2k2))*(M2 - tdip)*zim**2 - 4*(ab1*(-k1k2 + p1k1) + ab2*p2k1)**2*(ab2**2*M2 + (-1 + ab1)
     &   *ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*(M2 - tdip)*zim**2 + 4*k1k2*p1k1*sdip*(M2 - tdip)*zim*
     &   *2 - 4*k1k2*((-1 + ab1)*(-k1k2 + p1k1) + ab2*p2k1)*sdip*(M2 - tdip)*zim**2 + 2*p1k1*(ab1*(-k1k2 + p1k1) + ab2*p2k1)
     &   *sdip*(M2 - tdip)*zim**2 - 2*((-1 + ab1)*(-k1k2 + p1k1) + ab2*p2k1)*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*sdip*(M2 
     &   - tdip)*zim**2 - 4*k1k2*M2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*tdip*zim**2 + 10*k1k2*p1k2*(ab1*(-k1k2 + p1k1) 
     &   + ab2*p2k1)*tdip*zim**2 - 6*M2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)**2*tdip*zim**2 + 8*p1k2*(ab1*(-k1k2 + p1k1) + ab2
     &   *p2k1)**2*tdip*zim**2 + 2*k1k2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*tdip**2*zim**2 + 3*(ab1*(-k1k2 + p1k1) + ab2*p2k1)
     &   **2*tdip**2*zim**2 - 4*k2k3*p1k2**2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*zim*zjm - 16*k1k2*k2k3*p1k2*(ab2**2*M2 
     &   + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*zim*zjm - 12*k2k3*p1k2*(ab1*(-k1k2 + p1k1) 
     &   + ab2*p2k1)*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*zim*zjm + 16*k1k2*k2k3
     &   *(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))**2*zim*zjm + 8*k2k3*(ab1*(-k1k2 + p1k1) 
     &   + ab2*p2k1)*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))**2*zim*zjm - 2*k1k2
     &   *M2**2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*zim*zjm + 10*k1k2*M2*p1k2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*zim*zjm - 4
     &   *k1k2*p1k2**2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*zim*zjm - 6*M2**2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 
     &   + p1k3) + ab2*p2k3)*zim*zjm + 16*M2*p1k2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*zim*zjm 
     &   - 8*p1k2**2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*zim*zjm - 12*k1k2*p1k2*(ab2**2*M2 
     &   + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*zim*zjm - 16
     &   *p1k2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))
     &   *(ab1*(-k2k3 + p1k3) + ab2*p2k3)*zim*zjm + 8*k1k2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)
     &   *ab2*(p1p2 - p2k2))**2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*zim*zjm - 16*k1k2*k2k3*M2*sdip*zim*zjm + 4*k2k3*p1k1*p1k2
     &   *sdip*zim*zjm - 2*k1k3*p1k2**2*sdip*zim*zjm + 4*k1k2*p1k2*p1k3*sdip*zim*zjm - 4*k2k3*p1k2*((-1 + ab1)*(-k1k2 
     &   + p1k1) + ab2*p2k1)*sdip*zim*zjm - 8*k2k3*M2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*sdip*zim*zjm + 2*p1k2*p1k3*(ab1
     &   *(-k1k2 + p1k1) + ab2*p2k1)*sdip*zim*zjm + 4*k2k3*p1k1*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)
     &   *ab2*(p1p2 - p2k2))*sdip*zim*zjm - 4*k1k3*p1k2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 
     &   - p2k2))*sdip*zim*zjm + 4*k1k2*p1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))
     &   *sdip*zim*zjm - 4*k2k3*((-1 + ab1)*(-k1k2 + p1k1) + ab2*p2k1)*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) 
     &   + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*sdip*zim*zjm + 2*p1k3*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab2**2*M2 + (-1 + ab1)
     &   *ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*sdip*zim*zjm - 2*k1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2
     &   *p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))**2*sdip*zim*zjm - 4*k1k2*p1k2*((-1 + ab1)*(-k2k3 + p1k3) + ab2*p2k3)
     &   *sdip*zim*zjm - 2*p1k2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*((-1 + ab1)*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm - 4
     &   *k1k2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*((-1 + ab1)*(-k2k3 + p1k3) + ab2
     &   *p2k3)*sdip*zim*zjm - 2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2
     &   *ab1)*ab2*(p1p2 - p2k2))*((-1 + ab1)*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm - 16*k1k2*M2*(ab1*(-k2k3 + p1k3) 
     &   + ab2*p2k3)*sdip*zim*zjm - 2*M2*p1k1*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm + 2*p1k1*p1k2*(ab1*(-k2k3 + p1k3) 
     &   + ab2*p2k3)*sdip*zim*zjm + 2*M2*((-1 + ab1)*(-k1k2 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip
     &   *zim*zjm - 2*p1k2*((-1 + ab1)*(-k1k2 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm - 8*M2
     &   *(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm + 2*p1k1*(ab2**2*M2 + (-1 + ab1)
     &   *ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm - 2*((-1 + ab1)
     &   *(-k1k2 + p1k1) + ab2*p2k1)*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*(ab1
     &   *(-k2k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm + 4*k1k2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip**2*zim*zjm + 2*(ab1*(
     &   -k1k2 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip**2*zim*zjm + 8*k1k2*k2k3*sdip*(-2*M2 + sdip)*zim
     &   *zjm + 4*k2k3*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*sdip*(-2*M2 + sdip)*zim*zjm + 8*k1k2*k2k3*p1k2*(M2 - tdip)*zim
     &   *zjm + 10*k2k3*p1k2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(M2 - tdip)*zim*zjm - 8*k1k2*k2k3*(ab2**2*M2 + (-1 + ab1)*ab1
     &   *(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*(M2 - tdip)*zim*zjm - 2*k2k3*(ab1*(-k1k2 + p1k1) + ab2*p2k1)
     &   *(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*(M2 - tdip)*zim*zjm - 2*k1k2*(ab2
     &   **2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*(M2 
     &   - tdip)*zim*zjm + 8*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2
     &   *(p1p2 - p2k2))*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*(M2 - tdip)*zim*zjm - 4*k2k3*p1k1*sdip*(M2 - tdip)*zim*zjm + 4
     &   *k1k3*p1k2*sdip*(M2 - tdip)*zim*zjm - 4*k1k2*p1k3*sdip*(M2 - tdip)*zim*zjm + 4*k2k3*((-1 + ab1)*(-k1k2 + p1k1) 
     &   + ab2*p2k1)*sdip*(M2 - tdip)*zim*zjm - 2*p1k3*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*sdip*(M2 - tdip)*zim*zjm + 4*k1k3
     &   *(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*sdip*(M2 - tdip)*zim*zjm + 4*k1k2
     &   *((-1 + ab1)*(-k2k3 + p1k3) + ab2*p2k3)*sdip*(M2 - tdip)*zim*zjm + 2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*((-1 
     &   + ab1)*(-k2k3 + p1k3) + ab2*p2k3)*sdip*(M2 - tdip)*zim*zjm - 2*k2k3*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(M2 - tdip)
     &   **2*zim*zjm - 2*k1k3*sdip*(M2 - tdip)**2*zim*zjm + 4*k1k2*M2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*tdip*zim*zjm - 10
     &   *k1k2*p1k2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*tdip*zim*zjm + 12*M2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k3) 
     &   + ab2*p2k3)*tdip*zim*zjm - 16*p1k2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*tdip
     &   *zim*zjm + 2*p1k1*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip*tdip*zim*zjm - 2*((-1 + ab1)*(-k1k2 + p1k1) + ab2*p2k1)
     &   *(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip*tdip*zim*zjm - 2*k1k2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*tdip**2*zim*zjm - 6
     &   *(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*tdip**2*zim*zjm + 8*k2k3**2*p1k2*(ab2**2*M2 
     &   + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*zjm**2 - 8*k2k3**2*(ab2**2*M2 + (-1 + ab1)*ab1
     &   *(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))**2*zjm**2 + 2*k2k3*M2**2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*zjm
     &   **2 - 10*k2k3*M2*p1k2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*zjm**2 + 4*k2k3*p1k2**2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*zjm
     &   **2 + 12*k2k3*p1k2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*(ab1*(-k2k3 
     &   + p1k3) + ab2*p2k3)*zjm**2 - 8*k2k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))
     &   **2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*zjm**2 + 3*M2**2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)**2*zjm**2 - 8*M2*p1k2*(ab1
     &   *(-k2k3 + p1k3) + ab2*p2k3)**2*zjm**2 + 4*p1k2**2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)**2*zjm**2 + 8*p1k2*(ab2**2
     &   *M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*(ab1*(-k2k3 + p1k3) + ab2*p2k3)**2*zjm**2 
     &   + 8*k2k3**2*M2*sdip*zjm**2 - 4*k2k3*p1k2*p1k3*sdip*zjm**2 - 4*k2k3*p1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) 
     &   + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*sdip*zjm**2 + 4*k2k3*p1k2*((-1 + ab1)*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zjm
     &   **2 + 4*k2k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*((-1 + ab1)*(-k2k3 
     &   + p1k3) + ab2*p2k3)*sdip*zjm**2 + 16*k2k3*M2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zjm**2 + 2*M2*p1k3*(ab1*(-k2k3 
     &   + p1k3) + ab2*p2k3)*sdip*zjm**2 - 2*p1k2*p1k3*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zjm**2 - 2*p1k3*(ab2**2*M2 
     &   + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zjm**2 
     &   - 2*M2*((-1 + ab1)*(-k2k3 + p1k3) + ab2*p2k3)*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zjm**2 + 2*p1k2*((-1 + ab1)
     &   *(-k2k3 + p1k3) + ab2*p2k3)*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zjm**2 + 2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 
     &   - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*((-1 + ab1)*(-k2k3 + p1k3) + ab2*p2k3)*(ab1*(-k2k3 + p1k3) + ab2*p2k3)
     &   *sdip*zjm**2 + 4*M2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)**2*sdip*zjm**2 - 4*k2k3*(ab1*(-k2k3 + p1k3) + ab2*p2k3)
     &   *sdip**2*zjm**2 - (ab1*(-k2k3 + p1k3) + ab2*p2k3)**2*sdip**2*zjm**2 - 4*k2k3**2*sdip*(-2*M2 + sdip)*zjm**2 - 4
     &   *k2k3**2*p1k2*(M2 - tdip)*zjm**2 + 4*k2k3**2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 
     &   - p2k2))*(M2 - tdip)*zjm**2 + 2*k2k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))
     &   *(ab1*(-k2k3 + p1k3) + ab2*p2k3)*(M2 - tdip)*zjm**2 - 4*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 
     &   + 2*ab1)*ab2*(p1p2 - p2k2))*(ab1*(-k2k3 + p1k3) + ab2*p2k3)**2*(M2 - tdip)*zjm**2 + 4*k2k3*p1k3*sdip*(M2 - tdip)
     &   *zjm**2 - 4*k2k3*((-1 + ab1)*(-k2k3 + p1k3) + ab2*p2k3)*sdip*(M2 - tdip)*zjm**2 - 4*k2k3*M2*(ab1*(-k2k3 + p1k3) 
     &   + ab2*p2k3)*tdip*zjm**2 + 10*k2k3*p1k2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*tdip*zjm**2 - 6*M2*(ab1*(-k2k3 + p1k3) 
     &   + ab2*p2k3)**2*tdip*zjm**2 + 8*p1k2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)**2*tdip*zjm**2 - 2*p1k3*(ab1*(-k2k3 + p1k3) 
     &   + ab2*p2k3)*sdip*tdip*zjm**2 + 2*((-1 + ab1)*(-k2k3 + p1k3) + ab2*p2k3)*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip*tdip*zjm**2 
     &   + 2*k2k3*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*tdip**2*zjm**2 + 3*(ab1*(-k2k3 + p1k3) + ab2*p2k3)**2*tdip**2*zjm**2)
     &   *SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(8*sdip**3))

!************************************************   2 Re(Ms Vgg Mt^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_b = stsT2xx_M2Dipole_SpinTensor_D13_b + 2*dble(0)

!************************************************   2 Re(Ms Vgg Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_b = stsT2xx_M2Dipole_SpinTensor_D13_b + 2*dble((GGG*Nc*(-1 + Nc
     &   **2)*((-2 + ab1)*k1k2*zim - (-2 + ab1)*p1k1*zim - ab2*p2k1*zim + 2*k2k3*zjm - ab1*k2k3*zjm - 2*p1k3*zjm + ab1
     &   *p1k3*zjm + ab2*p2k3*zjm)*(8*ab1**3*(M2 - 2*p1k2)*p1k2*(k1k2*zim - p1k1*zim + (-k2k3 + p1k3)*zjm) + 8*ab2**3*M2
     &   *p1k2*(-(p2k1*zim) + p2k3*zjm) + 2*(p1k1*p1k2*sdip*zim + p1k1*sdip*tdip*zim + k1k2*(3*M2*p1k2 - 4*M2*sdip + sdip
     &   **2 - 3*p1k2*tdip)*zim - p1k2*p1k3*sdip*zjm - k2k3*sdip**2*zjm + 3*k2k3*p1k2*tdip*zjm - p1k3*sdip*tdip*zjm 
     &   + M2*(-(p1k1*sdip*zim) - 3*k2k3*p1k2*zjm + 4*k2k3*sdip*zjm + p1k3*sdip*zjm)) - 2*ab2**2*(k1k2*M2*(M2 + 4*p1k2 
     &   - tdip)*zim - 4*p1k2*(p1p2*p2k1*zim - p2k1*p2k2*zim + k2k3*M2*zjm - p1p2*p2k3*zjm + p2k2*p2k3*zjm) + M2*(-(p1k1
     &   *sdip*zim) + p1k3*sdip*zjm + k2k3*(-M2 + tdip)*zjm)) + ab2*(-4*p1k2**2*p2k1*zim - 2*p1k1*p1p2*sdip*zim + 2*p1k1
     &   *p2k2*sdip*zim + p2k1*sdip**2*zim + 2*k1k2*(p1p2 - p2k2)*(M2 + 4*p1k2 - tdip)*zim - 8*p1k2*p2k1*tdip*zim - p2k1
     &   *tdip**2*zim - 8*k2k3*p1k2*p1p2*zjm + 8*k2k3*p1k2*p2k2*zjm + 4*p1k2**2*p2k3*zjm + 2*p1k3*p1p2*sdip*zjm - 2
     &   *p1k3*p2k2*sdip*zjm - p2k3*sdip**2*zjm + 2*k2k3*p1p2*tdip*zjm - 2*k2k3*p2k2*tdip*zjm + 8*p1k2*p2k3*tdip*zjm + p2k3
     &   *tdip**2*zjm + M2**2*(-(p2k1*zim) + p2k3*zjm) + 2*M2*(p2k1*(-2*sdip + tdip)*zim + (-(k2k3*p1p2) + k2k3*p2k2 
     &   + 2*p2k3*sdip - p2k3*tdip)*zjm + 4*p1k2*(p2k1*zim - p2k3*zjm))) - 2*ab1**2*(k1k2*(M2**2 + M2*(6*p1k2 - tdip) 
     &   + 2*p1k2*(-8*p1k2 - 4*ab2*p1p2 + 4*ab2*p2k2 + tdip))*zim - k2k3*M2**2*zjm + 2*p1k2*(p1k1*(4*p1k2 + 4*ab2*p1p2 
     &   - 4*ab2*p2k2 + sdip)*zim + 4*ab2*(k2k3 - p1k3)*(p1p2 - p2k2)*zjm - (p1k3*(4*p1k2 + sdip) + k2k3*(-8*p1k2 + tdip))
     &   *zjm - 4*ab2*p1k2*(p2k1*zim - p2k3*zjm)) + M2*(-(p1k1*(4*p1k2 + sdip)*zim) + (-6*k2k3*p1k2 + 4*p1k2*p1k3 + p1k3
     &   *sdip + k2k3*tdip)*zjm + 4*ab2*p1k2*(p2k1*zim - p2k3*zjm))) + ab1*(-4*p1k1*p1k2**2*zim + 8*ab2*p1k1*p1k2*p1p2
     &   *zim - 16*ab2*p1k2**2*p2k1*zim - 16*ab2**2*p1k2*p1p2*p2k1*zim - 8*ab2*p1k1*p1k2*p2k2*zim + 16*ab2**2*p1k2*p2k1
     &   *p2k2*zim + 4*p1k1*p1k2*sdip*zim + 4*ab2*p1k1*p1p2*sdip*zim - 4*ab2*p1k1*p2k2*sdip*zim + p1k1*sdip**2*zim 
     &   - 8*p1k1*p1k2*tdip*zim - p1k1*tdip**2*zim + k1k2*(3*M2**2 - 12*p1k2**2 - sdip**2 + 4*M2*((-1 + 2*ab2**2)*p1k2 
     &   - ab2*p1p2 + ab2*p2k2 + sdip - tdip) + 4*ab2*p1p2*tdip - 4*ab2*p2k2*tdip + tdip**2 + 12*p1k2*(-2*ab2*p1p2 + 2*ab2
     &   *p2k2 + tdip))*zim + 12*k2k3*p1k2**2*zjm + 4*p1k2**2*p1k3*zjm + 24*ab2*k2k3*p1k2*p1p2*zjm - 8*ab2*p1k2*p1k3
     &   *p1p2*zjm - 24*ab2*k2k3*p1k2*p2k2*zjm + 8*ab2*p1k2*p1k3*p2k2*zjm + 16*ab2*p1k2**2*p2k3*zjm + 16*ab2**2*p1k2*p1p2
     &   *p2k3*zjm - 16*ab2**2*p1k2*p2k2*p2k3*zjm - 4*p1k2*p1k3*sdip*zjm - 4*ab2*p1k3*p1p2*sdip*zjm + 4*ab2*p1k3*p2k2
     &   *sdip*zjm + k2k3*sdip**2*zjm - p1k3*sdip**2*zjm - 12*k2k3*p1k2*tdip*zjm + 8*p1k2*p1k3*tdip*zjm - 4*ab2*k2k3*p1p2
     &   *tdip*zjm + 4*ab2*k2k3*p2k2*tdip*zjm - k2k3*tdip**2*zjm + p1k3*tdip**2*zjm + M2**2*(-(p1k1*zim) + (-3*k2k3 + p1k3)
     &   *zjm) + M2*(2*p1k1*(-4*(-1 + ab2**2)*p1k2 - 3*sdip + tdip)*zim + 8*ab2**2*p1k2*(-k2k3 + p1k3)*zjm + 4*ab2*k2k3*(p1p2 
     &   - p2k2)*zjm - 2*p1k3*(4*p1k2 - 3*sdip + tdip)*zjm + 4*k2k3*(p1k2 - sdip + tdip)*zjm + 8*ab2*p1k2*(p2k1*zim - p2k3*zjm))))
     &   *SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(8*sdip**2*(M2 - udip)))

!************************************************   2 Re(Ms Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_b = stsT2xx_M2Dipole_SpinTensor_D13_b + 2*dble((GGG*Nc*(-1 + Nc
     &   **2)*(3*k1k2*M2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*zim**2 - 4*k1k2**2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) 
     &   + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*zim**2 - 6*k1k2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab2**2*M2 + (-1 + ab1)*ab1
     &   *(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*zim**2 - 2*k1k2*p1k1*sdip*zim**2 + 2*k1k2*((-1 + ab1)*(-k1k2 
     &   + p1k1) + ab2*p2k1)*sdip*zim**2 - p1k1*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*sdip*zim**2 + ((-1 + ab1)*(-k1k2 + p1k1) 
     &   + ab2*p2k1)*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*sdip*zim**2 + 2*k1k2**2*(M2 - tdip)*zim**2 - (ab1*(-k1k2 + p1k1) 
     &   + ab2*p2k1)**2*(M2 - 2*p1k2 - tdip)*zim**2 - 3*k1k2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*tdip*zim**2 - 3*k2k3*M2*(ab1
     &   *(-k1k2 + p1k1) + ab2*p2k1)*zim*zjm + 8*k1k2*k2k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2
     &   *(p1p2 - p2k2))*zim*zjm + 6*k2k3*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) 
     &   + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*zim*zjm - 3*k1k2*M2*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*zim*zjm + 2*M2*(ab1*(-k1k2 
     &   + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*zim*zjm - 4*p1k2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab1*(
     &   -k2k3 + p1k3) + ab2*p2k3)*zim*zjm + 6*k1k2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 
     &   - p2k2))*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*zim*zjm + 2*k2k3*p1k1*sdip*zim*zjm - 2*k1k3*p1k2*sdip*zim*zjm + 2*k1k2
     &   *p1k3*sdip*zim*zjm - 2*k2k3*((-1 + ab1)*(-k1k2 + p1k1) + ab2*p2k1)*sdip*zim*zjm + p1k3*(ab1*(-k1k2 + p1k1) + ab2
     &   *p2k1)*sdip*zim*zjm - 2*k1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*sdip
     &   *zim*zjm - 2*k1k2*((-1 + ab1)*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm - (ab1*(-k1k2 + p1k1) + ab2*p2k1)*((-1 
     &   + ab1)*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm + p1k1*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm - ((-1 + ab1)
     &   *(-k1k2 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm - 4*k1k2*k2k3*(M2 - tdip)*zim*zjm 
     &   + 2*k1k3*sdip*(M2 - tdip)*zim*zjm + 3*k2k3*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*tdip*zim*zjm + 3*k1k2*(ab1*(-k2k3 + p1k3) 
     &   + ab2*p2k3)*tdip*zim*zjm - 2*(ab1*(-k1k2 + p1k1) + ab2*p2k1)*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*tdip*zim*zjm 
     &   - 4*k2k3**2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*zjm**2 + 3*k2k3*M2*(ab1
     &   *(-k2k3 + p1k3) + ab2*p2k3)*zjm**2 - 6*k2k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)*ab2*(p1p2 
     &   - p2k2))*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*zjm**2 - 2*k2k3*p1k3*sdip*zjm**2 + 2*k2k3*((-1 + ab1)*(-k2k3 + p1k3) 
     &   + ab2*p2k3)*sdip*zjm**2 - p1k3*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zjm**2 + ((-1 + ab1)*(-k2k3 + p1k3) + ab2
     &   *p2k3)*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*sdip*zjm**2 + 2*k2k3**2*(M2 - tdip)*zjm**2 - (ab1*(-k2k3 + p1k3) 
     &   + ab2*p2k3)**2*(M2 - 2*p1k2 - tdip)*zjm**2 - 3*k2k3*(ab1*(-k2k3 + p1k3) + ab2*p2k3)*tdip*zjm**2)
     &   *SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(8*sdip**2))

!************************************************   Mt^ast Vgg Mt   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_b = stsT2xx_M2Dipole_SpinTensor_D13_b +(((-1 + Nc**2)*(4*ab1*
     &   *3*(M2 - 2*p1k2)*p1k2 + 4*ab2**3*M2*p2k2 + (-M2 + p1k2)*sdip + ab2*(2*p1p2 - p2k2)*(2*p1k2 + sdip) - ab2**2*(4
     &   *M2*p1k2 + 4*p1p2*p2k2 - 2*p2k2**2 + M2*sdip) + ab1**2*(2*p1k2*(7*p1k2 + 4*ab2*p1p2 - 8*ab2*p2k2 + sdip) - M2*(8
     &   *p1k2 - 4*ab2*p2k2 + sdip)) + ab1*(-6*p1k2**2 + p1k2*(-12*ab2*p1p2 + 16*ab2*p2k2 - 3*sdip) + 2*ab2*(p1p2 - p2k2)
     &   *(4*ab2*p2k2 - sdip) + 2*M2*(2*(1 + ab2**2)*p1k2 - 2*ab2*p2k2 + sdip)))*(2*p1k1*zim - ab2*p2k1*zim - 2*p1k3
     &   *zjm + ab2*p2k3*zjm + ab1*(k1k2*zim - p1k1*zim - k2k3*zjm + p1k3*zjm))**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(2*Nc*sdip*(M2 - tdip)**2))

!************************************************   2 Re(Mt Vgg Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_b = stsT2xx_M2Dipole_SpinTensor_D13_b + 2*dble(-0.25d0*((-1 + Nc
     &   **2)*(2*ab1**2*(4*M2 - 7*p1k2)*p1k2 + ab2**2*(4*M2*p1k2 + 4*p1p2*p2k2 - 2*p2k2**2) + 2*(M2 - p1k2)*sdip - ab2
     &   *(2*p1p2 - p2k2)*(4*p1k2 + sdip) - 2*ab1*M2*(4*p1k2 - 2*ab2*p2k2 + sdip) + ab1*p1k2*(12*p1k2 + 4*ab2*(3*p1p2 - 4
     &   *p2k2) + 3*sdip))*((-2 + ab1)*k1k2*zim - (-2 + ab1)*p1k1*zim - ab2*p2k1*zim + 2*k2k3*zjm - ab1*k2k3*zjm - 2*p1k3
     &   *zjm + ab1*p1k3*zjm + ab2*p2k3*zjm)*(2*p1k1*zim - ab2*p2k1*zim - 2*p1k3*zjm + ab2*p2k3*zjm + ab1*(k1k2*zim 
     &   - p1k1*zim - k2k3*zjm + p1k3*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(Nc*sdip*(M2 - tdip)*(M2 - udip)))

!************************************************   2 Re(Mt Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_b = stsT2xx_M2Dipole_SpinTensor_D13_b + 2*dble(((-1 + Nc**2)*(ab1
     &   *M2*(-k1k2 + p1k1)*zim + ab2*M2*p2k1*zim + 2*k1k2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) + (-1 + 2*ab1)
     &   *ab2*(p1p2 - p2k2))*zim + ab1*(k1k2 - p1k1)*sdip*zim + (-k1k2 + p1k1)*sdip*zim - ab2*p2k1*sdip*zim + ab1*(k1k2 
     &   - p1k1)*tdip*zim - ab2*p2k1*tdip*zim + ab1*M2*(k2k3 - p1k3)*zjm - 2*k2k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k2) 
     &   + (-1 + 2*ab1)*ab2*(p1p2 - p2k2))*zjm - ab2*M2*p2k3*zjm + (k2k3 - p1k3)*sdip*zjm + ab1*(-k2k3 + p1k3)*sdip
     &   *zjm + ab2*p2k3*sdip*zjm + ab1*(-k2k3 + p1k3)*tdip*zjm + ab2*p2k3*tdip*zjm)*(-2*p1k1*zim + ab2*p2k1*zim + 2*p1k3
     &   *zjm - ab2*p2k3*zjm + ab1*(-(k1k2*zim) + p1k1*zim + k2k3*zjm - p1k3*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(2*Nc*sdip*(M2 - tdip)))

!************************************************   Mu^ast Vgg Mu   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_b = stsT2xx_M2Dipole_SpinTensor_D13_b +(((-1 + Nc**2)**2*(ab1
     &   *(-4*M2*p1k2 + 6*p1k2**2) + 2*ab2*p1k2*(-2*p1p2 + p2k2) + (M2 - p1k2)*sdip)*(-(ab2*p2k1*zim) + ab2*p2k3*zjm + (
     &   -2 + ab1)*(k1k2*zim - p1k1*zim + (-k2k3 + p1k3)*zjm))**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(2*Nc*sdip*(M2 - udip)**2))

!************************************************   2 Re(Mu Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_b = stsT2xx_M2Dipole_SpinTensor_D13_b + 2*dble(((2 - 3*Nc**2 
     &   + Nc**4)*(2*ab2*p1k2*p2k1*zim - p1k1*sdip*zim + k1k2*(M2 - 2*ab1*p1k2 - tdip)*zim - k2k3*M2*zjm - 2*ab2*p1k2*p2k3
     &   *zjm + p1k3*sdip*zjm + k2k3*tdip*zjm + 2*ab1*p1k2*(p1k1*zim + k2k3*zjm - p1k3*zjm))*(ab2*p2k1*zim - ab2*p2k3
     &   *zjm - (-2 + ab1)*(k1k2*zim - p1k1*zim + (-k2k3 + p1k3)*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(4.
     &   *Nc*sdip*(M2 - udip)))

!************************************************   Mv^ast Vgg Mv   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D13_b = stsT2xx_M2Dipole_SpinTensor_D13_b +(((4 - 5*Nc**2 + Nc**4)
     &   *(-(k1k3*sdip*zim*zjm) + 2*ab1*(k1k2*zim - k2k3*zjm)*(k1k2*zim - p1k1*zim - k2k3*zjm + p1k3*zjm) - 2*ab2*(k1k2
     &   *zim - k2k3*zjm)*(p2k1*zim - p2k3*zjm))*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(4*Nc*sdip))



      end


      double precision function stsT2xx_M2Dipole_SpinTensor_D23_b(sdip,tdip,udip)
      implicit none
#include "common/stsT2xx_Kinematics.h"
#include "common/stsT2xx_Model.h"
#include "common/stsT2xx_GenCouplings.h"

      double precision sdip,tdip,udip,M2,zi,zM,zim,zjm
      double precision xija,Rdip,kappa,ab1,ab2

      M2=MSf2(iSf1,itt1,igen1)  

      xija = (-k2k3 + p2k2 + p2k3)/(p2k2 + p2k3)

      Rdip = dsqrt(-4*M2*udip*xija**2 + (-M2 + udip + 2*M2*xija)**2)/kappa(udip,0d0,M2)

      zi =  p2k2/(p2k2 + p2k3)

      zM = (1-Rdip)/2d0

      ab1 = 1 - (M2 + udip)/(2*udip) + (udip - M2 + 2*M2*xija)/(2*Rdip*udip)
      ab2 = xija/Rdip

      zim = zi - zM
      zjm = 1 - zi - zM

      
!************************************************   Ms^ast Vgg Ms   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_b = ((GGG**2*Nc*(-1 + Nc**2)
     &   *(8*k1k2**2*p1k1*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*zim**2 - 8*k1k2
     &   **2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))**2*zim**2 + 2*k1k2*M2**2*(ab1*(
     &   -k1k2 + p1k2) + ab2*p2k2)*zim**2 - 10*k1k2*M2*p1k1*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*zim**2 + 4*k1k2*p1k1**2*(ab1
     &   *(-k1k2 + p1k2) + ab2*p2k2)*zim**2 + 12*k1k2*p1k1*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2
     &   *(p1p2 - p2k1))*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*zim**2 - 8*k1k2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (
     &   -1 + 2*ab1)*ab2*(p1p2 - p2k1))**2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*zim**2 + 3*M2**2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)
     &   **2*zim**2 - 8*M2*p1k1*(ab1*(-k1k2 + p1k2) + ab2*p2k2)**2*zim**2 + 4*p1k1**2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)
     &   **2*zim**2 + 8*p1k1*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*(ab1*(-k1k2 
     &   + p1k2) + ab2*p2k2)**2*zim**2 + 8*k1k2**2*M2*sdip*zim**2 - 4*k1k2*p1k1*p1k2*sdip*zim**2 - 4*k1k2*p1k2*(ab2**2
     &   *M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*sdip*zim**2 + 4*k1k2*p1k1*((-1 + ab1)*(-k1k2 
     &   + p1k2) + ab2*p2k2)*sdip*zim**2 + 4*k1k2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 
     &   - p2k1))*((-1 + ab1)*(-k1k2 + p1k2) + ab2*p2k2)*sdip*zim**2 + 16*k1k2*M2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*sdip
     &   *zim**2 - 2*p1k1*p1k2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*sdip*zim**2 - 2*p1k2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2
     &   *p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*sdip*zim**2 + 2*p1k1*((-1 + ab1)*(-k1k2 
     &   + p1k2) + ab2*p2k2)*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*sdip*zim**2 + 2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) 
     &   + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*((-1 + ab1)*(-k1k2 + p1k2) + ab2*p2k2)*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*sdip
     &   *zim**2 + 4*M2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)**2*sdip*zim**2 - 4*k1k2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*sdip
     &   **2*zim**2 - (ab1*(-k1k2 + p1k2) + ab2*p2k2)**2*sdip**2*zim**2 - 4*k1k2**2*sdip*(-2*M2 + sdip)*zim**2 - 4*k1k2
     &   **2*p1k1*(M2 - tdip)*zim**2 + 4*k1k2**2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))
     &   *(M2 - tdip)*zim**2 + 2*k1k2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*(ab1
     &   *(-k1k2 + p1k2) + ab2*p2k2)*(M2 - tdip)*zim**2 - 4*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)
     &   *ab2*(p1p2 - p2k1))*(ab1*(-k1k2 + p1k2) + ab2*p2k2)**2*(M2 - tdip)*zim**2 + 4*k1k2*p1k2*sdip*(M2 - tdip)*zim*
     &   *2 - 4*k1k2*((-1 + ab1)*(-k1k2 + p1k2) + ab2*p2k2)*sdip*(M2 - tdip)*zim**2 + 2*p1k2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)
     &   *sdip*(M2 - tdip)*zim**2 - 2*((-1 + ab1)*(-k1k2 + p1k2) + ab2*p2k2)*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*sdip*(M2 
     &   - tdip)*zim**2 - 4*k1k2*M2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*tdip*zim**2 + 10*k1k2*p1k1*(ab1*(-k1k2 + p1k2) 
     &   + ab2*p2k2)*tdip*zim**2 - 6*M2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)**2*tdip*zim**2 + 8*p1k1*(ab1*(-k1k2 + p1k2) + ab2
     &   *p2k2)**2*tdip*zim**2 + 2*k1k2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*tdip**2*zim**2 + 3*(ab1*(-k1k2 + p1k2) + ab2*p2k2)
     &   **2*tdip**2*zim**2 - 16*k1k2*k1k3*p1k1*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 
     &   - p2k1))*zim*zjm + 16*k1k2*k1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))**2
     &   *zim*zjm - 4*k1k3*p1k1**2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*zim*zjm - 12*k1k3*p1k1*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 
     &   - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*zim*zjm + 8*k1k3*(ab2**2*M2 + (
     &   -1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))**2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*zim*zjm - 2*k1k2
     &   *M2**2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zim*zjm + 10*k1k2*M2*p1k1*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zim*zjm - 4
     &   *k1k2*p1k1**2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zim*zjm - 12*k1k2*p1k1*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) 
     &   + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zim*zjm + 8*k1k2*(ab2**2*M2 + (-1 + ab1)*ab1
     &   *(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))**2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zim*zjm - 6*M2**2*(ab1*(
     &   -k1k2 + p1k2) + ab2*p2k2)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zim*zjm + 16*M2*p1k1*(ab1*(-k1k2 + p1k2) + ab2*p2k2)
     &   *(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zim*zjm - 8*p1k1**2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*(ab1*(-k1k3 + p1k3) + ab2
     &   *p2k3)*zim*zjm - 16*p1k1*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*(ab1*(-k1k2 
     &   + p1k2) + ab2*p2k2)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zim*zjm - 16*k1k2*k1k3*M2*sdip*zim*zjm - 2*k2k3*p1k1**2
     &   *sdip*zim*zjm + 4*k1k3*p1k1*p1k2*sdip*zim*zjm + 4*k1k2*p1k1*p1k3*sdip*zim*zjm - 4*k2k3*p1k1*(ab2**2*M2 + (-1 
     &   + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*sdip*zim*zjm + 4*k1k3*p1k2*(ab2**2*M2 + (-1 + ab1)*ab1
     &   *(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*sdip*zim*zjm + 4*k1k2*p1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 
     &   - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*sdip*zim*zjm - 2*k2k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (
     &   -1 + 2*ab1)*ab2*(p1p2 - p2k1))**2*sdip*zim*zjm - 4*k1k3*p1k1*((-1 + ab1)*(-k1k2 + p1k2) + ab2*p2k2)*sdip*zim
     &   *zjm - 4*k1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*((-1 + ab1)*(-k1k2 + p1k2) 
     &   + ab2*p2k2)*sdip*zim*zjm - 8*k1k3*M2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*sdip*zim*zjm + 2*p1k1*p1k3*(ab1*(-k1k2 
     &   + p1k2) + ab2*p2k2)*sdip*zim*zjm + 2*p1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 
     &   - p2k1))*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*sdip*zim*zjm - 4*k1k2*p1k1*((-1 + ab1)*(-k1k3 + p1k3) + ab2*p2k3)
     &   *sdip*zim*zjm - 4*k1k2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*((-1 + ab1)
     &   *(-k1k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm - 2*p1k1*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*((-1 + ab1)*(-k1k3 + p1k3) + ab2
     &   *p2k3)*sdip*zim*zjm - 2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*(ab1*(
     &   -k1k2 + p1k2) + ab2*p2k2)*((-1 + ab1)*(-k1k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm - 16*k1k2*M2*(ab1*(-k1k3 + p1k3) 
     &   + ab2*p2k3)*sdip*zim*zjm - 2*M2*p1k2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm + 2*p1k1*p1k2*(ab1*(-k1k3 + p1k3) 
     &   + ab2*p2k3)*sdip*zim*zjm + 2*p1k2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))
     &   *(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm + 2*M2*((-1 + ab1)*(-k1k2 + p1k2) + ab2*p2k2)*(ab1*(-k1k3 
     &   + p1k3) + ab2*p2k3)*sdip*zim*zjm - 2*p1k1*((-1 + ab1)*(-k1k2 + p1k2) + ab2*p2k2)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)
     &   *sdip*zim*zjm - 2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*((-1 + ab1)*(-k1k2 
     &   + p1k2) + ab2*p2k2)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm - 8*M2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*(ab1
     &   *(-k1k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm + 4*k1k2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip**2*zim*zjm + 2*(ab1*(
     &   -k1k2 + p1k2) + ab2*p2k2)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip**2*zim*zjm + 8*k1k2*k1k3*sdip*(-2*M2 + sdip)*zim
     &   *zjm + 4*k1k3*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*sdip*(-2*M2 + sdip)*zim*zjm + 8*k1k2*k1k3*p1k1*(M2 - tdip)*zim
     &   *zjm - 8*k1k2*k1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*(M2 - tdip)*zim
     &   *zjm + 10*k1k3*p1k1*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*(M2 - tdip)*zim*zjm - 2*k1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 
     &   - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*(M2 - tdip)*zim*zjm - 2*k1k2*(ab2
     &   **2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*(M2 
     &   - tdip)*zim*zjm + 8*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*(ab1*(-k1k2 + p1k2) 
     &   + ab2*p2k2)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*(M2 - tdip)*zim*zjm + 4*k2k3*p1k1*sdip*(M2 - tdip)*zim*zjm - 4
     &   *k1k3*p1k2*sdip*(M2 - tdip)*zim*zjm - 4*k1k2*p1k3*sdip*(M2 - tdip)*zim*zjm + 4*k2k3*(ab2**2*M2 + (-1 + ab1)*ab1
     &   *(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*sdip*(M2 - tdip)*zim*zjm + 4*k1k3*((-1 + ab1)*(-k1k2 + p1k2) 
     &   + ab2*p2k2)*sdip*(M2 - tdip)*zim*zjm - 2*p1k3*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*sdip*(M2 - tdip)*zim*zjm + 4*k1k2
     &   *((-1 + ab1)*(-k1k3 + p1k3) + ab2*p2k3)*sdip*(M2 - tdip)*zim*zjm + 2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*((-1 
     &   + ab1)*(-k1k3 + p1k3) + ab2*p2k3)*sdip*(M2 - tdip)*zim*zjm - 2*k1k3*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*(M2 - tdip)
     &   **2*zim*zjm - 2*k2k3*sdip*(M2 - tdip)**2*zim*zjm + 4*k1k2*M2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*tdip*zim*zjm - 10
     &   *k1k2*p1k1*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*tdip*zim*zjm + 12*M2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*(ab1*(-k1k3 + p1k3) 
     &   + ab2*p2k3)*tdip*zim*zjm - 16*p1k1*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*tdip
     &   *zim*zjm + 2*p1k2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip*tdip*zim*zjm - 2*((-1 + ab1)*(-k1k2 + p1k2) + ab2*p2k2)
     &   *(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip*tdip*zim*zjm - 2*k1k2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*tdip**2*zim*zjm - 6
     &   *(ab1*(-k1k2 + p1k2) + ab2*p2k2)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*tdip**2*zim*zjm + 8*k1k3**2*p1k1*(ab2**2*M2 
     &   + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*zjm**2 - 8*k1k3**2*(ab2**2*M2 + (-1 + ab1)*ab1
     &   *(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))**2*zjm**2 + 2*k1k3*M2**2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zjm
     &   **2 - 10*k1k3*M2*p1k1*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zjm**2 + 4*k1k3*p1k1**2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zjm
     &   **2 + 12*k1k3*p1k1*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*(ab1*(-k1k3 
     &   + p1k3) + ab2*p2k3)*zjm**2 - 8*k1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))
     &   **2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zjm**2 + 3*M2**2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)**2*zjm**2 - 8*M2*p1k1*(ab1
     &   *(-k1k3 + p1k3) + ab2*p2k3)**2*zjm**2 + 4*p1k1**2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)**2*zjm**2 + 8*p1k1*(ab2**2
     &   *M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*(ab1*(-k1k3 + p1k3) + ab2*p2k3)**2*zjm**2 
     &   + 8*k1k3**2*M2*sdip*zjm**2 - 4*k1k3*p1k1*p1k3*sdip*zjm**2 - 4*k1k3*p1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) 
     &   + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*sdip*zjm**2 + 4*k1k3*p1k1*((-1 + ab1)*(-k1k3 + p1k3) + ab2*p2k3)*sdip*zjm
     &   **2 + 4*k1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*((-1 + ab1)*(-k1k3 
     &   + p1k3) + ab2*p2k3)*sdip*zjm**2 + 16*k1k3*M2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip*zjm**2 + 2*M2*p1k3*(ab1*(-k1k3 
     &   + p1k3) + ab2*p2k3)*sdip*zjm**2 - 2*p1k1*p1k3*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip*zjm**2 - 2*p1k3*(ab2**2*M2 
     &   + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip*zjm**2 
     &   - 2*M2*((-1 + ab1)*(-k1k3 + p1k3) + ab2*p2k3)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip*zjm**2 + 2*p1k1*((-1 + ab1)
     &   *(-k1k3 + p1k3) + ab2*p2k3)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip*zjm**2 + 2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 
     &   - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*((-1 + ab1)*(-k1k3 + p1k3) + ab2*p2k3)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)
     &   *sdip*zjm**2 + 4*M2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)**2*sdip*zjm**2 - 4*k1k3*(ab1*(-k1k3 + p1k3) + ab2*p2k3)
     &   *sdip**2*zjm**2 - (ab1*(-k1k3 + p1k3) + ab2*p2k3)**2*sdip**2*zjm**2 - 4*k1k3**2*sdip*(-2*M2 + sdip)*zjm**2 - 4
     &   *k1k3**2*p1k1*(M2 - tdip)*zjm**2 + 4*k1k3**2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 
     &   - p2k1))*(M2 - tdip)*zjm**2 + 2*k1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))
     &   *(ab1*(-k1k3 + p1k3) + ab2*p2k3)*(M2 - tdip)*zjm**2 - 4*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)
     &   *ab2*(p1p2 - p2k1))*(ab1*(-k1k3 + p1k3) + ab2*p2k3)**2*(M2 - tdip)*zjm**2 + 4*k1k3*p1k3*sdip*(M2 - tdip)*zjm**2 
     &   - 4*k1k3*((-1 + ab1)*(-k1k3 + p1k3) + ab2*p2k3)*sdip*(M2 - tdip)*zjm**2 - 4*k1k3*M2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)
     &   *tdip*zjm**2 + 10*k1k3*p1k1*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*tdip*zjm**2 - 6*M2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)**2
     &   *tdip*zjm**2 + 8*p1k1*(ab1*(-k1k3 + p1k3) + ab2*p2k3)**2*tdip*zjm**2 - 2*p1k3*(ab1*(-k1k3 + p1k3) + ab2*p2k3)
     &   *sdip*tdip*zjm**2 + 2*((-1 + ab1)*(-k1k3 + p1k3) + ab2*p2k3)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip*tdip*zjm**2 
     &   + 2*k1k3*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*tdip**2*zjm**2 + 3*(ab1*(-k1k3 + p1k3) + ab2*p2k3)**2*tdip**2*zjm**2)
     &   *SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(8*sdip**3))

!************************************************   2 Re(Ms Vgg Mt^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_b = stsT2xx_M2Dipole_SpinTensor_D23_b + 2*dble(0)

!************************************************   2 Re(Ms Vgg Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_b = stsT2xx_M2Dipole_SpinTensor_D23_b + 2*dble((GGG*Nc*(-1 + Nc
     &   **2)*((-2 + ab1)*k1k2*zim - (-2 + ab1)*p1k2*zim - ab2*p2k2*zim + 2*k1k3*zjm - ab1*k1k3*zjm - 2*p1k3*zjm + ab1
     &   *p1k3*zjm + ab2*p2k3*zjm)*(8*ab1**3*(M2 - 2*p1k1)*p1k1*(k1k2*zim - p1k2*zim + (-k1k3 + p1k3)*zjm) + 8*ab2**3*M2
     &   *p1k1*(-(p2k2*zim) + p2k3*zjm) + 2*(p1k1*p1k2*sdip*zim + p1k2*sdip*tdip*zim + k1k2*(3*M2*p1k1 - 4*M2*sdip + sdip
     &   **2 - 3*p1k1*tdip)*zim - p1k1*p1k3*sdip*zjm - k1k3*sdip**2*zjm + 3*k1k3*p1k1*tdip*zjm - p1k3*sdip*tdip*zjm 
     &   + M2*(-(p1k2*sdip*zim) - 3*k1k3*p1k1*zjm + 4*k1k3*sdip*zjm + p1k3*sdip*zjm)) - 2*ab2**2*(k1k2*M2*(M2 + 4*p1k1 
     &   - tdip)*zim - 4*p1k1*(p1p2*p2k2*zim - p2k1*p2k2*zim + k1k3*M2*zjm - p1p2*p2k3*zjm + p2k1*p2k3*zjm) + M2*(-(p1k2
     &   *sdip*zim) + p1k3*sdip*zjm + k1k3*(-M2 + tdip)*zjm)) + ab2*(-4*p1k1**2*p2k2*zim - 2*p1k2*p1p2*sdip*zim + 2*p1k2
     &   *p2k1*sdip*zim + p2k2*sdip**2*zim + 2*k1k2*(p1p2 - p2k1)*(M2 + 4*p1k1 - tdip)*zim - 8*p1k1*p2k2*tdip*zim - p2k2
     &   *tdip**2*zim - 8*k1k3*p1k1*p1p2*zjm + 8*k1k3*p1k1*p2k1*zjm + 4*p1k1**2*p2k3*zjm + 2*p1k3*p1p2*sdip*zjm - 2
     &   *p1k3*p2k1*sdip*zjm - p2k3*sdip**2*zjm + 2*k1k3*p1p2*tdip*zjm - 2*k1k3*p2k1*tdip*zjm + 8*p1k1*p2k3*tdip*zjm + p2k3
     &   *tdip**2*zjm + M2**2*(-(p2k2*zim) + p2k3*zjm) + 2*M2*(p2k2*(-2*sdip + tdip)*zim + (-(k1k3*p1p2) + k1k3*p2k1 
     &   + 2*p2k3*sdip - p2k3*tdip)*zjm + 4*p1k1*(p2k2*zim - p2k3*zjm))) - 2*ab1**2*(k1k2*(M2**2 + M2*(6*p1k1 - tdip) 
     &   + 2*p1k1*(-8*p1k1 - 4*ab2*p1p2 + 4*ab2*p2k1 + tdip))*zim - k1k3*M2**2*zjm + M2*(-(p1k2*sdip*zim) + p1k3*sdip*zjm 
     &   + k1k3*tdip*zjm + p1k1*(-4*p1k2*zim + 4*ab2*p2k2*zim - 6*k1k3*zjm + 4*p1k3*zjm - 4*ab2*p2k3*zjm)) + 2*p1k1
     &   *(p1k2*sdip*zim - p1k3*sdip*zjm - k1k3*tdip*zjm + 4*ab2*(p1p2 - p2k1)*(p1k2*zim + k1k3*zjm - p1k3*zjm) + 4*p1k1
     &   *(p1k2*zim - ab2*p2k2*zim + 2*k1k3*zjm - p1k3*zjm + ab2*p2k3*zjm))) + ab1*(-4*p1k1**2*p1k2*zim + 8*ab2*p1k1*p1k2
     &   *p1p2*zim - 8*ab2*p1k1*p1k2*p2k1*zim - 16*ab2*p1k1**2*p2k2*zim - 16*ab2**2*p1k1*p1p2*p2k2*zim + 16*ab2**2*p1k1
     &   *p2k1*p2k2*zim + 4*p1k1*p1k2*sdip*zim + 4*ab2*p1k2*p1p2*sdip*zim - 4*ab2*p1k2*p2k1*sdip*zim + p1k2*sdip**2*zim 
     &   - 8*p1k1*p1k2*tdip*zim - p1k2*tdip**2*zim + k1k2*(3*M2**2 - 12*p1k1**2 - sdip**2 + 4*M2*((-1 + 2*ab2**2)*p1k1 
     &   - ab2*p1p2 + ab2*p2k1 + sdip - tdip) + 4*ab2*p1p2*tdip - 4*ab2*p2k1*tdip + tdip**2 + 12*p1k1*(-2*ab2*p1p2 
     &   + 2*ab2*p2k1 + tdip))*zim + 12*k1k3*p1k1**2*zjm + 4*p1k1**2*p1k3*zjm + 24*ab2*k1k3*p1k1*p1p2*zjm - 8*ab2*p1k1*p1k3
     &   *p1p2*zjm - 24*ab2*k1k3*p1k1*p2k1*zjm + 8*ab2*p1k1*p1k3*p2k1*zjm + 16*ab2*p1k1**2*p2k3*zjm + 16*ab2**2*p1k1
     &   *p1p2*p2k3*zjm - 16*ab2**2*p1k1*p2k1*p2k3*zjm - 4*p1k1*p1k3*sdip*zjm - 4*ab2*p1k3*p1p2*sdip*zjm + 4*ab2*p1k3*p2k1
     &   *sdip*zjm + k1k3*sdip**2*zjm - p1k3*sdip**2*zjm - 12*k1k3*p1k1*tdip*zjm + 8*p1k1*p1k3*tdip*zjm - 4*ab2*k1k3
     &   *p1p2*tdip*zjm + 4*ab2*k1k3*p2k1*tdip*zjm - k1k3*tdip**2*zjm + p1k3*tdip**2*zjm + M2**2*(-(p1k2*zim) + (-3*k1k3 
     &   + p1k3)*zjm) + M2*(2*p1k2*(-3*sdip + tdip)*zim + 2*(2*ab2*k1k3*(p1p2 - p2k1) - 2*k1k3*sdip + 3*p1k3*sdip + 2*k1k3*tdip
     &   - p1k3*tdip)*zjm + p1k1*(-8*(-1 + ab2**2)*p1k2*zim + 4*(k1k3 - 2*p1k3)*zjm + 8*ab2**2*(-k1k3 + p1k3)*zjm + 8*ab2
     &  *(p2k2*zim - p2k3*zjm)))))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**3)/(8*sdip**2*(M2 - udip)))

!************************************************   2 Re(Ms Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_b = stsT2xx_M2Dipole_SpinTensor_D23_b + 2*dble((GGG*Nc*(-1 + Nc
     &   **2)*(-4*k1k2**2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*zim**2 + 3*k1k2
     &   *M2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*zim**2 - 6*k1k2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2
     &   *(p1p2 - p2k1))*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*zim**2 - 2*k1k2*p1k2*sdip*zim**2 + 2*k1k2*((-1 + ab1)*(-k1k2 
     &   + p1k2) + ab2*p2k2)*sdip*zim**2 - p1k2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*sdip*zim**2 + ((-1 + ab1)*(-k1k2 + p1k2) 
     &   + ab2*p2k2)*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*sdip*zim**2 + 2*k1k2**2*(M2 - tdip)*zim**2 - (ab1*(-k1k2 + p1k2) 
     &   + ab2*p2k2)**2*(M2 - 2*p1k1 - tdip)*zim**2 - 3*k1k2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*tdip*zim**2 + 8*k1k2*k1k3
     &   *(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*zim*zjm - 3*k1k3*M2*(ab1*(-k1k2 
     &   + p1k2) + ab2*p2k2)*zim*zjm + 6*k1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))
     &   *(ab1*(-k1k2 + p1k2) + ab2*p2k2)*zim*zjm - 3*k1k2*M2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zim*zjm + 6*k1k2*(ab2**2
     &   *M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zim*zjm 
     &   + 2*M2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zim*zjm - 4*p1k1*(ab1*(-k1k2 + p1k2) + ab2
     &   *p2k2)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zim*zjm - 2*k2k3*p1k1*sdip*zim*zjm + 2*k1k3*p1k2*sdip*zim*zjm + 2*k1k2
     &   *p1k3*sdip*zim*zjm - 2*k2k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*sdip
     &   *zim*zjm - 2*k1k3*((-1 + ab1)*(-k1k2 + p1k2) + ab2*p2k2)*sdip*zim*zjm + p1k3*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*sdip
     &   *zim*zjm - 2*k1k2*((-1 + ab1)*(-k1k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm - (ab1*(-k1k2 + p1k2) + ab2*p2k2)*((-1 
     &   + ab1)*(-k1k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm + p1k2*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm - ((-1 + ab1)
     &   *(-k1k2 + p1k2) + ab2*p2k2)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip*zim*zjm - 4*k1k2*k1k3*(M2 - tdip)*zim*zjm 
     &   + 2*k2k3*sdip*(M2 - tdip)*zim*zjm + 3*k1k3*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*tdip*zim*zjm + 3*k1k2*(ab1*(-k1k3 
     &   + p1k3) + ab2*p2k3)*tdip*zim*zjm - 2*(ab1*(-k1k2 + p1k2) + ab2*p2k2)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*tdip*zim*zjm 
     &   - 4*k1k3**2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*zjm**2 + 3*k1k3*M2
     &   *(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zjm**2 - 6*k1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 
     &   - p2k1))*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*zjm**2 - 2*k1k3*p1k3*sdip*zjm**2 + 2*k1k3*((-1 + ab1)*(-k1k3 + p1k3) 
     &   + ab2*p2k3)*sdip*zjm**2 - p1k3*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip*zjm**2 + ((-1 + ab1)*(-k1k3 + p1k3) 
     &   + ab2*p2k3)*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*sdip*zjm**2 + 2*k1k3**2*(M2 - tdip)*zjm**2 - (ab1*(-k1k3 + p1k3) 
     &   + ab2*p2k3)**2*(M2 - 2*p1k1 - tdip)*zjm**2 - 3*k1k3*(ab1*(-k1k3 + p1k3) + ab2*p2k3)*tdip*zjm**2)
     &   *SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(8*sdip**2))

!************************************************   Mt^ast Vgg Mt   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_b = stsT2xx_M2Dipole_SpinTensor_D23_b +(((-1 + Nc**2)*(4*ab1*
     &   *3*(M2 - 2*p1k1)*p1k1 + 4*ab2**3*M2*p2k1 + (-M2 + p1k1)*sdip + ab2*(2*p1p2 - p2k1)*(2*p1k1 + sdip) - ab2**2*(4
     &   *M2*p1k1 + 4*p1p2*p2k1 - 2*p2k1**2 + M2*sdip) + ab1**2*(2*p1k1*(7*p1k1 + 4*ab2*p1p2 - 8*ab2*p2k1 + sdip) - M2*(8
     &   *p1k1 - 4*ab2*p2k1 + sdip)) + ab1*(-6*p1k1**2 + p1k1*(-12*ab2*p1p2 + 16*ab2*p2k1 - 3*sdip) + 2*ab2*(p1p2 - p2k1)
     &   *(4*ab2*p2k1 - sdip) + 2*M2*(2*(1 + ab2**2)*p1k1 - 2*ab2*p2k1 + sdip)))*(2*p1k2*zim - ab2*p2k2*zim - 2*p1k3
     &   *zjm + ab2*p2k3*zjm + ab1*(k1k2*zim - p1k2*zim - k1k3*zjm + p1k3*zjm))**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(2*Nc*sdip*(M2 - tdip)**2))

!************************************************   2 Re(Mt Vgg Mu^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_b = stsT2xx_M2Dipole_SpinTensor_D23_b + 2*dble(-0.25d0*((-1 + Nc
     &   **2)*(2*ab1**2*(4*M2 - 7*p1k1)*p1k1 + ab2**2*(4*M2*p1k1 + 4*p1p2*p2k1 - 2*p2k1**2) + 2*(M2 - p1k1)*sdip - ab2
     &   *(2*p1p2 - p2k1)*(4*p1k1 + sdip) - 2*ab1*M2*(4*p1k1 - 2*ab2*p2k1 + sdip) + ab1*p1k1*(12*p1k1 + 4*ab2*(3*p1p2 - 4
     &   *p2k1) + 3*sdip))*((-2 + ab1)*k1k2*zim - (-2 + ab1)*p1k2*zim - ab2*p2k2*zim + 2*k1k3*zjm - ab1*k1k3*zjm - 2*p1k3
     &   *zjm + ab1*p1k3*zjm + ab2*p2k3*zjm)*(2*p1k2*zim - ab2*p2k2*zim - 2*p1k3*zjm + ab2*p2k3*zjm + ab1*(k1k2*zim 
     &   - p1k2*zim - k1k3*zjm + p1k3*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(Nc*sdip*(M2 - tdip)*(M2 - udip)))

!************************************************   2 Re(Mt Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_b = stsT2xx_M2Dipole_SpinTensor_D23_b + 2*dble(((-1 + Nc**2)*(ab1
     &   *M2*(-k1k2 + p1k2)*zim + 2*k1k2*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))
     &   *zim + ab2*M2*p2k2*zim + ab1*(k1k2 - p1k2)*sdip*zim + (-k1k2 + p1k2)*sdip*zim - ab2*p2k2*sdip*zim + ab1*(k1k2 
     &   - p1k2)*tdip*zim - ab2*p2k2*tdip*zim + ab1*M2*(k1k3 - p1k3)*zjm - 2*k1k3*(ab2**2*M2 + (-1 + ab1)*ab1*(M2 - 2*p1k1) 
     &   + (-1 + 2*ab1)*ab2*(p1p2 - p2k1))*zjm - ab2*M2*p2k3*zjm + (k1k3 - p1k3)*sdip*zjm + ab1*(-k1k3 + p1k3)*sdip
     &   *zjm + ab2*p2k3*sdip*zjm + ab1*(-k1k3 + p1k3)*tdip*zjm + ab2*p2k3*tdip*zjm)*(-2*p1k2*zim + ab2*p2k2*zim + 2*p1k3
     &   *zjm - ab2*p2k3*zjm + ab1*(-(k1k2*zim) + p1k2*zim + k1k3*zjm - p1k3*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)
     &   **2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(2*Nc*sdip*(M2 - tdip)))

!************************************************   Mu^ast Vgg Mu   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_b = stsT2xx_M2Dipole_SpinTensor_D23_b +(((-1 + Nc**2)**2*(ab1
     &   *(-4*M2*p1k1 + 6*p1k1**2) + 2*ab2*p1k1*(-2*p1p2 + p2k1) + (M2 - p1k1)*sdip)*(-(ab2*p2k2*zim) + ab2*p2k3*zjm + (
     &   -2 + ab1)*(k1k2*zim - p1k2*zim + (-k1k3 + p1k3)*zjm))**2*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**4)/(2*Nc*sdip*(M2 - udip)**2))

!************************************************   2 Re(Mu Vgg Mv^ast)     ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_b = stsT2xx_M2Dipole_SpinTensor_D23_b + 2*dble(((2 - 3*Nc**2 
     &   + Nc**4)*(2*ab2*p1k1*p2k2*zim - p1k2*sdip*zim + k1k2*(M2 - 2*ab1*p1k1 - tdip)*zim - k1k3*M2*zjm - 2*ab2*p1k1*p2k3
     &   *zjm + p1k3*sdip*zjm + k1k3*tdip*zjm + 2*ab1*p1k1*(p1k2*zim + k1k3*zjm - p1k3*zjm))*(ab2*p2k2*zim - ab2*p2k3
     &   *zjm - (-2 + ab1)*(k1k2*zim - p1k2*zim + (-k1k3 + p1k3)*zjm))*SfSfg(isf1,itt1,igen1,isf2,itt2,igen2)**2*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2))/(4.
     &   *Nc*sdip*(M2 - udip)))

!************************************************   Mv^ast Vgg Mv   ************************************************

      stsT2xx_M2Dipole_SpinTensor_D23_b = stsT2xx_M2Dipole_SpinTensor_D23_b +(((4 - 5*Nc**2 + Nc**4)
     &   *(-(k2k3*sdip*zim*zjm) + 2*ab1*(k1k2*zim - k1k3*zjm)*(k1k2*zim - p1k2*zim - k1k3*zjm + p1k3*zjm) - 2*ab2*(k1k2
     &   *zim - k1k3*zjm)*(p2k2*zim - p2k3*zjm))*SfSfgg(isf1,itt1,igen1,isf2,itt2,igen2)**2)/(4*Nc*sdip))

      end

      


