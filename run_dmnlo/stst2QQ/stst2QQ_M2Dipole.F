      double precision function stst2QQ_M2Dipole()
      implicit none

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"

      double precision kappa,vijktilde,rho,sqrtLambda,Adip,Qbar2,Bdip,vDip,rhoJ,rhoK,etaJ
      double precision treeT1T3,treeT1T2,tree,treeT1T4,stst2qq_M2TreeDipoleT1T3
      double precision sbar,SqrtLambdaAB,d2,Cdip,d1,colorTree,stst2qq_M2TreeDipoleT1T4
      double complex cspen

! stst2QQ_M2Tree() has to be called before stst2QQ_M2Dipole() to ensure that stst2QQ_M2Tree3 
! and stst2QQ_M2Tree6 are set correctly

      treeT1T3 = stst2qq_M2TreeDipoleT1T3(sqrtS**2,tman,uman)
      treeT1T2 = -1d0/(2*Nc)*((Nc+1)*stst2QQ_M2Tree3 - (Nc-1)*stst2QQ_M2Tree6)
      tree = stst2QQ_M2Tree3 + stst2QQ_M2Tree6
      treeT1T4 = -CF*tree - treeT1T2 - treeT1T3


! Dipole variables for FE-FS
      vijktilde = kappa(1d0,muu3**2,muu4**2)/(1-muu3**2-muu4**2)
      rho = dsqrt((1-vijktilde)/(1+vijktilde))
   
! Dipole variables for IE-IS
      sbar = sqrtS**2 - m1**2 - m2**2
      SqrtLambdaAB = kappa(sqrtS**2,m1**2,m2**2)
      d2 = sbar/SqrtLambdaAB

! Definition of Mandelstam variables: t=(k1-p1)^2=(k2-p2)^2, u=(k2-p1)^2=(k1-p2)^2
! with p_1^2 = m_1^2,  p_2^2 = m_2^2,  k_1^2 = m_3^2,  k_2^2 = m_4^2 
! ------------ double poles -----------
      if(EpsPole.eq.-2d0) then

      stst2QQ_M2Dipole = 0d0
        
! ------------ simple poles -----------
      elseif(EpsPole.eq.-1d0) then
! ****************************************** FE-FS (from hep-ph/0201036 section 5.1) ***************
! I_{13,2}   ----> k1 as emitter and k2 as spectator
      colorTree = treeT1T2/CF ! =T3T4/CF
      stst2QQ_M2Dipole = -(CF*colorTree*gs**2*(vijktilde + dlog(rho)))/(8*Pi**2*vijktilde)

! I_{23,1}   ----> k2 as emitter and k1 as spectator
      colorTree = treeT1T2/CF ! =T3T4/CF
      stst2QQ_M2Dipole = stst2QQ_M2Dipole -(CF*colorTree*gs**2*(vijktilde + dlog(rho)))/(8*Pi**2*vijktilde)

! ****************************************** FE-IS *************************************************
! I_{13}^a   ----> k1 as emitter and p1 as spectator (Q-> t)
      sqrtLambda = kappa(tman,m3**2,m1**2)
      Qbar2 = tman - m1**2 - m3**2
      Adip = 2*sqrtLambda/(2*m3**2-Qbar2-sqrtLambda)
      vdip = sqrtLambda/(-Qbar2)
      
      colorTree = treeT1T3/CF
      stst2QQ_M2Dipole = stst2QQ_M2Dipole + (CF*colorTree*gs**2*(-1 + dlog(1 + Adip)/vdip))/(8*Pi**2)

! I_{13}^b   ----> k1 as emitter and p2 as spectator (Q-> u)
      sqrtLambda = kappa(uman,m3**2,m2**2)
      Qbar2 = uman - m2**2 - m3**2
      Adip = 2*sqrtLambda/(2*m3**2-Qbar2-sqrtLambda)
      vdip = sqrtLambda/(-Qbar2)

      colorTree = treeT1T4/CF
      stst2QQ_M2Dipole = stst2QQ_M2Dipole + (CF*colorTree*gs**2*(-1 + dlog(1 + Adip)/vdip))/(8*Pi**2)

! I_{23}^a   ----> k2 as emitter and p1 as spectator (Q-> u)
      sqrtLambda = kappa(uman,m4**2,m1**2)
      Qbar2 = uman - m1**2 - m4**2
      Adip = 2*sqrtLambda/(2*m4**2-Qbar2-sqrtLambda)
      vdip = sqrtLambda/(-Qbar2)

      colorTree = treeT1T4/CF ! =T2T3
      stst2QQ_M2Dipole = stst2QQ_M2Dipole + (CF*colorTree*gs**2*(-1 + dlog(1 + Adip)/vdip))/(8*Pi**2)

! I_{23}^b   ----> k2 as emitter and p2 as spectator (Q-> t)
      sqrtLambda = kappa(tman,m4**2,m2**2)
      Qbar2 = tman - m2**2 - m4**2
      Adip = 2*sqrtLambda/(2*m4**2-Qbar2-sqrtLambda)
      vdip = sqrtLambda/(-Qbar2)

      colorTree = treeT1T3/CF ! T2T4=T1T3
      stst2QQ_M2Dipole = stst2QQ_M2Dipole + (CF*colorTree*gs**2*(-1 + dlog(1 + Adip)/vdip))/(8*Pi**2)

! ****************************************** IE-FS *************************************************
! I^{a3}_1   ----> p1 as emitter and k1 as spectator (Q-> t)
      sqrtLambda = kappa(tman,m3**2,m1**2)
      Qbar2 = tman - m1**2 - m3**2
      Adip = 2*sqrtLambda/(2*m3**2-Qbar2-sqrtLambda)
      vdip = sqrtLambda/(-Qbar2)
      Bdip = 2*vdip/(1d0 - vdip)

      colorTree = treeT1T3/CF 
      stst2QQ_M2Dipole = stst2QQ_M2Dipole - (CF*colorTree*gs**2*(vDip - dlog((1 + Bdip)/(1 + Adip))))/(8*Pi**2*vDip)

! I^{b3}_1   ----> p2 as emitter and k1 as spectator (Q-> u)
      sqrtLambda = kappa(uman,m3**2,m2**2)
      Qbar2 = uman - m2**2 - m3**2
      Adip = 2*sqrtLambda/(2*m3**2-Qbar2-sqrtLambda)
      vdip = sqrtLambda/(-Qbar2)
      Bdip = 2*vdip/(1d0- vdip)

      colorTree = treeT1T4/CF 
      stst2QQ_M2Dipole = stst2QQ_M2Dipole - (CF*colorTree*gs**2*(vDip - dlog((1 + Bdip)/(1 + Adip))))/(8*Pi**2*vDip)

! I^{a3}_2   ----> p1 as emitter and k2 as spectator (Q-> u)
      sqrtLambda = kappa(uman,m4**2,m1**2)
      Qbar2 = uman - m1**2 - m4**2
      Adip = 2*sqrtLambda/(2*m4**2-Qbar2-sqrtLambda)
      vdip = sqrtLambda/(-Qbar2)
      Bdip = 2*vdip/(1d0 - vdip)

      colorTree = treeT1T4/CF 
      stst2QQ_M2Dipole = stst2QQ_M2Dipole - (CF*colorTree*gs**2*(vDip - dlog((1 + Bdip)/(1 + Adip))))/(8*Pi**2*vDip)

! I^{b3}_2   ----> p2 as emitter and k2 as spectator (Q-> t)
      sqrtLambda = kappa(tman,m4**2,m2**2)
      Qbar2 = tman - m2**2 - m4**2
      Adip = 2*sqrtLambda/(2*m4**2-Qbar2-sqrtLambda)
      vdip = sqrtLambda/(-Qbar2)
      Bdip = 2*vdip/(1d0 - vdip)

      colorTree = treeT1T3/CF !T2T4=T1T3
      stst2QQ_M2Dipole = stst2QQ_M2Dipole - (CF*colorTree*gs**2*(vDip - dlog((1 + Bdip)/(1 + Adip))))/(8*Pi**2*vDip)

! ****************************************** IE-IS  ***************
! I^{a3,b}   ----> p1 as emitter and p2 as spectator
      d1 = sbar + 2*m1**2 - SqrtLambdaAB
      Cdip = 2*SqrtLambdaAB/d1
      colorTree = treeT1T2/CF
      stst2QQ_M2Dipole = stst2QQ_M2Dipole - (CF*colorTree*gs**2*(1 - d2*dlog(1 + Cdip)))/(8*Pi**2)

! I^{b3,a}   ----> p2 as emitter and p1 as spectator
      d1 = sbar + 2*m2**2 - SqrtLambdaAB
      Cdip = 2*SqrtLambdaAB/d1
      colorTree = treeT1T2/CF
      stst2QQ_M2Dipole = stst2QQ_M2Dipole - (CF*colorTree*gs**2*(1 - d2*dlog(1 + Cdip)))/(8*Pi**2)

! ------------ finite pieces ---------
      elseif(EpsPole.eq.0d0) then

      ! ****************************************** FE-FS (from hep-ph/0201036 section 5.1) ***************
! I_{13,2}   ----> k1 as emitter and k2 as spectator
      rhoJ = dsqrt((1 + (2*muu3**2)/(1 - muu3**2 - muu4**2) - vijktilde)/(1 + (2*muu3**2)/(1 - muu3**2 - muu4**2) + vijktilde))
      rhoK = dsqrt((1 + (2*muu4**2)/(1 - muu3**2 - muu4**2) - vijktilde)/(1 + (2*muu4**2)/(1 - muu3**2 - muu4**2) + vijktilde))

      colorTree = treeT1T2/CF ! =T3T4/CF
      stst2QQ_M2Dipole = -IRdiv*(CF*colorTree*gs**2*(vijktilde + dlog(rho)))/(8*Pi**2*vijktilde)
     &    + (CF*colorTree*gs**2*(-3 + muu4/(1 - muu4) + (2*(1 - 2*muu4)*muu4)/(1 - muu3**2 - muu4**2) - dlog(muu3) + 
     &      2*dlog(-muu3**2 + (1 - muu4)**2) + (2*muu3**2*dlog(muu3/(1 -muu4)))/(1 - muu3**2 - muu4**2) - dlog(1 - muu4) - 
     &      dlog(muSc**2/sqrtS**2)*(1 + dlog(rho)/vijktilde) - 
     &      (Pi**2 - 12*dble(cspen(dcmplx(1 - rho))) + 12*dble(cspen(dcmplx(-rho))) - 3*dble(cspen(dcmplx(1 - rhoJ**2))) - 
     &         3*dble(cspen(dcmplx(1 - rhoK**2))) - 6*dlog(1 - (muu3 + muu4)**2)*dlog(rho) - 3*dlog(rhoJ)**2 - 3*dlog(rhoK)**2)/
     &       (3*vijktilde)))/(8*Pi**2)


! I_{23,1}   ----> k2 as emitter and k1 as spectator
      rhoJ = dsqrt((1 + (2*muu4**2)/(1 - muu3**2 - muu4**2) - vijktilde)/(1 + (2*muu4**2)/(1 - muu3**2 - muu4**2) + vijktilde))
      rhoK = dsqrt((1 + (2*muu3**2)/(1 - muu3**2 - muu4**2) - vijktilde)/(1 + (2*muu3**2)/(1 - muu3**2 - muu4**2) + vijktilde))

      colorTree = treeT1T2/CF ! =T3T4/CF
      stst2QQ_M2Dipole = stst2QQ_M2Dipole -IRdiv*(CF*colorTree*gs**2*(vijktilde + dlog(rho)))/(8*Pi**2*vijktilde)
     &   + (CF*colorTree*gs**2*(-3 + muu3/(1 - muu3) + (2*(1 - 2*muu3)*muu3)/(1 - muu3**2 - muu4**2) - dlog(1 - muu3) - dlog(muu4) + 
     &      (2*muu4**2*dlog(muu4/(1 - muu3)))/(1 - muu3**2 - muu4**2) + 2*dlog((1 - muu3)**2 - muu4**2) - 
     &      dlog(muSc**2/sqrtS**2)*(1 + dlog(rho)/vijktilde) - 
     &      (Pi**2 - 12*dble(cspen(dcmplx(1 - rho))) + 12*dble(cspen(dcmplx(-rho))) - 3*dble(cspen(dcmplx(1 - rhoJ**2))) - 
     &         3*dble(cspen(dcmplx(1 - rhoK**2))) - 6*dlog(1 - (muu3 + muu4)**2)*dlog(rho) - 3*dlog(rhoJ)**2 - 3*dlog(rhoK)**2)/
     &       (3*vijktilde)))/(8*Pi**2)
     
! ****************************************** FE-IS *************************************************
! I_{13}^a   ----> k1 as emitter and p1 as spectator (Q-> t)
      sqrtLambda = kappa(tman,m3**2,m1**2)
      Qbar2 = tman - m1**2 - m3**2
      Adip = 2*sqrtLambda/(2*m3**2-Qbar2-sqrtLambda)
      vdip = sqrtLambda/(-Qbar2)
      etaJ = m3**2/(-Qbar2)
      
      colorTree = treeT1T3/CF
      stst2QQ_M2Dipole = stst2QQ_M2Dipole - IRdiv*(CF*colorTree*gs**2*(vDip - dlog(1 + Adip)))/(8*Pi**2*vDip) + 
     &  - (CF*colorTree*gs**2*(4*dble(cspen(dcmplx(-Adip))) + dlog(1 + Adip)**2- 2*dlog(1 + Adip)*(dlog(-(muSc**2/Qbar2)) + dlog(etaJ/(-1 + x0FEIS)**2)) + 
     &       2*vDip*(2 + dlog(-(muSc**2/Qbar2)) + dlog(etaJ/(-1 + x0FEIS)**2))))/(16*Pi**2*vDip)
 
! I_{13}^b   ----> k1 as emitter and p2 as spectator (Q-> u)
      sqrtLambda = kappa(uman,m3**2,m2**2)
      Qbar2 = uman - m2**2 - m3**2
      Adip = 2*sqrtLambda/(2*m3**2-Qbar2-sqrtLambda)
      vdip = sqrtLambda/(-Qbar2)
      etaJ = m3**2/(-Qbar2)

      colorTree = treeT1T4/CF
      stst2QQ_M2Dipole = stst2QQ_M2Dipole - IRdiv*(CF*colorTree*gs**2*(vDip - dlog(1 + Adip)))/(8*Pi**2*vDip) + 
     &  - (CF*colorTree*gs**2*(4*dble(cspen(dcmplx(-Adip))) + dlog(1 + Adip)**2- 2*dlog(1 + Adip)*(dlog(-(muSc**2/Qbar2)) + dlog(etaJ/(-1 + x0FEIS)**2)) + 
     &       2*vDip*(2 + dlog(-(muSc**2/Qbar2)) + dlog(etaJ/(-1 + x0FEIS)**2))))/(16*Pi**2*vDip)
 

! I_{23}^a   ----> k2 as emitter and p1 as spectator (Q-> u)
      sqrtLambda = kappa(uman,m4**2,m1**2)
      Qbar2 = uman - m1**2 - m4**2
      Adip = 2*sqrtLambda/(2*m4**2-Qbar2-sqrtLambda)
      vdip = sqrtLambda/(-Qbar2)
      etaJ = m4**2/(-Qbar2)

      colorTree = treeT1T4/CF ! =T2T3
      stst2QQ_M2Dipole = stst2QQ_M2Dipole - IRdiv*(CF*colorTree*gs**2*(vDip - dlog(1 + Adip)))/(8*Pi**2*vDip) + 
     &  - (CF*colorTree*gs**2*(4*dble(cspen(dcmplx(-Adip))) + dlog(1 + Adip)**2- 2*dlog(1 + Adip)*(dlog(-(muSc**2/Qbar2)) + dlog(etaJ/(-1 + x0FEIS)**2)) + 
     &       2*vDip*(2 + dlog(-(muSc**2/Qbar2)) + dlog(etaJ/(-1 + x0FEIS)**2))))/(16*Pi**2*vDip)
 

! I_{23}^b   ----> k2 as emitter and p2 as spectator (Q-> t)
      sqrtLambda = kappa(tman,m4**2,m2**2)
      Qbar2 = tman - m2**2 - m4**2
      Adip = 2*sqrtLambda/(2*m4**2-Qbar2-sqrtLambda)
      vdip = sqrtLambda/(-Qbar2)
      etaJ = m4**2/(-Qbar2)

      colorTree = treeT1T3/CF ! T2T4=T1T3
      stst2QQ_M2Dipole = stst2QQ_M2Dipole - IRdiv*(CF*colorTree*gs**2*(vDip - dlog(1 + Adip)))/(8*Pi**2*vDip) + 
     &  - (CF*colorTree*gs**2*(4*dble(cspen(dcmplx(-Adip))) + dlog(1 + Adip)**2- 2*dlog(1 + Adip)*(dlog(-(muSc**2/Qbar2)) + dlog(etaJ/(-1 + x0FEIS)**2)) + 
     &       2*vDip*(2 + dlog(-(muSc**2/Qbar2)) + dlog(etaJ/(-1 + x0FEIS)**2))))/(16*Pi**2*vDip)

! ****************************************** IE-FS *************************************************
! I^{a3}_1   ----> p1 as emitter and k1 as spectator (Q-> t)
      sqrtLambda = kappa(tman,m3**2,m1**2)
      Qbar2 = tman - m1**2 - m3**2
      Adip = 2*sqrtLambda/(2*m3**2-Qbar2-sqrtLambda)
      vdip = sqrtLambda/(-Qbar2)
      Bdip = 2*vdip/(1d0 - vdip)
      etaJ = m3**2/(-Qbar2)

      colorTree = treeT1T3/CF 
      stst2QQ_M2Dipole = stst2QQ_M2Dipole  - IRdiv*(CF*colorTree*gs**2*(vDip - dlog((1 + Bdip)/(1 + Adip))))/(8*Pi**2*vDip) + 
     &  - (CF*colorTree*gs**2*((1 - dlog((1 + Bdip)/(1 + Adip))/vDip)*dlog(-(muSc**2/Qbar2)) + 
     &       (-2*dble(cspen(dcmplx(-Adip))) + 2*dble(cspen(dcmplx(-Bdip)))+ dlog(1 + Bdip) + (-dlog(1 + Adip)**2 + dlog(1 + Bdip)**2)/2d0 + 
     &          (vDip + dlog((1 + Adip)/(1 + Bdip)))*dlog(etaJ/(1 - x0FEIS)**2))/vDip))/(8*Pi**2)
    
! I^{b3}_1   ----> p2 as emitter and k1 as spectator (Q-> u)
      sqrtLambda = kappa(uman,m3**2,m2**2)
      Qbar2 = uman - m2**2 - m3**2
      Adip = 2*sqrtLambda/(2*m3**2-Qbar2-sqrtLambda)
      vdip = sqrtLambda/(-Qbar2)
      Bdip = 2*vdip/(1d0- vdip)
      etaJ = m3**2/(-Qbar2)

      colorTree = treeT1T4/CF 
      stst2QQ_M2Dipole = stst2QQ_M2Dipole  - IRdiv*(CF*colorTree*gs**2*(vDip - dlog((1 + Bdip)/(1 + Adip))))/(8*Pi**2*vDip) + 
     &  - (CF*colorTree*gs**2*((1 - dlog((1 + Bdip)/(1 + Adip))/vDip)*dlog(-(muSc**2/Qbar2)) + 
     &       (-2*dble(cspen(dcmplx(-Adip))) + 2*dble(cspen(dcmplx(-Bdip)))+ dlog(1 + Bdip) + (-dlog(1 + Adip)**2 + dlog(1 + Bdip)**2)/2d0 + 
     &          (vDip + dlog((1 + Adip)/(1 + Bdip)))*dlog(etaJ/(1 - x0FEIS)**2))/vDip))/(8*Pi**2)

! I^{a3}_2   ----> p1 as emitter and k2 as spectator (Q-> u)
      sqrtLambda = kappa(uman,m4**2,m1**2)
      Qbar2 = uman - m1**2 - m4**2
      Adip = 2*sqrtLambda/(2*m4**2-Qbar2-sqrtLambda)
      vdip = sqrtLambda/(-Qbar2)
      Bdip = 2*vdip/(1d0 - vdip)
      etaJ = m4**2/(-Qbar2)

      colorTree = treeT1T4/CF 
      stst2QQ_M2Dipole = stst2QQ_M2Dipole  - IRdiv*(CF*colorTree*gs**2*(vDip - dlog((1 + Bdip)/(1 + Adip))))/(8*Pi**2*vDip) + 
     &  - (CF*colorTree*gs**2*((1 - dlog((1 + Bdip)/(1 + Adip))/vDip)*dlog(-(muSc**2/Qbar2)) + 
     &       (-2*dble(cspen(dcmplx(-Adip))) + 2*dble(cspen(dcmplx(-Bdip)))+ dlog(1 + Bdip) + (-dlog(1 + Adip)**2 + dlog(1 + Bdip)**2)/2d0 + 
     &          (vDip + dlog((1 + Adip)/(1 + Bdip)))*dlog(etaJ/(1 - x0FEIS)**2))/vDip))/(8*Pi**2)

! I^{b3}_2   ----> p2 as emitter and k2 as spectator (Q-> t)
      sqrtLambda = kappa(tman,m4**2,m2**2)
      Qbar2 = tman - m2**2 - m4**2
      Adip = 2*sqrtLambda/(2*m4**2-Qbar2-sqrtLambda)
      vdip = sqrtLambda/(-Qbar2)
      Bdip = 2*vdip/(1d0 - vdip)
      etaJ = m4**2/(-Qbar2)

      colorTree = treeT1T3/CF !T2T4=T1T3
      stst2QQ_M2Dipole = stst2QQ_M2Dipole  - IRdiv*(CF*colorTree*gs**2*(vDip - dlog((1 + Bdip)/(1 + Adip))))/(8*Pi**2*vDip) + 
     &  - (CF*colorTree*gs**2*((1 - dlog((1 + Bdip)/(1 + Adip))/vDip)*dlog(-(muSc**2/Qbar2)) + 
     &       (-2*dble(cspen(dcmplx(-Adip))) + 2*dble(cspen(dcmplx(-Bdip)))+ dlog(1 + Bdip) + (-dlog(1 + Adip)**2 + dlog(1 + Bdip)**2)/2d0 + 
     &          (vDip + dlog((1 + Adip)/(1 + Bdip)))*dlog(etaJ/(1 - x0FEIS)**2))/vDip))/(8*Pi**2)

! ****************************************** IE-IS  ***************
! I^{a3,b}   ----> p1 as emitter and p2 as spectator
      d1 = sbar + 2*m1**2 - SqrtLambdaAB
      Cdip = 2*SqrtLambdaAB/d1
      
      colorTree = treeT1T2/CF
      stst2QQ_M2Dipole = stst2QQ_M2Dipole + (CF*colorTree*gs**2*IRdiv*(-1 + d2*dlog(1 + Cdip)))/(8*Pi**2) + 
     &  (CF*colorTree*gs**2*(-(((2 + Cdip)*dlog(1 + Cdip))/Cdip) - (d2*(4*dble(cspen(dcmplx(-Cdip))) + dlog(1 + Cdip)**2))/2d0 - 
     &       (1 - d2*dlog(1 + Cdip))*dlog(muSc**2/sbar) - (1- d2*dlog(1 + Cdip))*dlog(sqrtS**2/(sbar*(1 - x0IEIS)**2))))/(8*Pi**2)
      
! I^{b3,a}   ----> p2 as emitter and p1 as spectator
      d1 = sbar + 2*m2**2 - SqrtLambdaAB
      Cdip = 2*SqrtLambdaAB/d1

      colorTree = treeT1T2/CF
      stst2QQ_M2Dipole = stst2QQ_M2Dipole + (CF*colorTree*gs**2*IRdiv*(-1 + d2*dlog(1 + Cdip)))/(8*Pi**2) + 
     &  (CF*colorTree*gs**2*(-(((2 + Cdip)*dlog(1 + Cdip))/Cdip) - (d2*(4*dble(cspen(dcmplx(-Cdip))) + dlog(1 + Cdip)**2))/2d0 - 
     &       (1 - d2*dlog(1 + Cdip))*dlog(muSc**2/sbar) - (1- d2*dlog(1 + Cdip))*dlog(sqrtS**2/(sbar*(1 - x0IEIS)**2))))/(8*Pi**2)


      endif
       

      end


! Heaviside theta function
      double precision function stst2QQ_unit_step(x)
      implicit none

      double precision x

      if (x .ge. 0d0) then
       stst2QQ_unit_step=1d0
      else
       stst2QQ_unit_step=0d0
      endif

      end


      double precision function stst2QQ_M2DipolePlus_Quark(colorTreeNum,Q2IsTman)
      implicit none
#include "stst2QQ_DipoleKinematics.h"
#include "stst2QQ_Model.h"

      double precision tdip1,udip1,tdip,udip,colorDiptree1,colorDiptreeX
      double precision stst2qq_M2TreeDipoleT1T4,stst2qq_M2TreeDipoleT1T3
      integer Q2IsTman,colorTreeNum
      
! Q2IsTman should take the value 1 if Q^2 corresponds to the Mandelstam variable t and 
! take 0 if Q^2 corresponds to the Mandelstam variable u
      if (Q2IsTman.eq.1) then
      tdip1 = q2
      udip1 = m2Sum - sdip1 - q2

      tdip = q2X
      udip = m2Sum - sdip - q2X
      else 
      udip1 = q2
      tdip1 = m2Sum - sdip1  - q2

      udip = q2X
      tdip = m2Sum - sdip - q2X
      endif
 
! Decide on colour correlated tree level matrix element which is determined through the variable "colorTreeNum"
      if (colorTreeNum.eq.1) then
      colorDiptree1 = stst2qq_M2TreeDipoleT1T3(sdip1,tdip1,udip1)/CF
      colorDiptreeX = stst2qq_M2TreeDipoleT1T3(sdip,tdip,udip)/CF

      else if (colorTreeNum.eq.2) then
      colorDiptree1 = stst2qq_M2TreeDipoleT1T4(sdip1,tdip1,udip1)/CF
      colorDiptreeX = stst2qq_M2TreeDipoleT1T4(sdip,tdip,udip)/CF
      else
      write(*,*) "colorTreeNum can only take the values 1 and 2!"
      endif
       
! FE-IS
      stst2QQ_M2DipolePlus_Quark = CF*dQ2*2/(xdip**2*(1-xdip))*((xdip-1)**2/(4*(xdip*(etaJX-1)+1)**2) -1 + dlog(1+AdipX)/(Rdip*vdipX))*colorDiptreeX
     &          - 2*CF/(1-xdip)*(-1 + dlog(1+Adip1)/vdip1)*colorDiptree1

! IE-FS
      stst2QQ_M2DipolePlus_Quark = stst2QQ_M2DipolePlus_Quark + CF*dQ2*colorDiptreeX/(1-xdip)*2/(vDipX*Rdip*xdip**2)*(xdip*dlog(1+Bdipx) 
     &       - dlog(1+Adipx) - 4*etaAX*xdip**2*((etaJX-1)*xdip+1)/(1-vdipX**2*Rdip**2)*vdipX*Rdip)
     &          -  CF*colorDiptree1/(1-xdip)*2/vDip1*(dlog(1+Bdip1) - dlog(1+Adip1) - 4*etaA1*etaJ1/(1-vDip1**2)*vDip1)

      stst2QQ_M2DipolePlus_Quark = -gs**2/(8*Pi**2)*stst2QQ_M2DipolePlus_Quark

      end

      double precision function stst2QQ_M2DipolePlus_IEIS()
      implicit none
#include "stst2QQ_DipoleKinematics.h"
#include "stst2QQ_Model.h"

      double precision tdip,udip,tdip1,udip1,colorDiptree1,colorDiptreeX
      double precision stst2qq_M2TreeDipoleT1T2

! Mandelstam variables
      tdip = q2X
      udip = m2Sum - sdip - q2X

      tdip1 = q2
      udip1 = m2Sum - sdip1 - q2

! colour correlated leading order matrix elements squared
      colorDiptreeX = stst2qq_M2TreeDipoleT1T2(sdip,tdip,udip)/CF
      colorDiptree1 = stst2qq_M2TreeDipoleT1T2(sdip1,tdip1,udip1)/CF

      stst2QQ_M2DipolePlus_IEIS = -gs**2*CF/(4*Pi**2*(1-xdip))*(d2*dlog(1+Cdip)-1)*(colorDiptreeX*dQ2/sqrtlambdaSdip - colorDiptree1/SqrtLambdaAB) 

      end

      subroutine stst2QQ_SetDipoleKinematicsFEIS(xi,x0,sman,ma2,mmb2,mk2,mj2,intfac)
      implicit none 
! For the evaluation of the "plus"-distribution one has to distinguish between three different dependencies $f(x,Q^2)$. 
! For this purpose, we introduce two different variables for $Q^2$. These are q2 which is generated from the integration limits
! of $Q^2$ evaluated at x=1 $Q^2_\pm(1)$ and q2X which is generated from the $x$-dependent integration limits $Q^2_\pm(x)$ instead. 
! Variables with the suffix "X" refer to f(x,q2X), variables with the suffix "1" refer to f(1,q2) and those with the suffix "plus" correspond to f(x,q2).
#include "stst2QQ_Model.h"
#include "stst2QQ_DipoleKinematics.h"
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

      subroutine stst2QQ_SetDipoleKinematicsIEIS(xi, x0,sman,ma2,mmb2,m12,m22,intfac)
      implicit none 
#include "stst2QQ_Model.h"
#include "stst2QQ_DipoleKinematics.h"
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