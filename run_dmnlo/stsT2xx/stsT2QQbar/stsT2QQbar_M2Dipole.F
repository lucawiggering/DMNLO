      double precision function stsT2QQbar_M2Dipole()
      implicit none

#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      double precision kappa,Nf
      double precision treeT1T3,treeT1T2,tree,treeT1T4,stsT2QQbar_M2Tree
      double precision sbar,SqrtLambdaAB,d2,Cdip,d1,colorTree,gamQ,gammaQ,Kq
      double precision tbar,etaAt,rhox0t,Ax0t,Bx0t,ubar,etaAu,rhox0u,Ax0u,Bx0u
      double complex cspen

! Color correlated Born matrix element squared
      tree = stsT2QQbar_M2Tree()
      treeT1T3 = (2-Nc**2)/(2*Nc)*tree
      treeT1T2 = tree/(2*Nc)
      treeT1T4 = tree/(-Nc)    

! Useful variables
      gamQ = (3*CF)/2d0
      gammaQ = gamQ*IRdiv - CF/2d0 
      Kq = CF*(7d0/2d0 - Pi**2/6d0)

      tbar = tman - m1**2
      ubar = uman - m1**2
      etaAt = m1**2/(-tbar)
      etaAu = m1**2/(-ubar)

      rhox0t = dsqrt(1+4*etaAt*(x0FEISQ-1)*x0FEISQ)
      Ax0t = rhox0t*(2*x0FEISQ-3-rhox0t)/(2*(1-x0FEISQ)*((1-etaAt)*x0FEISQ-2))
      Bx0t = (rhox0t**2+rhox0t)/(2*(1-x0FEISQ)*etaAt*x0FEISQ)

      rhox0u = dsqrt(1+4*etaAu*(x0FEISQ-1)*x0FEISQ)
      Ax0u = rhox0u*(2*x0FEISQ-3-rhox0u)/(2*(1-x0FEISQ)*((1-etaAu)*x0FEISQ-2))      
      Bx0u = (rhox0u**2+rhox0u)/(2*(1-x0FEISQ)*etaAu*x0FEISQ)
   
! Dipole variables for IE-IS
      sbar = sman - 2*m1**2 
      SqrtLambdaAB = kappa(sman,m1**2,m1**2)
      d2 = sbar/SqrtLambdaAB
      d1 = sbar + 2*m1**2 - SqrtLambdaAB
      Cdip = 2*SqrtLambdaAB/d1

! Number of effectively massless quarks
      Nf = 4d0

! Definition of Mandelstam variables: t=(k1-p1)^2=(k2-p2)^2, u=(k2-p1)^2=(k1-p2)^2

! Don't forget to include c_epsilon in the poles for the integrated dipoles that contain double poles (gg + Qg) -----> extra finite terms
! ------------ double poles -----------
      if(EpsPole.eq.-2d0) then

! ****************************************** FE-FS *************************************************
! D_{13,2} and D_{23,1}
      colorTree = treeT1T2/CF ! T3 T4 = T1 T2
      stsT2QQbar_M2Dipole = - 2*(CF*colorTree*gs**2)/(8*Pi**2)

! ****************************************** FE-IS *************************************************
! D_{31}^a and  D_{32}^b    
      colorTree = treeT1T3/CF
      stsT2QQbar_M2Dipole = stsT2QQbar_M2Dipole - 2*(CF*colorTree*gs**2)/(8*Pi**2)

! D_{32}^a and D_{31}^b
      colorTree = treeT1T4/CF 
      stsT2QQbar_M2Dipole = stsT2QQbar_M2Dipole - 2*(CF*colorTree*gs**2)/(8*Pi**2)

! Rsult for Nf massless quarks
      stsT2QQbar_M2Dipole = Nf*stsT2QQbar_M2Dipole

! ------------ simple poles -----------
      elseif(EpsPole.eq.-1d0) then

! ****************************************** FE-FS *************************************************
! D_{13,2} and D_{23,1}
      colorTree = treeT1T2/CF ! T3 T4 = T1 T2
      stsT2QQbar_M2Dipole = - 2*(colorTree*gs**2*(gamQ + CF*dlog(muSc**2/sman)))/(8*Pi**2)

! ****************************************** FE-IS *************************************************
! D_{31}^a and  D_{32}^b: (Q^2 -> tman)      
      colorTree = treeT1T3/CF
      stsT2QQbar_M2Dipole = stsT2QQbar_M2Dipole - 2*(colorTree*gs**2*(gamQ+ CF*dlog(1 + etaAt) + CF*dlog(-(muSc**2/tbar))))/(8*Pi**2)

! D_{32}^a and D_{31}^b:  (Q^2 -> uman)
      colorTree = treeT1T4/CF 
      stsT2QQbar_M2Dipole = stsT2QQbar_M2Dipole - 2*(colorTree*gs**2*(gamQ+ CF*dlog(1 + etaAu) + CF*dlog(-(muSc**2/ubar))))/(8*Pi**2)

! ****************************************** IE-FS *************************************************
! D^{a3}_1 and D^{b3}_2: (Q^2 -> tman)
      colorTree = treeT1T3/CF 
      stsT2QQbar_M2Dipole = stsT2QQbar_M2Dipole + (CF*colorTree*gs**2*(-1 + dlog((1 + etaAt)/etaAt)))/(4*Pi**2) 

! D^{a3}_2 and D^{b3}_1: (Q^2 -> uman)
      colorTree = treeT1T4/CF 
      stsT2QQbar_M2Dipole = stsT2QQbar_M2Dipole + (CF*colorTree*gs**2*(-1 + dlog((1 + etaAu)/etaAu)))/(4*Pi**2)

! ****************************************** IE-IS *************************************************
! D^{a3,b} and D^{b3,a}
      colorTree = treeT1T2/CF
      stsT2QQbar_M2Dipole = stsT2QQbar_M2Dipole + 2*(CF*colorTree*gs**2*(-1 + d2*dlog(1 + Cdip)))/(8*Pi**2)

! Rsult for Nf massless quarks
      stsT2QQbar_M2Dipole = Nf*stsT2QQbar_M2Dipole

! ------------ finite pieces ---------
      elseif(EpsPole.eq.0d0) then

! ****************************************** FE-FS *************************************************
! D_{13,2} and D_{23,1}
      colorTree = treeT1T2/CF ! T3 T4 = T1 T2
      stsT2QQbar_M2Dipole = - 2*(CF*colorTree*gs**2*IRdiv**2)/(8*Pi**2) - 2*(CF*colorTree*gs**2*IRdiv*dlog(muSc**2/sman))/(8*Pi**2) + 
     &  2*(colorTree*gs**2*(-2*GammaQ - 2*gamQ - 2*Kq + CF*Pi**2 - 2*gamQ*dlog(muSc**2/sman) - CF*dlog(muSc**2/sman)**2))/(16*Pi**2) ! CS short eq. 6.16
          
! ****************************************** FE-IS *************************************************
! D_{31}^a and  D_{32}^b: (Q^2 -> tman)      
      colorTree = treeT1T3/CF
      stsT2QQbar_M2Dipole = stsT2QQbar_M2Dipole - 2*(CF*colorTree*gs**2*IRdiv**2)/(8*Pi**2) 
     & - 2*(colorTree*gs**2*IRdiv*(CF*dlog(1 + etaAt) + CF*dlog(-(muSc**2/tbar))))/(8*Pi**2) + 
     &  2*(colorTree*gs**2*(-6*GammaQ - 6*gamQ*dlog(-(muSc**2/tbar)) - 6*Kq + CF*Pi**2 - 12*CF*dble(cspen(dcmplx(-Ax0t))) - 3*CF*dlog(1 + etaAt)**2 - 6*CF*dlog(1 + etaAt)*dlog(-(muSc**2/tbar)) - 
     &       3*CF*dlog(-(muSc**2/tbar))**2 + 9*CF*dlog(1 - x0FEISQ)))/(48*Pi**2)

! D_{32}^a and D_{31}^b:  (Q^2 -> uman)
      colorTree = treeT1T4/CF 
      stsT2QQbar_M2Dipole = stsT2QQbar_M2Dipole  - 2*(CF*colorTree*gs**2*IRdiv**2)/(8*Pi**2) 
     &   - 2*(colorTree*gs**2*IRdiv*(CF*dlog(1 + etaAu) + CF*dlog(-(muSc**2/ubar))))/(8*Pi**2) + 
     &  2*(colorTree*gs**2*(-6*GammaQ - 6*gamQ*dlog(-(muSc**2/ubar)) - 6*Kq + CF*Pi**2 - 12*CF*dble(cspen(dcmplx(-Ax0u))) - 3*CF*dlog(1 + etaAu)**2 - 6*CF*dlog(1 + etaAu)*dlog(-(muSc**2/ubar)) - 
     &       3*CF*dlog(-(muSc**2/ubar))**2 + 9*CF*dlog(1 - x0FEISQ)))/(48*Pi**2)

! ****************************************** IE-FS *************************************************
! D^{a3}_1 and D^{b3}_2: (Q^2 -> tman)
      colorTree = treeT1T3/CF 
      stsT2QQbar_M2Dipole = stsT2QQbar_M2Dipole + (CF*colorTree*gs**2*IRdiv*(-1 + dlog((1 + etaAt)/etaAt)))/(4*Pi**2) + 
     &  (CF*colorTree*gs**2*(2*dble(cspen(dcmplx(-Ax0t))) - 2*dble(cspen(dcmplx(-Bx0t))) - 2*dlog(1 + Bx0t) - dlog(etaAt) - 
     &       dlog(etaAt)**2/2d0 + dlog(1 + etaAt)**2/2d0 - (1 - dlog((1 + etaAt)/etaAt))*dlog(-(muSc**2/tbar))))/(4*Pi**2)
 
! D^{a3}_2 and D^{b3}_1: (Q^2 -> uman)
      colorTree = treeT1T4/CF 
      stsT2QQbar_M2Dipole = stsT2QQbar_M2Dipole + (CF*colorTree*gs**2*IRdiv*(-1 + dlog((1 + etaAu)/etaAu)))/(4*Pi**2) + 
     &  (CF*colorTree*gs**2*(2*dble(cspen(dcmplx(-Ax0u))) - 2*dble(cspen(dcmplx(-Bx0u))) - 2*dlog(1 + Bx0u) - dlog(etaAu) - 
     &       dlog(etaAu)**2/2d0 + dlog(1 + etaAu)**2/2d0 - (1 - dlog((1 + etaAu)/etaAu))*dlog(-(muSc**2/ubar))))/(4*Pi**2)
 
! ****************************************** IE-IS *************************************************
! D^{a3,b} and D^{b3,a}       
      colorTree = treeT1T2/CF
      stsT2QQbar_M2Dipole = stsT2QQbar_M2Dipole + 2*(CF*colorTree*gs**2*IRdiv*(-1 + d2*dlog(1 + Cdip)))/(8*Pi**2) + 
     &  2*(CF*colorTree*gs**2*(-(((2 + Cdip)*dlog(1 + Cdip))/Cdip) - (d2*(4*dble(cspen(dcmplx(-Cdip))) + dlog(1 + Cdip)**2))/2d0 - 
     &       (1 - d2*dlog(1 + Cdip))*dlog(muSc**2/sbar) - (1- d2*dlog(1 + Cdip))*dlog(sman/(sbar*(1 - x0IEISQ)**2))))/(8*Pi**2)

! Rsult for Nf massless quarks
      stsT2QQbar_M2Dipole = Nf*stsT2QQbar_M2Dipole


      endif

      end
