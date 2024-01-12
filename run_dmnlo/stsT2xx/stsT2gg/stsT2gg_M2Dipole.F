! The first long paper by Catani and Seymour (hep-ph/9605323) is called "CS".
      
      double precision function stst2gg_M2Dipole()
      implicit none

#include "../common/stsT2xx_Kinematics.h"
#include "../common/stsT2xx_Model.h"
#include "../../util/DMNLO_LoopIntegrals.h"

      double precision M2_S,M2_8S,M2_8A,tree,treeT1T2,treeT1T3,treeT3T4,treeT1T4,treeT2T3,colorTree
      double precision d1,d2,sbar,SqrtlambdaAB,kappa,Cdip,Nf,M2,GammaG,Kg,gamG
      double complex cspen
      double precision stsT2gg_M2TreeDipole,stsT2gg_M2TreeDipoleT1T2,stsT2gg_M2TreeDipoleT1T3
      double precision tbar,etaAt,rhox0t,Ax0t,Bx0t,ubar,etaAu,rhox0u,Ax0u,Bx0u

!  number of massless quarks
      Nf = 4d0

! color correlated tree matrix elements 
      tree = stsT2gg_M2TreeDipole(sman,tman,uman)
      treeT1T3 = stsT2gg_M2TreeDipoleT1T3(sman,tman,uman) ! T2 T4 = T1 T3
      treeT1T2 = stsT2gg_M2TreeDipoleT1T2(sman,tman,uman)

      treeT3T4 = ((CF-CV)*tree + treeT1T2)
      treeT1T4 = -(CF*tree + treeT1T2 + treeT1T3) ! T2 T3 = T1 T4

! Squared squark mass 
      M2=MSf2(iSf1,itt1,igen1)

! Variables for FE-IS and IE-FS
! The following two variables correspond to \bar{Q}^2
      tbar = tman - M2
      ubar = uman - M2

      etaAt = M2/(-tbar)
      etaAu = M2/(-ubar)

      rhox0t = dsqrt(1+4*etaAt*(x0FEISg-1)*x0FEISg)
      Ax0t = rhox0t*(2*x0FEISg-3-rhox0t)/(2*(1-x0FEISg)*((1-etaAt)*x0FEISg-2))
      Bx0t = (rhox0t**2+rhox0t)/(2*(1-x0FEISg)*etaAt*x0FEISg)

      rhox0u = dsqrt(1+4*etaAu*(x0FEISg-1)*x0FEISg)
      Ax0u = rhox0u*(2*x0FEISg-3-rhox0u)/(2*(1-x0FEISg)*((1-etaAu)*x0FEISg-2))      
      Bx0u = (rhox0u**2+rhox0u)/(2*(1-x0FEISg)*etaAu*x0FEISg)


      gamG = 11/6d0*CV-2/3d0*TF*Nf
      GammaG = IRdiv*gamG - CV/6d0    
      Kg = CV*(67/18d0-Pi**2/6d0)-10/9d0*TF*Nf                         

! variables for IE-IS
      sbar = sman - 2*M2
      SqrtLambdaAB = kappa(sman,M2,M2)
      d2 = sbar/SqrtLambdaAB
      d1 = sbar + 2*M2 - SqrtLambdaAB
      Cdip = 2*SqrtLambdaAB/d1


! ------------ double poles -----------
      if(EpsPole.eq.-2d0) then

! ************* FE-FS *************
      colorTree = treeT3T4/CV    ! T3 T4/CA
      stst2gg_M2Dipole = -(CV*colorTree*gs**2)/(4*Pi**2)   ! CS eq. 5.34

! ************* FE-IS *************
! Squark as spectator (p1/pa)
      colorTree =  treeT1T3/CV ! T1 T3/CV
      stst2gg_M2Dipole = stst2gg_M2Dipole - (colorTree*CV*gs**2)/(4*Pi**2)   

! Antisquark as spectator (p2/pb)     
      colorTree = treeT1T4/CV ! T2 T3/CV
      stst2gg_M2Dipole = stst2gg_M2Dipole - (colorTree*CV*gs**2)/(4*Pi**2)

! ------------ simple poles -----------
      elseif(EpsPole.eq.-1d0) then
! --------------------------------------------  dipoles for 3 gluons in the final state ------------------------------------------------------
! ************* FE-FS *************
      colorTree = treeT3T4/CV    ! T3 T4/CA
      stst2gg_M2Dipole = - (colorTree*gs**2*(gamG + CV*dlog(muSc**2/sman)))/(4*Pi**2)   

! ************* FE-IS *************
! Squark as spectator (p1/pa)
      colorTree =  treeT1T3/CV ! T1 T3/CV
      stst2gg_M2Dipole = stst2gg_M2Dipole - (colorTree*gs**2*(gamG+ CV*dlog(1 + etaAt) + CV*dlog(-(muSc**2/tbar))))/(4*Pi**2)

! Antisquark as spectator (p2/pb)     
      colorTree =  treeT1T4/CV ! T2 T3/CV
      stst2gg_M2Dipole = stst2gg_M2Dipole - (colorTree*gs**2*(gamG+ CV*dlog(1 + etaAu) + CV*dlog(-(muSc**2/ubar))))/(4*Pi**2)

! ************* IE-FS *************
! emission from squark (p1/pa) with k1 and k2 as spectator (can be pulled together since T1 T3 = T2 T4)
      colorTree =  treeT1T3/CF  ! T1 T3/CF
      stst2gg_M2Dipole = stst2gg_M2Dipole + (CF*colorTree*gs**2*(-1 + dlog((1 + etaAt)/etaAt)))/(4*Pi**2)

! emission from squark (p1/pa) with k1 and k2 as spectator (can be pulled together since T1 T4 = T2 T3)
      colorTree =  treeT1T4/CF  ! T1 T4/CF
      stst2gg_M2Dipole = stst2gg_M2Dipole  + (CF*colorTree*gs**2*(-1 + dlog((1 + etaAu)/etaAu)))/(4*Pi**2)

! ************* IE-IS *************
! emission from squark and antisquark as spectator 
      colorTree = treeT1T2/CF  ! T1 T2/CF
      stst2gg_M2Dipole = stst2gg_M2Dipole + (CF*colorTree*gs**2*(-1 + d2*dlog(1 + Cdip)))/(8*Pi**2)

! emission from antisquark (p2/pb)     
      colorTree = treeT1T2/CF  ! T1 T2/CF
      stst2gg_M2Dipole = stst2gg_M2Dipole + (CF*colorTree*gs**2*(-1 + d2*dlog(1 + Cdip)))/(8*Pi**2)

! ------------ finite pieces ---------
      elseif(EpsPole.eq.0d0) then
! ************* FE-FS *************
      colorTree = treeT3T4/CV    ! T3 T4/CA
      stst2gg_M2Dipole = -(colorTree*CV*gs**2*IRdiv**2)/(4*Pi**2) - (colorTree*CV*gs**2*IRdiv*dlog(muSc**2/sman))/(4*Pi**2) + 
     &  (colorTree*gs**2*(-2*gamG - 2*GammaG - 2*Kg + CV*Pi**2 - 2*gamG*dlog(muSc**2/sman) - CV*dlog(muSc**2/sman)**2))/(8*Pi**2)

! ************* FE-IS *************
! Squark as spectator (p1/pa)
      colorTree =  treeT1T3/CV ! T1 T3/CV
      stst2gg_M2Dipole = stst2gg_M2Dipole - (colorTree*CV*gs**2*IRdiv**2)/(4*Pi**2) - (colorTree*gs**2*IRdiv*(CV*dlog(1 + etaAt) + CV*dlog(-(muSc**2/tbar))))/(4*Pi**2) + 
     &  (colorTree*gs**2*(-6*GammaG -6*dlog(muSc**2/(-tbar))*gamG - 6*Kg+ CV*Pi**2 - 12*CV*dble(cspen(dcmplx(-Ax0t))) - 3*CV*dlog(1 + etaAt)**2 + 
     &       11*CV*dlog(1 - x0FEISg) - 4*Nf*TF*dlog(1 - x0FEISg) - 6*CV*dlog(1 + etaAt)*dlog(-(muSc**2/tbar)) - 3*CV*dlog(-(muSc**2/tbar))**2))/(24*Pi**2)
 
! Antisquark as spectator (p2/pb)     
      colorTree =  treeT1T4/CV ! T2 T3/CV
      stst2gg_M2Dipole = stst2gg_M2Dipole - (colorTree*CV*gs**2*IRdiv**2)/(4*Pi**2) - (colorTree*gs**2*IRdiv*(CV*dlog(1 + etaAu) + CV*dlog(-(muSc**2/ubar))))/(4*Pi**2) + 
     &  (colorTree*gs**2*(-6*GammaG -6*dlog(muSc**2/(-ubar))*gamG - 6*Kg+ CV*Pi**2 - 12*CV*dble(cspen(dcmplx(-Ax0u))) - 3*CV*dlog(1 + etaAu)**2 + 
     &       11*CV*dlog(1 - x0FEISg) - 4*Nf*TF*dlog(1 - x0FEISg) - 6*CV*dlog(1 + etaAu)*dlog(-(muSc**2/ubar)) - 3*CV*dlog(-(muSc**2/ubar))**2))/(24*Pi**2)

! ************* IE-FS *************
! emission from squark (p1/pa) with k1 and k2 as spectator (can be pulled together since T1 T3 = T2 T4)
      colorTree =  treeT1T3/CF  ! T1 T3/CF
      stst2gg_M2Dipole = stst2gg_M2Dipole + (CF*colorTree*gs**2*IRdiv*(-1 + dlog((1 + etaAt)/etaAt)))/(4*Pi**2) + 
     &  (CF*colorTree*gs**2*(2*dble(cspen(dcmplx(-Ax0t))) - 2*dble(cspen(dcmplx(-Bx0t))) - 2*dlog(1 + Bx0t) - dlog(etaAt) - 
     &       dlog(etaAt)**2/2d0 + dlog(1 + etaAt)**2/2d0 - (1 - dlog((1 + etaAt)/etaAt))*dlog(-(muSc**2/tbar))))/(4*Pi**2)

! emission from squark (p1/pa) with k1 and k2 as spectator (can be pulled together since T1 T4 = T2 T3)
      colorTree =  treeT1T4/CF  ! T1 T4/CF
      stst2gg_M2Dipole = stst2gg_M2Dipole + (CF*colorTree*gs**2*IRdiv*(-1 + dlog((1 + etaAu)/etaAu)))/(4*Pi**2) + 
     &  (CF*colorTree*gs**2*(2*dble(cspen(dcmplx(-Ax0u))) - 2*dble(cspen(dcmplx(-Bx0u))) - 2*dlog(1 + Bx0u) - dlog(etaAu) - 
     &       dlog(etaAu)**2/2d0 + dlog(1 + etaAu)**2/2d0 - (1 - dlog((1 + etaAu)/etaAu))*dlog(-(muSc**2/ubar))))/(4*Pi**2)

! ************* IE-IS *************
! p1 as emitter and p2 as spectator and vice versa   
      colorTree = treeT1T2/CF
      stst2gg_M2Dipole = stst2gg_M2Dipole + 2*(CF*colorTree*gs**2*IRdiv*(-1 + d2*dlog(1 + Cdip)))/(8*Pi**2) + 
     &  2*(CF*colorTree*gs**2*(-(((2 + Cdip)*dlog(1 + Cdip))/Cdip) - (d2*(4*dble(cspen(dcmplx(-Cdip))) + dlog(1 + Cdip)**2))/2d0 - 
     &       (1 - d2*dlog(1 + Cdip))*dlog(muSc**2/sbar) - (1- d2*dlog(1 + Cdip))*dlog(sman/(sbar*(1 - x0IEISg)**2))))/(8*Pi**2)
      
      endif


      end