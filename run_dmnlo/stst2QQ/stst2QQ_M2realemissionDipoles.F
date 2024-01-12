      double precision function stst2QQ_M2realemissionDipoles()
      implicit none
#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"

      double precision vijkTilde,sdip,tdip,udip,yijk,vijk,zi,zj,colorTree,Q2
      double precision kappa,Qbar2,xdip,Rdip,P2ab,SqrtlambdaAB,xiab,stst2qq_M2TreeDipoleT1T4
      double precision stst2qq_M2TreeDipoleT1T2, stst2qq_M2TreeDipoleT1T3, stst2qq_M2TreeDipole

! ******************************************************************************************************************************
! ****************************************** FE-FS (from hep-ph/0201036 section 5.1) *******************************************
! ******************************************************************************************************************************

      Q2 = sqrtS**2
      vijkTilde = kappa(1d0,muu3**2,muu4**2)/(1-muu3**2-muu4**2)
      sdip = Q2

! D_{13,2}   ----> k1 as emitter and k2 as spectator
      tdip = m1**2 + m3**2 - 2*(m1**2 + p1p2 - ((m1**2 + p1p2)*(-m3**2+ m4**2 + Q2))/(2*Q2) - 
     &  ((p1k2 - ((m1**2 + p1p2)*(p1k2 + p2k2))/Q2)*kappa(Q2,m3**2,m4**2))/kappa(Q2,2*k1k3 + m3**2,m4**2)) ! tdip = (p1-\tilde{k}_{13})^2
      udip = m2**2 + m3**2 - 2*(m2**2 + p1p2 - ((m2**2 + p1p2)*(-m3**2+ m4**2 + Q2))/(2*Q2) - 
     &  ((p2k2 - ((m2**2 + p1p2)*(p1k2 + p2k2))/Q2)*kappa(Q2,m3**2,m4**2))/kappa(Q2,2*k1k3 + m3**2,m4**2)) ! udip = (p2-\tilde{k}_{13})^2

      yijk = k1k3/(k1k2 + k1k3 + k2k3)
      vijk = dsqrt(-4*muu4**2 + (2*muu4**2 + (1 - muu3**2 - muu4**2 - muu5**2)*(1 - yijk))**2)/((1 - muu3**2 - muu4**2 - muu5**2)*(1 - yijk))

      zj = k1k2/(k1k2 + k2k3)
      
      colorTree = stst2qq_M2TreeDipoleT1T2(sdip,tdip,udip)/CF
      stst2QQ_M2realemissionDipoles = (CF*colorTree*gs**2*(-((vijkTilde*(1 + m3**2/k1k3 + zj))/vijk) + 2/(1 - (1 - yijk)*zj)))/k1k3
      
! D_{23,1}   ----> k2 as emitter and k1 as spectator
      tdip = m1**2 + m3**2 - ((m1**2 + p1p2)*(m3**2 - m4**2 + Q2))/Q2 - 
     &  (2*(p1k1 - ((m1**2 + p1p2)*(p1k1 + p2k1))/Q2)*kappa(Q2,m4**2,m3**2))/kappa(Q2,2*k2k3 + m4**2,m3**2) ! tdip = (p1-\tilde{k}_1)^2
      udip =   m2**2 + m3**2 - ((m2**2 + p1p2)*(m3**2 - m4**2 + Q2))/Q2 - 
     &  (2*(p2k1 - ((m2**2 + p1p2)*(p1k1 + p2k1))/Q2)*kappa(Q2,m4**2,m3**2))/kappa(Q2,2*k2k3 + m4**2,m3**2) ! udip = (p2-\tilde{k}_1)^2

      yijk = k2k3/(k1k2 + k1k3 + k2k3)
      vijk = dsqrt(-4*muu3**2 + (2*muu3**2 + (1 - muu3**2 - muu4**2 - muu5**2)*(1 - yijk))**2)/((1 - muu3**2 - muu4**2 - muu5**2)*(1 - yijk))    

      zj = k1k2/(k1k2 + k1k3)
      
      colorTree = stst2qq_M2TreeDipoleT1T2(sdip,tdip,udip)/CF
      stst2QQ_M2realemissionDipoles = stst2QQ_M2realemissionDipoles + (CF*colorTree*gs**2*(-((vijkTilde*(1 + m4**2/k2k3 + zj))/vijk) + 2/(1 - (1 - yijk)*zj)))/k2k3
    
      
! ******************************************************************************************************************************
! ****************************************************** FE-IS and IE-FS *******************************************************
! ******************************************************************************************************************************
! dipoles with the same kinematical configuration are treated in one go
!------------------------------------------------------------------------------
      xdip = (-k1k3 + p1k1 + p1k3)/(p1k1 + p1k3)

      if (xdip.gt.x0FEIS) then

      tdip = m2**2 + m4**2 - 2*p2k2 !Q^2
      Qbar2 = tdip - m1**2 -m3**2
      Rdip = dsqrt(-4*m1**2*tdip*xdip**2 + (Qbar2 + 2*m1**2*xdip)**2)/kappa(tdip,m3**2,m1**2)
      sdip = m1**2 + m2**2 + 2*((p1p2*xdip)/Rdip + (m2**2 - p2k2)*(-(m1**2 - m3**2 + tdip)/(2*tdip) + (Qbar2 + 2*m1**2*xdip)/(2*Rdip*tdip)))
      udip = m1**2 + m4**2 - 2*((p1k2*xdip)/Rdip + (-m4**2 + p2k2)*(-(m1**2 - m3**2 + tdip)/(2*tdip) + (Qbar2 + 2*m1**2*xdip)/(2*Rdip*tdip)))

      zi = p1k3/(p1k1 + p1k3)
      zj = 1 - zi

      colorTree = stst2qq_M2TreeDipoleT1T3(sdip,tdip,udip)/CF

! D_{13}^a   ----> k1 as emitter and p1 as spectator
      stst2QQ_M2realemissionDipoles = stst2QQ_M2realemissionDipoles + (CF*colorTree*gs**2*(-2 + zi - m3**2/k1k3 + 2/(2 - xdip - zj)))/(k1k3*xdip)

! D^{a3}_1   ----> p1 as emitter and k1 as spectator 
      stst2QQ_M2realemissionDipoles = stst2QQ_M2realemissionDipoles + (CF*colorTree*gs**2*(-2 - (m1**2*xdip)/p1k3 + 2/(2 - xdip - zj)))/(p1k3*xdip)

      endif 
!------------------------------------------------------------------------------
      xdip = (-k1k3 + p2k1 + p2k3)/(p2k1 + p2k3)

      if (xdip.gt.x0FEIS) then

      udip = m1**2 + m4**2 - 2*p1k2 !Q^2
      Qbar2 = udip - m2**2 -m3**2
      Rdip = dsqrt(-4*m2**2*udip*xdip**2 + (Qbar2 + 2*m2**2*xdip)**2)/kappa(udip,m3**2,m2**2)
      sdip = m1**2 + m2**2 + 2*((p1p2*xdip)/Rdip + (m1**2 - p1k2)*(-(m2**2 - m3**2 + udip)/(2*udip) + (Qbar2 + 2*m2**2*xdip)/(2*Rdip*udip)))
      tdip = m2**2 + m4**2 - 2*((p2k2*xdip)/Rdip + (-m4**2 + p1k2)*(-(m2**2 - m3**2 + udip)/(2*udip) + (Qbar2 + 2*m2**2*xdip)/(2*Rdip*udip)))

      zi = p2k3/(p2k1 + p2k3)
      zj = 1 - zi

      colorTree = stst2qq_M2TreeDipoleT1T4(sdip,tdip,udip)/CF

! D_{13}^b   ----> k1 as emitter and p2 as spectator
      stst2QQ_M2realemissionDipoles = stst2QQ_M2realemissionDipoles + (CF*colorTree*gs**2*(-2 + zi - m3**2/k1k3  + 2/(2 - xdip - zj)))/(k1k3*xdip)

! D^{b3}_1   ----> p2 as emitter and k1 as spectator
      stst2QQ_M2realemissionDipoles = stst2QQ_M2realemissionDipoles + (CF*colorTree*gs**2*(-2 - (m2**2*xdip)/p2k3 + 2/(2 - xdip - zj)))/(p2k3*xdip)

      endif 
!------------------------------------------------------------------------------
      xdip = (-k2k3 + p1k2 + p1k3)/(p1k2 + p1k3)

      if (xdip.gt.x0FEIS) then

      udip = m2**2 + m3**2 - 2*p2k1 !Q^2
      Qbar2 = udip - m1**2 -m4**2
      Rdip = dsqrt(-4*m1**2*udip*xdip**2 + (Qbar2 + 2*m1**2*xdip)**2)/kappa(udip,m4**2,m1**2)
      sdip = m1**2 + m2**2 + 2*((p1p2*xdip)/Rdip + (m2**2 - p2k1)*(-(m1**2 - m4**2 + udip)/(2*udip) + (Qbar2 + 2*m1**2*xdip)/(2*Rdip*udip)))
      tdip = m1**2 + m3**2 - 2*((p1k1*xdip)/Rdip + (-m3**2 + p2k1)*(-(m1**2 - m4**2 + udip)/(2*udip) + (Qbar2 + 2*m1**2*xdip)/(2*Rdip*udip)))

      zi = p1k3/(p1k2 + p1k3)
      zj = 1 - zi

      colorTree = stst2qq_M2TreeDipoleT1T4(sdip,tdip,udip)/CF

! D_{23}^a   ----> k2 as emitter and p1 as spectator
      stst2QQ_M2realemissionDipoles = stst2QQ_M2realemissionDipoles + (CF*colorTree*gs**2*(-2 + zi - m4**2/k2k3 + 2/(2 - xdip - zj)))/(k2k3*xdip)

! D^{a3}_2   ----> p1 as emitter and k2 as spectator 
      stst2QQ_M2realemissionDipoles = stst2QQ_M2realemissionDipoles + (CF*colorTree*gs**2*(-2 - (m1**2*xdip)/p1k3 + 2/(2 - xdip - zj)))/(p1k3*xdip)

      endif 
!------------------------------------------------------------------------------
      xdip = (-k2k3 + p2k2 + p2k3)/(p2k2 + p2k3)

      if (xdip.gt.x0FEIS) then

      tdip = m1**2 + m3**2 - 2*p1k1 !Q^2
      Qbar2 = tdip - m2**2 -m4**2
      Rdip = dsqrt(-4*m2**2*tdip*xdip**2 + (Qbar2 + 2*m2**2*xdip)**2)/kappa(tdip,m4**2,m2**2)
      sdip = m1**2 + m2**2 + 2*((p1p2*xdip)/Rdip + (m1**2 - p1k1)*(-(m2**2 - m4**2 + tdip)/(2*tdip) + (Qbar2 + 2*m2**2*xdip)/(2*Rdip*tdip)))
      udip = m2**2 + m3**2 - 2*((p2k1*xdip)/Rdip + (-m3**2 + p1k1)*(-(m2**2 - m4**2 + tdip)/(2*tdip) + (Qbar2 + 2*m2**2*xdip)/(2*Rdip*tdip)))

      zi = p2k3/(p2k2 + p2k3)
      zj = 1 - zi

      colorTree = stst2qq_M2TreeDipoleT1T3(sdip,tdip,udip)/CF

! D_{23}^b   ----> k2 as emitter and p2 as spectator
      stst2QQ_M2realemissionDipoles = stst2QQ_M2realemissionDipoles + (CF*colorTree*gs**2*(-2 + zi - m4**2/k2k3 + 2/(2 - xdip - zj)))/(k2k3*xdip)

! D^{b3}_2   ----> p2 as emitter and k2 as spectator 
      stst2QQ_M2realemissionDipoles = stst2QQ_M2realemissionDipoles + (CF*colorTree*gs**2*(-2 - (m2**2*xdip)/p2k3 + 2/(2 - xdip - zj)))/(p2k3*xdip)

      endif 

! ******************************************************************************************************************************
! ****************************************** IE-IS  ****************************************************************************
! ******************************************************************************************************************************
      P2ab = m1**2 + m2**2 - 2*p1k3 + 2*p1p2 - 2*p2k3
      xiab = (-p1k3 + p1p2 - p2k3)/p1p2
      SqrtlambdaAB = kappa(sqrtS**2,m1**2,m2**2)

! D^{a3,b}   ----> p1 as emitter and p2 as spectator
      if (xiab.gt.x0IEIS) then
     
      sdip = P2ab
      tdip =  m2**2 + m4**2 + (2*(-m1**2 + m2**2 + P2ab)*(k2k3 - p1k2 - p2k2))/P2ab - 2*p2k2 - 
     &  ((m1**2 - 3*m2**2 - 2*p1p2 - P2ab + 2*p2k3)*((2*k2k3*m2**2 + (m1**2 - P2ab)*p2k2 - m2**2*(2*p1k2 + 3*p2k2))*SqrtlambdaAB + 
     &       (-2*m2**2*p1k2 + 2*p1p2*p2k2)*kappa(P2ab,m1**2,m2**2)))/
     &   ((-m2**4 + m1**2*(m2**2 + p1p2 - p2k3) + P2ab*(-p1p2 + p2k3) + m2**2*(-p1p2 - 3*P2ab + p2k3))*SqrtlambdaAB - 
     &     2*(m1**2*m2**2 - m2**2*p1k3 + p1p2*(-p1p2 +p2k3))*kappa(P2ab,m1**2,m2**2))        ! tdip = (p2-\tilde{k}_2)^2 = m2^2 + m4^2 - 2 p2 Lambda k2
      udip =   m2**2 + m3**2 + (2*(-m1**2 + m2**2 + P2ab)*(k1k3 -p1k1 - p2k1))/P2ab - 2*p2k1 - 
     &  ((m1**2 - 3*m2**2 - 2*p1p2 - P2ab + 2*p2k3)*((2*k1k3*m2**2 + (m1**2 - P2ab)*p2k1 - m2**2*(2*p1k1 + 3*p2k1))*SqrtlambdaAB + 
     &       (-2*m2**2*p1k1 + 2*p1p2*p2k1)*kappa(P2ab,m1**2,m2**2)))/
     &   ((-m2**4 + m1**2*(m2**2 + p1p2 - p2k3) + P2ab*(-p1p2 + p2k3) + m2**2*(-p1p2 - 3*P2ab + p2k3))*SqrtlambdaAB - 
     &     2*(m1**2*m2**2 - m2**2*p1k3 + p1p2*(-p1p2 +p2k3))*kappa(P2ab,m1**2,m2**2))        ! udip = (p2-\tilde{k}_1)^2 = m2^2 + m3^2 - 2 p2 Lambda k1


      colorTree = stst2qq_M2TreeDipoleT1T2(sdip,tdip,udip)/CF

      stst2QQ_M2realemissionDipoles = stst2QQ_M2realemissionDipoles + (CF*colorTree*gs**2*(-2 + 2/(1 - xiab) - (m1**2*xiab)/p1k3))/(p1k3*xiab)

      endif

! D^{b3,a}   ----> p2 as emitter and p1 as spectator
      if (xiab.gt.x0IEIS) then
     
      sdip = P2ab
      tdip =  m1**2 + m3**2 - 2*p1k1 + (2*(m1**2 - m2**2 + P2ab)*(k1k3- p1k1 - p2k1))/P2ab - 
     &  ((3*m1**2 - m2**2 - 2*p1k3 + 2*p1p2 + P2ab)*((2*k1k3*m1**2 + p1k1*(m2**2 - P2ab) - m1**2*(3*p1k1 + 2*p2k1))*SqrtlambdaAB + 
     &       2*(p1k1*p1p2 - m1**2*p2k1)*kappa(P2ab,m2**2,m1**2)))/
     &   ((m1**4 - m1**2*(m2**2 + p1k3 - p1p2 - 3*P2ab) + (p1k3 - p1p2)*(m2**2 - P2ab))*SqrtlambdaAB + 
     &     2*((p1k3 - p1p2)*p1p2 + m1**2*(m2**2 - p2k3))*kappa(P2ab,m2**2,m1**2))        ! tdip = (p1-\tilde{k}_1)^2 = m1^2 + m3^2 - 2 p1 Lambda k1
      udip =  m1**2 + m4**2 - 2*p1k2 + (2*(m1**2 - m2**2 + P2ab)*(k2k3- p1k2 - p2k2))/P2ab - 
     &  ((3*m1**2 - m2**2 - 2*p1k3 + 2*p1p2 + P2ab)*((2*k2k3*m1**2 + p1k2*(m2**2 - P2ab) - m1**2*(3*p1k2 + 2*p2k2))*SqrtlambdaAB + 
     &       2*(p1k2*p1p2 - m1**2*p2k2)*kappa(P2ab,m2**2,m1**2)))/
     &   ((m1**4 - m1**2*(m2**2 + p1k3 - p1p2 - 3*P2ab) + (p1k3 - p1p2)*(m2**2 - P2ab))*SqrtlambdaAB + 
     &     2*((p1k3 - p1p2)*p1p2 + m1**2*(m2**2 - p2k3))*kappa(P2ab,m2**2,m1**2))        ! udip = (p1-\tilde{k}_2)^2 = m2^2 + m4^2 - 2 p1 Lambda k2

      colorTree = stst2qq_M2TreeDipoleT1T2(sdip,tdip,udip)/CF

      stst2QQ_M2realemissionDipoles = stst2QQ_M2realemissionDipoles + (CF*colorTree*gs**2*(-2 + 2/(1 - xiab) - (m2**2*xiab)/p2k3))/(p2k3*xiab)

      endif

      end


      subroutine stst2QQ_setGauginoCoupling(k,n,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu)
      implicit none
#include "stst2QQ_Kinematics.h"
#include "stst2QQ_GenCouplings.h"
! k denotes neutralino/chargion in "normal" matrix element
! n denotes neutralino/chargino in complex conjugated matrix element
      double complex g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu  
      integer k,n

      g3Lt = ALt1(k,isf1,itt1,ftt1,igen1,fgen1)
      g3Rt = ARt1(k,isf1,itt1,ftt1,igen1,fgen1)
      g3Lct = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
      g3Rct = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))

      g4Lt = ALt1(k,isf2,itt2,ftt2,igen2,fgen2)
      g4Rt = ARt1(k,isf2,itt2,ftt2,igen2,fgen2)
      g4Lct = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
      g4Rct = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))

      g3Lu = ALt1(k,isf2,itt2,ftt1,igen2,fgen1)
      g3Ru = ARt1(k,isf2,itt2,ftt1,igen2,fgen1)
      g3Lcu = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
      g3Rcu = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))

      g4Lu = ALt1(k,isf1,itt1,ftt2,igen1,fgen2)
      g4Ru = ARt1(k,isf1,itt1,ftt2,igen1,fgen2)
      g4Lcu = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
      g4Rcu = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))

      end

      double precision function stst2qq_M2TreeDipoleT1T4(sdip,tdip,udip)
      implicit none
#include "stst2QQ_Model.h"

      double precision sdip,tdip,udip
      double precision stst2qq_M2TreeDipole,stst2qq_M2TreeDipoleT1T2,stst2qq_M2TreeDipoleT1T3

      stst2qq_M2TreeDipoleT1T4 = (-CF*stst2qq_M2TreeDipole(sdip,tdip,udip) - stst2qq_M2TreeDipoleT1T2(sdip,tdip,udip) -stst2qq_M2TreeDipoleT1T3(sdip,tdip,udip))

      end

      double precision function stst2qq_M2TreeDipoleT1T3(sdip,tdip,udip)
      implicit none

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
#include "stst2QQ_GenCouplings.h"

      integer i,k,n
      double precision sdip,tdip,udip,sqrtSdip,tredDip
      double complex g1Lt,g1Rt,g1Lct,g1Rct,g2Lt,g2Rt,g2Lct,g2Rct,g1Lu,g1Ru,g1Lcu,g1Rcu,g2Lu,g2Ru
      double complex g2Lcu,g2Rcu,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu 

    
c Dipole Mandelstam variables
      sqrtSdip = dsqrt(sdip)
      tredDip = tdip/sdip

c 	Propagator--------------------------------------------------------------------------------------

      do i=1,4 !hier gluino und neutralino unterscheiden
      DenGauT(i)=(tdip - MNeu(i)**2)! + cI*MNeu(i)*wNeu(i)
      enddo
      DenGauT(5)=(tdip - MCha(1)**2)! + cI*MCha(1)*wCha(1)
      DenGauT(6)=(tdip - MCha(2)**2)! + cI*MCha(2)*wCha(2)

      do i=1,4 !hier gluino und neutralino unterscheiden
      DenGauU(i)=(udip - MNeu(i)**2)! + cI*MNeu(i)*wNeu(i)
      enddo
      DenGauU(5)=(udip - MCha(1)**2)! + cI*MCha(1)*wCha(1)
      DenGauU(6)=(udip - MCha(2)**2)! + cI*MCha(2)*wCha(2)


      DenGlT = (tdip-MGl**2)
      DenGlU = (udip-MGl**2)

  	
      g1Lt = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
      g1Rt = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
      g1Lct = dconjg(g1Lt)
      g1Rct = dconjg(g1Rt)

      g2Lt = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
      g2Rt = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
      g2Lct = dconjg(g2Lt)
      g2Rct = dconjg(g2Rt)

      g1Lu = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
      g1Ru = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
      g1Lcu = dconjg(g1Lu)
      g1Rcu = dconjg(g1Ru)

      g2Lu = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
      g2Ru = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
      g2Lcu = dconjg(g2Lu)
      g2Rcu = dconjg(g2Ru)

      stst2qq_M2TreeDipoleT1T3 = 0d0

!************************************************   |M^{gaugino}_t|^2   ************************************************

      do k=1,6 
      do n=1,6 

      call stst2QQ_setGauginoCoupling(k,n,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu)

      stst2qq_M2TreeDipoleT1T3=stst2qq_M2TreeDipoleT1T3+((CF*Nc**2*sqrtSDip**4*(g3Rct*g3Rt*g4Lct*g4Lt*muu1
     &**2*muu2**2 + g3Lct*g3Lt*g4Rct*g4Rt*muu1**2*muu2**2 - g3Rct*g3Rt*g4Lct*g4Lt*muu2**2*muu3**2 
     &- g3Lct*g3Lt*g4Rct*g4Rt*muu2**2*muu3**2 - g3Rct*g3Rt*g4Lct*g4Lt*muu1**2*muu4**2 - g3Lct*g3Lt*g4Rct
     &*g4Rt*muu1**2*muu4**2 + g3Rct*g3Rt*g4Lct*g4Lt*muu3**2*muu4**2 + g3Lct*g3Lt*g4Rct*g4Rt*muu3
     &**2*muu4**2 + g3Rct*g3Rt*g4Lct*g4Lt*tredDip + g3Lct*g3Lt*g4Rct*g4Rt*tredDip - g3Rct*g3Rt*g4Lct
     &*g4Lt*muu1**2*tredDip - g3Lct*g3Lt*g4Rct*g4Rt*muu1**2*tredDip - g3Rct*g3Rt*g4Lct*g4Lt*muu2
     &**2*tredDip - g3Lct*g3Lt*g4Rct*g4Rt*muu2**2*tredDip + 2*g3Lct*g3Rt*g4Lt*g4Rct*muu3*muu4*tredDip 
     &+ 2*g3Lt*g3Rct*g4Lct*g4Rt*muu3*muu4*tredDip + g3Rct*g3Rt*g4Lct*g4Lt*tredDip**2 + g3Lct*g3Lt
     &*g4Rct*g4Rt*tredDip**2 + (g3Lt*g4Lt*(g3Lct*g4Rct*muu4*(-muu1**2 + muu3**2 + tredDip) + g3Rct
     &*g4Lct*muu3*(-muu2**2 + muu4**2 + tredDip)) + g3Rt*g4Rt*(g3Rct*g4Lct*muu4*(-muu1**2 + muu3**2 
     &+ tredDip) + g3Lct*g4Rct*muu3*(-muu2**2 + muu4**2 + tredDip)))*muuqt(k) + muupt(n)*(g3Lct*g4Lct
     &*(g3Lt*g4Rt*muu4*(-muu1**2 + muu3**2 + tredDip) + g3Rt*g4Lt*muu3*(-muu2**2 + muu4**2 + tredDip)) 
     &+ g3Rct*g4Rct*(g3Rt*g4Lt*muu4*(-muu1**2 + muu3**2 + tredDip) + g3Lt*g4Rt*muu3*(-muu2**2 
     &+ muu4**2 + tredDip)) + (g3Lct*g4Lct*(2*g3Rt*g4Rt*muu3*muu4 + g3Lt*g4Lt*(-1 + muu3**2 + muu4
     &**2)) + g3Rct*g4Rct*(2*g3Lt*g4Lt*muu3*muu4 + g3Rt*g4Rt*(-1 + muu3**2 + muu4**2)))*muuqt(k))))/(dconjg(DenGauT(n))*DenGauT(k)))

      enddo
      enddo


!************************************************   |M^{gluino}_t|^2   ************************************************

      stst2qq_M2TreeDipoleT1T3=stst2qq_M2TreeDipoleT1T3+(-(CF*sqrtSDip**4*(g1Rct*(g1Rt*g2Rct*muugl*(g2Rt*(
     &-1 + muu3**2 + muu4**2)*muugl + g2Lt*muu4*(-muu1**2 + muu3**2 + tredDip)) + g1Lt*g2Rct*muu3*muugl
     &*(2*g2Lt*muu4*muugl + g2Rt*(-muu2**2 + muu4**2 + tredDip)) + g1Lt*g2Lct*muu3*(2*g2Rt*muu4
     &*tredDip + g2Lt*muugl*(-muu2**2 + muu4**2 + tredDip)) + g1Rt*g2Lct*(g2Rt*muu4*muugl*(-muu1**2 
     &+ muu3**2 + tredDip) + g2Lt*(muu3**2*muu4**2 + muu1**2*(muu2**2 - muu4**2 - tredDip) + tredDip 
     &+ tredDip**2 - muu2**2*(muu3**2 + tredDip)))) + g1Lct*(g1Lt*g2Lct*muugl*(g2Lt*(-1 + muu3*
     &*2 + muu4**2)*muugl + g2Rt*muu4*(-muu1**2 + muu3**2 + tredDip)) + g1Rt*g2Lct*muu3*muugl*(2*g2Rt
     &*muu4*muugl + g2Lt*(-muu2**2 + muu4**2 + tredDip)) + g1Rt*g2Rct*muu3*(2*g2Lt*muu4*tredDip + g2Rt
     &*muugl*(-muu2**2 + muu4**2 + tredDip)) + g1Lt*g2Rct*(g2Lt*muu4*muugl*(-muu1**2 + muu3**2 
     &+ tredDip) + g2Rt*(muu3**2*muu4**2 + muu1**2*(muu2**2 - muu4**2 - tredDip) + tredDip + tredDip**2 - muu2**2*(muu3**2 + tredDip))))))/(4.*DenGlT*dconjg(DenGlT)))


!************************************************   2 Re(M^{gluino}_t M^{gaugino ast}_t)   ************************************************

      do n=1,6 

      call stst2QQ_setGauginoCoupling(1,n,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu)

      stst2qq_M2TreeDipoleT1T3=stst2qq_M2TreeDipoleT1T3+dble(0)

      enddo


!************************************************   |M^{gaugino}_u|^2   ************************************************

      do k=1,6 
      do n=1,6 

      call stst2QQ_setGauginoCoupling(k,n,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu)

      stst2qq_M2TreeDipoleT1T3=stst2qq_M2TreeDipoleT1T3+(0)

      enddo
      enddo


!************************************************   |M^{gluino}_u|^2   ************************************************

      stst2qq_M2TreeDipoleT1T3=stst2qq_M2TreeDipoleT1T3+((CF*(-2 + Nc**2)*sqrtSDip**4*(g1Rcu*(g1Ru*g2Rcu*muugl
     &*(g2Ru*(-1 + muu3**2 + muu4**2)*muugl + g2Lu*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip)) 
     &+ g1Lu*g2Rcu*muu3*muugl*(2*g2Lu*muu4*muugl + g2Ru*(-1 + muu2**2 + muu3**2 + 2*muu4*
     &*2 - tredDip)) + g1Lu*g2Lcu*muu3*(2*g2Ru*muu4*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tredDip) 
     &+ g2Lu*muugl*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + g1Ru*g2Lcu*(g2Ru*muu4*muugl
     &*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + g2Lu*(muu3**4 - muu4**2 + muu4**4 + muu3
     &**2*(-1 + muu2**2 + 3*muu4**2 - 2*tredDip) + muu1**2*(muu2**2 + muu4**2 - tredDip) + tredDip 
     &- muu2**2*tredDip - 2*muu4**2*tredDip + tredDip**2))) + g1Lcu*(g1Lu*g2Lcu*muugl*(g2Lu*(-1 + muu3
     &**2 + muu4**2)*muugl + g2Ru*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip)) + g1Ru*g2Ru
     &*muu3*muugl*(2*g2Lcu*muu4*muugl + g2Rcu*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + g1Ru
     &*g2Lu*muu3*(2*g2Rcu*muu4*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tredDip) + g2Lcu*muugl
     &*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + g1Lu*g2Rcu*(g2Lu*muu4*muugl*(-1 + muu1
     &**2 + 2*muu3**2 + muu4**2 - tredDip) + g2Ru*(muu3**4 - muu4**2 + muu4**4 + muu3**2*(-1 + muu2
     &**2 + 3*muu4**2 - 2*tredDip) + muu1**2*(muu2**2 + muu4**2 - tredDip) + tredDip - muu2**2*tredDip - 2*muu4**2*tredDip + tredDip**2)))))/(4.*DenGlU*dconjg(DenGlU)))


!************************************************   2 Re(M^{gluino}_u M^{gaugino ast}_u)   ************************************************

      do n=1,6 

      call stst2QQ_setGauginoCoupling(1,n,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu)

      stst2qq_M2TreeDipoleT1T3=stst2qq_M2TreeDipoleT1T3+dble((CF*Nc*sqrtSDip**4*(g1Ru*(2*g2Lu*g3Lcu*g4Rcu
     &*muu3*muu4*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tredDip) + g2Ru*g3Rcu*g4Lcu*muu4*muugl
     &*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + g2Ru*g3Lcu*g4Rcu*muu3*muugl*(-1 + muu2**2 
     &+ muu3**2 + 2*muu4**2 - tredDip) + g2Lu*g3Rcu*g4Lcu*(muu3**4 - muu4**2 + muu4**4 + muu3**2*(
     &-1 + muu2**2 + 3*muu4**2 - 2*tredDip) + muu1**2*(muu2**2 + muu4**2 - tredDip) + tredDip - muu2
     &**2*tredDip - 2*muu4**2*tredDip + tredDip**2)) + g1Lu*(2*g2Ru*g3Rcu*g4Lcu*muu3*muu4*(-1 + muu1
     &**2 + muu2**2 + muu3**2 + muu4**2 - tredDip) + g2Lu*g3Lcu*g4Rcu*muu4*muugl*(-1 + muu1**2 + 2
     &*muu3**2 + muu4**2 - tredDip) + g2Lu*g3Rcu*g4Lcu*muu3*muugl*(-1 + muu2**2 + muu3**2 + 2*muu4
     &**2 - tredDip) + g2Ru*g3Lcu*g4Rcu*(muu3**4 - muu4**2 + muu4**4 + muu3**2*(-1 + muu2**2 + 3*muu4
     &**2 - 2*tredDip) + muu1**2*(muu2**2 + muu4**2 - tredDip) + tredDip - muu2**2*tredDip - 2*muu4
     &**2*tredDip + tredDip**2)) + (g1Ru*(2*g2Ru*g3Lcu*g4Lcu*muu3*muu4*muugl + g2Ru*g3Rcu*g4Rcu*(
     &-1 + muu3**2 + muu4**2)*muugl + g2Lu*g3Rcu*g4Rcu*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) 
     &+ g2Lu*g3Lcu*g4Lcu*muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + g1Lu*(2*g2Lu
     &*g3Rcu*g4Rcu*muu3*muu4*muugl + g2Lu*g3Lcu*g4Lcu*(-1 + muu3**2 + muu4**2)*muugl + g2Ru*g3Lcu*g4Lcu
     &*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + g2Ru*g3Rcu*g4Rcu*muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)))*muupu(n)))/(DenGlU*dconjg(DenGauU(n))))

      enddo


!************************************************   2 Re(M^{gaugino}_t M^{gaugino ast}_u)   ************************************************

      do k=1,6 
      do n=1,6 

      call stst2QQ_setGauginoCoupling(k,n,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu)

      stst2qq_M2TreeDipoleT1T3=stst2qq_M2TreeDipoleT1T3+dble((-2*CF*Nc*sqrtSDip**4*(g3Rcu*g3Rt*g4Lcu*g4Lt
     &*muu1**2*muu2**2 + g3Lcu*g3Lt*g4Rcu*g4Rt*muu1**2*muu2**2 + g3Lcu*g3Rt*g4Lt*g4Rcu*muu1**2*muu3
     &*muu4 + g3Lt*g3Rcu*g4Lcu*g4Rt*muu1**2*muu3*muu4 + g3Lcu*g3Rt*g4Lt*g4Rcu*muu2**2*muu3*muu4 + g3Lt
     &*g3Rcu*g4Lcu*g4Rt*muu2**2*muu3*muu4 - g3Lcu*g3Rt*g4Lt*g4Rcu*muu3**3*muu4 - g3Lt*g3Rcu*g4Lcu
     &*g4Rt*muu3**3*muu4 - g3Rcu*g3Rt*g4Lcu*g4Lt*muu3**2*muu4**2 - g3Lcu*g3Lt*g4Rcu*g4Rt*muu3**2*muu4
     &**2 - g3Lcu*g3Rt*g4Lt*g4Rcu*muu3*muu4**3 - g3Lt*g3Rcu*g4Lcu*g4Rt*muu3*muu4**3 + g3Rcu*g3Rt
     &*g4Lcu*g4Lt*tredDip + g3Lcu*g3Lt*g4Rcu*g4Rt*tredDip - g3Rcu*g3Rt*g4Lcu*g4Lt*muu1**2*tredDip - g3Lcu
     &*g3Lt*g4Rcu*g4Rt*muu1**2*tredDip - g3Rcu*g3Rt*g4Lcu*g4Lt*muu2**2*tredDip - g3Lcu*g3Lt*g4Rcu
     &*g4Rt*muu2**2*tredDip - g3Rcu*g3Rt*g4Lcu*g4Lt*muu3**2*tredDip - g3Lcu*g3Lt*g4Rcu*g4Rt*muu3
     &**2*tredDip - g3Rcu*g3Rt*g4Lcu*g4Lt*muu4**2*tredDip - g3Lcu*g3Lt*g4Rcu*g4Rt*muu4**2*tredDip 
     &+ g3Rcu*g3Rt*g4Lcu*g4Lt*tredDip**2 + g3Lcu*g3Lt*g4Rcu*g4Rt*tredDip**2 - (g3Lt*g4Lt*(g3Lcu*g4Rcu
     &*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + g3Rcu*g4Lcu*muu3*(-1 + muu2**2 + muu3
     &**2 + 2*muu4**2 - tredDip)) + g3Rt*g4Rt*(g3Rcu*g4Lcu*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 
     &- tredDip) + g3Lcu*g4Rcu*muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)))*muuqt(k) + muupu(n)
     &*(g3Lcu*g4Lcu*(g3Lt*g4Rt*muu4*(muu1**2 - muu3**2 - tredDip) + g3Rt*g4Lt*muu3*(muu2**2 
     &- muu4**2 - tredDip)) + g3Rcu*g4Rcu*(g3Rt*g4Lt*muu4*(muu1**2 - muu3**2 - tredDip) + g3Lt*g4Rt
     &*muu3*(muu2**2 - muu4**2 - tredDip)) - (g3Lcu*g4Lcu*(2*g3Rt*g4Rt*muu3*muu4 + g3Lt*g4Lt*(-1 + muu3
     &**2 + muu4**2)) + g3Rcu*g4Rcu*(2*g3Lt*g4Lt*muu3*muu4 + g3Rt*g4Rt*(-1 + muu3**2 + muu4**2)))*muuqt(k))))/(dconjg(DenGauU(n))*DenGauT(k)))

      enddo
      enddo


!************************************************   2 Re(M^{gluino}_t M^{gluino ast}_u)   ************************************************

      stst2qq_M2TreeDipoleT1T3=stst2qq_M2TreeDipoleT1T3+((CF*(-2*CF + Nc)*sqrtSDip**4*(g1Rcu*(g1Lt*g2Lcu*muu3
     &*(g2Rt*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2) + g2Lt*muugl*(-1 + muu2**2 + muu3**2 
     &+ 2*muu4**2 - tredDip)) + g1Rt*g2Rcu*muugl*(g2Rt*(-1 + muu3**2 + muu4**2)*muugl + g2Lt*muu4*(
     &-muu1**2 + muu3**2 + tredDip)) + g1Lt*g2Rcu*muu3*muugl*(2*g2Lt*muu4*muugl + g2Rt*(-muu2**2 + muu4
     &**2 + tredDip)) + g1Rt*g2Lcu*(g2Rt*muu4*muugl*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) 
     &+ g2Lt*((-1 + muu2**2 + muu4**2 - tredDip)*tredDip + muu1**2*(-muu2**2 + tredDip) + muu3
     &**2*(muu4**2 + tredDip)))) + g1Lcu*(g1Rt*g2Rt*muu3*muugl*(2*g2Lcu*muu4*muugl + g2Rcu*(-1 + muu2
     &**2 + muu3**2 + 2*muu4**2 - tredDip)) + g1Lt*g2Lcu*muugl*(g2Lt*(-1 + muu3**2 + muu4**2)*muugl 
     &+ g2Rt*muu4*(-muu1**2 + muu3**2 + tredDip)) + g1Rt*g2Lt*muu3*(g2Rcu*muu4*(-muu1**2 - muu2*
     &*2 + muu3**2 + muu4**2) + g2Lcu*muugl*(-muu2**2 + muu4**2 + tredDip)) + g1Lt*g2Rcu*(g2Lt*muu4
     &*muugl*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + g2Rt*((-1 + muu2**2 + muu4**2 - tredDip)
     &*tredDip + muu1**2*(-muu2**2 + tredDip) + muu3**2*(muu4**2 + tredDip))))))/(2.*DenGlT*dconjg(DenGlU)))


!************************************************   2 Re(M^{gluino}_t M^{gaugino ast}_u)   ************************************************

      do n=1,6 

      call stst2QQ_setGauginoCoupling(1,n,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu)

      stst2qq_M2TreeDipoleT1T3=stst2qq_M2TreeDipoleT1T3+dble(-((CF*sqrtSDip**4*(g1Rt*(g2Lt*g3Lcu*g4Rcu*muu3
     &*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2) + g2Rt*g3Rcu*g4Lcu*muu4*muugl*(-1 + muu1**2 
     &+ 2*muu3**2 + muu4**2 - tredDip) + g2Rt*g3Lcu*g4Rcu*muu3*muugl*(-1 + muu2**2 + muu3**2 + 2*muu4
     &**2 - tredDip) + g2Lt*g3Rcu*g4Lcu*((-1 + muu2**2 + muu4**2 - tredDip)*tredDip + muu1**2*(-muu2
     &**2 + tredDip) + muu3**2*(muu4**2 + tredDip))) + g1Lt*(g2Rt*g3Rcu*g4Lcu*muu3*muu4*(-muu1**2 
     &- muu2**2 + muu3**2 + muu4**2) + g2Lt*g3Lcu*g4Rcu*muu4*muugl*(-1 + muu1**2 + 2*muu3**2 + muu4
     &**2 - tredDip) + g2Lt*g3Rcu*g4Lcu*muu3*muugl*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip) 
     &+ g2Rt*g3Lcu*g4Rcu*((-1 + muu2**2 + muu4**2 - tredDip)*tredDip + muu1**2*(-muu2**2 + tredDip) 
     &+ muu3**2*(muu4**2 + tredDip))) + (g1Rt*(2*g2Rt*g3Lcu*g4Lcu*muu3*muu4*muugl + g2Rt*g3Rcu*g4Rcu
     &*(-1 + muu3**2 + muu4**2)*muugl + g2Lt*g3Rcu*g4Rcu*muu4*(-muu1**2 + muu3**2 + tredDip) + g2Lt
     &*g3Lcu*g4Lcu*muu3*(-muu2**2 + muu4**2 + tredDip)) + g1Lt*(2*g2Lt*g3Rcu*g4Rcu*muu3*muu4*muugl 
     &+ g2Lt*g3Lcu*g4Lcu*(-1 + muu3**2 + muu4**2)*muugl + g2Rt*g3Lcu*g4Lcu*muu4*(-muu1**2 + muu3
     &**2 + tredDip) + g2Rt*g3Rcu*g4Rcu*muu3*(-muu2**2 + muu4**2 + tredDip)))*muupu(n)))/(DenGlT*dconjg(DenGauU(n)))))

      enddo


!************************************************   2 Re(M^{gluino}_u M^{gaugino ast}_t)   ************************************************

      do n=1,6 

      call stst2QQ_setGauginoCoupling(1,n,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu)

      stst2qq_M2TreeDipoleT1T3=stst2qq_M2TreeDipoleT1T3+dble((2*CF**2*Nc*sqrtSDip**4*(g1Ru*(g2Lu*g3Lct*g4Rct
     &*muu3*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2) + g2Ru*g3Rct*g4Lct*muu4*muugl*(-muu1**2 
     &+ muu3**2 + tredDip) + g2Ru*g3Lct*g4Rct*muu3*muugl*(-muu2**2 + muu4**2 + tredDip) + g2Lu*g3Rct
     &*g4Lct*((-1 + muu2**2 + muu4**2 - tredDip)*tredDip + muu1**2*(-muu2**2 + tredDip) + muu3**2
     &*(muu4**2 + tredDip))) + g1Lu*(g2Ru*g3Rct*g4Lct*muu3*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4
     &**2) + g2Lu*g3Lct*g4Rct*muu4*muugl*(-muu1**2 + muu3**2 + tredDip) + g2Lu*g3Rct*g4Lct*muu3*muugl
     &*(-muu2**2 + muu4**2 + tredDip) + g2Ru*g3Lct*g4Rct*((-1 + muu2**2 + muu4**2 - tredDip)*tredDip 
     &+ muu1**2*(-muu2**2 + tredDip) + muu3**2*(muu4**2 + tredDip))) + (g1Ru*(2*g2Ru*g3Lct*g4Lct
     &*muu3*muu4*muugl + g2Ru*g3Rct*g4Rct*(-1 + muu3**2 + muu4**2)*muugl + g2Lu*g3Rct*g4Rct*muu4
     &*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + g2Lu*g3Lct*g4Lct*muu3*(-1 + muu2**2 + muu3*
     &*2 + 2*muu4**2 - tredDip)) + g1Lu*(2*g2Lu*g3Rct*g4Rct*muu3*muu4*muugl + g2Lu*g3Lct*g4Lct*(-1 
     &+ muu3**2 + muu4**2)*muugl + g2Ru*g3Lct*g4Lct*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) 
     &+ g2Ru*g3Rct*g4Rct*muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)))*muupt(n)))/(DenGlU*dconjg(DenGauT(n))))

      enddo
      end

      double precision function stst2qq_M2TreeDipoleT1T2(sdip,tdip,udip)
      implicit none

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
#include "stst2QQ_GenCouplings.h"

      integer i,k,n
      double precision sdip,tdip,udip,sqrtSdip,tredDip
      double complex g1Lt,g1Rt,g1Lct,g1Rct,g2Lt,g2Rt,g2Lct,g2Rct,g1Lu,g1Ru,g1Lcu,g1Rcu,g2Lu,g2Ru
      double complex g2Lcu,g2Rcu,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu 

    
c Dipole Mandelstam variables
      sqrtSdip = dsqrt(sdip)
      tredDip = tdip/sdip

c 	Propagator--------------------------------------------------------------------------------------

      do i=1,4 !hier gluino und neutralino unterscheiden
      DenGauT(i)=(tdip - MNeu(i)**2)! + cI*MNeu(i)*wNeu(i)
      enddo
      DenGauT(5)=(tdip - MCha(1)**2)! + cI*MCha(1)*wCha(1)
      DenGauT(6)=(tdip - MCha(2)**2)! + cI*MCha(2)*wCha(2)

      do i=1,4 !hier gluino und neutralino unterscheiden
      DenGauU(i)=(udip - MNeu(i)**2)! + cI*MNeu(i)*wNeu(i)
      enddo
      DenGauU(5)=(udip - MCha(1)**2)! + cI*MCha(1)*wCha(1)
      DenGauU(6)=(udip - MCha(2)**2)! + cI*MCha(2)*wCha(2)


      DenGlT = (tdip-MGl**2)
      DenGlU = (udip-MGl**2)

  	
      g1Lt = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
      g1Rt = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
      g1Lct = dconjg(g1Lt)
      g1Rct = dconjg(g1Rt)

      g2Lt = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
      g2Rt = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
      g2Lct = dconjg(g2Lt)
      g2Rct = dconjg(g2Rt)

      g1Lu = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
      g1Ru = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
      g1Lcu = dconjg(g1Lu)
      g1Rcu = dconjg(g1Ru)

      g2Lu = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
      g2Ru = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
      g2Lcu = dconjg(g2Lu)
      g2Rcu = dconjg(g2Ru)

      stst2qq_M2TreeDipoleT1T2 = 0d0

!************************************************   |M^{gaugino}_t|^2   ************************************************

      do k=1,6 
      do n=1,6 

      call stst2QQ_setGauginoCoupling(k,n,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu)

      stst2qq_M2TreeDipoleT1T2=stst2qq_M2TreeDipoleT1T2+(0)

      enddo
      enddo


!************************************************   |M^{gluino}_t|^2   ************************************************

      stst2qq_M2TreeDipoleT1T2=stst2qq_M2TreeDipoleT1T2+((CF*sqrtSDip**4*(g1Rct*(g1Rt*g2Rct*muugl*(g2Rt*(
     &-1 + muu3**2 + muu4**2)*muugl + g2Lt*muu4*(-muu1**2 + muu3**2 + tredDip)) + g1Lt*g2Rct*muu3*muugl
     &*(2*g2Lt*muu4*muugl + g2Rt*(-muu2**2 + muu4**2 + tredDip)) + g1Lt*g2Lct*muu3*(2*g2Rt*muu4
     &*tredDip + g2Lt*muugl*(-muu2**2 + muu4**2 + tredDip)) + g1Rt*g2Lct*(g2Rt*muu4*muugl*(-muu1**2 
     &+ muu3**2 + tredDip) + g2Lt*(muu3**2*muu4**2 + muu1**2*(muu2**2 - muu4**2 - tredDip) + tredDip 
     &+ tredDip**2 - muu2**2*(muu3**2 + tredDip)))) + g1Lct*(g1Lt*g2Lct*muugl*(g2Lt*(-1 + muu3**2 
     &+ muu4**2)*muugl + g2Rt*muu4*(-muu1**2 + muu3**2 + tredDip)) + g1Rt*g2Lct*muu3*muugl*(2*g2Rt
     &*muu4*muugl + g2Lt*(-muu2**2 + muu4**2 + tredDip)) + g1Rt*g2Rct*muu3*(2*g2Lt*muu4*tredDip + g2Rt
     &*muugl*(-muu2**2 + muu4**2 + tredDip)) + g1Lt*g2Rct*(g2Lt*muu4*muugl*(-muu1**2 + muu3**2 
     &+ tredDip) + g2Rt*(muu3**2*muu4**2 + muu1**2*(muu2**2 - muu4**2 - tredDip) + tredDip + tredDip**2 - muu2**2*(muu3**2 + tredDip))))))/(2.*DenGlT*dconjg(DenGlT)))


!************************************************   2 Re(M^{gluino}_t M^{gaugino ast}_t)   ************************************************

      do n=1,6 

      call stst2QQ_setGauginoCoupling(1,n,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu)

      stst2qq_M2TreeDipoleT1T2=stst2qq_M2TreeDipoleT1T2+dble(-((CF*Nc*sqrtSDip**4*(g1Rt*(2*g2Lt*g3Lct*g4Rct
     &*muu3*muu4*tredDip + g2Rt*g3Rct*g4Lct*muu4*muugl*(-muu1**2 + muu3**2 + tredDip) + g2Rt*g3Lct
     &*g4Rct*muu3*muugl*(-muu2**2 + muu4**2 + tredDip) + g2Lt*g3Rct*g4Lct*(muu3**2*muu4**2 + muu1*
     &*2*(muu2**2 - muu4**2 - tredDip) + tredDip + tredDip**2 - muu2**2*(muu3**2 + tredDip))) + g1Lt
     &*(2*g2Rt*g3Rct*g4Lct*muu3*muu4*tredDip + g2Lt*g3Lct*g4Rct*muu4*muugl*(-muu1**2 + muu3**2 + tredDip) 
     &+ g2Lt*g3Rct*g4Lct*muu3*muugl*(-muu2**2 + muu4**2 + tredDip) + g2Rt*g3Lct*g4Rct*(muu3
     &**2*muu4**2 + muu1**2*(muu2**2 - muu4**2 - tredDip) + tredDip + tredDip**2 - muu2**2*(muu3**2 
     &+ tredDip))) + (g1Rt*(2*g2Rt*g3Lct*g4Lct*muu3*muu4*muugl + g2Rt*g3Rct*g4Rct*(-1 + muu3**2 + muu4
     &**2)*muugl + g2Lt*g3Rct*g4Rct*muu4*(-muu1**2 + muu3**2 + tredDip) + g2Lt*g3Lct*g4Lct*muu3
     &*(-muu2**2 + muu4**2 + tredDip)) + g1Lt*(2*g2Lt*g3Rct*g4Rct*muu3*muu4*muugl + g2Lt*g3Lct*g4Lct
     &*(-1 + muu3**2 + muu4**2)*muugl + g2Rt*g3Lct*g4Lct*muu4*(-muu1**2 + muu3**2 + tredDip) + g2Rt
     &*g3Rct*g4Rct*muu3*(-muu2**2 + muu4**2 + tredDip)))*muupt(n)))/(DenGlT*dconjg(DenGauT(n)))))

      enddo


!************************************************   |M^{gaugino}_u|^2   ************************************************

      do k=1,6 
      do n=1,6 

      call stst2QQ_setGauginoCoupling(k,n,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu)

      stst2qq_M2TreeDipoleT1T2=stst2qq_M2TreeDipoleT1T2+(0)

      enddo
      enddo


!************************************************   |M^{gluino}_u|^2   ************************************************

      stst2qq_M2TreeDipoleT1T2=stst2qq_M2TreeDipoleT1T2+((CF*sqrtSDip**4*(g1Rcu*(g1Ru*g2Rcu*muugl*(g2Ru*(
     &-1 + muu3**2 + muu4**2)*muugl + g2Lu*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip)) + g1Lu
     &*g2Rcu*muu3*muugl*(2*g2Lu*muu4*muugl + g2Ru*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) 
     &+ g1Lu*g2Lcu*muu3*(2*g2Ru*muu4*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tredDip) + g2Lu
     &*muugl*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + g1Ru*g2Lcu*(g2Ru*muu4*muugl*(-1 + muu1
     &**2 + 2*muu3**2 + muu4**2 - tredDip) + g2Lu*(muu3**4 - muu4**2 + muu4**4 + muu3**2*(-1 + muu2
     &**2 + 3*muu4**2 - 2*tredDip) + muu1**2*(muu2**2 + muu4**2 - tredDip) + tredDip - muu2**2*tredDip 
     &- 2*muu4**2*tredDip + tredDip**2))) + g1Lcu*(g1Lu*g2Lcu*muugl*(g2Lu*(-1 + muu3**2 + muu4
     &**2)*muugl + g2Ru*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip)) + g1Ru*g2Ru*muu3*muugl
     &*(2*g2Lcu*muu4*muugl + g2Rcu*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + g1Ru*g2Lu*muu3
     &*(2*g2Rcu*muu4*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tredDip) + g2Lcu*muugl*(-1 + muu2
     &**2 + muu3**2 + 2*muu4**2 - tredDip)) + g1Lu*g2Rcu*(g2Lu*muu4*muugl*(-1 + muu1**2 + 2*muu3
     &**2 + muu4**2 - tredDip) + g2Ru*(muu3**4 - muu4**2 + muu4**4 + muu3**2*(-1 + muu2**2 + 3*muu4
     &**2 - 2*tredDip) + muu1**2*(muu2**2 + muu4**2 - tredDip) + tredDip - muu2**2*tredDip - 2*muu4**2*tredDip + tredDip**2)))))/(2.*DenGlU*dconjg(DenGlU)))


!************************************************   2 Re(M^{gluino}_u M^{gaugino ast}_u)   ************************************************

      do n=1,6 

      call stst2QQ_setGauginoCoupling(1,n,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu)

      stst2qq_M2TreeDipoleT1T2=stst2qq_M2TreeDipoleT1T2+dble(-((CF*Nc*sqrtSDip**4*(g1Ru*(2*g2Lu*g3Lcu*g4Rcu
     &*muu3*muu4*(-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tredDip) + g2Ru*g3Rcu*g4Lcu*muu4*muugl
     &*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + g2Ru*g3Lcu*g4Rcu*muu3*muugl*(-1 + muu2*
     &*2 + muu3**2 + 2*muu4**2 - tredDip) + g2Lu*g3Rcu*g4Lcu*(muu3**4 - muu4**2 + muu4**4 + muu3**2
     &*(-1 + muu2**2 + 3*muu4**2 - 2*tredDip) + muu1**2*(muu2**2 + muu4**2 - tredDip) + tredDip - muu2
     &**2*tredDip - 2*muu4**2*tredDip + tredDip**2)) + g1Lu*(2*g2Ru*g3Rcu*g4Lcu*muu3*muu4*(-1 + muu1
     &**2 + muu2**2 + muu3**2 + muu4**2 - tredDip) + g2Lu*g3Lcu*g4Rcu*muu4*muugl*(-1 + muu1**2 
     &+ 2*muu3**2 + muu4**2 - tredDip) + g2Lu*g3Rcu*g4Lcu*muu3*muugl*(-1 + muu2**2 + muu3**2 + 2*muu4
     &**2 - tredDip) + g2Ru*g3Lcu*g4Rcu*(muu3**4 - muu4**2 + muu4**4 + muu3**2*(-1 + muu2**2 + 3*muu4
     &**2 - 2*tredDip) + muu1**2*(muu2**2 + muu4**2 - tredDip) + tredDip - muu2**2*tredDip - 2*muu4
     &**2*tredDip + tredDip**2)) + (g1Ru*(2*g2Ru*g3Lcu*g4Lcu*muu3*muu4*muugl + g2Ru*g3Rcu*g4Rcu
     &*(-1 + muu3**2 + muu4**2)*muugl + g2Lu*g3Rcu*g4Rcu*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) 
     &+ g2Lu*g3Lcu*g4Lcu*muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + g1Lu*(2*g2Lu
     &*g3Rcu*g4Rcu*muu3*muu4*muugl + g2Lu*g3Lcu*g4Lcu*(-1 + muu3**2 + muu4**2)*muugl + g2Ru*g3Lcu
     &*g4Lcu*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + g2Ru*g3Rcu*g4Rcu*muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)))*muupu(n)))/(DenGlU*dconjg(DenGauU(n)))))

      enddo


!************************************************   2 Re(M^{gaugino}_t M^{gaugino ast}_u)   ************************************************

      do k=1,6 
      do n=1,6 

      call stst2QQ_setGauginoCoupling(k,n,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu)

      stst2qq_M2TreeDipoleT1T2=stst2qq_M2TreeDipoleT1T2+dble((2*CF*Nc*sqrtSDip**4*(g3Rcu*g3Rt*g4Lcu*g4Lt*muu1
     &**2*muu2**2 + g3Lcu*g3Lt*g4Rcu*g4Rt*muu1**2*muu2**2 + g3Lcu*g3Rt*g4Lt*g4Rcu*muu1**2*muu3*muu4 
     &+ g3Lt*g3Rcu*g4Lcu*g4Rt*muu1**2*muu3*muu4 + g3Lcu*g3Rt*g4Lt*g4Rcu*muu2**2*muu3*muu4 + g3Lt
     &*g3Rcu*g4Lcu*g4Rt*muu2**2*muu3*muu4 - g3Lcu*g3Rt*g4Lt*g4Rcu*muu3**3*muu4 - g3Lt*g3Rcu*g4Lcu
     &*g4Rt*muu3**3*muu4 - g3Rcu*g3Rt*g4Lcu*g4Lt*muu3**2*muu4**2 - g3Lcu*g3Lt*g4Rcu*g4Rt*muu3**2*muu4
     &**2 - g3Lcu*g3Rt*g4Lt*g4Rcu*muu3*muu4**3 - g3Lt*g3Rcu*g4Lcu*g4Rt*muu3*muu4**3 + g3Rcu*g3Rt*g4Lcu
     &*g4Lt*tredDip + g3Lcu*g3Lt*g4Rcu*g4Rt*tredDip - g3Rcu*g3Rt*g4Lcu*g4Lt*muu1**2*tredDip - g3Lcu
     &*g3Lt*g4Rcu*g4Rt*muu1**2*tredDip - g3Rcu*g3Rt*g4Lcu*g4Lt*muu2**2*tredDip - g3Lcu*g3Lt*g4Rcu
     &*g4Rt*muu2**2*tredDip - g3Rcu*g3Rt*g4Lcu*g4Lt*muu3**2*tredDip - g3Lcu*g3Lt*g4Rcu*g4Rt*muu3
     &**2*tredDip - g3Rcu*g3Rt*g4Lcu*g4Lt*muu4**2*tredDip - g3Lcu*g3Lt*g4Rcu*g4Rt*muu4**2*tredDip + g3Rcu
     &*g3Rt*g4Lcu*g4Lt*tredDip**2 + g3Lcu*g3Lt*g4Rcu*g4Rt*tredDip**2 - (g3Lt*g4Lt*(g3Lcu*g4Rcu
     &*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + g3Rcu*g4Lcu*muu3*(-1 + muu2**2 + muu3
     &**2 + 2*muu4**2 - tredDip)) + g3Rt*g4Rt*(g3Rcu*g4Lcu*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 
     &- tredDip) + g3Lcu*g4Rcu*muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)))*muuqt(k) + muupu(n)
     &*(g3Lcu*g4Lcu*(g3Lt*g4Rt*muu4*(muu1**2 - muu3**2 - tredDip) + g3Rt*g4Lt*muu3*(muu2**2 - muu4
     &**2 - tredDip)) + g3Rcu*g4Rcu*(g3Rt*g4Lt*muu4*(muu1**2 - muu3**2 - tredDip) + g3Lt*g4Rt*muu3
     &*(muu2**2 - muu4**2 - tredDip)) - (g3Lcu*g4Lcu*(2*g3Rt*g4Rt*muu3*muu4 + g3Lt*g4Lt*(-1 + muu3
     &**2 + muu4**2)) + g3Rcu*g4Rcu*(2*g3Lt*g4Lt*muu3*muu4 + g3Rt*g4Rt*(-1 + muu3**2 + muu4**2)))*muuqt(k))))/(dconjg(DenGauU(n))*DenGauT(k)))

      enddo
      enddo


!************************************************   2 Re(M^{gluino}_t M^{gluino ast}_u)   ************************************************

      stst2qq_M2TreeDipoleT1T2=stst2qq_M2TreeDipoleT1T2+(-(CF*(-2*CF + Nc)*(1 + Nc**2)*sqrtSDip**4*(g1Rcu
     &*(g1Lt*g2Lcu*muu3*(g2Rt*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2) + g2Lt*muugl*(-1 + muu2
     &**2 + muu3**2 + 2*muu4**2 - tredDip)) + g1Rt*g2Rcu*muugl*(g2Rt*(-1 + muu3**2 + muu4**2)*muugl 
     &+ g2Lt*muu4*(-muu1**2 + muu3**2 + tredDip)) + g1Lt*g2Rcu*muu3*muugl*(2*g2Lt*muu4*muugl + g2Rt
     &*(-muu2**2 + muu4**2 + tredDip)) + g1Rt*g2Lcu*(g2Rt*muu4*muugl*(-1 + muu1**2 + 2*muu3**2 + muu4
     &**2 - tredDip) + g2Lt*((-1 + muu2**2 + muu4**2 - tredDip)*tredDip + muu1**2*(-muu2**2 + tredDip) 
     &+ muu3**2*(muu4**2 + tredDip)))) + g1Lcu*(g1Rt*g2Rt*muu3*muugl*(2*g2Lcu*muu4*muugl + g2Rcu
     &*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + g1Lt*g2Lcu*muugl*(g2Lt*(-1 + muu3**2 + muu4
     &**2)*muugl + g2Rt*muu4*(-muu1**2 + muu3**2 + tredDip)) + g1Rt*g2Lt*muu3*(g2Rcu*muu4*(-muu1
     &**2 - muu2**2 + muu3**2 + muu4**2) + g2Lcu*muugl*(-muu2**2 + muu4**2 + tredDip)) + g1Lt*g2Rcu
     &*(g2Lt*muu4*muugl*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + g2Rt*((-1 + muu2**2 + muu4
     &**2 - tredDip)*tredDip + muu1**2*(-muu2**2 + tredDip) + muu3**2*(muu4**2 + tredDip))))))/(2.*DenGlT*dconjg(DenGlU)))


!************************************************   2 Re(M^{gluino}_t M^{gaugino ast}_u)   ************************************************

      do n=1,6 

      call stst2QQ_setGauginoCoupling(1,n,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu)

      stst2qq_M2TreeDipoleT1T2=stst2qq_M2TreeDipoleT1T2+dble((CF*sqrtSDip**4*(g1Rt*(g2Lt*g3Lcu*g4Rcu*muu3
     &*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2) + g2Rt*g3Rcu*g4Lcu*muu4*muugl*(-1 + muu1**2 + 2
     &*muu3**2 + muu4**2 - tredDip) + g2Rt*g3Lcu*g4Rcu*muu3*muugl*(-1 + muu2**2 + muu3**2 + 2*muu4
     &**2 - tredDip) + g2Lt*g3Rcu*g4Lcu*((-1 + muu2**2 + muu4**2 - tredDip)*tredDip + muu1**2*(-muu2
     &**2 + tredDip) + muu3**2*(muu4**2 + tredDip))) + g1Lt*(g2Rt*g3Rcu*g4Lcu*muu3*muu4*(-muu1**2 
     &- muu2**2 + muu3**2 + muu4**2) + g2Lt*g3Lcu*g4Rcu*muu4*muugl*(-1 + muu1**2 + 2*muu3**2 + muu4
     &**2 - tredDip) + g2Lt*g3Rcu*g4Lcu*muu3*muugl*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip) + g2Rt
     &*g3Lcu*g4Rcu*((-1 + muu2**2 + muu4**2 - tredDip)*tredDip + muu1**2*(-muu2**2 + tredDip) 
     &+ muu3**2*(muu4**2 + tredDip))) + (g1Rt*(2*g2Rt*g3Lcu*g4Lcu*muu3*muu4*muugl + g2Rt*g3Rcu*g4Rcu
     &*(-1 + muu3**2 + muu4**2)*muugl + g2Lt*g3Rcu*g4Rcu*muu4*(-muu1**2 + muu3**2 + tredDip) + g2Lt
     &*g3Lcu*g4Lcu*muu3*(-muu2**2 + muu4**2 + tredDip)) + g1Lt*(2*g2Lt*g3Rcu*g4Rcu*muu3*muu4*muugl 
     &+ g2Lt*g3Lcu*g4Lcu*(-1 + muu3**2 + muu4**2)*muugl + g2Rt*g3Lcu*g4Lcu*muu4*(-muu1**2 + muu3**2 
     &+ tredDip) + g2Rt*g3Rcu*g4Rcu*muu3*(-muu2**2 + muu4**2 + tredDip)))*muupu(n)))/(DenGlT*dconjg(DenGauU(n))))

      enddo


!************************************************   2 Re(M^{gluino}_u M^{gaugino ast}_t)   ************************************************

      do n=1,6 

      call stst2QQ_setGauginoCoupling(1,n,g3Lt,g3Rt,g3Lct,g3Rct,g4Lt,g4Rt,g4Lct,g4Rct,g3Lu,g3Ru,g3Lcu,g3Rcu,g4Lu,g4Ru,g4Lcu,g4Rcu)

      stst2qq_M2TreeDipoleT1T2=stst2qq_M2TreeDipoleT1T2+dble((CF*sqrtSDip**4*(g1Ru*(g2Lu*g3Lct*g4Rct*muu3
     &*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2) + g2Ru*g3Rct*g4Lct*muu4*muugl*(-muu1**2 + muu3
     &**2 + tredDip) + g2Ru*g3Lct*g4Rct*muu3*muugl*(-muu2**2 + muu4**2 + tredDip) + g2Lu*g3Rct*g4Lct
     &*((-1 + muu2**2 + muu4**2 - tredDip)*tredDip + muu1**2*(-muu2**2 + tredDip) + muu3**2*(muu4*
     &*2 + tredDip))) + g1Lu*(g2Ru*g3Rct*g4Lct*muu3*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2) + g2Lu
     &*g3Lct*g4Rct*muu4*muugl*(-muu1**2 + muu3**2 + tredDip) + g2Lu*g3Rct*g4Lct*muu3*muugl*(-muu2
     &**2 + muu4**2 + tredDip) + g2Ru*g3Lct*g4Rct*((-1 + muu2**2 + muu4**2 - tredDip)*tredDip + muu1
     &**2*(-muu2**2 + tredDip) + muu3**2*(muu4**2 + tredDip))) + (g1Ru*(2*g2Ru*g3Lct*g4Lct*muu3
     &*muu4*muugl + g2Ru*g3Rct*g4Rct*(-1 + muu3**2 + muu4**2)*muugl + g2Lu*g3Rct*g4Rct*muu4*(-1 + muu1
     &**2 + 2*muu3**2 + muu4**2 - tredDip) + g2Lu*g3Lct*g4Lct*muu3*(-1 + muu2**2 + muu3**2 + 2*muu4
     &**2 - tredDip)) + g1Lu*(2*g2Lu*g3Rct*g4Rct*muu3*muu4*muugl + g2Lu*g3Lct*g4Lct*(-1 + muu3**2 
     &+ muu4**2)*muugl + g2Ru*g3Lct*g4Lct*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + g2Ru
     &*g3Rct*g4Rct*muu3*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)))*muupt(n)))/(DenGlU*dconjg(DenGauT(n))))

      enddo

      end

      double precision function stst2qq_M2TreeDipole(sdip,tdip,udip)
      implicit none

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"
#include "stst2QQ_GenCouplings.h"

      integer i,j,k,l,m,n
      double precision sdip,tdip,udip,sqrtSdip,tredDip
      double complex g1L, g1Lc, g2L, g2Lc, g3L, g3Lc, g4L, g4Lc
      double complex g1R, g1Rc, g2R, g2Rc, g3R, g3Rc, g4R, g4Rc
      double complex test, GltGlt,GluGlu,GltGlu,NeutNeut,NeuuNeuu,NeutNeuu,GltNeut,GltNeuu,GluNeut,GluNeuu
      
      test = 0d0
      GltGlt = 0d0
      GluGlu = 0d0
      GltGlu = 0d0
      NeutNeut = 0d0
      NeuuNeuu = 0d0
      NeutNeuu = 0d0
      GltNeut = 0d0
      GltNeuu = 0d0
      GluNeut = 0d0
      GluNeuu = 0d0

c Dipole Mandelstam variables
      sqrtSdip = dsqrt(sdip)
      tredDip = tdip/sdip
    
c 	Propagator--------------------------------------------------------------------------------------

      do i=1,4 !hier gluino und neutralino unterscheiden
          DenGauT(i)=(tdip - MNeu(i)**2)! + cI*MNeu(i)*wNeu(i)
      enddo
      DenGauT(5)=(tdip - MCha(1)**2)! + cI*MCha(1)*wCha(1)
      DenGauT(6)=(tdip - MCha(2)**2)! + cI*MCha(2)*wCha(2)
    
      do i=1,4 !hier gluino und neutralino unterscheiden
          DenGauU(i)=(udip - MNeu(i)**2)! + cI*MNeu(i)*wNeu(i)
      enddo
      DenGauU(5)=(udip - MCha(1)**2)! + cI*MCha(1)*wCha(1)
      DenGauU(6)=(udip - MCha(2)**2)! + cI*MCha(2)*wCha(2)


      DenGlT = (tdip-MGl**2)
      DenGlU = (udip-MGl**2)
      
! 	  print*,'prop',DenGlT,DenGlU

c-------------------------------------------------------------------------------------------------------
c------------------------------------t-channel Gluino squared-------------------------------------------
c-------------------------------------------------------------------------------------------------------

c 	  colfac = Nc/2d0*CF
    
      g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
      g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
      g1Lc = dconjg(g1L)
      g1Rc = dconjg(g1R)

      g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
      g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
      g2Lc = dconjg(g2L)
      g2Rc = dconjg(g2R)

      GltGlt = (-((sqrtSDip**4*(g1R*
     -            (g1Rc*g2R*muugl*(g2Rc*(-1 + muu3**2 + muu4**2)* 
     -              muugl + g2Lc*muu4*(-muu1**2 + muu3**2 + tredDip)) + 
     -           g1Lc*g2R*muu3*muugl*
     -            (2*g2Lc*muu4*muugl + 
     -              g2Rc*(-muu2**2 + muu4**2 + tredDip)) + 
     -           g1Lc*g2L*muu3*
     -            (2*g2Rc*muu4*tredDip + 
     -              g2Lc*muugl*(-muu2**2 + muu4**2 + tredDip)) + 
     -           g1Rc*g2L*(g2Rc*muu4*muugl*
     -               (-muu1**2 + muu3**2 + tredDip) + 
     -              g2Lc*(muu3**2*muu4**2 + 
     -                 muu1**2*(muu2**2 - muu4**2 - tredDip) + tredDip + 
     -                 tredDip**2 - muu2**2*(muu3**2 + tredDip)))) + 
     -        g1L*(g1Lc*g2L*muugl*
     -            (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -              g2Rc*muu4*(-muu1**2 + muu3**2 + tredDip)) + 
     -           g1Rc*g2L*muu3*muugl*
     -            (2*g2Rc*muu4*muugl + 
     -              g2Lc*(-muu2**2 + muu4**2 + tredDip)) + 
     -           g1Rc*g2R*muu3*
     -            (2*g2Lc*muu4*tredDip + 
     -              g2Rc*muugl*(-muu2**2 + muu4**2 + tredDip)) + 
     -           g1Lc*g2R*(g2Lc*muu4*muugl*
     -               (-muu1**2 + muu3**2 + tredDip) + 
     -              g2Rc*(muu3**2*muu4**2 + 
     -                 muu1**2*(muu2**2 - muu4**2 - tredDip) + tredDip + 
     -                 tredDip**2 - muu2**2*(muu3**2 + tredDip))))))/
     -    (DenGlT*dconjg(DenGlT))))
   
      
c-------------------------------------------------------------------------------------------------------
c------------------------------------u-channel Gluino squared-------------------------------------------
c-------------------------------------------------------------------------------------------------------
c 	  colfac = Nc/2d0*CF
  
    
      g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
      g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)
      g1Lc = dconjg(g1L)
      g1Rc = dconjg(g1R)

      g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
      g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)
      g2Lc = dconjg(g2L)
      g2Rc = dconjg(g2R)


      GluGlu = (-((sqrtSDip**4*(g1R*(g1Rc*g2R*muugl*
     -            (g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -              g2Lc*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip)) + 
     -           g1Lc*g2R*muu3*muugl*
     -            (2*g2Lc*muu4*muugl + 
     -              g2Rc*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) +
     -             g1Lc*g2L*muu3*
     -            (2*g2Rc*muu4*
     -               (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                 tredDip) + 
     -              g2Lc*muugl*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + 
     -           g1Rc*g2L*(g2Rc*muu4*muugl*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + 
     -              g2Lc*(muu3**4 - muu4**2 + muu4**4 + 
     -                 muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tredDip) + 
     -                 muu1**2*(muu2**2 + muu4**2 - tredDip) + tredDip - 
     -                 muu2**2*tredDip - 2*muu4**2*tredDip + tredDip**2))) + 
     -        g1L*(g1Lc*g2L*muugl*
     -            (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -              g2Rc*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip)) + 
     -           g1Rc*g2Rc*muu3*muugl*
     -            (2*g2L*muu4*muugl + 
     -              g2R*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) +
     -             g1Rc*g2Lc*muu3*
     -            (2*g2R*muu4*
     -               (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - 
     -                 tredDip) + 
     -              g2L*muugl*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + 
     -           g1Lc*g2R*(g2Lc*muu4*muugl*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + 
     -              g2Rc*(muu3**4 - muu4**2 + muu4**4 + 
     -                 muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tredDip) + 
     -                 muu1**2*(muu2**2 + muu4**2 - tredDip) + tredDip - 
     -                 muu2**2*tredDip - 2*muu4**2*tredDip + tredDip**2)))))/
     -    (DenGlU*dconjg(DenGlU))))

c-------------------------------------------------------------------------------------------------------
c------------------------------------t-und-u-channel Gluino Interferenc-------------------------------------------
c------------------------------------------------------------------------------------------------------- 
c 	  colfac = -1d0/2d0*CF
    
      g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
      g1Lc = dconjg(GluLt1(isf2,itt2,ftt1,igen2,fgen1))
      g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)
      g1Rc = dconjg(GluRt1(isf2,itt2,ftt1,igen2,fgen1))

      g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
      g2Lc = dconjg(GluLt1(isf1,itt1,ftt2,igen1,fgen2))
      g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)
      g2Rc = dconjg(GluRt1(isf1,itt1,ftt2,igen1,fgen2))

      GltGlu = dreal((-2*sqrtSDip**4*(g1R*(g1Rc*g2R*muugl*
     -          (g2Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2Lc*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -               tredDip)) + 
     -         g1Lc*g2R*muu3*muugl*
     -          (2*g2Lc*muu4*muugl + 
     -            g2Rc*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + 
     -         g1Lc*g2L*muu3*
     -          (g2Rc*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -            g2Lc*muugl*(-muu2**2 + muu4**2 + tredDip)) + 
     -         g1Rc*g2L*(g2Rc*muu4*muugl*
     -             (-muu1**2 + muu3**2 + tredDip) + 
     -            g2Lc*((-1 + muu2**2 + muu4**2 - tredDip)*tredDip + 
     -               muu1**2*(-muu2**2 + tredDip) + 
     -               muu3**2*(muu4**2 + tredDip)))) + 
     -      g1L*(g1Lc*g2L*muugl*
     -          (g2Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2Rc*muu4*(-1 + muu1**2 + 2*muu3**2 + muu4**2 - 
     -               tredDip)) + 
     -         g1Rc*g2Lc*muu3*
     -          (g2R*muu4*(-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -            g2L*muugl*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - 
     -               tredDip)) + 
     -         g1Rc*g2Rc*muu3*muugl*
     -          (2*g2L*muu4*muugl + g2R*(-muu2**2 + muu4**2 + tredDip)) +
     -          g1Lc*g2R*(g2Lc*muu4*muugl*
     -             (-muu1**2 + muu3**2 + tredDip) + 
     -            g2Rc*((-1 + muu2**2 + muu4**2 - tredDip)*tredDip + 
     -               muu1**2*(-muu2**2 + tredDip) + 
     -               muu3**2*(muu4**2 + tredDip))))))/
     -  (DenGlT*dconjg(DenGlU)))

c-------------------------------------------------------------------------------
c-------------------------------Neutralino squared------------------------------
c-------------------------------------------------------------------------------

c-------------------------------------------------------------------------------------------------------
c--------------------------------t-channel Neutralino squared-------------------------------------------
c-------------------------------------------------------------------------------------------------------      
      do n=1,6
      do k=1,6

c 	  colfac = Nc**2

      g3L = ALt1(k,isf1,itt1,ftt1,igen1,fgen1)
      g3R = ARt1(k,isf1,itt1,ftt1,igen1,fgen1)
      g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
      g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))

      g4L = ALt1(k,isf2,itt2,ftt2,igen2,fgen2)
      g4R = ARt1(k,isf2,itt2,ftt2,igen2,fgen2)
      g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
      g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))


      NeutNeut = NeutNeut + ( -((sqrtSDip**4*
     -        (g3R*g3Rc*g4L*g4Lc*muu1**2*muu2**2 + 
     -        g3L*g3Lc*g4R*g4Rc*muu1**2*muu2**2 - 
     -        g3R*g3Rc*g4L*g4Lc*muu2**2*muu3**2 - 
     -        g3L*g3Lc*g4R*g4Rc*muu2**2*muu3**2 - 
     -        g3R*g3Rc*g4L*g4Lc*muu1**2*muu4**2 - 
     -        g3L*g3Lc*g4R*g4Rc*muu1**2*muu4**2 + 
     -        g3R*g3Rc*g4L*g4Lc*muu3**2*muu4**2 + 
     -        g3L*g3Lc*g4R*g4Rc*muu3**2*muu4**2 + 
     -        g3R*g3Rc*g4L*g4Lc*tredDip + g3L*g3Lc*g4R*g4Rc*tredDip - 
     -        g3R*g3Rc*g4L*g4Lc*muu1**2*tredDip - 
     -        g3L*g3Lc*g4R*g4Rc*muu1**2*tredDip - 
     -        g3R*g3Rc*g4L*g4Lc*muu2**2*tredDip - 
     -        g3L*g3Lc*g4R*g4Rc*muu2**2*tredDip + 
     -        2*g3L*g3Rc*g4Lc*g4R*muu3*muu4*tredDip + 
     -        2*g3Lc*g3R*g4L*g4Rc*muu3*muu4*tredDip + 
     -        g3R*g3Rc*g4L*g4Lc*tredDip**2 + g3L*g3Lc*g4R*g4Rc*tredDip**2 + 
     -        (g3L*g4L*(g3Lc*g4Rc*muu4*(-muu1**2 + muu3**2 + tredDip) + 
     -              g3Rc*g4Lc*muu3*(-muu2**2 + muu4**2 + tredDip)) + 
     -           g3R*g4R*(g3Rc*g4Lc*muu4*
     -               (-muu1**2 + muu3**2 + tredDip) + 
     -              g3Lc*g4Rc*muu3*(-muu2**2 + muu4**2 + tredDip)))*
     -         muuqt(k) + muupt(n)*
     -         (g3Lc*g4Lc*(g3L*g4R*muu4*(-muu1**2 + muu3**2 + tredDip) + 
     -              g3R*g4L*muu3*(-muu2**2 + muu4**2 + tredDip)) + 
     -           g3Rc*g4Rc*(g3R*g4L*muu4*
     -               (-muu1**2 + muu3**2 + tredDip) + 
     -              g3L*g4R*muu3*(-muu2**2 + muu4**2 + tredDip)) + 
     -           (g3L*g4L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -                 g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)) + 
     -              g3R*g4R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -                 g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)))*muuqt(k)))
     -      )/(dconjg(DenGauT(n))*DenGauT(k))))

   
c-------------------------------------------------------------------------------------------------------
c--------------------------------u-channel Neutralino squared-------------------------------------------
c-------------------------------------------------------------------------------------------------------

c 	  colfac = Nc**2

      g3L = ALt1(k,isf2,itt2,ftt1,igen2,fgen1)
      g3R = ARt1(k,isf2,itt2,ftt1,igen2,fgen1)
      g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
      g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))

      g4L = ALt1(k,isf1,itt1,ftt2,igen1,fgen2)
      g4R = ARt1(k,isf1,itt1,ftt2,igen1,fgen2)
      g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
      g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))

      NeuuNeuu = NeuuNeuu + (-((sqrtSDip**4*(g3R*g3Rc*
     -        g4L*g4Lc*muu1**2*muu2**2 + 
     -        g3L*g3Lc*g4R*g4Rc*muu1**2*muu2**2 - 
     -        g3R*g3Rc*g4L*g4Lc*muu3**2 - g3L*g3Lc*g4R*g4Rc*muu3**2 + 
     -        g3R*g3Rc*g4L*g4Lc*muu2**2*muu3**2 + 
     -        g3L*g3Lc*g4R*g4Rc*muu2**2*muu3**2 + 
     -        g3R*g3Rc*g4L*g4Lc*muu3**4 + g3L*g3Lc*g4R*g4Rc*muu3**4 - 
     -        2*g3L*g3Rc*g4Lc*g4R*muu3*muu4 - 
     -        2*g3Lc*g3R*g4L*g4Rc*muu3*muu4 + 
     -        2*g3L*g3Rc*g4Lc*g4R*muu1**2*muu3*muu4 + 
     -        2*g3Lc*g3R*g4L*g4Rc*muu1**2*muu3*muu4 + 
     -        2*g3L*g3Rc*g4Lc*g4R*muu2**2*muu3*muu4 + 
     -        2*g3Lc*g3R*g4L*g4Rc*muu2**2*muu3*muu4 + 
     -        2*g3L*g3Rc*g4Lc*g4R*muu3**3*muu4 + 
     -        2*g3Lc*g3R*g4L*g4Rc*muu3**3*muu4 - 
     -        g3R*g3Rc*g4L*g4Lc*muu4**2 - g3L*g3Lc*g4R*g4Rc*muu4**2 + 
     -        g3R*g3Rc*g4L*g4Lc*muu1**2*muu4**2 + 
     -        g3L*g3Lc*g4R*g4Rc*muu1**2*muu4**2 + 
     -        3*g3R*g3Rc*g4L*g4Lc*muu3**2*muu4**2 + 
     -        3*g3L*g3Lc*g4R*g4Rc*muu3**2*muu4**2 + 
     -        2*g3L*g3Rc*g4Lc*g4R*muu3*muu4**3 + 
     -        2*g3Lc*g3R*g4L*g4Rc*muu3*muu4**3 + 
     -        g3R*g3Rc*g4L*g4Lc*muu4**4 + g3L*g3Lc*g4R*g4Rc*muu4**4 + 
     -        g3R*g3Rc*g4L*g4Lc*tredDip + g3L*g3Lc*g4R*g4Rc*tredDip - 
     -        g3R*g3Rc*g4L*g4Lc*muu1**2*tredDip - 
     -        g3L*g3Lc*g4R*g4Rc*muu1**2*tredDip - 
     -        g3R*g3Rc*g4L*g4Lc*muu2**2*tredDip - 
     -        g3L*g3Lc*g4R*g4Rc*muu2**2*tredDip - 
     -        2*g3R*g3Rc*g4L*g4Lc*muu3**2*tredDip - 
     -        2*g3L*g3Lc*g4R*g4Rc*muu3**2*tredDip - 
     -        2*g3L*g3Rc*g4Lc*g4R*muu3*muu4*tredDip - 
     -        2*g3Lc*g3R*g4L*g4Rc*muu3*muu4*tredDip - 
     -        2*g3R*g3Rc*g4L*g4Lc*muu4**2*tredDip - 
     -        2*g3L*g3Lc*g4R*g4Rc*muu4**2*tredDip + 
     -        g3R*g3Rc*g4L*g4Lc*tredDip**2 + g3L*g3Lc*g4R*g4Rc*tredDip**2 + 
     -        (g3L*g4L*(g3Lc*g4Rc*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + 
     -              g3Rc*g4Lc*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + 
     -           g3R*g4R*(g3Rc*g4Lc*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + 
     -              g3Lc*g4Rc*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)))*
     -         muuqu(k) + muupu(n)*
     -         (g3Lc*g4Lc*(g3L*g4R*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + 
     -              g3R*g4L*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + 
     -           g3Rc*g4Rc*(g3R*g4L*muu4*
     -               (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + 
     -              g3L*g4R*muu3*
     -               (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + 
     -           (g3L*g4L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -                 g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)) + 
     -              g3R*g4R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -                 g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)))*muuqu(k)))
     -      )/(dconjg(DenGauU(n))*DenGauU(k))))   

c-------------------------------------------------------------------------------------------------------
c-----------------------------t-und-u-channel Neutralino Inteference------------------------------------
c-------------------------------------------------------------------------------------------------------
c	  colfac = Nc
    
      g3L = ALt1(k,isf1,itt1,ftt1,igen1,fgen1)
      g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))

      g3R = ARt1(k,isf1,itt1,ftt1,igen1,fgen1)
      g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))

      g4L = ALt1(k,isf2,itt2,ftt2,igen2,fgen2)
      g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))

      g4R = ARt1(k,isf2,itt2,ftt2,igen2,fgen2)
      g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))

      NeutNeuu = NeutNeuu + dreal((-2*sqrtSDip**4*
     -      (-(g3R*g3Rc*g4L*g4Lc*muu1**2*muu2**2) - 
     -      g3L*g3Lc*g4R*g4Rc*muu1**2*muu2**2 - 
     -      g3L*g3Rc*g4Lc*g4R*muu1**2*muu3*muu4 - 
     -      g3Lc*g3R*g4L*g4Rc*muu1**2*muu3*muu4 - 
     -      g3L*g3Rc*g4Lc*g4R*muu2**2*muu3*muu4 - 
     -      g3Lc*g3R*g4L*g4Rc*muu2**2*muu3*muu4 + 
     -      g3L*g3Rc*g4Lc*g4R*muu3**3*muu4 + 
     -      g3Lc*g3R*g4L*g4Rc*muu3**3*muu4 + 
     -      g3R*g3Rc*g4L*g4Lc*muu3**2*muu4**2 + 
     -      g3L*g3Lc*g4R*g4Rc*muu3**2*muu4**2 + 
     -      g3L*g3Rc*g4Lc*g4R*muu3*muu4**3 + 
     -      g3Lc*g3R*g4L*g4Rc*muu3*muu4**3 - g3R*g3Rc*g4L*g4Lc*tredDip - 
     -      g3L*g3Lc*g4R*g4Rc*tredDip + g3R*g3Rc*g4L*g4Lc*muu1**2*tredDip + 
     -      g3L*g3Lc*g4R*g4Rc*muu1**2*tredDip + 
     -      g3R*g3Rc*g4L*g4Lc*muu2**2*tredDip + 
     -      g3L*g3Lc*g4R*g4Rc*muu2**2*tredDip + 
     -      g3R*g3Rc*g4L*g4Lc*muu3**2*tredDip + 
     -      g3L*g3Lc*g4R*g4Rc*muu3**2*tredDip + 
     -      g3R*g3Rc*g4L*g4Lc*muu4**2*tredDip + 
     -      g3L*g3Lc*g4R*g4Rc*muu4**2*tredDip - 
     -      g3R*g3Rc*g4L*g4Lc*tredDip**2 - g3L*g3Lc*g4R*g4Rc*tredDip**2 + 
     -      (g3L*g4L*(g3Lc*g4Rc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + 
     -            g3Rc*g4Lc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + 
     -         g3R*g4R*(g3Rc*g4Lc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + 
     -            g3Lc*g4Rc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)))*
     -       muuqt(k) + muupu(n)*
     -       (g3Lc*g4Lc*(g3L*g4R*muu4*(-muu1**2 + muu3**2 + tredDip) + 
     -            g3R*g4L*muu3*(-muu2**2 + muu4**2 + tredDip)) + 
     -         g3Rc*g4Rc*(g3R*g4L*muu4*(-muu1**2 + muu3**2 + tredDip) + 
     -            g3L*g4R*muu3*(-muu2**2 + muu4**2 + tredDip)) + 
     -         (g3L*g4L*(2*g3Rc*g4Rc*muu3*muu4 + 
     -               g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)) + 
     -            g3R*g4R*(2*g3Lc*g4Lc*muu3*muu4 + 
     -               g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)))*muuqt(k))))/
     -  (dconjg(DenGauU(n))*DenGauT(k)))
  
      enddo

c------------------------------------------------------------------------------------------------------
c----------------------------------------INTEFERENZ-PARTS---------------------------------------------- 
c------------------------------------------------------------------------------------------------------     

c-------------------------------------------------------------------------------------------------------
c----------------------------t-t-channel Gluino-Neutralino Inteference----------------------------------
c-------------------------------------------------------------------------------------------------------
c     colfac = 0d0
      
      g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
      g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)

      g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
      g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)

      g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
      g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))

      g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
      g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))

      GltNeut = GltNeut + dreal((-2*sqrtSDip**4*(g1R*(2*g2L*g3Lc*g4Rc*muu3*muu4*tredDip + 
     -         g2R*g3Rc*g4Lc*muu4*muugl*(-muu1**2 + muu3**2 + tredDip) + 
     -         g2R*g3Lc*g4Rc*muu3*muugl*(-muu2**2 + muu4**2 + tredDip) + 
     -         g2L*g3Rc*g4Lc*
     -          (muu3**2*muu4**2 + 
     -            muu1**2*(muu2**2 - muu4**2 - tredDip) + tredDip + 
     -            tredDip**2 - muu2**2*(muu3**2 + tredDip))) + 
     -      g1L*(2*g2R*g3Rc*g4Lc*muu3*muu4*tredDip + 
     -         g2L*g3Lc*g4Rc*muu4*muugl*(-muu1**2 + muu3**2 + tredDip) + 
     -         g2L*g3Rc*g4Lc*muu3*muugl*(-muu2**2 + muu4**2 + tredDip) + 
     -         g2R*g3Lc*g4Rc*
     -          (muu3**2*muu4**2 + 
     -            muu1**2*(muu2**2 - muu4**2 - tredDip) + tredDip + 
     -            tredDip**2 - muu2**2*(muu3**2 + tredDip))) + 
     -      (g1R*(2*g2R*g3Lc*g4Lc*muu3*muu4*muugl + 
     -            g2R*g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2L*g3Rc*g4Rc*muu4*(-muu1**2 + muu3**2 + tredDip) + 
     -            g2L*g3Lc*g4Lc*muu3*(-muu2**2 + muu4**2 + tredDip)) + 
     -         g1L*(2*g2L*g3Rc*g4Rc*muu3*muu4*muugl + 
     -            g2L*g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2R*g3Lc*g4Lc*muu4*(-muu1**2 + muu3**2 + tredDip) + 
     -            g2R*g3Rc*g4Rc*muu3*(-muu2**2 + muu4**2 + tredDip)))*
     -       muupt(n)))/(DenGlT*dconjg(DenGauT(n))))        
   
  
c-------------------------------------------------------------------------------------------------------
c----------------------------u-u-channel Gluino-Neutralino Inteference----------------------------------
c-------------------------------------------------------------------------------------------------------
c 	  colfac = 0d0
      
      g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
      g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)

      g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
      g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)

      g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
      g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))

      g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
      g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))

      GluNeuu = GluNeuu + dreal((-2*sqrtSDip**4*(g1R*(2*g2L*g3Lc*g4Rc*muu3*muu4*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tredDip) + 
     -         g2R*g3Rc*g4Lc*muu4*muugl*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + 
     -         g2R*g3Lc*g4Rc*muu3*muugl*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip) + 
     -         g2L*g3Rc*g4Lc*
     -          (muu3**4 - muu4**2 + muu4**4 + 
     -            muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tredDip) + 
     -            muu1**2*(muu2**2 + muu4**2 - tredDip) + tredDip - 
     -            muu2**2*tredDip - 2*muu4**2*tredDip + tredDip**2)) + 
     -      g1L*(2*g2R*g3Rc*g4Lc*muu3*muu4*
     -          (-1 + muu1**2 + muu2**2 + muu3**2 + muu4**2 - tredDip) + 
     -         g2L*g3Lc*g4Rc*muu4*muugl*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + 
     -         g2L*g3Rc*g4Lc*muu3*muugl*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip) + 
     -         g2R*g3Lc*g4Rc*
     -          (muu3**4 - muu4**2 + muu4**4 + 
     -            muu3**2*(-1 + muu2**2 + 3*muu4**2 - 2*tredDip) + 
     -            muu1**2*(muu2**2 + muu4**2 - tredDip) + tredDip - 
     -            muu2**2*tredDip - 2*muu4**2*tredDip + tredDip**2)) + 
     -      (g1R*(2*g2R*g3Lc*g4Lc*muu3*muu4*muugl + 
     -            g2R*g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2L*g3Rc*g4Rc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + 
     -            g2L*g3Lc*g4Lc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + 
     -         g1L*(2*g2L*g3Rc*g4Rc*muu3*muu4*muugl + 
     -            g2L*g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2R*g3Lc*g4Lc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + 
     -            g2R*g3Rc*g4Rc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)))*
     -       muupu(n)))/(DenGlU*dconjg(DenGauU(n))))        
      
!         
c-------------------------------------------------------------------------------------------------------
c----------------------------t-u-channel Gluino-Neutralino Inteference----------------------------------
c-------------------------------------------------------------------------------------------------------
c     colfac = CF*Nc
      
      g1L = GluLt1(isf1,itt1,ftt1,igen1,fgen1)
      g1R = GluRt1(isf1,itt1,ftt1,igen1,fgen1)

      g2L = GluLt1(isf2,itt2,ftt2,igen2,fgen2)
      g2R = GluRt1(isf2,itt2,ftt2,igen2,fgen2)

      g3Lc = dconjg(ALt1(n,isf2,itt2,ftt1,igen2,fgen1))
      g3Rc = dconjg(ARt1(n,isf2,itt2,ftt1,igen2,fgen1))

      g4Lc = dconjg(ALt1(n,isf1,itt1,ftt2,igen1,fgen2))
      g4Rc = dconjg(ARt1(n,isf1,itt1,ftt2,igen1,fgen2))
    
    
      GltNeuu = GltNeuu + dreal((-2*sqrtSDip**4*(g1R*(g2L*g3Lc*g4Rc*muu3*muu4*
     -          (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -         g2R*g3Rc*g4Lc*muu4*muugl*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + 
     -         g2R*g3Lc*g4Rc*muu3*muugl*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip) + 
     -         g2L*g3Rc*g4Lc*
     -          ((-1 + muu2**2 + muu4**2 - tredDip)*tredDip + 
     -            muu1**2*(-muu2**2 + tredDip) + muu3**2*(muu4**2 + tredDip)
     -            )) + g1L*(g2R*g3Rc*g4Lc*muu3*muu4*
     -          (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -         g2L*g3Lc*g4Rc*muu4*muugl*
     -          (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + 
     -         g2L*g3Rc*g4Lc*muu3*muugl*
     -          (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip) + 
     -         g2R*g3Lc*g4Rc*
     -          ((-1 + muu2**2 + muu4**2 - tredDip)*tredDip + 
     -            muu1**2*(-muu2**2 + tredDip) + muu3**2*(muu4**2 + tredDip)
     -            )) + (g1R*(2*g2R*g3Lc*g4Lc*muu3*muu4*muugl + 
     -            g2R*g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2L*g3Rc*g4Rc*muu4*(-muu1**2 + muu3**2 + tredDip) + 
     -            g2L*g3Lc*g4Lc*muu3*(-muu2**2 + muu4**2 + tredDip)) + 
     -         g1L*(2*g2L*g3Rc*g4Rc*muu3*muu4*muugl + 
     -            g2L*g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2R*g3Lc*g4Lc*muu4*(-muu1**2 + muu3**2 + tredDip) + 
     -            g2R*g3Rc*g4Rc*muu3*(-muu2**2 + muu4**2 + tredDip)))*
     -       muupu(n)))/(DenGlT*dconjg(DenGauU(n))))   
   
  
!         
c-------------------------------------------------------------------------------------------------------
c----------------------------u-t-channel Gluino-Neutralino Inteference----------------------------------
c-------------------------------------------------------------------------------------------------------
c     colfac = CF*Nc
      
      g1L = GluLt1(isf2,itt2,ftt1,igen2,fgen1)
      g1R = GluRt1(isf2,itt2,ftt1,igen2,fgen1)

      g2L = GluLt1(isf1,itt1,ftt2,igen1,fgen2)
      g2R = GluRt1(isf1,itt1,ftt2,igen1,fgen2)

      g3Lc = dconjg(ALt1(n,isf1,itt1,ftt1,igen1,fgen1))
      g3Rc = dconjg(ARt1(n,isf1,itt1,ftt1,igen1,fgen1))

      g4Lc = dconjg(ALt1(n,isf2,itt2,ftt2,igen2,fgen2))
      g4Rc = dconjg(ARt1(n,isf2,itt2,ftt2,igen2,fgen2))

      GluNeut = GluNeut + dreal((-2*sqrtSDip**4*(g1R*(g2L*g3Lc*g4Rc*muu3*muu4*
     -          (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -         g2R*g3Rc*g4Lc*muu4*muugl*(-muu1**2 + muu3**2 + tredDip) + 
     -         g2R*g3Lc*g4Rc*muu3*muugl*(-muu2**2 + muu4**2 + tredDip) + 
     -         g2L*g3Rc*g4Lc*
     -          ((-1 + muu2**2 + muu4**2 - tredDip)*tredDip + 
     -            muu1**2*(-muu2**2 + tredDip) + muu3**2*(muu4**2 + tredDip)
     -            )) + g1L*(g2R*g3Rc*g4Lc*muu3*muu4*
     -          (-muu1**2 - muu2**2 + muu3**2 + muu4**2) + 
     -         g2L*g3Lc*g4Rc*muu4*muugl*(-muu1**2 + muu3**2 + tredDip) + 
     -         g2L*g3Rc*g4Lc*muu3*muugl*(-muu2**2 + muu4**2 + tredDip) + 
     -         g2R*g3Lc*g4Rc*
     -          ((-1 + muu2**2 + muu4**2 - tredDip)*tredDip + 
     -            muu1**2*(-muu2**2 + tredDip) + muu3**2*(muu4**2 + tredDip)
     -            )) + (g1R*(2*g2R*g3Lc*g4Lc*muu3*muu4*muugl + 
     -            g2R*g3Rc*g4Rc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2L*g3Rc*g4Rc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + 
     -            g2L*g3Lc*g4Lc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)) + 
     -         g1L*(2*g2L*g3Rc*g4Rc*muu3*muu4*muugl + 
     -            g2L*g3Lc*g4Lc*(-1 + muu3**2 + muu4**2)*muugl + 
     -            g2R*g3Lc*g4Lc*muu4*
     -             (-1 + muu1**2 + 2*muu3**2 + muu4**2 - tredDip) + 
     -            g2R*g3Rc*g4Rc*muu3*
     -             (-1 + muu2**2 + muu3**2 + 2*muu4**2 - tredDip)))*
     -       muupt(n)))/(DenGlU*dconjg(DenGauT(n))))


      enddo
      
c---------------------Calculation of complete Tree-level-------------------------------------------        
! 
      stst2qq_M2TreeDipole = (Nc*CF/2d0)*(GltGlt+GluGlu)-CF/2d0*GltGlu+
     -      NC**2*(NeutNeut+NeuuNeuu) + Nc*NeutNeuu+
     -      0d0*(GltNeut+GluNeuu) + CF*Nc*(GltNeuu + GluNeut)


      end


