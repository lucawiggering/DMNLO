      double precision function NeuQ2qx_M2RadiationDipoleG()
! p1 (quark), p2 (gluon), p3 (gluon)

      implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

      double precision dipcolor,tree,NeuQ2qx_M2TreeDipole,muQ2,zi,zj,Q2,xdip
      double precision sdip,tdip,udip,Rdip,kappa,yijk,vijk,viji,zP,zM,vijkTilde
      double precision NeuQ2qx_M2_FEFS_gg,NeuQ2qx_M2_FEIS_gg

! ********************************************************** FE-FS  (from hep-ph/0201036v1 Section 5.1.1+5.1.2 ***********************************************************

      Q2 = sman
      muQ2 = m3**2/Q2

! -----------------------------------------------    Dipoles with quark as emitter     ----------------------------------------------------------------------------------

! (1.1) emission off quark (q->q(p1)+g(p3)) with gluon p2 as spectator

      sdip = Q2
      tdip = m2**2 - 2*(((-m3**2 + Q2)*(papb + m2**2))/(2*Q2) + (kappa(Q2,m3**2,0d0)*(pbp2 - ((pap2 + pbp2)*(papb + m2**2))/Q2))/kappa(Q2,2*p1p3 + m3**2,0d0))
      udip = m1**2 - 2*(((-m3**2 + Q2)*(papb + m1**2))/(2*Q2) + (kappa(Q2,m3**2,0d0)*(pap2 - ((pap2 + pbp2)*(papb + m1**2))/Q2))/kappa(Q2,2*p1p3 + m3**2,0d0))

      tree = NeuQ2qx_M2TreeDipole(sdip,tdip,udip)

      yijk = p1p3/(p1p2 + p1p3 + p2p3)
      vijk = 1d0
      viji = 1d0
      zP = ((1-muQ2)*yijk)/(muQ2+(1-muQ2)*yijk)*(1+viji*vijk)/2d0
      zM = ((1-muQ2)*yijk)/(muQ2+(1-muQ2)*yijk)*(1-viji*vijk)/2d0
      vijkTilde = kappa(1d0,muQ2,0d0)/(1-muQ2)

      zj = p1p2/(p1p2 + p2p3)

      dipcolor = -CV/(2*CF)
      NeuQ2qx_M2RadiationDipoleG =  0.5d0*(CF*dipcolor*gs**2*tree*(-((vijkTilde*(1 + m3**2/p1p3 + zj))/vijk) + 2/(1 - (1 - yijk)*zj)))/p1p3


! (1.2) emission off quark (q->q(p1)+g(p2)) with gluon p3 as spectator

      sdip = Q2
      tdip = m2**2 - 2*(((-m3**2 + Q2)*(papb + m2**2))/(2*Q2) + (kappa(Q2,m3**2,0d0)*(pbp3 - ((pap3 + pbp3)*(papb + m2**2))/Q2))/kappa(Q2,2*p1p2 + m3**2,0d0))
      udip = m1**2 - 2*(((-m3**2 + Q2)*(papb + m1**2))/(2*Q2) + (kappa(Q2,m3**2,0d0)*(pap3 - ((pap3 + pbp3)*(papb + m1**2))/Q2))/kappa(Q2,2*p1p2 + m3**2,0d0))

      tree = NeuQ2qx_M2TreeDipole(sdip,tdip,udip)

      yijk = p1p2/(p1p2 + p1p3 + p2p3)
      vijk = 1d0
      viji = 1d0
      zP = ((1-muQ2)*yijk)/(muQ2+(1-muQ2)*yijk)*(1+viji*vijk)/2d0
      zM = ((1-muQ2)*yijk)/(muQ2+(1-muQ2)*yijk)*(1-viji*vijk)/2d0
      vijkTilde = kappa(1d0,muQ2,0d0)/(1-muQ2)

      zj = p1p3/(p1p3 + p2p3)

      dipcolor = -CV/(2*CF)
      NeuQ2qx_M2RadiationDipoleG = NeuQ2qx_M2RadiationDipoleG + 0.5d0*(CF*dipcolor*gs**2*tree*(-((vijkTilde*(1 + m3**2/p1p2 + zj))/vijk) + 2/(1 - (1 - yijk)*zj)))/p1p2

! -----------------------------------------------    Dipoles with gluon as emitter     ----------------------------------------------------------------------------------

! (2) emission off gluon (g->g(p2)+g(p3))
      NeuQ2qx_M2RadiationDipoleG = NeuQ2qx_M2RadiationDipoleG +0.5d0*NeuQ2qx_M2_FEFS_gg()

! ****************************************************** FE-IS and IE-FS ***************************************************************************************

! --------------------------------------------    Dipoles with gluon as emitter/spectator  ----------------------------------------------------------------------

      xdip = (-p2p3 + pbp2 + pbp3)/(pbp2 + pbp3) 

      if (xdip.gt.x0) then

      tdip = m1**2 + m3**2 - 2*pap1 
      Rdip = dsqrt((tdip - m2**2 + 2*m2**2*xdip)**2 - 4*m2**2*tdip*xdip**2)/(-tdip+m2**2) 
      sdip = m1**2 + m2**2 + ((xdip*(sman - m1**2 - m2**2))/Rdip +(tdip-m2**2 +2*m2**2*xdip)/(tdip*Rdip)*(m1**2 - pap1)-((tdip + m2**2)/tdip)*(m1**2 -pap1)) 
      udip = m1**2 + m2**2  + m3**2 - sdip - tdip 

      zP = (1+Rdip)/2d0
      zM = (1-Rdip)/2d0

      zi = pbp3/(pbp2 + pbp3)
      zj = 1 - zi

      tree = NeuQ2qx_M2TreeDipole(sdip,tdip,udip)

! (5) FE-IS: emission off gluon (g->g(p2)+g(p3)) 
      NeuQ2qx_M2RadiationDipoleG = NeuQ2qx_M2RadiationDipoleG + 0.5d0*NeuQ2qx_M2_FEIS_gg()

! (8.1) IE-FS: gluon p2 is spectator  (\tilde{q}-> \tilde{q}(pb)+g(p3)), 
      dipcolor =  -CV/(2*CF)
      NeuQ2qx_M2RadiationDipoleG = NeuQ2qx_M2RadiationDipoleG +  0.5d0*(CF*dipcolor*gs**2*tree*(-2 - (m2**2*xdip)/pbp3 + 2/(2 - xdip - zj)))/(pbp3*xdip)

! (8.2) IE-FS: gluon p3 is spectator  (\tilde{q}-> \tilde{q}(pb)+g(p2)), 
      zi = pbp2/(pbp2 + pbp3)
      zj = 1 - zi

      dipcolor =  -CV/(2*CF)
      NeuQ2qx_M2RadiationDipoleG = NeuQ2qx_M2RadiationDipoleG +  0.5d0*(CF*dipcolor*gs**2*tree*(-2 - (m2**2*xdip)/pbp2 + 2/(2 - xdip - zj)))/(pbp2*xdip)

      endif


! --------------------------------------------     Dipoles with quark as emitter/spectator and with gluon p3 as i ----------------------------------------------------------------------

      xdip = (-p1p3 + pbp1 + pbp3)/(pbp1 + pbp3) 

      if (xdip.gt.x0) then

      udip = m1**2 - 2*pap2
      Rdip = dsqrt(((udip - m2**2 - m3**2 + 2*m2**2*xdip)**2 - 4*m2**2*udip*xdip**2)/(udip**2 + m2**4 + m3**4 - 2*udip*m2**2 - 2*udip*m3**2 - 2*m2**2*m3**2))
      sdip = m1**2 + m2**2 + ((xdip*(sman - m1**2 - m2**2))/Rdip + (udip- m3**2 - m2**2 + 2*m2**2*xdip)/(udip*Rdip)*(m1**2 - pap2) - ((udip + m2**2 - m3**2)/udip)*(m1**2 - pap2))
      tdip = m1**2 + m2**2 + m3**2 + m4**2 - sdip - udip

      tree = NeuQ2qx_M2TreeDipole(sdip,tdip,udip)

      zi = pbp3/(pbp1 + pbp3)
      zj = 1 - zi

! (4.1) FE-IS: emission off quark (q->q(p1)+g(p3)) 
      dipcolor =  (CV - 2*CF)/(2*CF)
      NeuQ2qx_M2RadiationDipoleG = NeuQ2qx_M2RadiationDipoleG +  0.5d0*(CF*dipcolor*gs**2*tree*(-2 - m3**2/p1p3 + zi + 2/(2 - xdip - zj)))/(p1p3*xdip)


! (7.1) IE-FS: quark is spectator (\tilde{q}-> \tilde{q}(pb)+g(p3))
      dipcolor =  (CV - 2*CF)/(2*CF)
      NeuQ2qx_M2RadiationDipoleG = NeuQ2qx_M2RadiationDipoleG +  0.5d0*(CF*dipcolor*gs**2*tree*(-2 - (m2**2*xdip)/pbp3 + 2/(2 - xdip - zj)))/(pbp3*xdip)

      endif

! --------------------------------------------     Dipoles with quark as emitter/spectator and with gluon p2 as i ----------------------------------------------------------------------

      xdip = (-p1p2 + pbp1 + pbp2)/(pbp1 + pbp2)

      if (xdip.gt.x0) then

      udip = m1**2 - 2*pap3
      Rdip = dsqrt((-4*m2**2*udip*xdip**2 + (-m2**2 - m3**2 + udip + 2*m2**2*xdip)**2)/(-4*m2**2*m3**2 + (-m2**2 - m3**2 + udip)**2))
      sdip = m1**2 + m2**2 + 2*((papb*xdip)/Rdip + (m1**2 - pap3)*(-(m2**2 - m3**2 + udip)/(2*udip) + (-m2**2 - m3**2 + udip + 2*m2**2*xdip)/(2*Rdip*udip)))
      tdip = m2**2 - 2*((pbp3*xdip)/Rdip + pap3*(-(m2**2 - m3**2 + udip)/(2*udip) + (-m2**2 - m3**2 + udip + 2*m2**2*xdip)/(2*Rdip*udip)))

      tree = NeuQ2qx_M2TreeDipole(sdip,tdip,udip)

      zi = pbp2/(pbp1 + pbp2)
      zj = 1 - zi

! (4.2) FE-IS: emission off quark (q->q(p1)+g(p2)) 
      dipcolor =  (CV - 2*CF)/(2*CF)
      NeuQ2qx_M2RadiationDipoleG = NeuQ2qx_M2RadiationDipoleG + 0.5d0*(CF*dipcolor*gs**2*tree*(-2 - m3**2/p1p2 + zi + 2/(2 - xdip - zj)))/(p1p2*xdip)


! (7.2) IE-FS: quark is spectator (\tilde{q}-> \tilde{q}(pb)+g(p2))
      dipcolor =  (CV - 2*CF)/(2*CF)
      NeuQ2qx_M2RadiationDipoleG = NeuQ2qx_M2RadiationDipoleG +  0.5d0*(CF*dipcolor*gs**2*tree*(-2 - (m2**2*xdip)/pbp2 + 2/(2 - xdip - zj)))/(pbp2*xdip)

      endif

      end
      
      double precision function NeuQ2qx_M2RadiationDipoleEW()
! pa -> neutralino, pb -> squark, p1 -> quark, p2-> Higgs, p3 -> gluon
      implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

      double precision udip, sdip, tdip, Qbar2, Rdip,xdip, zi, zj, colorTree, NeuQ2qx_M2TreeDipole, kappa

      xdip = (pbp3 + pbp1 - p1p3)/(pbp3 + pbp1)

      if (xdip.gt.x0) then

      udip = m1**2 + m4**2 - 2*pap2 !Q^2
      Qbar2 = udip - m2**2 -m3**2
      Rdip = dsqrt(-4*m2**2*udip*xdip**2 + (Qbar2 + 2*m2**2*xdip)**2)/kappa(udip,m3**2,m2**2)
      sdip = m1**2 + m2**2 + (2*papb*xdip)/Rdip + ((m1**2 - pap2)*(Qbar2- Rdip*(m2**2 - m3**2 + udip) + 2*m2**2*xdip))/(Rdip*udip)
      tdip = m2**2 + m4**2 - (2*pbp2*xdip)/Rdip + ((m4**2 - pap2)*(Qbar2- Rdip*(m2**2 - m3**2 + udip) + 2*m2**2*xdip))/(Rdip*udip)

      zi = pbp3/(pbp1 + pbp3)
      zj = 1 - zi

      colorTree = -NeuQ2qx_M2TreeDipole(sdip,tdip,udip)

! D_{31}^b   
      NeuQ2qx_M2RadiationDipoleEW = (CF*colorTree*gs**2*(-2 + zi - m3**2/p1p3 + 2/(2 - xdip - zj)))/(p1p3*xdip)

! D^{b3}_1   
      NeuQ2qx_M2RadiationDipoleEW = NeuQ2qx_M2RadiationDipoleEW + (CF*colorTree*gs**2*(-2 - (m2**2*xdip)/pbp3 + 2/(2 - xdip - zj)))/(pbp3*xdip)

      else 

      NeuQ2qx_M2RadiationDipoleEW = 0d0

      endif

      end

      double precision function NeuQ2qx_M2TreeH_Dipole(sdip,tdip,udip)

      implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

      complex*16 DenSmZ, DenSmH(4)
      double precision colorfactor
      double precision DenSmZ1, DenSmH1(4), M2test
      double precision treeSS, treeTT, treeUU, treeST, treeTS, treeSU, treeUS, treeTU, treeUT
      integer sumU, sumUc, sumT, sumTc
      integer sumUend, sumUcend, sumTend, sumTcend
      double precision sdip, tdip,udip

      NeuQ2qx_M2TreeH_Dipole = 0d0
      treeSS = 0d0
      treeTT = 0d0
      treeUU = 0d0
      treeST = 0d0
      treeTS = 0d0
      treeSU = 0d0
      treeUS = 0d0
      treeUT = 0d0
      treeTU = 0d0

      sumU = 0
      sumUc = 0
      sumT = 0
      sumTc = 0


      if ((xtype.eq.5).or.(xtype.eq.6)) then

      sumTend = 1
      sumTcend = 1
      sumUend = 3
      sumUcend = 3

      elseif (xtype.eq.7) then

      sumTend = 0
      sumTcend = 0
      sumUend = 3
      sumUcend = 3

      elseif (xtype.eq.8) then

      sumTend = 1
      sumTcend = 1
      sumUend = 1
      sumUcend = 1

      endif



      ! ==================================================================================================================================	
      !                             ---                     Higgs final states           ---
      ! ==================================================================================================================================



! --- Mss  ---	

      treeSS = NC*(FR(xtype,pstype,psgen)*(dconjg(FL(xtype,pstype,psgen))*m3*(AL(chimeigen,sfmeigen,sftype,sfgen)*(2*m1*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*sdip +
     &  dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*Mf(pstype,psgen)*(m1**2 - m2**2 + sdip)) +
     &  Mf(pstype,psgen)*AR(chimeigen,sfmeigen,sftype,sfgen)*(2*dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*m1*Mf(pstype,psgen) +
     &  dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m1**2 - m2**2 + sdip))) +
     &  dconjg(FR(xtype,pstype,psgen))*(Mf(pstype,psgen)*AR(chimeigen,sfmeigen,sftype,sfgen)*(Mf(pstype,psgen)*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m1**2 + m3**2 - tdip) +
     &  dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*m1*(m1**2 + m2**2 + 2*m3**2 - tdip - udip)) +
     &  AL(chimeigen,sfmeigen,sftype,sfgen)*(m1*Mf(pstype,psgen)*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m1**2 + m2**2 + 2*m3**2 - tdip - udip) +
     &  dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(m1**4 - m2**4 + sdip*(m3**2 - udip) + m1**2*(2*m3**2 - tdip - udip) + m2**2*(-2*m3**2 + sdip + tdip + udip))))) +
     &  FL(xtype,pstype,psgen)*(m3*dconjg(FR(xtype,pstype,psgen))*(AR(chimeigen,sfmeigen,sftype,sfgen)*(m1**2*Mf(pstype,psgen)*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen)) +
     &  2*dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*m1*sdip + Mf(pstype,psgen)*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(-m2**2 + sdip)) +
     &  AL(chimeigen,sfmeigen,sftype,sfgen)*Mf(pstype,psgen)*(2*m1*Mf(pstype,psgen)*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen)) + dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(m1**2 - m2**2 + sdip))) +
     &  dconjg(FL(xtype,pstype,psgen))*(AL(chimeigen,sfmeigen,sftype,sfgen)*Mf(pstype,psgen)*(dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*Mf(pstype,psgen)*(m1**2 + m3**2 - tdip) +
     &  m1*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m1**2 + m2**2 + 2*m3**2 - tdip - udip)) + AR(chimeigen,sfmeigen,sftype,sfgen)*(dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*m1*Mf(pstype,psgen)*(m1**2 + m2**2 + 2*m3**2 - tdip - udip) +
     &  dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m1**4 - m2**4 + sdip*(m3**2 - udip) + m1**2*(2*m3**2 - tdip - udip) + m2**2*(-2*m3**2 + sdip + tdip + udip))))))/(Mf(pstype,psgen)**2 - sdip + (cI*Mf(pstype,psgen)*width(pstype)))**2


      NeuQ2qx_M2TreeH_Dipole = NeuQ2qx_M2TreeH_Dipole + treeSS



! --- Mtt  ---
      do sumT = 0,sumTend
      do sumTc = 0,sumTcend

      treeTT = treeTT + NC*(G(xtype,sfmeigen,ptmeigen+sumT,sftype,sfgen)*dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))*(AL(chimeigen,ptmeigen+sumT,pttype,ptgen)*(2*m1*m3*dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen)) +
     &	dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))*(m1**2 + m3**2 - tdip)) + AR(chimeigen,ptmeigen+sumT,pttype,ptgen)*(2*dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))*m1*m3 +
     &  dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))*(m1**2 + m3**2 - tdip))))/((MSf(ptmeigen+sumT,pttype,ptgen)**2 - tdip)*(MSf(ptmeigen+sumTc,pttype,ptgen)**2 - tdip))  


      enddo
      enddo


      NeuQ2qx_M2TreeH_Dipole = NeuQ2qx_M2TreeH_Dipole + treeTT


! --- Muu  ---
      do sumU = 0,sumUend
      do sumUc = 0,sumUcend

      treeUU = treeUU + NC*(HR(chimeigen,pumeigen+sumU,xtype)*(dconjg(HL(chimeigen,pumeigen+sumUc,xtype))*m1*(MNeuCha(pumeigen+sumU)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(2*dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*MNeuCha(pumeigen+sumUc) +
     &	dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sdip - tdip)) + AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*MNeuCha(pumeigen+sumUc)*(m1**2 + 2*m3**2 + m4**2 - sdip - tdip) +
     &  2*AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*m3*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*udip) + dconjg(HR(chimeigen,pumeigen+sumUc, xtype))*(AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1 - m4)*(m1 + m4)*(m3*MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen)) +
     &  dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sdip - tdip)) + AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m3*MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen)) + dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m3**2 + m4**2 - sdip))*udip +
     &  MNeuCha(pumeigen+sumU)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + m3**2 - tdip) + dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*(m1**2 - m4**2 + udip)))) +
     &  HL(chimeigen,pumeigen+sumU, xtype)*(m1*dconjg(HR(chimeigen,pumeigen+sumUc,xtype))*(AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*MNeuCha(pumeigen+sumU)*(2*m3*MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen)) + dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sdip - tdip)) +
     &  MNeuCha(pumeigen+sumUc)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sdip - tdip) +
     &  2*dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*udip) + dconjg(HL(chimeigen,pumeigen+sumUc, xtype))*((m1 - m4)*(m1 + m4)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*MNeuCha(pumeigen+sumUc) +
     &  dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sdip - tdip)) + AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*MNeuCha(pumeigen+sumUc) + dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m3**2 + m4**2 - sdip))*udip +
     &  AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*MNeuCha(pumeigen+sumU)*(dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*MNeuCha(pumeigen+sumUc)*(m1**2 + m3**2 - tdip) + m3*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 - m4**2 + udip)))))/((MNeuCha(pumeigen+sumU)**2 - udip)*(MNeuCha(pumeigen+sumUc)**2 - udip))



      enddo
      enddo


      NeuQ2qx_M2TreeH_Dipole = NeuQ2qx_M2TreeH_Dipole + treeUU



!--- Mst  ---

      do sumTc = 0,sumTcend

      treeST = treeST + NC*(dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))*(FL(xtype,pstype,psgen)*(AL(chimeigen,sfmeigen, sftype,sfgen)*Mf(pstype,psgen)*(2*m1*m3*dconjg(AR(chimeigen,ptmeigen+sumTc, pttype,ptgen)) +
     &	dconjg(AL(chimeigen,ptmeigen+sumTc, pttype,ptgen))*(m1**2 + m3**2 - tdip)) + AR(chimeigen,sfmeigen, sftype,sfgen)*(m3*dconjg(AR(chimeigen,ptmeigen+sumTc, pttype,ptgen))*(m1**2 - m2**2 + sdip) + 
     &  dconjg(AL(chimeigen,ptmeigen+sumTc, pttype,ptgen))*m1*(m1**2 + m2**2 + 2*m3**2 - tdip - udip))) + 
     &  FR(xtype,pstype,psgen)*(Mf(pstype,psgen)*AR(chimeigen,sfmeigen,sftype,sfgen)*(2*dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))*m1*m3 + 
     &  dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))*(m1**2 + m3**2 - tdip)) + AL(chimeigen,sfmeigen, sftype,sfgen)*(dconjg(AL(chimeigen,ptmeigen+sumTc, pttype,ptgen))*m3*(m1**2 - m2**2 + sdip) +
     &  m1*dconjg(AR(chimeigen,ptmeigen+sumTc, pttype,ptgen))*(m1**2 + m2**2 + 2*m3**2 - tdip - udip)))))/((Mf(pstype,psgen)**2 - sdip + (cI*Mf(pstype,psgen)*width(pstype)))*(MSf(ptmeigen+sumTc,pttype,ptgen)**2 - tdip))

      enddo

      NeuQ2qx_M2TreeH_Dipole = NeuQ2qx_M2TreeH_Dipole + treeST


! --- Mts  ---

      do sumT = 0,sumTend

      treeTS = treeTS + NC*(G(xtype,sfmeigen,ptmeigen+sumT,sftype,sfgen)*(dconjg(FL(xtype,pstype,psgen))*(dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*Mf(pstype,psgen)*(2*m1*m3*AR(chimeigen,ptmeigen+sumT,pttype,ptgen) +
     &	AL(chimeigen,ptmeigen+sumT,pttype,ptgen)*(m1**2 + m3**2 - tdip)) + dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m3*AR(chimeigen,ptmeigen+sumT,pttype,ptgen)*(m1**2 - m2**2 + sdip) +
     &  AL(chimeigen,ptmeigen+sumT,pttype,ptgen)*m1*(m1**2 + m2**2 + 2*m3**2 - tdip - udip))) +
     &  dconjg(FR(xtype,pstype,psgen))*(Mf(pstype,psgen)*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(2*AL(chimeigen,ptmeigen+sumT,pttype,ptgen)*m1*m3 + AR(chimeigen,ptmeigen+sumT,pttype,ptgen)*(m1**2 + m3**2 - tdip)) +
     &  dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(AL(chimeigen,ptmeigen+sumT,pttype,ptgen)*m3*(m1**2 - m2**2 + sdip) +
     &  m1*AR(chimeigen,ptmeigen+sumT,pttype,ptgen)*(m1**2 + m2**2 + 2*m3**2 - tdip - udip)))))/((Mf(pstype,psgen)**2 - sdip + (cI*Mf(pstype,psgen)*width(pstype)))*(MSf(ptmeigen+sumT,pttype,ptgen)**2 - tdip))

      enddo


      NeuQ2qx_M2TreeH_Dipole = NeuQ2qx_M2TreeH_Dipole + treeTS

! --- Msu  ---

      do sumUc = 0,sumUcend

      treeSU = treeSU + NC*(dconjg(HL(chimeigen,pumeigen+sumUc,xtype))*(2*FL(xtype,pstype,psgen)*(AL(chimeigen,sfmeigen,sftype,sfgen)*Mf(pstype,psgen)*(dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*MNeuCha(pumeigen+sumUc)*(m1**2 + m3**2 - tdip) +
     &	m3*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 - m4**2 + udip)) + m1*AR(chimeigen,sfmeigen,sftype,sfgen)*(dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*MNeuCha(pumeigen+sumUc)*(m1**2 + m2**2 + 2*m3**2 - tdip - udip) +
     &  m3*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(-2*m2**2 - 2*m4**2 + sdip + tdip + udip))) + FR(xtype,pstype,psgen)*(2*m1*Mf(pstype,psgen)*AR(chimeigen,sfmeigen,sftype,sfgen)*(2*dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*MNeuCha(pumeigen+sumUc) +
     &  dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sdip - tdip)) + AL(chimeigen,sfmeigen,sftype,sfgen)*(2*dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*MNeuCha(pumeigen+sumUc)*(m1**2 - m2**2 + sdip) +
     &  dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(2*m1**4 + (m3**2 + m4**2 - sdip - tdip)*(sdip - tdip) + m1**2*(2*m2**2 + 4*m3**2 + 2*m4**2 - sdip - 3*tdip - udip) + (m3**2 + m4**2)*udip - udip**2 + m2**2*(-2*m4**2 + sdip - tdip + udip))))) +
     &  dconjg(HR(chimeigen,pumeigen+sumUc,xtype))*(2*FR(xtype,pstype,psgen)*(Mf(pstype,psgen)*AR(chimeigen,sfmeigen,sftype,sfgen)*(MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + m3**2 - tdip) +
     &  dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*(m1**2 - m4**2 + udip)) + AL(chimeigen,sfmeigen,sftype,sfgen)*m1*(MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + m2**2 + 2*m3**2 - tdip - udip) +
     &  dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*(-2*m2**2 - 2*m4**2 + sdip + tdip + udip))) + FL(xtype,pstype,psgen)*(2*AL(chimeigen,sfmeigen,sftype,sfgen)*m1*Mf(pstype,psgen)*(2*m3*MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen)) +
     &  dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sdip - tdip)) + AR(chimeigen,sfmeigen,sftype,sfgen)*(2*m3*MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 - m2**2 + sdip) +
     &  dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(2*m1**4 + (m3**2 + m4**2 - sdip - tdip)*(sdip - tdip) + m1**2*(2*m2**2 + 4*m3**2 + 2*m4**2 - sdip - 3*tdip - udip) +
     &  (m3**2 + m4**2)*udip - udip**2 + m2**2*(-2*m4**2 + sdip - tdip + udip))))))/(2.*(Mf(pstype,psgen)**2 - sdip + (cI*Mf(pstype,psgen)*width(pstype)))*(MNeuCha(pumeigen+sumUc)**2 - udip))

      enddo


      NeuQ2qx_M2TreeH_Dipole = NeuQ2qx_M2TreeH_Dipole + treeSU

! --- Mus  ---

      do sumU = 0,sumUend

      treeUS = treeUS + NC*(HR(chimeigen,pumeigen+sumU,xtype)*(2*dconjg(FR(xtype,pstype,psgen))*(Mf(pstype,psgen)*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(MNeuCha(pumeigen+sumU)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 + m3**2 - tdip) + 
     &	AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*m3*(m1**2 - m4**2 + udip)) + dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*m1*(MNeuCha(pumeigen+sumU)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 + m2**2 + 2*m3**2 - tdip - udip) +
     &  AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*m3*(-2*m2**2 - 2*m4**2 + sdip + tdip + udip))) + dconjg(FL(xtype,pstype,psgen))*(2*dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*m1*Mf(pstype,psgen)*(2*m3*MNeuCha(pumeigen+sumU)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen) +
     &  AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 + 2*m3**2 + m4**2 - sdip - tdip)) + dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(2*m3*MNeuCha(pumeigen+sumU)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 - m2**2 + sdip) +
     &  AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*(2*m1**4 + (m3**2 + m4**2 - sdip - tdip)*(sdip - tdip) + m1**2*(2*m2**2 + 4*m3**2 + 2*m4**2 - sdip - 3*tdip - udip) + (m3**2 + m4**2)*udip - udip**2 + m2**2*(-2*m4**2 + sdip - tdip + udip))))) +
     &  HL(chimeigen,pumeigen+sumU,xtype)*(2*dconjg(FL(xtype,pstype,psgen))*(dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*Mf(pstype,psgen)*(AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*MNeuCha(pumeigen+sumU)*(m1**2 + m3**2 - tdip) +
     &  m3*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 - m4**2 + udip)) + m1*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*MNeuCha(pumeigen+sumU)*(m1**2 + m2**2 + 2*m3**2 - tdip - udip) +
     &  m3*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(-2*m2**2 - 2*m4**2 + sdip + tdip + udip))) + dconjg(FR(xtype,pstype,psgen))*(2*m1*Mf(pstype,psgen)*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(2*AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*m3*MNeuCha(pumeigen+sumU) +
     &  AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 + 2*m3**2 + m4**2 - sdip - tdip)) + dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(2*AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*m3*MNeuCha(pumeigen+sumU)*(m1**2 - m2**2 + sdip) +
     &  AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(2*m1**4 + (m3**2 + m4**2 - sdip - tdip)*(sdip - tdip) + m1**2*(2*m2**2 + 4*m3**2 + 2*m4**2 - sdip - 3*tdip - udip) +
     &  (m3**2 + m4**2)*udip - udip**2 + m2**2*(-2*m4**2 + sdip - tdip + udip))))))/(2.*(Mf(pstype,psgen)**2 - sdip + (cI*Mf(pstype,psgen)*width(pstype)))*(MNeuCha(pumeigen+sumU)**2 - udip))

      enddo


      NeuQ2qx_M2TreeH_Dipole = NeuQ2qx_M2TreeH_Dipole + treeUS

      do sumT = 0,sumTend
      do sumUc = 0,sumUcend

      treeTU = treeTU + NC*(G(xtype,sfmeigen,ptmeigen+sumT,sftype,sfgen)*(dconjg(HR(chimeigen,pumeigen+sumUc,xtype))*(AL(chimeigen,ptmeigen+sumT,pttype,ptgen)*m1*(2*m3*MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen)) +
     &  dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sdip - tdip)) + AR(chimeigen,ptmeigen+sumT,pttype,ptgen)*(MNeuCha(pumeigen+sumUc)*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + m3**2 - tdip) +
     &  dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*(m1**2 - m4**2 + udip))) + dconjg(HL(chimeigen,pumeigen+sumUc,xtype))*(m1*AR(chimeigen,ptmeigen+sumT,pttype,ptgen)*(2*dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*m3*MNeuCha(pumeigen+sumUc) +
     &  dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 + 2*m3**2 + m4**2 - sdip - tdip)) + AL(chimeigen,ptmeigen+sumT,pttype,ptgen)*(dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))*MNeuCha(pumeigen+sumUc)*(m1**2 + m3**2 - tdip) +
     &  m3*dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))*(m1**2 - m4**2 + udip)))))/((MSf(ptmeigen+sumT,pttype,ptgen)**2 - tdip)*(MNeuCha(pumeigen+sumUc)**2 - udip))

      enddo
      enddo


      NeuQ2qx_M2TreeH_Dipole = NeuQ2qx_M2TreeH_Dipole + treeTU



! --- Mut  ---

      do sumTc = 0,sumTcend
      do sumU = 0,sumUend

      treeUT = treeUT + NC*(dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))*(HR(chimeigen,pumeigen+sumU,xtype)*(dconjg(AL(chimeigen, ptmeigen+sumTc,pttype,ptgen))*m1*(2*m3*MNeuCha(pumeigen+sumU)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen) +
     &	AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 + 2*m3**2 + m4**2 - sdip - tdip)) + dconjg(AR(chimeigen, ptmeigen+sumTc,pttype,ptgen))*(MNeuCha(pumeigen+sumU)*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 + m3**2 - tdip) +
     &  AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*m3*(m1**2 - m4**2 + udip))) + HL(chimeigen,pumeigen+sumU,xtype)*(m1*dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))*(2*AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*m3*MNeuCha(pumeigen+sumU) +
     &  AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 + 2*m3**2 + m4**2 - sdip - tdip)) + dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))*(AL(pumeigen+sumU,sfmeigen,sftype,sfgen)*MNeuCha(pumeigen+sumU)*(m1**2 + m3**2 - tdip) +
     &  m3*AR(pumeigen+sumU,sfmeigen,sftype,sfgen)*(m1**2 - m4**2 + udip)))))/((MSf(ptmeigen+sumTc,pttype,ptgen)**2 - tdip)*(MNeuCha(pumeigen+sumU)**2 - udip))

      enddo
      enddo


      NeuQ2qx_M2TreeH_Dipole = NeuQ2qx_M2TreeH_Dipole + treeUT


      end


      
      
      
      double precision function NeuQ2qx_M2RadiationDipoleGqq()

      implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

      double precision NeuQ2qx_M2_FEIS_qqbar,NeuQ2qx_M2_FEFS_qqbar,Nf,xdip


! Number of massless quarks
      Nf = 4d0
! ********************************************************** FE-FS  (from hep-ph/0201036v1 Section 5.1.1+5.1.2 ***********************************************************

! (3) gluon decay into massless qqbar pair (g->q+qbar)
      NeuQ2qx_M2RadiationDipoleGqq = Nf*NeuQ2qx_M2_FEFS_qqbar()


! ****************************************************** FE-IS  ***************************************************************************************

! --------------------------------------------    Dipole with gluon as emitter ----------------------------------------------------------------------
      
      xdip = (-p2p3 + pbp2 + pbp3)/(pbp2 + pbp3) 

      if (xdip.gt.x0) then

! (6) FE-IS: gluon decay into massless qqbar pair (g->q+qbar) 
      NeuQ2qx_M2RadiationDipoleGqq = NeuQ2qx_M2RadiationDipoleGqq + Nf*NeuQ2qx_M2_FEIS_qqbar()
      endif

      end

      double precision function NeuQ2qx_M2TreeGA_dipole(sdip,tdip,udip)
! return the matrix element for a gluon or a photon in the final state
      implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"



      integer i,j
      double complex DenSmZ, DenSmH(4)
      double precision colorfactor, DenSmZ1, DenSmH1(4),treeSS, treeTT, treeST, treeTS
      double precision sdip,tdip,udip

      NeuQ2qx_M2TreeGA_dipole	= 0d0
      treeSS = 0d0
      treeTT = 0d0
      treeST = 0d0
      treeTS = 0d0
  

! ==================================================================================================================================	
!                             ---                     gluon or photon final state           ---
! ==================================================================================================================================


! --- 	setting the specific couplings B,C for gluon or photon final state   ---
      if (xtype.eq.1) then
            colorfactor = 4./3.	 
      else if (xtype.eq.2) then
            colorfactor = 1.0
      endif

! --- Mss for gluon ---	
      treeSS = colorfactor*NC*((2*B(xtype,pstype)**2*(AL(chimeigen,sfmeigen,sftype,sfgen)*(2*m1*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m1**2*Mf(pstype,psgen) +
     & m2**2*Mf(pstype,psgen) + 2*m3**2*Mf(pstype,psgen) - 2*m3*Mf(pstype,psgen)**2 - 2*m3*sdip - Mf(pstype,psgen)*tdip - Mf(pstype,psgen)*udip) + 
     & dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(m1**4 - m2**4 + m3**2*Mf(pstype,psgen)**2 + m3**2*sdip - 4*m3*Mf(pstype,psgen)*sdip -
     & Mf(pstype,psgen)**2*tdip + m1**2*(2*m3**2 - 4*m3*Mf(pstype,psgen) + Mf(pstype,psgen)**2 - tdip - udip) - sdip*udip + m2**2*(-2*m3**2 + 4*m3*Mf(pstype,psgen) +
     & sdip + tdip + udip))) + AR(chimeigen,sfmeigen,sftype,sfgen)*(2*dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*m1*(m1**2*Mf(pstype,psgen) +
     & m2**2*Mf(pstype,psgen) + 2*m3**2*Mf(pstype,psgen) - 2*m3*Mf(pstype,psgen)**2 - 2*m3*sdip - Mf(pstype,psgen)*tdip - Mf(pstype,psgen)*udip) +
     & dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(m1**4 - m2**4 + m3**2*Mf(pstype,psgen)**2 + m3**2*sdip - 4*m3*Mf(pstype,psgen)*sdip - Mf(pstype,psgen)**2*tdip +
     & m1**2*(2*m3**2 - 4*m3*Mf(pstype,psgen) + Mf(pstype,psgen)**2 - tdip - udip) - sdip*udip + m2**2*(-2*m3**2 + 4*m3*Mf(pstype,psgen) + sdip + tdip + udip)))))/((Mf(pstype,psgen)**2 - sdip + (cI*Mf(pstype,psgen)*width(pstype)))**2))

      NeuQ2qx_M2TreeGA_dipole = NeuQ2qx_M2TreeGA_dipole + treeSS

! --- Mtt for gluon ---
      treeTT = colorfactor*NC*(-1.*C(xtype,pstype)**2*(AL(chimeigen,ptmeigen,pttype,ptgen)*(2*m1*m3*dconjg(AR(chimeigen,ptmeigen,pttype,ptgen)) + 
     & dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))*(m1**2 + m3**2 - tdip)) + 
     & AR(chimeigen,ptmeigen,pttype,ptgen)*(2*dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))*m1*m3 +
     & dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))*(m1**2 + m3**2 - tdip)))*(2*m2**2 - m4**2 + 2*tdip))/((MSf(ptmeigen,pttype,ptgen)**2 - tdip)**2)

      NeuQ2qx_M2TreeGA_dipole = NeuQ2qx_M2TreeGA_dipole + treeTT

! --- Mst for gluon ---
      treeST = colorfactor*NC*((-1.*B(xtype,pstype)*C(xtype,pstype)*(AL(chimeigen,sfmeigen,sftype,sfgen)*(2*m1*dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))*(-3*m1**2*m3 -
     & 2*m3*m4**2 + 3*m3**2*Mf(pstype,psgen) + m4**2*Mf(pstype,psgen) + m2**2*(m3 + 2*Mf(pstype,psgen)) + 2*m3*sdip - Mf(pstype,psgen)*sdip + m3*tdip + m3*udip -
     & 2*Mf(pstype,psgen)*udip) + dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))*(-2*m3*m4**2*Mf(pstype,psgen) + m3**2*sdip + m4**2*sdip +
     & 4*m3*Mf(pstype,psgen)*sdip - sdip**2 + m3**2*tdip + m4**2*tdip - tdip**2 + m1**2*(8*m2**2 + 6*m3*(m3 - Mf(pstype,psgen)) - sdip + tdip - 5*udip) +
     & m2**2*(2*m3*(m3 - 2*Mf(pstype,psgen)) + sdip - 3*tdip - udip) - (m3**2 + m4**2 - 2*m3*Mf(pstype,psgen))*udip + udip**2)) +
     & AR(chimeigen,sfmeigen,sftype,sfgen)*(2*dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))*m1*(-3*m1**2*m3 - 2*m3*m4**2 +
     & 3*m3**2*Mf(pstype,psgen) + m4**2*Mf(pstype,psgen) + m2**2*(m3 + 2*Mf(pstype,psgen)) + 2*m3*sdip - Mf(pstype,psgen)*sdip + m3*tdip + m3*udip -
     & 2*Mf(pstype,psgen)*udip) + dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))*(-2*m3*m4**2*Mf(pstype,psgen) + m3**2*sdip + m4**2*sdip +
     & 4*m3*Mf(pstype,psgen)*sdip - sdip**2 + m3**2*tdip + m4**2*tdip - tdip**2 + m1**2*(8*m2**2 + 6*m3*(m3 - Mf(pstype,psgen)) - sdip + tdip - 5*udip) +
     & m2**2*(2*m3*(m3 - 2*Mf(pstype,psgen)) + sdip - 3*tdip - udip) - (m3**2 + m4**2 - 2*m3*Mf(pstype,psgen))*udip + udip**2))))/(2.*(Mf(pstype,psgen)**2 -
     & sdip + (cI*Mf(pstype,psgen)*width(pstype)))*(MSf(ptmeigen,pttype,ptgen)**2 - tdip))) 

      NeuQ2qx_M2TreeGA_dipole = NeuQ2qx_M2TreeGA_dipole + treeST

 
! --- Mts for gluon ---
      treeTS = colorfactor*NC*((-1.*B(xtype,pstype)*C(xtype,pstype)*(dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(2*m1*AR(chimeigen,ptmeigen,pttype,ptgen)*(-3*m1**2*m3 -
     & 2*m3*m4**2 + 3*m3**2*Mf(pstype,psgen) + m4**2*Mf(pstype,psgen) + m2**2*(m3 + 2*Mf(pstype,psgen)) + 2*m3*sdip - Mf(pstype,psgen)*sdip + m3*tdip + m3*udip -
     & 2*Mf(pstype,psgen)*udip) + AL(chimeigen,ptmeigen,pttype,ptgen)*(-2*m3*m4**2*Mf(pstype,psgen) + m3**2*sdip + m4**2*sdip + 4*m3*Mf(pstype,psgen)*sdip -
     & sdip**2 + m3**2*tdip + m4**2*tdip - tdip**2 + m1**2*(8*m2**2 + 6*m3*(m3 - Mf(pstype,psgen)) - sdip + tdip - 5*udip) + m2**2*(2*m3*(m3 - 2*Mf(pstype,psgen)) +
     & sdip - 3*tdip - udip) - (m3**2 + m4**2 - 2*m3*Mf(pstype,psgen))*udip + udip**2)) +
     & dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(2*AL(chimeigen,ptmeigen,pttype,ptgen)*m1*(-3*m1**2*m3 - 2*m3*m4**2 +
     & 3*m3**2*Mf(pstype,psgen) + m4**2*Mf(pstype,psgen) + m2**2*(m3 + 2*Mf(pstype,psgen)) + 2*m3*sdip - Mf(pstype,psgen)*sdip + m3*tdip + m3*udip -
     & 2*Mf(pstype,psgen)*udip) + AR(chimeigen,ptmeigen,pttype,ptgen)*(-2*m3*m4**2*Mf(pstype,psgen) + m3**2*sdip + m4**2*sdip + 4*m3*Mf(pstype,psgen)*sdip -
     & sdip**2 + m3**2*tdip + m4**2*tdip - tdip**2 + m1**2*(8*m2**2 + 6*m3*(m3 - Mf(pstype,psgen)) - sdip + tdip - 5*udip) +  m2**2*(2*m3*(m3 - 2*Mf(pstype,psgen)) +
     & sdip - 3*tdip - udip) - (m3**2 + m4**2 - 2*m3*Mf(pstype,psgen))*udip + udip**2))))/(2.*(Mf(pstype,psgen)**2 - sdip + (cI*Mf(pstype,psgen)*width(pstype)))*(MSf(ptmeigen,pttype,ptgen)**2 - tdip)))

      NeuQ2qx_M2TreeGA_dipole = NeuQ2qx_M2TreeGA_dipole + treeTS

      end


      double precision function NeuQ2qx_M2TreeWZ_dipole(sdip,tdip,udip)
! returns tree level matrix element with a W or Z boson in the final state
      implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

      double precision colorfactor,sdip,tdip,udip
      double precision treeSS, treeTT, treeUU, treeST, treeTS, treeSU, treeUS, treeTU, treeUT
      integer sumT, sumTc, sumU, sumUc,k
        

      if (xtype.eq.3) then
            k = 3
      else if (xtype.eq.4) then
            k = 1
      endif


      colorfactor = 1d0

      NeuQ2qx_M2TreeWZ_dipole = 0d0
      treeSS = 0d0
      treeTT = 0d0
      treeUU = 0d0
      treeST = 0d0
      treeTS = 0d0
      treeSU = 0d0
      treeUS = 0d0
      treeTU = 0d0
      treeUT = 0d0



! ==================================================================================================================================	
!                             ---                     vector boson final states           ---
! ==================================================================================================================================

	
! --- Mss  ---	


        treeSS = (-1.*colorfactor*NC*(AL(chimeigen,sfmeigen,sftype,sfgen)*(dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*
     &          (BR(xtype,ftype)*((m1**4*(m3**2 - 1.*sdip) + m2**4*(-1.*m3**2 + sdip) + sdip*(m4**4 + (-1.*m3**2 + sdip)*tdip - 1.*m4**2*(sdip + tdip - 1.*udip)) + 
     &                  m1**2*(2.*m4**4 - 2.*m4**2*sdip - 1.*(m3**2 - 1.*sdip)*(tdip + udip)) + m2**2*(-2.*m4**4 + 2.*m4**2*sdip + (m3**2 - 1.*sdip)*(sdip + tdip + udip)))*
     &                dconjg(BR(xtype,ftype)) + 3.*m3*m4**2*(m1**2 - 1.*m2**2 + sdip)*dconjg(BL(xtype,ftype))*Mf(pstype,psgen)) + 
     &            BL(xtype,ftype)*Mf(pstype,psgen)*(3.*m3*m4**2*(m1**2 - 1.*m2**2 + sdip)*dconjg(BR(xtype,ftype)) + 
     &               (m4**4 + m1**2*(m3**2 - 1.*sdip) + (-1.*m3**2 + sdip)*udip - 1.*m4**2*(sdip - 1.*tdip + udip))*dconjg(BL(xtype,ftype))*Mf(pstype,psgen))) + 
     &         m1*dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(BL(xtype,ftype)*Mf(pstype,psgen)*
     &             ((2.*m4**4 + m1**2*(m3**2 - 1.*sdip) + m2**2*(m3**2 - 1.*sdip) - 2.*m4**2*sdip - 1.*m3**2*tdip + sdip*tdip - 1.*m3**2*udip + sdip*udip)*
     &                dconjg(BL(xtype,ftype)) + 6.*m3*m4**2*dconjg(BR(xtype,ftype))*Mf(pstype,psgen)) + 
     &            BR(xtype,ftype)*(6.*m3*m4**2*sdip*dconjg(BL(xtype,ftype)) + 
     &               (2.*m4**4 + m1**2*(m3**2 - 1.*sdip) + m2**2*(m3**2 - 1.*sdip) - 2.*m4**2*sdip - 1.*m3**2*tdip + sdip*tdip - 1.*m3**2*udip + sdip*udip)*dconjg(BR(xtype,ftype))*
     &                Mf(pstype,psgen)))) + AR(chimeigen,sfmeigen,sftype,sfgen)*
     &       (m1*dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(BR(xtype,ftype)*Mf(pstype,psgen)*
     &             ((2.*m4**4 + m1**2*(m3**2 - 1.*sdip) + m2**2*(m3**2 - 1.*sdip) - 2.*m4**2*sdip - 1.*m3**2*tdip + sdip*tdip - 1.*m3**2*udip + sdip*udip)*
     &                dconjg(BR(xtype,ftype)) + 6.*m3*m4**2*dconjg(BL(xtype,ftype))*Mf(pstype,psgen)) + 
     &            BL(xtype,ftype)*(6.*m3*m4**2*sdip*dconjg(BR(xtype,ftype)) + 
     &               (2.*m4**4 + m1**2*(m3**2 - 1.*sdip) + m2**2*(m3**2 - 1.*sdip) - 2.*m4**2*sdip - 1.*m3**2*tdip + sdip*tdip - 1.*m3**2*udip + sdip*udip)*dconjg(BL(xtype,ftype))*
     &                Mf(pstype,psgen))) + dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*
     &          (BL(xtype,ftype)*((m1**4*(m3**2 - 1.*sdip) + m2**4*(-1.*m3**2 + sdip) + sdip*(m4**4 + (-1.*m3**2 + sdip)*tdip - 1.*m4**2*(sdip + tdip - 1.*udip)) + 
     &                  m1**2*(2.*m4**4 - 2.*m4**2*sdip - 1.*(m3**2 - 1.*sdip)*(tdip + udip)) + m2**2*(-2.*m4**4 + 2.*m4**2*sdip + (m3**2 - 1.*sdip)*(sdip + tdip + udip)))*
     &                dconjg(BL(xtype,ftype)) + 3.*m3*m4**2*(m1**2 - 1.*m2**2 + sdip)*dconjg(BR(xtype,ftype))*Mf(pstype,psgen)) + 
     &            BR(xtype,ftype)*Mf(pstype,psgen)*(3.*m3*m4**2*(m1**2 - 1.*m2**2 + sdip)*dconjg(BL(xtype,ftype)) + 
     &               (m4**4 + m1**2*(m3**2 - 1.*sdip) + (-1.*m3**2 + sdip)*udip - 1.*m4**2*(sdip - 1.*tdip + udip))*dconjg(BR(xtype,ftype))*Mf(pstype,psgen))))))/
     &  (m4**2*(sdip - 1.*Mf(pstype,psgen)**2)**2)

        NeuQ2qx_M2TreeWZ_dipole = NeuQ2qx_M2TreeWZ_dipole + treeSS
	
	
! --- Mtt  ---
 
        do sumT = 0, 1
          do sumTc = 0, 1
     
        treeTT = (-1.*colorfactor*NC*(m2**4 + (m4**2 - 1.*tdip)**2 - 2.*m2**2*(m4**2 + tdip))*CL(ptmeigen + sumT,sfmeigen,sftype,sfgen,xtype)*
     &    (-1.*AL(chimeigen,ptmeigen + sumT,pttype,fgen)*((m1**2 + m3**2 - 1.*tdip)*dconjg(AL(chimeigen,ptmeigen + sumTc,pttype,fgen)) + 
     &         2.*m1*m3*dconjg(AR(chimeigen,ptmeigen + sumTc,pttype,fgen))) - 
     &      1.*AR(chimeigen,ptmeigen + sumT,pttype,fgen)*(2.*m1*m3*dconjg(AL(chimeigen,ptmeigen + sumTc,pttype,fgen)) + 
     &         (m1**2 + m3**2 - 1.*tdip)*dconjg(AR(chimeigen,ptmeigen + sumTc,pttype,fgen))))*dconjg(CL(ptmeigen + sumTc,sfmeigen,sftype,sfgen,xtype)))/
     &  (m4**2*(tdip - 1.*MSf(ptmeigen + sumT,pttype,ptgen)**2)*(tdip - 1.*MSf(ptmeigen + sumTc,pttype,ptgen)**2))

        
        NeuQ2qx_M2TreeWZ_dipole = NeuQ2qx_M2TreeWZ_dipole + treeTT


          enddo
        enddo
        

! --- Muu  ---

        do sumU = 0, k
          do sumUc = 0, k
 
     
        treeUU = NC*(-1.*(AR(pumeigen + sumU,sfmeigen,sftype,sfgen)*(dconjg(AR(pumeigen + sumUc,sfmeigen,sftype,sfgen))*
     &          (DR(chimeigen,pumeigen + sumU,sftype,xtype)*(((m1**6 - 2.*m4**6 + 2.*m4**4*sdip + 2.*m4**4*tdip + m1**4*(2.*m3**2 + 2.*m4**2 - 1.*sdip - 1.*tdip - 2.*udip) + 
     &                    2.*m4**4*udip - 2.*m4**2*sdip*udip - 1.*tdip*udip**2 + m3**2*(-4.*m4**4 + 2.*m4**2*udip + udip**2) - 
     &                    1.*m1**2*(m4**4 + m4**2*(sdip + tdip + udip) - 1.*udip*(sdip + 2.*tdip + udip) + m3**2*(-2.*m4**2 + 3.*udip)))*
     &                  dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)))/m4**2 - 
     &               3.*m1*(m1**2 + 2.*m3**2 + m4**2 - 1.*sdip - 1.*tdip)*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)) + 
     &            DL(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)*
     &             (-3.*m1*(m1**2 + 2.*m3**2 + m4**2 - 1.*sdip - 1.*tdip)*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &               (1.*(m4**4 + m1**2*(m3**2 - 1.*sdip) + (-1.*m3**2 + sdip)*udip - 1.*m4**2*(sdip - 1.*tdip + udip))*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*
     &                  MNeuCha(pumeigen + sumUc))/m4**2)) + (m3*dconjg(AL(pumeigen + sumUc,sfmeigen,sftype,sfgen))*
     &            (DL(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)*
     &               ((m1**4 - 2.*m4**4 + m1**2*(m4**2 - 2.*udip) + m4**2*udip + udip**2)*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &                 6.*m1*m4**2*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)) + 
     &              DR(chimeigen,pumeigen + sumU,sftype,xtype)*(-6.*m1*m4**2*udip*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)) + 
     &                 (m1**4 - 2.*m4**4 + m1**2*(m4**2 - 2.*udip) + m4**2*udip + udip**2)*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc))))/m4**2) + 
     &      AL(pumeigen + sumU,sfmeigen,sftype,sfgen)*((m3*dconjg(AR(pumeigen + sumUc,sfmeigen,sftype,sfgen))*
     &            (DR(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)*
     &               ((m1**4 - 2.*m4**4 + m1**2*(m4**2 - 2.*udip) + m4**2*udip + udip**2)*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &                 6.*m1*m4**2*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)) + 
     &              DL(chimeigen,pumeigen + sumU,sftype,xtype)*(-6.*m1*m4**2*udip*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)) + 
     &                 (m1**4 - 2.*m4**4 + m1**2*(m4**2 - 2.*udip) + m4**2*udip + udip**2)*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc))))/m4**2 + 
     &         dconjg(AL(pumeigen + sumUc,sfmeigen,sftype,sfgen))*(DL(chimeigen,pumeigen + sumU,sftype,xtype)*
     &             (((m1**6 - 2.*m4**6 + 2.*m4**4*sdip + 2.*m4**4*tdip + m1**4*(2.*m3**2 + 2.*m4**2 - 1.*sdip - 1.*tdip - 2.*udip) + 2.*m4**4*udip - 2.*m4**2*sdip*udip - 
     &                    1.*tdip*udip**2 + m3**2*(-4.*m4**4 + 2.*m4**2*udip + udip**2) - 
     &                    1.*m1**2*(m4**4 + m4**2*(sdip + tdip + udip) - 1.*udip*(sdip + 2.*tdip + udip) + m3**2*(-2.*m4**2 + 3.*udip)))*
     &                  dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)))/m4**2 - 
     &               3.*m1*(m1**2 + 2.*m3**2 + m4**2 - 1.*sdip - 1.*tdip)*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)) + 
     &            DR(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)*
     &             (-3.*m1*(m1**2 + 2.*m3**2 + m4**2 - 1.*sdip - 1.*tdip)*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &               (1.*(m4**4 + m1**2*(m3**2 - 1.*sdip) + (-1.*m3**2 + sdip)*udip - 1.*m4**2*(sdip - 1.*tdip + udip))*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*
     &                  MNeuCha(pumeigen + sumUc))/m4**2)))))/((-1.*udip + MNeuCha(pumeigen + sumU)**2)*(udip - 1.*MNeuCha(pumeigen + sumUc)**2))
        
     
 	NeuQ2qx_M2TreeWZ_dipole = NeuQ2qx_M2TreeWZ_dipole + treeUU

          enddo
        enddo


     	
! --- Mst  ---

        do sumTc = 0, 1

        treeST =  (0.5*colorfactor*NC*dconjg(CL(ptmeigen + sumTc,sfmeigen,sftype,sfgen,xtype))*
     &    (AL(chimeigen,sfmeigen,sftype,sfgen)*(2.*m1*dconjg(AR(chimeigen,ptmeigen + sumTc,pttype,ptgen))*
     &          (m3*(m2**4 + 2.*m4**4 - 2.*m4**2*sdip + m1**2*(m2**2 + 3.*m4**2 - 1.*tdip) - 3.*m4**2*tdip + tdip**2 + m2**2*(m4**2 - 2.*tdip - 1.*udip) - 1.*m4**2*udip + 
     &               tdip*udip)*BR(xtype,ftype) + (-1.*m4**4 + m4**2*sdip + m2**2*(-1.*m3**2 - 3.*m4**2 + sdip) + m4**2*tdip - 1.*sdip*tdip + m3**2*(-3.*m4**2 + tdip) + 
     &               2.*m4**2*udip)*BL(xtype,ftype)*Mf(pstype,psgen)) - 
     &         1.*dconjg(AL(chimeigen,ptmeigen + sumTc,pttype,ptgen))*((m2**4*(-2.*m3**2 + sdip + tdip - 1.*udip) + 
     &               m1**2*(-1.*m4**2*sdip + m3**2*(6.*m4**2 - 2.*tdip) + m4**2*tdip + sdip*tdip - 1.*tdip**2 + m2**2*(2.*m3**2 + 8.*m4**2 - 1.*sdip + tdip - 1.*udip) - 
     &                  5.*m4**2*udip + tdip*udip) + (m4**2 - 1.*tdip)*(-1.*sdip**2 - 1.*tdip**2 + m3**2*(sdip + tdip - 1.*udip) + m4**2*(sdip + tdip - 1.*udip) + udip**2) + 
     &               m2**2*(-1.*sdip**2 - 1.*sdip*tdip - 2.*tdip**2 + 2.*m4**2*(sdip - 1.*tdip - 1.*udip) + m3**2*(2.*m4**2 + sdip + 3.*tdip - 1.*udip) + tdip*udip + udip**2))*
     &             BR(xtype,ftype) + 2.*m3*(-1.*m4**4 + 2.*m4**2*sdip + m4**2*tdip + m1**2*(-1.*m2**2 - 3.*m4**2 + tdip) + m4**2*udip - 1.*tdip*udip + m2**2*(-3.*m4**2 + udip))*
     &             BL(xtype,ftype)*Mf(pstype,psgen))) + AR(chimeigen,sfmeigen,sftype,sfgen)*
     &       (2.*m1*dconjg(AL(chimeigen,ptmeigen + sumTc,pttype,ptgen))*(m3*
     &             (m2**4 + 2.*m4**4 - 2.*m4**2*sdip + m1**2*(m2**2 + 3.*m4**2 - 1.*tdip) - 3.*m4**2*tdip + tdip**2 + m2**2*(m4**2 - 2.*tdip - 1.*udip) - 1.*m4**2*udip + 
     &               tdip*udip)*BL(xtype,ftype) + (-1.*m4**4 + m4**2*sdip + m2**2*(-1.*m3**2 - 3.*m4**2 + sdip) + m4**2*tdip - 1.*sdip*tdip + m3**2*(-3.*m4**2 + tdip) + 
     &               2.*m4**2*udip)*BR(xtype,ftype)*Mf(pstype,psgen)) - 
     &         1.*dconjg(AR(chimeigen,ptmeigen + sumTc,pttype,ptgen))*((m2**4*(-2.*m3**2 + sdip + tdip - 1.*udip) + 
     &               m1**2*(-1.*m4**2*sdip + m3**2*(6.*m4**2 - 2.*tdip) + m4**2*tdip + sdip*tdip - 1.*tdip**2 + m2**2*(2.*m3**2 + 8.*m4**2 - 1.*sdip + tdip - 1.*udip) - 
     &                  5.*m4**2*udip + tdip*udip) + (m4**2 - 1.*tdip)*(-1.*sdip**2 - 1.*tdip**2 + m3**2*(sdip + tdip - 1.*udip) + m4**2*(sdip + tdip - 1.*udip) + udip**2) + 
     &               m2**2*(-1.*sdip**2 - 1.*sdip*tdip - 2.*tdip**2 + 2.*m4**2*(sdip - 1.*tdip - 1.*udip) + m3**2*(2.*m4**2 + sdip + 3.*tdip - 1.*udip) + tdip*udip + udip**2))*
     &             BL(xtype,ftype) + 2.*m3*(-1.*m4**4 + 2.*m4**2*sdip + m4**2*tdip + m1**2*(-1.*m2**2 - 3.*m4**2 + tdip) + m4**2*udip - 1.*tdip*udip + m2**2*(-3.*m4**2 + udip))*
     &             BR(xtype,ftype)*Mf(pstype,psgen)))))/(m4**2*(sdip - 1.*Mf(pstype,psgen)**2)*(tdip - 1.*MSf(ptmeigen + sumTc,pttype,ptgen)**2))

   	
	NeuQ2qx_M2TreeWZ_dipole = NeuQ2qx_M2TreeWZ_dipole + treeST

        enddo



! --- Mts  ---

        do sumT = 0, 1

        treeTS = (-0.5*colorfactor*NC*CL(ptmeigen + sumT,sfmeigen,sftype,sfgen,xtype)*
     &    (dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(-2.*m1*AR(chimeigen,ptmeigen + sumT,pttype,ptgen)*
     &          (m3*(m2**4 + 2.*m4**4 - 2.*m4**2*sdip + m1**2*(m2**2 + 3.*m4**2 - 1.*tdip) - 3.*m4**2*tdip + tdip**2 + m2**2*(m4**2 - 2.*tdip - 1.*udip) - 1.*m4**2*udip + 
     &               tdip*udip)*dconjg(BR(xtype,ftype)) + (-1.*m4**4 + m4**2*sdip + m2**2*(-1.*m3**2 - 3.*m4**2 + sdip) + m4**2*tdip - 1.*sdip*tdip + m3**2*(-3.*m4**2 + tdip) + 
     &               2.*m4**2*udip)*dconjg(BL(xtype,ftype))*Mf(pstype,psgen)) + 
     &         AL(chimeigen,ptmeigen + sumT,pttype,ptgen)*((m2**4*(-2.*m3**2 + sdip + tdip - 1.*udip) + 
     &               m1**2*(-1.*m4**2*sdip + m3**2*(6.*m4**2 - 2.*tdip) + m4**2*tdip + sdip*tdip - 1.*tdip**2 + m2**2*(2.*m3**2 + 8.*m4**2 - 1.*sdip + tdip - 1.*udip) - 
     &                  5.*m4**2*udip + tdip*udip) + (m4**2 - 1.*tdip)*(-1.*sdip**2 - 1.*tdip**2 + m3**2*(sdip + tdip - 1.*udip) + m4**2*(sdip + tdip - 1.*udip) + udip**2) + 
     &               m2**2*(-1.*sdip**2 - 1.*sdip*tdip - 2.*tdip**2 + 2.*m4**2*(sdip - 1.*tdip - 1.*udip) + m3**2*(2.*m4**2 + sdip + 3.*tdip - 1.*udip) + tdip*udip + udip**2))*
     &             dconjg(BR(xtype,ftype)) + 2.*m3*(-1.*m4**4 + 2.*m4**2*sdip + m4**2*tdip + m1**2*(-1.*m2**2 - 3.*m4**2 + tdip) + m4**2*udip - 1.*tdip*udip + 
     &               m2**2*(-3.*m4**2 + udip))*dconjg(BL(xtype,ftype))*Mf(pstype,psgen))) + 
     &      dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*(-2.*m1*AL(chimeigen,ptmeigen + sumT,pttype,ptgen)*
     &          (m3*(m2**4 + 2.*m4**4 - 2.*m4**2*sdip + m1**2*(m2**2 + 3.*m4**2 - 1.*tdip) - 3.*m4**2*tdip + tdip**2 + m2**2*(m4**2 - 2.*tdip - 1.*udip) - 1.*m4**2*udip + 
     &               tdip*udip)*dconjg(BL(xtype,ftype)) + (-1.*m4**4 + m4**2*sdip + m2**2*(-1.*m3**2 - 3.*m4**2 + sdip) + m4**2*tdip - 1.*sdip*tdip + m3**2*(-3.*m4**2 + tdip) + 
     &               2.*m4**2*udip)*dconjg(BR(xtype,ftype))*Mf(pstype,psgen)) + 
     &         AR(chimeigen,ptmeigen + sumT,pttype,ptgen)*((m2**4*(-2.*m3**2 + sdip + tdip - 1.*udip) + 
     &               m1**2*(-1.*m4**2*sdip + m3**2*(6.*m4**2 - 2.*tdip) + m4**2*tdip + sdip*tdip - 1.*tdip**2 + m2**2*(2.*m3**2 + 8.*m4**2 - 1.*sdip + tdip - 1.*udip) - 
     &                  5.*m4**2*udip + tdip*udip) + (m4**2 - 1.*tdip)*(-1.*sdip**2 - 1.*tdip**2 + m3**2*(sdip + tdip - 1.*udip) + m4**2*(sdip + tdip - 1.*udip) + udip**2) + 
     &               m2**2*(-1.*sdip**2 - 1.*sdip*tdip - 2.*tdip**2 + 2.*m4**2*(sdip - 1.*tdip - 1.*udip) + m3**2*(2.*m4**2 + sdip + 3.*tdip - 1.*udip) + tdip*udip + udip**2))*
     &             dconjg(BL(xtype,ftype)) + 2.*m3*(-1.*m4**4 + 2.*m4**2*sdip + m4**2*tdip + m1**2*(-1.*m2**2 - 3.*m4**2 + tdip) + m4**2*udip - 1.*tdip*udip + 
     &               m2**2*(-3.*m4**2 + udip))*dconjg(BR(xtype,ftype))*Mf(pstype,psgen)))))/
     &  (m4**2*(sdip - 1.*Mf(pstype,psgen)**2)*(tdip - 1.*MSf(ptmeigen + sumT,pttype,ptgen)**2))

  
       NeuQ2qx_M2TreeWZ_dipole = NeuQ2qx_M2TreeWZ_dipole + treeTS

        enddo


! --- Msu  ---

        do sumUc = 0, k

        treeSU =  NC*(0.5*(AL(chimeigen,sfmeigen,sftype,sfgen)*(2.*dconjg(AR(pumeigen + sumUc,sfmeigen,sftype,sfgen))*
     &          ((m3*BR(xtype,ftype)*(m1*(m1**4 - 4.*m4**4 + m4**2*sdip + 2.*m4**2*tdip + m1**2*(m2**2 + m4**2 - 1.*tdip - 2.*udip) + tdip*udip + udip**2 - 
     &                    1.*m2**2*(3.*m4**2 + udip))*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &                 3.*m4**2*(m1**2 - 1.*m2**2 + sdip)*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)))/m4**2 + 
     &            BL(xtype,ftype)*Mf(pstype,psgen)*(-3.*m1*(m1**2 + 2.*m3**2 + m4**2 - 1.*sdip - 1.*tdip)*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &               (1.*(m4**4 + m1**2*(m3**2 - 1.*sdip) + (-1.*m3**2 + sdip)*udip - 1.*m4**2*(sdip - 1.*tdip + udip))*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*
     &                  MNeuCha(pumeigen + sumUc))/m4**2)) + (dconjg(AL(pumeigen + sumUc,sfmeigen,sftype,sfgen))*
     &            (2.*m3*BL(xtype,ftype)*Mf(pstype,psgen)*((m1**4 - 2.*m4**4 + m1**2*(m4**2 - 2.*udip) + m4**2*udip + udip**2)*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &                 6.*m1*m4**2*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)) + 
     &              BR(xtype,ftype)*((-8.*m3**2*m4**2*sdip - 4.*m4**4*sdip + 4.*m4**2*sdip**2 + 4.*m4**2*sdip*tdip - 1.*m1**4*(sdip + tdip - 1.*udip) + m3**2*sdip*udip + 
     &                    m4**2*sdip*udip - 1.*sdip**2*udip + m3**2*tdip*udip + m4**2*tdip*udip - 1.*tdip**2*udip - 1.*m3**2*udip**2 - 1.*m4**2*udip**2 + udip**3 + 
     &                    m2**2*(4.*m4**4 - 4.*m4**2*(sdip + tdip) + m3**2*(8.*m4**2 - 2.*udip) + (sdip + tdip - 1.*udip)*udip) + 
     &                    m1**2*(-5.*m4**2*sdip + sdip**2 - 1.*m4**2*tdip + 2.*sdip*tdip + tdip**2 - 3.*m4**2*udip + sdip*udip + tdip*udip - 2.*udip**2 + 
     &                       m3**2*(-1.*sdip - 3.*tdip + udip) + m2**2*(4.*m3**2 + 8.*m4**2 - 3.*sdip - 1.*tdip + udip)))*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &                 2.*m1*(2.*m4**4 + m1**2*(m3**2 - 1.*sdip) + m2**2*(m3**2 - 1.*sdip) - 2.*m4**2*sdip - 1.*m3**2*tdip + sdip*tdip - 1.*m3**2*udip + sdip*udip)*
     &                  dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc))))/m4**2) + 
     &      AR(chimeigen,sfmeigen,sftype,sfgen)*((dconjg(AR(pumeigen + sumUc,sfmeigen,sftype,sfgen))*
     &            (2.*m3*BR(xtype,ftype)*Mf(pstype,psgen)*((m1**4 - 2.*m4**4 + m1**2*(m4**2 - 2.*udip) + m4**2*udip + udip**2)*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &                 6.*m1*m4**2*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)) + 
     &              BL(xtype,ftype)*((-8.*m3**2*m4**2*sdip - 4.*m4**4*sdip + 4.*m4**2*sdip**2 + 4.*m4**2*sdip*tdip - 1.*m1**4*(sdip + tdip - 1.*udip) + m3**2*sdip*udip + 
     &                    m4**2*sdip*udip - 1.*sdip**2*udip + m3**2*tdip*udip + m4**2*tdip*udip - 1.*tdip**2*udip - 1.*m3**2*udip**2 - 1.*m4**2*udip**2 + udip**3 + 
     &                    m2**2*(4.*m4**4 - 4.*m4**2*(sdip + tdip) + m3**2*(8.*m4**2 - 2.*udip) + (sdip + tdip - 1.*udip)*udip) + 
     &                    m1**2*(-5.*m4**2*sdip + sdip**2 - 1.*m4**2*tdip + 2.*sdip*tdip + tdip**2 - 3.*m4**2*udip + sdip*udip + tdip*udip - 2.*udip**2 + 
     &                       m3**2*(-1.*sdip - 3.*tdip + udip) + m2**2*(4.*m3**2 + 8.*m4**2 - 3.*sdip - 1.*tdip + udip)))*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &                 2.*m1*(2.*m4**4 + m1**2*(m3**2 - 1.*sdip) + m2**2*(m3**2 - 1.*sdip) - 2.*m4**2*sdip - 1.*m3**2*tdip + sdip*tdip - 1.*m3**2*udip + sdip*udip)*
     &                  dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc))))/m4**2 + 
     &         2.*dconjg(AL(pumeigen + sumUc,sfmeigen,sftype,sfgen))*((m3*BL(xtype,ftype)*
     &               (m1*(m1**4 - 4.*m4**4 + m4**2*sdip + 2.*m4**2*tdip + m1**2*(m2**2 + m4**2 - 1.*tdip - 2.*udip) + tdip*udip + udip**2 - 1.*m2**2*(3.*m4**2 + udip))*
     &                  dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)) - 3.*m4**2*(m1**2 - 1.*m2**2 + sdip)*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)))/
     &             m4**2 + BR(xtype,ftype)*Mf(pstype,psgen)*(-3.*m1*(m1**2 + 2.*m3**2 + m4**2 - 1.*sdip - 1.*tdip)*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)) - 
     &               (1.*(m4**4 + m1**2*(m3**2 - 1.*sdip) + (-1.*m3**2 + sdip)*udip - 1.*m4**2*(sdip - 1.*tdip + udip))*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*
     &                  MNeuCha(pumeigen + sumUc))/m4**2)))))/((sdip - 1.*Mf(pstype,psgen)**2)*(udip - 1.*MNeuCha(pumeigen + sumUc)**2))
  
	NeuQ2qx_M2TreeWZ_dipole = NeuQ2qx_M2TreeWZ_dipole + treeSU
     
        enddo



! --- Mus  ---

  	do sumU = 0, k     

        treeUS = NC*(0.5*(dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))*(2.*AR(pumeigen + sumU,sfmeigen,sftype,sfgen)*
     &          ((m3*dconjg(BR(xtype,ftype))*(m1*(m1**4 - 4.*m4**4 + m4**2*sdip + 2.*m4**2*tdip + m1**2*(m2**2 + m4**2 - 1.*tdip - 2.*udip) + tdip*udip + udip**2 - 
     &                    1.*m2**2*(3.*m4**2 + udip))*DR(chimeigen,pumeigen + sumU,sftype,xtype) - 
     &                 3.*m4**2*(m1**2 - 1.*m2**2 + sdip)*DL(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)))/m4**2 + 
     &            dconjg(BL(xtype,ftype))*Mf(pstype,psgen)*(-3.*m1*(m1**2 + 2.*m3**2 + m4**2 - 1.*sdip - 1.*tdip)*DR(chimeigen,pumeigen + sumU,sftype,xtype) - 
     &               (1.*(m4**4 + m1**2*(m3**2 - 1.*sdip) + (-1.*m3**2 + sdip)*udip - 1.*m4**2*(sdip - 1.*tdip + udip))*DL(chimeigen,pumeigen + sumU,sftype,xtype)*
     &                  MNeuCha(pumeigen + sumU))/m4**2)) + (AL(pumeigen + sumU,sfmeigen,sftype,sfgen)*
     &            (2.*m3*dconjg(BL(xtype,ftype))*Mf(pstype,psgen)*((m1**4 - 2.*m4**4 + m1**2*(m4**2 - 2.*udip) + m4**2*udip + udip**2)*DL(chimeigen,pumeigen + sumU,sftype,xtype) - 
     &                 6.*m1*m4**2*DR(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)) + 
     &              dconjg(BR(xtype,ftype))*((-8.*m3**2*m4**2*sdip - 4.*m4**4*sdip + 4.*m4**2*sdip**2 + 4.*m4**2*sdip*tdip - 1.*m1**4*(sdip + tdip - 1.*udip) + m3**2*sdip*udip + 
     &                    m4**2*sdip*udip - 1.*sdip**2*udip + m3**2*tdip*udip + m4**2*tdip*udip - 1.*tdip**2*udip - 1.*m3**2*udip**2 - 1.*m4**2*udip**2 + udip**3 + 
     &                    m2**2*(4.*m4**4 - 4.*m4**2*(sdip + tdip) + m3**2*(8.*m4**2 - 2.*udip) + (sdip + tdip - 1.*udip)*udip) + 
     &                    m1**2*(-5.*m4**2*sdip + sdip**2 - 1.*m4**2*tdip + 2.*sdip*tdip + tdip**2 - 3.*m4**2*udip + sdip*udip + tdip*udip - 2.*udip**2 + 
     &                       m3**2*(-1.*sdip - 3.*tdip + udip) + m2**2*(4.*m3**2 + 8.*m4**2 - 3.*sdip - 1.*tdip + udip)))*DL(chimeigen,pumeigen + sumU,sftype,xtype) - 
     &                 2.*m1*(2.*m4**4 + m1**2*(m3**2 - 1.*sdip) + m2**2*(m3**2 - 1.*sdip) - 2.*m4**2*sdip - 1.*m3**2*tdip + sdip*tdip - 1.*m3**2*udip + sdip*udip)*
     &                  DR(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU))))/m4**2) + 
     &      dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))*((AR(pumeigen + sumU,sfmeigen,sftype,sfgen)*
     &            (2.*m3*dconjg(BR(xtype,ftype))*Mf(pstype,psgen)*((m1**4 - 2.*m4**4 + m1**2*(m4**2 - 2.*udip) + m4**2*udip + udip**2)*DR(chimeigen,pumeigen + sumU,sftype,xtype) - 
     &                 6.*m1*m4**2*DL(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)) + 
     &              dconjg(BL(xtype,ftype))*((-8.*m3**2*m4**2*sdip - 4.*m4**4*sdip + 4.*m4**2*sdip**2 + 4.*m4**2*sdip*tdip - 1.*m1**4*(sdip + tdip - 1.*udip) + m3**2*sdip*udip + 
     &                    m4**2*sdip*udip - 1.*sdip**2*udip + m3**2*tdip*udip + m4**2*tdip*udip - 1.*tdip**2*udip - 1.*m3**2*udip**2 - 1.*m4**2*udip**2 + udip**3 + 
     &                    m2**2*(4.*m4**4 - 4.*m4**2*(sdip + tdip) + m3**2*(8.*m4**2 - 2.*udip) + (sdip + tdip - 1.*udip)*udip) + 
     &                    m1**2*(-5.*m4**2*sdip + sdip**2 - 1.*m4**2*tdip + 2.*sdip*tdip + tdip**2 - 3.*m4**2*udip + sdip*udip + tdip*udip - 2.*udip**2 + 
     &                       m3**2*(-1.*sdip - 3.*tdip + udip) + m2**2*(4.*m3**2 + 8.*m4**2 - 3.*sdip - 1.*tdip + udip)))*DR(chimeigen,pumeigen + sumU,sftype,xtype) - 
     &                 2.*m1*(2.*m4**4 + m1**2*(m3**2 - 1.*sdip) + m2**2*(m3**2 - 1.*sdip) - 2.*m4**2*sdip - 1.*m3**2*tdip + sdip*tdip - 1.*m3**2*udip + sdip*udip)*
     &                  DL(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU))))/m4**2 + 
     &         2.*AL(pumeigen + sumU,sfmeigen,sftype,sfgen)*((m3*dconjg(BL(xtype,ftype))*
     &               (m1*(m1**4 - 4.*m4**4 + m4**2*sdip + 2.*m4**2*tdip + m1**2*(m2**2 + m4**2 - 1.*tdip - 2.*udip) + tdip*udip + udip**2 - 1.*m2**2*(3.*m4**2 + udip))*
     &                  DL(chimeigen,pumeigen + sumU,sftype,xtype) - 3.*m4**2*(m1**2 - 1.*m2**2 + sdip)*DR(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)))/m4**2 + 
     &            dconjg(BR(xtype,ftype))*Mf(pstype,psgen)*(-3.*m1*(m1**2 + 2.*m3**2 + m4**2 - 1.*sdip - 1.*tdip)*DL(chimeigen,pumeigen + sumU,sftype,xtype) - 
     &               (1.*(m4**4 + m1**2*(m3**2 - 1.*sdip) + (-1.*m3**2 + sdip)*udip - 1.*m4**2*(sdip - 1.*tdip + udip))*DR(chimeigen,pumeigen + sumU,sftype,xtype)*
     &                  MNeuCha(pumeigen + sumU))/m4**2)))))/((sdip - 1.*Mf(pstype,psgen)**2)*(udip - 1.*MNeuCha(pumeigen + sumU)**2))

        NeuQ2qx_M2TreeWZ_dipole = NeuQ2qx_M2TreeWZ_dipole + treeUS

        enddo


! --- Mtu  ---

        do sumT = 0, 1
          do sumUc = 0, k

        treeTU =  NC*(CL(ptmeigen + sumT,sfmeigen,sftype,sfgen,xtype)*(AR(chimeigen,ptmeigen + sumT,pttype,ptgen)*
     &       (dconjg(AR(pumeigen + sumUc,sfmeigen,sftype,sfgen))*(((-3.*m3**2*m4**2*sdip - 1.*m4**4*sdip + m4**2*sdip**2 + m1**4*(m2**2 + 3.*m4**2 - 1.*tdip) - 
     &                 1.*m3**2*m4**2*tdip - 1.*m4**4*tdip + 2.*m4**2*sdip*tdip + m4**2*tdip**2 + m4**4*udip + m3**2*tdip*udip + m4**2*tdip*udip - 1.*tdip**2*udip - 
     &                 1.*m4**2*udip**2 + m1**2*(3.*m4**4 - 4.*m4**2*sdip + m3**2*(6.*m4**2 - 2.*tdip) - 5.*m4**2*tdip + sdip*tdip + tdip**2 + 
     &                    m2**2*(2.*m3**2 + 5.*m4**2 - 1.*sdip - 1.*tdip - 1.*udip) - 2.*m4**2*udip + tdip*udip) + 
     &                 m2**2*(m3**2*(4.*m4**2 - 1.*udip) + tdip*udip + m4**2*(-1.*sdip - 3.*tdip + udip)))*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)))/m4**2 + 
     &            m1*(-1.*m4**2 - (1.*m2**2*(m3**2 + 3.*m4**2 - 1.*sdip))/m4**2 + sdip + tdip - (1.*sdip*tdip)/m4**2 + m3**2*(-3. + tdip/m4**2) + 2.*udip)*
     &             dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)) + 
     &         (m3*(m4**4 - 2.*m4**2*sdip + m1**2*(m2**2 + 3.*m4**2 - 1.*tdip) - 1.*m4**2*tdip + m2**2*(3.*m4**2 - 1.*udip) - 1.*m4**2*udip + tdip*udip)*
     &            dconjg(AL(pumeigen + sumUc,sfmeigen,sftype,sfgen))*(m1*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)) + 
     &              dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)))/m4**2) + 
     &      AL(chimeigen,ptmeigen + sumT,pttype,ptgen)*((m3*(m4**4 - 2.*m4**2*sdip + m1**2*(m2**2 + 3.*m4**2 - 1.*tdip) - 1.*m4**2*tdip + m2**2*(3.*m4**2 - 1.*udip) - 
     &              1.*m4**2*udip + tdip*udip)*dconjg(AR(pumeigen + sumUc,sfmeigen,sftype,sfgen))*
     &            (m1*dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype)) + dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)))/m4**2 + 
     &         dconjg(AL(pumeigen + sumUc,sfmeigen,sftype,sfgen))*(((-3.*m3**2*m4**2*sdip - 1.*m4**4*sdip + m4**2*sdip**2 + m1**4*(m2**2 + 3.*m4**2 - 1.*tdip) - 
     &                 1.*m3**2*m4**2*tdip - 1.*m4**4*tdip + 2.*m4**2*sdip*tdip + m4**2*tdip**2 + m4**4*udip + m3**2*tdip*udip + m4**2*tdip*udip - 1.*tdip**2*udip - 
     &                 1.*m4**2*udip**2 + m1**2*(3.*m4**4 - 4.*m4**2*sdip + m3**2*(6.*m4**2 - 2.*tdip) - 5.*m4**2*tdip + sdip*tdip + tdip**2 + 
     &                    m2**2*(2.*m3**2 + 5.*m4**2 - 1.*sdip - 1.*tdip - 1.*udip) - 2.*m4**2*udip + tdip*udip) + 
     &                 m2**2*(m3**2*(4.*m4**2 - 1.*udip) + tdip*udip + m4**2*(-1.*sdip - 3.*tdip + udip)))*dconjg(DL(chimeigen,pumeigen + sumUc,sftype,xtype)))/m4**2 + 
     &            m1*(-1.*m4**2 - (1.*m2**2*(m3**2 + 3.*m4**2 - 1.*sdip))/m4**2 + sdip + tdip - (1.*sdip*tdip)/m4**2 + m3**2*(-3. + tdip/m4**2) + 2.*udip)*
     &             dconjg(DR(chimeigen,pumeigen + sumUc,sftype,xtype))*MNeuCha(pumeigen + sumUc)))))/
     &  ((udip - 1.*MNeuCha(pumeigen + sumUc)**2)*(tdip - 1.*MSf(ptmeigen + sumT,pttype,ptgen)**2))
   
	NeuQ2qx_M2TreeWZ_dipole = NeuQ2qx_M2TreeWZ_dipole + treeTU

          enddo
        enddo



! --- Mut  ---

        do sumTc = 0, 1
          do sumU = 0, k     

        treeUT = NC*(dconjg(CL(ptmeigen + sumTc,sfmeigen,sftype,sfgen,xtype))*(dconjg(AR(chimeigen,ptmeigen + sumTc,pttype,ptgen))*
     &       (AR(pumeigen + sumU,sfmeigen,sftype,sfgen)*(((-3.*m3**2*m4**2*sdip - 1.*m4**4*sdip + m4**2*sdip**2 + m1**4*(m2**2 + 3.*m4**2 - 1.*tdip) - 1.*m3**2*m4**2*tdip - 
     &                 1.*m4**4*tdip + 2.*m4**2*sdip*tdip + m4**2*tdip**2 + m4**4*udip + m3**2*tdip*udip + m4**2*tdip*udip - 1.*tdip**2*udip - 1.*m4**2*udip**2 + 
     &                 m1**2*(3.*m4**4 - 4.*m4**2*sdip + m3**2*(6.*m4**2 - 2.*tdip) - 5.*m4**2*tdip + sdip*tdip + tdip**2 + 
     &                    m2**2*(2.*m3**2 + 5.*m4**2 - 1.*sdip - 1.*tdip - 1.*udip) - 2.*m4**2*udip + tdip*udip) + 
     &                 m2**2*(m3**2*(4.*m4**2 - 1.*udip) + tdip*udip + m4**2*(-1.*sdip - 3.*tdip + udip)))*DR(chimeigen,pumeigen + sumU,sftype,xtype))/m4**2 + 
     &            m1*(-1.*m4**2 - (1.*m2**2*(m3**2 + 3.*m4**2 - 1.*sdip))/m4**2 + sdip + tdip - (1.*sdip*tdip)/m4**2 + m3**2*(-3. + tdip/m4**2) + 2.*udip)*
     &             DL(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)) + 
     &         (m3*(m4**4 - 2.*m4**2*sdip + m1**2*(m2**2 + 3.*m4**2 - 1.*tdip) - 1.*m4**2*tdip + m2**2*(3.*m4**2 - 1.*udip) - 1.*m4**2*udip + tdip*udip)*
     &            AL(pumeigen + sumU,sfmeigen,sftype,sfgen)*(m1*DL(chimeigen,pumeigen + sumU,sftype,xtype) + DR(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)))/m4**2) + 
     &      dconjg(AL(chimeigen,ptmeigen + sumTc,pttype,ptgen))*((m3*(m4**4 - 2.*m4**2*sdip + m1**2*(m2**2 + 3.*m4**2 - 1.*tdip) - 1.*m4**2*tdip + m2**2*(3.*m4**2 - 1.*udip) - 
     &              1.*m4**2*udip + tdip*udip)*AR(pumeigen + sumU,sfmeigen,sftype,sfgen)*
     &            (m1*DR(chimeigen,pumeigen + sumU,sftype,xtype) + DL(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)))/m4**2 + 
     &         AL(pumeigen + sumU,sfmeigen,sftype,sfgen)*(((-3.*m3**2*m4**2*sdip - 1.*m4**4*sdip + m4**2*sdip**2 + m1**4*(m2**2 + 3.*m4**2 - 1.*tdip) - 1.*m3**2*m4**2*tdip - 
     &                 1.*m4**4*tdip + 2.*m4**2*sdip*tdip + m4**2*tdip**2 + m4**4*udip + m3**2*tdip*udip + m4**2*tdip*udip - 1.*tdip**2*udip - 1.*m4**2*udip**2 + 
     &                 m1**2*(3.*m4**4 - 4.*m4**2*sdip + m3**2*(6.*m4**2 - 2.*tdip) - 5.*m4**2*tdip + sdip*tdip + tdip**2 + 
     &                    m2**2*(2.*m3**2 + 5.*m4**2 - 1.*sdip - 1.*tdip - 1.*udip) - 2.*m4**2*udip + tdip*udip) + 
     &                 m2**2*(m3**2*(4.*m4**2 - 1.*udip) + tdip*udip + m4**2*(-1.*sdip - 3.*tdip + udip)))*DL(chimeigen,pumeigen + sumU,sftype,xtype))/m4**2 + 
     &            m1*(-1.*m4**2 - (1.*m2**2*(m3**2 + 3.*m4**2 - 1.*sdip))/m4**2 + sdip + tdip - (1.*sdip*tdip)/m4**2 + m3**2*(-3. + tdip/m4**2) + 2.*udip)*
     &             DR(chimeigen,pumeigen + sumU,sftype,xtype)*MNeuCha(pumeigen + sumU)))))/((udip - 1.*MNeuCha(pumeigen + sumU)**2)*(tdip - 1.*MSf(ptmeigen + sumTc,pttype,ptgen)**2))

   
        NeuQ2qx_M2TreeWZ_dipole = NeuQ2qx_M2TreeWZ_dipole + treeUT

          enddo
        enddo

  
       end



      double precision function NeuQ2qx_M2_FEFS_gg()

      implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

      double precision dipcolor,muQ2,zi,zj,A1L,A1R,A1Lc,A1Rc,ab1,ab2, ab3, abb5
      double precision kappa,yijk,vijk,viji,zP,zM,vijkTilde,zim,zjm,Q2

      Q2 = sman
      muQ2 = m3**2/Q2


      A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
      A1R = AR(chimeigen,sfmeigen,sftype,sfgen)
      A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
      A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))

      yijk = p2p3/(p1p2 + p1p3 + p2p3) 
      vijk = dsqrt((2*muQ2+(1-muQ2)*(1-yijk))**2-4*muQ2)/((1-muQ2)*(1-yijk)) 
      viji = 1d0 
      zP = (1+viji*vijk)/2d0 
      zM = (1-viji*vijk)/2d0 

      zi = p1p3/(p1p2 + p1p3)
      zj = 1 - zi

      zim = zi - zM
      zjm = zj - zM

      ab1 = kappa(Q2, 0d0, m3**2)/kappa(Q2, 2*p2p3, m3**2)
      ab2 = (pap1+pbp1)/sman
      ab3 = (Q2 + m3**2)/(2*Q2)
      abb5 = 1/(1 - (1 - yijk)*zi) + 1/(1 - (1 - yijk)*zj) + (-2 + kappaDip*zM*zP)/vijk

      dipcolor =  -1/2d0

! Ms^2
      NeuQ2qx_M2_FEFS_gg = dipcolor/(2*p2p3)*((4*CV*gs**4*(A1Lc*A1R*abb5*m1*m3**3*p2p3*vijk + A1L*A1Rc*abb5*m1*m3**3*p2p3*vijk - A1Lc*A1R*abb5*m1*m3**3*Nc**2*p2p3*vijk - A1L*A1Rc*abb5*m1*m3**3*Nc**2*p2p3*vijk + A1L*A1Lc*abb5*m3**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*vijk + 
     -      A1R*A1Rc*abb5*m3**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*vijk - A1L*A1Lc*abb5*m3**2*Nc**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*vijk - 
     -      A1R*A1Rc*abb5*m3**2*Nc**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*vijk + A1L*A1Lc*abb5*m3**2*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*vijk + 
     -      A1R*A1Rc*abb5*m3**2*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*vijk - A1L*A1Lc*abb5*m3**2*Nc**2*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*vijk - A1R*A1Rc*abb5*m3**2*Nc**2*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*vijk - 
     -      A1Lc*A1R*abb5*m1*m3*p2p3*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*vijk - A1L*A1Rc*abb5*m1*m3*p2p3*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*vijk + 
     -      A1Lc*A1R*abb5*m1*m3*Nc**2*p2p3*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*vijk + A1L*A1Rc*abb5*m1*m3*Nc**2*p2p3*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*vijk + 
     -      A1L*A1Lc*abb5*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*vijk + 
     -      A1R*A1Rc*abb5*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*vijk - 
     -      A1L*A1Lc*abb5*Nc**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*vijk - 
     -      A1R*A1Rc*abb5*Nc**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*vijk - 
     -      A1L*A1Lc*pap3*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      A1R*A1Rc*pap3*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      A1L*A1Lc*Nc**2*pap3*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      A1R*A1Rc*Nc**2*pap3*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      A1Lc*A1R*m1*m3*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      A1L*A1Rc*m1*m3*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      A1Lc*A1R*m1*m3*Nc**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      A1L*A1Rc*m1*m3*Nc**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      A1Lc*A1R*m1*m3*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 - A1L*A1Rc*m1*m3*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 + A1Lc*A1R*m1*m3*Nc**2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 + 
     -      A1L*A1Rc*m1*m3*Nc**2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 - A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 - 
     -      A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 + 
     -      A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 + 
     -      A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 - 
     -      A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 - A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 + 
     -      A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 + A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 + 
     -      A1L*A1Lc*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*zim*zjm + 
     -      A1R*A1Rc*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*zim*zjm - 
     -      A1L*A1Lc*Nc**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*zim*zjm - 
     -      A1R*A1Rc*Nc**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*zim*zjm + 
     -      A1L*A1Lc*pap3*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zim*zjm + 
     -      A1R*A1Rc*pap3*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zim*zjm - 
     -      A1L*A1Lc*Nc**2*pap3*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zim*zjm - 
     -      A1R*A1Rc*Nc**2*pap3*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zim*zjm + 
     -      A1Lc*A1R*m1*m3*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      A1L*A1Rc*m1*m3*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      A1Lc*A1R*m1*m3*Nc**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      A1L*A1Rc*m1*m3*Nc**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      A1L*A1Lc*pap2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      A1R*A1Rc*pap2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      A1L*A1Lc*Nc**2*pap2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      A1R*A1Rc*Nc**2*pap2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      A1Lc*A1R*m1*m3*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      A1L*A1Rc*m1*m3*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      A1Lc*A1R*m1*m3*Nc**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      A1L*A1Rc*m1*m3*Nc**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      2*A1Lc*A1R*m1*m3*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 2*A1L*A1Rc*m1*m3*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      2*A1Lc*A1R*m1*m3*Nc**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 2*A1L*A1Rc*m1*m3*Nc**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      2*A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      2*A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      2*A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      2*A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      2*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      2*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      2*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      2*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      A1L*A1Lc*pap2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      A1R*A1Rc*pap2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      A1L*A1Lc*Nc**2*pap2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      A1R*A1Rc*Nc**2*pap2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      A1Lc*A1R*m1*m3*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      A1L*A1Rc*m1*m3*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      A1Lc*A1R*m1*m3*Nc**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      A1L*A1Rc*m1*m3*Nc**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      A1Lc*A1R*m1*m3*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 - A1L*A1Rc*m1*m3*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 + A1Lc*A1R*m1*m3*Nc**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 + 
     -      A1L*A1Rc*m1*m3*Nc**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 - A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 - 
     -      A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 + 
     -      A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 + 
     -      A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 - 
     -      A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 - A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 + 
     -      A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 + A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2))/
     -  (p2p3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*vijk))

! Mt2

      NeuQ2qx_M2_FEFS_gg = NeuQ2qx_M2_FEFS_gg + dipcolor/(2*p2p3)*((CV*gs**4*(-1 + Nc**2)*(-(A1Lc*A1R*m1*m3) - A1L*A1Rc*m1*m3 + A1L*A1Lc*(-(ab1*pap1) + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb)) + A1R*A1Rc*(-(ab1*pap1) + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb)))*
     -    (abb5*p2p3*((1 + ab1*ab2 - ab3)**2*m1**2 + (1 - ab1*ab2 + ab3)**2*m2**2 + ab1**2*m3**2 - 2*ab1*pap1 - 2*ab1**2*ab2*pap1 + 2*ab1*ab3*pap1 - 2*papb + 2*ab1**2*ab2**2*papb - 4*ab1*ab2*ab3*papb + 2*ab3**2*papb + 2*ab1*pbp1 - 2*ab1**2*ab2*pbp1 + 2*ab1*ab3*pbp1)*vijk - 
     -      ((-1 + ab3)*pap3*zim + pbp3*zim + ab3*pbp3*zim + pap2*zjm - ab3*pap2*zjm - pbp2*zjm - ab3*pbp2*zjm + ab1*(p1p3*zim - p1p2*zjm + ab2*(-(pap3*zim) - pbp3*zim + (pap2 + pbp2)*zjm)))**2))/(p2p3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*vijk))


! 2 MsMt
      NeuQ2qx_M2_FEFS_gg = NeuQ2qx_M2_FEFS_gg + 2*dipcolor/(2*p2p3)*((CV*gs**4*(-1 + Nc**2)*(2*A1L*A1Lc*ab1**2*ab2**2*abb5*m1**4*p2p3*vijk + 2*A1R*A1Rc*ab1**2*ab2**2*abb5*m1**4*p2p3*vijk - 4*A1L*A1Lc*ab1*ab2*ab3*abb5*m1**4*p2p3*vijk - 4*A1R*A1Rc*ab1*ab2*ab3*abb5*m1**4*p2p3*vijk + 2*A1L*A1Lc*ab3**2*abb5*m1**4*p2p3*vijk + 
     -      2*A1R*A1Rc*ab3**2*abb5*m1**4*p2p3*vijk + 4*A1L*A1Lc*ab1*ab2*abb5*m1**2*m2**2*p2p3*vijk + 4*A1R*A1Rc*ab1*ab2*abb5*m1**2*m2**2*p2p3*vijk + 2*A1L*A1Lc*ab1**2*ab2**2*abb5*m1**2*m2**2*p2p3*vijk + 2*A1R*A1Rc*ab1**2*ab2**2*abb5*m1**2*m2**2*p2p3*vijk - 
     -      4*A1L*A1Lc*ab3*abb5*m1**2*m2**2*p2p3*vijk - 4*A1R*A1Rc*ab3*abb5*m1**2*m2**2*p2p3*vijk - 4*A1L*A1Lc*ab1*ab2*ab3*abb5*m1**2*m2**2*p2p3*vijk - 4*A1R*A1Rc*ab1*ab2*ab3*abb5*m1**2*m2**2*p2p3*vijk + 2*A1L*A1Lc*ab3**2*abb5*m1**2*m2**2*p2p3*vijk + 
     -      2*A1R*A1Rc*ab3**2*abb5*m1**2*m2**2*p2p3*vijk - 2*A1Lc*A1R*ab1*ab2*abb5*m1**3*m3*p2p3*vijk - 2*A1L*A1Rc*ab1*ab2*abb5*m1**3*m3*p2p3*vijk + 2*A1Lc*A1R*ab3*abb5*m1**3*m3*p2p3*vijk + 2*A1L*A1Rc*ab3*abb5*m1**3*m3*p2p3*vijk - 2*A1Lc*A1R*abb5*m1*m2**2*m3*p2p3*vijk - 
     -      2*A1L*A1Rc*abb5*m1*m2**2*m3*p2p3*vijk + 2*A1L*A1Lc*ab1*ab2*abb5*m1**2*m3**2*p2p3*vijk + 2*A1R*A1Rc*ab1*ab2*abb5*m1**2*m3**2*p2p3*vijk - 2*A1L*A1Lc*ab3*abb5*m1**2*m3**2*p2p3*vijk - 2*A1R*A1Rc*ab3*abb5*m1**2*m3**2*p2p3*vijk - 2*A1Lc*A1R*abb5*m1*m3**3*p2p3*vijk - 
     -      2*A1L*A1Rc*abb5*m1*m3**3*p2p3*vijk - 4*A1L*A1Lc*ab1**2*ab2*abb5*m1**2*p2p3*pap1*vijk - 4*A1R*A1Rc*ab1**2*ab2*abb5*m1**2*p2p3*pap1*vijk + 4*A1L*A1Lc*ab1*ab3*abb5*m1**2*p2p3*pap1*vijk + 4*A1R*A1Rc*ab1*ab3*abb5*m1**2*p2p3*pap1*vijk - 2*A1L*A1Lc*ab1*abb5*m2**2*p2p3*pap1*vijk - 
     -      2*A1R*A1Rc*ab1*abb5*m2**2*p2p3*pap1*vijk - 2*A1L*A1Lc*ab1**2*ab2*abb5*m2**2*p2p3*pap1*vijk - 2*A1R*A1Rc*ab1**2*ab2*abb5*m2**2*p2p3*pap1*vijk + 2*A1L*A1Lc*ab1*ab3*abb5*m2**2*p2p3*pap1*vijk + 2*A1R*A1Rc*ab1*ab3*abb5*m2**2*p2p3*pap1*vijk + 2*A1Lc*A1R*ab1*abb5*m1*m3*p2p3*pap1*vijk + 
     -      2*A1L*A1Rc*ab1*abb5*m1*m3*p2p3*pap1*vijk - 2*A1L*A1Lc*ab1*abb5*m3**2*p2p3*pap1*vijk - 2*A1R*A1Rc*ab1*abb5*m3**2*p2p3*pap1*vijk + 2*A1L*A1Lc*ab1**2*abb5*p2p3*pap1**2*vijk + 2*A1R*A1Rc*ab1**2*abb5*p2p3*pap1**2*vijk + 2*A1L*A1Lc*ab1*ab2*abb5*m1**2*p2p3*papb*vijk + 
     -      2*A1R*A1Rc*ab1*ab2*abb5*m1**2*p2p3*papb*vijk + 6*A1L*A1Lc*ab1**2*ab2**2*abb5*m1**2*p2p3*papb*vijk + 6*A1R*A1Rc*ab1**2*ab2**2*abb5*m1**2*p2p3*papb*vijk - 2*A1L*A1Lc*ab3*abb5*m1**2*p2p3*papb*vijk - 2*A1R*A1Rc*ab3*abb5*m1**2*p2p3*papb*vijk - 
     -      12*A1L*A1Lc*ab1*ab2*ab3*abb5*m1**2*p2p3*papb*vijk - 12*A1R*A1Rc*ab1*ab2*ab3*abb5*m1**2*p2p3*papb*vijk + 6*A1L*A1Lc*ab3**2*abb5*m1**2*p2p3*papb*vijk + 6*A1R*A1Rc*ab3**2*abb5*m1**2*p2p3*papb*vijk + 2*A1L*A1Lc*ab1*ab2*abb5*m2**2*p2p3*papb*vijk + 
     -      2*A1R*A1Rc*ab1*ab2*abb5*m2**2*p2p3*papb*vijk + 2*A1L*A1Lc*ab1**2*ab2**2*abb5*m2**2*p2p3*papb*vijk + 2*A1R*A1Rc*ab1**2*ab2**2*abb5*m2**2*p2p3*papb*vijk - 2*A1L*A1Lc*ab3*abb5*m2**2*p2p3*papb*vijk - 2*A1R*A1Rc*ab3*abb5*m2**2*p2p3*papb*vijk - 
     -      4*A1L*A1Lc*ab1*ab2*ab3*abb5*m2**2*p2p3*papb*vijk - 4*A1R*A1Rc*ab1*ab2*ab3*abb5*m2**2*p2p3*papb*vijk + 2*A1L*A1Lc*ab3**2*abb5*m2**2*p2p3*papb*vijk + 2*A1R*A1Rc*ab3**2*abb5*m2**2*p2p3*papb*vijk - 2*A1Lc*A1R*abb5*m1*m3*p2p3*papb*vijk - 2*A1L*A1Rc*abb5*m1*m3*p2p3*papb*vijk - 
     -      2*A1Lc*A1R*ab1*ab2*abb5*m1*m3*p2p3*papb*vijk - 2*A1L*A1Rc*ab1*ab2*abb5*m1*m3*p2p3*papb*vijk + 2*A1Lc*A1R*ab3*abb5*m1*m3*p2p3*papb*vijk + 2*A1L*A1Rc*ab3*abb5*m1*m3*p2p3*papb*vijk - 2*A1L*A1Lc*abb5*m3**2*p2p3*papb*vijk - 2*A1R*A1Rc*abb5*m3**2*p2p3*papb*vijk + 
     -      2*A1L*A1Lc*ab1*ab2*abb5*m3**2*p2p3*papb*vijk + 2*A1R*A1Rc*ab1*ab2*abb5*m3**2*p2p3*papb*vijk - 2*A1L*A1Lc*ab3*abb5*m3**2*p2p3*papb*vijk - 2*A1R*A1Rc*ab3*abb5*m3**2*p2p3*papb*vijk - 6*A1L*A1Lc*ab1**2*ab2*abb5*p2p3*pap1*papb*vijk - 6*A1R*A1Rc*ab1**2*ab2*abb5*p2p3*pap1*papb*vijk + 
     -      6*A1L*A1Lc*ab1*ab3*abb5*p2p3*pap1*papb*vijk + 6*A1R*A1Rc*ab1*ab3*abb5*p2p3*pap1*papb*vijk + 4*A1L*A1Lc*ab1**2*ab2**2*abb5*p2p3*papb**2*vijk + 4*A1R*A1Rc*ab1**2*ab2**2*abb5*p2p3*papb**2*vijk - 8*A1L*A1Lc*ab1*ab2*ab3*abb5*p2p3*papb**2*vijk - 
     -      8*A1R*A1Rc*ab1*ab2*ab3*abb5*p2p3*papb**2*vijk + 4*A1L*A1Lc*ab3**2*abb5*p2p3*papb**2*vijk + 4*A1R*A1Rc*ab3**2*abb5*p2p3*papb**2*vijk - 2*A1L*A1Lc*ab1*abb5*m1**2*p2p3*pbp1*vijk - 2*A1R*A1Rc*ab1*abb5*m1**2*p2p3*pbp1*vijk - 2*A1L*A1Lc*ab1**2*ab2*abb5*m1**2*p2p3*pbp1*vijk - 
     -      2*A1R*A1Rc*ab1**2*ab2*abb5*m1**2*p2p3*pbp1*vijk + 2*A1L*A1Lc*ab1*ab3*abb5*m1**2*p2p3*pbp1*vijk + 2*A1R*A1Rc*ab1*ab3*abb5*m1**2*p2p3*pbp1*vijk + 2*A1L*A1Lc*ab1**2*abb5*p2p3*pap1*pbp1*vijk + 2*A1R*A1Rc*ab1**2*abb5*p2p3*pap1*pbp1*vijk - 
     -      2*A1L*A1Lc*ab1**2*ab2*abb5*p2p3*papb*pbp1*vijk - 2*A1R*A1Rc*ab1**2*ab2*abb5*p2p3*papb*pbp1*vijk + 2*A1L*A1Lc*ab1*ab3*abb5*p2p3*papb*pbp1*vijk + 2*A1R*A1Rc*ab1*ab3*abb5*p2p3*papb*pbp1*vijk + A1L*A1Lc*ab1**2*m1**2*p1p3**2*zim**2 + A1R*A1Rc*ab1**2*m1**2*p1p3**2*zim**2 + 
     -      A1Lc*A1R*ab1**2*m1*m3*p1p3**2*zim**2 + A1L*A1Rc*ab1**2*m1*m3*p1p3**2*zim**2 - A1L*A1Lc*ab1*m1**2*p1p3*pap3*zim**2 - A1R*A1Rc*ab1*m1**2*p1p3*pap3*zim**2 - 2*A1L*A1Lc*ab1**2*ab2*m1**2*p1p3*pap3*zim**2 - 2*A1R*A1Rc*ab1**2*ab2*m1**2*p1p3*pap3*zim**2 + 
     -      2*A1L*A1Lc*ab1*ab3*m1**2*p1p3*pap3*zim**2 + 2*A1R*A1Rc*ab1*ab3*m1**2*p1p3*pap3*zim**2 + A1L*A1Lc*ab1**2*ab2*m2**2*p1p3*pap3*zim**2 + A1R*A1Rc*ab1**2*ab2*m2**2*p1p3*pap3*zim**2 - A1L*A1Lc*ab1*ab3*m2**2*p1p3*pap3*zim**2 - A1R*A1Rc*ab1*ab3*m2**2*p1p3*pap3*zim**2 - 
     -      2*A1Lc*A1R*ab1**2*ab2*m1*m3*p1p3*pap3*zim**2 - 2*A1L*A1Rc*ab1**2*ab2*m1*m3*p1p3*pap3*zim**2 + 2*A1Lc*A1R*ab1*ab3*m1*m3*p1p3*pap3*zim**2 + 2*A1L*A1Rc*ab1*ab3*m1*m3*p1p3*pap3*zim**2 + A1L*A1Lc*ab1*m3**2*p1p3*pap3*zim**2 + A1R*A1Rc*ab1*m3**2*p1p3*pap3*zim**2 + 
     -      A1L*A1Lc*ab1*ab2*m1**2*pap3**2*zim**2 + A1R*A1Rc*ab1*ab2*m1**2*pap3**2*zim**2 + A1L*A1Lc*ab1**2*ab2**2*m1**2*pap3**2*zim**2 + A1R*A1Rc*ab1**2*ab2**2*m1**2*pap3**2*zim**2 - A1L*A1Lc*ab3*m1**2*pap3**2*zim**2 - A1R*A1Rc*ab3*m1**2*pap3**2*zim**2 - 
     -      2*A1L*A1Lc*ab1*ab2*ab3*m1**2*pap3**2*zim**2 - 2*A1R*A1Rc*ab1*ab2*ab3*m1**2*pap3**2*zim**2 + A1L*A1Lc*ab3**2*m1**2*pap3**2*zim**2 + A1R*A1Rc*ab3**2*m1**2*pap3**2*zim**2 - A1L*A1Lc*ab1*ab2*m2**2*pap3**2*zim**2 - A1R*A1Rc*ab1*ab2*m2**2*pap3**2*zim**2 - 
     -      A1L*A1Lc*ab1**2*ab2**2*m2**2*pap3**2*zim**2 - A1R*A1Rc*ab1**2*ab2**2*m2**2*pap3**2*zim**2 + A1L*A1Lc*ab3*m2**2*pap3**2*zim**2 + A1R*A1Rc*ab3*m2**2*pap3**2*zim**2 + 2*A1L*A1Lc*ab1*ab2*ab3*m2**2*pap3**2*zim**2 + 2*A1R*A1Rc*ab1*ab2*ab3*m2**2*pap3**2*zim**2 - 
     -      A1L*A1Lc*ab3**2*m2**2*pap3**2*zim**2 - A1R*A1Rc*ab3**2*m2**2*pap3**2*zim**2 - A1Lc*A1R*m1*m3*pap3**2*zim**2 - A1L*A1Rc*m1*m3*pap3**2*zim**2 + A1Lc*A1R*ab1**2*ab2**2*m1*m3*pap3**2*zim**2 + A1L*A1Rc*ab1**2*ab2**2*m1*m3*pap3**2*zim**2 - 2*A1Lc*A1R*ab1*ab2*ab3*m1*m3*pap3**2*zim**2 - 
     -      2*A1L*A1Rc*ab1*ab2*ab3*m1*m3*pap3**2*zim**2 + A1Lc*A1R*ab3**2*m1*m3*pap3**2*zim**2 + A1L*A1Rc*ab3**2*m1*m3*pap3**2*zim**2 - A1L*A1Lc*m3**2*pap3**2*zim**2 - A1R*A1Rc*m3**2*pap3**2*zim**2 - A1L*A1Lc*ab1*ab2*m3**2*pap3**2*zim**2 - A1R*A1Rc*ab1*ab2*m3**2*pap3**2*zim**2 + 
     -      A1L*A1Lc*ab3*m3**2*pap3**2*zim**2 + A1R*A1Rc*ab3*m3**2*pap3**2*zim**2 + A1L*A1Lc*ab1**2*p1p3**2*papb*zim**2 + A1R*A1Rc*ab1**2*p1p3**2*papb*zim**2 - A1L*A1Lc*ab1*p1p3*pap3*papb*zim**2 - A1R*A1Rc*ab1*p1p3*pap3*papb*zim**2 - A1L*A1Lc*ab1**2*ab2*p1p3*pap3*papb*zim**2 - 
     -      A1R*A1Rc*ab1**2*ab2*p1p3*pap3*papb*zim**2 + A1L*A1Lc*ab1*ab3*p1p3*pap3*papb*zim**2 + A1R*A1Rc*ab1*ab3*p1p3*pap3*papb*zim**2 - A1L*A1Lc*ab1**2*p1p3*pap3*pbp1*zim**2 - A1R*A1Rc*ab1**2*p1p3*pap3*pbp1*zim**2 + A1L*A1Lc*ab1*pap3**2*pbp1*zim**2 + A1R*A1Rc*ab1*pap3**2*pbp1*zim**2 + 
     -      A1L*A1Lc*ab1**2*ab2*pap3**2*pbp1*zim**2 + A1R*A1Rc*ab1**2*ab2*pap3**2*pbp1*zim**2 - A1L*A1Lc*ab1*ab3*pap3**2*pbp1*zim**2 - A1R*A1Rc*ab1*ab3*pap3**2*pbp1*zim**2 + A1L*A1Lc*ab1*m1**2*p1p3*pbp3*zim**2 + A1R*A1Rc*ab1*m1**2*p1p3*pbp3*zim**2 - 
     -      3*A1L*A1Lc*ab1**2*ab2*m1**2*p1p3*pbp3*zim**2 - 3*A1R*A1Rc*ab1**2*ab2*m1**2*p1p3*pbp3*zim**2 + 3*A1L*A1Lc*ab1*ab3*m1**2*p1p3*pbp3*zim**2 + 3*A1R*A1Rc*ab1*ab3*m1**2*p1p3*pbp3*zim**2 + 2*A1Lc*A1R*ab1*m1*m3*p1p3*pbp3*zim**2 + 2*A1L*A1Rc*ab1*m1*m3*p1p3*pbp3*zim**2 - 
     -      2*A1Lc*A1R*ab1**2*ab2*m1*m3*p1p3*pbp3*zim**2 - 2*A1L*A1Rc*ab1**2*ab2*m1*m3*p1p3*pbp3*zim**2 + 2*A1Lc*A1R*ab1*ab3*m1*m3*p1p3*pbp3*zim**2 + 2*A1L*A1Rc*ab1*ab3*m1*m3*p1p3*pbp3*zim**2 + A1L*A1Lc*ab1**2*p1p3*pap1*pbp3*zim**2 + A1R*A1Rc*ab1**2*p1p3*pap1*pbp3*zim**2 + 
     -      A1L*A1Lc*ab1*ab2*m1**2*pap3*pbp3*zim**2 + A1R*A1Rc*ab1*ab2*m1**2*pap3*pbp3*zim**2 + 3*A1L*A1Lc*ab1**2*ab2**2*m1**2*pap3*pbp3*zim**2 + 3*A1R*A1Rc*ab1**2*ab2**2*m1**2*pap3*pbp3*zim**2 - A1L*A1Lc*ab3*m1**2*pap3*pbp3*zim**2 - A1R*A1Rc*ab3*m1**2*pap3*pbp3*zim**2 - 
     -      6*A1L*A1Lc*ab1*ab2*ab3*m1**2*pap3*pbp3*zim**2 - 6*A1R*A1Rc*ab1*ab2*ab3*m1**2*pap3*pbp3*zim**2 + 3*A1L*A1Lc*ab3**2*m1**2*pap3*pbp3*zim**2 + 3*A1R*A1Rc*ab3**2*m1**2*pap3*pbp3*zim**2 + A1L*A1Lc*ab1*ab2*m2**2*pap3*pbp3*zim**2 + A1R*A1Rc*ab1*ab2*m2**2*pap3*pbp3*zim**2 - 
     -      A1L*A1Lc*ab1**2*ab2**2*m2**2*pap3*pbp3*zim**2 - A1R*A1Rc*ab1**2*ab2**2*m2**2*pap3*pbp3*zim**2 - A1L*A1Lc*ab3*m2**2*pap3*pbp3*zim**2 - A1R*A1Rc*ab3*m2**2*pap3*pbp3*zim**2 + 2*A1L*A1Lc*ab1*ab2*ab3*m2**2*pap3*pbp3*zim**2 + 2*A1R*A1Rc*ab1*ab2*ab3*m2**2*pap3*pbp3*zim**2 - 
     -      A1L*A1Lc*ab3**2*m2**2*pap3*pbp3*zim**2 - A1R*A1Rc*ab3**2*m2**2*pap3*pbp3*zim**2 - 2*A1Lc*A1R*ab1*ab2*m1*m3*pap3*pbp3*zim**2 - 2*A1L*A1Rc*ab1*ab2*m1*m3*pap3*pbp3*zim**2 + 2*A1Lc*A1R*ab1**2*ab2**2*m1*m3*pap3*pbp3*zim**2 + 2*A1L*A1Rc*ab1**2*ab2**2*m1*m3*pap3*pbp3*zim**2 + 
     -      2*A1Lc*A1R*ab3*m1*m3*pap3*pbp3*zim**2 + 2*A1L*A1Rc*ab3*m1*m3*pap3*pbp3*zim**2 - 4*A1Lc*A1R*ab1*ab2*ab3*m1*m3*pap3*pbp3*zim**2 - 4*A1L*A1Rc*ab1*ab2*ab3*m1*m3*pap3*pbp3*zim**2 + 2*A1Lc*A1R*ab3**2*m1*m3*pap3*pbp3*zim**2 + 2*A1L*A1Rc*ab3**2*m1*m3*pap3*pbp3*zim**2 + 
     -      A1L*A1Lc*m3**2*pap3*pbp3*zim**2 + A1R*A1Rc*m3**2*pap3*pbp3*zim**2 - A1L*A1Lc*ab1*ab2*m3**2*pap3*pbp3*zim**2 - A1R*A1Rc*ab1*ab2*m3**2*pap3*pbp3*zim**2 + A1L*A1Lc*ab3*m3**2*pap3*pbp3*zim**2 + A1R*A1Rc*ab3*m3**2*pap3*pbp3*zim**2 - A1L*A1Lc*ab1*pap1*pap3*pbp3*zim**2 - 
     -      A1R*A1Rc*ab1*pap1*pap3*pbp3*zim**2 - A1L*A1Lc*ab1**2*ab2*pap1*pap3*pbp3*zim**2 - A1R*A1Rc*ab1**2*ab2*pap1*pap3*pbp3*zim**2 + A1L*A1Lc*ab1*ab3*pap1*pap3*pbp3*zim**2 + A1R*A1Rc*ab1*ab3*pap1*pap3*pbp3*zim**2 + A1L*A1Lc*ab1*p1p3*papb*pbp3*zim**2 + 
     -      A1R*A1Rc*ab1*p1p3*papb*pbp3*zim**2 - 3*A1L*A1Lc*ab1**2*ab2*p1p3*papb*pbp3*zim**2 - 3*A1R*A1Rc*ab1**2*ab2*p1p3*papb*pbp3*zim**2 + 3*A1L*A1Lc*ab1*ab3*p1p3*papb*pbp3*zim**2 + 3*A1R*A1Rc*ab1*ab3*p1p3*papb*pbp3*zim**2 + 2*A1L*A1Lc*ab1*ab2*pap3*papb*pbp3*zim**2 + 
     -      2*A1R*A1Rc*ab1*ab2*pap3*papb*pbp3*zim**2 + 2*A1L*A1Lc*ab1**2*ab2**2*pap3*papb*pbp3*zim**2 + 2*A1R*A1Rc*ab1**2*ab2**2*pap3*papb*pbp3*zim**2 - 2*A1L*A1Lc*ab3*pap3*papb*pbp3*zim**2 - 2*A1R*A1Rc*ab3*pap3*papb*pbp3*zim**2 - 4*A1L*A1Lc*ab1*ab2*ab3*pap3*papb*pbp3*zim**2 - 
     -      4*A1R*A1Rc*ab1*ab2*ab3*pap3*papb*pbp3*zim**2 + 2*A1L*A1Lc*ab3**2*pap3*papb*pbp3*zim**2 + 2*A1R*A1Rc*ab3**2*pap3*papb*pbp3*zim**2 - A1L*A1Lc*ab1*pap3*pbp1*pbp3*zim**2 - A1R*A1Rc*ab1*pap3*pbp1*pbp3*zim**2 + A1L*A1Lc*ab1**2*ab2*pap3*pbp1*pbp3*zim**2 + 
     -      A1R*A1Rc*ab1**2*ab2*pap3*pbp1*pbp3*zim**2 - A1L*A1Lc*ab1*ab3*pap3*pbp1*pbp3*zim**2 - A1R*A1Rc*ab1*ab3*pap3*pbp1*pbp3*zim**2 - 2*A1L*A1Lc*ab1*ab2*m1**2*pbp3**2*zim**2 - 2*A1R*A1Rc*ab1*ab2*m1**2*pbp3**2*zim**2 + 2*A1L*A1Lc*ab1**2*ab2**2*m1**2*pbp3**2*zim**2 + 
     -      2*A1R*A1Rc*ab1**2*ab2**2*m1**2*pbp3**2*zim**2 + 2*A1L*A1Lc*ab3*m1**2*pbp3**2*zim**2 + 2*A1R*A1Rc*ab3*m1**2*pbp3**2*zim**2 - 4*A1L*A1Lc*ab1*ab2*ab3*m1**2*pbp3**2*zim**2 - 4*A1R*A1Rc*ab1*ab2*ab3*m1**2*pbp3**2*zim**2 + 2*A1L*A1Lc*ab3**2*m1**2*pbp3**2*zim**2 + 
     -      2*A1R*A1Rc*ab3**2*m1**2*pbp3**2*zim**2 + A1Lc*A1R*m1*m3*pbp3**2*zim**2 + A1L*A1Rc*m1*m3*pbp3**2*zim**2 - 2*A1Lc*A1R*ab1*ab2*m1*m3*pbp3**2*zim**2 - 2*A1L*A1Rc*ab1*ab2*m1*m3*pbp3**2*zim**2 + A1Lc*A1R*ab1**2*ab2**2*m1*m3*pbp3**2*zim**2 + 
     -      A1L*A1Rc*ab1**2*ab2**2*m1*m3*pbp3**2*zim**2 + 2*A1Lc*A1R*ab3*m1*m3*pbp3**2*zim**2 + 2*A1L*A1Rc*ab3*m1*m3*pbp3**2*zim**2 - 2*A1Lc*A1R*ab1*ab2*ab3*m1*m3*pbp3**2*zim**2 - 2*A1L*A1Rc*ab1*ab2*ab3*m1*m3*pbp3**2*zim**2 + A1Lc*A1R*ab3**2*m1*m3*pbp3**2*zim**2 + 
     -      A1L*A1Rc*ab3**2*m1*m3*pbp3**2*zim**2 + A1L*A1Lc*ab1*pap1*pbp3**2*zim**2 + A1R*A1Rc*ab1*pap1*pbp3**2*zim**2 - A1L*A1Lc*ab1**2*ab2*pap1*pbp3**2*zim**2 - A1R*A1Rc*ab1**2*ab2*pap1*pbp3**2*zim**2 + A1L*A1Lc*ab1*ab3*pap1*pbp3**2*zim**2 + A1R*A1Rc*ab1*ab3*pap1*pbp3**2*zim**2 - 
     -      2*A1L*A1Lc*ab1*ab2*papb*pbp3**2*zim**2 - 2*A1R*A1Rc*ab1*ab2*papb*pbp3**2*zim**2 + 2*A1L*A1Lc*ab1**2*ab2**2*papb*pbp3**2*zim**2 + 2*A1R*A1Rc*ab1**2*ab2**2*papb*pbp3**2*zim**2 + 2*A1L*A1Lc*ab3*papb*pbp3**2*zim**2 + 2*A1R*A1Rc*ab3*papb*pbp3**2*zim**2 - 
     -      4*A1L*A1Lc*ab1*ab2*ab3*papb*pbp3**2*zim**2 - 4*A1R*A1Rc*ab1*ab2*ab3*papb*pbp3**2*zim**2 + 2*A1L*A1Lc*ab3**2*papb*pbp3**2*zim**2 + 2*A1R*A1Rc*ab3**2*papb*pbp3**2*zim**2 - 2*A1L*A1Lc*ab1**2*m1**2*p1p2*p1p3*zim*zjm - 2*A1R*A1Rc*ab1**2*m1**2*p1p2*p1p3*zim*zjm - 
     -      2*A1Lc*A1R*ab1**2*m1*m3*p1p2*p1p3*zim*zjm - 2*A1L*A1Rc*ab1**2*m1*m3*p1p2*p1p3*zim*zjm + A1L*A1Lc*ab1*m1**2*p1p3*pap2*zim*zjm + A1R*A1Rc*ab1*m1**2*p1p3*pap2*zim*zjm + 2*A1L*A1Lc*ab1**2*ab2*m1**2*p1p3*pap2*zim*zjm + 2*A1R*A1Rc*ab1**2*ab2*m1**2*p1p3*pap2*zim*zjm - 
     -      2*A1L*A1Lc*ab1*ab3*m1**2*p1p3*pap2*zim*zjm - 2*A1R*A1Rc*ab1*ab3*m1**2*p1p3*pap2*zim*zjm - A1L*A1Lc*ab1**2*ab2*m2**2*p1p3*pap2*zim*zjm - A1R*A1Rc*ab1**2*ab2*m2**2*p1p3*pap2*zim*zjm + A1L*A1Lc*ab1*ab3*m2**2*p1p3*pap2*zim*zjm + A1R*A1Rc*ab1*ab3*m2**2*p1p3*pap2*zim*zjm + 
     -      2*A1Lc*A1R*ab1**2*ab2*m1*m3*p1p3*pap2*zim*zjm + 2*A1L*A1Rc*ab1**2*ab2*m1*m3*p1p3*pap2*zim*zjm - 2*A1Lc*A1R*ab1*ab3*m1*m3*p1p3*pap2*zim*zjm - 2*A1L*A1Rc*ab1*ab3*m1*m3*p1p3*pap2*zim*zjm - A1L*A1Lc*ab1*m3**2*p1p3*pap2*zim*zjm - A1R*A1Rc*ab1*m3**2*p1p3*pap2*zim*zjm + 
     -      A1L*A1Lc*ab1*m1**2*p1p2*pap3*zim*zjm + A1R*A1Rc*ab1*m1**2*p1p2*pap3*zim*zjm + 2*A1L*A1Lc*ab1**2*ab2*m1**2*p1p2*pap3*zim*zjm + 2*A1R*A1Rc*ab1**2*ab2*m1**2*p1p2*pap3*zim*zjm - 2*A1L*A1Lc*ab1*ab3*m1**2*p1p2*pap3*zim*zjm - 2*A1R*A1Rc*ab1*ab3*m1**2*p1p2*pap3*zim*zjm - 
     -      A1L*A1Lc*ab1**2*ab2*m2**2*p1p2*pap3*zim*zjm - A1R*A1Rc*ab1**2*ab2*m2**2*p1p2*pap3*zim*zjm + A1L*A1Lc*ab1*ab3*m2**2*p1p2*pap3*zim*zjm + A1R*A1Rc*ab1*ab3*m2**2*p1p2*pap3*zim*zjm + 2*A1Lc*A1R*ab1**2*ab2*m1*m3*p1p2*pap3*zim*zjm + 2*A1L*A1Rc*ab1**2*ab2*m1*m3*p1p2*pap3*zim*zjm - 
     -      2*A1Lc*A1R*ab1*ab3*m1*m3*p1p2*pap3*zim*zjm - 2*A1L*A1Rc*ab1*ab3*m1*m3*p1p2*pap3*zim*zjm - A1L*A1Lc*ab1*m3**2*p1p2*pap3*zim*zjm - A1R*A1Rc*ab1*m3**2*p1p2*pap3*zim*zjm - 2*A1L*A1Lc*ab1*ab2*m1**2*pap2*pap3*zim*zjm - 2*A1R*A1Rc*ab1*ab2*m1**2*pap2*pap3*zim*zjm - 
     -      2*A1L*A1Lc*ab1**2*ab2**2*m1**2*pap2*pap3*zim*zjm - 2*A1R*A1Rc*ab1**2*ab2**2*m1**2*pap2*pap3*zim*zjm + 2*A1L*A1Lc*ab3*m1**2*pap2*pap3*zim*zjm + 2*A1R*A1Rc*ab3*m1**2*pap2*pap3*zim*zjm + 4*A1L*A1Lc*ab1*ab2*ab3*m1**2*pap2*pap3*zim*zjm + 
     -      4*A1R*A1Rc*ab1*ab2*ab3*m1**2*pap2*pap3*zim*zjm - 2*A1L*A1Lc*ab3**2*m1**2*pap2*pap3*zim*zjm - 2*A1R*A1Rc*ab3**2*m1**2*pap2*pap3*zim*zjm + 2*A1L*A1Lc*ab1*ab2*m2**2*pap2*pap3*zim*zjm + 2*A1R*A1Rc*ab1*ab2*m2**2*pap2*pap3*zim*zjm + 2*A1L*A1Lc*ab1**2*ab2**2*m2**2*pap2*pap3*zim*zjm + 
     -      2*A1R*A1Rc*ab1**2*ab2**2*m2**2*pap2*pap3*zim*zjm - 2*A1L*A1Lc*ab3*m2**2*pap2*pap3*zim*zjm - 2*A1R*A1Rc*ab3*m2**2*pap2*pap3*zim*zjm - 4*A1L*A1Lc*ab1*ab2*ab3*m2**2*pap2*pap3*zim*zjm - 4*A1R*A1Rc*ab1*ab2*ab3*m2**2*pap2*pap3*zim*zjm + 2*A1L*A1Lc*ab3**2*m2**2*pap2*pap3*zim*zjm + 
     -      2*A1R*A1Rc*ab3**2*m2**2*pap2*pap3*zim*zjm + 2*A1Lc*A1R*m1*m3*pap2*pap3*zim*zjm + 2*A1L*A1Rc*m1*m3*pap2*pap3*zim*zjm - 2*A1Lc*A1R*ab1**2*ab2**2*m1*m3*pap2*pap3*zim*zjm - 2*A1L*A1Rc*ab1**2*ab2**2*m1*m3*pap2*pap3*zim*zjm + 4*A1Lc*A1R*ab1*ab2*ab3*m1*m3*pap2*pap3*zim*zjm + 
     -      4*A1L*A1Rc*ab1*ab2*ab3*m1*m3*pap2*pap3*zim*zjm - 2*A1Lc*A1R*ab3**2*m1*m3*pap2*pap3*zim*zjm - 2*A1L*A1Rc*ab3**2*m1*m3*pap2*pap3*zim*zjm + 2*A1L*A1Lc*m3**2*pap2*pap3*zim*zjm + 2*A1R*A1Rc*m3**2*pap2*pap3*zim*zjm + 2*A1L*A1Lc*ab1*ab2*m3**2*pap2*pap3*zim*zjm + 
     -      2*A1R*A1Rc*ab1*ab2*m3**2*pap2*pap3*zim*zjm - 2*A1L*A1Lc*ab3*m3**2*pap2*pap3*zim*zjm - 2*A1R*A1Rc*ab3*m3**2*pap2*pap3*zim*zjm - 2*A1L*A1Lc*ab1**2*p1p2*p1p3*papb*zim*zjm - 2*A1R*A1Rc*ab1**2*p1p2*p1p3*papb*zim*zjm + A1L*A1Lc*ab1*p1p3*pap2*papb*zim*zjm + 
     -      A1R*A1Rc*ab1*p1p3*pap2*papb*zim*zjm + A1L*A1Lc*ab1**2*ab2*p1p3*pap2*papb*zim*zjm + A1R*A1Rc*ab1**2*ab2*p1p3*pap2*papb*zim*zjm - A1L*A1Lc*ab1*ab3*p1p3*pap2*papb*zim*zjm - A1R*A1Rc*ab1*ab3*p1p3*pap2*papb*zim*zjm + A1L*A1Lc*ab1*p1p2*pap3*papb*zim*zjm + 
     -      A1R*A1Rc*ab1*p1p2*pap3*papb*zim*zjm + A1L*A1Lc*ab1**2*ab2*p1p2*pap3*papb*zim*zjm + A1R*A1Rc*ab1**2*ab2*p1p2*pap3*papb*zim*zjm - A1L*A1Lc*ab1*ab3*p1p2*pap3*papb*zim*zjm - A1R*A1Rc*ab1*ab3*p1p2*pap3*papb*zim*zjm + A1L*A1Lc*ab1**2*p1p3*pap2*pbp1*zim*zjm + 
     -      A1R*A1Rc*ab1**2*p1p3*pap2*pbp1*zim*zjm + A1L*A1Lc*ab1**2*p1p2*pap3*pbp1*zim*zjm + A1R*A1Rc*ab1**2*p1p2*pap3*pbp1*zim*zjm - 2*A1L*A1Lc*ab1*pap2*pap3*pbp1*zim*zjm - 2*A1R*A1Rc*ab1*pap2*pap3*pbp1*zim*zjm - 2*A1L*A1Lc*ab1**2*ab2*pap2*pap3*pbp1*zim*zjm - 
     -      2*A1R*A1Rc*ab1**2*ab2*pap2*pap3*pbp1*zim*zjm + 2*A1L*A1Lc*ab1*ab3*pap2*pap3*pbp1*zim*zjm + 2*A1R*A1Rc*ab1*ab3*pap2*pap3*pbp1*zim*zjm - A1L*A1Lc*ab1*m1**2*p1p3*pbp2*zim*zjm - A1R*A1Rc*ab1*m1**2*p1p3*pbp2*zim*zjm + 3*A1L*A1Lc*ab1**2*ab2*m1**2*p1p3*pbp2*zim*zjm + 
     -      3*A1R*A1Rc*ab1**2*ab2*m1**2*p1p3*pbp2*zim*zjm - 3*A1L*A1Lc*ab1*ab3*m1**2*p1p3*pbp2*zim*zjm - 3*A1R*A1Rc*ab1*ab3*m1**2*p1p3*pbp2*zim*zjm - 2*A1Lc*A1R*ab1*m1*m3*p1p3*pbp2*zim*zjm - 2*A1L*A1Rc*ab1*m1*m3*p1p3*pbp2*zim*zjm + 2*A1Lc*A1R*ab1**2*ab2*m1*m3*p1p3*pbp2*zim*zjm + 
     -      2*A1L*A1Rc*ab1**2*ab2*m1*m3*p1p3*pbp2*zim*zjm - 2*A1Lc*A1R*ab1*ab3*m1*m3*p1p3*pbp2*zim*zjm - 2*A1L*A1Rc*ab1*ab3*m1*m3*p1p3*pbp2*zim*zjm - A1L*A1Lc*ab1**2*p1p3*pap1*pbp2*zim*zjm - A1R*A1Rc*ab1**2*p1p3*pap1*pbp2*zim*zjm - A1L*A1Lc*ab1*ab2*m1**2*pap3*pbp2*zim*zjm - 
     -      A1R*A1Rc*ab1*ab2*m1**2*pap3*pbp2*zim*zjm - 3*A1L*A1Lc*ab1**2*ab2**2*m1**2*pap3*pbp2*zim*zjm - 3*A1R*A1Rc*ab1**2*ab2**2*m1**2*pap3*pbp2*zim*zjm + A1L*A1Lc*ab3*m1**2*pap3*pbp2*zim*zjm + A1R*A1Rc*ab3*m1**2*pap3*pbp2*zim*zjm + 6*A1L*A1Lc*ab1*ab2*ab3*m1**2*pap3*pbp2*zim*zjm + 
     -      6*A1R*A1Rc*ab1*ab2*ab3*m1**2*pap3*pbp2*zim*zjm - 3*A1L*A1Lc*ab3**2*m1**2*pap3*pbp2*zim*zjm - 3*A1R*A1Rc*ab3**2*m1**2*pap3*pbp2*zim*zjm - A1L*A1Lc*ab1*ab2*m2**2*pap3*pbp2*zim*zjm - A1R*A1Rc*ab1*ab2*m2**2*pap3*pbp2*zim*zjm + A1L*A1Lc*ab1**2*ab2**2*m2**2*pap3*pbp2*zim*zjm + 
     -      A1R*A1Rc*ab1**2*ab2**2*m2**2*pap3*pbp2*zim*zjm + A1L*A1Lc*ab3*m2**2*pap3*pbp2*zim*zjm + A1R*A1Rc*ab3*m2**2*pap3*pbp2*zim*zjm - 2*A1L*A1Lc*ab1*ab2*ab3*m2**2*pap3*pbp2*zim*zjm - 2*A1R*A1Rc*ab1*ab2*ab3*m2**2*pap3*pbp2*zim*zjm + A1L*A1Lc*ab3**2*m2**2*pap3*pbp2*zim*zjm + 
     -      A1R*A1Rc*ab3**2*m2**2*pap3*pbp2*zim*zjm + 2*A1Lc*A1R*ab1*ab2*m1*m3*pap3*pbp2*zim*zjm + 2*A1L*A1Rc*ab1*ab2*m1*m3*pap3*pbp2*zim*zjm - 2*A1Lc*A1R*ab1**2*ab2**2*m1*m3*pap3*pbp2*zim*zjm - 2*A1L*A1Rc*ab1**2*ab2**2*m1*m3*pap3*pbp2*zim*zjm - 2*A1Lc*A1R*ab3*m1*m3*pap3*pbp2*zim*zjm - 
     -      2*A1L*A1Rc*ab3*m1*m3*pap3*pbp2*zim*zjm + 4*A1Lc*A1R*ab1*ab2*ab3*m1*m3*pap3*pbp2*zim*zjm + 4*A1L*A1Rc*ab1*ab2*ab3*m1*m3*pap3*pbp2*zim*zjm - 2*A1Lc*A1R*ab3**2*m1*m3*pap3*pbp2*zim*zjm - 2*A1L*A1Rc*ab3**2*m1*m3*pap3*pbp2*zim*zjm - A1L*A1Lc*m3**2*pap3*pbp2*zim*zjm - 
     -      A1R*A1Rc*m3**2*pap3*pbp2*zim*zjm + A1L*A1Lc*ab1*ab2*m3**2*pap3*pbp2*zim*zjm + A1R*A1Rc*ab1*ab2*m3**2*pap3*pbp2*zim*zjm - A1L*A1Lc*ab3*m3**2*pap3*pbp2*zim*zjm - A1R*A1Rc*ab3*m3**2*pap3*pbp2*zim*zjm + A1L*A1Lc*ab1*pap1*pap3*pbp2*zim*zjm + A1R*A1Rc*ab1*pap1*pap3*pbp2*zim*zjm + 
     -      A1L*A1Lc*ab1**2*ab2*pap1*pap3*pbp2*zim*zjm + A1R*A1Rc*ab1**2*ab2*pap1*pap3*pbp2*zim*zjm - A1L*A1Lc*ab1*ab3*pap1*pap3*pbp2*zim*zjm - A1R*A1Rc*ab1*ab3*pap1*pap3*pbp2*zim*zjm - A1L*A1Lc*ab1*p1p3*papb*pbp2*zim*zjm - A1R*A1Rc*ab1*p1p3*papb*pbp2*zim*zjm + 
     -      3*A1L*A1Lc*ab1**2*ab2*p1p3*papb*pbp2*zim*zjm + 3*A1R*A1Rc*ab1**2*ab2*p1p3*papb*pbp2*zim*zjm - 3*A1L*A1Lc*ab1*ab3*p1p3*papb*pbp2*zim*zjm - 3*A1R*A1Rc*ab1*ab3*p1p3*papb*pbp2*zim*zjm - 2*A1L*A1Lc*ab1*ab2*pap3*papb*pbp2*zim*zjm - 2*A1R*A1Rc*ab1*ab2*pap3*papb*pbp2*zim*zjm - 
     -      2*A1L*A1Lc*ab1**2*ab2**2*pap3*papb*pbp2*zim*zjm - 2*A1R*A1Rc*ab1**2*ab2**2*pap3*papb*pbp2*zim*zjm + 2*A1L*A1Lc*ab3*pap3*papb*pbp2*zim*zjm + 2*A1R*A1Rc*ab3*pap3*papb*pbp2*zim*zjm + 4*A1L*A1Lc*ab1*ab2*ab3*pap3*papb*pbp2*zim*zjm + 4*A1R*A1Rc*ab1*ab2*ab3*pap3*papb*pbp2*zim*zjm - 
     -      2*A1L*A1Lc*ab3**2*pap3*papb*pbp2*zim*zjm - 2*A1R*A1Rc*ab3**2*pap3*papb*pbp2*zim*zjm + A1L*A1Lc*ab1*pap3*pbp1*pbp2*zim*zjm + A1R*A1Rc*ab1*pap3*pbp1*pbp2*zim*zjm - A1L*A1Lc*ab1**2*ab2*pap3*pbp1*pbp2*zim*zjm - A1R*A1Rc*ab1**2*ab2*pap3*pbp1*pbp2*zim*zjm + 
     -      A1L*A1Lc*ab1*ab3*pap3*pbp1*pbp2*zim*zjm + A1R*A1Rc*ab1*ab3*pap3*pbp1*pbp2*zim*zjm - A1L*A1Lc*ab1*m1**2*p1p2*pbp3*zim*zjm - A1R*A1Rc*ab1*m1**2*p1p2*pbp3*zim*zjm + 3*A1L*A1Lc*ab1**2*ab2*m1**2*p1p2*pbp3*zim*zjm + 3*A1R*A1Rc*ab1**2*ab2*m1**2*p1p2*pbp3*zim*zjm - 
     -      3*A1L*A1Lc*ab1*ab3*m1**2*p1p2*pbp3*zim*zjm - 3*A1R*A1Rc*ab1*ab3*m1**2*p1p2*pbp3*zim*zjm - 2*A1Lc*A1R*ab1*m1*m3*p1p2*pbp3*zim*zjm - 2*A1L*A1Rc*ab1*m1*m3*p1p2*pbp3*zim*zjm + 2*A1Lc*A1R*ab1**2*ab2*m1*m3*p1p2*pbp3*zim*zjm + 2*A1L*A1Rc*ab1**2*ab2*m1*m3*p1p2*pbp3*zim*zjm - 
     -      2*A1Lc*A1R*ab1*ab3*m1*m3*p1p2*pbp3*zim*zjm - 2*A1L*A1Rc*ab1*ab3*m1*m3*p1p2*pbp3*zim*zjm - A1L*A1Lc*ab1**2*p1p2*pap1*pbp3*zim*zjm - A1R*A1Rc*ab1**2*p1p2*pap1*pbp3*zim*zjm - A1L*A1Lc*ab1*ab2*m1**2*pap2*pbp3*zim*zjm - A1R*A1Rc*ab1*ab2*m1**2*pap2*pbp3*zim*zjm - 
     -      3*A1L*A1Lc*ab1**2*ab2**2*m1**2*pap2*pbp3*zim*zjm - 3*A1R*A1Rc*ab1**2*ab2**2*m1**2*pap2*pbp3*zim*zjm + A1L*A1Lc*ab3*m1**2*pap2*pbp3*zim*zjm + A1R*A1Rc*ab3*m1**2*pap2*pbp3*zim*zjm + 6*A1L*A1Lc*ab1*ab2*ab3*m1**2*pap2*pbp3*zim*zjm + 6*A1R*A1Rc*ab1*ab2*ab3*m1**2*pap2*pbp3*zim*zjm - 
     -      3*A1L*A1Lc*ab3**2*m1**2*pap2*pbp3*zim*zjm - 3*A1R*A1Rc*ab3**2*m1**2*pap2*pbp3*zim*zjm - A1L*A1Lc*ab1*ab2*m2**2*pap2*pbp3*zim*zjm - A1R*A1Rc*ab1*ab2*m2**2*pap2*pbp3*zim*zjm + A1L*A1Lc*ab1**2*ab2**2*m2**2*pap2*pbp3*zim*zjm + A1R*A1Rc*ab1**2*ab2**2*m2**2*pap2*pbp3*zim*zjm + 
     -      A1L*A1Lc*ab3*m2**2*pap2*pbp3*zim*zjm + A1R*A1Rc*ab3*m2**2*pap2*pbp3*zim*zjm - 2*A1L*A1Lc*ab1*ab2*ab3*m2**2*pap2*pbp3*zim*zjm - 2*A1R*A1Rc*ab1*ab2*ab3*m2**2*pap2*pbp3*zim*zjm + A1L*A1Lc*ab3**2*m2**2*pap2*pbp3*zim*zjm + A1R*A1Rc*ab3**2*m2**2*pap2*pbp3*zim*zjm + 
     -      2*A1Lc*A1R*ab1*ab2*m1*m3*pap2*pbp3*zim*zjm + 2*A1L*A1Rc*ab1*ab2*m1*m3*pap2*pbp3*zim*zjm - 2*A1Lc*A1R*ab1**2*ab2**2*m1*m3*pap2*pbp3*zim*zjm - 2*A1L*A1Rc*ab1**2*ab2**2*m1*m3*pap2*pbp3*zim*zjm - 2*A1Lc*A1R*ab3*m1*m3*pap2*pbp3*zim*zjm - 2*A1L*A1Rc*ab3*m1*m3*pap2*pbp3*zim*zjm + 
     -      4*A1Lc*A1R*ab1*ab2*ab3*m1*m3*pap2*pbp3*zim*zjm + 4*A1L*A1Rc*ab1*ab2*ab3*m1*m3*pap2*pbp3*zim*zjm - 2*A1Lc*A1R*ab3**2*m1*m3*pap2*pbp3*zim*zjm - 2*A1L*A1Rc*ab3**2*m1*m3*pap2*pbp3*zim*zjm - A1L*A1Lc*m3**2*pap2*pbp3*zim*zjm - A1R*A1Rc*m3**2*pap2*pbp3*zim*zjm + 
     -      A1L*A1Lc*ab1*ab2*m3**2*pap2*pbp3*zim*zjm + A1R*A1Rc*ab1*ab2*m3**2*pap2*pbp3*zim*zjm - A1L*A1Lc*ab3*m3**2*pap2*pbp3*zim*zjm - A1R*A1Rc*ab3*m3**2*pap2*pbp3*zim*zjm + A1L*A1Lc*ab1*pap1*pap2*pbp3*zim*zjm + A1R*A1Rc*ab1*pap1*pap2*pbp3*zim*zjm + 
     -      A1L*A1Lc*ab1**2*ab2*pap1*pap2*pbp3*zim*zjm + A1R*A1Rc*ab1**2*ab2*pap1*pap2*pbp3*zim*zjm - A1L*A1Lc*ab1*ab3*pap1*pap2*pbp3*zim*zjm - A1R*A1Rc*ab1*ab3*pap1*pap2*pbp3*zim*zjm - A1L*A1Lc*ab1*p1p2*papb*pbp3*zim*zjm - A1R*A1Rc*ab1*p1p2*papb*pbp3*zim*zjm + 
     -      3*A1L*A1Lc*ab1**2*ab2*p1p2*papb*pbp3*zim*zjm + 3*A1R*A1Rc*ab1**2*ab2*p1p2*papb*pbp3*zim*zjm - 3*A1L*A1Lc*ab1*ab3*p1p2*papb*pbp3*zim*zjm - 3*A1R*A1Rc*ab1*ab3*p1p2*papb*pbp3*zim*zjm - 2*A1L*A1Lc*ab1*ab2*pap2*papb*pbp3*zim*zjm - 2*A1R*A1Rc*ab1*ab2*pap2*papb*pbp3*zim*zjm - 
     -      2*A1L*A1Lc*ab1**2*ab2**2*pap2*papb*pbp3*zim*zjm - 2*A1R*A1Rc*ab1**2*ab2**2*pap2*papb*pbp3*zim*zjm + 2*A1L*A1Lc*ab3*pap2*papb*pbp3*zim*zjm + 2*A1R*A1Rc*ab3*pap2*papb*pbp3*zim*zjm + 4*A1L*A1Lc*ab1*ab2*ab3*pap2*papb*pbp3*zim*zjm + 4*A1R*A1Rc*ab1*ab2*ab3*pap2*papb*pbp3*zim*zjm - 
     -      2*A1L*A1Lc*ab3**2*pap2*papb*pbp3*zim*zjm - 2*A1R*A1Rc*ab3**2*pap2*papb*pbp3*zim*zjm + A1L*A1Lc*ab1*pap2*pbp1*pbp3*zim*zjm + A1R*A1Rc*ab1*pap2*pbp1*pbp3*zim*zjm - A1L*A1Lc*ab1**2*ab2*pap2*pbp1*pbp3*zim*zjm - A1R*A1Rc*ab1**2*ab2*pap2*pbp1*pbp3*zim*zjm + 
     -      A1L*A1Lc*ab1*ab3*pap2*pbp1*pbp3*zim*zjm + A1R*A1Rc*ab1*ab3*pap2*pbp1*pbp3*zim*zjm + 4*A1L*A1Lc*ab1*ab2*m1**2*pbp2*pbp3*zim*zjm + 4*A1R*A1Rc*ab1*ab2*m1**2*pbp2*pbp3*zim*zjm - 4*A1L*A1Lc*ab1**2*ab2**2*m1**2*pbp2*pbp3*zim*zjm - 4*A1R*A1Rc*ab1**2*ab2**2*m1**2*pbp2*pbp3*zim*zjm - 
     -      4*A1L*A1Lc*ab3*m1**2*pbp2*pbp3*zim*zjm - 4*A1R*A1Rc*ab3*m1**2*pbp2*pbp3*zim*zjm + 8*A1L*A1Lc*ab1*ab2*ab3*m1**2*pbp2*pbp3*zim*zjm + 8*A1R*A1Rc*ab1*ab2*ab3*m1**2*pbp2*pbp3*zim*zjm - 4*A1L*A1Lc*ab3**2*m1**2*pbp2*pbp3*zim*zjm - 4*A1R*A1Rc*ab3**2*m1**2*pbp2*pbp3*zim*zjm - 
     -      2*A1Lc*A1R*m1*m3*pbp2*pbp3*zim*zjm - 2*A1L*A1Rc*m1*m3*pbp2*pbp3*zim*zjm + 4*A1Lc*A1R*ab1*ab2*m1*m3*pbp2*pbp3*zim*zjm + 4*A1L*A1Rc*ab1*ab2*m1*m3*pbp2*pbp3*zim*zjm - 2*A1Lc*A1R*ab1**2*ab2**2*m1*m3*pbp2*pbp3*zim*zjm - 2*A1L*A1Rc*ab1**2*ab2**2*m1*m3*pbp2*pbp3*zim*zjm - 
     -      4*A1Lc*A1R*ab3*m1*m3*pbp2*pbp3*zim*zjm - 4*A1L*A1Rc*ab3*m1*m3*pbp2*pbp3*zim*zjm + 4*A1Lc*A1R*ab1*ab2*ab3*m1*m3*pbp2*pbp3*zim*zjm + 4*A1L*A1Rc*ab1*ab2*ab3*m1*m3*pbp2*pbp3*zim*zjm - 2*A1Lc*A1R*ab3**2*m1*m3*pbp2*pbp3*zim*zjm - 2*A1L*A1Rc*ab3**2*m1*m3*pbp2*pbp3*zim*zjm - 
     -      2*A1L*A1Lc*ab1*pap1*pbp2*pbp3*zim*zjm - 2*A1R*A1Rc*ab1*pap1*pbp2*pbp3*zim*zjm + 2*A1L*A1Lc*ab1**2*ab2*pap1*pbp2*pbp3*zim*zjm + 2*A1R*A1Rc*ab1**2*ab2*pap1*pbp2*pbp3*zim*zjm - 2*A1L*A1Lc*ab1*ab3*pap1*pbp2*pbp3*zim*zjm - 2*A1R*A1Rc*ab1*ab3*pap1*pbp2*pbp3*zim*zjm + 
     -      4*A1L*A1Lc*ab1*ab2*papb*pbp2*pbp3*zim*zjm + 4*A1R*A1Rc*ab1*ab2*papb*pbp2*pbp3*zim*zjm - 4*A1L*A1Lc*ab1**2*ab2**2*papb*pbp2*pbp3*zim*zjm - 4*A1R*A1Rc*ab1**2*ab2**2*papb*pbp2*pbp3*zim*zjm - 4*A1L*A1Lc*ab3*papb*pbp2*pbp3*zim*zjm - 4*A1R*A1Rc*ab3*papb*pbp2*pbp3*zim*zjm + 
     -      8*A1L*A1Lc*ab1*ab2*ab3*papb*pbp2*pbp3*zim*zjm + 8*A1R*A1Rc*ab1*ab2*ab3*papb*pbp2*pbp3*zim*zjm - 4*A1L*A1Lc*ab3**2*papb*pbp2*pbp3*zim*zjm - 4*A1R*A1Rc*ab3**2*papb*pbp2*pbp3*zim*zjm + A1L*A1Lc*ab1**2*m1**2*p1p2**2*zjm**2 + A1R*A1Rc*ab1**2*m1**2*p1p2**2*zjm**2 + 
     -      A1Lc*A1R*ab1**2*m1*m3*p1p2**2*zjm**2 + A1L*A1Rc*ab1**2*m1*m3*p1p2**2*zjm**2 - A1L*A1Lc*ab1*m1**2*p1p2*pap2*zjm**2 - A1R*A1Rc*ab1*m1**2*p1p2*pap2*zjm**2 - 2*A1L*A1Lc*ab1**2*ab2*m1**2*p1p2*pap2*zjm**2 - 2*A1R*A1Rc*ab1**2*ab2*m1**2*p1p2*pap2*zjm**2 + 
     -      2*A1L*A1Lc*ab1*ab3*m1**2*p1p2*pap2*zjm**2 + 2*A1R*A1Rc*ab1*ab3*m1**2*p1p2*pap2*zjm**2 + A1L*A1Lc*ab1**2*ab2*m2**2*p1p2*pap2*zjm**2 + A1R*A1Rc*ab1**2*ab2*m2**2*p1p2*pap2*zjm**2 - A1L*A1Lc*ab1*ab3*m2**2*p1p2*pap2*zjm**2 - A1R*A1Rc*ab1*ab3*m2**2*p1p2*pap2*zjm**2 - 
     -      2*A1Lc*A1R*ab1**2*ab2*m1*m3*p1p2*pap2*zjm**2 - 2*A1L*A1Rc*ab1**2*ab2*m1*m3*p1p2*pap2*zjm**2 + 2*A1Lc*A1R*ab1*ab3*m1*m3*p1p2*pap2*zjm**2 + 2*A1L*A1Rc*ab1*ab3*m1*m3*p1p2*pap2*zjm**2 + A1L*A1Lc*ab1*m3**2*p1p2*pap2*zjm**2 + A1R*A1Rc*ab1*m3**2*p1p2*pap2*zjm**2 + 
     -      A1L*A1Lc*ab1*ab2*m1**2*pap2**2*zjm**2 + A1R*A1Rc*ab1*ab2*m1**2*pap2**2*zjm**2 + A1L*A1Lc*ab1**2*ab2**2*m1**2*pap2**2*zjm**2 + A1R*A1Rc*ab1**2*ab2**2*m1**2*pap2**2*zjm**2 - A1L*A1Lc*ab3*m1**2*pap2**2*zjm**2 - A1R*A1Rc*ab3*m1**2*pap2**2*zjm**2 - 
     -      2*A1L*A1Lc*ab1*ab2*ab3*m1**2*pap2**2*zjm**2 - 2*A1R*A1Rc*ab1*ab2*ab3*m1**2*pap2**2*zjm**2 + A1L*A1Lc*ab3**2*m1**2*pap2**2*zjm**2 + A1R*A1Rc*ab3**2*m1**2*pap2**2*zjm**2 - A1L*A1Lc*ab1*ab2*m2**2*pap2**2*zjm**2 - A1R*A1Rc*ab1*ab2*m2**2*pap2**2*zjm**2 - 
     -      A1L*A1Lc*ab1**2*ab2**2*m2**2*pap2**2*zjm**2 - A1R*A1Rc*ab1**2*ab2**2*m2**2*pap2**2*zjm**2 + A1L*A1Lc*ab3*m2**2*pap2**2*zjm**2 + A1R*A1Rc*ab3*m2**2*pap2**2*zjm**2 + 2*A1L*A1Lc*ab1*ab2*ab3*m2**2*pap2**2*zjm**2 + 2*A1R*A1Rc*ab1*ab2*ab3*m2**2*pap2**2*zjm**2 - 
     -      A1L*A1Lc*ab3**2*m2**2*pap2**2*zjm**2 - A1R*A1Rc*ab3**2*m2**2*pap2**2*zjm**2 - A1Lc*A1R*m1*m3*pap2**2*zjm**2 - A1L*A1Rc*m1*m3*pap2**2*zjm**2 + A1Lc*A1R*ab1**2*ab2**2*m1*m3*pap2**2*zjm**2 + A1L*A1Rc*ab1**2*ab2**2*m1*m3*pap2**2*zjm**2 - 2*A1Lc*A1R*ab1*ab2*ab3*m1*m3*pap2**2*zjm**2 - 
     -      2*A1L*A1Rc*ab1*ab2*ab3*m1*m3*pap2**2*zjm**2 + A1Lc*A1R*ab3**2*m1*m3*pap2**2*zjm**2 + A1L*A1Rc*ab3**2*m1*m3*pap2**2*zjm**2 - A1L*A1Lc*m3**2*pap2**2*zjm**2 - A1R*A1Rc*m3**2*pap2**2*zjm**2 - A1L*A1Lc*ab1*ab2*m3**2*pap2**2*zjm**2 - A1R*A1Rc*ab1*ab2*m3**2*pap2**2*zjm**2 + 
     -      A1L*A1Lc*ab3*m3**2*pap2**2*zjm**2 + A1R*A1Rc*ab3*m3**2*pap2**2*zjm**2 + A1L*A1Lc*ab1**2*p1p2**2*papb*zjm**2 + A1R*A1Rc*ab1**2*p1p2**2*papb*zjm**2 - A1L*A1Lc*ab1*p1p2*pap2*papb*zjm**2 - A1R*A1Rc*ab1*p1p2*pap2*papb*zjm**2 - A1L*A1Lc*ab1**2*ab2*p1p2*pap2*papb*zjm**2 - 
     -      A1R*A1Rc*ab1**2*ab2*p1p2*pap2*papb*zjm**2 + A1L*A1Lc*ab1*ab3*p1p2*pap2*papb*zjm**2 + A1R*A1Rc*ab1*ab3*p1p2*pap2*papb*zjm**2 - A1L*A1Lc*ab1**2*p1p2*pap2*pbp1*zjm**2 - A1R*A1Rc*ab1**2*p1p2*pap2*pbp1*zjm**2 + A1L*A1Lc*ab1*pap2**2*pbp1*zjm**2 + A1R*A1Rc*ab1*pap2**2*pbp1*zjm**2 + 
     -      A1L*A1Lc*ab1**2*ab2*pap2**2*pbp1*zjm**2 + A1R*A1Rc*ab1**2*ab2*pap2**2*pbp1*zjm**2 - A1L*A1Lc*ab1*ab3*pap2**2*pbp1*zjm**2 - A1R*A1Rc*ab1*ab3*pap2**2*pbp1*zjm**2 + A1L*A1Lc*ab1*m1**2*p1p2*pbp2*zjm**2 + A1R*A1Rc*ab1*m1**2*p1p2*pbp2*zjm**2 - 
     -      3*A1L*A1Lc*ab1**2*ab2*m1**2*p1p2*pbp2*zjm**2 - 3*A1R*A1Rc*ab1**2*ab2*m1**2*p1p2*pbp2*zjm**2 + 3*A1L*A1Lc*ab1*ab3*m1**2*p1p2*pbp2*zjm**2 + 3*A1R*A1Rc*ab1*ab3*m1**2*p1p2*pbp2*zjm**2 + 2*A1Lc*A1R*ab1*m1*m3*p1p2*pbp2*zjm**2 + 2*A1L*A1Rc*ab1*m1*m3*p1p2*pbp2*zjm**2 - 
     -      2*A1Lc*A1R*ab1**2*ab2*m1*m3*p1p2*pbp2*zjm**2 - 2*A1L*A1Rc*ab1**2*ab2*m1*m3*p1p2*pbp2*zjm**2 + 2*A1Lc*A1R*ab1*ab3*m1*m3*p1p2*pbp2*zjm**2 + 2*A1L*A1Rc*ab1*ab3*m1*m3*p1p2*pbp2*zjm**2 + A1L*A1Lc*ab1**2*p1p2*pap1*pbp2*zjm**2 + A1R*A1Rc*ab1**2*p1p2*pap1*pbp2*zjm**2 + 
     -      A1L*A1Lc*ab1*ab2*m1**2*pap2*pbp2*zjm**2 + A1R*A1Rc*ab1*ab2*m1**2*pap2*pbp2*zjm**2 + 3*A1L*A1Lc*ab1**2*ab2**2*m1**2*pap2*pbp2*zjm**2 + 3*A1R*A1Rc*ab1**2*ab2**2*m1**2*pap2*pbp2*zjm**2 - A1L*A1Lc*ab3*m1**2*pap2*pbp2*zjm**2 - A1R*A1Rc*ab3*m1**2*pap2*pbp2*zjm**2 - 
     -      6*A1L*A1Lc*ab1*ab2*ab3*m1**2*pap2*pbp2*zjm**2 - 6*A1R*A1Rc*ab1*ab2*ab3*m1**2*pap2*pbp2*zjm**2 + 3*A1L*A1Lc*ab3**2*m1**2*pap2*pbp2*zjm**2 + 3*A1R*A1Rc*ab3**2*m1**2*pap2*pbp2*zjm**2 + A1L*A1Lc*ab1*ab2*m2**2*pap2*pbp2*zjm**2 + A1R*A1Rc*ab1*ab2*m2**2*pap2*pbp2*zjm**2 - 
     -      A1L*A1Lc*ab1**2*ab2**2*m2**2*pap2*pbp2*zjm**2 - A1R*A1Rc*ab1**2*ab2**2*m2**2*pap2*pbp2*zjm**2 - A1L*A1Lc*ab3*m2**2*pap2*pbp2*zjm**2 - A1R*A1Rc*ab3*m2**2*pap2*pbp2*zjm**2 + 2*A1L*A1Lc*ab1*ab2*ab3*m2**2*pap2*pbp2*zjm**2 + 2*A1R*A1Rc*ab1*ab2*ab3*m2**2*pap2*pbp2*zjm**2 - 
     -      A1L*A1Lc*ab3**2*m2**2*pap2*pbp2*zjm**2 - A1R*A1Rc*ab3**2*m2**2*pap2*pbp2*zjm**2 - 2*A1Lc*A1R*ab1*ab2*m1*m3*pap2*pbp2*zjm**2 - 2*A1L*A1Rc*ab1*ab2*m1*m3*pap2*pbp2*zjm**2 + 2*A1Lc*A1R*ab1**2*ab2**2*m1*m3*pap2*pbp2*zjm**2 + 2*A1L*A1Rc*ab1**2*ab2**2*m1*m3*pap2*pbp2*zjm**2 + 
     -      2*A1Lc*A1R*ab3*m1*m3*pap2*pbp2*zjm**2 + 2*A1L*A1Rc*ab3*m1*m3*pap2*pbp2*zjm**2 - 4*A1Lc*A1R*ab1*ab2*ab3*m1*m3*pap2*pbp2*zjm**2 - 4*A1L*A1Rc*ab1*ab2*ab3*m1*m3*pap2*pbp2*zjm**2 + 2*A1Lc*A1R*ab3**2*m1*m3*pap2*pbp2*zjm**2 + 2*A1L*A1Rc*ab3**2*m1*m3*pap2*pbp2*zjm**2 + 
     -      A1L*A1Lc*m3**2*pap2*pbp2*zjm**2 + A1R*A1Rc*m3**2*pap2*pbp2*zjm**2 - A1L*A1Lc*ab1*ab2*m3**2*pap2*pbp2*zjm**2 - A1R*A1Rc*ab1*ab2*m3**2*pap2*pbp2*zjm**2 + A1L*A1Lc*ab3*m3**2*pap2*pbp2*zjm**2 + A1R*A1Rc*ab3*m3**2*pap2*pbp2*zjm**2 - A1L*A1Lc*ab1*pap1*pap2*pbp2*zjm**2 - 
     -      A1R*A1Rc*ab1*pap1*pap2*pbp2*zjm**2 - A1L*A1Lc*ab1**2*ab2*pap1*pap2*pbp2*zjm**2 - A1R*A1Rc*ab1**2*ab2*pap1*pap2*pbp2*zjm**2 + A1L*A1Lc*ab1*ab3*pap1*pap2*pbp2*zjm**2 + A1R*A1Rc*ab1*ab3*pap1*pap2*pbp2*zjm**2 + A1L*A1Lc*ab1*p1p2*papb*pbp2*zjm**2 + 
     -      A1R*A1Rc*ab1*p1p2*papb*pbp2*zjm**2 - 3*A1L*A1Lc*ab1**2*ab2*p1p2*papb*pbp2*zjm**2 - 3*A1R*A1Rc*ab1**2*ab2*p1p2*papb*pbp2*zjm**2 + 3*A1L*A1Lc*ab1*ab3*p1p2*papb*pbp2*zjm**2 + 3*A1R*A1Rc*ab1*ab3*p1p2*papb*pbp2*zjm**2 + 2*A1L*A1Lc*ab1*ab2*pap2*papb*pbp2*zjm**2 + 
     -      2*A1R*A1Rc*ab1*ab2*pap2*papb*pbp2*zjm**2 + 2*A1L*A1Lc*ab1**2*ab2**2*pap2*papb*pbp2*zjm**2 + 2*A1R*A1Rc*ab1**2*ab2**2*pap2*papb*pbp2*zjm**2 - 2*A1L*A1Lc*ab3*pap2*papb*pbp2*zjm**2 - 2*A1R*A1Rc*ab3*pap2*papb*pbp2*zjm**2 - 4*A1L*A1Lc*ab1*ab2*ab3*pap2*papb*pbp2*zjm**2 - 
     -      4*A1R*A1Rc*ab1*ab2*ab3*pap2*papb*pbp2*zjm**2 + 2*A1L*A1Lc*ab3**2*pap2*papb*pbp2*zjm**2 + 2*A1R*A1Rc*ab3**2*pap2*papb*pbp2*zjm**2 - A1L*A1Lc*ab1*pap2*pbp1*pbp2*zjm**2 - A1R*A1Rc*ab1*pap2*pbp1*pbp2*zjm**2 + A1L*A1Lc*ab1**2*ab2*pap2*pbp1*pbp2*zjm**2 + 
     -      A1R*A1Rc*ab1**2*ab2*pap2*pbp1*pbp2*zjm**2 - A1L*A1Lc*ab1*ab3*pap2*pbp1*pbp2*zjm**2 - A1R*A1Rc*ab1*ab3*pap2*pbp1*pbp2*zjm**2 - 2*A1L*A1Lc*ab1*ab2*m1**2*pbp2**2*zjm**2 - 2*A1R*A1Rc*ab1*ab2*m1**2*pbp2**2*zjm**2 + 2*A1L*A1Lc*ab1**2*ab2**2*m1**2*pbp2**2*zjm**2 + 
     -      2*A1R*A1Rc*ab1**2*ab2**2*m1**2*pbp2**2*zjm**2 + 2*A1L*A1Lc*ab3*m1**2*pbp2**2*zjm**2 + 2*A1R*A1Rc*ab3*m1**2*pbp2**2*zjm**2 - 4*A1L*A1Lc*ab1*ab2*ab3*m1**2*pbp2**2*zjm**2 - 4*A1R*A1Rc*ab1*ab2*ab3*m1**2*pbp2**2*zjm**2 + 2*A1L*A1Lc*ab3**2*m1**2*pbp2**2*zjm**2 + 
     -      2*A1R*A1Rc*ab3**2*m1**2*pbp2**2*zjm**2 + A1Lc*A1R*m1*m3*pbp2**2*zjm**2 + A1L*A1Rc*m1*m3*pbp2**2*zjm**2 - 2*A1Lc*A1R*ab1*ab2*m1*m3*pbp2**2*zjm**2 - 2*A1L*A1Rc*ab1*ab2*m1*m3*pbp2**2*zjm**2 + A1Lc*A1R*ab1**2*ab2**2*m1*m3*pbp2**2*zjm**2 + 
     -      A1L*A1Rc*ab1**2*ab2**2*m1*m3*pbp2**2*zjm**2 + 2*A1Lc*A1R*ab3*m1*m3*pbp2**2*zjm**2 + 2*A1L*A1Rc*ab3*m1*m3*pbp2**2*zjm**2 - 2*A1Lc*A1R*ab1*ab2*ab3*m1*m3*pbp2**2*zjm**2 - 2*A1L*A1Rc*ab1*ab2*ab3*m1*m3*pbp2**2*zjm**2 + A1Lc*A1R*ab3**2*m1*m3*pbp2**2*zjm**2 + 
     -      A1L*A1Rc*ab3**2*m1*m3*pbp2**2*zjm**2 + A1L*A1Lc*ab1*pap1*pbp2**2*zjm**2 + A1R*A1Rc*ab1*pap1*pbp2**2*zjm**2 - A1L*A1Lc*ab1**2*ab2*pap1*pbp2**2*zjm**2 - A1R*A1Rc*ab1**2*ab2*pap1*pbp2**2*zjm**2 + A1L*A1Lc*ab1*ab3*pap1*pbp2**2*zjm**2 + A1R*A1Rc*ab1*ab3*pap1*pbp2**2*zjm**2 - 
     -      2*A1L*A1Lc*ab1*ab2*papb*pbp2**2*zjm**2 - 2*A1R*A1Rc*ab1*ab2*papb*pbp2**2*zjm**2 + 2*A1L*A1Lc*ab1**2*ab2**2*papb*pbp2**2*zjm**2 + 2*A1R*A1Rc*ab1**2*ab2**2*papb*pbp2**2*zjm**2 + 2*A1L*A1Lc*ab3*papb*pbp2**2*zjm**2 + 2*A1R*A1Rc*ab3*papb*pbp2**2*zjm**2 - 
     -      4*A1L*A1Lc*ab1*ab2*ab3*papb*pbp2**2*zjm**2 - 4*A1R*A1Rc*ab1*ab2*ab3*papb*pbp2**2*zjm**2 + 2*A1L*A1Lc*ab3**2*papb*pbp2**2*zjm**2 + 2*A1R*A1Rc*ab3**2*papb*pbp2**2*zjm**2))/
     -  (p2p3*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) + ab1*(-pap1 + ab2*(m1**2 + m2**2 + 2*papb) - pbp1))*((1 + ab1*ab2 - ab3)*m2**2 + papb + ab1*ab2*papb - ab3*papb - ab1*pbp1)*vijk))


      end

      double precision function NeuQ2qx_M2_FEFS_qqbar()
      implicit none
    
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
      
      double precision zi,zj,A1L,A1R,A1Lc,A1Rc,zP,zM,dipcolor,zim,zjm
      double precision ab1,ab2,ab3,abb5,muQ2,yijk,vijk,viji,kappa,Q2

      A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
      A1R = AR(chimeigen,sfmeigen,sftype,sfgen)
      A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
      A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))

      Q2 = sman
      muQ2 = m3**2/Q2

      yijk = p2p3/(p1p2 + p1p3 + p2p3) 
      vijk = dsqrt((2*muQ2+(1-muQ2)*(1-yijk))**2-4*muQ2)/((1-muQ2)*(1-yijk)) 
      viji = 1d0 
      zP = (1+viji*vijk)/2d0 
      zM = (1-viji*vijk)/2d0 


      zi = p1p3/(p1p2 + p1p3)
      zj = 1-zi

      zim = zi - zM
      zjm = zj - zM

      ab1 = kappa(Q2, 0d0, m3**2)/kappa(Q2, 2*p2p3, m3**2)
      ab2 = (pap1+pbp1)/sman
      ab3 = (Q2 + m3**2)/(2*Q2)

      abb5 = 1-2*kappaDip*zP*zM

      dipcolor =  -1/2d0


      NeuQ2qx_M2_FEFS_qqbar =  dipcolor/(2*p2p3*vijk)*((gs**4*TF*(2*A1Lc*A1R*abb5*m1*m3**3*p2p3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2 + 2*A1L*A1Rc*abb5*m1*m3**3*p2p3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2 - 
     -      2*A1Lc*A1R*abb5*m1*m3**3*Nc**2*p2p3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2 - 2*A1L*A1Rc*abb5*m1*m3**3*Nc**2*p2p3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2 + 
     -      2*A1L*A1Lc*abb5*m3**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2 + 
     -      2*A1R*A1Rc*abb5*m3**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2 - 
     -      2*A1L*A1Lc*abb5*m3**2*Nc**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2 - 
     -      2*A1R*A1Rc*abb5*m3**2*Nc**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2 + 
     -      2*A1L*A1Lc*abb5*m3**2*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2 + 
     -      2*A1R*A1Rc*abb5*m3**2*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2 - 
     -      2*A1L*A1Lc*abb5*m3**2*Nc**2*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2 - 
     -      2*A1R*A1Rc*abb5*m3**2*Nc**2*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2 + 
     -      2*A1L*A1Lc*abb5*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) + 
     -      2*A1R*A1Rc*abb5*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) - 
     -      2*A1L*A1Lc*abb5*Nc**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) - 
     -      2*A1R*A1Rc*abb5*Nc**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) + 
     -      2*A1L*A1Lc*abb5*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) + 
     -      2*A1R*A1Rc*abb5*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) - 
     -      2*A1L*A1Lc*abb5*Nc**2*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) - 
     -      2*A1R*A1Rc*abb5*Nc**2*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) + 
     -      4*A1Lc*A1R*abb5*m1*m3*p2p3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(ab3*(m2**2 + papb) + ab1*(-(ab2*(m2**2 + papb)) + pbp1))*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) + 
     -      4*A1L*A1Rc*abb5*m1*m3*p2p3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(ab3*(m2**2 + papb) + ab1*(-(ab2*(m2**2 + papb)) + pbp1))*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) - 
     -      4*A1Lc*A1R*abb5*m1*m3*Nc**2*p2p3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(ab3*(m2**2 + papb) + ab1*(-(ab2*(m2**2 + papb)) + pbp1))*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) - 
     -      4*A1L*A1Rc*abb5*m1*m3*Nc**2*p2p3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(ab3*(m2**2 + papb) + ab1*(-(ab2*(m2**2 + papb)) + pbp1))*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) + 
     -      2*A1L*A1Lc*abb5*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(ab3*(m2**2 + papb) + ab1*(-(ab2*(m2**2 + papb)) + pbp1))*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) + 2*A1R*A1Rc*abb5*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (ab3*(m2**2 + papb) + ab1*(-(ab2*(m2**2 + papb)) + pbp1))*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) - 
     -      2*A1L*A1Lc*abb5*Nc**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(ab3*(m2**2 + papb) + ab1*(-(ab2*(m2**2 + papb)) + pbp1))*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) - 2*A1R*A1Rc*abb5*Nc**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (ab3*(m2**2 + papb) + ab1*(-(ab2*(m2**2 + papb)) + pbp1))*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) + 
     -      4*A1L*A1Lc*abb5*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(ab3*(m2**2 + papb) + ab1*(-(ab2*(m2**2 + papb)) + pbp1))*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) + 4*A1R*A1Rc*abb5*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (ab3*(m2**2 + papb) + ab1*(-(ab2*(m2**2 + papb)) + pbp1))*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) - 
     -      4*A1L*A1Lc*abb5*Nc**2*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(ab3*(m2**2 + papb) + ab1*(-(ab2*(m2**2 + papb)) + pbp1))*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) - 4*A1R*A1Rc*abb5*Nc**2*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (ab3*(m2**2 + papb) + ab1*(-(ab2*(m2**2 + papb)) + pbp1))*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1)) + 2*A1Lc*A1R*abb5*m1*m2**2*m3*p2p3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2 + 
     -      2*A1L*A1Rc*abb5*m1*m2**2*m3*p2p3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2 - 2*A1Lc*A1R*abb5*m1*m2**2*m3*Nc**2*p2p3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2 - 
     -      2*A1L*A1Rc*abb5*m1*m2**2*m3*Nc**2*p2p3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2 + 
     -      2*A1L*A1Lc*abb5*m2**2*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2 + 
     -      2*A1R*A1Rc*abb5*m2**2*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2 - 
     -      2*A1L*A1Lc*abb5*m2**2*Nc**2*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2 - 
     -      2*A1R*A1Rc*abb5*m2**2*Nc**2*p2p3*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2 + 
     -      2*A1L*A1Lc*abb5*p2p3*papb*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2 + 
     -      2*A1R*A1Rc*abb5*p2p3*papb*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2 - 
     -      2*A1L*A1Lc*abb5*Nc**2*p2p3*papb*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2 - 
     -      2*A1R*A1Rc*abb5*Nc**2*p2p3*papb*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2 + 
     -      4*A1L*A1Lc*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3*zim**2 + 
     -      4*A1R*A1Rc*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3*zim**2 - 
     -      4*A1L*A1Lc*Nc**2*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3*zim**2 - 
     -      4*A1R*A1Rc*Nc**2*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3*zim**2 + 
     -      4*A1Lc*A1R*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3**2*zim**2 + 4*A1L*A1Rc*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3**2*zim**2 - 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3**2*zim**2 - 4*A1L*A1Rc*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3**2*zim**2 + 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3**2*zim**2 + 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3**2*zim**2 - 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3**2*zim**2 - 
     -      4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3**2*zim**2 - 
     -      2*A1L*A1Lc*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 - 
     -      2*A1R*A1Rc*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 + 
     -      2*A1L*A1Lc*Nc**2*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 + 
     -      2*A1R*A1Rc*Nc**2*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 + 
     -      4*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 + 
     -      4*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 - 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 - 
     -      4*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 + 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 + 4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 - 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 - 4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 - 
     -      4*A1Lc*A1R*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 - 
     -      4*A1L*A1Rc*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 + 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 + 
     -      4*A1L*A1Rc*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 - 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 - 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 + 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 + 
     -      4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim**2 - 
     -      2*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))**2*zim**2 - 
     -      2*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))**2*zim**2 + 
     -      2*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))**2*zim**2 + 
     -      2*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))**2*zim**2 - 
     -      2*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))**2*zim**2 - 2*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))**2*zim**2 + 
     -      2*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))**2*zim**2 + 2*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))**2*zim**2 + 
     -      A1Lc*A1R*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))**2*zim**2 + 
     -      A1L*A1Rc*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))**2*zim**2 - 
     -      A1Lc*A1R*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))**2*zim**2 - 
     -      A1L*A1Rc*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))**2*zim**2 + 
     -      A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))**2*zim**2 + 
     -      A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))**2*zim**2 - 
     -      A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))**2*zim**2 - 
     -      A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))**2*zim**2 + 
     -      4*A1L*A1Lc*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      4*A1R*A1Rc*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      4*A1L*A1Lc*Nc**2*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      4*A1R*A1Rc*Nc**2*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      8*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      8*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      8*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      8*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      4*A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 4*A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      4*A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 4*A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      8*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 8*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      8*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 8*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp3*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      4*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      4*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      4*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      4*A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 4*A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      4*A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 4*A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*
     -       zim**2 - 4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 4*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      4*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 4*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      4*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 2*A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      2*A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      2*A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      2*A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 - 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim**2 + 
     -      4*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 + 
     -      4*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 - 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 - 
     -      4*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 + 
     -      4*A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 + 
     -      4*A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 - 
     -      4*A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 - 
     -      4*A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 + 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 + 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 - 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 - 
     -      4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))**2*zim**2 - 
     -      4*A1L*A1Lc*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*zim*zjm - 
     -      4*A1R*A1Rc*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*zim*zjm + 
     -      4*A1L*A1Lc*Nc**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*zim*zjm + 
     -      4*A1R*A1Rc*Nc**2*p2p3*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*zim*zjm - 
     -      4*A1L*A1Lc*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*zim*zjm - 
     -      4*A1R*A1Rc*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*zim*zjm + 
     -      4*A1L*A1Lc*Nc**2*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*zim*zjm + 
     -      4*A1R*A1Rc*Nc**2*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*zim*zjm + 
     -      2*A1L*A1Lc*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zim*zjm + 
     -      2*A1R*A1Rc*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zim*zjm - 
     -      2*A1L*A1Lc*Nc**2*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zim*zjm - 
     -      2*A1R*A1Rc*Nc**2*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zim*zjm - 
     -      4*A1L*A1Lc*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zim*zjm - 
     -      4*A1R*A1Rc*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zim*zjm + 
     -      4*A1L*A1Lc*Nc**2*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zim*zjm + 
     -      4*A1R*A1Rc*Nc**2*pap3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zim*zjm - 
     -      4*A1L*A1Lc*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3*zim*zjm - 
     -      4*A1R*A1Rc*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3*zim*zjm + 
     -      4*A1L*A1Lc*Nc**2*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3*zim*zjm + 
     -      4*A1R*A1Rc*Nc**2*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp3*zim*zjm - 
     -      8*A1Lc*A1R*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*pbp3*zim*zjm - 8*A1L*A1Rc*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*pbp3*zim*zjm + 
     -      8*A1Lc*A1R*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*pbp3*zim*zjm + 8*A1L*A1Rc*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*pbp3*zim*zjm - 
     -      8*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*pbp3*zim*zjm - 
     -      8*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*pbp3*zim*zjm + 
     -      8*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*pbp3*zim*zjm + 
     -      8*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*pbp3*zim*zjm - 
     -      4*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*pbp3*zim*zjm - 
     -      4*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*pbp3*zim*zjm + 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*pbp3*zim*zjm + 
     -      4*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*pbp3*zim*zjm - 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*pbp3*zim*zjm - 4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*pbp3*zim*zjm + 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*pbp3*zim*zjm + 4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*pbp3*zim*zjm + 
     -      4*A1Lc*A1R*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*pbp3*zim*zjm + 
     -      4*A1L*A1Rc*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*pbp3*zim*zjm - 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*pbp3*zim*zjm - 
     -      4*A1L*A1Rc*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*pbp3*zim*zjm + 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*pbp3*zim*zjm + 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*pbp3*zim*zjm - 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*pbp3*zim*zjm - 
     -      4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*pbp3*zim*zjm - 
     -      8*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*pbp3*zim*zjm - 
     -      8*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*pbp3*zim*zjm + 
     -      8*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*pbp3*zim*zjm + 
     -      8*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*pbp3*zim*zjm - 
     -      4*A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*pbp3*zim*zjm - 4*A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*pbp3*zim*zjm + 
     -      4*A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*pbp3*zim*zjm + 4*A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*pbp3*zim*zjm - 
     -      8*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*pbp3*zim*zjm - 8*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*pbp3*zim*zjm + 
     -      8*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*pbp3*zim*zjm + 8*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*pbp3*zim*zjm + 
     -      2*A1L*A1Lc*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      2*A1R*A1Rc*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      2*A1L*A1Lc*Nc**2*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      2*A1R*A1Rc*Nc**2*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      4*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      4*A1Lc*A1R*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      4*A1L*A1Rc*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1L*A1Rc*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      4*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 4*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 4*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      2*A1Lc*A1R*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      2*A1L*A1Rc*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      2*A1Lc*A1R*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      2*A1L*A1Rc*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      2*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 2*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      2*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 2*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      4*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 4*A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      4*A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 4*A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*
     -       zjm + 4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*
     -       zim*zjm + 4*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 4*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*
     -       (-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 4*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      2*A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      2*A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      2*A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      2*A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm + 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(-(ab1*p1p3) + pap3 + pbp3 + ab1*ab2*(pap3 + pbp3) - ab3*(pap3 + pbp3))*zim*zjm - 
     -      4*A1L*A1Lc*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      4*A1R*A1Rc*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      4*A1L*A1Lc*Nc**2*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      4*A1R*A1Rc*Nc**2*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      8*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      8*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      8*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      8*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      4*A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 4*A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      4*A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 4*A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      8*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 8*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      8*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 8*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      4*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      4*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      4*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      4*A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 4*A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      4*A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 4*A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*
     -       zjm + 4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*
     -       zim*zjm + 4*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 4*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*
     -       (ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 4*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      2*A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      2*A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      2*A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      2*A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      8*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      8*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      8*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      8*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      8*A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      8*A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      8*A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      8*A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      8*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm - 
     -      8*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      8*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      8*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*(ab3*(pap3 + pbp3) + ab1*(p1p3 - ab2*(pap3 + pbp3)))*zim*zjm + 
     -      4*A1L*A1Lc*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*zjm**2 + 
     -      4*A1R*A1Rc*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*zjm**2 - 
     -      4*A1L*A1Lc*Nc**2*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*zjm**2 - 
     -      4*A1R*A1Rc*Nc**2*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*zjm**2 + 
     -      4*A1Lc*A1R*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2**2*zjm**2 + 4*A1L*A1Rc*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2**2*zjm**2 - 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2**2*zjm**2 - 4*A1L*A1Rc*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2**2*zjm**2 + 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2**2*zjm**2 + 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2**2*zjm**2 - 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2**2*zjm**2 - 
     -      4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2**2*zjm**2 - 
     -      2*A1L*A1Lc*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 - 
     -      2*A1R*A1Rc*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 + 
     -      2*A1L*A1Lc*Nc**2*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 + 
     -      2*A1R*A1Rc*Nc**2*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 + 
     -      4*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 + 
     -      4*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 - 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 - 
     -      4*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 + 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 + 4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 - 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 - 4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 - 
     -      4*A1Lc*A1R*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 - 
     -      4*A1L*A1Rc*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 + 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 + 
     -      4*A1L*A1Rc*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 - 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 - 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 + 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 + 
     -      4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*pbp2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*zjm**2 - 
     -      2*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))**2*zjm**2 - 
     -      2*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))**2*zjm**2 + 
     -      2*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))**2*zjm**2 + 
     -      2*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))**2*zjm**2 - 
     -      2*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))**2*zjm**2 - 2*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))**2*zjm**2 + 
     -      2*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))**2*zjm**2 + 2*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))**2*zjm**2 + 
     -      A1Lc*A1R*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))**2*zjm**2 + 
     -      A1L*A1Rc*m1*m3*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))**2*zjm**2 - 
     -      A1Lc*A1R*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))**2*zjm**2 - 
     -      A1L*A1Rc*m1*m3*Nc**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))**2*zjm**2 + 
     -      A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))**2*zjm**2 + 
     -      A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))**2*zjm**2 - 
     -      A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))**2*zjm**2 - 
     -      A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))**2*zjm**2 + 
     -      4*A1L*A1Lc*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      4*A1R*A1Rc*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      4*A1L*A1Lc*Nc**2*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      4*A1R*A1Rc*Nc**2*pap2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      8*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      8*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      8*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      8*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      4*A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 4*A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      4*A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 4*A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      8*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 8*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      8*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 8*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*pbp2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      4*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      4*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      4*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      4*A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 4*A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      4*A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 4*A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*
     -       zjm**2 - 4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 4*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      4*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 4*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      4*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*
     -       (ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 2*A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*
     -       (m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*(-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      2*A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      2*A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      2*A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 - 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)*(m3**2 - ab3*(m1**2 + m2**2 + 2*papb) - ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))*
     -       (-(ab1*p1p2) + pap2 + pbp2 + ab1*ab2*(pap2 + pbp2) - ab3*(pap2 + pbp2))*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))*zjm**2 + 
     -      4*A1Lc*A1R*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 + 
     -      4*A1L*A1Rc*m1*m3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 - 
     -      4*A1Lc*A1R*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 - 
     -      4*A1L*A1Rc*m1*m3*Nc**2*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 + 
     -      4*A1L*A1Lc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 + 
     -      4*A1R*A1Rc*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 - 
     -      4*A1L*A1Lc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 - 
     -      4*A1R*A1Rc*Nc**2*(m1**2 - ab1*pap1 + papb + ab1*ab2*(m1**2 + papb) - ab3*(m1**2 + papb))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 + 
     -      4*A1L*A1Lc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 + 
     -      4*A1R*A1Rc*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 - 
     -      4*A1L*A1Lc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2 - 
     -      4*A1R*A1Rc*Nc**2*(ab3*(m1**2 + papb) + ab1*(pap1 - ab2*(m1**2 + papb)))*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(ab3*(pap2 + pbp2) + ab1*(p1p2 - ab2*(pap2 + pbp2)))**2*zjm**2))/
     -  (p2p3*(m2**2 + papb + ab1*ab2*(m2**2 + papb) - ab3*(m2**2 + papb) - ab1*pbp1)**2*(-m3**2 + ab3*(m1**2 + m2**2 + 2*papb) + ab1*(pap1 - ab2*(m1**2 + m2**2 + 2*papb) + pbp1))**2))

      end



      double precision function NeuQ2qx_M2_FEIS_qqbar()
      implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

      double precision zi,zj,A1L,A1R,A1Lc,A1Rc
      double precision zP,zM,Rdip,ab1,zim,zjm,dipcolor,xdip,Q2


      A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
      A1R = AR(chimeigen,sfmeigen,sftype,sfgen)
      A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
      A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))

      xdip = (-p2p3 + pbp2 + pbp3)/(pbp2 + pbp3) 

      Q2 = m1**2 + m3**2 - 2*pap1 
      Rdip = dsqrt((Q2 - m2**2 + 2*m2**2*xdip)**2 - 4*m2**2*Q2*xdip**2)/(-Q2+m2**2) 

      zP = (1+Rdip)/2d0
      zM = (1-Rdip)/2d0

      zi = pbp3/(pbp2 + pbp3)
      zj = 1-zi

      zim = zi-zM
      zjm = zj-zM

      ab1 = (1d0/Rdip*(Q2-m2**2 + 2*m2**2*xdip)/(2*Q2) - (Q2 + m2**2)/(2*Q2) + 1)

      dipcolor =  -1/2d0

      NeuQ2qx_M2_FEIS_qqbar = dipcolor/(2*p2p3*xdip)*((2*A1Lc*A1R*gs**4*m1*m2**2*m3*TF)/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2 + (2*A1L*A1Rc*gs**4*m1*m2**2*m3*TF)/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2 - 
     -  (2*A1Lc*A1R*gs**4*m1*m2**2*m3*Nc**2*TF)/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2 - (2*A1L*A1Rc*gs**4*m1*m2**2*m3*Nc**2*TF)/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2 + 
     -  (2*A1L*A1Lc*gs**4*m2**2*pap1*TF)/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2 + (2*A1R*A1Rc*gs**4*m2**2*pap1*TF)/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2 - 
     -  (2*A1L*A1Lc*gs**4*m2**2*Nc**2*pap1*TF)/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2 - (2*A1R*A1Rc*gs**4*m2**2*Nc**2*pap1*TF)/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2 + 
     -  (2*A1L*A1Lc*gs**4*TF*(-m1**2 + ab1*(m1**2 - pap1) + pap1 + (papb*xdip)/Rdip))/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip) + 
     -  (2*A1R*A1Rc*gs**4*TF*(-m1**2 + ab1*(m1**2 - pap1) + pap1 + (papb*xdip)/Rdip))/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip) - 
     -  (2*A1L*A1Lc*gs**4*Nc**2*TF*(-m1**2 + ab1*(m1**2 - pap1) + pap1 + (papb*xdip)/Rdip))/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip) - 
     -  (2*A1R*A1Rc*gs**4*Nc**2*TF*(-m1**2 + ab1*(m1**2 - pap1) + pap1 + (papb*xdip)/Rdip))/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip) + (2*A1Lc*A1R*gs**4*m1*m3**3*TF)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 + 
     -  (2*A1L*A1Rc*gs**4*m1*m3**3*TF)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 - (2*A1Lc*A1R*gs**4*m1*m3**3*Nc**2*TF)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 - (2*A1L*A1Rc*gs**4*m1*m3**3*Nc**2*TF)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 + 
     -  (2*A1L*A1Lc*gs**4*m3**2*pap1*TF)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 + (2*A1R*A1Rc*gs**4*m3**2*pap1*TF)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 - (2*A1L*A1Lc*gs**4*m3**2*Nc**2*pap1*TF)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 - 
     -  (2*A1R*A1Rc*gs**4*m3**2*Nc**2*pap1*TF)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 + (2*A1L*A1Lc*gs**4*m3**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip))/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 + 
     -  (2*A1R*A1Rc*gs**4*m3**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip))/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 - (2*A1L*A1Lc*gs**4*m3**2*Nc**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip))/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 - 
     -  (2*A1R*A1Rc*gs**4*m3**2*Nc**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip))/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 - (2*A1L*A1Lc*gs**4*pap1*TF)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) - (2*A1R*A1Rc*gs**4*pap1*TF)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) + 
     -  (2*A1L*A1Lc*gs**4*Nc**2*pap1*TF)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) + (2*A1R*A1Rc*gs**4*Nc**2*pap1*TF)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) - (2*A1L*A1Lc*gs**4*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip))/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) - 
     -  (2*A1R*A1Rc*gs**4*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip))/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) + (2*A1L*A1Lc*gs**4*Nc**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip))/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) + 
     -  (2*A1R*A1Rc*gs**4*Nc**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip))/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*TF*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*gs**4*m1*m3*TF*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*pap1*TF*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*pap1*TF*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*Nc**2*pap1*TF*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*Nc**2*pap1*TF*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*gs**4*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*gs**4*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*gs**4*Nc**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*gs**4*Nc**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*p1p3**2*TF*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1L*A1Rc*gs**4*m1*m3*p1p3**2*TF*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*p1p3**2*TF*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*p1p3**2*TF*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*p1p3**2*pap1*TF*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1R*A1Rc*gs**4*p1p3**2*pap1*TF*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p3**2*pap1*TF*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1R*A1Rc*gs**4*Nc**2*p1p3**2*pap1*TF*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*p1p3**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1R*A1Rc*gs**4*p1p3**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p3**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1R*A1Rc*gs**4*Nc**2*p1p3**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*p1p3*pap3*TF*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - (4*A1R*A1Rc*gs**4*p1p3*pap3*TF*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + (4*A1L*A1Lc*gs**4*Nc**2*p1p3*pap3*TF*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*Nc**2*p1p3*pap3*TF*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - (2*A1L*A1Lc*gs**4*pap3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*gs**4*pap3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*gs**4*Nc**2*pap3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*gs**4*Nc**2*pap3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*p1p3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1L*A1Rc*gs**4*m1*m3*p1p3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*p1p3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*p1p3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*p1p3*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1R*A1Rc*gs**4*p1p3*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p3*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1R*A1Rc*gs**4*Nc**2*p1p3*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*gs**4*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*gs**4*Nc**2*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*p1p3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*gs**4*m1*m3*p1p3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*p1p3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*p1p3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*p1p3*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*p1p3*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p3*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*Nc**2*p1p3*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*gs**4*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*gs**4*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*gs**4*Nc**2*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*gs**4*Nc**2*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (A1Lc*A1R*gs**4*m1*m3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (A1L*A1Rc*gs**4*m1*m3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (A1L*A1Lc*gs**4*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (A1R*A1Rc*gs**4*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (A1L*A1Lc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (A1R*A1Rc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (2*A1Lc*A1R*gs**4*m1*m3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Rc*gs**4*m1*m3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*gs**4*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*gs**4*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*pap3*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*pap3*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*pap3*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*Nc**2*pap3*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (8*A1Lc*A1R*gs**4*m1*m3*p1p3*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Rc*gs**4*m1*m3*p1p3*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1Lc*A1R*gs**4*m1*m3*Nc**2*p1p3*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Rc*gs**4*m1*m3*Nc**2*p1p3*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Lc*gs**4*p1p3*pap1*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1R*A1Rc*gs**4*p1p3*pap1*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Lc*gs**4*Nc**2*p1p3*pap1*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1R*A1Rc*gs**4*Nc**2*p1p3*pap1*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*Nc**2*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Rc*gs**4*m1*m3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*gs**4*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*gs**4*m1*m3*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Rc*gs**4*m1*m3*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*pap1*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*gs**4*pap1*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*pap1*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*gs**4*Nc**2*pap1*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (8*A1Lc*A1R*gs**4*m1*m3*p1p2*p1p3*TF*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (8*A1L*A1Rc*gs**4*m1*m3*p1p2*p1p3*TF*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (8*A1Lc*A1R*gs**4*m1*m3*Nc**2*p1p2*p1p3*TF*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (8*A1L*A1Rc*gs**4*m1*m3*Nc**2*p1p2*p1p3*TF*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (8*A1L*A1Lc*gs**4*p1p2*p1p3*pap1*TF*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (8*A1R*A1Rc*gs**4*p1p2*p1p3*pap1*TF*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (8*A1L*A1Lc*gs**4*Nc**2*p1p2*p1p3*pap1*TF*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (8*A1R*A1Rc*gs**4*Nc**2*p1p2*p1p3*pap1*TF*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (8*A1L*A1Lc*gs**4*p1p2*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (8*A1R*A1Rc*gs**4*p1p2*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (8*A1L*A1Lc*gs**4*Nc**2*p1p2*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (8*A1R*A1Rc*gs**4*Nc**2*p1p2*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*p1p3*pap2*TF*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + (4*A1R*A1Rc*gs**4*p1p3*pap2*TF*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - (4*A1L*A1Lc*gs**4*Nc**2*p1p3*pap2*TF*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*Nc**2*p1p3*pap2*TF*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + (4*A1L*A1Lc*gs**4*p1p2*pap3*TF*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + (4*A1R*A1Rc*gs**4*p1p2*pap3*TF*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p2*pap3*TF*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - (4*A1R*A1Rc*gs**4*Nc**2*p1p2*pap3*TF*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim*zjm)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) + (4*A1R*A1Rc*gs**4*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim*zjm)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim*zjm)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) - (4*A1R*A1Rc*gs**4*Nc**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim*zjm)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) + 
     -  (2*A1L*A1Lc*gs**4*pap3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*gs**4*pap3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*gs**4*Nc**2*pap3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*gs**4*Nc**2*pap3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*p1p3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1L*A1Rc*gs**4*m1*m3*p1p3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*p1p3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*p1p3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*p1p3*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1R*A1Rc*gs**4*p1p3*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p3*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1R*A1Rc*gs**4*Nc**2*p1p3*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*gs**4*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*gs**4*Nc**2*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*p1p3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*gs**4*m1*m3*p1p3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*p1p3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*p1p3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*p1p3*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*p1p3*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p3*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*Nc**2*p1p3*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*gs**4*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*gs**4*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*gs**4*Nc**2*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*gs**4*Nc**2*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*pap3*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*pap3*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*Nc**2*pap3*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*Nc**2*pap3*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (8*A1Lc*A1R*gs**4*m1*m3*p1p3*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Rc*gs**4*m1*m3*p1p3*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1Lc*A1R*gs**4*m1*m3*Nc**2*p1p3*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Rc*gs**4*m1*m3*Nc**2*p1p3*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Lc*gs**4*p1p3*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1R*A1Rc*gs**4*p1p3*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Lc*gs**4*Nc**2*p1p3*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1R*A1Rc*gs**4*Nc**2*p1p3*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*Nc**2*p1p3*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*gs**4*pap2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*gs**4*pap2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*gs**4*Nc**2*pap2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*gs**4*Nc**2*pap2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*p1p2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1L*A1Rc*gs**4*m1*m3*p1p2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*p1p2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*p1p2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*p1p2*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1R*A1Rc*gs**4*p1p2*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p2*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1R*A1Rc*gs**4*Nc**2*p1p2*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*gs**4*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*gs**4*Nc**2*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*p1p2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*gs**4*m1*m3*p1p2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*p1p2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*p1p2*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*p1p2*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*p1p2*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p2*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*Nc**2*p1p2*pap1*TF*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*gs**4*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*gs**4*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*gs**4*Nc**2*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*gs**4*Nc**2*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1Lc*A1R*gs**4*m1*m3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (2*A1L*A1Rc*gs**4*m1*m3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (2*A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (2*A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (2*A1L*A1Lc*gs**4*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (2*A1R*A1Rc*gs**4*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (2*A1L*A1Lc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (2*A1R*A1Rc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*gs**4*m1*m3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Rc*gs**4*m1*m3*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*gs**4*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*gs**4*Nc**2*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*gs**4*m1*m3*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*Nc**2*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*pap2*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*pap2*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*Nc**2*pap2*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*Nc**2*pap2*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (8*A1Lc*A1R*gs**4*m1*m3*p1p2*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Rc*gs**4*m1*m3*p1p2*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1Lc*A1R*gs**4*m1*m3*Nc**2*p1p2*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Rc*gs**4*m1*m3*Nc**2*p1p2*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Lc*gs**4*p1p2*pap1*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1R*A1Rc*gs**4*p1p2*pap1*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Lc*gs**4*Nc**2*p1p2*pap1*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1R*A1Rc*gs**4*Nc**2*p1p2*pap1*TF*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*Nc**2*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Rc*gs**4*m1*m3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*gs**4*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*gs**4*m1*m3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1Lc*A1R*gs**4*m1*m3*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (8*A1L*A1Rc*gs**4*m1*m3*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (8*A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (8*A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (8*A1L*A1Lc*gs**4*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (8*A1R*A1Rc*gs**4*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (8*A1L*A1Lc*gs**4*Nc**2*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (8*A1R*A1Rc*gs**4*Nc**2*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*p1p2**2*TF*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1L*A1Rc*gs**4*m1*m3*p1p2**2*TF*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*p1p2**2*TF*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*p1p2**2*TF*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*p1p2**2*pap1*TF*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1R*A1Rc*gs**4*p1p2**2*pap1*TF*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p2**2*pap1*TF*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1R*A1Rc*gs**4*Nc**2*p1p2**2*pap1*TF*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*p1p2**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1R*A1Rc*gs**4*p1p2**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p2**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1R*A1Rc*gs**4*Nc**2*p1p2**2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*p1p2*pap2*TF*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - (4*A1R*A1Rc*gs**4*p1p2*pap2*TF*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + (4*A1L*A1Lc*gs**4*Nc**2*p1p2*pap2*TF*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*Nc**2*p1p2*pap2*TF*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - (2*A1L*A1Lc*gs**4*pap2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*gs**4*pap2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*gs**4*Nc**2*pap2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*gs**4*Nc**2*pap2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*p1p2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1L*A1Rc*gs**4*m1*m3*p1p2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*p1p2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*p1p2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*p1p2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1R*A1Rc*gs**4*p1p2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1R*A1Rc*gs**4*Nc**2*p1p2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*gs**4*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*gs**4*Nc**2*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*p1p2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*gs**4*m1*m3*p1p2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*p1p2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*p1p2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*p1p2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*p1p2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*Nc**2*p1p2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*gs**4*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*gs**4*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*gs**4*Nc**2*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*gs**4*Nc**2*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (A1Lc*A1R*gs**4*m1*m3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (A1L*A1Rc*gs**4*m1*m3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (A1L*A1Lc*gs**4*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (A1R*A1Rc*gs**4*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (A1L*A1Lc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (A1R*A1Rc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (2*A1Lc*A1R*gs**4*m1*m3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Rc*gs**4*m1*m3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*gs**4*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*gs**4*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*pap2*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*pap2*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*pap2*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*Nc**2*pap2*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (8*A1Lc*A1R*gs**4*m1*m3*p1p2*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Rc*gs**4*m1*m3*p1p2*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1Lc*A1R*gs**4*m1*m3*Nc**2*p1p2*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Rc*gs**4*m1*m3*Nc**2*p1p2*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Lc*gs**4*p1p2*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1R*A1Rc*gs**4*p1p2*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Lc*gs**4*Nc**2*p1p2*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1R*A1Rc*gs**4*Nc**2*p1p2*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*Nc**2*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*Nc**2*p1p2*TF*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Rc*gs**4*m1*m3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*gs**4*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*gs**4*m1*m3*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*gs**4*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*gs**4*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*gs**4*Nc**2*pap1*TF*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*gs**4*m1*m3*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Rc*gs**4*m1*m3*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*gs**4*m1*m3*Nc**2*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Rc*gs**4*m1*m3*Nc**2*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*gs**4*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*gs**4*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*gs**4*Nc**2*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*gs**4*Nc**2*pap1*TF*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2))

      end  

      double precision function NeuQ2qx_M2_FEIS_gg()
      implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

      double precision zi,zj,A1L,A1R,A1Lc,A1Rc,Q2,xdip
      double precision zP,zM,Rdip,ab1,ab2,zim,zjm,dipcolor


      A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
      A1R = AR(chimeigen,sfmeigen,sftype,sfgen)
      A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
      A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))

      xdip = (-p2p3 + pbp2 + pbp3)/(pbp2 + pbp3) 
      Q2 = m1**2 + m3**2 - 2*pap1 
      Rdip = dsqrt((Q2 - m2**2 + 2*m2**2*xdip)**2 - 4*m2**2*Q2*xdip**2)/(-Q2+m2**2) 

      zP = (1+Rdip)/2d0
      zM = (1-Rdip)/2d0

      zi = pbp3/(pbp2 + pbp3)
      zj = 1 - zi

      zim = zi-zM
      zjm = zj-zM

      ab1 = (1d0/Rdip*(Q2-m2**2 + 2*m2**2*xdip)/(2*Q2) - (Q2 + m2**2)/(2*Q2) + 1)
      ab2 = -2 + 1/(1 - xdip + zi) + 1/(2 - xdip - zi)

      dipcolor =  -1/2d0
   
! Ms^2+Mt^2+2MsMt
      NeuQ2qx_M2_FEIS_gg = dipcolor/(2*p2p3*xdip)*((4*A1Lc*A1R*ab2*CV*gs**4*m1*m2**2*m3)/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2 + (4*A1L*A1Rc*ab2*CV*gs**4*m1*m2**2*m3)/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2 - 
     -  (4*A1Lc*A1R*ab2*CV*gs**4*m1*m2**2*m3*Nc**2)/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2 - (4*A1L*A1Rc*ab2*CV*gs**4*m1*m2**2*m3*Nc**2)/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2 + 
     -  (4*A1L*A1Lc*ab2*CV*gs**4*m2**2*pap1)/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2 + (4*A1R*A1Rc*ab2*CV*gs**4*m2**2*pap1)/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2 - 
     -  (4*A1L*A1Lc*ab2*CV*gs**4*m2**2*Nc**2*pap1)/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2 - (4*A1R*A1Rc*ab2*CV*gs**4*m2**2*Nc**2*pap1)/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2 + 
     -  (4*A1L*A1Lc*ab2*CV*gs**4*(-m1**2 + ab1*(m1**2 - pap1) + pap1 + (papb*xdip)/Rdip))/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip) + 
     -  (4*A1R*A1Rc*ab2*CV*gs**4*(-m1**2 + ab1*(m1**2 - pap1) + pap1 + (papb*xdip)/Rdip))/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip) - 
     -  (4*A1L*A1Lc*ab2*CV*gs**4*Nc**2*(-m1**2 + ab1*(m1**2 - pap1) + pap1 + (papb*xdip)/Rdip))/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip) - 
     -  (4*A1R*A1Rc*ab2*CV*gs**4*Nc**2*(-m1**2 + ab1*(m1**2 - pap1) + pap1 + (papb*xdip)/Rdip))/(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip) + (4*A1Lc*A1R*ab2*CV*gs**4*m1*m3**3)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 + 
     -  (4*A1L*A1Rc*ab2*CV*gs**4*m1*m3**3)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 - (4*A1Lc*A1R*ab2*CV*gs**4*m1*m3**3*Nc**2)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 - (4*A1L*A1Rc*ab2*CV*gs**4*m1*m3**3*Nc**2)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 + 
     -  (4*A1L*A1Lc*ab2*CV*gs**4*m3**2*pap1)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 + (4*A1R*A1Rc*ab2*CV*gs**4*m3**2*pap1)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 - (4*A1L*A1Lc*ab2*CV*gs**4*m3**2*Nc**2*pap1)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 - 
     -  (4*A1R*A1Rc*ab2*CV*gs**4*m3**2*Nc**2*pap1)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 + (4*A1L*A1Lc*ab2*CV*gs**4*m3**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip))/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 + 
     -  (4*A1R*A1Rc*ab2*CV*gs**4*m3**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip))/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 - (4*A1L*A1Lc*ab2*CV*gs**4*m3**2*Nc**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip))/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 - 
     -  (4*A1R*A1Rc*ab2*CV*gs**4*m3**2*Nc**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip))/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2 - (4*A1L*A1Lc*ab2*CV*gs**4*pap1)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) - (4*A1R*A1Rc*ab2*CV*gs**4*pap1)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) + 
     -  (4*A1L*A1Lc*ab2*CV*gs**4*Nc**2*pap1)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) + (4*A1R*A1Rc*ab2*CV*gs**4*Nc**2*pap1)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) - (4*A1L*A1Lc*ab2*CV*gs**4*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip))/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) - 
     -  (4*A1R*A1Rc*ab2*CV*gs**4*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip))/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) + (4*A1L*A1Lc*ab2*CV*gs**4*Nc**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip))/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) + 
     -  (4*A1R*A1Rc*ab2*CV*gs**4*Nc**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip))/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) - 
     -  (8*A1Lc*A1R*ab2*CV*gs**4*m1*m3*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Rc*ab2*CV*gs**4*m1*m3*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1Lc*A1R*ab2*CV*gs**4*m1*m3*Nc**2*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Rc*ab2*CV*gs**4*m1*m3*Nc**2*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Lc*ab2*CV*gs**4*pap1*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1R*A1Rc*ab2*CV*gs**4*pap1*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Lc*ab2*CV*gs**4*Nc**2*pap1*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1R*A1Rc*ab2*CV*gs**4*Nc**2*pap1*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*ab2*CV*gs**4*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*ab2*CV*gs**4*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*ab2*CV*gs**4*Nc**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*ab2*CV*gs**4*Nc**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(m3**2 - pap1 + ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip))/((-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*p1p3**2*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1L*A1Rc*CV*gs**4*m1*m3*p1p3**2*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*p1p3**2*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*p1p3**2*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*p1p3**2*pap1*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1R*A1Rc*CV*gs**4*p1p3**2*pap1*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p3**2*pap1*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p3**2*pap1*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*p1p3**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1R*A1Rc*CV*gs**4*p1p3**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p3**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p3**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*p1p3*pap3*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + (4*A1R*A1Rc*CV*gs**4*p1p3*pap3*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p3*pap3*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p3*pap3*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + (2*A1L*A1Lc*CV*gs**4*pap3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*CV*gs**4*pap3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*CV*gs**4*Nc**2*pap3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*CV*gs**4*Nc**2*pap3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*p1p3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1L*A1Rc*CV*gs**4*m1*m3*p1p3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*p1p3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*p1p3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*p1p3*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1R*A1Rc*CV*gs**4*p1p3*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p3*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p3*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*CV*gs**4*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*p1p3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*p1p3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*p1p3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*p1p3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*p1p3*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*p1p3*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p3*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p3*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*CV*gs**4*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*CV*gs**4*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*CV*gs**4*Nc**2*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*CV*gs**4*Nc**2*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (A1Lc*A1R*CV*gs**4*m1*m3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (A1L*A1Rc*CV*gs**4*m1*m3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (A1L*A1Lc*CV*gs**4*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (A1R*A1Rc*CV*gs**4*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (A1L*A1Lc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (A1R*A1Rc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (2*A1Lc*A1R*CV*gs**4*m1*m3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Rc*CV*gs**4*m1*m3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*CV*gs**4*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*CV*gs**4*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*pap3*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*pap3*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*pap3*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*pap3*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (8*A1Lc*A1R*CV*gs**4*m1*m3*p1p3*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Rc*CV*gs**4*m1*m3*p1p3*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*p1p3*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*p1p3*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Lc*CV*gs**4*p1p3*pap1*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1R*A1Rc*CV*gs**4*p1p3*pap1*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Lc*CV*gs**4*Nc**2*p1p3*pap1*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1R*A1Rc*CV*gs**4*Nc**2*p1p3*pap1*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*CV*gs**4*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*pap1*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*CV*gs**4*pap1*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*pap1*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*pap1*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)**2*zim**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (8*A1Lc*A1R*CV*gs**4*m1*m3*p1p2*p1p3*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (8*A1L*A1Rc*CV*gs**4*m1*m3*p1p2*p1p3*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (8*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*p1p2*p1p3*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (8*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*p1p2*p1p3*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (8*A1L*A1Lc*CV*gs**4*p1p2*p1p3*pap1*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (8*A1R*A1Rc*CV*gs**4*p1p2*p1p3*pap1*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (8*A1L*A1Lc*CV*gs**4*Nc**2*p1p2*p1p3*pap1*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (8*A1R*A1Rc*CV*gs**4*Nc**2*p1p2*p1p3*pap1*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (8*A1L*A1Lc*CV*gs**4*p1p2*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (8*A1R*A1Rc*CV*gs**4*p1p2*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (8*A1L*A1Lc*CV*gs**4*Nc**2*p1p2*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (8*A1R*A1Rc*CV*gs**4*Nc**2*p1p2*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*p1p3*pap2*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - (4*A1R*A1Rc*CV*gs**4*p1p3*pap2*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p3*pap2*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p3*pap2*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - (4*A1L*A1Lc*CV*gs**4*p1p2*pap3*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - (4*A1R*A1Rc*CV*gs**4*p1p2*pap3*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p2*pap3*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p2*pap3*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim*zjm)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) - (4*A1R*A1Rc*CV*gs**4*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim*zjm)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim*zjm)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) + (4*A1R*A1Rc*CV*gs**4*Nc**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zim*zjm)/(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip) - 
     -  (2*A1L*A1Lc*CV*gs**4*pap3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*CV*gs**4*pap3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*CV*gs**4*Nc**2*pap3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*CV*gs**4*Nc**2*pap3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*p1p3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1L*A1Rc*CV*gs**4*m1*m3*p1p3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*p1p3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*p1p3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*p1p3*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1R*A1Rc*CV*gs**4*p1p3*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p3*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p3*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*CV*gs**4*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*p1p3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*p1p3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*p1p3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*p1p3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*p1p3*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*p1p3*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p3*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p3*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*CV*gs**4*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*CV*gs**4*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*CV*gs**4*Nc**2*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*CV*gs**4*Nc**2*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*pap3*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*pap3*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*pap3*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*pap3*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (8*A1Lc*A1R*CV*gs**4*m1*m3*p1p3*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Rc*CV*gs**4*m1*m3*p1p3*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*p1p3*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*p1p3*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Lc*CV*gs**4*p1p3*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1R*A1Rc*CV*gs**4*p1p3*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Lc*CV*gs**4*Nc**2*p1p3*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1R*A1Rc*CV*gs**4*Nc**2*p1p3*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p3*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*CV*gs**4*pap2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*CV*gs**4*pap2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*CV*gs**4*Nc**2*pap2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*CV*gs**4*Nc**2*pap2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*p1p2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1L*A1Rc*CV*gs**4*m1*m3*p1p2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*p1p2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*p1p2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*p1p2*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1R*A1Rc*CV*gs**4*p1p2*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p2*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p2*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*CV*gs**4*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*p1p2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*p1p2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*p1p2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*p1p2*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*p1p2*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*p1p2*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p2*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p2*pap1*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*CV*gs**4*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*CV*gs**4*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*CV*gs**4*Nc**2*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*CV*gs**4*Nc**2*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1Lc*A1R*CV*gs**4*m1*m3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (2*A1L*A1Rc*CV*gs**4*m1*m3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (2*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (2*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (2*A1L*A1Lc*CV*gs**4*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (2*A1R*A1Rc*CV*gs**4*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (2*A1L*A1Lc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (2*A1R*A1Rc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*CV*gs**4*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*pap2*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*pap2*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*pap2*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*pap2*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (8*A1Lc*A1R*CV*gs**4*m1*m3*p1p2*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Rc*CV*gs**4*m1*m3*p1p2*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*p1p2*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*p1p2*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Lc*CV*gs**4*p1p2*pap1*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1R*A1Rc*CV*gs**4*p1p2*pap1*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Lc*CV*gs**4*Nc**2*p1p2*pap1*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1R*A1Rc*CV*gs**4*Nc**2*p1p2*pap1*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*CV*gs**4*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1Lc*A1R*CV*gs**4*m1*m3*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (8*A1L*A1Rc*CV*gs**4*m1*m3*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (8*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (8*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (8*A1L*A1Lc*CV*gs**4*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (8*A1R*A1Rc*CV*gs**4*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (8*A1L*A1Lc*CV*gs**4*Nc**2*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (8*A1R*A1Rc*CV*gs**4*Nc**2*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p3 - pap3 + ab1*(-p1p3 + pap3) + (pbp3*xdip)/Rdip)*zim*zjm)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*p1p2**2*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1L*A1Rc*CV*gs**4*m1*m3*p1p2**2*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*p1p2**2*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*p1p2**2*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*p1p2**2*pap1*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1R*A1Rc*CV*gs**4*p1p2**2*pap1*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p2**2*pap1*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p2**2*pap1*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*p1p2**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1R*A1Rc*CV*gs**4*p1p2**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p2**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p2**2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*p1p2*pap2*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + (4*A1R*A1Rc*CV*gs**4*p1p2*pap2*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p2*pap2*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p2*pap2*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + (2*A1L*A1Lc*CV*gs**4*pap2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*CV*gs**4*pap2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*CV*gs**4*Nc**2*pap2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*CV*gs**4*Nc**2*pap2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*p1p2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1L*A1Rc*CV*gs**4*m1*m3*p1p2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*p1p2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*p1p2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*p1p2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - (4*A1R*A1Rc*CV*gs**4*p1p2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*CV*gs**4*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*p1p2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*p1p2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*p1p2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*p1p2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*p1p2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*p1p2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*CV*gs**4*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*CV*gs**4*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*CV*gs**4*Nc**2*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*CV*gs**4*Nc**2*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (A1Lc*A1R*CV*gs**4*m1*m3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (A1L*A1Rc*CV*gs**4*m1*m3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (A1L*A1Lc*CV*gs**4*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (A1R*A1Rc*CV*gs**4*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (A1L*A1Lc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (A1R*A1Rc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (2*A1Lc*A1R*CV*gs**4*m1*m3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Rc*CV*gs**4*m1*m3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1L*A1Lc*CV*gs**4*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (2*A1R*A1Rc*CV*gs**4*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1L*A1Lc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (2*A1R*A1Rc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*pap2*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*pap2*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*pap2*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*pap2*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)) + 
     -  (8*A1Lc*A1R*CV*gs**4*m1*m3*p1p2*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Rc*CV*gs**4*m1*m3*p1p2*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*p1p2*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*p1p2*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1L*A1Lc*CV*gs**4*p1p2*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (8*A1R*A1Rc*CV*gs**4*p1p2*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1L*A1Lc*CV*gs**4*Nc**2*p1p2*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (8*A1R*A1Rc*CV*gs**4*Nc**2*p1p2*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*p1p2*(ab1*(m1**2 - pap1) + (papb*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*CV*gs**4*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1L*A1Lc*CV*gs**4*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) + 
     -  (4*A1R*A1Rc*CV*gs**4*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*pap1*(ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)*(ab1*(-m3**2 + pap1) + (pbp1*xdip)/Rdip)) - 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1Lc*A1R*CV*gs**4*m1*m3*Nc**2*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Rc*CV*gs**4*m1*m3*Nc**2*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1L*A1Lc*CV*gs**4*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) - 
     -  (4*A1R*A1Rc*CV*gs**4*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1L*A1Lc*CV*gs**4*Nc**2*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2) + 
     -  (4*A1R*A1Rc*CV*gs**4*Nc**2*pap1*(p1p2 - pap2 + ab1*(-p1p2 + pap2) + (pbp2*xdip)/Rdip)**2*zjm**2)/(p2p3*(-m1**2 + m2**2 - m3**2 + ab1*(m1**2 + m3**2 - 2*pap1) + 2*pap1 + ((papb - pbp1)*xdip)/Rdip)**2))

      end

