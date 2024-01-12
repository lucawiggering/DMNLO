      double precision function NeuQ2qx_M2DipoleEW()
      implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

      double precision Viesingular,Viefinite,Vfesingular,Vfefinite
      double precision NeuQ2qx_M2Tree,v, aie, bfe, qbar
      double complex cspen

!======================================================================================
! Returns the x=1 dipole contribution for an EW particle in the final state (gamma,Z,W,h,H,A,H+)
!================ Soft contribution and finite terms ================
 
 ! IE-FS ( Gluon emmission from an initial state squark with a final state quark as spectator ): described by VQsingular and VQfinite.
 ! FE-IS ( Gluon emmission from a final state quark with an initial state squark as spectator ): described by Vqsingular and Vqfinite.



 ! ************* x-independant part (proportional to the delta function delta(1-x) ******************   
 ! In this part we always have x=1 

      qbar = uman - m3**2 - m2**2     ! this is qbar at point x = 1
       
      v = dsqrt(qbar**2 - 4d0*m2**2*m3**2)/(-qbar)
    
      aie= 2d0*v/(1d0 + 2d0*(-m3**2/qbar)  - v)
      bfe= 2d0*v/(1d0 - v)

! ========================================================================================================================================================

! ========================================================================================================================================================

      if(EpsPole.eq.-2d0) then       ! deals with double poles 

!     *singular part

      Viesingular = 0d0     ! Squark dipole ! We have no collinear singularity in this case -> no double pole 

      Vfesingular = 0d0     ! Quark dipole

      NeuQ2qx_M2DipoleEW = NeuQ2qx_M2Tree()*CF*GS**2/(8d0*pi**2)*(Viesingular + Vfesingular)          !  alphaS = GS**2/(4*pi)

      elseif(EpsPole.eq.-1d0) then

!    *singular part
     
      Viesingular = 1d0 - 1d0/v*dlog((1d0+bfe)/(1d0+aie))  
                                                            
      Vfesingular = 1d0 - 1d0/v*dlog(1d0+aie)

      NeuQ2qx_M2DipoleEW = NeuQ2qx_M2Tree()*CF*GS**2/(8d0*pi**2)*(Viesingular + Vfesingular)  

      elseif(EpsPole.eq.0d0) then

!    *singular part

      Viesingular = IRdiv*( 1d0 - 1d0/v*dlog(  (1d0+bfe)/(1d0+aie) ) ) 

      Vfesingular = IRdiv*( 1d0 - 1d0/v*dlog(1d0+aie) )                      

!    * non-singular part

      Viefinite = 1d0/v*((v+dlog((1d0+aie)/(1d0+bfe)))*dlog((muSc*m3)**2/(qbar**2*(1d0-x0)**2))+dlog(1d0+bfe)
     &       +(dlog(1d0+bfe)**2-dlog(1d0+aie)**2)/2d0+2*dble(cspen(dcmplx(-bfe)))-2*dble(cspen(dcmplx(-aie))))

      Vfefinite = 2d0/v*((v-dlog(aie+1d0))*dlog((m3*muSc)**2/(qbar**2*(1d0-x0)**2))/2d0+v+dlog(1d0+aie)**2/4d0+dble(cspen(dcmplx(-aie))))


      NeuQ2qx_M2DipoleEW = NeuQ2qx_M2Tree()*CF*GS**2/(8d0*pi**2)*(Viesingular + Vfesingular +  Viefinite + Vfefinite) 

      else
      

      NeuQ2qx_M2DipoleEW = 0d0

      endif
                                                                                              

      end



      double precision function NeuQ2qx_M2DipoleG()

      implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

      double precision Isingular,Ifinite,dipcolor,s12,rhox0t,Ax0t,Bx0t,etaAt,etaJu
      double precision  NeuQ2qx_M2Tree,vu,vt, ubar,tbar,Nf,Au,Bu
      double complex cspen
      double precision Vsymm,VQasymm,GammaQ,Kq,VGasymm,gamG,GammaG,Kg


! Mandelstams: t=(pQuark-pNeutralino)^2=(pSquark-pGluon)^2,  u=(pQuark-pSquark)^2=(pGluon-pNeutralino)^2
! The following two variables correspond to \bar{Q}^2
      ubar = uman - m3**2 - m2**2     
      tbar = tman - m2**2

! number of massless quarks (relevant for g-> q+qbar splitting
      Nf=4d0   

! Definition of some other useful variables
      vu = dsqrt(ubar**2 - 4*m2**2*m3**2)/(-ubar) 
      vt = dsqrt(tbar**2 - 4*m2**2*m3**2)/(-tbar) 
      etaAt = m2**2/(-tbar)
      etaJu = m3**2/(-ubar)
    
      Au = 2*vu/(1 + 2*etaJu - vu)
      Bu = 2*vu/(1-vu)

      gamG = 11/6d0*CV-2/3d0*TF*Nf                               !\gamma_g  , Catani-Seymour (short) eq. 5.91  
      Kg = CV*(67/18d0-Pi**2/6d0)-10/9d0*TF*Nf                   !Catani-Seymour (short) eq. 6.17
      Kq = CF*(7/2d0-Pi**2/6d0)
      GammaG = IRdiv*gamG - CV/6d0                                ! hep-ph/0011222  eq. 20   
      
      s12 = sman - m3**2

      rhox0t = dsqrt(1+4*etaAt*(x0-1)*x0)
      Ax0t = rhox0t*(2*x0-3-rhox0t)/(2*(1-x0)*((1-etaAt)*x0-2))
      Bx0t = (rhox0t**2+rhox0t)/(2*(1-x0)*etaAt*x0)

      if(EpsPole.eq.-2d0) then       ! gives double poles 

! FE-IS (from hep-ph/0201036v1 Section 6.2)
      Isingular = (CV*gs**2)/(16*Pi**2)
! FE-IS
      Isingular = Isingular + (CV*gs**2)/(16*Pi**2)

      NeuQ2qx_M2DipoleG =  NeuQ2qx_M2Tree()*Isingular

      elseif(EpsPole.eq.-1d0) then ! gives simple poles 

! --------------------------------------------- FE-FS  (from hep-ph/0201036v1 Section 6.2)-----------------------------------------------------
      Isingular = (gs**2*(CV + gamG + CV*(dlog(m3**2/s12) + dlog(muSc**2/s12))))/(16*Pi**2)

! --------------------------------------------- FE-IS -------------------------------------------------------------------------------------------
      Isingular = Isingular + (gs**2*(gamG +((2*CF - CV)*(vu - dlog(1 + Au)))/vu + CV*dlog(1 + etaAt) + CV*dlog(-(muSc**2/tbar))))/(16*Pi**2)

! --------------------------------------------- IE-FS -------------------------------------------------------------------------------------------
      Isingular = Isingular + (gs**2*(((2*CF - CV)*(vu - dlog((1 + Bu)/(1 + Au))))/vu - CV*(-1 + dlog((1 + etaAt)/etaAt))))/(16*Pi**2)

      NeuQ2qx_M2DipoleG =  NeuQ2qx_M2Tree()*Isingular

      elseif(EpsPole.eq.0d0) then! non-singular part

! --------------------------------------------- FE-FS  (from hep-ph/0201036v1 Section 5.1.3)-----------------------------------------------------
      Vsymm = (IRdiv**2 + IRdiv*dlog(muSc**2/s12) + dlog(muSc**2/s12)**2/2d0 - pi**2/6d0)/2d0   !aber woher kommt der  dlog(s12)-Term->aus Entwicklung des Vorfaktors nach epsilon
     &             + (IRdiv + dlog(muSc**2/s12))/2d0*dlog(m3**2/s12)- dlog(m3**2/s12)**2/4d0 - pi**2/12d0 
     &             - dlog(m3**2/s12)*dlog(s12/sqrtS**2)/2d0 - dlog(m3**2/sqrtS**2)*dlog(s12/sqrtS**2)/2d0         !Vergl. Catani-Seymour (kurz) eq. 6.20 

! (1) emission off quark (q->q+g) (from ChiChi2QQ_M2Dipole.F)
      dipcolor = -CV/(2*CF)

      VQasymm = 3d0/2d0*dlog(s12/sman) + pi**2/6d0 - dble(cspen(dcmplx(s12/sman))) - 2*dlog(s12/sman) - m3**2/s12*dlog(m3**2/sman) !Catani-Seymour (short) eq. 6.22

      GammaQ = CF*(IRdiv - 1d0/2d0*dlog(muSc**2) + 1d0/2d0*dlog(m3**2) - 2d0)

      Ifinite = -gs**2*dipcolor/(8*pi**2)*(CF*(Vsymm + VQasymm - pi**2/3d0) + GammaQ + CF*3/2d0*dlog(muSc**2/s12) + CF*3/2d0 + Kq)  !Vergl. Catani-Seymour (kurz) eq.6.16

! (2) emission off gluon (g->g+g) and (3) gluon decay into massless qqbar pair (g->q+qbar)
      dipcolor =  -1/2d0
                                                                              
      VGasymm = gamG/CV*(dlog(s12/sman)-2*dlog((sqrtS-m3)/sqrtS)-2*m3/(sqrtS+m3)) + pi**2/6d0 -  dble(cspen(dcmplx(s12/sman))) 
     &           +(kappaDip-2/3d0)*m3**2/s12*((2*TF/CV*Nf - 1)*dlog(2*m3/(sqrtS+m3)))                                               !Catani-Seymour (short) eq. 6.24                                                              

      Ifinite =  Ifinite-gs**2*dipcolor/(8*pi**2)*(CV*(Vsymm + VGasymm -pi**2/3d0) + GammaG + gamG*dlog(muSc**2/s12) + gamG + Kg)    !Vergl. Catani-Seymour (kurz) Glg.6.16
! --------------------------------------------- FE-IS -------------------------------------------------------------------------------------------
! emission off quark (q->q+g) 
      Ifinite = Ifinite + ((2*CF - CV)*gs**2*IRdiv*(vu - dlog(1 + Au)))/(16*Pi**2*vu) + 
     &  ((2*CF - CV)*gs**2*(4*dble(cspen(dcmplx(-Au))) + dlog(1 + Au)**2- 2*dlog(1 + Au)*(dlog(-(muSc**2/uBar)) + dlog(etaJu/(-1 + x0)**2)) + 
     &       2*vu*(2 + dlog(-(muSc**2/uBar)) + dlog(etaJu/(-1 + x0)**2))))/(32*Pi**2*vu)
 
! emission off gluon (g->g+g) and gluon decay into massless qqbar pair (g->q+qbar) 
      Ifinite = Ifinite +  (CV*gs**2*IRdiv**2)/(16*Pi**2) + (CV*gs**2*IRdiv*(dlog(1 + etaAt) + dlog(-(muSc**2/tBar))))/(16*Pi**2) + 
     &  (gs**2*(6*GammaG + 6*dlog(muSc**2/(-tbar))*gamG + 6*Kg - CV*Pi**2 + 12*CV*dble(cspen(dcmplx(-Ax0t))) + 3*CV*dlog(1 + etaAt)**2 + 6*CV*dlog(1 + etaAt)*dlog(-(muSc**2/tBar)) + 
     &       3*CV*dlog(-(muSc**2/tBar))**2 + (-11*CV + 4*Nf*TF)*dlog(1- x0)))/(96*Pi**2)

! --------------------------------------------- IE-FS -------------------------------------------------------------------------------------------
! gluon is spectator  (\tilde{q}-> \tilde{q}+g)
      Ifinite = Ifinite  - (CV*gs**2*IRdiv*(-1 + dlog((1 + etaAt)/etaAt)))/(16*Pi**2) + 
     &  (CV*gs**2*(-2*dble(cspen(dcmplx(-Ax0t))) + 2*dble(cspen(dcmplx(-Bx0t))) + 2*dlog(1d0 + Bx0t) + dlog(etaAt) + dlog(etaAt)**2/2d0 - 
     &       dlog(1 + etaAt)**2/2d0 + (1 - dlog((1 + etaAt)/etaAt))*dlog(-(muSc**2/tBar))))/(16*Pi**2)

! quark is spectator (\tilde{q}-> \tilde{q}+g) 
      Ifinite = Ifinite +  ((2*CF - CV)*gs**2*IRdiv*(vu - dlog((1 + Bu)/(1 + Au))))/(16*Pi**2*vu) + 
     &  (CF*(2 - CV/CF)*gs**2*((1 - dlog((1 + Bu)/(1 + Au))/vu)*dlog(-(muSc**2/uBar)) + 
     &       (-2*dble(cspen(dcmplx(-Au))) + 2*dble(cspen(dcmplx(-Bu)))+ dlog(1 + Bu) + (-dlog(1 + Au)**2 + dlog(1 + Bu)**2)/2d0 + 
     &          (vu + dlog((1 + Au)/(1 + Bu)))*dlog(etaJu/(1 - x0)**2))/vu))/(16*Pi**2)

      NeuQ2qx_M2DipoleG =  NeuQ2qx_M2Tree()*Ifinite

     
      else
      
      NeuQ2qx_M2DipoleG = 0d0

      endif
                                                                                              

      end
