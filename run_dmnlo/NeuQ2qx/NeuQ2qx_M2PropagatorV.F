************************************************
        double precision function NeuQ2qx_M2PropagatorV()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

        double precision MssV, MstV, MsuV, MtsV, MttV, MtuV, MusV, MutV, MuuV
        double precision propSS, propTT, propUU, propST, propTS, propSU, propUS, propTU, propUT

        double complex PvL,PvR,PsL,PsR
        double complex CvL,CvR,CsL,CsR
        double complex SEsquark

        double complex A1vL,A1vR,A1sL,A1sR,A2sL,A2sR,A2vL,A2vR,A3vL,A3vR,A3sL,A3sR
        double complex Bsm1L,Bsm1R,Bsm2L,Bsm2R,BsgL,BsgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,BvgL,BvgR
        double complex CLoop1,CLoop2
        double complex A1Lc,A1Rc,A2Lc,A2Rc,A3Lc,A3Rc,BLc,BRc,CLc,DLc,DRc
    
        double precision mq1,mq1c,mq2,mq2c,mq3,mq3c,wquark
        
        double precision xv, colorfactor
        integer k 
        double precision Krondelta
        

        double complex g0, g0L, g0R, g1, g1L, g1R, g2, g2L, g2R
        double precision Mloop0, Mloop1, Mloop2

        integer sumU, sumUc, sumT, sumTc
        integer sumUend, sumUcend
        integer sumT1, sumT2
        
        sumU = 0
        sumUc = 0
        sumT = 0
        sumTc = 0
        sumT1 = 0
        sumT2 = 0
                
        if (xtype.eq.1) then
        
           xv = 0d0
           colorfactor = 4d0/3d0
           
           sumUend = -1
           sumUcend = -1
        
        elseif (xtype.eq.2) then
        
           xv = 0d0
           colorfactor = 1d0
           
           sumUend = -1
           sumUcend = -1
           
        elseif (xtype.eq.3) then
        
           xv = 1/(m4**2)
           colorfactor = 1d0
        
           sumUend = 3
           sumUcend = 3
           
        elseif (xtype.eq.4) then
        
           xv = 1/(m4**2)
           colorfactor = 1d0
               
           sumUend = 1
           sumUcend = 1
           
         endif
        
        g0 = 0d0
        g0L = 0d0
        g0R = 0d0
        g1L = 0d0
        g1R = 0d0
        g2 = 0d0
        g2L = 0d0
        g2R = 0d0
        
        Mloop0 = 0d0
        Mloop1 = 0d0
        Mloop2 = 0d0
        
        mq1 = 0d0
        mq1c = 0d0
        mq2 = 0d0
        mq2c = 0d0
        mq3 = 0d0
        mq3c = 0d0
                
        NeuQ2qx_M2PropagatorV = 0d0
        propSS = 0d0
        propTT = 0d0
        propUU = 0d0
        propST = 0d0
        propTS = 0d0
        propSU = 0d0
        propUS = 0d0
        propUT = 0d0
        propTU = 0d0
        
! For comparison with Julia (= higgs) results
c        m4 = Mh0
c        print*, "tman", tman

! ----- Correction to fermion propagator in the S channel and Interference with the Tree-level S, T, U Channels -----

        CsL = 0d0
        CsR = 0d0
        CvL = 0d0
        CvR = 0d0
        PsL = 0d0
        PsR = 0d0
        PvL = 0d0
        PvR = 0d0

! Set A and B couplings to tree-level

        A1sL = AL(chimeigen,sfmeigen,sftype,sfgen)
        A1sR = AR(chimeigen,sfmeigen,sftype,sfgen)
        A1vL = 0d0
        A1vR = 0d0
        
        BsgL = BL(xtype,ftype)
        BsgR = BR(xtype,ftype)
        Bsm1L = 0d0
        Bsm1R = 0d0
        Bsm2L = 0d0
        Bsm2R = 0d0
        BvgL = 0d0
        BvgR = 0d0
        Bvm1L = 0d0
        Bvm1R = 0d0
        Bvm2L = 0d0
        Bvm2R = 0d0
        
c       Expressions in MV have a general propagator of the form ( (PsL.PL + PsR.PR) + (PvL.PL + PvR.PR).slash(p) )/Den
c       which is then expressed as (slash(p)+m).(CsL.PL + CsR.PR + CvL.PL.slash(p) + CvR.PR.slash(p)).(slash(p)+m)/Den^2.
c       First, values are given to CpL, etc. Then, PsL, etc. are expressed in function of CpL, etc. 
c       This basis is chosen so that the tree-level can be recovered by setting PsL = Mf(pstype,psgen), PsR = Mf(pstype,psgen), PvL = 1d0, PvR = 1d0.

! --- Gluon loop ---
        
        g0 = B(1,pstype)
        g1 = B(1,pstype)

        Mloop0 = 0
        Mloop1 = Mf(pstype,psgen) 

        mq1 = Mf(pstype,psgen) 
c        wquark = width(pstype)
        wquark = 0d0

        call myBget(sman,0d0,Mf(pstype,psgen)**2)

        CsL = CsL - (g0*g1*Mloop1*xb0)/(4.*Pi**2)
        CsR = CsR - (g0*g1*Mloop1*xb0)/(4.*Pi**2)
        CvL = CvL + (g0*g1*(xb0 + xb1))/(8.*Pi**2)
        CvR = CvR + (g0*g1*(xb0 + xb1))/(8.*Pi**2)
        
! --- Gluino loop ---

        do k=1,2

        g0L = -sqrt2*GS*USf(k,1,sftype,sfgen)
        g0R =  sqrt2*GS*USf(k,2,sftype,sfgen)
        g1L =  sqrt2*GS*USfC(k,2,sftype,sfgen)
        g1R = -sqrt2*GS*USfC(k,1,sftype,sfgen)       
                
        Mloop0 = MSf(k,sftype,sfgen)
        Mloop1 = MGl

        mq1 = Mf(pstype,psgen) 
c        wquark = width(pstype)
        wquark = 0d0

        call myBget(sman,MSf(k,sftype,sfgen)**2,MGl**2)

        CsL = CsL + (g0L*g1L*Mloop1*xb0)/(16.*Pi**2)
        CsR = CsR + (g0R*g1R*Mloop1*xb0)/(16.*Pi**2)
        CvL = CvL + (g0R*g1L*(xb0 + xb1))/(16.*Pi**2)
        CvR = CvR + (g0L*g1R*(xb0 + xb1))/(16.*Pi**2)
        
        enddo

        PsL = (CsL*Mf(pstype,psgen)**2 + (CsR + (CvL + CvR)*Mf(pstype,psgen))*sman)/(Mf(pstype,psgen)**2 - sman - cI*mq1*wquark)
        PsR = (CsR*Mf(pstype,psgen)**2 + (CsL + (CvL + CvR)*Mf(pstype,psgen))*sman)/(Mf(pstype,psgen)**2 - sman - cI*mq1*wquark)
        PvL = (Mf(pstype,psgen)*(CsL + CsR + CvL*Mf(pstype,psgen)) + CvR*sman)/(Mf(pstype,psgen)**2 - sman - cI*mq1*wquark)
        PvR = (Mf(pstype,psgen)*(CsL + CsR + CvR*Mf(pstype,psgen)) + CvL*sman)/(Mf(pstype,psgen)**2 - sman - cI*mq1*wquark)

! --- Mss  ---            

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))

        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)

        propSS = colorfactor*CF*NC*2d0*MssV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A1Lc,A1Rc,BLc,BRc,mq1,mq1c,xv,wquark)

        NeuQ2qx_M2PropagatorV = NeuQ2qx_M2PropagatorV + propSS

! --- Mst  ---
        
        do sumTc = 0, 1
                
        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
        mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
        propST = colorfactor*CF*NC*2*MstV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A2Lc,A2Rc,CLc,mq1,mq2c,xv,wquark)
        
        NeuQ2qx_M2PropagatorV = NeuQ2qx_M2PropagatorV + propST

        enddo

! --- Msu  ---

        do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

        propSU = colorfactor*CF*NC*2*MsuV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A3Lc,A3Rc,DLc,DRc,mq1,mq3c,xv,wquark)

        NeuQ2qx_M2PropagatorV = NeuQ2qx_M2PropagatorV + propSU
                                                
        enddo
        
! ----- Correction to Sfermion propagator in the T channel and Interference with the Tree-level S, T, U Channels -----

        do sumT1 = 0, 1
                do sumT2 = 0, 1                

        SEsquark = 0d0

! Set A and C couplings to tree-level

        A2sL = AL(chimeigen,ptmeigen+sumT1,pttype,ptgen)
        A2sR = AR(chimeigen,ptmeigen+sumT1,pttype,ptgen)
        A2vL = 0d0
        A2vR = 0d0

        CLoop1 = CL(ptmeigen+sumT2,sfmeigen,sftype,sfgen,xtype)
        CLoop2 = 0d0

        mq2 = MSf(ptmeigen+sumT1,pttype,ptgen) 

! --- Gluon loop ---

        g0 = CL(ptmeigen+sumT2,ptmeigen+sumT2,pttype,ptgen,1)*Krondelta(ptmeigen+sumT1,ptmeigen+sumT2)
        g1 = CL(ptmeigen+sumT1,ptmeigen+sumT1,pttype,ptgen,1)*Krondelta(ptmeigen+sumT1,ptmeigen+sumT2)
        
        Mloop0 = 0
        Mloop1 = MSf(ptmeigen+sumT1,pttype,ptgen) 

        call myAget(MSf(ptmeigen+sumT1,pttype,ptgen)**2)
        call myBget(tman,0d0,MSf(ptmeigen+sumT1,pttype,ptgen)**2)

        SEsquark = SEsquark - g0*g1*(xa0 + Mloop0**2*xb0 + 4*tman*(xb0 + xb1))/(16.*Pi**2*(MSf(ptmeigen+sumT2,pttype,ptgen)**2 - tman))
        
! --- Gluino loop ---

        g0L =  sqrt2*GS*USfC(ptmeigen+sumT2,2,pttype,ptgen)
        g0R = -sqrt2*GS*USfC(ptmeigen+sumT2,1,pttype,ptgen)        
        g1L = -sqrt2*GS*USf(ptmeigen+sumT1,1,pttype,ptgen)
        g1R =  sqrt2*GS*USf(ptmeigen+sumT1,2,pttype,ptgen)
        
        Mloop0 = MGl
        Mloop1 = Mf(pttype,ptgen)

        mq2 = MSf(ptmeigen+sumT1,pttype,ptgen) 

        call myAget(Mf(pttype,ptgen)**2)
        call myBget(tman,MGl**2,Mf(pttype,ptgen)**2)

        SEsquark = SEsquark - (-2*(g0L*g1L + g0R*g1R)*Mloop0*Mloop1*xb0 - 2*(g0R*g1L + g0L*g1R)*(xa0 + Mloop0**2*xb0 + 
     &               tman*xb1))/(16.*Pi**2*(-MSf(ptmeigen+sumT2,pttype,ptgen)**2 + tman))

! --- Squark loop ---

        do k=1,2

        mq2 = MSf(ptmeigen+sumT1,pttype,ptgen) 

        call myAget(MSf(k,pttype,ptgen)**2)

c        g0 = - GS**2*(Krondelta(ptmeigen+sumT2,k) - 2*USfC(1,ptmeigen+sumT2,pttype,ptgen)*USf(1,k,pttype,ptgen))*(Krondelta(k,ptmeigen+sumT1) -  
c     & 2*USfC(1,k,pttype,ptgen)*USf(1,ptmeigen+sumT1,pttype,ptgen))

        g0 = S(ptmeigen+sumT2,k,k,ptmeigen+sumT1,pttype,pttype,ptgen)
        
        SEsquark = SEsquark - g0*xa0/((16.*Pi**2)*(MSf(ptmeigen+sumT2,pttype,ptgen)**2 - tman))

        enddo

! --- Mts  ---            

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))

        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0

        propTS = SEsquark*colorfactor*CF*NC*2*MtsV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A1Lc,A1Rc,BLc,BRc,mq2,mq1c,xv,wquark)

        NeuQ2qx_M2PropagatorV = NeuQ2qx_M2PropagatorV + propTS        

! --- Mtt  ---
        
                        do sumTc = 0, 1
                
        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
        mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
        propTT = SEsquark*colorfactor*CF*NC*2*MttV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A2Lc,A2Rc,CLc,mq2,mq2c,xv)
        
        NeuQ2qx_M2PropagatorV = NeuQ2qx_M2PropagatorV + propTT

                        enddo

! --- Mtu  ---

                        do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

        propTU = SEsquark*colorfactor*CF*NC*2*MtuV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A3Lc,A3Rc,DLc,DRc,mq2,mq3c,xv)

        NeuQ2qx_M2PropagatorV = NeuQ2qx_M2PropagatorV + propTU
                                        
                        enddo

                enddo
        enddo

        end
