! =======================================================================
	  double precision function NeuQ2qx_M2PropagatorCounterV()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "NeuQ2qx_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"


        double precision MssV, MstV, MsuV, MtsV, MttV, MtuV, MusV, MutV, MuuV
	double precision propCTSS, propCTST, propCTSU, propCTTS, propCTTT, propCTTU

        double complex CvL,CvR,CsL,CsR
        double complex PvL,PvR,PsL,PsR
        double complex A1vL,A1vR,A1sL,A1sR,A2vL,A2vR,A2sL,A2sR,A3vL,A3vR,A3sL,A3sR
        double complex Bsm1L,Bsm1R,Bsm2L,Bsm2R,BsgL,BsgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,BvgL,BvgR
        double complex CLoop1,CLoop2
        double complex DLLoop,DRLoop
        double complex A1Lc,A1Rc,A2Lc,A2Rc,A3Lc,A3Rc,BLc,BRc,CLc,DLc,DRc
        double precision mq1,mq1c,mq2,mq2c,mq3,mq3c,wquark
        
        double precision xv, colorfactor 

        integer sumUc, sumTc
        integer sumUend, sumUcend
        integer sumT1, sumT2

        sumUc = 0
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

	NeuQ2qx_M2PropagatorCounterV = 0d0
        propCTSS = 0d0
        propCTST = 0d0
        propCTSU = 0d0
        propCTTS = 0d0
        propCTTT = 0d0
        propCTTU = 0d0

! S channel propagator counter-term

        mq1 = Mf(pstype,psgen) 
c        wquark = width(pstype)
        wquark = 0d0

        call NeuQ2qx_CalcPropQuark(pstype,psgen)

        CsL = PropSL(pstype,psgen)/CF
        CsR = PropSR(pstype,psgen)/CF
        CvL = PropVL(pstype,psgen)/CF
        CvR = PropVR(pstype,psgen)/CF
        
        PsL = (CsL*Mf(pstype,psgen)**2 + (CsR + (CvL + CvR)*Mf(pstype,psgen))*sman)/(Mf(pstype,psgen)**2 - sman - cI*mq1*wquark)
        PsR = (CsR*Mf(pstype,psgen)**2 + (CsL + (CvL + CvR)*Mf(pstype,psgen))*sman)/(Mf(pstype,psgen)**2 - sman - cI*mq1*wquark)
        PvL = (Mf(pstype,psgen)*(CsL + CsR + CvL*Mf(pstype,psgen)) + CvR*sman)/(Mf(pstype,psgen)**2 - sman - cI*mq1*wquark)
        PvR = (Mf(pstype,psgen)*(CsL + CsR + CvR*Mf(pstype,psgen)) + CvL*sman)/(Mf(pstype,psgen)**2 - sman - cI*mq1*wquark)

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

! --- Mss  ---
        
        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
 
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)
        
        propCTSS = colorfactor*CF*NC*2d0*MssV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A1Lc,A1Rc,BLc,BRc,mq1,mq1c,xv,wquark)

        NeuQ2qx_M2PropagatorCounterV = NeuQ2qx_M2PropagatorCounterV + propCTSS
        

! --- Mst  ---

        do sumTc = 0, 1

        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
        propCTST = colorfactor*CF*NC*2*MstV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A2Lc,A2Rc,CLc,mq1,mq2c,xv,wquark)
        
        NeuQ2qx_M2PropagatorCounterV = NeuQ2qx_M2PropagatorCounterV + propCTST

        enddo

! --- Msu  ---

        do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

        propCTSU = colorfactor*CF*NC*2*MsuV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A3Lc,A3Rc,DLc,DRc,mq1,mq3c,xv,wquark)

        NeuQ2qx_M2PropagatorCounterV = NeuQ2qx_M2PropagatorCounterV + propCTSU
                                        
        enddo

! T channel propagator counter-term

        do sumT1 = 0, 1
                do sumT2 = 0, 1
 
        call NeuQ2qx_CalcPropSquark(tman,pttype,ptgen)
        
        A2sL = AL(chimeigen,ptmeigen+sumT1,pttype,ptgen)
        A2sR = AR(chimeigen,ptmeigen+sumT1,pttype,ptgen)
        A2vL = 0d0
        A2vR = 0d0

        CLoop1 = CL(ptmeigen+sumT2,sfmeigen,sftype,sfgen,xtype)
        CLoop2 = 0d0

        mq2 = MSf(ptmeigen+sumT1,pttype,ptgen) 
        
! --- Mts  ---

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
 
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0

        propCTTS = (Prop(ptmeigen+sumT1,ptmeigen+sumT2,pttype,ptgen)/(MSf(ptmeigen+sumT2,pttype,ptgen)**2 - tman))*colorfactor*NC*2*MtsV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A1Lc,A1Rc,BLc,BRc,mq2,mq1c,xv,wquark)
        
        NeuQ2qx_M2PropagatorCounterV = NeuQ2qx_M2PropagatorCounterV + propCTTS

! --- Mtt  ---

                        do sumTc = 0, 1

        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
        propCTTT = (Prop(ptmeigen+sumT1,ptmeigen+sumT2,pttype,ptgen)/(MSf(ptmeigen+sumT2,pttype,ptgen)**2 - tman))*colorfactor*NC*2*MttV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A2Lc,A2Rc,CLc,mq2,mq2c,xv)
       
        NeuQ2qx_M2PropagatorCounterV = NeuQ2qx_M2PropagatorCounterV + propCTTT

                        enddo

! --- Mtu  ---

                        do sumUc = 0,sumUcend
	
        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

        propCTTU = (Prop(ptmeigen+sumT1,ptmeigen+sumT2,pttype,ptgen)/(MSf(ptmeigen+sumT2,pttype,ptgen)**2 - tman))*colorfactor*NC*2*MtuV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A3Lc,A3Rc,DLc,DRc,mq2,mq3c,xv)

        NeuQ2qx_M2PropagatorCounterV = NeuQ2qx_M2PropagatorCounterV + propCTTU
                                        
                        enddo
                
                enddo
        enddo
        
!         NeuQ2qx_M2PropagatorCounterV =  propCTSS

	  end
