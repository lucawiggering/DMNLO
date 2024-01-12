! =======================================================================
	  double precision function NeuQ2qx_M2VertexCounterV()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "NeuQ2qx_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"

        double precision MssV, MstV, MsuV, MtsV, MttV, MtuV, MusV, MutV, MuuV
    	double precision vertCTSS, vertCTST, vertCTSU, vertCTTS, vertCTTT, vertCTTU, vertCTUS, vertCTUT, vertCTUU

        double complex PvL,PvR,PsL,PsR
        double complex A1vL,A1vR,A1sL,A1sR,A2vL,A2vR,A2sL,A2sR,A3vL,A3vR,A3sL,A3sR
        double complex Bsm1L,Bsm1R,Bsm2L,Bsm2R,BsgL,BsgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,BvgL,BvgR
        double complex CLoop1,CLoop2
        double complex DLLoop,DRLoop
        double complex A1Lc,A1Rc,A2Lc,A2Rc,A3Lc,A3Rc,BLc,BRc,CLc,DLc,DRc
        double precision mq1,mq1c,mq2,mq2c,mq3,mq3c,wquark
        
        double precision xv, colorfactor 

        integer sumU, sumUc, sumT, sumTc
        integer sumUend, sumUcend

        sumU = 0
        sumUc = 0
        sumT = 0
        sumTc = 0

        if (xtype.eq.2) then
        
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

    
    	NeuQ2qx_M2VertexCounterV = 0d0
        vertCTSS = 0d0
        vertCTST = 0d0
        vertCTSU = 0d0
        vertCTTS = 0d0
        vertCTTT = 0d0
        vertCTTU = 0d0
        vertCTUS = 0d0
        vertCTUT = 0d0
        vertCTUU = 0d0


! ----- Counter-term for Neutralino - Sfermion - Fermion Coupling in S channel (A1 Coupling) -----

! Set propagator to tree-level

        PsL = Mf(pstype,psgen)
        PsR = Mf(pstype,psgen)
        PvL = 1d0
        PvR = 1d0        

! Set A coupling to counter-term

c        A1sL = dAL(chimeigen,sfmeigen,sftype,sfgen)/CF
c        A1sR = dAR(chimeigen,sfmeigen,sftype,sfgen)/CF
        A1sL = dAL(chimeigen,sfmeigen,sftype,sfgen)
        A1sR = dAR(chimeigen,sfmeigen,sftype,sfgen)
        A1vL = 0d0
        A1vR = 0d0
        
! Set B coupling to tree-level
        
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

        mq1 = Mf(pstype,psgen) 
c        wquark = width(pstype)
        wquark = 0D0

! --- Mss  ---
        
        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
 
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)
        
c        vertCTSS = colorfactor*CF*NC*2*MssV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A1Lc,A1Rc,BLc,BRc,mq1,mq1c,xv,wquark)
        vertCTSS = NC*2*MssV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A1Lc,A1Rc,BLc,BRc,mq1,mq1c,xv,wquark)

        NeuQ2qx_M2VertexCounterV = NeuQ2qx_M2VertexCounterV + vertCTSS

! --- Mst  ---

        do sumTc = 0, 1

        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
c        vertCTST = colorfactor*CF*NC*2*MstV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A2Lc,A2Rc,CLc,mq1,mq2c,xv,wquark)
        vertCTST = NC*2*MstV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A2Lc,A2Rc,CLc,mq1,mq2c,xv,wquark)

        NeuQ2qx_M2VertexCounterV = NeuQ2qx_M2VertexCounterV + vertCTST

        enddo

! --- Msu  ---

        do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

c        vertCTSU = colorfactor*CF*NC*2*MsuV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A3Lc,A3Rc,DLc,DRc,mq1,mq3c,xv,wquark)
        vertCTSU = NC*2*MsuV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A3Lc,A3Rc,DLc,DRc,mq1,mq3c,xv,wquark)
        
        NeuQ2qx_M2VertexCounterV = NeuQ2qx_M2VertexCounterV + vertCTSU
                                        
        enddo

! ----- Counter-term for Neutralino - Sfermion - Fermion Coupling in T channel (A2 Coupling) -----

        do sumT = 0,1
        
! Set A coupling to counter-term

c        A2sL = dAL(chimeigen,ptmeigen+sumT,pttype,ptgen)/CF
c        A2sR = dAR(chimeigen,ptmeigen+sumT,pttype,ptgen)/CF
        A2sL = dAL(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2sR = dAR(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2vL = 0d0
        A2vR = 0d0

! Set C coupling to tree-level

        CLoop1 = CL(ptmeigen+sumT,sfmeigen,sftype,sfgen,xtype)        
        CLoop2 = 0d0

        mq2 = MSf(ptmeigen+sumT,pttype,ptgen)

! --- Mts  ---

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
 
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
         wquark = 0D0
       
c        vertCTTS = colorfactor*CF*NC*2*MtsV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A1Lc,A1Rc,BLc,BRc,mq2,mq1c,xv,wquark)
        vertCTTS = NC*2*MtsV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A1Lc,A1Rc,BLc,BRc,mq2,mq1c,xv,wquark)
        
        NeuQ2qx_M2VertexCounterV = NeuQ2qx_M2VertexCounterV + vertCTTS

! --- Mtt  ---

                do sumTc = 0, 1

        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
c        vertCTTT = colorfactor*CF*NC*2*MttV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A2Lc,A2Rc,CLc,mq2,mq2c,xv)
        vertCTTT = NC*2*MttV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A2Lc,A2Rc,CLc,mq2,mq2c,xv)
        
        NeuQ2qx_M2VertexCounterV = NeuQ2qx_M2VertexCounterV + vertCTTT

                enddo

! --- Mtu  ---

                do sumUc = 0,sumUcend
	
        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

c        vertCTTU = colorfactor*CF*NC*2*MtuV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A3Lc,A3Rc,DLc,DRc,mq2,mq3c,xv)
        vertCTTU = NC*2*MtuV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A3Lc,A3Rc,DLc,DRc,mq2,mq3c,xv)

        NeuQ2qx_M2VertexCounterV = NeuQ2qx_M2VertexCounterV + vertCTTU
                                        
                enddo

        enddo

! ----- Counter-term for Neutralino - Sfermion - Fermion Coupling in U channel (A3 Coupling) -----

        do sumU = 0,sumUend

! Set A coupling to counter-term
        
c        A3sL = dAL(pumeigen+sumU,sfmeigen,sftype,sfgen)/CF
c        A3sR = dAR(pumeigen+sumU,sfmeigen,sftype,sfgen)/CF
        A3sL = dAL(pumeigen+sumU,sfmeigen,sftype,sfgen)
        A3sR = dAR(pumeigen+sumU,sfmeigen,sftype,sfgen)
        A3vL = 0d0
        A3vR = 0d0

! Set D coupling to tree-level
             
        DLLoop = DL(chimeigen,pumeigen + sumU,sftype,xtype)
        DRLoop = DR(chimeigen,pumeigen + sumU,sftype,xtype)     

        mq3 = MNeuCha(pumeigen+sumU)

! --- Mus  ---

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
 
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
         wquark = 0D0

c        vertCTUS = colorfactor*CF*NC*2*MusV(A3sL,A3sR,A3vL,A3vR,DLLoop,DRLoop,A1Lc,A1Rc,BLc,BRc,mq3,mq1c,xv,wquark)
        vertCTUS = NC*2*MusV(A3sL,A3sR,A3vL,A3vR,DLLoop,DRLoop,A1Lc,A1Rc,BLc,BRc,mq3,mq1c,xv,wquark)

        NeuQ2qx_M2VertexCounterV = NeuQ2qx_M2VertexCounterV + vertCTUS

! --- Mut  ---

                do sumTc = 0, 1
                
        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
c        vertCTUT = colorfactor*CF*NC*2*MutV(A3sL,A3sR,A3vL,A3vR,DLLoop,DRLoop,A2Lc,A2Rc,CLc,mq3,mq2c,xv)
        vertCTUT = NC*2*MutV(A3sL,A3sR,A3vL,A3vR,DLLoop,DRLoop,A2Lc,A2Rc,CLc,mq3,mq2c,xv)
        
        NeuQ2qx_M2VertexCounterV = NeuQ2qx_M2VertexCounterV + vertCTUT

                enddo

! --- Muu  ---

                do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

c        vertCTUU = colorfactor*CF*NC*2*MuuV(A3sL,A3sR,A3vL,A3vR,DLLoop,DRLoop,A3Lc,A3Rc,DLc,DRc,mq3,mq3c,xv)
        vertCTUU = NC*2*MuuV(A3sL,A3sR,A3vL,A3vR,DLLoop,DRLoop,A3Lc,A3Rc,DLc,DRc,mq3,mq3c,xv)

        NeuQ2qx_M2VertexCounterV = NeuQ2qx_M2VertexCounterV + vertCTUU
                                        
                enddo

        enddo

! ----- Counter-term for Vector - Fermion - Fermion Coupling in S channel (B Coupling) -----

! Set propagator to tree-level

        PsL = Mf(pstype,psgen)
        PsR = Mf(pstype,psgen)
        PvL = 1d0
        PvR = 1d0        

! Set A coupling to tree-level

        A1sL = AL(chimeigen,sfmeigen,sftype,sfgen)
        A1sR = AR(chimeigen,sfmeigen,sftype,sfgen)
        A1vL = 0d0
        A1vR = 0d0
        
! Set B coupling to counter-term
        
c        BsgL = dBL(xtype,ftype,sfgen)/CF
c        BsgR = dBR(xtype,ftype,sfgen)/CF
        BsgL = dBL(xtype,ftype,sfgen)
        BsgR = dBR(xtype,ftype,sfgen)
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

        mq1 = Mf(pstype,psgen) 
c        wquark = width(pstype)
        wquark = 0D0

! --- Mss  ---
        
        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
 
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)
        
c        vertCTSS = colorfactor*CF*NC*2*MssV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A1Lc,A1Rc,BLc,BRc,mq1,mq1c,xv,wquark)
        vertCTSS = NC*2*MssV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A1Lc,A1Rc,BLc,BRc,mq1,mq1c,xv,wquark)

        NeuQ2qx_M2VertexCounterV = NeuQ2qx_M2VertexCounterV + vertCTSS

! --- Mst  ---

        do sumTc = 0, 1

        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
c        vertCTST = colorfactor*CF*NC*2*MstV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A2Lc,A2Rc,CLc,mq1,mq2c,xv,wquark)
        vertCTST = NC*2*MstV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A2Lc,A2Rc,CLc,mq1,mq2c,xv,wquark)
        
        NeuQ2qx_M2VertexCounterV = NeuQ2qx_M2VertexCounterV + vertCTST

        enddo

! --- Msu  ---

        do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

c        vertCTSU = colorfactor*CF*NC*2*MsuV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A3Lc,A3Rc,DLc,DRc,mq1,mq3c,xv,wquark)
        vertCTSU = NC*2*MsuV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A3Lc,A3Rc,DLc,DRc,mq1,mq3c,xv,wquark)
        
        NeuQ2qx_M2VertexCounterV = NeuQ2qx_M2VertexCounterV + vertCTSU
                                        
        enddo        

! ----- Counter-term for Vector - Sfermion - Sfermion Coupling in T channel (C Coupling) -----

        do sumT = 0,1
        
! Set A coupling to tree-level

        A2sL = AL(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2sR = AR(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2vL = 0d0
        A2vR = 0d0

! Set C coupling to counter-term

c        CLoop1 = dCL(ptmeigen+sumT,sfmeigen,sftype,sfgen,xtype)/CF      
        CLoop1 = dCL(ptmeigen+sumT,sfmeigen,sftype,sfgen,xtype)      
        CLoop2 = 0d0

        mq2 = MSf(ptmeigen+sumT,pttype,ptgen)

! --- Mts  ---

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
 
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0D0
        
c        vertCTTS = colorfactor*CF*NC*2*MtsV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A1Lc,A1Rc,BLc,BRc,mq2,mq1c,xv,wquark)
        vertCTTS = NC*2*MtsV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A1Lc,A1Rc,BLc,BRc,mq2,mq1c,xv,wquark)
        
        NeuQ2qx_M2VertexCounterV = NeuQ2qx_M2VertexCounterV + vertCTTS

! --- Mtt  ---

                do sumTc = 0, 1

        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
c        vertCTTT = colorfactor*CF*NC*2*MttV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A2Lc,A2Rc,CLc,mq2,mq2c,xv)
        vertCTTT = NC*2*MttV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A2Lc,A2Rc,CLc,mq2,mq2c,xv)
        
        NeuQ2qx_M2VertexCounterV = NeuQ2qx_M2VertexCounterV + vertCTTT

                enddo

! --- Mtu  ---

                do sumUc = 0,sumUcend
	
        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

c        vertCTTU = colorfactor*CF*NC*2*MtuV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A3Lc,A3Rc,DLc,DRc,mq2,mq3c,xv)
        vertCTTU = NC*2*MtuV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A3Lc,A3Rc,DLc,DRc,mq2,mq3c,xv)

        NeuQ2qx_M2VertexCounterV = NeuQ2qx_M2VertexCounterV + vertCTTU
                                        
                enddo

        enddo        

	  end
