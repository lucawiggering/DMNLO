************************************************
        double precision function NeuQ2qx_M2TreeVgen()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double precision MssV, MstV, MsuV, MtsV, MttV, MtuV, MusV, MutV, MuuV
        double precision treeSS, treeTT, treeUU, treeST, treeTS, treeSU, treeUS, treeTU, treeUT
        
        double complex PvL,PvR,PsL,PsR
        double complex A1vL,A1vR,A1sL,A1sR,A2vL,A2vR,A2sL,A2sR,A3vL,A3vR,A3sL,A3sR
        double complex Bsm1L,Bsm1R,Bsm2L,Bsm2R,BsgL,BsgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,BvgL,BvgR
        double complex CLoop1,CLoop2
        double complex DLLoop,DRLoop
        double complex A1Lc,A1Rc,A2Lc,A2Rc,A3Lc,A3Rc,BLc,BRc,CLc,DLc,DRc
        double precision mq1,mq1c,mq2,mq2c,mq3,mq3c,wquark
        
        double precision xv, colorfactor 

        integer sumU, sumUc, sumT, sumTc
        integer sumUend, sumUcend, sumTend, sumTcend

        sumU = 0
        sumUc = 0
        sumT = 0
        sumTc = 0

        if (xtype.eq.1) then
        
           xv = 0d0
           colorfactor = 4./3.

           sumTend = 0
           sumTcend = 0
           sumUend = -1
           sumUcend = -1

        elseif (xtype.eq.2) then
        
           xv = 0d0
           colorfactor = 1d0
           
           sumTend = 1
           sumTcend = 1
           sumUend = -1
           sumUcend = -1
           
        elseif (xtype.eq.3) then
        
c Feynman Gauge
c           xv = 0d0
c Unitarity Gauge
           xv = 1/(m4**2)

           colorfactor = 1d0
        
           sumTend = 1
           sumTcend = 1
           sumUend = 3
           sumUcend = 3
           
        elseif (xtype.eq.4) then

c Feynman Gauge
c           xv = 0d0
c Unitarity Gauge
           xv = 1/(m4**2)
           colorfactor = 1d0
               
           sumTend = 1
           sumTcend = 1
           sumUend = 1
           sumUcend = 1
           
         endif

    
        NeuQ2qx_M2TreeVgen = 0d0
        treeSS = 0d0
        treeTT = 0d0
        treeUU = 0d0
        treeST = 0d0
        treeTS = 0d0
        treeSU = 0d0
        treeUS = 0d0
        treeUT = 0d0
        treeTU = 0d0
        
        

! S channel General Amplitude set to Tree Level

        PsL = Mf(pstype,psgen)
        PsR = Mf(pstype,psgen)
        PvL = 1d0
        PvR = 1d0        
        
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

        mq1 = Mf(pstype,psgen) 
!         wquark = width(pstype)
        wquark = 0d0

! --- Mss  ---
        
        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
 
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)
        
        treeSS = colorfactor*NC*MssV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A1Lc,A1Rc,BLc,BRc,mq1,mq1c,xv,wquark)

        NeuQ2qx_M2TreeVgen = NeuQ2qx_M2TreeVgen + treeSS


! --- Mst  ---

        do sumTc = 0, sumTcend

        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
        mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
        treeST = treeST + colorfactor*NC*MstV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A2Lc,A2Rc,CLc,mq1,mq2c,xv,wquark)

        enddo
        
        NeuQ2qx_M2TreeVgen = NeuQ2qx_M2TreeVgen + treeST

! --- Msu  ---

        do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

        treeSU = treeSU + colorfactor*NC*MsuV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A3Lc,A3Rc,DLc,DRc,mq1,mq3c,xv,wquark)
                                        
        enddo
        
        NeuQ2qx_M2TreeVgen = NeuQ2qx_M2TreeVgen + treeSU

        
        
        
! T channel General Amplitude set to Tree Level

! --- Mts  ---
        do sumT = 0, sumTend
        
        A2sL = AL(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2sR = AR(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2vL = 0d0
        A2vR = 0d0
                
        CLoop1 = CL(ptmeigen+sumT,sfmeigen,sftype,sfgen,xtype)
        CLoop2 = 0d0

        mq2 = MSf(ptmeigen+sumT,pttype,ptgen)

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
 
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0
        
        treeTS = treeTS + colorfactor*NC*MtsV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A1Lc,A1Rc,BLc,BRc,mq2,mq1c,xv,wquark)
        
        enddo
        
        NeuQ2qx_M2TreeVgen = NeuQ2qx_M2TreeVgen + treeTS

! --- Mtt  ---

        do sumT = 0, sumTend
              do sumTc = 0, sumTcend
        
        A2sL = AL(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2sR = AR(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2vL = 0d0
        A2vR = 0d0
                
        CLoop1 = CL(ptmeigen+sumT,sfmeigen,sftype,sfgen,xtype)
        CLoop2 = 0d0

        mq2 = MSf(ptmeigen+sumT,pttype,ptgen)

        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 

        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))

        mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
        treeTT = treeTT + colorfactor*NC*MttV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A2Lc,A2Rc,CLc,mq2,mq2c,xv)

                enddo
        enddo
                
        NeuQ2qx_M2TreeVgen = NeuQ2qx_M2TreeVgen + treeTT
        
        
! --- Mtu  ---

        do sumT = 0, sumTend
                do sumUc = 0,sumUcend
        
        A2sL = AL(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2sR = AR(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2vL = 0d0
        A2vR = 0d0
                
        CLoop1 = CL(ptmeigen+sumT,sfmeigen,sftype,sfgen,xtype)
        CLoop2 = 0d0

        mq2 = MSf(ptmeigen+sumT,pttype,ptgen)
        
        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

        treeTU = treeTU + colorfactor*NC*MtuV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A3Lc,A3Rc,DLc,DRc,mq2,mq3c,xv)
                                        
                                        
                enddo
        enddo
                
        NeuQ2qx_M2TreeVgen = NeuQ2qx_M2TreeVgen + treeTU

! U channel General Amplitude set to Tree Level

! --- Mus  ---

        do sumU = 0,sumUend
        
        A3sL = AL(pumeigen+sumU,sfmeigen,sftype,sfgen)
        A3sR = AR(pumeigen+sumU,sfmeigen,sftype,sfgen)
        A3vL = 0d0
        A3vR = 0d0
             
        DLLoop = DL(chimeigen,pumeigen + sumU,sftype,xtype)
        DRLoop = DR(chimeigen,pumeigen + sumU,sftype,xtype)     

        mq3 = MNeuCha(pumeigen+sumU)

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
 
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0
        
        treeUS = treeUS + colorfactor*NC*MusV(A3sL,A3sR,A3vL,A3vR,DLLoop,DRLoop,A1Lc,A1Rc,BLc,BRc,mq3,mq1c,xv,wquark)

        enddo
        
        NeuQ2qx_M2TreeVgen = NeuQ2qx_M2TreeVgen + treeUS

! --- Mut  ---

        do sumU = 0,sumUend

                do sumTc = 0, sumTcend

        A3sL = AL(pumeigen+sumU,sfmeigen,sftype,sfgen)
        A3sR = AR(pumeigen+sumU,sfmeigen,sftype,sfgen)
        A3vL = 0d0
        A3vR = 0d0
             
        DLLoop = DL(chimeigen,pumeigen + sumU,sftype,xtype)
        DRLoop = DR(chimeigen,pumeigen + sumU,sftype,xtype)     

        mq3 = MNeuCha(pumeigen+sumU)        
                
        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
        mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
        treeUT = treeUT + colorfactor*NC*MutV(A3sL,A3sR,A3vL,A3vR,DLLoop,DRLoop,A2Lc,A2Rc,CLc,mq3,mq2c,xv)
        
             enddo
        enddo
        
        NeuQ2qx_M2TreeVgen = NeuQ2qx_M2TreeVgen + treeUT

! --- Muu  ---

        do sumU = 0,sumUend
                do sumUc = 0,sumUcend

        A3sL = AL(pumeigen+sumU,sfmeigen,sftype,sfgen)
        A3sR = AR(pumeigen+sumU,sfmeigen,sftype,sfgen)
        A3vL = 0d0
        A3vR = 0d0
             
        DLLoop = DL(chimeigen,pumeigen + sumU,sftype,xtype)
        DRLoop = DR(chimeigen,pumeigen + sumU,sftype,xtype)     

        mq3 = MNeuCha(pumeigen+sumU)        
                
        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

        treeUU = treeUU + colorfactor*NC*MuuV(A3sL,A3sR,A3vL,A3vR,DLLoop,DRLoop,A3Lc,A3Rc,DLc,DRc,mq3,mq3c,xv)
                                        
                enddo

        enddo
        NeuQ2qx_M2TreeVgen = NeuQ2qx_M2TreeVgen + treeUU
        
        
        
        treeperc(1) = treeSS
        treeperc(2) = treeTT
        treeperc(3) = treeUU
        treeperc(4) = treeST + treeTS
        treeperc(5) = treeSU + treeUS
        treeperc(6) = treeTU + treeUT
        treeperc(7) = treeperc(1) + treeperc(2) + treeperc(3) + treeperc(4) + treeperc(5) + treeperc(6)
        
        
        NeuQ2qx_M2TreeVgen = treeSS + treeTT + treeTS + treeST + treeUU + treeSU + treeUS + treeTU + treeUT

        
        end