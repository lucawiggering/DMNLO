************************************************
	double precision function NeuQ2qx_M2VertexV()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

        double precision MssV, MstV, MsuV, MtsV, MttV, MtuV, MusV, MutV, MuuV
	double precision vertSS, vertTT, vertUU, vertST, vertTS, vertSU, vertUS, vertTU, vertUT

        double complex PvL,PvR,PsL,PsR
        
        double complex A1sL,A1sR,A1vL,A1vR,A2sL,A2sR,A2vL,A2vR,A3sL,A3sR,A3vL,A3vR
        double complex BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R
        double complex CLoop1,CLoop2
        double complex DLLoop,DRLoop
        double complex A1Lc,A1Rc,A2Lc,A2Rc,A3Lc,A3Rc,BLc,BRc,CLc,DLc,DRc
        
        double precision mq1,mq1c,mq2,mq2c,mq3,mq3c,wquark
        
        double precision xv, colorfactor 
        integer k,l
        double precision Krondelta
        
        double complex g0, g0L, g0R, g1, g1L, g1R, g2, g2L, g2R
        double precision Mloop0, Mloop1, Mloop2

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
    
    	NeuQ2qx_M2VertexV = 0d0
        vertSS = 0d0
        vertTT = 0d0
        vertUU = 0d0
        vertST = 0d0
        vertTS = 0d0
        vertSU = 0d0
        vertUS = 0d0
        vertUT = 0d0
        vertTU = 0d0

! ----- Correction to Neutralino - Sfermion - Fermion Coupling in S channel (A1 Coupling) -----

        A1sL = 0d0
        A1sR = 0d0
        A1vL = 0d0
        A1vR = 0d0

        mq1 = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0

! Set propagator to tree-level

        PsL = Mf(pstype,psgen)
        PsR = Mf(pstype,psgen)
        PvL = 1d0
        PvR = 1d0        

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

! --- Gluon loop ---

        g0 = - GS
        g2 = - GS
        g1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        g1R = AR(chimeigen,sfmeigen,sftype,sfgen)
        
        Mloop0 = Mf(pstype,psgen) 

    	call myCget(m1**2,m2**2,sqrtS**2,Mf(pstype,psgen)**2,m2**2,0d0)
    	call myBget(m2**2,0d0,m2**2)

        A1sL = A1sL + (g0*g2*(-(g1R*m1*Mloop0*(2*xc0 + xc1)) + g1L*(xb0 + Mloop0**2*xc0 + (-2*m2**2 + sman)*xc2 + 2*m1**2*(xc1 + xc2)) ))/(16.*Pi**2)
        A1sR = A1sR + (g0*g2*(-(g1L*m1*Mloop0*(2*xc0 + xc1)) + g1R*(xb0 + Mloop0**2*xc0 + (-2*m2**2 + sman)*xc2 + 2*m1**2*(xc1 + xc2))))/(16.*Pi**2)
        A1vL = A1vL + (g0*g2*(g1L*Mloop0*(xc0 - xc2) - g1R*m1*(xc1 + 2*xc2)))/(16.*Pi**2)
        A1vR = A1vR + (g0*g2*(g1R*Mloop0*(xc0 - xc2) - g1L*m1*(xc1 + 2*xc2)))/(16.*Pi**2)

! --- Gluino loop ---

        do k=1,2

        g0L = JR(sfmeigen,sftype,sfgen)
        g0R = JL(sfmeigen,sftype,sfgen)
        g1L = dconjg(AR(chimeigen,k,sftype,sfgen))
        g1R = dconjg(AL(chimeigen,k,sftype,sfgen))
        g2L = JR(k,sftype,sfgen)
        g2R = JL(k,sftype,sfgen)

        Mloop0 = MSf(k,sftype,sfgen)
        Mloop1 = Mf(pstype,psgen)
        Mloop2 = MGl

        call myCget(m1**2,m2**2,sqrtS**2,MSf(k,sftype,sfgen)**2,Mf(pstype,psgen)**2,MGl**2)
        call myBget(m2**2,Mf(pstype,psgen)**2,MGl**2)

        A1sL = A1sL - (g2L*(g0L*Mloop2*(g1L*Mloop1*xc0 + g1R*m1*(xc0 + xc1)) + 
     &                g0R*(g1R*m1*Mloop1*xc1 +g1L*(xb0 + Mloop0**2*xc0 + m1**2*xc1 + sman*xc2))))/(16.*Pi**2)
        A1sR = A1sR - (g2R*(g0R*Mloop2*(g1R*Mloop1*xc0 + g1L*m1*(xc0 + xc1)) + 
     &                g0L*(g1L*m1*Mloop1*xc1 + g1R*(xb0 + Mloop0**2*xc0 + m1**2*xc1 + sman*xc2))))/(16.*Pi**2)
        A1vL = A1vL - (g2R*(g0R*g1L*Mloop2*xc2 + g0L*(g1L*Mloop1*(xc0 + xc2) + g1R*m1*(xc0 + xc1 + xc2))))/(16.*Pi**2)
        A1vR = A1vR - (g2L*(g0L*g1R*Mloop2*xc2 + g0R*(g1R*Mloop1*(xc0 + xc2) + g1L*m1*(xc0 + xc1 + xc2))))/(16.*Pi**2)

        enddo
                
! --- Mss  ---            
                        
        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
 
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))
        
        mq1c = Mf(pstype,psgen)

        vertSS = colorfactor*CF*NC*2*MssV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A1Lc,A1Rc,BLc,BRc,mq1,mq1c,xv,wquark)

        NeuQ2qx_M2VertexV = NeuQ2qx_M2VertexV + vertSS

! --- Mst  ---
        
        do sumTc = 0, 1
                
        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))

        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
        vertST = colorfactor*CF*NC*2*MstV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A2Lc,A2Rc,CLc,mq1,mq2c,xv,wquark)
        
        NeuQ2qx_M2VertexV = NeuQ2qx_M2VertexV + vertST

        enddo

! --- Msu  ---

        do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

        vertSU = colorfactor*CF*NC*2*MsuV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A3Lc,A3Rc,DLc,DRc,mq1,mq3c,xv,wquark)

        NeuQ2qx_M2VertexV = NeuQ2qx_M2VertexV + vertSU
                                        
        enddo

! ----- Correction to Neutralino - Sfermion - Fermion Coupling in T channel (A2 Coupling) -----

        do sumT = 0, 1

        A2sL = 0d0
        A2sR = 0d0
        A2vL = 0d0
        A2vR = 0d0

        mq2 = MSf(ptmeigen+sumT,pttype,ptgen)

! Set C coupling to tree-level
        CLoop1 = CL(ptmeigen+sumT,sfmeigen,sftype,sfgen,xtype)        
        CLoop2 = 0d0

! --- Gluon loop ---

        g0 = - GS
        g2 = - GS
        g1L = AL(chimeigen,ptmeigen+sumT,pttype,ptgen)
        g1R = AR(chimeigen,ptmeigen+sumT,pttype,ptgen)
        
        Mloop0 = m3
       	Mloop1 = MSf(ptmeigen+sumT,pttype,ptgen)
  
   	    call myCget(m1**2,tman,m3**2,m3**2,MSf(ptmeigen+sumT,pttype,ptgen)**2,0d0)
   	    call myBget(tman,MSf(ptmeigen+sumT,pttype,ptgen)**2,0d0)
        
        A2sL = A2sL + (g0*g2*(g1L*(xb0 + Mloop0*(m3 + Mloop0)*xc0 + (m3**2 - m3*Mloop0 - 2*tman)*xc2 + 2*m1**2*(xc1 + xc2)) -  g1R*m1*(Mloop0*(2*xc0 + xc1) + m3*(xc1 + 2*xc2))))/(16.*Pi**2)
        A2sR = A2sR + (g0*g2*(g1R*(xb0 + Mloop0*(m3 + Mloop0)*xc0 + (m3**2 - m3*Mloop0 - 2*tman)*xc2 + 2*m1**2*(xc1 + xc2)) -  g1L*m1*(Mloop0*(2*xc0 + xc1) + m3*(xc1 + 2*xc2))))/(16.*Pi**2)
        A2vL = A2vL + 0d0
        A2vR = A2vR + 0d0

! --- Gluino loop ---

                do k=1,2

        g0L = JR(ptmeigen+SumT,pttype,ptgen)
        g0R = JL(ptmeigen+SumT,pttype,ptgen)
        g1L = dconjg(AR(chimeigen,k,pttype,ptgen))
        g1R = dconjg(AL(chimeigen,k,pttype,ptgen))
        g2L = JR(k,pttype,ptgen)
        g2R = JL(k,pttype,ptgen)
        
        Mloop0 = MSf(k,pttype,ptgen)
        Mloop1 = m3
        Mloop2 = MGl
     
        call myCget(m1**2,tman,m3**2,MSf(k,pttype,ptgen)**2,m3**2,MGl**2)
        call myBget(tman,m3**2,MGl**2)

        A2sL = A2sL - (g0L*(g2R*m3 + g2L*Mloop2)*(g1L*Mloop1*xc0 + g1R*m1*(xc0 + xc1)) + g0R*g2L*(g1R*m1*Mloop1*xc1 + g1L*(xb0 + Mloop0**2*xc0 + m1**2*xc1)) + 
     &                g0L*g2R*m3*(g1R*m1 + g1L*Mloop1)*xc2 + g0R*g1L*m3*(g2L*m3 + g2R*Mloop2)*xc2)/(16.*Pi**2)
        A2sR = A2sR - (g0R*(g2L*m3 + g2R*Mloop2)*(g1R*Mloop1*xc0 + g1L*m1*(xc0 + xc1)) + g0L*g2R*(g1L*m1*Mloop1*xc1 + g1R*(xb0 + Mloop0**2*xc0 + m1**2*xc1)) +
     &                g0R*g2L*m3*(g1L*m1 + g1R*Mloop1)*xc2 + g0L*g1R*m3*(g2R*m3 + g2L*Mloop2)*xc2)/(16.*Pi**2)
        A2vL = A2vL + 0d0
        A2vR = A2vR + 0d0

                enddo
        
! --- Mts  ---

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0

        vertTS = colorfactor*CF*NC*2*MtsV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A1Lc,A1Rc,BLc,BRc,mq2,mq1c,xv,wquark)

        NeuQ2qx_M2VertexV = NeuQ2qx_M2VertexV + vertTS
                
! --- Mtt  ---

                do sumTc = 0, 1

        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
    	mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
        vertTT = colorfactor*CF*NC*2*MttV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A2Lc,A2Rc,CLc,mq2,mq2c,xv)

        NeuQ2qx_M2VertexV = NeuQ2qx_M2VertexV + vertTT
                
                enddo

! --- Mtu  ---

                do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

        vertTU = colorfactor*CF*NC*2*MtuV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A3Lc,A3Rc,DLc,DRc,mq2,mq3c,xv)

        NeuQ2qx_M2VertexV = NeuQ2qx_M2VertexV + vertTU

                enddo

        enddo


! ----- Correction to Neutralino - Sfermion - Fermion Coupling in U channel (A3 Coupling) -----

        do sumU = 0,sumUend

        A3sL = 0d0
        A3sR = 0d0
        A3vL = 0d0
        A3vR = 0d0

        mq3 = MNeuCha(pumeigen+sumU)

! Set D coupling to tree-level
        
        DLLoop = DL(chimeigen,pumeigen + sumU,sftype,xtype)
        DRLoop = DR(chimeigen,pumeigen + sumU,sftype,xtype)

! --- Gluon loop ---

        g0 = - GS
        g2 = - GS
        g1L = AL(pumeigen + sumU,sfmeigen,sftype,sfgen)
        g1R = AR(pumeigen + sumU,sfmeigen,sftype,sfgen)

        Mloop0 = m3 

        call myCget(uman,m2**2,m3**2,m3**2,m2**2,0d0)
	    call myBget(m2**2,m2**2,0d0)

        A3sL = A3sL + (g0*g1L*g2*(xb0 + Mloop0*(m3 + Mloop0)*xc0 + 2*uman*xc1 - 2*m2**2*xc2 +  (m3**2 - m3*Mloop0 + 2*uman)*xc2))/(16.*Pi**2)
        A3sR = A3sR + (g0*g1R*g2*(xb0 + Mloop0*(m3 + Mloop0)*xc0 + 2*uman*xc1 - 2*m2**2*xc2 + (m3**2 - m3*Mloop0 + 2*uman)*xc2))/(16.*Pi**2)
        A3vL = A3vL - (g0*g1L*g2*(Mloop0*(2*xc0 + xc1) + m3*(xc1 + 2*xc2)))/(16.*Pi**2)
        A3vR = A3vR - (g0*g1R*g2*(Mloop0*(2*xc0 + xc1) + m3*(xc1 + 2*xc2)))/(16.*Pi**2)
       
! --- Gluino loop ---

                do k=1,2        

        g0L = JR(sfmeigen,sftype,sfgen)
        g0R = JL(sfmeigen,sftype,sfgen)
        g1L = dconjg(AR(pumeigen + sumU,k,ftype,fgen))
        g1R = dconjg(AL(pumeigen + sumU,k,ftype,fgen))
        g2L = JR(k,ftype,fgen)
        g2R = JL(k,ftype,fgen)

        Mloop0 = MSf(k,ftype,fgen)
        Mloop1 = Mf(sftype,sfgen)
        Mloop2 = MGl
                                
        call myCget(uman,m2**2,m3**2,MSf(k,ftype,fgen)**2,Mf(sftype,sfgen)**2,MGl**2)
        call myBget(m2**2,Mf(sftype,sfgen)**2,MGl**2)

        A3sL = A3sL - (g1L*(g0R*g2L*(xb0 + Mloop0**2*xc0 + uman*xc1) + g0R*m3*(g2L*m3 + g2R*Mloop2)*xc2 + g0L*Mloop1*(g2L*Mloop2*xc0 + g2R*m3*(xc0 + xc2))))/(16.*Pi**2)
        A3sR = A3sR - (g1R*(g0L*g2R*(xb0 + Mloop0**2*xc0 + uman*xc1) +  g0L*m3*(g2R*m3 + g2L*Mloop2)*xc2 + g0R*Mloop1*(g2R*Mloop2*xc0 + g2L*m3*(xc0 + xc2))))/(16.*Pi**2)
        A3vL = A3vL - (g1L*(g0L*g2R*Mloop1*xc1 + g0R*(g2R*Mloop2*(xc0 + xc1) + g2L*m3*(xc0 + xc1 + xc2))))/(16.*Pi**2)
        A3vR = A3vR - (g1R*(g0R*g2L*Mloop1*xc1 + g0L*(g2L*Mloop2*(xc0 + xc1) + g2R*m3*(xc0 + xc1 + xc2))))/(16.*Pi**2)

                enddo

! --- Mus  ---

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
 
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
         wquark = 0d0

        vertUS = colorfactor*CF*NC*2*MusV(A3sL,A3sR,A3vL,A3vR,DLLoop,DRLoop,A1Lc,A1Rc,BLc,BRc,mq3,mq1c,xv,wquark)

        NeuQ2qx_M2VertexV = NeuQ2qx_M2VertexV + vertUS

! --- Mut  ---
        
                do sumTc = 0, 1
                
        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
        vertUT = colorfactor*CF*NC*2*MutV(A3sL,A3sR,A3vL,A3vR,DLLoop,DRLoop,A2Lc,A2Rc,CLc,mq3,mq2c,xv)
        
        NeuQ2qx_M2VertexV = NeuQ2qx_M2VertexV + vertUT

                enddo

! --- Muu  ---

                do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

        vertUU = colorfactor*CF*NC*2*MuuV(A3sL,A3sR,A3vL,A3vR,DLLoop,DRLoop,A3Lc,A3Rc,DLc,DRc,mq3,mq3c,xv)

        NeuQ2qx_M2VertexV = NeuQ2qx_M2VertexV + vertUU
                                        
                enddo
        enddo

! ----- Correction to Fermion - Fermion - Vector Coupling in S channel (B Coupling) -----

        BsgL  = 0d0
        BsgR  = 0d0
        Bsm1L = 0d0
        Bsm1R = 0d0
        Bsm2L = 0d0
        Bsm2R = 0d0
        BvgL  = 0d0
        BvgR  = 0d0
        Bvm1L = 0d0
        Bvm1R = 0d0
        Bvm2L = 0d0
        Bvm2R = 0d0

        mq1 = Mf(pstype,psgen) 
c        wquark = width(pstype)
        wquark = 0d0

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
 
! --- Gluon loop ---

        g0L = BL(xtype,ftype)
        g0R = BR(xtype,ftype)
        g1L = BL(1,sftype)      
        g1R = BR(1,sftype)        
        g2L = BL(1,ftype)         
        g2R = BR(1,ftype)        

        Mloop0 = 0d0
        Mloop1 = Mf(pstype,psgen)
        Mloop2 = Mf(ftype,fgen)

        call myCget(sman,m4**2,m3**2,0d0,Mf(pstype,psgen)**2,Mf(ftype,fgen)**2)
        call myBget(m4**2,Mf(pstype,psgen)**2,Mf(ftype,fgen)**2)

        BsgL = BsgL + ( - (g0R*g1L*g2L*Mloop1*Mloop2*xc0) + g0R*g1R*g2R*m3*mq1*(xc0 + xc1 + xc2) + g0L*g1L*g2L*(m1**2 + m2**2 + 2*m3**2 - tman - uman)*(xc0 + xc1 + xc2) + 
     &              g0L*g1L*g2L*(xb0 + Mloop0**2*xc0 - 2*xc00 + sman*xc1 + m3**2*xc2))/(8.*Pi**2)
        BsgR = BsgR + ( - (g0L*g1R*g2R*Mloop1*Mloop2*xc0) + g0L*g1L*g2L*m3*mq1*(xc0 + xc1 + xc2) + g0R*g1R*g2R*(m1**2 + m2**2 + 2*m3**2 - tman - uman)*(xc0 + xc1 + xc2) + 
     &              g0R*g1R*g2R*(xb0 + Mloop0**2*xc0 - 2*xc00 + sman*xc1 + m3**2*xc2))/(8.*Pi**2)
        Bsm1L = Bsm1L - (g0R*g2R*(-(g1L*Mloop1*(xc0 + xc1 + xc2)) + g1R*mq1*(xc0 + 2*xc1 + xc11 + xc12 + xc2)) + 
     &             g0L*g1L*(-(g2R*Mloop2*(xc0 + xc1 + xc2)) + g2L*m3*(xc0 + xc1 + xc12 + 2*xc2 + xc22)))/(8.*Pi**2)
        Bsm1R = Bsm1R - (g0L*g2L*(-(g1R*Mloop1*(xc0 + xc1 + xc2)) + g1L*mq1*(xc0 + 2*xc1 + xc11 + xc12 + xc2)) + 
     &             g0R*g1R*(-(g2L*Mloop2*(xc0 + xc1 + xc2)) + g2R*m3*(xc0 + xc1 + xc12 + 2*xc2 + xc22)))/(8.*Pi**2)
        Bsm2L = Bsm2L + (g0R*g2R*(-(g1L*Mloop1*(xc0 - xc1 + xc2)) + g1R*mq1*(xc0 - xc11 + xc12 + xc2)) - 
     &           g0L*g1L*(-(g2R*Mloop2*(xc0 + xc1 - xc2)) + g2L*m3*(xc0 + xc1 + xc12 - xc22)))/(8.*Pi**2)
        Bsm2R = Bsm2R + (g0L*g2L*(-(g1R*Mloop1*(xc0 - xc1 + xc2)) + g1L*mq1*(xc0 - xc11 + xc12 + xc2)) - 
     &           g0R*g1R*(-(g2L*Mloop2*(xc0 + xc1 - xc2)) + g2R*m3*(xc0 + xc1 + xc12 - xc22)))/(8.*Pi**2)
        BvgL = BvgL + (g0R*g1R*g2R*m3*(xc0 + xc1 + xc2))/(8.*Pi**2)
        BvgR = BvgR + (g0L*g1L*g2L*m3*(xc0 + xc1 + xc2))/(8.*Pi**2)
        Bvm1L = Bvm1L - (g0R*g1R*g2R*(xc0 + 2*xc1 + xc11 + xc12 + xc2))/(8.*Pi**2)
        Bvm1R = Bvm1R - (g0L*g1L*g2L*(xc0 + 2*xc1 + xc11 + xc12 + xc2))/(8.*Pi**2)
        Bvm2L = Bvm2L + (g0R*g1R*g2R*(xc0 - xc11 + xc12 + xc2))/(8.*Pi**2)
        Bvm2R = Bvm2R + (g0L*g1L*g2L*(xc0 - xc11 + xc12 + xc2))/(8.*Pi**2)
            
! --- Gluino loop ---

        do k=1,2
                do l=1,2
                
        g0  = - CL(l,k,pstype,psgen,xtype)
        g1L = JL(k,pstype,psgen)
        g1R = JR(k,pstype,psgen)
        g2L = JR(l,ftype,fgen)
        g2R = JL(l,ftype,fgen)

        Mloop0 = MGl
        Mloop1 = MSf(k,pstype,psgen)
        Mloop2 = MSf(l,ftype,fgen)

        call myCget(sman,m4**2,m3**2,MGl**2,MSf(k,pstype,psgen)**2,MSf(l,ftype,fgen)**2)

        BsgL = BsgL - (g0*g1L*g2R*xc00)/(8.*Pi**2)
        BsgR = BsgR - (g0*g1R*g2L*xc00)/(8.*Pi**2)
        Bsm1L = Bsm1L - (g0*(g1R*g2L*mq1*(xc1 + xc11 + xc12) + g1L*(-(g2L*Mloop0*(xc0 + xc1 + xc2)) + g2R*m3*(xc12 + xc2 + xc22))))/(16.*Pi**2)
        Bsm1R = Bsm1R - (g0*(g1L*g2R*mq1*(xc1 + xc11 + xc12) + g1R*(-(g2R*Mloop0*(xc0 + xc1 + xc2)) + g2L*m3*(xc12 + xc2 + xc22))))/(16.*Pi**2)
        Bsm2L = Bsm2L + (g0*(g1R*g2L*mq1*(-xc11 + xc12) + g1L*(g2L*Mloop0*(xc1 - xc2) + g2R*m3*(-xc12 + xc22))))/(16.*Pi**2)
        Bsm2R = Bsm2R + (g0*(g1L*g2R*mq1*(-xc11 + xc12) + g1R*(g2R*Mloop0*(xc1 - xc2) + g2L*m3*(-xc12 + xc22))))/(16.*Pi**2)
        BvgL = BvgL + 0d0
        BvgR = BvgR + 0d0
        Bvm1L = Bvm1L - (g0*g1R*g2L*(xc1 + xc11 + xc12))/(16.*Pi**2)
        Bvm1R = Bvm1R - (g0*g1L*g2R*(xc1 + xc11 + xc12))/(16.*Pi**2)
        Bvm2L = Bvm2L + (g0*g1R*g2L*(-xc11 + xc12))/(16.*Pi**2)
        Bvm2R = Bvm2R + (g0*g1L*g2R*(-xc11 + xc12))/(16.*Pi**2)

                enddo
        enddo
  
! --- Mss  ---            

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))

        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)

        vertSS = colorfactor*CF*NC*2*MssV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A1Lc,A1Rc,BLc,BRc,mq1,mq1c,xv,wquark)

        NeuQ2qx_M2VertexV = NeuQ2qx_M2VertexV + vertSS

! --- Mst  ---

        do sumTc = 0, 1

        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))

        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))

        mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

        vertST = colorfactor*CF*NC*2*MstV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A2Lc,A2Rc,CLc,mq1,mq2c,xv,wquark)

        NeuQ2qx_M2VertexV = NeuQ2qx_M2VertexV + vertST

        enddo

! --- Msu  ---

        do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

        vertSU = colorfactor*CF*NC*2*MsuV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A3Lc,A3Rc,DLc,DRc,mq1,mq3c,xv,wquark)

       NeuQ2qx_M2VertexV = NeuQ2qx_M2VertexV + vertSU

        enddo

! ----- Correction to Sfermion - Sfermion - Vector Coupling in T channel (C Coupling) -----

        do sumT = 0,1

        CLoop1 = 0d0
        CLoop2 = 0d0

        mq2 = MSf(ptmeigen+sumT,pttype,ptgen)
c        wquark = width(pstype)
        wquark = 0d0

! Set A coupling to tree-level

        A2sL = AL(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2sR = AR(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2vL = 0d0
        A2vR = 0d0

! --- Gluon loop ---

        g0 = - CL(ptmeigen+sumT,sfmeigen,sftype,sfgen,xtype)
        g2 = - CL(ptmeigen+sumT,ptmeigen+sumT,pttype,ptgen,1)
        g1 = - CL(sfmeigen,sfmeigen,sftype,sfgen,1)
        
        Mloop0 = 0d0
        Mloop1 = MSf(sfmeigen,sftype,sfgen)
       	Mloop2 = MSf(ptmeigen+sumT,pttype,ptgen)

       	call myCget(m2**2,m4**2,tman,0d0,MSf(sfmeigen,sftype,sfgen)**2,MSf(ptmeigen+sumT,pttype,ptgen)**2)
        call myBget(m4**2,Mloop1**2,Mloop2**2)

        CLoop1 = CLoop1 - (g0*g1*g2*(Mloop0**2*xc0 + 2*tman*xc0 + 4*xc00 + Mloop0**2*xc1 + 3*tman*xc1 + tman*xc11 + 4*tman*xc12 + Mloop0**2*xc2 + 5*tman*xc2 + 
     &             3*tman*xc22 - m4**2*(2*xc0 + 3*xc1 + xc11 + 2*xc12 + 3*xc2 + xc22) + m2**2*(2*xc0 + 5*xc1 + 3*xc11 + 4*xc12 + 3*xc2 + xc22)))/(16.*Pi**2)
        CLoop2 = CLoop2 + (g0*g1*g2*(xb0 + 2*xb1 - Mloop0**2*xc1 - 2*tman*xc1 - tman*xc11 - 2*tman*xc12 - Mloop0**2*xc2 + 2*tman*xc2 + 
     &             m4**2*(2*xc1 + xc11 - 2*xc2 - xc22) + 3*tman*xc22 +  m2**2*(-2*xc1 - 3*xc11 + 2*(xc12 + xc2) + xc22)))/(16.*Pi**2)


! --- Gluino loop ---
        
        g0L = BR(xtype,ftype)
        g0R = BL(xtype,ftype)
        g1L = JR(sfmeigen,sftype,sfgen)
        g1R = JL(sfmeigen,sftype,sfgen)
        g2L = JL(ptmeigen+sumT,pttype,ptgen)
        g2R = JR(ptmeigen+sumT,pttype,ptgen)

        Mloop0 = MGl
        Mloop1 = Mf(sftype,sfgen)
        Mloop2 = Mf(pttype,ptgen)

        call myCget(m2**2,m4**2,tman,MGl**2,Mf(sftype,sfgen)**2,Mf(pttype,ptgen)**2)
        call myBget(m4**2,Mf(sftype,sfgen)**2,Mf(pttype,ptgen)**2)

        CLoop1 = CLoop1 + (g0R*(g1L*Mloop1*(g2R*Mloop2*(xc1 + xc2) + g2L*Mloop0*(xc0 + xc1 + xc2)) + g1R*(g2R*Mloop0*Mloop2*(xc0 + xc1 + xc2) + 
     &             g2L*(xb0 + m2**2*xc1 + tman*xc2 + Mloop0**2*(2*xc0 + xc1 + xc2)))) + g0L*(g1R*Mloop1*(g2L*Mloop2*(xc1 + xc2) + g2R*Mloop0*(xc0 + xc1 + xc2)) + 
     &             g1L*(g2L*Mloop0*Mloop2*(xc0 + xc1 + xc2) + g2R*(xb0 + m2**2*xc1 + tman*xc2 + Mloop0**2*(2*xc0 + xc1 + xc2)))))/(16.*Pi**2)
        
        CLoop2 = CLoop2 - (g0R*(g1L*Mloop1*(g2R*Mloop2*(-xc1 + xc2) +  g2L*Mloop0*(xc0 - xc1 + xc2)) + g1R*(-(g2R*Mloop0*Mloop2*(xc0 + xc1 - xc2)) + 
     &             g2L*(xb0 + 2*xb1 + m2**2*xc1 - tman*xc2 + Mloop0**2*(-xc1 + xc2)))) + g0L*(g1R*Mloop1*(g2L*Mloop2*(-xc1 + xc2) +  g2R*Mloop0*(xc0 - xc1 + xc2)) + 
     &             g1L*(-(g2L*Mloop0*Mloop2*(xc0 + xc1 - xc2)) + g2R*(xb0 + 2*xb1 + m2**2*xc1 - tman*xc2 + Mloop0**2*(-xc1 + xc2)))))/(16.*Pi**2)
        
! --- Gluon loop 3 (with quartic coupling) ---

        g0 = CLL(ptmeigen+sumT,sfmeigen,sftype,sfgen,xtype)
        g1 = CL(sfmeigen,sfmeigen,sftype,sfgen,1)
        
        Mloop0 = 0d0
        Mloop1 = MSf(sfmeigen,sftype,sfgen)

        call myBget(m2**2,0d0,MSf(sfmeigen,sftype,sfgen)**2)

        CLoop1 = CLoop1 - (g0*g1*(2*xb0 + xb1))/(32.*Pi**2)
        CLoop2 = CLoop2 - (g0*g1*(2*xb0 + xb1))/(32.*Pi**2)

! --- Gluon loop 4 (with quartic coupling) ---

        g0 = CLL(ptmeigen+sumT,sfmeigen,sftype,sfgen,xtype)
        g1 = CL(ptmeigen+sumT,ptmeigen+sumT,pttype,ptgen,1)
        
        Mloop0 = 0d0
        Mloop1 = MSf(ptmeigen+sumT,pttype,ptgen)

        call myBget(tman,0d0,MSf(ptmeigen+sumT,pttype,ptgen)**2)

        CLoop1 = CLoop1 - (g0*g1*(2*xb0 + xb1))/(32.*Pi**2)
        CLoop2 = CLoop2 + (g0*g1*(2*xb0 + xb1))/(32.*Pi**2) 

! --- 4-Squark loop ---

        do k=1,2
                do l=1,2

        g0 = CL(k,l,sftype,sfgen,xtype) 
c        g1 = - GS**2*(Krondelta(sfmeigen,l) - 2*USfC(sfmeigen,1,sftype,sfgen)*USf(l,1,sftype,sfgen))*(Krondelta(k,ptmeigen+sumT) -  
c     &          2*USfC(k,1,pttype,ptgen)*USf(ptmeigen+sumT,1,pttype,ptgen))
        g1 = S(sfmeigen,l,k,ptmeigen+sumT,sftype,ftype,sfgen)
        
c        Mloop0 = MSf(l,sftype,sfgen)
c        Mloop1 = MSf(k,sftype,sfgen)

        call myBget(m4**2,MSf(k,ftype,fgen)**2,MSf(l,sftype,sfgen)**2)

        CLoop1 = CLoop1 + 0d0
        CLoop2 = CLoop2 + (g0*g1*(xb0 + 2*xb1))/(16.*Pi**2)   

                enddo
        enddo


! --- Mts  ---

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)

        vertTS = colorfactor*CF*NC*2*MtsV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A1Lc,A1Rc,BLc,BRc,mq2,mq1c,xv,wquark)
 
        NeuQ2qx_M2VertexV = NeuQ2qx_M2VertexV + vertTS

! --- Mtt  ---

                do sumTc = 0, 1
                
        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
    	mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
        vertTT = colorfactor*CF*NC*2*MttV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A2Lc,A2Rc,CLc,mq2,mq2c,xv)
        
        NeuQ2qx_M2VertexV = NeuQ2qx_M2VertexV + vertTT
                
                enddo

! --- Mtu  ---

                do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

        vertTU = colorfactor*CF*NC*2*MtuV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A3Lc,A3Rc,DLc,DRc,mq2,mq3c,xv)

        NeuQ2qx_M2VertexV = NeuQ2qx_M2VertexV + vertTU

                enddo

        enddo

	end
