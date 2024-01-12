************************************************
	double precision function NeuQ2qx_M2VertexG()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

        double precision MssV, MstV, MtsV, MttV
	double precision vertAs, vertAt, vertBs, vertCt

        double complex PvL,PvR,PsL,PsR
        
        double complex A1sL,A1sR,A1vL,A1vR,A2sL,A2sR,A2vL,A2vR
        double complex BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R
        double complex CLoop1,CLoop2
        double complex A1Lc,A1Rc,A2Lc,A2Rc,BLc,BRc,CLc,CRc
        
        double precision mq1,mq1c,mq2,mq2c,wquark
        
        double precision xv
        integer k,l
        double precision Krondelta
        
        double complex g0, g0L, g0R, g1, g1L, g1R, g2, g2L, g2R
        double precision Mloop0, Mloop1, Mloop2

        integer sumT, sumTc


       	NeuQ2qx_M2VertexG = 0d0
        
        xv = 0d0

        sumT = 0
               
           
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

    
        vertAs = 0d0
        vertAt = 0d0
        vertBs = 0d0
        vertCt = 0d0
        

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

        g0 = -GS
        g2 = -GS
        g1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        g1R = AR(chimeigen,sfmeigen,sftype,sfgen)
        
        Mloop0 = Mf(pstype,psgen) 

    	call myCget(m1**2,m2**2,sqrtS**2,Mf(pstype,psgen)**2,m2**2,0d0)
    	call myBget(m2**2,0d0,m2**2)

        A1sL = CF*(g0*g2*(-(g1R*m1*Mloop0*(2*xc0 + xc1)) + g1L*(xb0 + Mloop0**2*xc0 + (-2*m2**2 + sman)*xc2 + 2*m1**2*(xc1 + xc2)) ))/(16.*Pi**2)
        A1sR = CF*(g0*g2*(-(g1L*m1*Mloop0*(2*xc0 + xc1)) + g1R*(xb0 + Mloop0**2*xc0 + (-2*m2**2 + sman)*xc2 + 2*m1**2*(xc1 + xc2))))/(16.*Pi**2)
        A1vL = CF*(g0*g2*(g1L*Mloop0*(xc0 - xc2) - g1R*m1*(xc1 + 2*xc2)))/(16.*Pi**2)
        A1vR = CF*(g0*g2*(g1R*Mloop0*(xc0 - xc2) - g1L*m1*(xc1 + 2*xc2)))/(16.*Pi**2)

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

        A1sL = A1sL - CF*(g2L*(g0L*Mloop2*(g1L*Mloop1*xc0 + g1R*m1*(xc0 + xc1)) + 
     &                g0R*(g1R*m1*Mloop1*xc1 +g1L*(xb0 + Mloop0**2*xc0 + m1**2*xc1 + sman*xc2))))/(16.*Pi**2)
        A1sR = A1sR - CF*(g2R*(g0R*Mloop2*(g1R*Mloop1*xc0 + g1L*m1*(xc0 + xc1)) + 
     &                g0L*(g1L*m1*Mloop1*xc1 + g1R*(xb0 + Mloop0**2*xc0 + m1**2*xc1 + sman*xc2))))/(16.*Pi**2)
        A1vL = A1vL - CF*(g2R*(g0R*g1L*Mloop2*xc2 + g0L*(g1L*Mloop1*(xc0 + xc2) + g1R*m1*(xc0 + xc1 + xc2))))/(16.*Pi**2)
        A1vR = A1vR - CF*(g2L*(g0L*g1R*Mloop2*xc2 + g0R*(g1R*Mloop1*(xc0 + xc2) + g1L*m1*(xc0 + xc1 + xc2))))/(16.*Pi**2)

        enddo
                
                
! --- tree-level - s-channel ---   
                        
        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
 
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))
        
        mq1c = Mf(pstype,psgen)

        vertAs = CF*NC*2*MssV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A1Lc,A1Rc,BLc,BRc,mq1,mq1c,xv,wquark)

! --- tree-level - t-channel ---     
                
        A2Lc = dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))

        CLc = dconjg(CL(ptmeigen,sfmeigen,sftype,sfgen,xtype))
 
        mq2c = MSf(ptmeigen,pttype,ptgen)
 
        vertAs = vertAs + CF*NC*2*MstV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A2Lc,A2Rc,CLc,mq1,mq2c,xv,wquark)
       
        
        NeuQ2qx_M2VertexG = NeuQ2qx_M2VertexG + vertAs
        

! ----- Correction to Neutralino - Sfermion - Fermion Coupling in T channel (A2 Coupling) -----

        A2sL = 0d0
        A2sR = 0d0
        A2vL = 0d0
        A2vR = 0d0

        mq2 = MSf(ptmeigen,pttype,ptgen)
        
        PsL = 1d0
        PsR = 1d0
        PvL = 1d0
        PvR = 1d0       

! Set C coupling to tree-level
        CLoop1 = CL(ptmeigen,sfmeigen,sftype,sfgen,xtype)        
        CLoop2 = 0d0

! --- Gluon loop ---

        g0 = -GS
        g2 = -GS
        g1L = AL(chimeigen,ptmeigen,pttype,ptgen)
        g1R = AR(chimeigen,ptmeigen,pttype,ptgen)
        
        Mloop0 = m3
       	Mloop1 = MSf(ptmeigen,pttype,ptgen)
  
   	    call myCget(m1**2,tman,m3**2,m3**2,MSf(ptmeigen,pttype,ptgen)**2,0d0)
   	    call myBget(tman,MSf(ptmeigen,pttype,ptgen)**2,0d0)
        
        A2sL = CF*(g0*g2*(g1L*(xb0 + Mloop0*(m3 + Mloop0)*xc0 + (m3**2 - m3*Mloop0 - 2*tman)*xc2 + 2*m1**2*(xc1 + xc2)) -  g1R*m1*(Mloop0*(2*xc0 + xc1) + m3*(xc1 + 2*xc2))))/(16.*Pi**2)
        A2sR = CF*(g0*g2*(g1R*(xb0 + Mloop0*(m3 + Mloop0)*xc0 + (m3**2 - m3*Mloop0 - 2*tman)*xc2 + 2*m1**2*(xc1 + xc2)) -  g1L*m1*(Mloop0*(2*xc0 + xc1) + m3*(xc1 + 2*xc2))))/(16.*Pi**2)
        A2vL = 0d0
        A2vR = 0d0

! --- Gluino loop ---

                do k=1,2

        g0L = JR(ptmeigen,pttype,ptgen)
        g0R = JL(ptmeigen,pttype,ptgen)
        g1L = dconjg(AR(chimeigen,k,pttype,ptgen))
        g1R = dconjg(AL(chimeigen,k,pttype,ptgen))
        g2L = JR(k,pttype,ptgen)
        g2R = JL(k,pttype,ptgen)
        
        Mloop0 = MSf(k,pttype,ptgen)
        Mloop1 = Mf(pttype,ptgen)
        Mloop2 = MGl
     
        call myCget(m1**2,tman,Mf(pttype,ptgen)**2,MSf(k,pttype,ptgen)**2,Mf(pttype,ptgen)**2,MGl**2)
        call myBget(tman,Mf(pttype,ptgen)**2,MGl**2)
        

        A2sL = A2sL - CF*(g0L*(g2R*m3 + g2L*Mloop2)*(g1L*Mloop1*xc0 + g1R*m1*(xc0 + xc1)) + g0R*g2L*(g1R*m1*Mloop1*xc1 + g1L*(xb0 + Mloop0**2*xc0 + m1**2*xc1)) + 
     &                g0L*g2R*m3*(g1R*m1 + g1L*Mloop1)*xc2 + g0R*g1L*m3*(g2L*m3 + g2R*Mloop2)*xc2)/(16.*Pi**2)
        A2sR = A2sR - CF*(g0R*(g2L*m3 + g2R*Mloop2)*(g1R*Mloop1*xc0 + g1L*m1*(xc0 + xc1)) + g0L*g2R*(g1L*m1*Mloop1*xc1 + g1R*(xb0 + Mloop0**2*xc0 + m1**2*xc1)) +
     &                g0R*g2L*m3*(g1L*m1 + g1R*Mloop1)*xc2 + g0L*g1R*m3*(g2R*m3 + g2L*Mloop2)*xc2)/(16.*Pi**2)
        A2vL = A2vL + 0d0
        A2vR = A2vR + 0d0

                enddo
        
! --- tree-level - s-channel ---  

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0

        vertAt = CF*NC*2*MtsV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A1Lc,A1Rc,BLc,BRc,mq2,mq1c,xv,wquark)
 
 
! --- tree-level - t-channel ---  


        A2Lc = dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen,sfmeigen,sftype,sfgen,xtype))
 
    	mq2c = MSf(ptmeigen,pttype,ptgen)
 
        vertAt = vertAt + CF*NC*2*MttV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A2Lc,A2Rc,CLc,mq2,mq2c,xv)
        
        
        NeuQ2qx_M2VertexG = NeuQ2qx_M2VertexG + vertAt


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
 
! --- Gluon-Quark-Quark-loop (I)---

        g0L = BL(1,ftype)
        g0R = BR(1,ftype)
        g1L = BL(1,ftype)      
        g1R = BR(1,ftype)        
        g2L = BL(1,ftype)         
        g2R = BR(1,ftype)        

        Mloop0 = 0d0
        Mloop1 = m3
        Mloop2 = m3


        call myCget(sman,0d0,m3**2,0d0,m3**2,m3**2)
        call myBget(0d0,m3**2,m3**2)

        BsgL = (CF-(CV/2d0))*( - (g0R*g1L*g2L*Mloop1*Mloop2*xc0) + g0R*g1R*g2R*m3*mq1*(xc0 + xc1 + xc2) + g0L*g1L*g2L*(m1**2 + m2**2 + 2*m3**2 - tman - uman)*(xc0 + xc1 + xc2) + 
     &              g0L*g1L*g2L*(xb0 + Mloop0**2*xc0 - 2*xc00 + sman*xc1 + m3**2*xc2))/(8.*Pi**2)
        BsgR = (CF-(CV/2d0))*( - (g0L*g1R*g2R*Mloop1*Mloop2*xc0) + g0L*g1L*g2L*m3*mq1*(xc0 + xc1 + xc2) + g0R*g1R*g2R*(m1**2 + m2**2 + 2*m3**2 - tman - uman)*(xc0 + xc1 + xc2) + 
     &              g0R*g1R*g2R*(xb0 + Mloop0**2*xc0 - 2*xc00 + sman*xc1 + m3**2*xc2))/(8.*Pi**2)
        Bsm1L = -(CF-(CV/2d0))*(g0R*g2R*(-(g1L*Mloop1*(xc0 + xc1 + xc2)) + g1R*mq1*(xc0 + 2*xc1 + xc11 + xc12 + xc2)) + 
     &             g0L*g1L*(-(g2R*Mloop2*(xc0 + xc1 + xc2)) + g2L*m3*(xc0 + xc1 + xc12 + 2*xc2 + xc22)))/(8.*Pi**2)
        Bsm1R = -(CF-(CV/2d0))*(g0L*g2L*(-(g1R*Mloop1*(xc0 + xc1 + xc2)) + g1L*mq1*(xc0 + 2*xc1 + xc11 + xc12 + xc2)) + 
     &             g0R*g1R*(-(g2L*Mloop2*(xc0 + xc1 + xc2)) + g2R*m3*(xc0 + xc1 + xc12 + 2*xc2 + xc22)))/(8.*Pi**2)
        Bsm2L = (CF-(CV/2d0))*(g0R*g2R*(-(g1L*Mloop1*(xc0 - xc1 + xc2)) + g1R*mq1*(xc0 - xc11 + xc12 + xc2)) - 
     &           g0L*g1L*(-(g2R*Mloop2*(xc0 + xc1 - xc2)) + g2L*m3*(xc0 + xc1 + xc12 - xc22)))/(8.*Pi**2)
        Bsm2R = (CF-(CV/2d0))*(g0L*g2L*(-(g1R*Mloop1*(xc0 - xc1 + xc2)) + g1L*mq1*(xc0 - xc11 + xc12 + xc2)) - 
     &           g0R*g1R*(-(g2L*Mloop2*(xc0 + xc1 - xc2)) + g2R*m3*(xc0 + xc1 + xc12 - xc22)))/(8.*Pi**2)
        BvgL = (CF-(CV/2d0))*(g0R*g1R*g2R*m3*(xc0 + xc1 + xc2))/(8.*Pi**2)
        BvgR = (CF-(CV/2d0))*(g0L*g1L*g2L*m3*(xc0 + xc1 + xc2))/(8.*Pi**2)
        Bvm1L = -(CF-(CV/2d0))*(g0R*g1R*g2R*(xc0 + 2*xc1 + xc11 + xc12 + xc2))/(8.*Pi**2)
        Bvm1R = -(CF-(CV/2d0))*(g0L*g1L*g2L*(xc0 + 2*xc1 + xc11 + xc12 + xc2))/(8.*Pi**2)
        Bvm2L = (CF-(CV/2d0))*(g0R*g1R*g2R*(xc0 - xc11 + xc12 + xc2))/(8.*Pi**2)
        Bvm2R = (CF-(CV/2d0))*(g0L*g1L*g2L*(xc0 - xc11 + xc12 + xc2))/(8.*Pi**2)
!             
            
            
! --- Gluino-Squark-Squark-loop (II)---

        do k=1,2
                
        g0  = GS
        g1L = JL(k,pstype,psgen)
        g1R = JR(k,pstype,psgen)
        g2L = JR(k,ftype,fgen)
        g2R = JL(k,ftype,fgen)

        Mloop0 = MGl
        Mloop1 = MSf(k,pstype,psgen)
        Mloop2 = MSf(k,ftype,fgen)

        call myCget(sman,0d0,m3**2,MGl**2,MSf(k,pstype,psgen)**2,MSf(k,ftype,fgen)**2)
        
        BsgL = BsgL - (CF-(CV/2d0))*(g0*g1L*g2R*xc00)/(8.*Pi**2)
        BsgR = BsgR - (CF-(CV/2d0))*(g0*g1R*g2L*xc00)/(8.*Pi**2)
        Bsm1L = Bsm1L - (CF-(CV/2d0))*(g0*(g1R*g2L*mq1*(xc1 + xc11 + xc12) + g1L*(-(g2L*Mloop0*(xc0 + xc1 + xc2)) + g2R*m3*(xc12 + xc2 + xc22))))/(16.*Pi**2)
        Bsm1R = Bsm1R - (CF-(CV/2d0))*(g0*(g1L*g2R*mq1*(xc1 + xc11 + xc12) + g1R*(-(g2R*Mloop0*(xc0 + xc1 + xc2)) + g2L*m3*(xc12 + xc2 + xc22))))/(16.*Pi**2)
        Bsm2L = Bsm2L + (CF-(CV/2d0))*(g0*(g1R*g2L*mq1*(-xc11 + xc12) + g1L*(g2L*Mloop0*(xc1 - xc2) + g2R*m3*(-xc12 + xc22))))/(16.*Pi**2)
        Bsm2R = Bsm2R + (CF-(CV/2d0))*(g0*(g1L*g2R*mq1*(-xc11 + xc12) + g1R*(g2R*Mloop0*(xc1 - xc2) + g2L*m3*(-xc12 + xc22))))/(16.*Pi**2)
        BvgL = BvgL + 0d0
        BvgR = BvgR + 0d0
        Bvm1L = Bvm1L - (CF-(CV/2d0))*(g0*g1R*g2L*(xc1 + xc11 + xc12))/(16.*Pi**2)
        Bvm1R = Bvm1R - (CF-(CV/2d0))*(g0*g1L*g2R*(xc1 + xc11 + xc12))/(16.*Pi**2)
        Bvm2L = Bvm2L + (CF-(CV/2d0))*(g0*g1R*g2L*(-xc11 + xc12))/(16.*Pi**2)
        Bvm2R = Bvm2R + (CF-(CV/2d0))*(g0*g1L*g2R*(-xc11 + xc12))/(16.*Pi**2)

        enddo
        

! --- Gluon-Gluon-Quark-loop (1)---
                
        g0  = -GS
        g1L = -GS
        g1R = -GS
        g2L = -GS
        g2R = -GS

        Mloop0 = m3
        Mloop1 = 0d0
        Mloop2 = 0d0

        call myCget(sman,0d0,m3**2,m3**2,0d0,0d0)

        BsgL = BsgL - (-CV/2d0)*(g0*(3d0*g1R*mq1*(g2L*Mloop0*xc0 + g2R*m3*(xc1 + xc2)) + 
     &       g1L*(3d0*g2R*Mloop0*m3*xc0 + g2L*(12d0*xc00 + sman*(3d0*xc1 + 2*(xc11 + xc12 + xc2)) + 
     &             m3**2*(2d0*xc1 + 2d0*xc12 + 3d0*xc2 + 2d0*xc22)))))/(16.*Pi**2)
        BsgR = BsgR -(-CV/2d0)*(g0*(3d0*g1L*mq1*(g2R*Mloop0*xc0 + g2L*m3*(xc1 + xc2)) + 
     &       g1R*(3d0*g2L*Mloop0*m3*xc0 + g2R*(12d0*xc00 + sman*(3d0*xc1 + 2d0*(xc11 + xc12 + xc2)) + 
     &             m3**2*(2d0*xc1 + 2d0*xc12 + 3d0*xc2 + 2d0*xc22)))))/(16.*Pi**2)
     
        Bsm1L = Bsm1L - (-CV/2d0)*(g0*(g1R*g2R*mq1*(2d0*xc11 + 2d0*xc12 - xc2) + 
     &       g1L*(3d0*g2R*Mloop0*(xc1 + xc2) + g2L*m3*(-xc1 + 2d0*(xc12 + xc22)))))/(16.*Pi**2)
        Bsm1R = Bsm1R - (-CV/2d0)*(g0*(g1L*g2L*mq1*(2d0*xc11 + 2d0*xc12 - xc2) + 
     &       g1R*(3d0*g2L*Mloop0*(xc1 + xc2) + g2R*m3*(-xc1 + 2d0*(xc12 + xc22)))))/(16.*Pi**2)
     
        Bsm2L = Bsm2L + (-CV/2d0)*(g0*(-(g1R*g2R*mq1*(2d0*xc1 + 2d0*xc11 - 2d0*xc12 + xc2)) + 
     &      g1L*(3d0*g2R*Mloop0*(-xc1 + xc2) + g2L*m3*(xc1 + 2d0*(-xc12 + xc2 + xc22)))))/(16.*Pi**2)
        Bsm2R = Bsm2R + (-CV/2d0)*(g0*(-(g1L*g2L*mq1*(2d0*xc1 + 2d0*xc11 - 2d0*xc12 + xc2)) + 
     &      g1R*(3d0*g2L*Mloop0*(-xc1 + xc2) + g2R*m3*(xc1 + 2d0*(-xc12 + xc2 + xc22)))))/(16.*Pi**2)
     
        BvgL = BvgL + (-CV/2d0)*(-3d0*g0*g1R*(g2L*Mloop0*xc0 + g2R*m3*(xc1 + xc2)))/(16.*Pi**2)
        BvgR = BvgR + (-CV/2d0)*(-3d0*g0*g1R*(g2L*Mloop0*xc0 + g2R*m3*(xc1 + xc2)))/(16.*Pi**2)
        
        Bvm1L = Bvm1L + (-CV/2d0)*(g0*g1R*g2R*(-2d0*xc11 - 2d0*xc12 + xc2))/(16.*Pi**2)
        Bvm1R = Bvm1R + (-CV/2d0)*(g0*g1L*g2L*(-2d0*xc11 - 2d0*xc12 + xc2))/(16.*Pi**2)
        
        Bvm2L = Bvm2L - (-CV/2d0)*(g0*g1R*g2R*(2d0*xc1 + 2d0*xc11 - 2d0*xc12 + xc2))/(16.*Pi**2)
        Bvm2R = Bvm2R - (-CV/2d0)*(g0*g1L*g2L*(2d0*xc1 + 2d0*xc11 - 2d0*xc12 + xc2))/(16.*Pi**2)
        

        
! --- Gluino-Gluino-Squark-loop (2)--- 

        do k=1,2
        
        g0  = GS
        g1L = JL(k,pstype,psgen)
        g1R = JR(k,pstype,psgen)
        g2L = JR(k,pstype,psgen)
        g2R = JL(k,pstype,psgen)

        Mloop0 = MSf(k,pstype,psgen)
        Mloop1 = MGl
        Mloop2 = MGl

        call myCget(sman, 0d0, m3**2, MSf(k,pstype,psgen)**2, MGl**2, MGl**2)
      
        BsgL = BsgL + (-CV/2d0)*(g0*(-(g1R*mq1*(g2R*(-(Mloop1*xc1) + Mloop2*(xc0 + xc1)) + g2L*m3*(xc0 + xc1 + xc2))) + 
     &      g1L*(-(g2L*m3*(-(Mloop2*xc2) + Mloop1*(xc0 + xc2))) + g2R*(-(Mloop1*Mloop2*xc0) + 2*xc00 + sman*(xc1 + xc11 + xc12) + m3**2*(xc12 + xc2 + xc22)))))/(16.*Pi**2)
        BsgR = BsgR + (-CV/2d0)*(g0*(-(g1L*mq1*(g2L*(-(Mloop1*xc1) + Mloop2*(xc0 + xc1)) + g2R*m3*(xc0 + xc1 + xc2))) + 
     &      g1R*(-(g2R*m3*(-(Mloop2*xc2) + Mloop1*(xc0 + xc2))) + g2L*(-(Mloop1*Mloop2*xc0) + 2*xc00 + sman*xc1 + sman*xc11 + m3**2*xc12 + sman*xc12 + m3**2*xc2 + 
     &            m3**2*xc22))))/(16.*Pi**2)
        
        Bsm1L = Bsm1L - (-CV/2d0)*(g0*(g1R*g2L*mq1*(xc1 + xc11 + xc12) + g1L*(g2L*(Mloop1*xc1 + Mloop2*xc2) + g2R*m3*(xc12 + xc2 + xc22))))/(16.*Pi**2)
        Bsm1R = Bsm1R - (-CV/2d0)*(g0*(g1L*g2R*mq1*(xc1 + xc11 + xc12) + g1R*(g2R*(Mloop1*xc1 + Mloop2*xc2) + g2L*m3*(xc12 + xc2 + xc22))))/(16.*Pi**2)
     
        Bsm2L = Bsm2L + (-CV/2d0)*(g0*(-(g1R*g2L*mq1*(xc1 + xc11 - xc12)) + g1L*(g2L*(-(Mloop1*xc1) + Mloop2*xc2) + g2R*m3*(-xc12 + xc2 + xc22))))/(16.*Pi**2)
        Bsm2R = Bsm2R + (-CV/2d0)*(g0*(-(g1L*g2R*mq1*(xc1 + xc11 - xc12)) + g1R*(g2R*(-(Mloop1*xc1) + Mloop2*xc2) + g2L*m3*(-xc12 + xc2 + xc22))))/(16.*Pi**2)
     
        BvgL = BvgL - (-CV/2d0)*(g0*g1R*(g2R*(-(Mloop1*xc1) + Mloop2*(xc0 + xc1)) + g2L*m3*(xc0 + xc1 + xc2)))/(16.*Pi**2)
        BvgR = BvgR - (-CV/2d0)*(g0*g1L*(g2L*(-(Mloop1*xc1) + Mloop2*(xc0 + xc1)) + g2R*m3*(xc0 + xc1 + xc2)))/(16.*Pi**2)
        
        Bvm1L = Bvm1L - (-CV/2d0)*(g0*g1R*g2L*(xc1 + xc11 + xc12))/(16.*Pi**2)
        Bvm1R = Bvm1R - (-CV/2d0)*(g0*g1L*g2R*(xc1 + xc11 + xc12))/(16.*Pi**2)
        
        Bvm2L = Bvm2L - (-CV/2d0)*(g0*g1R*g2L*(xc1 + xc11 - xc12))/(16.*Pi**2)
        Bvm2R = Bvm2R - (-CV/2d0)*(g0*g1L*g2R*(xc1 + xc11 - xc12))/(16.*Pi**2)
!         
        enddo
          
  
! --- tree-level - s-channel ---  
   

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))

        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)

        vertBs = CF*NC*2*MssV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A1Lc,A1Rc,BLc,BRc,mq1,mq1c,xv,wquark)
        

! --- tree-level - t-channel ---  

        A2Lc = dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))

        CLc = dconjg(CL(ptmeigen,sfmeigen,sftype,sfgen,xtype))

        mq2c = MSf(ptmeigen,pttype,ptgen)

        vertBs = vertBs + CF*NC*2*MstV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A2Lc,A2Rc,CLc,mq1,mq2c,xv,wquark)


        NeuQ2qx_M2VertexG = NeuQ2qx_M2VertexG + vertBs
        
        
! ----- Correction to Sfermion - Sfermion - Vector Coupling in T channel (C Coupling) ----


        do sumT=0,1

        CLoop1 = 0d0
        CLoop2 = 0d0
        

                
! --- Gluon-Squark-Squark-loop (III)---

        g0 = GS*Krondelta(sfmeigen,1+sumT)
        g2 = GS*Krondelta(sfmeigen,1+sumT)
        g1 = GS*Krondelta(sfmeigen,1+sumT)
        
        Mloop0 = 0d0
        Mloop1 = m2
       	Mloop2 = m2

       	call myCget(m2**2,0d0,tman,0d0,m2**2,m2**2)
        call myBget(0d0,m2**2,m2**2)

        CLoop1 = CLoop1 - (CF-(CV/2d0))*(g0*g1*g2*(Mloop0**2*xc0 + 2*tman*xc0 + 4*xc00 + Mloop0**2*xc1 + 3*tman*xc1 + tman*xc11 + 4*tman*xc12 + Mloop0**2*xc2 + 5*tman*xc2 + 
     &             3*tman*xc22 - m4**2*(2*xc0 + 3*xc1 + xc11 + 2*xc12 + 3*xc2 + xc22) + m2**2*(2*xc0 + 5*xc1 + 3*xc11 + 4*xc12 + 3*xc2 + xc22)))/(16.*Pi**2)
     
        CLoop2 = CLoop2 + (CF-(CV/2d0))*(g0*g1*g2*(xb0 + 2*xb1 - Mloop0**2*xc1 - 2*tman*xc1 - tman*xc11 - 2*tman*xc12 - Mloop0**2*xc2 + 2*tman*xc2 + 
     &             m4**2*(2*xc1 + xc11 - 2*xc2 - xc22) + 3*tman*xc22 +  m2**2*(-2*xc1 - 3*xc11 + 2*(xc12 + xc2) + xc22)))/(16.*Pi**2)



! --- Gluino-Quark-Quark-loop (IV)---
      
        
        g0L = -GS
        g0R = -GS
        g1L = JR(sfmeigen,sftype,sfgen)
        g1R = JL(sfmeigen,sftype,sfgen)
        g2L = JL(1+sumT,pttype,ptgen)
        g2R = JR(1+sumT,pttype,ptgen)

        Mloop0 = MGl
        Mloop1 = Mf(sftype,sfgen)
        Mloop2 = Mf(pttype,ptgen)

        call myCget(m2**2,0d0,tman,MGl**2,Mf(sftype,sfgen)**2,Mf(pttype,ptgen)**2)
        call myBget(0d0,Mf(sftype,sfgen)**2,Mf(pttype,ptgen)**2)

        CLoop1 = CLoop1 + (CF-(CV/2d0))*(g0R*(g1L*Mloop1*(g2R*Mloop2*(xc1 + xc2) + g2L*Mloop0*(xc0 + xc1 + xc2)) + g1R*(g2R*Mloop0*Mloop2*(xc0 + xc1 + xc2) + 
     &             g2L*(xb0 + m2**2*xc1 + tman*xc2 + Mloop0**2*(2*xc0 + xc1 + xc2)))) + g0L*(g1R*Mloop1*(g2L*Mloop2*(xc1 + xc2) + g2R*Mloop0*(xc0 + xc1 + xc2)) + 
     &             g1L*(g2L*Mloop0*Mloop2*(xc0 + xc1 + xc2) + g2R*(xb0 + m2**2*xc1 + tman*xc2 + Mloop0**2*(2*xc0 + xc1 + xc2)))))/(16.*Pi**2)
        
        CLoop2 = CLoop2 - (CF-(CV/2d0))*(g0R*(g1L*Mloop1*(g2R*Mloop2*(-xc1 + xc2) +  g2L*Mloop0*(xc0 - xc1 + xc2)) + g1R*(-(g2R*Mloop0*Mloop2*(xc0 + xc1 - xc2)) + 
     &             g2L*(xb0 + 2*xb1 + m2**2*xc1 - tman*xc2 + Mloop0**2*(-xc1 + xc2)))) + g0L*(g1R*Mloop1*(g2L*Mloop2*(-xc1 + xc2) +  g2R*Mloop0*(xc0 - xc1 + xc2)) + 
     &             g1L*(-(g2L*Mloop0*Mloop2*(xc0 + xc1 - xc2)) + g2R*(xb0 + 2*xb1 + m2**2*xc1 - tman*xc2 + Mloop0**2*(-xc1 + xc2)))))/(16.*Pi**2)
        
        
! --- Gluon-Gluon-Squark-loop (3) ---     
! check minus sign!!!!!!

        g0 = GS
        g1 = GS*Krondelta(sfmeigen,1+sumT)
        g2 = GS*Krondelta(sfmeigen,1+sumT)

        Mloop0 = m2
        Mloop1 = 0d0
        Mloop2 = 0d0


        call myCget(m2**2,0d0,tman,m2**2,0d0,0d0)
     
        CLoop1 =  CLoop1 + (-CV/2d0)*(g0*g1*g2*(12*xc00 - (m2**2 - tman)*(xc1 - xc11 - xc2 + xc22)))/(32d0*Pi**2)
        
        CLoop2 = CLoop2 + (-CV/2d0)*(g0*g1*g2*(m2**2*(-xc0 + 4*xc1 - 3*xc11 + 2*xc12 + xc22) + tman*(xc0 - xc11 - 2*xc12 - 4*xc2 + 3*xc22)))/(32d0*Pi**2)



! --- Gluino-Gluino-Quark-loop (4 ^= IV)---        

! ! check minus sign!!!!!!
        g0L = GS
        g0R = GS
        g1L = JR(sfmeigen,sftype,sfgen)
        g1R = JL(sfmeigen,sftype,sfgen)
        g2L = JL(1+sumT,pttype,ptgen)
        g2R = JR(1+sumT,pttype,ptgen)

        Mloop0 = Mf(sftype,sfgen)
        Mloop1 = MGl
        Mloop2 = MGl

        call myCget(m2**2,0d0,tman,Mf(sftype,sfgen)**2,MGl**2,MGl**2)
        call myBget(0d0,MGl**2,MGl**2)

        CLoop1 = CLoop1 + (-CV/2d0)*(g0R*(g1L*Mloop1*(g2R*Mloop2*(xc1 + xc2) + g2L*Mloop0*(xc0 + xc1 + xc2)) + g1R*(g2R*Mloop0*Mloop2*(xc0 + xc1 + xc2) + 
     &             g2L*(xb0 + m2**2*xc1 + tman*xc2 + Mloop0**2*(2*xc0 + xc1 + xc2)))) + g0L*(g1R*Mloop1*(g2L*Mloop2*(xc1 + xc2) + g2R*Mloop0*(xc0 + xc1 + xc2)) + 
     &             g1L*(g2L*Mloop0*Mloop2*(xc0 + xc1 + xc2) + g2R*(xb0 + m2**2*xc1 + tman*xc2 + Mloop0**2*(2*xc0 + xc1 + xc2)))))/(16.*Pi**2)
        
        CLoop2 = CLoop2 - (-CV/2d0)*(g0R*(g1L*Mloop1*(g2R*Mloop2*(-xc1 + xc2) +  g2L*Mloop0*(xc0 - xc1 + xc2)) + g1R*(-(g2R*Mloop0*Mloop2*(xc0 + xc1 - xc2)) + 
     &             g2L*(xb0 + 2*xb1 + m2**2*xc1 - tman*xc2 + Mloop0**2*(-xc1 + xc2)))) + g0L*(g1R*Mloop1*(g2L*Mloop2*(-xc1 + xc2) +  g2R*Mloop0*(xc0 - xc1 + xc2)) + 
     &             g1L*(-(g2L*Mloop0*Mloop2*(xc0 + xc1 - xc2)) + g2R*(xb0 + 2*xb1 + m2**2*xc1 - tman*xc2 + Mloop0**2*(-xc1 + xc2)))))/(16.*Pi**2)

!      
! --- Gluon-Squark-bubble - quartic (5) ---   

        g0 = GS**2*Krondelta(sfmeigen,1+sumT)
        g1 = - GS*Krondelta(sfmeigen,1+sumT)

        Mloop0 = 0d0
        Mloop1 = m2
   
        call myBget(m2**2,0d0,m2**2)

        CLoop1 = CLoop1 - 2d0*(CF-(CV/2d0))*(g0*g1*(2*xb0 + xb1))/(32.*Pi**2)
        CLoop1 = CLoop1 + (-CV/2d0)*(g0*g1*(2*xb0 + xb1))/(32.*Pi**2)
        
        CLoop2 = CLoop2 + 2d0*(CF-(CV/2d0))*(g0*g1*(2*xb0 + xb1))/(32.*Pi**2)
        CLoop2 = CLoop2 - (-CV/2d0)*(g0*g1*(2*xb0 + xb1))/(32.*Pi**2)

        

! --- Gluon-Squark-bubble - quartic (5bis) ---   

        g0 = GS**2*Krondelta(sfmeigen,1+sumT)
        g1 = - GS*Krondelta(sfmeigen,1+sumT)

        Mloop0 = 0d0
        Mloop1 = m2
   
        call myBget(tman,0d0,m2**2)

           CLoop1 = CLoop1 - 2d0*(CF-(CV/2d0))*(g0*g1*(2*xb0 + xb1))/(32.*Pi**2)
         CLoop1 = CLoop1 + (-CV/2d0)*(g0*g1*(2*xb0 + xb1))/(32.*Pi**2)
          
          CLoop2 = CLoop2 + 2d0*(CF-(CV/2d0))*(g0*g1*(2*xb0 + xb1))/(32.*Pi**2)
         CLoop2 = CLoop2 - (-CV/2d0)*(g0*g1*(2*xb0 + xb1))/(32.*Pi**2)
        

! ! --- Gluon-Gluon-bubble - quartic (6) --- 
! ! zero because of color factor!!!
! 
!         g0 = -GS
!         g1 = GS**2*Krondelta(sfmeigen,1+sumT)
! 
! 
!         Mloop0 = 0d0
!         Mloop1 = 0d0
! 
!         call myBget(0d0,0d0,0d0)
! 
!         CLoop1 = CLoop1 + 0d0
!         
!         CLoop2 = CLoop2 - (3d0*(xb0 + 2*xb1)*g0*g1)/(16.*Pi**2)
!         
       

! --- 4-Squark loop ---


	    do k=1,2
	      do l=1,2
	      
	  call myBget(0d0,MSf(k,pttype,ptgen)**2,MSf(l,sftype,sfgen)**2)

        g0 = -GS**2*Aij(1+sumT,k,sftype,sfgen)*Aij(k,sfmeigen,sftype,sfgen)
	  g1 = -GS

	  CLoop1 = CLoop1 + 0d0
          CLoop2 = CLoop2 + (-1d0/(2d0*NC))*(g0*g1*(xb0 + 2*xb1))/(16.*Pi**2)  

	      enddo
	    enddo



! input for general amplitude        
        
        mq2 = MSf(1+sumT,pttype,ptgen)
c        wquark = width(pstype)
        wquark = 0d0
        
       PvL = 1d0
       PvR = 1d0
       PsL = 1d0
       PsR = 1d0

! Set A coupling to tree-level

        A2sL = AL(chimeigen,1+sumT,pttype,ptgen)
        A2sR = AR(chimeigen,1+sumT,pttype,ptgen)
        A2vL = 0d0
        A2vR = 0d0

        
        
! --- tree-level - s-channel ---  


        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)

        vertCt = vertCt + CF*NC*2*MtsV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A1Lc,A1Rc,BLc,BRc,mq2,mq1c,xv,wquark)


! --- tree-level - t-channel ---  

                
        A2Lc = dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen,sfmeigen,sftype,sfgen,xtype))
 
    	mq2c = MSf(ptmeigen,pttype,ptgen)
 
        vertCt = vertCt + CF*NC*2*MttV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A2Lc,A2Rc,CLc,mq2,mq2c,xv)
                
        enddo

        
        NeuQ2qx_M2VertexG = NeuQ2qx_M2VertexG + vertCt
!         NeuQ2qx_M2VertexG = vertCt
        
        
	end
