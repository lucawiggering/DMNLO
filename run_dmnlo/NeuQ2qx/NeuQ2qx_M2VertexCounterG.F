! =======================================================================
	  double precision function NeuQ2qx_M2VertexCounterG()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "NeuQ2qx_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"

        double precision MssV, MstV, MtsV, MttV
    	double precision vertAs, vertAt, vertBs, vertCt

        double complex PvL,PvR,PsL,PsR
        double complex A1vL,A1vR,A1sL,A1sR,A2vL,A2vR,A2sL,A2sR
        double complex Bsm1L,Bsm1R,Bsm2L,Bsm2R,BsgL,BsgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,BvgL,BvgR
        double complex CLoop1,CLoop2

        double complex A1Lc,A1Rc,A2Lc,A2Rc,BLc,BRc,CLc
        double precision mq1,mq1c,mq2,mq2c,wquark
        
        double precision xv, colorfactor 
        
        integer sumT, sumTc, k

        
        xv = 0d0
        colorfactor = 4d0/3d0
           
        sumT = 0
        sumTc = 0
        
        vertAs = 0d0
        vertAt = 0d0
        vertBs = 0d0
        vertCt = 0d0

    	NeuQ2qx_M2VertexCounterG = 0d0



! ----- Counter-term for Neutralino - Sfermion - Fermion Coupling in S channel (A1 Coupling) -----


! Set A coupling to counter-term

        A1sL = dAL(chimeigen,sfmeigen,sftype,sfgen)
        A1sR = dAR(chimeigen,sfmeigen,sftype,sfgen)
        A1vL = 0d0
        A1vR = 0d0
        

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

        mq1 = Mf(pstype,psgen) 
c        wquark = width(pstype)
        wquark = 0D0

! --- Mss  ---
        
        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
 
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)
        
        vertAs = CF*NC*2*MssV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A1Lc,A1Rc,BLc,BRc,mq1,mq1c,xv,wquark)

        testvar2 = vertAs
! --- Mst  ---


        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
	mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
        vertAs = vertAs + CF*NC*2*MstV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A2Lc,A2Rc,CLc,mq1,mq2c,xv,wquark)

        NeuQ2qx_M2VertexCounterG = NeuQ2qx_M2VertexCounterG + vertAs
        

! ----- Counter-term for Neutralino - Sfermion - Fermion Coupling in T channel (A2 Coupling) -----

        
! Set A coupling to counter-term

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
       
        vertAt = CF*NC*2*MtsV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A1Lc,A1Rc,BLc,BRc,mq2,mq1c,xv,wquark)
        

! --- Mtt  ---

        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
	mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
        vertAt = vertAt + CF*NC*2*MttV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A2Lc,A2Rc,CLc,mq2,mq2c,xv)
        
        NeuQ2qx_M2VertexCounterG = NeuQ2qx_M2VertexCounterG + vertAt


!         print*,"------------B---------------"
        
      
! ----- Counter-term for Vector - Fermion - Fermion Coupling in S channel (B Coupling) -----

! Set B coupling to counter-term
        
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
        
        mq1 = Mf(pstype,psgen) 
c        wquark = width(pstype)
        wquark = 0d0

! --- Mss  ---
        
        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
 
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))
        

        mq1c = Mf(pstype,psgen)
        
        vertBs = CF*NC*2d0*MssV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A1Lc,A1Rc,BLc,BRc,mq1,mq1c,xv,wquark)
        


! --- Mst  ---

        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
	 mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
 
        vertBs = vertBs + CF*NC*2d0*MstV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A2Lc,A2Rc,CLc,mq1,mq2c,xv,wquark)
   
        NeuQ2qx_M2VertexCounterG = NeuQ2qx_M2VertexCounterG + vertBs


! ----- Counter-term for Vector - Sfermion - Sfermion Coupling in T channel (C Coupling) -----

!         print*,"------------C---------------"


        do k = 0, 1

        CLoop1 = 0d0
        CLoop2 = 0d0

! Set C coupling to counter-term

        CLoop1 = dCL(ptmeigen+k,sfmeigen,sftype,sfgen,xtype)      
        CLoop2 = 0d0
        
        
! Set A coupling to tree-level

        A2sL = AL(chimeigen,ptmeigen+k,pttype,ptgen)
        A2sR = AR(chimeigen,ptmeigen+k,pttype,ptgen)
        A2vL = 0d0
        A2vR = 0d0

        mq2 = MSf(ptmeigen+k,pttype,ptgen)
        
! --- Mts  ---

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
 
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))

        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0D0
        
        vertCt = vertCt + CF*NC*2d0*MtsV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A1Lc,A1Rc,BLc,BRc,mq2,mq1c,xv,wquark)



! --- Mtt  ---


        A2Lc = dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen,sfmeigen,sftype,sfgen,xtype))
 
	mq2c = MSf(ptmeigen,pttype,ptgen)
 
        vertCt = vertCt + CF*NC*2d0*MttV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A2Lc,A2Rc,CLc,mq2,mq2c,xv)

        enddo
        
        
        NeuQ2qx_M2VertexCounterG = NeuQ2qx_M2VertexCounterG + vertCt
        


	  end
