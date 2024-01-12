! =========================================================================
        double precision function NeuQ2qx_M2BoxG()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

        double precision MsBoxV, MtBoxV
        double precision boxGluonS, boxGluinoS, boxGluonT, boxGluinoT, boxGluonST, boxGluinoST, boxGluonT4
        double precision boxFVS, boxFVT
        
        double complex Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R
        double complex Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R

        double complex Bs0LK,Bs0RK,Bs1LK,Bs1RK,Bs2LK,Bs2RK,Bs3LK,Bs3RK
        double complex Bv0LK,Bv0RK,Bv1LK,Bv1RK,Bv2LK,Bv2RK,Bv3LK,Bv3RK

        double complex A1Lc,A1Rc,A2Lc,A2Rc,BLc,BRc,CLc,DLc
        
        double precision mq1,mq1c,mq2,mq2c,wquark
        
        double precision xv, colorfactor 
        double precision Krondelta

        double complex g0, g0L, g0R, g1, g1L, g1R, g2, g2L, g2R, g3, g3L, g3R
        double precision Mloop0, Mloop1, Mloop2, Mloop3
        
        double complex xc0e, xc1e, xc2e, xc0g, xc1g, xc2g, xc0h, xc1h, xc2h
        double complex xc0a, xc1a, xc2a, xc0c, xc1c, xc2c

        integer k
	
        
        xv = 0d0
          

        g0 = 0d0
        g0L = 0d0
        g0R = 0d0
        g1L = 0d0
        g1R = 0d0
        g2 = 0d0
        g2L = 0d0
        g2R = 0d0
        g3 = 0d0
        g3L = 0d0
        g3R = 0d0
        
        Mloop0 = 0d0
        Mloop1 = 0d0
        Mloop2 = 0d0
        Mloop3 = 0d0

        mq1c = 0d0
        mq2c = 0d0

        boxGluonS = 0d0
        boxGluinoS = 0d0
        boxGluonT = 0d0
        boxGluinoT = 0d0
        boxGluonST = 0d0
        boxGluinoST = 0d0
        boxGluonT4 = 0d0
        
        Bs0L = 0d0
        Bs0R = 0d0
        Bs1L = 0d0
        Bs1R = 0d0
        Bs2L = 0d0
        Bs2R = 0d0
        Bs3L = 0d0
        Bs3R = 0d0
        
        Bv0L = 0d0
        Bv0R = 0d0
        Bv1L = 0d0
        Bv1R = 0d0
        Bv2L = 0d0
        Bv2R = 0d0
        Bv3L = 0d0
        Bv3R = 0d0
        
        xc0a = 0d0
        xc1a = 0d0
        xc1a = 0d0
        
        xc0c = 0d0
        xc1c = 0d0
        xc2c = 0d0
        
        xc0e = 0d0
        xc1e = 0d0
        xc2e = 0d0
        
        xc0g = 0d0
        xc1g = 0d0
        xc2g = 0d0
        
        xc0h = 0d0
        xc1h = 0d0
        xc2h = 0d0

        
        NeuQ2qx_M2BoxG = 0d0
        

! ----- s channel box with gluon (I)-----

        g0 = -GS
        g1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        g1R = AR(chimeigen,sfmeigen,sftype,sfgen)
        g2L = -GS  
        g2R = -GS
        g3 = -GS
        
        
        Mloop1 = MSf(sfmeigen,sftype,sfgen)
        Mloop2 = Mf(sftype,sfgen)
        Mloop3 = Mf(ftype,fgen)

        call myDget(m2**2,m1**2,0d0,m3**2,sman,uman,0d0,Mloop1**2,Mloop2**2,Mloop3**2)
        call myCget(m1**2,0d0,uman,Mloop1**2,Mloop2**2,Mloop3**2)

        Bs0L = -1.*g0*g3*(g1R*m1*(g2R*Mloop2*(m3 - 1.*Mloop3)*xd2 + g2L*(xc1 + 4.*xd00 + sman*xd12 - 1.*tman*xd2 + sman*xd22 - 1.*tman*xd22 - 1.*tman*xd23 + 
     &          m1**2*(-1.*xd12 + xd2 + xd22 + xd23) + m3*Mloop3*xd3 - 1.*uman*(2.*xd0 + xd1 + xd12 + xd13 + 3.*xd2 + xd22 + xd23 + 2.*xd3) + 
     &          m2**2*(2.*xd0 + 3.*xd1 + xd11 + 2.*xd12 + xd13 + 3.*xd2 + xd22 + xd23 + 2.*xd3) + m3**2*(2.*xd0 + xd1 + xd12 + xd13 + 4.*xd2 + 2.*xd22 + 4.*xd23 + 3.*xd3 + xd33)))
     &     + g1L*(g2L*Mloop2*(4.*xd00 + sman*(xd12 + xd22) - 1.*tman*(xd2 + xd23) + m1**2*(-1.*xd12 + xd2 + xd23) + m3*Mloop3*xd3 - 
     &          1.*uman*(2.*xd0 + xd1 + xd13 + xd2 + xd23 + 2.*xd3) + m2**2*(2.*xd0 + 3.*xd1 + xd11 + xd12 + xd13 + 3.*xd2 + xd23 + 2.*xd3) + 
     &          m3**2*(2.*xd0 + xd1 + xd13 + 2.*xd2 + 2.*xd23 + 3.*xd3 + xd33)) - 
     &       1.*g2R*(Mloop3*(2.*xd00 + sman*(xd12 + xd2 + xd22) - 1.*tman*xd23 + m1**2*(-1.*xd12 + xd23) - 1.*uman*(xd13 + xd23 + 2.*xd3) + 
     &             m2**2*(2.*xd0 + 3.*xd1 + xd11 + xd12 + xd13 + 2.*xd2 + xd23 + 2.*xd3)) + 
     &          m3*(xc2 + m1**2*xd12 - 1.*sman*(xd12 + xd2 + xd22) - 1.*m2**2*(2.*xd0 + 3.*xd1 + xd11 + xd12 + 2.*(xd2 + xd3))) + m3**3*xd33 + 
     &          m3**2*Mloop3*(xd13 + 2.*(xd23 + xd3) + xd33))))
        Bs0R = -1.*g0*g3*(g1L*m1*(g2L*Mloop2*(m3 - 1.*Mloop3)*xd2 + g2R*(xc1 + 4.*xd00 + sman*xd12 - 1.*tman*xd2 + sman*xd22 - 1.*tman*xd22 - 1.*tman*xd23 + 
     &          m1**2*(-1.*xd12 + xd2 + xd22 + xd23) + m3*Mloop3*xd3 - 1.*uman*(2.*xd0 + xd1 + xd12 + xd13 + 3.*xd2 + xd22 + xd23 + 2.*xd3) + 
     &          m2**2*(2.*xd0 + 3.*xd1 + xd11 + 2.*xd12 + xd13 + 3.*xd2 + xd22 + xd23 + 2.*xd3) + m3**2*(2.*xd0 + xd1 + xd12 + xd13 + 4.*xd2 + 2.*xd22 + 4.*xd23 + 3.*xd3 + xd33)))
     &     + g1R*(g2R*Mloop2*(4.*xd00 + sman*(xd12 + xd22) - 1.*tman*(xd2 + xd23) + m1**2*(-1.*xd12 + xd2 + xd23) + m3*Mloop3*xd3 - 
     &          1.*uman*(2.*xd0 + xd1 + xd13 + xd2 + xd23 + 2.*xd3) + m2**2*(2.*xd0 + 3.*xd1 + xd11 + xd12 + xd13 + 3.*xd2 + xd23 + 2.*xd3) + 
     &          m3**2*(2.*xd0 + xd1 + xd13 + 2.*xd2 + 2.*xd23 + 3.*xd3 + xd33)) - 
     &       1.*g2L*(Mloop3*(2.*xd00 + sman*(xd12 + xd2 + xd22) - 1.*tman*xd23 + m1**2*(-1.*xd12 + xd23) - 1.*uman*(xd13 + xd23 + 2.*xd3) + 
     &             m2**2*(2.*xd0 + 3.*xd1 + xd11 + xd12 + xd13 + 2.*xd2 + xd23 + 2.*xd3)) + 
     &          m3*(xc2 + m1**2*xd12 - 1.*sman*(xd12 + xd2 + xd22) - 1.*m2**2*(2.*xd0 + 3.*xd1 + xd11 + xd12 + 2.*(xd2 + xd3))) + m3**3*xd33 + 
     &          m3**2*Mloop3*(xd13 + 2.*(xd23 + xd3) + xd33))))
     
        Bs1L = -2.*g0*g3*(-1.*g1R*g2L*m1*(m3 - 1.*Mloop3)*(xd2 + xd22) + g1L*(g2L*Mloop2*(-1.*m3 + Mloop3)*xd2 + g2R*((m3*(m3 + Mloop3) - 2.*uman)*xd23 + 2.*m2**2*(xd12 + xd2 + xd22 + xd23))))
        Bs1R = -2.*g0*g3*(-1.*g1L*g2R*m1*(m3 - 1.*Mloop3)*(xd2 + xd22) + g1R*(g2R*Mloop2*(-1.*m3 + Mloop3)*xd2 + g2L*((m3*(m3 + Mloop3) - 2.*uman)*xd23 + 2.*m2**2*(xd12 + xd2 + xd22 + xd23))))
        
        Bs2L = 2.*g0*g3*(g1R*m1*(2.*g2R*Mloop2*xd2 + g2L*(m3 - 1.*Mloop3)*(2.*xd0 + xd1 + xd12 + 3.*xd2 + xd22) + 2.*g2L*m3*xd3) + 
     &    g1L*(g2L*Mloop2*((m3 - 1.*Mloop3)*(2.*xd0 + xd1 + xd2) + 2.*m3*xd3) + 
     &       g2R*(4.*xd00 - 2.*(m2**2 - 1.*sman)*(xd12 + xd2 + xd22) - 2.*tman*xd23 + 2.*m1**2*(-1.*xd12 + xd23) - 1.*m3*Mloop3*(xd13 + xd23 + 2.*xd3) + 
     &          m3**2*(xd13 + 3.*xd23 + 2.*(xd3 + xd33)))))
        Bs2R = 2.*g0*g3*(g1L*m1*(2.*g2L*Mloop2*xd2 + g2R*(m3 - 1.*Mloop3)*(2.*xd0 + xd1 + xd12 + 3.*xd2 + xd22) + 2.*g2R*m3*xd3) + 
     &    g1R*(g2R*Mloop2*((m3 - 1.*Mloop3)*(2.*xd0 + xd1 + xd2) + 2.*m3*xd3) + 
     &       g2L*(4.*xd00 - 2.*(m2**2 - 1.*sman)*(xd12 + xd2 + xd22) - 2.*tman*xd23 + 2.*m1**2*(-1.*xd12 + xd23) - 1.*m3*Mloop3*(xd13 + xd23 + 2.*xd3) + 
     &          m3**2*(xd13 + 3.*xd23 + 2.*(xd3 + xd33)))))
     
        Bs3L = 2.*g0*g3*(g1R*g2L*m1*(m3 - 1.*Mloop3)*xd23 - 1.*g1L*g2R*(xc2 + 2.*xd00 + (m1**2 - 1.*tman)*xd23 + m3*Mloop3*xd33 + m3**2*(xd13 + 2.*(xd23 + xd3) + 3.*xd33) - 
     &       1.*uman*(xd13 + xd23 + 2.*(xd3 + xd33)) + m2**2*(3.*xd13 + 3.*xd23 + 2.*(xd3 + xd33))))
        Bs3R = 2.*g0*g3*(g1L*g2R*m1*(m3 - 1.*Mloop3)*xd23 - 1.*g1R*g2L*(xc2 + 2.*xd00 + (m1**2 - 1.*tman)*xd23 + m3*Mloop3*xd33 + m3**2*(xd13 + 2.*(xd23 + xd3) + 3.*xd33) - 
     &       1.*uman*(xd13 + xd23 + 2.*(xd3 + xd33)) + m2**2*(3.*xd13 + 3.*xd23 + 2.*(xd3 + xd33))))
     
        Bv0L = g0*g3*(-1.*g1R*m1*(2.*g2R*Mloop2*xd2 + g2L*(m3 - 1.*Mloop3)*(2.*xd0 + xd1 + 2.*xd2) + 2.*g2L*m3*xd3) + 
     &    g1L*(g2L*Mloop2*((-1.*m3 + Mloop3)*(2.*xd0 + xd1 + xd2) - 2.*m3*xd3) + 
     &       g2R*(xc0 + xc2 - 8.*xd00 - 3.*sman*xd12 + tman*xd12 - 2.*sman*xd2 + tman*xd2 - 3.*sman*xd22 + tman*xd22 + 3.*tman*xd23 - 1.*m1**2*(-2.*xd12 + xd2 + xd22 + 3.*xd23) + 
     &          m3*Mloop3*xd3 + uman*(2.*xd0 + 3.*xd1 + xd11 + 2.*xd12 + 3.*xd13 + 3.*xd2 + xd22 + 3.*xd23 + 2.*xd3) - 
     &          1.*m2**2*(2.*xd0 + 5.*xd1 + 4.*xd11 + 5.*xd12 + 3.*xd13 + 3.*xd2 + xd22 + 3.*xd23 + 2.*xd3) - 
     &          1.*m3**2*(2.*xd0 + 3.*xd1 + xd11 + 3.*xd12 + 5.*xd13 + 4.*xd2 + 2.*xd22 + 8.*xd23 + 5.*xd3 + 3.*xd33))))
        Bv0R = g0*g3*(-1.*g1L*m1*(2.*g2L*Mloop2*xd2 + g2R*(m3 - 1.*Mloop3)*(2.*xd0 + xd1 + 2.*xd2) + 2.*g2R*m3*xd3) + 
     &    g1R*(g2R*Mloop2*((-1.*m3 + Mloop3)*(2.*xd0 + xd1 + xd2) - 2.*m3*xd3) + 
     &       g2L*(xc0 + xc2 - 8.*xd00 - 3.*sman*xd12 + tman*xd12 - 2.*sman*xd2 + tman*xd2 - 3.*sman*xd22 + tman*xd22 + 3.*tman*xd23 - 1.*m1**2*(-2.*xd12 + xd2 + xd22 + 3.*xd23) + 
     &          m3*Mloop3*xd3 + uman*(2.*xd0 + 3.*xd1 + xd11 + 2.*xd12 + 3.*xd13 + 3.*xd2 + xd22 + 3.*xd23 + 2.*xd3) - 
     &          1.*m2**2*(2.*xd0 + 5.*xd1 + 4.*xd11 + 5.*xd12 + 3.*xd13 + 3.*xd2 + xd22 + 3.*xd23 + 2.*xd3) - 
     &          1.*m3**2*(2.*xd0 + 3.*xd1 + xd11 + 3.*xd12 + 5.*xd13 + 4.*xd2 + 2.*xd22 + 8.*xd23 + 5.*xd3 + 3.*xd33))))
     
        Bv1L = -2.*g0*g3*(2.*g1R*g2L*m1*(xd2 + xd22) + g1L*(2.*g2L*Mloop2*xd2 + g2R*(-1.*m3 + Mloop3)*(xd12 + xd2 + xd22) - 2.*g2R*m3*xd23))
        Bv1R = -2.*g0*g3*(2.*g1L*g2R*m1*(xd2 + xd22) + g1R*(2.*g2R*Mloop2*xd2 + g2L*(-1.*m3 + Mloop3)*(xd12 + xd2 + xd22) - 2.*g2L*m3*xd23))
        
        Bv2L = 2.*g0*g3*(-2.*g1R*g2L*m1*(xd12 + xd2 + xd22) + g1L*g2R*(-1.*Mloop3*(2.*xd0 + 3.*xd1 + xd11 + 2.*xd12 + 3.*xd2 + xd22) + 
     &       m3*(2.*xd0 + 3.*xd1 + xd11 + 2.*xd12 + 2.*xd13 + 3.*xd2 + xd22 + 2.*(xd23 + xd3))))
        Bv2R = 2.*g0*g3*(-2.*g1L*g2R*m1*(xd12 + xd2 + xd22) + g1R*g2L*(-1.*Mloop3*(2.*xd0 + 3.*xd1 + xd11 + 2.*xd12 + 3.*xd2 + xd22) + 
     &       m3*(2.*xd0 + 3.*xd1 + xd11 + 2.*xd12 + 2.*xd13 + 3.*xd2 + xd22 + 2.*(xd23 + xd3))))
     
        Bv3L = 2.*g0*g3*(-2.*g1R*g2L*m1*xd23 + g1L*g2R*((m3 - 1.*Mloop3)*(xd13 + xd23 + 2.*xd3) + 2.*m3*xd33))
        Bv3R = 2.*g0*g3*(-2.*g1L*g2R*m1*xd23 + g1R*g2L*((m3 - 1.*Mloop3)*(xd13 + xd23 + 2.*xd3) + 2.*m3*xd33))



! --- Interference with s channel tree-level --

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))
        
        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0
        
        boxGluonS = CF*NC*2*(CF-(CV/2d0))*MsBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A1Lc,A1Rc,BLc,BRc,mq1c,xv)
        
! --- Interference with t channel tree-level --

                
        A2Lc = dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen,sfmeigen,sftype,sfgen,xtype))
 
    	mq2c = MSf(ptmeigen,pttype,ptgen)
c	    wquark = width(pstype)
         wquark = 0d0

        boxGluonS = boxGluonS + CF*NC*2*(CF-(CV/2d0))*MtBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A2Lc,A2Rc,CLc,mq2c,xv)


         NeuQ2qx_M2BoxG = NeuQ2qx_M2BoxG + boxGluonS
        
                
                
                

! ----- s channel box with gluino (II) -----

        Bs0L = 0d0
        Bs0R = 0d0
        Bs1L = 0d0
        Bs1R = 0d0
        Bs2L = 0d0
        Bs2R = 0d0
        Bs3L = 0d0
        Bs3R = 0d0
        
        Bv0L = 0d0
        Bv0R = 0d0
        Bv1L = 0d0
        Bv1R = 0d0
        Bv2L = 0d0
        Bv2R = 0d0
        Bv3L = 0d0
        Bv3R = 0d0

        do k = 1, 2
            

        g0L =  sqrt2*GS*USfC(sfmeigen,2,sftype,sfgen)
        g0R = -sqrt2*GS*USfC(sfmeigen,1,sftype,sfgen)
        g1L = AR(chimeigen,k,sftype,sfgen)
        g1R = AL(chimeigen,k,sftype,sfgen)
        g2 =  -GS
        g3L =  sqrt2*GS*USfC(k,2,ftype,fgen)
        g3R = -sqrt2*GS*USfC(k,1,ftype,fgen)

        Mloop0 = MGl 
        Mloop1 = Mf(sftype,sfgen) 
        Mloop2 = MSf(k,sftype,sfgen)
        Mloop3 = MSf(k,ftype,fgen)
 
        call myDget(m2**2,m1**2,0d0,m3**2,sman,uman,Mloop0**2,Mloop1**2,Mloop2**2,Mloop3**2)
        call myCget(m1**2,0d0,uman,Mloop1**2,Mloop2**2,Mloop3**2)
        
        Bs0L = Bs0L - 2.*g1L*g2*g3R*(g0R*Mloop0 + g0L*Mloop1)*xd00
        Bs0R = Bs0R - 2.*g1R*g2*g3L*(g0L*Mloop0 + g0R*Mloop1)*xd00
        
        Bs1L = Bs1L + g2*(g0L*(-1.*g1R*g3L*m1*Mloop0*(xd2 + 2.*xd22) + g1L*Mloop1*(g3L*Mloop0*(xd0 + 2.*xd2) - 1.*g3R*m3*(2.*xd23 + xd3))) + 
     &    g0R*(-1.*g1R*g3L*m1*Mloop1*(xd2 + 2.*xd22) + g1L*(g3L*(xc0 + 2.*xc1 + m2**2*(xd1 + 2.*xd12) + Mloop0**2*(xd0 + 2.*xd2) - 1.*(m1**2 - 1.*sman)*(xd2 + 2.*xd22)) - 
     &          1.*g3R*m3*Mloop0*(2.*xd23 + xd3))))
        Bs1R = Bs1R + g2*(g0R*(-1.*g1L*g3R*m1*Mloop0*(xd2 + 2.*xd22) + g1R*Mloop1*(g3R*Mloop0*(xd0 + 2.*xd2) - 1.*g3L*m3*(2.*xd23 + xd3))) + 
     &    g0L*(-1.*g1L*g3R*m1*Mloop1*(xd2 + 2.*xd22) + g1R*(g3R*(xc0 + 2.*xc1 + m2**2*(xd1 + 2.*xd12) + Mloop0**2*(xd0 + 2.*xd2) - 1.*(m1**2 - 1.*sman)*(xd2 + 2.*xd22)) - 
     &          1.*g3L*m3*Mloop0*(2.*xd23 + xd3))))
     
        Bs2L = Bs2L + g2*(g0L*(-1.*g1R*g3L*m1*Mloop0*(2.*xd12 + xd2 + 2.*xd22) + g1L*Mloop1*(g3L*Mloop0*(xd0 + 2.*(xd1 + xd2)) - 1.*g3R*m3*(2.*(xd13 + xd23) + xd3))) - 
     &    1.*g0R*(g1R*g3L*m1*Mloop1*(2.*xd12 + xd2 + 2.*xd22) + g1L*(g3L*
     &           (xc0 + 2.*xc2 - 1.*m2**2*(xd1 + 2.*(xd11 + xd12)) - 1.*Mloop0**2*(xd0 + 2.*(xd1 + xd2)) + (m1**2 - 1.*sman)*(2.*xd12 + xd2 + 2.*xd22)) + 
     &          g3R*m3*Mloop0*(2.*(xd13 + xd23) + xd3))))
        Bs2R = Bs2R + g2*(g0R*(-1.*g1L*g3R*m1*Mloop0*(2.*xd12 + xd2 + 2.*xd22) + g1R*Mloop1*(g3R*Mloop0*(xd0 + 2.*(xd1 + xd2)) - 1.*g3L*m3*(2.*(xd13 + xd23) + xd3))) - 
     &    1.*g0L*(g1L*g3R*m1*Mloop1*(2.*xd12 + xd2 + 2.*xd22) + g1R*(g3R*
     &           (xc0 + 2.*xc2 - 1.*m2**2*(xd1 + 2.*(xd11 + xd12)) - 1.*Mloop0**2*(xd0 + 2.*(xd1 + xd2)) + (m1**2 - 1.*sman)*(2.*xd12 + xd2 + 2.*xd22)) + 
     &          g3L*m3*Mloop0*(2.*(xd13 + xd23) + xd3))))
     
        Bs3L = Bs3L + g2*(g0L*(-1.*g1R*g3L*m1*Mloop0*(xd2 + 2.*xd23) + g1L*Mloop1*(g3L*Mloop0*(xd0 + 2.*xd3) - 1.*g3R*m3*(xd3 + 2.*xd33))) + 
     &    g0R*(-1.*g1R*g3L*m1*Mloop1*(xd2 + 2.*xd23) + g1L*(g3L*(xc0 + 2.*xc2 + m2**2*(xd1 + 2.*xd13) - 1.*(m1**2 - 1.*sman)*(xd2 + 2.*xd23) + Mloop0**2*(xd0 + 2.*xd3)) - 
     &          1.*g3R*m3*Mloop0*(xd3 + 2.*xd33))))
        Bs3R = Bs3R + g2*(g0R*(-1.*g1L*g3R*m1*Mloop0*(xd2 + 2.*xd23) + g1R*Mloop1*(g3R*Mloop0*(xd0 + 2.*xd3) - 1.*g3L*m3*(xd3 + 2.*xd33))) + 
     &    g0L*(-1.*g1L*g3R*m1*Mloop1*(xd2 + 2.*xd23) + g1R*(g3R*(xc0 + 2.*xc2 + m2**2*(xd1 + 2.*xd13) - 1.*(m1**2 - 1.*sman)*(xd2 + 2.*xd23) + Mloop0**2*(xd0 + 2.*xd3)) - 
     &          1.*g3L*m3*Mloop0*(xd3 + 2.*xd33))))
     
        Bv0L = Bv0L + 2.*g0R*g1L*g2*g3L*xd00
        Bv0R = Bv0R + 2.*g0L*g1R*g2*g3R*xd00
        
        Bv1L = Bv1L - 1.*g2*(g0L*g3R*(g1R*m1*(xd2 + 2.*xd22) + g1L*Mloop1*(xd1 + 2.*xd12 + xd2 + 2.*xd22)) + g0R*g1L*(g3R*Mloop0*(xd0 + xd1 + 2.*xd12 + 3.*xd2 + 2.*xd22) - 1.*g3L*m3*(2.*xd23 + xd3)))
        Bv1R = Bv1R - 1.*g2*(g0R*g3L*(g1L*m1*(xd2 + 2.*xd22) + g1R*Mloop1*(xd1 + 2.*xd12 + xd2 + 2.*xd22)) + g0L*g1R*(g3L*Mloop0*(xd0 + xd1 + 2.*xd12 + 3.*xd2 + 2.*xd22) - 1.*g3R*m3*(2.*xd23 + xd3)))
        
        Bv2L = Bv2L - 1.*g2*(g0L*g3R*(g1R*m1*(2.*xd12 + xd2 + 2.*xd22) + g1L*Mloop1*(xd1 + 2.*xd11 + 4.*xd12 + xd2 + 2.*xd22)) + 
     &    g0R*g1L*(g3R*Mloop0*(xd0 + 3.*xd1 + 2.*xd11 + 4.*xd12 + 3.*xd2 + 2.*xd22) - 1.*g3L*m3*(2.*(xd13 + xd23) + xd3)))
        Bv2R = Bv2R - 1.*g2*(g0R*g3L*(g1L*m1*(2.*xd12 + xd2 + 2.*xd22) + g1R*Mloop1*(xd1 + 2.*xd11 + 4.*xd12 + xd2 + 2.*xd22)) + 
     &    g0L*g1R*(g3L*Mloop0*(xd0 + 3.*xd1 + 2.*xd11 + 4.*xd12 + 3.*xd2 + 2.*xd22) - 1.*g3R*m3*(2.*(xd13 + xd23) + xd3)))
     
        Bv3L = Bv3L - 1.*g2*(g0L*g3R*(g1R*m1*(xd2 + 2.*xd23) + g1L*Mloop1*(xd1 + 2.*xd13 + xd2 + 2.*xd23)) + 
     &    g0R*g1L*(g3R*Mloop0*(xd0 + xd1 + 2.*xd13 + xd2 + 2.*(xd23 + xd3)) - 1.*g3L*m3*(xd3 + 2.*xd33)))
        Bv3R = Bv3R - 1.*g2*(g0R*g3L*(g1L*m1*(xd2 + 2.*xd23) + g1R*Mloop1*(xd1 + 2.*xd13 + xd2 + 2.*xd23)) + 
     &    g0L*g1R*(g3L*Mloop0*(xd0 + xd1 + 2.*xd13 + xd2 + 2.*(xd23 + xd3)) - 1.*g3R*m3*(xd3 + 2.*xd33)))

     
        enddo
     
! --- Interference with s channel tree-level --

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))
        
        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0
        
        boxGluinoS = CF*NC*2*(CF-(CV/2d0))*MsBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A1Lc,A1Rc,BLc,BRc,mq1c,xv)
        
        
! --- Interference with t channel tree-level --

                
        A2Lc = dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen,sfmeigen,sftype,sfgen,xtype))
 
    	mq2c = MSf(ptmeigen,pttype,ptgen)
c	    wquark = width(pstype)
         wquark = 0d0

        boxGluinoS = boxGluinoS +  CF*NC*2*(CF-(CV/2d0))*MtBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A2Lc,A2Rc,CLc,mq2c,xv)
   
         NeuQ2qx_M2BoxG = NeuQ2qx_M2BoxG + boxGluinoS
        


! ----- t channel box with gluon (III)-----


        g0 = -GS
        g1 = -GS
        g2L = AL(chimeigen,ptmeigen,ftype,fgen)
        g2R = AR(chimeigen,ptmeigen,ftype,fgen)
        g3 = -GS
        
        Mloop0 = 0d0
        Mloop1 = Mf(ftype,fgen)
        Mloop2 = MSf(sfmeigen,sftype,sfgen)
        Mloop3 = MSf(sfmeigen,sftype,sfgen)
     
        call myDget(m3**2,m1**2,0d0,m2**2,tman,uman,0d0,Mloop1**2,Mloop2**2,Mloop3**2)
        call myCget(m1**2,0d0,uman,Mloop1**2,Mloop2**2,Mloop3**2)


        Bs0L = 2.*g0*g1*g2L*g3*(m3 - 1.*Mloop1)*xd00
        Bs0R = 2.*g0*g1*g2R*g3*(m3 - 1.*Mloop1)*xd00
        
        Bs1L = g0*g1*g3*(g2R*m1*(m3 - 1.*Mloop1)*(xd2 + 2.*xd22) + g2L*(xc0 + 2.*xc1 - 1.*(m1**2 - 1.*tman)*(xd2 + 2.*xd22) + m3*Mloop1*(xd1 + 2.*xd12 + xd2 + 2.*xd22) + 
     &       m3**2*(2.*xd0 + 3.*xd1 + 6.*xd12 + 6.*xd2 + 4.*xd22 + 2.*xd23 + xd3) - 1.*uman*(2.*(xd0 + xd1 + 2.*xd12 + 3.*xd2 + 2.*xd22 + xd23) + xd3) + 
     &       m2**2*(2.*(xd0 + xd1 + 2.*xd12 + 3.*xd2 + 2.*xd22 + 3.*xd23) + 3.*xd3)))
     
        Bs1R = g0*g1*g3*(g2L*m1*(m3 - 1.*Mloop1)*(xd2 + 2.*xd22) + g2R*(xc0 + 2.*xc1 - 1.*(m1**2 - 1.*tman)*(xd2 + 2.*xd22) + m3*Mloop1*(xd1 + 2.*xd12 + xd2 + 2.*xd22) + 
     &       m3**2*(2.*xd0 + 3.*xd1 + 6.*xd12 + 6.*xd2 + 4.*xd22 + 2.*xd23 + xd3) - 1.*uman*(2.*(xd0 + xd1 + 2.*xd12 + 3.*xd2 + 2.*xd22 + xd23) + xd3) + 
     &       m2**2*(2.*(xd0 + xd1 + 2.*xd12 + 3.*xd2 + 2.*xd22 + 3.*xd23) + 3.*xd3)))
     
        Bs2L = -1.*g0*g1*g3*(g2R*m1*(m3 - 1.*Mloop1)*(xd2 + 2.*xd23) + g2L*(xc0 + 2.*xc2 + 8.*xd00 - 1.*(m1**2 - 1.*tman)*(xd2 + 2.*xd23) + m3*Mloop1*(xd1 + 2.*xd13 + xd2 + 2.*xd23) - 
     &       1.*uman*(2.*xd0 + 2.*xd1 + 4.*xd13 + 2.*xd2 + 4.*xd23 + 5.*xd3 + 2.*xd33) + m3**2*(2.*xd0 + 3.*xd1 + 6.*xd13 + 2.*xd2 + 4.*xd23 + 5.*xd3 + 2.*xd33) + 
     &       m2**2*(2.*xd0 + 2.*xd1 + 4.*xd13 + 2.*xd2 + 4.*xd23 + 7.*xd3 + 6.*xd33)))
        Bs2R = -1.*g0*g1*g3*(g2L*m1*(m3 - 1.*Mloop1)*(xd2 + 2.*xd23) + g2R*(xc0 + 2.*xc2 + 8.*xd00 - 1.*(m1**2 - 1.*tman)*(xd2 + 2.*xd23) + m3*Mloop1*(xd1 + 2.*xd13 + xd2 + 2.*xd23) - 
     &       1.*uman*(2.*xd0 + 2.*xd1 + 4.*xd13 + 2.*xd2 + 4.*xd23 + 5.*xd3 + 2.*xd33) + m3**2*(2.*xd0 + 3.*xd1 + 6.*xd13 + 2.*xd2 + 4.*xd23 + 5.*xd3 + 2.*xd33) + 
     &       m2**2*(2.*xd0 + 2.*xd1 + 4.*xd13 + 2.*xd2 + 4.*xd23 + 7.*xd3 + 6.*xd33)))
     
        Bs3L = g0*g1*g3*(-1.*g2R*m1*(m3 - 1.*Mloop1)*(2.*xd12 + xd2 + 2.*xd22) + 
     &    g2L*(xc0 + 2.*xc2 - 4.*xd00 + (m1**2 - 1.*tman)*(2.*xd12 + xd2 + 2.*xd22) - 1.*m3*Mloop1*(xd1 + 2.*xd11 + 4.*xd12 + xd2 + 2.*xd22) + 
     &       uman*(2.*(xd0 + 3.*xd1 + 2.*xd11 + 4.*xd12 + xd13 + 3.*xd2 + 2.*xd22 + xd23) + xd3) - 
     &       1.*m3**2*(2.*xd0 + 7.*xd1 + 2.*(3.*xd11 + 5.*xd12 + xd13 + 3.*xd2 + 2.*xd22 + xd23) + xd3) - 
     &       1.*m2**2*(2.*xd0 + 6.*xd1 + 4.*xd11 + 8.*xd12 + 6.*xd13 + 6.*xd2 + 4.*xd22 + 6.*xd23 + 3.*xd3)))
        Bs3R = g0*g1*g3*(-1.*g2L*m1*(m3 - 1.*Mloop1)*(2.*xd12 + xd2 + 2.*xd22) + 
     &    g2R*(xc0 + 2.*xc2 - 4.*xd00 + (m1**2 - 1.*tman)*(2.*xd12 + xd2 + 2.*xd22) - 1.*m3*Mloop1*(xd1 + 2.*xd11 + 4.*xd12 + xd2 + 2.*xd22) + 
     &       uman*(2.*(xd0 + 3.*xd1 + 2.*xd11 + 4.*xd12 + xd13 + 3.*xd2 + 2.*xd22 + xd23) + xd3) - 
     &       1.*m3**2*(2.*xd0 + 7.*xd1 + 2.*(3.*xd11 + 5.*xd12 + xd13 + 3.*xd2 + 2.*xd22 + xd23) + xd3) - 
     &       1.*m2**2*(2.*xd0 + 6.*xd1 + 4.*xd11 + 8.*xd12 + 6.*xd13 + 6.*xd2 + 4.*xd22 + 6.*xd23 + 3.*xd3)))
     
        Bv0L = 4.*g0*g1*g2L*g3*xd00
        Bv0R = 4.*g0*g1*g2R*g3*xd00
        
        Bv1L = -1.*g0*g1*g3*(2.*g2R*m1*(xd2 + 2.*xd22) + g2L*(-1.*Mloop1*(2.*(xd0 + 2.*xd2 + xd23) + xd3) + m3*(2.*(xd0 + xd1 + 2.*xd12 + 3.*xd2 + 2.*xd22 + xd23) + xd3)))
        Bv1R = -1.*g0*g1*g3*(2.*g2L*m1*(xd2 + 2.*xd22) + g2R*(-1.*Mloop1*(2.*(xd0 + 2.*xd2 + xd23) + xd3) + m3*(2.*(xd0 + xd1 + 2.*xd12 + 3.*xd2 + 2.*xd22 + xd23) + xd3)))
        
        Bv2L = g0*g1*g3*(2.*g2R*m1*(xd2 + 2.*xd23) + g2L*(-1.*Mloop1*(2.*xd0 + 5.*xd3 + 2.*xd33) + m3*(2.*xd0 + 2.*xd1 + 4.*xd13 + 2.*xd2 + 4.*xd23 + 5.*xd3 + 2.*xd33)))
        Bv2R = g0*g1*g3*(2.*g2L*m1*(xd2 + 2.*xd23) + g2R*(-1.*Mloop1*(2.*xd0 + 5.*xd3 + 2.*xd33) + m3*(2.*xd0 + 2.*xd1 + 4.*xd13 + 2.*xd2 + 4.*xd23 + 5.*xd3 + 2.*xd33)))
        
        Bv3L = g0*g1*g3*(2.*g2R*m1*(2.*xd12 + xd2 + 2.*xd22) + g2L*(-2.*Mloop1*(xd0 + 2.*xd1 + xd13 + 2.*xd2 + xd23) - 1.*Mloop1*xd3 + 
     &       m3*(2.*(xd0 + 3.*xd1 + 2.*xd11 + 4.*xd12 + xd13 + 3.*xd2 + 2.*xd22 + xd23) + xd3)))
        Bv3R = g0*g1*g3*(2.*g2L*m1*(2.*xd12 + xd2 + 2.*xd22) + g2R*(-2.*Mloop1*(xd0 + 2.*xd1 + xd13 + 2.*xd2 + xd23) - 1.*Mloop1*xd3 + 
     &       m3*(2.*(xd0 + 3.*xd1 + 2.*xd11 + 4.*xd12 + xd13 + 3.*xd2 + 2.*xd22 + xd23) + xd3)))
      

! --- Interference with s channel tree-level --

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))
        
        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0
        
        boxGluonT = CF*NC*2*(CF-(CV/2d0))*MsBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A1Lc,A1Rc,BLc,BRc,mq1c,xv)



! --- Interference with t channel tree-level --

                
        A2Lc = dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen,sfmeigen,sftype,sfgen,xtype))
 
    	mq2c = MSf(ptmeigen,pttype,ptgen)
c	    wquark = width(pstype)
        wquark = 0d0
 
        boxGluonT  =  boxGluonT + CF*NC*2*(CF-(CV/2d0))*MtBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A2Lc,A2Rc,CLc,mq2c,xv)
     
         NeuQ2qx_M2BoxG = NeuQ2qx_M2BoxG + boxGluonT



! ----- t channel box with gluino  (IV)-----

        Bs0L = 0d0
        Bs0R = 0d0
        Bs1L = 0d0
        Bs1R = 0d0
        Bs2L = 0d0
        Bs2R = 0d0
        Bs3L = 0d0
        Bs3R = 0d0
        
        Bv0L = 0d0
        Bv0R = 0d0
        Bv1L = 0d0
        Bv1R = 0d0
        Bv2L = 0d0
        Bv2R = 0d0
        Bv3L = 0d0
        Bv3R = 0d0

        do k = 1, 2

        g0L =  JR(sfmeigen,sftype,sfgen)
        g0R =  JL(sfmeigen,sftype,sfgen)
        g1L =  JR(k,ftype,fgen)
        g1R =  JL(k,ftype,fgen)
        g2L = AR(chimeigen,k,ftype,fgen)
        g2R = AL(chimeigen,k,ftype,fgen)
        g3L = -GS
        g3R = -GS

        Mloop0 = MGl
        Mloop1 = MSf(k,ftype,fgen)
        Mloop2 = Mf(ftype,fgen)
        Mloop3 = Mf(sftype,sfgen)

        call myDget(m3**2,m1**2,0d0,m2**2,tman,uman,Mloop0**2,Mloop1**2,Mloop2**2,Mloop3**2)
        call myCget(m1**2,0d0,uman,Mloop1**2,Mloop2**2,Mloop3**2) 


        Bs0L = Bs0L + g0L*(g1L*m3*Mloop0*(-1.*g2R*g3L*m1*xd1 + g2L*(-1.*g3L*Mloop2*(xd1 + xd2) + g3R*Mloop3*(xd0 + xd1 + xd2))) + 
     &     g1R*(g2L*(-1.*g3R*Mloop3*(xc0 + Mloop0**2*xd0 - 2.*xd00 + tman*xd2 + m2**2*xd3 - 1.*uman*xd3 + m3**2*(xd1 + xd3)) + 
     &           g3L*Mloop2*(xc0 + Mloop0**2*xd0 + m1**2*xd2 - 1.*sman*xd2 + m2**2*(xd2 + xd3))) + 
     &        g2R*m1*(-1.*g3R*Mloop2*Mloop3*xd2 + g3L*(xc0 + xc1 - 1.*uman*xd12 + Mloop0**2*(xd0 + xd2) + m1**2*(xd2 + xd22) - 1.*(sman + uman)*(xd2 + xd22) + 
     &              m3**2*(xd12 + xd2 + xd22) + m2**2*(xd12 + 2.*(xd2 + xd22 + xd23) + xd3))))) + 
     &  g0R*(g1L*m3*(-1.*g2R*g3L*m1*Mloop3*xd1 + g2L*(-1.*g3L*Mloop2*Mloop3*(xd1 + xd2) + 
     &           g3R*(-1.*xc2 - 1.*sman*xd12 + Mloop0**2*(xd0 + xd1 + xd2) - 1.*sman*(xd2 + xd22) + m1**2*(xd12 + xd2 + xd22) + m3**2*(xd1 + xd11 + 2.*xd12 + xd2 + xd22) - 
     &              1.*uman*(xd1 + xd11 + 2.*xd12 + xd2 + xd22) + m2**2*(xd1 + xd11 + 3.*xd12 + 2.*(xd13 + xd2 + xd22 + xd23) + xd3)))) - 
     &     1.*g1R*Mloop0*(-1.*g2R*m1*(-1.*g3R*Mloop2*xd2 + g3L*Mloop3*(xd0 + xd2)) + 
     &        g2L*(-1.*g3L*Mloop2*Mloop3*xd0 + g3R*(xc0 + Mloop0**2*xd0 - 1.*uman*xd0 - 2.*xd00 - 1.*uman*xd1 + tman*xd2 - 1.*uman*xd2 - 1.*uman*xd3 + 
     &              m3**2*(xd0 + 2.*xd1 + xd2 + xd3) + m2**2*(xd0 + xd1 + xd2 + 2.*xd3)))))
        Bs0R = Bs0R + g0R*(g1R*m3*Mloop0*(-1.*g2L*g3R*m1*xd1 + g2R*(-1.*g3R*Mloop2*(xd1 + xd2) + g3L*Mloop3*(xd0 + xd1 + xd2))) + 
     &     g1L*(g2R*(-1.*g3L*Mloop3*(xc0 + Mloop0**2*xd0 - 2.*xd00 + tman*xd2 + m2**2*xd3 - 1.*uman*xd3 + m3**2*(xd1 + xd3)) + 
     &           g3R*Mloop2*(xc0 + Mloop0**2*xd0 + m1**2*xd2 - 1.*sman*xd2 + m2**2*(xd2 + xd3))) + 
     &        g2L*m1*(-1.*g3L*Mloop2*Mloop3*xd2 + g3R*(xc0 + xc1 - 1.*uman*xd12 + Mloop0**2*(xd0 + xd2) + m1**2*(xd2 + xd22) - 1.*(sman + uman)*(xd2 + xd22) + 
     &              m3**2*(xd12 + xd2 + xd22) + m2**2*(xd12 + 2.*(xd2 + xd22 + xd23) + xd3))))) + 
     &  g0L*(g1R*m3*(-1.*g2L*g3R*m1*Mloop3*xd1 + g2R*(-1.*g3R*Mloop2*Mloop3*(xd1 + xd2) + 
     &           g3L*(-1.*xc2 - 1.*sman*xd12 + Mloop0**2*(xd0 + xd1 + xd2) - 1.*sman*(xd2 + xd22) + m1**2*(xd12 + xd2 + xd22) + m3**2*(xd1 + xd11 + 2.*xd12 + xd2 + xd22) - 
     &              1.*uman*(xd1 + xd11 + 2.*xd12 + xd2 + xd22) + m2**2*(xd1 + xd11 + 3.*xd12 + 2.*(xd13 + xd2 + xd22 + xd23) + xd3)))) - 
     &     1.*g1L*Mloop0*(-1.*g2L*m1*(-1.*g3L*Mloop2*xd2 + g3R*Mloop3*(xd0 + xd2)) + 
     &        g2R*(-1.*g3R*Mloop2*Mloop3*xd0 + g3L*(xc0 + Mloop0**2*xd0 - 1.*uman*xd0 - 2.*xd00 - 1.*uman*xd1 + tman*xd2 - 1.*uman*xd2 - 1.*uman*xd3 + 
     &              m3**2*(xd0 + 2.*xd1 + xd2 + xd3) + m2**2*(xd0 + xd1 + xd2 + 2.*xd3)))))
     
        Bs1L = Bs1L - 2.*(g0L*(g1R*g2L*g3R*m3*Mloop3*(xd12 + xd2 + xd22) - 1.*g1L*g3L*Mloop0*(g2L*Mloop2*xd2 + g2R*m1*(xd2 + xd22))) + 
     &    g0R*(g1R*g2L*g3R*m3*Mloop0*(xd12 + xd2 + xd22) + g1L*(-1.*g2R*g3L*m1*Mloop3*(xd2 + xd22) + 
     &          g2L*(-1.*g3L*Mloop2*Mloop3*xd2 + g3R*(m2**2 + m3**2 - 1.*uman)*(xd12 + xd2 + xd22) + g3R*m2**2*xd23))))
        Bs1R = Bs1R - 2.*(g0R*(g1L*g2R*g3L*m3*Mloop3*(xd12 + xd2 + xd22) - 1.*g1R*g3R*Mloop0*(g2R*Mloop2*xd2 + g2L*m1*(xd2 + xd22))) + 
     &    g0L*(g1L*g2R*g3L*m3*Mloop0*(xd12 + xd2 + xd22) + g1R*(-1.*g2L*g3R*m1*Mloop3*(xd2 + xd22) + 
     &          g2R*(-1.*g3R*Mloop2*Mloop3*xd2 + g3L*(m2**2 + m3**2 - 1.*uman)*(xd12 + xd2 + xd22) + g3L*m2**2*xd23))))
     
        Bs2L = Bs2L - 2.*(g0L*(g1L*g3L*Mloop0*(g2L*Mloop2*(xd0 + xd3) + g2R*m1*(xd0 + xd2 + xd23 + xd3)) - 
     &       1.*g1R*m3*(g2R*g3L*m1*xd1 + g2L*(g3L*Mloop2*(xd1 + xd2) + g3R*Mloop3*(xd13 + xd23 + xd3)))) + 
     &    g0R*(-1.*g1R*g2L*g3R*m3*Mloop0*(xd0 + xd1 + xd13 + xd2 + xd23 + xd3) + 
     &       g1L*(g2R*m1*(g3R*Mloop2*xd2 + g3L*Mloop3*(xd23 + xd3)) + g2L*
     &           (g3L*Mloop2*Mloop3*xd3 + g3R*(xc0 + Mloop0**2*xd0 - 2.*xd00 + tman*xd2 + m3**2*(xd1 - 1.*xd13 - 1.*xd23) + uman*(xd13 + xd23) - 1.*m2**2*(xd13 + xd23 + xd33))))))
        Bs2R = Bs2R - 2.*(g0R*(g1R*g3R*Mloop0*(g2R*Mloop2*(xd0 + xd3) + g2L*m1*(xd0 + xd2 + xd23 + xd3)) - 
     &       1.*g1L*m3*(g2L*g3R*m1*xd1 + g2R*(g3R*Mloop2*(xd1 + xd2) + g3L*Mloop3*(xd13 + xd23 + xd3)))) + 
     &    g0L*(-1.*g1L*g2R*g3L*m3*Mloop0*(xd0 + xd1 + xd13 + xd2 + xd23 + xd3) + 
     &       g1R*(g2L*m1*(g3L*Mloop2*xd2 + g3R*Mloop3*(xd23 + xd3)) + g2R*
     &           (g3R*Mloop2*Mloop3*xd3 + g3L*(xc0 + Mloop0**2*xd0 - 2.*xd00 + tman*xd2 + m3**2*(xd1 - 1.*xd13 - 1.*xd23) + uman*(xd13 + xd23) - 1.*m2**2*(xd13 + xd23 + xd33))))))
     
        Bs3L = Bs3L - 2.*(g0L*(-1.*g1R*g2L*g3R*m3*Mloop3*(xd1 + xd11 + 2.*xd12 + xd2 + xd22) + g1L*Mloop0*(g2L*g3R*Mloop3*(xd0 + xd1 + xd2) + g2R*g3L*m1*(xd12 + xd2 + xd22))) + 
     &    g0R*(-1.*g1R*g2L*g3R*m3*Mloop0*(xd1 + xd11 + 2.*xd12 + xd2 + xd22) + 
     &       g1L*(g2R*g3L*m1*Mloop3*(xd12 + xd2 + xd22) + g2L*g3R*(-1.*xc2 - 1.*sman*xd12 + Mloop0**2*(xd0 + xd1 + xd2) - 1.*sman*(xd2 + xd22) + m1**2*(xd12 + xd2 + xd22) + 
     &             m2**2*(xd12 + xd13 + xd2 + xd22 + xd23 + xd3)))))
        Bs3R = Bs3R - 2.*(g0R*(-1.*g1L*g2R*g3L*m3*Mloop3*(xd1 + xd11 + 2.*xd12 + xd2 + xd22) + g1R*Mloop0*(g2R*g3L*Mloop3*(xd0 + xd1 + xd2) + g2L*g3R*m1*(xd12 + xd2 + xd22))) + 
     &    g0L*(-1.*g1L*g2R*g3L*m3*Mloop0*(xd1 + xd11 + 2.*xd12 + xd2 + xd22) + 
     &       g1R*(g2L*g3R*m1*Mloop3*(xd12 + xd2 + xd22) + g2R*g3L*(-1.*xc2 - 1.*sman*xd12 + Mloop0**2*(xd0 + xd1 + xd2) - 1.*sman*(xd2 + xd22) + m1**2*(xd12 + xd2 + xd22) + 
     &             m2**2*(xd12 + xd13 + xd2 + xd22 + xd23 + xd3)))))
     
        Bv0L = Bv0L + g0L*(-1.*g1R*m3*(g2R*g3L*m1*xd1 + g2L*g3L*Mloop2*(xd1 + xd2) + g2L*g3R*Mloop3*xd3) + 
     &     g1L*Mloop0*(g2R*g3L*m1*(xd0 + xd2 + xd3) + g2L*(-1.*g3R*Mloop3*xd3 + g3L*Mloop2*(xd0 + xd3)))) + 
     &  g0R*(-1.*g1R*g2L*g3R*m3*Mloop0*(xd0 + xd1 + xd2 + xd3) + g1L*(g2R*m1*(g3R*Mloop2*xd2 + g3L*Mloop3*xd3) + 
     &        g2L*(g3L*Mloop2*Mloop3*xd3 + g3R*(xc0 - 1.*xc2 - 4.*xd00 + uman*xd13 + tman*xd2 - 1.*m1**2*xd23 + (sman + uman)*xd23 - 1.*m3**2*(-1.*xd1 + xd13 + xd23) + 
     &              Mloop0**2*(xd0 - 1.*xd3) - 1.*m2**2*(xd13 + 2.*(xd23 + xd33))))))
        Bv0R = Bv0R + g0R*(-1.*g1L*m3*(g2L*g3R*m1*xd1 + g2R*g3R*Mloop2*(xd1 + xd2) + g2R*g3L*Mloop3*xd3) + 
     &     g1R*Mloop0*(g2L*g3R*m1*(xd0 + xd2 + xd3) + g2R*(-1.*g3L*Mloop3*xd3 + g3R*Mloop2*(xd0 + xd3)))) + 
     &  g0L*(-1.*g1L*g2R*g3L*m3*Mloop0*(xd0 + xd1 + xd2 + xd3) + g1R*(g2L*m1*(g3L*Mloop2*xd2 + g3R*Mloop3*xd3) + 
     &        g2R*(g3R*Mloop2*Mloop3*xd3 + g3L*(xc0 - 1.*xc2 - 4.*xd00 + uman*xd13 + tman*xd2 - 1.*m1**2*xd23 + (sman + uman)*xd23 - 1.*m3**2*(-1.*xd1 + xd13 + xd23) + 
     &              Mloop0**2*(xd0 - 1.*xd3) - 1.*m2**2*(xd13 + 2.*(xd23 + xd33))))))
     
        Bv1L = Bv1L + 2.*(g0R*g2L*g3R*(g1L*m3*(xd12 + xd2 + xd22) - 1.*g1R*Mloop0*xd23) + g0L*g1R*(g2R*g3L*m1*(xd2 + xd22) + g2L*(g3L*Mloop2*xd2 - 1.*g3R*Mloop3*xd23)))
        Bv1R = Bv1R + 2.*(g0L*g2R*g3L*(g1R*m3*(xd12 + xd2 + xd22) - 1.*g1L*Mloop0*xd23) + g0R*g1L*(g2L*g3R*m1*(xd2 + xd22) + g2R*(g3R*Mloop2*xd2 - 1.*g3L*Mloop3*xd23)))
        
        Bv2L = Bv2L + 2.*(g0L*g1R*(-1.*g2R*g3L*m1*xd23 + g2L*g3R*Mloop3*xd33) + g0R*g2L*g3R*(-1.*g1L*m3*(xd13 + xd23) + g1R*Mloop0*(xd3 + xd33)))
        Bv2R = Bv2R + 2.*(g0R*g1L*(-1.*g2L*g3R*m1*xd23 + g2R*g3L*Mloop3*xd33) + g0L*g2R*g3L*(-1.*g1R*m3*(xd13 + xd23) + g1L*Mloop0*(xd3 + xd33)))
        
        Bv3L = Bv3L + 2.*(g0L*g1R*(-1.*g2R*g3L*m1*(xd12 + xd2 + xd22) + g2L*g3R*Mloop3*(xd13 + xd23 + xd3)) + 
     &    g0R*g2L*g3R*(-1.*g1L*m3*(xd1 + xd11 + 2.*xd12 + xd2 + xd22) + g1R*Mloop0*(xd0 + xd1 + xd13 + xd2 + xd23 + xd3)))
        Bv3R = Bv3R + 2.*(g0R*g1L*(-1.*g2L*g3R*m1*(xd12 + xd2 + xd22) + g2R*g3L*Mloop3*(xd13 + xd23 + xd3)) + 
     &    g0L*g2R*g3L*(-1.*g1R*m3*(xd1 + xd11 + 2.*xd12 + xd2 + xd22) + g1L*Mloop0*(xd0 + xd1 + xd13 + xd2 + xd23 + xd3)))


        enddo

! --- Interference with s channel tree-level --

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))
        
        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0
        
        boxGluinoT = CF*NC*2*(CF-(CV/2d0))*MsBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A1Lc,A1Rc,BLc,BRc,mq1c,xv)



! --- Interference with t channel tree-level --
                
        A2Lc = dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen,sfmeigen,sftype,sfgen,xtype))
 
        mq2c = MSf(ptmeigen,pttype,ptgen)
c        wquark = width(pstype)
        wquark = 0d0
 
        boxGluinoT = boxGluinoT + CF*NC*2*(CF-(CV/2d0))*MtBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A2Lc,A2Rc,CLc,mq2c,xv)
     
         NeuQ2qx_M2BoxG = NeuQ2qx_M2BoxG + boxGluinoT

        
        

! ----- Four-Vertex box with gluon (V)-----


        g0 = GS**2
        g1L = AL(chimeigen,ptmeigen,ftype,fgen)
        g1R = AR(chimeigen,ptmeigen,ftype,fgen)
        g2 = -GS

        Mloop0 = 0d0
        Mloop1 = MSf(ptmeigen,ftype,fgen)
        Mloop2 = Mf(ftype,fgen)

        call myCget(tman,m1**2,m3**2,0d0,Mloop1**2,Mloop2**2)
 
        Bs0L = -1.*g0*g2*(g1R*m1*xc1 + g1L*(-1.*Mloop2*xc0 + m3*(xc0 + xc1 + xc2)))
        Bs0R = -1.*g0*g2*(g1L*m1*xc1 + g1R*(-1.*Mloop2*xc0 + m3*(xc0 + xc1 + xc2)))
        Bs1L = 0.
        Bs1R = 0.
        Bs2L = 0.
        Bs2R = 0.
        Bs3L = 2.*g0*g1L*g2*(xc0 + xc1 + xc2)
        Bs3R = 2.*g0*g1R*g2*(xc0 + xc1 + xc2)
        Bv0L = 0.
        Bv0R = 0.
        Bv1L = 0.
        Bv1R = 0.
        Bv2L = 0.
        Bv2R = 0.
        Bv3L = 0.
        Bv3R = 0.

! --- Interference with s channel tree-level --

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))
        
        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0
        
        boxGluonT4 = CF*NC*2*(2d0*CF-(CV/2d0))*MsBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A1Lc,A1Rc,BLc,BRc,mq1c,xv)


! --- Interference with t channel tree-level --

               
        A2Lc = dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen,sfmeigen,sftype,sfgen,xtype))
 
        mq2c = MSf(ptmeigen,pttype,ptgen)
c        wquark = width(pstype)
        wquark = 0d0
 
        boxGluonT4 = boxGluonT4 +  CF*NC*2*(2d0*CF-(CV/2d0))*MtBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A2Lc,A2Rc,CLc,mq2c,xv)
     
     
         NeuQ2qx_M2BoxG = NeuQ2qx_M2BoxG + boxGluonT4
        
	
! ----- ST box with gluon (1)-----

        Bs0L = 0d0
        Bs0R = 0d0
        Bs1L = 0d0
        Bs1R = 0d0
        Bs2L = 0d0
        Bs2R = 0d0
        Bs3L = 0d0
        Bs3R = 0d0
        
        Bv0L = 0d0
        Bv0R = 0d0
        Bv1L = 0d0
        Bv1R = 0d0
        Bv2L = 0d0
        Bv2R = 0d0
        Bv3L = 0d0
        Bv3R = 0d0


        Bs0LK = 0d0
        Bs0RK = 0d0
        Bs1LK = 0d0
        Bs1RK = 0d0
        Bs2LK = 0d0
        Bs2RK = 0d0
        Bs3LK = 0d0
        Bs3RK = 0d0
        
        Bv0LK = 0d0
        Bv0RK = 0d0
        Bv1LK = 0d0
        Bv1RK = 0d0
        Bv2LK = 0d0
        Bv2RK = 0d0
        Bv3LK = 0d0
        Bv3RK = 0d0
        


        g0 = -GS
        g1L = AL(chimeigen,ptmeigen,ftype,fgen)
        g1R = AR(chimeigen,ptmeigen,ftype,fgen)
        g2 = -GS
        g3 = -GS

        Mloop0 = 0d0
        Mloop1 = Mf(ftype,fgen)
        Mloop2 = MSf(sfmeigen,sftype,sfgen)
        Mloop3 = 0d0

        
        call myDget(m3**2, m1**2, m2**2, 0d0, tman, sman, 0d0, m3**2, m2**2, 0d0)
        call myCget(m1**2,m2**2,sman,Mf(ftype,fgen)**2,MSf(sfmeigen,sftype,sfgen)**2,0d0)
        

c        Bs0L = (g0*g2*g3*(-(g1R*m1*(xc1 + xc2 - tman*(xd2 + xd3) + m2**2*(xd12 + xd13 + xd2 + 2*(xd22 + xd23) + xd3) + (m1**2 - sman)*(xd22 + 2*xd23 + xd33) + m3**2*(xd12 + xd13 + xd22 + 2*xd23 + xd33) - 
c     &           uman*(xd12 + xd13 + xd22 + 2*xd23 + xd33))) - g1L*m3*(-xc0 - tman*(xd1 + xd2 + xd3) + m2**2*(xd1 + xd11 + 3*xd12 + xd13 + xd2 + 2*(xd22 + xd23) + xd3) + (m1**2 - sman)*(xd12 + xd13 + xd22 + 2*xd23 + xd33) + 
c     &         m3**2*(xd11 + 2*xd12 + 2*xd13 + xd22 + 2*xd23 + xd33) - uman*(xd11 + 2*xd12 + 2*xd13 + xd22 + 2*xd23 + xd33))))
c        Bs0R = (g0*g2*g3*(-(g1L*m1*(xc1 + xc2 - tman*(xd2 + xd3) + m2**2*(xd12 + xd13 + xd2 + 2*(xd22 + xd23) + xd3) + (m1**2 - sman)*(xd22 + 2*xd23 + xd33) + m3**2*(xd12 + xd13 + xd22 + 2*xd23 + xd33) - 
c     &           uman*(xd12 + xd13 + xd22 + 2*xd23 + xd33))) - g1R*m3*(-xc0 - tman*(xd1 + xd2 + xd3) + m2**2*(xd1 + xd11 + 3*xd12 + xd13 + xd2 + 2*(xd22 + xd23) + xd3) + (m1**2 - sman)*(xd12 + xd13 + xd22 + 2*xd23 + xd33) + 
c     &         m3**2*(xd11 + 2*xd12 + 2*xd13 + xd22 + 2*xd23 + xd33) - uman*(xd11 + 2*xd12 + 2*xd13 + xd22 + 2*xd23 + xd33))))
c        
c        Bs1L =  (g0*g2*g3*(-(g1R*m1*m3*xd1) + g1L*(xc1 + xc2 + 2*m3**2*xd0 - tman*xd0 - 3*uman*xd0 + 3*m3**2*xd1 - tman*xd1 - 3*uman*xd1 + 6*m3**2*xd12 - 4*uman*xd12 + 6*m3**2*xd13 - 4*uman*xd13 + 7*m3**2*xd2 - tman*xd2 - 7*uman*xd2 + 
c     &         5*m3**2*xd22 + tman*xd22 - 4*uman*xd22 + 9*m3**2*xd23 + 2*tman*xd23 - 7*uman*xd23 + 7*m3**2*xd3 - sman*xd3 - tman*xd3 - 7*uman*xd3 + m1**2*(xd0 + xd1 - xd22 - 2*xd23 + xd3 - xd33) + 
c     &         m2**2*(3*xd0 + 3*xd1 + 4*xd12 + 4*xd13 + 7*xd2 + 4*xd22 + 3*xd23 + 5*xd3 - xd33) + (4*m3**2 + tman - 3*uman)*xd33)))
c        Bs1R = (g0*g2*g3*(-(g1L*m1*m3*xd1) + g1R*(xc1 + xc2 + 2*m3**2*xd0 - tman*xd0 - 3*uman*xd0 + 3*m3**2*xd1 - tman*xd1 - 3*uman*xd1 + 6*m3**2*xd12 - 4*uman*xd12 + 6*m3**2*xd13 - 4*uman*xd13 + 7*m3**2*xd2 - tman*xd2 - 7*uman*xd2 + 
c     &         5*m3**2*xd22 + tman*xd22 - 4*uman*xd22 + 9*m3**2*xd23 + 2*tman*xd23 - 7*uman*xd23 + 7*m3**2*xd3 - sman*xd3 - tman*xd3 - 7*uman*xd3 + m1**2*(xd0 + xd1 - xd22 - 2*xd23 + xd3 - xd33) + 
c     &         m2**2*(3*xd0 + 3*xd1 + 4*xd12 + 4*xd13 + 7*xd2 + 4*xd22 + 3*xd23 + 5*xd3 - xd33) + (4*m3**2 + tman - 3*uman)*xd33)))
c        
c        Bs2L =   (g0*g2*g3*(3*g1R*m1*m3*xd1 + g1L*(2*xc0 + xc2 + 2*m3**2*xd0 + 3*tman*xd0 + uman*xd0 + 3*m3**2*xd1 + 3*tman*xd1 + uman*xd1 + 6*m3**2*xd13 - 4*uman*xd13 + m3**2*xd2 + 5*tman*xd2 + uman*xd2 + 5*m3**2*xd23 + tman*xd23 - 
c     &         4*uman*xd23 + 3*m3**2*xd3 + 3*sman*xd3 + 5*tman*xd3 - uman*xd3 + (4*m3**2 + tman - 3*uman)*xd33 - m2**2*(xd0 + xd1 - 4*xd13 + xd2 - 4*xd23 + xd3 + xd33) - m1**2*(3*xd0 + 3*xd1 + 2*xd2 + xd23 + 5*xd3 + xd33))))
c
c        Bs2R = (g0*g2*g3*(3*g1L*m1*m3*xd1 + g1R*(2*xc0 + xc2 + 2*m3**2*xd0 + 3*tman*xd0 + uman*xd0 + 3*m3**2*xd1 + 3*tman*xd1 + uman*xd1 + 6*m3**2*xd13 - 4*uman*xd13 + m3**2*xd2 + 5*tman*xd2 + uman*xd2 + 5*m3**2*xd23 + tman*xd23 - 
c     &         4*uman*xd23 + 3*m3**2*xd3 + 3*sman*xd3 + 5*tman*xd3 - uman*xd3 + (4*m3**2 + tman - 3*uman)*xd33 - m2**2*(xd0 + xd1 - 4*xd13 + xd2 - 4*xd23 + xd3 + xd33) - m1**2*(3*xd0 + 3*xd1 + 2*xd2 + xd23 + 5*xd3 + xd33))))
c
c        
c        Bs3L = (g0*g2*g3*(g1R*m1*m3*xd1 + g1L*(xc0 - 2*m3**2*xd0 - tman*xd0 + 3*uman*xd0 - 2*xd00 - 5*m3**2*xd1 - tman*xd1 + 5*uman*xd1 - 4*m3**2*xd11 + 2*uman*xd11 - 7*m3**2*xd12 - 2*sman*xd12 - tman*xd12 + 4*uman*xd12 - 
c     &         6*m3**2*xd13 - 2*sman*xd13 - tman*xd13 + 3*uman*xd13 - 5*m3**2*xd2 - 2*sman*xd2 - tman*xd2 + 5*uman*xd2 - 3*m3**2*xd22 - 2*sman*xd22 - tman*xd22 + 2*uman*xd22 - 5*m3**2*xd23 - 4*sman*xd23 - 2*tman*xd23 + 
c     &         3*uman*xd23 - 5*m3**2*xd3 - sman*xd3 - tman*xd3 + 5*uman*xd3 - m2**2*(xd0 + 3*xd1 + 2*xd11 + 2*xd12 + xd13 + xd2 - xd23 + 3*xd3 - xd33) - (2*m3**2 + 2*sman + tman - uman)*xd33 + 
c     &         m1**2*(-xd0 - xd1 + 3*xd12 + 3*xd13 + 2*xd2 + 3*xd22 + 6*xd23 + xd3 + 3*xd33))))
c        Bs3R = (g0*g2*g3*(g1L*m1*m3*xd1 + g1R*(xc0 - 2*m3**2*xd0 - tman*xd0 + 3*uman*xd0 - 2*xd00 - 5*m3**2*xd1 - tman*xd1 + 5*uman*xd1 - 4*m3**2*xd11 + 2*uman*xd11 - 7*m3**2*xd12 - 2*sman*xd12 - tman*xd12 + 4*uman*xd12 - 
c     &         6*m3**2*xd13 - 2*sman*xd13 - tman*xd13 + 3*uman*xd13 - 5*m3**2*xd2 - 2*sman*xd2 - tman*xd2 + 5*uman*xd2 - 3*m3**2*xd22 - 2*sman*xd22 - tman*xd22 + 2*uman*xd22 - 5*m3**2*xd23 - 4*sman*xd23 - 2*tman*xd23 + 
c     &         3*uman*xd23 - 5*m3**2*xd3 - sman*xd3 - tman*xd3 + 5*uman*xd3 - m2**2*(xd0 + 3*xd1 + 2*xd11 + 2*xd12 + xd13 + xd2 - xd23 + 3*xd3 - xd33) - (2*m3**2 + 2*sman + tman - uman)*xd33 + 
c     &         m1**2*(-xd0 - xd1 + 3*xd12 + 3*xd13 + 2*xd2 + 3*xd22 + 6*xd23 + xd3 + 3*xd33))))
c
c        
c        Bv0L = -(g1L*g0*g2*g3*(xc2 - 2*xd00 - tman*xd3 + m2**2*(xd13 + 2*xd23 + xd3) + (m1**2 - sman)*(xd23 + xd33) + m3**2*(xd13 + xd23 + xd33) - uman*(xd13 + xd23 + xd33)))
c        Bv0R = -(g1R*g0*g2*g3*(xc2 - 2*xd00 - tman*xd3 + m2**2*(xd13 + 2*xd23 + xd3) + (m1**2 - sman)*(xd23 + xd33) + m3**2*(xd13 + xd23 + xd33) - uman*(xd13 + xd23 + xd33)))
c        
c        Bv1L = -(g0*g2*g3*(g1R*m1*(3*xd2 + 4*xd22 + 8*xd23 + 2*xd3 + 4*xd33) + g1L*m3*(3*xd1 + 4*xd12 + 4*xd13 + 3*xd2 + 4*xd22 + 8*xd23 + 2*xd3 + 4*xd33)))
c        Bv1R =  -(g0*g2*g3*(g1L*m1*(3*xd2 + 4*xd22 + 8*xd23 + 2*xd3 + 4*xd33) + g1R*m3*(3*xd1 + 4*xd12 + 4*xd13 + 3*xd2 + 4*xd22 + 8*xd23 + 2*xd3 + 4*xd33)))
c        
c        Bv2L =   (g0*g2*g3*(g1R*m1*(xd2 - 2*(2*xd23 + xd3 + 2*xd33)) + g1L*m3*(xd1 - 4*xd13 + xd2 - 2*(2*xd23 + xd3 + 2*xd33))))
c        Bv2R =   (g0*g2*g3*(g1L*m1*(xd2 - 2*(2*xd23 + xd3 + 2*xd33)) + g1R*m3*(xd1 - 4*xd13 + xd2 - 2*(2*xd23 + xd3 + 2*xd33))))
c        
c        Bv3L = (g0*g2*g3*(g1R*m1*(4*xd12 + 4*xd13 + 3*xd2 + 4*xd22 + 8*xd23 + 2*xd3 + 4*xd33) + g1L*m3*(3*xd1 + 4*xd11 + 8*xd12 + 8*xd13 + 3*xd2 + 4*xd22 + 8*xd23 + 2*xd3 + 4*xd33)))
c        Bv3R = (g0*g2*g3*(g1L*m1*(4*xd12 + 4*xd13 + 3*xd2 + 4*xd22 + 8*xd23 + 2*xd3 + 4*xd33) + g1R*m3*(3*xd1 + 4*xd11 + 8*xd12 + 8*xd13 + 3*xd2 + 4*xd22 + 8*xd23 + 2*xd3 + 4*xd33)))


        Bs0L = g0*g2*g3*(-(g1R*m1*(xc1 + xc2 + m2**2*xd12 + m3**2*xd12 - uman*xd12 + m2**2*xd13 + m3**2*xd13 - uman*xd13 + m2**2*xd2 - tman*xd2 + m1**2*xd22 + 2*m2**2*xd22 + m3**2*xd22 - 
     &         sman*xd22 - uman*xd22 + 2*m1**2*xd23 + 2*m2**2*xd23 + 2*m3**2*xd23 - 2*sman*xd23 - 2*uman*xd23 + m2**2*xd3 - tman*xd3 + m1**2*xd33 + m3**2*xd33 - sman*xd33 - 
     &         uman*xd33)) - g1L*m3*(-xc0 - tman*xd1 + m3**2*xd11 - uman*xd11 + m1**2*xd12 + 2*m3**2*xd12 - sman*xd12 - 2*uman*xd12 + m1**2*xd13 + 2*m3**2*xd13 - sman*xd13 - 
     &       2*uman*xd13 - tman*xd2 + m1**2*xd22 + m3**2*xd22 - sman*xd22 - uman*xd22 + 2*m1**2*xd23 + 2*m3**2*xd23 - 2*sman*xd23 - 2*uman*xd23 - tman*xd3 + 
     &       m2**2*(xd1 + xd11 + 3*xd12 + xd13 + xd2 + 2*xd22 + 2*xd23 + xd3) + m1**2*xd33 + m3**2*xd33 - sman*xd33 - uman*xd33))      

        Bs0R = g0*g2*g3*(-(g1L*m1*(xc1 + xc2 + m2**2*xd12 + m3**2*xd12 - uman*xd12 + m2**2*xd13 + m3**2*xd13 - uman*xd13 + m2**2*xd2 - tman*xd2 + m1**2*xd22 + 2*m2**2*xd22 + m3**2*xd22 - 
     &         sman*xd22 - uman*xd22 + 2*m1**2*xd23 + 2*m2**2*xd23 + 2*m3**2*xd23 - 2*sman*xd23 - 2*uman*xd23 + m2**2*xd3 - tman*xd3 + m1**2*xd33 + m3**2*xd33 - sman*xd33 - 
     &         uman*xd33)) - g1R*m3*(-xc0 - tman*xd1 + m3**2*xd11 - uman*xd11 + m1**2*xd12 + 2*m3**2*xd12 - sman*xd12 - 2*uman*xd12 + m1**2*xd13 + 2*m3**2*xd13 - sman*xd13 - 
     &       2*uman*xd13 - tman*xd2 + m1**2*xd22 + m3**2*xd22 - sman*xd22 - uman*xd22 + 2*m1**2*xd23 + 2*m3**2*xd23 - 2*sman*xd23 - 2*uman*xd23 - tman*xd3 + 
     &       m2**2*(xd1 + xd11 + 3*xd12 + xd13 + xd2 + 2*xd22 + 2*xd23 + xd3) + m1**2*xd33 + m3**2*xd33 - sman*xd33 - uman*xd33))

        Bs1L = g0*g2*g3*(-(g1R*m1*m3*xd1) + g1L*(xc1 + xc2 + m1**2*xd0 + 3*m2**2*xd0 + 2*m3**2*xd0 - tman*xd0 - 3*uman*xd0 + m1**2*xd1 + 3*m2**2*xd1 + 3*m3**2*xd1 - tman*xd1 - 
     &       3*uman*xd1 + 4*m2**2*xd12 + 6*m3**2*xd12 - 4*uman*xd12 + 4*m2**2*xd13 + 6*m3**2*xd13 - 4*uman*xd13 + 7*m2**2*xd2 + 7*m3**2*xd2 - tman*xd2 - 7*uman*xd2 - m1**2*xd22 + 
     &       4*m2**2*xd22 + 5*m3**2*xd22 + tman*xd22 - 4*uman*xd22 - 2*m1**2*xd23 + 3*m2**2*xd23 + 9*m3**2*xd23 + 2*tman*xd23 - 7*uman*xd23 + m1**2*xd3 + 5*m2**2*xd3 + 
     &       7*m3**2*xd3 - sman*xd3 - tman*xd3 - 7*uman*xd3 - m1**2*xd33 - m2**2*xd33 + 4*m3**2*xd33 + tman*xd33 - 3*uman*xd33))

        Bs1R = g0*g2*g3*(-(g1L*m1*m3*xd1) + g1R*(xc1 + xc2 + m1**2*xd0 + 3*m2**2*xd0 + 2*m3**2*xd0 - tman*xd0 - 3*uman*xd0 + m1**2*xd1 + 3*m2**2*xd1 + 3*m3**2*xd1 - tman*xd1 - 
     &       3*uman*xd1 + 4*m2**2*xd12 + 6*m3**2*xd12 - 4*uman*xd12 + 4*m2**2*xd13 + 6*m3**2*xd13 - 4*uman*xd13 + 7*m2**2*xd2 + 7*m3**2*xd2 - tman*xd2 - 7*uman*xd2 - m1**2*xd22 + 
     &       4*m2**2*xd22 + 5*m3**2*xd22 + tman*xd22 - 4*uman*xd22 - 2*m1**2*xd23 + 3*m2**2*xd23 + 9*m3**2*xd23 + 2*tman*xd23 - 7*uman*xd23 + m1**2*xd3 + 5*m2**2*xd3 + 
     &       7*m3**2*xd3 - sman*xd3 - tman*xd3 - 7*uman*xd3 - m1**2*xd33 - m2**2*xd33 + 4*m3**2*xd33 + tman*xd33 - 3*uman*xd33))

        Bs2L = g0*g2*g3*(3*g1R*m1*m3*xd1 + g1L*(2*xc0 + xc2 - 3*m1**2*xd0 - m2**2*xd0 + 2*m3**2*xd0 + 3*tman*xd0 + uman*xd0 - 8*xd00 - 3*m1**2*xd1 - m2**2*xd1 + 3*m3**2*xd1 + 3*tman*xd1 + 
     &       uman*xd1 + 4*m2**2*xd13 + 6*m3**2*xd13 - 4*uman*xd13 - 2*m1**2*xd2 - m2**2*xd2 + m3**2*xd2 + 5*tman*xd2 + uman*xd2 - m1**2*xd23 + 4*m2**2*xd23 + 5*m3**2*xd23 + 
     &       tman*xd23 - 4*uman*xd23 - 5*m1**2*xd3 - m2**2*xd3 + 3*m3**2*xd3 + 3*sman*xd3 + 5*tman*xd3 - uman*xd3 - m1**2*xd33 - m2**2*xd33 + 4*m3**2*xd33 + tman*xd33 - 3*uman*xd33))

        Bs2R = g0*g2*g3*(3*g1L*m1*m3*xd1 + g1R*(2*xc0 + xc2 - 3*m1**2*xd0 - m2**2*xd0 + 2*m3**2*xd0 + 3*tman*xd0 + uman*xd0 - 8*xd00 - 3*m1**2*xd1 - m2**2*xd1 + 3*m3**2*xd1 + 3*tman*xd1 + 
     &       uman*xd1 + 4*m2**2*xd13 + 6*m3**2*xd13 - 4*uman*xd13 - 2*m1**2*xd2 - m2**2*xd2 + m3**2*xd2 + 5*tman*xd2 + uman*xd2 - m1**2*xd23 + 4*m2**2*xd23 + 5*m3**2*xd23 + 
     &       tman*xd23 - 4*uman*xd23 - 5*m1**2*xd3 - m2**2*xd3 + 3*m3**2*xd3 + 3*sman*xd3 + 5*tman*xd3 - uman*xd3 - m1**2*xd33 - m2**2*xd33 + 4*m3**2*xd33 + tman*xd33 - 3*uman*xd33))

        Bs3L = g0*g2*g3*(g1R*m1*m3*xd1 + g1L*(xc0 - m1**2*xd0 - m2**2*xd0 - 2*m3**2*xd0 - tman*xd0 + 3*uman*xd0 - 2*xd00 - m1**2*xd1 - 3*m2**2*xd1 - 5*m3**2*xd1 - tman*xd1 + 5*uman*xd1 - 
     &       2*m2**2*xd11 - 4*m3**2*xd11 + 2*uman*xd11 + 3*m1**2*xd12 - 2*m2**2*xd12 - 7*m3**2*xd12 - 2*sman*xd12 - tman*xd12 + 4*uman*xd12 + 3*m1**2*xd13 - m2**2*xd13 - 
     &       6*m3**2*xd13 - 2*sman*xd13 - tman*xd13 + 3*uman*xd13 + 2*m1**2*xd2 - m2**2*xd2 - 5*m3**2*xd2 - 2*sman*xd2 - tman*xd2 + 5*uman*xd2 + 3*m1**2*xd22 - 3*m3**2*xd22 - 
     &       2*sman*xd22 - tman*xd22 + 2*uman*xd22 + 6*m1**2*xd23 + m2**2*xd23 - 5*m3**2*xd23 - 4*sman*xd23 - 2*tman*xd23 + 3*uman*xd23 + m1**2*xd3 - 3*m2**2*xd3 - 5*m3**2*xd3 - 
     &       sman*xd3 - tman*xd3 + 5*uman*xd3 + 3*m1**2*xd33 + m2**2*xd33 - 2*m3**2*xd33 - 2*sman*xd33 - tman*xd33 + uman*xd33))

        Bs3R = g0*g2*g3*(g1L*m1*m3*xd1 + g1R*(xc0 - m1**2*xd0 - m2**2*xd0 - 2*m3**2*xd0 - tman*xd0 + 3*uman*xd0 - 2*xd00 - m1**2*xd1 - 3*m2**2*xd1 - 5*m3**2*xd1 - tman*xd1 + 5*uman*xd1 - 
     &       2*m2**2*xd11 - 4*m3**2*xd11 + 2*uman*xd11 + 3*m1**2*xd12 - 2*m2**2*xd12 - 7*m3**2*xd12 - 2*sman*xd12 - tman*xd12 + 4*uman*xd12 + 3*m1**2*xd13 - m2**2*xd13 - 
     &       6*m3**2*xd13 - 2*sman*xd13 - tman*xd13 + 3*uman*xd13 + 2*m1**2*xd2 - m2**2*xd2 - 5*m3**2*xd2 - 2*sman*xd2 - tman*xd2 + 5*uman*xd2 + 3*m1**2*xd22 - 3*m3**2*xd22 - 
     &       2*sman*xd22 - tman*xd22 + 2*uman*xd22 + 6*m1**2*xd23 + m2**2*xd23 - 5*m3**2*xd23 - 4*sman*xd23 - 2*tman*xd23 + 3*uman*xd23 + m1**2*xd3 - 3*m2**2*xd3 - 5*m3**2*xd3 - 
     &       sman*xd3 - tman*xd3 + 5*uman*xd3 + 3*m1**2*xd33 + m2**2*xd33 - 2*m3**2*xd33 - 2*sman*xd33 - tman*xd33 + uman*xd33))     

        Bv0L = -(g0*g1L*g2*g3*(xc2 - 6*xd00 + m2**2*xd13 + m3**2*xd13 - uman*xd13 + m1**2*xd23 + 2*m2**2*xd23 + m3**2*xd23 - sman*xd23 - uman*xd23 + m2**2*xd3 - tman*xd3 + m1**2*xd33 + 
     &      m3**2*xd33 - sman*xd33 - uman*xd33))

        Bv0R = -(g0*g1R*g2*g3*(xc2 - 6*xd00 + m2**2*xd13 + m3**2*xd13 - uman*xd13 + m1**2*xd23 + 2*m2**2*xd23 + m3**2*xd23 - sman*xd23 - uman*xd23 + m2**2*xd3 - tman*xd3 + m1**2*xd33 + 
     &      m3**2*xd33 - sman*xd33 - uman*xd33))

        Bv1L = g0*g2*g3*(-(g1R*m1*(3*xd2 + 4*xd22 + 8*xd23 + 2*xd3 + 4*xd33)) - g1L*m3*(3*xd1 + 4*xd12 + 4*xd13 + 3*xd2 + 4*xd22 + 8*xd23 + 2*xd3 + 4*xd33))     

        Bv1R = g0*g2*g3*(-(g1L*m1*(3*xd2 + 4*xd22 + 8*xd23 + 2*xd3 + 4*xd33)) - g1R*m3*(3*xd1 + 4*xd12 + 4*xd13 + 3*xd2 + 4*xd22 + 8*xd23 + 2*xd3 + 4*xd33))

        Bv2L = g0*g2*g3*(g1L*m3*(xd1 - 4*xd13 + xd2 - 4*xd23 - 2*xd3 - 4*xd33) + g1R*m1*(xd2 - 2*(2*xd23 + xd3 + 2*xd33)))

        Bv2R = g0*g2*g3*(g1R*m3*(xd1 - 4*xd13 + xd2 - 4*xd23 - 2*xd3 - 4*xd33) + g1L*m1*(xd2 - 2*(2*xd23 + xd3 + 2*xd33)))
        
        Bv3L = g0*g2*g3*(g1R*m1*(4*xd12 + 4*xd13 + 3*xd2 + 4*xd22 + 8*xd23 + 2*xd3 + 4*xd33) + g1L*m3*(3*xd1 + 4*xd11 + 8*xd12 + 8*xd13 + 3*xd2 + 4*xd22 + 8*xd23 + 2*xd3 + 4*xd33))

        Bv3R = g0*g2*g3*(g1L*m1*(4*xd12 + 4*xd13 + 3*xd2 + 4*xd22 + 8*xd23 + 2*xd3 + 4*xd33) + g1R*m3*(3*xd1 + 4*xd11 + 8*xd12 + 8*xd13 + 3*xd2 + 4*xd22 + 8*xd23 + 2*xd3 + 4*xd33))        

c        print*,'------------ Julia Gluon Box -----------------------------'
c        print*,'Bs0L=',Bs0L
c        print*,'Bs0R=',Bs0R
c        print*,'Bs1L=',Bs1L
c        print*,'Bs1R=',Bs1R
c        print*,'Bs2L=',Bs2L
c        print*,'Bs2R=',Bs2R
c        print*,'Bs3L=',Bs3L
c        print*,'Bs3R=',Bs3R
c        print*,'------------'
c        print*,'Bv0L=',Bv0L
c        print*,'Bv0R=',Bv0R
c        print*,'Bv1L=',Bv1L
c        print*,'Bv1R=',Bv1R
c        print*,'Bv2L=',Bv2L
c        print*,'Bv2R=',Bv2R
c        print*,'Bv3L=',Bv3L
c        print*,'Bv3R=',Bv3R
c
c        print*,'------------ Karol Gluon Box -----------------------------'
c        print*,'Bs0L=',Bs0LK
c        print*,'Bs0R=',Bs0RK
c        print*,'Bs1L=',Bs1LK
c        print*,'Bs1R=',Bs1RK
c        print*,'Bs2L=',Bs2LK
c        print*,'Bs2R=',Bs2RK
c        print*,'Bs3L=',Bs3LK
c        print*,'Bs3R=',Bs3RK
c        print*,'------------'
c        print*,'Bv0L=',Bv0LK
c        print*,'Bv0R=',Bv0RK
c        print*,'Bv1L=',Bv1LK
c        print*,'Bv1R=',Bv1RK
c        print*,'Bv2L=',Bv2LK
c        print*,'Bv2R=',Bv2RK
c        print*,'Bv3L=',Bv3LK
c        print*,'Bv3R=',Bv3RK

        
! --- Interference with s channel tree-level --

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))
        
        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0
        
        boxGluonST = CF*NC*(-CV/2d0)*2d0*MsBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A1Lc,A1Rc,BLc,BRc,mq1c,xv)

! 
! ! --- Interference with t channel tree-level --
               
        A2Lc = dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen,sfmeigen,sftype,sfgen,xtype))
 
        mq2c = MSf(ptmeigen,pttype,ptgen)
c        wquark = width(pstype)
        wquark = 0d0
 
        boxGluonST = boxGluonST +  CF*NC*(-CV/2d0)*2d0*MtBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A2Lc,A2Rc,CLc,mq2c,xv)

     
         NeuQ2qx_M2BoxG = NeuQ2qx_M2BoxG + boxGluonST

        
        
! ----- ST box with gluino (2)-----

c        Bs0L = 0d0
c        Bs0R = 0d0
c        Bs1L = 0d0
c        Bs1R = 0d0
c        Bs2L = 0d0
c        Bs2R = 0d0
c        Bs3L = 0d0
c        Bs3R = 0d0
c        
c        Bv0L = 0d0
c        Bv0R = 0d0
c        Bv1L = 0d0
c        Bv1R = 0d0
c        Bv2L = 0d0
c        Bv2R = 0d0
c        Bv3L = 0d0
c        Bv3R = 0d0
c
c        do k=1,2
c      
c
c        g0L = JR(k,ftype,fgen)
c        g0R = JL(k,ftype,fgen)
c        g1L = AL(chimeigen,k,sftype,sfgen)
c        g1R = AR(chimeigen,k,sftype,sfgen)
c        g2L = JL(sfmeigen,sftype,sfgen)
c        g2R = JR(sfmeigen,sftype,sfgen)
c        g3 = -GS
c
c        Mloop0 = Mf(ftype,fgen)
c        Mloop1 = MSf(k,ftype,fgen)
c        Mloop2 = MGl
c        Mloop3 = MGl
c
c     
c       call myDget(m2**2,sman,m3**2,tman,m1**2,0d0,Mf(ftype,fgen)**2,MGl**2,MSf(k,ftype,fgen)**2,MGl**2)
c       
c
c       call myCget(m2**2,tman,0d0,MGl**2,Mf(ftype,fgen)**2,MGl**2)
c       xc0a = xc0
c       xc1a = xc1
c       xc2a = xc2
c
c       call myCget(m3**2,m1**2,tman,MGl**2,MSf(k,ftype,fgen)**2,Mf(ftype,fgen)**2)
c       xc0c = xc0
c       xc1c = xc1
c       xc2c = xc2
c       
c       call myCget(sman,0d0,m3**2,MSf(k,ftype,fgen)**2,MGl**2,MGl**2) 
c
c 
c        Bs0L =  Bs0L + (g3*(g0L*m3*(g1L*(g2L*Mloop0*(-(Mloop2*xd1) + Mloop3*(xd0 + xd1)) + 
c     &            g2R*(xc0 - xc0c - xc1c - xc2c + Mloop0**2*xd0 + 
c     &               (m2**2 - m4**2 - Mloop2*Mloop3 + tman)*xd1 + Mloop3**2*xd2)) - 
c     &         g1R*m1*(g2R*Mloop0 + g2L*Mloop3)*xd3) + 
c     &      g0R*(g1R*m1*(g2L*(xc0 - xc0c - xc2c + Mloop0**2*xd0 - m4**2*xd1 + 
c     &               (Mloop3**2 + tman)*xd1 + Mloop3**2*xd2 + (m1**2 + m2**2 - sman)*xd3 - 
c     &               Mloop2*Mloop3*(xd1 + xd3)) + 
c     &            g2R*Mloop0*(-(Mloop2*(xd1 + xd3)) + Mloop3*(xd0 + xd1 + xd3))) + 
c     &         g1L*(g2L*Mloop0*(xc0 - xc0a - xc0c + 
c     &               (Mloop0**2 - Mloop1**2 + (Mloop2 - Mloop3)*Mloop3 + sman)*xd0 + 2*xd00 - 
c     &               uman*xd1 + m2**2*(xd1 + xd2) + m1**2*(xd1 + xd3)) + 
c     &            g2R*(Mloop2*(xc0 + Mloop0**2*xd0 + (m3**2 - uman)*xd1 + m2**2*(xd1 + xd2)) + 
c     &               Mloop3*(-xc0a - Mloop1**2*xd0 + 2*xd00 - m3**2*xd1 + m1**2*(xd0 + xd1 + xd3))))))) 
c        Bs0R = Bs0R + (g3*(g0R*m3*(g1R*(g2R*Mloop0*(-(Mloop2*xd1) + Mloop3*(xd0 + xd1)) + 
c     &            g2L*(xc0 - xc0c - xc1c - xc2c + Mloop0**2*xd0 + 
c     &               (m2**2 - m4**2 - Mloop2*Mloop3 + tman)*xd1 + Mloop3**2*xd2)) - 
c     &         g1L*m1*(g2L*Mloop0 + g2R*Mloop3)*xd3) + 
c     &      g0L*(g1L*m1*(g2R*(xc0 - xc0c - xc2c + Mloop0**2*xd0 - m4**2*xd1 + 
c     &               (Mloop3**2 + tman)*xd1 + Mloop3**2*xd2 + (m1**2 + m2**2 - sman)*xd3 - 
c     &               Mloop2*Mloop3*(xd1 + xd3)) + 
c     &            g2L*Mloop0*(-(Mloop2*(xd1 + xd3)) + Mloop3*(xd0 + xd1 + xd3))) + 
c     &         g1R*(g2R*Mloop0*(xc0 - xc0a - xc0c + 
c     &               (Mloop0**2 - Mloop1**2 + (Mloop2 - Mloop3)*Mloop3 + sman)*xd0 + 2*xd00 - 
c     &               uman*xd1 + m2**2*(xd1 + xd2) + m1**2*(xd1 + xd3)) + 
c     &            g2L*(Mloop2*(xc0 + Mloop0**2*xd0 + (m3**2 - uman)*xd1 + m2**2*(xd1 + xd2)) + 
c     &               Mloop3*(-xc0a - Mloop1**2*xd0 + 2*xd00 - m3**2*xd1 + m1**2*(xd0 + xd1 + xd3)))))))
c        
c        Bs1L =  Bs1L - (g3*(g0R*g1L*m3*(g2L*Mloop0 + g2R*Mloop3)*(xd2 + xd23 + xd33) + 
c     &       g0L*(-(g1R*m1*(g2R*Mloop0 + g2L*Mloop3)*(xd2 + xd22 + 2*xd23 + xd3 + xd33)) + 
c     &          g1L*(g2L*Mloop0*(Mloop2*(-xd1 + xd2) + Mloop3*(xd0 + xd1 + xd3)) + 
c     &             g2R*(-xc1 + (-m4**2 - Mloop2*Mloop3 + Mloop3**2 + tman)*xd1 - m1**2*xd2 + 
c     &                Mloop2*Mloop3*xd2 + sman*xd2 - Mloop3**2*xd3 + Mloop0**2*(xd0 + xd2 + xd3) + 
c     &                (m3**2 - uman)*(xd23 + xd33) + m2**2*(xd12 + xd13 + xd23 + xd3 + xd33))))))*2d0
c        Bs1R = Bs1R - (g3*(g0L*g1R*m3*(g2R*Mloop0 + g2L*Mloop3)*(xd2 + xd23 + xd33) + 
c     &       g0R*(-(g1L*m1*(g2L*Mloop0 + g2R*Mloop3)*(xd2 + xd22 + 2*xd23 + xd3 + xd33)) + 
c     &          g1R*(g2R*Mloop0*(Mloop2*(-xd1 + xd2) + Mloop3*(xd0 + xd1 + xd3)) + 
c     &             g2L*(-xc1 + (-m4**2 - Mloop2*Mloop3 + Mloop3**2 + tman)*xd1 - m1**2*xd2 + 
c     &                Mloop2*Mloop3*xd2 + sman*xd2 - Mloop3**2*xd3 + Mloop0**2*(xd0 + xd2 + xd3) + 
c     &                (m3**2 - uman)*(xd23 + xd33) + m2**2*(xd12 + xd13 + xd23 + xd3 + xd33))))))*2d0
c        
c        Bs2L = Bs2L + (g3*(g0R*g1L*m3*(g2L*Mloop0 + g2R*Mloop3)*(xd1 + xd13 - xd2) + 
c     &      g0L*(-(g1R*m1*(g2R*Mloop0 + g2L*Mloop3)*(xd1 + xd12 + xd13 - xd2)) + 
c     &         g1L*(g2L*Mloop0*(Mloop2*(xd1 - xd2) + Mloop3*xd2) + 
c     &            g2R*(xc1 + 2*xd00 + Mloop0**2*xd1 + Mloop2*Mloop3*xd1 - Mloop3**2*xd1 + sman*xd1 + 
c     &               (m3**2 - uman)*xd13 - Mloop2*Mloop3*xd2 - sman*xd2 + m1**2*(-xd1 + xd2) + 
c     &               m2**2*(xd11 + xd13 + xd2 - xd3) + Mloop3**2*xd3)))))*2d0
c        Bs2R = Bs2R + (g3*(g0L*g1R*m3*(g2R*Mloop0 + g2L*Mloop3)*(xd1 + xd13 - xd2) + 
c     &      g0R*(-(g1L*m1*(g2L*Mloop0 + g2R*Mloop3)*(xd1 + xd12 + xd13 - xd2)) + 
c     &         g1R*(g2R*Mloop0*(Mloop2*(xd1 - xd2) + Mloop3*xd2) + 
c     &            g2L*(xc1 + 2*xd00 + Mloop0**2*xd1 + Mloop2*Mloop3*xd1 - Mloop3**2*xd1 + sman*xd1 + 
c     &               (m3**2 - uman)*xd13 - Mloop2*Mloop3*xd2 - sman*xd2 + m1**2*(-xd1 + xd2) + 
c     &               m2**2*(xd11 + xd13 + xd2 - xd3) + Mloop3**2*xd3)))))*2d0
c        
c        Bs3L =  Bs3L + (g3*(g0R*g1L*m3*(g2L*Mloop0 + g2R*Mloop3)*(xd3 + xd33) + 
c     &      g0L*(-(g1R*m1*(g2R*Mloop0 + g2L*Mloop3)*(xd23 + xd3 + xd33)) + 
c     &         g1L*(g2L*Mloop0*Mloop2*xd3 + 
c     &            g2R*(xc0c + xc1c + xc2 + xc2c + 
c     &               (-m1**2 + Mloop0**2 + (Mloop2 - Mloop3)*Mloop3 + sman)*xd3 + 
c     &               (m3**2 - uman)*xd33 + m2**2*(xd13 + xd33))))))*2d0
c        Bs3R =  Bs3R + (g3*(g0L*g1R*m3*(g2R*Mloop0 + g2L*Mloop3)*(xd3 + xd33) + 
c     &      g0R*(-(g1L*m1*(g2L*Mloop0 + g2R*Mloop3)*(xd23 + xd3 + xd33)) + 
c     &         g1R*(g2R*Mloop0*Mloop2*xd3 + 
c     &            g2L*(xc0c + xc1c + xc2 + xc2c + 
c     &               (-m1**2 + Mloop0**2 + (Mloop2 - Mloop3)*Mloop3 + sman)*xd3 + 
c     &               (m3**2 - uman)*xd33 + m2**2*(xd13 + xd33))))))*2d0
c        
c        Bv0L = Bv0L + (g3*(g0R*m3*(g1L*(g2L*Mloop0*(xd0 + xd1 + xd2) + g2R*(Mloop2*xd1 + Mloop3*xd2)) - 
c     &         g1R*g2L*m1*xd3) + g0L*(-(g1R*m1*
c     &            (g2R*Mloop0*(xd0 + xd1 + xd2 + xd3) + g2L*(Mloop3*xd2 + Mloop2*(xd1 + xd3)))) + 
c     &         g1L*(g2L*Mloop0*(-(Mloop3*xd2) + Mloop2*(xd0 + xd2)) + 
c     &            g2R*(xc0a + Mloop1**2*xd0 - 2*xd00 + m3**2*xd1 + Mloop2*Mloop3*xd2 + sman*xd2 - 
c     &               Mloop3**2*xd3 + m2**2*(-xd2 + xd3) - m1**2*(xd0 + xd1 + xd2 + xd3))))))
c        Bv0R =  Bv0R + (g3*(g0L*m3*(g1R*(g2R*Mloop0*(xd0 + xd1 + xd2) + g2L*(Mloop2*xd1 + Mloop3*xd2)) - 
c     &         g1L*g2R*m1*xd3) + g0R*(-(g1L*m1*
c     &            (g2L*Mloop0*(xd0 + xd1 + xd2 + xd3) + g2R*(Mloop3*xd2 + Mloop2*(xd1 + xd3)))) + 
c     &         g1R*(g2R*Mloop0*(-(Mloop3*xd2) + Mloop2*(xd0 + xd2)) + 
c     &            g2L*(xc0a + Mloop1**2*xd0 - 2*xd00 + m3**2*xd1 + Mloop2*Mloop3*xd2 + sman*xd2 - 
c     &               Mloop3**2*xd3 + m2**2*(-xd2 + xd3) - m1**2*(xd0 + xd1 + xd2 + xd3))))))
c
c        
c        Bv1L =   Bv1L - (g3*(-(g0L*g1L*g2R*m3*(xd2 + xd23 + xd33)) + 
c     &       g0R*(g1L*(g2R*(Mloop2*(xd1 - xd2) + Mloop3*(xd12 + xd13 + xd2)) + 
c     &             g2L*Mloop0*(xd0 + xd1 + xd12 + xd13 + xd2 + xd3)) - 
c     &          g1R*g2L*m1*(xd2 + xd22 + 2*xd23 + xd3 + xd33))))*2d0
c        Bv1R =  Bv1R - (g3*(-(g0R*g1R*g2L*m3*(xd2 + xd23 + xd33)) + 
c     &       g0L*(g1R*(g2L*(Mloop2*(xd1 - xd2) + Mloop3*(xd12 + xd13 + xd2)) + 
c     &             g2R*Mloop0*(xd0 + xd1 + xd12 + xd13 + xd2 + xd3)) - 
c     &          g1L*g2R*m1*(xd2 + xd22 + 2*xd23 + xd3 + xd33))))*2d0
c        
c        Bv2L = Bv2L - (g3*(g0L*g1L*g2R*m3*(xd1 + xd13 - xd2) + 
c     &       g0R*(g1R*g2L*m1*(xd1 + xd12 + xd13 - xd2) - 
c     &          g1L*(g2L*Mloop0*(xd11 + xd2) + g2R*(-(Mloop2*xd1) + Mloop3*xd11 + Mloop2*xd2)))))*2d0
c        Bv2R = Bv2R - (g3*(g0R*g1R*g2L*m3*(xd1 + xd13 - xd2) + 
c     &       g0L*(g1L*g2R*m1*(xd1 + xd12 + xd13 - xd2) - 
c     &          g1R*(g2R*Mloop0*(xd11 + xd2) + g2L*(-(Mloop2*xd1) + Mloop3*xd11 + Mloop2*xd2)))))*2d0
c        
c        Bv3L = Bv3L - (g3*(g0L*g1L*g2R*m3*(xd3 + xd33) + 
c     &       g0R*(-(g1L*(g2L*Mloop0 + g2R*Mloop3)*xd13) + g1L*g2R*Mloop2*xd3 + 
c     &          g1R*g2L*m1*(xd23 + xd3 + xd33))))*2d0
c        Bv3R = Bv3R - (g3*(g0R*g1R*g2L*m3*(xd3 + xd33) + 
c     &       g0L*(-(g1R*(g2R*Mloop0 + g2L*Mloop3)*xd13) + g1R*g2L*Mloop2*xd3 + 
c     &          g1L*g2R*m1*(xd23 + xd3 + xd33))))*2d0
c
c
c        enddo
c
c      print*,'------------ Julia Gluino Box -----------------------------'
c      print*,'Bs0L=',Bs0L
c      print*,'Bs0R=',Bs0R
c      print*,'Bs1L=',Bs1L
c      print*,'Bs1R=',Bs1R
c      print*,'Bs2L=',Bs2L
c      print*,'Bs2R=',Bs2R
c      print*,'Bs3L=',Bs3L
c      print*,'Bs3R=',Bs3R
c      print*,'------------'
c      print*,'Bv0L=',Bv0L
c      print*,'Bv0R=',Bv0R
c      print*,'Bv1L=',Bv1L
c      print*,'Bv1R=',Bv1R
c      print*,'Bv2L=',Bv2L
c      print*,'Bv2R=',Bv2R
c      print*,'Bv3L=',Bv3L
c      print*,'Bv3R=',Bv3R

c Karol's version
        Bs0L = 0d0
        Bs0R = 0d0
        Bs1L = 0d0
        Bs1R = 0d0
        Bs2L = 0d0
        Bs2R = 0d0
        Bs3L = 0d0
        Bs3R = 0d0
        
        Bv0L = 0d0
        Bv0R = 0d0
        Bv1L = 0d0
        Bv1R = 0d0
        Bv2L = 0d0
        Bv2R = 0d0
        Bv3L = 0d0
        Bv3R = 0d0
        
        do k=1,2
        
            g0L = JR(k,ftype,fgen)
            g0R = JL(k,ftype,fgen)
            g1L = dconjg(AR(chimeigen,k,sftype,sfgen))
            g1R = dconjg(AL(chimeigen,k,sftype,sfgen))
            g2L = JR(sfmeigen,sftype,sfgen)
            g2R = JL(sfmeigen,sftype,sfgen)
            g3L = GS
            g3R = GS
            
            call myDget(m3**2, m1**2, m2**2, 0d0, tman, sman, MGl2, MSf(k,ftype,fgen)**2, Mf(sftype,sfgen)**2, MGl2)
            call myCget(m1**2,m2**2,sman,MSf(k,ftype,fgen)**2, Mf(sftype,sfgen)**2, MGl2)

            Bs0L = Bs0L + g0L*m3*(-(g1R*g3L*m1*(g2R*m3 + g2L*MGl)*xd1) + g1L*(g2L*m3*MGl*(-(g3L*(xd1 + xd2 + xd3)) + g3R*(xd0 + xd1 + xd2 + xd3)) + 
     &            g2R*(g3R*MGl**2*(xd0 + xd1 + xd2 + xd3) + g3L*(xc0 - 2*xd00 + m2**2*xd1 + m3**2*xd1 - tman*xd1 - uman*xd1 + m1**2*xd11 + m2**2*xd11 - tman*xd11 - uman*xd11 + 
     &              2*m1**2*xd12 + 3*m2**2*xd12 + m3**2*xd12 - sman*xd12 - 3*tman*xd12 - 2*uman*xd12 + 3*m1**2*xd13 + 3*m2**2*xd13 + 2*m3**2*xd13 - 2*sman*xd13 - 3*tman*xd13 - 
     &              3*uman*xd13 + m1**2*xd2 + 2*m2**2*xd2 + m3**2*xd2 - sman*xd2 - tman*xd2 - uman*xd2 + m1**2*xd22 + 2*m2**2*xd22 + m3**2*xd22 - sman*xd22 - 2*tman*xd22 - uman*xd22 + 
     &              3*m1**2*xd23 + 4*m2**2*xd23 + 3*m3**2*xd23 - 3*sman*xd23 - 4*tman*xd23 - 3*uman*xd23 + m1**2*xd3 + m2**2*xd3 + m3**2*xd3 - sman*xd3 - tman*xd3 - uman*xd3 + 
     &              2*m1**2*xd33 + 2*m2**2*xd33 + 2*m3**2*xd33 - 2*sman*xd33 - 2*tman*xd33 - 2*uman*xd33)))) - 
     &            g0R*(g1L*(g2L*m3*(-(g3L*MGl**2*xd0) + g3R*(xc0 - 2*xd00 + m3**2*xd1 - m1**2*xd2 - m2**2*xd2 + sman*xd2 + tman*xd2 - m1**2*xd3 - m2**2*xd3 - m3**2*xd3 + sman*xd3 + tman*xd3 + 
     &              uman*xd3)) + g2R*MGl*(g3R*(xc0 - 2*xd00 + m3**2*xd1 + tman*xd2 - m2**2*xd3 - m3**2*xd3 + tman*xd3 + uman*xd3) - 
     &           g3L*(xc0 - m3**2*xd0 + tman*xd0 + uman*xd0 - m1**2*xd1 + tman*xd1 + uman*xd1 - m3**2*xd2 + 2*tman*xd2 + uman*xd2 - m1**2*xd3 - 2*m3**2*xd3 + sman*xd3 + 2*tman*xd3 + 
     &              2*uman*xd3 - m2**2*(xd0 + xd1 + xd2 + 2*xd3)))) - g1R*m1*
     &          (g2R*m3*MGl*(-(g3R*(xd2 + xd3)) + g3L*(xd0 + xd2 + xd3)) + g2L*
     &          (g3L*MGl**2*(xd0 + xd2 + xd3) - g3R*(xc1 + xc2 + 2*xd00 - m1**2*xd12 - m2**2*xd12 + tman*xd12 + uman*xd12 - m1**2*xd13 - m2**2*xd13 + tman*xd13 + uman*xd13 - 
     &              m1**2*xd2 - 2*m2**2*xd2 - m3**2*xd2 + sman*xd2 + tman*xd2 + uman*xd2 - m1**2*xd22 - 2*m2**2*xd22 - m3**2*xd22 + sman*xd22 + 2*tman*xd22 + uman*xd22 - 
     &              3*m1**2*xd23 - 4*m2**2*xd23 - 3*m3**2*xd23 + 3*sman*xd23 + 4*tman*xd23 + 3*uman*xd23 - m1**2*xd3 - m2**2*xd3 - m3**2*xd3 + sman*xd3 + tman*xd3 + uman*xd3 - 
     &              2*m1**2*xd33 - 2*m2**2*xd33 - 2*m3**2*xd33 + 2*sman*xd33 + 2*tman*xd33 + 2*uman*xd33))))

            Bs0R = Bs0R + g0R*m3*(-(g1L*g3R*m1*(g2L*m3 + g2R*MGl)*xd1) + g1R*(g2R*m3*MGl*(-(g3R*(xd1 + xd2 + xd3)) + g3L*(xd0 + xd1 + xd2 + xd3)) + 
     &           g2L*(g3L*MGl**2*(xd0 + xd1 + xd2 + xd3) + g3R*(xc0 - 2*xd00 + m2**2*xd1 + m3**2*xd1 - tman*xd1 - uman*xd1 + m1**2*xd11 + m2**2*xd11 - tman*xd11 - uman*xd11 + 
     &              2*m1**2*xd12 + 3*m2**2*xd12 + m3**2*xd12 - sman*xd12 - 3*tman*xd12 - 2*uman*xd12 + 3*m1**2*xd13 + 3*m2**2*xd13 + 2*m3**2*xd13 - 2*sman*xd13 - 3*tman*xd13 - 
     &              3*uman*xd13 + m1**2*xd2 + 2*m2**2*xd2 + m3**2*xd2 - sman*xd2 - tman*xd2 - uman*xd2 + m1**2*xd22 + 2*m2**2*xd22 + m3**2*xd22 - sman*xd22 - 2*tman*xd22 - uman*xd22 + 
     &              3*m1**2*xd23 + 4*m2**2*xd23 + 3*m3**2*xd23 - 3*sman*xd23 - 4*tman*xd23 - 3*uman*xd23 + m1**2*xd3 + m2**2*xd3 + m3**2*xd3 - sman*xd3 - tman*xd3 - uman*xd3 + 
     &              2*m1**2*xd33 + 2*m2**2*xd33 + 2*m3**2*xd33 - 2*sman*xd33 - 2*tman*xd33 - 2*uman*xd33)))) - 
     &           g0L*(g1R*(g2R*m3*(-(g3R*MGl**2*xd0) + g3L*(xc0 - 2*xd00 + m3**2*xd1 - m1**2*xd2 - m2**2*xd2 + sman*xd2 + tman*xd2 - m1**2*xd3 - m2**2*xd3 - m3**2*xd3 + sman*xd3 + tman*xd3 + 
     &              uman*xd3)) + g2L*MGl*(g3L*(xc0 - 2*xd00 + m3**2*xd1 + tman*xd2 - m2**2*xd3 - m3**2*xd3 + tman*xd3 + uman*xd3) - 
     &           g3R*(xc0 - m3**2*xd0 + tman*xd0 + uman*xd0 - m1**2*xd1 + tman*xd1 + uman*xd1 - m3**2*xd2 + 2*tman*xd2 + uman*xd2 - m1**2*xd3 - 2*m3**2*xd3 + sman*xd3 + 2*tman*xd3 + 
     &              2*uman*xd3 - m2**2*(xd0 + xd1 + xd2 + 2*xd3)))) - g1L*m1*
     &         (g2L*m3*MGl*(-(g3L*(xd2 + xd3)) + g3R*(xd0 + xd2 + xd3)) + g2R*
     &         (g3R*MGl**2*(xd0 + xd2 + xd3) - g3L*(xc1 + xc2 + 2*xd00 - m1**2*xd12 - m2**2*xd12 + tman*xd12 + uman*xd12 - m1**2*xd13 - m2**2*xd13 + tman*xd13 + uman*xd13 - 
     &              m1**2*xd2 - 2*m2**2*xd2 - m3**2*xd2 + sman*xd2 + tman*xd2 + uman*xd2 - m1**2*xd22 - 2*m2**2*xd22 - m3**2*xd22 + sman*xd22 + 2*tman*xd22 + uman*xd22 - 
     &              3*m1**2*xd23 - 4*m2**2*xd23 - 3*m3**2*xd23 + 3*sman*xd23 + 4*tman*xd23 + 3*uman*xd23 - m1**2*xd3 - m2**2*xd3 - m3**2*xd3 + sman*xd3 + tman*xd3 + uman*xd3 - 
     &              2*m1**2*xd33 - 2*m2**2*xd33 - 2*m3**2*xd33 + 2*sman*xd33 + 2*tman*xd33 + 2*uman*xd33))))

           Bs1L = Bs1L + 2*(-(g0R*g1L*g3R*m3*(g2L*m3 + g2R*MGl)*(xd12 + xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33)) + 
     &          g0L*g3L*(g1R*m1*(g2R*m3 + g2L*MGl)*(xd2 + xd22 + 2*xd23 + xd3 + xd33) + 
     &          g1L*(g2L*m3*MGl*(xd2 + xd3) + g2R*(xc1 + xc2 + 2*xd00 - m1**2*xd12 - m2**2*xd12 + tman*xd12 + uman*xd12 - m1**2*xd13 - m2**2*xd13 + tman*xd13 + uman*xd13 - m2**2*xd2 - 
     &             m3**2*xd2 + tman*xd2 + uman*xd2 - m2**2*xd22 - m3**2*xd22 + 2*tman*xd22 + uman*xd22 - m1**2*xd23 - 3*m2**2*xd23 - 3*m3**2*xd23 + sman*xd23 + 4*tman*xd23 + 
     &             3*uman*xd23 - m2**2*xd3 - m3**2*xd3 + tman*xd3 + uman*xd3 - m1**2*xd33 - 2*m2**2*xd33 - 2*m3**2*xd33 + sman*xd33 + 2*tman*xd33 + 2*uman*xd33))))

           Bs1R = Bs1R + 2*(-(g0L*g1R*g3L*m3*(g2R*m3 + g2L*MGl)*(xd12 + xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33)) + 
     &          g0R*g3R*(g1L*m1*(g2L*m3 + g2R*MGl)*(xd2 + xd22 + 2*xd23 + xd3 + xd33) + 
     &          g1R*(g2R*m3*MGl*(xd2 + xd3) + g2L*(xc1 + xc2 + 2*xd00 - m1**2*xd12 - m2**2*xd12 + tman*xd12 + uman*xd12 - m1**2*xd13 - m2**2*xd13 + tman*xd13 + uman*xd13 - m2**2*xd2 - 
     &             m3**2*xd2 + tman*xd2 + uman*xd2 - m2**2*xd22 - m3**2*xd22 + 2*tman*xd22 + uman*xd22 - m1**2*xd23 - 3*m2**2*xd23 - 3*m3**2*xd23 + sman*xd23 + 4*tman*xd23 + 
     &             3*uman*xd23 - m2**2*xd3 - m3**2*xd3 + tman*xd3 + uman*xd3 - m1**2*xd33 - 2*m2**2*xd33 - 2*m3**2*xd33 + sman*xd33 + 2*tman*xd33 + 2*uman*xd33))))

           Bs2L = Bs2L + 2*(-(g0R*g1L*g3R*m3*(g2L*m3 + g2R*MGl)*(xd13 + xd23 + xd3 + xd33)) + 
     &          g0L*g3L*(g1R*m1*(g2R*m3 + g2L*MGl)*(xd23 + xd3 + xd33) + g1L*(g2L*m3*MGl*xd3 + 
     &          g2R*(xc2 + 2*xd00 - m1**2*xd13 - m2**2*xd13 + tman*xd13 + uman*xd13 - m2**2*xd23 - m3**2*xd23 + 2*tman*xd23 + uman*xd23 - m2**2*xd3 - m3**2*xd3 + tman*xd3 + uman*xd3 - 
     &             m1**2*xd33 - 2*m2**2*xd33 - 2*m3**2*xd33 + sman*xd33 + 2*tman*xd33 + 2*uman*xd33))))

           Bs2R = Bs2R + 2*(-(g0L*g1R*g3L*m3*(g2R*m3 + g2L*MGl)*(xd13 + xd23 + xd3 + xd33)) + 
     &          g0R*g3R*(g1L*m1*(g2L*m3 + g2R*MGl)*(xd23 + xd3 + xd33) + g1R*(g2R*m3*MGl*xd3 + 
     &          g2L*(xc2 + 2*xd00 - m1**2*xd13 - m2**2*xd13 + tman*xd13 + uman*xd13 - m2**2*xd23 - m3**2*xd23 + 2*tman*xd23 + uman*xd23 - m2**2*xd3 - m3**2*xd3 + tman*xd3 + uman*xd3 - 
     &             m1**2*xd33 - 2*m2**2*xd33 - 2*m3**2*xd33 + sman*xd33 + 2*tman*xd33 + 2*uman*xd33))))

           Bs3L = Bs3L - 2*(-(g0R*g1L*g3R*m3*(g2L*m3 + g2R*MGl)*(xd1 + xd11 + 2*xd12 + 2*xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33)) + 
     &          g0L*(g1R*g3L*m1*(g2R*m3 + g2L*MGl)*(xd12 + xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33) + 
     &          g1L*(g2L*g3R*m3*MGl*(xd0 + xd1 + xd2 + xd3) + g2R*(g3R*MGl**2*(xd0 + xd1 + xd2 + xd3) + 
     &             g3L*(m2**2*(xd12 + xd2 + xd22 + xd23) + m1**2*(xd12 + xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33) - sman*(xd12 + xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33))))))
     
           Bs3R = Bs3R - 2*(-(g0L*g1R*g3L*m3*(g2R*m3 + g2L*MGl)*(xd1 + xd11 + 2*xd12 + 2*xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33)) + 
     &          g0R*(g1L*g3R*m1*(g2L*m3 + g2R*MGl)*(xd12 + xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33) + 
     &          g1R*(g2R*g3L*m3*MGl*(xd0 + xd1 + xd2 + xd3) + g2L*(g3L*MGl**2*(xd0 + xd1 + xd2 + xd3) + 
     &             g3R*(m2**2*(xd12 + xd2 + xd22 + xd23) + m1**2*(xd12 + xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33) - sman*(xd12 + xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33)))))) 

          Bv0L = Bv0L - (g0R*m3*(g1R*g2L*g3R*m1*xd1 + g1L*(g2L*g3R*m3*(xd1 + xd2) + g2R*MGl*(-(g3R*xd3) + g3L*(xd0 + xd1 + xd2 + xd3))))) + 
     &          g0L*(g1R*m1*(g2R*g3L*m3*xd2 + g2L*MGl*(-(g3L*xd3) + g3R*(xd0 + xd2 + xd3))) + 
     &          g1L*(g2L*m3*MGl*(-(g3L*xd3) + g3R*(xd0 + xd3)) + g2R*(g3R*MGl**2*xd3 + 
     &           g3L*(xc0 - xc2 - 4*xd00 + m3**2*xd1 + m1**2*xd13 + m2**2*xd13 - tman*xd13 - uman*xd13 + tman*xd2 + m1**2*xd23 + 2*m2**2*xd23 + m3**2*xd23 - sman*xd23 - 2*tman*xd23 - 
     &              uman*xd23 + 2*m1**2*xd33 + 2*m2**2*xd33 + 2*m3**2*xd33 - 2*sman*xd33 - 2*tman*xd33 - 2*uman*xd33))))     
     
          Bv0R = Bv0R - (g0L*m3*(g1L*g2R*g3L*m1*xd1 + g1R*(g2R*g3L*m3*(xd1 + xd2) + g2L*MGl*(-(g3L*xd3) + g3R*(xd0 + xd1 + xd2 + xd3))))) + 
     &          g0R*(g1L*m1*(g2L*g3R*m3*xd2 + g2R*MGl*(-(g3R*xd3) + g3L*(xd0 + xd2 + xd3))) + 
     &          g1R*(g2R*m3*MGl*(-(g3R*xd3) + g3L*(xd0 + xd3)) + g2L*(g3L*MGl**2*xd3 + 
     &           g3R*(xc0 - xc2 - 4*xd00 + m3**2*xd1 + m1**2*xd13 + m2**2*xd13 - tman*xd13 - uman*xd13 + tman*xd2 + m1**2*xd23 + 2*m2**2*xd23 + m3**2*xd23 - sman*xd23 - 2*tman*xd23 - 
     &              uman*xd23 + 2*m1**2*xd33 + 2*m2**2*xd33 + 2*m3**2*xd33 - 2*sman*xd33 - 2*tman*xd33 - 2*uman*xd33))))     

          Bv1L = Bv1L + 2*(g0L*g1L*g2R*g3L*m3*(xd12 + xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33) + 
     &          g0R*g3R*(g1R*g2L*m1*(xd2 + xd22 + 2*xd23 + xd3 + xd33) + g1L*(g2R*MGl*(xd23 + xd33) + g2L*m3*(xd2 + xd23 + xd3 + xd33))))     

          Bv1R = Bv1R + 2*(g0R*g1R*g2L*g3R*m3*(xd12 + xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33) + 
     &          g0L*g3L*(g1L*g2R*m1*(xd2 + xd22 + 2*xd23 + xd3 + xd33) + g1R*(g2L*MGl*(xd23 + xd33) + g2R*m3*(xd2 + xd23 + xd3 + xd33))))     

          Bv2L = Bv2L + 2*(g0L*g1L*g2R*g3L*m3*(xd13 + xd23 + xd3 + xd33) + g0R*g3R*(g1R*g2L*m1*(xd23 + xd3 + xd33) + g1L*(g2R*MGl*xd33 + g2L*m3*(xd3 + xd33))))
          
          Bv2R = Bv2R + 2*(g0R*g1R*g2L*g3R*m3*(xd13 + xd23 + xd3 + xd33) + g0L*g3L*(g1L*g2R*m1*(xd23 + xd3 + xd33) + g1R*(g2L*MGl*xd33 + g2R*m3*(xd3 + xd33))))
          
          Bv3L = Bv3L - 2*(g0L*g1L*g2R*g3L*m3*(xd1 + xd11 + 2*xd12 + 2*xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33) + 
     &          g0R*(g1R*g2L*g3R*m1*(xd12 + xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33) + 
     &          g1L*(g2L*g3R*m3*(xd13 + xd23 + xd3 + xd33) + g2R*MGl*(-(g3L*(xd0 + xd1 + xd2 + xd3)) + g3R*(xd13 + xd23 + xd3 + xd33)))))

          Bv3R = Bv3R - 2*(g0R*g1R*g2L*g3R*m3*(xd1 + xd11 + 2*xd12 + 2*xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33) + 
     &          g0L*(g1L*g2R*g3L*m1*(xd12 + xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33) + 
     &          g1R*(g2R*g3L*m3*(xd13 + xd23 + xd3 + xd33) + g2L*MGl*(-(g3R*(xd0 + xd1 + xd2 + xd3)) + g3L*(xd13 + xd23 + xd3 + xd33)))))     
        end do

!         print*,'------------ Karol Gluino Box -----------------------------'
!         print*,'Bs0L=',Bs0L
!         print*,'Bs0R=',Bs0R
!         print*,'Bs1L=',Bs1L
!         print*,'Bs1R=',Bs1R
!         print*,'Bs2L=',Bs2L
!         print*,'Bs2R=',Bs2R
!         print*,'Bs3L=',Bs3L
!         print*,'Bs3R=',Bs3R
!         print*,'------------'
!         print*,'Bv0L=',Bv0L
!         print*,'Bv0R=',Bv0R
!         print*,'Bv1L=',Bv1L
!         print*,'Bv1R=',Bv1R
!         print*,'Bv2L=',Bv2L
!         print*,'Bv2R=',Bv2R
!         print*,'Bv3L=',Bv3L
!         print*,'Bv3R=',Bv3R

        
! --- Interference with s channel tree-level --

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))
        
        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0
        
        boxGluinoST = CF*NC*2d0*(-CV/2d0)*MsBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A1Lc,A1Rc,BLc,BRc,mq1c,xv)


! --- Interference with t channel tree-level --

               
        A2Lc = dconjg(AL(chimeigen,ptmeigen,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen,sfmeigen,sftype,sfgen,xtype))
 
        mq2c = MSf(ptmeigen,pttype,ptgen)
c        wquark = width(pstype)
        wquark = 0d0
 
        boxGluinoST = boxGluinoST +  CF*NC*2d0*(-CV/2d0)*MtBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A2Lc,A2Rc,CLc,mq2c,xv)
     
     
         NeuQ2qx_M2BoxG = NeuQ2qx_M2BoxG + boxGluinoST
        
!         print *,"boxGluonS",boxGluonS
!         print *,"boxGluinoS",boxGluinoS
!         print *,"boxGluonT",boxGluonT
!         print *,"boxGluonT4",boxGluonT4
!         print *,"boxGluinoT",boxGluinoT
!         print *,"boxGluonST",boxGluonST
!         print *,"boxGluinoST",boxGluinoST

!           NeuQ2qx_M2BoxG =  boxGluonST

        
        end
