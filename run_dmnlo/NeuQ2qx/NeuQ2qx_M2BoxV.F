! =========================================================================
        double precision function NeuQ2qx_M2BoxV()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

        double precision MsBoxV, MtBoxV, MuBoxV
        double precision boxSSgluon, boxSSgluino, boxSTgluon, boxSTgluino, boxSUgluon, boxSUgluino
        double precision boxTSgluon, boxTSgluino, boxTTgluon, boxTTgluino, boxTUgluon, boxTUgluino
        double precision boxFVS, boxFVT, boxFVU
        
        double complex Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R
        double complex Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R

c =====================
c B.H.: Box test (to be deleted)
c		double complex bjBs0L,bjBs0R,bjBs1L,bjBs1R,bjBs2L,bjBs2R,bjBs3L,bjBs3R
c		double complex bjBv0L,bjBv0R,bjBv1L,bjBv1R,bjBv2L,bjBv2R,bjBv3L,bjBv3R
c		complex*16 xd31,xd32,xd21		
c =====================

        double complex A1Lc,A1Rc,A2Lc,A2Rc,A3Lc,A3Rc,BLc,BRc,CLc,DLc,DRc
        
        double precision mq1,mq1c,mq2,mq2c,mq3,mq3c,wquark
        
        double precision xv, colorfactor 
        integer k,l
        double precision Krondelta

        double complex g0, g0L, g0R, g1, g1L, g1R, g2, g2L, g2R, g3, g3L, g3R
        double precision Mloop0, Mloop1, Mloop2, Mloop3
        
    	integer sumUc, sumTc
        integer sumUcend
        integer sumS1, sumS2, sumT
	
        sumUc = 0
        sumTc = 0
        sumS1 = 0d0
        sumS2 = 0d0
        sumT = 0d0
                
        if (xtype.eq.2) then
        
           xv = 0d0
           colorfactor = 1d0
           
           sumUcend = -1
           
        elseif (xtype.eq.3) then
        
           xv = 1/(m4**2)
c            xv = 0d0
           colorfactor = 1d0
        
           sumUcend = 3
           
        elseif (xtype.eq.4) then

           xv = 1/(m4**2)
           colorfactor = 1d0
               
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
        g3 = 0d0
        g3L = 0d0
        g3R = 0d0
        
        Mloop0 = 0d0
        Mloop1 = 0d0
        Mloop2 = 0d0
        Mloop3 = 0d0

        mq1c = 0d0
        mq2c = 0d0
        mq3c = 0d0

	
    	boxSSgluon = 0d0
    	boxSSgluino = 0d0
    	boxSTgluon = 0d0
    	boxSTgluino = 0d0
        boxSUgluon = 0d0
        boxSUgluino = 0d0
        boxTSgluon = 0d0
        boxTSgluino = 0d0
        boxTTgluon = 0d0
        boxTTgluino = 0d0
        boxTUgluon = 0d0
        boxTUgluino = 0d0
        NeuQ2qx_M2BoxV = 0d0

! ----- s channel box with gluon -----

        g0 = CL(sfmeigen,sfmeigen,sftype,sfgen,1)
        g1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        g1R = AR(chimeigen,sfmeigen,sftype,sfgen)
        g2L = BL(xtype,ftype)   
        g2R = BR(xtype,ftype)
        g3 = B(1,ftype)
        
        Mloop1 = MSf(sfmeigen,sftype,sfgen)
        Mloop2 = Mf(sftype,sfgen)
        Mloop3 = Mf(ftype,fgen)

        call myCget(m1**2,m4**2,uman,Mloop1**2,Mloop2**2,Mloop3**2)
        call myDget(m2**2,m1**2,m4**2,m3**2,sman,uman,0d0,Mloop1**2,Mloop2**2,Mloop3**2)

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


c ==========================
c B.H.: Box test s-channel with gluon (to be deleted)

c       xd31 = xd13
c       xd32 = xd23
c       xd21 = xd12
c   
c       bjBs0L = g0*g3*(-1.*g1R*m1*(g2R*MLoop2*(m3 - 1.*MLoop3)*xd2 + 
c     &       g2L*(xc0 + xc1 + MLoop0**2*xd0 - 2.*uman*xd0 - 1.*uman*xd1 + MLoop0**2*xd2 - 1.*tman*xd2 - 3.*uman*xd2 - 
c     &          1.*uman*xd21 - 1.*tman*xd22 - 1.*uman*xd22 + m1**2*(xd2 + xd22) + m3*MLoop3*xd3 - 2.*uman*xd3 + 
c     &          m2**2*(2.*xd0 + 3.*xd1 + 3.*xd2 + xd21 + xd22 + 2.*xd3) + 
c     &          m3**2*(2.*xd0 + xd1 + 4.*xd2 + xd21 + 2.*(xd22 + xd23) + 3.*xd3))) + 
c     &    g1L*(-1.*g2L*MLoop2*(xc0 + MLoop0**2*xd0 - 2.*uman*xd0 - 1.*uman*xd1 + m1**2*xd2 - 1.*tman*xd2 - 
c     &          1.*uman*xd2 + m3*MLoop3*xd3 - 2.*uman*xd3 + m2**2*(2.*xd0 + 3.*(xd1 + xd2) + 2.*xd3) + 
c     &          m3**2*(2.*xd0 + xd1 + 2.*xd2 + 3.*xd3)) + 
c     &       g2R*(2.*m3**2*MLoop3*xd3 + MLoop3*(xc0 + MLoop0**2*xd0 - 2.*xd00 + sman*xd2 - 2.*uman*xd3 + 
c     &             m2**2*(2.*xd0 + 3.*xd1 + 2.*(xd2 + xd3))) - 
c     &          1.*m3*(xc0 - 1.*xc2 - 4.*xd00 + sman*xd2 + MLoop0**2*(xd0 - 1.*xd3) + 
c     &             m2**2*(2.*xd0 + 3.*xd1 + 2.*(xd2 + xd3) - 1.*xd31) + uman*xd31) + 
c     &          m3*(m1**2 + m2**2 - 1.*tman - 1.*uman)*xd32 + m3**3*(xd31 + 2.*(xd32 + xd33)))))
c        bjBs0R = g0*g3*(-1.*g1L*m1*(g2L*MLoop2*(m3 - 1.*MLoop3)*xd2 + 
c     &       g2R*(xc0 + xc1 + MLoop0**2*xd0 - 2.*uman*xd0 - 1.*uman*xd1 + MLoop0**2*xd2 - 1.*tman*xd2 - 3.*uman*xd2 - 
c     &          1.*uman*xd21 - 1.*tman*xd22 - 1.*uman*xd22 + m1**2*(xd2 + xd22) + m3*MLoop3*xd3 - 2.*uman*xd3 + 
c     &          m2**2*(2.*xd0 + 3.*xd1 + 3.*xd2 + xd21 + xd22 + 2.*xd3) + 
c     &          m3**2*(2.*xd0 + xd1 + 4.*xd2 + xd21 + 2.*(xd22 + xd23) + 3.*xd3))) + 
c     &    g1R*(-1.*g2R*MLoop2*(xc0 + MLoop0**2*xd0 - 2.*uman*xd0 - 1.*uman*xd1 + m1**2*xd2 - 1.*tman*xd2 - 
c     &          1.*uman*xd2 + m3*MLoop3*xd3 - 2.*uman*xd3 + m2**2*(2.*xd0 + 3.*(xd1 + xd2) + 2.*xd3) + 
c     &          m3**2*(2.*xd0 + xd1 + 2.*xd2 + 3.*xd3)) + 
c     &       g2L*(2.*m3**2*MLoop3*xd3 + MLoop3*(xc0 + MLoop0**2*xd0 - 2.*xd00 + sman*xd2 - 2.*uman*xd3 + 
c     &             m2**2*(2.*xd0 + 3.*xd1 + 2.*(xd2 + xd3))) - 
c     &          1.*m3*(xc0 - 1.*xc2 - 4.*xd00 + sman*xd2 + MLoop0**2*(xd0 - 1.*xd3) + 
c     &             m2**2*(2.*xd0 + 3.*xd1 + 2.*(xd2 + xd3) - 1.*xd31) + uman*xd31) + 
c     &          m3*(m1**2 + m2**2 - 1.*tman - 1.*uman)*xd32 + m3**3*(xd31 + 2.*(xd32 + xd33)))))
c        bjBs1L = -2.*g0*g3*(-1.*g1R*g2L*m1*(m3 - 1.*MLoop3)*(xd2 + xd22) + 
c     &    g1L*(g2L*MLoop2*(-1.*m3 + MLoop3)*xd2 + 
c     &       g2R*(m3*MLoop3*xd23 + m3**2*xd32 - 2.*uman*xd32 + 2.*m2**2*(xd12 + xd2 + xd22 + xd32))))
c        bjBs1R = -2.*g0*g3*(-1.*g1L*g2R*m1*(m3 - 1.*MLoop3)*(xd2 + xd22) + 
c     &    g1R*(g2R*MLoop2*(-1.*m3 + MLoop3)*xd2 + 
c     &       g2L*(m3*MLoop3*xd23 + m3**2*xd32 - 2.*uman*xd32 + 2.*m2**2*(xd12 + xd2 + xd22 + xd32))))
c        bjBs2L = 2.*g0*g3*(g1R*m1*(2.*g2R*MLoop2*xd2 + g2L*
c     &        (-1.*MLoop3*(2.*xd0 + xd1 + xd12 + 3.*xd2 + xd22) + m3*(2.*xd0 + xd1 + 3.*xd2 + xd21 + xd22 + 2.*xd3))) + 
c     &    g1L*(g2L*MLoop2*((m3 - 1.*MLoop3)*(2.*xd0 + xd1 + xd2) + 2.*m3*xd3) + 
c     &       g2R*(2.*xc0 + 2.*MLoop0**2*xd0 - 4.*xd00 + 2.*sman*xd2 + 
c     &          m3*(2.*m3*xd3 - 1.*MLoop3*(xd13 + xd23 + 2.*xd3)) - 1.*(m3**2 - 2.*uman)*(xd31 + xd32) - 
c     &          2.*m2**2*(xd11 + xd12 + xd2 + xd21 + xd22 + xd31 + xd32))))
c        bjBs2R = 2.*g0*g3*(g1L*m1*(2.*g2L*MLoop2*xd2 + g2R*
c     &        (-1.*MLoop3*(2.*xd0 + xd1 + xd12 + 3.*xd2 + xd22) + m3*(2.*xd0 + xd1 + 3.*xd2 + xd21 + xd22 + 2.*xd3))) + 
c     &    g1R*(g2R*MLoop2*((m3 - 1.*MLoop3)*(2.*xd0 + xd1 + xd2) + 2.*m3*xd3) + 
c     &       g2L*(2.*xc0 + 2.*MLoop0**2*xd0 - 4.*xd00 + 2.*sman*xd2 + 
c     &          m3*(2.*m3*xd3 - 1.*MLoop3*(xd13 + xd23 + 2.*xd3)) - 1.*(m3**2 - 2.*uman)*(xd31 + xd32) - 
c     &          2.*m2**2*(xd11 + xd12 + xd2 + xd21 + xd22 + xd31 + xd32))))
c        bjBs3L = 2.*g0*g3*(g1R*g2L*m1*(m3*xd23 - 1.*MLoop3*xd32) - 
c     &    1.*g1L*g2R*(xc2 + 2.*xd00 + MLoop0**2*xd3 - 2.*uman*xd3 - 1.*uman*xd31 + m1**2*xd32 - 1.*tman*xd32 - 
c     &       1.*uman*xd32 + m3*MLoop3*xd33 - 2.*uman*xd33 + 
c     &       m2**2*(2.*xd13 + 2.*xd23 + 2.*xd3 + xd31 + xd32 + 2.*xd33) + m3**2*(2.*xd3 + xd31 + 2.*xd32 + 3.*xd33)))
c        bjBs3R = 2.*g0*g3*(g1L*g2R*m1*(m3*xd23 - 1.*MLoop3*xd32) - 
c     &    1.*g1R*g2L*(xc2 + 2.*xd00 + MLoop0**2*xd3 - 2.*uman*xd3 - 1.*uman*xd31 + m1**2*xd32 - 1.*tman*xd32 - 
c     &       1.*uman*xd32 + m3*MLoop3*xd33 - 2.*uman*xd33 + 
c     &       m2**2*(2.*xd13 + 2.*xd23 + 2.*xd3 + xd31 + xd32 + 2.*xd33) + m3**2*(2.*xd3 + xd31 + 2.*xd32 + 3.*xd33)))
c        bjBv0L = g0*g3*(-1.*g1R*m1*(2.*g2R*MLoop2*xd2 + g2L*(m3 - 1.*MLoop3)*(2.*xd0 + xd1 + 2.*xd2) + 2.*g2L*m3*xd3) - 
c     &    1.*g1L*(g2L*MLoop2*((m3 - 1.*MLoop3)*(2.*xd0 + xd1 + xd2) + 2.*m3*xd3) + 
c     &       g2R*(2.*xc0 - 1.*xc2 + 3.*MLoop0**2*xd0 - 2.*uman*xd0 - 4.*xd00 + MLoop0**2*xd1 - 3.*uman*xd1 - 
c     &          1.*uman*xd11 - 1.*tman*xd12 - 1.*uman*xd12 + MLoop0**2*xd2 + 2.*sman*xd2 - 1.*tman*xd2 - 3.*uman*xd2 - 
c     &          1.*uman*xd21 - 1.*tman*xd22 - 1.*uman*xd22 + m1**2*(xd12 + xd2 + xd22) - 1.*m3*MLoop3*xd3 - 
c     &          2.*uman*xd3 + m2**2*(2.*xd0 + 5.*xd1 + xd11 + xd12 + 3.*xd2 + xd21 + xd22 + 2.*xd3) + 
c     &          m3**2*(2.*xd0 + 3.*xd1 + xd11 + 2.*xd12 + 2.*xd13 + 4.*xd2 + xd21 + 2.*(xd22 + xd23) + 5.*xd3))))
c        bjBv0R = g0*g3*(-1.*g1L*m1*(2.*g2L*MLoop2*xd2 + g2R*(m3 - 1.*MLoop3)*(2.*xd0 + xd1 + 2.*xd2) + 2.*g2R*m3*xd3) - 
c     &    1.*g1R*(g2R*MLoop2*((m3 - 1.*MLoop3)*(2.*xd0 + xd1 + xd2) + 2.*m3*xd3) + 
c     &       g2L*(2.*xc0 - 1.*xc2 + 3.*MLoop0**2*xd0 - 2.*uman*xd0 - 4.*xd00 + MLoop0**2*xd1 - 3.*uman*xd1 - 
c     &          1.*uman*xd11 - 1.*tman*xd12 - 1.*uman*xd12 + MLoop0**2*xd2 + 2.*sman*xd2 - 1.*tman*xd2 - 3.*uman*xd2 - 
c     &          1.*uman*xd21 - 1.*tman*xd22 - 1.*uman*xd22 + m1**2*(xd12 + xd2 + xd22) - 1.*m3*MLoop3*xd3 - 
c     &          2.*uman*xd3 + m2**2*(2.*xd0 + 5.*xd1 + xd11 + xd12 + 3.*xd2 + xd21 + xd22 + 2.*xd3) + 
c     &          m3**2*(2.*xd0 + 3.*xd1 + xd11 + 2.*xd12 + 2.*xd13 + 4.*xd2 + xd21 + 2.*(xd22 + xd23) + 5.*xd3))))
c        bjBv1L = 2.*g0*g3*(-2.*g1R*g2L*m1*(xd2 + xd22) + g1L*
c     &     (-2.*g2L*MLoop2*xd2 + g2R*(-1.*MLoop3*(xd2 + xd21 + xd22) + m3*(xd12 + xd2 + xd22 + 2.*xd32))))
c        bjBv1R = 2.*g0*g3*(-2.*g1L*g2R*m1*(xd2 + xd22) + g1R*
c     &     (-2.*g2R*MLoop2*xd2 + g2L*(-1.*MLoop3*(xd2 + xd21 + xd22) + m3*(xd12 + xd2 + xd22 + 2.*xd32))))
c        bjBv2L = 2.*g0*g3*(-2.*g1R*g2L*m1*(xd2 + xd21 + xd22) + 
c     &    g1L*g2R*(-1.*MLoop3*(2.*xd0 + 3.*xd1 + xd11 + xd12 + 3.*xd2 + xd21 + xd22) + 
c     &       m3*(2.*xd0 + 3.*xd1 + xd11 + xd12 + 3.*xd2 + xd21 + xd22 + 2.*(xd3 + xd31 + xd32))))
c        bjBv2R = 2.*g0*g3*(-2.*g1L*g2R*m1*(xd2 + xd21 + xd22) + 
c     &    g1R*g2L*(-1.*MLoop3*(2.*xd0 + 3.*xd1 + xd11 + xd12 + 3.*xd2 + xd21 + xd22) + 
c     &       m3*(2.*xd0 + 3.*xd1 + xd11 + xd12 + 3.*xd2 + xd21 + xd22 + 2.*(xd3 + xd31 + xd32))))
c        bjBv3L = 2.*g0*g3*(-2.*g1R*g2L*m1*xd23 + g1L*g2R*(-1.*MLoop3*(2.*xd3 + xd31 + xd32) + m3*(xd13 + xd23 + 2.*(xd3 + xd33))))
c        bjBv3R = 2.*g0*g3*(-2.*g1L*g2R*m1*xd23 + g1R*g2L*(-1.*MLoop3*(2.*xd3 + xd31 + xd32) + m3*(xd13 + xd23 + 2.*(xd3 + xd33))))
    
   	
c   	write(*,*) ' '
c   	write(*,*) 'S-GLUON'
c   	write(*,*) ' '
c
c   	write(*,*) 'dBs0L = ',(Bs0L - bjBs0L)
c   	write(*,*) 'dBs0R = ',(Bs0R - bjBs0R)
c   	write(*,*) 'dBs1L = ',(Bs1L - bjBs1L)
c   	write(*,*) 'dBs1R = ',(Bs1R - bjBs1R)
c   	write(*,*) 'dBs2L = ',(Bs2L - bjBs2L)
c   	write(*,*) 'dBs2R = ',(Bs2R - bjBs2R)
c   	write(*,*) 'dBs3L = ',(Bs3L - bjBs3L)
c   	write(*,*) 'dBs3R = ',(Bs3R - bjBs3R)
c   	                                      
c   	write(*,*) 'dBv0L = ',(Bv0L - bjBv0L)
c   	write(*,*) 'dBv0L = ',(Bv0R - bjBv0R)
c   	write(*,*) 'dBv1L = ',(Bv1L - bjBv1L)
c   	write(*,*) 'dBv1L = ',(Bv1R - bjBv1R)
c   	write(*,*) 'dBv2L = ',(Bv2L - bjBv2L)
c   	write(*,*) 'dBv2L = ',(Bv2R - bjBv2R)
c   	write(*,*) 'dBv3L = ',(Bv3L - bjBv3L)
c   	write(*,*) 'dBv3L = ',(Bv3R - bjBv3R)
c
c       write(*,*)
c
c   	write(*,*) 'dBs0L = ',(Bs0L - bjBs0L)/Bs0L
c   	write(*,*) 'dBs0R = ',(Bs0R - bjBs0R)/Bs0R
c   	write(*,*) 'dBs1L = ',(Bs1L - bjBs1L)/Bs1L
c   	write(*,*) 'dBs1R = ',(Bs1R - bjBs1R)/Bs1R
c   	write(*,*) 'dBs2L = ',(Bs2L - bjBs2L)/Bs2L
c   	write(*,*) 'dBs2R = ',(Bs2R - bjBs2R)/Bs2R
c   	write(*,*) 'dBs3L = ',(Bs3L - bjBs3L)/Bs3L
c   	write(*,*) 'dBs3R = ',(Bs3R - bjBs3R)/Bs3R
c
c   	write(*,*) 'dBv0L = ',(Bv0L - bjBv0L)/Bv0L
c   	write(*,*) 'dBv0R = ',(Bv0R - bjBv0R)/Bv0R
c   	write(*,*) 'dBv1L = ',(Bv1L - bjBv1L)/Bv1L
c   	write(*,*) 'dBv1R = ',(Bv1R - bjBv1R)/Bv1R
c   	write(*,*) 'dBv2L = ',(Bv2L - bjBv2L)/Bv2L
c   	write(*,*) 'dBv2R = ',(Bv2R - bjBv2R)/Bv2R
c   	write(*,*) 'dBv3L = ',(Bv3L - bjBv3L)/Bv3L
c   	write(*,*) 'dBv3R = ',(Bv3R - bjBv3R)/Bv3R
c
c   	write(*,*) ' '
c   	write(*,*) ' '



! --- Interference with s channel tree-level --

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))
        
        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0
        
        boxSSgluon = colorfactor*CF*NC*2*MsBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A1Lc,A1Rc,BLc,BRc,mq1c,xv)
        
        NeuQ2qx_M2BoxV = NeuQ2qx_M2BoxV + boxSSgluon
        
! --- Interference with t channel tree-level --

        do sumTc = 0, 1
                
        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
    	mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
c	    wquark = width(pstype)
         wquark = 0d0

        boxSTgluon = colorfactor*CF*NC*2*MtBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A2Lc,A2Rc,CLc,mq2c,xv)
     
        NeuQ2qx_M2BoxV = NeuQ2qx_M2BoxV + boxSTgluon

        enddo

! --- Interference with u channel tree-level --

        do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)
c        wquark = width(pstype)
        wquark = 0d0

        boxSUgluon = colorfactor*CF*NC*2*MuBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A3Lc,A3Rc,DLc,DRc,mq3c,xv)

        NeuQ2qx_M2BoxV = NeuQ2qx_M2BoxV + boxSUgluon

        enddo

! ----- s channel box with gluino -----

        do sumS1 = 1, 2
                do sumS2 = 1, 2

        g0L =  sqrt2*GS*USfC(sfmeigen,2,sftype,sfgen)
        g0R = -sqrt2*GS*USfC(sfmeigen,1,sftype,sfgen)
c        g1L = AL(chimeigen,sumS1,sftype,sfgen)
c        g1R = AR(chimeigen,sumS1,sftype,sfgen)
c        Has been corrected (now consistent with Higgs implementation):
        g1L = AR(chimeigen,sumS1,sftype,sfgen)
        g1R = AL(chimeigen,sumS1,sftype,sfgen)
        g2 = CL(sumS2,sumS1,ftype,fgen,xtype)
        g3L =  sqrt2*GS*USfC(sumS2,2,ftype,fgen)
        g3R = -sqrt2*GS*USfC(sumS2,1,ftype,fgen)


        Mloop0 = MGl 
        Mloop1 = Mf(sftype,sfgen) 
        Mloop2 = MSf(sumS1,sftype,sfgen)
        Mloop3 = MSf(sumS2,ftype,fgen)
        
        call myCget(m1**2,m4**2,uman,Mloop1**2,Mloop2**2,Mloop3**2)
        call myDget(m2**2,m1**2,m4**2,m3**2,sman,uman,Mloop0**2,Mloop1**2,Mloop2**2,Mloop3**2)
        
        Bs0L = -2.*g1L*g2*g3R*(g0R*Mloop0 + g0L*Mloop1)*xd00
        Bs0R = -2.*g1R*g2*g3L*(g0L*Mloop0 + g0R*Mloop1)*xd00
        Bs1L = g2*(g0L*(-1.*g1R*g3L*m1*Mloop0*(xd2 + 2.*xd22) + g1L*Mloop1*(g3L*Mloop0*(xd0 + 2.*xd2) - 1.*g3R*m3*(2.*xd23 + xd3))) + 
     &    g0R*(-1.*g1R*g3L*m1*Mloop1*(xd2 + 2.*xd22) + g1L*(g3L*(xc0 + 2.*xc1 + m2**2*(xd1 + 2.*xd12) + Mloop0**2*(xd0 + 2.*xd2) - 1.*(m1**2 - 1.*sman)*(xd2 + 2.*xd22)) - 
     &          1.*g3R*m3*Mloop0*(2.*xd23 + xd3))))
        Bs1R = g2*(g0R*(-1.*g1L*g3R*m1*Mloop0*(xd2 + 2.*xd22) + g1R*Mloop1*(g3R*Mloop0*(xd0 + 2.*xd2) - 1.*g3L*m3*(2.*xd23 + xd3))) + 
     &    g0L*(-1.*g1L*g3R*m1*Mloop1*(xd2 + 2.*xd22) + g1R*(g3R*(xc0 + 2.*xc1 + m2**2*(xd1 + 2.*xd12) + Mloop0**2*(xd0 + 2.*xd2) - 1.*(m1**2 - 1.*sman)*(xd2 + 2.*xd22)) - 
     &          1.*g3L*m3*Mloop0*(2.*xd23 + xd3))))
        Bs2L = g2*(g0L*(-1.*g1R*g3L*m1*Mloop0*(2.*xd12 + xd2 + 2.*xd22) + g1L*Mloop1*(g3L*Mloop0*(xd0 + 2.*(xd1 + xd2)) - 1.*g3R*m3*(2.*(xd13 + xd23) + xd3))) - 
     &    1.*g0R*(g1R*g3L*m1*Mloop1*(2.*xd12 + xd2 + 2.*xd22) + g1L*(g3L*
     &           (xc0 + 2.*xc2 - 1.*m2**2*(xd1 + 2.*(xd11 + xd12)) - 1.*Mloop0**2*(xd0 + 2.*(xd1 + xd2)) + (m1**2 - 1.*sman)*(2.*xd12 + xd2 + 2.*xd22)) + 
     &          g3R*m3*Mloop0*(2.*(xd13 + xd23) + xd3))))
        Bs2R = g2*(g0R*(-1.*g1L*g3R*m1*Mloop0*(2.*xd12 + xd2 + 2.*xd22) + g1R*Mloop1*(g3R*Mloop0*(xd0 + 2.*(xd1 + xd2)) - 1.*g3L*m3*(2.*(xd13 + xd23) + xd3))) - 
     &    1.*g0L*(g1L*g3R*m1*Mloop1*(2.*xd12 + xd2 + 2.*xd22) + g1R*(g3R*
     &           (xc0 + 2.*xc2 - 1.*m2**2*(xd1 + 2.*(xd11 + xd12)) - 1.*Mloop0**2*(xd0 + 2.*(xd1 + xd2)) + (m1**2 - 1.*sman)*(2.*xd12 + xd2 + 2.*xd22)) + 
     &          g3L*m3*Mloop0*(2.*(xd13 + xd23) + xd3))))
        Bs3L = g2*(g0L*(-1.*g1R*g3L*m1*Mloop0*(xd2 + 2.*xd23) + g1L*Mloop1*(g3L*Mloop0*(xd0 + 2.*xd3) - 1.*g3R*m3*(xd3 + 2.*xd33))) + 
     &    g0R*(-1.*g1R*g3L*m1*Mloop1*(xd2 + 2.*xd23) + g1L*(g3L*(xc0 + 2.*xc2 + m2**2*(xd1 + 2.*xd13) - 1.*(m1**2 - 1.*sman)*(xd2 + 2.*xd23) + Mloop0**2*(xd0 + 2.*xd3)) - 
     &          1.*g3R*m3*Mloop0*(xd3 + 2.*xd33))))
        Bs3R = g2*(g0R*(-1.*g1L*g3R*m1*Mloop0*(xd2 + 2.*xd23) + g1R*Mloop1*(g3R*Mloop0*(xd0 + 2.*xd3) - 1.*g3L*m3*(xd3 + 2.*xd33))) + 
     &    g0L*(-1.*g1L*g3R*m1*Mloop1*(xd2 + 2.*xd23) + g1R*(g3R*(xc0 + 2.*xc2 + m2**2*(xd1 + 2.*xd13) - 1.*(m1**2 - 1.*sman)*(xd2 + 2.*xd23) + Mloop0**2*(xd0 + 2.*xd3)) - 
     &          1.*g3L*m3*Mloop0*(xd3 + 2.*xd33))))
        Bv0L = 2.*g0R*g1L*g2*g3L*xd00
        Bv0R = 2.*g0L*g1R*g2*g3R*xd00
        Bv1L = -1.*g2*(g0L*g3R*(g1R*m1*(xd2 + 2.*xd22) + g1L*Mloop1*(xd1 + 2.*xd12 + xd2 + 2.*xd22)) + g0R*g1L*(g3R*Mloop0*(xd0 + xd1 + 2.*xd12 + 3.*xd2 + 2.*xd22) - 1.*g3L*m3*(2.*xd23 + xd3)))
        Bv1R = -1.*g2*(g0R*g3L*(g1L*m1*(xd2 + 2.*xd22) + g1R*Mloop1*(xd1 + 2.*xd12 + xd2 + 2.*xd22)) + g0L*g1R*(g3L*Mloop0*(xd0 + xd1 + 2.*xd12 + 3.*xd2 + 2.*xd22) - 1.*g3R*m3*(2.*xd23 + xd3)))
        Bv2L = -1.*g2*(g0L*g3R*(g1R*m1*(2.*xd12 + xd2 + 2.*xd22) + g1L*Mloop1*(xd1 + 2.*xd11 + 4.*xd12 + xd2 + 2.*xd22)) + 
     &    g0R*g1L*(g3R*Mloop0*(xd0 + 3.*xd1 + 2.*xd11 + 4.*xd12 + 3.*xd2 + 2.*xd22) - 1.*g3L*m3*(2.*(xd13 + xd23) + xd3)))
        Bv2R = -1.*g2*(g0R*g3L*(g1L*m1*(2.*xd12 + xd2 + 2.*xd22) + g1R*Mloop1*(xd1 + 2.*xd11 + 4.*xd12 + xd2 + 2.*xd22)) + 
     &    g0L*g1R*(g3L*Mloop0*(xd0 + 3.*xd1 + 2.*xd11 + 4.*xd12 + 3.*xd2 + 2.*xd22) - 1.*g3R*m3*(2.*(xd13 + xd23) + xd3)))
        Bv3L = -1.*g2*(g0L*g3R*(g1R*m1*(xd2 + 2.*xd23) + g1L*Mloop1*(xd1 + 2.*xd13 + xd2 + 2.*xd23)) + 
     &    g0R*g1L*(g3R*Mloop0*(xd0 + xd1 + 2.*xd13 + xd2 + 2.*(xd23 + xd3)) - 1.*g3L*m3*(xd3 + 2.*xd33)))
        Bv3R = -1.*g2*(g0R*g3L*(g1L*m1*(xd2 + 2.*xd23) + g1R*Mloop1*(xd1 + 2.*xd13 + xd2 + 2.*xd23)) + 
     &    g0L*g1R*(g3L*Mloop0*(xd0 + xd1 + 2.*xd13 + xd2 + 2.*(xd23 + xd3)) - 1.*g3R*m3*(xd3 + 2.*xd33)))

c ==========================
c B.H.: Box test s-channel with gluino (to be deleted)

c       xd31 = xd13
c       xd32 = xd23
c       xd21 = xd12
c   
c       bjBs0L = -2.*g1L*g2*g3R*(g0R*MLoop0 + g0L*MLoop1)*xd00
c       bjBs0R = -2.*g1R*g2*g3L*(g0L*MLoop0 + g0R*MLoop1)*xd00
c       bjBs1L = g2*(g0L*(-1.*g1R*g3L*m1*MLoop0*(xd2 + 2.*xd22) + 
c    &       g1L*MLoop1*(g3L*MLoop0*(xd0 + 2.*xd2) - 1.*g3R*m3*(2.*xd23 + xd3))) + 
c    &    g0R*(-1.*g1R*g3L*m1*MLoop1*(xd2 + 2.*xd22) + 
c    &       g1L*(g3L*(xc0 + 2.*xc1 + m2**2*(xd1 + 2.*xd12) + MLoop0**2*(xd0 + 2.*xd2) - 
c    &             1.*(m1**2 - 1.*sman)*(xd2 + 2.*xd22)) - 1.*g3R*m3*MLoop0*(2.*xd23 + xd3))))
c       bjBs1R = g2*(g0R*(-1.*g1L*g3R*m1*MLoop0*(xd2 + 2.*xd22) + 
c    &       g1R*MLoop1*(g3R*MLoop0*(xd0 + 2.*xd2) - 1.*g3L*m3*(2.*xd23 + xd3))) + 
c    &    g0L*(-1.*g1L*g3R*m1*MLoop1*(xd2 + 2.*xd22) + 
c    &       g1R*(g3R*(xc0 + 2.*xc1 + m2**2*(xd1 + 2.*xd12) + MLoop0**2*(xd0 + 2.*xd2) - 
c    &             1.*(m1**2 - 1.*sman)*(xd2 + 2.*xd22)) - 1.*g3L*m3*MLoop0*(2.*xd23 + xd3))))
c       bjBs2L = g2*(g0L*(-1.*g1R*g3L*m1*MLoop0*(2.*xd12 + xd2 + 2.*xd22) + 
c    &       g1L*MLoop1*(g3L*MLoop0*(xd0 + 2.*(xd1 + xd2)) - 1.*g3R*m3*(2.*(xd13 + xd23) + xd3))) - 
c    &    1.*g0R*(g1R*g3L*m1*MLoop1*(2.*xd12 + xd2 + 2.*xd22) + 
c    &       g1L*(g3L*(xc0 + 2.*xc2 - 1.*m2**2*(xd1 + 2.*(xd11 + xd12)) - 1.*MLoop0**2*(xd0 + 2.*(xd1 + xd2)) + 
c    &             (m1**2 - 1.*sman)*(xd2 + 2.*(xd21 + xd22))) + g3R*m3*MLoop0*(2.*(xd13 + xd23) + xd3))))
c       bjBs2R = g2*(g0R*(-1.*g1L*g3R*m1*MLoop0*(2.*xd12 + xd2 + 2.*xd22) + 
c    &       g1R*MLoop1*(g3R*MLoop0*(xd0 + 2.*(xd1 + xd2)) - 1.*g3L*m3*(2.*(xd13 + xd23) + xd3))) - 
c    &    1.*g0L*(g1L*g3R*m1*MLoop1*(2.*xd12 + xd2 + 2.*xd22) + 
c    &       g1R*(g3R*(xc0 + 2.*xc2 - 1.*m2**2*(xd1 + 2.*(xd11 + xd12)) - 1.*MLoop0**2*(xd0 + 2.*(xd1 + xd2)) + 
c    &             (m1**2 - 1.*sman)*(xd2 + 2.*(xd21 + xd22))) + g3L*m3*MLoop0*(2.*(xd13 + xd23) + xd3))))
c       bjBs3L = g2*(g0L*(-1.*g1R*g3L*m1*MLoop0*(xd2 + 2.*xd32) + 
c    &       g1L*MLoop1*(g3L*MLoop0*(xd0 + 2.*xd3) - 1.*g3R*m3*(xd3 + 2.*xd33))) + 
c    &    g0R*(-1.*g1R*g3L*m1*MLoop1*(xd2 + 2.*xd32) + 
c    &       g1L*(g3L*(xc0 + 2.*xc2 + m2**2*(xd1 + 2.*xd13) - 1.*(m1**2 - 1.*sman)*(xd2 + 2.*xd23) + 
c    &             MLoop0**2*(xd0 + 2.*xd3)) - 1.*g3R*m3*MLoop0*(xd3 + 2.*xd33))))
c       bjBs3R = g2*(g0R*(-1.*g1L*g3R*m1*MLoop0*(xd2 + 2.*xd32) + 
c    &       g1R*MLoop1*(g3R*MLoop0*(xd0 + 2.*xd3) - 1.*g3L*m3*(xd3 + 2.*xd33))) + 
c    &    g0L*(-1.*g1L*g3R*m1*MLoop1*(xd2 + 2.*xd32) + 
c    &       g1R*(g3R*(xc0 + 2.*xc2 + m2**2*(xd1 + 2.*xd13) - 1.*(m1**2 - 1.*sman)*(xd2 + 2.*xd23) + 
c    &             MLoop0**2*(xd0 + 2.*xd3)) - 1.*g3L*m3*MLoop0*(xd3 + 2.*xd33))))
c       bjBv0L = 2.*g0R*g1L*g2*g3L*xd00
c       bjBv0R = 2.*g0L*g1R*g2*g3R*xd00
c       bjBv1L = g2*(-1.*g0L*g3R*(g1R*m1*(xd2 + 2.*xd22) + g1L*MLoop1*(xd1 + xd2 + 2.*(xd21 + xd22))) - 
c    &    1.*g0R*g1L*(g3R*MLoop0*(xd0 + xd1 + 3.*xd2 + 2.*(xd21 + xd22)) - 1.*g3L*m3*(xd3 + 2.*xd32)))
c       bjBv1R = g2*(-1.*g0R*g3L*(g1L*m1*(xd2 + 2.*xd22) + g1R*MLoop1*(xd1 + xd2 + 2.*(xd21 + xd22))) - 
c    &    1.*g0L*g1R*(g3L*MLoop0*(xd0 + xd1 + 3.*xd2 + 2.*(xd21 + xd22)) - 1.*g3R*m3*(xd3 + 2.*xd32)))
c       bjBv2L = g2*(-1.*g0L*g3R*(g1R*m1*(xd2 + 2.*(xd21 + xd22)) + 
c    &       g1L*MLoop1*(xd1 + 2.*xd11 + 2.*xd12 + xd2 + 2.*(xd21 + xd22))) - 
c    &    1.*g0R*g1L*(g3R*MLoop0*(xd0 + 3.*xd1 + 2.*xd11 + 2.*xd12 + 3.*xd2 + 2.*(xd21 + xd22)) - 
c    &       1.*g3L*m3*(xd3 + 2.*(xd31 + xd32))))
c       bjBv2R = g2*(-1.*g0R*g3L*(g1L*m1*(xd2 + 2.*(xd21 + xd22)) + 
c    &       g1R*MLoop1*(xd1 + 2.*xd11 + 2.*xd12 + xd2 + 2.*(xd21 + xd22))) - 
c    &    1.*g0L*g1R*(g3L*MLoop0*(xd0 + 3.*xd1 + 2.*xd11 + 2.*xd12 + 3.*xd2 + 2.*(xd21 + xd22)) - 
c    &       1.*g3R*m3*(xd3 + 2.*(xd31 + xd32))))
c       bjBv3L = g2*(-1.*g0L*g3R*(g1R*m1*(xd2 + 2.*xd23) + g1L*MLoop1*(xd1 + xd2 + 2.*(xd31 + xd32))) - 
c    &    1.*g0R*g1L*(g3R*MLoop0*(xd0 + xd1 + xd2 + 2.*(xd3 + xd31 + xd32)) - 1.*g3L*m3*(xd3 + 2.*xd33)))
c       bjBv3R = g2*(-1.*g0R*g3L*(g1L*m1*(xd2 + 2.*xd23) + g1R*MLoop1*(xd1 + xd2 + 2.*(xd31 + xd32))) - 
c    &    1.*g0L*g1R*(g3L*MLoop0*(xd0 + xd1 + xd2 + 2.*(xd3 + xd31 + xd32)) - 1.*g3R*m3*(xd3 + 2.*xd33)))
c    
c   	
c   	write(*,*) ' '
c   	write(*,*) 'S-GLUINO'
c   	write(*,*) ' '
c   	write(*,*) 'dBs0L = ',(Bs0L - bjBs0L)
c   	write(*,*) 'dBs0R = ',(Bs0R - bjBs0R)
c   	write(*,*) 'dBs1L = ',(Bs1L - bjBs1L)
c   	write(*,*) 'dBs1R = ',(Bs1R - bjBs1R)
c   	write(*,*) 'dBs2L = ',(Bs2L - bjBs2L)
c   	write(*,*) 'dBs2R = ',(Bs2R - bjBs2R)
c   	write(*,*) 'dBs3L = ',(Bs3L - bjBs3L)
c   	write(*,*) 'dBs3R = ',(Bs3R - bjBs3R)
c   	                                      
c   	write(*,*) 'dBv0L = ',(Bv0L - bjBv0L)
c   	write(*,*) 'dBv0L = ',(Bv0R - bjBv0R)
c   	write(*,*) 'dBv1L = ',(Bv1L - bjBv1L)
c   	write(*,*) 'dBv1L = ',(Bv1R - bjBv1R)
c   	write(*,*) 'dBv2L = ',(Bv2L - bjBv2L)
c   	write(*,*) 'dBv2L = ',(Bv2R - bjBv2R)
c   	write(*,*) 'dBv3L = ',(Bv3L - bjBv3L)
c   	write(*,*) 'dBv3L = ',(Bv3R - bjBv3R)
c
c       write(*,*)
c
c   	write(*,*) 'dBs0L = ',(Bs0L - bjBs0L)/Bs0L
c   	write(*,*) 'dBs0R = ',(Bs0R - bjBs0R)/Bs0R
c   	write(*,*) 'dBs1L = ',(Bs1L - bjBs1L)/Bs1L
c   	write(*,*) 'dBs1R = ',(Bs1R - bjBs1R)/Bs1R
c   	write(*,*) 'dBs2L = ',(Bs2L - bjBs2L)/Bs2L
c   	write(*,*) 'dBs2R = ',(Bs2R - bjBs2R)/Bs2R
c   	write(*,*) 'dBs3L = ',(Bs3L - bjBs3L)/Bs3L
c   	write(*,*) 'dBs3R = ',(Bs3R - bjBs3R)/Bs3R
c
c   	write(*,*) 'dBv0L = ',(Bv0L - bjBv0L)/Bv0L
c   	write(*,*) 'dBv0R = ',(Bv0R - bjBv0R)/Bv0R
c   	write(*,*) 'dBv1L = ',(Bv1L - bjBv1L)/Bv1L
c   	write(*,*) 'dBv1R = ',(Bv1R - bjBv1R)/Bv1R
c   	write(*,*) 'dBv2L = ',(Bv2L - bjBv2L)/Bv2L
c   	write(*,*) 'dBv2R = ',(Bv2R - bjBv2R)/Bv2R
c   	write(*,*) 'dBv3L = ',(Bv3L - bjBv3L)/Bv3L
c   	write(*,*) 'dBv3R = ',(Bv3R - bjBv3R)/Bv3R
c
c   	
c   	write(*,*) ' '
c   	write(*,*) ' '

c ==========================

! --- Interference with s channel tree-level --

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))
        
        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0
        
        boxSSgluino = colorfactor*CF*NC*2*MsBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A1Lc,A1Rc,BLc,BRc,mq1c,xv)
        
        NeuQ2qx_M2BoxV = NeuQ2qx_M2BoxV + boxSSgluino
        
! --- Interference with t channel tree-level --

                        do sumTc = 0, 1
                
        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
    	mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
c	    wquark = width(pstype)
         wquark = 0d0

        boxSTgluino = colorfactor*CF*NC*2*MtBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A2Lc,A2Rc,CLc,mq2c,xv)
   
        NeuQ2qx_M2BoxV = NeuQ2qx_M2BoxV + boxSTgluino

                        enddo


! --- Interference with u channel tree-level --

                        do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)
c        wquark = width(pstype)
        wquark = 0d0

        boxSUgluino = colorfactor*CF*NC*2* MuBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A3Lc,A3Rc,DLc,DRc,mq3c,xv)

        NeuQ2qx_M2BoxV = NeuQ2qx_M2BoxV + boxSUgluino


                        enddo

                enddo
        enddo


! ----- t channel box with gluon -----

        do sumT = 1, 2

        g0 = CL(sfmeigen,sfmeigen,sftype,sfgen,1)
        g1 = B(1,ftype)
        g2L = AL(chimeigen,sumT,ftype,fgen)
        g2R = AR(chimeigen,sumT,ftype,fgen)
        g3 = CL(sumT,sfmeigen,sftype,sfgen,xtype)      
        
        Mloop0 = 0d0
        Mloop1 = Mf(ftype,fgen)
        Mloop2 = MSf(sumT,ftype,fgen)
        Mloop3 = MSf(sfmeigen,sftype,sfgen)

        call myCget(m1**2,m4**2,uman,Mloop1**2,Mloop2**2,Mloop3**2)
        call myDget(m3**2,m1**2,m4**2,m2**2,tman,uman,0d0,Mloop1**2,Mloop2**2,Mloop3**2)

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
      
c ==========================
c B.H.: Box test t-channel with gluon (to be deleted)

c       xd31 = xd13
c       xd32 = xd23
c       xd21 = xd12
c   
c       bjBs0L = 2.*g0*g1*g2L*g3*(m3 - 1.*MLoop1)*xd00
c       bjBs0R = 2.*g0*g1*g2R*g3*(m3 - 1.*MLoop1)*xd00
c       bjBs1L = g0*g1*g3*(g2R*m1*(m3 - 1.*MLoop1)*(xd2 + 2.*xd22) + 
c    &    g2L*(xc0 + 2.*xc1 + MLoop0**2*xd0 - 2.*uman*xd0 - 2.*uman*xd1 - 4.*uman*xd12 + 2.*MLoop0**2*xd2 + tman*xd2 - 
c    &       6.*uman*xd2 + 2.*tman*xd22 - 4.*uman*xd22 - 1.*m1**2*(xd2 + 2.*xd22) + 
c    &       m3*MLoop1*(xd1 + xd2 + 2.*(xd21 + xd22)) - 1.*uman*xd3 - 2.*uman*xd32 + 
c    &       m3**2*(2.*xd0 + 3.*xd1 + 6.*xd12 + 6.*xd2 + 4.*xd22 + xd3 + 2.*xd32) + 
c    &       m2**2*(2.*xd0 + 2.*xd1 + 4.*xd12 + 6.*xd2 + 4.*xd22 + 3.*xd3 + 6.*xd32)))
c       bjBs1R = g0*g1*g3*(g2L*m1*(m3 - 1.*MLoop1)*(xd2 + 2.*xd22) + 
c    &    g2R*(xc0 + 2.*xc1 + MLoop0**2*xd0 - 2.*uman*xd0 - 2.*uman*xd1 - 4.*uman*xd12 + 2.*MLoop0**2*xd2 + tman*xd2 - 
c    &       6.*uman*xd2 + 2.*tman*xd22 - 4.*uman*xd22 - 1.*m1**2*(xd2 + 2.*xd22) + 
c    &       m3*MLoop1*(xd1 + xd2 + 2.*(xd21 + xd22)) - 1.*uman*xd3 - 2.*uman*xd32 + 
c    &       m3**2*(2.*xd0 + 3.*xd1 + 6.*xd12 + 6.*xd2 + 4.*xd22 + xd3 + 2.*xd32) + 
c    &       m2**2*(2.*xd0 + 2.*xd1 + 4.*xd12 + 6.*xd2 + 4.*xd22 + 3.*xd3 + 6.*xd32)))
c       bjBs2L = g0*g1*g3*(g2R*m1*(-1.*m3*(xd2 + 2.*xd23) + MLoop1*(xd2 + 2.*xd32)) - 
c    &    1.*g2L*(xc0 + 2.*xc2 + MLoop0**2*xd0 - 2.*uman*xd0 + 8.*xd00 - 2.*uman*xd1 - 4.*uman*xd13 - 2.*uman*xd2 - 
c    &       4.*uman*xd23 - 1.*(m1**2 - 1.*tman)*(xd2 + 2.*xd23) + 2.*MLoop0**2*xd3 - 5.*uman*xd3 + 
c    &       m3*MLoop1*(xd1 + xd2 + 2.*(xd31 + xd32)) - 2.*uman*xd33 + 
c    &       m3**2*(2.*xd0 + 3.*xd1 + 6.*xd13 + 2.*xd2 + 4.*xd23 + 5.*xd3 + 2.*xd33) + 
c    &       m2**2*(2.*xd0 + 2.*xd1 + 4.*xd13 + 2.*xd2 + 4.*xd23 + 7.*xd3 + 6.*xd33)))
c       bjBs2R = g0*g1*g3*(g2L*m1*(-1.*m3*(xd2 + 2.*xd23) + MLoop1*(xd2 + 2.*xd32)) - 
c    &    1.*g2R*(xc0 + 2.*xc2 + MLoop0**2*xd0 - 2.*uman*xd0 + 8.*xd00 - 2.*uman*xd1 - 4.*uman*xd13 - 2.*uman*xd2 - 
c    &       4.*uman*xd23 - 1.*(m1**2 - 1.*tman)*(xd2 + 2.*xd23) + 2.*MLoop0**2*xd3 - 5.*uman*xd3 + 
c    &       m3*MLoop1*(xd1 + xd2 + 2.*(xd31 + xd32)) - 2.*uman*xd33 + 
c    &       m3**2*(2.*xd0 + 3.*xd1 + 6.*xd13 + 2.*xd2 + 4.*xd23 + 5.*xd3 + 2.*xd33) + 
c    &       m2**2*(2.*xd0 + 2.*xd1 + 4.*xd13 + 2.*xd2 + 4.*xd23 + 7.*xd3 + 6.*xd33)))
c       bjBs3L = g0*g1*g3*(g2R*m1*(MLoop1*(2.*xd12 + xd2 + 2.*xd22) - 1.*m3*(xd2 + 2.*(xd21 + xd22))) + 
c    &    g2L*(xc0 + 2.*xc2 - 1.*MLoop0**2*xd0 + 2.*uman*xd0 - 4.*xd00 - 2.*MLoop0**2*xd1 + 6.*uman*xd1 + 
c    &       4.*uman*xd11 + 4.*uman*xd12 - 2.*MLoop0**2*xd2 - 1.*tman*xd2 + 6.*uman*xd2 - 2.*tman*xd21 + 4.*uman*xd21 - 
c    &       2.*tman*xd22 + 4.*uman*xd22 + m1**2*(xd2 + 2.*(xd21 + xd22)) - 
c    &       1.*m3*MLoop1*(xd1 + 2.*xd11 + 2.*xd12 + xd2 + 2.*(xd21 + xd22)) + uman*xd3 + 2.*uman*xd31 + 2.*uman*xd32 - 
c    &       1.*m2**2*(2.*xd0 + 6.*xd1 + 4.*xd11 + 4.*xd12 + 6.*xd2 + 4.*xd21 + 4.*xd22 + 3.*xd3 + 6.*xd31 + 6.*xd32) - 
c    &       1.*m3**2*(2.*xd0 + 7.*xd1 + 6.*xd11 + 6.*xd12 + 6.*xd2 + 4.*xd21 + 4.*xd22 + xd3 + 2.*(xd31 + xd32))))
c       bjBs3R = g0*g1*g3*(g2L*m1*(MLoop1*(2.*xd12 + xd2 + 2.*xd22) - 1.*m3*(xd2 + 2.*(xd21 + xd22))) + 
c    &    g2R*(xc0 + 2.*xc2 - 1.*MLoop0**2*xd0 + 2.*uman*xd0 - 4.*xd00 - 2.*MLoop0**2*xd1 + 6.*uman*xd1 + 
c    &       4.*uman*xd11 + 4.*uman*xd12 - 2.*MLoop0**2*xd2 - 1.*tman*xd2 + 6.*uman*xd2 - 2.*tman*xd21 + 4.*uman*xd21 - 
c    &       2.*tman*xd22 + 4.*uman*xd22 + m1**2*(xd2 + 2.*(xd21 + xd22)) - 
c    &       1.*m3*MLoop1*(xd1 + 2.*xd11 + 2.*xd12 + xd2 + 2.*(xd21 + xd22)) + uman*xd3 + 2.*uman*xd31 + 2.*uman*xd32 - 
c    &       1.*m2**2*(2.*xd0 + 6.*xd1 + 4.*xd11 + 4.*xd12 + 6.*xd2 + 4.*xd21 + 4.*xd22 + 3.*xd3 + 6.*xd31 + 6.*xd32) - 
c    &       1.*m3**2*(2.*xd0 + 7.*xd1 + 6.*xd11 + 6.*xd12 + 6.*xd2 + 4.*xd21 + 4.*xd22 + xd3 + 2.*(xd31 + xd32))))
c       bjBv0L = 4.*g0*g1*g2L*g3*xd00
c       bjBv0R = 4.*g0*g1*g2R*g3*xd00
c       bjBv1L = -1.*g0*g1*g3*(2.*g2R*m1*(xd2 + 2.*xd22) + 
c    &    g2L*(-1.*MLoop1*(2.*(xd0 + 2.*xd2 + xd23) + xd3) + 
c    &       m3*(2.*xd0 + 2.*xd1 + 4.*xd12 + 6.*xd2 + 4.*xd22 + xd3 + 2.*xd32)))
c       bjBv1R = -1.*g0*g1*g3*(2.*g2L*m1*(xd2 + 2.*xd22) + 
c    &    g2R*(-1.*MLoop1*(2.*(xd0 + 2.*xd2 + xd23) + xd3) + 
c    &       m3*(2.*xd0 + 2.*xd1 + 4.*xd12 + 6.*xd2 + 4.*xd22 + xd3 + 2.*xd32)))
c       bjBv2L = g0*g1*g3*(2.*g2R*m1*(xd2 + 2.*xd23) + g2L*
c    &     (-1.*MLoop1*(2.*xd0 + 5.*xd3 + 2.*xd33) + 
c    &       m3*(2.*xd0 + 2.*xd1 + 4.*xd13 + 2.*xd2 + 4.*xd23 + 5.*xd3 + 2.*xd33)))
c       bjBv2R = g0*g1*g3*(2.*g2L*m1*(xd2 + 2.*xd23) + g2R*
c    &     (-1.*MLoop1*(2.*xd0 + 5.*xd3 + 2.*xd33) + 
c    &       m3*(2.*xd0 + 2.*xd1 + 4.*xd13 + 2.*xd2 + 4.*xd23 + 5.*xd3 + 2.*xd33)))
c       bjBv3L = g0*g1*g3*(2.*g2R*m1*(xd2 + 2.*(xd21 + xd22)) + 
c    &    g2L*(-2.*MLoop1*(xd0 + 2.*xd1 + xd13 + 2.*xd2 + xd23) - 1.*MLoop1*xd3 + 
c    &       m3*(2.*xd0 + 6.*xd1 + 4.*xd11 + 4.*xd12 + 6.*xd2 + 4.*xd21 + 4.*xd22 + xd3 + 2.*(xd31 + xd32))))
c       bjBv3R = g0*g1*g3*(2.*g2L*m1*(xd2 + 2.*(xd21 + xd22)) + 
c    &    g2R*(-2.*MLoop1*(xd0 + 2.*xd1 + xd13 + 2.*xd2 + xd23) - 1.*MLoop1*xd3 + 
c    &       m3*(2.*xd0 + 6.*xd1 + 4.*xd11 + 4.*xd12 + 6.*xd2 + 4.*xd21 + 4.*xd22 + xd3 + 2.*(xd31 + xd32))))
c    
c   	
c   	write(*,*) ' '
c   	write(*,*) 'T-GLUON'
c   	write(*,*) ' m3, Mloop1 = ',m3,Mloop1
c   	write(*,*) ' '
c   	write(*,*) 'dBs0L = ',(Bs0L - bjBs0L),Bs0L,bjBs0L
c   	write(*,*) 'dBs0R = ',(Bs0R - bjBs0R),Bs0R,bjBs0R
c   	write(*,*) 'dBs1L = ',(Bs1L - bjBs1L),Bs1L,bjBs1L
c   	write(*,*) 'dBs1R = ',(Bs1R - bjBs1R),Bs1R,bjBs1R
c   	write(*,*) 'dBs2L = ',(Bs2L - bjBs2L),Bs2L,bjBs2L
c   	write(*,*) 'dBs2R = ',(Bs2R - bjBs2R),Bs2R,bjBs2R
c   	write(*,*) 'dBs3L = ',(Bs3L - bjBs3L),Bs3L,bjBs3L
c   	write(*,*) 'dBs3R = ',(Bs3R - bjBs3R),Bs3R,bjBs3R
c   	                                      
c   	write(*,*) 'dBv0L = ',(Bv0L - bjBv0L),Bv0L,bjBv0L
c   	write(*,*) 'dBv0L = ',(Bv0R - bjBv0R),Bv0R,bjBv0R
c   	write(*,*) 'dBv1L = ',(Bv1L - bjBv1L),Bv1L,bjBv1L
c   	write(*,*) 'dBv1L = ',(Bv1R - bjBv1R),Bv1R,bjBv1R
c   	write(*,*) 'dBv2L = ',(Bv2L - bjBv2L),Bv2L,bjBv2L
c   	write(*,*) 'dBv2L = ',(Bv2R - bjBv2R),Bv2R,bjBv2R
c   	write(*,*) 'dBv3L = ',(Bv3L - bjBv3L),Bv3L,bjBv3L
c   	write(*,*) 'dBv3L = ',(Bv3R - bjBv3R),Bv3R,bjBv3R
c
c       write(*,*)
c
c   	write(*,*) 'rel. dBs0L = ',(Bs0L - bjBs0L)/Bs0L
c   	write(*,*) 'rel. dBs0R = ',(Bs0R - bjBs0R)/Bs0R
c   	write(*,*) 'rel. dBs1L = ',(Bs1L - bjBs1L)/Bs1L
c   	write(*,*) 'rel. dBs1R = ',(Bs1R - bjBs1R)/Bs1R
c   	write(*,*) 'rel. dBs2L = ',(Bs2L - bjBs2L)/Bs2L
c   	write(*,*) 'rel. dBs2R = ',(Bs2R - bjBs2R)/Bs2R
c   	write(*,*) 'rel. dBs3L = ',(Bs3L - bjBs3L)/Bs3L
c   	write(*,*) 'rel. dBs3R = ',(Bs3R - bjBs3R)/Bs3R
c                     
c   	write(*,*) 'rel. dBv0L = ',(Bv0L - bjBv0L)/Bv0L
c   	write(*,*) 'rel. dBv0R = ',(Bv0R - bjBv0R)/Bv0R
c   	write(*,*) 'rel. dBv1L = ',(Bv1L - bjBv1L)/Bv1L
c   	write(*,*) 'rel. dBv1R = ',(Bv1R - bjBv1R)/Bv1R
c   	write(*,*) 'rel. dBv2L = ',(Bv2L - bjBv2L)/Bv2L
c   	write(*,*) 'rel. dBv2R = ',(Bv2R - bjBv2R)/Bv2R
c   	write(*,*) 'rel. dBv3L = ',(Bv3L - bjBv3L)/Bv3L
c   	write(*,*) 'rel. dBv3R = ',(Bv3R - bjBv3R)/Bv3R
c
c   	write(*,*) ' '
c   	write(*,*) ' '

c ==========================

! --- Interference with s channel tree-level --

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))
        
        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0
        
        boxTSgluon = colorfactor*CF*NC*2*MsBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A1Lc,A1Rc,BLc,BRc,mq1c,xv)

        NeuQ2qx_M2BoxV = NeuQ2qx_M2BoxV + boxTSgluon


! --- Interference with t channel tree-level --

                do sumTc = 0, 1
                
        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
    	mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
c	    wquark = width(pstype)
        wquark = 0d0
 
        boxTTgluon = colorfactor*CF*NC*2*MtBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A2Lc,A2Rc,CLc,mq2c,xv)
     
        NeuQ2qx_M2BoxV = NeuQ2qx_M2BoxV + boxTTgluon

                enddo
        

! --- Interference with u channel tree-level --

                do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)
c        wquark = width(pstype)
        wquark = 0d0

        boxTUgluon = colorfactor*CF*NC*2*MuBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A3Lc,A3Rc,DLc,DRc,mq3c,xv)

        NeuQ2qx_M2BoxV = NeuQ2qx_M2BoxV + boxTUgluon


                enddo

        enddo


! ----- t channel box with gluino -----

        do sumT = 1, 2

        g0L =  sqrt2*GS*USfC(sfmeigen,2,sftype,sfgen)
        g0R = -sqrt2*GS*USfC(sfmeigen,1,sftype,sfgen)
        g1L =  sqrt2*GS*USfC(sumT,2,ftype,fgen)
        g1R = -sqrt2*GS*USfC(sumT,1,ftype,fgen)
c        g2L = AL(chimeigen,sumT,ftype,fgen)
c        g2R = AR(chimeigen,sumT,ftype,fgen)
c        Has been corrected (now consistent with Higgs implementation):
        g2L = AR(chimeigen,sumT,ftype,fgen)
        g2R = AL(chimeigen,sumT,ftype,fgen)
        g3L = BR(xtype,ftype)
        g3R = BL(xtype,ftype)



        Mloop0 = MGl
        Mloop1 = MSf(sumT,ftype,fgen)
        Mloop2 = Mf(ftype,fgen)
        Mloop3 = Mf(sftype,sfgen)
       
        call myCget(m1**2,m4**2,uman,Mloop1**2,Mloop2**2,Mloop3**2)
        call myDget(m3**2,m1**2,m4**2,m2**2,tman,uman,Mloop0**2,Mloop1**2,Mloop2**2,Mloop3**2)

        Bs0L = g0L*(g1L*m3*Mloop0*(-1.*g2R*g3L*m1*xd1 + g2L*(-1.*g3L*Mloop2*(xd1 + xd2) + g3R*Mloop3*(xd0 + xd1 + xd2))) + 
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
        Bs0R = g0R*(g1R*m3*Mloop0*(-1.*g2L*g3R*m1*xd1 + g2R*(-1.*g3R*Mloop2*(xd1 + xd2) + g3L*Mloop3*(xd0 + xd1 + xd2))) + 
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
        Bs1L = -2.*(g0L*(g1R*g2L*g3R*m3*Mloop3*(xd12 + xd2 + xd22) - 1.*g1L*g3L*Mloop0*(g2L*Mloop2*xd2 + g2R*m1*(xd2 + xd22))) + 
     &    g0R*(g1R*g2L*g3R*m3*Mloop0*(xd12 + xd2 + xd22) + g1L*(-1.*g2R*g3L*m1*Mloop3*(xd2 + xd22) + 
     &          g2L*(-1.*g3L*Mloop2*Mloop3*xd2 + g3R*(m2**2 + m3**2 - 1.*uman)*(xd12 + xd2 + xd22) + g3R*m2**2*xd23))))
        Bs1R = -2.*(g0R*(g1L*g2R*g3L*m3*Mloop3*(xd12 + xd2 + xd22) - 1.*g1R*g3R*Mloop0*(g2R*Mloop2*xd2 + g2L*m1*(xd2 + xd22))) + 
     &    g0L*(g1L*g2R*g3L*m3*Mloop0*(xd12 + xd2 + xd22) + g1R*(-1.*g2L*g3R*m1*Mloop3*(xd2 + xd22) + 
     &          g2R*(-1.*g3R*Mloop2*Mloop3*xd2 + g3L*(m2**2 + m3**2 - 1.*uman)*(xd12 + xd2 + xd22) + g3L*m2**2*xd23))))
        Bs2L = -2.*(g0L*(g1L*g3L*Mloop0*(g2L*Mloop2*(xd0 + xd3) + g2R*m1*(xd0 + xd2 + xd23 + xd3)) - 
     &       1.*g1R*m3*(g2R*g3L*m1*xd1 + g2L*(g3L*Mloop2*(xd1 + xd2) + g3R*Mloop3*(xd13 + xd23 + xd3)))) + 
     &    g0R*(-1.*g1R*g2L*g3R*m3*Mloop0*(xd0 + xd1 + xd13 + xd2 + xd23 + xd3) + 
     &       g1L*(g2R*m1*(g3R*Mloop2*xd2 + g3L*Mloop3*(xd23 + xd3)) + g2L*
     &           (g3L*Mloop2*Mloop3*xd3 + g3R*(xc0 + Mloop0**2*xd0 - 2.*xd00 + tman*xd2 + m3**2*(xd1 - 1.*xd13 - 1.*xd23) + uman*(xd13 + xd23) - 1.*m2**2*(xd13 + xd23 + xd33))))))
        Bs2R = -2.*(g0R*(g1R*g3R*Mloop0*(g2R*Mloop2*(xd0 + xd3) + g2L*m1*(xd0 + xd2 + xd23 + xd3)) - 
     &       1.*g1L*m3*(g2L*g3R*m1*xd1 + g2R*(g3R*Mloop2*(xd1 + xd2) + g3L*Mloop3*(xd13 + xd23 + xd3)))) + 
     &    g0L*(-1.*g1L*g2R*g3L*m3*Mloop0*(xd0 + xd1 + xd13 + xd2 + xd23 + xd3) + 
     &       g1R*(g2L*m1*(g3L*Mloop2*xd2 + g3R*Mloop3*(xd23 + xd3)) + g2R*
     &           (g3R*Mloop2*Mloop3*xd3 + g3L*(xc0 + Mloop0**2*xd0 - 2.*xd00 + tman*xd2 + m3**2*(xd1 - 1.*xd13 - 1.*xd23) + uman*(xd13 + xd23) - 1.*m2**2*(xd13 + xd23 + xd33))))))
        Bs3L = -2.*(g0L*(-1.*g1R*g2L*g3R*m3*Mloop3*(xd1 + xd11 + 2.*xd12 + xd2 + xd22) + g1L*Mloop0*(g2L*g3R*Mloop3*(xd0 + xd1 + xd2) + g2R*g3L*m1*(xd12 + xd2 + xd22))) + 
     &    g0R*(-1.*g1R*g2L*g3R*m3*Mloop0*(xd1 + xd11 + 2.*xd12 + xd2 + xd22) + 
     &       g1L*(g2R*g3L*m1*Mloop3*(xd12 + xd2 + xd22) + g2L*g3R*(-1.*xc2 - 1.*sman*xd12 + Mloop0**2*(xd0 + xd1 + xd2) - 1.*sman*(xd2 + xd22) + m1**2*(xd12 + xd2 + xd22) + 
     &             m2**2*(xd12 + xd13 + xd2 + xd22 + xd23 + xd3)))))
        Bs3R = -2.*(g0R*(-1.*g1L*g2R*g3L*m3*Mloop3*(xd1 + xd11 + 2.*xd12 + xd2 + xd22) + g1R*Mloop0*(g2R*g3L*Mloop3*(xd0 + xd1 + xd2) + g2L*g3R*m1*(xd12 + xd2 + xd22))) + 
     &    g0L*(-1.*g1L*g2R*g3L*m3*Mloop0*(xd1 + xd11 + 2.*xd12 + xd2 + xd22) + 
     &       g1R*(g2L*g3R*m1*Mloop3*(xd12 + xd2 + xd22) + g2R*g3L*(-1.*xc2 - 1.*sman*xd12 + Mloop0**2*(xd0 + xd1 + xd2) - 1.*sman*(xd2 + xd22) + m1**2*(xd12 + xd2 + xd22) + 
     &             m2**2*(xd12 + xd13 + xd2 + xd22 + xd23 + xd3)))))
        Bv0L = g0L*(-1.*g1R*m3*(g2R*g3L*m1*xd1 + g2L*g3L*Mloop2*(xd1 + xd2) + g2L*g3R*Mloop3*xd3) + 
     &     g1L*Mloop0*(g2R*g3L*m1*(xd0 + xd2 + xd3) + g2L*(-1.*g3R*Mloop3*xd3 + g3L*Mloop2*(xd0 + xd3)))) + 
     &  g0R*(-1.*g1R*g2L*g3R*m3*Mloop0*(xd0 + xd1 + xd2 + xd3) + g1L*(g2R*m1*(g3R*Mloop2*xd2 + g3L*Mloop3*xd3) + 
     &        g2L*(g3L*Mloop2*Mloop3*xd3 + g3R*(xc0 - 1.*xc2 - 4.*xd00 + uman*xd13 + tman*xd2 - 1.*m1**2*xd23 + (sman + uman)*xd23 - 1.*m3**2*(-1.*xd1 + xd13 + xd23) + 
     &              Mloop0**2*(xd0 - 1.*xd3) - 1.*m2**2*(xd13 + 2.*(xd23 + xd33))))))
        Bv0R = g0R*(-1.*g1L*m3*(g2L*g3R*m1*xd1 + g2R*g3R*Mloop2*(xd1 + xd2) + g2R*g3L*Mloop3*xd3) + 
     &     g1R*Mloop0*(g2L*g3R*m1*(xd0 + xd2 + xd3) + g2R*(-1.*g3L*Mloop3*xd3 + g3R*Mloop2*(xd0 + xd3)))) + 
     &  g0L*(-1.*g1L*g2R*g3L*m3*Mloop0*(xd0 + xd1 + xd2 + xd3) + g1R*(g2L*m1*(g3L*Mloop2*xd2 + g3R*Mloop3*xd3) + 
     &        g2R*(g3R*Mloop2*Mloop3*xd3 + g3L*(xc0 - 1.*xc2 - 4.*xd00 + uman*xd13 + tman*xd2 - 1.*m1**2*xd23 + (sman + uman)*xd23 - 1.*m3**2*(-1.*xd1 + xd13 + xd23) + 
     &              Mloop0**2*(xd0 - 1.*xd3) - 1.*m2**2*(xd13 + 2.*(xd23 + xd33))))))
        Bv1L = 2.*(g0R*g2L*g3R*(g1L*m3*(xd12 + xd2 + xd22) - 1.*g1R*Mloop0*xd23) + g0L*g1R*(g2R*g3L*m1*(xd2 + xd22) + g2L*(g3L*Mloop2*xd2 - 1.*g3R*Mloop3*xd23)))
        Bv1R = 2.*(g0L*g2R*g3L*(g1R*m3*(xd12 + xd2 + xd22) - 1.*g1L*Mloop0*xd23) + g0R*g1L*(g2L*g3R*m1*(xd2 + xd22) + g2R*(g3R*Mloop2*xd2 - 1.*g3L*Mloop3*xd23)))
        Bv2L = 2.*(g0L*g1R*(-1.*g2R*g3L*m1*xd23 + g2L*g3R*Mloop3*xd33) + g0R*g2L*g3R*(-1.*g1L*m3*(xd13 + xd23) + g1R*Mloop0*(xd3 + xd33)))
        Bv2R = 2.*(g0R*g1L*(-1.*g2L*g3R*m1*xd23 + g2R*g3L*Mloop3*xd33) + g0L*g2R*g3L*(-1.*g1R*m3*(xd13 + xd23) + g1L*Mloop0*(xd3 + xd33)))
        Bv3L = 2.*(g0L*g1R*(-1.*g2R*g3L*m1*(xd12 + xd2 + xd22) + g2L*g3R*Mloop3*(xd13 + xd23 + xd3)) + 
     &    g0R*g2L*g3R*(-1.*g1L*m3*(xd1 + xd11 + 2.*xd12 + xd2 + xd22) + g1R*Mloop0*(xd0 + xd1 + xd13 + xd2 + xd23 + xd3)))
        Bv3R = 2.*(g0R*g1L*(-1.*g2L*g3R*m1*(xd12 + xd2 + xd22) + g2R*g3L*Mloop3*(xd13 + xd23 + xd3)) + 
     &    g0L*g2R*g3L*(-1.*g1R*m3*(xd1 + xd11 + 2.*xd12 + xd2 + xd22) + g1L*Mloop0*(xd0 + xd1 + xd13 + xd2 + xd23 + xd3)))


c ==========================
c B.H.: Box test t-channel with gluino (to be deleted)

c       xd31 = xd13
c       xd32 = xd23
c       xd21 = xd12
c   
c       bjBs0L = g0L*(g1L*m3*MLoop0*(-1.*g2R*g3L*m1*xd1 + g2L*(-1.*g3L*MLoop2*(xd1 + xd2) + g3R*MLoop3*(xd0 + xd1 + xd2))) + 
c    &     g1R*(g2L*(-1.*g3R*MLoop3*(xc0 + MLoop0**2*xd0 - 2.*xd00 + tman*xd2 + m2**2*xd3 - 1.*uman*xd3 + 
c    &              m3**2*(xd1 + xd3)) + g3L*MLoop2*(xc0 + MLoop0**2*xd0 + m1**2*xd2 - 1.*sman*xd2 + m2**2*(xd2 + xd3)))
c    &          + g2R*m1*(-1.*g3R*MLoop2*MLoop3*xd2 + 
c    &           g3L*(xc0 + xc1 + MLoop0**2*(xd0 + xd2) + m1**2*(xd2 + xd22) - 1.*sman*(xd2 + xd22) + 
c    &              m3**2*(xd2 + xd21 + xd22) - 1.*uman*(xd2 + xd21 + xd22) + 
c    &              m2**2*(2.*xd2 + xd21 + 2.*(xd22 + xd23) + xd3))))) + 
c    &  g0R*(g1L*m3*(-1.*g2R*g3L*m1*MLoop3*xd1 + 
c    &        g2L*(-1.*g3L*MLoop2*MLoop3*(xd1 + xd2) + 
c    &           g3R*(-1.*xc2 - 1.*sman*xd12 + MLoop0**2*(xd0 + xd1 + xd2) - 1.*sman*(xd2 + xd22) + 
c    &              m1**2*(xd12 + xd2 + xd22) + m3**2*(xd1 + xd11 + xd12 + xd2 + xd21 + xd22) - 
c    &              1.*uman*(xd1 + xd11 + xd12 + xd2 + xd21 + xd22) + 
c    &              m2**2*(xd1 + xd11 + 2.*xd12 + 2.*xd13 + 2.*xd2 + xd21 + 2.*(xd22 + xd23) + xd3)))) - 
c    &     1.*g1R*MLoop0*(-1.*g2R*m1*(-1.*g3R*MLoop2*xd2 + g3L*MLoop3*(xd0 + xd2)) + 
c    &        g2L*(-1.*g3L*MLoop2*MLoop3*xd0 + g3R*
c    &            (xc0 + MLoop0**2*xd0 - 1.*uman*xd0 - 2.*xd00 - 1.*uman*xd1 + tman*xd2 - 1.*uman*xd2 - 1.*uman*xd3 + 
c    &              m3**2*(xd0 + 2.*xd1 + xd2 + xd3) + m2**2*(xd0 + xd1 + xd2 + 2.*xd3)))))
c       bjBs0R = g0R*(g1R*m3*MLoop0*(-1.*g2L*g3R*m1*xd1 + g2R*(-1.*g3R*MLoop2*(xd1 + xd2) + g3L*MLoop3*(xd0 + xd1 + xd2))) + 
c    &     g1L*(g2R*(-1.*g3L*MLoop3*(xc0 + MLoop0**2*xd0 - 2.*xd00 + tman*xd2 + m2**2*xd3 - 1.*uman*xd3 + 
c    &              m3**2*(xd1 + xd3)) + g3R*MLoop2*(xc0 + MLoop0**2*xd0 + m1**2*xd2 - 1.*sman*xd2 + m2**2*(xd2 + xd3)))
c    &          + g2L*m1*(-1.*g3L*MLoop2*MLoop3*xd2 + 
c    &           g3R*(xc0 + xc1 + MLoop0**2*(xd0 + xd2) + m1**2*(xd2 + xd22) - 1.*sman*(xd2 + xd22) + 
c    &              m3**2*(xd2 + xd21 + xd22) - 1.*uman*(xd2 + xd21 + xd22) + 
c    &              m2**2*(2.*xd2 + xd21 + 2.*(xd22 + xd23) + xd3))))) + 
c    &  g0L*(g1R*m3*(-1.*g2L*g3R*m1*MLoop3*xd1 + 
c    &        g2R*(-1.*g3R*MLoop2*MLoop3*(xd1 + xd2) + 
c    &           g3L*(-1.*xc2 - 1.*sman*xd12 + MLoop0**2*(xd0 + xd1 + xd2) - 1.*sman*(xd2 + xd22) + 
c    &              m1**2*(xd12 + xd2 + xd22) + m3**2*(xd1 + xd11 + xd12 + xd2 + xd21 + xd22) - 
c    &              1.*uman*(xd1 + xd11 + xd12 + xd2 + xd21 + xd22) + 
c    &              m2**2*(xd1 + xd11 + 2.*xd12 + 2.*xd13 + 2.*xd2 + xd21 + 2.*(xd22 + xd23) + xd3)))) - 
c    &     1.*g1L*MLoop0*(-1.*g2L*m1*(-1.*g3L*MLoop2*xd2 + g3R*MLoop3*(xd0 + xd2)) + 
c    &        g2R*(-1.*g3R*MLoop2*MLoop3*xd0 + g3L*
c    &            (xc0 + MLoop0**2*xd0 - 1.*uman*xd0 - 2.*xd00 - 1.*uman*xd1 + tman*xd2 - 1.*uman*xd2 - 1.*uman*xd3 + 
c    &              m3**2*(xd0 + 2.*xd1 + xd2 + xd3) + m2**2*(xd0 + xd1 + xd2 + 2.*xd3)))))
c       bjBs1L = 2.*g0L*(-1.*g1R*g2L*g3R*m3*MLoop3*(xd2 + xd21 + xd22) + 
c    &     g1L*g3L*MLoop0*(g2L*MLoop2*xd2 + g2R*m1*(xd2 + xd22))) - 
c    &  2.*g0R*(g1R*g2L*g3R*m3*MLoop0*(xd2 + xd21 + xd22) + 
c    &     g1L*(-1.*g2R*g3L*m1*MLoop3*(xd2 + xd22) + 
c    &        g2L*(-1.*g3L*MLoop2*MLoop3*xd2 + g3R*(m2**2 + m3**2 - 1.*uman)*(xd12 + xd2 + xd22) + g3R*m2**2*xd32)))
c       bjBs1R = 2.*g0R*(-1.*g1L*g2R*g3L*m3*MLoop3*(xd2 + xd21 + xd22) + 
c    &     g1R*g3R*MLoop0*(g2R*MLoop2*xd2 + g2L*m1*(xd2 + xd22))) - 
c    &  2.*g0L*(g1L*g2R*g3L*m3*MLoop0*(xd2 + xd21 + xd22) + 
c    &     g1R*(-1.*g2L*g3R*m1*MLoop3*(xd2 + xd22) + 
c    &        g2R*(-1.*g3R*MLoop2*MLoop3*xd2 + g3L*(m2**2 + m3**2 - 1.*uman)*(xd12 + xd2 + xd22) + g3L*m2**2*xd32)))
c       bjBs2L = 2.*(g0L*(-1.*g1L*g3L*MLoop0*(g2L*MLoop2*(xd0 + xd3) + g2R*m1*(xd0 + xd2 + xd3 + xd32)) + 
c    &       g1R*m3*(g2R*g3L*m1*xd1 + g2L*(g3L*MLoop2*(xd1 + xd2) + g3R*MLoop3*(xd3 + xd31 + xd32)))) + 
c    &    g0R*(g1R*g2L*g3R*m3*MLoop0*(xd0 + xd1 + xd2 + xd3 + xd31 + xd32) - 
c    &       1.*g1L*(g2R*m1*(g3R*MLoop2*xd2 + g3L*MLoop3*(xd3 + xd32)) + 
c    &          g2L*(g3L*MLoop2*MLoop3*xd3 + g3R*
c    &              (xc0 + MLoop0**2*xd0 - 2.*xd00 + tman*xd2 + m3**2*(xd1 - 1.*xd13 - 1.*xd23) + uman*(xd13 + xd23) - 
c    &                1.*m2**2*(xd13 + xd23 + xd33))))))
c       bjBs2R = 2.*(g0R*(-1.*g1R*g3R*MLoop0*(g2R*MLoop2*(xd0 + xd3) + g2L*m1*(xd0 + xd2 + xd3 + xd32)) + 
c    &       g1L*m3*(g2L*g3R*m1*xd1 + g2R*(g3R*MLoop2*(xd1 + xd2) + g3L*MLoop3*(xd3 + xd31 + xd32)))) + 
c    &    g0L*(g1L*g2R*g3L*m3*MLoop0*(xd0 + xd1 + xd2 + xd3 + xd31 + xd32) - 
c    &       1.*g1R*(g2L*m1*(g3L*MLoop2*xd2 + g3R*MLoop3*(xd3 + xd32)) + 
c    &          g2R*(g3R*MLoop2*MLoop3*xd3 + g3L*
c    &              (xc0 + MLoop0**2*xd0 - 2.*xd00 + tman*xd2 + m3**2*(xd1 - 1.*xd13 - 1.*xd23) + uman*(xd13 + xd23) - 
c    &                1.*m2**2*(xd13 + xd23 + xd33))))))
c       bjBs3L = 2.*(g0L*(g1R*g2L*g3R*m3*MLoop3*(xd1 + xd11 + xd12 + xd2 + xd21 + xd22) - 
c    &       1.*g1L*MLoop0*(g2L*g3R*MLoop3*(xd0 + xd1 + xd2) + g2R*g3L*m1*(xd12 + xd2 + xd22))) + 
c    &    g0R*(g1R*g2L*g3R*m3*MLoop0*(xd1 + xd11 + xd12 + xd2 + xd21 + xd22) + 
c    &       g1L*(-1.*g2R*g3L*m1*MLoop3*(xd12 + xd2 + xd22) + 
c    &          g2L*g3R*(xc2 + sman*xd12 - 1.*MLoop0**2*(xd0 + xd1 + xd2) + sman*(xd2 + xd22) - 
c    &             1.*m1**2*(xd12 + xd2 + xd22) - 
c    &             1.*m2**2*(xd12 + 2.*xd13 + xd2 + xd22 + 2.*xd23 + xd3 - 1.*xd31 - 1.*xd32)))))
c       bjBs3R = 2.*(g0R*(g1L*g2R*g3L*m3*MLoop3*(xd1 + xd11 + xd12 + xd2 + xd21 + xd22) - 
c    &       1.*g1R*MLoop0*(g2R*g3L*MLoop3*(xd0 + xd1 + xd2) + g2L*g3R*m1*(xd12 + xd2 + xd22))) + 
c    &    g0L*(g1L*g2R*g3L*m3*MLoop0*(xd1 + xd11 + xd12 + xd2 + xd21 + xd22) + 
c    &       g1R*(-1.*g2L*g3R*m1*MLoop3*(xd12 + xd2 + xd22) + 
c    &          g2R*g3L*(xc2 + sman*xd12 - 1.*MLoop0**2*(xd0 + xd1 + xd2) + sman*(xd2 + xd22) - 
c    &             1.*m1**2*(xd12 + xd2 + xd22) - 
c    &             1.*m2**2*(xd12 + 2.*xd13 + xd2 + xd22 + 2.*xd23 + xd3 - 1.*xd31 - 1.*xd32)))))
c       bjBv0L = g0L*(-1.*g1R*m3*(g2R*g3L*m1*xd1 + g2L*g3L*MLoop2*(xd1 + xd2) + g2L*g3R*MLoop3*xd3) + 
c    &     g1L*MLoop0*(g2R*g3L*m1*(xd0 + xd2 + xd3) + g2L*(-1.*g3R*MLoop3*xd3 + g3L*MLoop2*(xd0 + xd3)))) + 
c    &  g0R*(-1.*g1R*g2L*g3R*m3*MLoop0*(xd0 + xd1 + xd2 + xd3) + 
c    &     g1L*(g2R*m1*(g3R*MLoop2*xd2 + g3L*MLoop3*xd3) + 
c    &        g2L*(g3L*MLoop2*MLoop3*xd3 + g3R*(xc0 - 1.*xc2 - 4.*xd00 + tman*xd2 + MLoop0**2*(xd0 - 1.*xd3) + 
c    &              uman*xd31 + m3**2*(xd1 - 1.*xd31 - 1.*xd32) + (-1.*m1**2 + sman + uman)*xd32 - 
c    &              1.*m2**2*(xd31 + 2.*(xd32 + xd33))))))
c       bjBv0R = g0R*(-1.*g1L*m3*(g2L*g3R*m1*xd1 + g2R*g3R*MLoop2*(xd1 + xd2) + g2R*g3L*MLoop3*xd3) + 
c    &     g1R*MLoop0*(g2L*g3R*m1*(xd0 + xd2 + xd3) + g2R*(-1.*g3L*MLoop3*xd3 + g3R*MLoop2*(xd0 + xd3)))) + 
c    &  g0L*(-1.*g1L*g2R*g3L*m3*MLoop0*(xd0 + xd1 + xd2 + xd3) + 
c    &     g1R*(g2L*m1*(g3L*MLoop2*xd2 + g3R*MLoop3*xd3) + 
c    &        g2R*(g3R*MLoop2*MLoop3*xd3 + g3L*(xc0 - 1.*xc2 - 4.*xd00 + tman*xd2 + MLoop0**2*(xd0 - 1.*xd3) + 
c    &              uman*xd31 + m3**2*(xd1 - 1.*xd31 - 1.*xd32) + (-1.*m1**2 + sman + uman)*xd32 - 
c    &              1.*m2**2*(xd31 + 2.*(xd32 + xd33))))))
c       bjBv1L = 2.*(g0R*g2L*g3R*(g1L*m3*(xd12 + xd2 + xd22) - 1.*g1R*MLoop0*xd23) + 
c    &    g0L*g1R*(g2R*g3L*m1*(xd2 + xd22) + g2L*(g3L*MLoop2*xd2 - 1.*g3R*MLoop3*xd23)))
c       bjBv1R = 2.*(g0L*g2R*g3L*(g1R*m3*(xd12 + xd2 + xd22) - 1.*g1L*MLoop0*xd23) + 
c    &    g0R*g1L*(g2L*g3R*m1*(xd2 + xd22) + g2R*(g3R*MLoop2*xd2 - 1.*g3L*MLoop3*xd23)))
c       bjBv2L = -2.*g0L*g1R*g2R*g3L*m1*xd23 + 2.*g2L*g3R*
c    &   (g0L*g1R*MLoop3*xd33 + g0R*(-1.*g1L*m3*(xd13 + xd23) + g1R*MLoop0*(xd3 + xd33)))
c       bjBv2R = -2.*g0R*g1L*g2L*g3R*m1*xd23 + 2.*g2R*g3L*
c    &   (g0R*g1L*MLoop3*xd33 + g0L*(-1.*g1R*m3*(xd13 + xd23) + g1L*MLoop0*(xd3 + xd33)))
c       bjBv3L = 2.*(g0L*g1R*(-1.*g2R*g3L*m1*(xd2 + xd21 + xd22) + g2L*g3R*MLoop3*(xd13 + xd23 + xd3)) + 
c    &    g0R*g2L*g3R*(-1.*g1L*m3*(xd1 + xd11 + xd12 + xd2 + xd21 + xd22) + 
c    &       g1R*MLoop0*(xd0 + xd1 + xd13 + xd2 + xd23 + xd3)))
c       bjBv3R = 2.*(g0R*g1L*(-1.*g2L*g3R*m1*(xd2 + xd21 + xd22) + g2R*g3L*MLoop3*(xd13 + xd23 + xd3)) + 
c    &    g0L*g2R*g3L*(-1.*g1R*m3*(xd1 + xd11 + xd12 + xd2 + xd21 + xd22) + 
c    &       g1L*MLoop0*(xd0 + xd1 + xd13 + xd2 + xd23 + xd3)))
c    
c   	
c   	write(*,*) ' '
c   	write(*,*) 'T-GLUINO'
c   	write(*,*) ' '
c   	write(*,*) 'dBs0L = ',(Bs0L - bjBs0L)
c   	write(*,*) 'dBs0R = ',(Bs0R - bjBs0R)
c   	write(*,*) 'dBs1L = ',(Bs1L - bjBs1L)
c   	write(*,*) 'dBs1R = ',(Bs1R - bjBs1R)
c   	write(*,*) 'dBs2L = ',(Bs2L - bjBs2L)
c   	write(*,*) 'dBs2R = ',(Bs2R - bjBs2R)
c   	write(*,*) 'dBs3L = ',(Bs3L - bjBs3L)
c   	write(*,*) 'dBs3R = ',(Bs3R - bjBs3R)
c   	                                      
c   	write(*,*) 'dBv0L = ',(Bv0L - bjBv0L)
c   	write(*,*) 'dBv0L = ',(Bv0R - bjBv0R)
c   	write(*,*) 'dBv1L = ',(Bv1L - bjBv1L)
c   	write(*,*) 'dBv1L = ',(Bv1R - bjBv1R)
c   	write(*,*) 'dBv2L = ',(Bv2L - bjBv2L)
c   	write(*,*) 'dBv2L = ',(Bv2R - bjBv2R)
c   	write(*,*) 'dBv3L = ',(Bv3L - bjBv3L)
c   	write(*,*) 'dBv3L = ',(Bv3R - bjBv3R)
c
c       write(*,*)
c
c   	write(*,*) 'dBs0L = ',(Bs0L - bjBs0L)/Bs0L
c   	write(*,*) 'dBs0R = ',(Bs0R - bjBs0R)/Bs0R
c   	write(*,*) 'dBs1L = ',(Bs1L - bjBs1L)/Bs1L
c   	write(*,*) 'dBs1R = ',(Bs1R - bjBs1R)/Bs1R
c   	write(*,*) 'dBs2L = ',(Bs2L - bjBs2L)/Bs2L
c   	write(*,*) 'dBs2R = ',(Bs2R - bjBs2R)/Bs2R
c   	write(*,*) 'dBs3L = ',(Bs3L - bjBs3L)/Bs3L
c   	write(*,*) 'dBs3R = ',(Bs3R - bjBs3R)/Bs3R
c
c   	write(*,*) 'dBv0L = ',(Bv0L - bjBv0L)/Bv0L
c   	write(*,*) 'dBv0R = ',(Bv0R - bjBv0R)/Bv0R
c   	write(*,*) 'dBv1L = ',(Bv1L - bjBv1L)/Bv1L
c   	write(*,*) 'dBv1R = ',(Bv1R - bjBv1R)/Bv1R
c   	write(*,*) 'dBv2L = ',(Bv2L - bjBv2L)/Bv2L
c   	write(*,*) 'dBv2R = ',(Bv2R - bjBv2R)/Bv2R
c   	write(*,*) 'dBv3L = ',(Bv3L - bjBv3L)/Bv3L
c   	write(*,*) 'dBv3R = ',(Bv3R - bjBv3R)/Bv3R
c
c   	write(*,*) ' '
c   	write(*,*) ' '

c ==========================

! --- Interference with s channel tree-level --

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        BLc = dconjg(BL(xtype,ftype))
        BRc = dconjg(BR(xtype,ftype))
        
        mq1c = Mf(pstype,psgen)
c        wquark = width(pstype)
        wquark = 0d0
        
        boxTSgluino = colorfactor*CF*NC*2*MsBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A1Lc,A1Rc,BLc,BRc,mq1c,xv)

        NeuQ2qx_M2BoxV = NeuQ2qx_M2BoxV + boxTSgluino


! --- Interference with t channel tree-level --

                do sumTc = 0, 1
                
        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
        mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
c        wquark = width(pstype)
        wquark = 0d0
 
        boxTTgluino = colorfactor*CF*NC*2*MtBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A2Lc,A2Rc,CLc,mq2c,xv)
     
        NeuQ2qx_M2BoxV = NeuQ2qx_M2BoxV + boxTTgluino

                enddo
        
! --- Interference with u channel tree-level --

                do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)
c        wquark = width(pstype)
        wquark = 0d0

        boxTUgluino = colorfactor*CF*NC*2* MuBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A3Lc,A3Rc,DLc,DRc,mq3c,xv)

        NeuQ2qx_M2BoxV = NeuQ2qx_M2BoxV + boxTUgluino


                enddo

        enddo


! ----- Four-Vertex box with gluon -----

        do sumT = 1, 2

        g0 = CLL(sumT,sfmeigen,sftype,sfgen,xtype)
        g1L = AL(chimeigen,sumT,ftype,fgen)
        g1R = AR(chimeigen,sumT,ftype,fgen)
        g2 = B(1,ftype)

        Mloop0 = 0d0
        Mloop1 = MSf(sumT,ftype,fgen)
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
        
        boxFVS = colorfactor*CF*NC*2*MsBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A1Lc,A1Rc,BLc,BRc,mq1c,xv)

        NeuQ2qx_M2BoxV = NeuQ2qx_M2BoxV + boxFVS


! --- Interference with t channel tree-level --

                do sumTc = 0, 1
                
        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        CLc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))
 
        mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)
c        wquark = width(pstype)
        wquark = 0d0
 
        boxFVT = colorfactor*CF*NC*2*MtBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A2Lc,A2Rc,CLc,mq2c,xv)
     
        NeuQ2qx_M2BoxV = NeuQ2qx_M2BoxV + boxFVT

                enddo
        

! --- Interference with u channel tree-level --

                do sumUc = 0,sumUcend

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        DLc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        DRc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)
c        wquark = width(pstype)
        wquark = 0d0

        boxFVU = colorfactor*CF*NC*2*MuBoxV(Bs0L,Bs0R,Bs1L,Bs1R,Bs2L,Bs2R,Bs3L,Bs3R,Bv0L,Bv0R,Bv1L,Bv1R,Bv2L,Bv2R,Bv3L,Bv3R,A3Lc,A3Rc,DLc,DRc,mq3c,xv)

        NeuQ2qx_M2BoxV = NeuQ2qx_M2BoxV + boxFVU

                enddo

        enddo
        
	
        end
