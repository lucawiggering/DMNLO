	double precision function MssV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A1Lc,A1Rc,BLc,BRc,mq1,mq1c,xv,wquark)
        
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

        double complex A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A1Lc,A1Rc,BLc,BRc
        double precision mq1,mq1c,xv,wquark
              
        MssV = 0d0
        MssV = dreal((A1Rc*A1vL*m1*(BLc*((-2*BvgL*m3**2 - Bvm1L*m3**3 + Bvm2L*m3**3 + 2*BvgL*m4**2 + 3*Bvm1L*m3*m4**2 + Bvm2L*m3*m4**2 + Bsm1R*(-3*m3**2 + m4**2) + 
     &           Bsm2R*(m3**2 + m4**2) + ((3*Bvm1R - Bvm2R)*m3**2 - (Bvm1R + Bvm2R)*m4**2)*mq1)*PsR + 
     &        (2*BsgL*(-m3**2 + m4**2) + Bsm2L*m3*(m3**2 + m4**2) - Bsm1L*(m3**3 - 3*m3*m4**2) + 
     &           (m3**2*(2*BvgL + (Bvm1L - Bvm2L)*m3) - (2*BvgL + (3*Bvm1L + Bvm2L)*m3)*m4**2)*mq1)*PvL)*sman + 
     &     BRc*mq1c*((Bsm2R*m3*(m3**2 + m4**2) + (m3**2*(2*BvgR + (Bvm1R - Bvm2R)*m3) - (2*BvgR + (3*Bvm1R + Bvm2R)*m3)*m4**2)*mq1)*PsR + 
     &        ((-(Bsm2R*m3) + Bvm1L*(-3*m3**2 + m4**2) + Bvm2L*(m3**2 + m4**2) + 2*BvgR*mq1 + (3*Bvm1R + Bvm2R)*m3*mq1)*PsR + 
     &           (-2*BvgR*m3**2 - Bvm1R*m3**3 + Bvm2R*m3**3 + 2*BvgR*m4**2 + 3*Bvm1R*m3*m4**2 + Bvm2R*m3*m4**2 + Bsm1L*(-3*m3**2 + m4**2) + 
     &              Bsm2L*(m3**2 + m4**2) + ((3*Bvm1L - Bvm2L)*m3**2 - (Bvm1L + Bvm2L)*m4**2)*mq1)*PvL)*sman - 2*BsgR*PsR*(m3**2 - m4**2 + sman) - 
     &        Bsm1R*m3*PsR*(m3**2 - 3*m4**2 + 3*sman)) - sman**2*
     &      (BRc*mq1c*((Bsm1L + Bsm2L + 2*BvgR + 3*Bvm1R*m3 + Bvm2R*m3)*PvL + (Bvm1L + Bvm2L)*(PsR - mq1*PvL)) + 
     &        BLc*(Bsm1R*PsR + Bsm2R*PsR + 2*BvgL*PsR + 3*Bvm1L*m3*PsR + Bvm2L*m3*PsR - Bvm1R*mq1*PsR - Bvm2R*mq1*PsR + 2*BsgL*PvL + 3*Bsm1L*m3*PvL + 
     &           Bsm2L*m3*PvL + 3*Bvm1R*m3**2*PvL - Bvm2R*m3**2*PvL - Bvm1R*m4**2*PvL - Bvm2R*m4**2*PvL - 2*BvgL*mq1*PvL - 3*Bvm1L*m3*mq1*PvL - 
     &           Bvm2L*m3*mq1*PvL + (Bvm1R + Bvm2R)*PvL*sman)) + 
     &     (m3**2 + m4**2 - sman)*(BRc*(-(Bsm2R*m3*m4**2) + BsgR*(m3 - m4)*(m3 + m4) + Bsm1R*(m3**3 - 2*m3*m4**2) - BvgR*m3**2*mq1 - Bvm1R*m3**3*mq1 + 
     &           BvgR*m4**2*mq1 + 2*Bvm1R*m3*m4**2*mq1 + Bvm2R*m3*m4**2*mq1)*mq1c*PsR + 
     &        (BLc*(Bsm1R*m3**2 + Bvm1L*m3**3 - Bsm2R*m4**2 - 2*Bvm1L*m3*m4**2 - Bvm2L*m3*m4**2 + BvgL*(m3 - m4)*(m3 + m4) - Bvm1R*m3**2*mq1 + 
     &              Bvm2R*m4**2*mq1)*PsR + BLc*(-(Bsm2L*m3*m4**2) + BsgL*(m3 - m4)*(m3 + m4) + Bsm1L*(m3**3 - 2*m3*m4**2) - BvgL*m3**2*mq1 - Bvm1L*m3**3*mq1 + 
     &              BvgL*m4**2*mq1 + 2*Bvm1L*m3*m4**2*mq1 + Bvm2L*m3*m4**2*mq1)*PvL + 
     &           BRc*mq1c*((-BsgR - Bvm2L*m4**2 + BvgR*mq1 + m3*(-Bsm1R + Bvm1L*m3 + Bvm1R*mq1))*PsR + 
     &              (Bsm1L*m3**2 + Bvm1R*m3**3 - Bsm2L*m4**2 - 2*Bvm1R*m3*m4**2 - Bvm2R*m3*m4**2 + BvgR*(m3 - m4)*(m3 + m4) - Bvm1L*m3**2*mq1 + 
     &                 Bvm2L*m4**2*mq1)*PvL))*sman - (BLc*(Bsm1R + BvgL + Bvm1L*m3 - Bvm1R*mq1)*PsR + 
     &           BLc*(BsgL + Bsm1L*m3 - Bvm1R*m3**2 + Bvm2R*m4**2 - (BvgL + Bvm1L*m3)*mq1)*PvL + 
     &           BRc*mq1c*(Bvm1L*PsR + (Bsm1L + BvgR + Bvm1R*m3 - Bvm1L*mq1)*PvL))*sman**2 - BLc*Bvm1R*PvL*sman**3)*xv - 
     &     BRc*(m3 - m4)*(m3 + m4)*(BsgR - BvgR*mq1)*mq1c*PsR*(-4 + m4**2*xv) + 
     &     BRc*mq1c*(-2*m3*m4**2*(Bsm1R + Bsm2R - (Bvm1R + Bvm2R)*mq1)*PsR - 
     &        2*(-2*BsgR*PsR + 2*BvgR*mq1*PsR + (4*BsgL*m3 - 2*BvgR*m3**2 + (2*BvgR + (Bvm1R + Bvm2R)*m3)*m4**2)*PvL + 4*BvgL*m3*(PsR - mq1*PvL))*sman + 
     &        4*BvgR*PvL*sman**2 + m4**2*(2*m3*m4**2*(Bsm1R + Bsm2R - (Bvm1R + Bvm2R)*mq1)*PsR + 
     &           (-(BsgR*PsR) + BvgR*mq1*PsR + (m3*(2*BsgL - BvgR*m3) + (BvgR + 2*(Bvm1R + Bvm2R)*m3)*m4**2)*PvL + 2*BvgL*m3*(PsR - mq1*PvL))*sman - 
     &           BvgR*PvL*sman**2)*xv) + BLc*sman*(-2*(-4*BvgR*m3*mq1*PsR - 2*BsgL*m3**2*PvL + 2*BsgL*m4**2*PvL + Bsm1L*m3*m4**2*PvL + Bsm2L*m3*m4**2*PvL + 
     &           Bvm1L*m3*m4**2*(PsR - mq1*PvL) + Bvm2L*m3*m4**2*(PsR - mq1*PvL) - 2*BsgL*PvL*sman + 4*BvgR*m3*PvL*sman) + 
     &        m4**2*(-2*BvgR*m3*mq1*PsR - BsgL*m3**2*PvL + BsgL*m4**2*PvL + 2*Bsm1L*m3*m4**2*PvL + 2*Bsm2L*m3*m4**2*PvL + 2*Bvm1L*m3*m4**2*(PsR - mq1*PvL) + 
     &           2*Bvm2L*m3*m4**2*(PsR - mq1*PvL) - BsgL*PvL*sman + 2*BvgR*m3*PvL*sman)*xv + 2*BsgR*m3*PsR*(-4 + m4**2*xv) - 
     &        BvgL*(PsR - mq1*PvL)*(m3**2 - m4**2 + sman)*(-4 + m4**2*xv))) + 
     &  A1Lc*A1vR*m1*(BRc*((-2*BvgR*m3**2 - Bvm1R*m3**3 + Bvm2R*m3**3 + 2*BvgR*m4**2 + 3*Bvm1R*m3*m4**2 + Bvm2R*m3*m4**2 + Bsm1L*(-3*m3**2 + m4**2) + 
     &           Bsm2L*(m3**2 + m4**2) + ((3*Bvm1L - Bvm2L)*m3**2 - (Bvm1L + Bvm2L)*m4**2)*mq1)*PsL + 
     &        (2*BsgR*(-m3**2 + m4**2) + Bsm2R*m3*(m3**2 + m4**2) - Bsm1R*(m3**3 - 3*m3*m4**2) + 
     &           (m3**2*(2*BvgR + (Bvm1R - Bvm2R)*m3) - (2*BvgR + (3*Bvm1R + Bvm2R)*m3)*m4**2)*mq1)*PvR)*sman - 
     &     BRc*sman**2*(Bsm1L*PsL + Bsm2L*PsL + 2*BvgR*PsL + 3*Bvm1R*m3*PsL + Bvm2R*m3*PsL - Bvm1L*mq1*PsL - Bvm2L*mq1*PsL + 2*BsgR*PvR + 3*Bsm1R*m3*PvR + 
     &        Bsm2R*m3*PvR + 3*Bvm1L*m3**2*PvR - Bvm2L*m3**2*PvR - Bvm1L*m4**2*PvR - Bvm2L*m4**2*PvR - 2*BvgR*mq1*PvR - 3*Bvm1R*m3*mq1*PvR - 
     &        Bvm2R*m3*mq1*PvR + (Bvm1L + Bvm2L)*PvR*sman) + BLc*mq1c*
     &      ((Bsm2L*m3*(m3**2 + m4**2) + (m3**2*(2*BvgL + (Bvm1L - Bvm2L)*m3) - (2*BvgL + (3*Bvm1L + Bvm2L)*m3)*m4**2)*mq1)*PsL + 
     &        (-(Bsm2L*m3*PsL) + Bvm1R*(-3*m3**2 + m4**2)*PsL + (Bvm2R*(m3**2 + m4**2) + 2*BvgL*mq1 + (3*Bvm1L + Bvm2L)*m3*mq1)*PsL + 
     &           (m3**2*(-3*Bsm1R + Bsm2R - 2*BvgL - Bvm1L*m3 + Bvm2L*m3) + Bsm1R*m4**2)*PvR)*sman - 2*BsgL*PsL*(m3**2 - m4**2 + sman) - 
     &        Bsm1L*m3*PsL*(m3**2 - 3*m4**2 + 3*sman)) + BRc*Bvm1L*PvR*sman**4*xv + 
     &     (BLc*(-(Bsm2L*m3**3*m4**2) - (BsgL + Bsm2L*m3)*m4**4 + Bsm1L*m3*(m3**2 - 2*m4**2)*(m3**2 + m4**2) + 
     &           (m3**2 + m4**2)*(-(m3**2*(BvgL + Bvm1L*m3)) + (BvgL + (2*Bvm1L + Bvm2L)*m3)*m4**2)*mq1)*mq1c*PsL + 
     &        (BLc*mq1c*((-2*BsgL*m3**2 + Bvm1R*m3**4 + Bsm2L*m3*m4**2 + Bvm1R*m3**2*m4**2 - Bvm2R*m3**2*m4**2 - Bvm2R*m4**4 + Bsm1L*m3*(-2*m3**2 + m4**2) + 
     &                 m3*(2*m3*(BvgL + Bvm1L*m3) - (Bvm1L + Bvm2L)*m4**2)*mq1)*PsL + 
     &              (m3**2 + m4**2)*(Bsm1R*m3**2 + Bvm1L*m3**3 - Bsm2R*m4**2 - 2*Bvm1L*m3*m4**2 - Bvm2L*m3*m4**2 + BvgL*(m3 - m4)*(m3 + m4) - 
     &                 Bvm1R*m3**2*mq1 + Bvm2R*m4**2*mq1)*PvR) + 
     &           BRc*(m3**2 + m4**2)*(Bsm1L*m3**2*PsL + (-((Bsm2L + Bvm2R*m3)*m4**2) + Bvm1R*(m3**3 - 2*m3*m4**2) - Bvm1L*m3**2*mq1 + Bvm2L*m4**2*mq1)*PsL + 
     &              (BsgR*(m3 - m4)*(m3 + m4) + m3*(Bsm1R*(m3**2 - 2*m4**2) - Bvm1R*m3**2*mq1 + m4**2*(-Bsm2R + (2*Bvm1R + Bvm2R)*mq1)))*PvR + 
     &              BvgR*(m3 - m4)*(m3 + m4)*(PsL - mq1*PvR)))*sman + 
     &        (BRc*(-2*BvgR*m3**2 - 2*Bvm1R*m3**3 + Bsm2L*m4**2 + Bvm1R*m3*m4**2 + Bvm2R*m3*m4**2 - Bsm1L*(2*m3**2 + m4**2) + 2*Bvm1L*m3**2*mq1 + 
     &              Bvm1L*m4**2*mq1 - Bvm2L*m4**2*mq1)*PsL + BRc*
     &            (-2*BsgR*m3**2 - 2*Bsm1R*m3**3 + Bvm1L*m3**4 + Bsm1R*m3*m4**2 + Bsm2R*m3*m4**2 + Bvm1L*m3**2*m4**2 - Bvm2L*m3**2*m4**2 - Bvm2L*m4**4 + 
     &              m3*(2*m3*(BvgR + Bvm1R*m3) - (Bvm1R + Bvm2R)*m4**2)*mq1)*PvR + 
     &           BLc*mq1c*((BsgL + Bsm1L*m3 - 2*Bvm1R*m3**2 - Bvm1R*m4**2 + Bvm2R*m4**2 - (BvgL + Bvm1L*m3)*mq1)*PsL + 
     &              (-2*BvgL*m3**2 - 2*Bvm1L*m3**3 + Bsm2R*m4**2 + Bvm1L*m3*m4**2 + Bvm2L*m3*m4**2 - Bsm1R*(2*m3**2 + m4**2) + 2*Bvm1R*m3**2*mq1 + 
     &                 Bvm1R*m4**2*mq1 - Bvm2R*m4**2*mq1)*PvR))*sman**2 + 
     &        (BRc*(Bsm1L + BvgR + Bvm1R*m3 - Bvm1L*mq1)*PsL + BRc*(BsgR + Bsm1R*m3 - 2*Bvm1L*m3**2 - Bvm1L*m4**2 + Bvm2L*m4**2 - (BvgR + Bvm1R*m3)*mq1)*PvR + 
     &           BLc*mq1c*(Bvm1R*PsL + (Bsm1R + BvgL + Bvm1L*m3 - Bvm1R*mq1)*PvR))*sman**3)*xv - 
     &     BLc*(m3 - m4)*(m3 + m4)*(BsgL - BvgL*mq1)*mq1c*PsL*(-4 + m4**2*xv) + 
     &     BLc*mq1c*(sman*(((Bsm2R + 2*BvgL + (3*Bvm1L + Bvm2L)*m3)*m4**2 - ((-3*Bvm1R + Bvm2R)*m3**2 + (Bvm1R + Bvm2R)*m4**2)*mq1)*PvR - 
     &           ((Bsm1R + Bsm2R + 2*BvgL + 3*Bvm1L*m3 + Bvm2L*m3)*PvR + (Bvm1R + Bvm2R)*(PsL - mq1*PvR))*sman) + BsgL*m3**4*PsL*xv) + 
     &     BLc*mq1c*(-2*m3*m4**2*(Bsm1L + Bsm2L - (Bvm1L + Bvm2L)*mq1)*PsL + 
     &        2*(2*BsgL*PsL - 2*BvgL*mq1*PsL - (4*BsgR*m3 - 2*BvgL*m3**2 + (2*BvgL + (Bvm1L + Bvm2L)*m3)*m4**2)*PvR - 4*BvgR*m3*(PsL - mq1*PvR))*sman + 
     &        4*BvgL*PvR*sman**2 + m4**2*(2*m3*m4**2*(Bsm1L + Bsm2L - (Bvm1L + Bvm2L)*mq1)*PsL + 
     &           (-(BsgL*PsL) + BvgL*mq1*PsL + (m3*(2*BsgR - BvgL*m3) + (BvgL + 2*(Bvm1L + Bvm2L)*m3)*m4**2)*PvR + 2*BvgR*m3*(PsL - mq1*PvR))*sman - 
     &           BvgL*PvR*sman**2)*xv) + BRc*sman*(-2*(-4*BvgL*m3*mq1*PsL - 2*BsgR*m3**2*PvR + 2*BsgR*m4**2*PvR + Bsm1R*m3*m4**2*PvR + Bsm2R*m3*m4**2*PvR + 
     &           Bvm1R*m3*m4**2*(PsL - mq1*PvR) + Bvm2R*m3*m4**2*(PsL - mq1*PvR) - 2*BsgR*PvR*sman + 4*BvgL*m3*PvR*sman) + 
     &        m4**2*(-2*BvgL*m3*mq1*PsL - BsgR*m3**2*PvR + BsgR*m4**2*PvR + 2*Bsm1R*m3*m4**2*PvR + 2*Bsm2R*m3*m4**2*PvR + 2*Bvm1R*m3*m4**2*(PsL - mq1*PvR) + 
     &           2*Bvm2R*m3*m4**2*(PsL - mq1*PvR) - BsgR*PvR*sman + 2*BvgL*m3*PvR*sman)*xv + 2*BsgL*m3*PsL*(-4 + m4**2*xv) - 
     &        BvgR*(PsL - mq1*PvR)*(m3**2 - m4**2 + sman)*(-4 + m4**2*xv))) + 
     &  A1Lc*A1sR*m1*(-((BRc*(Bvm1L + Bvm2L)*PsR + BRc*(Bsm1L + Bsm2L + 2*BvgR + 3*Bvm1R*m3 + Bvm2R*m3 - (Bvm1L + Bvm2L)*mq1)*PvL + 
     &          BLc*(Bvm1R + Bvm2R)*mq1c*PvL)*sman**2) + BLc*mq1c*
     &      (-3*Bsm1R*m3**2*PsR + Bsm2R*m3**2*PsR - 2*BvgL*m3**2*PsR - Bvm1L*m3**3*PsR + Bvm2L*m3**3*PsR + Bsm1R*m4**2*PsR + Bsm2R*m4**2*PsR + 
     &        2*BvgL*m4**2*PsR + 3*Bvm1L*m3*m4**2*PsR + Bvm2L*m3*m4**2*PsR + 3*Bvm1R*m3**2*mq1*PsR - Bvm2R*m3**2*mq1*PsR - Bvm1R*m4**2*mq1*PsR - 
     &        Bvm2R*m4**2*mq1*PsR - 2*BsgL*m3**2*PvL - Bsm1L*m3**3*PvL + Bsm2L*m3**3*PvL + 2*BsgL*m4**2*PvL + 3*Bsm1L*m3*m4**2*PvL + Bsm2L*m3*m4**2*PvL + 
     &        2*BvgL*m3**2*mq1*PvL + Bvm1L*m3**3*mq1*PvL - Bvm2L*m3**3*mq1*PvL - 2*BvgL*m4**2*mq1*PvL - 3*Bvm1L*m3*m4**2*mq1*PvL - Bvm2L*m3*m4**2*mq1*PvL + 
     &        (-((Bsm1R + Bsm2R + 2*BvgL + 3*Bvm1L*m3 + Bvm2L*m3 - (Bvm1R + Bvm2R)*mq1)*PsR) + 
     &           (-2*BsgL - 3*Bsm1L*m3 - Bsm2L*m3 - 3*Bvm1R*m3**2 + Bvm2R*m3**2 + Bvm1R*m4**2 + Bvm2R*m4**2 + 2*BvgL*mq1 + (3*Bvm1L + Bvm2L)*m3*mq1)*PvL)*sman)
     &       + 2*BRc*(-(m3*m4**2*(Bsm1R + Bsm2R - (Bvm1R + Bvm2R)*mq1)*PsR) + 
     &        (2*BsgR*PsR - 2*BvgR*mq1*PsR - (4*BsgL*m3 - 2*BvgR*m3**2 + (2*BvgR + (Bvm1R + Bvm2R)*m3)*m4**2)*PvL - 4*BvgL*m3*(PsR - mq1*PvL))*sman + 
     &        2*BvgR*PvL*sman**2) + BRc*((Bsm2R*m3*(m3**2 + m4**2) + (m3**2*(2*BvgR + (Bvm1R - Bvm2R)*m3) - (2*BvgR + (3*Bvm1R + Bvm2R)*m3)*m4**2)*mq1)*PsR + 
     &        ((-(Bsm2R*m3) + Bvm1L*(-3*m3**2 + m4**2) + Bvm2L*(m3**2 + m4**2) + 2*BvgR*mq1 + (3*Bvm1R + Bvm2R)*m3*mq1)*PsR + 
     &           (-2*BvgR*m3**2 - Bvm1R*m3**3 + Bvm2R*m3**3 + 2*BvgR*m4**2 + 3*Bvm1R*m3*m4**2 + Bvm2R*m3*m4**2 + Bsm1L*(-3*m3**2 + m4**2) + 
     &              Bsm2L*(m3**2 + m4**2) + ((3*Bvm1L - Bvm2L)*m3**2 - (Bvm1L + Bvm2L)*m4**2)*mq1)*PvL)*sman - 2*BsgR*PsR*(m3**2 - m4**2 + sman) - 
     &        Bsm1R*m3*PsR*(m3**2 - 3*m4**2 + 3*sman)) + BRc*m4**2*
     &      (2*m3*m4**2*(Bsm1R + Bsm2R - (Bvm1R + Bvm2R)*mq1)*PsR + 
     &        (-(BsgR*PsR) + BvgR*mq1*PsR + (m3*(2*BsgL - BvgR*m3) + (BvgR + 2*(Bvm1R + Bvm2R)*m3)*m4**2)*PvL + 2*BvgL*m3*(PsR - mq1*PvL))*sman - 
     &        BvgR*PvL*sman**2)*xv + (m3**2 + m4**2 - sman)*(BLc*mq1c*
     &         (Bvm1L*m3**3*PsR - Bsm2R*m4**2*PsR - 2*Bvm1L*m3*m4**2*PsR - Bvm2L*m3*m4**2*PsR - Bvm1R*m3**2*mq1*PsR + Bvm2R*m4**2*mq1*PsR + BsgL*m3**2*PvL + 
     &           Bsm1L*m3**3*PvL - BsgL*m4**2*PvL - 2*Bsm1L*m3*m4**2*PvL - Bsm2L*m3*m4**2*PvL - Bvm1L*m3**3*mq1*PvL + 2*Bvm1L*m3*m4**2*mq1*PvL + 
     &           Bvm2L*m3*m4**2*mq1*PvL + Bsm1R*PsR*(m3**2 - sman) + BvgL*(PsR - mq1*PvL)*(m3**2 - m4**2 - sman) - 
     &           (-(Bvm1R*mq1*PsR) + (BsgL + m3*(Bsm1L - Bvm1R*m3) + Bvm2R*m4**2)*PvL + Bvm1L*m3*(PsR - mq1*PvL))*sman - Bvm1R*PvL*sman**2) + 
     &        BRc*((-(Bsm2R*m3*m4**2) + BsgR*(m3 - m4)*(m3 + m4) + Bsm1R*(m3**3 - 2*m3*m4**2) - BvgR*m3**2*mq1 - Bvm1R*m3**3*mq1 + BvgR*m4**2*mq1 + 
     &              2*Bvm1R*m3*m4**2*mq1 + Bvm2R*m3*m4**2*mq1)*PsR + 
     &           ((-BsgR - Bvm2L*m4**2 + BvgR*mq1 + m3*(-Bsm1R + Bvm1L*m3 + Bvm1R*mq1))*PsR + 
     &              (Bsm1L*m3**2 + Bvm1R*m3**3 - Bsm2L*m4**2 - 2*Bvm1R*m3*m4**2 - Bvm2R*m3*m4**2 + BvgR*(m3 - m4)*(m3 + m4) - Bvm1L*m3**2*mq1 + 
     &                 Bvm2L*m4**2*mq1)*PvL)*sman - (Bvm1L*PsR + (Bsm1L + BvgR + Bvm1R*m3 - Bvm1L*mq1)*PvL)*sman**2))*xv + 
     &     (-(BRc*(m3 - m4)*(m3 + m4)*(BsgR - BvgR*mq1)) + BLc*(2*BsgR*m3 + BvgL*(-m3**2 + m4**2) - 2*BvgR*m3*mq1)*mq1c)*PsR*(-4 + m4**2*xv) + 
     &     BLc*mq1c*(2*(-((Bsm1L + Bsm2L)*m3*m4**2) + 2*BsgL*(m3 - m4)*(m3 + m4) + 2*BvgL*(-m3**2 + m4**2)*mq1)*PvL + 
     &        4*(BvgL*PsR + (BsgL - 2*BvgR*m3 - BvgL*mq1)*PvL)*sman + 
     &        m4**2*((2*(Bsm1L + Bsm2L)*m3*m4**2 + BvgL*(m3 - m4)*(m3 + m4)*mq1)*PvL + (-(BvgL*PsR) + 2*BvgR*m3*PvL + BvgL*mq1*PvL)*sman - 
     &           BsgL*PvL*(m3**2 - m4**2 + sman))*xv + 2*Bvm1L*m3*m4**2*(PsR - mq1*PvL)*(-1 + m4**2*xv) + 2*Bvm2L*m3*m4**2*(PsR - mq1*PvL)*(-1 + m4**2*xv))) + 
     &  A1Rc*A1sL*m1*(-(BLc*((Bsm1R + Bsm2R + 2*BvgL + 3*Bvm1L*m3 + Bvm2L*m3)*PvR + (Bvm1R + Bvm2R)*(PsL - mq1*PvR))*sman**2) + 
     &     2*BLc*(-(m3*m4**2*(Bsm1L + Bsm2L - (Bvm1L + Bvm2L)*mq1)*PsL) + 
     &        (2*BsgL*PsL - 2*BvgL*mq1*PsL - (4*BsgR*m3 - 2*BvgL*m3**2 + (2*BvgL + (Bvm1L + Bvm2L)*m3)*m4**2)*PvR - 4*BvgR*m3*(PsL - mq1*PvR))*sman + 
     &        2*BvgL*PvR*sman**2) - BRc*mq1c*(2*BvgR*m3**2*PsL + Bvm1R*m3**3*PsL - Bvm2R*m3**3*PsL - 2*BvgR*m4**2*PsL - 3*Bvm1R*m3*m4**2*PsL - 
     &        Bvm2R*m3*m4**2*PsL - 3*Bvm1L*m3**2*mq1*PsL + Bvm2L*m3**2*mq1*PsL + Bvm1L*m4**2*mq1*PsL + Bvm2L*m4**2*mq1*PsL + 2*BsgR*m3**2*PvR + 
     &        Bsm1R*m3**3*PvR - Bsm2R*m3**3*PvR - 2*BsgR*m4**2*PvR - 3*Bsm1R*m3*m4**2*PvR - Bsm2R*m3*m4**2*PvR - 2*BvgR*m3**2*mq1*PvR - Bvm1R*m3**3*mq1*PvR + 
     &        Bvm2R*m3**3*mq1*PvR + 2*BvgR*m4**2*mq1*PvR + 3*Bvm1R*m3*m4**2*mq1*PvR + Bvm2R*m3*m4**2*mq1*PvR - Bsm2L*PsL*(m3**2 + m4**2 - sman) + 
     &        ((2*BvgR + 3*Bvm1R*m3 + Bvm2R*m3 - (Bvm1L + Bvm2L)*mq1)*PsL + 2*BsgR*PvR + m3*(3*Bsm1R + Bsm2R + 3*Bvm1L*m3 - Bvm2L*m3)*PvR)*sman + 
     &        Bsm1L*PsL*(3*m3**2 - m4**2 + sman)) + BLc*((Bsm2L*m3*(m3**2 + m4**2) + 
     &           (m3**2*(2*BvgL + (Bvm1L - Bvm2L)*m3) - (2*BvgL + (3*Bvm1L + Bvm2L)*m3)*m4**2)*mq1)*PsL + 
     &        ((-(Bsm2L*m3) + Bvm1R*(-3*m3**2 + m4**2) + Bvm2R*(m3**2 + m4**2) + 2*BvgL*mq1 + (3*Bvm1L + Bvm2L)*m3*mq1)*PsL + 
     &           (-2*BvgL*m3**2 - Bvm1L*m3**3 + Bvm2L*m3**3 + 2*BvgL*m4**2 + 3*Bvm1L*m3*m4**2 + Bvm2L*m3*m4**2 + Bsm1R*(-3*m3**2 + m4**2) + 
     &              Bsm2R*(m3**2 + m4**2) + ((3*Bvm1R - Bvm2R)*m3**2 - (Bvm1R + Bvm2R)*m4**2)*mq1)*PvR)*sman - 2*BsgL*PsL*(m3**2 - m4**2 + sman) - 
     &        Bsm1L*m3*PsL*(m3**2 - 3*m4**2 + 3*sman)) + BLc*(m3**2 + m4**2)*
     &      (Bsm1L*m3**3 - 2*Bsm1L*m3*m4**2 - Bsm2L*m3*m4**2 + BsgL*(m3 - m4)*(m3 + m4) - BvgL*m3**2*mq1 - Bvm1L*m3**3*mq1 + BvgL*m4**2*mq1 + 
     &        2*Bvm1L*m3*m4**2*mq1 + Bvm2L*m3*m4**2*mq1)*PsL*xv + 
     &     BLc*m4**2*(2*m3*m4**2*(Bsm1L + Bsm2L - (Bvm1L + Bvm2L)*mq1)*PsL + 
     &        (-(BsgL*PsL) + BvgL*mq1*PsL + (m3*(2*BsgR - BvgL*m3) + (BvgL + 2*(Bvm1L + Bvm2L)*m3)*m4**2)*PvR + 2*BvgR*m3*(PsL - mq1*PvR))*sman - 
     &        BvgL*PvR*sman**2)*xv + (BRc*mq1c*(Bsm1L*m3**2*m4**2*PsL - Bsm2L*m3**2*m4**2*PsL - Bvm2R*m3**3*m4**2*PsL - Bsm2L*m4**4*PsL - Bvm2R*m3*m4**4*PsL - 
     &           Bvm1L*m3**4*mq1*PsL - Bvm1L*m3**2*m4**2*mq1*PsL + Bvm2L*m3**2*m4**2*mq1*PsL + Bvm2L*m4**4*mq1*PsL + BsgR*m3**4*PvR + Bsm1R*m3**5*PvR - 
     &           Bsm1R*m3**3*m4**2*PvR - Bsm2R*m3**3*m4**2*PvR - BsgR*m4**4*PvR - 2*Bsm1R*m3*m4**4*PvR - Bsm2R*m3*m4**4*PvR + Bvm2R*m3**3*m4**2*mq1*PvR + 
     &           Bvm2R*m3*m4**4*mq1*PvR + BvgR*(PsL - mq1*PvR)*(-m4**4 + (m3**2 - sman)**2) + 
     &           Bvm1R*m3*(PsL - mq1*PvR)*(m3**2 - 2*m4**2 - sman)*(m3**2 + m4**2 - sman) + 
     &           (((Bsm2L + Bvm2R*m3)*m4**2 - Bsm1L*(2*m3**2 + m4**2) + (2*Bvm1L*m3**2 + (Bvm1L - Bvm2L)*m4**2)*mq1)*PsL + 
     &              (-2*BsgR*m3**2 - 2*Bsm1R*m3**3 + Bvm1L*m3**4 + Bsm1R*m3*m4**2 + Bsm2R*m3*m4**2 + Bvm1L*m3**2*m4**2 - Bvm2L*m3**2*m4**2 - Bvm2L*m4**4 - 
     &                 Bvm2R*m3*m4**2*mq1)*PvR)*sman + (Bsm1L*PsL - Bvm1L*mq1*PsL + (BsgR + m3*(Bsm1R - 2*Bvm1L*m3) + (-Bvm1L + Bvm2L)*m4**2)*PvR)*sman**2 + 
     &           Bvm1L*PvR*sman**3) + BLc*sman*(Bvm1R*m3**4*PsL + Bsm2L*m3*m4**2*PsL + Bvm1R*m3**2*m4**2*PsL - Bvm2R*m3**2*m4**2*PsL - Bvm2R*m4**4*PsL + 
     &           2*BvgL*m3**2*mq1*PsL + 2*Bvm1L*m3**3*mq1*PsL - Bvm1L*m3*m4**2*mq1*PsL - Bvm2L*m3*m4**2*mq1*PsL + Bsm1R*m3**4*PvR + BvgL*m3**4*PvR + 
     &           Bvm1L*m3**5*PvR + Bsm1R*m3**2*m4**2*PvR - Bsm2R*m3**2*m4**2*PvR - Bvm1L*m3**3*m4**2*PvR - Bvm2L*m3**3*m4**2*PvR - Bsm2R*m4**4*PvR - 
     &           BvgL*m4**4*PvR - 2*Bvm1L*m3*m4**4*PvR - Bvm2L*m3*m4**4*PvR - Bvm1R*m3**4*mq1*PvR - Bvm1R*m3**2*m4**2*mq1*PvR + Bvm2R*m3**2*m4**2*mq1*PvR + 
     &           Bvm2R*m4**4*mq1*PvR - ((BvgL + Bvm1L*m3)*mq1*PsL + (2*m3**2*(Bsm1R + BvgL + Bvm1L*m3) + (Bsm1R - Bsm2R - (Bvm1L + Bvm2L)*m3)*m4**2)*PvR + 
     &              Bvm1R*(2*m3**2 + m4**2)*(PsL - mq1*PvR) + Bvm2R*m4**2*(-PsL + mq1*PvR))*sman + 
     &           (Bvm1R*PsL + (Bsm1R + BvgL + Bvm1L*m3 - Bvm1R*mq1)*PvR)*sman**2 + BsgL*PsL*(-2*m3**2 + sman) + Bsm1L*m3*PsL*(-2*m3**2 + m4**2 + sman)))*xv + 
     &     (-(BLc*(m3 - m4)*(m3 + m4)*(BsgL - BvgL*mq1)) + BRc*(2*BsgL*m3 + BvgR*(-m3**2 + m4**2) - 2*BvgL*m3*mq1)*mq1c)*PsL*(-4 + m4**2*xv) + 
     &     BRc*mq1c*(PvR*((2*BvgR + (3*Bvm1R + Bvm2R)*m3)*mq1 + (Bvm1L + Bvm2L)*(m4**2 - sman))*sman + Bsm1L*m3**4*PsL*xv) + 
     &     BRc*mq1c*(2*(-((Bsm1R + Bsm2R)*m3*m4**2) + 2*BsgR*(m3 - m4)*(m3 + m4) + 2*BvgR*(-m3**2 + m4**2)*mq1)*PvR + 
     &        4*(BvgR*PsL + (BsgR - 2*BvgL*m3 - BvgR*mq1)*PvR)*sman + 
     &        m4**2*((2*(Bsm1R + Bsm2R)*m3*m4**2 + BvgR*(m3 - m4)*(m3 + m4)*mq1)*PvR + (-(BvgR*PsL) + 2*BvgL*m3*PvR + BvgR*mq1*PvR)*sman - 
     &           BsgR*PvR*(m3**2 - m4**2 + sman))*xv + 2*Bvm1R*m3*m4**2*(PsL - mq1*PvR)*(-1 + m4**2*xv) + 2*Bvm2R*m3*m4**2*(PsL - mq1*PvR)*(-1 + m4**2*xv))) + 
     &  A1Rc*A1vR*(BLc*(-2*Bsm1L*m1**2*m3**3*PsL + 2*Bsm2L*m1**2*m3**3*PsL - Bsm1L*m3**5*PsL + Bsm2L*m3**5*PsL + 6*Bsm1L*m1**2*m3*m4**2*PsL + 
     &        2*Bsm2L*m1**2*m3*m4**2*PsL + 2*Bsm1L*m3**3*m4**2*PsL + 2*Bsm2L*m3**3*m4**2*PsL + 3*Bsm1L*m3*m4**4*PsL + Bsm2L*m3*m4**4*PsL - 
     &        2*BsgL*(m3 - m4)*(m3 + m4)*(2*m1**2 + m3**2 + m4**2)*PsL + 4*BvgL*m1**2*m3**2*mq1*PsL + 2*Bvm1L*m1**2*m3**3*mq1*PsL - 
     &        2*Bvm2L*m1**2*m3**3*mq1*PsL + 2*BvgL*m3**4*mq1*PsL + Bvm1L*m3**5*mq1*PsL - Bvm2L*m3**5*mq1*PsL - 4*BvgL*m1**2*m4**2*mq1*PsL - 
     &        6*Bvm1L*m1**2*m3*m4**2*mq1*PsL - 2*Bvm2L*m1**2*m3*m4**2*mq1*PsL - 2*Bvm1L*m3**3*m4**2*mq1*PsL - 2*Bvm2L*m3**3*m4**2*mq1*PsL - 
     &        2*BvgL*m4**4*mq1*PsL - 3*Bvm1L*m3*m4**4*mq1*PsL - Bvm2L*m3*m4**4*mq1*PsL - 2*BsgL*(m1**2 + m4**2)*PsL*sman + 
     &        (-(Bsm1L*m3*(3*m1**2 + m3**2 + 2*m4**2)*PsL) + (-(Bsm2L*m3*(m1**2 + m3**2)) - Bvm1R*(m3**4 + m1**2*(2*m3**2 - m4**2)) + 
     &              Bvm2R*(m3**4 + 2*m3**2*m4**2 + m1**2*(2*m3**2 + m4**2)) + 
     &              (2*BvgL*(m1**2 + m4**2) + m3*(Bvm2L*(m1**2 + m3**2) + Bvm1L*(3*m1**2 + m3**2 + 2*m4**2)))*mq1)*PsL + 
     &           (-((2*m1**2 + m3**2 + m4**2)*(m3**2*(2*BvgL + (Bvm1L - Bvm2L)*m3) - (2*BvgL + (3*Bvm1L + Bvm2L)*m3)*m4**2)) - 
     &              Bsm1R*(m3**4 + m1**2*(2*m3**2 - m4**2)) + Bsm2R*(m3**4 + 2*m3**2*m4**2 + m1**2*(2*m3**2 + m4**2)) + 
     &              ((Bvm1R - Bvm2R)*m3**2*(2*m1**2 + m3**2) - ((Bvm1R + Bvm2R)*m1**2 + 2*Bvm2R*m3**2)*m4**2)*mq1)*PvR)*sman) - 
     &     BRc*mq1c*((-2*Bvm1L*m3**2*(2*m1**2 + m3**2) + ((Bvm1L + Bvm2L)*m1**2 - 2*Bvm1L*m3**2)*m4**2 + (Bvm1L + Bvm2L)*m4**4)*mq1*PsL - 
     &        Bsm2L*(m1**2 + m4**2)*PsL*(m4**2 - sman) + (2*BvgR*(m1**2 + m3**2)*PsL + Bvm1R*m3*(3*m1**2 + 2*m3**2 + m4**2)*PsL + 
     &           (m1**2 + m4**2)*(Bvm2R*m3 - (Bvm1L + Bvm2L)*mq1)*PsL + 
     &           (2*BsgR*m1**2 + (3*Bsm1R + Bsm2R)*m1**2*m3 + 2*(BsgR + 2*Bvm1L*m1**2)*m3**2 + 2*Bsm1R*m3**3 + 2*Bvm1L*m3**4)*PvR)*sman + 
     &        Bsm1L*PsL*(2*m3**4 + 2*m3**2*m4**2 - m4**4 + m4**2*sman + m1**2*(4*m3**2 - m4**2 + sman))) + 
     &     BRc*mq1c*(-(Bsm1R*m3*m4**2*PvR*sman) - Bsm2R*m3*m4**2*PvR*sman + 2*BvgR*m1**2*mq1*PvR*sman + 3*Bvm1R*m1**2*m3*mq1*PvR*sman + 
     &        Bvm2R*m1**2*m3*mq1*PvR*sman + 2*BvgR*m3**2*mq1*PvR*sman + 2*Bvm1R*m3**3*mq1*PvR*sman + Bvm1R*m3*m4**2*mq1*PvR*sman + 
     &        Bvm2R*m3*m4**2*mq1*PvR*sman + 2*Bsm1L*m3**2*PsL*tman + 2*BvgR*PsL*sman*tman + 2*Bvm1R*m3*PsL*sman*tman + 2*BsgR*PvR*sman*tman + 
     &        2*Bsm1R*m3*PvR*sman*tman - 2*BvgR*mq1*PvR*sman*tman - 2*Bvm1R*m3*mq1*PvR*sman*tman + 
     &        ((2*Bsm1L*m3**2 - (Bsm1L + Bsm2L)*m4**2)*PsL + ((Bsm1L + Bsm2L + (Bvm1R + Bvm2R)*m3)*PsL + m3*(Bsm1R + Bsm2R - (Bvm1R + Bvm2R)*mq1)*PvR)*sman)*
     &         uman + Bvm2L*(m4**2 - sman)*(PvR*sman*(m1**2 + m4**2 - uman) + mq1*PsL*uman) + 
     &        Bvm1L*(m1**2*PvR*(m4**2 - sman)*sman + (m4**2 - sman)*(PvR*sman*(m4**2 - uman) + mq1*PsL*uman) - 
     &           2*m3**2*(m4**2*PvR*sman + (mq1*PsL - PvR*sman)*(tman + uman)))) + 
     &     (BRc*mq1c*(m4**2*(-(Bsm2L*m4**2*(m1**2 + m4**2)) + Bsm1L*(m3**4 + 2*m3**2*m4**2 + m1**2*(2*m3**2 + m4**2)) - 
     &              (-(Bvm2L*m4**2*(m1**2 + m4**2)) + Bvm1L*(m3**4 + 2*m3**2*m4**2 + m1**2*(2*m3**2 + m4**2)))*mq1)*PsL - 
     &           (Bsm1L*(m1**2*m3**2 + 2*(m1**2 + m3**2)*m4**2 + m4**4)*PsL + 
     &              ((m1**2 + m4**2)*(m3**2*(BvgR + Bvm1R*m3) - (Bsm2L - BvgR + Bvm2R*m3)*m4**2) - 
     &                 (-(Bvm2L*m4**2*(m1**2 + m4**2)) + Bvm1L*(m1**2*m3**2 + 2*(m1**2 + m3**2)*m4**2 + m4**4))*mq1)*PsL + 
     &              (-(m3*(Bsm2R*m1**2 + Bvm1L*m3*(2*m1**2 + m3**2))*m4**2) - ((Bvm1L - Bvm2L)*m1**2 + Bsm2R*m3 + 2*Bvm1L*m3**2)*m4**4 + Bvm2L*m4**6 + 
     &                 Bsm1R*m3**3*(m1**2 + m4**2) + BsgR*(m1**2 + m4**2)*(m3**2 + m4**2) - 
     &                 (m1**2 + m4**2)*(m3**2*(BvgR + Bvm1R*m3) + (BvgR - Bvm2R*m3)*m4**2)*mq1)*PvR)*sman + 
     &           ((m1**2 + m4**2)*(Bsm1L + BvgR + Bvm1R*m3 - Bvm1L*mq1)*PsL + 
     &              (-(Bvm1L*m1**2*m3**2) - 2*Bvm1L*m1**2*m4**2 + Bvm2L*m1**2*m4**2 - 2*Bvm1L*m3**2*m4**2 - Bvm1L*m4**4 + Bvm2L*m4**4 + BsgR*(m1**2 + m4**2) + 
     &                 Bsm1R*m3*(m1**2 + m4**2) - (BvgR + Bvm1R*m3)*(m1**2 + m4**2)*mq1)*PvR)*sman**2 + Bvm1L*(m1**2 + m4**2)*PvR*sman**3) + 
     &        BLc*sman*(-3*Bsm1L*m1**2*m3**3*PsL + 2*Bvm1R*m1**2*m3**4*PsL - 2*Bsm1L*m3**5*PsL + Bvm1R*m3**6*PsL + 2*Bsm1L*m1**2*m3*m4**2*PsL + 
     &           Bsm2L*m1**2*m3*m4**2*PsL - 2*Bvm2R*m1**2*m3**2*m4**2*PsL + Bsm2L*m3**3*m4**2*PsL + Bvm1R*m3**4*m4**2*PsL - Bvm2R*m3**4*m4**2*PsL - 
     &           Bvm1R*m1**2*m4**4*PsL - Bvm2R*m1**2*m4**4*PsL + Bsm1L*m3*m4**4*PsL - Bvm1R*m3**2*m4**4*PsL - 2*Bvm2R*m3**2*m4**4*PsL + 
     &           3*BvgL*m1**2*m3**2*mq1*PsL + 3*Bvm1L*m1**2*m3**3*mq1*PsL + 2*BvgL*m3**4*mq1*PsL + 2*Bvm1L*m3**5*mq1*PsL - BvgL*m1**2*m4**2*mq1*PsL - 
     &           2*Bvm1L*m1**2*m3*m4**2*mq1*PsL - Bvm2L*m1**2*m3*m4**2*mq1*PsL + BvgL*m3**2*m4**2*mq1*PsL - Bvm2L*m3**3*m4**2*mq1*PsL - BvgL*m4**4*mq1*PsL - 
     &           Bvm1L*m3*m4**4*mq1*PsL + 2*Bsm1R*m1**2*m3**4*PvR + 2*BvgL*m1**2*m3**4*PvR + 2*Bvm1L*m1**2*m3**5*PvR + Bsm1R*m3**6*PvR + BvgL*m3**6*PvR + 
     &           Bvm1L*m3**7*PvR - 2*Bsm2R*m1**2*m3**2*m4**2*PvR - 2*Bvm1L*m1**2*m3**3*m4**2*PvR - 2*Bvm2L*m1**2*m3**3*m4**2*PvR + Bsm1R*m3**4*m4**2*PvR - 
     &           Bsm2R*m3**4*m4**2*PvR + BvgL*m3**4*m4**2*PvR - Bvm2L*m3**5*m4**2*PvR - Bsm1R*m1**2*m4**4*PvR - Bsm2R*m1**2*m4**4*PvR - 
     &           2*BvgL*m1**2*m4**4*PvR - 4*Bvm1L*m1**2*m3*m4**4*PvR - 2*Bvm2L*m1**2*m3*m4**4*PvR - Bsm1R*m3**2*m4**4*PvR - 2*Bsm2R*m3**2*m4**4*PvR - 
     &           BvgL*m3**2*m4**4*PvR - 3*Bvm1L*m3**3*m4**4*PvR - 2*Bvm2L*m3**3*m4**4*PvR - BvgL*m4**6*PvR - 2*Bvm1L*m3*m4**6*PvR - Bvm2L*m3*m4**6*PvR - 
     &           2*Bvm1R*m1**2*m3**4*mq1*PvR - Bvm1R*m3**6*mq1*PvR + 2*Bvm2R*m1**2*m3**2*m4**2*mq1*PvR - Bvm1R*m3**4*m4**2*mq1*PvR + 
     &           Bvm2R*m3**4*m4**2*mq1*PvR + Bvm1R*m1**2*m4**4*mq1*PvR + Bvm2R*m1**2*m4**4*mq1*PvR + Bvm1R*m3**2*m4**4*mq1*PvR + 2*Bvm2R*m3**2*m4**4*mq1*PvR + 
     &           (Bsm1L*m3*(m1**2 + m3**2)*PsL - (m1**2 + m3**2)*(-(Bvm2R*m4**2) + (BvgL + Bvm1L*m3)*mq1)*PsL + 
     &              (-3*Bvm1L*m1**2*m3**3 - 2*Bvm1L*m3**5 + Bsm2R*m1**2*m4**2 + 2*Bvm1L*m1**2*m3*m4**2 + Bvm2L*m1**2*m3*m4**2 + Bsm2R*m3**2*m4**2 + 
     &                 Bvm2L*m3**3*m4**2 + Bvm1L*m3*m4**4 - Bsm1R*m3**2*(3*m1**2 + 2*m3**2 + m4**2) + 
     &                 BvgL*(-2*m3**4 - m3**2*m4**2 + m4**4 + m1**2*(-3*m3**2 + m4**2)) - Bvm2R*(m1**2 + m3**2)*m4**2*mq1)*PvR - 
     &              Bvm1R*m3**2*(3*m1**2 + 2*m3**2 + m4**2)*(PsL - mq1*PvR))*sman + 
     &           (m1**2 + m3**2)*(Bvm1R*PsL + (Bsm1R + BvgL + Bvm1L*m3 - Bvm1R*mq1)*PvR)*sman**2 + 
     &           BsgL*PsL*(-2*m3**4 + m4**4 + m3**2*(-m4**2 + sman) + m1**2*(-3*m3**2 + m4**2 + sman))))*xv - 
     &     BLc*(Bsm1R*m1**2*PvR*sman**2 + Bsm2R*m1**2*PvR*sman**2 + 2*BvgL*m1**2*PvR*sman**2 + 3*Bvm1L*m1**2*m3*PvR*sman**2 + Bvm2L*m1**2*m3*PvR*sman**2 + 
     &        Bsm1R*m3**2*PvR*sman**2 + Bsm2R*m3**2*PvR*sman**2 + Bvm1L*m3**3*PvR*sman**2 + Bvm2L*m3**3*PvR*sman**2 + 2*BvgL*m4**2*PvR*sman**2 + 
     &        2*Bvm1L*m3*m4**2*PvR*sman**2 - 2*BsgL*m3**2*PsL*tman - Bsm1L*m3**3*PsL*tman + Bsm2L*m3**3*PsL*tman + 2*BsgL*m4**2*PsL*tman + 
     &        3*Bsm1L*m3*m4**2*PsL*tman + Bsm2L*m3*m4**2*PsL*tman + 2*BvgL*m3**2*mq1*PsL*tman + Bvm1L*m3**3*mq1*PsL*tman - Bvm2L*m3**3*mq1*PsL*tman - 
     &        2*BvgL*m4**2*mq1*PsL*tman - 3*Bvm1L*m3*m4**2*mq1*PsL*tman - Bvm2L*m3*m4**2*mq1*PsL*tman - Bsm1L*m3*PsL*sman*tman - Bsm2L*m3*PsL*sman*tman + 
     &        Bvm1L*m3*mq1*PsL*sman*tman + Bvm2L*m3*mq1*PsL*sman*tman - Bsm1R*m3**2*PvR*sman*tman + Bsm2R*m3**2*PvR*sman*tman - 2*BvgL*m3**2*PvR*sman*tman - 
     &        Bvm1L*m3**3*PvR*sman*tman + Bvm2L*m3**3*PvR*sman*tman + Bsm1R*m4**2*PvR*sman*tman + Bsm2R*m4**2*PvR*sman*tman + 2*BvgL*m4**2*PvR*sman*tman + 
     &        3*Bvm1L*m3*m4**2*PvR*sman*tman + Bvm2L*m3*m4**2*PvR*sman*tman - Bsm1R*PvR*sman**2*tman - Bsm2R*PvR*sman**2*tman - Bvm1L*m3*PvR*sman**2*tman - 
     &        Bvm2L*m3*PvR*sman**2*tman + Bvm1R*(PsL - mq1*PvR)*sman*(m1**2*sman + (m4**2 - sman)*tman + m3**2*(sman - tman - uman)) - 2*BsgL*m3**2*PsL*uman - 
     &        Bsm1L*m3**3*PsL*uman + Bsm2L*m3**3*PsL*uman + 2*BsgL*m4**2*PsL*uman + 3*Bsm1L*m3*m4**2*PsL*uman + Bsm2L*m3*m4**2*PsL*uman + 
     &        2*BvgL*m3**2*mq1*PsL*uman + Bvm1L*m3**3*mq1*PsL*uman - Bvm2L*m3**3*mq1*PsL*uman - 2*BvgL*m4**2*mq1*PsL*uman - 3*Bvm1L*m3*m4**2*mq1*PsL*uman - 
     &        Bvm2L*m3*m4**2*mq1*PsL*uman - 2*BsgL*PsL*sman*uman - 2*Bsm1L*m3*PsL*sman*uman + 2*BvgL*mq1*PsL*sman*uman + 2*Bvm1L*m3*mq1*PsL*sman*uman - 
     &        Bsm1R*m3**2*PvR*sman*uman + Bsm2R*m3**2*PvR*sman*uman - 2*BvgL*m3**2*PvR*sman*uman - Bvm1L*m3**3*PvR*sman*uman + Bvm2L*m3**3*PvR*sman*uman + 
     &        2*BvgL*m4**2*PvR*sman*uman + 3*Bvm1L*m3*m4**2*PvR*sman*uman + Bvm2L*m3*m4**2*PvR*sman*uman - 2*BvgL*PvR*sman**2*uman - 
     &        2*Bvm1L*m3*PvR*sman**2*uman + Bvm2R*(PsL - mq1*PvR)*sman*(m1**2*sman + (m4**2 - sman)*tman + m3**2*(sman + tman + uman)) - 
     &        (m3**2 + m4**2)*(2*m1**2 + m3**2 + m4**2)*(-(Bsm2L*m3*m4**2) + BsgL*(m3 - m4)*(m3 + m4) + Bsm1L*(m3**3 - 2*m3*m4**2) - BvgL*m3**2*mq1 - 
     &           Bvm1L*m3**3*mq1 + BvgL*m4**2*mq1 + 2*Bvm1L*m3*m4**2*mq1 + Bvm2L*m3*m4**2*mq1)*PsL*xv) + 
     &     BRc*mq1c*(-(((m4**4*(-Bsm2L + Bvm2L*mq1)*PsL - (-(m4**2*(Bsm2L + Bvm2R*m3 - Bvm2L*mq1)*PsL) + 
     &                   (Bsm1R*m3**3 + BsgR*(m3**2 + m4**2) + m4**2*(-(Bsm2R*m3) + Bvm2L*m4**2 + Bvm2R*m3*mq1))*PvR + Bvm1R*m3**3*(PsL - mq1*PvR) + 
     &                   BvgR*(m3**2 + m4**2)*(PsL - mq1*PvR))*sman + 
     &                ((BvgR + Bvm1R*m3)*PsL + (BsgR + Bsm1R*m3 + Bvm2L*m4**2 - (BvgR + Bvm1R*m3)*mq1)*PvR)*sman**2)*uman + 
     &             Bsm1L*PsL*(m4**2*(m3**2 + m4**2 - sman)*tman + (m3**2 - sman)*(m4**2 - sman)*uman) + 
     &             Bvm1L*(-(mq1*PsL) + PvR*sman)*(m4**2*(m3**2 + m4**2 - sman)*tman + (m3**2 - sman)*(m4**2 - sman)*uman))*xv) + 
     &        BsgL*m3*(2*m1**2 + m3**2 + m4**2)*PsL*(-4 + m4**2*xv) - BvgL*m3*(2*m1**2 + m3**2 + m4**2)*mq1*PsL*(-4 + m4**2*xv)) + 
     &     BRc*mq1c*(-(Bvm1L*m1**2*m4**2*mq1*PsL) - Bvm2L*m1**2*m4**2*mq1*PsL - Bvm1L*m3**2*m4**2*mq1*PsL - Bvm2L*m3**2*m4**2*mq1*PsL + 
     &        4*BvgR*m1**2*PsL*sman + 4*BvgR*m3**2*PsL*sman + 4*BsgR*m1**2*PvR*sman - 8*BvgL*m1**2*m3*PvR*sman + 4*BsgR*m3**2*PvR*sman - 
     &        4*BvgL*m3**3*PvR*sman + Bvm1L*m1**2*m4**2*PvR*sman + Bvm2L*m1**2*m4**2*PvR*sman - 4*BvgL*m3*m4**2*PvR*sman + Bvm1L*m3**2*m4**2*PvR*sman + 
     &        Bvm2L*m3**2*m4**2*PvR*sman - 4*BvgR*m1**2*mq1*PvR*sman - 4*BvgR*m3**2*mq1*PvR*sman + 4*BsgL*m3*PsL*tman - 4*BvgL*m3*mq1*PsL*tman + 
     &        Bvm1L*m4**2*mq1*PsL*tman + Bvm2L*m4**2*mq1*PsL*tman - 4*BvgR*PsL*sman*tman - 4*BsgR*PvR*sman*tman + 4*BvgL*m3*PvR*sman*tman - 
     &        Bvm1L*m4**2*PvR*sman*tman - Bvm2L*m4**2*PvR*sman*tman + 4*BvgR*mq1*PvR*sman*tman + 4*BsgL*m3*PsL*uman - 4*BvgL*m3*mq1*PsL*uman + 
     &        4*BvgL*m3*PvR*sman*uman + m4**2*(-(BvgR*m1**2*PsL*sman) - BvgR*m3**2*PsL*sman - BsgR*m1**2*PvR*sman + 2*BvgL*m1**2*m3*PvR*sman - 
     &           BsgR*m3**2*PvR*sman + BvgL*m3**3*PvR*sman + BvgL*m3*m4**2*PvR*sman + BvgR*m1**2*mq1*PvR*sman + BvgR*m3**2*mq1*PvR*sman + 
     &           Bvm1L*m4**2*(mq1*PsL - PvR*sman)*(m1**2 + m3**2 - tman) + Bvm2L*m4**2*(mq1*PsL - PvR*sman)*(m1**2 + m3**2 - tman) - BsgL*m3*PsL*tman + 
     &           BvgL*m3*mq1*PsL*tman + BvgR*PsL*sman*tman + BsgR*PvR*sman*tman - BvgL*m3*PvR*sman*tman - BvgR*mq1*PvR*sman*tman - 
     &           m3*(BsgL*PsL - BvgL*mq1*PsL + BvgL*PvR*sman)*uman)*xv - Bsm1L*m4**2*PsL*(m1**2 + m3**2 - tman)*(-1 + m4**2*xv) - 
     &        Bsm2L*m4**2*PsL*(m1**2 + m3**2 - tman)*(-1 + m4**2*xv)) + 
     &     BLc*(4*Bvm1L*m1**2*m3*m4**2*mq1*PsL + 4*Bvm2L*m1**2*m3*m4**2*mq1*PsL + 2*Bvm1L*m3**3*m4**2*mq1*PsL + 2*Bvm2L*m3**3*m4**2*mq1*PsL + 
     &        2*Bvm1L*m3*m4**4*mq1*PsL + 2*Bvm2L*m3*m4**4*mq1*PsL + 4*BsgL*m1**2*PsL*sman - 8*BvgR*m1**2*m3*PsL*sman - 4*BvgR*m3**3*PsL*sman + 
     &        4*BsgL*m4**2*PsL*sman - Bvm1R*m1**2*m4**2*PsL*sman - Bvm2R*m1**2*m4**2*PsL*sman - 4*BvgR*m3*m4**2*PsL*sman - Bvm1R*m3**2*m4**2*PsL*sman - 
     &        Bvm2R*m3**2*m4**2*PsL*sman - 4*BvgL*m1**2*mq1*PsL*sman - 4*BvgL*m4**2*mq1*PsL*sman - 8*BsgR*m1**2*m3*PvR*sman + 8*BvgL*m1**2*m3**2*PvR*sman - 
     &        4*BsgR*m3**3*PvR*sman + 4*BvgL*m3**4*PvR*sman - Bsm1R*m1**2*m4**2*PvR*sman - Bsm2R*m1**2*m4**2*PvR*sman - 8*BvgL*m1**2*m4**2*PvR*sman - 
     &        4*BsgR*m3*m4**2*PvR*sman - 4*Bvm1L*m1**2*m3*m4**2*PvR*sman - 4*Bvm2L*m1**2*m3*m4**2*PvR*sman - Bsm1R*m3**2*m4**2*PvR*sman - 
     &        Bsm2R*m3**2*m4**2*PvR*sman - 2*Bvm1L*m3**3*m4**2*PvR*sman - 2*Bvm2L*m3**3*m4**2*PvR*sman - 4*BvgL*m4**4*PvR*sman - 2*Bvm1L*m3*m4**4*PvR*sman - 
     &        2*Bvm2L*m3*m4**4*PvR*sman + 8*BvgR*m1**2*m3*mq1*PvR*sman + 4*BvgR*m3**3*mq1*PvR*sman + Bvm1R*m1**2*m4**2*mq1*PvR*sman + 
     &        Bvm2R*m1**2*m4**2*mq1*PvR*sman + 4*BvgR*m3*m4**2*mq1*PvR*sman + Bvm1R*m3**2*m4**2*mq1*PvR*sman + Bvm2R*m3**2*m4**2*mq1*PvR*sman + 
     &        4*BvgL*m1**2*PvR*sman**2 + 4*BvgL*m4**2*PvR*sman**2 - 4*BsgL*m3**2*PsL*tman + 4*BsgL*m4**2*PsL*tman + 4*BvgL*m3**2*mq1*PsL*tman - 
     &        4*BvgL*m4**2*mq1*PsL*tman - 2*Bvm1L*m3*m4**2*mq1*PsL*tman - 2*Bvm2L*m3*m4**2*mq1*PsL*tman + 4*BvgR*m3*PsL*sman*tman + 
     &        Bvm1R*m4**2*PsL*sman*tman + Bvm2R*m4**2*PsL*sman*tman + 4*BsgR*m3*PvR*sman*tman - 4*BvgL*m3**2*PvR*sman*tman + Bsm1R*m4**2*PvR*sman*tman + 
     &        Bsm2R*m4**2*PvR*sman*tman + 4*BvgL*m4**2*PvR*sman*tman + 2*Bvm1L*m3*m4**2*PvR*sman*tman + 2*Bvm2L*m3*m4**2*PvR*sman*tman - 
     &        4*BvgR*m3*mq1*PvR*sman*tman - Bvm1R*m4**2*mq1*PvR*sman*tman - Bvm2R*m4**2*mq1*PvR*sman*tman - 4*BsgL*m3**2*PsL*uman + 4*BsgL*m4**2*PsL*uman + 
     &        4*BvgL*m3**2*mq1*PsL*uman - 4*BvgL*m4**2*mq1*PsL*uman - 2*Bvm1L*m3*m4**2*mq1*PsL*uman - 2*Bvm2L*m3*m4**2*mq1*PsL*uman - 4*BsgL*PsL*sman*uman + 
     &        4*BvgR*m3*PsL*sman*uman + 4*BvgL*mq1*PsL*sman*uman + 4*BsgR*m3*PvR*sman*uman - 4*BvgL*m3**2*PvR*sman*uman + 4*BvgL*m4**2*PvR*sman*uman + 
     &        2*Bvm1L*m3*m4**2*PvR*sman*uman + 2*Bvm2L*m3*m4**2*PvR*sman*uman - 4*BvgR*m3*mq1*PvR*sman*uman - 4*BvgL*PvR*sman**2*uman + 
     &        m4**2*(-(BsgL*m1**2*PsL*sman) + 2*BvgR*m1**2*m3*PsL*sman + BvgR*m3**3*PsL*sman - BsgL*m4**2*PsL*sman + Bvm1R*m1**2*m4**2*PsL*sman + 
     &           Bvm2R*m1**2*m4**2*PsL*sman + BvgR*m3*m4**2*PsL*sman + Bvm1R*m3**2*m4**2*PsL*sman + Bvm2R*m3**2*m4**2*PsL*sman + BvgL*m1**2*mq1*PsL*sman + 
     &           BvgL*m4**2*mq1*PsL*sman + 2*BsgR*m1**2*m3*PvR*sman - 2*BvgL*m1**2*m3**2*PvR*sman + BsgR*m3**3*PvR*sman - BvgL*m3**4*PvR*sman + 
     &           Bsm1R*m1**2*m4**2*PvR*sman + Bsm2R*m1**2*m4**2*PvR*sman + 2*BvgL*m1**2*m4**2*PvR*sman + BsgR*m3*m4**2*PvR*sman + Bsm1R*m3**2*m4**2*PvR*sman + 
     &           Bsm2R*m3**2*m4**2*PvR*sman + BvgL*m4**4*PvR*sman - 2*BvgR*m1**2*m3*mq1*PvR*sman - BvgR*m3**3*mq1*PvR*sman - Bvm1R*m1**2*m4**2*mq1*PvR*sman - 
     &           Bvm2R*m1**2*m4**2*mq1*PvR*sman - BvgR*m3*m4**2*mq1*PvR*sman - Bvm1R*m3**2*m4**2*mq1*PvR*sman - Bvm2R*m3**2*m4**2*mq1*PvR*sman - 
     &           BvgL*m1**2*PvR*sman**2 - BvgL*m4**2*PvR*sman**2 + BsgL*m3**2*PsL*tman - BsgL*m4**2*PsL*tman - BvgL*m3**2*mq1*PsL*tman + 
     &           BvgL*m4**2*mq1*PsL*tman - BvgR*m3*PsL*sman*tman - Bvm1R*m4**2*PsL*sman*tman - Bvm2R*m4**2*PsL*sman*tman - BsgR*m3*PvR*sman*tman + 
     &           BvgL*m3**2*PvR*sman*tman - Bsm1R*m4**2*PvR*sman*tman - Bsm2R*m4**2*PvR*sman*tman - BvgL*m4**2*PvR*sman*tman + BvgR*m3*mq1*PvR*sman*tman + 
     &           Bvm1R*m4**2*mq1*PvR*sman*tman + Bvm2R*m4**2*mq1*PvR*sman*tman - 
     &           2*Bvm1L*m3*m4**2*(mq1*PsL - PvR*sman)*(2*m1**2 + m3**2 + m4**2 - tman - uman) - 
     &           2*Bvm2L*m3*m4**2*(mq1*PsL - PvR*sman)*(2*m1**2 + m3**2 + m4**2 - tman - uman) + 
     &           (-(m3*(BvgR*PsL + BsgR*PvR - BvgR*mq1*PvR)*sman) + BsgL*PsL*(m3**2 - m4**2 + sman) + BvgL*(m3**2 - m4**2 + sman)*(-(mq1*PsL) + PvR*sman))*uman
     &           )*xv + 2*Bsm1L*m3*m4**2*PsL*(2*m1**2 + m3**2 + m4**2 - tman - uman)*(-1 + m4**2*xv) + 
     &        2*Bsm2L*m3*m4**2*PsL*(2*m1**2 + m3**2 + m4**2 - tman - uman)*(-1 + m4**2*xv)) + 
     &     BLc*((Bvm1L*m3**5*mq1*PsL*tman - Bvm1L*m3**3*m4**2*mq1*PsL*tman - Bvm2L*m3**3*m4**2*mq1*PsL*tman - 2*Bvm1L*m3*m4**4*mq1*PsL*tman - 
     &           Bvm2L*m3*m4**4*mq1*PsL*tman - Bvm1R*m3**4*PsL*sman*tman + Bvm2R*m3**2*m4**2*PsL*sman*tman + Bvm1R*m4**4*PsL*sman*tman + 
     &           Bvm2R*m4**4*PsL*sman*tman - 2*Bvm1L*m3**3*mq1*PsL*sman*tman + Bvm1L*m3*m4**2*mq1*PsL*sman*tman + Bvm2L*m3*m4**2*mq1*PsL*sman*tman - 
     &           Bsm1R*m3**4*PvR*sman*tman - Bvm1L*m3**5*PvR*sman*tman + Bsm2R*m3**2*m4**2*PvR*sman*tman + Bvm1L*m3**3*m4**2*PvR*sman*tman + 
     &           Bvm2L*m3**3*m4**2*PvR*sman*tman + Bsm1R*m4**4*PvR*sman*tman + Bsm2R*m4**4*PvR*sman*tman + 2*Bvm1L*m3*m4**4*PvR*sman*tman + 
     &           Bvm2L*m3*m4**4*PvR*sman*tman + Bvm1R*m3**4*mq1*PvR*sman*tman - Bvm2R*m3**2*m4**2*mq1*PvR*sman*tman - Bvm1R*m4**4*mq1*PvR*sman*tman - 
     &           Bvm2R*m4**4*mq1*PvR*sman*tman + 2*Bvm1R*m3**2*PsL*sman**2*tman - Bvm2R*m4**2*PsL*sman**2*tman + Bvm1L*m3*mq1*PsL*sman**2*tman + 
     &           2*Bsm1R*m3**2*PvR*sman**2*tman + 2*Bvm1L*m3**3*PvR*sman**2*tman - Bsm2R*m4**2*PvR*sman**2*tman - Bvm1L*m3*m4**2*PvR*sman**2*tman - 
     &           Bvm2L*m3*m4**2*PvR*sman**2*tman - 2*Bvm1R*m3**2*mq1*PvR*sman**2*tman + Bvm2R*m4**2*mq1*PvR*sman**2*tman - Bvm1R*PsL*sman**3*tman - 
     &           Bsm1R*PvR*sman**3*tman - Bvm1L*m3*PvR*sman**3*tman + Bvm1R*mq1*PvR*sman**3*tman + 
     &           m3*(-(Bvm2L*m4**2*(m3**2 + m4**2)*(mq1*PsL - PvR*sman)) + 
     &              m3*sman*(-(Bsm1R*m3**2*PvR) + m4**2*(Bsm2R*PvR + Bvm2R*(PsL - mq1*PvR)) - Bvm1R*(PsL - mq1*PvR)*(m3**2 - sman) + Bsm1R*PvR*sman) + 
     &              Bvm1L*(mq1*PsL - PvR*sman)*(m3**4 - 2*m4**4 + m4**2*sman - m3**2*(m4**2 + sman)))*uman + 
     &           Bsm2L*m3*m4**2*PsL*(-(sman*tman) + m3**2*(tman + uman) + m4**2*(tman + uman)) + 
     &           Bsm1L*m3*PsL*(-(sman**2*tman) - m3**4*(tman + uman) + 2*m4**4*(tman + uman) - m4**2*sman*(tman + uman) + 
     &              m3**2*(m4**2*(tman + uman) + sman*(2*tman + uman))))*xv + 
     &        BvgL*(4*(-m3**4 + m4**4)*mq1*PsL + (-(m4**6*mq1*PsL) + sman**2*(mq1*PsL - PvR*sman)*tman + m4**2*sman*(mq1*PsL - PvR*sman)*uman - 
     &              m4**4*(mq1*PsL - PvR*sman)*(tman + uman) + m3**2*sman*(-(mq1*PsL) + PvR*sman)*(2*tman + uman) + 
     &              m3**4*(m4**2*mq1*PsL + (mq1*PsL - PvR*sman)*(tman + uman)))*xv + 2*m1**2*(m3 - m4)*(m3 + m4)*mq1*PsL*(-4 + m4**2*xv)) + 
     &        BsgL*PsL*(-4*m4**4 + m3**2*sman*(2*tman + uman)*xv + (m4**6 - sman**2*tman - m4**2*sman*uman + m4**4*(tman + uman))*xv - 
     &           2*m1**2*(m3 - m4)*(m3 + m4)*(-4 + m4**2*xv) + m3**4*(4 - (m4**2 + tman + uman)*xv)))) + 
     &  A1Lc*A1vL*(BRc*(-2*Bsm1R*m1**2*m3**3*PsR + 2*Bsm2R*m1**2*m3**3*PsR - Bsm1R*m3**5*PsR + Bsm2R*m3**5*PsR + 6*Bsm1R*m1**2*m3*m4**2*PsR + 
     &        2*Bsm2R*m1**2*m3*m4**2*PsR + 2*Bsm1R*m3**3*m4**2*PsR + 2*Bsm2R*m3**3*m4**2*PsR + 3*Bsm1R*m3*m4**4*PsR + Bsm2R*m3*m4**4*PsR - 
     &        2*BsgR*(m3 - m4)*(m3 + m4)*(2*m1**2 + m3**2 + m4**2)*PsR + 4*BvgR*m1**2*m3**2*mq1*PsR + 2*Bvm1R*m1**2*m3**3*mq1*PsR - 
     &        2*Bvm2R*m1**2*m3**3*mq1*PsR + 2*BvgR*m3**4*mq1*PsR + Bvm1R*m3**5*mq1*PsR - Bvm2R*m3**5*mq1*PsR - 4*BvgR*m1**2*m4**2*mq1*PsR - 
     &        6*Bvm1R*m1**2*m3*m4**2*mq1*PsR - 2*Bvm2R*m1**2*m3*m4**2*mq1*PsR - 2*Bvm1R*m3**3*m4**2*mq1*PsR - 2*Bvm2R*m3**3*m4**2*mq1*PsR - 
     &        2*BvgR*m4**4*mq1*PsR - 3*Bvm1R*m3*m4**4*mq1*PsR - Bvm2R*m3*m4**4*mq1*PsR - 2*BsgR*(m1**2 + m4**2)*PsR*sman + 
     &        (-(Bsm1R*m3*(3*m1**2 + m3**2 + 2*m4**2)*PsR) + (-(Bsm2R*m3*(m1**2 + m3**2)) - Bvm1L*(m3**4 + m1**2*(2*m3**2 - m4**2)) + 
     &              Bvm2L*(m3**4 + 2*m3**2*m4**2 + m1**2*(2*m3**2 + m4**2)) + 
     &              (2*BvgR*(m1**2 + m4**2) + m3*(Bvm2R*(m1**2 + m3**2) + Bvm1R*(3*m1**2 + m3**2 + 2*m4**2)))*mq1)*PsR + 
     &           (-((2*m1**2 + m3**2 + m4**2)*(m3**2*(2*BvgR + (Bvm1R - Bvm2R)*m3) - (2*BvgR + (3*Bvm1R + Bvm2R)*m3)*m4**2)) - 
     &              Bsm1L*(m3**4 + m1**2*(2*m3**2 - m4**2)) + Bsm2L*(m3**4 + 2*m3**2*m4**2 + m1**2*(2*m3**2 + m4**2)) + 
     &              ((Bvm1L - Bvm2L)*m3**2*(2*m1**2 + m3**2) - ((Bvm1L + Bvm2L)*m1**2 + 2*Bvm2L*m3**2)*m4**2)*mq1)*PvL)*sman) - 
     &     BLc*mq1c*((-2*Bvm1R*m3**2*(2*m1**2 + m3**2) + ((Bvm1R + Bvm2R)*m1**2 - 2*Bvm1R*m3**2)*m4**2 + (Bvm1R + Bvm2R)*m4**4)*mq1*PsR - 
     &        Bsm2R*(m1**2 + m4**2)*PsR*(m4**2 - sman) + ((m1**2 + m4**2)*(Bvm2L*m3 - (Bvm1R + Bvm2R)*mq1)*PsR + 
     &           (Bsm2L*m1**2*m3 + 4*Bvm1R*m1**2*m3**2 + 2*Bvm1R*m3**4 + 2*BsgL*(m1**2 + m3**2) - Bvm1R*m1**2*m4**2 - Bvm2R*m1**2*m4**2 + Bsm2L*m3*m4**2 + 
     &              2*Bvm1R*m3**2*m4**2 - Bvm1R*m4**4 - Bvm2R*m4**4 + Bsm1L*m3*(3*m1**2 + 2*m3**2 + m4**2) - Bvm2L*m3*(m1**2 + m4**2)*mq1)*PvL + 
     &           2*BvgL*(m1**2 + m3**2)*(PsR - mq1*PvL) + Bvm1L*m3*(3*m1**2 + 2*m3**2 + m4**2)*(PsR - mq1*PvL))*sman + 
     &        Bsm1R*PsR*(2*m3**4 + 2*m3**2*m4**2 - m4**4 + m4**2*sman + m1**2*(4*m3**2 - m4**2 + sman))) - 
     &     BRc*(Bsm1L*m1**2*PvL*sman**2 + Bsm2L*m1**2*PvL*sman**2 + 2*BvgR*m1**2*PvL*sman**2 + 3*Bvm1R*m1**2*m3*PvL*sman**2 + Bvm2R*m1**2*m3*PvL*sman**2 + 
     &        Bsm1L*m3**2*PvL*sman**2 + Bsm2L*m3**2*PvL*sman**2 + Bvm1R*m3**3*PvL*sman**2 + Bvm2R*m3**3*PvL*sman**2 + 2*BvgR*m4**2*PvL*sman**2 + 
     &        2*Bvm1R*m3*m4**2*PvL*sman**2 - 2*BsgR*m3**2*PsR*tman - Bsm1R*m3**3*PsR*tman + Bsm2R*m3**3*PsR*tman + 2*BsgR*m4**2*PsR*tman + 
     &        3*Bsm1R*m3*m4**2*PsR*tman + Bsm2R*m3*m4**2*PsR*tman + 2*BvgR*m3**2*mq1*PsR*tman + Bvm1R*m3**3*mq1*PsR*tman - Bvm2R*m3**3*mq1*PsR*tman - 
     &        2*BvgR*m4**2*mq1*PsR*tman - 3*Bvm1R*m3*m4**2*mq1*PsR*tman - Bvm2R*m3*m4**2*mq1*PsR*tman - Bsm1R*m3*PsR*sman*tman - Bsm2R*m3*PsR*sman*tman + 
     &        Bvm1R*m3*mq1*PsR*sman*tman + Bvm2R*m3*mq1*PsR*sman*tman - Bsm1L*m3**2*PvL*sman*tman + Bsm2L*m3**2*PvL*sman*tman - 2*BvgR*m3**2*PvL*sman*tman - 
     &        Bvm1R*m3**3*PvL*sman*tman + Bvm2R*m3**3*PvL*sman*tman + Bsm1L*m4**2*PvL*sman*tman + Bsm2L*m4**2*PvL*sman*tman + 2*BvgR*m4**2*PvL*sman*tman + 
     &        3*Bvm1R*m3*m4**2*PvL*sman*tman + Bvm2R*m3*m4**2*PvL*sman*tman - Bsm1L*PvL*sman**2*tman - Bsm2L*PvL*sman**2*tman - Bvm1R*m3*PvL*sman**2*tman - 
     &        Bvm2R*m3*PvL*sman**2*tman + Bvm1L*(PsR - mq1*PvL)*sman*(m1**2*sman + (m4**2 - sman)*tman + m3**2*(sman - tman - uman)) + 
     &        ((Bsm2R*m3*(m3**2 + m4**2) + (m3**2*(2*BvgR + (Bvm1R - Bvm2R)*m3) - (2*BvgR + (3*Bvm1R + Bvm2R)*m3)*m4**2)*mq1)*PsR + 
     &           (2*(BvgR + Bvm1R*m3)*mq1*PsR + (m3**2*(-Bsm1L + Bsm2L - 2*BvgR - Bvm1R*m3 + Bvm2R*m3) + (2*BvgR + (3*Bvm1R + Bvm2R)*m3)*m4**2)*PvL)*sman - 
     &           2*(BvgR + Bvm1R*m3)*PvL*sman**2 - 2*BsgR*PsR*(m3**2 - m4**2 + sman) - Bsm1R*m3*PsR*(m3**2 - 3*m4**2 + 2*sman))*uman + 
     &        Bvm2L*(PsR - mq1*PvL)*sman*(m1**2*sman + (m4**2 - sman)*tman + m3**2*(sman + tman + uman))) - 
     &     BLc*mq1c*(Bvm2R*(m1**2 + m4**2)*PvL*sman**2 - 2*(Bsm1R*m3**2*PsR + ((BvgL + Bvm1L*m3)*PsR + (BsgL + Bsm1L*m3 - (BvgL + Bvm1L*m3)*mq1)*PvL)*sman)*
     &         tman + Bvm2R*(m4**2 - sman)*(-(mq1*PsR) + PvL*sman)*uman - 
     &        ((2*Bsm1R*m3**2 - (Bsm1R + Bsm2R)*m4**2)*PsR + ((Bsm1R + Bsm2R + (Bvm1L + Bvm2L)*m3)*PsR + m3*(Bsm1L + Bsm2L - (Bvm1L + Bvm2L)*mq1)*PvL)*sman)*
     &         uman + Bvm1R*(m1**2*PvL*sman**2 + m4**2*(-(mq1*PsR*uman) + PvL*sman*(sman + uman)) + (mq1*PsR - PvL*sman)*(sman*uman + 2*m3**2*(tman + uman)))) +
     &      (BLc*mq1c*(m4**2*(-(Bsm2R*m4**2*(m1**2 + m4**2)) + Bsm1R*(m3**4 + 2*m3**2*m4**2 + m1**2*(2*m3**2 + m4**2)) - 
     &              (-(Bvm2R*m4**2*(m1**2 + m4**2)) + Bvm1R*(m3**4 + 2*m3**2*m4**2 + m1**2*(2*m3**2 + m4**2)))*mq1)*PsR - 
     &           (Bsm1R*(m1**2*m3**2 + 2*(m1**2 + m3**2)*m4**2 + m4**4)*PsR + 
     &              ((m1**2 + m4**2)*(m3**2*(BvgL + Bvm1L*m3) - (Bsm2R - BvgL + Bvm2L*m3)*m4**2) - 
     &                 (-(Bvm2R*m4**2*(m1**2 + m4**2)) + Bvm1R*(m1**2*m3**2 + 2*(m1**2 + m3**2)*m4**2 + m4**4))*mq1)*PsR + 
     &              (-(m3*(Bsm2L*m1**2 + Bvm1R*m3*(2*m1**2 + m3**2))*m4**2) - ((Bvm1R - Bvm2R)*m1**2 + Bsm2L*m3 + 2*Bvm1R*m3**2)*m4**4 + Bvm2R*m4**6 + 
     &                 Bsm1L*m3**3*(m1**2 + m4**2) + BsgL*(m1**2 + m4**2)*(m3**2 + m4**2) - 
     &                 (m1**2 + m4**2)*(m3**2*(BvgL + Bvm1L*m3) + (BvgL - Bvm2L*m3)*m4**2)*mq1)*PvL)*sman + 
     &           ((m1**2 + m4**2)*(Bsm1R + BvgL + Bvm1L*m3 - Bvm1R*mq1)*PsR + 
     &              (-(Bvm1R*m1**2*m3**2) - 2*Bvm1R*m1**2*m4**2 + Bvm2R*m1**2*m4**2 - 2*Bvm1R*m3**2*m4**2 - Bvm1R*m4**4 + Bvm2R*m4**4 + BsgL*(m1**2 + m4**2) + 
     &                 Bsm1L*m3*(m1**2 + m4**2) - (BvgL + Bvm1L*m3)*(m1**2 + m4**2)*mq1)*PvL)*sman**2 + Bvm1R*(m1**2 + m4**2)*PvL*sman**3) + 
     &        BRc*(2*Bsm1R*m1**2*m3**5*PsR + Bsm1R*m3**7*PsR - 2*Bsm1R*m1**2*m3**3*m4**2*PsR - 2*Bsm2R*m1**2*m3**3*m4**2*PsR - Bsm2R*m3**5*m4**2*PsR - 
     &           4*Bsm1R*m1**2*m3*m4**4*PsR - 2*Bsm2R*m1**2*m3*m4**4*PsR - 3*Bsm1R*m3**3*m4**4*PsR - 2*Bsm2R*m3**3*m4**4*PsR - 2*Bsm1R*m3*m4**6*PsR - 
     &           Bsm2R*m3*m4**6*PsR - 2*BvgR*m1**2*m3**4*mq1*PsR - 2*Bvm1R*m1**2*m3**5*mq1*PsR - BvgR*m3**6*mq1*PsR - Bvm1R*m3**7*mq1*PsR + 
     &           2*Bvm1R*m1**2*m3**3*m4**2*mq1*PsR + 2*Bvm2R*m1**2*m3**3*m4**2*mq1*PsR - BvgR*m3**4*m4**2*mq1*PsR + Bvm2R*m3**5*m4**2*mq1*PsR + 
     &           2*BvgR*m1**2*m4**4*mq1*PsR + 4*Bvm1R*m1**2*m3*m4**4*mq1*PsR + 2*Bvm2R*m1**2*m3*m4**4*mq1*PsR + BvgR*m3**2*m4**4*mq1*PsR + 
     &           3*Bvm1R*m3**3*m4**4*mq1*PsR + 2*Bvm2R*m3**3*m4**4*mq1*PsR + BvgR*m4**6*mq1*PsR + 2*Bvm1R*m3*m4**6*mq1*PsR + Bvm2R*m3*m4**6*mq1*PsR - 
     &           (-(Bsm2R*m1**2*m3*m4**2*PsR) + 2*Bvm2L*m1**2*m3**2*m4**2*PsR - Bsm2R*m3**3*m4**2*PsR + Bvm2L*m3**4*m4**2*PsR + Bvm2L*m1**2*m4**4*PsR + 
     &              2*Bvm2L*m3**2*m4**4*PsR + Bsm1R*m3*(2*m3**4 - m4**4 + m1**2*(3*m3**2 - 2*m4**2))*PsR - 3*BvgR*m1**2*m3**2*mq1*PsR - 
     &              3*Bvm1R*m1**2*m3**3*mq1*PsR - 2*BvgR*m3**4*mq1*PsR - 2*Bvm1R*m3**5*mq1*PsR + BvgR*m1**2*m4**2*mq1*PsR + 2*Bvm1R*m1**2*m3*m4**2*mq1*PsR + 
     &              Bvm2R*m1**2*m3*m4**2*mq1*PsR - BvgR*m3**2*m4**2*mq1*PsR + Bvm2R*m3**3*m4**2*mq1*PsR + BvgR*m4**4*mq1*PsR + Bvm1R*m3*m4**4*mq1*PsR + 
     &              (-2*Bvm1R*m1**2*m3**5 - Bvm1R*m3**7 + 2*Bsm2L*m1**2*m3**2*m4**2 + 2*Bvm1R*m1**2*m3**3*m4**2 + 2*Bvm2R*m1**2*m3**3*m4**2 + 
     &                 Bsm2L*m3**4*m4**2 + Bvm2R*m3**5*m4**2 + Bsm2L*m1**2*m4**4 + 4*Bvm1R*m1**2*m3*m4**4 + 2*Bvm2R*m1**2*m3*m4**4 + 2*Bsm2L*m3**2*m4**4 + 
     &                 3*Bvm1R*m3**3*m4**4 + 2*Bvm2R*m3**3*m4**4 + 2*Bvm1R*m3*m4**6 + Bvm2R*m3*m4**6 - BvgR*(2*m1**2 + m3**2 + m4**2)*(m3**4 - m4**4) + 
     &                 Bsm1L*(-(m3**4*(2*m1**2 + m3**2)) - m3**4*m4**2 + (m1**2 + m3**2)*m4**4) - 
     &                 Bvm2L*m4**2*(m3**4 + 2*m3**2*m4**2 + m1**2*(2*m3**2 + m4**2))*mq1)*PvL - 
     &              Bvm1L*(2*m1**2*m3**4 + m3**6 + m3**4*m4**2 - (m1**2 + m3**2)*m4**4)*(PsR - mq1*PvL))*sman + 
     &           (Bsm1R*m3*(m1**2 + m3**2)*PsR - (m1**2 + m3**2)*(-(Bvm2L*m4**2) + (BvgR + Bvm1R*m3)*mq1)*PsR + 
     &              (-3*Bvm1R*m1**2*m3**3 - 2*Bvm1R*m3**5 + Bsm2L*m1**2*m4**2 + 2*Bvm1R*m1**2*m3*m4**2 + Bvm2R*m1**2*m3*m4**2 + Bsm2L*m3**2*m4**2 + 
     &                 Bvm2R*m3**3*m4**2 + Bvm1R*m3*m4**4 - Bsm1L*m3**2*(3*m1**2 + 2*m3**2 + m4**2) + 
     &                 BvgR*(-2*m3**4 - m3**2*m4**2 + m4**4 + m1**2*(-3*m3**2 + m4**2)) - Bvm2L*(m1**2 + m3**2)*m4**2*mq1)*PvL - 
     &              Bvm1L*m3**2*(3*m1**2 + 2*m3**2 + m4**2)*(PsR - mq1*PvL))*sman**2 + 
     &           (m1**2 + m3**2)*(Bvm1L*PsR + (Bsm1L + BvgR + Bvm1R*m3 - Bvm1L*mq1)*PvL)*sman**3 + 
     &           BsgR*PsR*(m3**2 + m4**2 - sman)*((m3 - m4)*(m3 + m4)*(2*m1**2 + m3**2 + m4**2) - (m1**2 + m3**2)*sman)))*xv + 
     &     BLc*mq1c*(-(((m4**4*(-Bsm2R + Bvm2R*mq1)*PsR - (-(m4**2*(Bsm2R + Bvm2L*m3 - Bvm2R*mq1)*PsR) + 
     &                   (Bsm1L*m3**3 + BsgL*(m3**2 + m4**2) + m4**2*(-(Bsm2L*m3) + Bvm2R*m4**2 + Bvm2L*m3*mq1))*PvL + Bvm1L*m3**3*(PsR - mq1*PvL) + 
     &                   BvgL*(m3**2 + m4**2)*(PsR - mq1*PvL))*sman + 
     &                ((BvgL + Bvm1L*m3)*PsR + (BsgL + Bsm1L*m3 + Bvm2R*m4**2 - (BvgL + Bvm1L*m3)*mq1)*PvL)*sman**2)*uman + 
     &             Bsm1R*PsR*(m4**2*(m3**2 + m4**2 - sman)*tman + (m3**2 - sman)*(m4**2 - sman)*uman) + 
     &             Bvm1R*(-(mq1*PsR) + PvL*sman)*(m4**2*(m3**2 + m4**2 - sman)*tman + (m3**2 - sman)*(m4**2 - sman)*uman))*xv) + 
     &        BsgR*m3*(2*m1**2 + m3**2 + m4**2)*PsR*(-4 + m4**2*xv) - BvgR*m3*(2*m1**2 + m3**2 + m4**2)*mq1*PsR*(-4 + m4**2*xv)) + 
     &     BLc*mq1c*(-(Bvm1R*m1**2*m4**2*mq1*PsR) - Bvm2R*m1**2*m4**2*mq1*PsR - Bvm1R*m3**2*m4**2*mq1*PsR - Bvm2R*m3**2*m4**2*mq1*PsR + 
     &        4*BvgL*m1**2*PsR*sman + 4*BvgL*m3**2*PsR*sman + 4*BsgL*m1**2*PvL*sman - 8*BvgR*m1**2*m3*PvL*sman + 4*BsgL*m3**2*PvL*sman - 
     &        4*BvgR*m3**3*PvL*sman + Bvm1R*m1**2*m4**2*PvL*sman + Bvm2R*m1**2*m4**2*PvL*sman - 4*BvgR*m3*m4**2*PvL*sman + Bvm1R*m3**2*m4**2*PvL*sman + 
     &        Bvm2R*m3**2*m4**2*PvL*sman - 4*BvgL*m1**2*mq1*PvL*sman - 4*BvgL*m3**2*mq1*PvL*sman + 4*BsgR*m3*PsR*tman - 4*BvgR*m3*mq1*PsR*tman + 
     &        Bvm1R*m4**2*mq1*PsR*tman + Bvm2R*m4**2*mq1*PsR*tman - 4*BvgL*PsR*sman*tman - 4*BsgL*PvL*sman*tman + 4*BvgR*m3*PvL*sman*tman - 
     &        Bvm1R*m4**2*PvL*sman*tman - Bvm2R*m4**2*PvL*sman*tman + 4*BvgL*mq1*PvL*sman*tman + 4*BsgR*m3*PsR*uman - 4*BvgR*m3*mq1*PsR*uman + 
     &        4*BvgR*m3*PvL*sman*uman + m4**2*(-(BvgL*m1**2*PsR*sman) - BvgL*m3**2*PsR*sman - BsgL*m1**2*PvL*sman + 2*BvgR*m1**2*m3*PvL*sman - 
     &           BsgL*m3**2*PvL*sman + BvgR*m3**3*PvL*sman + BvgR*m3*m4**2*PvL*sman + BvgL*m1**2*mq1*PvL*sman + BvgL*m3**2*mq1*PvL*sman + 
     &           Bvm1R*m4**2*(mq1*PsR - PvL*sman)*(m1**2 + m3**2 - tman) + Bvm2R*m4**2*(mq1*PsR - PvL*sman)*(m1**2 + m3**2 - tman) - BsgR*m3*PsR*tman + 
     &           BvgR*m3*mq1*PsR*tman + BvgL*PsR*sman*tman + BsgL*PvL*sman*tman - BvgR*m3*PvL*sman*tman - BvgL*mq1*PvL*sman*tman - 
     &           m3*(BsgR*PsR - BvgR*mq1*PsR + BvgR*PvL*sman)*uman)*xv - Bsm1R*m4**2*PsR*(m1**2 + m3**2 - tman)*(-1 + m4**2*xv) - 
     &        Bsm2R*m4**2*PsR*(m1**2 + m3**2 - tman)*(-1 + m4**2*xv)) + 
     &     BRc*(4*Bvm1R*m1**2*m3*m4**2*mq1*PsR + 4*Bvm2R*m1**2*m3*m4**2*mq1*PsR + 2*Bvm1R*m3**3*m4**2*mq1*PsR + 2*Bvm2R*m3**3*m4**2*mq1*PsR + 
     &        2*Bvm1R*m3*m4**4*mq1*PsR + 2*Bvm2R*m3*m4**4*mq1*PsR + 4*BsgR*m1**2*PsR*sman - 8*BvgL*m1**2*m3*PsR*sman - 4*BvgL*m3**3*PsR*sman + 
     &        4*BsgR*m4**2*PsR*sman - Bvm1L*m1**2*m4**2*PsR*sman - Bvm2L*m1**2*m4**2*PsR*sman - 4*BvgL*m3*m4**2*PsR*sman - Bvm1L*m3**2*m4**2*PsR*sman - 
     &        Bvm2L*m3**2*m4**2*PsR*sman - 4*BvgR*m1**2*mq1*PsR*sman - 4*BvgR*m4**2*mq1*PsR*sman - 8*BsgL*m1**2*m3*PvL*sman + 8*BvgR*m1**2*m3**2*PvL*sman - 
     &        4*BsgL*m3**3*PvL*sman + 4*BvgR*m3**4*PvL*sman - Bsm1L*m1**2*m4**2*PvL*sman - Bsm2L*m1**2*m4**2*PvL*sman - 8*BvgR*m1**2*m4**2*PvL*sman - 
     &        4*BsgL*m3*m4**2*PvL*sman - 4*Bvm1R*m1**2*m3*m4**2*PvL*sman - 4*Bvm2R*m1**2*m3*m4**2*PvL*sman - Bsm1L*m3**2*m4**2*PvL*sman - 
     &        Bsm2L*m3**2*m4**2*PvL*sman - 2*Bvm1R*m3**3*m4**2*PvL*sman - 2*Bvm2R*m3**3*m4**2*PvL*sman - 4*BvgR*m4**4*PvL*sman - 2*Bvm1R*m3*m4**4*PvL*sman - 
     &        2*Bvm2R*m3*m4**4*PvL*sman + 8*BvgL*m1**2*m3*mq1*PvL*sman + 4*BvgL*m3**3*mq1*PvL*sman + Bvm1L*m1**2*m4**2*mq1*PvL*sman + 
     &        Bvm2L*m1**2*m4**2*mq1*PvL*sman + 4*BvgL*m3*m4**2*mq1*PvL*sman + Bvm1L*m3**2*m4**2*mq1*PvL*sman + Bvm2L*m3**2*m4**2*mq1*PvL*sman + 
     &        4*BvgR*m1**2*PvL*sman**2 + 4*BvgR*m4**2*PvL*sman**2 - 4*BsgR*m3**2*PsR*tman + 4*BsgR*m4**2*PsR*tman + 4*BvgR*m3**2*mq1*PsR*tman - 
     &        4*BvgR*m4**2*mq1*PsR*tman - 2*Bvm1R*m3*m4**2*mq1*PsR*tman - 2*Bvm2R*m3*m4**2*mq1*PsR*tman + 4*BvgL*m3*PsR*sman*tman + 
     &        Bvm1L*m4**2*PsR*sman*tman + Bvm2L*m4**2*PsR*sman*tman + 4*BsgL*m3*PvL*sman*tman - 4*BvgR*m3**2*PvL*sman*tman + Bsm1L*m4**2*PvL*sman*tman + 
     &        Bsm2L*m4**2*PvL*sman*tman + 4*BvgR*m4**2*PvL*sman*tman + 2*Bvm1R*m3*m4**2*PvL*sman*tman + 2*Bvm2R*m3*m4**2*PvL*sman*tman - 
     &        4*BvgL*m3*mq1*PvL*sman*tman - Bvm1L*m4**2*mq1*PvL*sman*tman - Bvm2L*m4**2*mq1*PvL*sman*tman - 4*BsgR*m3**2*PsR*uman + 4*BsgR*m4**2*PsR*uman + 
     &        4*BvgR*m3**2*mq1*PsR*uman - 4*BvgR*m4**2*mq1*PsR*uman - 2*Bvm1R*m3*m4**2*mq1*PsR*uman - 2*Bvm2R*m3*m4**2*mq1*PsR*uman - 4*BsgR*PsR*sman*uman + 
     &        4*BvgL*m3*PsR*sman*uman + 4*BvgR*mq1*PsR*sman*uman + 4*BsgL*m3*PvL*sman*uman - 4*BvgR*m3**2*PvL*sman*uman + 4*BvgR*m4**2*PvL*sman*uman + 
     &        2*Bvm1R*m3*m4**2*PvL*sman*uman + 2*Bvm2R*m3*m4**2*PvL*sman*uman - 4*BvgL*m3*mq1*PvL*sman*uman - 4*BvgR*PvL*sman**2*uman + 
     &        m4**2*(-(BsgR*m1**2*PsR*sman) + 2*BvgL*m1**2*m3*PsR*sman + BvgL*m3**3*PsR*sman - BsgR*m4**2*PsR*sman + Bvm1L*m1**2*m4**2*PsR*sman + 
     &           Bvm2L*m1**2*m4**2*PsR*sman + BvgL*m3*m4**2*PsR*sman + Bvm1L*m3**2*m4**2*PsR*sman + Bvm2L*m3**2*m4**2*PsR*sman + BvgR*m1**2*mq1*PsR*sman + 
     &           BvgR*m4**2*mq1*PsR*sman + 2*BsgL*m1**2*m3*PvL*sman - 2*BvgR*m1**2*m3**2*PvL*sman + BsgL*m3**3*PvL*sman - BvgR*m3**4*PvL*sman + 
     &           Bsm1L*m1**2*m4**2*PvL*sman + Bsm2L*m1**2*m4**2*PvL*sman + 2*BvgR*m1**2*m4**2*PvL*sman + BsgL*m3*m4**2*PvL*sman + Bsm1L*m3**2*m4**2*PvL*sman + 
     &           Bsm2L*m3**2*m4**2*PvL*sman + BvgR*m4**4*PvL*sman - 2*BvgL*m1**2*m3*mq1*PvL*sman - BvgL*m3**3*mq1*PvL*sman - Bvm1L*m1**2*m4**2*mq1*PvL*sman - 
     &           Bvm2L*m1**2*m4**2*mq1*PvL*sman - BvgL*m3*m4**2*mq1*PvL*sman - Bvm1L*m3**2*m4**2*mq1*PvL*sman - Bvm2L*m3**2*m4**2*mq1*PvL*sman - 
     &           BvgR*m1**2*PvL*sman**2 - BvgR*m4**2*PvL*sman**2 + BsgR*m3**2*PsR*tman - BsgR*m4**2*PsR*tman - BvgR*m3**2*mq1*PsR*tman + 
     &           BvgR*m4**2*mq1*PsR*tman - BvgL*m3*PsR*sman*tman - Bvm1L*m4**2*PsR*sman*tman - Bvm2L*m4**2*PsR*sman*tman - BsgL*m3*PvL*sman*tman + 
     &           BvgR*m3**2*PvL*sman*tman - Bsm1L*m4**2*PvL*sman*tman - Bsm2L*m4**2*PvL*sman*tman - BvgR*m4**2*PvL*sman*tman + BvgL*m3*mq1*PvL*sman*tman + 
     &           Bvm1L*m4**2*mq1*PvL*sman*tman + Bvm2L*m4**2*mq1*PvL*sman*tman - 
     &           2*Bvm1R*m3*m4**2*(mq1*PsR - PvL*sman)*(2*m1**2 + m3**2 + m4**2 - tman - uman) - 
     &           2*Bvm2R*m3*m4**2*(mq1*PsR - PvL*sman)*(2*m1**2 + m3**2 + m4**2 - tman - uman) + 
     &           (-(m3*(BvgL*PsR + BsgL*PvL - BvgL*mq1*PvL)*sman) + BsgR*PsR*(m3**2 - m4**2 + sman) + BvgR*(m3**2 - m4**2 + sman)*(-(mq1*PsR) + PvL*sman))*uman
     &           )*xv + 2*Bsm1R*m3*m4**2*PsR*(2*m1**2 + m3**2 + m4**2 - tman - uman)*(-1 + m4**2*xv) + 
     &        2*Bsm2R*m3*m4**2*PsR*(2*m1**2 + m3**2 + m4**2 - tman - uman)*(-1 + m4**2*xv)) + 
     &     BRc*((Bvm1R*m3**5*mq1*PsR*tman - Bvm1R*m3**3*m4**2*mq1*PsR*tman - Bvm2R*m3**3*m4**2*mq1*PsR*tman - 2*Bvm1R*m3*m4**4*mq1*PsR*tman - 
     &           Bvm2R*m3*m4**4*mq1*PsR*tman - Bvm1L*m3**4*PsR*sman*tman + Bvm2L*m3**2*m4**2*PsR*sman*tman + Bvm1L*m4**4*PsR*sman*tman + 
     &           Bvm2L*m4**4*PsR*sman*tman - 2*Bvm1R*m3**3*mq1*PsR*sman*tman + Bvm1R*m3*m4**2*mq1*PsR*sman*tman + Bvm2R*m3*m4**2*mq1*PsR*sman*tman - 
     &           Bsm1L*m3**4*PvL*sman*tman - Bvm1R*m3**5*PvL*sman*tman + Bsm2L*m3**2*m4**2*PvL*sman*tman + Bvm1R*m3**3*m4**2*PvL*sman*tman + 
     &           Bvm2R*m3**3*m4**2*PvL*sman*tman + Bsm1L*m4**4*PvL*sman*tman + Bsm2L*m4**4*PvL*sman*tman + 2*Bvm1R*m3*m4**4*PvL*sman*tman + 
     &           Bvm2R*m3*m4**4*PvL*sman*tman + Bvm1L*m3**4*mq1*PvL*sman*tman - Bvm2L*m3**2*m4**2*mq1*PvL*sman*tman - Bvm1L*m4**4*mq1*PvL*sman*tman - 
     &           Bvm2L*m4**4*mq1*PvL*sman*tman + 2*Bvm1L*m3**2*PsR*sman**2*tman - Bvm2L*m4**2*PsR*sman**2*tman + Bvm1R*m3*mq1*PsR*sman**2*tman + 
     &           2*Bsm1L*m3**2*PvL*sman**2*tman + 2*Bvm1R*m3**3*PvL*sman**2*tman - Bsm2L*m4**2*PvL*sman**2*tman - Bvm1R*m3*m4**2*PvL*sman**2*tman - 
     &           Bvm2R*m3*m4**2*PvL*sman**2*tman - 2*Bvm1L*m3**2*mq1*PvL*sman**2*tman + Bvm2L*m4**2*mq1*PvL*sman**2*tman - Bvm1L*PsR*sman**3*tman - 
     &           Bsm1L*PvL*sman**3*tman - Bvm1R*m3*PvL*sman**3*tman + Bvm1L*mq1*PvL*sman**3*tman + 
     &           m3*(-(Bvm2R*m4**2*(m3**2 + m4**2)*(mq1*PsR - PvL*sman)) + 
     &              m3*sman*(-(Bsm1L*m3**2*PvL) + m4**2*(Bsm2L*PvL + Bvm2L*(PsR - mq1*PvL)) - Bvm1L*(PsR - mq1*PvL)*(m3**2 - sman) + Bsm1L*PvL*sman) + 
     &              Bvm1R*(mq1*PsR - PvL*sman)*(m3**4 - 2*m4**4 + m4**2*sman - m3**2*(m4**2 + sman)))*uman + 
     &           Bsm2R*m3*m4**2*PsR*(-(sman*tman) + m3**2*(tman + uman) + m4**2*(tman + uman)) + 
     &           Bsm1R*m3*PsR*(-(sman**2*tman) - m3**4*(tman + uman) + 2*m4**4*(tman + uman) - m4**2*sman*(tman + uman) + 
     &              m3**2*(m4**2*(tman + uman) + sman*(2*tman + uman))))*xv + 
     &        BvgR*(4*(-m3**4 + m4**4)*mq1*PsR + (-(m4**6*mq1*PsR) + sman**2*(mq1*PsR - PvL*sman)*tman + m4**2*sman*(mq1*PsR - PvL*sman)*uman - 
     &              m4**4*(mq1*PsR - PvL*sman)*(tman + uman) + m3**2*sman*(-(mq1*PsR) + PvL*sman)*(2*tman + uman) + 
     &              m3**4*(m4**2*mq1*PsR + (mq1*PsR - PvL*sman)*(tman + uman)))*xv + 2*m1**2*(m3 - m4)*(m3 + m4)*mq1*PsR*(-4 + m4**2*xv)) + 
     &        BsgR*PsR*(-4*m4**4 + m3**2*sman*(2*tman + uman)*xv + (m4**6 - sman**2*tman - m4**2*sman*uman + m4**4*(tman + uman))*xv - 
     &           2*m1**2*(m3 - m4)*(m3 + m4)*(-4 + m4**2*xv) + m3**4*(4 - (m4**2 + tman + uman)*xv)))) + 
     &  A1Lc*A1sL*(BLc*mq1c*(-(Bsm2L*m1**2*m3*PsL) - 8*BvgR*m1**2*m3*PsL - 4*Bvm1R*m1**2*m3**2*PsL - 4*BvgR*m3**3*PsL - 2*Bvm1R*m3**4*PsL + 
     &        2*Bvm1R*m1**2*m4**2*PsL + 2*Bvm2R*m1**2*m4**2*PsL - Bsm2L*m3*m4**2*PsL - 4*BvgR*m3*m4**2*PsL - Bvm1R*m3**2*m4**2*PsL + Bvm2R*m3**2*m4**2*PsL + 
     &        Bvm1R*m4**4*PsL + Bvm2R*m4**4*PsL - 2*BvgL*m1**2*mq1*PsL + 3*Bvm1L*m1**2*m3*mq1*PsL + Bvm2L*m1**2*m3*mq1*PsL - 2*BvgL*m3**2*mq1*PsL + 
     &        2*Bvm1L*m3**3*mq1*PsL + Bvm1L*m3*m4**2*mq1*PsL + Bvm2L*m3*m4**2*mq1*PsL - 8*BsgR*m1**2*m3*PvR - 4*Bsm1R*m1**2*m3**2*PvR - 4*BsgR*m3**3*PvR - 
     &        2*Bsm1R*m3**4*PvR + 2*Bsm1R*m1**2*m4**2*PvR + 2*Bsm2R*m1**2*m4**2*PvR - 4*BsgR*m3*m4**2*PvR - Bsm1R*m3**2*m4**2*PvR + Bsm2R*m3**2*m4**2*PvR + 
     &        Bsm1R*m4**4*PvR + Bsm2R*m4**4*PvR + 8*BvgR*m1**2*m3*mq1*PvR + 4*Bvm1R*m1**2*m3**2*mq1*PvR + 4*BvgR*m3**3*mq1*PvR + 2*Bvm1R*m3**4*mq1*PvR - 
     &        2*Bvm1R*m1**2*m4**2*mq1*PvR - 2*Bvm2R*m1**2*m4**2*mq1*PvR + 4*BvgR*m3*m4**2*mq1*PvR + Bvm1R*m3**2*m4**2*mq1*PvR - Bvm2R*m3**2*m4**2*mq1*PvR - 
     &        Bvm1R*m4**4*mq1*PvR - Bvm2R*m4**4*mq1*PvR - Bvm1R*m1**2*PsL*sman - Bvm2R*m1**2*PsL*sman - Bvm1R*m4**2*PsL*sman - Bvm2R*m4**2*PsL*sman - 
     &        Bsm1R*m1**2*PvR*sman - Bsm2R*m1**2*PvR*sman + 2*BvgL*m1**2*PvR*sman - 3*Bvm1L*m1**2*m3*PvR*sman - Bvm2L*m1**2*m3*PvR*sman + 
     &        2*BvgL*m3**2*PvR*sman - 2*Bvm1L*m3**3*PvR*sman - Bsm1R*m4**2*PvR*sman - Bsm2R*m4**2*PvR*sman - Bvm1L*m3*m4**2*PvR*sman - 
     &        Bvm2L*m3*m4**2*PvR*sman + Bvm1R*m1**2*mq1*PvR*sman + Bvm2R*m1**2*mq1*PvR*sman + Bvm1R*m4**2*mq1*PvR*sman + Bvm2R*m4**2*mq1*PvR*sman + 
     &        4*BvgR*m3*PsL*tman + 2*Bvm1R*m3**2*PsL*tman - Bvm1R*m4**2*PsL*tman - Bvm2R*m4**2*PsL*tman + 2*BvgL*mq1*PsL*tman - 2*Bvm1L*m3*mq1*PsL*tman + 
     &        4*BsgR*m3*PvR*tman + 2*Bsm1R*m3**2*PvR*tman - Bsm1R*m4**2*PvR*tman - Bsm2R*m4**2*PvR*tman - 4*BvgR*m3*mq1*PvR*tman - 
     &        2*Bvm1R*m3**2*mq1*PvR*tman + Bvm1R*m4**2*mq1*PvR*tman + Bvm2R*m4**2*mq1*PvR*tman - 2*BvgL*PvR*sman*tman + 2*Bvm1L*m3*PvR*sman*tman + 
     &        Bsm2L*m3*PsL*uman + 4*BvgR*m3*PsL*uman + 2*Bvm1R*m3**2*PsL*uman - Bvm1R*m4**2*PsL*uman - Bvm2R*m4**2*PsL*uman - Bvm1L*m3*mq1*PsL*uman - 
     &        Bvm2L*m3*mq1*PsL*uman + 4*BsgR*m3*PvR*uman + 2*Bsm1R*m3**2*PvR*uman - Bsm1R*m4**2*PvR*uman - Bsm2R*m4**2*PvR*uman - 4*BvgR*m3*mq1*PvR*uman - 
     &        2*Bvm1R*m3**2*mq1*PvR*uman + Bvm1R*m4**2*mq1*PvR*uman + Bvm2R*m4**2*mq1*PvR*uman + Bvm1R*PsL*sman*uman + Bvm2R*PsL*sman*uman + 
     &        Bsm1R*PvR*sman*uman + Bsm2R*PvR*sman*uman + Bvm1L*m3*PvR*sman*uman + Bvm2L*m3*PvR*sman*uman - Bvm1R*mq1*PvR*sman*uman - 
     &        Bvm2R*mq1*PvR*sman*uman + Bsm2L*m1**2*m3*m4**2*PsL*xv + 2*BvgR*m1**2*m3*m4**2*PsL*xv + 2*Bvm1R*m1**2*m3**2*m4**2*PsL*xv + 
     &        BvgR*m3**3*m4**2*PsL*xv + Bvm1R*m3**4*m4**2*PsL*xv - 2*Bvm2R*m1**2*m4**4*PsL*xv + Bsm2L*m3*m4**4*PsL*xv + BvgR*m3*m4**4*PsL*xv + 
     &        Bvm1R*m3**2*m4**4*PsL*xv - Bvm2R*m3**2*m4**4*PsL*xv - Bvm2R*m4**6*PsL*xv + BvgL*m1**2*m3**2*mq1*PsL*xv + Bvm1L*m1**2*m3**3*mq1*PsL*xv + 
     &        2*BvgL*m1**2*m4**2*mq1*PsL*xv - Bvm2L*m1**2*m3*m4**2*mq1*PsL*xv + 2*BvgL*m3**2*m4**2*mq1*PsL*xv + Bvm1L*m3**3*m4**2*mq1*PsL*xv + 
     &        BvgL*m4**4*mq1*PsL*xv - Bvm2L*m3*m4**4*mq1*PsL*xv + 2*BsgR*m1**2*m3*m4**2*PvR*xv + 2*Bsm1R*m1**2*m3**2*m4**2*PvR*xv + BsgR*m3**3*m4**2*PvR*xv + 
     &        Bsm1R*m3**4*m4**2*PvR*xv - 2*Bsm2R*m1**2*m4**4*PvR*xv + BsgR*m3*m4**4*PvR*xv + Bsm1R*m3**2*m4**4*PvR*xv - Bsm2R*m3**2*m4**4*PvR*xv - 
     &        Bsm2R*m4**6*PvR*xv - 2*BvgR*m1**2*m3*m4**2*mq1*PvR*xv - 2*Bvm1R*m1**2*m3**2*m4**2*mq1*PvR*xv - BvgR*m3**3*m4**2*mq1*PvR*xv - 
     &        Bvm1R*m3**4*m4**2*mq1*PvR*xv + 2*Bvm2R*m1**2*m4**4*mq1*PvR*xv - BvgR*m3*m4**4*mq1*PvR*xv - Bvm1R*m3**2*m4**4*mq1*PvR*xv + 
     &        Bvm2R*m3**2*m4**4*mq1*PvR*xv + Bvm2R*m4**6*mq1*PvR*xv - Bvm1R*m1**2*m3**2*PsL*sman*xv - 2*Bvm1R*m1**2*m4**2*PsL*sman*xv + 
     &        Bvm2R*m1**2*m4**2*PsL*sman*xv - 2*Bvm1R*m3**2*m4**2*PsL*sman*xv - Bvm1R*m4**4*PsL*sman*xv + Bvm2R*m4**4*PsL*sman*xv - 
     &        BvgL*m1**2*mq1*PsL*sman*xv - Bvm1L*m1**2*m3*mq1*PsL*sman*xv - BvgL*m4**2*mq1*PsL*sman*xv - Bvm1L*m3*m4**2*mq1*PsL*sman*xv - 
     &        Bsm1R*m1**2*m3**2*PvR*sman*xv - BvgL*m1**2*m3**2*PvR*sman*xv - Bvm1L*m1**2*m3**3*PvR*sman*xv - 2*Bsm1R*m1**2*m4**2*PvR*sman*xv + 
     &        Bsm2R*m1**2*m4**2*PvR*sman*xv - 2*BvgL*m1**2*m4**2*PvR*sman*xv + Bvm2L*m1**2*m3*m4**2*PvR*sman*xv - 2*Bsm1R*m3**2*m4**2*PvR*sman*xv - 
     &        2*BvgL*m3**2*m4**2*PvR*sman*xv - Bvm1L*m3**3*m4**2*PvR*sman*xv - Bsm1R*m4**4*PvR*sman*xv + Bsm2R*m4**4*PvR*sman*xv - BvgL*m4**4*PvR*sman*xv + 
     &        Bvm2L*m3*m4**4*PvR*sman*xv + Bvm1R*m1**2*m3**2*mq1*PvR*sman*xv + 2*Bvm1R*m1**2*m4**2*mq1*PvR*sman*xv - Bvm2R*m1**2*m4**2*mq1*PvR*sman*xv + 
     &        2*Bvm1R*m3**2*m4**2*mq1*PvR*sman*xv + Bvm1R*m4**4*mq1*PvR*sman*xv - Bvm2R*m4**4*mq1*PvR*sman*xv + Bvm1R*m1**2*PsL*sman**2*xv + 
     &        Bvm1R*m4**2*PsL*sman**2*xv + Bsm1R*m1**2*PvR*sman**2*xv + BvgL*m1**2*PvR*sman**2*xv + Bvm1L*m1**2*m3*PvR*sman**2*xv + 
     &        Bsm1R*m4**2*PvR*sman**2*xv + BvgL*m4**2*PvR*sman**2*xv + Bvm1L*m3*m4**2*PvR*sman**2*xv - Bvm1R*m1**2*mq1*PvR*sman**2*xv - 
     &        Bvm1R*m4**2*mq1*PvR*sman**2*xv - BvgR*m3*m4**2*PsL*tman*xv - Bvm1R*m3**2*m4**2*PsL*tman*xv + Bvm2R*m4**4*PsL*tman*xv - 
     &        BvgL*m4**2*mq1*PsL*tman*xv - BsgR*m3*m4**2*PvR*tman*xv - Bsm1R*m3**2*m4**2*PvR*tman*xv + Bsm2R*m4**4*PvR*tman*xv + 
     &        BvgR*m3*m4**2*mq1*PvR*tman*xv + Bvm1R*m3**2*m4**2*mq1*PvR*tman*xv - Bvm2R*m4**4*mq1*PvR*tman*xv + Bvm1R*m4**2*PsL*sman*tman*xv + 
     &        Bsm1R*m4**2*PvR*sman*tman*xv + BvgL*m4**2*PvR*sman*tman*xv - Bvm1R*m4**2*mq1*PvR*sman*tman*xv - Bsm2L*m3*m4**2*PsL*uman*xv - 
     &        BvgR*m3*m4**2*PsL*uman*xv - Bvm1R*m3**2*m4**2*PsL*uman*xv + Bvm2R*m4**4*PsL*uman*xv - BvgL*m3**2*mq1*PsL*uman*xv - Bvm1L*m3**3*mq1*PsL*uman*xv - 
     &        BvgL*m4**2*mq1*PsL*uman*xv + Bvm2L*m3*m4**2*mq1*PsL*uman*xv - BsgR*m3*m4**2*PvR*uman*xv - Bsm1R*m3**2*m4**2*PvR*uman*xv + 
     &        Bsm2R*m4**4*PvR*uman*xv + BvgR*m3*m4**2*mq1*PvR*uman*xv + Bvm1R*m3**2*m4**2*mq1*PvR*uman*xv - Bvm2R*m4**4*mq1*PvR*uman*xv + 
     &        Bvm1R*m3**2*PsL*sman*uman*xv + Bvm1R*m4**2*PsL*sman*uman*xv - Bvm2R*m4**2*PsL*sman*uman*xv + BvgL*mq1*PsL*sman*uman*xv + 
     &        Bvm1L*m3*mq1*PsL*sman*uman*xv + Bsm1R*m3**2*PvR*sman*uman*xv + BvgL*m3**2*PvR*sman*uman*xv + Bvm1L*m3**3*PvR*sman*uman*xv + 
     &        Bsm1R*m4**2*PvR*sman*uman*xv - Bsm2R*m4**2*PvR*sman*uman*xv + BvgL*m4**2*PvR*sman*uman*xv - Bvm2L*m3*m4**2*PvR*sman*uman*xv - 
     &        Bvm1R*m3**2*mq1*PvR*sman*uman*xv - Bvm1R*m4**2*mq1*PvR*sman*uman*xv + Bvm2R*m4**2*mq1*PvR*sman*uman*xv - Bvm1R*PsL*sman**2*uman*xv - 
     &        Bsm1R*PvR*sman**2*uman*xv - BvgL*PvR*sman**2*uman*xv - Bvm1L*m3*PvR*sman**2*uman*xv + Bvm1R*mq1*PvR*sman**2*uman*xv - 
     &        Bsm1L*m3*PsL*(m4**2 - 2*tman - uman - m4**2*sman*xv + sman*uman*xv + m1**2*(3 + m3**2*xv - sman*xv) + m3**2*(2 + m4**2*xv - uman*xv)) + 
     &        BsgL*PsL*(-2*tman - m4**4*xv + m4**2*sman*xv + m4**2*tman*xv + m4**2*uman*xv - sman*uman*xv + m1**2*(2 - m3**2*xv - 2*m4**2*xv + sman*xv) + 
     &           m3**2*(2 - 2*m4**2*xv + uman*xv))) + BRc*(2*Bsm2L*m1**2*m3**2*PsL + 4*BvgR*m1**2*m3**2*PsL - 2*Bvm1R*m1**2*m3**3*PsL + 
     &        2*Bvm2R*m1**2*m3**3*PsL + Bsm2L*m3**4*PsL + 2*BvgR*m3**4*PsL - Bvm1R*m3**5*PsL + Bvm2R*m3**5*PsL - 4*BvgR*m1**2*m4**2*PsL + 
     &        2*Bvm1R*m1**2*m3*m4**2*PsL - 2*Bvm2R*m1**2*m3*m4**2*PsL + Bsm2L*m3**2*m4**2*PsL - 2*BvgR*m4**4*PsL + Bvm1R*m3*m4**4*PsL - Bvm2R*m3*m4**4*PsL + 
     &        8*BvgL*m1**2*m3*mq1*PsL + 2*Bvm1L*m1**2*m3**2*mq1*PsL - 2*Bvm2L*m1**2*m3**2*mq1*PsL + 4*BvgL*m3**3*mq1*PsL + Bvm1L*m3**4*mq1*PsL - 
     &        Bvm2L*m3**4*mq1*PsL + 4*BvgL*m3*m4**2*mq1*PsL + Bvm1L*m3**2*m4**2*mq1*PsL - Bvm2L*m3**2*m4**2*mq1*PsL + 4*BsgR*m1**2*m3**2*PvR - 
     &        2*Bsm1R*m1**2*m3**3*PvR + 2*Bsm2R*m1**2*m3**3*PvR + 2*BsgR*m3**4*PvR - Bsm1R*m3**5*PvR + Bsm2R*m3**5*PvR - 4*BsgR*m1**2*m4**2*PvR + 
     &        2*Bsm1R*m1**2*m3*m4**2*PvR - 2*Bsm2R*m1**2*m3*m4**2*PvR - 2*BsgR*m4**4*PvR + Bsm1R*m3*m4**4*PvR - Bsm2R*m3*m4**4*PvR - 
     &        4*BvgR*m1**2*m3**2*mq1*PvR + 2*Bvm1R*m1**2*m3**3*mq1*PvR - 2*Bvm2R*m1**2*m3**3*mq1*PvR - 2*BvgR*m3**4*mq1*PvR + Bvm1R*m3**5*mq1*PvR - 
     &        Bvm2R*m3**5*mq1*PvR + 4*BvgR*m1**2*m4**2*mq1*PvR - 2*Bvm1R*m1**2*m3*m4**2*mq1*PvR + 2*Bvm2R*m1**2*m3*m4**2*mq1*PvR + 2*BvgR*m4**4*mq1*PvR - 
     &        Bvm1R*m3*m4**4*mq1*PvR + Bvm2R*m3*m4**4*mq1*PvR - Bsm2L*m1**2*PsL*sman + 2*BvgR*m1**2*PsL*sman - 3*Bvm1R*m1**2*m3*PsL*sman - 
     &        Bvm2R*m1**2*m3*PsL*sman - Bsm2L*m3**2*PsL*sman - Bvm1R*m3**3*PsL*sman - Bvm2R*m3**3*PsL*sman + 2*BvgR*m4**2*PsL*sman - 
     &        2*Bvm1R*m3*m4**2*PsL*sman + Bvm1L*m1**2*mq1*PsL*sman + Bvm2L*m1**2*mq1*PsL*sman + Bvm1L*m3**2*mq1*PsL*sman + Bvm2L*m3**2*mq1*PsL*sman + 
     &        2*BsgR*m1**2*PvR*sman - 3*Bsm1R*m1**2*m3*PvR*sman - Bsm2R*m1**2*m3*PvR*sman - 8*BvgL*m1**2*m3*PvR*sman - 2*Bvm1L*m1**2*m3**2*PvR*sman + 
     &        2*Bvm2L*m1**2*m3**2*PvR*sman - Bsm1R*m3**3*PvR*sman - Bsm2R*m3**3*PvR*sman - 4*BvgL*m3**3*PvR*sman - Bvm1L*m3**4*PvR*sman + 
     &        Bvm2L*m3**4*PvR*sman + 2*BsgR*m4**2*PvR*sman - 2*Bsm1R*m3*m4**2*PvR*sman - 4*BvgL*m3*m4**2*PvR*sman - Bvm1L*m3**2*m4**2*PvR*sman + 
     &        Bvm2L*m3**2*m4**2*PvR*sman - 2*BvgR*m1**2*mq1*PvR*sman + 3*Bvm1R*m1**2*m3*mq1*PvR*sman + Bvm2R*m1**2*m3*mq1*PvR*sman + 
     &        Bvm1R*m3**3*mq1*PvR*sman + Bvm2R*m3**3*mq1*PvR*sman - 2*BvgR*m4**2*mq1*PvR*sman + 2*Bvm1R*m3*m4**2*mq1*PvR*sman - Bvm1L*m1**2*PvR*sman**2 - 
     &        Bvm2L*m1**2*PvR*sman**2 - Bvm1L*m3**2*PvR*sman**2 - Bvm2L*m3**2*PvR*sman**2 - Bsm2L*m3**2*PsL*tman - 2*BvgR*m3**2*PsL*tman + 
     &        Bvm1R*m3**3*PsL*tman - Bvm2R*m3**3*PsL*tman + 2*BvgR*m4**2*PsL*tman - Bvm1R*m3*m4**2*PsL*tman + Bvm2R*m3*m4**2*PsL*tman - 
     &        4*BvgL*m3*mq1*PsL*tman - Bvm1L*m3**2*mq1*PsL*tman + Bvm2L*m3**2*mq1*PsL*tman - 2*BsgR*m3**2*PvR*tman + Bsm1R*m3**3*PvR*tman - 
     &        Bsm2R*m3**3*PvR*tman + 2*BsgR*m4**2*PvR*tman - Bsm1R*m3*m4**2*PvR*tman + Bsm2R*m3*m4**2*PvR*tman + 2*BvgR*m3**2*mq1*PvR*tman - 
     &        Bvm1R*m3**3*mq1*PvR*tman + Bvm2R*m3**3*mq1*PvR*tman - 2*BvgR*m4**2*mq1*PvR*tman + Bvm1R*m3*m4**2*mq1*PvR*tman - Bvm2R*m3*m4**2*mq1*PvR*tman + 
     &        Bsm2L*PsL*sman*tman + Bvm1R*m3*PsL*sman*tman + Bvm2R*m3*PsL*sman*tman - Bvm1L*mq1*PsL*sman*tman - Bvm2L*mq1*PsL*sman*tman + 
     &        Bsm1R*m3*PvR*sman*tman + Bsm2R*m3*PvR*sman*tman + 4*BvgL*m3*PvR*sman*tman + Bvm1L*m3**2*PvR*sman*tman - Bvm2L*m3**2*PvR*sman*tman - 
     &        Bvm1R*m3*mq1*PvR*sman*tman - Bvm2R*m3*mq1*PvR*sman*tman + Bvm1L*PvR*sman**2*tman + Bvm2L*PvR*sman**2*tman - Bsm2L*m3**2*PsL*uman - 
     &        2*BvgR*m3**2*PsL*uman + Bvm1R*m3**3*PsL*uman - Bvm2R*m3**3*PsL*uman + 2*BvgR*m4**2*PsL*uman - Bvm1R*m3*m4**2*PsL*uman + 
     &        Bvm2R*m3*m4**2*PsL*uman - 4*BvgL*m3*mq1*PsL*uman - Bvm1L*m3**2*mq1*PsL*uman + Bvm2L*m3**2*mq1*PsL*uman - 2*BsgR*m3**2*PvR*uman + 
     &        Bsm1R*m3**3*PvR*uman - Bsm2R*m3**3*PvR*uman + 2*BsgR*m4**2*PvR*uman - Bsm1R*m3*m4**2*PvR*uman + Bsm2R*m3*m4**2*PvR*uman + 
     &        2*BvgR*m3**2*mq1*PvR*uman - Bvm1R*m3**3*mq1*PvR*uman + Bvm2R*m3**3*mq1*PvR*uman - 2*BvgR*m4**2*mq1*PvR*uman + Bvm1R*m3*m4**2*mq1*PvR*uman - 
     &        Bvm2R*m3*m4**2*mq1*PvR*uman - 2*BvgR*PsL*sman*uman + 2*Bvm1R*m3*PsL*sman*uman - 2*BsgR*PvR*sman*uman + 2*Bsm1R*m3*PvR*sman*uman + 
     &        4*BvgL*m3*PvR*sman*uman + Bvm1L*m3**2*PvR*sman*uman - Bvm2L*m3**2*PvR*sman*uman + 2*BvgR*mq1*PvR*sman*uman - 2*Bvm1R*m3*mq1*PvR*sman*uman + 
     &        2*BvgR*m1**2*m3**4*PsL*xv + 2*Bvm1R*m1**2*m3**5*PsL*xv + BvgR*m3**6*PsL*xv + Bvm1R*m3**7*PsL*xv - 2*Bsm2L*m1**2*m3**2*m4**2*PsL*xv - 
     &        2*BvgR*m1**2*m3**2*m4**2*PsL*xv - 2*Bvm1R*m1**2*m3**3*m4**2*PsL*xv - 2*Bvm2R*m1**2*m3**3*m4**2*PsL*xv - Bsm2L*m3**4*m4**2*PsL*xv - 
     &        Bvm2R*m3**5*m4**2*PsL*xv + 2*Bvm2R*m1**2*m3*m4**4*PsL*xv - Bsm2L*m3**2*m4**4*PsL*xv - BvgR*m3**2*m4**4*PsL*xv - Bvm1R*m3**3*m4**4*PsL*xv + 
     &        Bvm2R*m3*m4**6*PsL*xv - 2*Bvm1L*m1**2*m3**4*mq1*PsL*xv - Bvm1L*m3**6*mq1*PsL*xv - 2*BvgL*m1**2*m3*m4**2*mq1*PsL*xv + 
     &        2*Bvm2L*m1**2*m3**2*m4**2*mq1*PsL*xv - BvgL*m3**3*m4**2*mq1*PsL*xv - Bvm1L*m3**4*m4**2*mq1*PsL*xv + Bvm2L*m3**4*m4**2*mq1*PsL*xv - 
     &        BvgL*m3*m4**4*mq1*PsL*xv + Bvm2L*m3**2*m4**4*mq1*PsL*xv + 2*BsgR*m1**2*m3**4*PvR*xv + 2*Bsm1R*m1**2*m3**5*PvR*xv + BsgR*m3**6*PvR*xv + 
     &        Bsm1R*m3**7*PvR*xv - 2*BsgR*m1**2*m3**2*m4**2*PvR*xv - 2*Bsm1R*m1**2*m3**3*m4**2*PvR*xv - 2*Bsm2R*m1**2*m3**3*m4**2*PvR*xv - 
     &        Bsm2R*m3**5*m4**2*PvR*xv + 2*Bsm2R*m1**2*m3*m4**4*PvR*xv - BsgR*m3**2*m4**4*PvR*xv - Bsm1R*m3**3*m4**4*PvR*xv + Bsm2R*m3*m4**6*PvR*xv - 
     &        2*BvgR*m1**2*m3**4*mq1*PvR*xv - 2*Bvm1R*m1**2*m3**5*mq1*PvR*xv - BvgR*m3**6*mq1*PvR*xv - Bvm1R*m3**7*mq1*PvR*xv + 
     &        2*BvgR*m1**2*m3**2*m4**2*mq1*PvR*xv + 2*Bvm1R*m1**2*m3**3*m4**2*mq1*PvR*xv + 2*Bvm2R*m1**2*m3**3*m4**2*mq1*PvR*xv + 
     &        Bvm2R*m3**5*m4**2*mq1*PvR*xv - 2*Bvm2R*m1**2*m3*m4**4*mq1*PvR*xv + BvgR*m3**2*m4**4*mq1*PvR*xv + Bvm1R*m3**3*m4**4*mq1*PvR*xv - 
     &        Bvm2R*m3*m4**6*mq1*PvR*xv - 3*BvgR*m1**2*m3**2*PsL*sman*xv - 3*Bvm1R*m1**2*m3**3*PsL*sman*xv - 2*BvgR*m3**4*PsL*sman*xv - 
     &        2*Bvm1R*m3**5*PsL*sman*xv + Bsm2L*m1**2*m4**2*PsL*sman*xv + 2*Bvm1R*m1**2*m3*m4**2*PsL*sman*xv + Bvm2R*m1**2*m3*m4**2*PsL*sman*xv + 
     &        Bsm2L*m3**2*m4**2*PsL*sman*xv - BvgR*m3**2*m4**2*PsL*sman*xv + Bvm2R*m3**3*m4**2*PsL*sman*xv + Bvm1R*m3*m4**4*PsL*sman*xv + 
     &        3*Bvm1L*m1**2*m3**2*mq1*PsL*sman*xv + 2*Bvm1L*m3**4*mq1*PsL*sman*xv - Bvm2L*m1**2*m4**2*mq1*PsL*sman*xv + Bvm1L*m3**2*m4**2*mq1*PsL*sman*xv - 
     &        Bvm2L*m3**2*m4**2*mq1*PsL*sman*xv - 3*BsgR*m1**2*m3**2*PvR*sman*xv - 3*Bsm1R*m1**2*m3**3*PvR*sman*xv - 2*BsgR*m3**4*PvR*sman*xv + 
     &        2*Bvm1L*m1**2*m3**4*PvR*sman*xv - 2*Bsm1R*m3**5*PvR*sman*xv + Bvm1L*m3**6*PvR*sman*xv + 2*Bsm1R*m1**2*m3*m4**2*PvR*sman*xv + 
     &        Bsm2R*m1**2*m3*m4**2*PvR*sman*xv + 2*BvgL*m1**2*m3*m4**2*PvR*sman*xv - BsgR*m3**2*m4**2*PvR*sman*xv - 2*Bvm2L*m1**2*m3**2*m4**2*PvR*sman*xv + 
     &        Bsm2R*m3**3*m4**2*PvR*sman*xv + BvgL*m3**3*m4**2*PvR*sman*xv + Bvm1L*m3**4*m4**2*PvR*sman*xv - Bvm2L*m3**4*m4**2*PvR*sman*xv + 
     &        Bsm1R*m3*m4**4*PvR*sman*xv + BvgL*m3*m4**4*PvR*sman*xv - Bvm2L*m3**2*m4**4*PvR*sman*xv + 3*BvgR*m1**2*m3**2*mq1*PvR*sman*xv + 
     &        3*Bvm1R*m1**2*m3**3*mq1*PvR*sman*xv + 2*BvgR*m3**4*mq1*PvR*sman*xv + 2*Bvm1R*m3**5*mq1*PvR*sman*xv - 2*Bvm1R*m1**2*m3*m4**2*mq1*PvR*sman*xv - 
     &        Bvm2R*m1**2*m3*m4**2*mq1*PvR*sman*xv + BvgR*m3**2*m4**2*mq1*PvR*sman*xv - Bvm2R*m3**3*m4**2*mq1*PvR*sman*xv - Bvm1R*m3*m4**4*mq1*PvR*sman*xv + 
     &        BvgR*m1**2*PsL*sman**2*xv + Bvm1R*m1**2*m3*PsL*sman**2*xv + BvgR*m3**2*PsL*sman**2*xv + Bvm1R*m3**3*PsL*sman**2*xv - 
     &        Bvm1L*m1**2*mq1*PsL*sman**2*xv - Bvm1L*m3**2*mq1*PsL*sman**2*xv + BsgR*m1**2*PvR*sman**2*xv + Bsm1R*m1**2*m3*PvR*sman**2*xv + 
     &        BsgR*m3**2*PvR*sman**2*xv - 3*Bvm1L*m1**2*m3**2*PvR*sman**2*xv + Bsm1R*m3**3*PvR*sman**2*xv - 2*Bvm1L*m3**4*PvR*sman**2*xv + 
     &        Bvm2L*m1**2*m4**2*PvR*sman**2*xv - Bvm1L*m3**2*m4**2*PvR*sman**2*xv + Bvm2L*m3**2*m4**2*PvR*sman**2*xv - BvgR*m1**2*mq1*PvR*sman**2*xv - 
     &        Bvm1R*m1**2*m3*mq1*PvR*sman**2*xv - BvgR*m3**2*mq1*PvR*sman**2*xv - Bvm1R*m3**3*mq1*PvR*sman**2*xv + Bvm1L*m1**2*PvR*sman**3*xv + 
     &        Bvm1L*m3**2*PvR*sman**3*xv - BvgR*m3**4*PsL*tman*xv - Bvm1R*m3**5*PsL*tman*xv + Bsm2L*m3**2*m4**2*PsL*tman*xv + BvgR*m3**2*m4**2*PsL*tman*xv + 
     &        Bvm1R*m3**3*m4**2*PsL*tman*xv + Bvm2R*m3**3*m4**2*PsL*tman*xv - Bvm2R*m3*m4**4*PsL*tman*xv + Bvm1L*m3**4*mq1*PsL*tman*xv + 
     &        BvgL*m3*m4**2*mq1*PsL*tman*xv - Bvm2L*m3**2*m4**2*mq1*PsL*tman*xv - BsgR*m3**4*PvR*tman*xv - Bsm1R*m3**5*PvR*tman*xv + 
     &        BsgR*m3**2*m4**2*PvR*tman*xv + Bsm1R*m3**3*m4**2*PvR*tman*xv + Bsm2R*m3**3*m4**2*PvR*tman*xv - Bsm2R*m3*m4**4*PvR*tman*xv + 
     &        BvgR*m3**4*mq1*PvR*tman*xv + Bvm1R*m3**5*mq1*PvR*tman*xv - BvgR*m3**2*m4**2*mq1*PvR*tman*xv - Bvm1R*m3**3*m4**2*mq1*PvR*tman*xv - 
     &        Bvm2R*m3**3*m4**2*mq1*PvR*tman*xv + Bvm2R*m3*m4**4*mq1*PvR*tman*xv + 2*BvgR*m3**2*PsL*sman*tman*xv + 2*Bvm1R*m3**3*PsL*sman*tman*xv - 
     &        Bsm2L*m4**2*PsL*sman*tman*xv - Bvm1R*m3*m4**2*PsL*sman*tman*xv - Bvm2R*m3*m4**2*PsL*sman*tman*xv - 2*Bvm1L*m3**2*mq1*PsL*sman*tman*xv + 
     &        Bvm2L*m4**2*mq1*PsL*sman*tman*xv + 2*BsgR*m3**2*PvR*sman*tman*xv + 2*Bsm1R*m3**3*PvR*sman*tman*xv - Bvm1L*m3**4*PvR*sman*tman*xv - 
     &        Bsm1R*m3*m4**2*PvR*sman*tman*xv - Bsm2R*m3*m4**2*PvR*sman*tman*xv - BvgL*m3*m4**2*PvR*sman*tman*xv + Bvm2L*m3**2*m4**2*PvR*sman*tman*xv - 
     &        2*BvgR*m3**2*mq1*PvR*sman*tman*xv - 2*Bvm1R*m3**3*mq1*PvR*sman*tman*xv + Bvm1R*m3*m4**2*mq1*PvR*sman*tman*xv + 
     &        Bvm2R*m3*m4**2*mq1*PvR*sman*tman*xv - BvgR*PsL*sman**2*tman*xv - Bvm1R*m3*PsL*sman**2*tman*xv + Bvm1L*mq1*PsL*sman**2*tman*xv - 
     &        BsgR*PvR*sman**2*tman*xv - Bsm1R*m3*PvR*sman**2*tman*xv + 2*Bvm1L*m3**2*PvR*sman**2*tman*xv - Bvm2L*m4**2*PvR*sman**2*tman*xv + 
     &        BvgR*mq1*PvR*sman**2*tman*xv + Bvm1R*m3*mq1*PvR*sman**2*tman*xv - Bvm1L*PvR*sman**3*tman*xv - BvgR*m3**4*PsL*uman*xv - Bvm1R*m3**5*PsL*uman*xv + 
     &        Bsm2L*m3**2*m4**2*PsL*uman*xv + BvgR*m3**2*m4**2*PsL*uman*xv + Bvm1R*m3**3*m4**2*PsL*uman*xv + Bvm2R*m3**3*m4**2*PsL*uman*xv - 
     &        Bvm2R*m3*m4**4*PsL*uman*xv + Bvm1L*m3**4*mq1*PsL*uman*xv + BvgL*m3*m4**2*mq1*PsL*uman*xv - Bvm2L*m3**2*m4**2*mq1*PsL*uman*xv - 
     &        BsgR*m3**4*PvR*uman*xv - Bsm1R*m3**5*PvR*uman*xv + BsgR*m3**2*m4**2*PvR*uman*xv + Bsm1R*m3**3*m4**2*PvR*uman*xv + 
     &        Bsm2R*m3**3*m4**2*PvR*uman*xv - Bsm2R*m3*m4**4*PvR*uman*xv + BvgR*m3**4*mq1*PvR*uman*xv + Bvm1R*m3**5*mq1*PvR*uman*xv - 
     &        BvgR*m3**2*m4**2*mq1*PvR*uman*xv - Bvm1R*m3**3*m4**2*mq1*PvR*uman*xv - Bvm2R*m3**3*m4**2*mq1*PvR*uman*xv + Bvm2R*m3*m4**4*mq1*PvR*uman*xv + 
     &        BvgR*m3**2*PsL*sman*uman*xv + Bvm1R*m3**3*PsL*sman*uman*xv - Bvm1R*m3*m4**2*PsL*sman*uman*xv - Bvm1L*m3**2*mq1*PsL*sman*uman*xv + 
     &        BsgR*m3**2*PvR*sman*uman*xv + Bsm1R*m3**3*PvR*sman*uman*xv - Bvm1L*m3**4*PvR*sman*uman*xv - Bsm1R*m3*m4**2*PvR*sman*uman*xv - 
     &        BvgL*m3*m4**2*PvR*sman*uman*xv + Bvm2L*m3**2*m4**2*PvR*sman*uman*xv - BvgR*m3**2*mq1*PvR*sman*uman*xv - Bvm1R*m3**3*mq1*PvR*sman*uman*xv + 
     &        Bvm1R*m3*m4**2*mq1*PvR*sman*uman*xv + Bvm1L*m3**2*PvR*sman**2*uman*xv + BsgL*m3*PsL*(2*m1**2 + m3**2 + m4**2 - tman - uman)*(-4 + m4**2*xv) + 
     &        Bsm1L*PsL*(m3**6*xv + sman*tman*(1 - sman*xv) + m3**4*(-1 + m4**2*xv - 2*sman*xv - tman*xv - uman*xv) + 
     &           m1**2*(2*m3**4*xv + sman*(-1 + sman*xv) - m3**2*(2 + 3*sman*xv)) + 
     &           m3**2*(tman + uman + sman**2*xv - m4**2*(1 + sman*xv) + sman*(-1 + 2*tman*xv + uman*xv))))) + 
     &  A1Rc*A1sR*(BRc*mq1c*(-(Bsm2R*m1**2*m3*PsR) - 8*BvgL*m1**2*m3*PsR - 4*Bvm1L*m1**2*m3**2*PsR - 4*BvgL*m3**3*PsR - 2*Bvm1L*m3**4*PsR + 
     &        2*Bvm1L*m1**2*m4**2*PsR + 2*Bvm2L*m1**2*m4**2*PsR - Bsm2R*m3*m4**2*PsR - 4*BvgL*m3*m4**2*PsR - Bvm1L*m3**2*m4**2*PsR + Bvm2L*m3**2*m4**2*PsR + 
     &        Bvm1L*m4**4*PsR + Bvm2L*m4**4*PsR - 2*BvgR*m1**2*mq1*PsR + 3*Bvm1R*m1**2*m3*mq1*PsR + Bvm2R*m1**2*m3*mq1*PsR - 2*BvgR*m3**2*mq1*PsR + 
     &        2*Bvm1R*m3**3*mq1*PsR + Bvm1R*m3*m4**2*mq1*PsR + Bvm2R*m3*m4**2*mq1*PsR - 8*BsgL*m1**2*m3*PvL - 4*Bsm1L*m1**2*m3**2*PvL - 4*BsgL*m3**3*PvL - 
     &        2*Bsm1L*m3**4*PvL + 2*Bsm1L*m1**2*m4**2*PvL + 2*Bsm2L*m1**2*m4**2*PvL - 4*BsgL*m3*m4**2*PvL - Bsm1L*m3**2*m4**2*PvL + Bsm2L*m3**2*m4**2*PvL + 
     &        Bsm1L*m4**4*PvL + Bsm2L*m4**4*PvL + 8*BvgL*m1**2*m3*mq1*PvL + 4*Bvm1L*m1**2*m3**2*mq1*PvL + 4*BvgL*m3**3*mq1*PvL + 2*Bvm1L*m3**4*mq1*PvL - 
     &        2*Bvm1L*m1**2*m4**2*mq1*PvL - 2*Bvm2L*m1**2*m4**2*mq1*PvL + 4*BvgL*m3*m4**2*mq1*PvL + Bvm1L*m3**2*m4**2*mq1*PvL - Bvm2L*m3**2*m4**2*mq1*PvL - 
     &        Bvm1L*m4**4*mq1*PvL - Bvm2L*m4**4*mq1*PvL - Bvm1L*m1**2*PsR*sman - Bvm2L*m1**2*PsR*sman - Bvm1L*m4**2*PsR*sman - Bvm2L*m4**2*PsR*sman - 
     &        Bsm1L*m1**2*PvL*sman - Bsm2L*m1**2*PvL*sman + 2*BvgR*m1**2*PvL*sman - 3*Bvm1R*m1**2*m3*PvL*sman - Bvm2R*m1**2*m3*PvL*sman + 
     &        2*BvgR*m3**2*PvL*sman - 2*Bvm1R*m3**3*PvL*sman - Bsm1L*m4**2*PvL*sman - Bsm2L*m4**2*PvL*sman - Bvm1R*m3*m4**2*PvL*sman - 
     &        Bvm2R*m3*m4**2*PvL*sman + Bvm1L*m1**2*mq1*PvL*sman + Bvm2L*m1**2*mq1*PvL*sman + Bvm1L*m4**2*mq1*PvL*sman + Bvm2L*m4**2*mq1*PvL*sman + 
     &        4*BvgL*m3*PsR*tman + 2*Bvm1L*m3**2*PsR*tman - Bvm1L*m4**2*PsR*tman - Bvm2L*m4**2*PsR*tman + 2*BvgR*mq1*PsR*tman - 2*Bvm1R*m3*mq1*PsR*tman + 
     &        4*BsgL*m3*PvL*tman + 2*Bsm1L*m3**2*PvL*tman - Bsm1L*m4**2*PvL*tman - Bsm2L*m4**2*PvL*tman - 4*BvgL*m3*mq1*PvL*tman - 
     &        2*Bvm1L*m3**2*mq1*PvL*tman + Bvm1L*m4**2*mq1*PvL*tman + Bvm2L*m4**2*mq1*PvL*tman - 2*BvgR*PvL*sman*tman + 2*Bvm1R*m3*PvL*sman*tman + 
     &        Bsm2R*m3*PsR*uman + 4*BvgL*m3*PsR*uman + 2*Bvm1L*m3**2*PsR*uman - Bvm1L*m4**2*PsR*uman - Bvm2L*m4**2*PsR*uman - Bvm1R*m3*mq1*PsR*uman - 
     &        Bvm2R*m3*mq1*PsR*uman + 4*BsgL*m3*PvL*uman + 2*Bsm1L*m3**2*PvL*uman - Bsm1L*m4**2*PvL*uman - Bsm2L*m4**2*PvL*uman - 4*BvgL*m3*mq1*PvL*uman - 
     &        2*Bvm1L*m3**2*mq1*PvL*uman + Bvm1L*m4**2*mq1*PvL*uman + Bvm2L*m4**2*mq1*PvL*uman + Bvm1L*PsR*sman*uman + Bvm2L*PsR*sman*uman + 
     &        Bsm1L*PvL*sman*uman + Bsm2L*PvL*sman*uman + Bvm1R*m3*PvL*sman*uman + Bvm2R*m3*PvL*sman*uman - Bvm1L*mq1*PvL*sman*uman - 
     &        Bvm2L*mq1*PvL*sman*uman + Bsm2R*m1**2*m3*m4**2*PsR*xv + 2*BvgL*m1**2*m3*m4**2*PsR*xv + 2*Bvm1L*m1**2*m3**2*m4**2*PsR*xv + 
     &        BvgL*m3**3*m4**2*PsR*xv + Bvm1L*m3**4*m4**2*PsR*xv - 2*Bvm2L*m1**2*m4**4*PsR*xv + Bsm2R*m3*m4**4*PsR*xv + BvgL*m3*m4**4*PsR*xv + 
     &        Bvm1L*m3**2*m4**4*PsR*xv - Bvm2L*m3**2*m4**4*PsR*xv - Bvm2L*m4**6*PsR*xv + BvgR*m1**2*m3**2*mq1*PsR*xv + Bvm1R*m1**2*m3**3*mq1*PsR*xv + 
     &        2*BvgR*m1**2*m4**2*mq1*PsR*xv - Bvm2R*m1**2*m3*m4**2*mq1*PsR*xv + 2*BvgR*m3**2*m4**2*mq1*PsR*xv + Bvm1R*m3**3*m4**2*mq1*PsR*xv + 
     &        BvgR*m4**4*mq1*PsR*xv - Bvm2R*m3*m4**4*mq1*PsR*xv + 2*BsgL*m1**2*m3*m4**2*PvL*xv + 2*Bsm1L*m1**2*m3**2*m4**2*PvL*xv + BsgL*m3**3*m4**2*PvL*xv + 
     &        Bsm1L*m3**4*m4**2*PvL*xv - 2*Bsm2L*m1**2*m4**4*PvL*xv + BsgL*m3*m4**4*PvL*xv + Bsm1L*m3**2*m4**4*PvL*xv - Bsm2L*m3**2*m4**4*PvL*xv - 
     &        Bsm2L*m4**6*PvL*xv - 2*BvgL*m1**2*m3*m4**2*mq1*PvL*xv - 2*Bvm1L*m1**2*m3**2*m4**2*mq1*PvL*xv - BvgL*m3**3*m4**2*mq1*PvL*xv - 
     &        Bvm1L*m3**4*m4**2*mq1*PvL*xv + 2*Bvm2L*m1**2*m4**4*mq1*PvL*xv - BvgL*m3*m4**4*mq1*PvL*xv - Bvm1L*m3**2*m4**4*mq1*PvL*xv + 
     &        Bvm2L*m3**2*m4**4*mq1*PvL*xv + Bvm2L*m4**6*mq1*PvL*xv - Bvm1L*m1**2*m3**2*PsR*sman*xv - 2*Bvm1L*m1**2*m4**2*PsR*sman*xv + 
     &        Bvm2L*m1**2*m4**2*PsR*sman*xv - 2*Bvm1L*m3**2*m4**2*PsR*sman*xv - Bvm1L*m4**4*PsR*sman*xv + Bvm2L*m4**4*PsR*sman*xv - 
     &        BvgR*m1**2*mq1*PsR*sman*xv - Bvm1R*m1**2*m3*mq1*PsR*sman*xv - BvgR*m4**2*mq1*PsR*sman*xv - Bvm1R*m3*m4**2*mq1*PsR*sman*xv - 
     &        Bsm1L*m1**2*m3**2*PvL*sman*xv - BvgR*m1**2*m3**2*PvL*sman*xv - Bvm1R*m1**2*m3**3*PvL*sman*xv - 2*Bsm1L*m1**2*m4**2*PvL*sman*xv + 
     &        Bsm2L*m1**2*m4**2*PvL*sman*xv - 2*BvgR*m1**2*m4**2*PvL*sman*xv + Bvm2R*m1**2*m3*m4**2*PvL*sman*xv - 2*Bsm1L*m3**2*m4**2*PvL*sman*xv - 
     &        2*BvgR*m3**2*m4**2*PvL*sman*xv - Bvm1R*m3**3*m4**2*PvL*sman*xv - Bsm1L*m4**4*PvL*sman*xv + Bsm2L*m4**4*PvL*sman*xv - BvgR*m4**4*PvL*sman*xv + 
     &        Bvm2R*m3*m4**4*PvL*sman*xv + Bvm1L*m1**2*m3**2*mq1*PvL*sman*xv + 2*Bvm1L*m1**2*m4**2*mq1*PvL*sman*xv - Bvm2L*m1**2*m4**2*mq1*PvL*sman*xv + 
     &        2*Bvm1L*m3**2*m4**2*mq1*PvL*sman*xv + Bvm1L*m4**4*mq1*PvL*sman*xv - Bvm2L*m4**4*mq1*PvL*sman*xv + Bvm1L*m1**2*PsR*sman**2*xv + 
     &        Bvm1L*m4**2*PsR*sman**2*xv + Bsm1L*m1**2*PvL*sman**2*xv + BvgR*m1**2*PvL*sman**2*xv + Bvm1R*m1**2*m3*PvL*sman**2*xv + 
     &        Bsm1L*m4**2*PvL*sman**2*xv + BvgR*m4**2*PvL*sman**2*xv + Bvm1R*m3*m4**2*PvL*sman**2*xv - Bvm1L*m1**2*mq1*PvL*sman**2*xv - 
     &        Bvm1L*m4**2*mq1*PvL*sman**2*xv - BvgL*m3*m4**2*PsR*tman*xv - Bvm1L*m3**2*m4**2*PsR*tman*xv + Bvm2L*m4**4*PsR*tman*xv - 
     &        BvgR*m4**2*mq1*PsR*tman*xv - BsgL*m3*m4**2*PvL*tman*xv - Bsm1L*m3**2*m4**2*PvL*tman*xv + Bsm2L*m4**4*PvL*tman*xv + 
     &        BvgL*m3*m4**2*mq1*PvL*tman*xv + Bvm1L*m3**2*m4**2*mq1*PvL*tman*xv - Bvm2L*m4**4*mq1*PvL*tman*xv + Bvm1L*m4**2*PsR*sman*tman*xv + 
     &        Bsm1L*m4**2*PvL*sman*tman*xv + BvgR*m4**2*PvL*sman*tman*xv - Bvm1L*m4**2*mq1*PvL*sman*tman*xv - Bsm2R*m3*m4**2*PsR*uman*xv - 
     &        BvgL*m3*m4**2*PsR*uman*xv - Bvm1L*m3**2*m4**2*PsR*uman*xv + Bvm2L*m4**4*PsR*uman*xv - BvgR*m3**2*mq1*PsR*uman*xv - Bvm1R*m3**3*mq1*PsR*uman*xv - 
     &        BvgR*m4**2*mq1*PsR*uman*xv + Bvm2R*m3*m4**2*mq1*PsR*uman*xv - BsgL*m3*m4**2*PvL*uman*xv - Bsm1L*m3**2*m4**2*PvL*uman*xv + 
     &        Bsm2L*m4**4*PvL*uman*xv + BvgL*m3*m4**2*mq1*PvL*uman*xv + Bvm1L*m3**2*m4**2*mq1*PvL*uman*xv - Bvm2L*m4**4*mq1*PvL*uman*xv + 
     &        Bvm1L*m3**2*PsR*sman*uman*xv + Bvm1L*m4**2*PsR*sman*uman*xv - Bvm2L*m4**2*PsR*sman*uman*xv + BvgR*mq1*PsR*sman*uman*xv + 
     &        Bvm1R*m3*mq1*PsR*sman*uman*xv + Bsm1L*m3**2*PvL*sman*uman*xv + BvgR*m3**2*PvL*sman*uman*xv + Bvm1R*m3**3*PvL*sman*uman*xv + 
     &        Bsm1L*m4**2*PvL*sman*uman*xv - Bsm2L*m4**2*PvL*sman*uman*xv + BvgR*m4**2*PvL*sman*uman*xv - Bvm2R*m3*m4**2*PvL*sman*uman*xv - 
     &        Bvm1L*m3**2*mq1*PvL*sman*uman*xv - Bvm1L*m4**2*mq1*PvL*sman*uman*xv + Bvm2L*m4**2*mq1*PvL*sman*uman*xv - Bvm1L*PsR*sman**2*uman*xv - 
     &        Bsm1L*PvL*sman**2*uman*xv - BvgR*PvL*sman**2*uman*xv - Bvm1R*m3*PvL*sman**2*uman*xv + Bvm1L*mq1*PvL*sman**2*uman*xv - 
     &        Bsm1R*m3*PsR*(m4**2 - 2*tman - uman - m4**2*sman*xv + sman*uman*xv + m1**2*(3 + m3**2*xv - sman*xv) + m3**2*(2 + m4**2*xv - uman*xv)) + 
     &        BsgR*PsR*(-2*tman - m4**4*xv + m4**2*sman*xv + m4**2*tman*xv + m4**2*uman*xv - sman*uman*xv + m1**2*(2 - m3**2*xv - 2*m4**2*xv + sman*xv) + 
     &           m3**2*(2 - 2*m4**2*xv + uman*xv))) + BLc*(2*Bsm2R*m1**2*m3**2*PsR + 4*BvgL*m1**2*m3**2*PsR - 2*Bvm1L*m1**2*m3**3*PsR + 
     &        2*Bvm2L*m1**2*m3**3*PsR + Bsm2R*m3**4*PsR + 2*BvgL*m3**4*PsR - Bvm1L*m3**5*PsR + Bvm2L*m3**5*PsR - 4*BvgL*m1**2*m4**2*PsR + 
     &        2*Bvm1L*m1**2*m3*m4**2*PsR - 2*Bvm2L*m1**2*m3*m4**2*PsR + Bsm2R*m3**2*m4**2*PsR - 2*BvgL*m4**4*PsR + Bvm1L*m3*m4**4*PsR - Bvm2L*m3*m4**4*PsR + 
     &        8*BvgR*m1**2*m3*mq1*PsR + 2*Bvm1R*m1**2*m3**2*mq1*PsR - 2*Bvm2R*m1**2*m3**2*mq1*PsR + 4*BvgR*m3**3*mq1*PsR + Bvm1R*m3**4*mq1*PsR - 
     &        Bvm2R*m3**4*mq1*PsR + 4*BvgR*m3*m4**2*mq1*PsR + Bvm1R*m3**2*m4**2*mq1*PsR - Bvm2R*m3**2*m4**2*mq1*PsR + 4*BsgL*m1**2*m3**2*PvL - 
     &        2*Bsm1L*m1**2*m3**3*PvL + 2*Bsm2L*m1**2*m3**3*PvL + 2*BsgL*m3**4*PvL - Bsm1L*m3**5*PvL + Bsm2L*m3**5*PvL - 4*BsgL*m1**2*m4**2*PvL + 
     &        2*Bsm1L*m1**2*m3*m4**2*PvL - 2*Bsm2L*m1**2*m3*m4**2*PvL - 2*BsgL*m4**4*PvL + Bsm1L*m3*m4**4*PvL - Bsm2L*m3*m4**4*PvL - 
     &        4*BvgL*m1**2*m3**2*mq1*PvL + 2*Bvm1L*m1**2*m3**3*mq1*PvL - 2*Bvm2L*m1**2*m3**3*mq1*PvL - 2*BvgL*m3**4*mq1*PvL + Bvm1L*m3**5*mq1*PvL - 
     &        Bvm2L*m3**5*mq1*PvL + 4*BvgL*m1**2*m4**2*mq1*PvL - 2*Bvm1L*m1**2*m3*m4**2*mq1*PvL + 2*Bvm2L*m1**2*m3*m4**2*mq1*PvL + 2*BvgL*m4**4*mq1*PvL - 
     &        Bvm1L*m3*m4**4*mq1*PvL + Bvm2L*m3*m4**4*mq1*PvL - Bsm2R*m1**2*PsR*sman + 2*BvgL*m1**2*PsR*sman - 3*Bvm1L*m1**2*m3*PsR*sman - 
     &        Bvm2L*m1**2*m3*PsR*sman - Bsm2R*m3**2*PsR*sman - Bvm1L*m3**3*PsR*sman - Bvm2L*m3**3*PsR*sman + 2*BvgL*m4**2*PsR*sman - 
     &        2*Bvm1L*m3*m4**2*PsR*sman + Bvm1R*m1**2*mq1*PsR*sman + Bvm2R*m1**2*mq1*PsR*sman + Bvm1R*m3**2*mq1*PsR*sman + Bvm2R*m3**2*mq1*PsR*sman + 
     &        2*BsgL*m1**2*PvL*sman - 3*Bsm1L*m1**2*m3*PvL*sman - Bsm2L*m1**2*m3*PvL*sman - 8*BvgR*m1**2*m3*PvL*sman - 2*Bvm1R*m1**2*m3**2*PvL*sman + 
     &        2*Bvm2R*m1**2*m3**2*PvL*sman - Bsm1L*m3**3*PvL*sman - Bsm2L*m3**3*PvL*sman - 4*BvgR*m3**3*PvL*sman - Bvm1R*m3**4*PvL*sman + 
     &        Bvm2R*m3**4*PvL*sman + 2*BsgL*m4**2*PvL*sman - 2*Bsm1L*m3*m4**2*PvL*sman - 4*BvgR*m3*m4**2*PvL*sman - Bvm1R*m3**2*m4**2*PvL*sman + 
     &        Bvm2R*m3**2*m4**2*PvL*sman - 2*BvgL*m1**2*mq1*PvL*sman + 3*Bvm1L*m1**2*m3*mq1*PvL*sman + Bvm2L*m1**2*m3*mq1*PvL*sman + 
     &        Bvm1L*m3**3*mq1*PvL*sman + Bvm2L*m3**3*mq1*PvL*sman - 2*BvgL*m4**2*mq1*PvL*sman + 2*Bvm1L*m3*m4**2*mq1*PvL*sman - Bvm1R*m1**2*PvL*sman**2 - 
     &        Bvm2R*m1**2*PvL*sman**2 - Bvm1R*m3**2*PvL*sman**2 - Bvm2R*m3**2*PvL*sman**2 - Bsm2R*m3**2*PsR*tman - 2*BvgL*m3**2*PsR*tman + 
     &        Bvm1L*m3**3*PsR*tman - Bvm2L*m3**3*PsR*tman + 2*BvgL*m4**2*PsR*tman - Bvm1L*m3*m4**2*PsR*tman + Bvm2L*m3*m4**2*PsR*tman - 
     &        4*BvgR*m3*mq1*PsR*tman - Bvm1R*m3**2*mq1*PsR*tman + Bvm2R*m3**2*mq1*PsR*tman - 2*BsgL*m3**2*PvL*tman + Bsm1L*m3**3*PvL*tman - 
     &        Bsm2L*m3**3*PvL*tman + 2*BsgL*m4**2*PvL*tman - Bsm1L*m3*m4**2*PvL*tman + Bsm2L*m3*m4**2*PvL*tman + 2*BvgL*m3**2*mq1*PvL*tman - 
     &        Bvm1L*m3**3*mq1*PvL*tman + Bvm2L*m3**3*mq1*PvL*tman - 2*BvgL*m4**2*mq1*PvL*tman + Bvm1L*m3*m4**2*mq1*PvL*tman - Bvm2L*m3*m4**2*mq1*PvL*tman + 
     &        Bsm2R*PsR*sman*tman + Bvm1L*m3*PsR*sman*tman + Bvm2L*m3*PsR*sman*tman - Bvm1R*mq1*PsR*sman*tman - Bvm2R*mq1*PsR*sman*tman + 
     &        Bsm1L*m3*PvL*sman*tman + Bsm2L*m3*PvL*sman*tman + 4*BvgR*m3*PvL*sman*tman + Bvm1R*m3**2*PvL*sman*tman - Bvm2R*m3**2*PvL*sman*tman - 
     &        Bvm1L*m3*mq1*PvL*sman*tman - Bvm2L*m3*mq1*PvL*sman*tman + Bvm1R*PvL*sman**2*tman + Bvm2R*PvL*sman**2*tman - Bsm2R*m3**2*PsR*uman - 
     &        2*BvgL*m3**2*PsR*uman + Bvm1L*m3**3*PsR*uman - Bvm2L*m3**3*PsR*uman + 2*BvgL*m4**2*PsR*uman - Bvm1L*m3*m4**2*PsR*uman + 
     &        Bvm2L*m3*m4**2*PsR*uman - 4*BvgR*m3*mq1*PsR*uman - Bvm1R*m3**2*mq1*PsR*uman + Bvm2R*m3**2*mq1*PsR*uman - 2*BsgL*m3**2*PvL*uman + 
     &        Bsm1L*m3**3*PvL*uman - Bsm2L*m3**3*PvL*uman + 2*BsgL*m4**2*PvL*uman - Bsm1L*m3*m4**2*PvL*uman + Bsm2L*m3*m4**2*PvL*uman + 
     &        2*BvgL*m3**2*mq1*PvL*uman - Bvm1L*m3**3*mq1*PvL*uman + Bvm2L*m3**3*mq1*PvL*uman - 2*BvgL*m4**2*mq1*PvL*uman + Bvm1L*m3*m4**2*mq1*PvL*uman - 
     &        Bvm2L*m3*m4**2*mq1*PvL*uman - 2*BvgL*PsR*sman*uman + 2*Bvm1L*m3*PsR*sman*uman - 2*BsgL*PvL*sman*uman + 2*Bsm1L*m3*PvL*sman*uman + 
     &        4*BvgR*m3*PvL*sman*uman + Bvm1R*m3**2*PvL*sman*uman - Bvm2R*m3**2*PvL*sman*uman + 2*BvgL*mq1*PvL*sman*uman - 2*Bvm1L*m3*mq1*PvL*sman*uman + 
     &        2*BvgL*m1**2*m3**4*PsR*xv + 2*Bvm1L*m1**2*m3**5*PsR*xv + BvgL*m3**6*PsR*xv + Bvm1L*m3**7*PsR*xv - 2*Bsm2R*m1**2*m3**2*m4**2*PsR*xv - 
     &        2*BvgL*m1**2*m3**2*m4**2*PsR*xv - 2*Bvm1L*m1**2*m3**3*m4**2*PsR*xv - 2*Bvm2L*m1**2*m3**3*m4**2*PsR*xv - Bsm2R*m3**4*m4**2*PsR*xv - 
     &        Bvm2L*m3**5*m4**2*PsR*xv + 2*Bvm2L*m1**2*m3*m4**4*PsR*xv - Bsm2R*m3**2*m4**4*PsR*xv - BvgL*m3**2*m4**4*PsR*xv - Bvm1L*m3**3*m4**4*PsR*xv + 
     &        Bvm2L*m3*m4**6*PsR*xv - 2*Bvm1R*m1**2*m3**4*mq1*PsR*xv - Bvm1R*m3**6*mq1*PsR*xv - 2*BvgR*m1**2*m3*m4**2*mq1*PsR*xv + 
     &        2*Bvm2R*m1**2*m3**2*m4**2*mq1*PsR*xv - BvgR*m3**3*m4**2*mq1*PsR*xv - Bvm1R*m3**4*m4**2*mq1*PsR*xv + Bvm2R*m3**4*m4**2*mq1*PsR*xv - 
     &        BvgR*m3*m4**4*mq1*PsR*xv + Bvm2R*m3**2*m4**4*mq1*PsR*xv + 2*BsgL*m1**2*m3**4*PvL*xv + 2*Bsm1L*m1**2*m3**5*PvL*xv + BsgL*m3**6*PvL*xv + 
     &        Bsm1L*m3**7*PvL*xv - 2*BsgL*m1**2*m3**2*m4**2*PvL*xv - 2*Bsm1L*m1**2*m3**3*m4**2*PvL*xv - 2*Bsm2L*m1**2*m3**3*m4**2*PvL*xv - 
     &        Bsm2L*m3**5*m4**2*PvL*xv + 2*Bsm2L*m1**2*m3*m4**4*PvL*xv - BsgL*m3**2*m4**4*PvL*xv - Bsm1L*m3**3*m4**4*PvL*xv + Bsm2L*m3*m4**6*PvL*xv - 
     &        2*BvgL*m1**2*m3**4*mq1*PvL*xv - 2*Bvm1L*m1**2*m3**5*mq1*PvL*xv - BvgL*m3**6*mq1*PvL*xv - Bvm1L*m3**7*mq1*PvL*xv + 
     &        2*BvgL*m1**2*m3**2*m4**2*mq1*PvL*xv + 2*Bvm1L*m1**2*m3**3*m4**2*mq1*PvL*xv + 2*Bvm2L*m1**2*m3**3*m4**2*mq1*PvL*xv + 
     &        Bvm2L*m3**5*m4**2*mq1*PvL*xv - 2*Bvm2L*m1**2*m3*m4**4*mq1*PvL*xv + BvgL*m3**2*m4**4*mq1*PvL*xv + Bvm1L*m3**3*m4**4*mq1*PvL*xv - 
     &        Bvm2L*m3*m4**6*mq1*PvL*xv - 3*BvgL*m1**2*m3**2*PsR*sman*xv - 3*Bvm1L*m1**2*m3**3*PsR*sman*xv - 2*BvgL*m3**4*PsR*sman*xv - 
     &        2*Bvm1L*m3**5*PsR*sman*xv + Bsm2R*m1**2*m4**2*PsR*sman*xv + 2*Bvm1L*m1**2*m3*m4**2*PsR*sman*xv + Bvm2L*m1**2*m3*m4**2*PsR*sman*xv + 
     &        Bsm2R*m3**2*m4**2*PsR*sman*xv - BvgL*m3**2*m4**2*PsR*sman*xv + Bvm2L*m3**3*m4**2*PsR*sman*xv + Bvm1L*m3*m4**4*PsR*sman*xv + 
     &        3*Bvm1R*m1**2*m3**2*mq1*PsR*sman*xv + 2*Bvm1R*m3**4*mq1*PsR*sman*xv - Bvm2R*m1**2*m4**2*mq1*PsR*sman*xv + Bvm1R*m3**2*m4**2*mq1*PsR*sman*xv - 
     &        Bvm2R*m3**2*m4**2*mq1*PsR*sman*xv - 3*BsgL*m1**2*m3**2*PvL*sman*xv - 3*Bsm1L*m1**2*m3**3*PvL*sman*xv - 2*BsgL*m3**4*PvL*sman*xv + 
     &        2*Bvm1R*m1**2*m3**4*PvL*sman*xv - 2*Bsm1L*m3**5*PvL*sman*xv + Bvm1R*m3**6*PvL*sman*xv + 2*Bsm1L*m1**2*m3*m4**2*PvL*sman*xv + 
     &        Bsm2L*m1**2*m3*m4**2*PvL*sman*xv + 2*BvgR*m1**2*m3*m4**2*PvL*sman*xv - BsgL*m3**2*m4**2*PvL*sman*xv - 2*Bvm2R*m1**2*m3**2*m4**2*PvL*sman*xv + 
     &        Bsm2L*m3**3*m4**2*PvL*sman*xv + BvgR*m3**3*m4**2*PvL*sman*xv + Bvm1R*m3**4*m4**2*PvL*sman*xv - Bvm2R*m3**4*m4**2*PvL*sman*xv + 
     &        Bsm1L*m3*m4**4*PvL*sman*xv + BvgR*m3*m4**4*PvL*sman*xv - Bvm2R*m3**2*m4**4*PvL*sman*xv + 3*BvgL*m1**2*m3**2*mq1*PvL*sman*xv + 
     &        3*Bvm1L*m1**2*m3**3*mq1*PvL*sman*xv + 2*BvgL*m3**4*mq1*PvL*sman*xv + 2*Bvm1L*m3**5*mq1*PvL*sman*xv - 2*Bvm1L*m1**2*m3*m4**2*mq1*PvL*sman*xv - 
     &        Bvm2L*m1**2*m3*m4**2*mq1*PvL*sman*xv + BvgL*m3**2*m4**2*mq1*PvL*sman*xv - Bvm2L*m3**3*m4**2*mq1*PvL*sman*xv - Bvm1L*m3*m4**4*mq1*PvL*sman*xv + 
     &        BvgL*m1**2*PsR*sman**2*xv + Bvm1L*m1**2*m3*PsR*sman**2*xv + BvgL*m3**2*PsR*sman**2*xv + Bvm1L*m3**3*PsR*sman**2*xv - 
     &        Bvm1R*m1**2*mq1*PsR*sman**2*xv - Bvm1R*m3**2*mq1*PsR*sman**2*xv + BsgL*m1**2*PvL*sman**2*xv + Bsm1L*m1**2*m3*PvL*sman**2*xv + 
     &        BsgL*m3**2*PvL*sman**2*xv - 3*Bvm1R*m1**2*m3**2*PvL*sman**2*xv + Bsm1L*m3**3*PvL*sman**2*xv - 2*Bvm1R*m3**4*PvL*sman**2*xv + 
     &        Bvm2R*m1**2*m4**2*PvL*sman**2*xv - Bvm1R*m3**2*m4**2*PvL*sman**2*xv + Bvm2R*m3**2*m4**2*PvL*sman**2*xv - BvgL*m1**2*mq1*PvL*sman**2*xv - 
     &        Bvm1L*m1**2*m3*mq1*PvL*sman**2*xv - BvgL*m3**2*mq1*PvL*sman**2*xv - Bvm1L*m3**3*mq1*PvL*sman**2*xv + Bvm1R*m1**2*PvL*sman**3*xv + 
     &        Bvm1R*m3**2*PvL*sman**3*xv - BvgL*m3**4*PsR*tman*xv - Bvm1L*m3**5*PsR*tman*xv + Bsm2R*m3**2*m4**2*PsR*tman*xv + BvgL*m3**2*m4**2*PsR*tman*xv + 
     &        Bvm1L*m3**3*m4**2*PsR*tman*xv + Bvm2L*m3**3*m4**2*PsR*tman*xv - Bvm2L*m3*m4**4*PsR*tman*xv + Bvm1R*m3**4*mq1*PsR*tman*xv + 
     &        BvgR*m3*m4**2*mq1*PsR*tman*xv - Bvm2R*m3**2*m4**2*mq1*PsR*tman*xv - BsgL*m3**4*PvL*tman*xv - Bsm1L*m3**5*PvL*tman*xv + 
     &        BsgL*m3**2*m4**2*PvL*tman*xv + Bsm1L*m3**3*m4**2*PvL*tman*xv + Bsm2L*m3**3*m4**2*PvL*tman*xv - Bsm2L*m3*m4**4*PvL*tman*xv + 
     &        BvgL*m3**4*mq1*PvL*tman*xv + Bvm1L*m3**5*mq1*PvL*tman*xv - BvgL*m3**2*m4**2*mq1*PvL*tman*xv - Bvm1L*m3**3*m4**2*mq1*PvL*tman*xv - 
     &        Bvm2L*m3**3*m4**2*mq1*PvL*tman*xv + Bvm2L*m3*m4**4*mq1*PvL*tman*xv + 2*BvgL*m3**2*PsR*sman*tman*xv + 2*Bvm1L*m3**3*PsR*sman*tman*xv - 
     &        Bsm2R*m4**2*PsR*sman*tman*xv - Bvm1L*m3*m4**2*PsR*sman*tman*xv - Bvm2L*m3*m4**2*PsR*sman*tman*xv - 2*Bvm1R*m3**2*mq1*PsR*sman*tman*xv + 
     &        Bvm2R*m4**2*mq1*PsR*sman*tman*xv + 2*BsgL*m3**2*PvL*sman*tman*xv + 2*Bsm1L*m3**3*PvL*sman*tman*xv - Bvm1R*m3**4*PvL*sman*tman*xv - 
     &        Bsm1L*m3*m4**2*PvL*sman*tman*xv - Bsm2L*m3*m4**2*PvL*sman*tman*xv - BvgR*m3*m4**2*PvL*sman*tman*xv + Bvm2R*m3**2*m4**2*PvL*sman*tman*xv - 
     &        2*BvgL*m3**2*mq1*PvL*sman*tman*xv - 2*Bvm1L*m3**3*mq1*PvL*sman*tman*xv + Bvm1L*m3*m4**2*mq1*PvL*sman*tman*xv + 
     &        Bvm2L*m3*m4**2*mq1*PvL*sman*tman*xv - BvgL*PsR*sman**2*tman*xv - Bvm1L*m3*PsR*sman**2*tman*xv + Bvm1R*mq1*PsR*sman**2*tman*xv - 
     &        BsgL*PvL*sman**2*tman*xv - Bsm1L*m3*PvL*sman**2*tman*xv + 2*Bvm1R*m3**2*PvL*sman**2*tman*xv - Bvm2R*m4**2*PvL*sman**2*tman*xv + 
     &        BvgL*mq1*PvL*sman**2*tman*xv + Bvm1L*m3*mq1*PvL*sman**2*tman*xv - Bvm1R*PvL*sman**3*tman*xv - BvgL*m3**4*PsR*uman*xv - Bvm1L*m3**5*PsR*uman*xv + 
     &        Bsm2R*m3**2*m4**2*PsR*uman*xv + BvgL*m3**2*m4**2*PsR*uman*xv + Bvm1L*m3**3*m4**2*PsR*uman*xv + Bvm2L*m3**3*m4**2*PsR*uman*xv - 
     &        Bvm2L*m3*m4**4*PsR*uman*xv + Bvm1R*m3**4*mq1*PsR*uman*xv + BvgR*m3*m4**2*mq1*PsR*uman*xv - Bvm2R*m3**2*m4**2*mq1*PsR*uman*xv - 
     &        BsgL*m3**4*PvL*uman*xv - Bsm1L*m3**5*PvL*uman*xv + BsgL*m3**2*m4**2*PvL*uman*xv + Bsm1L*m3**3*m4**2*PvL*uman*xv + 
     &        Bsm2L*m3**3*m4**2*PvL*uman*xv - Bsm2L*m3*m4**4*PvL*uman*xv + BvgL*m3**4*mq1*PvL*uman*xv + Bvm1L*m3**5*mq1*PvL*uman*xv - 
     &        BvgL*m3**2*m4**2*mq1*PvL*uman*xv - Bvm1L*m3**3*m4**2*mq1*PvL*uman*xv - Bvm2L*m3**3*m4**2*mq1*PvL*uman*xv + Bvm2L*m3*m4**4*mq1*PvL*uman*xv + 
     &        BvgL*m3**2*PsR*sman*uman*xv + Bvm1L*m3**3*PsR*sman*uman*xv - Bvm1L*m3*m4**2*PsR*sman*uman*xv - Bvm1R*m3**2*mq1*PsR*sman*uman*xv + 
     &        BsgL*m3**2*PvL*sman*uman*xv + Bsm1L*m3**3*PvL*sman*uman*xv - Bvm1R*m3**4*PvL*sman*uman*xv - Bsm1L*m3*m4**2*PvL*sman*uman*xv - 
     &        BvgR*m3*m4**2*PvL*sman*uman*xv + Bvm2R*m3**2*m4**2*PvL*sman*uman*xv - BvgL*m3**2*mq1*PvL*sman*uman*xv - Bvm1L*m3**3*mq1*PvL*sman*uman*xv + 
     &        Bvm1L*m3*m4**2*mq1*PvL*sman*uman*xv + Bvm1R*m3**2*PvL*sman**2*uman*xv + BsgR*m3*PsR*(2*m1**2 + m3**2 + m4**2 - tman - uman)*(-4 + m4**2*xv) + 
     &        Bsm1R*PsR*(m3**6*xv + sman*tman*(1 - sman*xv) + m3**4*(-1 + m4**2*xv - 2*sman*xv - tman*xv - uman*xv) + 
     &           m1**2*(2*m3**4*xv + sman*(-1 + sman*xv) - m3**2*(2 + 3*sman*xv)) + 
     &           m3**2*(tman + uman + sman**2*xv - m4**2*(1 + sman*xv) + sman*(-1 + 2*tman*xv + uman*xv))))) )/
     &     (( sman - mq1**2 + cI*mq1*wquark)*(sman - mq1c**2 - cI*mq1c*wquark)))
     
           end
          
          ! --- Mst  ---

        double precision function MstV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A2Lc,A2Rc,CLc,mq1,mq2c,xv,wquark)

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

        double complex A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A2Lc,A2Rc,CLc
        double precision mq1,mq2c,xv,wquark

        MstV = 0d0
        MstV = dreal((0.5*(2*A1vR*A2Lc*CLc*m1*(-(Bsm2L*m2**2*m3**2*PsL) + Bsm2L*m2**2*m4**2*PsL + 4*BvgL*m2**2*m3*mq1*PsL + 3*Bvm1L*m2**2*m3**2*mq1*PsL + 
     &     Bvm2L*m2**2*m3**2*mq1*PsL + 3*BvgL*m3**3*mq1*PsL + 3*Bvm1L*m3**4*mq1*PsL - 3*Bvm1L*m2**2*m4**2*mq1*PsL - Bvm2L*m2**2*m4**2*mq1*PsL + 
     &     BvgL*m3*m4**2*mq1*PsL - 2*Bvm1L*m3**2*m4**2*mq1*PsL - Bvm1L*m4**4*mq1*PsL - Bsm2L*m2**2*PsL*sman - 2*BvgR*m2**2*PsL*sman - 
     &     6*Bvm1R*m2**2*m3*PsL*sman - 2*Bvm2R*m2**2*m3*PsL*sman - 3*BvgR*m3**2*PsL*sman - 6*Bvm1R*m3**3*PsL*sman - BvgR*m4**2*PsL*sman - 
     &     2*Bvm1R*m3*m4**2*PsL*sman + 3*Bvm1L*m2**2*mq1*PsL*sman + Bvm2L*m2**2*mq1*PsL*sman - BvgL*m3*mq1*PsL*sman + 2*Bvm1L*m3**2*mq1*PsL*sman + 
     &     2*Bvm1L*m4**2*mq1*PsL*sman - 2*BsgR*m2**2*PvR*sman - 6*Bsm1R*m2**2*m3*PvR*sman - 2*Bsm2R*m2**2*m3*PvR*sman - 4*BvgL*m2**2*m3*PvR*sman - 
     &     3*BsgR*m3**2*PvR*sman - 3*Bvm1L*m2**2*m3**2*PvR*sman - Bvm2L*m2**2*m3**2*PvR*sman - 6*Bsm1R*m3**3*PvR*sman - 3*BvgL*m3**3*PvR*sman - 
     &     3*Bvm1L*m3**4*PvR*sman - BsgR*m4**2*PvR*sman + 3*Bvm1L*m2**2*m4**2*PvR*sman + Bvm2L*m2**2*m4**2*PvR*sman - 2*Bsm1R*m3*m4**2*PvR*sman - 
     &     BvgL*m3*m4**2*PvR*sman + 2*Bvm1L*m3**2*m4**2*PvR*sman + Bvm1L*m4**4*PvR*sman + 2*BvgR*m2**2*mq1*PvR*sman + 6*Bvm1R*m2**2*m3*mq1*PvR*sman + 
     &     2*Bvm2R*m2**2*m3*mq1*PvR*sman + 3*BvgR*m3**2*mq1*PvR*sman + 6*Bvm1R*m3**3*mq1*PvR*sman + BvgR*m4**2*mq1*PvR*sman + 2*Bvm1R*m3*m4**2*mq1*PvR*sman + 
     &     BvgR*PsL*sman**2 + 2*Bvm1R*m3*PsL*sman**2 - Bvm1L*mq1*PsL*sman**2 + BsgR*PvR*sman**2 - 3*Bvm1L*m2**2*PvR*sman**2 - Bvm2L*m2**2*PvR*sman**2 + 
     &     2*Bsm1R*m3*PvR*sman**2 + BvgL*m3*PvR*sman**2 - 2*Bvm1L*m3**2*PvR*sman**2 - 2*Bvm1L*m4**2*PvR*sman**2 - BvgR*mq1*PvR*sman**2 - 
     &     2*Bvm1R*m3*mq1*PvR*sman**2 + Bvm1L*PvR*sman**3 + Bsm2L*m3**2*PsL*tman - Bsm2L*m4**2*PsL*tman - 2*BvgL*m3*mq1*PsL*tman - Bvm1L*m3**2*mq1*PsL*tman - 
     &     Bvm2L*m3**2*mq1*PsL*tman + Bvm1L*m4**2*mq1*PsL*tman + Bvm2L*m4**2*mq1*PsL*tman + Bsm2L*PsL*sman*tman + 2*Bvm1R*m3*PsL*sman*tman + 
     &     2*Bvm2R*m3*PsL*sman*tman - Bvm1L*mq1*PsL*sman*tman - Bvm2L*mq1*PsL*sman*tman + 2*Bsm1R*m3*PvR*sman*tman + 2*Bsm2R*m3*PvR*sman*tman + 
     &     2*BvgL*m3*PvR*sman*tman + Bvm1L*m3**2*PvR*sman*tman + Bvm2L*m3**2*PvR*sman*tman - Bvm1L*m4**2*PvR*sman*tman - Bvm2L*m4**2*PvR*sman*tman - 
     &     2*Bvm1R*m3*mq1*PvR*sman*tman - 2*Bvm2R*m3*mq1*PvR*sman*tman + Bvm1L*PvR*sman**2*tman + Bvm2L*PvR*sman**2*tman - 2*BvgL*m3*mq1*PsL*uman - 
     &     2*Bvm1L*m3**2*mq1*PsL*uman + 2*Bvm1L*m4**2*mq1*PsL*uman + 2*BvgR*PsL*sman*uman + 4*Bvm1R*m3*PsL*sman*uman - 2*Bvm1L*mq1*PsL*sman*uman + 
     &     2*BsgR*PvR*sman*uman + 4*Bsm1R*m3*PvR*sman*uman + 2*BvgL*m3*PvR*sman*uman + 2*Bvm1L*m3**2*PvR*sman*uman - 2*Bvm1L*m4**2*PvR*sman*uman - 
     &     2*BvgR*mq1*PvR*sman*uman - 4*Bvm1R*m3*mq1*PvR*sman*uman + 2*Bvm1L*PvR*sman**2*uman + Bsm2L*m2**2*m3**2*m4**2*PsL*xv - Bsm2L*m2**2*m4**4*PsL*xv + 
     &     BvgL*m2**2*m3**3*mq1*PsL*xv + Bvm1L*m2**2*m3**4*mq1*PsL*xv - BvgL*m2**2*m3*m4**2*mq1*PsL*xv - Bvm1L*m2**2*m3**2*m4**2*mq1*PsL*xv - 
     &     Bvm2L*m2**2*m3**2*m4**2*mq1*PsL*xv + Bvm2L*m2**2*m4**4*mq1*PsL*xv - BvgR*m2**2*m3**2*PsL*sman*xv - 2*Bvm1R*m2**2*m3**3*PsL*sman*xv + 
     &     Bsm2L*m2**2*m4**2*PsL*sman*xv - BvgR*m2**2*m4**2*PsL*sman*xv + 2*Bvm2R*m2**2*m3*m4**2*PsL*sman*xv - BvgL*m2**2*m3*mq1*PsL*sman*xv + 
     &     Bvm1L*m2**2*m4**2*mq1*PsL*sman*xv - Bvm2L*m2**2*m4**2*mq1*PsL*sman*xv - BsgR*m2**2*m3**2*PvR*sman*xv - 2*Bsm1R*m2**2*m3**3*PvR*sman*xv - 
     &     BvgL*m2**2*m3**3*PvR*sman*xv - Bvm1L*m2**2*m3**4*PvR*sman*xv - BsgR*m2**2*m4**2*PvR*sman*xv + 2*Bsm2R*m2**2*m3*m4**2*PvR*sman*xv + 
     &     BvgL*m2**2*m3*m4**2*PvR*sman*xv + Bvm1L*m2**2*m3**2*m4**2*PvR*sman*xv + Bvm2L*m2**2*m3**2*m4**2*PvR*sman*xv - Bvm2L*m2**2*m4**4*PvR*sman*xv + 
     &     BvgR*m2**2*m3**2*mq1*PvR*sman*xv + 2*Bvm1R*m2**2*m3**3*mq1*PvR*sman*xv + BvgR*m2**2*m4**2*mq1*PvR*sman*xv - 
     &     2*Bvm2R*m2**2*m3*m4**2*mq1*PvR*sman*xv + BvgR*m2**2*PsL*sman**2*xv + 2*Bvm1R*m2**2*m3*PsL*sman**2*xv - Bvm1L*m2**2*mq1*PsL*sman**2*xv + 
     &     BsgR*m2**2*PvR*sman**2*xv + 2*Bsm1R*m2**2*m3*PvR*sman**2*xv + BvgL*m2**2*m3*PvR*sman**2*xv - Bvm1L*m2**2*m4**2*PvR*sman**2*xv + 
     &     Bvm2L*m2**2*m4**2*PvR*sman**2*xv - BvgR*m2**2*mq1*PvR*sman**2*xv - 2*Bvm1R*m2**2*m3*mq1*PvR*sman**2*xv + Bvm1L*m2**2*PvR*sman**3*xv - 
     &     Bsm2L*m3**2*m4**2*PsL*tman*xv + Bsm2L*m4**4*PsL*tman*xv - BvgL*m3**3*mq1*PsL*tman*xv - Bvm1L*m3**4*mq1*PsL*tman*xv + 
     &     BvgL*m3*m4**2*mq1*PsL*tman*xv + Bvm1L*m3**2*m4**2*mq1*PsL*tman*xv + Bvm2L*m3**2*m4**2*mq1*PsL*tman*xv - Bvm2L*m4**4*mq1*PsL*tman*xv + 
     &     BvgR*m3**2*PsL*sman*tman*xv + 2*Bvm1R*m3**3*PsL*sman*tman*xv - Bsm2L*m4**2*PsL*sman*tman*xv + BvgR*m4**2*PsL*sman*tman*xv - 
     &     2*Bvm2R*m3*m4**2*PsL*sman*tman*xv + BvgL*m3*mq1*PsL*sman*tman*xv - Bvm1L*m4**2*mq1*PsL*sman*tman*xv + Bvm2L*m4**2*mq1*PsL*sman*tman*xv + 
     &     BsgR*m3**2*PvR*sman*tman*xv + 2*Bsm1R*m3**3*PvR*sman*tman*xv + BvgL*m3**3*PvR*sman*tman*xv + Bvm1L*m3**4*PvR*sman*tman*xv + 
     &     BsgR*m4**2*PvR*sman*tman*xv - 2*Bsm2R*m3*m4**2*PvR*sman*tman*xv - BvgL*m3*m4**2*PvR*sman*tman*xv - Bvm1L*m3**2*m4**2*PvR*sman*tman*xv - 
     &     Bvm2L*m3**2*m4**2*PvR*sman*tman*xv + Bvm2L*m4**4*PvR*sman*tman*xv - BvgR*m3**2*mq1*PvR*sman*tman*xv - 2*Bvm1R*m3**3*mq1*PvR*sman*tman*xv - 
     &     BvgR*m4**2*mq1*PvR*sman*tman*xv + 2*Bvm2R*m3*m4**2*mq1*PvR*sman*tman*xv - BvgR*PsL*sman**2*tman*xv - 2*Bvm1R*m3*PsL*sman**2*tman*xv + 
     &     Bvm1L*mq1*PsL*sman**2*tman*xv - BsgR*PvR*sman**2*tman*xv - 2*Bsm1R*m3*PvR*sman**2*tman*xv - BvgL*m3*PvR*sman**2*tman*xv + 
     &     Bvm1L*m4**2*PvR*sman**2*tman*xv - Bvm2L*m4**2*PvR*sman**2*tman*xv + BvgR*mq1*PvR*sman**2*tman*xv + 2*Bvm1R*m3*mq1*PvR*sman**2*tman*xv - 
     &     Bvm1L*PvR*sman**3*tman*xv - Bsm1L*PsL*(m3**2 - m4**2 + sman)*
     &      (m4**2 - sman - tman - 2*uman + sman*tman*xv + m2**2*(3 + m3**2*xv - sman*xv) + m3**2*(3 - tman*xv)) - 
     &     BsgL*m3*PsL*(m4**2 - sman - 2*tman - 2*uman + m4**2*tman*xv + sman*tman*xv + m2**2*(4 + m3**2*xv - m4**2*xv - sman*xv) + m3**2*(3 - tman*xv))) + 
     &  2*A1vL*A2Rc*CLc*m1*(-(Bsm2R*m2**2*m3**2*PsR) + Bsm2R*m2**2*m4**2*PsR + 4*BvgR*m2**2*m3*mq1*PsR + 3*Bvm1R*m2**2*m3**2*mq1*PsR + 
     &     Bvm2R*m2**2*m3**2*mq1*PsR + 3*BvgR*m3**3*mq1*PsR + 3*Bvm1R*m3**4*mq1*PsR - 3*Bvm1R*m2**2*m4**2*mq1*PsR - Bvm2R*m2**2*m4**2*mq1*PsR + 
     &     BvgR*m3*m4**2*mq1*PsR - 2*Bvm1R*m3**2*m4**2*mq1*PsR - Bvm1R*m4**4*mq1*PsR - Bsm2R*m2**2*PsR*sman - 2*BvgL*m2**2*PsR*sman - 
     &     6*Bvm1L*m2**2*m3*PsR*sman - 2*Bvm2L*m2**2*m3*PsR*sman - 3*BvgL*m3**2*PsR*sman - 6*Bvm1L*m3**3*PsR*sman - BvgL*m4**2*PsR*sman - 
     &     2*Bvm1L*m3*m4**2*PsR*sman + 3*Bvm1R*m2**2*mq1*PsR*sman + Bvm2R*m2**2*mq1*PsR*sman - BvgR*m3*mq1*PsR*sman + 2*Bvm1R*m3**2*mq1*PsR*sman + 
     &     2*Bvm1R*m4**2*mq1*PsR*sman - 2*BsgL*m2**2*PvL*sman - 6*Bsm1L*m2**2*m3*PvL*sman - 2*Bsm2L*m2**2*m3*PvL*sman - 4*BvgR*m2**2*m3*PvL*sman - 
     &     3*BsgL*m3**2*PvL*sman - 3*Bvm1R*m2**2*m3**2*PvL*sman - Bvm2R*m2**2*m3**2*PvL*sman - 6*Bsm1L*m3**3*PvL*sman - 3*BvgR*m3**3*PvL*sman - 
     &     3*Bvm1R*m3**4*PvL*sman - BsgL*m4**2*PvL*sman + 3*Bvm1R*m2**2*m4**2*PvL*sman + Bvm2R*m2**2*m4**2*PvL*sman - 2*Bsm1L*m3*m4**2*PvL*sman - 
     &     BvgR*m3*m4**2*PvL*sman + 2*Bvm1R*m3**2*m4**2*PvL*sman + Bvm1R*m4**4*PvL*sman + 2*BvgL*m2**2*mq1*PvL*sman + 6*Bvm1L*m2**2*m3*mq1*PvL*sman + 
     &     2*Bvm2L*m2**2*m3*mq1*PvL*sman + 3*BvgL*m3**2*mq1*PvL*sman + 6*Bvm1L*m3**3*mq1*PvL*sman + BvgL*m4**2*mq1*PvL*sman + 2*Bvm1L*m3*m4**2*mq1*PvL*sman + 
     &     BvgL*PsR*sman**2 + 2*Bvm1L*m3*PsR*sman**2 - Bvm1R*mq1*PsR*sman**2 + BsgL*PvL*sman**2 - 3*Bvm1R*m2**2*PvL*sman**2 - Bvm2R*m2**2*PvL*sman**2 + 
     &     2*Bsm1L*m3*PvL*sman**2 + BvgR*m3*PvL*sman**2 - 2*Bvm1R*m3**2*PvL*sman**2 - 2*Bvm1R*m4**2*PvL*sman**2 - BvgL*mq1*PvL*sman**2 - 
     &     2*Bvm1L*m3*mq1*PvL*sman**2 + Bvm1R*PvL*sman**3 + Bsm2R*m3**2*PsR*tman - Bsm2R*m4**2*PsR*tman - 2*BvgR*m3*mq1*PsR*tman - Bvm1R*m3**2*mq1*PsR*tman - 
     &     Bvm2R*m3**2*mq1*PsR*tman + Bvm1R*m4**2*mq1*PsR*tman + Bvm2R*m4**2*mq1*PsR*tman + Bsm2R*PsR*sman*tman + 2*Bvm1L*m3*PsR*sman*tman + 
     &     2*Bvm2L*m3*PsR*sman*tman - Bvm1R*mq1*PsR*sman*tman - Bvm2R*mq1*PsR*sman*tman + 2*Bsm1L*m3*PvL*sman*tman + 2*Bsm2L*m3*PvL*sman*tman + 
     &     2*BvgR*m3*PvL*sman*tman + Bvm1R*m3**2*PvL*sman*tman + Bvm2R*m3**2*PvL*sman*tman - Bvm1R*m4**2*PvL*sman*tman - Bvm2R*m4**2*PvL*sman*tman - 
     &     2*Bvm1L*m3*mq1*PvL*sman*tman - 2*Bvm2L*m3*mq1*PvL*sman*tman + Bvm1R*PvL*sman**2*tman + Bvm2R*PvL*sman**2*tman - 2*BvgR*m3*mq1*PsR*uman - 
     &     2*Bvm1R*m3**2*mq1*PsR*uman + 2*Bvm1R*m4**2*mq1*PsR*uman + 2*BvgL*PsR*sman*uman + 4*Bvm1L*m3*PsR*sman*uman - 2*Bvm1R*mq1*PsR*sman*uman + 
     &     2*BsgL*PvL*sman*uman + 4*Bsm1L*m3*PvL*sman*uman + 2*BvgR*m3*PvL*sman*uman + 2*Bvm1R*m3**2*PvL*sman*uman - 2*Bvm1R*m4**2*PvL*sman*uman - 
     &     2*BvgL*mq1*PvL*sman*uman - 4*Bvm1L*m3*mq1*PvL*sman*uman + 2*Bvm1R*PvL*sman**2*uman + Bsm2R*m2**2*m3**2*m4**2*PsR*xv - Bsm2R*m2**2*m4**4*PsR*xv + 
     &     BvgR*m2**2*m3**3*mq1*PsR*xv + Bvm1R*m2**2*m3**4*mq1*PsR*xv - BvgR*m2**2*m3*m4**2*mq1*PsR*xv - Bvm1R*m2**2*m3**2*m4**2*mq1*PsR*xv - 
     &     Bvm2R*m2**2*m3**2*m4**2*mq1*PsR*xv + Bvm2R*m2**2*m4**4*mq1*PsR*xv - BvgL*m2**2*m3**2*PsR*sman*xv - 2*Bvm1L*m2**2*m3**3*PsR*sman*xv + 
     &     Bsm2R*m2**2*m4**2*PsR*sman*xv - BvgL*m2**2*m4**2*PsR*sman*xv + 2*Bvm2L*m2**2*m3*m4**2*PsR*sman*xv - BvgR*m2**2*m3*mq1*PsR*sman*xv + 
     &     Bvm1R*m2**2*m4**2*mq1*PsR*sman*xv - Bvm2R*m2**2*m4**2*mq1*PsR*sman*xv - BsgL*m2**2*m3**2*PvL*sman*xv - 2*Bsm1L*m2**2*m3**3*PvL*sman*xv - 
     &     BvgR*m2**2*m3**3*PvL*sman*xv - Bvm1R*m2**2*m3**4*PvL*sman*xv - BsgL*m2**2*m4**2*PvL*sman*xv + 2*Bsm2L*m2**2*m3*m4**2*PvL*sman*xv + 
     &     BvgR*m2**2*m3*m4**2*PvL*sman*xv + Bvm1R*m2**2*m3**2*m4**2*PvL*sman*xv + Bvm2R*m2**2*m3**2*m4**2*PvL*sman*xv - Bvm2R*m2**2*m4**4*PvL*sman*xv + 
     &     BvgL*m2**2*m3**2*mq1*PvL*sman*xv + 2*Bvm1L*m2**2*m3**3*mq1*PvL*sman*xv + BvgL*m2**2*m4**2*mq1*PvL*sman*xv - 
     &     2*Bvm2L*m2**2*m3*m4**2*mq1*PvL*sman*xv + BvgL*m2**2*PsR*sman**2*xv + 2*Bvm1L*m2**2*m3*PsR*sman**2*xv - Bvm1R*m2**2*mq1*PsR*sman**2*xv + 
     &     BsgL*m2**2*PvL*sman**2*xv + 2*Bsm1L*m2**2*m3*PvL*sman**2*xv + BvgR*m2**2*m3*PvL*sman**2*xv - Bvm1R*m2**2*m4**2*PvL*sman**2*xv + 
     &     Bvm2R*m2**2*m4**2*PvL*sman**2*xv - BvgL*m2**2*mq1*PvL*sman**2*xv - 2*Bvm1L*m2**2*m3*mq1*PvL*sman**2*xv + Bvm1R*m2**2*PvL*sman**3*xv - 
     &     Bsm2R*m3**2*m4**2*PsR*tman*xv + Bsm2R*m4**4*PsR*tman*xv - BvgR*m3**3*mq1*PsR*tman*xv - Bvm1R*m3**4*mq1*PsR*tman*xv + 
     &     BvgR*m3*m4**2*mq1*PsR*tman*xv + Bvm1R*m3**2*m4**2*mq1*PsR*tman*xv + Bvm2R*m3**2*m4**2*mq1*PsR*tman*xv - Bvm2R*m4**4*mq1*PsR*tman*xv + 
     &     BvgL*m3**2*PsR*sman*tman*xv + 2*Bvm1L*m3**3*PsR*sman*tman*xv - Bsm2R*m4**2*PsR*sman*tman*xv + BvgL*m4**2*PsR*sman*tman*xv - 
     &     2*Bvm2L*m3*m4**2*PsR*sman*tman*xv + BvgR*m3*mq1*PsR*sman*tman*xv - Bvm1R*m4**2*mq1*PsR*sman*tman*xv + Bvm2R*m4**2*mq1*PsR*sman*tman*xv + 
     &     BsgL*m3**2*PvL*sman*tman*xv + 2*Bsm1L*m3**3*PvL*sman*tman*xv + BvgR*m3**3*PvL*sman*tman*xv + Bvm1R*m3**4*PvL*sman*tman*xv + 
     &     BsgL*m4**2*PvL*sman*tman*xv - 2*Bsm2L*m3*m4**2*PvL*sman*tman*xv - BvgR*m3*m4**2*PvL*sman*tman*xv - Bvm1R*m3**2*m4**2*PvL*sman*tman*xv - 
     &     Bvm2R*m3**2*m4**2*PvL*sman*tman*xv + Bvm2R*m4**4*PvL*sman*tman*xv - BvgL*m3**2*mq1*PvL*sman*tman*xv - 2*Bvm1L*m3**3*mq1*PvL*sman*tman*xv - 
     &     BvgL*m4**2*mq1*PvL*sman*tman*xv + 2*Bvm2L*m3*m4**2*mq1*PvL*sman*tman*xv - BvgL*PsR*sman**2*tman*xv - 2*Bvm1L*m3*PsR*sman**2*tman*xv + 
     &     Bvm1R*mq1*PsR*sman**2*tman*xv - BsgL*PvL*sman**2*tman*xv - 2*Bsm1L*m3*PvL*sman**2*tman*xv - BvgR*m3*PvL*sman**2*tman*xv + 
     &     Bvm1R*m4**2*PvL*sman**2*tman*xv - Bvm2R*m4**2*PvL*sman**2*tman*xv + BvgL*mq1*PvL*sman**2*tman*xv + 2*Bvm1L*m3*mq1*PvL*sman**2*tman*xv - 
     &     Bvm1R*PvL*sman**3*tman*xv - Bsm1R*PsR*(m3**2 - m4**2 + sman)*
     &      (m4**2 - sman - tman - 2*uman + sman*tman*xv + m2**2*(3 + m3**2*xv - sman*xv) + m3**2*(3 - tman*xv)) - 
     &     BsgR*m3*PsR*(m4**2 - sman - 2*tman - 2*uman + m4**2*tman*xv + sman*tman*xv + m2**2*(4 + m3**2*xv - m4**2*xv - sman*xv) + m3**2*(3 - tman*xv))) + 
     &  2*A1sL*A2Rc*CLc*m1*(-2*Bsm2L*m2**2*m3*PsL - 4*BvgR*m2**2*m3*PsL - 3*Bvm1R*m2**2*m3**2*PsL - Bvm2R*m2**2*m3**2*PsL - 3*BvgR*m3**3*PsL - 
     &     3*Bvm1R*m3**4*PsL + 3*Bvm1R*m2**2*m4**2*PsL + Bvm2R*m2**2*m4**2*PsL - BvgR*m3*m4**2*PsL + 2*Bvm1R*m3**2*m4**2*PsL + Bvm1R*m4**4*PsL + 
     &     2*BvgL*m2**2*mq1*PsL + 6*Bvm1L*m2**2*m3*mq1*PsL + 2*Bvm2L*m2**2*m3*mq1*PsL + 3*BvgL*m3**2*mq1*PsL + 6*Bvm1L*m3**3*mq1*PsL + BvgL*m4**2*mq1*PsL + 
     &     2*Bvm1L*m3*m4**2*mq1*PsL - 4*BsgR*m2**2*m3*PvR - 3*Bsm1R*m2**2*m3**2*PvR - Bsm2R*m2**2*m3**2*PvR - 3*BsgR*m3**3*PvR - 3*Bsm1R*m3**4*PvR + 
     &     3*Bsm1R*m2**2*m4**2*PvR + Bsm2R*m2**2*m4**2*PvR - BsgR*m3*m4**2*PvR + 2*Bsm1R*m3**2*m4**2*PvR + Bsm1R*m4**4*PvR + 4*BvgR*m2**2*m3*mq1*PvR + 
     &     3*Bvm1R*m2**2*m3**2*mq1*PvR + Bvm2R*m2**2*m3**2*mq1*PvR + 3*BvgR*m3**3*mq1*PvR + 3*Bvm1R*m3**4*mq1*PvR - 3*Bvm1R*m2**2*m4**2*mq1*PvR - 
     &     Bvm2R*m2**2*m4**2*mq1*PvR + BvgR*m3*m4**2*mq1*PvR - 2*Bvm1R*m3**2*m4**2*mq1*PvR - Bvm1R*m4**4*mq1*PvR - 3*Bvm1R*m2**2*PsL*sman - 
     &     Bvm2R*m2**2*PsL*sman + BvgR*m3*PsL*sman - 2*Bvm1R*m3**2*PsL*sman - 2*Bvm1R*m4**2*PsL*sman - BvgL*mq1*PsL*sman - 2*Bvm1L*m3*mq1*PsL*sman - 
     &     3*Bsm1R*m2**2*PvR*sman - Bsm2R*m2**2*PvR*sman - 2*BvgL*m2**2*PvR*sman + BsgR*m3*PvR*sman - 6*Bvm1L*m2**2*m3*PvR*sman - 2*Bvm2L*m2**2*m3*PvR*sman - 
     &     2*Bsm1R*m3**2*PvR*sman - 3*BvgL*m3**2*PvR*sman - 6*Bvm1L*m3**3*PvR*sman - 2*Bsm1R*m4**2*PvR*sman - BvgL*m4**2*PvR*sman - 
     &     2*Bvm1L*m3*m4**2*PvR*sman + 3*Bvm1R*m2**2*mq1*PvR*sman + Bvm2R*m2**2*mq1*PvR*sman - BvgR*m3*mq1*PvR*sman + 2*Bvm1R*m3**2*mq1*PvR*sman + 
     &     2*Bvm1R*m4**2*mq1*PvR*sman + Bvm1R*PsL*sman**2 + Bsm1R*PvR*sman**2 + BvgL*PvR*sman**2 + 2*Bvm1L*m3*PvR*sman**2 - Bvm1R*mq1*PvR*sman**2 + 
     &     2*Bsm2L*m3*PsL*tman + 2*BvgR*m3*PsL*tman + Bvm1R*m3**2*PsL*tman + Bvm2R*m3**2*PsL*tman - Bvm1R*m4**2*PsL*tman - Bvm2R*m4**2*PsL*tman - 
     &     2*Bvm1L*m3*mq1*PsL*tman - 2*Bvm2L*m3*mq1*PsL*tman + 2*BsgR*m3*PvR*tman + Bsm1R*m3**2*PvR*tman + Bsm2R*m3**2*PvR*tman - Bsm1R*m4**2*PvR*tman - 
     &     Bsm2R*m4**2*PvR*tman - 2*BvgR*m3*mq1*PvR*tman - Bvm1R*m3**2*mq1*PvR*tman - Bvm2R*m3**2*mq1*PvR*tman + Bvm1R*m4**2*mq1*PvR*tman + 
     &     Bvm2R*m4**2*mq1*PvR*tman + Bvm1R*PsL*sman*tman + Bvm2R*PsL*sman*tman + Bsm1R*PvR*sman*tman + Bsm2R*PvR*sman*tman + 2*Bvm1L*m3*PvR*sman*tman + 
     &     2*Bvm2L*m3*PvR*sman*tman - Bvm1R*mq1*PvR*sman*tman - Bvm2R*mq1*PvR*sman*tman + 2*BvgR*m3*PsL*uman + 2*Bvm1R*m3**2*PsL*uman - 
     &     2*Bvm1R*m4**2*PsL*uman - 2*BvgL*mq1*PsL*uman - 4*Bvm1L*m3*mq1*PsL*uman + 2*BsgR*m3*PvR*uman + 2*Bsm1R*m3**2*PvR*uman - 2*Bsm1R*m4**2*PvR*uman - 
     &     2*BvgR*m3*mq1*PvR*uman - 2*Bvm1R*m3**2*mq1*PvR*uman + 2*Bvm1R*m4**2*mq1*PvR*uman + 2*Bvm1R*PsL*sman*uman + 2*Bsm1R*PvR*sman*uman + 
     &     2*BvgL*PvR*sman*uman + 4*Bvm1L*m3*PvR*sman*uman - 2*Bvm1R*mq1*PvR*sman*uman - BvgR*m2**2*m3**3*PsL*xv - Bvm1R*m2**2*m3**4*PsL*xv + 
     &     2*Bsm2L*m2**2*m3*m4**2*PsL*xv + BvgR*m2**2*m3*m4**2*PsL*xv + Bvm1R*m2**2*m3**2*m4**2*PsL*xv + Bvm2R*m2**2*m3**2*m4**2*PsL*xv - 
     &     Bvm2R*m2**2*m4**4*PsL*xv + BvgL*m2**2*m3**2*mq1*PsL*xv + 2*Bvm1L*m2**2*m3**3*mq1*PsL*xv + BvgL*m2**2*m4**2*mq1*PsL*xv - 
     &     2*Bvm2L*m2**2*m3*m4**2*mq1*PsL*xv - BsgR*m2**2*m3**3*PvR*xv - Bsm1R*m2**2*m3**4*PvR*xv + BsgR*m2**2*m3*m4**2*PvR*xv + 
     &     Bsm1R*m2**2*m3**2*m4**2*PvR*xv + Bsm2R*m2**2*m3**2*m4**2*PvR*xv - Bsm2R*m2**2*m4**4*PvR*xv + BvgR*m2**2*m3**3*mq1*PvR*xv + 
     &     Bvm1R*m2**2*m3**4*mq1*PvR*xv - BvgR*m2**2*m3*m4**2*mq1*PvR*xv - Bvm1R*m2**2*m3**2*m4**2*mq1*PvR*xv - Bvm2R*m2**2*m3**2*m4**2*mq1*PvR*xv + 
     &     Bvm2R*m2**2*m4**4*mq1*PvR*xv + BvgR*m2**2*m3*PsL*sman*xv - Bvm1R*m2**2*m4**2*PsL*sman*xv + Bvm2R*m2**2*m4**2*PsL*sman*xv - 
     &     BvgL*m2**2*mq1*PsL*sman*xv - 2*Bvm1L*m2**2*m3*mq1*PsL*sman*xv + BsgR*m2**2*m3*PvR*sman*xv - BvgL*m2**2*m3**2*PvR*sman*xv - 
     &     2*Bvm1L*m2**2*m3**3*PvR*sman*xv - Bsm1R*m2**2*m4**2*PvR*sman*xv + Bsm2R*m2**2*m4**2*PvR*sman*xv - BvgL*m2**2*m4**2*PvR*sman*xv + 
     &     2*Bvm2L*m2**2*m3*m4**2*PvR*sman*xv - BvgR*m2**2*m3*mq1*PvR*sman*xv + Bvm1R*m2**2*m4**2*mq1*PvR*sman*xv - Bvm2R*m2**2*m4**2*mq1*PvR*sman*xv + 
     &     Bvm1R*m2**2*PsL*sman**2*xv + Bsm1R*m2**2*PvR*sman**2*xv + BvgL*m2**2*PvR*sman**2*xv + 2*Bvm1L*m2**2*m3*PvR*sman**2*xv - 
     &     Bvm1R*m2**2*mq1*PvR*sman**2*xv + BvgR*m3**3*PsL*tman*xv + Bvm1R*m3**4*PsL*tman*xv - 2*Bsm2L*m3*m4**2*PsL*tman*xv - BvgR*m3*m4**2*PsL*tman*xv - 
     &     Bvm1R*m3**2*m4**2*PsL*tman*xv - Bvm2R*m3**2*m4**2*PsL*tman*xv + Bvm2R*m4**4*PsL*tman*xv - BvgL*m3**2*mq1*PsL*tman*xv - 
     &     2*Bvm1L*m3**3*mq1*PsL*tman*xv - BvgL*m4**2*mq1*PsL*tman*xv + 2*Bvm2L*m3*m4**2*mq1*PsL*tman*xv + BsgR*m3**3*PvR*tman*xv + Bsm1R*m3**4*PvR*tman*xv - 
     &     BsgR*m3*m4**2*PvR*tman*xv - Bsm1R*m3**2*m4**2*PvR*tman*xv - Bsm2R*m3**2*m4**2*PvR*tman*xv + Bsm2R*m4**4*PvR*tman*xv - BvgR*m3**3*mq1*PvR*tman*xv - 
     &     Bvm1R*m3**4*mq1*PvR*tman*xv + BvgR*m3*m4**2*mq1*PvR*tman*xv + Bvm1R*m3**2*m4**2*mq1*PvR*tman*xv + Bvm2R*m3**2*m4**2*mq1*PvR*tman*xv - 
     &     Bvm2R*m4**4*mq1*PvR*tman*xv - BvgR*m3*PsL*sman*tman*xv + Bvm1R*m4**2*PsL*sman*tman*xv - Bvm2R*m4**2*PsL*sman*tman*xv + BvgL*mq1*PsL*sman*tman*xv + 
     &     2*Bvm1L*m3*mq1*PsL*sman*tman*xv - BsgR*m3*PvR*sman*tman*xv + BvgL*m3**2*PvR*sman*tman*xv + 2*Bvm1L*m3**3*PvR*sman*tman*xv + 
     &     Bsm1R*m4**2*PvR*sman*tman*xv - Bsm2R*m4**2*PvR*sman*tman*xv + BvgL*m4**2*PvR*sman*tman*xv - 2*Bvm2L*m3*m4**2*PvR*sman*tman*xv + 
     &     BvgR*m3*mq1*PvR*sman*tman*xv - Bvm1R*m4**2*mq1*PvR*sman*tman*xv + Bvm2R*m4**2*mq1*PvR*sman*tman*xv - Bvm1R*PsL*sman**2*tman*xv - 
     &     Bsm1R*PvR*sman**2*tman*xv - BvgL*PvR*sman**2*tman*xv - 2*Bvm1L*m3*PvR*sman**2*tman*xv + Bvm1R*mq1*PvR*sman**2*tman*xv - 
     &     2*Bsm1L*m3*PsL*(m4**2 - sman - tman - 2*uman + sman*tman*xv + m2**2*(3 + m3**2*xv - sman*xv) + m3**2*(3 - tman*xv)) + 
     &     BsgL*PsL*(-m4**2 + sman + 2*uman + m4**2*tman*xv - sman*tman*xv - m2**2*(2 + m3**2*xv + m4**2*xv - sman*xv) + m3**2*(-3 + tman*xv))) + 
     &  2*A1sR*A2Lc*CLc*m1*(-2*Bsm2R*m2**2*m3*PsR - 4*BvgL*m2**2*m3*PsR - 3*Bvm1L*m2**2*m3**2*PsR - Bvm2L*m2**2*m3**2*PsR - 3*BvgL*m3**3*PsR - 
     &     3*Bvm1L*m3**4*PsR + 3*Bvm1L*m2**2*m4**2*PsR + Bvm2L*m2**2*m4**2*PsR - BvgL*m3*m4**2*PsR + 2*Bvm1L*m3**2*m4**2*PsR + Bvm1L*m4**4*PsR + 
     &     2*BvgR*m2**2*mq1*PsR + 6*Bvm1R*m2**2*m3*mq1*PsR + 2*Bvm2R*m2**2*m3*mq1*PsR + 3*BvgR*m3**2*mq1*PsR + 6*Bvm1R*m3**3*mq1*PsR + BvgR*m4**2*mq1*PsR + 
     &     2*Bvm1R*m3*m4**2*mq1*PsR - 4*BsgL*m2**2*m3*PvL - 3*Bsm1L*m2**2*m3**2*PvL - Bsm2L*m2**2*m3**2*PvL - 3*BsgL*m3**3*PvL - 3*Bsm1L*m3**4*PvL + 
     &     3*Bsm1L*m2**2*m4**2*PvL + Bsm2L*m2**2*m4**2*PvL - BsgL*m3*m4**2*PvL + 2*Bsm1L*m3**2*m4**2*PvL + Bsm1L*m4**4*PvL + 4*BvgL*m2**2*m3*mq1*PvL + 
     &     3*Bvm1L*m2**2*m3**2*mq1*PvL + Bvm2L*m2**2*m3**2*mq1*PvL + 3*BvgL*m3**3*mq1*PvL + 3*Bvm1L*m3**4*mq1*PvL - 3*Bvm1L*m2**2*m4**2*mq1*PvL - 
     &     Bvm2L*m2**2*m4**2*mq1*PvL + BvgL*m3*m4**2*mq1*PvL - 2*Bvm1L*m3**2*m4**2*mq1*PvL - Bvm1L*m4**4*mq1*PvL - 3*Bvm1L*m2**2*PsR*sman - 
     &     Bvm2L*m2**2*PsR*sman + BvgL*m3*PsR*sman - 2*Bvm1L*m3**2*PsR*sman - 2*Bvm1L*m4**2*PsR*sman - BvgR*mq1*PsR*sman - 2*Bvm1R*m3*mq1*PsR*sman - 
     &     3*Bsm1L*m2**2*PvL*sman - Bsm2L*m2**2*PvL*sman - 2*BvgR*m2**2*PvL*sman + BsgL*m3*PvL*sman - 6*Bvm1R*m2**2*m3*PvL*sman - 2*Bvm2R*m2**2*m3*PvL*sman - 
     &     2*Bsm1L*m3**2*PvL*sman - 3*BvgR*m3**2*PvL*sman - 6*Bvm1R*m3**3*PvL*sman - 2*Bsm1L*m4**2*PvL*sman - BvgR*m4**2*PvL*sman - 
     &     2*Bvm1R*m3*m4**2*PvL*sman + 3*Bvm1L*m2**2*mq1*PvL*sman + Bvm2L*m2**2*mq1*PvL*sman - BvgL*m3*mq1*PvL*sman + 2*Bvm1L*m3**2*mq1*PvL*sman + 
     &     2*Bvm1L*m4**2*mq1*PvL*sman + Bvm1L*PsR*sman**2 + Bsm1L*PvL*sman**2 + BvgR*PvL*sman**2 + 2*Bvm1R*m3*PvL*sman**2 - Bvm1L*mq1*PvL*sman**2 + 
     &     2*Bsm2R*m3*PsR*tman + 2*BvgL*m3*PsR*tman + Bvm1L*m3**2*PsR*tman + Bvm2L*m3**2*PsR*tman - Bvm1L*m4**2*PsR*tman - Bvm2L*m4**2*PsR*tman - 
     &     2*Bvm1R*m3*mq1*PsR*tman - 2*Bvm2R*m3*mq1*PsR*tman + 2*BsgL*m3*PvL*tman + Bsm1L*m3**2*PvL*tman + Bsm2L*m3**2*PvL*tman - Bsm1L*m4**2*PvL*tman - 
     &     Bsm2L*m4**2*PvL*tman - 2*BvgL*m3*mq1*PvL*tman - Bvm1L*m3**2*mq1*PvL*tman - Bvm2L*m3**2*mq1*PvL*tman + Bvm1L*m4**2*mq1*PvL*tman + 
     &     Bvm2L*m4**2*mq1*PvL*tman + Bvm1L*PsR*sman*tman + Bvm2L*PsR*sman*tman + Bsm1L*PvL*sman*tman + Bsm2L*PvL*sman*tman + 2*Bvm1R*m3*PvL*sman*tman + 
     &     2*Bvm2R*m3*PvL*sman*tman - Bvm1L*mq1*PvL*sman*tman - Bvm2L*mq1*PvL*sman*tman + 2*BvgL*m3*PsR*uman + 2*Bvm1L*m3**2*PsR*uman - 
     &     2*Bvm1L*m4**2*PsR*uman - 2*BvgR*mq1*PsR*uman - 4*Bvm1R*m3*mq1*PsR*uman + 2*BsgL*m3*PvL*uman + 2*Bsm1L*m3**2*PvL*uman - 2*Bsm1L*m4**2*PvL*uman - 
     &     2*BvgL*m3*mq1*PvL*uman - 2*Bvm1L*m3**2*mq1*PvL*uman + 2*Bvm1L*m4**2*mq1*PvL*uman + 2*Bvm1L*PsR*sman*uman + 2*Bsm1L*PvL*sman*uman + 
     &     2*BvgR*PvL*sman*uman + 4*Bvm1R*m3*PvL*sman*uman - 2*Bvm1L*mq1*PvL*sman*uman - BvgL*m2**2*m3**3*PsR*xv - Bvm1L*m2**2*m3**4*PsR*xv + 
     &     2*Bsm2R*m2**2*m3*m4**2*PsR*xv + BvgL*m2**2*m3*m4**2*PsR*xv + Bvm1L*m2**2*m3**2*m4**2*PsR*xv + Bvm2L*m2**2*m3**2*m4**2*PsR*xv - 
     &     Bvm2L*m2**2*m4**4*PsR*xv + BvgR*m2**2*m3**2*mq1*PsR*xv + 2*Bvm1R*m2**2*m3**3*mq1*PsR*xv + BvgR*m2**2*m4**2*mq1*PsR*xv - 
     &     2*Bvm2R*m2**2*m3*m4**2*mq1*PsR*xv - BsgL*m2**2*m3**3*PvL*xv - Bsm1L*m2**2*m3**4*PvL*xv + BsgL*m2**2*m3*m4**2*PvL*xv + 
     &     Bsm1L*m2**2*m3**2*m4**2*PvL*xv + Bsm2L*m2**2*m3**2*m4**2*PvL*xv - Bsm2L*m2**2*m4**4*PvL*xv + BvgL*m2**2*m3**3*mq1*PvL*xv + 
     &     Bvm1L*m2**2*m3**4*mq1*PvL*xv - BvgL*m2**2*m3*m4**2*mq1*PvL*xv - Bvm1L*m2**2*m3**2*m4**2*mq1*PvL*xv - Bvm2L*m2**2*m3**2*m4**2*mq1*PvL*xv + 
     &     Bvm2L*m2**2*m4**4*mq1*PvL*xv + BvgL*m2**2*m3*PsR*sman*xv - Bvm1L*m2**2*m4**2*PsR*sman*xv + Bvm2L*m2**2*m4**2*PsR*sman*xv - 
     &     BvgR*m2**2*mq1*PsR*sman*xv - 2*Bvm1R*m2**2*m3*mq1*PsR*sman*xv + BsgL*m2**2*m3*PvL*sman*xv - BvgR*m2**2*m3**2*PvL*sman*xv - 
     &     2*Bvm1R*m2**2*m3**3*PvL*sman*xv - Bsm1L*m2**2*m4**2*PvL*sman*xv + Bsm2L*m2**2*m4**2*PvL*sman*xv - BvgR*m2**2*m4**2*PvL*sman*xv + 
     &     2*Bvm2R*m2**2*m3*m4**2*PvL*sman*xv - BvgL*m2**2*m3*mq1*PvL*sman*xv + Bvm1L*m2**2*m4**2*mq1*PvL*sman*xv - Bvm2L*m2**2*m4**2*mq1*PvL*sman*xv + 
     &     Bvm1L*m2**2*PsR*sman**2*xv + Bsm1L*m2**2*PvL*sman**2*xv + BvgR*m2**2*PvL*sman**2*xv + 2*Bvm1R*m2**2*m3*PvL*sman**2*xv - 
     &     Bvm1L*m2**2*mq1*PvL*sman**2*xv + BvgL*m3**3*PsR*tman*xv + Bvm1L*m3**4*PsR*tman*xv - 2*Bsm2R*m3*m4**2*PsR*tman*xv - BvgL*m3*m4**2*PsR*tman*xv - 
     &     Bvm1L*m3**2*m4**2*PsR*tman*xv - Bvm2L*m3**2*m4**2*PsR*tman*xv + Bvm2L*m4**4*PsR*tman*xv - BvgR*m3**2*mq1*PsR*tman*xv - 
     &     2*Bvm1R*m3**3*mq1*PsR*tman*xv - BvgR*m4**2*mq1*PsR*tman*xv + 2*Bvm2R*m3*m4**2*mq1*PsR*tman*xv + BsgL*m3**3*PvL*tman*xv + Bsm1L*m3**4*PvL*tman*xv - 
     &     BsgL*m3*m4**2*PvL*tman*xv - Bsm1L*m3**2*m4**2*PvL*tman*xv - Bsm2L*m3**2*m4**2*PvL*tman*xv + Bsm2L*m4**4*PvL*tman*xv - BvgL*m3**3*mq1*PvL*tman*xv - 
     &     Bvm1L*m3**4*mq1*PvL*tman*xv + BvgL*m3*m4**2*mq1*PvL*tman*xv + Bvm1L*m3**2*m4**2*mq1*PvL*tman*xv + Bvm2L*m3**2*m4**2*mq1*PvL*tman*xv - 
     &     Bvm2L*m4**4*mq1*PvL*tman*xv - BvgL*m3*PsR*sman*tman*xv + Bvm1L*m4**2*PsR*sman*tman*xv - Bvm2L*m4**2*PsR*sman*tman*xv + BvgR*mq1*PsR*sman*tman*xv + 
     &     2*Bvm1R*m3*mq1*PsR*sman*tman*xv - BsgL*m3*PvL*sman*tman*xv + BvgR*m3**2*PvL*sman*tman*xv + 2*Bvm1R*m3**3*PvL*sman*tman*xv + 
     &     Bsm1L*m4**2*PvL*sman*tman*xv - Bsm2L*m4**2*PvL*sman*tman*xv + BvgR*m4**2*PvL*sman*tman*xv - 2*Bvm2R*m3*m4**2*PvL*sman*tman*xv + 
     &     BvgL*m3*mq1*PvL*sman*tman*xv - Bvm1L*m4**2*mq1*PvL*sman*tman*xv + Bvm2L*m4**2*mq1*PvL*sman*tman*xv - Bvm1L*PsR*sman**2*tman*xv - 
     &     Bsm1L*PvL*sman**2*tman*xv - BvgR*PvL*sman**2*tman*xv - 2*Bvm1R*m3*PvL*sman**2*tman*xv + Bvm1L*mq1*PvL*sman**2*tman*xv - 
     &     2*Bsm1R*m3*PsR*(m4**2 - sman - tman - 2*uman + sman*tman*xv + m2**2*(3 + m3**2*xv - sman*xv) + m3**2*(3 - tman*xv)) + 
     &     BsgR*PsR*(-m4**2 + sman + 2*uman + m4**2*tman*xv - sman*tman*xv - m2**2*(2 + m3**2*xv + m4**2*xv - sman*xv) + m3**2*(-3 + tman*xv))) - 
     &  2*A1sL*A2Lc*CLc*(4*BvgR*m1**2*m2**2*PsL - 3*BsgL*m1**2*m3*PsL - 2*BsgL*m2**2*m3*PsL + 6*Bvm1R*m1**2*m2**2*m3*PsL + 2*Bvm2R*m1**2*m2**2*m3*PsL + 
     &     6*BvgR*m1**2*m3**2*PsL + 4*BvgR*m2**2*m3**2*PsL + 6*Bvm1R*m1**2*m3**3*PsL + 3*Bvm1R*m2**2*m3**3*PsL + Bvm2R*m2**2*m3**3*PsL + 3*BvgR*m3**4*PsL + 
     &     3*Bvm1R*m3**5*PsL - BsgL*m3*m4**2*PsL + 2*Bvm1R*m1**2*m3*m4**2*PsL + 3*Bvm1R*m2**2*m3*m4**2*PsL + Bvm2R*m2**2*m3*m4**2*PsL + 
     &     3*BvgR*m3**2*m4**2*PsL + 4*Bvm1R*m3**3*m4**2*PsL + Bvm1R*m3*m4**4*PsL - 3*Bvm1L*m1**2*m2**2*mq1*PsL - Bvm2L*m1**2*m2**2*mq1*PsL + 
     &     3*BvgL*m1**2*m3*mq1*PsL + 2*BvgL*m2**2*m3*mq1*PsL - 3*Bvm1L*m1**2*m3**2*mq1*PsL - 3*Bvm1L*m2**2*m3**2*mq1*PsL - Bvm2L*m2**2*m3**2*mq1*PsL - 
     &     3*Bvm1L*m3**4*mq1*PsL - Bvm1L*m1**2*m4**2*mq1*PsL + BvgL*m3*m4**2*mq1*PsL - Bvm1L*m3**2*m4**2*mq1*PsL + 4*BsgR*m1**2*m2**2*PvR + 
     &     6*Bsm1R*m1**2*m2**2*m3*PvR + 2*Bsm2R*m1**2*m2**2*m3*PvR + 6*BsgR*m1**2*m3**2*PvR + 4*BsgR*m2**2*m3**2*PvR + 6*Bsm1R*m1**2*m3**3*PvR + 
     &     3*Bsm1R*m2**2*m3**3*PvR + Bsm2R*m2**2*m3**3*PvR + 3*BsgR*m3**4*PvR + 3*Bsm1R*m3**5*PvR + 2*Bsm1R*m1**2*m3*m4**2*PvR + 3*Bsm1R*m2**2*m3*m4**2*PvR + 
     &     Bsm2R*m2**2*m3*m4**2*PvR + 3*BsgR*m3**2*m4**2*PvR + 4*Bsm1R*m3**3*m4**2*PvR + Bsm1R*m3*m4**4*PvR - 4*BvgR*m1**2*m2**2*mq1*PvR - 
     &     6*Bvm1R*m1**2*m2**2*m3*mq1*PvR - 2*Bvm2R*m1**2*m2**2*m3*mq1*PvR - 6*BvgR*m1**2*m3**2*mq1*PvR - 4*BvgR*m2**2*m3**2*mq1*PvR - 
     &     6*Bvm1R*m1**2*m3**3*mq1*PvR - 3*Bvm1R*m2**2*m3**3*mq1*PvR - Bvm2R*m2**2*m3**3*mq1*PvR - 3*BvgR*m3**4*mq1*PvR - 3*Bvm1R*m3**5*mq1*PvR - 
     &     2*Bvm1R*m1**2*m3*m4**2*mq1*PvR - 3*Bvm1R*m2**2*m3*m4**2*mq1*PvR - Bvm2R*m2**2*m3*m4**2*mq1*PvR - 3*BvgR*m3**2*m4**2*mq1*PvR - 
     &     4*Bvm1R*m3**3*m4**2*mq1*PvR - Bvm1R*m3*m4**4*mq1*PvR + BvgR*m2**2*PsL*sman + 2*BsgL*m3*PsL*sman - 2*Bvm1R*m1**2*m3*PsL*sman - 
     &     2*BvgR*m3**2*PsL*sman - Bvm1R*m3**3*PsL*sman + BvgR*m4**2*PsL*sman - Bvm1R*m3*m4**2*PsL*sman + Bvm1L*m1**2*mq1*PsL*sman - 2*BvgL*m3*mq1*PsL*sman + 
     &     Bvm1L*m3**2*mq1*PsL*sman + BsgR*m2**2*PvR*sman + 3*Bvm1L*m1**2*m2**2*PvR*sman + Bvm2L*m1**2*m2**2*PvR*sman - 2*Bsm1R*m1**2*m3*PvR*sman - 
     &     3*BvgL*m1**2*m3*PvR*sman - 2*BvgL*m2**2*m3*PvR*sman - 2*BsgR*m3**2*PvR*sman + 3*Bvm1L*m1**2*m3**2*PvR*sman + 3*Bvm1L*m2**2*m3**2*PvR*sman + 
     &     Bvm2L*m2**2*m3**2*PvR*sman - Bsm1R*m3**3*PvR*sman + 3*Bvm1L*m3**4*PvR*sman + BsgR*m4**2*PvR*sman + Bvm1L*m1**2*m4**2*PvR*sman - 
     &     Bsm1R*m3*m4**2*PvR*sman - BvgL*m3*m4**2*PvR*sman + Bvm1L*m3**2*m4**2*PvR*sman - BvgR*m2**2*mq1*PvR*sman + 2*Bvm1R*m1**2*m3*mq1*PvR*sman + 
     &     2*BvgR*m3**2*mq1*PvR*sman + Bvm1R*m3**3*mq1*PvR*sman - BvgR*m4**2*mq1*PvR*sman + Bvm1R*m3*m4**2*mq1*PvR*sman - BvgR*PsL*sman**2 - 
     &     BsgR*PvR*sman**2 - Bvm1L*m1**2*PvR*sman**2 + 2*BvgL*m3*PvR*sman**2 - Bvm1L*m3**2*PvR*sman**2 + BvgR*mq1*PvR*sman**2 - BvgR*m1**2*PsL*tman - 
     &     3*BvgR*m2**2*PsL*tman - 2*Bvm1R*m1**2*m3*PsL*tman - 2*Bvm2R*m1**2*m3*PsL*tman - 3*Bvm1R*m2**2*m3*PsL*tman - Bvm2R*m2**2*m3*PsL*tman - 
     &     4*BvgR*m3**2*PsL*tman - 4*Bvm1R*m3**3*PsL*tman - Bvm2R*m3**3*PsL*tman - BvgR*m4**2*PsL*tman - 2*Bvm1R*m3*m4**2*PsL*tman - Bvm2R*m3*m4**2*PsL*tman + 
     &     Bvm1L*m1**2*mq1*PsL*tman + Bvm2L*m1**2*mq1*PsL*tman + 3*Bvm1L*m2**2*mq1*PsL*tman + Bvm2L*m2**2*mq1*PsL*tman + 4*Bvm1L*m3**2*mq1*PsL*tman + 
     &     Bvm2L*m3**2*mq1*PsL*tman + Bvm1L*m4**2*mq1*PsL*tman - BsgR*m1**2*PvR*tman - 3*BsgR*m2**2*PvR*tman - 2*Bsm1R*m1**2*m3*PvR*tman - 
     &     2*Bsm2R*m1**2*m3*PvR*tman - 3*Bsm1R*m2**2*m3*PvR*tman - Bsm2R*m2**2*m3*PvR*tman - 4*BsgR*m3**2*PvR*tman - 4*Bsm1R*m3**3*PvR*tman - 
     &     Bsm2R*m3**3*PvR*tman - BsgR*m4**2*PvR*tman - 2*Bsm1R*m3*m4**2*PvR*tman - Bsm2R*m3*m4**2*PvR*tman + BvgR*m1**2*mq1*PvR*tman + 
     &     3*BvgR*m2**2*mq1*PvR*tman + 2*Bvm1R*m1**2*m3*mq1*PvR*tman + 2*Bvm2R*m1**2*m3*mq1*PvR*tman + 3*Bvm1R*m2**2*m3*mq1*PvR*tman + 
     &     Bvm2R*m2**2*m3*mq1*PvR*tman + 4*BvgR*m3**2*mq1*PvR*tman + 4*Bvm1R*m3**3*mq1*PvR*tman + Bvm2R*m3**3*mq1*PvR*tman + BvgR*m4**2*mq1*PvR*tman + 
     &     2*Bvm1R*m3*m4**2*mq1*PvR*tman + Bvm2R*m3*m4**2*mq1*PvR*tman + BvgR*PsL*sman*tman + Bvm1R*m3*PsL*sman*tman - Bvm1L*mq1*PsL*sman*tman + 
     &     BsgR*PvR*sman*tman - Bvm1L*m1**2*PvR*sman*tman - Bvm2L*m1**2*PvR*sman*tman - 3*Bvm1L*m2**2*PvR*sman*tman - Bvm2L*m2**2*PvR*sman*tman + 
     &     Bsm1R*m3*PvR*sman*tman - 4*Bvm1L*m3**2*PvR*sman*tman - Bvm2L*m3**2*PvR*sman*tman - Bvm1L*m4**2*PvR*sman*tman - BvgR*mq1*PvR*sman*tman - 
     &     Bvm1R*m3*mq1*PvR*sman*tman + Bvm1L*PvR*sman**2*tman + BvgR*PsL*tman**2 + Bvm1R*m3*PsL*tman**2 + Bvm2R*m3*PsL*tman**2 - Bvm1L*mq1*PsL*tman**2 - 
     &     Bvm2L*mq1*PsL*tman**2 + BsgR*PvR*tman**2 + Bsm1R*m3*PvR*tman**2 + Bsm2R*m3*PvR*tman**2 - BvgR*mq1*PvR*tman**2 - Bvm1R*m3*mq1*PvR*tman**2 - 
     &     Bvm2R*m3*mq1*PvR*tman**2 + Bvm1L*PvR*sman*tman**2 + Bvm2L*PvR*sman*tman**2 - 3*BvgR*m1**2*PsL*uman - BvgR*m2**2*PsL*uman + BsgL*m3*PsL*uman - 
     &     4*Bvm1R*m1**2*m3*PsL*uman - 3*Bvm1R*m2**2*m3*PsL*uman - Bvm2R*m2**2*m3*PsL*uman - 4*BvgR*m3**2*PsL*uman - 5*Bvm1R*m3**3*PsL*uman - 
     &     BvgR*m4**2*PsL*uman - 3*Bvm1R*m3*m4**2*PsL*uman + 2*Bvm1L*m1**2*mq1*PsL*uman - BvgL*m3*mq1*PsL*uman + 2*Bvm1L*m3**2*mq1*PsL*uman - 
     &     3*BsgR*m1**2*PvR*uman - BsgR*m2**2*PvR*uman - 4*Bsm1R*m1**2*m3*PvR*uman - 3*Bsm1R*m2**2*m3*PvR*uman - Bsm2R*m2**2*m3*PvR*uman - 
     &     4*BsgR*m3**2*PvR*uman - 5*Bsm1R*m3**3*PvR*uman - BsgR*m4**2*PvR*uman - 3*Bsm1R*m3*m4**2*PvR*uman + 3*BvgR*m1**2*mq1*PvR*uman + 
     &     BvgR*m2**2*mq1*PvR*uman + 4*Bvm1R*m1**2*m3*mq1*PvR*uman + 3*Bvm1R*m2**2*m3*mq1*PvR*uman + Bvm2R*m2**2*m3*mq1*PvR*uman + 4*BvgR*m3**2*mq1*PvR*uman + 
     &     5*Bvm1R*m3**3*mq1*PvR*uman + BvgR*m4**2*mq1*PvR*uman + 3*Bvm1R*m3*m4**2*mq1*PvR*uman + Bvm1R*m3*PsL*sman*uman - 2*Bvm1L*m1**2*PvR*sman*uman + 
     &     Bsm1R*m3*PvR*sman*uman + BvgL*m3*PvR*sman*uman - 2*Bvm1L*m3**2*PvR*sman*uman - Bvm1R*m3*mq1*PvR*sman*uman + 2*BvgR*PsL*tman*uman + 
     &     3*Bvm1R*m3*PsL*tman*uman + Bvm2R*m3*PsL*tman*uman - 2*Bvm1L*mq1*PsL*tman*uman + 2*BsgR*PvR*tman*uman + 3*Bsm1R*m3*PvR*tman*uman + 
     &     Bsm2R*m3*PvR*tman*uman - 2*BvgR*mq1*PvR*tman*uman - 3*Bvm1R*m3*mq1*PvR*tman*uman - Bvm2R*m3*mq1*PvR*tman*uman + 2*Bvm1L*PvR*sman*tman*uman + 
     &     BvgR*PsL*uman**2 + 2*Bvm1R*m3*PsL*uman**2 + BsgR*PvR*uman**2 + 2*Bsm1R*m3*PvR*uman**2 - BvgR*mq1*PvR*uman**2 - 2*Bvm1R*m3*mq1*PvR*uman**2 - 
     &     BsgL*m1**2*m2**2*m3*PsL*xv + 2*BvgR*m1**2*m2**2*m3**2*PsL*xv + 2*Bvm1R*m1**2*m2**2*m3**3*PsL*xv + BvgR*m2**2*m3**4*PsL*xv + 
     &     Bvm1R*m2**2*m3**5*PsL*xv - BsgL*m2**2*m3*m4**2*PsL*xv - 2*Bvm2R*m1**2*m2**2*m3*m4**2*PsL*xv + BvgR*m2**2*m3**2*m4**2*PsL*xv + 
     &     Bvm1R*m2**2*m3**3*m4**2*PsL*xv - Bvm2R*m2**2*m3**3*m4**2*PsL*xv - Bvm2R*m2**2*m3*m4**4*PsL*xv + BvgL*m1**2*m2**2*m3*mq1*PsL*xv - 
     &     Bvm1L*m1**2*m2**2*m3**2*mq1*PsL*xv - Bvm1L*m2**2*m3**4*mq1*PsL*xv + Bvm2L*m1**2*m2**2*m4**2*mq1*PsL*xv + BvgL*m2**2*m3*m4**2*mq1*PsL*xv + 
     &     Bvm2L*m2**2*m3**2*m4**2*mq1*PsL*xv + 2*BsgR*m1**2*m2**2*m3**2*PvR*xv + 2*Bsm1R*m1**2*m2**2*m3**3*PvR*xv + BsgR*m2**2*m3**4*PvR*xv + 
     &     Bsm1R*m2**2*m3**5*PvR*xv - 2*Bsm2R*m1**2*m2**2*m3*m4**2*PvR*xv + BsgR*m2**2*m3**2*m4**2*PvR*xv + Bsm1R*m2**2*m3**3*m4**2*PvR*xv - 
     &     Bsm2R*m2**2*m3**3*m4**2*PvR*xv - Bsm2R*m2**2*m3*m4**4*PvR*xv - 2*BvgR*m1**2*m2**2*m3**2*mq1*PvR*xv - 2*Bvm1R*m1**2*m2**2*m3**3*mq1*PvR*xv - 
     &     BvgR*m2**2*m3**4*mq1*PvR*xv - Bvm1R*m2**2*m3**5*mq1*PvR*xv + 2*Bvm2R*m1**2*m2**2*m3*m4**2*mq1*PvR*xv - BvgR*m2**2*m3**2*m4**2*mq1*PvR*xv - 
     &     Bvm1R*m2**2*m3**3*m4**2*mq1*PvR*xv + Bvm2R*m2**2*m3**3*m4**2*mq1*PvR*xv + Bvm2R*m2**2*m3*m4**4*mq1*PvR*xv - BvgR*m1**2*m2**2*PsL*sman*xv - 
     &     2*Bvm1R*m1**2*m2**2*m3*PsL*sman*xv - BvgR*m2**2*m3**2*PsL*sman*xv - Bvm1R*m2**2*m3**3*PsL*sman*xv - Bvm1R*m2**2*m3*m4**2*PsL*sman*xv + 
     &     Bvm1L*m1**2*m2**2*mq1*PsL*sman*xv + Bvm1L*m2**2*m3**2*mq1*PsL*sman*xv - BsgR*m1**2*m2**2*PvR*sman*xv - 2*Bsm1R*m1**2*m2**2*m3*PvR*sman*xv - 
     &     BvgL*m1**2*m2**2*m3*PvR*sman*xv - BsgR*m2**2*m3**2*PvR*sman*xv + Bvm1L*m1**2*m2**2*m3**2*PvR*sman*xv - Bsm1R*m2**2*m3**3*PvR*sman*xv + 
     &     Bvm1L*m2**2*m3**4*PvR*sman*xv - Bvm2L*m1**2*m2**2*m4**2*PvR*sman*xv - Bsm1R*m2**2*m3*m4**2*PvR*sman*xv - BvgL*m2**2*m3*m4**2*PvR*sman*xv - 
     &     Bvm2L*m2**2*m3**2*m4**2*PvR*sman*xv + BvgR*m1**2*m2**2*mq1*PvR*sman*xv + 2*Bvm1R*m1**2*m2**2*m3*mq1*PvR*sman*xv + 
     &     BvgR*m2**2*m3**2*mq1*PvR*sman*xv + Bvm1R*m2**2*m3**3*mq1*PvR*sman*xv + Bvm1R*m2**2*m3*m4**2*mq1*PvR*sman*xv - Bvm1L*m1**2*m2**2*PvR*sman**2*xv - 
     &     Bvm1L*m2**2*m3**2*PvR*sman**2*xv + BsgL*m1**2*m3*PsL*tman*xv - 2*BvgR*m1**2*m3**2*PsL*tman*xv - BvgR*m2**2*m3**2*PsL*tman*xv - 
     &     2*Bvm1R*m1**2*m3**3*PsL*tman*xv - Bvm1R*m2**2*m3**3*PsL*tman*xv - BvgR*m3**4*PsL*tman*xv - Bvm1R*m3**5*PsL*tman*xv + BsgL*m3*m4**2*PsL*tman*xv + 
     &     2*Bvm2R*m1**2*m3*m4**2*PsL*tman*xv + Bvm2R*m2**2*m3*m4**2*PsL*tman*xv - BvgR*m3**2*m4**2*PsL*tman*xv - Bvm1R*m3**3*m4**2*PsL*tman*xv + 
     &     Bvm2R*m3**3*m4**2*PsL*tman*xv + Bvm2R*m3*m4**4*PsL*tman*xv - BvgL*m1**2*m3*mq1*PsL*tman*xv + Bvm1L*m1**2*m3**2*mq1*PsL*tman*xv + 
     &     Bvm1L*m2**2*m3**2*mq1*PsL*tman*xv + Bvm1L*m3**4*mq1*PsL*tman*xv - Bvm2L*m1**2*m4**2*mq1*PsL*tman*xv - Bvm2L*m2**2*m4**2*mq1*PsL*tman*xv - 
     &     BvgL*m3*m4**2*mq1*PsL*tman*xv - Bvm2L*m3**2*m4**2*mq1*PsL*tman*xv - 2*BsgR*m1**2*m3**2*PvR*tman*xv - BsgR*m2**2*m3**2*PvR*tman*xv - 
     &     2*Bsm1R*m1**2*m3**3*PvR*tman*xv - Bsm1R*m2**2*m3**3*PvR*tman*xv - BsgR*m3**4*PvR*tman*xv - Bsm1R*m3**5*PvR*tman*xv + 
     &     2*Bsm2R*m1**2*m3*m4**2*PvR*tman*xv + Bsm2R*m2**2*m3*m4**2*PvR*tman*xv - BsgR*m3**2*m4**2*PvR*tman*xv - Bsm1R*m3**3*m4**2*PvR*tman*xv + 
     &     Bsm2R*m3**3*m4**2*PvR*tman*xv + Bsm2R*m3*m4**4*PvR*tman*xv + 2*BvgR*m1**2*m3**2*mq1*PvR*tman*xv + BvgR*m2**2*m3**2*mq1*PvR*tman*xv + 
     &     2*Bvm1R*m1**2*m3**3*mq1*PvR*tman*xv + Bvm1R*m2**2*m3**3*mq1*PvR*tman*xv + BvgR*m3**4*mq1*PvR*tman*xv + Bvm1R*m3**5*mq1*PvR*tman*xv - 
     &     2*Bvm2R*m1**2*m3*m4**2*mq1*PvR*tman*xv - Bvm2R*m2**2*m3*m4**2*mq1*PvR*tman*xv + BvgR*m3**2*m4**2*mq1*PvR*tman*xv + 
     &     Bvm1R*m3**3*m4**2*mq1*PvR*tman*xv - Bvm2R*m3**3*m4**2*mq1*PvR*tman*xv - Bvm2R*m3*m4**4*mq1*PvR*tman*xv + BvgR*m1**2*PsL*sman*tman*xv + 
     &     BvgR*m2**2*PsL*sman*tman*xv + 2*Bvm1R*m1**2*m3*PsL*sman*tman*xv + Bvm1R*m2**2*m3*PsL*sman*tman*xv + BvgR*m3**2*PsL*sman*tman*xv + 
     &     Bvm1R*m3**3*PsL*sman*tman*xv + Bvm1R*m3*m4**2*PsL*sman*tman*xv - Bvm1L*m1**2*mq1*PsL*sman*tman*xv - Bvm1L*m2**2*mq1*PsL*sman*tman*xv - 
     &     Bvm1L*m3**2*mq1*PsL*sman*tman*xv + BsgR*m1**2*PvR*sman*tman*xv + BsgR*m2**2*PvR*sman*tman*xv + 2*Bsm1R*m1**2*m3*PvR*sman*tman*xv + 
     &     BvgL*m1**2*m3*PvR*sman*tman*xv + Bsm1R*m2**2*m3*PvR*sman*tman*xv + BsgR*m3**2*PvR*sman*tman*xv - Bvm1L*m1**2*m3**2*PvR*sman*tman*xv - 
     &     Bvm1L*m2**2*m3**2*PvR*sman*tman*xv + Bsm1R*m3**3*PvR*sman*tman*xv - Bvm1L*m3**4*PvR*sman*tman*xv + Bvm2L*m1**2*m4**2*PvR*sman*tman*xv + 
     &     Bvm2L*m2**2*m4**2*PvR*sman*tman*xv + Bsm1R*m3*m4**2*PvR*sman*tman*xv + BvgL*m3*m4**2*PvR*sman*tman*xv + Bvm2L*m3**2*m4**2*PvR*sman*tman*xv - 
     &     BvgR*m1**2*mq1*PvR*sman*tman*xv - BvgR*m2**2*mq1*PvR*sman*tman*xv - 2*Bvm1R*m1**2*m3*mq1*PvR*sman*tman*xv - Bvm1R*m2**2*m3*mq1*PvR*sman*tman*xv - 
     &     BvgR*m3**2*mq1*PvR*sman*tman*xv - Bvm1R*m3**3*mq1*PvR*sman*tman*xv - Bvm1R*m3*m4**2*mq1*PvR*sman*tman*xv + Bvm1L*m1**2*PvR*sman**2*tman*xv + 
     &     Bvm1L*m2**2*PvR*sman**2*tman*xv + Bvm1L*m3**2*PvR*sman**2*tman*xv + BvgR*m3**2*PsL*tman**2*xv + Bvm1R*m3**3*PsL*tman**2*xv - 
     &     Bvm2R*m3*m4**2*PsL*tman**2*xv - Bvm1L*m3**2*mq1*PsL*tman**2*xv + Bvm2L*m4**2*mq1*PsL*tman**2*xv + BsgR*m3**2*PvR*tman**2*xv + 
     &     Bsm1R*m3**3*PvR*tman**2*xv - Bsm2R*m3*m4**2*PvR*tman**2*xv - BvgR*m3**2*mq1*PvR*tman**2*xv - Bvm1R*m3**3*mq1*PvR*tman**2*xv + 
     &     Bvm2R*m3*m4**2*mq1*PvR*tman**2*xv - BvgR*PsL*sman*tman**2*xv - Bvm1R*m3*PsL*sman*tman**2*xv + Bvm1L*mq1*PsL*sman*tman**2*xv - 
     &     BsgR*PvR*sman*tman**2*xv - Bsm1R*m3*PvR*sman*tman**2*xv + Bvm1L*m3**2*PvR*sman*tman**2*xv - Bvm2L*m4**2*PvR*sman*tman**2*xv + 
     &     BvgR*mq1*PvR*sman*tman**2*xv + Bvm1R*m3*mq1*PvR*sman*tman**2*xv - Bvm1L*PvR*sman**2*tman**2*xv + BsgL*m2**2*m3*PsL*uman*xv - 
     &     BvgR*m2**2*m3**2*PsL*uman*xv - Bvm1R*m2**2*m3**3*PsL*uman*xv + Bvm2R*m2**2*m3*m4**2*PsL*uman*xv - BvgL*m2**2*m3*mq1*PsL*uman*xv - 
     &     BsgR*m2**2*m3**2*PvR*uman*xv - Bsm1R*m2**2*m3**3*PvR*uman*xv + Bsm2R*m2**2*m3*m4**2*PvR*uman*xv + BvgR*m2**2*m3**2*mq1*PvR*uman*xv + 
     &     Bvm1R*m2**2*m3**3*mq1*PvR*uman*xv - Bvm2R*m2**2*m3*m4**2*mq1*PvR*uman*xv + Bvm1R*m2**2*m3*PsL*sman*uman*xv + Bsm1R*m2**2*m3*PvR*sman*uman*xv + 
     &     BvgL*m2**2*m3*PvR*sman*uman*xv - Bvm1R*m2**2*m3*mq1*PvR*sman*uman*xv - BsgL*m3*PsL*tman*uman*xv + BvgR*m3**2*PsL*tman*uman*xv + 
     &     Bvm1R*m3**3*PsL*tman*uman*xv - Bvm2R*m3*m4**2*PsL*tman*uman*xv + BvgL*m3*mq1*PsL*tman*uman*xv + BsgR*m3**2*PvR*tman*uman*xv + 
     &     Bsm1R*m3**3*PvR*tman*uman*xv - Bsm2R*m3*m4**2*PvR*tman*uman*xv - BvgR*m3**2*mq1*PvR*tman*uman*xv - Bvm1R*m3**3*mq1*PvR*tman*uman*xv + 
     &     Bvm2R*m3*m4**2*mq1*PvR*tman*uman*xv - Bvm1R*m3*PsL*sman*tman*uman*xv - Bsm1R*m3*PvR*sman*tman*uman*xv - BvgL*m3*PvR*sman*tman*uman*xv + 
     &     Bvm1R*m3*mq1*PvR*sman*tman*uman*xv - Bsm2L*PsL*(m2**2 - tman)*(m1**2 + m3**2 - tman)*(-1 + m4**2*xv) - 
     &     Bsm1L*PsL*(m1**2 + m3**2 - tman)*(-m4**2 + sman + tman + 2*uman - sman*tman*xv + m2**2*(-3 - m3**2*xv + sman*xv) + m3**2*(-3 + tman*xv))) - 
     &  2*A1sR*A2Rc*CLc*(4*BvgL*m1**2*m2**2*PsR - 3*BsgR*m1**2*m3*PsR - 2*BsgR*m2**2*m3*PsR + 6*Bvm1L*m1**2*m2**2*m3*PsR + 2*Bvm2L*m1**2*m2**2*m3*PsR + 
     &     6*BvgL*m1**2*m3**2*PsR + 4*BvgL*m2**2*m3**2*PsR + 6*Bvm1L*m1**2*m3**3*PsR + 3*Bvm1L*m2**2*m3**3*PsR + Bvm2L*m2**2*m3**3*PsR + 3*BvgL*m3**4*PsR + 
     &     3*Bvm1L*m3**5*PsR - BsgR*m3*m4**2*PsR + 2*Bvm1L*m1**2*m3*m4**2*PsR + 3*Bvm1L*m2**2*m3*m4**2*PsR + Bvm2L*m2**2*m3*m4**2*PsR + 
     &     3*BvgL*m3**2*m4**2*PsR + 4*Bvm1L*m3**3*m4**2*PsR + Bvm1L*m3*m4**4*PsR - 3*Bvm1R*m1**2*m2**2*mq1*PsR - Bvm2R*m1**2*m2**2*mq1*PsR + 
     &     3*BvgR*m1**2*m3*mq1*PsR + 2*BvgR*m2**2*m3*mq1*PsR - 3*Bvm1R*m1**2*m3**2*mq1*PsR - 3*Bvm1R*m2**2*m3**2*mq1*PsR - Bvm2R*m2**2*m3**2*mq1*PsR - 
     &     3*Bvm1R*m3**4*mq1*PsR - Bvm1R*m1**2*m4**2*mq1*PsR + BvgR*m3*m4**2*mq1*PsR - Bvm1R*m3**2*m4**2*mq1*PsR + 4*BsgL*m1**2*m2**2*PvL + 
     &     6*Bsm1L*m1**2*m2**2*m3*PvL + 2*Bsm2L*m1**2*m2**2*m3*PvL + 6*BsgL*m1**2*m3**2*PvL + 4*BsgL*m2**2*m3**2*PvL + 6*Bsm1L*m1**2*m3**3*PvL + 
     &     3*Bsm1L*m2**2*m3**3*PvL + Bsm2L*m2**2*m3**3*PvL + 3*BsgL*m3**4*PvL + 3*Bsm1L*m3**5*PvL + 2*Bsm1L*m1**2*m3*m4**2*PvL + 3*Bsm1L*m2**2*m3*m4**2*PvL + 
     &     Bsm2L*m2**2*m3*m4**2*PvL + 3*BsgL*m3**2*m4**2*PvL + 4*Bsm1L*m3**3*m4**2*PvL + Bsm1L*m3*m4**4*PvL - 4*BvgL*m1**2*m2**2*mq1*PvL - 
     &     6*Bvm1L*m1**2*m2**2*m3*mq1*PvL - 2*Bvm2L*m1**2*m2**2*m3*mq1*PvL - 6*BvgL*m1**2*m3**2*mq1*PvL - 4*BvgL*m2**2*m3**2*mq1*PvL - 
     &     6*Bvm1L*m1**2*m3**3*mq1*PvL - 3*Bvm1L*m2**2*m3**3*mq1*PvL - Bvm2L*m2**2*m3**3*mq1*PvL - 3*BvgL*m3**4*mq1*PvL - 3*Bvm1L*m3**5*mq1*PvL - 
     &     2*Bvm1L*m1**2*m3*m4**2*mq1*PvL - 3*Bvm1L*m2**2*m3*m4**2*mq1*PvL - Bvm2L*m2**2*m3*m4**2*mq1*PvL - 3*BvgL*m3**2*m4**2*mq1*PvL - 
     &     4*Bvm1L*m3**3*m4**2*mq1*PvL - Bvm1L*m3*m4**4*mq1*PvL + BvgL*m2**2*PsR*sman + 2*BsgR*m3*PsR*sman - 2*Bvm1L*m1**2*m3*PsR*sman - 
     &     2*BvgL*m3**2*PsR*sman - Bvm1L*m3**3*PsR*sman + BvgL*m4**2*PsR*sman - Bvm1L*m3*m4**2*PsR*sman + Bvm1R*m1**2*mq1*PsR*sman - 2*BvgR*m3*mq1*PsR*sman + 
     &     Bvm1R*m3**2*mq1*PsR*sman + BsgL*m2**2*PvL*sman + 3*Bvm1R*m1**2*m2**2*PvL*sman + Bvm2R*m1**2*m2**2*PvL*sman - 2*Bsm1L*m1**2*m3*PvL*sman - 
     &     3*BvgR*m1**2*m3*PvL*sman - 2*BvgR*m2**2*m3*PvL*sman - 2*BsgL*m3**2*PvL*sman + 3*Bvm1R*m1**2*m3**2*PvL*sman + 3*Bvm1R*m2**2*m3**2*PvL*sman + 
     &     Bvm2R*m2**2*m3**2*PvL*sman - Bsm1L*m3**3*PvL*sman + 3*Bvm1R*m3**4*PvL*sman + BsgL*m4**2*PvL*sman + Bvm1R*m1**2*m4**2*PvL*sman - 
     &     Bsm1L*m3*m4**2*PvL*sman - BvgR*m3*m4**2*PvL*sman + Bvm1R*m3**2*m4**2*PvL*sman - BvgL*m2**2*mq1*PvL*sman + 2*Bvm1L*m1**2*m3*mq1*PvL*sman + 
     &     2*BvgL*m3**2*mq1*PvL*sman + Bvm1L*m3**3*mq1*PvL*sman - BvgL*m4**2*mq1*PvL*sman + Bvm1L*m3*m4**2*mq1*PvL*sman - BvgL*PsR*sman**2 - 
     &     BsgL*PvL*sman**2 - Bvm1R*m1**2*PvL*sman**2 + 2*BvgR*m3*PvL*sman**2 - Bvm1R*m3**2*PvL*sman**2 + BvgL*mq1*PvL*sman**2 - BvgL*m1**2*PsR*tman - 
     &     3*BvgL*m2**2*PsR*tman - 2*Bvm1L*m1**2*m3*PsR*tman - 2*Bvm2L*m1**2*m3*PsR*tman - 3*Bvm1L*m2**2*m3*PsR*tman - Bvm2L*m2**2*m3*PsR*tman - 
     &     4*BvgL*m3**2*PsR*tman - 4*Bvm1L*m3**3*PsR*tman - Bvm2L*m3**3*PsR*tman - BvgL*m4**2*PsR*tman - 2*Bvm1L*m3*m4**2*PsR*tman - Bvm2L*m3*m4**2*PsR*tman + 
     &     Bvm1R*m1**2*mq1*PsR*tman + Bvm2R*m1**2*mq1*PsR*tman + 3*Bvm1R*m2**2*mq1*PsR*tman + Bvm2R*m2**2*mq1*PsR*tman + 4*Bvm1R*m3**2*mq1*PsR*tman + 
     &     Bvm2R*m3**2*mq1*PsR*tman + Bvm1R*m4**2*mq1*PsR*tman - BsgL*m1**2*PvL*tman - 3*BsgL*m2**2*PvL*tman - 2*Bsm1L*m1**2*m3*PvL*tman - 
     &     2*Bsm2L*m1**2*m3*PvL*tman - 3*Bsm1L*m2**2*m3*PvL*tman - Bsm2L*m2**2*m3*PvL*tman - 4*BsgL*m3**2*PvL*tman - 4*Bsm1L*m3**3*PvL*tman - 
     &     Bsm2L*m3**3*PvL*tman - BsgL*m4**2*PvL*tman - 2*Bsm1L*m3*m4**2*PvL*tman - Bsm2L*m3*m4**2*PvL*tman + BvgL*m1**2*mq1*PvL*tman + 
     &     3*BvgL*m2**2*mq1*PvL*tman + 2*Bvm1L*m1**2*m3*mq1*PvL*tman + 2*Bvm2L*m1**2*m3*mq1*PvL*tman + 3*Bvm1L*m2**2*m3*mq1*PvL*tman + 
     &     Bvm2L*m2**2*m3*mq1*PvL*tman + 4*BvgL*m3**2*mq1*PvL*tman + 4*Bvm1L*m3**3*mq1*PvL*tman + Bvm2L*m3**3*mq1*PvL*tman + BvgL*m4**2*mq1*PvL*tman + 
     &     2*Bvm1L*m3*m4**2*mq1*PvL*tman + Bvm2L*m3*m4**2*mq1*PvL*tman + BvgL*PsR*sman*tman + Bvm1L*m3*PsR*sman*tman - Bvm1R*mq1*PsR*sman*tman + 
     &     BsgL*PvL*sman*tman - Bvm1R*m1**2*PvL*sman*tman - Bvm2R*m1**2*PvL*sman*tman - 3*Bvm1R*m2**2*PvL*sman*tman - Bvm2R*m2**2*PvL*sman*tman + 
     &     Bsm1L*m3*PvL*sman*tman - 4*Bvm1R*m3**2*PvL*sman*tman - Bvm2R*m3**2*PvL*sman*tman - Bvm1R*m4**2*PvL*sman*tman - BvgL*mq1*PvL*sman*tman - 
     &     Bvm1L*m3*mq1*PvL*sman*tman + Bvm1R*PvL*sman**2*tman + BvgL*PsR*tman**2 + Bvm1L*m3*PsR*tman**2 + Bvm2L*m3*PsR*tman**2 - Bvm1R*mq1*PsR*tman**2 - 
     &     Bvm2R*mq1*PsR*tman**2 + BsgL*PvL*tman**2 + Bsm1L*m3*PvL*tman**2 + Bsm2L*m3*PvL*tman**2 - BvgL*mq1*PvL*tman**2 - Bvm1L*m3*mq1*PvL*tman**2 - 
     &     Bvm2L*m3*mq1*PvL*tman**2 + Bvm1R*PvL*sman*tman**2 + Bvm2R*PvL*sman*tman**2 - 3*BvgL*m1**2*PsR*uman - BvgL*m2**2*PsR*uman + BsgR*m3*PsR*uman - 
     &     4*Bvm1L*m1**2*m3*PsR*uman - 3*Bvm1L*m2**2*m3*PsR*uman - Bvm2L*m2**2*m3*PsR*uman - 4*BvgL*m3**2*PsR*uman - 5*Bvm1L*m3**3*PsR*uman - 
     &     BvgL*m4**2*PsR*uman - 3*Bvm1L*m3*m4**2*PsR*uman + 2*Bvm1R*m1**2*mq1*PsR*uman - BvgR*m3*mq1*PsR*uman + 2*Bvm1R*m3**2*mq1*PsR*uman - 
     &     3*BsgL*m1**2*PvL*uman - BsgL*m2**2*PvL*uman - 4*Bsm1L*m1**2*m3*PvL*uman - 3*Bsm1L*m2**2*m3*PvL*uman - Bsm2L*m2**2*m3*PvL*uman - 
     &     4*BsgL*m3**2*PvL*uman - 5*Bsm1L*m3**3*PvL*uman - BsgL*m4**2*PvL*uman - 3*Bsm1L*m3*m4**2*PvL*uman + 3*BvgL*m1**2*mq1*PvL*uman + 
     &     BvgL*m2**2*mq1*PvL*uman + 4*Bvm1L*m1**2*m3*mq1*PvL*uman + 3*Bvm1L*m2**2*m3*mq1*PvL*uman + Bvm2L*m2**2*m3*mq1*PvL*uman + 4*BvgL*m3**2*mq1*PvL*uman + 
     &     5*Bvm1L*m3**3*mq1*PvL*uman + BvgL*m4**2*mq1*PvL*uman + 3*Bvm1L*m3*m4**2*mq1*PvL*uman + Bvm1L*m3*PsR*sman*uman - 2*Bvm1R*m1**2*PvL*sman*uman + 
     &     Bsm1L*m3*PvL*sman*uman + BvgR*m3*PvL*sman*uman - 2*Bvm1R*m3**2*PvL*sman*uman - Bvm1L*m3*mq1*PvL*sman*uman + 2*BvgL*PsR*tman*uman + 
     &     3*Bvm1L*m3*PsR*tman*uman + Bvm2L*m3*PsR*tman*uman - 2*Bvm1R*mq1*PsR*tman*uman + 2*BsgL*PvL*tman*uman + 3*Bsm1L*m3*PvL*tman*uman + 
     &     Bsm2L*m3*PvL*tman*uman - 2*BvgL*mq1*PvL*tman*uman - 3*Bvm1L*m3*mq1*PvL*tman*uman - Bvm2L*m3*mq1*PvL*tman*uman + 2*Bvm1R*PvL*sman*tman*uman + 
     &     BvgL*PsR*uman**2 + 2*Bvm1L*m3*PsR*uman**2 + BsgL*PvL*uman**2 + 2*Bsm1L*m3*PvL*uman**2 - BvgL*mq1*PvL*uman**2 - 2*Bvm1L*m3*mq1*PvL*uman**2 - 
     &     BsgR*m1**2*m2**2*m3*PsR*xv + 2*BvgL*m1**2*m2**2*m3**2*PsR*xv + 2*Bvm1L*m1**2*m2**2*m3**3*PsR*xv + BvgL*m2**2*m3**4*PsR*xv + 
     &     Bvm1L*m2**2*m3**5*PsR*xv - BsgR*m2**2*m3*m4**2*PsR*xv - 2*Bvm2L*m1**2*m2**2*m3*m4**2*PsR*xv + BvgL*m2**2*m3**2*m4**2*PsR*xv + 
     &     Bvm1L*m2**2*m3**3*m4**2*PsR*xv - Bvm2L*m2**2*m3**3*m4**2*PsR*xv - Bvm2L*m2**2*m3*m4**4*PsR*xv + BvgR*m1**2*m2**2*m3*mq1*PsR*xv - 
     &     Bvm1R*m1**2*m2**2*m3**2*mq1*PsR*xv - Bvm1R*m2**2*m3**4*mq1*PsR*xv + Bvm2R*m1**2*m2**2*m4**2*mq1*PsR*xv + BvgR*m2**2*m3*m4**2*mq1*PsR*xv + 
     &     Bvm2R*m2**2*m3**2*m4**2*mq1*PsR*xv + 2*BsgL*m1**2*m2**2*m3**2*PvL*xv + 2*Bsm1L*m1**2*m2**2*m3**3*PvL*xv + BsgL*m2**2*m3**4*PvL*xv + 
     &     Bsm1L*m2**2*m3**5*PvL*xv - 2*Bsm2L*m1**2*m2**2*m3*m4**2*PvL*xv + BsgL*m2**2*m3**2*m4**2*PvL*xv + Bsm1L*m2**2*m3**3*m4**2*PvL*xv - 
     &     Bsm2L*m2**2*m3**3*m4**2*PvL*xv - Bsm2L*m2**2*m3*m4**4*PvL*xv - 2*BvgL*m1**2*m2**2*m3**2*mq1*PvL*xv - 2*Bvm1L*m1**2*m2**2*m3**3*mq1*PvL*xv - 
     &     BvgL*m2**2*m3**4*mq1*PvL*xv - Bvm1L*m2**2*m3**5*mq1*PvL*xv + 2*Bvm2L*m1**2*m2**2*m3*m4**2*mq1*PvL*xv - BvgL*m2**2*m3**2*m4**2*mq1*PvL*xv - 
     &     Bvm1L*m2**2*m3**3*m4**2*mq1*PvL*xv + Bvm2L*m2**2*m3**3*m4**2*mq1*PvL*xv + Bvm2L*m2**2*m3*m4**4*mq1*PvL*xv - BvgL*m1**2*m2**2*PsR*sman*xv - 
     &     2*Bvm1L*m1**2*m2**2*m3*PsR*sman*xv - BvgL*m2**2*m3**2*PsR*sman*xv - Bvm1L*m2**2*m3**3*PsR*sman*xv - Bvm1L*m2**2*m3*m4**2*PsR*sman*xv + 
     &     Bvm1R*m1**2*m2**2*mq1*PsR*sman*xv + Bvm1R*m2**2*m3**2*mq1*PsR*sman*xv - BsgL*m1**2*m2**2*PvL*sman*xv - 2*Bsm1L*m1**2*m2**2*m3*PvL*sman*xv - 
     &     BvgR*m1**2*m2**2*m3*PvL*sman*xv - BsgL*m2**2*m3**2*PvL*sman*xv + Bvm1R*m1**2*m2**2*m3**2*PvL*sman*xv - Bsm1L*m2**2*m3**3*PvL*sman*xv + 
     &     Bvm1R*m2**2*m3**4*PvL*sman*xv - Bvm2R*m1**2*m2**2*m4**2*PvL*sman*xv - Bsm1L*m2**2*m3*m4**2*PvL*sman*xv - BvgR*m2**2*m3*m4**2*PvL*sman*xv - 
     &     Bvm2R*m2**2*m3**2*m4**2*PvL*sman*xv + BvgL*m1**2*m2**2*mq1*PvL*sman*xv + 2*Bvm1L*m1**2*m2**2*m3*mq1*PvL*sman*xv + 
     &     BvgL*m2**2*m3**2*mq1*PvL*sman*xv + Bvm1L*m2**2*m3**3*mq1*PvL*sman*xv + Bvm1L*m2**2*m3*m4**2*mq1*PvL*sman*xv - Bvm1R*m1**2*m2**2*PvL*sman**2*xv - 
     &     Bvm1R*m2**2*m3**2*PvL*sman**2*xv + BsgR*m1**2*m3*PsR*tman*xv - 2*BvgL*m1**2*m3**2*PsR*tman*xv - BvgL*m2**2*m3**2*PsR*tman*xv - 
     &     2*Bvm1L*m1**2*m3**3*PsR*tman*xv - Bvm1L*m2**2*m3**3*PsR*tman*xv - BvgL*m3**4*PsR*tman*xv - Bvm1L*m3**5*PsR*tman*xv + BsgR*m3*m4**2*PsR*tman*xv + 
     &     2*Bvm2L*m1**2*m3*m4**2*PsR*tman*xv + Bvm2L*m2**2*m3*m4**2*PsR*tman*xv - BvgL*m3**2*m4**2*PsR*tman*xv - Bvm1L*m3**3*m4**2*PsR*tman*xv + 
     &     Bvm2L*m3**3*m4**2*PsR*tman*xv + Bvm2L*m3*m4**4*PsR*tman*xv - BvgR*m1**2*m3*mq1*PsR*tman*xv + Bvm1R*m1**2*m3**2*mq1*PsR*tman*xv + 
     &     Bvm1R*m2**2*m3**2*mq1*PsR*tman*xv + Bvm1R*m3**4*mq1*PsR*tman*xv - Bvm2R*m1**2*m4**2*mq1*PsR*tman*xv - Bvm2R*m2**2*m4**2*mq1*PsR*tman*xv - 
     &     BvgR*m3*m4**2*mq1*PsR*tman*xv - Bvm2R*m3**2*m4**2*mq1*PsR*tman*xv - 2*BsgL*m1**2*m3**2*PvL*tman*xv - BsgL*m2**2*m3**2*PvL*tman*xv - 
     &     2*Bsm1L*m1**2*m3**3*PvL*tman*xv - Bsm1L*m2**2*m3**3*PvL*tman*xv - BsgL*m3**4*PvL*tman*xv - Bsm1L*m3**5*PvL*tman*xv + 
     &     2*Bsm2L*m1**2*m3*m4**2*PvL*tman*xv + Bsm2L*m2**2*m3*m4**2*PvL*tman*xv - BsgL*m3**2*m4**2*PvL*tman*xv - Bsm1L*m3**3*m4**2*PvL*tman*xv + 
     &     Bsm2L*m3**3*m4**2*PvL*tman*xv + Bsm2L*m3*m4**4*PvL*tman*xv + 2*BvgL*m1**2*m3**2*mq1*PvL*tman*xv + BvgL*m2**2*m3**2*mq1*PvL*tman*xv + 
     &     2*Bvm1L*m1**2*m3**3*mq1*PvL*tman*xv + Bvm1L*m2**2*m3**3*mq1*PvL*tman*xv + BvgL*m3**4*mq1*PvL*tman*xv + Bvm1L*m3**5*mq1*PvL*tman*xv - 
     &     2*Bvm2L*m1**2*m3*m4**2*mq1*PvL*tman*xv - Bvm2L*m2**2*m3*m4**2*mq1*PvL*tman*xv + BvgL*m3**2*m4**2*mq1*PvL*tman*xv + 
     &     Bvm1L*m3**3*m4**2*mq1*PvL*tman*xv - Bvm2L*m3**3*m4**2*mq1*PvL*tman*xv - Bvm2L*m3*m4**4*mq1*PvL*tman*xv + BvgL*m1**2*PsR*sman*tman*xv + 
     &     BvgL*m2**2*PsR*sman*tman*xv + 2*Bvm1L*m1**2*m3*PsR*sman*tman*xv + Bvm1L*m2**2*m3*PsR*sman*tman*xv + BvgL*m3**2*PsR*sman*tman*xv + 
     &     Bvm1L*m3**3*PsR*sman*tman*xv + Bvm1L*m3*m4**2*PsR*sman*tman*xv - Bvm1R*m1**2*mq1*PsR*sman*tman*xv - Bvm1R*m2**2*mq1*PsR*sman*tman*xv - 
     &     Bvm1R*m3**2*mq1*PsR*sman*tman*xv + BsgL*m1**2*PvL*sman*tman*xv + BsgL*m2**2*PvL*sman*tman*xv + 2*Bsm1L*m1**2*m3*PvL*sman*tman*xv + 
     &     BvgR*m1**2*m3*PvL*sman*tman*xv + Bsm1L*m2**2*m3*PvL*sman*tman*xv + BsgL*m3**2*PvL*sman*tman*xv - Bvm1R*m1**2*m3**2*PvL*sman*tman*xv - 
     &     Bvm1R*m2**2*m3**2*PvL*sman*tman*xv + Bsm1L*m3**3*PvL*sman*tman*xv - Bvm1R*m3**4*PvL*sman*tman*xv + Bvm2R*m1**2*m4**2*PvL*sman*tman*xv + 
     &     Bvm2R*m2**2*m4**2*PvL*sman*tman*xv + Bsm1L*m3*m4**2*PvL*sman*tman*xv + BvgR*m3*m4**2*PvL*sman*tman*xv + Bvm2R*m3**2*m4**2*PvL*sman*tman*xv - 
     &     BvgL*m1**2*mq1*PvL*sman*tman*xv - BvgL*m2**2*mq1*PvL*sman*tman*xv - 2*Bvm1L*m1**2*m3*mq1*PvL*sman*tman*xv - Bvm1L*m2**2*m3*mq1*PvL*sman*tman*xv - 
     &     BvgL*m3**2*mq1*PvL*sman*tman*xv - Bvm1L*m3**3*mq1*PvL*sman*tman*xv - Bvm1L*m3*m4**2*mq1*PvL*sman*tman*xv + Bvm1R*m1**2*PvL*sman**2*tman*xv + 
     &     Bvm1R*m2**2*PvL*sman**2*tman*xv + Bvm1R*m3**2*PvL*sman**2*tman*xv + BvgL*m3**2*PsR*tman**2*xv + Bvm1L*m3**3*PsR*tman**2*xv - 
     &     Bvm2L*m3*m4**2*PsR*tman**2*xv - Bvm1R*m3**2*mq1*PsR*tman**2*xv + Bvm2R*m4**2*mq1*PsR*tman**2*xv + BsgL*m3**2*PvL*tman**2*xv + 
     &     Bsm1L*m3**3*PvL*tman**2*xv - Bsm2L*m3*m4**2*PvL*tman**2*xv - BvgL*m3**2*mq1*PvL*tman**2*xv - Bvm1L*m3**3*mq1*PvL*tman**2*xv + 
     &     Bvm2L*m3*m4**2*mq1*PvL*tman**2*xv - BvgL*PsR*sman*tman**2*xv - Bvm1L*m3*PsR*sman*tman**2*xv + Bvm1R*mq1*PsR*sman*tman**2*xv - 
     &     BsgL*PvL*sman*tman**2*xv - Bsm1L*m3*PvL*sman*tman**2*xv + Bvm1R*m3**2*PvL*sman*tman**2*xv - Bvm2R*m4**2*PvL*sman*tman**2*xv + 
     &     BvgL*mq1*PvL*sman*tman**2*xv + Bvm1L*m3*mq1*PvL*sman*tman**2*xv - Bvm1R*PvL*sman**2*tman**2*xv + BsgR*m2**2*m3*PsR*uman*xv - 
     &     BvgL*m2**2*m3**2*PsR*uman*xv - Bvm1L*m2**2*m3**3*PsR*uman*xv + Bvm2L*m2**2*m3*m4**2*PsR*uman*xv - BvgR*m2**2*m3*mq1*PsR*uman*xv - 
     &     BsgL*m2**2*m3**2*PvL*uman*xv - Bsm1L*m2**2*m3**3*PvL*uman*xv + Bsm2L*m2**2*m3*m4**2*PvL*uman*xv + BvgL*m2**2*m3**2*mq1*PvL*uman*xv + 
     &     Bvm1L*m2**2*m3**3*mq1*PvL*uman*xv - Bvm2L*m2**2*m3*m4**2*mq1*PvL*uman*xv + Bvm1L*m2**2*m3*PsR*sman*uman*xv + Bsm1L*m2**2*m3*PvL*sman*uman*xv + 
     &     BvgR*m2**2*m3*PvL*sman*uman*xv - Bvm1L*m2**2*m3*mq1*PvL*sman*uman*xv - BsgR*m3*PsR*tman*uman*xv + BvgL*m3**2*PsR*tman*uman*xv + 
     &     Bvm1L*m3**3*PsR*tman*uman*xv - Bvm2L*m3*m4**2*PsR*tman*uman*xv + BvgR*m3*mq1*PsR*tman*uman*xv + BsgL*m3**2*PvL*tman*uman*xv + 
     &     Bsm1L*m3**3*PvL*tman*uman*xv - Bsm2L*m3*m4**2*PvL*tman*uman*xv - BvgL*m3**2*mq1*PvL*tman*uman*xv - Bvm1L*m3**3*mq1*PvL*tman*uman*xv + 
     &     Bvm2L*m3*m4**2*mq1*PvL*tman*uman*xv - Bvm1L*m3*PsR*sman*tman*uman*xv - Bsm1L*m3*PvL*sman*tman*uman*xv - BvgR*m3*PvL*sman*tman*uman*xv + 
     &     Bvm1L*m3*mq1*PvL*sman*tman*uman*xv - Bsm2R*PsR*(m2**2 - tman)*(m1**2 + m3**2 - tman)*(-1 + m4**2*xv) - 
     &     Bsm1R*PsR*(m1**2 + m3**2 - tman)*(-m4**2 + sman + tman + 2*uman - sman*tman*xv + m2**2*(-3 - m3**2*xv + sman*xv) + m3**2*(-3 + tman*xv))) - 
     &  2*A1vR*A2Rc*CLc*(2*Bsm2L*m1**2*m2**2*m3*PsL + Bsm2L*m2**2*m3**3*PsL + Bsm2L*m2**2*m3*m4**2*PsL - 4*BvgL*m1**2*m2**2*mq1*PsL - 
     &     6*Bvm1L*m1**2*m2**2*m3*mq1*PsL - 2*Bvm2L*m1**2*m2**2*m3*mq1*PsL - 6*BvgL*m1**2*m3**2*mq1*PsL - 4*BvgL*m2**2*m3**2*mq1*PsL - 
     &     6*Bvm1L*m1**2*m3**3*mq1*PsL - 3*Bvm1L*m2**2*m3**3*mq1*PsL - Bvm2L*m2**2*m3**3*mq1*PsL - 3*BvgL*m3**4*mq1*PsL - 3*Bvm1L*m3**5*mq1*PsL - 
     &     2*Bvm1L*m1**2*m3*m4**2*mq1*PsL - 3*Bvm1L*m2**2*m3*m4**2*mq1*PsL - Bvm2L*m2**2*m3*m4**2*mq1*PsL - 3*BvgL*m3**2*m4**2*mq1*PsL - 
     &     4*Bvm1L*m3**3*m4**2*mq1*PsL - Bvm1L*m3*m4**4*mq1*PsL + 3*Bvm1R*m1**2*m2**2*PsL*sman + Bvm2R*m1**2*m2**2*PsL*sman - 3*BvgR*m1**2*m3*PsL*sman - 
     &     2*BvgR*m2**2*m3*PsL*sman + 3*Bvm1R*m1**2*m3**2*PsL*sman + 3*Bvm1R*m2**2*m3**2*PsL*sman + Bvm2R*m2**2*m3**2*PsL*sman + 3*Bvm1R*m3**4*PsL*sman + 
     &     Bvm1R*m1**2*m4**2*PsL*sman - BvgR*m3*m4**2*PsL*sman + Bvm1R*m3**2*m4**2*PsL*sman - BvgL*m2**2*mq1*PsL*sman + 2*Bvm1L*m1**2*m3*mq1*PsL*sman + 
     &     2*BvgL*m3**2*mq1*PsL*sman + Bvm1L*m3**3*mq1*PsL*sman - BvgL*m4**2*mq1*PsL*sman + Bvm1L*m3*m4**2*mq1*PsL*sman + 3*Bsm1R*m1**2*m2**2*PvR*sman + 
     &     Bsm2R*m1**2*m2**2*PvR*sman + 4*BvgL*m1**2*m2**2*PvR*sman - 3*BsgR*m1**2*m3*PvR*sman - 2*BsgR*m2**2*m3*PvR*sman + 6*Bvm1L*m1**2*m2**2*m3*PvR*sman + 
     &     2*Bvm2L*m1**2*m2**2*m3*PvR*sman + 3*Bsm1R*m1**2*m3**2*PvR*sman + 6*BvgL*m1**2*m3**2*PvR*sman + 3*Bsm1R*m2**2*m3**2*PvR*sman + 
     &     Bsm2R*m2**2*m3**2*PvR*sman + 4*BvgL*m2**2*m3**2*PvR*sman + 6*Bvm1L*m1**2*m3**3*PvR*sman + 3*Bvm1L*m2**2*m3**3*PvR*sman + 
     &     Bvm2L*m2**2*m3**3*PvR*sman + 3*Bsm1R*m3**4*PvR*sman + 3*BvgL*m3**4*PvR*sman + 3*Bvm1L*m3**5*PvR*sman + Bsm1R*m1**2*m4**2*PvR*sman - 
     &     BsgR*m3*m4**2*PvR*sman + 2*Bvm1L*m1**2*m3*m4**2*PvR*sman + 3*Bvm1L*m2**2*m3*m4**2*PvR*sman + Bvm2L*m2**2*m3*m4**2*PvR*sman + 
     &     Bsm1R*m3**2*m4**2*PvR*sman + 3*BvgL*m3**2*m4**2*PvR*sman + 4*Bvm1L*m3**3*m4**2*PvR*sman + Bvm1L*m3*m4**4*PvR*sman - 
     &     3*Bvm1R*m1**2*m2**2*mq1*PvR*sman - Bvm2R*m1**2*m2**2*mq1*PvR*sman + 3*BvgR*m1**2*m3*mq1*PvR*sman + 2*BvgR*m2**2*m3*mq1*PvR*sman - 
     &     3*Bvm1R*m1**2*m3**2*mq1*PvR*sman - 3*Bvm1R*m2**2*m3**2*mq1*PvR*sman - Bvm2R*m2**2*m3**2*mq1*PvR*sman - 3*Bvm1R*m3**4*mq1*PvR*sman - 
     &     Bvm1R*m1**2*m4**2*mq1*PvR*sman + BvgR*m3*m4**2*mq1*PvR*sman - Bvm1R*m3**2*m4**2*mq1*PvR*sman - Bvm1R*m1**2*PsL*sman**2 + 2*BvgR*m3*PsL*sman**2 - 
     &     Bvm1R*m3**2*PsL*sman**2 + BvgL*mq1*PsL*sman**2 - Bsm1R*m1**2*PvR*sman**2 + BvgL*m2**2*PvR*sman**2 + 2*BsgR*m3*PvR*sman**2 - 
     &     2*Bvm1L*m1**2*m3*PvR*sman**2 - Bsm1R*m3**2*PvR*sman**2 - 2*BvgL*m3**2*PvR*sman**2 - Bvm1L*m3**3*PvR*sman**2 + BvgL*m4**2*PvR*sman**2 - 
     &     Bvm1L*m3*m4**2*PvR*sman**2 + Bvm1R*m1**2*mq1*PvR*sman**2 - 2*BvgR*m3*mq1*PvR*sman**2 + Bvm1R*m3**2*mq1*PvR*sman**2 - BvgL*PvR*sman**3 - 
     &     2*Bsm2L*m1**2*m3*PsL*tman - Bsm2L*m2**2*m3*PsL*tman - Bsm2L*m3**3*PsL*tman - Bsm2L*m3*m4**2*PsL*tman + BvgL*m1**2*mq1*PsL*tman + 
     &     3*BvgL*m2**2*mq1*PsL*tman + 2*Bvm1L*m1**2*m3*mq1*PsL*tman + 2*Bvm2L*m1**2*m3*mq1*PsL*tman + 3*Bvm1L*m2**2*m3*mq1*PsL*tman + 
     &     Bvm2L*m2**2*m3*mq1*PsL*tman + 4*BvgL*m3**2*mq1*PsL*tman + 4*Bvm1L*m3**3*mq1*PsL*tman + Bvm2L*m3**3*mq1*PsL*tman + BvgL*m4**2*mq1*PsL*tman + 
     &     2*Bvm1L*m3*m4**2*mq1*PsL*tman + Bvm2L*m3*m4**2*mq1*PsL*tman - Bvm1R*m1**2*PsL*sman*tman - Bvm2R*m1**2*PsL*sman*tman - 3*Bvm1R*m2**2*PsL*sman*tman - 
     &     Bvm2R*m2**2*PsL*sman*tman - 4*Bvm1R*m3**2*PsL*sman*tman - Bvm2R*m3**2*PsL*sman*tman - Bvm1R*m4**2*PsL*sman*tman - BvgL*mq1*PsL*sman*tman - 
     &     Bvm1L*m3*mq1*PsL*sman*tman - Bsm1R*m1**2*PvR*sman*tman - Bsm2R*m1**2*PvR*sman*tman - BvgL*m1**2*PvR*sman*tman - 3*Bsm1R*m2**2*PvR*sman*tman - 
     &     Bsm2R*m2**2*PvR*sman*tman - 3*BvgL*m2**2*PvR*sman*tman - 2*Bvm1L*m1**2*m3*PvR*sman*tman - 2*Bvm2L*m1**2*m3*PvR*sman*tman - 
     &     3*Bvm1L*m2**2*m3*PvR*sman*tman - Bvm2L*m2**2*m3*PvR*sman*tman - 4*Bsm1R*m3**2*PvR*sman*tman - Bsm2R*m3**2*PvR*sman*tman - 
     &     4*BvgL*m3**2*PvR*sman*tman - 4*Bvm1L*m3**3*PvR*sman*tman - Bvm2L*m3**3*PvR*sman*tman - Bsm1R*m4**2*PvR*sman*tman - BvgL*m4**2*PvR*sman*tman - 
     &     2*Bvm1L*m3*m4**2*PvR*sman*tman - Bvm2L*m3*m4**2*PvR*sman*tman + Bvm1R*m1**2*mq1*PvR*sman*tman + Bvm2R*m1**2*mq1*PvR*sman*tman + 
     &     3*Bvm1R*m2**2*mq1*PvR*sman*tman + Bvm2R*m2**2*mq1*PvR*sman*tman + 4*Bvm1R*m3**2*mq1*PvR*sman*tman + Bvm2R*m3**2*mq1*PvR*sman*tman + 
     &     Bvm1R*m4**2*mq1*PvR*sman*tman + Bvm1R*PsL*sman**2*tman + Bsm1R*PvR*sman**2*tman + BvgL*PvR*sman**2*tman + Bvm1L*m3*PvR*sman**2*tman - 
     &     Bvm1R*mq1*PvR*sman**2*tman + Bsm2L*m3*PsL*tman**2 - BvgL*mq1*PsL*tman**2 - Bvm1L*m3*mq1*PsL*tman**2 - Bvm2L*m3*mq1*PsL*tman**2 + 
     &     Bvm1R*PsL*sman*tman**2 + Bvm2R*PsL*sman*tman**2 + Bsm1R*PvR*sman*tman**2 + Bsm2R*PvR*sman*tman**2 + BvgL*PvR*sman*tman**2 + 
     &     Bvm1L*m3*PvR*sman*tman**2 + Bvm2L*m3*PvR*sman*tman**2 - Bvm1R*mq1*PvR*sman*tman**2 - Bvm2R*mq1*PvR*sman*tman**2 - Bsm2L*m2**2*m3*PsL*uman + 
     &     3*BvgL*m1**2*mq1*PsL*uman + BvgL*m2**2*mq1*PsL*uman + 4*Bvm1L*m1**2*m3*mq1*PsL*uman + 3*Bvm1L*m2**2*m3*mq1*PsL*uman + Bvm2L*m2**2*m3*mq1*PsL*uman + 
     &     4*BvgL*m3**2*mq1*PsL*uman + 5*Bvm1L*m3**3*mq1*PsL*uman + BvgL*m4**2*mq1*PsL*uman + 3*Bvm1L*m3*m4**2*mq1*PsL*uman - 2*Bvm1R*m1**2*PsL*sman*uman + 
     &     BvgR*m3*PsL*sman*uman - 2*Bvm1R*m3**2*PsL*sman*uman - Bvm1L*m3*mq1*PsL*sman*uman - 2*Bsm1R*m1**2*PvR*sman*uman - 3*BvgL*m1**2*PvR*sman*uman - 
     &     BvgL*m2**2*PvR*sman*uman + BsgR*m3*PvR*sman*uman - 4*Bvm1L*m1**2*m3*PvR*sman*uman - 3*Bvm1L*m2**2*m3*PvR*sman*uman - Bvm2L*m2**2*m3*PvR*sman*uman - 
     &     2*Bsm1R*m3**2*PvR*sman*uman - 4*BvgL*m3**2*PvR*sman*uman - 5*Bvm1L*m3**3*PvR*sman*uman - BvgL*m4**2*PvR*sman*uman - 
     &     3*Bvm1L*m3*m4**2*PvR*sman*uman + 2*Bvm1R*m1**2*mq1*PvR*sman*uman - BvgR*m3*mq1*PvR*sman*uman + 2*Bvm1R*m3**2*mq1*PvR*sman*uman + 
     &     Bvm1L*m3*PvR*sman**2*uman + Bsm2L*m3*PsL*tman*uman - 2*BvgL*mq1*PsL*tman*uman - 3*Bvm1L*m3*mq1*PsL*tman*uman - Bvm2L*m3*mq1*PsL*tman*uman + 
     &     2*Bvm1R*PsL*sman*tman*uman + 2*Bsm1R*PvR*sman*tman*uman + 2*BvgL*PvR*sman*tman*uman + 3*Bvm1L*m3*PvR*sman*tman*uman + Bvm2L*m3*PvR*sman*tman*uman - 
     &     2*Bvm1R*mq1*PvR*sman*tman*uman - BvgL*mq1*PsL*uman**2 - 2*Bvm1L*m3*mq1*PsL*uman**2 + BvgL*PvR*sman*uman**2 + 2*Bvm1L*m3*PvR*sman*uman**2 - 
     &     2*Bsm2L*m1**2*m2**2*m3*m4**2*PsL*xv - Bsm2L*m2**2*m3**3*m4**2*PsL*xv - Bsm2L*m2**2*m3*m4**4*PsL*xv - 2*BvgL*m1**2*m2**2*m3**2*mq1*PsL*xv - 
     &     2*Bvm1L*m1**2*m2**2*m3**3*mq1*PsL*xv - BvgL*m2**2*m3**4*mq1*PsL*xv - Bvm1L*m2**2*m3**5*mq1*PsL*xv + 2*Bvm2L*m1**2*m2**2*m3*m4**2*mq1*PsL*xv - 
     &     BvgL*m2**2*m3**2*m4**2*mq1*PsL*xv - Bvm1L*m2**2*m3**3*m4**2*mq1*PsL*xv + Bvm2L*m2**2*m3**3*m4**2*mq1*PsL*xv + Bvm2L*m2**2*m3*m4**4*mq1*PsL*xv - 
     &     BvgR*m1**2*m2**2*m3*PsL*sman*xv + Bvm1R*m1**2*m2**2*m3**2*PsL*sman*xv + Bvm1R*m2**2*m3**4*PsL*sman*xv - Bvm2R*m1**2*m2**2*m4**2*PsL*sman*xv - 
     &     BvgR*m2**2*m3*m4**2*PsL*sman*xv - Bvm2R*m2**2*m3**2*m4**2*PsL*sman*xv + BvgL*m1**2*m2**2*mq1*PsL*sman*xv + 2*Bvm1L*m1**2*m2**2*m3*mq1*PsL*sman*xv + 
     &     BvgL*m2**2*m3**2*mq1*PsL*sman*xv + Bvm1L*m2**2*m3**3*mq1*PsL*sman*xv + Bvm1L*m2**2*m3*m4**2*mq1*PsL*sman*xv - BsgR*m1**2*m2**2*m3*PvR*sman*xv + 
     &     Bsm1R*m1**2*m2**2*m3**2*PvR*sman*xv + 2*BvgL*m1**2*m2**2*m3**2*PvR*sman*xv + 2*Bvm1L*m1**2*m2**2*m3**3*PvR*sman*xv + 
     &     Bsm1R*m2**2*m3**4*PvR*sman*xv + BvgL*m2**2*m3**4*PvR*sman*xv + Bvm1L*m2**2*m3**5*PvR*sman*xv - Bsm2R*m1**2*m2**2*m4**2*PvR*sman*xv - 
     &     BsgR*m2**2*m3*m4**2*PvR*sman*xv - 2*Bvm2L*m1**2*m2**2*m3*m4**2*PvR*sman*xv - Bsm2R*m2**2*m3**2*m4**2*PvR*sman*xv + 
     &     BvgL*m2**2*m3**2*m4**2*PvR*sman*xv + Bvm1L*m2**2*m3**3*m4**2*PvR*sman*xv - Bvm2L*m2**2*m3**3*m4**2*PvR*sman*xv - Bvm2L*m2**2*m3*m4**4*PvR*sman*xv + 
     &     BvgR*m1**2*m2**2*m3*mq1*PvR*sman*xv - Bvm1R*m1**2*m2**2*m3**2*mq1*PvR*sman*xv - Bvm1R*m2**2*m3**4*mq1*PvR*sman*xv + 
     &     Bvm2R*m1**2*m2**2*m4**2*mq1*PvR*sman*xv + BvgR*m2**2*m3*m4**2*mq1*PvR*sman*xv + Bvm2R*m2**2*m3**2*m4**2*mq1*PvR*sman*xv - 
     &     Bvm1R*m1**2*m2**2*PsL*sman**2*xv - Bvm1R*m2**2*m3**2*PsL*sman**2*xv - Bsm1R*m1**2*m2**2*PvR*sman**2*xv - BvgL*m1**2*m2**2*PvR*sman**2*xv - 
     &     2*Bvm1L*m1**2*m2**2*m3*PvR*sman**2*xv - Bsm1R*m2**2*m3**2*PvR*sman**2*xv - BvgL*m2**2*m3**2*PvR*sman**2*xv - Bvm1L*m2**2*m3**3*PvR*sman**2*xv - 
     &     Bvm1L*m2**2*m3*m4**2*PvR*sman**2*xv + Bvm1R*m1**2*m2**2*mq1*PvR*sman**2*xv + Bvm1R*m2**2*m3**2*mq1*PvR*sman**2*xv + 
     &     2*Bsm2L*m1**2*m3*m4**2*PsL*tman*xv + Bsm2L*m2**2*m3*m4**2*PsL*tman*xv + Bsm2L*m3**3*m4**2*PsL*tman*xv + Bsm2L*m3*m4**4*PsL*tman*xv + 
     &     2*BvgL*m1**2*m3**2*mq1*PsL*tman*xv + BvgL*m2**2*m3**2*mq1*PsL*tman*xv + 2*Bvm1L*m1**2*m3**3*mq1*PsL*tman*xv + Bvm1L*m2**2*m3**3*mq1*PsL*tman*xv + 
     &     BvgL*m3**4*mq1*PsL*tman*xv + Bvm1L*m3**5*mq1*PsL*tman*xv - 2*Bvm2L*m1**2*m3*m4**2*mq1*PsL*tman*xv - Bvm2L*m2**2*m3*m4**2*mq1*PsL*tman*xv + 
     &     BvgL*m3**2*m4**2*mq1*PsL*tman*xv + Bvm1L*m3**3*m4**2*mq1*PsL*tman*xv - Bvm2L*m3**3*m4**2*mq1*PsL*tman*xv - Bvm2L*m3*m4**4*mq1*PsL*tman*xv + 
     &     BvgR*m1**2*m3*PsL*sman*tman*xv - Bvm1R*m1**2*m3**2*PsL*sman*tman*xv - Bvm1R*m2**2*m3**2*PsL*sman*tman*xv - Bvm1R*m3**4*PsL*sman*tman*xv + 
     &     Bvm2R*m1**2*m4**2*PsL*sman*tman*xv + Bvm2R*m2**2*m4**2*PsL*sman*tman*xv + BvgR*m3*m4**2*PsL*sman*tman*xv + Bvm2R*m3**2*m4**2*PsL*sman*tman*xv - 
     &     BvgL*m1**2*mq1*PsL*sman*tman*xv - BvgL*m2**2*mq1*PsL*sman*tman*xv - 2*Bvm1L*m1**2*m3*mq1*PsL*sman*tman*xv - Bvm1L*m2**2*m3*mq1*PsL*sman*tman*xv - 
     &     BvgL*m3**2*mq1*PsL*sman*tman*xv - Bvm1L*m3**3*mq1*PsL*sman*tman*xv - Bvm1L*m3*m4**2*mq1*PsL*sman*tman*xv + BsgR*m1**2*m3*PvR*sman*tman*xv - 
     &     Bsm1R*m1**2*m3**2*PvR*sman*tman*xv - 2*BvgL*m1**2*m3**2*PvR*sman*tman*xv - Bsm1R*m2**2*m3**2*PvR*sman*tman*xv - BvgL*m2**2*m3**2*PvR*sman*tman*xv - 
     &     2*Bvm1L*m1**2*m3**3*PvR*sman*tman*xv - Bvm1L*m2**2*m3**3*PvR*sman*tman*xv - Bsm1R*m3**4*PvR*sman*tman*xv - BvgL*m3**4*PvR*sman*tman*xv - 
     &     Bvm1L*m3**5*PvR*sman*tman*xv + Bsm2R*m1**2*m4**2*PvR*sman*tman*xv + Bsm2R*m2**2*m4**2*PvR*sman*tman*xv + BsgR*m3*m4**2*PvR*sman*tman*xv + 
     &     2*Bvm2L*m1**2*m3*m4**2*PvR*sman*tman*xv + Bvm2L*m2**2*m3*m4**2*PvR*sman*tman*xv + Bsm2R*m3**2*m4**2*PvR*sman*tman*xv - 
     &     BvgL*m3**2*m4**2*PvR*sman*tman*xv - Bvm1L*m3**3*m4**2*PvR*sman*tman*xv + Bvm2L*m3**3*m4**2*PvR*sman*tman*xv + Bvm2L*m3*m4**4*PvR*sman*tman*xv - 
     &     BvgR*m1**2*m3*mq1*PvR*sman*tman*xv + Bvm1R*m1**2*m3**2*mq1*PvR*sman*tman*xv + Bvm1R*m2**2*m3**2*mq1*PvR*sman*tman*xv + 
     &     Bvm1R*m3**4*mq1*PvR*sman*tman*xv - Bvm2R*m1**2*m4**2*mq1*PvR*sman*tman*xv - Bvm2R*m2**2*m4**2*mq1*PvR*sman*tman*xv - 
     &     BvgR*m3*m4**2*mq1*PvR*sman*tman*xv - Bvm2R*m3**2*m4**2*mq1*PvR*sman*tman*xv + Bvm1R*m1**2*PsL*sman**2*tman*xv + Bvm1R*m2**2*PsL*sman**2*tman*xv + 
     &     Bvm1R*m3**2*PsL*sman**2*tman*xv + Bsm1R*m1**2*PvR*sman**2*tman*xv + BvgL*m1**2*PvR*sman**2*tman*xv + Bsm1R*m2**2*PvR*sman**2*tman*xv + 
     &     BvgL*m2**2*PvR*sman**2*tman*xv + 2*Bvm1L*m1**2*m3*PvR*sman**2*tman*xv + Bvm1L*m2**2*m3*PvR*sman**2*tman*xv + Bsm1R*m3**2*PvR*sman**2*tman*xv + 
     &     BvgL*m3**2*PvR*sman**2*tman*xv + Bvm1L*m3**3*PvR*sman**2*tman*xv + Bvm1L*m3*m4**2*PvR*sman**2*tman*xv - Bvm1R*m1**2*mq1*PvR*sman**2*tman*xv - 
     &     Bvm1R*m2**2*mq1*PvR*sman**2*tman*xv - Bvm1R*m3**2*mq1*PvR*sman**2*tman*xv - Bsm2L*m3*m4**2*PsL*tman**2*xv - BvgL*m3**2*mq1*PsL*tman**2*xv - 
     &     Bvm1L*m3**3*mq1*PsL*tman**2*xv + Bvm2L*m3*m4**2*mq1*PsL*tman**2*xv + Bvm1R*m3**2*PsL*sman*tman**2*xv - Bvm2R*m4**2*PsL*sman*tman**2*xv + 
     &     BvgL*mq1*PsL*sman*tman**2*xv + Bvm1L*m3*mq1*PsL*sman*tman**2*xv + Bsm1R*m3**2*PvR*sman*tman**2*xv + BvgL*m3**2*PvR*sman*tman**2*xv + 
     &     Bvm1L*m3**3*PvR*sman*tman**2*xv - Bsm2R*m4**2*PvR*sman*tman**2*xv - Bvm2L*m3*m4**2*PvR*sman*tman**2*xv - Bvm1R*m3**2*mq1*PvR*sman*tman**2*xv + 
     &     Bvm2R*m4**2*mq1*PvR*sman*tman**2*xv - Bvm1R*PsL*sman**2*tman**2*xv - Bsm1R*PvR*sman**2*tman**2*xv - BvgL*PvR*sman**2*tman**2*xv - 
     &     Bvm1L*m3*PvR*sman**2*tman**2*xv + Bvm1R*mq1*PvR*sman**2*tman**2*xv + Bsm2L*m2**2*m3*m4**2*PsL*uman*xv + BvgL*m2**2*m3**2*mq1*PsL*uman*xv + 
     &     Bvm1L*m2**2*m3**3*mq1*PsL*uman*xv - Bvm2L*m2**2*m3*m4**2*mq1*PsL*uman*xv + BvgR*m2**2*m3*PsL*sman*uman*xv - Bvm1L*m2**2*m3*mq1*PsL*sman*uman*xv + 
     &     BsgR*m2**2*m3*PvR*sman*uman*xv - BvgL*m2**2*m3**2*PvR*sman*uman*xv - Bvm1L*m2**2*m3**3*PvR*sman*uman*xv + Bvm2L*m2**2*m3*m4**2*PvR*sman*uman*xv - 
     &     BvgR*m2**2*m3*mq1*PvR*sman*uman*xv + Bvm1L*m2**2*m3*PvR*sman**2*uman*xv - Bsm2L*m3*m4**2*PsL*tman*uman*xv - BvgL*m3**2*mq1*PsL*tman*uman*xv - 
     &     Bvm1L*m3**3*mq1*PsL*tman*uman*xv + Bvm2L*m3*m4**2*mq1*PsL*tman*uman*xv - BvgR*m3*PsL*sman*tman*uman*xv + Bvm1L*m3*mq1*PsL*sman*tman*uman*xv - 
     &     BsgR*m3*PvR*sman*tman*uman*xv + BvgL*m3**2*PvR*sman*tman*uman*xv + Bvm1L*m3**3*PvR*sman*tman*uman*xv - Bvm2L*m3*m4**2*PvR*sman*tman*uman*xv + 
     &     BvgR*m3*mq1*PvR*sman*tman*uman*xv - Bvm1L*m3*PvR*sman**2*tman*uman*xv - 
     &     Bsm1L*m3*PsL*(2*m1**2 + m3**2 + m4**2 - tman - uman)*
     &      (-m4**2 + sman + tman + 2*uman - sman*tman*xv + m2**2*(-3 - m3**2*xv + sman*xv) + m3**2*(-3 + tman*xv)) + 
     &     BsgL*PsL*(3*m3**4 + 3*m3**2*m4**2 - 2*m3**2*sman + m4**2*sman - sman**2 - 4*m3**2*tman - m4**2*tman + sman*tman + tman**2 - 4*m3**2*uman - 
     &        m4**2*uman + 2*tman*uman + uman**2 - m3**4*tman*xv - m3**2*m4**2*tman*xv + m3**2*sman*tman*xv + m3**2*tman**2*xv - sman*tman**2*xv + 
     &        m3**2*tman*uman*xv + m1**2*(-tman - 3*uman + sman*tman*xv + m2**2*(4 + 2*m3**2*xv - sman*xv) + m3**2*(6 - 2*tman*xv)) + 
     &        m2**2*(sman - 3*tman - uman + m3**4*xv + sman*tman*xv + m3**2*(4 + m4**2*xv - sman*xv - tman*xv - uman*xv)))) - 
     &  2*A1vL*A2Lc*CLc*(2*Bsm2R*m1**2*m2**2*m3*PsR + Bsm2R*m2**2*m3**3*PsR + Bsm2R*m2**2*m3*m4**2*PsR - 4*BvgR*m1**2*m2**2*mq1*PsR - 
     &     6*Bvm1R*m1**2*m2**2*m3*mq1*PsR - 2*Bvm2R*m1**2*m2**2*m3*mq1*PsR - 6*BvgR*m1**2*m3**2*mq1*PsR - 4*BvgR*m2**2*m3**2*mq1*PsR - 
     &     6*Bvm1R*m1**2*m3**3*mq1*PsR - 3*Bvm1R*m2**2*m3**3*mq1*PsR - Bvm2R*m2**2*m3**3*mq1*PsR - 3*BvgR*m3**4*mq1*PsR - 3*Bvm1R*m3**5*mq1*PsR - 
     &     2*Bvm1R*m1**2*m3*m4**2*mq1*PsR - 3*Bvm1R*m2**2*m3*m4**2*mq1*PsR - Bvm2R*m2**2*m3*m4**2*mq1*PsR - 3*BvgR*m3**2*m4**2*mq1*PsR - 
     &     4*Bvm1R*m3**3*m4**2*mq1*PsR - Bvm1R*m3*m4**4*mq1*PsR + 3*Bvm1L*m1**2*m2**2*PsR*sman + Bvm2L*m1**2*m2**2*PsR*sman - 3*BvgL*m1**2*m3*PsR*sman - 
     &     2*BvgL*m2**2*m3*PsR*sman + 3*Bvm1L*m1**2*m3**2*PsR*sman + 3*Bvm1L*m2**2*m3**2*PsR*sman + Bvm2L*m2**2*m3**2*PsR*sman + 3*Bvm1L*m3**4*PsR*sman + 
     &     Bvm1L*m1**2*m4**2*PsR*sman - BvgL*m3*m4**2*PsR*sman + Bvm1L*m3**2*m4**2*PsR*sman - BvgR*m2**2*mq1*PsR*sman + 2*Bvm1R*m1**2*m3*mq1*PsR*sman + 
     &     2*BvgR*m3**2*mq1*PsR*sman + Bvm1R*m3**3*mq1*PsR*sman - BvgR*m4**2*mq1*PsR*sman + Bvm1R*m3*m4**2*mq1*PsR*sman + 3*Bsm1L*m1**2*m2**2*PvL*sman + 
     &     Bsm2L*m1**2*m2**2*PvL*sman + 4*BvgR*m1**2*m2**2*PvL*sman - 3*BsgL*m1**2*m3*PvL*sman - 2*BsgL*m2**2*m3*PvL*sman + 6*Bvm1R*m1**2*m2**2*m3*PvL*sman + 
     &     2*Bvm2R*m1**2*m2**2*m3*PvL*sman + 3*Bsm1L*m1**2*m3**2*PvL*sman + 6*BvgR*m1**2*m3**2*PvL*sman + 3*Bsm1L*m2**2*m3**2*PvL*sman + 
     &     Bsm2L*m2**2*m3**2*PvL*sman + 4*BvgR*m2**2*m3**2*PvL*sman + 6*Bvm1R*m1**2*m3**3*PvL*sman + 3*Bvm1R*m2**2*m3**3*PvL*sman + 
     &     Bvm2R*m2**2*m3**3*PvL*sman + 3*Bsm1L*m3**4*PvL*sman + 3*BvgR*m3**4*PvL*sman + 3*Bvm1R*m3**5*PvL*sman + Bsm1L*m1**2*m4**2*PvL*sman - 
     &     BsgL*m3*m4**2*PvL*sman + 2*Bvm1R*m1**2*m3*m4**2*PvL*sman + 3*Bvm1R*m2**2*m3*m4**2*PvL*sman + Bvm2R*m2**2*m3*m4**2*PvL*sman + 
     &     Bsm1L*m3**2*m4**2*PvL*sman + 3*BvgR*m3**2*m4**2*PvL*sman + 4*Bvm1R*m3**3*m4**2*PvL*sman + Bvm1R*m3*m4**4*PvL*sman - 
     &     3*Bvm1L*m1**2*m2**2*mq1*PvL*sman - Bvm2L*m1**2*m2**2*mq1*PvL*sman + 3*BvgL*m1**2*m3*mq1*PvL*sman + 2*BvgL*m2**2*m3*mq1*PvL*sman - 
     &     3*Bvm1L*m1**2*m3**2*mq1*PvL*sman - 3*Bvm1L*m2**2*m3**2*mq1*PvL*sman - Bvm2L*m2**2*m3**2*mq1*PvL*sman - 3*Bvm1L*m3**4*mq1*PvL*sman - 
     &     Bvm1L*m1**2*m4**2*mq1*PvL*sman + BvgL*m3*m4**2*mq1*PvL*sman - Bvm1L*m3**2*m4**2*mq1*PvL*sman - Bvm1L*m1**2*PsR*sman**2 + 2*BvgL*m3*PsR*sman**2 - 
     &     Bvm1L*m3**2*PsR*sman**2 + BvgR*mq1*PsR*sman**2 - Bsm1L*m1**2*PvL*sman**2 + BvgR*m2**2*PvL*sman**2 + 2*BsgL*m3*PvL*sman**2 - 
     &     2*Bvm1R*m1**2*m3*PvL*sman**2 - Bsm1L*m3**2*PvL*sman**2 - 2*BvgR*m3**2*PvL*sman**2 - Bvm1R*m3**3*PvL*sman**2 + BvgR*m4**2*PvL*sman**2 - 
     &     Bvm1R*m3*m4**2*PvL*sman**2 + Bvm1L*m1**2*mq1*PvL*sman**2 - 2*BvgL*m3*mq1*PvL*sman**2 + Bvm1L*m3**2*mq1*PvL*sman**2 - BvgR*PvL*sman**3 - 
     &     2*Bsm2R*m1**2*m3*PsR*tman - Bsm2R*m2**2*m3*PsR*tman - Bsm2R*m3**3*PsR*tman - Bsm2R*m3*m4**2*PsR*tman + BvgR*m1**2*mq1*PsR*tman + 
     &     3*BvgR*m2**2*mq1*PsR*tman + 2*Bvm1R*m1**2*m3*mq1*PsR*tman + 2*Bvm2R*m1**2*m3*mq1*PsR*tman + 3*Bvm1R*m2**2*m3*mq1*PsR*tman + 
     &     Bvm2R*m2**2*m3*mq1*PsR*tman + 4*BvgR*m3**2*mq1*PsR*tman + 4*Bvm1R*m3**3*mq1*PsR*tman + Bvm2R*m3**3*mq1*PsR*tman + BvgR*m4**2*mq1*PsR*tman + 
     &     2*Bvm1R*m3*m4**2*mq1*PsR*tman + Bvm2R*m3*m4**2*mq1*PsR*tman - Bvm1L*m1**2*PsR*sman*tman - Bvm2L*m1**2*PsR*sman*tman - 3*Bvm1L*m2**2*PsR*sman*tman - 
     &     Bvm2L*m2**2*PsR*sman*tman - 4*Bvm1L*m3**2*PsR*sman*tman - Bvm2L*m3**2*PsR*sman*tman - Bvm1L*m4**2*PsR*sman*tman - BvgR*mq1*PsR*sman*tman - 
     &     Bvm1R*m3*mq1*PsR*sman*tman - Bsm1L*m1**2*PvL*sman*tman - Bsm2L*m1**2*PvL*sman*tman - BvgR*m1**2*PvL*sman*tman - 3*Bsm1L*m2**2*PvL*sman*tman - 
     &     Bsm2L*m2**2*PvL*sman*tman - 3*BvgR*m2**2*PvL*sman*tman - 2*Bvm1R*m1**2*m3*PvL*sman*tman - 2*Bvm2R*m1**2*m3*PvL*sman*tman - 
     &     3*Bvm1R*m2**2*m3*PvL*sman*tman - Bvm2R*m2**2*m3*PvL*sman*tman - 4*Bsm1L*m3**2*PvL*sman*tman - Bsm2L*m3**2*PvL*sman*tman - 
     &     4*BvgR*m3**2*PvL*sman*tman - 4*Bvm1R*m3**3*PvL*sman*tman - Bvm2R*m3**3*PvL*sman*tman - Bsm1L*m4**2*PvL*sman*tman - BvgR*m4**2*PvL*sman*tman - 
     &     2*Bvm1R*m3*m4**2*PvL*sman*tman - Bvm2R*m3*m4**2*PvL*sman*tman + Bvm1L*m1**2*mq1*PvL*sman*tman + Bvm2L*m1**2*mq1*PvL*sman*tman + 
     &     3*Bvm1L*m2**2*mq1*PvL*sman*tman + Bvm2L*m2**2*mq1*PvL*sman*tman + 4*Bvm1L*m3**2*mq1*PvL*sman*tman + Bvm2L*m3**2*mq1*PvL*sman*tman + 
     &     Bvm1L*m4**2*mq1*PvL*sman*tman + Bvm1L*PsR*sman**2*tman + Bsm1L*PvL*sman**2*tman + BvgR*PvL*sman**2*tman + Bvm1R*m3*PvL*sman**2*tman - 
     &     Bvm1L*mq1*PvL*sman**2*tman + Bsm2R*m3*PsR*tman**2 - BvgR*mq1*PsR*tman**2 - Bvm1R*m3*mq1*PsR*tman**2 - Bvm2R*m3*mq1*PsR*tman**2 + 
     &     Bvm1L*PsR*sman*tman**2 + Bvm2L*PsR*sman*tman**2 + Bsm1L*PvL*sman*tman**2 + Bsm2L*PvL*sman*tman**2 + BvgR*PvL*sman*tman**2 + 
     &     Bvm1R*m3*PvL*sman*tman**2 + Bvm2R*m3*PvL*sman*tman**2 - Bvm1L*mq1*PvL*sman*tman**2 - Bvm2L*mq1*PvL*sman*tman**2 - Bsm2R*m2**2*m3*PsR*uman + 
     &     3*BvgR*m1**2*mq1*PsR*uman + BvgR*m2**2*mq1*PsR*uman + 4*Bvm1R*m1**2*m3*mq1*PsR*uman + 3*Bvm1R*m2**2*m3*mq1*PsR*uman + Bvm2R*m2**2*m3*mq1*PsR*uman + 
     &     4*BvgR*m3**2*mq1*PsR*uman + 5*Bvm1R*m3**3*mq1*PsR*uman + BvgR*m4**2*mq1*PsR*uman + 3*Bvm1R*m3*m4**2*mq1*PsR*uman - 2*Bvm1L*m1**2*PsR*sman*uman + 
     &     BvgL*m3*PsR*sman*uman - 2*Bvm1L*m3**2*PsR*sman*uman - Bvm1R*m3*mq1*PsR*sman*uman - 2*Bsm1L*m1**2*PvL*sman*uman - 3*BvgR*m1**2*PvL*sman*uman - 
     &     BvgR*m2**2*PvL*sman*uman + BsgL*m3*PvL*sman*uman - 4*Bvm1R*m1**2*m3*PvL*sman*uman - 3*Bvm1R*m2**2*m3*PvL*sman*uman - Bvm2R*m2**2*m3*PvL*sman*uman - 
     &     2*Bsm1L*m3**2*PvL*sman*uman - 4*BvgR*m3**2*PvL*sman*uman - 5*Bvm1R*m3**3*PvL*sman*uman - BvgR*m4**2*PvL*sman*uman - 
     &     3*Bvm1R*m3*m4**2*PvL*sman*uman + 2*Bvm1L*m1**2*mq1*PvL*sman*uman - BvgL*m3*mq1*PvL*sman*uman + 2*Bvm1L*m3**2*mq1*PvL*sman*uman + 
     &     Bvm1R*m3*PvL*sman**2*uman + Bsm2R*m3*PsR*tman*uman - 2*BvgR*mq1*PsR*tman*uman - 3*Bvm1R*m3*mq1*PsR*tman*uman - Bvm2R*m3*mq1*PsR*tman*uman + 
     &     2*Bvm1L*PsR*sman*tman*uman + 2*Bsm1L*PvL*sman*tman*uman + 2*BvgR*PvL*sman*tman*uman + 3*Bvm1R*m3*PvL*sman*tman*uman + Bvm2R*m3*PvL*sman*tman*uman - 
     &     2*Bvm1L*mq1*PvL*sman*tman*uman - BvgR*mq1*PsR*uman**2 - 2*Bvm1R*m3*mq1*PsR*uman**2 + BvgR*PvL*sman*uman**2 + 2*Bvm1R*m3*PvL*sman*uman**2 - 
     &     2*Bsm2R*m1**2*m2**2*m3*m4**2*PsR*xv - Bsm2R*m2**2*m3**3*m4**2*PsR*xv - Bsm2R*m2**2*m3*m4**4*PsR*xv - 2*BvgR*m1**2*m2**2*m3**2*mq1*PsR*xv - 
     &     2*Bvm1R*m1**2*m2**2*m3**3*mq1*PsR*xv - BvgR*m2**2*m3**4*mq1*PsR*xv - Bvm1R*m2**2*m3**5*mq1*PsR*xv + 2*Bvm2R*m1**2*m2**2*m3*m4**2*mq1*PsR*xv - 
     &     BvgR*m2**2*m3**2*m4**2*mq1*PsR*xv - Bvm1R*m2**2*m3**3*m4**2*mq1*PsR*xv + Bvm2R*m2**2*m3**3*m4**2*mq1*PsR*xv + Bvm2R*m2**2*m3*m4**4*mq1*PsR*xv - 
     &     BvgL*m1**2*m2**2*m3*PsR*sman*xv + Bvm1L*m1**2*m2**2*m3**2*PsR*sman*xv + Bvm1L*m2**2*m3**4*PsR*sman*xv - Bvm2L*m1**2*m2**2*m4**2*PsR*sman*xv - 
     &     BvgL*m2**2*m3*m4**2*PsR*sman*xv - Bvm2L*m2**2*m3**2*m4**2*PsR*sman*xv + BvgR*m1**2*m2**2*mq1*PsR*sman*xv + 2*Bvm1R*m1**2*m2**2*m3*mq1*PsR*sman*xv + 
     &     BvgR*m2**2*m3**2*mq1*PsR*sman*xv + Bvm1R*m2**2*m3**3*mq1*PsR*sman*xv + Bvm1R*m2**2*m3*m4**2*mq1*PsR*sman*xv - BsgL*m1**2*m2**2*m3*PvL*sman*xv + 
     &     Bsm1L*m1**2*m2**2*m3**2*PvL*sman*xv + 2*BvgR*m1**2*m2**2*m3**2*PvL*sman*xv + 2*Bvm1R*m1**2*m2**2*m3**3*PvL*sman*xv + 
     &     Bsm1L*m2**2*m3**4*PvL*sman*xv + BvgR*m2**2*m3**4*PvL*sman*xv + Bvm1R*m2**2*m3**5*PvL*sman*xv - Bsm2L*m1**2*m2**2*m4**2*PvL*sman*xv - 
     &     BsgL*m2**2*m3*m4**2*PvL*sman*xv - 2*Bvm2R*m1**2*m2**2*m3*m4**2*PvL*sman*xv - Bsm2L*m2**2*m3**2*m4**2*PvL*sman*xv + 
     &     BvgR*m2**2*m3**2*m4**2*PvL*sman*xv + Bvm1R*m2**2*m3**3*m4**2*PvL*sman*xv - Bvm2R*m2**2*m3**3*m4**2*PvL*sman*xv - Bvm2R*m2**2*m3*m4**4*PvL*sman*xv + 
     &     BvgL*m1**2*m2**2*m3*mq1*PvL*sman*xv - Bvm1L*m1**2*m2**2*m3**2*mq1*PvL*sman*xv - Bvm1L*m2**2*m3**4*mq1*PvL*sman*xv + 
     &     Bvm2L*m1**2*m2**2*m4**2*mq1*PvL*sman*xv + BvgL*m2**2*m3*m4**2*mq1*PvL*sman*xv + Bvm2L*m2**2*m3**2*m4**2*mq1*PvL*sman*xv - 
     &     Bvm1L*m1**2*m2**2*PsR*sman**2*xv - Bvm1L*m2**2*m3**2*PsR*sman**2*xv - Bsm1L*m1**2*m2**2*PvL*sman**2*xv - BvgR*m1**2*m2**2*PvL*sman**2*xv - 
     &     2*Bvm1R*m1**2*m2**2*m3*PvL*sman**2*xv - Bsm1L*m2**2*m3**2*PvL*sman**2*xv - BvgR*m2**2*m3**2*PvL*sman**2*xv - Bvm1R*m2**2*m3**3*PvL*sman**2*xv - 
     &     Bvm1R*m2**2*m3*m4**2*PvL*sman**2*xv + Bvm1L*m1**2*m2**2*mq1*PvL*sman**2*xv + Bvm1L*m2**2*m3**2*mq1*PvL*sman**2*xv + 
     &     2*Bsm2R*m1**2*m3*m4**2*PsR*tman*xv + Bsm2R*m2**2*m3*m4**2*PsR*tman*xv + Bsm2R*m3**3*m4**2*PsR*tman*xv + Bsm2R*m3*m4**4*PsR*tman*xv + 
     &     2*BvgR*m1**2*m3**2*mq1*PsR*tman*xv + BvgR*m2**2*m3**2*mq1*PsR*tman*xv + 2*Bvm1R*m1**2*m3**3*mq1*PsR*tman*xv + Bvm1R*m2**2*m3**3*mq1*PsR*tman*xv + 
     &     BvgR*m3**4*mq1*PsR*tman*xv + Bvm1R*m3**5*mq1*PsR*tman*xv - 2*Bvm2R*m1**2*m3*m4**2*mq1*PsR*tman*xv - Bvm2R*m2**2*m3*m4**2*mq1*PsR*tman*xv + 
     &     BvgR*m3**2*m4**2*mq1*PsR*tman*xv + Bvm1R*m3**3*m4**2*mq1*PsR*tman*xv - Bvm2R*m3**3*m4**2*mq1*PsR*tman*xv - Bvm2R*m3*m4**4*mq1*PsR*tman*xv + 
     &     BvgL*m1**2*m3*PsR*sman*tman*xv - Bvm1L*m1**2*m3**2*PsR*sman*tman*xv - Bvm1L*m2**2*m3**2*PsR*sman*tman*xv - Bvm1L*m3**4*PsR*sman*tman*xv + 
     &     Bvm2L*m1**2*m4**2*PsR*sman*tman*xv + Bvm2L*m2**2*m4**2*PsR*sman*tman*xv + BvgL*m3*m4**2*PsR*sman*tman*xv + Bvm2L*m3**2*m4**2*PsR*sman*tman*xv - 
     &     BvgR*m1**2*mq1*PsR*sman*tman*xv - BvgR*m2**2*mq1*PsR*sman*tman*xv - 2*Bvm1R*m1**2*m3*mq1*PsR*sman*tman*xv - Bvm1R*m2**2*m3*mq1*PsR*sman*tman*xv - 
     &     BvgR*m3**2*mq1*PsR*sman*tman*xv - Bvm1R*m3**3*mq1*PsR*sman*tman*xv - Bvm1R*m3*m4**2*mq1*PsR*sman*tman*xv + BsgL*m1**2*m3*PvL*sman*tman*xv - 
     &     Bsm1L*m1**2*m3**2*PvL*sman*tman*xv - 2*BvgR*m1**2*m3**2*PvL*sman*tman*xv - Bsm1L*m2**2*m3**2*PvL*sman*tman*xv - BvgR*m2**2*m3**2*PvL*sman*tman*xv - 
     &     2*Bvm1R*m1**2*m3**3*PvL*sman*tman*xv - Bvm1R*m2**2*m3**3*PvL*sman*tman*xv - Bsm1L*m3**4*PvL*sman*tman*xv - BvgR*m3**4*PvL*sman*tman*xv - 
     &     Bvm1R*m3**5*PvL*sman*tman*xv + Bsm2L*m1**2*m4**2*PvL*sman*tman*xv + Bsm2L*m2**2*m4**2*PvL*sman*tman*xv + BsgL*m3*m4**2*PvL*sman*tman*xv + 
     &     2*Bvm2R*m1**2*m3*m4**2*PvL*sman*tman*xv + Bvm2R*m2**2*m3*m4**2*PvL*sman*tman*xv + Bsm2L*m3**2*m4**2*PvL*sman*tman*xv - 
     &     BvgR*m3**2*m4**2*PvL*sman*tman*xv - Bvm1R*m3**3*m4**2*PvL*sman*tman*xv + Bvm2R*m3**3*m4**2*PvL*sman*tman*xv + Bvm2R*m3*m4**4*PvL*sman*tman*xv - 
     &     BvgL*m1**2*m3*mq1*PvL*sman*tman*xv + Bvm1L*m1**2*m3**2*mq1*PvL*sman*tman*xv + Bvm1L*m2**2*m3**2*mq1*PvL*sman*tman*xv + 
     &     Bvm1L*m3**4*mq1*PvL*sman*tman*xv - Bvm2L*m1**2*m4**2*mq1*PvL*sman*tman*xv - Bvm2L*m2**2*m4**2*mq1*PvL*sman*tman*xv - 
     &     BvgL*m3*m4**2*mq1*PvL*sman*tman*xv - Bvm2L*m3**2*m4**2*mq1*PvL*sman*tman*xv + Bvm1L*m1**2*PsR*sman**2*tman*xv + Bvm1L*m2**2*PsR*sman**2*tman*xv + 
     &     Bvm1L*m3**2*PsR*sman**2*tman*xv + Bsm1L*m1**2*PvL*sman**2*tman*xv + BvgR*m1**2*PvL*sman**2*tman*xv + Bsm1L*m2**2*PvL*sman**2*tman*xv + 
     &     BvgR*m2**2*PvL*sman**2*tman*xv + 2*Bvm1R*m1**2*m3*PvL*sman**2*tman*xv + Bvm1R*m2**2*m3*PvL*sman**2*tman*xv + Bsm1L*m3**2*PvL*sman**2*tman*xv + 
     &     BvgR*m3**2*PvL*sman**2*tman*xv + Bvm1R*m3**3*PvL*sman**2*tman*xv + Bvm1R*m3*m4**2*PvL*sman**2*tman*xv - Bvm1L*m1**2*mq1*PvL*sman**2*tman*xv - 
     &     Bvm1L*m2**2*mq1*PvL*sman**2*tman*xv - Bvm1L*m3**2*mq1*PvL*sman**2*tman*xv - Bsm2R*m3*m4**2*PsR*tman**2*xv - BvgR*m3**2*mq1*PsR*tman**2*xv - 
     &     Bvm1R*m3**3*mq1*PsR*tman**2*xv + Bvm2R*m3*m4**2*mq1*PsR*tman**2*xv + Bvm1L*m3**2*PsR*sman*tman**2*xv - Bvm2L*m4**2*PsR*sman*tman**2*xv + 
     &     BvgR*mq1*PsR*sman*tman**2*xv + Bvm1R*m3*mq1*PsR*sman*tman**2*xv + Bsm1L*m3**2*PvL*sman*tman**2*xv + BvgR*m3**2*PvL*sman*tman**2*xv + 
     &     Bvm1R*m3**3*PvL*sman*tman**2*xv - Bsm2L*m4**2*PvL*sman*tman**2*xv - Bvm2R*m3*m4**2*PvL*sman*tman**2*xv - Bvm1L*m3**2*mq1*PvL*sman*tman**2*xv + 
     &     Bvm2L*m4**2*mq1*PvL*sman*tman**2*xv - Bvm1L*PsR*sman**2*tman**2*xv - Bsm1L*PvL*sman**2*tman**2*xv - BvgR*PvL*sman**2*tman**2*xv - 
     &     Bvm1R*m3*PvL*sman**2*tman**2*xv + Bvm1L*mq1*PvL*sman**2*tman**2*xv + Bsm2R*m2**2*m3*m4**2*PsR*uman*xv + BvgR*m2**2*m3**2*mq1*PsR*uman*xv + 
     &     Bvm1R*m2**2*m3**3*mq1*PsR*uman*xv - Bvm2R*m2**2*m3*m4**2*mq1*PsR*uman*xv + BvgL*m2**2*m3*PsR*sman*uman*xv - Bvm1R*m2**2*m3*mq1*PsR*sman*uman*xv + 
     &     BsgL*m2**2*m3*PvL*sman*uman*xv - BvgR*m2**2*m3**2*PvL*sman*uman*xv - Bvm1R*m2**2*m3**3*PvL*sman*uman*xv + Bvm2R*m2**2*m3*m4**2*PvL*sman*uman*xv - 
     &     BvgL*m2**2*m3*mq1*PvL*sman*uman*xv + Bvm1R*m2**2*m3*PvL*sman**2*uman*xv - Bsm2R*m3*m4**2*PsR*tman*uman*xv - BvgR*m3**2*mq1*PsR*tman*uman*xv - 
     &     Bvm1R*m3**3*mq1*PsR*tman*uman*xv + Bvm2R*m3*m4**2*mq1*PsR*tman*uman*xv - BvgL*m3*PsR*sman*tman*uman*xv + Bvm1R*m3*mq1*PsR*sman*tman*uman*xv - 
     &     BsgL*m3*PvL*sman*tman*uman*xv + BvgR*m3**2*PvL*sman*tman*uman*xv + Bvm1R*m3**3*PvL*sman*tman*uman*xv - Bvm2R*m3*m4**2*PvL*sman*tman*uman*xv + 
     &     BvgL*m3*mq1*PvL*sman*tman*uman*xv - Bvm1R*m3*PvL*sman**2*tman*uman*xv - 
     &     Bsm1R*m3*PsR*(2*m1**2 + m3**2 + m4**2 - tman - uman)*
     &      (-m4**2 + sman + tman + 2*uman - sman*tman*xv + m2**2*(-3 - m3**2*xv + sman*xv) + m3**2*(-3 + tman*xv)) + 
     &     BsgR*PsR*(3*m3**4 + 3*m3**2*m4**2 - 2*m3**2*sman + m4**2*sman - sman**2 - 4*m3**2*tman - m4**2*tman + sman*tman + tman**2 - 4*m3**2*uman - 
     &        m4**2*uman + 2*tman*uman + uman**2 - m3**4*tman*xv - m3**2*m4**2*tman*xv + m3**2*sman*tman*xv + m3**2*tman**2*xv - sman*tman**2*xv + 
     &        m3**2*tman*uman*xv + m1**2*(-tman - 3*uman + sman*tman*xv + m2**2*(4 + 2*m3**2*xv - sman*xv) + m3**2*(6 - 2*tman*xv)) + 
     &        m2**2*(sman - 3*tman - uman + m3**4*xv + sman*tman*xv + m3**2*(4 + m4**2*xv - sman*xv - tman*xv - uman*xv)))) 
     &      ))/((sman - mq1**2 + cI*mq1*wquark)*(tman - mq2c**2)))


        end

! --- Msu  ---

        double precision function MsuV(A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A3Lc,A3Rc,DLc,DRc,mq1,mq3c,xv,wquark)

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

        double complex A1sL,A1sR,A1vL,A1vR,PsL,PsR,PvL,PvR,BsgL,BsgR,Bsm1L,Bsm1R,Bsm2L,Bsm2R,BvgL,BvgR,Bvm1L,Bvm1R,Bvm2L,Bvm2R,A3Lc,A3Rc,DLc,DRc
        double precision mq1,mq3c,xv,wquark

        MsuV = 0d0

        MsuV = dreal((0.5*(A1sR*A3Rc*(8*BvgL*DRc*m1**2*m2**2*PsR - 8*BsgR*DRc*m1**2*m3*PsR - 4*BsgR*DRc*m2**2*m3*PsR + 12*Bvm1L*DRc*m1**2*m2**2*m3*PsR + 
     &     4*Bvm2L*DRc*m1**2*m2**2*m3*PsR + 8*BvgL*DRc*m1**2*m3**2*PsR + 12*BvgL*DRc*m2**2*m3**2*PsR - 4*BsgR*DRc*m3**3*PsR + 8*Bvm1L*DRc*m2**2*m3**3*PsR + 
     &     4*BvgL*DRc*m3**4*PsR - 4*BvgL*DRc*m2**2*m4**2*PsR + 8*Bvm1L*DRc*m1**2*m3*m4**2*PsR + 8*Bvm2L*DRc*m1**2*m3*m4**2*PsR + 
     &     4*Bvm1L*DRc*m2**2*m3*m4**2*PsR + 4*Bvm2L*DRc*m2**2*m3*m4**2*PsR + 4*BvgL*DRc*m3**2*m4**2*PsR + 4*Bvm1L*DRc*m3**3*m4**2*PsR + 
     &     4*Bvm2L*DRc*m3**3*m4**2*PsR + 4*Bvm1L*DRc*m3*m4**4*PsR + 4*Bvm2L*DRc*m3*m4**4*PsR - 6*Bvm1R*DRc*m1**2*m2**2*mq1*PsR - 
     &     2*Bvm2R*DRc*m1**2*m2**2*mq1*PsR - 8*BvgR*DRc*m1**2*m3*mq1*PsR - 4*BvgR*DRc*m2**2*m3*mq1*PsR - 2*Bvm1R*DRc*m1**2*m3**2*mq1*PsR + 
     &     2*Bvm2R*DRc*m1**2*m3**2*mq1*PsR - 8*Bvm1R*DRc*m2**2*m3**2*mq1*PsR - 4*BvgR*DRc*m3**3*mq1*PsR - 4*BvgL*DLc*m1*m3**2*mq3c*PsR - 
     &     2*Bvm1L*DLc*m1*m3**3*mq3c*PsR + 2*Bvm2L*DLc*m1*m3**3*mq3c*PsR + 4*BvgL*DLc*m1*m4**2*mq3c*PsR + 2*Bvm1L*DLc*m1*m3*m4**2*mq3c*PsR - 
     &     2*Bvm2L*DLc*m1*m3*m4**2*mq3c*PsR + 16*BvgR*DLc*m1*m3*mq1*mq3c*PsR + 6*Bvm1R*DLc*m1*m3**2*mq1*mq3c*PsR - 2*Bvm2R*DLc*m1*m3**2*mq1*mq3c*PsR - 
     &     2*Bvm1R*DLc*m1*m4**2*mq1*mq3c*PsR - 2*Bvm2R*DLc*m1*m4**2*mq1*mq3c*PsR + 8*BsgL*DRc*m1**2*m2**2*PvL + 12*Bsm1L*DRc*m1**2*m2**2*m3*PvL + 
     &     4*Bsm2L*DRc*m1**2*m2**2*m3*PvL - 8*BsgL*DRc*m1**2*m3**2*PvL + 4*BsgL*DRc*m2**2*m3**2*PvL + 8*Bsm1L*DRc*m2**2*m3**3*PvL - 4*BsgL*DRc*m3**4*PvL + 
     &     4*BsgL*DRc*m2**2*m4**2*PvL + 8*Bsm1L*DRc*m1**2*m3*m4**2*PvL + 8*Bsm2L*DRc*m1**2*m3*m4**2*PvL + 4*Bsm1L*DRc*m2**2*m3*m4**2*PvL + 
     &     4*Bsm2L*DRc*m2**2*m3*m4**2*PvL - 4*BsgL*DRc*m3**2*m4**2*PvL + 4*Bsm1L*DRc*m3**3*m4**2*PvL + 4*Bsm2L*DRc*m3**3*m4**2*PvL + 
     &     4*Bsm1L*DRc*m3*m4**4*PvL + 4*Bsm2L*DRc*m3*m4**4*PvL - 8*BvgL*DRc*m1**2*m2**2*mq1*PvL - 12*Bvm1L*DRc*m1**2*m2**2*m3*mq1*PvL - 
     &     4*Bvm2L*DRc*m1**2*m2**2*m3*mq1*PvL - 8*BvgL*DRc*m1**2*m3**2*mq1*PvL - 12*BvgL*DRc*m2**2*m3**2*mq1*PvL - 8*Bvm1L*DRc*m2**2*m3**3*mq1*PvL - 
     &     4*BvgL*DRc*m3**4*mq1*PvL + 4*BvgL*DRc*m2**2*m4**2*mq1*PvL - 8*Bvm1L*DRc*m1**2*m3*m4**2*mq1*PvL - 8*Bvm2L*DRc*m1**2*m3*m4**2*mq1*PvL - 
     &     4*Bvm1L*DRc*m2**2*m3*m4**2*mq1*PvL - 4*Bvm2L*DRc*m2**2*m3*m4**2*mq1*PvL - 4*BvgL*DRc*m3**2*m4**2*mq1*PvL - 4*Bvm1L*DRc*m3**3*m4**2*mq1*PvL - 
     &     4*Bvm2L*DRc*m3**3*m4**2*mq1*PvL - 4*Bvm1L*DRc*m3*m4**4*mq1*PvL - 4*Bvm2L*DRc*m3*m4**4*mq1*PvL + 4*BsgL*DLc*m1*m3**2*mq3c*PvL - 
     &     2*Bsm1L*DLc*m1*m3**3*mq3c*PvL + 2*Bsm2L*DLc*m1*m3**3*mq3c*PvL - 4*BsgL*DLc*m1*m4**2*mq3c*PvL + 2*Bsm1L*DLc*m1*m3*m4**2*mq3c*PvL - 
     &     2*Bsm2L*DLc*m1*m3*m4**2*mq3c*PvL + 4*BvgL*DLc*m1*m3**2*mq1*mq3c*PvL + 2*Bvm1L*DLc*m1*m3**3*mq1*mq3c*PvL - 2*Bvm2L*DLc*m1*m3**3*mq1*mq3c*PvL - 
     &     4*BvgL*DLc*m1*m4**2*mq1*mq3c*PvL - 2*Bvm1L*DLc*m1*m3*m4**2*mq1*mq3c*PvL + 2*Bvm2L*DLc*m1*m3*m4**2*mq1*mq3c*PvL + 4*BvgL*DRc*m1**2*PsR*sman + 
     &     4*BvgL*DRc*m2**2*PsR*sman + 4*BsgR*DRc*m3*PsR*sman + Bvm1L*DRc*m1**2*m3*PsR*sman - Bvm2L*DRc*m1**2*m3*PsR*sman + 3*Bvm1L*DRc*m2**2*m3*PsR*sman + 
     &     Bvm2L*DRc*m2**2*m3*PsR*sman - 4*BvgL*DRc*m3**2*PsR*sman - Bvm1L*DRc*m3**3*PsR*sman + Bvm2L*DRc*m3**3*PsR*sman + 4*BvgL*DRc*m4**2*PsR*sman - 
     &     Bvm1L*DRc*m3*m4**2*PsR*sman - 3*Bvm2L*DRc*m3*m4**2*PsR*sman - Bvm1R*DRc*m1**2*mq1*PsR*sman - Bvm2R*DRc*m1**2*mq1*PsR*sman - 
     &     Bvm1R*DRc*m2**2*mq1*PsR*sman - Bvm2R*DRc*m2**2*mq1*PsR*sman + 4*BvgR*DRc*m3*mq1*PsR*sman + 3*Bvm1R*DRc*m3**2*mq1*PsR*sman - 
     &     Bvm2R*DRc*m3**2*mq1*PsR*sman - Bvm1R*DRc*m4**2*mq1*PsR*sman - Bvm2R*DRc*m4**2*mq1*PsR*sman - 4*BvgL*DLc*m1*mq3c*PsR*sman - 
     &     6*Bvm1L*DLc*m1*m3*mq3c*PsR*sman - 2*Bvm2L*DLc*m1*m3*mq3c*PsR*sman + 2*Bvm1R*DLc*m1*mq1*mq3c*PsR*sman + 2*Bvm2R*DLc*m1*mq1*mq3c*PsR*sman + 
     &     6*Bvm1R*DRc*m1**2*m2**2*PvL*sman + 2*Bvm2R*DRc*m1**2*m2**2*PvL*sman + Bsm1L*DRc*m1**2*m3*PvL*sman - Bsm2L*DRc*m1**2*m3*PvL*sman + 
     &     8*BvgR*DRc*m1**2*m3*PvL*sman + 3*Bsm1L*DRc*m2**2*m3*PvL*sman + Bsm2L*DRc*m2**2*m3*PvL*sman + 4*BvgR*DRc*m2**2*m3*PvL*sman + 
     &     2*Bvm1R*DRc*m1**2*m3**2*PvL*sman - 2*Bvm2R*DRc*m1**2*m3**2*PvL*sman + 8*Bvm1R*DRc*m2**2*m3**2*PvL*sman - Bsm1L*DRc*m3**3*PvL*sman + 
     &     Bsm2L*DRc*m3**3*PvL*sman + 4*BvgR*DRc*m3**3*PvL*sman - Bsm1L*DRc*m3*m4**2*PvL*sman - 3*Bsm2L*DRc*m3*m4**2*PvL*sman - 
     &     4*BvgL*DRc*m1**2*mq1*PvL*sman - 4*BvgL*DRc*m2**2*mq1*PvL*sman - Bvm1L*DRc*m1**2*m3*mq1*PvL*sman + Bvm2L*DRc*m1**2*m3*mq1*PvL*sman - 
     &     3*Bvm1L*DRc*m2**2*m3*mq1*PvL*sman - Bvm2L*DRc*m2**2*m3*mq1*PvL*sman + 4*BvgL*DRc*m3**2*mq1*PvL*sman + Bvm1L*DRc*m3**3*mq1*PvL*sman - 
     &     Bvm2L*DRc*m3**3*mq1*PvL*sman - 4*BvgL*DRc*m4**2*mq1*PvL*sman + Bvm1L*DRc*m3*m4**2*mq1*PvL*sman + 3*Bvm2L*DRc*m3*m4**2*mq1*PvL*sman + 
     &     4*BsgL*DLc*m1*mq3c*PvL*sman - 6*Bsm1L*DLc*m1*m3*mq3c*PvL*sman - 2*Bsm2L*DLc*m1*m3*mq3c*PvL*sman - 16*BvgR*DLc*m1*m3*mq3c*PvL*sman - 
     &     6*Bvm1R*DLc*m1*m3**2*mq3c*PvL*sman + 2*Bvm2R*DLc*m1*m3**2*mq3c*PvL*sman + 2*Bvm1R*DLc*m1*m4**2*mq3c*PvL*sman + 2*Bvm2R*DLc*m1*m4**2*mq3c*PvL*sman + 
     &     4*BvgL*DLc*m1*mq1*mq3c*PvL*sman + 6*Bvm1L*DLc*m1*m3*mq1*mq3c*PvL*sman + 2*Bvm2L*DLc*m1*m3*mq1*mq3c*PvL*sman - 4*BvgL*DRc*PsR*sman**2 - 
     &     3*Bvm1L*DRc*m3*PsR*sman**2 - Bvm2L*DRc*m3*PsR*sman**2 + Bvm1R*DRc*mq1*PsR*sman**2 + Bvm2R*DRc*mq1*PsR*sman**2 + Bvm1R*DRc*m1**2*PvL*sman**2 + 
     &     Bvm2R*DRc*m1**2*PvL*sman**2 + Bvm1R*DRc*m2**2*PvL*sman**2 + Bvm2R*DRc*m2**2*PvL*sman**2 - 3*Bsm1L*DRc*m3*PvL*sman**2 - Bsm2L*DRc*m3*PvL*sman**2 - 
     &     4*BvgR*DRc*m3*PvL*sman**2 - 3*Bvm1R*DRc*m3**2*PvL*sman**2 + Bvm2R*DRc*m3**2*PvL*sman**2 + Bvm1R*DRc*m4**2*PvL*sman**2 + 
     &     Bvm2R*DRc*m4**2*PvL*sman**2 + 4*BvgL*DRc*mq1*PvL*sman**2 + 3*Bvm1L*DRc*m3*mq1*PvL*sman**2 + Bvm2L*DRc*m3*mq1*PvL*sman**2 - 
     &     2*Bvm1R*DLc*m1*mq3c*PvL*sman**2 - 2*Bvm2R*DLc*m1*mq3c*PvL*sman**2 - Bvm1R*DRc*PvL*sman**3 - Bvm2R*DRc*PvL*sman**3 - 
     &     DRc*((-4*BsgR + (5*Bvm1L + 3*Bvm2L)*m1**2 + (7*Bvm1L + Bvm2L)*m2**2)*m3 + 4*BvgL*(m1**2 + 2*(m2**2 + m3**2)))*PsR*tman - 
     &     2*Bvm1L*DRc*m3*m4**2*PsR*tman - 2*Bvm2L*DRc*m3*m4**2*PsR*tman + 8*BvgR*DRc*m3*mq1*PsR*tman + 4*BsgL*DRc*m1**2*PvL*tman + 
     &     4*BsgL*DRc*m2**2*PvL*tman + 12*BsgL*DRc*m3**2*PvL*tman + 4*BsgL*DRc*m4**2*PvL*tman - 2*Bsm1L*DRc*m3*m4**2*PvL*tman - 
     &     2*Bsm2L*DRc*m3*m4**2*PvL*tman + 2*Bvm1L*DRc*m3*m4**2*mq1*PvL*tman + 2*Bvm2L*DRc*m3*m4**2*mq1*PvL*tman - 8*BvgR*DRc*m3*PvL*sman*tman - 
     &     4*BsgL*DRc*PvL*tman**2 + Bsm1R*PsR*(-2*DLc*m1*mq3c*(3*m3**2 - m4**2 + sman) + 
     &        DRc*((-3*m3**2 + m4**2 - sman)*sman + m2**2*(8*m3**2 + sman - 5*tman) + m1**2*(6*m2**2 + 2*m3**2 + sman - tman) - m3**2*tman)) + 
     &     Bsm2R*PsR*(2*DLc*m1*mq3c*(m3**2 + m4**2 - sman) + DRc*
     &         ((m2**2 + m3**2 + m4**2 - sman)*sman + m1**2*(2*m2**2 - 2*m3**2 + sman - tman) - (m2**2 + m3**2)*tman)) - 4*BsgL*DRc*m1**2*PvL*uman - 
     &     4*BsgL*DRc*m2**2*PvL*uman + 4*BsgL*DRc*m3**2*PvL*uman - 4*BsgL*DRc*m4**2*PvL*uman + 4*BsgL*DRc*PvL*uman**2 - 
     &     DRc*(Bsm1R*m4**2*PsR*tman + Bsm2R*m4**2*PsR*tman + 4*BvgL*m4**2*PsR*tman - Bvm1R*m1**2*mq1*PsR*tman - Bvm2R*m1**2*mq1*PsR*tman - 
     &        5*Bvm1R*m2**2*mq1*PsR*tman - Bvm2R*m2**2*mq1*PsR*tman + 4*BvgR*m3*mq1*PsR*tman - Bvm1R*m3**2*mq1*PsR*tman - Bvm2R*m3**2*mq1*PsR*tman - 
     &        Bvm1R*m4**2*mq1*PsR*tman - Bvm2R*m4**2*mq1*PsR*tman + 4*BsgL*m1**2*PvL*tman + 8*BsgL*m2**2*PvL*tman + 5*Bsm1L*m1**2*m3*PvL*tman + 
     &        3*Bsm2L*m1**2*m3*PvL*tman + 7*Bsm1L*m2**2*m3*PvL*tman + Bsm2L*m2**2*m3*PvL*tman + 8*BsgL*m3**2*PvL*tman + 3*Bsm1L*m3**3*PvL*tman + 
     &        Bsm2L*m3**3*PvL*tman + 4*BsgL*m4**2*PvL*tman + 5*Bsm1L*m3*m4**2*PvL*tman + 3*Bsm2L*m3*m4**2*PvL*tman - 4*BvgL*m1**2*mq1*PvL*tman - 
     &        8*BvgL*m2**2*mq1*PvL*tman - 8*BvgL*m3**2*mq1*PvL*tman - 4*BvgL*m4**2*mq1*PvL*tman + Bvm1R*m1**2*PvL*sman*tman + Bvm2R*m1**2*PvL*sman*tman + 
     &        5*Bvm1R*m2**2*PvL*sman*tman + Bvm2R*m2**2*PvL*sman*tman - 4*BvgR*m3*PvL*sman*tman + Bvm1R*m3**2*PvL*sman*tman + Bvm2R*m3**2*PvL*sman*tman + 
     &        Bvm1R*m4**2*PvL*sman*tman + Bvm2R*m4**2*PvL*sman*tman - Bsm1R*PsR*tman**2 - Bsm2R*PsR*tman**2 - 4*BvgL*PsR*tman**2 + Bvm1R*mq1*PsR*tman**2 + 
     &        Bvm2R*mq1*PsR*tman**2 - 4*BsgL*PvL*tman**2 - 3*Bsm1L*m3*PvL*tman**2 - Bsm2L*m3*PvL*tman**2 + 4*BvgL*mq1*PvL*tman**2 - Bvm1R*PvL*sman*tman**2 - 
     &        Bvm2R*PvL*sman*tman**2 + (4*BvgL*m1**2*PsR + 4*BvgL*m3**2*PsR + Bsm2R*(m1**2 + m2**2 - m3**2 + m4**2)*PsR - 5*Bvm1R*m1**2*mq1*PsR - 
     &           Bvm2R*m1**2*mq1*PsR - Bvm1R*m2**2*mq1*PsR - Bvm2R*m2**2*mq1*PsR - 3*Bvm1R*m3**2*mq1*PsR + Bvm2R*m3**2*mq1*PsR - Bvm1R*m4**2*mq1*PsR - 
     &           Bvm2R*m4**2*mq1*PsR + 4*BsgL*m1**2*PvL + 7*Bsm1L*m1**2*m3*PvL + Bsm2L*m1**2*m3*PvL + 5*Bsm1L*m2**2*m3*PvL + 3*Bsm2L*m2**2*m3*PvL + 
     &           4*BsgL*m3**2*PvL + 3*Bsm1L*m3**3*PvL + Bsm2L*m3**3*PvL + 7*Bsm1L*m3*m4**2*PvL + 5*Bsm2L*m3*m4**2*PvL - 4*BvgL*m1**2*mq1*PvL - 
     &           4*BvgL*m3**2*mq1*PvL + 5*Bvm1R*m1**2*PvL*sman + Bvm2R*m1**2*PvL*sman + Bvm1R*m2**2*PvL*sman + Bvm2R*m2**2*PvL*sman - 2*Bsm1L*m3*PvL*sman - 
     &           2*Bsm2L*m3*PvL*sman + 3*Bvm1R*m3**2*PvL*sman - Bvm2R*m3**2*PvL*sman + Bvm1R*m4**2*PvL*sman + Bvm2R*m4**2*PvL*sman + 
     &           Bsm1R*PsR*(5*m1**2 + m2**2 + 3*m3**2 + m4**2 - 4*tman) - 
     &           2*(-2*Bvm1R*mq1*PsR + 2*BsgL*PvL + 3*Bsm1L*m3*PvL + Bsm2L*m3*PvL + 2*BvgL*(PsR - mq1*PvL) + 2*Bvm1R*PvL*sman)*tman)*uman - 
     &        ((Bsm1R + Bsm2R - (Bvm1R + Bvm2R)*mq1)*PsR + (3*Bsm1L + Bsm2L)*m3*PvL + (Bvm1R + Bvm2R)*PvL*sman)*uman**2 + 
     &        Bvm1L*m3*(-(((5*m1**2 + 7*m2**2)*mq1*PvL - (PsR - mq1*PvL)*(3*m3**2 + 5*m4**2 - 3*tman))*tman) + 
     &           (PsR - mq1*PvL)*(7*m1**2 + 5*m2**2 + 3*m3**2 + 7*m4**2 - 2*sman - 6*tman)*uman - 3*(PsR - mq1*PvL)*uman**2) + 
     &        Bvm2L*m3*(-(((3*m1**2 + m2**2)*mq1*PvL - (PsR - mq1*PvL)*(m3**2 + 3*m4**2 - tman))*tman) + 
     &           (PsR - mq1*PvL)*(m1**2 + 3*m2**2 + m3**2 + 5*m4**2 - 2*(sman + tman))*uman + (-PsR + mq1*PvL)*uman**2)) - 
     &     (-4*BvgL*DRc*m1**2*m2**2*m3**2*PsR - 4*Bvm1L*DRc*m1**2*m2**2*m3**3*PsR + 2*Bsm2R*DRc*m1**2*m2**2*m4**2*PsR + 4*Bvm2L*DRc*m1**2*m2**2*m3*m4**2*PsR - 
     &        2*Bsm2R*DRc*m1**2*m3**2*m4**2*PsR - 4*BvgL*DRc*m1**2*m3**2*m4**2*PsR - 2*BvgL*DRc*m2**2*m3**2*m4**2*PsR - 8*Bvm1L*DRc*m1**2*m3**3*m4**2*PsR - 
     &        4*Bvm1L*DRc*m2**2*m3**3*m4**2*PsR - 2*BvgL*DRc*m3**4*m4**2*PsR - 4*Bvm1L*DRc*m3**5*m4**2*PsR - 2*BvgL*DRc*m2**2*m4**4*PsR - 
     &        4*Bvm1L*DRc*m1**2*m3*m4**4*PsR + 4*Bvm2L*DRc*m1**2*m3*m4**4*PsR - 2*Bvm1L*DRc*m2**2*m3*m4**4*PsR + 2*Bvm2L*DRc*m2**2*m3*m4**4*PsR - 
     &        2*BvgL*DRc*m3**2*m4**4*PsR - 6*Bvm1L*DRc*m3**3*m4**4*PsR + 2*Bvm2L*DRc*m3**3*m4**4*PsR + 4*Bvm2L*DRc*m3*m4**6*PsR - 
     &        2*BvgR*DRc*m1**2*m2**2*m3*mq1*PsR + 2*Bvm1R*DRc*m1**2*m2**2*m3**2*mq1*PsR - 2*BvgR*DRc*m1**2*m3**3*mq1*PsR - 2*Bvm1R*DRc*m1**2*m3**4*mq1*PsR - 
     &        2*Bvm2R*DRc*m1**2*m2**2*m4**2*mq1*PsR - 4*BvgR*DRc*m1**2*m3*m4**2*mq1*PsR - 2*BvgR*DRc*m2**2*m3*m4**2*mq1*PsR + 
     &        2*Bvm2R*DRc*m1**2*m3**2*m4**2*mq1*PsR - 2*BvgR*DRc*m3**3*m4**2*mq1*PsR - 4*BvgR*DRc*m3*m4**4*mq1*PsR - 2*BvgL*DLc*m1*m3**4*mq3c*PsR - 
     &        2*Bvm1L*DLc*m1*m3**5*mq3c*PsR + 2*Bsm2R*DLc*m1*m3**2*m4**2*mq3c*PsR + 2*Bvm1L*DLc*m1*m3**3*m4**2*mq3c*PsR + 
     &        2*Bvm2L*DLc*m1*m3**3*m4**2*mq3c*PsR + 2*Bsm2R*DLc*m1*m4**4*mq3c*PsR + 2*BvgL*DLc*m1*m4**4*mq3c*PsR + 4*Bvm1L*DLc*m1*m3*m4**4*mq3c*PsR + 
     &        2*Bvm2L*DLc*m1*m3*m4**4*mq3c*PsR + 2*Bvm1R*DLc*m1*m3**4*mq1*mq3c*PsR + 2*Bvm1R*DLc*m1*m3**2*m4**2*mq1*mq3c*PsR - 
     &        2*Bvm2R*DLc*m1*m3**2*m4**2*mq1*mq3c*PsR - 2*Bvm2R*DLc*m1*m4**4*mq1*mq3c*PsR - 4*BsgL*DRc*m1**2*m2**2*m3**2*PvL - 
     &        4*Bsm1L*DRc*m1**2*m2**2*m3**3*PvL + 4*Bsm2L*DRc*m1**2*m2**2*m3*m4**2*PvL - 4*BsgL*DRc*m1**2*m3**2*m4**2*PvL - 2*BsgL*DRc*m2**2*m3**2*m4**2*PvL - 
     &        8*Bsm1L*DRc*m1**2*m3**3*m4**2*PvL - 4*Bsm1L*DRc*m2**2*m3**3*m4**2*PvL - 2*BsgL*DRc*m3**4*m4**2*PvL - 4*Bsm1L*DRc*m3**5*m4**2*PvL - 
     &        2*BsgL*DRc*m2**2*m4**4*PvL - 4*Bsm1L*DRc*m1**2*m3*m4**4*PvL + 4*Bsm2L*DRc*m1**2*m3*m4**4*PvL - 2*Bsm1L*DRc*m2**2*m3*m4**4*PvL + 
     &        2*Bsm2L*DRc*m2**2*m3*m4**4*PvL - 2*BsgL*DRc*m3**2*m4**4*PvL - 6*Bsm1L*DRc*m3**3*m4**4*PvL + 2*Bsm2L*DRc*m3**3*m4**4*PvL + 
     &        4*Bsm2L*DRc*m3*m4**6*PvL + 4*BvgL*DRc*m1**2*m2**2*m3**2*mq1*PvL + 4*Bvm1L*DRc*m1**2*m2**2*m3**3*mq1*PvL - 
     &        4*Bvm2L*DRc*m1**2*m2**2*m3*m4**2*mq1*PvL + 4*BvgL*DRc*m1**2*m3**2*m4**2*mq1*PvL + 2*BvgL*DRc*m2**2*m3**2*m4**2*mq1*PvL + 
     &        8*Bvm1L*DRc*m1**2*m3**3*m4**2*mq1*PvL + 4*Bvm1L*DRc*m2**2*m3**3*m4**2*mq1*PvL + 2*BvgL*DRc*m3**4*m4**2*mq1*PvL + 
     &        4*Bvm1L*DRc*m3**5*m4**2*mq1*PvL + 2*BvgL*DRc*m2**2*m4**4*mq1*PvL + 4*Bvm1L*DRc*m1**2*m3*m4**4*mq1*PvL - 4*Bvm2L*DRc*m1**2*m3*m4**4*mq1*PvL + 
     &        2*Bvm1L*DRc*m2**2*m3*m4**4*mq1*PvL - 2*Bvm2L*DRc*m2**2*m3*m4**4*mq1*PvL + 2*BvgL*DRc*m3**2*m4**4*mq1*PvL + 6*Bvm1L*DRc*m3**3*m4**4*mq1*PvL - 
     &        2*Bvm2L*DRc*m3**3*m4**4*mq1*PvL - 4*Bvm2L*DRc*m3*m4**6*mq1*PvL - 2*BsgL*DLc*m1*m3**4*mq3c*PvL - 2*Bsm1L*DLc*m1*m3**5*mq3c*PvL + 
     &        2*Bsm1L*DLc*m1*m3**3*m4**2*mq3c*PvL + 2*Bsm2L*DLc*m1*m3**3*m4**2*mq3c*PvL + 2*BsgL*DLc*m1*m4**4*mq3c*PvL + 4*Bsm1L*DLc*m1*m3*m4**4*mq3c*PvL + 
     &        2*Bsm2L*DLc*m1*m3*m4**4*mq3c*PvL + 2*BvgL*DLc*m1*m3**4*mq1*mq3c*PvL + 2*Bvm1L*DLc*m1*m3**5*mq1*mq3c*PvL - 
     &        2*Bvm1L*DLc*m1*m3**3*m4**2*mq1*mq3c*PvL - 2*Bvm2L*DLc*m1*m3**3*m4**2*mq1*mq3c*PvL - 2*BvgL*DLc*m1*m4**4*mq1*mq3c*PvL - 
     &        4*Bvm1L*DLc*m1*m3*m4**4*mq1*mq3c*PvL - 2*Bvm2L*DLc*m1*m3*m4**4*mq1*mq3c*PvL + 
     &        2*BsgR*DRc*m3*(m1**2 + m4**2)*PsR*(m2**2 + m3**2 + 2*m4**2 - sman) + 
     &        (Bsm2R*DRc*m1**2*m4**2*PsR + Bsm2R*DRc*m2**2*m4**2*PsR - Bvm2L*DRc*m1**2*m3*m4**2*PsR + Bvm2L*DRc*m2**2*m3*m4**2*PsR + 
     &           Bsm2R*DRc*m3**2*m4**2*PsR + Bvm2L*DRc*m3**3*m4**2*PsR + Bsm2R*DRc*m4**4*PsR - Bvm2L*DRc*m3*m4**4*PsR - 2*Bvm1R*DRc*m1**2*m2**2*mq1*PsR + 
     &           2*BvgR*DRc*m1**2*m3*mq1*PsR + 3*Bvm1R*DRc*m1**2*m3**2*mq1*PsR + Bvm1R*DRc*m2**2*m3**2*mq1*PsR + Bvm1R*DRc*m3**4*mq1*PsR - 
     &           Bvm2R*DRc*m1**2*m4**2*mq1*PsR - Bvm2R*DRc*m2**2*m4**2*mq1*PsR + 2*BvgR*DRc*m3*m4**2*mq1*PsR + Bvm1R*DRc*m3**2*m4**2*mq1*PsR - 
     &           Bvm2R*DRc*m3**2*m4**2*mq1*PsR - Bvm2R*DRc*m4**4*mq1*PsR - 2*Bsm2R*DLc*m1*m4**2*mq3c*PsR - 2*Bvm2L*DLc*m1*m3*m4**2*mq3c*PsR - 
     &           4*Bvm1R*DLc*m1*m3**2*mq1*mq3c*PsR - 2*Bvm1R*DLc*m1*m4**2*mq1*mq3c*PsR + 2*Bvm2R*DLc*m1*m4**2*mq1*mq3c*PsR + 
     &           (2*BvgR*DRc*m1**2*m2**2*m3 - 2*Bvm1R*DRc*m1**2*m2**2*m3**2 + 2*BvgR*DRc*m1**2*m3**3 + 2*Bvm1R*DRc*m1**2*m3**4 + 
     &              2*Bvm2R*DRc*m1**2*m2**2*m4**2 - Bsm2L*DRc*m1**2*m3*m4**2 + 4*BvgR*DRc*m1**2*m3*m4**2 + Bsm2L*DRc*m2**2*m3*m4**2 + 
     &              2*BvgR*DRc*m2**2*m3*m4**2 - 2*Bvm2R*DRc*m1**2*m3**2*m4**2 + Bsm2L*DRc*m3**3*m4**2 + 2*BvgR*DRc*m3**3*m4**2 - Bsm2L*DRc*m3*m4**4 + 
     &              4*BvgR*DRc*m3*m4**4 + BsgL*DRc*(-(m3**2*(m2**2 + m3**2)) + (m2**2 + 2*m3**2)*m4**2 + m4**4 + m1**2*(2*m2**2 - m3**2 + m4**2)) + 
     &              Bvm2L*DRc*m1**2*m3*m4**2*mq1 - Bvm2L*DRc*m2**2*m3*m4**2*mq1 - Bvm2L*DRc*m3**3*m4**2*mq1 + Bvm2L*DRc*m3*m4**4*mq1 + 
     &              4*BsgL*DLc*m1*m3**2*mq3c + 2*DLc*m1*(-(Bvm1R*m3**2*(m3**2 + m4**2)) + m4**2*(-(Bsm2L*m3) + Bvm2R*(m3**2 + m4**2) + Bvm2L*m3*mq1))*mq3c + 
     &              Bsm1L*m3*(DRc*(-(m3**2*(m2**2 + m3**2)) + (4*m2**2 + 7*m3**2)*m4**2 + 6*m4**4 + m1**2*(4*m2**2 + m3**2 + 8*m4**2)) + 
     &                 2*DLc*m1*(2*m3**2 - m4**2)*mq3c))*PvL + BvgL*
     &            (DRc*(-(m3**2*(m2**2 + m3**2)) + (m2**2 + 2*m3**2)*m4**2 + m4**4 + m1**2*(2*m2**2 - m3**2 + m4**2)) + 4*DLc*m1*m3**2*mq3c)*(PsR - mq1*PvL) + 
     &           Bvm1L*m3*(DRc*(-(m3**2*(m2**2 + m3**2)) + (4*m2**2 + 7*m3**2)*m4**2 + 6*m4**4 + m1**2*(4*m2**2 + m3**2 + 8*m4**2)) + 
     &              2*DLc*m1*(2*m3**2 - m4**2)*mq3c)*(PsR - mq1*PvL))*sman + 
     &        Bsm1R*PsR*(-2*DLc*m1*mq3c*(m3**4 + m3**2*(m4**2 - 2*sman) - m4**2*sman) - 
     &           DRc*(m3**2*(m2**2 + m3**2 + m4**2)*sman + m1**2*(-2*m3**4 + 2*m2**2*(m3**2 - sman) + 3*m3**2*sman))))*xv + 
     &     DRc*tman*(-2*BvgL*m1**2*m3**2*PsR - 2*Bvm1L*m1**2*m3**3*PsR - BvgL*m3**4*PsR - Bvm1L*m3**5*PsR + Bsm2R*m1**2*m4**2*PsR + 2*BvgL*m1**2*m4**2*PsR + 
     &        2*Bvm2L*m1**2*m3*m4**2*PsR + Bsm2R*m3**2*m4**2*PsR - Bvm1L*m3**3*m4**2*PsR + Bvm2L*m3**3*m4**2*PsR + BvgL*m4**4*PsR + Bvm2L*m3*m4**4*PsR + 
     &        2*BsgR*m3*(m1**2 + m4**2)*PsR - 2*BvgR*m1**2*m3*mq1*PsR + Bvm1R*m1**2*m3**2*mq1*PsR + Bvm1R*m3**4*mq1*PsR - Bvm2R*m1**2*m4**2*mq1*PsR - 
     &        2*BvgR*m3*m4**2*mq1*PsR - Bvm2R*m3**2*m4**2*mq1*PsR - 2*BsgL*m1**2*m3**2*PvL - 2*Bsm1L*m1**2*m3**3*PvL - BsgL*m3**4*PvL - Bsm1L*m3**5*PvL + 
     &        2*BsgL*m1**2*m4**2*PvL + 2*Bsm2L*m1**2*m3*m4**2*PvL - Bsm1L*m3**3*m4**2*PvL + Bsm2L*m3**3*m4**2*PvL + BsgL*m4**4*PvL + Bsm2L*m3*m4**4*PvL + 
     &        2*BvgL*m1**2*m3**2*mq1*PvL + 2*Bvm1L*m1**2*m3**3*mq1*PvL + BvgL*m3**4*mq1*PvL + Bvm1L*m3**5*mq1*PvL - 2*BvgL*m1**2*m4**2*mq1*PvL - 
     &        2*Bvm2L*m1**2*m3*m4**2*mq1*PvL + Bvm1L*m3**3*m4**2*mq1*PvL - Bvm2L*m3**3*m4**2*mq1*PvL - BvgL*m4**4*mq1*PvL - Bvm2L*m3*m4**4*mq1*PvL + 
     &        2*Bvm1L*m1**2*m3*PsR*sman + BvgL*m3**2*PsR*sman + Bvm1L*m3**3*PsR*sman - BvgL*m4**2*PsR*sman + Bvm1L*m3*m4**2*PsR*sman - 
     &        Bvm1R*m1**2*mq1*PsR*sman - Bvm1R*m3**2*mq1*PsR*sman + 2*Bsm1L*m1**2*m3*PvL*sman + 2*BvgR*m1**2*m3*PvL*sman + BsgL*m3**2*PvL*sman - 
     &        Bvm1R*m1**2*m3**2*PvL*sman + Bsm1L*m3**3*PvL*sman - Bvm1R*m3**4*PvL*sman - BsgL*m4**2*PvL*sman + Bvm2R*m1**2*m4**2*PvL*sman + 
     &        Bsm1L*m3*m4**2*PvL*sman + 2*BvgR*m3*m4**2*PvL*sman + Bvm2R*m3**2*m4**2*PvL*sman - 2*Bvm1L*m1**2*m3*mq1*PvL*sman - BvgL*m3**2*mq1*PvL*sman - 
     &        Bvm1L*m3**3*mq1*PvL*sman + BvgL*m4**2*mq1*PvL*sman - Bvm1L*m3*m4**2*mq1*PvL*sman + Bvm1R*m1**2*PvL*sman**2 + Bvm1R*m3**2*PvL*sman**2 - 
     &        Bsm1R*PsR*(m3**2 - sman)*(m1**2 + m3**2 - tman) + 
     &        (-(Bsm2R*m4**2*PsR) - Bvm2L*m3*m4**2*PsR - Bvm1R*m3**2*mq1*PsR + Bvm2R*m4**2*mq1*PsR + BsgL*m3**2*PvL + Bsm1L*m3**3*PvL - BsgL*m4**2*PvL - 
     &           Bsm2L*m3*m4**2*PvL + Bvm2L*m3*m4**2*mq1*PvL + Bvm1L*m3*(PsR - mq1*PvL)*(m3**2 - sman) + BvgL*(PsR - mq1*PvL)*(m3**2 - m4**2 - sman) + 
     &           Bvm1R*mq1*PsR*sman - BsgL*PvL*sman - Bsm1L*m3*PvL*sman + Bvm1R*m3**2*PvL*sman - Bvm2R*m4**2*PvL*sman - Bvm1R*PvL*sman**2)*tman)*xv + 
     &     (Bvm1L*DRc*m1**2*m3*PsR*sman**2 - Bvm1L*DRc*m2**2*m3*PsR*sman**2 - 2*Bvm1L*DRc*m3**3*PsR*sman**2 + Bsm2R*DRc*m4**2*PsR*sman**2 + 
     &        3*Bvm1L*DRc*m3*m4**2*PsR*sman**2 + Bvm2L*DRc*m3*m4**2*PsR*sman**2 + Bvm1R*DRc*m1**2*mq1*PsR*sman**2 + Bvm1R*DRc*m2**2*mq1*PsR*sman**2 + 
     &        2*Bvm1R*DRc*m3**2*mq1*PsR*sman**2 + Bvm1R*DRc*m4**2*mq1*PsR*sman**2 - Bvm2R*DRc*m4**2*mq1*PsR*sman**2 + 2*Bvm1L*DLc*m1*m3*mq3c*PsR*sman**2 - 
     &        2*Bvm1R*DLc*m1*mq1*mq3c*PsR*sman**2 - BsgL*DRc*m1**2*PvL*sman**2 - BsgL*DRc*m2**2*PvL*sman**2 - 2*Bvm1R*DRc*m1**2*m2**2*PvL*sman**2 + 
     &        Bsm1L*DRc*m1**2*m3*PvL*sman**2 + 2*BvgR*DRc*m1**2*m3*PvL*sman**2 - Bsm1L*DRc*m2**2*m3*PvL*sman**2 - 2*BsgL*DRc*m3**2*PvL*sman**2 + 
     &        3*Bvm1R*DRc*m1**2*m3**2*PvL*sman**2 + Bvm1R*DRc*m2**2*m3**2*PvL*sman**2 - 2*Bsm1L*DRc*m3**3*PvL*sman**2 + Bvm1R*DRc*m3**4*PvL*sman**2 - 
     &        Bvm2R*DRc*m1**2*m4**2*PvL*sman**2 - Bvm2R*DRc*m2**2*m4**2*PvL*sman**2 + 3*Bsm1L*DRc*m3*m4**2*PvL*sman**2 + Bsm2L*DRc*m3*m4**2*PvL*sman**2 + 
     &        2*BvgR*DRc*m3*m4**2*PvL*sman**2 + Bvm1R*DRc*m3**2*m4**2*PvL*sman**2 - Bvm2R*DRc*m3**2*m4**2*PvL*sman**2 - Bvm2R*DRc*m4**4*PvL*sman**2 - 
     &        Bvm1L*DRc*m1**2*m3*mq1*PvL*sman**2 + Bvm1L*DRc*m2**2*m3*mq1*PvL*sman**2 + 2*Bvm1L*DRc*m3**3*mq1*PvL*sman**2 - 
     &        3*Bvm1L*DRc*m3*m4**2*mq1*PvL*sman**2 - Bvm2L*DRc*m3*m4**2*mq1*PvL*sman**2 + 2*BsgL*DLc*m1*mq3c*PvL*sman**2 + 
     &        2*Bsm1L*DLc*m1*m3*mq3c*PvL*sman**2 - 4*Bvm1R*DLc*m1*m3**2*mq3c*PvL*sman**2 - 2*Bvm1R*DLc*m1*m4**2*mq3c*PvL*sman**2 + 
     &        2*Bvm2R*DLc*m1*m4**2*mq3c*PvL*sman**2 - 2*Bvm1L*DLc*m1*m3*mq1*mq3c*PvL*sman**2 + 
     &        BvgL*(PsR - mq1*PvL)*(2*DLc*m1*mq3c - DRc*(m1**2 + m2**2 + 2*m3**2 - sman))*sman**2 + Bvm1L*DRc*m3*PsR*sman**3 - Bvm1R*DRc*mq1*PsR*sman**3 + 
     &        BsgL*DRc*PvL*sman**3 - Bvm1R*DRc*m1**2*PvL*sman**3 - Bvm1R*DRc*m2**2*PvL*sman**3 + Bsm1L*DRc*m3*PvL*sman**3 - 2*Bvm1R*DRc*m3**2*PvL*sman**3 - 
     &        Bvm1R*DRc*m4**2*PvL*sman**3 + Bvm2R*DRc*m4**2*PvL*sman**3 - Bvm1L*DRc*m3*mq1*PvL*sman**3 + 2*Bvm1R*DLc*m1*mq3c*PvL*sman**3 + 
     &        Bvm1R*DRc*PvL*sman**4 + BvgL*DRc*(-(((m1**2 + m2**2)*m3**2 + (m1**2 - m2**2 + 2*m3**2)*m4**2)*(PsR - mq1*PvL)) + 
     &           (m1**2 + m2**2 + 2*m4**2)*PsR*sman)*tman + DRc*
     &         (Bsm2R*m2**2*m4**2*PsR + Bvm2L*m1**2*m3*m4**2*PsR + Bvm2L*m2**2*m3*m4**2*PsR + Bsm2R*m4**4*PsR + 2*Bvm2L*m3*m4**4*PsR + 
     &           Bvm1R*m2**2*m3**2*mq1*PsR - Bvm2R*m2**2*m4**2*mq1*PsR + Bvm1R*m3**2*m4**2*mq1*PsR - Bvm2R*m4**4*mq1*PsR - BsgL*m1**2*m3**2*PvL - 
     &           BsgL*m2**2*m3**2*PvL - Bsm1L*m1**2*m3**3*PvL - Bsm1L*m2**2*m3**3*PvL - BsgL*m1**2*m4**2*PvL + BsgL*m2**2*m4**2*PvL + 
     &           Bsm2L*m1**2*m3*m4**2*PvL + Bsm2L*m2**2*m3*m4**2*PvL - 2*BsgL*m3**2*m4**2*PvL - 4*Bsm1L*m3**3*m4**2*PvL - 2*Bsm1L*m3*m4**4*PvL + 
     &           2*Bsm2L*m3*m4**4*PvL - Bvm2L*m1**2*m3*m4**2*mq1*PvL - Bvm2L*m2**2*m3*m4**2*mq1*PvL - 2*Bvm2L*m3*m4**4*mq1*PvL + 
     &           (((m1**2 + m2**2)*(BsgL + Bsm1L*m3) + (2*BsgL + Bvm2R*m2**2 + 4*Bsm1L*m3)*m4**2)*PvL - Bvm1R*(m2**2 + m4**2)*(mq1*PsR + m3**2*PvL))*sman + 
     &           Bvm1L*m3*(-(((m1**2 + m2**2)*m3**2 + 4*m3**2*m4**2 + 2*m4**4)*(PsR - mq1*PvL)) + (m1**2 + m2**2 + 4*m4**2)*PsR*sman))*tman + 
     &        Bsm1R*PsR*(2*DLc*m1*mq3c*sman**2 - DRc*(m1**2 + m2**2 + 2*m3**2 + m4**2 - sman)*sman**2 + DRc*(m2**2 + m4**2)*(-m3**2 + sman)*tman))*xv + 
     &     (4*BvgR*DRc*m1**2*m3*m4**2*mq1*PsR + 2*BvgR*DRc*m2**2*m3*m4**2*mq1*PsR + 2*BvgR*DRc*m3**3*m4**2*mq1*PsR - 4*BvgR*DLc*m1*m3*m4**2*mq1*mq3c*PsR + 
     &        4*BsgL*DRc*m1**2*m3**2*m4**2*PvL + 2*BsgL*DRc*m2**2*m3**2*m4**2*PvL + 2*BsgL*DRc*m3**4*m4**2*PvL - 2*BsgL*DRc*m2**2*m4**4*PvL - 
     &        4*Bsm1L*DRc*m1**2*m3*m4**4*PvL - 4*Bsm2L*DRc*m1**2*m3*m4**4*PvL - 2*Bsm1L*DRc*m2**2*m3*m4**4*PvL - 2*Bsm2L*DRc*m2**2*m3*m4**4*PvL + 
     &        2*BsgL*DRc*m3**2*m4**4*PvL - 2*Bsm1L*DRc*m3**3*m4**4*PvL - 2*Bsm2L*DRc*m3**3*m4**4*PvL - 2*BsgL*DLc*m1*m3**2*m4**2*mq3c*PvL + 
     &        2*BsgL*DLc*m1*m4**4*mq3c*PvL + 4*Bsm1L*DLc*m1*m3*m4**4*mq3c*PvL + 4*Bsm2L*DLc*m1*m3*m4**4*mq3c*PvL - 2*BvgR*DRc*m3*m4**2*mq1*PsR*sman + 
     &        BsgL*DRc*m1**2*m4**2*PvL*sman + BsgL*DRc*m2**2*m4**2*PvL*sman - 4*BvgR*DRc*m1**2*m3*m4**2*PvL*sman - 2*BvgR*DRc*m2**2*m3*m4**2*PvL*sman - 
     &        BsgL*DRc*m3**2*m4**2*PvL*sman - 2*BvgR*DRc*m3**3*m4**2*PvL*sman + BsgL*DRc*m4**4*PvL*sman + 2*Bsm1L*DRc*m3*m4**4*PvL*sman + 
     &        2*Bsm2L*DRc*m3*m4**4*PvL*sman - 2*BsgL*DLc*m1*m4**2*mq3c*PvL*sman + 4*BvgR*DLc*m1*m3*m4**2*mq3c*PvL*sman - BsgL*DRc*m4**2*PvL*sman**2 + 
     &        2*BvgR*DRc*m3*m4**2*PvL*sman**2 - 2*BvgR*DRc*m3*m4**2*mq1*PsR*tman - BsgL*DRc*m1**2*m4**2*PvL*tman - BsgL*DRc*m2**2*m4**2*PvL*tman - 
     &        3*BsgL*DRc*m3**2*m4**2*PvL*tman - BsgL*DRc*m4**4*PvL*tman + 2*Bsm1L*DRc*m3*m4**4*PvL*tman + 2*Bsm2L*DRc*m3*m4**4*PvL*tman + 
     &        2*BvgR*DRc*m3*m4**2*PvL*sman*tman + BsgL*DRc*m4**2*PvL*tman**2 + 
     &        DRc*(m3*(2*BvgR*mq1*PsR + Bsm1L*m3**2*PvL - Bsm2L*m4**2*PvL - (Bsm1L + 2*BvgR)*PvL*sman)*tman + 
     &           BsgL*PvL*(m4**2*(m1**2 + m2**2 - m3**2 + m4**2) + (m3**2 - m4**2 + sman)*tman))*uman - BsgL*DRc*m4**2*PvL*uman**2 + 
     &        Bvm2L*m3*m4**2*(PsR - mq1*PvL)*(2*m4**2*(2*DLc*m1*mq3c + DRc*(-2*m1**2 - m2**2 - m3**2 + sman + tman)) - DRc*tman*uman) + 
     &        Bvm1L*m3*(PsR - mq1*PvL)*(2*m4**4*(2*DLc*m1*mq3c + DRc*(-2*m1**2 - m2**2 - m3**2 + sman + tman)) + DRc*(m3**2 - sman)*tman*uman))*xv + 
     &     DRc*(PvL*sman*(-(Bvm1L*m3*(m1**2 + m2**2 + 4*m4**2)*mq1) + Bvm1R*(m2**2 + m4**2)*sman)*tman + 
     &        (-(Bvm1L*m1**2*m3**3*PsR) - 3*Bvm1L*m2**2*m3**3*PsR - Bvm1L*m3**5*PsR + Bsm2R*m1**2*m4**2*PsR + Bsm2R*m2**2*m4**2*PsR + 
     &           Bvm2L*m1**2*m3*m4**2*PsR + 3*Bvm2L*m2**2*m3*m4**2*PsR - Bsm2R*m3**2*m4**2*PsR - 5*Bvm1L*m3**3*m4**2*PsR + Bvm2L*m3**3*m4**2*PsR + 
     &           Bsm2R*m4**4*PsR + 5*Bvm2L*m3*m4**4*PsR - 2*BvgR*m2**2*m3*mq1*PsR + Bvm1R*m1**2*m3**2*mq1*PsR + Bvm1R*m2**2*m3**2*mq1*PsR - 
     &           2*BvgR*m3**3*mq1*PsR - Bvm1R*m3**4*mq1*PsR - 4*BvgR*m3*m4**2*mq1*PsR + Bvm1R*m3**2*m4**2*mq1*PsR - BsgL*m1**2*m3**2*PvL - 
     &           3*BsgL*m2**2*m3**2*PvL - Bsm1L*m1**2*m3**3*PvL - 3*Bsm1L*m2**2*m3**3*PvL - BsgL*m3**4*PvL - Bsm1L*m3**5*PvL - BsgL*m1**2*m4**2*PvL - 
     &           BsgL*m2**2*m4**2*PvL + Bsm2L*m1**2*m3*m4**2*PvL + 3*Bsm2L*m2**2*m3*m4**2*PvL - 4*BsgL*m3**2*m4**2*PvL - 5*Bsm1L*m3**3*m4**2*PvL + 
     &           Bsm2L*m3**3*m4**2*PvL - BsgL*m4**4*PvL + 5*Bsm2L*m3*m4**4*PvL + Bvm1L*m1**2*m3**3*mq1*PvL + 3*Bvm1L*m2**2*m3**3*mq1*PvL + 
     &           Bvm1L*m3**5*mq1*PvL - Bvm2L*m1**2*m3*m4**2*mq1*PvL - 3*Bvm2L*m2**2*m3*m4**2*mq1*PvL + 5*Bvm1L*m3**3*m4**2*mq1*PvL - 
     &           Bvm2L*m3**3*m4**2*mq1*PvL - 5*Bvm2L*m3*m4**4*mq1*PvL - Bsm1R*(m1**2 + m2**2 - m3**2 + m4**2)*PsR*(m3**2 - sman) + 
     &           2*BsgR*m3*PsR*(m2**2 + m3**2 + 2*m4**2 - sman) + Bvm1L*m1**2*m3*PsR*sman + 3*Bvm1L*m2**2*m3*PsR*sman + 3*Bvm1L*m3**3*PsR*sman + 
     &           5*Bvm1L*m3*m4**2*PsR*sman - 2*Bvm2L*m3*m4**2*PsR*sman - Bvm1R*m1**2*mq1*PsR*sman - Bvm1R*m2**2*mq1*PsR*sman + 2*BvgR*m3*mq1*PsR*sman + 
     &           Bvm1R*m3**2*mq1*PsR*sman - Bvm1R*m4**2*mq1*PsR*sman + BsgL*m1**2*PvL*sman + BsgL*m2**2*PvL*sman + Bsm1L*m1**2*m3*PvL*sman + 
     &           3*Bsm1L*m2**2*m3*PvL*sman + 2*BvgR*m2**2*m3*PvL*sman + BsgL*m3**2*PvL*sman - Bvm1R*m1**2*m3**2*PvL*sman - Bvm1R*m2**2*m3**2*PvL*sman + 
     &           3*Bsm1L*m3**3*PvL*sman + 2*BvgR*m3**3*PvL*sman + Bvm1R*m3**4*PvL*sman + BsgL*m4**2*PvL*sman + 5*Bsm1L*m3*m4**2*PvL*sman - 
     &           2*Bsm2L*m3*m4**2*PvL*sman + 4*BvgR*m3*m4**2*PvL*sman - Bvm1R*m3**2*m4**2*PvL*sman - Bvm1L*m1**2*m3*mq1*PvL*sman - 
     &           3*Bvm1L*m2**2*m3*mq1*PvL*sman - 3*Bvm1L*m3**3*mq1*PvL*sman - 5*Bvm1L*m3*m4**2*mq1*PvL*sman + 2*Bvm2L*m3*m4**2*mq1*PvL*sman - 
     &           2*Bvm1L*m3*PsR*sman**2 + Bvm1R*m1**2*PvL*sman**2 + Bvm1R*m2**2*PvL*sman**2 - 2*Bsm1L*m3*PvL*sman**2 - 2*BvgR*m3*PvL*sman**2 - 
     &           Bvm1R*m3**2*PvL*sman**2 + Bvm1R*m4**2*PvL*sman**2 + 2*Bvm1L*m3*mq1*PvL*sman**2 + 
     &           (Bvm2L*m3*m4**2*(-PsR + mq1*PvL) + Bvm1L*m3*(PsR - mq1*PvL)*(m3**2 - sman) + 
     &              PvL*(BsgL*m3**2 + Bsm1L*m3**3 + BsgL*m4**2 - Bsm2L*m3*m4**2 - (BsgL + Bsm1L*m3)*sman))*tman)*uman + 
     &        (-(Bsm2R*m4**2*PsR) - Bvm2L*m3*m4**2*PsR - Bvm1R*m3**2*mq1*PsR + BsgL*m3**2*PvL + Bsm1L*m3**3*PvL + BsgL*m4**2*PvL - Bsm2L*m3*m4**2*PvL + 
     &           Bvm2L*m3*m4**2*mq1*PvL + Bsm1R*PsR*(m3**2 - sman) + Bvm1L*m3*(PsR - mq1*PvL)*(m3**2 - sman) + 
     &           (Bvm1R*mq1*PsR - (BsgL + m3*(Bsm1L - Bvm1R*m3))*PvL)*sman - Bvm1R*PvL*sman**2)*uman**2 + 
     &        BvgL*(-((m1**2 + m2**2 + 2*m4**2)*mq1*PvL*sman*tman) - 
     &           (PsR - mq1*PvL)*(m3**4 + 4*m3**2*m4**2 + m4**4 + m1**2*(m3**2 + m4**2 - sman) + m2**2*(3*m3**2 + m4**2 - sman) - m3**2*sman - m4**2*sman - 
     &              (m3**2 + m4**2 - sman)*tman)*uman + (PsR - mq1*PvL)*(m3**2 + m4**2 - sman)*uman**2) + 
     &        Bvm2R*m4**2*((mq1*PsR - PvL*sman)*uman*(-m1**2 - m2**2 + m3**2 + uman) + m4**2*(-(mq1*PsR*uman) + PvL*sman*(tman + uman))))*xv + 
     &     2*BsgR*m3*PsR*(4*DRc*(2*m1**2 + m2**2 + m3**2 - sman - tman) - DRc*(m4**2*(2*m1**2 + m2**2 + m3**2 - sman - tman) + tman*uman)*xv + 
     &        2*DLc*m1*mq3c*(-4 + m4**2*xv)) + BvgL*(PsR - mq1*PvL)*
     &      (-4*DRc*((m3**2 + m4**2 - sman - tman)*(2*m3**2 + sman - tman) + (-m3**2 + m4**2)*uman - uman**2 + m1**2*(4*m3**2 + sman - tman + uman) + 
     &           m2**2*(2*m3**2 - 2*m4**2 + sman - tman + uman)) + 
     &        DRc*(m4**2*(m3**2 + m4**2 - sman - tman)*(2*m3**2 + sman - tman) + (m4**4 - m4**2*tman + sman*tman + m3**2*(-m4**2 + tman))*uman - 
     &           m4**2*uman**2 + m1**2*m4**2*(4*m3**2 + sman - tman + uman) + m2**2*m4**2*(2*m3**2 - 2*m4**2 + sman - tman + uman))*xv - 
     &        2*DLc*m1*mq3c*(m3**2 - m4**2 + sman)*(-4 + m4**2*xv))) + 
     &  A1sL*A3Lc*(8*BvgR*DLc*m1**2*m2**2*PsL - 8*BsgL*DLc*m1**2*m3*PsL - 4*BsgL*DLc*m2**2*m3*PsL + 12*Bvm1R*DLc*m1**2*m2**2*m3*PsL + 
     &     4*Bvm2R*DLc*m1**2*m2**2*m3*PsL + 8*BvgR*DLc*m1**2*m3**2*PsL + 12*BvgR*DLc*m2**2*m3**2*PsL - 4*BsgL*DLc*m3**3*PsL + 8*Bvm1R*DLc*m2**2*m3**3*PsL + 
     &     4*BvgR*DLc*m3**4*PsL - 4*BvgR*DLc*m2**2*m4**2*PsL + 8*Bvm1R*DLc*m1**2*m3*m4**2*PsL + 8*Bvm2R*DLc*m1**2*m3*m4**2*PsL + 
     &     4*Bvm1R*DLc*m2**2*m3*m4**2*PsL + 4*Bvm2R*DLc*m2**2*m3*m4**2*PsL + 4*BvgR*DLc*m3**2*m4**2*PsL + 4*Bvm1R*DLc*m3**3*m4**2*PsL + 
     &     4*Bvm2R*DLc*m3**3*m4**2*PsL + 4*Bvm1R*DLc*m3*m4**4*PsL + 4*Bvm2R*DLc*m3*m4**4*PsL - 6*Bvm1L*DLc*m1**2*m2**2*mq1*PsL - 
     &     2*Bvm2L*DLc*m1**2*m2**2*mq1*PsL - 8*BvgL*DLc*m1**2*m3*mq1*PsL - 4*BvgL*DLc*m2**2*m3*mq1*PsL - 2*Bvm1L*DLc*m1**2*m3**2*mq1*PsL + 
     &     2*Bvm2L*DLc*m1**2*m3**2*mq1*PsL - 8*Bvm1L*DLc*m2**2*m3**2*mq1*PsL - 4*BvgL*DLc*m3**3*mq1*PsL - 4*BvgR*DRc*m1*m3**2*mq3c*PsL - 
     &     2*Bvm1R*DRc*m1*m3**3*mq3c*PsL + 2*Bvm2R*DRc*m1*m3**3*mq3c*PsL + 4*BvgR*DRc*m1*m4**2*mq3c*PsL + 2*Bvm1R*DRc*m1*m3*m4**2*mq3c*PsL - 
     &     2*Bvm2R*DRc*m1*m3*m4**2*mq3c*PsL + 16*BvgL*DRc*m1*m3*mq1*mq3c*PsL + 6*Bvm1L*DRc*m1*m3**2*mq1*mq3c*PsL - 2*Bvm2L*DRc*m1*m3**2*mq1*mq3c*PsL - 
     &     2*Bvm1L*DRc*m1*m4**2*mq1*mq3c*PsL - 2*Bvm2L*DRc*m1*m4**2*mq1*mq3c*PsL + 8*BsgR*DLc*m1**2*m2**2*PvR + 12*Bsm1R*DLc*m1**2*m2**2*m3*PvR + 
     &     4*Bsm2R*DLc*m1**2*m2**2*m3*PvR - 8*BsgR*DLc*m1**2*m3**2*PvR + 4*BsgR*DLc*m2**2*m3**2*PvR + 8*Bsm1R*DLc*m2**2*m3**3*PvR - 4*BsgR*DLc*m3**4*PvR + 
     &     4*BsgR*DLc*m2**2*m4**2*PvR + 8*Bsm1R*DLc*m1**2*m3*m4**2*PvR + 8*Bsm2R*DLc*m1**2*m3*m4**2*PvR + 4*Bsm1R*DLc*m2**2*m3*m4**2*PvR + 
     &     4*Bsm2R*DLc*m2**2*m3*m4**2*PvR - 4*BsgR*DLc*m3**2*m4**2*PvR + 4*Bsm1R*DLc*m3**3*m4**2*PvR + 4*Bsm2R*DLc*m3**3*m4**2*PvR + 
     &     4*Bsm1R*DLc*m3*m4**4*PvR + 4*Bsm2R*DLc*m3*m4**4*PvR - 8*BvgR*DLc*m1**2*m2**2*mq1*PvR - 12*Bvm1R*DLc*m1**2*m2**2*m3*mq1*PvR - 
     &     4*Bvm2R*DLc*m1**2*m2**2*m3*mq1*PvR - 8*BvgR*DLc*m1**2*m3**2*mq1*PvR - 12*BvgR*DLc*m2**2*m3**2*mq1*PvR - 8*Bvm1R*DLc*m2**2*m3**3*mq1*PvR - 
     &     4*BvgR*DLc*m3**4*mq1*PvR + 4*BvgR*DLc*m2**2*m4**2*mq1*PvR - 8*Bvm1R*DLc*m1**2*m3*m4**2*mq1*PvR - 8*Bvm2R*DLc*m1**2*m3*m4**2*mq1*PvR - 
     &     4*Bvm1R*DLc*m2**2*m3*m4**2*mq1*PvR - 4*Bvm2R*DLc*m2**2*m3*m4**2*mq1*PvR - 4*BvgR*DLc*m3**2*m4**2*mq1*PvR - 4*Bvm1R*DLc*m3**3*m4**2*mq1*PvR - 
     &     4*Bvm2R*DLc*m3**3*m4**2*mq1*PvR - 4*Bvm1R*DLc*m3*m4**4*mq1*PvR - 4*Bvm2R*DLc*m3*m4**4*mq1*PvR + 4*BsgR*DRc*m1*m3**2*mq3c*PvR - 
     &     2*Bsm1R*DRc*m1*m3**3*mq3c*PvR + 2*Bsm2R*DRc*m1*m3**3*mq3c*PvR - 4*BsgR*DRc*m1*m4**2*mq3c*PvR + 2*Bsm1R*DRc*m1*m3*m4**2*mq3c*PvR - 
     &     2*Bsm2R*DRc*m1*m3*m4**2*mq3c*PvR + 4*BvgR*DRc*m1*m3**2*mq1*mq3c*PvR + 2*Bvm1R*DRc*m1*m3**3*mq1*mq3c*PvR - 2*Bvm2R*DRc*m1*m3**3*mq1*mq3c*PvR - 
     &     4*BvgR*DRc*m1*m4**2*mq1*mq3c*PvR - 2*Bvm1R*DRc*m1*m3*m4**2*mq1*mq3c*PvR + 2*Bvm2R*DRc*m1*m3*m4**2*mq1*mq3c*PvR + 4*BvgR*DLc*m1**2*PsL*sman + 
     &     4*BvgR*DLc*m2**2*PsL*sman + 4*BsgL*DLc*m3*PsL*sman + Bvm1R*DLc*m1**2*m3*PsL*sman - Bvm2R*DLc*m1**2*m3*PsL*sman + 3*Bvm1R*DLc*m2**2*m3*PsL*sman + 
     &     Bvm2R*DLc*m2**2*m3*PsL*sman - 4*BvgR*DLc*m3**2*PsL*sman - Bvm1R*DLc*m3**3*PsL*sman + Bvm2R*DLc*m3**3*PsL*sman + 4*BvgR*DLc*m4**2*PsL*sman - 
     &     Bvm1R*DLc*m3*m4**2*PsL*sman - 3*Bvm2R*DLc*m3*m4**2*PsL*sman - Bvm1L*DLc*m1**2*mq1*PsL*sman - Bvm2L*DLc*m1**2*mq1*PsL*sman - 
     &     Bvm1L*DLc*m2**2*mq1*PsL*sman - Bvm2L*DLc*m2**2*mq1*PsL*sman + 4*BvgL*DLc*m3*mq1*PsL*sman + 3*Bvm1L*DLc*m3**2*mq1*PsL*sman - 
     &     Bvm2L*DLc*m3**2*mq1*PsL*sman - Bvm1L*DLc*m4**2*mq1*PsL*sman - Bvm2L*DLc*m4**2*mq1*PsL*sman - 4*BvgR*DRc*m1*mq3c*PsL*sman - 
     &     6*Bvm1R*DRc*m1*m3*mq3c*PsL*sman - 2*Bvm2R*DRc*m1*m3*mq3c*PsL*sman + 2*Bvm1L*DRc*m1*mq1*mq3c*PsL*sman + 2*Bvm2L*DRc*m1*mq1*mq3c*PsL*sman + 
     &     6*Bvm1L*DLc*m1**2*m2**2*PvR*sman + 2*Bvm2L*DLc*m1**2*m2**2*PvR*sman + Bsm1R*DLc*m1**2*m3*PvR*sman - Bsm2R*DLc*m1**2*m3*PvR*sman + 
     &     8*BvgL*DLc*m1**2*m3*PvR*sman + 3*Bsm1R*DLc*m2**2*m3*PvR*sman + Bsm2R*DLc*m2**2*m3*PvR*sman + 4*BvgL*DLc*m2**2*m3*PvR*sman + 
     &     2*Bvm1L*DLc*m1**2*m3**2*PvR*sman - 2*Bvm2L*DLc*m1**2*m3**2*PvR*sman + 8*Bvm1L*DLc*m2**2*m3**2*PvR*sman - Bsm1R*DLc*m3**3*PvR*sman + 
     &     Bsm2R*DLc*m3**3*PvR*sman + 4*BvgL*DLc*m3**3*PvR*sman - Bsm1R*DLc*m3*m4**2*PvR*sman - 3*Bsm2R*DLc*m3*m4**2*PvR*sman - 
     &     4*BvgR*DLc*m1**2*mq1*PvR*sman - 4*BvgR*DLc*m2**2*mq1*PvR*sman - Bvm1R*DLc*m1**2*m3*mq1*PvR*sman + Bvm2R*DLc*m1**2*m3*mq1*PvR*sman - 
     &     3*Bvm1R*DLc*m2**2*m3*mq1*PvR*sman - Bvm2R*DLc*m2**2*m3*mq1*PvR*sman + 4*BvgR*DLc*m3**2*mq1*PvR*sman + Bvm1R*DLc*m3**3*mq1*PvR*sman - 
     &     Bvm2R*DLc*m3**3*mq1*PvR*sman - 4*BvgR*DLc*m4**2*mq1*PvR*sman + Bvm1R*DLc*m3*m4**2*mq1*PvR*sman + 3*Bvm2R*DLc*m3*m4**2*mq1*PvR*sman + 
     &     4*BsgR*DRc*m1*mq3c*PvR*sman - 6*Bsm1R*DRc*m1*m3*mq3c*PvR*sman - 2*Bsm2R*DRc*m1*m3*mq3c*PvR*sman - 16*BvgL*DRc*m1*m3*mq3c*PvR*sman - 
     &     6*Bvm1L*DRc*m1*m3**2*mq3c*PvR*sman + 2*Bvm2L*DRc*m1*m3**2*mq3c*PvR*sman + 2*Bvm1L*DRc*m1*m4**2*mq3c*PvR*sman + 2*Bvm2L*DRc*m1*m4**2*mq3c*PvR*sman + 
     &     4*BvgR*DRc*m1*mq1*mq3c*PvR*sman + 6*Bvm1R*DRc*m1*m3*mq1*mq3c*PvR*sman + 2*Bvm2R*DRc*m1*m3*mq1*mq3c*PvR*sman - 4*BvgR*DLc*PsL*sman**2 - 
     &     3*Bvm1R*DLc*m3*PsL*sman**2 - Bvm2R*DLc*m3*PsL*sman**2 + Bvm1L*DLc*mq1*PsL*sman**2 + Bvm2L*DLc*mq1*PsL*sman**2 + Bvm1L*DLc*m1**2*PvR*sman**2 + 
     &     Bvm2L*DLc*m1**2*PvR*sman**2 + Bvm1L*DLc*m2**2*PvR*sman**2 + Bvm2L*DLc*m2**2*PvR*sman**2 - 3*Bsm1R*DLc*m3*PvR*sman**2 - Bsm2R*DLc*m3*PvR*sman**2 - 
     &     4*BvgL*DLc*m3*PvR*sman**2 - 3*Bvm1L*DLc*m3**2*PvR*sman**2 + Bvm2L*DLc*m3**2*PvR*sman**2 + Bvm1L*DLc*m4**2*PvR*sman**2 + 
     &     Bvm2L*DLc*m4**2*PvR*sman**2 + 4*BvgR*DLc*mq1*PvR*sman**2 + 3*Bvm1R*DLc*m3*mq1*PvR*sman**2 + Bvm2R*DLc*m3*mq1*PvR*sman**2 - 
     &     2*Bvm1L*DRc*m1*mq3c*PvR*sman**2 - 2*Bvm2L*DRc*m1*mq3c*PvR*sman**2 - Bvm1L*DLc*PvR*sman**3 - Bvm2L*DLc*PvR*sman**3 - 2*Bvm1R*DLc*m3*m4**2*PsL*tman - 
     &     2*Bvm2R*DLc*m3*m4**2*PsL*tman + 8*BvgL*DLc*m3*mq1*PsL*tman + 
     &     DLc*(-4*BvgR*(m1**2 + 2*(m2**2 + m3**2) + m4**2) - m3*
     &         (-4*BsgL + Bvm2R*(3*m1**2 + m2**2 + m3**2 + 3*m4**2) + Bvm1R*(5*m1**2 + 7*m2**2 + 3*m3**2 + 5*m4**2)) + 
     &        (-4*BvgL*m3 + Bvm2L*(m1**2 + m2**2 + m3**2 + m4**2) + Bvm1L*(m1**2 + 5*m2**2 + m3**2 + m4**2))*mq1)*PsL*tman + 4*BsgR*DLc*m1**2*PvR*tman + 
     &     4*BsgR*DLc*m2**2*PvR*tman + 12*BsgR*DLc*m3**2*PvR*tman + 4*BsgR*DLc*m4**2*PvR*tman - 2*Bsm1R*DLc*m3*m4**2*PvR*tman - 
     &     2*Bsm2R*DLc*m3*m4**2*PvR*tman + 2*Bvm1R*DLc*m3*m4**2*mq1*PvR*tman + 2*Bvm2R*DLc*m3*m4**2*mq1*PvR*tman - 8*BvgL*DLc*m3*PvR*sman*tman - 
     &     4*BsgR*DLc*PvR*tman**2 + Bsm1L*PsL*(-2*DRc*m1*mq3c*(3*m3**2 - m4**2 + sman) + 
     &        DLc*((-3*m3**2 + m4**2 - sman)*sman + m2**2*(8*m3**2 + sman - 5*tman) + m1**2*(6*m2**2 + 2*m3**2 + sman - tman) - (m3**2 + m4**2)*tman)) + 
     &     Bsm2L*PsL*(2*DRc*m1*mq3c*(m3**2 + m4**2 - sman) + DLc*
     &         ((m2**2 + m3**2 + m4**2 - sman)*sman + m1**2*(2*m2**2 - 2*m3**2 + sman - tman) - (m2**2 + m3**2 + m4**2)*tman)) - 4*BsgR*DLc*m1**2*PvR*uman - 
     &     4*BsgR*DLc*m2**2*PvR*uman + 4*BsgR*DLc*m3**2*PvR*uman - 4*BsgR*DLc*m4**2*PvR*uman + 4*BsgR*DLc*PvR*uman**2 + 
     &     (-2*Bsm2L*DLc*m1**2*m2**2*m4**2*PsL - 4*BsgL*DLc*m1**2*m3*m4**2*PsL - 2*BsgL*DLc*m2**2*m3*m4**2*PsL - 4*Bvm2R*DLc*m1**2*m2**2*m3*m4**2*PsL + 
     &        2*Bsm2L*DLc*m1**2*m3**2*m4**2*PsL + 4*BvgR*DLc*m1**2*m3**2*m4**2*PsL + 2*BvgR*DLc*m2**2*m3**2*m4**2*PsL - 2*BsgL*DLc*m3**3*m4**2*PsL + 
     &        2*BvgR*DLc*m3**4*m4**2*PsL + 2*BvgR*DLc*m2**2*m4**4*PsL - 4*BsgL*DLc*m3*m4**4*PsL - 4*Bvm2R*DLc*m1**2*m3*m4**4*PsL - 
     &        2*Bvm2R*DLc*m2**2*m3*m4**4*PsL + 2*BvgR*DLc*m3**2*m4**4*PsL - 2*Bvm2R*DLc*m3**3*m4**4*PsL - 4*Bvm2R*DLc*m3*m4**6*PsL + 
     &        2*BvgL*DLc*m1**2*m2**2*m3*mq1*PsL - 2*Bvm1L*DLc*m1**2*m2**2*m3**2*mq1*PsL + 2*BvgL*DLc*m1**2*m3**3*mq1*PsL + 2*Bvm1L*DLc*m1**2*m3**4*mq1*PsL + 
     &        2*Bvm2L*DLc*m1**2*m2**2*m4**2*mq1*PsL + 4*BvgL*DLc*m1**2*m3*m4**2*mq1*PsL + 2*BvgL*DLc*m2**2*m3*m4**2*mq1*PsL - 
     &        2*Bvm2L*DLc*m1**2*m3**2*m4**2*mq1*PsL + 2*BvgL*DLc*m3**3*m4**2*mq1*PsL + 4*BvgL*DLc*m3*m4**4*mq1*PsL + 2*BvgR*DRc*m1*m3**4*mq3c*PsL - 
     &        2*Bsm2L*DRc*m1*m3**2*m4**2*mq3c*PsL - 2*Bvm2R*DRc*m1*m3**3*m4**2*mq3c*PsL - 2*Bsm2L*DRc*m1*m4**4*mq3c*PsL - 2*BvgR*DRc*m1*m4**4*mq3c*PsL - 
     &        2*Bvm2R*DRc*m1*m3*m4**4*mq3c*PsL - 2*Bvm1L*DRc*m1*m3**4*mq1*mq3c*PsL - 2*Bvm1L*DRc*m1*m3**2*m4**2*mq1*mq3c*PsL + 
     &        2*Bvm2L*DRc*m1*m3**2*m4**2*mq1*mq3c*PsL + 2*Bvm2L*DRc*m1*m4**4*mq1*mq3c*PsL + 4*BsgR*DLc*m1**2*m2**2*m3**2*PvR + 
     &        4*Bsm1R*DLc*m1**2*m2**2*m3**3*PvR - 4*Bsm2R*DLc*m1**2*m2**2*m3*m4**2*PvR + 4*BsgR*DLc*m1**2*m3**2*m4**2*PvR + 2*BsgR*DLc*m2**2*m3**2*m4**2*PvR + 
     &        8*Bsm1R*DLc*m1**2*m3**3*m4**2*PvR + 4*Bsm1R*DLc*m2**2*m3**3*m4**2*PvR + 2*BsgR*DLc*m3**4*m4**2*PvR + 4*Bsm1R*DLc*m3**5*m4**2*PvR + 
     &        2*BsgR*DLc*m2**2*m4**4*PvR + 4*Bsm1R*DLc*m1**2*m3*m4**4*PvR - 4*Bsm2R*DLc*m1**2*m3*m4**4*PvR + 2*Bsm1R*DLc*m2**2*m3*m4**4*PvR - 
     &        2*Bsm2R*DLc*m2**2*m3*m4**4*PvR + 2*BsgR*DLc*m3**2*m4**4*PvR + 6*Bsm1R*DLc*m3**3*m4**4*PvR - 2*Bsm2R*DLc*m3**3*m4**4*PvR - 
     &        4*Bsm2R*DLc*m3*m4**6*PvR - 4*BvgR*DLc*m1**2*m2**2*m3**2*mq1*PvR + 4*Bvm2R*DLc*m1**2*m2**2*m3*m4**2*mq1*PvR - 
     &        4*BvgR*DLc*m1**2*m3**2*m4**2*mq1*PvR - 2*BvgR*DLc*m2**2*m3**2*m4**2*mq1*PvR - 2*BvgR*DLc*m3**4*m4**2*mq1*PvR - 2*BvgR*DLc*m2**2*m4**4*mq1*PvR + 
     &        4*Bvm2R*DLc*m1**2*m3*m4**4*mq1*PvR + 2*Bvm2R*DLc*m2**2*m3*m4**4*mq1*PvR - 2*BvgR*DLc*m3**2*m4**4*mq1*PvR + 2*Bvm2R*DLc*m3**3*m4**4*mq1*PvR + 
     &        4*Bvm2R*DLc*m3*m4**6*mq1*PvR + 2*BsgR*DRc*m1*m3**4*mq3c*PvR + 2*Bsm1R*DRc*m1*m3**5*mq3c*PvR - 2*Bsm1R*DRc*m1*m3**3*m4**2*mq3c*PvR - 
     &        2*Bsm2R*DRc*m1*m3**3*m4**2*mq3c*PvR - 2*BsgR*DRc*m1*m4**4*mq3c*PvR - 4*Bsm1R*DRc*m1*m3*m4**4*mq3c*PvR - 2*Bsm2R*DRc*m1*m3*m4**4*mq3c*PvR - 
     &        2*BvgR*DRc*m1*m3**4*mq1*mq3c*PvR + 2*Bvm2R*DRc*m1*m3**3*m4**2*mq1*mq3c*PvR + 2*BvgR*DRc*m1*m4**4*mq1*mq3c*PvR + 
     &        2*Bvm2R*DRc*m1*m3*m4**4*mq1*mq3c*PvR - (-2*BsgL*DLc*m1**2*m3*PsL + Bsm2L*DLc*m1**2*m4**2*PsL + Bsm2L*DLc*m2**2*m4**2*PsL - 
     &           2*BsgL*DLc*m3*m4**2*PsL - Bvm2R*DLc*m1**2*m3*m4**2*PsL + Bvm2R*DLc*m2**2*m3*m4**2*PsL + Bsm2L*DLc*m3**2*m4**2*PsL + 
     &           Bvm2R*DLc*m3**3*m4**2*PsL + Bsm2L*DLc*m4**4*PsL - Bvm2R*DLc*m3*m4**4*PsL - 2*Bvm1L*DLc*m1**2*m2**2*mq1*PsL + 2*BvgL*DLc*m1**2*m3*mq1*PsL + 
     &           3*Bvm1L*DLc*m1**2*m3**2*mq1*PsL + Bvm1L*DLc*m2**2*m3**2*mq1*PsL + Bvm1L*DLc*m3**4*mq1*PsL - Bvm2L*DLc*m1**2*m4**2*mq1*PsL - 
     &           Bvm2L*DLc*m2**2*m4**2*mq1*PsL + 2*BvgL*DLc*m3*m4**2*mq1*PsL + Bvm1L*DLc*m3**2*m4**2*mq1*PsL - Bvm2L*DLc*m3**2*m4**2*mq1*PsL - 
     &           Bvm2L*DLc*m4**4*mq1*PsL - 2*Bsm2L*DRc*m1*m4**2*mq3c*PsL - 2*Bvm2R*DRc*m1*m3*m4**2*mq3c*PsL - 4*Bvm1L*DRc*m1*m3**2*mq1*mq3c*PsL - 
     &           2*Bvm1L*DRc*m1*m4**2*mq1*mq3c*PsL + 2*Bvm2L*DRc*m1*m4**2*mq1*mq3c*PsL + 
     &           (2*BvgL*DLc*m1**2*m2**2*m3 - 2*Bvm1L*DLc*m1**2*m2**2*m3**2 + 2*BvgL*DLc*m1**2*m3**3 + 2*Bvm1L*DLc*m1**2*m3**4 + 
     &              2*Bvm2L*DLc*m1**2*m2**2*m4**2 - Bsm2R*DLc*m1**2*m3*m4**2 + 4*BvgL*DLc*m1**2*m3*m4**2 + Bsm2R*DLc*m2**2*m3*m4**2 + 
     &              2*BvgL*DLc*m2**2*m3*m4**2 - 2*Bvm2L*DLc*m1**2*m3**2*m4**2 + Bsm2R*DLc*m3**3*m4**2 + 2*BvgL*DLc*m3**3*m4**2 - Bsm2R*DLc*m3*m4**4 + 
     &              4*BvgL*DLc*m3*m4**4 + BsgR*DLc*(-(m3**2*(m2**2 + m3**2)) + (m2**2 + 2*m3**2)*m4**2 + m4**4 + m1**2*(2*m2**2 - m3**2 + m4**2)) + 
     &              Bvm2R*DLc*m1**2*m3*m4**2*mq1 - Bvm2R*DLc*m2**2*m3*m4**2*mq1 - Bvm2R*DLc*m3**3*m4**2*mq1 + Bvm2R*DLc*m3*m4**4*mq1 + 
     &              4*BsgR*DRc*m1*m3**2*mq3c + 2*DRc*m1*(-(Bvm1L*m3**2*(m3**2 + m4**2)) + m4**2*(-(Bsm2R*m3) + Bvm2L*(m3**2 + m4**2) + Bvm2R*m3*mq1))*mq3c + 
     &              Bsm1R*m3*(DLc*(-(m3**2*(m2**2 + m3**2)) + (4*m2**2 + 7*m3**2)*m4**2 + 6*m4**4 + m1**2*(4*m2**2 + m3**2 + 8*m4**2)) + 
     &                 2*DRc*m1*(2*m3**2 - m4**2)*mq3c))*PvR + BvgR*
     &            (DLc*(-(m3**2*(m2**2 + m3**2)) + (m2**2 + 2*m3**2)*m4**2 + m4**4 + m1**2*(2*m2**2 - m3**2 + m4**2)) + 4*DRc*m1*m3**2*mq3c)*(PsL - mq1*PvR))*
     &         sman + Bvm1R*m3*(PsL - mq1*PvR)*(2*DLc*(2*m1**2*m2**2*m3**2 + 2*m3**2*(2*m1**2 + m2**2 + m3**2)*m4**2 + (2*m1**2 + m2**2 + 3*m3**2)*m4**4) - 
     &           DLc*(-(m3**2*(m2**2 + m3**2)) + (4*m2**2 + 7*m3**2)*m4**2 + 6*m4**4 + m1**2*(4*m2**2 + m3**2 + 8*m4**2))*sman + 
     &           2*DRc*m1*mq3c*(m3**4 - 2*m4**4 + m4**2*sman - m3**2*(m4**2 + 2*sman))) + 
     &        Bsm1L*PsL*(2*DRc*m1*mq3c*(m3**4 + m3**2*(m4**2 - 2*sman) - m4**2*sman) + 
     &           DLc*(m3**2*(m2**2 + m3**2 + m4**2)*sman + m1**2*(-2*m3**4 - 2*m2**2*sman + 3*m3**2*sman))))*xv + 
     &     DLc*tman*(-(Bvm2L*m1**2*m4**2*mq1*PsL) - Bvm2L*m3**2*m4**2*mq1*PsL - 2*BsgR*m1**2*m3**2*PvR - 2*Bsm1R*m1**2*m3**3*PvR - BsgR*m3**4*PvR - 
     &        Bsm1R*m3**5*PvR + 2*BsgR*m1**2*m4**2*PvR + 2*Bsm2R*m1**2*m3*m4**2*PvR - Bsm1R*m3**3*m4**2*PvR + Bsm2R*m3**3*m4**2*PvR + BsgR*m4**4*PvR + 
     &        Bsm2R*m3*m4**4*PvR + 2*BvgR*m1**2*m3**2*mq1*PvR + 2*Bvm1R*m1**2*m3**3*mq1*PvR + BvgR*m3**4*mq1*PvR + Bvm1R*m3**5*mq1*PvR - 
     &        2*BvgR*m1**2*m4**2*mq1*PvR - 2*Bvm2R*m1**2*m3*m4**2*mq1*PvR + Bvm1R*m3**3*m4**2*mq1*PvR - Bvm2R*m3**3*m4**2*mq1*PvR - BvgR*m4**4*mq1*PvR - 
     &        Bvm2R*m3*m4**4*mq1*PvR + Bsm1L*m1**2*PsL*sman + 2*Bvm1R*m1**2*m3*PsL*sman + Bsm1L*m3**2*PsL*sman + BvgR*m3**2*PsL*sman + Bvm1R*m3**3*PsL*sman - 
     &        BvgR*m4**2*PsL*sman + Bvm1R*m3*m4**2*PsL*sman + 2*Bsm1R*m1**2*m3*PvR*sman + BsgR*m3**2*PvR*sman + Bsm1R*m3**3*PvR*sman - BsgR*m4**2*PvR*sman + 
     &        Bvm2L*m1**2*m4**2*PvR*sman + Bsm1R*m3*m4**2*PvR*sman + Bvm2L*m3**2*m4**2*PvR*sman - 2*Bvm1R*m1**2*m3*mq1*PvR*sman - BvgR*m3**2*mq1*PvR*sman - 
     &        Bvm1R*m3**3*mq1*PvR*sman + BvgR*m4**2*mq1*PvR*sman - Bvm1R*m3*m4**2*mq1*PvR*sman - 2*BvgL*m3*(m1**2 + m4**2)*(mq1*PsL - PvR*sman) + 
     &        Bvm1L*(m3**2 - sman)*(mq1*PsL - PvR*sman)*(m1**2 + m3**2 - tman) + 
     &        (Bvm1R*m3**3*PsL - Bsm2L*m4**2*PsL - Bvm2R*m3*m4**2*PsL + Bvm2L*m4**2*mq1*PsL + BsgR*m3**2*PvR + Bsm1R*m3**3*PvR - BsgR*m4**2*PvR - 
     &           Bsm2R*m3*m4**2*PvR - Bvm1R*m3**3*mq1*PvR + Bvm2R*m3*m4**2*mq1*PvR + Bsm1L*PsL*(m3**2 - sman) + BvgR*(PsL - mq1*PvR)*(m3**2 - m4**2 - sman) - 
     &           ((BsgR + Bsm1R*m3 + Bvm2L*m4**2)*PvR + Bvm1R*m3*(PsL - mq1*PvR))*sman)*tman)*xv + 
     &     (Bvm1R*DLc*m1**2*m3*PsL*sman**2 - Bvm1R*DLc*m2**2*m3*PsL*sman**2 - 2*Bvm1R*DLc*m3**3*PsL*sman**2 + Bsm2L*DLc*m4**2*PsL*sman**2 + 
     &        3*Bvm1R*DLc*m3*m4**2*PsL*sman**2 + Bvm2R*DLc*m3*m4**2*PsL*sman**2 + Bvm1L*DLc*m1**2*mq1*PsL*sman**2 + Bvm1L*DLc*m2**2*mq1*PsL*sman**2 + 
     &        2*Bvm1L*DLc*m3**2*mq1*PsL*sman**2 + Bvm1L*DLc*m4**2*mq1*PsL*sman**2 - Bvm2L*DLc*m4**2*mq1*PsL*sman**2 + 2*Bvm1R*DRc*m1*m3*mq3c*PsL*sman**2 - 
     &        2*Bvm1L*DRc*m1*mq1*mq3c*PsL*sman**2 - BsgR*DLc*m1**2*PvR*sman**2 - BsgR*DLc*m2**2*PvR*sman**2 - 2*Bvm1L*DLc*m1**2*m2**2*PvR*sman**2 + 
     &        Bsm1R*DLc*m1**2*m3*PvR*sman**2 + 2*BvgL*DLc*m1**2*m3*PvR*sman**2 - Bsm1R*DLc*m2**2*m3*PvR*sman**2 - 2*BsgR*DLc*m3**2*PvR*sman**2 + 
     &        3*Bvm1L*DLc*m1**2*m3**2*PvR*sman**2 + Bvm1L*DLc*m2**2*m3**2*PvR*sman**2 - 2*Bsm1R*DLc*m3**3*PvR*sman**2 + Bvm1L*DLc*m3**4*PvR*sman**2 - 
     &        Bvm2L*DLc*m1**2*m4**2*PvR*sman**2 - Bvm2L*DLc*m2**2*m4**2*PvR*sman**2 + 3*Bsm1R*DLc*m3*m4**2*PvR*sman**2 + Bsm2R*DLc*m3*m4**2*PvR*sman**2 + 
     &        2*BvgL*DLc*m3*m4**2*PvR*sman**2 + Bvm1L*DLc*m3**2*m4**2*PvR*sman**2 - Bvm2L*DLc*m3**2*m4**2*PvR*sman**2 - Bvm2L*DLc*m4**4*PvR*sman**2 - 
     &        Bvm1R*DLc*m1**2*m3*mq1*PvR*sman**2 + Bvm1R*DLc*m2**2*m3*mq1*PvR*sman**2 + 2*Bvm1R*DLc*m3**3*mq1*PvR*sman**2 - 
     &        3*Bvm1R*DLc*m3*m4**2*mq1*PvR*sman**2 - Bvm2R*DLc*m3*m4**2*mq1*PvR*sman**2 + 2*BsgR*DRc*m1*mq3c*PvR*sman**2 + 
     &        2*Bsm1R*DRc*m1*m3*mq3c*PvR*sman**2 - 4*Bvm1L*DRc*m1*m3**2*mq3c*PvR*sman**2 - 2*Bvm1L*DRc*m1*m4**2*mq3c*PvR*sman**2 + 
     &        2*Bvm2L*DRc*m1*m4**2*mq3c*PvR*sman**2 - 2*Bvm1R*DRc*m1*m3*mq1*mq3c*PvR*sman**2 + 
     &        BvgR*(PsL - mq1*PvR)*(2*DRc*m1*mq3c - DLc*(m1**2 + m2**2 + 2*m3**2 - sman))*sman**2 + Bvm1R*DLc*m3*PsL*sman**3 - Bvm1L*DLc*mq1*PsL*sman**3 + 
     &        BsgR*DLc*PvR*sman**3 - Bvm1L*DLc*m1**2*PvR*sman**3 - Bvm1L*DLc*m2**2*PvR*sman**3 + Bsm1R*DLc*m3*PvR*sman**3 - 2*Bvm1L*DLc*m3**2*PvR*sman**3 - 
     &        Bvm1L*DLc*m4**2*PvR*sman**3 + Bvm2L*DLc*m4**2*PvR*sman**3 - Bvm1R*DLc*m3*mq1*PvR*sman**3 + 2*Bvm1L*DRc*m1*mq3c*PvR*sman**3 + 
     &        Bvm1L*DLc*PvR*sman**4 + BvgR*DLc*(-(((m1**2 + m2**2)*m3**2 + (m1**2 - m2**2 + 2*m3**2)*m4**2)*(PsL - mq1*PvR)) + 
     &           (m1**2 + m2**2 + 2*m4**2)*PsL*sman)*tman + DLc*
     &         (Bvm2R*m1**2*m3*m4**2*PsL + Bvm2R*m2**2*m3*m4**2*PsL + 2*Bvm2R*m3*m4**4*PsL + Bsm2L*m4**2*(m2**2 + m4**2)*PsL + Bvm1L*m2**2*m3**2*mq1*PsL - 
     &           Bvm2L*m2**2*m4**2*mq1*PsL + Bvm1L*m3**2*m4**2*mq1*PsL - Bvm2L*m4**4*mq1*PsL - BsgR*m1**2*m3**2*PvR - BsgR*m2**2*m3**2*PvR - 
     &           Bsm1R*m1**2*m3**3*PvR - Bsm1R*m2**2*m3**3*PvR - BsgR*m1**2*m4**2*PvR + BsgR*m2**2*m4**2*PvR + Bsm2R*m1**2*m3*m4**2*PvR + 
     &           Bsm2R*m2**2*m3*m4**2*PvR - 2*BsgR*m3**2*m4**2*PvR - 4*Bsm1R*m3**3*m4**2*PvR - 2*Bsm1R*m3*m4**4*PvR + 2*Bsm2R*m3*m4**4*PvR - 
     &           Bvm2R*m1**2*m3*m4**2*mq1*PvR - Bvm2R*m2**2*m3*m4**2*mq1*PvR - 2*Bvm2R*m3*m4**4*mq1*PvR - Bvm1L*(m2**2 + m4**2)*mq1*PsL*sman + 
     &           Bvm1R*m3*(-(((m1**2 + m2**2)*m3**2 + 4*m3**2*m4**2 + 2*m4**4)*(PsL - mq1*PvR)) + (m1**2 + m2**2 + 4*m4**2)*PsL*sman))*tman + 
     &        Bsm1L*PsL*(2*DRc*m1*mq3c*sman**2 - DLc*(m1**2 + m2**2 + 2*m3**2 + m4**2 - sman)*sman**2 + DLc*(m2**2 + m4**2)*(-m3**2 + sman)*tman))*xv + 
     &     DLc*(2*BsgL*m1**2*m3*PsL*tman - Bsm1L*m1**2*m3**2*PsL*tman - 2*BvgR*m1**2*m3**2*PsL*tman - 2*Bvm1R*m1**2*m3**3*PsL*tman - Bsm1L*m3**4*PsL*tman - 
     &        BvgR*m3**4*PsL*tman - Bvm1R*m3**5*PsL*tman + Bsm2L*m1**2*m4**2*PsL*tman + 2*BvgR*m1**2*m4**2*PsL*tman + 2*BsgL*m3*m4**2*PsL*tman + 
     &        2*Bvm2R*m1**2*m3*m4**2*PsL*tman + Bsm2L*m3**2*m4**2*PsL*tman - Bvm1R*m3**3*m4**2*PsL*tman + Bvm2R*m3**3*m4**2*PsL*tman + BvgR*m4**4*PsL*tman + 
     &        Bvm2R*m3*m4**4*PsL*tman + BsgR*m1**2*PvR*sman*tman + BsgR*m2**2*PvR*sman*tman + Bsm1R*m1**2*m3*PvR*sman*tman + Bsm1R*m2**2*m3*PvR*sman*tman - 
     &        Bvm1L*m2**2*m3**2*PvR*sman*tman + 2*BsgR*m4**2*PvR*sman*tman + Bvm2L*m2**2*m4**2*PvR*sman*tman + 4*Bsm1R*m3*m4**2*PvR*sman*tman - 
     &        Bvm1L*m3**2*m4**2*PvR*sman*tman + Bvm2L*m4**4*PvR*sman*tman - BvgR*m1**2*mq1*PvR*sman*tman - BvgR*m2**2*mq1*PvR*sman*tman - 
     &        Bvm1R*m1**2*m3*mq1*PvR*sman*tman - Bvm1R*m2**2*m3*mq1*PvR*sman*tman - 2*BvgR*m4**2*mq1*PvR*sman*tman - 4*Bvm1R*m3*m4**2*mq1*PvR*sman*tman + 
     &        Bvm1L*m2**2*PvR*sman**2*tman + Bvm1L*m4**2*PvR*sman**2*tman + 
     &        (-(BvgR*m1**2*m3**2*PsL) - 3*BvgR*m2**2*m3**2*PsL - Bvm1R*m1**2*m3**3*PsL - 3*Bvm1R*m2**2*m3**3*PsL - BvgR*m3**4*PsL - Bvm1R*m3**5*PsL + 
     &           Bsm2L*m1**2*m4**2*PsL - BvgR*m1**2*m4**2*PsL + Bsm2L*m2**2*m4**2*PsL - BvgR*m2**2*m4**2*PsL + Bvm2R*m1**2*m3*m4**2*PsL + 
     &           3*Bvm2R*m2**2*m3*m4**2*PsL - Bsm2L*m3**2*m4**2*PsL - 4*BvgR*m3**2*m4**2*PsL - 5*Bvm1R*m3**3*m4**2*PsL + Bvm2R*m3**3*m4**2*PsL + 
     &           Bsm2L*m4**4*PsL - BvgR*m4**4*PsL + 5*Bvm2R*m3*m4**4*PsL - 2*BvgL*m2**2*m3*mq1*PsL + Bvm1L*m1**2*m3**2*mq1*PsL + Bvm1L*m2**2*m3**2*mq1*PsL - 
     &           2*BvgL*m3**3*mq1*PsL - Bvm1L*m3**4*mq1*PsL - Bvm2L*m1**2*m4**2*mq1*PsL - Bvm2L*m2**2*m4**2*mq1*PsL - 4*BvgL*m3*m4**2*mq1*PsL + 
     &           Bvm1L*m3**2*m4**2*mq1*PsL + Bvm2L*m3**2*m4**2*mq1*PsL - Bvm2L*m4**4*mq1*PsL - BsgR*m1**2*m3**2*PvR - 3*BsgR*m2**2*m3**2*PvR - 
     &           Bsm1R*m1**2*m3**3*PvR - 3*Bsm1R*m2**2*m3**3*PvR - BsgR*m3**4*PvR - Bsm1R*m3**5*PvR - BsgR*m1**2*m4**2*PvR - BsgR*m2**2*m4**2*PvR + 
     &           Bsm2R*m1**2*m3*m4**2*PvR + 3*Bsm2R*m2**2*m3*m4**2*PvR - 4*BsgR*m3**2*m4**2*PvR - 5*Bsm1R*m3**3*m4**2*PvR + Bsm2R*m3**3*m4**2*PvR - 
     &           BsgR*m4**4*PvR + 5*Bsm2R*m3*m4**4*PvR + BvgR*m1**2*m3**2*mq1*PvR + 3*BvgR*m2**2*m3**2*mq1*PvR + Bvm1R*m1**2*m3**3*mq1*PvR + 
     &           3*Bvm1R*m2**2*m3**3*mq1*PvR + BvgR*m3**4*mq1*PvR + Bvm1R*m3**5*mq1*PvR + BvgR*m1**2*m4**2*mq1*PvR + BvgR*m2**2*m4**2*mq1*PvR - 
     &           Bvm2R*m1**2*m3*m4**2*mq1*PvR - 3*Bvm2R*m2**2*m3*m4**2*mq1*PvR + 4*BvgR*m3**2*m4**2*mq1*PvR + 5*Bvm1R*m3**3*m4**2*mq1*PvR - 
     &           Bvm2R*m3**3*m4**2*mq1*PvR + BvgR*m4**4*mq1*PvR - 5*Bvm2R*m3*m4**4*mq1*PvR - Bsm1L*(m1**2 + m2**2 - m3**2 + m4**2)*PsL*(m3**2 - sman) + 
     &           2*BsgL*m3*PsL*(m2**2 + m3**2 + 2*m4**2 - sman) + BvgR*m1**2*PsL*sman + BvgR*m2**2*PsL*sman + Bvm1R*m1**2*m3*PsL*sman + 
     &           3*Bvm1R*m2**2*m3*PsL*sman + BvgR*m3**2*PsL*sman + 3*Bvm1R*m3**3*PsL*sman + BvgR*m4**2*PsL*sman + 5*Bvm1R*m3*m4**2*PsL*sman - 
     &           2*Bvm2R*m3*m4**2*PsL*sman - Bvm1L*m1**2*mq1*PsL*sman - Bvm1L*m2**2*mq1*PsL*sman + 2*BvgL*m3*mq1*PsL*sman + Bvm1L*m3**2*mq1*PsL*sman - 
     &           Bvm1L*m4**2*mq1*PsL*sman + BsgR*m1**2*PvR*sman + BsgR*m2**2*PvR*sman + Bsm1R*m1**2*m3*PvR*sman + 3*Bsm1R*m2**2*m3*PvR*sman + 
     &           2*BvgL*m2**2*m3*PvR*sman + BsgR*m3**2*PvR*sman - Bvm1L*m1**2*m3**2*PvR*sman - Bvm1L*m2**2*m3**2*PvR*sman + 3*Bsm1R*m3**3*PvR*sman + 
     &           2*BvgL*m3**3*PvR*sman + Bvm1L*m3**4*PvR*sman + BsgR*m4**2*PvR*sman + Bvm2L*m1**2*m4**2*PvR*sman + Bvm2L*m2**2*m4**2*PvR*sman + 
     &           5*Bsm1R*m3*m4**2*PvR*sman - 2*Bsm2R*m3*m4**2*PvR*sman + 4*BvgL*m3*m4**2*PvR*sman - Bvm1L*m3**2*m4**2*PvR*sman - Bvm2L*m3**2*m4**2*PvR*sman + 
     &           Bvm2L*m4**4*PvR*sman - BvgR*m1**2*mq1*PvR*sman - BvgR*m2**2*mq1*PvR*sman - Bvm1R*m1**2*m3*mq1*PvR*sman - 3*Bvm1R*m2**2*m3*mq1*PvR*sman - 
     &           BvgR*m3**2*mq1*PvR*sman - 3*Bvm1R*m3**3*mq1*PvR*sman - BvgR*m4**2*mq1*PvR*sman - 5*Bvm1R*m3*m4**2*mq1*PvR*sman + 
     &           2*Bvm2R*m3*m4**2*mq1*PvR*sman - 2*Bvm1R*m3*PsL*sman**2 + Bvm1L*m1**2*PvR*sman**2 + Bvm1L*m2**2*PvR*sman**2 - 2*Bsm1R*m3*PvR*sman**2 - 
     &           2*BvgL*m3*PvR*sman**2 - Bvm1L*m3**2*PvR*sman**2 + Bvm1L*m4**2*PvR*sman**2 + 2*Bvm1R*m3*mq1*PvR*sman**2 + 
     &           (-(Bvm2R*m3*m4**2*PsL) + BsgR*m3**2*PvR + Bsm1R*m3**3*PvR + BsgR*m4**2*PvR - Bsm2R*m3*m4**2*PvR + Bvm2R*m3*m4**2*mq1*PvR + 
     &              Bvm1R*m3*(PsL - mq1*PvR)*(m3**2 - sman) + BvgR*(PsL - mq1*PvR)*(m3**2 + m4**2 - sman) - (BsgR + Bsm1R*m3)*PvR*sman)*tman)*uman + 
     &        (Bvm1R*m3**3*PsL - Bsm2L*m4**2*PsL - Bvm2R*m3*m4**2*PsL - Bvm1L*m3**2*mq1*PsL + Bvm2L*m4**2*mq1*PsL + BsgR*m3**2*PvR + Bsm1R*m3**3*PvR + 
     &           BsgR*m4**2*PvR - Bsm2R*m3*m4**2*PvR - Bvm1R*m3**3*mq1*PvR + Bvm2R*m3*m4**2*mq1*PvR + Bsm1L*PsL*(m3**2 - sman) + 
     &           BvgR*(PsL - mq1*PvR)*(m3**2 + m4**2 - sman) - (-(Bvm1L*mq1*PsL) + (BsgR + m3*(Bsm1R - Bvm1L*m3) + Bvm2L*m4**2)*PvR + 
     &              Bvm1R*m3*(PsL - mq1*PvR))*sman - Bvm1L*PvR*sman**2)*uman**2)*xv + 
     &     (4*BvgL*DLc*m1**2*m3*m4**2*mq1*PsL + 2*BvgL*DLc*m2**2*m3*m4**2*mq1*PsL + 2*BvgL*DLc*m3**3*m4**2*mq1*PsL - 4*BvgL*DRc*m1*m3*m4**2*mq1*mq3c*PsL + 
     &        4*BsgR*DLc*m1**2*m3**2*m4**2*PvR + 2*BsgR*DLc*m2**2*m3**2*m4**2*PvR + 2*BsgR*DLc*m3**4*m4**2*PvR - 2*BsgR*DLc*m2**2*m4**4*PvR - 
     &        4*Bsm1R*DLc*m1**2*m3*m4**4*PvR - 4*Bsm2R*DLc*m1**2*m3*m4**4*PvR - 2*Bsm1R*DLc*m2**2*m3*m4**4*PvR - 2*Bsm2R*DLc*m2**2*m3*m4**4*PvR + 
     &        2*BsgR*DLc*m3**2*m4**4*PvR - 2*Bsm1R*DLc*m3**3*m4**4*PvR - 2*Bsm2R*DLc*m3**3*m4**4*PvR - 2*BsgR*DRc*m1*m3**2*m4**2*mq3c*PvR + 
     &        2*BsgR*DRc*m1*m4**4*mq3c*PvR + 4*Bsm1R*DRc*m1*m3*m4**4*mq3c*PvR + 4*Bsm2R*DRc*m1*m3*m4**4*mq3c*PvR - 2*BvgL*DLc*m3*m4**2*mq1*PsL*sman + 
     &        BsgR*DLc*m1**2*m4**2*PvR*sman + BsgR*DLc*m2**2*m4**2*PvR*sman - 4*BvgL*DLc*m1**2*m3*m4**2*PvR*sman - 2*BvgL*DLc*m2**2*m3*m4**2*PvR*sman - 
     &        BsgR*DLc*m3**2*m4**2*PvR*sman - 2*BvgL*DLc*m3**3*m4**2*PvR*sman + BsgR*DLc*m4**4*PvR*sman + 2*Bsm1R*DLc*m3*m4**4*PvR*sman + 
     &        2*Bsm2R*DLc*m3*m4**4*PvR*sman - 2*BsgR*DRc*m1*m4**2*mq3c*PvR*sman + 4*BvgL*DRc*m1*m3*m4**2*mq3c*PvR*sman - BsgR*DLc*m4**2*PvR*sman**2 + 
     &        2*BvgL*DLc*m3*m4**2*PvR*sman**2 - 2*BvgL*DLc*m3*m4**2*mq1*PsL*tman - BsgR*DLc*m1**2*m4**2*PvR*tman - BsgR*DLc*m2**2*m4**2*PvR*tman - 
     &        3*BsgR*DLc*m3**2*m4**2*PvR*tman - BsgR*DLc*m4**4*PvR*tman + 2*Bsm1R*DLc*m3*m4**4*PvR*tman + 2*Bsm2R*DLc*m3*m4**4*PvR*tman + 
     &        2*BvgL*DLc*m3*m4**2*PvR*sman*tman + BsgR*DLc*m4**2*PvR*tman**2 + 
     &        DLc*(m3*(2*BvgL*mq1*PsL + Bsm1R*m3**2*PvR - Bsm2R*m4**2*PvR - (Bsm1R + 2*BvgL)*PvR*sman)*tman + 
     &           BsgR*PvR*(m4**2*(m1**2 + m2**2 - m3**2 + m4**2) + (m3**2 - m4**2 + sman)*tman))*uman - BsgR*DLc*m4**2*PvR*uman**2 + 
     &        Bvm2R*m3*m4**2*(PsL - mq1*PvR)*(2*m4**2*(2*DRc*m1*mq3c + DLc*(-2*m1**2 - m2**2 - m3**2 + sman + tman)) - DLc*tman*uman) + 
     &        Bvm1R*m3*(PsL - mq1*PvR)*(2*m4**4*(2*DRc*m1*mq3c + DLc*(-2*m1**2 - m2**2 - m3**2 + sman + tman)) + DLc*(m3**2 - sman)*tman*uman))*xv - 
     &     DLc*(3*Bsm2R*m1**2*m3*PvR*tman + Bsm2R*m2**2*m3*PvR*tman + Bsm2R*m3**3*PvR*tman + 3*Bsm2R*m3*m4**2*PvR*tman - 4*BvgR*m1**2*mq1*PvR*tman - 
     &        8*BvgR*m2**2*mq1*PvR*tman - 5*Bvm1R*m1**2*m3*mq1*PvR*tman - 3*Bvm2R*m1**2*m3*mq1*PvR*tman - 7*Bvm1R*m2**2*m3*mq1*PvR*tman - 
     &        Bvm2R*m2**2*m3*mq1*PvR*tman - 8*BvgR*m3**2*mq1*PvR*tman - 3*Bvm1R*m3**3*mq1*PvR*tman - Bvm2R*m3**3*mq1*PvR*tman - 4*BvgR*m4**2*mq1*PvR*tman - 
     &        5*Bvm1R*m3*m4**2*mq1*PvR*tman - 3*Bvm2R*m3*m4**2*mq1*PvR*tman + Bvm1L*m1**2*PvR*sman*tman + Bvm2L*m1**2*PvR*sman*tman + 
     &        5*Bvm1L*m2**2*PvR*sman*tman + Bvm2L*m2**2*PvR*sman*tman - 4*BvgL*m3*PvR*sman*tman + Bvm1L*m3**2*PvR*sman*tman + Bvm2L*m3**2*PvR*sman*tman + 
     &        Bvm1L*m4**2*PvR*sman*tman + Bvm2L*m4**2*PvR*sman*tman - Bsm1L*PsL*tman**2 - Bsm2L*PsL*tman**2 - 4*BvgR*PsL*tman**2 - 3*Bvm1R*m3*PsL*tman**2 - 
     &        Bvm2R*m3*PsL*tman**2 + Bvm1L*mq1*PsL*tman**2 + Bvm2L*mq1*PsL*tman**2 - Bsm2R*m3*PvR*tman**2 + 4*BvgR*mq1*PvR*tman**2 + 
     &        3*Bvm1R*m3*mq1*PvR*tman**2 + Bvm2R*m3*mq1*PvR*tman**2 - Bvm1L*PvR*sman*tman**2 - Bvm2L*PvR*sman*tman**2 + 5*Bsm1L*m1**2*PsL*uman + 
     &        Bsm2L*m1**2*PsL*uman + 4*BvgR*m1**2*PsL*uman + Bsm1L*m2**2*PsL*uman + Bsm2L*m2**2*PsL*uman + 7*Bvm1R*m1**2*m3*PsL*uman + 
     &        Bvm2R*m1**2*m3*PsL*uman + 5*Bvm1R*m2**2*m3*PsL*uman + 3*Bvm2R*m2**2*m3*PsL*uman + 3*Bsm1L*m3**2*PsL*uman - Bsm2L*m3**2*PsL*uman + 
     &        4*BvgR*m3**2*PsL*uman + 3*Bvm1R*m3**3*PsL*uman + Bvm2R*m3**3*PsL*uman + Bsm1L*m4**2*PsL*uman + Bsm2L*m4**2*PsL*uman + 
     &        7*Bvm1R*m3*m4**2*PsL*uman + 5*Bvm2R*m3*m4**2*PsL*uman - 5*Bvm1L*m1**2*mq1*PsL*uman - Bvm2L*m1**2*mq1*PsL*uman - Bvm1L*m2**2*mq1*PsL*uman - 
     &        Bvm2L*m2**2*mq1*PsL*uman - 3*Bvm1L*m3**2*mq1*PsL*uman + Bvm2L*m3**2*mq1*PsL*uman - Bvm1L*m4**2*mq1*PsL*uman - Bvm2L*m4**2*mq1*PsL*uman + 
     &        Bsm2R*m1**2*m3*PvR*uman + 3*Bsm2R*m2**2*m3*PvR*uman + Bsm2R*m3**3*PvR*uman + 5*Bsm2R*m3*m4**2*PvR*uman - 4*BvgR*m1**2*mq1*PvR*uman - 
     &        7*Bvm1R*m1**2*m3*mq1*PvR*uman - Bvm2R*m1**2*m3*mq1*PvR*uman - 5*Bvm1R*m2**2*m3*mq1*PvR*uman - 3*Bvm2R*m2**2*m3*mq1*PvR*uman - 
     &        4*BvgR*m3**2*mq1*PvR*uman - 3*Bvm1R*m3**3*mq1*PvR*uman - Bvm2R*m3**3*mq1*PvR*uman - 7*Bvm1R*m3*m4**2*mq1*PvR*uman - 
     &        5*Bvm2R*m3*m4**2*mq1*PvR*uman - 2*Bvm1R*m3*PsL*sman*uman - 2*Bvm2R*m3*PsL*sman*uman + 5*Bvm1L*m1**2*PvR*sman*uman + Bvm2L*m1**2*PvR*sman*uman + 
     &        Bvm1L*m2**2*PvR*sman*uman + Bvm2L*m2**2*PvR*sman*uman - 2*Bsm2R*m3*PvR*sman*uman + 3*Bvm1L*m3**2*PvR*sman*uman - Bvm2L*m3**2*PvR*sman*uman + 
     &        Bvm1L*m4**2*PvR*sman*uman + Bvm2L*m4**2*PvR*sman*uman + 2*Bvm1R*m3*mq1*PvR*sman*uman + 2*Bvm2R*m3*mq1*PvR*sman*uman - 4*Bsm1L*PsL*tman*uman - 
     &        4*BvgR*PsL*tman*uman - 6*Bvm1R*m3*PsL*tman*uman - 2*Bvm2R*m3*PsL*tman*uman + 4*Bvm1L*mq1*PsL*tman*uman - 2*Bsm2R*m3*PvR*tman*uman + 
     &        4*BvgR*mq1*PvR*tman*uman + 6*Bvm1R*m3*mq1*PvR*tman*uman + 2*Bvm2R*m3*mq1*PvR*tman*uman - 4*Bvm1L*PvR*sman*tman*uman - Bsm1L*PsL*uman**2 - 
     &        Bsm2L*PsL*uman**2 - 3*Bvm1R*m3*PsL*uman**2 - Bvm2R*m3*PsL*uman**2 + Bvm1L*mq1*PsL*uman**2 + Bvm2L*mq1*PsL*uman**2 - Bsm2R*m3*PvR*uman**2 + 
     &        3*Bvm1R*m3*mq1*PvR*uman**2 + Bvm2R*m3*mq1*PvR*uman**2 - Bvm1L*PvR*sman*uman**2 - Bvm2L*PvR*sman*uman**2 + 
     &        4*BsgR*PvR*((m1**2 + 2*m2**2 + 2*m3**2 + m4**2 - tman)*tman + (m1**2 + m3**2 - tman)*uman) + 
     &        Bsm1R*m3*PvR*((5*m1**2 + 7*m2**2 + 3*m3**2 + 5*m4**2 - 3*tman)*tman + (7*m1**2 + 5*m2**2 + 3*m3**2 + 7*m4**2 - 2*sman - 6*tman)*uman - 
     &           3*uman**2) + 2*m1**2*m3*(-((Bsm1L + 2*BvgR)*m2**2*m3) + BsgL*(m2**2 + m3**2))*PsL*xv) + 
     &     2*BsgL*m3*PsL*(4*DLc*(2*m1**2 + m2**2 + m3**2 - sman - tman) - DLc*(m4**2*(2*m1**2 + m2**2 + m3**2 - sman - tman) + tman*uman)*xv + 
     &        2*DRc*m1*mq3c*(-4 + m4**2*xv)) + BvgR*(PsL - mq1*PvR)*
     &      (-4*DLc*((m3**2 + m4**2 - sman - tman)*(2*m3**2 + sman - tman) + (-m3**2 + m4**2)*uman - uman**2 + m1**2*(4*m3**2 + sman - tman + uman) + 
     &           m2**2*(2*m3**2 - 2*m4**2 + sman - tman + uman)) + 
     &        DLc*(m4**2*(m3**2 + m4**2 - sman - tman)*(2*m3**2 + sman - tman) + (m4**4 - m4**2*tman + sman*tman + m3**2*(-m4**2 + tman))*uman - 
     &           m4**2*uman**2 + m1**2*m4**2*(4*m3**2 + sman - tman + uman) + m2**2*m4**2*(2*m3**2 - 2*m4**2 + sman - tman + uman))*xv - 
     &        2*DRc*m1*mq3c*(m3**2 - m4**2 + sman)*(-4 + m4**2*xv))) + 
     &  A1sR*A3Lc*(2*(4*BvgL*DLc*m1*m2**2*m3*PsR - Bvm1L*DLc*m1*m2**2*m3**2*PsR + Bvm2L*DLc*m1*m2**2*m3**2*PsR - Bvm1L*DLc*m1*m3**4*PsR + 
     &        Bvm2L*DLc*m1*m3**4*PsR - 2*Bvm1L*DLc*m1*m2**2*m4**2*PsR + 4*BvgL*DLc*m1*m3*m4**2*PsR + 2*Bvm2L*DLc*m1*m3**2*m4**2*PsR - 
     &        3*Bvm1R*DLc*m1*m2**2*m3*mq1*PsR - Bvm2R*DLc*m1*m2**2*m3*mq1*PsR + Bvm1R*DLc*m1*m3**3*mq1*PsR - Bvm2R*DLc*m1*m3**3*mq1*PsR - 
     &        2*Bvm1R*DLc*m1*m3*m4**2*mq1*PsR - 2*Bvm2R*DLc*m1*m3*m4**2*mq1*PsR - 2*BsgR*DRc*m1**2*mq3c*PsR - 2*BsgR*DRc*m3**2*mq3c*PsR - 
     &        4*Bvm1L*DRc*m1**2*m3**2*mq3c*PsR - 2*Bvm1L*DRc*m3**4*mq3c*PsR + Bvm1L*DRc*m1**2*m4**2*mq3c*PsR + Bvm2L*DRc*m1**2*m4**2*mq3c*PsR - 
     &        2*Bvm1L*DRc*m3**2*m4**2*mq3c*PsR + Bvm1L*DRc*m4**4*mq3c*PsR + Bvm2L*DRc*m4**4*mq3c*PsR + 2*BvgR*DRc*m1**2*mq1*mq3c*PsR + 
     &        3*Bvm1R*DRc*m1**2*m3*mq1*mq3c*PsR + Bvm2R*DRc*m1**2*m3*mq1*mq3c*PsR + 2*BvgR*DRc*m3**2*mq1*mq3c*PsR + 2*Bvm1R*DRc*m3**3*mq1*mq3c*PsR + 
     &        Bvm1R*DRc*m3*m4**2*mq1*mq3c*PsR + Bvm2R*DRc*m3*m4**2*mq1*mq3c*PsR + 4*BsgL*DLc*m1*m2**2*m3*PvL - Bsm1L*DLc*m1*m2**2*m3**2*PvL + 
     &        Bsm2L*DLc*m1*m2**2*m3**2*PvL - Bsm1L*DLc*m1*m3**4*PvL + Bsm2L*DLc*m1*m3**4*PvL - 2*Bsm1L*DLc*m1*m2**2*m4**2*PvL + 4*BsgL*DLc*m1*m3*m4**2*PvL + 
     &        2*Bsm2L*DLc*m1*m3**2*m4**2*PvL - 4*BvgL*DLc*m1*m2**2*m3*mq1*PvL + Bvm1L*DLc*m1*m2**2*m3**2*mq1*PvL - Bvm2L*DLc*m1*m2**2*m3**2*mq1*PvL + 
     &        Bvm1L*DLc*m1*m3**4*mq1*PvL - Bvm2L*DLc*m1*m3**4*mq1*PvL + 2*Bvm1L*DLc*m1*m2**2*m4**2*mq1*PvL - 4*BvgL*DLc*m1*m3*m4**2*mq1*PvL - 
     &        2*Bvm2L*DLc*m1*m3**2*m4**2*mq1*PvL - 4*Bsm1L*DRc*m1**2*m3**2*mq3c*PvL - 2*Bsm1L*DRc*m3**4*mq3c*PvL + Bsm1L*DRc*m1**2*m4**2*mq3c*PvL + 
     &        Bsm2L*DRc*m1**2*m4**2*mq3c*PvL - 2*Bsm1L*DRc*m3**2*m4**2*mq3c*PvL + Bsm1L*DRc*m4**4*mq3c*PvL + Bsm2L*DRc*m4**4*mq3c*PvL + 
     &        4*Bvm1L*DRc*m1**2*m3**2*mq1*mq3c*PvL + 2*Bvm1L*DRc*m3**4*mq1*mq3c*PvL - Bvm1L*DRc*m1**2*m4**2*mq1*mq3c*PvL - 
     &        Bvm2L*DRc*m1**2*m4**2*mq1*mq3c*PvL + 2*Bvm1L*DRc*m3**2*m4**2*mq1*mq3c*PvL - Bvm1L*DRc*m4**4*mq1*mq3c*PvL - Bvm2L*DRc*m4**4*mq1*mq3c*PvL + 
     &        2*Bvm1L*DLc*m1*m2**2*PsR*sman - 2*BvgL*DLc*m1*m3*PsR*sman - Bvm1L*DLc*m1*m3**2*PsR*sman - Bvm2L*DLc*m1*m3**2*PsR*sman + 
     &        Bvm1R*DLc*m1*m3*mq1*PsR*sman + Bvm2R*DLc*m1*m3*mq1*PsR*sman - Bvm1L*DRc*m1**2*mq3c*PsR*sman - Bvm2L*DRc*m1**2*mq3c*PsR*sman - 
     &        Bvm1L*DRc*m4**2*mq3c*PsR*sman - Bvm2L*DRc*m4**2*mq3c*PsR*sman + 2*Bsm1L*DLc*m1*m2**2*PvL*sman - 2*BsgL*DLc*m1*m3*PvL*sman + 
     &        3*Bvm1R*DLc*m1*m2**2*m3*PvL*sman + Bvm2R*DLc*m1*m2**2*m3*PvL*sman - Bsm1L*DLc*m1*m3**2*PvL*sman - Bsm2L*DLc*m1*m3**2*PvL*sman - 
     &        Bvm1R*DLc*m1*m3**3*PvL*sman + Bvm2R*DLc*m1*m3**3*PvL*sman + 2*Bvm1R*DLc*m1*m3*m4**2*PvL*sman + 2*Bvm2R*DLc*m1*m3*m4**2*PvL*sman - 
     &        2*Bvm1L*DLc*m1*m2**2*mq1*PvL*sman + 2*BvgL*DLc*m1*m3*mq1*PvL*sman + Bvm1L*DLc*m1*m3**2*mq1*PvL*sman + Bvm2L*DLc*m1*m3**2*mq1*PvL*sman - 
     &        Bsm1L*DRc*m1**2*mq3c*PvL*sman - Bsm2L*DRc*m1**2*mq3c*PvL*sman - 2*BvgR*DRc*m1**2*mq3c*PvL*sman - 3*Bvm1R*DRc*m1**2*m3*mq3c*PvL*sman - 
     &        Bvm2R*DRc*m1**2*m3*mq3c*PvL*sman - 2*BvgR*DRc*m3**2*mq3c*PvL*sman - 2*Bvm1R*DRc*m3**3*mq3c*PvL*sman - Bsm1L*DRc*m4**2*mq3c*PvL*sman - 
     &        Bsm2L*DRc*m4**2*mq3c*PvL*sman - Bvm1R*DRc*m3*m4**2*mq3c*PvL*sman - Bvm2R*DRc*m3*m4**2*mq3c*PvL*sman + Bvm1L*DRc*m1**2*mq1*mq3c*PvL*sman + 
     &        Bvm2L*DRc*m1**2*mq1*mq3c*PvL*sman + Bvm1L*DRc*m4**2*mq1*mq3c*PvL*sman + Bvm2L*DRc*m4**2*mq1*mq3c*PvL*sman - Bvm1R*DLc*m1*m3*PvL*sman**2 - 
     &        Bvm2R*DLc*m1*m3*PvL*sman**2 + Bsm2R*m3*PsR*(-(DRc*(m1**2 + m4**2)*mq3c) + DLc*m1*(m2**2 + m3**2 + 2*m4**2 - sman - tman)) + 
     &        DLc*m1*m3*(-2*BvgL + (Bvm1L - Bvm2L)*m3)*PsR*tman - 
     &        Bsm1R*m3*PsR*(DRc*(3*m1**2 + 2*m3**2 + m4**2)*mq3c + DLc*m1*(-3*m2**2 + m3**2 - 2*m4**2 + sman + tman))) - 
     &     2*(-((2*BsgR*DRc*mq3c*PsR + 2*Bsm1R*DRc*m3*mq3c*PsR + 2*Bvm1L*DRc*m3**2*mq3c*PsR - 2*BvgR*DRc*mq1*mq3c*PsR - 2*BsgL*DLc*m1*m3*PvL + 
     &             Bsm1L*DLc*m1*m3**2*PvL - Bsm2L*DLc*m1*m3**2*PvL + 2*BvgL*DLc*m1*m3*mq1*PvL - Bvm1L*DLc*m1*m3**2*mq1*PvL + Bvm2L*DLc*m1*m3**2*mq1*PvL + 
     &             2*Bsm1L*DRc*m3**2*mq3c*PvL - 2*Bvm1L*DRc*m3**2*mq1*mq3c*PvL + 2*BvgR*DRc*mq3c*PvL*sman)*tman) + 
     &        (-2*Bvm1L*DLc*m1*m4**2*PsR - Bsm2R*DRc*m3*mq3c*PsR - 2*Bvm1L*DRc*m3**2*mq3c*PsR + Bvm1L*DRc*m4**2*mq3c*PsR + Bvm2L*DRc*m4**2*mq3c*PsR + 
     &           Bsm1R*m3*(2*DLc*m1 - DRc*mq3c)*PsR + 2*BsgL*DLc*m1*m3*PvL - 2*Bsm1L*DLc*m1*m4**2*PvL + 2*Bvm1L*DLc*m1*m4**2*mq1*PvL - 
     &           2*Bsm1L*DRc*m3**2*mq3c*PvL + Bsm1L*DRc*m4**2*mq3c*PvL + Bsm2L*DRc*m4**2*mq3c*PvL + 2*Bvm1L*DRc*m3**2*mq1*mq3c*PvL - 
     &           Bvm1L*DRc*m4**2*mq1*mq3c*PvL - Bvm2L*DRc*m4**2*mq1*mq3c*PvL + 2*BvgL*DLc*m1*m3*(PsR - mq1*PvL) + 
     &           (2*Bsm1L*DLc*m1*PvL + Bvm1L*(2*DLc*m1 - DRc*mq3c)*(PsR - mq1*PvL) - DRc*mq3c*(Bvm2L*PsR + (Bsm1L + Bsm2L - Bvm2L*mq1)*PvL))*sman)*uman - 
     &        Bvm2R*m3*(mq1*PsR - PvL*sman)*(DLc*m1*tman - DRc*mq3c*uman) - 
     &        Bvm1R*m3*(mq1*PsR - PvL*sman)*(-(DRc*mq3c*(2*tman + uman)) + DLc*m1*(tman + 2*uman))) + 
     &     2*(Bvm1L*DLc*m1*m2**2*m3**4*PsR + Bvm1L*DLc*m1*m3**6*PsR - Bsm2R*DLc*m1*m2**2*m3*m4**2*PsR + Bvm1L*DLc*m1*m2**2*m3**2*m4**2*PsR - 
     &        Bvm2L*DLc*m1*m2**2*m3**2*m4**2*PsR - Bsm2R*DLc*m1*m3**3*m4**2*PsR + Bvm1L*DLc*m1*m3**4*m4**2*PsR - Bvm2L*DLc*m1*m3**4*m4**2*PsR + 
     &        Bvm1L*DLc*m1*m2**2*m4**4*PsR - 2*Bsm2R*DLc*m1*m3*m4**4*PsR - Bvm1L*DLc*m1*m3**2*m4**4*PsR - 2*Bvm2L*DLc*m1*m3**2*m4**4*PsR - 
     &        Bvm1R*DLc*m1*m2**2*m3**3*mq1*PsR - Bvm1R*DLc*m1*m3**5*mq1*PsR + Bvm2R*DLc*m1*m2**2*m3*m4**2*mq1*PsR - 2*Bvm1R*DLc*m1*m3**3*m4**2*mq1*PsR + 
     &        Bvm2R*DLc*m1*m3**3*m4**2*mq1*PsR + 2*Bvm2R*DLc*m1*m3*m4**4*mq1*PsR - BsgR*DRc*m1**2*m3**2*mq3c*PsR - BsgR*DRc*m1**2*m4**2*mq3c*PsR + 
     &        Bsm2R*DRc*m1**2*m3*m4**2*mq3c*PsR - BsgR*DRc*m3**2*m4**2*mq3c*PsR + 2*Bvm1L*DRc*m1**2*m3**2*m4**2*mq3c*PsR + Bvm1L*DRc*m3**4*m4**2*mq3c*PsR - 
     &        BsgR*DRc*m4**4*mq3c*PsR + Bvm1L*DRc*m1**2*m4**4*mq3c*PsR - Bvm2L*DRc*m1**2*m4**4*mq3c*PsR + Bsm2R*DRc*m3*m4**4*mq3c*PsR + 
     &        2*Bvm1L*DRc*m3**2*m4**4*mq3c*PsR - Bvm2L*DRc*m4**6*mq3c*PsR + BvgR*DRc*m1**2*m3**2*mq1*mq3c*PsR + Bvm1R*DRc*m1**2*m3**3*mq1*mq3c*PsR + 
     &        BvgR*DRc*m1**2*m4**2*mq1*mq3c*PsR - Bvm2R*DRc*m1**2*m3*m4**2*mq1*mq3c*PsR + BvgR*DRc*m3**2*m4**2*mq1*mq3c*PsR + 
     &        Bvm1R*DRc*m3**3*m4**2*mq1*mq3c*PsR + BvgR*DRc*m4**4*mq1*mq3c*PsR - Bvm2R*DRc*m3*m4**4*mq1*mq3c*PsR + BsgL*DLc*m1*m2**2*m3**3*PvL + 
     &        Bsm1L*DLc*m1*m2**2*m3**4*PvL + BsgL*DLc*m1*m3**5*PvL + Bsm1L*DLc*m1*m3**6*PvL - BsgL*DLc*m1*m2**2*m3*m4**2*PvL + 
     &        Bsm1L*DLc*m1*m2**2*m3**2*m4**2*PvL - Bsm2L*DLc*m1*m2**2*m3**2*m4**2*PvL + BsgL*DLc*m1*m3**3*m4**2*PvL + Bsm1L*DLc*m1*m3**4*m4**2*PvL - 
     &        Bsm2L*DLc*m1*m3**4*m4**2*PvL + Bsm1L*DLc*m1*m2**2*m4**4*PvL - 2*BsgL*DLc*m1*m3*m4**4*PvL - Bsm1L*DLc*m1*m3**2*m4**4*PvL - 
     &        2*Bsm2L*DLc*m1*m3**2*m4**4*PvL - Bvm1L*DLc*m1*m2**2*m3**4*mq1*PvL - Bvm1L*DLc*m1*m3**6*mq1*PvL - Bvm1L*DLc*m1*m2**2*m3**2*m4**2*mq1*PvL + 
     &        Bvm2L*DLc*m1*m2**2*m3**2*m4**2*mq1*PvL - Bvm1L*DLc*m1*m3**4*m4**2*mq1*PvL + Bvm2L*DLc*m1*m3**4*m4**2*mq1*PvL - 
     &        Bvm1L*DLc*m1*m2**2*m4**4*mq1*PvL + Bvm1L*DLc*m1*m3**2*m4**4*mq1*PvL + 2*Bvm2L*DLc*m1*m3**2*m4**4*mq1*PvL + 
     &        2*Bsm1L*DRc*m1**2*m3**2*m4**2*mq3c*PvL + Bsm1L*DRc*m3**4*m4**2*mq3c*PvL + Bsm1L*DRc*m1**2*m4**4*mq3c*PvL - Bsm2L*DRc*m1**2*m4**4*mq3c*PvL + 
     &        2*Bsm1L*DRc*m3**2*m4**4*mq3c*PvL - Bsm2L*DRc*m4**6*mq3c*PvL - 2*Bvm1L*DRc*m1**2*m3**2*m4**2*mq1*mq3c*PvL - Bvm1L*DRc*m3**4*m4**2*mq1*mq3c*PvL - 
     &        Bvm1L*DRc*m1**2*m4**4*mq1*mq3c*PvL + Bvm2L*DRc*m1**2*m4**4*mq1*mq3c*PvL - 2*Bvm1L*DRc*m3**2*m4**4*mq1*mq3c*PvL + Bvm2L*DRc*m4**6*mq1*mq3c*PvL - 
     &        (-(Bsm2R*DLc*m1*m3*m4**2*PsR) - (DLc*m1*m3*(Bvm2L*m3*m4**2 + (-(Bvm2R*m4**2) + Bvm1R*(m2**2 + 2*(m3**2 + m4**2)))*mq1) + 
     &              DRc*(m1**2 + m4**2)*(BsgR + Bvm2L*m4**2 - (BvgR + Bvm1R*m3)*mq1)*mq3c)*PsR + 
     &           (DLc*m1*(Bsm1L*m3**2*(m2**2 + 2*m3**2) + Bsm1L*(m2**2 + m3**2)*m4**2 + BsgL*m3*(m2**2 + 2*m3**2 + m4**2) + 
     &                 m3*(-(Bvm1R*m3**2*(m2**2 + m3**2 + 2*m4**2)) + m4**2*(-(Bsm2L*m3) + Bvm2R*(m2**2 + m3**2 + 2*m4**2) + Bvm2L*m3*mq1))) + 
     &              DRc*(Bsm1L*(m1**2*m3**2 + 2*(m1**2 + m3**2)*m4**2 + m4**4) + 
     &                 (m1**2 + m4**2)*(Bvm1R*m3**3 + BvgR*(m3**2 + m4**2) - m4**2*(Bsm2L + Bvm2R*m3 - Bvm2L*mq1)))*mq3c)*PvL + 
     &           Bvm1L*(DLc*m1*(m2**2*(m3**2 + m4**2) + m3**2*(2*m3**2 + m4**2)) + DRc*(m1**2*m3**2 + 2*(m1**2 + m3**2)*m4**2 + m4**4)*mq3c)*(PsR - mq1*PvL))*
     &         sman + BvgL*DLc*m1*m3*(PsR - mq1*PvL)*((m3 - m4)*(m3 + m4)*(m2**2 + m3**2 + 2*m4**2) - (m2**2 + 2*m3**2 + m4**2)*sman) + 
     &        Bsm1R*m3*PsR*(-(DRc*(m1**2 + m4**2)*mq3c*(m3**2 - sman)) + DLc*m1*(m3**2*(m2**2 + m3**2 + 2*m4**2) - (m2**2 + 2*(m3**2 + m4**2))*sman)))*xv + 
     &     DLc*m1*(-(Bvm1L*m3**4*PsR) + 2*Bsm2R*m3*m4**2*PsR + Bvm1L*m3**2*m4**2*PsR + Bvm2L*m3**2*m4**2*PsR - Bvm2L*m4**4*PsR + 2*Bvm1R*m3**3*mq1*PsR - 
     &        2*Bvm2R*m3*m4**2*mq1*PsR - 2*BsgL*m3**3*PvL - Bsm1L*m3**4*PvL + 2*BsgL*m3*m4**2*PvL + Bsm1L*m3**2*m4**2*PvL + Bsm2L*m3**2*m4**2*PvL - 
     &        Bsm2L*m4**4*PvL + Bvm1L*m3**4*mq1*PvL - Bvm1L*m3**2*m4**2*mq1*PvL - Bvm2L*m3**2*m4**2*mq1*PvL + Bvm2L*m4**4*mq1*PvL - 
     &        2*BvgL*m3*(PsR - mq1*PvL)*(m3**2 - m4**2 - sman) + 
     &        ((Bvm2L*m4**2 - 2*Bvm1R*m3*mq1)*PsR + (2*BsgL*m3 - 2*Bvm1R*m3**3 + m4**2*(-Bsm1L + Bsm2L + 2*Bvm2R*m3 - Bvm2L*mq1))*PvL + 
     &           Bvm1L*m4**2*(-PsR + mq1*PvL))*sman + (Bvm1L*PsR + (Bsm1L + 2*Bvm1R*m3 - Bvm1L*mq1)*PvL)*sman**2 + 2*Bsm1R*m3*PsR*(-m3**2 + sman))*tman*xv + 
     &     (2*Bsm1R*DLc*m1*m3*PsR*sman**2 + 2*Bvm1L*DLc*m1*m3**2*PsR*sman**2 - 2*Bvm1R*DLc*m1*m3*mq1*PsR*sman**2 + 2*Bvm1L*DRc*m1**2*mq3c*PsR*sman**2 + 
     &        2*Bvm1L*DRc*m4**2*mq3c*PsR*sman**2 + 2*BsgL*DLc*m1*m3*PvL*sman**2 - 2*Bvm1R*DLc*m1*m2**2*m3*PvL*sman**2 + 2*Bsm1L*DLc*m1*m3**2*PvL*sman**2 - 
     &        4*Bvm1R*DLc*m1*m3**3*PvL*sman**2 - 4*Bvm1R*DLc*m1*m3*m4**2*PvL*sman**2 + 2*Bvm2R*DLc*m1*m3*m4**2*PvL*sman**2 - 
     &        2*Bvm1L*DLc*m1*m3**2*mq1*PvL*sman**2 + 2*Bsm1L*DRc*m1**2*mq3c*PvL*sman**2 + 2*BvgR*DRc*m1**2*mq3c*PvL*sman**2 + 
     &        2*Bvm1R*DRc*m1**2*m3*mq3c*PvL*sman**2 + 2*Bsm1L*DRc*m4**2*mq3c*PvL*sman**2 + 2*BvgR*DRc*m4**2*mq3c*PvL*sman**2 + 
     &        2*Bvm1R*DRc*m3*m4**2*mq3c*PvL*sman**2 - 2*Bvm1L*DRc*m1**2*mq1*mq3c*PvL*sman**2 - 2*Bvm1L*DRc*m4**2*mq1*mq3c*PvL*sman**2 + 
     &        2*BvgL*DLc*m1*m3*(PsR - mq1*PvL)*sman**2 + 2*Bvm1R*DLc*m1*m3*PvL*sman**3 + 
     &        (-(Bvm1L*(m3**2 + m4**2)*(DLc*m1*m3**2 + 2*DRc*m4**2*mq3c)*(PsR - mq1*PvL)) + 
     &           (m3**2 + m4**2)*((Bsm2L*DLc*m1*m4**2 - Bsm1L*(DLc*m1*m3**2 + 2*DRc*m4**2*mq3c))*PvL + Bvm2L*DLc*m1*m4**2*(PsR - mq1*PvL)) + 
     &           Bvm1L*(DLc*m1*(2*m3**2 + m4**2) + 2*DRc*m4**2*mq3c)*PsR*sman + DLc*m1*(2*Bsm1L*m3**2*PvL - m4**2*(Bvm2L*PsR + (-Bsm1L + Bsm2L)*PvL))*sman)*
     &         tman)*xv - (Bsm1L*PvL*sman*(-2*DRc*m4**2*mq3c + DLc*m1*sman)*tman + 
     &        Bvm1L*sman*(2*DRc*m4**2*mq1*mq3c*PvL + DLc*m1*(mq1*PvL*(2*m3**2 + m4**2 - sman) + PsR*sman))*tman + 
     &        2*Bvm1L*(PsR - mq1*PvL)*(DRc*mq3c*(m3**2 - sman)*(m4**2 - sman) + DLc*m1*m4**2*(m3**2 + m4**2 - sman))*uman + 
     &        2*(Bsm2R*DRc*m3*m4**2*mq3c*PsR + BvgR*DRc*m3**2*mq1*mq3c*PsR + Bvm1R*DRc*m3**3*mq1*mq3c*PsR + BvgR*DRc*m4**2*mq1*mq3c*PsR - 
     &           Bvm2R*DRc*m3*m4**2*mq1*mq3c*PsR + Bsm1L*DLc*m1*m3**2*m4**2*PvL + Bsm1L*DLc*m1*m4**4*PvL + Bsm1L*DRc*m3**2*m4**2*mq3c*PvL - 
     &           Bsm2L*DRc*m4**4*mq3c*PvL - BsgR*DRc*mq3c*PsR*(m3**2 + m4**2 - sman) - 
     &           (Bsm1L*DLc*m1*m4**2*PvL + DRc*mq3c*((BvgR + Bvm1R*m3)*mq1*PsR + 
     &                 (m3**2*(Bsm1L + BvgR + Bvm1R*m3) + (Bsm1L - Bsm2L + BvgR - Bvm2R*m3)*m4**2)*PvL))*sman + 
     &           DRc*(Bsm1L + BvgR + Bvm1R*m3)*mq3c*PvL*sman**2 + Bsm1R*DRc*m3*mq3c*PsR*(-m3**2 + sman))*uman + 
     &        Bvm2L*m4**2*(-(DLc*m1*mq1*PvL*sman*tman) - 2*DRc*mq3c*(PsR - mq1*PvL)*(m4**2 - sman)*uman))*xv + 
     &     2*(Bvm1L*DLc*m1*m2**2*m4**2*PsR + Bvm2L*DLc*m1*m2**2*m4**2*PsR - Bvm1L*DLc*m1*m3**2*m4**2*PsR - Bvm2L*DLc*m1*m3**2*m4**2*PsR - 
     &        4*BvgR*DLc*m1*m2**2*mq1*PsR + 4*BvgR*DLc*m1*m3**2*mq1*PsR + Bvm1L*DRc*m1**2*m4**2*mq3c*PsR + Bvm2L*DRc*m1**2*m4**2*mq3c*PsR + 
     &        Bvm1L*DRc*m3**2*m4**2*mq3c*PsR + Bvm2L*DRc*m3**2*m4**2*mq3c*PsR - 4*BvgR*DRc*m1**2*mq1*mq3c*PsR - 4*BvgR*DRc*m3**2*mq1*mq3c*PsR - 
     &        8*BsgL*DLc*m1*m2**2*m3*PvL + Bsm1L*DLc*m1*m2**2*m4**2*PvL + Bsm2L*DLc*m1*m2**2*m4**2*PvL - 8*BsgL*DLc*m1*m3*m4**2*PvL - 
     &        Bsm1L*DLc*m1*m3**2*m4**2*PvL - Bsm2L*DLc*m1*m3**2*m4**2*PvL - Bvm1L*DLc*m1*m2**2*m4**2*mq1*PvL - Bvm2L*DLc*m1*m2**2*m4**2*mq1*PvL + 
     &        Bvm1L*DLc*m1*m3**2*m4**2*mq1*PvL + Bvm2L*DLc*m1*m3**2*m4**2*mq1*PvL - 8*BsgL*DRc*m1**2*m3*mq3c*PvL - 4*BsgL*DRc*m3**3*mq3c*PvL + 
     &        Bsm1L*DRc*m1**2*m4**2*mq3c*PvL + Bsm2L*DRc*m1**2*m4**2*mq3c*PvL - 4*BsgL*DRc*m3*m4**2*mq3c*PvL + Bsm1L*DRc*m3**2*m4**2*mq3c*PvL + 
     &        Bsm2L*DRc*m3**2*m4**2*mq3c*PvL - Bvm1L*DRc*m1**2*m4**2*mq1*mq3c*PvL - Bvm2L*DRc*m1**2*m4**2*mq1*mq3c*PvL - Bvm1L*DRc*m3**2*m4**2*mq1*mq3c*PvL - 
     &        Bvm2L*DRc*m3**2*m4**2*mq1*mq3c*PvL + 4*BvgR*DLc*m1*m2**2*PvL*sman + 4*BsgL*DLc*m1*m3*PvL*sman - 4*BvgR*DLc*m1*m3**2*PvL*sman + 
     &        4*BvgR*DRc*m1**2*mq3c*PvL*sman + 4*BvgR*DRc*m3**2*mq3c*PvL*sman - Bvm1L*DRc*m4**2*mq3c*PsR*tman - Bvm2L*DRc*m4**2*mq3c*PsR*tman + 
     &        4*BvgR*DRc*mq1*mq3c*PsR*tman + 4*BsgL*DLc*m1*m3*PvL*tman + 4*BsgL*DRc*m3*mq3c*PvL*tman - Bsm1L*DRc*m4**2*mq3c*PvL*tman - 
     &        Bsm2L*DRc*m4**2*mq3c*PvL*tman + Bvm1L*DRc*m4**2*mq1*mq3c*PvL*tman + Bvm2L*DRc*m4**2*mq1*mq3c*PvL*tman - 4*BvgR*DRc*mq3c*PvL*sman*tman - 
     &        Bvm1L*DLc*m1*m4**2*PsR*uman - Bvm2L*DLc*m1*m4**2*PsR*uman + 4*BvgR*DLc*m1*mq1*PsR*uman + 4*BsgL*DLc*m1*m3*PvL*uman - 
     &        Bsm1L*DLc*m1*m4**2*PvL*uman - Bsm2L*DLc*m1*m4**2*PvL*uman + Bvm1L*DLc*m1*m4**2*mq1*PvL*uman + Bvm2L*DLc*m1*m4**2*mq1*PvL*uman + 
     &        4*BsgL*DRc*m3*mq3c*PvL*uman - 4*BvgR*DLc*m1*PvL*sman*uman - 
     &        m4**2*(-(BvgR*DLc*m1*m2**2*mq1*PsR) + BvgR*DLc*m1*m3**2*mq1*PsR - BvgR*DRc*m1**2*mq1*mq3c*PsR - BvgR*DRc*m3**2*mq1*mq3c*PsR - 
     &           2*BsgL*DLc*m1*m2**2*m3*PvL + Bsm1L*DLc*m1*m2**2*m4**2*PvL + Bsm2L*DLc*m1*m2**2*m4**2*PvL - 2*BsgL*DLc*m1*m3*m4**2*PvL - 
     &           Bsm1L*DLc*m1*m3**2*m4**2*PvL - Bsm2L*DLc*m1*m3**2*m4**2*PvL - 2*BsgL*DRc*m1**2*m3*mq3c*PvL - BsgL*DRc*m3**3*mq3c*PvL + 
     &           Bsm1L*DRc*m1**2*m4**2*mq3c*PvL + Bsm2L*DRc*m1**2*m4**2*mq3c*PvL - BsgL*DRc*m3*m4**2*mq3c*PvL + Bsm1L*DRc*m3**2*m4**2*mq3c*PvL + 
     &           Bsm2L*DRc*m3**2*m4**2*mq3c*PvL + BvgR*DLc*m1*m2**2*PvL*sman + BsgL*DLc*m1*m3*PvL*sman - BvgR*DLc*m1*m3**2*PvL*sman + 
     &           BvgR*DRc*m1**2*mq3c*PvL*sman + BvgR*DRc*m3**2*mq3c*PvL*sman + BvgR*DRc*mq1*mq3c*PsR*tman + BsgL*DLc*m1*m3*PvL*tman + 
     &           BsgL*DRc*m3*mq3c*PvL*tman - Bsm1L*DRc*m4**2*mq3c*PvL*tman - Bsm2L*DRc*m4**2*mq3c*PvL*tman - BvgR*DRc*mq3c*PvL*sman*tman + 
     &           Bvm1L*m4**2*(PsR - mq1*PvL)*(DRc*mq3c*(m1**2 + m3**2 - tman) + DLc*m1*(m2**2 - m3**2 - uman)) + 
     &           Bvm2L*m4**2*(PsR - mq1*PvL)*(DRc*mq3c*(m1**2 + m3**2 - tman) + DLc*m1*(m2**2 - m3**2 - uman)) + 
     &           ((-((Bsm1L + Bsm2L)*DLc*m1*m4**2) + BsgL*m3*(DLc*m1 + DRc*mq3c))*PvL + BvgR*DLc*m1*(mq1*PsR - PvL*sman))*uman)*xv - 
     &        BsgR*PsR*(DRc*mq3c*(m1**2 + m3**2 - tman) + DLc*m1*(m2**2 - m3**2 - uman))*(-4 + m4**2*xv) + 
     &        BvgL*m3*(PsR - mq1*PvL)*(DRc*mq3c*(2*m1**2 + m3**2 + m4**2 - tman - uman) + DLc*m1*(2*m2**2 + 2*m4**2 - sman - tman - uman))*(-4 + m4**2*xv))) + 
     &  A1sL*A3Rc*(-2*(-4*BvgR*DRc*m1*m2**2*m3*PsL + Bvm1R*DRc*m1*m2**2*m3**2*PsL - Bvm2R*DRc*m1*m2**2*m3**2*PsL + Bvm1R*DRc*m1*m3**4*PsL - 
     &        Bvm2R*DRc*m1*m3**4*PsL + 2*Bvm1R*DRc*m1*m2**2*m4**2*PsL - 4*BvgR*DRc*m1*m3*m4**2*PsL - 2*Bvm2R*DRc*m1*m3**2*m4**2*PsL + 
     &        3*Bvm1L*DRc*m1*m2**2*m3*mq1*PsL + Bvm2L*DRc*m1*m2**2*m3*mq1*PsL - Bvm1L*DRc*m1*m3**3*mq1*PsL + Bvm2L*DRc*m1*m3**3*mq1*PsL + 
     &        2*Bvm1L*DRc*m1*m3*m4**2*mq1*PsL + 2*Bvm2L*DRc*m1*m3*m4**2*mq1*PsL + 2*BsgL*DLc*m1**2*mq3c*PsL + 2*BsgL*DLc*m3**2*mq3c*PsL + 
     &        4*Bvm1R*DLc*m1**2*m3**2*mq3c*PsL + 2*Bvm1R*DLc*m3**4*mq3c*PsL - Bvm1R*DLc*m1**2*m4**2*mq3c*PsL - Bvm2R*DLc*m1**2*m4**2*mq3c*PsL + 
     &        2*Bvm1R*DLc*m3**2*m4**2*mq3c*PsL - Bvm1R*DLc*m4**4*mq3c*PsL - Bvm2R*DLc*m4**4*mq3c*PsL - 2*BvgL*DLc*m1**2*mq1*mq3c*PsL - 
     &        3*Bvm1L*DLc*m1**2*m3*mq1*mq3c*PsL - Bvm2L*DLc*m1**2*m3*mq1*mq3c*PsL - 2*BvgL*DLc*m3**2*mq1*mq3c*PsL - 2*Bvm1L*DLc*m3**3*mq1*mq3c*PsL - 
     &        Bvm1L*DLc*m3*m4**2*mq1*mq3c*PsL - Bvm2L*DLc*m3*m4**2*mq1*mq3c*PsL - 4*BsgR*DRc*m1*m2**2*m3*PvR + Bsm1R*DRc*m1*m2**2*m3**2*PvR - 
     &        Bsm2R*DRc*m1*m2**2*m3**2*PvR + Bsm1R*DRc*m1*m3**4*PvR - Bsm2R*DRc*m1*m3**4*PvR + 2*Bsm1R*DRc*m1*m2**2*m4**2*PvR - 4*BsgR*DRc*m1*m3*m4**2*PvR - 
     &        2*Bsm2R*DRc*m1*m3**2*m4**2*PvR + 4*BvgR*DRc*m1*m2**2*m3*mq1*PvR - Bvm1R*DRc*m1*m2**2*m3**2*mq1*PvR + Bvm2R*DRc*m1*m2**2*m3**2*mq1*PvR - 
     &        Bvm1R*DRc*m1*m3**4*mq1*PvR + Bvm2R*DRc*m1*m3**4*mq1*PvR - 2*Bvm1R*DRc*m1*m2**2*m4**2*mq1*PvR + 4*BvgR*DRc*m1*m3*m4**2*mq1*PvR + 
     &        2*Bvm2R*DRc*m1*m3**2*m4**2*mq1*PvR + 4*Bsm1R*DLc*m1**2*m3**2*mq3c*PvR + 2*Bsm1R*DLc*m3**4*mq3c*PvR - Bsm1R*DLc*m1**2*m4**2*mq3c*PvR - 
     &        Bsm2R*DLc*m1**2*m4**2*mq3c*PvR + 2*Bsm1R*DLc*m3**2*m4**2*mq3c*PvR - Bsm1R*DLc*m4**4*mq3c*PvR - Bsm2R*DLc*m4**4*mq3c*PvR - 
     &        4*Bvm1R*DLc*m1**2*m3**2*mq1*mq3c*PvR - 2*Bvm1R*DLc*m3**4*mq1*mq3c*PvR + Bvm1R*DLc*m1**2*m4**2*mq1*mq3c*PvR + 
     &        Bvm2R*DLc*m1**2*m4**2*mq1*mq3c*PvR - 2*Bvm1R*DLc*m3**2*m4**2*mq1*mq3c*PvR + Bvm1R*DLc*m4**4*mq1*mq3c*PvR + Bvm2R*DLc*m4**4*mq1*mq3c*PvR - 
     &        2*Bvm1R*DRc*m1*m2**2*PsL*sman + 2*BvgR*DRc*m1*m3*PsL*sman + Bvm1R*DRc*m1*m3**2*PsL*sman + Bvm2R*DRc*m1*m3**2*PsL*sman - 
     &        Bvm1L*DRc*m1*m3*mq1*PsL*sman - Bvm2L*DRc*m1*m3*mq1*PsL*sman + Bvm1R*DLc*m1**2*mq3c*PsL*sman + Bvm2R*DLc*m1**2*mq3c*PsL*sman + 
     &        Bvm1R*DLc*m4**2*mq3c*PsL*sman + Bvm2R*DLc*m4**2*mq3c*PsL*sman - 2*Bsm1R*DRc*m1*m2**2*PvR*sman + 2*BsgR*DRc*m1*m3*PvR*sman - 
     &        3*Bvm1L*DRc*m1*m2**2*m3*PvR*sman - Bvm2L*DRc*m1*m2**2*m3*PvR*sman + Bsm1R*DRc*m1*m3**2*PvR*sman + Bsm2R*DRc*m1*m3**2*PvR*sman + 
     &        Bvm1L*DRc*m1*m3**3*PvR*sman - Bvm2L*DRc*m1*m3**3*PvR*sman - 2*Bvm1L*DRc*m1*m3*m4**2*PvR*sman - 2*Bvm2L*DRc*m1*m3*m4**2*PvR*sman + 
     &        2*Bvm1R*DRc*m1*m2**2*mq1*PvR*sman - 2*BvgR*DRc*m1*m3*mq1*PvR*sman - Bvm1R*DRc*m1*m3**2*mq1*PvR*sman - Bvm2R*DRc*m1*m3**2*mq1*PvR*sman + 
     &        Bsm1R*DLc*m1**2*mq3c*PvR*sman + Bsm2R*DLc*m1**2*mq3c*PvR*sman + 2*BvgL*DLc*m1**2*mq3c*PvR*sman + 3*Bvm1L*DLc*m1**2*m3*mq3c*PvR*sman + 
     &        Bvm2L*DLc*m1**2*m3*mq3c*PvR*sman + 2*BvgL*DLc*m3**2*mq3c*PvR*sman + 2*Bvm1L*DLc*m3**3*mq3c*PvR*sman + Bsm1R*DLc*m4**2*mq3c*PvR*sman + 
     &        Bsm2R*DLc*m4**2*mq3c*PvR*sman + Bvm1L*DLc*m3*m4**2*mq3c*PvR*sman + Bvm2L*DLc*m3*m4**2*mq3c*PvR*sman - Bvm1R*DLc*m1**2*mq1*mq3c*PvR*sman - 
     &        Bvm2R*DLc*m1**2*mq1*mq3c*PvR*sman - Bvm1R*DLc*m4**2*mq1*mq3c*PvR*sman - Bvm2R*DLc*m4**2*mq1*mq3c*PvR*sman + Bvm1L*DRc*m1*m3*PvR*sman**2 + 
     &        Bvm2L*DRc*m1*m3*PvR*sman**2 + (-(DRc*m1*m3*(-2*BvgR + (Bvm1R - Bvm2R)*m3 + (Bvm1L + Bvm2L)*mq1)) - 
     &           2*DLc*(BsgL + Bvm1R*m3**2 - (BvgL + Bvm1L*m3)*mq1)*mq3c)*PsL*tman + 
     &        Bsm2L*m3*PsL*(DLc*(m1**2 + m4**2)*mq3c + DRc*m1*(-m2**2 - m3**2 - 2*m4**2 + sman + tman)) + 
     &        Bsm1L*m3*PsL*(DLc*mq3c*(3*m1**2 + 2*m3**2 + m4**2 - 2*tman) + DRc*m1*(-3*m2**2 + m3**2 - 2*m4**2 + sman + tman))) - 
     &     2*(-(Bsm1L*DRc*m1*m3**5*PsL) + Bsm2L*DRc*m1*m2**2*m3*m4**2*PsL + Bvm2R*DRc*m1*m2**2*m3**2*m4**2*PsL - 2*Bsm1L*DRc*m1*m3**3*m4**2*PsL + 
     &        Bsm2L*DRc*m1*m3**3*m4**2*PsL + Bvm2R*DRc*m1*m3**4*m4**2*PsL + 2*Bsm2L*DRc*m1*m3*m4**4*PsL + 2*Bvm2R*DRc*m1*m3**2*m4**4*PsL + 
     &        Bvm1L*DRc*m1*m2**2*m3**3*mq1*PsL + Bvm1L*DRc*m1*m3**5*mq1*PsL - Bvm2L*DRc*m1*m2**2*m3*m4**2*mq1*PsL + 2*Bvm1L*DRc*m1*m3**3*m4**2*mq1*PsL - 
     &        Bvm2L*DRc*m1*m3**3*m4**2*mq1*PsL - 2*Bvm2L*DRc*m1*m3*m4**4*mq1*PsL + BsgL*DLc*m1**2*m3**2*mq3c*PsL + Bsm1L*DLc*m1**2*m3**3*mq3c*PsL + 
     &        BsgL*DLc*m1**2*m4**2*mq3c*PsL - Bsm2L*DLc*m1**2*m3*m4**2*mq3c*PsL + BsgL*DLc*m3**2*m4**2*mq3c*PsL + Bsm1L*DLc*m3**3*m4**2*mq3c*PsL + 
     &        BsgL*DLc*m4**4*mq3c*PsL + Bvm2R*DLc*m1**2*m4**4*mq3c*PsL - Bsm2L*DLc*m3*m4**4*mq3c*PsL + Bvm2R*DLc*m4**6*mq3c*PsL - 
     &        BvgL*DLc*m1**2*m3**2*mq1*mq3c*PsL - Bvm1L*DLc*m1**2*m3**3*mq1*mq3c*PsL - BvgL*DLc*m1**2*m4**2*mq1*mq3c*PsL + 
     &        Bvm2L*DLc*m1**2*m3*m4**2*mq1*mq3c*PsL - BvgL*DLc*m3**2*m4**2*mq1*mq3c*PsL - Bvm1L*DLc*m3**3*m4**2*mq1*mq3c*PsL - BvgL*DLc*m4**4*mq1*mq3c*PsL + 
     &        Bvm2L*DLc*m3*m4**4*mq1*mq3c*PsL - BsgR*DRc*m1*m2**2*m3**3*PvR - Bsm1R*DRc*m1*m2**2*m3**4*PvR - BsgR*DRc*m1*m3**5*PvR - Bsm1R*DRc*m1*m3**6*PvR + 
     &        BsgR*DRc*m1*m2**2*m3*m4**2*PvR - Bsm1R*DRc*m1*m2**2*m3**2*m4**2*PvR + Bsm2R*DRc*m1*m2**2*m3**2*m4**2*PvR - BsgR*DRc*m1*m3**3*m4**2*PvR - 
     &        Bsm1R*DRc*m1*m3**4*m4**2*PvR + Bsm2R*DRc*m1*m3**4*m4**2*PvR - Bsm1R*DRc*m1*m2**2*m4**4*PvR + 2*BsgR*DRc*m1*m3*m4**4*PvR + 
     &        Bsm1R*DRc*m1*m3**2*m4**4*PvR + 2*Bsm2R*DRc*m1*m3**2*m4**4*PvR - Bvm2R*DRc*m1*m2**2*m3**2*m4**2*mq1*PvR - Bvm2R*DRc*m1*m3**4*m4**2*mq1*PvR - 
     &        2*Bvm2R*DRc*m1*m3**2*m4**4*mq1*PvR - 2*Bsm1R*DLc*m1**2*m3**2*m4**2*mq3c*PvR - Bsm1R*DLc*m3**4*m4**2*mq3c*PvR - Bsm1R*DLc*m1**2*m4**4*mq3c*PvR + 
     &        Bsm2R*DLc*m1**2*m4**4*mq3c*PvR - 2*Bsm1R*DLc*m3**2*m4**4*mq3c*PvR + Bsm2R*DLc*m4**6*mq3c*PvR - Bvm2R*DLc*m1**2*m4**4*mq1*mq3c*PvR - 
     &        Bvm2R*DLc*m4**6*mq1*mq3c*PvR + (-(Bsm2L*DRc*m1*m3*m4**2*PsL) + Bsm1L*m3*(DRc*m1*(m2**2 + 2*(m3**2 + m4**2)) - DLc*(m1**2 + m4**2)*mq3c)*PsL - 
     &           (DRc*m1*m3*(Bvm2R*m3*m4**2 + (-(Bvm2L*m4**2) + Bvm1L*(m2**2 + 2*(m3**2 + m4**2)))*mq1) + 
     &              DLc*(m1**2 + m4**2)*(BsgL + Bvm2R*m4**2 - (BvgL + Bvm1L*m3)*mq1)*mq3c)*PsL + 
     &           (DRc*m1*(Bsm1R*m3**2*(m2**2 + 2*m3**2) + Bsm1R*(m2**2 + m3**2)*m4**2 + BsgR*m3*(m2**2 + 2*m3**2 + m4**2) + 
     &                 m3*(-(Bvm1L*m3**2*(m2**2 + m3**2 + 2*m4**2)) + m4**2*(-(Bsm2R*m3) + Bvm2L*(m2**2 + m3**2 + 2*m4**2) + Bvm2R*m3*mq1))) + 
     &              DLc*((m1**2 + m4**2)*(m3**2*(BvgL + Bvm1L*m3) - (Bsm2R - BvgL + Bvm2L*m3)*m4**2) + Bsm1R*(m1**2*m3**2 + 2*(m1**2 + m3**2)*m4**2 + m4**4) + 
     &                 Bvm2R*m1**2*m4**2*mq1)*mq3c)*PvR)*sman + 
     &        BvgR*DRc*m1*m3*(PsL - mq1*PvR)*(-((m3 - m4)*(m3 + m4)*(m2**2 + m3**2 + 2*m4**2)) + (m2**2 + 2*m3**2 + m4**2)*sman) + 
     &        Bvm1R*(DLc*mq3c*(-(m1**2*(PsL - mq1*PvR)*(m4**4 + m3**2*(2*m4**2 - sman) - 2*m4**2*sman)) + 
     &              m4**2*(m3**4*(-PsL + mq1*PvR) - 2*m3**2*(PsL - mq1*PvR)*(m4**2 - sman) + m4**2*PsL*sman)) - 
     &           DRc*m1*(PsL - mq1*PvR)*(m2**2*(m3**4 + m3**2*m4**2 + m4**4 - (m3**2 + m4**2)*sman) + 
     &              m3**2*(m3**4 + m3**2*(m4**2 - 2*sman) - m4**2*(m4**2 + sman)))))*xv + 
     &     DRc*m1*(-((m3 - m4)*(m3 + m4)*(2*BsgR*m3 + Bsm1R*m3**2 - m3*(2*BvgR + Bvm1R*m3)*mq1 + m4**2*(-Bsm2R + Bvm2R*mq1))*PvR) + 
     &        (2*Bsm1L*m3*PsL + 2*BsgR*m3*PvR + 2*BvgR*m3*(PsL - mq1*PvR) + m4**2*((-Bsm1R + Bsm2R)*PvR - (Bvm1R - Bvm2R)*(PsL - mq1*PvR)))*sman + 
     &        (Bsm1R*PvR + Bvm1R*(PsL - mq1*PvR))*sman**2 + 2*Bvm1L*m3*(m3**2 - sman)*(mq1*PsL - PvR*sman) + 2*Bvm2L*m3*m4**2*(-(mq1*PsL) + PvR*sman))*tman*xv
     &      + (-2*Bvm2R*DLc*m4**4*mq1*mq3c*PvR*sman + 2*Bsm1L*DRc*m1*m3*PsL*sman**2 + 2*BvgR*DRc*m1*m3*PsL*sman**2 - 2*Bvm1L*DRc*m1*m3*mq1*PsL*sman**2 + 
     &        2*BsgR*DRc*m1*m3*PvR*sman**2 - 2*Bvm1L*DRc*m1*m2**2*m3*PvR*sman**2 + 2*Bsm1R*DRc*m1*m3**2*PvR*sman**2 - 4*Bvm1L*DRc*m1*m3**3*PvR*sman**2 - 
     &        4*Bvm1L*DRc*m1*m3*m4**2*PvR*sman**2 + 2*Bvm2L*DRc*m1*m3*m4**2*PvR*sman**2 - 2*BvgR*DRc*m1*m3*mq1*PvR*sman**2 + 
     &        2*Bsm1R*DLc*m1**2*mq3c*PvR*sman**2 + 2*BvgL*DLc*m1**2*mq3c*PvR*sman**2 + 2*Bvm1L*DLc*m1**2*m3*mq3c*PvR*sman**2 + 
     &        2*Bsm1R*DLc*m4**2*mq3c*PvR*sman**2 + 2*BvgL*DLc*m4**2*mq3c*PvR*sman**2 + 2*Bvm1L*DLc*m3*m4**2*mq3c*PvR*sman**2 + 2*Bvm1L*DRc*m1*m3*PvR*sman**3 + 
     &        2*Bvm1R*sman*(DLc*m4**4*mq1*mq3c*PvR + (DRc*m1*m3**2 + DLc*(m1**2 + m4**2)*mq3c)*(PsL - mq1*PvR)*sman) - 
     &        ((m3**2 + m4**2)*((Bsm1R*DRc*m1*m3**2 - Bsm2R*DRc*m1*m4**2 + 2*Bsm1R*DLc*m4**2*mq3c)*PvR + Bvm2R*DRc*m1*m4**2*(-PsL + mq1*PvR)) + 
     &           Bvm2R*DRc*m1*m4**2*PsL*sman)*tman + Bvm1R*(-((m3**2 + m4**2)*(DRc*m1*m3**2 + 2*DLc*m4**2*mq3c)*(PsL - mq1*PvR)) + 
     &           (DRc*m1*(2*m3**2 + m4**2) + 2*DLc*m4**2*mq3c)*PsL*sman)*tman)*xv - 
     &     (2*Bsm1L*DRc*m1*m3**3*PsL*tman + Bvm1R*DRc*m1*m3**4*PsL*tman - 2*Bsm2L*DRc*m1*m3*m4**2*PsL*tman - Bvm1R*DRc*m1*m3**2*m4**2*PsL*tman - 
     &        Bvm2R*DRc*m1*m3**2*m4**2*PsL*tman + Bvm2R*DRc*m1*m4**4*PsL*tman + 2*BvgR*DRc*m1*m3*(m3 - m4)*(m3 + m4)*PsL*tman - 
     &        2*Bsm1R*DRc*m1*m3**2*PvR*sman*tman - Bsm1R*DRc*m1*m4**2*PvR*sman*tman + Bsm2R*DRc*m1*m4**2*PvR*sman*tman + 
     &        2*Bvm1R*DRc*m1*m3**2*mq1*PvR*sman*tman + Bvm1R*DRc*m1*m4**2*mq1*PvR*sman*tman - Bvm2R*DRc*m1*m4**2*mq1*PvR*sman*tman - 
     &        2*Bsm1R*DLc*m4**2*mq3c*PvR*sman*tman + 2*Bvm1R*DLc*m4**2*mq1*mq3c*PvR*sman*tman + Bvm1R*DRc*m1*PsL*sman**2*tman + 
     &        Bsm1R*DRc*m1*PvR*sman**2*tman - Bvm1R*DRc*m1*mq1*PvR*sman**2*tman + 
     &        2*(-(Bsm1L*DLc*m3**3*mq3c*PsL) + Bsm2L*DLc*m3*m4**2*mq3c*PsL - Bvm2R*DLc*m4**4*mq3c*PsL + BvgL*DLc*m3**2*mq1*mq3c*PsL + 
     &           Bvm1L*DLc*m3**3*mq1*mq3c*PsL + BvgL*DLc*m4**2*mq1*mq3c*PsL - Bvm2L*DLc*m3*m4**2*mq1*mq3c*PsL + Bsm1R*DRc*m1*m3**2*m4**2*PvR + 
     &           Bsm1R*DRc*m1*m4**4*PvR + Bsm1R*DLc*m3**2*m4**2*mq3c*PvR - Bsm2R*DLc*m4**4*mq3c*PvR + Bvm2R*DLc*m4**4*mq1*mq3c*PvR + 
     &           Bvm1R*(PsL - mq1*PvR)*(DLc*mq3c*(m3**2 - sman)*(m4**2 - sman) + DRc*m1*m4**2*(m3**2 + m4**2 - sman)) - 
     &           BsgL*DLc*mq3c*PsL*(m3**2 + m4**2 - sman) + (Bsm1L*DLc*m3*mq3c*PsL - DLc*(BvgL + Bvm1L*m3)*mq1*mq3c*PsL - 
     &              (Bsm1R*DRc*m1*m4**2 + DLc*(m3**2*(Bsm1R + BvgL + Bvm1L*m3) + (Bsm1R - Bsm2R + BvgL - Bvm2L*m3)*m4**2)*mq3c)*PvR + 
     &              Bvm2R*DLc*m4**2*mq3c*(PsL - mq1*PvR))*sman + DLc*(Bsm1R + BvgL + Bvm1L*m3)*mq3c*PvR*sman**2)*uman)*xv - 
     &     2*(Bsm2R*DRc*m1*m3**2*PvR*tman - 2*BvgR*DRc*m1*m3*mq1*PvR*tman + Bvm1R*DRc*m1*m3**2*mq1*PvR*tman - Bvm2R*DRc*m1*m3**2*mq1*PvR*tman + 
     &        2*Bvm1R*DLc*m3**2*mq1*mq3c*PvR*tman + Bvm1L*DRc*m1*m3*PvR*sman*tman + Bvm2L*DRc*m1*m3*PvR*sman*tman - 2*BvgL*DLc*mq3c*PvR*sman*tman - 
     &        2*Bvm1L*DLc*m3*mq3c*PvR*sman*tman + 2*Bsm1L*DRc*m1*m3*PsL*uman + 2*BvgR*DRc*m1*m3*PsL*uman - 2*Bvm1R*DRc*m1*m4**2*PsL*uman - 
     &        2*Bvm1L*DRc*m1*m3*mq1*PsL*uman - Bsm1L*DLc*m3*mq3c*PsL*uman - Bsm2L*DLc*m3*mq3c*PsL*uman - 2*Bvm1R*DLc*m3**2*mq3c*PsL*uman + 
     &        Bvm1R*DLc*m4**2*mq3c*PsL*uman + Bvm2R*DLc*m4**2*mq3c*PsL*uman + Bvm1L*DLc*m3*mq1*mq3c*PsL*uman + Bvm2L*DLc*m3*mq1*mq3c*PsL*uman - 
     &        2*BvgR*DRc*m1*m3*mq1*PvR*uman + 2*Bvm1R*DRc*m1*m4**2*mq1*PvR*uman + Bsm2R*DLc*m4**2*mq3c*PvR*uman + 2*Bvm1R*DLc*m3**2*mq1*mq3c*PvR*uman - 
     &        Bvm1R*DLc*m4**2*mq1*mq3c*PvR*uman - Bvm2R*DLc*m4**2*mq1*mq3c*PvR*uman + 2*Bvm1R*DRc*m1*PsL*sman*uman - Bvm1R*DLc*mq3c*PsL*sman*uman - 
     &        Bvm2R*DLc*mq3c*PsL*sman*uman + 2*Bvm1L*DRc*m1*m3*PvR*sman*uman - 2*Bvm1R*DRc*m1*mq1*PvR*sman*uman - Bsm2R*DLc*mq3c*PvR*sman*uman - 
     &        Bvm1L*DLc*m3*mq3c*PvR*sman*uman - Bvm2L*DLc*m3*mq3c*PvR*sman*uman + Bvm1R*DLc*mq1*mq3c*PvR*sman*uman + Bvm2R*DLc*mq1*mq3c*PvR*sman*uman + 
     &        2*BsgR*DRc*m1*m3*PvR*(tman + uman) - Bsm1R*PvR*(DRc*m1*(m3**2*tman + 2*(m4**2 - sman)*uman) + 
     &           DLc*mq3c*((-m4**2 + sman)*uman + 2*m3**2*(tman + uman))) - Bsm1L*DRc*m1*m2**2*m3**3*PsL*xv) + 
     &     2*(Bvm1R*DRc*m1*m2**2*m4**2*PsL + Bvm2R*DRc*m1*m2**2*m4**2*PsL - Bvm1R*DRc*m1*m3**2*m4**2*PsL - Bvm2R*DRc*m1*m3**2*m4**2*PsL - 
     &        4*BvgL*DRc*m1*m2**2*mq1*PsL + 4*BvgL*DRc*m1*m3**2*mq1*PsL + Bvm1R*DLc*m1**2*m4**2*mq3c*PsL + Bvm2R*DLc*m1**2*m4**2*mq3c*PsL + 
     &        Bvm1R*DLc*m3**2*m4**2*mq3c*PsL + Bvm2R*DLc*m3**2*m4**2*mq3c*PsL - 4*BvgL*DLc*m1**2*mq1*mq3c*PsL - 4*BvgL*DLc*m3**2*mq1*mq3c*PsL - 
     &        8*BsgR*DRc*m1*m2**2*m3*PvR + Bsm1R*DRc*m1*m2**2*m4**2*PvR + Bsm2R*DRc*m1*m2**2*m4**2*PvR - 8*BsgR*DRc*m1*m3*m4**2*PvR - 
     &        Bsm1R*DRc*m1*m3**2*m4**2*PvR - Bsm2R*DRc*m1*m3**2*m4**2*PvR - Bvm1R*DRc*m1*m2**2*m4**2*mq1*PvR - Bvm2R*DRc*m1*m2**2*m4**2*mq1*PvR + 
     &        Bvm1R*DRc*m1*m3**2*m4**2*mq1*PvR + Bvm2R*DRc*m1*m3**2*m4**2*mq1*PvR - 8*BsgR*DLc*m1**2*m3*mq3c*PvR - 4*BsgR*DLc*m3**3*mq3c*PvR + 
     &        Bsm1R*DLc*m1**2*m4**2*mq3c*PvR + Bsm2R*DLc*m1**2*m4**2*mq3c*PvR - 4*BsgR*DLc*m3*m4**2*mq3c*PvR + Bsm1R*DLc*m3**2*m4**2*mq3c*PvR + 
     &        Bsm2R*DLc*m3**2*m4**2*mq3c*PvR - Bvm1R*DLc*m1**2*m4**2*mq1*mq3c*PvR - Bvm2R*DLc*m1**2*m4**2*mq1*mq3c*PvR - Bvm1R*DLc*m3**2*m4**2*mq1*mq3c*PvR - 
     &        Bvm2R*DLc*m3**2*m4**2*mq1*mq3c*PvR + 4*BvgL*DRc*m1*m2**2*PvR*sman + 4*BsgR*DRc*m1*m3*PvR*sman - 4*BvgL*DRc*m1*m3**2*PvR*sman + 
     &        4*BvgL*DLc*m1**2*mq3c*PvR*sman + 4*BvgL*DLc*m3**2*mq3c*PvR*sman - Bvm1R*DLc*m4**2*mq3c*PsL*tman - Bvm2R*DLc*m4**2*mq3c*PsL*tman + 
     &        4*BvgL*DLc*mq1*mq3c*PsL*tman + 4*BsgR*DRc*m1*m3*PvR*tman + 4*BsgR*DLc*m3*mq3c*PvR*tman - Bsm1R*DLc*m4**2*mq3c*PvR*tman - 
     &        Bsm2R*DLc*m4**2*mq3c*PvR*tman + Bvm1R*DLc*m4**2*mq1*mq3c*PvR*tman + Bvm2R*DLc*m4**2*mq1*mq3c*PvR*tman - 4*BvgL*DLc*mq3c*PvR*sman*tman - 
     &        Bvm1R*DRc*m1*m4**2*PsL*uman - Bvm2R*DRc*m1*m4**2*PsL*uman + 4*BvgL*DRc*m1*mq1*PsL*uman + 4*BsgR*DRc*m1*m3*PvR*uman - 
     &        Bsm1R*DRc*m1*m4**2*PvR*uman - Bsm2R*DRc*m1*m4**2*PvR*uman + Bvm1R*DRc*m1*m4**2*mq1*PvR*uman + Bvm2R*DRc*m1*m4**2*mq1*PvR*uman + 
     &        4*BsgR*DLc*m3*mq3c*PvR*uman - 4*BvgL*DRc*m1*PvR*sman*uman - 
     &        m4**2*(-(BvgL*DRc*m1*m2**2*mq1*PsL) + BvgL*DRc*m1*m3**2*mq1*PsL - BvgL*DLc*m1**2*mq1*mq3c*PsL - BvgL*DLc*m3**2*mq1*mq3c*PsL - 
     &           2*BsgR*DRc*m1*m2**2*m3*PvR + Bsm1R*DRc*m1*m2**2*m4**2*PvR + Bsm2R*DRc*m1*m2**2*m4**2*PvR - 2*BsgR*DRc*m1*m3*m4**2*PvR - 
     &           Bsm1R*DRc*m1*m3**2*m4**2*PvR - Bsm2R*DRc*m1*m3**2*m4**2*PvR - 2*BsgR*DLc*m1**2*m3*mq3c*PvR - BsgR*DLc*m3**3*mq3c*PvR + 
     &           Bsm1R*DLc*m1**2*m4**2*mq3c*PvR + Bsm2R*DLc*m1**2*m4**2*mq3c*PvR - BsgR*DLc*m3*m4**2*mq3c*PvR + Bsm1R*DLc*m3**2*m4**2*mq3c*PvR + 
     &           Bsm2R*DLc*m3**2*m4**2*mq3c*PvR + BvgL*DRc*m1*m2**2*PvR*sman + BsgR*DRc*m1*m3*PvR*sman - BvgL*DRc*m1*m3**2*PvR*sman + 
     &           BvgL*DLc*m1**2*mq3c*PvR*sman + BvgL*DLc*m3**2*mq3c*PvR*sman + BvgL*DLc*mq1*mq3c*PsL*tman + BsgR*DRc*m1*m3*PvR*tman + 
     &           BsgR*DLc*m3*mq3c*PvR*tman - Bsm1R*DLc*m4**2*mq3c*PvR*tman - Bsm2R*DLc*m4**2*mq3c*PvR*tman - BvgL*DLc*mq3c*PvR*sman*tman + 
     &           Bvm1R*m4**2*(PsL - mq1*PvR)*(DLc*mq3c*(m1**2 + m3**2 - tman) + DRc*m1*(m2**2 - m3**2 - uman)) + 
     &           Bvm2R*m4**2*(PsL - mq1*PvR)*(DLc*mq3c*(m1**2 + m3**2 - tman) + DRc*m1*(m2**2 - m3**2 - uman)) + 
     &           ((-((Bsm1R + Bsm2R)*DRc*m1*m4**2) + BsgR*m3*(DRc*m1 + DLc*mq3c))*PvR + BvgL*DRc*m1*(mq1*PsL - PvR*sman))*uman)*xv - 
     &        BsgL*PsL*(DLc*mq3c*(m1**2 + m3**2 - tman) + DRc*m1*(m2**2 - m3**2 - uman))*(-4 + m4**2*xv) + 
     &        BvgR*m3*(PsL - mq1*PvR)*(DLc*mq3c*(2*m1**2 + m3**2 + m4**2 - tman - uman) + DRc*m1*(2*m2**2 + 2*m4**2 - sman - tman - uman))*(-4 + m4**2*xv))) + 
     &  A1vR*A3Rc*(4*Bsm2L*DRc*m1**2*m2**2*m3*PsL + 8*Bsm2L*DRc*m1**2*m3*m4**2*PsL + 4*Bsm2L*DRc*m2**2*m3*m4**2*PsL + 4*Bsm2L*DRc*m3**3*m4**2*PsL + 
     &     4*Bsm2L*DRc*m3*m4**4*PsL - 8*BvgL*DRc*m1**2*m2**2*mq1*PsL - 12*Bvm1L*DRc*m1**2*m2**2*m3*mq1*PsL - 4*Bvm2L*DRc*m1**2*m2**2*m3*mq1*PsL + 
     &     8*BvgL*DRc*m1**2*m3**2*mq1*PsL - 4*BvgL*DRc*m2**2*m3**2*mq1*PsL - 8*Bvm1L*DRc*m2**2*m3**3*mq1*PsL + 4*BvgL*DRc*m3**4*mq1*PsL - 
     &     4*BvgL*DRc*m2**2*m4**2*mq1*PsL - 8*Bvm1L*DRc*m1**2*m3*m4**2*mq1*PsL - 8*Bvm2L*DRc*m1**2*m3*m4**2*mq1*PsL - 4*Bvm1L*DRc*m2**2*m3*m4**2*mq1*PsL - 
     &     4*Bvm2L*DRc*m2**2*m3*m4**2*mq1*PsL + 4*BvgL*DRc*m3**2*m4**2*mq1*PsL - 4*Bvm1L*DRc*m3**3*m4**2*mq1*PsL - 4*Bvm2L*DRc*m3**3*m4**2*mq1*PsL - 
     &     4*Bvm1L*DRc*m3*m4**4*mq1*PsL - 4*Bvm2L*DRc*m3*m4**4*mq1*PsL + 2*Bsm2L*DLc*m1*m3**3*mq3c*PsL - 2*Bsm2L*DLc*m1*m3*m4**2*mq3c*PsL - 
     &     4*BvgL*DLc*m1*m3**2*mq1*mq3c*PsL + 2*Bvm1L*DLc*m1*m3**3*mq1*mq3c*PsL - 2*Bvm2L*DLc*m1*m3**3*mq1*mq3c*PsL + 4*BvgL*DLc*m1*m4**2*mq1*mq3c*PsL - 
     &     2*Bvm1L*DLc*m1*m3*m4**2*mq1*mq3c*PsL + 2*Bvm2L*DLc*m1*m3*m4**2*mq1*mq3c*PsL + 6*Bvm1R*DRc*m1**2*m2**2*PsL*sman + 2*Bvm2R*DRc*m1**2*m2**2*PsL*sman - 
     &     Bsm2L*DRc*m1**2*m3*PsL*sman + 8*BvgR*DRc*m1**2*m3*PsL*sman + Bsm2L*DRc*m2**2*m3*PsL*sman + 4*BvgR*DRc*m2**2*m3*PsL*sman + 
     &     2*Bvm1R*DRc*m1**2*m3**2*PsL*sman - 2*Bvm2R*DRc*m1**2*m3**2*PsL*sman + 8*Bvm1R*DRc*m2**2*m3**2*PsL*sman + Bsm2L*DRc*m3**3*PsL*sman + 
     &     4*BvgR*DRc*m3**3*PsL*sman - 3*Bsm2L*DRc*m3*m4**2*PsL*sman - Bvm1L*DRc*m1**2*m3*mq1*PsL*sman + Bvm2L*DRc*m1**2*m3*mq1*PsL*sman - 
     &     3*Bvm1L*DRc*m2**2*m3*mq1*PsL*sman - Bvm2L*DRc*m2**2*m3*mq1*PsL*sman + Bvm1L*DRc*m3**3*mq1*PsL*sman - Bvm2L*DRc*m3**3*mq1*PsL*sman + 
     &     Bvm1L*DRc*m3*m4**2*mq1*PsL*sman + 3*Bvm2L*DRc*m3*m4**2*mq1*PsL*sman - 2*Bsm2L*DLc*m1*m3*mq3c*PsL*sman - 16*BvgR*DLc*m1*m3*mq3c*PsL*sman - 
     &     6*Bvm1R*DLc*m1*m3**2*mq3c*PsL*sman + 2*Bvm2R*DLc*m1*m3**2*mq3c*PsL*sman + 2*Bvm1R*DLc*m1*m4**2*mq3c*PsL*sman + 2*Bvm2R*DLc*m1*m4**2*mq3c*PsL*sman - 
     &     4*BvgL*DLc*m1*mq1*mq3c*PsL*sman + 6*Bvm1L*DLc*m1*m3*mq1*mq3c*PsL*sman + 2*Bvm2L*DLc*m1*m3*mq1*mq3c*PsL*sman + 6*Bsm1R*DRc*m1**2*m2**2*PvR*sman + 
     &     2*Bsm2R*DRc*m1**2*m2**2*PvR*sman + 8*BvgL*DRc*m1**2*m2**2*PvR*sman + 8*BsgR*DRc*m1**2*m3*PvR*sman + 4*BsgR*DRc*m2**2*m3*PvR*sman + 
     &     12*Bvm1L*DRc*m1**2*m2**2*m3*PvR*sman + 4*Bvm2L*DRc*m1**2*m2**2*m3*PvR*sman + 2*Bsm1R*DRc*m1**2*m3**2*PvR*sman - 2*Bsm2R*DRc*m1**2*m3**2*PvR*sman - 
     &     8*BvgL*DRc*m1**2*m3**2*PvR*sman + 8*Bsm1R*DRc*m2**2*m3**2*PvR*sman + 4*BvgL*DRc*m2**2*m3**2*PvR*sman + 4*BsgR*DRc*m3**3*PvR*sman + 
     &     8*Bvm1L*DRc*m2**2*m3**3*PvR*sman - 4*BvgL*DRc*m3**4*PvR*sman + 4*BvgL*DRc*m2**2*m4**2*PvR*sman + 8*Bvm1L*DRc*m1**2*m3*m4**2*PvR*sman + 
     &     8*Bvm2L*DRc*m1**2*m3*m4**2*PvR*sman + 4*Bvm1L*DRc*m2**2*m3*m4**2*PvR*sman + 4*Bvm2L*DRc*m2**2*m3*m4**2*PvR*sman - 4*BvgL*DRc*m3**2*m4**2*PvR*sman + 
     &     4*Bvm1L*DRc*m3**3*m4**2*PvR*sman + 4*Bvm2L*DRc*m3**3*m4**2*PvR*sman + 4*Bvm1L*DRc*m3*m4**4*PvR*sman + 4*Bvm2L*DRc*m3*m4**4*PvR*sman - 
     &     6*Bvm1R*DRc*m1**2*m2**2*mq1*PvR*sman - 2*Bvm2R*DRc*m1**2*m2**2*mq1*PvR*sman - 8*BvgR*DRc*m1**2*m3*mq1*PvR*sman - 4*BvgR*DRc*m2**2*m3*mq1*PvR*sman - 
     &     2*Bvm1R*DRc*m1**2*m3**2*mq1*PvR*sman + 2*Bvm2R*DRc*m1**2*m3**2*mq1*PvR*sman - 8*Bvm1R*DRc*m2**2*m3**2*mq1*PvR*sman - 
     &     4*BvgR*DRc*m3**3*mq1*PvR*sman - 16*BsgR*DLc*m1*m3*mq3c*PvR*sman - 6*Bsm1R*DLc*m1*m3**2*mq3c*PvR*sman + 2*Bsm2R*DLc*m1*m3**2*mq3c*PvR*sman + 
     &     4*BvgL*DLc*m1*m3**2*mq3c*PvR*sman - 2*Bvm1L*DLc*m1*m3**3*mq3c*PvR*sman + 2*Bvm2L*DLc*m1*m3**3*mq3c*PvR*sman + 2*Bsm1R*DLc*m1*m4**2*mq3c*PvR*sman + 
     &     2*Bsm2R*DLc*m1*m4**2*mq3c*PvR*sman - 4*BvgL*DLc*m1*m4**2*mq3c*PvR*sman + 2*Bvm1L*DLc*m1*m3*m4**2*mq3c*PvR*sman - 
     &     2*Bvm2L*DLc*m1*m3*m4**2*mq3c*PvR*sman + 16*BvgR*DLc*m1*m3*mq1*mq3c*PvR*sman + 6*Bvm1R*DLc*m1*m3**2*mq1*mq3c*PvR*sman - 
     &     2*Bvm2R*DLc*m1*m3**2*mq1*mq3c*PvR*sman - 2*Bvm1R*DLc*m1*m4**2*mq1*mq3c*PvR*sman - 2*Bvm2R*DLc*m1*m4**2*mq1*mq3c*PvR*sman + 
     &     Bvm1R*DRc*m1**2*PsL*sman**2 + Bvm2R*DRc*m1**2*PsL*sman**2 + Bvm1R*DRc*m2**2*PsL*sman**2 + Bvm2R*DRc*m2**2*PsL*sman**2 - Bsm2L*DRc*m3*PsL*sman**2 - 
     &     4*BvgR*DRc*m3*PsL*sman**2 - 3*Bvm1R*DRc*m3**2*PsL*sman**2 + Bvm2R*DRc*m3**2*PsL*sman**2 + Bvm1R*DRc*m4**2*PsL*sman**2 + 
     &     Bvm2R*DRc*m4**2*PsL*sman**2 + 3*Bvm1L*DRc*m3*mq1*PsL*sman**2 + Bvm2L*DRc*m3*mq1*PsL*sman**2 - 2*Bvm1R*DLc*m1*mq3c*PsL*sman**2 - 
     &     2*Bvm2R*DLc*m1*mq3c*PsL*sman**2 + Bsm1R*DRc*m1**2*PvR*sman**2 + Bsm2R*DRc*m1**2*PvR*sman**2 + Bsm1R*DRc*m2**2*PvR*sman**2 + 
     &     Bsm2R*DRc*m2**2*PvR*sman**2 - 4*BsgR*DRc*m3*PvR*sman**2 + Bvm1L*DRc*m1**2*m3*PvR*sman**2 - Bvm2L*DRc*m1**2*m3*PvR*sman**2 + 
     &     3*Bvm1L*DRc*m2**2*m3*PvR*sman**2 + Bvm2L*DRc*m2**2*m3*PvR*sman**2 - 3*Bsm1R*DRc*m3**2*PvR*sman**2 + Bsm2R*DRc*m3**2*PvR*sman**2 - 
     &     Bvm1L*DRc*m3**3*PvR*sman**2 + Bvm2L*DRc*m3**3*PvR*sman**2 + Bsm1R*DRc*m4**2*PvR*sman**2 + Bsm2R*DRc*m4**2*PvR*sman**2 - 
     &     Bvm1L*DRc*m3*m4**2*PvR*sman**2 - 3*Bvm2L*DRc*m3*m4**2*PvR*sman**2 - Bvm1R*DRc*m1**2*mq1*PvR*sman**2 - Bvm2R*DRc*m1**2*mq1*PvR*sman**2 - 
     &     Bvm1R*DRc*m2**2*mq1*PvR*sman**2 - Bvm2R*DRc*m2**2*mq1*PvR*sman**2 + 4*BvgR*DRc*m3*mq1*PvR*sman**2 + 3*Bvm1R*DRc*m3**2*mq1*PvR*sman**2 - 
     &     Bvm2R*DRc*m3**2*mq1*PvR*sman**2 - Bvm1R*DRc*m4**2*mq1*PvR*sman**2 - Bvm2R*DRc*m4**2*mq1*PvR*sman**2 - 2*Bsm1R*DLc*m1*mq3c*PvR*sman**2 - 
     &     2*Bsm2R*DLc*m1*mq3c*PvR*sman**2 + 4*BvgL*DLc*m1*mq3c*PvR*sman**2 - 6*Bvm1L*DLc*m1*m3*mq3c*PvR*sman**2 - 2*Bvm2L*DLc*m1*m3*mq3c*PvR*sman**2 + 
     &     2*Bvm1R*DLc*m1*mq1*mq3c*PvR*sman**2 + 2*Bvm2R*DLc*m1*mq1*mq3c*PvR*sman**2 - Bvm1R*DRc*PsL*sman**3 - Bvm2R*DRc*PsL*sman**3 - Bsm1R*DRc*PvR*sman**3 - 
     &     Bsm2R*DRc*PvR*sman**3 - 3*Bvm1L*DRc*m3*PvR*sman**3 - Bvm2L*DRc*m3*PvR*sman**3 + Bvm1R*DRc*mq1*PvR*sman**3 + Bvm2R*DRc*mq1*PvR*sman**3 - 
     &     2*Bsm2L*DRc*m3*m4**2*PsL*tman - 4*BvgL*DRc*m1**2*mq1*PsL*tman - 4*BvgL*DRc*m2**2*mq1*PsL*tman - 12*BvgL*DRc*m3**2*mq1*PsL*tman - 
     &     4*BvgL*DRc*m4**2*mq1*PsL*tman + 2*Bvm1L*DRc*m3*m4**2*mq1*PsL*tman + 2*Bvm2L*DRc*m3*m4**2*mq1*PsL*tman + 
     &     DRc*(-(Bsm2L*m3*(3*m1**2 + m2**2 + m3**2 + 3*m4**2)) + 
     &        (4*BvgL*(m1**2 + 2*(m2**2 + m3**2) + m4**2) + m3*(Bvm2L*(3*m1**2 + m2**2 + m3**2 + 3*m4**2) + Bvm1L*(5*m1**2 + 7*m2**2 + 3*m3**2 + 5*m4**2)))*mq1
     &        )*PsL*tman - 8*BvgR*DRc*m3*PsL*sman*tman + 4*BvgL*DRc*m1**2*PvR*sman*tman + 4*BvgL*DRc*m2**2*PvR*sman*tman - 8*BsgR*DRc*m3*PvR*sman*tman + 
     &     12*BvgL*DRc*m3**2*PvR*sman*tman + 4*BvgL*DRc*m4**2*PvR*sman*tman - 2*Bvm1L*DRc*m3*m4**2*PvR*sman*tman - 2*Bvm2L*DRc*m3*m4**2*PvR*sman*tman + 
     &     8*BvgR*DRc*m3*mq1*PvR*sman*tman + 4*BvgL*DRc*mq1*PsL*tman**2 - 4*BvgL*DRc*PvR*sman*tman**2 + 
     &     Bsm1L*m3*PsL*(-2*DLc*m1*mq3c*(m3**2 - 3*m4**2 + 3*sman) + 
     &        DRc*(2*m3**2*m4**2 + 4*m4**4 - m3**2*sman + m4**2*sman - 3*sman**2 + m2**2*(8*m3**2 + 2*m4**2 + 3*sman - 7*tman) + 
     &           m1**2*(12*m2**2 + 4*m4**2 + sman - 5*tman) - 3*m3**2*tman - 5*m4**2*tman)) + 
     &     4*BsgL*PsL*(-(DLc*m1*mq3c*(m3**2 - m4**2 + sman)) + DRc*
     &         (m3**4 + m3**2*m4**2 - m3**2*sman + m4**2*sman - sman**2 + m2**2*(3*m3**2 - m4**2 + sman - 2*tman) + m1**2*(2*m2**2 + 2*m3**2 + sman - tman) - 
     &           (2*m3**2 + m4**2)*tman)) + 4*BvgL*DRc*m1**2*mq1*PsL*uman + 4*BvgL*DRc*m2**2*mq1*PsL*uman - 4*BvgL*DRc*m3**2*mq1*PsL*uman + 
     &     4*BvgL*DRc*m4**2*mq1*PsL*uman - 4*BvgL*DRc*m1**2*PvR*sman*uman - 4*BvgL*DRc*m2**2*PvR*sman*uman + 4*BvgL*DRc*m3**2*PvR*sman*uman - 
     &     4*BvgL*DRc*m4**2*PvR*sman*uman - 4*BvgL*DRc*mq1*PsL*uman**2 + 4*BvgL*DRc*PvR*sman*uman**2 + 
     &     (4*BsgL*DRc*m1**2*m3**2*m4**2*PsL + 2*BsgL*DRc*m2**2*m3**2*m4**2*PsL + 2*BsgL*DRc*m3**4*m4**2*PsL + 2*BsgL*DRc*m2**2*m4**4*PsL + 
     &        2*BsgL*DRc*m3**2*m4**4*PsL - 4*BvgL*DRc*m1**2*m2**2*m3**2*mq1*PsL - 4*Bvm1L*DRc*m1**2*m2**2*m3**3*mq1*PsL + 
     &        4*Bvm2L*DRc*m1**2*m2**2*m3*m4**2*mq1*PsL - 4*BvgL*DRc*m1**2*m3**2*m4**2*mq1*PsL - 2*BvgL*DRc*m2**2*m3**2*m4**2*mq1*PsL - 
     &        8*Bvm1L*DRc*m1**2*m3**3*m4**2*mq1*PsL - 4*Bvm1L*DRc*m2**2*m3**3*m4**2*mq1*PsL - 2*BvgL*DRc*m3**4*m4**2*mq1*PsL - 
     &        4*Bvm1L*DRc*m3**5*m4**2*mq1*PsL - 2*BvgL*DRc*m2**2*m4**4*mq1*PsL - 4*Bvm1L*DRc*m1**2*m3*m4**4*mq1*PsL + 4*Bvm2L*DRc*m1**2*m3*m4**4*mq1*PsL - 
     &        2*Bvm1L*DRc*m2**2*m3*m4**4*mq1*PsL + 2*Bvm2L*DRc*m2**2*m3*m4**4*mq1*PsL - 2*BvgL*DRc*m3**2*m4**4*mq1*PsL - 6*Bvm1L*DRc*m3**3*m4**4*mq1*PsL + 
     &        2*Bvm2L*DRc*m3**3*m4**4*mq1*PsL + 4*Bvm2L*DRc*m3*m4**6*mq1*PsL + 2*BsgL*DLc*m1*m3**4*mq3c*PsL - 2*BsgL*DLc*m1*m4**4*mq3c*PsL - 
     &        2*BvgL*DLc*m1*m3**4*mq1*mq3c*PsL - 2*Bvm1L*DLc*m1*m3**5*mq1*mq3c*PsL + 2*Bvm1L*DLc*m1*m3**3*m4**2*mq1*mq3c*PsL + 
     &        2*Bvm2L*DLc*m1*m3**3*m4**2*mq1*mq3c*PsL + 2*BvgL*DLc*m1*m4**4*mq1*mq3c*PsL + 4*Bvm1L*DLc*m1*m3*m4**4*mq1*mq3c*PsL + 
     &        2*Bvm2L*DLc*m1*m3*m4**4*mq1*mq3c*PsL + (2*Bvm1R*DRc*m1**2*m2**2*m3**2*PsL - 2*Bvm1R*DRc*m1**2*m3**4*PsL - 2*Bvm2R*DRc*m1**2*m2**2*m4**2*PsL + 
     &           2*Bvm2R*DRc*m1**2*m3**2*m4**2*PsL + 2*BvgL*DRc*m1**2*m2**2*mq1*PsL + 4*Bvm1L*DRc*m1**2*m2**2*m3*mq1*PsL - BvgL*DRc*m1**2*m3**2*mq1*PsL - 
     &           BvgL*DRc*m2**2*m3**2*mq1*PsL + Bvm1L*DRc*m1**2*m3**3*mq1*PsL - Bvm1L*DRc*m2**2*m3**3*mq1*PsL - BvgL*DRc*m3**4*mq1*PsL - 
     &           Bvm1L*DRc*m3**5*mq1*PsL + BvgL*DRc*m1**2*m4**2*mq1*PsL + BvgL*DRc*m2**2*m4**2*mq1*PsL + 8*Bvm1L*DRc*m1**2*m3*m4**2*mq1*PsL - 
     &           Bvm2L*DRc*m1**2*m3*m4**2*mq1*PsL + 4*Bvm1L*DRc*m2**2*m3*m4**2*mq1*PsL + Bvm2L*DRc*m2**2*m3*m4**2*mq1*PsL + 2*BvgL*DRc*m3**2*m4**2*mq1*PsL + 
     &           7*Bvm1L*DRc*m3**3*m4**2*mq1*PsL + Bvm2L*DRc*m3**3*m4**2*mq1*PsL + BvgL*DRc*m4**4*mq1*PsL + 6*Bvm1L*DRc*m3*m4**4*mq1*PsL - 
     &           Bvm2L*DRc*m3*m4**4*mq1*PsL + 2*Bvm1R*DLc*m1*m3**4*mq3c*PsL + 2*Bvm1R*DLc*m1*m3**2*m4**2*mq3c*PsL - 2*Bvm2R*DLc*m1*m3**2*m4**2*mq3c*PsL - 
     &           2*Bvm2R*DLc*m1*m4**4*mq3c*PsL + 4*BvgL*DLc*m1*m3**2*mq1*mq3c*PsL + 4*Bvm1L*DLc*m1*m3**3*mq1*mq3c*PsL - 2*Bvm1L*DLc*m1*m3*m4**2*mq1*mq3c*PsL - 
     &           2*Bvm2L*DLc*m1*m3*m4**2*mq1*mq3c*PsL - BsgL*(DRc*
     &               (-(m3**2*(m2**2 + m3**2)) + (m2**2 + 2*m3**2)*m4**2 + m4**4 + m1**2*(2*m2**2 - m3**2 + m4**2)) + 4*DLc*m1*m3**2*mq3c)*PsL - 
     &           2*(-(Bsm1R*DRc*m1**2*m2**2*m3**2) - 2*BvgL*DRc*m1**2*m2**2*m3**2 - 2*Bvm1L*DRc*m1**2*m2**2*m3**3 + Bsm1R*DRc*m1**2*m3**4 + 
     &              Bsm2R*DRc*m1**2*m2**2*m4**2 + 2*Bvm2L*DRc*m1**2*m2**2*m3*m4**2 - Bsm2R*DRc*m1**2*m3**2*m4**2 - 2*BvgL*DRc*m1**2*m3**2*m4**2 - 
     &              BvgL*DRc*m2**2*m3**2*m4**2 - 4*Bvm1L*DRc*m1**2*m3**3*m4**2 - 2*Bvm1L*DRc*m2**2*m3**3*m4**2 - BvgL*DRc*m3**4*m4**2 - 
     &              2*Bvm1L*DRc*m3**5*m4**2 - BvgL*DRc*m2**2*m4**4 - 2*Bvm1L*DRc*m1**2*m3*m4**4 + 2*Bvm2L*DRc*m1**2*m3*m4**4 - Bvm1L*DRc*m2**2*m3*m4**4 + 
     &              Bvm2L*DRc*m2**2*m3*m4**4 - BvgL*DRc*m3**2*m4**4 - 3*Bvm1L*DRc*m3**3*m4**4 + Bvm2L*DRc*m3**3*m4**4 + 2*Bvm2L*DRc*m3*m4**6 + 
     &              BsgR*DRc*m3*(m1**2 + m4**2)*(m2**2 + m3**2 + 2*m4**2) + Bvm1R*DRc*m1**2*m2**2*m3**2*mq1 - Bvm1R*DRc*m1**2*m3**4*mq1 - 
     &              Bvm2R*DRc*m1**2*m2**2*m4**2*mq1 + Bvm2R*DRc*m1**2*m3**2*m4**2*mq1 + 
     &              DLc*m1*(-(m3**4*(Bsm1R + BvgL + Bvm1L*m3)) + m3**2*(-Bsm1R + Bsm2R + (Bvm1L + Bvm2L)*m3)*m4**2 + 
     &                 (Bsm2R + BvgL + (2*Bvm1L + Bvm2L)*m3)*m4**4 + m3**2*(Bvm1R*m3**2 + (Bvm1R - Bvm2R)*m4**2)*mq1)*mq3c)*PvR - 
     &           2*BvgR*DRc*m3*(m1**2 + m4**2)*(m2**2 + m3**2 + 2*m4**2)*(PsL - mq1*PvR))*sman - 
     &        Bsm2L*m3*m4**2*PsL*(2*DLc*m1*mq3c*(m3**2 + m4**2 - sman) + 
     &           DRc*(2*m4**2*(m2**2 + m3**2 + 2*m4**2) + m1**2*(4*(m2**2 + m4**2) - sman) + (m2**2 + m3**2 - m4**2)*sman)) + 
     &        Bsm1L*m3*PsL*(2*DRc*(2*m1**2*m2**2*m3**2 + 2*m3**2*(2*m1**2 + m2**2 + m3**2)*m4**2 + (2*m1**2 + m2**2 + 3*m3**2)*m4**4) - 
     &           DRc*(-(m3**2*(m2**2 + m3**2)) + (4*m2**2 + 7*m3**2)*m4**2 + 6*m4**4 + m1**2*(4*m2**2 + m3**2 + 8*m4**2))*sman + 
     &           2*DLc*m1*mq3c*(m3**4 - 2*m4**4 + m4**2*sman - m3**2*(m4**2 + 2*sman))))*xv + 
     &     (-(BsgL*DRc*m1**2*PsL*sman**2) - BsgL*DRc*m2**2*PsL*sman**2 - 2*Bvm1R*DRc*m1**2*m2**2*PsL*sman**2 + Bsm1L*DRc*m1**2*m3*PsL*sman**2 + 
     &        2*BvgR*DRc*m1**2*m3*PsL*sman**2 - Bsm1L*DRc*m2**2*m3*PsL*sman**2 - 2*BsgL*DRc*m3**2*PsL*sman**2 + 3*Bvm1R*DRc*m1**2*m3**2*PsL*sman**2 + 
     &        Bvm1R*DRc*m2**2*m3**2*PsL*sman**2 - 2*Bsm1L*DRc*m3**3*PsL*sman**2 + Bvm1R*DRc*m3**4*PsL*sman**2 + 3*Bsm1L*DRc*m3*m4**2*PsL*sman**2 + 
     &        Bsm2L*DRc*m3*m4**2*PsL*sman**2 + 2*BvgR*DRc*m3*m4**2*PsL*sman**2 + Bvm1R*DRc*m3**2*m4**2*PsL*sman**2 + BvgL*DRc*m1**2*mq1*PsL*sman**2 + 
     &        BvgL*DRc*m2**2*mq1*PsL*sman**2 - Bvm1L*DRc*m1**2*m3*mq1*PsL*sman**2 + Bvm1L*DRc*m2**2*m3*mq1*PsL*sman**2 + 2*BvgL*DRc*m3**2*mq1*PsL*sman**2 + 
     &        2*Bvm1L*DRc*m3**3*mq1*PsL*sman**2 - 3*Bvm1L*DRc*m3*m4**2*mq1*PsL*sman**2 - Bvm2L*DRc*m3*m4**2*mq1*PsL*sman**2 + 2*BsgL*DLc*m1*mq3c*PsL*sman**2 + 
     &        2*Bsm1L*DLc*m1*m3*mq3c*PsL*sman**2 - 4*Bvm1R*DLc*m1*m3**2*mq3c*PsL*sman**2 - 2*Bvm1R*DLc*m1*m4**2*mq3c*PsL*sman**2 - 
     &        2*BvgL*DLc*m1*mq1*mq3c*PsL*sman**2 - 2*Bvm1L*DLc*m1*m3*mq1*mq3c*PsL*sman**2 - 2*Bsm1R*DRc*m1**2*m2**2*PvR*sman**2 - 
     &        2*BvgL*DRc*m1**2*m2**2*PvR*sman**2 + 2*BsgR*DRc*m1**2*m3*PvR*sman**2 - 4*Bvm1L*DRc*m1**2*m2**2*m3*PvR*sman**2 + 
     &        3*Bsm1R*DRc*m1**2*m3**2*PvR*sman**2 + BvgL*DRc*m1**2*m3**2*PvR*sman**2 + Bsm1R*DRc*m2**2*m3**2*PvR*sman**2 + BvgL*DRc*m2**2*m3**2*PvR*sman**2 - 
     &        Bvm1L*DRc*m1**2*m3**3*PvR*sman**2 + Bvm1L*DRc*m2**2*m3**3*PvR*sman**2 + Bsm1R*DRc*m3**4*PvR*sman**2 + BvgL*DRc*m3**4*PvR*sman**2 + 
     &        Bvm1L*DRc*m3**5*PvR*sman**2 - Bsm2R*DRc*m1**2*m4**2*PvR*sman**2 - BvgL*DRc*m1**2*m4**2*PvR*sman**2 - Bsm2R*DRc*m2**2*m4**2*PvR*sman**2 - 
     &        BvgL*DRc*m2**2*m4**2*PvR*sman**2 + 2*BsgR*DRc*m3*m4**2*PvR*sman**2 - 8*Bvm1L*DRc*m1**2*m3*m4**2*PvR*sman**2 + 
     &        Bvm2L*DRc*m1**2*m3*m4**2*PvR*sman**2 - 4*Bvm1L*DRc*m2**2*m3*m4**2*PvR*sman**2 - Bvm2L*DRc*m2**2*m3*m4**2*PvR*sman**2 + 
     &        Bsm1R*DRc*m3**2*m4**2*PvR*sman**2 - Bsm2R*DRc*m3**2*m4**2*PvR*sman**2 - 2*BvgL*DRc*m3**2*m4**2*PvR*sman**2 - 
     &        7*Bvm1L*DRc*m3**3*m4**2*PvR*sman**2 - Bvm2L*DRc*m3**3*m4**2*PvR*sman**2 - Bsm2R*DRc*m4**4*PvR*sman**2 - BvgL*DRc*m4**4*PvR*sman**2 - 
     &        6*Bvm1L*DRc*m3*m4**4*PvR*sman**2 + Bvm2L*DRc*m3*m4**4*PvR*sman**2 + 2*Bvm1R*DRc*m1**2*m2**2*mq1*PvR*sman**2 - 
     &        2*BvgR*DRc*m1**2*m3*mq1*PvR*sman**2 - 3*Bvm1R*DRc*m1**2*m3**2*mq1*PvR*sman**2 - Bvm1R*DRc*m2**2*m3**2*mq1*PvR*sman**2 - 
     &        Bvm1R*DRc*m3**4*mq1*PvR*sman**2 - 2*BvgR*DRc*m3*m4**2*mq1*PvR*sman**2 - Bvm1R*DRc*m3**2*m4**2*mq1*PvR*sman**2 - 
     &        4*Bsm1R*DLc*m1*m3**2*mq3c*PvR*sman**2 - 4*BvgL*DLc*m1*m3**2*mq3c*PvR*sman**2 - 4*Bvm1L*DLc*m1*m3**3*mq3c*PvR*sman**2 - 
     &        2*Bsm1R*DLc*m1*m4**2*mq3c*PvR*sman**2 + 2*Bsm2R*DLc*m1*m4**2*mq3c*PvR*sman**2 + 2*Bvm1L*DLc*m1*m3*m4**2*mq3c*PvR*sman**2 + 
     &        2*Bvm2L*DLc*m1*m3*m4**2*mq3c*PvR*sman**2 + 4*Bvm1R*DLc*m1*m3**2*mq1*mq3c*PvR*sman**2 + 2*Bvm1R*DLc*m1*m4**2*mq1*mq3c*PvR*sman**2 + 
     &        BsgL*DRc*PsL*sman**3 - Bvm1R*DRc*m1**2*PsL*sman**3 - Bvm1R*DRc*m2**2*PsL*sman**3 + Bsm1L*DRc*m3*PsL*sman**3 - 2*Bvm1R*DRc*m3**2*PsL*sman**3 - 
     &        Bvm1R*DRc*m4**2*PsL*sman**3 - BvgL*DRc*mq1*PsL*sman**3 - Bvm1L*DRc*m3*mq1*PsL*sman**3 + 2*Bvm1R*DLc*m1*mq3c*PsL*sman**3 - 
     &        Bsm1R*DRc*m1**2*PvR*sman**3 - BvgL*DRc*m1**2*PvR*sman**3 - Bsm1R*DRc*m2**2*PvR*sman**3 - BvgL*DRc*m2**2*PvR*sman**3 + 
     &        Bvm1L*DRc*m1**2*m3*PvR*sman**3 - Bvm1L*DRc*m2**2*m3*PvR*sman**3 - 2*Bsm1R*DRc*m3**2*PvR*sman**3 - 2*BvgL*DRc*m3**2*PvR*sman**3 - 
     &        2*Bvm1L*DRc*m3**3*PvR*sman**3 - Bsm1R*DRc*m4**2*PvR*sman**3 + Bsm2R*DRc*m4**2*PvR*sman**3 + 3*Bvm1L*DRc*m3*m4**2*PvR*sman**3 + 
     &        Bvm2L*DRc*m3*m4**2*PvR*sman**3 + Bvm1R*DRc*m1**2*mq1*PvR*sman**3 + Bvm1R*DRc*m2**2*mq1*PvR*sman**3 + 2*Bvm1R*DRc*m3**2*mq1*PvR*sman**3 + 
     &        Bvm1R*DRc*m4**2*mq1*PvR*sman**3 + 2*Bsm1R*DLc*m1*mq3c*PvR*sman**3 + 2*BvgL*DLc*m1*mq3c*PvR*sman**3 + 2*Bvm1L*DLc*m1*m3*mq3c*PvR*sman**3 - 
     &        2*Bvm1R*DLc*m1*mq1*mq3c*PvR*sman**3 + Bvm1R*DRc*PsL*sman**4 + Bsm1R*DRc*PvR*sman**4 + BvgL*DRc*PvR*sman**4 + Bvm1L*DRc*m3*PvR*sman**4 - 
     &        Bvm1R*DRc*mq1*PvR*sman**4 - DRc*PsL*(-(Bsm2L*m1**2*m3*m4**2) - Bsm2L*m2**2*m3*m4**2 - 2*Bsm2L*m3*m4**4 - BvgL*m1**2*m3**2*mq1 - 
     &           BvgL*m2**2*m3**2*mq1 - Bvm1L*m1**2*m3**3*mq1 - Bvm1L*m2**2*m3**3*mq1 - BvgL*m1**2*m4**2*mq1 + BvgL*m2**2*m4**2*mq1 + 
     &           Bvm2L*m1**2*m3*m4**2*mq1 + Bvm2L*m2**2*m3*m4**2*mq1 - 2*BvgL*m3**2*m4**2*mq1 - 4*Bvm1L*m3**3*m4**2*mq1 - 2*Bvm1L*m3*m4**4*mq1 + 
     &           2*Bvm2L*m3*m4**4*mq1 + Bsm1L*m3*(2*m4**2*(2*m3**2 + m4**2 - 2*sman) + m1**2*(m3**2 - sman) + m2**2*(m3**2 - sman)) + 
     &           BsgL*(2*m4**2*(m3**2 - sman) + m2**2*(m3**2 - m4**2 - sman) + m1**2*(m3**2 + m4**2 - sman)) + 
     &           (Bvm1R*m3**2*(m2**2 + m4**2) + ((m1**2 + m2**2)*(BvgL + Bvm1L*m3) + 2*(BvgL + 2*Bvm1L*m3)*m4**2)*mq1)*sman)*tman + 
     &        Bvm2R*m4**2*sman*(-(DRc*(PsL - mq1*PvR)*(m1**2 + m2**2 + m3**2 + m4**2 - sman)*sman) + 2*DLc*m1*mq3c*(mq1*PvR*(m4**2 - sman) + PsL*sman) + 
     &           DRc*(m2**2 + m4**2)*PsL*tman))*xv + DRc*tman*(-2*Bvm2L*m1**2*m3*m4**2*mq1*PsL - Bvm2L*m3**3*m4**2*mq1*PsL - Bvm2L*m3*m4**4*mq1*PsL + 
     &        2*Bsm1L*m1**2*m3*PsL*sman + 2*BvgR*m1**2*m3*PsL*sman + BsgL*m3**2*PsL*sman - Bvm1R*m1**2*m3**2*PsL*sman + Bsm1L*m3**3*PsL*sman - 
     &        Bvm1R*m3**4*PsL*sman - BsgL*m4**2*PsL*sman + Bvm2R*m1**2*m4**2*PsL*sman + Bsm1L*m3*m4**2*PsL*sman + 2*BvgR*m3*m4**2*PsL*sman + 
     &        Bvm2R*m3**2*m4**2*PsL*sman + 2*BsgR*m1**2*m3*PvR*sman - Bsm1R*m1**2*m3**2*PvR*sman - Bsm1R*m3**4*PvR*sman + Bsm2R*m1**2*m4**2*PvR*sman + 
     &        2*BsgR*m3*m4**2*PvR*sman + 2*Bvm2L*m1**2*m3*m4**2*PvR*sman + Bsm2R*m3**2*m4**2*PvR*sman + Bvm2L*m3**3*m4**2*PvR*sman + Bvm2L*m3*m4**4*PvR*sman - 
     &        2*BvgR*m1**2*m3*mq1*PvR*sman + Bvm1R*m1**2*m3**2*mq1*PvR*sman + Bvm1R*m3**4*mq1*PvR*sman - Bvm2R*m1**2*m4**2*mq1*PvR*sman - 
     &        2*BvgR*m3*m4**2*mq1*PvR*sman - Bvm2R*m3**2*m4**2*mq1*PvR*sman + Bvm1R*m1**2*PsL*sman**2 + Bvm1R*m3**2*PsL*sman**2 + Bsm1R*m1**2*PvR*sman**2 + 
     &        Bsm1R*m3**2*PvR*sman**2 - Bvm1R*m1**2*mq1*PvR*sman**2 - Bvm1R*m3**2*mq1*PvR*sman**2 + 
     &        Bvm1L*m3*(m3**2 - sman)*(mq1*PsL - PvR*sman)*(2*m1**2 + m3**2 + m4**2 - tman) + 
     &        (m3*m4**2*(-Bsm2L + Bvm2L*mq1)*PsL + Bsm1L*m3*PsL*(m3**2 - sman) + BsgL*PsL*(m3**2 - m4**2 - sman) + 
     &           (Bsm1R*m3**2*PvR + Bvm1R*m3**2*(PsL - mq1*PvR) - m4**2*(Bvm2R*PsL + (Bsm2R + Bvm2L*m3 - Bvm2R*mq1)*PvR))*sman - 
     &           (Bsm1R*PvR + Bvm1R*(PsL - mq1*PvR))*sman**2)*tman + 
     &        BvgL*(mq1*PsL - PvR*sman)*((m3 - m4)*(m3 + m4)*(2*m1**2 + m3**2 + m4**2 - sman) + (-m3**2 + m4**2 + sman)*tman))*xv + 
     &     DRc*(-2*BsgL*m1**2*m3**2*PsL*tman - 2*Bsm1L*m1**2*m3**3*PsL*tman - BsgL*m3**4*PsL*tman - Bsm1L*m3**5*PsL*tman + 2*BsgL*m1**2*m4**2*PsL*tman + 
     &        2*Bsm2L*m1**2*m3*m4**2*PsL*tman - Bsm1L*m3**3*m4**2*PsL*tman + Bsm2L*m3**3*m4**2*PsL*tman + BsgL*m4**4*PsL*tman + Bsm2L*m3*m4**4*PsL*tman - 
     &        BvgL*m1**2*m3**2*PvR*sman*tman - Bsm1R*m2**2*m3**2*PvR*sman*tman - BvgL*m2**2*m3**2*PvR*sman*tman - Bvm1L*m1**2*m3**3*PvR*sman*tman - 
     &        Bvm1L*m2**2*m3**3*PvR*sman*tman - BvgL*m1**2*m4**2*PvR*sman*tman + Bsm2R*m2**2*m4**2*PvR*sman*tman + BvgL*m2**2*m4**2*PvR*sman*tman + 
     &        Bvm2L*m1**2*m3*m4**2*PvR*sman*tman + Bvm2L*m2**2*m3*m4**2*PvR*sman*tman - Bsm1R*m3**2*m4**2*PvR*sman*tman - 2*BvgL*m3**2*m4**2*PvR*sman*tman - 
     &        4*Bvm1L*m3**3*m4**2*PvR*sman*tman + Bsm2R*m4**4*PvR*sman*tman - 2*Bvm1L*m3*m4**4*PvR*sman*tman + 2*Bvm2L*m3*m4**4*PvR*sman*tman + 
     &        Bvm1R*m2**2*m3**2*mq1*PvR*sman*tman - Bvm2R*m2**2*m4**2*mq1*PvR*sman*tman + Bvm1R*m3**2*m4**2*mq1*PvR*sman*tman - 
     &        Bvm2R*m4**4*mq1*PvR*sman*tman + Bvm1R*m2**2*PsL*sman**2*tman + Bvm1R*m4**2*PsL*sman**2*tman + BvgL*m1**2*PvR*sman**2*tman + 
     &        Bsm1R*m2**2*PvR*sman**2*tman + BvgL*m2**2*PvR*sman**2*tman + Bvm1L*m1**2*m3*PvR*sman**2*tman + Bvm1L*m2**2*m3*PvR*sman**2*tman + 
     &        Bsm1R*m4**2*PvR*sman**2*tman + 2*BvgL*m4**2*PvR*sman**2*tman + 4*Bvm1L*m3*m4**2*PvR*sman**2*tman - Bvm1R*m2**2*mq1*PvR*sman**2*tman - 
     &        Bvm1R*m4**2*mq1*PvR*sman**2*tman + (Bsm2L*m1**2*m3*m4**2*PsL + 3*Bsm2L*m2**2*m3*m4**2*PsL + Bsm2L*m3**3*m4**2*PsL + 5*Bsm2L*m3*m4**4*PsL - 
     &           BsgL*(m3**2*(m1**2 + 3*m2**2 + m3**2) + (m1**2 + m2**2 + 4*m3**2)*m4**2 + m4**4)*PsL + BvgL*m1**2*m3**2*mq1*PsL + 
     &           3*BvgL*m2**2*m3**2*mq1*PsL + Bvm1L*m1**2*m3**3*mq1*PsL + 3*Bvm1L*m2**2*m3**3*mq1*PsL + BvgL*m3**4*mq1*PsL + Bvm1L*m3**5*mq1*PsL + 
     &           BvgL*m1**2*m4**2*mq1*PsL + BvgL*m2**2*m4**2*mq1*PsL - Bvm2L*m1**2*m3*m4**2*mq1*PsL - 3*Bvm2L*m2**2*m3*m4**2*mq1*PsL + 
     &           4*BvgL*m3**2*m4**2*mq1*PsL + 5*Bvm1L*m3**3*m4**2*mq1*PsL - Bvm2L*m3**3*m4**2*mq1*PsL + BvgL*m4**4*mq1*PsL - 5*Bvm2L*m3*m4**4*mq1*PsL + 
     &           2*BvgR*m2**2*m3*PsL*sman - Bvm1R*m1**2*m3**2*PsL*sman - Bvm1R*m2**2*m3**2*PsL*sman + 2*BvgR*m3**3*PsL*sman + Bvm1R*m3**4*PsL*sman + 
     &           Bvm2R*m1**2*m4**2*PsL*sman + Bvm2R*m2**2*m4**2*PsL*sman - 2*Bsm2L*m3*m4**2*PsL*sman + 4*BvgR*m3*m4**2*PsL*sman - Bvm1R*m3**2*m4**2*PsL*sman - 
     &           Bvm2R*m3**2*m4**2*PsL*sman + Bvm2R*m4**4*PsL*sman + BsgL*(m1**2 + m2**2 + m3**2 + m4**2)*PsL*sman - BvgL*m1**2*mq1*PsL*sman - 
     &           BvgL*m2**2*mq1*PsL*sman - Bvm1L*m1**2*m3*mq1*PsL*sman - 3*Bvm1L*m2**2*m3*mq1*PsL*sman - BvgL*m3**2*mq1*PsL*sman - 
     &           3*Bvm1L*m3**3*mq1*PsL*sman - BvgL*m4**2*mq1*PsL*sman - 5*Bvm1L*m3*m4**2*mq1*PsL*sman + 2*Bvm2L*m3*m4**2*mq1*PsL*sman + 
     &           2*BsgR*m2**2*m3*PvR*sman - Bsm1R*m1**2*m3**2*PvR*sman - BvgL*m1**2*m3**2*PvR*sman - Bsm1R*m2**2*m3**2*PvR*sman - 
     &           3*BvgL*m2**2*m3**2*PvR*sman + 2*BsgR*m3**3*PvR*sman - Bvm1L*m1**2*m3**3*PvR*sman - 3*Bvm1L*m2**2*m3**3*PvR*sman + Bsm1R*m3**4*PvR*sman - 
     &           BvgL*m3**4*PvR*sman - Bvm1L*m3**5*PvR*sman + Bsm2R*m1**2*m4**2*PvR*sman - BvgL*m1**2*m4**2*PvR*sman + Bsm2R*m2**2*m4**2*PvR*sman - 
     &           BvgL*m2**2*m4**2*PvR*sman + 4*BsgR*m3*m4**2*PvR*sman + Bvm2L*m1**2*m3*m4**2*PvR*sman + 3*Bvm2L*m2**2*m3*m4**2*PvR*sman - 
     &           Bsm1R*m3**2*m4**2*PvR*sman - Bsm2R*m3**2*m4**2*PvR*sman - 4*BvgL*m3**2*m4**2*PvR*sman - 5*Bvm1L*m3**3*m4**2*PvR*sman + 
     &           Bvm2L*m3**3*m4**2*PvR*sman + Bsm2R*m4**4*PvR*sman - BvgL*m4**4*PvR*sman + 5*Bvm2L*m3*m4**4*PvR*sman - 2*BvgR*m2**2*m3*mq1*PvR*sman + 
     &           Bvm1R*m1**2*m3**2*mq1*PvR*sman + Bvm1R*m2**2*m3**2*mq1*PvR*sman - 2*BvgR*m3**3*mq1*PvR*sman - Bvm1R*m3**4*mq1*PvR*sman - 
     &           Bvm2R*m1**2*m4**2*mq1*PvR*sman - Bvm2R*m2**2*m4**2*mq1*PvR*sman - 4*BvgR*m3*m4**2*mq1*PvR*sman + Bvm1R*m3**2*m4**2*mq1*PvR*sman + 
     &           Bvm2R*m3**2*m4**2*mq1*PvR*sman - Bvm2R*m4**4*mq1*PvR*sman + Bvm1R*m1**2*PsL*sman**2 + Bvm1R*m2**2*PsL*sman**2 - 2*BvgR*m3*PsL*sman**2 - 
     &           Bvm1R*m3**2*PsL*sman**2 + Bvm1R*m4**2*PsL*sman**2 + 2*Bvm1L*m3*mq1*PsL*sman**2 + Bsm1R*m1**2*PvR*sman**2 + BvgL*m1**2*PvR*sman**2 + 
     &           Bsm1R*m2**2*PvR*sman**2 + BvgL*m2**2*PvR*sman**2 - 2*BsgR*m3*PvR*sman**2 + Bvm1L*m1**2*m3*PvR*sman**2 + 3*Bvm1L*m2**2*m3*PvR*sman**2 - 
     &           Bsm1R*m3**2*PvR*sman**2 + BvgL*m3**2*PvR*sman**2 + 3*Bvm1L*m3**3*PvR*sman**2 + Bsm1R*m4**2*PvR*sman**2 + BvgL*m4**2*PvR*sman**2 + 
     &           5*Bvm1L*m3*m4**2*PvR*sman**2 - 2*Bvm2L*m3*m4**2*PvR*sman**2 - Bvm1R*m1**2*mq1*PvR*sman**2 - Bvm1R*m2**2*mq1*PvR*sman**2 + 
     &           2*BvgR*m3*mq1*PvR*sman**2 + Bvm1R*m3**2*mq1*PvR*sman**2 - Bvm1R*m4**2*mq1*PvR*sman**2 - 2*Bvm1L*m3*PvR*sman**3 - 
     &           Bsm1L*m3*PsL*(m3**2 - sman)*(m1**2 + 3*m2**2 + m3**2 + 5*m4**2 - 2*sman - tman) + BsgL*PsL*(m3**2 + m4**2 - sman)*tman - 
     &           (Bsm2L*m3*m4**2*PsL + (BvgL*m3**2 + Bvm1L*m3**3 + BvgL*m4**2 - Bvm2L*m3*m4**2 - (BvgL + Bvm1L*m3)*sman)*(mq1*PsL - PvR*sman))*tman)*uman + 
     &        ((Bsm1L*m3**3 - Bsm2L*m3*m4**2 + BsgL*(m3**2 + m4**2) - (m3**2*(BvgL + Bvm1L*m3) + (BvgL - Bvm2L*m3)*m4**2)*mq1)*PsL + 
     &           ((-BsgL - Bvm2R*m4**2 + BvgL*mq1 + m3*(-Bsm1L + Bvm1R*m3 + Bvm1L*mq1))*PsL + 
     &              (Bsm1R*m3**2 + Bvm1L*m3**3 - Bsm2R*m4**2 - Bvm2L*m3*m4**2 + BvgL*(m3**2 + m4**2) - Bvm1R*m3**2*mq1 + Bvm2R*m4**2*mq1)*PvR)*sman - 
     &           (Bvm1R*PsL + (Bsm1R + BvgL + Bvm1L*m3 - Bvm1R*mq1)*PvR)*sman**2)*uman**2)*xv - 
     &     (-2*m3*m4**2*((Bvm1L + Bvm2L)*m4**2*mq1*PsL - (BvgR*PsL + (BsgR + (Bvm1L + Bvm2L)*m4**2 - BvgR*mq1)*PvR)*sman)*
     &         (-2*DLc*m1*mq3c + DRc*(2*m1**2 + m2**2 + m3**2 - sman - tman)) + 
     &        DRc*m3*(2*(BvgR*PsL + BsgR*PvR - BvgR*mq1*PvR)*sman + Bvm1L*(m3**2 - sman)*(mq1*PsL - PvR*sman) + Bvm2L*m4**2*(-(mq1*PsL) + PvR*sman))*tman*
     &         uman + Bsm2L*m3*m4**2*PsL*(2*m4**2*(-2*DLc*m1*mq3c + DRc*(2*m1**2 + m2**2 + m3**2 - sman - tman)) + DRc*tman*uman) + 
     &        BvgL*(mq1*PsL - PvR*sman)*(2*DLc*m1*m4**2*mq3c*(-m3**2 + m4**2 - sman) + 
     &           DRc*(m4**2*(m3**2 + m4**2 - sman - tman)*(2*m3**2 + sman - tman) + (m4**4 - m4**2*tman + sman*tman + m3**2*(-m4**2 + tman))*uman - 
     &              m4**2*uman**2 + m1**2*m4**2*(4*m3**2 + sman - tman + uman) + m2**2*m4**2*(2*m3**2 - 2*m4**2 + sman - tman + uman))))*xv - 
     &     DRc*(-4*BvgR*m3*PsL*sman*tman + Bsm1R*m1**2*PvR*sman*tman + Bsm2R*m1**2*PvR*sman*tman + 4*BvgL*m1**2*PvR*sman*tman + 5*Bsm1R*m2**2*PvR*sman*tman + 
     &        Bsm2R*m2**2*PvR*sman*tman + 8*BvgL*m2**2*PvR*sman*tman - 4*BsgR*m3*PvR*sman*tman + 5*Bvm1L*m1**2*m3*PvR*sman*tman + 
     &        3*Bvm2L*m1**2*m3*PvR*sman*tman + 7*Bvm1L*m2**2*m3*PvR*sman*tman + Bvm2L*m2**2*m3*PvR*sman*tman + Bsm1R*m3**2*PvR*sman*tman + 
     &        Bsm2R*m3**2*PvR*sman*tman + 8*BvgL*m3**2*PvR*sman*tman + 3*Bvm1L*m3**3*PvR*sman*tman + Bvm2L*m3**3*PvR*sman*tman + Bsm1R*m4**2*PvR*sman*tman + 
     &        Bsm2R*m4**2*PvR*sman*tman + 4*BvgL*m4**2*PvR*sman*tman + 5*Bvm1L*m3*m4**2*PvR*sman*tman + 3*Bvm2L*m3*m4**2*PvR*sman*tman + 
     &        4*BvgR*m3*mq1*PvR*sman*tman - 4*BsgL*PsL*tman**2 - 3*Bsm1L*m3*PsL*tman**2 - Bsm2L*m3*PsL*tman**2 + 4*BvgL*mq1*PsL*tman**2 + 
     &        3*Bvm1L*m3*mq1*PsL*tman**2 + Bvm2L*m3*mq1*PsL*tman**2 - Bsm1R*PvR*sman*tman**2 - Bsm2R*PvR*sman*tman**2 - 4*BvgL*PvR*sman*tman**2 - 
     &        3*Bvm1L*m3*PvR*sman*tman**2 - Bvm2L*m3*PvR*sman*tman**2 + 4*BsgL*m1**2*PsL*uman + 7*Bsm1L*m1**2*m3*PsL*uman + Bsm2L*m1**2*m3*PsL*uman + 
     &        5*Bsm1L*m2**2*m3*PsL*uman + 3*Bsm2L*m2**2*m3*PsL*uman + 4*BsgL*m3**2*PsL*uman + 3*Bsm1L*m3**3*PsL*uman + Bsm2L*m3**3*PsL*uman + 
     &        7*Bsm1L*m3*m4**2*PsL*uman + 5*Bsm2L*m3*m4**2*PsL*uman - 4*BvgL*m1**2*mq1*PsL*uman - 7*Bvm1L*m1**2*m3*mq1*PsL*uman - 
     &        Bvm2L*m1**2*m3*mq1*PsL*uman - 5*Bvm1L*m2**2*m3*mq1*PsL*uman - 3*Bvm2L*m2**2*m3*mq1*PsL*uman - 4*BvgL*m3**2*mq1*PsL*uman - 
     &        3*Bvm1L*m3**3*mq1*PsL*uman - Bvm2L*m3**3*mq1*PsL*uman - 7*Bvm1L*m3*m4**2*mq1*PsL*uman - 5*Bvm2L*m3*m4**2*mq1*PsL*uman - 
     &        2*Bsm1L*m3*PsL*sman*uman - 2*Bsm2L*m3*PsL*sman*uman + 2*Bvm1L*m3*mq1*PsL*sman*uman + 2*Bvm2L*m3*mq1*PsL*sman*uman + 
     &        5*Bsm1R*m1**2*PvR*sman*uman + Bsm2R*m1**2*PvR*sman*uman + 4*BvgL*m1**2*PvR*sman*uman + Bsm1R*m2**2*PvR*sman*uman + Bsm2R*m2**2*PvR*sman*uman + 
     &        7*Bvm1L*m1**2*m3*PvR*sman*uman + Bvm2L*m1**2*m3*PvR*sman*uman + 5*Bvm1L*m2**2*m3*PvR*sman*uman + 3*Bvm2L*m2**2*m3*PvR*sman*uman + 
     &        3*Bsm1R*m3**2*PvR*sman*uman - Bsm2R*m3**2*PvR*sman*uman + 4*BvgL*m3**2*PvR*sman*uman + 3*Bvm1L*m3**3*PvR*sman*uman + Bvm2L*m3**3*PvR*sman*uman + 
     &        Bsm1R*m4**2*PvR*sman*uman + Bsm2R*m4**2*PvR*sman*uman + 7*Bvm1L*m3*m4**2*PvR*sman*uman + 5*Bvm2L*m3*m4**2*PvR*sman*uman - 
     &        2*Bvm1L*m3*PvR*sman**2*uman - 2*Bvm2L*m3*PvR*sman**2*uman - 4*BsgL*PsL*tman*uman - 6*Bsm1L*m3*PsL*tman*uman - 2*Bsm2L*m3*PsL*tman*uman + 
     &        4*BvgL*mq1*PsL*tman*uman + 6*Bvm1L*m3*mq1*PsL*tman*uman + 2*Bvm2L*m3*mq1*PsL*tman*uman - 4*Bsm1R*PvR*sman*tman*uman - 
     &        4*BvgL*PvR*sman*tman*uman - 6*Bvm1L*m3*PvR*sman*tman*uman - 2*Bvm2L*m3*PvR*sman*tman*uman - 3*Bsm1L*m3*PsL*uman**2 - Bsm2L*m3*PsL*uman**2 + 
     &        3*Bvm1L*m3*mq1*PsL*uman**2 + Bvm2L*m3*mq1*PsL*uman**2 - Bsm1R*PvR*sman*uman**2 - Bsm2R*PvR*sman*uman**2 - 3*Bvm1L*m3*PvR*sman*uman**2 - 
     &        Bvm2L*m3*PvR*sman*uman**2 + Bvm2R*(PsL - mq1*PvR)*sman*
     &         ((m1**2 + m2**2 + m3**2 + m4**2 - tman)*tman + (m1**2 + m2**2 - m3**2 + m4**2)*uman - uman**2) + 
     &        Bvm1R*(PsL - mq1*PvR)*sman*((m1**2 + 5*m2**2 + m3**2 + m4**2 - tman)*tman + (5*m1**2 + m2**2 + 3*m3**2 + m4**2 - 4*tman)*uman - uman**2) - 
     &        4*BsgL*m1**2*m2**2*m3**2*PsL*xv) + BsgL*PsL*(-4*DRc*
     &         ((m3**2 + m4**2 - sman - tman)*(2*m3**2 + sman - tman) + (-m3**2 + m4**2)*uman - uman**2 + m1**2*(4*m3**2 + sman - tman + uman) + 
     &           m2**2*(2*m3**2 - 2*m4**2 + sman - tman + uman)) + 
     &        DRc*(m4**2*(m3**2 + m4**2 - sman - tman)*(2*m3**2 + sman - tman) + (m4**4 - m4**2*tman + sman*tman + m3**2*(-m4**2 + tman))*uman - 
     &           m4**2*uman**2 + m1**2*m4**2*(4*m3**2 + sman - tman + uman) + m2**2*m4**2*(2*m3**2 - 2*m4**2 + sman - tman + uman))*xv - 
     &        2*DLc*m1*mq3c*(m3**2 - m4**2 + sman)*(-4 + m4**2*xv)) + 
     &     Bsm1L*m3*PsL*(2*DRc*m4**2*(2*m1**2 + m2**2 + m3**2 - sman - tman) + 
     &        DRc*(-2*m4**4*(2*m1**2 + m2**2 + m3**2 - sman - tman) + (m3**2 - sman)*tman*uman)*xv + 4*DLc*m1*m4**2*mq3c*(-1 + m4**2*xv))) + 
     &  A1vL*A3Lc*(4*Bsm2R*DLc*m1**2*m2**2*m3*PsR + 8*Bsm2R*DLc*m1**2*m3*m4**2*PsR + 4*Bsm2R*DLc*m2**2*m3*m4**2*PsR + 4*Bsm2R*DLc*m3**3*m4**2*PsR + 
     &     4*Bsm2R*DLc*m3*m4**4*PsR - 8*BvgR*DLc*m1**2*m2**2*mq1*PsR - 12*Bvm1R*DLc*m1**2*m2**2*m3*mq1*PsR - 4*Bvm2R*DLc*m1**2*m2**2*m3*mq1*PsR + 
     &     8*BvgR*DLc*m1**2*m3**2*mq1*PsR - 4*BvgR*DLc*m2**2*m3**2*mq1*PsR - 8*Bvm1R*DLc*m2**2*m3**3*mq1*PsR + 4*BvgR*DLc*m3**4*mq1*PsR - 
     &     4*BvgR*DLc*m2**2*m4**2*mq1*PsR - 8*Bvm1R*DLc*m1**2*m3*m4**2*mq1*PsR - 8*Bvm2R*DLc*m1**2*m3*m4**2*mq1*PsR - 4*Bvm1R*DLc*m2**2*m3*m4**2*mq1*PsR - 
     &     4*Bvm2R*DLc*m2**2*m3*m4**2*mq1*PsR + 4*BvgR*DLc*m3**2*m4**2*mq1*PsR - 4*Bvm1R*DLc*m3**3*m4**2*mq1*PsR - 4*Bvm2R*DLc*m3**3*m4**2*mq1*PsR - 
     &     4*Bvm1R*DLc*m3*m4**4*mq1*PsR - 4*Bvm2R*DLc*m3*m4**4*mq1*PsR + 2*Bsm2R*DRc*m1*m3**3*mq3c*PsR - 2*Bsm2R*DRc*m1*m3*m4**2*mq3c*PsR - 
     &     4*BvgR*DRc*m1*m3**2*mq1*mq3c*PsR + 2*Bvm1R*DRc*m1*m3**3*mq1*mq3c*PsR - 2*Bvm2R*DRc*m1*m3**3*mq1*mq3c*PsR + 4*BvgR*DRc*m1*m4**2*mq1*mq3c*PsR - 
     &     2*Bvm1R*DRc*m1*m3*m4**2*mq1*mq3c*PsR + 2*Bvm2R*DRc*m1*m3*m4**2*mq1*mq3c*PsR + 6*Bvm1L*DLc*m1**2*m2**2*PsR*sman + 2*Bvm2L*DLc*m1**2*m2**2*PsR*sman - 
     &     Bsm2R*DLc*m1**2*m3*PsR*sman + 8*BvgL*DLc*m1**2*m3*PsR*sman + Bsm2R*DLc*m2**2*m3*PsR*sman + 4*BvgL*DLc*m2**2*m3*PsR*sman + 
     &     2*Bvm1L*DLc*m1**2*m3**2*PsR*sman - 2*Bvm2L*DLc*m1**2*m3**2*PsR*sman + 8*Bvm1L*DLc*m2**2*m3**2*PsR*sman + Bsm2R*DLc*m3**3*PsR*sman + 
     &     4*BvgL*DLc*m3**3*PsR*sman - 3*Bsm2R*DLc*m3*m4**2*PsR*sman - Bvm1R*DLc*m1**2*m3*mq1*PsR*sman + Bvm2R*DLc*m1**2*m3*mq1*PsR*sman - 
     &     3*Bvm1R*DLc*m2**2*m3*mq1*PsR*sman - Bvm2R*DLc*m2**2*m3*mq1*PsR*sman + Bvm1R*DLc*m3**3*mq1*PsR*sman - Bvm2R*DLc*m3**3*mq1*PsR*sman + 
     &     Bvm1R*DLc*m3*m4**2*mq1*PsR*sman + 3*Bvm2R*DLc*m3*m4**2*mq1*PsR*sman - 2*Bsm2R*DRc*m1*m3*mq3c*PsR*sman - 16*BvgL*DRc*m1*m3*mq3c*PsR*sman - 
     &     6*Bvm1L*DRc*m1*m3**2*mq3c*PsR*sman + 2*Bvm2L*DRc*m1*m3**2*mq3c*PsR*sman + 2*Bvm1L*DRc*m1*m4**2*mq3c*PsR*sman + 2*Bvm2L*DRc*m1*m4**2*mq3c*PsR*sman - 
     &     4*BvgR*DRc*m1*mq1*mq3c*PsR*sman + 6*Bvm1R*DRc*m1*m3*mq1*mq3c*PsR*sman + 2*Bvm2R*DRc*m1*m3*mq1*mq3c*PsR*sman + 6*Bsm1L*DLc*m1**2*m2**2*PvL*sman + 
     &     2*Bsm2L*DLc*m1**2*m2**2*PvL*sman + 8*BvgR*DLc*m1**2*m2**2*PvL*sman + 8*BsgL*DLc*m1**2*m3*PvL*sman + 4*BsgL*DLc*m2**2*m3*PvL*sman + 
     &     12*Bvm1R*DLc*m1**2*m2**2*m3*PvL*sman + 4*Bvm2R*DLc*m1**2*m2**2*m3*PvL*sman + 2*Bsm1L*DLc*m1**2*m3**2*PvL*sman - 2*Bsm2L*DLc*m1**2*m3**2*PvL*sman - 
     &     8*BvgR*DLc*m1**2*m3**2*PvL*sman + 8*Bsm1L*DLc*m2**2*m3**2*PvL*sman + 4*BvgR*DLc*m2**2*m3**2*PvL*sman + 4*BsgL*DLc*m3**3*PvL*sman + 
     &     8*Bvm1R*DLc*m2**2*m3**3*PvL*sman - 4*BvgR*DLc*m3**4*PvL*sman + 4*BvgR*DLc*m2**2*m4**2*PvL*sman + 8*Bvm1R*DLc*m1**2*m3*m4**2*PvL*sman + 
     &     8*Bvm2R*DLc*m1**2*m3*m4**2*PvL*sman + 4*Bvm1R*DLc*m2**2*m3*m4**2*PvL*sman + 4*Bvm2R*DLc*m2**2*m3*m4**2*PvL*sman - 4*BvgR*DLc*m3**2*m4**2*PvL*sman + 
     &     4*Bvm1R*DLc*m3**3*m4**2*PvL*sman + 4*Bvm2R*DLc*m3**3*m4**2*PvL*sman + 4*Bvm1R*DLc*m3*m4**4*PvL*sman + 4*Bvm2R*DLc*m3*m4**4*PvL*sman - 
     &     6*Bvm1L*DLc*m1**2*m2**2*mq1*PvL*sman - 2*Bvm2L*DLc*m1**2*m2**2*mq1*PvL*sman - 8*BvgL*DLc*m1**2*m3*mq1*PvL*sman - 4*BvgL*DLc*m2**2*m3*mq1*PvL*sman - 
     &     2*Bvm1L*DLc*m1**2*m3**2*mq1*PvL*sman + 2*Bvm2L*DLc*m1**2*m3**2*mq1*PvL*sman - 8*Bvm1L*DLc*m2**2*m3**2*mq1*PvL*sman - 
     &     4*BvgL*DLc*m3**3*mq1*PvL*sman - 16*BsgL*DRc*m1*m3*mq3c*PvL*sman - 6*Bsm1L*DRc*m1*m3**2*mq3c*PvL*sman + 2*Bsm2L*DRc*m1*m3**2*mq3c*PvL*sman + 
     &     4*BvgR*DRc*m1*m3**2*mq3c*PvL*sman - 2*Bvm1R*DRc*m1*m3**3*mq3c*PvL*sman + 2*Bvm2R*DRc*m1*m3**3*mq3c*PvL*sman + 2*Bsm1L*DRc*m1*m4**2*mq3c*PvL*sman + 
     &     2*Bsm2L*DRc*m1*m4**2*mq3c*PvL*sman - 4*BvgR*DRc*m1*m4**2*mq3c*PvL*sman + 2*Bvm1R*DRc*m1*m3*m4**2*mq3c*PvL*sman - 
     &     2*Bvm2R*DRc*m1*m3*m4**2*mq3c*PvL*sman + 16*BvgL*DRc*m1*m3*mq1*mq3c*PvL*sman + 6*Bvm1L*DRc*m1*m3**2*mq1*mq3c*PvL*sman - 
     &     2*Bvm2L*DRc*m1*m3**2*mq1*mq3c*PvL*sman - 2*Bvm1L*DRc*m1*m4**2*mq1*mq3c*PvL*sman - 2*Bvm2L*DRc*m1*m4**2*mq1*mq3c*PvL*sman + 
     &     Bvm1L*DLc*m1**2*PsR*sman**2 + Bvm2L*DLc*m1**2*PsR*sman**2 + Bvm1L*DLc*m2**2*PsR*sman**2 + Bvm2L*DLc*m2**2*PsR*sman**2 - Bsm2R*DLc*m3*PsR*sman**2 - 
     &     4*BvgL*DLc*m3*PsR*sman**2 - 3*Bvm1L*DLc*m3**2*PsR*sman**2 + Bvm2L*DLc*m3**2*PsR*sman**2 + Bvm1L*DLc*m4**2*PsR*sman**2 + 
     &     Bvm2L*DLc*m4**2*PsR*sman**2 + 3*Bvm1R*DLc*m3*mq1*PsR*sman**2 + Bvm2R*DLc*m3*mq1*PsR*sman**2 - 2*Bvm1L*DRc*m1*mq3c*PsR*sman**2 - 
     &     2*Bvm2L*DRc*m1*mq3c*PsR*sman**2 + Bsm1L*DLc*m1**2*PvL*sman**2 + Bsm2L*DLc*m1**2*PvL*sman**2 + Bsm1L*DLc*m2**2*PvL*sman**2 + 
     &     Bsm2L*DLc*m2**2*PvL*sman**2 - 4*BsgL*DLc*m3*PvL*sman**2 + Bvm1R*DLc*m1**2*m3*PvL*sman**2 - Bvm2R*DLc*m1**2*m3*PvL*sman**2 + 
     &     3*Bvm1R*DLc*m2**2*m3*PvL*sman**2 + Bvm2R*DLc*m2**2*m3*PvL*sman**2 - 3*Bsm1L*DLc*m3**2*PvL*sman**2 + Bsm2L*DLc*m3**2*PvL*sman**2 - 
     &     Bvm1R*DLc*m3**3*PvL*sman**2 + Bvm2R*DLc*m3**3*PvL*sman**2 + Bsm1L*DLc*m4**2*PvL*sman**2 + Bsm2L*DLc*m4**2*PvL*sman**2 - 
     &     Bvm1R*DLc*m3*m4**2*PvL*sman**2 - 3*Bvm2R*DLc*m3*m4**2*PvL*sman**2 - Bvm1L*DLc*m1**2*mq1*PvL*sman**2 - Bvm2L*DLc*m1**2*mq1*PvL*sman**2 - 
     &     Bvm1L*DLc*m2**2*mq1*PvL*sman**2 - Bvm2L*DLc*m2**2*mq1*PvL*sman**2 + 4*BvgL*DLc*m3*mq1*PvL*sman**2 + 3*Bvm1L*DLc*m3**2*mq1*PvL*sman**2 - 
     &     Bvm2L*DLc*m3**2*mq1*PvL*sman**2 - Bvm1L*DLc*m4**2*mq1*PvL*sman**2 - Bvm2L*DLc*m4**2*mq1*PvL*sman**2 - 2*Bsm1L*DRc*m1*mq3c*PvL*sman**2 - 
     &     2*Bsm2L*DRc*m1*mq3c*PvL*sman**2 + 4*BvgR*DRc*m1*mq3c*PvL*sman**2 - 6*Bvm1R*DRc*m1*m3*mq3c*PvL*sman**2 - 2*Bvm2R*DRc*m1*m3*mq3c*PvL*sman**2 + 
     &     2*Bvm1L*DRc*m1*mq1*mq3c*PvL*sman**2 + 2*Bvm2L*DRc*m1*mq1*mq3c*PvL*sman**2 - Bvm1L*DLc*PsR*sman**3 - Bvm2L*DLc*PsR*sman**3 - Bsm1L*DLc*PvL*sman**3 - 
     &     Bsm2L*DLc*PvL*sman**3 - 3*Bvm1R*DLc*m3*PvL*sman**3 - Bvm2R*DLc*m3*PvL*sman**3 + Bvm1L*DLc*mq1*PvL*sman**3 + Bvm2L*DLc*mq1*PvL*sman**3 + 
     &     Bsm1R*m3*PsR*(-2*DRc*m1*mq3c*(m3**2 - 3*m4**2 + 3*sman) + 
     &        DLc*(2*m3**2*m4**2 + 4*m4**4 - m3**2*sman + m4**2*sman - 3*sman**2 + m2**2*(8*m3**2 + 2*m4**2 + 3*sman - 7*tman) + 
     &           m1**2*(12*m2**2 + 4*m4**2 + sman - 5*tman))) - Bsm2R*DLc*(3*m1**2 + m2**2)*m3*PsR*tman - 2*Bsm2R*DLc*m3*m4**2*PsR*tman - 
     &     4*BvgR*DLc*m1**2*mq1*PsR*tman - 4*BvgR*DLc*m2**2*mq1*PsR*tman - 12*BvgR*DLc*m3**2*mq1*PsR*tman - 4*BvgR*DLc*m4**2*mq1*PsR*tman + 
     &     2*Bvm1R*DLc*m3*m4**2*mq1*PsR*tman + 2*Bvm2R*DLc*m3*m4**2*mq1*PsR*tman - 8*BvgL*DLc*m3*PsR*sman*tman + 4*BvgR*DLc*m1**2*PvL*sman*tman + 
     &     4*BvgR*DLc*m2**2*PvL*sman*tman - 8*BsgL*DLc*m3*PvL*sman*tman + 12*BvgR*DLc*m3**2*PvL*sman*tman + 4*BvgR*DLc*m4**2*PvL*sman*tman - 
     &     2*Bvm1R*DLc*m3*m4**2*PvL*sman*tman - 2*Bvm2R*DLc*m3*m4**2*PvL*sman*tman + 8*BvgL*DLc*m3*mq1*PvL*sman*tman + 4*BvgR*DLc*mq1*PsR*tman**2 - 
     &     4*BvgR*DLc*PvL*sman*tman**2 + 4*BsgR*PsR*(-(DRc*m1*mq3c*(m3**2 - m4**2 + sman)) + 
     &        DLc*(m3**4 + m3**2*m4**2 - m3**2*sman + m4**2*sman - sman**2 + m2**2*(3*m3**2 - m4**2 + sman - 2*tman) + 
     &           m1**2*(2*m2**2 + 2*m3**2 + sman - tman) - 2*m3**2*tman)) + 4*BvgR*DLc*m1**2*mq1*PsR*uman + 4*BvgR*DLc*m2**2*mq1*PsR*uman - 
     &     4*BvgR*DLc*m3**2*mq1*PsR*uman + 4*BvgR*DLc*m4**2*mq1*PsR*uman - 4*BvgR*DLc*m1**2*PvL*sman*uman - 4*BvgR*DLc*m2**2*PvL*sman*uman + 
     &     4*BvgR*DLc*m3**2*PvL*sman*uman - 4*BvgR*DLc*m4**2*PvL*sman*uman - 4*BvgR*DLc*mq1*PsR*uman**2 + 4*BvgR*DLc*PvL*sman*uman**2 - 
     &     DLc*(4*BsgR*m4**2*PsR*tman - 4*BvgR*m1**2*mq1*PsR*tman - 8*BvgR*m2**2*mq1*PsR*tman - 5*Bvm1R*m1**2*m3*mq1*PsR*tman - 
     &        3*Bvm2R*m1**2*m3*mq1*PsR*tman - 7*Bvm1R*m2**2*m3*mq1*PsR*tman - Bvm2R*m2**2*m3*mq1*PsR*tman - 8*BvgR*m3**2*mq1*PsR*tman - 
     &        3*Bvm1R*m3**3*mq1*PsR*tman - Bvm2R*m3**3*mq1*PsR*tman - 4*BvgR*m4**2*mq1*PsR*tman - 5*Bvm1R*m3*m4**2*mq1*PsR*tman - 
     &        3*Bvm2R*m3*m4**2*mq1*PsR*tman + Bvm1L*m1**2*PsR*sman*tman + Bvm2L*m1**2*PsR*sman*tman + 5*Bvm1L*m2**2*PsR*sman*tman + 
     &        Bvm2L*m2**2*PsR*sman*tman - 4*BvgL*m3*PsR*sman*tman + Bvm1L*m3**2*PsR*sman*tman + Bvm2L*m3**2*PsR*sman*tman + Bvm1L*m4**2*PsR*sman*tman + 
     &        Bvm2L*m4**2*PsR*sman*tman + Bsm1L*m1**2*PvL*sman*tman + Bsm2L*m1**2*PvL*sman*tman + 4*BvgR*m1**2*PvL*sman*tman + 5*Bsm1L*m2**2*PvL*sman*tman + 
     &        Bsm2L*m2**2*PvL*sman*tman + 8*BvgR*m2**2*PvL*sman*tman - 4*BsgL*m3*PvL*sman*tman + 5*Bvm1R*m1**2*m3*PvL*sman*tman + 
     &        3*Bvm2R*m1**2*m3*PvL*sman*tman + 7*Bvm1R*m2**2*m3*PvL*sman*tman + Bvm2R*m2**2*m3*PvL*sman*tman + Bsm1L*m3**2*PvL*sman*tman + 
     &        Bsm2L*m3**2*PvL*sman*tman + 8*BvgR*m3**2*PvL*sman*tman + 3*Bvm1R*m3**3*PvL*sman*tman + Bvm2R*m3**3*PvL*sman*tman + Bsm1L*m4**2*PvL*sman*tman + 
     &        Bsm2L*m4**2*PvL*sman*tman + 4*BvgR*m4**2*PvL*sman*tman + 5*Bvm1R*m3*m4**2*PvL*sman*tman + 3*Bvm2R*m3*m4**2*PvL*sman*tman - 
     &        Bvm1L*m1**2*mq1*PvL*sman*tman - Bvm2L*m1**2*mq1*PvL*sman*tman - 5*Bvm1L*m2**2*mq1*PvL*sman*tman - Bvm2L*m2**2*mq1*PvL*sman*tman + 
     &        4*BvgL*m3*mq1*PvL*sman*tman - Bvm1L*m3**2*mq1*PvL*sman*tman - Bvm2L*m3**2*mq1*PvL*sman*tman - Bvm1L*m4**2*mq1*PvL*sman*tman - 
     &        Bvm2L*m4**2*mq1*PvL*sman*tman - 4*BsgR*PsR*tman**2 + 4*BvgR*mq1*PsR*tman**2 + 3*Bvm1R*m3*mq1*PsR*tman**2 + Bvm2R*m3*mq1*PsR*tman**2 - 
     &        Bvm1L*PsR*sman*tman**2 - Bvm2L*PsR*sman*tman**2 - Bsm1L*PvL*sman*tman**2 - Bsm2L*PvL*sman*tman**2 - 4*BvgR*PvL*sman*tman**2 - 
     &        3*Bvm1R*m3*PvL*sman*tman**2 - Bvm2R*m3*PvL*sman*tman**2 + Bvm1L*mq1*PvL*sman*tman**2 + Bvm2L*mq1*PvL*sman*tman**2 + 
     &        (-7*Bvm1R*m1**2*m3*mq1*PsR - Bvm2R*m1**2*m3*mq1*PsR - 5*Bvm1R*m2**2*m3*mq1*PsR - 3*Bvm2R*m2**2*m3*mq1*PsR - 3*Bvm1R*m3**3*mq1*PsR - 
     &           Bvm2R*m3**3*mq1*PsR - 7*Bvm1R*m3*m4**2*mq1*PsR - 5*Bvm2R*m3*m4**2*mq1*PsR + 5*Bvm1L*m1**2*PsR*sman + Bvm2L*m1**2*PsR*sman + 
     &           Bvm1L*m2**2*PsR*sman + Bvm2L*m2**2*PsR*sman + 3*Bvm1L*m3**2*PsR*sman - Bvm2L*m3**2*PsR*sman + Bvm1L*m4**2*PsR*sman + Bvm2L*m4**2*PsR*sman + 
     &           2*Bvm1R*m3*mq1*PsR*sman + 2*Bvm2R*m3*mq1*PsR*sman + 5*Bsm1L*m1**2*PvL*sman + Bsm2L*m1**2*PvL*sman + Bsm1L*m2**2*PvL*sman + 
     &           Bsm2L*m2**2*PvL*sman + 7*Bvm1R*m1**2*m3*PvL*sman + Bvm2R*m1**2*m3*PvL*sman + 5*Bvm1R*m2**2*m3*PvL*sman + 3*Bvm2R*m2**2*m3*PvL*sman + 
     &           3*Bsm1L*m3**2*PvL*sman - Bsm2L*m3**2*PvL*sman + 3*Bvm1R*m3**3*PvL*sman + Bvm2R*m3**3*PvL*sman + Bsm1L*m4**2*PvL*sman + Bsm2L*m4**2*PvL*sman + 
     &           7*Bvm1R*m3*m4**2*PvL*sman + 5*Bvm2R*m3*m4**2*PvL*sman - 5*Bvm1L*m1**2*mq1*PvL*sman - Bvm2L*m1**2*mq1*PvL*sman - Bvm1L*m2**2*mq1*PvL*sman - 
     &           Bvm2L*m2**2*mq1*PvL*sman - 3*Bvm1L*m3**2*mq1*PvL*sman + Bvm2L*m3**2*mq1*PvL*sman - Bvm1L*m4**2*mq1*PvL*sman - Bvm2L*m4**2*mq1*PvL*sman - 
     &           2*Bvm1R*m3*PvL*sman**2 - 2*Bvm2R*m3*PvL*sman**2 + 4*BsgR*PsR*(m1**2 + m3**2 - tman) - 4*BvgR*(mq1*PsR - PvL*sman)*(m1**2 + m3**2 - tman) + 
     &           2*((3*Bvm1R + Bvm2R)*m3*mq1*PsR - (2*Bvm1L*PsR + (2*Bsm1L + 3*Bvm1R*m3 + Bvm2R*m3 - 2*Bvm1L*mq1)*PvL)*sman)*tman)*uman + 
     &        ((3*Bvm1R + Bvm2R)*m3*mq1*PsR - ((Bsm1L + Bsm2L + (3*Bvm1R + Bvm2R)*m3)*PvL + (Bvm1L + Bvm2L)*(PsR - mq1*PvL))*sman)*uman**2 + 
     &        Bsm2R*m3*PsR*((m1**2 + 3*m2**2 - 2*sman)*uman + m3**2*(tman + uman) - (tman + uman)**2 + m4**2*(3*tman + 5*uman)) + 
     &        Bsm1R*m3*PsR*((7*m1**2 + 5*m2**2 - 2*sman)*uman + 3*m3**2*(tman + uman) - 3*(tman + uman)**2 + m4**2*(5*tman + 7*uman))) + 
     &     (-4*Bsm2R*DLc*m1**2*m2**2*m3*m4**2*PsR - 4*Bsm2R*DLc*m1**2*m3*m4**4*PsR - 2*Bsm2R*DLc*m2**2*m3*m4**4*PsR - 2*Bsm2R*DLc*m3**3*m4**4*PsR - 
     &        4*Bsm2R*DLc*m3*m4**6*PsR - 4*BvgR*DLc*m1**2*m2**2*m3**2*mq1*PsR - 4*Bvm1R*DLc*m1**2*m2**2*m3**3*mq1*PsR + 
     &        4*Bvm2R*DLc*m1**2*m2**2*m3*m4**2*mq1*PsR - 4*BvgR*DLc*m1**2*m3**2*m4**2*mq1*PsR - 2*BvgR*DLc*m2**2*m3**2*m4**2*mq1*PsR - 
     &        8*Bvm1R*DLc*m1**2*m3**3*m4**2*mq1*PsR - 4*Bvm1R*DLc*m2**2*m3**3*m4**2*mq1*PsR - 2*BvgR*DLc*m3**4*m4**2*mq1*PsR - 
     &        4*Bvm1R*DLc*m3**5*m4**2*mq1*PsR - 2*BvgR*DLc*m2**2*m4**4*mq1*PsR - 4*Bvm1R*DLc*m1**2*m3*m4**4*mq1*PsR + 4*Bvm2R*DLc*m1**2*m3*m4**4*mq1*PsR - 
     &        2*Bvm1R*DLc*m2**2*m3*m4**4*mq1*PsR + 2*Bvm2R*DLc*m2**2*m3*m4**4*mq1*PsR - 2*BvgR*DLc*m3**2*m4**4*mq1*PsR - 6*Bvm1R*DLc*m3**3*m4**4*mq1*PsR + 
     &        2*Bvm2R*DLc*m3**3*m4**4*mq1*PsR + 4*Bvm2R*DLc*m3*m4**6*mq1*PsR - 2*Bsm2R*DRc*m1*m3**3*m4**2*mq3c*PsR - 2*Bsm2R*DRc*m1*m3*m4**4*mq3c*PsR - 
     &        2*BvgR*DRc*m1*m3**4*mq1*mq3c*PsR - 2*Bvm1R*DRc*m1*m3**5*mq1*mq3c*PsR + 2*Bvm1R*DRc*m1*m3**3*m4**2*mq1*mq3c*PsR + 
     &        2*Bvm2R*DRc*m1*m3**3*m4**2*mq1*mq3c*PsR + 2*BvgR*DRc*m1*m4**4*mq1*mq3c*PsR + 4*Bvm1R*DRc*m1*m3*m4**4*mq1*mq3c*PsR + 
     &        2*Bvm2R*DRc*m1*m3*m4**4*mq1*mq3c*PsR + (-2*Bvm2L*DLc*m1**2*m2**2*m4**2*PsR + Bsm2R*DLc*m1**2*m3*m4**2*PsR - Bsm2R*DLc*m2**2*m3*m4**2*PsR + 
     &           2*Bvm2L*DLc*m1**2*m3**2*m4**2*PsR - Bsm2R*DLc*m3**3*m4**2*PsR + Bsm2R*DLc*m3*m4**4*PsR + 2*BvgR*DLc*m1**2*m2**2*mq1*PsR + 
     &           4*Bvm1R*DLc*m1**2*m2**2*m3*mq1*PsR - BvgR*DLc*m1**2*m3**2*mq1*PsR - BvgR*DLc*m2**2*m3**2*mq1*PsR + Bvm1R*DLc*m1**2*m3**3*mq1*PsR - 
     &           Bvm1R*DLc*m2**2*m3**3*mq1*PsR - BvgR*DLc*m3**4*mq1*PsR - Bvm1R*DLc*m3**5*mq1*PsR + BvgR*DLc*m1**2*m4**2*mq1*PsR + 
     &           BvgR*DLc*m2**2*m4**2*mq1*PsR + 8*Bvm1R*DLc*m1**2*m3*m4**2*mq1*PsR - Bvm2R*DLc*m1**2*m3*m4**2*mq1*PsR + 4*Bvm1R*DLc*m2**2*m3*m4**2*mq1*PsR + 
     &           Bvm2R*DLc*m2**2*m3*m4**2*mq1*PsR + 2*BvgR*DLc*m3**2*m4**2*mq1*PsR + 7*Bvm1R*DLc*m3**3*m4**2*mq1*PsR + Bvm2R*DLc*m3**3*m4**2*mq1*PsR + 
     &           BvgR*DLc*m4**4*mq1*PsR + 6*Bvm1R*DLc*m3*m4**4*mq1*PsR - Bvm2R*DLc*m3*m4**4*mq1*PsR + 2*Bsm2R*DRc*m1*m3*m4**2*mq3c*PsR - 
     &           2*Bvm2L*DRc*m1*m3**2*m4**2*mq3c*PsR - 2*Bvm2L*DRc*m1*m4**4*mq3c*PsR + 4*BvgR*DRc*m1*m3**2*mq1*mq3c*PsR + 4*Bvm1R*DRc*m1*m3**3*mq1*mq3c*PsR - 
     &           2*Bvm1R*DRc*m1*m3*m4**2*mq1*mq3c*PsR - 2*Bvm2R*DRc*m1*m3*m4**2*mq1*mq3c*PsR + 
     &           2*(Bsm1L*DLc*m1**2*m2**2*m3**2 + 2*BvgR*DLc*m1**2*m2**2*m3**2 + 2*Bvm1R*DLc*m1**2*m2**2*m3**3 - Bsm1L*DLc*m1**2*m3**4 - 
     &              Bsm2L*DLc*m1**2*m2**2*m4**2 - 2*Bvm2R*DLc*m1**2*m2**2*m3*m4**2 + Bsm2L*DLc*m1**2*m3**2*m4**2 + 2*BvgR*DLc*m1**2*m3**2*m4**2 + 
     &              BvgR*DLc*m2**2*m3**2*m4**2 + 4*Bvm1R*DLc*m1**2*m3**3*m4**2 + 2*Bvm1R*DLc*m2**2*m3**3*m4**2 + BvgR*DLc*m3**4*m4**2 + 
     &              2*Bvm1R*DLc*m3**5*m4**2 + BvgR*DLc*m2**2*m4**4 + 2*Bvm1R*DLc*m1**2*m3*m4**4 - 2*Bvm2R*DLc*m1**2*m3*m4**4 + Bvm1R*DLc*m2**2*m3*m4**4 - 
     &              Bvm2R*DLc*m2**2*m3*m4**4 + BvgR*DLc*m3**2*m4**4 + 3*Bvm1R*DLc*m3**3*m4**4 - Bvm2R*DLc*m3**3*m4**4 - 2*Bvm2R*DLc*m3*m4**6 - 
     &              BsgL*DLc*m3*(m1**2 + m4**2)*(m2**2 + m3**2 + 2*m4**2) + Bvm2L*DLc*m1**2*m2**2*m4**2*mq1 - Bvm2L*DLc*m1**2*m3**2*m4**2*mq1 + 
     &              DRc*m1*(m3**2 + m4**2)*(Bsm1L*m3**2 + Bvm1R*m3**3 - Bsm2L*m4**2 - 2*Bvm1R*m3*m4**2 - Bvm2R*m3*m4**2 + BvgR*(m3 - m4)*(m3 + m4) + 
     &                 Bvm2L*m4**2*mq1)*mq3c)*PvL - 2*BvgL*DLc*m3*(m1**2 + m4**2)*(m2**2 + m3**2 + 2*m4**2)*(PsR - mq1*PvL) + 
     &           2*Bvm1L*m1*m3**2*(DLc*m1*(m2 - m3)*(m2 + m3) + DRc*(m3**2 + m4**2)*mq3c)*(PsR - mq1*PvL))*sman + 
     &        BsgR*PsR*(2*DLc*(2*m1**2*m3**2*(m2**2 + m4**2) + (m2**2 + m3**2)*m4**2*(m3**2 + m4**2)) - 
     &           DLc*(-(m3**2*(m2**2 + m3**2)) + (m2**2 + 2*m3**2)*m4**2 + m4**4 + m1**2*(2*m2**2 - m3**2 + m4**2))*sman + 
     &           2*DRc*m1*mq3c*(m3**4 - m4**4 - 2*m3**2*sman)) + 
     &        Bsm1R*m3*PsR*(2*DLc*(2*m1**2*m2**2*m3**2 + 2*m3**2*(2*m1**2 + m2**2 + m3**2)*m4**2 + (2*m1**2 + m2**2 + 3*m3**2)*m4**4) - 
     &           DLc*(-(m3**2*(m2**2 + m3**2)) + (4*m2**2 + 7*m3**2)*m4**2 + 6*m4**4 + m1**2*(4*m2**2 + m3**2 + 8*m4**2))*sman + 
     &           2*DRc*m1*mq3c*(m3**4 - 2*m4**4 + m4**2*sman - m3**2*(m4**2 + 2*sman))))*xv + 
     &     DLc*tman*(2*Bsm2R*m1**2*m3*m4**2*PsR + Bsm2R*m3**3*m4**2*PsR + Bsm2R*m3*m4**4*PsR + 2*BvgR*m1**2*m3**2*mq1*PsR + 2*Bvm1R*m1**2*m3**3*mq1*PsR + 
     &        BvgR*m3**4*mq1*PsR + Bvm1R*m3**5*mq1*PsR - 2*BvgR*m1**2*m4**2*mq1*PsR - 2*Bvm2R*m1**2*m3*m4**2*mq1*PsR + Bvm1R*m3**3*m4**2*mq1*PsR - 
     &        Bvm2R*m3**3*m4**2*mq1*PsR - BvgR*m4**4*mq1*PsR - Bvm2R*m3*m4**4*mq1*PsR + 2*BvgL*m1**2*m3*PsR*sman - Bvm1L*m1**2*m3**2*PsR*sman - 
     &        Bvm1L*m3**4*PsR*sman + Bvm2L*m1**2*m4**2*PsR*sman + 2*BvgL*m3*m4**2*PsR*sman + Bvm2L*m3**2*m4**2*PsR*sman - 2*Bvm1R*m1**2*m3*mq1*PsR*sman - 
     &        BvgR*m3**2*mq1*PsR*sman - Bvm1R*m3**3*mq1*PsR*sman + BvgR*m4**2*mq1*PsR*sman - Bvm1R*m3*m4**2*mq1*PsR*sman + 2*BsgL*m1**2*m3*PvL*sman - 
     &        Bsm1L*m1**2*m3**2*PvL*sman - 2*BvgR*m1**2*m3**2*PvL*sman - 2*Bvm1R*m1**2*m3**3*PvL*sman - Bsm1L*m3**4*PvL*sman - BvgR*m3**4*PvL*sman - 
     &        Bvm1R*m3**5*PvL*sman + Bsm2L*m1**2*m4**2*PvL*sman + 2*BvgR*m1**2*m4**2*PvL*sman + 2*BsgL*m3*m4**2*PvL*sman + 2*Bvm2R*m1**2*m3*m4**2*PvL*sman + 
     &        Bsm2L*m3**2*m4**2*PvL*sman - Bvm1R*m3**3*m4**2*PvL*sman + Bvm2R*m3**3*m4**2*PvL*sman + BvgR*m4**4*PvL*sman + Bvm2R*m3*m4**4*PvL*sman - 
     &        2*BvgL*m1**2*m3*mq1*PvL*sman + Bvm1L*m1**2*m3**2*mq1*PvL*sman + Bvm1L*m3**4*mq1*PvL*sman - Bvm2L*m1**2*m4**2*mq1*PvL*sman - 
     &        2*BvgL*m3*m4**2*mq1*PvL*sman - Bvm2L*m3**2*m4**2*mq1*PvL*sman + Bvm1L*m1**2*PsR*sman**2 + Bvm1L*m3**2*PsR*sman**2 + Bsm1L*m1**2*PvL*sman**2 + 
     &        2*Bvm1R*m1**2*m3*PvL*sman**2 + Bsm1L*m3**2*PvL*sman**2 + BvgR*m3**2*PvL*sman**2 + Bvm1R*m3**3*PvL*sman**2 - BvgR*m4**2*PvL*sman**2 + 
     &        Bvm1R*m3*m4**2*PvL*sman**2 - Bvm1L*m1**2*mq1*PvL*sman**2 - Bvm1L*m3**2*mq1*PvL*sman**2 - 
     &        Bsm1R*m3*PsR*(m3**2 - sman)*(2*m1**2 + m3**2 + m4**2 - tman) - 
     &        (Bsm2R*m3*m4**2*PsR + (m3**2*(BvgR + Bvm1R*m3) - (BvgR + Bvm2R*m3)*m4**2)*mq1*PsR - 
     &           ((BvgR + Bvm1R*m3)*mq1*PsR + (m3**2*(Bsm1L + BvgR + Bvm1R*m3) - (Bsm2L + BvgR + Bvm2R*m3)*m4**2)*PvL + Bvm1L*m3**2*(PsR - mq1*PvL) + 
     &              Bvm2L*m4**2*(-PsR + mq1*PvL))*sman + (Bvm1L*PsR + (Bsm1L + BvgR + Bvm1R*m3 - Bvm1L*mq1)*PvL)*sman**2)*tman + 
     &        BsgR*PsR*(-((m3 - m4)*(m3 + m4)*(2*m1**2 + m3**2 + m4**2 - sman)) + (m3**2 - m4**2 - sman)*tman))*xv - 
     &     (-(Bsm1R*DLc*m1**2*m3*PsR*sman**2) - 2*BvgL*DLc*m1**2*m3*PsR*sman**2 + Bsm1R*DLc*m2**2*m3*PsR*sman**2 + 2*Bsm1R*DLc*m3**3*PsR*sman**2 + 
     &        Bvm2L*DLc*m1**2*m4**2*PsR*sman**2 + Bvm2L*DLc*m2**2*m4**2*PsR*sman**2 - 3*Bsm1R*DLc*m3*m4**2*PsR*sman**2 - Bsm2R*DLc*m3*m4**2*PsR*sman**2 - 
     &        2*BvgL*DLc*m3*m4**2*PsR*sman**2 + Bvm2L*DLc*m3**2*m4**2*PsR*sman**2 + Bvm2L*DLc*m4**4*PsR*sman**2 - BvgR*DLc*m1**2*mq1*PsR*sman**2 - 
     &        BvgR*DLc*m2**2*mq1*PsR*sman**2 + Bvm1R*DLc*m1**2*m3*mq1*PsR*sman**2 - Bvm1R*DLc*m2**2*m3*mq1*PsR*sman**2 - 2*BvgR*DLc*m3**2*mq1*PsR*sman**2 - 
     &        2*Bvm1R*DLc*m3**3*mq1*PsR*sman**2 + 3*Bvm1R*DLc*m3*m4**2*mq1*PsR*sman**2 + Bvm2R*DLc*m3*m4**2*mq1*PsR*sman**2 - 
     &        2*Bsm1R*DRc*m1*m3*mq3c*PsR*sman**2 - 2*Bvm2L*DRc*m1*m4**2*mq3c*PsR*sman**2 + 2*BvgR*DRc*m1*mq1*mq3c*PsR*sman**2 + 
     &        2*Bvm1R*DRc*m1*m3*mq1*mq3c*PsR*sman**2 + 2*Bsm1L*DLc*m1**2*m2**2*PvL*sman**2 + 2*BvgR*DLc*m1**2*m2**2*PvL*sman**2 - 
     &        2*BsgL*DLc*m1**2*m3*PvL*sman**2 + 4*Bvm1R*DLc*m1**2*m2**2*m3*PvL*sman**2 - 3*Bsm1L*DLc*m1**2*m3**2*PvL*sman**2 - 
     &        BvgR*DLc*m1**2*m3**2*PvL*sman**2 - Bsm1L*DLc*m2**2*m3**2*PvL*sman**2 - BvgR*DLc*m2**2*m3**2*PvL*sman**2 + Bvm1R*DLc*m1**2*m3**3*PvL*sman**2 - 
     &        Bvm1R*DLc*m2**2*m3**3*PvL*sman**2 - Bsm1L*DLc*m3**4*PvL*sman**2 - BvgR*DLc*m3**4*PvL*sman**2 - Bvm1R*DLc*m3**5*PvL*sman**2 + 
     &        Bsm2L*DLc*m1**2*m4**2*PvL*sman**2 + BvgR*DLc*m1**2*m4**2*PvL*sman**2 + Bsm2L*DLc*m2**2*m4**2*PvL*sman**2 + BvgR*DLc*m2**2*m4**2*PvL*sman**2 - 
     &        2*BsgL*DLc*m3*m4**2*PvL*sman**2 + 8*Bvm1R*DLc*m1**2*m3*m4**2*PvL*sman**2 - Bvm2R*DLc*m1**2*m3*m4**2*PvL*sman**2 + 
     &        4*Bvm1R*DLc*m2**2*m3*m4**2*PvL*sman**2 + Bvm2R*DLc*m2**2*m3*m4**2*PvL*sman**2 - Bsm1L*DLc*m3**2*m4**2*PvL*sman**2 + 
     &        Bsm2L*DLc*m3**2*m4**2*PvL*sman**2 + 2*BvgR*DLc*m3**2*m4**2*PvL*sman**2 + 7*Bvm1R*DLc*m3**3*m4**2*PvL*sman**2 + 
     &        Bvm2R*DLc*m3**3*m4**2*PvL*sman**2 + Bsm2L*DLc*m4**4*PvL*sman**2 + BvgR*DLc*m4**4*PvL*sman**2 + 6*Bvm1R*DLc*m3*m4**4*PvL*sman**2 - 
     &        Bvm2R*DLc*m3*m4**4*PvL*sman**2 + 2*BvgL*DLc*m1**2*m3*mq1*PvL*sman**2 - Bvm2L*DLc*m1**2*m4**2*mq1*PvL*sman**2 - 
     &        Bvm2L*DLc*m2**2*m4**2*mq1*PvL*sman**2 + 2*BvgL*DLc*m3*m4**2*mq1*PvL*sman**2 - Bvm2L*DLc*m3**2*m4**2*mq1*PvL*sman**2 - 
     &        Bvm2L*DLc*m4**4*mq1*PvL*sman**2 + 4*Bsm1L*DRc*m1*m3**2*mq3c*PvL*sman**2 + 4*BvgR*DRc*m1*m3**2*mq3c*PvL*sman**2 + 
     &        4*Bvm1R*DRc*m1*m3**3*mq3c*PvL*sman**2 + 2*Bsm1L*DRc*m1*m4**2*mq3c*PvL*sman**2 - 2*Bsm2L*DRc*m1*m4**2*mq3c*PvL*sman**2 - 
     &        2*Bvm1R*DRc*m1*m3*m4**2*mq3c*PvL*sman**2 - 2*Bvm2R*DRc*m1*m3*m4**2*mq3c*PvL*sman**2 + 2*Bvm2L*DRc*m1*m4**2*mq1*mq3c*PvL*sman**2 - 
     &        Bsm1R*DLc*m3*PsR*sman**3 - Bvm2L*DLc*m4**2*PsR*sman**3 + BvgR*DLc*mq1*PsR*sman**3 + Bvm1R*DLc*m3*mq1*PsR*sman**3 + Bsm1L*DLc*m1**2*PvL*sman**3 + 
     &        BvgR*DLc*m1**2*PvL*sman**3 + Bsm1L*DLc*m2**2*PvL*sman**3 + BvgR*DLc*m2**2*PvL*sman**3 - Bvm1R*DLc*m1**2*m3*PvL*sman**3 + 
     &        Bvm1R*DLc*m2**2*m3*PvL*sman**3 + 2*Bsm1L*DLc*m3**2*PvL*sman**3 + 2*BvgR*DLc*m3**2*PvL*sman**3 + 2*Bvm1R*DLc*m3**3*PvL*sman**3 + 
     &        Bsm1L*DLc*m4**2*PvL*sman**3 - Bsm2L*DLc*m4**2*PvL*sman**3 - 3*Bvm1R*DLc*m3*m4**2*PvL*sman**3 - Bvm2R*DLc*m3*m4**2*PvL*sman**3 + 
     &        Bvm2L*DLc*m4**2*mq1*PvL*sman**3 - 2*Bsm1L*DRc*m1*mq3c*PvL*sman**3 - 2*BvgR*DRc*m1*mq3c*PvL*sman**3 - 2*Bvm1R*DRc*m1*m3*mq3c*PvL*sman**3 - 
     &        Bsm1L*DLc*PvL*sman**4 - BvgR*DLc*PvL*sman**4 - Bvm1R*DLc*m3*PvL*sman**4 + 
     &        DLc*(-(Bsm2R*m3*m4**2*(m1**2 + m2**2 + 2*m4**2)*PsR) - 
     &           ((m1**2 + m2**2)*m3**2*(BvgR + Bvm1R*m3) + (-(Bvm2R*(m1**2 + m2**2)*m3) + 4*Bvm1R*m3**3 + BvgR*(m1**2 - m2**2 + 2*m3**2))*m4**2 + 
     &              2*(Bvm1R - Bvm2R)*m3*m4**4)*mq1*PsR + Bsm1R*m3*PsR*(2*m4**2*(2*m3**2 + m4**2 - 2*sman) + m1**2*(m3**2 - sman) + m2**2*(m3**2 - sman)) + 
     &           (-(Bvm2L*m4**2*(m2**2 + m4**2)*PsR) + ((m1**2 + m2**2)*(BvgR + Bvm1R*m3) + 2*(BvgR + 2*Bvm1R*m3)*m4**2)*mq1*PsR + 
     &              (m3**2*(Bsm1L*m2**2 + (m1**2 + m2**2)*(BvgR + Bvm1R*m3)) + 
     &                 (-(Bsm2L*m2**2) + BvgR*(m1 - m2)*(m1 + m2) - Bvm2R*(m1**2 + m2**2)*m3 + (Bsm1L + 2*BvgR)*m3**2)*m4**2)*PvL)*sman)*tman + 
     &        BsgR*PsR*((-2*DRc*m1*mq3c + DLc*(m1**2 + m2**2 + 2*m3**2 - sman))*sman**2 + 
     &           DLc*(2*m4**2*(m3**2 - sman) + m2**2*(m3**2 - m4**2 - sman) + m1**2*(m3**2 + m4**2 - sman))*tman) + 
     &        Bvm1L*sman*(2*DRc*m1*mq3c*(PsR - mq1*PvL)*(2*m3**2 + m4**2 - sman)*sman + 
     &           DLc*(-((PsR - mq1*PvL)*(m3**2 - sman)*(m2**2 + m3**2 + m4**2 - sman)*sman) + m1**2*(PsR - mq1*PvL)*sman*(2*m2**2 - 3*m3**2 + sman) + 
     &              m3**2*(m2**2 + m4**2)*PsR*tman)))*xv + DLc*(Bsm2L*m4**4*PvL*sman*tman + 2*Bvm2R*m3*m4**4*PvL*sman*tman + 
     &        Bvm1L*m2**2*m3**2*mq1*PvL*sman*tman - Bvm2L*m2**2*m4**2*mq1*PvL*sman*tman + Bvm1L*m3**2*m4**2*mq1*PvL*sman*tman - 
     &        Bvm2L*m4**4*mq1*PvL*sman*tman + Bvm1L*m2**2*PsR*sman**2*tman + Bvm1L*m4**2*PsR*sman**2*tman + BvgR*m1**2*PvL*sman**2*tman + 
     &        Bsm1L*m2**2*PvL*sman**2*tman + BvgR*m2**2*PvL*sman**2*tman + Bsm1L*m4**2*PvL*sman**2*tman + 2*BvgR*m4**2*PvL*sman**2*tman - 
     &        Bvm1L*m2**2*mq1*PvL*sman**2*tman - Bvm1L*m4**2*mq1*PvL*sman**2*tman + 
     &        (Bsm2R*m1**2*m3*m4**2*PsR + 3*Bsm2R*m2**2*m3*m4**2*PsR + Bsm2R*m3**3*m4**2*PsR + 5*Bsm2R*m3*m4**4*PsR - 
     &           BsgR*(m3**2*(m1**2 + 3*m2**2 + m3**2) + (m1**2 + m2**2 + 4*m3**2)*m4**2 + m4**4)*PsR + BvgR*m1**2*m3**2*mq1*PsR + 
     &           3*BvgR*m2**2*m3**2*mq1*PsR + BvgR*m3**4*mq1*PsR + BvgR*m1**2*m4**2*mq1*PsR + BvgR*m2**2*m4**2*mq1*PsR - Bvm2R*m1**2*m3*m4**2*mq1*PsR - 
     &           3*Bvm2R*m2**2*m3*m4**2*mq1*PsR + 4*BvgR*m3**2*m4**2*mq1*PsR - Bvm2R*m3**3*m4**2*mq1*PsR + BvgR*m4**4*mq1*PsR - 5*Bvm2R*m3*m4**4*mq1*PsR + 
     &           2*BvgL*m2**2*m3*PsR*sman - Bvm1L*m1**2*m3**2*PsR*sman - Bvm1L*m2**2*m3**2*PsR*sman + 2*BvgL*m3**3*PsR*sman + Bvm1L*m3**4*PsR*sman + 
     &           Bvm2L*m1**2*m4**2*PsR*sman + Bvm2L*m2**2*m4**2*PsR*sman - 2*Bsm2R*m3*m4**2*PsR*sman + 4*BvgL*m3*m4**2*PsR*sman - Bvm1L*m3**2*m4**2*PsR*sman - 
     &           Bvm2L*m3**2*m4**2*PsR*sman + Bvm2L*m4**4*PsR*sman + BsgR*(m1**2 + m2**2 + m3**2 + m4**2)*PsR*sman - BvgR*m1**2*mq1*PsR*sman - 
     &           BvgR*m2**2*mq1*PsR*sman - BvgR*m3**2*mq1*PsR*sman - BvgR*m4**2*mq1*PsR*sman + 2*Bvm2R*m3*m4**2*mq1*PsR*sman + 2*BsgL*m2**2*m3*PvL*sman - 
     &           Bsm1L*m1**2*m3**2*PvL*sman - BvgR*m1**2*m3**2*PvL*sman - Bsm1L*m2**2*m3**2*PvL*sman - 3*BvgR*m2**2*m3**2*PvL*sman + 2*BsgL*m3**3*PvL*sman + 
     &           Bsm1L*m3**4*PvL*sman - BvgR*m3**4*PvL*sman + Bsm2L*m1**2*m4**2*PvL*sman - BvgR*m1**2*m4**2*PvL*sman + Bsm2L*m2**2*m4**2*PvL*sman - 
     &           BvgR*m2**2*m4**2*PvL*sman + 4*BsgL*m3*m4**2*PvL*sman + Bvm2R*m1**2*m3*m4**2*PvL*sman + 3*Bvm2R*m2**2*m3*m4**2*PvL*sman - 
     &           Bsm1L*m3**2*m4**2*PvL*sman - Bsm2L*m3**2*m4**2*PvL*sman - 4*BvgR*m3**2*m4**2*PvL*sman + Bvm2R*m3**3*m4**2*PvL*sman + Bsm2L*m4**4*PvL*sman - 
     &           BvgR*m4**4*PvL*sman + 5*Bvm2R*m3*m4**4*PvL*sman - 2*BvgL*m2**2*m3*mq1*PvL*sman + Bvm1L*m1**2*m3**2*mq1*PvL*sman + 
     &           Bvm1L*m2**2*m3**2*mq1*PvL*sman - 2*BvgL*m3**3*mq1*PvL*sman - Bvm1L*m3**4*mq1*PvL*sman - Bvm2L*m1**2*m4**2*mq1*PvL*sman - 
     &           Bvm2L*m2**2*m4**2*mq1*PvL*sman - 4*BvgL*m3*m4**2*mq1*PvL*sman + Bvm1L*m3**2*m4**2*mq1*PvL*sman + Bvm2L*m3**2*m4**2*mq1*PvL*sman - 
     &           Bvm2L*m4**4*mq1*PvL*sman + Bvm1L*m1**2*PsR*sman**2 + Bvm1L*m2**2*PsR*sman**2 - 2*BvgL*m3*PsR*sman**2 - Bvm1L*m3**2*PsR*sman**2 + 
     &           Bvm1L*m4**2*PsR*sman**2 + Bsm1L*m1**2*PvL*sman**2 + BvgR*m1**2*PvL*sman**2 + Bsm1L*m2**2*PvL*sman**2 + BvgR*m2**2*PvL*sman**2 - 
     &           2*BsgL*m3*PvL*sman**2 - Bsm1L*m3**2*PvL*sman**2 + BvgR*m3**2*PvL*sman**2 + Bsm1L*m4**2*PvL*sman**2 + BvgR*m4**2*PvL*sman**2 - 
     &           2*Bvm2R*m3*m4**2*PvL*sman**2 - Bvm1L*m1**2*mq1*PvL*sman**2 - Bvm1L*m2**2*mq1*PvL*sman**2 + 2*BvgL*m3*mq1*PvL*sman**2 + 
     &           Bvm1L*m3**2*mq1*PvL*sman**2 - Bvm1L*m4**2*mq1*PvL*sman**2 - Bsm1R*m3*PsR*(m3**2 - sman)*(m1**2 + 3*m2**2 + m3**2 + 5*m4**2 - 2*sman - tman) + 
     &           BsgR*PsR*(m3**2 + m4**2 - sman)*tman - (Bsm2R*m3*m4**2*PsR + (-(Bvm2R*m3*m4**2) + BvgR*(m3**2 + m4**2 - sman))*(mq1*PsR - PvL*sman))*tman)*
     &         uman + ((Bsm1R*m3**3 - Bsm2R*m3*m4**2 + BsgR*(m3**2 + m4**2) - BvgR*m3**2*mq1 - BvgR*m4**2*mq1 + Bvm2R*m3*m4**2*mq1)*PsR + 
     &           (-((BsgR + m3*(Bsm1R - Bvm1L*m3) + Bvm2L*m4**2 - BvgR*mq1)*PsR) + 
     &              (Bsm1L*m3**2 - (Bsm2L + Bvm2R*m3)*m4**2 + BvgR*(m3**2 + m4**2) - Bvm1L*m3**2*mq1 + Bvm2L*m4**2*mq1)*PvL)*sman - 
     &           (Bvm1L*PsR + (Bsm1L + BvgR - Bvm1L*mq1)*PvL)*sman**2)*uman**2 + 
     &        Bvm1R*m3*(PvL*sman*(-2*m4**2*(2*m3**2 + m4**2 - 2*sman) + (m1**2 + m2**2)*sman)*tman + 
     &           (m3**2 - sman)*(mq1*PsR - PvL*sman)*(m1**2 + 3*m2**2 + m3**2 + 5*m4**2 - 2*sman - tman)*uman + (m3**2 - sman)*(-(mq1*PsR) + PvL*sman)*uman**2)
     &        )*xv - (2*m3*m4**2*((Bvm1R + Bvm2R)*m4**2*mq1*PsR - (BvgL*PsR + (BsgL + (Bvm1R + Bvm2R)*m4**2 - BvgL*mq1)*PvL)*sman)*
     &         (2*DRc*m1*mq3c + DLc*(-2*m1**2 - m2**2 - m3**2 + sman + tman)) + 
     &        DLc*m3*(2*(BvgL*PsR + BsgL*PvL - BvgL*mq1*PvL)*sman + Bvm1R*(m3**2 - sman)*(mq1*PsR - PvL*sman) + Bvm2R*m4**2*(-(mq1*PsR) + PvL*sman))*tman*
     &         uman + Bsm2R*m3*m4**2*PsR*(2*m4**2*(-2*DRc*m1*mq3c + DLc*(2*m1**2 + m2**2 + m3**2 - sman - tman)) + DLc*tman*uman) + 
     &        BvgR*(mq1*PsR - PvL*sman)*(2*DRc*m1*m4**2*mq3c*(-m3**2 + m4**2 - sman) + 
     &           DLc*(m4**2*(m3**2 + m4**2 - sman - tman)*(2*m3**2 + sman - tman) + (m4**4 - m4**2*tman + sman*tman + m3**2*(-m4**2 + tman))*uman - 
     &              m4**2*uman**2 + m1**2*m4**2*(4*m3**2 + sman - tman + uman) + m2**2*m4**2*(2*m3**2 - 2*m4**2 + sman - tman + uman))))*xv + 
     &     BsgR*PsR*(-4*DLc*((m3**2 + m4**2 - sman - tman)*(2*m3**2 + sman - tman) + (-m3**2 + m4**2)*uman - uman**2 + m1**2*(4*m3**2 + sman - tman + uman) + 
     &           m2**2*(2*m3**2 - 2*m4**2 + sman - tman + uman)) + 
     &        DLc*(m4**2*(m3**2 + m4**2 - sman - tman)*(2*m3**2 + sman - tman) + (m4**4 - m4**2*tman + sman*tman + m3**2*(-m4**2 + tman))*uman - 
     &           m4**2*uman**2 + m1**2*m4**2*(4*m3**2 + sman - tman + uman) + m2**2*m4**2*(2*m3**2 - 2*m4**2 + sman - tman + uman))*xv - 
     &        2*DRc*m1*mq3c*(m3**2 - m4**2 + sman)*(-4 + m4**2*xv)) + 
     &     Bsm1R*m3*PsR*(2*DLc*m4**2*(2*m1**2 + m2**2 + m3**2 - sman - tman) + 
     &        DLc*(-2*m4**4*(2*m1**2 + m2**2 + m3**2 - sman - tman) + (m3**2 - sman)*tman*uman)*xv + 4*DRc*m1*m4**2*mq3c*(-1 + m4**2*xv))) + 
     &  A1vR*A3Lc*(2*(Bsm2L*DLc*m1*m2**2*m3**2*PsL + Bsm2L*DLc*m1*m3**4*PsL + 2*Bsm2L*DLc*m1*m3**2*m4**2*PsL - 4*BvgL*DLc*m1*m2**2*m3*mq1*PsL + 
     &        Bvm1L*DLc*m1*m2**2*m3**2*mq1*PsL - Bvm2L*DLc*m1*m2**2*m3**2*mq1*PsL + Bvm1L*DLc*m1*m3**4*mq1*PsL - Bvm2L*DLc*m1*m3**4*mq1*PsL + 
     &        2*Bvm1L*DLc*m1*m2**2*m4**2*mq1*PsL - 4*BvgL*DLc*m1*m3*m4**2*mq1*PsL - 2*Bvm2L*DLc*m1*m3**2*m4**2*mq1*PsL + Bsm2L*DRc*m1**2*m4**2*mq3c*PsL + 
     &        Bsm2L*DRc*m4**4*mq3c*PsL + 4*Bvm1L*DRc*m1**2*m3**2*mq1*mq3c*PsL + 2*Bvm1L*DRc*m3**4*mq1*mq3c*PsL - Bvm1L*DRc*m1**2*m4**2*mq1*mq3c*PsL - 
     &        Bvm2L*DRc*m1**2*m4**2*mq1*mq3c*PsL + 2*Bvm1L*DRc*m3**2*m4**2*mq1*mq3c*PsL - Bvm1L*DRc*m4**4*mq1*mq3c*PsL - Bvm2L*DRc*m4**4*mq1*mq3c*PsL + 
     &        3*Bvm1R*DLc*m1*m2**2*m3*PsL*sman + Bvm2R*DLc*m1*m2**2*m3*PsL*sman - Bsm2L*DLc*m1*m3**2*PsL*sman - Bvm1R*DLc*m1*m3**3*PsL*sman + 
     &        Bvm2R*DLc*m1*m3**3*PsL*sman + 2*Bvm1R*DLc*m1*m3*m4**2*PsL*sman + 2*Bvm2R*DLc*m1*m3*m4**2*PsL*sman - 2*Bvm1L*DLc*m1*m2**2*mq1*PsL*sman + 
     &        2*BvgL*DLc*m1*m3*mq1*PsL*sman + Bvm1L*DLc*m1*m3**2*mq1*PsL*sman + Bvm2L*DLc*m1*m3**2*mq1*PsL*sman - Bsm2L*DRc*m1**2*mq3c*PsL*sman - 
     &        2*BvgR*DRc*m1**2*mq3c*PsL*sman - 3*Bvm1R*DRc*m1**2*m3*mq3c*PsL*sman - Bvm2R*DRc*m1**2*m3*mq3c*PsL*sman - 2*BvgR*DRc*m3**2*mq3c*PsL*sman - 
     &        2*Bvm1R*DRc*m3**3*mq3c*PsL*sman - Bsm2L*DRc*m4**2*mq3c*PsL*sman - Bvm1R*DRc*m3*m4**2*mq3c*PsL*sman - Bvm2R*DRc*m3*m4**2*mq3c*PsL*sman + 
     &        Bvm1L*DRc*m1**2*mq1*mq3c*PsL*sman + Bvm2L*DRc*m1**2*mq1*mq3c*PsL*sman + Bvm1L*DRc*m4**2*mq1*mq3c*PsL*sman + Bvm2L*DRc*m4**2*mq1*mq3c*PsL*sman + 
     &        3*Bsm1R*DLc*m1*m2**2*m3*PvR*sman + Bsm2R*DLc*m1*m2**2*m3*PvR*sman + 4*BvgL*DLc*m1*m2**2*m3*PvR*sman - Bvm1L*DLc*m1*m2**2*m3**2*PvR*sman + 
     &        Bvm2L*DLc*m1*m2**2*m3**2*PvR*sman - Bsm1R*DLc*m1*m3**3*PvR*sman + Bsm2R*DLc*m1*m3**3*PvR*sman - Bvm1L*DLc*m1*m3**4*PvR*sman + 
     &        Bvm2L*DLc*m1*m3**4*PvR*sman - 2*Bvm1L*DLc*m1*m2**2*m4**2*PvR*sman + 2*Bsm1R*DLc*m1*m3*m4**2*PvR*sman + 2*Bsm2R*DLc*m1*m3*m4**2*PvR*sman + 
     &        4*BvgL*DLc*m1*m3*m4**2*PvR*sman + 2*Bvm2L*DLc*m1*m3**2*m4**2*PvR*sman - 3*Bvm1R*DLc*m1*m2**2*m3*mq1*PvR*sman - 
     &        Bvm2R*DLc*m1*m2**2*m3*mq1*PvR*sman + Bvm1R*DLc*m1*m3**3*mq1*PvR*sman - Bvm2R*DLc*m1*m3**3*mq1*PvR*sman - 2*Bvm1R*DLc*m1*m3*m4**2*mq1*PvR*sman - 
     &        2*Bvm2R*DLc*m1*m3*m4**2*mq1*PvR*sman - 2*BsgR*DRc*m1**2*mq3c*PvR*sman - 3*Bsm1R*DRc*m1**2*m3*mq3c*PvR*sman - Bsm2R*DRc*m1**2*m3*mq3c*PvR*sman - 
     &        2*BsgR*DRc*m3**2*mq3c*PvR*sman - 4*Bvm1L*DRc*m1**2*m3**2*mq3c*PvR*sman - 2*Bsm1R*DRc*m3**3*mq3c*PvR*sman - 2*Bvm1L*DRc*m3**4*mq3c*PvR*sman + 
     &        Bvm1L*DRc*m1**2*m4**2*mq3c*PvR*sman + Bvm2L*DRc*m1**2*m4**2*mq3c*PvR*sman - Bsm1R*DRc*m3*m4**2*mq3c*PvR*sman - 
     &        Bsm2R*DRc*m3*m4**2*mq3c*PvR*sman - 2*Bvm1L*DRc*m3**2*m4**2*mq3c*PvR*sman + Bvm1L*DRc*m4**4*mq3c*PvR*sman + Bvm2L*DRc*m4**4*mq3c*PvR*sman + 
     &        2*BvgR*DRc*m1**2*mq1*mq3c*PvR*sman + 3*Bvm1R*DRc*m1**2*m3*mq1*mq3c*PvR*sman + Bvm2R*DRc*m1**2*m3*mq1*mq3c*PvR*sman + 
     &        2*BvgR*DRc*m3**2*mq1*mq3c*PvR*sman + 2*Bvm1R*DRc*m3**3*mq1*mq3c*PvR*sman + Bvm1R*DRc*m3*m4**2*mq1*mq3c*PvR*sman + 
     &        Bvm2R*DRc*m3*m4**2*mq1*mq3c*PvR*sman - Bvm1R*DLc*m1*m3*PsL*sman**2 - Bvm2R*DLc*m1*m3*PsL*sman**2 + 2*Bvm1L*DLc*m1*m2**2*PvR*sman**2 - 
     &        Bsm1R*DLc*m1*m3*PvR*sman**2 - Bsm2R*DLc*m1*m3*PvR*sman**2 - 2*BvgL*DLc*m1*m3*PvR*sman**2 - Bvm1L*DLc*m1*m3**2*PvR*sman**2 - 
     &        Bvm2L*DLc*m1*m3**2*PvR*sman**2 + Bvm1R*DLc*m1*m3*mq1*PvR*sman**2 + Bvm2R*DLc*m1*m3*mq1*PvR*sman**2 - Bvm1L*DRc*m1**2*mq3c*PvR*sman**2 - 
     &        Bvm2L*DRc*m1**2*mq3c*PvR*sman**2 - Bvm1L*DRc*m4**2*mq3c*PvR*sman**2 - Bvm2L*DRc*m4**2*mq3c*PvR*sman**2 - 
     &        Bsm1L*PsL*(DLc*m1*(m2**2*(m3**2 + 2*m4**2 - 2*sman) + m3**2*(m3**2 + sman)) + 
     &           DRc*mq3c*(2*m3**4 + 2*m3**2*m4**2 - m4**4 + m4**2*sman + m1**2*(4*m3**2 - m4**2 + sman))) + 
     &        2*BsgL*DLc*m1*m3*PsL*(2*m2**2 + 2*m4**2 - sman - tman) + Bsm1L*m3**2*(DLc*m1 + 2*DRc*mq3c)*PsL*tman - 
     &        m3*(DLc*m1*(-2*BvgL*mq1 + m3*(Bsm2L + (Bvm1L - Bvm2L)*mq1)) + 2*Bvm1L*DRc*m3*mq1*mq3c)*PsL*tman) - 
     &     2*(Bsm2L*DLc*m1*m2**2*m3**2*m4**2*PsL + Bsm2L*DLc*m1*m3**4*m4**2*PsL + 2*Bsm2L*DLc*m1*m3**2*m4**4*PsL + BvgL*DLc*m1*m2**2*m3**3*mq1*PsL + 
     &        Bvm1L*DLc*m1*m2**2*m3**4*mq1*PsL + BvgL*DLc*m1*m3**5*mq1*PsL + Bvm1L*DLc*m1*m3**6*mq1*PsL - BvgL*DLc*m1*m2**2*m3*m4**2*mq1*PsL + 
     &        Bvm1L*DLc*m1*m2**2*m3**2*m4**2*mq1*PsL - Bvm2L*DLc*m1*m2**2*m3**2*m4**2*mq1*PsL + BvgL*DLc*m1*m3**3*m4**2*mq1*PsL + 
     &        Bvm1L*DLc*m1*m3**4*m4**2*mq1*PsL - Bvm2L*DLc*m1*m3**4*m4**2*mq1*PsL + Bvm1L*DLc*m1*m2**2*m4**4*mq1*PsL - 2*BvgL*DLc*m1*m3*m4**4*mq1*PsL - 
     &        Bvm1L*DLc*m1*m3**2*m4**4*mq1*PsL - 2*Bvm2L*DLc*m1*m3**2*m4**4*mq1*PsL + Bsm2L*DRc*m1**2*m4**4*mq3c*PsL + Bsm2L*DRc*m4**6*mq3c*PsL + 
     &        2*Bvm1L*DRc*m1**2*m3**2*m4**2*mq1*mq3c*PsL + Bvm1L*DRc*m3**4*m4**2*mq1*mq3c*PsL + Bvm1L*DRc*m1**2*m4**4*mq1*mq3c*PsL - 
     &        Bvm2L*DRc*m1**2*m4**4*mq1*mq3c*PsL + 2*Bvm1L*DRc*m3**2*m4**4*mq1*mq3c*PsL - Bvm2L*DRc*m4**6*mq1*mq3c*PsL - 
     &        (Bsm2L*DLc*m1*m3**2*m4**2*PsL + BvgL*DLc*m1*m2**2*m3*mq1*PsL + Bvm1L*DLc*m1*m2**2*m3**2*mq1*PsL + 2*BvgL*DLc*m1*m3**3*mq1*PsL + 
     &           2*Bvm1L*DLc*m1*m3**4*mq1*PsL + Bvm1L*DLc*m1*m2**2*m4**2*mq1*PsL + BvgL*DLc*m1*m3*m4**2*mq1*PsL + Bvm1L*DLc*m1*m3**2*m4**2*mq1*PsL - 
     &           Bvm2L*DLc*m1*m3**2*m4**2*mq1*PsL - BvgR*DRc*m1**2*m3**2*mq3c*PsL + Bsm2L*DRc*m1**2*m4**2*mq3c*PsL - BvgR*DRc*m1**2*m4**2*mq3c*PsL - 
     &           BvgR*DRc*m3**2*m4**2*mq3c*PsL + Bsm2L*DRc*m4**4*mq3c*PsL - BvgR*DRc*m4**4*mq3c*PsL + Bvm1L*DRc*m1**2*m3**2*mq1*mq3c*PsL + 
     &           2*Bvm1L*DRc*m1**2*m4**2*mq1*mq3c*PsL - Bvm2L*DRc*m1**2*m4**2*mq1*mq3c*PsL + 2*Bvm1L*DRc*m3**2*m4**2*mq1*mq3c*PsL + 
     &           Bvm1L*DRc*m4**4*mq1*mq3c*PsL - Bvm2L*DRc*m4**4*mq1*mq3c*PsL + 
     &           (BvgL*DLc*m1*m3*(m3 - m4)*(m3 + m4)*(m2**2 + m3**2 + 2*m4**2) + 
     &              DLc*m1*(Bvm1L*m3**4*(m2**2 + m3**2) - m3*(Bsm2R + (-Bvm1L + Bvm2L)*m3)*(m2**2 + m3**2)*m4**2 + 
     &                 (Bvm1L*(m2 - m3)*(m2 + m3) - 2*m3*(Bsm2R + Bvm2L*m3))*m4**4) + 
     &              DRc*(m3*(Bsm2R*m1**2 + Bvm1L*m3*(2*m1**2 + m3**2))*m4**2 + ((Bvm1L - Bvm2L)*m1**2 + Bsm2R*m3 + 2*Bvm1L*m3**2)*m4**4 - Bvm2L*m4**6 - 
     &                 BsgR*(m1**2 + m4**2)*(m3**2 + m4**2) + BvgR*(m3**2*m4**2 + m1**2*(m3**2 + m4**2))*mq1)*mq3c + 
     &              Bsm1R*m3**3*(DLc*m1*(m2**2 + m3**2 + 2*m4**2) - DRc*(m1**2 + m4**2)*mq3c))*PvR + 
     &           Bvm1R*m3**3*(DLc*m1*(m2**2 + m3**2 + 2*m4**2)*(PsL - mq1*PvR) + DRc*mq3c*(-((m1**2 + m4**2)*PsL) + m1**2*mq1*PvR)) + 
     &           Bvm2R*m3*m4**2*(-(DLc*m1*(m2**2 + m3**2 + 2*m4**2)*(PsL - mq1*PvR)) + DRc*mq3c*(m4**2*PsL + m1**2*(PsL - mq1*PvR))))*sman + 
     &        BsgL*DLc*m1*m3*PsL*(-m3**4 - m3**2*(m4**2 - 2*sman) + m2**2*(m4**2 + sman) + m4**2*(2*m4**2 + sman)) + 
     &        Bsm1L*PsL*(DRc*mq3c*(m4**2*(-m3**4 - 2*m3**2*(m4**2 - sman) + m4**2*sman) + m1**2*(-m4**4 + 2*m4**2*sman + m3**2*(-2*m4**2 + sman))) - 
     &           DLc*m1*(m2**2*(m3**4 + m4**4 + m3**2*(m4**2 - sman) - m4**2*sman) + m3**2*(m3**4 + m3**2*(m4**2 - 2*sman) - m4**2*(m4**2 + sman)))))*xv + 
     &     DLc*m1*(Bvm2L*m4**2*(-m3**2 + m4**2)*mq1*PsL + (2*BsgL*m3*PsL + m4**2*(-Bsm1L + Bsm2L + 2*Bvm2R*m3 - Bvm2L*mq1)*PsL + 
     &           (-2*Bsm1R*m3**3 + m4**2*(2*Bsm2R*m3 + Bvm2L*(m3 - m4)*(m3 + m4) - 2*Bvm2R*m3*mq1))*PvR - 2*Bvm1R*m3**3*(PsL - mq1*PvR))*sman + 
     &        (Bsm1L*PsL + (2*Bsm1R*m3 + Bvm2L*m4**2)*PvR + 2*Bvm1R*m3*(PsL - mq1*PvR))*sman**2 + 2*BvgL*m3*(m3**2 - m4**2 - sman)*(mq1*PsL - PvR*sman) + 
     &        Bvm1L*(m3**2 - sman)*(m3**2 - m4**2 + sman)*(mq1*PsL - PvR*sman))*tman*xv + 
     &     (-2*Bvm2R*DRc*m3*m4**4*mq1*mq3c*PvR*sman + 2*BsgL*DLc*m1*m3*PsL*sman**2 + 2*Bsm1L*DLc*m1*m3**2*PsL*sman**2 + 2*Bvm2R*DLc*m1*m3*m4**2*PsL*sman**2 - 
     &        2*BvgL*DLc*m1*m3*mq1*PsL*sman**2 - 2*Bvm1L*DLc*m1*m3**2*mq1*PsL*sman**2 + 2*Bsm1L*DRc*m1**2*mq3c*PsL*sman**2 + 
     &        2*Bsm1L*DRc*m4**2*mq3c*PsL*sman**2 - 2*Bvm1L*DRc*m1**2*mq1*mq3c*PsL*sman**2 - 2*Bvm1L*DRc*m4**2*mq1*mq3c*PsL*sman**2 - 
     &        2*Bsm1R*DLc*m1*m2**2*m3*PvR*sman**2 - 2*BvgL*DLc*m1*m2**2*m3*PvR*sman**2 - 2*Bvm1L*DLc*m1*m2**2*m3**2*PvR*sman**2 - 
     &        4*Bsm1R*DLc*m1*m3**3*PvR*sman**2 - 4*BvgL*DLc*m1*m3**3*PvR*sman**2 - 4*Bvm1L*DLc*m1*m3**4*PvR*sman**2 - 2*Bvm1L*DLc*m1*m2**2*m4**2*PvR*sman**2 - 
     &        4*Bsm1R*DLc*m1*m3*m4**2*PvR*sman**2 + 2*Bsm2R*DLc*m1*m3*m4**2*PvR*sman**2 - 2*BvgL*DLc*m1*m3*m4**2*PvR*sman**2 - 
     &        2*Bvm1L*DLc*m1*m3**2*m4**2*PvR*sman**2 + 2*Bvm2L*DLc*m1*m3**2*m4**2*PvR*sman**2 - 2*Bvm2R*DLc*m1*m3*m4**2*mq1*PvR*sman**2 + 
     &        2*BsgR*DRc*m1**2*mq3c*PvR*sman**2 + 2*Bsm1R*DRc*m1**2*m3*mq3c*PvR*sman**2 - 2*Bvm1L*DRc*m1**2*m3**2*mq3c*PvR*sman**2 + 
     &        2*BsgR*DRc*m4**2*mq3c*PvR*sman**2 - 4*Bvm1L*DRc*m1**2*m4**2*mq3c*PvR*sman**2 + 2*Bvm2L*DRc*m1**2*m4**2*mq3c*PvR*sman**2 + 
     &        2*Bsm1R*DRc*m3*m4**2*mq3c*PvR*sman**2 - 4*Bvm1L*DRc*m3**2*m4**2*mq3c*PvR*sman**2 - 2*Bvm1L*DRc*m4**4*mq3c*PvR*sman**2 + 
     &        2*Bvm2L*DRc*m4**4*mq3c*PvR*sman**2 + 2*Bsm1R*DLc*m1*m3*PvR*sman**3 + 2*BvgL*DLc*m1*m3*PvR*sman**3 + 2*Bvm1L*DLc*m1*m3**2*PvR*sman**3 + 
     &        2*Bvm1L*DRc*m1**2*mq3c*PvR*sman**3 + 2*Bvm1L*DRc*m4**2*mq3c*PvR*sman**3 + 
     &        2*BvgR*DRc*mq3c*sman*(m4**4*mq1*PvR + (m1**2 + m4**2)*(PsL - mq1*PvR)*sman) + 
     &        2*Bvm1R*m3*sman*(-(DLc*m1*(PsL - mq1*PvR)*(m2**2 + 2*(m3**2 + m4**2) - sman)*sman) + 
     &           DRc*mq3c*(m3**2*m4**2*mq1*PvR + (m1**2 + m4**2)*(PsL - mq1*PvR)*sman)) + 
     &        PsL*((m3**2 + m4**2)*(DLc*m1*(-(Bsm1L*m3**2) + Bvm1L*m3**2*mq1 + m4**2*(Bsm2L - Bvm2L*mq1)) - 2*DRc*m4**2*(Bsm1L - Bvm1L*mq1)*mq3c) + 
     &           (DLc*m1*(Bsm1L*(2*m3**2 + m4**2) - 2*Bvm1L*m3**2*mq1 - m4**2*(Bsm2L + (Bvm1L - Bvm2L)*mq1)) + 2*DRc*m4**2*(Bsm1L - Bvm1L*mq1)*mq3c)*sman)*tman
     &        )*xv - (Bsm1L*DLc*m1*m3**4*PsL*tman - Bsm1L*DLc*m1*m3**2*m4**2*PsL*tman - Bsm2L*DLc*m1*m3**2*m4**2*PsL*tman + Bsm2L*DLc*m1*m4**4*PsL*tman + 
     &        2*BsgL*DLc*m1*m3*(m3 - m4)*(m3 + m4)*PsL*tman + Bvm1L*DLc*m1*m3**4*PvR*sman*tman + Bvm1L*DLc*m1*m3**2*m4**2*PvR*sman*tman - 
     &        Bvm2L*DLc*m1*m3**2*m4**2*PvR*sman*tman - Bvm2L*DLc*m1*m4**4*PvR*sman*tman + 2*Bvm1L*DRc*m3**2*m4**2*mq3c*PvR*sman*tman + 
     &        2*Bvm1L*DRc*m4**4*mq3c*PvR*sman*tman + Bsm1L*DLc*m1*PsL*sman**2*tman - Bvm1L*DLc*m1*mq1*PsL*sman**2*tman - 
     &        2*Bvm1L*DLc*m1*m3**2*PvR*sman**2*tman - Bvm1L*DLc*m1*m4**2*PvR*sman**2*tman + Bvm2L*DLc*m1*m4**2*PvR*sman**2*tman - 
     &        2*Bvm1L*DRc*m4**2*mq3c*PvR*sman**2*tman + Bvm1L*DLc*m1*PvR*sman**3*tman + 
     &        2*(Bsm1L*PsL*(DRc*mq3c*(m3**2 - sman)*(m4**2 - sman) + DLc*m1*m4**2*(m3**2 + m4**2 - sman)) - 
     &           Bvm1L*(DRc*mq3c*(m3**2 - sman)*(m4**2 - sman) + DLc*m1*m4**2*(m3**2 + m4**2 - sman))*(mq1*PsL - PvR*sman) + 
     &           DRc*mq3c*(m4**4*(-Bsm2L + Bvm2L*mq1)*PsL - (-(m4**2*(Bsm2L + Bvm2R*m3 - Bvm2L*mq1)*PsL) + 
     &                 (Bsm1R*m3**3 + BsgR*(m3**2 + m4**2) + m4**2*(-(Bsm2R*m3) + Bvm2L*m4**2 + Bvm2R*m3*mq1))*PvR + Bvm1R*m3**3*(PsL - mq1*PvR) + 
     &                 BvgR*(m3**2 + m4**2)*(PsL - mq1*PvR))*sman + 
     &              ((BvgR + Bvm1R*m3)*PsL + (BsgR + Bsm1R*m3 + Bvm2L*m4**2 - (BvgR + Bvm1R*m3)*mq1)*PvR)*sman**2))*uman)*xv + 
     &     2*(2*BvgR*DRc*mq3c*PsL*sman*tman - Bsm1R*DLc*m1*m3*PvR*sman*tman - Bsm2R*DLc*m1*m3*PvR*sman*tman - 2*BvgL*DLc*m1*m3*PvR*sman*tman + 
     &        Bvm1L*DLc*m1*m3**2*PvR*sman*tman - Bvm2L*DLc*m1*m3**2*PvR*sman*tman + 2*BsgR*DRc*mq3c*PvR*sman*tman + 2*Bsm1R*DRc*m3*mq3c*PvR*sman*tman + 
     &        2*Bvm1L*DRc*m3**2*mq3c*PvR*sman*tman - 2*BvgR*DRc*mq1*mq3c*PvR*sman*tman - 2*BsgL*DLc*m1*m3*PsL*uman + 2*Bsm1L*DLc*m1*m4**2*PsL*uman + 
     &        2*BvgL*DLc*m1*m3*mq1*PsL*uman - 2*Bvm1L*DLc*m1*m4**2*mq1*PsL*uman + 2*Bsm1L*DRc*m3**2*mq3c*PsL*uman - Bsm1L*DRc*m4**2*mq3c*PsL*uman - 
     &        Bsm2L*DRc*m4**2*mq3c*PsL*uman - 2*Bvm1L*DRc*m3**2*mq1*mq3c*PsL*uman + Bvm1L*DRc*m4**2*mq1*mq3c*PsL*uman + Bvm2L*DRc*m4**2*mq1*mq3c*PsL*uman - 
     &        2*Bsm1L*DLc*m1*PsL*sman*uman + 2*Bvm1L*DLc*m1*mq1*PsL*sman*uman + Bsm1L*DRc*mq3c*PsL*sman*uman + Bsm2L*DRc*mq3c*PsL*sman*uman - 
     &        Bvm1L*DRc*mq1*mq3c*PsL*sman*uman - Bvm2L*DRc*mq1*mq3c*PsL*sman*uman - 2*Bsm1R*DLc*m1*m3*PvR*sman*uman - 2*BvgL*DLc*m1*m3*PvR*sman*uman + 
     &        2*Bvm1L*DLc*m1*m4**2*PvR*sman*uman + Bsm1R*DRc*m3*mq3c*PvR*sman*uman + Bsm2R*DRc*m3*mq3c*PvR*sman*uman + 2*Bvm1L*DRc*m3**2*mq3c*PvR*sman*uman - 
     &        Bvm1L*DRc*m4**2*mq3c*PvR*sman*uman - Bvm2L*DRc*m4**2*mq3c*PvR*sman*uman - 2*Bvm1L*DLc*m1*PvR*sman**2*uman + Bvm1L*DRc*mq3c*PvR*sman**2*uman + 
     &        Bvm2L*DRc*mq3c*PvR*sman**2*uman - Bvm2R*m3*(PsL - mq1*PvR)*sman*(DLc*m1*tman - DRc*mq3c*uman) - 
     &        Bvm1R*m3*(PsL - mq1*PvR)*sman*(-(DRc*mq3c*(2*tman + uman)) + DLc*m1*(tman + 2*uman)) + BsgL*DLc*m1*m2**2*m3**3*PsL*xv) + 
     &     2*(Bsm2L*DLc*m1*m2**2*m4**2*PsL - Bsm2L*DLc*m1*m3**2*m4**2*PsL + 8*BvgL*DLc*m1*m2**2*m3*mq1*PsL - Bvm1L*DLc*m1*m2**2*m4**2*mq1*PsL - 
     &        Bvm2L*DLc*m1*m2**2*m4**2*mq1*PsL + 8*BvgL*DLc*m1*m3*m4**2*mq1*PsL + Bvm1L*DLc*m1*m3**2*m4**2*mq1*PsL + Bvm2L*DLc*m1*m3**2*m4**2*mq1*PsL + 
     &        Bsm2L*DRc*m1**2*m4**2*mq3c*PsL + Bsm2L*DRc*m3**2*m4**2*mq3c*PsL + 8*BvgL*DRc*m1**2*m3*mq1*mq3c*PsL + 4*BvgL*DRc*m3**3*mq1*mq3c*PsL - 
     &        Bvm1L*DRc*m1**2*m4**2*mq1*mq3c*PsL - Bvm2L*DRc*m1**2*m4**2*mq1*mq3c*PsL + 4*BvgL*DRc*m3*m4**2*mq1*mq3c*PsL - 
     &        Bvm1L*DRc*m3**2*m4**2*mq1*mq3c*PsL - Bvm2L*DRc*m3**2*m4**2*mq1*mq3c*PsL + 4*BvgR*DLc*m1*m2**2*PsL*sman - 4*BvgR*DLc*m1*m3**2*PsL*sman - 
     &        4*BvgL*DLc*m1*m3*mq1*PsL*sman + 4*BvgR*DRc*m1**2*mq3c*PsL*sman + 4*BvgR*DRc*m3**2*mq3c*PsL*sman + 4*BsgR*DLc*m1*m2**2*PvR*sman - 
     &        8*BvgL*DLc*m1*m2**2*m3*PvR*sman - 4*BsgR*DLc*m1*m3**2*PvR*sman + Bvm1L*DLc*m1*m2**2*m4**2*PvR*sman + Bvm2L*DLc*m1*m2**2*m4**2*PvR*sman - 
     &        8*BvgL*DLc*m1*m3*m4**2*PvR*sman - Bvm1L*DLc*m1*m3**2*m4**2*PvR*sman - Bvm2L*DLc*m1*m3**2*m4**2*PvR*sman - 4*BvgR*DLc*m1*m2**2*mq1*PvR*sman + 
     &        4*BvgR*DLc*m1*m3**2*mq1*PvR*sman + 4*BsgR*DRc*m1**2*mq3c*PvR*sman - 8*BvgL*DRc*m1**2*m3*mq3c*PvR*sman + 4*BsgR*DRc*m3**2*mq3c*PvR*sman - 
     &        4*BvgL*DRc*m3**3*mq3c*PvR*sman + Bvm1L*DRc*m1**2*m4**2*mq3c*PvR*sman + Bvm2L*DRc*m1**2*m4**2*mq3c*PvR*sman - 4*BvgL*DRc*m3*m4**2*mq3c*PvR*sman + 
     &        Bvm1L*DRc*m3**2*m4**2*mq3c*PvR*sman + Bvm2L*DRc*m3**2*m4**2*mq3c*PvR*sman - 4*BvgR*DRc*m1**2*mq1*mq3c*PvR*sman - 
     &        4*BvgR*DRc*m3**2*mq1*mq3c*PvR*sman + 4*BvgL*DLc*m1*m3*PvR*sman**2 - 4*BvgL*DLc*m1*m3*mq1*PsL*tman - Bsm2L*DRc*m4**2*mq3c*PsL*tman - 
     &        4*BvgL*DRc*m3*mq1*mq3c*PsL*tman + Bvm1L*DRc*m4**2*mq1*mq3c*PsL*tman + Bvm2L*DRc*m4**2*mq1*mq3c*PsL*tman - 4*BvgR*DRc*mq3c*PsL*sman*tman + 
     &        4*BvgL*DLc*m1*m3*PvR*sman*tman - 4*BsgR*DRc*mq3c*PvR*sman*tman + 4*BvgL*DRc*m3*mq3c*PvR*sman*tman - Bvm1L*DRc*m4**2*mq3c*PvR*sman*tman - 
     &        Bvm2L*DRc*m4**2*mq3c*PvR*sman*tman + 4*BvgR*DRc*mq1*mq3c*PvR*sman*tman - Bsm2L*DLc*m1*m4**2*PsL*uman - 4*BvgL*DLc*m1*m3*mq1*PsL*uman + 
     &        Bvm1L*DLc*m1*m4**2*mq1*PsL*uman + Bvm2L*DLc*m1*m4**2*mq1*PsL*uman - 4*BvgL*DRc*m3*mq1*mq3c*PsL*uman - 4*BvgR*DLc*m1*PsL*sman*uman - 
     &        4*BsgR*DLc*m1*PvR*sman*uman + 4*BvgL*DLc*m1*m3*PvR*sman*uman - Bvm1L*DLc*m1*m4**2*PvR*sman*uman - Bvm2L*DLc*m1*m4**2*PvR*sman*uman + 
     &        4*BvgR*DLc*m1*mq1*PvR*sman*uman + 4*BvgL*DRc*m3*mq3c*PvR*sman*uman + 
     &        m4**2*(((Bvm1L + Bvm2L)*m4**2*mq1*PsL - (BvgR*PsL + (BsgR + (Bvm1L + Bvm2L)*m4**2 - BvgR*mq1)*PvR)*sman)*
     &            (DRc*mq3c*(m1**2 + m3**2 - tman) + DLc*m1*(m2**2 - m3**2 - uman)) + 
     &           Bsm2L*m4**2*PsL*(-(DRc*mq3c*(m1**2 + m3**2 - tman)) + DLc*m1*(-m2**2 + m3**2 + uman)) + 
     &           BvgL*m3*(mq1*PsL - PvR*sman)*(DRc*mq3c*(-2*m1**2 - m3**2 - m4**2 + tman + uman) + DLc*m1*(-2*m2**2 - 2*m4**2 + sman + tman + uman)))*xv + 
     &        BsgL*m3*PsL*(DRc*mq3c*(2*m1**2 + m3**2 + m4**2 - tman - uman) + DLc*m1*(2*m2**2 + 2*m4**2 - sman - tman - uman))*(-4 + m4**2*xv) - 
     &        Bsm1L*m4**2*PsL*(DRc*mq3c*(m1**2 + m3**2 - tman) + DLc*m1*(m2**2 - m3**2 - uman))*(-1 + m4**2*xv))) + 
     &  A1vL*A3Rc*(-2*(-(Bsm2R*DRc*m1*m2**2*m3**2*PsR) - Bsm2R*DRc*m1*m3**4*PsR - 2*Bsm2R*DRc*m1*m3**2*m4**2*PsR + 4*BvgR*DRc*m1*m2**2*m3*mq1*PsR - 
     &        Bvm1R*DRc*m1*m2**2*m3**2*mq1*PsR + Bvm2R*DRc*m1*m2**2*m3**2*mq1*PsR - Bvm1R*DRc*m1*m3**4*mq1*PsR + Bvm2R*DRc*m1*m3**4*mq1*PsR - 
     &        2*Bvm1R*DRc*m1*m2**2*m4**2*mq1*PsR + 4*BvgR*DRc*m1*m3*m4**2*mq1*PsR + 2*Bvm2R*DRc*m1*m3**2*m4**2*mq1*PsR - Bsm2R*DLc*m1**2*m4**2*mq3c*PsR - 
     &        Bsm2R*DLc*m4**4*mq3c*PsR - 4*Bvm1R*DLc*m1**2*m3**2*mq1*mq3c*PsR - 2*Bvm1R*DLc*m3**4*mq1*mq3c*PsR + Bvm1R*DLc*m1**2*m4**2*mq1*mq3c*PsR + 
     &        Bvm2R*DLc*m1**2*m4**2*mq1*mq3c*PsR - 2*Bvm1R*DLc*m3**2*m4**2*mq1*mq3c*PsR + Bvm1R*DLc*m4**4*mq1*mq3c*PsR + Bvm2R*DLc*m4**4*mq1*mq3c*PsR - 
     &        ((DRc*m1*(-2*Bvm1R*m2**2*mq1 + m3*(-(Bsm2R*m3) + 2*BvgR*mq1 + (Bvm1R + Bvm2R)*m3*mq1)) + 
     &              DLc*(-2*BvgL*(m1**2 + m3**2) - Bsm2R*(m1**2 + m4**2) + (Bvm1R + Bvm2R)*(m1**2 + m4**2)*mq1)*mq3c)*PsR + 
     &           (DRc*m1*(m3*((3*Bsm1L + Bsm2L + 4*BvgR)*m2**2 + (-Bvm1R + Bvm2R)*m2**2*m3 + (-Bsm1L + Bsm2L)*m3**2 + (-Bvm1R + Bvm2R)*m3**3) + 
     &                 2*(-(Bvm1R*m2**2) + m3*(Bsm1L + Bsm2L + 2*BvgR + Bvm2R*m3))*m4**2) - 
     &              DLc*(Bsm2L*m1**2*m3 + 4*Bvm1R*m1**2*m3**2 + 2*Bvm1R*m3**4 + 2*BsgL*(m1**2 + m3**2) - Bvm1R*m1**2*m4**2 - Bvm2R*m1**2*m4**2 + 
     &                 Bsm2L*m3*m4**2 + 2*Bvm1R*m3**2*m4**2 - Bvm1R*m4**4 - Bvm2R*m4**4 + Bsm1L*m3*(3*m1**2 + 2*m3**2 + m4**2) - 2*BvgL*(m1**2 + m3**2)*mq1)*
     &               mq3c)*PvL + Bvm2L*m3*(DRc*m1*(m2**2 + m3**2 + 2*m4**2) - DLc*(m1**2 + m4**2)*mq3c)*(PsR - mq1*PvL) - 
     &           Bvm1L*m3*(DRc*m1*(-3*m2**2 + m3**2 - 2*m4**2) + DLc*(3*m1**2 + 2*m3**2 + m4**2)*mq3c)*(PsR - mq1*PvL))*sman + 
     &        ((Bvm1L + Bvm2L)*DRc*m1*m3*PsR + (DRc*m1*(Bvm1R*(-2*m2**2 + m3**2) + m3*(Bsm1L + Bsm2L + 2*BvgR + Bvm2R*m3 - (Bvm1L + Bvm2L)*mq1)) + 
     &              (Bvm1R + Bvm2R)*DLc*(m1**2 + m4**2)*mq3c)*PvL)*sman**2 + 
     &        Bsm1R*PsR*(DRc*m1*(m2**2*(m3**2 + 2*m4**2 - 2*sman) + m3**2*(m3**2 + sman)) + 
     &           DLc*mq3c*(2*m3**4 + 2*m3**2*m4**2 - m4**4 + m4**2*sman + m1**2*(4*m3**2 - m4**2 + sman))) + 
     &        2*BsgR*DRc*m1*m3*PsR*(-2*m2**2 - 2*m4**2 + sman + tman)) + 
     &     2*(Bsm1R*DRc*m1*m3**2*PsR*tman - Bsm2R*DRc*m1*m3**2*PsR*tman + 2*BvgR*DRc*m1*m3*mq1*PsR*tman - Bvm1R*DRc*m1*m3**2*mq1*PsR*tman + 
     &        Bvm2R*DRc*m1*m3**2*mq1*PsR*tman + 2*Bsm1R*DLc*m3**2*mq3c*PsR*tman - 2*Bvm1R*DLc*m3**2*mq1*mq3c*PsR*tman - Bvm1L*DRc*m1*m3*PsR*sman*tman - 
     &        Bvm2L*DRc*m1*m3*PsR*sman*tman + 2*BvgL*DLc*mq3c*PsR*sman*tman + 2*Bvm1L*DLc*m3*mq3c*PsR*sman*tman - Bsm1L*DRc*m1*m3*PvL*sman*tman - 
     &        Bsm2L*DRc*m1*m3*PvL*sman*tman - 2*BvgR*DRc*m1*m3*PvL*sman*tman + Bvm1R*DRc*m1*m3**2*PvL*sman*tman - Bvm2R*DRc*m1*m3**2*PvL*sman*tman + 
     &        Bvm1L*DRc*m1*m3*mq1*PvL*sman*tman + Bvm2L*DRc*m1*m3*mq1*PvL*sman*tman + 2*BsgL*DLc*mq3c*PvL*sman*tman + 2*Bsm1L*DLc*m3*mq3c*PvL*sman*tman + 
     &        2*Bvm1R*DLc*m3**2*mq3c*PvL*sman*tman - 2*BvgL*DLc*mq1*mq3c*PvL*sman*tman - 2*Bvm1L*DLc*m3*mq1*mq3c*PvL*sman*tman + 
     &        ((2*DRc*m1*(-(BsgR*m3) + BvgR*m3*mq1 + m4**2*(Bsm1R - Bvm1R*mq1)) + 
     &              DLc*(Bsm1R*(2*m3**2 - m4**2) - 2*Bvm1R*m3**2*mq1 + m4**2*(-Bsm2R + (Bvm1R + Bvm2R)*mq1))*mq3c)*PsR + 
     &           ((-2*DRc*m1*(Bsm1R + Bvm1L*m3 - Bvm1R*mq1) + DLc*(Bsm1R + Bsm2R + (Bvm1L + Bvm2L)*m3 - (Bvm1R + Bvm2R)*mq1)*mq3c)*PsR - 
     &              (2*DRc*m1*(-(Bvm1R*m4**2) + m3*(Bsm1L + BvgR - Bvm1L*mq1)) + 
     &                 DLc*(-(Bsm1L*m3) - Bsm2L*m3 - 2*Bvm1R*m3**2 + Bvm1R*m4**2 + Bvm2R*m4**2 + (Bvm1L + Bvm2L)*m3*mq1)*mq3c)*PvL)*sman + 
     &           (-2*Bvm1R*DRc*m1 + (Bvm1R + Bvm2R)*DLc*mq3c)*PvL*sman**2)*uman) - 
     &     2*(Bsm2R*DRc*m1*m2**2*m3**2*m4**2*PsR + Bsm2R*DRc*m1*m3**4*m4**2*PsR + 2*Bsm2R*DRc*m1*m3**2*m4**4*PsR + BvgR*DRc*m1*m2**2*m3**3*mq1*PsR + 
     &        Bvm1R*DRc*m1*m2**2*m3**4*mq1*PsR + BvgR*DRc*m1*m3**5*mq1*PsR + Bvm1R*DRc*m1*m3**6*mq1*PsR - BvgR*DRc*m1*m2**2*m3*m4**2*mq1*PsR + 
     &        Bvm1R*DRc*m1*m2**2*m3**2*m4**2*mq1*PsR - Bvm2R*DRc*m1*m2**2*m3**2*m4**2*mq1*PsR + BvgR*DRc*m1*m3**3*m4**2*mq1*PsR + 
     &        Bvm1R*DRc*m1*m3**4*m4**2*mq1*PsR - Bvm2R*DRc*m1*m3**4*m4**2*mq1*PsR + Bvm1R*DRc*m1*m2**2*m4**4*mq1*PsR - 2*BvgR*DRc*m1*m3*m4**4*mq1*PsR - 
     &        Bvm1R*DRc*m1*m3**2*m4**4*mq1*PsR - 2*Bvm2R*DRc*m1*m3**2*m4**4*mq1*PsR + Bsm2R*DLc*m1**2*m4**4*mq3c*PsR + Bsm2R*DLc*m4**6*mq3c*PsR + 
     &        2*Bvm1R*DLc*m1**2*m3**2*m4**2*mq1*mq3c*PsR + Bvm1R*DLc*m3**4*m4**2*mq1*mq3c*PsR + Bvm1R*DLc*m1**2*m4**4*mq1*mq3c*PsR - 
     &        Bvm2R*DLc*m1**2*m4**4*mq1*mq3c*PsR + 2*Bvm1R*DLc*m3**2*m4**4*mq1*mq3c*PsR - Bvm2R*DLc*m4**6*mq1*mq3c*PsR - 
     &        (Bsm2R*DRc*m1*m3**2*m4**2*PsR + BvgR*DRc*m1*m2**2*m3*mq1*PsR + Bvm1R*DRc*m1*m2**2*m3**2*mq1*PsR + 2*BvgR*DRc*m1*m3**3*mq1*PsR + 
     &           2*Bvm1R*DRc*m1*m3**4*mq1*PsR + Bvm1R*DRc*m1*m2**2*m4**2*mq1*PsR + BvgR*DRc*m1*m3*m4**2*mq1*PsR + Bvm1R*DRc*m1*m3**2*m4**2*mq1*PsR - 
     &           Bvm2R*DRc*m1*m3**2*m4**2*mq1*PsR - BvgL*DLc*m1**2*m3**2*mq3c*PsR + Bsm2R*DLc*m1**2*m4**2*mq3c*PsR - BvgL*DLc*m1**2*m4**2*mq3c*PsR - 
     &           BvgL*DLc*m3**2*m4**2*mq3c*PsR + Bsm2R*DLc*m4**4*mq3c*PsR - BvgL*DLc*m4**4*mq3c*PsR + Bvm1R*DLc*m1**2*m3**2*mq1*mq3c*PsR + 
     &           2*Bvm1R*DLc*m1**2*m4**2*mq1*mq3c*PsR - Bvm2R*DLc*m1**2*m4**2*mq1*mq3c*PsR + 2*Bvm1R*DLc*m3**2*m4**2*mq1*mq3c*PsR + 
     &           Bvm1R*DLc*m4**4*mq1*mq3c*PsR - Bvm2R*DLc*m4**4*mq1*mq3c*PsR + 
     &           (BvgR*DRc*m1*m3*(m3 - m4)*(m3 + m4)*(m2**2 + m3**2 + 2*m4**2) + 
     &              DRc*m1*(Bvm1R*m3**4*(m2**2 + m3**2) - m3*(Bsm2L + (-Bvm1R + Bvm2R)*m3)*(m2**2 + m3**2)*m4**2 + 
     &                 (Bvm1R*(m2 - m3)*(m2 + m3) - 2*m3*(Bsm2L + Bvm2R*m3))*m4**4) + 
     &              DLc*(m3*(Bsm2L*m1**2 + Bvm1R*m3*(2*m1**2 + m3**2))*m4**2 + ((Bvm1R - Bvm2R)*m1**2 + Bsm2L*m3 + 2*Bvm1R*m3**2)*m4**4 - Bvm2R*m4**6 - 
     &                 BsgL*(m1**2 + m4**2)*(m3**2 + m4**2) + BvgL*(m1**2 + m4**2)*(m3**2 + m4**2)*mq1)*mq3c + 
     &              Bsm1L*m3**3*(DRc*m1*(m2**2 + m3**2 + 2*m4**2) - DLc*(m1**2 + m4**2)*mq3c))*PvL + 
     &           Bvm1L*m3**3*(DRc*m1*(m2**2 + m3**2 + 2*m4**2) - DLc*(m1**2 + m4**2)*mq3c)*(PsR - mq1*PvL) - 
     &           Bvm2L*m3*m4**2*(DRc*m1*(m2**2 + m3**2 + 2*m4**2) - DLc*(m1**2 + m4**2)*mq3c)*(PsR - mq1*PvL))*sman + 
     &        BsgR*DRc*m1*m3*PsR*(-((m3 - m4)*(m3 + m4)*(m2**2 + m3**2 + 2*m4**2)) + (m2**2 + 2*m3**2 + m4**2)*sman) + 
     &        Bsm1R*PsR*(DLc*mq3c*(m4**2*(-m3**4 - 2*m3**2*(m4**2 - sman) + m4**2*sman) + m1**2*(-m4**4 + 2*m4**2*sman + m3**2*(-2*m4**2 + sman))) - 
     &           DRc*m1*(m2**2*(m3**4 + m4**4 + m3**2*(m4**2 - sman) - m4**2*sman) + m3**2*(m3**4 + m3**2*(m4**2 - 2*sman) - m4**2*(m4**2 + sman)))))*xv + 
     &     DRc*m1*((m3 - m4)*(m3 + m4)*(-2*BsgR*m3 - Bsm1R*m3**2 + m3*(2*BvgR + Bvm1R*m3)*mq1 + m4**2*(Bsm2R - Bvm2R*mq1))*PsR + 
     &        ((2*BsgR*m3 - 2*Bvm1L*m3**3 + (-Bsm1R + Bsm2R + 2*Bvm2L*m3)*m4**2 + (-2*BvgR*m3 + (Bvm1R - Bvm2R)*m4**2)*mq1)*PsR + 
     &           (-2*Bsm1L*m3**3 - 2*BvgR*m3**3 - Bvm1R*m3**4 + 2*Bsm2L*m3*m4**2 + 2*BvgR*m3*m4**2 + Bvm1R*m3**2*m4**2 + Bvm2R*m3**2*m4**2 - Bvm2R*m4**4 + 
     &              2*Bvm1L*m3**3*mq1 - 2*Bvm2L*m3*m4**2*mq1)*PvL)*sman + 
     &        (Bsm1R*PsR - Bvm1R*mq1*PsR + 2*(Bsm1L + BvgR)*m3*PvL + (-Bvm1R + Bvm2R)*m4**2*PvL + 2*Bvm1L*m3*(PsR - mq1*PvL))*sman**2 + Bvm1R*PvL*sman**3)*
     &      tman*xv + (2*BsgR*DRc*m1*m3*PsR*sman**2 + 2*Bsm1R*DRc*m1*m3**2*PsR*sman**2 + 2*Bvm2L*DRc*m1*m3*m4**2*PsR*sman**2 - 
     &        2*BvgR*DRc*m1*m3*mq1*PsR*sman**2 - 2*Bvm1R*DRc*m1*m3**2*mq1*PsR*sman**2 + 2*Bsm1R*DLc*m1**2*mq3c*PsR*sman**2 + 
     &        2*BvgL*DLc*m1**2*mq3c*PsR*sman**2 + 2*Bsm1R*DLc*m4**2*mq3c*PsR*sman**2 + 2*BvgL*DLc*m4**2*mq3c*PsR*sman**2 - 
     &        2*Bvm1R*DLc*m1**2*mq1*mq3c*PsR*sman**2 - 2*Bvm1R*DLc*m4**2*mq1*mq3c*PsR*sman**2 - 2*Bsm1L*DRc*m1*m2**2*m3*PvL*sman**2 - 
     &        2*BvgR*DRc*m1*m2**2*m3*PvL*sman**2 - 2*Bvm1R*DRc*m1*m2**2*m3**2*PvL*sman**2 - 4*Bsm1L*DRc*m1*m3**3*PvL*sman**2 - 
     &        4*BvgR*DRc*m1*m3**3*PvL*sman**2 - 4*Bvm1R*DRc*m1*m3**4*PvL*sman**2 - 2*Bvm1R*DRc*m1*m2**2*m4**2*PvL*sman**2 - 
     &        4*Bsm1L*DRc*m1*m3*m4**2*PvL*sman**2 + 2*Bsm2L*DRc*m1*m3*m4**2*PvL*sman**2 - 2*BvgR*DRc*m1*m3*m4**2*PvL*sman**2 - 
     &        2*Bvm1R*DRc*m1*m3**2*m4**2*PvL*sman**2 + 2*Bvm2R*DRc*m1*m3**2*m4**2*PvL*sman**2 - 2*Bvm2L*DRc*m1*m3*m4**2*mq1*PvL*sman**2 + 
     &        2*BsgL*DLc*m1**2*mq3c*PvL*sman**2 + 2*Bsm1L*DLc*m1**2*m3*mq3c*PvL*sman**2 - 2*Bvm1R*DLc*m1**2*m3**2*mq3c*PvL*sman**2 + 
     &        2*BsgL*DLc*m4**2*mq3c*PvL*sman**2 - 4*Bvm1R*DLc*m1**2*m4**2*mq3c*PvL*sman**2 + 2*Bvm2R*DLc*m1**2*m4**2*mq3c*PvL*sman**2 + 
     &        2*Bsm1L*DLc*m3*m4**2*mq3c*PvL*sman**2 - 4*Bvm1R*DLc*m3**2*m4**2*mq3c*PvL*sman**2 - 2*Bvm1R*DLc*m4**4*mq3c*PvL*sman**2 + 
     &        2*Bvm2R*DLc*m4**4*mq3c*PvL*sman**2 - 2*BvgL*DLc*m1**2*mq1*mq3c*PvL*sman**2 - 2*BvgL*DLc*m4**2*mq1*mq3c*PvL*sman**2 - 
     &        2*Bvm1L*m3*(PsR - mq1*PvL)*(-(DLc*(m1**2 + m4**2)*mq3c) + DRc*m1*(m2**2 + 2*(m3**2 + m4**2) - sman))*sman**2 + 2*Bsm1L*DRc*m1*m3*PvL*sman**3 + 
     &        2*BvgR*DRc*m1*m3*PvL*sman**3 + 2*Bvm1R*DRc*m1*m3**2*PvL*sman**3 + 2*Bvm1R*DLc*m1**2*mq3c*PvL*sman**3 + 2*Bvm1R*DLc*m4**2*mq3c*PvL*sman**3 - 
     &        ((m3**2 + m4**2)*(DRc*m1*(Bsm1R*m3**2 - Bvm1R*m3**2*mq1 + m4**2*(-Bsm2R + Bvm2R*mq1)) + 2*DLc*m4**2*(Bsm1R - Bvm1R*mq1)*mq3c)*PsR + 
     &           ((DRc*m1*(-(Bsm1R*(2*m3**2 + m4**2)) + 2*Bvm1R*m3**2*mq1 + m4**2*(Bsm2R + (Bvm1R - Bvm2R)*mq1)) - 2*DLc*m4**2*(Bsm1R - Bvm1R*mq1)*mq3c)*PsR + 
     &              Bvm1R*DRc*m1*m3**2*(m3**2 + m4**2)*PvL)*sman)*tman)*xv + 
     &     (sman*(Bvm2R*DRc*m1*m4**2*PvL*(m3**2 + m4**2 - sman) - Bsm1R*DRc*m1*PsR*sman + 
     &           Bvm1R*(-2*DLc*m4**2*mq3c*PvL*(m3**2 + m4**2 - sman) + DRc*m1*(mq1*PsR + PvL*(2*m3**2 + m4**2 - sman))*sman))*tman - 
     &        2*(Bsm1R*PsR*(DLc*mq3c*(m3**2 - sman)*(m4**2 - sman) + DRc*m1*m4**2*(m3**2 + m4**2 - sman)) - 
     &           Bvm1R*(DLc*mq3c*(m3**2 - sman)*(m4**2 - sman) + DRc*m1*m4**2*(m3**2 + m4**2 - sman))*(mq1*PsR - PvL*sman) + 
     &           DLc*mq3c*(m4**4*(-Bsm2R + Bvm2R*mq1)*PsR - (-(m4**2*(Bsm2R + Bvm2L*m3 - Bvm2R*mq1)*PsR) + 
     &                 (Bsm1L*m3**3 + BsgL*(m3**2 + m4**2) + m4**2*(-(Bsm2L*m3) + Bvm2R*m4**2 + Bvm2L*m3*mq1))*PvL + Bvm1L*m3**3*(PsR - mq1*PvL) + 
     &                 BvgL*(m3**2 + m4**2)*(PsR - mq1*PvL))*sman + 
     &              ((BvgL + Bvm1L*m3)*PsR + (BsgL + Bsm1L*m3 + Bvm2R*m4**2 - (BvgL + Bvm1L*m3)*mq1)*PvL)*sman**2))*uman)*xv + 
     &     2*(Bsm2R*DRc*m1*m2**2*m4**2*PsR - Bsm2R*DRc*m1*m3**2*m4**2*PsR + 8*BvgR*DRc*m1*m2**2*m3*mq1*PsR - Bvm1R*DRc*m1*m2**2*m4**2*mq1*PsR - 
     &        Bvm2R*DRc*m1*m2**2*m4**2*mq1*PsR + 8*BvgR*DRc*m1*m3*m4**2*mq1*PsR + Bvm1R*DRc*m1*m3**2*m4**2*mq1*PsR + Bvm2R*DRc*m1*m3**2*m4**2*mq1*PsR + 
     &        Bsm2R*DLc*m1**2*m4**2*mq3c*PsR + Bsm2R*DLc*m3**2*m4**2*mq3c*PsR + 8*BvgR*DLc*m1**2*m3*mq1*mq3c*PsR + 4*BvgR*DLc*m3**3*mq1*mq3c*PsR - 
     &        Bvm1R*DLc*m1**2*m4**2*mq1*mq3c*PsR - Bvm2R*DLc*m1**2*m4**2*mq1*mq3c*PsR + 4*BvgR*DLc*m3*m4**2*mq1*mq3c*PsR - 
     &        Bvm1R*DLc*m3**2*m4**2*mq1*mq3c*PsR - Bvm2R*DLc*m3**2*m4**2*mq1*mq3c*PsR + 4*BvgL*DRc*m1*m2**2*PsR*sman - 4*BvgL*DRc*m1*m3**2*PsR*sman - 
     &        4*BvgR*DRc*m1*m3*mq1*PsR*sman + 4*BvgL*DLc*m1**2*mq3c*PsR*sman + 4*BvgL*DLc*m3**2*mq3c*PsR*sman + 4*BsgL*DRc*m1*m2**2*PvL*sman - 
     &        8*BvgR*DRc*m1*m2**2*m3*PvL*sman - 4*BsgL*DRc*m1*m3**2*PvL*sman + Bvm1R*DRc*m1*m2**2*m4**2*PvL*sman + Bvm2R*DRc*m1*m2**2*m4**2*PvL*sman - 
     &        8*BvgR*DRc*m1*m3*m4**2*PvL*sman - Bvm1R*DRc*m1*m3**2*m4**2*PvL*sman - Bvm2R*DRc*m1*m3**2*m4**2*PvL*sman - 4*BvgL*DRc*m1*m2**2*mq1*PvL*sman + 
     &        4*BvgL*DRc*m1*m3**2*mq1*PvL*sman + 4*BsgL*DLc*m1**2*mq3c*PvL*sman - 8*BvgR*DLc*m1**2*m3*mq3c*PvL*sman + 4*BsgL*DLc*m3**2*mq3c*PvL*sman - 
     &        4*BvgR*DLc*m3**3*mq3c*PvL*sman + Bvm1R*DLc*m1**2*m4**2*mq3c*PvL*sman + Bvm2R*DLc*m1**2*m4**2*mq3c*PvL*sman - 4*BvgR*DLc*m3*m4**2*mq3c*PvL*sman + 
     &        Bvm1R*DLc*m3**2*m4**2*mq3c*PvL*sman + Bvm2R*DLc*m3**2*m4**2*mq3c*PvL*sman - 4*BvgL*DLc*m1**2*mq1*mq3c*PvL*sman - 
     &        4*BvgL*DLc*m3**2*mq1*mq3c*PvL*sman + 4*BvgR*DRc*m1*m3*PvL*sman**2 - 4*BvgR*DRc*m1*m3*mq1*PsR*tman - Bsm2R*DLc*m4**2*mq3c*PsR*tman - 
     &        4*BvgR*DLc*m3*mq1*mq3c*PsR*tman + Bvm1R*DLc*m4**2*mq1*mq3c*PsR*tman + Bvm2R*DLc*m4**2*mq1*mq3c*PsR*tman - 4*BvgL*DLc*mq3c*PsR*sman*tman + 
     &        4*BvgR*DRc*m1*m3*PvL*sman*tman - 4*BsgL*DLc*mq3c*PvL*sman*tman + 4*BvgR*DLc*m3*mq3c*PvL*sman*tman - Bvm1R*DLc*m4**2*mq3c*PvL*sman*tman - 
     &        Bvm2R*DLc*m4**2*mq3c*PvL*sman*tman + 4*BvgL*DLc*mq1*mq3c*PvL*sman*tman - Bsm2R*DRc*m1*m4**2*PsR*uman - 4*BvgR*DRc*m1*m3*mq1*PsR*uman + 
     &        Bvm1R*DRc*m1*m4**2*mq1*PsR*uman + Bvm2R*DRc*m1*m4**2*mq1*PsR*uman - 4*BvgR*DLc*m3*mq1*mq3c*PsR*uman - 4*BvgL*DRc*m1*PsR*sman*uman - 
     &        4*BsgL*DRc*m1*PvL*sman*uman + 4*BvgR*DRc*m1*m3*PvL*sman*uman - Bvm1R*DRc*m1*m4**2*PvL*sman*uman - Bvm2R*DRc*m1*m4**2*PvL*sman*uman + 
     &        4*BvgL*DRc*m1*mq1*PvL*sman*uman + 4*BvgR*DLc*m3*mq3c*PvL*sman*uman + 
     &        m4**2*(((Bvm1R + Bvm2R)*m4**2*mq1*PsR - (BvgL*PsR + (BsgL + (Bvm1R + Bvm2R)*m4**2 - BvgL*mq1)*PvL)*sman)*
     &            (DLc*mq3c*(m1**2 + m3**2 - tman) + DRc*m1*(m2**2 - m3**2 - uman)) + 
     &           Bsm2R*m4**2*PsR*(-(DLc*mq3c*(m1**2 + m3**2 - tman)) + DRc*m1*(-m2**2 + m3**2 + uman)) + 
     &           BvgR*m3*(mq1*PsR - PvL*sman)*(DLc*mq3c*(-2*m1**2 - m3**2 - m4**2 + tman + uman) + DRc*m1*(-2*m2**2 - 2*m4**2 + sman + tman + uman)))*xv + 
     &        BsgR*m3*PsR*(DLc*mq3c*(2*m1**2 + m3**2 + m4**2 - tman - uman) + DRc*m1*(2*m2**2 + 2*m4**2 - sman - tman - uman))*(-4 + m4**2*xv) - 
     &        Bsm1R*m4**2*PsR*(DLc*mq3c*(m1**2 + m3**2 - tman) + DRc*m1*(m2**2 - m3**2 - uman))*(-1 + m4**2*xv))) ))/
     &  ((sman - mq1**2 + cI*mq1*wquark)*(uman - mq3c**2)))
        
        end

        double precision function MtsV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A1Lc,A1Rc,BLc,BRc,mq2,mq1c,xv,wquark)

        implicit none
        
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

        double complex A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A1Lc,A1Rc,BLc,BRc
        double precision mq2,mq1c,xv,wquark

	MtsV = 0d0

        MtsV = dreal((A1Rc*(A2sL*m1*(BLc*m3*(-1.*CLoop2*(m3**2 - 1.*m4**2 - 1.*sman)*(-1. + m4**2*xv) + 
     &             CLoop1*(2.*m1**2 - 2.*m2**2 - 1.*m3**2 + m4**2 - 1.*sman - 1.*(m3**2 - 1.*m4**2 - 1.*sman)*(m2**2 - 1.*tman)*xv)) + 
     &          BRc*mq1c*(CLoop2*(m3**2 + m4**2 - 1.*sman)*(1. - 1.*m4**2*xv) + 
     &             CLoop1*(2.*m1**2 - 1.*m3**2 + m4**2 - 1.*sman - 2.*tman + (m3**2 + m4**2 - 1.*sman)*(-1.*m2**2 + tman)*xv))) + 
     &       A2sR*(BRc*m3*mq1c*(CLoop2*(m2**2 + m3**2 - 1.*sman - 1.*tman)*(1. - 1.*m4**2*xv) + 
     &             CLoop1*(2.*m1**2 + m2**2 - 1.*m3**2 - 1.*sman + tman - 1.*(m2**2 - 1.*tman)*(m2**2 + m3**2 - 1.*sman - 1.*tman)*xv)) + 
     &          BLc*(CLoop2*(m2**2*m3**2 + m1**2*(-1.*m3**2 + sman) - 1.*sman*tman)*(-1. + m4**2*xv) + 
     &             CLoop1*(2.*m1**4 - 1.*m1**2*(2.*m2**2 + m3**2 - 2.*m4**2 + sman + 2.*tman) + m2**4*m3**2*xv + m1**2*(m3**2 - 1.*sman)*(-1.*m2**2 + tman)*xv + 
     &                sman*tman*(-1. + tman*xv) - 1.*m2**2*(m3**2 - 2.*tman + (m3**2 + sman)*tman*xv))))) + 
     &    A1Lc*(A2sR*m1*(BRc*m3*(-1.*CLoop2*(m3**2 - 1.*m4**2 - 1.*sman)*(-1. + m4**2*xv) + 
     &             CLoop1*(2.*m1**2 - 2.*m2**2 - 1.*m3**2 + m4**2 - 1.*sman - 1.*(m3**2 - 1.*m4**2 - 1.*sman)*(m2**2 - 1.*tman)*xv)) + 
     &          BLc*mq1c*(CLoop2*(m3**2 + m4**2 - 1.*sman)*(1. - 1.*m4**2*xv) + 
     &             CLoop1*(2.*m1**2 - 1.*m3**2 + m4**2 - 1.*sman - 2.*tman + (m3**2 + m4**2 - 1.*sman)*(-1.*m2**2 + tman)*xv))) + 
     &       A2sL*(BLc*m3*mq1c*(CLoop2*(m2**2 + m3**2 - 1.*sman - 1.*tman)*(1. - 1.*m4**2*xv) + 
     &             CLoop1*(2.*m1**2 + m2**2 - 1.*m3**2 - 1.*sman + tman - 1.*(m2**2 - 1.*tman)*(m2**2 + m3**2 - 1.*sman - 1.*tman)*xv)) + 
     &          BRc*(CLoop2*(m2**2*m3**2 + m1**2*(-1.*m3**2 + sman) - 1.*sman*tman)*(-1. + m4**2*xv) + 
     &             CLoop1*(2.*m1**4 - 1.*m1**2*(2.*m2**2 + m3**2 - 2.*m4**2 + sman + 2.*tman) + m2**4*m3**2*xv + m1**2*(m3**2 - 1.*sman)*(-1.*m2**2 + tman)*xv + 
     &                sman*tman*(-1. + tman*xv) - 1.*m2**2*(m3**2 - 2.*tman + (m3**2 + sman)*tman*xv))))))/((mq1c**2 - 1.*sman + cI*mq1c*wquark)*(mq2**2 - 1.*tman)))
     

        end

! --- Mtt  ---
 
        double precision function MttV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A2Lc,A2Rc,CLc,mq2,mq2c,xv)

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

        double complex A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A2Lc,A2Rc,CLc
        double precision mq2,mq2c,xv

        MttV = 0d0
        MttV = dreal((CLc*(A2Lc*(2.*A2sR*m1*m3 + A2sL*(m1**2 + m3**2 - 1.*tman)) + A2Rc*(2.*A2sL*m1*m3 + A2sR*(m1**2 + m3**2 - 1.*tman)))*
     &    (CLoop2*(m2**2 - 1.*tman)*(-1. + m4**2*xv) + CLoop1*(-2.*m2**2 + m4**2 - 2.*tman + (m2**2 - 1.*tman)**2*xv)))/((mq2**2 - 1.*tman)*(mq2c**2 - 1.*tman)))


        end



! --- Mtu  ---

        double precision function MtuV(A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A3Lc,A3Rc,DLc,DRc,mq2,mq3c,xv)

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
        
        double complex A2sL,A2sR,A2vL,A2vR,PsL,PsR,PvL,PvR,CLoop1,CLoop2,A3Lc,A3Rc,DLc,DRc
        double precision mq2,mq3c,xv

        MtuV =  0d0
        MtuV = dreal((A2sL*(A3Rc*m3*(CLoop2*(DRc*m1*(m1**2 - 1.*m4**2 - 1.*uman) + DLc*mq3c*(m1**2 + m4**2 - 1.*uman))*(-1. + m4**2*xv) + 
     &          CLoop1*(DRc*m1*(m1**2 + 2.*m2**2 - 2.*m3**2 - 1.*m4**2 + uman + (m2**2 - 1.*tman)*(m1**2 - 1.*m4**2 - 1.*uman)*xv) + 
     &             DLc*mq3c*(m1**2 - 2.*m3**2 - 1.*m4**2 + 2.*tman + uman + (m2**2 - 1.*tman)*(m1**2 + m4**2 - 1.*uman)*xv))) + 
     &       A3Lc*(-1.*CLoop2*(DLc*m1**2*(m2 - 1.*m3)*(m2 + m3) + DLc*(m3**2 - 1.*tman)*uman + DRc*m1*mq3c*(-1.*m1**2 - 1.*m2**2 + tman + uman))*(-1. + m4**2*xv) + 
     &          CLoop1*(DRc*m1*mq3c*(m1**2 - 1.*m2**2 - 2.*m3**2 - 1.*tman + uman + (m2**2 - 1.*tman)*(m1**2 + m2**2 - 1.*tman - 1.*uman)*xv) + 
     &             DLc*(m1**2*m2**2 + m1**2*m3**2 + 2.*m2**2*m3**2 - 2.*m3**4 - 2.*m3**2*m4**2 - 2.*m2**2*tman + 2.*m3**2*tman + m3**2*uman + tman*uman + 
     &                (-1.*m2**2 + tman)*(m1**2*(m2 - 1.*m3)*(m2 + m3) + (m3**2 - 1.*tman)*uman)*xv)))) + 
     &    A2sR*(A3Lc*m3*(CLoop2*(DLc*m1*(m1**2 - 1.*m4**2 - 1.*uman) + DRc*mq3c*(m1**2 + m4**2 - 1.*uman))*(-1. + m4**2*xv) + 
     &          CLoop1*(DLc*m1*(m1**2 + 2.*m2**2 - 2.*m3**2 - 1.*m4**2 + uman + (m2**2 - 1.*tman)*(m1**2 - 1.*m4**2 - 1.*uman)*xv) + 
     &             DRc*mq3c*(m1**2 - 2.*m3**2 - 1.*m4**2 + 2.*tman + uman + (m2**2 - 1.*tman)*(m1**2 + m4**2 - 1.*uman)*xv))) + 
     &       A3Rc*(-1.*CLoop2*(DRc*m1**2*(m2 - 1.*m3)*(m2 + m3) + DRc*(m3**2 - 1.*tman)*uman + DLc*m1*mq3c*(-1.*m1**2 - 1.*m2**2 + tman + uman))*(-1. + m4**2*xv) + 
     &          CLoop1*(DLc*m1*mq3c*(m1**2 - 1.*m2**2 - 2.*m3**2 - 1.*tman + uman + (m2**2 - 1.*tman)*(m1**2 + m2**2 - 1.*tman - 1.*uman)*xv) + 
     &             DRc*(m1**2*m2**2 + m1**2*m3**2 + 2.*m2**2*m3**2 - 2.*m3**4 - 2.*m3**2*m4**2 - 2.*m2**2*tman + 2.*m3**2*tman + m3**2*uman + tman*uman + 
     &                (-1.*m2**2 + tman)*(m1**2*(m2 - 1.*m3)*(m2 + m3) + (m3**2 - 1.*tman)*uman)*xv)))))/((mq2**2 - 1.*tman)*(mq3c**2 - 1.*uman)))



        end


! --- Mus  ---

       double precision function MusV(A3sL,A3sR,A3vL,A3vR,DL,DR,A1Lc,A1Rc,BLc,BRc,mq3,mq1c,xv,wquark)

        implicit none
        
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

        double complex A3sL,A3sR,A3vL,A3vR,DL,DR,A1Lc,A1Rc,BLc,BRc
        double precision mq3,mq1c,xv,wquark

        MusV = 0d0
        MusV = dreal((A1Rc*(2.*A3sL*BLc*DL*m1**3*m3 - 2.*A3sL*BLc*DL*m1*m2**2*m3 + 2.*A3sL*BLc*DL*m1*m3**3 - 2.*A3sL*BLc*DL*m1*m3*m4**2 + 4.*A3sL*BRc*DL*m1*m2**2*mq1c - 
     &       4.*A3sL*BRc*DL*m1*m3**2*mq1c + 2.*A3vL*BLc*DR*m1**2*m2**2*mq3 - 2.*A3vL*BLc*DR*m2**4*mq3 - 4.*A3sL*BLc*DR*m1**2*m3*mq3 + 2.*A3vR*BLc*DL*m1**3*m3*mq3 + 
     &       4.*A3sL*BLc*DR*m2**2*m3*mq3 - 2.*A3vR*BLc*DL*m1*m2**2*m3*mq3 - 2.*A3vL*BLc*DR*m1**2*m3**2*mq3 + 2.*A3vL*BLc*DR*m2**2*m3**2*mq3 + 2.*A3vR*BLc*DL*m1*m3**3*mq3 - 
     &       2.*A3vR*BLc*DL*m1*m3*m4**2*mq3 - 2.*A3sL*BRc*DR*m2**2*mq1c*mq3 + 4.*A3vR*BRc*DL*m1*m2**2*mq1c*mq3 + 2.*A3vL*BRc*DR*m1**2*m3*mq1c*mq3 - 
     &       4.*A3vR*BRc*DL*m1*m3**2*mq1c*mq3 - 2.*A3sL*BRc*DR*m4**2*mq1c*mq3 - 2.*A3vL*BRc*DR*m3*m4**2*mq1c*mq3 + 2.*A3vL*BLc*DR*m2**2*mq3*sman - 4.*A3sL*BLc*DR*m3*mq3*sman - 
     &       2.*A3vL*BLc*DR*m3**2*mq3*sman + 2.*A3sL*BRc*DR*mq1c*mq3*sman - 4.*A3vR*BLc*DR*m1**2*m3*uman + 4.*A3vR*BLc*DR*m2**2*m3*uman + 2.*A3vL*BLc*DL*m1*m3**2*uman - 
     &       2.*A3vL*BLc*DL*m1*m4**2*uman - 4.*A3sL*BRc*DL*m1*mq1c*uman - 2.*A3vR*BRc*DR*m2**2*mq1c*uman - 8.*A3vL*BRc*DL*m1*m3*mq1c*uman - 2.*A3vR*BRc*DR*m4**2*mq1c*uman - 
     &       2.*A3vL*BLc*DR*m1**2*mq3*uman + 2.*A3vL*BLc*DR*m2**2*mq3*uman + 2.*A3sL*BRc*DR*mq1c*mq3*uman - 4.*A3vR*BRc*DL*m1*mq1c*mq3*uman + 2.*A3vL*BRc*DR*m3*mq1c*mq3*uman + 
     &       2.*A3vL*BLc*DL*m1*sman*uman - 4.*A3vR*BLc*DR*m3*sman*uman + 2.*A3vR*BRc*DR*mq1c*sman*uman - 2.*A3vL*BLc*DR*mq3*sman*uman + 2.*A3vR*BRc*DR*mq1c*uman**2 + 
     &       (A3sL*(BRc*mq1c*(DL*m1*m4**2*(-1.*m2**2 + m3**2 + uman) - 1.*DR*mq3*(-1.*m2**2*m4**2 + m3**2*m4**2 + m1**2*(m3**2 + m4**2 - 1.*sman) - 1.*m3**2*uman + sman*uman)) + 
     &             BLc*m3*(DR*m4**2*mq3*(m1**2 - 1.*m2**2 + sman) - 1.*DL*m1*(m1**2*(m3**2 - 1.*sman) + m4**2*(-1.*m2**2 + sman) + (-1.*m3**2 + m4**2 + sman)*uman))) + 
     &          A3vR*(BRc*mq1c*(DL*m1*m4**2*mq3*(-1.*m2**2 + m3**2 + uman) + 
     &                DR*uman*((m2 - 1.*m3)*(m2 + m3)*m4**2 - 1.*m1**2*(m3**2 + m4**2 - 1.*sman) + (m3**2 - 1.*sman)*uman)) + 
     &             BLc*m3*(DR*m4**2*(m1**2 - 1.*m2**2 + sman)*uman - 1.*DL*m1*mq3*(m1**2*(m3**2 - 1.*sman) + m4**2*(-1.*m2**2 + sman) + (-1.*m3**2 + m4**2 + sman)*uman))) + 
     &          A3vL*(BRc*m3*mq1c*(2.*DL*m1*m4**2*uman + DR*mq3*(m1**4 + uman*(-1.*m4**2 + uman) - 1.*m1**2*(m4**2 + 2.*uman))) + 
     &             BLc*(DL*m1*(m3**4 + sman*(-1.*m4**2 + sman) - 1.*m3**2*(m4**2 + 2.*sman))*uman + 
     &                DR*mq3*(-1.*m1**4*m3**2 + (m2**2*(-1.*m3**2 + sman) + sman*(m3**2 + m4**2 - 1.*sman - 1.*uman))*uman + 
     &                   m1**2*(-1.*m3**4 + m2**2*(m3**2 - 1.*sman) + sman*uman + m3**2*(m4**2 + sman + uman))))))*xv + 
     &       A3sR*(2.*BLc*(DL*m1*mq3*(m3**2 - 1.*m4**2 + sman) + DR*(m1**2 - 1.*m2**2 + sman)*(m2**2 - 1.*m3**2 - 1.*uman)) + 
     &          BLc*(DL*m1*mq3*(m3**4 + sman*(-1.*m4**2 + sman) - 1.*m3**2*(m4**2 + 2.*sman)) + 
     &             DR*(-1.*m1**4*m3**2 + (m2**2*(-1.*m3**2 + sman) + sman*(m3**2 + m4**2 - 1.*sman - 1.*uman))*uman + 
     &                m1**2*(-1.*m3**4 + m2**2*(m3**2 - 1.*sman) + sman*uman + m3**2*(m4**2 + sman + uman))))*xv + 
     &          BRc*m3*mq1c*(2.*DR*(m1**2 - 1.*m4**2 + uman) + DR*(m1**4 + uman*(-1.*m4**2 + uman) - 1.*m1**2*(m4**2 + 2.*uman))*xv + 2.*DL*m1*mq3*(-4. + m4**2*xv)))) + 
     &    A1Lc*(2.*A3sR*BRc*DR*m1**3*m3 - 2.*A3sR*BRc*DR*m1*m2**2*m3 + 2.*A3sR*BRc*DR*m1*m3**3 - 2.*A3sR*BRc*DR*m1*m3*m4**2 + 4.*A3sR*BLc*DR*m1*m2**2*mq1c - 
     &       4.*A3sR*BLc*DR*m1*m3**2*mq1c + 2.*A3vR*BRc*DL*m1**2*m2**2*mq3 - 2.*A3vR*BRc*DL*m2**4*mq3 - 4.*A3sR*BRc*DL*m1**2*m3*mq3 + 2.*A3vL*BRc*DR*m1**3*m3*mq3 + 
     &       4.*A3sR*BRc*DL*m2**2*m3*mq3 - 2.*A3vL*BRc*DR*m1*m2**2*m3*mq3 - 2.*A3vR*BRc*DL*m1**2*m3**2*mq3 + 2.*A3vR*BRc*DL*m2**2*m3**2*mq3 + 2.*A3vL*BRc*DR*m1*m3**3*mq3 - 
     &       2.*A3vL*BRc*DR*m1*m3*m4**2*mq3 - 2.*A3sR*BLc*DL*m2**2*mq1c*mq3 + 4.*A3vL*BLc*DR*m1*m2**2*mq1c*mq3 + 2.*A3vR*BLc*DL*m1**2*m3*mq1c*mq3 - 
     &       4.*A3vL*BLc*DR*m1*m3**2*mq1c*mq3 - 2.*A3sR*BLc*DL*m4**2*mq1c*mq3 - 2.*A3vR*BLc*DL*m3*m4**2*mq1c*mq3 + 2.*A3vR*BRc*DL*m2**2*mq3*sman - 4.*A3sR*BRc*DL*m3*mq3*sman - 
     &       2.*A3vR*BRc*DL*m3**2*mq3*sman + 2.*A3sR*BLc*DL*mq1c*mq3*sman - 4.*A3vL*BRc*DL*m1**2*m3*uman + 4.*A3vL*BRc*DL*m2**2*m3*uman + 2.*A3vR*BRc*DR*m1*m3**2*uman - 
     &       2.*A3vR*BRc*DR*m1*m4**2*uman - 4.*A3sR*BLc*DR*m1*mq1c*uman - 2.*A3vL*BLc*DL*m2**2*mq1c*uman - 8.*A3vR*BLc*DR*m1*m3*mq1c*uman - 2.*A3vL*BLc*DL*m4**2*mq1c*uman - 
     &       2.*A3vR*BRc*DL*m1**2*mq3*uman + 2.*A3vR*BRc*DL*m2**2*mq3*uman + 2.*A3sR*BLc*DL*mq1c*mq3*uman - 4.*A3vL*BLc*DR*m1*mq1c*mq3*uman + 2.*A3vR*BLc*DL*m3*mq1c*mq3*uman + 
     &       2.*A3vR*BRc*DR*m1*sman*uman - 4.*A3vL*BRc*DL*m3*sman*uman + 2.*A3vL*BLc*DL*mq1c*sman*uman - 2.*A3vR*BRc*DL*mq3*sman*uman + 2.*A3vL*BLc*DL*mq1c*uman**2 + 
     &       (A3sR*(BLc*mq1c*(DR*m1*m4**2*(-1.*m2**2 + m3**2 + uman) - 1.*DL*mq3*(-1.*m2**2*m4**2 + m3**2*m4**2 + m1**2*(m3**2 + m4**2 - 1.*sman) - 1.*m3**2*uman + sman*uman)) + 
     &             BRc*m3*(DL*m4**2*mq3*(m1**2 - 1.*m2**2 + sman) - 1.*DR*m1*(m1**2*(m3**2 - 1.*sman) + m4**2*(-1.*m2**2 + sman) + (-1.*m3**2 + m4**2 + sman)*uman))) + 
     &          A3vL*(BLc*mq1c*(DR*m1*m4**2*mq3*(-1.*m2**2 + m3**2 + uman) + 
     &                DL*uman*((m2 - 1.*m3)*(m2 + m3)*m4**2 - 1.*m1**2*(m3**2 + m4**2 - 1.*sman) + (m3**2 - 1.*sman)*uman)) + 
     &             BRc*m3*(DL*m4**2*(m1**2 - 1.*m2**2 + sman)*uman - 1.*DR*m1*mq3*(m1**2*(m3**2 - 1.*sman) + m4**2*(-1.*m2**2 + sman) + (-1.*m3**2 + m4**2 + sman)*uman))) + 
     &          A3vR*(BLc*m3*mq1c*(2.*DR*m1*m4**2*uman + DL*mq3*(m1**4 + uman*(-1.*m4**2 + uman) - 1.*m1**2*(m4**2 + 2.*uman))) + 
     &             BRc*(DR*m1*(m3**4 + sman*(-1.*m4**2 + sman) - 1.*m3**2*(m4**2 + 2.*sman))*uman + 
     &                DL*mq3*(-1.*m1**4*m3**2 + (m2**2*(-1.*m3**2 + sman) + sman*(m3**2 + m4**2 - 1.*sman - 1.*uman))*uman + 
     &                   m1**2*(-1.*m3**4 + m2**2*(m3**2 - 1.*sman) + sman*uman + m3**2*(m4**2 + sman + uman))))))*xv + 
     &       A3sL*(2.*BRc*(DR*m1*mq3*(m3**2 - 1.*m4**2 + sman) + DL*(m1**2 - 1.*m2**2 + sman)*(m2**2 - 1.*m3**2 - 1.*uman)) + 
     &          BRc*(DR*m1*mq3*(m3**4 + sman*(-1.*m4**2 + sman) - 1.*m3**2*(m4**2 + 2.*sman)) + 
     &             DL*(-1.*m1**4*m3**2 + (m2**2*(-1.*m3**2 + sman) + sman*(m3**2 + m4**2 - 1.*sman - 1.*uman))*uman + 
     &                m1**2*(-1.*m3**4 + m2**2*(m3**2 - 1.*sman) + sman*uman + m3**2*(m4**2 + sman + uman))))*xv + 
     &          BLc*m3*mq1c*(2.*DL*(m1**2 - 1.*m4**2 + uman) + DL*(m1**4 + uman*(-1.*m4**2 + uman) - 1.*m1**2*(m4**2 + 2.*uman))*xv + 2.*DR*m1*mq3*(-4. + m4**2*xv)))))/
     &  ((mq1c**2 - 1.*sman + cI*mq1c*wquark)*(mq3**2 - 1.*uman)))

        end




! --- Mut  ---

        double precision function MutV(A3sL,A3sR,A3vL,A3vR,DL,DR,A2Lc,A2Rc,CLc,mq3,mq2c,xv)

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

        double complex A3sL,A3sR,A3vL,A3vR,DL,DR,A2Lc,A2Rc,CLc
        double precision mq3,mq2c,xv

        MutV = 0d0
        MutV = dreal((CLc*(A2Lc*(A3sR*DR*m1**3*m3 + 2.*A3sR*DR*m1*m2**2*m3 - 2.*A3sR*DR*m1*m3**3 - 1.*A3sR*DR*m1*m3*m4**2 + A3vR*DL*m1**2*m2**2*mq3 + A3sR*DL*m1**2*m3*mq3 + 
     &         A3vL*DR*m1**3*m3*mq3 + 2.*A3vL*DR*m1*m2**2*m3*mq3 + A3vR*DL*m1**2*m3**2*mq3 + 2.*A3vR*DL*m2**2*m3**2*mq3 - 2.*A3sR*DL*m3**3*mq3 - 2.*A3vL*DR*m1*m3**3*mq3 - 
     &         2.*A3vR*DL*m3**4*mq3 - 1.*A3sR*DL*m3*m4**2*mq3 - 1.*A3vL*DR*m1*m3*m4**2*mq3 - 2.*A3vR*DL*m3**2*m4**2*mq3 - 2.*A3vR*DL*m2**2*mq3*tman + 2.*A3sR*DL*m3*mq3*tman + 
     &         2.*A3vR*DL*m3**2*mq3*tman + A3vR*DR*m1**3*uman - 1.*A3vR*DR*m1*m2**2*uman + A3sR*DR*m1*m3*uman + A3vL*DL*m1**2*m3*uman - 2.*A3vR*DR*m1*m3**2*uman - 
     &         2.*A3vL*DL*m3**3*uman - 1.*A3vL*DL*m3*m4**2*uman + A3sR*DL*m3*mq3*uman + A3vL*DR*m1*m3*mq3*uman + A3vR*DL*m3**2*mq3*uman - 1.*A3vR*DR*m1*tman*uman + 
     &         2.*A3vL*DL*m3*tman*uman + A3vR*DL*mq3*tman*uman + A3vR*DR*m1*uman**2 + A3vL*DL*m3*uman**2 + 
     &         (m2**2 - 1.*tman)*(m1*(A3vR*DL*m1*(-1.*m2**2 + m3**2) + A3vL*DR*m3*(m1 - 1.*m4)*(m1 + m4))*mq3 + 
     &            A3sR*m3*(DR*m1*(m1**2 - 1.*m4**2 - 1.*uman) + DL*mq3*(m1**2 + m4**2 - 1.*uman)) + 
     &            (A3vL*m3*(DL*(m1**2 + m4**2) - 1.*DR*m1*mq3) + A3vR*(DR*m1*(m1**2 + m2**2 - 1.*tman) + DL*mq3*(-1.*m3**2 + tman)))*uman - 1.*(A3vR*DR*m1 + A3vL*DL*m3)*uman**2)*xv
     &           + A3sL*(DR*m1*mq3*(m1**2 - 1.*m2**2 - 2.*m3**2 - 1.*tman + uman + (m2**2 - 1.*tman)*(m1**2 + m2**2 - 1.*tman - 1.*uman)*xv) + 
     &            DL*(-2.*m3**4 - 2.*m3**2*m4**2 + 2.*m3**2*tman + m3**2*uman + tman*uman + (m3**2 - 1.*tman)*tman*uman*xv + 
     &               m1**2*(m2**2 + m3**2 - 1.*(m2 - 1.*m3)*(m2 + m3)*(m2**2 - 1.*tman)*xv) - 1.*m2**2*(m3**2 - 1.*tman)*(-2. + uman*xv)))) + 
     &      A2Rc*(A3sL*DL*m1**3*m3 + 2.*A3sL*DL*m1*m2**2*m3 - 2.*A3sL*DL*m1*m3**3 - 1.*A3sL*DL*m1*m3*m4**2 + A3vL*DR*m1**2*m2**2*mq3 + A3sL*DR*m1**2*m3*mq3 + 
     &         A3vR*DL*m1**3*m3*mq3 + 2.*A3vR*DL*m1*m2**2*m3*mq3 + A3vL*DR*m1**2*m3**2*mq3 + 2.*A3vL*DR*m2**2*m3**2*mq3 - 2.*A3sL*DR*m3**3*mq3 - 2.*A3vR*DL*m1*m3**3*mq3 - 
     &         2.*A3vL*DR*m3**4*mq3 - 1.*A3sL*DR*m3*m4**2*mq3 - 1.*A3vR*DL*m1*m3*m4**2*mq3 - 2.*A3vL*DR*m3**2*m4**2*mq3 - 2.*A3vL*DR*m2**2*mq3*tman + 2.*A3sL*DR*m3*mq3*tman + 
     &         2.*A3vL*DR*m3**2*mq3*tman + A3vL*DL*m1**3*uman - 1.*A3vL*DL*m1*m2**2*uman + A3sL*DL*m1*m3*uman + A3vR*DR*m1**2*m3*uman - 2.*A3vL*DL*m1*m3**2*uman - 
     &         2.*A3vR*DR*m3**3*uman - 1.*A3vR*DR*m3*m4**2*uman + A3sL*DR*m3*mq3*uman + A3vR*DL*m1*m3*mq3*uman + A3vL*DR*m3**2*mq3*uman - 1.*A3vL*DL*m1*tman*uman + 
     &         2.*A3vR*DR*m3*tman*uman + A3vL*DR*mq3*tman*uman + A3vL*DL*m1*uman**2 + A3vR*DR*m3*uman**2 + 
     &         (m2**2 - 1.*tman)*(m1*(A3vL*DR*m1*(-1.*m2**2 + m3**2) + A3vR*DL*m3*(m1 - 1.*m4)*(m1 + m4))*mq3 + 
     &            A3sL*m3*(DL*m1*(m1**2 - 1.*m4**2 - 1.*uman) + DR*mq3*(m1**2 + m4**2 - 1.*uman)) + 
     &            (A3vR*m3*(DR*(m1**2 + m4**2) - 1.*DL*m1*mq3) + A3vL*(DL*m1*(m1**2 + m2**2 - 1.*tman) + DR*mq3*(-1.*m3**2 + tman)))*uman - 1.*(A3vL*DL*m1 + A3vR*DR*m3)*uman**2)*xv
     &           + A3sR*(DL*m1*mq3*(m1**2 - 1.*m2**2 - 2.*m3**2 - 1.*tman + uman + (m2**2 - 1.*tman)*(m1**2 + m2**2 - 1.*tman - 1.*uman)*xv) + 
     &            DR*(-2.*m3**4 - 2.*m3**2*m4**2 + 2.*m3**2*tman + m3**2*uman + tman*uman + (m3**2 - 1.*tman)*tman*uman*xv + 
     &               m1**2*(m2**2 + m3**2 - 1.*(m2 - 1.*m3)*(m2 + m3)*(m2**2 - 1.*tman)*xv) - 1.*m2**2*(m3**2 - 1.*tman)*(-2. + uman*xv))))))/
     &  ((mq2c**2 - 1.*tman)*(mq3**2 - 1.*uman)))



       end


! --- Muu  ---

        double precision function MuuV(A3sL,A3sR,A3vL,A3vR,DL,DR,A3Lc,A3Rc,DLc,DRc,mq3,mq3c,xv)

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"

        double complex A3sL,A3sR,A3vL,A3vR,DL,DR,A3Lc,A3Rc,DLc,DRc
        double precision mq3,mq3c,xv
        
        MuuV = 0d0
        MuuV = dreal((-1.*(A3Rc*(A3sR*(DR*(DLc*m1*mq3c*(m2**2 - 1.*m3**2 - 1.*uman)*(-4. + m4**2*xv) + 
     &               DRc*(2.*(m2 - 1.*m3)*(m2 + m3)*(m1 - 1.*m4)*(m1 + m4) - 2.*(m1**2 + m3**2 - 1.*sman)*uman + 
     &                  (m1**2*(m2 - 1.*m3)*(m2 + m3)*(m1 - 1.*m4)*(m1 + m4) + m1**2*(-1.*m1**2 - 2.*m2**2 + m3**2 + sman)*uman + (2.*m1**2 + m2**2 + m4**2 - 1.*sman)*uman**2 - 
     &                     1.*uman**3)*xv)) + DL*mq3*(DRc*m1*(m2**2 - 1.*m3**2 - 1.*uman)*(-4. + m4**2*xv) + 
     &               DLc*mq3c*(-2.*(sman + uman) + (m3**2 - 1.*sman)*(m1**2 - 1.*uman)*xv + m4**2*(2. + (m1**2 + m3**2)*xv) + m2**2*(2. - 1.*m4**2*xv)))) + 
     &         A3vL*(DR*mq3*(DLc*m1*mq3c*(m2**2 - 1.*m3**2 - 1.*uman)*(-4. + m4**2*xv) + 
     &               DRc*(2.*(m2 - 1.*m3)*(m2 + m3)*(m1 - 1.*m4)*(m1 + m4) - 2.*(m1**2 + m3**2 - 1.*sman)*uman + 
     &                  (m1**2*(m2 - 1.*m3)*(m2 + m3)*(m1 - 1.*m4)*(m1 + m4) + m1**2*(-1.*m1**2 - 2.*m2**2 + m3**2 + sman)*uman + (2.*m1**2 + m2**2 + m4**2 - 1.*sman)*uman**2 - 
     &                     1.*uman**3)*xv)) + DL*uman*(DRc*m1*(m2**2 - 1.*m3**2 - 1.*uman)*(-4. + m4**2*xv) + 
     &               DLc*mq3c*(-2.*(sman + uman) + (m3**2 - 1.*sman)*(m1**2 - 1.*uman)*xv + m4**2*(2. + (m1**2 + m3**2)*xv) + m2**2*(2. - 1.*m4**2*xv)))) + 
     &         m3*(A3sL*(DR*mq3*(-2.*DLc*m1*mq3c*(-4. + m4**2*xv) + DRc*(-1.*m1**4*xv + (m4**2 - 1.*uman)*(2. + uman*xv) + m1**2*(-2. + m4**2*xv + 2.*uman*xv))) + 
     &               DL*(-2.*DRc*m1*uman*(-4. + m4**2*xv) + DLc*mq3c*(-1.*m1**4*xv + (m4**2 - 1.*uman)*(2. + uman*xv) + m1**2*(-2. + m4**2*xv + 2.*uman*xv)))) + 
     &            A3vR*(DR*uman*(-2.*DLc*m1*mq3c*(-4. + m4**2*xv) + DRc*(-1.*m1**4*xv + (m4**2 - 1.*uman)*(2. + uman*xv) + m1**2*(-2. + m4**2*xv + 2.*uman*xv))) + 
     &               DL*mq3*(-2.*DRc*m1*uman*(-4. + m4**2*xv) + DLc*mq3c*(-1.*m1**4*xv + (m4**2 - 1.*uman)*(2. + uman*xv) + m1**2*(-2. + m4**2*xv + 2.*uman*xv)))))) + 
     &      A3Lc*(A3sL*(DL*(DRc*m1*mq3c*(m2**2 - 1.*m3**2 - 1.*uman)*(-4. + m4**2*xv) + 
     &               DLc*(2.*(m2 - 1.*m3)*(m2 + m3)*(m1 - 1.*m4)*(m1 + m4) - 2.*(m1**2 + m3**2 - 1.*sman)*uman + 
     &                  (m1**2*(m2 - 1.*m3)*(m2 + m3)*(m1 - 1.*m4)*(m1 + m4) + m1**2*(-1.*m1**2 - 2.*m2**2 + m3**2 + sman)*uman + (2.*m1**2 + m2**2 + m4**2 - 1.*sman)*uman**2 - 
     &                     1.*uman**3)*xv)) + DR*mq3*(DLc*m1*(m2**2 - 1.*m3**2 - 1.*uman)*(-4. + m4**2*xv) + 
     &               DRc*mq3c*(-2.*(sman + uman) + (m3**2 - 1.*sman)*(m1**2 - 1.*uman)*xv + m4**2*(2. + (m1**2 + m3**2)*xv) + m2**2*(2. - 1.*m4**2*xv)))) + 
     &         A3vR*(DL*mq3*(DRc*m1*mq3c*(m2**2 - 1.*m3**2 - 1.*uman)*(-4. + m4**2*xv) + 
     &               DLc*(2.*(m2 - 1.*m3)*(m2 + m3)*(m1 - 1.*m4)*(m1 + m4) - 2.*(m1**2 + m3**2 - 1.*sman)*uman + 
     &                  (m1**2*(m2 - 1.*m3)*(m2 + m3)*(m1 - 1.*m4)*(m1 + m4) + m1**2*(-1.*m1**2 - 2.*m2**2 + m3**2 + sman)*uman + (2.*m1**2 + m2**2 + m4**2 - 1.*sman)*uman**2 - 
     &                     1.*uman**3)*xv)) + DR*uman*(DLc*m1*(m2**2 - 1.*m3**2 - 1.*uman)*(-4. + m4**2*xv) + 
     &               DRc*mq3c*(-2.*(sman + uman) + (m3**2 - 1.*sman)*(m1**2 - 1.*uman)*xv + m4**2*(2. + (m1**2 + m3**2)*xv) + m2**2*(2. - 1.*m4**2*xv)))) + 
     &         m3*(A3sR*(DL*mq3*(-2.*DRc*m1*mq3c*(-4. + m4**2*xv) + DLc*(-1.*m1**4*xv + (m4**2 - 1.*uman)*(2. + uman*xv) + m1**2*(-2. + m4**2*xv + 2.*uman*xv))) + 
     &               DR*(-2.*DLc*m1*uman*(-4. + m4**2*xv) + DRc*mq3c*(-1.*m1**4*xv + (m4**2 - 1.*uman)*(2. + uman*xv) + m1**2*(-2. + m4**2*xv + 2.*uman*xv)))) + 
     &            A3vL*(DL*uman*(-2.*DRc*m1*mq3c*(-4. + m4**2*xv) + DLc*(-1.*m1**4*xv + (m4**2 - 1.*uman)*(2. + uman*xv) + m1**2*(-2. + m4**2*xv + 2.*uman*xv))) + 
     &               DR*mq3*(-2.*DLc*m1*uman*(-4. + m4**2*xv) + DRc*mq3c*(-1.*m1**4*xv + (m4**2 - 1.*uman)*(2. + uman*xv) + m1**2*(-2. + m4**2*xv + 2.*uman*xv))))))))/
     &  ((mq3**2 - 1.*uman)*(mq3c**2 - 1.*uman)))


        end

