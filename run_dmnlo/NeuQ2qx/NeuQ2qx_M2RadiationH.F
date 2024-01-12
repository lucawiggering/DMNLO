! =================================================================
	double precision function NeuQ2qx_M2RadiationHSS()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double precision radSiSi, radSiSp, radSiSf, radSpSp, radSpSf, radSfSf
        double complex A1L, A1Lc, A1R, A1Rc, F1L, F1R, F1Lc, F1Rc
        double complex g0, g1
        double precision mq1, mq1c


        NeuQ2qx_M2RadiationHSS = 0d0

        radSiSi = 0d0
        radSiSp = 0d0
        radSiSf = 0d0
        radSpSp = 0d0
        radSpSf = 0d0
        radSfSf = 0d0

        A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A1R = AR(chimeigen,sfmeigen,sftype,sfgen)
        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))

	F1L = FL(xtype,pstype,psgen)
	F1R = FR(xtype,pstype,psgen)
	F1Lc = dconjg(FL(xtype,pstype,psgen))
	F1Rc = dconjg(FR(xtype,pstype,psgen))

        g0 = -gs
	g1 = -gs
        
        mq1 = Mf(pstype,psgen) 
        mq1c = Mf(pstype,psgen) 


	radSiSi = dreal((-2*g0**2*(m1**2 + m2**2 + 2*(pbp1 + pbp2) - sman)*(A1R*(A1Lc*m1*(F1R*mq1*(2*F1Lc*m3*mq1c + F1Rc*(m3**2 - m4**2 + sman*(-1 + x1 + x2))) + 
     &            F1L*(2*F1Rc*m3*sman*(-1 + x1 + x2) + F1Lc*mq1c*(m3**2 - m4**2 + sman*(-1 + x1 + x2)))) + 
     &         2*A1Rc*(F1R*mq1*(F1Rc*mq1c*pap1 + F1Lc*m3*(pap1 + pap2)) + F1L*(F1Rc*m3*mq1c*(pap1 + pap2) + F1Lc*(m3**2*(pap1 + pap2) - m4**2*(pap1 + pap2) + pap2*sman*(-1 + x1 + x2))))) + 
     &      A1L*(A1Rc*m1*(F1L*mq1*(2*F1Rc*m3*mq1c + F1Lc*(m3**2 - m4**2 + sman*(-1 + x1 + x2))) + F1R*(2*F1Lc*m3*sman*(-1 + x1 + x2) + F1Rc*mq1c*(m3**2 - m4**2 + sman*(-1 + x1 + x2)))) + 
     &         2*A1Lc*(F1L*mq1*(F1Lc*mq1c*pap1 + F1Rc*m3*(pap1 + pap2)) + F1R*(F1Lc*m3*mq1c*(pap1 + pap2) + F1Rc*(m3**2*(pap1 + pap2) - m4**2*(pap1 + pap2) + pap2*sman*(-1 + x1 + x2)))))))/
     &  ((-m1**2 + m2**2 - 2*(pbp1 + pbp2) + sman)**2*(-mq1**2 + sman*(-1 + x1 + x2))*(-mq1c**2 + sman*(-1 + x1 + x2))))


        radSpSp = dreal((2*g1**2*(A1R*(-2*A1Lc*m1*(F1L*(F1Lc*mq1c*(-(mq1**2*pap1) - mq1**2*pbp1 + m4**2*(-mq1**2 + pap1 + pap2 + pbp1 + pbp2 - sman) - mq1**2*sman + pap2*sman + pbp2*sman - sman**2 + 
     &                  m3**2*(mq1**2 - pap1 - pap2 - pbp1 - pbp2 + sman) + mq1**2*sman*x1 - pap2*sman*x1 - pbp2*sman*x1 + sman**2*x1 + sman*(mq1**2 - pap2 - pbp2 + sman)*x2) - 
     &               F1Rc*m3*(mq1c**2*(pap1 + pap2 + pbp1 + pbp2) + mq1**2*(-2*mq1c**2 + pap1 + pap2 + pbp1 + pbp2) - 2*sman**2*(-1 + x1 + x2))) + 
     &            F1R*mq1*(F1Rc*(-(mq1c**2*pap1) - mq1c**2*pbp1 + m4**2*(-mq1c**2 + pap1 + pap2 + pbp1 + pbp2 - sman) - mq1c**2*sman + pap2*sman + pbp2*sman - sman**2 + m3**2*(mq1c**2 - pap1 - pap2 - pbp1 - pbp2 + sman) + 
     &                  mq1c**2*sman*x1 - pap2*sman*x1 - pbp2*sman*x1 + sman**2*x1 + sman*(mq1c**2 - pap2 - pbp2 + sman)*x2) - 2*F1Lc*m3*mq1c*(pap1 + pap2 + pbp1 + pbp2 - sman*(x1 + x2)))) + 
     &         A1Rc*(-2*F1R*mq1*(F1Rc*mq1c*(-(m3**2*pap1) + m4**2*pap1 - m3**2*pap2 + m4**2*pap2 + m3**2*sman - m4**2*sman + pap2*sman - pbp1*sman - sman**2 - pap2*sman*x1 + sman**2*x1 - pap2*sman*x2 + sman**2*x2 + 
     &                  m2**2*(-m3**2 + m4**2 + pap1 + pbp1 - sman*(-1 + x1 + x2)) + m1**2*(m3**2 - m4**2 - pap1 - pbp1 + sman*(-1 + x1 + x2))) + 
     &               F1Lc*m3*(-(mq1c**2*pap1) - mq1c**2*pap2 + mq1c**2*sman - pbp1*sman - pbp2*sman - sman**2 + sman**2*x1 + sman**2*x2 + m2**2*(-mq1c**2 + pap1 + pap2 + pbp1 + pbp2 - sman*(-1 + x1 + x2)) + 
     &                  m1**2*(mq1c**2 - pap1 - pap2 - pbp1 - pbp2 + sman*(-1 + x1 + x2)))) + 
     &            F1L*(-2*F1Rc*m3*mq1c*(-(mq1**2*pap1) - mq1**2*pap2 + mq1**2*sman - pbp1*sman - pbp2*sman - sman**2 + sman**2*x1 + sman**2*x2 + m2**2*(-mq1**2 + pap1 + pap2 + pbp1 + pbp2 - sman*(-1 + x1 + x2)) + 
     &                  m1**2*(mq1**2 - pap1 - pap2 - pbp1 - pbp2 + sman*(-1 + x1 + x2))) - 
     &               F1Lc*(-2*mq1**2*mq1c**2*pap1 + m3**2*mq1**2*sman - m4**2*mq1**2*sman + m3**2*mq1c**2*sman - m4**2*mq1c**2*sman - 2*m3**2*pbp1*sman + 2*m4**2*pbp1*sman - 2*m3**2*pbp2*sman + 2*m4**2*pbp2*sman - 
     &                  mq1**2*sman**2 - mq1c**2*sman**2 + 2*pbp2*sman**2 + mq1**2*sman**2*x1 + mq1c**2*sman**2*x1 - 2*pbp2*sman**2*x1 + (mq1**2 + mq1c**2 - 2*pbp2)*sman**2*x2 + 
     &                  m1**2*(m3**2*(mq1**2 + mq1c**2 - 2*(pap1 + pap2 + pbp1 + pbp2)) - m4**2*(mq1**2 + mq1c**2 - 2*(pap1 + pap2 + pbp1 + pbp2)) + (mq1**2 + mq1c**2 - 2*(pap2 + pbp2))*sman*(-1 + x1 + x2)) - 
     &                  m2**2*(m3**2*(mq1**2 + mq1c**2 - 2*(pap1 + pap2 + pbp1 + pbp2)) - m4**2*(mq1**2 + mq1c**2 - 2*(pap1 + pap2 + pbp1 + pbp2)) + (mq1**2 + mq1c**2 - 2*(pap2 + pbp2))*sman*(-1 + x1 + x2)))))) + 
     &      A1L*(-2*A1Rc*m1*(F1R*(F1Rc*mq1c*(-(mq1**2*pap1) - mq1**2*pbp1 + m4**2*(-mq1**2 + pap1 + pap2 + pbp1 + pbp2 - sman) - mq1**2*sman + pap2*sman + pbp2*sman - sman**2 + 
     &                  m3**2*(mq1**2 - pap1 - pap2 - pbp1 - pbp2 + sman) + mq1**2*sman*x1 - pap2*sman*x1 - pbp2*sman*x1 + sman**2*x1 + sman*(mq1**2 - pap2 - pbp2 + sman)*x2) - 
     &               F1Lc*m3*(mq1c**2*(pap1 + pap2 + pbp1 + pbp2) + mq1**2*(-2*mq1c**2 + pap1 + pap2 + pbp1 + pbp2) - 2*sman**2*(-1 + x1 + x2))) + 
     &            F1L*mq1*(F1Lc*(-(mq1c**2*pap1) - mq1c**2*pbp1 + m4**2*(-mq1c**2 + pap1 + pap2 + pbp1 + pbp2 - sman) - mq1c**2*sman + pap2*sman + pbp2*sman - sman**2 + m3**2*(mq1c**2 - pap1 - pap2 - pbp1 - pbp2 + sman) + 
     &                  mq1c**2*sman*x1 - pap2*sman*x1 - pbp2*sman*x1 + sman**2*x1 + sman*(mq1c**2 - pap2 - pbp2 + sman)*x2) - 2*F1Rc*m3*mq1c*(pap1 + pap2 + pbp1 + pbp2 - sman*(x1 + x2)))) + 
     &         A1Lc*(-2*F1L*mq1*(F1Lc*mq1c*(-(m3**2*pap1) + m4**2*pap1 - m3**2*pap2 + m4**2*pap2 + m3**2*sman - m4**2*sman + pap2*sman - pbp1*sman - sman**2 - pap2*sman*x1 + sman**2*x1 - pap2*sman*x2 + sman**2*x2 + 
     &                  m2**2*(-m3**2 + m4**2 + pap1 + pbp1 - sman*(-1 + x1 + x2)) + m1**2*(m3**2 - m4**2 - pap1 - pbp1 + sman*(-1 + x1 + x2))) + 
     &               F1Rc*m3*(-(mq1c**2*pap1) - mq1c**2*pap2 + mq1c**2*sman - pbp1*sman - pbp2*sman - sman**2 + sman**2*x1 + sman**2*x2 + m2**2*(-mq1c**2 + pap1 + pap2 + pbp1 + pbp2 - sman*(-1 + x1 + x2)) + 
     &                  m1**2*(mq1c**2 - pap1 - pap2 - pbp1 - pbp2 + sman*(-1 + x1 + x2)))) + 
     &            F1R*(-2*F1Lc*m3*mq1c*(-(mq1**2*pap1) - mq1**2*pap2 + mq1**2*sman - pbp1*sman - pbp2*sman - sman**2 + sman**2*x1 + sman**2*x2 + m2**2*(-mq1**2 + pap1 + pap2 + pbp1 + pbp2 - sman*(-1 + x1 + x2)) + 
     &                  m1**2*(mq1**2 - pap1 - pap2 - pbp1 - pbp2 + sman*(-1 + x1 + x2))) - 
     &               F1Rc*(-2*mq1**2*mq1c**2*pap1 + m3**2*mq1**2*sman - m4**2*mq1**2*sman + m3**2*mq1c**2*sman - m4**2*mq1c**2*sman - 2*m3**2*pbp1*sman + 2*m4**2*pbp1*sman - 2*m3**2*pbp2*sman + 2*m4**2*pbp2*sman - 
     &                  mq1**2*sman**2 - mq1c**2*sman**2 + 2*pbp2*sman**2 + mq1**2*sman**2*x1 + mq1c**2*sman**2*x1 - 2*pbp2*sman**2*x1 + (mq1**2 + mq1c**2 - 2*pbp2)*sman**2*x2 + 
     &                  m1**2*(m3**2*(mq1**2 + mq1c**2 - 2*(pap1 + pap2 + pbp1 + pbp2)) - m4**2*(mq1**2 + mq1c**2 - 2*(pap1 + pap2 + pbp1 + pbp2)) + (mq1**2 + mq1c**2 - 2*(pap2 + pbp2))*sman*(-1 + x1 + x2)) - 
     &                  m2**2*(m3**2*(mq1**2 + mq1c**2 - 2*(pap1 + pap2 + pbp1 + pbp2)) - m4**2*(mq1**2 + mq1c**2 - 2*(pap1 + pap2 + pbp1 + pbp2)) + (mq1**2 + mq1c**2 - 2*(pap2 + pbp2))*sman*(-1 + x1 + x2))))))))/
     &  ((mq1**2 - sman)*(-mq1c**2 + sman)*(mq1**2 - sman*(-1 + x1 + x2))*(-mq1c**2 + sman*(-1 + x1 + x2))))


        radSfSf = dreal((2*g1**2*(A1R*(2*A1Lc*m1*(F1Lc*mq1c*(F1R*m3*mq1*(m3**2 + m4**2 + sman - sman*x2) + F1L*(-(m3**2*(pap1 + 3*pap2 + pbp1 + 3*pbp2 - 3*sman)) + (pap1 + pap2 + pbp1 + pbp2 - sman)*(m4**2 + sman - sman*x2))) + 
     &            F1Rc*(F1L*m3*sman*(m3**2 + m4**2 + sman - sman*x2) + F1R*mq1*(-(m3**2*(pap1 + 3*pap2 + pbp1 + 3*pbp2 - 3*sman)) + (pap1 + pap2 + pbp1 + pbp2 - sman)*(m4**2 + sman - sman*x2)))) + 
     &         A1Rc*(F1Rc*mq1c*(F1L*m3*(m1**2 - m2**2 + sman)*(m3**2 + m4**2 + sman - sman*x2) + 
     &               F1R*mq1*(-2*m3**2*pap1 + 2*m4**2*pap1 - 6*m3**2*pap2 + 2*m4**2*pap2 + 3*m3**2*sman - m4**2*sman + 2*pap1*sman + 2*pap2*sman - sman**2 + m1**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) + 
     &                  sman*(-2*(pap1 + pap2) + sman)*x2 + m2**2*(-3*m3**2 + m4**2 + sman - sman*x2))) + 
     &            F1Lc*(F1R*m3*mq1*(m1**2 - m2**2 + sman)*(m3**2 + m4**2 + sman - sman*x2) - 
     &               F1L*(sman*(m3**2*(2*pbp1 + 6*pbp2 - 3*sman) - (2*(pbp1 + pbp2) - sman)*(m4**2 + sman - sman*x2)) + 
     &                  m1**2*(m3**2*(2*(pap1 + 3*pap2 + pbp1 + 3*pbp2) - 3*sman) - (2*(pap1 + pap2 + pbp1 + pbp2) - sman)*(m4**2 + sman - sman*x2)) + 
     &                  m2**2*(m3**2*(-2*(pap1 + 3*pap2 + pbp1 + 3*pbp2) + 3*sman) + (2*(pap1 + pap2 + pbp1 + pbp2) - sman)*(m4**2 + sman - sman*x2)))))) + 
     &      A1L*(2*A1Rc*m1*(F1L*mq1*(F1Rc*m3*mq1c*(m3**2 + m4**2 + sman - sman*x2) + F1Lc*(-(m3**2*(pap1 + 3*pap2 + pbp1 + 3*pbp2 - 3*sman)) + (pap1 + pap2 + pbp1 + pbp2 - sman)*(m4**2 + sman - sman*x2))) + 
     &            F1R*(F1Lc*m3*sman*(m3**2 + m4**2 + sman - sman*x2) + F1Rc*mq1c*(-(m3**2*(pap1 + 3*pap2 + pbp1 + 3*pbp2 - 3*sman)) + (pap1 + pap2 + pbp1 + pbp2 - sman)*(m4**2 + sman - sman*x2)))) + 
     &         A1Lc*(F1L*mq1*(F1Rc*m3*(m1**2 - m2**2 + sman)*(m3**2 + m4**2 + sman - sman*x2) + 
     &               F1Lc*mq1c*(-2*m3**2*pap1 + 2*m4**2*pap1 - 6*m3**2*pap2 + 2*m4**2*pap2 + 3*m3**2*sman - m4**2*sman + 2*pap1*sman + 2*pap2*sman - sman**2 + m1**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) + 
     &                  sman*(-2*(pap1 + pap2) + sman)*x2 + m2**2*(-3*m3**2 + m4**2 + sman - sman*x2))) + 
     &            F1R*(F1Lc*m3*mq1c*(m1**2 - m2**2 + sman)*(m3**2 + m4**2 + sman - sman*x2) - 
     &               F1Rc*(sman*(m3**2*(2*pbp1 + 6*pbp2 - 3*sman) - (2*(pbp1 + pbp2) - sman)*(m4**2 + sman - sman*x2)) + 
     &                  m1**2*(m3**2*(2*(pap1 + 3*pap2 + pbp1 + 3*pbp2) - 3*sman) - (2*(pap1 + pap2 + pbp1 + pbp2) - sman)*(m4**2 + sman - sman*x2)) + 
     &                  m2**2*(m3**2*(-2*(pap1 + 3*pap2 + pbp1 + 3*pbp2) + 3*sman) + (2*(pap1 + pap2 + pbp1 + pbp2) - sman)*(m4**2 + sman - sman*x2))))))))/
     &  ((mq1**2 - sman)*(-mq1c**2 + sman)*(m3**2 - m4**2 + sman*(-1 + x2))**2))


	radSiSp = dreal((g0*g1*(A1L*(2*A1Rc*m1*(F1L*mq1*(F1Rc*m3*mq1c*(-2*m1**2 + 2*m2**2 + 2*pap1 + 2*pap2 + 6*pbp1 + 6*pbp2 + sman*(-2 + x1 + x2)) + 
     &               F1Lc*(m3**2*mq1c**2 - m4**2*mq1c**2 + 4*mq1c**2*pbp1 - 4*pap2*pbp1 + 4*pap1*pbp2 - mq1c**2*sman - 2*pap1*sman - 2*pbp1*sman + 2*pap1*sman*x1 + pap2*sman*x1 + 2*pbp1*sman*x1 + pbp2*sman*x1 + 
     &                  (mq1c**2 + pap1 + pbp1)*sman*x2 + m1**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + m2**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2)))) + 
     &            F1R*(F1Lc*m3*(-2*(m1**2 - m2**2 - pap1 - pap2 - pbp1 - pbp2)*sman*(-1 + x1 + x2) + mq1c**2*(4*pbp1 + 4*pbp2 + sman*(-2 + x1 + x2))) + 
     &               F1Rc*mq1c*(2*m3**2*pap1 - 2*m4**2*pap1 + 2*m3**2*pap2 - 2*m4**2*pap2 + 6*m3**2*pbp1 - 6*m4**2*pbp1 + 4*pap2*pbp1 + 6*m3**2*pbp2 - 6*m4**2*pbp2 - 4*pap1*pbp2 - m3**2*sman + m4**2*sman - 2*pap2*sman - 
     &                  6*pbp2*sman + sman**2 + pap2*sman*x1 + 5*pbp2*sman*x1 - 2*sman**2*x1 + sman**2*x1**2 + sman*(pap1 + 2*pap2 + pbp1 + 6*pbp2 + sman*(-1 + x1))*x2 + m1**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + 
     &                  m2**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2))))) + A1Lc*(F1L*mq1*(F1Rc*m3*
     &                (4*mq1c**2*pap1 + 4*mq1c**2*pap2 + 2*mq1c**2*sman - 2*pap1*sman - 2*pap2*sman + 2*pbp1*sman + 2*pbp2*sman - 2*sman**2 + sman**2*x1 + sman**2*x2 - 
     &                  m2**2*(2*mq1c**2 - 6*pap1 - 6*pap2 + 2*pbp1 + 2*pbp2 + sman*(-2 + x1 + x2)) + m1**2*(-6*mq1c**2 + 2*pap1 + 2*pap2 + 10*pbp1 + 10*pbp2 + sman*(-2 + x1 + x2))) + 
     &               F1Lc*mq1c*(-12*pap2*pbp1 + 12*pap1*pbp2 + 2*(m3**2 - m4**2 - 3*pap1)*sman + 2*pbp1*sman - 2*sman**2 + 4*pap1*sman*x1 + 2*pap2*sman*x1 + sman**2*x1 + 
     &                  m1**2*(-2*m3**2 + 2*m4**2 + 2*pap1 + 10*pbp1 + sman*(2 - 3*x1 - 2*x2)) + 2*sman*(pap1 + sman)*x2 - m2**2*(2*m3**2 - 2*m4**2 - 6*pap1 + 2*pbp1 + sman*(-2 + x1 + 2*x2)))) + 
     &            F1R*(F1Rc*(4*m3**2*mq1c**2*pap1 - 4*m4**2*mq1c**2*pap1 + 4*m3**2*mq1c**2*pap2 - 4*m4**2*mq1c**2*pap2 + 8*mq1c**2*pap2*pbp1 - 8*mq1c**2*pap1*pbp2 + m3**2*mq1c**2*sman - m4**2*mq1c**2*sman - 
     &                  2*m3**2*pap1*sman + 2*m4**2*pap1*sman - 2*m3**2*pap2*sman + 2*m4**2*pap2*sman - 4*mq1c**2*pap2*sman + 2*m3**2*pbp1*sman - 2*m4**2*pbp1*sman - 4*pap2*pbp1*sman + 2*m3**2*pbp2*sman - 
     &                  2*m4**2*pbp2*sman + 4*pap1*pbp2*sman - m3**2*sman**2 + m4**2*sman**2 - mq1c**2*sman**2 + 2*pap2*sman**2 - 2*pbp2*sman**2 + sman**3 + 2*mq1c**2*pap2*sman*x1 + 4*pap2*pbp1*sman*x1 - 
     &                  4*pap1*pbp2*sman*x1 + mq1c**2*sman**2*x1 - 2*pap2*sman**2*x1 + 2*pbp2*sman**2*x1 - 2*sman**3*x1 + sman**3*x1**2 + 
     &                  sman*(-4*pap1*pbp2 + pap2*(4*pbp1 - 2*sman) + mq1c**2*(2*pap1 + 4*pap2 + sman) + sman*(2*pbp2 + sman*(-1 + x1)))*x2 - 
     &                  m2**2*(m3**2*(mq1c**2 - 6*pap1 - 6*pap2 + 2*(pbp1 + pbp2) - sman) + m4**2*(-mq1c**2 + 6*pap1 + 6*pap2 - 2*(pbp1 + pbp2) + sman) + sman*(mq1c**2 - 6*pap2 + 2*pbp2 + sman*(-1 + x1))*(-1 + x1 + x2)) + 
     &                  m1**2*(m3**2*(-3*mq1c**2 + 2*(pap1 + pap2 + 5*(pbp1 + pbp2)) - sman) + m4**2*(3*mq1c**2 - 2*(pap1 + pap2 + 5*(pbp1 + pbp2)) + sman) + 
     &                     sman*(-3*mq1c**2 + 2*pap2 + 10*pbp2 + sman*(-1 + x1))*(-1 + x1 + x2))) + 
     &               F1Lc*m3*mq1c*(m1**2*(2*pap1 + 2*pap2 + 10*pbp1 + 10*pbp2 + sman*(4 - 5*x1 - 5*x2)) + m2**2*(6*pap1 + 6*pap2 - 2*pbp1 - 2*pbp2 + sman*(4 - 3*x1 - 3*x2)) + 
     &                  sman*(2*(-3*pap2 + pbp1 + pbp2 - 2*sman) + (4*pap2 + 3*sman)*(x1 + x2) + pap1*(-6 + 4*x1 + 4*x2)))))) + 
     &      A1R*(2*A1Lc*m1*(F1R*mq1*(F1Lc*m3*mq1c*(-2*m1**2 + 2*m2**2 + 2*pap1 + 2*pap2 + 6*pbp1 + 6*pbp2 + sman*(-2 + x1 + x2)) + 
     &               F1Rc*(m3**2*mq1c**2 - m4**2*mq1c**2 + 4*mq1c**2*pbp1 - 4*pap2*pbp1 + 4*pap1*pbp2 - mq1c**2*sman - 2*pap1*sman - 2*pbp1*sman + 2*pap1*sman*x1 + pap2*sman*x1 + 2*pbp1*sman*x1 + pbp2*sman*x1 + 
     &                  (mq1c**2 + pap1 + pbp1)*sman*x2 + m1**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + m2**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2)))) + 
     &            F1L*(F1Rc*m3*(-2*(m1**2 - m2**2 - pap1 - pap2 - pbp1 - pbp2)*sman*(-1 + x1 + x2) + mq1c**2*(4*pbp1 + 4*pbp2 + sman*(-2 + x1 + x2))) + 
     &               F1Lc*mq1c*(2*m3**2*pap1 - 2*m4**2*pap1 + 2*m3**2*pap2 - 2*m4**2*pap2 + 6*m3**2*pbp1 - 6*m4**2*pbp1 + 4*pap2*pbp1 + 6*m3**2*pbp2 - 6*m4**2*pbp2 - 4*pap1*pbp2 - m3**2*sman + m4**2*sman - 2*pap2*sman - 
     &                  6*pbp2*sman + sman**2 + pap2*sman*x1 + 5*pbp2*sman*x1 - 2*sman**2*x1 + sman**2*x1**2 + sman*(pap1 + 2*pap2 + pbp1 + 6*pbp2 + sman*(-1 + x1))*x2 + m1**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + 
     &                  m2**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2))))) + A1Rc*(F1R*mq1*(F1Lc*m3*
     &                (4*mq1c**2*pap1 + 4*mq1c**2*pap2 + 2*mq1c**2*sman - 2*pap1*sman - 2*pap2*sman + 2*pbp1*sman + 2*pbp2*sman - 2*sman**2 + sman**2*x1 + sman**2*x2 - 
     &                  m2**2*(2*mq1c**2 - 6*pap1 - 6*pap2 + 2*pbp1 + 2*pbp2 + sman*(-2 + x1 + x2)) + m1**2*(-6*mq1c**2 + 2*pap1 + 2*pap2 + 10*pbp1 + 10*pbp2 + sman*(-2 + x1 + x2))) + 
     &               F1Rc*mq1c*(-12*pap2*pbp1 + 12*pap1*pbp2 + 2*(m3**2 - m4**2 - 3*pap1)*sman + 2*pbp1*sman - 2*sman**2 + 4*pap1*sman*x1 + 2*pap2*sman*x1 + sman**2*x1 + 
     &                  m1**2*(-2*m3**2 + 2*m4**2 + 2*pap1 + 10*pbp1 + sman*(2 - 3*x1 - 2*x2)) + 2*sman*(pap1 + sman)*x2 - m2**2*(2*m3**2 - 2*m4**2 - 6*pap1 + 2*pbp1 + sman*(-2 + x1 + 2*x2)))) + 
     &            F1L*(F1Lc*(4*m3**2*mq1c**2*pap1 - 4*m4**2*mq1c**2*pap1 + 4*m3**2*mq1c**2*pap2 - 4*m4**2*mq1c**2*pap2 + 8*mq1c**2*pap2*pbp1 - 8*mq1c**2*pap1*pbp2 + m3**2*mq1c**2*sman - m4**2*mq1c**2*sman - 
     &                  2*m3**2*pap1*sman + 2*m4**2*pap1*sman - 2*m3**2*pap2*sman + 2*m4**2*pap2*sman - 4*mq1c**2*pap2*sman + 2*m3**2*pbp1*sman - 2*m4**2*pbp1*sman - 4*pap2*pbp1*sman + 2*m3**2*pbp2*sman - 
     &                  2*m4**2*pbp2*sman + 4*pap1*pbp2*sman - m3**2*sman**2 + m4**2*sman**2 - mq1c**2*sman**2 + 2*pap2*sman**2 - 2*pbp2*sman**2 + sman**3 + 2*mq1c**2*pap2*sman*x1 + 4*pap2*pbp1*sman*x1 - 
     &                  4*pap1*pbp2*sman*x1 + mq1c**2*sman**2*x1 - 2*pap2*sman**2*x1 + 2*pbp2*sman**2*x1 - 2*sman**3*x1 + sman**3*x1**2 + 
     &                  sman*(-4*pap1*pbp2 + pap2*(4*pbp1 - 2*sman) + mq1c**2*(2*pap1 + 4*pap2 + sman) + sman*(2*pbp2 + sman*(-1 + x1)))*x2 - 
     &                  m2**2*(m3**2*(mq1c**2 - 6*pap1 - 6*pap2 + 2*(pbp1 + pbp2) - sman) + m4**2*(-mq1c**2 + 6*pap1 + 6*pap2 - 2*(pbp1 + pbp2) + sman) + sman*(mq1c**2 - 6*pap2 + 2*pbp2 + sman*(-1 + x1))*(-1 + x1 + x2)) + 
     &                  m1**2*(m3**2*(-3*mq1c**2 + 2*(pap1 + pap2 + 5*(pbp1 + pbp2)) - sman) + m4**2*(3*mq1c**2 - 2*(pap1 + pap2 + 5*(pbp1 + pbp2)) + sman) + 
     &                     sman*(-3*mq1c**2 + 2*pap2 + 10*pbp2 + sman*(-1 + x1))*(-1 + x1 + x2))) + 
     &               F1Rc*m3*mq1c*(m1**2*(2*pap1 + 2*pap2 + 10*pbp1 + 10*pbp2 + sman*(4 - 5*x1 - 5*x2)) + m2**2*(6*pap1 + 6*pap2 - 2*pbp1 - 2*pbp2 + sman*(4 - 3*x1 - 3*x2)) + 
     &                  sman*(2*(-3*pap2 + pbp1 + pbp2 - 2*sman) + (4*pap2 + 3*sman)*(x1 + x2) + pap1*(-6 + 4*x1 + 4*x2))))))))/
     &  (2.*(mq1c**2 - sman)*(-m1**2 + m2**2 - 2*(pbp1 + pbp2) + sman)*(mq1c**2 - sman*(-1 + x1 + x2))*(-mq1**2 + sman*(-1 + x1 + x2))))


        radSiSf = dreal((g0*g1*(A1R*(-2*A1Lc*m1*(F1R*mq1*(-2*F1Lc*m3*mq1c*(-m1**2 + m2**2 + m3**2 - m4**2 + 2*pbp1 - 2*pbp2 + sman*x2) + 
     &               F1Rc*(-2*m3**2*pap1 + 2*m4**2*pap1 - 2*m3**2*pbp1 + 2*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 8*pbp1*pbp2 - m3**2*sman + m4**2*sman - 4*pbp1*sman + sman**2 + 
     &                  m1**2*(m3**2 - m4**2 + 2*pap1 + 2*pbp1 + sman*(-1 + x2)) - sman*(2*(pap1 + pbp1) + sman)*x2 + m2**2*(-m3**2 + m4**2 - 2*pap1 - 2*pbp1 + sman - sman*x2))) + 
     &            F1L*(-(F1Rc*m3*(2*m3**2*pap1 - 2*m4**2*pap1 + 2*m3**2*pap2 - 2*m4**2*pap2 + 2*m3**2*pbp1 - 2*m4**2*pbp1 + 2*m3**2*pbp2 - 2*m4**2*pbp2 - 8*pap1*pbp2 - 8*pap2*pbp2 - 8*pbp1*pbp2 - 8*pbp2**2 + 4*pbp1*sman + 
     &                    4*pbp2*sman - 2*sman**2 + sman**2*x1 + sman*(2*(pap1 + pap2 + pbp1 + pbp2) + sman)*x2 - m1**2*(2*pap1 + 2*pap2 + 2*pbp1 + 2*pbp2 + sman*(-2 + x1 + x2)) + 
     &                    m2**2*(2*pap1 + 2*pap2 + 2*pbp1 + 2*pbp2 + sman*(-2 + x1 + x2)))) + 
     &               F1Lc*mq1c*(-m3**4 + 2*m3**2*m4**2 - m4**4 - 4*m3**2*pbp1 + 4*m4**2*pbp1 + m3**2*sman - m4**2*sman - m3**2*sman*x1 + m4**2*sman*x1 + 2*pbp2*sman*x1 + 
     &                  sman*(-2*m3**2 + 2*m4**2 - 2*pbp1 + sman - sman*x1)*x2 - sman**2*x2**2 + m2**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + m1**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2))))) + 
     &         A1Rc*(2*F1R*mq1*(F1Rc*mq1c*(2*m3**2*pap1 - 2*m4**2*pap1 - 4*pap2*pbp1 - 4*pap1*pbp2 + m3**2*sman - m4**2*sman + 2*pbp1*sman - sman**2 + sman*(2*pap1 + sman)*x2 + 
     &                  m2**2*(-m3**2 + m4**2 + 2*pap1 - 2*pbp1 + sman - sman*x2) + m1**2*(-m3**2 + m4**2 - 2*pap1 + 2*pbp1 + sman - sman*x2)) + 
     &               F1Lc*m3*(m1**2*(m3**2 - m4**2 + 4*pbp1 + sman*(-3 + x2)) + sman*(m3**2 - m4**2 - 4*pbp2 + sman + sman*x2) - m2**2*(m3**2 - m4**2 - 4*pap1 - 4*pap2 - 4*pbp2 + sman + sman*x2))) + 
     &            F1L*(2*F1Rc*m3*mq1c*(2*m3**2*pap1 - 2*m4**2*pap1 + 2*m3**2*pap2 - 2*m4**2*pap2 - 8*pap1*pbp2 - 8*pap2*pbp2 + 2*pbp1*sman + 2*pbp2*sman - 2*sman**2 + sman**2*x1 + sman*(2*(pap1 + pap2) + sman)*x2 + 
     &                  m2**2*(2*pap1 + 2*pap2 - 2*pbp1 - 2*pbp2 - sman*(-2 + x1 + x2)) - m1**2*(2*pap1 + 2*pap2 - 2*pbp1 - 2*pbp2 + sman*(-2 + x1 + x2))) + 
     &               F1Lc*((m3**2 - m4**2 - 4*pbp2 + sman + sman*x2)*(4*pap2*pbp1 - 4*pap1*pbp2 + sman*(m3**2 - m4**2 + sman*(-1 + x1 + x2))) + 
     &                  m1**2*(m3**4 + m4**4 - 8*pap2*pbp1 + 8*pap1*pbp2 - 4*pbp2*sman + 3*sman**2 - 3*sman**2*x1 + sman*(4*(pbp1 + pbp2 - sman) + sman*x1)*x2 + sman**2*x2**2 - 
     &                     m4**2*(8*pbp1 + 4*pbp2 + sman*(-4 + x1 + 2*x2)) + m3**2*(-2*m4**2 + 8*pbp1 + 4*pbp2 + sman*(-4 + x1 + 2*x2))) - 
     &                  m2**2*(m3**4 + m4**4 - 8*pap2*pbp1 + 8*pap1*pbp2 + 8*pap2*sman + 4*pbp2*sman - sman**2 - 4*pap2*sman*x1 - 4*pbp2*sman*x1 + sman**2*x1 + sman*(-4*(pap1 + 2*pap2 + pbp2) + sman*x1)*x2 + sman**2*x2**2 + 
     &                     m3**2*(-2*m4**2 - 8*pap1 - 8*pap2 - 4*pbp2 + sman*x1 + 2*sman*x2) + m4**2*(8*pap1 + 8*pap2 + 4*pbp2 - sman*(x1 + 2*x2))))))) + 
     &      A1L*(-2*A1Rc*m1*(F1L*mq1*(-2*F1Rc*m3*mq1c*(-m1**2 + m2**2 + m3**2 - m4**2 + 2*pbp1 - 2*pbp2 + sman*x2) + 
     &               F1Lc*(-2*m3**2*pap1 + 2*m4**2*pap1 - 2*m3**2*pbp1 + 2*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 8*pbp1*pbp2 - m3**2*sman + m4**2*sman - 4*pbp1*sman + sman**2 + 
     &                  m1**2*(m3**2 - m4**2 + 2*pap1 + 2*pbp1 + sman*(-1 + x2)) - sman*(2*(pap1 + pbp1) + sman)*x2 + m2**2*(-m3**2 + m4**2 - 2*pap1 - 2*pbp1 + sman - sman*x2))) + 
     &            F1R*(-(F1Lc*m3*(2*m3**2*pap1 - 2*m4**2*pap1 + 2*m3**2*pap2 - 2*m4**2*pap2 + 2*m3**2*pbp1 - 2*m4**2*pbp1 + 2*m3**2*pbp2 - 2*m4**2*pbp2 - 8*pap1*pbp2 - 8*pap2*pbp2 - 8*pbp1*pbp2 - 8*pbp2**2 + 4*pbp1*sman + 
     &                    4*pbp2*sman - 2*sman**2 + sman**2*x1 + sman*(2*(pap1 + pap2 + pbp1 + pbp2) + sman)*x2 - m1**2*(2*pap1 + 2*pap2 + 2*pbp1 + 2*pbp2 + sman*(-2 + x1 + x2)) + 
     &                    m2**2*(2*pap1 + 2*pap2 + 2*pbp1 + 2*pbp2 + sman*(-2 + x1 + x2)))) + 
     &               F1Rc*mq1c*(-m3**4 + 2*m3**2*m4**2 - m4**4 - 4*m3**2*pbp1 + 4*m4**2*pbp1 + m3**2*sman - m4**2*sman - m3**2*sman*x1 + m4**2*sman*x1 + 2*pbp2*sman*x1 + 
     &                  sman*(-2*m3**2 + 2*m4**2 - 2*pbp1 + sman - sman*x1)*x2 - sman**2*x2**2 + m2**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + m1**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2))))) + 
     &         A1Lc*(2*F1L*mq1*(F1Lc*mq1c*(2*m3**2*pap1 - 2*m4**2*pap1 - 4*pap2*pbp1 - 4*pap1*pbp2 + m3**2*sman - m4**2*sman + 2*pbp1*sman - sman**2 + sman*(2*pap1 + sman)*x2 + 
     &                  m2**2*(-m3**2 + m4**2 + 2*pap1 - 2*pbp1 + sman - sman*x2) + m1**2*(-m3**2 + m4**2 - 2*pap1 + 2*pbp1 + sman - sman*x2)) + 
     &               F1Rc*m3*(m1**2*(m3**2 - m4**2 + 4*pbp1 + sman*(-3 + x2)) + sman*(m3**2 - m4**2 - 4*pbp2 + sman + sman*x2) - m2**2*(m3**2 - m4**2 - 4*pap1 - 4*pap2 - 4*pbp2 + sman + sman*x2))) + 
     &            F1R*(2*F1Lc*m3*mq1c*(2*m3**2*pap1 - 2*m4**2*pap1 + 2*m3**2*pap2 - 2*m4**2*pap2 - 8*pap1*pbp2 - 8*pap2*pbp2 + 2*pbp1*sman + 2*pbp2*sman - 2*sman**2 + sman**2*x1 + sman*(2*(pap1 + pap2) + sman)*x2 + 
     &                  m2**2*(2*pap1 + 2*pap2 - 2*pbp1 - 2*pbp2 - sman*(-2 + x1 + x2)) - m1**2*(2*pap1 + 2*pap2 - 2*pbp1 - 2*pbp2 + sman*(-2 + x1 + x2))) + 
     &               F1Rc*((m3**2 - m4**2 - 4*pbp2 + sman + sman*x2)*(4*pap2*pbp1 - 4*pap1*pbp2 + sman*(m3**2 - m4**2 + sman*(-1 + x1 + x2))) + 
     &                  m1**2*(m3**4 + m4**4 - 8*pap2*pbp1 + 8*pap1*pbp2 - 4*pbp2*sman + 3*sman**2 - 3*sman**2*x1 + sman*(4*(pbp1 + pbp2 - sman) + sman*x1)*x2 + sman**2*x2**2 - 
     &                     m4**2*(8*pbp1 + 4*pbp2 + sman*(-4 + x1 + 2*x2)) + m3**2*(-2*m4**2 + 8*pbp1 + 4*pbp2 + sman*(-4 + x1 + 2*x2))) - 
     &                  m2**2*(m3**4 + m4**4 - 8*pap2*pbp1 + 8*pap1*pbp2 + 8*pap2*sman + 4*pbp2*sman - sman**2 - 4*pap2*sman*x1 - 4*pbp2*sman*x1 + sman**2*x1 + sman*(-4*(pap1 + 2*pap2 + pbp2) + sman*x1)*x2 + sman**2*x2**2 + 
     &                     m3**2*(-2*m4**2 - 8*pap1 - 8*pap2 - 4*pbp2 + sman*x1 + 2*sman*x2) + m4**2*(8*pap1 + 8*pap2 + 4*pbp2 - sman*(x1 + 2*x2)))))))))/
     &  (2.*(mq1c**2 - sman)*(-m1**2 + m2**2 - 2*(pbp1 + pbp2) + sman)*(m3**2 - m4**2 + sman*(-1 + x2))*(-mq1**2 + sman*(-1 + x1 + x2))))


        radSpSf = dreal((g1**2*(A1R*(2*A1Lc*m1*(-(F1Rc*(F1L*m3*(2*mq1**2*(pap1 + 2*pap2 + pbp1 + 2*pbp2 - 2*sman) + sman*(-m3**2 + m4**2 + sman*(3 - 2*x1 - 2*x2))) + 
     &                 F1R*mq1*(m3**2*(2*(pap1 + 2*pap2 + pbp1 + 2*pbp2) - 5*sman) + m4**2*(-2*(pap2 + pbp2) + 3*sman) + sman*(sman*(3 - 2*x1 - 3*x2) + 2*pap2*(-1 + x1 + x2) + 2*pbp2*(-1 + x1 + x2))))) + 
     &            F1Lc*mq1c*(F1R*m3*mq1*(m3**2 - m4**2 - 2*pap1 - 4*pap2 - 2*pbp1 - 4*pbp2 + sman + 2*sman*x1 + 2*sman*x2) + 
     &               F1L*(-(m4**2*(mq1**2 - 2*(pap2 + pbp2 - sman))) + m3**2*(3*mq1**2 - 2*(pap1 + 2*pap2 + pbp1 + 2*pbp2 - sman)) + sman*(mq1**2*(-1 + x2) - 2*(pap2 + pbp2 - sman)*(-1 + x1 + x2))))) + 
     &         A1Rc*(F1Lc*(F1L*(-4*mq1**2*pap2*pbp1 + 4*mq1**2*pap1*pbp2 + 3*m3**2*mq1**2*sman - m4**2*mq1**2*sman - 2*mq1**2*pap1*sman - 4*m3**2*pbp1*sman + 2*mq1**2*pbp1*sman - 8*m3**2*pbp2*sman + 4*m4**2*pbp2*sman + 
     &                  2*m3**2*sman**2 - 2*m4**2*sman**2 - mq1**2*sman**2 + 4*pbp2*sman**2 - 2*sman**3 - 4*pbp2*sman**2*x1 + 2*sman**3*x1 + sman**2*(mq1**2 - 4*pbp2 + 2*sman)*x2 + 
     &                  m2**2*(-2*mq1**2*pap1 - 2*mq1**2*pbp1 + m3**2*(-3*mq1**2 + 4*(pap1 + 2*pap2 + pbp1 + 2*pbp2) - 2*sman) + mq1**2*sman - 4*pap2*sman - 4*pbp2*sman + 2*sman**2 + 
     &                     m4**2*(mq1**2 - 4*(pap2 + pbp2) + 2*sman) + 4*pap2*sman*x1 + 4*pbp2*sman*x1 - 2*sman**2*x1 - sman*(mq1**2 - 4*(pap2 + pbp2) + 2*sman)*x2) + 
     &                  m1**2*(2*mq1**2*pap1 + 2*mq1**2*pbp1 - mq1**2*sman + 4*pap2*sman + 4*pbp2*sman - 2*sman**2 - m4**2*(mq1**2 - 4*(pap2 + pbp2) + 2*sman) + 
     &                     m3**2*(3*mq1**2 - 4*(pap1 + 2*pap2 + pbp1 + 2*pbp2) + 2*sman) - 4*pap2*sman*x1 - 4*pbp2*sman*x1 + 2*sman**2*x1 + sman*(mq1**2 - 4*(pap2 + pbp2) + 2*sman)*x2)) + 
     &               F1R*m3*mq1*(-4*pap2*pbp1 + 4*pap1*pbp2 + (m3 - m4)*(m3 + m4)*sman - 2*pap1*sman - 2*pap2*sman - 2*pbp1*sman - 6*pbp2*sman + sman**2 + 2*sman**2*x1 + 2*sman**2*x2 + 
     &                  m1**2*(m3**2 - m4**2 - 2*pap1 - 6*pap2 - 2*pbp1 - 6*pbp2 + sman + 2*sman*x1 + 2*sman*x2) + m2**2*(-m3**2 + m4**2 + 2*pap1 + 6*pap2 + 2*pbp1 + 6*pbp2 - sman*(1 + 2*x1 + 2*x2)))) + 
     &            F1Rc*mq1c*(F1L*m3*(-4*mq1**2*pap1 - 8*mq1**2*pap2 + 4*pap2*pbp1 - 4*pap1*pbp2 + (m3 - m4)*(m3 + m4)*sman + 4*mq1**2*sman + 2*pap1*sman + 2*pap2*sman - 2*pbp1*sman - 2*pbp2*sman - 3*sman**2 + 2*sman**2*x1 + 
     &                  m2**2*(-m3**2 + m4**2 - 4*mq1**2 + 2*pap1 + 2*pap2 + 2*pbp1 + 2*pbp2 + sman*(3 - 2*x1 - 2*x2)) + 2*sman**2*x2 + 
     &                  m1**2*(m3**2 - m4**2 + 4*mq1**2 - 2*pap1 - 2*pap2 - 2*pbp1 - 2*pbp2 + sman*(-3 + 2*x1 + 2*x2))) + 
     &               F1R*mq1*(-4*m3**2*pap1 - 8*m3**2*pap2 + 4*m4**2*pap2 + 4*pap2*pbp1 - 4*pap1*pbp2 + 5*m3**2*sman - 3*m4**2*sman + 2*pap1*sman + 4*pap2*sman - 2*pbp1*sman - 3*sman**2 - 4*pap2*sman*x1 + 2*sman**2*x1 + 
     &                  m2**2*(-5*m3**2 + 3*m4**2 + 2*pap1 + 2*pbp1 + sman*(3 - 2*x1 - 3*x2)) - 4*pap2*sman*x2 + 3*sman**2*x2 + m1**2*(5*m3**2 - 3*m4**2 - 2*pap1 - 2*pbp1 + sman*(-3 + 2*x1 + 3*x2)))))) + 
     &      A1L*(-2*A1Rc*m1*(F1L*mq1*(F1Rc*m3*mq1c*(-m3**2 + m4**2 + 2*pap1 + 4*pap2 + 2*pbp1 + 4*pbp2 - sman*(1 + 2*x1 + 2*x2)) + 
     &               F1Lc*(m3**2*(2*(pap1 + 2*pap2 + pbp1 + 2*pbp2) - 5*sman) + m4**2*(-2*(pap2 + pbp2) + 3*sman) + sman*(sman*(3 - 2*x1 - 3*x2) + 2*pap2*(-1 + x1 + x2) + 2*pbp2*(-1 + x1 + x2)))) + 
     &            F1R*(F1Lc*m3*(2*mq1**2*(pap1 + 2*pap2 + pbp1 + 2*pbp2 - 2*sman) + sman*(-m3**2 + m4**2 + sman*(3 - 2*x1 - 2*x2))) + 
     &               F1Rc*mq1c*(m4**2*(mq1**2 - 2*(pap2 + pbp2 - sman)) + m3**2*(-3*mq1**2 + 2*(pap1 + 2*pap2 + pbp1 + 2*pbp2 - sman)) + sman*(-(mq1**2*(-1 + x2)) + 2*(pap2 + pbp2 - sman)*(-1 + x1 + x2))))) + 
     &         A1Lc*(F1R*(F1Rc*(-4*mq1**2*pap2*pbp1 + 4*mq1**2*pap1*pbp2 + 3*m3**2*mq1**2*sman - m4**2*mq1**2*sman - 2*mq1**2*pap1*sman - 4*m3**2*pbp1*sman + 2*mq1**2*pbp1*sman - 8*m3**2*pbp2*sman + 4*m4**2*pbp2*sman + 
     &                  2*m3**2*sman**2 - 2*m4**2*sman**2 - mq1**2*sman**2 + 4*pbp2*sman**2 - 2*sman**3 - 4*pbp2*sman**2*x1 + 2*sman**3*x1 + sman**2*(mq1**2 - 4*pbp2 + 2*sman)*x2 + 
     &                  m2**2*(-2*mq1**2*pap1 - 2*mq1**2*pbp1 + m3**2*(-3*mq1**2 + 4*(pap1 + 2*pap2 + pbp1 + 2*pbp2) - 2*sman) + mq1**2*sman - 4*pap2*sman - 4*pbp2*sman + 2*sman**2 + 
     &                     m4**2*(mq1**2 - 4*(pap2 + pbp2) + 2*sman) + 4*pap2*sman*x1 + 4*pbp2*sman*x1 - 2*sman**2*x1 - sman*(mq1**2 - 4*(pap2 + pbp2) + 2*sman)*x2) + 
     &                  m1**2*(2*mq1**2*pap1 + 2*mq1**2*pbp1 - mq1**2*sman + 4*pap2*sman + 4*pbp2*sman - 2*sman**2 - m4**2*(mq1**2 - 4*(pap2 + pbp2) + 2*sman) + 
     &                     m3**2*(3*mq1**2 - 4*(pap1 + 2*pap2 + pbp1 + 2*pbp2) + 2*sman) - 4*pap2*sman*x1 - 4*pbp2*sman*x1 + 2*sman**2*x1 + sman*(mq1**2 - 4*(pap2 + pbp2) + 2*sman)*x2)) + 
     &               F1Lc*m3*mq1c*(-4*mq1**2*pap1 - 8*mq1**2*pap2 + 4*pap2*pbp1 - 4*pap1*pbp2 + (m3 - m4)*(m3 + m4)*sman + 4*mq1**2*sman + 2*pap1*sman + 2*pap2*sman - 2*pbp1*sman - 2*pbp2*sman - 3*sman**2 + 2*sman**2*x1 + 
     &                  m2**2*(-m3**2 + m4**2 - 4*mq1**2 + 2*pap1 + 2*pap2 + 2*pbp1 + 2*pbp2 + sman*(3 - 2*x1 - 2*x2)) + 2*sman**2*x2 + 
     &                  m1**2*(m3**2 - m4**2 + 4*mq1**2 - 2*pap1 - 2*pap2 - 2*pbp1 - 2*pbp2 + sman*(-3 + 2*x1 + 2*x2)))) + 
     &            F1L*mq1*(F1Rc*m3*(-4*pap2*pbp1 + 4*pap1*pbp2 + (m3 - m4)*(m3 + m4)*sman - 2*pap1*sman - 2*pap2*sman - 2*pbp1*sman - 6*pbp2*sman + sman**2 + 2*sman**2*x1 + 2*sman**2*x2 + 
     &                  m1**2*(m3**2 - m4**2 - 2*pap1 - 6*pap2 - 2*pbp1 - 6*pbp2 + sman + 2*sman*x1 + 2*sman*x2) + m2**2*(-m3**2 + m4**2 + 2*pap1 + 6*pap2 + 2*pbp1 + 6*pbp2 - sman*(1 + 2*x1 + 2*x2))) + 
     &               F1Lc*mq1c*(-4*m3**2*pap1 - 8*m3**2*pap2 + 4*m4**2*pap2 + 4*pap2*pbp1 - 4*pap1*pbp2 + 5*m3**2*sman - 3*m4**2*sman + 2*pap1*sman + 4*pap2*sman - 2*pbp1*sman - 3*sman**2 - 4*pap2*sman*x1 + 2*sman**2*x1 + 
     &                  m2**2*(-5*m3**2 + 3*m4**2 + 2*pap1 + 2*pbp1 + sman*(3 - 2*x1 - 3*x2)) - 4*pap2*sman*x2 + 3*sman**2*x2 + m1**2*(5*m3**2 - 3*m4**2 - 2*pap1 - 2*pbp1 + sman*(-3 + 2*x1 + 3*x2))))))))/
     &  ((mq1**2 - sman)*(-mq1c**2 + sman)*(m3**2 - m4**2 + sman*(-1 + x2))*(mq1**2 - sman*(-1 + x1 + x2))))

        NeuQ2qx_M2RadiationHSS = CF*NC*(radSiSi + radSpSp + radSfSf + 2d0*radSiSp + 2d0*radSiSf + 2d0*radSpSf)
!            print*,"NeuQ2qx_M2RadiationHSS",NeuQ2qx_M2RadiationHSS

	end

! =================================================================
	double precision function NeuQ2qx_M2RadiationHTT()

        implicit none


#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

	double precision radTiTi, radTiTp, radTiTf, radTpTp, radTpTf, radTfTf
        double complex A2L, A2Lc, A2R, A2Rc, G2, G2c
        double complex g0, g1
        double precision mq2, mq2c
	integer sumT, sumTc
	integer sumTend, sumTcend

        NeuQ2qx_M2RadiationHTT = 0d0

        radTiTi = 0d0        
        radTpTp = 0d0
        radTfTf = 0d0
        radTiTp = 0d0
        radTiTf = 0d0
        radTpTf = 0d0

	sumT = 0
	sumTc = 0

	if ((xtype.eq.5).or.(xtype.eq.6)) then

	sumTend = 1
	sumTcend = 1

	elseif (xtype.eq.7) then

	sumTend = 0
	sumTcend = 0

	elseif (xtype.eq.8) then

	sumTend = 1
	sumTcend = 1

        endif

	do sumT = 0,sumTend
	   do sumTc = 0,sumTcend

        A2L = AL(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2R = AR(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))

	G2 = G(xtype,sfmeigen,ptmeigen+sumT,sftype,sfgen)
	G2c = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

        g0 = -gs
	g1 = -gs
        
        mq2 = MSf(ptmeigen+sumT,pttype,ptgen)
        mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)


        radTiTi = dreal((-4*G2*g0**2*G2c*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*(m1**2 + m2**2 + 2*(pbp1 + pbp2) - sman))/
     &  ((m1**2 + m3**2 - mq2**2 - 2*pap1)*(m1**2 + m3**2 - mq2c**2 - 2*pap1)*(-m1**2 + m2**2 - 2*(pbp1 + pbp2) + sman)**2))


        radTpTp = dreal((-2*G2*g0**2*G2c*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*(2*m1**2 + 2*m2**2 + 2*m3**2 + 2*m4**2 - 2*pap1 + 2*pap2 + 2*pbp1 - 2*pbp2 - sman*(x1 + x2)))/
     &  ((m1**2 + m3**2 - mq2**2 - 2*pap1)*(m1**2 + m3**2 - mq2c**2 - 2*pap1)*(-m2**2 - m4**2 + mq2**2 + 2*pbp2)*(-m2**2 - m4**2 + mq2c**2 + 2*pbp2)))


        radTfTf = dreal((-2*G2*g1**2*G2c*(A2L*(2*A2Rc*m1*m3*(m3**2 + m4**2 + sman - sman*x2) + A2Lc*(-2*m3**2*pap1 + 2*m4**2*pap1 - 6*m3**2*pap2 + 2*m4**2*pap2 + 3*m3**2*sman - m4**2*sman + 2*pap1*sman + 2*pap2*sman - sman**2 + 
     &            m1**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) + sman*(-2*(pap1 + pap2) + sman)*x2 + m2**2*(-3*m3**2 + m4**2 + sman - sman*x2))) + 
     &      A2R*(2*A2Lc*m1*m3*(m3**2 + m4**2 + sman - sman*x2) + A2Rc*(-2*m3**2*pap1 + 2*m4**2*pap1 - 6*m3**2*pap2 + 2*m4**2*pap2 + 3*m3**2*sman - m4**2*sman + 2*pap1*sman + 2*pap2*sman - sman**2 + 
     &            m1**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) + sman*(-2*(pap1 + pap2) + sman)*x2 + m2**2*(-3*m3**2 + m4**2 + sman - sman*x2)))))/
     &  ((m2**2 + m4**2 - mq2**2 - 2*pbp2)*(m2**2 + m4**2 - mq2c**2 - 2*pbp2)*(-m3**2 + m4**2 + sman - sman*x2)**2))


        radTiTp = dreal((G2*g0**2*G2c*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*(4*m1**2 + 4*m2**2 + 2*m3**2 - 2*m4**2 - 2*pap1 - 2*pap2 + 6*pbp1 - 2*pbp2 + sman*(-2 - x1 + x2)))/
     &  ((m1**2 + m3**2 - mq2**2 - 2*pap1)*(m1**2 + m3**2 - mq2c**2 - 2*pap1)*(-m2**2 - m4**2 + mq2c**2 + 2*pbp2)*(m1**2 - m2**2 + 2*(pbp1 + pbp2) - sman)))


        radTiTf = dreal( (G2*g0*g1*G2c*(A2Lc*(2*A2R*m1*m3*(-m1**2 + m2**2 + m3**2 - m4**2 + 2*pbp1 - 2*pbp2 + sman*x2) + 
     &         A2L*(2*m3**2*pap1 - 2*m4**2*pap1 - 4*pap2*pbp1 - 4*pap1*pbp2 + m3**2*sman - m4**2*sman + 2*pbp1*sman - sman**2 + sman*(2*pap1 + sman)*x2 + m2**2*(-m3**2 + m4**2 + 2*pap1 - 2*pbp1 + sman - sman*x2) + 
     &            m1**2*(-m3**2 + m4**2 - 2*pap1 + 2*pbp1 + sman - sman*x2))) + A2Rc*(2*A2L*m1*m3*(-m1**2 + m2**2 + m3**2 - m4**2 + 2*pbp1 - 2*pbp2 + sman*x2) + 
     &         A2R*(2*m3**2*pap1 - 2*m4**2*pap1 - 4*pap2*pbp1 - 4*pap1*pbp2 + m3**2*sman - m4**2*sman + 2*pbp1*sman - sman**2 + sman*(2*pap1 + sman)*x2 + m2**2*(-m3**2 + m4**2 + 2*pap1 - 2*pbp1 + sman - sman*x2) + 
     &            m1**2*(-m3**2 + m4**2 - 2*pap1 + 2*pbp1 + sman - sman*x2)))))/((m1**2 + m3**2 - mq2**2 - 2*pap1)*(m2**2 + m4**2 - mq2c**2 - 2*pbp2)*(m1**2 - m2**2 + 2*(pbp1 + pbp2) - sman)*(m3**2 - m4**2 + sman*(-1 + x2))))


        radTpTf = dreal((G2*g0*g1*G2c*(A2Rc*(A2R*(-4*m4**2*pap1 - 8*pap1*pap2 + 4*pap2*pbp1 + 4*pap1*pbp2 + 2*pap1*sman - 2*pbp1*sman - 4*m3**2*(pap1 - 2*pap2 + sman) - 2*m1**2*(2*m4**2 - 3*pap1 + pbp1 + 2*sman) + m1**2*sman*x1 - 
     &            2*pap2*sman*x1 + sman**2*x1 + m2**2*(4*m3**2 - 6*pap1 + 2*pbp1 - sman*x1) + 2*(2*m1**2 + pap1)*sman*x2) + 
     &         2*A2L*m1*m3*(2*m1**2 - 2*m2**2 - 4*m3**2 - 4*m4**2 + 2*pap1 - 2*pap2 - 2*pbp1 + 2*pbp2 + sman*(-2 + x1 + 3*x2))) + 
     &      A2Lc*(A2L*(-4*m4**2*pap1 - 8*pap1*pap2 + 4*pap2*pbp1 + 4*pap1*pbp2 + 2*pap1*sman - 2*pbp1*sman - 4*m3**2*(pap1 - 2*pap2 + sman) - 2*m1**2*(2*m4**2 - 3*pap1 + pbp1 + 2*sman) + m1**2*sman*x1 - 2*pap2*sman*x1 + 
     &            sman**2*x1 + m2**2*(4*m3**2 - 6*pap1 + 2*pbp1 - sman*x1) + 2*(2*m1**2 + pap1)*sman*x2) + 2*A2R*m1*m3*(2*m1**2 - 2*m2**2 - 4*m3**2 - 4*m4**2 + 2*pap1 - 2*pap2 - 2*pbp1 + 2*pbp2 + sman*(-2 + x1 + 3*x2)))))/
     &  (2.*(m1**2 + m3**2 - mq2**2 - 2*pap1)*(m2**2 + m4**2 - mq2c**2 - 2*pbp2)*(-m2**2 - m4**2 + mq2**2 + 2*pbp2)*(m3**2 - m4**2 + sman*(-1 + x2))))


        NeuQ2qx_M2RadiationHTT = NeuQ2qx_M2RadiationHTT + CF*NC*(radTiTi + radTpTp + radTfTf + 2d0*radTiTp + 2d0*radTiTf + 2d0*radTpTf)
!            print*,"NeuQ2qx_M2RadiationHTT",NeuQ2qx_M2RadiationHTT

           enddo
        enddo


	end

! =================================================================
	double precision function NeuQ2qx_M2RadiationHUU()

        implicit none


#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double precision radUiUi, radUiUf, radUfUf
        double complex A3L, A3Lc, A3R, A3Rc, H3L, H3Lc, H3R, H3Rc
        double complex g0, g1
        double precision mq3, mq3c
        integer sumU, sumUc, sumUend, sumUcend

        NeuQ2qx_M2RadiationHUU = 0d0

        radUiUi = 0d0
        radUfUf = 0d0
        radUiUf = 0d0

        sumU = 0
        sumUc = 0

        if ((xtype.eq.5).or.(xtype.eq.6)) then
        
           sumUend = 3
           sumUcend = 3
           
        elseif (xtype.eq.7) then
        
           sumUend = 3
           sumUcend = 3
           
        elseif (xtype.eq.8) then
               
           sumUend = 1
           sumUcend = 1
           
         endif

        do sumU = 0, sumUend
           do sumUc = 0, sumUcend

           A3L = AL(pumeigen+sumU,sfmeigen,sftype,sfgen)
           A3R = AR(pumeigen+sumU,sfmeigen,sftype,sfgen)
           A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
           A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

           H3L = HL(chimeigen,pumeigen+sumU,xtype)
           H3R = HR(chimeigen,pumeigen+sumU,xtype)
           H3Lc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
           H3Rc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

           g0 = -gs
           g1 = -gs

	   mq3 =  MNeuCha(pumeigen+sumU)
	   mq3c = MNeuCha(pumeigen+sumUc)

           radUiUi = dreal((-2*g0**2*(m1**2 + m2**2 + 2*(pbp1 + pbp2) - sman)*(A3R*(2*A3Lc*m3*(H3R*mq3*(H3Lc*m1*mq3c + H3Rc*(m1**2 - pap2)) + H3L*(H3Rc*m1*(m1**2 + m4**2 - 2*pap2) + H3Lc*mq3c*(m1**2 - pap2))) + 
     &         A3Rc*(H3R*mq3*(2*H3Rc*mq3c*pap1 + H3Lc*m1*(m3**2 + m4**2 + 2*pap1 - sman*(-1 + x1 + x2))) + 
     &            H3L*(H3Rc*m1*mq3c*(m3**2 + m4**2 + 2*pap1 - sman*(-1 + x1 + x2)) + 2*H3Lc*(-(m3**2*pap2) - m4**2*(pap1 + pap2) + pap2*sman*(-1 + x1 + x2) + m1**2*(m3**2 + m4**2 + pap1 - sman*(-1 + x1 + x2)))))) + 
     &      A3L*(2*A3Rc*m3*(H3Rc*mq3c*(H3L*m1*mq3 + H3R*(m1**2 - pap2)) + H3Lc*(H3R*m1*(m1**2 + m4**2 - 2*pap2) + H3L*mq3*(m1**2 - pap2))) + 
     &         A3Lc*(H3L*mq3*(2*H3Lc*mq3c*pap1 + H3Rc*m1*(m3**2 + m4**2 + 2*pap1 - sman*(-1 + x1 + x2))) + 
     &            H3R*(H3Lc*m1*mq3c*(m3**2 + m4**2 + 2*pap1 - sman*(-1 + x1 + x2)) + 2*H3Rc*(-(m3**2*pap2) - m4**2*(pap1 + pap2) + pap2*sman*(-1 + x1 + x2) + m1**2*(m3**2 + m4**2 + pap1 - sman*(-1 + x1 + x2))))))))/
     &  ((m1**2 + m4**2 - mq3**2 - 2*pap2)*(m1**2 + m4**2 - mq3c**2 - 2*pap2)*(-m1**2 + m2**2 - 2*(pbp1 + pbp2) + sman)**2))


           radUfUf = dreal((2*g1**2*(A3R*(-2*A3Lc*m3*(H3R*mq3*(H3Lc*m1*mq3c + H3Rc*(m1**2 - pap2)) + H3L*(H3Rc*m1*(m1**2 + m4**2 - 2*pap2) + H3Lc*mq3c*(m1**2 - pap2)))*(m3**2 + m4**2 + sman - sman*x2) + 
     &         A3Rc*(H3R*mq3*(H3Rc*mq3c*(2*m3**2*pap1 - 2*m4**2*pap1 + 6*m3**2*pap2 - 2*m4**2*pap2 - 3*m3**2*sman + m4**2*sman - 2*pap1*sman - 2*pap2*sman + sman**2 + m2**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) + 
     &                  (2*(pap1 + pap2) - sman)*sman*x2 + m1**2*(-3*m3**2 + m4**2 + sman - sman*x2)) + 
     &               H3Lc*m1*(m3**4 - 6*m3**2*m4**2 + m4**4 + 2*m3**2*pap1 - 2*m4**2*pap1 + 6*m3**2*pap2 - 2*m4**2*pap2 - 3*m3**2*sman + m4**2*sman - 2*pap1*sman - 2*pap2*sman - m3**2*sman*x1 + m4**2*sman*x1 + sman**2*x1 + 
     &                  m2**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) - sman*(-3*m3**2 + m4**2 - 2*(pap1 + pap2) + sman*x1)*x2 + m1**2*(-3*m3**2 + m4**2 + sman - sman*x2))) + 
     &            H3L*(H3Rc*m1*mq3c*(m3**4 - 6*m3**2*m4**2 + m4**4 + 2*m3**2*pap1 - 2*m4**2*pap1 + 6*m3**2*pap2 - 2*m4**2*pap2 - 3*m3**2*sman + m4**2*sman - 2*pap1*sman - 2*pap2*sman - m3**2*sman*x1 + m4**2*sman*x1 + 
     &                  sman**2*x1 + m2**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) - sman*(-3*m3**2 + m4**2 - 2*(pap1 + pap2) + sman*x1)*x2 + m1**2*(-3*m3**2 + m4**2 + sman - sman*x2)) + 
     &               H3Lc*(-2*m3**2*m4**2*pap1 + 2*m4**4*pap1 - 2*m3**4*pap2 + 6*m3**2*m4**2*pap2 + 3*m3**2*m4**2*sman - m4**4*sman + 2*m4**2*pap1*sman + 2*m4**2*pap2*sman - m4**2*sman**2 + 2*pap2*sman**2 + 
     &                  2*m3**2*pap2*sman*x1 - 2*m4**2*pap2*sman*x1 - 2*pap2*sman**2*x1 + sman*(-6*m3**2*pap2 + m4**2*(-2*pap1 + sman) + 2*pap2*sman*(-1 + x1))*x2 + m1**4*(-3*m3**2 + m4**2 + sman - sman*x2) + 
     &                  m2**2*m4**2*(-3*m3**2 + m4**2 + sman - sman*x2) + m1**2*(2*m3**4 + m2**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) + (m4**2 - 2*pap1 - 2*pap2 - sman + 2*sman*x1)*(m4**2 + sman - sman*x2) + 
     &                     m3**2*(-9*m4**2 + 2*pap1 + 6*pap2 + sman*(-3 - 2*x1 + 6*x2))))))) + 
     &      A3L*(-2*A3Rc*m3*(H3Rc*mq3c*(H3L*m1*mq3 + H3R*(m1**2 - pap2)) + H3Lc*(H3R*m1*(m1**2 + m4**2 - 2*pap2) + H3L*mq3*(m1**2 - pap2)))*(m3**2 + m4**2 + sman - sman*x2) + 
     &         A3Lc*(H3L*mq3*(H3Lc*mq3c*(2*m3**2*pap1 - 2*m4**2*pap1 + 6*m3**2*pap2 - 2*m4**2*pap2 - 3*m3**2*sman + m4**2*sman - 2*pap1*sman - 2*pap2*sman + sman**2 + m2**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) + 
     &                  (2*(pap1 + pap2) - sman)*sman*x2 + m1**2*(-3*m3**2 + m4**2 + sman - sman*x2)) + 
     &               H3Rc*m1*(m3**4 - 6*m3**2*m4**2 + m4**4 + 2*m3**2*pap1 - 2*m4**2*pap1 + 6*m3**2*pap2 - 2*m4**2*pap2 - 3*m3**2*sman + m4**2*sman - 2*pap1*sman - 2*pap2*sman - m3**2*sman*x1 + m4**2*sman*x1 + sman**2*x1 + 
     &                  m2**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) - sman*(-3*m3**2 + m4**2 - 2*(pap1 + pap2) + sman*x1)*x2 + m1**2*(-3*m3**2 + m4**2 + sman - sman*x2))) + 
     &            H3R*(H3Lc*m1*mq3c*(m3**4 - 6*m3**2*m4**2 + m4**4 + 2*m3**2*pap1 - 2*m4**2*pap1 + 6*m3**2*pap2 - 2*m4**2*pap2 - 3*m3**2*sman + m4**2*sman - 2*pap1*sman - 2*pap2*sman - m3**2*sman*x1 + m4**2*sman*x1 + 
     &                  sman**2*x1 + m2**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) - sman*(-3*m3**2 + m4**2 - 2*(pap1 + pap2) + sman*x1)*x2 + m1**2*(-3*m3**2 + m4**2 + sman - sman*x2)) + 
     &               H3Rc*(-2*m3**2*m4**2*pap1 + 2*m4**4*pap1 - 2*m3**4*pap2 + 6*m3**2*m4**2*pap2 + 3*m3**2*m4**2*sman - m4**4*sman + 2*m4**2*pap1*sman + 2*m4**2*pap2*sman - m4**2*sman**2 + 2*pap2*sman**2 + 
     &                  2*m3**2*pap2*sman*x1 - 2*m4**2*pap2*sman*x1 - 2*pap2*sman**2*x1 + sman*(-6*m3**2*pap2 + m4**2*(-2*pap1 + sman) + 2*pap2*sman*(-1 + x1))*x2 + m1**4*(-3*m3**2 + m4**2 + sman - sman*x2) + 
     &                  m2**2*m4**2*(-3*m3**2 + m4**2 + sman - sman*x2) + m1**2*(2*m3**4 + m2**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) + (m4**2 - 2*pap1 - 2*pap2 - sman + 2*sman*x1)*(m4**2 + sman - sman*x2) + 
     &                     m3**2*(-9*m4**2 + 2*pap1 + 6*pap2 + sman*(-3 - 2*x1 + 6*x2)))))))))/((m1**2 + m4**2 - mq3**2 - 2*pap2)*(m1**2 + m4**2 - mq3c**2 - 2*pap2)*(-m3**2 + m4**2 + sman - sman*x2)**2))


           radUiUf = dreal((g0*g1*(A3R*(-(A3Lc*m3*(H3R*mq3*(2*H3Lc*m1*mq3c*(-m1**2 + m2**2 + m3**2 - m4**2 + 2*pbp1 - 2*pbp2 + sman*x2) - 
     &                 H3Rc*(2*m1**4 + 2*m3**2*pap2 - 2*m4**2*pap2 + 4*pap2*pbp1 + 4*pap1*pbp2 + m3**2*sman - m4**2*sman - 2*pbp2*sman + sman**2 + m2**2*(-m3**2 + m4**2 + 2*pap2 + 2*pbp2 + sman*(-1 + x1)) - sman**2*x1 + 
     &                    2*pap2*sman*x2 - m1**2*(2*m2**2 + 3*m3**2 - 3*m4**2 + 2*pap2 + 4*pbp1 - 2*pbp2 + sman - sman*x1 + 2*sman*x2))) + 
     &              H3L*(-2*H3Rc*m1*(m1**2 + m4**2 - 2*pap2)*(m1**2 - m2**2 - m3**2 + m4**2 - 2*pbp1 + 2*pbp2 - sman*x2) + 
     &                 H3Lc*mq3c*(-2*m1**4 - 2*m3**2*pap2 + 2*m4**2*pap2 - 4*pap2*pbp1 + 4*pap1*pbp2 + 8*pap2*pbp2 + m3**2*sman - m4**2*sman - 2*pbp2*sman + sman**2 + 
     &                    m2**2*(-m3**2 + m4**2 - 2*pap2 + 2*pbp2 + sman*(-1 + x1)) - sman**2*x1 - 2*pap2*sman*x2 + m1**2*(2*m2**2 + m3**2 - m4**2 + 2*pap2 + 4*pbp1 - 6*pbp2 + sman*(-1 + x1 + 2*x2)))))) + 
     &         A3Rc*(H3R*mq3*(H3Rc*mq3c*(-2*m3**2*pap1 + 2*m4**2*pap1 + 4*pap2*pbp1 + 4*pap1*pbp2 - m3**2*sman + m4**2*sman - 2*pbp1*sman + sman**2 + m1**2*(m3**2 - m4**2 + 2*pap1 - 2*pbp1 + sman*(-1 + x2)) + 
     &                  m2**2*(m3**2 - m4**2 - 2*pap1 + 2*pbp1 + sman*(-1 + x2)) - sman*(2*pap1 + sman)*x2) + 
     &               H3Lc*m1*(-m3**4 + m4**4 + 2*m4**2*pap1 - 4*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 2*m4**2*sman - 2*pbp1*sman + sman**2 - 2*m2**2*(m4**2 + pap1 - pbp1 + sman) + 
     &                  m3**2*(-2*pap1 + 4*pbp2 + sman*(-2 + x1)) + m2**2*sman*x1 - m4**2*sman*x1 - 2*pbp2*sman*x1 + m1**2*(2*(m3**2 + pap1 - pbp1) - sman*x1) + sman*(2*m2**2 - 2*(m4**2 + pap1 - pbp1 + sman) + sman*x1)*x2 + 
     &                  sman**2*x2**2)) + H3L*(H3Rc*m1*mq3c*(-m3**4 + m4**4 + 2*m4**2*pap1 - 4*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 2*m4**2*sman - 2*pbp1*sman + sman**2 - 2*m2**2*(m4**2 + pap1 - pbp1 + sman) + 
     &                  m3**2*(-2*pap1 + 4*pbp2 + sman*(-2 + x1)) + m2**2*sman*x1 - m4**2*sman*x1 - 2*pbp2*sman*x1 + m1**2*(2*(m3**2 + pap1 - pbp1) - sman*x1) + sman*(2*m2**2 - 2*(m4**2 + pap1 - pbp1 + sman) + sman*x1)*x2 + 
     &                  sman**2*x2**2) + H3Lc*(2*m3**2*m4**2*pap1 - 2*m4**4*pap1 + 2*m3**4*pap2 - 2*m4**4*pap2 + 4*m4**2*pap2*pbp1 - 4*m4**2*pap1*pbp2 - 8*m3**2*pap2*pbp2 + m3**2*m4**2*sman - m4**4*sman + 2*m3**2*pap2*sman - 
     &                  2*m4**2*pap2*sman + 2*m4**2*pbp1*sman - m4**2*sman**2 - 2*m3**2*pap2*sman*x1 + 2*m4**2*pap2*sman*x1 + 4*pap2*pbp2*sman*x1 + 
     &                  sman*(m4**2*(2*pap1 + 4*pap2 + sman) + 2*pap2*(-2*pbp1 + sman - sman*x1))*x2 - 2*pap2*sman**2*x2**2 + m1**4*(3*m3**2 + m4**2 + 2*pap1 - 2*pbp1 + sman - 2*sman*x1 - sman*x2) - 
     &                  m1**2*(2*m3**4 - 3*m4**4 + 2*m4**2*pap2 + 6*m4**2*pbp1 - 4*pap2*pbp1 - 4*pap1*pbp2 - 4*m4**2*sman + 2*pap2*sman + 2*pbp1*sman - sman**2 + 2*m4**2*sman*x1 - 2*pap2*sman*x1 + 4*pbp2*sman*x1 + 
     &                     m3**2*(m4**2 + 2*pap1 + 2*pap2 - 8*pbp2 + 3*sman - 2*sman*x1) + m2**2*(m3**2 + 3*m4**2 + 2*pap1 - 2*pbp1 + sman*(3 - 2*x1 - 3*x2)) + 
     &                     sman*(5*m4**2 + 2*pap1 - 2*pap2 - 4*pbp1 + 3*sman - 2*sman*x1)*x2 - 2*sman**2*x2**2) + 
     &                  m2**2*(m4**4 - m3**2*(m4**2 - 2*pap2) - 2*pap2*sman*(-1 + x1 + x2) + m4**2*(2*pap1 + 2*pap2 - 2*pbp1 + sman - sman*x2)))))) + 
     &      A3L*(-(A3Rc*m3*(H3Lc*(-2*H3R*m1*(m1**2 + m4**2 - 2*pap2)*(m1**2 - m2**2 - m3**2 + m4**2 - 2*pbp1 + 2*pbp2 - sman*x2) - 
     &                 H3L*mq3*(2*m1**4 + 2*m3**2*pap2 - 2*m4**2*pap2 + 4*pap2*pbp1 + 4*pap1*pbp2 + m3**2*sman - m4**2*sman - 2*pbp2*sman + sman**2 + m2**2*(-m3**2 + m4**2 + 2*pap2 + 2*pbp2 + sman*(-1 + x1)) - sman**2*x1 + 
     &                    2*pap2*sman*x2 - m1**2*(2*m2**2 + 3*m3**2 - 3*m4**2 + 2*pap2 + 4*pbp1 - 2*pbp2 + sman - sman*x1 + 2*sman*x2))) + 
     &              H3Rc*mq3c*(2*H3L*m1*mq3*(-m1**2 + m2**2 + m3**2 - m4**2 + 2*pbp1 - 2*pbp2 + sman*x2) + 
     &                 H3R*(-2*m1**4 - 2*m3**2*pap2 + 2*m4**2*pap2 - 4*pap2*pbp1 + 4*pap1*pbp2 + 8*pap2*pbp2 + m3**2*sman - m4**2*sman - 2*pbp2*sman + sman**2 + m2**2*(-m3**2 + m4**2 - 2*pap2 + 2*pbp2 + sman*(-1 + x1)) - 
     &                    sman**2*x1 - 2*pap2*sman*x2 + m1**2*(2*m2**2 + m3**2 - m4**2 + 2*pap2 + 4*pbp1 - 6*pbp2 + sman*(-1 + x1 + 2*x2)))))) + 
     &         A3Lc*(H3L*mq3*(H3Lc*mq3c*(-2*m3**2*pap1 + 2*m4**2*pap1 + 4*pap2*pbp1 + 4*pap1*pbp2 - m3**2*sman + m4**2*sman - 2*pbp1*sman + sman**2 + m1**2*(m3**2 - m4**2 + 2*pap1 - 2*pbp1 + sman*(-1 + x2)) + 
     &                  m2**2*(m3**2 - m4**2 - 2*pap1 + 2*pbp1 + sman*(-1 + x2)) - sman*(2*pap1 + sman)*x2) + 
     &               H3Rc*m1*(-m3**4 + m4**4 + 2*m4**2*pap1 - 4*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 2*m4**2*sman - 2*pbp1*sman + sman**2 - 2*m2**2*(m4**2 + pap1 - pbp1 + sman) + 
     &                  m3**2*(-2*pap1 + 4*pbp2 + sman*(-2 + x1)) + m2**2*sman*x1 - m4**2*sman*x1 - 2*pbp2*sman*x1 + m1**2*(2*(m3**2 + pap1 - pbp1) - sman*x1) + sman*(2*m2**2 - 2*(m4**2 + pap1 - pbp1 + sman) + sman*x1)*x2 + 
     &                  sman**2*x2**2)) + H3R*(H3Lc*m1*mq3c*(-m3**4 + m4**4 + 2*m4**2*pap1 - 4*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 2*m4**2*sman - 2*pbp1*sman + sman**2 - 2*m2**2*(m4**2 + pap1 - pbp1 + sman) + 
     &                  m3**2*(-2*pap1 + 4*pbp2 + sman*(-2 + x1)) + m2**2*sman*x1 - m4**2*sman*x1 - 2*pbp2*sman*x1 + m1**2*(2*(m3**2 + pap1 - pbp1) - sman*x1) + sman*(2*m2**2 - 2*(m4**2 + pap1 - pbp1 + sman) + sman*x1)*x2 + 
     &                  sman**2*x2**2) + H3Rc*(2*m3**2*m4**2*pap1 - 2*m4**4*pap1 + 2*m3**4*pap2 - 2*m4**4*pap2 + 4*m4**2*pap2*pbp1 - 4*m4**2*pap1*pbp2 - 8*m3**2*pap2*pbp2 + m3**2*m4**2*sman - m4**4*sman + 2*m3**2*pap2*sman - 
     &                  2*m4**2*pap2*sman + 2*m4**2*pbp1*sman - m4**2*sman**2 - 2*m3**2*pap2*sman*x1 + 2*m4**2*pap2*sman*x1 + 4*pap2*pbp2*sman*x1 + 
     &                  sman*(m4**2*(2*pap1 + 4*pap2 + sman) + 2*pap2*(-2*pbp1 + sman - sman*x1))*x2 - 2*pap2*sman**2*x2**2 + m1**4*(3*m3**2 + m4**2 + 2*pap1 - 2*pbp1 + sman - 2*sman*x1 - sman*x2) - 
     &                  m1**2*(2*m3**4 - 3*m4**4 + 2*m4**2*pap2 + 6*m4**2*pbp1 - 4*pap2*pbp1 - 4*pap1*pbp2 - 4*m4**2*sman + 2*pap2*sman + 2*pbp1*sman - sman**2 + 2*m4**2*sman*x1 - 2*pap2*sman*x1 + 4*pbp2*sman*x1 + 
     &                     m3**2*(m4**2 + 2*pap1 + 2*pap2 - 8*pbp2 + 3*sman - 2*sman*x1) + m2**2*(m3**2 + 3*m4**2 + 2*pap1 - 2*pbp1 + sman*(3 - 2*x1 - 3*x2)) + 
     &                     sman*(5*m4**2 + 2*pap1 - 2*pap2 - 4*pbp1 + 3*sman - 2*sman*x1)*x2 - 2*sman**2*x2**2) + 
     &                  m2**2*(m4**4 - m3**2*(m4**2 - 2*pap2) - 2*pap2*sman*(-1 + x1 + x2) + m4**2*(2*pap1 + 2*pap2 - 2*pbp1 + sman - sman*x2))))))))/
     &  ((m1**2 + m4**2 - mq3**2 - 2*pap2)*(m1**2 + m4**2 - mq3c**2 - 2*pap2)*(m1**2 - m2**2 + 2*(pbp1 + pbp2) - sman)*(-m3**2 + m4**2 + sman - sman*x2)))


           NeuQ2qx_M2RadiationHUU = NeuQ2qx_M2RadiationHUU + CF*NC*(radUiUi + radUfUf + 2d0*radUiUf)
!            print*,"NeuQ2qx_M2RadiationHUU",NeuQ2qx_M2RadiationHUU


           enddo
	enddo

	end


! =================================================================
	double precision function NeuQ2qx_M2RadiationHST()

        implicit none


#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

	double precision radSiTi, radSiTp, radSiTf, radSpTi, radSpTp, radSpTf, radSfTi, radSfTp, radSfTf 
        double complex A1L, A1R, F1L, F1R, A2Lc, A2Rc, G2c
        double complex g0, g1
        double precision mq1, mq2c
	integer sumTc, sumTcend

        NeuQ2qx_M2RadiationHST = 0d0

        radSiTi = 0d0       
        radSiTp = 0d0   
        radSiTf = 0d0    
        radSpTi = 0d0
        radSpTp = 0d0
        radSpTf = 0d0
        radSfTi = 0d0
        radSfTp = 0d0
        radSfTf = 0d0

	sumTc = 0

	if ((xtype.eq.5).or.(xtype.eq.6)) then

	sumTcend = 1

	elseif (xtype.eq.7) then

	sumTcend = 0

	elseif (xtype.eq.8) then

	sumTcend = 1

        endif

	do sumTc = 0,sumTcend

           A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
           A1R = AR(chimeigen,sfmeigen,sftype,sfgen)
           A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
           A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))

	   F1L = FL(xtype,pstype,psgen)
	   F1R = FR(xtype,pstype,psgen)
	   G2c = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

           g0 = -gs
	   g1 = -gs
        
           mq1 = Mf(pstype,psgen) 
           mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)


           radSiTi = dreal((-2*g0**2*G2c*(m1**2 + m2**2 + 2*(pbp1 + pbp2) - sman)*(A1R*(2*A2Rc*(F1R*mq1*pap1 + F1L*m3*(pap1 + pap2)) + A2Lc*m1*(2*F1R*m3*mq1 + F1L*(m3**2 - m4**2 + sman*(-1 + x1 + x2)))) + 
     &      A1L*(2*A2Lc*(F1L*mq1*pap1 + F1R*m3*(pap1 + pap2)) + A2Rc*m1*(2*F1L*m3*mq1 + F1R*(m3**2 - m4**2 + sman*(-1 + x1 + x2))))))/
     &  ((m1**2 + m3**2 - mq2c**2 - 2*pap1)*(-m1**2 + m2**2 - 2*(pbp1 + pbp2) + sman)**2*(-mq1**2 + sman*(-1 + x1 + x2))))


           radSiTp = dreal((g0**2*G2c*(4*m1**2 + 4*m2**2 + 2*m3**2 - 2*m4**2 - 2*pap1 - 2*pap2 + 6*pbp1 - 2*pbp2 + sman*(-2 - x1 + x2))*
     &    (A1R*(2*A2Rc*(F1R*mq1*pap1 + F1L*m3*(pap1 + pap2)) + A2Lc*m1*(2*F1R*m3*mq1 + F1L*(m3**2 - m4**2 + sman*(-1 + x1 + x2)))) + 
     &      A1L*(2*A2Lc*(F1L*mq1*pap1 + F1R*m3*(pap1 + pap2)) + A2Rc*m1*(2*F1L*m3*mq1 + F1R*(m3**2 - m4**2 + sman*(-1 + x1 + x2))))))/
     &  (2.*(m1**2 + m3**2 - mq2c**2 - 2*pap1)*(-m2**2 - m4**2 + mq2c**2 + 2*pbp2)*(m1**2 - m2**2 + 2*(pbp1 + pbp2) - sman)*(-mq1**2 + sman*(-1 + x1 + x2))))


           radSiTf = dreal((g0*g1*G2c*(A1R*(A2Rc*(F1R*mq1*(2*m3**2*pap1 - 2*m4**2*pap1 - 4*pap2*pbp1 - 4*pap1*pbp2 + m3**2*sman - m4**2*sman + 2*pbp1*sman - sman**2 + sman*(2*pap1 + sman)*x2 + 
     &               m2**2*(-m3**2 + m4**2 + 2*pap1 - 2*pbp1 + sman - sman*x2) + m1**2*(-m3**2 + m4**2 - 2*pap1 + 2*pbp1 + sman - sman*x2)) + 
     &            F1L*m3*(2*m3**2*pap1 - 2*m4**2*pap1 + 2*m3**2*pap2 - 2*m4**2*pap2 - 8*pap1*pbp2 - 8*pap2*pbp2 + 2*pbp1*sman + 2*pbp2*sman - 2*sman**2 + sman**2*x1 + sman*(2*(pap1 + pap2) + sman)*x2 + 
     &               m2**2*(2*pap1 + 2*pap2 - 2*pbp1 - 2*pbp2 - sman*(-2 + x1 + x2)) - m1**2*(2*pap1 + 2*pap2 - 2*pbp1 - 2*pbp2 + sman*(-2 + x1 + x2)))) + 
     &         A2Lc*m1*(2*F1R*m3*mq1*(-m1**2 + m2**2 + m3**2 - m4**2 + 2*pbp1 - 2*pbp2 + sman*x2) + 
     &            F1L*(m3**4 - 2*m3**2*m4**2 + m4**4 + 4*m3**2*pbp1 - 4*m4**2*pbp1 - m3**2*sman + m4**2*sman + m3**2*sman*x1 - m4**2*sman*x1 - 2*pbp2*sman*x1 + sman*(2*m3**2 - 2*m4**2 + 2*pbp1 + sman*(-1 + x1))*x2 + 
     &               sman**2*x2**2 + m1**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + m2**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2))))) + 
     &      A1L*(A2Lc*(F1L*mq1*(2*m3**2*pap1 - 2*m4**2*pap1 - 4*pap2*pbp1 - 4*pap1*pbp2 + m3**2*sman - m4**2*sman + 2*pbp1*sman - sman**2 + sman*(2*pap1 + sman)*x2 + 
     &               m2**2*(-m3**2 + m4**2 + 2*pap1 - 2*pbp1 + sman - sman*x2) + m1**2*(-m3**2 + m4**2 - 2*pap1 + 2*pbp1 + sman - sman*x2)) + 
     &            F1R*m3*(2*m3**2*pap1 - 2*m4**2*pap1 + 2*m3**2*pap2 - 2*m4**2*pap2 - 8*pap1*pbp2 - 8*pap2*pbp2 + 2*pbp1*sman + 2*pbp2*sman - 2*sman**2 + sman**2*x1 + sman*(2*(pap1 + pap2) + sman)*x2 + 
     &               m2**2*(2*pap1 + 2*pap2 - 2*pbp1 - 2*pbp2 - sman*(-2 + x1 + x2)) - m1**2*(2*pap1 + 2*pap2 - 2*pbp1 - 2*pbp2 + sman*(-2 + x1 + x2)))) + 
     &         A2Rc*m1*(2*F1L*m3*mq1*(-m1**2 + m2**2 + m3**2 - m4**2 + 2*pbp1 - 2*pbp2 + sman*x2) + 
     &            F1R*(m3**4 - 2*m3**2*m4**2 + m4**4 + 4*m3**2*pbp1 - 4*m4**2*pbp1 - m3**2*sman + m4**2*sman + m3**2*sman*x1 - m4**2*sman*x1 - 2*pbp2*sman*x1 + sman*(2*m3**2 - 2*m4**2 + 2*pbp1 + sman*(-1 + x1))*x2 + 
     &               sman**2*x2**2 + m1**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + m2**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2)))))))/
     &  ((m2**2 + m4**2 - mq2c**2 - 2*pbp2)*(-m1**2 + m2**2 - 2*(pbp1 + pbp2) + sman)*(-m3**2 + m4**2 + sman - sman*x2)*(-mq1**2 + sman*(-1 + x1 + x2))))


           radSpTi = dreal( (g0*g1*G2c*(A1L*(2*A2Rc*m1*(F1L*m3*mq1*(2*m1**2 - 2*m2**2 - 2*pap1 - 2*pap2 - 6*pbp1 - 6*pbp2 - sman*(-2 + x1 + x2)) + 
     &            F1R*(-(m3**2*mq1**2) + m4**2*mq1**2 - 4*mq1**2*pbp1 + 4*pap2*pbp1 - 4*pap1*pbp2 + mq1**2*sman + 2*pap1*sman + 2*pbp1*sman - 2*pap1*sman*x1 - pap2*sman*x1 - 2*pbp1*sman*x1 - pbp2*sman*x1 - 
     &               (mq1**2 + pap1 + pbp1)*sman*x2 + m2**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + m1**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2)))) + 
     &         A2Lc*(F1R*m3*(-4*mq1**2*pap1 - 4*mq1**2*pap2 - 2*mq1**2*sman + 2*pap1*sman + 2*pap2*sman - 2*pbp1*sman - 2*pbp2*sman + 2*sman**2 - sman**2*x1 - sman**2*x2 + 
     &               m1**2*(6*mq1**2 - 2*pap1 - 2*pap2 - 10*pbp1 - 10*pbp2 - sman*(-2 + x1 + x2)) + m2**2*(2*mq1**2 - 6*pap1 - 6*pap2 + 2*pbp1 + 2*pbp2 + sman*(-2 + x1 + x2))) + 
     &            F1L*mq1*(12*pap2*pbp1 - 12*pap1*pbp2 + 2*(-m3**2 + m4**2 + 3*pap1)*sman - 2*pbp1*sman + 2*sman**2 - 4*pap1*sman*x1 - 2*pap2*sman*x1 - sman**2*x1 - 2*sman*(pap1 + sman)*x2 + 
     &               m2**2*(2*m3**2 - 2*m4**2 - 6*pap1 + 2*pbp1 + sman*(-2 + x1 + 2*x2)) + m1**2*(2*m3**2 - 2*m4**2 - 2*pap1 - 10*pbp1 + sman*(-2 + 3*x1 + 2*x2))))) + 
     &      A1R*(2*A2Lc*m1*(F1R*m3*mq1*(2*m1**2 - 2*m2**2 - 2*pap1 - 2*pap2 - 6*pbp1 - 6*pbp2 - sman*(-2 + x1 + x2)) + 
     &            F1L*(-(m3**2*mq1**2) + m4**2*mq1**2 - 4*mq1**2*pbp1 + 4*pap2*pbp1 - 4*pap1*pbp2 + mq1**2*sman + 2*pap1*sman + 2*pbp1*sman - 2*pap1*sman*x1 - pap2*sman*x1 - 2*pbp1*sman*x1 - pbp2*sman*x1 - 
     &               (mq1**2 + pap1 + pbp1)*sman*x2 + m2**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + m1**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2)))) + 
     &         A2Rc*(F1L*m3*(-4*mq1**2*pap1 - 4*mq1**2*pap2 - 2*mq1**2*sman + 2*pap1*sman + 2*pap2*sman - 2*pbp1*sman - 2*pbp2*sman + 2*sman**2 - sman**2*x1 - sman**2*x2 + 
     &               m1**2*(6*mq1**2 - 2*pap1 - 2*pap2 - 10*pbp1 - 10*pbp2 - sman*(-2 + x1 + x2)) + m2**2*(2*mq1**2 - 6*pap1 - 6*pap2 + 2*pbp1 + 2*pbp2 + sman*(-2 + x1 + x2))) + 
     &            F1R*mq1*(12*pap2*pbp1 - 12*pap1*pbp2 + 2*(-m3**2 + m4**2 + 3*pap1)*sman - 2*pbp1*sman + 2*sman**2 - 4*pap1*sman*x1 - 2*pap2*sman*x1 - sman**2*x1 - 2*sman*(pap1 + sman)*x2 + 
     &               m2**2*(2*m3**2 - 2*m4**2 - 6*pap1 + 2*pbp1 + sman*(-2 + x1 + 2*x2)) + m1**2*(2*m3**2 - 2*m4**2 - 2*pap1 - 10*pbp1 + sman*(-2 + 3*x1 + 2*x2)))))))/
     &  (2.*(m1**2 + m3**2 - mq2c**2 - 2*pap1)*(mq1**2 - sman)*(m1**2 - m2**2 + 2*(pbp1 + pbp2) - sman)*(mq1**2 - sman*(-1 + x1 + x2))))


           radSpTp = dreal((g0*g1*G2c*(A1R*(A2Rc*(F1L*m3*(2*mq1**2*pap1 - 2*mq1**2*pap2 + m2**2*(-m3**2 + m4**2 - mq1**2 + 4*(pap1 + pap2)) + 4*pap2*pbp1 - 4*pap1*pbp2 + m1**2*(m3**2 - m4**2 - 3*mq1**2 + 4*(pbp1 + pbp2)) + 
     &               (m3**2 - m4**2 + mq1**2 - 2*(pap1 + pap2))*sman) + 2*F1R*mq1*(-(m2**2*(m3**2 - 2*pap1)) + m3**2*(-2*pap2 + sman) + pap1*(-2*m4**2 + sman*(-2 + x1 + x2)) + m1**2*(m4**2 + 2*pbp1 - sman*(-1 + x1 + x2)))) + 
     &         A2Lc*m1*(2*F1R*m3*mq1*(-m1**2 + m2**2 + m3**2 - m4**2 - 2*pap2 + 2*pbp1) + F1L*
     &             (3*m3**2*mq1**2 + m4**2*mq1**2 - 4*m4**2*pap1 - 2*mq1**2*pap1 - 4*m3**2*pap2 - 4*m4**2*pbp1 + 2*mq1**2*pbp1 - 4*pap2*pbp1 - 4*m3**2*pbp2 + 4*pap1*pbp2 + mq1**2*sman - 2*pap1*sman - 2*pbp1*sman - 
     &               mq1**2*sman*x1 + 2*pap1*sman*x1 + 2*pbp1*sman*x1 + (-mq1**2 + 2*(pap1 + pbp1))*sman*x2 + m1**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + m2**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2))))) + 
     &      A1L*(A2Lc*(F1R*m3*(2*mq1**2*pap1 - 2*mq1**2*pap2 + m2**2*(-m3**2 + m4**2 - mq1**2 + 4*(pap1 + pap2)) + 4*pap2*pbp1 - 4*pap1*pbp2 + m1**2*(m3**2 - m4**2 - 3*mq1**2 + 4*(pbp1 + pbp2)) + 
     &               (m3**2 - m4**2 + mq1**2 - 2*(pap1 + pap2))*sman) + 2*F1L*mq1*(-(m2**2*(m3**2 - 2*pap1)) + m3**2*(-2*pap2 + sman) + pap1*(-2*m4**2 + sman*(-2 + x1 + x2)) + m1**2*(m4**2 + 2*pbp1 - sman*(-1 + x1 + x2)))) + 
     &         A2Rc*m1*(2*F1L*m3*mq1*(-m1**2 + m2**2 + m3**2 - m4**2 - 2*pap2 + 2*pbp1) + F1R*
     &             (3*m3**2*mq1**2 + m4**2*mq1**2 - 4*m4**2*pap1 - 2*mq1**2*pap1 - 4*m3**2*pap2 - 4*m4**2*pbp1 + 2*mq1**2*pbp1 - 4*pap2*pbp1 - 4*m3**2*pbp2 + 4*pap1*pbp2 + mq1**2*sman - 2*pap1*sman - 2*pbp1*sman - 
     &               mq1**2*sman*x1 + 2*pap1*sman*x1 + 2*pbp1*sman*x1 + (-mq1**2 + 2*(pap1 + pbp1))*sman*x2 + m1**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + m2**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2)))))))/
     &  ((m1**2 + m3**2 - mq2c**2 - 2*pap1)*(-m2**2 - m4**2 + mq2c**2 + 2*pbp2)*(mq1**2 - sman)*(mq1**2 - sman*(-1 + x1 + x2))))


           radSpTf = dreal( (g1**2*G2c*(A1L*(2*A2Rc*m1*(F1L*m3*mq1*(-m3**2 + m4**2 + 2*pap1 + 4*pap2 + 2*pbp1 + 4*pbp2 - sman*(1 + 2*x1 + 2*x2)) + 
     &            F1R*(m4**2*(mq1**2 - 2*(pap2 + pbp2 - sman)) + m3**2*(-3*mq1**2 + 2*(pap1 + 2*pap2 + pbp1 + 2*pbp2 - sman)) + sman*(-(mq1**2*(-1 + x2)) + 2*(pap2 + pbp2 - sman)*(-1 + x1 + x2)))) + 
     &         A2Lc*(F1R*m3*(4*mq1**2*pap1 + 8*mq1**2*pap2 - 4*pap2*pbp1 + 4*pap1*pbp2 - m3**2*sman + m4**2*sman - 4*mq1**2*sman - 2*pap1*sman - 2*pap2*sman + 2*pbp1*sman + 2*pbp2*sman + 3*sman**2 - 2*sman**2*x1 + 
     &               m1**2*(-m3**2 + m4**2 - 4*mq1**2 + 2*pap1 + 2*pap2 + 2*pbp1 + 2*pbp2 + sman*(3 - 2*x1 - 2*x2)) - 2*sman**2*x2 + 
     &               m2**2*(m3**2 - m4**2 + 4*mq1**2 - 2*pap1 - 2*pap2 - 2*pbp1 - 2*pbp2 + sman*(-3 + 2*x1 + 2*x2))) + 
     &            F1L*mq1*(4*m3**2*pap1 + 8*m3**2*pap2 - 4*m4**2*pap2 - 4*pap2*pbp1 + 4*pap1*pbp2 - 5*m3**2*sman + 3*m4**2*sman - 2*pap1*sman - 4*pap2*sman + 2*pbp1*sman + 3*sman**2 + 4*pap2*sman*x1 - 2*sman**2*x1 + 
     &               m1**2*(-5*m3**2 + 3*m4**2 + 2*pap1 + 2*pbp1 + sman*(3 - 2*x1 - 3*x2)) + (4*pap2 - 3*sman)*sman*x2 + m2**2*(5*m3**2 - 3*m4**2 - 2*pap1 - 2*pbp1 + sman*(-3 + 2*x1 + 3*x2))))) + 
     &      A1R*(2*A2Lc*m1*(F1R*m3*mq1*(-m3**2 + m4**2 + 2*pap1 + 4*pap2 + 2*pbp1 + 4*pbp2 - sman*(1 + 2*x1 + 2*x2)) + 
     &            F1L*(m4**2*(mq1**2 - 2*(pap2 + pbp2 - sman)) + m3**2*(-3*mq1**2 + 2*(pap1 + 2*pap2 + pbp1 + 2*pbp2 - sman)) + sman*(-(mq1**2*(-1 + x2)) + 2*(pap2 + pbp2 - sman)*(-1 + x1 + x2)))) + 
     &         A2Rc*(F1L*m3*(4*mq1**2*pap1 + 8*mq1**2*pap2 - 4*pap2*pbp1 + 4*pap1*pbp2 - m3**2*sman + m4**2*sman - 4*mq1**2*sman - 2*pap1*sman - 2*pap2*sman + 2*pbp1*sman + 2*pbp2*sman + 3*sman**2 - 2*sman**2*x1 + 
     &               m1**2*(-m3**2 + m4**2 - 4*mq1**2 + 2*pap1 + 2*pap2 + 2*pbp1 + 2*pbp2 + sman*(3 - 2*x1 - 2*x2)) - 2*sman**2*x2 + 
     &               m2**2*(m3**2 - m4**2 + 4*mq1**2 - 2*pap1 - 2*pap2 - 2*pbp1 - 2*pbp2 + sman*(-3 + 2*x1 + 2*x2))) + 
     &            F1R*mq1*(4*m3**2*pap1 + 8*m3**2*pap2 - 4*m4**2*pap2 - 4*pap2*pbp1 + 4*pap1*pbp2 - 5*m3**2*sman + 3*m4**2*sman - 2*pap1*sman - 4*pap2*sman + 2*pbp1*sman + 3*sman**2 + 4*pap2*sman*x1 - 2*sman**2*x1 + 
     &               m1**2*(-5*m3**2 + 3*m4**2 + 2*pap1 + 2*pbp1 + sman*(3 - 2*x1 - 3*x2)) + (4*pap2 - 3*sman)*sman*x2 + m2**2*(5*m3**2 - 3*m4**2 - 2*pap1 - 2*pbp1 + sman*(-3 + 2*x1 + 3*x2)))))))/
     &  ((m2**2 + m4**2 - mq2c**2 - 2*pbp2)*(mq1**2 - sman)*(-m3**2 + m4**2 + sman - sman*x2)*(mq1**2 - sman*(-1 + x1 + x2))))


           radSfTi = dreal((g0*g1*G2c*(A1R*(A2Lc*m1*(-2*F1R*m3*mq1*(-m1**2 + m2**2 + m3**2 - m4**2 + 2*pbp1 - 2*pbp2 + sman*x2) + 
     &            F1L*(-2*m3**2*pap1 + 2*m4**2*pap1 - 2*m3**2*pbp1 + 2*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 8*pbp1*pbp2 - m3**2*sman + m4**2*sman - 4*pbp1*sman + sman**2 + 
     &               m1**2*(m3**2 - m4**2 + 2*pap1 + 2*pbp1 + sman*(-1 + x2)) - sman*(2*(pap1 + pbp1) + sman)*x2 + m2**2*(-m3**2 + m4**2 - 2*pap1 - 2*pbp1 + sman - sman*x2))) - 
     &         A2Rc*(F1R*mq1*(2*m3**2*pap1 - 2*m4**2*pap1 - 4*pap2*pbp1 - 4*pap1*pbp2 + m3**2*sman - m4**2*sman + 2*pbp1*sman - sman**2 + sman*(2*pap1 + sman)*x2 + m2**2*(-m3**2 + m4**2 + 2*pap1 - 2*pbp1 + sman - sman*x2) + 
     &               m1**2*(-m3**2 + m4**2 - 2*pap1 + 2*pbp1 + sman - sman*x2)) + F1L*m3*(m1**2*(m3**2 - m4**2 + 4*pbp1 + sman*(-3 + x2)) + sman*(m3**2 - m4**2 - 4*pbp2 + sman + sman*x2) - 
     &               m2**2*(m3**2 - m4**2 - 4*pap1 - 4*pap2 - 4*pbp2 + sman + sman*x2)))) + 
     &      A1L*(A2Rc*m1*(-2*F1L*m3*mq1*(-m1**2 + m2**2 + m3**2 - m4**2 + 2*pbp1 - 2*pbp2 + sman*x2) + 
     &            F1R*(-2*m3**2*pap1 + 2*m4**2*pap1 - 2*m3**2*pbp1 + 2*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 8*pbp1*pbp2 - m3**2*sman + m4**2*sman - 4*pbp1*sman + sman**2 + 
     &               m1**2*(m3**2 - m4**2 + 2*pap1 + 2*pbp1 + sman*(-1 + x2)) - sman*(2*(pap1 + pbp1) + sman)*x2 + m2**2*(-m3**2 + m4**2 - 2*pap1 - 2*pbp1 + sman - sman*x2))) - 
     &         A2Lc*(F1L*mq1*(2*m3**2*pap1 - 2*m4**2*pap1 - 4*pap2*pbp1 - 4*pap1*pbp2 + m3**2*sman - m4**2*sman + 2*pbp1*sman - sman**2 + sman*(2*pap1 + sman)*x2 + m2**2*(-m3**2 + m4**2 + 2*pap1 - 2*pbp1 + sman - sman*x2) + 
     &               m1**2*(-m3**2 + m4**2 - 2*pap1 + 2*pbp1 + sman - sman*x2)) + F1R*m3*(m1**2*(m3**2 - m4**2 + 4*pbp1 + sman*(-3 + x2)) + sman*(m3**2 - m4**2 - 4*pbp2 + sman + sman*x2) - 
     &               m2**2*(m3**2 - m4**2 - 4*pap1 - 4*pap2 - 4*pbp2 + sman + sman*x2))))))/((m1**2 + m3**2 - mq2c**2 - 2*pap1)*(mq1**2 - sman)*(m1**2 - m2**2 + 2*(pbp1 + pbp2) - sman)*(m3**2 - m4**2 + sman*(-1 + x2))))


           radSfTp = dreal((g0*g1*G2c*(A1R*(-2*A2Lc*m1*(F1R*m3*mq1*(2*m1**2 - 2*m2**2 - 4*m3**2 - 4*m4**2 + 2*pap1 - 2*pap2 - 2*pbp1 + 2*pbp2 + sman*(-2 + x1 + 3*x2)) + 
     &            F1L*(-2*m3**2*pap1 - 2*m4**2*pap1 + 4*m3**2*pap2 - 4*pap1*pap2 - 2*m3**2*pbp1 - 2*m4**2*pbp1 + 4*m3**2*pbp2 + 4*pbp1*pbp2 - 3*m3**2*sman - m4**2*sman + 2*pap1*sman - 2*pbp1*sman - sman**2 - pap2*sman*x1 - 
     &               pbp2*sman*x1 + sman**2*x1 + m1**2*(m3**2 - m4**2 + 2*pap1 + 2*pbp1 + sman*(-1 + x2)) + sman*(pap1 + pbp1 + sman)*x2 + m2**2*(-m3**2 + m4**2 - 2*pap1 - 2*pbp1 + sman - sman*x2))) + 
     &         A2Rc*(F1R*mq1*(4*m4**2*pap1 + 8*pap1*pap2 - 4*pap2*pbp1 - 4*pap1*pbp2 - 2*pap1*sman + 2*pbp1*sman + 4*m3**2*(pap1 - 2*pap2 + sman) + 2*m1**2*(2*m4**2 - 3*pap1 + pbp1 + 2*sman) - m1**2*sman*x1 + 
     &               2*pap2*sman*x1 - sman**2*x1 + m2**2*(-4*m3**2 + 6*pap1 - 2*pbp1 + sman*x1) - 2*(2*m1**2 + pap1)*sman*x2) + 
     &            F1L*m3*(8*pap2*pbp1 + 4*(m3**2 + m4**2)*sman - 2*pap2*sman - 2*pbp1*sman - 2*pbp2*sman + 4*sman**2 - 2*pap1*(4*pbp2 + sman) - sman**2*x1 - 3*sman**2*x2 + 
     &               m2**2*(-4*m3**2 - 4*m4**2 + 10*pap1 + 2*pap2 + 2*pbp1 + 2*pbp2 + sman*(-4 + x1 + 3*x2)) + m1**2*(4*m3**2 + 4*m4**2 - 2*pap1 + 6*pap2 + 6*pbp1 + 6*pbp2 - sman*(4 + x1 + 3*x2))))) + 
     &      A1L*(-2*A2Rc*m1*(F1L*m3*mq1*(2*m1**2 - 2*m2**2 - 4*m3**2 - 4*m4**2 + 2*pap1 - 2*pap2 - 2*pbp1 + 2*pbp2 + sman*(-2 + x1 + 3*x2)) + 
     &            F1R*(-2*m3**2*pap1 - 2*m4**2*pap1 + 4*m3**2*pap2 - 4*pap1*pap2 - 2*m3**2*pbp1 - 2*m4**2*pbp1 + 4*m3**2*pbp2 + 4*pbp1*pbp2 - 3*m3**2*sman - m4**2*sman + 2*pap1*sman - 2*pbp1*sman - sman**2 - pap2*sman*x1 - 
     &               pbp2*sman*x1 + sman**2*x1 + m1**2*(m3**2 - m4**2 + 2*pap1 + 2*pbp1 + sman*(-1 + x2)) + sman*(pap1 + pbp1 + sman)*x2 + m2**2*(-m3**2 + m4**2 - 2*pap1 - 2*pbp1 + sman - sman*x2))) + 
     &         A2Lc*(F1L*mq1*(4*m4**2*pap1 + 8*pap1*pap2 - 4*pap2*pbp1 - 4*pap1*pbp2 - 2*pap1*sman + 2*pbp1*sman + 4*m3**2*(pap1 - 2*pap2 + sman) + 2*m1**2*(2*m4**2 - 3*pap1 + pbp1 + 2*sman) - m1**2*sman*x1 + 
     &               2*pap2*sman*x1 - sman**2*x1 + m2**2*(-4*m3**2 + 6*pap1 - 2*pbp1 + sman*x1) - 2*(2*m1**2 + pap1)*sman*x2) + 
     &            F1R*m3*(8*pap2*pbp1 + 4*(m3**2 + m4**2)*sman - 2*pap2*sman - 2*pbp1*sman - 2*pbp2*sman + 4*sman**2 - 2*pap1*(4*pbp2 + sman) - sman**2*x1 - 3*sman**2*x2 + 
     &               m2**2*(-4*m3**2 - 4*m4**2 + 10*pap1 + 2*pap2 + 2*pbp1 + 2*pbp2 + sman*(-4 + x1 + 3*x2)) + m1**2*(4*m3**2 + 4*m4**2 - 2*pap1 + 6*pap2 + 6*pbp1 + 6*pbp2 - sman*(4 + x1 + 3*x2)))))))/
     &  (2.*(m1**2 + m3**2 - mq2c**2 - 2*pap1)*(-m2**2 - m4**2 + mq2c**2 + 2*pbp2)*(mq1**2 - sman)*(m3**2 - m4**2 + sman*(-1 + x2))))


           radSfTf = dreal((2*g1**2*G2c*(A1L*(2*A2Rc*m1*(F1L*m3*mq1*(m3**2 + m4**2 + sman - sman*x2) + F1R*(-(m3**2*(pap1 + 3*pap2 + pbp1 + 3*pbp2 - 3*sman)) + (pap1 + pap2 + pbp1 + pbp2 - sman)*(m4**2 + sman - sman*x2))) + 
     &         A2Lc*(F1R*m3*(m1**2 - m2**2 + sman)*(m3**2 + m4**2 + sman - sman*x2) + F1L*mq1*
     &             (-2*m3**2*pap1 + 2*m4**2*pap1 - 6*m3**2*pap2 + 2*m4**2*pap2 + 3*m3**2*sman - m4**2*sman + 2*pap1*sman + 2*pap2*sman - sman**2 + m1**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) + 
     &               sman*(-2*(pap1 + pap2) + sman)*x2 + m2**2*(-3*m3**2 + m4**2 + sman - sman*x2)))) + 
     &      A1R*(2*A2Lc*m1*(F1R*m3*mq1*(m3**2 + m4**2 + sman - sman*x2) + F1L*(-(m3**2*(pap1 + 3*pap2 + pbp1 + 3*pbp2 - 3*sman)) + (pap1 + pap2 + pbp1 + pbp2 - sman)*(m4**2 + sman - sman*x2))) + 
     &         A2Rc*(F1L*m3*(m1**2 - m2**2 + sman)*(m3**2 + m4**2 + sman - sman*x2) + F1R*mq1*
     &             (-2*m3**2*pap1 + 2*m4**2*pap1 - 6*m3**2*pap2 + 2*m4**2*pap2 + 3*m3**2*sman - m4**2*sman + 2*pap1*sman + 2*pap2*sman - sman**2 + m1**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) + 
     &               sman*(-2*(pap1 + pap2) + sman)*x2 + m2**2*(-3*m3**2 + m4**2 + sman - sman*x2))))))/((m2**2 + m4**2 - mq2c**2 - 2*pbp2)*(mq1**2 - sman)*(-m3**2 + m4**2 + sman - sman*x2)**2))

	   NeuQ2qx_M2RadiationHST = NeuQ2qx_M2RadiationHST + CF*NC*2d0*(radSiTi + radSiTp + radSiTf + radSpTi + radSpTp + radSpTf + radSfTi + radSfTp + radSfTf)
!            print*,"NeuQ2qx_M2RadiationHST",NeuQ2qx_M2RadiationHST

	enddo

	end

! =================================================================
	double precision function NeuQ2qx_M2RadiationHSU()

        implicit none


#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double precision radSiUi, radSiUf, radSpUi, radSpUf, radSfUi, radSfUf
        double complex A1L, A1R, F1L, F1R, A3Lc, A3Rc, H3Lc, H3Rc
        double complex g0, g1
        double precision mq1, mq3c
        integer sumUc, sumUcend


        NeuQ2qx_M2RadiationHSU = 0d0

        radSiUi = 0d0
        radSiUf = 0d0
        radSpUi = 0d0
        radSpUf = 0d0
        radSfUi = 0d0
        radSfUf = 0d0

        sumUc = 0

        if ((xtype.eq.5).or.(xtype.eq.6)) then
        
           sumUcend = 3
           
        elseif (xtype.eq.7) then
        
           sumUcend = 3
           
        elseif (xtype.eq.8) then
               
           sumUcend = 1
           
         endif


        do sumUc = 0, sumUcend

        A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A1R = AR(chimeigen,sfmeigen,sftype,sfgen)
        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

	F1L = FL(xtype,pstype,psgen)
	F1R = FR(xtype,pstype,psgen)

        H3Lc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
        H3Rc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

        g0 = -gs
	g1 = -gs
        
        mq1 = Mf(pstype,psgen) 
        mq3c = MNeuCha(pumeigen+sumUc)


        radSiUi = dreal((-2*g0**2*(m1**2 + m2**2 + 2*(pbp1 + pbp2) - sman)*(A1L*(A3Lc*(F1L*mq1*(2*H3Lc*mq3c*pap1 + H3Rc*m1*(m3**2 + m4**2 + 2*pap1 - sman*(-1 + x1 + x2))) + 
     &            F1R*m3*(2*H3Lc*mq3c*(pap1 + pap2) + H3Rc*m1*(m3**2 - m4**2 + 2*pap1 + 2*pap2 - sman*(-1 + x1 + x2)))) + 
     &         A3Rc*(2*F1L*m3*mq1*(H3Rc*m1*mq3c + H3Lc*(m1**2 - pap2)) + F1R*(H3Rc*m1*mq3c*(m3**2 - m4**2 + sman*(-1 + x1 + x2)) + 
     &               H3Lc*(m1**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2)) + 2*(m4**2*(pap1 + pap2) - pap2*sman*(-1 + x1 + x2)))))) + 
     &      A1R*(A3Rc*(F1R*mq1*(2*H3Rc*mq3c*pap1 + H3Lc*m1*(m3**2 + m4**2 + 2*pap1 - sman*(-1 + x1 + x2))) + F1L*m3*(2*H3Rc*mq3c*(pap1 + pap2) + H3Lc*m1*(m3**2 - m4**2 + 2*pap1 + 2*pap2 - sman*(-1 + x1 + x2)))) + 
     &         A3Lc*(2*F1R*m3*mq1*(H3Lc*m1*mq3c + H3Rc*(m1**2 - pap2)) + F1L*(H3Lc*m1*mq3c*(m3**2 - m4**2 + sman*(-1 + x1 + x2)) + 
     &               H3Rc*(m1**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2)) + 2*(m4**2*(pap1 + pap2) - pap2*sman*(-1 + x1 + x2))))))))/
     &  ((m1**2 + m4**2 - mq3c**2 - 2*pap2)*(-m1**2 + m2**2 - 2*(pbp1 + pbp2) + sman)**2*(-mq1**2 + sman*(-1 + x1 + x2))))


        radSiUf = dreal((g0*g1*(A1L*(A3Lc*(F1L*mq1*(H3Lc*mq3c*(-2*m3**2*pap1 + 2*m4**2*pap1 + 4*pap2*pbp1 + 4*pap1*pbp2 - m3**2*sman + m4**2*sman - 2*pbp1*sman + sman**2 + m1**2*(m3**2 - m4**2 + 2*pap1 - 2*pbp1 + sman*(-1 + x2)) + 
     &                  m2**2*(m3**2 - m4**2 - 2*pap1 + 2*pbp1 + sman*(-1 + x2)) - sman*(2*pap1 + sman)*x2) + 
     &               H3Rc*m1*(-m3**4 + m4**4 + 2*m4**2*pap1 - 4*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 2*m4**2*sman - 2*pbp1*sman + sman**2 - 2*m2**2*(m4**2 + pap1 - pbp1 + sman) + 
     &                  m3**2*(-2*pap1 + 4*pbp2 + sman*(-2 + x1)) + m2**2*sman*x1 - m4**2*sman*x1 - 2*pbp2*sman*x1 + m1**2*(2*(m3**2 + pap1 - pbp1) - sman*x1) + sman*(2*m2**2 - 2*(m4**2 + pap1 - pbp1 + sman) + sman*x1)*x2 + 
     &                  sman**2*x2**2)) - F1R*m3*(H3Rc*m1*(m3**4 - 2*m3**2*m4**2 + m4**4 + 2*m3**2*pap1 - 2*m4**2*pap1 + 2*m3**2*pap2 - 2*m4**2*pap2 - 4*m3**2*pbp2 + 4*m4**2*pbp2 - 8*pap1*pbp2 - 8*pap2*pbp2 + m3**2*sman - 
     &                  m4**2*sman + 2*pbp1*sman + 2*pbp2*sman - 2*sman**2 + m1**2*(-m3**2 + m4**2 - 2*pap1 - 2*pap2 + 2*(pbp1 + pbp2) + sman) - m3**2*sman*x1 + m4**2*sman*x1 + 2*pbp2*sman*x1 + sman**2*x1 + 
     &                  m2**2*(m3**2 - m4**2 + 2*pap1 + 2*pap2 - 2*pbp1 - 2*pbp2 + sman*(3 - 2*x1 - 2*x2)) + sman*(2*(pap1 + pap2 - pbp1 + sman) - sman*x1)*x2 - sman**2*x2**2) + 
     &               H3Lc*mq3c*(2*m3**2*pap1 - 2*m4**2*pap1 + 2*m3**2*pap2 - 2*m4**2*pap2 - 8*pap1*pbp2 - 8*pap2*pbp2 + 2*pbp1*sman + 2*pbp2*sman - 2*sman**2 + sman**2*x1 + sman*(2*(pap1 + pap2) + sman)*x2 + 
     &                  m2**2*(2*pap1 + 2*pap2 - 2*pbp1 - 2*pbp2 - sman*(-2 + x1 + x2)) - m1**2*(2*pap1 + 2*pap2 - 2*pbp1 - 2*pbp2 + sman*(-2 + x1 + x2))))) + 
     &         A3Rc*(F1L*m3*mq1*(2*H3Rc*m1*mq3c*(m1**2 - m2**2 - m3**2 + m4**2 - 2*pbp1 + 2*pbp2 - sman*x2) + 
     &               H3Lc*(2*m1**4 + 2*m3**2*pap2 - 2*m4**2*pap2 + 4*pap2*pbp1 + 4*pap1*pbp2 + m3**2*sman - m4**2*sman - 2*pbp2*sman + sman**2 + m2**2*(-m3**2 + m4**2 + 2*pap2 + 2*pbp2 + sman*(-1 + x1)) - sman**2*x1 + 
     &                  2*pap2*sman*x2 - m1**2*(2*m2**2 + 3*m3**2 - 3*m4**2 + 2*pap2 + 4*pbp1 - 2*pbp2 + sman - sman*x1 + 2*sman*x2))) + 
     &            F1R*(H3Rc*m1*mq3c*(-m3**4 + 2*m3**2*m4**2 - m4**4 - 4*m3**2*pbp1 + 4*m4**2*pbp1 + m3**2*sman - m4**2*sman - m3**2*sman*x1 + m4**2*sman*x1 + 2*pbp2*sman*x1 + 
     &                  sman*(-2*m3**2 + 2*m4**2 - 2*pbp1 + sman - sman*x1)*x2 - sman**2*x2**2 + m2**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + m1**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2))) + 
     &               H3Lc*(-2*m3**2*m4**2*pap1 + 2*m4**4*pap1 - 2*m3**2*m4**2*pap2 + 2*m4**4*pap2 + 4*m3**2*pap2*pbp1 - 4*m4**2*pap2*pbp1 + 4*m3**2*pap1*pbp2 + 4*m4**2*pap1*pbp2 + 8*m3**2*pap2*pbp2 + m3**4*sman - 
     &                  2*m3**2*m4**2*sman + m4**4*sman - 2*m3**2*pap2*sman + 2*m4**2*pap2*sman - 2*m4**2*pbp1*sman - 2*m3**2*pbp2*sman + m3**2*sman**2 + m4**2*sman**2 + 2*m3**2*pap2*sman*x1 - 2*m4**2*pap2*sman*x1 - 
     &                  4*pap2*pbp2*sman*x1 - m3**2*sman**2*x1 - sman*(m4**2*(2*pap1 + 4*pap2 + sman) - 2*pap2*(m3**2 + 2*pbp1 + sman*(-1 + x1)))*x2 + 2*pap2*sman**2*x2**2 + m1**4*(m3**2 - m4**2 + sman*(-1 + x1 + x2)) - 
     &                  m2**2*(m3**4 + m4**4 + m3**2*(-2*m4**2 - 2*pbp2 + sman - sman*x1) - 2*pap2*sman*(-1 + x1 + x2) + m4**2*(2*pap1 + 2*pap2 - 2*pbp1 + sman - sman*x2)) + 
     &                  m1**2*(-2*m3**4 - 2*m4**4 + 2*m4**2*pap1 + 2*m4**2*pap2 + 2*m4**2*pbp1 - 2*m4**2*sman + 2*pap2*sman + m4**2*sman*x1 - 2*pap2*sman*x1 + 2*pbp2*sman*x1 + 
     &                     sman*(3*m4**2 - 2*pap2 - 2*pbp1 + sman - sman*x1)*x2 - sman**2*x2**2 + m2**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + m3**2*(4*m4**2 - 2*(2*pbp1 + pbp2 + sman*x2))))))) + 
     &      A1R*(A3Rc*(F1R*mq1*(H3Rc*mq3c*(-2*m3**2*pap1 + 2*m4**2*pap1 + 4*pap2*pbp1 + 4*pap1*pbp2 - m3**2*sman + m4**2*sman - 2*pbp1*sman + sman**2 + m1**2*(m3**2 - m4**2 + 2*pap1 - 2*pbp1 + sman*(-1 + x2)) + 
     &                  m2**2*(m3**2 - m4**2 - 2*pap1 + 2*pbp1 + sman*(-1 + x2)) - sman*(2*pap1 + sman)*x2) + 
     &               H3Lc*m1*(-m3**4 + m4**4 + 2*m4**2*pap1 - 4*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 2*m4**2*sman - 2*pbp1*sman + sman**2 - 2*m2**2*(m4**2 + pap1 - pbp1 + sman) + 
     &                  m3**2*(-2*pap1 + 4*pbp2 + sman*(-2 + x1)) + m2**2*sman*x1 - m4**2*sman*x1 - 2*pbp2*sman*x1 + m1**2*(2*(m3**2 + pap1 - pbp1) - sman*x1) + sman*(2*m2**2 - 2*(m4**2 + pap1 - pbp1 + sman) + sman*x1)*x2 + 
     &                  sman**2*x2**2)) - F1L*m3*(H3Lc*m1*(m3**4 - 2*m3**2*m4**2 + m4**4 + 2*m3**2*pap1 - 2*m4**2*pap1 + 2*m3**2*pap2 - 2*m4**2*pap2 - 4*m3**2*pbp2 + 4*m4**2*pbp2 - 8*pap1*pbp2 - 8*pap2*pbp2 + m3**2*sman - 
     &                  m4**2*sman + 2*pbp1*sman + 2*pbp2*sman - 2*sman**2 + m1**2*(-m3**2 + m4**2 - 2*pap1 - 2*pap2 + 2*(pbp1 + pbp2) + sman) - m3**2*sman*x1 + m4**2*sman*x1 + 2*pbp2*sman*x1 + sman**2*x1 + 
     &                  m2**2*(m3**2 - m4**2 + 2*pap1 + 2*pap2 - 2*pbp1 - 2*pbp2 + sman*(3 - 2*x1 - 2*x2)) + sman*(2*(pap1 + pap2 - pbp1 + sman) - sman*x1)*x2 - sman**2*x2**2) + 
     &               H3Rc*mq3c*(2*m3**2*pap1 - 2*m4**2*pap1 + 2*m3**2*pap2 - 2*m4**2*pap2 - 8*pap1*pbp2 - 8*pap2*pbp2 + 2*pbp1*sman + 2*pbp2*sman - 2*sman**2 + sman**2*x1 + sman*(2*(pap1 + pap2) + sman)*x2 + 
     &                  m2**2*(2*pap1 + 2*pap2 - 2*pbp1 - 2*pbp2 - sman*(-2 + x1 + x2)) - m1**2*(2*pap1 + 2*pap2 - 2*pbp1 - 2*pbp2 + sman*(-2 + x1 + x2))))) + 
     &         A3Lc*(F1R*m3*mq1*(2*H3Lc*m1*mq3c*(m1**2 - m2**2 - m3**2 + m4**2 - 2*pbp1 + 2*pbp2 - sman*x2) + 
     &               H3Rc*(2*m1**4 + 2*m3**2*pap2 - 2*m4**2*pap2 + 4*pap2*pbp1 + 4*pap1*pbp2 + m3**2*sman - m4**2*sman - 2*pbp2*sman + sman**2 + m2**2*(-m3**2 + m4**2 + 2*pap2 + 2*pbp2 + sman*(-1 + x1)) - sman**2*x1 + 
     &                  2*pap2*sman*x2 - m1**2*(2*m2**2 + 3*m3**2 - 3*m4**2 + 2*pap2 + 4*pbp1 - 2*pbp2 + sman - sman*x1 + 2*sman*x2))) + 
     &            F1L*(H3Lc*m1*mq3c*(-m3**4 + 2*m3**2*m4**2 - m4**4 - 4*m3**2*pbp1 + 4*m4**2*pbp1 + m3**2*sman - m4**2*sman - m3**2*sman*x1 + m4**2*sman*x1 + 2*pbp2*sman*x1 + 
     &                  sman*(-2*m3**2 + 2*m4**2 - 2*pbp1 + sman - sman*x1)*x2 - sman**2*x2**2 + m2**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + m1**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2))) + 
     &               H3Rc*(-2*m3**2*m4**2*pap1 + 2*m4**4*pap1 - 2*m3**2*m4**2*pap2 + 2*m4**4*pap2 + 4*m3**2*pap2*pbp1 - 4*m4**2*pap2*pbp1 + 4*m3**2*pap1*pbp2 + 4*m4**2*pap1*pbp2 + 8*m3**2*pap2*pbp2 + m3**4*sman - 
     &                  2*m3**2*m4**2*sman + m4**4*sman - 2*m3**2*pap2*sman + 2*m4**2*pap2*sman - 2*m4**2*pbp1*sman - 2*m3**2*pbp2*sman + m3**2*sman**2 + m4**2*sman**2 + 2*m3**2*pap2*sman*x1 - 2*m4**2*pap2*sman*x1 - 
     &                  4*pap2*pbp2*sman*x1 - m3**2*sman**2*x1 - sman*(m4**2*(2*pap1 + 4*pap2 + sman) - 2*pap2*(m3**2 + 2*pbp1 + sman*(-1 + x1)))*x2 + 2*pap2*sman**2*x2**2 + m1**4*(m3**2 - m4**2 + sman*(-1 + x1 + x2)) - 
     &                  m2**2*(m3**4 + m4**4 + m3**2*(-2*m4**2 - 2*pbp2 + sman - sman*x1) - 2*pap2*sman*(-1 + x1 + x2) + m4**2*(2*pap1 + 2*pap2 - 2*pbp1 + sman - sman*x2)) + 
     &                  m1**2*(-2*m3**4 - 2*m4**4 + 2*m4**2*pap1 + 2*m4**2*pap2 + 2*m4**2*pbp1 - 2*m4**2*sman + 2*pap2*sman + m4**2*sman*x1 - 2*pap2*sman*x1 + 2*pbp2*sman*x1 + 
     &                     sman*(3*m4**2 - 2*pap2 - 2*pbp1 + sman - sman*x1)*x2 - sman**2*x2**2 + m2**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + m3**2*(4*m4**2 - 2*(2*pbp1 + pbp2 + sman*x2)))))))))/
     &  ((m1**2 + m4**2 - mq3c**2 - 2*pap2)*(m1**2 - m2**2 + 2*(pbp1 + pbp2) - sman)*(-m3**2 + m4**2 + sman - sman*x2)*(-mq1**2 + sman*(-1 + x1 + x2))))


        radSpUi = dreal((g0*g1*(A1L*(A3Lc*(F1R*m3*(H3Lc*mq3c*(-4*mq1**2*pap1 - 4*mq1**2*pap2 - 2*mq1**2*sman + 2*pap1*sman + 2*pap2*sman - 2*pbp1*sman - 2*pbp2*sman + 2*sman**2 - sman**2*x1 - sman**2*x2 + 
     &                  m1**2*(6*mq1**2 - 2*pap1 - 2*pap2 - 10*pbp1 - 10*pbp2 - sman*(-2 + x1 + x2)) + m2**2*(2*mq1**2 - 6*pap1 - 6*pap2 + 2*pbp1 + 2*pbp2 + sman*(-2 + x1 + x2))) + 
     &               H3Rc*m1*(-2*m3**2*mq1**2 + 2*m4**2*mq1**2 - 4*mq1**2*pap1 - 4*mq1**2*pap2 + 8*pap2*pbp1 + 8*mq1**2*pbp2 - 8*pap1*pbp2 - 4*mq1**2*sman + 2*pap1*sman - 2*pap2*sman - 2*pbp1*sman - 6*pbp2*sman + 2*sman**2 + 
     &                  2*mq1**2*sman*x1 + 2*pap2*sman*x1 + 2*pbp2*sman*x1 - sman**2*x1 + m1**2*(2*m3**2 - 2*m4**2 + 6*mq1**2 - 2*pap1 - 2*pap2 - 10*pbp1 - 10*pbp2 + sman*(4 - 3*x1 - 3*x2)) + 
     &                  (2*(pap1 + 2*pap2 + pbp1 + 2*pbp2) - sman)*sman*x2 + m2**2*(-2*m3**2 + 2*m4**2 + 2*mq1**2 - 6*pap1 - 6*pap2 + 2*pbp1 + 2*pbp2 + sman*(-4 + 3*x1 + 3*x2)))) + 
     &            F1L*mq1*(H3Lc*mq3c*(12*pap2*pbp1 - 12*pap1*pbp2 + 2*(-m3**2 + m4**2 + 3*pap1)*sman - 2*pbp1*sman + 2*sman**2 - 4*pap1*sman*x1 - 2*pap2*sman*x1 - sman**2*x1 - 2*sman*(pap1 + sman)*x2 + 
     &                  m2**2*(2*m3**2 - 2*m4**2 - 6*pap1 + 2*pbp1 + sman*(-2 + x1 + 2*x2)) + m1**2*(2*m3**2 - 2*m4**2 - 2*pap1 - 10*pbp1 + sman*(-2 + 3*x1 + 2*x2))) + 
     &               H3Rc*m1*(-4*m4**2*pap1 - 4*m4**2*pap2 - 12*m4**2*pbp1 + 20*pap2*pbp1 - 12*m4**2*pbp2 - 20*pap1*pbp2 + 4*m4**2*sman + 6*pap1*sman - 4*pap2*sman - 2*pbp1*sman - 12*pbp2*sman + 4*sman**2 - 2*m3**2*sman*x1 - 
     &                  4*pap1*sman*x1 + 10*pbp2*sman*x1 - 5*sman**2*x1 + 2*sman**2*x1**2 + m1**2*(4*m3**2 - 2*(pap1 + 5*pbp1) + sman*x1) + 2*sman*(-m3**2 + 2*pap2 + pbp1 + 6*pbp2 + sman*(-2 + x1))*x2 + 
     &                  m2**2*(-4*m4**2 - 6*pap1 + 2*pbp1 + sman*(-4 + 3*x1 + 4*x2))))) + 
     &         A3Rc*(F1L*m3*mq1*(2*H3Rc*m1*mq3c*(2*m1**2 - 2*m2**2 - 2*pap1 - 2*pap2 - 6*pbp1 - 6*pbp2 - sman*(-2 + x1 + x2)) + 
     &               H3Lc*(4*m1**4 + 2*(-6*pap1*pbp2 + pap2*(6*pbp1 + sman*(-3 + x1)) + sman*(-m3**2 + m4**2 + pbp2 + sman*(-1 + x1))) + sman*(2*pap1 + 4*pap2 + sman)*x2 + 
     &                  m1**2*(-4*m2**2 + 2*m3**2 - 2*(m4**2 + 2*pap1 + pap2 + 6*pbp1 + pbp2 - 3*sman + 2*sman*x1) - 5*sman*x2) + m2**2*(2*m3**2 - 2*(m4**2 - 3*pap2 + pbp2 + sman*(-1 + x1)) - sman*x2))) + 
     &            F1R*(2*H3Rc*m1*mq3c*(-(m3**2*mq1**2) + m4**2*mq1**2 - 4*mq1**2*pbp1 + 4*pap2*pbp1 - 4*pap1*pbp2 + mq1**2*sman + 2*pap1*sman + 2*pbp1*sman - 2*pap1*sman*x1 - pap2*sman*x1 - 2*pbp1*sman*x1 - pbp2*sman*x1 - 
     &                  (mq1**2 + pap1 + pbp1)*sman*x2 + m2**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + m1**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2))) + 
     &               H3Lc*(-4*m4**2*mq1**2*pap1 - 4*m4**2*mq1**2*pap2 + 8*mq1**2*pap2*pbp1 - 8*mq1**2*pap1*pbp2 - m3**2*mq1**2*sman - m4**2*mq1**2*sman + 2*m4**2*pap1*sman + 2*m4**2*pap2*sman - 4*mq1**2*pap2*sman - 
     &                  2*m4**2*pbp1*sman - 4*pap2*pbp1*sman - 2*m4**2*pbp2*sman + 4*pap1*pbp2*sman + m3**2*sman**2 + m4**2*sman**2 - mq1**2*sman**2 + 2*pap2*sman**2 - 2*pbp2*sman**2 + sman**3 + 2*mq1**2*pap2*sman*x1 + 
     &                  4*pap2*pbp1*sman*x1 - 4*pap1*pbp2*sman*x1 - m3**2*sman**2*x1 + mq1**2*sman**2*x1 - 2*pap2*sman**2*x1 + 2*pbp2*sman**2*x1 - 2*sman**3*x1 + sman**3*x1**2 + 
     &                  sman*(-4*pap1*pbp2 + pap2*(4*pbp1 - 2*sman) - m3**2*sman + mq1**2*(2*pap1 + 4*pap2 + sman) + sman*(2*pbp2 + sman*(-1 + x1)))*x2 + 2*m1**4*(m3**2 - m4**2 + sman*(-1 + x1 + x2)) + 
     &                  m1**2*(5*m4**2*mq1**2 - 2*m4**2*pap1 - 2*m4**2*pap2 - 10*m4**2*pbp1 - 8*mq1**2*pbp1 + 8*pap2*pbp1 - 10*m4**2*pbp2 - 8*pap1*pbp2 + m4**2*sman + 5*mq1**2*sman + 4*pap1*sman - 2*pap2*sman + 
     &                     4*pbp1*sman - 10*pbp2*sman + sman**2 - 3*mq1**2*sman*x1 - 4*pap1*sman*x1 - 4*pbp1*sman*x1 + 8*pbp2*sman*x1 - 2*sman**2*x1 + sman**2*x1**2 + 
     &                     sman*(-5*mq1**2 - 2*pap1 + 2*pap2 - 2*pbp1 + 10*pbp2 + sman*(-1 + x1))*x2 + m3**2*(mq1**2 - sman*(-1 + x1 + x2)) - 2*m2**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2))) + 
     &                  m2**2*(m4**2*(mq1**2 - 6*pap1 - 6*pap2 + 2*(pbp1 + pbp2) - sman) - sman*(mq1**2 - 6*pap2 + 2*pbp2 + sman*(-1 + x1))*(-1 + x1 + x2) + m3**2*(mq1**2 + sman*(-1 + x1 + x2))))))) + 
     &      A1R*(A3Rc*(F1L*m3*(H3Rc*mq3c*(-4*mq1**2*pap1 - 4*mq1**2*pap2 - 2*mq1**2*sman + 2*pap1*sman + 2*pap2*sman - 2*pbp1*sman - 2*pbp2*sman + 2*sman**2 - sman**2*x1 - sman**2*x2 + 
     &                  m1**2*(6*mq1**2 - 2*pap1 - 2*pap2 - 10*pbp1 - 10*pbp2 - sman*(-2 + x1 + x2)) + m2**2*(2*mq1**2 - 6*pap1 - 6*pap2 + 2*pbp1 + 2*pbp2 + sman*(-2 + x1 + x2))) + 
     &               H3Lc*m1*(-2*m3**2*mq1**2 + 2*m4**2*mq1**2 - 4*mq1**2*pap1 - 4*mq1**2*pap2 + 8*pap2*pbp1 + 8*mq1**2*pbp2 - 8*pap1*pbp2 - 4*mq1**2*sman + 2*pap1*sman - 2*pap2*sman - 2*pbp1*sman - 6*pbp2*sman + 2*sman**2 + 
     &                  2*mq1**2*sman*x1 + 2*pap2*sman*x1 + 2*pbp2*sman*x1 - sman**2*x1 + m1**2*(2*m3**2 - 2*m4**2 + 6*mq1**2 - 2*pap1 - 2*pap2 - 10*pbp1 - 10*pbp2 + sman*(4 - 3*x1 - 3*x2)) + 
     &                  (2*(pap1 + 2*pap2 + pbp1 + 2*pbp2) - sman)*sman*x2 + m2**2*(-2*m3**2 + 2*m4**2 + 2*mq1**2 - 6*pap1 - 6*pap2 + 2*pbp1 + 2*pbp2 + sman*(-4 + 3*x1 + 3*x2)))) + 
     &            F1R*mq1*(H3Rc*mq3c*(12*pap2*pbp1 - 12*pap1*pbp2 + 2*(-m3**2 + m4**2 + 3*pap1)*sman - 2*pbp1*sman + 2*sman**2 - 4*pap1*sman*x1 - 2*pap2*sman*x1 - sman**2*x1 - 2*sman*(pap1 + sman)*x2 + 
     &                  m2**2*(2*m3**2 - 2*m4**2 - 6*pap1 + 2*pbp1 + sman*(-2 + x1 + 2*x2)) + m1**2*(2*m3**2 - 2*m4**2 - 2*pap1 - 10*pbp1 + sman*(-2 + 3*x1 + 2*x2))) + 
     &               H3Lc*m1*(-4*m4**2*pap1 - 4*m4**2*pap2 - 12*m4**2*pbp1 + 20*pap2*pbp1 - 12*m4**2*pbp2 - 20*pap1*pbp2 + 4*m4**2*sman + 6*pap1*sman - 4*pap2*sman - 2*pbp1*sman - 12*pbp2*sman + 4*sman**2 - 2*m3**2*sman*x1 - 
     &                  4*pap1*sman*x1 + 10*pbp2*sman*x1 - 5*sman**2*x1 + 2*sman**2*x1**2 + m1**2*(4*m3**2 - 2*(pap1 + 5*pbp1) + sman*x1) + 2*sman*(-m3**2 + 2*pap2 + pbp1 + 6*pbp2 + sman*(-2 + x1))*x2 + 
     &                  m2**2*(-4*m4**2 - 6*pap1 + 2*pbp1 + sman*(-4 + 3*x1 + 4*x2))))) + 
     &         A3Lc*(F1R*m3*mq1*(2*H3Lc*m1*mq3c*(2*m1**2 - 2*m2**2 - 2*pap1 - 2*pap2 - 6*pbp1 - 6*pbp2 - sman*(-2 + x1 + x2)) + 
     &               H3Rc*(4*m1**4 + 2*(-6*pap1*pbp2 + pap2*(6*pbp1 + sman*(-3 + x1)) + sman*(-m3**2 + m4**2 + pbp2 + sman*(-1 + x1))) + sman*(2*pap1 + 4*pap2 + sman)*x2 + 
     &                  m1**2*(-4*m2**2 + 2*m3**2 - 2*(m4**2 + 2*pap1 + pap2 + 6*pbp1 + pbp2 - 3*sman + 2*sman*x1) - 5*sman*x2) + m2**2*(2*m3**2 - 2*(m4**2 - 3*pap2 + pbp2 + sman*(-1 + x1)) - sman*x2))) + 
     &            F1L*(2*H3Lc*m1*mq3c*(-(m3**2*mq1**2) + m4**2*mq1**2 - 4*mq1**2*pbp1 + 4*pap2*pbp1 - 4*pap1*pbp2 + mq1**2*sman + 2*pap1*sman + 2*pbp1*sman - 2*pap1*sman*x1 - pap2*sman*x1 - 2*pbp1*sman*x1 - pbp2*sman*x1 - 
     &                  (mq1**2 + pap1 + pbp1)*sman*x2 + m2**2*(-m3**2 + m4**2 - sman*(-1 + x1 + x2)) + m1**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2))) + 
     &               H3Rc*(-4*m4**2*mq1**2*pap1 - 4*m4**2*mq1**2*pap2 + 8*mq1**2*pap2*pbp1 - 8*mq1**2*pap1*pbp2 - m3**2*mq1**2*sman - m4**2*mq1**2*sman + 2*m4**2*pap1*sman + 2*m4**2*pap2*sman - 4*mq1**2*pap2*sman - 
     &                  2*m4**2*pbp1*sman - 4*pap2*pbp1*sman - 2*m4**2*pbp2*sman + 4*pap1*pbp2*sman + m3**2*sman**2 + m4**2*sman**2 - mq1**2*sman**2 + 2*pap2*sman**2 - 2*pbp2*sman**2 + sman**3 + 2*mq1**2*pap2*sman*x1 + 
     &                  4*pap2*pbp1*sman*x1 - 4*pap1*pbp2*sman*x1 - m3**2*sman**2*x1 + mq1**2*sman**2*x1 - 2*pap2*sman**2*x1 + 2*pbp2*sman**2*x1 - 2*sman**3*x1 + sman**3*x1**2 + 
     &                  sman*(-4*pap1*pbp2 + pap2*(4*pbp1 - 2*sman) - m3**2*sman + mq1**2*(2*pap1 + 4*pap2 + sman) + sman*(2*pbp2 + sman*(-1 + x1)))*x2 + 2*m1**4*(m3**2 - m4**2 + sman*(-1 + x1 + x2)) + 
     &                  m1**2*(5*m4**2*mq1**2 - 2*m4**2*pap1 - 2*m4**2*pap2 - 10*m4**2*pbp1 - 8*mq1**2*pbp1 + 8*pap2*pbp1 - 10*m4**2*pbp2 - 8*pap1*pbp2 + m4**2*sman + 5*mq1**2*sman + 4*pap1*sman - 2*pap2*sman + 
     &                     4*pbp1*sman - 10*pbp2*sman + sman**2 - 3*mq1**2*sman*x1 - 4*pap1*sman*x1 - 4*pbp1*sman*x1 + 8*pbp2*sman*x1 - 2*sman**2*x1 + sman**2*x1**2 + 
     &                     sman*(-5*mq1**2 - 2*pap1 + 2*pap2 - 2*pbp1 + 10*pbp2 + sman*(-1 + x1))*x2 + m3**2*(mq1**2 - sman*(-1 + x1 + x2)) - 2*m2**2*(m3**2 - m4**2 + sman*(-1 + x1 + x2))) + 
     &                  m2**2*(m4**2*(mq1**2 - 6*pap1 - 6*pap2 + 2*(pbp1 + pbp2) - sman) - sman*(mq1**2 - 6*pap2 + 2*pbp2 + sman*(-1 + x1))*(-1 + x1 + x2) + m3**2*(mq1**2 + sman*(-1 + x1 + x2)))))))))/
     &  (2.*(m1**2 + m4**2 - mq3c**2 - 2*pap2)*(mq1**2 - sman)*(m1**2 - m2**2 + 2*(pbp1 + pbp2) - sman)*(mq1**2 - sman*(-1 + x1 + x2))))


        radSpUf = dreal((g1**2*(A1L*(A3Lc*(F1R*m3*(H3Rc*m1*(2*m3**2*mq1**2 - 6*m4**2*mq1**2 + 4*m4**2*pap1 + 4*mq1**2*pap1 + 4*m3**2*pap2 + 8*mq1**2*pap2 + 4*m4**2*pbp1 - 4*pap2*pbp1 + 4*m3**2*pbp2 + 4*pap1*pbp2 - 3*m3**2*sman + 
     &                  3*m4**2*sman - 2*mq1**2*sman - 2*pap1*sman - 6*pap2*sman + 2*pbp1*sman - 2*pbp2*sman + sman**2 - 2*mq1**2*sman*x1 + 2*pap2*sman*x1 + 2*pbp2*sman*x1 + 
     &                  m1**2*(-m3**2 + m4**2 - 4*mq1**2 + 2*pap1 + 2*pap2 + 2*pbp1 + 2*pbp2 + sman*(3 - 2*x1 - 2*x2)) + 2*(mq1**2 - pap1 - pbp1)*sman*x2 + 
     &                  m2**2*(m3**2 - m4**2 + 4*mq1**2 - 2*pap1 - 2*pap2 - 2*pbp1 - 2*pbp2 + sman*(-3 + 2*x1 + 2*x2))) + 
     &               H3Lc*mq3c*(4*mq1**2*pap1 + 8*mq1**2*pap2 - 4*pap2*pbp1 + 4*pap1*pbp2 - m3**2*sman + m4**2*sman - 4*mq1**2*sman - 2*pap1*sman - 2*pap2*sman + 2*pbp1*sman + 2*pbp2*sman + 3*sman**2 - 2*sman**2*x1 + 
     &                  m1**2*(-m3**2 + m4**2 - 4*mq1**2 + 2*pap1 + 2*pap2 + 2*pbp1 + 2*pbp2 + sman*(3 - 2*x1 - 2*x2)) - 2*sman**2*x2 + 
     &                  m2**2*(m3**2 - m4**2 + 4*mq1**2 - 2*pap1 - 2*pap2 - 2*pbp1 - 2*pbp2 + sman*(-3 + 2*x1 + 2*x2)))) + 
     &            F1L*mq1*(H3Lc*mq3c*(4*m3**2*pap1 + 8*m3**2*pap2 - 4*m4**2*pap2 - 4*pap2*pbp1 + 4*pap1*pbp2 - 5*m3**2*sman + 3*m4**2*sman - 2*pap1*sman - 4*pap2*sman + 2*pbp1*sman + 3*sman**2 + 4*pap2*sman*x1 - 
     &                  2*sman**2*x1 + m1**2*(-5*m3**2 + 3*m4**2 + 2*pap1 + 2*pbp1 + sman*(3 - 2*x1 - 3*x2)) + (4*pap2 - 3*sman)*sman*x2 + m2**2*(5*m3**2 - 3*m4**2 - 2*pap1 - 2*pbp1 + sman*(-3 + 2*x1 + 3*x2))) + 
     &               H3Rc*m1*(2*m3**4 - 6*m3**2*m4**2 + 4*m4**4 + 4*m3**2*pap1 + 4*m4**2*pap1 + 16*m3**2*pap2 - 4*m4**2*pap2 + 4*m4**2*pbp1 - 4*pap2*pbp1 + 8*m3**2*pbp2 + 4*pap1*pbp2 - 5*m3**2*sman + 5*m4**2*sman - 
     &                  2*pap1*sman - 4*pap2*sman + 2*pbp1*sman + sman**2 - 2*m3**2*sman*x1 - 4*m4**2*sman*x1 + 2*pap2*sman*x1 - 2*pbp2*sman*x1 + m1**2*(-5*m3**2 + 3*m4**2 + 2*pap1 + 2*pbp1 + sman*(3 - 2*x1 - 3*x2)) + 
     &                  sman*(-6*m4**2 - 2*pap1 + 4*pap2 - 2*pbp1 - 3*sman + 2*sman*x1)*x2 + 2*sman**2*x2**2 + m2**2*(5*m3**2 - 3*m4**2 - 2*pap1 - 2*pbp1 + sman*(-3 + 2*x1 + 3*x2))))) + 
     &         A3Rc*(-2*F1L*m3*mq1*(H3Rc*m1*mq3c*(m3**2 - m4**2 - 2*pap1 - 4*pap2 - 2*pbp1 - 4*pbp2 + sman + 2*sman*x1 + 2*sman*x2) + 
     &               H3Lc*(-(m2**2*m4**2) + 2*pap2*(-m3**2 + pbp1) - 2*pap1*pbp2 - 2*m2**2*(pap2 + pbp2) + 2*pbp2*sman + m4**2*(-2*pap1 + sman) - pap2*sman*x1 + (m2**2 + pap1 - sman)*sman*x2 + 
     &                  m1**2*(m3**2 - 2*pap1 - 2*pap2 - 2*pbp1 - 2*pbp2 + sman + 2*sman*x1 + sman*x2))) + 
     &            F1R*(2*H3Rc*m1*mq3c*(m4**2*(mq1**2 - 2*(pap2 + pbp2 - sman)) + m3**2*(-3*mq1**2 + 2*(pap1 + 2*pap2 + pbp1 + 2*pbp2 - sman)) + sman*(-(mq1**2*(-1 + x2)) + 2*(pap2 + pbp2 - sman)*(-1 + x1 + x2))) + 
     &               H3Lc*(4*m4**2*mq1**2*pap1 + 8*m3**2*mq1**2*pap2 - 4*m3**2*pap2*pbp1 + 4*m3**2*pap1*pbp2 + m3**4*sman - 3*m3**2*m4**2*sman + 2*m4**4*sman - m3**2*mq1**2*sman - m4**2*mq1**2*sman + 2*m3**2*pap2*sman - 
     &                  2*m4**2*pap2*sman - 2*m3**2*pbp2*sman + 2*m4**2*pbp2*sman + m3**2*sman**2 + 2*m4**2*sman**2 - mq1**2*sman**2 - 2*pap2*sman**2 + 2*pbp2*sman**2 - 2*mq1**2*pap2*sman*x1 - m3**2*sman**2*x1 - 
     &                  2*m4**2*sman**2*x1 + mq1**2*sman**2*x1 + 2*pap2*sman**2*x1 - 2*pbp2*sman**2*x1 + sman*(mq1**2*(-2*pap1 + sman) + sman*(-3*m4**2 + 2*pap2 - 2*pbp2 + sman*(-1 + x1)))*x2 + sman**3*x2**2 + 
     &                  m2**2*(-m3**4 + m3**2*(3*m4**2 + mq1**2 + 2*pap2 + 2*pbp2 + sman*(-1 + x1)) - (2*m4**2 - mq1**2 + 2*(pap2 + pbp2) - sman*x2)*(m4**2 - sman*(-1 + x1 + x2))) + 
     &                  m1**2*(m3**4 + 2*m4**4 + m3**2*(-3*m4**2 - 7*mq1**2 + 4*pap1 + 6*pap2 + 4*pbp1 + 6*pbp2 - sman*(3 + x1)) + m4**2*(mq1**2 - 2*(pap2 + pbp2 + sman*(-3 + x1)) - 3*sman*x2) + 
     &                     sman*(mq1**2*(1 + x1 - x2) + (-1 + x1 + x2)*(2*(pap2 + pbp2 - 2*sman) + sman*x2))))))) + 
     &      A1R*(A3Rc*(F1L*m3*(H3Lc*m1*(2*m3**2*mq1**2 - 6*m4**2*mq1**2 + 4*m4**2*pap1 + 4*mq1**2*pap1 + 4*m3**2*pap2 + 8*mq1**2*pap2 + 4*m4**2*pbp1 - 4*pap2*pbp1 + 4*m3**2*pbp2 + 4*pap1*pbp2 - 3*m3**2*sman + 3*m4**2*sman - 
     &                  2*mq1**2*sman - 2*pap1*sman - 6*pap2*sman + 2*pbp1*sman - 2*pbp2*sman + sman**2 - 2*mq1**2*sman*x1 + 2*pap2*sman*x1 + 2*pbp2*sman*x1 + 
     &                  m1**2*(-m3**2 + m4**2 - 4*mq1**2 + 2*pap1 + 2*pap2 + 2*pbp1 + 2*pbp2 + sman*(3 - 2*x1 - 2*x2)) + 2*(mq1**2 - pap1 - pbp1)*sman*x2 + 
     &                  m2**2*(m3**2 - m4**2 + 4*mq1**2 - 2*pap1 - 2*pap2 - 2*pbp1 - 2*pbp2 + sman*(-3 + 2*x1 + 2*x2))) + 
     &               H3Rc*mq3c*(4*mq1**2*pap1 + 8*mq1**2*pap2 - 4*pap2*pbp1 + 4*pap1*pbp2 - m3**2*sman + m4**2*sman - 4*mq1**2*sman - 2*pap1*sman - 2*pap2*sman + 2*pbp1*sman + 2*pbp2*sman + 3*sman**2 - 2*sman**2*x1 + 
     &                  m1**2*(-m3**2 + m4**2 - 4*mq1**2 + 2*pap1 + 2*pap2 + 2*pbp1 + 2*pbp2 + sman*(3 - 2*x1 - 2*x2)) - 2*sman**2*x2 + 
     &                  m2**2*(m3**2 - m4**2 + 4*mq1**2 - 2*pap1 - 2*pap2 - 2*pbp1 - 2*pbp2 + sman*(-3 + 2*x1 + 2*x2)))) + 
     &            F1R*mq1*(H3Rc*mq3c*(4*m3**2*pap1 + 8*m3**2*pap2 - 4*m4**2*pap2 - 4*pap2*pbp1 + 4*pap1*pbp2 - 5*m3**2*sman + 3*m4**2*sman - 2*pap1*sman - 4*pap2*sman + 2*pbp1*sman + 3*sman**2 + 4*pap2*sman*x1 - 
     &                  2*sman**2*x1 + m1**2*(-5*m3**2 + 3*m4**2 + 2*pap1 + 2*pbp1 + sman*(3 - 2*x1 - 3*x2)) + (4*pap2 - 3*sman)*sman*x2 + m2**2*(5*m3**2 - 3*m4**2 - 2*pap1 - 2*pbp1 + sman*(-3 + 2*x1 + 3*x2))) + 
     &               H3Lc*m1*(2*m3**4 - 6*m3**2*m4**2 + 4*m4**4 + 4*m3**2*pap1 + 4*m4**2*pap1 + 16*m3**2*pap2 - 4*m4**2*pap2 + 4*m4**2*pbp1 - 4*pap2*pbp1 + 8*m3**2*pbp2 + 4*pap1*pbp2 - 5*m3**2*sman + 5*m4**2*sman - 
     &                  2*pap1*sman - 4*pap2*sman + 2*pbp1*sman + sman**2 - 2*m3**2*sman*x1 - 4*m4**2*sman*x1 + 2*pap2*sman*x1 - 2*pbp2*sman*x1 + m1**2*(-5*m3**2 + 3*m4**2 + 2*pap1 + 2*pbp1 + sman*(3 - 2*x1 - 3*x2)) + 
     &                  sman*(-6*m4**2 - 2*pap1 + 4*pap2 - 2*pbp1 - 3*sman + 2*sman*x1)*x2 + 2*sman**2*x2**2 + m2**2*(5*m3**2 - 3*m4**2 - 2*pap1 - 2*pbp1 + sman*(-3 + 2*x1 + 3*x2))))) + 
     &         A3Lc*(-2*F1R*m3*mq1*(H3Lc*m1*mq3c*(m3**2 - m4**2 - 2*pap1 - 4*pap2 - 2*pbp1 - 4*pbp2 + sman + 2*sman*x1 + 2*sman*x2) + 
     &               H3Rc*(-(m2**2*m4**2) + 2*pap2*(-m3**2 + pbp1) - 2*pap1*pbp2 - 2*m2**2*(pap2 + pbp2) + 2*pbp2*sman + m4**2*(-2*pap1 + sman) - pap2*sman*x1 + (m2**2 + pap1 - sman)*sman*x2 + 
     &                  m1**2*(m3**2 - 2*pap1 - 2*pap2 - 2*pbp1 - 2*pbp2 + sman + 2*sman*x1 + sman*x2))) + 
     &            F1L*(2*H3Lc*m1*mq3c*(m4**2*(mq1**2 - 2*(pap2 + pbp2 - sman)) + m3**2*(-3*mq1**2 + 2*(pap1 + 2*pap2 + pbp1 + 2*pbp2 - sman)) + sman*(-(mq1**2*(-1 + x2)) + 2*(pap2 + pbp2 - sman)*(-1 + x1 + x2))) + 
     &               H3Rc*(4*m4**2*mq1**2*pap1 + 8*m3**2*mq1**2*pap2 - 4*m3**2*pap2*pbp1 + 4*m3**2*pap1*pbp2 + m3**4*sman - 3*m3**2*m4**2*sman + 2*m4**4*sman - m3**2*mq1**2*sman - m4**2*mq1**2*sman + 2*m3**2*pap2*sman - 
     &                  2*m4**2*pap2*sman - 2*m3**2*pbp2*sman + 2*m4**2*pbp2*sman + m3**2*sman**2 + 2*m4**2*sman**2 - mq1**2*sman**2 - 2*pap2*sman**2 + 2*pbp2*sman**2 - 2*mq1**2*pap2*sman*x1 - m3**2*sman**2*x1 - 
     &                  2*m4**2*sman**2*x1 + mq1**2*sman**2*x1 + 2*pap2*sman**2*x1 - 2*pbp2*sman**2*x1 + sman*(mq1**2*(-2*pap1 + sman) + sman*(-3*m4**2 + 2*pap2 - 2*pbp2 + sman*(-1 + x1)))*x2 + sman**3*x2**2 + 
     &                  m2**2*(-m3**4 + m3**2*(3*m4**2 + mq1**2 + 2*pap2 + 2*pbp2 + sman*(-1 + x1)) - (2*m4**2 - mq1**2 + 2*(pap2 + pbp2) - sman*x2)*(m4**2 - sman*(-1 + x1 + x2))) + 
     &                  m1**2*(m3**4 + 2*m4**4 + m3**2*(-3*m4**2 - 7*mq1**2 + 4*pap1 + 6*pap2 + 4*pbp1 + 6*pbp2 - sman*(3 + x1)) + m4**2*(mq1**2 - 2*(pap2 + pbp2 + sman*(-3 + x1)) - 3*sman*x2) + 
     &                     sman*(mq1**2*(1 + x1 - x2) + (-1 + x1 + x2)*(2*(pap2 + pbp2 - 2*sman) + sman*x2)))))))))/
     &  ((m1**2 + m4**2 - mq3c**2 - 2*pap2)*(mq1**2 - sman)*(-m3**2 + m4**2 + sman - sman*x2)*(mq1**2 - sman*(-1 + x1 + x2))))


        radSfUi = dreal((g0*g1*(A1L*(-2*A3Lc*(F1L*mq1*(H3Lc*mq3c*(-2*m3**2*pap1 + 2*m4**2*pap1 + 4*pap2*pbp1 + 4*pap1*pbp2 - m3**2*sman + m4**2*sman - 2*pbp1*sman + sman**2 + m1**2*(m3**2 - m4**2 + 2*pap1 - 2*pbp1 + sman*(-1 + x2)) + 
     &                  m2**2*(m3**2 - m4**2 - 2*pap1 + 2*pbp1 + sman*(-1 + x2)) - sman*(2*pap1 + sman)*x2) + 
     &               H3Rc*m1*(-m3**4 + m4**4 + 2*m4**2*pap1 - 4*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 2*m4**2*sman - 2*pbp1*sman + sman**2 - 2*m2**2*(m4**2 + pap1 - pbp1 + sman) + 
     &                  m3**2*(-2*pap1 + 4*pbp2 + sman*(-2 + x1)) + m2**2*sman*x1 - m4**2*sman*x1 - 2*pbp2*sman*x1 + m1**2*(2*(m3**2 + pap1 - pbp1) - sman*x1) + sman*(2*m2**2 - 2*(m4**2 + pap1 - pbp1 + sman) + sman*x1)*x2 + 
     &                  sman**2*x2**2)) - F1R*m3*(H3Lc*mq3c*(m1**2*(m3**2 - m4**2 + 4*pbp1 + sman*(-3 + x2)) + sman*(m3**2 - m4**2 - 4*pbp2 + sman + sman*x2) - 
     &                  m2**2*(m3**2 - m4**2 - 4*pap1 - 4*pap2 - 4*pbp2 + sman + sman*x2)) + 
     &               H3Rc*m1*(2*m4**2*pap2 - 4*pap2*pbp1 + 2*m4**2*pbp2 + 4*pap1*pbp2 + 8*pap2*pbp2 + 8*pbp2**2 + 2*m1**2*(pap2 + 2*pbp1 + pbp2 - 2*sman) - 2*m3**2*(pap2 + pbp2 - sman) - 2*m4**2*sman - 8*pbp2*sman + 
     &                  2*sman**2 + m1**2*sman*x1 - sman**2*x1 + sman*(m1**2 - 2*(pap2 + pbp2) + sman)*x2 + m2**2*(4*pap1 + 2*pap2 + 2*pbp2 - sman*(x1 + x2))))) + 
     &         A3Rc*(2*F1L*m3*mq1*(2*H3Rc*m1*mq3c*(-m1**2 + m2**2 + m3**2 - m4**2 + 2*pbp1 - 2*pbp2 + sman*x2) + 
     &               H3Lc*(-2*m1**4 - 2*m3**2*pap2 + 2*m4**2*pap2 - 4*pap2*pbp1 + 4*pap1*pbp2 + 8*pap2*pbp2 + m3**2*sman - m4**2*sman - 2*pbp2*sman + sman**2 + m2**2*(-m3**2 + m4**2 - 2*pap2 + 2*pbp2 + sman*(-1 + x1)) - 
     &                  sman**2*x1 - 2*pap2*sman*x2 + m1**2*(2*m2**2 + m3**2 - m4**2 + 2*pap2 + 4*pbp1 - 6*pbp2 + sman*(-1 + x1 + 2*x2)))) + 
     &            F1R*(-2*H3Rc*m1*mq3c*(-2*m3**2*pap1 + 2*m4**2*pap1 - 2*m3**2*pbp1 + 2*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 8*pbp1*pbp2 - m3**2*sman + m4**2*sman - 4*pbp1*sman + sman**2 + 
     &                  m1**2*(m3**2 - m4**2 + 2*pap1 + 2*pbp1 + sman*(-1 + x2)) - sman*(2*(pap1 + pbp1) + sman)*x2 + m2**2*(-m3**2 + m4**2 - 2*pap1 - 2*pbp1 + sman - sman*x2)) + 
     &               H3Lc*(-2*m1**4*(m3**2 - m4**2 + 2*pap1 + 2*pbp1 + sman*(-1 + x2)) + m1**2*
     &                   (m3**4 - m4**4 - 4*m4**2*pap1 + 4*m4**2*pbp1 + 4*m4**2*pbp2 - 16*pap1*pbp2 - 16*pbp1*pbp2 - 6*m4**2*sman + 8*pbp1*sman + 4*pbp2*sman - 5*sman**2 + m4**2*sman*x1 + 3*sman**2*x1 + 
     &                     m3**2*(4*(pap1 + pbp1 - pbp2) - sman*x1) + 2*m2**2*(m3**2 - m4**2 + 2*pap1 + 2*pbp1 + sman*(-1 + x2)) + sman*(2*m4**2 + 4*pap1 - 4*pbp2 + 6*sman - sman*x1)*x2 - sman**2*x2**2) + 
     &                  (m3**2 - m4**2 - 4*pbp2 + sman + sman*x2)*(-4*pap2*pbp1 + 4*pap1*pbp2 + sman*(m3**2 + m4**2 - sman*(-1 + x1 + x2))) + 
     &                  m2**2*(-m3**4 + m4**4 - 8*pap2*pbp1 + 8*pap1*pbp2 + 8*pap2*sman + 4*pbp2*sman - sman**2 + m3**2*(4*pbp2 + sman*(-2 + x1)) - 4*pap2*sman*x1 - 4*pbp2*sman*x1 + sman**2*x1 + 
     &                     sman*(-4*(pap1 + 2*pap2 + pbp2) + sman*x1)*x2 + sman**2*x2**2 + m4**2*(8*pap1 + 8*pap2 + 4*pbp2 - sman*(x1 + 2*x2))))))) + 
     &      A1R*(2*A3Rc*(-(F1R*mq1*(H3Rc*mq3c*(-2*m3**2*pap1 + 2*m4**2*pap1 + 4*pap2*pbp1 + 4*pap1*pbp2 - m3**2*sman + m4**2*sman - 2*pbp1*sman + sman**2 + m1**2*(m3**2 - m4**2 + 2*pap1 - 2*pbp1 + sman*(-1 + x2)) + 
     &                    m2**2*(m3**2 - m4**2 - 2*pap1 + 2*pbp1 + sman*(-1 + x2)) - sman*(2*pap1 + sman)*x2) + 
     &                 H3Lc*m1*(-m3**4 + m4**4 + 2*m4**2*pap1 - 4*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 2*m4**2*sman - 2*pbp1*sman + sman**2 - 2*m2**2*(m4**2 + pap1 - pbp1 + sman) + 
     &                    m3**2*(-2*pap1 + 4*pbp2 + sman*(-2 + x1)) + m2**2*sman*x1 - m4**2*sman*x1 - 2*pbp2*sman*x1 + m1**2*(2*(m3**2 + pap1 - pbp1) - sman*x1) + 
     &                    sman*(2*m2**2 - 2*(m4**2 + pap1 - pbp1 + sman) + sman*x1)*x2 + sman**2*x2**2))) + 
     &            F1L*m3*(H3Rc*mq3c*(m1**2*(m3**2 - m4**2 + 4*pbp1 + sman*(-3 + x2)) + sman*(m3**2 - m4**2 - 4*pbp2 + sman + sman*x2) - m2**2*(m3**2 - m4**2 - 4*pap1 - 4*pap2 - 4*pbp2 + sman + sman*x2)) + 
     &               H3Lc*m1*(2*m4**2*pap2 - 4*pap2*pbp1 + 2*m4**2*pbp2 + 4*pap1*pbp2 + 8*pap2*pbp2 + 8*pbp2**2 + 2*m1**2*(pap2 + 2*pbp1 + pbp2 - 2*sman) - 2*m3**2*(pap2 + pbp2 - sman) - 2*m4**2*sman - 8*pbp2*sman + 
     &                  2*sman**2 + m1**2*sman*x1 - sman**2*x1 + sman*(m1**2 - 2*(pap2 + pbp2) + sman)*x2 + m2**2*(4*pap1 + 2*pap2 + 2*pbp2 - sman*(x1 + x2))))) + 
     &         A3Lc*(2*F1R*m3*mq1*(2*H3Lc*m1*mq3c*(-m1**2 + m2**2 + m3**2 - m4**2 + 2*pbp1 - 2*pbp2 + sman*x2) + 
     &               H3Rc*(-2*m1**4 - 2*m3**2*pap2 + 2*m4**2*pap2 - 4*pap2*pbp1 + 4*pap1*pbp2 + 8*pap2*pbp2 + m3**2*sman - m4**2*sman - 2*pbp2*sman + sman**2 + m2**2*(-m3**2 + m4**2 - 2*pap2 + 2*pbp2 + sman*(-1 + x1)) - 
     &                  sman**2*x1 - 2*pap2*sman*x2 + m1**2*(2*m2**2 + m3**2 - m4**2 + 2*pap2 + 4*pbp1 - 6*pbp2 + sman*(-1 + x1 + 2*x2)))) + 
     &            F1L*(-2*H3Lc*m1*mq3c*(-2*m3**2*pap1 + 2*m4**2*pap1 - 2*m3**2*pbp1 + 2*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 8*pbp1*pbp2 - m3**2*sman + m4**2*sman - 4*pbp1*sman + sman**2 + 
     &                  m1**2*(m3**2 - m4**2 + 2*pap1 + 2*pbp1 + sman*(-1 + x2)) - sman*(2*(pap1 + pbp1) + sman)*x2 + m2**2*(-m3**2 + m4**2 - 2*pap1 - 2*pbp1 + sman - sman*x2)) + 
     &               H3Rc*(-2*m1**4*(m3**2 - m4**2 + 2*pap1 + 2*pbp1 + sman*(-1 + x2)) + m1**2*
     &                   (m3**4 - m4**4 - 4*m4**2*pap1 + 4*m4**2*pbp1 + 4*m4**2*pbp2 - 16*pap1*pbp2 - 16*pbp1*pbp2 - 6*m4**2*sman + 8*pbp1*sman + 4*pbp2*sman - 5*sman**2 + m4**2*sman*x1 + 3*sman**2*x1 + 
     &                     m3**2*(4*(pap1 + pbp1 - pbp2) - sman*x1) + 2*m2**2*(m3**2 - m4**2 + 2*pap1 + 2*pbp1 + sman*(-1 + x2)) + sman*(2*m4**2 + 4*pap1 - 4*pbp2 + 6*sman - sman*x1)*x2 - sman**2*x2**2) + 
     &                  (m3**2 - m4**2 - 4*pbp2 + sman + sman*x2)*(-4*pap2*pbp1 + 4*pap1*pbp2 + sman*(m3**2 + m4**2 - sman*(-1 + x1 + x2))) + 
     &                  m2**2*(-m3**4 + m4**4 - 8*pap2*pbp1 + 8*pap1*pbp2 + 8*pap2*sman + 4*pbp2*sman - sman**2 + m3**2*(4*pbp2 + sman*(-2 + x1)) - 4*pap2*sman*x1 - 4*pbp2*sman*x1 + sman**2*x1 + 
     &                     sman*(-4*(pap1 + 2*pap2 + pbp2) + sman*x1)*x2 + sman**2*x2**2 + m4**2*(8*pap1 + 8*pap2 + 4*pbp2 - sman*(x1 + 2*x2)))))))))/
     &  (2.*(m1**2 + m4**2 - mq3c**2 - 2*pap2)*(mq1**2 - sman)*(m1**2 - m2**2 + 2*(pbp1 + pbp2) - sman)*(-m3**2 + m4**2 + sman - sman*x2)))


        radSfUf = dreal((g1**2*(A1L*(2*A3Lc*(F1R*m3*(H3Lc*mq3c*(m1**2 - m2**2 + sman) + H3Rc*m1*(m1**2 - m2**2 - 2*(pap2 + pbp2) + sman))*(m3**2 + m4**2 + sman - sman*x2) + 
     &            F1L*mq1*(H3Lc*mq3c*(-2*m3**2*pap1 + 2*m4**2*pap1 - 6*m3**2*pap2 + 2*m4**2*pap2 + 3*m3**2*sman - m4**2*sman + 2*pap1*sman + 2*pap2*sman - sman**2 + m1**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) + 
     &                  sman*(-2*(pap1 + pap2) + sman)*x2 + m2**2*(-3*m3**2 + m4**2 + sman - sman*x2)) + 
     &               H3Rc*m1*(-m3**4 + 6*m3**2*m4**2 - m4**4 - 2*m3**2*pap1 + 2*m4**2*pap1 - 6*m3**2*pap2 + 2*m4**2*pap2 + 3*m3**2*sman - m4**2*sman + 2*pap1*sman + 2*pap2*sman + m3**2*sman*x1 - m4**2*sman*x1 - sman**2*x1 + 
     &                  m1**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) + sman*(-3*m3**2 + m4**2 - 2*(pap1 + pap2) + sman*x1)*x2 + m2**2*(-3*m3**2 + m4**2 + sman - sman*x2)))) + 
     &         A3Rc*(4*F1L*m3*mq1*(H3Rc*m1*mq3c + H3Lc*(m1**2 - pap2))*(m3**2 + m4**2 + sman - sman*x2) - 
     &            F1R*(4*H3Rc*m1*mq3c*(m3**2*(pap1 + 3*pap2 + pbp1 + 3*pbp2 - 3*sman) - (pap1 + pap2 + pbp1 + pbp2 - sman)*(m4**2 + sman - sman*x2)) + 
     &               H3Lc*(-4*m3**2*pap2*pbp1 + 4*m4**2*pap2*pbp1 + 4*m3**2*pap1*pbp2 - 4*m4**2*pap1*pbp2 + m3**4*sman - 6*m3**2*m4**2*sman + m4**4*sman + 6*m3**2*pap2*sman - 2*m4**2*pap2*sman + 4*pap2*pbp1*sman - 
     &                  6*m3**2*pbp2*sman + 2*m4**2*pbp2*sman - 4*pap1*pbp2*sman - 2*pap2*sman**2 + 2*pbp2*sman**2 - sman**3 - m3**2*sman**2*x1 + m4**2*sman**2*x1 + sman**3*x1 - 
     &                  m2**2*(m3**4 + (m4**2 + sman)*(m4**2 + 2*pap2 + 2*pbp2 + sman*(-1 + x1)) - m3**2*(6*(m4**2 + pap2 + pbp2) + sman*x1)) + 
     &                  sman*(4*pap1*pbp2 + 2*pap2*(-2*pbp1 + sman) - sman*(-3*m3**2 + m4**2 + 2*pbp2 + sman*(-1 + x1)))*x2 + m2**2*sman*(-3*m3**2 + m4**2 + 2*pap2 + 2*pbp2 + sman*(-1 + x1))*x2 + 
     &                  m1**2*(m3**4 + m3**2*(-6*m4**2 + 4*pap1 + 6*pap2 + 4*pbp1 + 6*pbp2 - sman*(12 + x1 - 3*x2)) + (m4**2 - 4*pap1 - 2*pap2 - 4*pbp1 - 2*pbp2 + sman*(3 + x1))*(m4**2 + sman - sman*x2)))))) + 
     &      A1R*(2*A3Rc*(F1L*m3*(H3Rc*mq3c*(m1**2 - m2**2 + sman) + H3Lc*m1*(m1**2 - m2**2 - 2*(pap2 + pbp2) + sman))*(m3**2 + m4**2 + sman - sman*x2) + 
     &            F1R*mq1*(H3Rc*mq3c*(-2*m3**2*pap1 + 2*m4**2*pap1 - 6*m3**2*pap2 + 2*m4**2*pap2 + 3*m3**2*sman - m4**2*sman + 2*pap1*sman + 2*pap2*sman - sman**2 + m1**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) + 
     &                  sman*(-2*(pap1 + pap2) + sman)*x2 + m2**2*(-3*m3**2 + m4**2 + sman - sman*x2)) + 
     &               H3Lc*m1*(-m3**4 + 6*m3**2*m4**2 - m4**4 - 2*m3**2*pap1 + 2*m4**2*pap1 - 6*m3**2*pap2 + 2*m4**2*pap2 + 3*m3**2*sman - m4**2*sman + 2*pap1*sman + 2*pap2*sman + m3**2*sman*x1 - m4**2*sman*x1 - sman**2*x1 + 
     &                  m1**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) + sman*(-3*m3**2 + m4**2 - 2*(pap1 + pap2) + sman*x1)*x2 + m2**2*(-3*m3**2 + m4**2 + sman - sman*x2)))) + 
     &         A3Lc*(4*F1R*m3*mq1*(H3Lc*m1*mq3c + H3Rc*(m1**2 - pap2))*(m3**2 + m4**2 + sman - sman*x2) - 
     &            F1L*(4*H3Lc*m1*mq3c*(m3**2*(pap1 + 3*pap2 + pbp1 + 3*pbp2 - 3*sman) - (pap1 + pap2 + pbp1 + pbp2 - sman)*(m4**2 + sman - sman*x2)) + 
     &               H3Rc*(-4*m3**2*pap2*pbp1 + 4*m4**2*pap2*pbp1 + 4*m3**2*pap1*pbp2 - 4*m4**2*pap1*pbp2 + m3**4*sman - 6*m3**2*m4**2*sman + m4**4*sman + 6*m3**2*pap2*sman - 2*m4**2*pap2*sman + 4*pap2*pbp1*sman - 
     &                  6*m3**2*pbp2*sman + 2*m4**2*pbp2*sman - 4*pap1*pbp2*sman - 2*pap2*sman**2 + 2*pbp2*sman**2 - sman**3 - m3**2*sman**2*x1 + m4**2*sman**2*x1 + sman**3*x1 - 
     &                  m2**2*(m3**4 + (m4**2 + sman)*(m4**2 + 2*pap2 + 2*pbp2 + sman*(-1 + x1)) - m3**2*(6*(m4**2 + pap2 + pbp2) + sman*x1)) + 
     &                  sman*(4*pap1*pbp2 + 2*pap2*(-2*pbp1 + sman) - sman*(-3*m3**2 + m4**2 + 2*pbp2 + sman*(-1 + x1)))*x2 + m2**2*sman*(-3*m3**2 + m4**2 + 2*pap2 + 2*pbp2 + sman*(-1 + x1))*x2 + 
     &                  m1**2*(m3**4 + m3**2*(-6*m4**2 + 4*pap1 + 6*pap2 + 4*pbp1 + 6*pbp2 - sman*(12 + x1 - 3*x2)) + (m4**2 - 4*pap1 - 2*pap2 - 4*pbp1 - 2*pbp2 + sman*(3 + x1))*(m4**2 + sman - sman*x2))))))))/
     &  ((m1**2 + m4**2 - mq3c**2 - 2*pap2)*(mq1**2 - sman)*(-m3**2 + m4**2 + sman - sman*x2)**2))


	NeuQ2qx_M2RadiationHSU = NeuQ2qx_M2RadiationHSU + CF*NC*2d0*(radSiUi + radSiUf + radSpUi + radSpUf + radSfUi + radSfUf)
!            print*,"NeuQ2qx_M2RadiationHSU",NeuQ2qx_M2RadiationHSU
 
        enddo

	end

! =================================================================
	double precision function NeuQ2qx_M2RadiationHTU()

        implicit none


#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double precision radTiUi, radTiUf, radTpUi, radTpUf, radTfUi, radTfUf
        double complex A2L, A2R, G2, A3Lc, A3Rc, H3Lc, H3Rc
        double complex g0, g1
        double precision mq2, mq3c
        integer sumT, sumTend, sumUc, sumUcend


        NeuQ2qx_M2RadiationHTU = 0d0

        radTiUi = 0d0
        radTiUf = 0d0
        radTpUi = 0d0
        radTpUf = 0d0
        radTfUi = 0d0
        radTfUf = 0d0

        sumT = 0        
        sumUc = 0

        if ((xtype.eq.5).or.(xtype.eq.6)) then
        
           sumTend = 1
           sumUcend = 3
           
        elseif (xtype.eq.7) then
        
           sumTend = 0
           sumUcend = 3
           
        elseif (xtype.eq.8) then
               
           sumTend = 1
           sumUcend = 1
           
        endif

        do sumT = 0, sumTend
           do sumUc = 0, sumUcend

           A2L = AL(chimeigen,ptmeigen+sumT,pttype,ptgen)
           A2R = AR(chimeigen,ptmeigen+sumT,pttype,ptgen)
           A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
           A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

           G2 = G(xtype,sfmeigen,ptmeigen+sumT,sftype,sfgen)

           H3Lc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
           H3Rc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

           g0 = -gs
	   g1 = -gs
        
           mq2 = MSf(ptmeigen+sumT,pttype,ptgen)
           mq3c = MNeuCha(pumeigen+sumUc)

           radTiUi = dreal((-2*G2*g0**2*(m1**2 + m2**2 + 2*(pbp1 + pbp2) - sman)*(A2R*(2*A3Lc*m3*(H3Lc*m1*mq3c + H3Rc*(m1**2 - pap2)) + A3Rc*(2*H3Rc*mq3c*pap1 + H3Lc*m1*(m3**2 + m4**2 + 2*pap1 - sman*(-1 + x1 + x2)))) + 
     &      A2L*(2*A3Rc*m3*(H3Rc*m1*mq3c + H3Lc*(m1**2 - pap2)) + A3Lc*(2*H3Lc*mq3c*pap1 + H3Rc*m1*(m3**2 + m4**2 + 2*pap1 - sman*(-1 + x1 + x2))))))/
     &  ((m1**2 + m3**2 - mq2**2 - 2*pap1)*(m1**2 + m4**2 - mq3c**2 - 2*pap2)*(-m1**2 + m2**2 - 2*(pbp1 + pbp2) + sman)**2))


           radTiUf = dreal(-((G2*g0*g1*(A2L*(A3Lc*(H3Lc*mq3c*(-2*m3**2*pap1 + 2*m4**2*pap1 + 4*pap2*pbp1 + 4*pap1*pbp2 - m3**2*sman + m4**2*sman - 2*pbp1*sman + sman**2 + m1**2*(m3**2 - m4**2 + 2*pap1 - 2*pbp1 + sman*(-1 + x2)) + 
     &                 m2**2*(m3**2 - m4**2 - 2*pap1 + 2*pbp1 + sman*(-1 + x2)) - sman*(2*pap1 + sman)*x2) + 
     &              H3Rc*m1*(-m3**4 + m4**4 + 2*m4**2*pap1 - 4*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 2*m4**2*sman - 2*pbp1*sman + sman**2 - 2*m2**2*(m4**2 + pap1 - pbp1 + sman) + 
     &                 m3**2*(-2*pap1 + 4*pbp2 + sman*(-2 + x1)) + m2**2*sman*x1 - m4**2*sman*x1 - 2*pbp2*sman*x1 + m1**2*(2*(m3**2 + pap1 - pbp1) - sman*x1) + sman*(2*m2**2 - 2*(m4**2 + pap1 - pbp1 + sman) + sman*x1)*x2 + 
     &                 sman**2*x2**2)) + A3Rc*m3*(2*H3Rc*m1*mq3c*(m1**2 - m2**2 - m3**2 + m4**2 - 2*pbp1 + 2*pbp2 - sman*x2) + 
     &              H3Lc*(2*m1**4 + 2*m3**2*pap2 - 2*m4**2*pap2 + 4*pap2*pbp1 + 4*pap1*pbp2 + m3**2*sman - m4**2*sman - 2*pbp2*sman + sman**2 + m2**2*(-m3**2 + m4**2 + 2*pap2 + 2*pbp2 + sman*(-1 + x1)) - sman**2*x1 + 
     &                 2*pap2*sman*x2 - m1**2*(2*m2**2 + 3*m3**2 - 3*m4**2 + 2*pap2 + 4*pbp1 - 2*pbp2 + sman - sman*x1 + 2*sman*x2)))) + 
     &        A2R*(A3Rc*(H3Rc*mq3c*(-2*m3**2*pap1 + 2*m4**2*pap1 + 4*pap2*pbp1 + 4*pap1*pbp2 - m3**2*sman + m4**2*sman - 2*pbp1*sman + sman**2 + m1**2*(m3**2 - m4**2 + 2*pap1 - 2*pbp1 + sman*(-1 + x2)) + 
     &                 m2**2*(m3**2 - m4**2 - 2*pap1 + 2*pbp1 + sman*(-1 + x2)) - sman*(2*pap1 + sman)*x2) + 
     &              H3Lc*m1*(-m3**4 + m4**4 + 2*m4**2*pap1 - 4*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 2*m4**2*sman - 2*pbp1*sman + sman**2 - 2*m2**2*(m4**2 + pap1 - pbp1 + sman) + 
     &                 m3**2*(-2*pap1 + 4*pbp2 + sman*(-2 + x1)) + m2**2*sman*x1 - m4**2*sman*x1 - 2*pbp2*sman*x1 + m1**2*(2*(m3**2 + pap1 - pbp1) - sman*x1) + sman*(2*m2**2 - 2*(m4**2 + pap1 - pbp1 + sman) + sman*x1)*x2 + 
     &                 sman**2*x2**2)) + A3Lc*m3*(2*H3Lc*m1*mq3c*(m1**2 - m2**2 - m3**2 + m4**2 - 2*pbp1 + 2*pbp2 - sman*x2) + 
     &              H3Rc*(2*m1**4 + 2*m3**2*pap2 - 2*m4**2*pap2 + 4*pap2*pbp1 + 4*pap1*pbp2 + m3**2*sman - m4**2*sman - 2*pbp2*sman + sman**2 + m2**2*(-m3**2 + m4**2 + 2*pap2 + 2*pbp2 + sman*(-1 + x1)) - sman**2*x1 + 
     &                 2*pap2*sman*x2 - m1**2*(2*m2**2 + 3*m3**2 - 3*m4**2 + 2*pap2 + 4*pbp1 - 2*pbp2 + sman - sman*x1 + 2*sman*x2))))))/
     &    ((m1**2 + m3**2 - mq2**2 - 2*pap1)*(m1**2 + m4**2 - mq3c**2 - 2*pap2)*(m1**2 - m2**2 + 2*(pbp1 + pbp2) - sman)*(m3**2 - m4**2 + sman*(-1 + x2)))))


           radTpUi = dreal((G2*g0**2*(4*m1**2 + 4*m2**2 + 2*m3**2 - 2*m4**2 - 2*pap1 - 2*pap2 + 6*pbp1 - 2*pbp2 + sman*(-2 - x1 + x2))*
     &    (A2R*(2*A3Lc*m3*(H3Lc*m1*mq3c + H3Rc*(m1**2 - pap2)) + A3Rc*(2*H3Rc*mq3c*pap1 + H3Lc*m1*(m3**2 + m4**2 + 2*pap1 - sman*(-1 + x1 + x2)))) + 
     &      A2L*(2*A3Rc*m3*(H3Rc*m1*mq3c + H3Lc*(m1**2 - pap2)) + A3Lc*(2*H3Lc*mq3c*pap1 + H3Rc*m1*(m3**2 + m4**2 + 2*pap1 - sman*(-1 + x1 + x2))))))/
     &  (2.*(m1**2 + m3**2 - mq2**2 - 2*pap1)*(m1**2 + m4**2 - mq3c**2 - 2*pap2)*(-m2**2 - m4**2 + mq2**2 + 2*pbp2)*(m1**2 - m2**2 + 2*(pbp1 + pbp2) - sman)))


           radTpUf = dreal((G2*g0*g1*(A2R*(A3Lc*m3*(H3Rc*(2*(2*m1**4 - 4*m4**2*pap1 + 2*m4**2*pap2 - 4*pap1*pap2 + 2*pap2*pbp1 + 2*pap1*pbp2 + m4**2*sman + pap2*sman - pbp2*sman + sman**2 + m3**2*(2*pap2 + sman) + 
     &                  m2**2*(-m3**2 - m4**2 + pap2 + pbp2 + sman*(-1 + x1)) + (pap2 - sman)*sman*x1 - m1**2*(2*m2**2 + 5*m3**2 + m4**2 - 2*pap1 + 3*pap2 + 2*pbp1 - pbp2 + 3*sman - 2*sman*x1)) + 
     &               (5*m1**2 + m2**2 + 2*pap1 - 4*pap2 - sman)*sman*x2) + 2*H3Lc*m1*mq3c*(2*m1**2 - 2*m2**2 - 4*m3**2 - 4*m4**2 + 2*pap1 - 2*pap2 - 2*pbp1 + 2*pbp2 + sman*(-2 + x1 + 3*x2))) + 
     &         A3Rc*(H3Rc*mq3c*(-4*m4**2*pap1 - 8*pap1*pap2 + 4*pap2*pbp1 + 4*pap1*pbp2 + 2*pap1*sman - 2*pbp1*sman - 4*m3**2*(pap1 - 2*pap2 + sman) - 2*m1**2*(2*m4**2 - 3*pap1 + pbp1 + 2*sman) + m1**2*sman*x1 - 
     &               2*pap2*sman*x1 + sman**2*x1 + m2**2*(4*m3**2 - 6*pap1 + 2*pbp1 - sman*x1) + 2*(2*m1**2 + pap1)*sman*x2) + 
     &            H3Lc*m1*(-2*m3**4 - 12*m3**2*m4**2 - 2*m4**4 - 4*m3**2*pap1 + 4*m3**2*pap2 - 8*pap1*pap2 - 4*m4**2*pbp1 + 4*pap2*pbp1 + 4*m3**2*pbp2 + 4*pap1*pbp2 - 6*m3**2*sman - 2*m4**2*sman + 2*pap1*sman - 2*pbp1*sman + 
     &               2*m3**2*sman*x1 + 4*m4**2*sman*x1 - 2*pbp2*sman*x1 + sman**2*x1 + m1**2*(2*m3**2 - 2*m4**2 + 6*pap1 - 2*pbp1 - sman*(2 + x1 - 2*x2)) + 2*sman*(3*m3**2 + 2*m4**2 + pbp1 + sman - sman*x1)*x2 - 
     &               2*sman**2*x2**2 + m2**2*(2*m3**2 - 2*m4**2 - 6*pap1 + 2*pbp1 + sman*(-2 + x1 + 2*x2))))) + 
     &      A2L*(A3Rc*m3*(H3Lc*(2*(2*m1**4 - 4*m4**2*pap1 + 2*m4**2*pap2 - 4*pap1*pap2 + 2*pap2*pbp1 + 2*pap1*pbp2 + m4**2*sman + pap2*sman - pbp2*sman + sman**2 + m3**2*(2*pap2 + sman) + 
     &                  m2**2*(-m3**2 - m4**2 + pap2 + pbp2 + sman*(-1 + x1)) + (pap2 - sman)*sman*x1 - m1**2*(2*m2**2 + 5*m3**2 + m4**2 - 2*pap1 + 3*pap2 + 2*pbp1 - pbp2 + 3*sman - 2*sman*x1)) + 
     &               (5*m1**2 + m2**2 + 2*pap1 - 4*pap2 - sman)*sman*x2) + 2*H3Rc*m1*mq3c*(2*m1**2 - 2*m2**2 - 4*m3**2 - 4*m4**2 + 2*pap1 - 2*pap2 - 2*pbp1 + 2*pbp2 + sman*(-2 + x1 + 3*x2))) + 
     &         A3Lc*(H3Lc*mq3c*(-4*m4**2*pap1 - 8*pap1*pap2 + 4*pap2*pbp1 + 4*pap1*pbp2 + 2*pap1*sman - 2*pbp1*sman - 4*m3**2*(pap1 - 2*pap2 + sman) - 2*m1**2*(2*m4**2 - 3*pap1 + pbp1 + 2*sman) + m1**2*sman*x1 - 
     &               2*pap2*sman*x1 + sman**2*x1 + m2**2*(4*m3**2 - 6*pap1 + 2*pbp1 - sman*x1) + 2*(2*m1**2 + pap1)*sman*x2) + 
     &            H3Rc*m1*(-2*m3**4 - 12*m3**2*m4**2 - 2*m4**4 - 4*m3**2*pap1 + 4*m3**2*pap2 - 8*pap1*pap2 - 4*m4**2*pbp1 + 4*pap2*pbp1 + 4*m3**2*pbp2 + 4*pap1*pbp2 - 6*m3**2*sman - 2*m4**2*sman + 2*pap1*sman - 2*pbp1*sman + 
     &               2*m3**2*sman*x1 + 4*m4**2*sman*x1 - 2*pbp2*sman*x1 + sman**2*x1 + m1**2*(2*m3**2 - 2*m4**2 + 6*pap1 - 2*pbp1 - sman*(2 + x1 - 2*x2)) + 2*sman*(3*m3**2 + 2*m4**2 + pbp1 + sman - sman*x1)*x2 - 
     &               2*sman**2*x2**2 + m2**2*(2*m3**2 - 2*m4**2 - 6*pap1 + 2*pbp1 + sman*(-2 + x1 + 2*x2)))))))/
     &  (2.*(m1**2 + m3**2 - mq2**2 - 2*pap1)*(m1**2 + m4**2 - mq3c**2 - 2*pap2)*(-m2**2 - m4**2 + mq2**2 + 2*pbp2)*(m3**2 - m4**2 + sman*(-1 + x2))))


           radTfUi = dreal((G2*g0*g1*(A2L*(A3Lc*(H3Lc*mq3c*(-2*m3**2*pap1 + 2*m4**2*pap1 + 4*pap2*pbp1 + 4*pap1*pbp2 - m3**2*sman + m4**2*sman - 2*pbp1*sman + sman**2 + m1**2*(m3**2 - m4**2 + 2*pap1 - 2*pbp1 + sman*(-1 + x2)) + 
     &               m2**2*(m3**2 - m4**2 - 2*pap1 + 2*pbp1 + sman*(-1 + x2)) - sman*(2*pap1 + sman)*x2) + 
     &            H3Rc*m1*(-m3**4 + m4**4 + 2*m4**2*pap1 - 4*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 2*m4**2*sman - 2*pbp1*sman + sman**2 - 2*m2**2*(m4**2 + pap1 - pbp1 + sman) + m3**2*(-2*pap1 + 4*pbp2 + sman*(-2 + x1)) + 
     &               m2**2*sman*x1 - m4**2*sman*x1 - 2*pbp2*sman*x1 + m1**2*(2*(m3**2 + pap1 - pbp1) - sman*x1) + sman*(2*m2**2 - 2*(m4**2 + pap1 - pbp1 + sman) + sman*x1)*x2 + sman**2*x2**2)) - 
     &         A3Rc*m3*(2*H3Rc*m1*mq3c*(-m1**2 + m2**2 + m3**2 - m4**2 + 2*pbp1 - 2*pbp2 + sman*x2) + 
     &            H3Lc*(-2*m1**4 - 2*m3**2*pap2 + 2*m4**2*pap2 - 4*pap2*pbp1 + 4*pap1*pbp2 + 8*pap2*pbp2 + m3**2*sman - m4**2*sman - 2*pbp2*sman + sman**2 + m2**2*(-m3**2 + m4**2 - 2*pap2 + 2*pbp2 + sman*(-1 + x1)) - 
     &               sman**2*x1 - 2*pap2*sman*x2 + m1**2*(2*m2**2 + m3**2 - m4**2 + 2*pap2 + 4*pbp1 - 6*pbp2 + sman*(-1 + x1 + 2*x2))))) + 
     &      A2R*(A3Rc*(H3Rc*mq3c*(-2*m3**2*pap1 + 2*m4**2*pap1 + 4*pap2*pbp1 + 4*pap1*pbp2 - m3**2*sman + m4**2*sman - 2*pbp1*sman + sman**2 + m1**2*(m3**2 - m4**2 + 2*pap1 - 2*pbp1 + sman*(-1 + x2)) + 
     &               m2**2*(m3**2 - m4**2 - 2*pap1 + 2*pbp1 + sman*(-1 + x2)) - sman*(2*pap1 + sman)*x2) + 
     &            H3Lc*m1*(-m3**4 + m4**4 + 2*m4**2*pap1 - 4*m4**2*pbp1 + 4*pap2*pbp1 + 4*pap1*pbp2 + 2*m4**2*sman - 2*pbp1*sman + sman**2 - 2*m2**2*(m4**2 + pap1 - pbp1 + sman) + m3**2*(-2*pap1 + 4*pbp2 + sman*(-2 + x1)) + 
     &               m2**2*sman*x1 - m4**2*sman*x1 - 2*pbp2*sman*x1 + m1**2*(2*(m3**2 + pap1 - pbp1) - sman*x1) + sman*(2*m2**2 - 2*(m4**2 + pap1 - pbp1 + sman) + sman*x1)*x2 + sman**2*x2**2)) - 
     &         A3Lc*m3*(2*H3Lc*m1*mq3c*(-m1**2 + m2**2 + m3**2 - m4**2 + 2*pbp1 - 2*pbp2 + sman*x2) + 
     &            H3Rc*(-2*m1**4 - 2*m3**2*pap2 + 2*m4**2*pap2 - 4*pap2*pbp1 + 4*pap1*pbp2 + 8*pap2*pbp2 + m3**2*sman - m4**2*sman - 2*pbp2*sman + sman**2 + m2**2*(-m3**2 + m4**2 - 2*pap2 + 2*pbp2 + sman*(-1 + x1)) - 
     &               sman**2*x1 - 2*pap2*sman*x2 + m1**2*(2*m2**2 + m3**2 - m4**2 + 2*pap2 + 4*pbp1 - 6*pbp2 + sman*(-1 + x1 + 2*x2)))))))/
     &  ((m1**2 + m4**2 - mq3c**2 - 2*pap2)*(m2**2 + m4**2 - mq2**2 - 2*pbp2)*(m1**2 - m2**2 + 2*(pbp1 + pbp2) - sman)*(-m3**2 + m4**2 + sman - sman*x2)))


           radTfUf = dreal((2*G2*g1**2*(A2R*(-2*A3Lc*m3*(H3Lc*m1*mq3c + H3Rc*(m1**2 - pap2))*(m3**2 + m4**2 + sman - sman*x2) + 
     &         A3Rc*(H3Rc*mq3c*(2*m3**2*pap1 - 2*m4**2*pap1 + 6*m3**2*pap2 - 2*m4**2*pap2 - 3*m3**2*sman + m4**2*sman - 2*pap1*sman - 2*pap2*sman + sman**2 + m2**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) + 
     &               (2*(pap1 + pap2) - sman)*sman*x2 + m1**2*(-3*m3**2 + m4**2 + sman - sman*x2)) + 
     &            H3Lc*m1*(m3**4 - 6*m3**2*m4**2 + m4**4 + 2*m3**2*pap1 - 2*m4**2*pap1 + 6*m3**2*pap2 - 2*m4**2*pap2 - 3*m3**2*sman + m4**2*sman - 2*pap1*sman - 2*pap2*sman - m3**2*sman*x1 + m4**2*sman*x1 + sman**2*x1 + 
     &               m2**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) - sman*(-3*m3**2 + m4**2 - 2*(pap1 + pap2) + sman*x1)*x2 + m1**2*(-3*m3**2 + m4**2 + sman - sman*x2)))) + 
     &      A2L*(-2*A3Rc*m3*(H3Rc*m1*mq3c + H3Lc*(m1**2 - pap2))*(m3**2 + m4**2 + sman - sman*x2) + 
     &         A3Lc*(H3Lc*mq3c*(2*m3**2*pap1 - 2*m4**2*pap1 + 6*m3**2*pap2 - 2*m4**2*pap2 - 3*m3**2*sman + m4**2*sman - 2*pap1*sman - 2*pap2*sman + sman**2 + m2**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) + 
     &               (2*(pap1 + pap2) - sman)*sman*x2 + m1**2*(-3*m3**2 + m4**2 + sman - sman*x2)) + 
     &            H3Rc*m1*(m3**4 - 6*m3**2*m4**2 + m4**4 + 2*m3**2*pap1 - 2*m4**2*pap1 + 6*m3**2*pap2 - 2*m4**2*pap2 - 3*m3**2*sman + m4**2*sman - 2*pap1*sman - 2*pap2*sman - m3**2*sman*x1 + m4**2*sman*x1 + sman**2*x1 + 
     &               m2**2*(3*m3**2 - m4**2 + sman*(-1 + x2)) - sman*(-3*m3**2 + m4**2 - 2*(pap1 + pap2) + sman*x1)*x2 + m1**2*(-3*m3**2 + m4**2 + sman - sman*x2))))))/
     &  ((m1**2 + m4**2 - mq3c**2 - 2*pap2)*(m2**2 + m4**2 - mq2**2 - 2*pbp2)*(-m3**2 + m4**2 + sman - sman*x2)**2))

           NeuQ2qx_M2RadiationHTU = NeuQ2qx_M2RadiationHTU + CF*NC*2d0*(radTiUi + radTiUf + radTpUi + radTpUf + radTfUi + radTfUf)
!            print*,"NeuQ2qx_M2RadiationHTU",NeuQ2qx_M2RadiationHTU

           enddo
        enddo
 
	end

