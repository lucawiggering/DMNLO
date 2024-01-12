! =================================================================
!
!
! =================================================================

        double precision function NeuQ2qx_M2RadiationGSSsoft()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double complex radgSiSi, radgSpSp, radgSfSf, radgSpSi, radgSfSi, radgSfSp
        double complex A1L, A1Lc, A1R, A1Rc
        double precision wquark
           
        NeuQ2qx_M2RadiationGSSsoft = 0d0

        radgSiSi = 0d0
        radgSpSp = 0d0
        radgSfSf = 0d0
        
        radgSpSi = 0d0
        radgSfSi = 0d0
        radgSfSp = 0d0

        A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A1R = AR(chimeigen,sfmeigen,sftype,sfgen)
        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        

        wquark = 0d0

! colour class CF**2 CA
        radgSiSi = (CF**2*CV)*((GS**4*(m2**2 - pbp3)*(A1Lc*(2*A1R*m1*m3*(m1**2 + m2**2 + m3**2 + p1p3 - pap1 - 2*pap3 + 2*papb - pbp1 - 2*pbp3) + 
     &         A1L*(m2**2*pap1 - m3**2*pap1 - 4*m3**2*pap3 - 2*p1p3*pap3 + 4*m3**2*papb + 2*p1p3*papb + 
     &            m1**2*(4*m3**2 + 2*p1p3 - pap1 - 2*pbp1) + 2*pap3*pbp1 - 2*papb*pbp1 - 2*pap1*pbp3)) + 
     &      A1Rc*(2*A1L*m1*m3*(m1**2 + m2**2 + m3**2 + p1p3 - pap1 - 2*pap3 + 2*papb - pbp1 - 2*pbp3) + 
     &         A1R*(m2**2*pap1 - m3**2*pap1 - 4*m3**2*pap3 - 2*p1p3*pap3 + 4*m3**2*papb + 2*p1p3*papb + 
     &            m1**2*(4*m3**2 + 2*p1p3 - pap1 - 2*pbp1) + 2*pap3*pbp1 - 2*papb*pbp1 - 2*pap1*pbp3))))/(p1p2**2*pbp3**2))
     
        radgSpSp = (CF**2*CV)*((4*GS**4*(A1Lc*(2*A1R*m1*m3*(m3**4 + m3**2*p1p2 + m1**2*(m3**2 + p1p2) + m2**2*(m3**2 + p1p2) - m3**2*pap1 - 
     &            m3**2*pap2 + p1p2*pap2 + 2*m3**2*papb + 2*p1p2*papb - m3**2*pbp1 - m3**2*pbp2 + p1p2*pbp2) + 
     &         A1L*(-(m3**4*pap1) - m3**4*pap2 + m3**2*p1p2*pap2 + m2**2*(-(p1p2*pap2) + m3**2*(pap1 + pap2)) + 
     &            4*m3**4*papb + 4*m3**2*p1p2*papb - 2*m3**2*papb*pbp1 - 2*m3**2*papb*pbp2 + 2*p1p2*papb*pbp2 + 
     &            m1**2*(4*m3**4 + m3**2*(4*p1p2 - pap1 - pap2 - 2*pbp1 - 2*pbp2) + p1p2*(pap2 + 2*pbp2)))) + 
     &      A1Rc*(2*A1L*m1*m3*(m3**4 + m3**2*p1p2 + m1**2*(m3**2 + p1p2) + m2**2*(m3**2 + p1p2) - m3**2*pap1 - m3**2*pap2 + 
     &            p1p2*pap2 + 2*m3**2*papb + 2*p1p2*papb - m3**2*pbp1 - m3**2*pbp2 + p1p2*pbp2) + 
     &         A1R*(-(m3**4*pap1) - m3**4*pap2 + m3**2*p1p2*pap2 + m2**2*(-(p1p2*pap2) + m3**2*(pap1 + pap2)) + 
     &            4*m3**4*papb + 4*m3**2*p1p2*papb - 2*m3**2*papb*pbp1 - 2*m3**2*papb*pbp2 + 2*p1p2*papb*pbp2 + 
     &            m1**2*(4*m3**4 + m3**2*(4*p1p2 - pap1 - pap2 - 2*pbp1 - 2*pbp2) + p1p2*(pap2 + 2*pbp2))))))/
     &  (p1p2**2*(m1**2 + m2**2 - m3**2 + 2*papb)**2))
     
        radgSfSf = (CF**2*CV)*((4*GS**4*(A1Lc*(2*A1R*m1*m3*(m3**4 + m3**2*p1p3 + m1**2*(m3**2 + p1p3) + m2**2*(m3**2 + p1p3) - m3**2*pap1 - 
     &            m3**2*pap3 + p1p3*pap3 + 2*m3**2*papb + 2*p1p3*papb - m3**2*pbp1 - m3**2*pbp3 + p1p3*pbp3) + 
     &         A1L*(-(m3**4*pap1) - m3**4*pap3 + m3**2*p1p3*pap3 + m2**2*(-(p1p3*pap3) + m3**2*(pap1 + pap3)) + 
     &            4*m3**4*papb + 4*m3**2*p1p3*papb - 2*m3**2*papb*pbp1 - 2*m3**2*papb*pbp3 + 2*p1p3*papb*pbp3 + 
     &            m1**2*(4*m3**4 + m3**2*(4*p1p3 - pap1 - pap3 - 2*pbp1 - 2*pbp3) + p1p3*(pap3 + 2*pbp3)))) + 
     &      A1Rc*(2*A1L*m1*m3*(m3**4 + m3**2*p1p3 + m1**2*(m3**2 + p1p3) + m2**2*(m3**2 + p1p3) - m3**2*pap1 - m3**2*pap3 + 
     &            p1p3*pap3 + 2*m3**2*papb + 2*p1p3*papb - m3**2*pbp1 - m3**2*pbp3 + p1p3*pbp3) + 
     &         A1R*(-(m3**4*pap1) - m3**4*pap3 + m3**2*p1p3*pap3 + m2**2*(-(p1p3*pap3) + m3**2*(pap1 + pap3)) + 
     &            4*m3**4*papb + 4*m3**2*p1p3*papb - 2*m3**2*papb*pbp1 - 2*m3**2*papb*pbp3 + 2*p1p3*papb*pbp3 + 
     &            m1**2*(4*m3**4 + m3**2*(4*p1p3 - pap1 - pap3 - 2*pbp1 - 2*pbp3) + p1p3*(pap3 + 2*pbp3))))))/
     &  (p1p3**2*(m1**2 + m2**2 - m3**2 + 2*papb)**2))
     
        radgSpSi = (CF**2*CV)*(-(GS**4*(A1L*(A1Rc*m1*m3*(-(m3**2*p1p3) - 2*m3**2*p2p3 + 2*p2p3*pap1 - 2*p1p3*pap2 - 2*m3**2*pap3 + 2*p1p2*pap3 + 
     &             2*p1p3*pap3 + 4*p2p3*pap3 + 4*m3**2*papb - 4*p1p2*papb - 4*p1p3*papb - 8*p2p3*papb + 2*m3**2*pbp1 + 
     &             4*pap2*pbp1 - 2*pap3*pbp1 + 4*papb*pbp1 - m1**2*(p1p3 + 2*p2p3 - 2*pbp1 - 4*pbp2) + 4*m3**2*pbp2 - 
     &             4*pap1*pbp2 - 4*pap3*pbp2 + 8*papb*pbp2 + m2**2*(4*m3**2 - 4*p1p2 - 3*p1p3 - 6*p2p3 + 2*pbp1 + 4*pbp2) - 
     &             4*m3**2*pbp3 + 4*p1p2*pbp3 + 2*p1p3*pbp3 + 4*p2p3*pbp3 + 2*pap1*pbp3 + 4*pap2*pbp3) + 
     &          A1Lc*(m3**2*p2p3*pap1 - m3**2*p1p3*pap2 - m3**4*pap3 + m3**2*p1p2*pap3 + 2*m3**2*p1p3*pap3 + 
     &             4*m3**2*p2p3*pap3 + 2*m3**4*papb - 2*m3**2*p1p2*papb - 4*m3**2*p1p3*papb - 8*m3**2*p2p3*papb + 
     &             2*m3**2*pap2*pbp1 - 2*m3**2*pap3*pbp1 - 2*p2p3*pap3*pbp1 + 4*m3**2*papb*pbp1 + 2*p2p3*papb*pbp1 - 
     &             2*m3**2*pap1*pbp2 - 4*m3**2*pap3*pbp2 + 2*p1p3*pap3*pbp2 + 8*m3**2*papb*pbp2 - 2*p1p3*papb*pbp2 + 
     &             m2**2*(-3*p2p3*pap1 + 3*p1p3*pap2 - m3**2*pap3 + p1p2*pap3 + 2*m3**2*papb - 2*p1p2*papb - 2*pap2*pbp1 + 
     &                2*pap1*pbp2) + 2*m3**2*pap1*pbp3 + 2*p2p3*pap1*pbp3 + 4*m3**2*pap2*pbp3 - 2*p1p3*pap2*pbp3 - 
     &             2*m3**2*papb*pbp3 + 2*p1p2*papb*pbp3 + 
     &             m1**2*(4*m2**2*(m3**2 - p1p2) + p2p3*pap1 - p1p3*pap2 + p1p2*pap3 - 2*p1p2*papb + 2*pap2*pbp1 - 
     &                2*pap1*pbp2 + 4*p1p2*pbp3 - m3**2*(2*p1p3 + 4*p2p3 + pap3 - 2*papb - 4*pbp1 - 8*pbp2 + 4*pbp3)))) + 
     &       A1R*(A1Lc*m1*m3*(-(m3**2*p1p3) - 2*m3**2*p2p3 + 2*p2p3*pap1 - 2*p1p3*pap2 - 2*m3**2*pap3 + 2*p1p2*pap3 + 
     &             2*p1p3*pap3 + 4*p2p3*pap3 + 4*m3**2*papb - 4*p1p2*papb - 4*p1p3*papb - 8*p2p3*papb + 2*m3**2*pbp1 + 
     &             4*pap2*pbp1 - 2*pap3*pbp1 + 4*papb*pbp1 - m1**2*(p1p3 + 2*p2p3 - 2*pbp1 - 4*pbp2) + 4*m3**2*pbp2 - 
     &             4*pap1*pbp2 - 4*pap3*pbp2 + 8*papb*pbp2 + m2**2*(4*m3**2 - 4*p1p2 - 3*p1p3 - 6*p2p3 + 2*pbp1 + 4*pbp2) - 
     &             4*m3**2*pbp3 + 4*p1p2*pbp3 + 2*p1p3*pbp3 + 4*p2p3*pbp3 + 2*pap1*pbp3 + 4*pap2*pbp3) + 
     &          A1Rc*(m3**2*p2p3*pap1 - m3**2*p1p3*pap2 - m3**4*pap3 + m3**2*p1p2*pap3 + 2*m3**2*p1p3*pap3 + 
     &             4*m3**2*p2p3*pap3 + 2*m3**4*papb - 2*m3**2*p1p2*papb - 4*m3**2*p1p3*papb - 8*m3**2*p2p3*papb + 
     &             2*m3**2*pap2*pbp1 - 2*m3**2*pap3*pbp1 - 2*p2p3*pap3*pbp1 + 4*m3**2*papb*pbp1 + 2*p2p3*papb*pbp1 - 
     &             2*m3**2*pap1*pbp2 - 4*m3**2*pap3*pbp2 + 2*p1p3*pap3*pbp2 + 8*m3**2*papb*pbp2 - 2*p1p3*papb*pbp2 + 
     &             m2**2*(-3*p2p3*pap1 + 3*p1p3*pap2 - m3**2*pap3 + p1p2*pap3 + 2*m3**2*papb - 2*p1p2*papb - 2*pap2*pbp1 + 
     &                2*pap1*pbp2) + 2*m3**2*pap1*pbp3 + 2*p2p3*pap1*pbp3 + 4*m3**2*pap2*pbp3 - 2*p1p3*pap2*pbp3 - 
     &             2*m3**2*papb*pbp3 + 2*p1p2*papb*pbp3 + 
     &             m1**2*(4*m2**2*(m3**2 - p1p2) + p2p3*pap1 - p1p3*pap2 + p1p2*pap3 - 2*p1p2*papb + 2*pap2*pbp1 - 
     &                2*pap1*pbp2 + 4*p1p2*pbp3 - m3**2*(2*p1p3 + 4*p2p3 + pap3 - 2*papb - 4*pbp1 - 8*pbp2 + 4*pbp3))))))/
     &  (2.*p1p2**2*(m1**2 + m2**2 - m3**2 + 2*papb)*pbp3))
     
! colour class -CF/2
        radgSfSi = (-CF/2d0)*((GS**4*(A1L*(A1Rc*m1*m3*(p1p3**2 - 2*p1p3*pap1 - 2*p1p3*pap3 + 2*p1p3*papb - 4*m2**2*pbp1 - 4*p1p3*pbp1 + 
     &            4*pap1*pbp1 + 6*pap3*pbp1 - 8*papb*pbp1 + 4*pbp1**2 + 2*m1**2*(p1p3 - 2*pbp1 - pbp3) + 
     &            2*m3**2*(p1p3 - 2*pbp1 - pbp3) - 2*m2**2*pbp3 - 2*p1p3*pbp3 + 2*pap1*pbp3 + 2*pap3*pbp3 - 4*papb*pbp3 + 
     &            8*pbp1*pbp3 + 2*pbp3**2) + A1Lc*(-(m3**2*p1p3*pap1) - 2*m3**2*p1p3*pap3 + 3*m3**2*p1p3*papb + 
     &            p1p3**2*papb + 2*m3**2*pap1*pbp1 + 5*m3**2*pap3*pbp1 + p1p3*pap3*pbp1 - 8*m3**2*papb*pbp1 - 
     &            4*p1p3*papb*pbp1 - 2*pap3*pbp1**2 + 4*papb*pbp1**2 + m3**2*pap1*pbp3 - p1p3*pap1*pbp3 + 2*m3**2*pap3*pbp3 - 
     &            4*m3**2*papb*pbp3 + 2*pap1*pbp1*pbp3 + 4*papb*pbp1*pbp3 + 
     &            m1**2*(-2*m2**2*p1p3 + p1p3**2 - p1p3*pap1 - p1p3*papb - 4*p1p3*pbp1 + 2*pap1*pbp1 + pap3*pbp1 + 
     &               4*pbp1**2 + 4*m3**2*(p1p3 - 2*pbp1 - pbp3) + pap1*pbp3 + 4*pbp1*pbp3) + 
     &            m2**2*(p1p3*(pap1 - papb) - pap3*pbp1 - pap1*(2*pbp1 + pbp3)))) + 
     &      A1R*(A1Lc*m1*m3*(p1p3**2 - 2*p1p3*pap1 - 2*p1p3*pap3 + 2*p1p3*papb - 4*m2**2*pbp1 - 4*p1p3*pbp1 + 4*pap1*pbp1 + 
     &            6*pap3*pbp1 - 8*papb*pbp1 + 4*pbp1**2 + 2*m1**2*(p1p3 - 2*pbp1 - pbp3) + 2*m3**2*(p1p3 - 2*pbp1 - pbp3) - 
     &            2*m2**2*pbp3 - 2*p1p3*pbp3 + 2*pap1*pbp3 + 2*pap3*pbp3 - 4*papb*pbp3 + 8*pbp1*pbp3 + 2*pbp3**2) + 
     &         A1Rc*(-(m3**2*p1p3*pap1) - 2*m3**2*p1p3*pap3 + 3*m3**2*p1p3*papb + p1p3**2*papb + 2*m3**2*pap1*pbp1 + 
     &            5*m3**2*pap3*pbp1 + p1p3*pap3*pbp1 - 8*m3**2*papb*pbp1 - 4*p1p3*papb*pbp1 - 2*pap3*pbp1**2 + 
     &            4*papb*pbp1**2 + m3**2*pap1*pbp3 - p1p3*pap1*pbp3 + 2*m3**2*pap3*pbp3 - 4*m3**2*papb*pbp3 + 
     &            2*pap1*pbp1*pbp3 + 4*papb*pbp1*pbp3 + 
     &            m1**2*(-2*m2**2*p1p3 + p1p3**2 - p1p3*pap1 - p1p3*papb - 4*p1p3*pbp1 + 2*pap1*pbp1 + pap3*pbp1 + 
     &               4*pbp1**2 + 4*m3**2*(p1p3 - 2*pbp1 - pbp3) + pap1*pbp3 + 4*pbp1*pbp3) + 
     &            m2**2*(p1p3*(pap1 - papb) - pap3*pbp1 - pap1*(2*pbp1 + pbp3))))))/
     &  (p1p2*p1p3*(m1**2 + m2**2 - m3**2 + 2*papb)*pbp3))
     
        radgSfSp = (-CF/2d0)*((2*GS**4*(A1Lc*(A1R*m1*m3*(m3**4 + 2*m3**2*p1p2 + 2*m3**2*p1p3 + m3**2*p2p3 + 
     &            m1**2*(m3**2 + 2*p1p2 + 2*p1p3 + p2p3) + m2**2*(m3**2 + 2*p1p2 + 2*p1p3 + p2p3) + 2*m3**2*pap1 - 
     &            4*p1p2*pap1 - 4*p1p3*pap1 - 4*p2p3*pap1 + 2*m3**2*pap2 + 2*m3**2*pap3 + 2*m3**2*papb + 4*p1p2*papb + 
     &            4*p1p3*papb + 2*p2p3*papb + 2*m3**2*pbp1 - 4*p1p2*pbp1 - 4*p1p3*pbp1 - 4*p2p3*pbp1 + 2*m3**2*pbp2 + 
     &            2*m3**2*pbp3) + A1L*(m3**4*pap1 - 2*m3**2*p1p2*pap1 - 2*m3**2*p1p3*pap1 - 2*m3**2*p2p3*pap1 + m3**4*pap2 + 
     &            m3**4*pap3 - m2**2*(-2*(p1p2 + p1p3 + p2p3)*pap1 + m3**2*(pap1 + pap2 + pap3)) + 2*m3**4*papb + 
     &            4*m3**2*p1p2*papb + 4*m3**2*p1p3*papb + 2*m3**2*p2p3*papb + 2*m3**2*papb*pbp1 - 4*p1p2*papb*pbp1 - 
     &            4*p1p3*papb*pbp1 - 4*p2p3*papb*pbp1 + 2*m3**2*papb*pbp2 + 2*m3**2*papb*pbp3 + 
     &            m1**2*(2*m3**4 - 2*(p1p2 + p1p3 + p2p3)*(pap1 + 2*pbp1) + 
     &               m3**2*(4*p1p2 + 4*p1p3 + 2*p2p3 + pap1 + pap2 + pap3 + 2*pbp1 + 2*pbp2 + 2*pbp3)))) + 
     &      A1Rc*(A1L*m1*m3*(m3**4 + 2*m3**2*p1p2 + 2*m3**2*p1p3 + m3**2*p2p3 + m1**2*(m3**2 + 2*p1p2 + 2*p1p3 + p2p3) + 
     &            m2**2*(m3**2 + 2*p1p2 + 2*p1p3 + p2p3) + 2*m3**2*pap1 - 4*p1p2*pap1 - 4*p1p3*pap1 - 4*p2p3*pap1 + 
     &            2*m3**2*pap2 + 2*m3**2*pap3 + 2*m3**2*papb + 4*p1p2*papb + 4*p1p3*papb + 2*p2p3*papb + 2*m3**2*pbp1 - 
     &            4*p1p2*pbp1 - 4*p1p3*pbp1 - 4*p2p3*pbp1 + 2*m3**2*pbp2 + 2*m3**2*pbp3) + 
     &         A1R*(m3**4*pap1 - 2*m3**2*p1p2*pap1 - 2*m3**2*p1p3*pap1 - 2*m3**2*p2p3*pap1 + m3**4*pap2 + m3**4*pap3 - 
     &            m2**2*(-2*(p1p2 + p1p3 + p2p3)*pap1 + m3**2*(pap1 + pap2 + pap3)) + 2*m3**4*papb + 4*m3**2*p1p2*papb + 
     &            4*m3**2*p1p3*papb + 2*m3**2*p2p3*papb + 2*m3**2*papb*pbp1 - 4*p1p2*papb*pbp1 - 4*p1p3*papb*pbp1 - 
     &            4*p2p3*papb*pbp1 + 2*m3**2*papb*pbp2 + 2*m3**2*papb*pbp3 + 
     &            m1**2*(2*m3**4 - 2*(p1p2 + p1p3 + p2p3)*(pap1 + 2*pbp1) + 
     &               m3**2*(4*p1p2 + 4*p1p3 + 2*p2p3 + pap1 + pap2 + pap3 + 2*pbp1 + 2*pbp2 + 2*pbp3))))))/
     &  (p1p2*p1p3*(m1**2 + m2**2 - m3**2 + 2*papb)**2))
     

        NeuQ2qx_M2RadiationGSSsoft = 0.5*(radgSiSi + radgSpSp + radgSfSf + 2d0*(radgSpSi + radgSfSi + radgSfSp))
        
!        print*,"NeuQ2qx_M2RadiationGSSsoft",NeuQ2qx_M2RadiationGSSsoft

        
       end

        double precision function NeuQ2qx_M2RadiationGSSsoftcoll()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double complex radgSfSg, radgSiSg, radgSpSg
        double complex A1L, A1Lc, A1R, A1Rc
        double precision wquark
           
        NeuQ2qx_M2RadiationGSSsoftcoll = 0d0

        radgSfSg = 0d0
        radgSiSg = 0d0
        radgSpSg = 0d0


        A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A1R = AR(chimeigen,sfmeigen,sftype,sfgen)
        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        wquark = 0d0
     
! colour class -1/2 CA**2 CF
        radgSfSg = (-CV**2*CF/2d0)*((GS**4*(A1Lc*(A1R*m1*m3*(6*m3**2*p1p2 + 3*m3**2*p1p3 + 3*m3**2*p2p3 + 3*m1**2*(2*p1p2 + p1p3 + p2p3) + 
     &            3*m2**2*(2*p1p2 + p1p3 + p2p3) - 8*p1p2*pap1 - 4*p1p3*pap1 - 4*p2p3*pap1 - 4*m3**2*pap2 + 8*p1p3*pap2 - 
     &            2*m3**2*pap3 - 4*p1p2*pap3 - 4*p1p3*pap3 + 12*p1p2*papb + 6*p1p3*papb + 6*p2p3*papb - 8*p1p2*pbp1 - 
     &            4*p1p3*pbp1 - 4*p2p3*pbp1 - 4*m3**2*pbp2 + 8*p1p3*pbp2 - 2*m3**2*pbp3 - 4*p1p2*pbp3 - 4*p1p3*pbp3) + 
     &         A1L*(-4*m3**2*p1p2*pap1 - 2*m3**2*p1p3*pap1 - 2*m3**2*p2p3*pap1 - 2*m3**4*pap2 + 4*m3**2*p1p3*pap2 - 
     &            m3**4*pap3 - 2*m3**2*p1p2*pap3 - 2*m3**2*p1p3*pap3 + 
     &            m2**2*(2*p2p3*pap1 + 2*m3**2*pap2 + m3**2*pap3 + 2*p1p2*(2*pap1 + pap3) + 2*p1p3*(pap1 - 2*pap2 + pap3)) + 
     &            12*m3**2*p1p2*papb + 6*m3**2*p1p3*papb + 6*m3**2*p2p3*papb - 8*p1p2*papb*pbp1 - 4*p1p3*papb*pbp1 - 
     &            4*p2p3*papb*pbp1 - 4*m3**2*papb*pbp2 + 8*p1p3*papb*pbp2 - 2*m3**2*papb*pbp3 - 4*p1p2*papb*pbp3 - 
     &            4*p1p3*papb*pbp3 + m1**2*(m3**2*(12*p1p2 + 6*p1p3 + 6*p2p3 - 2*pap2 - pap3 - 4*pbp2 - 2*pbp3) - 
     &               2*(p2p3*(pap1 + 2*pbp1) + p1p2*(2*pap1 + pap3 + 4*pbp1 + 2*pbp3) + 
     &                  p1p3*(pap1 - 2*pap2 + pap3 + 2*pbp1 - 4*pbp2 + 2*pbp3))))) + 
     &      A1Rc*(A1L*m1*m3*(6*m3**2*p1p2 + 3*m3**2*p1p3 + 3*m3**2*p2p3 + 3*m1**2*(2*p1p2 + p1p3 + p2p3) + 
     &            3*m2**2*(2*p1p2 + p1p3 + p2p3) - 8*p1p2*pap1 - 4*p1p3*pap1 - 4*p2p3*pap1 - 4*m3**2*pap2 + 8*p1p3*pap2 - 
     &            2*m3**2*pap3 - 4*p1p2*pap3 - 4*p1p3*pap3 + 12*p1p2*papb + 6*p1p3*papb + 6*p2p3*papb - 8*p1p2*pbp1 - 
     &            4*p1p3*pbp1 - 4*p2p3*pbp1 - 4*m3**2*pbp2 + 8*p1p3*pbp2 - 2*m3**2*pbp3 - 4*p1p2*pbp3 - 4*p1p3*pbp3) + 
     &         A1R*(-4*m3**2*p1p2*pap1 - 2*m3**2*p1p3*pap1 - 2*m3**2*p2p3*pap1 - 2*m3**4*pap2 + 4*m3**2*p1p3*pap2 - 
     &            m3**4*pap3 - 2*m3**2*p1p2*pap3 - 2*m3**2*p1p3*pap3 + 
     &            m2**2*(2*p2p3*pap1 + 2*m3**2*pap2 + m3**2*pap3 + 2*p1p2*(2*pap1 + pap3) + 2*p1p3*(pap1 - 2*pap2 + pap3)) + 
     &            12*m3**2*p1p2*papb + 6*m3**2*p1p3*papb + 6*m3**2*p2p3*papb - 8*p1p2*papb*pbp1 - 4*p1p3*papb*pbp1 - 
     &            4*p2p3*papb*pbp1 - 4*m3**2*papb*pbp2 + 8*p1p3*papb*pbp2 - 2*m3**2*papb*pbp3 - 4*p1p2*papb*pbp3 - 
     &            4*p1p3*papb*pbp3 + m1**2*(m3**2*(12*p1p2 + 6*p1p3 + 6*p2p3 - 2*pap2 - pap3 - 4*pbp2 - 2*pbp3) - 
     &               2*(p2p3*(pap1 + 2*pbp1) + p1p2*(2*pap1 + pap3 + 4*pbp1 + 2*pbp3) + 
     &                  p1p3*(pap1 - 2*pap2 + pap3 + 2*pbp1 - 4*pbp2 + 2*pbp3)))))))/
     &  (p1p3*p2p3*(m1**2 + m2**2 - m3**2 + 2*papb)**2))
     
! colour class +1/2 CA**2 CF
        radgSiSg = (+CV**2*CF/2d0)*((GS**4*(A1L*(A1Rc*m1*m3*(3*m3**2*p2p3 + 2*p1p3*p2p3 - 2*p2p3*pap1 - 2*p1p3*pap2 - 2*p1p2*pap3 - 2*p1p3*pap3 - 
     &            3*p2p3*pap3 + 4*p1p2*papb + 2*p1p3*papb + 6*p2p3*papb - 4*p2p3*pbp1 + 4*pap2*pbp1 + 2*pap3*pbp1 - 
     &            6*m3**2*pbp2 - 4*p1p3*pbp2 + 4*pap1*pbp2 + 6*pap3*pbp2 - 12*papb*pbp2 + 8*pbp1*pbp2 + 
     &            m2**2*(4*p1p2 + 2*p1p3 + 3*p2p3 - 6*pbp2 - 3*pbp3) + 3*m1**2*(p2p3 - 2*pbp2 - pbp3) - 3*m3**2*pbp3 - 
     &            4*p1p2*pbp3 - 4*p1p3*pbp3 - 3*p2p3*pbp3 + 2*pap1*pbp3 + 6*pap3*pbp3 - 6*papb*pbp3 + 4*pbp1*pbp3 + 
     &            6*pbp2*pbp3 + 6*pbp3**2) + A1Lc*(-(m3**2*p2p3*pap1) - m3**2*p1p3*pap2 - m3**2*p1p2*pap3 - m3**2*p1p3*pap3 - 
     &            3*m3**2*p2p3*pap3 + 2*m3**2*p1p2*papb + m3**2*p1p3*papb + 6*m3**2*p2p3*papb + 2*p1p3*p2p3*papb + 
     &            2*m3**2*pap2*pbp1 + m3**2*pap3*pbp1 + p2p3*pap3*pbp1 - 4*p2p3*papb*pbp1 + 2*m3**2*pap1*pbp2 + 
     &            6*m3**2*pap3*pbp2 + p1p3*pap3*pbp2 - 12*m3**2*papb*pbp2 - 4*p1p3*papb*pbp2 - 4*pap3*pbp1*pbp2 + 
     &            8*papb*pbp1*pbp2 + m3**2*pap1*pbp3 - p2p3*pap1*pbp3 - p1p3*pap2*pbp3 + 6*m3**2*pap3*pbp3 - 
     &            6*m3**2*papb*pbp3 - 2*p1p2*papb*pbp3 - 6*p1p3*papb*pbp3 + 2*pap2*pbp1*pbp3 - 4*pap3*pbp1*pbp3 + 
     &            4*papb*pbp1*pbp3 + 2*pap1*pbp2*pbp3 + 4*pap1*pbp3**2 + 
     &            m2**2*(p2p3*pap1 - p1p2*pap3 + 2*p1p2*papb + p1p3*(pap2 + 3*pap3 + papb) - 2*pap2*pbp1 - pap3*pbp1 - 
     &               2*pap1*pbp2 - pap1*pbp3) + m1**2*
     &             (2*m2**2*(2*p1p2 + p1p3) + 2*p1p3*p2p3 - p2p3*pap1 - p1p3*pap2 - p1p2*pap3 - p1p3*pap3 + 2*p1p2*papb + 
     &               p1p3*papb - 4*p2p3*pbp1 + 2*pap2*pbp1 + pap3*pbp1 - 4*p1p3*pbp2 + 2*pap1*pbp2 + 8*pbp1*pbp2 + 
     &               6*m3**2*(p2p3 - 2*pbp2 - pbp3) - 4*p1p2*pbp3 - 4*p1p3*pbp3 + pap1*pbp3 + 4*pbp1*pbp3))) + 
     &      A1R*(A1Lc*m1*m3*(3*m3**2*p2p3 + 2*p1p3*p2p3 - 2*p2p3*pap1 - 2*p1p3*pap2 - 2*p1p2*pap3 - 2*p1p3*pap3 - 
     &            3*p2p3*pap3 + 4*p1p2*papb + 2*p1p3*papb + 6*p2p3*papb - 4*p2p3*pbp1 + 4*pap2*pbp1 + 2*pap3*pbp1 - 
     &            6*m3**2*pbp2 - 4*p1p3*pbp2 + 4*pap1*pbp2 + 6*pap3*pbp2 - 12*papb*pbp2 + 8*pbp1*pbp2 + 
     &            m2**2*(4*p1p2 + 2*p1p3 + 3*p2p3 - 6*pbp2 - 3*pbp3) + 3*m1**2*(p2p3 - 2*pbp2 - pbp3) - 3*m3**2*pbp3 - 
     &            4*p1p2*pbp3 - 4*p1p3*pbp3 - 3*p2p3*pbp3 + 2*pap1*pbp3 + 6*pap3*pbp3 - 6*papb*pbp3 + 4*pbp1*pbp3 + 
     &            6*pbp2*pbp3 + 6*pbp3**2) + A1Rc*(-(m3**2*p2p3*pap1) - m3**2*p1p3*pap2 - m3**2*p1p2*pap3 - m3**2*p1p3*pap3 - 
     &            3*m3**2*p2p3*pap3 + 2*m3**2*p1p2*papb + m3**2*p1p3*papb + 6*m3**2*p2p3*papb + 2*p1p3*p2p3*papb + 
     &            2*m3**2*pap2*pbp1 + m3**2*pap3*pbp1 + p2p3*pap3*pbp1 - 4*p2p3*papb*pbp1 + 2*m3**2*pap1*pbp2 + 
     &            6*m3**2*pap3*pbp2 + p1p3*pap3*pbp2 - 12*m3**2*papb*pbp2 - 4*p1p3*papb*pbp2 - 4*pap3*pbp1*pbp2 + 
     &            8*papb*pbp1*pbp2 + m3**2*pap1*pbp3 - p2p3*pap1*pbp3 - p1p3*pap2*pbp3 + 6*m3**2*pap3*pbp3 - 
     &            6*m3**2*papb*pbp3 - 2*p1p2*papb*pbp3 - 6*p1p3*papb*pbp3 + 2*pap2*pbp1*pbp3 - 4*pap3*pbp1*pbp3 + 
     &            4*papb*pbp1*pbp3 + 2*pap1*pbp2*pbp3 + 4*pap1*pbp3**2 + 
     &            m2**2*(p2p3*pap1 - p1p2*pap3 + 2*p1p2*papb + p1p3*(pap2 + 3*pap3 + papb) - 2*pap2*pbp1 - pap3*pbp1 - 
     &               2*pap1*pbp2 - pap1*pbp3) + m1**2*
     &             (2*m2**2*(2*p1p2 + p1p3) + 2*p1p3*p2p3 - p2p3*pap1 - p1p3*pap2 - p1p2*pap3 - p1p3*pap3 + 2*p1p2*papb + 
     &               p1p3*papb - 4*p2p3*pbp1 + 2*pap2*pbp1 + pap3*pbp1 - 4*p1p3*pbp2 + 2*pap1*pbp2 + 8*pbp1*pbp2 + 
     &               6*m3**2*(p2p3 - 2*pbp2 - pbp3) - 4*p1p2*pbp3 - 4*p1p3*pbp3 + pap1*pbp3 + 4*pbp1*pbp3)))))/
     &  (2.*p1p2*p2p3*(m1**2 + m2**2 - m3**2 + 2*papb)*pbp3))
     
        radgSpSg = (+CV**2*CF/2d0)*(-((GS**4*(A1Lc*(A1R*m1*m3*(3*m3**2*p1p2 + 6*m3**2*p1p3 + 3*m3**2*p2p3 + 3*m1**2*(p1p2 + 2*p1p3 + p2p3) + 
     &              3*m2**2*(p1p2 + 2*p1p3 + p2p3) - 4*p1p2*pap1 - 8*p1p3*pap1 - 4*p2p3*pap1 - 2*m3**2*pap2 - 4*p1p2*pap2 - 
     &              4*p1p3*pap2 - 4*m3**2*pap3 + 8*p1p2*pap3 + 6*p1p2*papb + 12*p1p3*papb + 6*p2p3*papb - 4*p1p2*pbp1 - 
     &              8*p1p3*pbp1 - 4*p2p3*pbp1 - 2*m3**2*pbp2 - 4*p1p2*pbp2 - 4*p1p3*pbp2 - 4*m3**2*pbp3 + 8*p1p2*pbp3) + 
     &           A1L*(-2*m3**2*p1p2*pap1 - 4*m3**2*p1p3*pap1 - 2*m3**2*p2p3*pap1 - m3**4*pap2 - 2*m3**2*p1p2*pap2 - 
     &              2*m3**2*p1p3*pap2 - 2*m3**4*pap3 + 4*m3**2*p1p2*pap3 + 
     &              m2**2*(2*p2p3*pap1 + m3**2*pap2 + 2*p1p3*(2*pap1 + pap2) + 2*p1p2*(pap1 + pap2 - 2*pap3) + 
     &                 2*m3**2*pap3) + 6*m3**2*p1p2*papb + 12*m3**2*p1p3*papb + 6*m3**2*p2p3*papb - 4*p1p2*papb*pbp1 - 
     &              8*p1p3*papb*pbp1 - 4*p2p3*papb*pbp1 - 2*m3**2*papb*pbp2 - 4*p1p2*papb*pbp2 - 4*p1p3*papb*pbp2 + 
     &              m1**2*(-2*(p2p3*(pap1 + 2*pbp1) + p1p3*(2*pap1 + pap2 + 4*pbp1 + 2*pbp2) + 
     &                    p1p2*(pap1 + pap2 + 2*(-pap3 + pbp1 + pbp2 - 2*pbp3))) + 
     &                 m3**2*(6*p1p2 + 12*p1p3 + 6*p2p3 - pap2 - 2*pap3 - 2*pbp2 - 4*pbp3)) - 4*m3**2*papb*pbp3 + 
     &              8*p1p2*papb*pbp3)) + A1Rc*(A1L*m1*m3*
     &            (3*m3**2*p1p2 + 6*m3**2*p1p3 + 3*m3**2*p2p3 + 3*m1**2*(p1p2 + 2*p1p3 + p2p3) + 
     &              3*m2**2*(p1p2 + 2*p1p3 + p2p3) - 4*p1p2*pap1 - 8*p1p3*pap1 - 4*p2p3*pap1 - 2*m3**2*pap2 - 4*p1p2*pap2 - 
     &              4*p1p3*pap2 - 4*m3**2*pap3 + 8*p1p2*pap3 + 6*p1p2*papb + 12*p1p3*papb + 6*p2p3*papb - 4*p1p2*pbp1 - 
     &              8*p1p3*pbp1 - 4*p2p3*pbp1 - 2*m3**2*pbp2 - 4*p1p2*pbp2 - 4*p1p3*pbp2 - 4*m3**2*pbp3 + 8*p1p2*pbp3) + 
     &           A1R*(-2*m3**2*p1p2*pap1 - 4*m3**2*p1p3*pap1 - 2*m3**2*p2p3*pap1 - m3**4*pap2 - 2*m3**2*p1p2*pap2 - 
     &              2*m3**2*p1p3*pap2 - 2*m3**4*pap3 + 4*m3**2*p1p2*pap3 + 
     &              m2**2*(2*p2p3*pap1 + m3**2*pap2 + 2*p1p3*(2*pap1 + pap2) + 2*p1p2*(pap1 + pap2 - 2*pap3) + 
     &                 2*m3**2*pap3) + 6*m3**2*p1p2*papb + 12*m3**2*p1p3*papb + 6*m3**2*p2p3*papb - 4*p1p2*papb*pbp1 - 
     &              8*p1p3*papb*pbp1 - 4*p2p3*papb*pbp1 - 2*m3**2*papb*pbp2 - 4*p1p2*papb*pbp2 - 4*p1p3*papb*pbp2 + 
     &              m1**2*(-2*(p2p3*(pap1 + 2*pbp1) + p1p3*(2*pap1 + pap2 + 4*pbp1 + 2*pbp2) + 
     &                    p1p2*(pap1 + pap2 + 2*(-pap3 + pbp1 + pbp2 - 2*pbp3))) + 
     &                 m3**2*(6*p1p2 + 12*p1p3 + 6*p2p3 - pap2 - 2*pap3 - 2*pbp2 - 4*pbp3)) - 4*m3**2*papb*pbp3 + 
     &              8*p1p2*papb*pbp3))))/(p1p2*p2p3*(m1**2 + m2**2 - m3**2 + 2*papb)**2)))



        NeuQ2qx_M2RadiationGSSsoftcoll = 0.5*(2d0*(radgSfSg + radgSiSg + radgSpSg))

  
!        print*,"NeuQ2qx_M2RadiationGSSsoftcoll",NeuQ2qx_M2RadiationGSSsoftcoll
        
       end

        double precision function NeuQ2qx_M2RadiationGSScoll()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double complex radgSgSg, radgSnSn1, radgSnSn2
        double complex A1L, A1Lc, A1R, A1Rc
        double precision wquark
           
        NeuQ2qx_M2RadiationGSScoll = 0d0

        radgSgSg = 0d0       
        radgSnSn1 = 0d0
        radgSnSn2 = 0d0

        A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A1R = AR(chimeigen,sfmeigen,sftype,sfgen)
        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        

        wquark = 0d0

     
! colour class CA**2 CF
        radgSgSg = (CV**2*CF)*((GS**4*(A1Lc*(A1R*m1*m3*(9*m1**2*p2p3 + 9*m2**2*p2p3 + 9*m3**2*p2p3 - 2*p2p3*pap1 - 4*p1p2*pap2 - 14*p1p3*pap2 - 
     &            14*p1p2*pap3 - 4*p1p3*pap3 + 18*p2p3*papb - 2*p2p3*pbp1 - 4*p1p2*pbp2 - 14*p1p3*pbp2 - 14*p1p2*pbp3 - 
     &            4*p1p3*pbp3) + A1L*(-(m3**2*p2p3*pap1) - 2*m3**2*p1p2*pap2 - 7*m3**2*p1p3*pap2 - 7*m3**2*p1p2*pap3 - 
     &            2*m3**2*p1p3*pap3 + m2**2*(p2p3*pap1 + 2*p1p2*pap2 + 7*p1p3*pap2 + 7*p1p2*pap3 + 2*p1p3*pap3) + 
     &            18*m3**2*p2p3*papb - 2*p2p3*papb*pbp1 - 4*p1p2*papb*pbp2 - 14*p1p3*papb*pbp2 - 14*p1p2*papb*pbp3 - 
     &            4*p1p3*papb*pbp3 + m1**2*(18*m3**2*p2p3 - 2*p1p2*pap2 - 7*p1p3*pap2 - 7*p1p2*pap3 - 2*p1p3*pap3 - 
     &               p2p3*(pap1 + 2*pbp1) - 4*p1p2*pbp2 - 14*p1p3*pbp2 - 14*p1p2*pbp3 - 4*p1p3*pbp3))) + 
     &      A1Rc*(A1L*m1*m3*(9*m1**2*p2p3 + 9*m2**2*p2p3 + 9*m3**2*p2p3 - 2*p2p3*pap1 - 4*p1p2*pap2 - 14*p1p3*pap2 - 
     &            14*p1p2*pap3 - 4*p1p3*pap3 + 18*p2p3*papb - 2*p2p3*pbp1 - 4*p1p2*pbp2 - 14*p1p3*pbp2 - 14*p1p2*pbp3 - 
     &            4*p1p3*pbp3) + A1R*(-(m3**2*p2p3*pap1) - 2*m3**2*p1p2*pap2 - 7*m3**2*p1p3*pap2 - 7*m3**2*p1p2*pap3 - 
     &            2*m3**2*p1p3*pap3 + m2**2*(p2p3*pap1 + 2*p1p2*pap2 + 7*p1p3*pap2 + 7*p1p2*pap3 + 2*p1p3*pap3) + 
     &            18*m3**2*p2p3*papb - 2*p2p3*papb*pbp1 - 4*p1p2*papb*pbp2 - 14*p1p3*papb*pbp2 - 14*p1p2*papb*pbp3 - 
     &            4*p1p3*papb*pbp3 + m1**2*(18*m3**2*p2p3 - 2*p1p2*pap2 - 7*p1p3*pap2 - 7*p1p2*pap3 - 2*p1p3*pap3 - 
     &               p2p3*(pap1 + 2*pbp1) - 4*p1p2*pbp2 - 14*p1p3*pbp2 - 14*p1p2*pbp3 - 4*p1p3*pbp3)))))/
     &  (p2p3**2*(m1**2 + m2**2 - m3**2 + 2*papb)**2))
     

! ghosts           
               radgSnSn1 = (CV**2*CF)*(GS**4*p1p2*(A1L*(2*A1Rc*m1*m3*(pap2 + pbp2) + A1Lc*(-(m2**2*pap2) + m3**2*pap2 + 2*papb*pbp2 + m1**2*(pap2 + 2*pbp2))) + 
     &      A1R*(2*A1Lc*m1*m3*(pap2 + pbp2) + A1Rc*(-(m2**2*pap2) + m3**2*pap2 + 2*papb*pbp2 + m1**2*(pap2 + 2*pbp2)))))/
     &  (p2p3**2*(m1**2 + m2**2 - m3**2 + 2*papb)**2)
     
               radgSnSn2 = (CV**2*CF)*(GS**4*p1p3*(A1L*(2*A1Rc*m1*m3*(pap3 + pbp3) + A1Lc*(-(m2**2*pap3) + m3**2*pap3 + 2*papb*pbp3 + m1**2*(pap3 + 2*pbp3))) + 
     &      A1R*(2*A1Lc*m1*m3*(pap3 + pbp3) + A1Rc*(-(m2**2*pap3) + m3**2*pap3 + 2*papb*pbp3 + m1**2*(pap3 + 2*pbp3)))))/
     &  (p2p3**2*(m1**2 + m2**2 - m3**2 + 2*papb)**2)


        NeuQ2qx_M2RadiationGSScoll = 0.5*(radgSgSg) - 0.5d0*(radgSnSn1  + radgSnSn2)

  
!        print*,"NeuQ2qx_M2RadiationGSScoll",NeuQ2qx_M2RadiationGSScoll
        
       end

! =================================================================

       double precision function NeuQ2qx_M2RadiationGSTsoft()
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double complex radgSiTi, radgSpTi, radgSfTp, radgSfTf, radgSiTp, radgSiTf, radgSpTp, radgSpTf, radgSfTi
        double complex A1L, A1R, A2Lc, A2Rc
        double precision wquark


        NeuQ2qx_M2RadiationGSTsoft = 0d0
    
        radgSfTp = 0d0
        radgSfTf = 0d0
        radgSfTi = 0d0
        
        radgSiTp = 0d0
        radgSiTf = 0d0
        radgSiTi = 0d0
        
        radgSpTp = 0d0
        radgSpTf = 0d0
        radgSpTi = 0d0

        
        A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A1R = AR(chimeigen,sfmeigen,sftype,sfgen)


        wquark = 0d0

        A2Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A2Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
      

! colour class CF**2 CA

        radgSiTi = (CF**2*CV)*((GS**4*(m2**2 - pbp3)*(A1L*(A2Rc*m1*m3*(-m1**2 + m2**2 + m3**2 - 2*p1p3 + 2*pbp1 - 2*pbp3) + 
     &         A2Lc*(m2**2*pap1 + m3**2*pap1 - 2*m3**2*pap3 + 2*pap1*pap3 + 2*m3**2*papb - 2*pap1*papb - 
     &            m1**2*(2*p1p3 + pap1 - 2*pbp1) - 2*pap1*pbp3)) + 
     &      A1R*(A2Lc*m1*m3*(-m1**2 + m2**2 + m3**2 - 2*p1p3 + 2*pbp1 - 2*pbp3) + 
     &         A2Rc*(m2**2*pap1 + m3**2*pap1 - 2*m3**2*pap3 + 2*pap1*pap3 + 2*m3**2*papb - 2*pap1*papb - 
     &            m1**2*(2*p1p3 + pap1 - 2*pbp1) - 2*pap1*pbp3))))/(p1p2*(m1**2 - m2**2 + m3**2 - 2*pap1)*pbp3**2))
     
        radgSpTi = (CF**2*CV)*(-(GS**4*(A2Rc*(A1R*(-(m3**2*p2p3*pap1) + m3**2*p1p3*pap2 - 2*m3**4*pap3 - m3**2*p1p2*pap3 + 2*m3**2*p1p3*pap3 + 
     &             2*m3**2*p2p3*pap3 + 2*m3**2*pap1*pap3 + 4*m3**4*papb + 2*m3**2*p1p2*papb - 6*m3**2*p1p3*papb + 
     &             2*p1p3**2*papb - 4*m3**2*p2p3*papb + 2*p1p3*p2p3*papb - 4*m3**2*pap1*papb + 2*p1p3*pap1*papb + 
     &             2*p2p3*pap1*papb - 2*m3**2*pap2*pbp1 - 2*p1p3*pap3*pbp1 - 2*p2p3*pap3*pbp1 - 2*pap1*pap3*pbp1 + 
     &             4*m3**2*papb*pbp1 - 2*p1p3*papb*pbp1 + 2*pap3*pbp1**2 + 2*m3**2*pap1*pbp2 - 2*m3**2*pap3*pbp2 - 
     &             2*pap1*pap3*pbp2 + 4*m3**2*papb*pbp2 - 2*p1p3*papb*pbp2 + 2*pap3*pbp1*pbp2 + 
     &             m2**2*(4*m3**2*pap1 - 3*p2p3*pap1 - 4*pap1**2 - 4*pap1*pap2 - p1p3*(4*pap1 + pap2) + p1p2*pap3 - 
     &                2*p1p2*papb + 4*pap1*pbp1 + 2*pap2*pbp1 + 2*pap1*pbp2) - 2*m3**2*pap1*pbp3 + 2*p1p3*pap1*pbp3 + 
     &             2*p2p3*pap1*pbp3 + 2*pap1**2*pbp3 + 2*m3**2*pap2*pbp3 + 2*pap1*pap2*pbp3 + 2*p1p2*papb*pbp3 - 
     &             2*pap1*pbp1*pbp3 - 2*pap2*pbp1*pbp3 + 
     &             m1**2*(2*p1p3**2 + p2p3*pap1 - p1p2*pap3 + 2*p1p2*papb - 2*m3**2*(p1p3 - 2*pbp1) - 4*p2p3*pbp1 - 
     &                4*pap1*pbp1 - 2*pap2*pbp1 + 4*pbp1**2 + p1p3*(2*p2p3 + 2*pap1 + pap2 - 6*pbp1 - 2*pbp2) - 2*pap1*pbp2 + 
     &                4*pbp1*pbp2 + 2*p1p2*pbp3)) + 
     &          A1L*m1*m3*(2*p1p3**2 + m1**2*p2p3 + 2*p1p3*p2p3 + 2*p1p3*pap1 + 2*p1p3*pap2 - 2*p1p2*pap3 + 2*p1p3*pap3 + 
     &             2*p2p3*pap3 + 2*pap1*pap3 + 4*p1p2*papb - 4*p1p3*papb - 2*p2p3*papb - 4*pap1*papb - 6*p1p3*pbp1 - 
     &             4*p2p3*pbp1 - 4*pap1*pbp1 - 4*pap2*pbp1 - 2*pap3*pbp1 + 4*papb*pbp1 + 4*pbp1**2 - 2*m1**2*pbp2 - 
     &             2*p1p3*pbp2 - 4*pap3*pbp2 + 4*papb*pbp2 + 4*pbp1*pbp2 + 
     &             m2**2*(4*m3**2 - 4*p1p3 - 3*p2p3 - 4*pap1 - 4*pap2 + 4*pbp1 + 2*pbp2) + 2*p1p2*pbp3 + 2*p1p3*pbp3 + 
     &             2*p2p3*pbp3 + 2*pap1*pbp3 + 4*pap2*pbp3 - 2*pbp1*pbp3 - 
     &             m3**2*(2*p1p3 + p2p3 + 2*(pap3 - 2*papb - 2*pbp1 - pbp2 + pbp3)))) + 
     &       A2Lc*(A1L*(-(m3**2*p2p3*pap1) + m3**2*p1p3*pap2 - 2*m3**4*pap3 - m3**2*p1p2*pap3 + 2*m3**2*p1p3*pap3 + 
     &             2*m3**2*p2p3*pap3 + 2*m3**2*pap1*pap3 + 4*m3**4*papb + 2*m3**2*p1p2*papb - 6*m3**2*p1p3*papb + 
     &             2*p1p3**2*papb - 4*m3**2*p2p3*papb + 2*p1p3*p2p3*papb - 4*m3**2*pap1*papb + 2*p1p3*pap1*papb + 
     &             2*p2p3*pap1*papb - 2*m3**2*pap2*pbp1 - 2*p1p3*pap3*pbp1 - 2*p2p3*pap3*pbp1 - 2*pap1*pap3*pbp1 + 
     &             4*m3**2*papb*pbp1 - 2*p1p3*papb*pbp1 + 2*pap3*pbp1**2 + 2*m3**2*pap1*pbp2 - 2*m3**2*pap3*pbp2 - 
     &             2*pap1*pap3*pbp2 + 4*m3**2*papb*pbp2 - 2*p1p3*papb*pbp2 + 2*pap3*pbp1*pbp2 + 
     &             m2**2*(4*m3**2*pap1 - 3*p2p3*pap1 - 4*pap1**2 - 4*pap1*pap2 - p1p3*(4*pap1 + pap2) + p1p2*pap3 - 
     &                2*p1p2*papb + 4*pap1*pbp1 + 2*pap2*pbp1 + 2*pap1*pbp2) - 2*m3**2*pap1*pbp3 + 2*p1p3*pap1*pbp3 + 
     &             2*p2p3*pap1*pbp3 + 2*pap1**2*pbp3 + 2*m3**2*pap2*pbp3 + 2*pap1*pap2*pbp3 + 2*p1p2*papb*pbp3 - 
     &             2*pap1*pbp1*pbp3 - 2*pap2*pbp1*pbp3 + 
     &             m1**2*(2*p1p3**2 + p2p3*pap1 - p1p2*pap3 + 2*p1p2*papb - 2*m3**2*(p1p3 - 2*pbp1) - 4*p2p3*pbp1 - 
     &                4*pap1*pbp1 - 2*pap2*pbp1 + 4*pbp1**2 + p1p3*(2*p2p3 + 2*pap1 + pap2 - 6*pbp1 - 2*pbp2) - 2*pap1*pbp2 + 
     &                4*pbp1*pbp2 + 2*p1p2*pbp3)) + 
     &          A1R*m1*m3*(2*p1p3**2 + m1**2*p2p3 + 2*p1p3*p2p3 + 2*p1p3*pap1 + 2*p1p3*pap2 - 2*p1p2*pap3 + 2*p1p3*pap3 + 
     &             2*p2p3*pap3 + 2*pap1*pap3 + 4*p1p2*papb - 4*p1p3*papb - 2*p2p3*papb - 4*pap1*papb - 6*p1p3*pbp1 - 
     &             4*p2p3*pbp1 - 4*pap1*pbp1 - 4*pap2*pbp1 - 2*pap3*pbp1 + 4*papb*pbp1 + 4*pbp1**2 - 2*m1**2*pbp2 - 
     &             2*p1p3*pbp2 - 4*pap3*pbp2 + 4*papb*pbp2 + 4*pbp1*pbp2 + 
     &             m2**2*(4*m3**2 - 4*p1p3 - 3*p2p3 - 4*pap1 - 4*pap2 + 4*pbp1 + 2*pbp2) + 2*p1p2*pbp3 + 2*p1p3*pbp3 + 
     &             2*p2p3*pbp3 + 2*pap1*pbp3 + 4*pap2*pbp3 - 2*pbp1*pbp3 - 
     &             m3**2*(2*p1p3 + p2p3 + 2*(pap3 - 2*papb - 2*pbp1 - pbp2 + pbp3))))))/
     &  (2.*p1p2*(m1**2 - m2**2 + m3**2 - 2*pap1)*(m1**2 + m2**2 - m3**2 + 2*papb)*pbp3))
     
        radgSfTp = (CF**2*CV)*(-(GS**4*(A2Lc*(A1R*m1*m3*(2*p1p2**2 + m1**2*p2p3 + 2*p1p2*p2p3 + 2*p1p2*pap1 + 2*p1p2*pap2 - 2*p1p3*pap2 + 
     &             2*p2p3*pap2 + 2*pap1*pap2 + 2*p1p2*pap3 - 4*p1p2*papb + 4*p1p3*papb - 2*p2p3*papb - 4*pap1*papb - 
     &             6*p1p2*pbp1 - 4*p2p3*pbp1 - 4*pap1*pbp1 - 2*pap2*pbp1 - 4*pap3*pbp1 + 4*papb*pbp1 + 4*pbp1**2 + 
     &             2*p1p2*pbp2 + 2*p1p3*pbp2 + 2*p2p3*pbp2 + 2*pap1*pbp2 + 4*pap3*pbp2 - 2*pbp1*pbp2 - 
     &             m3**2*(2*p1p2 + p2p3 + 2*(pap2 - 2*papb - 2*pbp1 + pbp2 - pbp3)) - 2*m1**2*pbp3 - 2*p1p2*pbp3 - 
     &             4*pap2*pbp3 + 4*papb*pbp3 + 4*pbp1*pbp3 + 
     &             m2**2*(4*m3**2 - 4*p1p2 - 3*p2p3 - 4*pap1 - 4*pap3 + 4*pbp1 + 2*pbp3)) + 
     &          A1L*(-(m3**2*p2p3*pap1) - 2*m3**4*pap2 + 2*m3**2*p1p2*pap2 - m3**2*p1p3*pap2 + 2*m3**2*p2p3*pap2 + 
     &             2*m3**2*pap1*pap2 + m3**2*p1p2*pap3 + 4*m3**4*papb - 6*m3**2*p1p2*papb + 2*p1p2**2*papb + 
     &             2*m3**2*p1p3*papb - 4*m3**2*p2p3*papb + 2*p1p2*p2p3*papb - 4*m3**2*pap1*papb + 2*p1p2*pap1*papb + 
     &             2*p2p3*pap1*papb - 2*p1p2*pap2*pbp1 - 2*p2p3*pap2*pbp1 - 2*pap1*pap2*pbp1 - 2*m3**2*pap3*pbp1 + 
     &             4*m3**2*papb*pbp1 - 2*p1p2*papb*pbp1 + 2*pap2*pbp1**2 - 2*m3**2*pap1*pbp2 + 2*p1p2*pap1*pbp2 + 
     &             2*p2p3*pap1*pbp2 + 2*pap1**2*pbp2 + 2*m3**2*pap3*pbp2 + 2*pap1*pap3*pbp2 + 2*p1p3*papb*pbp2 - 
     &             2*pap1*pbp1*pbp2 - 2*pap3*pbp1*pbp2 + 2*m3**2*pap1*pbp3 - 2*m3**2*pap2*pbp3 - 2*pap1*pap2*pbp3 + 
     &             4*m3**2*papb*pbp3 - 2*p1p2*papb*pbp3 + 2*pap2*pbp1*pbp3 + 
     &             m2**2*(4*m3**2*pap1 - 3*p2p3*pap1 - 4*pap1**2 + p1p3*pap2 - 4*pap1*pap3 - p1p2*(4*pap1 + pap3) - 
     &                2*p1p3*papb + 4*pap1*pbp1 + 2*pap3*pbp1 + 2*pap1*pbp3) + 
     &             m1**2*(2*p1p2**2 + p2p3*pap1 - p1p3*pap2 + 2*p1p3*papb - 2*m3**2*(p1p2 - 2*pbp1) - 4*p2p3*pbp1 - 
     &                4*pap1*pbp1 - 2*pap3*pbp1 + 4*pbp1**2 + 2*p1p3*pbp2 + p1p2*(2*p2p3 + 2*pap1 + pap3 - 6*pbp1 - 2*pbp3) - 
     &                2*pap1*pbp3 + 4*pbp1*pbp3))) + 
     &       A2Rc*(A1L*m1*m3*(2*p1p2**2 + m1**2*p2p3 + 2*p1p2*p2p3 + 2*p1p2*pap1 + 2*p1p2*pap2 - 2*p1p3*pap2 + 2*p2p3*pap2 + 
     &             2*pap1*pap2 + 2*p1p2*pap3 - 4*p1p2*papb + 4*p1p3*papb - 2*p2p3*papb - 4*pap1*papb - 6*p1p2*pbp1 - 
     &             4*p2p3*pbp1 - 4*pap1*pbp1 - 2*pap2*pbp1 - 4*pap3*pbp1 + 4*papb*pbp1 + 4*pbp1**2 + 2*p1p2*pbp2 + 
     &             2*p1p3*pbp2 + 2*p2p3*pbp2 + 2*pap1*pbp2 + 4*pap3*pbp2 - 2*pbp1*pbp2 - 
     &             m3**2*(2*p1p2 + p2p3 + 2*(pap2 - 2*papb - 2*pbp1 + pbp2 - pbp3)) - 2*m1**2*pbp3 - 2*p1p2*pbp3 - 
     &             4*pap2*pbp3 + 4*papb*pbp3 + 4*pbp1*pbp3 + 
     &             m2**2*(4*m3**2 - 4*p1p2 - 3*p2p3 - 4*pap1 - 4*pap3 + 4*pbp1 + 2*pbp3)) + 
     &          A1R*(-(m3**2*p2p3*pap1) - 2*m3**4*pap2 + 2*m3**2*p1p2*pap2 - m3**2*p1p3*pap2 + 2*m3**2*p2p3*pap2 + 
     &             2*m3**2*pap1*pap2 + m3**2*p1p2*pap3 + 4*m3**4*papb - 6*m3**2*p1p2*papb + 2*p1p2**2*papb + 
     &             2*m3**2*p1p3*papb - 4*m3**2*p2p3*papb + 2*p1p2*p2p3*papb - 4*m3**2*pap1*papb + 2*p1p2*pap1*papb + 
     &             2*p2p3*pap1*papb - 2*p1p2*pap2*pbp1 - 2*p2p3*pap2*pbp1 - 2*pap1*pap2*pbp1 - 2*m3**2*pap3*pbp1 + 
     &             4*m3**2*papb*pbp1 - 2*p1p2*papb*pbp1 + 2*pap2*pbp1**2 - 2*m3**2*pap1*pbp2 + 2*p1p2*pap1*pbp2 + 
     &             2*p2p3*pap1*pbp2 + 2*pap1**2*pbp2 + 2*m3**2*pap3*pbp2 + 2*pap1*pap3*pbp2 + 2*p1p3*papb*pbp2 - 
     &             2*pap1*pbp1*pbp2 - 2*pap3*pbp1*pbp2 + 2*m3**2*pap1*pbp3 - 2*m3**2*pap2*pbp3 - 2*pap1*pap2*pbp3 + 
     &             4*m3**2*papb*pbp3 - 2*p1p2*papb*pbp3 + 2*pap2*pbp1*pbp3 + 
     &             m2**2*(4*m3**2*pap1 - 3*p2p3*pap1 - 4*pap1**2 + p1p3*pap2 - 4*pap1*pap3 - p1p2*(4*pap1 + pap3) - 
     &                2*p1p3*papb + 4*pap1*pbp1 + 2*pap3*pbp1 + 2*pap1*pbp3) + 
     &             m1**2*(2*p1p2**2 + p2p3*pap1 - p1p3*pap2 + 2*p1p3*papb - 2*m3**2*(p1p2 - 2*pbp1) - 4*p2p3*pbp1 - 
     &                4*pap1*pbp1 - 2*pap3*pbp1 + 4*pbp1**2 + 2*p1p3*pbp2 + p1p2*(2*p2p3 + 2*pap1 + pap3 - 6*pbp1 - 2*pbp3) - 
     &                2*pap1*pbp3 + 4*pbp1*pbp3)))))/
     &  (2.*p1p3*(m1**2 - m2**2 + m3**2 - 2*pap1)*(m1**2 + m2**2 - m3**2 + 2*papb)*pbp2))
     
        radgSfTf = (CF**2*CV)*(-((GS**4*(A2Lc*(A1R*m1*m3*(2*m2**2*(m3**2 + p1p3) + p1p3*(p2p3 - pap2 + 2*papb - pbp2 - 2*pbp3) - 
     &              m3**2*(p1p2 + p2p3 + pap2 - 2*papb - 2*pbp1 + pbp2 - 2*pbp3)) - 
     &           A1L*(m3**4*pap2 + m3**2*p1p3*pap2 - 2*m2**2*(-(p1p3*pap3) + m3**2*(pap1 + pap3)) - 2*m3**4*papb + 
     &              m3**2*p1p2*papb - 2*m3**2*p1p3*papb + m3**2*p2p3*papb - p1p3*p2p3*papb - m3**2*pap2*pbp1 + 
     &              m3**2*pap1*pbp2 + m3**2*pap3*pbp2 - p1p3*pap3*pbp2 - m3**2*pap2*pbp3 + p1p3*pap2*pbp3 + 
     &              m1**2*(-(p1p3*(p2p3 - 2*pbp3)) + m3**2*(p1p2 + p2p3 - 2*(pbp1 + pbp3))))) + 
     &        A2Rc*(A1L*m1*m3*(2*m2**2*(m3**2 + p1p3) + p1p3*(p2p3 - pap2 + 2*papb - pbp2 - 2*pbp3) - 
     &              m3**2*(p1p2 + p2p3 + pap2 - 2*papb - 2*pbp1 + pbp2 - 2*pbp3)) - 
     &           A1R*(m3**4*pap2 + m3**2*p1p3*pap2 - 2*m2**2*(-(p1p3*pap3) + m3**2*(pap1 + pap3)) - 2*m3**4*papb + 
     &              m3**2*p1p2*papb - 2*m3**2*p1p3*papb + m3**2*p2p3*papb - p1p3*p2p3*papb - m3**2*pap2*pbp1 + 
     &              m3**2*pap1*pbp2 + m3**2*pap3*pbp2 - p1p3*pap3*pbp2 - m3**2*pap2*pbp3 + p1p3*pap2*pbp3 + 
     &              m1**2*(-(p1p3*(p2p3 - 2*pbp3)) + m3**2*(p1p2 + p2p3 - 2*(pbp1 + pbp3)))))))/
     &    (p1p3**2*(m1**2 + m2**2 - m3**2 + 2*papb)*pbp2)))
     
! colour class -CF/2
        radgSiTp = (-CF/2d0)*(-(GS**4*(2*m2**2 - p1p3 + p2p3 + pap3 - 2*papb + 2*pbp1 - 2*pbp2 - pbp3)*
     &     (A1L*(A2Rc*m1*m3*(-2*m2**2 + p1p2 - p2p3 + pap2 - 2*papb - 2*pbp1 + pbp2 + 2*pbp3) + 
     &          A2Lc*(-2*m2**2*pap1 - p2p3*pap1 + m3**2*pap2 + p1p3*pap2 - p1p2*pap3 - 2*m3**2*papb + p1p2*papb - 
     &             2*p1p3*papb + m1**2*(p1p2 - 2*pbp1) - pap2*pbp1 + 2*pap3*pbp1 + pap1*pbp2 + 2*pap1*pbp3)) + 
     &       A1R*(A2Lc*m1*m3*(-2*m2**2 + p1p2 - p2p3 + pap2 - 2*papb - 2*pbp1 + pbp2 + 2*pbp3) + 
     &          A2Rc*(-2*m2**2*pap1 - p2p3*pap1 + m3**2*pap2 + p1p3*pap2 - p1p2*pap3 - 2*m3**2*papb + p1p2*papb - 
     &             2*p1p3*papb + m1**2*(p1p2 - 2*pbp1) - pap2*pbp1 + 2*pap3*pbp1 + pap1*pbp2 + 2*pap1*pbp3))))/
     &  (4.*p1p2*(m1**2 - m2**2 + m3**2 - 2*pap1)*pbp2*pbp3))
     
        radgSiTf = (-CF/2d0)*(-(GS**4*(A1R*(A2Rc*(p1p3*p2p3*pap1 - m3**2*p1p3*pap2 - p1p3**2*pap2 + p1p2*p1p3*pap3 + 2*m3**2*p1p3*papb - 
     &             p1p2*p1p3*papb + 2*p1p3**2*papb - m3**2*p2p3*papb - 2*p2p3*pap1*pbp1 + 2*m3**2*pap2*pbp1 + 
     &             3*p1p3*pap2*pbp1 - 2*p1p2*pap3*pbp1 - 2*p1p3*pap3*pbp1 - 2*p2p3*pap3*pbp1 - 4*m3**2*papb*pbp1 + 
     &             2*p1p2*papb*pbp1 - 4*p1p3*papb*pbp1 + 2*p2p3*papb*pbp1 - 2*pap2*pbp1**2 + 4*pap3*pbp1**2 - 
     &             p1p3*pap1*pbp2 + m3**2*pap3*pbp2 + 2*p1p3*pap3*pbp2 - 2*p1p3*papb*pbp2 + 2*pap1*pbp1*pbp2 - 
     &             2*p1p3*pap1*pbp3 + m3**2*pap2*pbp3 - 2*p1p2*pap3*pbp3 + 4*pap1*pbp1*pbp3 - 2*pap2*pbp1*pbp3 + 
     &             8*pap3*pbp1*pbp3 + 2*pap1*pbp2*pbp3 + 
     &             m2**2*(-(p2p3*pap1) - 2*m3**2*pap3 + p1p2*pap3 + p1p3*(2*pap1 + pap2 - 4*pap3 + 2*papb) - 4*pap1*pbp1 - 
     &                2*pap3*pbp1 - 2*pap1*pbp3) + 
     &             m1**2*(2*m2**2*p1p3 + 2*p1p3*pbp1 + p2p3*pbp1 - 4*pbp1**2 - p1p3*pbp2 - 4*pbp1*pbp3 + 
     &                p1p2*(-p1p3 + 2*pbp1 + pbp3))) + 
     &          A2Lc*m1*m3*(-(p1p2*p1p3) + p1p3*p2p3 - p1p3*pap2 + 2*p1p3*papb - p2p3*papb + 2*p1p2*pbp1 + 2*p1p3*pbp1 - 
     &             p2p3*pbp1 + 2*pap2*pbp1 - 4*papb*pbp1 - 4*pbp1**2 - 2*p1p3*pbp2 + pap3*pbp2 + 2*pbp1*pbp2 + p1p2*pbp3 - 
     &             2*p1p3*pbp3 + pap2*pbp3 + 2*pbp2*pbp3 + m2**2*(4*p1p3 - p2p3 - 2*(pap3 + 2*pbp1 + pbp3)))) + 
     &       A1L*(A2Lc*(p1p3*p2p3*pap1 - m3**2*p1p3*pap2 - p1p3**2*pap2 + p1p2*p1p3*pap3 + 2*m3**2*p1p3*papb - 
     &             p1p2*p1p3*papb + 2*p1p3**2*papb - m3**2*p2p3*papb - 2*p2p3*pap1*pbp1 + 2*m3**2*pap2*pbp1 + 
     &             3*p1p3*pap2*pbp1 - 2*p1p2*pap3*pbp1 - 2*p1p3*pap3*pbp1 - 2*p2p3*pap3*pbp1 - 4*m3**2*papb*pbp1 + 
     &             2*p1p2*papb*pbp1 - 4*p1p3*papb*pbp1 + 2*p2p3*papb*pbp1 - 2*pap2*pbp1**2 + 4*pap3*pbp1**2 - 
     &             p1p3*pap1*pbp2 + m3**2*pap3*pbp2 + 2*p1p3*pap3*pbp2 - 2*p1p3*papb*pbp2 + 2*pap1*pbp1*pbp2 - 
     &             2*p1p3*pap1*pbp3 + m3**2*pap2*pbp3 - 2*p1p2*pap3*pbp3 + 4*pap1*pbp1*pbp3 - 2*pap2*pbp1*pbp3 + 
     &             8*pap3*pbp1*pbp3 + 2*pap1*pbp2*pbp3 + 
     &             m2**2*(-(p2p3*pap1) - 2*m3**2*pap3 + p1p2*pap3 + p1p3*(2*pap1 + pap2 - 4*pap3 + 2*papb) - 4*pap1*pbp1 - 
     &                2*pap3*pbp1 - 2*pap1*pbp3) + 
     &             m1**2*(2*m2**2*p1p3 + 2*p1p3*pbp1 + p2p3*pbp1 - 4*pbp1**2 - p1p3*pbp2 - 4*pbp1*pbp3 + 
     &                p1p2*(-p1p3 + 2*pbp1 + pbp3))) + 
     &          A2Rc*m1*m3*(-(p1p2*p1p3) + p1p3*p2p3 - p1p3*pap2 + 2*p1p3*papb - p2p3*papb + 2*p1p2*pbp1 + 2*p1p3*pbp1 - 
     &             p2p3*pbp1 + 2*pap2*pbp1 - 4*papb*pbp1 - 4*pbp1**2 - 2*p1p3*pbp2 + pap3*pbp2 + 2*pbp1*pbp2 + p1p2*pbp3 - 
     &             2*p1p3*pbp3 + pap2*pbp3 + 2*pbp2*pbp3 + m2**2*(4*p1p3 - p2p3 - 2*(pap3 + 2*pbp1 + pbp3))))))/
     &  (4.*p1p2*p1p3*pbp2*pbp3))
     
        radgSpTp = (-CF/2d0)*(-((GS**4*(A1L*(A2Rc*m1*m3*(-(m3**2*p1p2) + p1p2**2 + p1p2*pap2 + 2*p1p2*papb + 2*m3**2*pbp1 - 4*p1p2*pbp1 - 
     &              4*pap2*pbp1 + 4*pbp1**2 + m1**2*(p1p2 - 2*pbp1 - pbp2) + m3**2*pbp2 + p1p2*pbp2 + 2*papb*pbp2 + 
     &              m2**2*(-p1p2 - 2*pap2 + 2*pbp1 + pbp2)) + 
     &           A2Lc*(-(m3**2*p1p2*pap1) + m3**2*p1p2*pap2 - m3**2*p1p2*papb + p1p2**2*papb + 2*p1p2*pap1*papb + 
     &              2*p1p2*papb**2 + 2*m3**2*pap1*pbp1 - 3*m3**2*pap2*pbp1 - p1p2*pap2*pbp1 + 4*m3**2*papb*pbp1 - 
     &              4*p1p2*papb*pbp1 - 4*pap1*papb*pbp1 - 2*pap2*papb*pbp1 + 4*pap2*pbp1**2 + m3**2*pap1*pbp2 + 
     &              p1p2*pap1*pbp2 + 4*m3**2*papb*pbp2 - 2*pap1*papb*pbp2 - 4*pap1*pbp1*pbp2 + 
     &              m2**2*(-2*m3**2*pap2 + p1p2*(pap1 - papb) + 2*pap1*pbp1 + pap2*pbp1 + pap1*pbp2) + 
     &              m1**2*(-2*m2**2*p1p2 + p1p2**2 + p1p2*(pap1 + papb - 4*pbp1) - 2*pap1*pbp1 - pap2*pbp1 + 4*pbp1**2 - 
     &                 pap1*pbp2 + 4*pbp1*pbp2))) + 
     &        A1R*(A2Lc*m1*m3*(-(m3**2*p1p2) + p1p2**2 + p1p2*pap2 + 2*p1p2*papb + 2*m3**2*pbp1 - 4*p1p2*pbp1 - 4*pap2*pbp1 + 
     &              4*pbp1**2 + m1**2*(p1p2 - 2*pbp1 - pbp2) + m3**2*pbp2 + p1p2*pbp2 + 2*papb*pbp2 + 
     &              m2**2*(-p1p2 - 2*pap2 + 2*pbp1 + pbp2)) + 
     &           A2Rc*(-(m3**2*p1p2*pap1) + m3**2*p1p2*pap2 - m3**2*p1p2*papb + p1p2**2*papb + 2*p1p2*pap1*papb + 
     &              2*p1p2*papb**2 + 2*m3**2*pap1*pbp1 - 3*m3**2*pap2*pbp1 - p1p2*pap2*pbp1 + 4*m3**2*papb*pbp1 - 
     &              4*p1p2*papb*pbp1 - 4*pap1*papb*pbp1 - 2*pap2*papb*pbp1 + 4*pap2*pbp1**2 + m3**2*pap1*pbp2 + 
     &              p1p2*pap1*pbp2 + 4*m3**2*papb*pbp2 - 2*pap1*papb*pbp2 - 4*pap1*pbp1*pbp2 + 
     &              m2**2*(-2*m3**2*pap2 + p1p2*(pap1 - papb) + 2*pap1*pbp1 + pap2*pbp1 + pap1*pbp2) + 
     &              m1**2*(-2*m2**2*p1p2 + p1p2**2 + p1p2*(pap1 + papb - 4*pbp1) - 2*pap1*pbp1 - pap2*pbp1 + 4*pbp1**2 - 
     &                 pap1*pbp2 + 4*pbp1*pbp2)))))/
     &    (p1p2*(m1**2 - m2**2 + m3**2 - 2*pap1)*(m1**2 + m2**2 - m3**2 + 2*papb)*pbp2)))
     
        radgSpTf = (-CF/2d0)*(-((GS**4*(A2Lc*(A1R*m1*m3*(m2**2*p1p2 + p1p2*p1p3 - p1p2*pap1 - 2*p1p2*pap3 + p1p2*papb - p1p2*pbp1 - 2*p1p3*pbp1 - 
     &              p2p3*pbp1 + 2*pap1*pbp1 - pap2*pbp1 + 4*pap3*pbp1 + 2*pbp1**2 - p1p3*pbp2 + pap1*pbp2 + 2*pap3*pbp2 + 
     &              m3**2*(-p1p2 + 2*pbp1 + pbp2) - p1p2*pbp3 + 4*pbp1*pbp3 + 2*pbp2*pbp3) + 
     &           A1L*(-(m3**2*p1p2*pap1) - 2*m3**2*p1p2*pap3 + 
     &              m2**2*(p2p3*pap1 - (m3**2 + p1p3)*pap2 + p1p2*(2*pap1 + pap3)) + p1p2*p1p3*papb + 2*m3**2*pap1*pbp1 - 
     &              m3**2*pap2*pbp1 + 4*m3**2*pap3*pbp1 + p1p2*pap3*pbp1 + 2*m3**2*papb*pbp1 - 2*p1p2*papb*pbp1 - 
     &              2*p1p3*papb*pbp1 - 2*p2p3*papb*pbp1 + 2*pap2*pbp1**2 - 2*pap3*pbp1**2 + m3**2*pap1*pbp2 + 
     &              2*m3**2*pap3*pbp2 + 2*m3**2*papb*pbp2 - 2*pap1*pbp1*pbp2 - 2*pap3*pbp1*pbp2 - p1p2*pap1*pbp3 + 
     &              2*pap1*pbp1*pbp3 + 2*pap2*pbp1*pbp3 + 
     &              m1**2*(-2*p1p3*pbp1 - p2p3*pbp1 - p1p3*pbp2 + m3**2*(-p1p2 + 2*pbp1 + pbp2) + p1p2*(p1p3 + pbp3)))) + 
     &        A2Rc*(A1L*m1*m3*(m2**2*p1p2 + p1p2*p1p3 - p1p2*pap1 - 2*p1p2*pap3 + p1p2*papb - p1p2*pbp1 - 2*p1p3*pbp1 - 
     &              p2p3*pbp1 + 2*pap1*pbp1 - pap2*pbp1 + 4*pap3*pbp1 + 2*pbp1**2 - p1p3*pbp2 + pap1*pbp2 + 2*pap3*pbp2 + 
     &              m3**2*(-p1p2 + 2*pbp1 + pbp2) - p1p2*pbp3 + 4*pbp1*pbp3 + 2*pbp2*pbp3) + 
     &           A1R*(-(m3**2*p1p2*pap1) - 2*m3**2*p1p2*pap3 + 
     &              m2**2*(p2p3*pap1 - (m3**2 + p1p3)*pap2 + p1p2*(2*pap1 + pap3)) + p1p2*p1p3*papb + 2*m3**2*pap1*pbp1 - 
     &              m3**2*pap2*pbp1 + 4*m3**2*pap3*pbp1 + p1p2*pap3*pbp1 + 2*m3**2*papb*pbp1 - 2*p1p2*papb*pbp1 - 
     &              2*p1p3*papb*pbp1 - 2*p2p3*papb*pbp1 + 2*pap2*pbp1**2 - 2*pap3*pbp1**2 + m3**2*pap1*pbp2 + 
     &              2*m3**2*pap3*pbp2 + 2*m3**2*papb*pbp2 - 2*pap1*pbp1*pbp2 - 2*pap3*pbp1*pbp2 - p1p2*pap1*pbp3 + 
     &              2*pap1*pbp1*pbp3 + 2*pap2*pbp1*pbp3 + 
     &              m1**2*(-2*p1p3*pbp1 - p2p3*pbp1 - p1p3*pbp2 + m3**2*(-p1p2 + 2*pbp1 + pbp2) + p1p2*(p1p3 + pbp3))))))/
     &    (p1p2*p1p3*(m1**2 + m2**2 - m3**2 + 2*papb)*pbp2)))
     
        radgSfTi = (-CF/2d0)*(-((GS**4*(A1L*(A2Rc*m1*m3*(-(m3**2*p1p3) + p1p3**2 + p1p3*pap3 + 2*p1p3*papb + 2*m3**2*pbp1 - 4*p1p3*pbp1 - 
     &              4*pap3*pbp1 + 4*pbp1**2 + m1**2*(p1p3 - 2*pbp1 - pbp3) + m3**2*pbp3 + p1p3*pbp3 + 2*papb*pbp3 + 
     &              m2**2*(-p1p3 - 2*pap3 + 2*pbp1 + pbp3)) + 
     &           A2Lc*(-(m3**2*p1p3*pap1) + m3**2*p1p3*pap3 - m3**2*p1p3*papb + p1p3**2*papb + 2*p1p3*pap1*papb + 
     &              2*p1p3*papb**2 + 2*m3**2*pap1*pbp1 - 3*m3**2*pap3*pbp1 - p1p3*pap3*pbp1 + 4*m3**2*papb*pbp1 - 
     &              4*p1p3*papb*pbp1 - 4*pap1*papb*pbp1 - 2*pap3*papb*pbp1 + 4*pap3*pbp1**2 + m3**2*pap1*pbp3 + 
     &              p1p3*pap1*pbp3 + 4*m3**2*papb*pbp3 - 2*pap1*papb*pbp3 - 4*pap1*pbp1*pbp3 + 
     &              m2**2*(-2*m3**2*pap3 + p1p3*(pap1 - papb) + 2*pap1*pbp1 + pap3*pbp1 + pap1*pbp3) + 
     &              m1**2*(-2*m2**2*p1p3 + p1p3**2 + p1p3*(pap1 + papb - 4*pbp1) - 2*pap1*pbp1 - pap3*pbp1 + 4*pbp1**2 - 
     &                 pap1*pbp3 + 4*pbp1*pbp3))) + 
     &        A1R*(A2Lc*m1*m3*(-(m3**2*p1p3) + p1p3**2 + p1p3*pap3 + 2*p1p3*papb + 2*m3**2*pbp1 - 4*p1p3*pbp1 - 4*pap3*pbp1 + 
     &              4*pbp1**2 + m1**2*(p1p3 - 2*pbp1 - pbp3) + m3**2*pbp3 + p1p3*pbp3 + 2*papb*pbp3 + 
     &              m2**2*(-p1p3 - 2*pap3 + 2*pbp1 + pbp3)) + 
     &           A2Rc*(-(m3**2*p1p3*pap1) + m3**2*p1p3*pap3 - m3**2*p1p3*papb + p1p3**2*papb + 2*p1p3*pap1*papb + 
     &              2*p1p3*papb**2 + 2*m3**2*pap1*pbp1 - 3*m3**2*pap3*pbp1 - p1p3*pap3*pbp1 + 4*m3**2*papb*pbp1 - 
     &              4*p1p3*papb*pbp1 - 4*pap1*papb*pbp1 - 2*pap3*papb*pbp1 + 4*pap3*pbp1**2 + m3**2*pap1*pbp3 + 
     &              p1p3*pap1*pbp3 + 4*m3**2*papb*pbp3 - 2*pap1*papb*pbp3 - 4*pap1*pbp1*pbp3 + 
     &              m2**2*(-2*m3**2*pap3 + p1p3*(pap1 - papb) + 2*pap1*pbp1 + pap3*pbp1 + pap1*pbp3) + 
     &              m1**2*(-2*m2**2*p1p3 + p1p3**2 + p1p3*(pap1 + papb - 4*pbp1) - 2*pap1*pbp1 - pap3*pbp1 + 4*pbp1**2 - 
     &                 pap1*pbp3 + 4*pbp1*pbp3)))))/
     &    (p1p3*(m1**2 - m2**2 + m3**2 - 2*pap1)*(m1**2 + m2**2 - m3**2 + 2*papb)*pbp3)))

   
                  
        NeuQ2qx_M2RadiationGSTsoft = 2d0*(0.5*(radgSiTi + radgSpTi + radgSfTp + radgSfTf + radgSiTp + radgSiTf + radgSpTp + radgSpTf + radgSfTi))


!         print*,"NeuQ2qx_M2RadiationGSTsoft",NeuQ2qx_M2RadiationGSTsoft
        

    
        end

        double precision function NeuQ2qx_M2RadiationGSTsoftcoll()
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double complex radgSfTg, radgSiTg, radgSpTg
        double complex radgSgTp, radgSgTf, radgSgTi
        double complex A1L, A1R, A2Lc, A2Rc
        double precision wquark


        NeuQ2qx_M2RadiationGSTsoftcoll = 0d0

        
        radgSiTg = 0d0
        radgSpTg = 0d0
        radgSfTg = 0d0
        
        radgSgTi = 0d0
        radgSgTp = 0d0
        radgSgTf = 0d0
        


        A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A1R = AR(chimeigen,sfmeigen,sftype,sfgen)


        wquark = 0d0

        A2Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A2Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))

     
! colour class -CF/2
     
     
! colour class -1/2 CA**2 CF
        radgSfTg = (-CV**2*CF/2d0)*((GS**4*(A2Lc*(A1R*m1*m3*(-p1p2**2 + p1p3**2 + m2**2*(4*p1p2 + 2*p1p3 + p2p3) - p2p3*pap2 + p1p3*pap3 + 2*p1p3*papb + 
     &            p2p3*papb - 2*p1p3*pbp1 + p2p3*pbp1 - 2*pap2*pbp1 - 4*pap3*pbp1 - 2*m3**2*pbp2 + p1p3*pbp2 - p2p3*pbp2 - 
     &            2*pap1*pbp2 + pap3*pbp2 - 4*pbp1*pbp2 - p1p2*(p2p3 + pap2 - 4*papb - 2*pbp1 + pbp2 - pbp3) + 2*m3**2*pbp3 - 
     &            p1p3*pbp3 + 2*pap1*pbp3 + pap2*pbp3 - 2*pap3*pbp3 - 2*pbp1*pbp3 + 2*pbp2*pbp3 - 2*pbp3**2) + 
     &         A1L*(-(m3**2*p1p2*pap2) - m3**2*p2p3*pap2 + m3**2*p1p3*pap3 + 
     &            m2**2*(4*p1p2*pap1 + 2*p1p3*pap1 + p2p3*pap1 - 3*p1p3*pap2 + 3*p1p2*pap3) + 4*m3**2*p1p2*papb - 
     &            p1p2**2*papb + 2*m3**2*p1p3*papb + p1p3**2*papb + m3**2*p2p3*papb - p1p2*p2p3*papb - 2*m3**2*pap2*pbp1 + 
     &            p1p2*pap2*pbp1 + p2p3*pap2*pbp1 - 4*m3**2*pap3*pbp1 - p1p3*pap3*pbp1 - 2*p1p2*papb*pbp1 - 
     &            4*p1p3*papb*pbp1 + 2*pap2*pbp1**2 + 4*pap3*pbp1**2 - 2*m3**2*pap1*pbp2 - p1p2*pap1*pbp2 - p2p3*pap1*pbp2 + 
     &            m3**2*pap3*pbp2 - 2*m3**2*papb*pbp2 + 4*p1p3*papb*pbp2 - 2*pap1*pbp1*pbp2 - 4*pap3*pbp1*pbp2 + 
     &            2*m3**2*pap1*pbp3 + p1p3*pap1*pbp3 + m3**2*pap2*pbp3 - 2*m3**2*pap3*pbp3 + 2*m3**2*papb*pbp3 - 
     &            2*p1p2*papb*pbp3 - 2*p1p3*papb*pbp3 - 4*pap1*pbp1*pbp3 + 2*pap2*pbp1*pbp3 + 2*pap3*pbp1*pbp3 + 
     &            2*pap1*pbp2*pbp3 - 2*pap1*pbp3**2 + 
     &            m1**2*(-p1p2**2 + p1p3**2 + p2p3*pbp1 - 2*m3**2*pbp2 + p1p3*(-2*pbp1 + pbp2 - 2*pbp3) + 2*m3**2*pbp3 + 
     &               p1p2*(-p2p3 + 2*pbp1 + pbp3)))) + 
     &      A2Rc*(A1L*m1*m3*(-p1p2**2 + p1p3**2 + m2**2*(4*p1p2 + 2*p1p3 + p2p3) - p2p3*pap2 + p1p3*pap3 + 2*p1p3*papb + 
     &            p2p3*papb - 2*p1p3*pbp1 + p2p3*pbp1 - 2*pap2*pbp1 - 4*pap3*pbp1 - 2*m3**2*pbp2 + p1p3*pbp2 - p2p3*pbp2 - 
     &            2*pap1*pbp2 + pap3*pbp2 - 4*pbp1*pbp2 - p1p2*(p2p3 + pap2 - 4*papb - 2*pbp1 + pbp2 - pbp3) + 2*m3**2*pbp3 - 
     &            p1p3*pbp3 + 2*pap1*pbp3 + pap2*pbp3 - 2*pap3*pbp3 - 2*pbp1*pbp3 + 2*pbp2*pbp3 - 2*pbp3**2) + 
     &         A1R*(-(m3**2*p1p2*pap2) - m3**2*p2p3*pap2 + m3**2*p1p3*pap3 + 
     &            m2**2*(4*p1p2*pap1 + 2*p1p3*pap1 + p2p3*pap1 - 3*p1p3*pap2 + 3*p1p2*pap3) + 4*m3**2*p1p2*papb - 
     &            p1p2**2*papb + 2*m3**2*p1p3*papb + p1p3**2*papb + m3**2*p2p3*papb - p1p2*p2p3*papb - 2*m3**2*pap2*pbp1 + 
     &            p1p2*pap2*pbp1 + p2p3*pap2*pbp1 - 4*m3**2*pap3*pbp1 - p1p3*pap3*pbp1 - 2*p1p2*papb*pbp1 - 
     &            4*p1p3*papb*pbp1 + 2*pap2*pbp1**2 + 4*pap3*pbp1**2 - 2*m3**2*pap1*pbp2 - p1p2*pap1*pbp2 - p2p3*pap1*pbp2 + 
     &            m3**2*pap3*pbp2 - 2*m3**2*papb*pbp2 + 4*p1p3*papb*pbp2 - 2*pap1*pbp1*pbp2 - 4*pap3*pbp1*pbp2 + 
     &            2*m3**2*pap1*pbp3 + p1p3*pap1*pbp3 + m3**2*pap2*pbp3 - 2*m3**2*pap3*pbp3 + 2*m3**2*papb*pbp3 - 
     &            2*p1p2*papb*pbp3 - 2*p1p3*papb*pbp3 - 4*pap1*pbp1*pbp3 + 2*pap2*pbp1*pbp3 + 2*pap3*pbp1*pbp3 + 
     &            2*pap1*pbp2*pbp3 - 2*pap1*pbp3**2 + 
     &            m1**2*(-p1p2**2 + p1p3**2 + p2p3*pbp1 - 2*m3**2*pbp2 + p1p3*(-2*pbp1 + pbp2 - 2*pbp3) + 2*m3**2*pbp3 + 
     &               p1p2*(-p2p3 + 2*pbp1 + pbp3))))))/
     &  (p1p3*p2p3*(m1**2 - m2**2 + m3**2 - 2*pap1)*(m1**2 + m2**2 - m3**2 + 2*papb)))
     
        radgSgTp = (-CV**2*CF/2d0)*( -(GS**4*((A1L*(A2Rc*m1*m3*(2*m2**2 - p1p2 - pap2 + 2*papb + 2*pbp1 - pbp2) - 
     &             A2Lc*(-2*m2**2*pap1 + m3**2*pap2 - 2*m3**2*papb + p1p2*papb + m1**2*(p1p2 - 2*pbp1) - pap2*pbp1 + 
     &                pap1*pbp2)) + A1R*(A2Lc*m1*m3*(2*m2**2 - p1p2 - pap2 + 2*papb + 2*pbp1 - pbp2) - 
     &             A2Rc*(-2*m2**2*pap1 + m3**2*pap2 - 2*m3**2*papb + p1p2*papb + m1**2*(p1p2 - 2*pbp1) - pap2*pbp1 + 
     &                pap1*pbp2)))*(2*p1p2 + p1p3 - p2p3 - 2*pap2 - pap3 + 2*pbp2 + pbp3) + 
     &       2*(A1L*(A2Rc*m1*m3*(m1**2 - m2**2 - m3**2 + p1p2 + pap2 - 2*pbp1 + pbp2) + 
     &             A2Lc*(-(m2**2*pap1) - m3**2*pap1 + m3**2*pap2 - 2*m3**2*papb + p1p2*papb + 2*pap1*papb + 
     &                m1**2*(p1p2 + pap1 - 2*pbp1) - pap2*pbp1 + pap1*pbp2)) + 
     &          A1R*(A2Lc*m1*m3*(m1**2 - m2**2 - m3**2 + p1p2 + pap2 - 2*pbp1 + pbp2) + 
     &             A2Rc*(-(m2**2*pap1) - m3**2*pap1 + m3**2*pap2 - 2*m3**2*papb + p1p2*papb + 2*pap1*papb + 
     &                m1**2*(p1p2 + pap1 - 2*pbp1) - pap2*pbp1 + pap1*pbp2)))*(-p2p3 + pbp2 + 2*pbp3) + 
     &       (2*m2**2 - p1p2 + pap2 - 2*papb + 2*pbp1 - 3*pbp2)*
     &        (A1L*(A2Rc*m1*m3*(-p1p2 + p1p3 - pap2 + pap3 - pbp2 + pbp3) + 
     &             A2Lc*(m1**2*(-p1p2 + p1p3) + m3**2*(-pap2 + pap3) - p1p2*papb + p1p3*papb + pap2*pbp1 - pap3*pbp1 - 
     &                pap1*pbp2 + pap1*pbp3)) + A1R*
     &           (A2Lc*m1*m3*(-p1p2 + p1p3 - pap2 + pap3 - pbp2 + pbp3) + 
     &             A2Rc*(m1**2*(-p1p2 + p1p3) + m3**2*(-pap2 + pap3) - p1p2*papb + p1p3*papb + pap2*pbp1 - pap3*pbp1 - 
     &                pap1*pbp2 + pap1*pbp3)))))/
     &  (2.*p2p3*(m1**2 - m2**2 + m3**2 - 2*pap1)*(m1**2 + m2**2 - m3**2 + 2*papb)*pbp2))
     
        radgSgTf = (-CV**2*CF/2d0)*((GS**4*(A1L*(A2Rc*m1*m3*(p1p2**2 - 2*m3**2*p2p3 + 2*p1p3*p2p3 - m2**2*(4*p1p2 + 2*p1p3 + 3*p2p3) - 2*p2p3*pap1 + 
     &            p1p3*pap2 - 3*p2p3*pap3 - 2*p1p3*papb - 3*p2p3*papb - 4*p1p3*pbp1 - 3*p2p3*pbp1 + 2*pap2*pbp1 - 
     &            2*pap3*pbp1 + 2*m3**2*pbp2 + 2*pap1*pbp2 + 3*pap3*pbp2 + 4*pbp1*pbp2 + 
     &            p1p2*(2*p1p3 + p2p3 + pap2 + pap3 - 4*papb - 2*pbp1 + pbp2) + 4*m3**2*pbp3 - 4*p1p3*pbp3 - 3*p2p3*pbp3 + 
     &            4*pap1*pbp3 + 3*pap2*pbp3 + 6*pap3*pbp3 + 2*pbp1*pbp3 + 6*pbp2*pbp3 + 6*pbp3**2) + 
     &         A2Lc*(-2*m3**2*p2p3*pap1 + m3**2*p1p2*pap2 + m3**2*p1p3*pap2 + m3**2*p1p2*pap3 - 3*m3**2*p2p3*pap3 - 
     &            m2**2*(3*p2p3*pap1 + p1p2*(4*pap1 + pap3) + p1p3*(2*pap1 - 3*pap2 + 2*pap3)) - 4*m3**2*p1p2*papb + 
     &            p1p2**2*papb - 2*m3**2*p1p3*papb + 2*p1p2*p1p3*papb - 5*m3**2*p2p3*papb + p1p2*p2p3*papb + 
     &            2*p1p3*p2p3*papb + 2*m3**2*pap2*pbp1 - p1p2*pap2*pbp1 - p1p3*pap2*pbp1 - 2*m3**2*pap3*pbp1 - 
     &            p1p2*pap3*pbp1 + p2p3*pap3*pbp1 + 2*p1p2*papb*pbp1 - 2*p1p3*papb*pbp1 + 2*p2p3*papb*pbp1 - 2*pap2*pbp1**2 + 
     &            2*pap3*pbp1**2 + 2*m3**2*pap1*pbp2 + p1p2*pap1*pbp2 + p1p3*pap1*pbp2 + 3*m3**2*pap3*pbp2 + p1p2*pap3*pbp2 + 
     &            p1p3*pap3*pbp2 + 2*m3**2*papb*pbp2 - 4*p1p3*papb*pbp2 + 2*pap1*pbp1*pbp2 - 2*pap3*pbp1*pbp2 + 
     &            4*m3**2*pap1*pbp3 + p1p2*pap1*pbp3 - p2p3*pap1*pbp3 + 3*m3**2*pap2*pbp3 - p1p2*pap2*pbp3 - p1p3*pap2*pbp3 + 
     &            6*m3**2*pap3*pbp3 + 4*m3**2*papb*pbp3 - 2*p1p3*papb*pbp3 - 2*pap1*pbp1*pbp3 - 2*pap2*pbp1*pbp3 - 
     &            2*pap3*pbp1*pbp3 + 4*pap1*pbp2*pbp3 + 2*pap1*pbp3**2 + 
     &            m1**2*(p1p2**2 + 2*p1p3*p2p3 - 4*p1p3*pbp1 - p2p3*pbp1 - p1p3*pbp2 + p1p2*(2*p1p3 + p2p3 - 2*pbp1 - pbp3) - 
     &               4*p1p3*pbp3 + m3**2*(-2*p2p3 + 2*pbp2 + 4*pbp3)))) + 
     &      A1R*(A2Lc*m1*m3*(p1p2**2 - 2*m3**2*p2p3 + 2*p1p3*p2p3 - m2**2*(4*p1p2 + 2*p1p3 + 3*p2p3) - 2*p2p3*pap1 + 
     &            p1p3*pap2 - 3*p2p3*pap3 - 2*p1p3*papb - 3*p2p3*papb - 4*p1p3*pbp1 - 3*p2p3*pbp1 + 2*pap2*pbp1 - 
     &            2*pap3*pbp1 + 2*m3**2*pbp2 + 2*pap1*pbp2 + 3*pap3*pbp2 + 4*pbp1*pbp2 + 
     &            p1p2*(2*p1p3 + p2p3 + pap2 + pap3 - 4*papb - 2*pbp1 + pbp2) + 4*m3**2*pbp3 - 4*p1p3*pbp3 - 3*p2p3*pbp3 + 
     &            4*pap1*pbp3 + 3*pap2*pbp3 + 6*pap3*pbp3 + 2*pbp1*pbp3 + 6*pbp2*pbp3 + 6*pbp3**2) + 
     &         A2Rc*(-2*m3**2*p2p3*pap1 + m3**2*p1p2*pap2 + m3**2*p1p3*pap2 + m3**2*p1p2*pap3 - 3*m3**2*p2p3*pap3 - 
     &            m2**2*(3*p2p3*pap1 + p1p2*(4*pap1 + pap3) + p1p3*(2*pap1 - 3*pap2 + 2*pap3)) - 4*m3**2*p1p2*papb + 
     &            p1p2**2*papb - 2*m3**2*p1p3*papb + 2*p1p2*p1p3*papb - 5*m3**2*p2p3*papb + p1p2*p2p3*papb + 
     &            2*p1p3*p2p3*papb + 2*m3**2*pap2*pbp1 - p1p2*pap2*pbp1 - p1p3*pap2*pbp1 - 2*m3**2*pap3*pbp1 - 
     &            p1p2*pap3*pbp1 + p2p3*pap3*pbp1 + 2*p1p2*papb*pbp1 - 2*p1p3*papb*pbp1 + 2*p2p3*papb*pbp1 - 2*pap2*pbp1**2 + 
     &            2*pap3*pbp1**2 + 2*m3**2*pap1*pbp2 + p1p2*pap1*pbp2 + p1p3*pap1*pbp2 + 3*m3**2*pap3*pbp2 + p1p2*pap3*pbp2 + 
     &            p1p3*pap3*pbp2 + 2*m3**2*papb*pbp2 - 4*p1p3*papb*pbp2 + 2*pap1*pbp1*pbp2 - 2*pap3*pbp1*pbp2 + 
     &            4*m3**2*pap1*pbp3 + p1p2*pap1*pbp3 - p2p3*pap1*pbp3 + 3*m3**2*pap2*pbp3 - p1p2*pap2*pbp3 - p1p3*pap2*pbp3 + 
     &            6*m3**2*pap3*pbp3 + 4*m3**2*papb*pbp3 - 2*p1p3*papb*pbp3 - 2*pap1*pbp1*pbp3 - 2*pap2*pbp1*pbp3 - 
     &            2*pap3*pbp1*pbp3 + 4*pap1*pbp2*pbp3 + 2*pap1*pbp3**2 + 
     &            m1**2*(p1p2**2 + 2*p1p3*p2p3 - 4*p1p3*pbp1 - p2p3*pbp1 - p1p3*pbp2 + p1p2*(2*p1p3 + p2p3 - 2*pbp1 - pbp3) - 
     &               4*p1p3*pbp3 + m3**2*(-2*p2p3 + 2*pbp2 + 4*pbp3))))))/(2.*p1p3*p2p3*(m1**2 + m2**2 - m3**2 + 2*papb)*pbp2))
     
! colour class +1/2 CA**2 CF
        radgSiTg = (+CV**2*CF/2d0)*(-(GS**4*(A1L*(A2Rc*m1*m3*(-p2p3**2 + p2p3*pap2 - 4*m2**2*(p1p2 + 2*p1p3 + pap2 + 2*pap3) - 4*p2p3*papb - 
     &             4*p2p3*pbp1 + 2*p1p3*pbp2 + 3*p2p3*pbp2 - 2*pap2*pbp2 + 2*pap3*pbp2 + 4*papb*pbp2 + 4*pbp1*pbp2 - 
     &             2*pbp2**2 + 4*p1p3*pbp3 + 2*p2p3*pbp3 + 2*pap2*pbp3 + 4*pap3*pbp3 + 8*papb*pbp3 + 8*pbp1*pbp3 - 
     &             4*pbp3**2 + p1p2*(p2p3 - 2*pbp2 + 2*pbp3)) + 
     &          A2Lc*(-(p2p3**2*pap1) + p2p3*(-(p1p2*pap3) + m3**2*(pap2 - 4*papb) + p1p3*(pap2 - 4*papb) + p1p2*papb - 
     &                pap2*pbp1 + 4*pap3*pbp1 + 3*pap1*pbp2 + 2*pap1*pbp3) + 
     &             m1**2*(-4*m2**2*(p1p2 + 2*p1p3) - 4*p2p3*pbp1 + 2*p1p3*pbp2 + 4*pbp1*pbp2 + 4*p1p3*pbp3 + 8*pbp1*pbp3 + 
     &                p1p2*(p2p3 - 2*pbp2 + 2*pbp3)) - 
     &             2*(2*m2**2*(p1p3*pap2 - p1p2*pap3 + m3**2*(pap2 + 2*pap3) + p1p2*papb + 2*p1p3*papb - pap2*pbp1 - 
     &                   2*pap3*pbp1) + p1p3*pap2*pbp2 - p1p2*pap3*pbp2 + p1p2*papb*pbp2 - 3*p1p3*papb*pbp2 - 
     &                pap2*pbp1*pbp2 + 3*pap3*pbp1*pbp2 + pap1*pbp2**2 - p1p3*pap2*pbp3 + p1p2*pap3*pbp3 - p1p2*papb*pbp3 - 
     &                6*p1p3*papb*pbp3 + pap2*pbp1*pbp3 + 6*pap3*pbp1*pbp3 + 2*pap1*pbp3**2 + 
     &                m3**2*(pap2*(pbp2 - pbp3) - (pap3 + 2*papb)*(pbp2 + 2*pbp3))))) + 
     &       A1R*(A2Lc*m1*m3*(-p2p3**2 + p2p3*pap2 - 4*m2**2*(p1p2 + 2*p1p3 + pap2 + 2*pap3) - 4*p2p3*papb - 4*p2p3*pbp1 + 
     &             2*p1p3*pbp2 + 3*p2p3*pbp2 - 2*pap2*pbp2 + 2*pap3*pbp2 + 4*papb*pbp2 + 4*pbp1*pbp2 - 2*pbp2**2 + 
     &             4*p1p3*pbp3 + 2*p2p3*pbp3 + 2*pap2*pbp3 + 4*pap3*pbp3 + 8*papb*pbp3 + 8*pbp1*pbp3 - 4*pbp3**2 + 
     &             p1p2*(p2p3 - 2*pbp2 + 2*pbp3)) + 
     &          A2Rc*(-(p2p3**2*pap1) + p2p3*(-(p1p2*pap3) + m3**2*(pap2 - 4*papb) + p1p3*(pap2 - 4*papb) + p1p2*papb - 
     &                pap2*pbp1 + 4*pap3*pbp1 + 3*pap1*pbp2 + 2*pap1*pbp3) + 
     &             m1**2*(-4*m2**2*(p1p2 + 2*p1p3) - 4*p2p3*pbp1 + 2*p1p3*pbp2 + 4*pbp1*pbp2 + 4*p1p3*pbp3 + 8*pbp1*pbp3 + 
     &                p1p2*(p2p3 - 2*pbp2 + 2*pbp3)) - 
     &             2*(2*m2**2*(p1p3*pap2 - p1p2*pap3 + m3**2*(pap2 + 2*pap3) + p1p2*papb + 2*p1p3*papb - pap2*pbp1 - 
     &                   2*pap3*pbp1) + p1p3*pap2*pbp2 - p1p2*pap3*pbp2 + p1p2*papb*pbp2 - 3*p1p3*papb*pbp2 - 
     &                pap2*pbp1*pbp2 + 3*pap3*pbp1*pbp2 + pap1*pbp2**2 - p1p3*pap2*pbp3 + p1p2*pap3*pbp3 - p1p2*papb*pbp3 - 
     &                6*p1p3*papb*pbp3 + pap2*pbp1*pbp3 + 6*pap3*pbp1*pbp3 + 2*pap1*pbp3**2 + 
     &                m3**2*(pap2*(pbp2 - pbp3) - (pap3 + 2*papb)*(pbp2 + 2*pbp3)))))))/
     &  (4.*p1p2*p2p3*(m1**2 - m2**2 + m3**2 - 2*pap1)*pbp3))
     
        radgSpTg = (+CV**2*CF/2d0)*(-((GS**4*(A2Lc*(A1R*m1*m3*(p1p2**2 - p1p3**2 - p1p3*p2p3 + m2**2*(2*p1p2 + 4*p1p3 + p2p3) - p1p3*pap3 - p2p3*pap3 + 
     &              4*p1p3*papb + p2p3*papb + 2*p1p3*pbp1 + p2p3*pbp1 - 4*pap2*pbp1 - 2*pap3*pbp1 + 2*m3**2*pbp2 + 
     &              p1p3*pbp2 + 2*pap1*pbp2 - 2*pap2*pbp2 + pap3*pbp2 - 2*pbp1*pbp2 - 2*pbp2**2 - 2*m3**2*pbp3 - p1p3*pbp3 - 
     &              p2p3*pbp3 - 2*pap1*pbp3 + pap2*pbp3 - 4*pbp1*pbp3 + 2*pbp2*pbp3 + 
     &              p1p2*(pap2 + 2*papb - 2*pbp1 - pbp2 + pbp3)) + 
     &           A1L*(m3**2*p1p2*pap2 - m3**2*p1p3*pap3 - m3**2*p2p3*pap3 + 
     &              m2**2*(2*p1p2*pap1 + 4*p1p3*pap1 + p2p3*pap1 + 3*p1p3*pap2 - 3*p1p2*pap3) + 2*m3**2*p1p2*papb + 
     &              p1p2**2*papb + 4*m3**2*p1p3*papb - p1p3**2*papb + m3**2*p2p3*papb - p1p3*p2p3*papb - 4*m3**2*pap2*pbp1 - 
     &              p1p2*pap2*pbp1 - 2*m3**2*pap3*pbp1 + p1p3*pap3*pbp1 + p2p3*pap3*pbp1 - 4*p1p2*papb*pbp1 - 
     &              2*p1p3*papb*pbp1 + 4*pap2*pbp1**2 + 2*pap3*pbp1**2 + 2*m3**2*pap1*pbp2 + p1p2*pap1*pbp2 - 
     &              2*m3**2*pap2*pbp2 + m3**2*pap3*pbp2 + 2*m3**2*papb*pbp2 - 2*p1p2*papb*pbp2 - 2*p1p3*papb*pbp2 - 
     &              4*pap1*pbp1*pbp2 + 2*pap2*pbp1*pbp2 + 2*pap3*pbp1*pbp2 - 2*pap1*pbp2**2 - 2*m3**2*pap1*pbp3 - 
     &              p1p3*pap1*pbp3 - p2p3*pap1*pbp3 + m3**2*pap2*pbp3 - 2*m3**2*papb*pbp3 + 4*p1p2*papb*pbp3 - 
     &              2*pap1*pbp1*pbp3 - 4*pap2*pbp1*pbp3 + 2*pap1*pbp2*pbp3 + 
     &              m1**2*(p1p2**2 - p1p3**2 + p2p3*pbp1 + 2*m3**2*pbp2 + p1p3*(-p2p3 + 2*pbp1 + pbp2) - 2*m3**2*pbp3 + 
     &                 p1p2*(-2*pbp1 - 2*pbp2 + pbp3)))) + 
     &        A2Rc*(A1L*m1*m3*(p1p2**2 - p1p3**2 - p1p3*p2p3 + m2**2*(2*p1p2 + 4*p1p3 + p2p3) - p1p3*pap3 - p2p3*pap3 + 
     &              4*p1p3*papb + p2p3*papb + 2*p1p3*pbp1 + p2p3*pbp1 - 4*pap2*pbp1 - 2*pap3*pbp1 + 2*m3**2*pbp2 + 
     &              p1p3*pbp2 + 2*pap1*pbp2 - 2*pap2*pbp2 + pap3*pbp2 - 2*pbp1*pbp2 - 2*pbp2**2 - 2*m3**2*pbp3 - p1p3*pbp3 - 
     &              p2p3*pbp3 - 2*pap1*pbp3 + pap2*pbp3 - 4*pbp1*pbp3 + 2*pbp2*pbp3 + 
     &              p1p2*(pap2 + 2*papb - 2*pbp1 - pbp2 + pbp3)) + 
     &           A1R*(m3**2*p1p2*pap2 - m3**2*p1p3*pap3 - m3**2*p2p3*pap3 + 
     &              m2**2*(2*p1p2*pap1 + 4*p1p3*pap1 + p2p3*pap1 + 3*p1p3*pap2 - 3*p1p2*pap3) + 2*m3**2*p1p2*papb + 
     &              p1p2**2*papb + 4*m3**2*p1p3*papb - p1p3**2*papb + m3**2*p2p3*papb - p1p3*p2p3*papb - 4*m3**2*pap2*pbp1 - 
     &              p1p2*pap2*pbp1 - 2*m3**2*pap3*pbp1 + p1p3*pap3*pbp1 + p2p3*pap3*pbp1 - 4*p1p2*papb*pbp1 - 
     &              2*p1p3*papb*pbp1 + 4*pap2*pbp1**2 + 2*pap3*pbp1**2 + 2*m3**2*pap1*pbp2 + p1p2*pap1*pbp2 - 
     &              2*m3**2*pap2*pbp2 + m3**2*pap3*pbp2 + 2*m3**2*papb*pbp2 - 2*p1p2*papb*pbp2 - 2*p1p3*papb*pbp2 - 
     &              4*pap1*pbp1*pbp2 + 2*pap2*pbp1*pbp2 + 2*pap3*pbp1*pbp2 - 2*pap1*pbp2**2 - 2*m3**2*pap1*pbp3 - 
     &              p1p3*pap1*pbp3 - p2p3*pap1*pbp3 + m3**2*pap2*pbp3 - 2*m3**2*papb*pbp3 + 4*p1p2*papb*pbp3 - 
     &              2*pap1*pbp1*pbp3 - 4*pap2*pbp1*pbp3 + 2*pap1*pbp2*pbp3 + 
     &              m1**2*(p1p2**2 - p1p3**2 + p2p3*pbp1 + 2*m3**2*pbp2 + p1p3*(-p2p3 + 2*pbp1 + pbp2) - 2*m3**2*pbp3 + 
     &                 p1p2*(-2*pbp1 - 2*pbp2 + pbp3))))))/
     &    (p1p2*p2p3*(m1**2 - m2**2 + m3**2 - 2*pap1)*(m1**2 + m2**2 - m3**2 + 2*papb))))
     
        radgSgTi = (+CV**2*CF/2d0)*( -(GS**4*((p1p2 + 2*p1p3 - p2p3 - pap2 - 2*pap3 + pbp2 + 2*pbp3)*
     &        (A1L*(A2Rc*m1*m3*(-2*m2**2 + p1p3 + pap3 - 2*papb - 2*pbp1 + pbp3) + 
     &             A2Lc*(-2*m2**2*pap1 + m3**2*pap3 - 2*m3**2*papb + p1p3*papb + m1**2*(p1p3 - 2*pbp1) - pap3*pbp1 + 
     &                pap1*pbp3)) + A1R*(A2Lc*m1*m3*(-2*m2**2 + p1p3 + pap3 - 2*papb - 2*pbp1 + pbp3) + 
     &             A2Rc*(-2*m2**2*pap1 + m3**2*pap3 - 2*m3**2*papb + p1p3*papb + m1**2*(p1p3 - 2*pbp1) - pap3*pbp1 + 
     &                pap1*pbp3))) + 2*(p2p3 - 2*pbp2 - pbp3)*
     &        (A1L*(A2Rc*m1*m3*(m1**2 - m2**2 - m3**2 + p1p3 + pap3 - 2*pbp1 + pbp3) + 
     &             A2Lc*(-(m2**2*pap1) - m3**2*pap1 + m3**2*pap3 - 2*m3**2*papb + p1p3*papb + 2*pap1*papb + 
     &                m1**2*(p1p3 + pap1 - 2*pbp1) - pap3*pbp1 + pap1*pbp3)) + 
     &          A1R*(A2Lc*m1*m3*(m1**2 - m2**2 - m3**2 + p1p3 + pap3 - 2*pbp1 + pbp3) + 
     &             A2Rc*(-(m2**2*pap1) - m3**2*pap1 + m3**2*pap3 - 2*m3**2*papb + p1p3*papb + 2*pap1*papb + 
     &                m1**2*(p1p3 + pap1 - 2*pbp1) - pap3*pbp1 + pap1*pbp3))) + 
     &       (2*m2**2 - p1p3 + pap3 - 2*papb + 2*pbp1 - 3*pbp3)*
     &        (A1L*(A2Rc*m1*m3*(-p1p2 + p1p3 - pap2 + pap3 - pbp2 + pbp3) + 
     &             A2Lc*(m1**2*(-p1p2 + p1p3) + m3**2*(-pap2 + pap3) - p1p2*papb + p1p3*papb + pap2*pbp1 - pap3*pbp1 - 
     &                pap1*pbp2 + pap1*pbp3)) + A1R*
     &           (A2Lc*m1*m3*(-p1p2 + p1p3 - pap2 + pap3 - pbp2 + pbp3) + 
     &             A2Rc*(m1**2*(-p1p2 + p1p3) + m3**2*(-pap2 + pap3) - p1p2*papb + p1p3*papb + pap2*pbp1 - pap3*pbp1 - 
     &                pap1*pbp2 + pap1*pbp3)))))/
     &  (2.*p2p3*(m1**2 - m2**2 + m3**2 - 2*pap1)*(m1**2 + m2**2 - m3**2 + 2*papb)*pbp3))

                  
        NeuQ2qx_M2RadiationGSTsoftcoll = 2d0*(0.5*(radgSfTg + radgSgTp + radgSgTf + radgSiTg + radgSpTg + radgSgTi))


!         print*,"NeuQ2qx_M2RadiationGSTsoftcoll",NeuQ2qx_M2RadiationGSTsoftcoll
        

    
        end


        double precision function NeuQ2qx_M2RadiationGSTcoll()
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double complex radgSgTg, radgSnTn1, radgSnTn2
        double complex A1L, A1R, A2Lc, A2Rc
        double precision wquark


        NeuQ2qx_M2RadiationGSTcoll = 0d0
     

        radgSnTn1 = 0d0
        radgSnTn2 = 0d0

        radgSgTg = 0d0


        A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A1R = AR(chimeigen,sfmeigen,sftype,sfgen)


        wquark = 0d0

        A2Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A2Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))

     
! colour class CA**2 CF
        radgSgTg = (CV**2*CF)*((GS**4*(A2Lc*(A1R*m1*m3*(16*m2**2*p2p3 + p1p3*p2p3 + p2p3*pap2 + p2p3*pap3 + 16*p2p3*papb + 16*p2p3*pbp1 - 
     &            14*p1p3*pbp2 + p2p3*pbp2 - 4*pap2*pbp2 - 14*pap3*pbp2 - 4*pbp2**2 + p1p2*(p2p3 - 4*pbp2 - 14*pbp3) - 
     &            4*p1p3*pbp3 + p2p3*pbp3 - 14*pap2*pbp3 - 4*pap3*pbp3 - 28*pbp2*pbp3 - 4*pbp3**2) + 
     &         A1L*(16*m2**2*p2p3*pap1 + m3**2*p2p3*pap2 + m3**2*p2p3*pap3 + 16*m3**2*p2p3*papb + p1p2*p2p3*papb + 
     &            p1p3*p2p3*papb - p2p3*pap2*pbp1 - p2p3*pap3*pbp1 + p2p3*pap1*pbp2 - 4*m3**2*pap2*pbp2 - 
     &            14*m3**2*pap3*pbp2 - 4*p1p2*papb*pbp2 - 14*p1p3*papb*pbp2 + 4*pap2*pbp1*pbp2 + 14*pap3*pbp1*pbp2 - 
     &            4*pap1*pbp2**2 + m1**2*(16*p2p3*pbp1 + p1p2*(p2p3 - 4*pbp2 - 14*pbp3) + p1p3*(p2p3 - 14*pbp2 - 4*pbp3)) + 
     &            p2p3*pap1*pbp3 - 14*m3**2*pap2*pbp3 - 4*m3**2*pap3*pbp3 - 14*p1p2*papb*pbp3 - 4*p1p3*papb*pbp3 + 
     &            14*pap2*pbp1*pbp3 + 4*pap3*pbp1*pbp3 - 28*pap1*pbp2*pbp3 - 4*pap1*pbp3**2)) + 
     &      A2Rc*(A1L*m1*m3*(16*m2**2*p2p3 + p1p3*p2p3 + p2p3*pap2 + p2p3*pap3 + 16*p2p3*papb + 16*p2p3*pbp1 - 14*p1p3*pbp2 + 
     &            p2p3*pbp2 - 4*pap2*pbp2 - 14*pap3*pbp2 - 4*pbp2**2 + p1p2*(p2p3 - 4*pbp2 - 14*pbp3) - 4*p1p3*pbp3 + 
     &            p2p3*pbp3 - 14*pap2*pbp3 - 4*pap3*pbp3 - 28*pbp2*pbp3 - 4*pbp3**2) + 
     &         A1R*(16*m2**2*p2p3*pap1 + m3**2*p2p3*pap2 + m3**2*p2p3*pap3 + 16*m3**2*p2p3*papb + p1p2*p2p3*papb + 
     &            p1p3*p2p3*papb - p2p3*pap2*pbp1 - p2p3*pap3*pbp1 + p2p3*pap1*pbp2 - 4*m3**2*pap2*pbp2 - 
     &            14*m3**2*pap3*pbp2 - 4*p1p2*papb*pbp2 - 14*p1p3*papb*pbp2 + 4*pap2*pbp1*pbp2 + 14*pap3*pbp1*pbp2 - 
     &            4*pap1*pbp2**2 + m1**2*(16*p2p3*pbp1 + p1p2*(p2p3 - 4*pbp2 - 14*pbp3) + p1p3*(p2p3 - 14*pbp2 - 4*pbp3)) + 
     &            p2p3*pap1*pbp3 - 14*m3**2*pap2*pbp3 - 4*m3**2*pap3*pbp3 - 14*p1p2*papb*pbp3 - 4*p1p3*papb*pbp3 + 
     &            14*pap2*pbp1*pbp3 + 4*pap3*pbp1*pbp3 - 28*pap1*pbp2*pbp3 - 4*pap1*pbp3**2))))/
     &  (2.*p2p3**2*(m1**2 - m2**2 + m3**2 - 2*pap1)*(m1**2 + m2**2 - m3**2 + 2*papb)))

        
! ghosts      
         radgSnTn1 = (CV**2*CF)*(GS**4*(p1p2 - pap2 + pbp2)*(A1L*(A2Rc*m1*m3*(p1p2 + pap2 + pbp2) + A2Lc*(m1**2*p1p2 + m3**2*pap2 + p1p2*papb - pap2*pbp1 + pap1*pbp2)) + 
     &      A1R*(A2Lc*m1*m3*(p1p2 + pap2 + pbp2) + A2Rc*(m1**2*p1p2 + m3**2*pap2 + p1p2*papb - pap2*pbp1 + pap1*pbp2))))/
     &  (2.*p2p3**2*(m1**2 - m2**2 + m3**2 - 2*pap1)*(m1**2 + m2**2 - m3**2 + 2*papb))
         
         radgSnTn2 = (CV**2*CF)*(GS**4*(p1p3 - pap3 + pbp3)*(A1L*(A2Rc*m1*m3*(p1p3 + pap3 + pbp3) + A2Lc*(m1**2*p1p3 + m3**2*pap3 + p1p3*papb - pap3*pbp1 + pap1*pbp3)) + 
     &      A1R*(A2Lc*m1*m3*(p1p3 + pap3 + pbp3) + A2Rc*(m1**2*p1p3 + m3**2*pap3 + p1p3*papb - pap3*pbp1 + pap1*pbp3))))/
     &  (2.*p2p3**2*(m1**2 - m2**2 + m3**2 - 2*pap1)*(m1**2 + m2**2 - m3**2 + 2*papb))
   
                  
        NeuQ2qx_M2RadiationGSTcoll = 2d0*(0.5*(radgSgTg) - 0.5d0*(radgSnTn1 + radgSnTn2))


!         print*,"NeuQ2qx_M2RadiationGSTcoll",NeuQ2qx_M2RadiationGSTcoll
        

    
        end
        
! =================================================================
        double precision function NeuQ2qx_M2RadiationGSTcon()
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double complex radgSiTv, radgSpTv, radgSfTv
        double complex A1L, A1R, A2Lc, A2Rc
        double precision wquark


        NeuQ2qx_M2RadiationGSTcon = 0d0

        radgSiTv = 0d0
        radgSpTv = 0d0
        radgSfTv = 0d0
        

        A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A1R = AR(chimeigen,sfmeigen,sftype,sfgen)


        wquark = 0d0

        A2Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A2Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))        
        
! colour class (CF**2*CA - 1/2*CF)
        radgSiTv = (CF**2*CV - CF/2d0)*(-(GS**4*(A1L*(A2Rc*m1*m3*(-2*m2**2 + p1p3 + pap3 - 2*papb - 2*pbp1 + 3*pbp3) + 
     &          A2Lc*(-2*m2**2*pap1 + m3**2*pap3 - 2*m3**2*papb - p1p3*papb + m1**2*(p1p3 - 2*pbp1) + pap3*pbp1 + 3*pap1*pbp3)
     &          ) + A1R*(A2Lc*m1*m3*(-2*m2**2 + p1p3 + pap3 - 2*papb - 2*pbp1 + 3*pbp3) + 
     &          A2Rc*(-2*m2**2*pap1 + m3**2*pap3 - 2*m3**2*papb - p1p3*papb + m1**2*(p1p3 - 2*pbp1) + pap3*pbp1 + 3*pap1*pbp3)
     &          )))/(2.*p1p2*(m1**2 - m2**2 + m3**2 - 2*pap1)*pbp3))
     
        radgSpTv = (CF**2*CV - CF/2d0)*((2*GS**4*(A1L*(A2Rc*m1*m3*(-m3**2 + p1p2 - pap1 + pap2 - pbp1 + pbp2) + 
     &         A2Lc*(m1**2*(-m3**2 + p1p2) + p1p2*papb - m3**2*(pap1 - pap2 + papb) - pap2*pbp1 + pap1*pbp2)) + 
     &      A1R*(A2Lc*m1*m3*(-m3**2 + p1p2 - pap1 + pap2 - pbp1 + pbp2) + 
     &         A2Rc*(m1**2*(-m3**2 + p1p2) + p1p2*papb - m3**2*(pap1 - pap2 + papb) - pap2*pbp1 + pap1*pbp2))))/
     &  (p1p2*(m1**2 - m2**2 + m3**2 - 2*pap1)*(m1**2 + m2**2 - m3**2 + 2*papb)))
     
        radgSfTv = (CF**2*CV - CF/2d0)*((2*GS**4*(A1L*(A2Rc*m1*m3*(-m3**2 + p1p3 - pap1 + pap3 - pbp1 + pbp3) + 
     &         A2Lc*(m1**2*(-m3**2 + p1p3) + p1p3*papb - m3**2*(pap1 - pap3 + papb) - pap3*pbp1 + pap1*pbp3)) + 
     &      A1R*(A2Lc*m1*m3*(-m3**2 + p1p3 - pap1 + pap3 - pbp1 + pbp3) + 
     &         A2Rc*(m1**2*(-m3**2 + p1p3) + p1p3*papb - m3**2*(pap1 - pap3 + papb) - pap3*pbp1 + pap1*pbp3))))/
     &  (p1p3*(m1**2 - m2**2 + m3**2 - 2*pap1)*(m1**2 + m2**2 - m3**2 + 2*papb)))
        
        
        NeuQ2qx_M2RadiationGSTcon = 2d0*(0.5*(radgSiTv + radgSpTv + radgSfTv))

        
!          print*,"radgSiTv",radgSiTv
!          print*,"radgSpTv",radgSpTv
!          print*,"radgSfTv",radgSfTv
        
!         print*,"NeuQ2qx_M2RadiationGSTcon",NeuQ2qx_M2RadiationGSTcon


    
        end
        
        

! =================================================================
        double precision function NeuQ2qx_M2RadiationGTTsoft()
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double complex radgTiTi, radgTpTp, radgTfTf
        double complex radgTfTp, radgTfTi, radgTpTi

        double complex A2L, A2Lc, A2R, A2Rc
     
        NeuQ2qx_M2RadiationGTTsoft = 0d0


        radgTiTi = 0d0
        radgTpTp = 0d0
        radgTfTf = 0d0
        
        radgTfTi = 0d0
        radgTpTi = 0d0
        radgTfTp = 0d0

        A2L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A2R = AR(chimeigen,sfmeigen,sftype,sfgen)
        A2Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A2Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
! colour class CF**2 CA
        radgTiTi = (CF**2*CV)*((2*GS**4*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*
     &    (m1**2 + m2**2 + m3**2 - 2*p1p3 - 2*pap1 + 2*pap3 - 2*papb + 2*pbp1 - 2*pbp3)*(m2**2 - pbp3))/
     &  ((m1**2 - m2**2 + m3**2 - 2*pap1)**2*pbp3**2))
     
        radgTpTp = (CF**2*CV)*((2*GS**4*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*
     &    (m1**2 + m2**2 + m3**2 - 2*p1p2 - 2*pap1 + 2*pap2 - 2*papb + 2*pbp1 - 2*pbp2)*(m2**2 - pbp2))/
     &  ((m1**2 - m2**2 + m3**2 - 2*pap1)**2*pbp2**2))
     
        radgTfTf = (CF**2*CV)*((2*GS**4*(A2Lc*(A2R*m1*m3*(m3**2 + p1p3) + A2L*(-(p1p3*pap3) + m3**2*(pap1 + pap3))) + 
     &      A2Rc*(A2L*m1*m3*(m3**2 + p1p3) + A2R*(-(p1p3*pap3) + m3**2*(pap1 + pap3))))*(m2**2 - pbp2))/(p1p3**2*pbp2**2))
     
        radgTfTp = (CF**2*CV)*((GS**4*(m2**2 - pbp2)*(A2Lc*(A2R*m1*m3*(2*m3**2 - 2*p1p2 + p1p3 - p2p3 - 2*pap1 - pap3 + 2*pbp1 + pbp3) + 
     &         A2L*(m1**2*p1p3 - 2*p1p2*pap1 - p2p3*pap1 - 2*pap1**2 + p1p3*pap2 - p1p2*pap3 - 2*pap1*pap3 + 
     &            m3**2*(2*pap1 + pap3) - p1p3*papb + 2*pap1*pbp1 + pap3*pbp1 + pap1*pbp3)) + 
     &      A2Rc*(A2L*m1*m3*(2*m3**2 - 2*p1p2 + p1p3 - p2p3 - 2*pap1 - pap3 + 2*pbp1 + pbp3) + 
     &         A2R*(m1**2*p1p3 - 2*p1p2*pap1 - p2p3*pap1 - 2*pap1**2 + p1p3*pap2 - p1p2*pap3 - 2*pap1*pap3 + 
     &            m3**2*(2*pap1 + pap3) - p1p3*papb + 2*pap1*pbp1 + pap3*pbp1 + pap1*pbp3))))/
     &  (p1p3*(m1**2 - m2**2 + m3**2 - 2*pap1)*pbp2**2))
     
! colour class -CF/2
        radgTfTi = (-CF/2d0)*( -(GS**4*(2*m2**2 - p1p2 + p2p3 + pap2 - 2*papb + 2*pbp1 - pbp2 - 2*pbp3)*
     &     (A2Lc*(A2R*m1*m3*(p1p3 - 2*pbp1 - pbp3) + A2L*(p1p3*(pap1 + papb) - pap3*pbp1 - pap1*(2*pbp1 + pbp3))) + 
     &       A2Rc*(A2L*m1*m3*(p1p3 - 2*pbp1 - pbp3) + A2R*(p1p3*(pap1 + papb) - pap3*pbp1 - pap1*(2*pbp1 + pbp3)))))/
     &  (2.*p1p3*(m1**2 - m2**2 + m3**2 - 2*pap1)*pbp2*pbp3))
     
        radgTpTi = (-CF/2d0)*((GS**4*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*
     &    (2*m2**2 - p1p2 + p2p3 + pap2 - 2*papb + 2*pbp1 - pbp2 - 2*pbp3)*
     &    (2*m2**2 - p1p3 + p2p3 + pap3 - 2*papb + 2*pbp1 - 2*pbp2 - pbp3))/(2.*(m1**2 - m2**2 + m3**2 - 2*pap1)**2*pbp2*pbp3))
    
     
        NeuQ2qx_M2RadiationGTTsoft = 0.5*(radgTiTi + radgTpTp + radgTfTf + 2d0*(radgTfTi + radgTpTi + radgTfTp))
        
!         print*,"NeuQ2qx_M2RadiationGTTsoft",NeuQ2qx_M2RadiationGTTsoft


    
        end
        
! =================================================================
        double precision function NeuQ2qx_M2RadiationGTTsoftcoll()
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double complex radgTpTg, radgTfTg, radgTiTg

        double complex A2L, A2Lc, A2R, A2Rc
     
        NeuQ2qx_M2RadiationGTTsoftcoll = 0d0

        radgTpTg = 0d0
        radgTfTg = 0d0
        radgTiTg = 0d0

        A2L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A2R = AR(chimeigen,sfmeigen,sftype,sfgen)
        A2Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A2Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
     
! colour class -1/2 CA**2 CF
        radgTpTg = (-CV**2*CF/2d0)*( -(GS**4*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*
     &     (p2p3**2 + 4*m2**2*(2*p1p2 + p1p3 - 2*pap2 - pap3) + p2p3*pap3 - 4*p2p3*papb + 4*p2p3*pbp1 - 4*p1p2*pbp2 - 
     &       2*p2p3*pbp2 + 4*pap2*pbp2 + 2*pap3*pbp2 + 8*papb*pbp2 - 8*pbp1*pbp2 + 4*pbp2**2 - 
     &       p1p3*(p2p3 + 2*pbp2 - 2*pbp3) - 2*p1p2*pbp3 - 3*p2p3*pbp3 + 2*pap2*pbp3 - 2*pap3*pbp3 + 4*papb*pbp3 - 
     &       4*pbp1*pbp3 + 2*pbp3**2))/(2.*p2p3*(m1**2 - m2**2 + m3**2 - 2*pap1)**2*pbp2))
     
        radgTfTg = (-CV**2*CF/2d0)*((GS**4*(A2Lc*(A2R*m1*m3*(-4*m2**2*(2*p1p2 + p1p3 + p2p3) - 4*p2p3*pbp1 + 4*p1p2*pbp2 + 2*p2p3*pbp2 + 8*pbp1*pbp2 + 
     &            p1p3*(p2p3 + 2*pbp2 - 2*pbp3) + 2*p1p2*pbp3 - p2p3*pbp3 + 4*pbp1*pbp3 + 4*pbp2*pbp3 + 2*pbp3**2) + 
     &         A2L*(-4*m2**2*(p2p3*pap1 + p1p3*(pap1 - pap2) + p1p2*(2*pap1 + pap3)) - 4*p2p3*pap1*pbp1 - 2*p2p3*pap3*pbp1 + 
     &            4*p1p2*pap1*pbp2 + 2*p2p3*pap1*pbp2 + 2*p1p2*pap3*pbp2 + 8*pap1*pbp1*pbp2 + 4*pap3*pbp1*pbp2 + 
     &            2*p1p2*pap1*pbp3 - p2p3*pap1*pbp3 + p1p2*pap3*pbp3 + 4*pap1*pbp1*pbp3 + 2*pap3*pbp1*pbp3 + 
     &            4*pap1*pbp2*pbp3 + 2*pap1*pbp3**2 + 
     &            p1p3*(p2p3*(pap1 + 2*papb) + 2*pap1*(pbp2 - pbp3) - (pap2 + 2*papb)*(2*pbp2 + pbp3)))) + 
     &      A2Rc*(A2L*m1*m3*(-4*m2**2*(2*p1p2 + p1p3 + p2p3) - 4*p2p3*pbp1 + 4*p1p2*pbp2 + 2*p2p3*pbp2 + 8*pbp1*pbp2 + 
     &            p1p3*(p2p3 + 2*pbp2 - 2*pbp3) + 2*p1p2*pbp3 - p2p3*pbp3 + 4*pbp1*pbp3 + 4*pbp2*pbp3 + 2*pbp3**2) + 
     &         A2R*(-4*m2**2*(p2p3*pap1 + p1p3*(pap1 - pap2) + p1p2*(2*pap1 + pap3)) - 4*p2p3*pap1*pbp1 - 2*p2p3*pap3*pbp1 + 
     &            4*p1p2*pap1*pbp2 + 2*p2p3*pap1*pbp2 + 2*p1p2*pap3*pbp2 + 8*pap1*pbp1*pbp2 + 4*pap3*pbp1*pbp2 + 
     &            2*p1p2*pap1*pbp3 - p2p3*pap1*pbp3 + p1p2*pap3*pbp3 + 4*pap1*pbp1*pbp3 + 2*pap3*pbp1*pbp3 + 
     &            4*pap1*pbp2*pbp3 + 2*pap1*pbp3**2 + 
     &            p1p3*(p2p3*(pap1 + 2*papb) + 2*pap1*(pbp2 - pbp3) - (pap2 + 2*papb)*(2*pbp2 + pbp3))))))/
     &  (2.*p1p3*p2p3*(m1**2 - m2**2 + m3**2 - 2*pap1)*pbp2))
     
! colour class +1/2 CA**2 CF
        radgTiTg = (+CV**2*CF/2d0)*((GS**4*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*
     &    (p2p3**2 + p2p3*pap2 + 4*m2**2*(p1p2 + 2*p1p3 - pap2 - 2*pap3) - 4*p2p3*papb + 4*p2p3*pbp1 - 2*p1p3*pbp2 - 
     &      3*p2p3*pbp2 - 2*pap2*pbp2 + 2*pap3*pbp2 + 4*papb*pbp2 - 4*pbp1*pbp2 + 2*pbp2**2 - 4*p1p3*pbp3 - 2*p2p3*pbp3 + 
     &      2*pap2*pbp3 + 4*pap3*pbp3 + 8*papb*pbp3 - 8*pbp1*pbp3 + 4*pbp3**2 - p1p2*(p2p3 - 2*pbp2 + 2*pbp3)))/
     &  (2.*p2p3*(m1**2 - m2**2 + m3**2 - 2*pap1)**2*pbp3))
  
          
        NeuQ2qx_M2RadiationGTTsoftcoll = 0.5*(2d0*(radgTpTg + radgTfTg + radgTiTg))


!          print*,"fsoft",fsoft
!          print*,"spflag",spflag
!          print*,"radgTiTi",radgTiTi
!          print*,"radgTfTf",radgTfTf
!          print*,"radgTfTi",radgTfTi
        
!         print*,"NeuQ2qx_M2RadiationGTTsoftcoll",NeuQ2qx_M2RadiationGTTsoftcoll


        end

! =================================================================
        double precision function NeuQ2qx_M2RadiationGTTcoll()
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double complex radgTgTg, radgTnTn1, radgTnTn2

        double complex A2L, A2Lc, A2R, A2Rc
     
        NeuQ2qx_M2RadiationGTTcoll = 0d0

        radgTgTg = 0d0
        
        radgTnTn1 = 0d0
        radgTnTn2 = 0d0


        A2L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A2R = AR(chimeigen,sfmeigen,sftype,sfgen)
        A2Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A2Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
     
! colour class CA**2 CF
        radgTgTg = (CV**2*CF)*(-((GS**4*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*
     &      (-16*m2**2*p2p3 + p2p3**2 - 2*p2p3*(pbp2 + pbp3) + 4*(pbp2**2 + 7*pbp2*pbp3 + pbp3**2)))/
     &    (p2p3**2*(m1**2 - m2**2 + m3**2 - 2*pap1)**2)))
        
! ghosts 
        radgTnTn1 =  (CV**2*CF)*(GS**4*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*(p1p2 - pap2 + pbp2)**2)/(2.*p2p3**2*(m1**2 - m2**2 + m3**2 - 2*pap1)**2)
        
        radgTnTn2 =  (CV**2*CF)*(GS**4*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*(p1p3 - pap3 + pbp3)**2)/(2.*p2p3**2*(m1**2 - m2**2 + m3**2 - 2*pap1)**2)

     
        NeuQ2qx_M2RadiationGTTcoll = 0.5*(radgTgTg)  - 0.5d0*(radgTnTn1 + radgTnTn2)

!          print*,"fsoft",fsoft
!          print*,"spflag",spflag
!          print*,"radgTiTi",radgTiTi
!          print*,"radgTfTf",radgTfTf
!          print*,"radgTfTi",radgTfTi
        
!         print*,"NeuQ2qx_M2RadiationGTTcoll",NeuQ2qx_M2RadiationGTTcoll


    
        end
                
        
! =================================================================
        double precision function NeuQ2qx_M2RadiationGTTcon()
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double complex radgTvTv, radgTiTv, radgTpTv, radgTfTv

        double complex A2L, A2Lc, A2R, A2Rc
     
        NeuQ2qx_M2RadiationGTTcon = 0d0

        radgTvTv = 0d0
        
        radgTiTv = 0d0
        radgTpTv = 0d0
        radgTfTv = 0d0


        A2L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A2R = AR(chimeigen,sfmeigen,sftype,sfgen)
        A2Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A2Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))


     
! colour class (CF**2*CA - 1/2*CF)
        radgTiTv = (CF**2*CV - CF/2d0)*((GS**4*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*
     &    (2*m2**2 - p1p3 + pap3 - 2*papb + 2*pbp1 - 3*pbp3))/((m1**2 - m2**2 + m3**2 - 2*pap1)**2*pbp3))
     
        radgTpTv = (CF**2*CV - CF/2d0)*((GS**4*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*
     &    (2*m2**2 - p1p2 + pap2 - 2*papb + 2*pbp1 - 3*pbp2))/((m1**2 - m2**2 + m3**2 - 2*pap1)**2*pbp2))
     
        radgTfTv = (CF**2*CV - CF/2d0)*(-(GS**4*(A2Lc*(A2R*m1*m3*(2*p1p2 + p2p3 - 4*pbp1 - 2*pbp3) + 
     &          A2L*(2*p1p2*pap1 + p2p3*pap1 - p1p3*pap2 + p1p2*pap3 + 2*p1p3*papb - 4*pap1*pbp1 - 2*pap3*pbp1 - 2*pap1*pbp3))
     &         + A2Rc*(A2L*m1*m3*(2*p1p2 + p2p3 - 4*pbp1 - 2*pbp3) + 
     &          A2R*(2*p1p2*pap1 + p2p3*pap1 - p1p3*pap2 + p1p2*pap3 + 2*p1p3*papb - 4*pap1*pbp1 - 2*pap3*pbp1 - 2*pap1*pbp3))
     &       ))/(2.*p1p3*(m1**2 - m2**2 + m3**2 - 2*pap1)*pbp2))
     
! colour class 2*(CF**2*CA - 1/2*CF)
        radgTvTv = 2d0*(CF**2*CV - CF/2d0)*((8*GS**4*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1))/(m1**2 - m2**2 + m3**2 - 2*pap1)**2)
        
!         print*,"2d0*(CF**2*CV - CF/2d0)",(2d0*(CF**2*CV - CF/2d0))
!         print*,"8*GS**4",(8*GS**4)
!         print*,"A2Lc*A2R*m1*m3",(A2Lc*A2R*m1*m3)
!         print*,"A2L*A2Rc*m1*m3",(A2L*A2Rc*m1*m3)
!         print*,"A2L*A2Lc*pap1",(A2L*A2Lc*pap1)
!         print*,"A2R*A2Rc*pap1",(A2R*A2Rc*pap1)
!         print*,"(m1**2 - m2**2 + m3**2 - 2*pap1)**2",((m1**2 - m2**2 + m3**2 - 2*pap1)**2)
!         print*,"pap1",pap1
!         print*,"p1p2",p1p2
!         print*,"pap2",pap2
!         print*,"papb",papb
!         print*,"pbp1",pbp1
!         print*,"pbp2",pbp2

     
        NeuQ2qx_M2RadiationGTTcon = 0.5*(radgTvTv + 2d0*(radgTiTv + radgTpTv + radgTfTv))


!          print*,"fsoft",fsoft
!          print*,"spflag",spflag
!          print*,"radgTiTv",radgTiTv
!          print*,"radgTpTv",radgTpTv
!          print*,"radgTfTv",radgTfTv
!          print*,"radgTvTv",radgTvTv
!         
!         print*,"NeuQ2qx_M2RadiationGTTcon",NeuQ2qx_M2RadiationGTTcon


    
        end
        

! =================================================================
        double precision function NeuQ2qx_M2RadiationGqq()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double complex radgSqSq, radgTqTq, radgSqTq, radgTqgTqg, radgSqTqg, radgTqTqg
        double complex A1L, A1Lc, A1R, A1Rc, A2L, A2Lc, A2R, A2Rc, A3L, A3Lc, A3R, A3Rc
        double complex J3L, J3Lc, J3R, J3Rc, J5L, J5Lc, J5R, J5Rc
        double precision mq3, mq3c
        integer ttf,ttgen,k,l

! functions from DMNLO_Kinematics.F
	  double precision Krondelta


        NeuQ2qx_M2RadiationGqq = 0d0
       

        radgSqSq = 0d0
        radgTqTq = 0d0
        radgSqTq = 0d0
        radgTqgTqg = 0d0
        radgSqTqg = 0d0
        radgTqTqg = 0d0
        
        do ttf = 3, 4
           do ttgen = 1, 2
             

        A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A1R = AR(chimeigen,sfmeigen,sftype,sfgen)
        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        A2L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A2R = AR(chimeigen,sfmeigen,sftype,sfgen)
        A2Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A2Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
       
        

        radgSqSq = radgSqSq + ((CV**2 - 1)/4d0)*(4*GS**4*(A1Rc*(A1R*(m3**2*p1p3*pap2 + m3**2*p1p2*pap3 - m2**2*(p1p3*pap2 + p1p2*pap3) - 2*m3**2*p2p3*papb + 2*p1p3*papb*pbp2 + 2*p1p2*papb*pbp3 + 
     &            m1**2*(-2*m3**2*p2p3 + p1p3*pap2 + p1p2*pap3 + 2*p1p3*pbp2 + 2*p1p2*pbp3)) - 
     &         A1L*m1*m3*(p2p3*(m1**2 + m2**2 + m3**2 + 2*papb) - 2*(p1p3*(pap2 + pbp2) + p1p2*(pap3 + pbp3)))) + 
     &      A1Lc*(A1L*(m3**2*p1p3*pap2 + m3**2*p1p2*pap3 - m2**2*(p1p3*pap2 + p1p2*pap3) - 2*m3**2*p2p3*papb + 2*p1p3*papb*pbp2 + 2*p1p2*papb*pbp3 + 
     &            m1**2*(-2*m3**2*p2p3 + p1p3*pap2 + p1p2*pap3 + 2*p1p3*pbp2 + 2*p1p2*pbp3)) - 
     &         A1R*m1*m3*(p2p3*(m1**2 + m2**2 + m3**2 + 2*papb) - 2*(p1p3*(pap2 + pbp2) + p1p2*(pap3 + pbp3))))))/(p2p3**2*(m1**2 + m2**2 - m3**2 + 2*papb)**2)
            
        
        radgSqTq = radgSqTq + ((CV**2 - 1)/4d0)*(2*GS**4*(A1L*(A2Rc*m1*m3*(m1**2*p2p3 - 2*pap2*pap3 - p2p3*(m2**2 + m3**2 + 2*pbp1) + 2*(p1p2 + pbp2)*(p1p3 + pbp3)) + 
     &         A2Lc*(-(m2**2*p2p3*pap1) - m3**2*p2p3*pap1 + m3**2*p1p3*pap2 + m3**2*p1p2*pap3 - 2*m3**2*pap2*pap3 + 2*p1p2*p1p3*papb - 2*m3**2*p2p3*papb + 
     &            2*p2p3*pap1*papb - p1p3*pap2*papb - p1p2*pap3*papb - p1p3*pap2*pbp1 - p1p2*pap3*pbp1 + 2*pap2*pap3*pbp1 + p1p3*pap1*pbp2 + m3**2*pap3*pbp2 - 
     &            pap1*pap3*pbp2 + p1p3*papb*pbp2 - pap3*pbp1*pbp2 + (m3**2*pap2 + p1p2*(pap1 + papb) - pap2*(pap1 + pbp1) + 2*pap1*pbp2)*pbp3 + 
     &            m1**2*(p2p3*pap1 - p1p3*pap2 - 2*p2p3*pbp1 + p1p3*pbp2 + p1p2*(2*p1p3 - pap3 + pbp3)))) + 
     &      A1R*(A2Lc*m1*m3*(m1**2*p2p3 - 2*pap2*pap3 - p2p3*(m2**2 + m3**2 + 2*pbp1) + 2*(p1p2 + pbp2)*(p1p3 + pbp3)) + 
     &         A2Rc*(-(m2**2*p2p3*pap1) - m3**2*p2p3*pap1 + m3**2*p1p3*pap2 + m3**2*p1p2*pap3 - 2*m3**2*pap2*pap3 + 2*p1p2*p1p3*papb - 2*m3**2*p2p3*papb + 
     &            2*p2p3*pap1*papb - p1p3*pap2*papb - p1p2*pap3*papb - p1p3*pap2*pbp1 - p1p2*pap3*pbp1 + 2*pap2*pap3*pbp1 + p1p3*pap1*pbp2 + m3**2*pap3*pbp2 - 
     &            pap1*pap3*pbp2 + p1p3*papb*pbp2 - pap3*pbp1*pbp2 + (m3**2*pap2 + p1p2*(pap1 + papb) - pap2*(pap1 + pbp1) + 2*pap1*pbp2)*pbp3 + 
     &            m1**2*(p2p3*pap1 - p1p3*pap2 - 2*p2p3*pbp1 + p1p3*pbp2 + p1p2*(2*p1p3 - pap3 + pbp3))))))/
     &  (p2p3**2*(m1**2 - m2**2 + m3**2 - 2*pap1)*(m1**2 + m2**2 - m3**2 + 2*papb))
        
        
        radgTqTq = radgTqTq + ((CV**2 - 1)/4d0)*(-2*GS**4*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*
     &    (p2p3*(m1**2 + m2**2 + m3**2 - 2*(pap1 + papb - pbp1)) - 2*(p1p2 - pap2 + pbp2)*(p1p3 - pap3 + pbp3)))/(p2p3**2*(m1**2 - m2**2 + m3**2 - 2*pap1)**2)

         do k = 1, 2
              do l = 1, 2
              
              
        A3L = AL(chimeigen,k,ttf,ttgen)
        A3R = AR(chimeigen,k,ttf,ttgen)
        A3Lc = dconjg(AL(chimeigen,l,ttf,ttgen))
        A3Rc = dconjg(AR(chimeigen,l,ttf,ttgen))
        
        J3L = JL(sfmeigen,sftype,sfgen)
        J3R = JL(sfmeigen,sftype,sfgen)
        J3Lc = dconjg(JL(sfmeigen,sftype,sfgen))
        J3Rc = dconjg(JL(sfmeigen,sftype,sfgen))
        
        J5L = JL(k,ttf,ttgen)
        J5R = JR(k,ttf,ttgen)
        J5Lc = dconjg(JL(l,ttf,ttgen))
        J5Rc = dconjg(JR(l,ttf,ttgen))
      
        mq3 = MSf(k,ttf,ttgen)
        mq3c = MSf(l,ttf,ttgen)
                 
                 
                 
        radgTqgTqg = radgTqgTqg + ((CV**2 - 1)/4d0)*(-4*(A3L*A3Lc + A3R*A3Rc)*pap2*(((J3L*J3Lc*J5L*J5Lc + J3R*J3Rc*J5R*J5Rc)*(m2 - m3)*(m2 + m3) + (J3Lc*J3R + J3L*J3Rc)*(J5L*J5Lc + J5R*J5Rc)*m3*MGl - 
     &         (J3R*J3Rc*J5L*J5Lc + J3L*J3Lc*J5R*J5Rc)*MGl**2)*p1p3 + 
     &      (m3*(2*(J3L*J3Lc*J5L*J5Lc + J3R*J3Rc*J5R*J5Rc)*m3 - (J3Lc*J3R + J3L*J3Rc)*(J5L*J5Lc + J5R*J5Rc)*MGl) - 2*(J3L*J3Lc*J5L*J5Lc + J3R*J3Rc*J5R*J5Rc)*pbp1)*
     &       pbp3))/((m1**2 - mq3**2 - 2*pap2)*(m1**2 - mq3c**2 - 2*pap2)*(m2**2 + m3**2 - MGl**2 - 2*pbp1)**2)
        
        
        radgSqTqg = radgSqTqg + Krondelta(k,l)*((CV**2 - 1)/4d0)*(2*GS**2*(A1R*(A3Lc*J5Rc*m1*(J3Lc*MGl*(m3**2*p2p3 - 2*p1p3*(pap2 + pbp2)) + 
     &            J3Rc*m3*(-(m3**2*p2p3) + p2p3*pbp1 + p1p3*(2*pap2 + pbp2) + p1p2*pbp3 - 2*(pap2 + pbp2)*pbp3)) + 
     &         A3Rc*J5Lc*(J3Rc*m3*MGl*(-2*p1p3*pap2 + p2p3*(m1**2 + papb) - pap3*pbp2 + pap2*pbp3) + 
     &            J3Lc*(m2**2*(p2p3*pap1 - p1p3*pap2 - p1p2*pap3) + m3**2*(2*p1p3*pap2 - p2p3*papb + pap3*pbp2) + 
     &               (-3*m3**2*pap2 + 2*p1p2*papb + 2*pap2*pbp1 - 2*pap1*pbp2)*pbp3 + m1**2*(-(m3**2*p2p3) + p2p3*pbp1 - p1p3*pbp2 + p1p2*pbp3)))) + 
     &      A1L*(A3Rc*J5Lc*m1*(J3Rc*MGl*(m3**2*p2p3 - 2*p1p3*(pap2 + pbp2)) + 
     &            J3Lc*m3*(-(m3**2*p2p3) + p2p3*pbp1 + p1p3*(2*pap2 + pbp2) + p1p2*pbp3 - 2*(pap2 + pbp2)*pbp3)) + 
     &         A3Lc*J5Rc*(J3Lc*m3*MGl*(-2*p1p3*pap2 + p2p3*(m1**2 + papb) - pap3*pbp2 + pap2*pbp3) + 
     &            J3Rc*(m2**2*(p2p3*pap1 - p1p3*pap2 - p1p2*pap3) + m3**2*(2*p1p3*pap2 - p2p3*papb + pap3*pbp2) + 
     &               (-3*m3**2*pap2 + 2*p1p2*papb + 2*pap2*pbp1 - 2*pap1*pbp2)*pbp3 + m1**2*(-(m3**2*p2p3) + p2p3*pbp1 - p1p3*pbp2 + p1p2*pbp3))))))/
     &  (p2p3*(m1**2 - mq3c**2 - 2*pap2)*(m1**2 + m2**2 - m3**2 + 2*papb)*(m2**2 + m3**2 - MGl**2 - 2*pbp1))
        
        
        radgTqTqg =  radgTqTqg + Krondelta(k,l)*((CV**2 - 1)/4d0)*(-1d0)*((GS**2*(A2R*(A3Lc*J5Rc*m1*(J3Rc*m3*(p1p3*pap2 + p1p2*(-2*p1p3 + pap3) + p2p3*(-m2**2 + m3**2 - pap1 + papb) - pap3*pbp2 - pap2*pbp3 + 2*pbp2*pbp3) + 
     &              J3Lc*MGl*(-(p2p3*(m3**2 - pap1 + pbp1)) + p1p3*(-pap2 + pbp2) + p1p2*(2*p1p3 - pap3 + pbp3))) + 
     &           A3Rc*J5Lc*(J3Rc*m3*MGl*(m1**2*p2p3 + p1p3*pap2 - p2p3*(pap1 + papb) + pap3*(p1p2 - 2*pap2 + pbp2) + pap2*pbp3) - 
     &              J3Lc*(-(m3**2*p2p3*pap1) + m3**2*p1p3*pap2 + m3**2*p1p2*pap3 - 2*m3**2*pap2*pap3 + m2**2*(p2p3*pap1 + p1p3*pap2 - p1p2*pap3) + 
     &                 2*p1p2*p1p3*papb - 2*m3**2*p2p3*papb - 2*p1p3*pap2*papb + 2*p2p3*pap1*pbp1 - 2*p1p2*pap3*pbp1 + 2*pap2*pap3*pbp1 - 2*p1p3*pap1*pbp2 + 
     &                 2*m3**2*pap3*pbp2 + 2*(pap1*pap2 + p1p2*papb - pap2*pbp1 - pap1*pbp2)*pbp3 + m1**2*(m3**2*p2p3 - p2p3*pbp1 + p1p3*pbp2 - p1p2*pbp3)))) + 
     &        A2L*(A3Rc*J5Lc*m1*(J3Lc*m3*(p1p3*pap2 + p1p2*(-2*p1p3 + pap3) + p2p3*(-m2**2 + m3**2 - pap1 + papb) - pap3*pbp2 - pap2*pbp3 + 2*pbp2*pbp3) + 
     &              J3Rc*MGl*(-(p2p3*(m3**2 - pap1 + pbp1)) + p1p3*(-pap2 + pbp2) + p1p2*(2*p1p3 - pap3 + pbp3))) + 
     &           A3Lc*J5Rc*(J3Lc*m3*MGl*(m1**2*p2p3 + p1p3*pap2 - p2p3*(pap1 + papb) + pap3*(p1p2 - 2*pap2 + pbp2) + pap2*pbp3) - 
     &              J3Rc*(-(m3**2*p2p3*pap1) + m3**2*p1p3*pap2 + m3**2*p1p2*pap3 - 2*m3**2*pap2*pap3 + m2**2*(p2p3*pap1 + p1p3*pap2 - p1p2*pap3) + 
     &                 2*p1p2*p1p3*papb - 2*m3**2*p2p3*papb - 2*p1p3*pap2*papb + 2*p2p3*pap1*pbp1 - 2*p1p2*pap3*pbp1 + 2*pap2*pap3*pbp1 - 2*p1p3*pap1*pbp2 + 
     &                 2*m3**2*pap3*pbp2 + 2*(pap1*pap2 + p1p2*papb - pap2*pbp1 - pap1*pbp2)*pbp3 + m1**2*(m3**2*p2p3 - p2p3*pbp1 + p1p3*pbp2 - p1p2*pbp3))))))/
     &    (p2p3*(m1**2 - m2**2 + m3**2 - 2*pap1)*(m1**2 - mq3c**2 - 2*pap2)*(m2**2 + m3**2 - MGl**2 - 2*pbp1)))

!                    print*,'A1L',A1L
!                    print*,'A1R',A1R
!                    print*,'A1Lc',A1Lc
!                    print*,'A1Rc',A1Rc
!                    print*,'A2L',A2L
!                    print*,'A2R',A2R
!                    print*,'A2Lc',A2Lc
!                    print*,'A2Rc',A2Rc
!                    print*,'A3L',A3L
!                    print*,'A3R',A3R
!                    print*,'A3Lc',A3Lc
!                    print*,'A3Rc',A3Rc
!                    print*,'J3L',J3L
!                    print*,'J3R',J3R
!                    print*,'J3Lc',J3Lc
!                    print*,'J3Rc',J3Rc
!                    print*,'J5L',J5L
!                    print*,'J5R',J5R
!                    print*,'J5Lc',J5Lc
!                    print*,'J5Rc',J5Rc
!                    print*,'Mf(3,1)',Mf(3,1)
!                    print*,'Mf(3,2)',Mf(3,2)
!                    print*,'Mf(4,1)',Mf(4,1)
!                    print*,'Mf(4,2)',Mf(4,2)
!                    
!                    print*,'ttf',ttf
!                    print*,'ttgen',ttgen
!                    print*,'k',k
!                    print*,'l',l
!                    print*,'radgSqSq',radgSqSq
!                    print*,'radgSqTq',radgSqTq
!                    print*,'radgTqTq',radgTqTq
!                    print*,'radgTqgTqg',radgTqgTqg
!                    print*,'radgSqTqg',radgSqTqg
!                    print*,'radgTqTqg',radgTqTqg
     
                 enddo
              enddo
           enddo
        enddo


        NeuQ2qx_M2RadiationGqq = NeuQ2qx_M2RadiationGqq + dreal(radgSqSq + radgTqTq + radgTqgTqg + 2d0*(radgSqTq + radgSqTqg + radgTqTqg))
        
!         print *,"NeuQ2qx_M2RadiationGqq",NeuQ2qx_M2RadiationGqq



       end