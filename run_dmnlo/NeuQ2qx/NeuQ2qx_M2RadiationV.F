! =================================================================
        double precision function NeuQ2qx_M2RadiationVSS()

        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double precision radvSiSi, radvSiSp, radvSiSf, radvSpSp, radvSpSf, radvSfSf
        double complex A1L, A1Lc, A1R, A1Rc, B1L, B1Lc, B1R, B1Rc
        double precision mq1, mq1c
        double precision xv, colorfactor
        double precision wquark1, wquark2
        double precision radvSiSiOS, radvSiSpOS, radvSpSpOS

        if (xtype.eq.1) then
        
           xv = 0d0
           colorfactor = 4d0/3d0
           
        elseif (xtype.eq.2) then
        
           xv = 0d0
           colorfactor = 1d0
           
        elseif (xtype.eq.3) then
        
           xv = 1/(m4**2)
           colorfactor = 1d0
           
        elseif (xtype.eq.4) then

           xv = 1/(m4**2)
           colorfactor = 1d0
                   
        endif

        NeuQ2qx_M2RadiationVSS = 0d0

        radvSiSi = 0d0
        radvSiSp = 0d0
        radvSiSf = 0d0
        radvSpSp = 0d0
        radvSpSf = 0d0
        radvSfSf = 0d0

        A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A1R = AR(chimeigen,sfmeigen,sftype,sfgen)
        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        B1L = BL(xtype,ftype)
        B1R = BR(xtype,ftype)
        B1Lc = dconjg(BL(xtype,ftype))
        B1Rc = dconjg(BR(xtype,ftype))

        mq1 = Mf(pstype,psgen) 
        mq1c = Mf(pstype,psgen)



        if(xtype.eq.4) then

        wquark1 = mq1*1e-3
        wquark2 = width(pstype)

        else

        wquark1 = 0d0
        wquark2 = 0d0

        endif


        radvSiSi = dreal((2D0*GS**2*(kbk3 - 1D0*m2**2)*(A1R*(A1Lc*m1*(B1R*mq1*(2D0*B1Rc*(k1k2 + m3**2) + B1Rc*(2D0*k1k2**2 + (k1k2 - 1D0*m3**2)*m4**2)*xv + B1Lc*m3*mq1c*(-4D0 + m4**2*xv)) + 
     &             B1L*(B1Rc*m3*(2D0*k1k2 + m3**2 + m4**2)*(-4D0 + m4**2*xv) + B1Lc*mq1c*(2D0*(k1k2 + m3**2) + (2D0*k1k2**2 + (k1k2 - 1D0*m3**2)*m4**2)*xv))) + 
     &          A1Rc*(B1R*mq1*(B1Lc*m3*(pap1 + pap2)*(-4D0 + m4**2*xv) + B1Rc*mq1c*(2D0*pap1 - 1D0*m4**2*pap1*xv + 2D0*k1k2*pap2*xv)) + 
     &             B1L*(B1Rc*m3*mq1c*(pap1 + pap2)*(-4D0 + m4**2*xv) + B1Lc*
     &                 (-2D0*m4**2*pap1 + 4D0*k1k2*pap2 + (2D0*k1k2 + m4**2)**2*pap1*xv - 1D0*m3**2*(pap1*(-2D0 + m4**2*xv) + 2D0*pap2*(-2D0 + (k1k2 + m4**2)*xv)))))) + 
     &       A1L*(A1Rc*m1*(B1Rc*mq1c*(2D0*B1R*(k1k2 + m3**2) + B1R*(2D0*k1k2**2 + (k1k2 - 1D0*m3**2)*m4**2)*xv + B1L*m3*mq1*(-4D0 + m4**2*xv)) + 
     &             B1Lc*(B1R*m3*(2D0*k1k2 + m3**2 + m4**2)*(-4D0 + m4**2*xv) + B1L*mq1*(2D0*(k1k2 + m3**2) + (2D0*k1k2**2 + (k1k2 - 1D0*m3**2)*m4**2)*xv))) + 
     &          A1Lc*(B1L*mq1*(B1Rc*m3*(pap1 + pap2)*(-4D0 + m4**2*xv) + B1Lc*mq1c*(2D0*pap1 - 1D0*m4**2*pap1*xv + 2D0*k1k2*pap2*xv)) + 
     &             B1R*(B1Lc*m3*mq1c*(pap1 + pap2)*(-4D0 + m4**2*xv) + B1Rc*
     &                 (-2D0*m4**2*pap1 + 4D0*k1k2*pap2 + (2D0*k1k2 + m4**2)**2*pap1*xv - 1D0*m3**2*(pap1*(-2D0 + m4**2*xv) + 2D0*pap2*(-2D0 + (k1k2 + m4**2)*xv))))))))/
     &   (kbk3**2*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark2)*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1c**2 + cI*mq1c*wquark2)))


        radvSiSiOS = dreal((-1D0*GS**2*(A1Lc*(A1R*m1*mq1*((B1L*B1Lc + B1R*B1Rc)*(m3**4 + m3**2*m4**2 - 2D0*m4**4) - 6D0*(B1Lc*B1R + B1L*B1Rc)*m3*m4**2*mq1 - 
     &             1D0*(B1L*B1Lc + B1R*B1Rc)*(2D0*m3**2 - 1D0*m4**2)*mq1**2 + (B1L*B1Lc + B1R*B1Rc)*mq1**4) + 
     &          2D0*A1L*(B1R*(-3D0*B1Lc*m3*m4**2*mq1 + B1Rc*(m3**4 - 2D0*m4**4 + mq1**4 + m3**2*(m4**2 - 2D0*mq1**2)))*pap1 + 
     &             B1R*(-3D0*B1Lc*m3*m4**2*mq1 + B1Rc*(m3**4 + m3**2*(m4 - 1D0*mq1)*(m4 + mq1) - 2D0*m4**2*(m4 - 1D0*mq1)*(m4 + mq1)))*pap2 + 
     &             B1L*mq1*(-3D0*B1Rc*m3*m4**2*(pap1 + pap2) + B1Lc*mq1*(m4**2*pap1 - 1D0*(m3**2 + m4**2 - 1D0*mq1**2)*pap2)))) + 
     &       A1Rc*(A1L*m1*mq1*((B1L*B1Lc + B1R*B1Rc)*(m3**4 + m3**2*m4**2 - 2D0*m4**4) - 6D0*(B1Lc*B1R + B1L*B1Rc)*m3*m4**2*mq1 - 
     &             1D0*(B1L*B1Lc + B1R*B1Rc)*(2D0*m3**2 - 1D0*m4**2)*mq1**2 + (B1L*B1Lc + B1R*B1Rc)*mq1**4) + 
     &          2D0*A1R*(B1L*(-3D0*B1Rc*m3*m4**2*mq1 + B1Lc*(m3**4 - 2D0*m4**4 + mq1**4 + m3**2*(m4**2 - 2D0*mq1**2)))*pap1 + 
     &             B1L*(-3D0*B1Rc*m3*m4**2*mq1 + B1Lc*(m3**4 + m3**2*(m4 - 1D0*mq1)*(m4 + mq1) - 2D0*m4**2*(m4 - 1D0*mq1)*(m4 + mq1)))*pap2 + 
     &             B1R*mq1*(-3D0*B1Lc*m3*m4**2*(pap1 + pap2) + B1Rc*mq1*(m4**2*pap1 - 1D0*(m3**2 + m4**2 - 1D0*mq1**2)*pap2)))))*
     &     (-2D0*kakb + m1**2 + m2**2 + mq1**2 - 2D0*pap1 - 2D0*pap2 + 2D0*pbp1 + 2D0*pbp2))/
     &   (m4**2*(m1**2 - 1D0*m2**2 + mq1**2 - 2D0*(pap1 + pap2))**2*((2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2)**2 + mq1**2*wquark2**2)))


        radvSiSp = dreal((-1D0*GS**2*(A1R*(A1Lc*m1*(B1R*mq1*(B1Lc*m3*mq1c*(k1k3 + k2k3 + kak3 - 2D0*kakb + kbk3 - 2D0*(m2**2 + pbp1 + pbp2))*(-4D0 + m4**2*xv) + 
     &                B1Rc*(2D0*k1k2*(kak3 - 2D0*kakb + kbk3 - 2D0*m2**2) + 2D0*k1k2**2*(kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*xv + 
     &                   k1k2*((kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*m4**2 + 
     &                      2D0*(pap2*(k1k3 - 2D0*pbp1) + k2k3*(mq1c**2 - 1D0*pap1 - 1D0*pbp1) + (k1k3 - 2D0*mq1c**2 + 2D0*pap1)*pbp2))*xv - 
     &                   1D0*((kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*m3**2 + (mq1c**2 - 1D0*pap2)*(k1k3 - 2D0*pbp1) + k2k3*(pap1 + pbp1) - 1D0*(k1k3 + 2D0*pap1)*pbp2)*
     &                    (-2D0 + m4**2*xv))) + B1L*(B1Rc*m3*((kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*(2D0*k1k2 + m3**2 + m4**2) + mq1c**2*(k1k3 + k2k3 - 2D0*(pbp1 + pbp2)))*
     &                 (-4D0 + m4**2*xv) + B1Lc*mq1c*(2D0*k1k2*(2D0*k2k3 + kak3 - 2D0*kakb + kbk3 - 2D0*m2**2 - 4D0*pbp2) + 
     &                   2D0*k1k2**2*(2D0*k1k3 + kak3 - 2D0*kakb + kbk3 - 2D0*m2**2 - 4D0*pbp1)*xv + 
     &                   k1k2*((kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*m4**2 + 4D0*(-2D0*m4**2 + pap2)*pbp1 + 2D0*k2k3*(-1D0*m3**2 + pap1 + pbp1) + 4D0*(m3**2 - 1D0*pap1)*pbp2 + 
     &                      k1k3*(4D0*m4**2 - 2D0*(pap2 + pbp2)))*xv - 1D0*
     &                    (kak3*m3**2 - 2D0*kakb*m3**2 + kbk3*m3**2 - 2D0*m2**2*m3**2 + k2k3*(2D0*m3**2 - 1D0*pap1 - 1D0*pbp1) - 2D0*m3**2*pbp1 + 2D0*m4**2*pbp1 - 2D0*pap2*pbp1 - 
     &                      4D0*m3**2*pbp2 + 2D0*pap1*pbp2 + k1k3*(m3**2 - 1D0*m4**2 + pap2 + pbp2))*(-2D0 + m4**2*xv)))) + 
     &          A1Rc*(B1R*mq1*(B1Lc*m3*(k1k3*(kakb + m1**2) + k2k3*(kakb + m1**2) + (kbk3 - 2D0*m2**2)*(pap1 + pap2) + kak3*(mq1c**2 - 1D0*pbp1 - 1D0*pbp2) - 
     &                   2D0*(kakb*mq1c**2 + m1**2*(pbp1 + pbp2)))*(-4D0 + m4**2*xv) + 
     &                B1Rc*mq1c*(2D0*k1k2*(kak3 - 2D0*kakb) + 2D0*k1k2**2*(kak3 - 2D0*kakb)*xv + 
     &                   k1k2*((kak3 - 2D0*kakb)*m4**2 + 2D0*k2k3*(kakb + m1**2 - 1D0*pap1) + 2D0*pap2*(k1k3 + kbk3 - 2D0*(m2**2 + pbp1)) - 2D0*(kak3 + 2D0*m1**2 - 2D0*pap1)*pbp2)*
     &                    xv - 1D0*((kak3 - 2D0*kakb)*m3**2 + k1k3*(kakb + m1**2 - 1D0*pap2) - 1D0*(kak3 + 2D0*m1**2 - 2D0*pap2)*pbp1 + pap1*(k2k3 + kbk3 - 2D0*(m2**2 + pbp2)))*
     &                    (-2D0 + m4**2*xv))) + B1L*(B1Rc*m3*mq1c*(2D0*k1k2*(kak3 - 2D0*kakb) + k2k3*kakb + k2k3*m1**2 + k1k3*(kakb + m1**2) - 2D0*kakb*m3**2 - 2D0*kakb*m4**2 + 
     &                   (kbk3 - 2D0*m2**2)*(pap1 + pap2) - 2D0*m1**2*pbp1 + kak3*(m3**2 + m4**2 - 1D0*pbp1 - 1D0*pbp2) - 2D0*m1**2*pbp2)*(-4D0 + m4**2*xv) + 
     &                B1Lc*(2D0*k1k2**2*(2D0*k1k3*(kakb + m1**2) + (kak3 - 2D0*kakb)*mq1c**2 + 2D0*(kbk3 - 2D0*m2**2)*pap1 - 2D0*(kak3 + 2D0*m1**2)*pbp1)*xv - 
     &                   1D0*(2D0*k2k3*(kakb + m1**2)*m3**2 + k1k3*(kakb + m1**2)*(m3 - 1D0*m4)*(m3 + m4) + kak3*m3**2*mq1c**2 - 2D0*kakb*m3**2*mq1c**2 + kbk3*m3**2*pap1 - 
     &                      2D0*m2**2*m3**2*pap1 - 1D0*kbk3*m4**2*pap1 + 2D0*m2**2*m4**2*pap1 - 1D0*k2k3*mq1c**2*pap1 + 2D0*kbk3*m3**2*pap2 - 4D0*m2**2*m3**2*pap2 + 
     &                      k1k3*mq1c**2*pap2 - 1D0*kak3*m3**2*pbp1 - 2D0*m1**2*m3**2*pbp1 + kak3*m4**2*pbp1 + 2D0*m1**2*m4**2*pbp1 - 2D0*mq1c**2*pap2*pbp1 - 2D0*kak3*m3**2*pbp2 - 
     &                      4D0*m1**2*m3**2*pbp2 + 2D0*mq1c**2*pap1*pbp2)*(-2D0 + m4**2*xv) + 
     &                   k1k2*(4D0*k2k3*(kakb + m1**2) + 2D0*(kak3 - 2D0*kakb)*mq1c**2 + 4D0*(kbk3 - 2D0*m2**2)*pap2 - 4D0*(kak3 + 2D0*m1**2)*pbp2 + 
     &                      (-2D0*k2k3*(kakb + m1**2)*m3**2 + 4D0*k1k3*(kakb + m1**2)*m4**2 + kak3*m4**2*mq1c**2 - 2D0*kakb*m4**2*mq1c**2 + 4D0*kbk3*m4**2*pap1 - 
     &                         8D0*m2**2*m4**2*pap1 + 2D0*k2k3*mq1c**2*pap1 - 2D0*kbk3*m3**2*pap2 + 4D0*m2**2*m3**2*pap2 - 2D0*k1k3*mq1c**2*pap2 - 4D0*kak3*m4**2*pbp1 - 
     &                         8D0*m1**2*m4**2*pbp1 + 4D0*mq1c**2*pap2*pbp1 + 2D0*kak3*m3**2*pbp2 + 4D0*m1**2*m3**2*pbp2 - 4D0*mq1c**2*pap1*pbp2)*xv))))) + 
     &       A1L*(A1Rc*m1*(B1Lc*(B1R*m3*((kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*(2D0*k1k2 + m3**2 + m4**2) + mq1c**2*(k1k3 + k2k3 - 2D0*(pbp1 + pbp2)))*(-4D0 + m4**2*xv) + 
     &                B1L*mq1*(2D0*k1k2*(kak3 - 2D0*kakb + kbk3 - 2D0*m2**2) + 2D0*k1k2**2*(kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*xv + 
     &                   k1k2*((kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*m4**2 + 
     &                      2D0*(pap2*(k1k3 - 2D0*pbp1) + k2k3*(mq1c**2 - 1D0*pap1 - 1D0*pbp1) + (k1k3 - 2D0*mq1c**2 + 2D0*pap1)*pbp2))*xv - 
     &                   1D0*((kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*m3**2 + (mq1c**2 - 1D0*pap2)*(k1k3 - 2D0*pbp1) + k2k3*(pap1 + pbp1) - 1D0*(k1k3 + 2D0*pap1)*pbp2)*
     &                    (-2D0 + m4**2*xv))) + B1Rc*mq1c*(B1L*m3*mq1*(k1k3 + k2k3 + kak3 - 2D0*kakb + kbk3 - 2D0*(m2**2 + pbp1 + pbp2))*(-4D0 + m4**2*xv) + 
     &                B1R*(2D0*k1k2*(2D0*k2k3 + kak3 - 2D0*kakb + kbk3 - 2D0*m2**2 - 4D0*pbp2) + 2D0*k1k2**2*(2D0*k1k3 + kak3 - 2D0*kakb + kbk3 - 2D0*m2**2 - 4D0*pbp1)*xv + 
     &                   k1k2*((kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*m4**2 + 4D0*(-2D0*m4**2 + pap2)*pbp1 + 2D0*k2k3*(-1D0*m3**2 + pap1 + pbp1) + 4D0*(m3**2 - 1D0*pap1)*pbp2 + 
     &                      k1k3*(4D0*m4**2 - 2D0*(pap2 + pbp2)))*xv - 1D0*
     &                    (kak3*m3**2 - 2D0*kakb*m3**2 + kbk3*m3**2 - 2D0*m2**2*m3**2 + k2k3*(2D0*m3**2 - 1D0*pap1 - 1D0*pbp1) - 2D0*m3**2*pbp1 + 2D0*m4**2*pbp1 - 2D0*pap2*pbp1 - 
     &                      4D0*m3**2*pbp2 + 2D0*pap1*pbp2 + k1k3*(m3**2 - 1D0*m4**2 + pap2 + pbp2))*(-2D0 + m4**2*xv)))) + 
     &          A1Lc*(B1L*mq1*(B1Rc*m3*(k1k3*(kakb + m1**2) + k2k3*(kakb + m1**2) + (kbk3 - 2D0*m2**2)*(pap1 + pap2) + kak3*(mq1c**2 - 1D0*pbp1 - 1D0*pbp2) - 
     &                   2D0*(kakb*mq1c**2 + m1**2*(pbp1 + pbp2)))*(-4D0 + m4**2*xv) + 
     &                B1Lc*mq1c*(2D0*k1k2*(kak3 - 2D0*kakb) + 2D0*k1k2**2*(kak3 - 2D0*kakb)*xv + 
     &                   k1k2*((kak3 - 2D0*kakb)*m4**2 + 2D0*k2k3*(kakb + m1**2 - 1D0*pap1) + 2D0*pap2*(k1k3 + kbk3 - 2D0*(m2**2 + pbp1)) - 2D0*(kak3 + 2D0*m1**2 - 2D0*pap1)*pbp2)*
     &                    xv - 1D0*((kak3 - 2D0*kakb)*m3**2 + k1k3*(kakb + m1**2 - 1D0*pap2) - 1D0*(kak3 + 2D0*m1**2 - 2D0*pap2)*pbp1 + pap1*(k2k3 + kbk3 - 2D0*(m2**2 + pbp2)))*
     &                    (-2D0 + m4**2*xv))) + B1R*(B1Lc*m3*mq1c*(2D0*k1k2*(kak3 - 2D0*kakb) + k2k3*kakb + k2k3*m1**2 + k1k3*(kakb + m1**2) - 2D0*kakb*m3**2 - 2D0*kakb*m4**2 + 
     &                   (kbk3 - 2D0*m2**2)*(pap1 + pap2) - 2D0*m1**2*pbp1 + kak3*(m3**2 + m4**2 - 1D0*pbp1 - 1D0*pbp2) - 2D0*m1**2*pbp2)*(-4D0 + m4**2*xv) + 
     &                B1Rc*(2D0*k1k2**2*(2D0*k1k3*(kakb + m1**2) + (kak3 - 2D0*kakb)*mq1c**2 + 2D0*(kbk3 - 2D0*m2**2)*pap1 - 2D0*(kak3 + 2D0*m1**2)*pbp1)*xv - 
     &                   1D0*(2D0*k2k3*(kakb + m1**2)*m3**2 + k1k3*(kakb + m1**2)*(m3 - 1D0*m4)*(m3 + m4) + kak3*m3**2*mq1c**2 - 2D0*kakb*m3**2*mq1c**2 + kbk3*m3**2*pap1 - 
     &                      2D0*m2**2*m3**2*pap1 - 1D0*kbk3*m4**2*pap1 + 2D0*m2**2*m4**2*pap1 - 1D0*k2k3*mq1c**2*pap1 + 2D0*kbk3*m3**2*pap2 - 4D0*m2**2*m3**2*pap2 + 
     &                      k1k3*mq1c**2*pap2 - 1D0*kak3*m3**2*pbp1 - 2D0*m1**2*m3**2*pbp1 + kak3*m4**2*pbp1 + 2D0*m1**2*m4**2*pbp1 - 2D0*mq1c**2*pap2*pbp1 - 2D0*kak3*m3**2*pbp2 - 
     &                      4D0*m1**2*m3**2*pbp2 + 2D0*mq1c**2*pap1*pbp2)*(-2D0 + m4**2*xv) + 
     &                   k1k2*(4D0*k2k3*(kakb + m1**2) + 2D0*(kak3 - 2D0*kakb)*mq1c**2 + 4D0*(kbk3 - 2D0*m2**2)*pap2 - 4D0*(kak3 + 2D0*m1**2)*pbp2 + 
     &                      (-2D0*k2k3*(kakb + m1**2)*m3**2 + 4D0*k1k3*(kakb + m1**2)*m4**2 + kak3*m4**2*mq1c**2 - 2D0*kakb*m4**2*mq1c**2 + 4D0*kbk3*m4**2*pap1 - 
     &                         8D0*m2**2*m4**2*pap1 + 2D0*k2k3*mq1c**2*pap1 - 2D0*kbk3*m3**2*pap2 + 4D0*m2**2*m3**2*pap2 - 2D0*k1k3*mq1c**2*pap2 - 4D0*kak3*m4**2*pbp1 - 
     &                         8D0*m1**2*m4**2*pbp1 + 4D0*mq1c**2*pap2*pbp1 + 2D0*kak3*m3**2*pbp2 + 4D0*m1**2*m3**2*pbp2 - 4D0*mq1c**2*pap1*pbp2)*xv)))))))/
     &   (kbk3*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark2)*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1c**2)*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1c**2 + cI*mq1c*wquark2)))


        radvSiSpOS = dreal((-2D0*GS**2*(A1Rc*(A1L*m1*mq1*(B1L*(-6D0*B1Rc*m3*m4**2*mq1 + B1Lc*(m3**4 - 2D0*m4**4 + m4**2*mq1**2 + mq1**4 + m3**2*(m4**2 - 2D0*mq1**2))) + 
     &             B1R*(-6D0*B1Lc*m3*m4**2*mq1 + B1Rc*(m3**4 - 2D0*m4**4 + m4**2*mq1**2 + mq1**4 + m3**2*(m4**2 - 2D0*mq1**2))))*(-1D0*m1**2 + m2**2 + mq1**2 + 2D0*pbp1 + 2D0*pbp2)
     &            + 2D0*A1R*(B1R*mq1*(B1Rc*mq1*(m2**2*m4**2*pap1 + m4**2*mq1**2*pap1 - 1D0*m2**2*m3**2*pap2 - 1D0*m2**2*m4**2*pap2 + m2**2*mq1**2*pap2 - 1D0*m3**2*mq1**2*pap2 - 
     &                   1D0*m4**2*mq1**2*pap2 + mq1**4*pap2 + kakb*(m3**4 - 2D0*m4**4 + mq1**4 - 2D0*mq1**2*pap2 + m3**2*(m4**2 - 2D0*mq1**2 + 2D0*pap2) + 
     &                      m4**2*(mq1**2 - 2D0*pap1 + 2D0*pap2)) - 2D0*m3**2*pap2*pbp1 - 4D0*m4**2*pap2*pbp1 + 2D0*mq1**2*pap2*pbp1 + 
     &                   m1**2*((m3**2 - 1D0*mq1**2)*(pap2 - 2D0*pbp2) + m4**2*(-1D0*pap1 + pap2 + 2D0*pbp1 - 2D0*pbp2)) + 2D0*m3**2*pap1*pbp2 + 4D0*m4**2*pap1*pbp2 - 
     &                   2D0*mq1**2*pap1*pbp2) - 3D0*B1Lc*m3*m4**2*(2D0*kakb*(mq1**2 - 1D0*pap1 - 1D0*pap2) + (m2**2 + mq1**2)*(pap1 + pap2) - 
     &                   1D0*m1**2*(pap1 + pap2 - 2D0*(pbp1 + pbp2)))) + 
     &             B1L*(-3D0*B1Rc*m3*m4**2*mq1*(2D0*kakb*(mq1**2 - 1D0*pap1 - 1D0*pap2) + (m2**2 + mq1**2)*(pap1 + pap2) - 1D0*m1**2*(pap1 + pap2 - 2D0*(pbp1 + pbp2))) + 
     &                B1Lc*(m2**2*m3**4*pap1 + m2**2*m3**2*m4**2*pap1 - 2D0*m2**2*m4**4*pap1 - 2D0*m2**2*m3**2*mq1**2*pap1 + m3**4*mq1**2*pap1 + m3**2*m4**2*mq1**2*pap1 - 
     &                   2D0*m4**4*mq1**2*pap1 + m2**2*mq1**4*pap1 - 2D0*m3**2*mq1**4*pap1 + mq1**6*pap1 + m2**2*m3**4*pap2 + m2**2*m3**2*m4**2*pap2 - 2D0*m2**2*m4**4*pap2 - 
     &                   1D0*m2**2*m3**2*mq1**2*pap2 + m3**4*mq1**2*pap2 + 2D0*m2**2*m4**2*mq1**2*pap2 + m3**2*m4**2*mq1**2*pap2 - 2D0*m4**4*mq1**2*pap2 - 1D0*m3**2*mq1**4*pap2 + 
     &                   2D0*m4**2*mq1**4*pap2 + kakb*(mq1**6 - 2D0*mq1**4*pap1 + m4**2*(mq1**4 - 4D0*mq1**2*pap2) + m3**4*(mq1**2 - 2D0*(pap1 + pap2)) + 
     &                      m4**4*(-2D0*mq1**2 + 4D0*(pap1 + pap2)) + m3**2*(2D0*mq1**2*(-1D0*mq1**2 + 2D0*pap1 + pap2) + m4**2*(mq1**2 - 2D0*(pap1 + pap2)))) + 
     &                   2D0*m3**2*mq1**2*pap2*pbp1 + 4D0*m4**2*mq1**2*pap2*pbp1 - 2D0*mq1**4*pap2*pbp1 - 2D0*m3**2*mq1**2*pap1*pbp2 - 4D0*m4**2*mq1**2*pap1*pbp2 + 
     &                   2D0*mq1**4*pap1*pbp2 - 1D0*m1**2*(mq1**4*(pap1 - 2D0*pbp1) + 2D0*m4**2*mq1**2*(pap2 - 2D0*pbp2) + m3**4*(pap1 + pap2 - 2D0*(pbp1 + pbp2)) - 
     &                      2D0*m4**4*(pap1 + pap2 - 2D0*(pbp1 + pbp2)) + m3**2*(m4**2*(pap1 + pap2 - 2D0*(pbp1 + pbp2)) - 1D0*mq1**2*(2D0*pap1 + pap2 - 2D0*(2D0*pbp1 + pbp2)))))))
     &          ) + A1Lc*(A1R*m1*mq1*(B1L*(-6D0*B1Rc*m3*m4**2*mq1 + B1Lc*(m3**4 - 2D0*m4**4 + m4**2*mq1**2 + mq1**4 + m3**2*(m4**2 - 2D0*mq1**2))) + 
     &             B1R*(-6D0*B1Lc*m3*m4**2*mq1 + B1Rc*(m3**4 - 2D0*m4**4 + m4**2*mq1**2 + mq1**4 + m3**2*(m4**2 - 2D0*mq1**2))))*(-1D0*m1**2 + m2**2 + mq1**2 + 2D0*pbp1 + 2D0*pbp2)
     &            + 2D0*A1L*(B1L*mq1*(B1Lc*mq1*(m2**2*m4**2*pap1 + m4**2*mq1**2*pap1 - 1D0*m2**2*m3**2*pap2 - 1D0*m2**2*m4**2*pap2 + m2**2*mq1**2*pap2 - 1D0*m3**2*mq1**2*pap2 - 
     &                   1D0*m4**2*mq1**2*pap2 + mq1**4*pap2 + kakb*(m3**4 - 2D0*m4**4 + mq1**4 - 2D0*mq1**2*pap2 + m3**2*(m4**2 - 2D0*mq1**2 + 2D0*pap2) + 
     &                      m4**2*(mq1**2 - 2D0*pap1 + 2D0*pap2)) - 2D0*m3**2*pap2*pbp1 - 4D0*m4**2*pap2*pbp1 + 2D0*mq1**2*pap2*pbp1 + 
     &                   m1**2*((m3**2 - 1D0*mq1**2)*(pap2 - 2D0*pbp2) + m4**2*(-1D0*pap1 + pap2 + 2D0*pbp1 - 2D0*pbp2)) + 2D0*m3**2*pap1*pbp2 + 4D0*m4**2*pap1*pbp2 - 
     &                   2D0*mq1**2*pap1*pbp2) - 3D0*B1Rc*m3*m4**2*(2D0*kakb*(mq1**2 - 1D0*pap1 - 1D0*pap2) + (m2**2 + mq1**2)*(pap1 + pap2) - 
     &                   1D0*m1**2*(pap1 + pap2 - 2D0*(pbp1 + pbp2)))) + 
     &             B1R*(-3D0*B1Lc*m3*m4**2*mq1*(2D0*kakb*(mq1**2 - 1D0*pap1 - 1D0*pap2) + (m2**2 + mq1**2)*(pap1 + pap2) - 1D0*m1**2*(pap1 + pap2 - 2D0*(pbp1 + pbp2))) + 
     &                B1Rc*(m2**2*m3**4*pap1 + m2**2*m3**2*m4**2*pap1 - 2D0*m2**2*m4**4*pap1 - 2D0*m2**2*m3**2*mq1**2*pap1 + m3**4*mq1**2*pap1 + m3**2*m4**2*mq1**2*pap1 - 
     &                   2D0*m4**4*mq1**2*pap1 + m2**2*mq1**4*pap1 - 2D0*m3**2*mq1**4*pap1 + mq1**6*pap1 + m2**2*m3**4*pap2 + m2**2*m3**2*m4**2*pap2 - 2D0*m2**2*m4**4*pap2 - 
     &                   1D0*m2**2*m3**2*mq1**2*pap2 + m3**4*mq1**2*pap2 + 2D0*m2**2*m4**2*mq1**2*pap2 + m3**2*m4**2*mq1**2*pap2 - 2D0*m4**4*mq1**2*pap2 - 1D0*m3**2*mq1**4*pap2 + 
     &                   2D0*m4**2*mq1**4*pap2 + kakb*(mq1**6 - 2D0*mq1**4*pap1 + m4**2*(mq1**4 - 4D0*mq1**2*pap2) + m3**4*(mq1**2 - 2D0*(pap1 + pap2)) + 
     &                      m4**4*(-2D0*mq1**2 + 4D0*(pap1 + pap2)) + m3**2*(2D0*mq1**2*(-1D0*mq1**2 + 2D0*pap1 + pap2) + m4**2*(mq1**2 - 2D0*(pap1 + pap2)))) + 
     &                   2D0*m3**2*mq1**2*pap2*pbp1 + 4D0*m4**2*mq1**2*pap2*pbp1 - 2D0*mq1**4*pap2*pbp1 - 2D0*m3**2*mq1**2*pap1*pbp2 - 4D0*m4**2*mq1**2*pap1*pbp2 + 
     &                   2D0*mq1**4*pap1*pbp2 - 1D0*m1**2*(mq1**4*(pap1 - 2D0*pbp1) + 2D0*m4**2*mq1**2*(pap2 - 2D0*pbp2) + m3**4*(pap1 + pap2 - 2D0*(pbp1 + pbp2)) - 
     &                      2D0*m4**4*(pap1 + pap2 - 2D0*(pbp1 + pbp2)) + m3**2*(m4**2*(pap1 + pap2 - 2D0*(pbp1 + pbp2)) - 1D0*mq1**2*(2D0*pap1 + pap2 - 2D0*(2D0*pbp1 + pbp2)))))))
     &          )))/(m4**2*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2)*(m1**2 - 1D0*m2**2 + mq1**2 - 2D0*pap1 - 2D0*pap2)*
     &     (4D0*k1k2**2 + m3**4 + m4**4 - 2D0*m4**2*mq1**2 + mq1**4 + 2D0*m3**2*(m4**2 - 1D0*mq1**2) + 4D0*k1k2*(m3**2 + m4**2 - 1D0*mq1**2) + mq1**2*wquark2**2)))


        radvSiSf = dreal((GS**2*(A1L*(A1Rc*m1*(B1L*mq1*(-1D0*B1Rc*m3*mq1c*(k1k3 - 1D0*kbk3 - 2D0*pbp1)*(-4D0 + m4**2*xv) + 
     &                B1Lc*(-2D0*k1k3*(kakb + m2**2 + pap1 + pbp1) + 2D0*kbk3*(pap1 + 2D0*pbp1) + 2D0*pbp1*(kak3 + 2D0*(pap1 + pbp1)) + 
     &                   k1k3*(m4**2*(kakb + m2**2 + pap1 + pbp1) - 2D0*(k1k2 + pbp2)*(pap2 + pbp2))*xv + 
     &                   (kbk3*(-1D0*m4**2*(pap1 + 2D0*pbp1) + 2D0*k1k2*(pap2 + pbp2)) + pbp1*(-1D0*m4**2*(kak3 + 2D0*(pap1 + pbp1)) + 2D0*(2D0*k1k2 + k2k3)*(pap2 + pbp2)))*xv)) + 
     &             B1R*(B1Rc*mq1c*(2D0*(k2k3*pbp1 + m3**2*(kbk3 + 2D0*pbp1) + k1k2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) - 1D0*k1k3*(m3**2 + pbp2)) + 
     &                   (2D0*k1k2**2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) + k1k2*(2D0*k2k3*pbp1 + m4**2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) - 2D0*k1k3*pbp2) + 
     &                      m4**2*(m3**2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + k2k3*pbp1 - 1D0*k1k3*pbp2))*xv) + 
     &                B1Lc*m3*(4D0*(k1k3 - 1D0*kbk3 - 2D0*pbp1)*(pap1 + pap2 + pbp1 + pbp2) + 2D0*k1k2*k2k3*(kakb + m2**2)*xv + 
     &                   (m4**2*(kbk3*(pap1 + pap2 + 2D0*pbp1) + pbp1*(kak3 + 2D0*(pap1 + pap2 + pbp1))) - 1D0*(2D0*k1k2*(kak3 + kbk3) + m4**2*(kak3 - 2D0*pbp1))*pbp2 + 
     &                      k1k3*(-1D0*m4**2*(kakb + m2**2 + pap1 + pap2 + pbp1) - 1D0*(m4**2 - 2D0*pap2)*pbp2 + 2D0*pbp2**2) + k2k3*((kakb + m2**2)*m4**2 - 2D0*pbp1*(pap2 + pbp2))
     &                      )*xv))) + A1Lc*(B1R*(B1Lc*m3*mq1c*(4D0*(pap1 + pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + 2D0*k1k2*k2k3*kakb*xv + 
     &                   (-1D0*k1k3*m4**2*(kakb + pap1 + pap2) + k2k3*(kakb*m4**2 - 2D0*pap2*pbp1) + m4**2*(kak3*pbp1 + (pap1 + pap2)*(kbk3 + 2D0*pbp1)) - 
     &                      1D0*kak3*(2D0*k1k2 + m4**2)*pbp2 + 2D0*k1k3*pap2*pbp2)*xv) + 
     &                B1Rc*(2D0*k1k2*(kbk3*m1**2 - 1D0*k1k3*(kakb + m1**2) + kak3*m2**2 + 2D0*(kakb + m1**2)*pbp1) - 
     &                   2D0*k1k3*(m2**2*(2D0*pap1 + pap2) + kakb*(m3**2 - 2D0*pbp1) + pap2*pbp1 - 1D0*pap1*pbp2 + m1**2*(m3**2 + pbp2)) + 
     &                   2D0*(kbk3*m1**2*m3**2 - 1D0*k2k3*m2**2*pap1 + 2D0*k2k3*kakb*pbp1 + k2k3*m1**2*pbp1 + 2D0*kakb*m3**2*pbp1 + 2D0*m1**2*m3**2*pbp1 + 
     &                      2D0*kbk3*(pap1 + pap2)*pbp1 + 2D0*pap2*pbp1**2 - 2D0*pap1*pbp1*pbp2 + kak3*(m2**2*m3**2 - 2D0*pbp1*(pbp1 + pbp2))) - 
     &                   2D0*k1k2**2*(k1k3*(kakb + m1**2) + kak3*m2**2 - 2D0*kakb*(kbk3 + pbp1) - 1D0*m1**2*(kbk3 + 2D0*pbp1))*xv + 
     &                   k1k3*m4**2*((kakb + m1**2)*m3**2 + pap2*(m2**2 + pbp1) - 1D0*(2D0*kakb + m1**2 + pap1)*pbp2)*xv + 
     &                   (kak3*m2**2*m3**2*m4**2 + k2k3*m2**2*m4**2*pap1 - 2D0*k2k3*m2**2*m3**2*pap2 + k2k3*m1**2*m4**2*pbp1 - 2D0*m1**2*m3**2*m4**2*pbp1 - 
     &                      2D0*m4**2*pap2*pbp1**2 - 1D0*kbk3*m4**2*(m1**2*m3**2 + 2D0*pap2*pbp1) + 2D0*kbk3*m3**2*pap2*pbp2 + 2D0*m4**2*(kak3 + pap1)*pbp1*pbp2 - 
     &                      2D0*kak3*m3**2*pbp2**2 + 2D0*kakb*m3**2*(-1D0*m4**2*(kbk3 + pbp1) + k2k3*pbp2))*xv + 
     &                   k1k2*(kbk3*m1**2*m4**2 - 1D0*kak3*m2**2*m4**2 + 2D0*k2k3*m2**2*pap1 + 2D0*k2k3*m1**2*pbp1 + 2D0*m1**2*m4**2*pbp1 - 4D0*kbk3*pap2*pbp1 - 4D0*pap2*pbp1**2 + 
     &                      2D0*kakb*m4**2*(kbk3 + pbp1) + 4D0*(kak3 + pap1)*pbp1*pbp2 - 
     &                      1D0*k1k3*((kakb + m1**2)*m4**2 - 2D0*pap2*(m2**2 + pbp1) + 2D0*(2D0*kakb + m1**2 + pap1)*pbp2))*xv)) - 
     &             1D0*B1L*mq1*(B1Rc*m3*(4D0*(kakb + m1**2)*(kbk3 + 2D0*pbp1) - 
     &                   1D0*(kbk3*m1**2*m4**2 - 1D0*kak3*m2**2*m4**2 + 2D0*k2k3*m2**2*pap2 + 2D0*m1**2*m4**2*pbp1 + 2D0*kakb*m4**2*(kbk3 + pbp1) - 2D0*k2k3*kakb*pbp2 - 
     &                      2D0*kbk3*pap2*pbp2 + 2D0*kak3*pbp2**2)*xv + k1k3*(kakb + m1**2)*(-4D0 + m4**2*xv)) + 
     &                B1Lc*mq1c*(2D0*k1k3*(kakb + pap1) - 2D0*(kak3 + 2D0*pap1)*pbp1 + (m4**2*(kak3 + 2D0*pap1) - 2D0*(2D0*k1k2 + k2k3)*pap2)*pbp1*xv - 
     &                   1D0*k1k3*(m4**2*(kakb + pap1) - 2D0*pap2*(k1k2 + pbp2))*xv + kbk3*(-2D0*pap1 + m4**2*pap1*xv - 2D0*k1k2*pap2*xv))))) + 
     &       A1R*(A1Lc*m1*(B1R*mq1*(-1D0*B1Lc*m3*mq1c*(k1k3 - 1D0*kbk3 - 2D0*pbp1)*(-4D0 + m4**2*xv) + 
     &                B1Rc*(-2D0*k1k3*(kakb + m2**2 + pap1 + pbp1) + 2D0*kbk3*(pap1 + 2D0*pbp1) + 2D0*pbp1*(kak3 + 2D0*(pap1 + pbp1)) + 
     &                   k1k3*(m4**2*(kakb + m2**2 + pap1 + pbp1) - 2D0*(k1k2 + pbp2)*(pap2 + pbp2))*xv + 
     &                   (kbk3*(-1D0*m4**2*(pap1 + 2D0*pbp1) + 2D0*k1k2*(pap2 + pbp2)) + pbp1*(-1D0*m4**2*(kak3 + 2D0*(pap1 + pbp1)) + 2D0*(2D0*k1k2 + k2k3)*(pap2 + pbp2)))*xv)) + 
     &             B1L*(B1Lc*mq1c*(2D0*(k2k3*pbp1 + m3**2*(kbk3 + 2D0*pbp1) + k1k2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) - 1D0*k1k3*(m3**2 + pbp2)) + 
     &                   (2D0*k1k2**2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) + k1k2*(2D0*k2k3*pbp1 + m4**2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) - 2D0*k1k3*pbp2) + 
     &                      m4**2*(m3**2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + k2k3*pbp1 - 1D0*k1k3*pbp2))*xv) + 
     &                B1Rc*m3*(4D0*(k1k3 - 1D0*kbk3 - 2D0*pbp1)*(pap1 + pap2 + pbp1 + pbp2) + 2D0*k1k2*k2k3*(kakb + m2**2)*xv + 
     &                   (m4**2*(kbk3*(pap1 + pap2 + 2D0*pbp1) + pbp1*(kak3 + 2D0*(pap1 + pap2 + pbp1))) - 1D0*(2D0*k1k2*(kak3 + kbk3) + m4**2*(kak3 - 2D0*pbp1))*pbp2 + 
     &                      k1k3*(-1D0*m4**2*(kakb + m2**2 + pap1 + pap2 + pbp1) - 1D0*(m4**2 - 2D0*pap2)*pbp2 + 2D0*pbp2**2) + k2k3*((kakb + m2**2)*m4**2 - 2D0*pbp1*(pap2 + pbp2))
     &                      )*xv))) + A1Rc*(B1L*(B1Rc*m3*mq1c*(4D0*(pap1 + pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + 2D0*k1k2*k2k3*kakb*xv + 
     &                   (-1D0*k1k3*m4**2*(kakb + pap1 + pap2) + k2k3*(kakb*m4**2 - 2D0*pap2*pbp1) + m4**2*(kak3*pbp1 + (pap1 + pap2)*(kbk3 + 2D0*pbp1)) - 
     &                      1D0*kak3*(2D0*k1k2 + m4**2)*pbp2 + 2D0*k1k3*pap2*pbp2)*xv) + 
     &                B1Lc*(2D0*k1k2*(kbk3*m1**2 - 1D0*k1k3*(kakb + m1**2) + kak3*m2**2 + 2D0*(kakb + m1**2)*pbp1) - 
     &                   2D0*k1k3*(m2**2*(2D0*pap1 + pap2) + kakb*(m3**2 - 2D0*pbp1) + pap2*pbp1 - 1D0*pap1*pbp2 + m1**2*(m3**2 + pbp2)) + 
     &                   2D0*(kbk3*m1**2*m3**2 - 1D0*k2k3*m2**2*pap1 + 2D0*k2k3*kakb*pbp1 + k2k3*m1**2*pbp1 + 2D0*kakb*m3**2*pbp1 + 2D0*m1**2*m3**2*pbp1 + 
     &                      2D0*kbk3*(pap1 + pap2)*pbp1 + 2D0*pap2*pbp1**2 - 2D0*pap1*pbp1*pbp2 + kak3*(m2**2*m3**2 - 2D0*pbp1*(pbp1 + pbp2))) - 
     &                   2D0*k1k2**2*(k1k3*(kakb + m1**2) + kak3*m2**2 - 2D0*kakb*(kbk3 + pbp1) - 1D0*m1**2*(kbk3 + 2D0*pbp1))*xv + 
     &                   k1k3*m4**2*((kakb + m1**2)*m3**2 + pap2*(m2**2 + pbp1) - 1D0*(2D0*kakb + m1**2 + pap1)*pbp2)*xv + 
     &                   (kak3*m2**2*m3**2*m4**2 + k2k3*m2**2*m4**2*pap1 - 2D0*k2k3*m2**2*m3**2*pap2 + k2k3*m1**2*m4**2*pbp1 - 2D0*m1**2*m3**2*m4**2*pbp1 - 
     &                      2D0*m4**2*pap2*pbp1**2 - 1D0*kbk3*m4**2*(m1**2*m3**2 + 2D0*pap2*pbp1) + 2D0*kbk3*m3**2*pap2*pbp2 + 2D0*m4**2*(kak3 + pap1)*pbp1*pbp2 - 
     &                      2D0*kak3*m3**2*pbp2**2 + 2D0*kakb*m3**2*(-1D0*m4**2*(kbk3 + pbp1) + k2k3*pbp2))*xv + 
     &                   k1k2*(kbk3*m1**2*m4**2 - 1D0*kak3*m2**2*m4**2 + 2D0*k2k3*m2**2*pap1 + 2D0*k2k3*m1**2*pbp1 + 2D0*m1**2*m4**2*pbp1 - 4D0*kbk3*pap2*pbp1 - 4D0*pap2*pbp1**2 + 
     &                      2D0*kakb*m4**2*(kbk3 + pbp1) + 4D0*(kak3 + pap1)*pbp1*pbp2 - 
     &                      1D0*k1k3*((kakb + m1**2)*m4**2 - 2D0*pap2*(m2**2 + pbp1) + 2D0*(2D0*kakb + m1**2 + pap1)*pbp2))*xv)) - 
     &             1D0*B1R*mq1*(B1Lc*m3*(4D0*(kakb + m1**2)*(kbk3 + 2D0*pbp1) - 
     &                   1D0*(kbk3*m1**2*m4**2 - 1D0*kak3*m2**2*m4**2 + 2D0*k2k3*m2**2*pap2 + 2D0*m1**2*m4**2*pbp1 + 2D0*kakb*m4**2*(kbk3 + pbp1) - 2D0*k2k3*kakb*pbp2 - 
     &                      2D0*kbk3*pap2*pbp2 + 2D0*kak3*pbp2**2)*xv + k1k3*(kakb + m1**2)*(-4D0 + m4**2*xv)) + 
     &                B1Rc*mq1c*(2D0*k1k3*(kakb + pap1) - 2D0*(kak3 + 2D0*pap1)*pbp1 + (m4**2*(kak3 + 2D0*pap1) - 2D0*(2D0*k1k2 + k2k3)*pap2)*pbp1*xv - 
     &                   1D0*k1k3*(m4**2*(kakb + pap1) - 2D0*pap2*(k1k2 + pbp2))*xv + kbk3*(-2D0*pap1 + m4**2*pap1*xv - 2D0*k1k2*pap2*xv)))))))/
     &   (k1k3*kbk3*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark1)*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1c**2)))


        radvSpSp = dreal((-4D0*GS**2*(A1R*(A1Lc*m1*(B1L*(B1Rc*m3*(2D0*(2D0*kakb + m1**2 + m2**2)*(2D0*k1k2 + m3**2 + m4**2) + 
     &                   mq1**2*(2D0*mq1c**2 - 1D0*pap1 - 1D0*pap2 - 1D0*pbp1 - 1D0*pbp2) - 1D0*mq1c**2*(pap1 + pap2 + pbp1 + pbp2))*(-4D0 + m4**2*xv) + 
     &                B1Lc*mq1c*(4D0*k1k2*(2D0*kakb + m1**2 + m2**2 + mq1**2 - 1D0*pap2 - 1D0*pbp2) + 4D0*k1k2**2*(2D0*kakb + m1**2 + m2**2 + mq1**2 - 1D0*pap1 - 1D0*pbp1)*xv + 
     &                   2D0*k1k2*(2D0*kakb*m4**2 + m1**2*m4**2 + m2**2*m4**2 + m4**2*mq1**2 - 2D0*m4**2*pap1 - 1D0*mq1**2*pap2 - 2D0*m4**2*pbp1 - 1D0*mq1**2*pbp2 + 
     &                      m3**2*(pap2 + pbp2))*xv - 1D0*(4D0*kakb*m3**2 + 2D0*m1**2*m3**2 + 2D0*m2**2*m3**2 + 2D0*m3**2*mq1**2 - 1D0*m3**2*pap1 + m4**2*pap1 - 1D0*mq1**2*pap1 - 
     &                      2D0*m3**2*pap2 - 1D0*m3**2*pbp1 + m4**2*pbp1 - 1D0*mq1**2*pbp1 - 2D0*m3**2*pbp2)*(-2D0 + m4**2*xv))) + 
     &             B1R*mq1*(2D0*B1Lc*m3*mq1c*(2D0*k1k2 + 2D0*kakb + m1**2 + m2**2 + m3**2 + m4**2 - 1D0*pap1 - 1D0*pap2 - 1D0*pbp1 - 1D0*pbp2)*(-4D0 + m4**2*xv) + 
     &                B1Rc*(4D0*k1k2*(2D0*kakb + m1**2 + m2**2 + mq1c**2 - 1D0*pap2 - 1D0*pbp2) + 4D0*k1k2**2*(2D0*kakb + m1**2 + m2**2 + mq1c**2 - 1D0*pap1 - 1D0*pbp1)*xv + 
     &                   2D0*k1k2*(2D0*kakb*m4**2 + m1**2*m4**2 + m2**2*m4**2 + m4**2*mq1c**2 - 2D0*m4**2*pap1 - 1D0*mq1c**2*pap2 - 2D0*m4**2*pbp1 - 1D0*mq1c**2*pbp2 + 
     &                      m3**2*(pap2 + pbp2))*xv - 1D0*(4D0*kakb*m3**2 + 2D0*m1**2*m3**2 + 2D0*m2**2*m3**2 + 2D0*m3**2*mq1c**2 - 1D0*m3**2*pap1 + m4**2*pap1 - 
     &                      1D0*mq1c**2*pap1 - 2D0*m3**2*pap2 - 1D0*m3**2*pbp1 + m4**2*pbp1 - 1D0*mq1c**2*pbp1 - 2D0*m3**2*pbp2)*(-2D0 + m4**2*xv)))) + 
     &          A1Rc*(B1Lc*(B1R*m3*mq1*(4D0*k1k2*(kakb + m1**2) + 2D0*m1**2*m3**2 + 2D0*m1**2*m4**2 + 2D0*m1**2*mq1c**2 - 1D0*m1**2*pap1 + m2**2*pap1 - 1D0*mq1c**2*pap1 - 
     &                   1D0*m1**2*pap2 + m2**2*pap2 - 1D0*mq1c**2*pap2 - 2D0*m1**2*pbp1 + 2D0*kakb*(m3**2 + m4**2 + mq1c**2 - 1D0*pbp1 - 1D0*pbp2) - 2D0*m1**2*pbp2)*
     &                 (-4D0 + m4**2*xv) + B1L*(4D0*k1k2*(m2**2*pap2 + kakb*(mq1**2 + mq1c**2 - 2D0*pbp2) + m1**2*(mq1**2 + mq1c**2 - 1D0*pap2 - 2D0*pbp2)) + 
     &                   4D0*k1k2**2*(m2**2*pap1 + kakb*(mq1**2 + mq1c**2 - 2D0*pbp1) + m1**2*(mq1**2 + mq1c**2 - 1D0*pap1 - 2D0*pbp1))*xv + 
     &                   2D0*k1k2*(2D0*m2**2*m4**2*pap1 - 1D0*(m2**2*m3**2 + mq1**2*mq1c**2)*pap2 + kakb*m4**2*(mq1**2 + mq1c**2 - 4D0*pbp1) + 2D0*kakb*m3**2*pbp2 + 
     &                      m1**2*(m4**2*(mq1**2 + mq1c**2 - 2D0*pap1 - 4D0*pbp1) + m3**2*(pap2 + 2D0*pbp2)))*xv - 
     &                   1D0*(m2**2*m3**2*pap1 - 1D0*m2**2*m4**2*pap1 - 1D0*mq1**2*mq1c**2*pap1 + 2D0*m2**2*m3**2*pap2 + 
     &                      2D0*kakb*(m4**2*pbp1 + m3**2*(mq1**2 + mq1c**2 - 1D0*pbp1 - 2D0*pbp2)) + 
     &                      m1**2*(m4**2*(pap1 + 2D0*pbp1) + m3**2*(2D0*mq1**2 + 2D0*mq1c**2 - 1D0*pap1 - 2D0*(pap2 + pbp1 + 2D0*pbp2))))*(-2D0 + m4**2*xv))) + 
     &             B1Rc*mq1c*(B1L*m3*(4D0*k1k2*(kakb + m1**2) + 2D0*m1**2*m3**2 + 2D0*m1**2*m4**2 + 2D0*m1**2*mq1**2 - 1D0*m1**2*pap1 + m2**2*pap1 - 1D0*mq1**2*pap1 - 
     &                   1D0*m1**2*pap2 + m2**2*pap2 - 1D0*mq1**2*pap2 - 2D0*m1**2*pbp1 + 2D0*kakb*(m3**2 + m4**2 + mq1**2 - 1D0*pbp1 - 1D0*pbp2) - 2D0*m1**2*pbp2)*
     &                 (-4D0 + m4**2*xv) + B1R*mq1*(4D0*k1k2**2*(2D0*(kakb + m1**2) - 1D0*pap1)*xv - 
     &                   1D0*((m2**2 - 1D0*m3**2 + m4**2)*pap1 - 2D0*m3**2*pap2 + kakb*(4D0*m3**2 - 2D0*pbp1) + m1**2*(4D0*m3**2 - 1D0*pap1 - 2D0*pbp1))*(-2D0 + m4**2*xv) + 
     &                   2D0*k1k2*(-2D0*m4**2*pap1*xv + pap2*(-2D0 + (m2**2 + m3**2)*xv) + 2D0*kakb*(2D0 + m4**2*xv - 1D0*pbp2*xv) + 
     &                      m1**2*(4D0 + 2D0*m4**2*xv - 1D0*(pap2 + 2D0*pbp2)*xv)))))) + 
     &       A1L*(A1Rc*m1*(B1Rc*mq1c*(2D0*B1L*m3*mq1*(2D0*k1k2 + 2D0*kakb + m1**2 + m2**2 + m3**2 + m4**2 - 1D0*pap1 - 1D0*pap2 - 1D0*pbp1 - 1D0*pbp2)*(-4D0 + m4**2*xv) + 
     &                B1R*(4D0*k1k2*(2D0*kakb + m1**2 + m2**2 + mq1**2 - 1D0*pap2 - 1D0*pbp2) + 4D0*k1k2**2*(2D0*kakb + m1**2 + m2**2 + mq1**2 - 1D0*pap1 - 1D0*pbp1)*xv + 
     &                   2D0*k1k2*(2D0*kakb*m4**2 + m1**2*m4**2 + m2**2*m4**2 + m4**2*mq1**2 - 2D0*m4**2*pap1 - 1D0*mq1**2*pap2 - 2D0*m4**2*pbp1 - 1D0*mq1**2*pbp2 + 
     &                      m3**2*(pap2 + pbp2))*xv - 1D0*(4D0*kakb*m3**2 + 2D0*m1**2*m3**2 + 2D0*m2**2*m3**2 + 2D0*m3**2*mq1**2 - 1D0*m3**2*pap1 + m4**2*pap1 - 1D0*mq1**2*pap1 - 
     &                      2D0*m3**2*pap2 - 1D0*m3**2*pbp1 + m4**2*pbp1 - 1D0*mq1**2*pbp1 - 2D0*m3**2*pbp2)*(-2D0 + m4**2*xv))) + 
     &             B1Lc*(B1R*m3*(2D0*(2D0*kakb + m1**2 + m2**2)*(2D0*k1k2 + m3**2 + m4**2) + mq1**2*(2D0*mq1c**2 - 1D0*pap1 - 1D0*pap2 - 1D0*pbp1 - 1D0*pbp2) - 
     &                   1D0*mq1c**2*(pap1 + pap2 + pbp1 + pbp2))*(-4D0 + m4**2*xv) + 
     &                B1L*mq1*(4D0*k1k2*(2D0*kakb + m1**2 + m2**2 + mq1c**2 - 1D0*pap2 - 1D0*pbp2) + 4D0*k1k2**2*(2D0*kakb + m1**2 + m2**2 + mq1c**2 - 1D0*pap1 - 1D0*pbp1)*xv + 
     &                   2D0*k1k2*(2D0*kakb*m4**2 + m1**2*m4**2 + m2**2*m4**2 + m4**2*mq1c**2 - 2D0*m4**2*pap1 - 1D0*mq1c**2*pap2 - 2D0*m4**2*pbp1 - 1D0*mq1c**2*pbp2 + 
     &                      m3**2*(pap2 + pbp2))*xv - 1D0*(4D0*kakb*m3**2 + 2D0*m1**2*m3**2 + 2D0*m2**2*m3**2 + 2D0*m3**2*mq1c**2 - 1D0*m3**2*pap1 + m4**2*pap1 - 
     &                      1D0*mq1c**2*pap1 - 2D0*m3**2*pap2 - 1D0*m3**2*pbp1 + m4**2*pbp1 - 1D0*mq1c**2*pbp1 - 2D0*m3**2*pbp2)*(-2D0 + m4**2*xv)))) + 
     &          A1Lc*(B1R*(B1Lc*m3*mq1c*(4D0*k1k2*(kakb + m1**2) + 2D0*m1**2*m3**2 + 2D0*m1**2*m4**2 + 2D0*m1**2*mq1**2 - 1D0*m1**2*pap1 + m2**2*pap1 - 1D0*mq1**2*pap1 - 
     &                   1D0*m1**2*pap2 + m2**2*pap2 - 1D0*mq1**2*pap2 - 2D0*m1**2*pbp1 + 2D0*kakb*(m3**2 + m4**2 + mq1**2 - 1D0*pbp1 - 1D0*pbp2) - 2D0*m1**2*pbp2)*
     &                 (-4D0 + m4**2*xv) + B1Rc*(4D0*k1k2*(m2**2*pap2 + kakb*(mq1**2 + mq1c**2 - 2D0*pbp2) + m1**2*(mq1**2 + mq1c**2 - 1D0*pap2 - 2D0*pbp2)) + 
     &                   4D0*k1k2**2*(m2**2*pap1 + kakb*(mq1**2 + mq1c**2 - 2D0*pbp1) + m1**2*(mq1**2 + mq1c**2 - 1D0*pap1 - 2D0*pbp1))*xv + 
     &                   2D0*k1k2*(2D0*m2**2*m4**2*pap1 - 1D0*(m2**2*m3**2 + mq1**2*mq1c**2)*pap2 + kakb*m4**2*(mq1**2 + mq1c**2 - 4D0*pbp1) + 2D0*kakb*m3**2*pbp2 + 
     &                      m1**2*(m4**2*(mq1**2 + mq1c**2 - 2D0*pap1 - 4D0*pbp1) + m3**2*(pap2 + 2D0*pbp2)))*xv - 
     &                   1D0*(m2**2*m3**2*pap1 - 1D0*m2**2*m4**2*pap1 - 1D0*mq1**2*mq1c**2*pap1 + 2D0*m2**2*m3**2*pap2 + 
     &                      2D0*kakb*(m4**2*pbp1 + m3**2*(mq1**2 + mq1c**2 - 1D0*pbp1 - 2D0*pbp2)) + 
     &                      m1**2*(m4**2*(pap1 + 2D0*pbp1) + m3**2*(2D0*mq1**2 + 2D0*mq1c**2 - 1D0*pap1 - 2D0*(pap2 + pbp1 + 2D0*pbp2))))*(-2D0 + m4**2*xv))) + 
     &             B1L*mq1*(B1Rc*m3*(4D0*k1k2*(kakb + m1**2) + 2D0*m1**2*m3**2 + 2D0*m1**2*m4**2 + 2D0*m1**2*mq1c**2 - 1D0*m1**2*pap1 + m2**2*pap1 - 1D0*mq1c**2*pap1 - 
     &                   1D0*m1**2*pap2 + m2**2*pap2 - 1D0*mq1c**2*pap2 - 2D0*m1**2*pbp1 + 2D0*kakb*(m3**2 + m4**2 + mq1c**2 - 1D0*pbp1 - 1D0*pbp2) - 2D0*m1**2*pbp2)*
     &                 (-4D0 + m4**2*xv) + B1Lc*mq1c*(4D0*k1k2**2*(2D0*(kakb + m1**2) - 1D0*pap1)*xv - 
     &                   1D0*((m2**2 - 1D0*m3**2 + m4**2)*pap1 - 2D0*m3**2*pap2 + kakb*(4D0*m3**2 - 2D0*pbp1) + m1**2*(4D0*m3**2 - 1D0*pap1 - 2D0*pbp1))*(-2D0 + m4**2*xv) + 
     &                   2D0*k1k2*(-2D0*m4**2*pap1*xv + pap2*(-2D0 + (m2**2 + m3**2)*xv) + 2D0*kakb*(2D0 + m4**2*xv - 1D0*pbp2*xv) + 
     &                      m1**2*(4D0 + 2D0*m4**2*xv - 1D0*(pap2 + 2D0*pbp2)*xv))))))))/
     &   ((2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2)*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark2)*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1c**2)*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1c**2 + cI*mq1c*wquark2)))

        radvSpSpOS = dreal((4D0*GS**2*(A1Rc*(B1Rc*mq1*(-1D0*A1L*m1*(-6D0*B1L*m3*m4**2*mq1 + B1R*(m3**4 - 2D0*m4**4 + m4**2*mq1**2 + mq1**4 + m3**2*(m4**2 - 2D0*mq1**2)))*
     &              (2D0*kakb + m1**2 + m2**2 + mq1**2 - 1D0*pap1 - 1D0*pap2 - 1D0*pbp1 - 1D0*pbp2) - 
     &             1D0*A1R*(3D0*B1L*m3*m4**2*(-1D0*(m2**2 - 1D0*mq1**2)*(pap1 + pap2) + 2D0*kakb*(-2D0*mq1**2 + pbp1 + pbp2) + 
     &                   m1**2*(-4D0*mq1**2 + pap1 + pap2 + 2D0*pbp1 + 2D0*pbp2)) + 
     &                B1R*mq1*(-1D0*m3**4*pap1 + m2**2*m4**2*pap1 - 1D0*m3**2*m4**2*pap1 + 2D0*m4**4*pap1 + 2D0*m3**2*mq1**2*pap1 - 1D0*mq1**4*pap1 - 1D0*m2**2*m3**2*pap2 - 
     &                   1D0*m3**4*pap2 - 1D0*m2**2*m4**2*pap2 - 1D0*m3**2*m4**2*pap2 + 2D0*m4**4*pap2 + m2**2*mq1**2*pap2 + m3**2*mq1**2*pap2 - 2D0*m4**2*mq1**2*pap2 + 
     &                   2D0*kakb*(m3**4 - 2D0*m4**4 + mq1**4 - 1D0*mq1**2*pbp2 + m3**2*(m4**2 - 2D0*mq1**2 + pbp2) + m4**2*(mq1**2 - 1D0*pbp1 + pbp2)) + 
     &                   m1**2*(2D0*m3**4 - 4D0*m4**4 + mq1**2*(2D0*mq1**2 - 1D0*pap2 - 2D0*pbp2) + m3**2*(2D0*m4**2 - 4D0*mq1**2 + pap2 + 2D0*pbp2) + 
     &                      m4**2*(2D0*mq1**2 - 1D0*pap1 + pap2 - 2D0*pbp1 + 2D0*pbp2))))) + 
     &          B1Lc*(-1D0*A1L*m1*mq1*(-6D0*B1R*m3*m4**2*mq1 + B1L*(m3**4 - 2D0*m4**4 + m4**2*mq1**2 + mq1**4 + m3**2*(m4**2 - 2D0*mq1**2)))*
     &              (2D0*kakb + m1**2 + m2**2 + mq1**2 - 1D0*pap1 - 1D0*pap2 - 1D0*pbp1 - 1D0*pbp2) - 
     &             1D0*A1R*(3D0*B1R*m3*m4**2*mq1*(-1D0*(m2**2 - 1D0*mq1**2)*(pap1 + pap2) + 2D0*kakb*(-2D0*mq1**2 + pbp1 + pbp2) + 
     &                   m1**2*(-4D0*mq1**2 + pap1 + pap2 + 2D0*pbp1 + 2D0*pbp2)) + 
     &                B1L*(m2**2*m3**4*pap1 + m2**2*m3**2*m4**2*pap1 - 2D0*m2**2*m4**4*pap1 - 2D0*m2**2*m3**2*mq1**2*pap1 + m2**2*mq1**4*pap1 - 1D0*m4**2*mq1**4*pap1 + 
     &                   m2**2*m3**4*pap2 + m2**2*m3**2*m4**2*pap2 - 2D0*m2**2*m4**4*pap2 - 1D0*m2**2*m3**2*mq1**2*pap2 + 2D0*m2**2*m4**2*mq1**2*pap2 + m3**2*mq1**4*pap2 + 
     &                   m4**2*mq1**4*pap2 - 1D0*mq1**6*pap2 + 2D0*kakb*(mq1**6 - 1D0*mq1**4*pbp1 + m3**4*(mq1**2 - 1D0*pbp1 - 1D0*pbp2) + 2D0*m4**4*(-1D0*mq1**2 + pbp1 + pbp2) + 
     &                      m4**2*(mq1**4 - 2D0*mq1**2*pbp2) + m3**2*(m4**2*(mq1**2 - 1D0*pbp1 - 1D0*pbp2) + mq1**2*(-2D0*mq1**2 + 2D0*pbp1 + pbp2))) + 
     &                   m1**2*(mq1**4*(2D0*mq1**2 - 1D0*pap1 - 2D0*pbp1) + 2D0*m4**2*mq1**2*(mq1**2 - 1D0*pap2 - 2D0*pbp2) + 
     &                      m3**4*(2D0*mq1**2 - 1D0*pap1 - 1D0*pap2 - 2D0*pbp1 - 2D0*pbp2) + 2D0*m4**4*(-2D0*mq1**2 + pap1 + pap2 + 2D0*pbp1 + 2D0*pbp2) + 
     &                      m3**2*(m4**2*(2D0*mq1**2 - 1D0*pap1 - 1D0*pap2 - 2D0*pbp1 - 2D0*pbp2) + mq1**2*(-4D0*mq1**2 + 2D0*pap1 + pap2 + 4D0*pbp1 + 2D0*pbp2))))))) + 
     &       A1Lc*(B1Lc*mq1*(-1D0*A1R*m1*(-6D0*B1R*m3*m4**2*mq1 + B1L*(m3**4 - 2D0*m4**4 + m4**2*mq1**2 + mq1**4 + m3**2*(m4**2 - 2D0*mq1**2)))*
     &              (2D0*kakb + m1**2 + m2**2 + mq1**2 - 1D0*pap1 - 1D0*pap2 - 1D0*pbp1 - 1D0*pbp2) - 
     &             1D0*A1L*(3D0*B1R*m3*m4**2*(-1D0*(m2**2 - 1D0*mq1**2)*(pap1 + pap2) + 2D0*kakb*(-2D0*mq1**2 + pbp1 + pbp2) + 
     &                   m1**2*(-4D0*mq1**2 + pap1 + pap2 + 2D0*pbp1 + 2D0*pbp2)) + 
     &                B1L*mq1*(-1D0*m3**4*pap1 + m2**2*m4**2*pap1 - 1D0*m3**2*m4**2*pap1 + 2D0*m4**4*pap1 + 2D0*m3**2*mq1**2*pap1 - 1D0*mq1**4*pap1 - 1D0*m2**2*m3**2*pap2 - 
     &                   1D0*m3**4*pap2 - 1D0*m2**2*m4**2*pap2 - 1D0*m3**2*m4**2*pap2 + 2D0*m4**4*pap2 + m2**2*mq1**2*pap2 + m3**2*mq1**2*pap2 - 2D0*m4**2*mq1**2*pap2 + 
     &                   2D0*kakb*(m3**4 - 2D0*m4**4 + mq1**4 - 1D0*mq1**2*pbp2 + m3**2*(m4**2 - 2D0*mq1**2 + pbp2) + m4**2*(mq1**2 - 1D0*pbp1 + pbp2)) + 
     &                   m1**2*(2D0*m3**4 - 4D0*m4**4 + mq1**2*(2D0*mq1**2 - 1D0*pap2 - 2D0*pbp2) + m3**2*(2D0*m4**2 - 4D0*mq1**2 + pap2 + 2D0*pbp2) + 
     &                      m4**2*(2D0*mq1**2 - 1D0*pap1 + pap2 - 2D0*pbp1 + 2D0*pbp2))))) + 
     &          B1Rc*(-1D0*A1R*m1*mq1*(-6D0*B1L*m3*m4**2*mq1 + B1R*(m3**4 - 2D0*m4**4 + m4**2*mq1**2 + mq1**4 + m3**2*(m4**2 - 2D0*mq1**2)))*
     &              (2D0*kakb + m1**2 + m2**2 + mq1**2 - 1D0*pap1 - 1D0*pap2 - 1D0*pbp1 - 1D0*pbp2) - 
     &             1D0*A1L*(3D0*B1L*m3*m4**2*mq1*(-1D0*(m2**2 - 1D0*mq1**2)*(pap1 + pap2) + 2D0*kakb*(-2D0*mq1**2 + pbp1 + pbp2) + 
     &                   m1**2*(-4D0*mq1**2 + pap1 + pap2 + 2D0*pbp1 + 2D0*pbp2)) + 
     &                B1R*(m2**2*m3**4*pap1 + m2**2*m3**2*m4**2*pap1 - 2D0*m2**2*m4**4*pap1 - 2D0*m2**2*m3**2*mq1**2*pap1 + m2**2*mq1**4*pap1 - 1D0*m4**2*mq1**4*pap1 + 
     &                   m2**2*m3**4*pap2 + m2**2*m3**2*m4**2*pap2 - 2D0*m2**2*m4**4*pap2 - 1D0*m2**2*m3**2*mq1**2*pap2 + 2D0*m2**2*m4**2*mq1**2*pap2 + m3**2*mq1**4*pap2 + 
     &                   m4**2*mq1**4*pap2 - 1D0*mq1**6*pap2 + 2D0*kakb*(mq1**6 - 1D0*mq1**4*pbp1 + m3**4*(mq1**2 - 1D0*pbp1 - 1D0*pbp2) + 2D0*m4**4*(-1D0*mq1**2 + pbp1 + pbp2) + 
     &                      m4**2*(mq1**4 - 2D0*mq1**2*pbp2) + m3**2*(m4**2*(mq1**2 - 1D0*pbp1 - 1D0*pbp2) + mq1**2*(-2D0*mq1**2 + 2D0*pbp1 + pbp2))) + 
     &                   m1**2*(mq1**4*(2D0*mq1**2 - 1D0*pap1 - 2D0*pbp1) + 2D0*m4**2*mq1**2*(mq1**2 - 1D0*pap2 - 2D0*pbp2) + 
     &                      m3**4*(2D0*mq1**2 - 1D0*pap1 - 1D0*pap2 - 2D0*pbp1 - 2D0*pbp2) + 2D0*m4**4*(-2D0*mq1**2 + pap1 + pap2 + 2D0*pbp1 + 2D0*pbp2) + 
     &                      m3**2*(m4**2*(2D0*mq1**2 - 1D0*pap1 - 1D0*pap2 - 2D0*pbp1 - 2D0*pbp2) + mq1**2*(-4D0*mq1**2 + 2D0*pap1 + pap2 + 4D0*pbp1 + 2D0*pbp2)))))))))/
     &   (m4**2*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2)**2*((2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2)**2 + mq1**2*wquark2**2)))


        radvSpSf = dreal((-2D0*GS**2*(A1L*(A1Lc*(B1L*mq1*(B1Lc*mq1c*(-2D0*k1k3*(kakb + m1**2 - 2D0*pap1) + 
     &                   2D0*((2D0*(k1k2 + k2k3) - 1D0*kbk3)*pap1 + m3**2*(kakb + m1**2 + pap1 - 1D0*pap2) + kak3*pbp1) - 1D0*k1k3*(2D0*k1k2 + m4**2)*pap2*xv + 
     &                   (2D0*k1k2**2*(kak3 + kakb + m1**2) + m4**2*(-1D0*k2k3*pap1 - 1D0*m3**2*(kak3 + kakb + m1**2 + pap1)) + m3**2*(2D0*k2k3 + m4**2)*pap2 + 
     &                      k1k2*(kak3*m4**2 + 2D0*k2k3*(kakb + m1**2 - 1D0*pap1) + 2D0*pap2*(kbk3 + m3**2 + pbp1) - 2D0*(kak3 + pap1)*pbp2))*xv) + 
     &                B1Rc*m3*(-2D0*(k2k3*kakb + 2D0*kakb*kbk3 + k2k3*m1**2 + kak3*m1**2 + 2D0*kbk3*m1**2 + 2D0*k1k2*(kakb + m1**2) + k1k3*(kakb + m1**2) - 1D0*kak3*m2**2 + 
     &                      2D0*kakb*m3**2 + 2D0*m1**2*m3**2 - 1D0*kbk3*pap1 + 2D0*m1**2*pap1 - 2D0*m2**2*pap1 - 1D0*kbk3*pap2 + kak3*pbp1 + 4D0*kakb*pbp1 + 4D0*m1**2*pbp1 - 
     &                      2D0*pap2*pbp1 + kak3*pbp2 + 2D0*pap1*pbp2) + 
     &                   (2D0*k2k3*(m2**2*pap2 + kakb*(m4**2 - 2D0*pbp2) + m1**2*(m4**2 - 1D0*pap2 - 2D0*pbp2)) + 
     &                      k1k2*((kakb + m1**2)*(2D0*k2k3 + m4**2) - 2D0*kbk3*pap2 + 2D0*kak3*pbp2) + 
     &                      m4**2*(2D0*kbk3*m1**2 + m1**2*m3**2 + m1**2*pap1 - 1D0*m2**2*pap1 - 2D0*kbk3*pap2 + 2D0*m1**2*pbp1 - 1D0*pap2*pbp1 + 
     &                         kakb*(2D0*kbk3 + m3**2 + 2D0*pbp1) + pap1*pbp2 + kak3*(m1**2 - 1D0*m2**2 + 2D0*pbp2)))*xv)) + 
     &             B1R*(B1Lc*m3*mq1c*(-2D0*((kakb + m1**2)*(2D0*k1k2 + k1k3 + k2k3 + 2D0*m3**2) + mq1**2*(kak3 + 2D0*pap1) + kbk3*(pap1 + pap2) - 1D0*(kak3 - 2D0*pap2)*pbp1) + 
     &                   2D0*(kak3 + 2D0*pap1)*pbp2 + (2D0*k2k3*(kakb + m1**2)*m4**2 - 2D0*k2k3*mq1**2*pap2 + 
     &                      k1k2*((kakb + m1**2)*(2D0*k2k3 + m4**2) + 2D0*kbk3*pap2 - 2D0*kak3*pbp2) + 
     &                      m4**2*((kakb + m1**2)*m3**2 + mq1**2*(kak3 + pap1) + pap2*(2D0*kbk3 + pbp1) - 1D0*(2D0*kak3 + pap1)*pbp2))*xv) + 
     &                B1Rc*(2D0*(-2D0*k1k2*m2**2*pap1 - 2D0*k2k3*m2**2*pap1 - 1D0*m2**2*m3**2*pap1 + kbk3*mq1**2*pap1 + m2**2*m3**2*pap2 - 1D0*kak3*mq1**2*pbp1 + 
     &                      kakb*(4D0*(k1k2 + k2k3)*pbp1 + m3**2*(mq1**2 + 2D0*pbp1 - 2D0*pbp2)) + 
     &                      m1**2*(2D0*(k1k2 + k2k3)*(pap1 + 2D0*pbp1) + m3**2*(mq1**2 + pap1 - 1D0*pap2 + 2D0*pbp1 - 2D0*pbp2))) + 
     &                   (-1D0*kak3*m1**2*m3**2*m4**2 - 2D0*kbk3*m1**2*m3**2*m4**2 + kak3*m2**2*m3**2*m4**2 - 1D0*m1**2*m3**2*m4**2*mq1**2 + 
     &                      2D0*k1k2**2*(kak3*(m1 - 1D0*m2)*(m1 + m2) + (kakb + m1**2)*(2D0*kbk3 + mq1**2)) - 1D0*k2k3*m1**2*m4**2*pap1 + k2k3*m2**2*m4**2*pap1 - 
     &                      1D0*m1**2*m3**2*m4**2*pap1 + m2**2*m3**2*m4**2*pap1 + 2D0*k2k3*m1**2*m3**2*pap2 - 2D0*k2k3*m2**2*m3**2*pap2 + m1**2*m3**2*m4**2*pap2 - 
     &                      1D0*m2**2*m3**2*m4**2*pap2 - 2D0*k2k3*m1**2*m4**2*pbp1 - 2D0*m1**2*m3**2*m4**2*pbp1 - 
     &                      1D0*kakb*m4**2*(2D0*k2k3*pbp1 + m3**2*(2D0*kbk3 + mq1**2 + 2D0*pbp1)) + 2D0*(kakb + m1**2)*m3**2*(2D0*k2k3 + m4**2)*pbp2 + 
     &                      k1k2*(2D0*kakb*kbk3*m4**2 + kak3*m1**2*m4**2 + 2D0*kbk3*m1**2*m4**2 - 1D0*kak3*m2**2*m4**2 + 2D0*k2k3*kakb*mq1**2 + 2D0*k2k3*m1**2*mq1**2 - 
     &                         2D0*k2k3*m1**2*pap1 + 2D0*k2k3*m2**2*pap1 + 2D0*m1**2*m3**2*pap2 - 2D0*m2**2*m3**2*pap2 - 2D0*kbk3*mq1**2*pap2 - 4D0*k2k3*kakb*pbp1 - 
     &                         4D0*k2k3*m1**2*pbp1 - 2D0*mq1**2*pap2*pbp1 + 2D0*(2D0*(kakb + m1**2)*m3**2 + mq1**2*(kak3 + pap1))*pbp2))*xv + 
     &                   k1k3*(-2D0*(kakb + m1**2)*mq1**2 + 4D0*(m1 - 1D0*m2)*(m1 + m2)*pap1 + 8D0*(kakb + m1**2)*pbp1 - 
     &                      1D0*(2D0*k1k2 + m4**2)*((m1 - 1D0*m2)*(m1 + m2)*pap2 + 2D0*(kakb + m1**2)*pbp2)*xv)))) + 
     &          A1Rc*m1*(B1Lc*(B1L*mq1*(-2D0*k1k3*(2D0*kakb + m1**2 + m2**2 - 2D0*(pap1 + pbp1)) + 
     &                   2D0*(2D0*(k1k2 + k2k3)*(pap1 + pbp1) + m3**2*(2D0*kakb + m1**2 + m2**2 + pap1 - 1D0*pap2 + pbp1 - 1D0*pbp2)) - 
     &                   1D0*k1k3*(2D0*k1k2 + m4**2)*(pap2 + pbp2)*xv + (2D0*k1k2**2*(kak3 + 2D0*kakb + kbk3 + m1**2 + m2**2) - 1D0*kak3*m3**2*m4**2 - 2D0*kakb*m3**2*m4**2 - 
     &                      1D0*kbk3*m3**2*m4**2 - 1D0*m1**2*m3**2*m4**2 - 1D0*m2**2*m3**2*m4**2 - 1D0*k2k3*m4**2*pap1 - 1D0*m3**2*m4**2*pap1 + 2D0*k2k3*m3**2*pap2 + 
     &                      m3**2*m4**2*pap2 - 1D0*k2k3*m4**2*pbp1 - 1D0*m3**2*m4**2*pbp1 + 2D0*k2k3*m3**2*pbp2 + m3**2*m4**2*pbp2 + 
     &                      k1k2*((kak3 + kbk3)*m4**2 + 2D0*k2k3*(2D0*kakb + m1**2 + m2**2 - 1D0*pap1 - 1D0*pbp1) + 2D0*m3**2*(pap2 + pbp2)))*xv) + 
     &                B1R*m3*(-2D0*(2D0*kakb + m1**2 + m2**2)*(2D0*k1k2 + k1k3 + k2k3 + 2D0*m3**2) - 2D0*mq1**2*(kak3 + kbk3 + 2D0*(pap1 + pbp1)) + 
     &                   (k1k2*(2D0*kakb + m1**2 + m2**2)*(2D0*k2k3 + m4**2) + m4**2*((2D0*kakb + m1**2 + m2**2)*m3**2 + mq1**2*(kak3 + kbk3 + pap1 + pbp1)) + 
     &                      2D0*k2k3*((2D0*kakb + m1**2 + m2**2)*m4**2 - 1D0*mq1**2*(pap2 + pbp2)))*xv)) + 
     &             B1Rc*mq1c*(B1L*m3*mq1*(-2D0*(2D0*k1k2 + k1k3 + k2k3 + kak3 + kbk3 + 2D0*(m3**2 + pap1 + pbp1)) + 
     &                   (k1k2*(2D0*k2k3 + m4**2) + m4**2*(kak3 + kbk3 + m3**2 + pap1 + pbp1) + 2D0*k2k3*(m4**2 - 1D0*pap2 - 1D0*pbp2))*xv) + 
     &                B1R*(4D0*(k1k2 + k2k3)*(pap1 + pbp1) - 2D0*k1k3*(mq1**2 - 2D0*(pap1 + pbp1)) + 
     &                   (2D0*k1k2**2*(kak3 + kbk3 + mq1**2) + k1k2*((kak3 + kbk3)*m4**2 + 2D0*k2k3*(mq1**2 - 1D0*pap1 - 1D0*pbp1)) - 1D0*k2k3*m4**2*(pap1 + pbp1))*xv - 
     &                   1D0*k1k3*(2D0*k1k2 + m4**2)*(pap2 + pbp2)*xv + m3**2*
     &                    (2D0*(mq1**2 + pap1 - 1D0*pap2 + pbp1 - 1D0*pbp2) - 
     &                      1D0*(m4**2*(kak3 + kbk3 + mq1**2 + pap1 - 1D0*pap2 + pbp1 - 1D0*pbp2) - 2D0*(k1k2 + k2k3)*(pap2 + pbp2))*xv))))) + 
     &       A1R*(A1Rc*(B1Rc*mq1c*(B1R*mq1*(-2D0*k1k3*(kakb + m1**2 - 2D0*pap1) + 2D0*((2D0*(k1k2 + k2k3) - 1D0*kbk3)*pap1 + m3**2*(kakb + m1**2 + pap1 - 1D0*pap2) + kak3*pbp1) - 
     &                   1D0*k1k3*(2D0*k1k2 + m4**2)*pap2*xv + (2D0*k1k2**2*(kak3 + kakb + m1**2) + m4**2*(-1D0*k2k3*pap1 - 1D0*m3**2*(kak3 + kakb + m1**2 + pap1)) + 
     &                      m3**2*(2D0*k2k3 + m4**2)*pap2 + k1k2*(kak3*m4**2 + 2D0*k2k3*(kakb + m1**2 - 1D0*pap1) + 2D0*pap2*(kbk3 + m3**2 + pbp1) - 2D0*(kak3 + pap1)*pbp2))*xv) + 
     &                B1L*m3*(-2D0*((kakb + m1**2)*(2D0*k1k2 + k1k3 + k2k3 + 2D0*m3**2) + mq1**2*(kak3 + 2D0*pap1) + kbk3*(pap1 + pap2) - 1D0*(kak3 - 2D0*pap2)*pbp1) + 
     &                   2D0*(kak3 + 2D0*pap1)*pbp2 + (2D0*k2k3*(kakb + m1**2)*m4**2 - 2D0*k2k3*mq1**2*pap2 + 
     &                      k1k2*((kakb + m1**2)*(2D0*k2k3 + m4**2) + 2D0*kbk3*pap2 - 2D0*kak3*pbp2) + 
     &                      m4**2*((kakb + m1**2)*m3**2 + mq1**2*(kak3 + pap1) + pap2*(2D0*kbk3 + pbp1) - 1D0*(2D0*kak3 + pap1)*pbp2))*xv)) + 
     &             B1Lc*(B1R*m3*mq1*(-2D0*(k2k3*kakb + 2D0*kakb*kbk3 + k2k3*m1**2 + kak3*m1**2 + 2D0*kbk3*m1**2 + 2D0*k1k2*(kakb + m1**2) + k1k3*(kakb + m1**2) - 1D0*kak3*m2**2 + 
     &                      2D0*kakb*m3**2 + 2D0*m1**2*m3**2 - 1D0*kbk3*pap1 + 2D0*m1**2*pap1 - 2D0*m2**2*pap1 - 1D0*kbk3*pap2 + kak3*pbp1 + 4D0*kakb*pbp1 + 4D0*m1**2*pbp1 - 
     &                      2D0*pap2*pbp1 + kak3*pbp2 + 2D0*pap1*pbp2) + 
     &                   (2D0*k2k3*(m2**2*pap2 + kakb*(m4**2 - 2D0*pbp2) + m1**2*(m4**2 - 1D0*pap2 - 2D0*pbp2)) + 
     &                      k1k2*((kakb + m1**2)*(2D0*k2k3 + m4**2) - 2D0*kbk3*pap2 + 2D0*kak3*pbp2) + 
     &                      m4**2*(2D0*kbk3*m1**2 + m1**2*m3**2 + m1**2*pap1 - 1D0*m2**2*pap1 - 2D0*kbk3*pap2 + 2D0*m1**2*pbp1 - 1D0*pap2*pbp1 + 
     &                         kakb*(2D0*kbk3 + m3**2 + 2D0*pbp1) + pap1*pbp2 + kak3*(m1**2 - 1D0*m2**2 + 2D0*pbp2)))*xv) + 
     &                B1L*(2D0*(-2D0*k1k2*m2**2*pap1 - 2D0*k2k3*m2**2*pap1 - 1D0*m2**2*m3**2*pap1 + kbk3*mq1**2*pap1 + m2**2*m3**2*pap2 - 1D0*kak3*mq1**2*pbp1 + 
     &                      kakb*(4D0*(k1k2 + k2k3)*pbp1 + m3**2*(mq1**2 + 2D0*pbp1 - 2D0*pbp2)) + 
     &                      m1**2*(2D0*(k1k2 + k2k3)*(pap1 + 2D0*pbp1) + m3**2*(mq1**2 + pap1 - 1D0*pap2 + 2D0*pbp1 - 2D0*pbp2))) + 
     &                   (-1D0*kak3*m1**2*m3**2*m4**2 - 2D0*kbk3*m1**2*m3**2*m4**2 + kak3*m2**2*m3**2*m4**2 - 1D0*m1**2*m3**2*m4**2*mq1**2 + 
     &                      2D0*k1k2**2*(kak3*(m1 - 1D0*m2)*(m1 + m2) + (kakb + m1**2)*(2D0*kbk3 + mq1**2)) - 1D0*k2k3*m1**2*m4**2*pap1 + k2k3*m2**2*m4**2*pap1 - 
     &                      1D0*m1**2*m3**2*m4**2*pap1 + m2**2*m3**2*m4**2*pap1 + 2D0*k2k3*m1**2*m3**2*pap2 - 2D0*k2k3*m2**2*m3**2*pap2 + m1**2*m3**2*m4**2*pap2 - 
     &                      1D0*m2**2*m3**2*m4**2*pap2 - 2D0*k2k3*m1**2*m4**2*pbp1 - 2D0*m1**2*m3**2*m4**2*pbp1 - 
     &                      1D0*kakb*m4**2*(2D0*k2k3*pbp1 + m3**2*(2D0*kbk3 + mq1**2 + 2D0*pbp1)) + 2D0*(kakb + m1**2)*m3**2*(2D0*k2k3 + m4**2)*pbp2 + 
     &                      k1k2*(2D0*kakb*kbk3*m4**2 + kak3*m1**2*m4**2 + 2D0*kbk3*m1**2*m4**2 - 1D0*kak3*m2**2*m4**2 + 2D0*k2k3*kakb*mq1**2 + 2D0*k2k3*m1**2*mq1**2 - 
     &                         2D0*k2k3*m1**2*pap1 + 2D0*k2k3*m2**2*pap1 + 2D0*m1**2*m3**2*pap2 - 2D0*m2**2*m3**2*pap2 - 2D0*kbk3*mq1**2*pap2 - 4D0*k2k3*kakb*pbp1 - 
     &                         4D0*k2k3*m1**2*pbp1 - 2D0*mq1**2*pap2*pbp1 + 2D0*(2D0*(kakb + m1**2)*m3**2 + mq1**2*(kak3 + pap1))*pbp2))*xv + 
     &                   k1k3*(-2D0*(kakb + m1**2)*mq1**2 + 4D0*(m1 - 1D0*m2)*(m1 + m2)*pap1 + 8D0*(kakb + m1**2)*pbp1 - 
     &                      1D0*(2D0*k1k2 + m4**2)*((m1 - 1D0*m2)*(m1 + m2)*pap2 + 2D0*(kakb + m1**2)*pbp2)*xv)))) + 
     &          A1Lc*m1*(B1R*mq1*(B1Lc*m3*mq1c*(-2D0*(2D0*k1k2 + k1k3 + k2k3 + kak3 + kbk3 + 2D0*(m3**2 + pap1 + pbp1)) + 
     &                   (k1k2*(2D0*k2k3 + m4**2) + m4**2*(kak3 + kbk3 + m3**2 + pap1 + pbp1) + 2D0*k2k3*(m4**2 - 1D0*pap2 - 1D0*pbp2))*xv) + 
     &                B1Rc*(-2D0*k1k3*(2D0*kakb + m1**2 + m2**2 - 2D0*(pap1 + pbp1)) + 
     &                   2D0*(2D0*(k1k2 + k2k3)*(pap1 + pbp1) + m3**2*(2D0*kakb + m1**2 + m2**2 + pap1 - 1D0*pap2 + pbp1 - 1D0*pbp2)) - 
     &                   1D0*k1k3*(2D0*k1k2 + m4**2)*(pap2 + pbp2)*xv + (2D0*k1k2**2*(kak3 + 2D0*kakb + kbk3 + m1**2 + m2**2) - 1D0*kak3*m3**2*m4**2 - 2D0*kakb*m3**2*m4**2 - 
     &                      1D0*kbk3*m3**2*m4**2 - 1D0*m1**2*m3**2*m4**2 - 1D0*m2**2*m3**2*m4**2 - 1D0*k2k3*m4**2*pap1 - 1D0*m3**2*m4**2*pap1 + 2D0*k2k3*m3**2*pap2 + 
     &                      m3**2*m4**2*pap2 - 1D0*k2k3*m4**2*pbp1 - 1D0*m3**2*m4**2*pbp1 + 2D0*k2k3*m3**2*pbp2 + m3**2*m4**2*pbp2 + 
     &                      k1k2*((kak3 + kbk3)*m4**2 + 2D0*k2k3*(2D0*kakb + m1**2 + m2**2 - 1D0*pap1 - 1D0*pbp1) + 2D0*m3**2*(pap2 + pbp2)))*xv)) + 
     &             B1L*(B1Rc*m3*(-2D0*(2D0*kakb + m1**2 + m2**2)*(2D0*k1k2 + k1k3 + k2k3 + 2D0*m3**2) - 2D0*mq1**2*(kak3 + kbk3 + 2D0*(pap1 + pbp1)) + 
     &                   (k1k2*(2D0*kakb + m1**2 + m2**2)*(2D0*k2k3 + m4**2) + m4**2*((2D0*kakb + m1**2 + m2**2)*m3**2 + mq1**2*(kak3 + kbk3 + pap1 + pbp1)) + 
     &                      2D0*k2k3*((2D0*kakb + m1**2 + m2**2)*m4**2 - 1D0*mq1**2*(pap2 + pbp2)))*xv) + 
     &                B1Lc*mq1c*(4D0*(k1k2 + k2k3)*(pap1 + pbp1) - 2D0*k1k3*(mq1**2 - 2D0*(pap1 + pbp1)) + 
     &                   (2D0*k1k2**2*(kak3 + kbk3 + mq1**2) + k1k2*((kak3 + kbk3)*m4**2 + 2D0*k2k3*(mq1**2 - 1D0*pap1 - 1D0*pbp1)) - 1D0*k2k3*m4**2*(pap1 + pbp1))*xv - 
     &                   1D0*k1k3*(2D0*k1k2 + m4**2)*(pap2 + pbp2)*xv + m3**2*
     &                    (2D0*(mq1**2 + pap1 - 1D0*pap2 + pbp1 - 1D0*pbp2) - 
     &                      1D0*(m4**2*(kak3 + kbk3 + mq1**2 + pap1 - 1D0*pap2 + pbp1 - 1D0*pbp2) - 2D0*(k1k2 + k2k3)*(pap2 + pbp2))*xv)))))))/
     &   (k1k3*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2)*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark1)*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1c**2)))


        radvSfSf = dreal((-2D0*GS**2*(A1R*(A1Lc*m1*(B1R*mq1*(B1Lc*m3*(k1k3 + m3**2)*mq1c*(-4D0 + m4**2*xv) + 
     &                B1Rc*(-2D0*k1k3*(kak3 + kbk3) + k1k3*((kak3 + kbk3)*m4**2 - 2D0*k2k3*(pap2 + pbp2))*xv + 
     &                   m3**2*(2D0*(kak3 + kbk3 + pap1 + pbp1) - 1D0*(m4**2*(kak3 + kbk3 + pap1 + pbp1) - 2D0*(k1k2 + k2k3)*(pap2 + pbp2))*xv))) + 
     &             B1L*(B1Rc*(2D0*kakb + m1**2 + m2**2)*m3*(k1k3 + m3**2)*(-4D0 + m4**2*xv) + 
     &                B1Lc*mq1c*(-2D0*k1k3*(kak3 + kbk3) + k1k3*((kak3 + kbk3)*m4**2 - 2D0*k2k3*(pap2 + pbp2))*xv + 
     &                   m3**2*(2D0*(kak3 + kbk3 + pap1 + pbp1) - 1D0*(m4**2*(kak3 + kbk3 + pap1 + pbp1) - 2D0*(k1k2 + k2k3)*(pap2 + pbp2))*xv)))) + 
     &          A1Rc*(B1R*mq1*(B1Lc*(kakb + m1**2)*m3*(k1k3 + m3**2)*(-4D0 + m4**2*xv) + 
     &                B1Rc*mq1c*(k1k3*(-2D0*kak3 + kak3*m4**2*xv - 2D0*k2k3*pap2*xv) + m3**2*(2D0*(kak3 + pap1) - 1D0*(m4**2*(kak3 + pap1) - 2D0*(k1k2 + k2k3)*pap2)*xv))) + 
     &             B1L*(B1Rc*(kakb + m1**2)*m3*(k1k3 + m3**2)*mq1c*(-4D0 + m4**2*xv) + 
     &                B1Lc*(m3**2*(4D0*kbk3*m1**2 + 2D0*m1**2*pap1 - 2D0*m2**2*pap1 + 4D0*m1**2*pbp1 + 4D0*kakb*(kbk3 + pbp1) - 
     &                      2D0*kakb*(m4**2*(kbk3 + pbp1) - 2D0*(k1k2 + k2k3)*pbp2)*xv + 
     &                      (-2D0*kbk3*m1**2*m4**2 + m2**2*(m4**2*pap1 - 2D0*(k1k2 + k2k3)*pap2) + m1**2*(-1D0*m4**2*(pap1 + 2D0*pbp1) + 2D0*(k1k2 + k2k3)*(pap2 + 2D0*pbp2)))*xv - 
     &                      1D0*kak3*(m1 - 1D0*m2)*(m1 + m2)*(-2D0 + m4**2*xv)) + 
     &                   k1k3*(kak3*(m1 - 1D0*m2)*(m1 + m2)*(-2D0 + m4**2*xv) + 2D0*kakb*(-2D0*k2k3*pbp2*xv + kbk3*(-2D0 + m4**2*xv)) + 
     &                      2D0*(k2k3*(m2**2*pap2 - 1D0*m1**2*(pap2 + 2D0*pbp2))*xv + kbk3*m1**2*(-2D0 + m4**2*xv))))))) + 
     &       A1L*(A1Rc*m1*(B1Rc*mq1c*(B1L*m3*(k1k3 + m3**2)*mq1*(-4D0 + m4**2*xv) + 
     &                B1R*(-2D0*k1k3*(kak3 + kbk3) + k1k3*((kak3 + kbk3)*m4**2 - 2D0*k2k3*(pap2 + pbp2))*xv + 
     &                   m3**2*(2D0*(kak3 + kbk3 + pap1 + pbp1) - 1D0*(m4**2*(kak3 + kbk3 + pap1 + pbp1) - 2D0*(k1k2 + k2k3)*(pap2 + pbp2))*xv))) + 
     &             B1Lc*(B1R*(2D0*kakb + m1**2 + m2**2)*m3*(k1k3 + m3**2)*(-4D0 + m4**2*xv) + 
     &                B1L*mq1*(-2D0*k1k3*(kak3 + kbk3) + k1k3*((kak3 + kbk3)*m4**2 - 2D0*k2k3*(pap2 + pbp2))*xv + 
     &                   m3**2*(2D0*(kak3 + kbk3 + pap1 + pbp1) - 1D0*(m4**2*(kak3 + kbk3 + pap1 + pbp1) - 2D0*(k1k2 + k2k3)*(pap2 + pbp2))*xv)))) + 
     &          A1Lc*(B1L*mq1*(B1Rc*(kakb + m1**2)*m3*(k1k3 + m3**2)*(-4D0 + m4**2*xv) + 
     &                B1Lc*mq1c*(k1k3*(-2D0*kak3 + kak3*m4**2*xv - 2D0*k2k3*pap2*xv) + m3**2*(2D0*(kak3 + pap1) - 1D0*(m4**2*(kak3 + pap1) - 2D0*(k1k2 + k2k3)*pap2)*xv))) + 
     &             B1R*(B1Lc*(kakb + m1**2)*m3*(k1k3 + m3**2)*mq1c*(-4D0 + m4**2*xv) + 
     &                B1Rc*(m3**2*(4D0*kbk3*m1**2 + 2D0*m1**2*pap1 - 2D0*m2**2*pap1 + 4D0*m1**2*pbp1 + 4D0*kakb*(kbk3 + pbp1) - 
     &                      2D0*kakb*(m4**2*(kbk3 + pbp1) - 2D0*(k1k2 + k2k3)*pbp2)*xv + 
     &                      (-2D0*kbk3*m1**2*m4**2 + m2**2*(m4**2*pap1 - 2D0*(k1k2 + k2k3)*pap2) + m1**2*(-1D0*m4**2*(pap1 + 2D0*pbp1) + 2D0*(k1k2 + k2k3)*(pap2 + 2D0*pbp2)))*xv - 
     &                      1D0*kak3*(m1 - 1D0*m2)*(m1 + m2)*(-2D0 + m4**2*xv)) + 
     &                   k1k3*(kak3*(m1 - 1D0*m2)*(m1 + m2)*(-2D0 + m4**2*xv) + 2D0*kakb*(-2D0*k2k3*pbp2*xv + kbk3*(-2D0 + m4**2*xv)) + 
     &                      2D0*(k2k3*(m2**2*pap2 - 1D0*m1**2*(pap2 + 2D0*pbp2))*xv + kbk3*m1**2*(-2D0 + m4**2*xv)))))))))/
     &   (k1k3**2*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2)*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1c**2)))

        if(xtype.eq.4) then

        NeuQ2qx_M2RadiationVSS = CF*NC*(radvSiSi + 2.*radvSiSp + 2.*radvSiSf + radvSpSp + 2.*radvSpSf + radvSfSf - radvSiSiOS - radvSiSpOS - radvSpSpOS)

        else

        NeuQ2qx_M2RadiationVSS = CF*NC*(radvSiSi + 2.*radvSiSp + 2.*radvSiSf + radvSpSp + 2.*radvSpSf + radvSfSf)

        endif


       end

! =================================================================
        double precision function NeuQ2qx_M2RadiationVST()
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double precision radvSiTi, radvSiTp, radvSiTf, radvSiTv, radvSpTi, radvSpTp, radvSpTf, radvSpTv, radvSfTi, radvSfTp, radvSfTf, radvSfTv
        double complex A1L, A1R, A1Lc, A1Rc, A2L, A2R, A2Lc, A2Rc, B1L, B1R, B1Lc, B1Rc, C2L, C2Lc, CLL2, CLL2c

        double precision mq1, mq1c, mq2, mq2c
        double precision xv, colorfactor
        integer sumTc
        double precision wquark1, wquark2

        if (xtype.eq.1) then
        
           xv = 0d0
           colorfactor = 4d0/3d0
           
        elseif (xtype.eq.2) then
        
           xv = 0d0
           colorfactor = 1d0
           
        elseif (xtype.eq.3) then
        
           xv = 1/(m4**2)
           colorfactor = 1d0
           
        elseif (xtype.eq.4) then

           xv = 1/(m4**2)
           colorfactor = 1d0
                   
        endif

        NeuQ2qx_M2RadiationVST = 0d0

        do sumTc = 0, 1

        radvSiTi = 0d0
        radvSiTp = 0d0
        radvSiTf = 0d0
        radvSiTv = 0d0
        radvSpTi = 0d0
        radvSpTp = 0d0
        radvSpTf = 0d0
        radvSpTv = 0d0
        radvSfTi = 0d0
        radvSfTp = 0d0
        radvSfTf = 0d0
        radvSfTv = 0d0

        A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A1R = AR(chimeigen,sfmeigen,sftype,sfgen)

        A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
        A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
        
        B1L = BL(xtype,ftype)
        B1R = BR(xtype,ftype)

        B1Lc = dconjg(BL(xtype,ftype))
        B1Rc = dconjg(BR(xtype,ftype))

        mq1 = Mf(pstype,psgen) 
        mq1c = Mf(pstype,psgen) 

        A2L = AL(chimeigen,ptmeigen+sumTc,pttype,ptgen)
        A2R = AR(chimeigen,ptmeigen+sumTc,pttype,ptgen)

        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
 
        C2L = CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype)

        C2Lc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))

        CLL2 = CLL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype)

        CLL2c = dconjg(CLL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))

        mq2 = MSf(ptmeigen+sumTc,pttype,ptgen)

        mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

        wquark1 = 0d0
        wquark2 = 0d0


        radvSiTi = dreal((-4D0*C2Lc*GS**2*(kbk3 - 1D0*m2**2)*(A1R*(-1D0*A2Rc*(B1R*m3*mq1*(-2D0*m1**2 + 2D0*pap1 + pap2) + 
     &             B1L*(2D0*m1**2*m3**2 + (2D0*m3**2 + m4**2 - 4D0*pap1)*pap1 + 2D0*k1k2*(m1**2 + pap1) + (m3**2 - 4D0*pap1)*pap2)) + 
     &          A2Rc*(2D0*k1k2 + m4**2 - 2D0*pap2)*(B1L*(2D0*k1k2 + m4**2)*pap1 + m3*(-1D0*B1L*m3 + B1R*mq1)*pap2)*xv + 
     &          A2Lc*m1*(B1R*mq1*(-1D0*k1k2 - 2D0*m3**2 + 2D0*pap1 + k1k2*(2D0*k1k2 + m4**2 - 2D0*pap2)*xv) + 
     &             B1L*m3*(-3D0*k1k2 - 2D0*m3**2 - 1D0*m4**2 + 2D0*(pap1 + pap2) + (k1k2 + m4**2)*(2D0*k1k2 + m4**2 - 2D0*pap2)*xv))) + 
     &       A1L*(-1D0*A2Lc*(B1L*m3*mq1*(-2D0*m1**2 + 2D0*pap1 + pap2) + 
     &             B1R*(2D0*m1**2*m3**2 + (2D0*m3**2 + m4**2 - 4D0*pap1)*pap1 + 2D0*k1k2*(m1**2 + pap1) + (m3**2 - 4D0*pap1)*pap2)) + 
     &          A2Lc*(2D0*k1k2 + m4**2 - 2D0*pap2)*(B1R*(2D0*k1k2 + m4**2)*pap1 + m3*(-1D0*B1R*m3 + B1L*mq1)*pap2)*xv + 
     &          A2Rc*m1*(B1L*mq1*(-1D0*k1k2 - 2D0*m3**2 + 2D0*pap1 + k1k2*(2D0*k1k2 + m4**2 - 2D0*pap2)*xv) + 
     &             B1R*m3*(-3D0*k1k2 - 2D0*m3**2 - 1D0*m4**2 + 2D0*(pap1 + pap2) + (k1k2 + m4**2)*(2D0*k1k2 + m4**2 - 2D0*pap2)*xv)))))/
     &   (kbk3*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark1)*(m1**2 + m3**2 - 1D0*mq2c**2 - 2D0*pap1 + (cI*mq2c*wst(ptmeigen+sumTc)))*(2D0*k1k2 + m1**2 - 1D0*m2**2 + m3**2 + m4**2 - 2D0*(pap1 + pap2))))


        radvSiTp = dreal((C2Lc*GS**2*(k1k3 - 1D0*k2k3 - 1D0*kak3 + 2D0*kakb + kbk3 - 2D0*(m2**2 + pbp1 - 1D0*pbp2))*
     &     (A1R*(A2Lc*m1*(B1R*mq1*(-1D0*k1k2 + 2D0*pbp1 + k1k2*(m4**2 - 2D0*pbp2)*xv) + 
     &             B1L*m3*(-1D0*k1k2 - 1D0*m4**2 + 2D0*(pbp1 + pbp2) + (k1k2 + m4**2)*(m4**2 - 2D0*pbp2)*xv)) + 
     &          A2Rc*(B1R*m3*mq1*(2D0*kakb - 1D0*pap2 + pap2*(m4**2 - 2D0*pbp2)*xv) + 
     &             B1L*(-2D0*k1k2*kakb - 2D0*kakb*m3**2 - 2D0*k1k2*pap1 - 1D0*m4**2*pap1 + m3**2*pap2 + 4D0*pap1*pbp1 + 2D0*pap2*pbp1 + 2D0*pap1*pbp2 + 
     &                (2D0*k1k2*pap1 + m4**2*pap1 - 1D0*m3**2*pap2)*(m4**2 - 2D0*pbp2)*xv))) + 
     &       A1L*(A2Rc*m1*(B1L*mq1*(-1D0*k1k2 + 2D0*pbp1 + k1k2*(m4**2 - 2D0*pbp2)*xv) + 
     &             B1R*m3*(-1D0*k1k2 - 1D0*m4**2 + 2D0*(pbp1 + pbp2) + (k1k2 + m4**2)*(m4**2 - 2D0*pbp2)*xv)) + 
     &          A2Lc*(B1L*m3*mq1*(2D0*kakb - 1D0*pap2 + pap2*(m4**2 - 2D0*pbp2)*xv) + 
     &             B1R*(-2D0*k1k2*kakb - 2D0*kakb*m3**2 - 2D0*k1k2*pap1 - 1D0*m4**2*pap1 + m3**2*pap2 + 4D0*pap1*pbp1 + 2D0*pap2*pbp1 + 2D0*pap1*pbp2 + 
     &                (2D0*k1k2*pap1 + m4**2*pap1 - 1D0*m3**2*pap2)*(m4**2 - 2D0*pbp2)*xv)))))/
     &   (kbk3*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark2)*(m1**2 + m3**2 - 1D0*mq2c**2 - 2D0*pap1 + (cI*mq2c*wst(ptmeigen+sumTc)))*(m2**2 + m4**2 - 1D0*mq2c**2 - 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc)))))


        radvSiTf = dreal((C2Lc*GS**2*(A1R*(A2Lc*m1*(B1L*m3*(-1D0*k1k2*k1k3 + k1k2*kbk3 - 2D0*k1k3*m2**2 - 2D0*k2k3*m2**2 - 1D0*k1k3*m4**2 + kbk3*m4**2 + 2D0*k1k2*pbp1 + 2D0*k1k3*pbp1 - 
     &                1D0*k2k3*pbp1 + 2D0*m4**2*pbp1 - 4D0*pbp1**2 + 3D0*k1k3*pbp2 - 4D0*pbp1*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*(k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + m4**2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + k2k3*pbp1 - 1D0*k1k3*pbp2)*xv) + 
     &             B1R*mq1*(-1D0*k1k2*k1k3 + k1k2*kbk3 + 2D0*k1k3*m2**2 + 2D0*k1k2*pbp1 + 2D0*k1k3*pbp1 + k2k3*pbp1 - 4D0*kbk3*pbp1 - 4D0*pbp1**2 - 1D0*k1k3*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*(k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*k2k3*pbp1 + k1k3*pbp2)*xv)) + 
     &          A2Rc*(B1R*m3*mq1*(2D0*k1k3*kakb - 1D0*k2k3*kakb - 2D0*kak3*m2**2 - 1D0*k1k3*pap2 + kbk3*pap2 - 4D0*kakb*pbp1 + 2D0*pap2*pbp1 + kak3*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*(k2k3*kakb + pap2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*kak3*pbp2)*xv) + 
     &             B1L*(k2k3*kakb*m3**2 + 2D0*kak3*m2**2*m3**2 - 2D0*k2k3*m2**2*pap1 + kbk3*m4**2*pap1 - 1D0*kbk3*m3**2*pap2 + 4D0*k2k3*kakb*pbp1 + 4D0*kakb*m3**2*pbp1 + 
     &                kak3*m4**2*pbp1 - 4D0*kbk3*pap1*pbp1 + 2D0*m4**2*pap1*pbp1 - 4D0*kbk3*pap2*pbp1 - 2D0*m3**2*pap2*pbp1 - 4D0*kak3*pbp1**2 - 8D0*pap1*pbp1**2 - 
     &                4D0*pap2*pbp1**2 - 1D0*kak3*m3**2*pbp2 - 4D0*kak3*pbp1*pbp2 - 4D0*pap1*pbp1*pbp2 - 
     &                1D0*(m4**2 - 2D0*pbp2)*(m4**2*(kbk3*pap1 + kak3*pbp1 + 2D0*pap1*pbp1) + m3**2*(k2k3*kakb - 1D0*kbk3*pap2 - 2D0*pap2*pbp1 - 1D0*kak3*pbp2))*xv + 
     &                k1k3*(-2D0*kakb*m3**2 - 1D0*kakb*m4**2 - 1D0*m4**2*pap1 + 2D0*m2**2*pap2 + m3**2*pap2 + 4D0*kakb*pbp1 + 4D0*pap1*pbp1 + 2D0*pap2*pbp1 + 2D0*pap1*pbp2 + 
     &                   (m4**2*(kakb + pap1) - 1D0*m3**2*pap2)*(m4**2 - 2D0*pbp2)*xv) + 
     &                2D0*k1k2*(-2D0*k1k3*kakb - 1D0*k1k3*pap1 + kbk3*pap1 + 2D0*kakb*pbp1 + 2D0*pap1*pbp1 + kak3*(m2**2 + pbp1) + 
     &                   (k1k3*(kakb + pap1) - 1D0*kak3*pbp1 - 1D0*pap1*(kbk3 + 2D0*pbp1))*(m4**2 - 2D0*pbp2)*xv)))) + 
     &       A1L*(A2Rc*m1*(B1R*m3*(-1D0*k1k2*k1k3 + k1k2*kbk3 - 2D0*k1k3*m2**2 - 2D0*k2k3*m2**2 - 1D0*k1k3*m4**2 + kbk3*m4**2 + 2D0*k1k2*pbp1 + 2D0*k1k3*pbp1 - 1D0*k2k3*pbp1 + 
     &                2D0*m4**2*pbp1 - 4D0*pbp1**2 + 3D0*k1k3*pbp2 - 4D0*pbp1*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*(k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + m4**2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + k2k3*pbp1 - 1D0*k1k3*pbp2)*xv) + 
     &             B1L*mq1*(-1D0*k1k2*k1k3 + k1k2*kbk3 + 2D0*k1k3*m2**2 + 2D0*k1k2*pbp1 + 2D0*k1k3*pbp1 + k2k3*pbp1 - 4D0*kbk3*pbp1 - 4D0*pbp1**2 - 1D0*k1k3*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*(k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*k2k3*pbp1 + k1k3*pbp2)*xv)) + 
     &          A2Lc*(B1L*m3*mq1*(2D0*k1k3*kakb - 1D0*k2k3*kakb - 2D0*kak3*m2**2 - 1D0*k1k3*pap2 + kbk3*pap2 - 4D0*kakb*pbp1 + 2D0*pap2*pbp1 + kak3*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*(k2k3*kakb + pap2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*kak3*pbp2)*xv) + 
     &             B1R*(k2k3*kakb*m3**2 + 2D0*kak3*m2**2*m3**2 - 2D0*k2k3*m2**2*pap1 + kbk3*m4**2*pap1 - 1D0*kbk3*m3**2*pap2 + 4D0*k2k3*kakb*pbp1 + 4D0*kakb*m3**2*pbp1 + 
     &                kak3*m4**2*pbp1 - 4D0*kbk3*pap1*pbp1 + 2D0*m4**2*pap1*pbp1 - 4D0*kbk3*pap2*pbp1 - 2D0*m3**2*pap2*pbp1 - 4D0*kak3*pbp1**2 - 8D0*pap1*pbp1**2 - 
     &                4D0*pap2*pbp1**2 - 1D0*kak3*m3**2*pbp2 - 4D0*kak3*pbp1*pbp2 - 4D0*pap1*pbp1*pbp2 - 
     &                1D0*(m4**2 - 2D0*pbp2)*(m4**2*(kbk3*pap1 + kak3*pbp1 + 2D0*pap1*pbp1) + m3**2*(k2k3*kakb - 1D0*kbk3*pap2 - 2D0*pap2*pbp1 - 1D0*kak3*pbp2))*xv + 
     &                k1k3*(-2D0*kakb*m3**2 - 1D0*kakb*m4**2 - 1D0*m4**2*pap1 + 2D0*m2**2*pap2 + m3**2*pap2 + 4D0*kakb*pbp1 + 4D0*pap1*pbp1 + 2D0*pap2*pbp1 + 2D0*pap1*pbp2 + 
     &                   (m4**2*(kakb + pap1) - 1D0*m3**2*pap2)*(m4**2 - 2D0*pbp2)*xv) + 
     &                2D0*k1k2*(-2D0*k1k3*kakb - 1D0*k1k3*pap1 + kbk3*pap1 + 2D0*kakb*pbp1 + 2D0*pap1*pbp1 + kak3*(m2**2 + pbp1) + 
     &                   (k1k3*(kakb + pap1) - 1D0*kak3*pbp1 - 1D0*pap1*(kbk3 + 2D0*pbp1))*(m4**2 - 2D0*pbp2)*xv))))))/
     &   (k1k3*kbk3*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark2)*(m2**2 + m4**2 - 1D0*mq2c**2 - 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc)))))


        radvSiTv = dreal((CLL2c*GS*(A1R*(A2Lc*m1*(B1R*mq1*(-1D0*k1k3 + 2D0*pbp1 + k1k2*(k2k3 - 2D0*pbp2)*xv) + 
     &             B1L*m3*(-1D0*k1k3 - 1D0*k2k3 + 2D0*(pbp1 + pbp2) + (k1k2 + m4**2)*(k2k3 - 2D0*pbp2)*xv)) + 
     &          A2Rc*(B1R*m3*mq1*(-1D0*kak3 + 2D0*kakb + pap2*(k2k3 - 2D0*pbp2)*xv) + 
     &             B1L*((kak3 - 2D0*kakb)*(k1k2 + m3**2) - 1D0*k2k3*pap1 - 1D0*k1k3*(2D0*pap1 + pap2) + 4D0*pap1*pbp1 + 2D0*pap2*pbp1 + 2D0*pap1*pbp2 + 
     &                (2D0*k1k2*pap1 + m4**2*pap1 - 1D0*m3**2*pap2)*(k2k3 - 2D0*pbp2)*xv))) + 
     &       A1L*(A2Rc*m1*(B1L*mq1*(-1D0*k1k3 + 2D0*pbp1 + k1k2*(k2k3 - 2D0*pbp2)*xv) + B1R*m3*(-1D0*k1k3 - 1D0*k2k3 + 2D0*(pbp1 + pbp2) + (k1k2 + m4**2)*(k2k3 - 2D0*pbp2)*xv)) + 
     &          A2Lc*(B1L*m3*mq1*(-1D0*kak3 + 2D0*kakb + pap2*(k2k3 - 2D0*pbp2)*xv) + 
     &             B1R*((kak3 - 2D0*kakb)*(k1k2 + m3**2) - 1D0*k2k3*pap1 - 1D0*k1k3*(2D0*pap1 + pap2) + 4D0*pap1*pbp1 + 2D0*pap2*pbp1 + 2D0*pap1*pbp2 + 
     &                (2D0*k1k2*pap1 + m4**2*pap1 - 1D0*m3**2*pap2)*(k2k3 - 2D0*pbp2)*xv)))))/
     &   (kbk3*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark1)*(-2D0*k2k3 + 2D0*kbk3 - 1D0*m2**2 - 1D0*m4**2 + mq2c**2 + 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc)))))


        radvSpTi = dreal((2D0*C2Lc*GS**2*(-1D0*A1R*(A2Rc*(B1R*m3*mq1*(3D0*k1k2*kak3 - 6D0*k1k2*kakb + 2D0*k1k3*kakb + k2k3*kakb - 1D0*k2k3*m1**2 - 2D0*kak3*m1**2 + 4D0*kakb*m1**2 - 
     &                2D0*kbk3*m1**2 + 4D0*m1**2*m2**2 + 2D0*kak3*m3**2 - 4D0*kakb*m3**2 + kak3*m4**2 - 2D0*kakb*m4**2 + k2k3*pap1 + 2D0*kbk3*pap1 - 4D0*m2**2*pap1 - 
     &                1D0*k1k3*pap2 + kbk3*pap2 - 2D0*m2**2*pap2 - 2D0*kak3*pbp1 + 2D0*pap2*pbp1 - 1D0*(kak3 - 2D0*m1**2 + 2D0*pap1)*pbp2) + 
     &             B1L*(k2k3*kakb*m3**2 + k2k3*m1**2*m3**2 + 2D0*kak3*m1**2*m3**2 - 4D0*kakb*m1**2*m3**2 + 2D0*kbk3*m1**2*m3**2 - 4D0*m1**2*m2**2*m3**2 + 2D0*kak3*m3**2*mq1**2 - 
     &                4D0*kakb*m3**2*mq1**2 - 4D0*k2k3*kakb*pap1 - 2D0*k2k3*m1**2*pap1 + 2D0*kbk3*m3**2*pap1 - 4D0*m2**2*m3**2*pap1 + kbk3*m4**2*pap1 - 2D0*m2**2*m4**2*pap1 + 
     &                k2k3*mq1**2*pap1 - 4D0*kak3*mq1**2*pap1 + 8D0*kakb*mq1**2*pap1 - 4D0*kbk3*pap1**2 + 8D0*m2**2*pap1**2 + kbk3*m3**2*pap2 - 2D0*m2**2*m3**2*pap2 - 
     &                4D0*kbk3*pap1*pap2 + 8D0*m2**2*pap1*pap2 + 2D0*k2k3*m1**2*pbp1 - 2D0*kak3*m3**2*pbp1 - 4D0*m1**2*m3**2*pbp1 - 1D0*kak3*m4**2*pbp1 - 2D0*m1**2*m4**2*pbp1 - 
     &                4D0*m1**2*mq1**2*pbp1 + 4D0*kak3*pap1*pbp1 + 8D0*m1**2*pap1*pbp1 + 4D0*m1**2*pap2*pbp1 + 2D0*mq1**2*pap2*pbp1 + 
     &                k1k2*(2D0*k1k3*(kakb + m1**2) + kak3*(2D0*m1**2 + mq1**2 - 2D0*pbp1) - 
     &                   2D0*(kakb*(2D0*m1**2 + mq1**2) - 1D0*(kbk3 - 2D0*m2**2)*(m1**2 + pap1) + 2D0*m1**2*pbp1)) + 
     &                k1k3*(kakb*(2D0*m3**2 + m4**2 - 4D0*pap1) - 1D0*mq1**2*pap2 + m1**2*(2D0*m3**2 + m4**2 + 2D0*mq1**2 - 4D0*pap1 - 2D0*pap2 - 2D0*pbp2)) - 
     &                1D0*((kak3 + 2D0*m1**2)*m3**2 + 2D0*(-2D0*(kak3 + m1**2) + mq1**2)*pap1)*pbp2)) + 
     &          A2Lc*m1*(B1R*mq1*(k2k3*m3**2 + 4D0*kak3*m3**2 - 8D0*kakb*m3**2 + 2D0*kbk3*m3**2 - 4D0*m2**2*m3**2 - 1D0*k2k3*pap1 - 4D0*kak3*pap1 + 8D0*kakb*pap1 - 2D0*kbk3*pap1 + 
     &                4D0*m2**2*pap1 + k1k2*(2D0*k1k3 + 3D0*kak3 - 6D0*kakb + kbk3 - 2D0*m2**2 - 4D0*pbp1) + k2k3*pbp1 - 2D0*kak3*pbp1 - 4D0*m1**2*pbp1 - 4D0*m3**2*pbp1 - 
     &                2D0*m4**2*pbp1 + 8D0*pap1*pbp1 + 6D0*pap2*pbp1 + k1k3*(2D0*kakb + 2D0*m1**2 + 2D0*m3**2 + m4**2 - 4D0*pap1 - 3D0*pap2 - 1D0*pbp2) - 2D0*m3**2*pbp2 + 
     &                2D0*pap1*pbp2) + B1L*m3*(-2D0*k2k3*kakb - 2D0*k2k3*m1**2 + 3D0*k1k2*(kak3 - 2D0*kakb + kbk3 - 2D0*m2**2) + 2D0*kak3*m3**2 - 4D0*kakb*m3**2 + 2D0*kbk3*m3**2 - 
     &                4D0*m2**2*m3**2 + kak3*m4**2 - 2D0*kakb*m4**2 + kbk3*m4**2 - 2D0*m2**2*m4**2 + k2k3*mq1**2 - 2D0*kak3*mq1**2 + 4D0*kakb*mq1**2 + k2k3*pap1 - 2D0*kbk3*pap1 + 
     &                4D0*m2**2*pap1 - 2D0*kbk3*pap2 + 4D0*m2**2*pap2 + k2k3*pbp1 + 2D0*kak3*pbp1 + 4D0*m1**2*pbp1 - 4D0*mq1**2*pbp1 + 2D0*pap2*pbp1 + 
     &                2D0*(kak3 + 2D0*m1**2 - 1D0*mq1**2 - 1D0*pap1)*pbp2 - 1D0*k1k3*(2D0*kakb + 2D0*m1**2 - 2D0*mq1**2 + pap2 + pbp2)))) - 
     &       1D0*A1L*(A2Lc*(B1L*m3*mq1*(3D0*k1k2*kak3 - 6D0*k1k2*kakb + 2D0*k1k3*kakb + k2k3*kakb - 1D0*k2k3*m1**2 - 2D0*kak3*m1**2 + 4D0*kakb*m1**2 - 2D0*kbk3*m1**2 + 
     &                4D0*m1**2*m2**2 + 2D0*kak3*m3**2 - 4D0*kakb*m3**2 + kak3*m4**2 - 2D0*kakb*m4**2 + k2k3*pap1 + 2D0*kbk3*pap1 - 4D0*m2**2*pap1 - 1D0*k1k3*pap2 + kbk3*pap2 - 
     &                2D0*m2**2*pap2 - 2D0*kak3*pbp1 + 2D0*pap2*pbp1 - 1D0*(kak3 - 2D0*m1**2 + 2D0*pap1)*pbp2) + 
     &             B1R*(k2k3*kakb*m3**2 + k2k3*m1**2*m3**2 + 2D0*kak3*m1**2*m3**2 - 4D0*kakb*m1**2*m3**2 + 2D0*kbk3*m1**2*m3**2 - 4D0*m1**2*m2**2*m3**2 + 2D0*kak3*m3**2*mq1**2 - 
     &                4D0*kakb*m3**2*mq1**2 - 4D0*k2k3*kakb*pap1 - 2D0*k2k3*m1**2*pap1 + 2D0*kbk3*m3**2*pap1 - 4D0*m2**2*m3**2*pap1 + kbk3*m4**2*pap1 - 2D0*m2**2*m4**2*pap1 + 
     &                k2k3*mq1**2*pap1 - 4D0*kak3*mq1**2*pap1 + 8D0*kakb*mq1**2*pap1 - 4D0*kbk3*pap1**2 + 8D0*m2**2*pap1**2 + kbk3*m3**2*pap2 - 2D0*m2**2*m3**2*pap2 - 
     &                4D0*kbk3*pap1*pap2 + 8D0*m2**2*pap1*pap2 + 2D0*k2k3*m1**2*pbp1 - 2D0*kak3*m3**2*pbp1 - 4D0*m1**2*m3**2*pbp1 - 1D0*kak3*m4**2*pbp1 - 2D0*m1**2*m4**2*pbp1 - 
     &                4D0*m1**2*mq1**2*pbp1 + 4D0*kak3*pap1*pbp1 + 8D0*m1**2*pap1*pbp1 + 4D0*m1**2*pap2*pbp1 + 2D0*mq1**2*pap2*pbp1 + 
     &                k1k2*(2D0*k1k3*(kakb + m1**2) + kak3*(2D0*m1**2 + mq1**2 - 2D0*pbp1) - 
     &                   2D0*(kakb*(2D0*m1**2 + mq1**2) - 1D0*(kbk3 - 2D0*m2**2)*(m1**2 + pap1) + 2D0*m1**2*pbp1)) + 
     &                k1k3*(kakb*(2D0*m3**2 + m4**2 - 4D0*pap1) - 1D0*mq1**2*pap2 + m1**2*(2D0*m3**2 + m4**2 + 2D0*mq1**2 - 4D0*pap1 - 2D0*pap2 - 2D0*pbp2)) - 
     &                1D0*((kak3 + 2D0*m1**2)*m3**2 + 2D0*(-2D0*(kak3 + m1**2) + mq1**2)*pap1)*pbp2)) + 
     &          A2Rc*m1*(B1L*mq1*(k2k3*m3**2 + 4D0*kak3*m3**2 - 8D0*kakb*m3**2 + 2D0*kbk3*m3**2 - 4D0*m2**2*m3**2 - 1D0*k2k3*pap1 - 4D0*kak3*pap1 + 8D0*kakb*pap1 - 2D0*kbk3*pap1 + 
     &                4D0*m2**2*pap1 + k1k2*(2D0*k1k3 + 3D0*kak3 - 6D0*kakb + kbk3 - 2D0*m2**2 - 4D0*pbp1) + k2k3*pbp1 - 2D0*kak3*pbp1 - 4D0*m1**2*pbp1 - 4D0*m3**2*pbp1 - 
     &                2D0*m4**2*pbp1 + 8D0*pap1*pbp1 + 6D0*pap2*pbp1 + k1k3*(2D0*kakb + 2D0*m1**2 + 2D0*m3**2 + m4**2 - 4D0*pap1 - 3D0*pap2 - 1D0*pbp2) - 2D0*m3**2*pbp2 + 
     &                2D0*pap1*pbp2) + B1R*m3*(-2D0*k2k3*kakb - 2D0*k2k3*m1**2 + 3D0*k1k2*(kak3 - 2D0*kakb + kbk3 - 2D0*m2**2) + 2D0*kak3*m3**2 - 4D0*kakb*m3**2 + 2D0*kbk3*m3**2 - 
     &                4D0*m2**2*m3**2 + kak3*m4**2 - 2D0*kakb*m4**2 + kbk3*m4**2 - 2D0*m2**2*m4**2 + k2k3*mq1**2 - 2D0*kak3*mq1**2 + 4D0*kakb*mq1**2 + k2k3*pap1 - 2D0*kbk3*pap1 + 
     &                4D0*m2**2*pap1 - 2D0*kbk3*pap2 + 4D0*m2**2*pap2 + k2k3*pbp1 + 2D0*kak3*pbp1 + 4D0*m1**2*pbp1 - 4D0*mq1**2*pbp1 + 2D0*pap2*pbp1 + 
     &                2D0*(kak3 + 2D0*m1**2 - 1D0*mq1**2 - 1D0*pap1)*pbp2 - 1D0*k1k3*(2D0*kakb + 2D0*m1**2 - 2D0*mq1**2 + pap2 + pbp2)))) - 
     &       1D0*(2D0*k1k2 + m4**2 - 2D0*pap2)*(-1D0*A1R*(A2Rc*(B1R*m3*mq1*
     &                 (k1k2*(kak3 - 2D0*kakb) + (kak3 - 2D0*kakb)*m4**2 + k2k3*(kakb + m1**2 - 1D0*pap1) + pap2*(k1k3 + kbk3 - 2D0*(m2**2 + pbp1)) - 
     &                   1D0*(kak3 + 2D0*m1**2 - 2D0*pap1)*pbp2) + B1L*(-1D0*k2k3*(kakb + m1**2)*m3**2 + k1k3*kakb*m4**2 + k1k3*m1**2*m4**2 + kbk3*m4**2*pap1 - 
     &                   2D0*m2**2*m4**2*pap1 + k2k3*mq1**2*pap1 - 1D0*kbk3*m3**2*pap2 + 2D0*m2**2*m3**2*pap2 - 1D0*k1k3*mq1**2*pap2 - 1D0*kak3*m4**2*pbp1 - 2D0*m1**2*m4**2*pbp1 + 
     &                   2D0*mq1**2*pap2*pbp1 + k1k2*(2D0*k1k3*(kakb + m1**2) + (kak3 - 2D0*kakb)*mq1**2 + 2D0*(kbk3 - 2D0*m2**2)*pap1 - 2D0*(kak3 + 2D0*m1**2)*pbp1) + 
     &                   kak3*m3**2*pbp2 + 2D0*m1**2*m3**2*pbp2 - 2D0*mq1**2*pap1*pbp2)) + 
     &             A2Lc*m1*(B1R*mq1*(k1k2*(2D0*k1k3 + kak3 - 2D0*kakb + kbk3 - 2D0*m2**2 - 4D0*pbp1) + (m4**2 - 1D0*pap2)*(k1k3 - 2D0*pbp1) + k2k3*(-1D0*m3**2 + pap1 + pbp1) - 
     &                   1D0*(k1k3 - 2D0*m3**2 + 2D0*pap1)*pbp2) + B1L*m3*
     &                 ((kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*(k1k2 + m4**2) + pap2*(k1k3 - 2D0*pbp1) + k2k3*(mq1**2 - 1D0*pap1 - 1D0*pbp1) + (k1k3 - 2D0*mq1**2 + 2D0*pap1)*pbp2)))
     &           - 1D0*A1L*(A2Lc*(B1L*m3*mq1*(k1k2*(kak3 - 2D0*kakb) + (kak3 - 2D0*kakb)*m4**2 + k2k3*(kakb + m1**2 - 1D0*pap1) + pap2*(k1k3 + kbk3 - 2D0*(m2**2 + pbp1)) - 
     &                   1D0*(kak3 + 2D0*m1**2 - 2D0*pap1)*pbp2) + B1R*(-1D0*k2k3*(kakb + m1**2)*m3**2 + k1k3*kakb*m4**2 + k1k3*m1**2*m4**2 + kbk3*m4**2*pap1 - 
     &                   2D0*m2**2*m4**2*pap1 + k2k3*mq1**2*pap1 - 1D0*kbk3*m3**2*pap2 + 2D0*m2**2*m3**2*pap2 - 1D0*k1k3*mq1**2*pap2 - 1D0*kak3*m4**2*pbp1 - 2D0*m1**2*m4**2*pbp1 + 
     &                   2D0*mq1**2*pap2*pbp1 + k1k2*(2D0*k1k3*(kakb + m1**2) + (kak3 - 2D0*kakb)*mq1**2 + 2D0*(kbk3 - 2D0*m2**2)*pap1 - 2D0*(kak3 + 2D0*m1**2)*pbp1) + 
     &                   kak3*m3**2*pbp2 + 2D0*m1**2*m3**2*pbp2 - 2D0*mq1**2*pap1*pbp2)) + 
     &             A2Rc*m1*(B1L*mq1*(k1k2*(2D0*k1k3 + kak3 - 2D0*kakb + kbk3 - 2D0*m2**2 - 4D0*pbp1) + (m4**2 - 1D0*pap2)*(k1k3 - 2D0*pbp1) + k2k3*(-1D0*m3**2 + pap1 + pbp1) - 
     &                   1D0*(k1k3 - 2D0*m3**2 + 2D0*pap1)*pbp2) + B1R*m3*
     &                 ((kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*(k1k2 + m4**2) + pap2*(k1k3 - 2D0*pbp1) + k2k3*(mq1**2 - 1D0*pap1 - 1D0*pbp1) + (k1k3 - 2D0*mq1**2 + 2D0*pap1)*pbp2))))
     &         *xv))/((2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2 - cI*mq1*wquark1)*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark1)*(m1**2 + m3**2 - 1D0*mq2c**2 - 2D0*pap1 + (cI*mq2c*wst(ptmeigen+sumTc)))*
     &     (2D0*k1k2 + m1**2 - 1D0*m2**2 + m3**2 + m4**2 - 2D0*(pap1 + pap2))))


        radvSpTp = dreal((2D0*C2Lc*GS**2*(A1R*(A2Rc*(B1R*m3*mq1*(-2D0*k1k2*kakb - 4D0*kakb**2 - 2D0*kakb*m1**2 + 2D0*kakb*m2**2 + 4D0*m1**2*m2**2 + 2D0*kakb*m3**2 - 2D0*kakb*m4**2 + 
     &                2D0*m1**2*m4**2 - 4D0*m2**2*pap1 - 2D0*m4**2*pap1 + 2D0*k1k2*pap2 + 2D0*kakb*pap2 + m1**2*pap2 - 1D0*m2**2*pap2 - 1D0*m3**2*pap2 + m4**2*pap2 + 
     &                8D0*kakb*pbp1 - 6D0*pap2*pbp1 - 6D0*m1**2*pbp2 + 6D0*pap1*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*(2D0*m4**2*pap1 + 2D0*k1k2*(kakb - 1D0*pap2) + (-2D0*kakb + m2**2 + m3**2)*pap2 - 1D0*m4**2*pap2 + 2D0*pap2*pbp1 - 
     &                   1D0*m1**2*(2D0*m4**2 + pap2 - 2D0*pbp2) - 2D0*pap1*pbp2)*xv) + 
     &             B1L*(4D0*kakb**2*m3**2 - 4D0*m1**2*m2**2*m3**2 - 2D0*m1**2*m3**2*m4**2 + 2D0*m2**2*m3**2*pap1 + m1**2*m4**2*pap1 - 3D0*m2**2*m4**2*pap1 + 
     &                2D0*m2**2*mq1**2*pap1 + m4**2*mq1**2*pap1 - 1D0*m1**2*m3**2*pap2 - 3D0*m2**2*m3**2*pap2 + m3**2*mq1**2*pap2 + 2D0*m1**2*m3**2*pbp1 - 4D0*m1**2*m4**2*pbp1 - 
     &                2D0*m1**2*mq1**2*pbp1 - 4D0*m1**2*pap1*pbp1 + 4D0*m2**2*pap1*pbp1 + 4D0*mq1**2*pap1*pbp1 - 2D0*m1**2*pap2*pbp1 + 2D0*m2**2*pap2*pbp1 + m3**2*pap2*pbp1 - 
     &                1D0*m4**2*pap2*pbp1 - 1D0*mq1**2*pap2*pbp1 + 8D0*m1**2*pbp1**2 + 8D0*pap2*pbp1**2 + 6D0*m1**2*m3**2*pbp2 - 2D0*m1**2*pap1*pbp2 + 2D0*m2**2*pap1*pbp2 - 
     &                1D0*m3**2*pap1*pbp2 + m4**2*pap1*pbp2 - 3D0*mq1**2*pap1*pbp2 + 8D0*m1**2*pbp1*pbp2 - 8D0*pap1*pbp1*pbp2 - 
     &                1D0*(m4**2 - 2D0*pbp2)*(m4**2*mq1**2*pap1 + m3**2*mq1**2*pap2 + m2**2*(-1D0*m4**2*pap1 + m3**2*pap2) + m3**2*pap2*pbp1 - 1D0*m4**2*pap2*pbp1 + 
     &                   mq1**2*pap2*pbp1 + m1**2*(m4**2*(pap1 - 2D0*pbp1) - 1D0*m3**2*(2D0*m4**2 + pap2 - 2D0*pbp2)) - 1D0*(m3**2 - 1D0*m4**2 + mq1**2)*pap1*pbp2)*xv + 
     &                2D0*kakb*(m1**2*m3**2 - 1D0*m2**2*m3**2 - 1D0*m3**2*m4**2 - 1D0*m3**2*mq1**2 + m4**2*pap1 - 1D0*m3**2*pap2 - 4D0*pap1*pbp1 - 2D0*pap2*pbp1 + 4D0*m3**2*pbp2 - 
     &                   2D0*pap1*pbp2 + (-1D0*m4**2*pap1 + m3**2*(m4**2 + pap2))*(m4**2 - 2D0*pbp2)*xv) - 2D0*k1k2**2*(kakb + m1**2)*(-1D0 + m4**2*xv - 2D0*pbp2*xv) + 
     &                k1k2*(4D0*kakb**2 + kakb*(2D0*m1**2 - 2D0*m2**2 - 1D0*m3**2 + m4**2 + mq1**2 + 4D0*pap1 - 8D0*pbp1 + 
     &                      (m3**2 - 1D0*m4**2 + mq1**2 - 4D0*pap1)*(m4**2 - 2D0*pbp2)*xv) + 
     &                   m1**2*(-4D0*m2**2 - 1D0*m3**2 + m4**2 + mq1**2 + 2D0*pap1 - 8D0*pbp1 + (m3**2 - 1D0*m4**2 - 1D0*mq1**2 - 2D0*pap1 + 4D0*pbp1)*(m4**2 - 2D0*pbp2)*xv) + 
     &                   2D0*(-1D0*pap2*pbp1 + pap1*(m2**2 - 1D0*mq1**2 + pbp2) + (pap2*pbp1 + pap1*(m2**2 + mq1**2 - 1D0*pbp2))*(m4**2 - 2D0*pbp2)*xv)))) + 
     &          A2Lc*m1*(B1R*mq1*(2D0*k1k2**2 + 4D0*k1k2*kakb + k1k2*m1**2 - 1D0*k1k2*m2**2 - 1D0*k1k2*m3**2 - 4D0*m2**2*m3**2 + k1k2*m4**2 - 2D0*m3**2*m4**2 + 4D0*m2**2*pap1 + 
     &                2D0*m4**2*pap1 - 8D0*k1k2*pbp1 - 4D0*kakb*pbp1 - 2D0*m1**2*pbp1 + 2D0*m2**2*pbp1 + 2D0*m3**2*pbp1 - 2D0*m4**2*pbp1 - 2D0*pap2*pbp1 + 8D0*pbp1**2 + 
     &                6D0*m3**2*pbp2 - 6D0*pap1*pbp2 - 1D0*(2D0*k1k2**2 + k1k2*(m1**2 - 1D0*m2**2 - 1D0*m3**2 + m4**2 - 4D0*pbp1) + 
     &                   2D0*(pap2*pbp1 - 1D0*(m3**2 - 1D0*pap1)*(m4**2 - 1D0*pbp2)))*(m4**2 - 2D0*pbp2)*xv) + 
     &             B1L*m3*(2D0*kakb*(m3 - 1D0*m4)*(m3 + m4) + m4**2*mq1**2 - 2D0*m4**2*pap1 + m4**2*pap2 + mq1**2*pap2 - 2D0*m4**2*pbp1 + 2D0*mq1**2*pbp1 - 6D0*pap2*pbp1 + 
     &                m4**2*pbp2 - 3D0*mq1**2*pbp2 + 6D0*pap1*pbp2 - 1D0*m3**2*(pap2 + pbp2) + m1**2*(m4**2 - 2D0*(pbp1 + pbp2)) - 2D0*kakb*(mq1**2 - 2D0*(pbp1 + pbp2)) - 
     &                1D0*(m4**2 - 2D0*pbp2)*(m1**2*m4**2 + mq1**2*pap2 - 2D0*pap2*pbp1 - 1D0*mq1**2*pbp2 + 2D0*pap1*pbp2 - 1D0*m3**2*(pap2 + pbp2) + 
     &                   m4**2*(mq1**2 - 2D0*pap1 + pap2 - 2D0*pbp1 + pbp2))*xv + k1k2*(m1**2 - 1D0*m2**2 - 1D0*mq1**2 + 2D0*(pap2 + pbp2))*(1D0 - 1D0*m4**2*xv + 2D0*pbp2*xv) + 
     &                m2**2*(2D0*m3**2 + 2D0*(mq1**2 - 2D0*pap1 - 2D0*pap2 + pbp1 + pbp2) + m4**4*xv - 1D0*m4**2*(3D0 + 2D0*pbp2*xv))))) + 
     &       A1L*(A2Lc*(B1L*m3*mq1*(-2D0*k1k2*kakb - 4D0*kakb**2 - 2D0*kakb*m1**2 + 2D0*kakb*m2**2 + 4D0*m1**2*m2**2 + 2D0*kakb*m3**2 - 2D0*kakb*m4**2 + 2D0*m1**2*m4**2 - 
     &                4D0*m2**2*pap1 - 2D0*m4**2*pap1 + 2D0*k1k2*pap2 + 2D0*kakb*pap2 + m1**2*pap2 - 1D0*m2**2*pap2 - 1D0*m3**2*pap2 + m4**2*pap2 + 8D0*kakb*pbp1 - 6D0*pap2*pbp1 - 
     &                6D0*m1**2*pbp2 + 6D0*pap1*pbp2 + (m4**2 - 2D0*pbp2)*
     &                 (2D0*m4**2*pap1 + 2D0*k1k2*(kakb - 1D0*pap2) + (-2D0*kakb + m2**2 + m3**2)*pap2 - 1D0*m4**2*pap2 + 2D0*pap2*pbp1 - 1D0*m1**2*(2D0*m4**2 + pap2 - 2D0*pbp2) - 
     &                   2D0*pap1*pbp2)*xv) + B1R*(4D0*kakb**2*m3**2 - 4D0*m1**2*m2**2*m3**2 - 2D0*m1**2*m3**2*m4**2 + 2D0*m2**2*m3**2*pap1 + m1**2*m4**2*pap1 - 
     &                3D0*m2**2*m4**2*pap1 + 2D0*m2**2*mq1**2*pap1 + m4**2*mq1**2*pap1 - 1D0*m1**2*m3**2*pap2 - 3D0*m2**2*m3**2*pap2 + m3**2*mq1**2*pap2 + 2D0*m1**2*m3**2*pbp1 - 
     &                4D0*m1**2*m4**2*pbp1 - 2D0*m1**2*mq1**2*pbp1 - 4D0*m1**2*pap1*pbp1 + 4D0*m2**2*pap1*pbp1 + 4D0*mq1**2*pap1*pbp1 - 2D0*m1**2*pap2*pbp1 + 2D0*m2**2*pap2*pbp1 + 
     &                m3**2*pap2*pbp1 - 1D0*m4**2*pap2*pbp1 - 1D0*mq1**2*pap2*pbp1 + 8D0*m1**2*pbp1**2 + 8D0*pap2*pbp1**2 + 6D0*m1**2*m3**2*pbp2 - 2D0*m1**2*pap1*pbp2 + 
     &                2D0*m2**2*pap1*pbp2 - 1D0*m3**2*pap1*pbp2 + m4**2*pap1*pbp2 - 3D0*mq1**2*pap1*pbp2 + 8D0*m1**2*pbp1*pbp2 - 8D0*pap1*pbp1*pbp2 - 
     &                1D0*(m4**2 - 2D0*pbp2)*(m4**2*mq1**2*pap1 + m3**2*mq1**2*pap2 + m2**2*(-1D0*m4**2*pap1 + m3**2*pap2) + m3**2*pap2*pbp1 - 1D0*m4**2*pap2*pbp1 + 
     &                   mq1**2*pap2*pbp1 + m1**2*(m4**2*(pap1 - 2D0*pbp1) - 1D0*m3**2*(2D0*m4**2 + pap2 - 2D0*pbp2)) - 1D0*(m3**2 - 1D0*m4**2 + mq1**2)*pap1*pbp2)*xv + 
     &                2D0*kakb*(m1**2*m3**2 - 1D0*m2**2*m3**2 - 1D0*m3**2*m4**2 - 1D0*m3**2*mq1**2 + m4**2*pap1 - 1D0*m3**2*pap2 - 4D0*pap1*pbp1 - 2D0*pap2*pbp1 + 4D0*m3**2*pbp2 - 
     &                   2D0*pap1*pbp2 + (-1D0*m4**2*pap1 + m3**2*(m4**2 + pap2))*(m4**2 - 2D0*pbp2)*xv) - 2D0*k1k2**2*(kakb + m1**2)*(-1D0 + m4**2*xv - 2D0*pbp2*xv) + 
     &                k1k2*(4D0*kakb**2 + kakb*(2D0*m1**2 - 2D0*m2**2 - 1D0*m3**2 + m4**2 + mq1**2 + 4D0*pap1 - 8D0*pbp1 + 
     &                      (m3**2 - 1D0*m4**2 + mq1**2 - 4D0*pap1)*(m4**2 - 2D0*pbp2)*xv) + 
     &                   m1**2*(-4D0*m2**2 - 1D0*m3**2 + m4**2 + mq1**2 + 2D0*pap1 - 8D0*pbp1 + (m3**2 - 1D0*m4**2 - 1D0*mq1**2 - 2D0*pap1 + 4D0*pbp1)*(m4**2 - 2D0*pbp2)*xv) + 
     &                   2D0*(-1D0*pap2*pbp1 + pap1*(m2**2 - 1D0*mq1**2 + pbp2) + (pap2*pbp1 + pap1*(m2**2 + mq1**2 - 1D0*pbp2))*(m4**2 - 2D0*pbp2)*xv)))) + 
     &          A2Rc*m1*(B1L*mq1*(2D0*k1k2**2 + 4D0*k1k2*kakb + k1k2*m1**2 - 1D0*k1k2*m2**2 - 1D0*k1k2*m3**2 - 4D0*m2**2*m3**2 + k1k2*m4**2 - 2D0*m3**2*m4**2 + 4D0*m2**2*pap1 + 
     &                2D0*m4**2*pap1 - 8D0*k1k2*pbp1 - 4D0*kakb*pbp1 - 2D0*m1**2*pbp1 + 2D0*m2**2*pbp1 + 2D0*m3**2*pbp1 - 2D0*m4**2*pbp1 - 2D0*pap2*pbp1 + 8D0*pbp1**2 + 
     &                6D0*m3**2*pbp2 - 6D0*pap1*pbp2 - 1D0*(2D0*k1k2**2 + k1k2*(m1**2 - 1D0*m2**2 - 1D0*m3**2 + m4**2 - 4D0*pbp1) + 
     &                   2D0*(pap2*pbp1 - 1D0*(m3**2 - 1D0*pap1)*(m4**2 - 1D0*pbp2)))*(m4**2 - 2D0*pbp2)*xv) + 
     &             B1R*m3*(2D0*kakb*(m3 - 1D0*m4)*(m3 + m4) + m4**2*mq1**2 - 2D0*m4**2*pap1 + m4**2*pap2 + mq1**2*pap2 - 2D0*m4**2*pbp1 + 2D0*mq1**2*pbp1 - 6D0*pap2*pbp1 + 
     &                m4**2*pbp2 - 3D0*mq1**2*pbp2 + 6D0*pap1*pbp2 - 1D0*m3**2*(pap2 + pbp2) + m1**2*(m4**2 - 2D0*(pbp1 + pbp2)) - 2D0*kakb*(mq1**2 - 2D0*(pbp1 + pbp2)) - 
     &                1D0*(m4**2 - 2D0*pbp2)*(m1**2*m4**2 + mq1**2*pap2 - 2D0*pap2*pbp1 - 1D0*mq1**2*pbp2 + 2D0*pap1*pbp2 - 1D0*m3**2*(pap2 + pbp2) + 
     &                   m4**2*(mq1**2 - 2D0*pap1 + pap2 - 2D0*pbp1 + pbp2))*xv + k1k2*(m1**2 - 1D0*m2**2 - 1D0*mq1**2 + 2D0*(pap2 + pbp2))*(1D0 - 1D0*m4**2*xv + 2D0*pbp2*xv) + 
     &                m2**2*(2D0*m3**2 + 2D0*(mq1**2 - 2D0*pap1 - 2D0*pap2 + pbp1 + pbp2) + m4**4*xv - 1D0*m4**2*(3D0 + 2D0*pbp2*xv)))))))/
     &   ((2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2 - cI*mq1*wquark2)*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark2)*(m1**2 + m3**2 - 1D0*mq2c**2 - 2D0*pap1 + (cI*mq2c*wst(ptmeigen+sumTc)))*(m2**2 + m4**2 - 1D0*mq2c**2 - 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc)))))


        radvSpTf = dreal((2D0*C2Lc*GS**2*(A1R*(A2Lc*m1*(B1L*m3*(-2D0*k1k2*kak3 + 2D0*k1k2*kakb - 2D0*k1k2*kbk3 + 2D0*k1k2*m2**2 + 2D0*kakb*m3**2 + 2D0*m2**2*m3**2 - 2D0*kak3*m4**2 - 
     &                2D0*kbk3*m4**2 - 1D0*k1k2*mq1**2 + k2k3*mq1**2 - 2D0*kbk3*mq1**2 - 1D0*m4**2*pap1 - 1D0*m3**2*pap2 + 4D0*kak3*pbp1 + 4D0*kbk3*pbp1 - 1D0*m4**2*pbp1 + 
     &                2D0*mq1**2*pbp1 - 2D0*pap2*pbp1 + 4D0*kak3*pbp2 + 4D0*kbk3*pbp2 - 1D0*m3**2*pbp2 + 2D0*pap1*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*(-1D0*k2k3*mq1**2 + k1k2*(2D0*(kak3 + kbk3) + mq1**2) + m4**2*(2D0*kak3 + 2D0*kbk3 + pap1 + pbp1) + m3**2*(pap2 + pbp2))*xv) + 
     &             B1R*mq1*(-1D0*k2k3*m3**2 - 2D0*kakb*m3**2 + 2D0*kbk3*m3**2 - 2D0*m2**2*m3**2 - 1D0*m3**2*m4**2 + m3**2*pap2 - 2D0*k2k3*pbp1 + 4D0*kak3*pbp1 + 4D0*kbk3*pbp1 + 
     &                2D0*m3**2*pbp1 + 4D0*pap1*pbp1 + 4D0*pbp1**2 + 3D0*m3**2*pbp2 + k1k3*(m4**2 - 2D0*(2D0*pbp1 + pbp2)) + 
     &                (-1D0*k1k3*m4**2 + m3**2*(k2k3 + m4**2 - 1D0*pap2 - 1D0*pbp2))*(m4**2 - 2D0*pbp2)*xv + 
     &                k1k2*(2D0*k1k3 - 2D0*kak3 - 1D0*m3**2 - 2D0*pap1 - 2D0*pbp1 + (-2D0*k1k3 + 2D0*kak3 + 2D0*kbk3 + m3**2 + 2D0*(pap1 + pbp1))*(m4**2 - 2D0*pbp2)*xv))) + 
     &          A2Rc*(B1R*m3*mq1*(-2D0*k1k2*kak3 + k1k2*kakb + k2k3*kakb - 1D0*k1k2*m1**2 + k2k3*m1**2 - 2D0*kbk3*m1**2 - 2D0*kak3*m2**2 + 2D0*kakb*m3**2 - 2D0*kak3*m4**2 - 
     &                2D0*m2**2*pap1 - 1D0*m4**2*pap1 - 1D0*kbk3*pap2 - 1D0*m3**2*pap2 + 4D0*kak3*pbp1 + 4D0*kakb*pbp1 + 2D0*m1**2*pbp1 - 3D0*pap2*pbp1 + 5D0*kak3*pbp2 + 
     &                3D0*pap1*pbp2 + (m4**2 - 2D0*pbp2)*(-1D0*k2k3*(kakb + m1**2) + k1k2*(2D0*kak3 + kakb + m1**2) + m4**2*(2D0*kak3 + pap1) + pap2*(kbk3 + m3**2 + pbp1) - 
     &                   1D0*(kak3 + pap1)*pbp2)*xv) + B1L*(2D0*kbk3*m1**2*m3**2 + 2D0*kak3*m2**2*m3**2 + k1k3*kakb*m4**2 + k1k3*m1**2*m4**2 - 1D0*kakb*m3**2*m4**2 - 
     &                1D0*m1**2*m3**2*m4**2 - 2D0*kakb*m3**2*mq1**2 + 2D0*m2**2*m3**2*pap1 - 1D0*kbk3*m4**2*pap1 - 2D0*k1k3*m2**2*pap2 + kbk3*m3**2*pap2 - 2D0*m2**2*m3**2*pap2 + 
     &                m3**2*mq1**2*pap2 - 4D0*k1k3*kakb*pbp1 - 4D0*k1k3*m1**2*pbp1 + 2D0*m1**2*m3**2*pbp1 + kak3*m4**2*pbp1 + 4D0*kak3*mq1**2*pbp1 + 4D0*kbk3*pap1*pbp1 + 
     &                4D0*mq1**2*pap1*pbp1 + 4D0*kbk3*pap2*pbp1 + m3**2*pap2*pbp1 - 4D0*kak3*pbp1**2 + 4D0*pap2*pbp1**2 - 2D0*k1k3*m1**2*pbp2 - 1D0*kak3*m3**2*pbp2 + 
     &                4D0*kakb*m3**2*pbp2 + 2D0*m1**2*m3**2*pbp2 - 1D0*m3**2*pap1*pbp2 - 4D0*kak3*pbp1*pbp2 - 4D0*pap1*pbp1*pbp2 - 
     &                1D0*(m4**2 - 2D0*pbp2)*(k1k3*(kakb + m1**2)*m4**2 - 1D0*kakb*m3**2*m4**2 - 1D0*m1**2*m3**2*m4**2 - 1D0*kbk3*m4**2*pap1 + kbk3*m3**2*pap2 + 
     &                   m3**2*mq1**2*pap2 + kak3*m4**2*pbp1 + m3**2*pap2*pbp1 - 1D0*m3**2*(kak3 + pap1)*pbp2)*xv + 
     &                k2k3*(-1D0*(kakb + m1**2)*m3**2 + 2D0*m2**2*pap1 - 2D0*(2D0*kakb + m1**2)*pbp1 + (kakb + m1**2)*m3**2*(m4**2 - 2D0*pbp2)*xv) + 
     &                k1k2*(2D0*k1k3*kakb + 2D0*k1k3*m1**2 + 2D0*kbk3*m1**2 - 1D0*kakb*m3**2 - 1D0*m1**2*m3**2 - 2D0*kbk3*pap1 + 4D0*m2**2*pap1 - 2D0*mq1**2*pap1 - 4D0*kakb*pbp1 + 
     &                   2D0*kak3*(m2**2 - 1D0*mq1**2 + pbp1) - 1D0*((kakb + m1**2)*(2D0*k1k3 - 1D0*m3**2) - 2D0*(kbk3*pap1 + mq1**2*(kak3 + pap1) - 1D0*kak3*pbp1))*
     &                    (m4**2 - 2D0*pbp2)*xv)))) + A1L*(A2Rc*m1*(B1R*m3*
     &              (-2D0*k1k2*kak3 + 2D0*k1k2*kakb - 2D0*k1k2*kbk3 + 2D0*k1k2*m2**2 + 2D0*kakb*m3**2 + 2D0*m2**2*m3**2 - 2D0*kak3*m4**2 - 2D0*kbk3*m4**2 - 1D0*k1k2*mq1**2 + 
     &                k2k3*mq1**2 - 2D0*kbk3*mq1**2 - 1D0*m4**2*pap1 - 1D0*m3**2*pap2 + 4D0*kak3*pbp1 + 4D0*kbk3*pbp1 - 1D0*m4**2*pbp1 + 2D0*mq1**2*pbp1 - 2D0*pap2*pbp1 + 
     &                4D0*kak3*pbp2 + 4D0*kbk3*pbp2 - 1D0*m3**2*pbp2 + 2D0*pap1*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*(-1D0*k2k3*mq1**2 + k1k2*(2D0*(kak3 + kbk3) + mq1**2) + m4**2*(2D0*kak3 + 2D0*kbk3 + pap1 + pbp1) + m3**2*(pap2 + pbp2))*xv) + 
     &             B1L*mq1*(-1D0*k2k3*m3**2 - 2D0*kakb*m3**2 + 2D0*kbk3*m3**2 - 2D0*m2**2*m3**2 - 1D0*m3**2*m4**2 + m3**2*pap2 - 2D0*k2k3*pbp1 + 4D0*kak3*pbp1 + 4D0*kbk3*pbp1 + 
     &                2D0*m3**2*pbp1 + 4D0*pap1*pbp1 + 4D0*pbp1**2 + 3D0*m3**2*pbp2 + k1k3*(m4**2 - 2D0*(2D0*pbp1 + pbp2)) + 
     &                (-1D0*k1k3*m4**2 + m3**2*(k2k3 + m4**2 - 1D0*pap2 - 1D0*pbp2))*(m4**2 - 2D0*pbp2)*xv + 
     &                k1k2*(2D0*k1k3 - 2D0*kak3 - 1D0*m3**2 - 2D0*pap1 - 2D0*pbp1 + (-2D0*k1k3 + 2D0*kak3 + 2D0*kbk3 + m3**2 + 2D0*(pap1 + pbp1))*(m4**2 - 2D0*pbp2)*xv))) + 
     &          A2Lc*(B1L*m3*mq1*(-2D0*k1k2*kak3 + k1k2*kakb + k2k3*kakb - 1D0*k1k2*m1**2 + k2k3*m1**2 - 2D0*kbk3*m1**2 - 2D0*kak3*m2**2 + 2D0*kakb*m3**2 - 2D0*kak3*m4**2 - 
     &                2D0*m2**2*pap1 - 1D0*m4**2*pap1 - 1D0*kbk3*pap2 - 1D0*m3**2*pap2 + 4D0*kak3*pbp1 + 4D0*kakb*pbp1 + 2D0*m1**2*pbp1 - 3D0*pap2*pbp1 + 5D0*kak3*pbp2 + 
     &                3D0*pap1*pbp2 + (m4**2 - 2D0*pbp2)*(-1D0*k2k3*(kakb + m1**2) + k1k2*(2D0*kak3 + kakb + m1**2) + m4**2*(2D0*kak3 + pap1) + pap2*(kbk3 + m3**2 + pbp1) - 
     &                   1D0*(kak3 + pap1)*pbp2)*xv) + B1R*(2D0*kbk3*m1**2*m3**2 + 2D0*kak3*m2**2*m3**2 + k1k3*kakb*m4**2 + k1k3*m1**2*m4**2 - 1D0*kakb*m3**2*m4**2 - 
     &                1D0*m1**2*m3**2*m4**2 - 2D0*kakb*m3**2*mq1**2 + 2D0*m2**2*m3**2*pap1 - 1D0*kbk3*m4**2*pap1 - 2D0*k1k3*m2**2*pap2 + kbk3*m3**2*pap2 - 2D0*m2**2*m3**2*pap2 + 
     &                m3**2*mq1**2*pap2 - 4D0*k1k3*kakb*pbp1 - 4D0*k1k3*m1**2*pbp1 + 2D0*m1**2*m3**2*pbp1 + kak3*m4**2*pbp1 + 4D0*kak3*mq1**2*pbp1 + 4D0*kbk3*pap1*pbp1 + 
     &                4D0*mq1**2*pap1*pbp1 + 4D0*kbk3*pap2*pbp1 + m3**2*pap2*pbp1 - 4D0*kak3*pbp1**2 + 4D0*pap2*pbp1**2 - 2D0*k1k3*m1**2*pbp2 - 1D0*kak3*m3**2*pbp2 + 
     &                4D0*kakb*m3**2*pbp2 + 2D0*m1**2*m3**2*pbp2 - 1D0*m3**2*pap1*pbp2 - 4D0*kak3*pbp1*pbp2 - 4D0*pap1*pbp1*pbp2 - 
     &                1D0*(m4**2 - 2D0*pbp2)*(k1k3*(kakb + m1**2)*m4**2 - 1D0*kakb*m3**2*m4**2 - 1D0*m1**2*m3**2*m4**2 - 1D0*kbk3*m4**2*pap1 + kbk3*m3**2*pap2 + 
     &                   m3**2*mq1**2*pap2 + kak3*m4**2*pbp1 + m3**2*pap2*pbp1 - 1D0*m3**2*(kak3 + pap1)*pbp2)*xv + 
     &                k2k3*(-1D0*(kakb + m1**2)*m3**2 + 2D0*m2**2*pap1 - 2D0*(2D0*kakb + m1**2)*pbp1 + (kakb + m1**2)*m3**2*(m4**2 - 2D0*pbp2)*xv) + 
     &                k1k2*(2D0*k1k3*kakb + 2D0*k1k3*m1**2 + 2D0*kbk3*m1**2 - 1D0*kakb*m3**2 - 1D0*m1**2*m3**2 - 2D0*kbk3*pap1 + 4D0*m2**2*pap1 - 2D0*mq1**2*pap1 - 4D0*kakb*pbp1 + 
     &                   2D0*kak3*(m2**2 - 1D0*mq1**2 + pbp1) - 1D0*((kakb + m1**2)*(2D0*k1k3 - 1D0*m3**2) - 2D0*(kbk3*pap1 + mq1**2*(kak3 + pap1) - 1D0*kak3*pbp1))*
     &                    (m4**2 - 2D0*pbp2)*xv))))))/
     &   (k1k3*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2 - cI*mq1*wquark2)*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark2)*(m2**2 + m4**2 - 1D0*mq2c**2 - 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc)))))


        radvSpTv = dreal((-2D0*CLL2c*GS*(A1R*(2D0*A2Rc*(B1R*m3*mq1*(-2D0*kakb - 2D0*m1**2 + pap1 + pap2) + 
     &             B1L*((kakb + m1**2)*(k1k2 + m3**2) - 2D0*mq1**2*pap1 - 1D0*pap2*pbp1 + pap1*pbp2)) + 
     &          A2Rc*(B1R*m3*mq1*(2D0*k1k2*pap2 + m4**2*(kakb + m1**2 - 1D0*pap1 + pap2)) + 
     &             B1L*(2D0*k1k2**2*(kakb + m1**2) + k1k2*((kakb + m1**2)*m4**2 - 2D0*pap2*pbp1 + 2D0*pap1*pbp2) - 
     &                1D0*m4**2*((kakb + m1**2)*m3**2 + pap2*pbp1 - 1D0*pap1*(mq1**2 + pbp2))))*xv + 
     &          A2Lc*m1*(B1R*mq1*(2D0*(k1k2 + m3**2 - 2D0*(pap1 + pbp1)) + (2D0*k1k2**2 + k1k2*m4**2 + m4**2*(-1D0*m3**2 + pap1 + pbp1))*xv) + 
     &             B1L*m3*(-4D0*mq1**2 + 2D0*(pap1 + pap2 + pbp1 + pbp2) + (2D0*k1k2*(pap2 + pbp2) + m4**2*(mq1**2 - 1D0*pap1 + pap2 - 1D0*pbp1 + pbp2))*xv))) + 
     &       A1L*(2D0*A2Lc*(B1L*m3*mq1*(-2D0*kakb - 2D0*m1**2 + pap1 + pap2) + B1R*((kakb + m1**2)*(k1k2 + m3**2) - 2D0*mq1**2*pap1 - 1D0*pap2*pbp1 + pap1*pbp2)) + 
     &          A2Lc*(B1L*m3*mq1*(2D0*k1k2*pap2 + m4**2*(kakb + m1**2 - 1D0*pap1 + pap2)) + 
     &             B1R*(2D0*k1k2**2*(kakb + m1**2) + k1k2*((kakb + m1**2)*m4**2 - 2D0*pap2*pbp1 + 2D0*pap1*pbp2) - 
     &                1D0*m4**2*((kakb + m1**2)*m3**2 + pap2*pbp1 - 1D0*pap1*(mq1**2 + pbp2))))*xv + 
     &          A2Rc*m1*(B1L*mq1*(2D0*(k1k2 + m3**2 - 2D0*(pap1 + pbp1)) + (2D0*k1k2**2 + k1k2*m4**2 + m4**2*(-1D0*m3**2 + pap1 + pbp1))*xv) + 
     &             B1R*m3*(-4D0*mq1**2 + 2D0*(pap1 + pap2 + pbp1 + pbp2) + (2D0*k1k2*(pap2 + pbp2) + m4**2*(mq1**2 - 1D0*pap1 + pap2 - 1D0*pbp1 + pbp2))*xv)))))/
     &   ((2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2 - cI*mq1*wquark1)*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark1)*(2D0*k2k3 - 2D0*kbk3 + m2**2 + m4**2 - 1D0*mq2c**2 - 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc)))))


        radvSfTi = dreal((2D0*C2Lc*GS**2*(A1R*(A2Rc*(B1R*m3*mq1*(2D0*k1k3*kakb + k2k3*kakb + 2D0*k1k3*m1**2 - 2D0*kbk3*m1**2 - 2D0*k1k3*pap1 + 2D0*kbk3*pap1 - 1D0*k1k3*pap2 + kbk3*pap2 - 
     &                4D0*m1**2*pbp1 + 4D0*pap1*pbp1 + 2D0*pap2*pbp1 - 1D0*kak3*(2D0*pbp1 + pbp2) + 
     &                (2D0*k1k2 + m4**2 - 2D0*pap2)*(-1D0*k2k3*kakb + k1k3*pap2 - 1D0*pap2*(kbk3 + 2D0*pbp1) + kak3*pbp2)*xv) + 
     &             B1L*(-1D0*k1k2*k1k3*kakb + 4D0*k1k3*kakb**2 + 2D0*k1k2*kakb*kbk3 - 1D0*k1k2*k1k3*m1**2 + 2D0*k1k3*kakb*m1**2 + k1k2*kbk3*m1**2 - 1D0*k1k2*kak3*m2**2 - 
     &                2D0*k1k3*m1**2*m2**2 - 2D0*k1k3*kakb*m3**2 + 4D0*kakb*kbk3*m3**2 - 2D0*k1k3*m1**2*m3**2 + 2D0*kbk3*m1**2*m3**2 - 2D0*kak3*m2**2*m3**2 + 4D0*k1k3*kakb*pap1 - 
     &                4D0*kakb*kbk3*pap1 + 2D0*k1k3*m1**2*pap1 - 2D0*kbk3*m1**2*pap1 + 4D0*k1k3*m2**2*pap1 + k2k3*m2**2*pap1 + k1k3*m2**2*pap2 + 2D0*k1k2*kakb*pbp1 - 
     &                4D0*k1k3*kakb*pbp1 - 4D0*kak3*kakb*pbp1 + 2D0*k1k2*m1**2*pbp1 - 2D0*k1k3*m1**2*pbp1 + k2k3*m1**2*pbp1 - 2D0*kak3*m1**2*pbp1 + 4D0*kbk3*m1**2*pbp1 + 
     &                4D0*kakb*m3**2*pbp1 + 4D0*m1**2*m3**2*pbp1 - 8D0*kakb*pap1*pbp1 - 4D0*kbk3*pap1*pbp1 - 4D0*m1**2*pap1*pbp1 + k1k3*pap2*pbp1 - 2D0*kbk3*pap2*pbp1 + 
     &                4D0*kak3*pbp1**2 + 4D0*m1**2*pbp1**2 - 2D0*pap2*pbp1**2 - 2D0*k1k3*kakb*pbp2 - 1D0*k1k3*m1**2*pbp2 - 1D0*k1k3*pap1*pbp2 + 2D0*kak3*pbp1*pbp2 + 
     &                2D0*pap1*pbp1*pbp2 + (2D0*k1k2 + m4**2 - 2D0*pap2)*
     &                 (-1D0*k2k3*(m2**2*pap1 + m1**2*pbp1) - 1D0*pap2*(-2D0*pbp1*(kbk3 + pbp1) + k1k3*(m2**2 + pbp1)) + 
     &                   k1k2*(k1k3*(kakb + m1**2) + kak3*m2**2 - 2D0*kakb*(kbk3 + pbp1) - 1D0*m1**2*(kbk3 + 2D0*pbp1)) + 
     &                   (k1k3*(2D0*kakb + m1**2 + pap1) - 2D0*(kak3 + pap1)*pbp1)*pbp2)*xv)) + 
     &          A2Lc*m1*(B1R*mq1*(-1D0*k1k2*k1k3 + 2D0*k1k3*kakb + k1k2*kbk3 - 2D0*k1k3*m3**2 + 2D0*kbk3*m3**2 + 2D0*k1k3*pap1 - 2D0*kbk3*pap1 + 2D0*k1k2*pbp1 + k2k3*pbp1 - 
     &                2D0*kak3*pbp1 + 4D0*m3**2*pbp1 - 4D0*pap1*pbp1 - 1D0*k1k3*pbp2 + 
     &                (2D0*k1k2 + m4**2 - 2D0*pap2)*(k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*k2k3*pbp1 + k1k3*pbp2)*xv) + 
     &             B1L*m3*(4D0*k1k3*kakb + k2k3*kakb + 2D0*k1k3*m1**2 - 2D0*kbk3*m1**2 + 2D0*k1k3*m2**2 + k2k3*m2**2 - 2D0*kak3*m2**2 - 2D0*k1k3*pap1 + 2D0*kbk3*pap1 - 
     &                1D0*k1k3*pap2 + kbk3*pap2 - 2D0*k1k3*pbp1 - 2D0*kak3*pbp1 - 4D0*kakb*pbp1 - 4D0*m1**2*pbp1 + 4D0*pap1*pbp1 + 2D0*pap2*pbp1 + 4D0*pbp1**2 - 1D0*k1k3*pbp2 - 
     &                1D0*kak3*pbp2 + 2D0*pbp1*pbp2 + (2D0*k1k2 + m4**2 - 2D0*pap2)*
     &                 (-1D0*k2k3*(kakb + m2**2) - 1D0*kbk3*pap2 + kak3*pbp2 + k1k3*(pap2 + pbp2) - 2D0*pbp1*(pap2 + pbp2))*xv))) + 
     &       A1L*(A2Lc*(B1L*m3*mq1*(2D0*k1k3*kakb + k2k3*kakb + 2D0*k1k3*m1**2 - 2D0*kbk3*m1**2 - 2D0*k1k3*pap1 + 2D0*kbk3*pap1 - 1D0*k1k3*pap2 + kbk3*pap2 - 4D0*m1**2*pbp1 + 
     &                4D0*pap1*pbp1 + 2D0*pap2*pbp1 - 1D0*kak3*(2D0*pbp1 + pbp2) + 
     &                (2D0*k1k2 + m4**2 - 2D0*pap2)*(-1D0*k2k3*kakb + k1k3*pap2 - 1D0*pap2*(kbk3 + 2D0*pbp1) + kak3*pbp2)*xv) + 
     &             B1R*(-1D0*k1k2*k1k3*kakb + 4D0*k1k3*kakb**2 + 2D0*k1k2*kakb*kbk3 - 1D0*k1k2*k1k3*m1**2 + 2D0*k1k3*kakb*m1**2 + k1k2*kbk3*m1**2 - 1D0*k1k2*kak3*m2**2 - 
     &                2D0*k1k3*m1**2*m2**2 - 2D0*k1k3*kakb*m3**2 + 4D0*kakb*kbk3*m3**2 - 2D0*k1k3*m1**2*m3**2 + 2D0*kbk3*m1**2*m3**2 - 2D0*kak3*m2**2*m3**2 + 4D0*k1k3*kakb*pap1 - 
     &                4D0*kakb*kbk3*pap1 + 2D0*k1k3*m1**2*pap1 - 2D0*kbk3*m1**2*pap1 + 4D0*k1k3*m2**2*pap1 + k2k3*m2**2*pap1 + k1k3*m2**2*pap2 + 2D0*k1k2*kakb*pbp1 - 
     &                4D0*k1k3*kakb*pbp1 - 4D0*kak3*kakb*pbp1 + 2D0*k1k2*m1**2*pbp1 - 2D0*k1k3*m1**2*pbp1 + k2k3*m1**2*pbp1 - 2D0*kak3*m1**2*pbp1 + 4D0*kbk3*m1**2*pbp1 + 
     &                4D0*kakb*m3**2*pbp1 + 4D0*m1**2*m3**2*pbp1 - 8D0*kakb*pap1*pbp1 - 4D0*kbk3*pap1*pbp1 - 4D0*m1**2*pap1*pbp1 + k1k3*pap2*pbp1 - 2D0*kbk3*pap2*pbp1 + 
     &                4D0*kak3*pbp1**2 + 4D0*m1**2*pbp1**2 - 2D0*pap2*pbp1**2 - 2D0*k1k3*kakb*pbp2 - 1D0*k1k3*m1**2*pbp2 - 1D0*k1k3*pap1*pbp2 + 2D0*kak3*pbp1*pbp2 + 
     &                2D0*pap1*pbp1*pbp2 + (2D0*k1k2 + m4**2 - 2D0*pap2)*
     &                 (-1D0*k2k3*(m2**2*pap1 + m1**2*pbp1) - 1D0*pap2*(-2D0*pbp1*(kbk3 + pbp1) + k1k3*(m2**2 + pbp1)) + 
     &                   k1k2*(k1k3*(kakb + m1**2) + kak3*m2**2 - 2D0*kakb*(kbk3 + pbp1) - 1D0*m1**2*(kbk3 + 2D0*pbp1)) + 
     &                   (k1k3*(2D0*kakb + m1**2 + pap1) - 2D0*(kak3 + pap1)*pbp1)*pbp2)*xv)) + 
     &          A2Rc*m1*(B1L*mq1*(-1D0*k1k2*k1k3 + 2D0*k1k3*kakb + k1k2*kbk3 - 2D0*k1k3*m3**2 + 2D0*kbk3*m3**2 + 2D0*k1k3*pap1 - 2D0*kbk3*pap1 + 2D0*k1k2*pbp1 + k2k3*pbp1 - 
     &                2D0*kak3*pbp1 + 4D0*m3**2*pbp1 - 4D0*pap1*pbp1 - 1D0*k1k3*pbp2 + 
     &                (2D0*k1k2 + m4**2 - 2D0*pap2)*(k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*k2k3*pbp1 + k1k3*pbp2)*xv) + 
     &             B1R*m3*(4D0*k1k3*kakb + k2k3*kakb + 2D0*k1k3*m1**2 - 2D0*kbk3*m1**2 + 2D0*k1k3*m2**2 + k2k3*m2**2 - 2D0*kak3*m2**2 - 2D0*k1k3*pap1 + 2D0*kbk3*pap1 - 
     &                1D0*k1k3*pap2 + kbk3*pap2 - 2D0*k1k3*pbp1 - 2D0*kak3*pbp1 - 4D0*kakb*pbp1 - 4D0*m1**2*pbp1 + 4D0*pap1*pbp1 + 2D0*pap2*pbp1 + 4D0*pbp1**2 - 1D0*k1k3*pbp2 - 
     &                1D0*kak3*pbp2 + 2D0*pbp1*pbp2 + (2D0*k1k2 + m4**2 - 2D0*pap2)*
     &                 (-1D0*k2k3*(kakb + m2**2) - 1D0*kbk3*pap2 + kak3*pbp2 + k1k3*(pap2 + pbp2) - 2D0*pbp1*(pap2 + pbp2))*xv)))))/
     &   (k1k3*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2 - cI*mq1*wquark1)*(m1**2 + m3**2 - 1D0*mq2c**2 - 2D0*pap1 + (cI*mq2c*wst(ptmeigen+sumTc)))*(2D0*k1k2 + m1**2 - 1D0*m2**2 + m3**2 + m4**2 - 2D0*(pap1 + pap2))))


        radvSfTp = dreal((C2Lc*GS**2*(A1R*(A2Lc*m1*(B1R*mq1*(2D0*k1k2**2 + 2D0*k1k2*k2k3 + k1k2*kak3 - 3D0*k1k2*kbk3 - 2D0*k1k2*m3**2 - 1D0*k2k3*m3**2 + 2D0*kbk3*m3**2 + 2D0*k1k2*pap1 + 
     &                k2k3*pap1 - 2D0*kbk3*pap1 - 6D0*k1k2*pbp1 - 3D0*k2k3*pbp1 - 2D0*kak3*pbp1 + 4D0*kbk3*pbp1 + 4D0*m3**2*pbp1 - 4D0*pap1*pbp1 + 4D0*pbp1**2 + 
     &                k1k3*(2D0*kakb - 2D0*m2**2 - 1D0*m4**2 - 1D0*pap2 + 3D0*pbp2) + 
     &                (-2D0*k1k2**2 + k2k3*(m3**2 - 1D0*pap1 + pbp1) + k1k2*(-2D0*k2k3 - 1D0*kak3 + kbk3 + 2D0*(m3**2 - 1D0*pap1 + pbp1)) + k1k3*(m4**2 + pap2 - 1D0*pbp2))*
     &                 (m4**2 - 2D0*pbp2)*xv) + B1L*m3*(k1k2*kak3 - 4D0*k1k2*kakb + 2D0*k1k3*kakb - 2D0*k2k3*kakb + k1k2*kbk3 + 4D0*kakb*kbk3 + k2k3*m1**2 - 2D0*kbk3*m1**2 - 
     &                4D0*k1k2*m2**2 + 2D0*k1k3*m2**2 - 3D0*k2k3*m2**2 - 4D0*kak3*m2**2 + 2D0*kbk3*m2**2 + 4D0*kakb*m3**2 + 4D0*m2**2*m3**2 - 1D0*kak3*m4**2 - 1D0*kbk3*m4**2 - 
     &                1D0*k2k3*pap1 - 4D0*kakb*pap1 + 2D0*kbk3*pap1 - 4D0*m2**2*pap1 + 2D0*k1k2*pap2 - 1D0*k1k3*pap2 + 2D0*k2k3*pap2 - 4D0*kbk3*pap2 - 2D0*m3**2*pap2 + 
     &                2D0*pap1*pap2 - 1D0*k2k3*pbp1 - 2D0*kak3*pbp1 + 4D0*kakb*pbp1 + 4D0*m2**2*pbp1 - 2D0*pap2*pbp1 + 2D0*k1k2*pbp2 - 1D0*k1k3*pbp2 + 2D0*k2k3*pbp2 + 
     &                4D0*kak3*pbp2 - 2D0*m3**2*pbp2 + 2D0*pap1*pbp2 - 2D0*pbp1*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*((kak3 + kbk3)*m4**2 + 2D0*pap2*(kbk3 + m3**2 - 1D0*pap1 + pbp1) + k2k3*(-1D0*m1**2 + m2**2 + pap1 - 2D0*pap2 + pbp1 - 2D0*pbp2) + 
     &                   2D0*(-1D0*kak3 + m3**2 - 1D0*pap1 + pbp1)*pbp2 + k1k3*(pap2 + pbp2) - 1D0*k1k2*(kak3 + kbk3 + 2D0*(pap2 + pbp2)))*xv)) + 
     &          A2Rc*(B1R*m3*mq1*(k1k2*kak3 - 4D0*k1k2*kakb + 2D0*k1k3*kakb - 3D0*k2k3*kakb + 4D0*kakb*kbk3 + k2k3*m1**2 - 2D0*kbk3*m1**2 - 2D0*kak3*m2**2 + 4D0*kakb*m3**2 - 
     &                1D0*kak3*m4**2 - 1D0*k2k3*pap1 - 4D0*kakb*pap1 + 2D0*kbk3*pap1 + 2D0*k1k2*pap2 - 1D0*k1k3*pap2 + 2D0*k2k3*pap2 - 3D0*kbk3*pap2 - 2D0*m3**2*pap2 + 
     &                2D0*pap1*pap2 - 2D0*kak3*pbp1 + 4D0*kakb*pbp1 - 2D0*pap2*pbp1 + 3D0*kak3*pbp2 + 
     &                (k2k3*(kakb - 1D0*m1**2 + pap1 - 2D0*pap2) - 1D0*k1k2*(kak3 + 2D0*pap2) + pap2*(k1k3 + kbk3 + 2D0*(m3**2 - 1D0*pap1 + pbp1)) + kak3*(m4**2 - 1D0*pbp2))*
     &                 (m4**2 - 2D0*pbp2)*xv) + B1L*(2D0*k1k3*kakb*m1**2 - 2D0*k1k3*kakb*m2**2 - 1D0*k2k3*kakb*m3**2 - 1D0*k2k3*m1**2*m3**2 + 2D0*kbk3*m1**2*m3**2 + 
     &                2D0*kak3*m2**2*m3**2 - 1D0*k1k3*kakb*m4**2 - 1D0*k1k3*m1**2*m4**2 + 2D0*k2k3*kakb*pap1 + k2k3*m1**2*pap1 - 2D0*kbk3*m1**2*pap1 - 3D0*k2k3*m2**2*pap1 - 
     &                4D0*kak3*m2**2*pap1 + 2D0*kbk3*m2**2*pap1 + 4D0*m2**2*m3**2*pap1 - 1D0*kbk3*m4**2*pap1 - 4D0*m2**2*pap1**2 - 1D0*k1k3*m1**2*pap2 + k1k3*m2**2*pap2 + 
     &                kbk3*m3**2*pap2 - 2D0*kbk3*pap1*pap2 - 4D0*k2k3*m1**2*pbp1 - 2D0*kak3*m1**2*pbp1 + 4D0*kbk3*m1**2*pbp1 + 2D0*kak3*m2**2*pbp1 + 4D0*m1**2*m3**2*pbp1 + 
     &                kak3*m4**2*pbp1 - 4D0*m1**2*pap1*pbp1 + 4D0*m2**2*pap1*pbp1 - 2D0*k2k3*pap2*pbp1 + 2D0*kbk3*pap2*pbp1 + 2D0*m3**2*pap2*pbp1 - 2D0*pap1*pap2*pbp1 + 
     &                4D0*m1**2*pbp1**2 + 2D0*pap2*pbp1**2 + 2D0*k1k3*kakb*pbp2 + 2D0*k1k3*m1**2*pbp2 - 1D0*kak3*m3**2*pbp2 + 2D0*k2k3*pap1*pbp2 + 2D0*kak3*pap1*pbp2 - 
     &                2D0*m3**2*pap1*pbp2 + 2D0*pap1**2*pbp2 - 2D0*kak3*pbp1*pbp2 - 2D0*pap1*pbp1*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*(k2k3*kakb*m3**2 + k2k3*m1**2*m3**2 + k1k3*kakb*m4**2 + k1k3*m1**2*m4**2 - 2D0*k2k3*kakb*pap1 - 1D0*k2k3*m1**2*pap1 + k2k3*m2**2*pap1 + 
     &                   kbk3*m4**2*pap1 + k1k3*m1**2*pap2 + k1k3*m2**2*pap2 - 1D0*kbk3*m3**2*pap2 + 2D0*kbk3*pap1*pap2 + 2D0*k2k3*m1**2*pbp1 - 1D0*kak3*m4**2*pbp1 + 
     &                   2D0*k2k3*pap2*pbp1 - 2D0*kbk3*pap2*pbp1 - 2D0*m3**2*pap2*pbp1 + 2D0*pap1*pap2*pbp1 - 2D0*pap2*pbp1**2 + 
     &                   (-2D0*k1k3*kakb + 2D0*pap1*(-1D0*k2k3 + m3**2 - 1D0*pap1 + pbp1) + kak3*(m3**2 - 2D0*pap1 + 2D0*pbp1))*pbp2)*xv - 
     &                2D0*k1k2**2*(kakb + m1**2)*(-1D0 + m4**2*xv - 2D0*pbp2*xv) + 
     &                k1k2*(kak3*(m1 - 1D0*m2)*(m1 + m2) - 2D0*(kbk3*m1**2 + m1**2*m3**2 - 1D0*m1**2*pap1 + 2D0*m2**2*pap1 + 3D0*m1**2*pbp1 + pap2*pbp1 + 
     &                      kakb*(kbk3 + m3**2 - 1D0*pap1 + pbp1) - 1D0*pap1*pbp2) + 
     &                   (m4**2 - 2D0*pbp2)*(-1D0*kak3*(m1**2 + m2**2) + 2D0*(pap2*pbp1 + m1**2*(m3**2 - 1D0*pap1 + pbp1) + kakb*(kbk3 + m3**2 - 1D0*pap1 + pbp1) - 1D0*pap1*pbp2))*
     &                    xv - 2D0*k2k3*(kakb + m1**2)*(-1D0 + m4**2*xv - 2D0*pbp2*xv))))) + 
     &       A1L*(A2Rc*m1*(B1L*mq1*(2D0*k1k2**2 + 2D0*k1k2*k2k3 + k1k2*kak3 - 3D0*k1k2*kbk3 - 2D0*k1k2*m3**2 - 1D0*k2k3*m3**2 + 2D0*kbk3*m3**2 + 2D0*k1k2*pap1 + k2k3*pap1 - 
     &                2D0*kbk3*pap1 - 6D0*k1k2*pbp1 - 3D0*k2k3*pbp1 - 2D0*kak3*pbp1 + 4D0*kbk3*pbp1 + 4D0*m3**2*pbp1 - 4D0*pap1*pbp1 + 4D0*pbp1**2 + 
     &                k1k3*(2D0*kakb - 2D0*m2**2 - 1D0*m4**2 - 1D0*pap2 + 3D0*pbp2) + 
     &                (-2D0*k1k2**2 + k2k3*(m3**2 - 1D0*pap1 + pbp1) + k1k2*(-2D0*k2k3 - 1D0*kak3 + kbk3 + 2D0*(m3**2 - 1D0*pap1 + pbp1)) + k1k3*(m4**2 + pap2 - 1D0*pbp2))*
     &                 (m4**2 - 2D0*pbp2)*xv) + B1R*m3*(k1k2*kak3 - 4D0*k1k2*kakb + 2D0*k1k3*kakb - 2D0*k2k3*kakb + k1k2*kbk3 + 4D0*kakb*kbk3 + k2k3*m1**2 - 2D0*kbk3*m1**2 - 
     &                4D0*k1k2*m2**2 + 2D0*k1k3*m2**2 - 3D0*k2k3*m2**2 - 4D0*kak3*m2**2 + 2D0*kbk3*m2**2 + 4D0*kakb*m3**2 + 4D0*m2**2*m3**2 - 1D0*kak3*m4**2 - 1D0*kbk3*m4**2 - 
     &                1D0*k2k3*pap1 - 4D0*kakb*pap1 + 2D0*kbk3*pap1 - 4D0*m2**2*pap1 + 2D0*k1k2*pap2 - 1D0*k1k3*pap2 + 2D0*k2k3*pap2 - 4D0*kbk3*pap2 - 2D0*m3**2*pap2 + 
     &                2D0*pap1*pap2 - 1D0*k2k3*pbp1 - 2D0*kak3*pbp1 + 4D0*kakb*pbp1 + 4D0*m2**2*pbp1 - 2D0*pap2*pbp1 + 2D0*k1k2*pbp2 - 1D0*k1k3*pbp2 + 2D0*k2k3*pbp2 + 
     &                4D0*kak3*pbp2 - 2D0*m3**2*pbp2 + 2D0*pap1*pbp2 - 2D0*pbp1*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*((kak3 + kbk3)*m4**2 + 2D0*pap2*(kbk3 + m3**2 - 1D0*pap1 + pbp1) + k2k3*(-1D0*m1**2 + m2**2 + pap1 - 2D0*pap2 + pbp1 - 2D0*pbp2) + 
     &                   2D0*(-1D0*kak3 + m3**2 - 1D0*pap1 + pbp1)*pbp2 + k1k3*(pap2 + pbp2) - 1D0*k1k2*(kak3 + kbk3 + 2D0*(pap2 + pbp2)))*xv)) + 
     &          A2Lc*(B1L*m3*mq1*(k1k2*kak3 - 4D0*k1k2*kakb + 2D0*k1k3*kakb - 3D0*k2k3*kakb + 4D0*kakb*kbk3 + k2k3*m1**2 - 2D0*kbk3*m1**2 - 2D0*kak3*m2**2 + 4D0*kakb*m3**2 - 
     &                1D0*kak3*m4**2 - 1D0*k2k3*pap1 - 4D0*kakb*pap1 + 2D0*kbk3*pap1 + 2D0*k1k2*pap2 - 1D0*k1k3*pap2 + 2D0*k2k3*pap2 - 3D0*kbk3*pap2 - 2D0*m3**2*pap2 + 
     &                2D0*pap1*pap2 - 2D0*kak3*pbp1 + 4D0*kakb*pbp1 - 2D0*pap2*pbp1 + 3D0*kak3*pbp2 + 
     &                (k2k3*(kakb - 1D0*m1**2 + pap1 - 2D0*pap2) - 1D0*k1k2*(kak3 + 2D0*pap2) + pap2*(k1k3 + kbk3 + 2D0*(m3**2 - 1D0*pap1 + pbp1)) + kak3*(m4**2 - 1D0*pbp2))*
     &                 (m4**2 - 2D0*pbp2)*xv) + B1R*(2D0*k1k3*kakb*m1**2 - 2D0*k1k3*kakb*m2**2 - 1D0*k2k3*kakb*m3**2 - 1D0*k2k3*m1**2*m3**2 + 2D0*kbk3*m1**2*m3**2 + 
     &                2D0*kak3*m2**2*m3**2 - 1D0*k1k3*kakb*m4**2 - 1D0*k1k3*m1**2*m4**2 + 2D0*k2k3*kakb*pap1 + k2k3*m1**2*pap1 - 2D0*kbk3*m1**2*pap1 - 3D0*k2k3*m2**2*pap1 - 
     &                4D0*kak3*m2**2*pap1 + 2D0*kbk3*m2**2*pap1 + 4D0*m2**2*m3**2*pap1 - 1D0*kbk3*m4**2*pap1 - 4D0*m2**2*pap1**2 - 1D0*k1k3*m1**2*pap2 + k1k3*m2**2*pap2 + 
     &                kbk3*m3**2*pap2 - 2D0*kbk3*pap1*pap2 - 4D0*k2k3*m1**2*pbp1 - 2D0*kak3*m1**2*pbp1 + 4D0*kbk3*m1**2*pbp1 + 2D0*kak3*m2**2*pbp1 + 4D0*m1**2*m3**2*pbp1 + 
     &                kak3*m4**2*pbp1 - 4D0*m1**2*pap1*pbp1 + 4D0*m2**2*pap1*pbp1 - 2D0*k2k3*pap2*pbp1 + 2D0*kbk3*pap2*pbp1 + 2D0*m3**2*pap2*pbp1 - 2D0*pap1*pap2*pbp1 + 
     &                4D0*m1**2*pbp1**2 + 2D0*pap2*pbp1**2 + 2D0*k1k3*kakb*pbp2 + 2D0*k1k3*m1**2*pbp2 - 1D0*kak3*m3**2*pbp2 + 2D0*k2k3*pap1*pbp2 + 2D0*kak3*pap1*pbp2 - 
     &                2D0*m3**2*pap1*pbp2 + 2D0*pap1**2*pbp2 - 2D0*kak3*pbp1*pbp2 - 2D0*pap1*pbp1*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*(k2k3*kakb*m3**2 + k2k3*m1**2*m3**2 + k1k3*kakb*m4**2 + k1k3*m1**2*m4**2 - 2D0*k2k3*kakb*pap1 - 1D0*k2k3*m1**2*pap1 + k2k3*m2**2*pap1 + 
     &                   kbk3*m4**2*pap1 + k1k3*m1**2*pap2 + k1k3*m2**2*pap2 - 1D0*kbk3*m3**2*pap2 + 2D0*kbk3*pap1*pap2 + 2D0*k2k3*m1**2*pbp1 - 1D0*kak3*m4**2*pbp1 + 
     &                   2D0*k2k3*pap2*pbp1 - 2D0*kbk3*pap2*pbp1 - 2D0*m3**2*pap2*pbp1 + 2D0*pap1*pap2*pbp1 - 2D0*pap2*pbp1**2 + 
     &                   (-2D0*k1k3*kakb + 2D0*pap1*(-1D0*k2k3 + m3**2 - 1D0*pap1 + pbp1) + kak3*(m3**2 - 2D0*pap1 + 2D0*pbp1))*pbp2)*xv - 
     &                2D0*k1k2**2*(kakb + m1**2)*(-1D0 + m4**2*xv - 2D0*pbp2*xv) + 
     &                k1k2*(kak3*(m1 - 1D0*m2)*(m1 + m2) - 2D0*(kbk3*m1**2 + m1**2*m3**2 - 1D0*m1**2*pap1 + 2D0*m2**2*pap1 + 3D0*m1**2*pbp1 + pap2*pbp1 + 
     &                      kakb*(kbk3 + m3**2 - 1D0*pap1 + pbp1) - 1D0*pap1*pbp2) + 
     &                   (m4**2 - 2D0*pbp2)*(-1D0*kak3*(m1**2 + m2**2) + 2D0*(pap2*pbp1 + m1**2*(m3**2 - 1D0*pap1 + pbp1) + kakb*(kbk3 + m3**2 - 1D0*pap1 + pbp1) - 1D0*pap1*pbp2))*
     &                    xv - 2D0*k2k3*(kakb + m1**2)*(-1D0 + m4**2*xv - 2D0*pbp2*xv)))))))/
     &   (k1k3*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2 - cI*mq1*wquark1)*(m1**2 + m3**2 - 1D0*mq2c**2 - 2D0*pap1 + (cI*mq2c*wst(ptmeigen+sumTc)))*(m2**2 + m4**2 - 1D0*mq2c**2 - 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc)))))


        radvSfTf = dreal((2D0*C2Lc*GS**2*(A1R*(A2Lc*m1*(B1L*m3*(k1k3 + m3**2)*(2D0*kakb + 2D0*m2**2 - 1D0*(pap2 + pbp2)*(1D0 - 1D0*m4**2*xv + 2D0*pbp2*xv)) + 
     &             B1R*mq1*(m3**2*(-1D0*k1k2 - 1D0*k2k3 + 2D0*(kbk3 + pbp1) + (k1k2 + k2k3)*(m4**2 - 2D0*pbp2)*xv) + k1k3*(k2k3 - 2D0*kbk3 - 1D0*k2k3*m4**2*xv + 2D0*k2k3*pbp2*xv)))
     &            + A2Rc*(B1R*m3*(k1k3 + m3**2)*mq1*(2D0*kakb + pap2*(-1D0 + m4**2*xv - 2D0*pbp2*xv)) + 
     &             B1L*(k1k3*(k2k3*kakb + k2k3*m1**2 - 2D0*kbk3*m1**2 - 2D0*kak3*m2**2 - 1D0*kbk3*pap2 + kak3*pbp2 - 
     &                   1D0*(m4**2 - 2D0*pbp2)*(k2k3*(kakb + m1**2) - 1D0*kbk3*pap2 + kak3*pbp2)*xv) + 
     &                m3**2*(-1D0*k1k2*kakb - 1D0*k2k3*kakb - 1D0*k1k2*m1**2 - 1D0*k2k3*m1**2 + 2D0*kbk3*m1**2 + kbk3*pap2 + 2D0*m1**2*pbp1 + pap2*pbp1 + 
     &                   (kak3 + pap1)*(2D0*m2**2 - 1D0*pbp2) + (m4**2 - 2D0*pbp2)*((k1k2 + k2k3)*(kakb + m1**2) - 1D0*pap2*(kbk3 + pbp1) + (kak3 + pap1)*pbp2)*xv)))) + 
     &       A1L*(A2Rc*m1*(B1R*m3*(k1k3 + m3**2)*(2D0*kakb + 2D0*m2**2 - 1D0*(pap2 + pbp2)*(1D0 - 1D0*m4**2*xv + 2D0*pbp2*xv)) + 
     &             B1L*mq1*(m3**2*(-1D0*k1k2 - 1D0*k2k3 + 2D0*(kbk3 + pbp1) + (k1k2 + k2k3)*(m4**2 - 2D0*pbp2)*xv) + k1k3*(k2k3 - 2D0*kbk3 - 1D0*k2k3*m4**2*xv + 2D0*k2k3*pbp2*xv)))
     &            + A2Lc*(B1L*m3*(k1k3 + m3**2)*mq1*(2D0*kakb + pap2*(-1D0 + m4**2*xv - 2D0*pbp2*xv)) + 
     &             B1R*(k1k3*(k2k3*kakb + k2k3*m1**2 - 2D0*kbk3*m1**2 - 2D0*kak3*m2**2 - 1D0*kbk3*pap2 + kak3*pbp2 - 
     &                   1D0*(m4**2 - 2D0*pbp2)*(k2k3*(kakb + m1**2) - 1D0*kbk3*pap2 + kak3*pbp2)*xv) + 
     &                m3**2*(-1D0*k1k2*kakb - 1D0*k2k3*kakb - 1D0*k1k2*m1**2 - 1D0*k2k3*m1**2 + 2D0*kbk3*m1**2 + kbk3*pap2 + 2D0*m1**2*pbp1 + pap2*pbp1 + 
     &                   (kak3 + pap1)*(2D0*m2**2 - 1D0*pbp2) + (m4**2 - 2D0*pbp2)*((k1k2 + k2k3)*(kakb + m1**2) - 1D0*pap2*(kbk3 + pbp1) + (kak3 + pap1)*pbp2)*xv))))))/
     &   (k1k3**2*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2 - cI*mq1*wquark1)*(m2**2 + m4**2 - 1D0*mq2c**2 - 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc)))))


        radvSfTv = dreal((CLL2c*GS*(A1L*(-2D0*A2Lc*(B1L*m3*mq1*(kak3 - 1D0*pap1) + B1R*((kakb + m1**2)*(k1k3 - 1D0*m3**2) + kbk3*pap1 - 1D0*kak3*pbp1)) + 
     &          A2Lc*(B1L*m3*mq1*(kak3*m4**2 - 2D0*(k1k2 + k2k3)*pap2) + 
     &             B1R*(-2D0*k1k2**2*(kakb + m1**2) + k1k3*(kakb + m1**2)*m4**2 + kbk3*m4**2*pap1 - 1D0*kak3*m4**2*pbp1 + 2D0*k2k3*pap2*pbp1 - 2D0*k2k3*pap1*pbp2 - 
     &                2D0*k1k2*(k2k3*(kakb + m1**2) - 1D0*pap2*pbp1 + pap1*pbp2)))*xv + 
     &          A2Rc*m1*(B1R*m3*(-2D0*kak3 + 2D0*(-1D0*kbk3 + pap1 + pbp1) + ((kak3 + kbk3)*m4**2 - 2D0*(k1k2 + k2k3)*(pap2 + pbp2))*xv) + 
     &             B1L*mq1*(2D0*m3**2 - 2D0*k1k2*(k1k2 + k2k3)*xv + k1k3*(-2D0 + m4**2*xv)))) + 
     &       A1R*(-2D0*A2Rc*(B1R*m3*mq1*(kak3 - 1D0*pap1) + B1L*((kakb + m1**2)*(k1k3 - 1D0*m3**2) + kbk3*pap1 - 1D0*kak3*pbp1)) + 
     &          A2Rc*(B1R*m3*mq1*(kak3*m4**2 - 2D0*(k1k2 + k2k3)*pap2) + 
     &             B1L*(-2D0*k1k2**2*(kakb + m1**2) + k1k3*(kakb + m1**2)*m4**2 + kbk3*m4**2*pap1 - 1D0*kak3*m4**2*pbp1 + 2D0*k2k3*pap2*pbp1 - 2D0*k2k3*pap1*pbp2 - 
     &                2D0*k1k2*(k2k3*(kakb + m1**2) - 1D0*pap2*pbp1 + pap1*pbp2)))*xv + 
     &          A2Lc*m1*(B1L*m3*(-2D0*kak3 + 2D0*(-1D0*kbk3 + pap1 + pbp1) + ((kak3 + kbk3)*m4**2 - 2D0*(k1k2 + k2k3)*(pap2 + pbp2))*xv) + 
     &             B1R*mq1*(2D0*m3**2 - 2D0*k1k2*(k1k2 + k2k3)*xv + k1k3*(-2D0 + m4**2*xv))))))/
     &   (k1k3*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2 - cI*mq1*wquark1)*(2D0*k2k3 - 2D0*kbk3 + m2**2 + m4**2 - 1D0*mq2c**2 - 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc)))))


            NeuQ2qx_M2RadiationVST = NeuQ2qx_M2RadiationVST + 2.*CF*NC*(radvSiTi + radvSiTp + radvSiTf + radvSiTv + radvSpTi + radvSpTp + radvSpTf + radvSpTv + radvSfTi + radvSfTp + radvSfTf + radvSfTv)


        enddo

        end

! =================================================================
        double precision function NeuQ2qx_M2RadiationVSU()
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double precision radvSiUi, radvSiUf, radvSpUi, radvSpUf, radvSfUi, radvSfUf
        double complex A1L, A1R, A3Lc, A3Rc, B1L, B1R, D3Lc, D3Rc
        double precision mq1, mq3c
        double precision xv, colorfactor
        integer sumUc, sumUend
        double precision wquark


        if (xtype.eq.1) then
        
           xv = 0d0
           colorfactor = 4d0/3d0
           sumUend = -1
           
        elseif (xtype.eq.2) then
        
           xv = 0d0
           colorfactor = 1d0
           sumUend = -1
           
        elseif (xtype.eq.3) then
        
           xv = 1/(m4**2)
           colorfactor = 1d0
           sumUend = 3
    
           
        elseif (xtype.eq.4) then

           xv = 1/(m4**2)
           colorfactor = 1d0
           sumUend = 1
                   
        endif

        NeuQ2qx_M2RadiationVSU = 0d0

        do sumUc = 0, sumUend

        radvSiUi = 0d0
        radvSiUf = 0d0
        radvSpUi = 0d0
        radvSpUf = 0d0
        radvSfUi = 0d0
        radvSfUf = 0d0


        A1L = AL(chimeigen,sfmeigen,sftype,sfgen)
        A1R = AR(chimeigen,sfmeigen,sftype,sfgen)

        B1L = BL(xtype,ftype)
        B1R = BR(xtype,ftype)

        mq1 = Mf(pstype,psgen) 

c        wquark = width(pstype)
        wquark = 0d0

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        D3Lc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        D3Rc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)


        radvSiUi = dreal((-2D0*GS**2*(kbk3 - 1D0*m2**2)*(A1L*(A3Rc*(B1L*mq1*(D3Rc*m1*(k1k2 - 1D0*pap1)*(-4D0 + m4**2*xv) + D3Lc*mq3c*(-2D0*pap1 + m4**2*pap1*xv - 2D0*k1k2*pap2*xv)) + 
     &             B1R*m3*(-1D0*D3Lc*mq3c*(pap1 + pap2)*(-4D0 + m4**2*xv) + 
     &                D3Rc*m1*(-2D0*(pap1 + pap2) + m4**2*(2D0 + (m4**2 + pap1 - 1D0*pap2)*xv) + k1k2*(2D0 + m4**2*xv - 2D0*pap2*xv)))) + 
     &          A3Lc*(-1D0*B1R*(4D0*D3Lc*(k1k2 - 1D0*pap1)*(pap1 + pap2) + 
     &                D3Lc*(-2D0*k1k2**2*m1**2 + m1**2*m3**2*m4**2 - 1D0*k1k2*m4**2*(m1**2 + 2D0*pap1) + 4D0*k1k2*pap1*pap2 - 1D0*(m4**2 - 2D0*pap2)*(m4**2*pap1 - 1D0*m3**2*pap2))*
     &                 xv + D3Rc*m1*mq3c*(2D0*(k1k2 + m3**2) + (2D0*k1k2**2 + (k1k2 - 1D0*m3**2)*m4**2)*xv)) + 
     &             B1L*m3*mq1*(D3Rc*m1*mq3c*(4D0 - 1D0*m4**2*xv) + D3Lc*(m1**2*(-2D0 + m4**2*xv) + pap2*(2D0 + m4**2*xv - 2D0*pap2*xv))))) + 
     &       A1R*(A3Lc*(B1R*mq1*(D3Lc*m1*(k1k2 - 1D0*pap1)*(-4D0 + m4**2*xv) + D3Rc*mq3c*(-2D0*pap1 + m4**2*pap1*xv - 2D0*k1k2*pap2*xv)) + 
     &             B1L*m3*(-1D0*D3Rc*mq3c*(pap1 + pap2)*(-4D0 + m4**2*xv) + 
     &                D3Lc*m1*(-2D0*(pap1 + pap2) + m4**2*(2D0 + (m4**2 + pap1 - 1D0*pap2)*xv) + k1k2*(2D0 + m4**2*xv - 2D0*pap2*xv)))) + 
     &          A3Rc*(-1D0*B1L*(4D0*D3Rc*(k1k2 - 1D0*pap1)*(pap1 + pap2) + 
     &                D3Rc*(-2D0*k1k2**2*m1**2 + m1**2*m3**2*m4**2 - 1D0*k1k2*m4**2*(m1**2 + 2D0*pap1) + 4D0*k1k2*pap1*pap2 - 1D0*(m4**2 - 2D0*pap2)*(m4**2*pap1 - 1D0*m3**2*pap2))*
     &                 xv + D3Lc*m1*mq3c*(2D0*(k1k2 + m3**2) + (2D0*k1k2**2 + (k1k2 - 1D0*m3**2)*m4**2)*xv)) + 
     &             B1R*m3*mq1*(D3Lc*m1*mq3c*(4D0 - 1D0*m4**2*xv) + D3Rc*(m1**2*(-2D0 + m4**2*xv) + pap2*(2D0 + m4**2*xv - 2D0*pap2*xv)))))))/
     &   (kbk3**2*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark)*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)))


        radvSiUf = dreal((GS**2*(A1R*(A3Rc*(B1R*m3*mq1*(2D0*D3Rc*(-1D0*k2k3*kakb - 1D0*(m1**2 - 1D0*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + kak3*pbp2) + 
     &                D3Rc*(k2k3*kakb*(m4**2 - 2D0*pap2) + (m1**2*m4**2 + (m4**2 - 2D0*pap2)*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*kak3*(m4**2 - 2D0*pap2)*pbp2)*xv + 
     &                D3Lc*m1*mq3c*(k1k3 - 1D0*kbk3 - 2D0*pbp1)*(4D0 - 1D0*m4**2*xv)) + 
     &             B1L*(4D0*D3Rc*(pap1 + pap2)*((k2k3 - 1D0*kak3)*pbp1 - 1D0*pap1*(kbk3 + 2D0*pbp1) + k1k2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) + k1k3*(kakb + pap1 - 1D0*pbp2)) + 
     &                D3Rc*(-1D0*k1k3*m1**2*m3**2*m4**2 + kbk3*m1**2*m3**2*m4**2 + k1k3*kakb*m4**4 + k1k3*m4**4*pap1 - 1D0*kbk3*m4**4*pap1 - 
     &                   1D0*k2k3*kakb*m3**2*(m4**2 - 2D0*pap2) - 2D0*k1k3*kakb*m4**2*pap2 - 1D0*k1k3*m3**2*m4**2*pap2 + kbk3*m3**2*m4**2*pap2 - 2D0*k1k3*m4**2*pap1*pap2 + 
     &                   2D0*kbk3*m4**2*pap1*pap2 + 2D0*k1k3*m3**2*pap2**2 - 2D0*kbk3*m3**2*pap2**2 + 2D0*k1k2**2*m1**2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*k2k3*m1**2*m4**2*pbp1 + 
     &                   2D0*m1**2*m3**2*m4**2*pbp1 - 1D0*kak3*m4**4*pbp1 - 2D0*m4**4*pap1*pbp1 + 2D0*kak3*m4**2*pap2*pbp1 + 2D0*m3**2*m4**2*pap2*pbp1 + 4D0*m4**2*pap1*pap2*pbp1 - 
     &                   4D0*m3**2*pap2**2*pbp1 + k1k3*m1**2*m4**2*pbp2 + kak3*m3**2*m4**2*pbp2 - 2D0*kak3*m3**2*pap2*pbp2 + 
     &                   k1k2*(-1D0*kbk3*m4**2*(m1**2 + 2D0*pap1) + 4D0*kbk3*pap1*pap2 - 2D0*(k2k3*m1**2 + m4**2*(kak3 + m1**2 + 2D0*pap1) - 2D0*(kak3 + 2D0*pap1)*pap2)*pbp1 + 
     &                      k1k3*(2D0*kakb*(m4**2 - 2D0*pap2) + 2D0*pap1*(m4**2 - 2D0*pap2) + m1**2*(m4**2 + 2D0*pbp2))))*xv + 
     &                D3Lc*m1*mq3c*(2D0*(k2k3*pbp1 + m3**2*(kbk3 + 2D0*pbp1) + k1k2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) - 1D0*k1k3*(m3**2 + pbp2)) + 
     &                   (2D0*k1k2**2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) + k1k2*(2D0*k2k3*pbp1 + m4**2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) - 2D0*k1k3*pbp2) + 
     &                      m4**2*(m3**2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + k2k3*pbp1 - 1D0*k1k3*pbp2))*xv))) + 
     &          A3Lc*(B1R*mq1*(D3Lc*m1*(kbk3*pap1 + k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + (-1D0*k2k3 + kak3 + 2D0*pap1)*pbp1 - 1D0*k1k3*(kakb + pap1 - 1D0*pbp2))*(-4D0 + m4**2*xv) + 
     &                D3Rc*mq3c*(-2D0*k1k3*(kakb + pap1) + 2D0*kak3*pbp1 + 2D0*pap1*(kbk3 + 2D0*pbp1) + 
     &                   (kbk3*(-1D0*m4**2*pap1 + 2D0*k1k2*pap2) + (-1D0*m4**2*(kak3 + 2D0*pap1) + 2D0*(2D0*k1k2 + k2k3)*pap2)*pbp1 + 
     &                      k1k3*(m4**2*(kakb + pap1) - 2D0*pap2*(k1k2 + pbp2)))*xv)) + 
     &             B1L*m3*(D3Rc*mq3c*(4D0*(pap1 + pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + 2D0*k1k2*k2k3*kakb*xv + 
     &                   (-1D0*k1k3*m4**2*(kakb + pap1 + pap2) + k2k3*(kakb*m4**2 - 2D0*pap2*pbp1) + m4**2*(kak3*pbp1 + (pap1 + pap2)*(kbk3 + 2D0*pbp1)) - 
     &                      1D0*kak3*(2D0*k1k2 + m4**2)*pbp2 + 2D0*k1k3*pap2*pbp2)*xv) + 
     &                D3Lc*m1*(2D0*(-1D0*k2k3*kakb + kbk3*(-1D0*m4**2 + pap1 + pap2) + k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*k2k3*pbp1 + kak3*pbp1 - 2D0*m4**2*pbp1 + 
     &                      2D0*pap1*pbp1 + 2D0*pap2*pbp1 - 1D0*k1k3*(kakb - 1D0*m4**2 + pap1 + pap2 - 1D0*pbp2) + kak3*pbp2) + 
     &                   (k1k2*(-2D0*k2k3*kakb + (m4**2 - 2D0*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + 2D0*kak3*pbp2) + 
     &                      m4**2*(k2k3*(-1D0*kakb + pbp1) - 1D0*(m4**2 + pap1 - 1D0*pap2)*(kbk3 + 2D0*pbp1) + k1k3*(kakb + m4**2 + pap1 - 1D0*pap2 - 1D0*pbp2) + 
     &                         kak3*(-1D0*pbp1 + pbp2)))*xv)))) + A1L*(A3Lc*
     &           (B1L*m3*mq1*(2D0*D3Lc*(-1D0*k2k3*kakb - 1D0*(m1**2 - 1D0*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + kak3*pbp2) + 
     &                D3Lc*(k2k3*kakb*(m4**2 - 2D0*pap2) + (m1**2*m4**2 + (m4**2 - 2D0*pap2)*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*kak3*(m4**2 - 2D0*pap2)*pbp2)*xv + 
     &                D3Rc*m1*mq3c*(k1k3 - 1D0*kbk3 - 2D0*pbp1)*(4D0 - 1D0*m4**2*xv)) + 
     &             B1R*(4D0*D3Lc*(pap1 + pap2)*((k2k3 - 1D0*kak3)*pbp1 - 1D0*pap1*(kbk3 + 2D0*pbp1) + k1k2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) + k1k3*(kakb + pap1 - 1D0*pbp2)) + 
     &                D3Lc*(-1D0*k1k3*m1**2*m3**2*m4**2 + kbk3*m1**2*m3**2*m4**2 + k1k3*kakb*m4**4 + k1k3*m4**4*pap1 - 1D0*kbk3*m4**4*pap1 - 
     &                   1D0*k2k3*kakb*m3**2*(m4**2 - 2D0*pap2) - 2D0*k1k3*kakb*m4**2*pap2 - 1D0*k1k3*m3**2*m4**2*pap2 + kbk3*m3**2*m4**2*pap2 - 2D0*k1k3*m4**2*pap1*pap2 + 
     &                   2D0*kbk3*m4**2*pap1*pap2 + 2D0*k1k3*m3**2*pap2**2 - 2D0*kbk3*m3**2*pap2**2 + 2D0*k1k2**2*m1**2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*k2k3*m1**2*m4**2*pbp1 + 
     &                   2D0*m1**2*m3**2*m4**2*pbp1 - 1D0*kak3*m4**4*pbp1 - 2D0*m4**4*pap1*pbp1 + 2D0*kak3*m4**2*pap2*pbp1 + 2D0*m3**2*m4**2*pap2*pbp1 + 4D0*m4**2*pap1*pap2*pbp1 - 
     &                   4D0*m3**2*pap2**2*pbp1 + k1k3*m1**2*m4**2*pbp2 + kak3*m3**2*m4**2*pbp2 - 2D0*kak3*m3**2*pap2*pbp2 + 
     &                   k1k2*(-1D0*kbk3*m4**2*(m1**2 + 2D0*pap1) + 4D0*kbk3*pap1*pap2 - 2D0*(k2k3*m1**2 + m4**2*(kak3 + m1**2 + 2D0*pap1) - 2D0*(kak3 + 2D0*pap1)*pap2)*pbp1 + 
     &                      k1k3*(2D0*kakb*(m4**2 - 2D0*pap2) + 2D0*pap1*(m4**2 - 2D0*pap2) + m1**2*(m4**2 + 2D0*pbp2))))*xv + 
     &                D3Rc*m1*mq3c*(2D0*(k2k3*pbp1 + m3**2*(kbk3 + 2D0*pbp1) + k1k2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) - 1D0*k1k3*(m3**2 + pbp2)) + 
     &                   (2D0*k1k2**2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) + k1k2*(2D0*k2k3*pbp1 + m4**2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) - 2D0*k1k3*pbp2) + 
     &                      m4**2*(m3**2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + k2k3*pbp1 - 1D0*k1k3*pbp2))*xv))) + 
     &          A3Rc*(B1L*mq1*(D3Rc*m1*(kbk3*pap1 + k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + (-1D0*k2k3 + kak3 + 2D0*pap1)*pbp1 - 1D0*k1k3*(kakb + pap1 - 1D0*pbp2))*(-4D0 + m4**2*xv) + 
     &                D3Lc*mq3c*(-2D0*k1k3*(kakb + pap1) + 2D0*kak3*pbp1 + 2D0*pap1*(kbk3 + 2D0*pbp1) + 
     &                   (kbk3*(-1D0*m4**2*pap1 + 2D0*k1k2*pap2) + (-1D0*m4**2*(kak3 + 2D0*pap1) + 2D0*(2D0*k1k2 + k2k3)*pap2)*pbp1 + 
     &                      k1k3*(m4**2*(kakb + pap1) - 2D0*pap2*(k1k2 + pbp2)))*xv)) + 
     &             B1R*m3*(D3Lc*mq3c*(4D0*(pap1 + pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + 2D0*k1k2*k2k3*kakb*xv + 
     &                   (-1D0*k1k3*m4**2*(kakb + pap1 + pap2) + k2k3*(kakb*m4**2 - 2D0*pap2*pbp1) + m4**2*(kak3*pbp1 + (pap1 + pap2)*(kbk3 + 2D0*pbp1)) - 
     &                      1D0*kak3*(2D0*k1k2 + m4**2)*pbp2 + 2D0*k1k3*pap2*pbp2)*xv) + 
     &                D3Rc*m1*(2D0*(-1D0*k2k3*kakb + kbk3*(-1D0*m4**2 + pap1 + pap2) + k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*k2k3*pbp1 + kak3*pbp1 - 2D0*m4**2*pbp1 + 
     &                      2D0*pap1*pbp1 + 2D0*pap2*pbp1 - 1D0*k1k3*(kakb - 1D0*m4**2 + pap1 + pap2 - 1D0*pbp2) + kak3*pbp2) + 
     &                   (k1k2*(-2D0*k2k3*kakb + (m4**2 - 2D0*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + 2D0*kak3*pbp2) + 
     &                      m4**2*(k2k3*(-1D0*kakb + pbp1) - 1D0*(m4**2 + pap1 - 1D0*pap2)*(kbk3 + 2D0*pbp1) + k1k3*(kakb + m4**2 + pap1 - 1D0*pap2 - 1D0*pbp2) + 
     &                         kak3*(-1D0*pbp1 + pbp2)))*xv))))))/(k1k3*kbk3*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark)*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)))


        radvSpUi = dreal((-1D0*GS**2*(A1L*(A3Lc*(B1R*(D3Lc*(4D0*k1k2*(k1k3*(kakb + m1**2) + k2k3*(kakb + m1**2) + (kbk3 - 2D0*m2**2)*(pap1 + pap2) + kak3*(mq1**2 - 1D0*pbp1 - 1D0*pbp2) - 
     &                      2D0*(kakb*mq1**2 + m1**2*(pbp1 + pbp2))) - 4D0*pap1*
     &                    (k1k3*(kakb + m1**2) + k2k3*(kakb + m1**2) + (kbk3 - 2D0*m2**2)*(pap1 + pap2) + kak3*(mq1**2 - 1D0*pbp1 - 1D0*pbp2) - 
     &                      2D0*(kakb*mq1**2 + m1**2*(pbp1 + pbp2))) - 2D0*k1k2**2*m1**2*(kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*xv + 
     &                   (kak3*m1**2*m3**2*m4**2 - 2D0*kakb*m1**2*m3**2*m4**2 + kbk3*m1**2*m3**2*m4**2 - 2D0*m1**2*m2**2*m3**2*m4**2 - 1D0*k1k3*kakb*m4**4 - 1D0*k1k3*m1**2*m4**4 + 
     &                      k1k3*m1**2*m4**2*mq1**2 - 1D0*kbk3*m4**4*pap1 + 2D0*m2**2*m4**4*pap1 + 2D0*k1k3*kakb*m4**2*pap2 + k1k3*m1**2*m4**2*pap2 + kbk3*m3**2*m4**2*pap2 - 
     &                      2D0*m2**2*m3**2*m4**2*pap2 + k1k3*m4**2*mq1**2*pap2 + 2D0*kbk3*m4**2*pap1*pap2 - 4D0*m2**2*m4**2*pap1*pap2 - 2D0*kbk3*m3**2*pap2**2 + 
     &                      4D0*m2**2*m3**2*pap2**2 - 2D0*k1k3*mq1**2*pap2**2 + kak3*m4**4*pbp1 + 2D0*m1**2*m4**4*pbp1 - 2D0*m1**2*m4**2*mq1**2*pbp1 - 2D0*kak3*m4**2*pap2*pbp1 - 
     &                      2D0*m1**2*m4**2*pap2*pbp1 - 2D0*m4**2*mq1**2*pap2*pbp1 + 4D0*mq1**2*pap2**2*pbp1 + 
     &                      k2k3*(kakb*m3**2*(m4**2 - 2D0*pap2) - 1D0*mq1**2*pap1*(m4**2 - 2D0*pap2) + m1**2*(m3**2*(m4**2 - 2D0*pap2) + m4**2*(pap1 + pbp1))) - 
     &                      1D0*(m4**2*(k1k3*m1**2 + (kak3 + 2D0*m1**2)*m3**2 + 2D0*(m1 - 1D0*mq1)*(m1 + mq1)*pap1) - 2D0*((kak3 + 2D0*m1**2)*m3**2 - 2D0*mq1**2*pap1)*pap2)*pbp2)*
     &                    xv - 1D0*k1k2*(-2D0*kakb*m1**2*m4**2 + kbk3*m1**2*m4**2 - 2D0*m1**2*m2**2*m4**2 + 2D0*k2k3*m1**2*mq1**2 - 2D0*kakb*m4**2*mq1**2 - 2D0*k2k3*m1**2*pap1 + 
     &                      2D0*kbk3*m4**2*pap1 - 4D0*m2**2*m4**2*pap1 + 4D0*kakb*mq1**2*pap2 - 4D0*kbk3*pap1*pap2 + 8D0*m2**2*pap1*pap2 + 
     &                      kak3*(m1**2*m4**2 + (m4**2 - 2D0*pap2)*(mq1**2 - 2D0*pbp1)) - 2D0*k2k3*m1**2*pbp1 - 4D0*m1**2*m4**2*pbp1 + 4D0*m1**2*pap2*pbp1 + 
     &                      4D0*m1**2*(-1D0*mq1**2 + pap1)*pbp2 + 2D0*k1k3*(kakb*(m4**2 - 2D0*pap2) + m1**2*(m4**2 - 1D0*pap2 + pbp2)))*xv) + 
     &                D3Rc*m1*mq3c*(2D0*k1k2*(kak3 - 2D0*kakb + kbk3 - 2D0*m2**2) + 2D0*k1k2**2*(kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*xv + 
     &                   k1k2*((kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*m4**2 + 
     &                      2D0*(pap2*(k1k3 - 2D0*pbp1) + k2k3*(mq1**2 - 1D0*pap1 - 1D0*pbp1) + (k1k3 - 2D0*mq1**2 + 2D0*pap1)*pbp2))*xv - 
     &                   1D0*((kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*m3**2 + (mq1**2 - 1D0*pap2)*(k1k3 - 2D0*pbp1) + k2k3*(pap1 + pbp1) - 1D0*(k1k3 + 2D0*pap1)*pbp2)*
     &                    (-2D0 + m4**2*xv))) - 1D0*B1L*m3*mq1*(-1D0*D3Rc*m1*mq3c*(k1k3 + k2k3 + kak3 - 2D0*kakb + kbk3 - 2D0*(m2**2 + pbp1 + pbp2))*(-4D0 + m4**2*xv) + 
     &                D3Lc*(-2D0*(2D0*kakb*(-1D0*m1**2 + m4**2) + k2k3*(-1D0*kakb + pap1) + k1k3*(m1**2 - 1D0*pap2) + (m1**2 - 1D0*pap2)*(kbk3 - 2D0*(m2**2 + pbp1)) - 
     &                      2D0*pap1*pbp2 + kak3*(m1**2 - 1D0*m4**2 + pbp2)) + 
     &                   (kak3*m1**2*m4**2 - 2D0*kakb*m1**2*m4**2 + kbk3*m1**2*m4**2 - 2D0*m1**2*m2**2*m4**2 + kak3*m4**4 - 2D0*kakb*m4**4 + 
     &                      k2k3*m4**2*(kakb + 2D0*m1**2 - 1D0*pap1) - 2D0*kak3*m4**2*pap2 + 4D0*kakb*m4**2*pap2 + kbk3*m4**2*pap2 - 2D0*m2**2*m4**2*pap2 - 
     &                      2D0*k2k3*(kakb + m1**2 - 1D0*pap1)*pap2 - 2D0*kbk3*pap2**2 + 4D0*m2**2*pap2**2 + k1k3*(m1**2*m4**2 + (m4**2 - 2D0*pap2)*pap2) - 2D0*m1**2*m4**2*pbp1 - 
     &                      2D0*m4**2*pap2*pbp1 + 4D0*pap2**2*pbp1 + (-1D0*m4**2*(kak3 + 4D0*m1**2 - 2D0*pap1) + 2D0*(kak3 + 2D0*m1**2 - 2D0*pap1)*pap2)*pbp2)*xv + 
     &                   k1k2*(kak3 - 2D0*kakb)*(2D0 + m4**2*xv - 2D0*pap2*xv)))) - 
     &          1D0*A3Rc*(B1R*m3*(-1D0*D3Lc*mq3c*(k1k3*(kakb + m1**2) + k2k3*(kakb + m1**2) + (kbk3 - 2D0*m2**2)*(pap1 + pap2) + kak3*(mq1**2 - 1D0*pbp1 - 1D0*pbp2) - 
     &                   2D0*(kakb*mq1**2 + m1**2*(pbp1 + pbp2)))*(-4D0 + m4**2*xv) + 
     &                D3Rc*m1*(2D0*kak3*m4**2 - 2D0*(2D0*kakb*m4**2 - 1D0*kbk3*m4**2 + 2D0*m2**2*m4**2 + kak3*mq1**2 - 2D0*kakb*mq1**2 + kbk3*pap1 - 2D0*m2**2*pap1 + 
     &                      k1k3*(kakb + m1**2 - 1D0*pap2) + kbk3*pap2 - 2D0*m2**2*pap2 - 1D0*(kak3 + 2D0*m1**2 - 2D0*pap2)*pbp1 + k2k3*(kakb + m1**2 - 1D0*mq1**2 + pap1 + pbp1))
     &                    + 2D0*(k1k3 + kak3 + 2D0*(m1**2 - 1D0*mq1**2 + pap1))*pbp2 + 
     &                   (kak3*m4**4 - 2D0*kakb*m4**4 + kbk3*m4**4 - 2D0*m2**2*m4**4 + kak3*m4**2*mq1**2 - 2D0*kakb*m4**2*mq1**2 + kbk3*m4**2*pap1 - 2D0*m2**2*m4**2*pap1 - 
     &                      2D0*kak3*m4**2*pap2 + 4D0*kakb*m4**2*pap2 - 1D0*kbk3*m4**2*pap2 + 2D0*m2**2*m4**2*pap2 + 
     &                      k2k3*(kakb*m4**2 + m1**2*m4**2 + (m4**2 - 2D0*pap2)*(mq1**2 - 1D0*pap1 - 1D0*pbp1)) - 1D0*kak3*m4**2*pbp1 - 2D0*m1**2*m4**2*pbp1 - 
     &                      2D0*m4**2*pap2*pbp1 + 4D0*pap2**2*pbp1 - 1D0*(m4**2*(kak3 + 2D0*(m1**2 + mq1**2 - 1D0*pap1)) + 4D0*(-1D0*mq1**2 + pap1)*pap2)*pbp2 + 
     &                      k1k3*(kakb*m4**2 + m1**2*m4**2 + (m4**2 - 2D0*pap2)*(pap2 + pbp2)))*xv + k1k2*(kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*(2D0 + m4**2*xv - 2D0*pap2*xv))) + 
     &             B1L*mq1*(D3Rc*m1*(4D0*(-1D0*k1k2 + pap1)*(k1k3 + k2k3 + kak3 - 2D0*kakb + kbk3 - 2D0*(m2**2 + pbp1 + pbp2)) + 
     &                   (2D0*k1k2**2*(kak3 - 2D0*kakb) - 1D0*k2k3*m3**2*m4**2 - 1D0*kak3*m3**2*m4**2 + 2D0*kakb*m3**2*m4**2 - 1D0*kbk3*m4**2*pap1 + 2D0*m2**2*m4**2*pap1 + 
     &                      2D0*k2k3*m3**2*pap2 - 2D0*k2k3*pap1*pap2 + k2k3*m4**2*pbp1 + kak3*m4**2*pbp1 + 2D0*m1**2*m4**2*pbp1 - 2D0*m4**4*pbp1 - 2D0*k2k3*pap2*pbp1 + 
     &                      4D0*m4**2*pap2*pbp1 - 4D0*pap2**2*pbp1 + 2D0*m3**2*m4**2*pbp2 - 4D0*m3**2*pap2*pbp2 + 4D0*pap1*pap2*pbp2 + 
     &                      k1k2*((kbk3 - 2D0*m2**2)*m4**2 + 2D0*k2k3*(kakb + m1**2 - 1D0*pap1) + 2D0*(k1k3 + kak3)*(m4**2 - 1D0*pap2) + 4D0*kakb*(-1D0*m4**2 + pap2) + 
     &                         4D0*(-1D0*m4**2 + pap2)*pbp1 - 2D0*kak3*pbp2 + 4D0*(-1D0*m1**2 + pap1)*pbp2) - 
     &                      1D0*k1k3*(kakb*m4**2 + m1**2*m4**2 - 1D0*m4**4 + 2D0*m4**2*pap2 - 2D0*pap2**2 + m4**2*pbp2 - 2D0*pap2*pbp2))*xv) + 
     &                D3Lc*mq3c*(-2D0*k1k2**2*(kak3 - 2D0*kakb)*xv + ((kak3 - 2D0*kakb)*m3**2 + k1k3*(kakb + m1**2 - 1D0*pap2) - 1D0*(kak3 + 2D0*m1**2 - 2D0*pap2)*pbp1 + 
     &                      pap1*(k2k3 + kbk3 - 2D0*(m2**2 + pbp2)))*(-2D0 + m4**2*xv) + 
     &                   k1k2*(-2D0*(pap2*(k1k3 + kbk3 - 2D0*(m2**2 + pbp1)) + (m1**2 - 1D0*pap1)*(k2k3 - 2D0*pbp2))*xv + kakb*(4D0 - 2D0*k2k3*xv + 2D0*m4**2*xv) - 
     &                      1D0*kak3*(2D0 + m4**2*xv - 2D0*pbp2*xv)))))) + 
     &       A1R*(A3Rc*(B1L*(D3Rc*(4D0*k1k2*(k1k3*(kakb + m1**2) + k2k3*(kakb + m1**2) + (kbk3 - 2D0*m2**2)*(pap1 + pap2) + kak3*(mq1**2 - 1D0*pbp1 - 1D0*pbp2) - 
     &                      2D0*(kakb*mq1**2 + m1**2*(pbp1 + pbp2))) - 4D0*pap1*
     &                    (k1k3*(kakb + m1**2) + k2k3*(kakb + m1**2) + (kbk3 - 2D0*m2**2)*(pap1 + pap2) + kak3*(mq1**2 - 1D0*pbp1 - 1D0*pbp2) - 
     &                      2D0*(kakb*mq1**2 + m1**2*(pbp1 + pbp2))) - 2D0*k1k2**2*m1**2*(kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*xv + 
     &                   (kak3*m1**2*m3**2*m4**2 - 2D0*kakb*m1**2*m3**2*m4**2 + kbk3*m1**2*m3**2*m4**2 - 2D0*m1**2*m2**2*m3**2*m4**2 - 1D0*k1k3*kakb*m4**4 - 1D0*k1k3*m1**2*m4**4 + 
     &                      k1k3*m1**2*m4**2*mq1**2 - 1D0*kbk3*m4**4*pap1 + 2D0*m2**2*m4**4*pap1 + 2D0*k1k3*kakb*m4**2*pap2 + k1k3*m1**2*m4**2*pap2 + kbk3*m3**2*m4**2*pap2 - 
     &                      2D0*m2**2*m3**2*m4**2*pap2 + k1k3*m4**2*mq1**2*pap2 + 2D0*kbk3*m4**2*pap1*pap2 - 4D0*m2**2*m4**2*pap1*pap2 - 2D0*kbk3*m3**2*pap2**2 + 
     &                      4D0*m2**2*m3**2*pap2**2 - 2D0*k1k3*mq1**2*pap2**2 + kak3*m4**4*pbp1 + 2D0*m1**2*m4**4*pbp1 - 2D0*m1**2*m4**2*mq1**2*pbp1 - 2D0*kak3*m4**2*pap2*pbp1 - 
     &                      2D0*m1**2*m4**2*pap2*pbp1 - 2D0*m4**2*mq1**2*pap2*pbp1 + 4D0*mq1**2*pap2**2*pbp1 + 
     &                      k2k3*(kakb*m3**2*(m4**2 - 2D0*pap2) - 1D0*mq1**2*pap1*(m4**2 - 2D0*pap2) + m1**2*(m3**2*(m4**2 - 2D0*pap2) + m4**2*(pap1 + pbp1))) - 
     &                      1D0*(m4**2*(k1k3*m1**2 + (kak3 + 2D0*m1**2)*m3**2 + 2D0*(m1 - 1D0*mq1)*(m1 + mq1)*pap1) - 2D0*((kak3 + 2D0*m1**2)*m3**2 - 2D0*mq1**2*pap1)*pap2)*pbp2)*
     &                    xv - 1D0*k1k2*(-2D0*kakb*m1**2*m4**2 + kbk3*m1**2*m4**2 - 2D0*m1**2*m2**2*m4**2 + 2D0*k2k3*m1**2*mq1**2 - 2D0*kakb*m4**2*mq1**2 - 2D0*k2k3*m1**2*pap1 + 
     &                      2D0*kbk3*m4**2*pap1 - 4D0*m2**2*m4**2*pap1 + 4D0*kakb*mq1**2*pap2 - 4D0*kbk3*pap1*pap2 + 8D0*m2**2*pap1*pap2 + 
     &                      kak3*(m1**2*m4**2 + (m4**2 - 2D0*pap2)*(mq1**2 - 2D0*pbp1)) - 2D0*k2k3*m1**2*pbp1 - 4D0*m1**2*m4**2*pbp1 + 4D0*m1**2*pap2*pbp1 + 
     &                      4D0*m1**2*(-1D0*mq1**2 + pap1)*pbp2 + 2D0*k1k3*(kakb*(m4**2 - 2D0*pap2) + m1**2*(m4**2 - 1D0*pap2 + pbp2)))*xv) + 
     &                D3Lc*m1*mq3c*(2D0*k1k2*(kak3 - 2D0*kakb + kbk3 - 2D0*m2**2) + 2D0*k1k2**2*(kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*xv + 
     &                   k1k2*((kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*m4**2 + 
     &                      2D0*(pap2*(k1k3 - 2D0*pbp1) + k2k3*(mq1**2 - 1D0*pap1 - 1D0*pbp1) + (k1k3 - 2D0*mq1**2 + 2D0*pap1)*pbp2))*xv - 
     &                   1D0*((kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*m3**2 + (mq1**2 - 1D0*pap2)*(k1k3 - 2D0*pbp1) + k2k3*(pap1 + pbp1) - 1D0*(k1k3 + 2D0*pap1)*pbp2)*
     &                    (-2D0 + m4**2*xv))) - 1D0*B1R*m3*mq1*(-1D0*D3Lc*m1*mq3c*(k1k3 + k2k3 + kak3 - 2D0*kakb + kbk3 - 2D0*(m2**2 + pbp1 + pbp2))*(-4D0 + m4**2*xv) + 
     &                D3Rc*(-2D0*(2D0*kakb*(-1D0*m1**2 + m4**2) + k2k3*(-1D0*kakb + pap1) + k1k3*(m1**2 - 1D0*pap2) + (m1**2 - 1D0*pap2)*(kbk3 - 2D0*(m2**2 + pbp1)) - 
     &                      2D0*pap1*pbp2 + kak3*(m1**2 - 1D0*m4**2 + pbp2)) + 
     &                   (kak3*m1**2*m4**2 - 2D0*kakb*m1**2*m4**2 + kbk3*m1**2*m4**2 - 2D0*m1**2*m2**2*m4**2 + kak3*m4**4 - 2D0*kakb*m4**4 + 
     &                      k2k3*m4**2*(kakb + 2D0*m1**2 - 1D0*pap1) - 2D0*kak3*m4**2*pap2 + 4D0*kakb*m4**2*pap2 + kbk3*m4**2*pap2 - 2D0*m2**2*m4**2*pap2 - 
     &                      2D0*k2k3*(kakb + m1**2 - 1D0*pap1)*pap2 - 2D0*kbk3*pap2**2 + 4D0*m2**2*pap2**2 + k1k3*(m1**2*m4**2 + (m4**2 - 2D0*pap2)*pap2) - 2D0*m1**2*m4**2*pbp1 - 
     &                      2D0*m4**2*pap2*pbp1 + 4D0*pap2**2*pbp1 + (-1D0*m4**2*(kak3 + 4D0*m1**2 - 2D0*pap1) + 2D0*(kak3 + 2D0*m1**2 - 2D0*pap1)*pap2)*pbp2)*xv + 
     &                   k1k2*(kak3 - 2D0*kakb)*(2D0 + m4**2*xv - 2D0*pap2*xv)))) - 
     &          1D0*A3Lc*(B1L*m3*(-1D0*D3Rc*mq3c*(k1k3*(kakb + m1**2) + k2k3*(kakb + m1**2) + (kbk3 - 2D0*m2**2)*(pap1 + pap2) + kak3*(mq1**2 - 1D0*pbp1 - 1D0*pbp2) - 
     &                   2D0*(kakb*mq1**2 + m1**2*(pbp1 + pbp2)))*(-4D0 + m4**2*xv) + 
     &                D3Lc*m1*(2D0*kak3*m4**2 - 2D0*(2D0*kakb*m4**2 - 1D0*kbk3*m4**2 + 2D0*m2**2*m4**2 + kak3*mq1**2 - 2D0*kakb*mq1**2 + kbk3*pap1 - 2D0*m2**2*pap1 + 
     &                      k1k3*(kakb + m1**2 - 1D0*pap2) + kbk3*pap2 - 2D0*m2**2*pap2 - 1D0*(kak3 + 2D0*m1**2 - 2D0*pap2)*pbp1 + k2k3*(kakb + m1**2 - 1D0*mq1**2 + pap1 + pbp1))
     &                    + 2D0*(k1k3 + kak3 + 2D0*(m1**2 - 1D0*mq1**2 + pap1))*pbp2 + 
     &                   (kak3*m4**4 - 2D0*kakb*m4**4 + kbk3*m4**4 - 2D0*m2**2*m4**4 + kak3*m4**2*mq1**2 - 2D0*kakb*m4**2*mq1**2 + kbk3*m4**2*pap1 - 2D0*m2**2*m4**2*pap1 - 
     &                      2D0*kak3*m4**2*pap2 + 4D0*kakb*m4**2*pap2 - 1D0*kbk3*m4**2*pap2 + 2D0*m2**2*m4**2*pap2 + 
     &                      k2k3*(kakb*m4**2 + m1**2*m4**2 + (m4**2 - 2D0*pap2)*(mq1**2 - 1D0*pap1 - 1D0*pbp1)) - 1D0*kak3*m4**2*pbp1 - 2D0*m1**2*m4**2*pbp1 - 
     &                      2D0*m4**2*pap2*pbp1 + 4D0*pap2**2*pbp1 - 1D0*(m4**2*(kak3 + 2D0*(m1**2 + mq1**2 - 1D0*pap1)) + 4D0*(-1D0*mq1**2 + pap1)*pap2)*pbp2 + 
     &                      k1k3*(kakb*m4**2 + m1**2*m4**2 + (m4**2 - 2D0*pap2)*(pap2 + pbp2)))*xv + k1k2*(kak3 - 2D0*kakb + kbk3 - 2D0*m2**2)*(2D0 + m4**2*xv - 2D0*pap2*xv))) + 
     &             B1R*mq1*(D3Lc*m1*(4D0*(-1D0*k1k2 + pap1)*(k1k3 + k2k3 + kak3 - 2D0*kakb + kbk3 - 2D0*(m2**2 + pbp1 + pbp2)) + 
     &                   (2D0*k1k2**2*(kak3 - 2D0*kakb) - 1D0*k2k3*m3**2*m4**2 - 1D0*kak3*m3**2*m4**2 + 2D0*kakb*m3**2*m4**2 - 1D0*kbk3*m4**2*pap1 + 2D0*m2**2*m4**2*pap1 + 
     &                      2D0*k2k3*m3**2*pap2 - 2D0*k2k3*pap1*pap2 + k2k3*m4**2*pbp1 + kak3*m4**2*pbp1 + 2D0*m1**2*m4**2*pbp1 - 2D0*m4**4*pbp1 - 2D0*k2k3*pap2*pbp1 + 
     &                      4D0*m4**2*pap2*pbp1 - 4D0*pap2**2*pbp1 + 2D0*m3**2*m4**2*pbp2 - 4D0*m3**2*pap2*pbp2 + 4D0*pap1*pap2*pbp2 + 
     &                      k1k2*((kbk3 - 2D0*m2**2)*m4**2 + 2D0*k2k3*(kakb + m1**2 - 1D0*pap1) + 2D0*(k1k3 + kak3)*(m4**2 - 1D0*pap2) + 4D0*kakb*(-1D0*m4**2 + pap2) + 
     &                         4D0*(-1D0*m4**2 + pap2)*pbp1 - 2D0*kak3*pbp2 + 4D0*(-1D0*m1**2 + pap1)*pbp2) - 
     &                      1D0*k1k3*(kakb*m4**2 + m1**2*m4**2 - 1D0*m4**4 + 2D0*m4**2*pap2 - 2D0*pap2**2 + m4**2*pbp2 - 2D0*pap2*pbp2))*xv) + 
     &                D3Rc*mq3c*(-2D0*k1k2**2*(kak3 - 2D0*kakb)*xv + ((kak3 - 2D0*kakb)*m3**2 + k1k3*(kakb + m1**2 - 1D0*pap2) - 1D0*(kak3 + 2D0*m1**2 - 2D0*pap2)*pbp1 + 
     &                      pap1*(k2k3 + kbk3 - 2D0*(m2**2 + pbp2)))*(-2D0 + m4**2*xv) + 
     &                   k1k2*(-2D0*(pap2*(k1k3 + kbk3 - 2D0*(m2**2 + pbp1)) + (m1**2 - 1D0*pap1)*(k2k3 - 2D0*pbp2))*xv + kakb*(4D0 - 2D0*k2k3*xv + 2D0*m4**2*xv) - 
     &                      1D0*kak3*(2D0 + m4**2*xv - 2D0*pbp2*xv))))))))/
     &   (kbk3*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2 - cI*mq1*wquark)*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark)*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)))


        radvSpUf = dreal((2D0*GS**2*(A1L*(A3Rc*(B1R*m3*(D3Rc*m1*(-2D0*(k1k2*(kak3 + kakb + kbk3 + m1**2) + (kakb + m1**2)*(k1k3 + k2k3 + m3**2) + (kak3 + kbk3)*m4**2) + 
     &                   2D0*(k1k2 + 2D0*k2k3 - 2D0*kak3)*mq1**2 + 2D0*(-1D0*k2k3 + 2D0*kak3 + kbk3 + m4**2 - 1D0*mq1**2)*pap1 + 2D0*(k1k3 + 2D0*kak3 + kbk3 + m3**2)*pap2 + 
     &                   2D0*(-1D0*k2k3 + kak3 + m4**2 - 1D0*pap2)*pbp1 + 2D0*(k1k3 + kak3 + m3**2 + pap1)*pbp2 + 
     &                   (kakb*m3**2*m4**2 + m1**2*m3**2*m4**2 + 2D0*kak3*m4**4 + 2D0*kbk3*m4**4 + kak3*m4**2*mq1**2 + k2k3*m4**2*(2D0*(kakb + m1**2) - 1D0*mq1**2) + m4**4*pap1 + 
     &                      m4**2*mq1**2*pap1 - 4D0*kak3*m4**2*pap2 - 2D0*kbk3*m4**2*pap2 + m3**2*m4**2*pap2 - 2D0*m4**2*pap1*pap2 - 2D0*m3**2*pap2**2 + m4**4*pbp1 - 
     &                      1D0*m4**2*pap2*pbp1 + k1k2*(2D0*k2k3*(kakb + m1**2) + m4**2*(kakb + 2D0*kbk3 + m1**2 + mq1**2) - 2D0*(kbk3 + mq1**2)*pap2 + 
     &                         2D0*kak3*(m4**2 - 2D0*pap2 - 1D0*pbp2)) - 1D0*(m4**2*(2D0*kak3 - 1D0*m3**2 + pap1) + 2D0*m3**2*pap2)*pbp2)*xv) - 
     &                1D0*D3Lc*mq3c*(-2D0*((kakb + m1**2)*(2D0*k1k2 + k1k3 + k2k3 + 2D0*m3**2) + mq1**2*(kak3 + 2D0*pap1) + kbk3*(pap1 + pap2) - 1D0*(kak3 - 2D0*pap2)*pbp1) + 
     &                   2D0*(kak3 + 2D0*pap1)*pbp2 + (2D0*k2k3*(kakb + m1**2)*m4**2 - 2D0*k2k3*mq1**2*pap2 + 
     &                      k1k2*((kakb + m1**2)*(2D0*k2k3 + m4**2) + 2D0*kbk3*pap2 - 2D0*kak3*pbp2) + 
     &                      m4**2*((kakb + m1**2)*m3**2 + mq1**2*(kak3 + pap1) + pap2*(2D0*kbk3 + pbp1) - 1D0*(2D0*kak3 + pap1)*pbp2))*xv)) + 
     &             B1L*mq1*(D3Lc*mq3c*(2D0*k1k3*(kakb + m1**2 - 2D0*pap1) - 2D0*((2D0*(k1k2 + k2k3) - 1D0*kbk3)*pap1 + m3**2*(kakb + m1**2 + pap1 - 1D0*pap2) + kak3*pbp1) + 
     &                   k1k3*(2D0*k1k2 + m4**2)*pap2*xv + (-2D0*k1k2**2*(kak3 + kakb + m1**2) + m4**2*(k2k3*pap1 + m3**2*(kak3 + kakb + m1**2 + pap1)) - 
     &                      1D0*m3**2*(2D0*k2k3 + m4**2)*pap2 - 1D0*k1k2*(kak3*m4**2 + 2D0*k2k3*(kakb + m1**2 - 1D0*pap1) + 2D0*pap2*(kbk3 + m3**2 + pbp1)) + 
     &                      2D0*k1k2*(kak3 + pap1)*pbp2)*xv) + D3Rc*m1*(-2D0*k1k2*(2D0*k2k3 + kbk3 + 2D0*(m3**2 + pbp1)) + 
     &                   2D0*(-1D0*k2k3*(m3**2 + pbp1) + kak3*(m3**2 + 2D0*pap1 + pbp1) + pap1*(kbk3 + 2D0*(m3**2 + pap1 + pbp1)) + k1k3*(-1D0*kakb - 1D0*m1**2 + m4**2 + pbp2)) - 
     &                   1D0*(kak3*m3**2*m4**2 + kakb*m3**2*m4**2 + m1**2*m3**2*m4**2 + k1k3*m4**4 - 1D0*m3**2*m4**4 + m3**2*m4**2*pap1 + k2k3*m4**2*(-1D0*m3**2 + pap1) - 
     &                      1D0*k1k3*m4**2*pap2 + 2D0*m3**2*m4**2*pap2 - 2D0*m3**2*pap2**2 + m3**2*(m4**2 - 2D0*pap2)*pbp2)*xv + 
     &                   k1k2*((2D0*kbk3 + m3**2)*m4**2 - 2D0*kbk3*pap2 + 2D0*k1k3*(-1D0*m4**2 + pap2) + 
     &                      2D0*(k2k3*(kakb + m1**2 - 1D0*pap1) - 2D0*pap1*pap2 - 1D0*pap2*pbp1 + m4**2*(pap1 + pbp1) - 1D0*pap1*pbp2) + kak3*(3D0*m4**2 - 2D0*(2D0*pap2 + pbp2)))*
     &                    xv + 2D0*k1k2**2*(-2D0 + (kak3 + kakb + m1**2)*xv)))) + 
     &          A3Lc*(B1L*m3*mq1*(-2D0*D3Lc*(k1k2*(kak3 - 1D0*kakb) + 2D0*kak3*kakb + k1k3*m1**2 + 2D0*kak3*m1**2 + m1**2*m3**2 + kak3*m4**2 - 2D0*kak3*pap1 + m1**2*pap1 - 
     &                   1D0*m4**2*pap1 + k2k3*(-2D0*kakb - 1D0*m1**2 + pap1) - 1D0*k1k3*pap2 - 2D0*kak3*pap2 - 1D0*m3**2*pap2 + m1**2*pbp1 - 1D0*pap2*pbp1 + pap1*pbp2) + 
     &                D3Lc*(kak3*m1**2*m4**2 + kbk3*m1**2*m4**2 + m1**2*m3**2*m4**2 + 2D0*kak3*m4**4 + m1**2*m4**2*pap1 + m4**4*pap1 + 
     &                   k1k2*((2D0*kak3 + kakb)*(m4**2 - 2D0*pap2) + 2D0*m1**2*(k2k3 + m4**2 - 1D0*pap2)) - 4D0*kak3*m4**2*pap2 + kbk3*m4**2*pap2 + m3**2*m4**2*pap2 - 
     &                   2D0*m4**2*pap1*pap2 - 2D0*kbk3*pap2**2 - 2D0*m3**2*pap2**2 + m1**2*m4**2*pbp1 + m4**2*pap2*pbp1 - 2D0*pap2**2*pbp1 + 
     &                   k2k3*(-1D0*kakb*(m4**2 - 2D0*pap2) + m1**2*(m4**2 - 2D0*pbp2)) - 1D0*(kak3 + pap1)*(m4**2 - 2D0*pap2)*pbp2)*xv + 
     &                D3Rc*m1*mq3c*(2D0*(2D0*k1k2 + k1k3 + k2k3 + kak3 + kbk3 + 2D0*(m3**2 + pap1 + pbp1)) - 
     &                   1D0*(k1k2*(2D0*k2k3 + m4**2) + m4**2*(kak3 + kbk3 + m3**2 + pap1 + pbp1) + 2D0*k2k3*(m4**2 - 1D0*pap2 - 1D0*pbp2))*xv)) - 
     &             1D0*B1R*(D3Rc*m1*mq3c*(4D0*(k1k2 + k2k3)*(pap1 + pbp1) - 2D0*k1k3*(mq1**2 - 2D0*(pap1 + pbp1)) + 
     &                   (2D0*k1k2**2*(kak3 + kbk3 + mq1**2) + k1k2*((kak3 + kbk3)*m4**2 + 2D0*k2k3*(mq1**2 - 1D0*pap1 - 1D0*pbp1)) - 1D0*k2k3*m4**2*(pap1 + pbp1))*xv - 
     &                   1D0*k1k3*(2D0*k1k2 + m4**2)*(pap2 + pbp2)*xv + m3**2*
     &                    (2D0*(mq1**2 + pap1 - 1D0*pap2 + pbp1 - 1D0*pbp2) - 
     &                      1D0*(m4**2*(kak3 + kbk3 + mq1**2 + pap1 - 1D0*pap2 + pbp1 - 1D0*pbp2) - 2D0*(k1k2 + k2k3)*(pap2 + pbp2))*xv)) + 
     &                D3Lc*(-4D0*kak3*kakb*m3**2 - 2D0*kak3*m1**2*m3**2 + 2D0*kbk3*m1**2*m3**2 - 2D0*k1k3*kakb*m4**2 - 2D0*k1k3*m1**2*m4**2 + 2D0*k1k3*m1**2*mq1**2 + 
     &                   4D0*k1k3*kakb*pap1 - 4D0*kakb*m3**2*pap1 - 4D0*m1**2*m3**2*pap1 + 2D0*kbk3*m4**2*pap1 - 4D0*kak3*mq1**2*pap1 - 4D0*kbk3*pap1**2 - 4D0*mq1**2*pap1**2 + 
     &                   4D0*k1k3*kakb*pap2 + 2D0*k1k3*m1**2*pap2 - 2D0*kbk3*m3**2*pap2 - 2D0*k1k3*mq1**2*pap2 - 4D0*kbk3*pap1*pap2 - 4D0*k1k3*m1**2*pbp1 - 2D0*kak3*m4**2*pbp1 + 
     &                   4D0*kak3*pap1*pbp1 + 4D0*k1k3*pap2*pbp1 - 4D0*pap1*pap2*pbp1 - 2D0*k1k3*m1**2*pbp2 + 2D0*kak3*m3**2*pbp2 - 4D0*k1k3*pap1*pbp2 + 4D0*kak3*pap1*pbp2 + 
     &                   4D0*pap1**2*pbp2 + 2D0*k2k3*(kakb*m3**2 + mq1**2*pap1 + m1**2*(m3**2 - 1D0*pap1 - 1D0*pbp1) + 2D0*pap2*pbp1 - 2D0*pap1*pbp2) + 
     &                   (kak3*m1**2*m3**2*m4**2 + kbk3*m1**2*m3**2*m4**2 + k1k3*kakb*m4**4 + k1k3*m1**2*m4**4 - 1D0*kakb*m3**2*m4**4 - 1D0*m1**2*m3**2*m4**4 + 
     &                      m1**2*m3**2*m4**2*mq1**2 + m1**2*m3**2*m4**2*pap1 - 1D0*kbk3*m4**4*pap1 - 2D0*k1k3*kakb*m4**2*pap2 - 1D0*k1k3*m1**2*m4**2*pap2 + 
     &                      2D0*kakb*m3**2*m4**2*pap2 + kbk3*m3**2*m4**2*pap2 + m1**2*m3**2*m4**2*pap2 + m3**2*m4**2*mq1**2*pap2 + 2D0*kbk3*m4**2*pap1*pap2 - 
     &                      2D0*kbk3*m3**2*pap2**2 - 2D0*m3**2*mq1**2*pap2**2 + m1**2*m3**2*m4**2*pbp1 + kak3*m4**4*pbp1 - 2D0*kak3*m4**2*pap2*pbp1 + m3**2*m4**2*pap2*pbp1 - 
     &                      2D0*m3**2*pap2**2*pbp1 + (k1k3*m1**2*m4**2 + m3**2*(-1D0*m4**2*(kak3 + m1**2 + pap1) + 2D0*(kak3 + pap1)*pap2))*pbp2)*xv + 
     &                   k2k3*(-1D0*kakb*m3**2*(m4**2 - 2D0*pap2) + m1**2*(m4**2*(pap1 + pbp1) - 1D0*m3**2*(m4**2 + 2D0*pbp2)))*xv + 
     &                   k1k2*(2D0*kbk3*m1**2 + 4D0*k2k3*(kakb + m1**2) + 4D0*kakb*m3**2 + 4D0*m1**2*m3**2 - 4D0*kakb*pap1 + 4D0*kbk3*pap1 - 4D0*m1**2*pap1 + 4D0*mq1**2*pap1 + 
     &                      4D0*pap2*pbp1 - 2D0*kak3*(2D0*kakb + m1**2 - 1D0*mq1**2 + 2D0*pbp1) - 4D0*pap1*pbp2 + 2D0*k2k3*m1**2*(-1D0*mq1**2 + pap1 + pbp1)*xv + 
     &                      (-1D0*kbk3*m1**2*m4**2 - 1D0*kakb*m3**2*m4**2 - 1D0*m1**2*m3**2*m4**2 - 2D0*kbk3*m4**2*pap1 - 2D0*m4**2*mq1**2*pap1 + 2D0*kakb*m3**2*pap2 + 
     &                         4D0*kbk3*pap1*pap2 + 4D0*mq1**2*pap1*pap2 - 1D0*kak3*(m1**2*m4**2 + 2D0*(m4**2 - 2D0*pap2)*(mq1**2 - 1D0*pbp1)) - 2D0*m1**2*m3**2*pbp2 + 
     &                         2D0*k1k3*(kakb*(m4**2 - 2D0*pap2) + m1**2*(m4**2 - 1D0*pap2 + pbp2)))*xv) + 2D0*k1k2**2*(2D0*kakb + m1**2*(2D0 - 1D0*(kak3 + kbk3 + mq1**2)*xv))))))
     &        + A1R*(A3Lc*(B1L*m3*(D3Lc*m1*(-2D0*(k1k2*(kak3 + kakb + kbk3 + m1**2) + (kakb + m1**2)*(k1k3 + k2k3 + m3**2) + (kak3 + kbk3)*m4**2) + 
     &                   2D0*(k1k2 + 2D0*k2k3 - 2D0*kak3)*mq1**2 + 2D0*(-1D0*k2k3 + 2D0*kak3 + kbk3 + m4**2 - 1D0*mq1**2)*pap1 + 2D0*(k1k3 + 2D0*kak3 + kbk3 + m3**2)*pap2 + 
     &                   2D0*(-1D0*k2k3 + kak3 + m4**2 - 1D0*pap2)*pbp1 + 2D0*(k1k3 + kak3 + m3**2 + pap1)*pbp2 + 
     &                   (kakb*m3**2*m4**2 + m1**2*m3**2*m4**2 + 2D0*kak3*m4**4 + 2D0*kbk3*m4**4 + kak3*m4**2*mq1**2 + k2k3*m4**2*(2D0*(kakb + m1**2) - 1D0*mq1**2) + m4**4*pap1 + 
     &                      m4**2*mq1**2*pap1 - 4D0*kak3*m4**2*pap2 - 2D0*kbk3*m4**2*pap2 + m3**2*m4**2*pap2 - 2D0*m4**2*pap1*pap2 - 2D0*m3**2*pap2**2 + m4**4*pbp1 - 
     &                      1D0*m4**2*pap2*pbp1 + k1k2*(2D0*k2k3*(kakb + m1**2) + m4**2*(kakb + 2D0*kbk3 + m1**2 + mq1**2) - 2D0*(kbk3 + mq1**2)*pap2 + 
     &                         2D0*kak3*(m4**2 - 2D0*pap2 - 1D0*pbp2)) - 1D0*(m4**2*(2D0*kak3 - 1D0*m3**2 + pap1) + 2D0*m3**2*pap2)*pbp2)*xv) - 
     &                1D0*D3Rc*mq3c*(-2D0*((kakb + m1**2)*(2D0*k1k2 + k1k3 + k2k3 + 2D0*m3**2) + mq1**2*(kak3 + 2D0*pap1) + kbk3*(pap1 + pap2) - 1D0*(kak3 - 2D0*pap2)*pbp1) + 
     &                   2D0*(kak3 + 2D0*pap1)*pbp2 + (2D0*k2k3*(kakb + m1**2)*m4**2 - 2D0*k2k3*mq1**2*pap2 + 
     &                      k1k2*((kakb + m1**2)*(2D0*k2k3 + m4**2) + 2D0*kbk3*pap2 - 2D0*kak3*pbp2) + 
     &                      m4**2*((kakb + m1**2)*m3**2 + mq1**2*(kak3 + pap1) + pap2*(2D0*kbk3 + pbp1) - 1D0*(2D0*kak3 + pap1)*pbp2))*xv)) + 
     &             B1R*mq1*(D3Rc*mq3c*(2D0*k1k3*(kakb + m1**2 - 2D0*pap1) - 2D0*((2D0*(k1k2 + k2k3) - 1D0*kbk3)*pap1 + m3**2*(kakb + m1**2 + pap1 - 1D0*pap2) + kak3*pbp1) + 
     &                   k1k3*(2D0*k1k2 + m4**2)*pap2*xv + (-2D0*k1k2**2*(kak3 + kakb + m1**2) + m4**2*(k2k3*pap1 + m3**2*(kak3 + kakb + m1**2 + pap1)) - 
     &                      1D0*m3**2*(2D0*k2k3 + m4**2)*pap2 - 1D0*k1k2*(kak3*m4**2 + 2D0*k2k3*(kakb + m1**2 - 1D0*pap1) + 2D0*pap2*(kbk3 + m3**2 + pbp1)) + 
     &                      2D0*k1k2*(kak3 + pap1)*pbp2)*xv) + D3Lc*m1*(-2D0*k1k2*(2D0*k2k3 + kbk3 + 2D0*(m3**2 + pbp1)) + 
     &                   2D0*(-1D0*k2k3*(m3**2 + pbp1) + kak3*(m3**2 + 2D0*pap1 + pbp1) + pap1*(kbk3 + 2D0*(m3**2 + pap1 + pbp1)) + k1k3*(-1D0*kakb - 1D0*m1**2 + m4**2 + pbp2)) - 
     &                   1D0*(kak3*m3**2*m4**2 + kakb*m3**2*m4**2 + m1**2*m3**2*m4**2 + k1k3*m4**4 - 1D0*m3**2*m4**4 + m3**2*m4**2*pap1 + k2k3*m4**2*(-1D0*m3**2 + pap1) - 
     &                      1D0*k1k3*m4**2*pap2 + 2D0*m3**2*m4**2*pap2 - 2D0*m3**2*pap2**2 + m3**2*(m4**2 - 2D0*pap2)*pbp2)*xv + 
     &                   k1k2*((2D0*kbk3 + m3**2)*m4**2 - 2D0*kbk3*pap2 + 2D0*k1k3*(-1D0*m4**2 + pap2) + 
     &                      2D0*(k2k3*(kakb + m1**2 - 1D0*pap1) - 2D0*pap1*pap2 - 1D0*pap2*pbp1 + m4**2*(pap1 + pbp1) - 1D0*pap1*pbp2) + kak3*(3D0*m4**2 - 2D0*(2D0*pap2 + pbp2)))*
     &                    xv + 2D0*k1k2**2*(-2D0 + (kak3 + kakb + m1**2)*xv)))) + 
     &          A3Rc*(B1R*m3*mq1*(-2D0*D3Rc*(k1k2*(kak3 - 1D0*kakb) + 2D0*kak3*kakb + k1k3*m1**2 + 2D0*kak3*m1**2 + m1**2*m3**2 + kak3*m4**2 - 2D0*kak3*pap1 + m1**2*pap1 - 
     &                   1D0*m4**2*pap1 + k2k3*(-2D0*kakb - 1D0*m1**2 + pap1) - 1D0*k1k3*pap2 - 2D0*kak3*pap2 - 1D0*m3**2*pap2 + m1**2*pbp1 - 1D0*pap2*pbp1 + pap1*pbp2) + 
     &                D3Rc*(kak3*m1**2*m4**2 + kbk3*m1**2*m4**2 + m1**2*m3**2*m4**2 + 2D0*kak3*m4**4 + m1**2*m4**2*pap1 + m4**4*pap1 + 
     &                   k1k2*((2D0*kak3 + kakb)*(m4**2 - 2D0*pap2) + 2D0*m1**2*(k2k3 + m4**2 - 1D0*pap2)) - 4D0*kak3*m4**2*pap2 + kbk3*m4**2*pap2 + m3**2*m4**2*pap2 - 
     &                   2D0*m4**2*pap1*pap2 - 2D0*kbk3*pap2**2 - 2D0*m3**2*pap2**2 + m1**2*m4**2*pbp1 + m4**2*pap2*pbp1 - 2D0*pap2**2*pbp1 + 
     &                   k2k3*(-1D0*kakb*(m4**2 - 2D0*pap2) + m1**2*(m4**2 - 2D0*pbp2)) - 1D0*(kak3 + pap1)*(m4**2 - 2D0*pap2)*pbp2)*xv + 
     &                D3Lc*m1*mq3c*(2D0*(2D0*k1k2 + k1k3 + k2k3 + kak3 + kbk3 + 2D0*(m3**2 + pap1 + pbp1)) - 
     &                   1D0*(k1k2*(2D0*k2k3 + m4**2) + m4**2*(kak3 + kbk3 + m3**2 + pap1 + pbp1) + 2D0*k2k3*(m4**2 - 1D0*pap2 - 1D0*pbp2))*xv)) - 
     &             1D0*B1L*(D3Lc*m1*mq3c*(4D0*(k1k2 + k2k3)*(pap1 + pbp1) - 2D0*k1k3*(mq1**2 - 2D0*(pap1 + pbp1)) + 
     &                   (2D0*k1k2**2*(kak3 + kbk3 + mq1**2) + k1k2*((kak3 + kbk3)*m4**2 + 2D0*k2k3*(mq1**2 - 1D0*pap1 - 1D0*pbp1)) - 1D0*k2k3*m4**2*(pap1 + pbp1))*xv - 
     &                   1D0*k1k3*(2D0*k1k2 + m4**2)*(pap2 + pbp2)*xv + m3**2*
     &                    (2D0*(mq1**2 + pap1 - 1D0*pap2 + pbp1 - 1D0*pbp2) - 
     &                      1D0*(m4**2*(kak3 + kbk3 + mq1**2 + pap1 - 1D0*pap2 + pbp1 - 1D0*pbp2) - 2D0*(k1k2 + k2k3)*(pap2 + pbp2))*xv)) + 
     &                D3Rc*(-4D0*kak3*kakb*m3**2 - 2D0*kak3*m1**2*m3**2 + 2D0*kbk3*m1**2*m3**2 - 2D0*k1k3*kakb*m4**2 - 2D0*k1k3*m1**2*m4**2 + 2D0*k1k3*m1**2*mq1**2 + 
     &                   4D0*k1k3*kakb*pap1 - 4D0*kakb*m3**2*pap1 - 4D0*m1**2*m3**2*pap1 + 2D0*kbk3*m4**2*pap1 - 4D0*kak3*mq1**2*pap1 - 4D0*kbk3*pap1**2 - 4D0*mq1**2*pap1**2 + 
     &                   4D0*k1k3*kakb*pap2 + 2D0*k1k3*m1**2*pap2 - 2D0*kbk3*m3**2*pap2 - 2D0*k1k3*mq1**2*pap2 - 4D0*kbk3*pap1*pap2 - 4D0*k1k3*m1**2*pbp1 - 2D0*kak3*m4**2*pbp1 + 
     &                   4D0*kak3*pap1*pbp1 + 4D0*k1k3*pap2*pbp1 - 4D0*pap1*pap2*pbp1 - 2D0*k1k3*m1**2*pbp2 + 2D0*kak3*m3**2*pbp2 - 4D0*k1k3*pap1*pbp2 + 4D0*kak3*pap1*pbp2 + 
     &                   4D0*pap1**2*pbp2 + 2D0*k2k3*(kakb*m3**2 + mq1**2*pap1 + m1**2*(m3**2 - 1D0*pap1 - 1D0*pbp1) + 2D0*pap2*pbp1 - 2D0*pap1*pbp2) + 
     &                   (kak3*m1**2*m3**2*m4**2 + kbk3*m1**2*m3**2*m4**2 + k1k3*kakb*m4**4 + k1k3*m1**2*m4**4 - 1D0*kakb*m3**2*m4**4 - 1D0*m1**2*m3**2*m4**4 + 
     &                      m1**2*m3**2*m4**2*mq1**2 + m1**2*m3**2*m4**2*pap1 - 1D0*kbk3*m4**4*pap1 - 2D0*k1k3*kakb*m4**2*pap2 - 1D0*k1k3*m1**2*m4**2*pap2 + 
     &                      2D0*kakb*m3**2*m4**2*pap2 + kbk3*m3**2*m4**2*pap2 + m1**2*m3**2*m4**2*pap2 + m3**2*m4**2*mq1**2*pap2 + 2D0*kbk3*m4**2*pap1*pap2 - 
     &                      2D0*kbk3*m3**2*pap2**2 - 2D0*m3**2*mq1**2*pap2**2 + m1**2*m3**2*m4**2*pbp1 + kak3*m4**4*pbp1 - 2D0*kak3*m4**2*pap2*pbp1 + m3**2*m4**2*pap2*pbp1 - 
     &                      2D0*m3**2*pap2**2*pbp1 + (k1k3*m1**2*m4**2 + m3**2*(-1D0*m4**2*(kak3 + m1**2 + pap1) + 2D0*(kak3 + pap1)*pap2))*pbp2)*xv + 
     &                   k2k3*(-1D0*kakb*m3**2*(m4**2 - 2D0*pap2) + m1**2*(m4**2*(pap1 + pbp1) - 1D0*m3**2*(m4**2 + 2D0*pbp2)))*xv + 
     &                   k1k2*(2D0*kbk3*m1**2 + 4D0*k2k3*(kakb + m1**2) + 4D0*kakb*m3**2 + 4D0*m1**2*m3**2 - 4D0*kakb*pap1 + 4D0*kbk3*pap1 - 4D0*m1**2*pap1 + 4D0*mq1**2*pap1 + 
     &                      4D0*pap2*pbp1 - 2D0*kak3*(2D0*kakb + m1**2 - 1D0*mq1**2 + 2D0*pbp1) - 4D0*pap1*pbp2 + 2D0*k2k3*m1**2*(-1D0*mq1**2 + pap1 + pbp1)*xv + 
     &                      (-1D0*kbk3*m1**2*m4**2 - 1D0*kakb*m3**2*m4**2 - 1D0*m1**2*m3**2*m4**2 - 2D0*kbk3*m4**2*pap1 - 2D0*m4**2*mq1**2*pap1 + 2D0*kakb*m3**2*pap2 + 
     &                         4D0*kbk3*pap1*pap2 + 4D0*mq1**2*pap1*pap2 - 1D0*kak3*(m1**2*m4**2 + 2D0*(m4**2 - 2D0*pap2)*(mq1**2 - 1D0*pbp1)) - 2D0*m1**2*m3**2*pbp2 + 
     &                         2D0*k1k3*(kakb*(m4**2 - 2D0*pap2) + m1**2*(m4**2 - 1D0*pap2 + pbp2)))*xv) + 2D0*k1k2**2*(2D0*kakb + m1**2*(2D0 - 1D0*(kak3 + kbk3 + mq1**2)*xv)))))))
     &     )/(k1k3*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2 - cI*mq1*wquark)*(2D0*k1k2 + m3**2 + m4**2 - 1D0*mq1**2 - cI*mq1*wquark)*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)))


        radvSfUi = dreal((-1D0*GS**2*(A1R*(A3Rc*(B1R*m3*mq1*(2D0*D3Rc*(-1D0*k2k3*kakb + (m1**2 - 1D0*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + kak3*pbp2) + 
     &                D3Rc*(k2k3*kakb*(m4**2 - 2D0*pap2) - 1D0*(m1**2*m4**2 + (m4**2 - 2D0*pap2)*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*kak3*(m4**2 - 2D0*pap2)*pbp2)*xv + 
     &                D3Lc*m1*mq3c*(k1k3 - 1D0*kbk3 - 2D0*pbp1)*(-4D0 + m4**2*xv)) + 
     &             B1L*(-1D0*D3Lc*m1*mq3c*(-2D0*k1k3*(kakb + m2**2 + pap1 + pbp1) + 2D0*kbk3*(pap1 + 2D0*pbp1) + 2D0*pbp1*(kak3 + 2D0*(pap1 + pbp1)) + 
     &                   k1k3*(m4**2*(kakb + m2**2 + pap1 + pbp1) - 2D0*(k1k2 + pbp2)*(pap2 + pbp2))*xv + 
     &                   (kbk3*(-1D0*m4**2*(pap1 + 2D0*pbp1) + 2D0*k1k2*(pap2 + pbp2)) + pbp1*(-1D0*m4**2*(kak3 + 2D0*(pap1 + pbp1)) + 2D0*(2D0*k1k2 + k2k3)*(pap2 + pbp2)))*xv) + 
     &                D3Rc*(4D0*(kakb + m1**2)*(kbk3*pap1 + k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + (-1D0*k2k3 + kak3 + 2D0*pap1)*pbp1 - 1D0*k1k3*(kakb + pap1 - 1D0*pbp2)) + 
     &                   (-1D0*kbk3*m1**2*m4**2*pap1 + k2k3*m2**2*m4**2*pap1 - 2D0*k2k3*m2**2*pap1*pap2 + k2k3*m1**2*m4**2*pbp1 - 1D0*kak3*m1**2*m4**2*pbp1 - 
     &                      2D0*kbk3*m1**2*m4**2*pbp1 - 2D0*m1**2*m4**2*pap1*pbp1 - 2D0*kbk3*m4**2*pap2*pbp1 + 4D0*kbk3*pap2**2*pbp1 - 2D0*m1**2*m4**2*pbp1**2 - 
     &                      2D0*m4**2*pap2*pbp1**2 + 4D0*pap2**2*pbp1**2 + 2D0*(k2k3*m1**2 + (kak3 + pap1)*(m4**2 - 2D0*pap2))*pbp1*pbp2 + 
     &                      k1k3*((m4**2 - 2D0*pap2)*(pap2*(m2**2 + pbp1) - 1D0*pap1*pbp2) + kakb*(m1**2*m4**2 - 2D0*m4**2*pbp2 + 4D0*pap2*pbp2) + 
     &                         m1**2*(m4**2*(m2**2 + pap1 + pbp1) - 1D0*m4**2*pbp2 - 2D0*pbp2**2)) + 
     &                      k1k2*(kbk3*m1**2*m4**2 - 1D0*kak3*m2**2*m4**2 + 2D0*kak3*m2**2*pap2 + 2D0*m1**2*m4**2*pbp1 + 2D0*kakb*(m4**2 - 2D0*pap2)*(kbk3 + pbp1) + 
     &                         2D0*m1**2*(kbk3 + 2D0*pbp1)*pbp2 - 1D0*k1k3*(kakb*(m4**2 - 2D0*pap2) + m1**2*(m4**2 + 2D0*pbp2))))*xv))) + 
     &          A3Lc*(B1R*mq1*(D3Lc*m1*((k2k3 - 1D0*kak3)*pbp1 - 1D0*pap1*(kbk3 + 2D0*pbp1) + k1k2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) + k1k3*(kakb + pap1 - 1D0*pbp2))*
     &                 (-4D0 + m4**2*xv) + D3Rc*mq3c*(2D0*k1k3*(kakb + pap1) - 2D0*(kak3 + 2D0*pap1)*pbp1 + (m4**2*(kak3 + 2D0*pap1) - 2D0*(2D0*k1k2 + k2k3)*pap2)*pbp1*xv - 
     &                   1D0*k1k3*(m4**2*(kakb + pap1) - 2D0*pap2*(k1k2 + pbp2))*xv + kbk3*(-2D0*pap1 + m4**2*pap1*xv - 2D0*k1k2*pap2*xv))) + 
     &             B1L*m3*(D3Rc*mq3c*(4D0*(kakb + m1**2)*(kbk3 + 2D0*pbp1) - 
     &                   1D0*(kbk3*m1**2*m4**2 - 1D0*kak3*m2**2*m4**2 + 2D0*k2k3*m2**2*pap2 + 2D0*m1**2*m4**2*pbp1 + 2D0*kakb*m4**2*(kbk3 + pbp1) - 2D0*k2k3*kakb*pbp2 - 
     &                      2D0*kbk3*pap2*pbp2 + 2D0*kak3*pbp2**2)*xv + k1k3*(kakb + m1**2)*(-4D0 + m4**2*xv)) + 
     &                D3Lc*m1*(-2D0*k2k3*(kakb + m2**2) - 2D0*kbk3*(2D0*kakb + m1**2 - 1D0*pap2) - 4D0*(kakb + m1**2 - 1D0*pap2)*pbp1 + 4D0*(kbk3 + pbp1)*pbp2 + 
     &                   2D0*kak3*(m2**2 + pbp2) + k2k3*(m2**2*m4**2 + kakb*(m4**2 - 2D0*(pap2 + pbp2)))*xv + 
     &                   (-1D0*kak3*m2**2*m4**2 + 2D0*m1**2*m4**2*pbp1 + 2D0*m4**2*pap2*pbp1 - 4D0*pap2**2*pbp1 + 2D0*kakb*m4**2*(kbk3 + pbp1) - 
     &                      1D0*(m4**2 - 2D0*pap2)*(kak3 - 2D0*pbp1)*pbp2 + 2D0*kak3*pbp2**2 + kbk3*(m1**2*m4**2 + pap2*(m4**2 - 2D0*(pap2 + pbp2))))*xv - 
     &                   1D0*k1k3*(2D0*(pap2 + pbp2) + (m4**2 - 2D0*pap2)*(pap2 + pbp2)*xv + kakb*(-2D0 + m4**2*xv) + m1**2*(-2D0 + m4**2*xv)))))) + 
     &       A1L*(A3Lc*(B1L*m3*mq1*(2D0*D3Lc*(-1D0*k2k3*kakb + (m1**2 - 1D0*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + kak3*pbp2) + 
     &                D3Lc*(k2k3*kakb*(m4**2 - 2D0*pap2) - 1D0*(m1**2*m4**2 + (m4**2 - 2D0*pap2)*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*kak3*(m4**2 - 2D0*pap2)*pbp2)*xv + 
     &                D3Rc*m1*mq3c*(k1k3 - 1D0*kbk3 - 2D0*pbp1)*(-4D0 + m4**2*xv)) + 
     &             B1R*(-1D0*D3Rc*m1*mq3c*(-2D0*k1k3*(kakb + m2**2 + pap1 + pbp1) + 2D0*kbk3*(pap1 + 2D0*pbp1) + 2D0*pbp1*(kak3 + 2D0*(pap1 + pbp1)) + 
     &                   k1k3*(m4**2*(kakb + m2**2 + pap1 + pbp1) - 2D0*(k1k2 + pbp2)*(pap2 + pbp2))*xv + 
     &                   (kbk3*(-1D0*m4**2*(pap1 + 2D0*pbp1) + 2D0*k1k2*(pap2 + pbp2)) + pbp1*(-1D0*m4**2*(kak3 + 2D0*(pap1 + pbp1)) + 2D0*(2D0*k1k2 + k2k3)*(pap2 + pbp2)))*xv) + 
     &                D3Lc*(4D0*(kakb + m1**2)*(kbk3*pap1 + k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + (-1D0*k2k3 + kak3 + 2D0*pap1)*pbp1 - 1D0*k1k3*(kakb + pap1 - 1D0*pbp2)) + 
     &                   (-1D0*kbk3*m1**2*m4**2*pap1 + k2k3*m2**2*m4**2*pap1 - 2D0*k2k3*m2**2*pap1*pap2 + k2k3*m1**2*m4**2*pbp1 - 1D0*kak3*m1**2*m4**2*pbp1 - 
     &                      2D0*kbk3*m1**2*m4**2*pbp1 - 2D0*m1**2*m4**2*pap1*pbp1 - 2D0*kbk3*m4**2*pap2*pbp1 + 4D0*kbk3*pap2**2*pbp1 - 2D0*m1**2*m4**2*pbp1**2 - 
     &                      2D0*m4**2*pap2*pbp1**2 + 4D0*pap2**2*pbp1**2 + 2D0*(k2k3*m1**2 + (kak3 + pap1)*(m4**2 - 2D0*pap2))*pbp1*pbp2 + 
     &                      k1k3*((m4**2 - 2D0*pap2)*(pap2*(m2**2 + pbp1) - 1D0*pap1*pbp2) + kakb*(m1**2*m4**2 - 2D0*m4**2*pbp2 + 4D0*pap2*pbp2) + 
     &                         m1**2*(m4**2*(m2**2 + pap1 + pbp1) - 1D0*m4**2*pbp2 - 2D0*pbp2**2)) + 
     &                      k1k2*(kbk3*m1**2*m4**2 - 1D0*kak3*m2**2*m4**2 + 2D0*kak3*m2**2*pap2 + 2D0*m1**2*m4**2*pbp1 + 2D0*kakb*(m4**2 - 2D0*pap2)*(kbk3 + pbp1) + 
     &                         2D0*m1**2*(kbk3 + 2D0*pbp1)*pbp2 - 1D0*k1k3*(kakb*(m4**2 - 2D0*pap2) + m1**2*(m4**2 + 2D0*pbp2))))*xv))) + 
     &          A3Rc*(B1L*mq1*(D3Rc*m1*((k2k3 - 1D0*kak3)*pbp1 - 1D0*pap1*(kbk3 + 2D0*pbp1) + k1k2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) + k1k3*(kakb + pap1 - 1D0*pbp2))*
     &                 (-4D0 + m4**2*xv) + D3Lc*mq3c*(2D0*k1k3*(kakb + pap1) - 2D0*(kak3 + 2D0*pap1)*pbp1 + (m4**2*(kak3 + 2D0*pap1) - 2D0*(2D0*k1k2 + k2k3)*pap2)*pbp1*xv - 
     &                   1D0*k1k3*(m4**2*(kakb + pap1) - 2D0*pap2*(k1k2 + pbp2))*xv + kbk3*(-2D0*pap1 + m4**2*pap1*xv - 2D0*k1k2*pap2*xv))) + 
     &             B1R*m3*(D3Lc*mq3c*(4D0*(kakb + m1**2)*(kbk3 + 2D0*pbp1) - 
     &                   1D0*(kbk3*m1**2*m4**2 - 1D0*kak3*m2**2*m4**2 + 2D0*k2k3*m2**2*pap2 + 2D0*m1**2*m4**2*pbp1 + 2D0*kakb*m4**2*(kbk3 + pbp1) - 2D0*k2k3*kakb*pbp2 - 
     &                      2D0*kbk3*pap2*pbp2 + 2D0*kak3*pbp2**2)*xv + k1k3*(kakb + m1**2)*(-4D0 + m4**2*xv)) + 
     &                D3Rc*m1*(-2D0*k2k3*(kakb + m2**2) - 2D0*kbk3*(2D0*kakb + m1**2 - 1D0*pap2) - 4D0*(kakb + m1**2 - 1D0*pap2)*pbp1 + 4D0*(kbk3 + pbp1)*pbp2 + 
     &                   2D0*kak3*(m2**2 + pbp2) + k2k3*(m2**2*m4**2 + kakb*(m4**2 - 2D0*(pap2 + pbp2)))*xv + 
     &                   (-1D0*kak3*m2**2*m4**2 + 2D0*m1**2*m4**2*pbp1 + 2D0*m4**2*pap2*pbp1 - 4D0*pap2**2*pbp1 + 2D0*kakb*m4**2*(kbk3 + pbp1) - 
     &                      1D0*(m4**2 - 2D0*pap2)*(kak3 - 2D0*pbp1)*pbp2 + 2D0*kak3*pbp2**2 + kbk3*(m1**2*m4**2 + pap2*(m4**2 - 2D0*(pap2 + pbp2))))*xv - 
     &                   1D0*k1k3*(2D0*(pap2 + pbp2) + (m4**2 - 2D0*pap2)*(pap2 + pbp2)*xv + kakb*(-2D0 + m4**2*xv) + m1**2*(-2D0 + m4**2*xv))))))))/
     &   (k1k3*kbk3*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2 - cI*mq1*wquark)*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)))


        radvSfUf = dreal((2D0*GS**2*(A1R*(A3Lc*(B1L*m3*(k1k3 + m3**2)*(-1D0*D3Rc*(kakb + m1**2)*mq3c*(-4D0 + m4**2*xv) + 
     &                D3Lc*m1*(2D0*(pap2 + pbp2) + (m4**2 - 2D0*pap2)*(pap2 + pbp2)*xv + kakb*(-2D0 + m4**2*xv) + m1**2*(-2D0 + m4**2*xv))) + 
     &             B1R*mq1*(-1D0*D3Lc*m1*(k1k3*(k2k3 - 1D0*kak3) + m3**2*(-1D0*k1k2 - 1D0*k2k3 + kak3 + pap1))*(-4D0 + m4**2*xv) + 
     &                D3Rc*mq3c*(k1k3*(2D0*kak3 - 1D0*kak3*m4**2*xv + 2D0*k2k3*pap2*xv) + m3**2*(-2D0*(kak3 + pap1) + (m4**2*(kak3 + pap1) - 2D0*(k1k2 + k2k3)*pap2)*xv)))) + 
     &          A3Rc*(B1R*m3*(k1k3 + m3**2)*mq1*(D3Lc*m1*mq3c*(4D0 - 1D0*m4**2*xv) + D3Rc*(m1**2*(-2D0 + m4**2*xv) + pap2*(2D0 + m4**2*xv - 2D0*pap2*xv))) + 
     &             B1L*(4D0*D3Rc*(kakb + m1**2)*(k1k3*(k2k3 - 1D0*kak3) + m3**2*(-1D0*k1k2 - 1D0*k2k3 + kak3 + pap1)) + 
     &                D3Rc*(k1k3*(kbk3*(m1**2*m4**2 + (m4**2 - 2D0*pap2)*pap2) + kak3*(m1**2*m4**2 - 1D0*m4**2*pbp2 + 2D0*pap2*pbp2) - 
     &                      1D0*k2k3*(kakb*(m4**2 - 2D0*pap2) + m1**2*(m4**2 + 2D0*pbp2))) + 
     &                   m3**2*(-1D0*kak3*m1**2*m4**2 - 1D0*kbk3*m1**2*m4**2 - 1D0*m1**2*m4**2*pap1 - 1D0*kbk3*m4**2*pap2 + 2D0*kbk3*pap2**2 - 1D0*m1**2*m4**2*pbp1 - 
     &                      1D0*m4**2*pap2*pbp1 + 2D0*pap2**2*pbp1 + (kak3 + pap1)*(m4**2 - 2D0*pap2)*pbp2 + k1k2*(kakb*(m4**2 - 2D0*pap2) + m1**2*(m4**2 + 2D0*pbp2)) + 
     &                      k2k3*(kakb*(m4**2 - 2D0*pap2) + m1**2*(m4**2 + 2D0*pbp2))))*xv + 
     &                D3Lc*m1*mq3c*(2D0*k1k3*(kak3 + kbk3) - 1D0*k1k3*((kak3 + kbk3)*m4**2 - 2D0*k2k3*(pap2 + pbp2))*xv + 
     &                   m3**2*(-2D0*(kak3 + kbk3 + pap1 + pbp1) + (m4**2*(kak3 + kbk3 + pap1 + pbp1) - 2D0*(k1k2 + k2k3)*(pap2 + pbp2))*xv))))) + 
     &       A1L*(A3Rc*(B1R*m3*(k1k3 + m3**2)*(-1D0*D3Lc*(kakb + m1**2)*mq3c*(-4D0 + m4**2*xv) + 
     &                D3Rc*m1*(2D0*(pap2 + pbp2) + (m4**2 - 2D0*pap2)*(pap2 + pbp2)*xv + kakb*(-2D0 + m4**2*xv) + m1**2*(-2D0 + m4**2*xv))) + 
     &             B1L*mq1*(-1D0*D3Rc*m1*(k1k3*(k2k3 - 1D0*kak3) + m3**2*(-1D0*k1k2 - 1D0*k2k3 + kak3 + pap1))*(-4D0 + m4**2*xv) + 
     &                D3Lc*mq3c*(k1k3*(2D0*kak3 - 1D0*kak3*m4**2*xv + 2D0*k2k3*pap2*xv) + m3**2*(-2D0*(kak3 + pap1) + (m4**2*(kak3 + pap1) - 2D0*(k1k2 + k2k3)*pap2)*xv)))) + 
     &          A3Lc*(B1L*m3*(k1k3 + m3**2)*mq1*(D3Rc*m1*mq3c*(4D0 - 1D0*m4**2*xv) + D3Lc*(m1**2*(-2D0 + m4**2*xv) + pap2*(2D0 + m4**2*xv - 2D0*pap2*xv))) + 
     &             B1R*(4D0*D3Lc*(kakb + m1**2)*(k1k3*(k2k3 - 1D0*kak3) + m3**2*(-1D0*k1k2 - 1D0*k2k3 + kak3 + pap1)) + 
     &                D3Lc*(k1k3*(kbk3*(m1**2*m4**2 + (m4**2 - 2D0*pap2)*pap2) + kak3*(m1**2*m4**2 - 1D0*m4**2*pbp2 + 2D0*pap2*pbp2) - 
     &                      1D0*k2k3*(kakb*(m4**2 - 2D0*pap2) + m1**2*(m4**2 + 2D0*pbp2))) + 
     &                   m3**2*(-1D0*kak3*m1**2*m4**2 - 1D0*kbk3*m1**2*m4**2 - 1D0*m1**2*m4**2*pap1 - 1D0*kbk3*m4**2*pap2 + 2D0*kbk3*pap2**2 - 1D0*m1**2*m4**2*pbp1 - 
     &                      1D0*m4**2*pap2*pbp1 + 2D0*pap2**2*pbp1 + (kak3 + pap1)*(m4**2 - 2D0*pap2)*pbp2 + k1k2*(kakb*(m4**2 - 2D0*pap2) + m1**2*(m4**2 + 2D0*pbp2)) + 
     &                      k2k3*(kakb*(m4**2 - 2D0*pap2) + m1**2*(m4**2 + 2D0*pbp2))))*xv + 
     &                D3Rc*m1*mq3c*(2D0*k1k3*(kak3 + kbk3) - 1D0*k1k3*((kak3 + kbk3)*m4**2 - 2D0*k2k3*(pap2 + pbp2))*xv + 
     &                   m3**2*(-2D0*(kak3 + kbk3 + pap1 + pbp1) + (m4**2*(kak3 + kbk3 + pap1 + pbp1) - 2D0*(k1k2 + k2k3)*(pap2 + pbp2))*xv)))))))/
     &   (k1k3**2*(2D0*kakb + m1**2 + m2**2 - 1D0*mq1**2 - cI*mq1*wquark)*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)))
                  
         
        NeuQ2qx_M2RadiationVSU = NeuQ2qx_M2RadiationVSU + 2.*CF*NC*(radvSiUi + radvSiUf + radvSpUi + radvSpUf + radvSfUi + radvSfUf)

        enddo
    
    
        end

! =================================================================
        double precision function NeuQ2qx_M2RadiationVTT()
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double precision radvTiTi, radvTiTp, radvTiTf, radvTiTv, radvTpTp, radvTpTf, radvTpTv, radvTfTf, radvTfTv, radvTvTv
        double complex A2L, A2Lc, A2R, A2Rc, C2L, C2Lc, CLL2, CLL2c
        double precision mq2, mq2c
        double precision xv, colorfactor
        integer sumT, sumTc

        if (xtype.eq.1) then
        
           xv = 0d0
           colorfactor = 4d0/3d0
           
        elseif (xtype.eq.2) then
        
           xv = 0d0
           colorfactor = 1d0
           
        elseif (xtype.eq.3) then
        
           xv = 1/(m4**2)
           colorfactor = 1d0
           
        elseif (xtype.eq.4) then

           xv = 1/(m4**2)
           colorfactor = 1d0
                   
        endif

        NeuQ2qx_M2RadiationVTT = 0d0

        do sumT = 0, 1
        do sumTc = 0, 1

        radvTiTi = 0d0        
        radvTiTp = 0d0
        radvTiTf = 0d0
        radvTiTv = 0d0
        radvTpTp = 0d0
        radvTpTf = 0d0
        radvTpTv = 0d0
        radvTfTf = 0d0
        radvTfTv = 0d0
        radvTvTv = 0d0

        A2L = AL(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2R = AR(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
        A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))

        C2L = CL(ptmeigen+sumT,sfmeigen,sftype,sfgen,xtype)
        C2Lc = dconjg(CL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))

        CLL2 =  CLL(ptmeigen+sumT,sfmeigen,sftype,sfgen,xtype)
        CLL2c =  dconjg(CLL(ptmeigen+sumTc,sfmeigen,sftype,sfgen,xtype))

        mq2 = MSf(ptmeigen+sumT,pttype,ptgen)
        mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

        

        radvTiTi = dreal((-8D0*C2L*C2Lc*GS**2*(kbk3 - 1D0*m2**2)*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*
     &     (4D0*k1k2 + 4D0*m1**2 + 4D0*m3**2 + m4**2 - 8D0*pap1 - 4D0*pap2 - 1D0*(2D0*k1k2 + m4**2 - 2D0*pap2)**2*xv))/
     &   ((m1**2 + m3**2 - 1D0*mq2**2 - 2D0*pap1 - (cI*mq2*wst(ptmeigen+sumT)))*(m1**2 + m3**2 - 1D0*mq2c**2 - 2D0*pap1 + (cI*mq2c*wst(ptmeigen+sumTc)))*(2D0*k1k2 + m1**2 - 1D0*m2**2 + m3**2 + m4**2 - 2D0*(pap1 + pap2))**2))


        radvTiTp = dreal((-2D0*C2L*C2Lc*GS**2*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*(k1k3 - 1D0*k2k3 - 1D0*kak3 + 2D0*kakb + kbk3 - 2D0*(m2**2 + pbp1 - 1D0*pbp2))*
     &     (2D0*k1k2 + 4D0*kakb + m4**2 - 2D0*(pap2 + 2D0*pbp1 + pbp2) - 1D0*(2D0*k1k2 + m4**2 - 2D0*pap2)*(m4**2 - 2D0*pbp2)*xv))/
     &   ((m1**2 + m3**2 - 1D0*mq2**2 - 2D0*pap1 - (cI*mq2*wst(ptmeigen+sumT)))*(m1**2 + m3**2 - 1D0*mq2c**2 - 2D0*pap1 + (cI*mq2c*wst(ptmeigen+sumTc)))*(2D0*k1k2 + m1**2 - 1D0*m2**2 + m3**2 + m4**2 - 2D0*(pap1 + pap2))*
     &     (-1D0*m2**2 - 1D0*m4**2 + mq2c**2 + 2D0*pbp2 - (cI*mq2c*wst(ptmeigen+sumTc)))))


        radvTiTf = dreal((2D0*C2L*C2Lc*GS**2*(A2L*(A2Rc*m1*m3*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + A2Lc*(k1k3*(kakb + pap1) - 1D0*kak3*pbp1 - 1D0*pap1*(kbk3 + 2D0*pbp1))) + 
     &       A2R*(A2Lc*m1*m3*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + A2Rc*(k1k3*(kakb + pap1) - 1D0*kak3*pbp1 - 1D0*pap1*(kbk3 + 2D0*pbp1))))*
     &     (2D0*k1k2 + 4D0*kakb + m4**2 - 2D0*(pap2 + 2D0*pbp1 + pbp2) - 1D0*(2D0*k1k2 + m4**2 - 2D0*pap2)*(m4**2 - 2D0*pbp2)*xv))/
     &   (k1k3*(m1**2 + m3**2 - 1D0*mq2**2 - 2D0*pap1 - (cI*mq2*wst(ptmeigen+sumT)))*(2D0*k1k2 + m1**2 - 1D0*m2**2 + m3**2 + m4**2 - 2D0*(pap1 + pap2))*(m2**2 + m4**2 - 1D0*mq2c**2 - 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc)))))


        radvTiTv = dreal((-2D0*C2L*CLL2c*GS*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*
     &     (2D0*k1k3 + k2k3 - 2D0*(kak3 - 2D0*kakb + 2D0*pbp1 + pbp2) - 1D0*(2D0*k1k2 + m4**2 - 2D0*pap2)*(k2k3 - 2D0*pbp2)*xv))/
     &   ((m1**2 + m3**2 - 1D0*mq2**2 - 2D0*pap1 - (cI*mq2*wst(ptmeigen+sumT)))*(2D0*k1k2 + m1**2 - 1D0*m2**2 + m3**2 + m4**2 - 2D0*(pap1 + pap2))*(2D0*k2k3 - 2D0*kbk3 + m2**2 + m4**2 - 1D0*mq2c**2 - 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc))))
     &   )


        radvTpTp = dreal((-2D0*C2L*C2Lc*GS**2*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*
     &     (-2D0*k1k2 - 2D0*kakb + m1**2 + m2**2 + m3**2 + m4**2 - 2D0*pap1 + 2D0*pap2 + 2D0*pbp1 - 2D0*pbp2)*(-4D0*m2**2 - 1D0*m4**2 + 4D0*pbp2 + (m4**2 - 2D0*pbp2)**2*xv))/
     &   ((m1**2 + m3**2 - 1D0*mq2**2 - 2D0*pap1 - (cI*mq2*wst(ptmeigen+sumT)))*(m1**2 + m3**2 - 1D0*mq2c**2 - 2D0*pap1 + (cI*mq2c*wst(ptmeigen+sumTc)))*(m2**2 + m4**2 - 1D0*mq2**2 - 2D0*pbp2 - (cI*mq2*wst(ptmeigen+sumT)))*(m2**2 + m4**2 - 1D0*mq2c**2 - 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc)))))


        radvTpTf = dreal((-1D0*C2L*C2Lc*GS**2*(A2L*(A2Rc*m1*m3*(-2D0*k1k2 + k1k3 - 1D0*k2k3 - 1D0*kak3 + kbk3 + 2D0*(m3**2 - 1D0*pap1 + pbp1)) + 
     &          A2Lc*(-1D0*k1k2*(kak3 + 2D0*pap1) + k1k3*(-1D0*kakb + m1**2 + pap2) + kak3*(m3**2 - 2D0*pap1 + pbp1) + pap1*(-1D0*k2k3 + kbk3 + 2D0*(m3**2 - 1D0*pap1 + pbp1)))) + 
     &       A2R*(A2Lc*m1*m3*(-2D0*k1k2 + k1k3 - 1D0*k2k3 - 1D0*kak3 + kbk3 + 2D0*(m3**2 - 1D0*pap1 + pbp1)) + 
     &          A2Rc*(-1D0*k1k2*(kak3 + 2D0*pap1) + k1k3*(-1D0*kakb + m1**2 + pap2) + kak3*(m3**2 - 2D0*pap1 + pbp1) + pap1*(-1D0*k2k3 + kbk3 + 2D0*(m3**2 - 1D0*pap1 + pbp1)))))*
     &     (4D0*m2**2 + m4**2 - 4D0*pbp2 - 1D0*(m4**2 - 2D0*pbp2)**2*xv))/
     &   (k1k3*(m1**2 + m3**2 - 1D0*mq2**2 - 2D0*pap1 - (cI*mq2*wst(ptmeigen+sumT)))*(m2**2 + m4**2 - 1D0*mq2c**2 - 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc)))*(-1D0*m2**2 - 1D0*m4**2 + mq2**2 + 2D0*pbp2 + (cI*mq2*wst(ptmeigen+sumT)))))


        radvTpTv = dreal((-2D0*C2L*CLL2c*GS*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*
     &     (-1D0*k1k2 - 2D0*kakb + 2D0*m2**2 + m4**2 + pap2 + 2D0*pbp1 - 3D0*pbp2 + (m4**2 - 2D0*pbp2)*(k1k2 - 1D0*m4**2 - 1D0*pap2 + pbp2)*xv))/
     &   ((m1**2 + m3**2 - 1D0*mq2**2 - 2D0*pap1 - (cI*mq2*wst(ptmeigen+sumT)))*(2D0*k2k3 - 2D0*kbk3 + m2**2 + m4**2 - 1D0*mq2c**2 - 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc)))*(-1D0*m2**2 - 1D0*m4**2 + mq2**2 + 2D0*pbp2 + (cI*mq2*wst(ptmeigen+sumT)))))


        radvTfTf = dreal((2D0*C2L*C2Lc*GS**2*(A2L*(A2Rc*m1*m3*(k1k3 + m3**2) + A2Lc*(-1D0*k1k3*kak3 + m3**2*(kak3 + pap1))) + 
     &       A2R*(A2Lc*m1*m3*(k1k3 + m3**2) + A2Rc*(-1D0*k1k3*kak3 + m3**2*(kak3 + pap1))))*(4D0*m2**2 + m4**2 - 4D0*pbp2 - 1D0*(m4**2 - 2D0*pbp2)**2*xv))/
     &   (k1k3**2*(m2**2 + m4**2 - 1D0*mq2**2 - 2D0*pbp2 - (cI*mq2*wst(ptmeigen+sumT)))*(m2**2 + m4**2 - 1D0*mq2c**2 - 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc)))))


        radvTfTv = dreal((C2L*CLL2c*GS*(A2L*(A2Rc*m1*m3*(-2D0*k1k2 - 1D0*k2k3 + 2D0*kbk3 + 4D0*pbp1 + (2D0*k1k2 + k2k3)*(m4**2 - 2D0*pbp2)*xv) + 
     &          A2Lc*(-1D0*k1k2*kak3 - 2D0*k1k3*kakb - 2D0*k1k2*pap1 - 1D0*k2k3*pap1 + 2D0*kbk3*pap1 + k1k3*pap2 + 2D0*kak3*pbp1 + 4D0*pap1*pbp1 + 
     &             (k1k2*kak3 + 2D0*k1k2*pap1 + k2k3*pap1 - 1D0*k1k3*pap2)*(m4**2 - 2D0*pbp2)*xv)) + 
     &       A2R*(A2Lc*m1*m3*(-2D0*k1k2 - 1D0*k2k3 + 2D0*kbk3 + 4D0*pbp1 + (2D0*k1k2 + k2k3)*(m4**2 - 2D0*pbp2)*xv) + 
     &          A2Rc*(-1D0*k1k2*kak3 - 2D0*k1k3*kakb - 2D0*k1k2*pap1 - 1D0*k2k3*pap1 + 2D0*kbk3*pap1 + k1k3*pap2 + 2D0*kak3*pbp1 + 4D0*pap1*pbp1 + 
     &             (k1k2*kak3 + 2D0*k1k2*pap1 + k2k3*pap1 - 1D0*k1k3*pap2)*(m4**2 - 2D0*pbp2)*xv))))/
     &   (k1k3*(m2**2 + m4**2 - 1D0*mq2**2 - 2D0*pbp2 - (cI*mq2*wst(ptmeigen+sumT)))*(2D0*k2k3 - 2D0*kbk3 + m2**2 + m4**2 - 1D0*mq2c**2 - 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc)))))


        radvTvTv = dreal((-2D0*CLL2*CLL2c*(A2Lc*A2R*m1*m3 + A2L*A2Rc*m1*m3 + A2L*A2Lc*pap1 + A2R*A2Rc*pap1)*(-4D0 + m4**2*xv))/
     &   ((2D0*k2k3 - 2D0*kbk3 + m2**2 + m4**2 - 1D0*mq2**2 - 2D0*pbp2 - (cI*mq2*wst(ptmeigen+sumT)))*(2D0*k2k3 - 2D0*kbk3 + m2**2 + m4**2 - 1D0*mq2c**2 - 2D0*pbp2 + (cI*mq2c*wst(ptmeigen+sumTc)))))
                         

            NeuQ2qx_M2RadiationVTT = NeuQ2qx_M2RadiationVTT + CF*NC*(radvTiTi + 2.*radvTiTp + 2.*radvTiTf + 2.*radvTiTv + radvTpTp + 2.*radvTpTf + 2.*radvTpTv + radvTfTf + 2.*radvTfTv + radvTvTv)
   
            
        enddo
        enddo
        
    
        end

! =================================================================
        double precision function NeuQ2qx_M2RadiationVTU()
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double precision radvTiUi, radvTiUf, radvTpUi, radvTpUf, radvTfUi, radvTfUf, radvTvUi, radvTvUf
        double complex A2L, A2R, A3Lc, A3Rc, C2L, D3Lc, D3Rc, CLL2
        double precision mq2, mq3c
        double precision xv, colorfactor
        integer sumT, sumUc, sumUend

         if (xtype.eq.1) then
        
           xv = 0d0
           colorfactor = 4d0/3d0
           sumUend = -1
           
        elseif (xtype.eq.2) then
        
           xv = 0d0
           colorfactor = 1d0
           sumUend = -1
           
        elseif (xtype.eq.3) then
        
           xv = 1/(m4**2)
           colorfactor = 1d0
           sumUend = 3
    
           
        elseif (xtype.eq.4) then

           xv = 1/(m4**2)
           colorfactor = 1d0
           sumUend = 1
                   
        endif

        NeuQ2qx_M2RadiationVTU = 0d0

        do sumT = 0, 1
        do sumUc = 0, sumUend
    
        radvTiUi = 0d0
        radvTiUf = 0d0
        radvTpUi = 0d0
        radvTpUf = 0d0
        radvTfUi = 0d0
        radvTfUf = 0d0
        radvTvUi = 0d0
        radvTvUf = 0d0

        A2L = AL(chimeigen,ptmeigen+sumT,pttype,ptgen)
        A2R = AR(chimeigen,ptmeigen+sumT,pttype,ptgen)
 
        C2L = CL(ptmeigen+sumT,sfmeigen,sftype,sfgen,xtype)

        CLL2 =  CLL(ptmeigen+sumT,sfmeigen,sftype,sfgen,xtype)
 
        mq2 = MSf(ptmeigen+sumT,pttype,ptgen)

        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        D3Lc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        D3Rc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3c = MNeuCha(pumeigen+sumUc)

        radvTiUi = dreal((4D0*C2L*GS**2*(kbk3 - 1D0*m2**2)*(A2L*(A3Rc*m3*(D3Rc*m1*(-2D0*k1k2 - 2D0*m1**2 - 1D0*m4**2 + 2D0*pap1 + 3D0*pap2 + 
     &                (2D0*k1k2 + m4**2 - 2D0*pap2)*(m4**2 - 1D0*pap2)*xv) + D3Lc*mq3c*(-2D0*m1**2 + 2D0*pap1 + pap2 + pap2*(-2D0*k1k2 - 1D0*m4**2 + 2D0*pap2)*xv)) + 
     &          A3Lc*(D3Rc*m1*mq3c*(k1k2 + 2D0*m3**2 - 2D0*pap1 - 1D0*k1k2*(2D0*k1k2 + m4**2 - 2D0*pap2)*xv) + 
     &             D3Lc*(k1k2*m1**2 - 2D0*m1**2*m3**2 - 4D0*k1k2*pap1 - 2D0*m1**2*pap1 - 1D0*m4**2*pap1 + 4D0*pap1**2 + 2D0*m3**2*pap2 + 2D0*pap1*pap2 + 
     &                (2D0*k1k2 + m4**2 - 2D0*pap2)*(k1k2*m1**2 + m4**2*pap1 - 2D0*pap1*pap2)*xv))) + 
     &       A2R*(A3Lc*m3*(D3Lc*m1*(-2D0*k1k2 - 2D0*m1**2 - 1D0*m4**2 + 2D0*pap1 + 3D0*pap2 + (2D0*k1k2 + m4**2 - 2D0*pap2)*(m4**2 - 1D0*pap2)*xv) + 
     &             D3Rc*mq3c*(-2D0*m1**2 + 2D0*pap1 + pap2 + pap2*(-2D0*k1k2 - 1D0*m4**2 + 2D0*pap2)*xv)) + 
     &          A3Rc*(D3Lc*m1*mq3c*(k1k2 + 2D0*m3**2 - 2D0*pap1 - 1D0*k1k2*(2D0*k1k2 + m4**2 - 2D0*pap2)*xv) + 
     &             D3Rc*(k1k2*m1**2 - 2D0*m1**2*m3**2 - 4D0*k1k2*pap1 - 2D0*m1**2*pap1 - 1D0*m4**2*pap1 + 4D0*pap1**2 + 2D0*m3**2*pap2 + 2D0*pap1*pap2 + 
     &                (2D0*k1k2 + m4**2 - 2D0*pap2)*(k1k2*m1**2 + m4**2*pap1 - 2D0*pap1*pap2)*xv)))))/
     &   (kbk3*(m1**2 + m3**2 - 1D0*mq2**2 - 2D0*pap1 - (cI*mq2*wst(ptmeigen+sumT)))*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)*(2D0*k1k2 + m1**2 - 1D0*m2**2 + m3**2 + m4**2 - 2D0*(pap1 + pap2))))


        radvTiUf = dreal((2D0*C2L*GS**2*(A2L*(A3Rc*m3*(D3Rc*m1*(2D0*k1k2*k1k3 - 2D0*k1k3*kakb + k2k3*kakb - 2D0*k1k2*kbk3 + 2D0*k1k3*m1**2 - 2D0*kbk3*m1**2 + k1k3*m4**2 - 1D0*kbk3*m4**2 - 
     &                2D0*k1k3*pap1 + 2D0*kbk3*pap1 - 3D0*k1k3*pap2 + 3D0*kbk3*pap2 - 4D0*k1k2*pbp1 - 2D0*k2k3*pbp1 + 2D0*kak3*pbp1 - 4D0*m1**2*pbp1 - 2D0*m4**2*pbp1 + 
     &                4D0*pap1*pbp1 + 6D0*pap2*pbp1 + 2D0*k1k3*pbp2 - 1D0*kak3*pbp2 + 
     &                (2D0*k1k2 + m4**2 - 2D0*pap2)*(k2k3*kakb - 1D0*(m4**2 - 1D0*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*kak3*pbp2)*xv) + 
     &             D3Lc*mq3c*(2D0*k1k3*kakb + k2k3*kakb + 2D0*k1k3*m1**2 - 2D0*kbk3*m1**2 - 2D0*k1k3*pap1 + 2D0*kbk3*pap1 - 1D0*k1k3*pap2 + kbk3*pap2 - 4D0*m1**2*pbp1 + 
     &                4D0*pap1*pbp1 + 2D0*pap2*pbp1 - 1D0*kak3*(2D0*pbp1 + pbp2) + (2D0*k1k2 + m4**2 - 2D0*pap2)*(-1D0*k2k3*kakb + pap2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + kak3*pbp2)*xv
     &                )) + A3Lc*(D3Rc*m1*mq3c*(-1D0*k1k2*k1k3 + 2D0*k1k3*kakb + k1k2*kbk3 - 2D0*k1k3*m3**2 + 2D0*kbk3*m3**2 + 2D0*k1k3*pap1 - 2D0*kbk3*pap1 + 2D0*k1k2*pbp1 + 
     &                k2k3*pbp1 - 2D0*kak3*pbp1 + 4D0*m3**2*pbp1 - 4D0*pap1*pbp1 - 1D0*k1k3*pbp2 + 
     &                (2D0*k1k2 + m4**2 - 2D0*pap2)*(k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*k2k3*pbp1 + k1k3*pbp2)*xv) + 
     &             D3Lc*(-1D0*k1k2*k1k3*m1**2 + 2D0*k1k3*kakb*m1**2 + k1k2*kbk3*m1**2 + 2D0*k2k3*kakb*m3**2 + 2D0*k1k3*m1**2*m3**2 - 2D0*kbk3*m1**2*m3**2 + k1k3*kakb*m4**2 + 
     &                4D0*k1k2*k1k3*pap1 - 4D0*k1k3*kakb*pap1 - 4D0*k1k2*kbk3*pap1 + 2D0*k1k3*m1**2*pap1 - 2D0*kbk3*m1**2*pap1 + k1k3*m4**2*pap1 - 1D0*kbk3*m4**2*pap1 - 
     &                4D0*k1k3*pap1**2 + 4D0*kbk3*pap1**2 - 2D0*k1k3*kakb*pap2 - 2D0*k1k3*m3**2*pap2 + 2D0*kbk3*m3**2*pap2 - 2D0*k1k3*pap1*pap2 + 2D0*kbk3*pap1*pap2 + 
     &                2D0*k1k2*m1**2*pbp1 + k2k3*m1**2*pbp1 - 2D0*kak3*m1**2*pbp1 - 4D0*m1**2*m3**2*pbp1 - 1D0*kak3*m4**2*pbp1 - 8D0*k1k2*pap1*pbp1 - 4D0*k2k3*pap1*pbp1 + 
     &                4D0*kak3*pap1*pbp1 - 4D0*m1**2*pap1*pbp1 - 2D0*m4**2*pap1*pbp1 + 8D0*pap1**2*pbp1 + 2D0*kak3*pap2*pbp1 + 4D0*m3**2*pap2*pbp1 + 4D0*pap1*pap2*pbp1 - 
     &                1D0*k1k3*m1**2*pbp2 - 2D0*kak3*m3**2*pbp2 + 4D0*k1k3*pap1*pbp2 - 
     &                1D0*(2D0*k1k2 + m4**2 - 2D0*pap2)*(k1k2*m1**2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*k2k3*m1**2*pbp1 + 
     &                   (m4**2 - 2D0*pap2)*(k1k3*(kakb + pap1) - 1D0*kak3*pbp1 - 1D0*pap1*(kbk3 + 2D0*pbp1)) + k1k3*m1**2*pbp2)*xv))) + 
     &       A2R*(A3Lc*m3*(D3Lc*m1*(2D0*k1k2*k1k3 - 2D0*k1k3*kakb + k2k3*kakb - 2D0*k1k2*kbk3 + 2D0*k1k3*m1**2 - 2D0*kbk3*m1**2 + k1k3*m4**2 - 1D0*kbk3*m4**2 - 2D0*k1k3*pap1 + 
     &                2D0*kbk3*pap1 - 3D0*k1k3*pap2 + 3D0*kbk3*pap2 - 4D0*k1k2*pbp1 - 2D0*k2k3*pbp1 + 2D0*kak3*pbp1 - 4D0*m1**2*pbp1 - 2D0*m4**2*pbp1 + 4D0*pap1*pbp1 + 
     &                6D0*pap2*pbp1 + 2D0*k1k3*pbp2 - 1D0*kak3*pbp2 + (2D0*k1k2 + m4**2 - 2D0*pap2)*
     &                 (k2k3*kakb - 1D0*(m4**2 - 1D0*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*kak3*pbp2)*xv) + 
     &             D3Rc*mq3c*(2D0*k1k3*kakb + k2k3*kakb + 2D0*k1k3*m1**2 - 2D0*kbk3*m1**2 - 2D0*k1k3*pap1 + 2D0*kbk3*pap1 - 1D0*k1k3*pap2 + kbk3*pap2 - 4D0*m1**2*pbp1 + 
     &                4D0*pap1*pbp1 + 2D0*pap2*pbp1 - 1D0*kak3*(2D0*pbp1 + pbp2) + (2D0*k1k2 + m4**2 - 2D0*pap2)*(-1D0*k2k3*kakb + pap2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + kak3*pbp2)*xv
     &                )) + A3Rc*(D3Lc*m1*mq3c*(-1D0*k1k2*k1k3 + 2D0*k1k3*kakb + k1k2*kbk3 - 2D0*k1k3*m3**2 + 2D0*kbk3*m3**2 + 2D0*k1k3*pap1 - 2D0*kbk3*pap1 + 2D0*k1k2*pbp1 + 
     &                k2k3*pbp1 - 2D0*kak3*pbp1 + 4D0*m3**2*pbp1 - 4D0*pap1*pbp1 - 1D0*k1k3*pbp2 + 
     &                (2D0*k1k2 + m4**2 - 2D0*pap2)*(k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*k2k3*pbp1 + k1k3*pbp2)*xv) + 
     &             D3Rc*(-1D0*k1k2*k1k3*m1**2 + 2D0*k1k3*kakb*m1**2 + k1k2*kbk3*m1**2 + 2D0*k2k3*kakb*m3**2 + 2D0*k1k3*m1**2*m3**2 - 2D0*kbk3*m1**2*m3**2 + k1k3*kakb*m4**2 + 
     &                4D0*k1k2*k1k3*pap1 - 4D0*k1k3*kakb*pap1 - 4D0*k1k2*kbk3*pap1 + 2D0*k1k3*m1**2*pap1 - 2D0*kbk3*m1**2*pap1 + k1k3*m4**2*pap1 - 1D0*kbk3*m4**2*pap1 - 
     &                4D0*k1k3*pap1**2 + 4D0*kbk3*pap1**2 - 2D0*k1k3*kakb*pap2 - 2D0*k1k3*m3**2*pap2 + 2D0*kbk3*m3**2*pap2 - 2D0*k1k3*pap1*pap2 + 2D0*kbk3*pap1*pap2 + 
     &                2D0*k1k2*m1**2*pbp1 + k2k3*m1**2*pbp1 - 2D0*kak3*m1**2*pbp1 - 4D0*m1**2*m3**2*pbp1 - 1D0*kak3*m4**2*pbp1 - 8D0*k1k2*pap1*pbp1 - 4D0*k2k3*pap1*pbp1 + 
     &                4D0*kak3*pap1*pbp1 - 4D0*m1**2*pap1*pbp1 - 2D0*m4**2*pap1*pbp1 + 8D0*pap1**2*pbp1 + 2D0*kak3*pap2*pbp1 + 4D0*m3**2*pap2*pbp1 + 4D0*pap1*pap2*pbp1 - 
     &                1D0*k1k3*m1**2*pbp2 - 2D0*kak3*m3**2*pbp2 + 4D0*k1k3*pap1*pbp2 - 
     &                1D0*(2D0*k1k2 + m4**2 - 2D0*pap2)*(k1k2*m1**2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*k2k3*m1**2*pbp1 + 
     &                   (m4**2 - 2D0*pap2)*(k1k3*(kakb + pap1) - 1D0*kak3*pbp1 - 1D0*pap1*(kbk3 + 2D0*pbp1)) + k1k3*m1**2*pbp2)*xv)))))/
     &   (k1k3*(m1**2 + m3**2 - 1D0*mq2**2 - 2D0*pap1 - (cI*mq2*wst(ptmeigen+sumT)))*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)*(2D0*k1k2 + m1**2 - 1D0*m2**2 + m3**2 + m4**2 - 2D0*(pap1 + pap2))))


        radvTpUi = dreal((C2L*GS**2*(k1k3 - 1D0*k2k3 - 1D0*kak3 + 2D0*kakb + kbk3 - 2D0*(m2**2 + pbp1 - 1D0*pbp2))*
     &     (A2L*(A3Lc*(D3Lc*(2D0*k1k2*kakb - 1D0*k1k2*m1**2 - 4D0*kakb*pap1 - 1D0*m4**2*pap1 + 2D0*pap1*pap2 + 2D0*m1**2*pbp1 - 2D0*pap2*pbp1 + 2D0*pap1*pbp2 + 
     &                (k1k2*m1**2 + m4**2*pap1 - 2D0*pap1*pap2)*(m4**2 - 2D0*pbp2)*xv) + D3Rc*m1*mq3c*(k1k2 - 2D0*pbp1 - 1D0*k1k2*m4**2*xv + 2D0*k1k2*pbp2*xv)) + 
     &          A3Rc*m3*(D3Rc*m1*(-2D0*kakb - 1D0*m4**2 + pap2 + 2D0*pbp2 + (m4**2 - 1D0*pap2)*(m4**2 - 2D0*pbp2)*xv) + 
     &             D3Lc*mq3c*(-2D0*kakb + pap2 - 1D0*m4**2*pap2*xv + 2D0*pap2*pbp2*xv))) + 
     &       A2R*(A3Rc*(D3Rc*(2D0*k1k2*kakb - 1D0*k1k2*m1**2 - 4D0*kakb*pap1 - 1D0*m4**2*pap1 + 2D0*pap1*pap2 + 2D0*m1**2*pbp1 - 2D0*pap2*pbp1 + 2D0*pap1*pbp2 + 
     &                (k1k2*m1**2 + m4**2*pap1 - 2D0*pap1*pap2)*(m4**2 - 2D0*pbp2)*xv) + D3Lc*m1*mq3c*(k1k2 - 2D0*pbp1 - 1D0*k1k2*m4**2*xv + 2D0*k1k2*pbp2*xv)) + 
     &          A3Lc*m3*(D3Lc*m1*(-2D0*kakb - 1D0*m4**2 + pap2 + 2D0*pbp2 + (m4**2 - 1D0*pap2)*(m4**2 - 2D0*pbp2)*xv) + 
     &             D3Rc*mq3c*(-2D0*kakb + pap2 - 1D0*m4**2*pap2*xv + 2D0*pap2*pbp2*xv)))))/
     &   (kbk3*(m1**2 + m3**2 - 1D0*mq2**2 - 2D0*pap1 - (cI*mq2*wst(ptmeigen+sumT)))*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)*(-1D0*m2**2 - 1D0*m4**2 + mq2**2 + 2D0*pbp2 + (cI*mq2*wst(ptmeigen+sumT)))))


        radvTpUf = dreal((C2L*GS**2*(A2L*(A3Rc*m3*(D3Lc*mq3c*(-1D0*k1k2*kak3 + 4D0*k1k2*kakb - 2D0*k1k3*kakb + 3D0*k2k3*kakb - 4D0*kakb*kbk3 - 1D0*k2k3*m1**2 + 2D0*kbk3*m1**2 + 
     &                2D0*kak3*m2**2 - 4D0*kakb*m3**2 + kak3*m4**2 + k2k3*pap1 + 4D0*kakb*pap1 - 2D0*kbk3*pap1 - 2D0*k1k2*pap2 + k1k3*pap2 - 2D0*k2k3*pap2 + 3D0*kbk3*pap2 + 
     &                2D0*m3**2*pap2 - 2D0*pap1*pap2 + 2D0*kak3*pbp1 - 4D0*kakb*pbp1 + 2D0*pap2*pbp1 - 3D0*kak3*pbp2 - 
     &                1D0*(k2k3*(kakb - 1D0*m1**2 + pap1 - 2D0*pap2) - 1D0*k1k2*(kak3 + 2D0*pap2) + pap2*(k1k3 + kbk3 + 2D0*(m3**2 - 1D0*pap1 + pbp1)) + kak3*(m4**2 - 1D0*pbp2))*
     &                 (m4**2 - 2D0*pbp2)*xv) + D3Rc*m1*(k1k2*kak3 + 4D0*k1k2*kakb - 2D0*k1k3*kakb - 1D0*k2k3*kakb + 4D0*kak3*kakb - 2D0*k1k2*kbk3 + k2k3*m1**2 - 2D0*kbk3*m1**2 + 
     &                2D0*k2k3*m2**2 - 2D0*kak3*m2**2 - 4D0*kakb*m3**2 + 2D0*k1k2*m4**2 - 1D0*k1k3*m4**2 + k2k3*m4**2 + kbk3*m4**2 - 2D0*m3**2*m4**2 - 1D0*k2k3*pap1 + 
     &                4D0*kakb*pap1 + 2D0*kbk3*pap1 + 2D0*m4**2*pap1 - 2D0*k1k2*pap2 + k1k3*pap2 - 2D0*kak3*pap2 + kbk3*pap2 + 2D0*m3**2*pap2 - 2D0*pap1*pap2 + 2D0*k2k3*pbp1 - 
     &                2D0*kak3*pbp1 - 4D0*kakb*pbp1 - 2D0*m4**2*pbp1 + 2D0*pap2*pbp1 - 4D0*k1k2*pbp2 + 2D0*k1k3*pbp2 - 4D0*k2k3*pbp2 + kak3*pbp2 + 4D0*m3**2*pbp2 - 4D0*pap1*pbp2 + 
     &                4D0*pbp1*pbp2 - 1D0*(m4**2 - 2D0*pbp2)*(k2k3*(-1D0*kakb + m1**2 + m4**2 - 1D0*pap1) + k1k2*(kak3 + 2D0*m4**2 - 2D0*pap2) + k1k3*(-1D0*m4**2 + pap2) - 
     &                   1D0*(m4**2 - 1D0*pap2)*(kbk3 + 2D0*(m3**2 - 1D0*pap1 + pbp1)) + kak3*(-2D0*pap2 + pbp2))*xv)) + 
     &          A3Lc*(D3Rc*m1*mq3c*(-2D0*k1k2**2 - 2D0*k1k2*k2k3 - 1D0*k1k2*kak3 + 3D0*k1k2*kbk3 + 2D0*k1k2*m3**2 + k2k3*m3**2 - 2D0*kbk3*m3**2 - 2D0*k1k2*pap1 - 1D0*k2k3*pap1 + 
     &                2D0*kbk3*pap1 + 6D0*k1k2*pbp1 + 3D0*k2k3*pbp1 + 2D0*kak3*pbp1 - 4D0*kbk3*pbp1 - 4D0*m3**2*pbp1 + 4D0*pap1*pbp1 - 4D0*pbp1**2 + 
     &                k1k3*(-2D0*kakb + 2D0*m2**2 + m4**2 + pap2 - 3D0*pbp2) - 
     &                1D0*(-2D0*k1k2**2 + k2k3*(m3**2 - 1D0*pap1 + pbp1) + k1k2*(-2D0*k2k3 - 1D0*kak3 + kbk3 + 2D0*(m3**2 - 1D0*pap1 + pbp1)) + k1k3*(m4**2 + pap2 - 1D0*pbp2))*
     &                 (m4**2 - 2D0*pbp2)*xv) + D3Lc*(2D0*k2k3*kakb*m3**2 - 4D0*kak3*kakb*m3**2 - 1D0*k2k3*m1**2*m3**2 + 2D0*kbk3*m1**2*m3**2 - 1D0*kak3*m3**2*m4**2 + 
     &                8D0*kak3*kakb*pap1 - 4D0*kakb*kbk3*pap1 + k2k3*m1**2*pap1 - 2D0*kbk3*m1**2*pap1 + 2D0*k2k3*m2**2*pap1 - 8D0*kakb*m3**2*pap1 + k2k3*m4**2*pap1 + 
     &                2D0*kak3*m4**2*pap1 + kbk3*m4**2*pap1 - 2D0*m3**2*m4**2*pap1 + 8D0*kakb*pap1**2 + 2D0*m4**2*pap1**2 + 2D0*kak3*m3**2*pap2 - 2D0*kbk3*m3**2*pap2 - 
     &                2D0*k2k3*pap1*pap2 - 4D0*kak3*pap1*pap2 + 6D0*kbk3*pap1*pap2 + 4D0*m3**2*pap1*pap2 - 4D0*pap1**2*pap2 - 4D0*kak3*kakb*pbp1 - 1D0*k2k3*m1**2*pbp1 - 
     &                2D0*kak3*m1**2*pbp1 + 4D0*kbk3*m1**2*pbp1 + 4D0*m1**2*m3**2*pbp1 - 3D0*kak3*m4**2*pbp1 - 8D0*kakb*pap1*pbp1 - 4D0*m1**2*pap1*pbp1 - 2D0*m4**2*pap1*pbp1 + 
     &                4D0*k2k3*pap2*pbp1 + 2D0*kak3*pap2*pbp1 - 4D0*kbk3*pap2*pbp1 - 4D0*m3**2*pap2*pbp1 + 8D0*pap1*pap2*pbp1 + 4D0*m1**2*pbp1**2 - 4D0*pap2*pbp1**2 + 
     &                2D0*kak3*m3**2*pbp2 - 4D0*k2k3*pap1*pbp2 - 4D0*kak3*pap1*pbp2 + 4D0*m3**2*pap1*pbp2 - 4D0*pap1**2*pbp2 + 4D0*kak3*pbp1*pbp2 + 4D0*pap1*pbp1*pbp2 + 
     &                (k2k3*(-1D0*pap1*(m4**2 - 2D0*pap2) + m1**2*(m3**2 - 1D0*pap1 + pbp1)) + 
     &                   (m4**2 - 2D0*pap2)*(kak3*(m3**2 - 2D0*pap1 + pbp1) + pap1*(kbk3 + 2D0*(m3**2 - 1D0*pap1 + pbp1))))*(m4**2 - 2D0*pbp2)*xv + 
     &                k1k2*(-4D0*k2k3*kakb + 4D0*kak3*kakb + 4D0*kakb*kbk3 + 2D0*k2k3*m1**2 + kak3*m1**2 - 5D0*kbk3*m1**2 - 2D0*kak3*m2**2 + 4D0*kakb*m3**2 - 2D0*m1**2*m3**2 + 
     &                   kak3*m4**2 + 4D0*kakb*pap1 + 2D0*m1**2*pap1 + 2D0*m4**2*pap1 - 2D0*kak3*pap2 - 4D0*pap1*pap2 + 4D0*kakb*pbp1 - 6D0*m1**2*pbp1 + 4D0*pap2*pbp1 - 
     &                   4D0*pap1*pbp2 - 1D0*(2D0*k2k3*m1**2 - 1D0*kbk3*m1**2 - 2D0*m1**2*m3**2 + 2D0*m1**2*pap1 + 2D0*m4**2*pap1 + kak3*(m1**2 + m4**2 - 2D0*pap2) - 
     &                      4D0*pap1*pap2 - 2D0*m1**2*pbp1)*(m4**2 - 2D0*pbp2)*xv) + 2D0*k1k2**2*(-2D0*kakb + m1**2*(1D0 - 1D0*m4**2*xv + 2D0*pbp2*xv)) + 
     &                k1k3*(4D0*kakb**2 + kakb*(-2D0*m1**2 + 3D0*m4**2 - 6D0*pap2 - 4D0*pbp2 - 1D0*(m4**2 - 2D0*pap2)*(m4**2 - 2D0*pbp2)*xv) + 
     &                   m1**2*(-2D0*m2**2 - 2D0*m4**2 + pap2 + 5D0*pbp2 + (-2D0*m4**2 + pap2 + pbp2)*(-1D0*m4**2 + 2D0*pbp2)*xv) + 
     &                   pap2*(2D0*m2**2 + (m4**2 - 2D0*pap2)*(-1D0 + m4**2*xv - 2D0*pbp2*xv)))))) + 
     &       A2R*(A3Lc*m3*(D3Rc*mq3c*(-1D0*k1k2*kak3 + 4D0*k1k2*kakb - 2D0*k1k3*kakb + 3D0*k2k3*kakb - 4D0*kakb*kbk3 - 1D0*k2k3*m1**2 + 2D0*kbk3*m1**2 + 2D0*kak3*m2**2 - 
     &                4D0*kakb*m3**2 + kak3*m4**2 + k2k3*pap1 + 4D0*kakb*pap1 - 2D0*kbk3*pap1 - 2D0*k1k2*pap2 + k1k3*pap2 - 2D0*k2k3*pap2 + 3D0*kbk3*pap2 + 2D0*m3**2*pap2 - 
     &                2D0*pap1*pap2 + 2D0*kak3*pbp1 - 4D0*kakb*pbp1 + 2D0*pap2*pbp1 - 3D0*kak3*pbp2 - 
     &                1D0*(k2k3*(kakb - 1D0*m1**2 + pap1 - 2D0*pap2) - 1D0*k1k2*(kak3 + 2D0*pap2) + pap2*(k1k3 + kbk3 + 2D0*(m3**2 - 1D0*pap1 + pbp1)) + kak3*(m4**2 - 1D0*pbp2))*
     &                 (m4**2 - 2D0*pbp2)*xv) + D3Lc*m1*(k1k2*kak3 + 4D0*k1k2*kakb - 2D0*k1k3*kakb - 1D0*k2k3*kakb + 4D0*kak3*kakb - 2D0*k1k2*kbk3 + k2k3*m1**2 - 2D0*kbk3*m1**2 + 
     &                2D0*k2k3*m2**2 - 2D0*kak3*m2**2 - 4D0*kakb*m3**2 + 2D0*k1k2*m4**2 - 1D0*k1k3*m4**2 + k2k3*m4**2 + kbk3*m4**2 - 2D0*m3**2*m4**2 - 1D0*k2k3*pap1 + 
     &                4D0*kakb*pap1 + 2D0*kbk3*pap1 + 2D0*m4**2*pap1 - 2D0*k1k2*pap2 + k1k3*pap2 - 2D0*kak3*pap2 + kbk3*pap2 + 2D0*m3**2*pap2 - 2D0*pap1*pap2 + 2D0*k2k3*pbp1 - 
     &                2D0*kak3*pbp1 - 4D0*kakb*pbp1 - 2D0*m4**2*pbp1 + 2D0*pap2*pbp1 - 4D0*k1k2*pbp2 + 2D0*k1k3*pbp2 - 4D0*k2k3*pbp2 + kak3*pbp2 + 4D0*m3**2*pbp2 - 4D0*pap1*pbp2 + 
     &                4D0*pbp1*pbp2 - 1D0*(m4**2 - 2D0*pbp2)*(k2k3*(-1D0*kakb + m1**2 + m4**2 - 1D0*pap1) + k1k2*(kak3 + 2D0*m4**2 - 2D0*pap2) + k1k3*(-1D0*m4**2 + pap2) - 
     &                   1D0*(m4**2 - 1D0*pap2)*(kbk3 + 2D0*(m3**2 - 1D0*pap1 + pbp1)) + kak3*(-2D0*pap2 + pbp2))*xv)) + 
     &          A3Rc*(D3Lc*m1*mq3c*(-2D0*k1k2**2 - 2D0*k1k2*k2k3 - 1D0*k1k2*kak3 + 3D0*k1k2*kbk3 + 2D0*k1k2*m3**2 + k2k3*m3**2 - 2D0*kbk3*m3**2 - 2D0*k1k2*pap1 - 1D0*k2k3*pap1 + 
     &                2D0*kbk3*pap1 + 6D0*k1k2*pbp1 + 3D0*k2k3*pbp1 + 2D0*kak3*pbp1 - 4D0*kbk3*pbp1 - 4D0*m3**2*pbp1 + 4D0*pap1*pbp1 - 4D0*pbp1**2 + 
     &                k1k3*(-2D0*kakb + 2D0*m2**2 + m4**2 + pap2 - 3D0*pbp2) - 
     &                1D0*(-2D0*k1k2**2 + k2k3*(m3**2 - 1D0*pap1 + pbp1) + k1k2*(-2D0*k2k3 - 1D0*kak3 + kbk3 + 2D0*(m3**2 - 1D0*pap1 + pbp1)) + k1k3*(m4**2 + pap2 - 1D0*pbp2))*
     &                 (m4**2 - 2D0*pbp2)*xv) + D3Rc*(2D0*k2k3*kakb*m3**2 - 4D0*kak3*kakb*m3**2 - 1D0*k2k3*m1**2*m3**2 + 2D0*kbk3*m1**2*m3**2 - 1D0*kak3*m3**2*m4**2 + 
     &                8D0*kak3*kakb*pap1 - 4D0*kakb*kbk3*pap1 + k2k3*m1**2*pap1 - 2D0*kbk3*m1**2*pap1 + 2D0*k2k3*m2**2*pap1 - 8D0*kakb*m3**2*pap1 + k2k3*m4**2*pap1 + 
     &                2D0*kak3*m4**2*pap1 + kbk3*m4**2*pap1 - 2D0*m3**2*m4**2*pap1 + 8D0*kakb*pap1**2 + 2D0*m4**2*pap1**2 + 2D0*kak3*m3**2*pap2 - 2D0*kbk3*m3**2*pap2 - 
     &                2D0*k2k3*pap1*pap2 - 4D0*kak3*pap1*pap2 + 6D0*kbk3*pap1*pap2 + 4D0*m3**2*pap1*pap2 - 4D0*pap1**2*pap2 - 4D0*kak3*kakb*pbp1 - 1D0*k2k3*m1**2*pbp1 - 
     &                2D0*kak3*m1**2*pbp1 + 4D0*kbk3*m1**2*pbp1 + 4D0*m1**2*m3**2*pbp1 - 3D0*kak3*m4**2*pbp1 - 8D0*kakb*pap1*pbp1 - 4D0*m1**2*pap1*pbp1 - 2D0*m4**2*pap1*pbp1 + 
     &                4D0*k2k3*pap2*pbp1 + 2D0*kak3*pap2*pbp1 - 4D0*kbk3*pap2*pbp1 - 4D0*m3**2*pap2*pbp1 + 8D0*pap1*pap2*pbp1 + 4D0*m1**2*pbp1**2 - 4D0*pap2*pbp1**2 + 
     &                2D0*kak3*m3**2*pbp2 - 4D0*k2k3*pap1*pbp2 - 4D0*kak3*pap1*pbp2 + 4D0*m3**2*pap1*pbp2 - 4D0*pap1**2*pbp2 + 4D0*kak3*pbp1*pbp2 + 4D0*pap1*pbp1*pbp2 + 
     &                (k2k3*(-1D0*pap1*(m4**2 - 2D0*pap2) + m1**2*(m3**2 - 1D0*pap1 + pbp1)) + 
     &                   (m4**2 - 2D0*pap2)*(kak3*(m3**2 - 2D0*pap1 + pbp1) + pap1*(kbk3 + 2D0*(m3**2 - 1D0*pap1 + pbp1))))*(m4**2 - 2D0*pbp2)*xv + 
     &                k1k2*(-4D0*k2k3*kakb + 4D0*kak3*kakb + 4D0*kakb*kbk3 + 2D0*k2k3*m1**2 + kak3*m1**2 - 5D0*kbk3*m1**2 - 2D0*kak3*m2**2 + 4D0*kakb*m3**2 - 2D0*m1**2*m3**2 + 
     &                   kak3*m4**2 + 4D0*kakb*pap1 + 2D0*m1**2*pap1 + 2D0*m4**2*pap1 - 2D0*kak3*pap2 - 4D0*pap1*pap2 + 4D0*kakb*pbp1 - 6D0*m1**2*pbp1 + 4D0*pap2*pbp1 - 
     &                   4D0*pap1*pbp2 - 1D0*(2D0*k2k3*m1**2 - 1D0*kbk3*m1**2 - 2D0*m1**2*m3**2 + 2D0*m1**2*pap1 + 2D0*m4**2*pap1 + kak3*(m1**2 + m4**2 - 2D0*pap2) - 
     &                      4D0*pap1*pap2 - 2D0*m1**2*pbp1)*(m4**2 - 2D0*pbp2)*xv) + 2D0*k1k2**2*(-2D0*kakb + m1**2*(1D0 - 1D0*m4**2*xv + 2D0*pbp2*xv)) + 
     &                k1k3*(4D0*kakb**2 + kakb*(-2D0*m1**2 + 3D0*m4**2 - 6D0*pap2 - 4D0*pbp2 - 1D0*(m4**2 - 2D0*pap2)*(m4**2 - 2D0*pbp2)*xv) + 
     &                   m1**2*(-2D0*m2**2 - 2D0*m4**2 + pap2 + 5D0*pbp2 + (-2D0*m4**2 + pap2 + pbp2)*(-1D0*m4**2 + 2D0*pbp2)*xv) + 
     &                   pap2*(2D0*m2**2 + (m4**2 - 2D0*pap2)*(-1D0 + m4**2*xv - 2D0*pbp2*xv))))))))/
     &   (k1k3*(m1**2 + m3**2 - 1D0*mq2**2 - 2D0*pap1 - (cI*mq2*wst(ptmeigen+sumT)))*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)*(-1D0*m2**2 - 1D0*m4**2 + mq2**2 + 2D0*pbp2 + (cI*mq2*wst(ptmeigen+sumT)))))


        radvTfUi = dreal((C2L*GS**2*(A2L*(A3Rc*m3*(D3Rc*m1*(2D0*k1k3*kakb + k2k3*kakb - 4D0*kakb*kbk3 - 2D0*k2k3*m2**2 + 2D0*kak3*m2**2 + k1k3*m4**2 - 1D0*kbk3*m4**2 - 1D0*k1k3*pap2 + 
     &                kbk3*pap2 - 4D0*kakb*pbp1 - 2D0*m4**2*pbp1 + 2D0*pap2*pbp1 - 2D0*k1k3*pbp2 - 1D0*kak3*pbp2 + 4D0*kbk3*pbp2 + 4D0*pbp1*pbp2 - 
     &                1D0*(m4**2 - 2D0*pbp2)*(k2k3*kakb + (m4**2 - 1D0*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*kak3*pbp2)*xv) + 
     &             D3Lc*mq3c*(2D0*k1k3*kakb - 1D0*k2k3*kakb - 2D0*kak3*m2**2 - 1D0*k1k3*pap2 + kbk3*pap2 - 4D0*kakb*pbp1 + 2D0*pap2*pbp1 + kak3*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*(k2k3*kakb + pap2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*kak3*pbp2)*xv)) + 
     &          A3Lc*(D3Rc*m1*mq3c*(-1D0*k1k2*k1k3 + k1k2*kbk3 + 2D0*k1k3*m2**2 + 2D0*k1k2*pbp1 + 2D0*k1k3*pbp1 + k2k3*pbp1 - 4D0*kbk3*pbp1 - 4D0*pbp1**2 - 1D0*k1k3*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*(k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*k2k3*pbp1 + k1k3*pbp2)*xv) + 
     &             D3Lc*(-2D0*k1k2*k1k3*kakb + 4D0*k1k3*kakb**2 + k1k2*k1k3*m1**2 - 1D0*k1k2*kbk3*m1**2 + 2D0*k1k2*kak3*m2**2 - 2D0*k1k3*m1**2*m2**2 + k1k3*kakb*m4**2 + 
     &                4D0*k1k3*kakb*pap1 - 4D0*kakb*kbk3*pap1 - 2D0*k2k3*m2**2*pap1 + k1k3*m4**2*pap1 - 1D0*kbk3*m4**2*pap1 - 2D0*k1k3*kakb*pap2 + 2D0*k1k3*m2**2*pap2 - 
     &                2D0*k1k3*pap1*pap2 + 2D0*kbk3*pap1*pap2 + 4D0*k1k2*kakb*pbp1 + 4D0*k2k3*kakb*pbp1 - 4D0*kak3*kakb*pbp1 - 2D0*k1k2*m1**2*pbp1 - 2D0*k1k3*m1**2*pbp1 - 
     &                1D0*k2k3*m1**2*pbp1 + 4D0*kbk3*m1**2*pbp1 - 1D0*kak3*m4**2*pbp1 - 8D0*kakb*pap1*pbp1 - 2D0*m4**2*pap1*pbp1 + 2D0*k1k3*pap2*pbp1 + 2D0*kak3*pap2*pbp1 - 
     &                4D0*kbk3*pap2*pbp1 + 4D0*pap1*pap2*pbp1 + 4D0*m1**2*pbp1**2 - 4D0*pap2*pbp1**2 - 4D0*k1k3*kakb*pbp2 + k1k3*m1**2*pbp2 - 2D0*k1k3*pap1*pbp2 + 
     &                4D0*kbk3*pap1*pbp2 + 4D0*pap1*pbp1*pbp2 - 1D0*(m4**2 - 2D0*pbp2)*
     &                 (k1k2*m1**2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*k2k3*m1**2*pbp1 + (m4**2 - 2D0*pap2)*(k1k3*(kakb + pap1) - 1D0*kak3*pbp1 - 1D0*pap1*(kbk3 + 2D0*pbp1)) + 
     &                   k1k3*m1**2*pbp2)*xv))) + A2R*(A3Lc*m3*(D3Lc*m1*(2D0*k1k3*kakb + k2k3*kakb - 4D0*kakb*kbk3 - 2D0*k2k3*m2**2 + 2D0*kak3*m2**2 + k1k3*m4**2 - 
     &                1D0*kbk3*m4**2 - 1D0*k1k3*pap2 + kbk3*pap2 - 4D0*kakb*pbp1 - 2D0*m4**2*pbp1 + 2D0*pap2*pbp1 - 2D0*k1k3*pbp2 - 1D0*kak3*pbp2 + 4D0*kbk3*pbp2 + 4D0*pbp1*pbp2 - 
     &                1D0*(m4**2 - 2D0*pbp2)*(k2k3*kakb + (m4**2 - 1D0*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*kak3*pbp2)*xv) + 
     &             D3Rc*mq3c*(2D0*k1k3*kakb - 1D0*k2k3*kakb - 2D0*kak3*m2**2 - 1D0*k1k3*pap2 + kbk3*pap2 - 4D0*kakb*pbp1 + 2D0*pap2*pbp1 + kak3*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*(k2k3*kakb + pap2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*kak3*pbp2)*xv)) + 
     &          A3Rc*(D3Lc*m1*mq3c*(-1D0*k1k2*k1k3 + k1k2*kbk3 + 2D0*k1k3*m2**2 + 2D0*k1k2*pbp1 + 2D0*k1k3*pbp1 + k2k3*pbp1 - 4D0*kbk3*pbp1 - 4D0*pbp1**2 - 1D0*k1k3*pbp2 + 
     &                (m4**2 - 2D0*pbp2)*(k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*k2k3*pbp1 + k1k3*pbp2)*xv) + 
     &             D3Rc*(-2D0*k1k2*k1k3*kakb + 4D0*k1k3*kakb**2 + k1k2*k1k3*m1**2 - 1D0*k1k2*kbk3*m1**2 + 2D0*k1k2*kak3*m2**2 - 2D0*k1k3*m1**2*m2**2 + k1k3*kakb*m4**2 + 
     &                4D0*k1k3*kakb*pap1 - 4D0*kakb*kbk3*pap1 - 2D0*k2k3*m2**2*pap1 + k1k3*m4**2*pap1 - 1D0*kbk3*m4**2*pap1 - 2D0*k1k3*kakb*pap2 + 2D0*k1k3*m2**2*pap2 - 
     &                2D0*k1k3*pap1*pap2 + 2D0*kbk3*pap1*pap2 + 4D0*k1k2*kakb*pbp1 + 4D0*k2k3*kakb*pbp1 - 4D0*kak3*kakb*pbp1 - 2D0*k1k2*m1**2*pbp1 - 2D0*k1k3*m1**2*pbp1 - 
     &                1D0*k2k3*m1**2*pbp1 + 4D0*kbk3*m1**2*pbp1 - 1D0*kak3*m4**2*pbp1 - 8D0*kakb*pap1*pbp1 - 2D0*m4**2*pap1*pbp1 + 2D0*k1k3*pap2*pbp1 + 2D0*kak3*pap2*pbp1 - 
     &                4D0*kbk3*pap2*pbp1 + 4D0*pap1*pap2*pbp1 + 4D0*m1**2*pbp1**2 - 4D0*pap2*pbp1**2 - 4D0*k1k3*kakb*pbp2 + k1k3*m1**2*pbp2 - 2D0*k1k3*pap1*pbp2 + 
     &                4D0*kbk3*pap1*pbp2 + 4D0*pap1*pbp1*pbp2 - 1D0*(m4**2 - 2D0*pbp2)*
     &                 (k1k2*m1**2*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*k2k3*m1**2*pbp1 + (m4**2 - 2D0*pap2)*(k1k3*(kakb + pap1) - 1D0*kak3*pbp1 - 1D0*pap1*(kbk3 + 2D0*pbp1)) + 
     &                   k1k3*m1**2*pbp2)*xv)))))/(k1k3*kbk3*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)*(m2**2 + m4**2 - 1D0*mq2**2 - 2D0*pbp2 - (cI*mq2*wst(ptmeigen+sumT)))))


        radvTfUf = dreal((2D0*C2L*GS**2*(A2R*(A3Lc*m3*(k1k3 + m3**2)*(D3Lc*m1*(2D0*kakb + m4**2 - 1D0*pap2 - 2D0*pbp2 - 1D0*(m4**2 - 1D0*pap2)*(m4**2 - 2D0*pbp2)*xv) + 
     &             D3Rc*mq3c*(2D0*kakb - 1D0*pap2 + pap2*(m4**2 - 2D0*pbp2)*xv)) + 
     &          A3Rc*(D3Rc*(m3**2*(-2D0*k1k2*kakb - 2D0*k2k3*kakb + 4D0*kak3*kakb + k1k2*m1**2 + k2k3*m1**2 - 2D0*kbk3*m1**2 + kak3*m4**2 + 4D0*kakb*pap1 + m4**2*pap1 - 
     &                   2D0*kak3*pap2 + 2D0*kbk3*pap2 - 2D0*pap1*pap2 - 2D0*m1**2*pbp1 + 2D0*pap2*pbp1 - 2D0*kak3*pbp2 - 2D0*pap1*pbp2 - 
     &                   1D0*((k1k2 + k2k3)*m1**2 + (kak3 + pap1)*(m4**2 - 2D0*pap2))*(m4**2 - 2D0*pbp2)*xv) + 
     &                k1k3*(2D0*k2k3*kakb - 4D0*kak3*kakb - 1D0*k2k3*m1**2 + 2D0*kbk3*m1**2 - 1D0*kak3*m4**2 + 2D0*kak3*pap2 - 2D0*kbk3*pap2 + 2D0*kak3*pbp2 + 
     &                   (k2k3*m1**2 + kak3*m4**2 - 2D0*kak3*pap2)*(m4**2 - 2D0*pbp2)*xv)) + 
     &             D3Lc*m1*mq3c*(m3**2*(-1D0*k1k2 - 1D0*k2k3 + 2D0*(kbk3 + pbp1) + (k1k2 + k2k3)*(m4**2 - 2D0*pbp2)*xv) + 
     &                k1k3*(k2k3 - 2D0*kbk3 - 1D0*k2k3*m4**2*xv + 2D0*k2k3*pbp2*xv)))) + 
     &       A2L*(A3Rc*m3*(k1k3 + m3**2)*(D3Rc*m1*(2D0*kakb + m4**2 - 1D0*pap2 - 2D0*pbp2 - 1D0*(m4**2 - 1D0*pap2)*(m4**2 - 2D0*pbp2)*xv) + 
     &             D3Lc*mq3c*(2D0*kakb - 1D0*pap2 + pap2*(m4**2 - 2D0*pbp2)*xv)) + 
     &          A3Lc*(D3Lc*(m3**2*(-2D0*k1k2*kakb - 2D0*k2k3*kakb + 4D0*kak3*kakb + k1k2*m1**2 + k2k3*m1**2 - 2D0*kbk3*m1**2 + kak3*m4**2 + 4D0*kakb*pap1 + m4**2*pap1 - 
     &                   2D0*kak3*pap2 + 2D0*kbk3*pap2 - 2D0*pap1*pap2 - 2D0*m1**2*pbp1 + 2D0*pap2*pbp1 - 2D0*kak3*pbp2 - 2D0*pap1*pbp2 - 
     &                   1D0*((k1k2 + k2k3)*m1**2 + (kak3 + pap1)*(m4**2 - 2D0*pap2))*(m4**2 - 2D0*pbp2)*xv) + 
     &                k1k3*(2D0*k2k3*kakb - 4D0*kak3*kakb - 1D0*k2k3*m1**2 + 2D0*kbk3*m1**2 - 1D0*kak3*m4**2 + 2D0*kak3*pap2 - 2D0*kbk3*pap2 + 2D0*kak3*pbp2 + 
     &                   (k2k3*m1**2 + kak3*m4**2 - 2D0*kak3*pap2)*(m4**2 - 2D0*pbp2)*xv)) + 
     &             D3Rc*m1*mq3c*(m3**2*(-1D0*k1k2 - 1D0*k2k3 + 2D0*(kbk3 + pbp1) + (k1k2 + k2k3)*(m4**2 - 2D0*pbp2)*xv) + 
     &                k1k3*(k2k3 - 2D0*kbk3 - 1D0*k2k3*m4**2*xv + 2D0*k2k3*pbp2*xv))))))/(k1k3**2*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)*(m2**2 + m4**2 - 1D0*mq2**2 - 2D0*pbp2 - (cI*mq2*wst(ptmeigen+sumT)))))


        radvTvUi = dreal((CLL2*GS*(A2L*(A3Rc*m3*(D3Rc*m1*(k2k3 - 1D0*kak3 + 2D0*kakb - 2D0*pbp2 - 1D0*(m4**2 - 1D0*pap2)*(k2k3 - 2D0*pbp2)*xv) + 
     &             D3Lc*mq3c*(-1D0*kak3 + 2D0*kakb + pap2*(k2k3 - 2D0*pbp2)*xv)) + 
     &          A3Lc*(D3Rc*m1*mq3c*(-1D0*k1k3 + 2D0*pbp1 + k1k2*(k2k3 - 2D0*pbp2)*xv) + 
     &             D3Lc*(k1k2*kak3 - 2D0*k1k2*kakb + k2k3*pap1 - 2D0*kak3*pap1 + 4D0*kakb*pap1 + (m1**2 - 1D0*pap2)*(k1k3 - 2D0*pbp1) - 2D0*pap1*pbp2 - 
     &                1D0*(k1k2*m1**2 + m4**2*pap1 - 2D0*pap1*pap2)*(k2k3 - 2D0*pbp2)*xv))) + 
     &       A2R*(A3Lc*m3*(D3Lc*m1*(k2k3 - 1D0*kak3 + 2D0*kakb - 2D0*pbp2 - 1D0*(m4**2 - 1D0*pap2)*(k2k3 - 2D0*pbp2)*xv) + 
     &             D3Rc*mq3c*(-1D0*kak3 + 2D0*kakb + pap2*(k2k3 - 2D0*pbp2)*xv)) + 
     &          A3Rc*(D3Lc*m1*mq3c*(-1D0*k1k3 + 2D0*pbp1 + k1k2*(k2k3 - 2D0*pbp2)*xv) + 
     &             D3Rc*(k1k2*kak3 - 2D0*k1k2*kakb + k2k3*pap1 - 2D0*kak3*pap1 + 4D0*kakb*pap1 + (m1**2 - 1D0*pap2)*(k1k3 - 2D0*pbp1) - 2D0*pap1*pbp2 - 
     &                1D0*(k1k2*m1**2 + m4**2*pap1 - 2D0*pap1*pap2)*(k2k3 - 2D0*pbp2)*xv)))))/
     &   (kbk3*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)*(-2D0*k2k3 + 2D0*kbk3 - 1D0*m2**2 - 1D0*m4**2 + mq2**2 + 2D0*pbp2 + (cI*mq2*wst(ptmeigen+sumT)))))


        radvTvUf = dreal((CLL2*GS*(A2L*(2D0*A3Lc*(D3Rc*m1*(-1D0*k1k3 + m3**2)*mq3c + D3Lc*(2D0*pap1*(-1D0*k1k2 - 1D0*k2k3 + kak3 + pap1) + m3**2*(-1D0*m1**2 + pap2))) + 
     &          A3Lc*(D3Rc*m1*(-2D0*k1k2*(k1k2 + k2k3) + k1k3*m4**2)*mq3c + 
     &             D3Lc*(2D0*k1k2**2*m1**2 + k1k2*(2D0*k2k3*m1**2 + (kak3 + 2D0*pap1)*(m4**2 - 2D0*pap2)) + k2k3*pap1*(m4**2 - 2D0*pap2) - 
     &                1D0*k1k3*(m1**2*m4**2 + (m4**2 - 2D0*pap2)*pap2)))*xv + 
     &          A3Rc*m3*(D3Rc*m1*(-2D0*k1k2 - 4D0*k2k3 + 4D0*kak3 + 2D0*pap1 + ((2D0*k1k2 + k2k3 - 1D0*kak3)*m4**2 - 2D0*k1k2*pap2)*xv) + 
     &             D3Lc*mq3c*(2D0*pap1 - 2D0*(k1k2 + k2k3)*pap2*xv + kak3*(-2D0 + m4**2*xv)))) + 
     &       A2R*(2D0*A3Rc*(D3Lc*m1*(-1D0*k1k3 + m3**2)*mq3c + D3Rc*(2D0*pap1*(-1D0*k1k2 - 1D0*k2k3 + kak3 + pap1) + m3**2*(-1D0*m1**2 + pap2))) + 
     &          A3Rc*(D3Lc*m1*(-2D0*k1k2*(k1k2 + k2k3) + k1k3*m4**2)*mq3c + 
     &             D3Rc*(2D0*k1k2**2*m1**2 + k1k2*(2D0*k2k3*m1**2 + (kak3 + 2D0*pap1)*(m4**2 - 2D0*pap2)) + k2k3*pap1*(m4**2 - 2D0*pap2) - 
     &                1D0*k1k3*(m1**2*m4**2 + (m4**2 - 2D0*pap2)*pap2)))*xv + 
     &          A3Lc*m3*(D3Lc*m1*(-2D0*k1k2 - 4D0*k2k3 + 4D0*kak3 + 2D0*pap1 + ((2D0*k1k2 + k2k3 - 1D0*kak3)*m4**2 - 2D0*k1k2*pap2)*xv) + 
     &             D3Rc*mq3c*(2D0*pap1 - 2D0*(k1k2 + k2k3)*pap2*xv + kak3*(-2D0 + m4**2*xv))))))/
     &   (k1k3*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)*(2D0*k2k3 - 2D0*kbk3 + m2**2 + m4**2 - 1D0*mq2**2 - 2D0*pbp2 - (cI*mq2*wst(ptmeigen+sumT)))))


        NeuQ2qx_M2RadiationVTU = NeuQ2qx_M2RadiationVTU + 2.*CF*NC*(radvTiUi + radvTiUf + radvTpUi + radvTpUf + radvTfUi + radvTfUf + radvTvUi + radvTvUf)

        enddo
        enddo
    
    
        end

! =================================================================
        double precision function NeuQ2qx_M2RadiationVUU()
        implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"

        double precision radvUiUi, radvUiUf, radvUfUf
        double complex A3L, A3Lc, A3R, A3Rc, D3L, D3Lc, D3R, D3Rc
        double precision mq3, mq3c
        double precision xv, colorfactor
        integer sumU, sumUc, sumUend

         if (xtype.eq.1) then
        
           xv = 0d0
           colorfactor = 4d0/3d0
           sumUend = -1
           
        elseif (xtype.eq.2) then
        
           xv = 0d0
           colorfactor = 1d0
           sumUend = -1
           
        elseif (xtype.eq.3) then
        
           xv = 1/(m4**2)
           colorfactor = 1d0
           sumUend = 3
    
           
        elseif (xtype.eq.4) then

           xv = 1/(m4**2)
           colorfactor = 1d0
           sumUend = 1
                   
        endif

        NeuQ2qx_M2RadiationVUU = 0d0

        do sumU = 0, sumUend
        do sumUc = 0, sumUend
    
        radvUiUi = 0d0
        radvUiUf = 0d0
        radvUfUf = 0d0

        A3L = AL(pumeigen+sumU,sfmeigen,sftype,sfgen)
        A3R = AR(pumeigen+sumU,sfmeigen,sftype,sfgen)
        A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
        A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))

        D3L = DL(chimeigen,pumeigen+sumU,sftype,xtype)
        D3R = DR(chimeigen,pumeigen+sumU,sftype,xtype)
        D3Lc = dconjg(DL(chimeigen,pumeigen+sumUc,sftype,xtype))
        D3Rc = dconjg(DR(chimeigen,pumeigen+sumUc,sftype,xtype))

        mq3 = MNeuCha(pumeigen+sumU)
        mq3c = MNeuCha(pumeigen+sumUc)

        radvUiUi = dreal((2D0*GS**2*(kbk3 - 1D0*m2**2)*(A3L*(A3Rc*m3*(-1D0*D3R*mq3*
     &              (D3Lc*m1*mq3c*(4D0 - 1D0*m4**2*xv) + D3Rc*(m1**2*(-2D0 + m4**2*xv) + pap2*(2D0 + m4**2*xv - 2D0*pap2*xv))) + 
     &             D3L*(D3Rc*m1*(m1**2 + m4**2 - 2D0*pap2)*(-4D0 + m4**2*xv) - 1D0*D3Lc*mq3c*(m1**2*(-2D0 + m4**2*xv) + pap2*(2D0 + m4**2*xv - 2D0*pap2*xv)))) + 
     &          A3Lc*(D3R*mq3*(-1D0*D3Lc*m1*(k1k2 - 1D0*pap1)*(-4D0 + m4**2*xv) + D3Rc*mq3c*(2D0*k1k2*pap2*xv + pap1*(2D0 - 1D0*m4**2*xv))) + 
     &             D3L*(-1D0*D3Rc*m1*mq3c*(k1k2 - 1D0*pap1)*(-4D0 + m4**2*xv) + 
     &                D3Lc*(pap1*(-2D0*m4**2 + (m4**2 - 2D0*pap2)**2*xv + m1**2*(2D0 - 1D0*m4**2*xv)) + k1k2*(4D0*pap2 + 2D0*m1**2*(-2D0 + m4**2*xv - 1D0*pap2*xv)))))) + 
     &       A3R*(A3Lc*m3*(-1D0*D3Rc*mq3c*(D3L*m1*mq3*(4D0 - 1D0*m4**2*xv) + D3R*(m1**2*(-2D0 + m4**2*xv) + pap2*(2D0 + m4**2*xv - 2D0*pap2*xv))) + 
     &             D3Lc*(D3R*m1*(m1**2 + m4**2 - 2D0*pap2)*(-4D0 + m4**2*xv) - 1D0*D3L*mq3*(m1**2*(-2D0 + m4**2*xv) + pap2*(2D0 + m4**2*xv - 2D0*pap2*xv)))) + 
     &          A3Rc*(D3L*mq3*(-1D0*D3Rc*m1*(k1k2 - 1D0*pap1)*(-4D0 + m4**2*xv) + D3Lc*mq3c*(2D0*k1k2*pap2*xv + pap1*(2D0 - 1D0*m4**2*xv))) + 
     &             D3R*(-1D0*D3Lc*m1*mq3c*(k1k2 - 1D0*pap1)*(-4D0 + m4**2*xv) + 
     &                D3Rc*(pap1*(-2D0*m4**2 + (m4**2 - 2D0*pap2)**2*xv + m1**2*(2D0 - 1D0*m4**2*xv)) + k1k2*(4D0*pap2 + 2D0*m1**2*(-2D0 + m4**2*xv - 1D0*pap2*xv))))))))/
     &   (kbk3**2*(m1**2 + m4**2 - 1D0*mq3**2 - 2D0*pap2)*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)))


        radvUiUf = dreal((-1D0*GS**2*(A3L*(A3Rc*m3*(D3L*(D3Rc*m1*(m1**2 + m4**2 - 2D0*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1)*(-4D0 + m4**2*xv) + 
     &                D3Lc*mq3c*(-2D0*k2k3*kakb + 2D0*(m1**2 - 1D0*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + 2D0*kak3*pbp2 + 
     &                   (k2k3*kakb*(m4**2 - 2D0*pap2) - 1D0*(m1**2*m4**2 + (m4**2 - 2D0*pap2)*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*kak3*(m4**2 - 2D0*pap2)*pbp2)*xv)) + 
     &             D3R*mq3*(D3Lc*m1*mq3c*(k1k3 - 1D0*kbk3 - 2D0*pbp1)*(-4D0 + m4**2*xv) + 
     &                D3Rc*(2D0*(m1**2 - 1D0*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 2D0*kak3*pbp2 + 
     &                   (-1D0*(m1**2*m4**2 + (m4**2 - 2D0*pap2)*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + kak3*(m4**2 - 2D0*pap2)*pbp2)*xv + 
     &                   k2k3*kakb*(2D0 - 1D0*m4**2*xv + 2D0*pap2*xv)))) + 
     &          A3Lc*(D3R*mq3*(D3Lc*m1*((k2k3 - 1D0*kak3)*pbp1 - 1D0*pap1*(kbk3 + 2D0*pbp1) + k1k2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) + k1k3*(kakb + pap1 - 1D0*pbp2))*
     &                 (-4D0 + m4**2*xv) + D3Rc*mq3c*(2D0*k1k3*(kakb + pap1) - 2D0*(kak3 + 2D0*pap1)*pbp1 + (m4**2*(kak3 + 2D0*pap1) - 2D0*(2D0*k1k2 + k2k3)*pap2)*pbp1*xv - 
     &                   1D0*k1k3*(m4**2*(kakb + pap1) - 2D0*pap2*(k1k2 + pbp2))*xv + kbk3*(-2D0*pap1 + m4**2*pap1*xv - 2D0*k1k2*pap2*xv))) + 
     &             D3L*(D3Rc*m1*mq3c*((k2k3 - 1D0*kak3)*pbp1 - 1D0*pap1*(kbk3 + 2D0*pbp1) + k1k2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) + k1k3*(kakb + pap1 - 1D0*pbp2))*(-4D0 + m4**2*xv) + 
     &                D3Lc*(2D0*(kbk3*(-1D0*m1**2 + m4**2)*pap1 + k1k3*(m1 - 1D0*m4)*(m1 + m4)*(kakb + pap1) + 
     &                      (-1D0*(m1 - 1D0*m4)*(m1 + m4)*(kak3 + 2D0*pap1) + 2D0*k2k3*(m1**2 - 1D0*pap2))*pbp1 + 2D0*k1k3*(-1D0*m1**2 + pap2)*pbp2) + 
     &                   (k1k3*(kakb + pap1)*(-1D0*m1*m4 + m4**2 - 2D0*pap2)*(m4*(m1 + m4) - 2D0*pap2) + kbk3*pap1*(m4*(m1 + m4) - 2D0*pap2)*((m1 - 1D0*m4)*m4 + 2D0*pap2) + 
     &                      (2D0*k2k3*m1**2*(-1D0*m4**2 + pap2) + (kak3 + 2D0*pap1)*(m4*(m1 + m4) - 2D0*pap2)*((m1 - 1D0*m4)*m4 + 2D0*pap2))*pbp1 + 
     &                      2D0*k1k3*m1**2*(m4**2 - 1D0*pap2)*pbp2)*xv - 2D0*k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1)*(-2D0*pap2 + m1**2*(2D0 - 1D0*m4**2*xv + pap2*xv)))))) + 
     &       A3R*(A3Lc*m3*(D3Rc*mq3c*(2D0*D3R*(-1D0*k2k3*kakb + (m1**2 - 1D0*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + kak3*pbp2) + 
     &                D3R*(k2k3*kakb*(m4**2 - 2D0*pap2) - 1D0*(m1**2*m4**2 + (m4**2 - 2D0*pap2)*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 1D0*kak3*(m4**2 - 2D0*pap2)*pbp2)*xv + 
     &                D3L*m1*mq3*(k1k3 - 1D0*kbk3 - 2D0*pbp1)*(-4D0 + m4**2*xv)) + 
     &             D3Lc*(D3R*m1*(m1**2 + m4**2 - 2D0*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1)*(-4D0 + m4**2*xv) + 
     &                D3L*mq3*(2D0*(m1**2 - 1D0*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) - 2D0*kak3*pbp2 + 
     &                   (-1D0*(m1**2*m4**2 + (m4**2 - 2D0*pap2)*pap2)*(k1k3 - 1D0*kbk3 - 2D0*pbp1) + kak3*(m4**2 - 2D0*pap2)*pbp2)*xv + 
     &                   k2k3*kakb*(2D0 - 1D0*m4**2*xv + 2D0*pap2*xv)))) + 
     &          A3Rc*(D3L*mq3*(D3Rc*m1*((k2k3 - 1D0*kak3)*pbp1 - 1D0*pap1*(kbk3 + 2D0*pbp1) + k1k2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) + k1k3*(kakb + pap1 - 1D0*pbp2))*
     &                 (-4D0 + m4**2*xv) + D3Lc*mq3c*(2D0*k1k3*(kakb + pap1) - 2D0*(kak3 + 2D0*pap1)*pbp1 + (m4**2*(kak3 + 2D0*pap1) - 2D0*(2D0*k1k2 + k2k3)*pap2)*pbp1*xv - 
     &                   1D0*k1k3*(m4**2*(kakb + pap1) - 2D0*pap2*(k1k2 + pbp2))*xv + kbk3*(-2D0*pap1 + m4**2*pap1*xv - 2D0*k1k2*pap2*xv))) + 
     &             D3R*(D3Lc*m1*mq3c*((k2k3 - 1D0*kak3)*pbp1 - 1D0*pap1*(kbk3 + 2D0*pbp1) + k1k2*(-1D0*k1k3 + kbk3 + 2D0*pbp1) + k1k3*(kakb + pap1 - 1D0*pbp2))*(-4D0 + m4**2*xv) + 
     &                D3Rc*(2D0*(kbk3*(-1D0*m1**2 + m4**2)*pap1 + k1k3*(m1 - 1D0*m4)*(m1 + m4)*(kakb + pap1) + 
     &                      (-1D0*(m1 - 1D0*m4)*(m1 + m4)*(kak3 + 2D0*pap1) + 2D0*k2k3*(m1**2 - 1D0*pap2))*pbp1 + 2D0*k1k3*(-1D0*m1**2 + pap2)*pbp2) + 
     &                   (k1k3*(kakb + pap1)*(-1D0*m1*m4 + m4**2 - 2D0*pap2)*(m4*(m1 + m4) - 2D0*pap2) + kbk3*pap1*(m4*(m1 + m4) - 2D0*pap2)*((m1 - 1D0*m4)*m4 + 2D0*pap2) + 
     &                      (2D0*k2k3*m1**2*(-1D0*m4**2 + pap2) + (kak3 + 2D0*pap1)*(m4*(m1 + m4) - 2D0*pap2)*((m1 - 1D0*m4)*m4 + 2D0*pap2))*pbp1 + 
     &                      2D0*k1k3*m1**2*(m4**2 - 1D0*pap2)*pbp2)*xv - 2D0*k1k2*(k1k3 - 1D0*kbk3 - 2D0*pbp1)*(-2D0*pap2 + m1**2*(2D0 - 1D0*m4**2*xv + pap2*xv))))))))/
     &   (k1k3*kbk3*(m1**2 + m4**2 - 1D0*mq3**2 - 2D0*pap2)*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)))


        radvUfUf = dreal((2D0*GS**2*(A3L*(-1D0*A3Rc*m3*(k1k3 + m3**2)*(-1D0*D3R*mq3*
     &              (D3Lc*m1*mq3c*(4D0 - 1D0*m4**2*xv) + D3Rc*(m1**2*(-2D0 + m4**2*xv) + pap2*(2D0 + m4**2*xv - 2D0*pap2*xv))) + 
     &             D3L*(D3Rc*m1*(m1**2 + m4**2 - 2D0*pap2)*(-4D0 + m4**2*xv) - 1D0*D3Lc*mq3c*(m1**2*(-2D0 + m4**2*xv) + pap2*(2D0 + m4**2*xv - 2D0*pap2*xv)))) + 
     &          A3Lc*(D3R*mq3*(-1D0*D3Lc*m1*(k1k3*(k2k3 - 1D0*kak3) + m3**2*(-1D0*k1k2 - 1D0*k2k3 + kak3 + pap1))*(-4D0 + m4**2*xv) + 
     &                D3Rc*mq3c*(k1k3*(2D0*kak3 - 1D0*kak3*m4**2*xv + 2D0*k2k3*pap2*xv) + m3**2*(-2D0*(kak3 + pap1) + (m4**2*(kak3 + pap1) - 2D0*(k1k2 + k2k3)*pap2)*xv))) + 
     &             D3L*(-1D0*D3Rc*m1*mq3c*(k1k3*(k2k3 - 1D0*kak3) + m3**2*(-1D0*k1k2 - 1D0*k2k3 + kak3 + pap1))*(-4D0 + m4**2*xv) + 
     &                D3Lc*(k1k3*(kak3*(-2D0*m4**2 + (m4**2 - 2D0*pap2)**2*xv + m1**2*(2D0 - 1D0*m4**2*xv)) + k2k3*(4D0*pap2 + 2D0*m1**2*(-2D0 + m4**2*xv - 1D0*pap2*xv))) + 
     &                   m3**2*((kak3 + pap1)*(2D0*m4**2 - 1D0*(m4**2 - 2D0*pap2)**2*xv + m1**2*(-2D0 + m4**2*xv)) + 2D0*k1k2*(-2D0*pap2 + m1**2*(2D0 - 1D0*m4**2*xv + pap2*xv)) + 
     &                      2D0*k2k3*(-2D0*pap2 + m1**2*(2D0 - 1D0*m4**2*xv + pap2*xv))))))) + 
     &       A3R*(-1D0*A3Lc*m3*(k1k3 + m3**2)*(-1D0*D3Rc*mq3c*(D3L*m1*mq3*(4D0 - 1D0*m4**2*xv) + D3R*(m1**2*(-2D0 + m4**2*xv) + pap2*(2D0 + m4**2*xv - 2D0*pap2*xv))) + 
     &             D3Lc*(D3R*m1*(m1**2 + m4**2 - 2D0*pap2)*(-4D0 + m4**2*xv) - 1D0*D3L*mq3*(m1**2*(-2D0 + m4**2*xv) + pap2*(2D0 + m4**2*xv - 2D0*pap2*xv)))) + 
     &          A3Rc*(D3L*mq3*(-1D0*D3Rc*m1*(k1k3*(k2k3 - 1D0*kak3) + m3**2*(-1D0*k1k2 - 1D0*k2k3 + kak3 + pap1))*(-4D0 + m4**2*xv) + 
     &                D3Lc*mq3c*(k1k3*(2D0*kak3 - 1D0*kak3*m4**2*xv + 2D0*k2k3*pap2*xv) + m3**2*(-2D0*(kak3 + pap1) + (m4**2*(kak3 + pap1) - 2D0*(k1k2 + k2k3)*pap2)*xv))) + 
     &             D3R*(-1D0*D3Lc*m1*mq3c*(k1k3*(k2k3 - 1D0*kak3) + m3**2*(-1D0*k1k2 - 1D0*k2k3 + kak3 + pap1))*(-4D0 + m4**2*xv) + 
     &                D3Rc*(k1k3*(kak3*(-2D0*m4**2 + (m4**2 - 2D0*pap2)**2*xv + m1**2*(2D0 - 1D0*m4**2*xv)) + k2k3*(4D0*pap2 + 2D0*m1**2*(-2D0 + m4**2*xv - 1D0*pap2*xv))) + 
     &                   m3**2*((kak3 + pap1)*(2D0*m4**2 - 1D0*(m4**2 - 2D0*pap2)**2*xv + m1**2*(-2D0 + m4**2*xv)) + 2D0*k1k2*(-2D0*pap2 + m1**2*(2D0 - 1D0*m4**2*xv + pap2*xv)) + 
     &                      2D0*k2k3*(-2D0*pap2 + m1**2*(2D0 - 1D0*m4**2*xv + pap2*xv)))))))))/
     &   (k1k3**2*(m1**2 + m4**2 - 1D0*mq3**2 - 2D0*pap2)*(m1**2 + m4**2 - 1D0*mq3c**2 - 2D0*pap2)))
                  
    
        NeuQ2qx_M2RadiationVUU = NeuQ2qx_M2RadiationVUU + CF*NC*(radvUiUi + 2.*radvUiUf + radvUfUf)

        enddo
        enddo
    
       
        end
