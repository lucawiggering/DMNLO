*
* =====================================================================
*
*  Function returns Alpha_s/(4*pi) 
*  given in DRbar with 6 active quark flavours in the MSSM
*  
*  Using three-loop running for MSbar,SM with two-loop decoupling and 
*  MS-DR conversion
*
*  no running in MSSM is implemented -> decoupling happens at the final scale
*
*------------------------------------------------------------------------
*  All formulas for decoupling coefficients are taken from 

*  Harlander, Mihaila, Steinhauser - arXiv: 0706.2953 [hep-ph]
*  'Running of alpha_s and m_b in the MSSM'

*  Harlander, Mihaila, Steinhauser - hep-ph/0509048
*  'Two-loop matching coefficients for the strong coupling in the MSSM'

*  Bauer, Mihaila, Salomon         - arXiv: 0810.5101 [hep-ph]
*  'Matching coefficients for alpha_s and m_b to O(a^2) in the MSSM'

 
      FUNCTION ASDR (Qinp2)

#include "../util/DMNLO_Model.h"

      double precision ASDR,Qinp2
      double precision ASMS5Q,ASDR5Q
      double precision zetaDC(2),test,test1,test2
      
* external functions      
      double precision AS,mssmdec2l,mssmdec1l,mssmzeta2l
      
      external ASzetaDecCoe

* three-loop running of alpha_s^{MS,SM,nf=5} returns alpha/pi	
      ASMS5Q = 4d0*AS( Qinp2, MZ**2, AlfasMZ/(4d0*pi), 5, 3, 5)

* two-loop MS-DR conversion of alpha_s^{MS,SM,nf=5}/pi
* Eq.7 in hep-ph/0509048 (no active gluino, nf=5)
      ASDR5Q = ASMS5Q*(1d0 + ASMS5Q/4d0 + ASMS5Q**2*(11d0/8d0 - 5d0/12d0))

* two-loop SUSY decoupling of alpha_s^{DR,SM,nf=5}/pi
* adding top-quark & all relevant SUSY particles in one step

* decoupling coefficients
      call ASzetaDecCoe(Qinp2,zetaDC)

* running of alpha_s^{DR,MSSM,nf=6} returns alpha/(4*pi)	
      ASDR = 1d0/4d0*ASDR5Q*(1d0 - ASDR5Q*zetaDC(1) + ASDR5Q**2*(2d0*zetaDC(1)**2 - zetaDC(2)))

* in case of problematic stop,sbottom,gluino masses zetaDC(2) breaks down and causes ASDR<0
      if(ASDR.lt.0d0)then
      
        print*,'ASDR<0 - reducing decoupling to one-loop'
        ASDR = 1d0/4d0*ASDR5Q*(1d0 - ASDR5Q*zetaDC(1))
        
      endif
            
      RETURN
      END

* =====================================================================
*  Subroutine for SUSY two-loop decoupling coefficients for alpha_s
*
      SUBROUTINE ASzetaDecCoe(Qdec2,zetaDC)
      
#include "../util/DMNLO_Model.h"

      double precision zetaDC(2),Qdec2
      double precision logsquark,Mmass,Mtop,Mtop2
      double precision logT,logM,logTM

      double precision mssmzeta2l

* Common SUSY mass of colored particles (for two-loop matching)
      Mmass = (MSf(1,4,3)*MSf(1,3,3)*MSf(1,4,2)*MSf(1,3,2)*MGl)**(1d0/5d0)
      Mtop = MTos
      Mtop2 = MTos**2
      
* short forms for logarithms
      logT = dlog(Qdec2/Mtop2)
      logM = dlog(Qdec2/Mmass**2)
      logTM = dlog(Mtop2/Mmass**2)

* based on Eq.20,21 in hep-ph/0509048 and Eq.(4.1) in arXiv: 0810.5101      
      logsquark = dlog(Qdec2/MSf2(1,4,3)) + dlog(Qdec2/MSf2(2,4,3))
     &          + dlog(Qdec2/MSf2(1,3,3)) + dlog(Qdec2/MSf2(2,3,3))
     &          + dlog(Qdec2/MSf2(1,4,2)) + dlog(Qdec2/MSf2(2,4,2))
     &          + dlog(Qdec2/MSf2(1,3,2)) + dlog(Qdec2/MSf2(2,3,2))
     &          + dlog(Qdec2/MSf2(1,4,1)) + dlog(Qdec2/MSf2(2,4,1))
     &          + dlog(Qdec2/MSf2(1,3,1)) + dlog(Qdec2/MSf2(2,3,1))
      
* one-loop matching with full mass dependance 
      zetaDC(1) = -1d0/6d0*dlog(Qdec2/Mtop2) - 1d0/2d0*dlog(Qdec2/MGl**2) - 1d0/24d0*logsquark
      
* two-loop matching with one SUSY mass only (Eq.6 of arXiv: 0706.2953)
c      zetaDC(2) = -215d0/96d0 - 19d0/24d0*logT - 5d0/2d0*logM + (1d0/6d0*logT + logM)**2
c     &        + (Mtop2/Mmass**2)*(5d0/48d0 + 3d0/8d0*logTM) - 7d0*pi/36d0*(Mtop/Mmass)**3     
c     &        + (Mtop2/Mmass**2)**2*(881d0/7200d0 - 1d0/80d0*logTM) + 7d0*pi/288d0*(Mtop/Mmass)**5     

* two-loop matching - exact from L.Mihaila
      zetaDC(2) = -mssmzeta2l(dsqrt(Qdec2),MGl,MTOS,MSf(1,3,3),MSf(2,3,3),Xf(3,3),MSf(1,4,3),MSf(2,4,3),MSf(1,3,2),MSf(2,3,2))
      
      END


* ..The running coupling of QCD,  
*
*         AS  =  a_s  =  alpha_s(mu_r^2)/(4 pi),
*
*    obtained by integrating the evolution equation for a fixed number
*    of massless flavours  NF.  Except at leading order (LO),  AS  is 
*    obtained using a fourth-order Runge-Kutta integration. 
*
* ..The initial and final scales  R20  and  R2,  the value  AS0  at
*    R20, and  NF  are passed as function arguments.  The coefficients 
*    of the beta function up to  a_s^5 (N^3LO)  are provided by the 
*    common-block  BETA.  The order of the expansion  NAORD  (defined 
*    as the 'n' in N^nLO) and the number of steps  NASTPS  for the 
*    integration beyond LO are given by the common-block  ASPAR.
*
* =====================================================================
*
*
       FUNCTION AS (R2, R20, AS0, NF, NAORD, NASTPS)
*
       IMPLICIT DOUBLE PRECISION (A - Z)
       INTEGER NFMIN, NFMAX, NF, NAORD, NASTPS, K1
c       double precision FBETA3, FBETA1
       PARAMETER (NFMIN = 3, NFMAX = 6)
       PARAMETER ( SXTH = 0.16666 66666 66666 D0 )

	external BETAFCT
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks 
*
       COMMON / BETA   / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     ,                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
*
* ..The beta functions FBETAn at N^nLO for n = 1, 2, and 3
*
       FBETA2(A) = - A**2 * ( BETA0(NF) + A * ( BETA1(NF)
     ,                        + A * BETA2(NF) ) )
       FBETA3(A) = - A**2 * ( BETA0(NF) + A * ( BETA1(NF)
     ,                        + A * (BETA2(NF) + A * BETA3(NF)) ) )
*
* ---------------------------------------------------------------------
*
* ..Initial value, evolution distance and step size
*

       call BETAFCT

       AS = AS0
       LRRAT = LOG (R2/R20)
       DLR = LRRAT / NASTPS
*
* ..Solution of the evolution equation depending on  NAORD
*   (fourth-order Runge-Kutta beyond the leading order)
*
       IF (NAORD .EQ. 0) THEN
*
         AS = AS0 / (1.+ BETA0(NF) * AS0 * LRRAT)
*
       ELSE IF (NAORD .EQ. 1) THEN
*
       DO 2 K1 = 1, NASTPS
         XK0 = DLR * FBETA1 (AS,NF)
         XK1 = DLR * FBETA1 (AS + 0.5 * XK0,NF)
         XK2 = DLR * FBETA1 (AS + 0.5 * XK1,NF)
         XK3 = DLR * FBETA1 (AS + XK2,NF)
         AS = AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  2    CONTINUE
*
       ELSE IF (NAORD .EQ. 2) THEN
*
       DO 3 K1 = 1, NASTPS
         XK0 = DLR * FBETA2 (AS)
         XK1 = DLR * FBETA2 (AS + 0.5 * XK0)
         XK2 = DLR * FBETA2 (AS + 0.5 * XK1)
         XK3 = DLR * FBETA2 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  3    CONTINUE
*  
       ELSE IF (NAORD .EQ. 3) THEN
*
       DO 4 K1 = 1, NASTPS
         XK0 = DLR * FBETA3 (AS)
         XK1 = DLR * FBETA3 (AS + 0.5 * XK0)
         XK2 = DLR * FBETA3 (AS + 0.5 * XK1)
         XK3 = DLR * FBETA3 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  4    CONTINUE
       END IF
*
* ---------------------------------------------------------------------
*
       RETURN
       END
*
* =================================================================av==
*
* ..File: betafct.f
*
*
* ..The subroutine BETAFCT for the coefficients  BETA0...BETA3  of the 
*    beta function of QCD up to order alpha_s^5 (N^3LO), normalized by 
*
*        d a_s / d ln mu_r^2  =  - BETA0 a_s^2 - BETA1 a_s^3 - ... 
*
*    with  a_s = alpha_s/(4*pi). 
*
* ..The MSbar coefficients are written to the common-block  BETA  for 
*   NF = 3...6  (parameters NFMIN, NFMAX) quark flavours.
*
* ..The factors CF, CA and TF  are taken from the common-block  COLOUR. 
*    Beyond NLO the QCD colour factors are hard-wired in this routine,
*    and the numerical coefficients are truncated to six digits.
*
* =====================================================================
*
*
       SUBROUTINE BETAFCT
*
       IMPLICIT DOUBLE PRECISION (A - Z)
       INTEGER NFMIN, NFMAX, NF
       PARAMETER (NFMIN = 3, NFMAX = 6)
*
* ---------------------------------------------------------------------
*
* ..Input common-block
*
**       COMMON / COLOUR / CF, CA, TR
*
* ..Output common-block
*
       COMMON / BETA   / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     1                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)

	CA = 3d0
	CF = 4/3d0
	TR = 1/2d0

*
* ---------------------------------------------------------------------
*
* ..The full LO and NLO coefficients 
*
       B00 =  11./3.D0 * CA
       B01 =  -4./3.D0 * TR
       B10 =  34./3.D0 * CA**2
       B11 = -20./3.D0 * CA*TR - 4.* CF*TR
*
* ..Flavour-number loop and output to the array
*
       DO 1 NF = NFMIN, NFMAX
*
       BETA0(NF) = B00 + B01 * NF
       BETA1(NF) = B10 + B11 * NF
*
       BETA2(NF) = 1428.50 - 279.611 * NF + 6.01852 * NF**2
       BETA3(NF) = 29243.0 - 6946.30 * NF + 405.089 * NF**2 
     1             + 1.49931 * NF**3
*
* ---------------------------------------------------------------------
*
  1    CONTINUE
*
       RETURN
       END
*
* =================================================================av==

	FUNCTION FBETA1(A,NF)

        IMPLICIT DOUBLE PRECISION (B)
        INTEGER NFMIN, NFMAX, NF
        PARAMETER (NFMIN = 3, NFMAX = 6)
	
	double precision A,FBETA1

        COMMON / BETA   / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     1                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)

	FBETA1 = - A**2 * ( BETA0(NF) + A *   BETA1(NF) )

	end

* =================================================================av==

	FUNCTION FBETA3(A,NF)

        IMPLICIT DOUBLE PRECISION (B)
        INTEGER NFMIN, NFMAX, NF
        PARAMETER (NFMIN = 3, NFMAX = 6)
	
	double precision A,FBETA3

        COMMON / BETA   / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     1                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)

	FBETA3 = - A**2 *BETA0(NF) 

	FBETA3 = FBETA3 - A**3*BETA1(NF)
	FBETA3 = FBETA3 - A**4*BETA2(NF)
	FBETA3 = FBETA3 - A**5*BETA3(NF)

	FBETA3 = - A**2 * ( BETA0(NF) + A * ( BETA1(NF)
     ,                        + A * (BETA2(NF) + A * BETA3(NF)) ) )

	end


* ========================================================================
* == From here on below - auxillary functions for AS decoupling in MSSM ==
* ========================================================================

c==============================================================
c function b0fin

      real*8 function b0fin(q,m1,m2,mu)
c..
c..   This is the finite part of the B0 function.
c..   Only the real part is implemented; the imaginary part
c..   is set to zero.
c..   Note: not all cases are covered yet; corresponding 
c..   'if-endif' still missing.
c..   
c..   
c..   q2 == q^2
c..   m1 == m_1
c..   m2 == m_2
c..

      implicit real*8 (a-h,k-z)
      implicit  integer(i)
      
      complex*16 ieps

      ieps = dcmplx(0d0,1d-17)
     
      mu2 = mu*mu
      q2  = q*q

      if ((m2.eq.0d0).and.(m1.ne.0d0)) then
     	   b0fin = dreal(2d0 - dlog(m1**2/mu**2)
     &           - (q2-m1**2)/q2*cdlog((m1**2-q2-ieps)/m1**2))

      else if ((m1.eq.0d0).and.(m2.ne.0d0)) then
     	   b0fin = dreal(2d0 - dlog(m2**2/mu**2)
     &           - (q2-m2**2)/q2*cdlog((m2**2-q2-ieps)/m2**2))

      else

       if ( (q2.le.(m1-m2)**2) ) then
         b0fin = 
     .        (2.d0-dlog(m1*m2/mu2)
     .        +(m1**2-m2**2)/q2*dlog(m2/m1)
     .        + dsqrt((m1+m2)**2-q2)*dsqrt((m1-m2)**2-q2)/q2 *
     .        dlog(  
     .        (dsqrt((m1+m2)**2-q2)+dsqrt((m1-m2)**2-q2)) /
     .        (dsqrt((m1+m2)**2-q2)-dsqrt((m1-m2)**2-q2)) 
     .        ) )
       else if ( (q2.ge.(m1+m2)**2) ) then
         b0fin = 
     .        (2.d0-dlog(m1*m2/mu2)
     .        +(m1**2-m2**2)/q2*dlog(m2/m1)
     .        - dsqrt(q2-(m1+m2)**2)*dsqrt(q2-(m1-m2)**2)/q2 *
     .        dlog(-  
     .        (dsqrt(q2-(m1+m2)**2)+dsqrt(q2-(m1-m2)**2)) /
     .        (dsqrt(q2-(m1+m2)**2)-dsqrt(q2-(m1-m2)**2)) 
     .        ) )
       else if ( (q2.lt.(m1+m2)**2).and.(q2.gt.(m1-m2)**2) ) then
         b0fin = 
     .        (2.d0-dlog(m1*m2/mu2)
     .        +(m1**2-m2**2)/q2*dlog(m2/m1)
     .        - 2.d0/q2 *
     .        dsqrt( dabs(q2-(m1+m2)**2) )*
     .        dsqrt( dabs(q2-(m1-m2)**2) )*
     .        datan( dsqrt(dabs(q2-(m1-m2)**2)) / 
     .        dsqrt(dabs(q2-(m1+m2)**2)) )
     .        )
      else
         write(6,*) '<function b0fin>: ',
     .        'arguments (',q2,',',m1,',',m2,') not implemented.'
         stop
       endif
      endif
      return
      end


c =============================================================
c function Tfin

      real*8 function Tfin(M1,M2,M3)
      implicit real*8 (a-h,k-z)
      complex*16 dilog,dtphixyc
      real*8 M1,M2,M3

c Constants
      data pi/3.14159265358979323846264338328d0/
      data z2/1.6449340668482264d0/
      data z3/1.2020569031595942d0/

       xx=(M1/M3)**2
       yy=(M2/M3)**2
     
      Tfin=(-((M1**2 + M2**2 + M3**2)*(42.d0 + Pi**2)) +
     .     6.d0*M1**2*dlog(xx)**2 -
     .     6.d0*(M1**2 + M2**2 - M3**2)*dlog(xx)*dlog(yy) +
     .     6d0*M2**2*dlog(yy)**2)/12.d0
     .     -davytausk(xx,yy,M3)*0.5d0

       return
       end

c =============================================================
c function davytausk(xin,yin,m3in)

      real*8 function davytausk(xin,yin,m3in)
      implicit real*8 (a-h,k-z)
      complex*16 dilog,dtphixyc
      
c Constants
      data pi/3.14159265358979323846264338328d0/
      data z2/1.6449340668482264d0/
      data z3/1.2020569031595942d0/
   
     
  
      sqx = dsqrt(xin)
      sqy = dsqrt(yin)

      dtlam2 = (1.d0 - xin - yin)**2 - 4*xin*yin

c..   First the case lambda^2 >= 0, which is given by Eq.(4.10).
c..   But before, we need to do some mapping to cover all
c..   kinematical cases (see discussion below Eq.(4.10)).
      if (dtlam2.gt.0) then
c..   m1 + m2 <= m3:
         if     ( (sqx+sqy).le.1.d0 ) then
            xx = xin
            yy = yin
            m32 = m3in**2
c..   m2 + m3 <= m1  (m3 <-> m1):
         elseif ( (sqx-sqy).gt.1.d0 ) then
            xx = 1/xin
            yy = yin/xin
            m32 = xin*m3in**2
c..   m1 + m3 <= m2  (m2 <-> m3):
         elseif ( (sqy-sqx).gt.1.d0 ) then
            xx = xin/yin
            yy = 1/yin
            m32 = yin*m3in**2
         endif

c..   the sign of lambda^2 is invariant under the mappings above, but
c..   lambda itself is not:
         dtlam2xy = (1.d0 - xx - yy)**2 - 4*xx*yy
         dtlambda = dsqrt(dtlam2xy)

c..   Eq.(4.10) times lambda^2:
         dtphixyc = 
     &        dtlambda * m32 * ( 2*dlog( dtarg(xx,yy,dtlambda) )
     &        *dlog( dtarg(yy,xx,dtlambda) )
     &        - dlog(xx)*dlog(yy)
     &        - 2*dilog( dcmplx(dtarg(xx,yy,dtlambda) ))
     &        - 2*dilog( dcmplx(dtarg(yy,xx,dtlambda)) )
     &        + pi*pi/3.d0 )




c..   Now the case lambda^2 < 0:
      elseif (dtlam2.lt.0) then
c..   and again the mappings, but this time   
c..   m1 + m2 >= m3  is the default:
         if     ( (sqx+sqy).ge.1.d0 ) then
            xx = xin
            yy = yin
            m32 = m3in**2
         elseif ( (sqx-sqy).lt.1.d0 ) then
            xx = 1/xin
            yy = yin/xin
            m32 = xin*m3in**2
         elseif ( (sqy-sqx).lt.1.d0 ) then
            xx = xin/yin
            yy = 1/yin
            m32 = yin*m3in**2
         endif

         dtlam2xy = (1.d0 - xx - yy)**2 - 4*xx*yy
         dtlambda = dsqrt(-dtlam2xy)

c..   Eq.(4.15) times lambda^2:
         dtphixyc = - 2.d0 * m32 * dtlambda *
     &        ( dtclaus(1.d0,xx,yy)
     &        + dtclaus(yy,xx,1.d0)
     &        + dtclaus(xx,yy,1.d0) )

      else
         dtphixyc = 0.d0
      endif
         
      if (dabs(dimag(dtphixyc)).gt.1d-10) stop 1003
      davytausk = dreal(dtphixyc)

     

       return
       end

c ============================================================
      real*8 function dtclaus(arg1,arg2,arg3)
c..
c..   This is Clausen's function with the particular argument as 
c..   defined in Eq.(4.15) of
c..   Davydychev, Tausk,  NPB 397 (1993) 123.
c..
c..   dclaus(x) is taken from CERNLIB which has to be compiled in.
c..   
      implicit real*8 (a-z)
      external dclaus

      argu = 2*acos( (- arg1 + arg2 + arg3)/2.d0/dsqrt(arg2*arg3) )
      dtclaus = dclaus(argu)

      return
      end

c ============================================================


C  function dtarg:

      real*8 function dtarg(xin,yin,xylam)

      implicit real*8 (a-z)

      dtarg = (1 + xin - yin - xylam)/2.d0

      return
      end

c ============================================================
c   function Tfin0

      real*8 function Tfin0(M1,M2)
      implicit real*8 (a-h,k-z)
      complex*16 dtphix
      
   
      if ( M2.ge.M1 ) then  
       Tfin0 =   dreal(dtphix(M1,M2))
      else
       Tfin0 =   dreal(dtphix(M2,M1))
      endif

      return
      end
      
c ============================================================

c   function dtphix

      complex*16 function dtphix(m1,m2)
      complex*16 dilog,zz
      real*8 m1,m2

c Constants
      data pi/3.14159265358979323846264338328d0/
      data z2/1.6449340668482264d0/
      data z3/1.2020569031595942d0/
   
      
       zz=dcmplx(1.d0 - m1**2/m2**2)
     
       dtphix = (-7.d0*(m1**2 + m2**2))/2.d0
     .   - ((m1**2 + m2**2)*Pi**2)/12.d0 
     .   + (m1**2*(dlog(m1**2/m2**2))**2)/2.d0 
     .   +(m1**2 - m2**2)*dilog(zz)


      return
      end

c =============================================================


c =============================================================
c functions  to be linked as objects
c =============================================================

c function d123
      real*8 function d123(x,y,z)
      implicit  real*8(a-h,k-z)
      implicit  integer(i)


      d123= 1/((x - y - z)*(x + y - z)*(x - y + z)*(x + y + z))

      return 
      end


c =============================================================
c function denominator
      real*8 function den(x,y)
      implicit  real*8(a-h,k-z)
      implicit  integer(i)


      den=1/(x**2-y**2)
      return 
      end
c =============================================================
c logarithm as a function of two arguments

      real*8 function ln(x,y)
      implicit  real*8(a-h,k-z)
      implicit  integer(i)


      ln=dlog(x**2/y**2)
      return 
      end




c =============================================================

c function dilog

      complex*16 function dilog(xx)
c..
c..   Dilogarithm: dilog(x) = Li_2(x)
c..
      complex*16 xx,wgplg,CSPEN
      
      dilog = CSPEN(xx)
c      dilog = wgplg(1,1,xx)

      end
c =============================================================
C-{{{ FUNCTION WGPLG:

      COMPLEX*16 FUNCTION WGPLG(N,P,X)
c..
c..   Nielson's Generalized Polylogarithm
c..
c..   taken from zwprod.f (R. Hamberg, T. Matsuura and W.L. van Neerven)
c..   (originally from CERNLIB?)
c..
      INTEGER P,P1,NC(10),INDEX(31)
      DOUBLE PRECISION FCT(0:4),SGN(0:4),U(0:4),S1(4,4),C(4,4)
      DOUBLE PRECISION A(0:30,10)
      DOUBLE PRECISION X,X1,H,ALFA,R,Q,C1,C2,B0,B1,B2,ZERO,HALF

      COMPLEX*16 V(0:5),SK,SM

      DATA FCT /1.0D0,1.0D0,2.0D0,6.0D0,24.0D0/
      DATA SGN /1.0D0,-1.0D0,1.0D0,-1.0D0,1.0D0/
      DATA ZERO /0.0D0/, HALF /0.5D0/
      DATA C1 /1.33333 33333 333D0/, C2 /0.33333 33333 3333D0/

      DATA S1(1,1) /1.64493 40668 482D0/
      DATA S1(1,2) /1.20205 69031 596D0/
      DATA S1(1,3) /1.08232 32337 111D0/
      DATA S1(1,4) /1.03692 77551 434D0/
      DATA S1(2,1) /1.20205 69031 596D0/
      DATA S1(2,2) /2.70580 80842 778D-1/
      DATA S1(2,3) /9.65511 59989 444D-2/
      DATA S1(3,1) /1.08232 32337 111D0/
      DATA S1(3,2) /9.65511 59989 444D-2/
      DATA S1(4,1) /1.03692 77551 434D0/

      DATA C(1,1) / 1.64493 40668 482D0/
      DATA C(1,2) / 1.20205 69031 596D0/
      DATA C(1,3) / 1.08232 32337 111D0/
      DATA C(1,4) / 1.03692 77551 434D0/
      DATA C(2,1) / 0.00000 00000 000D0/
      DATA C(2,2) /-1.89406 56589 945D0/
      DATA C(2,3) /-3.01423 21054 407D0/
      DATA C(3,1) / 1.89406 56589 945D0/
      DATA C(3,2) / 3.01423 21054 407D0/
      DATA C(4,1) / 0.00000 00000 000D0/

      DATA INDEX /1,2,3,4,6*0,5,6,7,7*0,8,9,8*0,10/

      DATA NC /24,26,28,30,22,24,26,19,22,17/

      DATA A( 0,1) / .96753 21504 3498D0/
      DATA A( 1,1) / .16607 30329 2785D0/
      DATA A( 2,1) / .02487 93229 2423D0/
      DATA A( 3,1) / .00468 63619 5945D0/
      DATA A( 4,1) / .00100 16274 9616D0/
      DATA A( 5,1) / .00023 20021 9609D0/
      DATA A( 6,1) / .00005 68178 2272D0/
      DATA A( 7,1) / .00001 44963 0056D0/
      DATA A( 8,1) / .00000 38163 2946D0/
      DATA A( 9,1) / .00000 10299 0426D0/
      DATA A(10,1) / .00000 02835 7538D0/
      DATA A(11,1) / .00000 00793 8705D0/
      DATA A(12,1) / .00000 00225 3670D0/
      DATA A(13,1) / .00000 00064 7434D0/
      DATA A(14,1) / .00000 00018 7912D0/
      DATA A(15,1) / .00000 00005 5029D0/
      DATA A(16,1) / .00000 00001 6242D0/
      DATA A(17,1) / .00000 00000 4827D0/
      DATA A(18,1) / .00000 00000 1444D0/
      DATA A(19,1) / .00000 00000 0434D0/
      DATA A(20,1) / .00000 00000 0131D0/
      DATA A(21,1) / .00000 00000 0040D0/
      DATA A(22,1) / .00000 00000 0012D0/
      DATA A(23,1) / .00000 00000 0004D0/
      DATA A(24,1) / .00000 00000 0001D0/

      DATA A( 0,2) / .95180 88912 7832D0/
      DATA A( 1,2) / .43131 13184 6532D0/
      DATA A( 2,2) / .10002 25071 4905D0/
      DATA A( 3,2) / .02442 41559 5220D0/
      DATA A( 4,2) / .00622 51246 3724D0/
      DATA A( 5,2) / .00164 07883 1235D0/
      DATA A( 6,2) / .00044 40792 0265D0/
      DATA A( 7,2) / .00012 27749 4168D0/
      DATA A( 8,2) / .00003 45398 1284D0/
      DATA A( 9,2) / .00000 98586 9565D0/
      DATA A(10,2) / .00000 28485 6995D0/
      DATA A(11,2) / .00000 08317 0847D0/
      DATA A(12,2) / .00000 02450 3950D0/
      DATA A(13,2) / .00000 00727 6496D0/
      DATA A(14,2) / .00000 00217 5802D0/
      DATA A(15,2) / .00000 00065 4616D0/
      DATA A(16,2) / .00000 00019 8033D0/
      DATA A(17,2) / .00000 00006 0204D0/
      DATA A(18,2) / .00000 00001 8385D0/
      DATA A(19,2) / .00000 00000 5637D0/
      DATA A(20,2) / .00000 00000 1735D0/
      DATA A(21,2) / .00000 00000 0536D0/
      DATA A(22,2) / .00000 00000 0166D0/
      DATA A(23,2) / .00000 00000 0052D0/
      DATA A(24,2) / .00000 00000 0016D0/
      DATA A(25,2) / .00000 00000 0005D0/
      DATA A(26,2) / .00000 00000 0002D0/

      DATA A( 0,3) / .98161 02799 1365D0/
      DATA A( 1,3) / .72926 80632 0726D0/
      DATA A( 2,3) / .22774 71490 9321D0/
      DATA A( 3,3) / .06809 08329 6197D0/
      DATA A( 4,3) / .02013 70118 3064D0/
      DATA A( 5,3) / .00595 47848 0197D0/
      DATA A( 6,3) / .00176 76901 3959D0/
      DATA A( 7,3) / .00052 74821 8502D0/
      DATA A( 8,3) / .00015 82746 1460D0/
      DATA A( 9,3) / .00004 77492 2076D0/
      DATA A(10,3) / .00001 44792 0408D0/
      DATA A(11,3) / .00000 44115 4886D0/
      DATA A(12,3) / .00000 13500 3870D0/
      DATA A(13,3) / .00000 04148 1779D0/
      DATA A(14,3) / .00000 01279 3307D0/
      DATA A(15,3) / .00000 00395 9070D0/
      DATA A(16,3) / .00000 00122 9055D0/
      DATA A(17,3) / .00000 00038 2658D0/
      DATA A(18,3) / .00000 00011 9459D0/
      DATA A(19,3) / .00000 00003 7386D0/
      DATA A(20,3) / .00000 00001 1727D0/
      DATA A(21,3) / .00000 00000 3687D0/
      DATA A(22,3) / .00000 00000 1161D0/
      DATA A(23,3) / .00000 00000 0366D0/
      DATA A(24,3) / .00000 00000 0116D0/
      DATA A(25,3) / .00000 00000 0037D0/
      DATA A(26,3) / .00000 00000 0012D0/
      DATA A(27,3) / .00000 00000 0004D0/
      DATA A(28,3) / .00000 00000 0001D0/

      DATA A( 0,4) /1.06405 21184 614 D0/
      DATA A( 1,4) /1.06917 20744 981 D0/
      DATA A( 2,4) / .41527 19325 1768D0/
      DATA A( 3,4) / .14610 33293 6222D0/
      DATA A( 4,4) / .04904 73264 8784D0/
      DATA A( 5,4) / .01606 34086 0396D0/
      DATA A( 6,4) / .00518 88935 0790D0/
      DATA A( 7,4) / .00166 29871 7324D0/
      DATA A( 8,4) / .00053 05827 9969D0/
      DATA A( 9,4) / .00016 88702 9251D0/
      DATA A(10,4) / .00005 36832 8059D0/
      DATA A(11,4) / .00001 70592 3313D0/
      DATA A(12,4) / .00000 54217 4374D0/
      DATA A(13,4) / .00000 17239 4082D0/
      DATA A(14,4) / .00000 05485 3275D0/
      DATA A(15,4) / .00000 01746 7795D0/
      DATA A(16,4) / .00000 00556 7550D0/
      DATA A(17,4) / .00000 00177 6234D0/
      DATA A(18,4) / .00000 00056 7224D0/
      DATA A(19,4) / .00000 00018 1313D0/
      DATA A(20,4) / .00000 00005 8012D0/
      DATA A(21,4) / .00000 00001 8579D0/
      DATA A(22,4) / .00000 00000 5955D0/
      DATA A(23,4) / .00000 00000 1911D0/
      DATA A(24,4) / .00000 00000 0614D0/
      DATA A(25,4) / .00000 00000 0197D0/
      DATA A(26,4) / .00000 00000 0063D0/
      DATA A(27,4) / .00000 00000 0020D0/
      DATA A(28,4) / .00000 00000 0007D0/
      DATA A(29,4) / .00000 00000 0002D0/
      DATA A(30,4) / .00000 00000 0001D0/

      DATA A( 0,5) / .97920 86066 9175D0/
      DATA A( 1,5) / .08518 81314 8683D0/
      DATA A( 2,5) / .00855 98522 2013D0/
      DATA A( 3,5) / .00121 17721 4413D0/
      DATA A( 4,5) / .00020 72276 8531D0/
      DATA A( 5,5) / .00003 99695 8691D0/
      DATA A( 6,5) / .00000 83806 4065D0/
      DATA A( 7,5) / .00000 18684 8945D0/
      DATA A( 8,5) / .00000 04366 6087D0/
      DATA A( 9,5) / .00000 01059 1733D0/
      DATA A(10,5) / .00000 00264 7892D0/
      DATA A(11,5) / .00000 00067 8700D0/
      DATA A(12,5) / .00000 00017 7654D0/
      DATA A(13,5) / .00000 00004 7342D0/
      DATA A(14,5) / .00000 00001 2812D0/
      DATA A(15,5) / .00000 00000 3514D0/
      DATA A(16,5) / .00000 00000 0975D0/
      DATA A(17,5) / .00000 00000 0274D0/
      DATA A(18,5) / .00000 00000 0077D0/
      DATA A(19,5) / .00000 00000 0022D0/
      DATA A(20,5) / .00000 00000 0006D0/
      DATA A(21,5) / .00000 00000 0002D0/
      DATA A(22,5) / .00000 00000 0001D0/

      DATA A( 0,6) / .95021 85196 3952D0/
      DATA A( 1,6) / .29052 52916 1433D0/
      DATA A( 2,6) / .05081 77406 1716D0/
      DATA A( 3,6) / .00995 54376 7280D0/
      DATA A( 4,6) / .00211 73389 5031D0/
      DATA A( 5,6) / .00047 85947 0550D0/
      DATA A( 6,6) / .00011 33432 1308D0/
      DATA A( 7,6) / .00002 78473 3104D0/
      DATA A( 8,6) / .00000 70478 8108D0/
      DATA A( 9,6) / .00000 18278 8740D0/
      DATA A(10,6) / .00000 04838 7492D0/
      DATA A(11,6) / .00000 01303 3842D0/
      DATA A(12,6) / .00000 00356 3769D0/
      DATA A(13,6) / .00000 00098 7174D0/
      DATA A(14,6) / .00000 00027 6586D0/
      DATA A(15,6) / .00000 00007 8279D0/
      DATA A(16,6) / .00000 00002 2354D0/
      DATA A(17,6) / .00000 00000 6435D0/
      DATA A(18,6) / .00000 00000 1866D0/
      DATA A(19,6) / .00000 00000 0545D0/
      DATA A(20,6) / .00000 00000 0160D0/
      DATA A(21,6) / .00000 00000 0047D0/
      DATA A(22,6) / .00000 00000 0014D0/
      DATA A(23,6) / .00000 00000 0004D0/
      DATA A(24,6) / .00000 00000 0001D0/

      DATA A( 0,7) / .95064 03218 6777D0/
      DATA A( 1,7) / .54138 28546 5171D0/
      DATA A( 2,7) / .13649 97959 0321D0/
      DATA A( 3,7) / .03417 94232 8207D0/
      DATA A( 4,7) / .00869 02788 3583D0/
      DATA A( 5,7) / .00225 28408 4155D0/
      DATA A( 6,7) / .00059 51608 9806D0/
      DATA A( 7,7) / .00015 99561 7766D0/
      DATA A( 8,7) / .00004 36521 3096D0/
      DATA A( 9,7) / .00001 20747 4688D0/
      DATA A(10,7) / .00000 33801 8176D0/
      DATA A(11,7) / .00000 09563 2476D0/
      DATA A(12,7) / .00000 02731 3129D0/
      DATA A(13,7) / .00000 00786 6968D0/
      DATA A(14,7) / .00000 00228 3195D0/
      DATA A(15,7) / .00000 00066 7205D0/
      DATA A(16,7) / .00000 00019 6191D0/
      DATA A(17,7) / .00000 00005 8018D0/
      DATA A(18,7) / .00000 00001 7246D0/
      DATA A(19,7) / .00000 00000 5151D0/
      DATA A(20,7) / .00000 00000 1545D0/
      DATA A(21,7) / .00000 00000 0465D0/
      DATA A(22,7) / .00000 00000 0141D0/
      DATA A(23,7) / .00000 00000 0043D0/
      DATA A(24,7) / .00000 00000 0013D0/
      DATA A(25,7) / .00000 00000 0004D0/
      DATA A(26,7) / .00000 00000 0001D0/

      DATA A( 0,8) / .98800 01167 2229D0/
      DATA A( 1,8) / .04364 06760 9601D0/
      DATA A( 2,8) / .00295 09117 8278D0/
      DATA A( 3,8) / .00031 47780 9720D0/
      DATA A( 4,8) / .00004 31484 6029D0/
      DATA A( 5,8) / .00000 69381 8230D0/
      DATA A( 6,8) / .00000 12464 0350D0/
      DATA A( 7,8) / .00000 02429 3628D0/
      DATA A( 8,8) / .00000 00504 0827D0/
      DATA A( 9,8) / .00000 00109 9075D0/
      DATA A(10,8) / .00000 00024 9467D0/
      DATA A(11,8) / .00000 00005 8540D0/
      DATA A(12,8) / .00000 00001 4127D0/
      DATA A(13,8) / .00000 00000 3492D0/
      DATA A(14,8) / .00000 00000 0881D0/
      DATA A(15,8) / .00000 00000 0226D0/
      DATA A(16,8) / .00000 00000 0059D0/
      DATA A(17,8) / .00000 00000 0016D0/
      DATA A(18,8) / .00000 00000 0004D0/
      DATA A(19,8) / .00000 00000 0001D0/

      DATA A( 0,9) / .95768 50654 6350D0/
      DATA A( 1,9) / .19725 24967 9534D0/
      DATA A( 2,9) / .02603 37031 3918D0/
      DATA A( 3,9) / .00409 38216 8261D0/
      DATA A( 4,9) / .00072 68170 7110D0/
      DATA A( 5,9) / .00014 09187 9261D0/
      DATA A( 6,9) / .00002 92045 8914D0/
      DATA A( 7,9) / .00000 63763 1144D0/
      DATA A( 8,9) / .00000 14516 7850D0/
      DATA A( 9,9) / .00000 03420 5281D0/
      DATA A(10,9) / .00000 00829 4302D0/
      DATA A(11,9) / .00000 00206 0784D0/
      DATA A(12,9) / .00000 00052 2823D0/
      DATA A(13,9) / .00000 00013 5066D0/
      DATA A(14,9) / .00000 00003 5451D0/
      DATA A(15,9) / .00000 00000 9436D0/
      DATA A(16,9) / .00000 00000 2543D0/
      DATA A(17,9) / .00000 00000 0693D0/
      DATA A(18,9) / .00000 00000 0191D0/
      DATA A(19,9) / .00000 00000 0053D0/
      DATA A(20,9) / .00000 00000 0015D0/
      DATA A(21,9) / .00000 00000 0004D0/
      DATA A(22,9) / .00000 00000 0001D0/

      DATA A( 0,10) / .99343 65167 1347D0/
      DATA A( 1,10) / .02225 77012 6826D0/
      DATA A( 2,10) / .00101 47557 4703D0/
      DATA A( 3,10) / .00008 17515 6250D0/
      DATA A( 4,10) / .00000 89997 3547D0/
      DATA A( 5,10) / .00000 12082 3987D0/
      DATA A( 6,10) / .00000 01861 6913D0/
      DATA A( 7,10) / .00000 00317 4723D0/
      DATA A( 8,10) / .00000 00058 5215D0/
      DATA A( 9,10) / .00000 00011 4739D0/
      DATA A(10,10) / .00000 00002 3652D0/
      DATA A(11,10) / .00000 00000 5082D0/
      DATA A(12,10) / .00000 00000 1131D0/
      DATA A(13,10) / .00000 00000 0259D0/
      DATA A(14,10) / .00000 00000 0061D0/
      DATA A(15,10) / .00000 00000 0015D0/
      DATA A(16,10) / .00000 00000 0004D0/
      DATA A(17,10) / .00000 00000 0001D0/

      IF(N .LT. 1 .OR. N .GT. 4 .OR. P .LT. 1 .OR. P .GT. 4 .OR.
     1   N+P .GT. 5) THEN
       WGPLG=ZERO
       PRINT 1000, N,P
       RETURN
      END IF
      IF(X .EQ. SGN(0)) THEN
       WGPLG=S1(N,P)
       RETURN
      END IF

      IF(X .GT. FCT(2) .OR. X .LT. SGN(1)) THEN
       X1=SGN(0)/X
       H=C1*X1+C2
       ALFA=H+H
       V(0)=SGN(0)
       V(1)=LOG(DCMPLX(-X,ZERO))
       DO 33 L = 2,N+P
   33  V(L)=V(1)*V(L-1)/L
       SK=ZERO
       DO 34 K = 0,P-1
       P1=P-K
       R=X1**P1/(FCT(P1)*FCT(N-1))
       SM=ZERO
       DO 35 M = 0,K
       N1=N+K-M
       L=INDEX(10*N1+P1-10)
       B1=ZERO
       B2=ZERO
       DO 31 I = NC(L),0,-1
       B0=A(I,L)+ALFA*B1-B2
       B2=B1
   31  B1=B0
       Q=(FCT(N1-1)/FCT(K-M))*(B0-H*B2)*R/P1**N1
   35  SM=SM+V(M)*Q
   34  SK=SK+SGN(K)*SM
       SM=ZERO
       DO 36 M = 0,N-1
   36  SM=SM+V(M)*C(N-M,P)
       WGPLG=SGN(N)*SK+SGN(P)*(SM+V(N+P))
       RETURN
      END IF

      IF(X .GT. HALF) THEN
       X1=SGN(0)-X
       H=C1*X1+C2
       ALFA=H+H
       V(0)=SGN(0)
       U(0)=SGN(0)
       V(1)=LOG(DCMPLX(X1,ZERO))
       U(1)=LOG(X)
       DO 23 L = 2,P
   23  V(L)=V(1)*V(L-1)/L
       DO 26 L = 2,N
   26  U(L)=U(1)*U(L-1)/L
       SK=ZERO
       DO 24 K = 0,N-1
       P1=N-K
       R=X1**P1/FCT(P1)
       SM=ZERO
       DO 25 M = 0,P-1
       N1=P-M
       L=INDEX(10*N1+P1-10)
       B1=ZERO
       B2=ZERO
       DO 12 I = NC(L),0,-1
       B0=A(I,L)+ALFA*B1-B2
       B2=B1
   12  B1=B0
       Q=SGN(M)*(B0-H*B2)*R/P1**N1
   25  SM=SM+V(M)*Q
   24  SK=SK+U(K)*(S1(P1,P)-SM)
       WGPLG=SK+SGN(P)*U(N)*V(P)
       RETURN
      END IF

      L=INDEX(10*N+P-10)
      H=C1*X+C2
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO 11 I = NC(L),0,-1
      B0=A(I,L)+ALFA*B1-B2
      B2=B1
   11 B1=B0
      WGPLG=(B0-H*B2)*X**P/(FCT(P)*P**N)
      RETURN
 1000 FORMAT(/' ***** CERN SUBROUTINE WGPLG ... ILLEGAL VALUES',
     1        '   N = ',I3,'   P = ',I3)
      END

C-}}}


c =============================================================
c  -function DCLAUS:
c.. Clausen's function from CERNLIB

      FUNCTION DCLAUS(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      DIMENSION A(0:8),B(0:13)
 
      PARAMETER (R1 = 1, HF =R1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (PI2 = 2*PI, PIH = PI/2, RPIH = 2/PI)
 
      DATA A( 0) / 0.02795 28319 73575 6613D0/
      DATA A( 1) / 0.00017 63088 74389 8116D0/
      DATA A( 2) / 0.00000 12662 74146 1157D0/
      DATA A( 3) / 0.00000 00117 17181 8134D0/
      DATA A( 4) / 0.00000 00001 23006 4129D0/
      DATA A( 5) / 0.00000 00000 01395 2729D0/
      DATA A( 6) / 0.00000 00000 00016 6908D0/
      DATA A( 7) / 0.00000 00000 00000 2076D0/
      DATA A( 8) / 0.00000 00000 00000 0027D0/
 
      DATA B( 0) / 0.63909 70888 57265 341D0/
      DATA B( 1) /-0.05498 05693 01851 716D0/
      DATA B( 2) /-0.00096 12619 45950 606D0/
      DATA B( 3) /-0.00003 20546 86822 550D0/
      DATA B( 4) /-0.00000 13294 61695 426D0/
      DATA B( 5) /-0.00000 00620 93601 824D0/
      DATA B( 6) /-0.00000 00031 29600 656D0/
      DATA B( 7) /-0.00000 00001 66351 954D0/
      DATA B( 8) /-0.00000 00000 09196 527D0/
      DATA B( 9) /-0.00000 00000 00524 004D0/
      DATA B(10) /-0.00000 00000 00030 580D0/
      DATA B(11) /-0.00000 00000 00001 820D0/
      DATA B(12) /-0.00000 00000 00000 110D0/
      DATA B(13) /-0.00000 00000 00000 007D0/
 
      V=MOD(ABS(X),PI2)
      S=SIGN(R1,X)
      IF(V .GT. PI) THEN
       V=PI2-V
       S=-S
      ENDIF
      IF(V .EQ. 0 .OR. V .EQ. PI) THEN
       H=0
      ELSEIF(V .LT. PIH) THEN
       U=RPIH*V
       H=2*U**2-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 8,0,-1
       B0=A(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=V*(1-LOG(V)+HF*V**2*(B0-H*B2))
      ELSE
       U=RPIH*V-2
       H=2*U**2-1
       ALFA=H+H
       B1=0
       B2=0
       DO 2 I = 13,0,-1
       B0=B(I)+ALFA*B1-B2
       B2=B1
    2  B1=B0
       H=(PI-V)*(B0-H*B2)
      ENDIF
      DCLAUS=S*H
      RETURN
      END

* =====================================================================
*
*  Function returns Alpha_s/(4*pi) given in MSbar with 6 active quark 
*  flavours in the SM (needed in Coulomb corrections)
*
* (from Q_0 with alpha_s(Q_0,nf=5) do three loop running within MSbar+SM+nf=5
*  up to the relevant scale, then add top-Quark effects at one loop)
*  
*  Using three-loop running for MSbar,SM with two-loop decoupling 
*------------------------------------------------------------------------
*  All formulas for decoupling coefficients are taken from 

*  Harlander, Mihaila, Steinhauser - hep-ph/0509048
*  'Two-loop matching coefficients for the strong coupling in the MSSM'


 
      FUNCTION ASMS (Qinp2)

#include "../util/DMNLO_Model.h"

      double precision ASMS,Qinp2
      double precision ASMS5Q
      double precision zetaDC(2)
      
* external functions      
      double precision AS
      
      external ASzetaDecCoe_Mo

* three-loop running of alpha_s^{MS,SM,nf=5} returns alpha/pi	
      ASMS5Q = 4d0*AS( Qinp2, MZ**2, AlfasMZ/(4d0*pi), 5, 3, 5)

* one-loop SM decoupling of alpha_s^{MS,SM,nf=5}/pi
* adding top-quark in one step

* decoupling coefficients at 2-loop
      call ASzetaDecCoe_Mo(Qinp2,zetaDC)                                  

* running of alpha_s^{MS,SM,nf=6} returns alpha/(4*pi)	
      ASMS = 1d0/4d0*ASMS5Q*(1d0 - ASMS5Q*zetaDC(1)+ASMS5Q**2*(2d0*zetaDC(1)**2-zetaDC(2)))   !factor 2 from correctly solving Eq.19 in 0509048v1 for alpha_s(n_f)

      RETURN
      END


* =====================================================================
*  Subroutine for top quark two-loop decoupling coefficients for alpha_s
*
      SUBROUTINE ASzetaDecCoe_Mo(Qdec2,zetaDC)
      
#include "../util/DMNLO_Model.h"

      double precision zetaDC(2),Qdec2
      double precision Mtop,Mtop2

* Common mass of colored SM-particles (for one-loop matching)
      Mtop = MTos
      Mtop2 = MTos**2


*taken from 0509048 Eq.20 with switch (DRbar->MSbar)

* one-loop matching with full mass dependance (all quarks except for the top are set to zero mass)
      zetaDC(1) = -1d0/6d0*dlog(Qdec2/Mtop2) 
* two-loop matching with full mass dependance (all quarks except for the top are set to zero mass)
      zetaDC(2) = -7d0/24d0-19d0/24d0*dlog(Qdec2/Mtop2) + 1d0/36d0*dlog(Qdec2/Mtop2)**2      
      
      END
