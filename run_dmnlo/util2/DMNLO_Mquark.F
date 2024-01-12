! ======= Bottom quark pole mass from mb^MS_SM (mb) ========

      function MB2OS(mbms,Q,isch)

#include "../util/DMNLO_Model.h"
#include "../util/DMNLO_RenScheme.h"

* functions
	double precision AS,delOS
	double precision MB2OS

* variables
      integer isch
	double precision ASpiMB,LLog
	double precision mbms,Q,delqm
	double precision fOSMS, fOSDR
      double precision qmass(4)

	DATA qmass / 2d-3, 45d-4, 1d-1, 1.27d0/
      
!**** Transformation to MB_DR(Q) in the Standard Model
! from Baer, Ferrandis, Melnikov and Tata hep-ph/0207126

      delqm = 0d0
      do i=1,4
        delqm = delqm + delOS(qmass(i)/mbms)
      enddo
      
	ASpiMB = (4d0)*AS(Q**2, MZ**2, AlfasMZ/(4d0*pi), 5, 3, 5)
	LLog = dlog(mbms**2/Q**2)
      fOSMS = 1d0 + ASpiMB*(LLog - 4d0/3d0) + ASpiMB**2*(-11d0/24d0*LLog**2 + 197d0/72d0*LLog 
     &       -187d0/32d0 + zeta3/6d0 - pi**2/9d0*dlog(2d0) - pi**2/9d0- delqm)

      fOSDR = 1d0 + ASpiMB*(5d0/3d0 - LLog) + ASpiMB**2*( 35d0/24d0*LLog**2 - 395d0/72d0*LLog 
     &       +2159d0/288d0 - zeta3/6d0 + pi**2/9d0*dlog(2d0) + pi**2/9d0+ delqm)
           
      if (isch.eq.0) then      
        MB2OS = mbms/fOSMS
      else if (isch.eq.1) then
! doesn't work... need Alpha_s in DRbar (?)
        MB2OS = mbms*fOSDR
!        MB2OS = -99d0
      else 
        MB2OS = -99d0
      end if
      
      return
      end

! Auxillary function to compute OS bottom quark mass
      function delOS(mbrat)

      double precision delOS,mbrat

#include "../util/DMNLO_Model.h"

      delOS = pi**2/6d0*mbrat*(1d0+mbrat**2) - mbrat**2

      end


! ======= Top quark DR mass @ Q**2 in the MSSM ========

      function MTDR(mtosinp,Q)

#include "../util/DMNLO_Model.h"
#include "../util/DMNLO_RenScheme.h"

! functions
	double precision MTDR

! variables
	double precision mtosinp,Q,mtit,test

	mtit = mtosinp

	do i=1,3
	 
	 call DMNLO_CalcPifth(mtit**2,Q,3,3)
	 dMf(3,3) = (1d0/2d0)*(PifVLth(3,3) + PifVRth(3,3)) + (PifSLth(3,3) + PifSRth(3,3))/(2d0*mtit)
	 mtit = mtosinp/(1d0-dMf(3,3))

      enddo

      MTDR = mtit

      return
      end

! =====================================================================
!
!  Function MBDR returns MB^{DR,MSSM}(Q) from MB^{MS,SM}(Q0)
!  given in the Standard Model with inf-flavours
!  
! For the first step - MB^{MS,SM}(Q0) -> MB^{MS,SM}(Q)
! using up-to 4-loop running (should use 3-loop)
!
!  All formulas for beta function coefficients and anomalous dimensions are taken from
!
!  K.G. Chetyrkin, Phys. Lett. B 404 (1997) 161
!  S.A. Larin, T. van Ritbergen, and J.A.M. Vermaseren, Phys. Lett. B 405 (1997) 327
!  Chetyrkin, Kuhn, Steinhauser hep-ph/0004189
!
! For the second step - MB^{MS,SM}(Q0) -> MB^{MS,SM}(Q)
! using 2-loop decoupling
!
!  DR & MSSM decoupling references are
!
!  Harlander, Kant, Mihaila, Steinhauser         - hep-ph/0607240
!  'Dimensional Reduction applied to QCD at three loops'
!
!  Bauer, Mihaila, Salomon         - arXiv: 0810.5101 [hep-ph]
!  'Matching coefficients for alpha_s and m_b to O(a^2) in the MSSM'
 

! ---------------------------------------------------------------------

      FUNCTION MBDR(MBinp,Q02,Q2,inf,nloop)
      implicit none
#include "../util/DMNLO_Model.h"
      
      integer inf, nloop
      double precision MBinp, Q02, Q2
      double precision ASMS50, ASMS5Q,cQ,c0
      
      double precision Ltl,Lgl,Lqb1,Lqb2
      double precision MBMSQ,MBDRQ,MBDR
      double precision ASEV,zetaEVb,ASDRQ,MBzetaDC(2)

! external functions & subroutines      
      double precision Ccoef,AS,ASDR

! shortcuts
      Ltl = dlog(Q2/MT2)
      Lgl = dlog(Q2/MGl2)
      Lqb1 = dlog(Q2/MSf2(1,4,3))
      Lqb2 = dlog(Q2/MSf2(2,4,3))

! three-loop running in the SM      
      ASMS50 = 4d0*AS( Q02, MZ**2, AlfasMZ/(4d0*pi), inf, 3, 5)
      ASMS5Q = 4d0*AS( Q2,  MZ**2, AlfasMZ/(4d0*pi), inf, 3, 5)
      
      cQ = Ccoef(ASMS5Q,inf,nloop)
      c0 = Ccoef(ASMS50,inf,nloop)

! based on Eq.9 in hep-ph/0004189 (checked with RunDec function mMS2mMS)           
      MBMSQ = MBinp*cQ/c0

! MS to DR conversion (attention evanescent coupling needed)
! Eq.(4.2) from 0810.5101 
      ASDRQ = 4d0*ASDR(Q2)

      zetaEVb = 1d0 + ASDRQ*(-TF*Ltl/2d0 + CV/4d0*(2d0 + Lgl + MSf2(1,4,3)*(Lgl-Lqb1)/(MGl2 - MSf2(1,4,3))
     &   + MSf2(2,4,3)*(Lgl-Lqb2)/(MGl2 - MSf2(2,4,3)) + CF/4d0*( -3d0 - 2d0*Lgl +
     & (-1d0 - 2d0*Lgl + 2d0*Lqb1 + (Lqb1-Lgl)*MSf2(1,4,3)/(MGl2-MSf2(1,4,3)))*MSf2(1,4,3)/(MGl2-MSf2(1,4,3)) +
     & (-1d0 - 2d0*Lgl + 2d0*Lqb2 + (Lqb2-Lgl)*MSf2(2,4,3)/(MGl2-MSf2(2,4,3)))*MSf2(2,4,3)/(MGl2-MSf2(2,4,3)) )))

      ASEV = zetaEVb*ASDRQ

! Eq.(4.4) from hep-ph/0607240
      MBDRQ = MBMSQ*(1d0 - ASEV*CF/(4d0) + ASMS5Q**2*11d0/192d0*CV*CF - ASMS5Q*ASEV*(CF**2/4d0 + CV*CF*3d0/32d0) 
     &        + ASEV**2*(CF**2*3d0/32d0 + 1d0/32d0*CF*TF))

! MSSM decoupling coefficients (from 0810.5101 Mathematica file)
      call MBzetaDecCoe(Q2,MBDRQ,MBzetaDC)      

      MBDR = MBDRQ*(1d0 - ASDRQ*MBzetaDC(1) + ASDRQ**2*(2d0*MBzetaDC(1)**2 - MBzetaDC(2)))

      RETURN
      END

! =====================================================================
!  Subroutine for SUSY two-loop decoupling coefficients for m_b
!
      SUBROUTINE MBzetaDecCoe(Qdec2,mbinp,zetaDC)
      implicit none
#include "../util/DMNLO_Model.h"

      double precision zetaDC(2),Qdec2,mbinp
	double precision LogGl,LogB1,LogB2,msb12,msb22

      double precision mssmzeta2lmb

! shortcut logs & masses
      LogGl = dlog(Qdec2/MGl**2)
      LogB1 = dlog(Qdec2/MSf2(1,4,3))
      LogB2 = dlog(Qdec2/MSf2(2,4,3))

      msb12 = MSf2(1,4,3)
      msb22 = MSf2(2,4,3)

! one-loop matching with full mass dependance (Eq.4.3 in 0810.5101)
      zetaDC(1) = CF*(
     &      -(1d0+logB1)*msb12/(4d0*(msb12 - MGl2)) + 
     &      ((3d0+2d0*logB1)*msb12**2 - (3d0+2d0*logGl)*MGl2**2)/(16d0*(msb12 - MGl2)**2) 
     &      -(1d0+logB2)*msb22/(4d0*(msb22 - MGl2)) + 
     &      ((3d0+2d0*logB2)*msb22**2 - (3d0+2d0*logGl)*MGl2**2)/(16d0*(msb22 - MGl2)**2) 
     &      + USf(1,1,4,3)*USf(1,2,4,3)*MGl*(msb12*logB1 - MGl2*logGl)/(2d0*mbinp*(msb12 - MGl2)) 
     &      - USf(1,1,4,3)*USf(1,2,4,3)*MGl*(msb22*logB2 - MGl2*logGl)/(2d0*mbinp*(msb22 - MGl2)) 
     &      )

! two-loop matching - exact from L.Mihaila's Mathematica code
      zetaDC(2) =-mssmzeta2lmb(dsqrt(Qdec2),MGl,MTOS,MSf(1,3,3),MSf(2,3,3),Xf(3,3),MSf(1,4,3),MSf(2,4,3),Xf(4,3))

      END


! ---------------------------------------------------------------------
! Coefficient for MS-bar running mass from RunDec hep-ph/0004189 Eq.10
!
      FUNCTION Ccoef(xx,inf,nloop)
      implicit none
      integer inf,nloop,nf
      double precision xx,cf,ca,tf,N,zeta3,zeta4,zeta5
      
      double precision beta0, beta1, beta2, beta3, GM0, GM1, GM2, GM3
      double precision b1, b2, b3, c0, c1, c2, c3
      
      double precision Ccoef
      
      parameter (zeta3 = 1.2020569031595942D0)
      parameter (zeta4 = 1.0823232337111381D0)
      parameter (zeta5 = 1.03692775514337D0)

      N = 3d0
      CA = 3d0
      CF = 4/3d0
      TF = 1/2d0
      
      nf=inf
      
! beta function coefficients (difference of 1/4, 1/16 ... to AS function)      
      beta0 = 1/4d0*(11d0 - 2/3d0*nf)
      beta1 = 1/16d0*(102d0 - 38/3d0*nf)
      beta2 = 1/64d0*(2857d0/2d0 - 5033d0/18d0*nf + 325d0/54d0*nf**2)
      beta3 = 1/256d0*(149753d0/6d0 + 3564d0*zeta3 - (1078361d0/162d0 + 6508d0/27d0*zeta3)*nf
     &        + (50065d0/162d0 + 6472d0/81d0*zeta3)*nf**2 + 1093d0/729d0*nf**3)

! anom. dim beta coefficients (difference of 1/4, 1/16 ... to AS function)     
       GM0 = 1/4d0*3d0*CF
       GM1 = 1/16d0*(3d0/2d0*CF**2 + 97d0/6d0*CF*CA - 10d0/3d0*CF*TF*nf)
       GM2 = 1/64d0*(129d0/2d0*CF**3 - 129d0/4d0*CF**2*CA + 11413d0/108d0*CF*CA**2 
     &    + (-46d0 + 48d0*zeta3)*CF**2*TF*nf + (-556d0/27d0 - 48d0*zeta3)*CF*CA*TF*nf
     &    - 140d0/27d0*CF*TF**2*nf**2)
       GM3 = 1/256d0*((-1261d0/8d0 - 336d0*zeta3)*CF**4 + (15349d0/12d0 + 316d0*zeta3)*CA*CF**3
     &    + (-34045d0/36d0 - 152d0*zeta3 + 440d0*zeta5)*CF**2*CA**2  
     &    + (70055d0/72d0 + 1418d0/9d0*zeta3 - 440d0*zeta5)*CF*CA**3  
     &    + (-280d0/3d0 + 552d0*zeta3 - 480d0*zeta5)*CF**3*TF*nf
     &    + (-8819d0/27d0 + 368d0*zeta3 - 264d0*zeta4 + 80d0*zeta5)*CF**2*CA*TF*nf
     &    + (-65459d0/162d0 - 2684d0/3d0*zeta3 + 264d0*zeta4 + 400d0*zeta5)*CF*CA**2*TF*nf
     &    + (304d0/27d0 - 160d0*zeta3 + 96d0*zeta4)*CF**2*TF**2*nf**2
     &    + (1342d0/81d0 + 160d0*zeta3 - 96d0*zeta4)*CF*CA*TF**2*nf**2
     &    + (-664d0/81d0 + 128d0/9d0*zeta3)*CF*TF**3*nf**3
     &    + (-32d0 + 240d0*zeta3)*(N**2 - 1d0)*(N**2 + 6d0)/48d0
     &    + (64d0 - 480d0*zeta3)*nf*(N**2 - 1d0)*(N**4 - 6d0*N**2 + 18d0)/(96d0*N**3) )

! c & b coefficients from RunDec paper
      b1 = beta1/beta0      
      b2 = beta2/beta0      
      b3 = beta3/beta0      
            
      c0 = GM0/beta0      
      c1 = GM1/beta0      
      c2 = GM2/beta0      
      c3 = GM3/beta0      
      
      if (nloop.eq.1) then
            
         Ccoef = xx**c0
         
      else if (nloop.eq.2) then
            
         Ccoef = xx**c0*(1d0 + (c1 - b1*c0)*xx)

      else if (nloop.eq.3) then
            
         Ccoef = xx**c0*(1d0 + (c1 - b1*c0)*xx + 1d0/2d0*((c1 - b1*c0)**2 + c2 - b1*c1 + b1**2*c0 - b2*c0)*xx**2)

      else if (nloop.eq.4) then
            
         Ccoef = xx**c0*(1d0 + (c1 - b1*c0)*xx + 1d0/2d0*((c1 - b1*c0)**2 + c2 - b1*c1 + b1**2*c0 - b2*c0)*xx**2
     &      + (1/6d0*(c1 - b1*c0)**3 + 1/2d0*(c1 - b1*c0)*(c2 - b1*c1 + b1**2*c0 - b2*c0)
     &      + 1/3d0*(c3 - b1*c2 + b1**2*c1 - b2*c1 - b1**3*c0 + 2d0*b1*b2*c0 - b3*c0))*xx**3)

      endif

      RETURN
      END

! --- temp functions correcting arguments mismatch between Fortran and Mathematica from L.Mihaila

      function B0finTMP(M1,M2,M3,mu)
      implicit none
      double precision M1,M2,M3,mu
      double precision B0finTMP, B0fin
      
      B0finTMP = B0fin(dsqrt(M1),M2,M3,mu)
      
      return
      end


      function TfinTMP(M1,M2,M3,mu)
      implicit none
      double precision M1,M2,M3,mu
      double precision TfinTMP, Tfin

      TfinTMP = Tfin(M1,M2,M3)

      return
      end

      function TfinZTMP(M1,M2,mu)
      implicit none
      double precision M1,M2,mu
      double precision TfinZTMP, Tfin0
      
      TfinZTMP = Tfin0(M1,M2)
      
      return
      end


      function mMS2mMS(mmu0, alsmu0, alsmu, nf, nloops)
      implicit none
      double precision  mmu0, alsmu0, alsmu,c0,cmu,mMS2mMS,Ccoef,Pi
      integer nf, nloops
   
      pi = 3.14159265358d0
      cmu = Ccoef(alsmu/Pi,nf,nloops)
      c0 = Ccoef(alsmu0/Pi,nf,nloops)
   
! based on Eq.9 in hep-ph/0004189 (checked with RunDec function mMS2mMS)           
      mMS2mMS = mmu0*cmu/c0;

      end 
      