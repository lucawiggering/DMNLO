! =========================================================================
	  double precision function NeuQ2qx_M2VertexH()

	  implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "NeuQ2qx_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"


	  double complex A1vL,A1vR,A1sL,A1sR,A2sL,A2sR,A3vL,A3vR,A3sL,A3sR
	  double complex A1Lc,A1Rc,A2Lc,A2Rc,A3Lc,A3Rc
	  double complex FvL,FvR,FsL,FsR
	  double complex FLc,FRc
	  double complex Gscalar,Gc
	  double complex HtreeL,HtreeR,HLc,HRc

	  double complex PvL,PvR,PsL,PsR,P

	  double precision mq1,mq1c,mq2,mq2c,mq3,mq3c,wquark

	  double complex g0,g0L,g0R,g1,g1L,g1R,g2,g2L,g2R
	  double precision Mloop0,Mloop1,Mloop2

	  double precision vertAs, vertAt, vertAu, vertFs, vertGt
	  double precision MssH, MttH, MuuH, MsuH, MstH, MtsH, MtuH, MusH, MutH
	  

	  integer i,k,l

	  integer sumU, sumUc, sumT, sumTc
	  integer sumUend, sumUcend, sumTend, sumTcend

	  sumU = 0
	  sumUc = 0
	  sumT = 0
	  sumTc = 0

	  if ((xtype.eq.5).or.(xtype.eq.6)) then

	  sumTend = 1
	  sumTcend = 1
	  sumUend = 3
	  sumUcend = 3

	  elseif (xtype.eq.7) then

	  sumTend = 0
	  sumTcend = 0
	  sumUend = 3
	  sumUcend = 3

	  elseif (xtype.eq.8) then

	  sumTend = 1
	  sumTcend = 1
	  sumUend = 1
	  sumUcend = 1

	  endif


	  NeuQ2qx_M2VertexH = 0d0



! ------------------------- Correction A-coupling (Neutralino - Sfermion - Fermion) -------------------------
!                                                 s-channel

	  vertAs = 0d0

! --------- Gluon correction to Neutralino-Sfermion-Fermion vertex in S-channel ---------

	  call myCget(m1**2,m2**2,sqrtS**2,Mf(pstype,psgen)**2,m2**2,0d0)
	  call myBget(m2**2,m2**2,0d0)

	  g0 = -GS
	  g1L = AL(chimeigen,sfmeigen,sftype,sfgen)
	  g1R = AR(chimeigen,sfmeigen,sftype,sfgen)
	  g2L = -GS
	  g2R = -GS

	  Mloop0 = Mf(pstype,psgen)

	  A1sL = (g0*g1R*g2L*Mloop0*m1*(2*xc0 + xc1)*(-1d0) + g0*g1L*g2L*(xb0 + Mloop0**2*xc0 +
     & sman*xc2 + 2*m1**2*xc1 + 2*m1**2*xc2 - 2*m2**2*xc2))/(16.*Pi**2)

	  A1sR = (g0*g1L*g2R*Mloop0*m1*(2*xc0 + xc1)*(-1d0) + g0*g1R*g2R*(xb0 + Mloop0**2*xc0 +
     & sman*xc2 + 2*m1**2*xc1 + 2*m1**2*xc2 - 2*m2**2*xc2))/(16.*Pi**2)

	  A1vL = (g0*g1R*g2L*Mloop0*(xc0 - xc2) - g0*g1L*g2L*m1*(2*xc2 + xc1))/(16.*Pi**2)

	  A1vR = (g0*g1L*g2R*Mloop0*(xc0 - xc2) - g0*g1R*g2R*m1*(2*xc2 + xc1))/(16.*Pi**2)



! --------- Gluino correction to Neutralino-Sfermion-Fermion vertex in S-channel ---------

	  do k=1,2

	  call myCget(m1**2,m2**2,sqrtS**2,MSf(k,sftype,sfgen)**2,Mf(pstype,psgen)**2,MGl**2)
	  call myBget(m2**2,Mf(pstype,psgen)**2,MGl**2)

	  g0L =  sqrt2*GS*USf(sfmeigen,2,sftype,sfgen)
	  g0R = -sqrt2*GS*USf(sfmeigen,1,sftype,sfgen)
	  g1L = dconjg(AR(chimeigen,k,sftype,sfgen))
	  g1R = dconjg(AL(chimeigen,k,sftype,sfgen))
	  g2L =  sqrt2*GS*USf(k,2,sftype,sfgen)
	  g2R = -sqrt2*GS*USf(k,1,sftype,sfgen)

	  Mloop0 = MSf(k,sftype,sfgen)
  	  Mloop1 = Mf(pstype,psgen)
	  Mloop2 = MGl

	  A1sL = A1sL - ((g0L*g1L*g2L*Mloop1*Mloop2 + g0L*g1R*g2L*Mloop2*m1)*xc0 + g0R*g1L*g2L*(xb0 + Mloop0**2*xc0) +
     & g0L*g1R*g2L*Mloop2*m1*xc1 + g0R*g1R*g2L*Mloop1*m1*xc1 + g0R*g1L*g2L*(m1**2*xc1 + sman*xc2))/(16.*Pi**2)

	  A1sR = A1sR - ((g0R*g1R*g2R*Mloop1*Mloop2 + g0R*g1L*g2R*Mloop2*m1)*xc0 + g0L*g1R*g2R*(xb0 + Mloop0**2*xc0) +
     & g0R*g1L*g2R*Mloop2*m1*xc1 + g0L*g1L*g2R*Mloop1*m1*xc1 + g0L*g1R*g2R*(m1**2*xc1 + sman*xc2))/(16.*Pi**2)

	  A1vL = A1vL - ((g0R*g1R*g2L*Mloop1 + g0R*g1L*g2L*m1)*xc0 + g0L*g1R*g2L*Mloop2*xc2 + g0R*g1R*g2L*Mloop1*xc2 + g0R*g1L*g2L*m1*(xc1 + xc2))/(16.*Pi**2)
	  A1vR = A1vR - ((g0L*g1L*g2R*Mloop1 + g0L*g1R*g2R*m1)*xc0 + g0R*g1L*g2R*Mloop2*xc2 + g0L*g1L*g2R*Mloop1*xc2 + g0L*g1R*g2R*m1*(xc1 + xc2))/(16.*Pi**2)

	  enddo


! propagator
	  PvL = 1d0
	  PvR = 1d0
	  PsL = Mf(pstype,psgen)
	  PsR = Mf(pstype,psgen)

! Higgs-Fermion-Fermion
	  FvL = 0d0
	  FvR = 0d0
	  FsL = FL(xtype,pstype,psgen)
	  FsR = FR(xtype,pstype,psgen)

! mass of propagator & width
	  mq1 = Mf(pstype,psgen)
	  wquark = 0d0


! --- tree-level - s-channel ---

	    A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
	    A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
	    FLc = dconjg(FL(xtype,pstype,psgen))
	    FRc = dconjg(FR(xtype,pstype,psgen))

	    mq1c = Mf(pstype,psgen)

	    vertAs = vertAs + CF*NC*2d0*MssH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A1Lc,A1Rc,FLc,FRc,mq1,mq1c,wquark)


! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    vertAs = vertAs + CF*NC*2d0*MstH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A2Lc,A2Rc,Gc,mq1,mq2c,wquark)

	  enddo


! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    vertAs = vertAs + CF*NC*2d0*MsuH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A3Lc,A3Rc,HLc,HRc,mq1,mq3c,wquark)

	  enddo

!  	    testvar = vertAs
!  	    print *,"vertAs",vertAs, UVdiv

	  NeuQ2qx_M2VertexH = NeuQ2qx_M2VertexH + vertAs
!          print *,"vertAs",vertAs, IRdiv
!          print *,"vertAs",vertAs

! -----------------------------------------------------------------------------------------------------



! ------------------------- Correction F-coupling (Higgs - Fermion - Fermion) -------------------------
!                                             s-channel

	  vertFs = 0d0

! --------- Gluon correction to Higgs-Fermion-Fermion vertex in S-channel ---------

	  call myCget(sqrtS**2,m4**2,m3**2,0d0,Mf(pstype,psgen)**2,m3**2)
	  call myBget(m4**2,Mf(pstype,psgen)**2,m3**2)

	  g0L = -GS
	  g0R = -GS
	  g1L = -GS
	  g1R = -GS
	  g2L = FL(xtype,pstype,psgen)
	  g2R = FR(xtype,pstype,psgen)

	  Mloop0 = 0d0
	  Mloop1 = Mf(pstype,psgen)
	  Mloop2 = m3

	  FsL = (g0R*g1L*g2L*(2*xb0 + (sman + m3**2 - m4**2 + 2*Mloop0**2)*xc0 + (3*sman + m3**2 - m4**2)*xc1 + (sman + 3*m3**2 - m4**2)*xc2) +
     &  g0R*g1R*g2R*Mloop1*m3*(- xc0 - xc2) + 2*g0R*g1L*g2R*Mloop1*Mloop2*xc0 - g0R*g1R*g2L*Mloop2*m3*xc2)/(8.*Pi**2)

	  FsR = (g0L*g1R*g2R*(2*xb0 + (sman + m3**2 - m4**2)*xc0 + (3*sman + m3**2 - m4**2)*xc1 + (sman + 3*m3**2 - m4**2)*xc2) -
     &  g0L*g1L*g2L*Mloop1*m3*(xc0 + xc2) + 2*g0L*g1R*g2L*Mloop1*Mloop2*xc0 - g0L*g1L*g2R*Mloop2*m3*xc2)/(8.*Pi**2)

	  FvL = -(g0L*g1L*g2L*Mloop1*xc1 + g0L*g1L*g2R*Mloop2*(xc0 + xc1))/(8.*Pi**2)

	  FvR = -(g0R*g1R*g2R*Mloop1*xc1 + g0R*g1R*g2L*Mloop2*(xc0 + xc1))/(8.*Pi**2)


! --------- Gluino correction to Higgs-Fermion-Fermion vertex in S-channel ---------

	  do k=1,2
	      do l=1,2

	  call myCget(sqrtS**2,m4**2,m3**2,MGl**2,MSf(k,pstype,psgen)**2,MSf(l,pstype,psgen)**2)

	  g0L = -sqrt2*GS*USf(k,1,pstype,psgen)
	  g0R = sqrt2*GS*USf(k,2,pstype,psgen)
	  g1L = sqrt2*GS*USf(l,2,pstype,psgen)
	  g1R = -sqrt2*GS*USf(l,1,pstype,psgen)
	  g2 = G(xtype,k,l,pstype,psgen)

	  Mloop0 = MGl
	  Mloop1 = MSf(k,pstype,psgen)
	  Mloop2 = MSf(l,pstype,psgen)

	  FsL = FsL - (g0L*g1L*g2*Mloop0*xc0 - g0L*g1R*g2*m3*xc2)/(16.*Pi**2)

	  FsR = FsR - (g0R*g1R*g2*Mloop0*xc0 - g0R*g1L*g2*m3*xc2)/(16.*Pi**2)

	  FvL = FvL + (g0R*g1L*g2*xc1)/(16.*Pi**2)

	  FvR = FvR + (g0L*g1R*g2*xc1)/(16.*Pi**2)

	        enddo
	    enddo


! Neutralino-Fermion-Sfermion
	  A1vL = 0d0
	  A1vR = 0d0
	  A1sL = AL(chimeigen,sfmeigen,sftype,sfgen)
	  A1sR = AR(chimeigen,sfmeigen,sftype,sfgen)

! propagator
	  PvL = 1d0
	  PvR = 1d0
	  PsL = Mf(pstype,psgen)
	  PsR = Mf(pstype,psgen)

! Higgs-Fermion-Fermion - corrected vertex
!     FvL=0d0
!     FvR=0d0
!     FsL=FL(xtype,pstype,psgen)
!     FsR=FR(xtype,pstype,psgen)

! mass of propagator & width
	  mq1 = Mf(pstype,psgen)
	  wquark = 0d0


! --- tree-level - s-channel ---
	    A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
	    A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
	    FLc = dconjg(FL(xtype,pstype,psgen))
	    FRc = dconjg(FR(xtype,pstype,psgen))

	    mq1c = Mf(pstype,psgen)

	    vertFs = vertFs + CF*NC*2d0*MssH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A1Lc,A1Rc,FLc,FRc,mq1,mq1c,wquark)


! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    vertFs = vertFs + CF*NC*2d0*MstH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A2Lc,A2Rc,Gc,mq1,mq2c,wquark)

	  enddo


! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    vertFs = vertFs + CF*NC*2d0*MsuH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A3Lc,A3Rc,HLc,HRc,mq1,mq3c,wquark)

	  enddo
!  	    testvar1 = vertFs
!  	    print *,"vertFs",vertFs, UVdiv

	  NeuQ2qx_M2VertexH = NeuQ2qx_M2VertexH + vertFs
!          print *,"vertFs",vertFs, IRdiv
!          print *,"vertFs",vertFs

! -----------------------------------------------------------------------------------------------------------



! ------------------------- Correction A-coupling (Neutralino - Sfermion - Fermion) -------------------------
!                                                 t-channel

	  vertAt = 0d0

	  do sumT = 0,sumTend

! --------- Gluon correction to Neutralino-Sfermion-Fermion vertex in T-channel ---------

	  call myCget(m1**2,tman,m3**2,m3**2,MSf(ptmeigen+sumT,pttype,ptgen)**2,0d0)
	  call myBget(tman,MSf(ptmeigen+sumT,pttype,ptgen)**2,0d0)

	  g0L = -GS
	  g0R = -GS
	  g1L = AL(chimeigen,ptmeigen+sumT,pttype,ptgen)
	  g1R = AR(chimeigen,ptmeigen+sumT,pttype,ptgen)
	  g2L = -GS
	  g2R = -GS

	  Mloop0 = m3
	  Mloop1 = MSf(ptmeigen+sumT,pttype,ptgen)
	  Mloop2 = 0d0

	  A2sL = (g0L*g1L*g2L*(xb0 + Mloop0**2*xc0 - m3**2*xc2 + 2*m1**2*xc1 + 2*(m1**2 + m3**2 - tman)*xc2) -
     &  g0R*g1R*g2R*m1*m3*(xc1 + 2*xc2) - g0L*g1R*g2L*Mloop0*m1*(2*xc0 + xc1) + g0R*g1L*g2R*Mloop0*m3*(xc0 - xc2))/(16.*Pi**2)

	  A2sR = (g0R*g1R*g2R*(xb0 + Mloop0**2*xc0 - m3**2*xc2 + 2*m1**2*xc1 + 2*(m1**2 + m3**2 - tman)*xc2) -
     &  g0L*g1L*g2L*m1*m3*(xc1 + 2*xc2) - g0R*g1L*g2R*Mloop0*m1*(2*xc0 + xc1) + g0L*g1R*g2L*Mloop0*m3*(xc0 - xc2))/(16.*Pi**2)


! --------- Gluino correction to Neutralino-Sfermion-Fermion vertex in T-channel ---------

	  do k=1,2

	  call myCget(m1**2,tman,m3**2,MSf(k,pttype,ptgen)**2,Mf(pttype,ptgen)**2,MGl**2)
	  call myBget(tman,Mf(pttype,ptgen)**2,MGl**2)

	  g0L =  sqrt2*GS*USf(ptmeigen+sumT,2,pttype,ptgen)
	  g0R = -sqrt2*GS*USf(ptmeigen+sumT,1,pttype,ptgen)
	  g1L = dconjg(AR(chimeigen,k,pttype,ptgen))
	  g1R = dconjg(AL(chimeigen,k,pttype,ptgen))
	  g2L =  sqrt2*GS*USf(k,2,pttype,ptgen)
	  g2R = -sqrt2*GS*USf(k,1,pttype,ptgen)

	  Mloop0 = MSf(k,pttype,ptgen)
	  Mloop1 = Mf(pttype,ptgen)
	  Mloop2 = MGl


	  A2sL = A2sL - (g0L*g1L*g2L*Mloop1*Mloop2*xc0 + g0L*g1R*g2L*Mloop2*m1*(xc0 + xc1) + g0L*g1L*g2R*Mloop1*m3*(xc0 + xc2) +
     & g0L*g1R*g2R*m1*m3*(xc0 + xc1 + xc2) + g0R*g1L*g2R*Mloop2*m3*xc2 + g0R*g1R*g2L*Mloop1*m1*xc1 +
     & g0R*g1L*g2L*(xb0 + Mloop0**2*xc0 + m1**2*xc1 + m3**2*xc2))/(16.*Pi**2)

	  A2sR = A2sR - (g0R*g1R*g2R*Mloop1*Mloop2*xc0 + g0R*g1L*g2R*Mloop2*m1*(xc0 + xc1) + g0R*g1R*g2L*Mloop1*m3*(xc0 + xc2) +
     & g0R*g1L*g2L*m1*m3*(xc0 + xc1 + xc2) + g0L*g1R*g2L*Mloop2*m3*xc2 + g0L*g1L*g2R*Mloop1*m1*xc1 +
     & g0L*g1R*g2R*(xb0 + Mloop0**2*xc0 + m1**2*xc1 + m3**2*xc2))/(16.*Pi**2)

	  enddo



! Neutralino-Fermion-Sfermion - corrected Vertex
!	  A2sL = AL(chimeigen,ptmeigen+sumT,pttype,ptgen)
!	  A2sR = AR(chimeigen,ptmeigen+sumT,pttype,ptgen)

! propagator
	  P = 1d0

! Higgs-Sfermion-Sfermion
          Gscalar=G(xtype,sfmeigen,ptmeigen+sumT,sftype,sfgen)

! mass of propagator & width
	  mq2 = MSf(ptmeigen+sumT,pttype,ptgen)


! --- tree-level - s-channel ---
	    A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
	    A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
	    FLc = dconjg(FL(xtype,pstype,psgen))
	    FRc = dconjg(FR(xtype,pstype,psgen))

	    mq1c = Mf(pstype,psgen)
	    wquark = 0d0

	    vertAt = vertAt + CF*NC*2d0*MtsH(A2sL,A2sR,P,Gscalar,A1Lc,A1Rc,FLc,FRc,mq2,mq1c,wquark)


! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    vertAt = vertAt + CF*NC*2d0*MttH(A2sL,A2sR,P,Gscalar,A2Lc,A2Rc,Gc,mq2,mq2c)

	  enddo


! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    vertAt = vertAt + CF*NC*2d0*MtuH(A2sL,A2sR,P,Gscalar,A3Lc,A3Rc,HLc,HRc,mq2,mq3c)

	  enddo

	  enddo

!  	    testvar2 = vertAt
!  	    print *,"vertAt",vertAt, UVdiv

	  NeuQ2qx_M2VertexH = NeuQ2qx_M2VertexH + vertAt
!          print *,"vertAt",vertAt, IRdiv
!          print *,"vertAt",vertAt

! -------------------------------------------------------------------------------------------------------



! ------------------------- Correction G-coupling (Higgs - Sfermion - Sfermion) -------------------------
!                                              t-channel

	  vertGt = 0d0
	  Gscalar = 0d0
! because of 4-Squark H3-case...
	  do sumT = 0,1


! --------- Gluon correction to Sfermion-Sfermion-Higgs vertex in T-channel ---------

	  call myCget(m4**2,tman,m2**2,m2**2,MSf(1+sumT,pttype,ptgen)**2,0d0)
	  call myBget(tman,MSf(1+sumT,pttype,ptgen)**2,0d0)

	  g0 = G(xtype,sfmeigen,1+sumT,sftype,sfgen)
	  g1 = -GS
	  g2 = -GS

	  Mloop0 = m2
	  Mloop1 = MSf(1+sumT,pttype,ptgen)
	  Mloop2 = 0d0

	  Gscalar =  g0*g1*g2*(xb0 + (Mloop0**2 + tman - m4**2)*xc0 + (tman + m4**2 - m2**2)*xc1 + (m4**2 - m2**2 - tman)*xc2)/(16.*Pi**2)


! --------- Gluino correction to Sfermion-Sfermion-Higgs vertex in T-channel ---------

	  call myCget(m4**2,tman,m2**2,Mf(sftype,sfgen)**2,Mf(pttype,ptgen)**2,MGl**2)
	  call myBget(tman,Mf(pttype,ptgen)**2,MGl**2)

	  g0L = FL(xtype,sftype,sfgen)
	  g0R = FR(xtype,sftype,sfgen)
	  g1L = -sqrt2*GS*USf(1+sumT,1,pttype,ptgen)
	  g1R =  sqrt2*GS*USf(1+sumT,2,pttype,ptgen)
	  g2L = sqrt2*GS*USf(sfmeigen,2,sftype,sfgen)
	  g2R = -sqrt2*GS*USf(sfmeigen,1,sftype,sfgen)

	  Mloop0 = Mf(sftype,sfgen)
	  Mloop1 = Mf(pttype,ptgen)
	  Mloop2 = MGl

	  Gscalar = Gscalar + ((g0L*g1L*g2L + g0R*g1R*g2R)*2*Mloop0*Mloop1*Mloop2*xc0 +
     &  (g0R*g1R*g2L + g0L*g1L*g2R)*Mloop1*(2*(xb0 + Mloop0**2*xc0 + m2**2*xc2) + (m2**2 + m4**2 - tman)*xc1) +
     &  (g0L*g1R*g2R + g0R*g1L*g2L)*Mloop2*(2*(xb0 + Mloop0**2*xc0 + m4**2*xc1) + (m2**2 + m4**2 - tman)*xc2) +
     &  (g0L*g1R*g2L + g0R*g1L*g2R)*Mloop0*(2*(xb0 + Mloop0**2*xc0 + m4**2*xc1 + m2**2*xc2) + (m2**2 + m4**2 - tman)*(xc0 + xc1 + xc2)))/(16.*Pi**2)


! --------- 4-Squark correction to Sfermion-Sfermion-Higgs vertex in T-channel ---------

	    do l=1,2
	      do i=1,2
	  call myBget(m4**2,MSf(l,pttype,ptgen)**2,MSf(i,sftype,sfgen)**2)

	  g0 = S(sfmeigen,i,l,1+sumT,sftype,pttype,sfgen)
	  g1 = G(xtype,i,l,sftype,sfgen)

	  Gscalar = Gscalar + g0*g1*xb0/(16.*Pi**2)



	      enddo
	    enddo


! Neutralino-Fermion-Sfermion
	  A2sL = AL(chimeigen,1+sumT,pttype,ptgen)
	  A2sR = AR(chimeigen,1+sumT,pttype,ptgen)

! propagator
	  P = 1d0

! Higgs-Sfermion-Sfermion - corrected Vertex
!      Gscalar=G(xtype,sfmeigen,ptmeigen+sumT,sftype,sfgen)

! mass of propagator & width
	  mq2 = MSf(1+sumT,pttype,ptgen)


! --- tree-level - s-channel ---
	    A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
	    A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
	    FLc = dconjg(FL(xtype,pstype,psgen))
	    FRc = dconjg(FR(xtype,pstype,psgen))

	    mq1c = Mf(pstype,psgen)
	    wquark = 0d0
	    

	    vertGt = vertGt + CF*NC*2d0*MtsH(A2sL,A2sR,P,Gscalar,A1Lc,A1Rc,FLc,FRc,mq2,mq1c,wquark)


! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    vertGt = vertGt + CF*NC*2d0*MttH(A2sL,A2sR,P,Gscalar,A2Lc,A2Rc,Gc,mq2,mq2c)

	  enddo


! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    vertGt = vertGt + CF*NC*2d0*MtuH(A2sL,A2sR,P,Gscalar,A3Lc,A3Rc,HLc,HRc,mq2,mq3c)

	  enddo

	  enddo

!  	    testvar3 = vertGt
!  	    print *,"vertGt",vertGt, UVdiv

	  NeuQ2qx_M2VertexH = NeuQ2qx_M2VertexH + vertGt
!          print *,"vertGt",vertGt, IRdiv
!          print *,"vertGt",vertGt

! -----------------------------------------------------------------------------------------------------------



! ------------------------- Correction A-coupling (Neutralino - Sfermion - Fermion) -------------------------
!                                              u-channel

	  vertAu = 0d0

	  do sumU = 0,sumUend

! --------- Gluon correction to Neutralino-Sfermion-Fermion vertex in U-channel ---------

	  call myCget(uman,m2**2,m3**2,m3**2,MSf(sfmeigen,sftype,sfgen)**2,0d0)
	  call myBget(m2**2,MSf(sfmeigen,sftype,sfgen)**2,0d0)

	  g0L = -GS
	  g0R = -GS
	  g1L = AL(pumeigen+sumU,sfmeigen,sftype,sfgen)
	  g1R = AR(pumeigen+sumU,sfmeigen,sftype,sfgen)
	  g2L = -GS
	  g2R = -GS

	  Mloop0 = m3
	  Mloop1 = MSf(sfmeigen,sftype,sfgen)
	  Mloop2 = 0d0

	  A3sL = (g0L*g1L*g2L*(xb0 + Mloop0**2*xc0 + 2*uman*xc1 + (2*uman - 2*m2**2 + m3**2)*xc2) +
     & g0R*g1L*g2R*m3*Mloop0*(xc0 - xc2))/(16.*Pi**2)

	  A3sR = (g0R*g1R*g2R*(xb0 + Mloop0**2*xc0 + 2*uman*xc1 + (2*uman - 2*m2**2 + m3**2)*xc2) +
     & g0L*g1R*g2L*m3*Mloop0*(xc0 - xc2))/(16.*Pi**2)

	  A3vL = -(g0L*g1R*g2L*Mloop0*(2*xc0 + xc1) + g0R*g1R*g2R*m3*(xc1 + 2*xc2))/(16.*Pi**2)

	  A3vR = -(g0R*g1L*g2R*Mloop0*(2*xc0 + xc1) + g0L*g1L*g2L*m3*(xc1 + 2*xc2))/(16.*Pi**2)


! --------- Gluino correction to Neutralino-Sfermion-Fermion vertex in U-channel ---------

	  do k=1,2

! ftype, fgen are type & gen of the outgoing quark

	  call myCget(uman,m2**2,m3**2,MSf(k,ftype,fgen)**2,Mf(sftype,sfgen)**2,MGl**2)
	  call myBget(m2**2,Mf(sftype,sfgen)**2,MGl**2)

	  g0L =  sqrt2*GS*USf(sfmeigen,2,sftype,sfgen)
	  g0R = -sqrt2*GS*USf(sfmeigen,1,sftype,sfgen)
	  g1L = dconjg(AR(pumeigen+sumU,k,ftype,fgen))
	  g1R = dconjg(AL(pumeigen+sumU,k,ftype,fgen))
	  g2L =  sqrt2*GS*USf(k,2,ftype,fgen)
	  g2R = -sqrt2*GS*USf(k,1,ftype,fgen)

	  Mloop0 = MSf(k,ftype,fgen)
	  Mloop1 = Mf(sftype,sfgen)
	  Mloop2 = MGl


	  A3sL = A3sL - (g0L*g1L*g2L*Mloop1*Mloop2*xc0 + g0L*g1L*g2R*Mloop1*m3*(xc0 + xc2) + g0R*g1L*g2R*Mloop2*m3*xc2 +
     &  g0R*g1L*g2L*(xb0 + Mloop0**2*xc0 + uman*xc1 + m3**2*xc2))/(16.*Pi**2)

	  A3sR = A3sR - (g0R*g1R*g2R*Mloop1*Mloop2*xc0 + g0R*g1R*g2L*Mloop1*m3*(xc0 + xc2) + g0L*g1R*g2L*Mloop2*m3*xc2 +
     &  g0L*g1R*g2R*(xb0 + Mloop0**2*xc0 + uman*xc1 + m3**2*xc2))/(16.*Pi**2)

	  A3vL = A3vL - (g0R*g1R*g2L*Mloop1*xc1 + g0L*g1R*g2L*Mloop2*(xc0 + xc1) + g0L*g1R*g2R*m3*(xc0 + xc1 + xc2))/(16.*Pi**2)

	  A3vR = A3vR - (g0L*g1L*g2R*Mloop1*xc1 + g0R*g1L*g2R*Mloop2*(xc0 + xc1) + g0R*g1L*g2L*m3*(xc0 + xc1 + xc2))/(16.*Pi**2)

	  enddo


! Neutralino-Fermion-Sfermion - corrected Vertex
!	   A3sL = AL(pumeigen+sumU,sfmeigen,sftype,sfgen)
!	   A3sR = AR(pumeigen+sumU,sfmeigen,sftype,sfgen)
!	   A3vL = 0d0
!	   A3vR = 0d0

! Higgs-Neutralino/Chargino-Sfermion
	  HtreeL = HL(chimeigen,pumeigen+sumU,xtype)
	  HtreeR = HR(chimeigen,pumeigen+sumU,xtype)

! mass of propagator & width
	   mq3 =  MNeuCha(pumeigen+sumU)


! --- tree-level - s-channel ---
	    A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
	    A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
	    FLc = dconjg(FL(xtype,pstype,psgen))
	    FRc = dconjg(FR(xtype,pstype,psgen))

	    mq1c = Mf(pstype,psgen)
	    wquark = 0d0

	    vertAu = vertAu + CF*NC*2d0*MusH(A3vL,A3vR,A3sL,A3sR,HtreeL,HtreeR,A1Lc,A1Rc,FLc,FRc,mq3,mq1c,wquark)


! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    vertAu = vertAu + CF*NC*2d0*MutH(A3vL,A3vR,A3sL,A3sR,HtreeL,HtreeR,A2Lc,A2Rc,Gc,mq3,mq2c)

	  enddo


! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    vertAu = vertAu + CF*NC*2d0*MuuH(A3vL,A3vR,A3sL,A3sR,HtreeL,HtreeR,A3Lc,A3Rc,HLc,HRc,mq3,mq3c)

	  enddo

	  enddo

!  	    testvar4 = vertAu
!  	    print *,"vertAu",vertAu, UVdiv

	  NeuQ2qx_M2VertexH = NeuQ2qx_M2VertexH + vertAu
!          print *,"vertAu",vertAu, IRdiv
!          print *,"vertAu",vertAu

! -----------------------------------------------------------------------------------------------------------

	   end


