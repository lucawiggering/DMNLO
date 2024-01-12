! =======================================================================
	  double precision function NeuQ2qx_M2VertexCounterH()

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

	  integer k,l

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


	  NeuQ2qx_M2VertexCounterH = 0d0

! ------------------------- Correction A-coupling (Neutralino - Sfermion - Fermion) -------------------------
!                                                 s-channel

	  vertAs = 0d0


	  A1sL = dAL(chimeigen,sfmeigen,sftype,sfgen)
	  A1sR = dAR(chimeigen,sfmeigen,sftype,sfgen)
	  A1vL = 0d0
	  A1vR = 0d0

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

	    vertAs = vertAs + NC*2d0*MssH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A1Lc,A1Rc,FLc,FRc,mq1,mq1c,wquark)


! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    vertAs = vertAs + NC*2d0*MstH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A2Lc,A2Rc,Gc,mq1,mq2c,wquark)

	  enddo



! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    vertAs = vertAs + NC*2d0*MsuH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A3Lc,A3Rc,HLc,HRc,mq1,mq3c,wquark)

	  enddo


	  NeuQ2qx_M2VertexCounterH = NeuQ2qx_M2VertexCounterH + vertAs
	  
! -----------------------------------------------------------------------------------------------------



! ------------------------- Correction F-coupling (Higgs - Fermion - Fermion) -------------------------
!                                             s-channel

	  vertFs = 0d0


	  FsL = dFL(xtype,pstype,psgen)
	  FsR = dFR(xtype,pstype,psgen)
	  FvL = 0d0
	  FvR = 0d0


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

	    vertFs = vertFs + NC*2d0*MssH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A1Lc,A1Rc,FLc,FRc,mq1,mq1c,wquark)

! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    vertFs = vertFs + NC*2d0*MstH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A2Lc,A2Rc,Gc,mq1,mq2c,wquark)

	  enddo


! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    vertFs = vertFs + NC*2d0*MsuH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A3Lc,A3Rc,HLc,HRc,mq1,mq3c,wquark)

	  enddo


	  NeuQ2qx_M2VertexCounterH = NeuQ2qx_M2VertexCounterH + vertFs

! -----------------------------------------------------------------------------------------------------------



! ------------------------- Correction A-coupling (Neutralino - Sfermion - Fermion) -------------------------
!                                                 t-channel

	  vertAt = 0d0

	  do sumT = 0,sumTend


	  A2sL = dAL(chimeigen,ptmeigen+sumT,pttype,ptgen)
	  A2sR = dAR(chimeigen,ptmeigen+sumT,pttype,ptgen)


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

	    vertAt = vertAt + NC*2d0*MtsH(A2sL,A2sR,P,Gscalar,A1Lc,A1Rc,FLc,FRc,mq2,mq1c,wquark)


! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    vertAt = vertAt + NC*2d0*MttH(A2sL,A2sR,P,Gscalar,A2Lc,A2Rc,Gc,mq2,mq2c)

	  enddo


! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    vertAt = vertAt + NC*2d0*MtuH(A2sL,A2sR,P,Gscalar,A3Lc,A3Rc,HLc,HRc,mq2,mq3c)

	  enddo

	  enddo



	  NeuQ2qx_M2VertexCounterH = NeuQ2qx_M2VertexCounterH + vertAt
          
          
! -------------------------------------------------------------------------------------------------------



! ------------------------- Correction G-coupling (Higgs - Sfermion - Sfermion) -------------------------
!                                              t-channel

	  vertGt = 0d0
          Gscalar = 0d0

	  do sumT = 0,sumTend


	  Gscalar = dG(xtype,sfmeigen,ptmeigen+sumT,sftype,sfgen)


! Neutralino-Fermion-Sfermion
	  A2sL = AL(chimeigen,ptmeigen+sumT,pttype,ptgen)
	  A2sR = AR(chimeigen,ptmeigen+sumT,pttype,ptgen)

! propagator
	  P = 1d0

! Higgs-Sfermion-Sfermion - corrected Vertex
!      Gscalar=G(xtype,sfmeigen,ptmeigen+sumT,sftype,sfgen)

! mass of propagator & width
	  mq2 = MSf(ptmeigen+sumT,pttype,ptgen)


! --- tree-level - s-channel ---
	    A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
	    A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
	    FLc = dconjg(FL(xtype,pstype,psgen))
	    FRc = dconjg(FR(xtype,pstype,psgen))

	    mq1c = Mf(pstype,psgen)
	    wquark = 0d0

	    vertGt = vertGt + NC*2d0*MtsH(A2sL,A2sR,P,Gscalar,A1Lc,A1Rc,FLc,FRc,mq2,mq1c,wquark)


! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    vertGt = vertGt + NC*2d0*MttH(A2sL,A2sR,P,Gscalar,A2Lc,A2Rc,Gc,mq2,mq2c)

	  enddo


! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    vertGt = vertGt + NC*2d0*MtuH(A2sL,A2sR,P,Gscalar,A3Lc,A3Rc,HLc,HRc,mq2,mq3c)

	  enddo

	  enddo

	  NeuQ2qx_M2VertexCounterH = NeuQ2qx_M2VertexCounterH + vertGt  
          
! -----------------------------------------------------------------------------------------------------------



! ------------------------- Correction A-coupling (Neutralino - Sfermion - Fermion) -------------------------
!                                              u-channel

	  vertAu = 0d0

	  do sumU = 0,sumUend


	  A3sL = dAL(pumeigen+sumU,sfmeigen,sftype,sfgen)
	  A3sR = dAR(pumeigen+sumU,sfmeigen,sftype,sfgen)
	  A3vL = 0d0
	  A3vR = 0d0

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

	    vertAu = vertAu + NC*2d0*MusH(A3vL,A3vR,A3sL,A3sR,HtreeL,HtreeR,A1Lc,A1Rc,FLc,FRc,mq3,mq1c,wquark)


! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    vertAu = vertAu + NC*2d0*MutH(A3vL,A3vR,A3sL,A3sR,HtreeL,HtreeR,A2Lc,A2Rc,Gc,mq3,mq2c)

	  enddo


! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    vertAu = vertAu + NC*2d0*MuuH(A3vL,A3vR,A3sL,A3sR,HtreeL,HtreeR,A3Lc,A3Rc,HLc,HRc,mq3,mq3c)

	  enddo

	  enddo


	  NeuQ2qx_M2VertexCounterH = NeuQ2qx_M2VertexCounterH + vertAu
          

! -----------------------------------------------------------------------------------------------------------



	  end
