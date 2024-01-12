! ========================================================================
	  double precision function NeuQ2qx_M2PropagatorCounterH()

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

	  double complex PvL,PvR,PsL,PsR,P,ZvL,ZvR,ZsL,ZsR

	  double precision mq1,mq1c,mq2,mq2c,mq3,mq3c,wquark

	  double complex g0,g0L,g0R,g1,g1L,g1R
	  double precision Mloop0,Mloop1

	  double precision props, propt
	  double precision MssH, MttH, MuuH, MstH, MtsH, MsuH, MtuH, MusH, MutH


	  integer sumT1, sumT2, sumTc, sumUc
	  integer sumT1end, sumT2end, sumTcend, sumUcend


	  sumT1 = 0
	  sumT2 = 0
	  sumTc = 0
	  sumUc = 0

	  if ((xtype.eq.5).or.(xtype.eq.6)) then

	  sumT1end = 1
	  sumT2end = 1
	  sumTcend = 1
	  sumUcend = 3

	  elseif (xtype.eq.7) then

	  sumT1end = 0
	  sumT2end = 0
	  sumTcend = 0
	  sumUcend = 3

	  elseif (xtype.eq.8) then

	  sumT1end = 1
	  sumT2end = 1
	  sumTcend = 1
	  sumUcend = 1

	  endif


	  NeuQ2qx_M2PropagatorCounterH = 0d0
	  props = 0d0
	  propt = 0d0


! ------------------------- Counterterms t-channel propagator  -------------------------

	  P = 0d0

          call NeuQ2qx_CalcPropSquark(tman,pttype,ptgen)

	  do sumT1 = 0,sumT1end
	    do sumT2 = 0,sumT2end
! 	  tman=sman
! 	  uman=sman

! propagator - the counterterm
	  P = Prop(ptmeigen+sumT2,ptmeigen+sumT1,pttype,ptgen)

! Neutralino-Fermion-Sfermion
	  A2sL = AL(chimeigen,ptmeigen+sumT1,pttype,ptgen)
	  A2sR = AR(chimeigen,ptmeigen+sumT1,pttype,ptgen)

! Higgs-Sfermion-Sfermion
          Gscalar=G(xtype,sfmeigen,ptmeigen+sumT2,sftype,sfgen)

! mass of propagator & width
	  mq2 = MSf(ptmeigen+sumT1,pttype,ptgen)


! --- tree-level - s-channel ---
	    A1Lc = dconjg(AL(chimeigen,sfmeigen,sftype,sfgen))
	    A1Rc = dconjg(AR(chimeigen,sfmeigen,sftype,sfgen))
	    FLc = dconjg(FL(xtype,pstype,psgen))
	    FRc = dconjg(FR(xtype,pstype,psgen))

	    mq1c = Mf(pstype,psgen)
	    wquark = 0d0

	    propt = propt - NC*2d0*MtsH(A2sL,A2sR,P,Gscalar,A1Lc,A1Rc,FLc,FRc,mq2,mq1c,wquark)/(tman - MSf(ptmeigen+sumT2,pttype,ptgen)**2)



! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    propt = propt - NC*2d0*MttH(A2sL,A2sR,P,Gscalar,A2Lc,A2Rc,Gc,mq2,mq2c)/(tman - MSf(ptmeigen+sumT2,pttype,ptgen)**2)

	  enddo


! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    propt = propt - NC*2d0*MtuH(A2sL,A2sR,P,Gscalar,A3Lc,A3Rc,HLc,HRc,mq2,mq3c)/(tman - MSf(ptmeigen+sumT2,pttype,ptgen)**2)

	  enddo

	  enddo
	  enddo

	  NeuQ2qx_M2PropagatorCounterH = NeuQ2qx_M2PropagatorCounterH + propt



! ------------------------- Correction in s-channel propagator  -------------------------

          call NeuQ2qx_CalcPropQuark(pstype,psgen)

	  ZsL = 0d0
	  ZsR = 0d0
	  ZvL = 0d0
	  ZvR = 0d0


	  ZsL = PropSL(pstype,psgen)
	  ZsR = PropSR(pstype,psgen)
	  ZvL = PropVL(pstype,psgen)
	  ZvR = PropVR(pstype,psgen)

	  PsL = 0d0
	  PsR = 0d0
	  PvL = 0d0
	  PvR = 0d0

	  PsL = (ZvL + ZvR)*Mf(pstype,psgen)*sman + ZsL*Mf(pstype,psgen)**2 + ZsR*sman

	  PsR = (ZvL + ZvR)*Mf(pstype,psgen)*sman + ZsR*Mf(pstype,psgen)**2 + ZsL*sman

	  PvL = ZvL*Mf(pstype,psgen)**2 + ZvR*sman + (ZsL + ZsR)*Mf(pstype,psgen)

	  PvR = ZvR*Mf(pstype,psgen)**2 + ZvL*sman + (ZsL + ZsR)*Mf(pstype,psgen)


! Neutralino-Fermion-Sfermion
	  A1vL = 0d0
	  A1vR = 0d0
	  A1sL = AL(chimeigen,sfmeigen,sftype,sfgen)
	  A1sR = AR(chimeigen,sfmeigen,sftype,sfgen)

! propagator -  corrected
! 	  PvL = 1d0
! 	  PvR = 1d0
! 	  PsL = Mf(pstype,psgen)
! 	  PsR = Mf(pstype,psgen)

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

	    props = props - NC*2d0*MssH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A1Lc,A1Rc,FLc,FRc,mq1,mq1c,wquark)/(sman - mq1**2 -
     &  (cI*mq1*wquark))



! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    props = props - NC*2d0*MstH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A2Lc,A2Rc,Gc,mq1,mq2c,wquark)/(sman -
     & mq1c**2 - (cI*mq1c*wquark))

	  enddo



! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    props = props - NC*2d0*MsuH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A3Lc,A3Rc,HLc,HRc,mq1,mq3c,wquark)/(sman -
     & mq1c**2 - (cI*mq1c*wquark))


	  enddo


	  NeuQ2qx_M2PropagatorCounterH = NeuQ2qx_M2PropagatorCounterH + props

	  end
