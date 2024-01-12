! ==========================================================================
	  double precision function NeuQ2qx_M2PropagatorH()

	  implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
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

	  double complex NeuQ2qx_PropLoopSfermion

	  integer k,r
	  integer sumT1, sumT2, sumTc, sumUc
	  integer sumT1end, sumT2end, sumTcend, sumUcend

	  r = 0


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


	  NeuQ2qx_M2PropagatorH = 0d0
	  props = 0d0
	  propt = 0d0

! ------------------------- Correction in s-channel propagator (Quark) -------------------------


! --------- Gluon contribution ---------

	  call myBget(sman,0d0,Mf(pstype,psgen)**2)

	  g0 = -GS
	  g1 = -GS

	  Mloop0 = 0d0
	  Mloop1 = Mf(pstype,psgen)

	  ZsL = g0*g1*(Mloop1*(1d0*r - 2d0*xb0))/(8.*Pi**2)
	  ZsR = ZsL
	  ZvL = g0*g1*(2d0*xb0 + 2d0*xb1 - 1d0*r)/(16.*Pi**2)
	  ZvR = ZvL

! --------- Gluino contribution ---------

	  do k=1,2

	  call myBget(sman,MGl**2,MSf(k,pstype,psgen)**2)


	  g0L = -sqrt2*GS*USf(k,1,pstype,psgen)
	  g0R = sqrt2*GS*USf(k,2,pstype,psgen)
	  g1L = sqrt2*GS*USf(k,2,pstype,psgen)
	  g1R = -sqrt2*GS*USf(k,1,pstype,psgen)

	  Mloop0 = MGl
	  Mloop1 = MSf(k,pstype,psgen)


	  ZsL = ZsL + (Mloop0*g0L*g1L*xb0)/(16.*Pi**2)
	  ZsR = ZsR + (Mloop0*g0R*g1R*xb0)/(16.*Pi**2)
	  ZvL = ZvL - (g0R*g1L*xb1)/(16.*Pi**2)
	  ZvR = ZvR - (g0L*g1R*xb1)/(16.*Pi**2)


	  enddo


	  PsL = (ZvL + ZvR)*Mf(pstype,psgen)*sman + ZsL*Mf(pstype,psgen)**2 + ZsR*sman

	  PsR = (ZvL + ZvR)*Mf(pstype,psgen)*sman + ZsR*Mf(pstype,psgen)**2 + ZsL*sman

	  PvL = ZvL*Mf(pstype,psgen)**2 + ZvR*sman + (ZsL + ZsR)*Mf(pstype,psgen)

	  PvR = ZvR*Mf(pstype,psgen)**2 + ZvL*sman + (ZsL + ZsR)*Mf(pstype,psgen)



! Neutralino-Fermion-Sfermion
	  A1vL = 0d0
	  A1vR = 0d0
	  A1sL = AL(chimeigen,sfmeigen,sftype,sfgen)
	  A1sR = AR(chimeigen,sfmeigen,sftype,sfgen)

!! propagator -  corrected
!	  PvL = 1d0
!	  PvR = 1d0
!	  PsL = Mf(pstype,psgen)
!	  PsR = Mf(pstype,psgen)

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

	    props = props - CF*NC*2d0*MssH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A1Lc,A1Rc,FLc,FRc,mq1,mq1c,wquark)/(sman - mq1**2 +
     &         (cI*mq1*wquark))


! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    props = props - CF*NC*2d0*MstH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A2Lc,A2Rc,Gc,mq1,mq2c,wquark)/(sman -
     &         mq1c**2 + (cI*mq1c*wquark))


	  enddo



! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)


	    props = props - CF*NC*2d0*MsuH(A1vL,A1vR,A1sL,A1sR,PvL,PvR,PsL,PsR,FvL,FvR,FsL,FsR,A3Lc,A3Rc,HLc,HRc,mq1,mq3c,wquark)/(sman -
     & mq1c**2 + (cI*mq1c*wquark))


	  enddo



	  NeuQ2qx_M2PropagatorH = NeuQ2qx_M2PropagatorH + props



! ------------------------- Correction in t-channel propagator (Squark) -------------------------

	  do sumT1 = 0,sumT1end
	    do sumT2 = 0,sumT2end


! propagator - the corrected one
	  P = NeuQ2qx_PropLoopSfermion(tman,ptmeigen+sumT2,ptmeigen+sumT1,pttype,ptgen)

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


	    propt = propt - CF*NC*2d0*MtsH(A2sL,A2sR,P,Gscalar,A1Lc,A1Rc,FLc,FRc,mq2,mq1c,wquark)/(tman - MSf(ptmeigen+sumT2,pttype,ptgen)**2)



! --- tree-level - t-channel ---

	  do sumTc = 0,sumTcend

	    A2Lc = dconjg(AL(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    A2Rc = dconjg(AR(chimeigen,ptmeigen+sumTc,pttype,ptgen))
	    Gc = dconjg(G(xtype,sfmeigen,ptmeigen+sumTc,sftype,sfgen))

	    mq2c = MSf(ptmeigen+sumTc,pttype,ptgen)

	    propt = propt - CF*NC*2d0*MttH(A2sL,A2sR,P,Gscalar,A2Lc,A2Rc,Gc,mq2,mq2c)/(tman - MSf(ptmeigen+sumT2,pttype,ptgen)**2)

	  enddo


! --- tree-level - u-channel ---

	  do sumUc = 0,sumUcend

	    A3Lc = dconjg(AL(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    A3Rc = dconjg(AR(pumeigen+sumUc,sfmeigen,sftype,sfgen))
	    HLc = dconjg(HL(chimeigen,pumeigen+sumUc,xtype))
	    HRc = dconjg(HR(chimeigen,pumeigen+sumUc,xtype))

	    mq3c = MNeuCha(pumeigen+sumUc)

	    propt = propt - CF*NC*2d0*MtuH(A2sL,A2sR,P,Gscalar,A3Lc,A3Rc,HLc,HRc,mq2,mq3c)/(tman - MSf(ptmeigen+sumT2,pttype,ptgen)**2)

	  enddo

	  enddo
	  enddo

	  NeuQ2qx_M2PropagatorH = NeuQ2qx_M2PropagatorH + propt

	  end


! ------------------------- Function for propagator in t-channel (Sfermion) -------------------------

	  double complex function NeuQ2qx_PropLoopSfermion(p2,ineigen,outeigen,type,fam)

	  implicit none

#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	  double precision p2
	  integer ineigen,outeigen,type,fam,m
	  double complex g0,g0L,g0R,g1,g1L,g1R,g2,g2L,g2R
	  double precision Mloop0,Mloop1,Mloop2
	  double precision A01,Krondelta


	  NeuQ2qx_PropLoopSfermion = 0d0

! --------- Gluon contribution ----------

	  call myBget(p2,0d0,MSf(ineigen,type,fam)**2)
	  call myAget(MSf(ineigen,type,fam)**2)

	  g0 = -GS*Krondelta(ineigen,outeigen)
	  g1 = -GS*Krondelta(ineigen,outeigen)


	  Mloop0 = 0d0
	  Mloop1 = MSf(ineigen,type,fam)

	  NeuQ2qx_PropLoopSfermion = -(g0*g1*(2d0*(p2 + Mloop1**2)*xb0 - xa0))/(16.*Pi**2)


! --------- Gluino contribution ---------

	  call myBget(p2,MGl**2,Mf(type,fam)**2)

	  g0L = sqrt2*GS*USf(outeigen,2,type,fam)
	  g0R = -sqrt2*GS*USf(outeigen,1,type,fam)
	  g1L = -sqrt2*GS*USf(ineigen,1,type,fam)
	  g1R = sqrt2*GS*USf(ineigen,2,type,fam)

	  Mloop0 = MGl
	  Mloop1 = Mf(type,fam)


	  NeuQ2qx_PropLoopSfermion = NeuQ2qx_PropLoopSfermion - ((g0L*g1R + g0R*g1L)*(A01(Mloop1) + A01(Mloop0) + (Mloop0**2 + Mloop1**2 -p2)*xb0)
     &  + 2d0*(g0R*g1R + g0L*g1L)*Mloop0*Mloop1*xb0)/(16.*Pi**2)
     

! --------- Squark contribution ---------

	  do m=1,2

	      g0 = S(ineigen,m,m,outeigen,type,type,fam)

	    Mloop0 = MSf(m,type,fam)

	    NeuQ2qx_PropLoopSfermion = NeuQ2qx_PropLoopSfermion - (g0*A01(Mloop0))/(16.*Pi**2)

	  enddo


	  end
