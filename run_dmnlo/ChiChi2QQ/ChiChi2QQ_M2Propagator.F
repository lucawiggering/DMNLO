	function PropLoop(p2,k,l,type,gg)

! type,gg are type,generation of sfermion in the tree-level propagator to be corrected

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"
	
	double precision p2
	integer type,gg,m,k,l

	double complex g1L,g1R,g2L,g2R,ML0,ML1,g1,g2
	double precision A01, Krondelta
	double complex B02
	
	double complex PropLoop
	
	PropLoop = 0d0

! *** Gluon Loop ***
! Note: gluon one-point diagram is zero due to A0(mg=0)=0

	g1 = -GS*Krondelta(k,l)
	g2 = -GS*Krondelta(k,l)
	
	ML0 = 0d0 
	ML1 = MSf(k,type,gg)

	PropLoop = PropLoop - 
     &		(g1*g2*(2*A01(ML0) - A01(ML1) + (-ML0**2 + 2*ML1**2 + 2*p2)*B02(p2,ML0,ML1)))/(16.*Pi**2)
		
! *** Gluino Loop ***
    	g1L = sqrt2*GS*USf(k,2,type,gg)
    	g1R = -sqrt2*GS*USf(k,1,type,gg)
    	g2L = -sqrt2*GS*USf(l,1,type,gg)
    	g2R = sqrt2*GS*USf(l,2,type,gg)
      
    	ML0 = Mf(type,gg) 
    	ML1 = MGl
	
	PropLoop = PropLoop - 
     &		(2*(g1L*g2L + g1R*g2R)*ML0*ML1*B02(p2,ML0,ML1) + 
     &     	(g1R*g2L + g1L*g2R)*(A01(ML0) + A01(ML1) + (ML0**2 + ML1**2 - p2)*B02(p2,ML0,ML1)))/(16.*Pi**2)

! *** 4 Squark Loop ***
	do m=1,2
      
	g1 = -GS**2*(USf(m,1,type,gg)*USf(l,1,type,gg)-USf(m,2,type,gg)*USf(l,2,type,gg))*
     &              (USf(m,1,type,gg)*USf(k,1,type,gg)-USf(m,2,type,gg)*USf(k,2,type,gg))
      
	ML0 = MSf(m,type,gg)
	
	PropLoop = PropLoop - (g1*A01(ML0))/(16.*Pi**2)	
	
	enddo
	
	end
	
c ************************************************

	function ChiChi2QQ_M2Propagator()

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"


	integer i,j,k,l,m,tcount,imin,imax,jmin,jmax,itypeT,itypeU,att,antitt
	double precision swap

	complex*16 DenSmV, DenSmH(6)
	double precision DenSmV1, DenSmH1(6)

	double complex PropLoop
	double precision ChiChi2QQ_M2Propagator

	double precision PropZ, PropH, Propt, Propu, Proptemp, Proptemp2
	
	PropZ = 0d0
	PropH = 0d0
	Propt = 0d0
	Propu = 0d0
	Proptemp = 0d0
	Proptemp2 = 0d0

	ChiChi2QQ_M2Propagator = 0d0

! widths of Z-boson and Higgs included in propagators	

	if (icha.eq.1) then
		DenSmV = (sqrtS**2 - MW2) + cI*MW*wW
		DenSmV1 = DenSmV*dconjg(DenSmV)/(sqrtS**2 - MW2)
	else
		DenSmV = (sqrtS**2 - MZ2) + cI*MZ*wZ
		DenSmV1 = DenSmV*dconjg(DenSmV)/(sqrtS**2 - MZ2)
	endif

! itypeT,U follows the type of the sfermion exchanged in t,u-channel interactions
! att gives the type of the antiquark final state
! antitt is the opposite type to the final state quark type

	if (icha.eq.0) then
		att   = tt
		itypeU = tt
		itypeT = tt
		antitt = tt
	else if (icha.eq.1) then
! tt is defined as normally as up-type (quark) type in chargino-neutralino co-annihilations		
		if (tt.eq.3) then 
		      att    = tt+1
		      itypeU = tt
		      itypeT = tt+1
		      antitt = tt+1
		else if (tt.eq.4) then 
		      att    = tt-1
		      itypeU = tt
		      itypeT = tt-1
		      antitt = tt-1
		endif	
	else if (icha.eq.2) then
		if (tt.eq.3) then 
		      att   = tt
		      itypeU = tt
		      itypeT = tt+1
		      antitt = tt+1
		else if (tt.eq.4) then 
		      att   = tt
		      itypeU = tt-1
		      itypeT = tt
		      antitt = tt-1
		endif	
	endif	

	if (icha.eq.1) then
! only charged Higgs exchange		
	  imin = 5
	  imax = 6 
	  jmin = 5
	  jmax = 6
	else
! exchange of neutral Higgses 		
	  imin = 1
	  imax = 4 
	  jmin = 1
	  jmax = 4	  	
	endif

	do i = imin, imax
	  DenSmH(i) = (sqrtS**2 - MHk(i)**2) + cI*MHk(i)*wHk(i)
	  DenSmH1(i) = DenSmH(i)*dconjg(DenSmH(i))/(sqrtS**2 - MHk(i)**2)
	enddo

	
!*************************************
!************* U-Loop (old T-loop) ***
!*************************************

	do k=1,2
	 do l=1,2

! Interference with Photon

	if (icha.eq.2) then
		
		ChiChi2QQ_M2Propagator = ChiChi2QQ_M2Propagator + (-4*ApL*BpL*CF*NC*(GL(indi,l,tt,gen)*(2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*HL(indj,k,tt,gen) + 
     -         (muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*ured) + (muu3**2 - ured)**2)*HR(indj,k,tt,gen)) + 
     -      GR(indi,l,tt,gen)*((muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*ured) + (muu3**2 - ured)**2)*HL(indj,k,tt,gen) + 
     -         2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*HR(indj,k,tt,gen)))*PropLoop(sqrtS**2*ured,k,l,itypeU,gen))/
     -  (sqrtS**2*(-ured + MSf(k,itypeU,gen)**2/sqrtS**2)*(ured - MSf(l,itypeU,gen)**2/sqrtS**2))

	endif

! Interference with Z/W-boson

	 ChiChi2QQ_M2Propagator = ChiChi2QQ_M2Propagator + (-4*CF*NC*(dconjg(AR)*(BR*(2*muu1*muu4*GL(indi,l,tt,gen) + (muu1**2 + muu4**2 - ured)*GR(indi,l,tt,gen))*
     -          ((muu2**2 + muu3**2 - ured)*HL(indj,k,tt,gen) + 2*muu2*muu3*HR(indj,k,tt,gen)) + 
     -         BL*(muu1*GL(indi,l,tt,gen)*(muu3*(muu2**2 + muu4**2 - tred)*HL(indj,k,tt,gen) - 
     -               muu2*(-1 + muu3**2 + muu4**2)*HR(indj,k,tt,gen)) + 
     -            muu4*GR(indi,l,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu3*HL(indj,k,tt,gen)) + 
     -               muu2*(muu1**2 + muu3**2 - tred)*HR(indj,k,tt,gen)))) + 
     -      dconjg(AL)*(BL*((muu1**2 + muu4**2 - ured)*GL(indi,l,tt,gen) + 2*muu1*muu4*GR(indi,l,tt,gen))*
     -          (2*muu2*muu3*HL(indj,k,tt,gen) + (muu2**2 + muu3**2 - ured)*HR(indj,k,tt,gen)) + 
     -         BR*(muu4*GL(indi,l,tt,gen)*(muu2*(muu1**2 + muu3**2 - tred)*HL(indj,k,tt,gen) - 
     -               (-1 + muu1**2 + muu2**2)*muu3*HR(indj,k,tt,gen)) + 
     -            muu1*GR(indi,l,tt,gen)*(-(muu2*(-1 + muu3**2 + muu4**2)*HL(indj,k,tt,gen)) + 
     -               muu3*(muu2**2 + muu4**2 - tred)*HR(indj,k,tt,gen)))))*PropLoop(sqrtS**2*ured,k,l,itypeU,gen))/
     -  (DenSmV1*(ured - MSf(k,itypeU,gen)**2/sqrtS**2)*(-ured + MSf(l,itypeU,gen)**2/sqrtS**2))

	PropZ = PropZ + ChiChi2QQ_M2Propagator - Proptemp2
	Proptemp2 = ChiChi2QQ_M2Propagator

! Interference with Higgs-boson

	  do j=jmin,jmax
		
	 ChiChi2QQ_M2Propagator = ChiChi2QQ_M2Propagator + ((CF*NC*(dconjg(CL(j))*(2*dconjg(DR(j))*(muu1*GR(indi,l,tt,gen)*
     -               (muu3*(muu2**2 + muu4**2 - tred)*HL(indj,k,tt,gen) - muu2*(-1 + muu3**2 + muu4**2)*HR(indj,k,tt,gen)) + 
     -              muu4*GL(indi,l,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu3*HL(indj,k,tt,gen)) + 
     -                 muu2*(muu1**2 + muu3**2 - tred)*HR(indj,k,tt,gen))) + 
     -           dconjg(DL(j))*(-2*muu1*muu4*GR(indi,l,tt,gen)*
     -               ((muu2**2 + muu3**2 - ured)*HL(indj,k,tt,gen) + 2*muu2*muu3*HR(indj,k,tt,gen)) + 
     -              GL(indi,l,tt,gen)*((-1 + muu3**2 + muu4**2 + muu2**2*(1 - 2*muu4**2) - 
     -                    (muu2**2 + muu3**2 + muu4**2 - tred - ured)*(tred - ured) + muu1**2*(1 - 2*muu3**2 - tred + ured))*
     -                  HL(indj,k,tt,gen) - 2*muu2*muu3*(muu1**2 + muu4**2 - ured)*HR(indj,k,tt,gen)))) + 
     -        dconjg(CR(j))*(2*dconjg(DL(j))*(muu4*GR(indi,l,tt,gen)*
     -               (muu2*(muu1**2 + muu3**2 - tred)*HL(indj,k,tt,gen) - (-1 + muu1**2 + muu2**2)*muu3*HR(indj,k,tt,gen)) + 
     -              muu1*GL(indi,l,tt,gen)*(-(muu2*(-1 + muu3**2 + muu4**2)*HL(indj,k,tt,gen)) + 
     -                 muu3*(muu2**2 + muu4**2 - tred)*HR(indj,k,tt,gen))) + 
     -           dconjg(DR(j))*(-2*muu1*muu4*GL(indi,l,tt,gen)*
     -               (2*muu2*muu3*HL(indj,k,tt,gen) + (muu2**2 + muu3**2 - ured)*HR(indj,k,tt,gen)) + 
     -              GR(indi,l,tt,gen)*(-2*muu2*muu3*(muu1**2 + muu4**2 - ured)*HL(indj,k,tt,gen) + 
     -                 (-1 + muu3**2 + muu4**2 + muu2**2*(1 - 2*muu4**2) - 
     -                    (muu2**2 + muu3**2 + muu4**2 - tred - ured)*(tred - ured) + muu1**2*(1 - 2*muu3**2 - tred + ured))*
     -                  HR(indj,k,tt,gen)))))*PropLoop(sqrtS**2*ured,k,l,itypeU,gen))/
     -    (DenSmH1(j)*(-ured + MSf(k,itypeU,gen)**2/sqrtS**2)*(ured - MSf(l,itypeU,gen)**2/sqrtS**2)))

	  enddo

	PropH = PropH + ChiChi2QQ_M2Propagator - Proptemp2
	Proptemp2 = ChiChi2QQ_M2Propagator
	
	
! Interference with new U-channel (old T-channel)

	  do j=1,2
		
	 ChiChi2QQ_M2Propagator = ChiChi2QQ_M2Propagator + (2*CF*NC*(((muu1**2 + muu4**2 - ured)*dconjg(GL(indi,j,tt,gen)) + 2*muu1*muu4*dconjg(GR(indi,j,tt,gen)))*
     -       GL(indi,l,tt,gen) + (2*muu1*muu4*dconjg(GL(indi,j,tt,gen)) + (muu1**2 + muu4**2 - ured)*dconjg(GR(indi,j,tt,gen)))*
     -       GR(indi,l,tt,gen))*(((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,j,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,j,tt,gen)))*
     -       HL(indj,k,tt,gen) + (2*muu2*muu3*dconjg(HL(indj,j,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,j,tt,gen)))*
     -       HR(indj,k,tt,gen))*PropLoop(sqrtS**2*ured,k,l,itypeU,gen))/
     -  (sqrtS**2*(ured - MSf(j,itypeU,gen)**2/sqrtS**2)*(-ured + MSf(k,itypeU,gen)**2/sqrtS**2)*
     -    (ured - MSf(l,itypeU,gen)**2/sqrtS**2))

	  enddo

	Propu = Propu + ChiChi2QQ_M2Propagator - Proptemp2
	Proptemp2 = ChiChi2QQ_M2Propagator

! Interference with new T-channel (old T-channel)

	  do j=1,2
		
	 ChiChi2QQ_M2Propagator = ChiChi2QQ_M2Propagator + (CF*NC*(dconjg(ELc(indi,j,tt,gen))*(2*dconjg(FR(indj,j,tt,gen))*
     -          (muu4*GL(indi,l,tt,gen)*(muu2*(muu1**2 + muu3**2 - tred)*HL(indj,k,tt,gen) - 
     -               (-1 + muu1**2 + muu2**2)*muu3*HR(indj,k,tt,gen)) + 
     -            muu1*GR(indi,l,tt,gen)*(-(muu2*(-1 + muu3**2 + muu4**2)*HL(indj,k,tt,gen)) + 
     -               muu3*(muu2**2 + muu4**2 - tred)*HR(indj,k,tt,gen))) + 
     -         dconjg(FL(indj,j,tt,gen))*(2*muu1*muu4*GR(indi,l,tt,gen)*
     -             ((muu2**2 + muu3**2 - ured)*HL(indj,k,tt,gen) + 2*muu2*muu3*HR(indj,k,tt,gen)) + 
     -            GL(indi,l,tt,gen)*((-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - 
     -                  (muu2**2 + muu3**2 + muu4**2 - tred)*tred + muu1**2*(1 + 2*muu2**2 - tred - ured) - 
     -                  (muu2**2 + muu3**2 + muu4**2)*ured + ured**2)*HL(indj,k,tt,gen) + 
     -               2*muu2*muu3*(muu1**2 + muu4**2 - ured)*HR(indj,k,tt,gen)))) + 
     -      dconjg(ER(indi,j,tt,gen))*(2*dconjg(FL(indj,j,tt,gen))*
     -          (muu1*GL(indi,l,tt,gen)*(muu3*(muu2**2 + muu4**2 - tred)*HL(indj,k,tt,gen) - 
     -               muu2*(-1 + muu3**2 + muu4**2)*HR(indj,k,tt,gen)) + 
     -            muu4*GR(indi,l,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu3*HL(indj,k,tt,gen)) + 
     -               muu2*(muu1**2 + muu3**2 - tred)*HR(indj,k,tt,gen))) + 
     -         dconjg(FR(indj,j,tt,gen))*(2*muu1*muu4*GL(indi,l,tt,gen)*
     -             (2*muu2*muu3*HL(indj,k,tt,gen) + (muu2**2 + muu3**2 - ured)*HR(indj,k,tt,gen)) + 
     -            GR(indi,l,tt,gen)*(2*muu2*muu3*(muu1**2 + muu4**2 - ured)*HL(indj,k,tt,gen) + 
     -               (-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - (muu2**2 + muu3**2 + muu4**2 - tred)*tred + 
     -                  muu1**2*(1 + 2*muu2**2 - tred - ured) - (muu2**2 + muu3**2 + muu4**2)*ured + ured**2)*HR(indj,k,tt,gen)))))*
     -    PropLoop(sqrtS**2*ured,k,l,itypeU,gen))/
     -  (sqrtS**2*(-tred + MSf(j,itypeT,gen)**2/sqrtS**2)*(-ured + MSf(k,itypeU,gen)**2/sqrtS**2)*
     -    (ured - MSf(l,itypeU,gen)**2/sqrtS**2))

	  enddo
	
	Propt = Propt + ChiChi2QQ_M2Propagator - Proptemp2
	Proptemp2 = ChiChi2QQ_M2Propagator
	
	enddo
	enddo
	
cc	write(*,*)'ChiChi U-Channel Propagator = ', ChiChi2QQ_M2Propagator
cc	write(*,*)'ChiChi U-Channel Propagator mit Z = ', PropZ
cc	write(*,*)'ChiChi U-Channel Propagator mit H = ', PropH
cc	write(*,*)'ChiChi U-Channel Propagator mit t = ', Propt
cc	write(*,*)'ChiChi U-Channel Propagator mit u = ', Propu
	Proptemp = ChiChi2QQ_M2Propagator
	PropZ = 0d0
	PropH = 0d0
	Propt = 0d0 
	Propu = 0d0
	

!*************************************
!************* T-Loop (old U-loop) ***
!*************************************

	do k=1,2
	 do l=1,2

! Interference with photon

	if (icha.eq.2) then

	 ChiChi2QQ_M2Propagator = ChiChi2QQ_M2Propagator + (4*ApL*BpL*CF*NC*(ELc(indi,k,tt,gen)*(2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*FL(indj,l,tt,gen) + 
     -         (muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*tred) + (muu3**2 - tred)**2)*FR(indj,l,tt,gen)) + 
     -      ER(indi,k,tt,gen)*((muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*tred) + (muu3**2 - tred)**2)*FL(indj,l,tt,gen) + 
     -         2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*FR(indj,l,tt,gen)))*PropLoop(sqrtS**2*tred,k,l,itypeT,gen))/
     -  (sqrtS**2*(-tred + MSf(k,itypeT,gen)**2/sqrtS**2)*(tred - MSf(l,itypeT,gen)**2/sqrtS**2))

	endif

! Interference with Z boson

	 ChiChi2QQ_M2Propagator = ChiChi2QQ_M2Propagator + (4*CF*NC*(dconjg(AL)*(BR*((muu1**2 + muu3**2 - tred)*ELc(indi,k,tt,gen) + 2*muu1*muu3*ER(indi,k,tt,gen))*
     -          (2*muu2*muu4*FL(indj,l,tt,gen) + (muu2**2 + muu4**2 - tred)*FR(indj,l,tt,gen)) + 
     -         BL*(muu3*ELc(indi,k,tt,gen)*(muu2*(muu1**2 + muu4**2 - ured)*FL(indj,l,tt,gen) - 
     -               (-1 + muu1**2 + muu2**2)*muu4*FR(indj,l,tt,gen)) + 
     -            muu1*ER(indi,k,tt,gen)*(-(muu2*(-1 + muu3**2 + muu4**2)*FL(indj,l,tt,gen)) + 
     -               muu4*(muu2**2 + muu3**2 - ured)*FR(indj,l,tt,gen)))) + 
     -      dconjg(AR)*(BL*(2*muu1*muu3*ELc(indi,k,tt,gen) + (muu1**2 + muu3**2 - tred)*ER(indi,k,tt,gen))*
     -          ((muu2**2 + muu4**2 - tred)*FL(indj,l,tt,gen) + 2*muu2*muu4*FR(indj,l,tt,gen)) + 
     -         BR*(muu1*ELc(indi,k,tt,gen)*(muu4*(muu2**2 + muu3**2 - ured)*FL(indj,l,tt,gen) - 
     -               muu2*(-1 + muu3**2 + muu4**2)*FR(indj,l,tt,gen)) + 
     -            muu3*ER(indi,k,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu4*FL(indj,l,tt,gen)) + 
     -               muu2*(muu1**2 + muu4**2 - ured)*FR(indj,l,tt,gen)))))*PropLoop(sqrtS**2*tred,k,l,itypeT,gen))/
     -  (DenSmV1*(tred - MSf(k,itypeT,gen)**2/sqrtS**2)*(-tred + MSf(l,itypeT,gen)**2/sqrtS**2))	

	PropZ = PropZ + ChiChi2QQ_M2Propagator - Proptemp2
	Proptemp2 = ChiChi2QQ_M2Propagator
	
! Interference with Higgs-boson

	  do j=jmin,jmax

	 ChiChi2QQ_M2Propagator = ChiChi2QQ_M2Propagator + ((CF*NC*(dconjg(CR(j))*(dconjg(DR(j))*(-2*muu1*muu3*ELc(indi,k,tt,gen)*
     -               (2*muu2*muu4*FL(indj,l,tt,gen) + (muu2**2 + muu4**2 - tred)*FR(indj,l,tt,gen)) + 
     -              ER(indi,k,tt,gen)*(-2*muu2*muu4*(muu1**2 + muu3**2 - tred)*FL(indj,l,tt,gen) + 
     -                 (-1 + muu3**2 + muu2**2*(1 - 2*muu3**2) + muu4**2 + 
     -                    (muu2**2 + muu3**2 + muu4**2 - tred - ured)*(tred - ured) + muu1**2*(1 - 2*muu4**2 + tred - ured))*
     -                  FR(indj,l,tt,gen))) + 2*dconjg(DL(j))*
     -            (muu3*ER(indi,k,tt,gen)*(muu2*(muu1**2 + muu4**2 - ured)*FL(indj,l,tt,gen) - 
     -                 (-1 + muu1**2 + muu2**2)*muu4*FR(indj,l,tt,gen)) + 
     -              muu1*ELc(indi,k,tt,gen)*(-(muu2*(-1 + muu3**2 + muu4**2)*FL(indj,l,tt,gen)) + 
     -                 muu4*(muu2**2 + muu3**2 - ured)*FR(indj,l,tt,gen)))) + 
     -        dconjg(CL(j))*(dconjg(DL(j))*(-2*muu1*muu3*ER(indi,k,tt,gen)*
     -               ((muu2**2 + muu4**2 - tred)*FL(indj,l,tt,gen) + 2*muu2*muu4*FR(indj,l,tt,gen)) + 
     -              ELc(indi,k,tt,gen)*((-1 + muu3**2 + muu2**2*(1 - 2*muu3**2) + muu4**2 + 
     -                    (muu2**2 + muu3**2 + muu4**2 - tred - ured)*(tred - ured) + muu1**2*(1 - 2*muu4**2 + tred - ured))*
     -                  FL(indj,l,tt,gen) - 2*muu2*muu4*(muu1**2 + muu3**2 - tred)*FR(indj,l,tt,gen))) + 
     -           2*dconjg(DR(j))*(muu1*ER(indi,k,tt,gen)*
     -               (muu4*(muu2**2 + muu3**2 - ured)*FL(indj,l,tt,gen) - muu2*(-1 + muu3**2 + muu4**2)*FR(indj,l,tt,gen)) + 
     -              muu3*ELc(indi,k,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu4*FL(indj,l,tt,gen)) + 
     -                 muu2*(muu1**2 + muu4**2 - ured)*FR(indj,l,tt,gen)))))*PropLoop(sqrtS**2*tred,k,l,itypeT,gen))/
     -    (DenSmH1(j)*(-tred + MSf(k,itypeT,gen)**2/sqrtS**2)*(tred - MSf(l,itypeT,gen)**2/sqrtS**2)))		

	  enddo
	
	PropH = PropH + ChiChi2QQ_M2Propagator - Proptemp2
	Proptemp2 = ChiChi2QQ_M2Propagator
	
! Interference with new U-channel (old T-channel)

	  do j=1,2
		
	 ChiChi2QQ_M2Propagator = ChiChi2QQ_M2Propagator - (CF*NC*(ELc(indi,k,tt,gen)*((2*muu1*muu4*dconjg(GR(indi,j,tt,gen))*
     -             ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,j,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,j,tt,gen))) + 
     -            dconjg(GL(indi,j,tt,gen))*((-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - 
     -                  (muu2**2 + muu3**2 + muu4**2 - tred)*tred + muu1**2*(1 + 2*muu2**2 - tred - ured) - 
     -                  (muu2**2 + muu3**2 + muu4**2)*ured + ured**2)*dconjg(HL(indj,j,tt,gen)) + 
     -               2*muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,j,tt,gen))))*FL(indj,l,tt,gen) + 
     -         2*(muu4*dconjg(GL(indi,j,tt,gen))*(muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,j,tt,gen)) - 
     -               (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,j,tt,gen))) + 
     -            muu1*dconjg(GR(indi,j,tt,gen))*(-(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,j,tt,gen))) + 
     -               muu3*(muu2**2 + muu4**2 - tred)*dconjg(HR(indj,j,tt,gen))))*FR(indj,l,tt,gen)) + 
     -      ER(indi,k,tt,gen)*(2*(muu1*dconjg(GL(indi,j,tt,gen))*
     -             (muu3*(muu2**2 + muu4**2 - tred)*dconjg(HL(indj,j,tt,gen)) - 
     -               muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,j,tt,gen))) + 
     -            muu4*dconjg(GR(indi,j,tt,gen))*(-((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,j,tt,gen))) + 
     -               muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,j,tt,gen))))*FL(indj,l,tt,gen) + 
     -         (2*muu1*muu4*dconjg(GL(indi,j,tt,gen))*
     -             (2*muu2*muu3*dconjg(HL(indj,j,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,j,tt,gen))) + 
     -            dconjg(GR(indi,j,tt,gen))*(2*muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,j,tt,gen)) + 
     -               (-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - (muu2**2 + muu3**2 + muu4**2 - tred)*tred + 
     -                  muu1**2*(1 + 2*muu2**2 - tred - ured) - (muu2**2 + muu3**2 + muu4**2)*ured + ured**2)*
     -                dconjg(HR(indj,j,tt,gen))))*FR(indj,l,tt,gen)))*PropLoop(sqrtS**2*tred,k,l,itypeT,gen))/
     -  (sqrtS**2*(-ured + MSf(j,itypeU,gen)**2/sqrtS**2)*(-tred + MSf(k,itypeT,gen)**2/sqrtS**2)*
     -    (-tred + MSf(l,itypeT,gen)**2/sqrtS**2))		

	  enddo
	
	Propu = Propu + ChiChi2QQ_M2Propagator - Proptemp2
	Proptemp2 = ChiChi2QQ_M2Propagator
	
! Interference with new T-channel (old U-channel)

	  do j=1,2
		
	 ChiChi2QQ_M2Propagator = ChiChi2QQ_M2Propagator + (2*CF*NC*(((muu1**2 + muu3**2 - tred)*dconjg(ELc(indi,j,tt,gen)) + 2*muu1*muu3*dconjg(ER(indi,j,tt,gen)))*
     -       ELc(indi,k,tt,gen) + (2*muu1*muu3*dconjg(ELc(indi,j,tt,gen)) + (muu1**2 + muu3**2 - tred)*dconjg(ER(indi,j,tt,gen)))*
     -       ER(indi,k,tt,gen))*(((muu2**2 + muu4**2 - tred)*dconjg(FL(indj,j,tt,gen)) + 2*muu2*muu4*dconjg(FR(indj,j,tt,gen)))*
     -       FL(indj,l,tt,gen) + (2*muu2*muu4*dconjg(FL(indj,j,tt,gen)) + (muu2**2 + muu4**2 - tred)*dconjg(FR(indj,j,tt,gen)))*
     -       FR(indj,l,tt,gen))*PropLoop(sqrtS**2*tred,k,l,itypeT,gen))/
     -  (sqrtS**2*(-tred + MSf(j,itypeT,gen)**2/sqrtS**2)*(tred - MSf(k,itypeT,gen)**2/sqrtS**2)*
     -    (tred - MSf(l,itypeT,gen)**2/sqrtS**2))		

	  enddo

	Propt = Propt + ChiChi2QQ_M2Propagator - Proptemp2
	Proptemp2 = ChiChi2QQ_M2Propagator

	 enddo
	enddo
	
cc	write(*,*)'ChiChi T-Channel Propagator = ', ChiChi2QQ_M2Propagator - Proptemp
cc	write(*,*)'ChiChi T-Channel Propagator mit Z = ', PropZ
cc	write(*,*)'ChiChi T-Channel Propagator mit H = ', PropH
cc	write(*,*)'ChiChi T-Channel Propagator mit t = ', Propt
cc	write(*,*)'ChiChi T-Channel Propagator mit u = ', Propu

	
	end	
