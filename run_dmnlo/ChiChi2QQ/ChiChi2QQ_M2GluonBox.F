************************************************

	function ChiChi2QQ_M2BoxGluon()

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"


	integer i,j,k,l,imin,imax,jmin,jmax,itypeT,itypeU,att
	double complex g1L,g1R,g2L,g2R,g3L,g4L
	double precision swap, tempBox, tempBox2, tempBoxZ, tempBoxH, tempBoxu, tempBoxt
	
	complex*16 DenSmV, DenSmH(6)
	double precision DenSmV1, DenSmH1(6)

	double precision ChiChi2QQ_M2BoxGluon
	external myDget

	ChiChi2QQ_M2BoxGluon = 0d0
	tempBox2 = 0d0
	tempBoxZ = 0d0
	tempBoxH = 0d0
	tempBoxu = 0d0
	tempBoxt = 0d0

! widths of Z-boson and Higgs included in propagators	

	if (icha.eq.1) then
		DenSmV = (sqrtS**2 - MW2) + cI*MW*wW
		DenSmV1 = DenSmV*dconjg(DenSmV)/(sqrtS**2 - MW2)
	else
		DenSmV = (sqrtS**2 - MZ2) + cI*MZ*wZ
		DenSmV1 = DenSmV*dconjg(DenSmV)/(sqrtS**2 - MZ2)
	endif

! itypeU,U follows the type of the sfermion exchanged in t,u-channel interactions
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

!***********************************
!******* Gluon U-channel Box *******
!***********************************
	
	do k = 1,2
		
	g1L = GL(indi,k,tt,gen)
	g1R = GR(indi,k,tt,gen)
	g2L = HL(indj,k,tt,gen) 
	g2R = HR(indj,k,tt,gen)
	g3L = -GS
	g4L = -GS
	
	call myDget(muu4**2*sqrtS**2,muu3**2*sqrtS**2,muu2**2*sqrtS**2,muu1**2*sqrtS**2,sqrtS**2,ured*sqrtS**2,muu4**2*sqrtS**2,0d0,muu3**2*sqrtS**2,MSf2(k,itypeU,gen))
	
! Photon contr.

	if (icha.eq.2) then

	
	ChiChi2QQ_M2BoxGluon = ChiChi2QQ_M2BoxGluon + (ApL*BpL*CF*g3L*g4L*NC*sqrtS**4*(2*muu3*(g1L*(g2R*muu3*
     &             (muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*ured) + (1 + ured)**2 - muu3**2*(1 + 2*ured)) + 
     &            g2L*muu1*(1 - 4*muu3**4 + muu3**2*(4 - 4*muu1**2 + 4*ured))) + 
     &         g1R*(g2L*muu3*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*ured) + (1 + ured)**2 - muu3**2*(1 + 2*ured)) + 
     &            g2R*muu1*(1 - 4*muu3**4 + muu3**2*(4 - 4*muu1**2 + 4*ured))))*xd0 - 
     &      (4*(g1R*(4*g2R*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured) + 
     &              g2L*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + ured**2)) + 
     &           g1L*(4*g2L*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured) + 
     &              g2R*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + ured**2)))*xd00)/sqrtS**2 - 
     &      2*(g1R*(g2R*muu1*muu3*(1 + 4*muu3**4 + muu1**2*(2 + 4*muu3**2) - 2*ured - 4*muu3**2*(2 + ured)) + 
     &            g2L*(-4*muu3**6 + muu1**4*(1 - 4*muu3**2) + ured**2 - 2*muu3**2*ured*(3 + 2*ured) + muu3**4*(1 + 8*ured) + 
     &               muu1**2*(1 - 8*muu3**4 - 2*ured + muu3**2*(4 + 8*ured)))) + 
     &         g1L*(g2L*muu1*muu3*(1 + 4*muu3**4 + muu1**2*(2 + 4*muu3**2) - 2*ured - 4*muu3**2*(2 + ured)) + 
     &            g2R*(-4*muu3**6 + muu1**4*(1 - 4*muu3**2) + ured**2 - 2*muu3**2*ured*(3 + 2*ured) + muu3**4*(1 + 8*ured) + 
     &               muu1**2*(1 - 8*muu3**4 - 2*ured + muu3**2*(4 + 8*ured)))))*xd1 + 
     &      2*muu3**2*(g1R*(2*g2R*muu1*muu3*(1 - 2*muu1**2 - 2*muu3**2 + 2*ured) + 
     &            g2L*(muu1**4 + muu3**4 + muu1**2*(-3 + 2*muu3**2 - 2*ured) + ured*(2 + ured) - muu3**2*(1 + 2*ured))) + 
     &         g1L*(2*g2L*muu1*muu3*(1 - 2*muu1**2 - 2*muu3**2 + 2*ured) + 
     &            g2R*(muu1**4 + muu3**4 + muu1**2*(-3 + 2*muu3**2 - 2*ured) + ured*(2 + ured) - muu3**2*(1 + 2*ured))))*xd11 - 
     &      2*(g1R*(2*g2R*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured) + 
     &            g2L*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + ured**2)) + 
     &         g1L*(2*g2L*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured) + 
     &            g2R*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + ured**2)))*xd12 - 
     &      2*(g1R*(2*g2R*muu1*muu3*(2*muu1**4 + 2*muu3**4 + muu1**2*(-1 + 4*muu3**2 - 4*ured) + ured*(1 + 2*ured) - muu3**2*(1 + 4*ured)) + 
     &            g2L*(muu1**4*(1 - 4*muu3**2) + (muu3**2 - ured)**2 - 2*muu1**2*(2*muu3**4 + ured - muu3**2*(3 + 2*ured)))) + 
     &         g1L*(2*g2L*muu1*muu3*(2*muu1**4 + 2*muu3**4 + muu1**2*(-1 + 4*muu3**2 - 4*ured) + ured*(1 + 2*ured) - muu3**2*(1 + 4*ured)) + 
     &            g2R*(muu1**4*(1 - 4*muu3**2) + (muu3**2 - ured)**2 - 2*muu1**2*(2*muu3**4 + ured - muu3**2*(3 + 2*ured)))))*xd13 - 
     &      2*(g1R*(2*g2R*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured) + 
     &            g2L*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + ured**2)) + 
     &         g1L*(2*g2L*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured) + 
     &            g2R*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + ured**2)))*xd2 - 
     &      2*(g1R*(2*g2R*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured) + 
     &            g2L*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + ured**2)) + 
     &         g1L*(2*g2L*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured) + 
     &            g2R*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + ured**2)))*xd22 - 
     &      2*(g1R*(2*g2R*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured) + 
     &            g2L*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + ured**2)) + 
     &         g1L*(2*g2L*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured) + 
     &            g2R*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + ured**2)))*xd23 - 
     &      2*muu1*(g1R*(g2R*muu3*(1 + 2*muu1**4 + 2*muu3**4 + muu1**2*(-2 + 4*muu3**2 - 4*ured) + 4*ured + 2*ured**2 - 4*muu3**2*(1 + ured)) + 
     &            g2L*muu1*(1 - 8*muu3**4 + muu3**2*(6 - 8*muu1**2 + 8*ured))) + 
     &         g1L*(g2L*muu3*(1 + 2*muu1**4 + 2*muu3**4 + muu1**2*(-2 + 4*muu3**2 - 4*ured) + 4*ured + 2*ured**2 - 4*muu3**2*(1 + ured)) + 
     &            g2R*muu1*(1 - 8*muu3**4 + muu3**2*(6 - 8*muu1**2 + 8*ured))))*xd3 + 
     &      2*muu1**2*(g1R*(2*g2R*muu1*muu3*(1 - 2*muu1**2 - 2*muu3**2 + 2*ured) + 
     &            g2L*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*ured) + ured*(2 + ured) - muu3**2*(3 + 2*ured))) + 
     &         g1L*(2*g2L*muu1*muu3*(1 - 2*muu1**2 - 2*muu3**2 + 2*ured) + 
     &            g2R*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*ured) + ured*(2 + ured) - muu3**2*(3 + 2*ured))))*xd33))/(4.*Pi**2)
	
	endif

! Zboson contr.
	
	ChiChi2QQ_M2BoxGluon = ChiChi2QQ_M2BoxGluon -(CF*g3L*g4L*NC*sqrtS**6*(2*muu1*xd33*((-(BR*(g1L*muu1*muu4*(2*g2R*(-1 + muu1**2 + muu2**2)*muu3 - g2L*muu2*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                  g1R*(-(g2R*muu1**2*muu3*(-1 + muu1**2 + muu3**2 - ured)) + g2L*muu2*(muu1**2 + muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))))
     &              + BL*muu1*(g1R*muu1*muu4*(2*g2L*muu2*muu3 + g2R*(muu2**2 + muu3**2 - ured)) + 
     &                g1L*(g2L*muu2*muu3*(muu1**2 + muu4**2 - ured) + g2R*(muu1**2*muu3**2 + muu2**2*muu4**2 - ured))))*dconjg(AL) + 
     &          (-(BL*(g1R*muu1*muu4*(2*g2L*(-1 + muu1**2 + muu2**2)*muu3 - g2R*muu2*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                  g1L*(-(g2L*muu1**2*muu3*(-1 + muu1**2 + muu3**2 - ured)) + g2R*muu2*(muu1**2 + muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))))
     &              + BR*muu1*(g1L*muu1*muu4*(2*g2R*muu2*muu3 + g2L*(muu2**2 + muu3**2 - ured)) + 
     &                g1R*(g2R*muu2*muu3*(muu1**2 + muu4**2 - ured) + g2L*(muu1**2*muu3**2 + muu2**2*muu4**2 - ured))))*dconjg(AR)) + 
     &       2*muu4*xd11*((-(BR*(g1R*muu1*muu4*(2*g2L*muu2*(-1 + muu3**2 + muu4**2) - g2R*muu3*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                  g1L*(g2R*muu3*(-1 + muu1**2 + muu3**2 - ured)*(muu1**2 + muu4**2 - ured) - g2L*muu2*muu4**2*(-1 + muu2**2 + muu4**2 - ured)))) + 
     &             BL*muu4*(g1R*muu1*muu4*(2*g2L*muu2*muu3 + g2R*(muu2**2 + muu3**2 - ured)) + 
     &                g1L*(g2L*muu2*muu3*(muu1**2 + muu4**2 - ured) + g2R*(muu1**2*muu3**2 + muu2**2*muu4**2 - ured))))*dconjg(AL) + 
     &          (-(BL*(g1L*muu1*muu4*(2*g2R*muu2*(-1 + muu3**2 + muu4**2) - g2L*muu3*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                  g1R*(g2L*muu3*(-1 + muu1**2 + muu3**2 - ured)*(muu1**2 + muu4**2 - ured) - g2R*muu2*muu4**2*(-1 + muu2**2 + muu4**2 - ured)))) + 
     &             BR*muu4*(g1L*muu1*muu4*(2*g2R*muu2*muu3 + g2L*(muu2**2 + muu3**2 - ured)) + 
     &                g1R*(g2R*muu2*muu3*(muu1**2 + muu4**2 - ured) + g2L*(muu1**2*muu3**2 + muu2**2*muu4**2 - ured))))*dconjg(AR)) + 
     &       (4*xd00*((-(BR*(g1R*muu1*(g2L*muu2*(-1 + muu3**2 + muu4**2) - 2*g2R*muu3*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                    g1L*muu4*(g2R*(-1 + muu1**2 + muu2**2)*muu3 - 2*g2L*muu2*(-1 + muu2**2 + muu4**2 - ured)))) + 
     &               BL*(2*g2L*muu2*muu3 + g2R*(muu2**2 + muu3**2 - ured))*(2*g1R*muu1*muu4 + g1L*(muu1**2 + muu4**2 - ured)))*dconjg(AL) + 
     &            (-(BL*(g1L*muu1*(g2R*muu2*(-1 + muu3**2 + muu4**2) - 2*g2L*muu3*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                    g1R*muu4*(g2L*(-1 + muu1**2 + muu2**2)*muu3 - 2*g2R*muu2*(-1 + muu2**2 + muu4**2 - ured)))) + 
     &               BR*(2*g2R*muu2*muu3 + g2L*(muu2**2 + muu3**2 - ured))*(2*g1L*muu1*muu4 + g1R*(muu1**2 + muu4**2 - ured)))*dconjg(AR)))/sqrtS**2 + 
     &       2*xd13*((BL*(g1R*muu1*muu4*(2*g2L*muu2*muu3 + g2R*(muu2**2 + muu3**2 - ured)) + 
     &                g1L*(g2L*muu2*muu3*(muu1**2 + muu4**2 - ured) + g2R*(muu1**2*muu3**2 + muu2**2*muu4**2 - ured)))*(muu1**2 + muu4**2 - ured) - 
     &             BR*(g1L*muu4*(g2R*muu3*(3*muu1**4 + muu1**2*(-3 + muu2**2 + 2*muu3**2 + muu4**2 - 3*ured) + (-1 + muu2**2)*(muu4**2 - ured)) - 
     &                   g2L*muu2*(muu1**2 + muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                g1R*muu1*(-(g2R*muu3*(-1 + muu1**2 + muu3**2 - ured)*(muu1**2 + muu4**2 - ured)) + 
     &                   g2L*muu2*(3*muu4**4 + muu1**2*(-1 + muu3**2 + muu4**2) + muu4**2*(-3 + 2*muu2**2 + muu3**2 - 3*ured) + ured - muu3**2*ured))))*
     &           dconjg(AL) + (BR*(g1L*muu1*muu4*(2*g2R*muu2*muu3 + g2L*(muu2**2 + muu3**2 - ured)) + 
     &                g1R*(g2R*muu2*muu3*(muu1**2 + muu4**2 - ured) + g2L*(muu1**2*muu3**2 + muu2**2*muu4**2 - ured)))*(muu1**2 + muu4**2 - ured) - 
     &             BL*(g1R*muu4*(g2L*muu3*(3*muu1**4 + muu1**2*(-3 + muu2**2 + 2*muu3**2 + muu4**2 - 3*ured) + (-1 + muu2**2)*(muu4**2 - ured)) - 
     &                   g2R*muu2*(muu1**2 + muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                g1L*muu1*(-(g2L*muu3*(-1 + muu1**2 + muu3**2 - ured)*(muu1**2 + muu4**2 - ured)) + 
     &                   g2R*muu2*(3*muu4**4 + muu1**2*(-1 + muu3**2 + muu4**2) + muu4**2*(-3 + 2*muu2**2 + muu3**2 - 3*ured) + ured - muu3**2*ured))))*
     &           dconjg(AR)) - 2*xd3*((BR*(-(g1L*muu4*(g2L*muu2*
     &                      ((muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured) + muu1**2*(-2 + muu2**2 + 3*muu3**2 + 2*muu4**2 - ured)) + 
     &                     g2R*muu1**2*muu3*(3 - 3*muu1**2 + muu2**2 - 2*muu3**2 + 2*ured))) - 
     &                g1R*muu1*(g2L*muu2*(1 - (1 + muu1**2)*muu3**2 + muu4**2 - 2*muu4**4 + muu2**2*(-1 + muu3**2 - 2*muu4**2) + ured + 2*muu4**2*ured) + 
     &                   g2R*muu3*(1 + muu1**4 - muu3**2 - 3*muu4**2 + muu3**2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) + 
     &                      muu1**2*(-2 + muu3**2 + 3*muu4**2 - 2*ured) + 2*ured - muu3**2*ured - muu4**2*ured + ured**2))) + 
     &             BL*muu1*(g1L*muu1*(-1 + muu3**2 + muu4**2)*(2*g2L*muu2*muu3 + g2R*(muu2**2 + muu3**2 - ured)) + 
     &                g1R*muu4*(2*g2L*muu2*muu3*(-1 + muu2**2 + muu4**2 - ured) + 
     &                   g2R*(muu2**4 + (-1 + muu1**2)*(muu3**2 + ured) - muu2**2*(1 + muu1**2 - muu3**2 - 2*muu4**2 + ured)))))*dconjg(AL) + 
     &          (BL*(-(g1R*muu4*(g2R*muu2*((muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured) + muu1**2*(-2 + muu2**2 + 3*muu3**2 + 2*muu4**2 - ured)) + 
     &                     g2L*muu1**2*muu3*(3 - 3*muu1**2 + muu2**2 - 2*muu3**2 + 2*ured))) - 
     &                g1L*muu1*(g2R*muu2*(1 - (1 + muu1**2)*muu3**2 + muu4**2 - 2*muu4**4 + muu2**2*(-1 + muu3**2 - 2*muu4**2) + ured + 2*muu4**2*ured) + 
     &                   g2L*muu3*(1 + muu1**4 - muu3**2 - 3*muu4**2 + muu3**2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) + 
     &                      muu1**2*(-2 + muu3**2 + 3*muu4**2 - 2*ured) + 2*ured - muu3**2*ured - muu4**2*ured + ured**2))) + 
     &             BR*muu1*(g1R*muu1*(-1 + muu3**2 + muu4**2)*(2*g2R*muu2*muu3 + g2L*(muu2**2 + muu3**2 - ured)) + 
     &                g1L*muu4*(2*g2R*muu2*muu3*(-1 + muu2**2 + muu4**2 - ured) + 
     &                   g2L*(muu2**4 + (-1 + muu1**2)*(muu3**2 + ured) - muu2**2*(1 + muu1**2 - muu3**2 - 2*muu4**2 + ured)))))*dconjg(AR)) - 
     &       2*muu4*xd0*((BR*(g1R*muu1*muu4*(g2L*muu2*(-1 + 3*muu3**2 + muu4**2) + g2R*muu3*(2 - 2*muu1**2 + muu2**2 - muu3**2 + ured)) + 
     &                g1L*(g2R*muu3*(muu1**4 + muu2**2*(-1 + muu3**2 - muu4**2) + muu1**2*(-2 + muu3**2 + muu4**2 - 2*ured) + 
     &                      (-1 + muu3**2 - ured)*(-1 + muu4**2 - ured)) + 
     &                   g2L*muu2*(1 - (1 + muu1**2)*muu3**2 - muu4**4 + muu2**2*(-1 + muu3**2 - muu4**2) + ured + muu4**2*ured))) + 
     &             BL*(-(g1R*muu1*(-1 + muu3**2)*(2*g2L*muu2*muu3 + g2R*(muu2**2 + muu3**2 - ured))) + 
     &                g1L*muu4*(g2L*muu2*muu3*(2 + muu1**2 - 2*muu2**2 - muu4**2 + ured) + 
     &                   g2R*(-muu2**4 + muu3**2 - muu1**2*ured + muu2**2*(1 + muu1**2 - muu3**2 - muu4**2 + ured)))))*dconjg(AL) + 
     &          (BL*(g1L*muu1*muu4*(g2R*muu2*(-1 + 3*muu3**2 + muu4**2) + g2L*muu3*(2 - 2*muu1**2 + muu2**2 - muu3**2 + ured)) + 
     &                g1R*(g2L*muu3*(muu1**4 + muu2**2*(-1 + muu3**2 - muu4**2) + muu1**2*(-2 + muu3**2 + muu4**2 - 2*ured) + 
     &                      (-1 + muu3**2 - ured)*(-1 + muu4**2 - ured)) + 
     &                   g2R*muu2*(1 - (1 + muu1**2)*muu3**2 - muu4**4 + muu2**2*(-1 + muu3**2 - muu4**2) + ured + muu4**2*ured))) + 
     &             BR*(-(g1L*muu1*(-1 + muu3**2)*(2*g2R*muu2*muu3 + g2L*(muu2**2 + muu3**2 - ured))) + 
     &                g1R*muu4*(g2R*muu2*muu3*(2 + muu1**2 - 2*muu2**2 - muu4**2 + ured) + 
     &                   g2L*(-muu2**4 + muu3**2 - muu1**2*ured + muu2**2*(1 + muu1**2 - muu3**2 - muu4**2 + ured)))))*dconjg(AR)) - 
     &       2*xd2*(-((BR*(g1R*muu1*(g2L*muu2*(1 + muu4**2 - 2*muu4**4 + muu3**2*(-1 + 2*muu4**2)) + 
     &                     g2R*muu3*(-muu3**4 + muu4**2*(-3 + 3*muu1**2 - 2*muu2**2 - ured) + muu3**2*(1 - muu1**2 + muu4**2 + ured))) + 
     &                  g1L*muu4*(g2L*muu2*(muu2**2*(2 - 3*muu3**2 + muu4**2) + (2 + muu4**2)*(-1 + muu4**2 - ured) + 
     &                        muu3**2*(3 + 2*muu1**2 - muu4**2 + ured)) + 
     &                     g2R*muu3*(1 - 2*muu1**4 + muu2**2*(-1 + 2*muu4**2 - 2*ured) + muu1**2*(1 + 2*muu2**2 - 2*muu3**2 + 2*ured)))) + 
     &               BL*(2*g1R*muu1*muu3*muu4*(g2R*(1 - muu1**2 + muu2**2)*muu3 + g2L*muu2*(1 + muu3**2 - muu4**2)) - 
     &                  g1L*(2*g2L*muu2*muu3*(-(muu2**2*muu4**2) + muu1**2*(-1 + muu4**2) + ured) + 
     &                     g2R*(-2*muu2**4*muu4**2 + muu1**2*(-1 + 2*muu4**2)*(muu2**2 + muu3**2 - ured) + (muu3**2 - ured)*(-muu4**2 + ured) + 
     &                        muu2**2*(muu4**2 - 2*muu4**4 + ured + 2*muu4**2*ured)))))*dconjg(AL)) - 
     &          (BL*(g1L*muu1*(g2R*muu2*(1 + muu4**2 - 2*muu4**4 + muu3**2*(-1 + 2*muu4**2)) + 
     &                   g2L*muu3*(-muu3**4 + muu4**2*(-3 + 3*muu1**2 - 2*muu2**2 - ured) + muu3**2*(1 - muu1**2 + muu4**2 + ured))) + 
     &                g1R*muu4*(g2R*muu2*(muu2**2*(2 - 3*muu3**2 + muu4**2) + (2 + muu4**2)*(-1 + muu4**2 - ured) + 
     &                      muu3**2*(3 + 2*muu1**2 - muu4**2 + ured)) + 
     &                   g2L*muu3*(1 - 2*muu1**4 + muu2**2*(-1 + 2*muu4**2 - 2*ured) + muu1**2*(1 + 2*muu2**2 - 2*muu3**2 + 2*ured)))) + 
     &             BR*(2*g1L*muu1*muu3*muu4*(g2L*(1 - muu1**2 + muu2**2)*muu3 + g2R*muu2*(1 + muu3**2 - muu4**2)) - 
     &                g1R*(2*g2R*muu2*muu3*(-(muu2**2*muu4**2) + muu1**2*(-1 + muu4**2) + ured) + 
     &                   g2L*(-2*muu2**4*muu4**2 + muu1**2*(-1 + 2*muu4**2)*(muu2**2 + muu3**2 - ured) + (muu3**2 - ured)*(-muu4**2 + ured) + 
     &                      muu2**2*(muu4**2 - 2*muu4**4 + ured + 2*muu4**2*ured)))))*dconjg(AR)) + 
     &       2*xd12*((BR*(g1R*muu1*(g2L*muu2*(1 + muu3**4 + 2*muu4**2 - 3*muu4**4 + 2*muu3**2*(-1 + muu4**2)) - 
     &                   g2R*muu3*(-1 + muu3**2 - muu4**2)*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                g1L*muu4*(-(g2L*muu2*(-1 + muu3**2 - muu4**2)*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                   g2R*muu3*(1 - 2*muu1**4 - muu3**2 + muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2 - 2*ured) + 
     &                      muu1**2*(1 + 2*muu2**2 - muu3**2 - muu4**2 + 2*ured)))) + 
     &             BL*(g1R*muu1*muu4*(1 - muu3**2 + muu4**2)*(2*g2L*muu2*muu3 + g2R*(muu2**2 + muu3**2 - ured)) + 
     &                g1L*(-(g2L*muu2*muu3*(-1 + muu3**2 - muu4**2)*(muu1**2 + muu4**2 - ured)) + 
     &                   g2R*(muu1**4*muu3**2 + muu2**4*muu4**2 + muu3**2*muu4**2 - 2*muu4**2*ured + ured**2 - 
     &                      muu2**2*(-2*muu4**4 + muu3**2*(2*muu4**2 - ured) + ured + muu4**2*ured) - 
     &                      muu1**2*(muu2**2*(-1 + muu3**2 + muu4**2) + (1 + muu3**2 - muu4**2)*ured)))))*dconjg(AL) + 
     &          (BL*(g1L*muu1*(g2R*muu2*(1 + muu3**4 + 2*muu4**2 - 3*muu4**4 + 2*muu3**2*(-1 + muu4**2)) - 
     &                   g2L*muu3*(-1 + muu3**2 - muu4**2)*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                g1R*muu4*(-(g2R*muu2*(-1 + muu3**2 - muu4**2)*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                   g2L*muu3*(1 - 2*muu1**4 - muu3**2 + muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2 - 2*ured) + 
     &                      muu1**2*(1 + 2*muu2**2 - muu3**2 - muu4**2 + 2*ured)))) + 
     &             BR*(g1L*muu1*muu4*(1 - muu3**2 + muu4**2)*(2*g2R*muu2*muu3 + g2L*(muu2**2 + muu3**2 - ured)) + 
     &                g1R*(-(g2R*muu2*muu3*(-1 + muu3**2 - muu4**2)*(muu1**2 + muu4**2 - ured)) + 
     &                   g2L*(muu1**4*muu3**2 + muu2**4*muu4**2 + muu3**2*muu4**2 - 2*muu4**2*ured + ured**2 - 
     &                      muu2**2*(-2*muu4**4 + muu3**2*(2*muu4**2 - ured) + ured + muu4**2*ured) - 
     &                      muu1**2*(muu2**2*(-1 + muu3**2 + muu4**2) + (1 + muu3**2 - muu4**2)*ured)))))*dconjg(AR)) + 
     &       2*xd1*((BR*(g1R*muu1*(g2L*muu2*(1 + muu3**4 + muu4**2 - 2*muu4**4 - muu3**2*(2 + 3*muu4**2)) + 
     &                   g2R*muu3*(1 + muu3**4 - 4*muu4**2 - muu2**2*muu4**2 + muu1**2*(-1 + muu3**2 + 4*muu4**2) + muu3**2*(-2 + 3*muu4**2 - ured) + 
     &                      ured - 3*muu4**2*ured)) + g1L*muu4*
     &                 (g2L*muu2*(muu1**2*muu3**2 + 3*muu4**4 + muu4**2*(-3 + 3*muu2**2 + muu3**2 - 3*ured) - muu3**2*ured) + 
     &                   g2R*muu3*(-2*muu1**4 - muu1**2*(-2 + muu3**2 + muu4**2 - 4*ured) + (-3 + 2*muu3**2 - 2*ured)*ured + 
     &                      muu4**2*(1 + 2*muu2**2 - 2*muu3**2 + 2*ured)))) - 
     &             BL*(g1R*muu1*muu4*(-1 + muu3**2 + muu4**2)*(2*g2L*muu2*muu3 + g2R*(muu2**2 + muu3**2 - ured)) + 
     &                g1L*(2*g2L*muu2*muu3*(-(muu2**2*muu4**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + ured - muu3**2*ured) + 
     &                   g2R*(-(muu2**4*muu4**2) + muu3**4*muu4**2 + muu3**2*muu4**4 + muu3**2*ured - muu3**4*ured + 2*muu4**2*ured - 
     &                      2*muu3**2*muu4**2*ured - muu4**4*ured - ured**2 + muu3**2*ured**2 + muu4**2*ured**2 - 
     &                      muu2**2*(muu4**4 + (-1 + muu3**2)*ured) + 
     &                      muu1**2*(muu3**4 + muu2**2*(-1 + muu3**2 + 2*muu4**2) + ured - 2*muu4**2*ured - muu3**2*(1 + ured))))))*dconjg(AL) + 
     &          (BL*(g1L*muu1*(g2R*muu2*(1 + muu3**4 + muu4**2 - 2*muu4**4 - muu3**2*(2 + 3*muu4**2)) + 
     &                   g2L*muu3*(1 + muu3**4 - 4*muu4**2 - muu2**2*muu4**2 + muu1**2*(-1 + muu3**2 + 4*muu4**2) + muu3**2*(-2 + 3*muu4**2 - ured) + 
     &                      ured - 3*muu4**2*ured)) + g1R*muu4*
     &                 (g2R*muu2*(muu1**2*muu3**2 + 3*muu4**4 + muu4**2*(-3 + 3*muu2**2 + muu3**2 - 3*ured) - muu3**2*ured) + 
     &                   g2L*muu3*(-2*muu1**4 - muu1**2*(-2 + muu3**2 + muu4**2 - 4*ured) + (-3 + 2*muu3**2 - 2*ured)*ured + 
     &                      muu4**2*(1 + 2*muu2**2 - 2*muu3**2 + 2*ured)))) - 
     &             BR*(g1L*muu1*muu4*(-1 + muu3**2 + muu4**2)*(2*g2R*muu2*muu3 + g2L*(muu2**2 + muu3**2 - ured)) + 
     &                g1R*(2*g2R*muu2*muu3*(-(muu2**2*muu4**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + ured - muu3**2*ured) + 
     &                   g2L*(-(muu2**4*muu4**2) + muu3**4*muu4**2 + muu3**2*muu4**4 + muu3**2*ured - muu3**4*ured + 2*muu4**2*ured - 
     &                      2*muu3**2*muu4**2*ured - muu4**4*ured - ured**2 + muu3**2*ured**2 + muu4**2*ured**2 - 
     &                      muu2**2*(muu4**4 + (-1 + muu3**2)*ured) + 
     &                      muu1**2*(muu3**4 + muu2**2*(-1 + muu3**2 + 2*muu4**2) + ured - 2*muu4**2*ured - muu3**2*(1 + ured))))))*dconjg(AR)) + 
     &       2*xd23*((BR*(g1L*muu4*(g2R*(-3*muu1**4 + (-1 + muu2**2)**2 + 2*muu1**2*(1 + muu2**2))*muu3 - 
     &                   g2L*muu2*(-1 - muu1**2 + muu2**2)*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                g1R*muu1*(g2R*(1 + muu1**2 - muu2**2)*muu3*(-1 + muu1**2 + muu3**2 - ured) + 
     &                   g2L*muu2*(1 - muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - 2*muu4**4 + muu2**2*(-1 + muu3**2 - muu4**2) + 
     &                      muu1**2*(1 + muu3**2 - muu4**2) - 2*muu3**2*ured + 2*muu4**2*ured))) + 
     &             BL*(g1R*muu1*(1 + muu1**2 - muu2**2)*muu4*(2*g2L*muu2*muu3 + g2R*(muu2**2 + muu3**2 - ured)) + 
     &                g1L*(g2L*muu2*(1 + muu1**2 - muu2**2)*muu3*(muu1**2 + muu4**2 - ured) + 
     &                   g2R*(2*muu1**4*muu3**2 + (muu4**2 - ured)*(-((-1 + muu2**2)*muu3**2) + muu2**2*muu4**2 - ured) + 
     &                      muu1**2*(muu3**4 + muu2**2*(1 - 2*muu3**2) + (-2 + muu4**2)*ured - muu3**2*(muu4**2 + ured))))))*dconjg(AL) + 
     &          (BL*(g1R*muu4*(g2L*(-3*muu1**4 + (-1 + muu2**2)**2 + 2*muu1**2*(1 + muu2**2))*muu3 - 
     &                   g2R*muu2*(-1 - muu1**2 + muu2**2)*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                g1L*muu1*(g2L*(1 + muu1**2 - muu2**2)*muu3*(-1 + muu1**2 + muu3**2 - ured) + 
     &                   g2R*muu2*(1 - muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - 2*muu4**4 + muu2**2*(-1 + muu3**2 - muu4**2) + 
     &                      muu1**2*(1 + muu3**2 - muu4**2) - 2*muu3**2*ured + 2*muu4**2*ured))) + 
     &             BR*(g1L*muu1*(1 + muu1**2 - muu2**2)*muu4*(2*g2R*muu2*muu3 + g2L*(muu2**2 + muu3**2 - ured)) + 
     &                g1R*(g2R*muu2*(1 + muu1**2 - muu2**2)*muu3*(muu1**2 + muu4**2 - ured) + 
     &                   g2L*(2*muu1**4*muu3**2 + (muu4**2 - ured)*(-((-1 + muu2**2)*muu3**2) + muu2**2*muu4**2 - ured) + 
     &                      muu1**2*(muu3**4 + muu2**2*(1 - 2*muu3**2) + (-2 + muu4**2)*ured - muu3**2*(muu4**2 + ured))))))*dconjg(AR)) + 
     &       2*xd22*((-(BR*(g1R*muu1*(g2L*muu2*(-1 + muu3**4 - 2*muu3**2*muu4**2 + muu4**4) - g2R*muu3*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                  g1L*muu4*(g2R*(-1 + muu1**4 - 2*muu1**2*muu2**2 + muu2**4)*muu3 - g2L*muu2*(-1 + muu2**2 + muu4**2 - ured)))) + 
     &             BL*(g1R*muu1*muu4*(2*g2L*muu2*muu3 + g2R*(muu2**2 + muu3**2 - ured)) + 
     &                g1L*(g2L*muu2*muu3*(muu1**2 + muu4**2 - ured) + 
     &                   g2R*(muu1**4*muu3**2 + muu2**4*muu4**2 + (muu3**2 - ured)*(muu4**2 - ured) - 
     &                      muu2**2*(-muu4**4 + muu3**2*(muu4**2 - ured) + ured + muu4**2*ured) - 
     &                      muu1**2*(-muu3**4 + muu2**2*(-1 + muu3**2 + muu4**2) + ured - muu4**2*ured + muu3**2*(muu4**2 + ured))))))*dconjg(AL) + 
     &          (-(BL*(g1L*muu1*(g2R*muu2*(-1 + muu3**4 - 2*muu3**2*muu4**2 + muu4**4) - g2L*muu3*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                  g1R*muu4*(g2L*(-1 + muu1**4 - 2*muu1**2*muu2**2 + muu2**4)*muu3 - g2R*muu2*(-1 + muu2**2 + muu4**2 - ured)))) + 
     &             BR*(g1L*muu1*muu4*(2*g2R*muu2*muu3 + g2L*(muu2**2 + muu3**2 - ured)) + 
     &                g1R*(g2R*muu2*muu3*(muu1**2 + muu4**2 - ured) + 
     &                   g2L*(muu1**4*muu3**2 + muu2**4*muu4**2 + (muu3**2 - ured)*(muu4**2 - ured) - 
     &                      muu2**2*(-muu4**4 + muu3**2*(muu4**2 - ured) + ured + muu4**2*ured) - 
     &                      muu1**2*(-muu3**4 + muu2**2*(-1 + muu3**2 + muu4**2) + ured - muu4**2*ured + muu3**2*(muu4**2 + ured))))))*dconjg(AR))))/
     &  (4.*DenSmV1*Pi**2) 
	
	
	tempBoxZ = tempBoxZ + ChiChi2QQ_M2BoxGluon - tempBox2
	tempBox2 = ChiChi2QQ_M2BoxGluon

! Hboson contr.

	
	do i = imin,imax

	ChiChi2QQ_M2BoxGluon = ChiChi2QQ_M2BoxGluon -(CF*g3L*g4L*NC*sqrtS**6*(muu1*xd33*(dconjg(CR(i))*((g1L*
     &                 (2*g2L*muu1**2*muu2*(-1 + muu3**2 + muu4**2) + 
     &                   g2R*muu3*(-(muu1**2*(muu2**2 - muu3**2 + muu4**2)) - (-1 + muu2**2)*(muu4**2 - ured))) + 
     &                g1R*muu1*muu4*(2*g2R*(-1 + muu1**2 + muu2**2)*muu3 - g2L*muu2*(-1 + muu2**2 + muu4**2 - ured)))*dconjg(DL(i)) + 
     &             (g1L*muu4*(4*g2L*muu1**2*muu2*muu3 + g2R*(muu1**2*(1 + muu3**2 - muu4**2) - (-1 + muu2**2)*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                g1R*muu1*(g2R*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2) - g2L*muu2*muu3*(muu1**2 + muu4**2 - ured)))*dconjg(DR(i))) + 
     &          dconjg(CL(i))*((g1R*muu4*(4*g2R*muu1**2*muu2*muu3 + 
     &                   g2L*(muu1**2*(1 + muu3**2 - muu4**2) - (-1 + muu2**2)*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                g1L*muu1*(g2L*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2) - g2R*muu2*muu3*(muu1**2 + muu4**2 - ured)))*dconjg(DL(i)) + 
     &             (g1R*(2*g2R*muu1**2*muu2*(-1 + muu3**2 + muu4**2) + 
     &                   g2L*muu3*(-(muu1**2*(muu2**2 - muu3**2 + muu4**2)) - (-1 + muu2**2)*(muu4**2 - ured))) + 
     &                g1L*muu1*muu4*(2*g2L*(-1 + muu1**2 + muu2**2)*muu3 - g2R*muu2*(-1 + muu2**2 + muu4**2 - ured)))*dconjg(DR(i)))) + 
     &       (2*xd00*(dconjg(CR(i))*((g1L*muu1*(4*g2L*muu2*(-1 + muu3**2 + muu4**2) + g2R*muu3*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                  g1R*muu4*(4*g2R*(-1 + muu1**2 + muu2**2)*muu3 + g2L*muu2*(-1 + muu2**2 + muu4**2 - ured)))*dconjg(DL(i)) + 
     &               (g1L*muu1*muu4*(8*g2L*muu2*muu3 + g2R*(muu2**2 + muu3**2 - ured)) + 
     &                  g1R*(2*g2R*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2) + g2L*muu2*muu3*(muu1**2 + muu4**2 - ured)))*dconjg(DR(i))) + 
     &            dconjg(CL(i))*((g1R*muu1*muu4*(8*g2R*muu2*muu3 + g2L*(muu2**2 + muu3**2 - ured)) + 
     &                  g1L*(2*g2L*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2) + g2R*muu2*muu3*(muu1**2 + muu4**2 - ured)))*dconjg(DL(i)) + 
     &               (g1R*muu1*(4*g2R*muu2*(-1 + muu3**2 + muu4**2) + g2L*muu3*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                  g1L*muu4*(4*g2L*(-1 + muu1**2 + muu2**2)*muu3 + g2R*muu2*(-1 + muu2**2 + muu4**2 - ured)))*dconjg(DR(i)))))/sqrtS**2 + 
     &       xd23*(dconjg(CL(i))*((g1L*(2*g2R*muu1**2*muu2*muu3*(-1 + muu3**2 - muu4**2) + g2L*(muu1**4 - (-1 + muu2**2)**2)*(-1 + muu3**2 + muu4**2)) + 
     &                2*g1R*muu1*muu4*(2*g2R*muu2*(1 + muu1**2 - muu2**2)*muu3 + g2L*(-1 + muu2**2 + muu4**2 + muu1**2*(1 + muu3**2 - muu4**2) - ured)))*
     &              dconjg(DL(i)) + 2*(g1L*muu4*(g2L*(muu1**4 - (-1 + muu2**2)**2)*muu3 + g2R*muu1**2*muu2*(1 + muu3**2 - muu4**2)) + 
     &                g1R*muu1*(g2R*muu2*(1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + g2L*muu3*(muu4**2 + muu1**2*(muu3**2 - muu4**2) - ured)))*
     &              dconjg(DR(i))) + dconjg(CR(i))*(2*(g1R*muu4*(g2R*(muu1**4 - (-1 + muu2**2)**2)*muu3 + g2L*muu1**2*muu2*(1 + muu3**2 - muu4**2)) + 
     &                g1L*muu1*(g2L*muu2*(1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + g2R*muu3*(muu4**2 + muu1**2*(muu3**2 - muu4**2) - ured)))*
     &              dconjg(DL(i)) + (g1R*(2*g2L*muu1**2*muu2*muu3*(-1 + muu3**2 - muu4**2) + g2R*(muu1**4 - (-1 + muu2**2)**2)*(-1 + muu3**2 + muu4**2)) + 
     &                2*g1L*muu1*muu4*(2*g2L*muu2*(1 + muu1**2 - muu2**2)*muu3 + g2R*(-1 + muu2**2 + muu4**2 + muu1**2*(1 + muu3**2 - muu4**2) - ured)))*
     &              dconjg(DR(i)))) + xd13*(dconjg(CL(i))*((g1L*
     &                 (-4*g2R*muu1**2*muu2*muu3*muu4**2 + g2L*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2)*(muu1**2 + muu4**2 - ured)) + 
     &                2*g1R*muu1*muu4*(2*g2R*muu2*muu3*(muu1**2 + muu4**2 - ured) + 
     &                   g2L*(-1 + muu2**2 + muu3**2 - muu2**2*muu3**2 + muu4**2 - muu1**2*(-1 + muu4**2) - ured)))*dconjg(DL(i)) + 
     &             2*(g1L*muu4*(-(g2R*muu1**2*muu2*(-1 + muu3**2 + muu4**2)) + g2L*(-1 + muu1**2 + muu2**2)*muu3*(muu1**2 + muu4**2 - ured)) + 
     &                g1R*muu1*(-(g2L*(-1 + muu1**2 + muu2**2)*muu3*muu4**2) + g2R*muu2*(-1 + muu3**2 + muu4**2)*(muu1**2 + muu4**2 - ured)))*dconjg(DR(i))
     &             ) + dconjg(CR(i))*(2*(g1R*muu4*(-(g2L*muu1**2*muu2*(-1 + muu3**2 + muu4**2)) + 
     &                   g2R*(-1 + muu1**2 + muu2**2)*muu3*(muu1**2 + muu4**2 - ured)) + 
     &                g1L*muu1*(-(g2R*(-1 + muu1**2 + muu2**2)*muu3*muu4**2) + g2L*muu2*(-1 + muu3**2 + muu4**2)*(muu1**2 + muu4**2 - ured)))*dconjg(DL(i))
     &               + (g1R*(-4*g2L*muu1**2*muu2*muu3*muu4**2 + g2R*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2)*(muu1**2 + muu4**2 - ured)) + 
     &                2*g1L*muu1*muu4*(2*g2L*muu2*muu3*(muu1**2 + muu4**2 - ured) + 
     &                   g2R*(-1 + muu2**2 + muu3**2 - muu2**2*muu3**2 + muu4**2 - muu1**2*(-1 + muu4**2) - ured)))*dconjg(DR(i)))) + 
     &       xd22*(dconjg(CR(i))*((g1R*muu4*(2*g2R*(-1 + muu1**2 + muu2**2)*muu3 + 
     &                   g2L*muu2*(muu3**2 - muu2**2*muu3**2 + muu2**2*muu4**2 + muu1**2*(1 + muu3**2 - muu4**2) - ured)) + 
     &                g1L*muu1*(2*g2L*muu2*(-1 + muu3**2 + muu4**2) + 
     &                   g2R*muu3*(muu4**2 + muu1**2*(muu3**2 - muu4**2) + muu2**2*(1 - muu3**2 + muu4**2) - ured)))*dconjg(DL(i)) + 
     &             (g1L*muu1*muu4*(4*g2L*muu2*muu3 + g2R*(-1 + muu4**2 + muu1**2*(1 + muu3**2 - muu4**2) + muu2**2*(-muu3**2 + muu4**2) - ured)) + 
     &                g1R*(g2R*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   g2L*muu2*muu3*(-1 + (1 + muu1**2)*muu3**2 - muu1**2*muu4**2 + muu2**2*(1 - muu3**2 + muu4**2) - ured)))*dconjg(DR(i))) + 
     &          dconjg(CL(i))*((g1R*muu1*muu4*(4*g2R*muu2*muu3 + 
     &                   g2L*(-1 + muu4**2 + muu1**2*(1 + muu3**2 - muu4**2) + muu2**2*(-muu3**2 + muu4**2) - ured)) + 
     &                g1L*(g2L*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   g2R*muu2*muu3*(-1 + (1 + muu1**2)*muu3**2 - muu1**2*muu4**2 + muu2**2*(1 - muu3**2 + muu4**2) - ured)))*dconjg(DL(i)) + 
     &             (g1L*muu4*(2*g2L*(-1 + muu1**2 + muu2**2)*muu3 + 
     &                   g2R*muu2*(muu3**2 - muu2**2*muu3**2 + muu2**2*muu4**2 + muu1**2*(1 + muu3**2 - muu4**2) - ured)) + 
     &                g1R*muu1*(2*g2R*muu2*(-1 + muu3**2 + muu4**2) + 
     &                   g2L*muu3*(muu4**2 + muu1**2*(muu3**2 - muu4**2) + muu2**2*(1 - muu3**2 + muu4**2) - ured)))*dconjg(DR(i)))) + 
     &       xd12*(dconjg(CR(i))*(-2*(g1L*muu1*(g2R*(-1 + muu1**2 - muu2**2)*muu3*muu4**2 + g2L*muu2*(1 - 2*muu3**2 + muu3**4 - muu4**4)) + 
     &                g1R*muu4*(g2R*(-1 + muu1**2 + muu2**2)*muu3*(-1 + muu3**2 - muu4**2) + g2L*muu2*(-(muu2**2*muu4**2) + muu1**2*(-1 + muu4**2) + ured))
     &                )*dconjg(DL(i)) + (g1R*(2*g2L*muu2*(-1 - muu1**2 + muu2**2)*muu3*muu4**2 - 
     &                   g2R*(-1 + muu1**2 + muu2**2)*(1 - 2*muu3**2 + muu3**4 - muu4**4)) + 
     &                2*g1L*muu1*muu4*(2*g2L*muu2*muu3*(1 - muu3**2 + muu4**2) + 
     &                   g2R*(-1 + muu3**2 + muu4**2 + muu2**2*muu4**2 - muu1**2*(-1 + muu4**2) - ured)))*dconjg(DR(i))) + 
     &          dconjg(CL(i))*((g1L*(2*g2R*muu2*(-1 - muu1**2 + muu2**2)*muu3*muu4**2 - 
     &                   g2L*(-1 + muu1**2 + muu2**2)*(1 - 2*muu3**2 + muu3**4 - muu4**4)) + 
     &                2*g1R*muu1*muu4*(2*g2R*muu2*muu3*(1 - muu3**2 + muu4**2) + 
     &                   g2L*(-1 + muu3**2 + muu4**2 + muu2**2*muu4**2 - muu1**2*(-1 + muu4**2) - ured)))*dconjg(DL(i)) - 
     &             2*(g1R*muu1*(g2L*(-1 + muu1**2 - muu2**2)*muu3*muu4**2 + g2R*muu2*(1 - 2*muu3**2 + muu3**4 - muu4**4)) + 
     &                g1L*muu4*(g2L*(-1 + muu1**2 + muu2**2)*muu3*(-1 + muu3**2 - muu4**2) + g2R*muu2*(-(muu2**2*muu4**2) + muu1**2*(-1 + muu4**2) + ured))
     &                )*dconjg(DR(i)))) + muu4*xd0*(dconjg(CR(i))*
     &           ((g1L*muu1*muu4*(g2L*muu2*(-1 + 3*muu3**2 + muu4**2) + g2R*muu3*(2 - 2*muu1**2 + muu2**2 - muu3**2 + ured)) + 
     &                g1R*(g2R*muu3*(1 - muu3**2 - 2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) + muu1**2*(-1 + 2*muu4**2) + ured) + 
     &                   g2L*muu2*(1 - muu3**2 - muu4**2 - muu3**2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) - muu1**2*(-1 + 2*muu3**2 + muu4**2) + 
     &                      muu3**2*ured)))*dconjg(DL(i)) + (g1R*muu4*
     &                 (g2L*muu2*muu3*(-1 - 2*muu1**2 + muu2**2 - muu4**2 + ured) + 
     &                   g2R*(1 - 2*muu3**2 + muu2**2*(-1 + 2*muu3**2) - muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + ured)) + 
     &                g1L*muu1*(g2L*muu2*muu3*(-1 + muu3**2 + 3*muu4**2) + 
     &                   g2R*(-1 + 2*muu3**2 - muu3**4 + muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) - muu1**2*(-1 + 2*muu3**2 + muu4**2) + muu3**2*ured)))*
     &              dconjg(DR(i))) + dconjg(CL(i))*((g1L*muu4*(g2R*muu2*muu3*(-1 - 2*muu1**2 + muu2**2 - muu4**2 + ured) + 
     &                   g2L*(1 - 2*muu3**2 + muu2**2*(-1 + 2*muu3**2) - muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + ured)) + 
     &                g1R*muu1*(g2R*muu2*muu3*(-1 + muu3**2 + 3*muu4**2) + 
     &                   g2L*(-1 + 2*muu3**2 - muu3**4 + muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) - muu1**2*(-1 + 2*muu3**2 + muu4**2) + muu3**2*ured)))*
     &              dconjg(DL(i)) + (g1R*muu1*muu4*(g2R*muu2*(-1 + 3*muu3**2 + muu4**2) + g2L*muu3*(2 - 2*muu1**2 + muu2**2 - muu3**2 + ured)) + 
     &                g1L*(g2L*muu3*(1 - muu3**2 - 2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) + muu1**2*(-1 + 2*muu4**2) + ured) + 
     &                   g2R*muu2*(1 - muu3**2 - muu4**2 - muu3**2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) - muu1**2*(-1 + 2*muu3**2 + muu4**2) + 
     &                      muu3**2*ured)))*dconjg(DR(i)))) + xd2*
     &        (-(dconjg(CL(i))*((g1R*muu1*muu4*(2*g2R*muu2*muu3*(-2 + muu3**2 - muu4**2) + 
     &                     g2L*(2 - 3*muu3**2 - 2*muu4**2 + muu2**2*(1 - 2*muu3**2 - 2*muu4**2) + 2*muu1**2*(-1 + muu3**2 + muu4**2) + ured)) + 
     &                  g1L*(g2R*muu2*muu3*((3 - 4*muu2**2)*muu4**2 + muu1**2*(-1 + 4*muu4**2) + ured) + 
     &                     g2L*(-1 + 2*muu3**2 - muu3**4 + muu4**4 + muu1**2*(1 - muu4**4 + muu3**2*(-2 + muu4**2)) + 
     &                        muu2**2*(1 + muu3**4 - muu3**2*(2 + muu4**2)) + muu3**2*ured - muu4**2*ured)))*dconjg(DL(i)) + 
     &               (g1R*muu1*(g2R*muu2*(2 - 3*muu3**2 + muu3**4 - muu4**2 - muu4**4) - 
     &                     g2L*muu3*(-1 + muu3**2 + 4*muu4**2 + 4*muu2**2*muu4**2 + muu1**2*(1 - 4*muu4**2) - ured)) + 
     &                  g1L*muu4*(g2L*(-1 + muu1**2 + muu2**2)*muu3*(-2 + muu3**2 - muu4**2) + 
     &                     g2R*muu2*(-1 + 2*muu3**2 + muu4**2 + muu2**2*(1 - 2*muu3**2 - 2*muu4**2) + 2*muu1**2*(-1 + muu3**2 + muu4**2) + ured)))*
     &                dconjg(DR(i)))) - dconjg(CR(i))*((g1L*muu1*
     &                 (g2L*muu2*(2 - 3*muu3**2 + muu3**4 - muu4**2 - muu4**4) - 
     &                   g2R*muu3*(-1 + muu3**2 + 4*muu4**2 + 4*muu2**2*muu4**2 + muu1**2*(1 - 4*muu4**2) - ured)) + 
     &                g1R*muu4*(g2R*(-1 + muu1**2 + muu2**2)*muu3*(-2 + muu3**2 - muu4**2) + 
     &                   g2L*muu2*(-1 + 2*muu3**2 + muu4**2 + muu2**2*(1 - 2*muu3**2 - 2*muu4**2) + 2*muu1**2*(-1 + muu3**2 + muu4**2) + ured)))*
     &              dconjg(DL(i)) + (g1L*muu1*muu4*(2*g2L*muu2*muu3*(-2 + muu3**2 - muu4**2) + 
     &                   g2R*(2 - 3*muu3**2 - 2*muu4**2 + muu2**2*(1 - 2*muu3**2 - 2*muu4**2) + 2*muu1**2*(-1 + muu3**2 + muu4**2) + ured)) + 
     &                g1R*(g2L*muu2*muu3*((3 - 4*muu2**2)*muu4**2 + muu1**2*(-1 + 4*muu4**2) + ured) + 
     &                   g2R*(-1 + 2*muu3**2 - muu3**4 + muu4**4 + muu1**2*(1 - muu4**4 + muu3**2*(-2 + muu4**2)) + 
     &                      muu2**2*(1 + muu3**4 - muu3**2*(2 + muu4**2)) + muu3**2*ured - muu4**2*ured)))*dconjg(DR(i)))) + 
     &       muu4*xd11*(dconjg(CL(i))*((g1L*muu4*(g2L*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2) - g2R*muu2*muu3*(muu1**2 + muu4**2 - ured)) + 
     &                g1R*muu1*(4*g2R*muu2*muu3*muu4**2 + g2L*(-1 - muu3**4 + muu4**2 + muu2**2*muu4**2 - muu1**2*(-1 + muu3**2 + muu4**2) - ured + 
     &                      muu3**2*(2 + ured))))*dconjg(DL(i)) + 
     &             (g1R*muu1*muu4*(2*g2R*muu2*(-1 + muu3**2 + muu4**2) - g2L*muu3*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                g1L*(2*g2L*(-1 + muu1**2 + muu2**2)*muu3*muu4**2 - 
     &                   g2R*muu2*(-(muu2**2*muu4**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + ured - muu3**2*ured)))*dconjg(DR(i))) + 
     &          dconjg(CR(i))*((g1L*muu1*muu4*(2*g2L*muu2*(-1 + muu3**2 + muu4**2) - g2R*muu3*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                g1R*(2*g2R*(-1 + muu1**2 + muu2**2)*muu3*muu4**2 - 
     &                   g2L*muu2*(-(muu2**2*muu4**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + ured - muu3**2*ured)))*dconjg(DL(i)) + 
     &             (g1R*muu4*(g2R*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2) - g2L*muu2*muu3*(muu1**2 + muu4**2 - ured)) + 
     &                g1L*muu1*(4*g2L*muu2*muu3*muu4**2 + g2R*(-1 - muu3**4 + muu4**2 + muu2**2*muu4**2 - muu1**2*(-1 + muu3**2 + muu4**2) - ured + 
     &                      muu3**2*(2 + ured))))*dconjg(DR(i)))) + 
     &       xd1*(dconjg(CL(i))*((g1L*(g2R*muu2*muu3*(muu2**2*muu4**2 - muu3**2*muu4**2 - 3*muu4**4 - muu1**2*(-1 + muu3**2 + 4*muu4**2) - ured + 
     &                      muu3**2*ured + 3*muu4**2*ured) + g2L*
     &                    (-2*muu4**4 + muu1**2*(-muu3**4 + 2*muu4**2*(-1 + muu4**2) + muu3**2*(1 + muu4**2)) + (-1 + muu3**2)*ured + 
     &                      muu4**2*(2 - 3*muu3**2 + muu2**2*(-1 + 2*muu3**2) + 2*ured))) + 
     &                g1R*muu1*muu4*(g2R*muu2*muu3*(1 - muu3**2 + 5*muu4**2) + 
     &                   g2L*(-2 - 3*muu3**4 + 2*muu4**2 + muu2**2*muu4**2 + muu1**2*(2 - 3*muu3**2 - 2*muu4**2) - 2*ured + muu4**2*ured + 
     &                      muu3**2*(5 - muu4**2 + 3*ured))))*dconjg(DL(i)) + 
     &             (-(g1R*muu1*(g2R*muu2*(1 + muu3**4 + muu4**2 - 2*muu4**4 - muu3**2*(2 + 3*muu4**2)) + 
     &                     g2L*muu3*(1 + muu3**4 - 4*muu4**2 - muu2**2*muu4**2 + muu1**2*(-1 + muu3**2 + 4*muu4**2) + muu3**2*(-2 + 3*muu4**2 - ured) + 
     &                        ured - 3*muu4**2*ured))) + g1L*muu4*
     &                 (g2L*muu3*((-3 + 2*muu2**2)*muu4**2 - muu1**2*(muu3**2 - 3*muu4**2) + ured) + 
     &                   g2R*muu2*(-muu4**4 + muu1**2*(2 - 3*muu3**2 - 2*muu4**2) + (-2 + 3*muu3**2)*ured + muu4**2*(1 + muu2**2 - 3*muu3**2 + ured))))*
     &              dconjg(DR(i))) + dconjg(CR(i))*((-(g1L*muu1*
     &                   (g2L*muu2*(1 + muu3**4 + muu4**2 - 2*muu4**4 - muu3**2*(2 + 3*muu4**2)) + 
     &                     g2R*muu3*(1 + muu3**4 - 4*muu4**2 - muu2**2*muu4**2 + muu1**2*(-1 + muu3**2 + 4*muu4**2) + muu3**2*(-2 + 3*muu4**2 - ured) + 
     &                        ured - 3*muu4**2*ured))) + g1R*muu4*
     &                 (g2R*muu3*((-3 + 2*muu2**2)*muu4**2 - muu1**2*(muu3**2 - 3*muu4**2) + ured) + 
     &                   g2L*muu2*(-muu4**4 + muu1**2*(2 - 3*muu3**2 - 2*muu4**2) + (-2 + 3*muu3**2)*ured + muu4**2*(1 + muu2**2 - 3*muu3**2 + ured))))*
     &              dconjg(DL(i)) + (g1R*(g2L*muu2*muu3*(muu2**2*muu4**2 - muu3**2*muu4**2 - 3*muu4**4 - muu1**2*(-1 + muu3**2 + 4*muu4**2) - ured + 
     &                      muu3**2*ured + 3*muu4**2*ured) + g2R*
     &                    (-2*muu4**4 + muu1**2*(-muu3**4 + 2*muu4**2*(-1 + muu4**2) + muu3**2*(1 + muu4**2)) + (-1 + muu3**2)*ured + 
     &                      muu4**2*(2 - 3*muu3**2 + muu2**2*(-1 + 2*muu3**2) + 2*ured))) + 
     &                g1L*muu1*muu4*(g2L*muu2*muu3*(1 - muu3**2 + 5*muu4**2) + 
     &                   g2R*(-2 - 3*muu3**4 + 2*muu4**2 + muu2**2*muu4**2 + muu1**2*(2 - 3*muu3**2 - 2*muu4**2) - 2*ured + muu4**2*ured + 
     &                      muu3**2*(5 - muu4**2 + 3*ured))))*dconjg(DR(i)))) + 
     &       xd3*(dconjg(CL(i))*((-(g1R*muu1*muu4*(g2L*(2 - 3*muu3**2 - 2*muu4**2 + 2*muu1**2*(-1 + muu3**2 + muu4**2) + 
     &                        muu2**2*(-2 + 3*muu3**2 + muu4**2) + ured) + g2R*muu2*muu3*(-1 - 3*muu1**2 + muu2**2 - 2*muu4**2 + 2*ured))) + 
     &                g1L*(-(g2R*muu1**2*muu2*muu3*(-1 + muu3**2 + 7*muu4**2)) + 
     &                   g2L*(muu1**4*(-1 + 2*muu3**2 + muu4**2) + (-1 + muu2**2)*muu4**2*(-1 + muu2**2 + muu4**2 - ured) + 
     &                      muu1**2*(muu3**4 + (-1 + muu4**2)**2 - muu2**2*(-1 + muu3**2 + muu4**2) - muu3**2*(2 + muu4**2 + ured)))))*dconjg(DL(i)) + 
     &             (g1L*muu4*(g2R*muu1**2*muu2*(3 - 5*muu3**2 - 3*muu4**2) + 
     &                   g2L*muu3*(2*muu1**4 + muu1**2*(-2 + muu4**2 - ured) + (-1 + muu2**2)*(muu4**2 - ured))) + 
     &                g1R*muu1*(g2L*muu3*(-1 + muu3**2 + 4*muu4**2 + muu1**2*(1 - 4*muu4**2) - muu2**2*(-1 + muu3**2 + 3*muu4**2) - ured) + 
     &                   g2R*muu2*(-1 + muu3**2 - muu2**2*(-1 + muu3**2) + muu3**2*muu4**2 + muu4**4 + muu1**2*(-1 + 2*muu3**2 + muu4**2) - muu3**2*ured - 
     &                      muu4**2*ured)))*dconjg(DR(i))) + dconjg(CR(i))*
     &           ((g1R*muu4*(g2L*muu1**2*muu2*(3 - 5*muu3**2 - 3*muu4**2) + 
     &                   g2R*muu3*(2*muu1**4 + muu1**2*(-2 + muu4**2 - ured) + (-1 + muu2**2)*(muu4**2 - ured))) + 
     &                g1L*muu1*(g2R*muu3*(-1 + muu3**2 + 4*muu4**2 + muu1**2*(1 - 4*muu4**2) - muu2**2*(-1 + muu3**2 + 3*muu4**2) - ured) + 
     &                   g2L*muu2*(-1 + muu3**2 - muu2**2*(-1 + muu3**2) + muu3**2*muu4**2 + muu4**4 + muu1**2*(-1 + 2*muu3**2 + muu4**2) - muu3**2*ured - 
     &                      muu4**2*ured)))*dconjg(DL(i)) + (-(g1L*muu1*muu4*
     &                   (g2R*(2 - 3*muu3**2 - 2*muu4**2 + 2*muu1**2*(-1 + muu3**2 + muu4**2) + muu2**2*(-2 + 3*muu3**2 + muu4**2) + ured) + 
     &                     g2L*muu2*muu3*(-1 - 3*muu1**2 + muu2**2 - 2*muu4**2 + 2*ured))) + 
     &                g1R*(-(g2L*muu1**2*muu2*muu3*(-1 + muu3**2 + 7*muu4**2)) + 
     &                   g2R*(muu1**4*(-1 + 2*muu3**2 + muu4**2) + (-1 + muu2**2)*muu4**2*(-1 + muu2**2 + muu4**2 - ured) + 
     &                      muu1**2*(muu3**4 + (-1 + muu4**2)**2 - muu2**2*(-1 + muu3**2 + muu4**2) - muu3**2*(2 + muu4**2 + ured)))))*dconjg(DR(i))))))/
     &  (4.*Pi**2*DenSmH1(i)) 
	
	enddo

	tempBoxH = tempBoxH + ChiChi2QQ_M2BoxGluon - tempBox2
	tempBox2 = ChiChi2QQ_M2BoxGluon

	do i = 1,2

! U-channel contr.

	 ChiChi2QQ_M2BoxGluon = ChiChi2QQ_M2BoxGluon + (CF*g3L*g4L*NC*sqrtS**4*(muu4**2*xd11*(-(g1L*(g2L*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                  (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*((1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                     dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g2R*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*(2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                    (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                    (1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*dconjg(HR(indj,i,tt,gen)))))) - 
     &         g1R*(g2L*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*(-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                  2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*((1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                   dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*(2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                  (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                  (1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*dconjg(HR(indj,i,tt,gen)))))) + 
     &      muu1**2*xd33*(-(g1L*(g2L*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                  (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*((1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                     dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g2R*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*(2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                    (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                    (1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*dconjg(HR(indj,i,tt,gen)))))) - 
     &         g1R*(g2L*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*(-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                  2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*((1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                   dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*(2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                  (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                  (1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*dconjg(HR(indj,i,tt,gen)))))) + 
     &      (-1 + muu2**2 + muu3**2 - tred)*xd13*(g1L*(g2L*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*((1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                   dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*(2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                  (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                  (1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &         g1R*(g2L*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*(-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                  2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*((1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                   dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*(2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                  (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                  (1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*dconjg(HR(indj,i,tt,gen)))))) - 
     &      (2*xd00*(g1L*(-(g2R*((-1 + muu2**2 + muu3**2 - tred)*dconjg(GL(indi,i,tt,gen)) - 2*muu1*muu4*dconjg(GR(indi,i,tt,gen)))*
     &                 (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &              2*g2L*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                  (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*((1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                     dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &           g1R*(-(g2L*(-2*muu1*muu4*dconjg(GL(indi,i,tt,gen)) + (-1 + muu2**2 + muu3**2 - tred)*dconjg(GR(indi,i,tt,gen)))*
     &                 (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen)))) + 
     &              2*g2R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                  (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                    (1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*dconjg(HR(indj,i,tt,gen)))))))/sqrtS**2 + 
     &      muu4*xd0*(g1L*(g2L*(muu1*(-1 + muu3**2)*dconjg(GR(indi,i,tt,gen))*
     &                ((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               muu4*dconjg(GL(indi,i,tt,gen))*((-muu1**4 + muu3**2 - muu2**2*tred + muu1**2*(1 + muu2**2 - muu3**2 - muu4**2 + tred))*
     &                   dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(1 + 2*muu1**2 - muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(muu1*(-1 + muu3**2)*dconjg(GR(indi,i,tt,gen))*
     &                (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               muu4*dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(1 + 2*muu1**2 - muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  (-muu1**4 + muu3**2 - muu2**2*tred + muu1**2*(1 + muu2**2 - muu3**2 - muu4**2 + tred))*dconjg(HR(indj,i,tt,gen))))) + 
     &         g1R*(g2L*(muu1*(-1 + muu3**2)*dconjg(GL(indi,i,tt,gen))*
     &                ((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               muu4*dconjg(GR(indi,i,tt,gen))*((-muu1**4 + muu3**2 - muu2**2*tred + muu1**2*(1 + muu2**2 - muu3**2 - muu4**2 + tred))*
     &                   dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(1 + 2*muu1**2 - muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(muu1*(-1 + muu3**2)*dconjg(GL(indi,i,tt,gen))*
     &                (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               muu4*dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(1 + 2*muu1**2 - muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  (-muu1**4 + muu3**2 - muu2**2*tred + muu1**2*(1 + muu2**2 - muu3**2 - muu4**2 + tred))*dconjg(HR(indj,i,tt,gen)))))) + 
     &      muu1*xd3*(g1R*(g2L*(-(muu1*(-1 + muu3**2 + muu4**2)*dconjg(GR(indi,i,tt,gen))*
     &                  ((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*dconjg(HR(indj,i,tt,gen)))) + 
     &               muu4*dconjg(GL(indi,i,tt,gen))*((1 + muu1**4 - 2*muu3**2 - muu4**2 + tred - muu1**2*(2 + muu2**2 - 2*muu3**2 - muu4**2 + tred) + 
     &                     muu2**2*(-1 + muu4**2 + tred))*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(-(muu1*(-1 + muu3**2 + muu4**2)*dconjg(GR(indi,i,tt,gen))*
     &                  (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &               muu4*dconjg(GL(indi,i,tt,gen))*(-2*muu2*muu3*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  (1 + muu1**4 - 2*muu3**2 - muu4**2 + tred - muu1**2*(2 + muu2**2 - 2*muu3**2 - muu4**2 + tred) + muu2**2*(-1 + muu4**2 + tred))*
     &                   dconjg(HR(indj,i,tt,gen))))) + g1L*(g2L*
     &             (-(muu1*(-1 + muu3**2 + muu4**2)*dconjg(GL(indi,i,tt,gen))*
     &                  ((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*dconjg(HR(indj,i,tt,gen)))) + 
     &               muu4*dconjg(GR(indi,i,tt,gen))*((1 + muu1**4 - 2*muu3**2 - muu4**2 + tred - muu1**2*(2 + muu2**2 - 2*muu3**2 - muu4**2 + tred) + 
     &                     muu2**2*(-1 + muu4**2 + tred))*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(-(muu1*(-1 + muu3**2 + muu4**2)*dconjg(GL(indi,i,tt,gen))*
     &                  (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &               muu4*dconjg(GR(indi,i,tt,gen))*(-2*muu2*muu3*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  (1 + muu1**4 - 2*muu3**2 - muu4**2 + tred - muu1**2*(2 + muu2**2 - 2*muu3**2 - muu4**2 + tred) + muu2**2*(-1 + muu4**2 + tred))*
     &                   dconjg(HR(indj,i,tt,gen)))))) + xd12*(-(g1R*
     &            (g2L*(-(muu1*muu4*(-1 + muu3**2 - muu4**2)*dconjg(GL(indi,i,tt,gen))*
     &                    (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen)))) + 
     &                 dconjg(GR(indi,i,tt,gen))*((1 - 2*muu3**2 + muu2**4*muu3**2 + muu3**4 + muu1**4*muu4**2 + muu3**2*muu4**2 - muu4**4 + 
     &                       muu2**2*(1 + muu3**2 - muu4**2)*(-1 + muu3**2 - muu4**2 - tred) + 2*tred - 2*muu3**2*tred + tred**2 - 
     &                       muu1**2*(1 + muu3**4 + 2*muu4**2 - muu4**4 + muu2**2*(-1 + muu3**2 + muu4**2) + tred + muu4**2*tred - muu3**2*(2 + tred)))*
     &                     dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(-1 + muu3**2 - muu4**2)*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g2R*(-1 + muu3**2 - muu4**2)*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                  (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                    (-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))))) + 
     &         g1L*(g2L*(-1 + muu3**2 - muu4**2)*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*((1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                   dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) - 
     &            g2R*(-(muu1*muu4*(-1 + muu3**2 - muu4**2)*dconjg(GR(indi,i,tt,gen))*
     &                  (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &               dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(-1 + muu3**2 - muu4**2)*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  (1 - 2*muu3**2 + muu2**4*muu3**2 + muu3**4 + muu1**4*muu4**2 + muu3**2*muu4**2 - muu4**4 + 
     &                     muu2**2*(1 + muu3**2 - muu4**2)*(-1 + muu3**2 - muu4**2 - tred) + 2*tred - 2*muu3**2*tred + tred**2 - 
     &                     muu1**2*(1 + muu3**4 + 2*muu4**2 - muu4**4 + muu2**2*(-1 + muu3**2 + muu4**2) + tred + muu4**2*tred - muu3**2*(2 + tred)))*
     &                   dconjg(HR(indj,i,tt,gen)))))) + xd2*(g1L*
     &          (g2L*(2*muu1*muu3*muu4*dconjg(GR(indi,i,tt,gen))*
     &                ((-1 + muu1**2 - muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) - muu2*(1 + muu3**2 - muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*((-1 + muu3**2 + muu2**4*muu3**2 + muu4**2 - muu1**4*muu4**2 - tred + 
     &                     muu2**2*(1 + muu3**4 + muu3**2*(-1 + muu4**2 - tred) - muu4**2*(1 + tred)) + 
     &                     muu1**2*(1 - muu4**4 - muu3**2*(1 + muu2**2 + muu4**2 - tred) + muu4**2*(1 + muu2**2 + tred)))*dconjg(HL(indj,i,tt,gen)) + 
     &                  2*muu2*muu3*(-1 + muu3**2 + muu4**2 + muu1**2*muu4**2 - muu2**2*(-1 + muu4**2) - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(-2*muu1*muu3*muu4*dconjg(GR(indi,i,tt,gen))*
     &                (muu2*(1 + muu3**2 - muu4**2)*dconjg(HL(indj,i,tt,gen)) + (1 - muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*(2*muu2*muu3*(-1 + muu3**2 + muu4**2 + muu1**2*muu4**2 - muu2**2*(-1 + muu4**2) - tred)*
     &                   dconjg(HL(indj,i,tt,gen)) + (-2*muu1**4*muu4**2 + (-1 + muu3**2 - tred)*(1 - muu4**2 + tred) + 
     &                     muu2**2*(1 + muu4**2*(-1 + 2*muu3**2 - 2*tred) + tred) + 
     &                     muu1**2*(1 - muu3**2 + 2*muu4**2 - 2*muu4**4 + muu2**2*(-1 + 2*muu4**2) + tred + 2*muu4**2*tred))*dconjg(HR(indj,i,tt,gen)))))
     &          + g1R*(g2L*(2*muu1*muu3*muu4*dconjg(GL(indi,i,tt,gen))*
     &                ((-1 + muu1**2 - muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) - muu2*(1 + muu3**2 - muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*((-2*muu1**4*muu4**2 + (-1 + muu3**2 - tred)*(1 - muu4**2 + tred) + 
     &                     muu2**2*(1 + muu4**2*(-1 + 2*muu3**2 - 2*tred) + tred) + 
     &                     muu1**2*(1 - muu3**2 + 2*muu4**2 - 2*muu4**4 + muu2**2*(-1 + 2*muu4**2) + tred + 2*muu4**2*tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                  2*muu2*muu3*(-1 + muu3**2 + muu4**2 + muu1**2*muu4**2 - muu2**2*(-1 + muu4**2) - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(-2*muu1*muu3*muu4*dconjg(GL(indi,i,tt,gen))*
     &                (muu2*(1 + muu3**2 - muu4**2)*dconjg(HL(indj,i,tt,gen)) + (1 - muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*(2*muu2*muu3*(-1 + muu3**2 + muu4**2 + muu1**2*muu4**2 - muu2**2*(-1 + muu4**2) - tred)*
     &                   dconjg(HL(indj,i,tt,gen)) + (-1 + muu3**2 + muu2**4*muu3**2 + muu4**2 - muu1**4*muu4**2 - tred + 
     &                     muu2**2*(1 + muu3**4 + muu3**2*(-1 + muu4**2 - tred) - muu4**2*(1 + tred)) + 
     &                     muu1**2*(1 - muu4**4 - muu3**2*(1 + muu2**2 + muu4**2 - tred) + muu4**2*(1 + muu2**2 + tred)))*dconjg(HR(indj,i,tt,gen)))))) + 
     &      xd23*(-(g1R*(g2L*(-(muu1*(1 + muu1**2 - muu2**2)*muu4*dconjg(GL(indi,i,tt,gen))*
     &                    ((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*dconjg(HR(indj,i,tt,gen)))) + 
     &                 dconjg(GR(indi,i,tt,gen))*((muu1**4*(-1 + muu3**2 + muu4**2) + 
     &                       (-1 + muu2**2 + muu3**2 - tred)*(-1 + muu4**2 + muu2**2*(1 + muu3**2 - muu4**2) - tred) + 
     &                       muu1**2*(muu2**2*(1 - 2*muu3**2) + muu4**4 + muu3**2*tred - muu4**2*(2 + muu3**2 + tred)))*dconjg(HL(indj,i,tt,gen)) + 
     &                    muu2*(-1 - muu1**2 + muu2**2)*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g2R*(1 + muu1**2 - muu2**2)*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                  (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                    (1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*dconjg(HR(indj,i,tt,gen)))))) - 
     &         g1L*(g2L*(1 + muu1**2 - muu2**2)*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*((1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                   dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(-(muu1*(1 + muu1**2 - muu2**2)*muu4*dconjg(GR(indi,i,tt,gen))*
     &                  (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &               dconjg(GL(indi,i,tt,gen))*(muu2*(-1 - muu1**2 + muu2**2)*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  (muu1**4*(-1 + muu3**2 + muu4**2) + (-1 + muu2**2 + muu3**2 - tred)*(-1 + muu4**2 + muu2**2*(1 + muu3**2 - muu4**2) - tred) + 
     &                     muu1**2*(muu2**2*(1 - 2*muu3**2) + muu4**4 + muu3**2*tred - muu4**2*(2 + muu3**2 + tred)))*dconjg(HR(indj,i,tt,gen)))))) + 
     &      xd22*(g1R*(g2L*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                ((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*((-(muu2**4*muu3**2) - muu1**4*muu4**2 + (-1 + muu3**2 - tred)*(1 - muu4**2 + tred) + 
     &                     muu1**2*(1 + muu4**2 - muu4**4 + muu2**2*(-1 + muu3**2 + muu4**2) + muu3**2*(-1 + muu4**2 - tred) + tred + muu4**2*tred) + 
     &                     muu2**2*(-muu3**4 - (-1 + muu4**2)*(1 + tred) + muu3**2*(1 + muu4**2 + tred)))*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*(-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + 
     &                  (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  (-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &         g1L*(g2L*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 
     &                  2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*
     &                   dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*(-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + 
     &                  (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  (-(muu2**4*muu3**2) - muu1**4*muu4**2 + (-1 + muu3**2 - tred)*(1 - muu4**2 + tred) + 
     &                     muu1**2*(1 + muu4**2 - muu4**4 + muu2**2*(-1 + muu3**2 + muu4**2) + muu3**2*(-1 + muu4**2 - tred) + tred + muu4**2*tred) + 
     &                     muu2**2*(-muu3**4 - (-1 + muu4**2)*(1 + tred) + muu3**2*(1 + muu4**2 + tred)))*dconjg(HR(indj,i,tt,gen)))))) + 
     &      xd1*(-(g1L*(g2L*(muu1*muu4*(-1 + muu3**2 + muu4**2)*dconjg(GR(indi,i,tt,gen))*
     &                  ((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*((1 - 2*muu3**2 + muu3**4 - muu4**2 + muu1**4*muu4**2 + muu3**2*muu4**2 - muu3**4*muu4**2 - 
     &                       muu3**2*muu4**4 + 2*tred - 3*muu3**2*tred + muu3**4*tred - 2*muu4**2*tred + 2*muu3**2*muu4**2*tred + muu4**4*tred + tred**2 - 
     &                       muu3**2*tred**2 - muu4**2*tred**2 - 
     &                       muu1**2*(1 + muu3**4 + muu4**2 - muu4**4 + muu2**2*(-1 + muu3**2 + 2*muu4**2) + tred - muu3**2*(2 + muu4**2 + tred)) + 
     &                       muu2**2*(-1 - tred + muu3**2*(1 - muu4**2 + tred) + muu4**2*(1 + 2*tred)))*dconjg(HL(indj,i,tt,gen)) + 
     &                    2*muu2*muu3*(1 + muu3**4 - muu4**2 - muu1**2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) + tred - muu3**2*(2 + tred))*
     &                     dconjg(HR(indj,i,tt,gen)))) + g2R*(muu1*muu4*(-1 + muu3**2 + muu4**2)*dconjg(GR(indi,i,tt,gen))*
     &                  (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(2*muu2*muu3*(1 + muu3**4 - muu4**2 - muu1**2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) + tred - 
     &                       muu3**2*(2 + tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                    (1 - 2*muu3**2 + muu3**4 - muu4**2 + muu1**4*muu4**2 + muu3**2*muu4**2 - muu3**4*muu4**2 - muu3**2*muu4**4 + 2*tred - 
     &                       3*muu3**2*tred + muu3**4*tred - 2*muu4**2*tred + 2*muu3**2*muu4**2*tred + muu4**4*tred + tred**2 - muu3**2*tred**2 - 
     &                       muu4**2*tred**2 - muu1**2*(1 + muu3**4 + muu4**2 - muu4**4 + muu2**2*(-1 + muu3**2 + 2*muu4**2) + tred - 
     &                          muu3**2*(2 + muu4**2 + tred)) + muu2**2*(-1 - tred + muu3**2*(1 - muu4**2 + tred) + muu4**2*(1 + 2*tred)))*
     &                     dconjg(HR(indj,i,tt,gen)))))) - g1R*
     &          (g2L*(muu1*muu4*(-1 + muu3**2 + muu4**2)*dconjg(GL(indi,i,tt,gen))*
     &                ((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*((1 - 2*muu3**2 + muu3**4 - muu4**2 + muu1**4*muu4**2 + muu3**2*muu4**2 - muu3**4*muu4**2 - 
     &                     muu3**2*muu4**4 + 2*tred - 3*muu3**2*tred + muu3**4*tred - 2*muu4**2*tred + 2*muu3**2*muu4**2*tred + muu4**4*tred + tred**2 - 
     &                     muu3**2*tred**2 - muu4**2*tred**2 - 
     &                     muu1**2*(1 + muu3**4 + muu4**2 - muu4**4 + muu2**2*(-1 + muu3**2 + 2*muu4**2) + tred - muu3**2*(2 + muu4**2 + tred)) + 
     &                     muu2**2*(-1 - tred + muu3**2*(1 - muu4**2 + tred) + muu4**2*(1 + 2*tred)))*dconjg(HL(indj,i,tt,gen)) + 
     &                  2*muu2*muu3*(1 + muu3**4 - muu4**2 - muu1**2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) + tred - muu3**2*(2 + tred))*
     &                   dconjg(HR(indj,i,tt,gen)))) + g2R*(muu1*muu4*(-1 + muu3**2 + muu4**2)*dconjg(GL(indi,i,tt,gen))*
     &                (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*(2*muu2*muu3*(1 + muu3**4 - muu4**2 - muu1**2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) + tred - 
     &                     muu3**2*(2 + tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                  (1 - 2*muu3**2 + muu3**4 - muu4**2 + muu1**4*muu4**2 + muu3**2*muu4**2 - muu3**4*muu4**2 - muu3**2*muu4**4 + 2*tred - 
     &                     3*muu3**2*tred + muu3**4*tred - 2*muu4**2*tred + 2*muu3**2*muu4**2*tred + muu4**4*tred + tred**2 - muu3**2*tred**2 - 
     &                     muu4**2*tred**2 - muu1**2*(1 + muu3**4 + muu4**2 - muu4**4 + muu2**2*(-1 + muu3**2 + 2*muu4**2) + tred - 
     &                        muu3**2*(2 + muu4**2 + tred)) + muu2**2*(-1 - tred + muu3**2*(1 - muu4**2 + tred) + muu4**2*(1 + 2*tred)))*
     &                   dconjg(HR(indj,i,tt,gen))))))))/(4.*Pi**2*(-ured + MSf(i,itypeU,gen)**2/sqrtS**2)) 

	
	tempBoxu = tempBoxu + ChiChi2QQ_M2BoxGluon - tempBox2
	tempBox2 = ChiChi2QQ_M2BoxGluon
	
! T-channel contr.

	ChiChi2QQ_M2BoxGluon = ChiChi2QQ_M2BoxGluon + (CF*g3L*g4L*NC*sqrtS**4*(xd22*(-(dconjg(ELc(indi,i,tt,gen))*
     &            ((-(g1R*muu1*muu4*(-2*g2R*muu2*muu3 - g2L*(-1 + muu1**2 - muu2**2)*(1 + muu3**2 - muu4**2))) + 
     &                 g1L*(-(g2R*muu2*(-1 - muu1**2 + muu2**2)*muu3*(-1 + muu3**2 - muu4**2)) + 
     &                    g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &              (g1L*muu4*(g2R*(-1 + muu1**4 - 2*muu1**2*muu2**2 + muu2**4)*muu3 + g2L*muu2*(muu1**2 + muu3**2 - tred)) + 
     &                 g1R*muu1*(g2L*muu2*(-1 + muu3**4 - 2*muu3**2*muu4**2 + muu4**4) + g2R*muu3*(muu2**2 + muu4**2 - tred)))*
     &               dconjg(FR(indj,i,tt,gen)))) - dconjg(ER(indi,i,tt,gen))*
     &          ((g1R*muu4*(g2L*(-1 + muu1**4 - 2*muu1**2*muu2**2 + muu2**4)*muu3 + g2R*muu2*(muu1**2 + muu3**2 - tred)) + 
     &               g1L*muu1*(g2R*muu2*(-1 + muu3**4 - 2*muu3**2*muu4**2 + muu4**4) + g2L*muu3*(muu2**2 + muu4**2 - tred)))*
     &             dconjg(FL(indj,i,tt,gen)) + (-(g1L*muu1*muu4*
     &                  (-2*g2L*muu2*muu3 - g2R*(-1 + muu1**2 - muu2**2)*(1 + muu3**2 - muu4**2))) + 
     &               g1R*(-(g2L*muu2*(-1 - muu1**2 + muu2**2)*muu3*(-1 + muu3**2 - muu4**2)) + 
     &                  g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FR(indj,i,tt,gen)))) - 
     &      muu1*xd33*(dconjg(ELc(indi,i,tt,gen))*((g1R*muu4*
     &                (2*g2R*muu1**2*muu2*muu3 + g2L*(-1 + muu1**2 + muu2**2)*(muu1**2 + muu3**2 - tred)) + 
     &               g1L*muu1*(2*g2R*muu2*muu3*(-1 + muu2**2 + muu3**2 - tred) + g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*
     &             dconjg(FL(indj,i,tt,gen)) + (g1L*muu1*muu4*
     &                (2*g2R*(-1 + muu1**2 + muu2**2)*muu3 + g2L*muu2*(muu1**2 + muu3**2 - tred)) + 
     &               g1R*(g2L*muu2*(muu1**2 + muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + 
     &                  g2R*muu1**2*muu3*(muu2**2 + muu4**2 - tred)))*dconjg(FR(indj,i,tt,gen))) + 
     &         dconjg(ER(indi,i,tt,gen))*((g1R*muu1*muu4*
     &                (2*g2L*(-1 + muu1**2 + muu2**2)*muu3 + g2R*muu2*(muu1**2 + muu3**2 - tred)) + 
     &               g1L*(g2R*muu2*(muu1**2 + muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + 
     &                  g2L*muu1**2*muu3*(muu2**2 + muu4**2 - tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1L*muu4*(2*g2L*muu1**2*muu2*muu3 + g2R*(-1 + muu1**2 + muu2**2)*(muu1**2 + muu3**2 - tred)) + 
     &               g1R*muu1*(2*g2L*muu2*muu3*(-1 + muu2**2 + muu3**2 - tred) + g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*
     &             dconjg(FR(indj,i,tt,gen)))) - muu4*xd11*
     &       (dconjg(ELc(indi,i,tt,gen))*((g1R*muu1*(2*g2R*muu2*muu3*muu4**2 + 
     &                  g2L*(-1 + muu3**2 + muu4**2)*(muu2**2 + muu4**2 - tred)) + 
     &               g1L*muu4*(2*g2R*muu2*muu3*(-1 + muu2**2 + muu3**2 - tred) + g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*
     &             dconjg(FL(indj,i,tt,gen)) + (g1R*muu1*muu4*
     &                (2*g2L*muu2*(-1 + muu3**2 + muu4**2) + g2R*muu3*(muu2**2 + muu4**2 - tred)) + 
     &               g1L*(g2L*muu2*muu4**2*(muu1**2 + muu3**2 - tred) + 
     &                  g2R*muu3*(-1 + muu2**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred)))*dconjg(FR(indj,i,tt,gen))) + 
     &         dconjg(ER(indi,i,tt,gen))*((g1L*muu1*muu4*
     &                (2*g2R*muu2*(-1 + muu3**2 + muu4**2) + g2L*muu3*(muu2**2 + muu4**2 - tred)) + 
     &               g1R*(g2R*muu2*muu4**2*(muu1**2 + muu3**2 - tred) + 
     &                  g2L*muu3*(-1 + muu2**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1L*muu1*(2*g2L*muu2*muu3*muu4**2 + g2R*(-1 + muu3**2 + muu4**2)*(muu2**2 + muu4**2 - tred)) + 
     &               g1R*muu4*(2*g2L*muu2*muu3*(-1 + muu2**2 + muu3**2 - tred) + g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*
     &             dconjg(FR(indj,i,tt,gen)))) - (2*xd00*
     &         (dconjg(ELc(indi,i,tt,gen))*((g1R*muu1*muu4*(4*g2R*muu2*muu3 + g2L*(-1 + muu1**2 + muu4**2 - tred)) + 
     &                 g1L*(g2R*muu2*muu3*(-1 + muu2**2 + muu3**2 - tred) + 2*g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*
     &               dconjg(FL(indj,i,tt,gen)) + (g1L*muu4*
     &                  (g2R*(-1 + muu1**2 + muu2**2)*muu3 + 2*g2L*muu2*(muu1**2 + muu3**2 - tred)) + 
     &                 g1R*muu1*(g2L*muu2*(-1 + muu3**2 + muu4**2) + 2*g2R*muu3*(muu2**2 + muu4**2 - tred)))*
     &               dconjg(FR(indj,i,tt,gen))) + dconjg(ER(indi,i,tt,gen))*
     &            ((g1R*muu4*(g2L*(-1 + muu1**2 + muu2**2)*muu3 + 2*g2R*muu2*(muu1**2 + muu3**2 - tred)) + 
     &                 g1L*muu1*(g2R*muu2*(-1 + muu3**2 + muu4**2) + 2*g2L*muu3*(muu2**2 + muu4**2 - tred)))*
     &               dconjg(FL(indj,i,tt,gen)) + (g1L*muu1*muu4*(4*g2L*muu2*muu3 + g2R*(-1 + muu1**2 + muu4**2 - tred)) + 
     &                 g1R*(g2L*muu2*muu3*(-1 + muu2**2 + muu3**2 - tred) + 2*g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*
     &               dconjg(FR(indj,i,tt,gen)))))/sqrtS**2 + 
     &      xd13*(dconjg(ELc(indi,i,tt,gen))*((g1R*muu1*muu4*
     &                (2*g2R*muu2*muu3*(-1 + muu2**2 + muu3**2 - tred) + 
     &                  g2L*(1 - muu3**2 - muu4**2 + muu3**2*muu4**2 + muu2**2*(-1 + 2*muu3**2 + muu4**2 - tred) + 
     &                     muu1**2*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred) - muu3**2*tred - muu4**2*tred + tred**2)) + 
     &               g1L*(g2L*(-1 + muu2**2 + muu3**2 - tred)*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred) + 
     &                  g2R*muu2*muu3*(1 + muu2**4 + muu3**4 + 4*muu1**2*muu4**2 + 2*muu2**2*(-1 + muu3**2 - tred) + 2*tred + 
     &                     tred**2 - 2*muu3**2*(1 + tred))))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1L*muu4*(g2R*muu3*(muu1**2*(-1 + 3*muu2**2 + muu3**2 + 2*muu4**2 - 3*tred) + 
     &                     (-1 + muu2**2)*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                  g2L*muu2*(muu1**2 + muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred)) + 
     &               g1R*muu1*(g2R*muu3*(-1 + muu2**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred) + 
     &                  g2L*muu2*(1 + muu3**4 - muu4**2 + 2*muu1**2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) + 
     &                     muu3**2*(-2 + 3*muu4**2 - tred) + tred - 3*muu4**2*tred)))*dconjg(FR(indj,i,tt,gen))) + 
     &         dconjg(ER(indi,i,tt,gen))*((g1R*muu4*(g2L*muu3*
     &                   (muu1**2*(-1 + 3*muu2**2 + muu3**2 + 2*muu4**2 - 3*tred) + (-1 + muu2**2)*(-1 + muu2**2 + muu3**2 - tred))
     &                    + g2R*muu2*(muu1**2 + muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred)) + 
     &               g1L*muu1*(g2L*muu3*(-1 + muu2**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred) + 
     &                  g2R*muu2*(1 + muu3**4 - muu4**2 + 2*muu1**2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) + 
     &                     muu3**2*(-2 + 3*muu4**2 - tred) + tred - 3*muu4**2*tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1L*muu1*muu4*(2*g2L*muu2*muu3*(-1 + muu2**2 + muu3**2 - tred) + 
     &                  g2R*(1 - muu3**2 - muu4**2 + muu3**2*muu4**2 + muu2**2*(-1 + 2*muu3**2 + muu4**2 - tred) + 
     &                     muu1**2*(-1 + muu2**2 + muu3**2 + 2*muu4**2 - tred) - muu3**2*tred - muu4**2*tred + tred**2)) + 
     &               g1R*(g2R*(-1 + muu2**2 + muu3**2 - tred)*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred) + 
     &                  g2L*muu2*muu3*(1 + muu2**4 + muu3**4 + 4*muu1**2*muu4**2 + 2*muu2**2*(-1 + muu3**2 - tred) + 2*tred + 
     &                     tred**2 - 2*muu3**2*(1 + tred))))*dconjg(FR(indj,i,tt,gen)))) - 
     &      muu4*xd0*(dconjg(ELc(indi,i,tt,gen))*((g1R*muu1*
     &                (g2R*muu2*muu3*(1 - muu3**2 + 3*muu4**2) + 
     &                  g2L*(-(muu1**2*muu4**2) + 2*muu3**2*muu4**2 + muu2**2*(-1 + 2*muu3**2 + muu4**2) + tred - 2*muu3**2*tred))
     &                + g1L*muu4*(g2R*muu2*muu3*(-1 + muu1**2 + muu2**2 + 2*muu3**2 - 2*tred) + 
     &                  g2L*(-muu1**4 + muu2**2*(2*muu3**2 - tred) - tred + muu1**2*(1 + muu2**2 + tred))))*
     &             dconjg(FL(indj,i,tt,gen)) + (g1R*muu1*muu4*
     &                (g2L*muu2*(-1 + 3*muu3**2 + muu4**2) + g2R*muu3*(1 - muu1**2 + 2*muu2**2 + muu4**2 - tred)) + 
     &               g1L*(g2R*muu3*(muu2**4 + muu2**2*(-1 + 2*muu3**2 - 2*tred) + (muu3**2 - tred)*(muu4**2 - tred)) + 
     &                  g2L*muu2*(muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(1 - muu3**2 + muu4**2) - tred - muu4**2*tred)))*
     &             dconjg(FR(indj,i,tt,gen))) + dconjg(ER(indi,i,tt,gen))*
     &          ((g1L*muu1*muu4*(g2R*muu2*(-1 + 3*muu3**2 + muu4**2) + g2L*muu3*(1 - muu1**2 + 2*muu2**2 + muu4**2 - tred)) + 
     &               g1R*(g2L*muu3*(muu2**4 + muu2**2*(-1 + 2*muu3**2 - 2*tred) + (muu3**2 - tred)*(muu4**2 - tred)) + 
     &                  g2R*muu2*(muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(1 - muu3**2 + muu4**2) - tred - muu4**2*tred)))*
     &             dconjg(FL(indj,i,tt,gen)) + (g1L*muu1*
     &                (g2L*muu2*muu3*(1 - muu3**2 + 3*muu4**2) + 
     &                  g2R*(-(muu1**2*muu4**2) + 2*muu3**2*muu4**2 + muu2**2*(-1 + 2*muu3**2 + muu4**2) + tred - 2*muu3**2*tred))
     &                + g1R*muu4*(g2L*muu2*muu3*(-1 + muu1**2 + muu2**2 + 2*muu3**2 - 2*tred) + 
     &                  g2R*(-muu1**4 + muu2**2*(2*muu3**2 - tred) - tred + muu1**2*(1 + muu2**2 + tred))))*
     &             dconjg(FR(indj,i,tt,gen)))) + xd12*
     &       (dconjg(ELc(indi,i,tt,gen))*((g1R*muu1*muu4*
     &                (2*g2R*muu2*muu3*(-1 + muu3**2 - muu4**2) + 
     &                  g2L*(1 + 2*muu2**2 - muu3**2 - 2*muu2**2*muu4**2 + muu3**2*muu4**2 - muu4**4 + 
     &                     muu1**2*(-1 + muu3**2 + muu4**2) - tred - muu3**2*tred + muu4**2*tred)) + 
     &               g1L*(g2L*(-1 + muu3**2 - muu4**2)*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred) + 
     &                  g2R*muu2*muu3*(1 + muu3**4 + 3*muu4**2 + 2*muu1**2*muu4**2 + muu2**2*(-1 + muu3**2 - 3*muu4**2) + tred + 
     &                     muu4**2*tred - muu3**2*(2 + muu4**2 + tred))))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1L*muu4*(g2R*muu3*(1 - 2*muu2**4 - muu3**2 + muu4**2 - muu2**2*(-1 + muu3**2 + muu4**2 - 2*tred) + 
     &                     muu1**2*(-1 + 2*muu2**2 + muu3**2 + muu4**2 - 2*tred)) + 
     &                  g2L*muu2*(-1 + muu3**2 - muu4**2)*(muu1**2 + muu3**2 - tred)) + 
     &               g1R*muu1*(g2L*muu2*(1 + muu3**4 + 2*muu4**2 - 3*muu4**4 + 2*muu3**2*(-1 + muu4**2)) + 
     &                  g2R*muu3*(-1 + muu3**2 - muu4**2)*(muu2**2 + muu4**2 - tred)))*dconjg(FR(indj,i,tt,gen))) + 
     &         dconjg(ER(indi,i,tt,gen))*((g1R*muu4*(g2L*muu3*
     &                   (1 - 2*muu2**4 - muu3**2 + muu4**2 - muu2**2*(-1 + muu3**2 + muu4**2 - 2*tred) + 
     &                     muu1**2*(-1 + 2*muu2**2 + muu3**2 + muu4**2 - 2*tred)) + 
     &                  g2R*muu2*(-1 + muu3**2 - muu4**2)*(muu1**2 + muu3**2 - tred)) + 
     &               g1L*muu1*(g2R*muu2*(1 + muu3**4 + 2*muu4**2 - 3*muu4**4 + 2*muu3**2*(-1 + muu4**2)) + 
     &                  g2L*muu3*(-1 + muu3**2 - muu4**2)*(muu2**2 + muu4**2 - tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1L*muu1*muu4*(2*g2L*muu2*muu3*(-1 + muu3**2 - muu4**2) + 
     &                  g2R*(1 + 2*muu2**2 - muu3**2 - 2*muu2**2*muu4**2 + muu3**2*muu4**2 - muu4**4 + 
     &                     muu1**2*(-1 + muu3**2 + muu4**2) - tred - muu3**2*tred + muu4**2*tred)) + 
     &               g1R*(g2R*(-1 + muu3**2 - muu4**2)*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred) + 
     &                  g2L*muu2*muu3*(1 + muu3**4 + 3*muu4**2 + 2*muu1**2*muu4**2 + muu2**2*(-1 + muu3**2 - 3*muu4**2) + tred + 
     &                     muu4**2*tred - muu3**2*(2 + muu4**2 + tred))))*dconjg(FR(indj,i,tt,gen)))) + 
     &      xd23*(dconjg(ELc(indi,i,tt,gen))*((g1R*muu1*muu4*
     &                (2*g2R*muu2*(-1 - muu1**2 + muu2**2)*muu3 - 
     &                  g2L*(-1 + muu1**4 - 2*muu3**2 + muu4**2 + tred + muu2**2*(1 - muu4**2 + tred) - 
     &                     muu1**2*(muu2**2 - 2*muu3**2 + muu4**2 + tred))) + 
     &               g1L*(g2L*(-1 - muu1**2 + muu2**2)*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred) + 
     &                  g2R*muu2*muu3*((-1 + muu2**2)*(-1 + muu2**2 + muu3**2 - tred) + 
     &                     muu1**2*(3 - muu2**2 - 3*muu3**2 + 2*muu4**2 + tred))))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1L*muu4*(g2R*(-3*muu1**4 + (-1 + muu2**2)**2 + 2*muu1**2*(1 + muu2**2))*muu3 + 
     &                  g2L*muu2*(-1 - muu1**2 + muu2**2)*(muu1**2 + muu3**2 - tred)) + 
     &               g1R*muu1*(g2R*(-1 - muu1**2 + muu2**2)*muu3*(muu2**2 + muu4**2 - tred) + 
     &                  g2L*muu2*(1 + muu3**2 - 2*muu3**4 - muu4**2 + 2*muu3**2*muu4**2 + muu2**2*(-1 - muu3**2 + muu4**2) + 
     &                     muu1**2*(1 - muu3**2 + muu4**2) + 2*muu3**2*tred - 2*muu4**2*tred)))*dconjg(FR(indj,i,tt,gen))) + 
     &         dconjg(ER(indi,i,tt,gen))*((g1R*muu4*(g2L*(-3*muu1**4 + (-1 + muu2**2)**2 + 2*muu1**2*(1 + muu2**2))*muu3 + 
     &                  g2R*muu2*(-1 - muu1**2 + muu2**2)*(muu1**2 + muu3**2 - tred)) + 
     &               g1L*muu1*(g2L*(-1 - muu1**2 + muu2**2)*muu3*(muu2**2 + muu4**2 - tred) + 
     &                  g2R*muu2*(1 + muu3**2 - 2*muu3**4 - muu4**2 + 2*muu3**2*muu4**2 + muu2**2*(-1 - muu3**2 + muu4**2) + 
     &                     muu1**2*(1 - muu3**2 + muu4**2) + 2*muu3**2*tred - 2*muu4**2*tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1L*muu1*muu4*(2*g2L*muu2*(-1 - muu1**2 + muu2**2)*muu3 - 
     &                  g2R*(-1 + muu1**4 - 2*muu3**2 + muu4**2 + tred + muu2**2*(1 - muu4**2 + tred) - 
     &                     muu1**2*(muu2**2 - 2*muu3**2 + muu4**2 + tred))) + 
     &               g1R*(g2R*(-1 - muu1**2 + muu2**2)*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred) + 
     &                  g2L*muu2*muu3*((-1 + muu2**2)*(-1 + muu2**2 + muu3**2 - tred) + 
     &                     muu1**2*(3 - muu2**2 - 3*muu3**2 + 2*muu4**2 + tred))))*dconjg(FR(indj,i,tt,gen)))) + 
     &      xd3*(dconjg(ELc(indi,i,tt,gen))*((g1R*muu1*muu4*
     &                (g2R*muu2*muu3*(-3 + muu1**2 + 3*muu2**2 + 4*muu3**2 - 4*tred) + 
     &                  g2L*(-muu1**4 + muu2**2*(3*muu3**2 - tred) - 2*tred + muu1**2*(1 + muu2**2 + muu4**2 + tred))) + 
     &               g1L*(g2R*muu1**2*muu2*muu3*(1 - muu3**2 + 5*muu4**2) + 
     &                  g2L*((-1 + muu2**2)*muu4**2*(muu3**2 - tred) + 
     &                     muu1**2*(muu2**2*(-1 + 2*muu3**2 + 2*muu4**2) + muu4**2*(-1 + 3*muu3**2 - tred) + tred - 2*muu3**2*tred)
     &                     )))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1L*muu4*(g2L*muu2*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + 
     &                     muu1**2*(-2 + muu2**2 + 4*muu3**2 + muu4**2 - tred)) + 
     &                  g2R*muu1**2*muu3*(1 - muu1**2 + 3*muu2**2 + 2*muu4**2 - 2*tred)) + 
     &               g1R*muu1*(g2L*muu2*(muu2**2*muu3**2 + 2*muu3**2*muu4**2 + muu1**2*(1 - muu3**2 + 2*muu4**2) - tred - 
     &                     2*muu4**2*tred) + g2R*muu3*
     &                   (muu2**4 + muu2**2*(-1 + 2*muu3**2 + 2*muu4**2 - 2*tred) + muu4**2*(-2 + 2*muu1**2 + muu3**2 - tred) + 
     &                     tred*(-muu3**2 + tred))))*dconjg(FR(indj,i,tt,gen))) + 
     &         dconjg(ER(indi,i,tt,gen))*((g1R*muu4*(g2R*muu2*
     &                   ((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(-2 + muu2**2 + 4*muu3**2 + muu4**2 - tred))
     &                   + g2L*muu1**2*muu3*(1 - muu1**2 + 3*muu2**2 + 2*muu4**2 - 2*tred)) + 
     &               g1L*muu1*(g2R*muu2*(muu2**2*muu3**2 + 2*muu3**2*muu4**2 + muu1**2*(1 - muu3**2 + 2*muu4**2) - tred - 
     &                     2*muu4**2*tred) + g2L*muu3*
     &                   (muu2**4 + muu2**2*(-1 + 2*muu3**2 + 2*muu4**2 - 2*tred) + muu4**2*(-2 + 2*muu1**2 + muu3**2 - tred) + 
     &                     tred*(-muu3**2 + tred))))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1L*muu1*muu4*(g2L*muu2*muu3*(-3 + muu1**2 + 3*muu2**2 + 4*muu3**2 - 4*tred) + 
     &                  g2R*(-muu1**4 + muu2**2*(3*muu3**2 - tred) - 2*tred + muu1**2*(1 + muu2**2 + muu4**2 + tred))) + 
     &               g1R*(g2L*muu1**2*muu2*muu3*(1 - muu3**2 + 5*muu4**2) + 
     &                  g2R*((-1 + muu2**2)*muu4**2*(muu3**2 - tred) + 
     &                     muu1**2*(muu2**2*(-1 + 2*muu3**2 + 2*muu4**2) + muu4**2*(-1 + 3*muu3**2 - tred) + tred - 2*muu3**2*tred)
     &                     )))*dconjg(FR(indj,i,tt,gen)))) + 
     &      xd1*(dconjg(ELc(indi,i,tt,gen))*((g1R*muu1*muu4*
     &                (-(g2R*muu2*muu3*(-1 + muu3**2 + 7*muu4**2)) - 
     &                  g2L*(-1 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - muu1**2*(-1 + muu3**2 + 2*muu4**2) + 
     &                     muu2**2*(-2 + 3*muu3**2 + 2*muu4**2) + tred - 2*muu3**2*tred)) + 
     &               g1L*(g2R*muu2*muu3*(1 + muu3**4 + 2*muu4**2 - muu1**2*muu4**2 + muu2**2*(-1 + muu3**2 - 2*muu4**2) + tred + 
     &                     3*muu4**2*tred - muu3**2*(2 + 3*muu4**2 + tred)) + 
     &                  g2L*(muu1**4*muu4**2 + muu3**4*(-muu4**2 + tred) - 
     &                     muu1**2*(muu4**2 + muu4**4 + muu2**2*(-1 + muu3**2 + 2*muu4**2) + tred - muu3**2*tred) + 
     &                     tred*(1 + muu4**4 + muu2**2*(-1 + 2*muu4**2) + tred - muu4**2*tred) - 
     &                     muu3**2*(muu4**4 + muu4**2*(-1 + 3*muu2**2 - 2*tred) + tred*(2 - muu2**2 + tred)))))*
     &             dconjg(FL(indj,i,tt,gen)) + (g1R*muu1*
     &                (g2L*muu2*(1 + muu3**4 + muu4**2 - 2*muu4**4 - muu3**2*(2 + 3*muu4**2)) + 
     &                  g2R*muu3*(muu1**2*muu4**2 - muu3**2*muu4**2 - 3*muu4**4 - muu2**2*(-1 + muu3**2 + 4*muu4**2) - tred + 
     &                     muu3**2*tred + 3*muu4**2*tred)) + 
     &               g1L*muu4*(g2L*muu2*(-muu3**4 + 3*muu4**2*(-muu1**2 + tred) + muu3**2*(1 - muu2**2 - 3*muu4**2 + tred)) + 
     &                  g2R*muu3*(1 - 2*muu2**4 - muu3**2 - 2*muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) - tred + 
     &                     2*muu3**2*tred + 2*muu4**2*tred - 2*tred**2 + muu2**2*(1 - 2*muu3**2 + 4*tred))))*
     &             dconjg(FR(indj,i,tt,gen))) + dconjg(ER(indi,i,tt,gen))*
     &          ((g1L*muu1*(g2R*muu2*(1 + muu3**4 + muu4**2 - 2*muu4**4 - muu3**2*(2 + 3*muu4**2)) + 
     &                  g2L*muu3*(muu1**2*muu4**2 - muu3**2*muu4**2 - 3*muu4**4 - muu2**2*(-1 + muu3**2 + 4*muu4**2) - tred + 
     &                     muu3**2*tred + 3*muu4**2*tred)) + 
     &               g1R*muu4*(g2R*muu2*(-muu3**4 + 3*muu4**2*(-muu1**2 + tred) + muu3**2*(1 - muu2**2 - 3*muu4**2 + tred)) + 
     &                  g2L*muu3*(1 - 2*muu2**4 - muu3**2 - 2*muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) - tred + 
     &                     2*muu3**2*tred + 2*muu4**2*tred - 2*tred**2 + muu2**2*(1 - 2*muu3**2 + 4*tred))))*
     &             dconjg(FL(indj,i,tt,gen)) + (g1L*muu1*muu4*
     &                (-(g2L*muu2*muu3*(-1 + muu3**2 + 7*muu4**2)) - 
     &                  g2R*(-1 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - muu1**2*(-1 + muu3**2 + 2*muu4**2) + 
     &                     muu2**2*(-2 + 3*muu3**2 + 2*muu4**2) + tred - 2*muu3**2*tred)) + 
     &               g1R*(g2L*muu2*muu3*(1 + muu3**4 + 2*muu4**2 - muu1**2*muu4**2 + muu2**2*(-1 + muu3**2 - 2*muu4**2) + tred + 
     &                     3*muu4**2*tred - muu3**2*(2 + 3*muu4**2 + tred)) + 
     &                  g2R*(muu1**4*muu4**2 + muu3**4*(-muu4**2 + tred) - 
     &                     muu1**2*(muu4**2 + muu4**4 + muu2**2*(-1 + muu3**2 + 2*muu4**2) + tred - muu3**2*tred) + 
     &                     tred*(1 + muu4**4 + muu2**2*(-1 + 2*muu4**2) + tred - muu4**2*tred) - 
     &                     muu3**2*(muu4**4 + muu4**2*(-1 + 3*muu2**2 - 2*tred) + tred*(2 - muu2**2 + tred)))))*
     &             dconjg(FR(indj,i,tt,gen)))) + xd2*
     &       (dconjg(ELc(indi,i,tt,gen))*(-((g1R*muu1*muu4*
     &                  (2*g2R*muu2*muu3*(1 - 2*muu3**2 + 2*muu4**2) + 
     &                    g2L*(-1 + muu4**2 + muu1**2*(1 - 2*muu4**2) + 2*muu2**2*(-1 + muu4**2) + tred)) + 
     &                 g1L*(g2R*muu2*muu3*(-1 + muu3**2 - 2*muu4**2 - 2*muu1**2*muu4**2 + muu2**2*(1 + 2*muu4**2) - tred) - 
     &                    g2L*(-(muu2**4*muu3**2) + muu1**4*muu4**2 + (1 - muu3**2 + muu4**2)*tred + 
     &                       muu1**2*(muu2**2*(muu3**2 - muu4**2) + muu4**2*(-2 + 2*muu3**2 - tred) - muu3**2*tred) + 
     &                       muu2**2*(muu4**2*tred + muu3**2*(-2*muu4**2 + tred)))))*dconjg(FL(indj,i,tt,gen))) + 
     &            (g1R*muu1*(g2L*muu2*(1 + muu4**2 - 2*muu4**4 + muu3**2*(-1 + 2*muu4**2)) + 
     &                  g2R*muu3*(-muu4**4 + muu2**2*(muu3**2 - 3*muu4**2) - muu3**2*tred + 
     &                     muu4**2*(-2 + 2*muu1**2 + muu3**2 + tred))) + 
     &               g1L*muu4*(g2L*muu2*(-2*muu2**2*muu3**2 + muu3**4 - muu3**2*muu4**2 + muu1**2*(-2 + 3*muu3**2 - muu4**2) + 
     &                     2*tred - muu3**2*tred + muu4**2*tred) + 
     &                  g2R*muu3*(1 - 2*muu2**4 + muu1**2*(-1 + 2*muu2**2 + 2*muu4**2 - 2*tred) + muu2**2*(1 - 2*muu3**2 + 2*tred))
     &                  ))*dconjg(FR(indj,i,tt,gen))) + 
     &         dconjg(ER(indi,i,tt,gen))*((g1L*muu1*(g2R*muu2*(1 + muu4**2 - 2*muu4**4 + muu3**2*(-1 + 2*muu4**2)) + 
     &                  g2L*muu3*(-muu4**4 + muu2**2*(muu3**2 - 3*muu4**2) - muu3**2*tred + 
     &                     muu4**2*(-2 + 2*muu1**2 + muu3**2 + tred))) + 
     &               g1R*muu4*(g2R*muu2*(-2*muu2**2*muu3**2 + muu3**4 - muu3**2*muu4**2 + muu1**2*(-2 + 3*muu3**2 - muu4**2) + 
     &                     2*tred - muu3**2*tred + muu4**2*tred) + 
     &                  g2L*muu3*(1 - 2*muu2**4 + muu1**2*(-1 + 2*muu2**2 + 2*muu4**2 - 2*tred) + muu2**2*(1 - 2*muu3**2 + 2*tred))
     &                  ))*dconjg(FL(indj,i,tt,gen)) - 
     &            (g1L*muu1*muu4*(2*g2L*muu2*muu3*(1 - 2*muu3**2 + 2*muu4**2) + 
     &                  g2R*(-1 + muu4**2 + muu1**2*(1 - 2*muu4**2) + 2*muu2**2*(-1 + muu4**2) + tred)) + 
     &               g1R*(g2L*muu2*muu3*(-1 + muu3**2 - 2*muu4**2 - 2*muu1**2*muu4**2 + muu2**2*(1 + 2*muu4**2) - tred) - 
     &                  g2R*(-(muu2**4*muu3**2) + muu1**4*muu4**2 + (1 - muu3**2 + muu4**2)*tred + 
     &                     muu1**2*(muu2**2*(muu3**2 - muu4**2) + muu4**2*(-2 + 2*muu3**2 - tred) - muu3**2*tred) + 
     &                     muu2**2*(muu4**2*tred + muu3**2*(-2*muu4**2 + tred)))))*dconjg(FR(indj,i,tt,gen))))))/
     &  (4.*Pi**2*(-tred + MSf(i,itypeT,gen)**2/sqrtS**2))

	
	tempBoxt = tempBoxt + ChiChi2QQ_M2BoxGluon - tempBox2
	tempBox2 = ChiChi2QQ_M2BoxGluon

	enddo

	enddo

c	write(*,*)'ChiChi U-Channel Gluon Box (Neu) = ', ChiChi2QQ_M2BoxGluon
c	write(*,*)'ChiChi U-Channel Gluon Box (Neu) mit Z = ', tempBoxZ
c	write(*,*)'ChiChi U-Channel Gluon Box (Neu) mit H = ', tempBoxH
c	write(*,*)'ChiChi U-Channel Gluon Box (Neu) mit t = ', tempBoxt
c	write(*,*)'ChiChi U-Channel Gluon Box (Neu) mit u = ', tempBoxu
	tempBox = ChiChi2QQ_M2BoxGluon
	tempBoxZ = 0d0
	tempBoxH = 0d0
	tempBoxt = 0d0 
	tempBoxu = 0d0


!**********************************
!****** Gluon T-channel Box *******
!**********************************

	do k = 1,2

	g1L = ELc(indi,k,tt,gen)
	g1R = ER(indi,k,tt,gen)
	g2L = FL(indj,k,tt,gen) 
	g2R = FR(indj,k,tt,gen)
	g3L = -GS
	g4L = -GS

	call myDget(muu3**2*sqrtS**2,muu4**2*sqrtS**2,muu2**2*sqrtS**2,muu1**2*sqrtS**2,sqrtS**2,tred*sqrtS**2,muu3**2*sqrtS**2,0d0,muu4**2*sqrtS**2,MSf2(k,itypeT,gen))
	  
	  
! Photon contr.

	if (icha.eq.2) then
	
	ChiChi2QQ_M2BoxGluon = ChiChi2QQ_M2BoxGluon + (ApL*BpL*CF*g3L*g4L*NC*sqrtS**4*(-2*muu3*(g1L*(g2R*muu3*
     &             (muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*tred) + (1 + tred)**2 - muu3**2*(1 + 2*tred)) + 
     &            g2L*muu1*(1 - 4*muu3**4 + muu3**2*(4 - 4*muu1**2 + 4*tred))) + 
     &         g1R*(g2L*muu3*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*tred) + (1 + tred)**2 - muu3**2*(1 + 2*tred)) + 
     &            g2R*muu1*(1 - 4*muu3**4 + muu3**2*(4 - 4*muu1**2 + 4*tred))))*xd0 + 
     &      (4*(g1R*(4*g2R*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred) + 
     &              g2L*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)) + 
     &           g1L*(4*g2L*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred) + 
     &              g2R*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)))*xd00)/sqrtS**2 + 
     &      2*(g1R*(g2R*muu1*muu3*(1 + 4*muu3**4 + muu1**2*(2 + 4*muu3**2) - 2*tred - 4*muu3**2*(2 + tred)) + 
     &            g2L*(-4*muu3**6 + muu1**4*(1 - 4*muu3**2) + tred**2 - 2*muu3**2*tred*(3 + 2*tred) + muu3**4*(1 + 8*tred) + 
     &               muu1**2*(1 - 8*muu3**4 - 2*tred + muu3**2*(4 + 8*tred)))) + 
     &         g1L*(g2L*muu1*muu3*(1 + 4*muu3**4 + muu1**2*(2 + 4*muu3**2) - 2*tred - 4*muu3**2*(2 + tred)) + 
     &            g2R*(-4*muu3**6 + muu1**4*(1 - 4*muu3**2) + tred**2 - 2*muu3**2*tred*(3 + 2*tred) + muu3**4*(1 + 8*tred) + 
     &               muu1**2*(1 - 8*muu3**4 - 2*tred + muu3**2*(4 + 8*tred)))))*xd1 - 
     &      2*muu3**2*(g1R*(2*g2R*muu1*muu3*(1 - 2*muu1**2 - 2*muu3**2 + 2*tred) + 
     &            g2L*(muu1**4 + muu3**4 + muu1**2*(-3 + 2*muu3**2 - 2*tred) + tred*(2 + tred) - muu3**2*(1 + 2*tred))) + 
     &         g1L*(2*g2L*muu1*muu3*(1 - 2*muu1**2 - 2*muu3**2 + 2*tred) + 
     &            g2R*(muu1**4 + muu3**4 + muu1**2*(-3 + 2*muu3**2 - 2*tred) + tred*(2 + tred) - muu3**2*(1 + 2*tred))))*xd11 + 
     &      2*(g1R*(2*g2R*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred) + 
     &            g2L*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)) + 
     &         g1L*(2*g2L*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred) + 
     &            g2R*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)))*xd12 + 
     &      2*(g1R*(2*g2R*muu1*muu3*(2*muu1**4 + 2*muu3**4 + muu1**2*(-1 + 4*muu3**2 - 4*tred) + tred*(1 + 2*tred) - muu3**2*(1 + 4*tred)) + 
     &            g2L*(muu1**4*(1 - 4*muu3**2) + (muu3**2 - tred)**2 - 2*muu1**2*(2*muu3**4 + tred - muu3**2*(3 + 2*tred)))) + 
     &         g1L*(2*g2L*muu1*muu3*(2*muu1**4 + 2*muu3**4 + muu1**2*(-1 + 4*muu3**2 - 4*tred) + tred*(1 + 2*tred) - muu3**2*(1 + 4*tred)) + 
     &            g2R*(muu1**4*(1 - 4*muu3**2) + (muu3**2 - tred)**2 - 2*muu1**2*(2*muu3**4 + tred - muu3**2*(3 + 2*tred)))))*xd13 + 
     &      2*(g1R*(2*g2R*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred) + 
     &            g2L*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)) + 
     &         g1L*(2*g2L*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred) + 
     &            g2R*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)))*xd2 + 
     &      2*(g1R*(2*g2R*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred) + 
     &            g2L*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)) + 
     &         g1L*(2*g2L*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred) + 
     &            g2R*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)))*xd22 + 
     &      2*(g1R*(2*g2R*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred) + 
     &            g2L*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)) + 
     &         g1L*(2*g2L*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred) + 
     &            g2R*(muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)))*xd23 + 
     &      2*muu1*(g1R*(g2R*muu3*(1 + 2*muu1**4 + 2*muu3**4 + muu1**2*(-2 + 4*muu3**2 - 4*tred) + 4*tred + 2*tred**2 - 4*muu3**2*(1 + tred)) + 
     &            g2L*muu1*(1 - 8*muu3**4 + muu3**2*(6 - 8*muu1**2 + 8*tred))) + 
     &         g1L*(g2L*muu3*(1 + 2*muu1**4 + 2*muu3**4 + muu1**2*(-2 + 4*muu3**2 - 4*tred) + 4*tred + 2*tred**2 - 4*muu3**2*(1 + tred)) + 
     &            g2R*muu1*(1 - 8*muu3**4 + muu3**2*(6 - 8*muu1**2 + 8*tred))))*xd3 - 
     &      2*muu1**2*(g1R*(2*g2R*muu1*muu3*(1 - 2*muu1**2 - 2*muu3**2 + 2*tred) + 
     &            g2L*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*tred) + tred*(2 + tred) - muu3**2*(3 + 2*tred))) + 
     &         g1L*(2*g2L*muu1*muu3*(1 - 2*muu1**2 - 2*muu3**2 + 2*tred) + 
     &            g2R*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*tred) + tred*(2 + tred) - muu3**2*(3 + 2*tred))))*xd33))/(4.*Pi**2)
	
	endif

! Zboson contr.
	
	ChiChi2QQ_M2BoxGluon = ChiChi2QQ_M2BoxGluon -(CF*g3L*g4L*NC*sqrtS**6*(-2*muu1*xd33*((-(BL*(g1L*muu1*muu3*(2*g2R*(-1 + muu1**2 + muu2**2)*muu4 - g2L*muu2*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                  g1R*(g2L*muu2*(muu1**2 + muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) - g2R*muu1**2*muu4*(-1 + muu1**2 + muu4**2 - tred)))) + 
     &             BR*muu1*(g1R*muu1*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &                g1L*(g2L*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred))))*dconjg(AL) + 
     &          (-(BR*(g1R*muu1*muu3*(2*g2L*(-1 + muu1**2 + muu2**2)*muu4 - g2R*muu2*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                  g1L*(g2R*muu2*(muu1**2 + muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) - g2L*muu1**2*muu4*(-1 + muu1**2 + muu4**2 - tred)))) + 
     &             BL*muu1*(g1L*muu1*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &                g1R*(g2R*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred))))*dconjg(AR)) - 
     &       2*muu3*xd11*((-(BL*(g1R*muu1*muu3*(2*g2L*muu2*(-1 + muu3**2 + muu4**2) - g2R*muu4*(-1 + muu1**2 + muu4**2 - tred)) + 
     &                  g1L*(-(g2L*muu2*muu3**2*(-1 + muu2**2 + muu3**2 - tred)) + g2R*muu4*(muu1**2 + muu3**2 - tred)*(-1 + muu1**2 + muu4**2 - tred))))
     &              + BR*muu3*(g1R*muu1*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &                g1L*(g2L*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred))))*dconjg(AL) + 
     &          (-(BR*(g1L*muu1*muu3*(2*g2R*muu2*(-1 + muu3**2 + muu4**2) - g2L*muu4*(-1 + muu1**2 + muu4**2 - tred)) + 
     &                  g1R*(-(g2R*muu2*muu3**2*(-1 + muu2**2 + muu3**2 - tred)) + g2L*muu4*(muu1**2 + muu3**2 - tred)*(-1 + muu1**2 + muu4**2 - tred))))
     &              + BL*muu3*(g1L*muu1*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &                g1R*(g2R*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred))))*dconjg(AR)) - 
     &       (4*xd00*((-(BL*(g1L*muu3*(g2R*(-1 + muu1**2 + muu2**2)*muu4 - 2*g2L*muu2*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                    g1R*muu1*(g2L*muu2*(-1 + muu3**2 + muu4**2) - 2*g2R*muu4*(-1 + muu1**2 + muu4**2 - tred)))) + 
     &               BR*(2*g1R*muu1*muu3 + g1L*(muu1**2 + muu3**2 - tred))*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)))*dconjg(AL) + 
     &            (-(BR*(g1R*muu3*(g2L*(-1 + muu1**2 + muu2**2)*muu4 - 2*g2R*muu2*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                    g1L*muu1*(g2R*muu2*(-1 + muu3**2 + muu4**2) - 2*g2L*muu4*(-1 + muu1**2 + muu4**2 - tred)))) + 
     &               BL*(2*g1L*muu1*muu3 + g1R*(muu1**2 + muu3**2 - tred))*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)))*dconjg(AR)))/sqrtS**2 - 
     &       2*xd13*((BR*(g1R*muu1*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &                g1L*(g2L*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*(muu1**2 + muu3**2 - tred) - 
     &             BL*(g1L*muu3*(g2R*muu4*(3*muu1**4 + muu1**2*(-3 + muu2**2 + muu3**2 + 2*muu4**2 - 3*tred) + (-1 + muu2**2)*(muu3**2 - tred)) - 
     &                   g2L*muu2*(muu1**2 + muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                g1R*muu1*(-(g2R*muu4*(muu1**2 + muu3**2 - tred)*(-1 + muu1**2 + muu4**2 - tred)) + 
     &                   g2L*muu2*(3*muu3**4 + muu1**2*(-1 + muu3**2 + muu4**2) + muu3**2*(-3 + 2*muu2**2 + muu4**2 - 3*tred) + tred - muu4**2*tred))))*
     &           dconjg(AL) + (BL*(g1L*muu1*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &                g1R*(g2R*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*(muu1**2 + muu3**2 - tred) - 
     &             BR*(g1R*muu3*(g2L*muu4*(3*muu1**4 + muu1**2*(-3 + muu2**2 + muu3**2 + 2*muu4**2 - 3*tred) + (-1 + muu2**2)*(muu3**2 - tred)) - 
     &                   g2R*muu2*(muu1**2 + muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                g1L*muu1*(-(g2L*muu4*(muu1**2 + muu3**2 - tred)*(-1 + muu1**2 + muu4**2 - tred)) + 
     &                   g2R*muu2*(3*muu3**4 + muu1**2*(-1 + muu3**2 + muu4**2) + muu3**2*(-3 + 2*muu2**2 + muu4**2 - 3*tred) + tred - muu4**2*tred))))*
     &           dconjg(AR)) + 2*muu3*xd0*((BL*(g1R*muu1*muu3*
     &                 (g2L*muu2*(-1 + muu3**2 + 3*muu4**2) + g2R*muu4*(2 - 2*muu1**2 + muu2**2 - muu4**2 + tred)) + 
     &                g1L*(g2R*muu4*(muu1**4 + muu2**2*(-1 - muu3**2 + muu4**2) + muu1**2*(-2 + muu3**2 + muu4**2 - 2*tred) + 
     &                      (-1 + muu3**2 - tred)*(-1 + muu4**2 - tred)) - 
     &                   g2L*muu2*(-1 + muu3**4 + muu4**2 + muu1**2*muu4**2 + muu2**2*(1 + muu3**2 - muu4**2) - tred - muu3**2*tred))) + 
     &             BR*(-(g1R*muu1*(-1 + muu4**2)*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred))) + 
     &                g1L*muu3*(g2L*muu2*muu4*(2 + muu1**2 - 2*muu2**2 - muu3**2 + tred) + 
     &                   g2R*(-muu2**4 + muu4**2 - muu1**2*tred + muu2**2*(1 + muu1**2 - muu3**2 - muu4**2 + tred)))))*dconjg(AL) + 
     &          (BR*(g1L*muu1*muu3*(g2R*muu2*(-1 + muu3**2 + 3*muu4**2) + g2L*muu4*(2 - 2*muu1**2 + muu2**2 - muu4**2 + tred)) + 
     &                g1R*(g2L*muu4*(muu1**4 + muu2**2*(-1 - muu3**2 + muu4**2) + muu1**2*(-2 + muu3**2 + muu4**2 - 2*tred) + 
     &                      (-1 + muu3**2 - tred)*(-1 + muu4**2 - tred)) - 
     &                   g2R*muu2*(-1 + muu3**4 + muu4**2 + muu1**2*muu4**2 + muu2**2*(1 + muu3**2 - muu4**2) - tred - muu3**2*tred))) + 
     &             BL*(-(g1L*muu1*(-1 + muu4**2)*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred))) + 
     &                g1R*muu3*(g2R*muu2*muu4*(2 + muu1**2 - 2*muu2**2 - muu3**2 + tred) + 
     &                   g2L*(-muu2**4 + muu4**2 - muu1**2*tred + muu2**2*(1 + muu1**2 - muu3**2 - muu4**2 + tred)))))*dconjg(AR)) - 
     &       2*xd3*((BR*muu1*(-(g1L*muu1*(-1 + muu3**2 + muu4**2)*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred))) + 
     &                g1R*muu3*(-2*g2L*muu2*muu4*(-1 + muu2**2 + muu3**2 - tred) + 
     &                   g2R*(-muu2**4 + muu2**2*(1 + muu1**2 - 2*muu3**2 - muu4**2 + tred) - (-1 + muu1**2)*(muu4**2 + tred)))) + 
     &             BL*(g1L*muu3*(g2L*muu2*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(-2 + muu2**2 + 2*muu3**2 + 3*muu4**2 - tred)) + 
     &                   g2R*muu1**2*muu4*(3 - 3*muu1**2 + muu2**2 - 2*muu4**2 + 2*tred)) + 
     &                g1R*muu1*(g2R*muu4*(1 + muu1**4 - 3*muu3**2 - muu4**2 + muu3**2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) + 
     &                      muu1**2*(-2 + 3*muu3**2 + muu4**2 - 2*tred) + 2*tred - muu3**2*tred - muu4**2*tred + tred**2) + 
     &                   g2L*muu2*(1 - 2*muu3**4 - muu4**2 - muu1**2*muu4**2 + muu2**2*(-1 - 2*muu3**2 + muu4**2) + tred + muu3**2*(1 + 2*tred)))))*
     &           dconjg(AL) + (BL*muu1*(-(g1R*muu1*(-1 + muu3**2 + muu4**2)*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred))) + 
     &                g1L*muu3*(-2*g2R*muu2*muu4*(-1 + muu2**2 + muu3**2 - tred) + 
     &                   g2L*(-muu2**4 + muu2**2*(1 + muu1**2 - 2*muu3**2 - muu4**2 + tred) - (-1 + muu1**2)*(muu4**2 + tred)))) + 
     &             BR*(g1R*muu3*(g2R*muu2*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(-2 + muu2**2 + 2*muu3**2 + 3*muu4**2 - tred)) + 
     &                   g2L*muu1**2*muu4*(3 - 3*muu1**2 + muu2**2 - 2*muu4**2 + 2*tred)) + 
     &                g1L*muu1*(g2L*muu4*(1 + muu1**4 - 3*muu3**2 - muu4**2 + muu3**2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) + 
     &                      muu1**2*(-2 + 3*muu3**2 + muu4**2 - 2*tred) + 2*tred - muu3**2*tred - muu4**2*tred + tred**2) + 
     &                   g2R*muu2*(1 - 2*muu3**4 - muu4**2 - muu1**2*muu4**2 + muu2**2*(-1 - 2*muu3**2 + muu4**2) + tred + muu3**2*(1 + 2*tred)))))*
     &           dconjg(AR)) - 2*xd2*((BL*(g1R*muu1*(g2L*muu2*(1 - 2*muu3**4 - muu4**2 + muu3**2*(1 + 2*muu4**2)) + 
     &                   g2R*muu4*(-(muu4**2*(-1 + muu1**2 + muu4**2 - tred)) + muu3**2*(-3 + 3*muu1**2 - 2*muu2**2 + muu4**2 - tred))) + 
     &                g1L*muu3*(g2L*muu2*(-2 + muu3**4 + 3*muu4**2 + 2*muu1**2*muu4**2 + muu2**2*(2 + muu3**2 - 3*muu4**2) - 2*tred + muu4**2*tred - 
     &                      muu3**2*(-1 + muu4**2 + tred)) + g2R*muu4*
     &                    (1 - 2*muu1**4 + muu2**2*(-1 + 2*muu3**2 - 2*tred) + muu1**2*(1 + 2*muu2**2 - 2*muu4**2 + 2*tred)))) + 
     &             BR*(2*g1R*muu1*muu3*muu4*(g2R*(1 - muu1**2 + muu2**2)*muu4 + g2L*muu2*(1 - muu3**2 + muu4**2)) - 
     &                g1L*(2*g2L*muu2*muu4*(-(muu2**2*muu3**2) + muu1**2*(-1 + muu3**2) + tred) + 
     &                   g2R*(-2*muu2**4*muu3**2 + muu1**2*(-1 + 2*muu3**2)*(muu2**2 + muu4**2 - tred) + (muu3**2 - tred)*(-muu4**2 + tred) + 
     &                      muu2**2*(muu3**2 - 2*muu3**4 + tred + 2*muu3**2*tred)))))*dconjg(AL) + 
     &          (BR*(g1L*muu1*(g2R*muu2*(1 - 2*muu3**4 - muu4**2 + muu3**2*(1 + 2*muu4**2)) + 
     &                   g2L*muu4*(-(muu4**2*(-1 + muu1**2 + muu4**2 - tred)) + muu3**2*(-3 + 3*muu1**2 - 2*muu2**2 + muu4**2 - tred))) + 
     &                g1R*muu3*(g2R*muu2*(-2 + muu3**4 + 3*muu4**2 + 2*muu1**2*muu4**2 + muu2**2*(2 + muu3**2 - 3*muu4**2) - 2*tred + muu4**2*tred - 
     &                      muu3**2*(-1 + muu4**2 + tred)) + g2L*muu4*
     &                    (1 - 2*muu1**4 + muu2**2*(-1 + 2*muu3**2 - 2*tred) + muu1**2*(1 + 2*muu2**2 - 2*muu4**2 + 2*tred)))) + 
     &             BL*(2*g1L*muu1*muu3*muu4*(g2L*(1 - muu1**2 + muu2**2)*muu4 + g2R*muu2*(1 - muu3**2 + muu4**2)) - 
     &                g1R*(2*g2R*muu2*muu4*(-(muu2**2*muu3**2) + muu1**2*(-1 + muu3**2) + tred) + 
     &                   g2L*(-2*muu2**4*muu3**2 + muu1**2*(-1 + 2*muu3**2)*(muu2**2 + muu4**2 - tred) + (muu3**2 - tred)*(-muu4**2 + tred) + 
     &                      muu2**2*(muu3**2 - 2*muu3**4 + tred + 2*muu3**2*tred)))))*dconjg(AR)) - 
     &       2*xd1*((BL*(g1L*muu3*(g2L*muu2*(3*muu3**4 + muu3**2*(-3 + 3*muu2**2 + muu4**2 - 3*tred) + muu4**2*(muu1**2 - tred)) + 
     &                   g2R*muu4*(-2*muu1**4 - muu1**2*(-2 + muu3**2 + muu4**2 - 4*tred) + (-3 + 2*muu4**2 - 2*tred)*tred + 
     &                      muu3**2*(1 + 2*muu2**2 - 2*muu4**2 + 2*tred))) + 
     &                g1R*muu1*(g2L*muu2*(-2*muu3**4 + muu3**2*(1 - 3*muu4**2) + (-1 + muu4**2)**2) + 
     &                   g2R*muu4*(muu1**2*(-1 + 4*muu3**2 + muu4**2) + (-1 + muu4**2)*(-1 + muu4**2 - tred) - muu3**2*(4 + muu2**2 - 3*muu4**2 + 3*tred)))
     &                ) - BR*(g1R*muu1*muu3*(-1 + muu3**2 + muu4**2)*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &                g1L*(2*g2L*muu2*muu4*(-(muu2**2*muu3**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + tred - muu4**2*tred) + 
     &                   g2R*(-(muu2**4*muu3**2) + muu3**4*muu4**2 + muu3**2*muu4**4 + 2*muu3**2*tred - muu3**4*tred + muu4**2*tred - 
     &                      2*muu3**2*muu4**2*tred - muu4**4*tred - tred**2 + muu3**2*tred**2 + muu4**2*tred**2 - 
     &                      muu2**2*(muu3**4 + (-1 + muu4**2)*tred) + 
     &                      muu1**2*(muu4**4 + muu2**2*(-1 + 2*muu3**2 + muu4**2) + tred - 2*muu3**2*tred - muu4**2*(1 + tred))))))*dconjg(AL) + 
     &          (BR*(g1R*muu3*(g2R*muu2*(3*muu3**4 + muu3**2*(-3 + 3*muu2**2 + muu4**2 - 3*tred) + muu4**2*(muu1**2 - tred)) + 
     &                   g2L*muu4*(-2*muu1**4 - muu1**2*(-2 + muu3**2 + muu4**2 - 4*tred) + (-3 + 2*muu4**2 - 2*tred)*tred + 
     &                      muu3**2*(1 + 2*muu2**2 - 2*muu4**2 + 2*tred))) + 
     &                g1L*muu1*(g2R*muu2*(-2*muu3**4 + muu3**2*(1 - 3*muu4**2) + (-1 + muu4**2)**2) + 
     &                   g2L*muu4*(muu1**2*(-1 + 4*muu3**2 + muu4**2) + (-1 + muu4**2)*(-1 + muu4**2 - tred) - muu3**2*(4 + muu2**2 - 3*muu4**2 + 3*tred)))
     &                ) - BL*(g1L*muu1*muu3*(-1 + muu3**2 + muu4**2)*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &                g1R*(2*g2R*muu2*muu4*(-(muu2**2*muu3**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + tred - muu4**2*tred) + 
     &                   g2L*(-(muu2**4*muu3**2) + muu3**4*muu4**2 + muu3**2*muu4**4 + 2*muu3**2*tred - muu3**4*tred + muu4**2*tred - 
     &                      2*muu3**2*muu4**2*tred - muu4**4*tred - tred**2 + muu3**2*tred**2 + muu4**2*tred**2 - 
     &                      muu2**2*(muu3**4 + (-1 + muu4**2)*tred) + 
     &                      muu1**2*(muu4**4 + muu2**2*(-1 + 2*muu3**2 + muu4**2) + tred - 2*muu3**2*tred - muu4**2*(1 + tred))))))*dconjg(AR)) - 
     &       2*xd23*((BL*(g1L*muu3*(g2R*(-3*muu1**4 + (-1 + muu2**2)**2 + 2*muu1**2*(1 + muu2**2))*muu4 - 
     &                   g2L*muu2*(-1 - muu1**2 + muu2**2)*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                g1R*muu1*(g2R*(1 + muu1**2 - muu2**2)*muu4*(-1 + muu1**2 + muu4**2 - tred) + 
     &                   g2L*muu2*(1 + muu3**2 - 2*muu3**4 - muu4**2 + 2*muu3**2*muu4**2 + muu2**2*(-1 - muu3**2 + muu4**2) + 
     &                      muu1**2*(1 - muu3**2 + muu4**2) + 2*muu3**2*tred - 2*muu4**2*tred))) + 
     &             BR*(g1R*muu1*(1 + muu1**2 - muu2**2)*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &                g1L*(g2L*muu2*(1 + muu1**2 - muu2**2)*muu4*(muu1**2 + muu3**2 - tred) + 
     &                   g2R*(2*muu1**4*muu4**2 + (muu3**2 - tred)*(muu4**2 + muu2**2*(muu3**2 - muu4**2) - tred) + 
     &                      muu1**2*(muu4**4 + muu2**2*(1 - 2*muu4**2) - 2*tred - muu4**2*tred + muu3**2*(-muu4**2 + tred))))))*dconjg(AL) + 
     &          (BR*(g1R*muu3*(g2L*(-3*muu1**4 + (-1 + muu2**2)**2 + 2*muu1**2*(1 + muu2**2))*muu4 - 
     &                   g2R*muu2*(-1 - muu1**2 + muu2**2)*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                g1L*muu1*(g2L*(1 + muu1**2 - muu2**2)*muu4*(-1 + muu1**2 + muu4**2 - tred) + 
     &                   g2R*muu2*(1 + muu3**2 - 2*muu3**4 - muu4**2 + 2*muu3**2*muu4**2 + muu2**2*(-1 - muu3**2 + muu4**2) + 
     &                      muu1**2*(1 - muu3**2 + muu4**2) + 2*muu3**2*tred - 2*muu4**2*tred))) + 
     &             BL*(g1L*muu1*(1 + muu1**2 - muu2**2)*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &                g1R*(g2R*muu2*(1 + muu1**2 - muu2**2)*muu4*(muu1**2 + muu3**2 - tred) + 
     &                   g2L*(2*muu1**4*muu4**2 + (muu3**2 - tred)*(muu4**2 + muu2**2*(muu3**2 - muu4**2) - tred) + 
     &                      muu1**2*(muu4**4 + muu2**2*(1 - 2*muu4**2) - 2*tred - muu4**2*tred + muu3**2*(-muu4**2 + tred))))))*dconjg(AR)) - 
     &       2*xd22*((-(BL*(g1L*muu3*(g2R*(-1 + muu1**4 - 2*muu1**2*muu2**2 + muu2**4)*muu4 - g2L*muu2*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                  g1R*muu1*(g2L*muu2*(-1 + muu3**4 - 2*muu3**2*muu4**2 + muu4**4) - g2R*muu4*(-1 + muu1**2 + muu4**2 - tred)))) + 
     &             BR*(g1R*muu1*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &                g1L*(g2L*muu2*muu4*(muu1**2 + muu3**2 - tred) + 
     &                   g2R*(muu2**4*muu3**2 + muu1**4*muu4**2 + (muu3**2 - tred)*(muu4**2 - tred) - 
     &                      muu1**2*(-muu4**4 + muu2**2*(-1 + muu3**2 + muu4**2) + muu3**2*(muu4**2 - tred) + tred + muu4**2*tred) + 
     &                      muu2**2*(muu3**4 + (-1 + muu4**2)*tred - muu3**2*(muu4**2 + tred))))))*dconjg(AL) + 
     &          (-(BR*(g1R*muu3*(g2L*(-1 + muu1**4 - 2*muu1**2*muu2**2 + muu2**4)*muu4 - g2R*muu2*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                  g1L*muu1*(g2R*muu2*(-1 + muu3**4 - 2*muu3**2*muu4**2 + muu4**4) - g2L*muu4*(-1 + muu1**2 + muu4**2 - tred)))) + 
     &             BL*(g1L*muu1*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &                g1R*(g2R*muu2*muu4*(muu1**2 + muu3**2 - tred) + 
     &                   g2L*(muu2**4*muu3**2 + muu1**4*muu4**2 + (muu3**2 - tred)*(muu4**2 - tred) - 
     &                      muu1**2*(-muu4**4 + muu2**2*(-1 + muu3**2 + muu4**2) + muu3**2*(muu4**2 - tred) + tred + muu4**2*tred) + 
     &                      muu2**2*(muu3**4 + (-1 + muu4**2)*tred - muu3**2*(muu4**2 + tred))))))*dconjg(AR)) - 
     &       2*xd12*((BL*(g1R*muu1*(g2L*muu2*(-3*muu3**4 + (-1 + muu4**2)**2 + 2*muu3**2*(1 + muu4**2)) + 
     &                   g2R*muu4*(1 + muu3**2 - muu4**2)*(-1 + muu1**2 + muu4**2 - tred)) + 
     &                g1L*muu3*(g2L*muu2*(1 + muu3**2 - muu4**2)*(-1 + muu2**2 + muu3**2 - tred) + 
     &                   g2R*muu4*(1 - 2*muu1**4 + muu3**2 - muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2 - 2*tred) + 
     &                      muu1**2*(1 + 2*muu2**2 - muu3**2 - muu4**2 + 2*tred)))) + 
     &             BR*(g1R*muu1*muu3*(1 + muu3**2 - muu4**2)*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &                g1L*(g2L*muu2*muu4*(1 + muu3**2 - muu4**2)*(muu1**2 + muu3**2 - tred) + 
     &                   g2R*(muu2**4*muu3**2 + muu1**4*muu4**2 + muu3**2*muu4**2 - 2*muu3**2*tred + tred**2 - 
     &                      muu1**2*(muu2**2*(-1 + muu3**2 + muu4**2) + (1 - muu3**2 + muu4**2)*tred) + 
     &                      muu2**2*(2*muu3**4 + (-1 + muu4**2)*tred - muu3**2*(2*muu4**2 + tred))))))*dconjg(AL) + 
     &          (BR*(g1L*muu1*(g2R*muu2*(-3*muu3**4 + (-1 + muu4**2)**2 + 2*muu3**2*(1 + muu4**2)) + 
     &                   g2L*muu4*(1 + muu3**2 - muu4**2)*(-1 + muu1**2 + muu4**2 - tred)) + 
     &                g1R*muu3*(g2R*muu2*(1 + muu3**2 - muu4**2)*(-1 + muu2**2 + muu3**2 - tred) + 
     &                   g2L*muu4*(1 - 2*muu1**4 + muu3**2 - muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2 - 2*tred) + 
     &                      muu1**2*(1 + 2*muu2**2 - muu3**2 - muu4**2 + 2*tred)))) + 
     &             BL*(g1L*muu1*muu3*(1 + muu3**2 - muu4**2)*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &                g1R*(g2R*muu2*muu4*(1 + muu3**2 - muu4**2)*(muu1**2 + muu3**2 - tred) + 
     &                   g2L*(muu2**4*muu3**2 + muu1**4*muu4**2 + muu3**2*muu4**2 - 2*muu3**2*tred + tred**2 - 
     &                      muu1**2*(muu2**2*(-1 + muu3**2 + muu4**2) + (1 - muu3**2 + muu4**2)*tred) + 
     &                      muu2**2*(2*muu3**4 + (-1 + muu4**2)*tred - muu3**2*(2*muu4**2 + tred))))))*dconjg(AR))))/(4.*DenSmV1*Pi**2)
	
	
	tempBoxZ = tempBoxZ + ChiChi2QQ_M2BoxGluon - tempBox2
	tempBox2 = ChiChi2QQ_M2BoxGluon
	
! Hboson contr.

	do i= imin,imax
	
	ChiChi2QQ_M2BoxGluon = ChiChi2QQ_M2BoxGluon -(CF*g3L*g4L*NC*sqrtS**6*(muu1*xd33*(dconjg(CR(i))*((g1L*
     &                 (2*g2L*muu1**2*muu2*(-1 + muu3**2 + muu4**2) + 
     &                   g2R*muu4*(-(muu1**2*(muu2**2 + muu3**2 - muu4**2)) - (-1 + muu2**2)*(muu3**2 - tred))) + 
     &                g1R*muu1*muu3*(2*g2R*(-1 + muu1**2 + muu2**2)*muu4 - g2L*muu2*(-1 + muu2**2 + muu3**2 - tred)))*dconjg(DL(i)) + 
     &             (g1L*muu3*(4*g2L*muu1**2*muu2*muu4 + g2R*(muu1**2*(1 - muu3**2 + muu4**2) - (-1 + muu2**2)*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                g1R*muu1*(g2R*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2) - g2L*muu2*muu4*(muu1**2 + muu3**2 - tred)))*dconjg(DR(i))) + 
     &          dconjg(CL(i))*((g1R*muu3*(4*g2R*muu1**2*muu2*muu4 + 
     &                   g2L*(muu1**2*(1 - muu3**2 + muu4**2) - (-1 + muu2**2)*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                g1L*muu1*(g2L*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2) - g2R*muu2*muu4*(muu1**2 + muu3**2 - tred)))*dconjg(DL(i)) + 
     &             (g1R*(2*g2R*muu1**2*muu2*(-1 + muu3**2 + muu4**2) + 
     &                   g2L*muu4*(-(muu1**2*(muu2**2 + muu3**2 - muu4**2)) - (-1 + muu2**2)*(muu3**2 - tred))) + 
     &                g1L*muu1*muu3*(2*g2L*(-1 + muu1**2 + muu2**2)*muu4 - g2R*muu2*(-1 + muu2**2 + muu3**2 - tred)))*dconjg(DR(i)))) + 
     &       (2*xd00*(dconjg(CL(i))*((g1L*(2*g2L*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2) + g2R*muu2*muu4*(muu1**2 + muu3**2 - tred)) + 
     &                  g1R*muu1*muu3*(8*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)))*dconjg(DL(i)) + 
     &               (g1L*muu3*(4*g2L*(-1 + muu1**2 + muu2**2)*muu4 + g2R*muu2*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                  g1R*muu1*(4*g2R*muu2*(-1 + muu3**2 + muu4**2) + g2L*muu4*(-1 + muu1**2 + muu4**2 - tred)))*dconjg(DR(i))) + 
     &            dconjg(CR(i))*((g1R*muu3*(4*g2R*(-1 + muu1**2 + muu2**2)*muu4 + g2L*muu2*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                  g1L*muu1*(4*g2L*muu2*(-1 + muu3**2 + muu4**2) + g2R*muu4*(-1 + muu1**2 + muu4**2 - tred)))*dconjg(DL(i)) + 
     &               (g1R*(2*g2R*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2) + g2L*muu2*muu4*(muu1**2 + muu3**2 - tred)) + 
     &                  g1L*muu1*muu3*(8*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)))*dconjg(DR(i)))))/sqrtS**2 + 
     &       xd12*(dconjg(CL(i))*((g1L*(2*g2R*muu2*(-1 - muu1**2 + muu2**2)*muu3**2*muu4 + g2L*(-1 + muu1**2 + muu2**2)*(muu3**4 - (-1 + muu4**2)**2)) + 
     &                2*g1R*muu1*muu3*(2*g2R*muu2*muu4*(1 + muu3**2 - muu4**2) + 
     &                   g2L*(-1 + (1 + muu2**2)*muu3**2 - muu1**2*(-1 + muu3**2) + muu4**2 - tred)))*dconjg(DL(i)) + 
     &             2*(g1R*muu1*(g2L*(1 - muu1**2 + muu2**2)*muu3**2*muu4 + g2R*muu2*(muu3**4 - (-1 + muu4**2)**2)) + 
     &                g1L*muu3*(g2L*(-1 + muu1**2 + muu2**2)*muu4*(1 + muu3**2 - muu4**2) + g2R*muu2*(muu1**2 - muu1**2*muu3**2 + muu2**2*muu3**2 - tred)))
     &               *dconjg(DR(i))) + dconjg(CR(i))*(2*(g1L*muu1*(g2R*(1 - muu1**2 + muu2**2)*muu3**2*muu4 + g2L*muu2*(muu3**4 - (-1 + muu4**2)**2)) + 
     &                g1R*muu3*(g2R*(-1 + muu1**2 + muu2**2)*muu4*(1 + muu3**2 - muu4**2) + g2L*muu2*(muu1**2 - muu1**2*muu3**2 + muu2**2*muu3**2 - tred)))
     &               *dconjg(DL(i)) + (g1R*(2*g2L*muu2*(-1 - muu1**2 + muu2**2)*muu3**2*muu4 + 
     &                   g2R*(-1 + muu1**2 + muu2**2)*(muu3**4 - (-1 + muu4**2)**2)) + 
     &                2*g1L*muu1*muu3*(2*g2L*muu2*muu4*(1 + muu3**2 - muu4**2) + 
     &                   g2R*(-1 + (1 + muu2**2)*muu3**2 - muu1**2*(-1 + muu3**2) + muu4**2 - tred)))*dconjg(DR(i)))) + 
     &       xd13*(dconjg(CL(i))*((g1L*(-4*g2R*muu1**2*muu2*muu3**2*muu4 + 
     &                   g2L*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2)*(muu1**2 + muu3**2 - tred)) + 
     &                2*g1R*muu1*muu3*(2*g2R*muu2*muu4*(muu1**2 + muu3**2 - tred) + 
     &                   g2L*(-1 + muu2**2 + muu3**2 - muu1**2*(-1 + muu3**2) + muu4**2 - muu2**2*muu4**2 - tred)))*dconjg(DL(i)) + 
     &             2*(g1L*muu3*(-(g2R*muu1**2*muu2*(-1 + muu3**2 + muu4**2)) + g2L*(-1 + muu1**2 + muu2**2)*muu4*(muu1**2 + muu3**2 - tred)) + 
     &                g1R*muu1*(-(g2L*(-1 + muu1**2 + muu2**2)*muu3**2*muu4) + g2R*muu2*(-1 + muu3**2 + muu4**2)*(muu1**2 + muu3**2 - tred)))*dconjg(DR(i))
     &             ) + dconjg(CR(i))*(2*(g1R*muu3*(-(g2L*muu1**2*muu2*(-1 + muu3**2 + muu4**2)) + 
     &                   g2R*(-1 + muu1**2 + muu2**2)*muu4*(muu1**2 + muu3**2 - tred)) + 
     &                g1L*muu1*(-(g2R*(-1 + muu1**2 + muu2**2)*muu3**2*muu4) + g2L*muu2*(-1 + muu3**2 + muu4**2)*(muu1**2 + muu3**2 - tred)))*dconjg(DL(i))
     &               + (g1R*(-4*g2L*muu1**2*muu2*muu3**2*muu4 + g2R*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2)*(muu1**2 + muu3**2 - tred)) + 
     &                2*g1L*muu1*muu3*(2*g2L*muu2*muu4*(muu1**2 + muu3**2 - tred) + 
     &                   g2R*(-1 + muu2**2 + muu3**2 - muu1**2*(-1 + muu3**2) + muu4**2 - muu2**2*muu4**2 - tred)))*dconjg(DR(i)))) + 
     &       xd23*(dconjg(CL(i))*((g1L*(2*g2R*muu1**2*muu2*muu4*(-1 - muu3**2 + muu4**2) + g2L*(muu1**4 - (-1 + muu2**2)**2)*(-1 + muu3**2 + muu4**2)) + 
     &                2*g1R*muu1*muu3*(2*g2R*muu2*(1 + muu1**2 - muu2**2)*muu4 + g2L*(-1 + muu2**2 + muu3**2 + muu1**2*(1 - muu3**2 + muu4**2) - tred)))*
     &              dconjg(DL(i)) + 2*(g1L*muu3*(g2L*(muu1**4 - (-1 + muu2**2)**2)*muu4 + g2R*muu1**2*muu2*(1 - muu3**2 + muu4**2)) + 
     &                g1R*muu1*(g2R*muu2*(1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + g2L*muu4*(muu3**2 - muu1**2*muu3**2 + muu1**2*muu4**2 - tred)))
     &               *dconjg(DR(i))) + dconjg(CR(i))*(2*(g1R*muu3*(g2R*(muu1**4 - (-1 + muu2**2)**2)*muu4 + g2L*muu1**2*muu2*(1 - muu3**2 + muu4**2)) + 
     &                g1L*muu1*(g2L*muu2*(1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + g2R*muu4*(muu3**2 - muu1**2*muu3**2 + muu1**2*muu4**2 - tred)))
     &               *dconjg(DL(i)) + (g1R*(2*g2L*muu1**2*muu2*muu4*(-1 - muu3**2 + muu4**2) + 
     &                   g2R*(muu1**4 - (-1 + muu2**2)**2)*(-1 + muu3**2 + muu4**2)) + 
     &                2*g1L*muu1*muu3*(2*g2L*muu2*(1 + muu1**2 - muu2**2)*muu4 + g2R*(-1 + muu2**2 + muu3**2 + muu1**2*(1 - muu3**2 + muu4**2) - tred)))*
     &              dconjg(DR(i)))) + xd22*(dconjg(CR(i))*((g1L*muu1*
     &                 (2*g2L*muu2*(-1 + muu3**2 + muu4**2) + g2R*muu4*
     &                    (muu3**2 - muu1**2*muu3**2 + muu1**2*muu4**2 + muu2**2*(1 + muu3**2 - muu4**2) - tred)) + 
     &                g1R*muu3*(2*g2R*(-1 + muu1**2 + muu2**2)*muu4 + 
     &                   g2L*muu2*(muu4**2 + muu2**2*(muu3**2 - muu4**2) + muu1**2*(1 - muu3**2 + muu4**2) - tred)))*dconjg(DL(i)) + 
     &             (g1R*(g2R*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   g2L*muu2*muu4*(-1 + muu4**2 + muu2**2*(1 + muu3**2 - muu4**2) + muu1**2*(-muu3**2 + muu4**2) - tred)) + 
     &                g1L*muu1*muu3*(4*g2L*muu2*muu4 + g2R*(-1 + (1 + muu2**2)*muu3**2 - muu2**2*muu4**2 + muu1**2*(1 - muu3**2 + muu4**2) - tred)))*
     &              dconjg(DR(i))) + dconjg(CL(i))*((g1L*(g2L*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   g2R*muu2*muu4*(-1 + muu4**2 + muu2**2*(1 + muu3**2 - muu4**2) + muu1**2*(-muu3**2 + muu4**2) - tred)) + 
     &                g1R*muu1*muu3*(4*g2R*muu2*muu4 + g2L*(-1 + (1 + muu2**2)*muu3**2 - muu2**2*muu4**2 + muu1**2*(1 - muu3**2 + muu4**2) - tred)))*
     &              dconjg(DL(i)) + (g1R*muu1*(2*g2R*muu2*(-1 + muu3**2 + muu4**2) + 
     &                   g2L*muu4*(muu3**2 - muu1**2*muu3**2 + muu1**2*muu4**2 + muu2**2*(1 + muu3**2 - muu4**2) - tred)) + 
     &                g1L*muu3*(2*g2L*(-1 + muu1**2 + muu2**2)*muu4 + 
     &                   g2R*muu2*(muu4**2 + muu2**2*(muu3**2 - muu4**2) + muu1**2*(1 - muu3**2 + muu4**2) - tred)))*dconjg(DR(i)))) + 
     &       muu3*xd11*(dconjg(CR(i))*((g1L*muu1*muu3*(2*g2L*muu2*(-1 + muu3**2 + muu4**2) - g2R*muu4*(-1 + muu1**2 + muu4**2 - tred)) + 
     &                g1R*(2*g2R*(-1 + muu1**2 + muu2**2)*muu3**2*muu4 - 
     &                   g2L*muu2*(-(muu2**2*muu3**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + tred - muu4**2*tred)))*dconjg(DL(i)) + 
     &             (g1L*muu1*(4*g2L*muu2*muu3**2*muu4 + g2R*((1 + muu2**2)*muu3**2 - muu1**2*(-1 + muu3**2 + muu4**2) - 
     &                      (-1 + muu4**2)*(-1 + muu4**2 - tred))) + 
     &                g1R*muu3*(g2R*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2) - g2L*muu2*muu4*(muu1**2 + muu3**2 - tred)))*dconjg(DR(i))) + 
     &          dconjg(CL(i))*((g1R*muu1*(4*g2R*muu2*muu3**2*muu4 + 
     &                   g2L*((1 + muu2**2)*muu3**2 - muu1**2*(-1 + muu3**2 + muu4**2) - (-1 + muu4**2)*(-1 + muu4**2 - tred))) + 
     &                g1L*muu3*(g2L*(-1 + muu1**2 + muu2**2)*(-1 + muu3**2 + muu4**2) - g2R*muu2*muu4*(muu1**2 + muu3**2 - tred)))*dconjg(DL(i)) + 
     &             (g1R*muu1*muu3*(2*g2R*muu2*(-1 + muu3**2 + muu4**2) - g2L*muu4*(-1 + muu1**2 + muu4**2 - tred)) + 
     &                g1L*(2*g2L*(-1 + muu1**2 + muu2**2)*muu3**2*muu4 - 
     &                   g2R*muu2*(-(muu2**2*muu3**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + tred - muu4**2*tred)))*dconjg(DR(i)))) + 
     &       xd2*(dconjg(CL(i))*((g1R*muu1*muu3*(2*g2R*muu2*muu4*(2 + muu3**2 - muu4**2) - 
     &                   g2L*(2 - 2*muu3**2 - 3*muu4**2 + muu2**2*(1 - 2*muu3**2 - 2*muu4**2) + 2*muu1**2*(-1 + muu3**2 + muu4**2) + tred)) + 
     &                g1L*(g2R*muu2*muu4*((-3 + 4*muu2**2)*muu3**2 + muu1**2*(1 - 4*muu3**2) - tred) + 
     &                   g2L*(1 - muu3**4 - 2*muu4**2 + muu4**4 + muu1**2*(-1 + muu3**4 + 2*muu4**2 - muu3**2*muu4**2) + 
     &                      muu2**2*(-1 + (2 + muu3**2)*muu4**2 - muu4**4) + muu3**2*tred - muu4**2*tred)))*dconjg(DL(i)) + 
     &             (g1R*muu1*(g2R*muu2*(-2 + muu3**2 + muu3**4 + 3*muu4**2 - muu4**4) + 
     &                   g2L*muu4*(-1 + 4*(1 + muu2**2)*muu3**2 + muu1**2*(1 - 4*muu3**2) + muu4**2 - tred)) + 
     &                g1L*muu3*(g2L*(-1 + muu1**2 + muu2**2)*muu4*(2 + muu3**2 - muu4**2) + 
     &                   g2R*muu2*(1 - muu3**2 - 2*muu4**2 - 2*muu1**2*(-1 + muu3**2 + muu4**2) + muu2**2*(-1 + 2*muu3**2 + 2*muu4**2) - tred)))*
     &              dconjg(DR(i))) + dconjg(CR(i))*((g1L*muu1*(g2L*muu2*(-2 + muu3**2 + muu3**4 + 3*muu4**2 - muu4**4) + 
     &                   g2R*muu4*(-1 + 4*(1 + muu2**2)*muu3**2 + muu1**2*(1 - 4*muu3**2) + muu4**2 - tred)) + 
     &                g1R*muu3*(g2R*(-1 + muu1**2 + muu2**2)*muu4*(2 + muu3**2 - muu4**2) + 
     &                   g2L*muu2*(1 - muu3**2 - 2*muu4**2 - 2*muu1**2*(-1 + muu3**2 + muu4**2) + muu2**2*(-1 + 2*muu3**2 + 2*muu4**2) - tred)))*
     &              dconjg(DL(i)) + (g1L*muu1*muu3*(2*g2L*muu2*muu4*(2 + muu3**2 - muu4**2) - 
     &                   g2R*(2 - 2*muu3**2 - 3*muu4**2 + muu2**2*(1 - 2*muu3**2 - 2*muu4**2) + 2*muu1**2*(-1 + muu3**2 + muu4**2) + tred)) + 
     &                g1R*(g2L*muu2*muu4*((-3 + 4*muu2**2)*muu3**2 + muu1**2*(1 - 4*muu3**2) - tred) + 
     &                   g2R*(1 - muu3**4 - 2*muu4**2 + muu4**4 + muu1**2*(-1 + muu3**4 + 2*muu4**2 - muu3**2*muu4**2) + 
     &                      muu2**2*(-1 + (2 + muu3**2)*muu4**2 - muu4**4) + muu3**2*tred - muu4**2*tred)))*dconjg(DR(i)))) + 
     &       muu3*xd0*(dconjg(CL(i))*((g1L*muu3*(g2R*muu2*muu4*(-1 - 2*muu1**2 + muu2**2 - muu3**2 + tred) + 
     &                   g2L*(1 - muu3**2 - 2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + muu2**2*(-1 + 2*muu4**2) + tred)) + 
     &                g1R*muu1*(g2R*muu2*muu4*(-1 + 3*muu3**2 + muu4**2) + 
     &                   g2L*(-1 + muu3**2 + 2*muu4**2 - muu4**4 + muu2**2*(-1 + muu3**2 + muu4**2) - muu1**2*(-1 + muu3**2 + 2*muu4**2) + muu4**2*tred)))*
     &              dconjg(DL(i)) + (g1R*muu1*muu3*(g2R*muu2*(-1 + muu3**2 + 3*muu4**2) + g2L*muu4*(2 - 2*muu1**2 + muu2**2 - muu4**2 + tred)) + 
     &                g1L*(g2L*muu4*(1 - 2*muu3**2 + muu1**2*(-1 + 2*muu3**2) - muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) + tred) + 
     &                   g2R*muu2*(1 - muu3**2 - muu4**2 - muu3**2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) - muu1**2*(-1 + muu3**2 + 2*muu4**2) + 
     &                      muu4**2*tred)))*dconjg(DR(i))) + dconjg(CR(i))*
     &           ((g1L*muu1*muu3*(g2L*muu2*(-1 + muu3**2 + 3*muu4**2) + g2R*muu4*(2 - 2*muu1**2 + muu2**2 - muu4**2 + tred)) + 
     &                g1R*(g2R*muu4*(1 - 2*muu3**2 + muu1**2*(-1 + 2*muu3**2) - muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) + tred) + 
     &                   g2L*muu2*(1 - muu3**2 - muu4**2 - muu3**2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) - muu1**2*(-1 + muu3**2 + 2*muu4**2) + 
     &                      muu4**2*tred)))*dconjg(DL(i)) + (g1R*muu3*
     &                 (g2L*muu2*muu4*(-1 - 2*muu1**2 + muu2**2 - muu3**2 + tred) + 
     &                   g2R*(1 - muu3**2 - 2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + muu2**2*(-1 + 2*muu4**2) + tred)) + 
     &                g1L*muu1*(g2L*muu2*muu4*(-1 + 3*muu3**2 + muu4**2) + 
     &                   g2R*(-1 + muu3**2 + 2*muu4**2 - muu4**4 + muu2**2*(-1 + muu3**2 + muu4**2) - muu1**2*(-1 + muu3**2 + 2*muu4**2) + muu4**2*tred)))*
     &              dconjg(DR(i)))) + xd1*(dconjg(CR(i))*((g1R*muu3*
     &                 (g2R*muu4*((-3 + 3*muu1**2 + 2*muu2**2)*muu3**2 - muu1**2*muu4**2 + tred) + 
     &                   g2L*muu2*(-muu3**4 + muu1**2*(2 - 2*muu3**2 - 3*muu4**2) + (-2 + 3*muu4**2)*tred + muu3**2*(1 + muu2**2 - 3*muu4**2 + tred))) - 
     &                g1L*muu1*(g2L*muu2*(-2*muu3**4 + muu3**2*(1 - 3*muu4**2) + (-1 + muu4**2)**2) + 
     &                   g2R*muu4*(muu1**2*(-1 + 4*muu3**2 + muu4**2) + (-1 + muu4**2)*(-1 + muu4**2 - tred) - muu3**2*(4 + muu2**2 - 3*muu4**2 + 3*tred)))
     &                )*dconjg(DL(i)) + (g1L*muu1*muu3*(g2L*muu2*muu4*(1 + 5*muu3**2 - muu4**2) + 
     &                   g2R*(muu1**2*(2 - 2*muu3**2 - 3*muu4**2) - (-2 + 3*muu4**2)*(-1 + muu4**2 - tred) + muu3**2*(2 + muu2**2 - muu4**2 + tred))) + 
     &                g1R*(g2L*muu2*muu4*(muu2**2*muu3**2 - 3*muu3**4 - muu3**2*muu4**2 - muu1**2*(-1 + 4*muu3**2 + muu4**2) - tred + 3*muu3**2*tred + 
     &                      muu4**2*tred) + g2R*(2*(-1 + muu1**2)*muu3**4 - (-1 + muu4**2)*(muu1**2*muu4**2 - tred) + 
     &                      muu3**2*(2 - 3*muu4**2 + muu1**2*(-2 + muu4**2) + muu2**2*(-1 + 2*muu4**2) + 2*tred))))*dconjg(DR(i))) + 
     &          dconjg(CL(i))*((g1R*muu1*muu3*(g2R*muu2*muu4*(1 + 5*muu3**2 - muu4**2) + 
     &                   g2L*(muu1**2*(2 - 2*muu3**2 - 3*muu4**2) - (-2 + 3*muu4**2)*(-1 + muu4**2 - tred) + muu3**2*(2 + muu2**2 - muu4**2 + tred))) + 
     &                g1L*(g2R*muu2*muu4*(muu2**2*muu3**2 - 3*muu3**4 - muu3**2*muu4**2 - muu1**2*(-1 + 4*muu3**2 + muu4**2) - tred + 3*muu3**2*tred + 
     &                      muu4**2*tred) + g2L*(2*(-1 + muu1**2)*muu3**4 - (-1 + muu4**2)*(muu1**2*muu4**2 - tred) + 
     &                      muu3**2*(2 - 3*muu4**2 + muu1**2*(-2 + muu4**2) + muu2**2*(-1 + 2*muu4**2) + 2*tred))))*dconjg(DL(i)) + 
     &             (g1L*muu3*(g2L*muu4*((-3 + 3*muu1**2 + 2*muu2**2)*muu3**2 - muu1**2*muu4**2 + tred) + 
     &                   g2R*muu2*(-muu3**4 + muu1**2*(2 - 2*muu3**2 - 3*muu4**2) + (-2 + 3*muu4**2)*tred + muu3**2*(1 + muu2**2 - 3*muu4**2 + tred))) - 
     &                g1R*muu1*(g2R*muu2*(-2*muu3**4 + muu3**2*(1 - 3*muu4**2) + (-1 + muu4**2)**2) + 
     &                   g2L*muu4*(muu1**2*(-1 + 4*muu3**2 + muu4**2) + (-1 + muu4**2)*(-1 + muu4**2 - tred) - muu3**2*(4 + muu2**2 - 3*muu4**2 + 3*tred)))
     &                )*dconjg(DR(i)))) + xd3*(dconjg(CL(i))*((-(g1R*muu1*muu3*
     &                   (g2L*(2 - 2*muu3**2 - 3*muu4**2 + 2*muu1**2*(-1 + muu3**2 + muu4**2) + muu2**2*(-2 + muu3**2 + 3*muu4**2) + tred) + 
     &                     g2R*muu2*muu4*(-1 - 3*muu1**2 + muu2**2 - 2*muu3**2 + 2*tred))) + 
     &                g1L*(-(g2R*muu1**2*muu2*muu4*(-1 + 7*muu3**2 + muu4**2)) + 
     &                   g2L*(muu1**4*(-1 + muu3**2 + 2*muu4**2) + (-1 + muu2**2)*muu3**2*(-1 + muu2**2 + muu3**2 - tred) + 
     &                      muu1**2*(1 + muu3**4 - 2*muu4**2 + muu4**4 - muu3**2*(2 + muu4**2) - muu2**2*(-1 + muu3**2 + muu4**2) - muu4**2*tred))))*
     &              dconjg(DL(i)) + (g1L*muu3*(g2R*muu1**2*muu2*(3 - 3*muu3**2 - 5*muu4**2) + 
     &                   g2L*muu4*(2*muu1**4 + muu1**2*(-2 + muu3**2 - tred) + (-1 + muu2**2)*(muu3**2 - tred))) + 
     &                g1R*muu1*(g2L*muu4*(-1 + 4*muu3**2 + muu1**2*(1 - 4*muu3**2) + muu4**2 - muu2**2*(-1 + 3*muu3**2 + muu4**2) - tred) + 
     &                   g2R*muu2*(-1 + muu3**4 + muu4**2 + muu3**2*muu4**2 - muu2**2*(-1 + muu4**2) + muu1**2*(-1 + muu3**2 + 2*muu4**2) - muu3**2*tred - 
     &                      muu4**2*tred)))*dconjg(DR(i))) + dconjg(CR(i))*
     &           ((g1R*muu3*(g2L*muu1**2*muu2*(3 - 3*muu3**2 - 5*muu4**2) + 
     &                   g2R*muu4*(2*muu1**4 + muu1**2*(-2 + muu3**2 - tred) + (-1 + muu2**2)*(muu3**2 - tred))) + 
     &                g1L*muu1*(g2R*muu4*(-1 + 4*muu3**2 + muu1**2*(1 - 4*muu3**2) + muu4**2 - muu2**2*(-1 + 3*muu3**2 + muu4**2) - tred) + 
     &                   g2L*muu2*(-1 + muu3**4 + muu4**2 + muu3**2*muu4**2 - muu2**2*(-1 + muu4**2) + muu1**2*(-1 + muu3**2 + 2*muu4**2) - muu3**2*tred - 
     &                      muu4**2*tred)))*dconjg(DL(i)) + (-(g1L*muu1*muu3*
     &                   (g2R*(2 - 2*muu3**2 - 3*muu4**2 + 2*muu1**2*(-1 + muu3**2 + muu4**2) + muu2**2*(-2 + muu3**2 + 3*muu4**2) + tred) + 
     &                     g2L*muu2*muu4*(-1 - 3*muu1**2 + muu2**2 - 2*muu3**2 + 2*tred))) + 
     &                g1R*(-(g2L*muu1**2*muu2*muu4*(-1 + 7*muu3**2 + muu4**2)) + 
     &                   g2R*(muu1**4*(-1 + muu3**2 + 2*muu4**2) + (-1 + muu2**2)*muu3**2*(-1 + muu2**2 + muu3**2 - tred) + 
     &                      muu1**2*(1 + muu3**4 - 2*muu4**2 + muu4**4 - muu3**2*(2 + muu4**2) - muu2**2*(-1 + muu3**2 + muu4**2) - muu4**2*tred))))*
     &              dconjg(DR(i))))))/(4.*Pi**2*DenSmH1(i))
	
	enddo

	tempBoxH = tempBoxH + ChiChi2QQ_M2BoxGluon - tempBox2
	tempBox2 = ChiChi2QQ_M2BoxGluon
	
! U-channel contr.

	do i= 1,2
	
	ChiChi2QQ_M2BoxGluon = ChiChi2QQ_M2BoxGluon + (CF*g3L*g4L*NC*sqrtS**4*(xd22*(g1L*(-(g2R*((-1 - muu1**2 + muu2**2)*muu4*dconjg(GL(indi,i,tt,gen)) + 
     &                 muu1*(-1 + muu3**2 - muu4**2)*dconjg(GR(indi,i,tt,gen)))*
     &               (muu2*(1 + muu3**2 - muu4**2)*dconjg(HL(indj,i,tt,gen)) + (1 - muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2L*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 
     &                  2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*
     &                   dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &         g1R*(-(g2L*(muu1*(-1 + muu3**2 - muu4**2)*dconjg(GL(indi,i,tt,gen)) + (-1 - muu1**2 + muu2**2)*muu4*dconjg(GR(indi,i,tt,gen)))*
     &               ((1 - muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + muu2*(1 + muu3**2 - muu4**2)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*(-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + 
     &                  (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  (-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))))) + 
     &      muu1*xd33*(g1L*(g2R*(-2*muu1*muu4*dconjg(GL(indi,i,tt,gen)) + (-1 + muu2**2 + muu3**2 - tred)*dconjg(GR(indi,i,tt,gen)))*
     &             (-(muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &            g2L*muu1*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                ((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*
     &                   dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &         g1R*(g2L*((-1 + muu2**2 + muu3**2 - tred)*dconjg(GL(indi,i,tt,gen)) - 2*muu1*muu4*dconjg(GR(indi,i,tt,gen)))*
     &             ((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) - muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &            g2R*muu1*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  (-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))))) + 
     &      muu3*xd11*(g1L*(g2R*(-(muu4*(muu1**2 + muu3**2 - tred)*dconjg(GL(indi,i,tt,gen))) + muu1*(-1 + muu3**2 + muu4**2)*dconjg(GR(indi,i,tt,gen)))*
     &             (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &            g2L*muu3*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                ((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*
     &                   dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &         g1R*(g2L*(muu1*(-1 + muu3**2 + muu4**2)*dconjg(GL(indi,i,tt,gen)) - muu4*(muu1**2 + muu3**2 - tred)*dconjg(GR(indi,i,tt,gen)))*
     &             ((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &            g2R*muu3*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  (-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))))) - 
     &      (2*xd00*(g1L*(2*g2L*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                  (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*((1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                     dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g2R*(muu4*dconjg(GL(indi,i,tt,gen))*(-(muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                    (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 muu1*dconjg(GR(indi,i,tt,gen))*(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) - 
     &                    muu3*(muu2**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &           g1R*(g2L*(muu1*dconjg(GL(indi,i,tt,gen))*(-(muu3*(muu2**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                    muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                 muu4*dconjg(GR(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                    muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &              2*g2R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                  (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                    (1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*dconjg(HR(indj,i,tt,gen)))))))/sqrtS**2 + 
     &      xd13*(g1R*(g2R*(muu1**2 + muu3**2 - tred)*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  (-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2L*(muu4*dconjg(GR(indi,i,tt,gen))*(-(muu3*(3*muu1**4 + muu1**2*(-3 + muu2**2 + muu3**2 + 2*muu4**2 - 3*tred) + 
     &                       (-1 + muu2**2)*(muu3**2 - tred))*dconjg(HL(indj,i,tt,gen))) + 
     &                  muu2*(muu1**4 + muu1**2*(6*muu3**2 - 2*tred) + (muu3**2 - tred)**2)*dconjg(HR(indj,i,tt,gen))) - 
     &               muu1*dconjg(GL(indi,i,tt,gen))*(-(muu3*(2 - 2*muu3**2 - 2*muu4**2 + muu3**2*muu4**2 + 
     &                       muu1**2*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred) + muu2**2*(-2 + muu3**2 + 2*muu4**2 - tred) + 2*tred - muu3**2*tred - 
     &                       muu4**2*tred + tred**2)*dconjg(HL(indj,i,tt,gen))) + 
     &                  muu2*(3*muu3**4 + muu1**2*(-1 + muu3**2 + muu4**2) + muu3**2*(-3 + 2*muu2**2 + muu4**2 - 3*tred) + tred - muu4**2*tred)*
     &                   dconjg(HR(indj,i,tt,gen))))) + g1L*(g2L*(muu1**2 + muu3**2 - tred)*
     &             (muu1*muu4*dconjg(GR(indi,i,tt,gen))*((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 
     &                  2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*
     &                   dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(muu4*dconjg(GL(indi,i,tt,gen))*(muu2*(muu1**4 + muu1**2*(6*muu3**2 - 2*tred) + (muu3**2 - tred)**2)*dconjg(HL(indj,i,tt,gen)) - 
     &                  muu3*(3*muu1**4 + muu1**2*(-3 + muu2**2 + muu3**2 + 2*muu4**2 - 3*tred) + (-1 + muu2**2)*(muu3**2 - tred))*
     &                   dconjg(HR(indj,i,tt,gen))) - muu1*dconjg(GR(indi,i,tt,gen))*
     &                (muu2*(3*muu3**4 + muu1**2*(-1 + muu3**2 + muu4**2) + muu3**2*(-3 + 2*muu2**2 + muu4**2 - 3*tred) + tred - muu4**2*tred)*
     &                   dconjg(HL(indj,i,tt,gen)) - muu3*(2 - 2*muu3**2 - 2*muu4**2 + muu3**2*muu4**2 + 
     &                     muu1**2*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred) + muu2**2*(-2 + muu3**2 + 2*muu4**2 - tred) + 2*tred - muu3**2*tred - 
     &                     muu4**2*tred + tred**2)*dconjg(HR(indj,i,tt,gen)))))) + 
     &      muu3*xd0*(g1R*(g2L*(muu1*dconjg(GL(indi,i,tt,gen))*
     &                ((-((1 + muu2**2)*muu3**2) + muu1**2*(-1 + muu3**2 + 2*muu4**2) + (-1 + 2*muu4**2)*(-1 + muu4**2 - tred))*
     &                   dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu3**2 + 3*muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &               muu4*dconjg(GR(indi,i,tt,gen))*((-muu1**4 + muu2**2*(1 + muu3**2 - muu4**2) - muu1**2*(-2 + muu3**2 + muu4**2 - 2*tred) + 
     &                     (-1 + muu3**2 - tred)*(1 - muu4**2 + tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu2*muu3*(-1 + muu1**2 + muu2**2 + 2*muu3**2 - 2*tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(dconjg(GR(indi,i,tt,gen))*(muu2*(-1 + muu3**4 + muu4**2 + muu1**2*muu4**2 + muu2**2*(1 + muu3**2 - muu4**2) - tred - muu3**2*tred)*
     &                   dconjg(HL(indj,i,tt,gen)) + muu3*(-1 + muu2**4 + muu3**2 + muu4**2 + muu2**2*(muu3**2 - muu4**2 - tred) - 
     &                     muu1**2*(-1 + muu2**2 + muu3**2 + muu4**2 - tred) - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               muu1*muu4*dconjg(GL(indi,i,tt,gen))*(-(muu2*(1 + 3*muu3**2 - muu4**2)*dconjg(HL(indj,i,tt,gen))) - 
     &                  muu3*(2 - 2*muu1**2 + muu2**2 - muu4**2 + tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &         g1L*(g2L*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*(-(muu3*(2 - 2*muu1**2 + muu2**2 - muu4**2 + tred)*dconjg(HL(indj,i,tt,gen))) - 
     &                  muu2*(1 + 3*muu3**2 - muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*(muu3*(-1 + muu2**4 + muu3**2 + muu4**2 + muu2**2*(muu3**2 - muu4**2 - tred) - 
     &                     muu1**2*(-1 + muu2**2 + muu3**2 + muu4**2 - tred) - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu2*(-1 + muu3**4 + muu4**2 + muu1**2*muu4**2 + muu2**2*(1 + muu3**2 - muu4**2) - tred - muu3**2*tred)*dconjg(HR(indj,i,tt,gen))))
     &              + g2R*(muu1*dconjg(GR(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu3**2 + 3*muu4**2)*dconjg(HL(indj,i,tt,gen))) + 
     &                  (-((1 + muu2**2)*muu3**2) + muu1**2*(-1 + muu3**2 + 2*muu4**2) + (-1 + 2*muu4**2)*(-1 + muu4**2 - tred))*dconjg(HR(indj,i,tt,gen)))
     &                 + muu4*dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(-1 + muu1**2 + muu2**2 + 2*muu3**2 - 2*tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  (-muu1**4 + muu2**2*(1 + muu3**2 - muu4**2) - muu1**2*(-2 + muu3**2 + muu4**2 - 2*tred) + 
     &                     (-1 + muu3**2 - tred)*(1 - muu4**2 + tred))*dconjg(HR(indj,i,tt,gen)))))) + 
     &      xd23*(g1L*(g2L*(1 + muu1**2 - muu2**2)*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                ((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*
     &                   dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(muu4*dconjg(GL(indi,i,tt,gen))*(muu2*(muu1**4 - (-1 + muu2**2)*(muu3**2 - tred) - 
     &                     muu1**2*(-3 + muu2**2 - 3*muu3**2 + 2*muu4**2 + tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                  (-3*muu1**4 + (-1 + muu2**2)**2 + 2*muu1**2*(1 + muu2**2))*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               muu1*dconjg(GR(indi,i,tt,gen))*(muu2*(1 + muu3**2 - 2*muu3**4 - muu4**2 + 2*muu3**2*muu4**2 + muu2**2*(-1 - muu3**2 + muu4**2) + 
     &                     muu1**2*(1 - muu3**2 + muu4**2) + 2*muu3**2*tred - 2*muu4**2*tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu3*(2 - muu2**4 - 2*muu3**2 + muu4**2 + muu1**2*(-2 + muu2**2 + 2*muu3**2 - muu4**2 - tred) + tred + 
     &                     muu2**2*(-1 - muu4**2 + tred))*dconjg(HR(indj,i,tt,gen))))) + 
     &         g1R*(g2R*(1 + muu1**2 - muu2**2)*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  (-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2L*(muu1*dconjg(GL(indi,i,tt,gen))*(muu3*(2 - muu2**4 - 2*muu3**2 + muu4**2 + muu1**2*(-2 + muu2**2 + 2*muu3**2 - muu4**2 - tred) + 
     &                     tred + muu2**2*(-1 - muu4**2 + tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu2*(1 + muu3**2 - 2*muu3**4 - muu4**2 + 2*muu3**2*muu4**2 + muu2**2*(-1 - muu3**2 + muu4**2) + muu1**2*(1 - muu3**2 + muu4**2) + 
     &                     2*muu3**2*tred - 2*muu4**2*tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               muu4*dconjg(GR(indi,i,tt,gen))*((-3*muu1**4 + (-1 + muu2**2)**2 + 2*muu1**2*(1 + muu2**2))*muu3*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu2*(muu1**4 - (-1 + muu2**2)*(muu3**2 - tred) - muu1**2*(-3 + muu2**2 - 3*muu3**2 + 2*muu4**2 + tred))*dconjg(HR(indj,i,tt,gen)))
     &               ))) + xd12*(g1R*(g2R*(1 + muu3**2 - muu4**2)*
     &             (muu1*muu4*dconjg(GL(indi,i,tt,gen))*(-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + 
     &                  (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  (-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2L*(muu1*dconjg(GL(indi,i,tt,gen))*(muu3*(2 - 2*muu3**2 + 2*muu1**2*(-1 + muu3**2) - muu4**2 + muu3**2*muu4**2 - muu4**4 - 
     &                     muu2**2*(-1 + muu3**2 + muu4**2) + tred - muu3**2*tred + muu4**2*tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu2*(-3*muu3**4 + (-1 + muu4**2)**2 + 2*muu3**2*(1 + muu4**2))*dconjg(HR(indj,i,tt,gen))) + 
     &               muu4*dconjg(GR(indi,i,tt,gen))*(muu3*(1 - 2*muu1**4 + muu3**2 - muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2 - 2*tred) + 
     &                     muu1**2*(1 + 2*muu2**2 - muu3**2 - muu4**2 + 2*tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu2*(muu3**4 + muu1**2*(1 + 3*muu3**2 - muu4**2) + (-1 + muu4**2)*tred - muu3**2*(-3 + 2*muu2**2 + muu4**2 + tred))*
     &                   dconjg(HR(indj,i,tt,gen))))) + g1L*(g2L*(1 + muu3**2 - muu4**2)*
     &             (muu1*muu4*dconjg(GR(indi,i,tt,gen))*((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 
     &                  2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*
     &                   dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(muu1*dconjg(GR(indi,i,tt,gen))*(muu2*(-3*muu3**4 + (-1 + muu4**2)**2 + 2*muu3**2*(1 + muu4**2))*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu3*(2 - 2*muu3**2 + 2*muu1**2*(-1 + muu3**2) - muu4**2 + muu3**2*muu4**2 - muu4**4 - muu2**2*(-1 + muu3**2 + muu4**2) + tred - 
     &                     muu3**2*tred + muu4**2*tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               muu4*dconjg(GL(indi,i,tt,gen))*(muu2*(muu3**4 + muu1**2*(1 + 3*muu3**2 - muu4**2) + (-1 + muu4**2)*tred - 
     &                     muu3**2*(-3 + 2*muu2**2 + muu4**2 + tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu3*(1 - 2*muu1**4 + muu3**2 - muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2 - 2*tred) + 
     &                     muu1**2*(1 + 2*muu2**2 - muu3**2 - muu4**2 + 2*tred))*dconjg(HR(indj,i,tt,gen)))))) + 
     &      xd1*(g1L*(g2L*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                ((muu1**2*(-1 + 4*muu3**2 + muu4**2) + (-1 + muu4**2)*(-1 + muu4**2 - tred) - muu3**2*(4 + muu2**2 - 3*muu4**2 + 3*tred))*
     &                   dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + 7*muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*(((-1 + muu4**2)*(-1 + muu2**2 + muu4**2 - tred)*tred + muu3**4*(2 + muu2**2 - muu4**2 + tred) - 
     &                     muu1**2*(2*muu3**4 + muu2**2*(-1 + 2*muu3**2 + muu4**2) + muu3**2*(-2 + muu4**2 - 2*tred) + tred - muu4**2*tred) + 
     &                     muu3**2*(-2 + muu2**4 + 3*muu4**2 - muu4**4 + muu2**2*(1 - 2*muu4**2) - 4*tred + 2*muu4**2*tred - tred**2))*
     &                   dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(3*muu3**4 + muu3**2*(-3 + 3*muu2**2 + muu4**2 - 3*tred) + muu4**2*(muu1**2 - tred))*
     &                   dconjg(HR(indj,i,tt,gen)))) + g2R*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                (muu2*(-2*muu3**4 + muu3**2*(1 - 3*muu4**2) + (-1 + muu4**2)**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu3*(2 - 2*muu3**2 - 4*muu4**2 + 2*muu4**4 - muu2**2*(-1 + 2*muu3**2 + muu4**2) + muu1**2*(-2 + 2*muu3**2 + 3*muu4**2) + tred - 
     &                     2*muu4**2*tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               muu4*dconjg(GL(indi,i,tt,gen))*(muu2*(muu2**2*muu3**2 + 3*muu3**4 - muu3**2*muu4**2 + muu1**2*(1 + 2*muu3**2 - muu4**2) - tred - 
     &                     3*muu3**2*tred + muu4**2*tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu3*(-2*muu1**4 - muu1**2*(-2 + muu3**2 + muu4**2 - 4*tred) + (-3 + 2*muu4**2 - 2*tred)*tred + 
     &                     muu3**2*(1 + 2*muu2**2 - 2*muu4**2 + 2*tred))*dconjg(HR(indj,i,tt,gen))))) + 
     &         g1R*(g2L*(muu1*dconjg(GL(indi,i,tt,gen))*(muu3*(2 - 2*muu3**2 - 4*muu4**2 + 2*muu4**4 - muu2**2*(-1 + 2*muu3**2 + muu4**2) + 
     &                     muu1**2*(-2 + 2*muu3**2 + 3*muu4**2) + tred - 2*muu4**2*tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu2*(-2*muu3**4 + muu3**2*(1 - 3*muu4**2) + (-1 + muu4**2)**2)*dconjg(HR(indj,i,tt,gen))) + 
     &               muu4*dconjg(GR(indi,i,tt,gen))*(muu3*(-2*muu1**4 - muu1**2*(-2 + muu3**2 + muu4**2 - 4*tred) + (-3 + 2*muu4**2 - 2*tred)*tred + 
     &                     muu3**2*(1 + 2*muu2**2 - 2*muu4**2 + 2*tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu2*(muu2**2*muu3**2 + 3*muu3**4 - muu3**2*muu4**2 + muu1**2*(1 + 2*muu3**2 - muu4**2) - tred - 3*muu3**2*tred + muu4**2*tred)*
     &                   dconjg(HR(indj,i,tt,gen)))) + g2R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                (-(muu2*muu3*(-1 + 7*muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen))) + 
     &                  (muu1**2*(-1 + 4*muu3**2 + muu4**2) + (-1 + muu4**2)*(-1 + muu4**2 - tred) - muu3**2*(4 + muu2**2 - 3*muu4**2 + 3*tred))*
     &                   dconjg(HR(indj,i,tt,gen))) + dconjg(GR(indi,i,tt,gen))*
     &                (muu2*muu3*(3*muu3**4 + muu3**2*(-3 + 3*muu2**2 + muu4**2 - 3*tred) + muu4**2*(muu1**2 - tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                  ((-1 + muu4**2)*(-1 + muu2**2 + muu4**2 - tred)*tred + muu3**4*(2 + muu2**2 - muu4**2 + tred) - 
     &                     muu1**2*(2*muu3**4 + muu2**2*(-1 + 2*muu3**2 + muu4**2) + muu3**2*(-2 + muu4**2 - 2*tred) + tred - muu4**2*tred) + 
     &                     muu3**2*(-2 + muu2**4 + 3*muu4**2 - muu4**4 + muu2**2*(1 - 2*muu4**2) - 4*tred + 2*muu4**2*tred - tred**2))*
     &                   dconjg(HR(indj,i,tt,gen)))))) + xd3*(g1L*
     &          (g2L*(-(dconjg(GL(indi,i,tt,gen))*((muu1**4*(-1 + muu3**2 + 2*muu4**2) + (-1 + muu2**2)*muu3**2*(-1 + muu2**2 + muu3**2 - tred) + 
     &                       muu1**2*(muu3**4 + (-1 + 2*muu4**2)*(-1 + muu4**2 - tred) - muu3**2*(2 + tred)))*dconjg(HL(indj,i,tt,gen)) - 
     &                    muu2*muu3*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(-2 + muu2**2 + 2*muu3**2 + 3*muu4**2 - tred))*
     &                     dconjg(HR(indj,i,tt,gen)))) + muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                ((1 + muu1**4 - 3*muu3**2 - muu4**2 + muu3**2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) + 
     &                     muu1**2*(-2 + 3*muu3**2 + muu4**2 - 2*tred) + 2*tred - muu3**2*tred - muu4**2*tred + tred**2)*dconjg(HL(indj,i,tt,gen)) - 
     &                  muu2*muu3*(-1 + 3*muu1**2 + muu2**2 + 4*muu3**2 - 4*tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*muu1*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                (muu2*(1 + 5*muu3**2 - muu4**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu3*(3 - 3*muu1**2 + muu2**2 - 2*muu4**2 + 2*tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GR(indi,i,tt,gen))*(muu2*(1 - 2*muu3**4 - muu4**2 - muu1**2*muu4**2 + muu2**2*(-1 - 2*muu3**2 + muu4**2) + tred + 
     &                     muu3**2*(1 + 2*tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu3*(-muu2**4 - 2*(-1 + muu3**2 + muu4**2 - tred) + muu1**2*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred) + 
     &                     muu2**2*(-1 - muu3**2 + 2*muu4**2 + tred))*dconjg(HR(indj,i,tt,gen))))) + 
     &         g1R*(g2L*muu1*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                (muu3*(3 - 3*muu1**2 + muu2**2 - 2*muu4**2 + 2*tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu2*(1 + 5*muu3**2 - muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &               dconjg(GL(indi,i,tt,gen))*(muu3*(-muu2**4 - 2*(-1 + muu3**2 + muu4**2 - tred) + muu1**2*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred) + 
     &                     muu2**2*(-1 - muu3**2 + 2*muu4**2 + tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu2*(1 - 2*muu3**4 - muu4**2 - muu1**2*muu4**2 + muu2**2*(-1 - 2*muu3**2 + muu4**2) + tred + muu3**2*(1 + 2*tred))*
     &                   dconjg(HR(indj,i,tt,gen)))) + g2R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                (-(muu2*muu3*(-1 + 3*muu1**2 + muu2**2 + 4*muu3**2 - 4*tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                  (1 + muu1**4 - 3*muu3**2 - muu4**2 + muu3**2*muu4**2 + muu2**2*(-1 + muu3**2 + muu4**2) + 
     &                     muu1**2*(-2 + 3*muu3**2 + muu4**2 - 2*tred) + 2*tred - muu3**2*tred - muu4**2*tred + tred**2)*dconjg(HR(indj,i,tt,gen))) - 
     &               dconjg(GR(indi,i,tt,gen))*(-(muu2*muu3*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + 
     &                       muu1**2*(-2 + muu2**2 + 2*muu3**2 + 3*muu4**2 - tred))*dconjg(HL(indj,i,tt,gen))) + 
     &                  (muu1**4*(-1 + muu3**2 + 2*muu4**2) + (-1 + muu2**2)*muu3**2*(-1 + muu2**2 + muu3**2 - tred) + 
     &                     muu1**2*(muu3**4 + (-1 + 2*muu4**2)*(-1 + muu4**2 - tred) - muu3**2*(2 + tred)))*dconjg(HR(indj,i,tt,gen)))))) + 
     &      xd2*(g1L*(g2L*(-(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                  ((muu4**2*(-1 + muu1**2 + muu4**2 - tred) + muu3**2*(3 - 3*muu1**2 + 2*muu2**2 - muu4**2 + tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                    2*muu2*muu3*(1 + 2*muu3**2 - 2*muu4**2)*dconjg(HR(indj,i,tt,gen)))) + 
     &               dconjg(GL(indi,i,tt,gen))*((muu2**4*muu3**2 - muu1**4*muu4**2 + (1 + muu3**2 - muu4**2)*(-1 + muu3**2 + muu4**2 - tred) + 
     &                     muu1**2*(1 - muu3**4 - muu4**4 + muu2**2*(-muu3**2 + muu4**2) + muu3**2*tred + muu4**2*tred) + 
     &                     muu2**2*(1 + muu3**4 + muu4**4 - muu3**2*tred - muu4**2*(2 + tred)))*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu2*muu3*(-2 + muu3**4 + 3*muu4**2 + 2*muu1**2*muu4**2 + muu2**2*(2 + muu3**2 - 3*muu4**2) - 2*tred + muu4**2*tred - 
     &                     muu3**2*(-1 + muu4**2 + tred))*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(muu1*dconjg(GR(indi,i,tt,gen))*(muu2*(1 - 2*muu3**4 - muu4**2 + muu3**2*(1 + 2*muu4**2))*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu3*(2 + muu2**2 - 2*muu3**2 - 2*muu2**2*muu3**2 + 2*muu1**2*(-1 + muu3**2) - muu4**2 + tred)*dconjg(HR(indj,i,tt,gen))) + 
     &               muu4*dconjg(GL(indi,i,tt,gen))*(muu2*((3 - 2*muu2**2)*muu3**2 + muu1**2*(1 + 2*muu3**2) - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu3*(1 - 2*muu1**4 + muu2**2*(-1 + 2*muu3**2 - 2*tred) + muu1**2*(1 + 2*muu2**2 - 2*muu4**2 + 2*tred))*dconjg(HR(indj,i,tt,gen))))
     &            ) + g1R*(g2L*(muu1*dconjg(GL(indi,i,tt,gen))*
     &                (muu3*(2 + muu2**2 - 2*muu3**2 - 2*muu2**2*muu3**2 + 2*muu1**2*(-1 + muu3**2) - muu4**2 + tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu2*(1 - 2*muu3**4 - muu4**2 + muu3**2*(1 + 2*muu4**2))*dconjg(HR(indj,i,tt,gen))) + 
     &               muu4*dconjg(GR(indi,i,tt,gen))*(muu3*(1 - 2*muu1**4 + muu2**2*(-1 + 2*muu3**2 - 2*tred) + 
     &                     muu1**2*(1 + 2*muu2**2 - 2*muu4**2 + 2*tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                  muu2*((3 - 2*muu2**2)*muu3**2 + muu1**2*(1 + 2*muu3**2) - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &            g2R*(-(muu1*muu4*dconjg(GL(indi,i,tt,gen))*(2*muu2*muu3*(1 + 2*muu3**2 - 2*muu4**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                    (muu4**2*(-1 + muu1**2 + muu4**2 - tred) + muu3**2*(3 - 3*muu1**2 + 2*muu2**2 - muu4**2 + tred))*dconjg(HR(indj,i,tt,gen)))) + 
     &               dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(-2 + muu3**4 + 3*muu4**2 + 2*muu1**2*muu4**2 + muu2**2*(2 + muu3**2 - 3*muu4**2) - 2*tred + 
     &                     muu4**2*tred - muu3**2*(-1 + muu4**2 + tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                  (muu2**4*muu3**2 - muu1**4*muu4**2 + (1 + muu3**2 - muu4**2)*(-1 + muu3**2 + muu4**2 - tred) + 
     &                     muu1**2*(1 - muu3**4 - muu4**4 + muu2**2*(-muu3**2 + muu4**2) + muu3**2*tred + muu4**2*tred) + 
     &                     muu2**2*(1 + muu3**4 + muu4**4 - muu3**2*tred - muu4**2*(2 + tred)))*dconjg(HR(indj,i,tt,gen))))))))/
     &  (4.*Pi**2*(-ured + MSf(i,itypeU,gen)**2/sqrtS**2))

	tempBoxu = tempBoxu + ChiChi2QQ_M2BoxGluon - tempBox2
	tempBox2 = ChiChi2QQ_M2BoxGluon

	enddo

! T-channel contr.

	do i= 1,2
	
	ChiChi2QQ_M2BoxGluon = ChiChi2QQ_M2BoxGluon + (CF*g3L*g4L*NC*sqrtS**4*(-(muu3**2*xd11*(dconjg(ELc(indi,i,tt,gen))*
     &            ((g1R*muu1*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &                 g1L*(g2R*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &              (g1R*muu1*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &                 g1L*(g2L*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FR(indj,i,tt,gen))) + 
     &           dconjg(ER(indi,i,tt,gen))*((g1L*muu1*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &                 g1R*(g2R*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &              (g1L*muu1*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &                 g1R*(g2L*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FR(indj,i,tt,gen))))) - 
     &      (muu1**2 + muu3**2 - tred)*xd13*(dconjg(ELc(indi,i,tt,gen))*
     &          ((g1R*muu1*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &               g1L*(g2R*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1R*muu1*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &               g1L*(g2L*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FR(indj,i,tt,gen))) + 
     &         dconjg(ER(indi,i,tt,gen))*((g1L*muu1*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &               g1R*(g2R*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1L*muu1*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &               g1R*(g2L*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FR(indj,i,tt,gen)))) - 
     &      muu1**2*xd33*(dconjg(ELc(indi,i,tt,gen))*((g1R*muu1*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &               g1L*(g2R*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1R*muu1*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &               g1L*(g2L*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FR(indj,i,tt,gen))) + 
     &         dconjg(ER(indi,i,tt,gen))*((g1L*muu1*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &               g1R*(g2R*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1L*muu1*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &               g1R*(g2L*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FR(indj,i,tt,gen)))) - 
     &      (2*xd00*(dconjg(ER(indi,i,tt,gen))*((2*g1L*muu1*muu3 + g1R*(muu1**2 + muu3**2 - tred))*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred))*
     &               dconjg(FL(indj,i,tt,gen)) + 2*(g1L*muu1*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &                 g1R*(g2L*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FR(indj,i,tt,gen))) + 
     &           dconjg(ELc(indi,i,tt,gen))*(2*(g1R*muu1*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &                 g1L*(g2R*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &              (2*g1R*muu1*muu3 + g1L*(muu1**2 + muu3**2 - tred))*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred))*dconjg(FR(indj,i,tt,gen)))))/
     &       sqrtS**2 + muu3*xd0*(dconjg(ELc(indi,i,tt,gen))*((-(g1R*muu1*(-1 + muu4**2)*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred))) + 
     &               g1L*muu3*(g2R*muu2*muu4*(2 + muu1**2 - 2*muu2**2 - muu3**2 + tred) + 
     &                  g2L*(-muu2**4 + muu4**2 - muu1**2*tred + muu2**2*(1 + muu1**2 - muu3**2 - muu4**2 + tred))))*dconjg(FL(indj,i,tt,gen)) + 
     &            (-(g1R*muu1*(-1 + muu4**2)*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred))) + 
     &               g1L*muu3*(g2L*muu2*muu4*(2 + muu1**2 - 2*muu2**2 - muu3**2 + tred) + 
     &                  g2R*(-muu2**4 + muu4**2 - muu1**2*tred + muu2**2*(1 + muu1**2 - muu3**2 - muu4**2 + tred))))*dconjg(FR(indj,i,tt,gen))) + 
     &         dconjg(ER(indi,i,tt,gen))*((-(g1L*muu1*(-1 + muu4**2)*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred))) + 
     &               g1R*muu3*(g2R*muu2*muu4*(2 + muu1**2 - 2*muu2**2 - muu3**2 + tred) + 
     &                  g2L*(-muu2**4 + muu4**2 - muu1**2*tred + muu2**2*(1 + muu1**2 - muu3**2 - muu4**2 + tred))))*dconjg(FL(indj,i,tt,gen)) + 
     &            (-(g1L*muu1*(-1 + muu4**2)*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred))) + 
     &               g1R*muu3*(g2L*muu2*muu4*(2 + muu1**2 - 2*muu2**2 - muu3**2 + tred) + 
     &                  g2R*(-muu2**4 + muu4**2 - muu1**2*tred + muu2**2*(1 + muu1**2 - muu3**2 - muu4**2 + tred))))*dconjg(FR(indj,i,tt,gen)))) + 
     &      muu1*xd3*(dconjg(ER(indi,i,tt,gen))*((g1R*muu1*(-1 + muu3**2 + muu4**2)*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &               g1L*muu3*(2*g2R*muu2*muu4*(-1 + muu2**2 + muu3**2 - tred) + 
     &                  g2L*(muu2**4 - muu2**2*(1 + muu1**2 - 2*muu3**2 - muu4**2 + tred) + (-1 + muu1**2)*(muu4**2 + tred))))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1R*muu1*(-1 + muu3**2 + muu4**2)*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &               g1L*muu3*(2*g2L*muu2*muu4*(-1 + muu2**2 + muu3**2 - tred) + 
     &                  g2R*(muu2**4 - muu2**2*(1 + muu1**2 - 2*muu3**2 - muu4**2 + tred) + (-1 + muu1**2)*(muu4**2 + tred))))*dconjg(FR(indj,i,tt,gen)))
     &          + dconjg(ELc(indi,i,tt,gen))*((g1L*muu1*(-1 + muu3**2 + muu4**2)*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &               g1R*muu3*(2*g2R*muu2*muu4*(-1 + muu2**2 + muu3**2 - tred) + 
     &                  g2L*(muu2**4 - muu2**2*(1 + muu1**2 - 2*muu3**2 - muu4**2 + tred) + (-1 + muu1**2)*(muu4**2 + tred))))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1L*muu1*(-1 + muu3**2 + muu4**2)*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &               g1R*muu3*(2*g2L*muu2*muu4*(-1 + muu2**2 + muu3**2 - tred) + 
     &                  g2R*(muu2**4 - muu2**2*(1 + muu1**2 - 2*muu3**2 - muu4**2 + tred) + (-1 + muu1**2)*(muu4**2 + tred))))*dconjg(FR(indj,i,tt,gen))))
     &       + xd2*(dconjg(ER(indi,i,tt,gen))*((2*g1L*muu1*muu3*muu4*(g2L*(-1 + muu1**2 - muu2**2)*muu4 + g2R*muu2*(-1 + muu3**2 - muu4**2)) + 
     &               g1R*(2*g2R*muu2*muu4*(-(muu2**2*muu3**2) + muu1**2*(-1 + muu3**2) + tred) + 
     &                  g2L*(-2*muu2**4*muu3**2 + muu1**2*(-1 + 2*muu3**2)*(muu2**2 + muu4**2 - tred) + (muu3**2 - tred)*(-muu4**2 + tred) + 
     &                     muu2**2*(muu3**2 - 2*muu3**4 + tred + 2*muu3**2*tred))))*dconjg(FL(indj,i,tt,gen)) + 
     &            (2*g1L*muu1*muu3*muu4*(g2R*(-1 + muu1**2 - muu2**2)*muu4 + g2L*muu2*(-1 + muu3**2 - muu4**2)) + 
     &               g1R*(2*g2L*muu2*muu4*(-(muu2**2*muu3**2) + muu1**2*(-1 + muu3**2) + tred) + 
     &                  g2R*(-(muu2**4*muu3**2) + muu1**4*muu4**2 - muu2**2*(muu3**2 + muu4**2)*(muu3**2 - tred) + tred + 
     &                     muu1**2*(muu4**4 + muu2**2*(muu3**2 - muu4**2) + muu4**2*(-2 + muu3**2 - tred) - muu3**2*tred))))*dconjg(FR(indj,i,tt,gen))) + 
     &         dconjg(ELc(indi,i,tt,gen))*((2*g1R*muu1*muu3*muu4*(g2L*(-1 + muu1**2 - muu2**2)*muu4 + g2R*muu2*(-1 + muu3**2 - muu4**2)) + 
     &               g1L*(2*g2R*muu2*muu4*(-(muu2**2*muu3**2) + muu1**2*(-1 + muu3**2) + tred) + 
     &                  g2L*(-(muu2**4*muu3**2) + muu1**4*muu4**2 - muu2**2*(muu3**2 + muu4**2)*(muu3**2 - tred) + tred + 
     &                     muu1**2*(muu4**4 + muu2**2*(muu3**2 - muu4**2) + muu4**2*(-2 + muu3**2 - tred) - muu3**2*tred))))*dconjg(FL(indj,i,tt,gen)) + 
     &            (2*g1R*muu1*muu3*muu4*(g2R*(-1 + muu1**2 - muu2**2)*muu4 + g2L*muu2*(-1 + muu3**2 - muu4**2)) + 
     &               g1L*(2*g2L*muu2*muu4*(-(muu2**2*muu3**2) + muu1**2*(-1 + muu3**2) + tred) + 
     &                  g2R*(-2*muu2**4*muu3**2 + muu1**2*(-1 + 2*muu3**2)*(muu2**2 + muu4**2 - tred) + (muu3**2 - tred)*(-muu4**2 + tred) + 
     &                     muu2**2*(muu3**2 - 2*muu3**4 + tred + 2*muu3**2*tred))))*dconjg(FR(indj,i,tt,gen)))) + 
     &      xd1*(dconjg(ELc(indi,i,tt,gen))*((g1R*muu1*muu3*(-1 + muu3**2 + muu4**2)*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &               g1L*(2*g2R*muu2*muu4*(-(muu2**2*muu3**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + tred - muu4**2*tred) + 
     &                  g2L*(-(muu2**4*muu3**2) + muu3**4*muu4**2 + muu3**2*muu4**4 + 2*muu3**2*tred - muu3**4*tred + muu4**2*tred - 
     &                     2*muu3**2*muu4**2*tred - muu4**4*tred - tred**2 + muu3**2*tred**2 + muu4**2*tred**2 - muu2**2*(muu3**4 + (-1 + muu4**2)*tred) + 
     &                     muu1**2*(muu4**4 + muu2**2*(-1 + 2*muu3**2 + muu4**2) + tred - 2*muu3**2*tred - muu4**2*(1 + tred)))))*dconjg(FL(indj,i,tt,gen))
     &              + (g1R*muu1*muu3*(-1 + muu3**2 + muu4**2)*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &               g1L*(2*g2L*muu2*muu4*(-(muu2**2*muu3**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + tred - muu4**2*tred) + 
     &                  g2R*(-(muu2**4*muu3**2) + muu3**4*muu4**2 + muu3**2*muu4**4 + 2*muu3**2*tred - muu3**4*tred + muu4**2*tred - 
     &                     2*muu3**2*muu4**2*tred - muu4**4*tred - tred**2 + muu3**2*tred**2 + muu4**2*tred**2 - muu2**2*(muu3**4 + (-1 + muu4**2)*tred) + 
     &                     muu1**2*(muu4**4 + muu2**2*(-1 + 2*muu3**2 + muu4**2) + tred - 2*muu3**2*tred - muu4**2*(1 + tred)))))*dconjg(FR(indj,i,tt,gen))
     &            ) + dconjg(ER(indi,i,tt,gen))*((g1L*muu1*muu3*(-1 + muu3**2 + muu4**2)*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &               g1R*(2*g2R*muu2*muu4*(-(muu2**2*muu3**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + tred - muu4**2*tred) + 
     &                  g2L*(-(muu2**4*muu3**2) + muu3**4*muu4**2 + muu3**2*muu4**4 + 2*muu3**2*tred - muu3**4*tred + muu4**2*tred - 
     &                     2*muu3**2*muu4**2*tred - muu4**4*tred - tred**2 + muu3**2*tred**2 + muu4**2*tred**2 - muu2**2*(muu3**4 + (-1 + muu4**2)*tred) + 
     &                     muu1**2*(muu4**4 + muu2**2*(-1 + 2*muu3**2 + muu4**2) + tred - 2*muu3**2*tred - muu4**2*(1 + tred)))))*dconjg(FL(indj,i,tt,gen))
     &              + (g1L*muu1*muu3*(-1 + muu3**2 + muu4**2)*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &               g1R*(2*g2L*muu2*muu4*(-(muu2**2*muu3**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + tred - muu4**2*tred) + 
     &                  g2R*(-(muu2**4*muu3**2) + muu3**4*muu4**2 + muu3**2*muu4**4 + 2*muu3**2*tred - muu3**4*tred + muu4**2*tred - 
     &                     2*muu3**2*muu4**2*tred - muu4**4*tred - tred**2 + muu3**2*tred**2 + muu4**2*tred**2 - muu2**2*(muu3**4 + (-1 + muu4**2)*tred) + 
     &                     muu1**2*(muu4**4 + muu2**2*(-1 + 2*muu3**2 + muu4**2) + tred - 2*muu3**2*tred - muu4**2*(1 + tred)))))*dconjg(FR(indj,i,tt,gen))
     &            )) + xd23*(-(dconjg(ER(indi,i,tt,gen))*((g1L*muu1*(1 + muu1**2 - muu2**2)*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &                 g1R*(g2R*muu2*(1 + muu1**2 - muu2**2)*muu4*(muu1**2 + muu3**2 - tred) + 
     &                    g2L*(2*muu1**4*muu4**2 + (muu3**2 - tred)*(muu4**2 + muu2**2*(muu3**2 - muu4**2) - tred) + 
     &                       muu1**2*(muu4**4 + muu2**2*(1 - 2*muu4**2) - 2*tred - muu4**2*tred + muu3**2*(-muu4**2 + tred)))))*dconjg(FL(indj,i,tt,gen))
     &               + (1 + muu1**2 - muu2**2)*(g1L*muu1*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &                 g1R*(g2L*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FR(indj,i,tt,gen)))) - 
     &         dconjg(ELc(indi,i,tt,gen))*((1 + muu1**2 - muu2**2)*
     &             (g1R*muu1*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &               g1L*(g2R*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1R*muu1*(1 + muu1**2 - muu2**2)*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &               g1L*(g2L*muu2*(1 + muu1**2 - muu2**2)*muu4*(muu1**2 + muu3**2 - tred) + 
     &                  g2R*(2*muu1**4*muu4**2 + (muu3**2 - tred)*(muu4**2 + muu2**2*(muu3**2 - muu4**2) - tred) + 
     &                     muu1**2*(muu4**4 + muu2**2*(1 - 2*muu4**2) - 2*tred - muu4**2*tred + muu3**2*(-muu4**2 + tred)))))*dconjg(FR(indj,i,tt,gen))))
     &       + xd22*(-(dconjg(ER(indi,i,tt,gen))*((g1L*muu1*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &                 g1R*(g2R*muu2*muu4*(muu1**2 + muu3**2 - tred) + 
     &                    g2L*(muu2**4*muu3**2 + muu1**4*muu4**2 + (muu3**2 - tred)*(muu4**2 - tred) - 
     &                       muu1**2*(-muu4**4 + muu2**2*(-1 + muu3**2 + muu4**2) + muu3**2*(muu4**2 - tred) + tred + muu4**2*tred) + 
     &                       muu2**2*(muu3**4 + (-1 + muu4**2)*tred - muu3**2*(muu4**2 + tred)))))*dconjg(FL(indj,i,tt,gen)) + 
     &              (g1L*muu1*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &                 g1R*(g2L*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FR(indj,i,tt,gen)))) - 
     &         dconjg(ELc(indi,i,tt,gen))*((g1R*muu1*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &               g1L*(g2R*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1R*muu1*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &               g1L*(g2L*muu2*muu4*(muu1**2 + muu3**2 - tred) + 
     &                  g2R*(muu2**4*muu3**2 + muu1**4*muu4**2 + (muu3**2 - tred)*(muu4**2 - tred) - 
     &                     muu1**2*(-muu4**4 + muu2**2*(-1 + muu3**2 + muu4**2) + muu3**2*(muu4**2 - tred) + tred + muu4**2*tred) + 
     &                     muu2**2*(muu3**4 + (-1 + muu4**2)*tred - muu3**2*(muu4**2 + tred)))))*dconjg(FR(indj,i,tt,gen)))) + 
     &      xd12*(-(dconjg(ER(indi,i,tt,gen))*((g1L*muu1*muu3*(1 + muu3**2 - muu4**2)*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &                 g1R*(g2R*muu2*muu4*(1 + muu3**2 - muu4**2)*(muu1**2 + muu3**2 - tred) + 
     &                    g2L*(muu2**4*muu3**2 + muu1**4*muu4**2 + muu3**2*muu4**2 - 2*muu3**2*tred + tred**2 - 
     &                       muu1**2*(muu2**2*(-1 + muu3**2 + muu4**2) + (1 - muu3**2 + muu4**2)*tred) + 
     &                       muu2**2*(2*muu3**4 + (-1 + muu4**2)*tred - muu3**2*(2*muu4**2 + tred)))))*dconjg(FL(indj,i,tt,gen)) + 
     &              (1 + muu3**2 - muu4**2)*(g1L*muu1*muu3*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &                 g1R*(g2L*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FR(indj,i,tt,gen)))) - 
     &         dconjg(ELc(indi,i,tt,gen))*((1 + muu3**2 - muu4**2)*
     &             (g1R*muu1*muu3*(2*g2R*muu2*muu4 + g2L*(muu2**2 + muu4**2 - tred)) + 
     &               g1L*(g2R*muu2*muu4*(muu1**2 + muu3**2 - tred) + g2L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &            (g1R*muu1*muu3*(1 + muu3**2 - muu4**2)*(2*g2L*muu2*muu4 + g2R*(muu2**2 + muu4**2 - tred)) + 
     &               g1L*(g2L*muu2*muu4*(1 + muu3**2 - muu4**2)*(muu1**2 + muu3**2 - tred) + 
     &                  g2R*(muu2**4*muu3**2 + muu1**4*muu4**2 + muu3**2*muu4**2 - 2*muu3**2*tred + tred**2 - 
     &                     muu1**2*(muu2**2*(-1 + muu3**2 + muu4**2) + (1 - muu3**2 + muu4**2)*tred) + 
     &                     muu2**2*(2*muu3**4 + (-1 + muu4**2)*tred - muu3**2*(2*muu4**2 + tred)))))*dconjg(FR(indj,i,tt,gen))))))/
     &  (4.*Pi**2*(-tred + MSf(i,itypeT,gen)**2/sqrtS**2))
	
	enddo

	tempBoxt = tempBoxt + ChiChi2QQ_M2BoxGluon - tempBox2
	tempBox2 = ChiChi2QQ_M2BoxGluon
	
	enddo
	
	
c	write(*,*)'ChiChi T-Channel Gluon Box (Neu) = ', ChiChi2QQ_M2BoxGluon - tempBox
c	write(*,*)'ChiChi T-Channel Gluon Box (Neu) mit Z = ', tempBoxZ
c	write(*,*)'ChiChi T-Channel Gluon Box (Neu) mit H = ', tempBoxH
c	write(*,*)'ChiChi T-Channel Gluon Box (Neu) mit t = ', tempBoxt
c	write(*,*)'ChiChi T-Channel Gluon Box (Neu) mit u = ', tempBoxu

	
	end
