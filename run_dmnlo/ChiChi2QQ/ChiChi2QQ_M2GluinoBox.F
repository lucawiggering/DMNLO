************************************************

	function ChiChi2QQ_M2BoxGluino()

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"


	integer i,j,k,l,imin,imax,jmin,jmax,itypeT,itypeU,att,antitt
	double complex g1L,g1R,g2L,g2R,g3L,g3R,g4L,g4R
	double precision swap, M4Lr, M4L, tempBox, tempBox2, tempBoxZ, tempBoxH, tempBoxu, tempBoxt

	complex*16 DenSmV, DenSmH(6)
	double precision DenSmV1, DenSmH1(6)

	double precision ChiChi2QQ_M2BoxGluino
	external myDget

	ChiChi2QQ_M2BoxGluino = 0d0
	tempBoxZ = 0d0
	tempBoxH = 0d0
	tempBoxu = 0d0
	tempBoxt = 0d0
	tempBox2 = 0d0
	
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
!****** Gluino U-channel Box *******
!***********************************

	  do k = 1, 2
	    do l = 1, 2
			
		if(icha==0)			then
			g1L = HL(indi, l, tt, gen)
			g1R = HR(indi, l, tt, gen)
			g2L = GL(indj, k, tt, gen)
			g2R = GR(indj, k, tt, gen)
		else if(icha==1)	then
			g1L = ELc(indi, l, tt, gen)
			g1R = ER(indi, l, tt, gen)
			g2L = dconjg(HR(indj, k, tt, gen))
			g2R = dconjg(HL(indj, k, tt, gen))
		else if(icha==2)	then
! couplings set to E & F of the top final state		      
			if(tt==4)		then
				g1L =  EL*Mf(tt-1,gen)/(SW*sqrt2*MW*SB)*VChaC(indi,2)*USf(l,1,tt,gen) 
				g1R = -EL/SW*UCha(indi,1)*USf(l,1,tt,gen) + EL*Mf(tt,gen)/(SW*sqrt2*MW*CB)*UCha(indi,2)*USf(l,2,tt,gen) 
				g2L = -EL/SW*UChaC(indj,1)*USf(k,1,tt,gen) + EL*Mf(tt,gen)/(SW*sqrt2*MW*CB)*UChaC(indj,2)*USf(k,2,tt,gen)
				g2R =  EL*Mf(tt-1,gen)/(SW*sqrt2*MW*SB)*VCha(indj,2)*USf(k,1,tt,gen)
			else
				g1L = 0d0
				g1R = 0d0
				g2L = 0d0
				g2R = 0d0
			endif
		endif
		
		g3L = -sqrt2*GS*USf(l,1,att,gen)
		g3R = sqrt2*GS*USf(l,2,att,gen)
		g4L = sqrt2*GS*USf(k,2,tt,gen)
		g4R = -sqrt2*GS*USf(k,1,tt,gen)

		M4Lr = Mf(itypeU,gen)/sqrtS
		M4L = Mf(itypeU,gen)

		call myDget(muu4**2*sqrtS**2,muu3**2*sqrtS**2,muu2**2*sqrtS**2,muu1**2*sqrtS**2,sqrtS**2,ured*sqrtS**2,MSf2(l,att,gen),MGl2,MSf2(k,tt,gen),Mf2(itypeU,gen)) 
		
		
! Photon contr.

	if (icha.eq.2) then
		
	ChiChi2QQ_M2BoxGluino = ChiChi2QQ_M2BoxGluino + (ApL*BpL*CF*NC*sqrtS**4*(-2*muu1*(g1L*g2L*M4Lr + g1R*g2R*M4Lr + g1R*g2L*muu1 + g1L*g2R*muu1)*muu3*
     &       (g3L*g4L*mGlr + g3R*g4R*mGlr + g3R*g4L*muu3 + g3L*g4R*muu3)*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured)*xd0 - 
     &      (2*(g1L*g2R*(g3L*g4R*(muu1**4 + muu3**4 + muu3**2*(1 - 2*ured) + muu1**2*(1 + 2*muu3**2 - 2*ured) + ured**2) + 
     &              g3R*g4L*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*ured) + (1 + ured)**2 - muu3**2*(1 + 2*ured))) + 
     &           g1R*g2L*(g3R*g4L*(muu1**4 + muu3**4 + muu3**2*(1 - 2*ured) + muu1**2*(1 + 2*muu3**2 - 2*ured) + ured**2) + 
     &              g3L*g4R*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*ured) + (1 + ured)**2 - muu3**2*(1 + 2*ured))))*xd00)/sqrtS**2 - 
     &      2*muu3*(g1R*(g2R*(g3R*g4L + g3L*g4R)*M4Lr*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured) + 
     &            g2L*(g3L*(g4R*muu3*(3*muu1**4 + muu3**4 + muu1**2*(-2 + 4*muu3**2 - 4*ured) + ured - 2*muu3**2*ured + ured**2) + 
     &                  g4L*mGlr*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*ured) + ured - 2*muu3**2*ured + ured**2)) + 
     &               g3R*(g4L*muu3*(3*muu1**4 + muu3**4 + muu1**2*(-2 + 4*muu3**2 - 4*ured) + ured - 2*muu3**2*ured + ured**2) + 
     &                  g4R*mGlr*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*ured) + ured - 2*muu3**2*ured + ured**2)))) + 
     &         g1L*(g2L*(g3R*g4L + g3L*g4R)*M4Lr*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured) + 
     &            g2R*(g3L*(g4R*muu3*(3*muu1**4 + muu3**4 + muu1**2*(-2 + 4*muu3**2 - 4*ured) + ured - 2*muu3**2*ured + ured**2) + 
     &                  g4L*mGlr*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*ured) + ured - 2*muu3**2*ured + ured**2)) + 
     &               g3R*(g4L*muu3*(3*muu1**4 + muu3**4 + muu1**2*(-2 + 4*muu3**2 - 4*ured) + ured - 2*muu3**2*ured + ured**2) + 
     &                  g4R*mGlr*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*ured) + ured - 2*muu3**2*ured + ured**2)))))*xd1 - 
     &      2*(g1R*g2L + g1L*g2R)*(g3R*g4L + g3L*g4R)*muu3**2*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*ured) + ured - 2*muu3**2*ured + ured**2)*
     &       xd11 + (g1R*g2L*(g3L*g4R*(-muu1**6 - muu3**6 - 3*muu3**2*ured*(1 + ured) + ured*(1 + ured)**2 + muu3**4*(1 + 3*ured) + 
     &               muu1**4*(1 - 7*muu3**2 + 3*ured) + muu1**2*(-7*muu3**4 - 3*ured*(1 + ured) + 2*muu3**2*(1 + 5*ured))) - 
     &            g3R*g4L*(muu1**6 + muu3**6 + muu1**4*(7*muu3**2 - 3*ured) - 3*muu3**4*ured - ured**2*(1 + ured) + muu3**2*ured*(1 + 3*ured) + 
     &               muu1**2*(7*muu3**4 + ured*(1 + 3*ured) - 2*muu3**2*(2 + 5*ured)))) + 
     &         g1L*g2R*(g3R*g4L*(-muu1**6 - muu3**6 - 3*muu3**2*ured*(1 + ured) + ured*(1 + ured)**2 + muu3**4*(1 + 3*ured) + 
     &               muu1**4*(1 - 7*muu3**2 + 3*ured) + muu1**2*(-7*muu3**4 - 3*ured*(1 + ured) + 2*muu3**2*(1 + 5*ured))) - 
     &            g3L*g4R*(muu1**6 + muu3**6 + muu1**4*(7*muu3**2 - 3*ured) - 3*muu3**4*ured - ured**2*(1 + ured) + muu3**2*ured*(1 + 3*ured) + 
     &               muu1**2*(7*muu3**4 + ured*(1 + 3*ured) - 2*muu3**2*(2 + 5*ured)))))*xd13 - 
     &      2*muu1*(g1R*(g2R*(g3R*g4L + g3L*g4R)*M4Lr*(muu1**4 + muu3**4 + 2*muu1**2*(muu3**2 - ured) + ured*(1 + ured) - muu3**2*(1 + 2*ured)) + 
     &            g2L*muu1*(g3R*(g4R*mGlr*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured) + 
     &                  g4L*(muu1**4 + 3*muu3**4 + muu1**2*(4*muu3**2 - 2*ured) + ured*(1 + ured) - 2*muu3**2*(1 + 2*ured))) + 
     &               g3L*(g4L*mGlr*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured) + 
     &                  g4R*(muu1**4 + 3*muu3**4 + muu1**2*(4*muu3**2 - 2*ured) + ured*(1 + ured) - 2*muu3**2*(1 + 2*ured))))) + 
     &         g1L*(g2L*(g3R*g4L + g3L*g4R)*M4Lr*(muu1**4 + muu3**4 + 2*muu1**2*(muu3**2 - ured) + ured*(1 + ured) - muu3**2*(1 + 2*ured)) + 
     &            g2R*muu1*(g3R*(g4R*mGlr*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured) + 
     &                  g4L*(muu1**4 + 3*muu3**4 + muu1**2*(4*muu3**2 - 2*ured) + ured*(1 + ured) - 2*muu3**2*(1 + 2*ured))) + 
     &               g3L*(g4L*mGlr*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured) + 
     &                  g4R*(muu1**4 + 3*muu3**4 + muu1**2*(4*muu3**2 - 2*ured) + ured*(1 + ured) - 2*muu3**2*(1 + 2*ured))))))*xd3 - 
     &      2*(g1R*g2L + g1L*g2R)*(g3R*g4L + g3L*g4R)*muu1**2*(muu1**4 + muu3**4 + 2*muu1**2*(muu3**2 - ured) + ured*(1 + ured) - muu3**2*(1 + 2*ured))*
     &       xd33))/(4.*Pi**2)
		
	endif

! Zboson contr.

	ChiChi2QQ_M2BoxGluino = ChiChi2QQ_M2BoxGluino -(CF*NC*sqrtS**6*(xd0*((BR*(g1R*muu1*(g2R*M4Lr*(g3L*muu4*(g4R*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4L*mGlr*(muu2**2 + muu3**2 - ured)) + 
     &                      g3R*(g4R*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4L*muu4**2*(muu2**2 + muu3**2 - ured))) + 
     &                   g2L*muu2*(g3L*muu4*(g4R*muu3*(muu1**2 + muu4**2 - ured) + g4L*mGlr*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                      g3R*(g4R*mGlr*muu3*(muu1**2 + muu4**2 - ured) + g4L*muu4**2*(-1 + muu2**2 + muu4**2 - ured)))) + 
     &                g1L*(g2R*muu1**2*(g3L*muu4*(g4R*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4L*mGlr*(muu2**2 + muu3**2 - ured)) + 
     &                      g3R*(g4R*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4L*muu4**2*(muu2**2 + muu3**2 - ured))) + 
     &                   g2L*M4Lr*muu2*(g3L*muu4*(g4R*muu3*(muu1**2 + muu4**2 - ured) + g4L*mGlr*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                      g3R*(g4R*mGlr*muu3*(muu1**2 + muu4**2 - ured) + g4L*muu4**2*(-1 + muu2**2 + muu4**2 - ured))))) + 
     &             BL*(g1R*muu1*(g2R*M4Lr*(g3R*muu4*(g4L*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4R*mGlr*(muu2**2 + muu3**2 - ured)) + 
     &                      g3L*(g4L*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4R*muu4**2*(muu2**2 + muu3**2 - ured))) + 
     &                   g2L*muu2*(g3R*muu4*(g4L*muu3*(muu1**2 + muu4**2 - ured) + g4R*mGlr*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                      g3L*(g4L*mGlr*muu3*(muu1**2 + muu4**2 - ured) + g4R*muu4**2*(-1 + muu2**2 + muu4**2 - ured)))) + 
     &                g1L*(g2R*muu1**2*(g3R*muu4*(g4L*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4R*mGlr*(muu2**2 + muu3**2 - ured)) + 
     &                      g3L*(g4L*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4R*muu4**2*(muu2**2 + muu3**2 - ured))) + 
     &                   g2L*M4Lr*muu2*(g3R*muu4*(g4L*muu3*(muu1**2 + muu4**2 - ured) + g4R*mGlr*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                      g3L*(g4L*mGlr*muu3*(muu1**2 + muu4**2 - ured) + g4R*muu4**2*(-1 + muu2**2 + muu4**2 - ured))))))*dconjg(AL) + 
     &          (BR*(g1L*muu1*(g2L*M4Lr*(g3L*muu4*(g4R*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4L*mGlr*(muu2**2 + muu3**2 - ured)) + 
     &                      g3R*(g4R*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4L*muu4**2*(muu2**2 + muu3**2 - ured))) + 
     &                   g2R*muu2*(g3L*muu4*(g4R*muu3*(muu1**2 + muu4**2 - ured) + g4L*mGlr*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                      g3R*(g4R*mGlr*muu3*(muu1**2 + muu4**2 - ured) + g4L*muu4**2*(-1 + muu2**2 + muu4**2 - ured)))) + 
     &                g1R*(g2L*muu1**2*(g3L*muu4*(g4R*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4L*mGlr*(muu2**2 + muu3**2 - ured)) + 
     &                      g3R*(g4R*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4L*muu4**2*(muu2**2 + muu3**2 - ured))) + 
     &                   g2R*M4Lr*muu2*(g3L*muu4*(g4R*muu3*(muu1**2 + muu4**2 - ured) + g4L*mGlr*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                      g3R*(g4R*mGlr*muu3*(muu1**2 + muu4**2 - ured) + g4L*muu4**2*(-1 + muu2**2 + muu4**2 - ured))))) + 
     &             BL*(g1L*muu1*(g2L*M4Lr*(g3R*muu4*(g4L*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4R*mGlr*(muu2**2 + muu3**2 - ured)) + 
     &                      g3L*(g4L*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4R*muu4**2*(muu2**2 + muu3**2 - ured))) + 
     &                   g2R*muu2*(g3R*muu4*(g4L*muu3*(muu1**2 + muu4**2 - ured) + g4R*mGlr*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                      g3L*(g4L*mGlr*muu3*(muu1**2 + muu4**2 - ured) + g4R*muu4**2*(-1 + muu2**2 + muu4**2 - ured)))) + 
     &                g1R*(g2L*muu1**2*(g3R*muu4*(g4L*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4R*mGlr*(muu2**2 + muu3**2 - ured)) + 
     &                      g3L*(g4L*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4R*muu4**2*(muu2**2 + muu3**2 - ured))) + 
     &                   g2R*M4Lr*muu2*(g3R*muu4*(g4L*muu3*(muu1**2 + muu4**2 - ured) + g4R*mGlr*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                      g3L*(g4L*mGlr*muu3*(muu1**2 + muu4**2 - ured) + g4R*muu4**2*(-1 + muu2**2 + muu4**2 - ured))))))*dconjg(AR)) + 
     &       (2*xd00*((BL*(g1R*g2L*muu1*muu2*(4*g3R*g4L*muu3*muu4 - g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &                  g1L*g2R*(-(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4) + g3L*g4R*(muu2**2 + muu3**2 - ured)*(muu1**2 + muu4**2 - ured))) + 
     &               BR*(g1R*g2L*muu1*muu2*(4*g3L*g4R*muu3*muu4 - g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &                  g1L*g2R*(-(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4) + g3R*g4L*(-1 + muu1**2 + muu3**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))))
     &              *dconjg(AL) + (BR*(g1L*g2R*muu1*muu2*(4*g3L*g4R*muu3*muu4 - g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &                  g1R*g2L*(-(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4) + g3R*g4L*(muu2**2 + muu3**2 - ured)*(muu1**2 + muu4**2 - ured))) + 
     &               BL*(g1L*g2R*muu1*muu2*(4*g3R*g4L*muu3*muu4 - g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &                  g1R*g2L*(-(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4) + g3L*g4R*(-1 + muu1**2 + muu3**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))))
     &              *dconjg(AR)))/sqrtS**2 + muu1*xd33*((BR*(g1R*g2L*muu2*
     &                 (2*g3L*g4R*muu1**2*muu3*muu4 + g3R*g4L*(muu1**2*(muu2**2 - muu3**2 - ured) + (muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                g1L*g2R*muu1*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3R*g4L*(muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))) + 
     &             BL*(g1R*g2L*muu2*(2*g3R*g4L*muu1**2*muu3*muu4 + 
     &                   g3L*g4R*(muu1**2*(muu2**2 - muu3**2 - ured) + (muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                g1L*g2R*muu1*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3L*g4R*(muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))))*dconjg(AL) + 
     &          (BR*(g1L*g2R*muu2*(2*g3L*g4R*muu1**2*muu3*muu4 + 
     &                   g3R*g4L*(muu1**2*(muu2**2 - muu3**2 - ured) + (muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                g1R*g2L*muu1*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3R*g4L*(muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))) + 
     &             BL*(g1L*g2R*muu2*(2*g3R*g4L*muu1**2*muu3*muu4 + 
     &                   g3L*g4R*(muu1**2*(muu2**2 - muu3**2 - ured) + (muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                g1R*g2L*muu1*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3L*g4R*(muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))))*dconjg(AR)) + 
     &       muu4*xd11*((BR*(g1R*g2L*muu1*muu2*muu4*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &                g1L*g2R*(g3L*g4R*muu3*(muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 - 2*ured) + ured - muu3**2*ured - 
     &                      muu4**2*ured + ured**2) + g3R*g4L*muu4*
     &                    (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))) + 
     &             BL*(g1R*g2L*muu1*muu2*muu4*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &                g1L*g2R*(g3R*g4L*muu3*(muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 - 2*ured) + ured - muu3**2*ured - 
     &                      muu4**2*ured + ured**2) + g3L*g4R*muu4*
     &                    (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))))*dconjg(AL) + 
     &          (BR*(g1L*g2R*muu1*muu2*muu4*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &                g1R*g2L*(g3L*g4R*muu3*(muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 - 2*ured) + ured - muu3**2*ured - 
     &                      muu4**2*ured + ured**2) + g3R*g4L*muu4*
     &                    (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))) + 
     &             BL*(g1L*g2R*muu1*muu2*muu4*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &                g1R*g2L*(g3R*g4L*muu3*(muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 - 2*ured) + ured - muu3**2*ured - 
     &                      muu4**2*ured + ured**2) + g3L*g4R*muu4*
     &                    (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))))*dconjg(AR)) + 
     &       xd23*((BR*(2*g1R*g2L*muu1*muu2*(g3L*g4R*(1 + muu1**2 - muu2**2)*muu3*muu4 - 
     &                   g3R*g4L*(muu1**2*muu3**2 - muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured))) + 
     &                g1L*g2R*(2*g3L*g4R*muu1**2*(-1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                   g3R*g4L*(muu1**4*(muu2**2 - muu3**2 - ured) - muu2**2*(muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured) + 
     &                      muu1**2*(-muu2**4 + muu3**2 - muu3**4 + 2*muu3**2*muu4**2 + muu2**2*(muu3**2 + muu4**2) + ured - 2*muu4**2*ured + ured**2))))
     &              + BL*(-2*g1R*g2L*muu1*muu2*(g3R*g4L*(-1 - muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3L*g4R*(muu1**2*muu3**2 - muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured))) + 
     &                g1L*g2R*(2*g3R*g4L*muu1**2*(-1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                   g3L*g4R*(muu1**4*(muu2**2 - muu3**2 - ured) - muu2**2*(muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured) + 
     &                      muu1**2*(-muu2**4 + muu3**2 - muu3**4 + 2*muu3**2*muu4**2 + muu2**2*(muu3**2 + muu4**2) + ured - 2*muu4**2*ured + ured**2)))))*
     &           dconjg(AL) + (BR*(2*g1L*g2R*muu1*muu2*(g3L*g4R*(1 + muu1**2 - muu2**2)*muu3*muu4 - 
     &                   g3R*g4L*(muu1**2*muu3**2 - muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured))) + 
     &                g1R*g2L*(2*g3L*g4R*muu1**2*(-1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                   g3R*g4L*(muu1**4*(muu2**2 - muu3**2 - ured) - muu2**2*(muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured) + 
     &                      muu1**2*(-muu2**4 + muu3**2 - muu3**4 + 2*muu3**2*muu4**2 + muu2**2*(muu3**2 + muu4**2) + ured - 2*muu4**2*ured + ured**2))))
     &              + BL*(-2*g1L*g2R*muu1*muu2*(g3R*g4L*(-1 - muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3L*g4R*(muu1**2*muu3**2 - muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured))) + 
     &                g1R*g2L*(2*g3R*g4L*muu1**2*(-1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                   g3L*g4R*(muu1**4*(muu2**2 - muu3**2 - ured) - muu2**2*(muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured) + 
     &                      muu1**2*(-muu2**4 + muu3**2 - muu3**4 + 2*muu3**2*muu4**2 + muu2**2*(muu3**2 + muu4**2) + ured - 2*muu4**2*ured + ured**2)))))*
     &           dconjg(AR)) + xd22*((BR*(g1R*g2L*muu1*muu2*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(muu3**4 + muu4**2*(-1 + muu4**2) - muu3**2*(1 + 2*muu4**2))) + 
     &                g1L*g2R*(g3L*g4R*(muu1**4 + muu2**2*(-1 + muu2**2) - muu1**2*(1 + 2*muu2**2))*muu3*muu4 + 
     &                   g3R*g4L*(-(muu1**4*muu3**2) + muu2**2*(-muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured)) + 
     &                      muu1**2*(-muu3**4 + muu4**2*(muu2**2 - ured) + muu3**2*(1 + muu2**2 + muu4**2 + ured))))) + 
     &             BL*(g1R*g2L*muu1*muu2*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(muu3**4 + muu4**2*(-1 + muu4**2) - muu3**2*(1 + 2*muu4**2))) + 
     &                g1L*g2R*(g3R*g4L*(muu1**4 + muu2**2*(-1 + muu2**2) - muu1**2*(1 + 2*muu2**2))*muu3*muu4 + 
     &                   g3L*g4R*(-(muu1**4*muu3**2) + muu2**2*(-muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured)) + 
     &                      muu1**2*(-muu3**4 + muu4**2*(muu2**2 - ured) + muu3**2*(1 + muu2**2 + muu4**2 + ured))))))*dconjg(AL) + 
     &          (BR*(g1L*g2R*muu1*muu2*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(muu3**4 + muu4**2*(-1 + muu4**2) - muu3**2*(1 + 2*muu4**2))) + 
     &                g1R*g2L*(g3L*g4R*(muu1**4 + muu2**2*(-1 + muu2**2) - muu1**2*(1 + 2*muu2**2))*muu3*muu4 + 
     &                   g3R*g4L*(-(muu1**4*muu3**2) + muu2**2*(-muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured)) + 
     &                      muu1**2*(-muu3**4 + muu4**2*(muu2**2 - ured) + muu3**2*(1 + muu2**2 + muu4**2 + ured))))) + 
     &             BL*(g1L*g2R*muu1*muu2*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(muu3**4 + muu4**2*(-1 + muu4**2) - muu3**2*(1 + 2*muu4**2))) + 
     &                g1R*g2L*(g3R*g4L*(muu1**4 + muu2**2*(-1 + muu2**2) - muu1**2*(1 + 2*muu2**2))*muu3*muu4 + 
     &                   g3L*g4R*(-(muu1**4*muu3**2) + muu2**2*(-muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured)) + 
     &                      muu1**2*(-muu3**4 + muu4**2*(muu2**2 - ured) + muu3**2*(1 + muu2**2 + muu4**2 + ured))))))*dconjg(AR)) + 
     &       xd13*((BL*(2*g1R*g2L*muu1*muu2*muu4*(g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) + g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                g1L*g2R*(2*g3R*g4L*muu1**2*muu3*muu4*(-1 + muu1**2 + muu3**2 - ured) + 
     &                   g3L*g4R*(muu1**2 + muu4**2 - ured)*(muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))))
     &              + BR*(2*g1R*g2L*muu1*muu2*muu4*(g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) + g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                g1L*g2R*(2*g3L*g4R*muu1**2*muu3*muu4*(-1 + muu1**2 + muu3**2 - ured) + 
     &                   g3R*g4L*(muu1**4*(muu2**2 - muu3**2 - ured) - 
     &                      (-1 + muu2**2 + muu4**2 - ured)*(muu2**2*muu4**2 + (-1 + muu4**2 - ured)*ured + muu3**2*(-muu4**2 + ured)) + 
     &                      muu1**2*(-muu3**4 + muu3**2*(1 + 2*muu4**2) + muu2**2*(-1 + muu3**2 + 2*muu4**2 - 2*ured) + ured*(2 - 3*muu4**2 + 2*ured))))))*
     &           dconjg(AL) + (BR*(2*g1L*g2R*muu1*muu2*muu4*(g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) + g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                g1R*g2L*(2*g3L*g4R*muu1**2*muu3*muu4*(-1 + muu1**2 + muu3**2 - ured) + 
     &                   g3R*g4L*(muu1**2 + muu4**2 - ured)*(muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))))
     &              + BL*(2*g1L*g2R*muu1*muu2*muu4*(g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) + g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                g1R*g2L*(2*g3R*g4L*muu1**2*muu3*muu4*(-1 + muu1**2 + muu3**2 - ured) + 
     &                   g3L*g4R*(muu1**4*(muu2**2 - muu3**2 - ured) - 
     &                      (-1 + muu2**2 + muu4**2 - ured)*(muu2**2*muu4**2 + (-1 + muu4**2 - ured)*ured + muu3**2*(-muu4**2 + ured)) + 
     &                      muu1**2*(-muu3**4 + muu3**2*(1 + 2*muu4**2) + muu2**2*(-1 + muu3**2 + 2*muu4**2 - 2*ured) + ured*(2 - 3*muu4**2 + 2*ured))))))*
     &           dconjg(AR)) + xd12*((BR*(2*g1R*g2L*muu1*muu2*muu4*(g3R*g4L*muu4*(-1 - muu3**2 + muu4**2) + g3L*g4R*muu3*(1 - muu3**2 + muu4**2)) + 
     &                g1L*g2R*(2*g3L*g4R*muu3*muu4*(muu1**4 - muu1**2*(1 + muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                   g3R*g4L*(-(muu1**4*muu3**2) - muu2**4*muu4**2 + muu2**2*muu4**2*(1 + muu3**2 - muu4**2) - 
     &                      (muu3**2 - muu4**2)*(muu3**2*(muu4**2 - ured) + ured*(1 - muu4**2 + ured)) + 
     &                      muu1**2*(-muu3**4 + 2*muu4**2*(muu2**2 - ured) + muu3**2*(1 + muu4**2 + 2*ured))))) + 
     &             BL*(2*g1R*g2L*muu1*muu2*muu4*(g3L*g4R*muu4*(-1 - muu3**2 + muu4**2) + g3R*g4L*muu3*(1 - muu3**2 + muu4**2)) + 
     &                g1L*g2R*(2*g3R*g4L*muu3*muu4*(muu1**4 - muu1**2*(1 + muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                   g3L*g4R*(-(muu1**4*muu3**2) - muu2**4*muu4**2 + muu2**2*muu4**2*(1 + muu3**2 - muu4**2) - 
     &                      (muu3**2 - muu4**2)*(muu3**2*(muu4**2 - ured) + ured*(1 - muu4**2 + ured)) + 
     &                      muu1**2*(-muu3**4 + 2*muu4**2*(muu2**2 - ured) + muu3**2*(1 + muu4**2 + 2*ured))))))*dconjg(AL) + 
     &          (BR*(2*g1L*g2R*muu1*muu2*muu4*(g3R*g4L*muu4*(-1 - muu3**2 + muu4**2) + g3L*g4R*muu3*(1 - muu3**2 + muu4**2)) + 
     &                g1R*g2L*(2*g3L*g4R*muu3*muu4*(muu1**4 - muu1**2*(1 + muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                   g3R*g4L*(-(muu1**4*muu3**2) - muu2**4*muu4**2 + muu2**2*muu4**2*(1 + muu3**2 - muu4**2) - 
     &                      (muu3**2 - muu4**2)*(muu3**2*(muu4**2 - ured) + ured*(1 - muu4**2 + ured)) + 
     &                      muu1**2*(-muu3**4 + 2*muu4**2*(muu2**2 - ured) + muu3**2*(1 + muu4**2 + 2*ured))))) + 
     &             BL*(2*g1L*g2R*muu1*muu2*muu4*(g3L*g4R*muu4*(-1 - muu3**2 + muu4**2) + g3R*g4L*muu3*(1 - muu3**2 + muu4**2)) + 
     &                g1R*g2L*(2*g3R*g4L*muu3*muu4*(muu1**4 - muu1**2*(1 + muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                   g3L*g4R*(-(muu1**4*muu3**2) - muu2**4*muu4**2 + muu2**2*muu4**2*(1 + muu3**2 - muu4**2) - 
     &                      (muu3**2 - muu4**2)*(muu3**2*(muu4**2 - ured) + ured*(1 - muu4**2 + ured)) + 
     &                      muu1**2*(-muu3**4 + 2*muu4**2*(muu2**2 - ured) + muu3**2*(1 + muu4**2 + 2*ured))))))*dconjg(AR)) + 
     &       xd3*((BR*(g1R*muu1*(g2L*muu2*(g3R*(g4L*(muu1**2*(muu2**2 - muu3**2 - ured) + (-1 + muu2**2 + muu4**2 - ured)*(2*muu4**2 - ured)) + 
     &                         g4R*mGlr*muu3*(muu1**2 + muu4**2 - ured)) + 
     &                      g3L*muu4*(g4R*muu3*(3*muu1**2 + muu4**2 - ured) + g4L*mGlr*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                   g2R*M4Lr*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                      g3R*g4L*(muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))) + 
     &                g1L*(g2L*M4Lr*muu2*(2*g3L*g4R*muu1**2*muu3*muu4 + 
     &                      g3R*g4L*(muu1**2*(muu2**2 - muu3**2 - ured) + (muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                   g2R*muu1**2*(g3L*muu4*(g4L*mGlr*(muu2**2 + muu3**2 - ured) + g4R*muu3*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured)) + 
     &                      g3R*(g4R*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + 
     &                         g4L*(2*muu3**2*muu4**2 + muu1**2*(muu2**2 - ured) + muu2**2*(muu4**2 - ured) + ured - muu3**2*ured - 2*muu4**2*ured + 
     &                            ured**2))))) + BL*(g1R*muu1*(g2L*muu2*
     &                    (g3L*(g4R*(muu1**2*(muu2**2 - muu3**2 - ured) + (-1 + muu2**2 + muu4**2 - ured)*(2*muu4**2 - ured)) + 
     &                         g4L*mGlr*muu3*(muu1**2 + muu4**2 - ured)) + 
     &                      g3R*muu4*(g4L*muu3*(3*muu1**2 + muu4**2 - ured) + g4R*mGlr*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                   g2R*M4Lr*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                      g3L*g4R*(muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))) + 
     &                g1L*(g2L*M4Lr*muu2*(2*g3R*g4L*muu1**2*muu3*muu4 + 
     &                      g3L*g4R*(muu1**2*(muu2**2 - muu3**2 - ured) + (muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                   g2R*muu1**2*(g3R*muu4*(g4R*mGlr*(muu2**2 + muu3**2 - ured) + g4L*muu3*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured)) + 
     &                      g3L*(g4L*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + 
     &                         g4R*(2*muu3**2*muu4**2 + muu1**2*(muu2**2 - ured) + muu2**2*(muu4**2 - ured) + ured - muu3**2*ured - 2*muu4**2*ured + 
     &                            ured**2))))))*dconjg(AL) + (BR*
     &              (g1L*muu1*(g2R*muu2*(g3R*(g4L*(muu1**2*(muu2**2 - muu3**2 - ured) + (-1 + muu2**2 + muu4**2 - ured)*(2*muu4**2 - ured)) + 
     &                         g4R*mGlr*muu3*(muu1**2 + muu4**2 - ured)) + 
     &                      g3L*muu4*(g4R*muu3*(3*muu1**2 + muu4**2 - ured) + g4L*mGlr*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                   g2L*M4Lr*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                      g3R*g4L*(muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))) + 
     &                g1R*(g2R*M4Lr*muu2*(2*g3L*g4R*muu1**2*muu3*muu4 + 
     &                      g3R*g4L*(muu1**2*(muu2**2 - muu3**2 - ured) + (muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                   g2L*muu1**2*(g3L*muu4*(g4L*mGlr*(muu2**2 + muu3**2 - ured) + g4R*muu3*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured)) + 
     &                      g3R*(g4R*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + 
     &                         g4L*(2*muu3**2*muu4**2 + muu1**2*(muu2**2 - ured) + muu2**2*(muu4**2 - ured) + ured - muu3**2*ured - 2*muu4**2*ured + 
     &                            ured**2))))) + BL*(g1L*muu1*(g2R*muu2*
     &                    (g3L*(g4R*(muu1**2*(muu2**2 - muu3**2 - ured) + (-1 + muu2**2 + muu4**2 - ured)*(2*muu4**2 - ured)) + 
     &                         g4L*mGlr*muu3*(muu1**2 + muu4**2 - ured)) + 
     &                      g3R*muu4*(g4L*muu3*(3*muu1**2 + muu4**2 - ured) + g4R*mGlr*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                   g2L*M4Lr*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                      g3L*g4R*(muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))) + 
     &                g1R*(g2R*M4Lr*muu2*(2*g3R*g4L*muu1**2*muu3*muu4 + 
     &                      g3L*g4R*(muu1**2*(muu2**2 - muu3**2 - ured) + (muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                   g2L*muu1**2*(g3R*muu4*(g4R*mGlr*(muu2**2 + muu3**2 - ured) + g4L*muu3*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured)) + 
     &                      g3L*(g4L*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + 
     &                         g4R*(2*muu3**2*muu4**2 + muu1**2*(muu2**2 - ured) + muu2**2*(muu4**2 - ured) + ured - muu3**2*ured - 2*muu4**2*ured + 
     &                            ured**2))))))*dconjg(AR)) + xd1*
     &        ((BR*(g1R*muu1*muu4*(g2L*muu2*(g3R*muu4*(2*g4R*mGlr*muu3 + g4L*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)) + 
     &                      g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2) + g4R*muu3*(muu1**2 + 3*muu4**2 - ured))) + 
     &                   g2R*M4Lr*(g3L*g4R*muu3*(-1 + muu1**2 + muu3**2 - ured) + g3R*g4L*muu4*(muu2**2 + muu3**2 - ured))) + 
     &                g1L*(g2L*M4Lr*muu2*muu4*(g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) + g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                   g2R*(g3R*(g4R*mGlr*muu3*(muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 - 2*ured) + ured - muu3**2*ured - 
     &                            muu4**2*ured + ured**2) + g4L*muu4**2*
     &                          (muu1**2*(2*muu2**2 + muu3**2 - 2*ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))) + 
     &                      g3L*muu4*(g4R*muu3*(2*muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-2 + 2*muu3**2 - 3*ured) + ured - muu3**2*ured - 
     &                            muu4**2*ured + ured**2) + g4L*mGlr*
     &                          (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))))) + 
     &             BL*(g1R*muu1*muu4*(g2L*muu2*(g3L*muu4*(2*g4L*mGlr*muu3 + g4R*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)) + 
     &                      g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2) + g4L*muu3*(muu1**2 + 3*muu4**2 - ured))) + 
     &                   g2R*M4Lr*(g3R*g4L*muu3*(-1 + muu1**2 + muu3**2 - ured) + g3L*g4R*muu4*(muu2**2 + muu3**2 - ured))) + 
     &                g1L*(g2L*M4Lr*muu2*muu4*(g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) + g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                   g2R*(g3L*(g4L*mGlr*muu3*(muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 - 2*ured) + ured - muu3**2*ured - 
     &                            muu4**2*ured + ured**2) + g4R*muu4**2*
     &                          (muu1**2*(2*muu2**2 + muu3**2 - 2*ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))) + 
     &                      g3R*muu4*(g4L*muu3*(2*muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-2 + 2*muu3**2 - 3*ured) + ured - muu3**2*ured - 
     &                            muu4**2*ured + ured**2) + g4R*mGlr*
     &                          (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))))))*dconjg(AL) + 
     &          (BR*(g1L*muu1*muu4*(g2R*muu2*(g3R*muu4*(2*g4R*mGlr*muu3 + g4L*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)) + 
     &                      g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2) + g4R*muu3*(muu1**2 + 3*muu4**2 - ured))) + 
     &                   g2L*M4Lr*(g3L*g4R*muu3*(-1 + muu1**2 + muu3**2 - ured) + g3R*g4L*muu4*(muu2**2 + muu3**2 - ured))) + 
     &                g1R*(g2R*M4Lr*muu2*muu4*(g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) + g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                   g2L*(g3R*(g4R*mGlr*muu3*(muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 - 2*ured) + ured - muu3**2*ured - 
     &                            muu4**2*ured + ured**2) + g4L*muu4**2*
     &                          (muu1**2*(2*muu2**2 + muu3**2 - 2*ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))) + 
     &                      g3L*muu4*(g4R*muu3*(2*muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-2 + 2*muu3**2 - 3*ured) + ured - muu3**2*ured - 
     &                            muu4**2*ured + ured**2) + g4L*mGlr*
     &                          (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))))) + 
     &             BL*(g1L*muu1*muu4*(g2R*muu2*(g3L*muu4*(2*g4L*mGlr*muu3 + g4R*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)) + 
     &                      g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2) + g4L*muu3*(muu1**2 + 3*muu4**2 - ured))) + 
     &                   g2L*M4Lr*(g3R*g4L*muu3*(-1 + muu1**2 + muu3**2 - ured) + g3L*g4R*muu4*(muu2**2 + muu3**2 - ured))) + 
     &                g1R*(g2R*M4Lr*muu2*muu4*(g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) + g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                   g2L*(g3L*(g4L*mGlr*muu3*(muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 - 2*ured) + ured - muu3**2*ured - 
     &                            muu4**2*ured + ured**2) + g4R*muu4**2*
     &                          (muu1**2*(2*muu2**2 + muu3**2 - 2*ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))) + 
     &                      g3R*muu4*(g4L*muu3*(2*muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-2 + 2*muu3**2 - 3*ured) + ured - muu3**2*ured - 
     &                            muu4**2*ured + ured**2) + g4R*mGlr*
     &                          (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))))))*dconjg(AR)) + 
     &       xd2*((BR*(g1R*muu1*(g2R*M4Lr*(g3L*g4R*(-1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                      g3R*g4L*(-muu3**4 + muu4**2*(muu2**2 - ured) + muu3**2*(1 - muu1**2 + muu4**2 + ured))) + 
     &                   g2L*muu2*(g3L*muu4*(-(g4L*mGlr*(1 + muu3**2 - muu4**2)) + g4R*muu3*(2 + muu1**2 - muu2**2 - muu3**2 + muu4**2)) + 
     &                      g3R*(g4R*mGlr*muu3*(1 - muu3**2 + muu4**2) + 
     &                         g4L*(-(muu1**2*muu3**2) + 2*muu4**4 + muu4**2*(-2 + muu2**2 - 2*muu3**2 - ured) + muu3**2*ured)))) + 
     &                g1L*(g2L*M4Lr*muu2*(g3L*g4R*(1 + muu1**2 - muu2**2)*muu3*muu4 - 
     &                      g3R*g4L*(muu1**2*muu3**2 - muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured))) + 
     &                   g2R*(g3L*muu4*(g4R*muu3*(2*muu1**4 - muu1**2*(2 + 2*muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                         g4L*mGlr*(-muu2**4 + muu1**2*(muu3**2 - ured) + muu2**2*(1 + muu1**2 - muu4**2 + ured))) + 
     &                      g3R*(g4R*mGlr*muu3*(muu1**4 - muu1**2*(1 + muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                         g4L*(-(muu1**4*muu3**2) - muu2**2*muu4**2*(-1 + muu2**2 + muu4**2 - ured) + 
     &                            muu1**2*(-muu3**4 + 2*muu4**2*(muu2**2 - ured) + muu3**2*(1 + 2*muu4**2 + ured))))))) + 
     &             BL*(g1R*muu1*(g2R*M4Lr*(g3R*g4L*(-1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                      g3L*g4R*(-muu3**4 + muu4**2*(muu2**2 - ured) + muu3**2*(1 - muu1**2 + muu4**2 + ured))) + 
     &                   g2L*muu2*(g3R*muu4*(-(g4R*mGlr*(1 + muu3**2 - muu4**2)) + g4L*muu3*(2 + muu1**2 - muu2**2 - muu3**2 + muu4**2)) + 
     &                      g3L*(g4L*mGlr*muu3*(1 - muu3**2 + muu4**2) + 
     &                         g4R*(-(muu1**2*muu3**2) + 2*muu4**4 + muu4**2*(-2 + muu2**2 - 2*muu3**2 - ured) + muu3**2*ured)))) + 
     &                g1L*(g2L*M4Lr*muu2*(g3R*g4L*(1 + muu1**2 - muu2**2)*muu3*muu4 - 
     &                      g3L*g4R*(muu1**2*muu3**2 - muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured))) + 
     &                   g2R*(g3R*muu4*(g4L*muu3*(2*muu1**4 - muu1**2*(2 + 2*muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                         g4R*mGlr*(-muu2**4 + muu1**2*(muu3**2 - ured) + muu2**2*(1 + muu1**2 - muu4**2 + ured))) + 
     &                      g3L*(g4L*mGlr*muu3*(muu1**4 - muu1**2*(1 + muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                         g4R*(-(muu1**4*muu3**2) - muu2**2*muu4**2*(-1 + muu2**2 + muu4**2 - ured) + 
     &                            muu1**2*(-muu3**4 + 2*muu4**2*(muu2**2 - ured) + muu3**2*(1 + 2*muu4**2 + ured))))))))*dconjg(AL) + 
     &          (BR*(g1L*muu1*(g2L*M4Lr*(g3L*g4R*(-1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                      g3R*g4L*(-muu3**4 + muu4**2*(muu2**2 - ured) + muu3**2*(1 - muu1**2 + muu4**2 + ured))) + 
     &                   g2R*muu2*(g3L*muu4*(-(g4L*mGlr*(1 + muu3**2 - muu4**2)) + g4R*muu3*(2 + muu1**2 - muu2**2 - muu3**2 + muu4**2)) + 
     &                      g3R*(g4R*mGlr*muu3*(1 - muu3**2 + muu4**2) + 
     &                         g4L*(-(muu1**2*muu3**2) + 2*muu4**4 + muu4**2*(-2 + muu2**2 - 2*muu3**2 - ured) + muu3**2*ured)))) + 
     &                g1R*(g2R*M4Lr*muu2*(g3L*g4R*(1 + muu1**2 - muu2**2)*muu3*muu4 - 
     &                      g3R*g4L*(muu1**2*muu3**2 - muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured))) + 
     &                   g2L*(g3L*muu4*(g4R*muu3*(2*muu1**4 - muu1**2*(2 + 2*muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                         g4L*mGlr*(-muu2**4 + muu1**2*(muu3**2 - ured) + muu2**2*(1 + muu1**2 - muu4**2 + ured))) + 
     &                      g3R*(g4R*mGlr*muu3*(muu1**4 - muu1**2*(1 + muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                         g4L*(-(muu1**4*muu3**2) - muu2**2*muu4**2*(-1 + muu2**2 + muu4**2 - ured) + 
     &                            muu1**2*(-muu3**4 + 2*muu4**2*(muu2**2 - ured) + muu3**2*(1 + 2*muu4**2 + ured))))))) + 
     &             BL*(g1L*muu1*(g2L*M4Lr*(g3R*g4L*(-1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                      g3L*g4R*(-muu3**4 + muu4**2*(muu2**2 - ured) + muu3**2*(1 - muu1**2 + muu4**2 + ured))) + 
     &                   g2R*muu2*(g3R*muu4*(-(g4R*mGlr*(1 + muu3**2 - muu4**2)) + g4L*muu3*(2 + muu1**2 - muu2**2 - muu3**2 + muu4**2)) + 
     &                      g3L*(g4L*mGlr*muu3*(1 - muu3**2 + muu4**2) + 
     &                         g4R*(-(muu1**2*muu3**2) + 2*muu4**4 + muu4**2*(-2 + muu2**2 - 2*muu3**2 - ured) + muu3**2*ured)))) + 
     &                g1R*(g2R*M4Lr*muu2*(g3R*g4L*(1 + muu1**2 - muu2**2)*muu3*muu4 - 
     &                      g3L*g4R*(muu1**2*muu3**2 - muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured))) + 
     &                   g2L*(g3R*muu4*(g4L*muu3*(2*muu1**4 - muu1**2*(2 + 2*muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                         g4R*mGlr*(-muu2**4 + muu1**2*(muu3**2 - ured) + muu2**2*(1 + muu1**2 - muu4**2 + ured))) + 
     &                      g3L*(g4L*mGlr*muu3*(muu1**4 - muu1**2*(1 + muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                         g4R*(-(muu1**4*muu3**2) - muu2**2*muu4**2*(-1 + muu2**2 + muu4**2 - ured) + 
     &                            muu1**2*(-muu3**4 + 2*muu4**2*(muu2**2 - ured) + muu3**2*(1 + 2*muu4**2 + ured))))))))*dconjg(AR))))/(4.*DenSmV1*Pi**2)
	
		tempBoxZ = tempBoxZ + ChiChi2QQ_M2BoxGluino - tempBox2
		tempBox2 = ChiChi2QQ_M2BoxGluino


! Hboson contr.

	do i = imin,imax
	
	ChiChi2QQ_M2BoxGluino = ChiChi2QQ_M2BoxGluino -(CF*NC*sqrtS**6*((g1R*muu1*(2*g2R*M4Lr*muu2*muu4*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &             g2L*(g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2)*(-1 + muu1**2 + muu3**2 - ured) + 
     &                   2*g4R*muu3*muu4**2*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured)) + 
     &                g3R*muu4*(2*g4R*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4L*(-1 + muu3**2 + muu4**2)*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured))
     &                )) + g1L*(g2L*M4Lr*(-1 + muu1**2 + muu2**2)*muu4*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &             g2R*muu2*(g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2)*(muu1**2 + muu4**2 - ured) + 2*g4R*muu3*muu4**2*(3*muu1**2 + muu4**2 - ured)) + 
     &                g3R*muu4*(2*g4R*mGlr*muu3*(muu1**2 + muu4**2 - ured) + g4L*(-1 + muu3**2 + muu4**2)*(3*muu1**2 + muu4**2 - ured)))))*xd1*
     &        dconjg(CL(i))*dconjg(DL(i)) + (-(g1L*g2R*muu2*(-(g3R*g4L*muu4*
     &                  (-1 + muu3**2 - muu3**2*muu4**2 + muu4**4 + 2*muu1**2*(-1 + muu4**2) - muu2**2*(-1 + muu3**2 + muu4**2) + ured + muu3**2*ured - 
     &                    muu4**2*ured)) + g3L*g4R*muu3*(-muu4**4 + muu1**2*(-1 + muu3**2 - 3*muu4**2) + ured - muu3**2*ured + 
     &                  muu4**2*(-3 + 2*muu2**2 + muu3**2 + ured)))) - 
     &          g1R*g2L*muu1*(g3R*g4L*muu4*(-2 + muu3**2 + muu3**4 + 2*muu4**2 - muu3**2*muu4**2 - 2*muu1**2*(-1 + muu4**2) + 
     &                muu2**2*(-1 + muu3**2 + muu4**2) - ured - muu3**2*ured + muu4**2*ured) + 
     &             g3L*g4R*muu3*(1 + muu3**4 + 3*muu4**2 + 2*muu2**2*muu4**2 + muu1**2*(-1 + muu3**2 - 3*muu4**2) + ured + muu4**2*ured - 
     &                muu3**2*(2 + muu4**2 + ured))))*xd12*dconjg(CL(i))*dconjg(DL(i)) + 
     &       (g1L*g2R*muu2*(g3L*g4R*muu3*(muu1**4 + muu1**2*(6*muu4**2 - 2*ured) + (muu4**2 - ured)**2) + 
     &             g3R*g4L*muu4*((muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured) + muu1**2*(-3 + muu2**2 + 2*muu3**2 + 3*muu4**2 - ured))) + 
     &          g1R*g2L*muu1*(g3R*g4L*muu4*(2 - 2*muu3**2 - 2*muu4**2 + muu3**2*muu4**2 + muu2**2*(-2 + 2*muu3**2 + muu4**2 - ured) + 
     &                muu1**2*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured) + 2*ured - muu3**2*ured - muu4**2*ured + ured**2) + 
     &             g3L*g4R*muu3*(muu1**4 + muu1**2*(-1 + muu3**2 + 3*muu4**2 - 2*ured) + muu4**2*(-3 + 2*muu2**2 + muu3**2 - ured) + 
     &                ured*(1 - muu3**2 + ured))))*xd13*dconjg(CL(i))*dconjg(DL(i)) + 
     &       (g1L*(-(g2L*M4Lr*(-1 + muu1**2 + muu2**2)*(g3L*g4R*muu3*(-1 + muu3**2 - muu4**2) + g3R*g4L*muu4*(1 + muu3**2 - muu4**2))) + 
     &             g2R*muu2*(g3R*muu4*(2*g4R*mGlr*(1 + muu1**2 - muu2**2)*muu3 + 
     &                   g4L*(-(muu1**2*(3 + muu3**2 - 3*muu4**2)) - (-1 + muu2**2)*(-1 + muu3**2 + muu4**2))) + 
     &                g3L*(g4L*mGlr*(1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   2*g4R*muu3*(-((-1 + muu2**2)*muu4**2) + muu1**2*(1 - muu3**2 + 2*muu4**2))))) + 
     &          g1R*muu1*(2*g2R*M4Lr*muu2*(g3R*g4L*muu4*(-1 - muu3**2 + muu4**2) + g3L*g4R*muu3*(1 - muu3**2 + muu4**2)) + 
     &             g2L*(2*g3R*muu4*(g4R*mGlr*(-1 + muu1**2 - muu2**2)*muu3 - g4L*(-1 + muu2**2*muu3**2 + muu4**2 - muu1**2*(-1 + muu4**2))) + 
     &                g3L*(g4L*mGlr*(-1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   g4R*muu3*(-1 + muu3**2 - 3*muu4**2 - muu2**2*(-1 + muu3**2 + muu4**2) + muu1**2*(1 - muu3**2 + 3*muu4**2))))))*xd2*dconjg(CL(i))*
     &        dconjg(DL(i)) + (-(g1L*g2R*muu2*(g3R*g4L*muu4*((-1 + muu2**2)*(-1 + muu2**2 + muu4**2 - ured) + 
     &                  muu1**2*(3 - muu2**2 + 2*muu3**2 - 3*muu4**2 + ured)) + 
     &               g3L*g4R*muu3*(-muu1**4 + (-1 + muu2**2)*(muu4**2 - ured) + muu1**2*(-3 + muu2**2 + 2*muu3**2 - 3*muu4**2 + ured)))) + 
     &          g1R*g2L*muu1*(g3R*g4L*muu4*(2 - muu2**4 + muu3**2 - 2*muu4**2 + muu1**2*(-2 + muu2**2 - muu3**2 + 2*muu4**2 - ured) + ured + 
     &                muu2**2*(-1 - muu3**2 + ured)) + g3L*g4R*muu3*
     &              (-1 + muu1**4 + muu3**2 - 2*muu4**2 + ured + muu2**2*(1 - muu3**2 + ured) - muu1**2*(muu2**2 + muu3**2 - 2*muu4**2 + ured))))*xd23*
     &        dconjg(CL(i))*dconjg(DL(i)) + (g1R*muu1*(g2L*(-1 + muu1**2 + muu2**2)*
     &              (g3R*muu4*(2*g4R*mGlr*muu3 + g4L*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)) + 
     &                g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2) + g4R*muu3*(muu1**2 + 3*muu4**2 - ured))) + 
     &             2*g2R*M4Lr*muu2*(g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) + g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured))) + 
     &          g1L*(2*g2R*muu1**2*muu2*(g3R*muu4*(2*g4R*mGlr*muu3 + g4L*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)) + 
     &                g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2) + g4R*muu3*(muu1**2 + 3*muu4**2 - ured))) + 
     &             g2L*M4Lr*(-1 + muu1**2 + muu2**2)*(g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) + g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured))))*xd3*
     &        dconjg(CL(i))*dconjg(DL(i)) + (g1L*muu1*(2*g2L*M4Lr*muu2*muu4*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &             g2R*(g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2)*(-1 + muu1**2 + muu3**2 - ured) + 
     &                   2*g4R*muu3*muu4**2*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured)) + 
     &                g3R*muu4*(2*g4R*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4L*(-1 + muu3**2 + muu4**2)*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured))
     &                )) + g1R*(g2R*M4Lr*(-1 + muu1**2 + muu2**2)*muu4*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &             g2L*muu2*(g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2)*(muu1**2 + muu4**2 - ured) + 2*g4R*muu3*muu4**2*(3*muu1**2 + muu4**2 - ured)) + 
     &                g3R*muu4*(2*g4R*mGlr*muu3*(muu1**2 + muu4**2 - ured) + g4L*(-1 + muu3**2 + muu4**2)*(3*muu1**2 + muu4**2 - ured)))))*xd1*
     &        dconjg(CR(i))*dconjg(DL(i)) + (-(g1R*g2L*muu2*(-(g3R*g4L*muu4*
     &                  (-1 + muu3**2 - muu3**2*muu4**2 + muu4**4 + 2*muu1**2*(-1 + muu4**2) - muu2**2*(-1 + muu3**2 + muu4**2) + ured + muu3**2*ured - 
     &                    muu4**2*ured)) + g3L*g4R*muu3*(-muu4**4 + muu1**2*(-1 + muu3**2 - 3*muu4**2) + ured - muu3**2*ured + 
     &                  muu4**2*(-3 + 2*muu2**2 + muu3**2 + ured)))) - 
     &          g1L*g2R*muu1*(g3R*g4L*muu4*(-2 + muu3**2 + muu3**4 + 2*muu4**2 - muu3**2*muu4**2 - 2*muu1**2*(-1 + muu4**2) + 
     &                muu2**2*(-1 + muu3**2 + muu4**2) - ured - muu3**2*ured + muu4**2*ured) + 
     &             g3L*g4R*muu3*(1 + muu3**4 + 3*muu4**2 + 2*muu2**2*muu4**2 + muu1**2*(-1 + muu3**2 - 3*muu4**2) + ured + muu4**2*ured - 
     &                muu3**2*(2 + muu4**2 + ured))))*xd12*dconjg(CR(i))*dconjg(DL(i)) + 
     &       (g1R*g2L*muu2*(g3L*g4R*muu3*(muu1**4 + muu1**2*(6*muu4**2 - 2*ured) + (muu4**2 - ured)**2) + 
     &             g3R*g4L*muu4*((muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured) + muu1**2*(-3 + muu2**2 + 2*muu3**2 + 3*muu4**2 - ured))) + 
     &          g1L*g2R*muu1*(g3R*g4L*muu4*(2 - 2*muu3**2 - 2*muu4**2 + muu3**2*muu4**2 + muu2**2*(-2 + 2*muu3**2 + muu4**2 - ured) + 
     &                muu1**2*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured) + 2*ured - muu3**2*ured - muu4**2*ured + ured**2) + 
     &             g3L*g4R*muu3*(muu1**4 + muu1**2*(-1 + muu3**2 + 3*muu4**2 - 2*ured) + muu4**2*(-3 + 2*muu2**2 + muu3**2 - ured) + 
     &                ured*(1 - muu3**2 + ured))))*xd13*dconjg(CR(i))*dconjg(DL(i)) + 
     &       (g1R*(-(g2R*M4Lr*(-1 + muu1**2 + muu2**2)*(g3L*g4R*muu3*(-1 + muu3**2 - muu4**2) + g3R*g4L*muu4*(1 + muu3**2 - muu4**2))) + 
     &             g2L*muu2*(g3R*muu4*(2*g4R*mGlr*(1 + muu1**2 - muu2**2)*muu3 + 
     &                   g4L*(-(muu1**2*(3 + muu3**2 - 3*muu4**2)) - (-1 + muu2**2)*(-1 + muu3**2 + muu4**2))) + 
     &                g3L*(g4L*mGlr*(1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   2*g4R*muu3*(-((-1 + muu2**2)*muu4**2) + muu1**2*(1 - muu3**2 + 2*muu4**2))))) + 
     &          g1L*muu1*(2*g2L*M4Lr*muu2*(g3R*g4L*muu4*(-1 - muu3**2 + muu4**2) + g3L*g4R*muu3*(1 - muu3**2 + muu4**2)) + 
     &             g2R*(2*g3R*muu4*(g4R*mGlr*(-1 + muu1**2 - muu2**2)*muu3 - g4L*(-1 + muu2**2*muu3**2 + muu4**2 - muu1**2*(-1 + muu4**2))) + 
     &                g3L*(g4L*mGlr*(-1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   g4R*muu3*(-1 + muu3**2 - 3*muu4**2 - muu2**2*(-1 + muu3**2 + muu4**2) + muu1**2*(1 - muu3**2 + 3*muu4**2))))))*xd2*dconjg(CR(i))*
     &        dconjg(DL(i)) + (-(g1R*g2L*muu2*(g3R*g4L*muu4*((-1 + muu2**2)*(-1 + muu2**2 + muu4**2 - ured) + 
     &                  muu1**2*(3 - muu2**2 + 2*muu3**2 - 3*muu4**2 + ured)) + 
     &               g3L*g4R*muu3*(-muu1**4 + (-1 + muu2**2)*(muu4**2 - ured) + muu1**2*(-3 + muu2**2 + 2*muu3**2 - 3*muu4**2 + ured)))) + 
     &          g1L*g2R*muu1*(g3R*g4L*muu4*(2 - muu2**4 + muu3**2 - 2*muu4**2 + muu1**2*(-2 + muu2**2 - muu3**2 + 2*muu4**2 - ured) + ured + 
     &                muu2**2*(-1 - muu3**2 + ured)) + g3L*g4R*muu3*
     &              (-1 + muu1**4 + muu3**2 - 2*muu4**2 + ured + muu2**2*(1 - muu3**2 + ured) - muu1**2*(muu2**2 + muu3**2 - 2*muu4**2 + ured))))*xd23*
     &        dconjg(CR(i))*dconjg(DL(i)) + (g1L*muu1*(g2R*(-1 + muu1**2 + muu2**2)*
     &              (g3R*muu4*(2*g4R*mGlr*muu3 + g4L*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)) + 
     &                g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2) + g4R*muu3*(muu1**2 + 3*muu4**2 - ured))) + 
     &             2*g2L*M4Lr*muu2*(g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) + g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured))) + 
     &          g1R*(2*g2L*muu1**2*muu2*(g3R*muu4*(2*g4R*mGlr*muu3 + g4L*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)) + 
     &                g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2) + g4R*muu3*(muu1**2 + 3*muu4**2 - ured))) + 
     &             g2R*M4Lr*(-1 + muu1**2 + muu2**2)*(g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) + g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured))))*xd3*
     &        dconjg(CR(i))*dconjg(DL(i)) + (g1R*muu1*(2*g2R*M4Lr*muu2*muu4*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &             g2L*(g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2)*(-1 + muu1**2 + muu3**2 - ured) + 
     &                   2*g4L*muu3*muu4**2*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured)) + 
     &                g3L*muu4*(2*g4L*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4R*(-1 + muu3**2 + muu4**2)*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured))
     &                )) + g1L*(g2L*M4Lr*(-1 + muu1**2 + muu2**2)*muu4*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &             g2R*muu2*(g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2)*(muu1**2 + muu4**2 - ured) + 2*g4L*muu3*muu4**2*(3*muu1**2 + muu4**2 - ured)) + 
     &                g3L*muu4*(2*g4L*mGlr*muu3*(muu1**2 + muu4**2 - ured) + g4R*(-1 + muu3**2 + muu4**2)*(3*muu1**2 + muu4**2 - ured)))))*xd1*
     &        dconjg(CL(i))*dconjg(DR(i)) + (-(g1L*g2R*muu2*(-(g3L*g4R*muu4*
     &                  (-1 + muu3**2 - muu3**2*muu4**2 + muu4**4 + 2*muu1**2*(-1 + muu4**2) - muu2**2*(-1 + muu3**2 + muu4**2) + ured + muu3**2*ured - 
     &                    muu4**2*ured)) + g3R*g4L*muu3*(-muu4**4 + muu1**2*(-1 + muu3**2 - 3*muu4**2) + ured - muu3**2*ured + 
     &                  muu4**2*(-3 + 2*muu2**2 + muu3**2 + ured)))) - 
     &          g1R*g2L*muu1*(g3L*g4R*muu4*(-2 + muu3**2 + muu3**4 + 2*muu4**2 - muu3**2*muu4**2 - 2*muu1**2*(-1 + muu4**2) + 
     &                muu2**2*(-1 + muu3**2 + muu4**2) - ured - muu3**2*ured + muu4**2*ured) + 
     &             g3R*g4L*muu3*(1 + muu3**4 + 3*muu4**2 + 2*muu2**2*muu4**2 + muu1**2*(-1 + muu3**2 - 3*muu4**2) + ured + muu4**2*ured - 
     &                muu3**2*(2 + muu4**2 + ured))))*xd12*dconjg(CL(i))*dconjg(DR(i)) + 
     &       (g1L*g2R*muu2*(g3R*g4L*muu3*(muu1**4 + muu1**2*(6*muu4**2 - 2*ured) + (muu4**2 - ured)**2) + 
     &             g3L*g4R*muu4*((muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured) + muu1**2*(-3 + muu2**2 + 2*muu3**2 + 3*muu4**2 - ured))) + 
     &          g1R*g2L*muu1*(g3L*g4R*muu4*(2 - 2*muu3**2 - 2*muu4**2 + muu3**2*muu4**2 + muu2**2*(-2 + 2*muu3**2 + muu4**2 - ured) + 
     &                muu1**2*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured) + 2*ured - muu3**2*ured - muu4**2*ured + ured**2) + 
     &             g3R*g4L*muu3*(muu1**4 + muu1**2*(-1 + muu3**2 + 3*muu4**2 - 2*ured) + muu4**2*(-3 + 2*muu2**2 + muu3**2 - ured) + 
     &                ured*(1 - muu3**2 + ured))))*xd13*dconjg(CL(i))*dconjg(DR(i)) + 
     &       (g1L*(-(g2L*M4Lr*(-1 + muu1**2 + muu2**2)*(g3R*g4L*muu3*(-1 + muu3**2 - muu4**2) + g3L*g4R*muu4*(1 + muu3**2 - muu4**2))) + 
     &             g2R*muu2*(g3L*muu4*(2*g4L*mGlr*(1 + muu1**2 - muu2**2)*muu3 + 
     &                   g4R*(-(muu1**2*(3 + muu3**2 - 3*muu4**2)) - (-1 + muu2**2)*(-1 + muu3**2 + muu4**2))) + 
     &                g3R*(g4R*mGlr*(1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   2*g4L*muu3*(-((-1 + muu2**2)*muu4**2) + muu1**2*(1 - muu3**2 + 2*muu4**2))))) + 
     &          g1R*muu1*(2*g2R*M4Lr*muu2*(g3L*g4R*muu4*(-1 - muu3**2 + muu4**2) + g3R*g4L*muu3*(1 - muu3**2 + muu4**2)) + 
     &             g2L*(2*g3L*muu4*(g4L*mGlr*(-1 + muu1**2 - muu2**2)*muu3 - g4R*(-1 + muu2**2*muu3**2 + muu4**2 - muu1**2*(-1 + muu4**2))) + 
     &                g3R*(g4R*mGlr*(-1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   g4L*muu3*(-1 + muu3**2 - 3*muu4**2 - muu2**2*(-1 + muu3**2 + muu4**2) + muu1**2*(1 - muu3**2 + 3*muu4**2))))))*xd2*dconjg(CL(i))*
     &        dconjg(DR(i)) + (-(g1L*g2R*muu2*(g3L*g4R*muu4*((-1 + muu2**2)*(-1 + muu2**2 + muu4**2 - ured) + 
     &                  muu1**2*(3 - muu2**2 + 2*muu3**2 - 3*muu4**2 + ured)) + 
     &               g3R*g4L*muu3*(-muu1**4 + (-1 + muu2**2)*(muu4**2 - ured) + muu1**2*(-3 + muu2**2 + 2*muu3**2 - 3*muu4**2 + ured)))) + 
     &          g1R*g2L*muu1*(g3L*g4R*muu4*(2 - muu2**4 + muu3**2 - 2*muu4**2 + muu1**2*(-2 + muu2**2 - muu3**2 + 2*muu4**2 - ured) + ured + 
     &                muu2**2*(-1 - muu3**2 + ured)) + g3R*g4L*muu3*
     &              (-1 + muu1**4 + muu3**2 - 2*muu4**2 + ured + muu2**2*(1 - muu3**2 + ured) - muu1**2*(muu2**2 + muu3**2 - 2*muu4**2 + ured))))*xd23*
     &        dconjg(CL(i))*dconjg(DR(i)) + (g1R*muu1*(g2L*(-1 + muu1**2 + muu2**2)*
     &              (g3L*muu4*(2*g4L*mGlr*muu3 + g4R*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)) + 
     &                g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2) + g4L*muu3*(muu1**2 + 3*muu4**2 - ured))) + 
     &             2*g2R*M4Lr*muu2*(g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) + g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured))) + 
     &          g1L*(2*g2R*muu1**2*muu2*(g3L*muu4*(2*g4L*mGlr*muu3 + g4R*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)) + 
     &                g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2) + g4L*muu3*(muu1**2 + 3*muu4**2 - ured))) + 
     &             g2L*M4Lr*(-1 + muu1**2 + muu2**2)*(g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) + g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured))))*xd3*
     &        dconjg(CL(i))*dconjg(DR(i)) + (g1L*muu1*(2*g2L*M4Lr*muu2*muu4*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &             g2R*(g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2)*(-1 + muu1**2 + muu3**2 - ured) + 
     &                   2*g4L*muu3*muu4**2*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured)) + 
     &                g3L*muu4*(2*g4L*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4R*(-1 + muu3**2 + muu4**2)*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured))
     &                )) + g1R*(g2R*M4Lr*(-1 + muu1**2 + muu2**2)*muu4*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &             g2L*muu2*(g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2)*(muu1**2 + muu4**2 - ured) + 2*g4L*muu3*muu4**2*(3*muu1**2 + muu4**2 - ured)) + 
     &                g3L*muu4*(2*g4L*mGlr*muu3*(muu1**2 + muu4**2 - ured) + g4R*(-1 + muu3**2 + muu4**2)*(3*muu1**2 + muu4**2 - ured)))))*xd1*
     &        dconjg(CR(i))*dconjg(DR(i)) + (-(g1R*g2L*muu2*(-(g3L*g4R*muu4*
     &                  (-1 + muu3**2 - muu3**2*muu4**2 + muu4**4 + 2*muu1**2*(-1 + muu4**2) - muu2**2*(-1 + muu3**2 + muu4**2) + ured + muu3**2*ured - 
     &                    muu4**2*ured)) + g3R*g4L*muu3*(-muu4**4 + muu1**2*(-1 + muu3**2 - 3*muu4**2) + ured - muu3**2*ured + 
     &                  muu4**2*(-3 + 2*muu2**2 + muu3**2 + ured)))) - 
     &          g1L*g2R*muu1*(g3L*g4R*muu4*(-2 + muu3**2 + muu3**4 + 2*muu4**2 - muu3**2*muu4**2 - 2*muu1**2*(-1 + muu4**2) + 
     &                muu2**2*(-1 + muu3**2 + muu4**2) - ured - muu3**2*ured + muu4**2*ured) + 
     &             g3R*g4L*muu3*(1 + muu3**4 + 3*muu4**2 + 2*muu2**2*muu4**2 + muu1**2*(-1 + muu3**2 - 3*muu4**2) + ured + muu4**2*ured - 
     &                muu3**2*(2 + muu4**2 + ured))))*xd12*dconjg(CR(i))*dconjg(DR(i)) + 
     &       (g1R*g2L*muu2*(g3R*g4L*muu3*(muu1**4 + muu1**2*(6*muu4**2 - 2*ured) + (muu4**2 - ured)**2) + 
     &             g3L*g4R*muu4*((muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured) + muu1**2*(-3 + muu2**2 + 2*muu3**2 + 3*muu4**2 - ured))) + 
     &          g1L*g2R*muu1*(g3L*g4R*muu4*(2 - 2*muu3**2 - 2*muu4**2 + muu3**2*muu4**2 + muu2**2*(-2 + 2*muu3**2 + muu4**2 - ured) + 
     &                muu1**2*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured) + 2*ured - muu3**2*ured - muu4**2*ured + ured**2) + 
     &             g3R*g4L*muu3*(muu1**4 + muu1**2*(-1 + muu3**2 + 3*muu4**2 - 2*ured) + muu4**2*(-3 + 2*muu2**2 + muu3**2 - ured) + 
     &                ured*(1 - muu3**2 + ured))))*xd13*dconjg(CR(i))*dconjg(DR(i)) + 
     &       (g1R*(-(g2R*M4Lr*(-1 + muu1**2 + muu2**2)*(g3R*g4L*muu3*(-1 + muu3**2 - muu4**2) + g3L*g4R*muu4*(1 + muu3**2 - muu4**2))) + 
     &             g2L*muu2*(g3L*muu4*(2*g4L*mGlr*(1 + muu1**2 - muu2**2)*muu3 + 
     &                   g4R*(-(muu1**2*(3 + muu3**2 - 3*muu4**2)) - (-1 + muu2**2)*(-1 + muu3**2 + muu4**2))) + 
     &                g3R*(g4R*mGlr*(1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   2*g4L*muu3*(-((-1 + muu2**2)*muu4**2) + muu1**2*(1 - muu3**2 + 2*muu4**2))))) + 
     &          g1L*muu1*(2*g2L*M4Lr*muu2*(g3L*g4R*muu4*(-1 - muu3**2 + muu4**2) + g3R*g4L*muu3*(1 - muu3**2 + muu4**2)) + 
     &             g2R*(2*g3L*muu4*(g4L*mGlr*(-1 + muu1**2 - muu2**2)*muu3 - g4R*(-1 + muu2**2*muu3**2 + muu4**2 - muu1**2*(-1 + muu4**2))) + 
     &                g3R*(g4R*mGlr*(-1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   g4L*muu3*(-1 + muu3**2 - 3*muu4**2 - muu2**2*(-1 + muu3**2 + muu4**2) + muu1**2*(1 - muu3**2 + 3*muu4**2))))))*xd2*dconjg(CR(i))*
     &        dconjg(DR(i)) + (-(g1R*g2L*muu2*(g3L*g4R*muu4*((-1 + muu2**2)*(-1 + muu2**2 + muu4**2 - ured) + 
     &                  muu1**2*(3 - muu2**2 + 2*muu3**2 - 3*muu4**2 + ured)) + 
     &               g3R*g4L*muu3*(-muu1**4 + (-1 + muu2**2)*(muu4**2 - ured) + muu1**2*(-3 + muu2**2 + 2*muu3**2 - 3*muu4**2 + ured)))) + 
     &          g1L*g2R*muu1*(g3L*g4R*muu4*(2 - muu2**4 + muu3**2 - 2*muu4**2 + muu1**2*(-2 + muu2**2 - muu3**2 + 2*muu4**2 - ured) + ured + 
     &                muu2**2*(-1 - muu3**2 + ured)) + g3R*g4L*muu3*
     &              (-1 + muu1**4 + muu3**2 - 2*muu4**2 + ured + muu2**2*(1 - muu3**2 + ured) - muu1**2*(muu2**2 + muu3**2 - 2*muu4**2 + ured))))*xd23*
     &        dconjg(CR(i))*dconjg(DR(i)) + (g1L*muu1*(g2R*(-1 + muu1**2 + muu2**2)*
     &              (g3L*muu4*(2*g4L*mGlr*muu3 + g4R*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)) + 
     &                g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2) + g4L*muu3*(muu1**2 + 3*muu4**2 - ured))) + 
     &             2*g2L*M4Lr*muu2*(g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) + g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured))) + 
     &          g1R*(2*g2L*muu1**2*muu2*(g3L*muu4*(2*g4L*mGlr*muu3 + g4R*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)) + 
     &                g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2) + g4L*muu3*(muu1**2 + 3*muu4**2 - ured))) + 
     &             g2R*M4Lr*(-1 + muu1**2 + muu2**2)*(g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) + g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured))))*xd3*
     &        dconjg(CR(i))*dconjg(DR(i)) - xd22*((g1R*g2L*muu1*(-1 + muu1**2 - muu2**2) + g1L*g2R*muu2*(1 + muu1**2 - muu2**2))*dconjg(CL(i)) + 
     &          (g1L*g2R*muu1*(-1 + muu1**2 - muu2**2) + g1R*g2L*muu2*(1 + muu1**2 - muu2**2))*dconjg(CR(i)))*
     &        ((g3L*g4R*muu3*(-1 + muu3**2 - muu4**2) + g3R*g4L*muu4*(1 + muu3**2 - muu4**2))*dconjg(DL(i)) + 
     &          (g3R*g4L*muu3*(-1 + muu3**2 - muu4**2) + g3L*g4R*muu4*(1 + muu3**2 - muu4**2))*dconjg(DR(i))) + 
     &       muu4*xd11*((g1R*g2L*muu1*(-1 + muu1**2 + muu3**2 - ured) + g1L*g2R*muu2*(muu1**2 + muu4**2 - ured))*dconjg(CL(i)) + 
     &          (g1L*g2R*muu1*(-1 + muu1**2 + muu3**2 - ured) + g1R*g2L*muu2*(muu1**2 + muu4**2 - ured))*dconjg(CR(i)))*
     &        ((2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2))*dconjg(DL(i)) + 
     &          (2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2))*dconjg(DR(i))) + 
     &       xd0*((g1R*muu1*(2*g2R*M4Lr*muu2 + g2L*(-1 + muu1**2 + muu2**2)) + g1L*(2*g2R*muu1**2*muu2 + g2L*M4Lr*(-1 + muu1**2 + muu2**2)))*
     &           dconjg(CL(i)) + (g1L*muu1*(2*g2L*M4Lr*muu2 + g2R*(-1 + muu1**2 + muu2**2)) + 
     &             g1R*(2*g2L*muu1**2*muu2 + g2R*M4Lr*(-1 + muu1**2 + muu2**2)))*dconjg(CR(i)))*
     &        ((g3R*muu4*(2*g4R*mGlr*muu3 + g4L*(-1 + muu3**2 + muu4**2)) + g3L*(2*g4R*muu3*muu4**2 + g4L*mGlr*(-1 + muu3**2 + muu4**2)))*dconjg(DL(i)) + 
     &          (g3L*muu4*(2*g4L*mGlr*muu3 + g4R*(-1 + muu3**2 + muu4**2)) + g3R*(2*g4L*muu3*muu4**2 + g4R*mGlr*(-1 + muu3**2 + muu4**2)))*dconjg(DR(i)))
     &        + muu1*xd33*((2*g1L*g2R*muu1*muu2 + g1R*g2L*(-1 + muu1**2 + muu2**2))*dconjg(CL(i)) + 
     &          (2*g1R*g2L*muu1*muu2 + g1L*g2R*(-1 + muu1**2 + muu2**2))*dconjg(CR(i)))*
     &        ((g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) + g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured))*dconjg(DL(i)) + 
     &          (g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) + g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured))*dconjg(DR(i))) + 
     &       (2*xd00*(dconjg(CR(i))*((g1L*g2R*muu1*(g3L*g4R*muu3*(-1 + muu1**2 + muu3**2 - ured) + g3R*g4L*muu4*(muu2**2 + muu3**2 - ured)) + 
     &                  g1R*g2L*muu2*(g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) + g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured)))*dconjg(DL(i)) + 
     &               (g1L*g2R*muu1*(g3R*g4L*muu3*(-1 + muu1**2 + muu3**2 - ured) + g3L*g4R*muu4*(muu2**2 + muu3**2 - ured)) + 
     &                  g1R*g2L*muu2*(g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) + g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured)))*dconjg(DR(i))) + 
     &            dconjg(CL(i))*((g1R*g2L*muu1*(g3L*g4R*muu3*(-1 + muu1**2 + muu3**2 - ured) + g3R*g4L*muu4*(muu2**2 + muu3**2 - ured)) + 
     &                  g1L*g2R*muu2*(g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) + g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured)))*dconjg(DL(i)) + 
     &               (g1R*g2L*muu1*(g3R*g4L*muu3*(-1 + muu1**2 + muu3**2 - ured) + g3L*g4R*muu4*(muu2**2 + muu3**2 - ured)) + 
     &                  g1L*g2R*muu2*(g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) + g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured)))*dconjg(DR(i)))))/sqrtS**2))/
     &  (8.*Pi**2*DenSmH1(i))
	
	enddo

	tempBoxH = tempBoxH + ChiChi2QQ_M2BoxGluino - tempBox2
	tempBox2 = ChiChi2QQ_M2BoxGluino

	do i = 1, 2

! U-channel contr.
	
	ChiChi2QQ_M2BoxGluino = ChiChi2QQ_M2BoxGluino + (CF*NC*sqrtS**4*((xd00*(g1L*g2R*(-(g3L*g4R*((muu1**2 + muu4**2 - ured)*dconjg(GL(indi,i,tt,gen)) + 2*muu1*muu4*dconjg(GR(indi,i,tt,gen)))*
     &                 (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g3R*g4L*(muu4*dconjg(GL(indi,i,tt,gen))*(muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                    (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 muu1*dconjg(GR(indi,i,tt,gen))*(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                    muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))))) + 
     &           g1R*g2L*(-(g3R*g4L*(2*muu1*muu4*dconjg(GL(indi,i,tt,gen)) + (muu1**2 + muu4**2 - ured)*dconjg(GR(indi,i,tt,gen)))*
     &                 ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen)))) + 
     &              g3L*g4R*(muu1*dconjg(GL(indi,i,tt,gen))*(muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                    muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                 muu4*dconjg(GR(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + 
     &                    muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))))))/sqrtS**2 - 
     &      (muu1*xd33*(g1R*g2L*(g3R*g4L*muu1*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                  ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*((muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))*
     &                     dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g3L*g4R*(dconjg(GL(indi,i,tt,gen))*(muu3*(muu1**2*(muu2**2 - muu3**2 + muu4**2) + (-1 + muu2**2)*(muu4**2 - ured))*
     &                     dconjg(HL(indj,i,tt,gen)) + muu2*(muu1**2*(muu2**2 - muu3**2 - ured) + (muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))*
     &                     dconjg(HR(indj,i,tt,gen))) + muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                  ((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen))))) + 
     &           g1L*g2R*(g3R*g4L*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                  (muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*(muu2*(muu1**2*(muu2**2 - muu3**2 - ured) + (muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))*
     &                     dconjg(HL(indj,i,tt,gen)) + muu3*(muu1**2*(muu2**2 - muu3**2 + muu4**2) + (-1 + muu2**2)*(muu4**2 - ured))*
     &                     dconjg(HR(indj,i,tt,gen)))) + g3L*g4R*muu1*
     &               (muu1*muu4*dconjg(GR(indi,i,tt,gen))*(2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                    (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))*dconjg(HR(indj,i,tt,gen)))))))/2. - 
     &      (muu4*xd11*(g1R*g2L*(g3R*g4L*muu4*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                  ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*((muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))*
     &                     dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g3L*g4R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                  (muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*(muu3*(muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 - 2*ured) + ured - 
     &                       muu3**2*ured - muu4**2*ured + ured**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                    muu2*(-(muu2**2*muu4**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + ured - muu3**2*ured)*dconjg(HR(indj,i,tt,gen)))))
     &             + g1L*g2R*(g3R*g4L*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                  (muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) + muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(muu2*(-(muu2**2*muu4**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + ured - muu3**2*ured)*
     &                     dconjg(HL(indj,i,tt,gen)) + muu3*(muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 - 2*ured) + ured - 
     &                       muu3**2*ured - muu4**2*ured + ured**2)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g3L*g4R*muu4*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                  (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                    (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))*dconjg(HR(indj,i,tt,gen)))))))/2. + 
     &      (xd13*(-(g1R*g2L*(g3R*g4L*(muu1**2 + muu4**2 - ured)*
     &                 (muu1*muu4*dconjg(GL(indi,i,tt,gen))*((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                      2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                   dconjg(GR(indi,i,tt,gen))*((muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))*
     &                       dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &                2*g3L*g4R*muu1*muu4*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                    (muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                   muu4*dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + 
     &                      muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))))) - 
     &           g1L*g2R*(2*g3R*g4L*muu1*muu4*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                  (muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 muu1*dconjg(GL(indi,i,tt,gen))*(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                    muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g3L*g4R*(muu1**2 + muu4**2 - ured)*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                  (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                    (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))*dconjg(HR(indj,i,tt,gen)))))))/2. + 
     &      (xd23*(-(g1R*g2L*(g3R*g4L*(muu1*(1 + muu1**2 - muu2**2)*muu4*dconjg(GL(indi,i,tt,gen))*
     &                    ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                   dconjg(GR(indi,i,tt,gen))*((muu1**4*(muu2**2 - muu3**2 - ured) - muu2**2*(muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured) + 
     &                         muu1**2*(-muu2**4 + muu3**2 - muu3**4 + 2*muu3**2*muu4**2 + muu2**2*(muu3**2 + muu4**2) + ured - 2*muu4**2*ured + ured**2))*
     &                       dconjg(HL(indj,i,tt,gen)) + muu2*(1 + muu1**2 - muu2**2)*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))) - 
     &                2*g3L*g4R*muu1*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                    ((1 - muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + muu2*(1 + muu3**2 - muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                   dconjg(GL(indi,i,tt,gen))*(muu3*(muu4**2 + muu1**2*(muu3**2 - muu4**2) - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                      muu2*(muu1**2*muu3**2 - muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured))*dconjg(HR(indj,i,tt,gen)))))) - 
     &           g1L*g2R*(-2*g3R*g4L*muu1*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                  (muu2*(1 + muu3**2 - muu4**2)*dconjg(HL(indj,i,tt,gen)) + (1 - muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*(muu2*(muu1**2*muu3**2 - muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured))*
     &                     dconjg(HL(indj,i,tt,gen)) + muu3*(muu4**2 + muu1**2*(muu3**2 - muu4**2) - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g3L*g4R*(muu1*(1 + muu1**2 - muu2**2)*muu4*dconjg(GR(indi,i,tt,gen))*
     &                  (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(muu2*(1 + muu1**2 - muu2**2)*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                    (muu1**4*(muu2**2 - muu3**2 - ured) - muu2**2*(muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured) + 
     &                       muu1**2*(-muu2**4 + muu3**2 - muu3**4 + 2*muu3**2*muu4**2 + muu2**2*(muu3**2 + muu4**2) + ured - 2*muu4**2*ured + ured**2))*
     &                     dconjg(HR(indj,i,tt,gen)))))))/2. + 
     &      (xd22*(-(g1R*g2L*(g3R*g4L*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                    ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                   dconjg(GR(indi,i,tt,gen))*((-(muu1**4*muu3**2) + muu2**2*(-muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured)) + 
     &                         muu1**2*(-muu3**4 + muu4**2*(muu2**2 - ured) + muu3**2*(1 + muu2**2 + muu4**2 + ured)))*dconjg(HL(indj,i,tt,gen)) + 
     &                      muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))) - 
     &                g3L*g4R*(muu1*dconjg(GL(indi,i,tt,gen))*(muu3*(muu4**2 + muu1**2*(muu3**2 - muu4**2) + muu2**2*(1 - muu3**2 + muu4**2) - ured)*
     &                       dconjg(HL(indj,i,tt,gen)) + muu2*(-muu3**4 + muu4**2 - muu4**4 + muu3**2*(1 + 2*muu4**2))*dconjg(HR(indj,i,tt,gen))) + 
     &                   muu4*dconjg(GR(indi,i,tt,gen))*((-muu1**4 + muu2**2 - muu2**4 + muu1**2*(1 + 2*muu2**2))*muu3*dconjg(HL(indj,i,tt,gen)) + 
     &                      muu2*(muu3**2 - muu2**2*muu3**2 + muu2**2*muu4**2 + muu1**2*(1 + muu3**2 - muu4**2) - ured)*dconjg(HR(indj,i,tt,gen)))))) - 
     &           g1L*g2R*(-(g3R*g4L*(muu4*dconjg(GL(indi,i,tt,gen))*
     &                    (muu2*(muu3**2 - muu2**2*muu3**2 + muu2**2*muu4**2 + muu1**2*(1 + muu3**2 - muu4**2) - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                      (-muu1**4 + muu2**2 - muu2**4 + muu1**2*(1 + 2*muu2**2))*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                   muu1*dconjg(GR(indi,i,tt,gen))*(muu2*(-muu3**4 + muu4**2 - muu4**4 + muu3**2*(1 + 2*muu4**2))*dconjg(HL(indj,i,tt,gen)) + 
     &                      muu3*(muu4**2 + muu1**2*(muu3**2 - muu4**2) + muu2**2*(1 - muu3**2 + muu4**2) - ured)*dconjg(HR(indj,i,tt,gen))))) + 
     &              g3L*g4R*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                  (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                    (-(muu1**4*muu3**2) + muu2**2*(-muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured)) + 
     &                       muu1**2*(-muu3**4 + muu4**2*(muu2**2 - ured) + muu3**2*(1 + muu2**2 + muu4**2 + ured)))*dconjg(HR(indj,i,tt,gen)))))))/2. + 
     &      (xd12*(g1R*g2L*(g3R*g4L*(muu1*muu4*(-1 + muu3**2 - muu4**2)*dconjg(GL(indi,i,tt,gen))*
     &                  ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*((muu1**4*muu3**2 + muu2**4*muu4**2 + muu2**2*muu4**2*(-1 - muu3**2 + muu4**2) + 
     &                       (muu3**2 - muu4**2)*(muu3**2*(muu4**2 - ured) + ured*(1 - muu4**2 + ured)) + 
     &                       muu1**2*(muu3**4 + 2*muu4**2*(-muu2**2 + ured) - muu3**2*(1 + muu4**2 + 2*ured)))*dconjg(HL(indj,i,tt,gen)) + 
     &                    muu2*muu3*(-1 + muu3**2 - muu4**2)*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))) - 
     &              2*g3L*g4R*muu4*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                  ((-1 + muu1**2 - muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) - muu2*(1 + muu3**2 - muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*(muu3*(muu1**4 - muu1**2*(1 + muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured))*
     &                     dconjg(HL(indj,i,tt,gen)) + muu2*(-(muu2**2*muu4**2) + muu1**2*(-1 + muu4**2) + ured)*dconjg(HR(indj,i,tt,gen))))) + 
     &           g1L*g2R*(-2*g3R*g4L*muu4*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                  (-(muu2*(1 + muu3**2 - muu4**2)*dconjg(HL(indj,i,tt,gen))) + (-1 + muu1**2 - muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(muu2*(-(muu2**2*muu4**2) + muu1**2*(-1 + muu4**2) + ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                    muu3*(muu1**4 - muu1**2*(1 + muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured))*dconjg(HR(indj,i,tt,gen)))) + 
     &              g3L*g4R*(muu1*muu4*(-1 + muu3**2 - muu4**2)*dconjg(GR(indi,i,tt,gen))*
     &                  (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(-1 + muu3**2 - muu4**2)*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                    (muu1**4*muu3**2 + muu2**4*muu4**2 + muu2**2*muu4**2*(-1 - muu3**2 + muu4**2) + 
     &                       (muu3**2 - muu4**2)*(muu3**2*(muu4**2 - ured) + ured*(1 - muu4**2 + ured)) + 
     &                       muu1**2*(muu3**4 + 2*muu4**2*(-muu2**2 + ured) - muu3**2*(1 + muu4**2 + 2*ured)))*dconjg(HR(indj,i,tt,gen)))))))/2. + 
     &      (xd0*(-(g1R*(g2L*muu1*(g3R*(g4L*muu4*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                          ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                         dconjg(GL(indi,i,tt,gen))*((muu1**2*muu3**2 + muu2**2*muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &                      g4R*mGlr*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                          (muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen)))
     &                          + muu4*dconjg(GL(indi,i,tt,gen))*
     &                          ((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen))))
     &                      ) + g3L*(g4L*mGlr*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                          ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                         dconjg(GL(indi,i,tt,gen))*((muu1**2*muu3**2 + muu2**2*muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &                      g4R*muu4*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                          (muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen)))
     &                          + muu4*dconjg(GL(indi,i,tt,gen))*
     &                          ((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen))))
     &                      )) + g2R*M4Lr*(g3R*(g4L*muu4*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                          (muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen)))
     &                          + muu1*dconjg(GL(indi,i,tt,gen))*
     &                          (muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) + muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))))
     &                        + g4R*mGlr*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                          (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                         dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                            (muu1**2*muu3**2 + muu2**2*muu4**2 - ured)*dconjg(HR(indj,i,tt,gen))))) + 
     &                   g3L*(g4L*mGlr*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                          (muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen)))
     &                          + muu1*dconjg(GL(indi,i,tt,gen))*
     &                          (muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) + muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))))
     &                        + g4R*muu4*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                          (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                         dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                            (muu1**2*muu3**2 + muu2**2*muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))))))) - 
     &           g1L*(g2L*M4Lr*(g3R*(g4L*muu4*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                        ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GL(indi,i,tt,gen))*((muu1**2*muu3**2 + muu2**2*muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4R*mGlr*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                        (muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu4*dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen))))) + 
     &                 g3L*(g4L*mGlr*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                        ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GL(indi,i,tt,gen))*((muu1**2*muu3**2 + muu2**2*muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4R*muu4*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                        (muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu4*dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))))) + 
     &              g2R*muu1*(g3R*(g4L*muu4*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu1*dconjg(GL(indi,i,tt,gen))*(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4R*mGlr*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                        (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                          (muu1**2*muu3**2 + muu2**2*muu4**2 - ured)*dconjg(HR(indj,i,tt,gen))))) + 
     &                 g3L*(g4L*mGlr*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu1*dconjg(GL(indi,i,tt,gen))*(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4R*muu4*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                        (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                          (muu1**2*muu3**2 + muu2**2*muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))))))))/2. + 
     &      (xd3*(-(g1R*(g2R*M4Lr*(g3R*g4L*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                       (muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                      dconjg(GR(indi,i,tt,gen))*(muu2*(muu1**2*(muu2**2 - muu3**2 - ured) + (muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))*
     &                          dconjg(HL(indj,i,tt,gen)) + muu3*(muu1**2*(muu2**2 - muu3**2 + muu4**2) + (-1 + muu2**2)*(muu4**2 - ured))*
     &                          dconjg(HR(indj,i,tt,gen)))) + g3L*g4R*muu1*
     &                    (muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                       (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                      dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                         (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))*dconjg(HR(indj,i,tt,gen))))) + 
     &                g2L*muu1*(g3L*(g4L*mGlr*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                          ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                         dconjg(GL(indi,i,tt,gen))*((muu1**2*muu3**2 + muu2**2*muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &                      g4R*(dconjg(GL(indi,i,tt,gen))*(muu3*(muu1**2*(muu2**2 - muu3**2 + 2*muu4**2) + (-1 + muu2**2)*(2*muu4**2 - ured))*
     &                             dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*(muu1**2*(muu2**2 - muu3**2 - ured) + (-1 + muu2**2 + muu4**2 - ured)*(2*muu4**2 - ured))*dconjg(HR(indj,i,tt,gen)))
     &                           + muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                          (muu3*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)*dconjg(HR(indj,i,tt,gen))))) + 
     &                   g3R*(g4R*mGlr*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                          (muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen)))
     &                          + muu4*dconjg(GL(indi,i,tt,gen))*
     &                          ((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen))))
     &                        + g4L*(muu4*dconjg(GL(indi,i,tt,gen))*
     &                          ((muu2**2*muu4**2 + muu1**2*(muu2**2 + 2*muu3**2 - ured) - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*muu3*(3*muu1**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                         muu1*dconjg(GR(indi,i,tt,gen))*((2*muu3**2*muu4**2 + muu1**2*(muu2**2 - ured) + muu2**2*(muu4**2 - ured) + ured - 
     &                               muu3**2*ured - 2*muu4**2*ured + ured**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*muu3*(muu1**2 + 3*muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))))))) - 
     &           g1L*(g2L*M4Lr*(g3R*g4L*muu1*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                     ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                    dconjg(GR(indi,i,tt,gen))*((muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))*
     &                        dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &                 g3L*g4R*(dconjg(GL(indi,i,tt,gen))*(muu3*(muu1**2*(muu2**2 - muu3**2 + muu4**2) + (-1 + muu2**2)*(muu4**2 - ured))*
     &                        dconjg(HL(indj,i,tt,gen)) + muu2*(muu1**2*(muu2**2 - muu3**2 - ured) + (muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))*
     &                        dconjg(HR(indj,i,tt,gen))) + muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                     ((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen))))) + 
     &              g2R*muu1*(g3R*(g4L*(dconjg(GR(indi,i,tt,gen))*
     &                        (muu2*(muu1**2*(muu2**2 - muu3**2 - ured) + (-1 + muu2**2 + muu4**2 - ured)*(2*muu4**2 - ured))*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu3*(muu1**2*(muu2**2 - muu3**2 + 2*muu4**2) + (-1 + muu2**2)*(2*muu4**2 - ured))*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                        (muu2*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu3*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4R*mGlr*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                        (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                          (muu1**2*muu3**2 + muu2**2*muu4**2 - ured)*dconjg(HR(indj,i,tt,gen))))) + 
     &                 g3L*(g4L*mGlr*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu1*dconjg(GL(indi,i,tt,gen))*(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4R*(muu4*dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(3*muu1**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                          (muu2**2*muu4**2 + muu1**2*(muu2**2 + 2*muu3**2 - ured) - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu1*dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(muu1**2 + 3*muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                          (2*muu3**2*muu4**2 + muu1**2*(muu2**2 - ured) + muu2**2*(muu4**2 - ured) + ured - muu3**2*ured - 2*muu4**2*ured + ured**2)*
     &                           dconjg(HR(indj,i,tt,gen)))))))))/2. + 
     &      (xd1*(-(g1R*(g2R*M4Lr*muu4*(g3R*g4L*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                       (muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                      muu1*dconjg(GL(indi,i,tt,gen))*(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                         muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &                   g3L*g4R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                       (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                      dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                         (muu1**2*muu3**2 + muu2**2*muu4**2 - ured)*dconjg(HR(indj,i,tt,gen))))) + 
     &                g2L*(g3R*(g4L*muu4*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                          ((muu1**2*(2*muu2**2 + muu3**2 - 2*ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))*
     &                             dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(3*muu1**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                         muu1*dconjg(GL(indi,i,tt,gen))*((muu1**2*muu3**2 + 2*muu2**2*muu4**2 + muu3**2*muu4**2 - ured - muu4**2*ured)*
     &                             dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(muu1**2 + 3*muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &                      g4R*mGlr*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                          (muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen)))
     &                          + dconjg(GR(indi,i,tt,gen))*(muu3*
     &                             (muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 - 2*ured) + ured - muu3**2*ured - muu4**2*ured + 
     &                               ured**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*(-(muu2**2*muu4**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + ured - muu3**2*ured)*
     &                             dconjg(HR(indj,i,tt,gen))))) + 
     &                   g3L*muu4*(g4L*mGlr*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                          ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                         dconjg(GR(indi,i,tt,gen))*((muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))*
     &                             dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &                      g4R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                          (muu3*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                         dconjg(GR(indi,i,tt,gen))*(muu3*
     &                             (2*muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-2 + 2*muu3**2 - 3*ured) + ured - muu3**2*ured - 
     &                               muu4**2*ured + ured**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*(-(muu2**2*muu4**2) + muu3**2*muu4**2 + 2*muu1**2*(-1 + muu3**2 + muu4**2) + ured - muu3**2*ured)*
     &                             dconjg(HR(indj,i,tt,gen)))))))) - 
     &           g1L*(g2L*M4Lr*muu4*(g3R*g4L*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                     ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                    dconjg(GL(indi,i,tt,gen))*((muu1**2*muu3**2 + muu2**2*muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                       muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen)))) + 
     &                 g3L*g4R*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                     (muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                    muu4*dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + 
     &                       muu2*(-1 + muu2**2 + muu4**2 - ured)*dconjg(HR(indj,i,tt,gen))))) + 
     &              g2R*(g3L*(g4L*mGlr*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) + muu3*(-1 + muu1**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GL(indi,i,tt,gen))*(muu2*(-(muu2**2*muu4**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + ured - 
     &                             muu3**2*ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu3*(muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 - 2*ured) + ured - muu3**2*ured - muu4**2*ured + 
     &                             ured**2)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4R*muu4*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                        (muu2*muu3*(muu1**2 + 3*muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                          (muu1**2*muu3**2 + 2*muu2**2*muu4**2 + muu3**2*muu4**2 - ured - muu4**2*ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu4*dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(3*muu1**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                          (muu1**2*(2*muu2**2 + muu3**2 - 2*ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))*
     &                           dconjg(HR(indj,i,tt,gen))))) + 
     &                 g3R*muu4*(g4L*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (muu2*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu3*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GL(indi,i,tt,gen))*(muu2*(-(muu2**2*muu4**2) + muu3**2*muu4**2 + 2*muu1**2*(-1 + muu3**2 + muu4**2) + ured - 
     &                             muu3**2*ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu3*(2*muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-2 + 2*muu3**2 - 3*ured) + ured - muu3**2*ured - 
     &                             muu4**2*ured + ured**2)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4R*mGlr*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                          (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))*dconjg(HR(indj,i,tt,gen))))))))
     &         )/2. + (xd2*(-(g1R*(g2R*M4Lr*(-(g3R*g4L*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                         (muu2*(1 + muu3**2 - muu4**2)*dconjg(HL(indj,i,tt,gen)) + (1 - muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                        dconjg(GR(indi,i,tt,gen))*(muu2*(muu1**2*muu3**2 - muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured))*
     &                            dconjg(HL(indj,i,tt,gen)) + muu3*(muu4**2 + muu1**2*(muu3**2 - muu4**2) - ured)*dconjg(HR(indj,i,tt,gen))))) + 
     &                   g3L*g4R*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                       (muu2*(1 + muu1**2 - muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + 
     &                         (muu3**2 - muu2**2*muu3**2 + muu2**2*muu4**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                      muu1*dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(1 - muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                         (-muu3**4 + muu4**2*(muu2**2 - ured) + muu3**2*(1 - muu1**2 + muu4**2 + ured))*dconjg(HR(indj,i,tt,gen))))) + 
     &                g2L*(g3R*(g4L*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                          (((1 + muu1**2)*muu3**2 + muu2**2*(1 - 2*muu3**2 + muu4**2) - 2*ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*muu3*(2 + muu1**2 - muu2**2 - muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                         dconjg(GR(indi,i,tt,gen))*((-(muu1**4*muu3**2) - muu2**2*muu4**2*(-1 + muu2**2 + muu4**2 - ured) + 
     &                               muu1**2*(-muu3**4 + 2*muu4**2*(muu2**2 - ured) + muu3**2*(1 + 2*muu4**2 + ured)))*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*muu3*(-((-1 + muu2**2)*muu4**2) + muu1**2*(1 - muu3**2 + 2*muu4**2))*dconjg(HR(indj,i,tt,gen)))) + 
     &                      g4R*mGlr*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                          ((-1 + muu1**2 - muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) - muu2*(1 + muu3**2 - muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                         dconjg(GR(indi,i,tt,gen))*(muu3*(muu1**4 - muu1**2*(1 + muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured))*
     &                             dconjg(HL(indj,i,tt,gen)) + muu2*(-(muu2**2*muu4**2) + muu1**2*(-1 + muu4**2) + ured)*dconjg(HR(indj,i,tt,gen))))) + 
     &                   g3L*(g4L*mGlr*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                          ((-muu2**4 + muu1**2*(muu3**2 - ured) + muu2**2*(1 + muu1**2 - muu4**2 + ured))*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*(1 + muu1**2 - muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                         muu1*dconjg(GL(indi,i,tt,gen))*((muu2**2 + muu1**2*muu3**2 - muu2**2*muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*muu3*(1 - muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen)))) - 
     &                      g4R*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                          (muu3*(-2*muu1**4 + muu2**2*(muu4**2 - ured) + muu1**2*(2 + 2*muu2**2 - muu3**2 + ured))*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*(muu2**2*muu4**2 + muu1**2*(2 + muu3**2 - 2*muu4**2) - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                         muu1*dconjg(GL(indi,i,tt,gen))*(muu3*((2 + muu2**2)*muu4**2 + muu1**2*(muu3**2 - 2*muu4**2) - ured)*
     &                             dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*(muu1**2*muu3**2 - 2*muu4**4 - muu3**2*ured + muu4**2*(2 - muu2**2 + 2*muu3**2 + ured))*dconjg(HR(indj,i,tt,gen)))))
     &                   ))) - g1L*(g2L*M4Lr*(g3R*g4L*(muu4*dconjg(GL(indi,i,tt,gen))*
     &                     ((muu3**2 - muu2**2*muu3**2 + muu2**2*muu4**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                       muu2*(1 + muu1**2 - muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                    muu1*dconjg(GR(indi,i,tt,gen))*((-muu3**4 + muu4**2*(muu2**2 - ured) + muu3**2*(1 - muu1**2 + muu4**2 + ured))*
     &                        dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(1 - muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen)))) - 
     &                 g3L*g4R*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                     ((1 - muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + muu2*(1 + muu3**2 - muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                    dconjg(GL(indi,i,tt,gen))*(muu3*(muu4**2 + muu1**2*(muu3**2 - muu4**2) - ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                       muu2*(muu1**2*muu3**2 - muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured))*dconjg(HR(indj,i,tt,gen))))) + 
     &              g2R*(g3R*(-(g4L*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                          (muu2*(muu1**2*muu3**2 - 2*muu4**4 - muu3**2*ured + muu4**2*(2 - muu2**2 + 2*muu3**2 + ured))*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu3*((2 + muu2**2)*muu4**2 + muu1**2*(muu3**2 - 2*muu4**2) - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                         muu4*dconjg(GL(indi,i,tt,gen))*(muu2*(muu2**2*muu4**2 + muu1**2*(2 + muu3**2 - 2*muu4**2) - ured)*
     &                             dconjg(HL(indj,i,tt,gen)) + 
     &                            muu3*(-2*muu1**4 + muu2**2*(muu4**2 - ured) + muu1**2*(2 + 2*muu2**2 - muu3**2 + ured))*dconjg(HR(indj,i,tt,gen))))) + 
     &                    g4R*mGlr*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                        (muu2*muu3*(1 - muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                          (muu2**2 + muu1**2*muu3**2 - muu2**2*muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu4*dconjg(GL(indi,i,tt,gen))*(muu2*(1 + muu1**2 - muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + 
     &                          (-muu2**4 + muu1**2*(muu3**2 - ured) + muu2**2*(1 + muu1**2 - muu4**2 + ured))*dconjg(HR(indj,i,tt,gen))))) + 
     &                 g3L*(g4L*mGlr*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (-(muu2*(1 + muu3**2 - muu4**2)*dconjg(HL(indj,i,tt,gen))) + (-1 + muu1**2 - muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GL(indi,i,tt,gen))*(muu2*(-(muu2**2*muu4**2) + muu1**2*(-1 + muu4**2) + ured)*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu3*(muu1**4 - muu1**2*(1 + muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured))*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4R*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (muu2*muu3*(2 + muu1**2 - muu2**2 - muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                          ((1 + muu1**2)*muu3**2 + muu2**2*(1 - 2*muu3**2 + muu4**2) - 2*ured)*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(-((-1 + muu2**2)*muu4**2) + muu1**2*(1 - muu3**2 + 2*muu4**2))*
     &                           dconjg(HL(indj,i,tt,gen)) + (-(muu1**4*muu3**2) - muu2**2*muu4**2*(-1 + muu2**2 + muu4**2 - ured) + 
     &                             muu1**2*(-muu3**4 + 2*muu4**2*(muu2**2 - ured) + muu3**2*(1 + 2*muu4**2 + ured)))*dconjg(HR(indj,i,tt,gen)))))))))/2.))/
     &  (4.*Pi**2*(-ured + MSf(i,itypeU,gen)**2/sqrtS**2))

	tempBoxu = tempBoxu + ChiChi2QQ_M2BoxGluino - tempBox2
	tempBox2 = ChiChi2QQ_M2BoxGluino

! T-channel contr.

	ChiChi2QQ_M2BoxGluino = ChiChi2QQ_M2BoxGluino + (CF*NC*sqrtS**4*((xd00*(dconjg(ER(indi,i,tt,gen))*((g1L*g2R*muu1*muu2*(4*g3R*g4L*muu3*muu4 - g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &                 g1R*g2L*(-(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4) + g3L*g4R*(-1 + muu1**2 + muu3**2 - ured)*(-1 + muu2**2 + muu4**2 - ured)))*
     &               dconjg(FL(indj,i,tt,gen)) + (g1L*g2R*muu1*
     &                  (-2*g3R*g4L*muu3*(-1 + muu1**2 + muu3**2 - ured) + g3L*g4R*muu4*(muu2**2 + muu3**2 - ured)) + 
     &                 g1R*g2L*muu2*(g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) - 2*g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured)))*dconjg(FR(indj,i,tt,gen)))
     &            + dconjg(ELc(indi,i,tt,gen))*((g1R*g2L*muu1*
     &                  (-2*g3L*g4R*muu3*(-1 + muu1**2 + muu3**2 - ured) + g3R*g4L*muu4*(muu2**2 + muu3**2 - ured)) + 
     &                 g1L*g2R*muu2*(g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) - 2*g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured)))*dconjg(FL(indj,i,tt,gen))
     &               + (g1R*g2L*muu1*muu2*(4*g3L*g4R*muu3*muu4 - g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &                 g1L*g2R*(-(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4) + g3R*g4L*(-1 + muu1**2 + muu3**2 - ured)*(-1 + muu2**2 + muu4**2 - ured)))*
     &               dconjg(FR(indj,i,tt,gen)))))/sqrtS**2 + (muu1*xd33*
     &         (dconjg(ER(indi,i,tt,gen))*((g1L*g2R*muu2*(2*g3R*g4L*muu1**2*muu3*muu4 + 
     &                    g3L*g4R*(muu1**2*(muu2**2 - muu3**2 - ured) + (muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                 g1R*g2L*muu1*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                    g3L*g4R*(muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))))*dconjg(FL(indj,i,tt,gen))
     &               - (g1L*g2R*(g3L*g4R*muu4*(muu1**2*(-1 - muu3**2 + muu4**2) + (-1 + muu2**2)*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                    g3R*g4L*muu1**2*muu3*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                 g1R*g2L*muu1*muu2*(g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) + g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured)))*
     &               dconjg(FR(indj,i,tt,gen))) + dconjg(ELc(indi,i,tt,gen))*
     &            (-((g1R*g2L*(g3R*g4L*muu4*(muu1**2*(-1 - muu3**2 + muu4**2) + (-1 + muu2**2)*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                      g3L*g4R*muu1**2*muu3*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                   g1L*g2R*muu1*muu2*(g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) + g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured)))*
     &                 dconjg(FL(indj,i,tt,gen))) + (g1R*g2L*muu2*
     &                  (2*g3L*g4R*muu1**2*muu3*muu4 + g3R*g4L*(muu1**2*(muu2**2 - muu3**2 - ured) + (muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                 g1L*g2R*muu1*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                    g3R*g4L*(muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))))*dconjg(FR(indj,i,tt,gen))))
     &         )/2. + (muu4*xd11*(dconjg(ER(indi,i,tt,gen))*((g1L*g2R*muu1*muu2*muu4*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &                 g1R*g2L*(g3R*g4L*muu3*(muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 - 2*ured) + ured - muu3**2*ured - 
     &                       muu4**2*ured + ured**2) + g3L*g4R*muu4*
     &                     (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))))*dconjg(FL(indj,i,tt,gen)) - 
     &              (g1R*g2L*muu2*muu4*(g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) + g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                 g1L*g2R*muu1*(g3R*g4L*muu3*muu4*(-1 + muu1**2 + muu3**2 - ured) + 
     &                    g3L*g4R*(1 + muu3**4 - muu4**2 - muu2**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + ured - muu3**2*(2 + ured))))*
     &               dconjg(FR(indj,i,tt,gen))) + dconjg(ELc(indi,i,tt,gen))*
     &            (-((g1L*g2R*muu2*muu4*(g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) + g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                   g1R*g2L*muu1*(g3L*g4R*muu3*muu4*(-1 + muu1**2 + muu3**2 - ured) + 
     &                      g3R*g4L*(1 + muu3**4 - muu4**2 - muu2**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + ured - muu3**2*(2 + ured))))*
     &                 dconjg(FL(indj,i,tt,gen))) + (g1R*g2L*muu1*muu2*muu4*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &                 g1L*g2R*(g3L*g4R*muu3*(muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 - 2*ured) + ured - muu3**2*ured - 
     &                       muu4**2*ured + ured**2) + g3R*g4L*muu4*
     &                     (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))))*dconjg(FR(indj,i,tt,gen)))))/2. + 
     &      (xd23*(dconjg(ER(indi,i,tt,gen))*((-2*g1L*g2R*muu1*muu2*
     &                  (g3R*g4L*(-1 - muu1**2 + muu2**2)*muu3*muu4 + 
     &                    g3L*g4R*(muu1**2*muu3**2 - muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured))) + 
     &                 g1R*g2L*(2*g3R*g4L*muu1**2*(-1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                    g3L*g4R*(muu1**4*(muu2**2 - muu3**2 - ured) - muu2**2*(muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured) + 
     &                       muu1**2*(-muu2**4 + muu3**2 - muu3**4 + 2*muu3**2*muu4**2 + muu2**2*(muu3**2 + muu4**2) + ured - 2*muu4**2*ured + ured**2))))*
     &               dconjg(FL(indj,i,tt,gen)) + (g1R*g2L*muu2*
     &                  (2*g3R*g4L*muu1**2*muu3*(-1 + muu3**2 - muu4**2) - g3L*g4R*(1 + muu1**2 - muu2**2)*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                 g1L*g2R*muu1*(-(g3R*g4L*(1 + muu1**2 - muu2**2)*muu3*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                    2*g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 + muu1**2*(1 + muu3**2 - muu4**2) - ured)))*dconjg(FR(indj,i,tt,gen))) + 
     &           dconjg(ELc(indi,i,tt,gen))*((g1L*g2R*muu2*(2*g3L*g4R*muu1**2*muu3*(-1 + muu3**2 - muu4**2) - 
     &                    g3R*g4L*(1 + muu1**2 - muu2**2)*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                 g1R*g2L*muu1*(-(g3L*g4R*(1 + muu1**2 - muu2**2)*muu3*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                    2*g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 + muu1**2*(1 + muu3**2 - muu4**2) - ured)))*dconjg(FL(indj,i,tt,gen)) + 
     &              (2*g1R*g2L*muu1*muu2*(g3L*g4R*(1 + muu1**2 - muu2**2)*muu3*muu4 - 
     &                    g3R*g4L*(muu1**2*muu3**2 - muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured))) + 
     &                 g1L*g2R*(2*g3L*g4R*muu1**2*(-1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                    g3R*g4L*(muu1**4*(muu2**2 - muu3**2 - ured) - muu2**2*(muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured) + 
     &                       muu1**2*(-muu2**4 + muu3**2 - muu3**4 + 2*muu3**2*muu4**2 + muu2**2*(muu3**2 + muu4**2) + ured - 2*muu4**2*ured + ured**2))))*
     &               dconjg(FR(indj,i,tt,gen)))))/2. + (xd22*(dconjg(ER(indi,i,tt,gen))*
     &            ((g1L*g2R*muu1*muu2*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(muu3**4 + muu4**2*(-1 + muu4**2) - muu3**2*(1 + 2*muu4**2))) + 
     &                 g1R*g2L*(g3R*g4L*(muu1**4 + muu2**2*(-1 + muu2**2) - muu1**2*(1 + 2*muu2**2))*muu3*muu4 + 
     &                    g3L*g4R*(-(muu1**4*muu3**2) + muu2**2*(-muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured)) + 
     &                       muu1**2*(-muu3**4 + muu4**2*(muu2**2 - ured) + muu3**2*(1 + muu2**2 + muu4**2 + ured)))))*dconjg(FL(indj,i,tt,gen)) + 
     &              (g1L*g2R*muu1*(-(g3R*g4L*muu3*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                    g3L*g4R*muu4*(-1 + muu4**2 + muu1**2*(1 + muu3**2 - muu4**2) + muu2**2*(-muu3**2 + muu4**2) - ured)) + 
     &                 g1R*g2L*muu2*(-(g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                    g3R*g4L*muu3*(-1 + (1 + muu1**2)*muu3**2 - muu1**2*muu4**2 + muu2**2*(1 - muu3**2 + muu4**2) - ured)))*dconjg(FR(indj,i,tt,gen)))
     &             + dconjg(ELc(indi,i,tt,gen))*((g1R*g2L*muu1*
     &                  (-(g3L*g4R*muu3*(-1 + muu1**2 + muu3**2 - ured)) + 
     &                    g3R*g4L*muu4*(-1 + muu4**2 + muu1**2*(1 + muu3**2 - muu4**2) + muu2**2*(-muu3**2 + muu4**2) - ured)) + 
     &                 g1L*g2R*muu2*(-(g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                    g3L*g4R*muu3*(-1 + (1 + muu1**2)*muu3**2 - muu1**2*muu4**2 + muu2**2*(1 - muu3**2 + muu4**2) - ured)))*dconjg(FL(indj,i,tt,gen))
     &               + (g1R*g2L*muu1*muu2*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(muu3**4 + muu4**2*(-1 + muu4**2) - muu3**2*(1 + 2*muu4**2))) + 
     &                 g1L*g2R*(g3L*g4R*(muu1**4 + muu2**2*(-1 + muu2**2) - muu1**2*(1 + 2*muu2**2))*muu3*muu4 + 
     &                    g3R*g4L*(-(muu1**4*muu3**2) + muu2**2*(-muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured)) + 
     &                       muu1**2*(-muu3**4 + muu4**2*(muu2**2 - ured) + muu3**2*(1 + muu2**2 + muu4**2 + ured)))))*dconjg(FR(indj,i,tt,gen)))))/2. + 
     &      (xd0*(-(dconjg(ELc(indi,i,tt,gen))*((g1R*muu1*(g2R*M4Lr*muu2*
     &                       (g3R*muu4*(2*g4R*mGlr*muu3 + g4L*(-1 + muu3**2 + muu4**2)) + g3L*(2*g4R*muu3*muu4**2 + g4L*mGlr*(-1 + muu3**2 + muu4**2))) + 
     &                      g2L*(g3R*muu4*(g4R*mGlr*(-1 + muu1**2 + muu2**2)*muu3 + 
     &                            g4L*(1 - muu3**2 + muu2**2*(-1 + muu3**2) - muu4**2 + muu1**2*(-1 + muu4**2) + ured)) + 
     &                         g3L*(g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4**2 + 
     &                            g4L*mGlr*(1 - muu3**2 + muu2**2*(-1 + muu3**2) - muu4**2 + muu1**2*(-1 + muu4**2) + ured)))) + 
     &                   g1L*(g2R*muu1**2*muu2*(g3R*muu4*(2*g4R*mGlr*muu3 + g4L*(-1 + muu3**2 + muu4**2)) + 
     &                         g3L*(2*g4R*muu3*muu4**2 + g4L*mGlr*(-1 + muu3**2 + muu4**2))) + 
     &                      g2L*M4Lr*(g3R*muu4*(g4R*mGlr*(-1 + muu1**2 + muu2**2)*muu3 + 
     &                            g4L*(1 - muu3**2 + muu2**2*(-1 + muu3**2) - muu4**2 + muu1**2*(-1 + muu4**2) + ured)) + 
     &                         g3L*(g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4**2 + 
     &                            g4L*mGlr*(1 - muu3**2 + muu2**2*(-1 + muu3**2) - muu4**2 + muu1**2*(-1 + muu4**2) + ured)))))*dconjg(FL(indj,i,tt,gen))
     &                 - (g1R*muu1*(g2R*M4Lr*(g3L*muu4*(g4R*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4L*mGlr*(muu2**2 + muu3**2 - ured)) + 
     &                         g3R*(g4R*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4L*muu4**2*(muu2**2 + muu3**2 - ured))) + 
     &                      g2L*muu2*(g3L*muu4*(g4R*muu3*(muu1**2 + muu4**2 - ured) + g4L*mGlr*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                         g3R*(g4R*mGlr*muu3*(muu1**2 + muu4**2 - ured) + g4L*muu4**2*(-1 + muu2**2 + muu4**2 - ured)))) + 
     &                   g1L*(g2R*muu1**2*(g3L*muu4*(g4R*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4L*mGlr*(muu2**2 + muu3**2 - ured)) + 
     &                         g3R*(g4R*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4L*muu4**2*(muu2**2 + muu3**2 - ured))) + 
     &                      g2L*M4Lr*muu2*(g3L*muu4*(g4R*muu3*(muu1**2 + muu4**2 - ured) + g4L*mGlr*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                         g3R*(g4R*mGlr*muu3*(muu1**2 + muu4**2 - ured) + g4L*muu4**2*(-1 + muu2**2 + muu4**2 - ured)))))*dconjg(FR(indj,i,tt,gen))))
     &            - dconjg(ER(indi,i,tt,gen))*(-((g1L*muu1*(g2L*M4Lr*
     &                       (g3R*muu4*(g4L*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4R*mGlr*(muu2**2 + muu3**2 - ured)) + 
     &                         g3L*(g4L*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4R*muu4**2*(muu2**2 + muu3**2 - ured))) + 
     &                      g2R*muu2*(g3R*muu4*(g4L*muu3*(muu1**2 + muu4**2 - ured) + g4R*mGlr*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                         g3L*(g4L*mGlr*muu3*(muu1**2 + muu4**2 - ured) + g4R*muu4**2*(-1 + muu2**2 + muu4**2 - ured)))) + 
     &                   g1R*(g2L*muu1**2*(g3R*muu4*(g4L*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4R*mGlr*(muu2**2 + muu3**2 - ured)) + 
     &                         g3L*(g4L*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + g4R*muu4**2*(muu2**2 + muu3**2 - ured))) + 
     &                      g2R*M4Lr*muu2*(g3R*muu4*(g4L*muu3*(muu1**2 + muu4**2 - ured) + g4R*mGlr*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                         g3L*(g4L*mGlr*muu3*(muu1**2 + muu4**2 - ured) + g4R*muu4**2*(-1 + muu2**2 + muu4**2 - ured)))))*dconjg(FL(indj,i,tt,gen)))
     &               + (g1L*muu1*(g2L*M4Lr*muu2*(g3L*muu4*(2*g4L*mGlr*muu3 + g4R*(-1 + muu3**2 + muu4**2)) + 
     &                       g3R*(2*g4L*muu3*muu4**2 + g4R*mGlr*(-1 + muu3**2 + muu4**2))) + 
     &                    g2R*(g3L*muu4*(g4L*mGlr*(-1 + muu1**2 + muu2**2)*muu3 + 
     &                          g4R*(1 - muu3**2 + muu2**2*(-1 + muu3**2) - muu4**2 + muu1**2*(-1 + muu4**2) + ured)) + 
     &                       g3R*(g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4**2 + 
     &                          g4R*mGlr*(1 - muu3**2 + muu2**2*(-1 + muu3**2) - muu4**2 + muu1**2*(-1 + muu4**2) + ured)))) + 
     &                 g1R*(g2L*muu1**2*muu2*(g3L*muu4*(2*g4L*mGlr*muu3 + g4R*(-1 + muu3**2 + muu4**2)) + 
     &                       g3R*(2*g4L*muu3*muu4**2 + g4R*mGlr*(-1 + muu3**2 + muu4**2))) + 
     &                    g2R*M4Lr*(g3L*muu4*(g4L*mGlr*(-1 + muu1**2 + muu2**2)*muu3 + 
     &                          g4R*(1 - muu3**2 + muu2**2*(-1 + muu3**2) - muu4**2 + muu1**2*(-1 + muu4**2) + ured)) + 
     &                       g3R*(g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4**2 + 
     &                          g4R*mGlr*(1 - muu3**2 + muu2**2*(-1 + muu3**2) - muu4**2 + muu1**2*(-1 + muu4**2) + ured)))))*dconjg(FR(indj,i,tt,gen)))))/
     &       2. + (xd13*(dconjg(ER(indi,i,tt,gen))*((2*g1L*g2R*muu1*muu2*muu4*
     &                  (g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) + g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                 g1R*g2L*(2*g3R*g4L*muu1**2*muu3*muu4*(-1 + muu1**2 + muu3**2 - ured) + 
     &                    g3L*g4R*(muu1**4*(muu2**2 - muu3**2 - ured) - 
     &                       (-1 + muu2**2 + muu4**2 - ured)*(muu2**2*muu4**2 + (-1 + muu4**2 - ured)*ured + muu3**2*(-muu4**2 + ured)) + 
     &                       muu1**2*(-muu3**4 + muu3**2*(1 + 2*muu4**2) + muu2**2*(-1 + muu3**2 + 2*muu4**2 - 2*ured) + ured*(2 - 3*muu4**2 + 2*ured)))))*
     &               dconjg(FL(indj,i,tt,gen)) - (g1R*g2L*muu2*muu4*
     &                  (4*g3R*g4L*muu1**2*muu3*muu4 + g3L*g4R*(muu1**2 + muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                 g1L*g2R*muu1*(g3R*g4L*muu3*(-1 + muu1**2 + muu3**2 - ured)*(muu1**2 + muu4**2 - ured) + 
     &                    2*g3L*g4R*muu4*(1 - muu3**2 + muu2**2*(-1 + muu3**2) - muu4**2 + muu1**2*(-1 + muu4**2) + ured)))*dconjg(FR(indj,i,tt,gen))) + 
     &           dconjg(ELc(indi,i,tt,gen))*(-((g1L*g2R*muu2*muu4*
     &                    (4*g3L*g4R*muu1**2*muu3*muu4 + g3R*g4L*(muu1**2 + muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                   g1R*g2L*muu1*(g3L*g4R*muu3*(-1 + muu1**2 + muu3**2 - ured)*(muu1**2 + muu4**2 - ured) + 
     &                      2*g3R*g4L*muu4*(1 - muu3**2 + muu2**2*(-1 + muu3**2) - muu4**2 + muu1**2*(-1 + muu4**2) + ured)))*dconjg(FL(indj,i,tt,gen))) + 
     &              (2*g1R*g2L*muu1*muu2*muu4*(g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) + g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                 g1L*g2R*(2*g3L*g4R*muu1**2*muu3*muu4*(-1 + muu1**2 + muu3**2 - ured) + 
     &                    g3R*g4L*(muu1**4*(muu2**2 - muu3**2 - ured) - 
     &                       (-1 + muu2**2 + muu4**2 - ured)*(muu2**2*muu4**2 + (-1 + muu4**2 - ured)*ured + muu3**2*(-muu4**2 + ured)) + 
     &                       muu1**2*(-muu3**4 + muu3**2*(1 + 2*muu4**2) + muu2**2*(-1 + muu3**2 + 2*muu4**2 - 2*ured) + ured*(2 - 3*muu4**2 + 2*ured)))))*
     &               dconjg(FR(indj,i,tt,gen)))))/2. + (xd12*(dconjg(ER(indi,i,tt,gen))*
     &            ((2*g1L*g2R*muu1*muu2*muu4*(g3L*g4R*muu4*(-1 - muu3**2 + muu4**2) + g3R*g4L*muu3*(1 - muu3**2 + muu4**2)) + 
     &                 g1R*g2L*(2*g3R*g4L*muu3*muu4*(muu1**4 - muu1**2*(1 + muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                    g3L*g4R*(-(muu1**4*muu3**2) - muu2**4*muu4**2 + muu2**2*muu4**2*(1 + muu3**2 - muu4**2) - 
     &                       (muu3**2 - muu4**2)*(muu3**2*(muu4**2 - ured) + ured*(1 - muu4**2 + ured)) + 
     &                       muu1**2*(-muu3**4 + 2*muu4**2*(muu2**2 - ured) + muu3**2*(1 + muu4**2 + 2*ured)))))*dconjg(FL(indj,i,tt,gen)) + 
     &              (g1R*g2L*muu2*muu4*(-2*g3R*g4L*(1 + muu1**2 - muu2**2)*muu3*muu4 + g3L*g4R*(-1 + muu3**2 - muu4**2)*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                 g1L*g2R*muu1*(g3R*g4L*muu3*(-1 + muu3**2 - muu4**2)*(-1 + muu1**2 + muu3**2 - ured) + 
     &                    2*g3L*g4R*muu4*(-1 + muu3**2 + muu4**2 + muu2**2*muu4**2 - muu1**2*(-1 + muu4**2) - ured)))*dconjg(FR(indj,i,tt,gen))) + 
     &           dconjg(ELc(indi,i,tt,gen))*((g1L*g2R*muu2*muu4*
     &                  (-2*g3L*g4R*(1 + muu1**2 - muu2**2)*muu3*muu4 + g3R*g4L*(-1 + muu3**2 - muu4**2)*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                 g1R*g2L*muu1*(g3L*g4R*muu3*(-1 + muu3**2 - muu4**2)*(-1 + muu1**2 + muu3**2 - ured) + 
     &                    2*g3R*g4L*muu4*(-1 + muu3**2 + muu4**2 + muu2**2*muu4**2 - muu1**2*(-1 + muu4**2) - ured)))*dconjg(FL(indj,i,tt,gen)) + 
     &              (2*g1R*g2L*muu1*muu2*muu4*(g3R*g4L*muu4*(-1 - muu3**2 + muu4**2) + g3L*g4R*muu3*(1 - muu3**2 + muu4**2)) + 
     &                 g1L*g2R*(2*g3L*g4R*muu3*muu4*(muu1**4 - muu1**2*(1 + muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                    g3R*g4L*(-(muu1**4*muu3**2) - muu2**4*muu4**2 + muu2**2*muu4**2*(1 + muu3**2 - muu4**2) - 
     &                       (muu3**2 - muu4**2)*(muu3**2*(muu4**2 - ured) + ured*(1 - muu4**2 + ured)) + 
     &                       muu1**2*(-muu3**4 + 2*muu4**2*(muu2**2 - ured) + muu3**2*(1 + muu4**2 + 2*ured)))))*dconjg(FR(indj,i,tt,gen)))))/2. + 
     &      (xd3*(dconjg(ER(indi,i,tt,gen))*((g1L*muu1*(g2R*muu2*
     &                     (g3L*(g4R*(muu1**2*(muu2**2 - muu3**2 - ured) + (-1 + muu2**2 + muu4**2 - ured)*(2*muu4**2 - ured)) + 
     &                          g4L*mGlr*muu3*(muu1**2 + muu4**2 - ured)) + 
     &                       g3R*muu4*(g4L*muu3*(3*muu1**2 + muu4**2 - ured) + g4R*mGlr*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                    g2L*M4Lr*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                       g3L*g4R*(muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))) + 
     &                 g1R*(g2R*M4Lr*muu2*(2*g3R*g4L*muu1**2*muu3*muu4 + 
     &                       g3L*g4R*(muu1**2*(muu2**2 - muu3**2 - ured) + (muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                    g2L*muu1**2*(g3R*muu4*(g4R*mGlr*(muu2**2 + muu3**2 - ured) + g4L*muu3*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured)) + 
     &                       g3L*(g4L*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + 
     &                          g4R*(2*muu3**2*muu4**2 + muu1**2*(muu2**2 - ured) + muu2**2*(muu4**2 - ured) + ured - muu3**2*ured - 2*muu4**2*ured + 
     &                             ured**2)))))*dconjg(FL(indj,i,tt,gen)) - 
     &              (g1R*(g2L*muu1**2*muu2*(g3L*muu4*(2*g4L*mGlr*muu3 + g4R*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)) + 
     &                       g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2) + g4L*muu3*(muu1**2 + 3*muu4**2 - ured))) + 
     &                    g2R*M4Lr*(g3L*g4R*muu4*(muu1**2*(-1 - muu3**2 + muu4**2) + (-1 + muu2**2)*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                       g3R*g4L*muu1**2*muu3*(-1 + muu1**2 + muu3**2 - ured))) + 
     &                 g1L*muu1*(g2L*M4Lr*muu2*(g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) + g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                    g2R*(g3R*(g4L*muu3*(muu1**4 + (-1 + muu2**2)*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2 - ured)) + 
     &                          g4R*mGlr*(1 - muu3**2 + muu2**2*(-1 + muu3**2) - muu4**2 + muu1**2*(-1 + muu4**2) + ured)) + 
     &                       g3L*muu4*(g4L*mGlr*(-1 + muu1**2 + muu2**2)*muu3 + 
     &                          g4R*(2 + muu2**4 - muu3**2 - 2*muu4**2 - muu1**2*(2 + muu3**2 - 2*muu4**2) + muu2**2*(-3 + muu3**2 + muu4**2 - ured) + 
     &                             2*ured)))))*dconjg(FR(indj,i,tt,gen))) + 
     &           dconjg(ELc(indi,i,tt,gen))*(-((g1L*(g2R*muu1**2*muu2*
     &                       (g3R*muu4*(2*g4R*mGlr*muu3 + g4L*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)) + 
     &                         g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2) + g4R*muu3*(muu1**2 + 3*muu4**2 - ured))) + 
     &                      g2L*M4Lr*(g3R*g4L*muu4*(muu1**2*(-1 - muu3**2 + muu4**2) + (-1 + muu2**2)*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                         g3L*g4R*muu1**2*muu3*(-1 + muu1**2 + muu3**2 - ured))) + 
     &                   g1R*muu1*(g2R*M4Lr*muu2*(g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) + g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                      g2L*(g3L*(g4R*muu3*(muu1**4 + (-1 + muu2**2)*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2 - ured)) + 
     &                            g4L*mGlr*(1 - muu3**2 + muu2**2*(-1 + muu3**2) - muu4**2 + muu1**2*(-1 + muu4**2) + ured)) + 
     &                         g3R*muu4*(g4R*mGlr*(-1 + muu1**2 + muu2**2)*muu3 + 
     &                            g4L*(2 + muu2**4 - muu3**2 - 2*muu4**2 - muu1**2*(2 + muu3**2 - 2*muu4**2) + muu2**2*(-3 + muu3**2 + muu4**2 - ured) + 
     &                               2*ured)))))*dconjg(FL(indj,i,tt,gen))) + 
     &              (g1R*muu1*(g2L*muu2*(g3R*(g4L*(muu1**2*(muu2**2 - muu3**2 - ured) + (-1 + muu2**2 + muu4**2 - ured)*(2*muu4**2 - ured)) + 
     &                          g4R*mGlr*muu3*(muu1**2 + muu4**2 - ured)) + 
     &                       g3L*muu4*(g4R*muu3*(3*muu1**2 + muu4**2 - ured) + g4L*mGlr*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                    g2R*M4Lr*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                       g3R*g4L*(muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured)))) + 
     &                 g1L*(g2L*M4Lr*muu2*(2*g3L*g4R*muu1**2*muu3*muu4 + 
     &                       g3R*g4L*(muu1**2*(muu2**2 - muu3**2 - ured) + (muu4**2 - ured)*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                    g2R*muu1**2*(g3L*muu4*(g4L*mGlr*(muu2**2 + muu3**2 - ured) + g4R*muu3*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured)) + 
     &                       g3R*(g4R*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + 
     &                          g4L*(2*muu3**2*muu4**2 + muu1**2*(muu2**2 - ured) + muu2**2*(muu4**2 - ured) + ured - muu3**2*ured - 2*muu4**2*ured + 
     &                             ured**2)))))*dconjg(FR(indj,i,tt,gen)))))/2. + 
     &      (xd1*(dconjg(ER(indi,i,tt,gen))*((g1L*muu1*muu4*(g2R*muu2*
     &                     (g3L*muu4*(2*g4L*mGlr*muu3 + g4R*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)) + 
     &                       g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2) + g4L*muu3*(muu1**2 + 3*muu4**2 - ured))) + 
     &                    g2L*M4Lr*(g3R*g4L*muu3*(-1 + muu1**2 + muu3**2 - ured) + g3L*g4R*muu4*(muu2**2 + muu3**2 - ured))) + 
     &                 g1R*(g2R*M4Lr*muu2*muu4*(g3R*g4L*muu3*(muu1**2 + muu4**2 - ured) + g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                    g2L*(g3L*(g4L*mGlr*muu3*(muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 - 2*ured) + ured - muu3**2*ured - 
     &                             muu4**2*ured + ured**2) + g4R*muu4**2*
     &                           (muu1**2*(2*muu2**2 + muu3**2 - 2*ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))) + 
     &                       g3R*muu4*(g4L*muu3*(2*muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-2 + 2*muu3**2 - 3*ured) + ured - muu3**2*ured - 
     &                             muu4**2*ured + ured**2) + g4R*mGlr*
     &                           (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))))))*dconjg(FL(indj,i,tt,gen))
     &               - (g1R*muu4*(g2L*muu2*(g3L*(g4R*(muu1**2*(-1 + muu3**2 + muu4**2) + muu4**2*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                          g4L*mGlr*muu3*(muu1**2 + muu4**2 - ured)) + 
     &                       g3R*muu4*(g4L*muu3*(3*muu1**2 + muu4**2 - ured) + g4R*mGlr*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                    g2R*M4Lr*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                       g3L*g4R*(1 - muu3**2 + muu2**2*(-1 + muu3**2) - muu4**2 + muu1**2*(-1 + muu4**2) + ured))) + 
     &                 g1L*muu1*(g2L*M4Lr*muu2*muu4*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &                    g2R*(g3L*muu4*(g4L*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + 
     &                          g4R*(2 - 3*muu3**2 + muu3**4 - 2*muu4**2 + muu2**2*(-1 + muu3**2 - muu4**2) + muu1**2*(-2 + muu3**2 + 2*muu4**2) + 
     &                             2*ured - muu3**2*ured)) + g3R*
     &                        (g4L*muu3*muu4**2*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured) + 
     &                          g4R*mGlr*(1 + muu3**4 - muu4**2 - muu2**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + ured - muu3**2*(2 + ured))))))*
     &               dconjg(FR(indj,i,tt,gen))) + dconjg(ELc(indi,i,tt,gen))*
     &            (-((g1L*muu4*(g2R*muu2*(g3R*(g4L*(muu1**2*(-1 + muu3**2 + muu4**2) + muu4**2*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                            g4R*mGlr*muu3*(muu1**2 + muu4**2 - ured)) + 
     &                         g3L*muu4*(g4R*muu3*(3*muu1**2 + muu4**2 - ured) + g4L*mGlr*(-1 + muu2**2 + muu4**2 - ured))) + 
     &                      g2L*M4Lr*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                         g3R*g4L*(1 - muu3**2 + muu2**2*(-1 + muu3**2) - muu4**2 + muu1**2*(-1 + muu4**2) + ured))) + 
     &                   g1R*muu1*(g2R*M4Lr*muu2*muu4*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &                      g2L*(g3R*muu4*(g4R*mGlr*muu3*(-1 + muu1**2 + muu3**2 - ured) + 
     &                            g4L*(2 - 3*muu3**2 + muu3**4 - 2*muu4**2 + muu2**2*(-1 + muu3**2 - muu4**2) + muu1**2*(-2 + muu3**2 + 2*muu4**2) + 
     &                               2*ured - muu3**2*ured)) + 
     &                         g3L*(g4R*muu3*muu4**2*(-2 + 2*muu1**2 + muu2**2 + muu3**2 - ured) + 
     &                            g4L*mGlr*(1 + muu3**4 - muu4**2 - muu2**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + ured - muu3**2*(2 + ured))))))*
     &                 dconjg(FL(indj,i,tt,gen))) + (g1R*muu1*muu4*
     &                  (g2L*muu2*(g3R*muu4*(2*g4R*mGlr*muu3 + g4L*(-2 + muu2**2 + muu3**2 + 2*muu4**2 - ured)) + 
     &                       g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2) + g4R*muu3*(muu1**2 + 3*muu4**2 - ured))) + 
     &                    g2R*M4Lr*(g3L*g4R*muu3*(-1 + muu1**2 + muu3**2 - ured) + g3R*g4L*muu4*(muu2**2 + muu3**2 - ured))) + 
     &                 g1L*(g2L*M4Lr*muu2*muu4*(g3L*g4R*muu3*(muu1**2 + muu4**2 - ured) + g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 - ured)) + 
     &                    g2R*(g3R*(g4R*mGlr*muu3*(muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 - 2*ured) + ured - muu3**2*ured - 
     &                             muu4**2*ured + ured**2) + g4L*muu4**2*
     &                           (muu1**2*(2*muu2**2 + muu3**2 - 2*ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))) + 
     &                       g3L*muu4*(g4R*muu3*(2*muu1**4 - muu2**2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-2 + 2*muu3**2 - 3*ured) + ured - muu3**2*ured - 
     &                             muu4**2*ured + ured**2) + g4L*mGlr*
     &                           (muu1**2*(muu2**2 - ured) + muu3**2*(muu4**2 - ured) + ured*(1 - muu2**2 - muu4**2 + ured))))))*dconjg(FR(indj,i,tt,gen)))
     &           ))/2. + (xd2*(dconjg(ER(indi,i,tt,gen))*((g1L*muu1*
     &                  (g2L*M4Lr*(g3R*g4L*(-1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                       g3L*g4R*(-muu3**4 + muu4**2*(muu2**2 - ured) + muu3**2*(1 - muu1**2 + muu4**2 + ured))) + 
     &                    g2R*muu2*(g3R*muu4*(-(g4R*mGlr*(1 + muu3**2 - muu4**2)) + g4L*muu3*(2 + muu1**2 - muu2**2 - muu3**2 + muu4**2)) + 
     &                       g3L*(g4L*mGlr*muu3*(1 - muu3**2 + muu4**2) + 
     &                          g4R*(-(muu1**2*muu3**2) + 2*muu4**4 + muu4**2*(-2 + muu2**2 - 2*muu3**2 - ured) + muu3**2*ured)))) + 
     &                 g1R*(g2R*M4Lr*muu2*(g3R*g4L*(1 + muu1**2 - muu2**2)*muu3*muu4 - 
     &                       g3L*g4R*(muu1**2*muu3**2 - muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured))) + 
     &                    g2L*(g3R*muu4*(g4L*muu3*(2*muu1**4 - muu1**2*(2 + 2*muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                          g4R*mGlr*(-muu2**4 + muu1**2*(muu3**2 - ured) + muu2**2*(1 + muu1**2 - muu4**2 + ured))) + 
     &                       g3L*(g4L*mGlr*muu3*(muu1**4 - muu1**2*(1 + muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                          g4R*(-(muu1**4*muu3**2) - muu2**2*muu4**2*(-1 + muu2**2 + muu4**2 - ured) + 
     &                             muu1**2*(-muu3**4 + 2*muu4**2*(muu2**2 - ured) + muu3**2*(1 + 2*muu4**2 + ured)))))))*dconjg(FL(indj,i,tt,gen)) + 
     &              (g1R*(g2L*muu2*(g3R*(g4L*muu3*((-1 + muu2**2)*muu4**2 + muu1**2*(-1 + muu3**2 - 2*muu4**2)) - 
     &                          g4R*mGlr*(-1 + (1 + muu1**2)*muu3**2 - muu2**2*(-1 + muu3**2) + muu4**2 - ured)) + 
     &                       g3L*muu4*(g4L*mGlr*(-1 - muu1**2 + muu2**2)*muu3 - 
     &                          g4R*(-1 + muu2**2 + muu3**2 - muu2**2*muu3**2 + muu4**2 + muu1**2*(-1 + muu4**2) - ured))) + 
     &                    g2R*M4Lr*(g3L*g4R*muu4*(-1 + muu2**2 + muu4**2 + muu1**2*(1 + muu3**2 - muu4**2) - ured) + 
     &                       g3R*g4L*muu3*(1 - muu1**2 - muu3**2 + muu2**2*(-1 + muu3**2 - muu4**2) + ured))) + 
     &                 g1L*muu1*(g2L*M4Lr*muu2*(g3R*g4L*muu3*(-1 + muu3**2 - muu4**2) + g3L*g4R*muu4*(1 + muu3**2 - muu4**2)) + 
     &                    g2R*(g3L*muu4*(g4L*mGlr*(1 - muu1**2 + muu2**2)*muu3 + 
     &                          g4R*(-2 + muu3**2 + 2*muu4**2 + muu1**2*(2 + muu3**2 - 2*muu4**2) + muu2**2*(1 + muu4**2) - 2*ured)) + 
     &                       g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2 + muu2**2*muu4**2 - muu1**2*(-1 + muu4**2) - ured) + 
     &                          g4L*muu3*(1 - muu3**2 + muu2**2*(-1 + muu3**2) + muu4**2 - muu1**2*(1 + muu4**2) + ured)))))*dconjg(FR(indj,i,tt,gen))) + 
     &           dconjg(ELc(indi,i,tt,gen))*((g1L*(g2R*muu2*(g3L*
     &                        (g4R*muu3*((-1 + muu2**2)*muu4**2 + muu1**2*(-1 + muu3**2 - 2*muu4**2)) - 
     &                          g4L*mGlr*(-1 + (1 + muu1**2)*muu3**2 - muu2**2*(-1 + muu3**2) + muu4**2 - ured)) + 
     &                       g3R*muu4*(g4R*mGlr*(-1 - muu1**2 + muu2**2)*muu3 - 
     &                          g4L*(-1 + muu2**2 + muu3**2 - muu2**2*muu3**2 + muu4**2 + muu1**2*(-1 + muu4**2) - ured))) + 
     &                    g2L*M4Lr*(g3R*g4L*muu4*(-1 + muu2**2 + muu4**2 + muu1**2*(1 + muu3**2 - muu4**2) - ured) + 
     &                       g3L*g4R*muu3*(1 - muu1**2 - muu3**2 + muu2**2*(-1 + muu3**2 - muu4**2) + ured))) + 
     &                 g1R*muu1*(g2R*M4Lr*muu2*(g3L*g4R*muu3*(-1 + muu3**2 - muu4**2) + g3R*g4L*muu4*(1 + muu3**2 - muu4**2)) + 
     &                    g2L*(g3R*muu4*(g4R*mGlr*(1 - muu1**2 + muu2**2)*muu3 + 
     &                          g4L*(-2 + muu3**2 + 2*muu4**2 + muu1**2*(2 + muu3**2 - 2*muu4**2) + muu2**2*(1 + muu4**2) - 2*ured)) + 
     &                       g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2 + muu2**2*muu4**2 - muu1**2*(-1 + muu4**2) - ured) + 
     &                          g4R*muu3*(1 - muu3**2 + muu2**2*(-1 + muu3**2) + muu4**2 - muu1**2*(1 + muu4**2) + ured)))))*dconjg(FL(indj,i,tt,gen)) + 
     &              (g1R*muu1*(g2R*M4Lr*(g3L*g4R*(-1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                       g3R*g4L*(-muu3**4 + muu4**2*(muu2**2 - ured) + muu3**2*(1 - muu1**2 + muu4**2 + ured))) + 
     &                    g2L*muu2*(g3L*muu4*(-(g4L*mGlr*(1 + muu3**2 - muu4**2)) + g4R*muu3*(2 + muu1**2 - muu2**2 - muu3**2 + muu4**2)) + 
     &                       g3R*(g4R*mGlr*muu3*(1 - muu3**2 + muu4**2) + 
     &                          g4L*(-(muu1**2*muu3**2) + 2*muu4**4 + muu4**2*(-2 + muu2**2 - 2*muu3**2 - ured) + muu3**2*ured)))) + 
     &                 g1L*(g2L*M4Lr*muu2*(g3L*g4R*(1 + muu1**2 - muu2**2)*muu3*muu4 - 
     &                       g3R*g4L*(muu1**2*muu3**2 - muu4**4 - muu3**2*ured + muu4**2*(1 - muu2**2 + muu3**2 + ured))) + 
     &                    g2R*(g3L*muu4*(g4R*muu3*(2*muu1**4 - muu1**2*(2 + 2*muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                          g4L*mGlr*(-muu2**4 + muu1**2*(muu3**2 - ured) + muu2**2*(1 + muu1**2 - muu4**2 + ured))) + 
     &                       g3R*(g4R*mGlr*muu3*(muu1**4 - muu1**2*(1 + muu2**2 - muu3**2 + ured) + muu2**2*(-muu4**2 + ured)) + 
     &                          g4L*(-(muu1**4*muu3**2) - muu2**2*muu4**2*(-1 + muu2**2 + muu4**2 - ured) + 
     &                             muu1**2*(-muu3**4 + 2*muu4**2*(muu2**2 - ured) + muu3**2*(1 + 2*muu4**2 + ured)))))))*dconjg(FR(indj,i,tt,gen)))))/2.))/
     &  (4.*Pi**2*(-tred + MSf(i,itypeT,gen)**2/sqrtS**2))
	
	tempBoxt = tempBoxt + ChiChi2QQ_M2BoxGluino - tempBox2
	tempBox2 = ChiChi2QQ_M2BoxGluino
	
	enddo

	enddo
	enddo


!	write(*,*)'ChiChi U-Channel Gluino Box (Neu) = ', ChiChi2QQ_M2BoxGluino
!	write(*,*)'ChiChi U-Channel Gluino Box (Neu) mit Z = ', tempBoxZ
!	write(*,*)'ChiChi U-Channel Gluino Box (Neu) mit H = ', tempBoxH
!	write(*,*)'ChiChi U-Channel Gluino Box (Neu) mit t = ', tempBoxt
!	write(*,*)'ChiChi U-Channel Gluino Box (Neu) mit u = ', tempBoxu
	tempBox = ChiChi2QQ_M2BoxGluino
	tempBoxZ = 0d0
	tempBoxH = 0d0
	tempBoxt = 0d0 		
	tempBoxu = 0d0



!************************************
!******* Gluino T-channel Box *******
!************************************

	  do k = 1, 2
	    do l = 1, 2


		if(icha==0)		then
				g1L = FL(indi, l, tt, gen)
				g1R = FR(indi, l, tt, gen)
				g2L = ELc(indj, k, tt, gen)
				g2R = ER(indj, k, tt, gen)
		else if(icha==1) then
				g1L = GL(indi, l, tt, gen)
				g1R = GR(indi, l, tt, gen)
				g2L = dconjg(FR(indj, k, tt, gen))
				g2R = dconjg(FL(indj, k, tt, gen))
		else if(icha==2) then
! couplings set to G & H of the bottom final state		      
				if(tt==3)	then
				   g1L = -EL/SW*VChaC(indi,1)*USf(l,1,tt,gen) + EL*Mf(tt,gen)/(SW*sqrt2*MW*SB)*VChaC(indi,2)*USf(l,2,tt,gen)
				   g1R =  EL*Mf(tt+1,gen)/(SW*sqrt2*MW*CB)*UCha(indi,2)*USf(l,1,tt,gen)
				   g2L =  EL*Mf(tt+1,gen)/(SW*sqrt2*MW*CB)*UChaC(indj,2)*USf(k,1,tt,gen)
				   g2R = -EL/SW*VCha(indj,1)*USf(k,1,tt,gen) + EL*Mf(tt,gen)/(SW*sqrt2*MW*SB)*VCha(indj,2)*USf(k,2,tt,gen)
				else
					g1L = 0d0
					g1R = 0d0
					g2L = 0d0
					g2R = 0d0
				endif
		endif
		
	    g3L = sqrt2*GS*USf(l,2,tt,gen)
	    g3R = -sqrt2*GS*USf(l,1,tt,gen)
	    g4L = -sqrt2*GS*USf(k,1,att,gen)
	    g4R = sqrt2*GS*USf(k,2,att,gen)

	    M4Lr = Mf(itypeT,gen)/sqrtS
	    M4L  = Mf(itypeT,gen)

	    call myDget(muu3**2*sqrtS**2,muu4**2*sqrtS**2,muu2**2*sqrtS**2,muu1**2*sqrtS**2,sqrtS**2,tred*sqrtS**2,MSf2(l,tt,gen),MGl2,MSf2(k,att,gen),Mf2(itypeT,gen))
			

! Photon contr.

	if (icha.eq.2) then
	
	ChiChi2QQ_M2BoxGluino = ChiChi2QQ_M2BoxGluino + (BpL*CF*NC*sqrtS**4*(2*muu1*(g1L*g2L*M4Lr + g1R*g2R*M4Lr + g1R*g2L*muu1 + g1L*g2R*muu1)*muu3*
     -       (g3L*g4L*mGlr + g3R*g4R*mGlr + g3R*g4L*muu3 + g3L*g4R*muu3)*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred)*xd0 + 
     -      (2*(g1L*g2R*(g3L*g4R*(muu1**4 + muu3**4 + muu3**2*(1 - 2*tred) + muu1**2*(1 + 2*muu3**2 - 2*tred) + tred**2) + 
     -              g3R*g4L*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*tred) + (1 + tred)**2 - muu3**2*(1 + 2*tred))) + 
     -           g1R*g2L*(g3R*g4L*(muu1**4 + muu3**4 + muu3**2*(1 - 2*tred) + muu1**2*(1 + 2*muu3**2 - 2*tred) + tred**2) + 
     -              g3L*g4R*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*tred) + (1 + tred)**2 - muu3**2*(1 + 2*tred))))*xd00)/sqrtS**2 + 
     -      2*muu3*(g1R*(g2R*(g3R*g4L + g3L*g4R)*M4Lr*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred) + 
     -            g2L*(g3L*(g4R*muu3*(3*muu1**4 + muu3**4 + muu1**2*(-2 + 4*muu3**2 - 4*tred) + tred - 2*muu3**2*tred + tred**2) + 
     -                  g4L*mGlr*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*tred) + tred - 2*muu3**2*tred + tred**2)) + 
     -               g3R*(g4L*muu3*(3*muu1**4 + muu3**4 + muu1**2*(-2 + 4*muu3**2 - 4*tred) + tred - 2*muu3**2*tred + tred**2) + 
     -                  g4R*mGlr*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*tred) + tred - 2*muu3**2*tred + tred**2)))) + 
     -         g1L*(g2L*(g3R*g4L + g3L*g4R)*M4Lr*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred) + 
     -            g2R*(g3L*(g4R*muu3*(3*muu1**4 + muu3**4 + muu1**2*(-2 + 4*muu3**2 - 4*tred) + tred - 2*muu3**2*tred + tred**2) + 
     -                  g4L*mGlr*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*tred) + tred - 2*muu3**2*tred + tred**2)) + 
     -               g3R*(g4L*muu3*(3*muu1**4 + muu3**4 + muu1**2*(-2 + 4*muu3**2 - 4*tred) + tred - 2*muu3**2*tred + tred**2) + 
     -                  g4R*mGlr*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*tred) + tred - 2*muu3**2*tred + tred**2)))))*xd1 + 
     -      2*(g1R*g2L + g1L*g2R)*(g3R*g4L + g3L*g4R)*muu3**2*(muu1**4 + muu3**4 + muu1**2*(-1 + 2*muu3**2 - 2*tred) + tred - 2*muu3**2*tred + tred**2)*
     -       xd11 + (g1R*g2L*(g3L*g4R*(muu1**6 + muu3**6 + muu1**4*(-1 + 7*muu3**2 - 3*tred) + 3*muu3**2*tred*(1 + tred) - tred*(1 + tred)**2 - 
     -               muu3**4*(1 + 3*tred) + muu1**2*(7*muu3**4 + 3*tred*(1 + tred) - 2*muu3**2*(1 + 5*tred))) + 
     -            g3R*g4L*(muu1**6 + muu3**6 + muu1**4*(7*muu3**2 - 3*tred) - 3*muu3**4*tred - tred**2*(1 + tred) + muu3**2*tred*(1 + 3*tred) + 
     -               muu1**2*(7*muu3**4 + tred*(1 + 3*tred) - 2*muu3**2*(2 + 5*tred)))) + 
     -         g1L*g2R*(g3R*g4L*(muu1**6 + muu3**6 + muu1**4*(-1 + 7*muu3**2 - 3*tred) + 3*muu3**2*tred*(1 + tred) - tred*(1 + tred)**2 - 
     -               muu3**4*(1 + 3*tred) + muu1**2*(7*muu3**4 + 3*tred*(1 + tred) - 2*muu3**2*(1 + 5*tred))) + 
     -            g3L*g4R*(muu1**6 + muu3**6 + muu1**4*(7*muu3**2 - 3*tred) - 3*muu3**4*tred - tred**2*(1 + tred) + muu3**2*tred*(1 + 3*tred) + 
     -               muu1**2*(7*muu3**4 + tred*(1 + 3*tred) - 2*muu3**2*(2 + 5*tred)))))*xd13 + 
     -      2*muu1*(g1R*(g2R*(g3R*g4L + g3L*g4R)*M4Lr*(muu1**4 + muu3**4 + 2*muu1**2*(muu3**2 - tred) + tred*(1 + tred) - muu3**2*(1 + 2*tred)) + 
     -            g2L*muu1*(g3R*(g4R*mGlr*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred) + 
     -                  g4L*(muu1**4 + 3*muu3**4 + muu1**2*(4*muu3**2 - 2*tred) + tred*(1 + tred) - 2*muu3**2*(1 + 2*tred))) + 
     -               g3L*(g4L*mGlr*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred) + 
     -                  g4R*(muu1**4 + 3*muu3**4 + muu1**2*(4*muu3**2 - 2*tred) + tred*(1 + tred) - 2*muu3**2*(1 + 2*tred))))) + 
     -         g1L*(g2L*(g3R*g4L + g3L*g4R)*M4Lr*(muu1**4 + muu3**4 + 2*muu1**2*(muu3**2 - tred) + tred*(1 + tred) - muu3**2*(1 + 2*tred)) + 
     -            g2R*muu1*(g3R*(g4R*mGlr*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred) + 
     -                  g4L*(muu1**4 + 3*muu3**4 + muu1**2*(4*muu3**2 - 2*tred) + tred*(1 + tred) - 2*muu3**2*(1 + 2*tred))) + 
     -               g3L*(g4L*mGlr*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred) + 
     -                  g4R*(muu1**4 + 3*muu3**4 + muu1**2*(4*muu3**2 - 2*tred) + tred*(1 + tred) - 2*muu3**2*(1 + 2*tred))))))*xd3 + 
     -      2*(g1R*g2L + g1L*g2R)*(g3R*g4L + g3L*g4R)*muu1**2*(muu1**4 + muu3**4 + 2*muu1**2*(muu3**2 - tred) + tred*(1 + tred) - muu3**2*(1 + 2*tred))*
     -       xd33)*dconjg(ApL))/(4.*Pi**2)
	
	endif

! Zboson contr.
	
	ChiChi2QQ_M2BoxGluino = ChiChi2QQ_M2BoxGluino -(CF*NC*sqrtS**6*(xd0*(-((BL*(g1R*muu1*(g2L*muu2*(g3L*muu3*(g4R*muu4*(muu1**2 + muu3**2 - tred) + g4L*mGlr*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                        g3R*(g4R*mGlr*muu4*(muu1**2 + muu3**2 - tred) + g4L*muu3**2*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                     g2R*M4Lr*(g3L*muu3*(g4R*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4L*mGlr*(muu2**2 + muu4**2 - tred)) + 
     &                        g3R*(g4R*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4L*muu3**2*(muu2**2 + muu4**2 - tred)))) + 
     &                  g1L*(g2L*M4Lr*muu2*(g3L*muu3*(g4R*muu4*(muu1**2 + muu3**2 - tred) + g4L*mGlr*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                        g3R*(g4R*mGlr*muu4*(muu1**2 + muu3**2 - tred) + g4L*muu3**2*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                     g2R*muu1**2*(g3L*muu3*(g4R*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4L*mGlr*(muu2**2 + muu4**2 - tred)) + 
     &                        g3R*(g4R*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4L*muu3**2*(muu2**2 + muu4**2 - tred))))) + 
     &               BR*(g1R*muu1*(g2L*muu2*(g3R*muu3*(g4L*muu4*(muu1**2 + muu3**2 - tred) + g4R*mGlr*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                        g3L*(g4L*mGlr*muu4*(muu1**2 + muu3**2 - tred) + g4R*muu3**2*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                     g2R*M4Lr*(g3R*muu3*(g4L*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4R*mGlr*(muu2**2 + muu4**2 - tred)) + 
     &                        g3L*(g4L*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4R*muu3**2*(muu2**2 + muu4**2 - tred)))) + 
     &                  g1L*(g2L*M4Lr*muu2*(g3R*muu3*(g4L*muu4*(muu1**2 + muu3**2 - tred) + g4R*mGlr*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                        g3L*(g4L*mGlr*muu4*(muu1**2 + muu3**2 - tred) + g4R*muu3**2*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                     g2R*muu1**2*(g3R*muu3*(g4L*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4R*mGlr*(muu2**2 + muu4**2 - tred)) + 
     &                        g3L*(g4L*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4R*muu3**2*(muu2**2 + muu4**2 - tred))))))*dconjg(AL)) - 
     &          (BL*(g1L*muu1*(g2R*muu2*(g3L*muu3*(g4R*muu4*(muu1**2 + muu3**2 - tred) + g4L*mGlr*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                      g3R*(g4R*mGlr*muu4*(muu1**2 + muu3**2 - tred) + g4L*muu3**2*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                   g2L*M4Lr*(g3L*muu3*(g4R*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4L*mGlr*(muu2**2 + muu4**2 - tred)) + 
     &                      g3R*(g4R*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4L*muu3**2*(muu2**2 + muu4**2 - tred)))) + 
     &                g1R*(g2R*M4Lr*muu2*(g3L*muu3*(g4R*muu4*(muu1**2 + muu3**2 - tred) + g4L*mGlr*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                      g3R*(g4R*mGlr*muu4*(muu1**2 + muu3**2 - tred) + g4L*muu3**2*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                   g2L*muu1**2*(g3L*muu3*(g4R*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4L*mGlr*(muu2**2 + muu4**2 - tred)) + 
     &                      g3R*(g4R*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4L*muu3**2*(muu2**2 + muu4**2 - tred))))) + 
     &             BR*(g1L*muu1*(g2R*muu2*(g3R*muu3*(g4L*muu4*(muu1**2 + muu3**2 - tred) + g4R*mGlr*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                      g3L*(g4L*mGlr*muu4*(muu1**2 + muu3**2 - tred) + g4R*muu3**2*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                   g2L*M4Lr*(g3R*muu3*(g4L*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4R*mGlr*(muu2**2 + muu4**2 - tred)) + 
     &                      g3L*(g4L*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4R*muu3**2*(muu2**2 + muu4**2 - tred)))) + 
     &                g1R*(g2R*M4Lr*muu2*(g3R*muu3*(g4L*muu4*(muu1**2 + muu3**2 - tred) + g4R*mGlr*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                      g3L*(g4L*mGlr*muu4*(muu1**2 + muu3**2 - tred) + g4R*muu3**2*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                   g2L*muu1**2*(g3R*muu3*(g4L*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4R*mGlr*(muu2**2 + muu4**2 - tred)) + 
     &                      g3L*(g4L*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4R*muu3**2*(muu2**2 + muu4**2 - tred))))))*dconjg(AR)) + 
     &       (2*xd00*((BR*(g1R*g2L*muu1*muu2*(-4*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &                  g1L*g2R*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 - g3L*g4R*(muu1**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred))) + 
     &               BL*(g1R*g2L*muu1*muu2*(-4*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &                  g1L*g2R*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 + g3R*g4L*(-1 + muu2**2 + muu3**2 - tred)*(1 - muu1**2 - muu4**2 + tred))))*
     &             dconjg(AL) + (BL*(g1L*g2R*muu1*muu2*(-4*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &                  g1R*g2L*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 - g3R*g4L*(muu1**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred))) + 
     &               BR*(g1L*g2R*muu1*muu2*(-4*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &                  g1R*g2L*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 + g3L*g4R*(-1 + muu2**2 + muu3**2 - tred)*(1 - muu1**2 - muu4**2 + tred))))*
     &             dconjg(AR)))/sqrtS**2 - muu1*xd33*((BL*(g1R*g2L*muu2*
     &                 (2*g3L*g4R*muu1**2*muu3*muu4 + g3R*g4L*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))) + 
     &                g1L*g2R*muu1*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3R*g4L*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))) + 
     &             BR*(g1R*g2L*muu2*(2*g3R*g4L*muu1**2*muu3*muu4 + 
     &                   g3L*g4R*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))) + 
     &                g1L*g2R*muu1*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3L*g4R*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))))*dconjg(AL) + 
     &          (BL*(g1L*g2R*muu2*(2*g3L*g4R*muu1**2*muu3*muu4 + 
     &                   g3R*g4L*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))) + 
     &                g1R*g2L*muu1*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3R*g4L*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))) + 
     &             BR*(g1L*g2R*muu2*(2*g3R*g4L*muu1**2*muu3*muu4 + 
     &                   g3L*g4R*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))) + 
     &                g1R*g2L*muu1*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3L*g4R*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))))*dconjg(AR)) - 
     &       muu3*xd11*((BL*(g1R*g2L*muu1*muu2*muu3*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &                g1L*g2R*(g3L*g4R*muu4*(muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu4**2 - 2*tred) + tred - muu3**2*tred - 
     &                      muu4**2*tred + tred**2) + g3R*g4L*muu3*
     &                    (muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))) + 
     &             BR*(g1R*g2L*muu1*muu2*muu3*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &                g1L*g2R*(g3R*g4L*muu4*(muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu4**2 - 2*tred) + tred - muu3**2*tred - 
     &                      muu4**2*tred + tred**2) + g3L*g4R*muu3*
     &                    (muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))))*dconjg(AL) + 
     &          (BL*(g1L*g2R*muu1*muu2*muu3*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &                g1R*g2L*(g3L*g4R*muu4*(muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu4**2 - 2*tred) + tred - muu3**2*tred - 
     &                      muu4**2*tred + tred**2) + g3R*g4L*muu3*
     &                    (muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))) + 
     &             BR*(g1L*g2R*muu1*muu2*muu3*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &                g1R*g2L*(g3R*g4L*muu4*(muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu4**2 - 2*tred) + tred - muu3**2*tred - 
     &                      muu4**2*tred + tred**2) + g3L*g4R*muu3*
     &                    (muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))))*dconjg(AR)) + 
     &       xd23*((BL*(2*g1R*g2L*muu1*muu2*(g3L*g4R*(-1 - muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3R*g4L*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 - muu2**2 + muu4**2 + tred))) + 
     &                g1L*g2R*(2*g3L*g4R*muu1**2*(1 - muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3R*g4L*(muu2**2*(muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**4*(-muu2**2 + muu4**2 + tred) + 
     &                      muu1**2*(muu2**4 - (1 + 2*muu3**2)*muu4**2 + muu4**4 - muu2**2*(muu3**2 + muu4**2) + (-1 + 2*muu3**2 - tred)*tred)))) + 
     &             BR*(2*g1R*g2L*muu1*muu2*(g3R*g4L*(-1 - muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3L*g4R*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 - muu2**2 + muu4**2 + tred))) + 
     &                g1L*g2R*(2*g3R*g4L*muu1**2*(1 - muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3L*g4R*(muu2**2*(muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**4*(-muu2**2 + muu4**2 + tred) + 
     &                      muu1**2*(muu2**4 - (1 + 2*muu3**2)*muu4**2 + muu4**4 - muu2**2*(muu3**2 + muu4**2) + (-1 + 2*muu3**2 - tred)*tred)))))*
     &           dconjg(AL) + (BL*(2*g1L*g2R*muu1*muu2*(g3L*g4R*(-1 - muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3R*g4L*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 - muu2**2 + muu4**2 + tred))) + 
     &                g1R*g2L*(2*g3L*g4R*muu1**2*(1 - muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3R*g4L*(muu2**2*(muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**4*(-muu2**2 + muu4**2 + tred) + 
     &                      muu1**2*(muu2**4 - (1 + 2*muu3**2)*muu4**2 + muu4**4 - muu2**2*(muu3**2 + muu4**2) + (-1 + 2*muu3**2 - tred)*tred)))) + 
     &             BR*(2*g1L*g2R*muu1*muu2*(g3R*g4L*(-1 - muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3L*g4R*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 - muu2**2 + muu4**2 + tred))) + 
     &                g1R*g2L*(2*g3R*g4L*muu1**2*(1 - muu1**2 + muu2**2)*muu3*muu4 + 
     &                   g3L*g4R*(muu2**2*(muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**4*(-muu2**2 + muu4**2 + tred) + 
     &                      muu1**2*(muu2**4 - (1 + 2*muu3**2)*muu4**2 + muu4**4 - muu2**2*(muu3**2 + muu4**2) + (-1 + 2*muu3**2 - tred)*tred)))))*
     &           dconjg(AR)) + xd22*(-((BL*(g1R*g2L*muu1*muu2*
     &                   (2*g3L*g4R*muu3*muu4 + g3R*g4L*(muu3**4 + muu4**2*(-1 + muu4**2) - muu3**2*(1 + 2*muu4**2))) + 
     &                  g1L*g2R*(g3L*g4R*(muu1**4 + muu2**2*(-1 + muu2**2) - muu1**2*(1 + 2*muu2**2))*muu3*muu4 + 
     &                     g3R*g4L*(-(muu2**4*muu3**2) + muu1**2*(-muu4**4 - muu3**2*tred + muu4**2*(1 - muu1**2 + muu3**2 + tred)) + 
     &                        muu2**2*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 + muu1**2 + muu4**2 + tred))))) + 
     &               BR*(g1R*g2L*muu1*muu2*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(muu3**4 + muu4**2*(-1 + muu4**2) - muu3**2*(1 + 2*muu4**2))) + 
     &                  g1L*g2R*(g3R*g4L*(muu1**4 + muu2**2*(-1 + muu2**2) - muu1**2*(1 + 2*muu2**2))*muu3*muu4 + 
     &                     g3L*g4R*(-(muu2**4*muu3**2) + muu1**2*(-muu4**4 - muu3**2*tred + muu4**2*(1 - muu1**2 + muu3**2 + tred)) + 
     &                        muu2**2*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 + muu1**2 + muu4**2 + tred))))))*dconjg(AL)) - 
     &          (BL*(g1L*g2R*muu1*muu2*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(muu3**4 + muu4**2*(-1 + muu4**2) - muu3**2*(1 + 2*muu4**2))) + 
     &                g1R*g2L*(g3L*g4R*(muu1**4 + muu2**2*(-1 + muu2**2) - muu1**2*(1 + 2*muu2**2))*muu3*muu4 + 
     &                   g3R*g4L*(-(muu2**4*muu3**2) + muu1**2*(-muu4**4 - muu3**2*tred + muu4**2*(1 - muu1**2 + muu3**2 + tred)) + 
     &                      muu2**2*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 + muu1**2 + muu4**2 + tred))))) + 
     &             BR*(g1L*g2R*muu1*muu2*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(muu3**4 + muu4**2*(-1 + muu4**2) - muu3**2*(1 + 2*muu4**2))) + 
     &                g1R*g2L*(g3R*g4L*(muu1**4 + muu2**2*(-1 + muu2**2) - muu1**2*(1 + 2*muu2**2))*muu3*muu4 + 
     &                   g3L*g4R*(-(muu2**4*muu3**2) + muu1**2*(-muu4**4 - muu3**2*tred + muu4**2*(1 - muu1**2 + muu3**2 + tred)) + 
     &                      muu2**2*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 + muu1**2 + muu4**2 + tred))))))*dconjg(AR)) + 
     &       xd13*(-((BR*(2*g1R*g2L*muu1*muu2*muu3*(g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + g3L*g4R*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                  g1L*g2R*(2*g3R*g4L*muu1**2*muu3*muu4*(-1 + muu1**2 + muu4**2 - tred) + 
     &                     g3L*g4R*(muu1**2 + muu3**2 - tred)*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))))
     &                 + BL*(2*g1R*g2L*muu1*muu2*muu3*(g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + g3R*g4L*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                  g1L*g2R*(2*g3L*g4R*muu1**2*muu3*muu4*(-1 + muu1**2 + muu4**2 - tred) + 
     &                     g3R*g4L*(muu1**4*(muu2**2 - muu4**2 - tred) - 
     &                        (-1 + muu2**2 + muu3**2 - tred)*(muu2**2*muu3**2 + (-1 + muu4**2 - tred)*tred + muu3**2*(-muu4**2 + tred)) + 
     &                        muu1**2*((1 + 2*muu3**2)*muu4**2 - muu4**4 + muu2**2*(-1 + 2*muu3**2 + muu4**2 - 2*tred) + tred*(2 - 3*muu3**2 + 2*tred))))))
     &              *dconjg(AL)) - (BL*(2*g1L*g2R*muu1*muu2*muu3*
     &                 (g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + g3R*g4L*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                g1R*g2L*(2*g3L*g4R*muu1**2*muu3*muu4*(-1 + muu1**2 + muu4**2 - tred) + 
     &                   g3R*g4L*(muu1**2 + muu3**2 - tred)*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))))
     &              + BR*(2*g1L*g2R*muu1*muu2*muu3*(g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + g3L*g4R*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                g1R*g2L*(2*g3R*g4L*muu1**2*muu3*muu4*(-1 + muu1**2 + muu4**2 - tred) + 
     &                   g3L*g4R*(muu1**4*(muu2**2 - muu4**2 - tred) - 
     &                      (-1 + muu2**2 + muu3**2 - tred)*(muu2**2*muu3**2 + (-1 + muu4**2 - tred)*tred + muu3**2*(-muu4**2 + tred)) + 
     &                      muu1**2*((1 + 2*muu3**2)*muu4**2 - muu4**4 + muu2**2*(-1 + 2*muu3**2 + muu4**2 - 2*tred) + tred*(2 - 3*muu3**2 + 2*tred))))))*
     &           dconjg(AR)) + xd12*((BL*(-2*g1R*g2L*muu1*muu2*muu3*(g3R*g4L*muu3*(-1 + muu3**2 - muu4**2) + g3L*g4R*muu4*(1 + muu3**2 - muu4**2)) + 
     &                g1L*g2R*(2*g3L*g4R*muu3*muu4*(-muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(1 + muu2**2 - muu4**2 + tred)) + 
     &                   g3R*g4L*(muu2**4*muu3**2 + muu1**4*muu4**2 + muu2**2*muu3**2*(-1 - 2*muu1**2 + muu3**2 - muu4**2) - 
     &                      (muu3**2 - muu4**2)*(muu3**2*(muu4**2 - tred) + tred*(1 - muu4**2 + tred)) + 
     &                      muu1**2*(muu4**4 + 2*muu3**2*tred - muu4**2*(1 + muu3**2 + 2*tred))))) + 
     &             BR*(-2*g1R*g2L*muu1*muu2*muu3*(g3L*g4R*muu3*(-1 + muu3**2 - muu4**2) + g3R*g4L*muu4*(1 + muu3**2 - muu4**2)) + 
     &                g1L*g2R*(2*g3R*g4L*muu3*muu4*(-muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(1 + muu2**2 - muu4**2 + tred)) + 
     &                   g3L*g4R*(muu2**4*muu3**2 + muu1**4*muu4**2 + muu2**2*muu3**2*(-1 - 2*muu1**2 + muu3**2 - muu4**2) - 
     &                      (muu3**2 - muu4**2)*(muu3**2*(muu4**2 - tred) + tred*(1 - muu4**2 + tred)) + 
     &                      muu1**2*(muu4**4 + 2*muu3**2*tred - muu4**2*(1 + muu3**2 + 2*tred))))))*dconjg(AL) + 
     &          (BL*(-2*g1L*g2R*muu1*muu2*muu3*(g3R*g4L*muu3*(-1 + muu3**2 - muu4**2) + g3L*g4R*muu4*(1 + muu3**2 - muu4**2)) + 
     &                g1R*g2L*(2*g3L*g4R*muu3*muu4*(-muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(1 + muu2**2 - muu4**2 + tred)) + 
     &                   g3R*g4L*(muu2**4*muu3**2 + muu1**4*muu4**2 + muu2**2*muu3**2*(-1 - 2*muu1**2 + muu3**2 - muu4**2) - 
     &                      (muu3**2 - muu4**2)*(muu3**2*(muu4**2 - tred) + tred*(1 - muu4**2 + tred)) + 
     &                      muu1**2*(muu4**4 + 2*muu3**2*tred - muu4**2*(1 + muu3**2 + 2*tred))))) + 
     &             BR*(-2*g1L*g2R*muu1*muu2*muu3*(g3L*g4R*muu3*(-1 + muu3**2 - muu4**2) + g3R*g4L*muu4*(1 + muu3**2 - muu4**2)) + 
     &                g1R*g2L*(2*g3R*g4L*muu3*muu4*(-muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(1 + muu2**2 - muu4**2 + tred)) + 
     &                   g3L*g4R*(muu2**4*muu3**2 + muu1**4*muu4**2 + muu2**2*muu3**2*(-1 - 2*muu1**2 + muu3**2 - muu4**2) - 
     &                      (muu3**2 - muu4**2)*(muu3**2*(muu4**2 - tred) + tred*(1 - muu4**2 + tred)) + 
     &                      muu1**2*(muu4**4 + 2*muu3**2*tred - muu4**2*(1 + muu3**2 + 2*tred))))))*dconjg(AR)) + 
     &       xd3*(-((BL*(g1R*muu1*(g2L*muu2*(g3R*(g4L*((-1 + muu2**2 + muu3**2 - tred)*(2*muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred)) + 
     &                           g4R*mGlr*muu4*(muu1**2 + muu3**2 - tred)) + 
     &                        g3L*muu3*(g4R*muu4*(3*muu1**2 + muu3**2 - tred) + g4L*mGlr*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                     g2R*M4Lr*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                        g3R*g4L*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))) + 
     &                  g1L*(g2L*M4Lr*muu2*(2*g3L*g4R*muu1**2*muu3*muu4 + 
     &                        g3R*g4L*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))) + 
     &                     g2R*muu1**2*(g3L*muu3*(g4L*mGlr*(muu2**2 + muu4**2 - tred) + g4R*muu4*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred)) + 
     &                        g3R*(g4R*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + 
     &                           g4L*(2*muu3**2*muu4**2 + muu1**2*(muu2**2 - tred) + muu2**2*(muu3**2 - tred) + tred - 2*muu3**2*tred - muu4**2*tred + 
     &                              tred**2))))) + BR*(g1R*muu1*
     &                   (g2L*muu2*(g3L*(g4R*((-1 + muu2**2 + muu3**2 - tred)*(2*muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred)) + 
     &                           g4L*mGlr*muu4*(muu1**2 + muu3**2 - tred)) + 
     &                        g3R*muu3*(g4L*muu4*(3*muu1**2 + muu3**2 - tred) + g4R*mGlr*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                     g2R*M4Lr*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                        g3L*g4R*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))) + 
     &                  g1L*(g2L*M4Lr*muu2*(2*g3R*g4L*muu1**2*muu3*muu4 + 
     &                        g3L*g4R*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))) + 
     &                     g2R*muu1**2*(g3R*muu3*(g4R*mGlr*(muu2**2 + muu4**2 - tred) + g4L*muu4*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred)) + 
     &                        g3L*(g4L*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + 
     &                           g4R*(2*muu3**2*muu4**2 + muu1**2*(muu2**2 - tred) + muu2**2*(muu3**2 - tred) + tred - 2*muu3**2*tred - muu4**2*tred + 
     &                              tred**2))))))*dconjg(AL)) - 
     &          (BL*(g1L*muu1*(g2R*muu2*(g3R*(g4L*((-1 + muu2**2 + muu3**2 - tred)*(2*muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred)) + 
     &                         g4R*mGlr*muu4*(muu1**2 + muu3**2 - tred)) + 
     &                      g3L*muu3*(g4R*muu4*(3*muu1**2 + muu3**2 - tred) + g4L*mGlr*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                   g2L*M4Lr*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                      g3R*g4L*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))) + 
     &                g1R*(g2R*M4Lr*muu2*(2*g3L*g4R*muu1**2*muu3*muu4 + 
     &                      g3R*g4L*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))) + 
     &                   g2L*muu1**2*(g3L*muu3*(g4L*mGlr*(muu2**2 + muu4**2 - tred) + g4R*muu4*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred)) + 
     &                      g3R*(g4R*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + 
     &                         g4L*(2*muu3**2*muu4**2 + muu1**2*(muu2**2 - tred) + muu2**2*(muu3**2 - tred) + tred - 2*muu3**2*tred - muu4**2*tred + 
     &                            tred**2))))) + BR*(g1L*muu1*(g2R*muu2*
     &                    (g3L*(g4R*((-1 + muu2**2 + muu3**2 - tred)*(2*muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred)) + 
     &                         g4L*mGlr*muu4*(muu1**2 + muu3**2 - tred)) + 
     &                      g3R*muu3*(g4L*muu4*(3*muu1**2 + muu3**2 - tred) + g4R*mGlr*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                   g2L*M4Lr*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                      g3L*g4R*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))) + 
     &                g1R*(g2R*M4Lr*muu2*(2*g3R*g4L*muu1**2*muu3*muu4 + 
     &                      g3L*g4R*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))) + 
     &                   g2L*muu1**2*(g3R*muu3*(g4R*mGlr*(muu2**2 + muu4**2 - tred) + g4L*muu4*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred)) + 
     &                      g3L*(g4L*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + 
     &                         g4R*(2*muu3**2*muu4**2 + muu1**2*(muu2**2 - tred) + muu2**2*(muu3**2 - tred) + tred - 2*muu3**2*tred - muu4**2*tred + 
     &                            tred**2))))))*dconjg(AR)) + xd1*
     &        (-((BL*(g1R*muu1*muu3*(g2L*muu2*(g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2) + g4R*muu4*(muu1**2 + 3*muu3**2 - tred)) + 
     &                        g3R*muu3*(2*g4R*mGlr*muu4 + g4L*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred))) + 
     &                     g2R*M4Lr*(g3L*g4R*muu4*(-1 + muu1**2 + muu4**2 - tred) + g3R*g4L*muu3*(muu2**2 + muu4**2 - tred))) + 
     &                  g1L*(g2L*M4Lr*muu2*muu3*(g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + g3R*g4L*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                     g2R*(g3R*(g4R*mGlr*muu4*(muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu4**2 - 2*tred) + tred - muu3**2*tred - 
     &                              muu4**2*tred + tred**2) + g4L*muu3**2*
     &                            (muu1**2*(2*muu2**2 + muu4**2 - 2*tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))) + 
     &                        g3L*muu3*(g4R*muu4*(2*muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-2 + 2*muu4**2 - 3*tred) + tred - 
     &                              muu3**2*tred - muu4**2*tred + tred**2) + 
     &                           g4L*mGlr*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))))) + 
     &               BR*(g1R*muu1*muu3*(g2L*muu2*(g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2) + g4L*muu4*(muu1**2 + 3*muu3**2 - tred)) + 
     &                        g3L*muu3*(2*g4L*mGlr*muu4 + g4R*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred))) + 
     &                     g2R*M4Lr*(g3R*g4L*muu4*(-1 + muu1**2 + muu4**2 - tred) + g3L*g4R*muu3*(muu2**2 + muu4**2 - tred))) + 
     &                  g1L*(g2L*M4Lr*muu2*muu3*(g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + g3L*g4R*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                     g2R*(g3L*(g4L*mGlr*muu4*(muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu4**2 - 2*tred) + tred - muu3**2*tred - 
     &                              muu4**2*tred + tred**2) + g4R*muu3**2*
     &                            (muu1**2*(2*muu2**2 + muu4**2 - 2*tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))) + 
     &                        g3R*muu3*(g4L*muu4*(2*muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-2 + 2*muu4**2 - 3*tred) + tred - 
     &                              muu3**2*tred - muu4**2*tred + tred**2) + 
     &                           g4R*mGlr*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))))))*dconjg(AL)) - 
     &          (BL*(g1L*muu1*muu3*(g2R*muu2*(g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2) + g4R*muu4*(muu1**2 + 3*muu3**2 - tred)) + 
     &                      g3R*muu3*(2*g4R*mGlr*muu4 + g4L*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred))) + 
     &                   g2L*M4Lr*(g3L*g4R*muu4*(-1 + muu1**2 + muu4**2 - tred) + g3R*g4L*muu3*(muu2**2 + muu4**2 - tred))) + 
     &                g1R*(g2R*M4Lr*muu2*muu3*(g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + g3R*g4L*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                   g2L*(g3R*(g4R*mGlr*muu4*(muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu4**2 - 2*tred) + tred - muu3**2*tred - 
     &                            muu4**2*tred + tred**2) + g4L*muu3**2*
     &                          (muu1**2*(2*muu2**2 + muu4**2 - 2*tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))) + 
     &                      g3L*muu3*(g4R*muu4*(2*muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-2 + 2*muu4**2 - 3*tred) + tred - muu3**2*tred - 
     &                            muu4**2*tred + tred**2) + g4L*mGlr*
     &                          (muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))))) + 
     &             BR*(g1L*muu1*muu3*(g2R*muu2*(g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2) + g4L*muu4*(muu1**2 + 3*muu3**2 - tred)) + 
     &                      g3L*muu3*(2*g4L*mGlr*muu4 + g4R*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred))) + 
     &                   g2L*M4Lr*(g3R*g4L*muu4*(-1 + muu1**2 + muu4**2 - tred) + g3L*g4R*muu3*(muu2**2 + muu4**2 - tred))) + 
     &                g1R*(g2R*M4Lr*muu2*muu3*(g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + g3L*g4R*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                   g2L*(g3L*(g4L*mGlr*muu4*(muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu4**2 - 2*tred) + tred - muu3**2*tred - 
     &                            muu4**2*tred + tred**2) + g4R*muu3**2*
     &                          (muu1**2*(2*muu2**2 + muu4**2 - 2*tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))) + 
     &                      g3R*muu3*(g4L*muu4*(2*muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-2 + 2*muu4**2 - 3*tred) + tred - muu3**2*tred - 
     &                            muu4**2*tred + tred**2) + g4R*mGlr*
     &                          (muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))))))*dconjg(AR)) + 
     &       xd2*((BL*(g1R*muu1*(g2R*M4Lr*(g3L*g4R*(1 - muu1**2 + muu2**2)*muu3*muu4 - 
     &                      g3R*g4L*(muu2**2*muu3**2 - muu4**4 - muu3**2*tred + muu4**2*(1 - muu1**2 + muu3**2 + tred))) + 
     &                   g2L*muu2*(g3L*muu3*(g4L*mGlr*(1 - muu3**2 + muu4**2) + g4R*muu4*(-2 - muu1**2 + muu2**2 - muu3**2 + muu4**2)) + 
     &                      g3R*(g4R*mGlr*muu4*(-1 - muu3**2 + muu4**2) + 
     &                         g4L*(-2*muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(2 - muu2**2 + 2*muu4**2 + tred))))) + 
     &                g1L*(g2L*M4Lr*muu2*(g3L*g4R*(-1 - muu1**2 + muu2**2)*muu3*muu4 + 
     &                      g3R*g4L*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 - muu2**2 + muu4**2 + tred))) + 
     &                   g2R*(g3L*muu3*(g4L*mGlr*(muu2**4 - muu2**2*(1 + muu1**2 - muu3**2 + tred) + muu1**2*(-muu4**2 + tred)) + 
     &                         g4R*muu4*(-2*muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(2 + 2*muu2**2 - muu4**2 + tred))) + 
     &                      g3R*(g4R*mGlr*muu4*(-muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(1 + muu2**2 - muu4**2 + tred)) + 
     &                         g4L*(muu2**4*muu3**2 + muu2**2*muu3**2*(-1 - 2*muu1**2 + muu3**2 - tred) + 
     &                            muu1**2*(muu4**4 + muu4**2*(-1 + muu1**2 - 2*muu3**2 - tred) + 2*muu3**2*tred)))))) + 
     &             BR*(g1R*muu1*(g2R*M4Lr*(g3R*g4L*(1 - muu1**2 + muu2**2)*muu3*muu4 - 
     &                      g3L*g4R*(muu2**2*muu3**2 - muu4**4 - muu3**2*tred + muu4**2*(1 - muu1**2 + muu3**2 + tred))) + 
     &                   g2L*muu2*(g3R*muu3*(g4R*mGlr*(1 - muu3**2 + muu4**2) + g4L*muu4*(-2 - muu1**2 + muu2**2 - muu3**2 + muu4**2)) + 
     &                      g3L*(g4L*mGlr*muu4*(-1 - muu3**2 + muu4**2) + 
     &                         g4R*(-2*muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(2 - muu2**2 + 2*muu4**2 + tred))))) + 
     &                g1L*(g2L*M4Lr*muu2*(g3R*g4L*(-1 - muu1**2 + muu2**2)*muu3*muu4 + 
     &                      g3L*g4R*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 - muu2**2 + muu4**2 + tred))) + 
     &                   g2R*(g3R*muu3*(g4R*mGlr*(muu2**4 - muu2**2*(1 + muu1**2 - muu3**2 + tred) + muu1**2*(-muu4**2 + tred)) + 
     &                         g4L*muu4*(-2*muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(2 + 2*muu2**2 - muu4**2 + tred))) + 
     &                      g3L*(g4L*mGlr*muu4*(-muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(1 + muu2**2 - muu4**2 + tred)) + 
     &                         g4R*(muu2**4*muu3**2 + muu2**2*muu3**2*(-1 - 2*muu1**2 + muu3**2 - tred) + 
     &                            muu1**2*(muu4**4 + muu4**2*(-1 + muu1**2 - 2*muu3**2 - tred) + 2*muu3**2*tred)))))))*dconjg(AL) + 
     &          (BL*(g1L*muu1*(g2L*M4Lr*(g3L*g4R*(1 - muu1**2 + muu2**2)*muu3*muu4 - 
     &                      g3R*g4L*(muu2**2*muu3**2 - muu4**4 - muu3**2*tred + muu4**2*(1 - muu1**2 + muu3**2 + tred))) + 
     &                   g2R*muu2*(g3L*muu3*(g4L*mGlr*(1 - muu3**2 + muu4**2) + g4R*muu4*(-2 - muu1**2 + muu2**2 - muu3**2 + muu4**2)) + 
     &                      g3R*(g4R*mGlr*muu4*(-1 - muu3**2 + muu4**2) + 
     &                         g4L*(-2*muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(2 - muu2**2 + 2*muu4**2 + tred))))) + 
     &                g1R*(g2R*M4Lr*muu2*(g3L*g4R*(-1 - muu1**2 + muu2**2)*muu3*muu4 + 
     &                      g3R*g4L*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 - muu2**2 + muu4**2 + tred))) + 
     &                   g2L*(g3L*muu3*(g4L*mGlr*(muu2**4 - muu2**2*(1 + muu1**2 - muu3**2 + tred) + muu1**2*(-muu4**2 + tred)) + 
     &                         g4R*muu4*(-2*muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(2 + 2*muu2**2 - muu4**2 + tred))) + 
     &                      g3R*(g4R*mGlr*muu4*(-muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(1 + muu2**2 - muu4**2 + tred)) + 
     &                         g4L*(muu2**4*muu3**2 + muu2**2*muu3**2*(-1 - 2*muu1**2 + muu3**2 - tred) + 
     &                            muu1**2*(muu4**4 + muu4**2*(-1 + muu1**2 - 2*muu3**2 - tred) + 2*muu3**2*tred)))))) + 
     &             BR*(g1L*muu1*(g2L*M4Lr*(g3R*g4L*(1 - muu1**2 + muu2**2)*muu3*muu4 - 
     &                      g3L*g4R*(muu2**2*muu3**2 - muu4**4 - muu3**2*tred + muu4**2*(1 - muu1**2 + muu3**2 + tred))) + 
     &                   g2R*muu2*(g3R*muu3*(g4R*mGlr*(1 - muu3**2 + muu4**2) + g4L*muu4*(-2 - muu1**2 + muu2**2 - muu3**2 + muu4**2)) + 
     &                      g3L*(g4L*mGlr*muu4*(-1 - muu3**2 + muu4**2) + 
     &                         g4R*(-2*muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(2 - muu2**2 + 2*muu4**2 + tred))))) + 
     &                g1R*(g2R*M4Lr*muu2*(g3R*g4L*(-1 - muu1**2 + muu2**2)*muu3*muu4 + 
     &                      g3L*g4R*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 - muu2**2 + muu4**2 + tred))) + 
     &                   g2L*(g3R*muu3*(g4R*mGlr*(muu2**4 - muu2**2*(1 + muu1**2 - muu3**2 + tred) + muu1**2*(-muu4**2 + tred)) + 
     &                         g4L*muu4*(-2*muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(2 + 2*muu2**2 - muu4**2 + tred))) + 
     &                      g3L*(g4L*mGlr*muu4*(-muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(1 + muu2**2 - muu4**2 + tred)) + 
     &                         g4R*(muu2**4*muu3**2 + muu2**2*muu3**2*(-1 - 2*muu1**2 + muu3**2 - tred) + 
     &                            muu1**2*(muu4**4 + muu4**2*(-1 + muu1**2 - 2*muu3**2 - tred) + 2*muu3**2*tred)))))))*dconjg(AR))))/(4.*DenSmV1*Pi**2)
	

	tempBoxZ = tempBoxZ + ChiChi2QQ_M2BoxGluino - tempBox2
	tempBox2 = ChiChi2QQ_M2BoxGluino

! Hboson contr.


	do i = imin,imax
	
	ChiChi2QQ_M2BoxGluino = ChiChi2QQ_M2BoxGluino -(CF*NC*sqrtS**6*((g1L*(g2L*M4Lr*(-1 + muu1**2 + muu2**2)*muu3*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &             g2R*muu2*(g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2)*(muu1**2 + muu3**2 - tred) + 2*g4R*muu3**2*muu4*(3*muu1**2 + muu3**2 - tred)) + 
     &                g3R*muu3*(2*g4R*mGlr*muu4*(muu1**2 + muu3**2 - tred) + g4L*(-1 + muu3**2 + muu4**2)*(3*muu1**2 + muu3**2 - tred)))) + 
     &          g1R*muu1*(2*g2R*M4Lr*muu2*muu3*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &             g2L*(g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2)*(-1 + muu1**2 + muu4**2 - tred) + 
     &                   2*g4R*muu3**2*muu4*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred)) + 
     &                g3R*muu3*(2*g4R*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4L*(-1 + muu3**2 + muu4**2)*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred))
     &                )))*xd1*dconjg(CL(i))*dconjg(DL(i)) + (g1R*g2L*muu1*
     &           (g3R*g4L*muu3*(2 - 2*muu3**2 + 2*muu1**2*(-1 + muu3**2) - muu4**2 + muu3**2*muu4**2 - muu4**4 - muu2**2*(-1 + muu3**2 + muu4**2) + tred - 
     &                muu3**2*tred + muu4**2*tred) - g3L*g4R*muu4*
     &              (muu1**2*(-1 - 3*muu3**2 + muu4**2) + (-1 + muu4**2)*(-1 + muu4**2 - tred) + muu3**2*(3 + 2*muu2**2 - muu4**2 + tred))) + 
     &          g1L*g2R*muu2*(g3R*g4L*muu3*(-1 + muu3**4 + 2*muu1**2*(-1 + muu3**2) + muu4**2 - muu3**2*muu4**2 - muu2**2*(-1 + muu3**2 + muu4**2) + 
     &                tred - muu3**2*tred + muu4**2*tred) + g3L*g4R*muu4*
     &              (muu3**4 + muu1**2*(1 + 3*muu3**2 - muu4**2) + (-1 + muu4**2)*tred - muu3**2*(-3 + 2*muu2**2 + muu4**2 + tred))))*xd12*dconjg(CL(i))*
     &        dconjg(DL(i)) + (g1L*g2R*muu2*(g3L*g4R*muu4*(muu1**4 + muu1**2*(6*muu3**2 - 2*tred) + (muu3**2 - tred)**2) + 
     &             g3R*g4L*muu3*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(-3 + muu2**2 + 3*muu3**2 + 2*muu4**2 - tred))) + 
     &          g1R*g2L*muu1*(g3R*g4L*muu3*(2 - 2*muu3**2 - 2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred) + 
     &                muu2**2*(-2 + muu3**2 + 2*muu4**2 - tred) + 2*tred - muu3**2*tred - muu4**2*tred + tred**2) + 
     &             g3L*g4R*muu4*(muu1**4 + muu1**2*(-1 + 3*muu3**2 + muu4**2 - 2*tred) + muu3**2*(-3 + 2*muu2**2 + muu4**2 - tred) + 
     &                tred*(1 - muu4**2 + tred))))*xd13*dconjg(CL(i))*dconjg(DL(i)) + 
     &       (g1R*muu1*(2*g2R*M4Lr*muu2*(g3R*g4L*muu3*(-1 + muu3**2 - muu4**2) + g3L*g4R*muu4*(1 + muu3**2 - muu4**2)) + 
     &             g2L*(2*g3R*muu3*(g4R*mGlr*(-1 + muu1**2 - muu2**2)*muu4 + g4L*(1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu2**2*muu4**2)) + 
     &                g3L*(g4L*mGlr*(-1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   g4R*muu4*(-1 - 3*muu3**2 + muu4**2 + muu1**2*(1 + 3*muu3**2 - muu4**2) - muu2**2*(-1 + muu3**2 + muu4**2))))) + 
     &          g1L*(g2L*M4Lr*(-1 + muu1**2 + muu2**2)*(g3R*g4L*muu3*(-1 + muu3**2 - muu4**2) + g3L*g4R*muu4*(1 + muu3**2 - muu4**2)) + 
     &             g2R*muu2*(g3L*(g4L*mGlr*(1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   2*g4R*muu4*(-((-1 + muu2**2)*muu3**2) + muu1**2*(1 + 2*muu3**2 - muu4**2))) + 
     &                g3R*muu3*(2*g4R*mGlr*(1 + muu1**2 - muu2**2)*muu4 + 
     &                   g4L*(muu1**2*(-3 + 3*muu3**2 - muu4**2) - (-1 + muu2**2)*(-1 + muu3**2 + muu4**2))))))*xd2*dconjg(CL(i))*dconjg(DL(i)) + 
     &       (g1R*g2L*muu1*(g3R*g4L*muu3*(2 - muu2**4 - 2*muu3**2 + muu4**2 + muu1**2*(-2 + muu2**2 + 2*muu3**2 - muu4**2 - tred) + tred + 
     &                muu2**2*(-1 - muu4**2 + tred)) + g3L*g4R*muu4*
     &              (-1 + muu1**4 - 2*muu3**2 + muu4**2 + tred + muu2**2*(1 - muu4**2 + tred) - muu1**2*(muu2**2 - 2*muu3**2 + muu4**2 + tred))) - 
     &          g1L*g2R*muu2*(g3R*g4L*muu3*((-1 + muu2**2)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(3 - muu2**2 - 3*muu3**2 + 2*muu4**2 + tred)) + 
     &             g3L*g4R*muu4*(-muu1**4 + (-1 + muu2**2)*(muu3**2 - tred) + muu1**2*(-3 + muu2**2 - 3*muu3**2 + 2*muu4**2 + tred))))*xd23*dconjg(CL(i))*
     &        dconjg(DL(i)) + (g1R*muu1*(g2L*(-1 + muu1**2 + muu2**2)*
     &              (g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2) + g4R*muu4*(muu1**2 + 3*muu3**2 - tred)) + 
     &                g3R*muu3*(2*g4R*mGlr*muu4 + g4L*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred))) + 
     &             2*g2R*M4Lr*muu2*(g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + g3R*g4L*muu3*(-1 + muu2**2 + muu3**2 - tred))) + 
     &          g1L*(2*g2R*muu1**2*muu2*(g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2) + g4R*muu4*(muu1**2 + 3*muu3**2 - tred)) + 
     &                g3R*muu3*(2*g4R*mGlr*muu4 + g4L*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred))) + 
     &             g2L*M4Lr*(-1 + muu1**2 + muu2**2)*(g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + g3R*g4L*muu3*(-1 + muu2**2 + muu3**2 - tred))))*xd3*
     &        dconjg(CL(i))*dconjg(DL(i)) + (g1R*(g2R*M4Lr*(-1 + muu1**2 + muu2**2)*muu3*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &             g2L*muu2*(g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2)*(muu1**2 + muu3**2 - tred) + 2*g4R*muu3**2*muu4*(3*muu1**2 + muu3**2 - tred)) + 
     &                g3R*muu3*(2*g4R*mGlr*muu4*(muu1**2 + muu3**2 - tred) + g4L*(-1 + muu3**2 + muu4**2)*(3*muu1**2 + muu3**2 - tred)))) + 
     &          g1L*muu1*(2*g2L*M4Lr*muu2*muu3*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &             g2R*(g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2)*(-1 + muu1**2 + muu4**2 - tred) + 
     &                   2*g4R*muu3**2*muu4*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred)) + 
     &                g3R*muu3*(2*g4R*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4L*(-1 + muu3**2 + muu4**2)*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred))
     &                )))*xd1*dconjg(CR(i))*dconjg(DL(i)) + (g1L*g2R*muu1*
     &           (g3R*g4L*muu3*(2 - 2*muu3**2 + 2*muu1**2*(-1 + muu3**2) - muu4**2 + muu3**2*muu4**2 - muu4**4 - muu2**2*(-1 + muu3**2 + muu4**2) + tred - 
     &                muu3**2*tred + muu4**2*tred) - g3L*g4R*muu4*
     &              (muu1**2*(-1 - 3*muu3**2 + muu4**2) + (-1 + muu4**2)*(-1 + muu4**2 - tred) + muu3**2*(3 + 2*muu2**2 - muu4**2 + tred))) + 
     &          g1R*g2L*muu2*(g3R*g4L*muu3*(-1 + muu3**4 + 2*muu1**2*(-1 + muu3**2) + muu4**2 - muu3**2*muu4**2 - muu2**2*(-1 + muu3**2 + muu4**2) + 
     &                tred - muu3**2*tred + muu4**2*tred) + g3L*g4R*muu4*
     &              (muu3**4 + muu1**2*(1 + 3*muu3**2 - muu4**2) + (-1 + muu4**2)*tred - muu3**2*(-3 + 2*muu2**2 + muu4**2 + tred))))*xd12*dconjg(CR(i))*
     &        dconjg(DL(i)) + (g1R*g2L*muu2*(g3L*g4R*muu4*(muu1**4 + muu1**2*(6*muu3**2 - 2*tred) + (muu3**2 - tred)**2) + 
     &             g3R*g4L*muu3*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(-3 + muu2**2 + 3*muu3**2 + 2*muu4**2 - tred))) + 
     &          g1L*g2R*muu1*(g3R*g4L*muu3*(2 - 2*muu3**2 - 2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred) + 
     &                muu2**2*(-2 + muu3**2 + 2*muu4**2 - tred) + 2*tred - muu3**2*tred - muu4**2*tred + tred**2) + 
     &             g3L*g4R*muu4*(muu1**4 + muu1**2*(-1 + 3*muu3**2 + muu4**2 - 2*tred) + muu3**2*(-3 + 2*muu2**2 + muu4**2 - tred) + 
     &                tred*(1 - muu4**2 + tred))))*xd13*dconjg(CR(i))*dconjg(DL(i)) + 
     &       (g1L*muu1*(2*g2L*M4Lr*muu2*(g3R*g4L*muu3*(-1 + muu3**2 - muu4**2) + g3L*g4R*muu4*(1 + muu3**2 - muu4**2)) + 
     &             g2R*(2*g3R*muu3*(g4R*mGlr*(-1 + muu1**2 - muu2**2)*muu4 + g4L*(1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu2**2*muu4**2)) + 
     &                g3L*(g4L*mGlr*(-1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   g4R*muu4*(-1 - 3*muu3**2 + muu4**2 + muu1**2*(1 + 3*muu3**2 - muu4**2) - muu2**2*(-1 + muu3**2 + muu4**2))))) + 
     &          g1R*(g2R*M4Lr*(-1 + muu1**2 + muu2**2)*(g3R*g4L*muu3*(-1 + muu3**2 - muu4**2) + g3L*g4R*muu4*(1 + muu3**2 - muu4**2)) + 
     &             g2L*muu2*(g3L*(g4L*mGlr*(1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   2*g4R*muu4*(-((-1 + muu2**2)*muu3**2) + muu1**2*(1 + 2*muu3**2 - muu4**2))) + 
     &                g3R*muu3*(2*g4R*mGlr*(1 + muu1**2 - muu2**2)*muu4 + 
     &                   g4L*(muu1**2*(-3 + 3*muu3**2 - muu4**2) - (-1 + muu2**2)*(-1 + muu3**2 + muu4**2))))))*xd2*dconjg(CR(i))*dconjg(DL(i)) + 
     &       (g1L*g2R*muu1*(g3R*g4L*muu3*(2 - muu2**4 - 2*muu3**2 + muu4**2 + muu1**2*(-2 + muu2**2 + 2*muu3**2 - muu4**2 - tred) + tred + 
     &                muu2**2*(-1 - muu4**2 + tred)) + g3L*g4R*muu4*
     &              (-1 + muu1**4 - 2*muu3**2 + muu4**2 + tred + muu2**2*(1 - muu4**2 + tred) - muu1**2*(muu2**2 - 2*muu3**2 + muu4**2 + tred))) - 
     &          g1R*g2L*muu2*(g3R*g4L*muu3*((-1 + muu2**2)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(3 - muu2**2 - 3*muu3**2 + 2*muu4**2 + tred)) + 
     &             g3L*g4R*muu4*(-muu1**4 + (-1 + muu2**2)*(muu3**2 - tred) + muu1**2*(-3 + muu2**2 - 3*muu3**2 + 2*muu4**2 + tred))))*xd23*dconjg(CR(i))*
     &        dconjg(DL(i)) + (g1L*muu1*(g2R*(-1 + muu1**2 + muu2**2)*
     &              (g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2) + g4R*muu4*(muu1**2 + 3*muu3**2 - tred)) + 
     &                g3R*muu3*(2*g4R*mGlr*muu4 + g4L*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred))) + 
     &             2*g2L*M4Lr*muu2*(g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + g3R*g4L*muu3*(-1 + muu2**2 + muu3**2 - tred))) + 
     &          g1R*(2*g2L*muu1**2*muu2*(g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2) + g4R*muu4*(muu1**2 + 3*muu3**2 - tred)) + 
     &                g3R*muu3*(2*g4R*mGlr*muu4 + g4L*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred))) + 
     &             g2R*M4Lr*(-1 + muu1**2 + muu2**2)*(g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + g3R*g4L*muu3*(-1 + muu2**2 + muu3**2 - tred))))*xd3*
     &        dconjg(CR(i))*dconjg(DL(i)) + (g1L*(g2L*M4Lr*(-1 + muu1**2 + muu2**2)*muu3*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &             g2R*muu2*(g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2)*(muu1**2 + muu3**2 - tred) + 2*g4L*muu3**2*muu4*(3*muu1**2 + muu3**2 - tred)) + 
     &                g3L*muu3*(2*g4L*mGlr*muu4*(muu1**2 + muu3**2 - tred) + g4R*(-1 + muu3**2 + muu4**2)*(3*muu1**2 + muu3**2 - tred)))) + 
     &          g1R*muu1*(2*g2R*M4Lr*muu2*muu3*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &             g2L*(g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2)*(-1 + muu1**2 + muu4**2 - tred) + 
     &                   2*g4L*muu3**2*muu4*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred)) + 
     &                g3L*muu3*(2*g4L*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4R*(-1 + muu3**2 + muu4**2)*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred))
     &                )))*xd1*dconjg(CL(i))*dconjg(DR(i)) + (g1R*g2L*muu1*
     &           (g3L*g4R*muu3*(2 - 2*muu3**2 + 2*muu1**2*(-1 + muu3**2) - muu4**2 + muu3**2*muu4**2 - muu4**4 - muu2**2*(-1 + muu3**2 + muu4**2) + tred - 
     &                muu3**2*tred + muu4**2*tred) - g3R*g4L*muu4*
     &              (muu1**2*(-1 - 3*muu3**2 + muu4**2) + (-1 + muu4**2)*(-1 + muu4**2 - tred) + muu3**2*(3 + 2*muu2**2 - muu4**2 + tred))) + 
     &          g1L*g2R*muu2*(g3L*g4R*muu3*(-1 + muu3**4 + 2*muu1**2*(-1 + muu3**2) + muu4**2 - muu3**2*muu4**2 - muu2**2*(-1 + muu3**2 + muu4**2) + 
     &                tred - muu3**2*tred + muu4**2*tred) + g3R*g4L*muu4*
     &              (muu3**4 + muu1**2*(1 + 3*muu3**2 - muu4**2) + (-1 + muu4**2)*tred - muu3**2*(-3 + 2*muu2**2 + muu4**2 + tred))))*xd12*dconjg(CL(i))*
     &        dconjg(DR(i)) + (g1L*g2R*muu2*(g3R*g4L*muu4*(muu1**4 + muu1**2*(6*muu3**2 - 2*tred) + (muu3**2 - tred)**2) + 
     &             g3L*g4R*muu3*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(-3 + muu2**2 + 3*muu3**2 + 2*muu4**2 - tred))) + 
     &          g1R*g2L*muu1*(g3L*g4R*muu3*(2 - 2*muu3**2 - 2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred) + 
     &                muu2**2*(-2 + muu3**2 + 2*muu4**2 - tred) + 2*tred - muu3**2*tred - muu4**2*tred + tred**2) + 
     &             g3R*g4L*muu4*(muu1**4 + muu1**2*(-1 + 3*muu3**2 + muu4**2 - 2*tred) + muu3**2*(-3 + 2*muu2**2 + muu4**2 - tred) + 
     &                tred*(1 - muu4**2 + tred))))*xd13*dconjg(CL(i))*dconjg(DR(i)) + 
     &       (g1R*muu1*(2*g2R*M4Lr*muu2*(g3L*g4R*muu3*(-1 + muu3**2 - muu4**2) + g3R*g4L*muu4*(1 + muu3**2 - muu4**2)) + 
     &             g2L*(2*g3L*muu3*(g4L*mGlr*(-1 + muu1**2 - muu2**2)*muu4 + g4R*(1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu2**2*muu4**2)) + 
     &                g3R*(g4R*mGlr*(-1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   g4L*muu4*(-1 - 3*muu3**2 + muu4**2 + muu1**2*(1 + 3*muu3**2 - muu4**2) - muu2**2*(-1 + muu3**2 + muu4**2))))) + 
     &          g1L*(g2L*M4Lr*(-1 + muu1**2 + muu2**2)*(g3L*g4R*muu3*(-1 + muu3**2 - muu4**2) + g3R*g4L*muu4*(1 + muu3**2 - muu4**2)) + 
     &             g2R*muu2*(g3R*(g4R*mGlr*(1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   2*g4L*muu4*(-((-1 + muu2**2)*muu3**2) + muu1**2*(1 + 2*muu3**2 - muu4**2))) + 
     &                g3L*muu3*(2*g4L*mGlr*(1 + muu1**2 - muu2**2)*muu4 + 
     &                   g4R*(muu1**2*(-3 + 3*muu3**2 - muu4**2) - (-1 + muu2**2)*(-1 + muu3**2 + muu4**2))))))*xd2*dconjg(CL(i))*dconjg(DR(i)) + 
     &       (g1R*g2L*muu1*(g3L*g4R*muu3*(2 - muu2**4 - 2*muu3**2 + muu4**2 + muu1**2*(-2 + muu2**2 + 2*muu3**2 - muu4**2 - tred) + tred + 
     &                muu2**2*(-1 - muu4**2 + tred)) + g3R*g4L*muu4*
     &              (-1 + muu1**4 - 2*muu3**2 + muu4**2 + tred + muu2**2*(1 - muu4**2 + tred) - muu1**2*(muu2**2 - 2*muu3**2 + muu4**2 + tred))) - 
     &          g1L*g2R*muu2*(g3L*g4R*muu3*((-1 + muu2**2)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(3 - muu2**2 - 3*muu3**2 + 2*muu4**2 + tred)) + 
     &             g3R*g4L*muu4*(-muu1**4 + (-1 + muu2**2)*(muu3**2 - tred) + muu1**2*(-3 + muu2**2 - 3*muu3**2 + 2*muu4**2 + tred))))*xd23*dconjg(CL(i))*
     &        dconjg(DR(i)) + (g1R*muu1*(g2L*(-1 + muu1**2 + muu2**2)*
     &              (g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2) + g4L*muu4*(muu1**2 + 3*muu3**2 - tred)) + 
     &                g3L*muu3*(2*g4L*mGlr*muu4 + g4R*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred))) + 
     &             2*g2R*M4Lr*muu2*(g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + g3L*g4R*muu3*(-1 + muu2**2 + muu3**2 - tred))) + 
     &          g1L*(2*g2R*muu1**2*muu2*(g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2) + g4L*muu4*(muu1**2 + 3*muu3**2 - tred)) + 
     &                g3L*muu3*(2*g4L*mGlr*muu4 + g4R*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred))) + 
     &             g2L*M4Lr*(-1 + muu1**2 + muu2**2)*(g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + g3L*g4R*muu3*(-1 + muu2**2 + muu3**2 - tred))))*xd3*
     &        dconjg(CL(i))*dconjg(DR(i)) + (g1R*(g2R*M4Lr*(-1 + muu1**2 + muu2**2)*muu3*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &             g2L*muu2*(g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2)*(muu1**2 + muu3**2 - tred) + 2*g4L*muu3**2*muu4*(3*muu1**2 + muu3**2 - tred)) + 
     &                g3L*muu3*(2*g4L*mGlr*muu4*(muu1**2 + muu3**2 - tred) + g4R*(-1 + muu3**2 + muu4**2)*(3*muu1**2 + muu3**2 - tred)))) + 
     &          g1L*muu1*(2*g2L*M4Lr*muu2*muu3*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &             g2R*(g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2)*(-1 + muu1**2 + muu4**2 - tred) + 
     &                   2*g4L*muu3**2*muu4*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred)) + 
     &                g3L*muu3*(2*g4L*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4R*(-1 + muu3**2 + muu4**2)*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred))
     &                )))*xd1*dconjg(CR(i))*dconjg(DR(i)) + (g1L*g2R*muu1*
     &           (g3L*g4R*muu3*(2 - 2*muu3**2 + 2*muu1**2*(-1 + muu3**2) - muu4**2 + muu3**2*muu4**2 - muu4**4 - muu2**2*(-1 + muu3**2 + muu4**2) + tred - 
     &                muu3**2*tred + muu4**2*tred) - g3R*g4L*muu4*
     &              (muu1**2*(-1 - 3*muu3**2 + muu4**2) + (-1 + muu4**2)*(-1 + muu4**2 - tred) + muu3**2*(3 + 2*muu2**2 - muu4**2 + tred))) + 
     &          g1R*g2L*muu2*(g3L*g4R*muu3*(-1 + muu3**4 + 2*muu1**2*(-1 + muu3**2) + muu4**2 - muu3**2*muu4**2 - muu2**2*(-1 + muu3**2 + muu4**2) + 
     &                tred - muu3**2*tred + muu4**2*tred) + g3R*g4L*muu4*
     &              (muu3**4 + muu1**2*(1 + 3*muu3**2 - muu4**2) + (-1 + muu4**2)*tred - muu3**2*(-3 + 2*muu2**2 + muu4**2 + tred))))*xd12*dconjg(CR(i))*
     &        dconjg(DR(i)) + (g1R*g2L*muu2*(g3R*g4L*muu4*(muu1**4 + muu1**2*(6*muu3**2 - 2*tred) + (muu3**2 - tred)**2) + 
     &             g3L*g4R*muu3*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(-3 + muu2**2 + 3*muu3**2 + 2*muu4**2 - tred))) + 
     &          g1L*g2R*muu1*(g3L*g4R*muu3*(2 - 2*muu3**2 - 2*muu4**2 + muu3**2*muu4**2 + muu1**2*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred) + 
     &                muu2**2*(-2 + muu3**2 + 2*muu4**2 - tred) + 2*tred - muu3**2*tred - muu4**2*tred + tred**2) + 
     &             g3R*g4L*muu4*(muu1**4 + muu1**2*(-1 + 3*muu3**2 + muu4**2 - 2*tred) + muu3**2*(-3 + 2*muu2**2 + muu4**2 - tred) + 
     &                tred*(1 - muu4**2 + tred))))*xd13*dconjg(CR(i))*dconjg(DR(i)) + 
     &       (g1L*muu1*(2*g2L*M4Lr*muu2*(g3L*g4R*muu3*(-1 + muu3**2 - muu4**2) + g3R*g4L*muu4*(1 + muu3**2 - muu4**2)) + 
     &             g2R*(2*g3L*muu3*(g4L*mGlr*(-1 + muu1**2 - muu2**2)*muu4 + g4R*(1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu2**2*muu4**2)) + 
     &                g3R*(g4R*mGlr*(-1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   g4L*muu4*(-1 - 3*muu3**2 + muu4**2 + muu1**2*(1 + 3*muu3**2 - muu4**2) - muu2**2*(-1 + muu3**2 + muu4**2))))) + 
     &          g1R*(g2R*M4Lr*(-1 + muu1**2 + muu2**2)*(g3L*g4R*muu3*(-1 + muu3**2 - muu4**2) + g3R*g4L*muu4*(1 + muu3**2 - muu4**2)) + 
     &             g2L*muu2*(g3R*(g4R*mGlr*(1 + muu1**2 - muu2**2)*(-1 + muu3**2 + muu4**2) + 
     &                   2*g4L*muu4*(-((-1 + muu2**2)*muu3**2) + muu1**2*(1 + 2*muu3**2 - muu4**2))) + 
     &                g3L*muu3*(2*g4L*mGlr*(1 + muu1**2 - muu2**2)*muu4 + 
     &                   g4R*(muu1**2*(-3 + 3*muu3**2 - muu4**2) - (-1 + muu2**2)*(-1 + muu3**2 + muu4**2))))))*xd2*dconjg(CR(i))*dconjg(DR(i)) + 
     &       (g1L*g2R*muu1*(g3L*g4R*muu3*(2 - muu2**4 - 2*muu3**2 + muu4**2 + muu1**2*(-2 + muu2**2 + 2*muu3**2 - muu4**2 - tred) + tred + 
     &                muu2**2*(-1 - muu4**2 + tred)) + g3R*g4L*muu4*
     &              (-1 + muu1**4 - 2*muu3**2 + muu4**2 + tred + muu2**2*(1 - muu4**2 + tred) - muu1**2*(muu2**2 - 2*muu3**2 + muu4**2 + tred))) - 
     &          g1R*g2L*muu2*(g3L*g4R*muu3*((-1 + muu2**2)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(3 - muu2**2 - 3*muu3**2 + 2*muu4**2 + tred)) + 
     &             g3R*g4L*muu4*(-muu1**4 + (-1 + muu2**2)*(muu3**2 - tred) + muu1**2*(-3 + muu2**2 - 3*muu3**2 + 2*muu4**2 + tred))))*xd23*dconjg(CR(i))*
     &        dconjg(DR(i)) + (g1L*muu1*(g2R*(-1 + muu1**2 + muu2**2)*
     &              (g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2) + g4L*muu4*(muu1**2 + 3*muu3**2 - tred)) + 
     &                g3L*muu3*(2*g4L*mGlr*muu4 + g4R*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred))) + 
     &             2*g2L*M4Lr*muu2*(g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + g3L*g4R*muu3*(-1 + muu2**2 + muu3**2 - tred))) + 
     &          g1R*(2*g2L*muu1**2*muu2*(g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2) + g4L*muu4*(muu1**2 + 3*muu3**2 - tred)) + 
     &                g3L*muu3*(2*g4L*mGlr*muu4 + g4R*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred))) + 
     &             g2R*M4Lr*(-1 + muu1**2 + muu2**2)*(g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + g3L*g4R*muu3*(-1 + muu2**2 + muu3**2 - tred))))*xd3*
     &        dconjg(CR(i))*dconjg(DR(i)) + xd22*((g1R*g2L*muu1*(-1 + muu1**2 - muu2**2) + g1L*g2R*muu2*(1 + muu1**2 - muu2**2))*dconjg(CL(i)) + 
     &          (g1L*g2R*muu1*(-1 + muu1**2 - muu2**2) + g1R*g2L*muu2*(1 + muu1**2 - muu2**2))*dconjg(CR(i)))*
     &        ((g3R*g4L*muu3*(-1 + muu3**2 - muu4**2) + g3L*g4R*muu4*(1 + muu3**2 - muu4**2))*dconjg(DL(i)) + 
     &          (g3L*g4R*muu3*(-1 + muu3**2 - muu4**2) + g3R*g4L*muu4*(1 + muu3**2 - muu4**2))*dconjg(DR(i))) + 
     &       muu3*xd11*((g1L*g2R*muu2*(muu1**2 + muu3**2 - tred) + g1R*g2L*muu1*(-1 + muu1**2 + muu4**2 - tred))*dconjg(CL(i)) + 
     &          (g1R*g2L*muu2*(muu1**2 + muu3**2 - tred) + g1L*g2R*muu1*(-1 + muu1**2 + muu4**2 - tred))*dconjg(CR(i)))*
     &        ((2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2))*dconjg(DL(i)) + 
     &          (2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2))*dconjg(DR(i))) + 
     &       xd0*((g1R*muu1*(2*g2R*M4Lr*muu2 + g2L*(-1 + muu1**2 + muu2**2)) + g1L*(2*g2R*muu1**2*muu2 + g2L*M4Lr*(-1 + muu1**2 + muu2**2)))*
     &           dconjg(CL(i)) + (g1L*muu1*(2*g2L*M4Lr*muu2 + g2R*(-1 + muu1**2 + muu2**2)) + 
     &             g1R*(2*g2L*muu1**2*muu2 + g2R*M4Lr*(-1 + muu1**2 + muu2**2)))*dconjg(CR(i)))*
     &        ((g3R*muu3*(2*g4R*mGlr*muu4 + g4L*(-1 + muu3**2 + muu4**2)) + g3L*(2*g4R*muu3**2*muu4 + g4L*mGlr*(-1 + muu3**2 + muu4**2)))*dconjg(DL(i)) + 
     &          (g3L*muu3*(2*g4L*mGlr*muu4 + g4R*(-1 + muu3**2 + muu4**2)) + g3R*(2*g4L*muu3**2*muu4 + g4R*mGlr*(-1 + muu3**2 + muu4**2)))*dconjg(DR(i)))
     &        + muu1*xd33*((2*g1L*g2R*muu1*muu2 + g1R*g2L*(-1 + muu1**2 + muu2**2))*dconjg(CL(i)) + 
     &          (2*g1R*g2L*muu1*muu2 + g1L*g2R*(-1 + muu1**2 + muu2**2))*dconjg(CR(i)))*
     &        ((g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + g3R*g4L*muu3*(-1 + muu2**2 + muu3**2 - tred))*dconjg(DL(i)) + 
     &          (g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + g3L*g4R*muu3*(-1 + muu2**2 + muu3**2 - tred))*dconjg(DR(i))) + 
     &       (2*xd00*(dconjg(CL(i))*((g1L*g2R*muu2*(g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + g3R*g4L*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                  g1R*g2L*muu1*(g3L*g4R*muu4*(-1 + muu1**2 + muu4**2 - tred) + g3R*g4L*muu3*(muu2**2 + muu4**2 - tred)))*dconjg(DL(i)) + 
     &               (g1L*g2R*muu2*(g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + g3L*g4R*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                  g1R*g2L*muu1*(g3R*g4L*muu4*(-1 + muu1**2 + muu4**2 - tred) + g3L*g4R*muu3*(muu2**2 + muu4**2 - tred)))*dconjg(DR(i))) + 
     &            dconjg(CR(i))*((g1R*g2L*muu2*(g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + g3R*g4L*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                  g1L*g2R*muu1*(g3L*g4R*muu4*(-1 + muu1**2 + muu4**2 - tred) + g3R*g4L*muu3*(muu2**2 + muu4**2 - tred)))*dconjg(DL(i)) + 
     &               (g1R*g2L*muu2*(g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + g3L*g4R*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                  g1L*g2R*muu1*(g3R*g4L*muu4*(-1 + muu1**2 + muu4**2 - tred) + g3L*g4R*muu3*(muu2**2 + muu4**2 - tred)))*dconjg(DR(i)))))/sqrtS**2))/
     &  (8.*Pi**2*DenSmH1(i))
	
	enddo

	tempBoxH = tempBoxH + ChiChi2QQ_M2BoxGluino - tempBox2
	tempBox2 = ChiChi2QQ_M2BoxGluino
	
! U-channel contr.

	do i = 1,2
		
		ChiChi2QQ_M2BoxGluino = ChiChi2QQ_M2BoxGluino + (CF*NC*sqrtS**4*((xd00*(g1R*g2L*(-(g3L*g4R*(-2*muu1*muu4*dconjg(GL(indi,i,tt,gen)) + (-1 + muu2**2 + muu3**2 - tred)*dconjg(GR(indi,i,tt,gen)))*
     &                 (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen)))) + 
     &              g3R*g4L*(muu1*dconjg(GL(indi,i,tt,gen))*(muu3*(muu2**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 
     &                    muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                 muu4*dconjg(GR(indi,i,tt,gen))*(-((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen))) + 
     &                    muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &           g1L*g2R*(-(g3R*g4L*((-1 + muu2**2 + muu3**2 - tred)*dconjg(GL(indi,i,tt,gen)) - 2*muu1*muu4*dconjg(GR(indi,i,tt,gen)))*
     &                 (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g3L*g4R*(muu4*dconjg(GL(indi,i,tt,gen))*(muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 
     &                    (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 muu1*dconjg(GR(indi,i,tt,gen))*(-(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen))) + 
     &                    muu3*(muu2**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))))))/sqrtS**2 + 
     &      (muu1*xd33*(g1R*g2L*(g3R*g4L*(dconjg(GL(indi,i,tt,gen))*
     &                  (muu3*(muu1**2*(1 - muu3**2 + muu4**2) - (-1 + muu2**2)*(-1 + muu2**2 + muu3**2 - tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                    muu2*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))*dconjg(HR(indj,i,tt,gen))) + 
     &                 muu1*muu4*dconjg(GR(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                    muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g3L*g4R*muu1*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                  (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*((muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))*
     &                     dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &           g1L*g2R*(g3L*g4R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                  (-(muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*(muu2*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))*
     &                     dconjg(HL(indj,i,tt,gen)) + muu3*(muu1**2*(1 - muu3**2 + muu4**2) - (-1 + muu2**2)*(-1 + muu2**2 + muu3**2 - tred))*
     &                     dconjg(HR(indj,i,tt,gen)))) + g3R*g4L*muu1*
     &               (muu1*muu4*dconjg(GR(indi,i,tt,gen))*(2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                    (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                    (muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))*dconjg(HR(indj,i,tt,gen)))))))/2. + 
     &      (muu3*xd11*(g1R*g2L*(g3R*g4L*(muu1*dconjg(GL(indi,i,tt,gen))*
     &                  (((1 + muu2**2)*muu3**2 - muu1**2*(-1 + muu3**2 + muu4**2) - (-1 + muu4**2)*(-1 + muu4**2 - tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                    muu2*muu3*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                 muu4*dconjg(GR(indi,i,tt,gen))*((muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu4**2 - 2*tred) + tred - 
     &                       muu3**2*tred - muu4**2*tred + tred**2)*dconjg(HL(indj,i,tt,gen)) - 
     &                    muu2*muu3*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g3L*g4R*muu3*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                  (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*((muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))*
     &                     dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &           g1L*g2R*(g3L*g4R*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                  (muu2*muu3*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                    ((1 + muu2**2)*muu3**2 - muu1**2*(-1 + muu3**2 + muu4**2) - (-1 + muu4**2)*(-1 + muu4**2 - tred))*dconjg(HR(indj,i,tt,gen))) + 
     &                 muu4*dconjg(GL(indi,i,tt,gen))*(-(muu2*muu3*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                    (muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu4**2 - 2*tred) + tred - muu3**2*tred - muu4**2*tred + tred**2)*
     &                     dconjg(HR(indj,i,tt,gen)))) + g3R*g4L*muu3*
     &               (muu1*muu4*dconjg(GR(indi,i,tt,gen))*(2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                    (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                    (muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))*dconjg(HR(indj,i,tt,gen)))))))/2. + 
     &      (xd23*(g1R*g2L*(g3L*g4R*(-(muu1*(1 + muu1**2 - muu2**2)*muu4*dconjg(GL(indi,i,tt,gen))*
     &                    ((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*dconjg(HR(indj,i,tt,gen)))) + 
     &                 dconjg(GR(indi,i,tt,gen))*((-(muu2**2*(muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred)) + muu1**4*(muu2**2 - muu4**2 - tred) + 
     &                       muu1**2*(-muu2**4 + muu4**2 + 2*muu3**2*muu4**2 - muu4**4 + muu2**2*(muu3**2 + muu4**2) + tred - 2*muu3**2*tred + tred**2))*
     &                     dconjg(HL(indj,i,tt,gen)) + muu2*(-1 - muu1**2 + muu2**2)*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) - 
     &              2*g3R*g4L*muu1*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                  ((1 - muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + muu2*(1 + muu3**2 - muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(muu3*(1 - muu2**2 - muu3**2 + muu1**2*(-1 + muu3**2 - muu4**2) + tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                    muu2*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 - muu2**2 + muu4**2 + tred))*dconjg(HR(indj,i,tt,gen))))) + 
     &           g1L*g2R*(-2*g3L*g4R*muu1*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                  (muu2*(1 + muu3**2 - muu4**2)*dconjg(HL(indj,i,tt,gen)) + (1 - muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*(muu2*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 - muu2**2 + muu4**2 + tred))*
     &                     dconjg(HL(indj,i,tt,gen)) + muu3*(1 - muu2**2 - muu3**2 + muu1**2*(-1 + muu3**2 - muu4**2) + tred)*dconjg(HR(indj,i,tt,gen))))
     &               + g3R*g4L*(-(muu1*(1 + muu1**2 - muu2**2)*muu4*dconjg(GR(indi,i,tt,gen))*
     &                    (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(muu2*(-1 - muu1**2 + muu2**2)*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                    (-(muu2**2*(muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred)) + muu1**4*(muu2**2 - muu4**2 - tred) + 
     &                       muu1**2*(-muu2**4 + muu4**2 + 2*muu3**2*muu4**2 - muu4**4 + muu2**2*(muu3**2 + muu4**2) + tred - 2*muu3**2*tred + tred**2))*
     &                     dconjg(HR(indj,i,tt,gen)))))))/2. + 
     &      (xd22*(g1R*g2L*(g3L*g4R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                  (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*((-(muu2**4*muu3**2) + muu1**2*(-muu4**4 - muu3**2*tred + muu4**2*(1 - muu1**2 + muu3**2 + tred)) + 
     &                       muu2**2*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 + muu1**2 + muu4**2 + tred)))*dconjg(HL(indj,i,tt,gen)) - 
     &                    muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) - 
     &              g3R*g4L*(muu1*dconjg(GL(indi,i,tt,gen))*(muu3*(1 - (1 + muu2**2)*muu3**2 + muu2**2*muu4**2 + muu1**2*(-1 + muu3**2 - muu4**2) + tred)*
     &                     dconjg(HL(indj,i,tt,gen)) + muu2*(-muu3**4 + muu4**2 - muu4**4 + muu3**2*(1 + 2*muu4**2))*dconjg(HR(indj,i,tt,gen))) + 
     &                 muu4*dconjg(GR(indi,i,tt,gen))*((-muu1**4 + muu2**2 - muu2**4 + muu1**2*(1 + 2*muu2**2))*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                    muu2*(-1 + muu4**2 + muu2**2*(1 + muu3**2 - muu4**2) + muu1**2*(-muu3**2 + muu4**2) - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &           g1L*g2R*(-(g3L*g4R*(muu4*dconjg(GL(indi,i,tt,gen))*
     &                    (-(muu2*(-1 + muu4**2 + muu2**2*(1 + muu3**2 - muu4**2) + muu1**2*(-muu3**2 + muu4**2) - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                      (-muu1**4 + muu2**2 - muu2**4 + muu1**2*(1 + 2*muu2**2))*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                   muu1*dconjg(GR(indi,i,tt,gen))*(muu2*(-muu3**4 + muu4**2 - muu4**4 + muu3**2*(1 + 2*muu4**2))*dconjg(HL(indj,i,tt,gen)) + 
     &                      muu3*(1 - (1 + muu2**2)*muu3**2 + muu2**2*muu4**2 + muu1**2*(-1 + muu3**2 - muu4**2) + tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &              g3R*g4L*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                  (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                    (-(muu2**4*muu3**2) + muu1**2*(-muu4**4 - muu3**2*tred + muu4**2*(1 - muu1**2 + muu3**2 + tred)) + 
     &                       muu2**2*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 + muu1**2 + muu4**2 + tred)))*dconjg(HR(indj,i,tt,gen)))))))/2. + 
     &      (xd13*(g1R*g2L*(2*g3R*g4L*muu1*muu3*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                  ((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*
     &                     dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g3L*g4R*(-(muu1*muu4*(muu1**2 + muu3**2 - tred)*dconjg(GL(indi,i,tt,gen))*
     &                    ((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*dconjg(HR(indj,i,tt,gen)))) + 
     &                 dconjg(GR(indi,i,tt,gen))*((muu1**4*(muu2**2 - muu4**2 - tred) - 
     &                       (-1 + muu2**2 + muu3**2 - tred)*(muu2**2*muu3**2 + (-1 + muu4**2 - tred)*tred + muu3**2*(-muu4**2 + tred)) + 
     &                       muu1**2*((1 + 2*muu3**2)*muu4**2 - muu4**4 + muu2**2*(-1 + 2*muu3**2 + muu4**2 - 2*tred) + tred*(2 - 3*muu3**2 + 2*tred)))*
     &                     dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(muu1**2 + muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &           g1L*g2R*(2*g3L*g4R*muu1*muu3*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                  (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                    (-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g3R*g4L*(-(muu1*muu4*(muu1**2 + muu3**2 - tred)*dconjg(GR(indi,i,tt,gen))*
     &                    (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(-(muu2*muu3*(muu1**2 + muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                    (muu1**4*(muu2**2 - muu4**2 - tred) - 
     &                       (-1 + muu2**2 + muu3**2 - tred)*(muu2**2*muu3**2 + (-1 + muu4**2 - tred)*tred + muu3**2*(-muu4**2 + tred)) + 
     &                       muu1**2*((1 + 2*muu3**2)*muu4**2 - muu4**4 + muu2**2*(-1 + 2*muu3**2 + muu4**2 - 2*tred) + tred*(2 - 3*muu3**2 + 2*tred)))*
     &                     dconjg(HR(indj,i,tt,gen)))))))/2. + 
     &      (xd12*(g1R*g2L*(-2*g3R*g4L*muu3*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                  ((-muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(1 + muu2**2 - muu4**2 + tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                    muu2*(1 + muu1**2 - muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 muu1*dconjg(GL(indi,i,tt,gen))*(-((-1 + (1 + muu2**2)*muu3**2 - muu1**2*(-1 + muu3**2) + muu4**2 - tred)*
     &                       dconjg(HL(indj,i,tt,gen))) + muu2*muu3*(1 - muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen)))) + 
     &              g3L*g4R*(muu1*muu4*(1 + muu3**2 - muu4**2)*dconjg(GL(indi,i,tt,gen))*
     &                  (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GR(indi,i,tt,gen))*((-(muu2**4*muu3**2) - muu1**4*muu4**2 + muu2**2*muu3**2*(1 + 2*muu1**2 - muu3**2 + muu4**2) + 
     &                       (muu3**2 - muu4**2)*(muu3**2*(muu4**2 - tred) + tred*(1 - muu4**2 + tred)) + 
     &                       muu1**2*(-muu4**4 - 2*muu3**2*tred + muu4**2*(1 + muu3**2 + 2*tred)))*dconjg(HL(indj,i,tt,gen)) - 
     &                    muu2*muu3*(1 + muu3**2 - muu4**2)*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &           g1L*g2R*(-2*g3L*g4R*muu3*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                  (muu2*muu3*(1 - muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) - 
     &                    (-1 + (1 + muu2**2)*muu3**2 - muu1**2*(-1 + muu3**2) + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                 muu4*dconjg(GL(indi,i,tt,gen))*(muu2*(1 + muu1**2 - muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + 
     &                    (-muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(1 + muu2**2 - muu4**2 + tred))*dconjg(HR(indj,i,tt,gen)))) + 
     &              g3R*g4L*(muu1*muu4*(1 + muu3**2 - muu4**2)*dconjg(GR(indi,i,tt,gen))*
     &                  (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                 dconjg(GL(indi,i,tt,gen))*(-(muu2*muu3*(1 + muu3**2 - muu4**2)*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                    (-(muu2**4*muu3**2) - muu1**4*muu4**2 + muu2**2*muu3**2*(1 + 2*muu1**2 - muu3**2 + muu4**2) + 
     &                       (muu3**2 - muu4**2)*(muu3**2*(muu4**2 - tred) + tred*(1 - muu4**2 + tred)) + 
     &                       muu1**2*(-muu4**4 - 2*muu3**2*tred + muu4**2*(1 + muu3**2 + 2*tred)))*dconjg(HR(indj,i,tt,gen)))))))/2. + 
     &      (xd0*(-(g1R*(g2L*muu1*(g3L*(g4R*muu3*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                          (-(muu3*(muu2**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                         muu4*dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                            muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                      g4L*mGlr*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                          (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                         dconjg(GL(indi,i,tt,gen))*((1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                             dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &                   g3R*(g4R*mGlr*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                          (-(muu3*(muu2**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                         muu4*dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                            muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                      g4L*muu3*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                          (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                         dconjg(GL(indi,i,tt,gen))*((1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                             dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))))) + 
     &                g2R*M4Lr*(g3R*(g4L*muu3*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                          (-(muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                         muu1*dconjg(GL(indi,i,tt,gen))*(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) - 
     &                            muu3*(muu2**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                      g4R*mGlr*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                          (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                         dconjg(GR(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                            (1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &                   g3L*(g4L*mGlr*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                          (-(muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                         muu1*dconjg(GL(indi,i,tt,gen))*(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) - 
     &                            muu3*(muu2**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                      g4R*muu3*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                          (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                         dconjg(GR(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                            (1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*dconjg(HR(indj,i,tt,gen)))))))) - 
     &           g1L*(g2L*M4Lr*(g3L*(g4R*muu3*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                        (-(muu3*(muu2**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu4*dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                          muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4L*mGlr*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GL(indi,i,tt,gen))*((1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                           dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &                 g3R*(g4R*mGlr*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                        (-(muu3*(muu2**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu4*dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                          muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4L*muu3*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GL(indi,i,tt,gen))*((1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                           dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))))) + 
     &              g2R*muu1*(g3R*(g4L*muu3*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (-(muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu1*dconjg(GL(indi,i,tt,gen))*(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) - 
     &                          muu3*(muu2**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4R*mGlr*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                        (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GR(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                          (1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &                 g3L*(g4L*mGlr*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (-(muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu1*dconjg(GL(indi,i,tt,gen))*(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) - 
     &                          muu3*(muu2**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4R*muu3*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                        (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GR(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                          (1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*dconjg(HR(indj,i,tt,gen)))))))))/2. + 
     &      (xd3*(g1R*(g2R*M4Lr*(g3L*g4R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                     (-(muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                    dconjg(GR(indi,i,tt,gen))*(muu2*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))*
     &                        dconjg(HL(indj,i,tt,gen)) + muu3*(muu1**2*(1 - muu3**2 + muu4**2) - (-1 + muu2**2)*(-1 + muu2**2 + muu3**2 - tred))*
     &                        dconjg(HR(indj,i,tt,gen)))) + g3R*g4L*muu1*
     &                  (muu1*muu4*dconjg(GR(indi,i,tt,gen))*(2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                       (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                    dconjg(GL(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                       (muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))*dconjg(HR(indj,i,tt,gen))))) + 
     &              g2L*muu1*(g3R*(g4R*mGlr*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                        (muu3*(muu2**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu4*dconjg(GL(indi,i,tt,gen))*(-((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen))) + 
     &                          muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4L*(dconjg(GL(indi,i,tt,gen))*(muu3*
     &                           (-2 - muu2**4 + 2*muu3**2 + muu4**2 + muu1**2*(2 - 2*muu3**2 + muu4**2) - muu2**2*(-3 + muu3**2 + muu4**2 - tred) - 
     &                             2*tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu2*((-1 + muu2**2 + muu3**2 - tred)*(2*muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))*dconjg(HR(indj,i,tt,gen)))
     &                        + muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (muu3*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu2*(-muu1**2 - 3*muu3**2 + tred)*dconjg(HR(indj,i,tt,gen))))) - 
     &                 g3L*(g4L*mGlr*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GL(indi,i,tt,gen))*((1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                           dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4R*(muu1*dconjg(GR(indi,i,tt,gen))*(-((2*muu3**2*muu4**2 + muu1**2*(muu2**2 - tred) + muu2**2*(muu3**2 - tred) + tred - 
     &                               2*muu3**2*tred - muu4**2*tred + tred**2)*dconjg(HL(indj,i,tt,gen))) + 
     &                          muu2*muu3*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu4*dconjg(GL(indi,i,tt,gen))*((muu1**4 + (-1 + muu2**2)*muu3**2 + muu1**2*(-1 + muu3**2 + muu4**2 - tred))*
     &                           dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(-3*muu1**2 - muu3**2 + tred)*dconjg(HR(indj,i,tt,gen))))))) + 
     &           g1L*(g2L*M4Lr*(g3R*g4L*(dconjg(GL(indi,i,tt,gen))*
     &                     (muu3*(muu1**2*(1 - muu3**2 + muu4**2) - (-1 + muu2**2)*(-1 + muu2**2 + muu3**2 - tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                       muu2*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))*dconjg(HR(indj,i,tt,gen))) + 
     &                    muu1*muu4*dconjg(GR(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                       muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                 g3L*g4R*muu1*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                     (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                    dconjg(GR(indi,i,tt,gen))*((muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))*
     &                        dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &              g2R*muu1*(g3L*(g4L*mGlr*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu1*dconjg(GL(indi,i,tt,gen))*(-(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen))) + 
     &                          muu3*(muu2**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4R*(dconjg(GR(indi,i,tt,gen))*(muu2*((-1 + muu2**2 + muu3**2 - tred)*(2*muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))*
     &                           dconjg(HL(indj,i,tt,gen)) + muu3*
     &                           (-2 - muu2**4 + 2*muu3**2 + muu4**2 + muu1**2*(2 - 2*muu3**2 + muu4**2) - muu2**2*(-3 + muu3**2 + muu4**2 - tred) - 
     &                             2*tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                        (muu2*(-muu1**2 - 3*muu3**2 + tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu3*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))))) - 
     &                 g3R*(g4R*mGlr*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                        (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GR(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                          (1 - muu3**2 + muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4L*(muu4*dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(-3*muu1**2 - muu3**2 + tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                          (muu1**4 + (-1 + muu2**2)*muu3**2 + muu1**2*(-1 + muu3**2 + muu4**2 - tred))*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu1*dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 
     &                          (2*muu3**2*muu4**2 + muu1**2*(muu2**2 - tred) + muu2**2*(muu3**2 - tred) + tred - 2*muu3**2*tred - muu4**2*tred + tred**2)*
     &                           dconjg(HR(indj,i,tt,gen)))))))))/2. + 
     &      (xd1*(g1R*(g2R*M4Lr*muu3*(g3R*g4L*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                     (muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                    muu1*dconjg(GL(indi,i,tt,gen))*(-(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen))) + 
     &                       muu3*(muu2**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                 g3L*g4R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                     (-2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                    dconjg(GR(indi,i,tt,gen))*(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                       (-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &              g2L*(g3L*(g4L*mGlr*(muu1*dconjg(GL(indi,i,tt,gen))*
     &                        (((1 + muu2**2)*muu3**2 - muu1**2*(-1 + muu3**2 + muu4**2) - (-1 + muu4**2)*(-1 + muu4**2 - tred))*
     &                           dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu4*dconjg(GR(indi,i,tt,gen))*((muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu4**2 - 2*tred) + tred - 
     &                             muu3**2*tred - muu4**2*tred + tred**2)*dconjg(HL(indj,i,tt,gen)) - 
     &                          muu2*muu3*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4R*muu3*(dconjg(GR(indi,i,tt,gen))*(muu3*
     &                           (muu1**2*(2*muu2**2 + muu4**2 - 2*tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))*
     &                           dconjg(HL(indj,i,tt,gen)) - muu2*(muu1**2*(-1 + muu3**2 + muu4**2) + muu3**2*(-1 + muu2**2 + muu3**2 - tred))*
     &                           dconjg(HR(indj,i,tt,gen))) + muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                        (-(muu3*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                          muu2*(muu1**2 + 3*muu3**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &                 g3R*muu3*(g4R*mGlr*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                        (-((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GR(indi,i,tt,gen))*((muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))*
     &                           dconjg(HL(indj,i,tt,gen)) - muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4L*(muu1*dconjg(GL(indi,i,tt,gen))*((-2 + 2*muu3**2 + 3*muu4**2 - muu4**4 + muu2**2*(1 + muu3**2 - muu4**2) - 
     &                             muu1**2*(-2 + 2*muu3**2 + muu4**2) - 2*tred + muu4**2*tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu2*muu3*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu4*dconjg(GR(indi,i,tt,gen))*((2*muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-2 + 2*muu4**2 - 3*tred) + tred - 
     &                             muu3**2*tred - muu4**2*tred + tred**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu2*muu3*(-3*muu1**2 - muu3**2 + tred)*dconjg(HR(indj,i,tt,gen))))))) + 
     &           g1L*(g2L*M4Lr*muu3*(g3L*g4R*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                     (muu3*(muu2**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                    muu4*dconjg(GL(indi,i,tt,gen))*(-((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen))) + 
     &                       muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen)))) + 
     &                 g3R*g4L*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                     ((-1 + muu1**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                    dconjg(GL(indi,i,tt,gen))*((-1 + muu1**2 + muu2**2 + muu3**2 - muu1**2*muu3**2 + muu4**2 - muu2**2*muu4**2 - tred)*
     &                        dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HR(indj,i,tt,gen))))) + 
     &              g2R*(g3R*(g4R*mGlr*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                        (muu2*muu3*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                          ((1 + muu2**2)*muu3**2 - muu1**2*(-1 + muu3**2 + muu4**2) - (-1 + muu4**2)*(-1 + muu4**2 - tred))*dconjg(HR(indj,i,tt,gen))
     &                          ) + muu4*dconjg(GL(indi,i,tt,gen))*
     &                        (-(muu2*muu3*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                          (muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu4**2 - 2*tred) + tred - muu3**2*tred - muu4**2*tred + 
     &                             tred**2)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4L*muu3*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (muu2*(muu1**2 + 3*muu3**2 - tred)*dconjg(HL(indj,i,tt,gen)) - 
     &                          muu3*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GL(indi,i,tt,gen))*(-(muu2*(muu1**2*(-1 + muu3**2 + muu4**2) + muu3**2*(-1 + muu2**2 + muu3**2 - tred))*
     &                             dconjg(HL(indj,i,tt,gen))) + 
     &                          muu3*(muu1**2*(2*muu2**2 + muu4**2 - 2*tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))*
     &                           dconjg(HR(indj,i,tt,gen))))) + 
     &                 g3L*muu3*(g4R*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                        (muu2*muu3*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                          (-2 + 2*muu3**2 + 3*muu4**2 - muu4**4 + muu2**2*(1 + muu3**2 - muu4**2) - muu1**2*(-2 + 2*muu3**2 + muu4**2) - 2*tred + 
     &                             muu4**2*tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu4*dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(-3*muu1**2 - muu3**2 + tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                          (2*muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-2 + 2*muu4**2 - 3*tred) + tred - muu3**2*tred - muu4**2*tred + 
     &                             tred**2)*dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4L*mGlr*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) - (-1 + muu1**2 + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GL(indi,i,tt,gen))*(-(muu2*muu3*(-1 + muu2**2 + muu3**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                          (muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))*dconjg(HR(indj,i,tt,gen))))))))
     &         )/2. + (xd2*(g1R*(g2R*M4Lr*(-(g3L*g4R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                       (muu2*(1 + muu3**2 - muu4**2)*dconjg(HL(indj,i,tt,gen)) + (1 - muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                      dconjg(GR(indi,i,tt,gen))*(muu2*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 - muu2**2 + muu4**2 + tred))*
     &                          dconjg(HL(indj,i,tt,gen)) + muu3*(1 - muu2**2 - muu3**2 + muu1**2*(-1 + muu3**2 - muu4**2) + tred)*
     &                          dconjg(HR(indj,i,tt,gen))))) + 
     &                 g3R*g4L*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                     (muu2*(1 + muu1**2 - muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) - 
     &                       (-1 + muu1**2 + muu4**2 + muu2**2*(1 + muu3**2 - muu4**2) - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                    muu1*dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(1 - muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) + 
     &                       (muu2**2*muu3**2 - muu4**4 - muu3**2*tred + muu4**2*(1 - muu1**2 + muu3**2 + tred))*dconjg(HR(indj,i,tt,gen))))) + 
     &              g2L*(-(g3L*(g4L*mGlr*(muu4*dconjg(GR(indi,i,tt,gen))*
     &                          ((-muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(1 + muu2**2 - muu4**2 + tred))*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*(1 + muu1**2 - muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                         muu1*dconjg(GL(indi,i,tt,gen))*(-((-1 + (1 + muu2**2)*muu3**2 - muu1**2*(-1 + muu3**2) + muu4**2 - tred)*
     &                               dconjg(HL(indj,i,tt,gen))) + muu2*muu3*(1 - muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen)))) - 
     &                      g4R*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                          (-((-1 + muu2**2 - muu3**2 + muu1**2*(1 + muu3**2) + muu4**2 - muu2**2*muu4**2 - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                            muu2*muu3*(2 + muu1**2 - muu2**2 + muu3**2 - muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                         dconjg(GR(indi,i,tt,gen))*((-(muu2**4*muu3**2) + muu2**2*muu3**2*(1 + 2*muu1**2 - muu3**2 + tred) + 
     &                               muu1**2*(-muu4**4 - 2*muu3**2*tred + muu4**2*(1 - muu1**2 + 2*muu3**2 + tred)))*dconjg(HL(indj,i,tt,gen)) + 
     &                            muu2*muu3*(1 - muu3**2 - muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*dconjg(HR(indj,i,tt,gen))))))
     &                  + g3R*(g4R*mGlr*(muu1*muu4*dconjg(GL(indi,i,tt,gen))*
     &                        ((1 - muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + muu2*(1 + muu3**2 - muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GR(indi,i,tt,gen))*(muu3*(-muu2**4 + muu1**2*(muu4**2 - tred) + muu2**2*(1 + muu1**2 - muu3**2 + tred))*
     &                           dconjg(HL(indj,i,tt,gen)) + muu2*(1 - muu3**2 - muu4**2 - muu1**2*muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                           dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4L*(muu4*dconjg(GR(indi,i,tt,gen))*(muu3*(2*muu1**4 + muu2**2*(-muu3**2 + tred) - muu1**2*(2 + 2*muu2**2 - muu4**2 + tred))*
     &                           dconjg(HL(indj,i,tt,gen)) + muu2*((-1 + muu2**2)*muu3**2 + muu1**2*(-1 - 2*muu3**2 + muu4**2))*dconjg(HR(indj,i,tt,gen)))
     &                        + muu1*dconjg(GL(indi,i,tt,gen))*
     &                        (muu3*(-2 + 2*muu3**2 + muu2**2*(1 + muu3**2) + muu4**2 + muu1**2*(2 - 2*muu3**2 + muu4**2) - 2*tred)*
     &                           dconjg(HL(indj,i,tt,gen)) + muu2*(2*muu3**4 + muu3**2*(-2 + muu2**2 - 2*muu4**2 - tred) + muu4**2*(-muu1**2 + tred))*
     &                           dconjg(HR(indj,i,tt,gen))))))) + 
     &           g1L*(g2L*M4Lr*(g3L*g4R*(muu4*dconjg(GL(indi,i,tt,gen))*
     &                     (-((-1 + muu1**2 + muu4**2 + muu2**2*(1 + muu3**2 - muu4**2) - tred)*dconjg(HL(indj,i,tt,gen))) + 
     &                       muu2*(1 + muu1**2 - muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                    muu1*dconjg(GR(indi,i,tt,gen))*((muu2**2*muu3**2 - muu4**4 - muu3**2*tred + muu4**2*(1 - muu1**2 + muu3**2 + tred))*
     &                        dconjg(HL(indj,i,tt,gen)) + muu2*muu3*(1 - muu3**2 + muu4**2)*dconjg(HR(indj,i,tt,gen)))) - 
     &                 g3R*g4L*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                     ((1 - muu1**2 + muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + muu2*(1 + muu3**2 - muu4**2)*dconjg(HR(indj,i,tt,gen))) + 
     &                    dconjg(GL(indi,i,tt,gen))*(muu3*(1 - muu2**2 - muu3**2 + muu1**2*(-1 + muu3**2 - muu4**2) + tred)*dconjg(HL(indj,i,tt,gen)) + 
     &                       muu2*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 - muu2**2 + muu4**2 + tred))*dconjg(HR(indj,i,tt,gen))))) + 
     &              g2R*(g3L*(g4L*mGlr*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (muu2*(1 + muu3**2 - muu4**2)*dconjg(HL(indj,i,tt,gen)) + (1 - muu1**2 + muu2**2)*muu3*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GL(indi,i,tt,gen))*(muu2*(1 - muu3**2 - muu4**2 - muu1**2*muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                           dconjg(HL(indj,i,tt,gen)) + muu3*(-muu2**4 + muu1**2*(muu4**2 - tred) + muu2**2*(1 + muu1**2 - muu3**2 + tred))*
     &                           dconjg(HR(indj,i,tt,gen)))) + 
     &                    g4R*(muu1*dconjg(GR(indi,i,tt,gen))*(muu2*(2*muu3**4 + muu3**2*(-2 + muu2**2 - 2*muu4**2 - tred) + muu4**2*(-muu1**2 + tred))*
     &                           dconjg(HL(indj,i,tt,gen)) + muu3*
     &                           (-2 + 2*muu3**2 + muu2**2*(1 + muu3**2) + muu4**2 + muu1**2*(2 - 2*muu3**2 + muu4**2) - 2*tred)*dconjg(HR(indj,i,tt,gen)))
     &                         + muu4*dconjg(GL(indi,i,tt,gen))*
     &                        (muu2*((-1 + muu2**2)*muu3**2 + muu1**2*(-1 - 2*muu3**2 + muu4**2))*dconjg(HL(indj,i,tt,gen)) + 
     &                          muu3*(2*muu1**4 + muu2**2*(-muu3**2 + tred) - muu1**2*(2 + 2*muu2**2 - muu4**2 + tred))*dconjg(HR(indj,i,tt,gen))))) - 
     &                 g3R*(g4R*mGlr*(muu1*dconjg(GR(indi,i,tt,gen))*
     &                        (muu2*muu3*(1 - muu3**2 + muu4**2)*dconjg(HL(indj,i,tt,gen)) - 
     &                          (-1 + (1 + muu2**2)*muu3**2 - muu1**2*(-1 + muu3**2) + muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                       muu4*dconjg(GL(indi,i,tt,gen))*(muu2*(1 + muu1**2 - muu2**2)*muu3*dconjg(HL(indj,i,tt,gen)) + 
     &                          (-muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(1 + muu2**2 - muu4**2 + tred))*dconjg(HR(indj,i,tt,gen)))) - 
     &                    g4L*(muu1*muu4*dconjg(GR(indi,i,tt,gen))*
     &                        (muu2*muu3*(2 + muu1**2 - muu2**2 + muu3**2 - muu4**2)*dconjg(HL(indj,i,tt,gen)) - 
     &                          (-1 + muu2**2 - muu3**2 + muu1**2*(1 + muu3**2) + muu4**2 - muu2**2*muu4**2 - tred)*dconjg(HR(indj,i,tt,gen))) + 
     &                       dconjg(GL(indi,i,tt,gen))*(muu2*muu3*(1 - muu3**2 - muu1**2*(-1 + muu3**2) - muu4**2 + muu2**2*(-1 + muu4**2) + tred)*
     &                           dconjg(HL(indj,i,tt,gen)) + (-(muu2**4*muu3**2) + muu2**2*muu3**2*(1 + 2*muu1**2 - muu3**2 + tred) + 
     &                             muu1**2*(-muu4**4 - 2*muu3**2*tred + muu4**2*(1 - muu1**2 + 2*muu3**2 + tred)))*dconjg(HR(indj,i,tt,gen)))))))))/2.))/
     &  (4.*Pi**2*(-ured + MSf(i,itypeU,gen)**2/sqrtS**2))
	
	enddo

	tempBoxu = tempBoxu + ChiChi2QQ_M2BoxGluino - tempBox2
	tempBox2 = ChiChi2QQ_M2BoxGluino
	
! T-channel contr.

	do i = 1,2

	ChiChi2QQ_M2BoxGluino = ChiChi2QQ_M2BoxGluino + (CF*NC*sqrtS**4*((xd0*(-(dconjg(ELc(indi,i,tt,gen))*((g1R*muu1*
     &                    (g2R*M4Lr*muu2*(g3R*muu3*(2*g4R*mGlr*muu4 + g4L*(-1 + muu3**2 + muu4**2)) + 
     &                         g3L*(2*g4R*muu3**2*muu4 + g4L*mGlr*(-1 + muu3**2 + muu4**2))) + 
     &                      g2L*(g3R*muu3*(g4R*mGlr*(-1 + muu1**2 + muu2**2)*muu4 + g4L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)) + 
     &                         g3L*(g4R*(-1 + muu1**2 + muu2**2)*muu3**2*muu4 + g4L*mGlr*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))) + 
     &                   g1L*(g2R*muu1**2*muu2*(g3R*muu3*(2*g4R*mGlr*muu4 + g4L*(-1 + muu3**2 + muu4**2)) + 
     &                         g3L*(2*g4R*muu3**2*muu4 + g4L*mGlr*(-1 + muu3**2 + muu4**2))) + 
     &                      g2L*M4Lr*(g3R*muu3*(g4R*mGlr*(-1 + muu1**2 + muu2**2)*muu4 + g4L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)) + 
     &                         g3L*(g4R*(-1 + muu1**2 + muu2**2)*muu3**2*muu4 + g4L*mGlr*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))))*
     &                 dconjg(FL(indj,i,tt,gen)) + (g1R*muu1*(g2L*muu2*
     &                       (g3R*muu3*(g4L*muu4*(muu1**2 + muu3**2 - tred) + g4R*mGlr*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                         g3L*(g4L*mGlr*muu4*(muu1**2 + muu3**2 - tred) + g4R*muu3**2*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                      g2R*M4Lr*(g3R*muu3*(g4L*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4R*mGlr*(muu2**2 + muu4**2 - tred)) + 
     &                         g3L*(g4L*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4R*muu3**2*(muu2**2 + muu4**2 - tred)))) + 
     &                   g1L*(g2L*M4Lr*muu2*(g3R*muu3*(g4L*muu4*(muu1**2 + muu3**2 - tred) + g4R*mGlr*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                         g3L*(g4L*mGlr*muu4*(muu1**2 + muu3**2 - tred) + g4R*muu3**2*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                      g2R*muu1**2*(g3R*muu3*(g4L*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4R*mGlr*(muu2**2 + muu4**2 - tred)) + 
     &                         g3L*(g4L*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4R*muu3**2*(muu2**2 + muu4**2 - tred)))))*dconjg(FR(indj,i,tt,gen))))
     &            - dconjg(ER(indi,i,tt,gen))*((g1L*muu1*(g2R*muu2*
     &                     (g3L*muu3*(g4R*muu4*(muu1**2 + muu3**2 - tred) + g4L*mGlr*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                       g3R*(g4R*mGlr*muu4*(muu1**2 + muu3**2 - tred) + g4L*muu3**2*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                    g2L*M4Lr*(g3L*muu3*(g4R*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4L*mGlr*(muu2**2 + muu4**2 - tred)) + 
     &                       g3R*(g4R*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4L*muu3**2*(muu2**2 + muu4**2 - tred)))) + 
     &                 g1R*(g2R*M4Lr*muu2*(g3L*muu3*(g4R*muu4*(muu1**2 + muu3**2 - tred) + g4L*mGlr*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                       g3R*(g4R*mGlr*muu4*(muu1**2 + muu3**2 - tred) + g4L*muu3**2*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                    g2L*muu1**2*(g3L*muu3*(g4R*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4L*mGlr*(muu2**2 + muu4**2 - tred)) + 
     &                       g3R*(g4R*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + g4L*muu3**2*(muu2**2 + muu4**2 - tred)))))*dconjg(FL(indj,i,tt,gen)) + 
     &              (g1L*muu1*(g2L*M4Lr*muu2*(g3L*muu3*(2*g4L*mGlr*muu4 + g4R*(-1 + muu3**2 + muu4**2)) + 
     &                       g3R*(2*g4L*muu3**2*muu4 + g4R*mGlr*(-1 + muu3**2 + muu4**2))) + 
     &                    g2R*(g3L*muu3*(g4L*mGlr*(-1 + muu1**2 + muu2**2)*muu4 + g4R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)) + 
     &                       g3R*(g4L*(-1 + muu1**2 + muu2**2)*muu3**2*muu4 + g4R*mGlr*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))) + 
     &                 g1R*(g2L*muu1**2*muu2*(g3L*muu3*(2*g4L*mGlr*muu4 + g4R*(-1 + muu3**2 + muu4**2)) + 
     &                       g3R*(2*g4L*muu3**2*muu4 + g4R*mGlr*(-1 + muu3**2 + muu4**2))) + 
     &                    g2R*M4Lr*(g3L*muu3*(g4L*mGlr*(-1 + muu1**2 + muu2**2)*muu4 + g4R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)) + 
     &                       g3R*(g4L*(-1 + muu1**2 + muu2**2)*muu3**2*muu4 + g4R*mGlr*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)))))*
     &               dconjg(FR(indj,i,tt,gen)))))/2. + (xd00*(dconjg(ER(indi,i,tt,gen))*
     &            ((g1L*g2R*muu1*muu2*(-4*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &                 g1R*g2L*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 - g3R*g4L*(muu1**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred)))*
     &               dconjg(FL(indj,i,tt,gen)) + (g1R*g2L*muu2*
     &                  (-2*g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + g3L*g4R*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                 g1L*g2R*muu1*(g3R*g4L*muu4*(-1 + muu1**2 + muu4**2 - tred) - 2*g3L*g4R*muu3*(muu2**2 + muu4**2 - tred)))*dconjg(FR(indj,i,tt,gen)))
     &            + dconjg(ELc(indi,i,tt,gen))*((g1L*g2R*muu2*
     &                  (-2*g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + g3R*g4L*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                 g1R*g2L*muu1*(g3L*g4R*muu4*(-1 + muu1**2 + muu4**2 - tred) - 2*g3R*g4L*muu3*(muu2**2 + muu4**2 - tred)))*dconjg(FL(indj,i,tt,gen))
     &               + (g1R*g2L*muu1*muu2*(-4*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &                 g1L*g2R*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 - g3L*g4R*(muu1**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred)))*
     &               dconjg(FR(indj,i,tt,gen)))))/sqrtS**2 - (muu1*xd33*
     &         (dconjg(ER(indi,i,tt,gen))*((g1L*g2R*muu2*(2*g3L*g4R*muu1**2*muu3*muu4 + 
     &                    g3R*g4L*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))) + 
     &                 g1R*g2L*muu1*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                    g3R*g4L*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))))*dconjg(FL(indj,i,tt,gen))
     &               + (g1R*g2L*muu1*muu2*(g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + g3L*g4R*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                 g1L*g2R*(g3R*g4L*muu4*(muu1**2*(muu2**2 + muu3**2 - muu4**2) + (-1 + muu2**2)*(muu3**2 - tred)) + 
     &                    g3L*g4R*muu1**2*muu3*(muu2**2 + muu4**2 - tred)))*dconjg(FR(indj,i,tt,gen))) + 
     &           dconjg(ELc(indi,i,tt,gen))*((g1L*g2R*muu1*muu2*(g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + g3R*g4L*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                 g1R*g2L*(g3L*g4R*muu4*(muu1**2*(muu2**2 + muu3**2 - muu4**2) + (-1 + muu2**2)*(muu3**2 - tred)) + 
     &                    g3R*g4L*muu1**2*muu3*(muu2**2 + muu4**2 - tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &              (g1R*g2L*muu2*(2*g3R*g4L*muu1**2*muu3*muu4 + 
     &                    g3L*g4R*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))) + 
     &                 g1L*g2R*muu1*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                    g3L*g4R*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))))*dconjg(FR(indj,i,tt,gen))))
     &         )/2. - (muu3*xd11*(dconjg(ER(indi,i,tt,gen))*((g1L*g2R*muu1*muu2*muu3*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &                 g1R*g2L*(g3L*g4R*muu4*(muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu4**2 - 2*tred) + tred - muu3**2*tred - 
     &                       muu4**2*tred + tred**2) + g3R*g4L*muu3*
     &                     (muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))))*dconjg(FL(indj,i,tt,gen)) + 
     &              (g1L*g2R*muu1*muu3*(g3R*g4L*muu4*(-1 + muu1**2 + muu4**2 - tred) + g3L*g4R*muu3*(muu2**2 + muu4**2 - tred)) + 
     &                 g1R*g2L*muu2*(g3R*g4L*muu3*muu4*(muu1**2 + muu3**2 - tred) + 
     &                    g3L*g4R*(-(muu2**2*muu3**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + tred - muu4**2*tred)))*
     &               dconjg(FR(indj,i,tt,gen))) + dconjg(ELc(indi,i,tt,gen))*
     &            ((g1R*g2L*muu1*muu3*(g3L*g4R*muu4*(-1 + muu1**2 + muu4**2 - tred) + g3R*g4L*muu3*(muu2**2 + muu4**2 - tred)) + 
     &                 g1L*g2R*muu2*(g3L*g4R*muu3*muu4*(muu1**2 + muu3**2 - tred) + 
     &                    g3R*g4L*(-(muu2**2*muu3**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + tred - muu4**2*tred)))*
     &               dconjg(FL(indj,i,tt,gen)) + (g1R*g2L*muu1*muu2*muu3*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &                 g1L*g2R*(g3R*g4L*muu4*(muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu4**2 - 2*tred) + tred - muu3**2*tred - 
     &                       muu4**2*tred + tred**2) + g3L*g4R*muu3*
     &                     (muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))))*dconjg(FR(indj,i,tt,gen)))))/2. + 
     &      (xd13*(-(dconjg(ER(indi,i,tt,gen))*((2*g1L*g2R*muu1*muu2*muu3*
     &                    (g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + g3R*g4L*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                   g1R*g2L*(2*g3L*g4R*muu1**2*muu3*muu4*(-1 + muu1**2 + muu4**2 - tred) + 
     &                      g3R*g4L*(muu1**2 + muu3**2 - tred)*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))
     &                   )*dconjg(FL(indj,i,tt,gen)) + (g1R*g2L*muu2*
     &                    (2*g3L*g4R*muu1**2*muu3*(-1 + muu3**2 + muu4**2) + g3R*g4L*muu4*(muu1**2 + muu3**2 - tred)**2) + 
     &                   g1L*g2R*muu1*muu3*(2*g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 + g3L*g4R*(muu1**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred)))*
     &                 dconjg(FR(indj,i,tt,gen)))) - dconjg(ELc(indi,i,tt,gen))*
     &            ((g1L*g2R*muu2*(2*g3R*g4L*muu1**2*muu3*(-1 + muu3**2 + muu4**2) + g3L*g4R*muu4*(muu1**2 + muu3**2 - tred)**2) + 
     &                 g1R*g2L*muu1*muu3*(2*g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 + g3R*g4L*(muu1**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred)))*
     &               dconjg(FL(indj,i,tt,gen)) + (2*g1R*g2L*muu1*muu2*muu3*
     &                  (g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + g3L*g4R*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                 g1L*g2R*(2*g3R*g4L*muu1**2*muu3*muu4*(-1 + muu1**2 + muu4**2 - tred) + 
     &                    g3L*g4R*(muu1**2 + muu3**2 - tred)*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))))*
     &               dconjg(FR(indj,i,tt,gen)))))/2. + (xd23*(-(dconjg(ER(indi,i,tt,gen))*
     &              ((2*g1L*g2R*muu1*muu2*(g3L*g4R*(1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                      g3R*g4L*(muu3**4 + muu3**2*(-1 + muu2**2 - muu4**2 - tred) + muu4**2*(-muu1**2 + tred))) + 
     &                   g1R*g2L*(2*g3L*g4R*muu1**2*(-1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                      g3R*g4L*(-(muu2**2*(muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred)) + muu1**4*(muu2**2 - muu4**2 - tred) + 
     &                         muu1**2*(-muu2**4 + muu4**2 + 2*muu3**2*muu4**2 - muu4**4 + muu2**2*(muu3**2 + muu4**2) + tred - 2*muu3**2*tred + tred**2)))
     &                   )*dconjg(FL(indj,i,tt,gen)) + (g1R*g2L*muu2*
     &                    (2*g3L*g4R*muu1**2*muu3*(-1 + muu3**2 - muu4**2) + g3R*g4L*(1 + muu1**2 - muu2**2)*muu4*(muu1**2 + muu3**2 - tred)) + 
     &                   g1L*g2R*muu1*(g3L*g4R*(1 + muu1**2 - muu2**2)*muu3*(muu2**2 + muu4**2 - tred) + 
     &                      2*g3R*g4L*muu4*((-1 + muu1**2)*muu3**2 - muu1**2*muu4**2 + tred)))*dconjg(FR(indj,i,tt,gen)))) - 
     &           dconjg(ELc(indi,i,tt,gen))*((g1L*g2R*muu2*(2*g3R*g4L*muu1**2*muu3*(-1 + muu3**2 - muu4**2) + 
     &                    g3L*g4R*(1 + muu1**2 - muu2**2)*muu4*(muu1**2 + muu3**2 - tred)) + 
     &                 g1R*g2L*muu1*(g3R*g4L*(1 + muu1**2 - muu2**2)*muu3*(muu2**2 + muu4**2 - tred) + 
     &                    2*g3L*g4R*muu4*((-1 + muu1**2)*muu3**2 - muu1**2*muu4**2 + tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &              (-2*g1R*g2L*muu1*muu2*(g3R*g4L*(-1 - muu1**2 + muu2**2)*muu3*muu4 + 
     &                    g3L*g4R*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 - muu2**2 + muu4**2 + tred))) + 
     &                 g1L*g2R*(2*g3R*g4L*muu1**2*(-1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                    g3L*g4R*(-(muu2**2*(muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred)) + muu1**4*(muu2**2 - muu4**2 - tred) + 
     &                       muu1**2*(-muu2**4 + muu4**2 + 2*muu3**2*muu4**2 - muu4**4 + muu2**2*(muu3**2 + muu4**2) + tred - 2*muu3**2*tred + tred**2))))*
     &               dconjg(FR(indj,i,tt,gen)))))/2. + (xd22*(-(dconjg(ER(indi,i,tt,gen))*
     &              ((g1L*g2R*muu1*muu2*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(muu3**4 + muu4**2*(-1 + muu4**2) - muu3**2*(1 + 2*muu4**2))) + 
     &                   g1R*g2L*(g3L*g4R*(muu1**4 + muu2**2*(-1 + muu2**2) - muu1**2*(1 + 2*muu2**2))*muu3*muu4 + 
     &                      g3R*g4L*(-(muu2**4*muu3**2) + muu1**2*(-muu4**4 - muu3**2*tred + muu4**2*(1 - muu1**2 + muu3**2 + tred)) + 
     &                         muu2**2*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 + muu1**2 + muu4**2 + tred)))))*dconjg(FL(indj,i,tt,gen)) + 
     &                (g1L*g2R*muu1*(g3L*g4R*muu3*(muu2**2 + muu4**2 - tred) + 
     &                      g3R*g4L*muu4*((-1 + muu1**2)*muu3**2 - muu1**2*muu4**2 + muu2**2*(-1 - muu3**2 + muu4**2) + tred)) + 
     &                   g1R*g2L*muu2*(g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + 
     &                      g3L*g4R*muu3*(-muu4**2 + muu1**2*(-1 + muu3**2 - muu4**2) + muu2**2*(-muu3**2 + muu4**2) + tred)))*dconjg(FR(indj,i,tt,gen))))
     &            - dconjg(ELc(indi,i,tt,gen))*((g1R*g2L*muu1*(g3R*g4L*muu3*(muu2**2 + muu4**2 - tred) + 
     &                    g3L*g4R*muu4*((-1 + muu1**2)*muu3**2 - muu1**2*muu4**2 + muu2**2*(-1 - muu3**2 + muu4**2) + tred)) + 
     &                 g1L*g2R*muu2*(g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + 
     &                    g3R*g4L*muu3*(-muu4**2 + muu1**2*(-1 + muu3**2 - muu4**2) + muu2**2*(-muu3**2 + muu4**2) + tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &              (g1R*g2L*muu1*muu2*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(muu3**4 + muu4**2*(-1 + muu4**2) - muu3**2*(1 + 2*muu4**2))) + 
     &                 g1L*g2R*(g3R*g4L*(muu1**4 + muu2**2*(-1 + muu2**2) - muu1**2*(1 + 2*muu2**2))*muu3*muu4 + 
     &                    g3L*g4R*(-(muu2**4*muu3**2) + muu1**2*(-muu4**4 - muu3**2*tred + muu4**2*(1 - muu1**2 + muu3**2 + tred)) + 
     &                       muu2**2*(-muu3**4 + muu4**2*(muu1**2 - tred) + muu3**2*(1 + muu1**2 + muu4**2 + tred)))))*dconjg(FR(indj,i,tt,gen)))))/2. + 
     &      (xd12*(-(dconjg(ER(indi,i,tt,gen))*((2*g1L*g2R*muu1*muu2*muu3*(g3R*g4L*muu3*(-1 + muu3**2 - muu4**2) + g3L*g4R*muu4*(1 + muu3**2 - muu4**2)) + 
     &                   g1R*g2L*(-2*g3L*g4R*muu3*muu4*(-muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(1 + muu2**2 - muu4**2 + tred)) + 
     &                      g3R*g4L*(-(muu2**4*muu3**2) - muu1**4*muu4**2 + muu2**2*muu3**2*(1 + 2*muu1**2 - muu3**2 + muu4**2) + 
     &                         (muu3**2 - muu4**2)*(muu3**2*(muu4**2 - tred) + tred*(1 - muu4**2 + tred)) + 
     &                         muu1**2*(-muu4**4 - 2*muu3**2*tred + muu4**2*(1 + muu3**2 + 2*tred)))))*dconjg(FL(indj,i,tt,gen)) + 
     &                (g1L*g2R*muu1*muu3*(2*g3R*g4L*(-1 + muu1**2 - muu2**2)*muu3*muu4 + g3L*g4R*(1 + muu3**2 - muu4**2)*(muu2**2 + muu4**2 - tred)) + 
     &                   g1R*g2L*muu2*(g3R*g4L*muu4*(1 + muu3**2 - muu4**2)*(muu1**2 + muu3**2 - tred) + 
     &                      2*g3L*g4R*muu3*(-(muu2**2*muu3**2) + muu1**2*(-1 + muu3**2) + tred)))*dconjg(FR(indj,i,tt,gen)))) - 
     &           dconjg(ELc(indi,i,tt,gen))*((g1R*g2L*muu1*muu3*
     &                  (2*g3L*g4R*(-1 + muu1**2 - muu2**2)*muu3*muu4 + g3R*g4L*(1 + muu3**2 - muu4**2)*(muu2**2 + muu4**2 - tred)) + 
     &                 g1L*g2R*muu2*(g3L*g4R*muu4*(1 + muu3**2 - muu4**2)*(muu1**2 + muu3**2 - tred) + 
     &                    2*g3R*g4L*muu3*(-(muu2**2*muu3**2) + muu1**2*(-1 + muu3**2) + tred)))*dconjg(FL(indj,i,tt,gen)) + 
     &              (2*g1R*g2L*muu1*muu2*muu3*(g3L*g4R*muu3*(-1 + muu3**2 - muu4**2) + g3R*g4L*muu4*(1 + muu3**2 - muu4**2)) + 
     &                 g1L*g2R*(-2*g3R*g4L*muu3*muu4*(-muu1**4 + muu2**2*(muu3**2 - tred) + muu1**2*(1 + muu2**2 - muu4**2 + tred)) + 
     &                    g3L*g4R*(-(muu2**4*muu3**2) - muu1**4*muu4**2 + muu2**2*muu3**2*(1 + 2*muu1**2 - muu3**2 + muu4**2) + 
     &                       (muu3**2 - muu4**2)*(muu3**2*(muu4**2 - tred) + tred*(1 - muu4**2 + tred)) + 
     &                       muu1**2*(-muu4**4 - 2*muu3**2*tred + muu4**2*(1 + muu3**2 + 2*tred)))))*dconjg(FR(indj,i,tt,gen)))))/2. + 
     &      (xd3*(-(dconjg(ER(indi,i,tt,gen))*((g1L*muu1*(g2R*muu2*
     &                       (g3R*(g4L*((-1 + muu2**2 + muu3**2 - tred)*(2*muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred)) + 
     &                            g4R*mGlr*muu4*(muu1**2 + muu3**2 - tred)) + 
     &                         g3L*muu3*(g4R*muu4*(3*muu1**2 + muu3**2 - tred) + g4L*mGlr*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                      g2L*M4Lr*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                         g3R*g4L*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))) + 
     &                   g1R*(g2R*M4Lr*muu2*(2*g3L*g4R*muu1**2*muu3*muu4 + 
     &                         g3R*g4L*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))) + 
     &                      g2L*muu1**2*(g3L*muu3*(g4L*mGlr*(muu2**2 + muu4**2 - tred) + g4R*muu4*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred)) + 
     &                         g3R*(g4R*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + 
     &                            g4L*(2*muu3**2*muu4**2 + muu1**2*(muu2**2 - tred) + muu2**2*(muu3**2 - tred) + tred - 2*muu3**2*tred - muu4**2*tred + 
     &                               tred**2)))))*dconjg(FL(indj,i,tt,gen)) + 
     &                (g1L*muu1*(g2R*(g3R*(g4L*muu4*(muu1**2*(muu2**2 + 2*muu3**2 - muu4**2) + (-1 + muu2**2)*(2*muu3**2 - tred)) + 
     &                            g4R*mGlr*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)) + 
     &                         g3L*muu3*(g4L*mGlr*(-1 + muu1**2 + muu2**2)*muu4 + g4R*(muu2**2*muu3**2 + muu1**2*(muu2**2 + 2*muu4**2 - tred) - tred))) + 
     &                      g2L*M4Lr*muu2*(g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + g3L*g4R*muu3*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                   g1R*(g2L*muu1**2*muu2*(g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2) + g4L*muu4*(muu1**2 + 3*muu3**2 - tred)) + 
     &                         g3L*muu3*(2*g4L*mGlr*muu4 + g4R*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred))) + 
     &                      g2R*M4Lr*(g3R*g4L*muu4*(muu1**2*(muu2**2 + muu3**2 - muu4**2) + (-1 + muu2**2)*(muu3**2 - tred)) + 
     &                         g3L*g4R*muu1**2*muu3*(muu2**2 + muu4**2 - tred))))*dconjg(FR(indj,i,tt,gen)))) - 
     &           dconjg(ELc(indi,i,tt,gen))*((g1R*muu1*(g2L*(g3L*
     &                        (g4R*muu4*(muu1**2*(muu2**2 + 2*muu3**2 - muu4**2) + (-1 + muu2**2)*(2*muu3**2 - tred)) + 
     &                          g4L*mGlr*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)) + 
     &                       g3R*muu3*(g4R*mGlr*(-1 + muu1**2 + muu2**2)*muu4 + g4L*(muu2**2*muu3**2 + muu1**2*(muu2**2 + 2*muu4**2 - tred) - tred))) + 
     &                    g2R*M4Lr*muu2*(g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + g3R*g4L*muu3*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                 g1L*(g2R*muu1**2*muu2*(g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2) + g4R*muu4*(muu1**2 + 3*muu3**2 - tred)) + 
     &                       g3R*muu3*(2*g4R*mGlr*muu4 + g4L*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred))) + 
     &                    g2L*M4Lr*(g3L*g4R*muu4*(muu1**2*(muu2**2 + muu3**2 - muu4**2) + (-1 + muu2**2)*(muu3**2 - tred)) + 
     &                       g3R*g4L*muu1**2*muu3*(muu2**2 + muu4**2 - tred))))*dconjg(FL(indj,i,tt,gen)) + 
     &              (g1R*muu1*(g2L*muu2*(g3L*(g4R*((-1 + muu2**2 + muu3**2 - tred)*(2*muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred)) + 
     &                          g4L*mGlr*muu4*(muu1**2 + muu3**2 - tred)) + 
     &                       g3R*muu3*(g4L*muu4*(3*muu1**2 + muu3**2 - tred) + g4R*mGlr*(-1 + muu2**2 + muu3**2 - tred))) + 
     &                    g2R*M4Lr*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 + 
     &                       g3L*g4R*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred)))) + 
     &                 g1L*(g2L*M4Lr*muu2*(2*g3R*g4L*muu1**2*muu3*muu4 + 
     &                       g3L*g4R*((muu3**2 - tred)*(-1 + muu2**2 + muu3**2 - tred) + muu1**2*(muu2**2 - muu4**2 - tred))) + 
     &                    g2R*muu1**2*(g3R*muu3*(g4R*mGlr*(muu2**2 + muu4**2 - tred) + g4L*muu4*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred)) + 
     &                       g3L*(g4L*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + 
     &                          g4R*(2*muu3**2*muu4**2 + muu1**2*(muu2**2 - tred) + muu2**2*(muu3**2 - tred) + tred - 2*muu3**2*tred - muu4**2*tred + 
     &                             tred**2)))))*dconjg(FR(indj,i,tt,gen)))))/2. + 
     &      (xd1*(-(dconjg(ER(indi,i,tt,gen))*((g1L*muu1*muu3*(g2R*muu2*
     &                       (g3L*(g4L*mGlr*(-1 + muu3**2 + muu4**2) + g4R*muu4*(muu1**2 + 3*muu3**2 - tred)) + 
     &                         g3R*muu3*(2*g4R*mGlr*muu4 + g4L*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred))) + 
     &                      g2L*M4Lr*(g3L*g4R*muu4*(-1 + muu1**2 + muu4**2 - tred) + g3R*g4L*muu3*(muu2**2 + muu4**2 - tred))) + 
     &                   g1R*(g2R*M4Lr*muu2*muu3*(g3L*g4R*muu4*(muu1**2 + muu3**2 - tred) + g3R*g4L*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                      g2L*(g3R*(g4R*mGlr*muu4*(muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu4**2 - 2*tred) + tred - muu3**2*tred - 
     &                               muu4**2*tred + tred**2) + 
     &                            g4L*muu3**2*(muu1**2*(2*muu2**2 + muu4**2 - 2*tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))) + 
     &                         g3L*muu3*(g4R*muu4*(2*muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-2 + 2*muu4**2 - 3*tred) + tred - 
     &                               muu3**2*tred - muu4**2*tred + tred**2) + 
     &                            g4L*mGlr*(muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))))))*
     &                 dconjg(FL(indj,i,tt,gen)) + (g1L*muu1*muu3*
     &                    (g2L*M4Lr*muu2*(2*g3R*g4L*muu3*muu4 + g3L*g4R*(-1 + muu3**2 + muu4**2)) + 
     &                      g2R*(g3R*muu3*(g4R*mGlr*(muu2**2 + muu4**2 - tred) + g4L*muu4*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred)) + 
     &                         g3L*(g4L*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + 
     &                            g4R*(2*muu2**2*muu3**2 + muu1**2*muu4**2 + muu3**2*muu4**2 - tred - muu3**2*tred)))) + 
     &                   g1R*(g2R*M4Lr*muu3*(g3R*g4L*(-1 + muu1**2 + muu2**2)*muu3*muu4 + g3L*g4R*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)) + 
     &                      g2L*muu2*(g3R*(g4L*muu3**2*muu4*(3*muu1**2 + muu3**2 - tred) + 
     &                            g4R*mGlr*(-(muu2**2*muu3**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + tred - muu4**2*tred)) + 
     &                         g3L*muu3*(g4L*mGlr*muu4*(muu1**2 + muu3**2 - tred) + 
     &                            g4R*(-(muu2**2*muu3**2) + muu3**2*muu4**2 + 2*muu1**2*(-1 + muu3**2 + muu4**2) + tred - muu4**2*tred)))))*
     &                 dconjg(FR(indj,i,tt,gen)))) - dconjg(ELc(indi,i,tt,gen))*
     &            ((g1R*muu1*muu3*(g2R*M4Lr*muu2*(2*g3L*g4R*muu3*muu4 + g3R*g4L*(-1 + muu3**2 + muu4**2)) + 
     &                    g2L*(g3L*muu3*(g4L*mGlr*(muu2**2 + muu4**2 - tred) + g4R*muu4*(-2 + 2*muu1**2 + muu2**2 + muu4**2 - tred)) + 
     &                       g3R*(g4R*mGlr*muu4*(-1 + muu1**2 + muu4**2 - tred) + 
     &                          g4L*(2*muu2**2*muu3**2 + muu1**2*muu4**2 + muu3**2*muu4**2 - tred - muu3**2*tred)))) + 
     &                 g1L*(g2L*M4Lr*muu3*(g3L*g4R*(-1 + muu1**2 + muu2**2)*muu3*muu4 + g3R*g4L*(muu2**2*muu3**2 + muu1**2*muu4**2 - tred)) + 
     &                    g2R*muu2*(g3L*(g4R*muu3**2*muu4*(3*muu1**2 + muu3**2 - tred) + 
     &                          g4L*mGlr*(-(muu2**2*muu3**2) + muu3**2*muu4**2 + muu1**2*(-1 + muu3**2 + muu4**2) + tred - muu4**2*tred)) + 
     &                       g3R*muu3*(g4R*mGlr*muu4*(muu1**2 + muu3**2 - tred) + 
     &                          g4L*(-(muu2**2*muu3**2) + muu3**2*muu4**2 + 2*muu1**2*(-1 + muu3**2 + muu4**2) + tred - muu4**2*tred)))))*
     &               dconjg(FL(indj,i,tt,gen)) + (g1R*muu1*muu3*
     &                  (g2L*muu2*(g3R*(g4R*mGlr*(-1 + muu3**2 + muu4**2) + g4L*muu4*(muu1**2 + 3*muu3**2 - tred)) + 
     &                       g3L*muu3*(2*g4L*mGlr*muu4 + g4R*(-2 + muu2**2 + 2*muu3**2 + muu4**2 - tred))) + 
     &                    g2R*M4Lr*(g3R*g4L*muu4*(-1 + muu1**2 + muu4**2 - tred) + g3L*g4R*muu3*(muu2**2 + muu4**2 - tred))) + 
     &                 g1L*(g2L*M4Lr*muu2*muu3*(g3R*g4L*muu4*(muu1**2 + muu3**2 - tred) + g3L*g4R*muu3*(-1 + muu2**2 + muu3**2 - tred)) + 
     &                    g2R*(g3L*(g4L*mGlr*muu4*(muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-1 + muu4**2 - 2*tred) + tred - muu3**2*tred - 
     &                             muu4**2*tred + tred**2) + g4R*muu3**2*
     &                           (muu1**2*(2*muu2**2 + muu4**2 - 2*tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))) + 
     &                       g3R*muu3*(g4L*muu4*(2*muu1**4 - muu2**2*muu3**2 + muu3**2*muu4**2 + muu1**2*(-2 + 2*muu4**2 - 3*tred) + tred - muu3**2*tred - 
     &                             muu4**2*tred + tred**2) + g4R*mGlr*
     &                           (muu1**2*(muu2**2 - tred) + muu3**2*(muu4**2 - tred) + tred*(1 - muu2**2 - muu4**2 + tred))))))*dconjg(FR(indj,i,tt,gen)))
     &           ))/2. + (xd2*(-(dconjg(ER(indi,i,tt,gen))*((g1L*muu1*
     &                    (g2L*M4Lr*(g3L*g4R*(-1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                         g3R*g4L*(muu2**2*muu3**2 - muu4**4 - muu3**2*tred + muu4**2*(1 - muu1**2 + muu3**2 + tred))) + 
     &                      g2R*muu2*(g3L*muu3*(g4L*mGlr*(-1 + muu3**2 - muu4**2) + g4R*muu4*(2 + muu1**2 - muu2**2 + muu3**2 - muu4**2)) + 
     &                         g3R*(g4R*mGlr*muu4*(1 + muu3**2 - muu4**2) + 
     &                            g4L*(2*muu3**4 + muu3**2*(-2 + muu2**2 - 2*muu4**2 - tred) + muu4**2*(-muu1**2 + tred))))) + 
     &                   g1R*(g2R*M4Lr*muu2*(g3L*g4R*(1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                         g3R*g4L*(muu3**4 + muu3**2*(-1 + muu2**2 - muu4**2 - tred) + muu4**2*(-muu1**2 + tred))) + 
     &                      g2L*(g3L*muu3*(g4L*mGlr*(-muu2**4 + muu1**2*(muu4**2 - tred) + muu2**2*(1 + muu1**2 - muu3**2 + tred)) + 
     &                            g4R*muu4*(2*muu1**4 + muu2**2*(-muu3**2 + tred) - muu1**2*(2 + 2*muu2**2 - muu4**2 + tred))) + 
     &                         g3R*(g4R*mGlr*muu4*(muu1**4 + muu2**2*(-muu3**2 + tred) - muu1**2*(1 + muu2**2 - muu4**2 + tred)) + 
     &                            g4L*(-(muu2**4*muu3**2) + muu2**2*muu3**2*(1 + 2*muu1**2 - muu3**2 + tred) + 
     &                               muu1**2*(-muu4**4 - 2*muu3**2*tred + muu4**2*(1 - muu1**2 + 2*muu3**2 + tred)))))))*dconjg(FL(indj,i,tt,gen)) + 
     &                (g1R*(g2R*M4Lr*(g3L*g4R*muu3*(muu4**2 + muu2**2*(muu3**2 - muu4**2) - tred) + 
     &                         g3R*g4L*muu4*((-1 + muu1**2)*muu3**2 - muu1**2*muu4**2 + tred)) + 
     &                      g2L*muu2*(g3R*(g4L*muu4*(-((-1 + muu2**2)*muu3**2) + muu1**2*(1 + 2*muu3**2 - muu4**2)) + 
     &                            g4R*mGlr*(-(muu2**2*muu3**2) + muu1**2*(-1 + muu3**2) + tred)) + 
     &                         g3L*muu3*(g4L*mGlr*(1 + muu1**2 - muu2**2)*muu4 + g4R*(-(muu2**2*muu3**2) + muu1**2*(-2 + 2*muu3**2 - muu4**2) + tred)))) + 
     &                   g1L*muu1*(g2L*M4Lr*muu2*(g3L*g4R*muu3*(-1 + muu3**2 - muu4**2) + g3R*g4L*muu4*(1 + muu3**2 - muu4**2)) + 
     &                      g2R*(g3L*muu3*(g4L*mGlr*(-1 + muu1**2 - muu2**2)*muu4 + 
     &                            g4R*((1 + muu1**2)*muu4**2 + muu2**2*(1 + muu3**2 - 2*muu4**2) - 2*tred)) - 
     &                         g3R*(g4L*muu4*((2 - 2*muu1**2 + muu2**2)*muu3**2 + muu1**2*muu4**2 - tred) + 
     &                            g4R*mGlr*(-(muu1**2*muu4**2) + muu2**2*(-1 + muu4**2) + tred)))))*dconjg(FR(indj,i,tt,gen)))) - 
     &           dconjg(ELc(indi,i,tt,gen))*((g1L*(g2L*M4Lr*(g3R*g4L*muu3*(muu4**2 + muu2**2*(muu3**2 - muu4**2) - tred) + 
     &                       g3L*g4R*muu4*((-1 + muu1**2)*muu3**2 - muu1**2*muu4**2 + tred)) + 
     &                    g2R*muu2*(g3L*(g4R*muu4*(-((-1 + muu2**2)*muu3**2) + muu1**2*(1 + 2*muu3**2 - muu4**2)) + 
     &                          g4L*mGlr*(-(muu2**2*muu3**2) + muu1**2*(-1 + muu3**2) + tred)) + 
     &                       g3R*muu3*(g4R*mGlr*(1 + muu1**2 - muu2**2)*muu4 + g4L*(-(muu2**2*muu3**2) + muu1**2*(-2 + 2*muu3**2 - muu4**2) + tred)))) + 
     &                 g1R*muu1*(g2R*M4Lr*muu2*(g3R*g4L*muu3*(-1 + muu3**2 - muu4**2) + g3L*g4R*muu4*(1 + muu3**2 - muu4**2)) + 
     &                    g2L*(g3R*muu3*(g4R*mGlr*(-1 + muu1**2 - muu2**2)*muu4 + 
     &                          g4L*((1 + muu1**2)*muu4**2 + muu2**2*(1 + muu3**2 - 2*muu4**2) - 2*tred)) - 
     &                       g3L*(g4R*muu4*((2 - 2*muu1**2 + muu2**2)*muu3**2 + muu1**2*muu4**2 - tred) + 
     &                          g4L*mGlr*(-(muu1**2*muu4**2) + muu2**2*(-1 + muu4**2) + tred)))))*dconjg(FL(indj,i,tt,gen)) + 
     &              (g1R*muu1*(g2R*M4Lr*(g3R*g4L*(-1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                       g3L*g4R*(muu2**2*muu3**2 - muu4**4 - muu3**2*tred + muu4**2*(1 - muu1**2 + muu3**2 + tred))) + 
     &                    g2L*muu2*(g3R*muu3*(g4R*mGlr*(-1 + muu3**2 - muu4**2) + g4L*muu4*(2 + muu1**2 - muu2**2 + muu3**2 - muu4**2)) + 
     &                       g3L*(g4L*mGlr*muu4*(1 + muu3**2 - muu4**2) + 
     &                          g4R*(2*muu3**4 + muu3**2*(-2 + muu2**2 - 2*muu4**2 - tred) + muu4**2*(-muu1**2 + tred))))) + 
     &                 g1L*(g2L*M4Lr*muu2*(g3R*g4L*(1 + muu1**2 - muu2**2)*muu3*muu4 + 
     &                       g3L*g4R*(muu3**4 + muu3**2*(-1 + muu2**2 - muu4**2 - tred) + muu4**2*(-muu1**2 + tred))) + 
     &                    g2R*(g3R*muu3*(g4R*mGlr*(-muu2**4 + muu1**2*(muu4**2 - tred) + muu2**2*(1 + muu1**2 - muu3**2 + tred)) + 
     &                          g4L*muu4*(2*muu1**4 + muu2**2*(-muu3**2 + tred) - muu1**2*(2 + 2*muu2**2 - muu4**2 + tred))) + 
     &                       g3L*(g4L*mGlr*muu4*(muu1**4 + muu2**2*(-muu3**2 + tred) - muu1**2*(1 + muu2**2 - muu4**2 + tred)) + 
     &                          g4R*(-(muu2**4*muu3**2) + muu2**2*muu3**2*(1 + 2*muu1**2 - muu3**2 + tred) + 
     &                             muu1**2*(-muu4**4 - 2*muu3**2*tred + muu4**2*(1 - muu1**2 + 2*muu3**2 + tred)))))))*dconjg(FR(indj,i,tt,gen)))))/2.))/
     &  (4.*Pi**2*(-tred + MSf(i,itypeT,gen)**2/sqrtS**2))

	enddo

	tempBoxt = tempBoxt + ChiChi2QQ_M2BoxGluino - tempBox2
	tempBox2 = ChiChi2QQ_M2BoxGluino

	enddo
	enddo
	
	
!	write(*,*)'ChiChi T-Channel Gluino Box (Neu) = ', ChiChi2QQ_M2BoxGluino - tempBox
!	write(*,*)'ChiChi T-Channel Gluino Box (Neu) mit Z = ', tempBoxZ
!	write(*,*)'ChiChi T-Channel Gluino Box (Neu) mit H = ', tempBoxH
!	write(*,*)'ChiChi T-Channel Gluino Box (Neu) mit t = ', tempBoxt
!	write(*,*)'ChiChi T-Channel Gluino Box (Neu) mit u = ', tempBoxu

	
	end


