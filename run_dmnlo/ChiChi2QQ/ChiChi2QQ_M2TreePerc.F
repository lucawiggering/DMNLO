************************************************

      function ChiChi2QQ_M2TreePerc(ichannel)

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"

	integer i,j,imin,imax,jmin,jmax,itypeU,itypeT,att,antitt
	complex*16 DenSmV, DenSmH(6)
	double precision DenSmV1, DenSmH1(6), M2test, sign

	double precision ChiChi2QQ_M2TreePerc

	ChiChi2QQ_M2TreePerc = 0d0
	tempTree = 0d0


! itypeT,U follows the type of the sfermion exchanged in t,u-channel interactions
! att gives the type of the antiquark final state 
! antitt is the opposite type to the final state quark type

	if (icha.eq.0) then
		att   = tt
		itypeU = tt
		itypeT = tt
		antitt = tt
	else if (icha.eq.1) then
! tt is defined normally as up-type (quark) type in chargino-neutralino co-annihilations		
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

! widths of Z-boson and Higgs included in propagators	

	if (icha.eq.1) then
		DenSmV = (sqrtS**2 - MW2) + cI*MW*wW
		DenSmV1 = DenSmV*dconjg(DenSmV)/(sqrtS**2 - MW2)
	else
		DenSmV = (sqrtS**2 - MZ2) + cI*MZ*wZ
		DenSmV1 = DenSmV*dconjg(DenSmV)/(sqrtS**2 - MZ2)
	endif

	do i = imin, imax
	  DenSmH(i) = (sqrtS**2 - MHk(i)**2) + cI*MHk(i)*wHk(i)
	  DenSmH1(i) = DenSmH(i)*dconjg(DenSmH(i))/(sqrtS**2 - MHk(i)**2)
	enddo

! Photon contributions only in cha-cha case
! cha cha -> b B : no t-channel
! cha cha -> t T : no u-channel
! Use couplings E/F, G/H to switch on/off t-/u-channels	
	
! *********** Vector boson contr.	


        if(ichannel.eq.11) then

! Photon exchange in Chargino annihilation (m3=m4)
	if (icha.eq.2) then
		
! Photon
	ChiChi2QQ_M2TreePerc = ChiChi2QQ_M2TreePerc + 8*ApL**2*BpL**2*NC*(1 + 2*muu3**4 + 2*muu1*(muu2 + 2*muu2*muu3**2) + 
     &    muu1**2*(-1 + 2*muu2**2 - 2*ured) + 2*ured - 4*muu3**2*ured + 2*ured**2 - muu2**2*(1 + 2*ured))	
	
cc	write(*,*)'P Tree mit P = ',ChiChi2QQ_M2TreePerc - tempTree
	tempTree = ChiChi2QQ_M2TreePerc
		
! Photon-Z interference
	ChiChi2QQ_M2TreePerc = ChiChi2QQ_M2TreePerc + (8*ApL*BpL*NC*sqrtS**2*((AR*BL + AL*BR)*(muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*tred) + (muu3**2 - tred)**2) + 
     &      (AL*BL + AR*BR)*(muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*ured) + (muu3**2 - ured)**2)))/(dreal(dconjg(DenSmV)))
	
cc	write(*,*)'P Tree mit Z = ',ChiChi2QQ_M2TreePerc - tempTree
	tempTree = ChiChi2QQ_M2TreePerc
	
	endif


! Massive vector-boson (Z/W)
	ChiChi2QQ_M2TreePerc = ChiChi2QQ_M2TreePerc + (-4*NC*sqrtS**4*((AR*muu1*muu2*((BL**2 + BR**2)*(-1 + muu3**2) - 8*BL*BR*muu3*muu4 + (BL**2 + BR**2)*muu4**2) - 
     &         AL*(-2*BL*BR*(-1 + muu1**2 + muu2**2)*muu3*muu4 + BR**2*(muu1**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred) + 
     &            BL**2*(muu2**2 + muu3**2 - ured)*(muu1**2 + muu4**2 - ured)))*dconjg(AL) + 
     &      (AL*muu1*muu2*((BL**2 + BR**2)*(-1 + muu3**2) - 8*BL*BR*muu3*muu4 + (BL**2 + BR**2)*muu4**2) - 
     &         AR*(-2*BL*BR*(-1 + muu1**2 + muu2**2)*muu3*muu4 + BL**2*(muu1**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred) + 
     &            BR**2*(muu2**2 + muu3**2 - ured)*(muu1**2 + muu4**2 - ured)))*dconjg(AR)))/(DenSmV*dconjg(DenSmV))
	
cc	write(*,*)'Z Tree mit Z = ',ChiChi2QQ_M2TreePerc - tempTree
	tempTree = ChiChi2QQ_M2TreePerc

        
        elseif(ichannel.eq.22) then


! Higgs contr.
	do i = imin, imax
	  do j = jmin, jmax

	ChiChi2QQ_M2TreePerc = ChiChi2QQ_M2TreePerc + (NC*sqrtS**4*(CL(i)*((-1 + muu1**2 + muu2**2)*dconjg(CL(j)) + 2*muu1*muu2*dconjg(CR(j))) + 
     &      CR(i)*(2*muu1*muu2*dconjg(CL(j)) + (-1 + muu1**2 + muu2**2)*dconjg(CR(j))))*
     &    (((-1 + muu3**2 + muu4**2)*dconjg(DL(j)) + 2*muu3*muu4*dconjg(DR(j)))*DL(i) + 
     &      (2*muu3*muu4*dconjg(DL(j)) + (-1 + muu3**2 + muu4**2)*dconjg(DR(j)))*DR(i)))/(dconjg(DenSmH(j))*DenSmH(i))

	  enddo
	enddo
	
cc	write(*,*)'H Tree mit H = ',ChiChi2QQ_M2TreePerc - tempTree
	tempTree = ChiChi2QQ_M2TreePerc
	
	elseif(ichannel.eq.12) then

! Vector-Higgs contr.

! Photon exchange in Chargino annihilation (m3=m4)
	if (icha.eq.2) then

	do i = imin,imax
		
	ChiChi2QQ_M2TreePerc = ChiChi2QQ_M2TreePerc + (4*ApL*BpL*muu1*muu3*NC*sqrtS**2*(tred - ured)*
     &    (CL(i)*(DL(i) + DR(i)) + CR(i)*(DL(i) + DR(i)) + (dconjg(CL(i)) + dconjg(CR(i)))*(dconjg(DL(i)) +dconjg(DR(i)))))/(dreal(dconjg(DenSmH(i))))
		
	enddo	
	
cc	write(*,*)'P Tree mit H = ',ChiChi2QQ_M2TreePerc - tempTree
	tempTree = ChiChi2QQ_M2TreePerc
	
	endif
	
! Z/W-Higgs contr.	
	do i = imin,imax
	
	ChiChi2QQ_M2TreePerc = ChiChi2QQ_M2TreePerc + (-2*NC*sqrtS**4*(-(muu4*(muu1*(muu2**2 + muu3**2 - ured)*(AR*dconjg(CL(i)) + AL*dconjg(CR(i))) - 
     &           muu2*(muu1**2 + muu3**2 - tred)*(AL*dconjg(CL(i)) + AR*dconjg(CR(i))))*(BR*dconjg(DL(i)) + BL*dconjg(DR(i)))) + 
     &      muu3*(muu1*(muu2**2 + muu4**2 - tred)*(AR*dconjg(CL(i)) + AL*dconjg(CR(i))) - 
     &         muu2*(muu1**2 + muu4**2 - ured)*(AL*dconjg(CL(i)) + AR*dconjg(CR(i))))*(BL*dconjg(DL(i)) + BR*dconjg(DR(i))) - 
     &      muu4*(muu1*(muu2**2 + muu3**2 - ured)*(CR(i)*dconjg(AL) + CL(i)*dconjg(AR)) - 
     &         muu2*(muu1**2 + muu3**2 - tred)*(CL(i)*dconjg(AL) + CR(i)*dconjg(AR)))*(dconjg(BR)*DL(i) + dconjg(BL)*DR(i)) + 
     &      muu3*(muu1*(muu2**2 + muu4**2 - tred)*(CR(i)*dconjg(AL) + CL(i)*dconjg(AR)) - 
     &         muu2*(muu1**2 + muu4**2 - ured)*(CL(i)*dconjg(AL) + CR(i)*dconjg(AR)))*(dconjg(BL)*DL(i) + dconjg(BR)*DR(i))))/
     &  dreal((DenSmV*dconjg(DenSmH(i))))
	 
	enddo
	
cc	write(*,*)'Z Tree mit H = ',ChiChi2QQ_M2TreePerc - tempTree
	tempTree = ChiChi2QQ_M2TreePerc

        elseif(ichannel.eq.33) then

! T-channel 
! (zero for up-type quarks in chargino-chargino annihilation through zero couplings)
	do i = 1, 2
	  do j = 1, 2

	ChiChi2QQ_M2TreePerc = ChiChi2QQ_M2TreePerc - ((NC*(((muu1**2 + muu3**2 - tred)*dconjg(ELc(indi,j,tt,gen)) + 2*muu1*muu3*dconjg(ER(indi,j,tt,gen)))*
     &         ELc(indi,i,tt,gen) + (2*muu1*muu3*dconjg(ELc(indi,j,tt,gen)) + 
     &           (muu1**2 + muu3**2 - tred)*dconjg(ER(indi,j,tt,gen)))*ER(indi,i,tt,gen))*
     &      (-(((muu2**2 + muu4**2 - tred)*dconjg(FL(indj,j,tt,gen)) + 2*muu2*muu4*dconjg(FR(indj,j,tt,gen)))*
     &           FL(indj,i,tt,gen)) - (2*muu2*muu4*dconjg(FL(indj,j,tt,gen)) + 
     &           (muu2**2 + muu4**2 - tred)*dconjg(FR(indj,j,tt,gen)))*FR(indj,i,tt,gen)))/
     &    ((-tred + MSf(i,itypeT,gen)**2/sqrtS**2)*(-tred + MSf(j,itypeT,gen)**2/sqrtS**2)))
	
	  enddo
	enddo

cc	write(*,*)'t Tree mit t = ',ChiChi2QQ_M2TreePerc - tempTree
	tempTree = ChiChi2QQ_M2TreePerc


! TU-channel contr. 
! relative factor (-1) from U-ch.spinor order
	do i = 1, 2
	  do j = 1, 2
	  
	ChiChi2QQ_M2TreePerc = ChiChi2QQ_M2TreePerc - dreal((NC*(ELc(indi,i,tt,gen)*((2*muu1*muu4*dconjg(GR(indi,j,tt,gen))*
     &               ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,j,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,j,tt,gen))) + 
     &              dconjg(GL(indi,j,tt,gen))*((-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - 
     &                    (muu2**2 + muu3**2 + muu4**2 - tred)*tred + muu1**2*(1 + 2*muu2**2 - tred - ured) - 
     &                    (muu2**2 + muu3**2 + muu4**2)*ured + ured**2)*dconjg(HL(indj,j,tt,gen)) + 
     &                 2*muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,j,tt,gen))))*FL(indj,i,tt,gen) + 
     &           2*(muu4*dconjg(GL(indi,j,tt,gen))*(muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,j,tt,gen)) - 
     &                 (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,j,tt,gen))) + 
     &              muu1*dconjg(GR(indi,j,tt,gen))*(-(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,j,tt,gen))) + 
     &                 muu3*(muu2**2 + muu4**2 - tred)*dconjg(HR(indj,j,tt,gen))))*FR(indj,i,tt,gen)) + 
     &        ER(indi,i,tt,gen)*(2*(muu1*dconjg(GL(indi,j,tt,gen))*
     &               (muu3*(muu2**2 + muu4**2 - tred)*dconjg(HL(indj,j,tt,gen)) - 
     &                 muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,j,tt,gen))) + 
     &              muu4*dconjg(GR(indi,j,tt,gen))*(-((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,j,tt,gen))) + 
     &                 muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,j,tt,gen))))*FL(indj,i,tt,gen) + 
     &           (2*muu1*muu4*dconjg(GL(indi,j,tt,gen))*
     &               (2*muu2*muu3*dconjg(HL(indj,j,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,j,tt,gen))) + 
     &              dconjg(GR(indi,j,tt,gen))*(2*muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,j,tt,gen)) + 
     &                 (-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - (muu2**2 + muu3**2 + muu4**2 - tred)*tred + 
     &                    muu1**2*(1 + 2*muu2**2 - tred - ured) - (muu2**2 + muu3**2 + muu4**2)*ured + ured**2)*
     &                  dconjg(HR(indj,j,tt,gen))))*FR(indj,i,tt,gen))))/
     &    ((-tred + MSf(i,itypeT,gen)**2/sqrtS**2)*(-ured + MSf(j,itypeU,gen)**2/sqrtS**2)))

	  enddo
	enddo

cc	write(*,*)'t Tree mit u = ',ChiChi2QQ_M2TreePerc - tempTree
	tempTree = ChiChi2QQ_M2TreePerc

! U-channel 
! (zero for down-type quarks in chargino-chargino annihilation through zero couplings)
	do i = 1, 2
	  do j = 1, 2

	ChiChi2QQ_M2TreePerc = ChiChi2QQ_M2TreePerc -((NC*(((muu1**2 + muu4**2 - ured)*dconjg(GL(indi,j,tt,gen)) + 2*muu1*muu4*dconjg(GR(indi,j,tt,gen)))*GL(indi,i,tt,gen) + 
     &        (2*muu1*muu4*dconjg(GL(indi,j,tt,gen)) + (muu1**2 + muu4**2 - ured)*dconjg(GR(indi,j,tt,gen)))*GR(indi,i,tt,gen))*
     &      (-(((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,j,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,j,tt,gen)))*
     &           HL(indj,i,tt,gen)) - (2*muu2*muu3*dconjg(HL(indj,j,tt,gen)) + 
     &           (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,j,tt,gen)))*HR(indj,i,tt,gen)))/
     &    ((-ured + MSf(i,itypeU,gen)**2/sqrtS**2)*(-ured + MSf(j,itypeU,gen)**2/sqrtS**2)))
	  
	  enddo
	enddo
	
cc	write(*,*)'u Tree mit u = ',ChiChi2QQ_M2TreePerc - tempTree
	tempTree = ChiChi2QQ_M2TreePerc


        elseif(ichannel.eq.13) then
	
	
! TS-channel photon (for cha cha -> t tbar)
! relative factor (-1)(-1) from vector propagator, T-ch.spinor order
	if (icha.eq.2) then
	
	do i = 1, 2
	
	ChiChi2QQ_M2TreePerc = ChiChi2QQ_M2TreePerc + (-4*ApL*BpL*NC*(ELc(indi,i,tt,gen)*(2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*FL(indj,i,tt,gen) + 
     &         (muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)*FR(indj,i,tt,gen)) + 
     &      ER(indi,i,tt,gen)*((muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)*
     &          FL(indj,i,tt,gen) + 2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*FR(indj,i,tt,gen))))/
     &  (-tred + MSf(i,itypeT,gen)**2/sqrtS**2)
	 
	enddo
	
cc	write(*,*)'t Tree mit P = ',ChiChi2QQ_M2TreePerc - tempTree
	tempTree = ChiChi2QQ_M2TreePerc
	
	endif	

	
! TS-channel vector-boson (Z/W)
! relative factor (-1)(-1) from vector propagator, T-ch.spinor order
	do i = 1,2
	
	ChiChi2QQ_M2TreePerc = ChiChi2QQ_M2TreePerc + (-4*NC*sqrtS**2*(dconjg(AL)*(BR*((muu1**2 + muu3**2 - tred)*ELc(indi,i,tt,gen) + 2*muu1*muu3*ER(indi,i,tt,gen))*
     &          (2*muu2*muu4*FL(indj,i,tt,gen) + (muu2**2 + muu4**2 - tred)*FR(indj,i,tt,gen)) + 
     &         BL*(muu3*ELc(indi,i,tt,gen)*(muu2*(muu1**2 + muu4**2 - ured)*FL(indj,i,tt,gen) - 
     &               (-1 + muu1**2 + muu2**2)*muu4*FR(indj,i,tt,gen)) + 
     &            muu1*ER(indi,i,tt,gen)*(-(muu2*(-1 + muu3**2 + muu4**2)*FL(indj,i,tt,gen)) + 
     &               muu4*(muu2**2 + muu3**2 - ured)*FR(indj,i,tt,gen)))) + 
     &      dconjg(AR)*(BL*(2*muu1*muu3*ELc(indi,i,tt,gen) + (muu1**2 + muu3**2 - tred)*ER(indi,i,tt,gen))*
     &          ((muu2**2 + muu4**2 - tred)*FL(indj,i,tt,gen) + 2*muu2*muu4*FR(indj,i,tt,gen)) + 
     &         BR*(muu1*ELc(indi,i,tt,gen)*(muu4*(muu2**2 + muu3**2 - ured)*FL(indj,i,tt,gen) - 
     &               muu2*(-1 + muu3**2 + muu4**2)*FR(indj,i,tt,gen)) + 
     &            muu3*ER(indi,i,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu4*FL(indj,i,tt,gen)) + 
     &               muu2*(muu1**2 + muu4**2 - ured)*FR(indj,i,tt,gen))))))/(DenSmV1*(-tred + MSf(i,itypeT,gen)**2/sqrtS**2))

	enddo
	
cc	write(*,*)'t Tree mit Z = ',ChiChi2QQ_M2TreePerc - tempTree
	tempTree = ChiChi2QQ_M2TreePerc


! US-channel photon (for cha cha -> b bbar)
! relative factor (-1)(-1)(-1) from vector propagator, S-ch.spinor order, Charge conjugated vector coupling
	if (icha.eq.2) then
	
	do i = 1, 2
	
	ChiChi2QQ_M2TreePerc = ChiChi2QQ_M2TreePerc + (4*ApL*BpL*NC*(GR(indi,i,tt,gen)*((muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + 
     &            ured**2)*HL(indj,i,tt,gen) + 2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*HR(indj,i,tt,gen)) + 
     &      GL(indi,i,tt,gen)*(2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*HL(indj,i,tt,gen) + 
     &         (muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + ured**2)*HR(indj,i,tt,gen))))/
     &  (-ured + MSf(i,itypeU,gen)**2/sqrtS**2)
	
	enddo
	
cc	write(*,*)'u Tree mit P = ',ChiChi2QQ_M2TreePerc - tempTree
	tempTree = ChiChi2QQ_M2TreePerc
	
	endif
	
	
! US-channel vector-boson (Z/W)
! relative factor (-1)(-1)(-1) from vector propagator, S-ch.spinor order, Charge conjugated vector coupling
	do i = 1,2

	ChiChi2QQ_M2TreePerc = ChiChi2QQ_M2TreePerc + dreal(4*NC*sqrtS**2*(dconjg(AR)*(BR*(2*muu1*muu4*GL(indi,i,tt,gen) + (muu1**2 + muu4**2 - ured)*GR(indi,i,tt,gen))*
     &          ((muu2**2 + muu3**2 - ured)*HL(indj,i,tt,gen) + 2*muu2*muu3*HR(indj,i,tt,gen)) + 
     &         BL*(muu1*GL(indi,i,tt,gen)*(muu3*(muu2**2 + muu4**2 - tred)*HL(indj,i,tt,gen) - 
     &               muu2*(-1 + muu3**2 + muu4**2)*HR(indj,i,tt,gen)) + 
     &            muu4*GR(indi,i,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu3*HL(indj,i,tt,gen)) + 
     &               muu2*(muu1**2 + muu3**2 - tred)*HR(indj,i,tt,gen)))) + 
     &      dconjg(AL)*(BL*((muu1**2 + muu4**2 - ured)*GL(indi,i,tt,gen) + 2*muu1*muu4*GR(indi,i,tt,gen))*
     &          (2*muu2*muu3*HL(indj,i,tt,gen) + (muu2**2 + muu3**2 - ured)*HR(indj,i,tt,gen)) + 
     &         BR*(muu4*GL(indi,i,tt,gen)*(muu2*(muu1**2 + muu3**2 - tred)*HL(indj,i,tt,gen) - 
     &               (-1 + muu1**2 + muu2**2)*muu3*HR(indj,i,tt,gen)) + 
     &            muu1*GR(indi,i,tt,gen)*(-(muu2*(-1 + muu3**2 + muu4**2)*HL(indj,i,tt,gen)) + 
     &               muu3*(muu2**2 + muu4**2 - tred)*HR(indj,i,tt,gen))))))/(DenSmV1*(-ured + MSf(i,itypeU,gen)**2/sqrtS**2))
	
	enddo
	
cc	write(*,*)'u Tree mit Z = ',ChiChi2QQ_M2TreePerc - tempTree
	tempTree = ChiChi2QQ_M2TreePerc

      
        elseif(ichannel.eq.23) then


! TS-channel Higgs boson
! relative factor (-1) from T-ch.spinor order
	do i=1,2
	  do j = jmin, jmax
	
	ChiChi2QQ_M2TreePerc = ChiChi2QQ_M2TreePerc + dreal(NC*sqrtS**2*(dconjg(CL(j))*(dconjg(DL(j))*(2*muu1*muu3*ER(indi,i,tt,gen)*
     &             ((muu2**2 + muu4**2 - tred)*FL(indj,i,tt,gen) + 2*muu2*muu4*FR(indj,i,tt,gen)) + 
     &            ELc(indi,i,tt,gen)*((1 - muu3**2 - muu4**2 - muu3**2*tred - muu4**2*tred + tred**2 + muu3**2*ured + 
     &                  muu4**2*ured - ured**2 + muu2**2*(-1 + 2*muu3**2 - tred + ured) + muu1**2*(-1 + 2*muu4**2 - tred + ured))
     &                 *FL(indj,i,tt,gen) + 2*muu2*muu4*(muu1**2 + muu3**2 - tred)*FR(indj,i,tt,gen))) + 
     &         2*dconjg(DR(j))*(muu1*ER(indi,i,tt,gen)*
     &             (-(muu4*(muu2**2 + muu3**2 - ured)*FL(indj,i,tt,gen)) + muu2*(-1 + muu3**2 + muu4**2)*FR(indj,i,tt,gen)) + 
     &            muu3*ELc(indi,i,tt,gen)*((-1 + muu1**2 + muu2**2)*muu4*FL(indj,i,tt,gen) - 
     &               muu2*(muu1**2 + muu4**2 - ured)*FR(indj,i,tt,gen)))) + 
     &      dconjg(CR(j))*(2*dconjg(DL(j))*(muu3*ER(indi,i,tt,gen)*
     &             (-(muu2*(muu1**2 + muu4**2 - ured)*FL(indj,i,tt,gen)) + (-1 + muu1**2 + muu2**2)*muu4*FR(indj,i,tt,gen)) + 
     &            muu1*ELc(indi,i,tt,gen)*(muu2*(-1 + muu3**2 + muu4**2)*FL(indj,i,tt,gen) - 
     &               muu4*(muu2**2 + muu3**2 - ured)*FR(indj,i,tt,gen))) + 
     &         dconjg(DR(j))*(2*muu1*muu3*ELc(indi,i,tt,gen)*
     &             (2*muu2*muu4*FL(indj,i,tt,gen) + (muu2**2 + muu4**2 - tred)*FR(indj,i,tt,gen)) + 
     &            ER(indi,i,tt,gen)*(2*muu2*muu4*(muu1**2 + muu3**2 - tred)*FL(indj,i,tt,gen) + 
     &               (1 - muu3**2 - muu4**2 - muu3**2*tred - muu4**2*tred + tred**2 + muu3**2*ured + muu4**2*ured - ured**2 + 
     &                  muu2**2*(-1 + 2*muu3**2 - tred + ured) + muu1**2*(-1 + 2*muu4**2 - tred + ured))*FR(indj,i,tt,gen))))))/
     &  (DenSmH1(j)*(-tred + MSf(i,itypeT,gen)**2/sqrtS**2))

	  enddo
	enddo
	
cc	write(*,*)'t Tree mit H = ',ChiChi2QQ_M2TreePerc - tempTree
	tempTree = ChiChi2QQ_M2TreePerc


! US-channel Higgs boson
! relative factor (-1) from S-ch.spinor order
	do i=1,2
	  do j = jmin, jmax
	  
	ChiChi2QQ_M2TreePerc = ChiChi2QQ_M2TreePerc + dreal(NC*sqrtS**2*(dconjg(CL(j))*(2*dconjg(DR(j))*
     &          (muu1*GR(indi,i,tt,gen)*(-(muu3*(muu2**2 + muu4**2 - tred)*HL(indj,i,tt,gen)) + 
     &               muu2*(-1 + muu3**2 + muu4**2)*HR(indj,i,tt,gen)) + 
     &            muu4*GL(indi,i,tt,gen)*((-1 + muu1**2 + muu2**2)*muu3*HL(indj,i,tt,gen) - 
     &               muu2*(muu1**2 + muu3**2 - tred)*HR(indj,i,tt,gen))) + 
     &         dconjg(DL(j))*(2*muu1*muu4*GR(indi,i,tt,gen)*
     &             ((muu2**2 + muu3**2 - ured)*HL(indj,i,tt,gen) + 2*muu2*muu3*HR(indj,i,tt,gen)) + 
     &            GL(indi,i,tt,gen)*(((-1 + muu3**2 + muu4**2 - tred)*(-1 + tred) + muu1**2*(-1 + 2*muu3**2 + tred - ured) + 
     &                  muu2**2*(-1 + 2*muu4**2 + tred - ured) - (muu3**2 + muu4**2)*ured + ured**2)*HL(indj,i,tt,gen) + 
     &               2*muu2*muu3*(muu1**2 + muu4**2 - ured)*HR(indj,i,tt,gen)))) + 
     &      dconjg(CR(j))*(2*dconjg(DL(j))*(muu4*GR(indi,i,tt,gen)*
     &             (-(muu2*(muu1**2 + muu3**2 - tred)*HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu2**2)*muu3*HR(indj,i,tt,gen)) + 
     &            muu1*GL(indi,i,tt,gen)*(muu2*(-1 + muu3**2 + muu4**2)*HL(indj,i,tt,gen) - 
     &               muu3*(muu2**2 + muu4**2 - tred)*HR(indj,i,tt,gen))) + 
     &         dconjg(DR(j))*(2*muu1*muu4*GL(indi,i,tt,gen)*
     &             (2*muu2*muu3*HL(indj,i,tt,gen) + (muu2**2 + muu3**2 - ured)*HR(indj,i,tt,gen)) + 
     &            GR(indi,i,tt,gen)*(2*muu2*muu3*(muu1**2 + muu4**2 - ured)*HL(indj,i,tt,gen) + 
     &               ((-1 + muu3**2 + muu4**2 - tred)*(-1 + tred) + muu1**2*(-1 + 2*muu3**2 + tred - ured) + 
     &                  muu2**2*(-1 + 2*muu4**2 + tred - ured) - (muu3**2 + muu4**2)*ured + ured**2)*HR(indj,i,tt,gen))))))/
     &  (dconjg(DenSmH(j))*(-ured + MSf(i,itypeU,gen)**2/sqrtS**2))

	  enddo
	enddo
	
cc	write(*,*)'u Tree mit H = ',ChiChi2QQ_M2TreePerc - tempTree
	tempTree = ChiChi2QQ_M2TreePerc

        endif  !end of ichannel

	end
