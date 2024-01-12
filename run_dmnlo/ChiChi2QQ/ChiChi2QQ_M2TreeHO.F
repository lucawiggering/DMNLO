************************************************

	function ChiChi2QQ_M2TreeHO()

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"

	integer i,j,imin,imax,jmin,jmax,itypeU,itypeT,att,antitt
	complex*16 DenSmV, DenSmH(6)
	double precision DenSmV1, DenSmH1(6), M2test, sign

	double precision ChiChi2QQ_M2TreeHO

	ChiChi2QQ_M2TreeHO = 0d0


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
	
! *********** Vector boson contr.				c			


! PS: No higher order corrections to these processes!
! Photon exchange in Chargino annihilation (m3=m4)
!	if (icha.eq.2) then
!		
! Photon
!	ChiChi2QQ_M2TreeHO = ChiChi2QQ_M2TreeHO + 8*ApL**2*BpL**2*NC*(1 + 2*muu3**4 + 2*muu1*(muu2 + 2*muu2*muu3**2) + 
!     &    muu1**2*(-1 + 2*muu2**2 - 2*ured) + 2*ured - 4*muu3**2*ured + 2*ured**2 - muu2**2*(1 + 2*ured))	
!
! Photon-Z interference
!	ChiChi2QQ_M2TreeHO = ChiChi2QQ_M2TreeHO + (8*ApL*BpL*NC*sqrtS**2*((AR*BL + AL*BR)*(muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*tred) + (muu3**2 - tred)**2) + 
!     &      (AL*BL + AR*BR)*(muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*ured) + (muu3**2 - ured)**2)))/DenSmV
!
!	endif
!
! Massive vector-boson (Z/W)
!	ChiChi2QQ_M2TreeHO = ChiChi2QQ_M2TreeHO + (-4*NC*sqrtS**4*((AR*muu1*muu2*((BL**2 + BR**2)*(-1 + muu3**2) - 8*BL*BR*muu3*muu4 + (BL**2 + BR**2)*muu4**2) - 
!     &         AL*(-2*BL*BR*(-1 + muu1**2 + muu2**2)*muu3*muu4 + BR**2*(muu1**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred) + 
!     &            BL**2*(muu2**2 + muu3**2 - ured)*(muu1**2 + muu4**2 - ured)))*dconjg(AL) + 
!     &      (AL*muu1*muu2*((BL**2 + BR**2)*(-1 + muu3**2) - 8*BL*BR*muu3*muu4 + (BL**2 + BR**2)*muu4**2) - 
!     &         AR*(-2*BL*BR*(-1 + muu1**2 + muu2**2)*muu3*muu4 + BL**2*(muu1**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred) + 
!     &            BR**2*(muu2**2 + muu3**2 - ured)*(muu1**2 + muu4**2 - ured)))*dconjg(AR)))/(DenSmV*dconjg(DenSmV))


! Higgs contr.

	do i = imin, imax
	  do j = jmin, jmax

	
	ChiChi2QQ_M2TreeHO = ChiChi2QQ_M2TreeHO + (NC*sqrtS**4*(CL(i)*((-1 + muu1**2 + muu2**2)*dconjg(CL(j)) + 2*muu1*muu2*dconjg(CR(j))) + 
     &      CR(i)*(2*muu1*muu2*dconjg(CL(j)) + (-1 + muu1**2 + muu2**2)*dconjg(CR(j))))*
     &    (((-1 + muu3**2 + muu4**2)*dconjg(DLho(j)) + 2*muu3*muu4*dconjg(DRho(j)))*DLho(i) + 
     &      (2*muu3*muu4*dconjg(DLho(j)) + (-1 + muu3**2 + muu4**2)*dconjg(DRho(j)))*DRho(i)))/(dconjg(DenSmH(j))*DenSmH(i))

	
	  enddo

! Vector-Higgs contr.

! Photon exchange in Chargino annihilation (m3=m4)
	if (icha.eq.2) then
	
	ChiChi2QQ_M2TreeHO = ChiChi2QQ_M2TreeHO + (4*ApL*BpL*muu1*muu3*NC*sqrtS**2*(tred - ured)*
     &    (CL(i)*(DLho(i) + DRho(i)) + CR(i)*(DLho(i) + DRho(i)) + (CL(j) + CR(j))*(DLho(j) + DRho(j))))/(dreal(dconjg(DenSmH(i))))
		
	endif
	
	ChiChi2QQ_M2TreeHO = ChiChi2QQ_M2TreeHO + (-2*NC*sqrtS**4*(-(muu4*(muu1*(muu2**2 + muu3**2 - ured)*(AR*dconjg(CL(i)) + AL*dconjg(CR(i))) - 
     &           muu2*(muu1**2 + muu3**2 - tred)*(AL*dconjg(CL(i)) + AR*dconjg(CR(i))))*(BR*dconjg(DLho(i)) + BL*dconjg(DRho(i)))) + 
     &      muu3*(muu1*(muu2**2 + muu4**2 - tred)*(AR*dconjg(CL(i)) + AL*dconjg(CR(i))) - 
     &         muu2*(muu1**2 + muu4**2 - ured)*(AL*dconjg(CL(i)) + AR*dconjg(CR(i))))*(BL*dconjg(DLho(i)) + BR*dconjg(DRho(i))) - 
     &      muu4*(muu1*(muu2**2 + muu3**2 - ured)*(CR(i)*dconjg(AL) + CL(i)*dconjg(AR)) - 
     &         muu2*(muu1**2 + muu3**2 - tred)*(CL(i)*dconjg(AL) + CR(i)*dconjg(AR)))*(dconjg(BR)*DLho(i) + dconjg(BL)*DRho(i)) + 
     &      muu3*(muu1*(muu2**2 + muu4**2 - tred)*(CR(i)*dconjg(AL) + CL(i)*dconjg(AR)) - 
     &         muu2*(muu1**2 + muu4**2 - ured)*(CL(i)*dconjg(AL) + CR(i)*dconjg(AR)))*(dconjg(BL)*DLho(i) + dconjg(BR)*DRho(i))))/
     &  dreal((DenSmV*dconjg(DenSmH(i))))
	 
	enddo

! T-channel 
! (zero for up-type quarks in chargino-chargino annihilation through zero couplings)
		
	do i = 1, 2
	  do j = 1, 2

	ChiChi2QQ_M2TreeHO = ChiChi2QQ_M2TreeHO - ((NC*(((muu1**2 + muu3**2 - tred)*dconjg(ELcho(indi,j,tt,gen)) + 2*muu1*muu3*dconjg(ERho(indi,j,tt,gen)))*
     &         ELcho(indi,i,tt,gen) + (2*muu1*muu3*dconjg(ELcho(indi,j,tt,gen)) + 
     &           (muu1**2 + muu3**2 - tred)*dconjg(ERho(indi,j,tt,gen)))*ERho(indi,i,tt,gen))*
     &      (-(((muu2**2 + muu4**2 - tred)*dconjg(FLho(indj,j,tt,gen)) + 2*muu2*muu4*dconjg(FRho(indj,j,tt,gen)))*
     &           FLho(indj,i,tt,gen)) - (2*muu2*muu4*dconjg(FLho(indj,j,tt,gen)) + 
     &           (muu2**2 + muu4**2 - tred)*dconjg(FRho(indj,j,tt,gen)))*FRho(indj,i,tt,gen)))/
     &    ((-tred + MSf(i,itypeT,gen)**2/sqrtS**2)*(-tred + MSf(j,itypeT,gen)**2/sqrtS**2)))
	
	
	  enddo
	enddo




	do i = 1, 2
	  do j = 1, 2

! TU-channel contr. 
! relative factor (-1) from U-ch.spinor order

	ChiChi2QQ_M2TreeHO = ChiChi2QQ_M2TreeHO - dreal((NC*(ELcho(indi,i,tt,gen)*((2*muu1*muu4*dconjg(GRho(indi,j,tt,gen))*
     &               ((muu2**2 + muu3**2 - ured)*dconjg(HLho(indj,j,tt,gen)) + 2*muu2*muu3*dconjg(HRho(indj,j,tt,gen))) + 
     &              dconjg(GLho(indi,j,tt,gen))*((-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - 
     &                    (muu2**2 + muu3**2 + muu4**2 - tred)*tred + muu1**2*(1 + 2*muu2**2 - tred - ured) - 
     &                    (muu2**2 + muu3**2 + muu4**2)*ured + ured**2)*dconjg(HLho(indj,j,tt,gen)) + 
     &                 2*muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HRho(indj,j,tt,gen))))*FLho(indj,i,tt,gen) + 
     &           2*(muu4*dconjg(GLho(indi,j,tt,gen))*(muu2*(muu1**2 + muu3**2 - tred)*dconjg(HLho(indj,j,tt,gen)) - 
     &                 (-1 + muu1**2 + muu2**2)*muu3*dconjg(HRho(indj,j,tt,gen))) + 
     &              muu1*dconjg(GRho(indi,j,tt,gen))*(-(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HLho(indj,j,tt,gen))) + 
     &                 muu3*(muu2**2 + muu4**2 - tred)*dconjg(HRho(indj,j,tt,gen))))*FRho(indj,i,tt,gen)) + 
     &        ERho(indi,i,tt,gen)*(2*(muu1*dconjg(GLho(indi,j,tt,gen))*
     &               (muu3*(muu2**2 + muu4**2 - tred)*dconjg(HLho(indj,j,tt,gen)) - 
     &                 muu2*(-1 + muu3**2 + muu4**2)*dconjg(HRho(indj,j,tt,gen))) + 
     &              muu4*dconjg(GRho(indi,j,tt,gen))*(-((-1 + muu1**2 + muu2**2)*muu3*dconjg(HLho(indj,j,tt,gen))) + 
     &                 muu2*(muu1**2 + muu3**2 - tred)*dconjg(HRho(indj,j,tt,gen))))*FLho(indj,i,tt,gen) + 
     &           (2*muu1*muu4*dconjg(GLho(indi,j,tt,gen))*
     &               (2*muu2*muu3*dconjg(HLho(indj,j,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HRho(indj,j,tt,gen))) + 
     &              dconjg(GRho(indi,j,tt,gen))*(2*muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HLho(indj,j,tt,gen)) + 
     &                 (-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - (muu2**2 + muu3**2 + muu4**2 - tred)*tred + 
     &                    muu1**2*(1 + 2*muu2**2 - tred - ured) - (muu2**2 + muu3**2 + muu4**2)*ured + ured**2)*
     &                  dconjg(HRho(indj,j,tt,gen))))*FRho(indj,i,tt,gen))))/
     &    ((-tred + MSf(i,itypeT,gen)**2/sqrtS**2)*(-ured + MSf(j,itypeU,gen)**2/sqrtS**2)))

	  enddo
	enddo

	do i = 1, 2

	if (icha.eq.2) then
		
! TS-channel photon (for cha cha -> t tbar)
! relative factor (-1)(-1) from vector propagator, T-ch.spinor order

	ChiChi2QQ_M2TreeHO = ChiChi2QQ_M2TreeHO + (-4*ApL*BpL*NC*(ELcho(indi,i,tt,gen)*(2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*FLho(indj,i,tt,gen) + 
     &         (muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)*FRho(indj,i,tt,gen)) + 
     &      ERho(indi,i,tt,gen)*((muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)*
     &          FLho(indj,i,tt,gen) + 2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*FRho(indj,i,tt,gen))))/
     &  (-tred + MSf(i,itypeT,gen)**2/sqrtS**2)

	endif	

! TS-channel vector-boson (Z/W)
! relative factor (-1)(-1) from vector propagator, T-ch.spinor order

	ChiChi2QQ_M2TreeHO = ChiChi2QQ_M2TreeHO + (-4*NC*sqrtS**2*(dconjg(AL)*(BR*((muu1**2 + muu3**2 - tred)*ELcho(indi,i,tt,gen) + 2*muu1*muu3*ERho(indi,i,tt,gen))*
     &          (2*muu2*muu4*FLho(indj,i,tt,gen) + (muu2**2 + muu4**2 - tred)*FRho(indj,i,tt,gen)) + 
     &         BL*(muu3*ELcho(indi,i,tt,gen)*(muu2*(muu1**2 + muu4**2 - ured)*FLho(indj,i,tt,gen) - 
     &               (-1 + muu1**2 + muu2**2)*muu4*FRho(indj,i,tt,gen)) + 
     &            muu1*ERho(indi,i,tt,gen)*(-(muu2*(-1 + muu3**2 + muu4**2)*FLho(indj,i,tt,gen)) + 
     &               muu4*(muu2**2 + muu3**2 - ured)*FRho(indj,i,tt,gen)))) + 
     &      dconjg(AR)*(BL*(2*muu1*muu3*ELcho(indi,i,tt,gen) + (muu1**2 + muu3**2 - tred)*ERho(indi,i,tt,gen))*
     &          ((muu2**2 + muu4**2 - tred)*FLho(indj,i,tt,gen) + 2*muu2*muu4*FRho(indj,i,tt,gen)) + 
     &         BR*(muu1*ELcho(indi,i,tt,gen)*(muu4*(muu2**2 + muu3**2 - ured)*FLho(indj,i,tt,gen) - 
     &               muu2*(-1 + muu3**2 + muu4**2)*FRho(indj,i,tt,gen)) + 
     &            muu3*ERho(indi,i,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu4*FLho(indj,i,tt,gen)) + 
     &               muu2*(muu1**2 + muu4**2 - ured)*FRho(indj,i,tt,gen))))))/(DenSmV1*(-tred + MSf(i,itypeT,gen)**2/sqrtS**2))

	do j = jmin, jmax


! TS-channel Higgs boson
! relative factor (-1) from T-ch.spinor order

	ChiChi2QQ_M2TreeHO = ChiChi2QQ_M2TreeHO + dreal(NC*sqrtS**2*(dconjg(CL(j))*(dconjg(DLho(j))*(2*muu1*muu3*ERho(indi,i,tt,gen)*
     &             ((muu2**2 + muu4**2 - tred)*FLho(indj,i,tt,gen) + 2*muu2*muu4*FRho(indj,i,tt,gen)) + 
     &            ELcho(indi,i,tt,gen)*((1 - muu3**2 - muu4**2 - muu3**2*tred - muu4**2*tred + tred**2 + muu3**2*ured + 
     &                  muu4**2*ured - ured**2 + muu2**2*(-1 + 2*muu3**2 - tred + ured) + muu1**2*(-1 + 2*muu4**2 - tred + ured))
     &                 *FLho(indj,i,tt,gen) + 2*muu2*muu4*(muu1**2 + muu3**2 - tred)*FRho(indj,i,tt,gen))) + 
     &         2*dconjg(DRho(j))*(muu1*ERho(indi,i,tt,gen)*
     &             (-(muu4*(muu2**2 + muu3**2 - ured)*FLho(indj,i,tt,gen)) + muu2*(-1 + muu3**2 + muu4**2)*FRho(indj,i,tt,gen)) + 
     &            muu3*ELcho(indi,i,tt,gen)*((-1 + muu1**2 + muu2**2)*muu4*FLho(indj,i,tt,gen) - 
     &               muu2*(muu1**2 + muu4**2 - ured)*FRho(indj,i,tt,gen)))) + 
     &      dconjg(CR(j))*(2*dconjg(DLho(j))*(muu3*ERho(indi,i,tt,gen)*
     &             (-(muu2*(muu1**2 + muu4**2 - ured)*FLho(indj,i,tt,gen)) + (-1 + muu1**2 + muu2**2)*muu4*FRho(indj,i,tt,gen)) + 
     &            muu1*ELcho(indi,i,tt,gen)*(muu2*(-1 + muu3**2 + muu4**2)*FLho(indj,i,tt,gen) - 
     &               muu4*(muu2**2 + muu3**2 - ured)*FRho(indj,i,tt,gen))) + 
     &         dconjg(DRho(j))*(2*muu1*muu3*ELcho(indi,i,tt,gen)*
     &             (2*muu2*muu4*FLho(indj,i,tt,gen) + (muu2**2 + muu4**2 - tred)*FRho(indj,i,tt,gen)) + 
     &            ERho(indi,i,tt,gen)*(2*muu2*muu4*(muu1**2 + muu3**2 - tred)*FLho(indj,i,tt,gen) + 
     &               (1 - muu3**2 - muu4**2 - muu3**2*tred - muu4**2*tred + tred**2 + muu3**2*ured + muu4**2*ured - ured**2 + 
     &                  muu2**2*(-1 + 2*muu3**2 - tred + ured) + muu1**2*(-1 + 2*muu4**2 - tred + ured))*FRho(indj,i,tt,gen))))))/
     &  (DenSmH1(j)*(-tred + MSf(i,itypeT,gen)**2/sqrtS**2))

	  enddo
	enddo


! U-channel 
! (zero for down-type quarks in chargino-chargino annihilation through zero couplings)

	do i = 1, 2
	  do j = 1, 2

	ChiChi2QQ_M2TreeHO = ChiChi2QQ_M2TreeHO -((NC*(((muu1**2 + muu4**2 - ured)*dconjg(GLho(indi,j,tt,gen)) + 2*muu1*muu4*dconjg(GRho(indi,j,tt,gen)))*GLho(indi,i,tt,gen) + 
     &        (2*muu1*muu4*dconjg(GLho(indi,j,tt,gen)) + (muu1**2 + muu4**2 - ured)*dconjg(GRho(indi,j,tt,gen)))*GRho(indi,i,tt,gen))*
     &      (-(((muu2**2 + muu3**2 - ured)*dconjg(HLho(indj,j,tt,gen)) + 2*muu2*muu3*dconjg(HRho(indj,j,tt,gen)))*
     &           HLho(indj,i,tt,gen)) - (2*muu2*muu3*dconjg(HLho(indj,j,tt,gen)) + 
     &           (muu2**2 + muu3**2 - ured)*dconjg(HRho(indj,j,tt,gen)))*HRho(indj,i,tt,gen)))/
     &    ((-ured + MSf(i,itypeU,gen)**2/sqrtS**2)*(-ured + MSf(j,itypeU,gen)**2/sqrtS**2)))
	  
	  enddo
	enddo
	
	
      do i = 1, 2

	if (icha.eq.2) then

! US-channel photon (for cha cha -> b bbar)
! relative factor (-1)(-1)(-1) from vector propagator, S-ch.spinor order, Charge conjugated vector coupling

	ChiChi2QQ_M2TreeHO = ChiChi2QQ_M2TreeHO + (4*ApL*BpL*NC*(GRho(indi,i,tt,gen)*((muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + 
     &            ured**2)*HLho(indj,i,tt,gen) + 2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*HRho(indj,i,tt,gen)) + 
     &      GLho(indi,i,tt,gen)*(2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*HLho(indj,i,tt,gen) + 
     &         (muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + ured**2)*HRho(indj,i,tt,gen))))/
     &  (-ured + MSf(i,itypeU,gen)**2/sqrtS**2)
		
	endif

! US-channel vector-boson (Z/W)
! relative factor (-1)(-1)(-1) from vector propagator, S-ch.spinor order, Charge conjugated vector coupling

	ChiChi2QQ_M2TreeHO = ChiChi2QQ_M2TreeHO + dreal(4*NC*sqrtS**2*(dconjg(AR)*(BR*(2*muu1*muu4*GLho(indi,i,tt,gen) + (muu1**2 + muu4**2 - ured)*GRho(indi,i,tt,gen))*
     &          ((muu2**2 + muu3**2 - ured)*HLho(indj,i,tt,gen) + 2*muu2*muu3*HRho(indj,i,tt,gen)) + 
     &         BL*(muu1*GLho(indi,i,tt,gen)*(muu3*(muu2**2 + muu4**2 - tred)*HLho(indj,i,tt,gen) - 
     &               muu2*(-1 + muu3**2 + muu4**2)*HRho(indj,i,tt,gen)) + 
     &            muu4*GRho(indi,i,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu3*HLho(indj,i,tt,gen)) + 
     &               muu2*(muu1**2 + muu3**2 - tred)*HRho(indj,i,tt,gen)))) + 
     &      dconjg(AL)*(BL*((muu1**2 + muu4**2 - ured)*GLho(indi,i,tt,gen) + 2*muu1*muu4*GRho(indi,i,tt,gen))*
     &          (2*muu2*muu3*HLho(indj,i,tt,gen) + (muu2**2 + muu3**2 - ured)*HRho(indj,i,tt,gen)) + 
     &         BR*(muu4*GLho(indi,i,tt,gen)*(muu2*(muu1**2 + muu3**2 - tred)*HLho(indj,i,tt,gen) - 
     &               (-1 + muu1**2 + muu2**2)*muu3*HRho(indj,i,tt,gen)) + 
     &            muu1*GRho(indi,i,tt,gen)*(-(muu2*(-1 + muu3**2 + muu4**2)*HLho(indj,i,tt,gen)) + 
     &               muu3*(muu2**2 + muu4**2 - tred)*HRho(indj,i,tt,gen))))))/(DenSmV1*(-ured + MSf(i,itypeU,gen)**2/sqrtS**2))

	  do j = jmin, jmax


! US-channel Higgs boson
! relative factor (-1) from S-ch.spinor order

	ChiChi2QQ_M2TreeHO = ChiChi2QQ_M2TreeHO + dreal(NC*sqrtS**2*(dconjg(CL(j))*(2*dconjg(DRho(j))*
     &          (muu1*GRho(indi,i,tt,gen)*(-(muu3*(muu2**2 + muu4**2 - tred)*HLho(indj,i,tt,gen)) + 
     &               muu2*(-1 + muu3**2 + muu4**2)*HRho(indj,i,tt,gen)) + 
     &            muu4*GLho(indi,i,tt,gen)*((-1 + muu1**2 + muu2**2)*muu3*HLho(indj,i,tt,gen) - 
     &               muu2*(muu1**2 + muu3**2 - tred)*HRho(indj,i,tt,gen))) + 
     &         dconjg(DLho(j))*(2*muu1*muu4*GRho(indi,i,tt,gen)*
     &             ((muu2**2 + muu3**2 - ured)*HLho(indj,i,tt,gen) + 2*muu2*muu3*HRho(indj,i,tt,gen)) + 
     &            GLho(indi,i,tt,gen)*(((-1 + muu3**2 + muu4**2 - tred)*(-1 + tred) + muu1**2*(-1 + 2*muu3**2 + tred - ured) + 
     &                  muu2**2*(-1 + 2*muu4**2 + tred - ured) - (muu3**2 + muu4**2)*ured + ured**2)*HLho(indj,i,tt,gen) + 
     &               2*muu2*muu3*(muu1**2 + muu4**2 - ured)*HRho(indj,i,tt,gen)))) + 
     &      dconjg(CR(j))*(2*dconjg(DLho(j))*(muu4*GRho(indi,i,tt,gen)*
     &             (-(muu2*(muu1**2 + muu3**2 - tred)*HLho(indj,i,tt,gen)) + (-1 + muu1**2 + muu2**2)*muu3*HRho(indj,i,tt,gen)) + 
     &            muu1*GLho(indi,i,tt,gen)*(muu2*(-1 + muu3**2 + muu4**2)*HLho(indj,i,tt,gen) - 
     &               muu3*(muu2**2 + muu4**2 - tred)*HRho(indj,i,tt,gen))) + 
     &         dconjg(DRho(j))*(2*muu1*muu4*GLho(indi,i,tt,gen)*
     &             (2*muu2*muu3*HLho(indj,i,tt,gen) + (muu2**2 + muu3**2 - ured)*HRho(indj,i,tt,gen)) + 
     &            GRho(indi,i,tt,gen)*(2*muu2*muu3*(muu1**2 + muu4**2 - ured)*HLho(indj,i,tt,gen) + 
     &               ((-1 + muu3**2 + muu4**2 - tred)*(-1 + tred) + muu1**2*(-1 + 2*muu3**2 + tred - ured) + 
     &                  muu2**2*(-1 + 2*muu4**2 + tred - ured) - (muu3**2 + muu4**2)*ured + ured**2)*HRho(indj,i,tt,gen))))))/
     &  (dconjg(DenSmH(j))*(-ured + MSf(i,itypeU,gen)**2/sqrtS**2))

	  enddo
	enddo

c	print*,'ChiChi2QQ_M2TreeHO', ChiChi2QQ_M2TreeHO

	end
