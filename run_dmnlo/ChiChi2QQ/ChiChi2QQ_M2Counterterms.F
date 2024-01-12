!************************************************

	function ChiChi2QQ_M2VertexSCounter()

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"
#include "ChiChi2QQ_GenCounterterms.h"

	integer i,j,k,l,m,tcount,imin,imax,jmin,jmax,itypeT,itypeU,att,antitt
	complex*16 DenSmV, DenSmH(6)
	double precision DenSmV1, DenSmH1(6), tempCounter, tempH, tempZ, tempt, tempu

	double precision ChiChi2QQ_M2VertexSCounter

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


	ChiChi2QQ_M2VertexSCounter = 0d0
	tempCounter = 0d0
	tempH = 0d0
	tempZ = 0d0
	tempu = 0d0
	tempt = 0d0

**********************************
* Counterterms for photon vertex *
**********************************

! for photon & interference only diagonal case considered (muu1=muu2, muu3=muu4)

	if (icha.eq.2) then
		
! Photon counterterm interference with photon tree-level
	
		ChiChi2QQ_M2VertexSCounter = ChiChi2QQ_M2VertexSCounter + 8*ApL**2*BpL*(dBpL + dBpR)*NC*(1 + 2*(muu1**2 + muu3**2 - tred)**2 + 2*tred)	
	
cc	write(*,*)'ChiChi P Vertex SCounterterm mit P = ', ChiChi2QQ_M2VertexSCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexSCounter
	
! Photon counterterm interference with Z-boson tree-level

		ChiChi2QQ_M2VertexSCounter = ChiChi2QQ_M2VertexSCounter + (8*ApL*NC*sqrtS**2*((BR*(dBpL*(1 + 2*muu1**2)*muu3**2 + 
     -            dBpR*(muu1**2 + muu1**4 + muu3**4 - 2*(muu1**2 + muu3**2)*tred + tred**2)) + 
     -         BL*(dBpR*(1 + 2*muu1**2)*muu3**2 + dBpL*(muu1**4 + (1 - muu3**2 + tred)**2 - muu1**2*(1 + 2*tred))))*dconjg(AL) + 
     -      (BL*(dBpR*(1 + 2*muu1**2)*muu3**2 + dBpL*(muu1**2 + muu1**4 + muu3**4 - 2*(muu1**2 + muu3**2)*tred + tred**2)) + 
     -         BR*(dBpL*(1 + 2*muu1**2)*muu3**2 + dBpR*(muu1**4 + (1 - muu3**2 + tred)**2 - muu1**2*(1 + 2*tred))))*dconjg(AR)))/
     -  (dreal(dconjg(DenSmV)))
	
cc	write(*,*)'ChiChi P Vertex SCounterterm mit Z = ', ChiChi2QQ_M2VertexSCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexSCounter
	
! Photon counterterm interference with Higgs-boson tree-level

	do i= imin,imax

	 	ChiChi2QQ_M2VertexSCounter = ChiChi2QQ_M2VertexSCounter + (-4*ApL*(dBpL + dBpR)*muu1*muu3*NC*sqrtS**2*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred)*(dconjg(CL(i)) + dconjg(CR(i)))*
     -    (dconjg(DL(i)) + dconjg(DR(i))))/(dreal(dconjg(DenSmH(i))))
	
	enddo
	
cc	write(*,*)'ChiChi P Vertex SCounterterm mit H = ', ChiChi2QQ_M2VertexSCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexSCounter
	
	do j=1,2
		
! Photon counterterm interference with T-channel tree-level (old U-channel)

	 	ChiChi2QQ_M2VertexSCounter = ChiChi2QQ_M2VertexSCounter + (-4*ApL*NC*(dBpR*(dconjg(ELc(indi,j,tt,gen))*
     -          (muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*dconjg(FL(indj,j,tt,gen)) + 
     -            (muu1**2 + muu1**4 + muu3**4 - 2*(muu1**2 + muu3**2)*tred + tred**2)*dconjg(FR(indj,j,tt,gen))) + 
     -         muu3*dconjg(ER(indi,j,tt,gen))*((1 + 2*muu1**2)*muu3*dconjg(FL(indj,j,tt,gen)) + 
     -            muu1*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*dconjg(FR(indj,j,tt,gen)))) + 
     -      dBpL*(muu3*dconjg(ELc(indi,j,tt,gen))*(muu1*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*dconjg(FL(indj,j,tt,gen)) + 
     -            (1 + 2*muu1**2)*muu3*dconjg(FR(indj,j,tt,gen))) + 
     -         dconjg(ER(indi,j,tt,gen))*((muu1**2 + muu1**4 + muu3**4 - 2*(muu1**2 + muu3**2)*tred + tred**2)*
     -             dconjg(FL(indj,j,tt,gen)) + muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*dconjg(FR(indj,j,tt,gen))))))/
     -  (-tred + MSf(j,itypeT,gen)**2/sqrtS**2)
	
	enddo
	
cc	write(*,*)'ChiChi P Vertex SCounterterm mit t = ', ChiChi2QQ_M2VertexSCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexSCounter
	
	do j=1,2
	
! Photon counterterm interference with U-channel tree-level (old T-channel)

		ChiChi2QQ_M2VertexSCounter = ChiChi2QQ_M2VertexSCounter + (4*ApL*NC*(dBpR*(muu3*dconjg(GL(indi,j,tt,gen))*
     -          (muu1*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*dconjg(HL(indj,j,tt,gen)) + 
     -            (1 + 2*muu1**2)*muu3*dconjg(HR(indj,j,tt,gen))) + 
     -         dconjg(GR(indi,j,tt,gen))*((muu1**2 + muu1**4 + muu3**4 - 2*(muu1**2 + muu3**2)*ured + ured**2)*
     -             dconjg(HL(indj,j,tt,gen)) + muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*dconjg(HR(indj,j,tt,gen)))) + 
     -      dBpL*(muu3*dconjg(GR(indi,j,tt,gen))*((1 + 2*muu1**2)*muu3*dconjg(HL(indj,j,tt,gen)) + 
     -            muu1*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*dconjg(HR(indj,j,tt,gen))) + 
     -         dconjg(GL(indi,j,tt,gen))*(muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*dconjg(HL(indj,j,tt,gen)) + 
     -            (muu1**2 + muu1**4 + muu3**4 - 2*(muu1**2 + muu3**2)*ured + ured**2)*dconjg(HR(indj,j,tt,gen))))))/
     -  (-ured + MSf(j,itypeU,gen)**2/sqrtS**2)
	 
	enddo
	
cc	write(*,*)'ChiChi P Vertex SCounterterm mit u = ', ChiChi2QQ_M2VertexSCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexSCounter

	endif

*************************************
* Counterterms for Z/W-boson vertex *
*************************************

! Z-boson counterterm interference with photon tree-level

	if (icha.eq.2) then

		ChiChi2QQ_M2VertexSCounter = ChiChi2QQ_M2VertexSCounter + (8*ApL*BpL*NC*sqrtS**2*(AL*(dBR*(muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*tred) + (muu3**2 - tred)**2) + 
     -         dBL*(1 - muu3**2 + ((-1 + muu1)*muu1 + muu3**2 - tred)*(muu1 + muu1**2 + muu3**2 - tred) + 2*tred)) + 
     -      AR*(dBL*muu1**2 + dBL*(muu3**2 + (muu1**2 + muu3**2 - tred)**2) + 
     -         dBR*(1 - muu3**2 + ((-1 + muu1)*muu1 + muu3**2 - tred)*(muu1 + muu1**2 + muu3**2 - tred) + 2*tred))))/(dreal(dconjg(DenSmV)))

	endif

! Z-boson counterterm interference with Z-boson tree-level

    	ChiChi2QQ_M2VertexSCounter = ChiChi2QQ_M2VertexSCounter + (-8*NC*sqrtS**4*((AR*muu1*muu2*((BL*dBL + BR*dBR)*(-1 + muu3**2) - 4*(BR*dBL + BL*dBR)*muu3*muu4 + 
     -            (BL*dBL + BR*dBR)*muu4**2) + AL*(BR*
     -             (dBL*(-1 + muu1**2 + muu2**2)*muu3*muu4 - dBR*(muu1**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred)) + 
     -            BL*(dBR*(-1 + muu1**2 + muu2**2)*muu3*muu4 - dBL*(muu2**2 + muu3**2 - ured)*(muu1**2 + muu4**2 - ured))))*
     -       dconjg(AL) + (AL*muu1*muu2*((BL*dBL + BR*dBR)*(-1 + muu3**2) - 4*(BR*dBL + BL*dBR)*muu3*muu4 + 
     -            (BL*dBL + BR*dBR)*muu4**2) + AR*(BL*
     -             (dBR*(-1 + muu1**2 + muu2**2)*muu3*muu4 - dBL*(muu1**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred)) + 
     -            BR*(dBL*(-1 + muu1**2 + muu2**2)*muu3*muu4 - dBR*(muu2**2 + muu3**2 - ured)*(muu1**2 + muu4**2 - ured))))*
     -       dconjg(AR)))/(DenSmV*dconjg(DenSmV))

cc	write(*,*)'ChiChi Z Vertex SCounterterm mit Z = ', ChiChi2QQ_M2VertexSCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexSCounter

! Z-boson counterterm interference with Higgs-boson tree-level

	do i= imin,imax

		ChiChi2QQ_M2VertexSCounter = ChiChi2QQ_M2VertexSCounter + (-4*NC*sqrtS**4*(muu1*(AR*dconjg(CL(i)) + AL*dconjg(CR(i)))*
     -       (-(muu4*(muu2**2 + muu3**2 - ured)*(dBR*dconjg(DL(i)) + dBL*dconjg(DR(i)))) + 
     -         muu3*(muu2**2 + muu4**2 - tred)*(dBL*dconjg(DL(i)) + dBR*dconjg(DR(i)))) - 
     -      muu2*(AL*dconjg(CL(i)) + AR*dconjg(CR(i)))*
     -       (-(muu4*(muu1**2 + muu3**2 - tred)*(dBR*dconjg(DL(i)) + dBL*dconjg(DR(i)))) + 
     -         muu3*(muu1**2 + muu4**2 - ured)*(dBL*dconjg(DL(i)) + dBR*dconjg(DR(i))))))/(dreal(DenSmV*dconjg(DenSmH(i))))

	enddo

cc	write(*,*)'ChiChi Z Vertex SCounterterm mit H = ', ChiChi2QQ_M2VertexSCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexSCounter

	do j=1,2

! Z-boson counterterm interference with U-channel tree-level (old T-channel)

 	 ChiChi2QQ_M2VertexSCounter = ChiChi2QQ_M2VertexSCounter + (4*NC*sqrtS**2*(AR*(dBR*(2*muu1*muu4*dconjg(GL(indi,j,tt,gen)) + (muu1**2 + muu4**2 - ured)*dconjg(GR(indi,j,tt,gen)))*
     -          ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,j,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,j,tt,gen))) + 
     -         dBL*(muu1*dconjg(GL(indi,j,tt,gen))*(muu3*(muu2**2 + muu4**2 - tred)*dconjg(HL(indj,j,tt,gen)) - 
     -               muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,j,tt,gen))) + 
     -            muu4*dconjg(GR(indi,j,tt,gen))*(-((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,j,tt,gen))) + 
     -               muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,j,tt,gen))))) + 
     -      AL*(dBL*((muu1**2 + muu4**2 - ured)*dconjg(GL(indi,j,tt,gen)) + 2*muu1*muu4*dconjg(GR(indi,j,tt,gen)))*
     -          (2*muu2*muu3*dconjg(HL(indj,j,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,j,tt,gen))) + 
     -         dBR*(muu4*dconjg(GL(indi,j,tt,gen))*(muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,j,tt,gen)) - 
     -               (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,j,tt,gen))) + 
     -            muu1*dconjg(GR(indi,j,tt,gen))*(-(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,j,tt,gen))) + 
     -               muu3*(muu2**2 + muu4**2 - tred)*dconjg(HR(indj,j,tt,gen)))))))/
     -  (DenSmV1*(-ured + MSf(j,itypeU,gen)**2/sqrtS**2))
	
	enddo
	
cc	write(*,*)'ChiChi Z Vertex SCounterterm mit u = ', ChiChi2QQ_M2VertexSCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexSCounter

	do j=1,2
	
! Z-boson counterterm interference with T-channel tree-level (old U-channel)

 	 ChiChi2QQ_M2VertexSCounter = ChiChi2QQ_M2VertexSCounter + (-4*NC*sqrtS**2*(AL*(dBR*((muu1**2 + muu3**2 - tred)*dconjg(ELc(indi,j,tt,gen)) + 2*muu1*muu3*dconjg(ER(indi,j,tt,gen)))*
     -          (2*muu2*muu4*dconjg(FL(indj,j,tt,gen)) + (muu2**2 + muu4**2 - tred)*dconjg(FR(indj,j,tt,gen))) + 
     -         dBL*(muu3*dconjg(ELc(indi,j,tt,gen))*(muu2*(muu1**2 + muu4**2 - ured)*dconjg(FL(indj,j,tt,gen)) - 
     -               (-1 + muu1**2 + muu2**2)*muu4*dconjg(FR(indj,j,tt,gen))) + 
     -            muu1*dconjg(ER(indi,j,tt,gen))*(-(muu2*(-1 + muu3**2 + muu4**2)*dconjg(FL(indj,j,tt,gen))) + 
     -               muu4*(muu2**2 + muu3**2 - ured)*dconjg(FR(indj,j,tt,gen))))) + 
     -      AR*(dBL*(2*muu1*muu3*dconjg(ELc(indi,j,tt,gen)) + (muu1**2 + muu3**2 - tred)*dconjg(ER(indi,j,tt,gen)))*
     -          ((muu2**2 + muu4**2 - tred)*dconjg(FL(indj,j,tt,gen)) + 2*muu2*muu4*dconjg(FR(indj,j,tt,gen))) + 
     -         dBR*(muu1*dconjg(ELc(indi,j,tt,gen))*(muu4*(muu2**2 + muu3**2 - ured)*dconjg(FL(indj,j,tt,gen)) - 
     -               muu2*(-1 + muu3**2 + muu4**2)*dconjg(FR(indj,j,tt,gen))) + 
     -            muu3*dconjg(ER(indi,j,tt,gen))*(-((-1 + muu1**2 + muu2**2)*muu4*dconjg(FL(indj,j,tt,gen))) + 
     -               muu2*(muu1**2 + muu4**2 - ured)*dconjg(FR(indj,j,tt,gen)))))))/
     -  (DenSmV1*(-tred + MSf(j,itypeT,gen)**2/sqrtS**2))

	 enddo
	
cc	write(*,*)'ChiChi Z Vertex SCounterterm mit t = ', ChiChi2QQ_M2VertexSCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexSCounter

*************************
* H vertex counterterms *
*************************

	do i=imin,imax

! Higgs-boson counterterm interference with photon tree-level

	if (icha.eq.2) then

		ChiChi2QQ_M2VertexSCounter = ChiChi2QQ_M2VertexSCounter + (8*ApL*BpL*muu1*muu3*NC*sqrtS**2*(tred - ured)*(CL(i) + CR(i))*(dDL(i) + dDR(i)))/(dreal(dconjg(DenSmH(i))))

	endif

! Higgs-boson counterterm interference with Z/W-boson tree-level

	ChiChi2QQ_M2VertexSCounter = ChiChi2QQ_M2VertexSCounter + (-4*NC*sqrtS**4*(-(muu4*(muu1*(muu2**2 + muu3**2 - ured)*(CR(i)*dconjg(AL) + CL(i)*dconjg(AR)) - 
     -           muu2*(muu1**2 + muu3**2 - tred)*(CL(i)*dconjg(AL) + CR(i)*dconjg(AR)))*(dconjg(BR)*dDL(i) + dconjg(BL)*dDR(i)))
     -       + muu3*(muu1*(muu2**2 + muu4**2 - tred)*(CR(i)*dconjg(AL) + CL(i)*dconjg(AR)) - 
     -         muu2*(muu1**2 + muu4**2 - ured)*(CL(i)*dconjg(AL) + CR(i)*dconjg(AR)))*(dconjg(BL)*dDL(i) + dconjg(BR)*dDR(i))))/
     -  (dreal(DenSmV*dconjg(DenSmH(i))))
 
	tempZ = tempZ + ChiChi2QQ_M2VertexSCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexSCounter
	
! Higgs-boson counterterm interference with Higgs-boson tree-level

       do j= jmin,jmax

 	ChiChi2QQ_M2VertexSCounter = ChiChi2QQ_M2VertexSCounter + (2*NC*sqrtS**4*(2*muu1*muu2*(CR(i)*dconjg(CL(j)) + CL(i)*dconjg(CR(j))) + 
     -      (-1 + muu1**2)*(CL(i)*dconjg(CL(j)) + CR(i)*dconjg(CR(j))) + muu2**2*(CL(i)*dconjg(CL(j)) + CR(i)*dconjg(CR(j))))*
     -    (2*muu3*muu4*(dconjg(DR(j))*dDL(i) + dconjg(DL(j))*dDR(i)) + 
     -      (-1 + muu3**2)*(dconjg(DL(j))*dDL(i) + dconjg(DR(j))*dDR(i)) + muu4**2*(dconjg(DL(j))*dDL(i) + dconjg(DR(j))*dDR(i))))
     -   /(dconjg(DenSmH(j))*DenSmH(i))

       enddo

	tempH = tempH + ChiChi2QQ_M2VertexSCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexSCounter
	
! Higgs-boson counterterm interference with U-channel tree-level (old T-channel)

       do j=1,2
	
	 ChiChi2QQ_M2VertexSCounter = ChiChi2QQ_M2VertexSCounter -((NC*sqrtS**2*(CL(i)*((-2*muu1*muu4*dconjg(GR(indi,j,tt,gen))*
     -               ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,j,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,j,tt,gen))) + 
     -              dconjg(GL(indi,j,tt,gen))*((-1 + muu3**2 + muu4**2 + muu2**2*(1 - 2*muu4**2) - 
     -                    (muu2**2 + muu3**2 + muu4**2 - tred - ured)*(tred - ured) + muu1**2*(1 - 2*muu3**2 - tred + ured))*
     -                  dconjg(HL(indj,j,tt,gen)) - 2*muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,j,tt,gen))))*dDL(i) + 
     -           2*(muu1*dconjg(GR(indi,j,tt,gen))*(muu3*(muu2**2 + muu4**2 - tred)*dconjg(HL(indj,j,tt,gen)) - 
     -                 muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,j,tt,gen))) + 
     -              muu4*dconjg(GL(indi,j,tt,gen))*(-((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,j,tt,gen))) + 
     -                 muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,j,tt,gen))))*dDR(i)) + 
     -        CR(i)*(2*(muu4*dconjg(GR(indi,j,tt,gen))*
     -               (muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,j,tt,gen)) - 
     -                 (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,j,tt,gen))) + 
     -              muu1*dconjg(GL(indi,j,tt,gen))*(-(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,j,tt,gen))) + 
     -                 muu3*(muu2**2 + muu4**2 - tred)*dconjg(HR(indj,j,tt,gen))))*dDL(i) + 
     -           (-2*muu1*muu4*dconjg(GL(indi,j,tt,gen))*
     -               (2*muu2*muu3*dconjg(HL(indj,j,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,j,tt,gen))) + 
     -              dconjg(GR(indi,j,tt,gen))*(-2*muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,j,tt,gen)) + 
     -                 (-1 + muu3**2 + muu4**2 + muu2**2*(1 - 2*muu4**2) - 
     -                    (muu2**2 + muu3**2 + muu4**2 - tred - ured)*(tred - ured) + muu1**2*(1 - 2*muu3**2 - tred + ured))*
     -                  dconjg(HR(indj,j,tt,gen))))*dDR(i))))/(DenSmH1(i)*(-ured + MSf(j,itypeU,gen)**2/sqrtS**2)))
	
	enddo
	
	tempu = tempu + ChiChi2QQ_M2VertexSCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexSCounter
	
	do j=1,2
	
! Higgs-boson counterterm interference with T-channel tree-level (old U-channel)

	 ChiChi2QQ_M2VertexSCounter = ChiChi2QQ_M2VertexSCounter -((NC*sqrtS**2*(CR(i)*(2*(muu3*dconjg(ER(indi,j,tt,gen))*
     -               (muu2*(muu1**2 + muu4**2 - ured)*dconjg(FL(indj,j,tt,gen)) - 
     -                 (-1 + muu1**2 + muu2**2)*muu4*dconjg(FR(indj,j,tt,gen))) + 
     -              muu1*dconjg(ELc(indi,j,tt,gen))*(-(muu2*(-1 + muu3**2 + muu4**2)*dconjg(FL(indj,j,tt,gen))) + 
     -                 muu4*(muu2**2 + muu3**2 - ured)*dconjg(FR(indj,j,tt,gen))))*dDL(i) + 
     -           (-2*muu1*muu3*dconjg(ELc(indi,j,tt,gen))*
     -               (2*muu2*muu4*dconjg(FL(indj,j,tt,gen)) + (muu2**2 + muu4**2 - tred)*dconjg(FR(indj,j,tt,gen))) + 
     -              dconjg(ER(indi,j,tt,gen))*(-2*muu2*muu4*(muu1**2 + muu3**2 - tred)*dconjg(FL(indj,j,tt,gen)) + 
     -                 (-1 + muu3**2 + muu2**2*(1 - 2*muu3**2) + muu4**2 + 
     -                    (muu2**2 + muu3**2 + muu4**2 - tred - ured)*(tred - ured) + muu1**2*(1 - 2*muu4**2 + tred - ured))*
     -                  dconjg(FR(indj,j,tt,gen))))*dDR(i)) + 
     -        CL(i)*((-2*muu1*muu3*dconjg(ER(indi,j,tt,gen))*
     -               ((muu2**2 + muu4**2 - tred)*dconjg(FL(indj,j,tt,gen)) + 2*muu2*muu4*dconjg(FR(indj,j,tt,gen))) + 
     -              dconjg(ELc(indi,j,tt,gen))*((-1 + muu3**2 + muu2**2*(1 - 2*muu3**2) + muu4**2 + 
     -                    (muu2**2 + muu3**2 + muu4**2 - tred - ured)*(tred - ured) + muu1**2*(1 - 2*muu4**2 + tred - ured))*
     -                  dconjg(FL(indj,j,tt,gen)) - 2*muu2*muu4*(muu1**2 + muu3**2 - tred)*dconjg(FR(indj,j,tt,gen))))*dDL(i) + 
     -           2*(muu1*dconjg(ER(indi,j,tt,gen))*(muu4*(muu2**2 + muu3**2 - ured)*dconjg(FL(indj,j,tt,gen)) - 
     -                 muu2*(-1 + muu3**2 + muu4**2)*dconjg(FR(indj,j,tt,gen))) + 
     -              muu3*dconjg(ELc(indi,j,tt,gen))*(-((-1 + muu1**2 + muu2**2)*muu4*dconjg(FL(indj,j,tt,gen))) + 
     -                 muu2*(muu1**2 + muu4**2 - ured)*dconjg(FR(indj,j,tt,gen))))*dDR(i))))/
     -    (DenSmH1(i)*(-tred + MSf(j,itypeT,gen)**2/sqrtS**2)))

	enddo

	tempt = tempt + ChiChi2QQ_M2VertexSCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexSCounter
	
	enddo

cc	write(*,*)'ChiChi H Vertex SCounterterm mit Z = ', tempZ
cc	write(*,*)'ChiChi H Vertex SCounterterm mit H = ', tempH
cc	write(*,*)'ChiChi H Vertex SCounterterm mit t = ', tempt
cc	write(*,*)'ChiChi H Vertex SCounterterm mit u = ', tempu
	
	end

************************************************

	function ChiChi2QQ_M2VertexTUCounter()

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"
#include "ChiChi2QQ_GenCounterterms.h"

	integer i,j,k,l,m,tcount,imin,imax,jmin,jmax,itypeT,itypeU,att,antitt
	complex*16 DenSmV, DenSmH(6)
	double precision DenSmV1, DenSmH1(6), tempCounter, tempH, tempZ, tempt, tempu

	double precision ChiChi2QQ_M2VertexTUCounter

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

	ChiChi2QQ_M2VertexTUCounter = 0d0
	tempCounter = 0d0
	tempH = 0d0
	tempZ = 0d0
	tempu = 0d0
	tempt = 0d0

*************************************************
* U-channel vertex counterterms (old T-channel) *
*************************************************

	do i=1,2

! (p2k1)-vertex with photon
	if (icha.eq.2) then

!neu Version

                ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (-4*ApL*BpL*NC*sqrtS**2*(dHR(indj,i,tt,gen)*((muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + ured**2)*GL(indi,i,tt,gen) + 
     -         2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*GR(indi,i,tt,gen)) + 
     -      dHL(indj,i,tt,gen)*(2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*GL(indi,i,tt,gen) + 
     -         (muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + ured**2)*GR(indi,i,tt,gen))))/(sqrtS**2*ured - MSf(i,itypeU,gen)**2)

!alte Version

! 		ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (4*ApL*BpL*NC*(dHL(indj,i,tt,gen)*(2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*GL(indi,i,tt,gen) + 
!      -         (muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*ured) + (muu3**2 - ured)**2)*GR(indi,i,tt,gen)) + 
!      -      dHR(indj,i,tt,gen)*((muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*ured) + (muu3**2 - ured)**2)*GL(indi,i,tt,gen) + 
!      -         2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*GR(indi,i,tt,gen))))/(-ured + MSf(i,itypeU,gen)**2/sqrtS**2)


! (p1k2)-vertex with photon

!neu Version

                ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (-4*ApL*BpL*NC*sqrtS**2*(dGR(indi,i,tt,gen)*((muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + ured**2)*HL(indj,i,tt,gen) + 
     -         2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*HR(indj,i,tt,gen)) + 
     -      dGL(indi,i,tt,gen)*(2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*HL(indj,i,tt,gen) + 
     -         (muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*ured) - 2*muu3**2*ured + ured**2)*HR(indj,i,tt,gen))))/(sqrtS**2*ured - MSf(i,itypeU,gen)**2)

!alte Version

! 		ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (4*ApL*BpL*NC*(dGL(indi,i,tt,gen)*(2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*HL(indj,i,tt,gen) + 
!      -         (muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*ured) + (muu3**2 - ured)**2)*HR(indj,i,tt,gen)) + 
!      -      dGR(indi,i,tt,gen)*((muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*ured) + (muu3**2 - ured)**2)*HL(indj,i,tt,gen) + 
!      -         2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*HR(indj,i,tt,gen))))/(-ured + MSf(i,itypeU,gen)**2/sqrtS**2)

	endif

! (p2k1)-vertex with Z

	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (4*NC*sqrtS**2*(dconjg(AR)*(BR*((muu2**2 + muu3**2 - ured)*dHL(indj,i,tt,gen) + 2*muu2*muu3*dHR(indj,i,tt,gen))*
     -          (2*muu1*muu4*GL(indi,i,tt,gen) + (muu1**2 + muu4**2 - ured)*GR(indi,i,tt,gen)) + 
     -         BL*(muu3*dHL(indj,i,tt,gen)*(muu1*(muu2**2 + muu4**2 - tred)*GL(indi,i,tt,gen) - 
     -               (-1 + muu1**2 + muu2**2)*muu4*GR(indi,i,tt,gen)) + 
     -            muu2*dHR(indj,i,tt,gen)*(-(muu1*(-1 + muu3**2 + muu4**2)*GL(indi,i,tt,gen)) + 
     -               muu4*(muu1**2 + muu3**2 - tred)*GR(indi,i,tt,gen)))) + 
     -      dconjg(AL)*(BL*(2*muu2*muu3*dHL(indj,i,tt,gen) + (muu2**2 + muu3**2 - ured)*dHR(indj,i,tt,gen))*
     -          ((muu1**2 + muu4**2 - ured)*GL(indi,i,tt,gen) + 2*muu1*muu4*GR(indi,i,tt,gen)) + 
     -         BR*(muu2*dHL(indj,i,tt,gen)*(muu4*(muu1**2 + muu3**2 - tred)*GL(indi,i,tt,gen) - 
     -               muu1*(-1 + muu3**2 + muu4**2)*GR(indi,i,tt,gen)) + 
     -            muu3*dHR(indj,i,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu4*GL(indi,i,tt,gen)) + 
     -               muu1*(muu2**2 + muu4**2 - tred)*GR(indi,i,tt,gen))))))/(DenSmV1*(-ured + MSf(i,itypeU,gen)**2/sqrtS**2))
	
	  
! (p1k2)-vertex with Z

	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (4*NC*sqrtS**2*(dconjg(AR)*(BR*(2*muu1*muu4*dGL(indi,i,tt,gen) + (muu1**2 + muu4**2 - ured)*dGR(indi,i,tt,gen))*
     -          ((muu2**2 + muu3**2 - ured)*HL(indj,i,tt,gen) + 2*muu2*muu3*HR(indj,i,tt,gen)) + 
     -         BL*(muu1*dGL(indi,i,tt,gen)*(muu3*(muu2**2 + muu4**2 - tred)*HL(indj,i,tt,gen) - 
     -               muu2*(-1 + muu3**2 + muu4**2)*HR(indj,i,tt,gen)) + 
     -            muu4*dGR(indi,i,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu3*HL(indj,i,tt,gen)) + 
     -               muu2*(muu1**2 + muu3**2 - tred)*HR(indj,i,tt,gen)))) + 
     -      dconjg(AL)*(BL*((muu1**2 + muu4**2 - ured)*dGL(indi,i,tt,gen) + 2*muu1*muu4*dGR(indi,i,tt,gen))*
     -          (2*muu2*muu3*HL(indj,i,tt,gen) + (muu2**2 + muu3**2 - ured)*HR(indj,i,tt,gen)) + 
     -         BR*(muu4*dGL(indi,i,tt,gen)*(muu2*(muu1**2 + muu3**2 - tred)*HL(indj,i,tt,gen) - 
     -               (-1 + muu1**2 + muu2**2)*muu3*HR(indj,i,tt,gen)) + 
     -            muu1*dGR(indi,i,tt,gen)*(-(muu2*(-1 + muu3**2 + muu4**2)*HL(indj,i,tt,gen)) + 
     -               muu3*(muu2**2 + muu4**2 - tred)*HR(indj,i,tt,gen))))))/(DenSmV1*(-ured + MSf(i,itypeU,gen)**2/sqrtS**2))

	tempZ = tempZ + ChiChi2QQ_M2VertexTUCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexTUCounter
	
! (p2k1)-vertex with Higgs
	
	do j= jmin,jmax

	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter -((NC*sqrtS**2*(dconjg(CL(j))*(-2*muu2*dHR(indj,i,tt,gen)*
     -            (muu3*dconjg(DL(j))*((muu1**2 + muu4**2 - ured)*GL(indi,i,tt,gen) + 2*muu1*muu4*GR(indi,i,tt,gen)) + 
     -              dconjg(DR(j))*(-(muu4*(muu1**2 + muu3**2 - tred)*GL(indi,i,tt,gen)) + 
     -                 muu1*(-1 + muu3**2 + muu4**2)*GR(indi,i,tt,gen))) + 
     -           dHL(indj,i,tt,gen)*(2*muu3*dconjg(DR(j))*
     -               (-((-1 + muu1**2 + muu2**2)*muu4*GL(indi,i,tt,gen)) + muu1*(muu2**2 + muu4**2 - tred)*GR(indi,i,tt,gen)) + 
     -              dconjg(DL(j))*((-1 + muu3**2 + muu4**2 + muu2**2*(1 - 2*muu4**2) - 
     -                    (muu2**2 + muu3**2 + muu4**2 - tred - ured)*(tred - ured) + muu1**2*(1 - 2*muu3**2 - tred + ured))*
     -                  GL(indi,i,tt,gen) - 2*muu1*muu4*(muu2**2 + muu3**2 - ured)*GR(indi,i,tt,gen)))) + 
     -        dconjg(CR(j))*(-2*muu2*dHL(indj,i,tt,gen)*
     -            (dconjg(DL(j))*(muu1*(-1 + muu3**2 + muu4**2)*GL(indi,i,tt,gen) - 
     -                 muu4*(muu1**2 + muu3**2 - tred)*GR(indi,i,tt,gen)) + 
     -              muu3*dconjg(DR(j))*(2*muu1*muu4*GL(indi,i,tt,gen) + (muu1**2 + muu4**2 - ured)*GR(indi,i,tt,gen))) + 
     -           dHR(indj,i,tt,gen)*(2*muu3*dconjg(DL(j))*
     -               (muu1*(muu2**2 + muu4**2 - tred)*GL(indi,i,tt,gen) - (-1 + muu1**2 + muu2**2)*muu4*GR(indi,i,tt,gen)) + 
     -              dconjg(DR(j))*(-2*muu1*muu4*(muu2**2 + muu3**2 - ured)*GL(indi,i,tt,gen) + 
     -                 (-1 + muu3**2 + muu4**2 + muu2**2*(1 - 2*muu4**2) - 
     -                    (muu2**2 + muu3**2 + muu4**2 - tred - ured)*(tred - ured) + muu1**2*(1 - 2*muu3**2 - tred + ured))*
     -                  GR(indi,i,tt,gen))))))/(DenSmH1(j)*(-ured + MSf(i,itypeU,gen)**2/sqrtS**2)))
     
! (p1k2)-vertex with Higgs

	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter -((NC*sqrtS**2*(dconjg(CL(j))*(-2*muu1*dGR(indi,i,tt,gen)*
     -            (muu4*dconjg(DL(j))*((muu2**2 + muu3**2 - ured)*HL(indj,i,tt,gen) + 2*muu2*muu3*HR(indj,i,tt,gen)) + 
     -              dconjg(DR(j))*(-(muu3*(muu2**2 + muu4**2 - tred)*HL(indj,i,tt,gen)) + 
     -                 muu2*(-1 + muu3**2 + muu4**2)*HR(indj,i,tt,gen))) + 
     -           dGL(indi,i,tt,gen)*(2*muu4*dconjg(DR(j))*
     -               (-((-1 + muu1**2 + muu2**2)*muu3*HL(indj,i,tt,gen)) + muu2*(muu1**2 + muu3**2 - tred)*HR(indj,i,tt,gen)) + 
     -              dconjg(DL(j))*((-1 + muu3**2 + muu4**2 + muu2**2*(1 - 2*muu4**2) - 
     -                    (muu2**2 + muu3**2 + muu4**2 - tred - ured)*(tred - ured) + muu1**2*(1 - 2*muu3**2 - tred + ured))*
     -                  HL(indj,i,tt,gen) - 2*muu2*muu3*(muu1**2 + muu4**2 - ured)*HR(indj,i,tt,gen)))) + 
     -        dconjg(CR(j))*(-2*muu1*dGL(indi,i,tt,gen)*
     -            (dconjg(DL(j))*(muu2*(-1 + muu3**2 + muu4**2)*HL(indj,i,tt,gen) - 
     -                 muu3*(muu2**2 + muu4**2 - tred)*HR(indj,i,tt,gen)) + 
     -              muu4*dconjg(DR(j))*(2*muu2*muu3*HL(indj,i,tt,gen) + (muu2**2 + muu3**2 - ured)*HR(indj,i,tt,gen))) + 
     -           dGR(indi,i,tt,gen)*(2*muu4*dconjg(DL(j))*
     -               (muu2*(muu1**2 + muu3**2 - tred)*HL(indj,i,tt,gen) - (-1 + muu1**2 + muu2**2)*muu3*HR(indj,i,tt,gen)) + 
     -              dconjg(DR(j))*(-2*muu2*muu3*(muu1**2 + muu4**2 - ured)*HL(indj,i,tt,gen) + 
     -                 (-1 + muu3**2 + muu4**2 + muu2**2*(1 - 2*muu4**2) - 
     -                    (muu2**2 + muu3**2 + muu4**2 - tred - ured)*(tred - ured) + muu1**2*(1 - 2*muu3**2 - tred + ured))*
     -                  HR(indj,i,tt,gen))))))/(DenSmH1(j)*(-ured + MSf(i,itypeU,gen)**2/sqrtS**2)))

	enddo
	
	tempH = tempH + ChiChi2QQ_M2VertexTUCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexTUCounter
	
! (p2k1)-vertex with U-channel (with old T-channel)

	do j=1,2

	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (2*NC*(((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,j,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,j,tt,gen)))*
     -       dHL(indj,i,tt,gen) + (2*muu2*muu3*dconjg(HL(indj,j,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,j,tt,gen)))*
     -       dHR(indj,i,tt,gen))*(((muu1**2 + muu4**2 - ured)*dconjg(GL(indi,j,tt,gen)) + 2*muu1*muu4*dconjg(GR(indi,j,tt,gen)))*
     -       GL(indi,i,tt,gen) + (2*muu1*muu4*dconjg(GL(indi,j,tt,gen)) + (muu1**2 + muu4**2 - ured)*dconjg(GR(indi,j,tt,gen)))*
     -       GR(indi,i,tt,gen)))/((-ured + MSf(i,itypeU,gen)**2/sqrtS**2)*(-ured + MSf(j,itypeU,gen)**2/sqrtS**2))

! (p1k2)-vertex with U-channel (with old T-channel)

	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (2*NC*(((muu1**2 + muu4**2 - ured)*dconjg(GL(indi,j,tt,gen)) + 2*muu1*muu4*dconjg(GR(indi,j,tt,gen)))*
     -       dGL(indi,i,tt,gen) + (2*muu1*muu4*dconjg(GL(indi,j,tt,gen)) + (muu1**2 + muu4**2 - ured)*dconjg(GR(indi,j,tt,gen)))*
     -       dGR(indi,i,tt,gen))*(((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,j,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,j,tt,gen)))*
     -       HL(indj,i,tt,gen) + (2*muu2*muu3*dconjg(HL(indj,j,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,j,tt,gen)))*
     -       HR(indj,i,tt,gen)))/((-ured + MSf(i,itypeU,gen)**2/sqrtS**2)*(-ured + MSf(j,itypeU,gen)**2/sqrtS**2))

	tempu = tempu + ChiChi2QQ_M2VertexTUCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexTUCounter
	
! (p2k1)-vertex with T-channel (with old U-channel)

	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (NC*(dHL(indj,i,tt,gen)*(dconjg(ELc(indi,j,tt,gen))*
     -          (2*muu2*dconjg(FR(indj,j,tt,gen))*(-(muu4*(muu1**2 + muu3**2 - tred)*GL(indi,i,tt,gen)) + 
     -               muu1*(-1 + muu3**2 + muu4**2)*GR(indi,i,tt,gen)) + 
     -            dconjg(FL(indj,j,tt,gen))*(-((-1 + muu1**2 + muu2**2 + 2*muu1**2*muu2**2 + muu3**2 + muu4**2 + 
     -                    2*muu3**2*muu4**2)*GL(indi,i,tt,gen)) + 
     -               (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred*GL(indi,i,tt,gen) - tred**2*GL(indi,i,tt,gen) + 
     -               (muu1**2 + muu2**2 + muu3**2 + muu4**2)*ured*GL(indi,i,tt,gen) - ured**2*GL(indi,i,tt,gen) - 
     -               2*muu1*muu4*(muu2**2 + muu3**2 - ured)*GR(indi,i,tt,gen))) - 
     -         2*muu3*dconjg(ER(indi,j,tt,gen))*(dconjg(FL(indj,j,tt,gen))*
     -             (muu1*(muu2**2 + muu4**2 - tred)*GL(indi,i,tt,gen) - (-1 + muu1**2 + muu2**2)*muu4*GR(indi,i,tt,gen)) + 
     -            muu2*dconjg(FR(indj,j,tt,gen))*(2*muu1*muu4*GL(indi,i,tt,gen) + (muu1**2 + muu4**2 - ured)*GR(indi,i,tt,gen))))
     -       + dHR(indj,i,tt,gen)*(-2*muu3*dconjg(ELc(indi,j,tt,gen))*
     -          (muu2*dconjg(FL(indj,j,tt,gen))*((muu1**2 + muu4**2 - ured)*GL(indi,i,tt,gen) + 2*muu1*muu4*GR(indi,i,tt,gen)) + 
     -            dconjg(FR(indj,j,tt,gen))*(-((-1 + muu1**2 + muu2**2)*muu4*GL(indi,i,tt,gen)) + 
     -               muu1*(muu2**2 + muu4**2 - tred)*GR(indi,i,tt,gen))) + 
     -         dconjg(ER(indi,j,tt,gen))*(2*muu2*dconjg(FL(indj,j,tt,gen))*
     -             (muu1*(-1 + muu3**2 + muu4**2)*GL(indi,i,tt,gen) - muu4*(muu1**2 + muu3**2 - tred)*GR(indi,i,tt,gen)) + 
     -            dconjg(FR(indj,j,tt,gen))*(-2*muu1*muu4*(muu2**2 + muu3**2 - ured)*GL(indi,i,tt,gen) - 
     -               (-1 + muu1**2 + muu2**2 + 2*muu1**2*muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2)*GR(indi,i,tt,gen) + 
     -               (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred*GR(indi,i,tt,gen) - tred**2*GR(indi,i,tt,gen) + 
     -               (muu1**2 + muu2**2 + muu3**2 + muu4**2)*ured*GR(indi,i,tt,gen) - ured**2*GR(indi,i,tt,gen))))))/
     -  ((-ured + MSf(i,itypeU,gen)**2/sqrtS**2)*(-tred + MSf(j,itypeT,gen)**2/sqrtS**2))

! (p1k2)-vertex with T-channel (with old U-channel)

	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (NC*(dGL(indi,i,tt,gen)*(-2*muu1*dconjg(ER(indi,j,tt,gen))*
     -          (dconjg(FL(indj,j,tt,gen))*(muu3*(muu2**2 + muu4**2 - tred)*HL(indj,i,tt,gen) - 
     -               muu2*(-1 + muu3**2 + muu4**2)*HR(indj,i,tt,gen)) + 
     -            muu4*dconjg(FR(indj,j,tt,gen))*(2*muu2*muu3*HL(indj,i,tt,gen) + (muu2**2 + muu3**2 - ured)*HR(indj,i,tt,gen)))
     -          + dconjg(ELc(indi,j,tt,gen))*(2*muu4*dconjg(FR(indj,j,tt,gen))*
     -             (-(muu2*(muu1**2 + muu3**2 - tred)*HL(indj,i,tt,gen)) + (-1 + muu1**2 + muu2**2)*muu3*HR(indj,i,tt,gen)) + 
     -            dconjg(FL(indj,j,tt,gen))*(-((-1 + muu1**2 + muu2**2 + 2*muu1**2*muu2**2 + muu3**2 + muu4**2 + 
     -                    2*muu3**2*muu4**2)*HL(indj,i,tt,gen)) + 
     -               (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred*HL(indj,i,tt,gen) - tred**2*HL(indj,i,tt,gen) + 
     -               (muu1**2 + muu2**2 + muu3**2 + muu4**2)*ured*HL(indj,i,tt,gen) - ured**2*HL(indj,i,tt,gen) - 
     -               2*muu2*muu3*(muu1**2 + muu4**2 - ured)*HR(indj,i,tt,gen)))) + 
     -      dGR(indi,i,tt,gen)*(-2*muu1*dconjg(ELc(indi,j,tt,gen))*
     -          (muu4*dconjg(FL(indj,j,tt,gen))*((muu2**2 + muu3**2 - ured)*HL(indj,i,tt,gen) + 2*muu2*muu3*HR(indj,i,tt,gen)) + 
     -            dconjg(FR(indj,j,tt,gen))*(-(muu2*(-1 + muu3**2 + muu4**2)*HL(indj,i,tt,gen)) + 
     -               muu3*(muu2**2 + muu4**2 - tred)*HR(indj,i,tt,gen))) + 
     -         dconjg(ER(indi,j,tt,gen))*(2*muu4*dconjg(FL(indj,j,tt,gen))*
     -             ((-1 + muu1**2 + muu2**2)*muu3*HL(indj,i,tt,gen) - muu2*(muu1**2 + muu3**2 - tred)*HR(indj,i,tt,gen)) + 
     -            dconjg(FR(indj,j,tt,gen))*(-2*muu2*muu3*(muu1**2 + muu4**2 - ured)*HL(indj,i,tt,gen) - 
     -               (-1 + muu1**2 + muu2**2 + 2*muu1**2*muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2)*HR(indj,i,tt,gen) + 
     -               (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred*HR(indj,i,tt,gen) - tred**2*HR(indj,i,tt,gen) + 
     -               (muu1**2 + muu2**2 + muu3**2 + muu4**2)*ured*HR(indj,i,tt,gen) - ured**2*HR(indj,i,tt,gen))))))/
     -  ((-ured + MSf(i,itypeU,gen)**2/sqrtS**2)*(-tred + MSf(j,itypeT,gen)**2/sqrtS**2))
	
	tempt = tempt + ChiChi2QQ_M2VertexTUCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexTUCounter

        enddo
	enddo
	
cc	write(*,*)'ChiChi U Channel Counterterm mit Z = ', tempZ
cc	write(*,*)'ChiChi U Channel Counterterm mit H = ', tempH
cc	write(*,*)'ChiChi U Channel Counterterm mit t = ', tempt
cc	write(*,*)'ChiChi U Channel Counterterm mit u = ', tempu
	tempZ=0d0
	tempH=0d0
	tempt=0d0
	tempu=0d0
	
	do i=1,2
	
!*************************************************
!* T-channel vertex counterterms (old U-channel) *
!*************************************************

	
	if (icha.eq.2) then
	
! (p2k2)-vertex with photon

!neue Version

	 	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (4*ApL*BpL*NC*sqrtS**2*(dFL(indj,i,tt,gen)*(2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*ELc(indi,i,tt,gen) + 
     -         (muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)*ER(indi,i,tt,gen)) + 
     -      dFR(indj,i,tt,gen)*((muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)*ELc(indi,i,tt,gen) + 
     -         2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*ER(indi,i,tt,gen))))/(sqrtS**2*tred - MSf(i,itypeT,gen)**2)

!alte Version
		
! 	 	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (-4*ApL*BpL*NC*(dFL(indj,i,tt,gen)*(2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*ELc(indi,i,tt,gen) + 
!      -         (muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*tred) + (muu3**2 - tred)**2)*ER(indi,i,tt,gen)) + 
!      -      ((muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*tred) + (muu3**2 - tred)**2)*ELc(indi,i,tt,gen) + 
!      -         2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*ER(indi,i,tt,gen))*FR(indj,i,tt,gen)))/
!      -  (-tred + MSf(i,itypeT,gen)**2/sqrtS**2)

! (p1k1)-vertex with photon

!neue Version

               ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (4*ApL*BpL*NC*sqrtS**2*(dELc(indi,i,tt,gen)*(2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*FL(indj,i,tt,gen) + 
     -         (muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)*FR(indj,i,tt,gen)) + 
     -      dER(indi,i,tt,gen)*((muu1**4 + muu3**2 + muu3**4 + muu1**2*(1 + 2*muu3**2 - 2*tred) - 2*muu3**2*tred + tred**2)*FL(indj,i,tt,gen) + 
     -         2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*FR(indj,i,tt,gen))))/(sqrtS**2*tred - MSf(i,itypeT,gen)**2)

!alte Version

! 		ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (-4*ApL*BpL*NC*(dELc(indi,i,tt,gen)*(2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*FL(indj,i,tt,gen) + 
!      -         (muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*tred) + (muu3**2 - tred)**2)*FR(indj,i,tt,gen)) + 
!      -      dER(indi,i,tt,gen)*((muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*tred) + (muu3**2 - tred)**2)*FL(indj,i,tt,gen) + 
!      -         2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*FR(indj,i,tt,gen))))/(-tred + MSf(i,itypeT,gen)**2/sqrtS**2)

	endif

! (p2k2)-vertex with Z

	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (-4*NC*sqrtS**2*(dconjg(AL)*(BR*((muu1**2 + muu3**2 - tred)*ELc(indi,i,tt,gen) + 2*muu1*muu3*ER(indi,i,tt,gen))*
     -          (2*muu2*muu4*dFL(indj,i,tt,gen) + (muu2**2 + muu4**2 - tred)*dFR(indj,i,tt,gen)) + 
     -         BL*(muu2*dFL(indj,i,tt,gen)*(muu3*(muu1**2 + muu4**2 - ured)*ELc(indi,i,tt,gen) - 
     -               muu1*(-1 + muu3**2 + muu4**2)*ER(indi,i,tt,gen)) + 
     -            muu4*(-((-1 + muu1**2 + muu2**2)*muu3*ELc(indi,i,tt,gen)) + muu1*(muu2**2 + muu3**2 - ured)*ER(indi,i,tt,gen))*
     -             dFR(indj,i,tt,gen))) + dconjg(AR)*(BL*
     -          (2*muu1*muu3*ELc(indi,i,tt,gen) + (muu1**2 + muu3**2 - tred)*ER(indi,i,tt,gen))*
     -          ((muu2**2 + muu4**2 - tred)*dFL(indj,i,tt,gen) + 2*muu2*muu4*dFR(indj,i,tt,gen)) + 
     -         BR*(muu1*ELc(indi,i,tt,gen)*(muu4*(muu2**2 + muu3**2 - ured)*dFL(indj,i,tt,gen) - 
     -               muu2*(-1 + muu3**2 + muu4**2)*dFR(indj,i,tt,gen)) + 
     -            muu3*ER(indi,i,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu4*dFL(indj,i,tt,gen)) + 
     -               muu2*(muu1**2 + muu4**2 - ured)*dFR(indj,i,tt,gen))))))/(DenSmV1*(-tred + MSf(i,itypeT,gen)**2/sqrtS**2))

! (p1k1)-vertex with Z

	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (-4*NC*sqrtS**2*(dconjg(AL)*(BR*((muu1**2 + muu3**2 - tred)*dELc(indi,i,tt,gen) + 2*muu1*muu3*dER(indi,i,tt,gen))*
     -          (2*muu2*muu4*FL(indj,i,tt,gen) + (muu2**2 + muu4**2 - tred)*FR(indj,i,tt,gen)) + 
     -         BL*(muu3*dELc(indi,i,tt,gen)*(muu2*(muu1**2 + muu4**2 - ured)*FL(indj,i,tt,gen) - 
     -               (-1 + muu1**2 + muu2**2)*muu4*FR(indj,i,tt,gen)) + 
     -            muu1*dER(indi,i,tt,gen)*(-(muu2*(-1 + muu3**2 + muu4**2)*FL(indj,i,tt,gen)) + 
     -               muu4*(muu2**2 + muu3**2 - ured)*FR(indj,i,tt,gen)))) + 
     -      dconjg(AR)*(BL*(2*muu1*muu3*dELc(indi,i,tt,gen) + (muu1**2 + muu3**2 - tred)*dER(indi,i,tt,gen))*
     -          ((muu2**2 + muu4**2 - tred)*FL(indj,i,tt,gen) + 2*muu2*muu4*FR(indj,i,tt,gen)) + 
     -         BR*(muu1*dELc(indi,i,tt,gen)*(muu4*(muu2**2 + muu3**2 - ured)*FL(indj,i,tt,gen) - 
     -               muu2*(-1 + muu3**2 + muu4**2)*FR(indj,i,tt,gen)) + 
     -            muu3*dER(indi,i,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu4*FL(indj,i,tt,gen)) + 
     -               muu2*(muu1**2 + muu4**2 - ured)*FR(indj,i,tt,gen))))))/(DenSmV1*(-tred + MSf(i,itypeT,gen)**2/sqrtS**2))
	
	tempZ = tempZ + ChiChi2QQ_M2VertexTUCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexTUCounter
	
! (p2k2)-vertex with Higgs

	do j= jmin,jmax
	
	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter -((NC*sqrtS**2*(dconjg(CL(j))*(-2*muu2*dFR(indj,i,tt,gen)*
     -            (muu4*dconjg(DL(j))*((muu1**2 + muu3**2 - tred)*ELc(indi,i,tt,gen) + 2*muu1*muu3*ER(indi,i,tt,gen)) + 
     -              dconjg(DR(j))*(-(muu3*(muu1**2 + muu4**2 - ured)*ELc(indi,i,tt,gen)) + 
     -                 muu1*(-1 + muu3**2 + muu4**2)*ER(indi,i,tt,gen))) + 
     -           dFL(indj,i,tt,gen)*(dconjg(DL(j))*((-1 + muu3**2 + muu2**2*(1 - 2*muu3**2) + muu4**2 + 
     -                    (muu2**2 + muu3**2 + muu4**2 - tred - ured)*(tred - ured) + muu1**2*(1 - 2*muu4**2 + tred - ured))*
     -                  ELc(indi,i,tt,gen) - 2*muu1*muu3*(muu2**2 + muu4**2 - tred)*ER(indi,i,tt,gen)) + 
     -              2*muu4*dconjg(DR(j))*(-((-1 + muu1**2 + muu2**2)*muu3*ELc(indi,i,tt,gen)) + 
     -                 muu1*(muu2**2 + muu3**2 - ured)*ER(indi,i,tt,gen)))) + 
     -        dconjg(CR(j))*(dFR(indj,i,tt,gen)*(2*muu4*dconjg(DL(j))*
     -               (muu1*(muu2**2 + muu3**2 - ured)*ELc(indi,i,tt,gen) - (-1 + muu1**2 + muu2**2)*muu3*ER(indi,i,tt,gen)) + 
     -              dconjg(DR(j))*(-2*muu1*muu3*(muu2**2 + muu4**2 - tred)*ELc(indi,i,tt,gen) + 
     -                 (-1 + muu3**2 + muu2**2*(1 - 2*muu3**2) + muu4**2 + 
     -                    (muu2**2 + muu3**2 + muu4**2 - tred - ured)*(tred - ured) + muu1**2*(1 - 2*muu4**2 + tred - ured))*
     -                  ER(indi,i,tt,gen))) - 2*muu2*dFL(indj,i,tt,gen)*
     -            (muu4*dconjg(DR(j))*(2*muu1*muu3*ELc(indi,i,tt,gen) + (muu1**2 + muu3**2 - tred)*ER(indi,i,tt,gen)) + 
     -              dconjg(DL(j))*(muu1*(-1 + muu3**2 + muu4**2)*ELc(indi,i,tt,gen) - 
     -                 muu3*(muu1**2 + muu4**2 - ured)*ER(indi,i,tt,gen))))))/
     -    (DenSmH1(j)*(-tred + MSf(i,itypeT,gen)**2/sqrtS**2)))

! (p1k1)-vertex with Higgs

	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter -((NC*sqrtS**2*(dconjg(CR(j))*(dER(indi,i,tt,gen)*
     -            (2*muu3*dconjg(DL(j))*(muu2*(muu1**2 + muu4**2 - ured)*FL(indj,i,tt,gen) - 
     -                 (-1 + muu1**2 + muu2**2)*muu4*FR(indj,i,tt,gen)) + 
     -              dconjg(DR(j))*(-2*muu2*muu4*(muu1**2 + muu3**2 - tred)*FL(indj,i,tt,gen) + 
     -                 (-1 + muu3**2 + muu2**2*(1 - 2*muu3**2) + muu4**2 + 
     -                    (muu2**2 + muu3**2 + muu4**2 - tred - ured)*(tred - ured) + muu1**2*(1 - 2*muu4**2 + tred - ured))*
     -                  FR(indj,i,tt,gen))) - 2*muu1*dELc(indi,i,tt,gen)*
     -            (muu3*dconjg(DR(j))*(2*muu2*muu4*FL(indj,i,tt,gen) + (muu2**2 + muu4**2 - tred)*FR(indj,i,tt,gen)) + 
     -              dconjg(DL(j))*(muu2*(-1 + muu3**2 + muu4**2)*FL(indj,i,tt,gen) - 
     -                 muu4*(muu2**2 + muu3**2 - ured)*FR(indj,i,tt,gen)))) + 
     -        dconjg(CL(j))*(-2*muu1*dER(indi,i,tt,gen)*
     -            (muu3*dconjg(DL(j))*((muu2**2 + muu4**2 - tred)*FL(indj,i,tt,gen) + 2*muu2*muu4*FR(indj,i,tt,gen)) + 
     -              dconjg(DR(j))*(-(muu4*(muu2**2 + muu3**2 - ured)*FL(indj,i,tt,gen)) + 
     -                 muu2*(-1 + muu3**2 + muu4**2)*FR(indj,i,tt,gen))) + 
     -           dELc(indi,i,tt,gen)*(dconjg(DL(j))*((-1 + muu3**2 + muu2**2*(1 - 2*muu3**2) + muu4**2 + 
     -                    (muu2**2 + muu3**2 + muu4**2 - tred - ured)*(tred - ured) + muu1**2*(1 - 2*muu4**2 + tred - ured))*
     -                  FL(indj,i,tt,gen) - 2*muu2*muu4*(muu1**2 + muu3**2 - tred)*FR(indj,i,tt,gen)) + 
     -              2*muu3*dconjg(DR(j))*(-((-1 + muu1**2 + muu2**2)*muu4*FL(indj,i,tt,gen)) + 
     -                 muu2*(muu1**2 + muu4**2 - ured)*FR(indj,i,tt,gen))))))/
     -    (DenSmH1(j)*(-tred + MSf(i,itypeT,gen)**2/sqrtS**2)))

	enddo
	
	tempH = tempH + ChiChi2QQ_M2VertexTUCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexTUCounter

! (p2k2)-vertex with U-channel (old T-channel)

	do j=1,2

	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (NC*(dFL(indj,i,tt,gen)*((-2*muu1*muu4*dconjg(GR(indi,j,tt,gen))*
     -             ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,j,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,j,tt,gen))) + 
     -            dconjg(GL(indi,j,tt,gen))*(-((-1 + muu1**2 + muu2**2 + 2*muu1**2*muu2**2 + muu3**2 + muu4**2 + 
     -                    2*muu3**2*muu4**2)*dconjg(HL(indj,j,tt,gen))) + 
     -               (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred*dconjg(HL(indj,j,tt,gen)) - tred**2*dconjg(HL(indj,j,tt,gen)) + 
     -               (muu1**2 + muu2**2 + muu3**2 + muu4**2)*ured*dconjg(HL(indj,j,tt,gen)) - ured**2*dconjg(HL(indj,j,tt,gen)) - 
     -               2*muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,j,tt,gen))))*ELc(indi,i,tt,gen) + 
     -         2*(muu1*dconjg(GL(indi,j,tt,gen))*(-(muu3*(muu2**2 + muu4**2 - tred)*dconjg(HL(indj,j,tt,gen))) + 
     -               muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,j,tt,gen))) + 
     -            muu4*dconjg(GR(indi,j,tt,gen))*((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,j,tt,gen)) - 
     -               muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,j,tt,gen))))*ER(indi,i,tt,gen)) + 
     -      dFR(indj,i,tt,gen)*(2*(muu4*dconjg(GL(indi,j,tt,gen))*
     -             (-(muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,j,tt,gen))) + 
     -               (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,j,tt,gen))) + 
     -            muu1*dconjg(GR(indi,j,tt,gen))*(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,j,tt,gen)) - 
     -               muu3*(muu2**2 + muu4**2 - tred)*dconjg(HR(indj,j,tt,gen))))*ELc(indi,i,tt,gen) + 
     -         (-2*muu1*muu4*dconjg(GL(indi,j,tt,gen))*
     -             (2*muu2*muu3*dconjg(HL(indj,j,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,j,tt,gen))) + 
     -            dconjg(GR(indi,j,tt,gen))*(-2*muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,j,tt,gen)) - 
     -               (-1 + muu1**2 + muu2**2 + 2*muu1**2*muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2)*
     -                dconjg(HR(indj,j,tt,gen)) + (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred*dconjg(HR(indj,j,tt,gen)) - 
     -               tred**2*dconjg(HR(indj,j,tt,gen)) + (muu1**2 + muu2**2 + muu3**2 + muu4**2)*ured*dconjg(HR(indj,j,tt,gen)) - 
     -               ured**2*dconjg(HR(indj,j,tt,gen))))*ER(indi,i,tt,gen))))/
     -  ((-tred + MSf(i,itypeT,gen)**2/sqrtS**2)*(-ured + MSf(j,itypeU,gen)**2/sqrtS**2))

! (p1k1)-vertex with U-channel (old T-channel)

	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (NC*(dELc(indi,i,tt,gen)*((-2*muu1*muu4*dconjg(GR(indi,j,tt,gen))*
     -             ((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,j,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,j,tt,gen))) + 
     -            dconjg(GL(indi,j,tt,gen))*(-((-1 + muu1**2 + muu2**2 + 2*muu1**2*muu2**2 + muu3**2 + muu4**2 + 
     -                    2*muu3**2*muu4**2)*dconjg(HL(indj,j,tt,gen))) + 
     -               (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred*dconjg(HL(indj,j,tt,gen)) - tred**2*dconjg(HL(indj,j,tt,gen)) + 
     -               (muu1**2 + muu2**2 + muu3**2 + muu4**2)*ured*dconjg(HL(indj,j,tt,gen)) - ured**2*dconjg(HL(indj,j,tt,gen)) - 
     -               2*muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HR(indj,j,tt,gen))))*FL(indj,i,tt,gen) + 
     -         2*(muu4*dconjg(GL(indi,j,tt,gen))*(-(muu2*(muu1**2 + muu3**2 - tred)*dconjg(HL(indj,j,tt,gen))) + 
     -               (-1 + muu1**2 + muu2**2)*muu3*dconjg(HR(indj,j,tt,gen))) + 
     -            muu1*dconjg(GR(indi,j,tt,gen))*(muu2*(-1 + muu3**2 + muu4**2)*dconjg(HL(indj,j,tt,gen)) - 
     -               muu3*(muu2**2 + muu4**2 - tred)*dconjg(HR(indj,j,tt,gen))))*FR(indj,i,tt,gen)) + 
     -      dER(indi,i,tt,gen)*(2*(muu1*dconjg(GL(indi,j,tt,gen))*
     -             (-(muu3*(muu2**2 + muu4**2 - tred)*dconjg(HL(indj,j,tt,gen))) + 
     -               muu2*(-1 + muu3**2 + muu4**2)*dconjg(HR(indj,j,tt,gen))) + 
     -            muu4*dconjg(GR(indi,j,tt,gen))*((-1 + muu1**2 + muu2**2)*muu3*dconjg(HL(indj,j,tt,gen)) - 
     -               muu2*(muu1**2 + muu3**2 - tred)*dconjg(HR(indj,j,tt,gen))))*FL(indj,i,tt,gen) + 
     -         (-2*muu1*muu4*dconjg(GL(indi,j,tt,gen))*
     -             (2*muu2*muu3*dconjg(HL(indj,j,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,j,tt,gen))) + 
     -            dconjg(GR(indi,j,tt,gen))*(-2*muu2*muu3*(muu1**2 + muu4**2 - ured)*dconjg(HL(indj,j,tt,gen)) - 
     -               (-1 + muu1**2 + muu2**2 + 2*muu1**2*muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2)*
     -                dconjg(HR(indj,j,tt,gen)) + (muu1**2 + muu2**2 + muu3**2 + muu4**2)*tred*dconjg(HR(indj,j,tt,gen)) - 
     -               tred**2*dconjg(HR(indj,j,tt,gen)) + (muu1**2 + muu2**2 + muu3**2 + muu4**2)*ured*dconjg(HR(indj,j,tt,gen)) - 
     -               ured**2*dconjg(HR(indj,j,tt,gen))))*FR(indj,i,tt,gen))))/
     -  ((-tred + MSf(i,itypeT,gen)**2/sqrtS**2)*(-ured + MSf(j,itypeU,gen)**2/sqrtS**2))

	tempu = tempu + ChiChi2QQ_M2VertexTUCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexTUCounter
	
! (p2k2)-vertex with T-channel (old U-channel)

	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (2*NC*(((muu2**2 + muu4**2 - tred)*dconjg(FL(indj,j,tt,gen)) + 2*muu2*muu4*dconjg(FR(indj,j,tt,gen)))*
     -       dFL(indj,i,tt,gen) + (2*muu2*muu4*dconjg(FL(indj,j,tt,gen)) + (muu2**2 + muu4**2 - tred)*dconjg(FR(indj,j,tt,gen)))*
     -       dFR(indj,i,tt,gen))*(((muu1**2 + muu3**2 - tred)*dconjg(ELc(indi,j,tt,gen)) + 2*muu1*muu3*dconjg(ER(indi,j,tt,gen)))*
     -       ELc(indi,i,tt,gen) + (2*muu1*muu3*dconjg(ELc(indi,j,tt,gen)) + (muu1**2 + muu3**2 - tred)*dconjg(ER(indi,j,tt,gen)))*
     -       ER(indi,i,tt,gen)))/((-tred + MSf(i,itypeT,gen)**2/sqrtS**2)*(-tred + MSf(j,itypeT,gen)**2/sqrtS**2))

! (p1k1)-vertex with T-channel (old U-channel)

	ChiChi2QQ_M2VertexTUCounter = ChiChi2QQ_M2VertexTUCounter + (2*NC*(((muu1**2 + muu3**2 - tred)*dconjg(ELc(indi,j,tt,gen)) + 2*muu1*muu3*dconjg(ER(indi,j,tt,gen)))*
     -       dELc(indi,i,tt,gen) + (2*muu1*muu3*dconjg(ELc(indi,j,tt,gen)) + 
     -         (muu1**2 + muu3**2 - tred)*dconjg(ER(indi,j,tt,gen)))*dER(indi,i,tt,gen))*
     -    (((muu2**2 + muu4**2 - tred)*dconjg(FL(indj,j,tt,gen)) + 2*muu2*muu4*dconjg(FR(indj,j,tt,gen)))*FL(indj,i,tt,gen) + 
     -      (2*muu2*muu4*dconjg(FL(indj,j,tt,gen)) + (muu2**2 + muu4**2 - tred)*dconjg(FR(indj,j,tt,gen)))*FR(indj,i,tt,gen)))/
     -  ((-tred + MSf(i,itypeT,gen)**2/sqrtS**2)*(-tred + MSf(j,itypeT,gen)**2/sqrtS**2))
	
	tempt = tempt + ChiChi2QQ_M2VertexTUCounter - tempCounter
	tempCounter = ChiChi2QQ_M2VertexTUCounter
	
	
	  enddo
	  
	enddo
		
cc	write(*,*)'ChiChi T Channel Counterterm mit Z = ', tempZ
cc	write(*,*)'ChiChi T Channel Counterterm mit H = ', tempH
cc	write(*,*)'ChiChi T Channel Counterterm mit t = ', tempt
cc	write(*,*)'ChiChi T Channel Counterterm mit u = ', tempu
	
	end

!************************************************

	function ChiChi2QQ_M2PropagatorCounter()

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"
#include "ChiChi2QQ_GenCounterterms.h"

	integer i,j,k,l,m,tcount,imin,imax,jmin,jmax,itypeT,itypeU,att,antitt

	complex*16 DenSmV, DenSmH(6)
	double precision DenSmV1, DenSmH1(6), tempCounter, tempH, tempZ, tempt, tempu

	double precision Krondelta
	double precision ChiChi2QQ_M2PropagatorCounter
	external ChiChi2QQ_CalcProp

	ChiChi2QQ_M2PropagatorCounter = 0d0
	tempCounter = 0d0
	tempH = 0d0
	tempZ = 0d0
	tempu = 0d0
	tempt = 0d0

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

!*****************************************************
!*********** U-Loop counterterms (old T-channel) *****
!*****************************************************

	 call ChiChi2QQ_CalcProp(sqrtS**2*ured,itypeU,gen)

	 do k=1,2
	 do l=1,2
	 
! U-channel prop. counterterm interference with photon tree-level

	if (icha.eq.2) then
		
		ChiChi2QQ_M2PropagatorCounter = ChiChi2QQ_M2PropagatorCounter + (-4*ApL*BpL*NC*(GL(indi,l,tt,gen)*(2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*HL(indj,k,tt,gen) + 
     -         (muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*ured) + (muu3**2 - ured)**2)*HR(indj,k,tt,gen)) + 
     -      GR(indi,l,tt,gen)*((muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*ured) + (muu3**2 - ured)**2)*HL(indj,k,tt,gen) + 
     -         2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - tred - 2*ured)*HR(indj,k,tt,gen)))*Prop(k,l,itypeU,gen))/
     -  (sqrtS**2*(-ured + MSf(k,itypeU,gen)**2/sqrtS**2)*(ured - MSf(l,itypeU,gen)**2/sqrtS**2))

	endif
	
! U-channel prop. counterterm interference with Z/W-boson tree-level

	ChiChi2QQ_M2PropagatorCounter = ChiChi2QQ_M2PropagatorCounter + (-4*NC*(dconjg(AR)*(BR*(2*muu1*muu4*GL(indi,l,tt,gen) 
     -			+ (muu1**2 + muu4**2 - ured)*GR(indi,l,tt,gen))*
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
     -               muu3*(muu2**2 + muu4**2 - tred)*HR(indj,k,tt,gen)))))*Prop(k,l,itypeU,gen))/
     -  (DenSmV1*(ured - MSf(k,itypeU,gen)**2/sqrtS**2)*(-ured + MSf(l,itypeU,gen)**2/sqrtS**2))
	 
	tempZ = tempZ + ChiChi2QQ_M2PropagatorCounter - tempCounter
	tempCounter = ChiChi2QQ_M2PropagatorCounter
	
! U-channel prop. counterterm interference with Higgs-boson tree-level

	 do j=jmin,jmax
		
	 ChiChi2QQ_M2PropagatorCounter = ChiChi2QQ_M2PropagatorCounter + ((NC*(dconjg(CL(j))*(2*dconjg(DR(j))*(muu1*GR(indi,l,tt,gen)*
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
     -                  HR(indj,k,tt,gen)))))*Prop(k,l,itypeU,gen))/
     -    (DenSmH1(j)*(-ured + MSf(k,itypeU,gen)**2/sqrtS**2)*(ured - MSf(l,itypeU,gen)**2/sqrtS**2)))

	 enddo

	tempH = tempH + ChiChi2QQ_M2PropagatorCounter - tempCounter
	tempCounter = ChiChi2QQ_M2PropagatorCounter
	
! U-channel prop. counterterm interference with U-channel tree-level (with old T-channel)

	 do j=1,2
		
	 ChiChi2QQ_M2PropagatorCounter = ChiChi2QQ_M2PropagatorCounter + (2*NC*(((muu1**2 + muu4**2 - ured)*dconjg(GL(indi,j,tt,gen)) + 2*muu1*muu4*dconjg(GR(indi,j,tt,gen)))*
     -       GL(indi,l,tt,gen) + (2*muu1*muu4*dconjg(GL(indi,j,tt,gen)) + (muu1**2 + muu4**2 - ured)*dconjg(GR(indi,j,tt,gen)))*
     -       GR(indi,l,tt,gen))*(((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,j,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,j,tt,gen)))*
     -       HL(indj,k,tt,gen) + (2*muu2*muu3*dconjg(HL(indj,j,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,j,tt,gen)))*
     -       HR(indj,k,tt,gen))*Prop(k,l,itypeU,gen))/
     -  (sqrtS**2*(ured - MSf(j,itypeU,gen)**2/sqrtS**2)*(-ured + MSf(k,itypeU,gen)**2/sqrtS**2)*
     -    (ured - MSf(l,itypeU,gen)**2/sqrtS**2))

	 enddo

	tempu = tempu + ChiChi2QQ_M2PropagatorCounter - tempCounter
	tempCounter = ChiChi2QQ_M2PropagatorCounter
	
! U-channel prop. counterterm interference with T-channel tree-level (with old U-channel)

	 do j=1,2
		
	 ChiChi2QQ_M2PropagatorCounter = ChiChi2QQ_M2PropagatorCounter + (NC*(dconjg(ELc(indi,j,tt,gen))*(2*dconjg(FR(indj,j,tt,gen))*
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
     -    Prop(k,l,itypeU,gen))/
     -  (sqrtS**2*(-tred + MSf(j,itypeT,gen)**2/sqrtS**2)*(-ured + MSf(k,itypeU,gen)**2/sqrtS**2)*
     -    (ured - MSf(l,itypeU,gen)**2/sqrtS**2))
	 
	 enddo

	tempt = tempt + ChiChi2QQ_M2PropagatorCounter - tempCounter
	tempCounter = ChiChi2QQ_M2PropagatorCounter
	
	 enddo
	 enddo

cc	write(*,*)'ChiChi PropagatorCounterterm UChannel mit Z = ', tempZ
cc	write(*,*)'ChiChi PropagatorCounterterm UChannel mit H = ', tempH
cc	write(*,*)'ChiChi PropagatorCounterterm UChannel mit t = ', tempt
cc	write(*,*)'ChiChi PropagatorCounterterm UChannel mit u = ', tempu
	tempZ=0d0
	tempH=0d0
	tempt=0d0
	tempu=0d0
	
!*****************************************************
!*********** T-Loop counterterms (old U-channel) *****
!*****************************************************

	 call ChiChi2QQ_CalcProp(sqrtS**2*tred,itypeT,gen)

	 do k=1,2
	 do l=1,2

! T-channel prop. counterterm interference with photon tree-level

	 if (icha.eq.2) then

	 	ChiChi2QQ_M2PropagatorCounter = ChiChi2QQ_M2PropagatorCounter + (4*ApL*BpL*NC*(ELc(indi,k,tt,gen)*(2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*FL(indj,l,tt,gen) + 
     -         (muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*tred) + (muu3**2 - tred)**2)*FR(indj,l,tt,gen)) + 
     -      ER(indi,k,tt,gen)*((muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*tred) + (muu3**2 - tred)**2)*FL(indj,l,tt,gen) + 
     -         2*muu1*muu3*(3*muu1**2 + 3*muu3**2 - 2*tred - ured)*FR(indj,l,tt,gen)))*Prop(k,l,itypeT,gen))/
     -  (sqrtS**2*(-tred + MSf(k,itypeT,gen)**2/sqrtS**2)*(tred - MSf(l,itypeT,gen)**2/sqrtS**2))

	 endif

! T-channel prop. counterterm interference with Z/W-boson tree-level

	 ChiChi2QQ_M2PropagatorCounter = ChiChi2QQ_M2PropagatorCounter + (4*NC*(dconjg(AL)*(BR*((muu1**2 + muu3**2 - tred)*ELc(indi,k,tt,gen) + 2*muu1*muu3*ER(indi,k,tt,gen))*
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
     -               muu2*(muu1**2 + muu4**2 - ured)*FR(indj,l,tt,gen)))))*Prop(k,l,itypeT,gen))/
     -  (DenSmV1*(tred - MSf(k,itypeT,gen)**2/sqrtS**2)*(-tred + MSf(l,itypeT,gen)**2/sqrtS**2))	

	tempZ = tempZ + ChiChi2QQ_M2PropagatorCounter - tempCounter
	tempCounter = ChiChi2QQ_M2PropagatorCounter
	
! T-channel prop. counterterm interference with Higgs-boson tree-level

	 do j=jmin,jmax

	 ChiChi2QQ_M2PropagatorCounter = ChiChi2QQ_M2PropagatorCounter + ((NC*(dconjg(CR(j))*(dconjg(DR(j))*(-2*muu1*muu3*ELc(indi,k,tt,gen)*
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
     -                 muu2*(muu1**2 + muu4**2 - ured)*FR(indj,l,tt,gen)))))*Prop(k,l,itypeT,gen))/
     -    (DenSmH1(j)*(-tred + MSf(k,itypeT,gen)**2/sqrtS**2)*(tred - MSf(l,itypeT,gen)**2/sqrtS**2)))		

	 enddo
	
	tempH = tempH + ChiChi2QQ_M2PropagatorCounter - tempCounter
	tempCounter = ChiChi2QQ_M2PropagatorCounter
	
! T-channel prop. counterterm interference with U-channel tree-level (old T-channel)

	 do j=1,2
		
	 ChiChi2QQ_M2PropagatorCounter = ChiChi2QQ_M2PropagatorCounter - (NC*(ELc(indi,k,tt,gen)*((2*muu1*muu4*dconjg(GR(indi,j,tt,gen))*
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
     -                dconjg(HR(indj,j,tt,gen))))*FR(indj,l,tt,gen)))*Prop(k,l,itypeT,gen))/
     -  (sqrtS**2*(-ured + MSf(j,itypeU,gen)**2/sqrtS**2)*(-tred + MSf(k,itypeT,gen)**2/sqrtS**2)*
     -    (-tred + MSf(l,itypeT,gen)**2/sqrtS**2))
	 enddo
	
	tempu = tempu + ChiChi2QQ_M2PropagatorCounter - tempCounter
	tempCounter = ChiChi2QQ_M2PropagatorCounter
	
! T-channel prop. counterterm interference with T-channel tree-level (old U-channel)

	 do j=1,2
		
	 ChiChi2QQ_M2PropagatorCounter = ChiChi2QQ_M2PropagatorCounter + (2*NC*(((muu1**2 + muu3**2 - tred)*dconjg(ELc(indi,j,tt,gen)) + 2*muu1*muu3*dconjg(ER(indi,j,tt,gen)))*
     -       ELc(indi,k,tt,gen) + (2*muu1*muu3*dconjg(ELc(indi,j,tt,gen)) + (muu1**2 + muu3**2 - tred)*dconjg(ER(indi,j,tt,gen)))*
     -       ER(indi,k,tt,gen))*(((muu2**2 + muu4**2 - tred)*dconjg(FL(indj,j,tt,gen)) + 2*muu2*muu4*dconjg(FR(indj,j,tt,gen)))*
     -       FL(indj,l,tt,gen) + (2*muu2*muu4*dconjg(FL(indj,j,tt,gen)) + (muu2**2 + muu4**2 - tred)*dconjg(FR(indj,j,tt,gen)))*
     -       FR(indj,l,tt,gen))*Prop(k,l,itypeT,gen))/
     -  (sqrtS**2*(-tred + MSf(j,itypeT,gen)**2/sqrtS**2)*(tred - MSf(k,itypeT,gen)**2/sqrtS**2)*
     -    (tred - MSf(l,itypeT,gen)**2/sqrtS**2))
	 enddo
	
	tempt = tempt + ChiChi2QQ_M2PropagatorCounter - tempCounter
	tempCounter = ChiChi2QQ_M2PropagatorCounter
	
	 enddo
	 enddo
	
cc	write(*,*)'ChiChi PropagatorCounterterm TChannel mit Z = ', tempZ
cc	write(*,*)'ChiChi PropagatorCounterterm TChannel mit H = ', tempH
cc	write(*,*)'ChiChi PropagatorCounterterm TChannel mit t = ', tempt
cc	write(*,*)'ChiChi PropagatorCounterterm TChannel mit u = ', tempu
	
	 end
