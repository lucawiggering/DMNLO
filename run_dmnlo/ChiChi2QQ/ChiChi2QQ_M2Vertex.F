************************************************

	function ChiChi2QQ_M2Vertex()

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"
#include "ChiChi2QQ_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer i,j,k,l,tcount,imin,imax,jmin,jmax,itypeU,itypeT,att,antitt
	complex*16 g0,g0L,g0R,g1L,g1R,g2L,g2R,g3L,g3R,g4L,g4R
	double precision M0
	double precision Ml0,Ml1,tempCounter, tempH, tempZ, tempt, tempu
	double precision tempH11, tempZ11, tempt11, tempu11, tempH12, tempZ12, tempt12, tempu12
	double precision tempH21, tempZ21, tempt21, tempu21, tempH22, tempZ22, tempt22, tempu22
	complex*16 BL1,BL2,BL3,BR1,BR2,BR3,BpL1,BpL2,BpL3,BpR1,BpR2,BpR3,DL1,DR1,EL1,ER1,FL1,FR1,GL1,GR1,HL1,HR1
	double precision swap, test
	
	complex*16 BL1g,BL2g,BL3g,BR1g,BR2g,BR3g
	
	double complex nc0,nc1,nc2,nc00,nc11,nc12,nc22

	complex*16 DenSmV, DenSmH(6)
	double precision DenSmV1, DenSmH1(6)
	
	complex*16 tempGluonL, tempGluonR
	
* functions & subroutines
	double precision ChiChi2QQ_M2Vertex
	complex*16 B02,C03f
	external myCget, myBget

	ChiChi2QQ_M2Vertex = 0d0
	tempCounter = 0d0
	tempH = 0d0
	tempZ = 0d0
	tempu = 0d0
	tempt = 0d0
	
	tempH11 = 0d0
	tempZ11 = 0d0
	tempu11 = 0d0
	tempt11 = 0d0
	tempH12 = 0d0
	tempZ12 = 0d0
	tempu12 = 0d0
	tempt12 = 0d0
	tempH21 = 0d0
	tempZ21 = 0d0
	tempu21 = 0d0
	tempt21 = 0d0
	tempH22 = 0d0
	tempZ22 = 0d0
	tempu22 = 0d0
	tempt22 = 0d0

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

! Higgs contributions

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

! widths of vector-boson and Higgs included in propagators
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


!****************************************************
!************ Loop with Vector boson (photon, Z, W) *
!****************************************************

! Gluon loop 

! Gluon correction to Z/W diagram		 

	call myCget(muu4**2*sqrtS**2,sqrtS**2,muu3**2*sqrtS**2,0d0,muu4**2*sqrtS**2,muu3**2*sqrtS**2) 
	call myBget(sqrtS**2,muu4**2*sqrtS**2,muu3**2*sqrtS**2)
	
	g0L = BL
	g0R = BR
	g1L = -GS
	g1R = -GS
	g2L = -GS
	g2R = -GS
	
	Ml0 = 0d0
	Ml1 = muu4*sqrtS
	Ml2 = muu3*sqrtS
	
	BL1 = (g0L*g1L*g2L*(xb0 + Ml0**2*xc0 - sqrtS**2*xc0 + muu3**2*sqrtS**2*xc0 + muu4**2*sqrtS**2*xc0 - 2*xc00 - 
     &       sqrtS**2*xc1 + muu3**2*sqrtS**2*xc1 + 2*muu4**2*sqrtS**2*xc1 - sqrtS**2*xc2 + 2*muu3**2*sqrtS**2*xc2 + 
     &       muu4**2*sqrtS**2*xc2) + g0R*(-(g1L*g2L*Ml1*Ml2*xc0) + g1R*g2R*muu3*muu4*sqrtS**2*(xc0 + xc1 + xc2)))/
     &  (8.*Pi**2)

	BR1 = (g0R*g1R*g2R*(xb0 + Ml0**2*xc0 - sqrtS**2*xc0 + muu3**2*sqrtS**2*xc0 + muu4**2*sqrtS**2*xc0 - 2*xc00 - 
     &       sqrtS**2*xc1 + muu3**2*sqrtS**2*xc1 + 2*muu4**2*sqrtS**2*xc1 - sqrtS**2*xc2 + 2*muu3**2*sqrtS**2*xc2 + 
     &       muu4**2*sqrtS**2*xc2) + g0L*(-(g1R*g2R*Ml1*Ml2*xc0) + g1L*g2L*muu3*muu4*sqrtS**2*(xc0 + xc1 + xc2)))/
     &  (8.*Pi**2)

	BL2 = -(g0R*g1R*g2R*muu4*sqrtS*(xc0 + xc1 + xc12 + xc2) - g1L*g2R*(g0L*Ml2*xc2 + g0R*Ml1*(xc0 + xc2)) + 
     &     g0L*g1L*g2L*muu3*sqrtS*(xc2 + xc22))/(4.*Pi**2)
	
	BR2 = -(g0L*g1L*g2L*muu4*sqrtS*(xc0 + xc1 + xc12 + xc2) - g1R*g2L*(g0R*Ml2*xc2 + g0L*Ml1*(xc0 + xc2)) + 
     &     g0R*g1R*g2R*muu3*sqrtS*(xc2 + xc22))/(4.*Pi**2)

	BL3 = (-(g1L*g2R*(g0R*Ml1*xc1 + g0L*Ml2*(xc0 + xc1))) + g0R*g1R*g2R*muu4*sqrtS*(xc1 + xc11) + 
     &    g0L*g1L*g2L*muu3*sqrtS*(xc0 + xc1 + xc12 + xc2))/(4.*Pi**2)

	BR3 = (-(g1R*g2L*(g0L*Ml1*xc1 + g0R*Ml2*(xc0 + xc1))) + g0L*g1L*g2L*muu4*sqrtS*(xc1 + xc11) + 
     &    g0R*g1R*g2R*muu3*sqrtS*(xc0 + xc1 + xc12 + xc2))/(4.*Pi**2)
	
cc	write(*,*)' --- LOOP GLUON --- '
cc	write(*,*)'BL1 =',BL1
cc	write(*,*)'BR1 =',BR1
cc	write(*,*)'BL2 =',BL2
cc	write(*,*)'BR2 =',BR2
cc	write(*,*)'BL3 =',BL3
cc	write(*,*)'BR3 =',BR3
cc
cc      write(*,*)' --- Cget gluon--- '
cc      write(*,*)'xc0 =',xc0
cc      write(*,*)'xc1 =',xc1
cc      write(*,*)'xc2 =',xc2
cc      write(*,*)'xc00 =',xc00
cc      write(*,*)'xc11 =',xc11
cc      write(*,*)'xc12 =',xc12
cc      write(*,*)'xc22 =',xc22

	
! Gluon correction to photon diagram
	if (icha.eq.2) then
		g0L = BpL
		g0R = BpL
		g1L = -GS
		g1R = -GS
		g2L = -GS
		g2R = -GS
      	
		Ml0 = 0d0
		Ml1 = muu4*sqrtS
		Ml2 = muu3*sqrtS
      	
		BpL1 = (g0L*g1L*g2L*(xb0 + Ml0**2*xc0 - sqrtS**2*xc0 + muu3**2*sqrtS**2*xc0 + muu4**2*sqrtS**2*xc0 - 2*xc00 - 
     &	       sqrtS**2*xc1 + muu3**2*sqrtS**2*xc1 + 2*muu4**2*sqrtS**2*xc1 - sqrtS**2*xc2 + 2*muu3**2*sqrtS**2*xc2 + 
     &	       muu4**2*sqrtS**2*xc2) + g0R*(-(g1L*g2L*Ml1*Ml2*xc0) + g1R*g2R*muu3*muu4*sqrtS**2*(xc0 + xc1 + xc2)))/
     &	  (8.*Pi**2)
      	
		BpR1 = (g0R*g1R*g2R*(xb0 + Ml0**2*xc0 - sqrtS**2*xc0 + muu3**2*sqrtS**2*xc0 + muu4**2*sqrtS**2*xc0 - 2*xc00 - 
     &	       sqrtS**2*xc1 + muu3**2*sqrtS**2*xc1 + 2*muu4**2*sqrtS**2*xc1 - sqrtS**2*xc2 + 2*muu3**2*sqrtS**2*xc2 + 
     &	       muu4**2*sqrtS**2*xc2) + g0L*(-(g1R*g2R*Ml1*Ml2*xc0) + g1L*g2L*muu3*muu4*sqrtS**2*(xc0 + xc1 + xc2)))/
     &	  (8.*Pi**2)
      	
		BpL2 = -(g0R*g1R*g2R*muu4*sqrtS*(xc0 + xc1 + xc12 + xc2) - g1L*g2R*(g0L*Ml2*xc2 + g0R*Ml1*(xc0 + xc2)) + 
     &	     g0L*g1L*g2L*muu3*sqrtS*(xc2 + xc22))/(4.*Pi**2)
      	
		BpR2 = -(g0L*g1L*g2L*muu4*sqrtS*(xc0 + xc1 + xc12 + xc2) - g1R*g2L*(g0R*Ml2*xc2 + g0L*Ml1*(xc0 + xc2)) + 
     &	     g0R*g1R*g2R*muu3*sqrtS*(xc2 + xc22))/(4.*Pi**2)
      	
		BpL3 = (-(g1L*g2R*(g0R*Ml1*xc1 + g0L*Ml2*(xc0 + xc1))) + g0R*g1R*g2R*muu4*sqrtS*(xc1 + xc11) + 
     &	    g0L*g1L*g2L*muu3*sqrtS*(xc0 + xc1 + xc12 + xc2))/(4.*Pi**2)
      	
		BpR3 = (-(g1R*g2L*(g0L*Ml1*xc1 + g0R*Ml2*(xc0 + xc1))) + g0L*g1L*g2L*muu4*sqrtS*(xc1 + xc11) + 
     &	    g0R*g1R*g2R*muu3*sqrtS*(xc0 + xc1 + xc12 + xc2))/(4.*Pi**2)
		
	endif

! Gluino loop *******
		
	do k=1,2
	 do l=1,2

! Correction to Z/W diagram

	    g0  = -ZSfSf(k,l,tt,gen)
	    g1L = -sqrt2*GS*USf(l,1,att,gen)
	    g1R = sqrt2*GS*USf(l,2,att,gen)
	    g2L = sqrt2*GS*USf(k,2,tt,gen)
	    g2R = -sqrt2*GS*USf(k,1,tt,gen)

	    Ml0 = MGl
	    Ml1 = MSf(l,att,gen)
	    Ml2 = MSf(k,tt,gen)
          
	    call myCget(muu4**2*sqrtS**2,sqrtS**2,muu3**2*sqrtS**2,MGl**2,MSf2(l,att,gen),MSf2(k,tt,gen))
	    call myBget(sqrtS**2,MSf2(l,att,gen),MSf2(k,tt,gen))
	
	
	    BL1 = BL1 - (g0*g1L*g2R*xc00)/(8.*Pi**2)
	    BR1 = BR1 - (g0*g1R*g2L*xc00)/(8.*Pi**2)
	    
	    BL2 = BL2 - g0*(g1R*g2L*muu4*sqrtS*(xc1 + 2*xc12) + 
     &           g1L*(-(g2L*Ml0*(xc0 + 2*xc2)) + 
     &              g2R*muu3*sqrtS*(xc2 + 2*xc22)))/(16.*Pi**2)
          
	    BR2 = BR2 - (g0*(g1L*g2R*muu4*sqrtS*(xc1 + 2*xc12) + 
     &           g1R*(-(g2R*Ml0*(xc0 + 2*xc2)) + 
     &              g2L*muu3*sqrtS*(xc2 + 2*xc22))))/(16.*Pi**2)
          
	    BL3 = BL3 + (g0*(g1R*g2L*muu4*sqrtS*(xc1 + 2*xc11) + 
     &          g1L*(-(g2L*Ml0*(xc0 + 2*xc1)) + 
     &             g2R*muu3*sqrtS*(2*xc12 + xc2))))/(16.*Pi**2)
          
	    BR3 = BR3 + (g0*(g1L*g2R*muu4*sqrtS*(xc1 + 2*xc11) + 
     &          g1R*(-(g2R*Ml0*(xc0 + 2*xc1)) + 
     &             g2L*muu3*sqrtS*(2*xc12 + xc2))))/(16.*Pi**2)


! Correction to photon diagram
	 	if (icha.eq.2) then

	   	    g0  = -GSfSf(k,l,tt,gen)
		    g1L = -sqrt2*GS*USf(l,1,att,gen)
		    g1R = sqrt2*GS*USf(l,2,att,gen)
		    g2L = sqrt2*GS*USf(k,2,tt,gen)
		    g2R = -sqrt2*GS*USf(k,1,tt,gen)

	    		BpL1 = BpL1 - (g0*g1L*g2R*xc00)/(8.*Pi**2)
	    		BpR1 = BpR1 - (g0*g1R*g2L*xc00)/(8.*Pi**2)

     	            BpL2 = BpL2 - g0*(g1R*g2L*muu4*sqrtS*(xc1 + 2*xc12) + 
     &                   g1L*(-(g2L*Ml0*(xc0 + 2*xc2)) + 
     &                      g2R*muu3*sqrtS*(xc2 + 2*xc22)))/(16.*Pi**2)
                  
     	            BpR2 = BpR2 - (g0*(g1L*g2R*muu4*sqrtS*(xc1 + 2*xc12) + 
     &                   g1R*(-(g2R*Ml0*(xc0 + 2*xc2)) + 
     &                      g2L*muu3*sqrtS*(xc2 + 2*xc22))))/(16.*Pi**2)
                  
     	            BpL3 = BpL3 + (g0*(g1R*g2L*muu4*sqrtS*(xc1 + 2*xc11) + 
     &                  g1L*(-(g2L*Ml0*(xc0 + 2*xc1)) + 
     &                     g2R*muu3*sqrtS*(2*xc12 + xc2))))/(16.*Pi**2)
                  
     	            BpR3 = BpR3 + (g0*(g1L*g2R*muu4*sqrtS*(xc1 + 2*xc11) + 
     &                  g1R*(-(g2R*Ml0*(xc0 + 2*xc1)) + 
     &                     g2L*muu3*sqrtS*(2*xc12 + xc2))))/(16.*Pi**2)

	 	endif
	
	 enddo
	enddo
	
cc	write(*,*)' --- TREE --- '
cc	write(*,*)'BL =',BL
cc	write(*,*)'BR =',BR
cc	write(*,*)' --- LOOP --- '
cc	write(*,*)'BL1 =',BL1
cc	write(*,*)'BR1 =',BR1
cc	write(*,*)'BL2 =',BL2
cc	write(*,*)'BR2 =',BR2
cc	write(*,*)'BL3 =',BL3
cc	write(*,*)'BR3 =',BR3


! Photon contributions only in cha-cha case
! cha cha -> b B : no t-channel
! cha cha -> t T : no u-channel
! Use couplings E/F, G/H to switch on/off t-/u-channels	

! for photon & interference only diagonal case considered (muu1=muu2, muu3=muu4)
	if (icha.eq.2) then
		
! Photon-Loop interference with photon tree-level
		ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + CF*NC*ApL**2*BpL*(8*(BpL1 + BpR1)*(1 + 2*(muu1**2 + muu3**2 - ured)**2 + 2*ured) +
     &			16*(BpL2 - BpL3 + BpR2 - BpR3)*muu3*sqrtS*(((-1 + muu1)*muu1 + muu3**2 - ured)*(muu1 + muu1**2 + muu3**2 - ured) + ured))	
	
cc	write(*,*)'ChiChi P Vertex mit P = ', ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! Photon-Loop interference with Z-boson tree-level

		ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + ApL*CF*NC*sqrtS**2*(8*((BpL1*(BR*(1 + 2*muu1**2)*muu3**2 + 
     &            	BL*(muu1**2 + muu1**4 + muu3**4 - 2*(muu1**2 + muu3**2)*ured + ured**2)) + 
     &         		BpR1*(BL*(1 + 2*muu1**2)*muu3**2 + BR*(muu1**4 + (1 - muu3**2 + ured)**2 - muu1**2*(1 + 2*ured))))*dconjg(AL) + 
     &      		(BR*(BpL1*(1 + 2*muu1**2)*muu3**2 + BpR1*(muu1**2 + muu1**4 + muu3**4 - 2*(muu1**2 + muu3**2)*ured + ured**2)) + 
     &         		BL*(BpR1*(1 + 2*muu1**2)*muu3**2 + BpL1*(muu1**4 + (1 - muu3**2 + ured)**2 - muu1**2*(1 + 2*ured))))*dconjg(AR)) +
     &    		4*(BpL2 - BpL3 + BpR2 - BpR3)*(BL + BR)*(dconjg(AL) + dconjg(AR))*muu3*sqrtS*
     &         		(((-1 + muu1)*muu1 + muu3**2 - ured)*(muu1 + muu1**2 + muu3**2 - ured) + ured))/(dreal(dconjg(DenSmV)))
	
cc	write(*,*)'ChiChi P Vertex mit Z = ', ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! Z-Loop interference with photon tree-level

		ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + ApL*BpL*CF*NC*sqrtS**2*(8*(AR*(BR1*(muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*ured) + (muu3**2 - ured)**2) + 
     &         		BL1*(muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*tred) + (muu3**2 - tred)**2)) + 
     &      		AL*(BL1*(muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*ured) + (muu3**2 - ured)**2) + 
     &         		BR1*(muu1**4 + muu3**2 + muu1**2*(1 + 2*muu3**2 - 2*tred) + (muu3**2 - tred)**2))) +
     &			2*(AL + AR)*(BL2 - BL3 + BR2 - BR3)*muu3*sqrtS*(-1 + 4*muu3**2 + (ured - tred)**2))/(dreal(dconjg(DenSmV)))

cc	write(*,*)'ChiChi Z Vertex mit P = ', ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! Photon-Loop interference with Higgs-boson tree-level

	        do i= imin,imax

	 	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + CF*NC*ApL*sqrtS**2*(4*(BpL1 + BpR1)*muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured)*(dconjg(CL(i)) + dconjg(CR(i)))*
     &        		(dconjg(DL(i)) + dconjg(DR(i))) +	
     &   		2*muu1*sqrtS*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured)*(dconjg(CL(i)) + dconjg(CR(i)))*
     &       		(-(BpL2*dconjg(DL(i))) + BpL3*dconjg(DL(i)) - BpR2*dconjg(DR(i)) + BpR3*dconjg(DR(i)) + 
     &        		2*(BpL2 - BpL3 + BpR2 - BpR3)*muu3**2*(dconjg(DL(i)) + dconjg(DR(i)))))/(dreal(dconjg(DenSmH(i))))
	
	        enddo
	
cc	write(*,*)'ChiChi P Vertex mit H = ', ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! Photon-Loop interference with T-channel tree-level (old U-channel)
			
	do j =1,2
			
	 	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + CF*NC*ApL*(4*(BpR1*(muu3*dconjg(FL(indj,j,tt,gen))*(muu1*(1 + muu1**2 + muu3**2 - tred)*dconjg(ELc(indi,j,tt,gen)) + 
     &            	(1 + 2*muu1**2)*muu3*dconjg(ER(indi,j,tt,gen))) + 
     &         		dconjg(FR(indj,j,tt,gen))*((muu1**2 + muu1**4 + muu3**4 - 2*(muu1**2 + muu3**2)*tred + tred**2)*
     &             	dconjg(ELc(indi,j,tt,gen)) + muu1*muu3*(1 + muu1**2 + muu3**2 - tred)*dconjg(ER(indi,j,tt,gen)))) + 
     &      		BpL1*(muu3*dconjg(FR(indj,j,tt,gen))*((1 + 2*muu1**2)*muu3*dconjg(ELc(indi,j,tt,gen)) + 
     &            	muu1*(1 + muu1**2 + muu3**2 - tred)*dconjg(ER(indi,j,tt,gen))) + 
     &         		dconjg(FL(indj,j,tt,gen))*(muu1*muu3*(1 + muu1**2 + muu3**2 - tred)*dconjg(ELc(indi,j,tt,gen)) + 
     &            	(muu1**2 + muu1**4 + muu3**4 - 2*(muu1**2 + muu3**2)*tred + tred**2)*dconjg(ER(indi,j,tt,gen)))))
     &      		-2*sqrtS*(BpL2*(muu3*dconjg(FR(indj,j,tt,gen))*
     &          	((((-1 + muu1)*muu1 + muu3**2 - tred)*(muu1 + muu1**2 + muu3**2 - tred) + tred)*dconjg(ELc(indi,j,tt,gen)) + 
     &            	muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred)*dconjg(ER(indi,j,tt,gen))) + 
     &         		dconjg(FL(indj,j,tt,gen))*(muu1*(2*muu3**4 + muu1**2*(-1 + 2*muu3**2) + tred - 2*muu3**2*tred)*
     &             	dconjg(ELc(indi,j,tt,gen)) + muu3*(((-1 + muu1)*muu1 + muu3**2 - tred)*(muu1 + muu1**2 + muu3**2 - tred) + tred)*
     &             	dconjg(ER(indi,j,tt,gen)))) - BpL3*(muu3*dconjg(FR(indj,j,tt,gen))*
     &          	((((-1 + muu1)*muu1 + muu3**2 - tred)*(muu1 + muu1**2 + muu3**2 - tred) + tred)*dconjg(ELc(indi,j,tt,gen)) + 
     &            	muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred)*dconjg(ER(indi,j,tt,gen))) + 
     &         		dconjg(FL(indj,j,tt,gen))*(muu1*(2*muu3**4 + muu1**2*(-1 + 2*muu3**2) + tred - 2*muu3**2*tred)*
     &             	dconjg(ELc(indi,j,tt,gen)) + muu3*(((-1 + muu1)*muu1 + muu3**2 - tred)*(muu1 + muu1**2 + muu3**2 - tred) + tred)*
     &             	dconjg(ER(indi,j,tt,gen)))) + (BpR2 - BpR3)*
     &       		(muu3*dconjg(FL(indj,j,tt,gen))*(muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*tred)*dconjg(ELc(indi,j,tt,gen)) + 
     &            	(((-1 + muu1)*muu1 + muu3**2 - tred)*(muu1 + muu1**2 + muu3**2 - tred) + tred)*dconjg(ER(indi,j,tt,gen))) + 
     &         		dconjg(FR(indj,j,tt,gen))*(muu3*(((-1 + muu1)*muu1 + muu3**2 - tred)*(muu1 + muu1**2 + muu3**2 - tred) + tred)*
     &             	dconjg(ELc(indi,j,tt,gen)) + muu1*(2*muu3**4 + muu1**2*(-1 + 2*muu3**2) + tred - 2*muu3**2*tred)*
     &             	dconjg(ER(indi,j,tt,gen)))))
     &       		)/(tred - MSf(j,itypeT,gen)**2/sqrtS**2)
	 
	enddo
		
cc	write(*,*)'ChiChi P Vertex mit t = ', ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
	do j =1,2
	
! Photon-Loop interference with U-channel tree-level (old T-channel)

		ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + NC*CF*ApL*(4*(BpL1*(muu3*dconjg(HL(indj,j,tt,gen))*(muu1*(1 + muu1**2 + muu3**2 - ured)*dconjg(GL(indi,j,tt,gen)) +     
     &            (1 + 2*muu1**2)*muu3*dconjg(GR(indi,j,tt,gen))) +     
     &         dconjg(HR(indj,j,tt,gen))*((muu1**2 + muu1**4 + muu3**4 - 2*(muu1**2 + muu3**2)*ured + ured**2)*dconjg(GL(indi,j,tt,gen)) +     
     &            muu1*muu3*(1 + muu1**2 + muu3**2 - ured)*dconjg(GR(indi,j,tt,gen)))) +     
     &      BpR1*(muu3*dconjg(HR(indj,j,tt,gen))*((1 + 2*muu1**2)*muu3*dconjg(GL(indi,j,tt,gen)) +     
     &            muu1*(1 + muu1**2 + muu3**2 - ured)*dconjg(GR(indi,j,tt,gen))) +     
     &         dconjg(HL(indj,j,tt,gen))*(muu1*muu3*(1 + muu1**2 + muu3**2 - ured)*dconjg(GL(indi,j,tt,gen)) +     
     &            (muu1**2 + muu1**4 + muu3**4 - 2*(muu1**2 + muu3**2)*ured + ured**2)*dconjg(GR(indi,j,tt,gen)))))
     &     + 2*sqrtS*(BpL2*(muu3*dconjg(HR(indj,j,tt,gen))*((((-1 + muu1)*muu1 + muu3**2 - ured)*(muu1 + muu1**2 + muu3**2 - ured) + ured)*
     &             dconjg(GL(indi,j,tt,gen)) + muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured)*dconjg(GR(indi,j,tt,gen))) +     
     &         dconjg(HL(indj,j,tt,gen))*(muu1*(2*muu3**4 + muu1**2*(-1 + 2*muu3**2) + ured - 2*muu3**2*ured)*dconjg(GL(indi,j,tt,gen)) +     
     &            muu3*(((-1 + muu1)*muu1 + muu3**2 - ured)*(muu1 + muu1**2 + muu3**2 - ured) + ured)*dconjg(GR(indi,j,tt,gen)))) -     
     &      BpL3*(muu3*dconjg(HR(indj,j,tt,gen))*((((-1 + muu1)*muu1 + muu3**2 - ured)*(muu1 + muu1**2 + muu3**2 - ured) + ured)*
     &             dconjg(GL(indi,j,tt,gen)) + muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured)*dconjg(GR(indi,j,tt,gen))) +     
     &         dconjg(HL(indj,j,tt,gen))*(muu1*(2*muu3**4 + muu1**2*(-1 + 2*muu3**2) + ured - 2*muu3**2*ured)*dconjg(GL(indi,j,tt,gen)) +     
     &            muu3*(((-1 + muu1)*muu1 + muu3**2 - ured)*(muu1 + muu1**2 + muu3**2 - ured) + ured)*dconjg(GR(indi,j,tt,gen)))) +     
     &      (BpR2 - BpR3)*(muu3*dconjg(HL(indj,j,tt,gen))*
     &          (muu1*muu3*(-1 + 2*muu1**2 + 2*muu3**2 - 2*ured)*dconjg(GL(indi,j,tt,gen)) +     
     &            (((-1 + muu1)*muu1 + muu3**2 - ured)*(muu1 + muu1**2 + muu3**2 - ured) + ured)*dconjg(GR(indi,j,tt,gen))) +     
     &         dconjg(HR(indj,j,tt,gen))*(muu3*(((-1 + muu1)*muu1 + muu3**2 - ured)*(muu1 + muu1**2 + muu3**2 - ured) + ured)*
     &             dconjg(GL(indi,j,tt,gen)) + muu1*(2*muu3**4 + muu1**2*(-1 + 2*muu3**2) + ured - 2*muu3**2*ured)*dconjg(GR(indi,j,tt,gen))))))/
     &  	(-ured + MSf(j,itypeU,gen)**2/sqrtS**2)
			
	enddo
		
cc	write(*,*)'ChiChi P Vertex mit u = ', ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
	endif



! Z/W-loop interference with Z/W-boson tree-level

        ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + CF*NC*sqrtS**4*(-8*((AR*muu1*muu2*((BL*BL1 + BR*BR1)*(-1 + muu3**2) - 4*(BL1*BR + BL*BR1)*muu3*muu4 + (BL*BL1 + BR*BR1)*muu4**2) + 
     &    	AL*(BL1*(BR*(-1 + muu1**2 + muu2**2)*muu3*muu4 - BL*(muu2**2 + muu3**2 - ured)*(muu1**2 + muu4**2 - ured)) + 
     &          BR1*(BL*(-1 + muu1**2 + muu2**2)*muu3*muu4 - BR*(muu1**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred))))*dconjg(AL)
     &       	+ (AL*muu1*muu2*((BL*BL1 + BR*BR1)*(-1 + muu3**2) - 4*(BL1*BR + BL*BR1)*muu3*muu4 + (BL*BL1 + BR*BR1)*muu4**2) + 
     &         	AR*(BR1*(BL*(-1 + muu1**2 + muu2**2)*muu3*muu4 - BR*(muu2**2 + muu3**2 - ured)*(muu1**2 + muu4**2 - ured)) + 
     &          BL1*(BR*(-1 + muu1**2 + muu2**2)*muu3*muu4 - BL*(muu1**2 + muu3**2 - tred)*(muu2**2 + muu4**2 - tred))))*dconjg(AR))
     &     	+ 2*sqrtS*(AL*((BR2*(BR*muu3*(-1 + muu2**2 + muu4**2 + muu3**2*(1 + 2*muu4**2) - (muu2**2 + muu3**2 + muu4**2 - ured)*ured + 
     &          muu1**2*(1 + 2*muu2**2 - ured - tred) - (muu2**2 + muu3**2 + muu4**2)*tred + tred**2) + 
     &          2*BL*muu4*(-muu3**4 + muu1**2*(-muu2**2 - 2*muu3**2 + ured) + (muu2**2 - ured)*tred + 
     &          muu3**2*(1 - 2*muu2**2 + ured + tred))) + 
     &          BL2*(BL*muu3*(-1 + muu2**2 + muu4**2 + muu3**2*(1 + 2*muu4**2) - (muu2**2 + muu3**2 + muu4**2 - ured)*ured + 
     &          muu1**2*(1 + 2*muu2**2 - ured - tred) - (muu2**2 + muu3**2 + muu4**2)*tred + tred**2) + 
     &          2*BR*muu4*(-muu3**4 + muu1**2*(-muu2**2 - 2*muu3**2 + ured) + (muu2**2 - ured)*tred + 
     &          muu3**2*(1 - 2*muu2**2 + ured + tred))))*dconjg(AL) + 
     &         	2*muu1*muu2*muu3*((BL*BL2 + BR*BR2)*(-1 + muu3**2) + 2*(BL2*BR + BL*BR2)*muu3*muu4 + (BL*BL2 + BR*BR2)*muu4**2)*
     &          dconjg(AR)) + AR*(2*muu1*muu2*muu3*((BL*BL2 + BR*BR2)*(-1 + muu3**2) + 2*(BL2*BR + BL*BR2)*muu3*muu4 + 
     &          (BL*BL2 + BR*BR2)*muu4**2)*dconjg(AL) + 
     &         	(BR2*(BR*muu3*(-1 + muu2**2 + muu4**2 + muu3**2*(1 + 2*muu4**2) - (muu2**2 + muu3**2 + muu4**2 - ured)*ured + 
     &                  muu1**2*(1 + 2*muu2**2 - ured - tred) - (muu2**2 + muu3**2 + muu4**2)*tred + tred**2) + 
     &               2*BL*muu4*(-muu3**4 + muu1**2*(-muu2**2 - 2*muu3**2 + ured) + (muu2**2 - ured)*tred + 
     &                  muu3**2*(1 - 2*muu2**2 + ured + tred))) + 
     &            BL2*(BL*muu3*(-1 + muu2**2 + muu4**2 + muu3**2*(1 + 2*muu4**2) - (muu2**2 + muu3**2 + muu4**2 - ured)*ured + 
     &                  muu1**2*(1 + 2*muu2**2 - ured - tred) - (muu2**2 + muu3**2 + muu4**2)*tred + tred**2) + 
     &               2*BR*muu4*(-muu3**4 + muu1**2*(-muu2**2 - 2*muu3**2 + ured) + (muu2**2 - ured)*tred + 
     &                  muu3**2*(1 - 2*muu2**2 + ured + tred))))*dconjg(AR)))
     &    	+ 2*sqrtS*(AL*((BL3*(-(BR*muu4*(-1 + muu2**2 + muu4**2 + muu3**2*(1 + 2*muu4**2) - (muu2**2 + muu3**2 + muu4**2 - ured)*ured + 
     &                    muu1**2*(1 + 2*muu2**2 - ured - tred) - (muu2**2 + muu3**2 + muu4**2)*tred + tred**2)) + 
     &               2*BL*muu3*(muu4**4 + muu1**2*(muu2**2 + 2*muu4**2 - tred) + muu4**2*(-1 + 2*muu2**2 - ured - tred) + 
     &                  ured*(-muu2**2 + tred))) + BR3*
     &             (-(BL*muu4*(-1 + muu2**2 + muu4**2 + muu3**2*(1 + 2*muu4**2) - (muu2**2 + muu3**2 + muu4**2 - ured)*ured + 
     &                    muu1**2*(1 + 2*muu2**2 - ured - tred) - (muu2**2 + muu3**2 + muu4**2)*tred + tred**2)) + 
     &               2*BR*muu3*(muu4**4 + muu1**2*(muu2**2 + 2*muu4**2 - tred) + muu4**2*(-1 + 2*muu2**2 - ured - tred) + 
     &                  ured*(-muu2**2 + tred))))*dconjg(AL) - 
     &         	2*muu1*muu2*muu4*((BL3*BR + BL*BR3)*(-1 + muu3**2) + 2*(BL*BL3 + BR*BR3)*muu3*muu4 + (BL3*BR + BL*BR3)*muu4**2)*
     &          dconjg(AR)) + AR*(-2*muu1*muu2*muu4*((BL3*BR + BL*BR3)*(-1 + muu3**2) + 2*(BL*BL3 + BR*BR3)*muu3*muu4 + 
     &            (BL3*BR + BL*BR3)*muu4**2)*dconjg(AL) + 
     &         	(BL3*(-(BR*muu4*(-1 + muu2**2 + muu4**2 + muu3**2*(1 + 2*muu4**2) - (muu2**2 + muu3**2 + muu4**2 - ured)*ured + 
     &                    muu1**2*(1 + 2*muu2**2 - ured - tred) - (muu2**2 + muu3**2 + muu4**2)*tred + tred**2)) + 
     &               2*BL*muu3*(muu4**4 + muu1**2*(muu2**2 + 2*muu4**2 - tred) + muu4**2*(-1 + 2*muu2**2 - ured - tred) + 
     &                  ured*(-muu2**2 + tred))) + BR3*
     &             (-(BL*muu4*(-1 + muu2**2 + muu4**2 + muu3**2*(1 + 2*muu4**2) - (muu2**2 + muu3**2 + muu4**2 - ured)*ured + 
     &                    muu1**2*(1 + 2*muu2**2 - ured - tred) - (muu2**2 + muu3**2 + muu4**2)*tred + tred**2)) + 
     &               2*BR*muu3*(muu4**4 + muu1**2*(muu2**2 + 2*muu4**2 - tred) + muu4**2*(-1 + 2*muu2**2 - ured - tred) + 
     &                  ured*(-muu2**2 + tred))))*dconjg(AR)))
     &       	)/(DenSmV*dconjg(DenSmV))

cc	write(*,*)'ChiChi Z Vertex mit Z = ', ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex

	
! Z/W-loop interference with Higgs-boson tree-level

	do i= imin,imax

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + CF*NC*sqrtS**4*(-4*(-(muu2*(AL*dconjg(CL(i)) + AR*dconjg(CR(i)))*
     &         	(-(muu4*(muu1**2 + muu3**2 - tred)*(BR1*dconjg(DL(i)) + BL1*dconjg(DR(i)))) + 
     &           muu3*(muu1**2 + muu4**2 - ured)*(BL1*dconjg(DL(i)) + BR1*dconjg(DR(i))))) + 
     &      	muu1*(AR*dconjg(CL(i)) + AL*dconjg(CR(i)))*
     &       	(-(muu4*(muu2**2 + muu3**2 - ured)*(BR1*dconjg(DL(i)) + BL1*dconjg(DR(i)))) + 
     &         	muu3*(muu2**2 + muu4**2 - tred)*(BL1*dconjg(DL(i)) + BR1*dconjg(DR(i)))))
     &  	+ 2*sqrtS*(muu1*(muu2**2 + muu3**2 - ured)*(AR*dconjg(CL(i)) + AL*dconjg(CR(i))) - 
     &      	muu2*(muu1**2 + muu3**2 - tred)*(AL*dconjg(CL(i)) + AR*dconjg(CR(i))))*
     &    	(2*muu3*muu4*(BR2*dconjg(DL(i)) + BL2*dconjg(DR(i))) + (-1 + muu3**2)*(BL2*dconjg(DL(i)) + BR2*dconjg(DR(i))) + 
     &      	muu4**2*(BL2*dconjg(DL(i)) + BR2*dconjg(DR(i))))
     &  	+ 2*sqrtS*(muu1*(muu2**2 + muu4**2 - tred)*(AR*dconjg(CL(i)) + AL*dconjg(CR(i))) - 
     &      	muu2*(muu1**2 + muu4**2 - ured)*(AL*dconjg(CL(i)) + AR*dconjg(CR(i))))*
     &    	(2*muu3*muu4*(BR3*dconjg(DL(i)) + BL3*dconjg(DR(i))) + (-1 + muu3**2)*(BL3*dconjg(DL(i)) + BR3*dconjg(DR(i))) + 
     &      	muu4**2*(BL3*dconjg(DL(i)) + BR3*dconjg(DR(i)))))/(dreal(DenSmV*dconjg(DenSmH(i))))

	enddo
	
cc	write(*,*)'ChiChi Z Vertex mit H = ', ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex

	do j=1,2

! Z/W-loop interference with U-channel tree-level (old T-channel)

 	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + CF*NC*sqrtS**2*(
     &     	-4*(AR*(BR1*((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,j,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,j,tt,gen)))*
     &          (2*muu1*muu4*dconjg(GL(indi,j,tt,gen)) + (muu1**2 + muu4**2 - ured)*dconjg(GR(indi,j,tt,gen))) + 
     &         	BL1*(muu3*dconjg(HL(indj,j,tt,gen))*(muu1*(muu2**2 + muu4**2 - tred)*dconjg(GL(indi,j,tt,gen)) - 
     &               (-1 + muu1**2 + muu2**2)*muu4*dconjg(GR(indi,j,tt,gen))) + 
     &            muu2*dconjg(HR(indj,j,tt,gen))*(-(muu1*(-1 + muu3**2 + muu4**2)*dconjg(GL(indi,j,tt,gen))) + 
     &               muu4*(muu1**2 + muu3**2 - tred)*dconjg(GR(indi,j,tt,gen))))) + 
     &      	AL*(BL1*(2*muu2*muu3*dconjg(HL(indj,j,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,j,tt,gen)))*
     &          ((muu1**2 + muu4**2 - ured)*dconjg(GL(indi,j,tt,gen)) + 2*muu1*muu4*dconjg(GR(indi,j,tt,gen))) + 
     &         	BR1*(muu2*dconjg(HL(indj,j,tt,gen))*(muu4*(muu1**2 + muu3**2 - tred)*dconjg(GL(indi,j,tt,gen)) - 
     &               muu1*(-1 + muu3**2 + muu4**2)*dconjg(GR(indi,j,tt,gen))) + 
     &            muu3*dconjg(HR(indj,j,tt,gen))*(-((-1 + muu1**2 + muu2**2)*muu4*dconjg(GL(indi,j,tt,gen))) + 
     &               muu1*(muu2**2 + muu4**2 - tred)*dconjg(GR(indi,j,tt,gen))))))
     &    	+ sqrtS*(AL*(-((BL2*muu3 - BR3*muu4)*(2*muu2*muu3*dconjg(HL(indj,j,tt,gen))*
     &               ((muu1**2 + muu4**2 - ured)*dconjg(GL(indi,j,tt,gen)) + 2*muu1*muu4*dconjg(GR(indi,j,tt,gen))) + 
     &              dconjg(HR(indj,j,tt,gen))*((-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - 
     &                    (muu2**2 + muu3**2 + muu4**2 - ured)*ured + muu1**2*(1 + 2*muu2**2 - ured - tred) - 
     &                    (muu2**2 + muu3**2 + muu4**2)*tred + tred**2)*dconjg(GL(indi,j,tt,gen)) + 
     &                 2*muu1*muu4*(muu2**2 + muu3**2 - ured)*dconjg(GR(indi,j,tt,gen))))) - 
     &         	2*BR2*(muu2*muu3*dconjg(HL(indj,j,tt,gen))*
     &             (-(muu4*(muu1**2 + muu3**2 - tred)*dconjg(GL(indi,j,tt,gen))) + 
     &               muu1*(-1 + muu3**2 + muu4**2)*dconjg(GR(indi,j,tt,gen))) + 
     &            dconjg(HR(indj,j,tt,gen))*(muu4*(-muu3**4 + muu1**2*(-muu2**2 - 2*muu3**2 + ured) + (muu2**2 - ured)*tred + 
     &                  muu3**2*(1 - 2*muu2**2 + ured + tred))*dconjg(GL(indi,j,tt,gen)) + 
     &               muu1*(muu3**4 + muu2**2*(-1 + 2*muu3**2 + muu4**2) + ured - muu4**2*ured + muu3**2*(-1 + 2*muu4**2 - ured - tred))*
     &                dconjg(GR(indi,j,tt,gen)))) + 2*BL3*
     &          (muu2*dconjg(HL(indj,j,tt,gen))*((muu4**4 + muu1**2*(-1 + muu3**2 + 2*muu4**2) + ured - muu3**2*ured + 
     &                  muu4**2*(-1 + 2*muu3**2 - ured - tred))*dconjg(GL(indi,j,tt,gen)) + 
     &               muu1*muu4*(-1 + muu3**2 + muu4**2)*dconjg(GR(indi,j,tt,gen))) + 
     &            muu3*dconjg(HR(indj,j,tt,gen))*((-muu4**4 + ured*(muu2**2 - tred) + muu1**2*(-muu2**2 - 2*muu4**2 + tred) + 
     &                  muu4**2*(1 - 2*muu2**2 + ured + tred))*dconjg(GL(indi,j,tt,gen)) - 
     &               muu1*muu4*(muu2**2 + muu4**2 - tred)*dconjg(GR(indi,j,tt,gen))))) + 
     &      	AR*(-((BR2*muu3 - BL3*muu4)*(2*muu2*muu3*dconjg(HR(indj,j,tt,gen))*
     &               (2*muu1*muu4*dconjg(GL(indi,j,tt,gen)) + (muu1**2 + muu4**2 - ured)*dconjg(GR(indi,j,tt,gen))) + 
     &              dconjg(HL(indj,j,tt,gen))*(2*muu1*muu4*(muu2**2 + muu3**2 - ured)*dconjg(GL(indi,j,tt,gen)) + 
     &                 (-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - (muu2**2 + muu3**2 + muu4**2 - ured)*ured + 
     &                    muu1**2*(1 + 2*muu2**2 - ured - tred) - (muu2**2 + muu3**2 + muu4**2)*tred + tred**2)*dconjg(GR(indi,j,tt,gen))))) - 
     &              2*BL2*(muu2*muu3*dconjg(HR(indj,j,tt,gen))*
     &             (muu1*(-1 + muu3**2 + muu4**2)*dconjg(GL(indi,j,tt,gen)) - muu4*(muu1**2 + muu3**2 - tred)*dconjg(GR(indi,j,tt,gen)))
     &             + dconjg(HL(indj,j,tt,gen))*(muu1*(muu3**4 + muu2**2*(-1 + 2*muu3**2 + muu4**2) + ured - muu4**2*ured + 
     &                  muu3**2*(-1 + 2*muu4**2 - ured - tred))*dconjg(GL(indi,j,tt,gen)) + 
     &               muu4*(-muu3**4 + muu1**2*(-muu2**2 - 2*muu3**2 + ured) + (muu2**2 - ured)*tred + 
     &                  muu3**2*(1 - 2*muu2**2 + ured + tred))*dconjg(GR(indi,j,tt,gen)))) + 
     &         	2*BR3*(muu2*dconjg(HR(indj,j,tt,gen))*(muu1*muu4*(-1 + muu3**2 + muu4**2)*dconjg(GL(indi,j,tt,gen)) + 
     &               (muu4**4 + muu1**2*(-1 + muu3**2 + 2*muu4**2) + ured - muu3**2*ured + muu4**2*(-1 + 2*muu3**2 - ured - tred))*
     &                dconjg(GR(indi,j,tt,gen))) + muu3*dconjg(HL(indj,j,tt,gen))*
     &             (-(muu1*muu4*(muu2**2 + muu4**2 - tred)*dconjg(GL(indi,j,tt,gen))) + 
     &               (-muu4**4 + ured*(muu2**2 - tred) + muu1**2*(-muu2**2 - 2*muu4**2 + tred) + muu4**2*(1 - 2*muu2**2 + ured + tred))*
     &                dconjg(GR(indi,j,tt,gen)))))) 
     &          )/(DenSmV1*(ured - MSf(j,itypeU,gen)**2/sqrtS**2))
	
	enddo
	
cc	write(*,*)'ChiChi Z Vertex mit u = ', ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex

	do j=1,2
	
! Z-loop interference with T-channel tree-level (old U-channel)

 	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + CF*NC*sqrtS**2*(
     &		4*(AL*(BR1*(2*muu2*muu4*dconjg(FL(indj,j,tt,gen)) + (muu2**2 + muu4**2 - tred)*dconjg(FR(indj,j,tt,gen)))*
     &          ((muu1**2 + muu3**2 - tred)*dconjg(ELc(indi,j,tt,gen)) + 2*muu1*muu3*dconjg(ER(indi,j,tt,gen))) + 
     &         	BL1*(muu2*dconjg(FL(indj,j,tt,gen))*(muu3*(muu1**2 + muu4**2 - ured)*dconjg(ELc(indi,j,tt,gen)) - 
     &               muu1*(-1 + muu3**2 + muu4**2)*dconjg(ER(indi,j,tt,gen))) + 
     &            muu4*dconjg(FR(indj,j,tt,gen))*(-((-1 + muu1**2 + muu2**2)*muu3*dconjg(ELc(indi,j,tt,gen))) + 
     &               muu1*(muu2**2 + muu3**2 - ured)*dconjg(ER(indi,j,tt,gen))))) + 
     &      	AR*(BL1*((muu2**2 + muu4**2 - tred)*dconjg(FL(indj,j,tt,gen)) + 2*muu2*muu4*dconjg(FR(indj,j,tt,gen)))*
     &          (2*muu1*muu3*dconjg(ELc(indi,j,tt,gen)) + (muu1**2 + muu3**2 - tred)*dconjg(ER(indi,j,tt,gen))) + 
     &         	BR1*(muu4*dconjg(FL(indj,j,tt,gen))*(muu1*(muu2**2 + muu3**2 - ured)*dconjg(ELc(indi,j,tt,gen)) - 
     &               (-1 + muu1**2 + muu2**2)*muu3*dconjg(ER(indi,j,tt,gen))) + 
     &            muu2*dconjg(FR(indj,j,tt,gen))*(-(muu1*(-1 + muu3**2 + muu4**2)*dconjg(ELc(indi,j,tt,gen))) + 
     &               muu3*(muu1**2 + muu4**2 - ured)*dconjg(ER(indi,j,tt,gen))))))
     &     + 	sqrtS*(AL*(2*BL2*(muu2*dconjg(FL(indj,j,tt,gen))*
     &             ((muu1**2*(-1 + 2*muu3**2 + muu4**2) + muu3**2*(-1 + muu3**2 + 2*muu4**2 - ured) + tred - (muu3**2 + muu4**2)*tred)*
     &                dconjg(ELc(indi,j,tt,gen)) + muu1*muu3*(-1 + muu3**2 + muu4**2)*dconjg(ER(indi,j,tt,gen))) + 
     &            muu4*dconjg(FR(indj,j,tt,gen))*((-muu3**4 + muu1**2*(-muu2**2 - 2*muu3**2 + ured) + (muu2**2 - ured)*tred + 
     &                  muu3**2*(1 - 2*muu2**2 + ured + tred))*dconjg(ELc(indi,j,tt,gen)) - 
     &               muu1*muu3*(muu2**2 + muu3**2 - ured)*dconjg(ER(indi,j,tt,gen)))) + 
     &         	(BR2*muu3 - BL3*muu4)*(2*muu2*muu4*dconjg(FL(indj,j,tt,gen))*
     &             ((muu1**2 + muu3**2 - tred)*dconjg(ELc(indi,j,tt,gen)) + 2*muu1*muu3*dconjg(ER(indi,j,tt,gen))) + 
     &            dconjg(FR(indj,j,tt,gen))*((-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - 
     &                  (muu2**2 + muu3**2 + muu4**2 - ured)*ured + muu1**2*(1 + 2*muu2**2 - ured - tred) - 
     &                  (muu2**2 + muu3**2 + muu4**2)*tred + tred**2)*dconjg(ELc(indi,j,tt,gen)) + 
     &               2*muu1*muu3*(muu2**2 + muu4**2 - tred)*dconjg(ER(indi,j,tt,gen)))) - 
     &         	2*BR3*(muu2*muu4*dconjg(FL(indj,j,tt,gen))*
     &             (-(muu3*(muu1**2 + muu4**2 - ured)*dconjg(ELc(indi,j,tt,gen))) + 
     &               muu1*(-1 + muu3**2 + muu4**2)*dconjg(ER(indi,j,tt,gen))) + 
     &            dconjg(FR(indj,j,tt,gen))*(muu3*(-muu4**4 + ured*(muu2**2 - tred) + muu1**2*(-muu2**2 - 2*muu4**2 + tred) + 
     &                  muu4**2*(1 - 2*muu2**2 + ured + tred))*dconjg(ELc(indi,j,tt,gen)) + 
     &               muu1*(muu4**4 + muu2**2*(-1 + muu3**2 + 2*muu4**2) + muu4**2*(-1 + 2*muu3**2 - ured - tred) + tred - muu3**2*tred)*
     &                dconjg(ER(indi,j,tt,gen))))) + AR*
     &       	((BL2*muu3 - BR3*muu4)*(2*muu2*muu4*dconjg(FR(indj,j,tt,gen))*
     &             (2*muu1*muu3*dconjg(ELc(indi,j,tt,gen)) + (muu1**2 + muu3**2 - tred)*dconjg(ER(indi,j,tt,gen))) + 
     &            dconjg(FL(indj,j,tt,gen))*(2*muu1*muu3*(muu2**2 + muu4**2 - tred)*dconjg(ELc(indi,j,tt,gen)) + 
     &               (-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - (muu2**2 + muu3**2 + muu4**2 - ured)*ured + 
     &                  muu1**2*(1 + 2*muu2**2 - ured - tred) - (muu2**2 + muu3**2 + muu4**2)*tred + tred**2)*dconjg(ER(indi,j,tt,gen))))
     &           + 2*BR2*(muu2*dconjg(FR(indj,j,tt,gen))*
     &             (muu1*muu3*(-1 + muu3**2 + muu4**2)*dconjg(ELc(indi,j,tt,gen)) + 
     &               (muu1**2*(-1 + 2*muu3**2 + muu4**2) + muu3**2*(-1 + muu3**2 + 2*muu4**2 - ured) + tred - (muu3**2 + muu4**2)*tred)*
     &                dconjg(ER(indi,j,tt,gen))) + muu4*dconjg(FL(indj,j,tt,gen))*
     &             (-(muu1*muu3*(muu2**2 + muu3**2 - ured)*dconjg(ELc(indi,j,tt,gen))) + 
     &               (-muu3**4 + muu1**2*(-muu2**2 - 2*muu3**2 + ured) + (muu2**2 - ured)*tred + muu3**2*(1 - 2*muu2**2 + ured + tred))*
     &                dconjg(ER(indi,j,tt,gen)))) - 2*BL3*
     &          (muu2*muu4*dconjg(FR(indj,j,tt,gen))*(muu1*(-1 + muu3**2 + muu4**2)*dconjg(ELc(indi,j,tt,gen)) - 
     &               muu3*(muu1**2 + muu4**2 - ured)*dconjg(ER(indi,j,tt,gen))) + 
     &            dconjg(FL(indj,j,tt,gen))*(muu1*(muu4**4 + muu2**2*(-1 + muu3**2 + 2*muu4**2) + 
     &                  muu4**2*(-1 + 2*muu3**2 - ured - tred) + tred - muu3**2*tred)*dconjg(ELc(indi,j,tt,gen)) + 
     &               muu3*(-muu4**4 + ured*(muu2**2 - tred) + muu1**2*(-muu2**2 - 2*muu4**2 + tred) + 
     &                  muu4**2*(1 - 2*muu2**2 + ured + tred))*dconjg(ER(indi,j,tt,gen)))))) 
     &             )/(DenSmV1*(tred - MSf(j,itypeT,gen)**2/sqrtS**2))

	 enddo

cc	write(*,*)'ChiChi Z Vertex mit t = ', ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex

****************************************
*********** Loop with Higgs ************
****************************************

	do i= imin,imax

!*******************************************************************************************	
!************ Gluon loop *******
!*******************************

	call myCget(muu4**2*sqrtS**2,sqrtS**2,muu3**2*sqrtS**2,0d0,muu4**2*sqrtS**2,muu3**2*sqrtS**2)
	call myBget(sqrtS**2,muu4**2*sqrtS**2,muu3**2*sqrtS**2)
	
	g0L = DL(i)
	g0R = DR(i)
	g1L = -GS
	g1R = -GS
	g2L = -GS
	g2R = -GS
	
	Ml0 = 0d0
	Ml1 = muu4*sqrtS
	Ml2 = muu3*sqrtS


	DL1 = (g0R*(-(g1R*g2R*Ml2*muu4*sqrtS*(xc0 + xc1)) + g1L*Ml1*(2*g2R*Ml2*xc0 - g2L*muu3*sqrtS*(xc0 + xc2))) + 
     -    g0L*(-(g1R*g2R*Ml1*muu4*sqrtS*xc1) + g1L*(-(g2L*Ml2*muu3*sqrtS*xc2) + 
     -          g2R*(2*xb0 + 2*Ml0**2*xc0 + sqrtS**2*
     -              ((-1 + muu3**2 + muu4**2)*xc0 + (-1 + muu3**2 + 3*muu4**2)*xc1 + (-1 + 3*muu3**2 + muu4**2)*xc2)))))/
     -  (8.*Pi**2)

	DR1 = (g0L*(-(g1L*g2L*Ml2*muu4*sqrtS*(xc0 + xc1)) + g1R*Ml1*(2*g2L*Ml2*xc0 - g2R*muu3*sqrtS*(xc0 + xc2))) + 
     -    g0R*(-(g1L*g2L*Ml1*muu4*sqrtS*xc1) + g1R*(-(g2R*Ml2*muu3*sqrtS*xc2) + 
     -          g2L*(2*xb0 + 2*Ml0**2*xc0 + sqrtS**2*
     -              ((-1 + muu3**2 + muu4**2)*xc0 + (-1 + muu3**2 + 3*muu4**2)*xc1 + (-1 + 3*muu3**2 + muu4**2)*xc2)))))/
     -  (8.*Pi**2)
	
cc	write(*,*)' --- LOOP GLUON --- '
cc	write(*,*)'DL1 =',DL1
cc	write(*,*)'DR1 =',DR1
cc	
cc	write(*,*)' --- Cget gluon--- '
cc      write(*,*)'xc0 =',xc0
cc      write(*,*)'xc1 =',xc1
cc      write(*,*)'xc2 =',xc2
cc      write(*,*)'xc00 =',xc00
cc      write(*,*)'xc11 =',xc11
cc      write(*,*)'xc12 =',xc12
cc      write(*,*)'xc22 =',xc22
      
	
!******************************************************************************************	
!*********** Gluino loop *******
!*******************************
	
	
	do k=1,2
	 do l=1,2

	 g0  = HSfSf(i,k,l,tt,gen)
	 g1L = -sqrt2*GS*USf(l,1,att,gen)
	 g1R = sqrt2*GS*USf(l,2,att,gen)
	 g2L = sqrt2*GS*USf(k,2,tt,gen)
	 g2R = -sqrt2*GS*USf(k,1,tt,gen)
	 
	 Ml0 = MGl
	 Ml1 = MSf(l,att,gen)
	 Ml2 = MSf(k,tt,gen)
	 
	 call myCget(muu4**2*sqrtS**2,sqrtS**2,muu3**2*sqrtS**2,MGl**2,MSf2(l,att,gen),MSf2(k,tt,gen))
	 call myBget(sqrtS**2,MSf2(l,att,gen),MSf2(k,tt,gen))

	 DL1 = DL1 + (g0*(-(g1L*g2L*Ml0*xc0) + g1R*g2L*muu4*sqrtS*xc1 + g1L*g2R*muu3*sqrtS*xc2))/(16.*Pi**2)
	 DR1 = DR1 + (g0*(-(g1R*g2R*Ml0*xc0) + g1L*g2R*muu4*sqrtS*xc1 + g1R*g2L*muu3*sqrtS*xc2))/(16.*Pi**2)

	 enddo
	enddo

cc	write(*,*)' --- TREE --- '
cc	write(*,*)'DL =',DL(i)
cc	write(*,*)'DR =',DR(i)
cc	write(*,*)' --- LOOP --- '
cc	write(*,*)'DL1 =',DL1
cc	write(*,*)'DR1 =',DR1

	
! Higgs-loop interference with photon tree-level

	if (icha.eq.2) then

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex - 8*ApL*BpL*CF*(DL1 + DR1)*muu1*muu3*NC*sqrtS**2*(ured - tred)*(CL(i) + CR(i))/(dreal(dconjg(DenSmH(i))))

	endif

! Higgs-loop interference with Z/W tree-level

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex - 4*CF*NC*sqrtS**4*(-(
     &           muu4*(muu1*(muu2**2 + muu3**2 - ured)*(CR(i)*dconjg(AL) + CL(i)*dconjg(AR)) - 
     &           muu2*(muu1**2 + muu3**2 - tred)*(CL(i)*dconjg(AL) + CR(i)*dconjg(AR)))*(DR1*dconjg(BL) + DL1*dconjg(BR))) + 
     &           muu3*(muu1*(muu2**2 + muu4**2 - tred)*(CR(i)*dconjg(AL) + CL(i)*dconjg(AR)) - 
     &           muu2*(muu1**2 + muu4**2 - ured)*(CL(i)*dconjg(AL) + CR(i)*dconjg(AR)))*(DL1*dconjg(BL) + DR1*dconjg(BR)))
     &  /(dreal(DenSmV*dconjg(DenSmH(i))))
 
	tempZ = tempZ + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! Higgs-loop interference with Higgs tree-level

       do j= jmin,jmax

	 ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + 2*CF*NC*sqrtS**4*(2*muu1*muu2*(CR(i)*dconjg(CL(j)) + CL(i)*dconjg(CR(j))) + 
     &      (-1 + muu1**2)*(CL(i)*dconjg(CL(j)) + CR(i)*dconjg(CR(j))) + muu2**2*(CL(i)*dconjg(CL(j)) + CR(i)*dconjg(CR(j))))*
     &    (2*muu3*muu4*(DR1*dconjg(DL(j)) + DL1*dconjg(DR(j))) + (-1 + muu3**2)*(DL1*dconjg(DL(j)) + DR1*dconjg(DR(j))) + 
     &      muu4**2*(DL1*dconjg(DL(j)) + DR1*dconjg(DR(j))))/(DenSmH(i)*dconjg(DenSmH(j)))

       enddo
	
	tempH = tempH + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! Higgs-loop interference with T- and U-channel tree-level

	do j=1,2
	 
	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + dreal(NC*CF*sqrtS**2*(dconjg(CL(i))*(2*dconjg(DR1)*
     &          (muu1*GR(indi,j,tt,gen)*(-(muu3*(muu2**2 + muu4**2 - tred)*HL(indj,j,tt,gen)) + 
     &               muu2*(-1 + muu3**2 + muu4**2)*HR(indj,j,tt,gen)) + 
     &            muu4*GL(indi,j,tt,gen)*((-1 + muu1**2 + muu2**2)*muu3*HL(indj,j,tt,gen) - 
     &               muu2*(muu1**2 + muu3**2 - tred)*HR(indj,j,tt,gen))) + 
     &         dconjg(DL1)*(2*muu1*muu4*GR(indi,j,tt,gen)*
     &             ((muu2**2 + muu3**2 - ured)*HL(indj,j,tt,gen) + 2*muu2*muu3*HR(indj,j,tt,gen)) + 
     &            GL(indi,j,tt,gen)*(((-1 + muu3**2 + muu4**2 - tred)*(-1 + tred) + muu1**2*(-1 + 2*muu3**2 + tred - ured) + 
     &                  muu2**2*(-1 + 2*muu4**2 + tred - ured) - (muu3**2 + muu4**2)*ured + ured**2)*HL(indj,j,tt,gen) + 
     &               2*muu2*muu3*(muu1**2 + muu4**2 - ured)*HR(indj,j,tt,gen)))) + 
     &      dconjg(CR(i))*(2*dconjg(DL1)*(muu4*GR(indi,j,tt,gen)*
     &             (-(muu2*(muu1**2 + muu3**2 - tred)*HL(indj,j,tt,gen)) + (-1 + muu1**2 + muu2**2)*muu3*HR(indj,j,tt,gen)) + 
     &            muu1*GL(indi,j,tt,gen)*(muu2*(-1 + muu3**2 + muu4**2)*HL(indj,j,tt,gen) - 
     &               muu3*(muu2**2 + muu4**2 - tred)*HR(indj,j,tt,gen))) + 
     &         dconjg(DR1)*(2*muu1*muu4*GL(indi,j,tt,gen)*
     &             (2*muu2*muu3*HL(indj,j,tt,gen) + (muu2**2 + muu3**2 - ured)*HR(indj,j,tt,gen)) + 
     &            GR(indi,j,tt,gen)*(2*muu2*muu3*(muu1**2 + muu4**2 - ured)*HL(indj,j,tt,gen) + 
     &               ((-1 + muu3**2 + muu4**2 - tred)*(-1 + tred) + muu1**2*(-1 + 2*muu3**2 + tred - ured) + 
     &                  muu2**2*(-1 + 2*muu4**2 + tred - ured) - (muu3**2 + muu4**2)*ured + ured**2)*HR(indj,j,tt,gen))))))/
     &  (dconjg(DenSmH(i))*(-ured + MSf(j,itypeU,gen)**2/sqrtS**2))
	 
	tempu = tempu + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + dreal(CF*NC*sqrtS**2*(dconjg(CL(i))*(dconjg(DL1)*(2*muu1*muu3*ER(indi,j,tt,gen)*
     &             ((muu2**2 + muu4**2 - tred)*FL(indj,j,tt,gen) + 2*muu2*muu4*FR(indj,j,tt,gen)) + 
     &            ELc(indi,j,tt,gen)*((1 - muu3**2 - muu4**2 - muu3**2*tred - muu4**2*tred + tred**2 + muu3**2*ured + 
     &                  muu4**2*ured - ured**2 + muu2**2*(-1 + 2*muu3**2 - tred + ured) + muu1**2*(-1 + 2*muu4**2 - tred + ured))
     &                 *FL(indj,j,tt,gen) + 2*muu2*muu4*(muu1**2 + muu3**2 - tred)*FR(indj,j,tt,gen))) + 
     &         2*dconjg(DR1)*(muu1*ER(indi,j,tt,gen)*
     &             (-(muu4*(muu2**2 + muu3**2 - ured)*FL(indj,j,tt,gen)) + muu2*(-1 + muu3**2 + muu4**2)*FR(indj,j,tt,gen)) + 
     &            muu3*ELc(indi,j,tt,gen)*((-1 + muu1**2 + muu2**2)*muu4*FL(indj,j,tt,gen) - 
     &               muu2*(muu1**2 + muu4**2 - ured)*FR(indj,j,tt,gen)))) + 
     &      dconjg(CR(i))*(2*dconjg(DL1)*(muu3*ER(indi,j,tt,gen)*
     &             (-(muu2*(muu1**2 + muu4**2 - ured)*FL(indj,j,tt,gen)) + (-1 + muu1**2 + muu2**2)*muu4*FR(indj,j,tt,gen)) + 
     &            muu1*ELc(indi,j,tt,gen)*(muu2*(-1 + muu3**2 + muu4**2)*FL(indj,j,tt,gen) - 
     &               muu4*(muu2**2 + muu3**2 - ured)*FR(indj,j,tt,gen))) + 
     &         dconjg(DR1)*(2*muu1*muu3*ELc(indi,j,tt,gen)*
     &             (2*muu2*muu4*FL(indj,j,tt,gen) + (muu2**2 + muu4**2 - tred)*FR(indj,j,tt,gen)) + 
     &            ER(indi,j,tt,gen)*(2*muu2*muu4*(muu1**2 + muu3**2 - tred)*FL(indj,j,tt,gen) + 
     &               (1 - muu3**2 - muu4**2 - muu3**2*tred - muu4**2*tred + tred**2 + muu3**2*ured + muu4**2*ured - ured**2 + 
     &                  muu2**2*(-1 + 2*muu3**2 - tred + ured) + muu1**2*(-1 + 2*muu4**2 - tred + ured))*FR(indj,j,tt,gen))))))/
     &  (DenSmH1(i)*(-tred + MSf(j,itypeT,gen)**2/sqrtS**2))
	 
	tempt = tempt + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	 
	 enddo

	enddo

cc	write(*,*)'ChiChi H Vertex mit Z = ', tempZ
cc	write(*,*)'ChiChi H Vertex mit H = ', tempH
cc	write(*,*)'ChiChi H Vertex mit t = ', tempt
cc	write(*,*)'ChiChi H Vertex mit u = ', tempu


!*************************************************
!*********** Loop in the T&U-channels ************ 
!*************************************************
	
	
	do j=1,2

! old T-channel p2k1 vertex --> new U-channel p2k1 vertex

! Gluon
	call myCget(muu3**2*sqrtS**2,ured*sqrtS**2,muu2**2*sqrtS**2,muu3**2*sqrtS**2,0d0,MSf2(j,itypeU,gen))
	call myBget(ured*sqrtS**2,0d0,MSf2(j,itypeU,gen))

	 M0  = muu3*sqrtS
	 g0  = -GS
	 g3L = -GS
	 g4L = HL(indj,j,tt,gen)
	 g4R = HR(indj,j,tt,gen)

	 HL1 = (g0*g3L*(-2*g4R*muu2*muu3*sqrtS**2*(xc0 + xc1 + xc2) + g4L*(xb0 + 2*sqrtS**2*(muu3**2*xc0 - ured*xc1 + muu2**2*(xc1 + xc2)))))/(16.*Pi**2)
	 HR1 = (g0*g3L*(-2*g4L*muu2*muu3*sqrtS**2*(xc0 + xc1 + xc2) + g4R*(xb0 + 2*sqrtS**2*(muu3**2*xc0 - ured*xc1 + muu2**2*(xc1 + xc2)))))/(16.*Pi**2)
	 
	 tempGluonL = HL1
	 tempGluonR = HR1
cc	 write(*,*)'HL1 Gluon = ', HL1
cc	 write(*,*)'HR1 Gluon = ', HR1

cc      write(*,*)' --- Cget gluon p2k1--- '
cc      write(*,*)'xc0 =',xc0
cc      write(*,*)'xc1 =',xc1
cc      write(*,*)'xc2 =',xc2
cc      write(*,*)'xc00 =',xc00
cc      write(*,*)'xc11 =',xc11
cc      write(*,*)'xc12 =',xc12
cc      write(*,*)'xc22 =',xc22


	 do k=1,2

! Gluino

	  call myCget(muu3**2*sqrtS**2,ured*sqrtS**2,muu2**2*sqrtS**2,MSf2(k,tt,gen),MGl**2,Mf2(itypeU,gen))
	  call myBget(ured*sqrtS**2,MGl**2,Mf2(itypeU,gen))

	  g0L = sqrt2*GS*USf(j,2,itypeU,gen)
	  g0R = -sqrt2*GS*USf(j,1,itypeU,gen)

	  g3L = sqrt2*GS*USf(k,2,tt,gen)
	  g3R = -sqrt2*GS*USf(k,1,tt,gen)

	  g4L = dconjg(HR(indj,k,itypeU,gen))
	  g4R = dconjg(HL(indj,k,itypeU,gen))
	  
	  if(icha.eq.2) then
	    if (tt.eq.4) then
! g4L,g4R are set to FL,FR but these are zero in this case -> explicit forms
		
		 g4L = -EL/SW*UChaC(indj,1)*USf(k,1,itypeU+1,gen) + EL*Mf(itypeU+1,gen)/(SW*sqrt2*MW*CB)*UChaC(indj,2)*USf(k,2,itypeU+1,gen)
		 g4R = EL*Mf(itypeU,gen)/(SW*sqrt2*MW*SB)*VCha(indj,2)*USf(k,1,itypeU+1,gen)
		 
		else 

		 g4L = 0d0
		 g4R = 0d0

		endif
	  endif

	  HL1 = HL1 - (g0L*sqrtS**2*(g3L*MGlr*(g4R*muu2*(xc0 + xc2) + (g4L*xc0*Mf(itypeU,gen))/sqrtS) + 
     -        g3R*muu3*(g4R*muu2*(xc0 + xc1 + xc2) + (g4L*(xc0 + xc1)*Mf(itypeU,gen))/sqrtS)) + 
     -     g0R*(g3R*g4L*MGlr*muu3*sqrtS**2*xc1 + g3L*
     -         (g4L*(xb0 + sqrtS**2*(muu3**2*xc1 + muu2**2*xc2)) + g4R*muu2*sqrtS*xc2*Mf(itypeU,gen))) + 
     -     g0R*g3L*g4L*xc0*MSf2(k,tt,gen))/(16.*Pi**2)

	  HR1 = HR1 -(g0R*sqrtS**2*(g3R*MGlr*(g4L*muu2*(xc0 + xc2) + (g4R*xc0*Mf(itypeU,gen))/sqrtS) + 
     -        g3L*muu3*(g4L*muu2*(xc0 + xc1 + xc2) + (g4R*(xc0 + xc1)*Mf(itypeU,gen))/sqrtS)) + 
     -     g0L*(g3L*g4R*MGlr*muu3*sqrtS**2*xc1 + g3R*
     -         (g4R*(xb0 + sqrtS**2*(muu3**2*xc1 + muu2**2*xc2)) + g4L*muu2*sqrtS*xc2*Mf(itypeU,gen))) + 
     -     g0L*g3R*g4R*xc0*MSf2(k,tt,gen))/(16.*Pi**2)

	enddo
	
!	write(*,*)' --- TREE --- '
!	write(*,*)'HL =',HL(indj,j,tt,gen)
!	write(*,*)'HR =',HR(indj,j,tt,gen)
!	write(*,*)' --- LOOP --- '
!	write(*,*)'HL1 =',HL1
!	write(*,*)'HR1 =',HR1
	
cc	write(*,*)'HL1 Gluino = ', HL1 - tempGluonL
cc	write(*,*)'HR1 Gluino = ', HR1 - tempGluonR


! old T-channel p1k2 vertex ==> new U-channel p1k2 vertex

! Gluon

	call myCget(muu4**2*sqrtS**2,ured*sqrtS**2,muu1**2*sqrtS**2,muu4**2*sqrtS**2,0d0,MSf2(j,itypeU,gen))
	call myBget(ured*sqrtS**2,0d0,MSf2(j,itypeU,gen))

	 M0  = muu4*sqrtS
	 g0  = -GS
	 g3L = -GS
	 g4L = GL(indi,j,tt,gen)
	 g4R = GR(indi,j,tt,gen)

	 GL1 = (g0*g3L*(-2*g4R*muu1*muu4*sqrtS**2*(xc0 + xc1 + xc2) + g4L*(xb0 + 2*sqrtS**2*(muu4**2*xc0 - ured*xc1 + muu1**2*(xc1 + xc2)))))/(16.*Pi**2)
	 GR1 = (g0*g3L*(-2*g4L*muu1*muu4*sqrtS**2*(xc0 + xc1 + xc2) + g4R*(xb0 + 2*sqrtS**2*(muu4**2*xc0 - ured*xc1 + muu1**2*(xc1 + xc2)))))/(16.*Pi**2)	
	
	 tempGluonL = GL1
	 tempGluonR = GR1
cc	 write(*,*)'GL1 Gluon = ', GL1
cc	 write(*,*)'GR1 Gluon = ', GR1

cc      write(*,*)' --- Cget gluon p1k2--- '
cc      write(*,*)'xc0 =',xc0
cc      write(*,*)'xc1 =',xc1
cc      write(*,*)'xc2 =',xc2
cc      write(*,*)'xc00 =',xc00
cc      write(*,*)'xc11 =',xc11
cc      write(*,*)'xc12 =',xc12
cc      write(*,*)'xc22 =',xc22
	 
! Gluino

	 do k=1,2

	  call myCget(muu4**2*sqrtS**2,ured*sqrtS**2,muu1**2*sqrtS**2,MSf2(k,att,gen),MGl**2,Mf2(itypeU,gen))
	  call myBget(ured*sqrtS**2,MGl**2,Mf2(itypeU,gen))

	  g0R = sqrt2*GS*USf(j,2,itypeU,gen)
	  g0L = -sqrt2*GS*USf(j,1,itypeU,gen)

	  g3R = sqrt2*GS*USf(k,2,att,gen)
	  g3L = -sqrt2*GS*USf(k,1,att,gen)

	  g4L = ELc(indi,k,itypeU,gen)
	  g4R = ER(indi,k,itypeU,gen)

	  if(icha.eq.2) then
	    if (tt.eq.4) then
! g4L,g4R are set to ELc,ER but these are zero in this case -> explicit forms
		
	     g4L =  EL*Mf(itypeU,gen)/(SW*sqrt2*MW*SB)*VChaC(indi,2)*USf(k,1,itypeU+1,gen) 
	     g4R =  -EL/SW*UCha(indi,1)*USf(k,1,itypeU+1,gen) + EL*Mf(itypeU+1,gen)/(SW*sqrt2*MW*CB)*UCha(indi,2)*USf(k,2,itypeU+1,gen) 
		 
		else 

		 g4L = 0d0
		 g4R = 0d0

		endif
	  endif


	  GL1 = GL1 - (g0L*sqrtS**2*(g3L*MGlr*(g4R*muu1*(xc0 + xc2) + (g4L*xc0*Mf(itypeU,gen))/sqrtS) + 
     -        g3R*muu4*(g4R*muu1*(xc0 + xc1 + xc2) + (g4L*(xc0 + xc1)*Mf(itypeU,gen))/sqrtS)) + 
     -     g0R*(g3R*g4L*MGlr*muu4*sqrtS**2*xc1 + g3L*
     -         (g4L*(xb0 + sqrtS**2*(muu4**2*xc1 + muu1**2*xc2)) + g4R*muu1*sqrtS*xc2*Mf(itypeU,gen))) + 
     -     g0R*g3L*g4L*xc0*MSf2(k,att,gen))/(16.*Pi**2)

	  GR1 = GR1 - (g0R*sqrtS**2*(g3R*MGlr*(g4L*muu1*(xc0 + xc2) + (g4R*xc0*Mf(itypeU,gen))/sqrtS) + 
     -        g3L*muu4*(g4L*muu1*(xc0 + xc1 + xc2) + (g4R*(xc0 + xc1)*Mf(itypeU,gen))/sqrtS)) + 
     -     g0L*(g3L*g4R*MGlr*muu4*sqrtS**2*xc1 + g3R*
     -         (g4R*(xb0 + sqrtS**2*(muu4**2*xc1 + muu1**2*xc2)) + g4L*muu1*sqrtS*xc2*Mf(itypeU,gen))) + 
     -     g0L*g3R*g4R*xc0*MSf2(k,att,gen))/(16.*Pi**2)

	enddo
	
!	write(*,*)' --- TREE --- '
!	write(*,*)'GL =',GL(indi,j,tt,gen)
!	write(*,*)'GR =',GR(indi,j,tt,gen)
!	write(*,*)' --- LOOP --- '
!	write(*,*)'GL1 =',GL1
!	write(*,*)'GR1 =',GR1
	
cc	write(*,*)'GL1 Gluino = ', GL1 - tempGluonL
cc	write(*,*)'GR1 Gluino = ', GR1 - tempGluonR
	
! old U-channel p2k2 vertex ==> new T-channel p2k2 vertex

! Gluon
	call myCget(muu4**2*sqrtS**2,tred*sqrtS**2,muu2**2*sqrtS**2,muu4**2*sqrtS**2,0d0,MSf2(j,itypeT,gen))
	call myBget(tred*sqrtS**2,0d0,MSf2(j,itypeT,gen))
	  
	 M0  = muu4*sqrtS 
	 g0  = -GS
	 g3L = -GS
	 g4L = FL(indj,j,tt,gen)
	 g4R = FR(indj,j,tt,gen)

	 FL1 = (g0*g3L*(-2*g4R*muu2*muu4*sqrtS**2*(xc0 + xc1 + xc2) + g4L*(xb0 + 2*sqrtS**2*(muu4**2*xc0 - tred*xc1 + muu2**2*(xc1 + xc2)))))/(16.*Pi**2)
	 FR1 = (g0*g3L*(-2*g4L*muu2*muu4*sqrtS**2*(xc0 + xc1 + xc2) + g4R*(xb0 + 2*sqrtS**2*(muu4**2*xc0 - tred*xc1 + muu2**2*(xc1 + xc2)))))/(16.*Pi**2)
	
	 tempGluonL = FL1
	 tempGluonR = FR1
cc	 write(*,*)'FL1 Gluon = ', FL1
cc	 write(*,*)'FR1 Gluon = ', FR1

cc      write(*,*)' --- Cget gluon p2k2--- '
cc      write(*,*)'xc0 =',xc0
cc      write(*,*)'xc1 =',xc1
cc      write(*,*)'xc2 =',xc2
cc      write(*,*)'xc00 =',xc00
cc      write(*,*)'xc11 =',xc11
cc      write(*,*)'xc12 =',xc12
cc      write(*,*)'xc22 =',xc22

	 
! Gluino
	 do k=1,2

	  call myCget(muu4**2*sqrtS**2,tred*sqrtS**2,muu2**2*sqrtS**2,MSf2(k,att,gen),MGl**2,Mf2(itypeT,gen))
	  call myBget(tred*sqrtS**2,MGl**2,Mf2(itypeT,gen))

	  g0R = sqrt2*GS*USf(j,2,itypeT,gen)
	  g0L = -sqrt2*GS*USf(j,1,itypeT,gen)

	  g3R = sqrt2*GS*USf(k,2,att,gen)
	  g3L = -sqrt2*GS*USf(k,1,att,gen)

	  g4L = dconjg(FR(indj,k,tt,gen))
	  g4R = dconjg(FL(indj,k,tt,gen))
	  
	  if(icha.eq.2) then
	    if (tt.eq.3) then
! g4L, g4R are set to HL, HR but these are zero in this case -> explicit forms
		
		g4L =  EL*Mf(itypeT,gen)/(SW*sqrt2*MW*CB)*UChaC(indj,2)*USf(k,1,itypeT-1,gen)
		g4R  =  -EL/SW*VCha(indj,1)*USf(k,1,itypeT-1,gen) + EL*Mf(itypeT-1,gen)/(SW*sqrt2*MW*SB)*VCha(indj,2)*USf(k,2,itypeT-1,gen)
		
		else 

		 g4L = 0d0
		 g4R = 0d0

		endif
	  endif
	  

	  FL1 = FL1 - (g0L*sqrtS**2*(g3L*MGlr*(g4R*muu2*(xc0 + xc2) + (g4L*xc0*Mf(itypeT,gen))/sqrtS) + 
     -        g3R*muu4*(g4R*muu2*(xc0 + xc1 + xc2) + (g4L*(xc0 + xc1)*Mf(itypeT,gen))/sqrtS)) + 
     -     g0R*(g3R*g4L*MGlr*muu4*sqrtS**2*xc1 + g3L*
     -         (g4L*(xb0 + sqrtS**2*(muu4**2*xc1 + muu2**2*xc2)) + g4R*muu2*sqrtS*xc2*Mf(itypeT,gen))) + 
     -     g0R*g3L*g4L*xc0*MSf2(k,att,gen))/(16.*Pi**2)

	  FR1 = FR1 - (g0R*sqrtS**2*(g3R*MGlr*(g4L*muu2*(xc0 + xc2) + (g4R*xc0*Mf(itypeT,gen))/sqrtS) + 
     -        g3L*muu4*(g4L*muu2*(xc0 + xc1 + xc2) + (g4R*(xc0 + xc1)*Mf(itypeT,gen))/sqrtS)) + 
     -     g0L*(g3L*g4R*MGlr*muu4*sqrtS**2*xc1 + g3R*
     -         (g4R*(xb0 + sqrtS**2*(muu4**2*xc1 + muu2**2*xc2)) + g4L*muu2*sqrtS*xc2*Mf(itypeT,gen))) + 
     -     g0L*g3R*g4R*xc0*MSf2(k,att,gen))/(16.*Pi**2)

	enddo
	
!	write(*,*)' --- TREE --- '
!	write(*,*)'FL =',FL(indj,j,tt,gen)
!	write(*,*)'FR =',FR(indj,j,tt,gen)
!	write(*,*)' --- LOOP --- '
!	write(*,*)'FL1 =',FL1
!	write(*,*)'FR1 =',FR1
	
cc	write(*,*)'FL1 Gluino = ', FL1 - tempGluonL
cc	write(*,*)'FR1 Gluino = ', FR1 - tempGluonR
cc      write(*,*)' --- Cget gluon p1k1--- '
cc      write(*,*)'xc0 =',xc0
cc      write(*,*)'xc1 =',xc1
cc      write(*,*)'xc2 =',xc2
cc      write(*,*)'xc00 =',xc00
cc      write(*,*)'xc11 =',xc11
cc      write(*,*)'xc12 =',xc12
cc      write(*,*)'xc22 =',xc22

	
	
! old U-channel p1k1 vertex ==> new T-channel p1k1 vertex

! Gluon
	
	call myCget(muu1**2*sqrtS**2,tred*sqrtS**2,muu3**2*sqrtS**2,muu3**2*sqrtS**2,MSf2(j,itypeT,gen),0d0)
	call myBget(tred*sqrtS**2,MSf2(j,itypeT,gen),0d0)
	  
	 g0  = -GS
	 g1L = ELc(indi,j,tt,gen)
	 g1R = ER(indi,j,tt,gen)
	 g2L = -GS
	 g2R = -GS

	 Ml0  = muu3*sqrtS 
	 Ml1  = MSf(j,itypeT,gen) 
	 Ml2  = 0d0 	

	 EL1 = (g0*(-(g1R*muu1*sqrtS*(g2R*Ml0*(2*xc0 + xc1) + g2L*muu3*sqrtS*(xc1 + 2*xc2))) + 
     -      g1L*(g2L*Ml0*muu3*sqrtS*(xc0 - xc2) + g2R*
     -          (xb0 + Ml0**2*xc0 + sqrtS**2*((muu3**2 - 2*tred)*xc2 + 2*muu1**2*(xc1 + xc2))))))/(16.*Pi**2)
	
	 ER1 = (g0*(-(g1L*muu1*sqrtS*(g2L*Ml0*(2*xc0 + xc1) + g2R*muu3*sqrtS*(xc1 + 2*xc2))) + 
     -      g1R*(g2R*Ml0*muu3*sqrtS*(xc0 - xc2) + g2L*
     -          (xb0 + Ml0**2*xc0 + sqrtS**2*((muu3**2 - 2*tred)*xc2 + 2*muu1**2*(xc1 + xc2))))))/(16.*Pi**2)
	
	 tempGluonL = EL1
	 tempGluonR = ER1
cc	 write(*,*)'EL1 Gluon = ', EL1
cc	 write(*,*)'ER1 Gluon = ', ER1
	 
! Gluino

	 do k=1,2

	  call myCget(muu3**2*sqrtS**2,tred*sqrtS**2,muu1**2*sqrtS**2,MSf2(k,tt,gen),MGl**2,Mf2(itypeT,gen))
	  call myBget(tred*sqrtS**2,MGl**2,Mf2(itypeT,gen))

	  g0L = sqrt2*GS*USf(j,2,itypeT,gen)
	  g0R = -sqrt2*GS*USf(j,1,itypeT,gen)

	  g3L = sqrt2*GS*USf(k,2,tt,gen)
	  g3R = -sqrt2*GS*USf(k,1,tt,gen)

	  g4L = GL(indi,k,tt,gen)
	  g4R = GR(indi,k,tt,gen)
	  
	  if(icha.eq.2) then
	    if (tt.eq.3) then
! g4L, g4R are set to GL, GR but these are zero in this case -> explicit forms
		
	     g4L =  -EL/SW*VChaC(indi,1)*USf(k,1,itypeT-1,gen) + EL*Mf(itypeT-1,gen)/(SW*sqrt2*MW*SB)*VChaC(indi,2)*USf(k,2,itypeT-1,gen)
	     g4R =  EL*Mf(itypeT,gen)/(SW*sqrt2*MW*CB)*UCha(indi,2)*USf(k,1,itypeT-1,gen)

		else 

		 g4L = 0d0
		 g4R = 0d0

		endif
	  endif
	  

	  EL1 = EL1 - (g0L*sqrtS**2*(g3L*MGlr*(g4R*muu1*(xc0 + xc2) + (g4L*xc0*Mf(itypeT,gen))/sqrtS) + 
     -        g3R*muu3*(g4R*muu1*(xc0 + xc1 + xc2) + (g4L*(xc0 + xc1)*Mf(itypeT,gen))/sqrtS)) + 
     -     g0R*(g3R*g4L*MGlr*muu3*sqrtS**2*xc1 + g3L*
     -         (g4L*(xb0 + sqrtS**2*(muu3**2*xc1 + muu1**2*xc2)) + g4R*muu1*sqrtS*xc2*Mf(itypeT,gen))) + 
     -     g0R*g3L*g4L*xc0*MSf2(k,tt,gen))/(16.*Pi**2)

	  ER1 = ER1 - (g0R*sqrtS**2*(g3R*MGlr*(g4L*muu1*(xc0 + xc2) + (g4R*xc0*Mf(itypeT,gen))/sqrtS) + 
     -        g3L*muu3*(g4L*muu1*(xc0 + xc1 + xc2) + (g4R*(xc0 + xc1)*Mf(itypeT,gen))/sqrtS)) + 
     -     g0L*(g3L*g4R*MGlr*muu3*sqrtS**2*xc1 + g3R*
     -         (g4R*(xb0 + sqrtS**2*(muu3**2*xc1 + muu1**2*xc2)) + g4L*muu1*sqrtS*xc2*Mf(itypeT,gen))) + 
     -     g0L*g3R*g4R*xc0*MSf2(k,tt,gen))/(16.*Pi**2)

	enddo
	
cc	write(*,*)'EL1 Gluino = ', EL1 - tempGluonL
cc	write(*,*)'ER1 Gluino = ', ER1 - tempGluonR

cc	write(*,*)' --- TREE --- '
cc	write(*,*)'EL =',ELc(indi,j,tt,gen)
cc	write(*,*)'ER =',ER(indi,j,tt,gen)
cc	write(*,*)' --- LOOP --- '
cc	write(*,*)'EL1 =',EL1
cc	write(*,*)'ER1 =',ER1
	

! T-channel-loop interference photon tree-level

	if (icha.eq.2) then
	
! (p2k2)-vertex with photon

	 ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + 4*ApL*BpL*CF*NC*(tred**2*(FR1*ELc(indi,j,tt,gen) + FL1*ER(indi,j,tt,gen)) + 
     -      (1 + muu1**2 + muu3**2)*(muu1**2*(FR1*ELc(indi,j,tt,gen) + FL1*ER(indi,j,tt,gen)) + 
     -         muu3**2*(FR1*ELc(indi,j,tt,gen) + FL1*ER(indi,j,tt,gen)) + 2*muu1*muu3*(FL1*ELc(indi,j,tt,gen) + FR1*ER(indi,j,tt,gen))) - 
     -      2*tred*(FR1*((muu1**2 + muu3**2)*ELc(indi,j,tt,gen) + muu1*muu3*ER(indi,j,tt,gen)) + 
     -         FL1*(muu1*muu3*ELc(indi,j,tt,gen) + (muu1**2 + muu3**2)*ER(indi,j,tt,gen))))/(tred - MSf(j,itypeT,gen)**2/sqrtS**2)

! (p1k1)-vertex with photon

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + 4*ApL*BpL*CF*NC*(tred**2*(ER1*FL(indj,j,tt,gen) + EL1*FR(indj,j,tt,gen)) - 
     -      2*tred*((EL1*muu1*muu3 + ER1*(muu1**2 + muu3**2))*FL(indj,j,tt,gen) + 
     -         (ER1*muu1*muu3 + EL1*(muu1**2 + muu3**2))*FR(indj,j,tt,gen)) + 
     -      (1 + muu1**2 + muu3**2)*(muu1**2*(ER1*FL(indj,j,tt,gen) + EL1*FR(indj,j,tt,gen)) + 
     -         muu3**2*(ER1*FL(indj,j,tt,gen) + EL1*FR(indj,j,tt,gen)) + 2*muu1*muu3*(EL1*FL(indj,j,tt,gen) + ER1*FR(indj,j,tt,gen))))/
     -  (tred - MSf(j,itypeT,gen)**2/sqrtS**2)


! U-channel-loop interference photon tree-level

! (p2k1)-vertex with photon

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + (4*ApL*BpL*CF*NC*(ured**2*(HR1*GL(indi,j,tt,gen) + HL1*GR(indi,j,tt,gen)) +     
     &      (1 + muu1**2 + muu3**2)*(muu1**2*(HR1*GL(indi,j,tt,gen) + HL1*GR(indi,j,tt,gen)) +     
     &         muu3**2*(HR1*GL(indi,j,tt,gen) + HL1*GR(indi,j,tt,gen)) + 2*muu1*muu3*(HL1*GL(indi,j,tt,gen) + HR1*GR(indi,j,tt,gen))) -     
     &      2*ured*(HR1*((muu1**2 + muu3**2)*GL(indi,j,tt,gen) + muu1*muu3*GR(indi,j,tt,gen)) +     
     &         HL1*(muu1*muu3*GL(indi,j,tt,gen) + (muu1**2 + muu3**2)*GR(indi,j,tt,gen)))))/(-ured + MSf(j,itypeU,gen)**2/sqrtS**2)

! (p1k2)-vertex with photon

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + (4*ApL*BpL*CF*NC*(ured**2*(GR1*HL(indj,j,tt,gen) + GL1*HR(indj,j,tt,gen)) -     
     &      2*ured*((GL1*muu1*muu3 + GR1*(muu1**2 + muu3**2))*HL(indj,j,tt,gen) + (GR1*muu1*muu3 + GL1*(muu1**2 + muu3**2))*HR(indj,j,tt,gen)) +     
     &      (1 + muu1**2 + muu3**2)*(muu1**2*(GR1*HL(indj,j,tt,gen) + GL1*HR(indj,j,tt,gen)) +     
     &         muu3**2*(GR1*HL(indj,j,tt,gen) + GL1*HR(indj,j,tt,gen)) + 2*muu1*muu3*(GL1*HL(indj,j,tt,gen) + GR1*HR(indj,j,tt,gen)))))/
     &  	(-ured + MSf(j,itypeU,gen)**2/sqrtS**2)

	endif


! U-channel-loop interference Z/W tree-level

! (p2k1)-vertex with Z/W
	
	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + dreal(4*NC*CF*sqrtS**2*(dconjg(AR)*(BR*(2*muu1*muu4*GL(indi,j,tt,gen) + (muu1**2 + muu4**2 - ured)*GR(indi,j,tt,gen))*
     &          ((muu2**2 + muu3**2 - ured)*HL1 + 2*muu2*muu3*HR1) + 
     &         BL*(muu1*GL(indi,j,tt,gen)*(muu3*(muu2**2 + muu4**2 - tred)*HL1 - 
     &               muu2*(-1 + muu3**2 + muu4**2)*HR1) + 
     &            muu4*GR(indi,j,tt,gen)*(-((-1 + muu1**2 + muu2**2)*muu3*HL1) + 
     &               muu2*(muu1**2 + muu3**2 - tred)*HR1))) + 
     &      dconjg(AL)*(BL*((muu1**2 + muu4**2 - ured)*GL(indi,j,tt,gen) + 2*muu1*muu4*GR(indi,j,tt,gen))*
     &          (2*muu2*muu3*HL1 + (muu2**2 + muu3**2 - ured)*HR1) + 
     &         BR*(muu4*GL(indi,j,tt,gen)*(muu2*(muu1**2 + muu3**2 - tred)*HL1 - 
     &               (-1 + muu1**2 + muu2**2)*muu3*HR1) + 
     &            muu1*GR(indi,j,tt,gen)*(-(muu2*(-1 + muu3**2 + muu4**2)*HL1) + 
     &               muu3*(muu2**2 + muu4**2 - tred)*HR1)))))/(DenSmV1*(-ured + MSf(j,itypeU,gen)**2/sqrtS**2))
	 
	tempZ21 = tempZ21 + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex

! (p1k2)-vertex with Z/W
	
	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + dreal(4*NC*CF*sqrtS**2*(dconjg(AR)*(BR*(2*muu1*muu4*GL1 + (muu1**2 + muu4**2 - ured)*GR1)*
     &          ((muu2**2 + muu3**2 - ured)*HL(indj,j,tt,gen) + 2*muu2*muu3*HR(indj,j,tt,gen)) + 
     &         BL*(muu1*GL1*(muu3*(muu2**2 + muu4**2 - tred)*HL(indj,j,tt,gen) - 
     &               muu2*(-1 + muu3**2 + muu4**2)*HR(indj,j,tt,gen)) + 
     &            muu4*GR1*(-((-1 + muu1**2 + muu2**2)*muu3*HL(indj,j,tt,gen)) + 
     &               muu2*(muu1**2 + muu3**2 - tred)*HR(indj,j,tt,gen)))) + 
     &      dconjg(AL)*(BL*((muu1**2 + muu4**2 - ured)*GL1 + 2*muu1*muu4*GR1)*
     &          (2*muu2*muu3*HL(indj,j,tt,gen) + (muu2**2 + muu3**2 - ured)*HR(indj,j,tt,gen)) + 
     &         BR*(muu4*GL1*(muu2*(muu1**2 + muu3**2 - tred)*HL(indj,j,tt,gen) - 
     &               (-1 + muu1**2 + muu2**2)*muu3*HR(indj,j,tt,gen)) + 
     &            muu1*GR1*(-(muu2*(-1 + muu3**2 + muu4**2)*HL(indj,j,tt,gen)) + 
     &               muu3*(muu2**2 + muu4**2 - tred)*HR(indj,j,tt,gen))))))/(DenSmV1*(-ured + MSf(j,itypeU,gen)**2/sqrtS**2))
	 
	tempZ12 = tempZ12 + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! T-channel-loop interference Z/W tree-level

! (p2k2)-vertex with Z/W

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + 4*CF*NC*sqrtS**2*(dconjg(AL)*(BR*(2*FL1*muu2*muu4 + FR1*(muu2**2 + muu4**2 - tred))*
     -          ((muu1**2 + muu3**2 - tred)*ELc(indi,j,tt,gen) + 2*muu1*muu3*ER(indi,j,tt,gen)) + 
     -         BL*(FL1*muu2*(muu3*(muu1**2 + muu4**2 - ured)*ELc(indi,j,tt,gen) - muu1*(-1 + muu3**2 + muu4**2)*ER(indi,j,tt,gen)) + 
     -            FR1*muu4*(-((-1 + muu1**2 + muu2**2)*muu3*ELc(indi,j,tt,gen)) + muu1*(muu2**2 + muu3**2 - ured)*ER(indi,j,tt,gen)))) + 
     -      dconjg(AR)*(BL*(2*FR1*muu2*muu4 + FL1*(muu2**2 + muu4**2 - tred))*
     -          (2*muu1*muu3*ELc(indi,j,tt,gen) + (muu1**2 + muu3**2 - tred)*ER(indi,j,tt,gen)) + 
     -         BR*(FL1*muu4*(muu1*(muu2**2 + muu3**2 - ured)*ELc(indi,j,tt,gen) - (-1 + muu1**2 + muu2**2)*muu3*ER(indi,j,tt,gen)) + 
     -            FR1*muu2*(-(muu1*(-1 + muu3**2 + muu4**2)*ELc(indi,j,tt,gen)) + muu3*(muu1**2 + muu4**2 - ured)*ER(indi,j,tt,gen)))))/
     -  (DenSmV1*(tred - MSf(j,itypeT,gen)**2/sqrtS**2))

	tempZ22 = tempZ22 + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! (p1k1)-vertex with Z/W

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex + 4*CF*NC*sqrtS**2*(dconjg(AR)*(BL*(2*EL1*muu1*muu3 + ER1*(muu1**2 + muu3**2 - tred))*
     -          ((muu2**2 + muu4**2 - tred)*FL(indj,j,tt,gen) + 2*muu2*muu4*FR(indj,j,tt,gen)) + 
     -         BR*(muu4*(-(ER1*(-1 + muu1**2 + muu2**2)*muu3) + EL1*muu1*(muu2**2 + muu3**2 - ured))*FL(indj,j,tt,gen) + 
     -            muu2*(-(EL1*muu1*(-1 + muu3**2 + muu4**2)) + ER1*muu3*(muu1**2 + muu4**2 - ured))*FR(indj,j,tt,gen))) + 
     -      dconjg(AL)*(BL*(muu2*(-(ER1*muu1*(-1 + muu3**2 + muu4**2)) + EL1*muu3*(muu1**2 + muu4**2 - ured))*FL(indj,j,tt,gen) + 
     -            muu4*(-(EL1*(-1 + muu1**2 + muu2**2)*muu3) + ER1*muu1*(muu2**2 + muu3**2 - ured))*FR(indj,j,tt,gen)) + 
     -         BR*(2*ER1*muu1*muu3 + EL1*(muu1**2 + muu3**2 - tred))*
     -          (2*muu2*muu4*FL(indj,j,tt,gen) + (muu2**2 + muu4**2 - tred)*FR(indj,j,tt,gen))))/
     -  (DenSmV1*(tred - MSf(j,itypeT,gen)**2/sqrtS**2))

	tempZ11 = tempZ11 + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! U-channel-loop interference Higgs tree-level

	     do i= imin,imax

! (p2k1)-vertex with Higgs

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex - CF*NC*sqrtS**2*(dconjg(CL(i))*(dconjg(DL(i))*
     -            (2*HR1*muu2*muu3*((muu1**2 + muu4**2 - ured)*GL(indi,j,tt,gen) + 2*muu1*muu4*GR(indi,j,tt,gen)) + 
     -              HL1*((1 - muu2**2 - muu3**2 - muu4**2 + 2*muu2**2*muu4**2 - 
     -                    (muu2**2 + muu3**2 + muu4**2 - ured - tred)*(ured - tred) + muu1**2*(-1 + 2*muu3**2 - ured + tred))*
     -                  GL(indi,j,tt,gen) + 2*muu1*muu4*(muu2**2 + muu3**2 - ured)*GR(indi,j,tt,gen))) + 
     -           2*dconjg(DR(i))*(HR1*muu2*(-(muu4*(muu1**2 + muu3**2 - tred)*GL(indi,j,tt,gen)) + 
     -                 muu1*(-1 + muu3**2 + muu4**2)*GR(indi,j,tt,gen)) + 
     -              HL1*muu3*((-1 + muu1**2 + muu2**2)*muu4*GL(indi,j,tt,gen) - muu1*(muu2**2 + muu4**2 - tred)*GR(indi,j,tt,gen)))) + 
     -        dconjg(CR(i))*(2*dconjg(DL(i))*(HR1*muu3*(-(muu1*(muu2**2 + muu4**2 - tred)*GL(indi,j,tt,gen)) + 
     -                 (-1 + muu1**2 + muu2**2)*muu4*GR(indi,j,tt,gen)) + 
     -              HL1*muu2*(muu1*(-1 + muu3**2 + muu4**2)*GL(indi,j,tt,gen) - muu4*(muu1**2 + muu3**2 - tred)*GR(indi,j,tt,gen))) + 
     -           dconjg(DR(i))*(2*HL1*muu2*muu3*(2*muu1*muu4*GL(indi,j,tt,gen) + (muu1**2 + muu4**2 - ured)*GR(indi,j,tt,gen)) + 
     -              HR1*(2*muu1*muu4*(muu2**2 + muu3**2 - ured)*GL(indi,j,tt,gen) + 
     -                 (1 - muu2**2 - muu3**2 - muu4**2 + 2*muu2**2*muu4**2 - (muu2**2 + muu3**2 + muu4**2 - ured - tred)*(ured - tred) + 
     -                    muu1**2*(-1 + 2*muu3**2 - ured + tred))*GR(indi,j,tt,gen)))))/
     -    (DenSmH1(i)*(ured - MSf(j,itypeU,gen)**2/sqrtS**2))
	
	tempH21 = tempH21 + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! (p1k2)-vertex with Higgs

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex - CF*NC*sqrtS**2*(dconjg(CL(i))*(dconjg(DL(i))*
     -            ((2*GR1*muu1*muu4*(muu2**2 + muu3**2 - ured) + 
     -                 GL1*(1 - muu2**2 - muu3**2 - muu4**2 + 2*muu2**2*muu4**2 - 
     -                    (muu2**2 + muu3**2 + muu4**2 - ured - tred)*(ured - tred) + muu1**2*(-1 + 2*muu3**2 - ured + tred)))*
     -               HL(indj,j,tt,gen) + 2*muu2*muu3*(2*GR1*muu1*muu4 + GL1*(muu1**2 + muu4**2 - ured))*HR(indj,j,tt,gen)) + 
     -           2*dconjg(DR(i))*(muu3*(GL1*(-1 + muu1**2 + muu2**2)*muu4 - GR1*muu1*(muu2**2 + muu4**2 - tred))*HL(indj,j,tt,gen) + 
     -              muu2*(GR1*muu1*(-1 + muu3**2 + muu4**2) - GL1*muu4*(muu1**2 + muu3**2 - tred))*HR(indj,j,tt,gen))) + 
     -        dconjg(CR(i))*(2*dconjg(DL(i))*(muu2*(GL1*muu1*(-1 + muu3**2 + muu4**2) - GR1*muu4*(muu1**2 + muu3**2 - tred))*
     -               HL(indj,j,tt,gen) + muu3*(GR1*(-1 + muu1**2 + muu2**2)*muu4 - GL1*muu1*(muu2**2 + muu4**2 - tred))*HR(indj,j,tt,gen)
     -              ) + dconjg(DR(i))*(2*muu2*muu3*(2*GL1*muu1*muu4 + GR1*(muu1**2 + muu4**2 - ured))*HL(indj,j,tt,gen) + 
     -              (2*GL1*muu1*muu4*(muu2**2 + muu3**2 - ured) + 
     -                 GR1*(1 - muu2**2 - muu3**2 - muu4**2 + 2*muu2**2*muu4**2 - 
     -                    (muu2**2 + muu3**2 + muu4**2 - ured - tred)*(ured - tred) + muu1**2*(-1 + 2*muu3**2 - ured + tred)))*
     -               HR(indj,j,tt,gen))))/(DenSmH1(i)*(ured - MSf(j,itypeU,gen)**2/sqrtS**2))

	tempH12 = tempH12 + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! T-channel-loop interference Higgs tree-level

! (p2k2)-vertex with Higgs

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex - CF*NC*sqrtS**2*(dconjg(CL(i))*(2*dconjg(DR(i))*
     -            (FR1*muu2*(-(muu3*(muu1**2 + muu4**2 - ured)*ELc(indi,j,tt,gen)) + muu1*(-1 + muu3**2 + muu4**2)*ER(indi,j,tt,gen)) + 
     -              FL1*muu4*((-1 + muu1**2 + muu2**2)*muu3*ELc(indi,j,tt,gen) - muu1*(muu2**2 + muu3**2 - ured)*ER(indi,j,tt,gen))) + 
     -           dconjg(DL(i))*(2*FR1*muu2*muu4*((muu1**2 + muu3**2 - tred)*ELc(indi,j,tt,gen) + 2*muu1*muu3*ER(indi,j,tt,gen)) + 
     -              FL1*(((-1 + muu3**2 + muu4**2 - ured)*(-1 + ured) + muu2**2*(-1 + 2*muu3**2 + ured - tred) + 
     -                    muu1**2*(-1 + 2*muu4**2 + ured - tred) - (muu3**2 + muu4**2)*tred + tred**2)*ELc(indi,j,tt,gen) + 
     -                 2*muu1*muu3*(muu2**2 + muu4**2 - tred)*ER(indi,j,tt,gen)))) + 
     -        dconjg(CR(i))*(2*dconjg(DL(i))*(FR1*muu4*(-(muu1*(muu2**2 + muu3**2 - ured)*ELc(indi,j,tt,gen)) + 
     -                 (-1 + muu1**2 + muu2**2)*muu3*ER(indi,j,tt,gen)) + 
     -              FL1*muu2*(muu1*(-1 + muu3**2 + muu4**2)*ELc(indi,j,tt,gen) - muu3*(muu1**2 + muu4**2 - ured)*ER(indi,j,tt,gen))) + 
     -           dconjg(DR(i))*(2*FL1*muu2*muu4*(2*muu1*muu3*ELc(indi,j,tt,gen) + (muu1**2 + muu3**2 - tred)*ER(indi,j,tt,gen)) + 
     -              FR1*(2*muu1*muu3*(muu2**2 + muu4**2 - tred)*ELc(indi,j,tt,gen) + 
     -                 ((-1 + muu3**2 + muu4**2 - ured)*(-1 + ured) + muu2**2*(-1 + 2*muu3**2 + ured - tred) + 
     -                    muu1**2*(-1 + 2*muu4**2 + ured - tred) - (muu3**2 + muu4**2)*tred + tred**2)*ER(indi,j,tt,gen)))))/
     -    (DenSmH1(i)*(tred - MSf(j,itypeT,gen)**2/sqrtS**2))

	tempH22 = tempH22 + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! (p1k1)-vertex with Higgs

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex - CF*NC*sqrtS**2*(dconjg(CL(i))*(2*dconjg(DR(i))*
     -            (muu4*(EL1*(-1 + muu1**2 + muu2**2)*muu3 - ER1*muu1*(muu2**2 + muu3**2 - ured))*FL(indj,j,tt,gen) + 
     -              muu2*(ER1*muu1*(-1 + muu3**2 + muu4**2) - EL1*muu3*(muu1**2 + muu4**2 - ured))*FR(indj,j,tt,gen)) + 
     -           dconjg(DL(i))*((2*ER1*muu1*muu3*(muu2**2 + muu4**2 - tred) + 
     -                 EL1*((-1 + muu3**2 + muu4**2 - ured)*(-1 + ured) + muu2**2*(-1 + 2*muu3**2 + ured - tred) + 
     -                    muu1**2*(-1 + 2*muu4**2 + ured - tred) - (muu3**2 + muu4**2)*tred + tred**2))*FL(indj,j,tt,gen) + 
     -              2*muu2*muu4*(2*ER1*muu1*muu3 + EL1*(muu1**2 + muu3**2 - tred))*FR(indj,j,tt,gen))) + 
     -        dconjg(CR(i))*(2*dconjg(DL(i))*(muu2*(EL1*muu1*(-1 + muu3**2 + muu4**2) - ER1*muu3*(muu1**2 + muu4**2 - ured))*
     -               FL(indj,j,tt,gen) + muu4*(ER1*(-1 + muu1**2 + muu2**2)*muu3 - EL1*muu1*(muu2**2 + muu3**2 - ured))*FR(indj,j,tt,gen))
     -             + dconjg(DR(i))*(2*muu2*muu4*(2*EL1*muu1*muu3 + ER1*(muu1**2 + muu3**2 - tred))*FL(indj,j,tt,gen) + 
     -              (2*EL1*muu1*muu3*(muu2**2 + muu4**2 - tred) + 
     -                 ER1*((-1 + muu3**2 + muu4**2 - ured)*(-1 + ured) + muu2**2*(-1 + 2*muu3**2 + ured - tred) + 
     -                    muu1**2*(-1 + 2*muu4**2 + ured - tred) - (muu3**2 + muu4**2)*tred + tred**2))*FR(indj,j,tt,gen))))/
     -    (DenSmH1(i)*(tred - MSf(j,itypeT,gen)**2/sqrtS**2))

	 enddo
	
	tempH11 = tempH11 + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex

! U-channel-loop interference U-channel tree-level

	  do i=1,2
		
! (p2k1)-vertex with U

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex - 2d0*CF*NC*(HL1*((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen))) + 
     -        HR1*(2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen))))*
     -      (-(((muu1**2 + muu4**2 - ured)*dconjg(GL(indi,i,tt,gen)) + 2*muu1*muu4*dconjg(GR(indi,i,tt,gen)))*GL(indi,j,tt,gen)) - 
     -        (2*muu1*muu4*dconjg(GL(indi,i,tt,gen)) + (muu1**2 + muu4**2 - ured)*dconjg(GR(indi,i,tt,gen)))*GR(indi,j,tt,gen))/
     -    ((ured - MSf(i,itypeU,gen)**2/sqrtS**2)*(ured - MSf(j,itypeU,gen)**2/sqrtS**2))
	
	tempu21 = tempu21 + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! (p1k2)-vertex with U

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex - 2d0*CF*NC*(-(GL1*((muu1**2 + muu4**2 - ured)*dconjg(GL(indi,i,tt,gen)) + 2*muu1*muu4*dconjg(GR(indi,i,tt,gen)))) - 
     -        GR1*(2*muu1*muu4*dconjg(GL(indi,i,tt,gen)) + (muu1**2 + muu4**2 - ured)*dconjg(GR(indi,i,tt,gen))))*
     -      (((muu2**2 + muu3**2 - ured)*dconjg(HL(indj,i,tt,gen)) + 2*muu2*muu3*dconjg(HR(indj,i,tt,gen)))*HL(indj,j,tt,gen) + 
     -        (2*muu2*muu3*dconjg(HL(indj,i,tt,gen)) + (muu2**2 + muu3**2 - ured)*dconjg(HR(indj,i,tt,gen)))*HR(indj,j,tt,gen))/
     -    ((ured - MSf(i,itypeU,gen)**2/sqrtS**2)*(ured - MSf(j,itypeU,gen)**2/sqrtS**2))

	tempu12 = tempu12 + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! U-channel-loop interference T-channel tree-level

! (p2k1)-vertex with T

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex - CF*NC*(HL1*((2*muu2*muu4*dconjg(FR(indj,i,tt,gen))*
     -             ((muu1**2 + muu3**2 - tred)*dconjg(ELc(indi,i,tt,gen)) + 2*muu1*muu3*dconjg(ER(indi,i,tt,gen))) + 
     -            dconjg(FL(indj,i,tt,gen))*((-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - 
     -                  (muu2**2 + muu3**2 + muu4**2 - ured)*ured + muu1**2*(1 + 2*muu2**2 - ured - tred) - 
     -                  (muu2**2 + muu3**2 + muu4**2)*tred + tred**2)*dconjg(ELc(indi,i,tt,gen)) + 
     -               2*muu1*muu3*(muu2**2 + muu4**2 - tred)*dconjg(ER(indi,i,tt,gen))))*GL(indi,j,tt,gen) + 
     -         2*(muu4*dconjg(FL(indj,i,tt,gen))*(muu1*(muu2**2 + muu3**2 - ured)*dconjg(ELc(indi,i,tt,gen)) - 
     -               (-1 + muu1**2 + muu2**2)*muu3*dconjg(ER(indi,i,tt,gen))) + 
     -            muu2*dconjg(FR(indj,i,tt,gen))*(-(muu1*(-1 + muu3**2 + muu4**2)*dconjg(ELc(indi,i,tt,gen))) + 
     -               muu3*(muu1**2 + muu4**2 - ured)*dconjg(ER(indi,i,tt,gen))))*GR(indi,j,tt,gen)) + 
     -      HR1*(2*(muu2*dconjg(FL(indj,i,tt,gen))*(muu3*(muu1**2 + muu4**2 - ured)*dconjg(ELc(indi,i,tt,gen)) - 
     -               muu1*(-1 + muu3**2 + muu4**2)*dconjg(ER(indi,i,tt,gen))) + 
     -            muu4*dconjg(FR(indj,i,tt,gen))*(-((-1 + muu1**2 + muu2**2)*muu3*dconjg(ELc(indi,i,tt,gen))) + 
     -               muu1*(muu2**2 + muu3**2 - ured)*dconjg(ER(indi,i,tt,gen))))*GL(indi,j,tt,gen) + 
     -         (2*muu2*muu4*dconjg(FL(indj,i,tt,gen))*(2*muu1*muu3*dconjg(ELc(indi,i,tt,gen)) + 
     -               (muu1**2 + muu3**2 - tred)*dconjg(ER(indi,i,tt,gen))) + 
     -            dconjg(FR(indj,i,tt,gen))*(2*muu1*muu3*(muu2**2 + muu4**2 - tred)*dconjg(ELc(indi,i,tt,gen)) + 
     -               (-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - (muu2**2 + muu3**2 + muu4**2 - ured)*ured + 
     -                  muu1**2*(1 + 2*muu2**2 - ured - tred) - (muu2**2 + muu3**2 + muu4**2)*tred + tred**2)*dconjg(ER(indi,i,tt,gen))))*
     -          GR(indi,j,tt,gen)))/((tred - MSf(i,itypeT,gen)**2/sqrtS**2)*(ured - MSf(j,itypeU,gen)**2/sqrtS**2))
	
	tempt21 = tempt21 + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! (p1k2)-vertex with T

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex - CF*NC*((2*GR1*(muu4*dconjg(FL(indj,i,tt,gen))*(muu1*(muu2**2 + muu3**2 - ured)*dconjg(ELc(indi,i,tt,gen)) - 
     -               (-1 + muu1**2 + muu2**2)*muu3*dconjg(ER(indi,i,tt,gen))) + 
     -            muu2*dconjg(FR(indj,i,tt,gen))*(-(muu1*(-1 + muu3**2 + muu4**2)*dconjg(ELc(indi,i,tt,gen))) + 
     -               muu3*(muu1**2 + muu4**2 - ured)*dconjg(ER(indi,i,tt,gen)))) + 
     -         GL1*(2*muu2*muu4*dconjg(FR(indj,i,tt,gen))*
     -             ((muu1**2 + muu3**2 - tred)*dconjg(ELc(indi,i,tt,gen)) + 2*muu1*muu3*dconjg(ER(indi,i,tt,gen))) + 
     -            dconjg(FL(indj,i,tt,gen))*((-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - 
     -                  (muu2**2 + muu3**2 + muu4**2 - ured)*ured + muu1**2*(1 + 2*muu2**2 - ured - tred) - 
     -                  (muu2**2 + muu3**2 + muu4**2)*tred + tred**2)*dconjg(ELc(indi,i,tt,gen)) + 
     -               2*muu1*muu3*(muu2**2 + muu4**2 - tred)*dconjg(ER(indi,i,tt,gen)))))*HL(indj,j,tt,gen) + 
     -      (2*GL1*(muu2*dconjg(FL(indj,i,tt,gen))*(muu3*(muu1**2 + muu4**2 - ured)*dconjg(ELc(indi,i,tt,gen)) - 
     -               muu1*(-1 + muu3**2 + muu4**2)*dconjg(ER(indi,i,tt,gen))) + 
     -            muu4*dconjg(FR(indj,i,tt,gen))*(-((-1 + muu1**2 + muu2**2)*muu3*dconjg(ELc(indi,i,tt,gen))) + 
     -               muu1*(muu2**2 + muu3**2 - ured)*dconjg(ER(indi,i,tt,gen)))) + 
     -         GR1*(2*muu2*muu4*dconjg(FL(indj,i,tt,gen))*
     -             (2*muu1*muu3*dconjg(ELc(indi,i,tt,gen)) + (muu1**2 + muu3**2 - tred)*dconjg(ER(indi,i,tt,gen))) + 
     -            dconjg(FR(indj,i,tt,gen))*(2*muu1*muu3*(muu2**2 + muu4**2 - tred)*dconjg(ELc(indi,i,tt,gen)) + 
     -               (-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - (muu2**2 + muu3**2 + muu4**2 - ured)*ured + 
     -                  muu1**2*(1 + 2*muu2**2 - ured - tred) - (muu2**2 + muu3**2 + muu4**2)*tred + tred**2)*dconjg(ER(indi,i,tt,gen)))))*
     -       HR(indj,j,tt,gen))/((tred - MSf(i,itypeT,gen)**2/sqrtS**2)*(ured - MSf(j,itypeU,gen)**2/sqrtS**2))
	
	tempt12 = tempt12 + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! T-channel-loop interference U-channel tree-level

! (p2k2)-vertex with U

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex - CF*NC*(dconjg(HL(indj,i,tt,gen))*(2*dconjg(GR(indi,i,tt,gen))*
     -          (FL1*muu4*(muu1*(muu2**2 + muu3**2 - ured)*ELc(indi,j,tt,gen) - (-1 + muu1**2 + muu2**2)*muu3*ER(indi,j,tt,gen)) + 
     -            FR1*muu2*(-(muu1*(-1 + muu3**2 + muu4**2)*ELc(indi,j,tt,gen)) + muu3*(muu1**2 + muu4**2 - ured)*ER(indi,j,tt,gen))) + 
     -         dconjg(GL(indi,i,tt,gen))*(2*FR1*muu2*muu4*
     -             ((muu1**2 + muu3**2 - tred)*ELc(indi,j,tt,gen) + 2*muu1*muu3*ER(indi,j,tt,gen)) + 
     -            FL1*((-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - (muu2**2 + muu3**2 + muu4**2 - ured)*ured + 
     -                  muu1**2*(1 + 2*muu2**2 - ured - tred) - (muu2**2 + muu3**2 + muu4**2)*tred + tred**2)*ELc(indi,j,tt,gen) + 
     -               2*muu1*muu3*(muu2**2 + muu4**2 - tred)*ER(indi,j,tt,gen)))) + 
     -      dconjg(HR(indj,i,tt,gen))*(2*dconjg(GL(indi,i,tt,gen))*
     -          (FL1*muu2*(muu3*(muu1**2 + muu4**2 - ured)*ELc(indi,j,tt,gen) - muu1*(-1 + muu3**2 + muu4**2)*ER(indi,j,tt,gen)) + 
     -            FR1*muu4*(-((-1 + muu1**2 + muu2**2)*muu3*ELc(indi,j,tt,gen)) + muu1*(muu2**2 + muu3**2 - ured)*ER(indi,j,tt,gen))) + 
     -         dconjg(GR(indi,i,tt,gen))*(2*FL1*muu2*muu4*
     -             (2*muu1*muu3*ELc(indi,j,tt,gen) + (muu1**2 + muu3**2 - tred)*ER(indi,j,tt,gen)) + 
     -            FR1*(2*muu1*muu3*(muu2**2 + muu4**2 - tred)*ELc(indi,j,tt,gen) + 
     -               (-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - (muu2**2 + muu3**2 + muu4**2 - ured)*ured + 
     -                  muu1**2*(1 + 2*muu2**2 - ured - tred) - (muu2**2 + muu3**2 + muu4**2)*tred + tred**2)*ER(indi,j,tt,gen)))))/
     -  ((ured - MSf(i,itypeU,gen)**2/sqrtS**2)*(tred - MSf(j,itypeT,gen)**2/sqrtS**2))
	
	tempu22 = tempu22 + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! (p1k1)-vertex with U

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex - CF*NC*(dconjg(HL(indj,i,tt,gen))*(2*dconjg(GR(indi,i,tt,gen))*
     -          (muu4*(-(ER1*(-1 + muu1**2 + muu2**2)*muu3) + EL1*muu1*(muu2**2 + muu3**2 - ured))*FL(indj,j,tt,gen) + 
     -            muu2*(-(EL1*muu1*(-1 + muu3**2 + muu4**2)) + ER1*muu3*(muu1**2 + muu4**2 - ured))*FR(indj,j,tt,gen)) + 
     -         dconjg(GL(indi,i,tt,gen))*((2*ER1*muu1*muu3*(muu2**2 + muu4**2 - tred) + 
     -               EL1*(-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - (muu2**2 + muu3**2 + muu4**2 - ured)*ured + 
     -                  muu1**2*(1 + 2*muu2**2 - ured - tred) - (muu2**2 + muu3**2 + muu4**2)*tred + tred**2))*FL(indj,j,tt,gen) + 
     -            2*muu2*muu4*(2*ER1*muu1*muu3 + EL1*(muu1**2 + muu3**2 - tred))*FR(indj,j,tt,gen))) + 
     -      dconjg(HR(indj,i,tt,gen))*(2*dconjg(GL(indi,i,tt,gen))*
     -          (muu2*(-(ER1*muu1*(-1 + muu3**2 + muu4**2)) + EL1*muu3*(muu1**2 + muu4**2 - ured))*FL(indj,j,tt,gen) + 
     -            muu4*(-(EL1*(-1 + muu1**2 + muu2**2)*muu3) + ER1*muu1*(muu2**2 + muu3**2 - ured))*FR(indj,j,tt,gen)) + 
     -         dconjg(GR(indi,i,tt,gen))*(2*muu2*muu4*(2*EL1*muu1*muu3 + ER1*(muu1**2 + muu3**2 - tred))*FL(indj,j,tt,gen) + 
     -            (2*EL1*muu1*muu3*(muu2**2 + muu4**2 - tred) + 
     -               ER1*(-1 + muu2**2 + muu3**2 + muu4**2 + 2*muu3**2*muu4**2 - (muu2**2 + muu3**2 + muu4**2 - ured)*ured + 
     -                  muu1**2*(1 + 2*muu2**2 - ured - tred) - (muu2**2 + muu3**2 + muu4**2)*tred + tred**2))*FR(indj,j,tt,gen))))/
     -  ((ured - MSf(i,itypeU,gen)**2/sqrtS**2)*(tred - MSf(j,itypeT,gen)**2/sqrtS**2))
	
	tempu11 = tempu11 + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! T-channel-loop interference T-channel tree-level

! (p2k2)-vertex with T

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex - 2d0*CF*NC*(FL1*((muu2**2 + muu4**2 - tred)*dconjg(FL(indj,i,tt,gen)) + 2*muu2*muu4*dconjg(FR(indj,i,tt,gen))) + 
     -        FR1*(2*muu2*muu4*dconjg(FL(indj,i,tt,gen)) + (muu2**2 + muu4**2 - tred)*dconjg(FR(indj,i,tt,gen))))*
     -      (-(((muu1**2 + muu3**2 - tred)*dconjg(ELc(indi,i,tt,gen)) + 2*muu1*muu3*dconjg(ER(indi,i,tt,gen)))*ELc(indi,j,tt,gen)) - 
     -        (2*muu1*muu3*dconjg(ELc(indi,i,tt,gen)) + (muu1**2 + muu3**2 - tred)*dconjg(ER(indi,i,tt,gen)))*ER(indi,j,tt,gen))/
     -    ((tred - MSf(i,itypeT,gen)**2/sqrtS**2)*(tred - MSf(j,itypeT,gen)**2/sqrtS**2))
	
	tempt22 = tempt22 + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
! (p1k1)-vertex with T

	ChiChi2QQ_M2Vertex = ChiChi2QQ_M2Vertex - 2*CF*NC*(-(EL1*((muu1**2 + muu3**2 - tred)*dconjg(ELc(indi,i,tt,gen)) + 2*muu1*muu3*dconjg(ER(indi,i,tt,gen)))) - 
     -      ER1*(2*muu1*muu3*dconjg(ELc(indi,i,tt,gen)) + (muu1**2 + muu3**2 - tred)*dconjg(ER(indi,i,tt,gen))))*
     -    (((muu2**2 + muu4**2 - tred)*dconjg(FL(indj,i,tt,gen)) + 2*muu2*muu4*dconjg(FR(indj,i,tt,gen)))*FL(indj,j,tt,gen) + 
     -      (2*muu2*muu4*dconjg(FL(indj,i,tt,gen)) + (muu2**2 + muu4**2 - tred)*dconjg(FR(indj,i,tt,gen)))*FR(indj,j,tt,gen))/
     -  ((tred - MSf(i,itypeT,gen)**2/sqrtS**2)*(tred - MSf(j,itypeT,gen)**2/sqrtS**2))
	
	tempt11 = tempt11 + ChiChi2QQ_M2Vertex - tempCounter
	tempCounter = ChiChi2QQ_M2Vertex
	
	  enddo
	enddo
	
	
cc	write(*,*)'ChiChi TU11 Vertex mit Z = ', tempZ11
cc	write(*,*)'ChiChi TU11 Vertex mit H = ', tempH11
cc	write(*,*)'ChiChi TU11 Vertex mit t = ', tempt11
cc	write(*,*)'ChiChi TU11 Vertex mit u = ', tempu11
cc    
cc	write(*,*)'ChiChi TU12 Vertex mit Z = ', tempZ12
cc	write(*,*)'ChiChi TU12 Vertex mit H = ', tempH12
cc	write(*,*)'ChiChi TU12 Vertex mit t = ', tempt12
cc	write(*,*)'ChiChi TU12 Vertex mit u = ', tempu12
cc    
cc	write(*,*)'ChiChi TU21 Vertex mit Z = ', tempZ21
cc	write(*,*)'ChiChi TU21 Vertex mit H = ', tempH21
cc	write(*,*)'ChiChi TU21 Vertex mit t = ', tempt21
cc	write(*,*)'ChiChi TU21 Vertex mit u = ', tempu21
cc    
cc	write(*,*)'ChiChi TU22 Vertex mit Z = ', tempZ22
cc	write(*,*)'ChiChi TU22 Vertex mit H = ', tempH22
cc	write(*,*)'ChiChi TU22 Vertex mit t = ', tempt22
cc	write(*,*)'ChiChi TU22 Vertex mit u = ', tempu22
	
cc	write(*,*)'---------------------------'
cc	write(*,*)'ER1 = ', ER1
cc	write(*,*)'EL1 = ', EL1
cc	write(*,*)'FR1 = ', FR1
cc	write(*,*)'FL1 = ', FL1
cc	write(*,*)'GR1 = ', GR1
cc	write(*,*)'GL1 = ', GL1
cc	write(*,*)'HR1 = ', HR1
cc	write(*,*)'HL1 = ', HL1
    
cc	write(*,*)'ER = ', ER(indi,1,tt,gen)
cc	write(*,*)'ELc = ', ELc(indi,1,tt,gen)
cc	write(*,*)'FR = ', FR(indj,1,tt,gen)
cc	write(*,*)'FL = ', FL(indj,1,tt,gen)
cc	write(*,*)'GR = ', GR(indi,1,tt,gen)
cc	write(*,*)'GL = ', GL(indi,1,tt,gen)
cc	write(*,*)'HR = ', HR(indj,1,tt,gen)
cc	write(*,*)'HL = ', HL(indj,1,tt,gen)
cc	write(*,*)'---------------------------'
	
	end
