!==================================================================================================
!============= Tests the loopfunction in the case of a vanishing Gram determinant =================
!=======================       Author: Patrick Steppeler		  =================================
!==================================================================================================

	subroutine DD_Looptest
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_Model.h"
#include "../util/DMNLO_LoopIntegrals.h"
#include "DD_GenCouplings.h"
	
!----------------------- Declare variables --------------------------------------------------------
	integer tt, gen, i, j, k, l
	double precision p1p1pat, p2p2pat
	double precision p10p10pat, p20p20pat, p30p30pat, p21p21pat, p32p32pat, p31p31pat
	double precision rotp10p10pat, rotp20p20pat, rotp30p30pat, rotp21p21pat, rotp32p32pat, rotp31p31pat
	double precision t_Mandeltest, s_Mandeltest, s_min, Mandelvar
	double precision m0m0, m1m1, m2m2, m3m3, pcm, costh, kappa
	complex*16 results(18), resultstemp(18)
    
	open (unit=20,file="results.txt",action="write",status="replace")
    

!======================= ACTUAL CALCULATIONS ======================================================
	
!----------------------- Renormalization scales for loop integrals --------------------------------
	UVdiv = 0d1
	IRdiv = 0d1
 	EpsPole = 0d0
	
	l = 1

!----------------------- Squark Gluon Vertex ------------------------------------------------------
C 	tt = 3
C 	gen = 3
C 	i = 1
C
C 	m0m0 = Mf(tt,gen)**2
C 	m1m1 = MSf(i,tt,gen)**2
C 	m2m2 = 0d0
C
C 	p1p1pat = Mneu(1)**2
C 	p2p2pat = Mf(tt,gen)**2
C 	s_min = (Mneu(1) + Mf(tt,gen))**2
C 	s_Mandeltest = s_min + 10d0

	
!----------------------- Higgs Gluino Vertex ------------------------------------------------------
C 	tt = 3
C 	gen = 3
C 	i = 1
C 	j = 2
C
C 	m0m0 = MGl**2
C 	m1m1 = MSf(i,tt,gen)**2
C 	m2m2 = MSf(j,tt,gen)**2
C
C 	p1p1pat = Mf(tt,gen)**2
C 	p2p2pat = Mf(tt,gen)**2
C 	t_Mandeltest = - 1d0

!----------------------- Loop over pcm ------------------------------------------------------------
	pcm = 10.0d0
	costh = 0.5d0

	do while (pcm.ge.0d0)

!----------------------- Gluon Box B1 ------------------------------------------------------------
C 		tt = 3
C 		gen = 1
C 		i = 1
C
C 		s_min = (Mneu(1) + Mf(tt,gen))**2
C
C 		m0m0 = 0d0
C 		m1m1 = MSf(i,tt,gen)**2
C 		m2m2 = Mf(tt,gen)**2
C 		m3m3 = Mf(tt,gen)**2
C
C ! General kinematics
C 		if(pcm.ne.0d0) then
C 			s_Mandeltest = MNeu(1)**2 + 2.0d0*pcm**2 + Mf(tt,gen)**2 + 2.0d0*
C      &				dsqrt(MNeu(1)**2 + pcm**2)*dsqrt(Mf(tt,gen)**2 +pcm**2)
C 			t_Mandeltest = 2.0d0*MNeu(1)**2 - 0.5d0*(s_Mandeltest + MNeu(1)**2 - Mf(tt,gen)**2)**2/s_Mandeltest +
C      &				0.5d0*kappa(s_Mandeltest, MNeu(1)**2, Mf(tt,gen)**2)**2/s_Mandeltest*costh
C
C ! Nonrelativistic approximation for pcm = 0, this choice sets the whole code pcm independent
C 		else
C 			s_Mandeltest = (MNeu(1) + Mf(tt,gen))**2
C 			t_Mandeltest = 0d0
C 		endif
C
C 		p10p10pat = s_Mandeltest
C 		p20p20pat = Mf(tt,gen)**2
C 		p30p30pat = Mf(tt,gen)**2
C 		p21p21pat = MNeu(1)**2
C 		p32p32pat = t_Mandeltest
C 		p31p31pat = MNeu(1)**2
			
C !----------------------- Gluino Box B2 ------------------------------------------------------------
		tt = 3
		gen = 3
		i = 1
		j = 2

		s_min = (Mneu(1) + Mf(tt,gen))**2

		m0m0 = MGl**2
		m1m1 = MSf(i,tt,gen)**2
		m2m2 = MSf(j,tt,gen)**2
		m3m3 = Mf(tt,gen)**2
		
! General kinematics
		if(pcm.ne.0d0) then
			s_Mandeltest = MNeu(1)**2 + 2.0d0*pcm**2 + Mf(tt,gen)**2 + 2.0d0*
     &				dsqrt(MNeu(1)**2 + pcm**2)*dsqrt(Mf(tt,gen)**2 +pcm**2)
			t_Mandeltest = 2.0d0*MNeu(1)**2 - 0.5d0*(s_Mandeltest + MNeu(1)**2 - Mf(tt,gen)**2)**2/s_Mandeltest +
     &				0.5d0*kappa(s_Mandeltest, MNeu(1)**2, Mf(tt,gen)**2)**2/s_Mandeltest*costh

! Nonrelativistic approximation for pcm = 0, this choice sets the whole code pcm independent
		else
			s_Mandeltest = (MNeu(1) + Mf(tt,gen))**2
			t_Mandeltest = 0d0
		endif

		p10p10pat = Mf(tt,gen)**2
		p20p20pat = Mf(tt,gen)**2
		p30p30pat = s_Mandeltest
		p21p21pat = t_Mandeltest
		p32p32pat = Mneu(1)**2
		p31p31pat = Mneu(1)**2

!----------------------- Loop over Mandelstam variable --------------------------------------------
C 	do while (t_Mandeltest.le.0d0)
C 		Mandelvar = t_Mandeltest
C
C 	do while (s_Mandeltest.ge.s_min)
C  		Mandelvar = s_Mandeltest

!----------------------- Alternative & classic C calculation --------------------------------------
		call myCget(p1p1pat, Mandelvar, p2p2pat, m0m0, m1m1, m2m2)
		call AltCget(p1p1pat, Mandelvar, p2p2pat, m0m0, m1m1, m2m2, results)
		
!----------------------- Alternative & classic D Calculation --------------------------------------

! Check if rotation 1 <-> 3 is necessary
		if ((p30p30pat.eq.p20p20pat).and.(p32p32pat.eq.0d0).and.(p31p31pat.eq.p21p21pat).and.(p10p10pat.ne.p30p30pat)) then

			write(*,*)'Yes, rotate!'

			rotp10p10pat = p30p30pat
			rotp20p20pat = p20p20pat
			rotp30p30pat = p10p10pat

			rotp21p21pat = p32p32pat
			rotp32p32pat = p21p21pat
			rotp31p31pat = p31p31pat

			call AltDget(rotp10p10pat, rotp21p21pat, rotp32p32pat, rotp30p30pat, rotp20p20pat, rotp31p31pat, m0m0, m3m3, m2m2, m1m1, results)

			do k = 1, 18
				resultstemp(k) = results(k)
			enddo

			results(3) = resultstemp(7)
			results(4) = resultstemp(7)
			results(7) = resultstemp(3)
			results(8) = resultstemp(4)

! No rotation necessary
		else

			write(*,*)'No rotation.'
			call AltDget(p10p10pat, p21p21pat, p32p32pat, p30p30pat, p20p20pat, p31p31pat, m0m0, m1m1, m2m2, m3m3, results)

		endif

		call myDget(p10p10pat, p21p21pat, p32p32pat, p30p30pat, p20p20pat, p31p31pat, m0m0, m1m1, m2m2, m3m3)
		
!----------------------- Output for C calculation -------------------------------------------------
C 		write(*,*)' --- Loopergebnisse: --- '
C c 		write(*,*)' s = ', s_Mandeltest
C c 		write(*,*)' s_crit = ', s_Mandeltest-s_min
C 		write(*,*)' t = ', t_Mandeltest
C 		write(*,*)' C0 = ', xc0
C 		write(*,*)' C0_alt = ', results(1), results(2)
C C 		write(*,*)' C1 = ', xc1
C C 		write(*,*)' C1_alt = ', results(3), results(4)
C C 		write(*,*)' C2 = ', xc2
C C 		write(*,*)' C2_alt = ', results(5), results(6)
C C 		write(*,*)' C1 + C2 = ', xc1 + xc2
C C 		write(*,*)' C1_alt + C2_alt = ', results(3) + results(5), results(4) + results(6)
C C 		write(*,*)' Ca_alt = ', results(7), results(8)
C 		write(*,*)' C00 = ', xc00
C  		write(*,*)' C00_alt = ', results(9), results(10)
C  		write(*,*)' C11 = ', xc11
C c 		write(*,*)' C11_alt = ', results(13), results(14)
C  		write(*,*)' C12 = ', xc12
C c 		write(*,*)' C12_alt = ', results(15), results(16)
C  		write(*,*)' C22 = ', xc22
C c		write(*,*)' C22_alt = ', results(17), results(18)
C  		write(*,*)' C11 + 2*C12 + C22 = ', xc11 + 2d0*xc12 + xc22
C  		write(*,*)' C11_alt +2*C12_alt + C22_alt = ', results(11), results(12)
C 		write(*,*)' ----------------------- '
C 		write(*,*)''
C
C        write (20,*),l, t_Mandeltest,dreal(xc1), dreal(xc2),dreal(results(7)),dreal(xc00),dreal(xc11),dreal(xc12),dreal(xc22),dreal(results(9)),dreal(results(11))
C 		l = l + 1
		
!----------------------- Output for D calculation -------------------------------------------------
C 		write(*,*)' --- Loopergebnisse: --- '
C 		write(*,*)' pcm = ',pcm
C 		write(*,*)' s = ', s_Mandeltest
C 		write(*,*)' s_crit = ', s_Mandeltest - s_min
C 		write(*,*)' t = ', t_Mandeltest
C 		write(*,*)' D0 = ', xd0
C 		write(*,*)' D0_alt = ', results(1), results(2)
C 		write(*,*)' D1 = ', xd1
C 		write(*,*)' D1_alt = ', results(3), results(4)
C 		write(*,*)' D2 = ', xd2
C 		write(*,*)' D2_alt = ', results(5), results(6)
C 		write(*,*)' D3 = ', xd3
C 		write(*,*)' D3_alt = ', results(7), results(8)
C 		write(*,*)' D00 = ', xd00
C 		write(*,*)' D00_alt = ', results(9), results(10)
C 		write(*,*)' D11 = ', xd11
C 		write(*,*)' D11_alt = ', results(11), results(12)
C 		write(*,*)' D13 + D12 = ', xd13 + xd12
C 		write(*,*)' D3a_alt = ', results(13), results(14)
C 		write(*,*)' D22 + 2D23 + D33 = ', xd22 + 2d0*xd23 + xd33
C 		write(*,*)' Daa_alt = ', results(15), results(16)
C 		write(*,*)' ----------------------- '
C 		write(*,*)''
		
        write (20,*) l, pcm,dreal(xd1), dreal(xd2),dreal(xd3),dreal(xd00),dreal(xd11),dreal(xd12),dreal(xd13),dreal(xd22),dreal(xd23),dreal(xd33),dreal(results(3)),dreal(results(5)),dreal(results(7)),dreal(results(9)),dreal(results(11)),dreal(results(13)),dreal(results(15))
 		l = l + 1
		
!----------------------- Loop over Mandelstam variable --------------------------------------------
C 		t_Mandeltest = t_Mandeltest + 0.0001d0
C 		if(abs(t_Mandeltest).lt.1d-7) then
C 			t_Mandeltest = 0d0
C 		endif

C 		s_Mandeltest = s_Mandeltest - 0.1d0
C 		if(abs(s_Mandeltest-s_min).lt.1d-6) then
C 			s_Mandeltest = s_min
C 		endif

C	enddo

!----------------------- Loop over pcm --------------------------------------------
		pcm = pcm - 1d0
		if(abs(pcm).lt.1d-3) then
			pcm = 0d0
		endif

	enddo
	

!======================= END OF PROGRAM ===========================================================
	close (unit=20)
	end
	
!==================================================================================================
!============= Calculates D functions in the case of a vanishing Gram determinant =================
!=======================       Author: Patrick Steppeler		  =================================
!==================================================================================================

	subroutine AltDget(p10p10pat, p21p21pat, p32p32pat, p30p30pat, p20p20pat, p31p31pat, m0m0, m1m1, m2m2, m3m3, results)
	implicit none
	
!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "../util/DMNLO_LoopIntegrals.h"

!----------------------- Declare variables --------------------------------------------------------
	double precision p10p10pat, p20p20pat, p30p30pat, p21p21pat, p32p32pat, p31p31pat
	double precision p1p2pat, p1p3pat, p2p3pat
	double precision m0m0, m1m1, m2m2, m3m3
	double precision f1, f2, f3
	double precision x1(2), x2(2), Z11(2), Z12(2), Z21(2), Z22(2), DetZ(2), Y(2)
	complex*16 W1(2), W2(2)
	complex*16 D0_alt(2), D1_alt(2), D2_alt(2), D3_alt(2), D00_alt(2), Da_alt(2)
	complex*16 D11_alt(2), D12_alt(2), D13_alt(2), D22_alt(2), D23_alt(2), D33_alt(2)
	complex*16 D3a_alt(2), Daa_alt(2)
	complex*16 Da_test(2), D3_test(2)
	complex*16 C0012, C0013, C0023, C0123
	complex*16 C1012, C1013, C1023, C1123
	complex*16 C2012, C2013, C2023, C2123
	complex*16 C00012, C00013, C00023, C00123
	complex*16 C11012, C11013, C11023, C11123
	complex*16 C12012, C12013, C12023, C12123
	complex*16 C22012, C22013, C22023, C22123
	complex*16 results(18)
	integer i, k
			
!----------------------- Preparations -------------------------------------------------------------
	do i = 1,2
		D0_alt(i) = 0d0
		D1_alt(i) = 0d0
		D2_alt(i) = 0d0
		D3_alt(i) = 0d0
		Da_alt(i) = 0d0
		D00_alt(i) = 0d0
		D11_alt(i) = 0d0
		D12_alt(i) = 0d0
		D13_alt(i) = 0d0
		D22_alt(i) = 0d0
		D23_alt(i) = 0d0
		D33_alt(i) = 0d0
		Daa_alt(i) = 0d0
		D3a_alt(i) = 0d0
	enddo
	
	f1 = p10p10pat - m1m1 + m0m0
	f2 = p20p20pat - m2m2 + m0m0
	f3 = p30p30pat - m3m3 + m0m0
	
	p1p2pat = - 0.5d0*(p21p21pat - p10p10pat - p20p20pat)
	p1p3pat = - 0.5d0*(p31p31pat - p10p10pat - p30p30pat)
	p2p3pat = - 0.5d0*(p32p32pat - p20p20pat - p30p30pat)
	
	call myCget(p10p10pat,p21p21pat, p20p20pat, m0m0, m1m1, m2m2)
c	write(*,*)'Argument set 1: ', p10p10pat,p21p21pat, p20p20pat, m0m0, m1m1, m2m2
	C0012 = xc0
	C1012 = xc1
	C2012 = xc2
	C00012 = xc00
	C11012 = xc11
	C12012 = xc12
	C22012 = xc22
c	write(*,*)'C0, C1, C2 = ', xc0, xc1, xc2
c	write(*,*)'C00, C11, C12, C22 = ', xc00, xc11, xc12, xc22
	
	call myCget(p10p10pat,p31p31pat, p30p30pat, m0m0, m1m1, m3m3)
c	write(*,*)'Argument set 2: ', p10p10pat,p31p31pat, p30p30pat, m0m0, m1m1, m3m3
	C0013 = xc0
	C1013 = xc1
	C2013 = xc2
	C00013 = xc00
	C11013 = xc11
	C12013 = xc12
	C22013 = xc22
c	write(*,*)'C0, C1, C2 = ', xc0, xc1, xc2
c	write(*,*)'C00, C11, C12, C22 = ', xc00, xc11, xc12, xc22
	
	call myCget(p20p20pat,p32p32pat, p30p30pat, m0m0, m2m2, m3m3)
c	write(*,*)'Argument set 3: ', p20p20pat,p32p32pat, p30p30pat, m0m0, m2m2, m3m3
	C0023 = xc0
	C1023 = xc1
	C2023 = xc2
	C00023 = xc00
	C11023 = xc11
	C12023 = xc12
	C22023 = xc22
c	write(*,*)'C0, C1, C2 = ', xc0, xc1, xc2
c	write(*,*)'C00, C11, C12, C22 = ', xc00, xc11, xc12, xc22
	
	call myCget(p21p21pat,p32p32pat, p31p31pat, m1m1, m2m2, m3m3)
c	write(*,*)'Argument set 4: ', p21p21pat,p32p32pat, p31p31pat, m1m1, m2m2, m3m3
	C0123 = xc0
	C1123 = xc1
	C2123 = xc2
	C00123 = xc00
	C11123 = xc11
	C12123 = xc12
	C22123 = xc22
c	write(*,*)'C0, C1, C2 = ', xc0, xc1, xc2
c	write(*,*)'C00, C11, C12, C22 = ', xc00, xc11, xc12, xc22
	
	call myDget(p10p10pat, p21p21pat, p32p32pat, p30p30pat, p20p20pat, p31p31pat, m0m0, m1m1, m2m2, m3m3)
	
!----------------------- Alternative D0 determination ---------------------------------------------

! Box B2 or B4 with zero momentum transfer D(p, p, p3, m0**2, m1**2, m2**2, m3**2) = (p**2, 0, (p3-p)**2, p3**2, p**2, (p3-p)**2, m0**2, m1**2, m2**2, m3**2)
	if ((p10p10pat.eq.p20p20pat).and.(p21p21pat.eq.0d0).and.(p31p31pat.eq.p32p32pat).and.(p10p10pat.ne.p30p30pat)) then
				
		if(p10p10pat.eq.0d0) then
					
			D0_alt(1) = xd0
			D0_alt(2) = xd0
			
		else
			
			D0_alt(1) = 1d0/(p1p3pat*f3 - p30p30pat*f1)*(p1p3pat*(C0012 - C0123) - p30p30pat*(C0023 - C0123))
			D0_alt(2) = 1d0/(p1p3pat*f3 - p30p30pat*f2)*(p1p3pat*(C0012 - C0123) - p30p30pat*(C0013 - C0123))
			
		endif
		
	else
		
		write(*,*)'Error in alternative D tensor determination - Case not found'

	endif
	
!----------------------- Alternative D1, D2 & CD2 determination -----------------------------------

! Box B2 or B4 with zero momentum transfer D(p, p, p3, m0**2, m1**2, m2**2, m3**2) = (p**2, 0, (p3-p)**2, p3**2, p**2, (p3-p)**2, m0**2, m1**2, m2**2, m3**2)
	if ((p10p10pat.eq.p20p20pat).and.(p21p21pat.eq.0d0).and.(p31p31pat.eq.p32p32pat).and.(p10p10pat.ne.p30p30pat)) then
				
		if(p10p10pat.eq.0d0) then
					
			do i = 1,2
				D1_alt(i) = 0d0
				D2_alt(i) = 0d0
				D3_alt(i) = 1d0/(2d0*p30p30pat)*(C0012 - C0123 - f3*xd0)
				D11_alt(i) = 0d0
				D12_alt(i) = 0d0
				D13_alt(i) = 0d0
				D22_alt(i) = 0d0
				D23_alt(i) = 0d0
			enddo
			
			D00_alt(1) = 1d0/3d0*(C0123 + m0m0*xd0 - 0.5d0*(C2012 - C2123 - f3*D3_alt(1)))
			D00_alt(2) = 1d0/3d0*(C0023 + m1m1*xd0 - 0.5d0*(C2012 - C2123 - f3*D3_alt(2)))
			D33_alt(1) = 1d0/p30p30pat*(2d0/3d0*(C2012 -C2123 - f3*D3_alt(1)) - 1d0/3d0*(C0123 + m0m0*xd0))
			D33_alt(2) = 1d0/p30p30pat*(2d0/3d0*(C2012 -C2123 - f3*D3_alt(2)) - 1d0/3d0*(C0023 + m1m1*xd0))
			
		else
			
			x1(1) = p10p10pat
			x1(2) = - p1p3pat
			x2(1) = - p1p3pat
			x2(2) = p30p30pat
			
			do k = 1,2
			
				Y(k) = x1(k)*f3 + x2(k)*f1
		
				Z11(k) = Y(k) + x1(k)/2d0*f3
				Z12(k) = x1(k)/2d0*f1
				Z21(k) = x2(k)/2d0*f3
				Z22(k) = Y(k) + x2(k)/2d0*f1
				
				DetZ(k) = Z11(k)*Z22(k) - Z12(k)*Z21(k)
				
				W1(k) = x1(k)*(-C2123 - m0m0*D0_alt(k) + 0.5*(C1023 - C0123)) + x2(k)*(C2023 - C2123)
				W2(k) = x1(k)*(C1012 + C2012 + C0123 + C2123) + x2(k)*(3d0/2d0*C1023 + C2123 - m0m0*D0_alt(k) + 0.5d0*C0123)
				
				D3_alt(k) = 1d0/DetZ(k)*(Z22(k)*W1(k) - Z12(k)*W2(k))
				Da_alt(k) = 1d0/DetZ(k)*(-Z21(k)*W1(k) + Z11(k)*W2(k))
				
				D1_alt(k) = Da_alt(k)/2d0
				D2_alt(k) = Da_alt(k)/2d0
				D00_alt(k) = 0.5d0*(m0m0*D0_alt(k) + 0.5d0*(f3*D3_alt(k) + f1*Da_alt(k) - C1023 + C0123))
				
			enddo
			
! Mathematica results for D3 and Da				
C 			D3_test(1) = -((-((-2d0*f3*p1p3pat + 3d0*f1*p30p30pat)*(-(C0123*p10p10pat)/4d0 + (C1023*p10p10pat)/4d0 - (C2123*p10p10pat)/2d0 - (C2023*p1p3pat)/2d0 +
C      -     	(C2123*p1p3pat)/2d0 - (m0m0*p10p10pat*xd0)/2d0)) -
C      -        (f1*p10p10pat*(2d0*C0123*p1p3pat + 2d0*C1012*p1p3pat + 2d0*C2012*p1p3pat + 2d0*C2123*p1p3pat - C0123*p30p30pat - 3d0*C1023*p30p30pat - 2d0*C2123*p30p30pat +
C      -        2d0*m0m0*p30p30pat*xd0))/4d0)/((-3d0*f3**2*p10p10pat*p1p3pat)/2d0 + f1*f3*p1p3pat**2 + 2d0*f1*f3*p10p10pat*p30p30pat - (3d0*f1**2*p1p3pat*p30p30pat)/2d0))
C
C 			Da_test(1) = -((-3d0*C0123*f3*p10p10pat*p1p3pat - 3d0*C1012*f3*p10p10pat*p1p3pat - 3d0*C2012*f3*p10p10pat*p1p3pat - 3d0*C2123*f3*p10p10pat*p1p3pat + 2d0*C0123*f1*p1p3pat**2 +
C      -        2d0*C1012*f1*p1p3pat**2 + 2d0*C2012*f1*p1p3pat**2 + 2d0*C2123*f1*p1p3pat**2 + 2d0*C0123*f3*p10p10pat*p30p30pat + 4d0*C1023*f3*p10p10pat*p30p30pat +
C      -        4d0*C2123*f3*p10p10pat*p30p30pat - C0123*f1*p1p3pat*p30p30pat - 3d0*C1023*f1*p1p3pat*p30p30pat - 2d0*C2123*f1*p1p3pat*p30p30pat + C2023*f3*p1p3pat*p30p30pat -
C      -        C2123*f3*p1p3pat*p30p30pat - 2d0*f3*m0m0*p10p10pat*p30p30pat*xd0 + 2d0*f1*m0m0*p1p3pat*p30p30pat*xd0)/
C      -		(3d0*f3**2*p10p10pat*p1p3pat - 2d0*f1*f3*p1p3pat**2 - 4d0*f1*f3*p10p10pat*p30p30pat + 3d0*f1**2*p1p3pat*p30p30pat))
C
C 			D3_test(2) = -((-(f1*p1p3pat*(-2d0*C0123*p10p10pat - 2d0*C1012*p10p10pat - 2d0*C2012*p10p10pat - 2d0*C2123*p10p10pat + C0123*p1p3pat + 3d0*C1023*p1p3pat +
C      -        2d0*C2123*p1p3pat - 2d0*m0m0*p1p3pat*xd0)) - (2d0*f3*p10p10pat - 3d0*f1*p1p3pat)*
C      -        (-(C0123*p1p3pat) + C1023*p1p3pat - 2d0*C2123*p1p3pat - 2d0*C2023*p30p30pat + 2d0*C2123*p30p30pat - 2d0*m0m0*p1p3pat*xd0))/
C      -		(f1*f3*p1p3pat**2 - (2d0*f3*p10p10pat - 3d0*f1*p1p3pat)*(-3d0*f3*p1p3pat + 2d0*f1*p30p30pat)))
C
C 			Da_test(2) = -((-3d0*C0123*f3*p10p10pat*p1p3pat - 3d0*C1012*f3*p10p10pat*p1p3pat - 3d0*C2012*f3*p10p10pat*p1p3pat - 3d0*C2123*f3*p10p10pat*p1p3pat + 2d0*C0123*f3*p1p3pat**2 +
C      -        4d0*C1023*f3*p1p3pat**2 + 4d0*C2123*f3*p1p3pat**2 + 2d0*C0123*f1*p10p10pat*p30p30pat + 2d0*C1012*f1*p10p10pat*p30p30pat + 2d0*C2012*f1*p10p10pat*p30p30pat +
C      -        2d0*C2123*f1*p10p10pat*p30p30pat - C0123*f1*p1p3pat*p30p30pat - 3d0*C1023*f1*p1p3pat*p30p30pat - 2d0*C2123*f1*p1p3pat*p30p30pat + C2023*f3*p1p3pat*p30p30pat -
C      -        C2123*f3*p1p3pat*p30p30pat - 2d0*f3*m0m0*p1p3pat**2*xd0 + 2d0*f1*m0m0*p1p3pat*p30p30pat*xd0)/
C      -		(3d0*f3**2*p10p10pat*p1p3pat - 4d0*f1*f3*p1p3pat**2 - 2d0*f1*f3*p10p10pat*p30p30pat + 3d0*f1**2*p1p3pat*p30p30pat))
			
! Mathematica results for D33, D3a and Daa						
			D33_alt(1) = -((-(((C12023*p10p10pat)/3. - (C2123*p10p10pat)/3. - (C22123*p10p10pat)/2. - (2*D3_alt(1)*m0m0*p10p10pat)/3. - (C22023*p1p3pat)/2. + 
     -          (C22123*p1p3pat)/2.)*((f3**2*p1p3pat**2)/3. - (5*f1*f3*p1p3pat*p30p30pat)/6. + (5*f1**2*p30p30pat**2)/9.)) - 
     -         	(f1*p10p10pat*(-(((f3*p1p3pat)/2. - (5*f1*p30p30pat)/6.)*
     -          ((C0123*p1p3pat)/6. + (C11023*p1p3pat)/6. + (2*C2123*p1p3pat)/3. + (C22123*p1p3pat)/2. - (Da_alt(1)*m0m0*p1p3pat)/3. - (2*C12023*p30p30pat)/3. - 
     -          (C2123*p30p30pat)/3. - (C22123*p30p30pat)/2. + (D3_alt(1)*m0m0*p30p30pat)/3.)) - 
     -          (f1*p1p3pat*((C0123*p1p3pat)/2. - (C11012*p1p3pat)/2. - C12012*p1p3pat + C2123*p1p3pat - (C22012*p1p3pat)/2. + (C22123*p1p3pat)/2. - 
     -          (C0123*p30p30pat)/6. + (5*C11023*p30p30pat)/6. - (2*C2123*p30p30pat)/3. - (C22123*p30p30pat)/2. - (2*Da_alt(1)*m0m0*p30p30pat)/3.))/6.))/3.)/
     -       	(-(((-5*f3*p10p10pat)/6. + (f1*p1p3pat)/2.)*((f3**2*p1p3pat**2)/3. - (5*f1*f3*p1p3pat*p30p30pat)/6. + (5*f1**2*p30p30pat**2)/9.)) - 
     -         	(f1*p10p10pat*(-(f3**2*p1p3pat*p30p30pat)/12. + (5*f1*f3*p30p30pat**2)/36.))/3.))
	 
     		D3a_alt(1) = -(5*C0123*f1*f3*p10p10pat*p1p3pat**2 - 5*C11012*f1*f3*p10p10pat*p1p3pat**2 - 10*C12012*f1*f3*p10p10pat*p1p3pat**2 + 
     -        	10*C2123*f1*f3*p10p10pat*p1p3pat**2 - 5*C22012*f1*f3*p10p10pat*p1p3pat**2 + 5*C22123*f1*f3*p10p10pat*p1p3pat**2 + 5*C0123*f3**2*p10p10pat*p1p3pat**2 + 
     -        	5*C11023*f3**2*p10p10pat*p1p3pat**2 + 20*C2123*f3**2*p10p10pat*p1p3pat**2 + 15*C22123*f3**2*p10p10pat*p1p3pat**2 - 10*Da_alt(1)*f3**2*m0m0*p10p10pat*p1p3pat**2 - 
     -        	3*C0123*f1**2*p1p3pat**3 + 3*C11012*f1**2*p1p3pat**3 + 6*C12012*f1**2*p1p3pat**3 - 6*C2123*f1**2*p1p3pat**3 + 3*C22012*f1**2*p1p3pat**3 - 
     -        	3*C22123*f1**2*p1p3pat**3 - 3*C0123*f1*f3*p1p3pat**3 - 3*C11023*f1*f3*p1p3pat**3 - 12*C2123*f1*f3*p1p3pat**3 - 9*C22123*f1*f3*p1p3pat**3 + 
     -        	6*Da_alt(1)*f1*f3*m0m0*p1p3pat**3 - 10*C0123*f1*f3*p10p10pat*p1p3pat*p30p30pat - 40*C2123*f1*f3*p10p10pat*p1p3pat*p30p30pat - 
     -        	30*C22123*f1*f3*p10p10pat*p1p3pat*p30p30pat - 18*C12023*f3**2*p10p10pat*p1p3pat*p30p30pat - 12*C2123*f3**2*p10p10pat*p1p3pat*p30p30pat - 
     -        	18*C22123*f3**2*p10p10pat*p1p3pat*p30p30pat + 10*Da_alt(1)*f1*f3*m0m0*p10p10pat*p1p3pat*p30p30pat + 6*D3_alt(1)*f3**2*m0m0*p10p10pat*p1p3pat*p30p30pat + 
     -        	6*C0123*f1**2*p1p3pat**2*p30p30pat + 24*C2123*f1**2*p1p3pat**2*p30p30pat + 18*C22123*f1**2*p1p3pat**2*p30p30pat + 12*C12023*f1*f3*p1p3pat**2*p30p30pat + 
     -        	6*C2123*f1*f3*p1p3pat**2*p30p30pat + 9*C22123*f1*f3*p1p3pat**2*p30p30pat - 3*C22023*f3**2*p1p3pat**2*p30p30pat + 3*C22123*f3**2*p1p3pat**2*p30p30pat - 
     -        	6*Da_alt(1)*f1**2*m0m0*p1p3pat**2*p30p30pat - 6*D3_alt(1)*f1*f3*m0m0*p1p3pat**2*p30p30pat + 30*C12023*f1*f3*p10p10pat*p30p30pat**2 + 
     -        	20*C2123*f1*f3*p10p10pat*p30p30pat**2 + 30*C22123*f1*f3*p10p10pat*p30p30pat**2 - 10*D3_alt(1)*f1*f3*m0m0*p10p10pat*p30p30pat**2 - 
     -        	20*C12023*f1**2*p1p3pat*p30p30pat**2 - 10*C2123*f1**2*p1p3pat*p30p30pat**2 - 15*C22123*f1**2*p1p3pat*p30p30pat**2 + 5*C22023*f1*f3*p1p3pat*p30p30pat**2 - 
     -       	5*C22123*f1*f3*p1p3pat*p30p30pat**2 + 10*D3_alt(1)*f1**2*m0m0*p1p3pat*p30p30pat**2)/
     -     		(2.*(-10*f3**3*p10p10pat*p1p3pat**2 + 6*f1*f3**2*p1p3pat**3 + 24*f1*f3**2*p10p10pat*p1p3pat*p30p30pat - 15*f1**2*f3*p1p3pat**2*p30p30pat - 
     -         	15*f1**2*f3*p10p10pat*p30p30pat**2 + 10*f1**3*p1p3pat*p30p30pat**2))
	 
 			Daa_alt(1) = -((-10*C0123*f3**2*p10p10pat*p1p3pat**2 + 10*C11012*f3**2*p10p10pat*p1p3pat**2 + 20*C12012*f3**2*p10p10pat*p1p3pat**2 - 
     -         	20*C2123*f3**2*p10p10pat*p1p3pat**2 + 10*C22012*f3**2*p10p10pat*p1p3pat**2 - 10*C22123*f3**2*p10p10pat*p1p3pat**2 + 6*C0123*f1*f3*p1p3pat**3 - 
     -         	6*C11012*f1*f3*p1p3pat**3 - 12*C12012*f1*f3*p1p3pat**3 + 12*C2123*f1*f3*p1p3pat**3 - 6*C22012*f1*f3*p1p3pat**3 + 6*C22123*f1*f3*p1p3pat**3 + 
     -         	9*C0123*f1*f3*p10p10pat*p1p3pat*p30p30pat - 9*C11012*f1*f3*p10p10pat*p1p3pat*p30p30pat - 18*C12012*f1*f3*p10p10pat*p1p3pat*p30p30pat + 
     -         	18*C2123*f1*f3*p10p10pat*p1p3pat*p30p30pat - 9*C22012*f1*f3*p10p10pat*p1p3pat*p30p30pat + 9*C22123*f1*f3*p10p10pat*p1p3pat*p30p30pat + 
     -         	5*C0123*f3**2*p10p10pat*p1p3pat*p30p30pat - 15*C11023*f3**2*p10p10pat*p1p3pat*p30p30pat + 20*C2123*f3**2*p10p10pat*p1p3pat*p30p30pat + 
     -         	15*C22123*f3**2*p10p10pat*p1p3pat*p30p30pat + 10*Da_alt(1)*f3**2*m0m0*p10p10pat*p1p3pat*p30p30pat - 6*C0123*f1**2*p1p3pat**2*p30p30pat + 
     -         	6*C11012*f1**2*p1p3pat**2*p30p30pat + 12*C12012*f1**2*p1p3pat**2*p30p30pat - 12*C2123*f1**2*p1p3pat**2*p30p30pat + 6*C22012*f1**2*p1p3pat**2*p30p30pat - 
     -         	6*C22123*f1**2*p1p3pat**2*p30p30pat - 3*C0123*f1*f3*p1p3pat**2*p30p30pat + 9*C11023*f1*f3*p1p3pat**2*p30p30pat - 12*C2123*f1*f3*p1p3pat**2*p30p30pat - 
     -         	9*C22123*f1*f3*p1p3pat**2*p30p30pat - 6*Da_alt(1)*f1*f3*m0m0*p1p3pat**2*p30p30pat - 3*C0123*f1*f3*p10p10pat*p30p30pat**2 + 
     -         	15*C11023*f1*f3*p10p10pat*p30p30pat**2 - 12*C2123*f1*f3*p10p10pat*p30p30pat**2 - 9*C22123*f1*f3*p10p10pat*p30p30pat**2 - 
     -         	6*C12023*f3**2*p10p10pat*p30p30pat**2 - 4*C2123*f3**2*p10p10pat*p30p30pat**2 - 6*C22123*f3**2*p10p10pat*p30p30pat**2 - 
     -         	12*Da_alt(1)*f1*f3*m0m0*p10p10pat*p30p30pat**2 + 2*D3_alt(1)*f3**2*m0m0*p10p10pat*p30p30pat**2 + 2*C0123*f1**2*p1p3pat*p30p30pat**2 - 
     -         	10*C11023*f1**2*p1p3pat*p30p30pat**2 + 8*C2123*f1**2*p1p3pat*p30p30pat**2 + 6*C22123*f1**2*p1p3pat*p30p30pat**2 + 4*C12023*f1*f3*p1p3pat*p30p30pat**2 + 
     -         	2*C2123*f1*f3*p1p3pat*p30p30pat**2 + 3*C22123*f1*f3*p1p3pat*p30p30pat**2 - C22023*f3**2*p1p3pat*p30p30pat**2 + C22123*f3**2*p1p3pat*p30p30pat**2 + 
     -         	8*Da_alt(1)*f1**2*m0m0*p1p3pat*p30p30pat**2 - 2*D3_alt(1)*f1*f3*m0m0*p1p3pat*p30p30pat**2)/
     -       	(-10*f3**3*p10p10pat*p1p3pat**2 + 6*f1*f3**2*p1p3pat**3 + 24*f1*f3**2*p10p10pat*p1p3pat*p30p30pat - 15*f1**2*f3*p1p3pat**2*p30p30pat - 
     -         	15*f1**2*f3*p10p10pat*p30p30pat**2 + 10*f1**3*p1p3pat*p30p30pat**2))

	 		D33_alt(2) = -((-2*f1*p1p3pat*(-((3*f3*p10p10pat - 5*f1*p1p3pat)*
     -          (-(C0123*p10p10pat) - C11023*p10p10pat - 4*C2123*p10p10pat - 3*C22123*p10p10pat + 2*Da_alt(2)*m0m0*p10p10pat + 4*C12023*p1p3pat + 2*C2123*p1p3pat + 
     -          3*C22123*p1p3pat - 2*D3_alt(2)*m0m0*p1p3pat)) + f1*p10p10pat*
     -          (3*C0123*p10p10pat - 3*C11012*p10p10pat - 6*C12012*p10p10pat + 6*C2123*p10p10pat - 3*C22012*p10p10pat + 3*C22123*p10p10pat - C0123*p1p3pat + 
     -          5*C11023*p1p3pat - 4*C2123*p1p3pat - 3*C22123*p1p3pat - 4*Da_alt(2)*m0m0*p1p3pat)) - 
     -         	(-12*f3**2*p10p10pat**2 + 30*f1*f3*p10p10pat*p1p3pat - 20*f1**2*p1p3pat**2)*
     -          (2*C12023*p1p3pat - 2*C2123*p1p3pat - 3*C22123*p1p3pat - 4*D3_alt(2)*m0m0*p1p3pat - 3*C22023*p30p30pat + 3*C22123*p30p30pat))/
     -       	(-2*f1*f3*p1p3pat**2*(3*f3*p10p10pat - 5*f1*p1p3pat) - 
     -         	(-12*f3**2*p10p10pat**2 + 30*f1*f3*p10p10pat*p1p3pat - 20*f1**2*p1p3pat**2)*(-5*f3*p1p3pat + 3*f1*p30p30pat)))
	 
     		D3a_alt(2) = -(5*C0123*f1*f3*p10p10pat**2*p1p3pat - 5*C11012*f1*f3*p10p10pat**2*p1p3pat - 10*C12012*f1*f3*p10p10pat**2*p1p3pat + 
     -        	10*C2123*f1*f3*p10p10pat**2*p1p3pat - 5*C22012*f1*f3*p10p10pat**2*p1p3pat + 5*C22123*f1*f3*p10p10pat**2*p1p3pat + 5*C0123*f3**2*p10p10pat**2*p1p3pat + 
     -        	5*C11023*f3**2*p10p10pat**2*p1p3pat + 20*C2123*f3**2*p10p10pat**2*p1p3pat + 15*C22123*f3**2*p10p10pat**2*p1p3pat - 10*Da_alt(2)*f3**2*m0m0*p10p10pat**2*p1p3pat - 
     -        	10*C0123*f1*f3*p10p10pat*p1p3pat**2 - 40*C2123*f1*f3*p10p10pat*p1p3pat**2 - 30*C22123*f1*f3*p10p10pat*p1p3pat**2 - 18*C12023*f3**2*p10p10pat*p1p3pat**2 - 
     -        	12*C2123*f3**2*p10p10pat*p1p3pat**2 - 18*C22123*f3**2*p10p10pat*p1p3pat**2 + 10*Da_alt(2)*f1*f3*m0m0*p10p10pat*p1p3pat**2 + 
     -        	6*D3_alt(2)*f3**2*m0m0*p10p10pat*p1p3pat**2 + 30*C12023*f1*f3*p1p3pat**3 + 20*C2123*f1*f3*p1p3pat**3 + 30*C22123*f1*f3*p1p3pat**3 - 10*D3_alt(2)*f1*f3*m0m0*p1p3pat**3 - 
     -        	3*C0123*f1**2*p10p10pat**2*p30p30pat + 3*C11012*f1**2*p10p10pat**2*p30p30pat + 6*C12012*f1**2*p10p10pat**2*p30p30pat - 
     -        	6*C2123*f1**2*p10p10pat**2*p30p30pat + 3*C22012*f1**2*p10p10pat**2*p30p30pat - 3*C22123*f1**2*p10p10pat**2*p30p30pat - 
     -        	3*C0123*f1*f3*p10p10pat**2*p30p30pat - 3*C11023*f1*f3*p10p10pat**2*p30p30pat - 12*C2123*f1*f3*p10p10pat**2*p30p30pat - 
     -        	9*C22123*f1*f3*p10p10pat**2*p30p30pat + 6*Da_alt(2)*f1*f3*m0m0*p10p10pat**2*p30p30pat + 6*C0123*f1**2*p10p10pat*p1p3pat*p30p30pat + 
     -        	24*C2123*f1**2*p10p10pat*p1p3pat*p30p30pat + 18*C22123*f1**2*p10p10pat*p1p3pat*p30p30pat + 12*C12023*f1*f3*p10p10pat*p1p3pat*p30p30pat + 
     -        	6*C2123*f1*f3*p10p10pat*p1p3pat*p30p30pat + 9*C22123*f1*f3*p10p10pat*p1p3pat*p30p30pat - 3*C22023*f3**2*p10p10pat*p1p3pat*p30p30pat + 
     -        	3*C22123*f3**2*p10p10pat*p1p3pat*p30p30pat - 6*Da_alt(2)*f1**2*m0m0*p10p10pat*p1p3pat*p30p30pat - 6*D3_alt(2)*f1*f3*m0m0*p10p10pat*p1p3pat*p30p30pat - 
     -        	20*C12023*f1**2*p1p3pat**2*p30p30pat - 10*C2123*f1**2*p1p3pat**2*p30p30pat - 15*C22123*f1**2*p1p3pat**2*p30p30pat + 5*C22023*f1*f3*p1p3pat**2*p30p30pat - 
     -        	5*C22123*f1*f3*p1p3pat**2*p30p30pat + 10*D3_alt(2)*f1**2*m0m0*p1p3pat**2*p30p30pat)/
     -     		(2.*(-10*f3**3*p10p10pat**2*p1p3pat + 24*f1*f3**2*p10p10pat*p1p3pat**2 - 15*f1**2*f3*p1p3pat**3 + 6*f1*f3**2*p10p10pat**2*p30p30pat - 
     -         	15*f1**2*f3*p10p10pat*p1p3pat*p30p30pat + 10*f1**3*p1p3pat**2*p30p30pat))
	 
     		Daa_alt(2) = -((-10*C0123*f3**2*p10p10pat**2*p1p3pat + 10*C11012*f3**2*p10p10pat**2*p1p3pat + 20*C12012*f3**2*p10p10pat**2*p1p3pat - 
     -         	20*C2123*f3**2*p10p10pat**2*p1p3pat + 10*C22012*f3**2*p10p10pat**2*p1p3pat - 10*C22123*f3**2*p10p10pat**2*p1p3pat + 9*C0123*f1*f3*p10p10pat*p1p3pat**2 - 
     -         	9*C11012*f1*f3*p10p10pat*p1p3pat**2 - 18*C12012*f1*f3*p10p10pat*p1p3pat**2 + 18*C2123*f1*f3*p10p10pat*p1p3pat**2 - 9*C22012*f1*f3*p10p10pat*p1p3pat**2 + 
     -         	9*C22123*f1*f3*p10p10pat*p1p3pat**2 + 5*C0123*f3**2*p10p10pat*p1p3pat**2 - 15*C11023*f3**2*p10p10pat*p1p3pat**2 + 20*C2123*f3**2*p10p10pat*p1p3pat**2 + 
     -         	15*C22123*f3**2*p10p10pat*p1p3pat**2 + 10*Da_alt(2)*f3**2*m0m0*p10p10pat*p1p3pat**2 - 3*C0123*f1*f3*p1p3pat**3 + 15*C11023*f1*f3*p1p3pat**3 - 
     -         	12*C2123*f1*f3*p1p3pat**3 - 9*C22123*f1*f3*p1p3pat**3 - 6*C12023*f3**2*p1p3pat**3 - 4*C2123*f3**2*p1p3pat**3 - 6*C22123*f3**2*p1p3pat**3 - 
     -         	12*Da_alt(2)*f1*f3*m0m0*p1p3pat**3 + 2*D3_alt(2)*f3**2*m0m0*p1p3pat**3 + 6*C0123*f1*f3*p10p10pat**2*p30p30pat - 6*C11012*f1*f3*p10p10pat**2*p30p30pat - 
     -         	12*C12012*f1*f3*p10p10pat**2*p30p30pat + 12*C2123*f1*f3*p10p10pat**2*p30p30pat - 6*C22012*f1*f3*p10p10pat**2*p30p30pat + 
     -         	6*C22123*f1*f3*p10p10pat**2*p30p30pat - 6*C0123*f1**2*p10p10pat*p1p3pat*p30p30pat + 6*C11012*f1**2*p10p10pat*p1p3pat*p30p30pat + 
     -         	12*C12012*f1**2*p10p10pat*p1p3pat*p30p30pat - 12*C2123*f1**2*p10p10pat*p1p3pat*p30p30pat + 6*C22012*f1**2*p10p10pat*p1p3pat*p30p30pat - 
     -         	6*C22123*f1**2*p10p10pat*p1p3pat*p30p30pat - 3*C0123*f1*f3*p10p10pat*p1p3pat*p30p30pat + 9*C11023*f1*f3*p10p10pat*p1p3pat*p30p30pat - 
     -         	12*C2123*f1*f3*p10p10pat*p1p3pat*p30p30pat - 9*C22123*f1*f3*p10p10pat*p1p3pat*p30p30pat - 6*Da_alt(2)*f1*f3*m0m0*p10p10pat*p1p3pat*p30p30pat + 
     -         	2*C0123*f1**2*p1p3pat**2*p30p30pat - 10*C11023*f1**2*p1p3pat**2*p30p30pat + 8*C2123*f1**2*p1p3pat**2*p30p30pat + 6*C22123*f1**2*p1p3pat**2*p30p30pat + 
     -         	4*C12023*f1*f3*p1p3pat**2*p30p30pat + 2*C2123*f1*f3*p1p3pat**2*p30p30pat + 3*C22123*f1*f3*p1p3pat**2*p30p30pat - C22023*f3**2*p1p3pat**2*p30p30pat + 
     -         	C22123*f3**2*p1p3pat**2*p30p30pat + 8*Da_alt(2)*f1**2*m0m0*p1p3pat**2*p30p30pat - 2*D3_alt(2)*f1*f3*m0m0*p1p3pat**2*p30p30pat)/
     -       	(-10*f3**3*p10p10pat**2*p1p3pat + 24*f1*f3**2*p10p10pat*p1p3pat**2 - 15*f1**2*f3*p1p3pat**3 + 6*f1*f3**2*p10p10pat**2*p30p30pat - 
     -         	15*f1**2*f3*p10p10pat*p1p3pat*p30p30pat + 10*f1**3*p1p3pat**2*p30p30pat))
			
			do k = 1,2
				
				D11_alt(k) = Daa_alt(k)/4d0
				D12_alt(k) = Daa_alt(k)/4d0
				D22_alt(k) = Daa_alt(k)/4d0
				D13_alt(k) = D3a_alt(k)/2d0
				D23_alt(k) = D3a_alt(k)/2d0
				
			enddo
			
		endif
		
	else
	
		write(*,*)'Error in alternative D tensor determination - Case not found'
			
	endif
		
	results(1) = D0_alt(1)
	results(2) = D0_alt(2)
	results(3) = D1_alt(1)
	results(4) = D1_alt(2)
	results(5) = D2_alt(1)
	results(6) = D2_alt(2)
	results(7) = D3_alt(1)
	results(8) = D3_alt(2)
	results(9) = D00_alt(1)
	results(10) = D00_alt(2)
	results(11) = D33_alt(1)
	results(12) = D33_alt(2)
	results(13) = D3a_alt(1)
	results(14) = D3a_alt(2)
	results(15) = Daa_alt(1)
	results(16) = Daa_alt(2)
		
!======================= END OF PROGRAM ===========================================================

	end

	
!==================================================================================================
!============= Calculates C functions in the case of a vanishing Gram determinant =================
!=======================       Author: Patrick Steppeler		  =================================
!==================================================================================================

	subroutine AltCget(p1p1pat, Mandelvar, p2p2pat, m0m0, m1m1, m2m2, results)
	implicit none
	
!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "../util/DMNLO_LoopIntegrals.h"

!----------------------- Declare variables --------------------------------------------------------
	integer i, k
	double precision p1p1pat, p2p2pat, p1p2pat, m1m1, m2m2, m0m0, f1, f2
	double precision Mandelvar, s_min
	double precision x1(2), x2(2), Z11(2), Z12(2), Z21(2), Z22(2), DetZ(2), Y(2)
	complex *16 R7(2), R8(2)
	complex *16 C0_alt(2), B012, B112, C1_alt(2), C2_alt(2), C00_alt(2), Caa_alt(2)
	complex *16 C11_alt(2), C12_alt(2), C22_alt(2)
	complex *16 Ca_alt(2), B001, B002, B101, B102, results(18)
	complex *16 B0001, B0002, B0012, B1101, B1102, B1112
		
!----------------------- Preparations -------------------------------------------------------------
	do i = 1,2
		C0_alt(i) = 0d0
		C1_alt(i) = 0d0 
		C2_alt(i) = 0d0
		Ca_alt(i) = 0d0
		Caa_alt(i) = 0d0
		C00_alt(i) = 0d0
		C11_alt(i) = 0d0
		C12_alt(i) = 0d0
		C22_alt(i) = 0d0
	enddo
	
	call myBget(p1p1pat,m0m0,m1m1)
	B001 = xb0
	B101 = xb1
	B0001 = xb00
	B1101 = xb11
	
	call myBget(p2p2pat,m0m0,m2m2)
	B002 = xb0
	B102 = xb1
	B0002 = xb00
	B1102 = xb11
	
	call myBget(Mandelvar,m1m1,m2m2)
	B012 = xb0
	B112 = xb1
	B0012 = xb00
	B1112 = xb11
	
	f1 = p1p1pat - m1m1 + m0m0
	f2 = p2p2pat - m2m2 + m0m0
	p1p2pat = -0.5d0*(Mandelvar - p1p1pat - p2p2pat)
	s_min = (dsqrt(p1p1pat) + dsqrt(p2p2pat))**2
	
	
!----------------------- Alternative C0 determination ---------------------------------------------
	! For t-channel
	if((dabs(p1p1pat - p2p2pat).lt.1d-4).and.((dabs(p1p2pat - p2p2pat).lt.1d-4).or.((p2p2pat.ne.0d0).and.(dabs(p1p2pat/p2p2pat - 1d0).lt.1d-4)))) then
		if(dabs(f1 - f2).lt.1d-4) then
			write(*,*)'Can not determine C0 in an alternative fashion, redirect to classic C0'
			call myCget(p1p1pat, Mandelvar, p2p2pat, m0m0, m1m1, m2m2)
			C0_alt(1) = xc0
			C0_alt(2) = xc0
		else
			C0_alt(1) = (B002 - B001)/(f1 - f2)
			C0_alt(2) = (B001 - B002)/(f2 - f1)
		endif
	else
		C0_alt(1) = (p2p2pat*B002 + (p1p2pat-p2p2pat)*B012 - p1p2pat*B001)/(p2p2pat*f1 - p1p2pat*f2)
		C0_alt(2) = (p1p1pat*B001 + (p1p2pat-p1p1pat)*B012 - p1p2pat*B002)/(p1p1pat*f2 - p1p2pat*f1)
	endif
	
	! For s-channel
	if(p1p1pat.ne.p2p2pat) then
		if(p2p2pat.ne.0d0) then
			C0_alt(1) = (p2p2pat*B002 + (p1p2pat-p2p2pat)*B012 - p1p2pat*B001)/(p2p2pat*f1 - p1p2pat*f2)
			C0_alt(2) = (p1p1pat*B001 + (p1p2pat-p1p1pat)*B012 - p1p2pat*B002)/(p1p1pat*f2 - p1p2pat*f1)
		elseif(f2.ne.0d0) then
			C0_alt(1) = (B012 - B001)/( - f2)
			C0_alt(2) = (p1p1pat*B001 + (p1p2pat-p1p1pat)*B012 - p1p2pat*B002)/(p1p1pat*f2 - p1p2pat*f1)
		else
			write(*,*)'Can not determine C0 in an alternative fashion, redirect to classic C0'
			call myCget(p1p1pat, Mandelvar, p2p2pat, m0m0, m1m1, m2m2)
			C0_alt(1) = xc0
			C0_alt(2) = xc0
		endif
	endif

	
!----------------------- Select special case ------------------------------------------------------

! t-channel with both momenta equal and zero momentum transfer C(p**2, 0, p**2, m0**2, m1**2, m2**2)	
		if ((Mandelvar.eq.0d0).and.(p1p1pat.eq.p2p2pat)) then
				
			if(p1p1pat.eq.0d0) then
				Ca_alt(1) = 0d0
				Ca_alt(2) = 0d0
				Caa_alt(1) = 0d0
				Caa_alt(2) = 0d0
				C00_alt(1) = 1/4d0*(B012 + m0m0*xc0 + 0.5d0) 
				C00_alt(2) = 1/4d0*(B012 + m0m0*xc0 + 0.5d0)
			else
				Ca_alt(1) = 1/(2d0*p1p1pat)*(B001-B012-f2*xc0)
				Ca_alt(2) = 1/(2d0*p1p1pat)*(B002-B012-f1*xc0)
				Caa_alt(1) = 1d0/(3d0*p1p1pat)*(B012 - m0m0*xc0 + 2d0*B102 - 2d0*f1*Ca_alt(1) - 0.5d0) 
				Caa_alt(2) = 1d0/(3d0*p1p1pat)*(B012 - m0m0*xc0 + 2d0*B101 - 2d0*f2*Ca_alt(2) - 0.5d0)
				C00_alt(1) = 1/3d0*m0m0*xc0 + 1/6d0*(B012 - B102 + f1*Ca_alt(1) + 1d0)
				C00_alt(2) = 1/3d0*m0m0*xc0 + 1/6d0*(B012 - B101 + f2*Ca_alt(2) + 1d0)
			endif
		
! s-channel or u-channel with one momentum zero and zero momentum transfer C(p**2, p**2, 0, m0**2, m1**2, m2**2)
		elseif((Mandelvar.eq.p1p1pat).and.(p2p2pat.eq.0d0)) then
		
			if(p1p1pat.eq.0d0) then
				Ca_alt(1) = 0d0
				Ca_alt(2) = 0d0
				Caa_alt(1) = 0d0
				Caa_alt(2) = 0d0
				C00_alt(1) = 1/4d0*(B012 + m0m0*xc0 + 0.5d0) 
				C00_alt(2) = 1/4d0*(B012 + m0m0*xc0 + 0.5d0)
			else
				Ca_alt(1) = 1/(2d0*p1p1pat)*(B002-B012-f1*xc0)
				Ca_alt(2) = 1/(2d0*p1p1pat)*(B002-B012-f1*xc0)
				Caa_alt(1) = 1d0/(3d0*p1p1pat)*(B012 - m0m0*xc0 + 2d0*B112 - 2d0*f1*Ca_alt(1) - 0.5d0) 
				Caa_alt(2) = 1d0/(3d0*p1p1pat)*(B012 - m0m0*xc0 + 2d0*B112 - 2d0*f1*Ca_alt(2) - 0.5d0)
				C00_alt(1) = 1/3d0*m0m0*xc0 + 1/6d0*(B012 - B112 + f1*Ca_alt(1) + 1d0)
				C00_alt(2) = 1/3d0*m0m0*xc0 + 1/6d0*(B012 - B112 + f1*Ca_alt(2) + 1d0)
			endif
			
! s-channel or u-channel with one momentum zero and zero momentum transfer C(0, p**2, p**2, m0**2, m1**2, m2**2)
		elseif((Mandelvar.eq.p2p2pat).and.(p1p1pat.eq.0d0)) then
		
			if(p2p2pat.eq.0d0) then
				Ca_alt(1) = 0d0
				Ca_alt(2) = 0d0
				Caa_alt(1) = 0d0
				Caa_alt(2) = 0d0
				C00_alt(1) = 1/4d0*(B012 + m0m0*xc0 + 0.5d0) 
				C00_alt(2) = 1/4d0*(B012 + m0m0*xc0 + 0.5d0)
			else
				Ca_alt(1) = 1/(2d0*p2p2pat)*(B001-B012-f2*xc0)
				Ca_alt(2) = 1/(2d0*p2p2pat)*(B001-B012-f2*xc0) 
				Caa_alt(1) = - 1d0/(3d0*p2p2pat)*(B012 + m0m0*xc0 + 2d0*B112 + 2d0*f2*Ca_alt(1) + 0.5d0) 
				Caa_alt(2) = - 1d0/(3d0*p2p2pat)*(B012 + m0m0*xc0 + 2d0*B112 + 2d0*f2*Ca_alt(2) + 0.5d0)
				C00_alt(1) = 1/3d0*m0m0*xc0 + 1/6d0*(2d0*B012 + B112 + f2*Ca_alt(1) + 1d0)
				C00_alt(2) = 1/3d0*m0m0*xc0 + 1/6d0*(2d0*B012 + B112 + f2*Ca_alt(2) + 1d0)
			endif
				
! s-channel with no momentum zero but zero momentum transfer C(p1**2, (m(p1)+m(p2))**2, p2**2, m0**2, m1**2, m2**2)
		elseif((Mandelvar.eq.((dsqrt(p1p1pat) + dsqrt(p2p2pat))**2).or.(abs(Mandelvar - ((dsqrt(p1p1pat) + dsqrt(p2p2pat))**2)).lt.0.01d0)).and.(p1p1pat.ne.0d0).and.(p2p2pat.ne.0d0)) then
						
			x1(1) = p2p2pat
			x1(2) = -p1p2pat
			x2(1) = -p1p2pat
			x2(2) = p1p1pat
				
			do k=1,2
		
				Y(k) = x1(k)*f1 + x2(k)*f2
			
				Z11(k) = Y(k) + x1(k)/2d0*f1
				Z12(k) = x1(k)/2d0*f2
				Z21(k) = x2(k)/2d0*f1
				Z22(k) = Y(k) + x2(k)/2d0*f2
			
				DetZ(k) = Z11(k)*Z22(k) - Z12(k)*Z21(k)
			
				R7(k) = x1(k)*(B112 + B012 - m0m0*xc0 - (0.5d0*B012 + 0.5d0)) + x2(k)*(B101 + B112 + B012)
				R8(k) = x1(k)*(B102 - B112) + x2(k)*(-B112 - m0m0*xc0 - (0.5d0*B012 + 0.5d0))
		
				C1_alt(k) = (Z22(k)*R7(k) - Z12(k)*R8(k))/DetZ(k)
				C2_alt(k) = (-Z21(k)*R7(k) + Z11(k)*R8(k))/DetZ(k)
				
				C00_alt(k) = 0.5d0*m0m0*xc0 + 0.25d0*(B012 + f1*C1_alt(k) + f2*C2_alt(k) + 1d0)
							
			enddo
			
! Mathematica results
	        C11_alt(1) = -((-(((B012*p1p2pat)/2. - (B1101*p1p2pat)/2. + (B1112*p1p2pat)/2. + B112*p1p2pat + p2p2pat/9. - (B012*p2p2pat)/6. - (B1112*p2p2pat)/2. - 
     -          (2*B112*p2p2pat)/3. - (2*C1_alt(1)*m0m0*p2p2pat)/3.)*(-(f1*f2*p1p1pat*p2p2pat)/18. - 
     -          ((-5*f2*p1p1pat)/6. + (f1*p1p2pat)/2.)*((-2*f2*p1p2pat)/3. + (2*f1*p2p2pat)/3.))) - 
     -         	(f2*p2p2pat*((f2*(p1p1pat/9. - (B1112*p1p1pat)/2. - (B112*p1p1pat)/3. - (2*C2_alt(1)*m0m0*p1p1pat)/3. - (B1102*p1p2pat)/2. + (B1112*p1p2pat)/2.)*p2p2pat)/6. - 
     -          ((-5*f2*p1p1pat)/6. + (f1*p1p2pat)/2.)*(p1p2pat/18. + (B012*p1p2pat)/6. + (B1112*p1p2pat)/2. + (2*B112*p1p2pat)/3. - (C1_alt(1)*m0m0*p1p2pat)/3. - 
     -          p2p2pat/18. - (B1112*p2p2pat)/2. - (B112*p2p2pat)/3. + (C2_alt(1)*m0m0*p2p2pat)/3.)))/3.)/
     -       	((5*f2**3*p1p1pat*p1p2pat**2)/18. - (f1*f2**2*p1p2pat**3)/6. - (2*f1*f2**2*p1p1pat*p1p2pat*p2p2pat)/3. + (5*f1**2*f2*p1p2pat**2*p2p2pat)/12. + 
     -         	(5*f1**2*f2*p1p1pat*p2p2pat**2)/12. - (5*f1**3*p1p2pat*p2p2pat**2)/18.))
		 
	        C12_alt(1) = -(-15*B012*f1*f2*p1p1pat*p1p2pat**2 + 15*B1101*f1*f2*p1p1pat*p1p2pat**2 - 15*B1112*f1*f2*p1p1pat*p1p2pat**2 - 30*B112*f1*f2*p1p1pat*p1p2pat**2 - 
     -        	5*f2**2*p1p1pat*p1p2pat**2 - 15*B012*f2**2*p1p1pat*p1p2pat**2 - 45*B1112*f2**2*p1p1pat*p1p2pat**2 - 60*B112*f2**2*p1p1pat*p1p2pat**2 + 
     -        	30*C1_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat**2 + 9*B012*f1**2*p1p2pat**3 - 9*B1101*f1**2*p1p2pat**3 + 9*B1112*f1**2*p1p2pat**3 + 18*B112*f1**2*p1p2pat**3 + 
     -        	3*f1*f2*p1p2pat**3 + 9*B012*f1*f2*p1p2pat**3 + 27*B1112*f1*f2*p1p2pat**3 + 36*B112*f1*f2*p1p2pat**3 - 18*C1_alt(1)*f1*f2*m0m0*p1p2pat**3 + 
     -        	5*f1*f2*p1p1pat*p1p2pat*p2p2pat + 30*B012*f1*f2*p1p1pat*p1p2pat*p2p2pat + 90*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 
     -        	120*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 3*f2**2*p1p1pat*p1p2pat*p2p2pat + 54*B1112*f2**2*p1p1pat*p1p2pat*p2p2pat + 
     -        	36*B112*f2**2*p1p1pat*p1p2pat*p2p2pat - 30*C1_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat - 18*C2_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat*p2p2pat - 3*f1**2*p1p2pat**2*p2p2pat - 
     -        	18*B012*f1**2*p1p2pat**2*p2p2pat - 54*B1112*f1**2*p1p2pat**2*p2p2pat - 72*B112*f1**2*p1p2pat**2*p2p2pat - 3*f1*f2*p1p2pat**2*p2p2pat - 
     -        	27*B1112*f1*f2*p1p2pat**2*p2p2pat - 18*B112*f1*f2*p1p2pat**2*p2p2pat + 9*B1102*f2**2*p1p2pat**2*p2p2pat - 9*B1112*f2**2*p1p2pat**2*p2p2pat + 
     -        	18*C1_alt(1)*f1**2*m0m0*p1p2pat**2*p2p2pat + 18*C2_alt(1)*f1*f2*m0m0*p1p2pat**2*p2p2pat - 5*f1*f2*p1p1pat*p2p2pat**2 - 90*B1112*f1*f2*p1p1pat*p2p2pat**2 - 
     -        	60*B112*f1*f2*p1p1pat*p2p2pat**2 + 30*C2_alt(1)*f1*f2*m0m0*p1p1pat*p2p2pat**2 + 5*f1**2*p1p2pat*p2p2pat**2 + 45*B1112*f1**2*p1p2pat*p2p2pat**2 + 
     -        	30*B112*f1**2*p1p2pat*p2p2pat**2 - 15*B1102*f1*f2*p1p2pat*p2p2pat**2 + 15*B1112*f1*f2*p1p2pat*p2p2pat**2 - 30*C2_alt(1)*f1**2*m0m0*p1p2pat*p2p2pat**2)/
     -     		(6.*(10*f2**3*p1p1pat*p1p2pat**2 - 6*f1*f2**2*p1p2pat**3 - 24*f1*f2**2*p1p1pat*p1p2pat*p2p2pat + 15*f1**2*f2*p1p2pat**2*p2p2pat + 
     -         	15*f1**2*f2*p1p1pat*p2p2pat**2 - 10*f1**3*p1p2pat*p2p2pat**2))
		 
        	C22_alt(1) = -(3*B012*f1**2*p1p1pat*p1p2pat**2 - 3*B1101*f1**2*p1p1pat*p1p2pat**2 + 3*B1112*f1**2*p1p1pat*p1p2pat**2 + 6*B112*f1**2*p1p1pat*p1p2pat**2 + 
     -        	f1*f2*p1p1pat*p1p2pat**2 + 3*B012*f1*f2*p1p1pat*p1p2pat**2 + 9*B1112*f1*f2*p1p1pat*p1p2pat**2 + 12*B112*f1*f2*p1p1pat*p1p2pat**2 - 
     -       	4*f2**2*p1p1pat*p1p2pat**2 + 18*B1112*f2**2*p1p1pat*p1p2pat**2 + 12*B112*f2**2*p1p1pat*p1p2pat**2 - 6*C1_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat**2 + 
     -        	24*C2_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat**2 + 18*B1102*f2**2*p1p2pat**3 - 18*B1112*f2**2*p1p2pat**3 - f1**2*p1p1pat*p1p2pat*p2p2pat - 
     -        	6*B012*f1**2*p1p1pat*p1p2pat*p2p2pat - 18*B1112*f1**2*p1p1pat*p1p2pat*p2p2pat - 24*B112*f1**2*p1p1pat*p1p2pat*p2p2pat + 9*f1*f2*p1p1pat*p1p2pat*p2p2pat - 
     -        	54*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat - 36*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 6*C1_alt(1)*f1**2*m0m0*p1p1pat*p1p2pat*p2p2pat - 
     -        	54*C2_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat - 45*B1102*f1*f2*p1p2pat**2*p2p2pat + 45*B1112*f1*f2*p1p2pat**2*p2p2pat - 5*f1**2*p1p1pat*p2p2pat**2 + 
     -        	45*B1112*f1**2*p1p1pat*p2p2pat**2 + 30*B112*f1**2*p1p1pat*p2p2pat**2 + 30*C2_alt(1)*f1**2*m0m0*p1p1pat*p2p2pat**2 + 30*B1102*f1**2*p1p2pat*p2p2pat**2 - 
     -        	30*B1112*f1**2*p1p2pat*p2p2pat**2)/
     -     		(3.*(10*f2**3*p1p1pat*p1p2pat**2 - 6*f1*f2**2*p1p2pat**3 - 24*f1*f2**2*p1p1pat*p1p2pat*p2p2pat + 15*f1**2*f2*p1p2pat**2*p2p2pat + 
     -         	15*f1**2*f2*p1p1pat*p2p2pat**2 - 10*f1**3*p1p2pat*p2p2pat**2))
			

	 		C11_alt(2) = -((-((9*B012*p1p1pat - 9*B1101*p1p1pat + 9*B1112*p1p1pat + 18*B112*p1p1pat + 2*p1p2pat - 3*B012*p1p2pat - 9*B1112*p1p2pat - 12*B112*p1p2pat - 
     -          12*C1_alt(2)*m0m0*p1p2pat)*(18*f1*f2*p1p2pat**2 - (12*f2*p1p1pat - 12*f1*p1p2pat)*(-15*f2*p1p2pat + 9*f1*p2p2pat))) - 
     -         	6*f2*p1p2pat*(-3*f2*p1p2pat*(2*p1p2pat - 9*B1112*p1p2pat - 6*B112*p1p2pat - 12*C2_alt(2)*m0m0*p1p2pat - 9*B1102*p2p2pat + 9*B1112*p2p2pat) - 
     -          (-p1p1pat - 3*B012*p1p1pat - 9*B1112*p1p1pat - 12*B112*p1p1pat + 6*C1_alt(2)*m0m0*p1p1pat + p1p2pat + 9*B1112*p1p2pat + 6*B112*p1p2pat - 6*C2_alt(2)*m0m0*p1p2pat)*
     -          (-15*f2*p1p2pat + 9*f1*p2p2pat)))/
     -       	(18*f1*f2*p1p1pat*p1p2pat*(-15*f2*p1p2pat + 9*f1*p2p2pat) - 
     -         	(9*f2*p1p1pat - 15*f1*p1p2pat)*(18*f1*f2*p1p2pat**2 - (12*f2*p1p1pat - 12*f1*p1p2pat)*(-15*f2*p1p2pat + 9*f1*p2p2pat))))
	 
        	C12_alt(2) = -(-15*B012*f1*f2*p1p1pat**2*p1p2pat + 15*B1101*f1*f2*p1p1pat**2*p1p2pat - 15*B1112*f1*f2*p1p1pat**2*p1p2pat - 30*B112*f1*f2*p1p1pat**2*p1p2pat - 
     -        	5*f2**2*p1p1pat**2*p1p2pat - 15*B012*f2**2*p1p1pat**2*p1p2pat - 45*B1112*f2**2*p1p1pat**2*p1p2pat - 60*B112*f2**2*p1p1pat**2*p1p2pat + 
     -        	30*C1_alt(2)*f2**2*m0m0*p1p1pat**2*p1p2pat + 5*f1*f2*p1p1pat*p1p2pat**2 + 30*B012*f1*f2*p1p1pat*p1p2pat**2 + 90*B1112*f1*f2*p1p1pat*p1p2pat**2 + 
     -        	120*B112*f1*f2*p1p1pat*p1p2pat**2 + 3*f2**2*p1p1pat*p1p2pat**2 + 54*B1112*f2**2*p1p1pat*p1p2pat**2 + 36*B112*f2**2*p1p1pat*p1p2pat**2 - 
     -        	30*C1_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat**2 - 18*C2_alt(2)*f2**2*m0m0*p1p1pat*p1p2pat**2 - 5*f1*f2*p1p2pat**3 - 90*B1112*f1*f2*p1p2pat**3 - 60*B112*f1*f2*p1p2pat**3 + 
     -        	30*C2_alt(2)*f1*f2*m0m0*p1p2pat**3 + 9*B012*f1**2*p1p1pat**2*p2p2pat - 9*B1101*f1**2*p1p1pat**2*p2p2pat + 9*B1112*f1**2*p1p1pat**2*p2p2pat + 
     -        	18*B112*f1**2*p1p1pat**2*p2p2pat + 3*f1*f2*p1p1pat**2*p2p2pat + 9*B012*f1*f2*p1p1pat**2*p2p2pat + 27*B1112*f1*f2*p1p1pat**2*p2p2pat + 
     -        	36*B112*f1*f2*p1p1pat**2*p2p2pat - 18*C1_alt(2)*f1*f2*m0m0*p1p1pat**2*p2p2pat - 3*f1**2*p1p1pat*p1p2pat*p2p2pat - 18*B012*f1**2*p1p1pat*p1p2pat*p2p2pat - 
     -        	54*B1112*f1**2*p1p1pat*p1p2pat*p2p2pat - 72*B112*f1**2*p1p1pat*p1p2pat*p2p2pat - 3*f1*f2*p1p1pat*p1p2pat*p2p2pat - 
     -        	27*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat - 18*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 9*B1102*f2**2*p1p1pat*p1p2pat*p2p2pat - 
     -        	9*B1112*f2**2*p1p1pat*p1p2pat*p2p2pat + 18*C1_alt(2)*f1**2*m0m0*p1p1pat*p1p2pat*p2p2pat + 18*C2_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat + 5*f1**2*p1p2pat**2*p2p2pat + 
     -        	45*B1112*f1**2*p1p2pat**2*p2p2pat + 30*B112*f1**2*p1p2pat**2*p2p2pat - 15*B1102*f1*f2*p1p2pat**2*p2p2pat + 15*B1112*f1*f2*p1p2pat**2*p2p2pat - 
     -        	30*C2_alt(2)*f1**2*m0m0*p1p2pat**2*p2p2pat)/
     -     	  	(6.*(10*f2**3*p1p1pat**2*p1p2pat - 24*f1*f2**2*p1p1pat*p1p2pat**2 + 15*f1**2*f2*p1p2pat**3 - 6*f1*f2**2*p1p1pat**2*p2p2pat + 
     -        	15*f1**2*f2*p1p1pat*p1p2pat*p2p2pat - 10*f1**3*p1p2pat**2*p2p2pat))
	 
     	  	C22_alt(2) = -(3*B012*f1**2*p1p1pat**2*p1p2pat - 3*B1101*f1**2*p1p1pat**2*p1p2pat + 3*B1112*f1**2*p1p1pat**2*p1p2pat + 6*B112*f1**2*p1p1pat**2*p1p2pat + 
     -        	f1*f2*p1p1pat**2*p1p2pat + 3*B012*f1*f2*p1p1pat**2*p1p2pat + 9*B1112*f1*f2*p1p1pat**2*p1p2pat + 12*B112*f1*f2*p1p1pat**2*p1p2pat - 
     -        	4*f2**2*p1p1pat**2*p1p2pat + 18*B1112*f2**2*p1p1pat**2*p1p2pat + 12*B112*f2**2*p1p1pat**2*p1p2pat - 6*C1_alt(2)*f1*f2*m0m0*p1p1pat**2*p1p2pat + 
     -        	24*C2_alt(2)*f2**2*m0m0*p1p1pat**2*p1p2pat - f1**2*p1p1pat*p1p2pat**2 - 6*B012*f1**2*p1p1pat*p1p2pat**2 - 18*B1112*f1**2*p1p1pat*p1p2pat**2 - 
     -        	24*B112*f1**2*p1p1pat*p1p2pat**2 + 9*f1*f2*p1p1pat*p1p2pat**2 - 54*B1112*f1*f2*p1p1pat*p1p2pat**2 - 36*B112*f1*f2*p1p1pat*p1p2pat**2 + 
     -        	6*C1_alt(2)*f1**2*m0m0*p1p1pat*p1p2pat**2 - 54*C2_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat**2 - 5*f1**2*p1p2pat**3 + 45*B1112*f1**2*p1p2pat**3 + 30*B112*f1**2*p1p2pat**3 + 
     -        	30*C2_alt(2)*f1**2*m0m0*p1p2pat**3 + 18*B1102*f2**2*p1p1pat**2*p2p2pat - 18*B1112*f2**2*p1p1pat**2*p2p2pat - 45*B1102*f1*f2*p1p1pat*p1p2pat*p2p2pat + 
     -        	45*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 30*B1102*f1**2*p1p2pat**2*p2p2pat - 30*B1112*f1**2*p1p2pat**2*p2p2pat)/
     -     	  	(3.*(10*f2**3*p1p1pat**2*p1p2pat - 24*f1*f2**2*p1p1pat*p1p2pat**2 + 15*f1**2*f2*p1p2pat**3 - 6*f1*f2**2*p1p1pat**2*p2p2pat + 
     -        	15*f1**2*f2*p1p1pat*p1p2pat*p2p2pat - 10*f1**3*p1p2pat**2*p2p2pat))

			
! u-channel with no momentum zero but zero momentum transfer C(p1**2, (m(p1)-m(p2))**2, p2**2, m0**2, m1**2, m2**2)
		elseif((Mandelvar.eq.((dsqrt(p1p1pat) - dsqrt(p2p2pat))**2).or.(abs(Mandelvar - ((dsqrt(p1p1pat) - dsqrt(p2p2pat))**2)).lt.0.01d0)).and.(p1p1pat.ne.0d0).and.(p2p2pat.ne.0d0)) then
		
			x1(1) = p2p2pat
			x1(2) = -p1p2pat
			x2(1) = -p1p2pat
			x2(2) = p1p1pat
				
			do k=1,2
		
				Y(k) = x1(k)*f1 + x2(k)*f2
			
				Z11(k) = Y(k) + x1(k)/2d0*f1
				Z12(k) = x1(k)/2d0*f2
				Z21(k) = x2(k)/2d0*f1
				Z22(k) = Y(k) + x2(k)/2d0*f2
			
				DetZ(k) = Z11(k)*Z22(k) - Z12(k)*Z21(k)
			
				R7(k) = x1(k)*(B112 + B012 - m0m0*xc0 - (0.5d0*B012 + 0.5d0)) + x2(k)*(B101 + B112 + B012)
				R8(k) = x1(k)*(B102 - B112) + x2(k)*(-B112 - m0m0*xc0 - (0.5d0*B012 + 0.5d0))
		
				C1_alt(k) = (Z22(k)*R7(k) - Z12(k)*R8(k))/DetZ(k)
				C2_alt(k) = (-Z21(k)*R7(k) + Z11(k)*R8(k))/DetZ(k)
				
				C00_alt(k) = 0.5d0*m0m0*xc0 + 0.25d0*(B012 + f1*C1_alt(k) + f2*C2_alt(k) + 1d0) 
				
			enddo

! Mathematica results
	        C11_alt(1) = -((-(((B012*p1p2pat)/2. - (B1101*p1p2pat)/2. + (B1112*p1p2pat)/2. + B112*p1p2pat + p2p2pat/9. - (B012*p2p2pat)/6. - (B1112*p2p2pat)/2. - 
     -          (2*B112*p2p2pat)/3. - (2*C1_alt(1)*m0m0*p2p2pat)/3.)*(-(f1*f2*p1p1pat*p2p2pat)/18. - 
     -          ((-5*f2*p1p1pat)/6. + (f1*p1p2pat)/2.)*((-2*f2*p1p2pat)/3. + (2*f1*p2p2pat)/3.))) - 
     -         	(f2*p2p2pat*((f2*(p1p1pat/9. - (B1112*p1p1pat)/2. - (B112*p1p1pat)/3. - (2*C2_alt(1)*m0m0*p1p1pat)/3. - (B1102*p1p2pat)/2. + (B1112*p1p2pat)/2.)*p2p2pat)/6. - 
     -          ((-5*f2*p1p1pat)/6. + (f1*p1p2pat)/2.)*(p1p2pat/18. + (B012*p1p2pat)/6. + (B1112*p1p2pat)/2. + (2*B112*p1p2pat)/3. - (C1_alt(1)*m0m0*p1p2pat)/3. - 
     -          p2p2pat/18. - (B1112*p2p2pat)/2. - (B112*p2p2pat)/3. + (C2_alt(1)*m0m0*p2p2pat)/3.)))/3.)/
     -       	((5*f2**3*p1p1pat*p1p2pat**2)/18. - (f1*f2**2*p1p2pat**3)/6. - (2*f1*f2**2*p1p1pat*p1p2pat*p2p2pat)/3. + (5*f1**2*f2*p1p2pat**2*p2p2pat)/12. + 
     -         	(5*f1**2*f2*p1p1pat*p2p2pat**2)/12. - (5*f1**3*p1p2pat*p2p2pat**2)/18.))
		 
	        C12_alt(1) = -(-15*B012*f1*f2*p1p1pat*p1p2pat**2 + 15*B1101*f1*f2*p1p1pat*p1p2pat**2 - 15*B1112*f1*f2*p1p1pat*p1p2pat**2 - 30*B112*f1*f2*p1p1pat*p1p2pat**2 - 
     -        	5*f2**2*p1p1pat*p1p2pat**2 - 15*B012*f2**2*p1p1pat*p1p2pat**2 - 45*B1112*f2**2*p1p1pat*p1p2pat**2 - 60*B112*f2**2*p1p1pat*p1p2pat**2 + 
     -        	30*C1_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat**2 + 9*B012*f1**2*p1p2pat**3 - 9*B1101*f1**2*p1p2pat**3 + 9*B1112*f1**2*p1p2pat**3 + 18*B112*f1**2*p1p2pat**3 + 
     -        	3*f1*f2*p1p2pat**3 + 9*B012*f1*f2*p1p2pat**3 + 27*B1112*f1*f2*p1p2pat**3 + 36*B112*f1*f2*p1p2pat**3 - 18*C1_alt(1)*f1*f2*m0m0*p1p2pat**3 + 
     -        	5*f1*f2*p1p1pat*p1p2pat*p2p2pat + 30*B012*f1*f2*p1p1pat*p1p2pat*p2p2pat + 90*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 
     -        	120*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 3*f2**2*p1p1pat*p1p2pat*p2p2pat + 54*B1112*f2**2*p1p1pat*p1p2pat*p2p2pat + 
     -        	36*B112*f2**2*p1p1pat*p1p2pat*p2p2pat - 30*C1_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat - 18*C2_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat*p2p2pat - 3*f1**2*p1p2pat**2*p2p2pat - 
     -        	18*B012*f1**2*p1p2pat**2*p2p2pat - 54*B1112*f1**2*p1p2pat**2*p2p2pat - 72*B112*f1**2*p1p2pat**2*p2p2pat - 3*f1*f2*p1p2pat**2*p2p2pat - 
     -        	27*B1112*f1*f2*p1p2pat**2*p2p2pat - 18*B112*f1*f2*p1p2pat**2*p2p2pat + 9*B1102*f2**2*p1p2pat**2*p2p2pat - 9*B1112*f2**2*p1p2pat**2*p2p2pat + 
     -        	18*C1_alt(1)*f1**2*m0m0*p1p2pat**2*p2p2pat + 18*C2_alt(1)*f1*f2*m0m0*p1p2pat**2*p2p2pat - 5*f1*f2*p1p1pat*p2p2pat**2 - 90*B1112*f1*f2*p1p1pat*p2p2pat**2 - 
     -        	60*B112*f1*f2*p1p1pat*p2p2pat**2 + 30*C2_alt(1)*f1*f2*m0m0*p1p1pat*p2p2pat**2 + 5*f1**2*p1p2pat*p2p2pat**2 + 45*B1112*f1**2*p1p2pat*p2p2pat**2 + 
     -        	30*B112*f1**2*p1p2pat*p2p2pat**2 - 15*B1102*f1*f2*p1p2pat*p2p2pat**2 + 15*B1112*f1*f2*p1p2pat*p2p2pat**2 - 30*C2_alt(1)*f1**2*m0m0*p1p2pat*p2p2pat**2)/
     -     		(6.*(10*f2**3*p1p1pat*p1p2pat**2 - 6*f1*f2**2*p1p2pat**3 - 24*f1*f2**2*p1p1pat*p1p2pat*p2p2pat + 15*f1**2*f2*p1p2pat**2*p2p2pat + 
     -         	15*f1**2*f2*p1p1pat*p2p2pat**2 - 10*f1**3*p1p2pat*p2p2pat**2))

        	C22_alt(1) = -(3*B012*f1**2*p1p1pat*p1p2pat**2 - 3*B1101*f1**2*p1p1pat*p1p2pat**2 + 3*B1112*f1**2*p1p1pat*p1p2pat**2 + 6*B112*f1**2*p1p1pat*p1p2pat**2 + 
     -        	f1*f2*p1p1pat*p1p2pat**2 + 3*B012*f1*f2*p1p1pat*p1p2pat**2 + 9*B1112*f1*f2*p1p1pat*p1p2pat**2 + 12*B112*f1*f2*p1p1pat*p1p2pat**2 - 
     -       	4*f2**2*p1p1pat*p1p2pat**2 + 18*B1112*f2**2*p1p1pat*p1p2pat**2 + 12*B112*f2**2*p1p1pat*p1p2pat**2 - 6*C1_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat**2 + 
     -        	24*C2_alt(1)*f2**2*m0m0*p1p1pat*p1p2pat**2 + 18*B1102*f2**2*p1p2pat**3 - 18*B1112*f2**2*p1p2pat**3 - f1**2*p1p1pat*p1p2pat*p2p2pat - 
     -        	6*B012*f1**2*p1p1pat*p1p2pat*p2p2pat - 18*B1112*f1**2*p1p1pat*p1p2pat*p2p2pat - 24*B112*f1**2*p1p1pat*p1p2pat*p2p2pat + 9*f1*f2*p1p1pat*p1p2pat*p2p2pat - 
     -        	54*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat - 36*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 6*C1_alt(1)*f1**2*m0m0*p1p1pat*p1p2pat*p2p2pat - 
     -        	54*C2_alt(1)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat - 45*B1102*f1*f2*p1p2pat**2*p2p2pat + 45*B1112*f1*f2*p1p2pat**2*p2p2pat - 5*f1**2*p1p1pat*p2p2pat**2 + 
     -        	45*B1112*f1**2*p1p1pat*p2p2pat**2 + 30*B112*f1**2*p1p1pat*p2p2pat**2 + 30*C2_alt(1)*f1**2*m0m0*p1p1pat*p2p2pat**2 + 30*B1102*f1**2*p1p2pat*p2p2pat**2 - 
     -        	30*B1112*f1**2*p1p2pat*p2p2pat**2)/
     -     		(3.*(10*f2**3*p1p1pat*p1p2pat**2 - 6*f1*f2**2*p1p2pat**3 - 24*f1*f2**2*p1p1pat*p1p2pat*p2p2pat + 15*f1**2*f2*p1p2pat**2*p2p2pat + 
     -         	15*f1**2*f2*p1p1pat*p2p2pat**2 - 10*f1**3*p1p2pat*p2p2pat**2))


	 		C11_alt(2) = -((-((9*B012*p1p1pat - 9*B1101*p1p1pat + 9*B1112*p1p1pat + 18*B112*p1p1pat + 2*p1p2pat - 3*B012*p1p2pat - 9*B1112*p1p2pat - 12*B112*p1p2pat - 
     -          12*C1_alt(2)*m0m0*p1p2pat)*(18*f1*f2*p1p2pat**2 - (12*f2*p1p1pat - 12*f1*p1p2pat)*(-15*f2*p1p2pat + 9*f1*p2p2pat))) - 
     -         	6*f2*p1p2pat*(-3*f2*p1p2pat*(2*p1p2pat - 9*B1112*p1p2pat - 6*B112*p1p2pat - 12*C2_alt(2)*m0m0*p1p2pat - 9*B1102*p2p2pat + 9*B1112*p2p2pat) - 
     -          (-p1p1pat - 3*B012*p1p1pat - 9*B1112*p1p1pat - 12*B112*p1p1pat + 6*C1_alt(2)*m0m0*p1p1pat + p1p2pat + 9*B1112*p1p2pat + 6*B112*p1p2pat - 6*C2_alt(2)*m0m0*p1p2pat)*
     -          (-15*f2*p1p2pat + 9*f1*p2p2pat)))/
     -       	(18*f1*f2*p1p1pat*p1p2pat*(-15*f2*p1p2pat + 9*f1*p2p2pat) - 
     -         	(9*f2*p1p1pat - 15*f1*p1p2pat)*(18*f1*f2*p1p2pat**2 - (12*f2*p1p1pat - 12*f1*p1p2pat)*(-15*f2*p1p2pat + 9*f1*p2p2pat))))

        	C12_alt(2) = -(-15*B012*f1*f2*p1p1pat**2*p1p2pat + 15*B1101*f1*f2*p1p1pat**2*p1p2pat - 15*B1112*f1*f2*p1p1pat**2*p1p2pat - 30*B112*f1*f2*p1p1pat**2*p1p2pat - 
     -        	5*f2**2*p1p1pat**2*p1p2pat - 15*B012*f2**2*p1p1pat**2*p1p2pat - 45*B1112*f2**2*p1p1pat**2*p1p2pat - 60*B112*f2**2*p1p1pat**2*p1p2pat + 
     -        	30*C1_alt(2)*f2**2*m0m0*p1p1pat**2*p1p2pat + 5*f1*f2*p1p1pat*p1p2pat**2 + 30*B012*f1*f2*p1p1pat*p1p2pat**2 + 90*B1112*f1*f2*p1p1pat*p1p2pat**2 + 
     -        	120*B112*f1*f2*p1p1pat*p1p2pat**2 + 3*f2**2*p1p1pat*p1p2pat**2 + 54*B1112*f2**2*p1p1pat*p1p2pat**2 + 36*B112*f2**2*p1p1pat*p1p2pat**2 - 
     -        	30*C1_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat**2 - 18*C2_alt(2)*f2**2*m0m0*p1p1pat*p1p2pat**2 - 5*f1*f2*p1p2pat**3 - 90*B1112*f1*f2*p1p2pat**3 - 60*B112*f1*f2*p1p2pat**3 + 
     -        	30*C2_alt(2)*f1*f2*m0m0*p1p2pat**3 + 9*B012*f1**2*p1p1pat**2*p2p2pat - 9*B1101*f1**2*p1p1pat**2*p2p2pat + 9*B1112*f1**2*p1p1pat**2*p2p2pat + 
     -        	18*B112*f1**2*p1p1pat**2*p2p2pat + 3*f1*f2*p1p1pat**2*p2p2pat + 9*B012*f1*f2*p1p1pat**2*p2p2pat + 27*B1112*f1*f2*p1p1pat**2*p2p2pat + 
     -        	36*B112*f1*f2*p1p1pat**2*p2p2pat - 18*C1_alt(2)*f1*f2*m0m0*p1p1pat**2*p2p2pat - 3*f1**2*p1p1pat*p1p2pat*p2p2pat - 18*B012*f1**2*p1p1pat*p1p2pat*p2p2pat - 
     -        	54*B1112*f1**2*p1p1pat*p1p2pat*p2p2pat - 72*B112*f1**2*p1p1pat*p1p2pat*p2p2pat - 3*f1*f2*p1p1pat*p1p2pat*p2p2pat - 
     -        	27*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat - 18*B112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 9*B1102*f2**2*p1p1pat*p1p2pat*p2p2pat - 
     -        	9*B1112*f2**2*p1p1pat*p1p2pat*p2p2pat + 18*C1_alt(2)*f1**2*m0m0*p1p1pat*p1p2pat*p2p2pat + 18*C2_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat*p2p2pat + 5*f1**2*p1p2pat**2*p2p2pat + 
     -        	45*B1112*f1**2*p1p2pat**2*p2p2pat + 30*B112*f1**2*p1p2pat**2*p2p2pat - 15*B1102*f1*f2*p1p2pat**2*p2p2pat + 15*B1112*f1*f2*p1p2pat**2*p2p2pat - 
     -        	30*C2_alt(2)*f1**2*m0m0*p1p2pat**2*p2p2pat)/
     -     	  	(6.*(10*f2**3*p1p1pat**2*p1p2pat - 24*f1*f2**2*p1p1pat*p1p2pat**2 + 15*f1**2*f2*p1p2pat**3 - 6*f1*f2**2*p1p1pat**2*p2p2pat + 
     -        	15*f1**2*f2*p1p1pat*p1p2pat*p2p2pat - 10*f1**3*p1p2pat**2*p2p2pat))

     	  	C22_alt(2) = -(3*B012*f1**2*p1p1pat**2*p1p2pat - 3*B1101*f1**2*p1p1pat**2*p1p2pat + 3*B1112*f1**2*p1p1pat**2*p1p2pat + 6*B112*f1**2*p1p1pat**2*p1p2pat + 
     -        	f1*f2*p1p1pat**2*p1p2pat + 3*B012*f1*f2*p1p1pat**2*p1p2pat + 9*B1112*f1*f2*p1p1pat**2*p1p2pat + 12*B112*f1*f2*p1p1pat**2*p1p2pat - 
     -        	4*f2**2*p1p1pat**2*p1p2pat + 18*B1112*f2**2*p1p1pat**2*p1p2pat + 12*B112*f2**2*p1p1pat**2*p1p2pat - 6*C1_alt(2)*f1*f2*m0m0*p1p1pat**2*p1p2pat + 
     -        	24*C2_alt(2)*f2**2*m0m0*p1p1pat**2*p1p2pat - f1**2*p1p1pat*p1p2pat**2 - 6*B012*f1**2*p1p1pat*p1p2pat**2 - 18*B1112*f1**2*p1p1pat*p1p2pat**2 - 
     -        	24*B112*f1**2*p1p1pat*p1p2pat**2 + 9*f1*f2*p1p1pat*p1p2pat**2 - 54*B1112*f1*f2*p1p1pat*p1p2pat**2 - 36*B112*f1*f2*p1p1pat*p1p2pat**2 + 
     -        	6*C1_alt(2)*f1**2*m0m0*p1p1pat*p1p2pat**2 - 54*C2_alt(2)*f1*f2*m0m0*p1p1pat*p1p2pat**2 - 5*f1**2*p1p2pat**3 + 45*B1112*f1**2*p1p2pat**3 + 30*B112*f1**2*p1p2pat**3 + 
     -        	30*C2_alt(2)*f1**2*m0m0*p1p2pat**3 + 18*B1102*f2**2*p1p1pat**2*p2p2pat - 18*B1112*f2**2*p1p1pat**2*p2p2pat - 45*B1102*f1*f2*p1p1pat*p1p2pat*p2p2pat + 
     -        	45*B1112*f1*f2*p1p1pat*p1p2pat*p2p2pat + 30*B1102*f1**2*p1p2pat**2*p2p2pat - 30*B1112*f1**2*p1p2pat**2*p2p2pat)/
     -     	  	(3.*(10*f2**3*p1p1pat**2*p1p2pat - 24*f1*f2**2*p1p1pat*p1p2pat**2 + 15*f1**2*f2*p1p2pat**3 - 6*f1*f2**2*p1p1pat**2*p2p2pat + 
     -        	15*f1**2*f2*p1p1pat*p1p2pat*p2p2pat - 10*f1**3*p1p2pat**2*p2p2pat))

		else
	
			write(*,*)'Error in alternative C tensor determination - Case not found'

		endif
				
	results(1) = C0_alt(1)
	results(2) = C0_alt(2)
	results(3) = C1_alt(1)
	results(4) = C1_alt(2)
	results(5) = C2_alt(1)
	results(6) = C2_alt(2)
	results(7) = Ca_alt(1)
	results(8) = Ca_alt(2)
	results(9) = C00_alt(1)
	results(10) = C00_alt(2)
	results(11) = Caa_alt(1)
	results(12) = Caa_alt(2)
	results(13) = C11_alt(1)
	results(14) = C11_alt(2)
	results(15) = C12_alt(1)
	results(16) = C12_alt(2)
	results(17) = C22_alt(1)
	results(18) = C22_alt(2)
			
!======================= END OF PROGRAM ===========================================================

	end

