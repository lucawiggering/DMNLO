!==================================================================================================
!======================= Initialision for direct detection in DM@NLO ==============================
!=======================        Author: Patrick Steppeler		     ==============================
!==================================================================================================

	subroutine DD_Init
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_Model.h"
#include "DD_GenCouplings.h"
#include "DD_GenCounterterms.h"
	
!----------------------- Declare variables --------------------------------------------------------
	integer nucleon 			! 1 = proton,	2 = neutron
	integer uptype				! 1 = up type,	2 = down type
	integer gen					! Quark generation
	integer contrib 			! 1 = tree, 2 = vertex, 3 = vertex counter, 4 = prop, 5 = prop counter, 6 = box
								! 7 = treeHO, 8 = EFTVertex, 9 = EFTVertexCounter
	integer tt					! 1 = neutrinos, 2 = electrons (e, mu, tau), 3 = up quarks, 4 = down quarks
	double precision MBDR,MTDR
	double complex dMf_Temp		! To avoid that dMf(3,3) gets replaced when calling MTDR
	
!======================= ACTUAL CALCULATIONS ======================================================

!----------------------- Set nuclear factors ------------------------------------------------------
	do nucleon = 1,2
		do uptype = 1,2
			do gen = 1,3
				f_T(nucleon,uptype,gen) = 0d0
				Delta_q(nucleon,uptype,gen) = 0d0
			enddo
		enddo
	enddo

!----------------------- Spin dependent
! micrOMEGAs & DM@NLO setup 
	if((Nucl_Mode.eq.1).or.(Nucl_Mode.eq.2).or.(Nucl_Mode.eq.4)) then
		Delta_q(1,1,1) = 0.842d0
		Delta_q(1,2,1) = -0.427d0
		Delta_q(1,2,2) = -0.085d0

		Delta_q(2,1,1) = -0.427d0
		Delta_q(2,2,1) = 0.842d0
		Delta_q(2,2,2) = -0.085d0

! DarkSUSY setup
	elseif(Nucl_Mode.eq.3) then
		Delta_q(1,1,1) = 0.77d0
		Delta_q(1,2,1) = -0.40d0
		Delta_q(1,2,2) = -0.12d0

		Delta_q(2,1,1) = -0.40d0
		Delta_q(2,2,1) = 0.77d0
		Delta_q(2,2,2) = -0.12d0
	
	else
		write(*,*)'WARNING! Nuclear input not set.'
	endif
	
!	Delta_q(1,1,3) = 0.1d0	!To simulate SD for massive quarks

!----------------------- Spin independent
! micrOMEGAs setup 2.4.1
	if(Nucl_Mode.eq.1) then
		f_T(1,1,1) = 0.023d0
		f_T(1,2,1) = 0.033d0
		f_T(1,1,2) = 0.0507d0
		f_T(1,2,2) = 0.26d0
		f_T(1,1,3) = 0.0507d0
		f_T(1,2,3) = 0.0507d0

		f_T(2,1,1) = 0.018d0
		f_T(2,2,1) = 0.042d0
		f_T(2,1,2) = 0.0504d0
		f_T(2,2,2) = 0.26d0
		f_T(2,1,3) = 0.0504d0
		f_T(2,2,3) = 0.0504d0

! micrOMEGAs setup 4.0.3
	elseif(Nucl_Mode.eq.2) then
		f_T(1,1,1) = 0.0153d0
		f_T(1,2,1) = 0.0191d0
		f_T(1,1,2) = 0.0682d0
		f_T(1,2,2) = 0.0447d0
		f_T(1,1,3) = 0.0682d0
		f_T(1,2,3) = 0.0682d0

		f_T(2,1,1) = 0.0110d0
		f_T(2,2,1) = 0.0273d0
		f_T(2,1,2) = 0.0679d0
		f_T(2,2,2) = 0.0447d0
		f_T(2,1,3) = 0.0679d0
		f_T(2,2,3) = 0.0679d0
					
! DarkSUSY setup
	elseif(Nucl_Mode.eq.3) then
		f_T(1,1,1) = 0.023d0
		f_T(1,2,1) = 0.034d0
		f_T(1,1,2) = 0.0595d0
		f_T(1,2,2) = 0.14d0
		f_T(1,1,3) = 0.0595d0
		f_T(1,2,3) = 0.0595d0

		f_T(2,1,1) = 0.019d0
		f_T(2,2,1) = 0.041d0
		f_T(2,1,2) = 0.0592d0
		f_T(2,2,2) = 0.14d0
		f_T(2,1,3) = 0.0592d0
		f_T(2,2,3) = 0.0592d0

! Crivellin/Hoferichter setup
	elseif(Nucl_Mode.eq.4) then
		f_T(1,1,1) = 0.0208d0
		f_T(1,2,1) = 0.0411d0
		f_T(1,1,2) = 0.0663d0
		f_T(1,2,2) = 0.043d0
		f_T(1,1,3) = 0.0663d0
		f_T(1,2,3) = 0.0663d0

		f_T(2,1,1) = 0.0189d0
		f_T(2,2,1) = 0.0451d0
		f_T(2,1,2) = 0.0661d0
		f_T(2,2,2) = 0.043d0
		f_T(2,1,3) = 0.0661d0
		f_T(2,2,3) = 0.0661d0
		
	else
		write(*,*)'WARNING! Nuclear input not set.'
	endif
	
!----------------------- Set nuclear masses -------------------------------------------------------
	m_Proton = 0.938271998d0
	m_Neutron = 0.9396d0
	m_Nucleon(1) = m_Proton
	m_Nucleon(2) = m_Neutron
	
	do nucleon=1,2
		mu_Nucleon(nucleon) = m_Nucleon(nucleon)*MNeu(1)/(MNeu(1) + m_Nucleon(nucleon))
	enddo
!----------------------- Set four fermion couplings to zero ---------------------------------------
	do nucleon = 1,2
		g_4Fs(nucleon) = (0d0, 0d0)
		g_4Fa(nucleon) = (0d0, 0d0)
		do contrib = 1,9
			g_4FaContribs(nucleon, contrib) = (0d0, 0d0)
			g_4FsContribs(nucleon, contrib) = (0d0, 0d0)
		enddo
	enddo
	
	g_4Fsp = g_4Fs(1)
	g_4Fsn = g_4Fs(2)
	g_4Fap = g_4Fa(1)
	g_4Fan = g_4Fa(2)
	
!----------------------- Low scale for matching ---------------------------------------------------
	Q_low = 5d0

!----------------------- Masses of quarks ---------------------------------------------------------

! Masses info:
! Mf		- fermion mass at Q_high, zero for gen = 1, 2
! Mf_high	- fermion mass at Q_high, never zero (for m_t Mf_high might be to on-shell mass, 
!			  depending on the choice of the scheme)
! Mf_low	- fermion mass at Q_low, never zero

! CalcHEP values
! The code depends on the input values for the light quark masses. However, this dependence is neglible. Multiplying by a
! factor of 10 doesn't change the final result
	Mf_high(3,1) = 0.0056d0
	Mf_high(4,1) = 0.0099d0
	Mf_high(3,2) = 0.71914743546778404d0
	Mf_high(4,2) = 0.18736986174637255d0
	Mf_high(3,3) = Mf(3,3)
	Mf_high(4,3) = Mf(4,3)
	
!In principle low scale values for the first generation quarks are needed. In practice this value cancels out exactly anyhow.
	Mf_low(3,1) = 0.0056d0
	Mf_low(4,1) = 0.0099d0
	Mf_low(3,2) = 0.71914743546778404d0
	Mf_low(4,2) = 0.18736986174637255d0
	dMf_Temp = dMf(3,3)
	Mf_low(3,3) = MTDR(MTos,Q_low)
	dMf(3,3) = dMf_Temp
	Mf_low(4,3) = MBDR(MBSM,MBSM**2,Q_low**2,5,3)
		
C 	do tt = 3,4
C 		do gen = 1,3
C 			write(*,*)'tt, gen', tt, gen
C 			write(*,*)'Mf = ',Mf(tt,gen)
C 			write(*,*)'Mf_high = ',Mf_high(tt,gen)
C 			write(*,*)'Mf_low = ',Mf_low(tt,gen)
C 		enddo
C 	enddo
C
C 	write(*,*)'MBDR @ 1 TeV',MBDR(MBSM,MBSM**2,QScale**2,5,3)
C 	write(*,*)'MTDR @ 1 TeV',MTDR(MTos,QScale)
	
	
!======================= END OF PROGRAM ===========================================================

	end


