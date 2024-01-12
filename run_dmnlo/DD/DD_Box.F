!==================================================================================================
!======================= Calculates box contributions to four fermion couplings ===================
!=======================				Author: Patrick Steppeler	   	        ===================
!==================================================================================================

	subroutine DD_Box	
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_Model.h"
#include "DD_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"
	
!----------------------- Declare variables --------------------------------------------------------
	complex*16 axialfactor(2,3), scalarfactor(2,3)
	integer nucleon, uptype, gen, tt, i, j
	logical DD_udsCheck
		
!======================= ACTUAL CALCULATIONS ======================================================

!======================= UNSUPPRESSED CONTRIBUTIONS ===============================================

!----------------------- Unsuppressed scalar contributions ----------------------------------------
	do gen = 1,3
		do uptype = 1,2
			scalarfactor(uptype,gen) = 0.0d0
		enddo
	enddo

!----------------------- B1 contribution
! Note: The old box mode (0) activates the old amplitude which is wrong in the limit v -> 0. This
! option is kept for testing purposes.

	do gen = 1,3
!	gen = 3
		do uptype = 1,2
			do i = 1,2
			tt = uptype + 2

			call myDget(s_Mandel(tt,gen),MNeu(1)**2,t_Mandel(tt,gen),Mf(tt,gen)**2,Mf(tt,gen)**2,MNeu(1)**2,0d0,MSf(i,tt,gen)**2,Mf(tt,gen)**2,Mf(tt,gen)**2)

			if(Box_Mode.eq.0) then
				scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - (Alfas*CF*(-(Mf(tt,gen)*(FR(1,i,tt,gen)*((xd1 + xd11 + 2d0*xd12 + 2d0*xd13 + xd2 + xd22 + 2d0*xd23 + xd3 + xd33)*ER(1,i,tt,gen)*Mf(tt,gen) + 2d0*xd12*ELc(1,i,tt,gen)*MNeu(1)) +
     -          	FL(1,i,tt,gen)*((xd1 + xd11 + 2d0*xd12 + 2d0*xd13 + xd2 + xd22 + 2d0*xd23 + xd3 + xd33)*ELc(1,i,tt,gen)*Mf(tt,gen) + 2d0*xd12*ER(1,i,tt,gen)*MNeu(1)))) -
     -       		2d0*(xd0 + xd1 + xd12 + xd2 + xd23 + xd3)*(ELc(1,i,tt,gen)*FL(1,i,tt,gen) + ER(1,i,tt,gen)*FR(1,i,tt,gen))*pbp2(tt,gen)) - 8d0*(FR(1,i,tt,gen)*ER(1,i,tt,gen) + FL(1,i,tt,gen)*ELc(1,i,tt,gen))*xd00)/(16.*pi)
 			else
				scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - (Alfas*CF*(-2*(FL(1,i,tt,gen)*((-xd1 + xd11 + xd12 + xd13)*ER(1,i,tt,gen)*Mf(tt,gen)*MNeu(1) +
     -             	ELc(1,i,tt,gen)*(4*xd00 + (xd1 + xd11 + 2*xd12 + 2*xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33)*Mf(tt,gen)**2 + xd11*MNeu(1)**2)) +
     -          	FR(1,i,tt,gen)*((-xd1 + xd11 + xd12 + xd13)*ELc(1,i,tt,gen)*Mf(tt,gen)*MNeu(1) +
     -             	ER(1,i,tt,gen)*(4*xd00 + (xd1 + xd11 + 2*xd12 + 2*xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33)*Mf(tt,gen)**2 + xd11*MNeu(1)**2))) -
     -       		2*xd0*(ELc(1,i,tt,gen)*FL(1,i,tt,gen) + ER(1,i,tt,gen)*FR(1,i,tt,gen))*pbp2(tt,gen)))/(16.*pi)
 		   	endif

			enddo
		enddo
	enddo

!----------------------- B2 contribution
	do gen = 1,3
!	gen = 3
		do uptype = 1,2
			do i = 1,2
			do j = 1,2

			tt = uptype + 2

			call myDget(Mf(tt,gen)**2, t_Mandel(tt,gen), MNeu(1)**2, s_Mandel(tt,gen), Mf(tt,gen)**2, MNeu(1)**2, MGl**2, MSf(i,tt,gen)**2, MSf(j,tt,gen)**2, Mf(tt,gen)**2)

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) + CF*((-(ELc(1,j,tt,gen)*FL(1,i,tt,gen)*Mf(tt,gen)*(JR(i,tt,gen)*(-(MGl*xd0*IR(j,tt,gen)) + (xd1 + xd2 + xd3)*IL(j,tt,gen)*Mf(tt,gen)) +
     -          JL(i,tt,gen)*(-(MGl*xd0*IL(j,tt,gen)) + (xd1 + xd2 + xd3)*IR(j,tt,gen)*Mf(tt,gen)))) -
     -       	ER(1,j,tt,gen)*FR(1,i,tt,gen)*Mf(tt,gen)*(JR(i,tt,gen)*(-(MGl*xd0*IR(j,tt,gen)) + (xd1 + xd2 + xd3)*IL(j,tt,gen)*Mf(tt,gen)) +
     -          JL(i,tt,gen)*(-(MGl*xd0*IL(j,tt,gen)) + (xd1 + xd2 + xd3)*IR(j,tt,gen)*Mf(tt,gen))) + MGl*(xd0 + xd3)*ER(1,j,tt,gen)*FL(1,i,tt,gen)*IL(j,tt,gen)*JL(i,tt,gen)*MNeu(1) +
     -       	MGl*(xd0 + xd3)*ELc(1,j,tt,gen)*FR(1,i,tt,gen)*IL(j,tt,gen)*JL(i,tt,gen)*MNeu(1) + MGl*(xd0 + xd3)*ER(1,j,tt,gen)*FL(1,i,tt,gen)*IR(j,tt,gen)*JR(i,tt,gen)*MNeu(1) +
     -       	MGl*(xd0 + xd3)*ELc(1,j,tt,gen)*FR(1,i,tt,gen)*IR(j,tt,gen)*JR(i,tt,gen)*MNeu(1) -
     -       	(xd1 + xd13 + xd2 + xd23 + xd3 + xd33)*ER(1,j,tt,gen)*FL(1,i,tt,gen)*IR(j,tt,gen)*JL(i,tt,gen)*Mf(tt,gen)*MNeu(1) -
     -       	(xd1 + xd13 + xd2 + xd23 + xd3 + xd33)*ELc(1,j,tt,gen)*FR(1,i,tt,gen)*IR(j,tt,gen)*JL(i,tt,gen)*Mf(tt,gen)*MNeu(1) -
     -       	(xd1 + xd13 + xd2 + xd23 + xd3 + xd33)*ER(1,j,tt,gen)*FL(1,i,tt,gen)*IL(j,tt,gen)*JR(i,tt,gen)*Mf(tt,gen)*MNeu(1) -
     -       	(xd1 + xd13 + xd2 + xd23 + xd3 + xd33)*ELc(1,j,tt,gen)*FR(1,i,tt,gen)*IL(j,tt,gen)*JR(i,tt,gen)*Mf(tt,gen)*MNeu(1)))/(64.*pi**2)

			enddo
			enddo
		enddo
	enddo

!----------------------- B3 contribution
! Note: The old box mode (0) activates the old amplitude which is wrong in the limit v -> 0. This option is kept for testing purposes.

	do gen = 1,3
!	gen = 3
		do uptype = 1,2
			do i = 1,2
			tt = uptype + 2

			call myDget(u_Mandel(tt,gen), MNeu(1)**2, t_Mandel(tt,gen), Mf(tt,gen)**2, Mf(tt,gen)**2,MNeu(1)**2, 0d0, MSf(i,tt,gen)**2, Mf(tt,gen)**2, Mf(tt,gen)**2)

			if(Box_Mode.eq.0) then
				scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - (Alfas*CF*(-(Mf(tt,gen)*(ER(1,i,tt,gen)*((xd1 + xd11 + 2d0*xd12 + 2d0*xd13 + xd2 + xd22 + 2d0*xd23 + xd3 + xd33)*FR(1,i,tt,gen)*Mf(tt,gen) + 2d0*xd13*FL(1,i,tt,gen)*MNeu(1)) +
     -          	ELc(1,i,tt,gen)*((xd1 + xd11 + 2d0*xd12 + 2d0*xd13 + xd2 + xd22 + 2d0*xd23 + xd3 + xd33)*FL(1,i,tt,gen)*Mf(tt,gen) + 2d0*xd13*FR(1,i,tt,gen)*MNeu(1)))) -
     -       		2d0*(xd0 + xd1 + xd13 + xd2 + xd23 + xd3)*(ELc(1,i,tt,gen)*FL(1,i,tt,gen) + ER(1,i,tt,gen)*FR(1,i,tt,gen))*pbp2(tt,gen)) - 8d0*(FR(1,i,tt,gen)*ER(1,i,tt,gen) + FL(1,i,tt,gen)*ELc(1,i,tt,gen))*xd00)/(16.*pi)
	 		else
				scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - (Alfas*CF*(-2*(FL(1,i,tt,gen)*((-xd1 + xd11 + xd12 + xd13)*ER(1,i,tt,gen)*Mf(tt,gen)*MNeu(1) +
     -             	ELc(1,i,tt,gen)*(4*xd00 + (xd1 + xd11 + 2*xd12 + 2*xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33)*Mf(tt,gen)**2 + xd11*MNeu(1)**2)) +
     -          	FR(1,i,tt,gen)*((-xd1 + xd11 + xd12 + xd13)*ELc(1,i,tt,gen)*Mf(tt,gen)*MNeu(1) +
     -             	ER(1,i,tt,gen)*(4*xd00 + (xd1 + xd11 + 2*xd12 + 2*xd13 + xd2 + xd22 + 2*xd23 + xd3 + xd33)*Mf(tt,gen)**2 + xd11*MNeu(1)**2))) -
     -       		2*xd0*(ELc(1,i,tt,gen)*FL(1,i,tt,gen) + ER(1,i,tt,gen)*FR(1,i,tt,gen))*pbp2(tt,gen)))/(16.*pi)
 	   		endif

			enddo
		enddo
	enddo

!----------------------- B4 contribution
	do gen = 1,3
!	gen = 3
		do uptype = 1,2
!			uptype = 1
			do i = 1,2
			do j = 1,2
				
			tt = uptype + 2

			call myDget(Mf(tt,gen)**2, t_Mandel(tt,gen), MNeu(1)**2, u_Mandel(tt,gen), Mf(tt,gen)**2, MNeu(1)**2, MGl**2, MSf(i,tt,gen)**2, MSf(j,tt,gen)**2, Mf(tt,gen)**2)

C  			write(*,*)''
C  			write(*,*)'tt, i, j', tt, i, j
C  			write(*,*)''
C  			write(*,*)'D0', xd0
C  			write(*,*)'D1, D2, D3', xd1, xd2, xd3
C  			write(*,*)'D13, D23, D33', xd13, xd23, xd33

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) + CF*((-(ELc(1,j,tt,gen)*FL(1,i,tt,gen)*Mf(tt,gen)*(JR(i,tt,gen)*(-(MGl*xd0*IR(j,tt,gen)) + (xd1 + xd2 + xd3)*IL(j,tt,gen)*Mf(tt,gen)) +
     -          JL(i,tt,gen)*(-(MGl*xd0*IL(j,tt,gen)) + (xd1 + xd2 + xd3)*IR(j,tt,gen)*Mf(tt,gen)))) -
     -       	ER(1,j,tt,gen)*FR(1,i,tt,gen)*Mf(tt,gen)*(JR(i,tt,gen)*(-(MGl*xd0*IR(j,tt,gen)) + (xd1 + xd2 + xd3)*IL(j,tt,gen)*Mf(tt,gen)) +
     -          JL(i,tt,gen)*(-(MGl*xd0*IL(j,tt,gen)) + (xd1 + xd2 + xd3)*IR(j,tt,gen)*Mf(tt,gen))) + MGl*(xd0 + xd3)*ER(1,j,tt,gen)*FL(1,i,tt,gen)*IL(j,tt,gen)*JL(i,tt,gen)*MNeu(1) +
     -       	MGl*(xd0 + xd3)*ELc(1,j,tt,gen)*FR(1,i,tt,gen)*IL(j,tt,gen)*JL(i,tt,gen)*MNeu(1) + MGl*(xd0 + xd3)*ER(1,j,tt,gen)*FL(1,i,tt,gen)*IR(j,tt,gen)*JR(i,tt,gen)*MNeu(1) +
     -       	MGl*(xd0 + xd3)*ELc(1,j,tt,gen)*FR(1,i,tt,gen)*IR(j,tt,gen)*JR(i,tt,gen)*MNeu(1) -
     -       	(xd1 + xd13 + xd2 + xd23 + xd3 + xd33)*ER(1,j,tt,gen)*FL(1,i,tt,gen)*IR(j,tt,gen)*JL(i,tt,gen)*Mf(tt,gen)*MNeu(1) -
     -       	(xd1 + xd13 + xd2 + xd23 + xd3 + xd33)*ELc(1,j,tt,gen)*FR(1,i,tt,gen)*IR(j,tt,gen)*JL(i,tt,gen)*Mf(tt,gen)*MNeu(1) -
     -       	(xd1 + xd13 + xd2 + xd23 + xd3 + xd33)*ER(1,j,tt,gen)*FL(1,i,tt,gen)*IL(j,tt,gen)*JR(i,tt,gen)*Mf(tt,gen)*MNeu(1) -
     -       	(xd1 + xd13 + xd2 + xd23 + xd3 + xd33)*ELc(1,j,tt,gen)*FR(1,i,tt,gen)*IL(j,tt,gen)*JR(i,tt,gen)*Mf(tt,gen)*MNeu(1)))/(64.*pi**2)

			enddo
			enddo
		enddo
	enddo
	
!----------------------- Total unsuppressed scalar contribution
	do gen = 1,3
!	gen = 3	
		do uptype = 1,2
			do nucleon = 1,2
				tt = uptype + 2
				g_4Fs(nucleon) = g_4Fs(nucleon) + f_T(nucleon,uptype,gen)*m_Nucleon(nucleon)/Mf_low(tt,gen)*scalarfactor(uptype,gen)*Running_SI(tt,gen)
				g_4FsContribs(nucleon,6) = g_4FsContribs(nucleon,6) + f_T(nucleon,uptype,gen)*m_Nucleon(nucleon)/Mf_low(tt,gen)*scalarfactor(uptype,gen)*Running_SI(tt,gen)
			enddo
		enddo
	enddo

!----------------------- Unsuppressed axial contributions -----------------------------------------
	do gen = 1,3
		do uptype = 1,2
			axialfactor(uptype,gen) = 0.0d0
		enddo
	enddo

!Note: For simplicity we loop over all quark types in the following, but for spin dependent only u,d and s contribute. This is
!	   taken into account by the uds Check.

!----------------------- B1 contribution
! Note: The old box mode (0) activates the old amplitude which is wrong in the limit v -> 0. This option is kept for testing purposes.

	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
			tt = uptype+2

			if(DD_udsCheck(tt,gen)) then

			call myDget(s_Mandel(tt,gen),MNeu(1)**2,t_Mandel(tt,gen),Mf(tt,gen)**2,Mf(tt,gen)**2,MNeu(1)**2,0d0,MSf(i,tt,gen)**2,Mf(tt,gen)**2,Mf(tt,gen)**2)

			if(Box_Mode.eq.0) then
				axialfactor(uptype,gen) = axialfactor(uptype,gen) - (Alfas*CF*(Mf(tt,gen)*(FL(1,i,tt,gen)*((xd1 + xd11 + xd12 + 2d0*xd13 + xd2 + xd22 + xd23 + xd3 + xd33)*ER(1,i,tt,gen)*Mf(tt,gen) + xd12*ELc(1,i,tt,gen)*MNeu(1)) +
     -      		FR(1,i,tt,gen)*((xd1 + xd11 + xd12 + 2d0*xd13 + xd2 + xd22 + xd23 + xd3 + xd33)*ELc(1,i,tt,gen)*Mf(tt,gen) + xd12*ER(1,i,tt,gen)*MNeu(1))) +
     -        		2d0*(xd0 + xd1 + xd12 + xd2 + xd23 + xd3)*(ER(1,i,tt,gen)*FL(1,i,tt,gen) + ELc(1,i,tt,gen)*FR(1,i,tt,gen))*pbp2(tt,gen)) + 2d0*(FL(1,i,tt,gen)*ER(1,i,tt,gen) + FR(1,i,tt,gen)*ELc(1,i,tt,gen))*xd00)/(16.*pi)
 			else
				axialfactor(uptype,gen) = axialfactor(uptype,gen) - (Alfas*CF*(FR(1,i,tt,gen)*(-2*xd1*ER(1,i,tt,gen)*Mf(tt,gen)*MNeu(1) +
     -           	ELc(1,i,tt,gen)*(2*xd00 + (2*xd1 + xd11 + 2*xd12 + 2*xd13 + 2*xd2 + xd22 + 2*xd23 + 2*xd3 + xd33)*Mf(tt,gen)**2 + xd11*MNeu(1)**2)) +
     -        		FL(1,i,tt,gen)*(-2*xd1*ELc(1,i,tt,gen)*Mf(tt,gen)*MNeu(1) +
     -           	ER(1,i,tt,gen)*(2*xd00 + (2*xd1 + xd11 + 2*xd12 + 2*xd13 + 2*xd2 + xd22 + 2*xd23 + 2*xd3 + xd33)*Mf(tt,gen)**2 + xd11*MNeu(1)**2)) +
     -        		2*xd0*(ER(1,i,tt,gen)*FL(1,i,tt,gen) + ELc(1,i,tt,gen)*FR(1,i,tt,gen))*pbp2(tt,gen)))/(16.*pi)
 		   	endif

			endif

			enddo
		enddo
	enddo

!----------------------- B2 contribution
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
				do j = 1,2
				tt = uptype+2

				if(DD_udsCheck(tt,gen)) then

				call myDget(Mf(tt,gen)**2, t_Mandel(tt,gen), MNeu(1)**2, s_Mandel(tt,gen), Mf(tt,gen)**2, MNeu(1)**2, MGl**2, MSf(i,tt,gen)**2, MSf(j,tt,gen)**2, Mf(tt,gen)**2)

				axialfactor(uptype,gen) = axialfactor(uptype,gen) + CF*((-(xd00*ER(1,j,tt,gen)*FL(1,i,tt,gen)*IR(j,tt,gen)*JL(i,tt,gen)) + xd00*ELc(1,j,tt,gen)*FR(1,i,tt,gen)*IR(j,tt,gen)*JL(i,tt,gen) +
     -				xd00*ER(1,j,tt,gen)*FL(1,i,tt,gen)*IL(j,tt,gen)*JR(i,tt,gen) - xd00*ELc(1,j,tt,gen)*FR(1,i,tt,gen)*IL(j,tt,gen)*JR(i,tt,gen)))/(64.*pi**2)

				endif

				enddo
			enddo
		enddo
	enddo

!----------------------- B3 contribution
! Note: The old box mode (0) activates the old amplitude which is wrong in the limit v -> 0. This option is kept for testing purposes.

	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
			tt = uptype+2

			if(DD_udsCheck(tt,gen)) then

			call myDget(u_Mandel(tt,gen), MNeu(1)**2, t_Mandel(tt,gen), Mf(tt,gen)**2, Mf(tt,gen)**2,MNeu(1)**2, 0d0, MSf(i,tt,gen)**2, Mf(tt,gen)**2, Mf(tt,gen)**2)

			if(Box_Mode.eq.0) then
				axialfactor(uptype,gen) = axialfactor(uptype,gen) - (Alfas*CF*(Mf(tt,gen)*(ELc(1,i,tt,gen)*((xd1 + xd11 + 2d0*xd12 + xd13 + xd2 + xd22 + xd23 + xd3 + xd33)*FR(1,i,tt,gen)*Mf(tt,gen) + xd13*FL(1,i,tt,gen)*MNeu(1)) +
     -      		ER(1,i,tt,gen)*((xd1 + xd11 + 2d0*xd12 + xd13 + xd2 + xd22 + xd23 + xd3 + xd33)*FL(1,i,tt,gen)*Mf(tt,gen) + xd13*FR(1,i,tt,gen)*MNeu(1))) +
     -        		2d0*(xd0 + xd1 + xd13 + xd2 + xd23 + xd3)*(ER(1,i,tt,gen)*FL(1,i,tt,gen) + ELc(1,i,tt,gen)*FR(1,i,tt,gen))*pbp2(tt,gen)) + 2d0*(FL(1,i,tt,gen)*ER(1,i,tt,gen) + FR(1,i,tt,gen)*ELc(1,i,tt,gen))*xd00)/(16.*pi)
 			else
				axialfactor(uptype,gen) = axialfactor(uptype,gen) - (Alfas*CF*(FR(1,i,tt,gen)*(-2*xd1*ER(1,i,tt,gen)*Mf(tt,gen)*MNeu(1) +
     -           	ELc(1,i,tt,gen)*(2*xd00 + (2*xd1 + xd11 + 2*xd12 + 2*xd13 + 2*xd2 + xd22 + 2*xd23 + 2*xd3 + xd33)*Mf(tt,gen)**2 + xd11*MNeu(1)**2)) +
     -        		FL(1,i,tt,gen)*(-2*xd1*ELc(1,i,tt,gen)*Mf(tt,gen)*MNeu(1) +
     -           	ER(1,i,tt,gen)*(2*xd00 + (2*xd1 + xd11 + 2*xd12 + 2*xd13 + 2*xd2 + xd22 + 2*xd23 + 2*xd3 + xd33)*Mf(tt,gen)**2 + xd11*MNeu(1)**2)) +
     -        		2*xd0*(ER(1,i,tt,gen)*FL(1,i,tt,gen) + ELc(1,i,tt,gen)*FR(1,i,tt,gen))*pbp2(tt,gen)))/(16.*pi)
 		   	endif

			endif

			enddo
		enddo
	enddo

!----------------------- B4 contribution
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
				do j = 1,2
				tt = uptype+2

				if(DD_udsCheck(tt,gen)) then

				call myDget(Mf(tt,gen)**2, t_Mandel(tt,gen), MNeu(1)**2, u_Mandel(tt,gen), Mf(tt,gen)**2, MNeu(1)**2, MGl**2, MSf(i,tt,gen)**2, MSf(j,tt,gen)**2, Mf(tt,gen)**2)

				axialfactor(uptype,gen) = axialfactor(uptype,gen) + CF*((-(xd00*ER(1,j,tt,gen)*FL(1,i,tt,gen)*IR(j,tt,gen)*JL(i,tt,gen)) + xd00*ELc(1,j,tt,gen)*FR(1,i,tt,gen)*IR(j,tt,gen)*JL(i,tt,gen) +
     -				xd00*ER(1,j,tt,gen)*FL(1,i,tt,gen)*IL(j,tt,gen)*JR(i,tt,gen) - xd00*ELc(1,j,tt,gen)*FR(1,i,tt,gen)*IL(j,tt,gen)*JR(i,tt,gen)))/(64.*pi**2)

				endif

				enddo
			enddo
		enddo
	enddo

!----------------------- Total unsuppressed axial contribution
	do gen = 1,3
!	gen = 3
		do uptype = 1,2
			do nucleon = 1,2
				g_4Fa(nucleon) = g_4Fa(nucleon) + Delta_q(nucleon,uptype,gen)*axialfactor(uptype,gen)*Running_SD(tt,gen)
				g_4FaContribs(nucleon,6) = g_4FaContribs(nucleon,6) + Delta_q(nucleon,uptype,gen)*axialfactor(uptype,gen)*Running_SD(tt,gen)
			enddo
		enddo
	enddo
	
! ====================== SUPPRESSED CONTRIBUTIONS =================================================

	if(Supp_Mode==1.OR.Supp_Mode==2) then

!----------------------- Suppressed scalar contributions ------------------------------------------

!----------------------- Total suppressed scalar contribution

!----------------------- Suppressed axial contributions -------------------------------------------
		
!----------------------- Total suppressed axial contribution

	endif				! End of Supp_Mode if
	
! ====================== STRONGLY SUPPRESSED CONTRIBUTIONS ========================================

	if(Supp_Mode==2) then

!----------------------- Strongly suppressed scalar contributions ---------------------------------

!----------------------- Total strongly suppressed scalar contribution

!----------------------- Strongly suppressed axial contributions ----------------------------------
		
!----------------------- Total strongly suppressed axial contribution

	endif				! End of Supp_Mode if
	
!======================= END OF PROGRAM ===========================================================
	
	g_4Fsp = g_4Fs(1)
	g_4Fsn = g_4Fs(2)
	g_4Fap = g_4Fa(1)
	g_4Fan = g_4Fa(2)
	
	end


