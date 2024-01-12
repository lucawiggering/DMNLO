!==================================================================================================
!======================= Calculates vertex contributions to four fermion couplings ================
!=======================				Author: Patrick Steppeler	   	           ================
!==================================================================================================

	subroutine DD_Vertex
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_Model.h"
#include "DD_GenCouplings.h"
#include "../util/DMNLO_LoopIntegrals.h"

!----------------------- Declare variables --------------------------------------------------------
	complex*16 ALV1, ARV1, ALV2, ARV2, ALV5, ARV5, ALV6, ARV6, ALV7, ARV7, ALV8, ARV8
	complex*16 AL1V3, AR1V3, AL2V3, AR2V3, AL3V3, AR3V3 
	complex*16 AL1V4, AR1V4, AL2V4, AR2V4, AL3V4, AR3V4
	complex*16 axialfactor(2,3), scalarfactor(2,3)
	integer nucleon, uptype, gen, tt, i, j, a
	logical DD_udsCheck

!----------------------- Declare external subroutines ---------------------------------------------
	external myBget, myCget
	
!======================= ACTUAL CALCULATIONS ======================================================
	
!======================= UNSUPPRESSED CONTRIBUTIONS ===============================================

!----------------------- Unsuppressed scalar contributions ----------------------------------------
	do gen = 1,3
		do uptype = 1,2
			scalarfactor(uptype,gen) = 0.0d0
		enddo
	enddo
	
!----------------------- Higgs contribution (h0 & H0) V1: Gluon Vertex, t-Channel
! Note: This contribution is completely eliminated by the matching.
!
!	do gen = 1,3
!!	gen = 3
!		do uptype = 1,2
!			tt = uptype + 2
!
!			call myCget(Mf(tt,gen)**2, t_Mandel(tt,gen),Mf(tt,gen)**2, 0.0d0, Mf(tt,gen)**2, Mf(tt,gen)**2)
!			call myBget(t_Mandel(tt,gen), Mf(tt,gen)**2, Mf(tt,gen)**2)
!					
!			do i = 1,2
!				
!				ALV1 = CF*Alfas/(4d0*pi)*(-2.0d0*(xc1 + xc2)*DR(i,tt,gen)*Mf(tt,gen)**2 + 
!     &					DL(i,tt,gen)*(-2.0d0 + 4.0d0*xb0 + 4.0d0*pbp2(tt,gen)*(xc0 + xc1 + xc2) + 
!     &					2.0d0*(xc1 + xc2)*Mf(tt,gen)**2))
!
!				ARV1 = CF*Alfas/(4d0*pi)*(-2.0d0*(xc1 + xc2)*DL(i,tt,gen)*Mf(tt,gen)**2 + 
!     &					DR(i,tt,gen)*(-2.0d0 + 4.0d0*xb0 + 4.0d0*pbp2(tt,gen)*(xc0 + xc1 + xc2) + 
!     &					2.0d0*(xc1 + xc2)*Mf(tt,gen)**2))
!	
!				if (Denom_Mode==0) then
!				
!				scalarfactor(uptype,gen) = scalarfactor(uptype,gen) + CR(i,tt,gen)
!     &										*(ALV1 + ARV1)/(2.0d0*MHk(i)**2)
!				
!				else
!				
!				scalarfactor(uptype,gen) = scalarfactor(uptype,gen) + CR(i,tt,gen)
!     &										*(ALV1 + ARV1)/(2.0d0*(MHk(i)**2-t_Mandel(tt,gen)))
!				
!				endif
!
!			enddo
!		enddo
!	enddo

!----------------------- Higgs contribution (h0 & H0) V2: Gluino Vertex, t-Channel
	do gen = 1,3
!	gen = 3
		do uptype = 1,2
			tt = uptype + 2

			do i = 1,2
			do j = 1,2

			call myCget(Mf(tt,gen)**2,t_Mandel(tt,gen),Mf(tt,gen)**2,MGl**2,MSf(i,tt,gen)**2,MSf(j,tt,gen)**2)

			do a = 1,2

			ALV2 = -CF/(4d0*pi)**2*(HSfSf(a,i,j,tt,gen)*(MGl*xc0*IL(i,tt,gen)*JL(j,tt,gen) -
     &			xc1*IR(i,tt,gen)*JL(j,tt,gen)*Mf(tt,gen) -
     &			xc2*IL(i,tt,gen)*JR(j,tt,gen)*Mf(tt,gen)))

			ARV2 = -CF/(4d0*pi)**2*(HSfSf(a,i,j,tt,gen)*(MGl*xc0*IR(i,tt,gen)*JR(j,tt,gen) -
     &			xc2*IR(i,tt,gen)*JL(j,tt,gen)*Mf(tt,gen) -
     &			xc1*IL(i,tt,gen)*JR(j,tt,gen)*Mf(tt,gen)))

			if (Denom_Mode==0) then

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) + CR(a,tt,gen)
     &			*(ALV2 + ARV2)/(2.0d0*MHk(a)**2)

			else

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) + CR(a,tt,gen)
     &			*(ALV2 + ARV2)/(2.0d0*(MHk(a)**2-t_Mandel(tt,gen)))

			endif

			enddo
			enddo
			enddo
		enddo

	enddo
	
!----------------------- Squark contribution V5: Gluon Vertex, Outgoing Squark, s-Channel
	do gen = 1,3
!	gen = 3
		do uptype = 1,2
			tt = uptype + 2

			do i = 1,2

			call myCget(Mneu(1)**2, s_Mandel(tt,gen), Mf(tt,gen)**2, Mf(tt,gen)**2, MSf(i,tt,gen)**2, 0d0)
			call myBget(s_Mandel(tt,gen), MSf(i, tt, gen)**2, 0d0)

			ALV5 = -(CF*Alfas/(4d0*pi)*(2.0d0*(xc0 + xc1 + xc2)*FR(1,i,tt,gen)*Mf(tt,gen)*MNeu(1)
     &				- FL(1,i,tt,gen)*(xb0 - 4.0d0*papb(tt,gen)*xc2 + 2.0d0*(xc0 - xc2)*Mf(tt,gen)**2
     &				+ 2.0d0*xc1*MNeu(1)**2)))

			ARV5 = -(CF*Alfas/(4d0*pi)*(2.0d0*(xc0 + xc1 + xc2)*FL(1,i,tt,gen)*Mf(tt,gen)*MNeu(1)
     &				- FR(1,i,tt,gen)*(xb0 - 4.0d0*papb(tt,gen)*xc2 + 2.0d0*(xc0 - xc2)*Mf(tt,gen)**2
     &				+ 2.0d0*xc1*MNeu(1)**2)))

			if (Denom_Mode==0) then

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV5 + ALV5)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) +
     &			(ARV5 - ALV5)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/MSf(i,tt,gen)**2

			else

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV5 + ALV5)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) +
     &			(ARV5 - ALV5)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/(MSf(i,tt,gen)**2-s_Mandel(tt,gen))

			endif

			enddo
		enddo
	enddo

!----------------------- Squark contribution V6: Gluino Vertex, Outgoing Squark, s-Channel
	do gen = 1,3
!	gen = 3
		do uptype = 1,2
			tt = uptype + 2

			do i = 1,2
			do j = 1,2

			call myCget(Mf(tt,gen)**2,s_Mandel(tt,gen),MNeu(1)**2,MSf(j,tt,gen)**2,MGl**2,Mf(tt,gen)**2)
			call myBget(s_Mandel(tt,gen), MGl**2, Mf(tt,gen)**2)

			ALV6 = CF/(4d0*pi)**2*(-(ER(1,j,tt,gen)*(xc2*IL(j,tt,gen)*IR(i,tt,gen)*Mf(tt,gen) +
     &				IL(i,tt,gen)*(MGl*(xc0 + xc2)*IL(j,tt,gen) + (xc0 + xc1 + xc2)*IR(j,tt,gen)*
     &				Mf(tt,gen)))*MNeu(1)) - ELc(1,j,tt,gen)*(IR(j,tt,gen)*Mf(tt,gen)*(MGl*xc1*
     &				IR(i,tt,gen) + (xc0 + xc1)*IL(i,tt,gen)*Mf(tt,gen)) + IL(j,tt,gen)*(MGl*xc0*
     &				IL(i,tt,gen)*Mf(tt,gen) + IR(i,tt,gen)*(xb0 + xc1*Mf(tt,gen)**2 + xc2*
     &				MNeu(1)**2 + xc0*MSf(j,tt,gen)**2))))

			ARV6 = CF/(4d0*pi)**2*(-(ELc(1,j,tt,gen)*(xc2*IL(i,tt,gen)*IR(j,tt,gen)*Mf(tt,gen) +
     &				IR(i,tt,gen)*(MGl*(xc0 + xc2)*IR(j,tt,gen) + (xc0 + xc1 + xc2)*IL(j,tt,gen)*
     &				Mf(tt,gen)))*MNeu(1)) - ER(1,j,tt,gen)*(IR(i,tt,gen)*Mf(tt,gen)*(MGl*xc0*
     &				IR(j,tt,gen) + (xc0 + xc1)*IL(j,tt,gen)*Mf(tt,gen)) + IL(i,tt,gen)*(MGl*xc1*
     &				IL(j,tt,gen)*Mf(tt,gen) + IR(j,tt,gen)*(xb0 + xc1*Mf(tt,gen)**2 + xc2*
     &				MNeu(1)**2 + xc0*MSf(j,tt,gen)**2))))

			if (Denom_Mode==0) then

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV6 + ALV6)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) +
     &			(ARV6 - ALV6)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/MSf(i,tt,gen)**2

			else

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV6 + ALV6)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) +
     &			(ARV6 - ALV6)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/(MSf(i,tt,gen)**2-s_Mandel(tt,gen))

			endif

			enddo
			enddo
		enddo
	enddo

!----------------------- Squark contribution V7: Gluon Vertex, Ingoing Squark, s-Channel
	do gen = 1,3
!	gen = 3
		do uptype = 1,2
			tt = uptype + 2

			do i = 1,2

			call myCget(Mneu(1)**2, s_Mandel(tt,gen), Mf(tt,gen)**2, Mf(tt,gen)**2, MSf(i,tt,gen)**2, 0d0)
			call myBget(s_Mandel(tt,gen), MSf(i, tt, gen)**2, 0d0)

			ALV7 = -(CF*Alfas/(4d0*pi)*(2.0d0*(xc0 + xc1 + xc2)*ER(1,i,tt,gen)*Mf(tt,gen)*MNeu(1)
     &				- ELc(1,i,tt,gen)*(xb0 - 4.0d0*p1p2(tt,gen)*xc2 + 2.0d0*(xc0 - xc2)*Mf(tt,gen)**2
     &				+ 2.0d0*xc1*MNeu(1)**2)))

			ARV7 = -(CF*Alfas/(4d0*pi)*(2.0d0*(xc0 + xc1 + xc2)*ELc(1,i,tt,gen)*Mf(tt,gen)*MNeu(1)
     &				- ER(1,i,tt,gen)*(xb0 - 4.0d0*p1p2(tt,gen)*xc2 + 2.0d0*(xc0 - xc2)*Mf(tt,gen)**2
     &				+ 2.0d0*xc1*MNeu(1)**2)))

			if (Denom_Mode==0) then

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV7 + ALV7)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) +
     &			(ARV7 - ALV7)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/MSf(i,tt,gen)**2

			else

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV7 + ALV7)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) +
     &			(ARV7 - ALV7)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/(MSf(i,tt,gen)**2-s_Mandel(tt,gen))

			endif

			enddo
		enddo
	enddo

!----------------------- Squark contribution V8: Gluino Vertex, Ingoing Squark, s-Channel
	do gen = 1,3
!	gen = 3
		do uptype = 1,2
			tt = uptype + 2

			do i = 1,2
			do j = 1,2

			call myCget(MNeu(1)**2,s_Mandel(tt,gen),Mf(tt,gen)**2,MSf(j,tt,gen)**2,Mf(tt,gen)**2,MGl**2)
			call myBget(s_Mandel(tt,gen), Mf(tt,gen)**2, MGl**2)

			ALV8 = CF/(4d0*pi)**2*(-(FR(1,j,tt,gen)*(xc1*JL(j,tt,gen)*JR(i,tt,gen)*Mf(tt,gen) +
     &				JL(i,tt,gen)*(MGl*(xc0 + xc1)*JL(j,tt,gen) + (xc0 + xc1 + xc2)*JR(j,tt,gen)*
     &				Mf(tt,gen)))*MNeu(1)) - FL(1,j,tt,gen)*(JR(j,tt,gen)*Mf(tt,gen)*(MGl*xc2*
     &				JR(i,tt,gen) + (xc0 + xc2)*JL(i,tt,gen)*Mf(tt,gen)) + JL(j,tt,gen)*(MGl*xc0*
     &				JL(i,tt,gen)*Mf(tt,gen) + JR(i,tt,gen)*(xb0 + xc2*Mf(tt,gen)**2 + xc1*
     &				MNeu(1)**2 + xc0*MSf(j,tt,gen)**2))))

			ARV8 = CF/(4d0*pi)**2*(-(FL(1,j,tt,gen)*(xc1*JL(i,tt,gen)*JR(j,tt,gen)*Mf(tt,gen) +
     &				JR(i,tt,gen)*(MGl*(xc0 + xc1)*JR(j,tt,gen) + (xc0 + xc1 + xc2)*JL(j,tt,gen)*
     &				Mf(tt,gen)))*MNeu(1)) - FR(1,j,tt,gen)*(JR(i,tt,gen)*Mf(tt,gen)*(MGl*xc0*
     &				JR(j,tt,gen) + (xc0 + xc2)*JL(j,tt,gen)*Mf(tt,gen)) + JL(i,tt,gen)*(MGl*xc2*
     &				JL(j,tt,gen)*Mf(tt,gen) + JR(j,tt,gen)*(xb0 + xc2*Mf(tt,gen)**2 + xc1*
     &				MNeu(1)**2 + xc0*MSf(j,tt,gen)**2))))

			if (Denom_Mode==0) then

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV8 + ALV8)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) +
     &			(ARV8 - ALV8)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/MSf(i,tt,gen)**2

			else

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV8 + ALV8)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) +
     &			(ARV8 - ALV8)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/(MSf(i,tt,gen)**2-s_Mandel(tt,gen))

			endif

			enddo
			enddo
		enddo
	enddo

!----------------------- Squark contribution V5: Gluon Vertex, Outgoing Squark, u-Channel
! NOTE: In case of the u-Channels we have s <-> u because pb <-> - p2, so the calculations have the same form.
! Effectively s <-> u, papb = p1p2 <-> -pap2 = -p1pb

	do gen = 1,3
!	gen = 3
		do uptype = 1,2
			tt = uptype + 2

			do i = 1,2

			call myCget(Mneu(1)**2, u_Mandel(tt,gen), Mf(tt,gen)**2, Mf(tt,gen)**2, MSf(i,tt,gen)**2, 0d0)
			call myBget(u_Mandel(tt,gen), MSf(i, tt, gen)**2, 0d0)

			ALV5 = -(CF*Alfas/(4d0*pi)*(2.0d0*(xc0 + xc1 + xc2)*FR(1,i,tt,gen)*Mf(tt,gen)*MNeu(1)
     &				- FL(1,i,tt,gen)*(xb0 + 4.0d0*pap2(tt,gen)*xc2 + 2.0d0*(xc0 - xc2)*Mf(tt,gen)**2
     &				+ 2.0d0*xc1*MNeu(1)**2)))

			ARV5 = -(CF*Alfas/(4d0*pi)*(2.0d0*(xc0 + xc1 + xc2)*FL(1,i,tt,gen)*Mf(tt,gen)*MNeu(1)
     &				- FR(1,i,tt,gen)*(xb0 + 4.0d0*pap2(tt,gen)*xc2 + 2.0d0*(xc0 - xc2)*Mf(tt,gen)**2
     &				+ 2.0d0*xc1*MNeu(1)**2)))

			if (Denom_Mode==0) then

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV5 + ALV5)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) +
     &			(ARV5 - ALV5)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/MSf(i,tt,gen)**2

			else

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV5 + ALV5)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) +
     &			(ARV5 - ALV5)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/(MSf(i,tt,gen)**2-u_Mandel(tt,gen))

			endif

			enddo
		enddo
	enddo

!----------------------- Squark contribution V6: Gluino Vertex, Outgoing Squark, u-Channel
	do gen = 1,3
!	gen = 3
		do uptype = 1,2
			tt = uptype + 2

			do i = 1,2
			do j = 1,2

			call myCget(Mf(tt,gen)**2,u_Mandel(tt,gen),MNeu(1)**2,MSf(j,tt,gen)**2,MGl**2,Mf(tt,gen)**2)
			call myBget(u_Mandel(tt,gen), MGl**2, Mf(tt,gen)**2)

			ALV6 = CF/(4d0*pi)**2*(-(ER(1,j,tt,gen)*(xc2*IL(j,tt,gen)*IR(i,tt,gen)*Mf(tt,gen) +
     &				IL(i,tt,gen)*(MGl*(xc0 + xc2)*IL(j,tt,gen) + (xc0 + xc1 + xc2)*IR(j,tt,gen)*
     &				Mf(tt,gen)))*MNeu(1)) - ELc(1,j,tt,gen)*(IR(j,tt,gen)*Mf(tt,gen)*(MGl*xc1*
     &				IR(i,tt,gen) + (xc0 + xc1)*IL(i,tt,gen)*Mf(tt,gen)) + IL(j,tt,gen)*(MGl*xc0*
     &				IL(i,tt,gen)*Mf(tt,gen) + IR(i,tt,gen)*(xb0 + xc1*Mf(tt,gen)**2 + xc2*
     &				MNeu(1)**2 + xc0*MSf(j,tt,gen)**2))))

			ARV6 = CF/(4d0*pi)**2*(-(ELc(1,j,tt,gen)*(xc2*IL(i,tt,gen)*IR(j,tt,gen)*Mf(tt,gen) +
     &				IR(i,tt,gen)*(MGl*(xc0 + xc2)*IR(j,tt,gen) + (xc0 + xc1 + xc2)*IL(j,tt,gen)*
     &				Mf(tt,gen)))*MNeu(1)) - ER(1,j,tt,gen)*(IR(i,tt,gen)*Mf(tt,gen)*(MGl*xc0*
     &				IR(j,tt,gen) + (xc0 + xc1)*IL(j,tt,gen)*Mf(tt,gen)) + IL(i,tt,gen)*(MGl*xc1*
     &				IL(j,tt,gen)*Mf(tt,gen) + IR(j,tt,gen)*(xb0 + xc1*Mf(tt,gen)**2 + xc2*
     &				MNeu(1)**2 + xc0*MSf(j,tt,gen)**2))))

			if (Denom_Mode==0) then

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV6 + ALV6)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) +
     &			(ARV6 - ALV6)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/MSf(i,tt,gen)**2

			else

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV6 + ALV6)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) +
     &			(ARV6 - ALV6)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/(MSf(i,tt,gen)**2-u_Mandel(tt,gen))

			endif

			enddo
			enddo
		enddo
	enddo

!----------------------- Squark contribution V7: Gluon Vertex, Ingoing Squark, u-Channel
	do gen = 1,3
!	gen = 3
		do uptype = 1,2
			tt = uptype + 2

			do i = 1,2

			call myCget(Mneu(1)**2, u_Mandel(tt,gen), Mf(tt,gen)**2, Mf(tt,gen)**2, MSf(i,tt,gen)**2, 0d0)
			call myBget(u_Mandel(tt,gen), MSf(i, tt, gen)**2, 0d0)

			ALV7 = -(CF*Alfas/(4d0*pi)*(2.0d0*(xc0 + xc1 + xc2)*ER(1,i,tt,gen)*Mf(tt,gen)*MNeu(1)
     &				- ELc(1,i,tt,gen)*(xb0 + 4.0d0*pap2(tt,gen)*xc2 + 2.0d0*(xc0 - xc2)*Mf(tt,gen)**2
     &				+ 2.0d0*xc1*MNeu(1)**2)))

			ARV7 = -(CF*Alfas/(4d0*pi)*(2.0d0*(xc0 + xc1 + xc2)*ELc(1,i,tt,gen)*Mf(tt,gen)*MNeu(1)
     &				- ER(1,i,tt,gen)*(xb0 + 4.0d0*pap2(tt,gen)*xc2 + 2.0d0*(xc0 - xc2)*Mf(tt,gen)**2
     &				+ 2.0d0*xc1*MNeu(1)**2)))

			if (Denom_Mode==0) then

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV7 + ALV7)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) +
     &			(ARV7 - ALV7)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/MSf(i,tt,gen)**2

			else

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV7 + ALV7)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) +
     &			(ARV7 - ALV7)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/(MSf(i,tt,gen)**2-u_Mandel(tt,gen))

			endif

			enddo
		enddo
	enddo

!----------------------- Squark contribution V8: Gluino Vertex, Ingoing Squark, u-Channel
	do gen = 1,3
!	gen = 3
		do uptype = 1,2
			tt = uptype + 2

			do i = 1,2
			do j = 1,2

			call myCget(MNeu(1)**2,u_Mandel(tt,gen),Mf(tt,gen)**2,MSf(j,tt,gen)**2,Mf(tt,gen)**2,MGl**2)
			call myBget(u_Mandel(tt,gen), Mf(tt,gen)**2, MGl**2)

			ALV8 = CF/(4d0*pi)**2*(-(FR(1,j,tt,gen)*(xc1*JL(j,tt,gen)*JR(i,tt,gen)*Mf(tt,gen) +
     &				JL(i,tt,gen)*(MGl*(xc0 + xc1)*JL(j,tt,gen) + (xc0 + xc1 + xc2)*JR(j,tt,gen)*
     &				Mf(tt,gen)))*MNeu(1)) - FL(1,j,tt,gen)*(JR(j,tt,gen)*Mf(tt,gen)*(MGl*xc2*
     &				JR(i,tt,gen) + (xc0 + xc2)*JL(i,tt,gen)*Mf(tt,gen)) + JL(j,tt,gen)*(MGl*xc0*
     &				JL(i,tt,gen)*Mf(tt,gen) + JR(i,tt,gen)*(xb0 + xc2*Mf(tt,gen)**2 + xc1*
     &				MNeu(1)**2 + xc0*MSf(j,tt,gen)**2))))

			ARV8 = CF/(4d0*pi)**2*(-(FL(1,j,tt,gen)*(xc1*JL(i,tt,gen)*JR(j,tt,gen)*Mf(tt,gen) +
     &				JR(i,tt,gen)*(MGl*(xc0 + xc1)*JR(j,tt,gen) + (xc0 + xc1 + xc2)*JL(j,tt,gen)*
     &				Mf(tt,gen)))*MNeu(1)) - FR(1,j,tt,gen)*(JR(i,tt,gen)*Mf(tt,gen)*(MGl*xc0*
     &				JR(j,tt,gen) + (xc0 + xc2)*JL(j,tt,gen)*Mf(tt,gen)) + JL(i,tt,gen)*(MGl*xc2*
     &				JL(j,tt,gen)*Mf(tt,gen) + JR(j,tt,gen)*(xb0 + xc2*Mf(tt,gen)**2 + xc1*
     &				MNeu(1)**2 + xc0*MSf(j,tt,gen)**2))))

			if (Denom_Mode==0) then

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV8 + ALV8)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) +
     &			(ARV8 - ALV8)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/MSf(i,tt,gen)**2

			else

			scalarfactor(uptype,gen) = scalarfactor(uptype,gen) - 1.0d0/16.0d0*
     &			((ARV8 + ALV8)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) +
     &			(ARV8 - ALV8)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/(MSf(i,tt,gen)**2-u_Mandel(tt,gen))

			endif

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
				g_4Fs(nucleon) = g_4Fs(nucleon) + f_T(nucleon,uptype,gen)*m_Nucleon(nucleon)
     &							 /Mf_low(tt,gen)*scalarfactor(uptype,gen)*Running_SI(tt,gen)
				g_4FsContribs(nucleon,2) = g_4FsContribs(nucleon,2) + f_T(nucleon,uptype,gen)*m_Nucleon(nucleon)
     &							 /Mf_low(tt,gen)*scalarfactor(uptype,gen)*Running_SI(tt,gen)

			enddo
		enddo
	enddo
	
!----------------------- Unsuppressed axial contributions -----------------------------------------
	do gen = 1,3
		do uptype = 1,2
			axialfactor(uptype,gen) = 0.0d0
		enddo
	enddo

! Note: For simplicity we loop over all quark types in the following, but for spin dependent only u,d and s contribute.
! This is taken into account by the uds Check.
	
!----------------------- Z contribution V3: Gluon Vertex, t-Channel
! Note: This contribution is completely eliminated by the matching.

!	do gen = 1, 3
!		do uptype = 1,2
!		tt = uptype+2
!		
!		if (DD_udsCheck(tt,gen)) then
!		
!		call myCget(Mf(tt,gen)**2, t_Mandel(tt,gen), Mf(tt,gen)**2, 0.0d0, Mf(tt,gen)**2, Mf(tt,gen)**2) 
!		call myBget(t_Mandel(tt,gen), Mf(tt,gen)**2, Mf(tt,gen)**2)
!	
!		AL1V3 = CF*Alfas/(4d0*pi)*(2.0d0*(xc1 + xc2)*BR(tt,gen)*Mf(tt,gen)**2 + BL(tt,gen)*
!     &		(-2.0d0 + 2.0d0*xb0 + 4.0d0*pbp2(tt,gen)*xc0 + 4.0d0*xc00 - 8.0d0*xc00 + 1.0d0 + 4.0d0*
!     &		pbp2(tt,gen)*xc1 + 4.0d0*pbp2(tt,gen)*xc2 + 2.0d0*(xc1 + xc2)*Mf(tt,gen)**2))
!
!		AR1V3 = CF*Alfas/(4d0*pi)*(2.0d0*(xc1 + xc2)*BL(tt,gen)*Mf(tt,gen)**2 + BR(tt,gen)*
!     &		(-2.0d0 + 2.0d0*xb0 + 4.0d0*pbp2(tt,gen)*xc0 + 4.0d0*xc00 - 8.0d0*xc00 + 1.0d0 + 4.0d0*
!     &		pbp2(tt,gen)*xc1 + 4.0d0*pbp2(tt,gen)*xc2 + 2.0d0*(xc1 + xc2)*Mf(tt,gen)**2))
!	
!		if (Denom_Mode==0) then
!	
!			axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.5d0*AR(tt,gen)*(AL1V3-AR1V3)/MZ2
!	
!		else
!	
!			axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.5d0*AR(tt,gen)*(AL1V3-AR1V3)/(MZ2-t_Mandel(tt,gen))
!	
!		endif
!		
!		endif
!		
!		enddo
!	enddo
			
!----------------------- Z contribution V4: Gluino Vertex, t-Channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
				do j = 1,2
				tt = uptype+2

				if(DD_udsCheck(tt,gen)) then

				call myCget(Mf(tt,gen)**2, t_Mandel(tt,gen), Mf(tt,gen)**2, MGl**2, MSf(i,tt,gen)**2, MSf(j,tt,gen)**2)

				AL1V4 = 2.0d0*CF/(4d0*pi)**2*xc00*IL(i,tt,gen)*JR(j,tt,gen)*ZSfSf(i,j,tt,gen)
				AR1V4 = 2.0d0*CF/(4d0*pi)**2*xc00*IR(i,tt,gen)*JL(j,tt,gen)*ZSfSf(i,j,tt,gen)

				if (Denom_Mode==0) then

					axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.5d0*AR(tt,gen)*(AL1V4-AR1V4)/MZ2

				else

					axialfactor(uptype,gen) = axialfactor(uptype,gen) + 0.5d0*AR(tt,gen)*(AL1V4-AR1V4)/(MZ2-t_Mandel(tt,gen))

				endif

				endif

				enddo
			enddo
		enddo
	enddo

!----------------------- Squark contribution V5: Gluon Vertex, Outgoing Squark, s-Channel
	do gen = 1,3
!		gen = 3
		do uptype = 1,2
			do i = 1,2

			tt = uptype + 2
			if(DD_udsCheck(tt,gen)) then

			call myCget(Mneu(1)**2, s_Mandel(tt,gen), Mf(tt,gen)**2, Mf(tt,gen)**2, MSf(i,tt,gen)**2, 0d0)
			call myBget(s_Mandel(tt,gen), MSf(i, tt, gen)**2, 0d0)

			ALV5 = -(CF*Alfas/(4d0*pi)*(2.0d0*(xc0 + xc1 + xc2)*FR(1,i,tt,gen)*Mf(tt,gen)*MNeu(1)
     &			- FL(1,i,tt,gen)*(xb0 - 4.0d0*papb(tt,gen)*xc2 + 2.0d0*(xc0 - xc2)*Mf(tt,gen)**2
     &			+ 2.0d0*xc1*MNeu(1)**2)))

			ARV5 = -(CF*Alfas/(4d0*pi)*(2.0d0*(xc0 + xc1 + xc2)*FL(1,i,tt,gen)*Mf(tt,gen)*MNeu(1)
     &			- FR(1,i,tt,gen)*(xb0 - 4.0d0*papb(tt,gen)*xc2 + 2.0d0*(xc0 - xc2)*Mf(tt,gen)**2
     &			+ 2.0d0*xc1*MNeu(1)**2)))

			if (Denom_Mode==0) then

			axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &			(-(ARV5 + ALV5)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) +
     &			(ARV5 - ALV5)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/MSf(i,tt,gen)**2

			else

			axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &			(-(ARV5 + ALV5)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) +
     &			(ARV5 - ALV5)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/(MSf(i,tt,gen)**2-s_Mandel(tt,gen))

			endif

			endif

			enddo
		enddo
	enddo

!----------------------- Squark contribution V6: Gluino Vertex, Outgoing Squark, s-Channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
				do j = 1,2
				tt = uptype+2

				if(DD_udsCheck(tt,gen)) then

				call myCget(Mf(tt,gen)**2,s_Mandel(tt,gen),MNeu(1)**2,MSf(j,tt,gen)**2,MGl**2,Mf(tt,gen)**2)
				call myBget(s_Mandel(tt,gen), MGl**2, Mf(tt,gen)**2)

				ALV6 = CF/(4d0*pi)**2*(-(ER(1,j,tt,gen)*(xc2*IL(j,tt,gen)*IR(i,tt,gen)*Mf(tt,gen) +
     &				IL(i,tt,gen)*(MGl*(xc0 + xc2)*IL(j,tt,gen) + (xc0 + xc1 + xc2)*IR(j,tt,gen)*
     &				Mf(tt,gen)))*MNeu(1)) - ELc(1,j,tt,gen)*(IR(j,tt,gen)*Mf(tt,gen)*(MGl*xc1*
     &				IR(i,tt,gen) + (xc0 + xc1)*IL(i,tt,gen)*Mf(tt,gen)) + IL(j,tt,gen)*(MGl*xc0*
     &				IL(i,tt,gen)*Mf(tt,gen) + IR(i,tt,gen)*(xb0 + xc1*Mf(tt,gen)**2 + xc2*
     &				MNeu(1)**2 + xc0*MSf(j,tt,gen)**2))))

				ARV6 = CF/(4d0*pi)**2*(-(ELc(1,j,tt,gen)*(xc2*IL(i,tt,gen)*IR(j,tt,gen)*Mf(tt,gen) +
     &				IR(i,tt,gen)*(MGl*(xc0 + xc2)*IR(j,tt,gen) + (xc0 + xc1 + xc2)*IL(j,tt,gen)*
     &				Mf(tt,gen)))*MNeu(1)) - ER(1,j,tt,gen)*(IR(i,tt,gen)*Mf(tt,gen)*(MGl*xc0*
     &				IR(j,tt,gen) + (xc0 + xc1)*IL(j,tt,gen)*Mf(tt,gen)) + IL(i,tt,gen)*(MGl*xc1*
     &				IL(j,tt,gen)*Mf(tt,gen) + IR(j,tt,gen)*(xb0 + xc1*Mf(tt,gen)**2 + xc2*
     &				MNeu(1)**2 + xc0*MSf(j,tt,gen)**2))))

				if (Denom_Mode==0) then

				axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &				(-(ARV6 + ALV6)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) +
     &				(ARV6 - ALV6)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/MSf(i,tt,gen)**2

				else

				axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &				(-(ARV6 + ALV6)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) +
     &				(ARV6 - ALV6)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/(MSf(i,tt,gen)**2-s_Mandel(tt,gen))

				endif

				endif

				enddo
			enddo
		enddo
	enddo


!----------------------- Squark contribution V7: Gluon Vertex, Ingoing Squark, s-Channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2

			tt = uptype + 2
			if(DD_udsCheck(tt,gen)) then

			call myCget(Mneu(1)**2, s_Mandel(tt,gen), Mf(tt,gen)**2, Mf(tt,gen)**2, MSf(i,tt,gen)**2, 0d0)
			call myBget(s_Mandel(tt,gen), MSf(i, tt, gen)**2, 0d0)

			ALV7 = -(CF*Alfas/(4d0*pi)*(2.0d0*(xc0 + xc1 + xc2)*ER(1,i,tt,gen)*Mf(tt,gen)*MNeu(1)
     &			- ELc(1,i,tt,gen)*(xb0 - 4.0d0*p1p2(tt,gen)*xc2 + 2.0d0*(xc0 - xc2)*Mf(tt,gen)**2
     &			+ 2.0d0*xc1*MNeu(1)**2)))

			ARV7 = -(CF*Alfas/(4d0*pi)*(2.0d0*(xc0 + xc1 + xc2)*ELc(1,i,tt,gen)*Mf(tt,gen)*MNeu(1)
     &			- ER(1,i,tt,gen)*(xb0 - 4.0d0*p1p2(tt,gen)*xc2 + 2.0d0*(xc0 - xc2)*Mf(tt,gen)**2
     &			+ 2.0d0*xc1*MNeu(1)**2)))

			if (Denom_Mode==0) then

			axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &			(-(ARV7 + ALV7)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) +
     &			(ARV7 - ALV7)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/MSf(i,tt,gen)**2

			else

			axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &			(-(ARV7 + ALV7)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) +
     &			(ARV7 - ALV7)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/(MSf(i,tt,gen)**2-s_Mandel(tt,gen))

			endif

			endif

			enddo
		enddo
	enddo

!----------------------- Squark contribution V8: Gluino Vertex, Ingoing Squark, s-Channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
				do j = 1,2
				tt = uptype+2

				if(DD_udsCheck(tt,gen)) then

				call myCget(MNeu(1)**2,s_Mandel(tt,gen),Mf(tt,gen)**2,MSf(j,tt,gen)**2,Mf(tt,gen)**2,MGl**2)
				call myBget(s_Mandel(tt,gen), Mf(tt,gen)**2, MGl**2)

				ALV8 = CF/(4d0*pi)**2*(-(FR(1,j,tt,gen)*(xc1*JL(j,tt,gen)*JR(i,tt,gen)*Mf(tt,gen) +
     &				JL(i,tt,gen)*(MGl*(xc0 + xc1)*JL(j,tt,gen) + (xc0 + xc1 + xc2)*JR(j,tt,gen)*
     &				Mf(tt,gen)))*MNeu(1)) - FL(1,j,tt,gen)*(JR(j,tt,gen)*Mf(tt,gen)*(MGl*xc2*
     &				JR(i,tt,gen) + (xc0 + xc2)*JL(i,tt,gen)*Mf(tt,gen)) + JL(j,tt,gen)*(MGl*xc0*
     &				JL(i,tt,gen)*Mf(tt,gen) + JR(i,tt,gen)*(xb0 + xc2*Mf(tt,gen)**2 + xc1*
     &				MNeu(1)**2 + xc0*MSf(j,tt,gen)**2))))

				ARV8 = CF/(4d0*pi)**2*(-(FL(1,j,tt,gen)*(xc1*JL(i,tt,gen)*JR(j,tt,gen)*Mf(tt,gen) +
     &				JR(i,tt,gen)*(MGl*(xc0 + xc1)*JR(j,tt,gen) + (xc0 + xc1 + xc2)*JL(j,tt,gen)*
     &				Mf(tt,gen)))*MNeu(1)) - FR(1,j,tt,gen)*(JR(i,tt,gen)*Mf(tt,gen)*(MGl*xc0*
     &				JR(j,tt,gen) + (xc0 + xc2)*JL(j,tt,gen)*Mf(tt,gen)) + JL(i,tt,gen)*(MGl*xc2*
     &				JL(j,tt,gen)*Mf(tt,gen) + JR(j,tt,gen)*(xb0 + xc2*Mf(tt,gen)**2 + xc1*
     &				MNeu(1)**2 + xc0*MSf(j,tt,gen)**2))))

				if (Denom_Mode==0) then

				axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &				(-(ARV8 + ALV8)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) +
     &				(ARV8 - ALV8)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/MSf(i,tt,gen)**2

				else

				axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &				(-(ARV8 + ALV8)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) +
     &				(ARV8 - ALV8)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/(MSf(i,tt,gen)**2-s_Mandel(tt,gen))

				endif

				endif

				enddo
			enddo
		enddo
	enddo


!----------------------- Squark contribution V5: Gluon Vertex, Outgoing Squark, u-Channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2

			tt = uptype + 2
			if(DD_udsCheck(tt,gen)) then

			call myCget(Mneu(1)**2, u_Mandel(tt,gen), Mf(tt,gen)**2, Mf(tt,gen)**2, MSf(i,tt,gen)**2, 0d0)
			call myBget(u_Mandel(tt,gen), MSf(i, tt, gen)**2, 0d0)

			ALV5 = -(CF*Alfas/(4d0*pi)*(2.0d0*(xc0 + xc1 + xc2)*FR(1,i,tt,gen)*Mf(tt,gen)*MNeu(1)
     &			- FL(1,i,tt,gen)*(xb0 + 4.0d0*pap2(tt,gen)*xc2 + 2.0d0*(xc0 - xc2)*Mf(tt,gen)**2
     &			+ 2.0d0*xc1*MNeu(1)**2)))

			ARV5 = -(CF*Alfas/(4d0*pi)*(2.0d0*(xc0 + xc1 + xc2)*FL(1,i,tt,gen)*Mf(tt,gen)*MNeu(1)
     &			- FR(1,i,tt,gen)*(xb0 + 4.0d0*pap2(tt,gen)*xc2 + 2.0d0*(xc0 - xc2)*Mf(tt,gen)**2
     &			+ 2.0d0*xc1*MNeu(1)**2)))

			if (Denom_Mode==0) then

			axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &			(-(ARV5 + ALV5)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) +
     &			(ARV5 - ALV5)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/MSf(i,tt,gen)**2

			else

			axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &			(-(ARV5 + ALV5)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) +
     &			(ARV5 - ALV5)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/(MSf(i,tt,gen)**2-u_Mandel(tt,gen))

			endif

			endif

			enddo
		enddo
	enddo

!----------------------- Squark contribution V6: Gluino Vertex, Outgoing Squark, u-Channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
				do j = 1,2
				tt = uptype+2

				if(DD_udsCheck(tt,gen)) then

				call myCget(Mf(tt,gen)**2,u_Mandel(tt,gen),MNeu(1)**2,MSf(j,tt,gen)**2,MGl**2,Mf(tt,gen)**2)
				call myBget(u_Mandel(tt,gen), MGl**2, Mf(tt,gen)**2)

				ALV6 = CF/(4d0*pi)**2*(-(ER(1,j,tt,gen)*(xc2*IL(j,tt,gen)*IR(i,tt,gen)*Mf(tt,gen) +
     &				IL(i,tt,gen)*(MGl*(xc0 + xc2)*IL(j,tt,gen) + (xc0 + xc1 + xc2)*IR(j,tt,gen)*
     &				Mf(tt,gen)))*MNeu(1)) - ELc(1,j,tt,gen)*(IR(j,tt,gen)*Mf(tt,gen)*(MGl*xc1*
     &				IR(i,tt,gen) + (xc0 + xc1)*IL(i,tt,gen)*Mf(tt,gen)) + IL(j,tt,gen)*(MGl*xc0*
     &				IL(i,tt,gen)*Mf(tt,gen) + IR(i,tt,gen)*(xb0 + xc1*Mf(tt,gen)**2 + xc2*
     &				MNeu(1)**2 + xc0*MSf(j,tt,gen)**2))))

				ARV6 = CF/(4d0*pi)**2*(-(ELc(1,j,tt,gen)*(xc2*IL(i,tt,gen)*IR(j,tt,gen)*Mf(tt,gen) +
     &				IR(i,tt,gen)*(MGl*(xc0 + xc2)*IR(j,tt,gen) + (xc0 + xc1 + xc2)*IL(j,tt,gen)*
     &				Mf(tt,gen)))*MNeu(1)) - ER(1,j,tt,gen)*(IR(i,tt,gen)*Mf(tt,gen)*(MGl*xc0*
     &				IR(j,tt,gen) + (xc0 + xc1)*IL(j,tt,gen)*Mf(tt,gen)) + IL(i,tt,gen)*(MGl*xc1*
     &				IL(j,tt,gen)*Mf(tt,gen) + IR(j,tt,gen)*(xb0 + xc1*Mf(tt,gen)**2 + xc2*
     &				MNeu(1)**2 + xc0*MSf(j,tt,gen)**2))))

				if (Denom_Mode==0) then

				axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &				(-(ARV6 + ALV6)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) +
     &				(ARV6 - ALV6)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/MSf(i,tt,gen)**2

				else

				axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &				(-(ARV6 + ALV6)*(ER(1,i,tt,gen) + ELc(1,i,tt,gen)) +
     &				(ARV6 - ALV6)*(ER(1,i,tt,gen) - ELc(1,i,tt,gen)))/(MSf(i,tt,gen)**2-u_Mandel(tt,gen))

				endif

				endif

				enddo
			enddo
		enddo
	enddo

!----------------------- Squark contribution V7: Gluon Vertex, Ingoing Squark, u-Channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2

			tt = uptype + 2
			if(DD_udsCheck(tt,gen)) then

			call myCget(Mneu(1)**2, u_Mandel(tt,gen), Mf(tt,gen)**2, Mf(tt,gen)**2, MSf(i,tt,gen)**2, 0d0)
			call myBget(u_Mandel(tt,gen), MSf(i, tt, gen)**2, 0d0)

			ALV7 = -(CF*Alfas/(4d0*pi)*(2.0d0*(xc0 + xc1 + xc2)*ER(1,i,tt,gen)*Mf(tt,gen)*MNeu(1)
     &			- ELc(1,i,tt,gen)*(xb0 + 4.0d0*pap2(tt,gen)*xc2 + 2.0d0*(xc0 - xc2)*Mf(tt,gen)**2
     &			+ 2.0d0*xc1*MNeu(1)**2)))

			ARV7 = -(CF*Alfas/(4d0*pi)*(2.0d0*(xc0 + xc1 + xc2)*ELc(1,i,tt,gen)*Mf(tt,gen)*MNeu(1)
     &			- ER(1,i,tt,gen)*(xb0 + 4.0d0*pap2(tt,gen)*xc2 + 2.0d0*(xc0 - xc2)*Mf(tt,gen)**2
     &			+ 2.0d0*xc1*MNeu(1)**2)))

			if (Denom_Mode==0) then

			axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &			(-(ARV7 + ALV7)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) +
     &			(ARV7 - ALV7)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/MSf(i,tt,gen)**2

			else

			axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &			(-(ARV7 + ALV7)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) +
     &			(ARV7 - ALV7)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/(MSf(i,tt,gen)**2-u_Mandel(tt,gen))

			endif

			endif

			enddo
		enddo
	enddo

!----------------------- Squark contribution V8: Gluino Vertex, Ingoing Squark, u-Channel
	do gen = 1,3
		do uptype = 1,2
			do i = 1,2
				do j = 1,2
				tt = uptype+2

				if(DD_udsCheck(tt,gen)) then

				call myCget(MNeu(1)**2,u_Mandel(tt,gen),Mf(tt,gen)**2,MSf(j,tt,gen)**2,Mf(tt,gen)**2,MGl**2)
				call myBget(u_Mandel(tt,gen), Mf(tt,gen)**2, MGl**2)

				ALV8 = CF/(4d0*pi)**2*(-(FR(1,j,tt,gen)*(xc1*JL(j,tt,gen)*JR(i,tt,gen)*Mf(tt,gen) +
     &				JL(i,tt,gen)*(MGl*(xc0 + xc1)*JL(j,tt,gen) + (xc0 + xc1 + xc2)*JR(j,tt,gen)*
     &				Mf(tt,gen)))*MNeu(1)) - FL(1,j,tt,gen)*(JR(j,tt,gen)*Mf(tt,gen)*(MGl*xc2*
     &				JR(i,tt,gen) + (xc0 + xc2)*JL(i,tt,gen)*Mf(tt,gen)) + JL(j,tt,gen)*(MGl*xc0*
     &				JL(i,tt,gen)*Mf(tt,gen) + JR(i,tt,gen)*(xb0 + xc2*Mf(tt,gen)**2 + xc1*
     &				MNeu(1)**2 + xc0*MSf(j,tt,gen)**2))))

				ARV8 = CF/(4d0*pi)**2*(-(FL(1,j,tt,gen)*(xc1*JL(i,tt,gen)*JR(j,tt,gen)*Mf(tt,gen) +
     &				JR(i,tt,gen)*(MGl*(xc0 + xc1)*JR(j,tt,gen) + (xc0 + xc1 + xc2)*JL(j,tt,gen)*
     &				Mf(tt,gen)))*MNeu(1)) - FR(1,j,tt,gen)*(JR(i,tt,gen)*Mf(tt,gen)*(MGl*xc0*
     &				JR(j,tt,gen) + (xc0 + xc2)*JL(j,tt,gen)*Mf(tt,gen)) + JL(i,tt,gen)*(MGl*xc2*
     &				JL(j,tt,gen)*Mf(tt,gen) + JR(j,tt,gen)*(xb0 + xc2*Mf(tt,gen)**2 + xc1*
     &				MNeu(1)**2 + xc0*MSf(j,tt,gen)**2))))

				if (Denom_Mode==0) then

				axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &				(-(ARV8 + ALV8)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) +
     &				(ARV8 - ALV8)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/MSf(i,tt,gen)**2

				else

				axialfactor(uptype,gen) = axialfactor(uptype,gen) - 1.0d0/16.0d0*
     &				(-(ARV8 + ALV8)*(FR(1,i,tt,gen) + FL(1,i,tt,gen)) +
     &				(ARV8 - ALV8)*(FR(1,i,tt,gen) - FL(1,i,tt,gen)))/(MSf(i,tt,gen)**2-u_Mandel(tt,gen))

				endif

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
				g_4Fa(nucleon) = g_4Fa(nucleon) + Delta_q(nucleon,uptype,gen)
     &							 *axialfactor(uptype,gen)*Running_SD(tt,gen)
				g_4FaContribs(nucleon,2) = g_4FaContribs(nucleon,2) + Delta_q(nucleon,uptype,gen)
     &							 *axialfactor(uptype,gen)*Running_SD(tt,gen)
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


