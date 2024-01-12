!==================================================================================================
!======================= Kinematics for direct detection in DM@NLO  ===============================
!=======================        Author: Patrick Steppeler		    ===============================
!==================================================================================================

	subroutine DD_SetKinematics(p_cm, costh)	
	implicit none

!======================= INITIALISATION AND SETUP =================================================

!----------------------- Include external header files --------------------------------------------
#include "DD_Model.h"
	
!----------------------- Declare variables --------------------------------------------------------
	double precision p_cm, costh, kappa, q, P
	integer i, tt, gen
	
!======================= ACTUAL CALCULATIONS ======================================================	
	
!----------------------- Set Mandelstam variables -------------------------------------------------
	do tt = 3,4
		do gen = 1,3

! General kinematics
			if(p_cm.ne.0d0) then
				s_Mandel(tt,gen) = MNeu(1)**2 + 2.0d0*p_cm**2 + Mf(tt,gen)**2 + 2.0d0*
     &					dsqrt(MNeu(1)**2 + p_cm**2)*dsqrt(Mf(tt,gen)**2 +p_cm**2)
				t_Mandel(tt,gen) = 2.0d0*MNeu(1)**2 - 0.5d0*(s_Mandel(tt,gen) + MNeu(1)**2 - Mf(tt,gen)**2)**2/s_Mandel(tt,gen) + 
     &					0.5d0*kappa(s_Mandel(tt,gen), MNeu(1)**2, Mf(tt,gen)**2)**2/s_Mandel(tt,gen)*costh
				u_Mandel(tt,gen) = 2.0d0*MNeu(1)**2 + 2.0d0*Mf(tt,gen)**2 - t_Mandel(tt,gen) - s_Mandel(tt,gen)
				
! Nonrelativistic approximation for pcm = 0, this choice sets the whole code pcm independent
			else
				s_Mandel(tt,gen) = (MNeu(1) + Mf(tt,gen))**2
				t_Mandel(tt,gen) = 0d0
				u_Mandel(tt,gen) = (MNeu(1) - Mf(tt,gen))**2
			endif
	
cc			write(*,*)'tt,gen =',tt,gen
cc			write(*,*)'Mf(tt,gen)**2',Mf(tt,gen)**2
cc			write(*,*)'MNeu(1)**2',MNeu(1)**2
cc			write(*,*)'s_Mandel(tt,gen)',s_Mandel(tt,gen)
cc			write(*,*)'t_Mandel(tt,gen)',t_Mandel(tt,gen)
cc			write(*,*)'u_Mandel(tt,gen)',u_Mandel(tt,gen)
			
		enddo
	enddo
	
!----------------------- Set products of four momenta ---------------------------------------------
	do tt =3,4
		do gen = 1,3
			pbp2(tt,gen) = Mf(tt,gen)**2 - 0.5d0*t_Mandel(tt,gen)
			pap2(tt,gen) = 0.5d0*(Mf(tt,gen)**2 + MNeu(1)**2 - u_Mandel(tt,gen))
			p1p2(tt,gen) = 0.5d0*(s_Mandel(tt,gen) - MNeu(1)**2 - Mf(tt,gen)**2)
			papb(tt,gen) = 0.5d0*(s_Mandel(tt,gen) - MNeu(1)**2 - Mf(tt,gen)**2)
		enddo
	enddo
	
!----------------------- Calculate Suppresion Factors ---------------------------------------------
	
!======================= END OF PROGRAM ===========================================================

	end


