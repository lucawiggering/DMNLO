!**********************************
!**** subroutine to set the kinematics of the ingoing 4-momenta and magnitude of
!**** in/out-going momenta
!***********************************

	subroutine ChiChi2QQ_SetKinematics

#include "ChiChi2QQ_Kinematics.h"
#include "ChiChi2QQ_Model.h"

	integer iMO

! No charginos in initial state
	if (icha.eq.0) then
	   muu1 = MNeu(indi)/sqrtS
	   muu2 = MNeu(indj)/sqrtS
	   muu3 = Mf(tt,gen)/sqrtS
	   muu4 = Mf(tt,gen)/sqrtS

! One chargino in initial state
	else if (icha.eq.1) then
	   muu1 = MCha(indi)/sqrtS
	   muu2 = MNeu(indj)/sqrtS
	   muu3 = Mf(tt,gen)/sqrtS
	   muu4 = Mf(tt+1,gen)/sqrtS	

! Two charginos in initial state
	else if (icha.eq.2) then
	   muu1 = MCha(indi)/sqrtS
	   muu2 = MCha(indj)/sqrtS
	   muu3 = Mf(tt,gen)/sqrtS
	   muu4 = Mf(tt,gen)/sqrtS		

	endif

! Gluon and gluino masses	
	muu5 = 0d0
	MGlr = MGl/sqrtS

c masses for dipoles
      m1 = muu1*sqrtS
      m2 = muu2*sqrtS
      m3 = muu3*sqrtS
      m4 = muu4*sqrtS



	end
