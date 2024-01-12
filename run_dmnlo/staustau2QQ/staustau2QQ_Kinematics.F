!**********************************
!**** subroutine to set the kinematics of the ingoing 4-momenta and magnitude of
!**** in/out-going momenta
!***********************************

	subroutine staustau2QQ_SetKinematics

#include "staustau2QQ_Kinematics.h"
#include "staustau2QQ_Model.h"

	muu1 = MSf(isf1,itt1,igen1)/sqrtS
	muu2 = MSf(isf2,itt2,igen2)/sqrtS
	muu3 = Mf(ftt1,fgen1)/sqrtS
	muu4 = Mf(ftt2,fgen2)/sqrtS

! Gluon and gluino masses	
	muu5 = 0d0

	end
	
