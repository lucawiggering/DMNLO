!**********************************
!**** subroutine to set the kinematics of the ingoing 4-momenta and magnitude of
!**** in/out-going momenta
!***********************************

	subroutine stst2qq_SetKinematics

#include "stst2QQ_Kinematics.h"
#include "stst2QQ_Model.h"

	muu1 = MSf(iSf1,itt1,igen1)/sqrtS
	muu2 = MSf(iSf2,itt2,igen2)/sqrtS		
	muu3 = Mf(ftt1,fgen1)/sqrtS
	muu4 = Mf(ftt2,fgen2)/sqrtS
	muu5 = 0d0
	
	m1 = MSf(iSf1,itt1,igen1)
	m2 = MSf(iSf2,itt2,igen2)
	m3 = Mf(ftt1,fgen1)
	m4 = Mf(ftt2,fgen2)
	m5 = 0d0

!Reduced masses of propagator
	do i=1,4
		muuqt(i) = MNeu(i)/sqrtS
		muupt(i) = MNeu(i)/sqrtS
	enddo	
	muuqt(5) = MCha(1)/sqrtS
	muupt(5) = MCha(1)/sqrtS
	muuqt(6) = MCha(2)/sqrtS
	muupt(6) = MCha(2)/sqrtS
	
	do i=1,4
		muuqu(i) = MNeu(i)/sqrtS
		muupu(i) = MNeu(i)/sqrtS
	enddo	
	muuqu(5) = MCha(1)/sqrtS
	muupu(5) = MCha(1)/sqrtS
	muuqu(6) = MCha(2)/sqrtS
	muupu(6) = MCha(2)/sqrtS
	
	muugl = MGl/sqrtS	
	
	  do i=1,4
		MGau(i) = MNeu(i)
      enddo
    
      MGau(5) = MCha(1)
      MGau(6) = MCha(2)

	end
	
