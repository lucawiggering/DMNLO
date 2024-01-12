! =====================================================================
! In this file Self-Energies and Renormalization constants are taken from NeuQ2qx_SetCounterterms.F, while Propagator and Coupling Counter-terms
! are the one calculated by QLB for the Vector boson case (gluon not implemented yet)
! =====================================================================

	subroutine NeuQ2qx_SetCountertermsV()
	implicit none

#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_GenCouplings.h"
#include "NeuQ2qx_GenCounterterms.h"

	integer i,j,k,l,m,type,gen
        integer typex, sfmeigenin, sfmeigenout, sfgaugeeigenin, sfgaugeeigenout
	integer ttf, genf, ttsf, gensf, eigenchi, eigensf
	double precision mbmass,mtmass, mtit, MSf2temp(2,4,3)
        double precision CL1, CL2, Krondelta

        double precision CLR11,CLR22
        double complex dCLnew

	external NeuQ2qx_CalcPifV, NeuQ2qx_CalcdPifV, NeuQ2qx_CalcPiSfV, NeuQ2qx_CalcdPiSfV
	external NeuQ2qx_CalcPropSquarkV, NeuQ2qx_CalcPropQuarkV



! ===================================
! Coupling counterterms 
! ===================================

! Neutralino/Chargino - Sfermion - Fermion

        do eigensf = 1, 2
 	        do ttsf = 1, 4
 	          do gensf = 1, 3
 	            do eigenchi = 1, 4

! --- incoming up-type sfermion ---

          if ((ttsf.eq.1).or.(ttsf.eq.3)) then

 	              dAL(eigenchi,eigensf,ttsf,gensf) = dconjg(dZfR(ttsf,gensf))*AL(eigenchi,eigensf,ttsf,gensf)/2. + dZSf(1,eigensf,ttsf,gensf)*AL(eigenchi,1,ttsf,gensf)/2. +
     &          dZSf(2,eigensf,ttsf,gensf)*AL(eigenchi,2,ttsf,gensf)/2. - hf(ttsf,gensf)*ZNeuC(eigenchi,4)*(dUSf(eigensf,1,ttsf,gensf) + USf(eigensf,1,ttsf,gensf)*dMf(ttsf,gensf)/Mf(ttsf,gensf)) +
     &          (sqrt2*EL*ef(ttsf)/CW)*ZNeuC(eigenchi,1)*dUSf(eigensf,2,ttsf,gensf)
 	        
                dAR(eigenchi,eigensf,ttsf,gensf) = dconjg(dZfL(ttsf,gensf))*AR(eigenchi,eigensf,ttsf,gensf)/2. + dZSf(1,eigensf,ttsf,gensf)*AR(eigenchi,1,ttsf,gensf)/2. +
     &          dZSf(2,eigensf,ttsf,gensf)*AR(eigenchi,2,ttsf,gensf)/2. - hf(ttsf,gensf)*ZNeu(eigenchi,4)*(dUSf(eigensf,2,ttsf,gensf) + USf(eigensf,2,ttsf,gensf)*dMf(ttsf,gensf)/Mf(ttsf,gensf)) -
     &          (sqrt2*EL/SW)*((ef(ttsf) - I3f(ttsf))*ZNeu(eigenchi,1)*SW/CW + I3f(ttsf)*ZNeu(eigenchi,2))*dUSf(eigensf,1,ttsf,gensf)

     
! --- incoming down-type sfermion---     
	
          else

                dAL(eigenchi,eigensf,ttsf,gensf) = dconjg(dZfR(ttsf,gensf))*AL(eigenchi,eigensf,ttsf,gensf)/2. + dZSf(1,eigensf,ttsf,gensf)*AL(eigenchi,1,ttsf,gensf)/2. +
     &          dZSf(2,eigensf,ttsf,gensf)*AL(eigenchi,2,ttsf,gensf)/2. - hf(ttsf,gensf)*ZNeuC(eigenchi,3)*(dUSf(eigensf,1,ttsf,gensf) + USf(eigensf,1,ttsf,gensf)*dMf(ttsf,gensf)/Mf(ttsf,gensf)) +
     &          (sqrt2*EL*ef(ttsf)/CW)*ZNeuC(eigenchi,1)*dUSf(eigensf,2,ttsf,gensf)
        	        
                dAR(eigenchi,eigensf,ttsf,gensf) = dconjg(dZfL(ttsf,gensf))*AR(eigenchi,eigensf,ttsf,gensf)/2. + dZSf(1,eigensf,ttsf,gensf)*AR(eigenchi,1,ttsf,gensf)/2. +
     &          dZSf(2,eigensf,ttsf,gensf)*AR(eigenchi,2,ttsf,gensf)/2. - hf(ttsf,gensf)*ZNeu(eigenchi,3)*(dUSf(eigensf,2,ttsf,gensf) + USf(eigensf,2,ttsf,gensf)*dMf(ttsf,gensf)/Mf(ttsf,gensf)) -
     &          (sqrt2*EL/SW)*((ef(ttsf) - I3f(ttsf))*ZNeu(eigenchi,1)*SW/CW + I3f(ttsf)*ZNeu(eigenchi,2))*dUSf(eigensf,1,ttsf,gensf)
        
          endif
        
 	               
 	            enddo

              do eigenchi = 5, 6

! --- incoming up-type sfermion ---

          if ((ttsf.eq.1).or.(ttsf.eq.3)) then

            dAL(eigenchi,eigensf,ttsf,gensf) = dconjg(dZfR(ttsf+1,gensf))*AL(eigenchi,eigensf,ttsf,gensf)/2. +
     &                                         (dUSf(eigensf,1,ttsf,gensf)/USf(eigensf,1,ttsf,gensf))*AL(eigenchi,eigensf,ttsf,gensf) +
     &                                         (dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))*AL(eigenchi,eigensf,ttsf,gensf) +
     &                                         dZSf(1,eigensf,ttsf,gensf)*AL(eigenchi,1,ttsf,gensf)/2. + dZSf(2,eigensf,ttsf,gensf)*AL(eigenchi,2,ttsf,gensf)/2.

            dAR(eigenchi,eigensf,ttsf,gensf) = dconjg(dZfL(ttsf+1,gensf))*AR(eigenchi,eigensf,ttsf,gensf)/2. +
     &                                         dUSf(eigensf,1,ttsf,gensf)*(-EL/SW*VCha(eigenchi-4,1)) +
     &                                         dUSf(eigensf,2,ttsf,gensf)*(EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*SB)*VCha(eigenchi-4,2)) +
     &                                         dMf(ttsf,gensf)*(EL/(SW*sqrt2*MW*SB)*VCha(eigenchi-4,2)*USf(eigensf,2,ttsf,gensf)) +
     &                                         dZSf(1,eigensf,ttsf,gensf)*AR(eigenchi,1,ttsf,gensf)/2. + dZSf(2,eigensf,ttsf,gensf)*AR(eigenchi,2,ttsf,gensf)/2.

! --- incoming down-type sfermion---     
  
          else

            dAL(eigenchi,eigensf,ttsf,gensf) = dconjg(dZfR(ttsf-1,gensf))*AL(eigenchi,eigensf,ttsf,gensf)/2. +
     &                                         (dUSf(eigensf,1,ttsf,gensf)/USf(eigensf,1,ttsf,gensf))*AL(eigenchi,eigensf,ttsf,gensf) +
     &                                         (dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))*AL(eigenchi,eigensf,ttsf,gensf) +
     &                                         dZSf(1,eigensf,ttsf,gensf)*AL(eigenchi,1,ttsf,gensf)/2. + dZSf(2,eigensf,ttsf,gensf)*AL(eigenchi,2,ttsf,gensf)/2.

            dAR(eigenchi,eigensf,ttsf,gensf) = dconjg(dZfL(ttsf-1,gensf))*AR(eigenchi,eigensf,ttsf,gensf)/2. +
     &                                         dUSf(eigensf,1,ttsf,gensf)*(-EL/SW*UCha(eigenchi-4,1)) +
     &                                         dUSf(eigensf,2,ttsf,gensf)*(EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*CB)*UCha(eigenchi-4,2)) +
     &                                         dMf(ttsf,gensf)*(EL/(SW*sqrt2*MW*CB)*UCha(eigenchi-4,2)*USf(eigensf,2,ttsf,gensf)) +
     &                                         dZSf(1,eigensf,ttsf,gensf)*AR(eigenchi,1,ttsf,gensf)/2. + dZSf(2,eigensf,ttsf,gensf)*AR(eigenchi,2,ttsf,gensf)/2.
        
          endif
        
                 
              enddo
 	          enddo
 	        enddo
 	      enddo  

! Z/W/Photon - Fermion - Fermion

        do typex = 2, 4
          do ttf = 1, 4
            do genf = 1, 3

              if (typex.eq.4) then

! --- outgoing down-type fermion ---

          if ((ttf.eq.2).or.(ttf.eq.4)) then

          dBL(typex,ttf,genf) = (dZfL(ttf-1,genf) + dconjg(dZfL(ttf,genf)))*BL(typex,ttf)/2.
          dBR(typex,ttf,genf) = (dZfR(ttf-1,genf) + dconjg(dZfR(ttf,genf)))*BR(typex,ttf)/2.

! --- outgoing up-type fermion---     
  
          else

          dBL(typex,ttf,genf) = (dZfL(ttf+1,genf) + dconjg(dZfL(ttf,genf)))*BL(typex,ttf)/2.
          dBR(typex,ttf,genf) = (dZfR(ttf+1,genf) + dconjg(dZfR(ttf,genf)))*BR(typex,ttf)/2.

          endif   

              else

              dBL(typex,ttf,genf) = (dZfL(ttf,genf) + dconjg(dZfL(ttf,genf)))*BL(typex,ttf)/2.
              dBR(typex,ttf,genf) = (dZfR(ttf,genf) + dconjg(dZfR(ttf,genf)))*BR(typex,ttf)/2.

              endif     
        
            enddo
          enddo
        enddo

! Z - Sfermion - Sfermion

        do sfmeigenout = 1, 2
          do eigensf = 1, 2
            do ttsf = 3, 4
              do gensf = 3, 3

              CL1 = - g2h1*I3f(ttsf)/CW
              CL2 = g2h1*SW**2*ef(ttsf)/CW

        
              dCL(sfmeigenout,eigensf,ttsf,gensf,3) = CL1*(dUSf(eigensf,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf,gensf) + USf(eigensf,1,ttsf,gensf)*dUSf(sfmeigenout,1,ttsf,gensf)) +
     &        (CL1*USf(eigensf,1,ttsf,gensf)*USf(1,1,ttsf,gensf) + CL2*Krondelta(eigensf,1))*dconjg(dZSf(1,sfmeigenout,ttsf,gensf))/2. +
     &        (CL1*USf(eigensf,1,ttsf,gensf)*USf(2,1,ttsf,gensf) + CL2*Krondelta(eigensf,2))*dconjg(dZSf(2,sfmeigenout,ttsf,gensf))/2. +
     &        (CL1*USf(1,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf,gensf) + CL2*Krondelta(1,sfmeigenout))*dZSf(1,eigensf,ttsf,gensf)/2. +
     &        (CL1*USf(2,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf,gensf) + CL2*Krondelta(2,sfmeigenout))*dZSf(2,eigensf,ttsf,gensf)/2.

     		   CLR11= -g2h1/CW*(I3f(ttsf)-SW**2*ef(ttsf))
     		   CLR22= -g2h1/CW*(-SW**2*ef(ttsf))

     		  dCLnew = 0.5d0*CL(1,eigensf,ttsf,gensf,3)*dconjg(dZSf(1,sfmeigenout,ttsf,gensf))
     &		          +0.5d0*CL(2,eigensf,ttsf,gensf,3)*dconjg(dZSf(2,sfmeigenout,ttsf,gensf))
     &   		      +0.5d0*CL(sfmeigenout,1,ttsf,gensf,3)*dZSf(1,eigensf,ttsf,gensf)
     &   		      +0.5d0*CL(sfmeigenout,2,ttsf,gensf,3)*dZSf(2,eigensf,ttsf,gensf)
     &   		      +dUSf(sfmeigenout,1,ttsf,gensf)*CLR11*USf(eigensf,1,ttsf,gensf)
     &   		      +dUSf(sfmeigenout,2,ttsf,gensf)*CLR22*USf(eigensf,2,ttsf,gensf)
     &   		      +USf(sfmeigenout,1,ttsf,gensf)*CLR11*dUSf(eigensf,1,ttsf,gensf)
     &   		      +USf(sfmeigenout,2,ttsf,gensf)*CLR22*dUSf(eigensf,2,ttsf,gensf)

    
              enddo
            enddo
          enddo
        enddo        

! Photon - Sfermion - Sfermion


        do sfmeigenout = 1, 2
          do eigensf = 1, 2
            do ttsf = 1, 4
              do gensf = 1, 3

              CL1 = 0d0
              CL2 = -ef(ttsf)*EL
        
              dCL(sfmeigenout,eigensf,ttsf,gensf,2) = CL2*(dconjg(dZSf(eigensf,sfmeigenout,ttsf,gensf)) + dZSf(sfmeigenout,eigensf,ttsf,gensf))/2.
             
              enddo
            enddo
          enddo
        enddo

! W - Sfermion - Sfermion

        do sfmeigenout = 1, 2
          do eigensf = 1, 2
            do ttsf = 1, 4
              do gensf = 1, 3                      

              CL1 = - g2h1/sqrt2
              CL2 = 0d0

        
! --- incoming up-type sfermion ---

              if ((ttsf.eq.1).or.(ttsf.eq.3)) then

              dCL(sfmeigenout,eigensf,ttsf,gensf,4) = CL1*(dUSf(eigensf,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf+1,gensf) + USf(eigensf,1,ttsf,gensf)*dUSf(sfmeigenout,1,ttsf+1,gensf)) +
     &        CL1*USf(eigensf,1,ttsf,gensf)*USf(1,1,ttsf+1,gensf)*dconjg(dZSf(1,sfmeigenout,ttsf+1,gensf))/2. +
     &        CL1*USf(eigensf,1,ttsf,gensf)*USf(2,1,ttsf+1,gensf)*dconjg(dZSf(2,sfmeigenout,ttsf+1,gensf))/2. +
     &        CL1*USf(1,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf+1,gensf)*dZSf(1,eigensf,ttsf,gensf)/2. +
     &        CL1*USf(2,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf+1,gensf)*dZSf(2,eigensf,ttsf,gensf)/2.
        
! --- incoming down-type sfermion ---

              else

              dCL(sfmeigenout,eigensf,ttsf,gensf,4) = CL1*(dUSf(eigensf,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf-1,gensf) + USf(eigensf,1,ttsf,gensf)*dUSf(sfmeigenout,1,ttsf-1,gensf)) +
     &        CL1*USf(eigensf,1,ttsf,gensf)*USf(1,1,ttsf-1,gensf)*dconjg(dZSf(1,sfmeigenout,ttsf-1,gensf))/2. +
     &        CL1*USf(eigensf,1,ttsf,gensf)*USf(2,1,ttsf-1,gensf)*dconjg(dZSf(2,sfmeigenout,ttsf-1,gensf))/2. +
     &        CL1*USf(1,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf-1,gensf)*dconjg(dZSf(1,eigensf,ttsf,gensf))/2. +
     &        CL1*USf(2,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf-1,gensf)*dconjg(dZSf(2,eigensf,ttsf,gensf))/2.
        
              endif
        
              enddo
            enddo
          enddo
        enddo
        

 	end

! ! =======================================
! ! Propagator counterterms 
! ! =======================================
! 
! ! Squarks
! 	subroutine NeuQ2qx_CalcPropSquarkV(myP2,type,gen)
! 
! #include "NeuQ2qx_Model.h"
! #include "NeuQ2qx_GenCounterterms.h"
! #include "../util/DMNLO_LoopIntegrals.h"
! 
! 	integer type,gen,i,j
! 	double precision myP2, Krondelta
! 
!         do i=1,2
!                 do j=1,2
!                 Prop(i,j,type,gen) = (myP2 - Msf(i,type,gen)**2)*dZSf(i,j,type,gen)/2 + (myP2 - Msf(j,type,gen)**2)*dconjg(dZSf(j,i,type,gen))/2 - Krondelta(i,j)*dMSf2(i,type,gen)
!                 enddo
!         enddo
!         
!         
! 	end
! 
! 
! ! Quarks
! 	subroutine NeuQ2qx_CalcPropQuarkV(type,gen)
! 
! #include "NeuQ2qx_Model.h"
! #include "NeuQ2qx_GenCounterterms.h"
! #include "../util/DMNLO_LoopIntegrals.h"
! 
! 	integer type,gen
! 
!         PropSL(type,gen) = - Mf(type,gen)*(dZfL(type,gen) + dconjg(dZfR(type,gen)))/2 - dMf(type,gen)
!         PropSR(type,gen) = - Mf(type,gen)*(dZfR(type,gen) + dconjg(dZfL(type,gen)))/2 - dMf(type,gen)             
!         PropVL(type,gen) = dreal(dZfR(type,gen))
!         PropVR(type,gen) = dreal(dZfL(type,gen))
! 
! 
! 	end
! 
