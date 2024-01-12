! =====================================================================
	subroutine NeuQ2qx_SetCounterterms()
	implicit none

#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_GenCouplings.h"
#include "NeuQ2qx_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer i,j,k,l,m,type,gen
      integer typex, sfmeigenin, sfmeigenout, sfgaugeeigenin, sfgaugeeigenout
	integer ttf, genf, ttsf, gensf, eigenchi, eigensf
	double precision mbmass,mtmass, mtit, MSf2temp(2,4,3)
	double precision CL1, CL2
      double precision CLR11,CLR22
      double complex dCLnew

	double complex Gnew

! functions from DMNLO_Kinematics.F
      double precision Krondelta, ef, I3f, hf

	external NeuQ2qx_CalcPropSquark, NeuQ2qx_CalcPropQuark
        external NeuQ2qx_CalcPifV, NeuQ2qx_CalcdPifV, NeuQ2qx_CalcPiSfV, NeuQ2qx_CalcdPiSfV
	external NeuQ2qx_CalcPropSquarkV, NeuQ2qx_CalcPropQuarkV
	
! ===================================
! Coupling counterterms 
! ===================================

! ===== Neu-Sf-f -- A-coupling =====

!  e.g. AL/R(chimeigen, sfmeigen, sftype, sfgen)
	 
         do eigensf = 1, 2
 	   do ttsf = 1, 4
 	      do gensf = 1, 3
 	      

! ------------------------------ Vertex Neutralino - Sfermion (in) - Fermion (out) ------------------------------

 	         do eigenchi = 1, 4

! --- incoming up-type sfermion ---

 	               if ((ttsf.eq.1).or.(ttsf.eq.3)) then
 	               
 	               	  dAR(eigenchi,eigensf,ttsf,gensf) = - EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*SB)*ZNeu(eigenchi,4)*USf(eigensf,2,ttsf,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dUSf(eigensf,2,ttsf,gensf)/USf(eigensf,2,ttsf,gensf)) -
     &                    sqrt2*EL/SW*((ef(ttsf) - I3f(ttsf))*SW/CW*ZNeu(eigenchi,1) +  I3f(ttsf)*ZNeu(eigenchi,2))*USf(eigensf,1,ttsf,gensf)*(dUSf(eigensf,1,ttsf,gensf)/USf(eigensf,1,ttsf,gensf))+
     &                    AR(eigenchi,1,ttsf,gensf)*1/2d0*dZSf(1,eigensf,ttsf,gensf) + AR(eigenchi,2,ttsf,gensf)*1/2d0*dZSf(2,eigensf,ttsf,gensf) +
     &                    AR(eigenchi,eigensf,ttsf,gensf)*1/2d0*dconjg(dZfL(ttsf,gensf))
     
                          dAL(eigenchi,eigensf,ttsf,gensf) =  - EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*SB)*ZNeuC(eigenchi,4)*USf(eigensf,1,ttsf,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dUSf(eigensf,1,ttsf,gensf)/USf(eigensf,1,ttsf,gensf)) + 
     &                    sqrt2*EL*ef(ttsf)/CW*ZNeuC(eigenchi,1)*USf(eigensf,2,ttsf,gensf)*(dUSf(eigensf,2,ttsf,gensf)/USf(eigensf,2,ttsf,gensf))+
     & 	                  AL(eigenchi,1,ttsf,gensf)*1/2d0*dZSf(1,eigensf,ttsf,gensf) + AL(eigenchi,2,ttsf,gensf)*1/2d0*dZSf(2,eigensf,ttsf,gensf) + 
     &                    AL(eigenchi,eigensf,ttsf,gensf)*1/2d0*dconjg(dZfR(ttsf,gensf))	

! --- incoming down-type sfermion---     

      		       else
!       		       
                           dAR(eigenchi,eigensf,ttsf,gensf) =  - EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*CB)*ZNeu(eigenchi,3)*USf(eigensf,2,ttsf,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dUSf(eigensf,2,ttsf,gensf)/USf(eigensf,2,ttsf,gensf)) -
     &                     sqrt2*EL/SW*(I3f(ttsf)*ZNeu(eigenchi,2) + (ef(ttsf) - I3f(ttsf))*SW/CW*ZNeu(eigenchi,1))*USf(eigensf,1,ttsf,gensf)*(dUSf(eigensf,1,ttsf,gensf)/USf(eigensf,1,ttsf,gensf))+
     & 	                   AR(eigenchi,1,ttsf,gensf)*1/2d0*dZSf(1,eigensf,ttsf,gensf) + AR(eigenchi,2,ttsf,gensf)*1/2d0*dZSf(2,eigensf,ttsf,gensf) + 
     &                     AR(eigenchi,eigensf,ttsf,gensf)*1/2d0*dconjg(dZfL(ttsf,gensf))
! 
                           dAL(eigenchi,eigensf,ttsf,gensf) =  - EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*CB)*ZNeuC(eigenchi,3)*USf(eigensf,1,ttsf,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dUSf(eigensf,1,ttsf,gensf)/USf(eigensf,1,ttsf,gensf)) + 
     &                     sqrt2*EL*ef(ttsf)/CW*ZNeuC(eigenchi,1)*USf(eigensf,2,ttsf,gensf)*(dUSf(eigensf,2,ttsf,gensf)/USf(eigensf,2,ttsf,gensf))+
     & 	                   AL(eigenchi,1,ttsf,gensf)*1/2d0*dZSf(1,eigensf,ttsf,gensf) + AL(eigenchi,2,ttsf,gensf)*1/2d0*dZSf(2,eigensf,ttsf,gensf) + 
     &                     AL(eigenchi,eigensf,ttsf,gensf)*1/2d0*dconjg(dZfR(ttsf,gensf))
	               endif
	          enddo
! ------------------------------ Vertex Chargino - Sfermion (in) - Fermion (out) ------------------------------

	          do eigenchi= 5, 6   

! --- incoming up-type sfermion ---	          

	               if ((ttsf.eq.1).or.(ttsf.eq.3)) then
 	         	
	         	dAL(eigenchi,eigensf,ttsf,gensf) =  EL*Mf(ttsf+1,gensf)/(SW*sqrt2*MW*CB)*UChaC(eigenchi-4,2)*USf(eigensf,1,ttsf,gensf)*(dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf) + dUSf(eigensf,1,ttsf,gensf)/USf(eigensf,1,ttsf,gensf)) + 
     &                  AL(eigenchi,1,ttsf,gensf)*1/2d0*dZSf(1,eigensf,ttsf,gensf) +
     &                  AL(eigenchi,2,ttsf,gensf)*1/2d0*dZSf(2,eigensf,ttsf,gensf) +
     &                  AL(eigenchi,eigensf,ttsf,gensf)*1/2d0*dconjg(dZfR(ttsf+1,gensf))

 	         	dAR(eigenchi,eigensf,ttsf,gensf) = -EL/SW*VCha(eigenchi-4,1)*USf(eigensf,1,ttsf,gensf)*(dUSf(eigensf,1,ttsf,gensf)/USf(eigensf,1,ttsf,gensf)) + EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*SB)*VCha(eigenchi-4,2)*USf(eigensf,2,ttsf,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dUSf(eigensf,2,ttsf,gensf)/USf(eigensf,2,ttsf,gensf)) + 
     &                  AR(eigenchi,1,ttsf,gensf)*1/2d0*dZSf(1,eigensf,ttsf,gensf) +
     &                  AR(eigenchi,2,ttsf,gensf)*1/2d0*dZSf(2,eigensf,ttsf,gensf) +
     &                  AR(eigenchi,eigensf,ttsf,gensf)*1/2d0*dconjg(dZfL(ttsf+1,gensf))

! --- incoming down-type sfermion --- 	 

 	         	else

 	         	dAL(eigenchi,eigensf,ttsf,gensf) = EL*Mf(ttsf-1,gensf)/(SW*sqrt2*MW*SB)*VChaC(eigenchi-4,2)*USf(eigensf,1,ttsf,gensf)*(dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf) + (dUSf(eigensf,1,ttsf,gensf)/USf(eigensf,1,ttsf,gensf))) + 
     &                  AL(eigenchi,1,ttsf,gensf)*1/2d0*dZSf(1,eigensf,ttsf,gensf) +
     &                  AL(eigenchi,2,ttsf,gensf)*1/2d0*dZSf(2,eigensf,ttsf,gensf) +
     &                  AL(eigenchi,eigensf,ttsf,gensf)*1/2d0*dconjg(dZfR(ttsf-1,gensf))

 	         	dAR(eigenchi,eigensf,ttsf,gensf) = -EL/SW*UCha(eigenchi-4,1)*USf(eigensf,1,ttsf,gensf)*(dUSf(eigensf,1,ttsf,gensf)/USf(eigensf,1,ttsf,gensf)) + EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*CB)*UCha(eigenchi-4,2)*USf(eigensf,2,ttsf,gensf)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dUSf(eigensf,2,ttsf,gensf)/USf(eigensf,2,ttsf,gensf)) + 
     &                  AR(eigenchi,1,ttsf,gensf)*1/2d0*dZSf(1,eigensf,ttsf,gensf) +
     &                  AR(eigenchi,2,ttsf,gensf)*1/2d0*dZSf(2,eigensf,ttsf,gensf) +
     &                  AR(eigenchi,eigensf,ttsf,gensf)*1/2d0*dconjg(dZfL(ttsf-1,gensf))
 	         	
 	         	endif
	          enddo
 	      enddo
	   enddo
	enddo


! Z/W/Photon/Gluon - Fermion - Fermion

        do typex = 1, 4
          do ttf = 1, 4
            do genf = 1, 3
            
              if (typex.eq.1) then
                                 
                 dBL(typex,ttf,genf) = B(1,ttf)*((dZfL(ttf,genf) + dconjg(dZfL(ttf,genf)) + dZG)/2d0 + DeltaGs)
                 dBR(typex,ttf,genf) = B(1,ttf)*((dZfL(ttf,genf) + dconjg(dZfL(ttf,genf)) + dZG)/2d0 + DeltaGs)
                 
                 
                 



              elseif (typex.eq.4) then

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
     &		      +0.5d0*CL(2,eigensf,ttsf,gensf,3)*dconjg(dZSf(2,sfmeigenout,ttsf,gensf))
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


 ! Gluon - Sfermion - Sfermion


        do sfmeigenout = 1, 2
          do eigensf = 1, 2
            do ttsf = 1, 4
              do gensf = 1, 3

                  
                  dCL(sfmeigenout,eigensf,ttsf,gensf,1) = -GS*((dconjg(dZSf(eigensf,sfmeigenout,ttsf,gensf)) + dZSf(sfmeigenout,eigensf,ttsf,gensf) + dZG*Krondelta(sfmeigenout,eigensf))/2d0 + DeltaGs*Krondelta(sfmeigenout,eigensf))
                  
                        
             
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
        

 
! 
! ! ! ===== Z-f-f =====
! ! 	dBL = BL*(1/2d0*dZfL(ttf,genf) + 1/2d0*dconjg(dZfL(ttf,genf)))
! ! 	dBR = BR*(1/2d0*dZfR(ttf,genf) + 1/2d0*dconjg(dZfR(ttf,genf)))
! 
! 
! ===== H0-f-f -- F-coupling =====

! e.g.  FL/R(xtype,pstype,psgen) 
	      

        do ttf = 1, 4
           do genf = 1, 3
           
! --- incoming up-type fermion ---           
           
              if ((ttf.eq.1).or.(ttf.eq.3)) then
! h              
!	     	  call myBget(m4**2,Mfos(pstype,psgen)**2,m3**2)

                  dFL(5,ttf,genf) = -EL*Mf(ttf,genf)*CA/(SW*2d0*MW*SB)*(dMf(ttf,genf)/Mf(ttf,genf) + 1/2d0*dZfR(ttf,genf) + 1/2d0*dconjg(dZfL(ttf,genf)))
                  dFR(5,ttf,genf) = -EL*Mf(ttf,genf)*CA/(SW*2d0*MW*SB)*(dMf(ttf,genf)/Mf(ttf,genf) + 1/2d0*dZfL(ttf,genf) + 1/2d0*dconjg(dZfR(ttf,genf)))
! H                  
                  dFL(6,ttf,genf) = -EL*Mf(ttf,genf)*SA/(SW*2d0*MW*SB)*(dMf(ttf,genf)/Mf(ttf,genf) + 1/2d0*dZfR(ttf,genf) + 1/2d0*dconjg(dZfL(ttf,genf)))
                  dFR(6,ttf,genf) = -EL*Mf(ttf,genf)*SA/(SW*2d0*MW*SB)*(dMf(ttf,genf)/Mf(ttf,genf) + 1/2d0*dZfL(ttf,genf) + 1/2d0*dconjg(dZfR(ttf,genf)))
! H3                  
                  dFL(7,ttf,genf) = -cI*EL*Mf(ttf,genf)*CB/(SW*2d0*MW*SB)*(dMf(ttf,genf)/Mf(ttf,genf) + 1/2d0*dZfR(ttf,genf) + 1/2d0*dconjg(dZfL(ttf,genf)))
                  dFR(7,ttf,genf) =  cI*EL*Mf(ttf,genf)*CB/(SW*2d0*MW*SB)*(dMf(ttf,genf)/Mf(ttf,genf) + 1/2d0*dZfL(ttf,genf) + 1/2d0*dconjg(dZfR(ttf,genf)))
! H+                 
		  dFL(8,ttf,genf) =  EL*Mf(ttf+1,genf)*SB/(sqrt2*SW*MW*CB)*(dMf(ttf+1,genf)/Mf(ttf+1,genf) + 1/2d0*dZfR(ttf,genf) + 1/2d0*dconjg(dZfL(ttf+1,genf)))                
		  dFR(8,ttf,genf) =  EL*Mf(ttf,genf)*CB/(sqrt2*SW*MW*SB)*(dMf(ttf,genf)/Mf(ttf,genf) + 1/2d0*dZfL(ttf,genf) + 1/2d0*dconjg(dZfR(ttf+1,genf)))
 
! --- incoming down-type fermion --- 

 	             else
! h 	
 	          dFL(5,ttf,genf) = EL*Mf(ttf,genf)*SA/(SW*2d0*MW*CB)*(dMf(ttf,genf)/Mf(ttf,genf) + 1/2d0*dZfR(ttf,genf) + 1/2d0*dconjg(dZfL(ttf,genf)))   
 	          dFR(5,ttf,genf) = EL*Mf(ttf,genf)*SA/(SW*2d0*MW*CB)*(dMf(ttf,genf)/Mf(ttf,genf) + 1/2d0*dZfL(ttf,genf) + 1/2d0*dconjg(dZfR(ttf,genf)))
! H 	          
 	          dFL(6,ttf,genf) = -EL*Mf(ttf,genf)*CA/(SW*2d0*MW*CB)*(dMf(ttf,genf)/Mf(ttf,genf) + 1/2d0*dZfR(ttf,genf) + 1/2d0*dconjg(dZfL(ttf,genf))) 
 	          dFR(6,ttf,genf) = -EL*Mf(ttf,genf)*CA/(SW*2d0*MW*CB)*(dMf(ttf,genf)/Mf(ttf,genf) + 1/2d0*dZfL(ttf,genf) + 1/2d0*dconjg(dZfR(ttf,genf))) 
! H3 	          
 	          dFL(7,ttf,genf) = -cI*EL*Mf(ttf,genf)*SB/(SW*2d0*MW*CB)*(dMf(ttf,genf)/Mf(ttf,genf) + 1/2d0*dZfR(ttf,genf) + 1/2d0*dconjg(dZfL(ttf,genf))) 
 	          dFR(7,ttf,genf) =  cI*EL*Mf(ttf,genf)*SB/(SW*2d0*MW*CB)*(dMf(ttf,genf)/Mf(ttf,genf) + 1/2d0*dZfL(ttf,genf) + 1/2d0*dconjg(dZfR(ttf,genf)))
! H- 	          
 	          dFL(8,ttf,genf) = EL*Mf(ttf-1,genf)*CB/(sqrt2*SW*MW*SB)*(dMf(ttf-1,genf)/Mf(ttf-1,genf) + 1/2d0*dZfR(ttf,genf) + 1/2d0*dconjg(dZfL(ttf-1,genf)))
		  dFR(8,ttf,genf) = EL*Mf(ttf,genf)*SB/(sqrt2*SW*MW*CB)*(dMf(ttf,genf)/Mf(ttf,genf) + 1/2d0*dZfL(ttf,genf) + 1/2d0*dconjg(dZfR(ttf-1,genf)))
 	           	          
 	      endif
 	   enddo
 	enddo


! ===== H-Sf-Sf -- G-coupling =====

! =====================================================================================================		
! ================================== Vertex Squark - Squark - Higgs ===================================
! =====================================================================================================	

! e.g.  dG(xtype,sfgaugeeigenin,sfgaugeeigenout,sftype(in),sfgen(in)) 



     	do ttsf = 1, 4
           do gensf = 1, 3
           
! --- incoming up-type sfermion --- 

	if ((ttsf.eq.1).or.(ttsf.eq.3)) then 
! h
	dGLR(5,1,1,ttsf,gensf) = (-sqrt2*hf(ttsf,gensf)*Mf(ttsf,gensf)*CA)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(5,2,2,ttsf,gensf) = (-sqrt2*hf(ttsf,gensf)*Mf(ttsf,gensf)*CA)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(5,1,2,ttsf,gensf) = (-hf(ttsf,gensf)/sqrt2*Af(ttsf,gensf)*CA)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) - (hf(ttsf,gensf)/sqrt2*MUE*SA)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(5,2,1,ttsf,gensf) = dGLR(5,1,2,ttsf,gensf)

! H
	dGLR(6,1,1,ttsf,gensf) = (-sqrt2*hf(ttsf,gensf)*Mf(ttsf,gensf)*SA)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(6,2,2,ttsf,gensf) = (-sqrt2*hf(ttsf,gensf)*Mf(ttsf,gensf)*SA)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(6,1,2,ttsf,gensf) = (-hf(ttsf,gensf)/sqrt2*Af(ttsf,gensf)*SA)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) + (hf(ttsf,gensf)/sqrt2*MUE*CA)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(6,2,1,ttsf,gensf) =  dGLR(6,1,2,ttsf,gensf)
	
! H3
	dGLR(7,1,1,ttsf,gensf) = 0d0
	dGLR(7,2,2,ttsf,gensf) = 0d0
	dGLR(7,1,2,ttsf,gensf) = (cI*hf(ttsf,gensf)/sqrt2*Af(ttsf,gensf)*CB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) + (cI*hf(ttsf,gensf)/sqrt2*MUE*SB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(7,2,1,ttsf,gensf) = - dGLR(7,1,2,ttsf,gensf)

! H+
	dGLR(8,1,1,ttsf,gensf) = (hf(ttsf+1,gensf)*Mf(ttsf+1,gensf)*SB)*(2d0*dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf)) + (hf(ttsf,gensf)*Mf(ttsf,gensf)*CB)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(8,2,2,ttsf,gensf) = ((hf(ttsf,gensf)*Mf(ttsf+1,gensf)*CB) + hf(ttsf+1,gensf)*Mf(ttsf,gensf)*SB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))
	dGLR(8,1,2,ttsf,gensf) = (hf(ttsf,gensf)*Af(ttsf,gensf)*CB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) + (hf(ttsf,gensf)*MUE*SB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(8,2,1,ttsf,gensf) = (hf(ttsf+1,gensf)*Af(ttsf+1,gensf)*SB)*(dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf) + dAf(ttsf+1,gensf)/Af(ttsf+1,gensf)) + (hf(ttsf+1,gensf)*MUE*CB)*(dMf(ttsf+1,gensf)/Mf(ttsf+1,gensf))


! --- incoming down-type sfermion --- 

	else  
! h
	dGLR(5,1,1,ttsf,gensf) = (sqrt2*hf(ttsf,gensf)*Mf(ttsf,gensf)*SA)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(5,2,2,ttsf,gensf) = (sqrt2*hf(ttsf,gensf)*Mf(ttsf,gensf)*SA)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(5,1,2,ttsf,gensf) = (hf(ttsf,gensf)/sqrt2*Af(ttsf,gensf)*SA)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) + (hf(ttsf,gensf)/sqrt2*MUE*CA)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(5,2,1,ttsf,gensf) = dGLR(5,1,2,ttsf,gensf)
! H
	dGLR(6,1,1,ttsf,gensf) = (-sqrt2*hf(ttsf,gensf)*Mf(ttsf,gensf)*CA)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(6,2,2,ttsf,gensf) = (-sqrt2*hf(ttsf,gensf)*Mf(ttsf,gensf)*CA)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(6,1,2,ttsf,gensf) = (-hf(ttsf,gensf)/sqrt2*Af(ttsf,gensf)*CA)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) + (hf(ttsf,gensf)/sqrt2*MUE*SA)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(6,2,1,ttsf,gensf) = dGLR(6,1,2,ttsf,gensf)
! H3
	dGLR(7,1,1,ttsf,gensf) = 0d0
	dGLR(7,2,2,ttsf,gensf) = 0d0
	dGLR(7,1,2,ttsf,gensf) = (cI*hf(ttsf,gensf)/sqrt2*Af(ttsf,gensf)*SB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) + (cI*hf(ttsf,gensf)/sqrt2*MUE*CB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(7,2,1,ttsf,gensf) = -dGLR(7,1,2,ttsf,gensf)

! H-
	dGLR(8,1,1,ttsf,gensf) = (hf(ttsf,gensf)*Mf(ttsf,gensf)*SB)*(2d0*dMf(ttsf,gensf)/Mf(ttsf,gensf)) + (hf(ttsf-1,gensf)*Mf(ttsf-1,gensf)*CB)*(2d0*dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))
	dGLR(8,2,2,ttsf,gensf) = (hf(ttsf-1,gensf)*Mf(ttsf,gensf)*CB + hf(ttsf,gensf)*Mf(ttsf-1,gensf)*SB)*(dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf) + dMf(ttsf,gensf)/Mf(ttsf,gensf))
	
 	dGLR(8,1,2,ttsf,gensf) = (hf(ttsf,gensf)*Af(ttsf,gensf)*SB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf) + dAf(ttsf,gensf)/Af(ttsf,gensf)) + (hf(ttsf,gensf)*MUE*CB)*(dMf(ttsf,gensf)/Mf(ttsf,gensf))
	dGLR(8,2,1,ttsf,gensf) = (hf(ttsf-1,gensf)*Af(ttsf-1,gensf)*CB)*(dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf) + dAf(ttsf-1,gensf)/Af(ttsf-1,gensf)) + (hf(ttsf-1,gensf)*MUE*SB)*(dMf(ttsf-1,gensf)/Mf(ttsf-1,gensf))
	

	endif
	
		
! ------------------------------ Vertex Squark - Squark - neutral Higgs ------------------------------

	do typex=5,7
	   do sfmeigenin=1,2
	      do sfmeigenout=1,2

	         dG(typex,sfmeigenin,sfmeigenout,ttsf,gensf) = (0d0,0d0) 

	         do sfgaugeeigenin=1,2
	            do sfgaugeeigenout=1,2

    	 dG(typex,sfmeigenin,sfmeigenout,ttsf,gensf) = dG(typex,sfmeigenin,sfmeigenout,ttsf,gensf) +
     &   dUSf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)*GLR(typex,sfgaugeeigenout,sfgaugeeigenin,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf,gensf) +
     &   USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)*dGLR(typex,sfgaugeeigenout,sfgaugeeigenin,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf,gensf) +
     &   USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)*GLR(typex,sfgaugeeigenout,sfgaugeeigenin,ttsf,gensf)*dUSf(sfmeigenout,sfgaugeeigenout,ttsf,gensf)	               
   	            enddo

		dG(typex,sfmeigenin,sfmeigenout,ttsf,gensf) = dG(typex,sfmeigenin,sfmeigenout,ttsf,gensf) + 
     &   1/2d0*(dconjg(dZSf(sfgaugeeigenin,sfmeigenout,ttsf,gensf))*G(typex,sfmeigenin,sfgaugeeigenin,ttsf,gensf) +               
     &                 dZSf(sfgaugeeigenin,sfmeigenin,ttsf,gensf)*G(typex,sfgaugeeigenin,sfmeigenout,ttsf,gensf))
	         enddo
 	      enddo
	   enddo
	enddo






! ------------------------------ Vertex Squark - Squark - charged Higgs ------------------------------
	
        do sfmeigenin=1,2
	   do sfmeigenout=1,2
	      typex=8
	      dG(typex,sfmeigenin,sfmeigenout,ttsf,gensf) = (0d0,0d0) 
! --- incoming up-type sfermion ---	
	      if ((ttsf.eq.1).or.(ttsf.eq.3)) then
	
	      do sfgaugeeigenin=1,2
	         do sfgaugeeigenout=1,2 
	            dG(typex,sfmeigenin,sfmeigenout,ttsf,gensf) = dG(typex,sfmeigenin,sfmeigenout,ttsf,gensf) +
     & dGLR(typex,sfgaugeeigenout,sfgaugeeigenin,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf+1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf) +
     & GLR(typex,sfgaugeeigenout,sfgaugeeigenin,ttsf,gensf)*dUSf(sfmeigenout,sfgaugeeigenout,ttsf+1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf) +
     & GLR(typex,sfgaugeeigenout,sfgaugeeigenin,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf+1,gensf)*dUSf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)
                   enddo
		    dG(typex,sfmeigenin,sfmeigenout,ttsf,gensf) = dG(typex,sfmeigenin,sfmeigenout,ttsf,gensf) + 
     &   1/2d0*(dconjg(dZSf(sfgaugeeigenin,sfmeigenout,ttsf+1,gensf))*G(typex,sfmeigenin,sfgaugeeigenin,ttsf,gensf) +               
     &                 dZSf(sfgaugeeigenin,sfmeigenin,ttsf,gensf)*G(typex,sfgaugeeigenin,sfmeigenout,ttsf,gensf))
	         enddo

! --- incoming down-type sfermion ---
                    
                 else
                 do sfgaugeeigenin=1,2
	         do sfgaugeeigenout=1,2 
                    dG(typex,sfmeigenin,sfmeigenout,ttsf,gensf) = dG(typex,sfmeigenin,sfmeigenout,ttsf,gensf) +
     & dGLR(typex,sfgaugeeigenout,sfgaugeeigenin,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf-1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf) +
     & GLR(typex,sfgaugeeigenout,sfgaugeeigenin,ttsf,gensf)*dUSf(sfmeigenout,sfgaugeeigenout,ttsf-1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf) +
     & GLR(typex,sfgaugeeigenout,sfgaugeeigenin,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf-1,gensf)*dUSf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)
                enddo
		    dG(typex,sfmeigenin,sfmeigenout,ttsf,gensf) = dG(typex,sfmeigenin,sfmeigenout,ttsf,gensf) + 
     &   1/2d0*(dconjg(dZSf(sfgaugeeigenin,sfmeigenout,ttsf-1,gensf))*G(typex,sfmeigenin,sfgaugeeigenin,ttsf,gensf) +               
     &                 dZSf(sfgaugeeigenin,sfmeigenin,ttsf,gensf)*G(typex,sfgaugeeigenin,sfmeigenout,ttsf,gensf))
 		
	         enddo
                 endif

	   enddo
	enddo
	
	enddo
	enddo

 	end

! =======================================
! Propagator counterterms 
! =======================================

! Squarks
	subroutine NeuQ2qx_CalcPropSquark(myP2,type,gen)

#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer type,gen,i,j
	double precision myP2, Krondelta

	do i = 1, 2
	  do j = 1, 2
c                  myP2 = 10000d0
	    Prop(i,j,type,gen) = 1/2d0*dZSf(i,j,type,gen)*(myP2 - MSf2(i,type,gen)) + 1/2d0*dZSf(j,i,type,gen)*(myP2 - MSf2(j,type,gen)) - Krondelta(i,j)*dMSf2(i,type,gen)


	  enddo
 	enddo

	end


! Quarks
	subroutine NeuQ2qx_CalcPropQuark(type,gen)

#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_GenCounterterms.h"
#include "../util/DMNLO_LoopIntegrals.h"

	integer type,gen

! PropVR/L according to PL/PR before kslash...
! dZfl, ... etc. follows convention like ChiChi2QQ code...

	    PropVR(type,gen) = 1/2d0*(dZfL(type,gen) + dconjg(dZfL(type,gen)))
	    PropVL(type,gen) = 1/2d0*(dZfR(type,gen) + dconjg(dZfR(type,gen)))
	    PropSL(type,gen) = (-1)*1/2d0*(Mf(type,gen)*dZfL(type,gen) + Mf(type,gen)*dconjg(dZfR(type,gen)))- dMf(type,gen)
	    PropSR(type,gen) = (-1)*1/2d0*(Mf(type,gen)*dZfR(type,gen) + Mf(type,gen)*dconjg(dZfL(type,gen)))- dMf(type,gen)

	end