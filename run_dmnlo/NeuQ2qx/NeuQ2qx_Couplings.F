        subroutine NeuQ2qx_SetCouplings
	implicit none

#include "NeuQ2qx_Model.h"
#include "NeuQ2qx_Kinematics.h"
#include "NeuQ2qx_GenCouplings.h"
#include "NeuQ2qx_GenCounterterms.h"



	integer eigenchi,eigensf,ttsf,ttsf1,ttsf2,gensf,ttf,genf,typex
	integer sfmeigenout,sfmeigenin,sfmeigenout1,sfmeigenin1,sfmeigenout2,sfmeigenin2,sfgaugeeigenout,sfgaugeeigenin,eigenchiin,eigenchiout
	double precision yukawaswitch
      
      double complex BA(4),BV(4)
	
! functions from DMNLO_Kinematics.F
	double precision Krondelta, ef, I3f, hf
	
! calcHEP often neglects the yukawa-couplings for first and second generation
! especially for the second generation one can see a difference between this code and calcHEP
! for the ligthest higgs you can set "yukawaswitch = 0d0" for setting the yukawa couplings equal to zero

	yukawaswitch = 1d0  
	
	
! ===============================================================================================================	
! ============================== Vertex Neutralino - Sfermion (in) - Fermion (out) ==============================
! ==============================          Chargino - Sfermion (in) - Fermion (out) ==============================
! ===============================================================================================================


!  e.g. AL/R(chimeigen, sfmeigen, sftype, sfgen)
	 
        do eigensf = 1, 2
 	      do ttsf = 1, 4
 	      do gensf = 1, 3
 	      
 	      
! ------------------------------ Vertex Neutralino - Sfermion (in) - Fermion (out) ------------------------------

 	         do eigenchi = 1, 4

! --- incoming up-type sfermion ---

 	               if ((ttsf.eq.1).or.(ttsf.eq.3)) then
 	               
 	               	  AR(eigenchi,eigensf,ttsf,gensf) = - yukawaswitch*EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*SB)*ZNeu(eigenchi,4)*USf(eigensf,2,ttsf,gensf) -
     &                    sqrt2*EL/SW*((ef(ttsf) - I3f(ttsf))*SW/CW*ZNeu(eigenchi,1) +  I3f(ttsf)*ZNeu(eigenchi,2))*USf(eigensf,1,ttsf,gensf)
     
                          AL(eigenchi,eigensf,ttsf,gensf) =  - EL*yukawaswitch*Mf(ttsf,gensf)/(SW*sqrt2*MW*SB)*ZNeuC(eigenchi,4)*USf(eigensf,1,ttsf,gensf) + 
     &                    sqrt2*EL*ef(ttsf)/CW*ZNeuC(eigenchi,1)*USf(eigensf,2,ttsf,gensf)
     
!                           print*,'eigensf',eigensf
!                           print*,'ttsf',ttsf
!                           print*,'gensf',gensf
!                           
!                           print*,'AR(eigenchi,eigensf,ttsf,gensf)',AR(eigenchi,eigensf,ttsf,gensf)
!                           print*,'AL(eigenchi,eigensf,ttsf,gensf)',AL(eigenchi,eigensf,ttsf,gensf)
!                           
!                           print*,'R1',(- yukawaswitch*EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*SB)*ZNeu(eigenchi,4)*USf(eigensf,2,ttsf,gensf))
!                           print*,'R2',(-
!      &                    sqrt2*EL/SW*((ef(ttsf) - I3f(ttsf))*SW/CW*ZNeu(eigenchi,1) +  I3f(ttsf)*ZNeu(eigenchi,2))*USf(eigensf,1,ttsf,gensf))
!                           
!                           print*,'L1',(- EL*yukawaswitch*Mf(ttsf,gensf)/(SW*sqrt2*MW*SB)*ZNeuC(eigenchi,4)*USf(eigensf,1,ttsf,gensf))
!                           print*,'L2',(sqrt2*EL*ef(ttsf)/CW*ZNeuC(eigenchi,1)*USf(eigensf,2,ttsf,gensf))
     
! debugging version without mixing matrices     
!                           AR(eigenchi,eigensf,ttsf,gensf) = - EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*SB) -
!      &                    sqrt2*EL/SW*((ef(ttsf) - I3f(ttsf))*SW/CW +  I3f(ttsf))
!      
!                           AL(eigenchi,eigensf,ttsf,gensf) =  - EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*SB) + 
!      &                    sqrt2*EL*ef(ttsf)/CW
     

! --- incoming down-type sfermion---     

      		       else
      		       
                           AR(eigenchi,eigensf,ttsf,gensf) =  - EL*yukawaswitch*Mf(ttsf,gensf)/(SW*sqrt2*MW*CB)*ZNeu(eigenchi,3)*USf(eigensf,2,ttsf,gensf) -
     &                     sqrt2*EL/SW*(I3f(ttsf)*ZNeu(eigenchi,2) + (ef(ttsf) - I3f(ttsf))*SW/CW*ZNeu(eigenchi,1))*USf(eigensf,1,ttsf,gensf)

                           AL(eigenchi,eigensf,ttsf,gensf) =  - EL*yukawaswitch*Mf(ttsf,gensf)/(SW*sqrt2*MW*CB)*ZNeuC(eigenchi,3)*USf(eigensf,1,ttsf,gensf) + 
     &                     sqrt2*EL*ef(ttsf)/CW*ZNeuC(eigenchi,1)*USf(eigensf,2,ttsf,gensf)
	               endif
	          enddo

	          
! ------------------------------ Vertex Chargino - Sfermion (in) - Fermion (out) ------------------------------

	          do eigenchi= 5, 6   

! --- incoming up-type sfermion ---	          

	          if ((ttsf.eq.1).or.(ttsf.eq.3)) then

 	         	KL(eigenchi,eigensf,ttsf,gensf) =  -EL/SW*UChaC(eigenchi-4,1)*USf(eigensf,1,ttsf+1,gensf) + EL*Mf(ttsf+1,gensf)/(SW*sqrt2*MW*CB)*UChaC(eigenchi-4,2)*USf(eigensf,2,ttsf+1,gensf)
 	         	KR(eigenchi,eigensf,ttsf,gensf) =  EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*SB)*VCha(eigenchi-4,2)*USf(eigensf,1,ttsf+1,gensf)
 	         	
	         	AL(eigenchi,eigensf,ttsf,gensf) =  EL*Mf(ttsf+1,gensf)/(SW*sqrt2*MW*CB)*UChaC(eigenchi-4,2)*USf(eigensf,1,ttsf,gensf)  
 	         	AR(eigenchi,eigensf,ttsf,gensf) = -EL/SW*VCha(eigenchi-4,1)*USf(eigensf,1,ttsf,gensf) + EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*SB)*VCha(eigenchi-4,2)*USf(eigensf,2,ttsf,gensf) 

! --- incoming down-type sfermion --- 	 

 	         	else

 	         	KL(eigenchi,eigensf,ttsf,gensf) =  -EL/SW*VChaC(eigenchi-4,1)*USf(eigensf,1,ttsf-1,gensf) + EL*Mf(ttsf-1,gensf)/(SW*sqrt2*MW*SB)*VChaC(eigenchi-4,2)*USf(eigensf,2,ttsf-1,gensf)
 	         	KR(eigenchi,eigensf,ttsf,gensf) =  EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*CB)*UCha(eigenchi-4,2)*USf(eigensf,1,ttsf-1,gensf)

 	         	AL(eigenchi,eigensf,ttsf,gensf) = EL*Mf(ttsf-1,gensf)/(SW*sqrt2*MW*SB)*VChaC(eigenchi-4,2)*USf(eigensf,1,ttsf,gensf) 
 	         	AR(eigenchi,eigensf,ttsf,gensf) = -EL/SW*UCha(eigenchi-4,1)*USf(eigensf,1,ttsf,gensf) + EL*Mf(ttsf,gensf)/(SW*sqrt2*MW*CB)*UCha(eigenchi-4,2)*USf(eigensf,2,ttsf,gensf) 
 	         	
 	         	endif
	          enddo

 	      enddo
	      enddo
	      enddo
	 

       


! ===================================================================================================
! ============================== Vertex Fermion - Fermion - Gluon / Photon ==========================
! ===================================================================================================

! e.g.  B(xtype,pstype) 

	B(1,1) =  0d0
	B(1,2) =  0d0

	B(1,3) = -GS
	B(1,4) = -GS

	B(2,1) = -ef(1)*EL
	B(2,2) = -ef(2)*EL
	B(2,3) = -ef(3)*EL
	B(2,4) = -ef(4)*EL
	
	do typex = 1, 2
	        do ttf = 1, 4
	                BL(typex,ttf) = B(typex,ttf)
                    BR(typex,ttf) = B(typex,ttf)
                enddo
        enddo
	
! ===================================================================================================
! =================================== Vertex Fermion - Fermion - W / Z ==============================
! ===================================================================================================

        do ttf = 1, 4

                BV(ttf) = I3f(ttf) - 2d0*SW**2d0*ef(ttf)
                BA(ttf) = I3f(ttf)

                BL(3,ttf) = - (g2h1/(2d0*CW))*(BV(ttf) + BA(ttf))
                BR(3,ttf) = - (g2h1/(2d0*CW))*(BV(ttf) - BA(ttf))

                BL(4,ttf) = - g2h1/sqrt2
                BR(4,ttf) = 0d0

	enddo


	
	
! =====================================================================================================	
! ============================== Vertex Sfermion - Sfermion - Gluon / Photon ==========================
! =====================================================================================================

! e.g.  C(xtype,pstype) 


	C(1,1) =  0d0
	C(1,2) =  0d0

	C(1,3) = -GS
	C(1,4) = -GS

	C(2,1) = -ef(1)*EL
	C(2,2) = -ef(2)*EL
	C(2,3) = -ef(3)*EL
	C(2,4) = -ef(4)*EL

        do sfmeigenout = 1, 2
          do eigensf = 1, 2   
            do ttsf = 1, 4
              do gensf = 1, 3
                       CL(sfmeigenout,eigensf,ttsf,gensf,1) = C(1,ttsf)*Krondelta(eigensf,sfmeigenout)
                       CL(sfmeigenout,eigensf,ttsf,gensf,2) = C(2,ttsf)*Krondelta(eigensf,sfmeigenout)
              enddo
       	    enddo
       	  enddo
      	enddo
  
        



! ===================================================================================================     
! ====================================== Sfermion - Sfermion - W / Z ================================
! ===================================================================================================



        do sfmeigenout = 1, 2
          do eigensf = 1, 2   
            do ttsf = 1, 4
              do gensf = 1, 3                                 
     
                CL(sfmeigenout,eigensf,ttsf,gensf,3) = - (g2h1/CW)*(I3f(ttsf)*USf(eigensf,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf,gensf) -
     &          SW**2*ef(ttsf)*Krondelta(eigensf,sfmeigenout))
     
! --- incoming up-type sfermion ---
               
                if ((ttsf.eq.1).or.(ttsf.eq.3)) then

    	                CL(sfmeigenout,eigensf,ttsf,gensf,4) = - g2h1/(sqrt2)*USf(eigensf,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf+1,gensf)

! --- incoming down-type sfermion ---
               
                else
    	       
    	                CL(sfmeigenout,eigensf,ttsf,gensf,4) = - g2h1/(sqrt2)*USf(eigensf,1,ttsf,gensf)*USf(sfmeigenout,1,ttsf-1,gensf)
                                
                endif
                   
                enddo
 	      enddo
 	    enddo
	  enddo


! ===================================================================================================     
! ====================================== Sfermion - Sfermion - Gluon - Vector ============================
! ===================================================================================================


! Z
        do sfmeigenout = 1, 2
          do eigensf = 1, 2   
            do ttsf = 1, 4
              do gensf = 1, 3                                 
     
                CLL(sfmeigenout,eigensf,ttsf,gensf,3) = (2*g2h1*GS/CW)*(I3f(ttsf)*USf(sfmeigenout,1,ttsf,gensf)*USf(eigensf,1,ttsf,gensf)-ef(ttsf)*Krondelta(eigensf,sfmeigenout)*SW**2)
                   
          enddo
 	      enddo
 	    enddo
	  enddo


! Photon
        do sfmeigenout = 1, 2
          do eigensf = 1, 2   
            do ttsf = 1, 4
              do gensf = 1, 3                                 
     
        			CLL(sfmeigenout,eigensf,ttsf,gensf,2) = 2*ef(ttsf)*EL*GS*Krondelta(eigensf,sfmeigenout)
                   
          enddo
 	      enddo
 	    enddo
	  enddo

! W
        do sfmeigenout = 1, 2
          do eigensf = 1, 2   
            do ttsf = 1, 4
              do gensf = 1, 3                                 
     
! --- incoming up-type sfermion ---
              
               if ((ttsf.eq.1).or.(ttsf.eq.3)) then
    
                CLL(sfmeigenout,eigensf,ttsf,gensf,4) = g2h1*GS*sqrt2*USf(sfmeigenout,1,ttsf+1,gensf)*USf(eigensf,1,ttsf,gensf)
           
! --- incoming down-type sfermion ---
              
               else

                CLL(sfmeigenout,eigensf,ttsf,gensf,4) = g2h1*GS*sqrt2*USf(sfmeigenout,1,ttsf-1,gensf)*USf(eigensf,1,ttsf,gensf)

               endif
                   
          enddo
 	      enddo
 	    enddo
	  enddo


! Gluon
            do sfmeigenout = 1, 2
              do eigensf = 1, 2   
                do ttsf = 1, 4
                  do gensf = 1, 3                                 
             
                        CLL(sfmeigenout,eigensf,ttsf,gensf,1) = GS**2*Krondelta(eigensf,sfmeigenout)
                           
                  enddo
                enddo
              enddo
            enddo


! ===================================================================================================     
! ========= Vertex Neutralino (in) - Virtual Neutralino/Chargino (out) - Vector Boson (out) =========
! ===================================================================================================     
                
                
        do eigenchi = 1, 4
          do eigenchiout = 1, 4
            do typex = 1, 3
              do ttsf = 1, 3
                              
        if ((typex.eq.1).OR.(typex.eq.2)) then
              
              DL(eigenchi,eigenchiout,ttsf,typex) = 0d0
              DR(eigenchi,eigenchiout,ttsf,typex) = 0d0
   
        else if (typex.eq.3) then
              
              DL(eigenchi,eigenchiout,ttsf,typex) = (g2h1/(2d0*CW))*( - ZNeu(eigenchiout,3)*ZNeuC(eigenchi,3) + ZNeu(eigenchiout,4)*ZNeuC(eigenchi,4) )
              DR(eigenchi,eigenchiout,ttsf,typex) = (g2h1/(2d0*CW))*(   ZNeuC(eigenchiout,3)*ZNeu(eigenchi,3) - ZNeuC(eigenchiout,4)*ZNeu(eigenchi,4) )

        endif

              enddo
            enddo
          enddo
        enddo

        do eigenchi = 1, 4
          do eigenchiout = 5, 6
            do ttsf = 1, 3

! --- incoming up-type sfermion ---
               
                if ((ttsf.eq.1).or.(ttsf.eq.3)) then

                        DR(eigenchi,eigenchiout,ttsf,4) = - g2h1*( ZNeu(eigenchi,2)*VChaC(eigenchiout-4,1) - (1/sqrt2)*ZNeu(eigenchi,4)*VChaC(eigenchiout-4,2) )
                        DL(eigenchi,eigenchiout,ttsf,4) = - g2h1*( ZNeuC(eigenchi,2)*UCha(eigenchiout-4,1) + (1/sqrt2)*ZNeuC(eigenchi,3)*UCha(eigenchiout-4,2) )

! --- incoming down-type sfermion ---
               
                else

c !!!! apparently not consistent with Feynman rules but gives good agreement with CalcHEP !!!! -> To be checked

                        DL(eigenchi,eigenchiout,ttsf,4) = g2h1*( ZNeu(eigenchi,2)*VChaC(eigenchiout-4,1) - (1/sqrt2)*ZNeu(eigenchi,4)*VChaC(eigenchiout-4,2) )
                        DR(eigenchi,eigenchiout,ttsf,4) = g2h1*( ZNeuC(eigenchi,2)*UCha(eigenchiout-4,1) + (1/sqrt2)*ZNeuC(eigenchi,3)*UCha(eigenchiout-4,2) )
              
                endif
            
            enddo
          enddo
        enddo




! =====================================================================================================	
! ==================================== Vertex Quark - Quark - Higgs ===================================
! =====================================================================================================

! e.g.  FL/R(xtype,pstype,psgen) 
	      

        do ttf = 1, 4
           do genf = 1, 3
           
! --- incoming up-type fermion ---           
           
              if ((ttf.eq.1).or.(ttf.eq.3)) then
! h              
                  FL(5,ttf,genf) = -yukawaswitch*EL*Mf(ttf,genf)*CA/(SW*2d0*MW*SB) 
                  FR(5,ttf,genf) =  FL(5,ttf,genf)
! H                  
                  FL(6,ttf,genf) = -EL*Mf(ttf,genf)*SA/(SW*2d0*MW*SB) 
                  FR(6,ttf,genf) =  FL(6,ttf,genf)
! H3                  
                  FL(7,ttf,genf) = -cI*EL*Mf(ttf,genf)*CB/(SW*2d0*MW*SB) 
                  FR(7,ttf,genf) = -FL(7,ttf,genf)
! H+                 
		  FL(8,ttf,genf) =  EL*Mf(ttf+1,genf)*SB/(sqrt2*SW*MW*CB)                  
		  FR(8,ttf,genf) =  EL*Mf(ttf,genf)*CB/(sqrt2*SW*MW*SB) 
 
! --- incoming down-type fermion --- 

 	             else
! h 	
 	          FL(5,ttf,genf) = yukawaswitch*EL*Mf(ttf,genf)*SA/(SW*2d0*MW*CB)   
 	          FR(5,ttf,genf) = FL(5,ttf,genf)
! H 	          
 	          FL(6,ttf,genf) = -EL*Mf(ttf,genf)*CA/(SW*2d0*MW*CB) 
 	          FR(6,ttf,genf) =  FL(6,ttf,genf)
! H3 	          
 	          FL(7,ttf,genf) = -cI*EL*Mf(ttf,genf)*SB/(SW*2d0*MW*CB) 
 	          FR(7,ttf,genf) = -FL(7,ttf,genf)
! H- 	          
 	          FL(8,ttf,genf) = EL*Mf(ttf-1,genf)*CB/(sqrt2*SW*MW*SB)
		  FR(8,ttf,genf) = EL*Mf(ttf,genf)*SB/(sqrt2*SW*MW*CB) 
 	           	          
 	      endif
 	   enddo
 	enddo

 	
 	
 	
! =====================================================================================================		
! ================================== Vertex Squark - Squark - Higgs ===================================
! =====================================================================================================	

! e.g.  G(xtype,sfgaugeeigenin,sfgaugeeigenout,sftype(in),sfgen(in)) 


     	do ttsf = 1, 4
           do gensf = 1, 3
           
! --- incoming up-type sfermion --- 

	if ((ttsf.eq.1).or.(ttsf.eq.3)) then 
! h
	GLR(5,1,1,ttsf,gensf) = -yukawaswitch*sqrt2*hf(ttsf,gensf)*Mf(ttsf,gensf)*CA + EL*MZ*(I3f(ttsf)-ef(ttsf)*SW2)*SAB/(SW*CW) 
	GLR(5,2,2,ttsf,gensf) = -yukawaswitch*sqrt2*hf(ttsf,gensf)*Mf(ttsf,gensf)*CA + EL*MZ*ef(ttsf)*SW2*SAB/(SW*CW) 
	GLR(5,1,2,ttsf,gensf) = -yukawaswitch*hf(ttsf,gensf)/sqrt2*(Af(ttsf,gensf)*CA + MUE*SA) 
	GLR(5,2,1,ttsf,gensf) = GLR(5,1,2,ttsf,gensf)
! H
	GLR(6,1,1,ttsf,gensf) = -sqrt2*hf(ttsf,gensf)*Mf(ttsf,gensf)*SA - EL*MZ*(I3f(ttsf)-ef(ttsf)*SW2)*CAB/(SW*CW) 
	GLR(6,2,2,ttsf,gensf) = -sqrt2*hf(ttsf,gensf)*Mf(ttsf,gensf)*SA - EL*MZ*ef(ttsf)*SW2*CAB/(SW*CW) 
	GLR(6,1,2,ttsf,gensf) = -hf(ttsf,gensf)/sqrt2*(Af(ttsf,gensf)*SA - MUE*CA) 
	GLR(6,2,1,ttsf,gensf) =  GLR(6,1,2,ttsf,gensf)
	
! H3
	GLR(7,1,1,ttsf,gensf) = 0d0
	GLR(7,2,2,ttsf,gensf) = 0d0
	GLR(7,1,2,ttsf,gensf) = cI*hf(ttsf,gensf)/sqrt2*(Af(ttsf,gensf)*CB + MUE*SB) 
	GLR(7,2,1,ttsf,gensf) = - GLR(7,1,2,ttsf,gensf)

! H+
	GLR(8,1,1,ttsf,gensf) = hf(ttsf+1,gensf)*Mf(ttsf+1,gensf)*SB + hf(ttsf,gensf)*Mf(ttsf,gensf)*CB - EL*MW/(sqrt2*SW)*S2B
	GLR(8,2,2,ttsf,gensf) = hf(ttsf,gensf)*Mf(ttsf+1,gensf)*CB + hf(ttsf+1,gensf)*Mf(ttsf,gensf)*SB
	GLR(8,1,2,ttsf,gensf) = hf(ttsf,gensf)*(Af(ttsf,gensf)*CB + MUE*SB) 
	GLR(8,2,1,ttsf,gensf) = hf(ttsf+1,gensf)*(Af(ttsf+1,gensf)*SB + MUE*CB)


! --- incoming down-type sfermion --- 

	else  
! h
	GLR(5,1,1,ttsf,gensf) = yukawaswitch*sqrt2*hf(ttsf,gensf)*Mf(ttsf,gensf)*SA + EL*MZ*(I3f(ttsf)-ef(ttsf)*SW2)*SAB/(SW*CW) 
	GLR(5,2,2,ttsf,gensf) = yukawaswitch*sqrt2*hf(ttsf,gensf)*Mf(ttsf,gensf)*SA + EL*MZ*ef(ttsf)*SW2*SAB/(SW*CW) 
	GLR(5,1,2,ttsf,gensf) = yukawaswitch*hf(ttsf,gensf)/sqrt2*(Af(ttsf,gensf)*SA + MUE*CA) 
	GLR(5,2,1,ttsf,gensf) = GLR(5,1,2,ttsf,gensf)
! H
	GLR(6,1,1,ttsf,gensf) = -sqrt2*hf(ttsf,gensf)*Mf(ttsf,gensf)*CA - EL*MZ*(I3f(ttsf)-ef(ttsf)*SW2)*CAB/(SW*CW) 
	GLR(6,2,2,ttsf,gensf) = -sqrt2*hf(ttsf,gensf)*Mf(ttsf,gensf)*CA - EL*MZ*ef(ttsf)*SW2*CAB/(SW*CW) 
	GLR(6,1,2,ttsf,gensf) = -hf(ttsf,gensf)/sqrt2*(Af(ttsf,gensf)*CA - MUE*SA) 
	GLR(6,2,1,ttsf,gensf) = GLR(6,1,2,ttsf,gensf)
! H3
	GLR(7,1,1,ttsf,gensf) = 0d0
	GLR(7,2,2,ttsf,gensf) = 0d0
	GLR(7,1,2,ttsf,gensf) = cI*hf(ttsf,gensf)/sqrt2*(Af(ttsf,gensf)*SB + MUE*CB) 
	GLR(7,2,1,ttsf,gensf) = -GLR(7,1,2,ttsf,gensf)
! H-
	GLR(8,1,1,ttsf,gensf) = hf(ttsf,gensf)*Mf(ttsf,gensf)*SB + hf(ttsf-1,gensf)*Mf(ttsf-1,gensf)*CB - EL*MW/(sqrt2*SW)*S2B
	GLR(8,2,2,ttsf,gensf) = hf(ttsf-1,gensf)*Mf(ttsf,gensf)*CB + hf(ttsf,gensf)*Mf(ttsf-1,gensf)*SB
	
 	GLR(8,1,2,ttsf,gensf) = hf(ttsf,gensf)*(Af(ttsf,gensf)*SB + MUE*CB) 
	GLR(8,2,1,ttsf,gensf) = hf(ttsf-1,gensf)*(Af(ttsf-1,gensf)*CB + MUE*SB)
	

	endif
	
		
! ------------------------------ Vertex Squark - Squark - neutral Higgs ------------------------------

	do typex=5,7
	   do sfmeigenin=1,2
	      do sfmeigenout=1,2

	         G(typex,sfmeigenin,sfmeigenout,ttsf,gensf) = (0d0,0d0) 

	         do sfgaugeeigenin=1,2
	            do sfgaugeeigenout=1,2
	            	
    	               G(typex,sfmeigenin,sfmeigenout,ttsf,gensf) = G(typex,sfmeigenin,sfmeigenout,ttsf,gensf) + GLR(typex,sfgaugeeigenout,sfgaugeeigenin,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)
	               
   	            enddo
	         enddo
 	      enddo
	   enddo
	enddo


	
! ------------------------------ Vertex Squark - Squark - charged Higgs ------------------------------
	
        do sfmeigenin=1,2
	   do sfmeigenout=1,2
	      typex=8
	      G(typex,sfmeigenin,sfmeigenout,ttsf,gensf) = (0d0,0d0) 
	
	      do sfgaugeeigenin=1,2
	         do sfgaugeeigenout=1,2

! --- incoming up-type sfermion ---	         
	         
	            if ((ttsf.eq.1).or.(ttsf.eq.3)) then
	            
	            G(typex,sfmeigenin,sfmeigenout,ttsf,gensf) = G(typex,sfmeigenin,sfmeigenout,ttsf,gensf) + GLR(typex,sfgaugeeigenout,sfgaugeeigenin,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf+1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)

! --- incoming down-type sfermion ---
                    
                    else
                    
                    G(typex,sfmeigenin,sfmeigenout,ttsf,gensf) = G(typex,sfmeigenin,sfmeigenout,ttsf,gensf) + GLR(typex,sfgaugeeigenout,sfgaugeeigenin,ttsf,gensf)*USf(sfmeigenout,sfgaugeeigenout,ttsf-1,gensf)*USf(sfmeigenin,sfgaugeeigenin,ttsf,gensf)
                    
                    endif
                 enddo
	      enddo
	   enddo
	enddo
	
	enddo
	enddo


	
	
! ==============================================================================================================	
! ================================== Vertizes Neutralino - Neutralino - Higgs ==================================
! ==================================          Neutralino -  Chargino  - Higgs ==================================
! ==================================           Chargino  -  Chargino  - Higgs ==================================
! ==============================================================================================================	

! e.g. HR(chimeigenout,chimeigenin,xtype) 


! ------------------------------ Vertex Neutralino - Neutralino - Higgs ------------------------------

	do eigenchiin=1,4
	   do eigenchiout=1,4
	   
! h
	   	HL(eigenchiin, eigenchiout,5) = -EL/(2d0*SW)*((-SA)*(ZNeuC(eigenchiout,3)*ZNeuC(eigenchiin,2) + ZNeuC(eigenchiin,3)*ZNeuC(eigenchiout,2) - SW/CW*(ZNeuC(eigenchiout,3)*ZNeuC(eigenchiin,1) + ZNeuC(eigenchiin,3)*ZNeuC(eigenchiout,1)))
     &	        + (-CA)*(ZNeuC(eigenchiout,4)*ZNeuC(eigenchiin,2) + ZNeuC(eigenchiin,4)*ZNeuC(eigenchiout,2) - SW/CW*(ZNeuC(eigenchiout,4)*ZNeuC(eigenchiin,1) + ZNeuC(eigenchiin,4)*ZNeuC(eigenchiout,1))) )
		HR(eigenchiin, eigenchiout,5) = -EL/(2d0*SW)*((-SA)*(ZNeu(eigenchiout,3)*ZNeu(eigenchiin,2)  +  ZNeu(eigenchiin,3)*ZNeu(eigenchiout,2)  - SW/CW*(ZNeu(eigenchiout,3)*ZNeu(eigenchiin,1)  +  ZNeu(eigenchiin,3)*ZNeu(eigenchiout,1)))
     &	        + (-CA)*(ZNeu(eigenchiout,4)*ZNeu(eigenchiin,2)  +  ZNeu(eigenchiin,4)*ZNeu(eigenchiout,2)  - SW/CW*(ZNeu(eigenchiout,4)*ZNeu(eigenchiin,1)  +  ZNeu(eigenchiin,4)*ZNeu(eigenchiout,1))) )
! H     
		HL(eigenchiin, eigenchiout,6) = -EL/(2d0*SW)*((+CA)*(ZNeuC(eigenchiout,3)*ZNeuC(eigenchiin,2) + ZNeuC(eigenchiin,3)*ZNeuC(eigenchiout,2) - SW/CW*(ZNeuC(eigenchiout,3)*ZNeuC(eigenchiin,1) + ZNeuC(eigenchiin,3)*ZNeuC(eigenchiout,1)))
     &	        + (-SA)*(ZNeuC(eigenchiout,4)*ZNeuC(eigenchiin,2) + ZNeuC(eigenchiin,4)*ZNeuC(eigenchiout,2) - SW/CW*(ZNeuC(eigenchiout,4)*ZNeuC(eigenchiin,1) + ZNeuC(eigenchiin,4)*ZNeuC(eigenchiout,1))) )
		HR(eigenchiin, eigenchiout,6) = -EL/(2d0*SW)*((+CA)*(ZNeu(eigenchiout,3)*ZNeu(eigenchiin,2)  +  ZNeu(eigenchiin,3)*ZNeu(eigenchiout,2)  - SW/CW*(ZNeu(eigenchiout,3)*ZNeu(eigenchiin,1)  +  ZNeu(eigenchiin,3)*ZNeu(eigenchiout,1)))
     &	        + (-SA)*(ZNeu(eigenchiout,4)*ZNeu(eigenchiin,2)  +  ZNeu(eigenchiin,4)*ZNeu(eigenchiout,2)  - SW/CW*(ZNeu(eigenchiout,4)*ZNeu(eigenchiin,1)  +  ZNeu(eigenchiin,4)*ZNeu(eigenchiout,1))) )
! H3     
		HL(eigenchiin, eigenchiout,7) = -cI*EL/(2d0*SW)*((-SB)*(ZNeuC(eigenchiout,3)*ZNeuC(eigenchiin,2) + ZNeuC(eigenchiin,3)*ZNeuC(eigenchiout,2) - SW/CW*(ZNeuC(eigenchiout,3)*ZNeuC(eigenchiin,1) + ZNeuC(eigenchiin,3)*ZNeuC(eigenchiout,1)))
     &	           + (+CB)*(ZNeuC(eigenchiout,4)*ZNeuC(eigenchiin,2) + ZNeuC(eigenchiin,4)*ZNeuC(eigenchiout,2) - SW/CW*(ZNeuC(eigenchiout,4)*ZNeuC(eigenchiin,1) + ZNeuC(eigenchiin,4)*ZNeuC(eigenchiout,1))) )
		HR(eigenchiin, eigenchiout,7) = +cI*EL/(2d0*SW)*((-SB)*(ZNeu(eigenchiout,3)*ZNeu(eigenchiin,2)  +  ZNeu(eigenchiin,3)*ZNeu(eigenchiout,2)  - SW/CW*(ZNeu(eigenchiout,3)*ZNeu(eigenchiin,1)  +  ZNeu(eigenchiin,3)*ZNeu(eigenchiout,1)))
     &	           + (+CB)*(ZNeu(eigenchiout,4)*ZNeu(eigenchiin,2)  +  ZNeu(eigenchiin,4)*ZNeu(eigenchiout,2)  - SW/CW*(ZNeu(eigenchiout,4)*ZNeu(eigenchiin,1)  +  ZNeu(eigenchiin,4)*ZNeu(eigenchiout,1))) )
! H+		
		HL(eigenchiin, eigenchiout,8) = 0d0
		HR(eigenchiin, eigenchiout,8) = 0d0
		
	
	   enddo
     	enddo

! ------------------------------ Vertex Neutralino - Chargino - Higgs ------------------------------	
  
	do eigenchiin=1,4
	   do eigenchiout=5,6
	   
! h
	      HL(eigenchiin, eigenchiout,5) = 0d0
	      HR(eigenchiin, eigenchiout,5) = 0d0
! H	      
	      HL(eigenchiin, eigenchiout,6) = 0d0	
	      HR(eigenchiin, eigenchiout,6) = 0d0
! H3	      
	      HL(eigenchiin, eigenchiout,7) = 0d0
	      HR(eigenchiin, eigenchiout,7) = 0d0
! H+/-	      
	      if ((sftype.eq.1).or.(sftype.eq.3)) then

	      HL(eigenchiin, eigenchiout,8) = -EL/SW*CB*(VChaC(eigenchiout-4,1)*ZNeuC(eigenchiin,4) + 1/sqrt2*(ZNeuC(eigenchiin,2) + SW/CW*ZNeuC(eigenchiin,1))*VChaC(eigenchiout-4,2)) 
	      HR(eigenchiin, eigenchiout,8) = -EL/SW*SB*(UCha(eigenchiout-4,1)*ZNeu(eigenchiin,3) - 1/sqrt2*(ZNeu(eigenchiin,2) + SW/CW*ZNeu(eigenchiin,1))*UCha(eigenchiout-4,2))
	      
	      else
	      
	      HR(eigenchiin, eigenchiout,8) = -EL/SW*CB*(VChaC(eigenchiout-4,1)*ZNeuC(eigenchiin,4) + 1/sqrt2*(ZNeuC(eigenchiin,2) + SW/CW*ZNeuC(eigenchiin,1))*VChaC(eigenchiout-4,2)) 
	      HL(eigenchiin, eigenchiout,8) = -EL/SW*SB*(UCha(eigenchiout-4,1)*ZNeu(eigenchiin,3) - 1/sqrt2*(ZNeu(eigenchiin,2) + SW/CW*ZNeu(eigenchiin,1))*UCha(eigenchiout-4,2))
	      
	      endif
	     
	   enddo
     	enddo
     	
     	
     	do eigenchiin=5,6
	   do eigenchiout=1,4

! h	   
	      HL(eigenchiin, eigenchiout,5) = 0d0
	      HR(eigenchiin, eigenchiout,5) = 0d0
! H	      
	      HL(eigenchiin, eigenchiout,6) = 0d0	
	      HR(eigenchiin, eigenchiout,6) = 0d0
! H3	      
	      HL(eigenchiin, eigenchiout,7) = 0d0
	      HR(eigenchiin, eigenchiout,7) = 0d0
! H+	     
	      HL(eigenchiin, eigenchiout,8) = -EL/SW*SB*(UCha(eigenchiout,1)*ZNeu(eigenchiin-4,3) - 1/sqrt2*(ZNeu(eigenchiin-4,2) + SW/CW*ZNeu(eigenchiin-4,1))*UCha(eigenchiout,2))
	      HR(eigenchiin, eigenchiout,8) = -EL/SW*CB*(VChaC(eigenchiout,1)*ZNeuC(eigenchiin-4,4) + 1/sqrt2*(ZNeuC(eigenchiin-4,2) + SW/CW*ZNeuC(eigenchiin-4,1))*VChaC(eigenchiout,2)) 

	   enddo
     	enddo
	
! ------------------------------ Vertex Chargino - Chargino - Higgs (H-) ------------------------------
	

	do eigenchiin=5,6
	   do eigenchiout=5,6

! h
	HL(eigenchiin, eigenchiout,5) = -EL/(SW*sqrt2)*(-SA*VChaC(eigenchiout-4,1)*UChaC(eigenchiin-4,2) + CA*VChaC(eigenchiout-4,2)*UChaC(eigenchiin-4,1))
	HR(eigenchiin, eigenchiout,5) = -EL/(SW*sqrt2)*(-SA*VCha(eigenchiin-4,1)*UCha(eigenchiout-4,2)   + CA*VCha(eigenchiin-4,2)*UCha(eigenchiout-4,1))
! H	
	HL(eigenchiin, eigenchiout,6) = -EL/(SW*sqrt2)*(CA*VChaC(eigenchiout-4,1)*UChaC(eigenchiin-4,2) + SA*VChaC(eigenchiout-4,2)*UChaC(eigenchiin-4,1))
	HR(eigenchiin, eigenchiout,6) = -EL/(SW*sqrt2)*(CA*VCha(eigenchiin-4,1)*UCha(eigenchiout-4,2) + SA*VCha(eigenchiin-4,2)*UCha(eigenchiout-4,1))
! H3	
	HL(eigenchiin, eigenchiout,7) = -cI*EL/(SW*sqrt2)*(-SB*VChaC(eigenchiout-4,1)*UChaC(eigenchiin-4,2) - CB*VChaC(eigenchiout-4,2)*UChaC(eigenchiin-4,1))
	HR(eigenchiin, eigenchiout,7) = cI*EL/(SW*sqrt2)*(-SB*VCha(eigenchiin-4,1)*UCha(eigenchiout-4,2) - CB*VCha(eigenchiin-4,2)*UCha(eigenchiout-4,1))	
! H-
	HL(eigenchiin, eigenchiout,8) = 0d0
	HR(eigenchiin, eigenchiout,8) = 0d0

	
	   enddo
     	enddo
     	

! ===================================================================================================== 
! =============================== Vertex Squark - Squark - Squark - Squark ============================
! =====================================================================================================
	

	do sfmeigenin1=1,2
	   do sfmeigenout1=1,2
	      do ttsf1 = 1, 4
	         do gensf = 1, 3

          Aij(sfmeigenout1,sfmeigenin1,ttsf1,gensf) = USf(sfmeigenout1,1,ttsf1,gensf)*USf(sfmeigenin1,1,ttsf1,gensf) - USf(sfmeigenout1,2,ttsf1,gensf)*USf(sfmeigenin1,2,ttsf1,gensf)

        	     do sfmeigenin2=1,2
                   do sfmeigenout2=1,2
		         do ttsf2 = 1, 4
         
         S(sfmeigenin1,sfmeigenout1,sfmeigenin2,sfmeigenout2,ttsf1,ttsf2,gensf) = 
     &        -GS**2*(USf(sfmeigenout1,1,ttsf1,gensf)*USf(sfmeigenin1,1,ttsf1,gensf)-USf(sfmeigenout1,2,ttsf1,gensf)*USf(sfmeigenin1,2,ttsf1,gensf))*
     &               (USf(sfmeigenout2,1,ttsf2,gensf)*USf(sfmeigenin2,1,ttsf2,gensf)-USf(sfmeigenout2,2,ttsf2,gensf)*USf(sfmeigenin2,2,ttsf2,gensf))
			      enddo
			  enddo
		      enddo
		  enddo
	      enddo
	   enddo
     	enddo


! ===================================================================================================== 
! ==================== Vertex Gluino - Squark (outgoing) - Quark (ingoing) =================
! =====================================================================================================


          do eigensf=1,2
             do ttsf = 1,4
               do gensf = 1,3

        JL(eigensf,ttsf,gensf) = - GS*sqrt2*USf(eigensf,1,ttsf,gensf)
        JR(eigensf,ttsf,gensf) =   GS*sqrt2*USf(eigensf,2,ttsf,gensf)

              enddo
             enddo
           enddo
  







  	end
  